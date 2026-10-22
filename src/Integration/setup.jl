# ============================================================================
# From a PyrolysisProblem to the core's definition, workspace and state
# ============================================================================

"""
	CoreProblem

A `PyrolysisProblem` set up on the finite-volume core: the problem (with its
own copy of the mesh), the definition and workspace the residual reads, the
initial state, the relative tolerance and the absolute-tolerance vector, the
colored Jacobian, the initial temperatures and masses, and the
[`CellDeactivation`](@ref) state (`nothing` when cells are not deactivated).
"""
struct CoreProblem{P,D<:CoreDefinition,W<:CoreWorkspace,J<:ColoredJacobianCache,X}
    problem::P
    def::D
    ws::W
    u0::Vector{Float64}
    reltol::Float64
    abstol::Vector{Float64}
    jacobian::J
    T0::Vector{Float64}
    m0::Matrix{Float64}
    deactivation::X
end

_energy_state(s::Symbol) =
    s === :temperature ? TemperatureState() :
    s === :enthalpy ? EnthalpyState() :
    throw(ArgumentError("energy_state must be :temperature or :enthalpy, got :$s"))

# A copy of the problem with its own mesh, so applying the initial conditions
# never touches the caller's problem.
_own_problem(p::PyrolysisProblem) = PyrolysisProblem(
    deepcopy(p.mesh),
    p.material,
    p.bc_set,
    p.T_initial,
    p.ξ_initial,
    p.tspan,
    p.dt_initial,
    p.output_times,
)

# The physics this core implements. A problem that needs anything else is
# refused rather than solved without it.
function _require_supported(
    problem::PyrolysisProblem,
    use_ale::Bool,
    radiation_model,
    correction::Correction,
    coordinates::CoordinateSystem,
    deactivation::Bool,
)
    material = problem.material
    radiation_model in (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT) || throw(
        ArgumentError("the finite-volume core has no $(radiation_model) radiation model"),
    )
    if radiation_model === BEER_LAMBERT && !(coordinates isa Planar)
        # The intensity rows attenuate the intensity and the cell sources read
        # the face areas, which agree only where the areas are equal.
        throw(ArgumentError("in-depth radiation needs planar coordinates"))
    end
    if material.lateral_shrinkage_law !== nothing
        use_ale || throw(
            ArgumentError(
                "material $(material.name) has a lateral_shrinkage_law, which needs a " *
                "moving mesh: solve with use_ale = true",
            ),
        )
        coordinates isa Planar ||
            throw(ArgumentError("a lateral shrinkage law needs planar coordinates"))
    end
    if correction isa PrescribedCorrection
        use_ale || throw(ArgumentError("a prescribed mesh velocity needs use_ale = true"))
        n_reactions(material) == 0 || throw(
            ArgumentError(
                "a prescribed mesh velocity moves the mesh over a solid at rest, so the " *
                "material must not react (it has $(n_reactions(material)) reactions)",
            ),
        )
        deactivation && throw(
            ArgumentError("cells are not deactivated under a prescribed mesh velocity"),
        )
    elseif !(correction isa NoCorrection)
        use_ale || throw(
            ArgumentError(
                "a correction moves the mesh relative to the solid, which needs " *
                "use_ale = true",
            ),
        )
    end
    return nothing
end

# The interior nodes between two cells that start with different condensed
# components (`amounts`, `NC × n`, their initial masses or concentrations):
# the layer interfaces an equidistribution correction must not move across.
function _composition_interfaces(
    amounts::AbstractMatrix{Float64},
    material::Material{NC},
) where {NC}
    present(i) = [amounts[j, i] > 0 && !is_gas(material.components[j]) for j = 1:NC]
    return [nd for nd = 2:size(amounts, 2) if present(nd - 1) != present(nd)]
end

# The correction of a solve that names none: the equidistribution correction
# where the mesh moves, cells are deactivated and some layer can vanish
# ([`vanishes`](@ref)), so that its cells thin together and leave in one
# event; no correction otherwise. `ξ0` holds the cells' initial
# concentrations (`NC × n`).
function _installed_correction(
    material::Material{NC},
    ξ0::AbstractMatrix{Float64},
    use_ale::Bool,
    deactivation::Bool,
) where {NC}
    (use_ale && deactivation) || return NoCorrection()
    v = vanishing(material)
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    any(i -> vanishes(v, ξ0, i:i, condensed), axes(ξ0, 2)) || return NoCorrection()
    return EquidistributionCorrection()
end

# An equidistribution correction needs every layer interface among the
# Lagrangian nodes, or it mixes the layers' compositions.
function _require_layers(::Correction, m0, material, lagrangian) end
function _require_layers(::EquidistributionCorrection, m0, material, lagrangian)
    loose = setdiff(_composition_interfaces(m0, material), lagrangian)
    isempty(loose) || throw(
        ArgumentError(
            "the cells beside node(s) $(loose) start with different condensed " *
            "components: give the layer interfaces as `lagrangian` nodes, so the " *
            "correction does not move the mesh across them",
        ),
    )
    return nothing
end

"""
	CORRECTION_BAND

The default near-threshold band `ε_band` of a solve whose bound
[`EquidistributionCorrection`](@ref) moves some layer: cells next to the
exposed side whose condensed mass is below this fraction of their reference
go with a cell that crosses the threshold, so a layer that the correction
carries to burnout together leaves in one invocation.
"""
const CORRECTION_BAND = 0.1

# Whether the bound correction of `closure` moves some layer.
_moves(closure) =
    closure !== nothing &&
    closure.correction isa LayerCorrection &&
    any(!iszero, closure.correction.gain)

# The gas transport a problem needs: none without a gas component, Darcy–Fick
# when a condensed component is permeable, Fickian when gas can move at all
# (some component carries a gas-transfer coefficient, or a boundary condition
# passes mass), none otherwise.
function _gas_transport(problem::PyrolysisProblem)
    material = problem.material
    any(is_gas, material.components) || return NoGasTransport()
    has_darcy_flow(material) && return DarcyFickTransport()
    moves =
        any(c -> !identically_zero(c.gas_transfer_coeff), material.components) ||
        any(bc -> !(bc isa ImpermeableBC), values(problem.bc_set.mass))
    return moves ? FickTransport() : NoGasTransport()
end

# The kernels in assembler order: reactions, gas transport (whose boundary
# mass fluxes the surface balance reads), in-depth radiation, conduction, the
# thermal boundary, and under a correction the inventory fluxes relative to
# the mesh.
function _kernels(
    problem,
    topo,
    g0,
    laws,
    gas::GasTransport,
    radiation_model,
    relative::Tuple,
)
    material = problem.material
    reactions = n_reactions(material) > 0 ? (ReactionKernel(material),) : ()
    transport =
        gas isa NoGasTransport ? () :
        (GasTransportKernel(gas, topo, g0, problem.bc_set, material, laws),)
    radiation =
        radiation_model === BEER_LAMBERT ? (BeerLambertKernel(topo, g0, problem.bc_set),) :
        ()
    boundary = ThermalBoundaryKernel(
        topo,
        g0,
        problem.bc_set;
        radiation_model,
        laws = gas isa NoGasTransport ? nothing : laws,
    )
    return (
        reactions...,
        transport...,
        radiation...,
        ConductionKernel(topo, g0),
        boundary,
        relative...,
    )
end

# The largest incident intensity over the solve window [W/m²], at least 1: the
# scale of the face-intensity tolerances.
function _intensity_scale(kernels::Tuple, ds, tspan)
    kernel = find_kernel(BeerLambertKernel, kernels)
    kernel === nothing && return 1.0
    times = range(Float64(tspan[1]), Float64(tspan[2]); length = 101)
    return max(1.0, maximum(t -> abs(incident_intensity(kernel, ds, t)), times))
end

# The closure's unknowns, and under a variable cross section the area and its
# lateral rate, that satisfy their algebraic rows at time `t`: one kernel pass
# for the strain rates and history rates, the lateral rows from the column
# sums, a second pass at the consistent lateral rate, then the closure's
# unknowns. Returns `nothing`, or the residual and scale of a closure row left
# above round-off.
function _consistent_motion!(u, def::CoreDefinition, ws::CoreWorkspace, t)
    cl = def.closure
    cl === nothing && return nothing
    layout = def.layout
    du = similar(u)
    residual!(du, u, def, t, ws)
    if has_chi(layout)
        ΣG, ΣM, ΣĠ, ΣṀ = 0.0, 0.0, 0.0, 0.0
        for i = 1:n_cells(layout)
            kG, kM = history_index(layout, i, 1), history_index(layout, i, 2)
            ΣG += u[kG]
            ΣM += u[kM]
            ΣĠ += du[kG]
            ΣṀ += du[kM]
        end
        A, θ_A = consistent_lateral(cl.lateral, ΣG, ΣM, ΣĠ, ΣṀ)
        u[area_index(layout)] = A
        u[lateral_rate_index(layout)] = θ_A
        residual!(du, u, def, t, ws)
    end
    return _consistent_unknowns!(u, du, cl, def, ws, t)
end

# The node velocities by the closure's sweep.
function _consistent_unknowns!(u, du, cl::BalanceRowClosure, def, ws, t)
    layout = def.layout
    ds, c = ws.derived, ws.contributions
    D = dimension(layout)
    w = consistent_velocities!(copy(ds.w), cl, c, ds, def.setup.topology, t)
    for nd = 1:n_nodes(layout), k = 1:D
        u[velocity_index(layout, nd, k)] = w[D*(nd-1)+k]
    end
    residual!(du, u, def, t, ws)
    return _left_velocity(du, layout, ds, c, w)
end

# The first velocity row left above round-off at the consistent velocities, and
# its scale: on a line the terms of its balance, on a column the largest.
function _left_velocity(
    du,
    layout::StateLayout{NC,E,Mov,Rad,Chi,1},
    ds,
    c,
    w,
) where {NC,E,Mov,Rad,Chi}
    g = ds.geometry
    for nd = 1:n_nodes(layout)
        r = abs(du[velocity_index(layout, nd)])
        scale =
            sum(abs, c.θ .* ds.V_eval) +
            sum(abs, g.face_area .* w) +
            sum(abs, g.face_area .* ds.c) +
            abs(w[nd])
        r <= 1e-10 * scale || return (r, scale)
    end
    return nothing
end

function _left_velocity(du, layout::StateLayout, ds, c, w)
    D = dimension(layout)
    base = sum(abs, c.θ .* ds.V_eval) + sum(ds.geometry.face_area) * maximum(abs, w)
    for nd = 1:n_nodes(layout), k = 1:D
        r = abs(du[velocity_index(layout, nd, k)])
        scale = base + abs(w[D*(nd-1)+k])
        r <= 1e-10 * scale || return (r, scale)
    end
    return nothing
end

# The multipliers by one symmetric positive-definite solve.
function _consistent_unknowns!(u, du, cl::PotentialClosure, def, ws, t)
    layout = def.layout
    ds, c = ws.derived, ws.contributions
    φ = consistent_multipliers(cl, c, ds, def.setup.topology)
    for i = 1:n_cells(layout)
        u[multiplier_index(layout, i)] = φ[i]
    end
    residual!(du, u, def, t, ws)
    scale = sum(abs, c.θ .* ds.V_eval) + sum(ds.geometry.face_area) * maximum(abs, ds.w)
    for i = 1:n_cells(layout)
        r = abs(du[multiplier_index(layout, i)])
        r <= 1e-10 * scale || return (r, scale)
    end
    return nothing
end

"""
	consistent_state!(u, def, ws, t) -> u

Overwrite the algebraic unknowns of `u` with the values that satisfy their
rows at time `t`, given its differential unknowns: the face intensities by
the attenuation recurrence from the exposed face, then the cross section and
the node velocities or the cells' multipliers by the mesh closure. The
integrator starts from such a state; a closure row left above round-off is
an error.
"""
function consistent_state!(u::AbstractVector, def::CoreDefinition, ws::CoreWorkspace, t)
    _consistent_intensities!(u, def, ws, t)
    left = _consistent_motion!(u, def, ws, t)
    left === nothing || error(
        "the closure's consistent unknowns leave a row at $(left[1]) (scale $(left[2]))",
    )
    return u
end

# The face intensities of `u` that satisfy their algebraic rows at time `t`:
# the consistent values the integrator starts from.
function _consistent_intensities!(u, def::CoreDefinition, ws::CoreWorkspace, t)
    layout = def.layout
    kernel = find_kernel(BeerLambertKernel, def.kernels)
    kernel === nothing && return u
    ds = ws.derived
    derive!(ds, u, layout, def.setup)
    sweep_intensities!(ds.I, kernel, ds, def.setup.topology, t)
    for f in eachindex(ds.I)
        u[intensity_index(layout, f)] = ds.I[f]
    end
    return u
end

"""
	setup(problem; energy_state = :temperature, reltol = 1e-6, abstol = 1e-8,
	      threshold = 1e-2, use_ale = false, radiation_model = NO_RADIATION,
	      coordinates = Planar(), w_override = nothing,
	      correction = nothing, lagrangian = nothing,
	      time_scale = problem.dt_initial, sources = NoSources(),
	      deactivation = false, band = nothing) -> CoreProblem

Set up `problem` on the finite-volume core: apply its initial conditions to its
own copy of the mesh, build the 1D mesh from the node positions and the
cross-section area, the initial inventories `mᵢⱼ = ξⱼ(zᵢ) Vᵢ`, the energy
state from `T(zᵢ)`, the gas transport the material and the mass conditions
need, the kernels, the face intensities consistent with the incident flux
under `BEER_LAMBERT`, the tolerance vector and the colored Jacobian.

`radiation_model` is `NO_RADIATION` or `SURFACE_ABSORPTION` (the incident
flux enters the surface balance at the exposed cell's absorptivity, or the
condition's own) or `BEER_LAMBERT` (it is absorbed in depth, through the
face intensities of [`BeerLambertKernel`](@ref)).

`use_ale` moves the mesh with the solid: the back node is pinned and every
other node follows the cells' reaction strain rates through the balance rows
of [`BalanceRowClosure`](@ref), with the node velocities consistent at the
start; a material with a lateral shrinkage law then also carries the
conversion-history pair and its cross section's area and lateral rate.
`correction` moves the mesh relative to the solid
([`EquidistributionCorrection`](@ref); [`NoCorrection`](@ref) keeps it with
the solid), the balance rows then holding the solid velocity `w − c`, and the
inventories cross the faces at the relative velocity `−c`
([`RelativeFluxKernel`](@ref)). By default a solve gets the equidistribution
correction when the mesh moves, cells are deactivated and some layer can
vanish (its condensed components gasify completely, [`vanishing`](@ref)), and
no correction otherwise; a layer that cannot vanish keeps its mesh with the
solid under any equidistribution correction. `lagrangian` lists the mesh
nodes that always move with the solid, the layer interfaces; under an
equidistribution correction every interior node between cells that start
with different condensed components must be among them, and by default they
are those nodes. `w_override`, a function `w(x, t)`, prescribes the velocity
of every node instead
([`PrescribedCorrection`](@ref)), moving the mesh over a solid at rest
(verification; the material must not react, and no other correction is
given), the inventories crossing the faces at `−w`. `coordinates` reads the
node positions as
radii of a cylinder (unit axial length) or a sphere (`Cylindrical()`,
`Spherical()`); planar meshes take the problem's cross-section area.
`time_scale` [s] scales the velocity tolerances (the first step).

`deactivation` deactivates, during the solve, every cell whose condensed mass
falls below the fraction `threshold` (`ε_m`) of its initial condensed mass
([`CellDeactivation`](@ref), [`ConservativeTransfer`](@ref)), from the
exposed `:top` side, on a moving mesh together with the run of cells outward
of it whose condensed mass is below the fraction `band` (`ε_band`, in
`[threshold, 1)`; by default the larger of the threshold and
[`CORRECTION_BAND`](@ref) when the bound equidistribution correction moves
some layer, otherwise the threshold, so that only the cells below it go) of
theirs
([`deactivate!`](@ref)).
`threshold` also scales the mass tolerances; `sources` is the external-source
hook.
"""
function setup(
    problem::PyrolysisProblem;
    energy_state::Symbol = :temperature,
    reltol::Real = 1e-6,
    abstol::Real = 1e-8,
    threshold::Real = 1e-2,
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    coordinates::CoordinateSystem = Planar(),
    w_override = nothing,
    correction::Union{Nothing,Correction} = nothing,
    lagrangian::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    time_scale::Real = problem.dt_initial,
    sources = NoSources(),
    deactivation::Bool = false,
    band::Union{Nothing,Real} = nothing,
)
    if w_override !== nothing
        (correction === nothing || correction isa NoCorrection) || throw(
            ArgumentError("w_override prescribes the mesh velocity: give no correction"),
        )
        correction = PrescribedCorrection(w_override)
    end
    problem = _own_problem(problem)
    apply_initial_conditions!(problem)
    material = problem.material
    NC = n_components(material)
    old = problem.mesh
    ξ0 = [Float64(s.ξ[j]) for j = 1:NC, s in old.cell_states]
    correction =
        something(correction, _installed_correction(material, ξ0, use_ale, deactivation))
    lagrangian = something(
        lagrangian,
        correction isa EquidistributionCorrection ?
        _composition_interfaces(ξ0, material) : Int[],
    )
    _require_supported(
        problem,
        use_ale,
        radiation_model,
        correction,
        coordinates,
        deactivation,
    )
    x = [node.z for node in old.nodes]
    mesh =
        coordinates isa Planar ?
        line_mesh(
            x;
            cross_section = old.cross_section_area_initial,
            bottom = :bottom,
            top = :top,
            fixed = !use_ale,
            lagrangian,
        ) :
        line_mesh(
            x;
            coordinates,
            bottom = :bottom,
            top = :top,
            fixed = !use_ale,
            lagrangian,
        )
    topo = mesh.topology
    n = n_cells(topo)
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, x, topo.cross_section, topo)

    T0 = [Float64(old.cell_states[i].T) for i = 1:n]
    m0 = [Float64(old.cell_states[i].ξ[j]) * g0.cell_volume[i] for j = 1:NC, i = 1:n]
    _require_layers(correction, m0, material, topo.lagrangian)

    reference = CellReference(
        mesh,
        m0,
        material;
        reltol = Float64(reltol),
        threshold = Float64(threshold),
    )
    energy = _energy_state(energy_state)
    laws = CaloricLaws(material)
    caloric =
        energy isa EnthalpyState ?
        CaloricExtension(laws, reference.mass_scale, reference.capacity_floor) : nothing
    gas = _gas_transport(problem)
    state_setup = StateSetup(topo, material, laws, reference, caloric, gas)
    chi = material.lateral_shrinkage_law !== nothing
    layout = StateLayout(
        NC,
        n,
        energy;
        moving = use_ale,
        intensities = radiation_model === BEER_LAMBERT,
        chi,
    )

    M0 = [
        sum(m0[j, i] for j = 1:NC if is_solid(material.components[j]); init = 0.0) for
        i = 1:n
    ]
    condensed0 = [
        sum(m0[j, i] for j = 1:NC if !is_gas(material.components[j]); init = 0.0) for
        i = 1:n
    ]
    lateral = chi ? LateralShrinkage(material, topo.cross_section, sum(M0)) : nothing
    closure =
        use_ale ?
        BalanceRowClosure(
            topo;
            correction = bind_correction(correction, topo, reference, m0, material),
            lateral,
        ) : nothing
    closure === nothing ||
        rank_check(closure, topo, g0) ||
        throw(ArgumentError("the balance rows do not determine the node velocities"))
    relative =
        correction isa NoCorrection ? () :
        (
            RelativeFluxKernel(
                topo,
                g0,
                material,
                laws,
                reference;
                history = chi,
                temperature_scale = Float64(abstol),
                history_scale = Float64(reltol * threshold),
            ),
        )
    band = something(band, _moves(closure) ? max(threshold, CORRECTION_BAND) : threshold)
    cells =
        deactivation ?
        CellDeactivation(
            ConservativeTransfer(),
            threshold,
            band,
            condensed0,
            map(c -> !is_gas(c), material.components),
            something(
                caloric,
                CaloricExtension(laws, reference.mass_scale, reference.capacity_floor),
            ),
        ) : nothing
    kernels = _kernels(problem, topo, g0, laws, gas, radiation_model, relative)
    def = CoreDefinition(layout, state_setup, kernels, closure, sources)
    ws = CoreWorkspace(def, T0)

    u0 = zeros(state_length(layout))
    for i = 1:n
        mi = ntuple(j -> m0[j, i], NC)
        for j = 1:NC
            u0[mass_index(layout, i, j)] = mi[j]
        end
        u0[energy_index(layout, i)] = energy_state_value(energy, state_setup, T0[i], mi, i)
        chi && (u0[history_index(layout, i, 2)] = M0[i])
    end
    if use_ale
        for nd = 1:n_nodes(layout)
            u0[position_index(layout, nd)] = x[nd]
        end
        chi && (u0[area_index(layout)] = topo.cross_section)
    end
    t0 = problem.tspan[1]
    consistent_state!(u0, def, ws, t0)
    tol = abstol_vector(
        layout,
        state_setup,
        m0,
        T0;
        reltol = Float64(reltol),
        abstol = Float64(abstol),
        threshold = Float64(threshold),
        intensity_scale = _intensity_scale(kernels, ws.derived, problem.tspan),
        time_scale = Float64(time_scale),
    )
    cj = ColoredJacobianCache(def, ws, u0, t0)
    return CoreProblem(problem, def, ws, u0, Float64(reltol), tol, cj, T0, m0, cells)
end

# ============================================================================
# Meshes of more than one dimension
# ============================================================================

"""
	setup(problem::PyrolysisProblem{<:Mesh{D}}; energy_state = :temperature,
	      reltol = 1e-6, abstol = 1e-8, threshold = 1e-2,
	      radiation_model = NO_RADIATION, use_ale = false, closure = :potential,
	      slip = (), w_override = nothing, tensor = nothing, history = false,
	      deformation = false, time_scale = problem.dt_initial,
	      sources = NoSources()) -> CoreProblem    (D = 2, 3)

Set up a problem whose mesh is a mesh of quadrilaterals or hexahedra
([`box_mesh`](@ref), [`cell_mesh`](@ref)): the initial inventories
`mᵢⱼ = ξⱼ(xᵢ) Vᵢ` and the energy state from `T(xᵢ)` at the cell centroids
`xᵢ` (the problem's initial-condition functions take an `SVector{D}`),
the reaction, conduction and thermal-boundary kernels, the tolerance vector
and the colored Jacobian. The thermal conditions are keyed by the mesh's
boundary tags.

`tensor` is each cell's dimensionless conductivity tensor, an
`SMatrix{D,D}` for every cell or a function of the centroid returning one
(the identity by default): the conductivity is `k` times it.

`use_ale` moves the mesh with the solid. The closure is `:potential`, one
multiplier per cell with the node velocities of minimum lumped-volume norm
that satisfy every cell's volume balance ([`PotentialClosure`](@ref)), the
pinned nodes at rest and on each `tag => k` of `slip` the component `k`
zero; or `:balance`, the balance rows of a single column of cells
([`BalanceRowClosure`](@ref)). `w_override`, a function
`w(x::SVector, t) -> SVector`, prescribes every node's velocity instead,
moving the mesh over a solid at rest (the material must not react), the
inventories crossing the faces at `−w` ([`RelativeFluxKernel`](@ref)).

`history` carries the conversion-history pair (`G = 0`, `M` the cell's
initial solid mass) and `deformation` the deformation gradient (the
identity at the start) on a moving mesh.

Gas transport, in-depth radiation, the variable cross section, the
equidistribution correction and deactivation are a line's: a problem that
needs one is refused.
"""
function setup(
    problem::PyrolysisProblem{<:Mesh{D}};
    energy_state::Symbol = :temperature,
    reltol::Real = 1e-6,
    abstol::Real = 1e-8,
    threshold::Real = 1e-2,
    radiation_model = NO_RADIATION,
    use_ale::Bool = false,
    closure::Symbol = :potential,
    slip = (),
    w_override = nothing,
    tensor = nothing,
    history::Bool = false,
    deformation::Bool = false,
    time_scale::Real = problem.dt_initial,
    sources = NoSources(),
) where {D}
    D in (2, 3) || throw(ArgumentError("a line's mesh is set up from a `Mesh1D`"))
    material = problem.material
    NC = n_components(material)
    radiation_model in (NO_RADIATION, SURFACE_ABSORPTION) ||
        throw(ArgumentError("in-depth radiation is a line's: the mesh has $D dimensions"))
    _gas_transport(problem) isa NoGasTransport ||
        throw(ArgumentError("gas transport is a line's: the mesh has $D dimensions"))
    material.lateral_shrinkage_law === nothing || throw(
        ArgumentError("the variable cross section is a line's: the mesh has $D dimensions"),
    )
    closure in (:potential, :balance) ||
        throw(ArgumentError("closure must be :potential or :balance, got :$closure"))
    w_override === nothing ||
        use_ale ||
        throw(ArgumentError("a prescribed mesh velocity needs use_ale = true"))
    if w_override !== nothing && n_reactions(material) > 0
        throw(
            ArgumentError(
                "a prescribed mesh velocity moves the mesh over a solid at rest, so the " *
                "material must not react (it has $(n_reactions(material)) reactions)",
            ),
        )
    end
    mesh = problem.mesh
    topo = mesh.topology
    is_fixed(topo) == !use_ale || throw(
        ArgumentError(
            use_ale ? "use_ale = true needs a mesh whose nodes are not all pinned" :
            "a fixed solve needs a mesh with every node pinned (`fixed = true`)",
        ),
    )
    n = n_cells(topo)
    x = reduce(vcat, mesh.positions)
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, x, topo.cross_section, topo)
    T0 = [Float64(problem.T_initial(g0.cell_center[i])) for i = 1:n]
    m0 = [
        Float64(problem.ξ_initial[j](g0.cell_center[i])) * g0.cell_volume[i] for
        j = 1:NC, i = 1:n
    ]

    reference = CellReference(
        mesh,
        m0,
        material;
        reltol = Float64(reltol),
        threshold = Float64(threshold),
    )
    energy = _energy_state(energy_state)
    laws = CaloricLaws(material)
    caloric =
        energy isa EnthalpyState ?
        CaloricExtension(laws, reference.mass_scale, reference.capacity_floor) : nothing
    state_setup = StateSetup(topo, material, laws, reference, caloric, NoGasTransport())
    potential = use_ale && w_override === nothing && closure === :potential
    layout = StateLayout(
        NC,
        topo,
        energy;
        moving = use_ale,
        history,
        deformation,
        unknowns = potential ? :multipliers : :velocities,
    )
    cl =
        !use_ale ? nothing :
        w_override !== nothing ?
        BalanceRowClosure(topo; correction = PrescribedCorrection(w_override)) :
        potential ? PotentialClosure(topo; slip) : BalanceRowClosure(topo)
    cl === nothing ||
        rank_check(cl, topo, g0) ||
        throw(ArgumentError("the closure's rows do not determine the mesh's motion"))
    A = _tensors(tensor, g0, n, Val(D))
    reactions = n_reactions(material) > 0 ? (ReactionKernel(material),) : ()
    relative =
        w_override === nothing ? () :
        (
            RelativeFluxKernel(
                topo,
                g0,
                material,
                laws,
                reference;
                history,
                temperature_scale = Float64(abstol),
                history_scale = Float64(reltol * threshold),
            ),
        )
    kernels = (
        reactions...,
        ConductionKernel(topo, g0; tensor = A),
        ThermalBoundaryKernel(topo, g0, problem.bc_set; radiation_model, tensor = A),
        relative...,
    )
    def = CoreDefinition(layout, state_setup, kernels, cl, sources)
    ws = CoreWorkspace(def, T0)

    u0 = zeros(state_length(layout))
    for i = 1:n
        mi = ntuple(j -> m0[j, i], NC)
        for j = 1:NC
            u0[mass_index(layout, i, j)] = mi[j]
        end
        u0[energy_index(layout, i)] = energy_state_value(energy, state_setup, T0[i], mi, i)
        history && (
            u0[history_index(layout, i, 2)] = sum(
                m0[j, i] for j = 1:NC if is_solid(material.components[j]);
                init = 0.0,
            )
        )
        if deformation
            for a = 1:D
                u0[deformation_index(layout, i, a, a)] = 1.0
            end
        end
    end
    if use_ale
        for nd = 1:n_nodes(topo), k = 1:D
            u0[position_index(layout, nd, k)] = x[D*(nd-1)+k]
        end
    end
    t0 = problem.tspan[1]
    consistent_state!(u0, def, ws, t0)
    tol = abstol_vector(
        layout,
        state_setup,
        m0,
        T0;
        reltol = Float64(reltol),
        abstol = Float64(abstol),
        threshold = Float64(threshold),
        time_scale = Float64(time_scale),
    )
    cj = ColoredJacobianCache(def, ws, u0, t0)
    return CoreProblem(problem, def, ws, u0, Float64(reltol), tol, cj, T0, m0, nothing)
end

# Every cell's dimensionless conductivity tensor.
_tensors(::Nothing, g0, n, ::Val{D}) where {D} = fill(SMatrix{D,D,Float64}(I), n)
_tensors(A::SMatrix, g0, n, ::Val{D}) where {D} = fill(SMatrix{D,D,Float64}(A), n)
_tensors(A::AbstractVector, g0, n, ::Val{D}) where {D} =
    [SMatrix{D,D,Float64}(a) for a in A]
_tensors(A::Function, g0, n, ::Val{D}) where {D} =
    [SMatrix{D,D,Float64}(A(g0.cell_center[i])) for i = 1:n]
