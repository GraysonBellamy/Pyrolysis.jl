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
    w_override,
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
    if w_override !== nothing
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
    end
    return nothing
end

# Whether a gas-transfer coefficient is zero across the property bracket.
_vanishes(λ) = all(T -> λ(T) == 0, (T_PROPERTY_MIN, 300.0, 1000.0, T_PROPERTY_MAX))

# The gas transport a problem needs: none without a gas component, Darcy–Fick
# when a condensed component is permeable, Fickian when gas can move at all
# (some component carries a gas-transfer coefficient, or a boundary condition
# passes mass), none otherwise.
function _gas_transport(problem::PyrolysisProblem)
    material = problem.material
    any(is_gas, material.components) || return NoGasTransport()
    has_darcy_flow(material) && return DarcyFickTransport()
    moves =
        any(c -> !_vanishes(c.gas_transfer_coeff), material.components) ||
        any(bc -> !(bc isa ImpermeableBC), values(problem.bc_set.mass))
    return moves ? FickTransport() : NoGasTransport()
end

# The kernels in assembler order: reactions, gas transport (whose boundary
# mass fluxes the surface balance reads), in-depth radiation, conduction, the
# thermal boundary, and under a prescribed mesh velocity the inventory fluxes
# relative to the mesh.
function _kernels(
    problem,
    topo,
    g0,
    laws,
    gas::GasTransport,
    radiation_model,
    relative::Tuple,
    switch,
)
    material = problem.material
    reactions = n_reactions(material) > 0 ? (ReactionKernel(material; switch),) : ()
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

# The node velocities, and under a variable cross section the area and its
# lateral rate, that satisfy their algebraic rows at time `t`: one kernel pass
# for the strain rates and history rates, the lateral rows from the column
# sums, a second pass at the consistent lateral rate, then the sweep. Returns
# `nothing`, or the residual and scale of a velocity row left above round-off.
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
    ds, c = ws.derived, ws.contributions
    w = consistent_velocities!(copy(ds.w), cl, c, ds, def.setup.topology, t)
    for nd = 1:n_nodes(layout)
        u[velocity_index(layout, nd)] = w[nd]
    end
    residual!(du, u, def, t, ws)
    g = ds.geometry
    for nd = 1:n_nodes(layout)
        r = abs(du[velocity_index(layout, nd)])
        scale = sum(abs, c.θ .* ds.V_eval) + sum(abs, g.face_area .* w) + abs(w[nd])
        r <= 1e-10 * scale || return (r, scale)
    end
    return nothing
end

"""
	consistent_state!(u, def, ws, t) -> u

Overwrite the algebraic unknowns of `u` with the values that satisfy their
rows at time `t`, given its differential unknowns: the face intensities by
the attenuation recurrence from the exposed face, then the cross section and
the node velocities by the mesh closure. The integrator starts from such a
state; a velocity row left above round-off is an error.
"""
function consistent_state!(u::AbstractVector, def::CoreDefinition, ws::CoreWorkspace, t)
    _consistent_intensities!(u, def, ws, t)
    left = _consistent_motion!(u, def, ws, t)
    left === nothing || error(
        "the consistent node velocities leave a velocity row at $(left[1]) (scale $(left[2]))",
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
	      time_scale = problem.dt_initial, sources = NoSources(),
	      deactivation = false, band = threshold) -> CoreProblem

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
`w_override`, a function `w(x, t)`, prescribes the velocity of every node
instead, moving the mesh over a solid at rest (verification; the material
must not react), and the inventories cross the faces relative to the mesh
([`RelativeFluxKernel`](@ref)). `coordinates` reads the node positions as
radii of a cylinder (unit axial length) or a sphere (`Cylindrical()`,
`Spherical()`); planar meshes take the problem's cross-section area.
`time_scale` [s] scales the velocity tolerances (the first step).

`deactivation` deactivates, during the solve, every cell whose condensed mass
falls below the fraction `threshold` (`ε_m`) of its initial condensed mass
([`CellDeactivation`](@ref), [`ConservativeTransfer`](@ref)), from the
exposed `:top` side, on a moving mesh together with the run of cells outward
of it whose condensed mass is below the fraction `band` (`ε_band`, in
`[threshold, 1)`; at its default only the cells below the threshold) of
theirs ([`deactivate!`](@ref)); the reactions of a cell then switch off
smoothly below a few times `ε_m` ([`ReactivitySwitch`](@ref)).
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
    time_scale::Real = problem.dt_initial,
    sources = NoSources(),
    deactivation::Bool = false,
    band::Real = threshold,
)
    _require_supported(
        problem,
        use_ale,
        radiation_model,
        w_override,
        coordinates,
        deactivation,
    )
    problem = _own_problem(problem)
    apply_initial_conditions!(problem)
    material = problem.material
    NC = n_components(material)
    old = problem.mesh

    x = [node.z for node in old.nodes]
    mesh =
        coordinates isa Planar ?
        line_mesh(
            x;
            cross_section = old.cross_section_area_initial,
            bottom = :bottom,
            top = :top,
            fixed = !use_ale,
        ) : line_mesh(x; coordinates, bottom = :bottom, top = :top, fixed = !use_ale)
    topo = mesh.topology
    n = n_cells(topo)
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, x, topo.cross_section, topo)

    T0 = [Float64(old.cell_states[i].T) for i = 1:n]
    m0 = [Float64(old.cell_states[i].ξ[j]) * g0.cell_volume[i] for j = 1:NC, i = 1:n]

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
    closure = use_ale ? BalanceRowClosure(topo; override = w_override, lateral) : nothing
    closure === nothing ||
        rank_check(closure, topo, g0) ||
        throw(ArgumentError("the balance rows do not determine the node velocities"))
    relative =
        w_override === nothing ? () :
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
    switch =
        cells === nothing ? nothing :
        ReactivitySwitch(cells.threshold, cells.reference, cells.solid)
    kernels = _kernels(problem, topo, g0, laws, gas, radiation_model, relative, switch)
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
