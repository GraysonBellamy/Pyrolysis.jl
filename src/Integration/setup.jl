# ============================================================================
# From a PyrolysisProblem to its definition, workspace and state
# ============================================================================

"""
	PreparedProblem

A `PyrolysisProblem` set up for integration: the problem, the
definition and workspace the residual reads, the
initial state, the relative tolerance and the absolute-tolerance vector, the
colored Jacobian, the initial temperatures and masses, and the
[`CellDeactivation`](@ref) state (`nothing` when cells are not deactivated).
"""
struct PreparedProblem{
    P,
    D<:ResidualDefinition,
    W<:ResidualWorkspace,
    J<:ColoredJacobianCache,
    X,
}
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

# The physics implemented on the problem's mesh. A problem that needs anything
# else is refused rather than solved without it.
function _require_supported(
    problem::PyrolysisProblem{<:Mesh{D}},
    use_ale::Bool,
    radiation_model,
    correction::Correction,
    deactivation::Bool;
    closure::Symbol,
    slip,
    tensor,
    history::Bool,
    deformation::Bool,
) where {D}
    material = problem.material
    coordinates = problem.mesh.topology.coordinates
    if D == 1
        closure === :balance ||
            throw(ArgumentError("a line's mesh moves by its balance rows: give no closure"))
        for (name, given) in (
            :slip => !isempty(slip),
            :tensor => tensor !== nothing,
            :history => history,
            :deformation => deformation,
        )
            given &&
                throw(ArgumentError("`$name` applies to meshes of two or three dimensions"))
        end
    else
        closure in (:potential, :balance) ||
            throw(ArgumentError("closure must be :potential or :balance, got :$closure"))
        line = (
            "in-depth radiation" => radiation_model !== SURFACE_ABSORPTION,
            "gas transport" => gas_can_move(problem),
            "the variable cross section" => material.lateral_shrinkage_law !== nothing,
            "deactivation" => deactivation,
            "the equidistribution correction" =>
                !(correction isa Union{NoCorrection,PrescribedCorrection}),
            "a reaction order other than one" => needs_lineage(material),
        )
        for (what, needed) in line
            needed && throw(ArgumentError("$what is a line's: the mesh has $D dimensions"))
        end
    end
    if radiation_model !== SURFACE_ABSORPTION && !(coordinates isa Planar)
        # The radiation rows carry fluxes per unit face area and the cell
        # sources read the face areas, which agree only where the areas are equal.
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
    if needs_lineage(material)
        name = first(r for r in material.reactions if r.order != 1).name
        correction isa NoCorrection || throw(
            ArgumentError(
                "reaction :$name has an order other than one, whose cells must keep " *
                "their material, which a correction moves between them: give no correction",
            ),
        )
    end
    return nothing
end

# The mesh's topology with the nodes a solve moves: without `use_ale` every
# node pinned; with it the mesh's own pinned nodes, which must leave some
# free. `lagrangian` joins the mesh's Lagrangian nodes.
function _motion(
    topo::MeshTopology{D,C,FT,CT},
    use_ale::Bool,
    lagrangian::AbstractVector{<:Integer},
) where {D,C,FT,CT}
    nn = n_nodes(topo)
    use_ale &&
        is_fixed(topo) &&
        throw(ArgumentError("use_ale = true needs a mesh whose nodes are not all pinned"))
    lag = sort!(unique!(vcat(topo.lagrangian, Int.(lagrangian))))
    all(n -> 1 <= n <= nn, lag) ||
        throw(ArgumentError("lagrangian node ids must lie in 1:$nn"))
    return MeshTopology{D,C,FT,CT}(
        topo.coordinates,
        topo.faces,
        topo.cells,
        nn,
        topo.boundary_tags,
        use_ale ? topo.pinned : collect(1:nn),
        lag,
        topo.cross_section,
    )
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

# Whether the solve's sources add mass to the cells.
_external_mass(::NoSources) = false
_external_mass(s::VolumetricSources) = s.mass !== nothing

"""
	lineage_references(graph, m0, condensed0, sources) -> references

The reference of every reaction of order other than one, from the rules of the
[`LineageGraph`](@ref) `graph` and the initial fields: `references[r]` holds the
per-cell reference [kg] `c_{r,i} M_{0,i}` of reaction `r` (empty for a first-order
reaction), zero where the cell holds none of the lineage. `m0` holds the initial
masses (`NC × n`) and `condensed0` the initial condensed mass `M_{0,i}` of each
cell.

With `Y_{v,i} = m_{v,i}/M_{0,i}`, exactly zero in a cell that starts empty, the
coefficient is `c_{r,i} = Σ_v Y_{v,i} Φ(v)` over the reactant's ancestors, `Φ` the
graph's reach, or under a `basis` `b` the value `b` on the cells holding some
ancestor. Throws an `ArgumentError` when the sources add mass (every condensed
node is then externally fed), when a cell holds a reactant that is also formed
there and the reaction has no `basis`, and when a positive coefficient is below
`1e-12` or a `basis` below `1e-6`, so that every positive reference is a normal
float and exact zero means absence.
"""
function lineage_references(
    graph::LineageGraph,
    m0::AbstractMatrix{Float64},
    condensed0::AbstractVector{Float64},
    sources,
)
    NC, n = size(m0)
    converted = findall(!=(1.0), graph.order)
    (_external_mass(sources) && !isempty(converted)) && throw(
        ArgumentError(
            "reaction :$(graph.names[first(converted)]) has an order other than one " *
            "in a solve whose sources add mass: its lineage has no initial inventory",
        ),
    )
    Y = zeros(NC, n)
    for i = 1:n
        condensed0[i] > 0 || continue
        for j = 1:NC
            graph.condensed[j] && (Y[j, i] = m0[j, i] / condensed0[i])
        end
    end
    references = [Float64[] for _ in graph.order]
    for r in converted
        a, b = graph.reactant[r], graph.basis[r]
        name = graph.names[r]
        !isnan(b) &&
            b < 1e-6 &&
            throw(
                ArgumentError(
                    "reaction :$name has basis = $b: a basis below 1e-6 is refused",
                ),
            )
        c = zeros(n)
        for i = 1:n
            if isnan(b)
                Φ = graph.reach[r]
                formed = sum(Y[v, i] * Φ[v] for v = 1:NC if v != a; init = 0.0)
                (Y[a, i] > 0 && formed > 0) && throw(
                    ArgumentError(
                        "reaction :$name: cell $i holds :$(graph.components[a]) and the " *
                        "ancestors it is formed from, so its reference is ambiguous: " *
                        "give `basis`",
                    ),
                )
                c[i] = Y[a, i] + formed
            else
                c[i] = any(v -> graph.ancestors[r][v] && Y[v, i] > 0, 1:NC) ? b : 0.0
            end
            0 < c[i] < 1e-12 && throw(
                ArgumentError(
                    "reaction :$name: the lineage fraction $(c[i]) of cell $i is below " *
                    "1e-12, too small to be a normalisation",
                ),
            )
        end
        references[r] = c .* condensed0
    end
    return references
end

# The reference of every reaction in every cell [kg], reactions × cells: zero
# for a first-order reaction and throughout an all-first-order material.
function _reaction_references(
    material::Material{NC,NR},
    m0,
    condensed0,
    sources,
) where {NC,NR}
    table = zeros(NR, size(m0, 2))
    needs_lineage(material) || return table
    references = lineage_references(lineage_graph(material), m0, condensed0, sources)
    for (r, reference) in enumerate(references)
        isempty(reference) || (table[r, :] = reference)
    end
    return table
end

# The fraction of its initial volume a cell may keep once fully converted and
# still leave at an order other than one. A cell that keeps more is a layer the
# front burns beneath: merged into its neighbour it would stop resolving that layer,
# and the reactant then taken in would burn at the temperature of the middle of it.
const VANISHING_VOLUME = 1e-2

# The fraction of its initial volume each cell keeps once its components are
# converted completely ([`retained_volumes`](@ref)), `m0` its initial masses
# (`NC × n`) and `T0` its initial temperatures. Zero in a cell with no condensed mass.
function _retained_volume(material::Material{NC}, m0, T0) where {NC}
    return map(axes(m0, 2)) do i
        v = retained_volumes(material, T0[i])
        kept, initial = 0.0, 0.0
        for j = 1:NC
            is_gas(material.components[j]) && continue
            kept += m0[j, i] * v[j]
            initial += m0[j, i] / material.components[j].density(T0[i])
        end
        initial > 0 ? kept / initial : 0.0
    end
end

# What the deactivation threshold of `material` reads at an order other than one:
# the component of each reaction with a condensed reactant, and each cell's
# reactive reference, the sum over those reactions of the reaction's lineage
# reference or, at order one, where the lineage is not defined for every network,
# its reactant's initial mass, zero in a cell that keeps `VANISHING_VOLUME` or more
# of its volume, which is then never tested. Nothing at order one.
function _reactive_references(material::Material, references::Matrix{Float64}, m0, T0)
    needs_lineage(material) || return Int[], Float64[]
    reactants = Int[]
    reactive = zeros(size(m0, 2))
    for (r, reaction) in enumerate(material.reactions)
        a = reaction.reactant
        is_gas(material.components[a]) && continue
        push!(reactants, a)
        reactive .+= reaction.order == 1 ? view(m0, a, :) : view(references, r, :)
    end
    reactive[_retained_volume(material, m0, T0) .>= VANISHING_VOLUME] .= 0.0
    return reactants, reactive
end

# The deactivation threshold of a solve that names none: below order one the
# rate dies only where the conversion onset acts, so a cell should leave once
# its conversion is a decade inside it.
_default_threshold(material::Material, deactivation::Bool) =
    deactivation && needs_lineage(material) ?
    material.depletion_limiter.conversion_threshold / 10 : 1e-2

# The correction of a solve that names none: the equidistribution correction
# where the mesh moves, cells are deactivated and deactivation at `threshold`
# can remove some cell ([`vanishes`](@ref)), so that a layer's cells thin
# together and leave in one event; no correction otherwise, and none at a
# reaction order other than one, whose cells keep their material. `ξ0` holds the
# cells' initial concentrations (`NC × n`).
function _installed_correction(
    material::Material{NC},
    ξ0::AbstractMatrix{Float64},
    use_ale::Bool,
    deactivation::Bool,
    threshold::Real,
) where {NC}
    (use_ale && deactivation && !needs_lineage(material)) || return NoCorrection()
    r = residues(material)
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    any(i -> vanishes(r, ξ0, i:i, condensed, threshold), axes(ξ0, 2)) ||
        return NoCorrection()
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
[`EquidistributionCorrection`](@ref) moves some layer: cells beside a cell
that crosses the threshold whose condensed mass is below this fraction of their
reference go with it, so a layer that the correction carries to burnout together
leaves in one invocation.
"""
const CORRECTION_BAND = 0.1

# Whether the bound correction of `closure` moves some layer: a gain already
# set, or layers whose times the solve estimates.
_moves(closure) =
    closure isa BalanceRowClosure &&
    closure.correction isa LayerCorrection &&
    (any(!iszero, closure.correction.gain) || estimates_layer_times(closure.correction))

# Pore state exists for every gas-bearing material, independently of mobility.
_gas_transport(problem::PyrolysisProblem) =
    has_gas(problem.material) ? problem.material.gas_transport : NoGasTransport()

function gas_can_move(problem::PyrolysisProblem)
    mat=problem.material
    has_gas(mat) || return false
    any(
        bc->!(bc isa ImpermeableBC) &&
            !(bc isa MassFluxBC && all(x->x isa Real && iszero(x), values(bc.fluxes))),
        values(problem.bc_set.mass),
    ) && return true
    return !closed_transport_channel(mat)
end

# The kernels in assembler order: reactions, gas transport (whose boundary
# mass fluxes the surface balance reads), in-depth radiation, conduction, the
# thermal boundary, and under a correction the inventory fluxes relative to
# the mesh. `tensors` holds every cell's dimensionless conductivity tensor on
# more than one dimension, `nothing` on a line.
function _kernels(
    problem,
    topo,
    g0,
    laws,
    gas::GasTransport,
    radiation_model,
    relative::Tuple,
    tensors,
    pressure_tolerance,
    references::Matrix{Float64},
)
    material = problem.material
    reactions = n_reactions(material) > 0 ? (ReactionKernel(material, references),) : ()
    transport =
        !gas_can_move(problem) ? () :
        (GasTransportKernel(gas, topo, g0, problem.bc_set, material, laws),)
    radiation =
        radiation_model === RADIATIVE_TRANSFER ?
        (RadiativeTransferKernel(topo, g0, problem.bc_set, material),) : ()
    boundary = ThermalBoundaryKernel(
        topo,
        g0,
        problem.bc_set;
        radiation_model,
        laws = gas isa NoGasTransport ? nothing : laws,
        tensor = tensors,
        material,
        pressure_tolerance,
    )
    return (
        reactions...,
        transport...,
        radiation...,
        ConductionKernel(topo, g0; tensor = tensors),
        boundary,
        relative...,
    )
end

# The radiation rows per face of a radiation model.
_radiation_rows(model::RadiationModel, bc_set) =
    model === RADIATIVE_TRANSFER ? radiative_transfer_rows(bc_set) : 0

# The scale of the radiation rows' tolerances [W/m²], at least 1: under
# radiative transfer the larger of the largest radiation arriving at either
# face over the solve window (diffuse members, the surroundings' F σ T∞⁴ and
# the beams), sampled at 101 times, and the emissive power of the hottest
# initial cell. The window is read whole, so a later end can raise the scale.
function _transfer_scale(kernels::Tuple, ds, tspan)
    kernel = find_kernel(RadiativeTransferKernel, kernels)
    kernel === nothing && return 1.0
    times = range(Float64(tspan[1]), Float64(tspan[2]); length = 101)
    conditions =
        kernel.back isa OpenBack ? (kernel.top, kernel.back.condition) : (kernel.top,)
    arriving = 0.0
    for t in times, bc in conditions
        arriving = max(arriving, abs(arriving_radiation(bc, t)))
    end
    return max(1.0, arriving, STEFAN_BOLTZMANN * maximum(ds.T)^4)
end

# The area and its lateral rate under a variable cross section: one kernel
# pass for the history rates, the lateral rows from the column sums, then a
# second pass at them. `du` and the workspace hold the last pass.
function _consistent_lateral!(du, u, def::ResidualDefinition, ws::ResidualWorkspace, t)
    layout = def.layout
    residual!(du, u, def, t, ws)
    ΣG, ΣM, ΣĠ, ΣṀ = 0.0, 0.0, 0.0, 0.0
    for i = 1:n_cells(layout)
        kG, kM = history_index(layout, i, 1), history_index(layout, i, 2)
        ΣG += u[kG]
        ΣM += u[kM]
        ΣĠ += du[kG]
        ΣṀ += du[kM]
    end
    A, θ_A = consistent_lateral(def.closure.lateral, ΣG, ΣM, ΣĠ, ΣṀ)
    u[area_index(layout)] = A
    u[lateral_rate_index(layout)] = θ_A
    residual!(du, u, def, t, ws)
    return du
end

# The algebraic unknowns of `u` that satisfy their rows at time `t`, given its
# differential unknowns, in the order the rows read each other: under a
# variable cross section its rows first, twice, since the radiation rows'
# optical depths and the history rates the lateral rate reads both read the
# concentrations at the area; then the radiation rows; then one kernel pass
# for the strain rates and the closure's unknowns. Returns `nothing`, or the
# residual and scale of a closure row left above round-off.
function _consistent_algebraic!(u, def::ResidualDefinition, ws::ResidualWorkspace, t)
    cl = def.closure
    du = similar(u)
    if has_chi(def.layout)
        _consistent_lateral!(du, u, def, ws, t)
        _consistent_lateral!(du, u, def, ws, t)
    end
    _consistent_radiation!(u, def, ws, t)
    cl === nothing && return nothing
    residual!(du, u, def, t, ws)
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
rows at time `t`, given its differential unknowns: under a variable cross
section the area and its lateral rate first, since the radiation rows'
optical depths and the history rates read the concentrations at the area;
then the face fluxes of radiative transfer; then the node velocities or the
cells' multipliers by the mesh closure. The
integrator starts from such a state; a closure row left above round-off is
an error.
"""
function consistent_state!(
    u::AbstractVector,
    def::ResidualDefinition,
    ws::ResidualWorkspace,
    t,
)
    left = _consistent_algebraic!(u, def, ws, t)
    left === nothing || error(
        "the closure's consistent unknowns leave a row at $(left[1]) (scale $(left[2]))",
    )
    return u
end

# The face fluxes of `u` that satisfy the radiative transfer rows at time `t`.
function _consistent_radiation!(u, def::ResidualDefinition, ws::ResidualWorkspace, t)
    layout = def.layout
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    kernel === nothing && return u
    ds = ws.derived
    derive!(ds, u, layout, def.setup)
    consistent_radiation!(ds.q, kernel, ds, def.setup.topology, t)
    for f in axes(ds.q, 2), k in axes(ds.q, 1)
        u[radiation_index(layout, f, k)] = ds.q[k, f]
    end
    return u
end

# The array of the slots' initial cells: the one a bound correction reads.
_slot_origin(cl::BalanceRowClosure{<:LayerCorrection}, n) = cl.correction.origin
_slot_origin(cl, n) = collect(1:n)

_register_kernel!(t::SlotTables, kernel) = t
_register_kernel!(t::SlotTables, k::ReactionKernel) = register_lineage!(t, k.reference)
_register_kernel!(
    t::SlotTables,
    k::Union{ConductionKernel,GasTransportKernel,ThermalBoundaryKernel},
) = register_tables!(t, :face, k.d_floor)
_register_kernel!(t::SlotTables, k::RadiativeTransferKernel) =
    register_tables!(t, :matrix, k.d_floor)
function _register_kernel!(t::SlotTables, k::RelativeFluxKernel)
    register_tables!(t, :matrix, k.ε)
    register_tables!(t, :cell, k.η, k.mass_scale)
    return register_tables!(t, :sides, k.d_floor)
end

_register_closure!(t::SlotTables, closure) = t
function _register_closure!(t::SlotTables, cl::BalanceRowClosure{<:LayerCorrection})
    corr = cl.correction
    push!(t.gains, () -> slot_gains!(corr))
    return t
end

# Every table of the solve that is read by slot. The caloric extension reads
# the mass scale and capacity floor of `reference`, the same arrays.
function _register_tables!(
    t::SlotTables,
    cells::CellDeactivation,
    reference::CellReference,
    def::ResidualDefinition,
    tol::Vector{Float64};
    scale::Float64,
    time_scale::Float64,
)
    register_tables!(
        t,
        :cell,
        reference.volume,
        reference.solid_mass,
        reference.mass_scale,
        reference.volume_floor,
        reference.capacity_floor,
        reference.thickness,
        cells.reference,
    )
    isempty(cells.reactive) || register_tables!(t, :cell, cells.reactive)
    register_tables!(t, :node, reference.positions)
    register_rows!(t, tol, reference.thickness; scale, time_scale)
    foreach(k -> _register_kernel!(t, k), def.kernels)
    return _register_closure!(t, def.closure)
end

"""
	setup(problem; energy_state = :temperature, reltol = 1e-6, abstol = 1e-8,
	      threshold = nothing, pressure_tolerance = 1.0, use_ale = false,
	      radiation_model = SURFACE_ABSORPTION, correction = nothing, lagrangian = nothing,
	      deactivation = false, band = nothing, closure = nothing, slip = (),
	      tensor = nothing, history = false, deformation = false,
	      time_scale = problem.dt_initial, sources = NoSources()) -> PreparedProblem

Set up `problem`, whose mesh is a `Mesh{D}` ([`uniform_mesh`](@ref)
and its siblings on a line, [`box_mesh`](@ref) and [`cell_mesh`](@ref) in two
and three dimensions): the initial inventories
`mᵢⱼ = ξⱼ(xᵢ) Vᵢ` and the energy state from `T(xᵢ)` at the cell centres `xᵢ`
([`initial_fields`](@ref)), the gas transport the material and the mass
conditions need, the kernels, the radiation rows consistent with the rest of
the state, the tolerance vector and the colored Jacobian. The coordinate system and the cross-section measure are the mesh's;
the thermal and mass conditions are keyed by its boundary tags.

`radiation_model` is `SURFACE_ABSORPTION` (the incident flux enters the surface
balance at the exposed cell's absorptivity, or the condition's own) or
`RADIATIVE_TRANSFER` (every cell absorbs and
emits, through the face fluxes of [`RadiativeTransferKernel`](@ref); the
`RadiativeBC` and `RadiativeFluxBC` members of the `:top` and `:bottom`
conditions are the radiation arriving from outside, and leave the surface
balance; a `RadiativeFluxBC` on `:top` arrives diffuse or as the beams of its
`distribution`).

`use_ale` moves the mesh with the solid: the mesh's pinned nodes (on a line
the back node) stay at rest and every other node follows the cells' reaction
strain rates through the mesh closure, with its unknowns consistent at the
start; without it every node is pinned. On a line the closure is the
balance rows of [`BalanceRowClosure`](@ref), and a material with a lateral
shrinkage law also carries the conversion-history pair and its cross
section's area and lateral rate. `correction` moves a line's mesh relative
to the solid ([`EquidistributionCorrection`](@ref); [`NoCorrection`](@ref)
keeps it with the solid), the balance rows then holding the solid velocity
`w − c`, and the inventories cross the faces at the relative velocity `−c`
([`RelativeFluxKernel`](@ref)). By default a solve gets the equidistribution
correction when the mesh moves, cells are deactivated and deactivation can
remove some cell (its components, converted completely, leave less condensed
mass than `threshold` of it, [`residues`](@ref)), and no correction
otherwise; a layer deactivation cannot remove keeps its mesh with the solid
under any equidistribution correction. `lagrangian` lists the mesh
nodes that always move with the solid, the layer interfaces, beside the
mesh's own; under an equidistribution correction every interior node between
cells that start with different condensed components must be among them, and
by default they are those nodes. A [`PrescribedCorrection`](@ref) prescribes the
velocity of every node instead, moving the mesh over a solid at rest
(verification; the material must not react), the inventories crossing the faces
at `−w`. `time_scale` [s] scales the velocity tolerances (the first step).

A reaction of order other than one is evaluated on the conversion of its
reactant against its lineage inventory, a constant of each cell
([`lineage_references`](@ref)): the cells keep their material. Such a reaction
is refused on a mesh of more than one dimension and with a correction, which
would mix material at different conversions.

`deactivation` deactivates, during the solve, every cell whose condensed mass
falls below the fraction `threshold` (`ε_m`, 1e-2 by default) of its initial
condensed mass, or, at a reaction order other than one, whose reactive inventory
falls below that fraction of its reference ([`below_threshold`](@ref); the default
is then a tenth of the limiter's `conversion_threshold`, 1e-4 with the default
limiter) ([`CellDeactivation`](@ref)): on a moving
mesh wherever it sits, together with the contiguous cells around it below the
fraction `band` (`ε_band`, in `[threshold, 1)`; by default the larger of the
threshold and [`CORRECTION_BAND`](@ref) when the bound equidistribution correction
moves some layer, otherwise the threshold, so that only the cells below it go),
the live cells above moving down into their slots; on a fixed mesh from the exposed
`:top` side only ([`deactivate!`](@ref)).
`threshold` also scales the mass tolerances; `sources` is the external-source
hook.

On a mesh of quadrilaterals or hexahedra the closure is `:potential` (the
default), one multiplier per cell with the node velocities of minimum
lumped-volume norm that satisfy every cell's volume balance
([`PotentialClosure`](@ref)), the pinned nodes at rest and on each `tag => k`
of `slip` the component `k` zero; or `:balance`, the balance rows of a single
column of cells. `tensor` is each cell's dimensionless conductivity tensor,
an `SMatrix{D,D}` for every cell or a function of the centroid returning one
(the identity by default): the conductivity is `k` times it. `history`
carries the conversion-history pair (`G = 0`, `M` the cell's initial solid
mass) and `deformation` the deformation gradient (the identity at the
start) on a moving mesh. Gas transport, in-depth radiation, the variable
cross section, the equidistribution correction and deactivation are a
line's, and `closure`, `slip`, `tensor`, `history` and `deformation` are the
other meshes': a problem that names one where it does not apply is refused.
"""
function setup(
    problem::PyrolysisProblem{<:Mesh{D}};
    energy_state::Symbol = :temperature,
    reltol::Real = 1e-6,
    abstol::Real = 1e-8,
    threshold::Union{Nothing,Real} = nothing,
    pressure_tolerance::Real = 1.0,
    use_ale::Bool = false,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    correction::Union{Nothing,Correction} = nothing,
    lagrangian::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    deactivation::Bool = false,
    band::Union{Nothing,Real} = nothing,
    closure::Union{Nothing,Symbol} = nothing,
    slip = (),
    tensor = nothing,
    history::Bool = false,
    deformation::Bool = false,
    time_scale::Real = problem.dt_initial,
    sources = NoSources(),
) where {D}
    material = problem.material
    NC = n_components(material)
    mesh = problem.mesh
    threshold = something(threshold, _default_threshold(material, deactivation))
    T0, ξ0 = initial_fields(problem)
    isfinite(pressure_tolerance) && pressure_tolerance>0 ||
        throw(ArgumentError("pressure_tolerance must be finite and positive"))
    errors=Problem.validate_initial_pore_state(problem)
    isempty(errors) || throw(ArgumentError(join(errors, "\n")))
    for bc in values(problem.bc_set.mass)
        validate_gas_boundary(bc, material, problem.tspan[1])
    end
    correction = something(
        correction,
        _installed_correction(material, ξ0, use_ale, deactivation, threshold),
    )
    validate_transport_support(
        material,
        ξ0;
        moving_support = !(correction isa NoCorrection),
    )
    closure = something(closure, D == 1 ? :balance : :potential)
    _require_supported(
        problem,
        use_ale,
        radiation_model,
        correction,
        deactivation;
        closure,
        slip,
        tensor,
        history,
        deformation,
    )
    lagrangian = something(
        lagrangian,
        correction isa EquidistributionCorrection ?
        _composition_interfaces(ξ0, material) : Int[],
    )
    topo = _motion(mesh.topology, use_ale, lagrangian)
    n = n_cells(topo)
    x = reduce(vcat, mesh.positions)
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, x, topo.cross_section, topo)
    m0 = [ξ0[j, i] * g0.cell_volume[i] for j = 1:NC, i = 1:n]
    _require_layers(correction, m0, material, topo.lagrangian)
    condensed0 = [
        sum(m0[j, i] for j = 1:NC if !is_gas(material.components[j]); init = 0.0) for
        i = 1:n
    ]

    reference = CellReference(
        Mesh(topo, mesh.positions),
        m0,
        material;
        reltol = Float64(reltol),
        threshold = Float64(threshold),
        pressure_tolerance = Float64(pressure_tolerance),
    )
    energy = _energy_state(energy_state)
    laws = CaloricLaws(material)
    caloric =
        energy isa EnthalpyState ?
        CaloricExtension(laws, reference.mass_scale, reference.capacity_floor) : nothing
    gas = _gas_transport(problem)
    references = _reaction_references(material, m0, condensed0, sources)
    state_setup = StateSetup(topo, material, laws, reference, caloric, gas)
    chi = material.lateral_shrinkage_law !== nothing
    history = history || chi
    potential = use_ale && closure === :potential && !(correction isa PrescribedCorrection)
    layout = StateLayout(
        NC,
        topo,
        energy;
        moving = use_ale,
        radiation_rows = _radiation_rows(radiation_model, problem.bc_set),
        chi,
        history,
        deformation,
        unknowns = potential ? :multipliers : :velocities,
    )

    M0 = [
        sum(m0[j, i] for j = 1:NC if is_solid(material.components[j]); init = 0.0) for
        i = 1:n
    ]
    lateral = chi ? LateralShrinkage(material, topo.cross_section, sum(M0)) : nothing
    cl =
        !use_ale ? nothing :
        potential ? PotentialClosure(topo; slip) :
        BalanceRowClosure(
            topo;
            correction = bind_correction(
                correction,
                topo,
                reference,
                m0,
                material,
                threshold,
            ),
            lateral,
        )
    cl === nothing ||
        rank_check(cl, topo, g0) ||
        throw(ArgumentError("the closure's rows do not determine the mesh's motion"))
    relative =
        correction isa NoCorrection ? () :
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
    band = something(band, _moves(cl) ? max(threshold, CORRECTION_BAND) : threshold)
    reactants, reactive = _reactive_references(material, references, m0, T0)
    cells =
        deactivation ?
        CellDeactivation(
            threshold,
            band,
            condensed0,
            map(c -> !is_gas(c), material.components),
            something(
                caloric,
                CaloricExtension(laws, reference.mass_scale, reference.capacity_floor),
            );
            origin = _slot_origin(cl, n),
            reactants,
            reactive,
            lineage = references,
        ) : nothing
    tensors = D == 1 ? nothing : _tensors(tensor, g0, n, Val(D))
    kernels = _kernels(
        problem,
        topo,
        g0,
        laws,
        gas,
        radiation_model,
        relative,
        tensors,
        pressure_tolerance,
        references,
    )
    def = ResidualDefinition(layout, state_setup, kernels, cl, sources)
    ws = ResidualWorkspace(def, T0)

    u0 = zeros(state_length(layout))
    for i = 1:n
        mi = ntuple(j -> m0[j, i], NC)
        for j = 1:NC
            u0[mass_index(layout, i, j)] = mi[j]
        end
        u0[energy_index(layout, i)] = energy_state_value(energy, state_setup, T0[i], mi, i)
        history && (u0[history_index(layout, i, 2)] = M0[i])
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
        radiation_scale = _transfer_scale(kernels, ws.derived, problem.tspan),
        time_scale = Float64(time_scale),
    )
    cells === nothing || _register_tables!(
        cells.tables,
        cells,
        reference,
        def,
        tol;
        scale = reltol * threshold,
        time_scale = Float64(time_scale),
    )
    cj = ColoredJacobianCache(def, ws)
    cp=PreparedProblem(problem, def, ws, u0, Float64(reltol), tol, cj, T0, m0, cells)
    derive!(ws.derived, u0, layout, state_setup)
    _report_support(find_kernel(ConductionKernel, kernels), ws.derived, topo)
    check=DomainCheck(cp)
    _porosity_out(check, ws.derived, gas) && throw(ArgumentError(check.monitor.message))
    return cp
end

# The live cells whose conduction gradient is a regularized fit at the
# initial state.
function _report_support(conduction::ConductionKernel, ds, topo)
    s = gradient_support(conduction, ds, topo)
    s.cells > 0 &&
        @info "setup(): the conduction gradients of $(s.cells) cells are regularized fits " *
              "(smallest normalized determinant $(s.worst))"
    return nothing
end

# Every cell's dimensionless conductivity tensor.
_tensors(::Nothing, g0, n, ::Val{D}) where {D} = fill(SMatrix{D,D,Float64}(I), n)
_tensors(A::SMatrix, g0, n, ::Val{D}) where {D} = fill(SMatrix{D,D,Float64}(A), n)
_tensors(A::AbstractVector, g0, n, ::Val{D}) where {D} =
    [SMatrix{D,D,Float64}(a) for a in A]
_tensors(A::Function, g0, n, ::Val{D}) where {D} =
    [SMatrix{D,D,Float64}(A(g0.cell_center[i])) for i = 1:n]
