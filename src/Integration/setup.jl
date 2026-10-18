# ============================================================================
# From a PyrolysisProblem to the core's definition, workspace and state
# ============================================================================

"""
	CoreProblem

A `PyrolysisProblem` set up on the finite-volume core: the problem (with its
own copy of the mesh), the definition and workspace the residual reads, the
initial state, the absolute-tolerance vector, the colored Jacobian, and the
initial temperatures and masses.
"""
struct CoreProblem{P,D<:CoreDefinition,W<:CoreWorkspace,J<:ColoredJacobianCache}
    problem::P
    def::D
    ws::W
    u0::Vector{Float64}
    abstol::Vector{Float64}
    jacobian::J
    T0::Vector{Float64}
    m0::Matrix{Float64}
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
function _require_supported(problem::PyrolysisProblem, use_ale::Bool, radiation_model)
    material = problem.material
    use_ale && throw(ArgumentError("the finite-volume core does not move the mesh yet"))
    radiation_model === NO_RADIATION ||
        throw(ArgumentError("the finite-volume core does not model in-depth radiation yet"))
    material.lateral_shrinkage_law === nothing ||
        throw(ArgumentError("the finite-volume core has no variable cross section yet"))
    has_darcy_flow(material) &&
        throw(ArgumentError("the finite-volume core has no Darcy transport yet"))
    for comp in material.components
        λ = comp.gas_transfer_coeff
        all(T -> λ(T) == 0, (T_PROPERTY_MIN, 300.0, 1000.0, T_PROPERTY_MAX)) || throw(
            ArgumentError(
                "the finite-volume core has no gas transport yet (component :$(comp.name) has λ ≠ 0)",
            ),
        )
    end
    for (tag, bc) in pairs(problem.bc_set.mass)
        bc isa ImpermeableBC || throw(
            ArgumentError(
                "the finite-volume core has no mass boundary flux yet (tag :$tag)",
            ),
        )
    end
    return nothing
end

"""
	setup(problem; energy_state = :temperature, reltol = 1e-6, abstol = 1e-8,
	      threshold = 1e-2, use_ale = false, radiation_model = NO_RADIATION,
	      sources = NoSources()) -> CoreProblem

Set up `problem` on the finite-volume core: apply its initial conditions to its
own copy of the mesh, build the 1D mesh from the node positions and the
cross-section area, the initial inventories `mᵢⱼ = ξⱼ(zᵢ) Vᵢ`, the energy
state from `T(zᵢ)`, the kernels, the tolerance vector and the colored Jacobian.

`threshold` is the deactivation mass fraction `ε_m`, which also scales the
mass tolerances; `sources` is the external-source hook.
"""
function setup(
    problem::PyrolysisProblem;
    energy_state::Symbol = :temperature,
    reltol::Real = 1e-6,
    abstol::Real = 1e-8,
    threshold::Real = 1e-2,
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    sources = NoSources(),
)
    _require_supported(problem, use_ale, radiation_model)
    problem = _own_problem(problem)
    apply_initial_conditions!(problem)
    material = problem.material
    NC = n_components(material)
    old = problem.mesh

    x = [node.z for node in old.nodes]
    mesh = line_mesh(
        x;
        cross_section = old.cross_section_area_initial,
        bottom = :bottom,
        top = :top,
        fixed = !use_ale,
    )
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
    state_setup = StateSetup(topo, material, laws, reference, caloric)
    layout = StateLayout(NC, n, energy; moving = use_ale)

    kernels = (
        (n_reactions(material) > 0 ? (ReactionKernel(material),) : ())...,
        ConductionKernel(topo, g0),
        ThermalBoundaryKernel(topo, g0, problem.bc_set),
    )
    def = CoreDefinition(layout, state_setup, kernels, sources)
    ws = CoreWorkspace(def, T0)

    u0 = zeros(state_length(layout))
    for i = 1:n
        mi = ntuple(j -> m0[j, i], NC)
        for j = 1:NC
            u0[mass_index(layout, i, j)] = mi[j]
        end
        u0[energy_index(layout, i)] = energy_state_value(energy, state_setup, T0[i], mi, i)
    end
    tol = abstol_vector(
        layout,
        state_setup,
        m0,
        T0;
        reltol = Float64(reltol),
        abstol = Float64(abstol),
        threshold = Float64(threshold),
    )
    cj = ColoredJacobianCache(def, ws, u0, problem.tspan[1])
    return CoreProblem(problem, def, ws, u0, tol, cj, T0, m0)
end
