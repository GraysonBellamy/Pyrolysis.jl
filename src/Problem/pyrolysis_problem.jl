"""
PyrolysisProblem and PyrolysisSolution type definitions.
"""

# ============================================================================
# PyrolysisProblem Type
# ============================================================================

"""
	PyrolysisProblem{M,Mat,BCS,FT,Fξ} <: AbstractPyrolysisProblem

Complete specification of a pyrolysis simulation problem.

# Type Parameters
- `M`: Mesh type
- `Mat`: Material type
- `BCS`: BoundaryConditionSet type (type-stable BC storage)
- `FT`: Initial temperature function type
- `Fξ`: Initial concentration functions type

# Fields
- `mesh::M`: The computational mesh
- `material::Mat`: Material definition with components and reactions
- `bc_set::BCS`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial::FT`: Initial temperature function z -> T(z)
- `ξ_initial::Vector{Fξ}`: Initial concentration functions [z -> ξⱼ(z)]
- `tspan::Tuple{Float64,Float64}`: Simulation time span (t_start, t_end)
- `dt_initial::Float64`: Initial timestep [s]
- `output_times::Vector{Float64}`: Times at which to save output

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
	mass    = (top = ConvectiveMassBC(h_m = 0.01), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
	mesh = create_uniform_mesh(0.01, 100, 2),
	material = pmma,
	bc_set = bc_set,
	T_initial = z -> 300.0,
	ξ_initial = [z -> 1190.0, z -> 0.0],
	tspan = (0.0, 100.0)
)
```
"""
struct PyrolysisProblem{M,Mat,BCS,FT,Fξ} <: AbstractPyrolysisProblem
    mesh::M
    material::Mat
    bc_set::BCS
    T_initial::FT
    ξ_initial::Vector{Fξ}
    tspan::Tuple{Float64,Float64}
    dt_initial::Float64
    output_times::Vector{Float64}
end

# ============================================================================
# Constructors
# ============================================================================

"""
	PyrolysisProblem(; mesh, material, bc_set::BoundaryConditionSet,
					 T_initial=z->300.0, ξ_initial=nothing,
					 tspan, dt_initial=1e-4, output_times=Float64[])

Create a PyrolysisProblem with BoundaryConditionSet (preferred, type-stable).

Radiative heating, convective coupling, and ambient temperatures are carried
exclusively by `bc_set` (`RadiativeFluxBC`, `ConvectiveBC`, `RadiativeBC`).

# Arguments
- `mesh`: The computational mesh
- `material`: Material definition
- `bc_set`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial=z->300.0`: Initial temperature function
- `ξ_initial=nothing`: Initial concentration functions (defaults from material)
- `tspan`: Simulation time span (required)
- `dt_initial=1e-4`: Initial timestep
- `output_times=Float64[]`: Times to save output (auto-generated if empty)

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
	           bottom = AdiabaticBC()),
	mass    = (top = DiffusiveMassBC(Inf), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
	mesh = mesh,
	material = pmma,
	bc_set = bc_set,
	T_initial = z -> 300.0,
	tspan = (0.0, 100.0)
)
```
"""
function PyrolysisProblem(;
    mesh,
    material,
    bc_set::BoundaryConditionSet,
    T_initial = z -> 300.0,
    ξ_initial = nothing,
    tspan::Tuple{Real,Real},
    dt_initial::Real = 1e-4,
    output_times::Vector{<:Real} = Float64[],
)
    NC = n_components(material)

    # Create default initial concentration functions if not provided
    if ξ_initial === nothing
        # Default: first solid component at full density, others at zero
        ξ_funcs = Vector{Function}(undef, NC)
        for j = 1:NC
            comp = material.components[j]
            if j == 1 && is_solid(comp)
                # First solid at reference density
                ρ_ref = comp.density(300.0)
                ξ_funcs[j] = z -> ρ_ref
            else
                ξ_funcs[j] = z -> 0.0
            end
        end
        ξ_initial = ξ_funcs
    end

    # Convert ξ_initial to vector of functions if needed
    if ξ_initial isa Function
        ξ_initial = [ξ_initial for _ = 1:NC]
    end

    # Auto-generate output times if empty
    output_times_float = Float64.(output_times)
    if isempty(output_times_float)
        t_start, t_end = Float64.(tspan)
        n_outputs = min(101, max(11, round(Int, (t_end - t_start) / dt_initial / 100)))
        output_times_float = collect(range(t_start, t_end, length = n_outputs))
    end

    PyrolysisProblem(
        mesh,
        material,
        bc_set,
        T_initial,
        ξ_initial,
        (Float64(tspan[1]), Float64(tspan[2])),
        Float64(dt_initial),
        output_times_float,
    )
end

"""
	PyrolysisProblem(mesh, material, bc_set::BoundaryConditionSet, tspan; kwargs...)

Simplified constructor for common use cases with type-stable BoundaryConditionSet.

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
)
problem = PyrolysisProblem(mesh, pmma, bc_set, (0.0, 100.0))
```
"""
function PyrolysisProblem(
    mesh,
    material,
    bc_set::BoundaryConditionSet,
    tspan::Tuple;
    kwargs...,
)
    PyrolysisProblem(;
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        tspan = tspan,
        kwargs...,
    )
end

# ============================================================================
# PyrolysisSolution Type
# ============================================================================

"""
	PyrolysisSolution{P, ExtraT} <: AbstractPyrolysisSolution

Container for pyrolysis simulation results, organized as a fixed core schema
of concretely-typed fields plus a typed `extras` NamedTuple for opt-in
diagnostics.

# Core fields
- `problem::P`: Reference to the original problem
- `t::Vector{Float64}`: Solution times [s]
- `T::Matrix{Float64}`: Temperature profiles, shape `(n_cells, n_steps)` [K]
- `ξ::Array{Float64,3}`: Concentration profiles, shape `(NC, n_cells, n_steps)` [kg/m³]
- `z::Matrix{Float64}`: Node positions, shape `(n_nodes, n_steps)` [m]; an
  empty `0×0` matrix in Eulerian mode
- `surface_T::Vector{Float64}`: Surface temperature at each saved time [K]
- `mass_loss_rate::Vector{Float64}`: Surface MLR per current area `A(t)` [kg/(m²·s)]
- `energy_residual::Vector{Float64}`: Energy-conservation residual; zeros when
  the conservation tracker is disabled
- `retcode::Symbol`: Return code from ODE solver
- `solve_time::Float64`: Wall-clock solve time [s]
- `raw_solution`: Raw DifferentialEquations.jl solution object (interpolation source)
- `extras::ExtraT`: Typed NamedTuple of opt-in diagnostics (see below)

# Extras
`extras` aggregates secondary diagnostics keyed by symbolic name. Vectors
inside `extras` are mutable, so the per-step extraction code and the
diagnostics callback write through the NamedTuple in place. Names:
`gas_generation_rate`, `back_temperature`, `total_mass`, `thickness_history`,
`surface_mass_flux`, `surface_heat_flux_{net,absorbed,reradiation,convection,
conduction,transpiration}`, `surface_emissivity`, `surface_conductivity`,
`back_heat_flux`, `MLR_gauge`, `cross_section_area_history`, `diagnostics_log`.

# Accessors
```julia
state_t   = sol(t)            # interpolate raw solution at time t
T_profile = sol.T[:, i]       # temperature profile at saved time i
T_surface = sol.surface_T[i]  # surface temperature at saved time i
mlr       = sol.mass_loss_rate[i]
mlr_gauge = sol.extras.MLR_gauge[i]
```
"""
mutable struct PyrolysisSolution{P<:AbstractPyrolysisProblem,ExtraT<:NamedTuple} <:
               AbstractPyrolysisSolution
    problem::P
    t::Vector{Float64}
    T::Matrix{Float64}
    ξ::Array{Float64,3}
    z::Matrix{Float64}
    surface_T::Vector{Float64}
    mass_loss_rate::Vector{Float64}
    energy_residual::Vector{Float64}
    retcode::Symbol
    solve_time::Float64
    # Raw integrator solution; left untyped so reassignment from `nothing`
    # (set at construction) to the SciMLBase.ODESolution returned by solve()
    # does not collide with a concretized type parameter.
    raw_solution::Any
    extras::ExtraT
end

"""
	_empty_extras() -> NamedTuple

Build the canonical `extras` NamedTuple of empty mutable vectors. This pins
the extras schema for the fixed-shape mutable struct and lets the diagnostics
callback and extraction code both write through `sol.extras.x` without
re-typing the parent solution.
"""
function _empty_extras()
    return (
        gas_generation_rate = Float64[],
        back_temperature = Float64[],
        total_mass = Float64[],
        thickness_history = Float64[],
        surface_mass_flux = Float64[],
        surface_heat_flux_net = Float64[],
        surface_heat_flux_absorbed = Float64[],
        surface_heat_flux_reradiation = Float64[],
        surface_heat_flux_convection = Float64[],
        surface_heat_flux_conduction = Float64[],
        surface_heat_flux_transpiration = Float64[],
        volumetric_absorption = Float64[],
        radiation_transmitted = Float64[],
        surface_emissivity = Float64[],
        surface_conductivity = Float64[],
        back_heat_flux = Float64[],
        MLR_gauge = Float64[],
        cross_section_area_history = Float64[],
        live_cells = Int[],
        park_times = Float64[],
        diagnostics_log = Any[],
    )
end

"""
	PyrolysisSolution(problem)

Create an empty solution container for `problem`. The core arrays are
zero-sized until `extract_solution!` (or its ALE sibling) populates them, and
the `extras` NamedTuple of vectors is ready to receive in-place writes.
"""
function PyrolysisSolution(problem::PyrolysisProblem)
    PyrolysisSolution(
        problem,
        Float64[],                            # t
        Matrix{Float64}(undef, 0, 0),         # T
        Array{Float64,3}(undef, 0, 0, 0),    # ξ
        Matrix{Float64}(undef, 0, 0),         # z
        Float64[],                            # surface_T
        Float64[],                            # mass_loss_rate
        Float64[],                            # energy_residual
        :Default,                             # retcode
        0.0,                                  # solve_time
        nothing,                              # raw_solution
        _empty_extras(),                      # extras
    )
end

"""
	_allocate_extras!(sol::PyrolysisSolution, n_times::Int) -> sol

Replace `sol.extras` with `n_times`-sized zero vectors for every numeric
extras field, preserving any `diagnostics_log` rows the callback may have
already pushed. Type-stable: the resulting NamedTuple matches the schema
pinned at construction by `_empty_extras`.
"""
function _allocate_extras!(sol::PyrolysisSolution, n_times::Int)
    prev_log = sol.extras.diagnostics_log
    sol.extras = (
        gas_generation_rate = zeros(n_times),
        back_temperature = zeros(n_times),
        total_mass = zeros(n_times),
        thickness_history = zeros(n_times),
        surface_mass_flux = zeros(n_times),
        surface_heat_flux_net = zeros(n_times),
        surface_heat_flux_absorbed = zeros(n_times),
        surface_heat_flux_reradiation = zeros(n_times),
        surface_heat_flux_convection = zeros(n_times),
        surface_heat_flux_conduction = zeros(n_times),
        surface_heat_flux_transpiration = zeros(n_times),
        volumetric_absorption = zeros(n_times),
        radiation_transmitted = zeros(n_times),
        surface_emissivity = zeros(n_times),
        surface_conductivity = zeros(n_times),
        back_heat_flux = zeros(n_times),
        MLR_gauge = zeros(n_times),
        cross_section_area_history = zeros(n_times),
        live_cells = zeros(Int, n_times),
        park_times = Float64[],
        diagnostics_log = prev_log,
    )
    return sol
end

# ============================================================================
# Accessors for PyrolysisProblem
# ============================================================================

"""
	get_mesh(prob::PyrolysisProblem)

Get the mesh from a problem.
"""
get_mesh(prob::PyrolysisProblem) = prob.mesh

"""
	get_material(prob::PyrolysisProblem)

Get the material from a problem.
"""
get_material(prob::PyrolysisProblem) = prob.material

"""
	n_components(prob::PyrolysisProblem)

Get number of components in the problem's material.
"""
n_components(prob::PyrolysisProblem) = n_components(prob.material)

"""
	tspan(prob::PyrolysisProblem)

Get the time span of the problem.
"""
tspan(prob::PyrolysisProblem) = prob.tspan

# ============================================================================
# Accessors for PyrolysisSolution
# ============================================================================

"""
	get_times(sol::PyrolysisSolution)

Get the solution times.
"""
get_times(sol::PyrolysisSolution) = sol.t

"""
	get_temperatures(sol::PyrolysisSolution)

Get the temperature profile matrix `(n_cells, n_steps)`.
"""
get_temperatures(sol::PyrolysisSolution) = sol.T

"""
	get_concentrations(sol::PyrolysisSolution)

Get the concentration array `(NC, n_cells, n_steps)`.
"""
get_concentrations(sol::PyrolysisSolution) = sol.ξ

"""
	get_gas_generation_rate(sol::PyrolysisSolution)

Get volumetric gas generation rate history (from `extras`).

Note: This is NOT the actual mass loss rate (surface flux). The canonical MLR
lives in `sol.mass_loss_rate`.
"""
get_gas_generation_rate(sol::PyrolysisSolution) = sol.extras.gas_generation_rate

"""
	n_timesteps(sol::PyrolysisSolution)

Get number of saved timesteps.
"""
n_timesteps(sol::PyrolysisSolution) = length(sol.t)

"""
	succeeded(sol::PyrolysisSolution)

Check if simulation completed successfully.
"""
succeeded(sol::PyrolysisSolution) = sol.retcode == :Success

# ============================================================================
# Interpolation Interface
# ============================================================================

"""
	(sol::PyrolysisSolution)(t::Real)

Interpolate solution at time t. Returns the state vector.
"""
function (sol::PyrolysisSolution)(t::Real)
    if sol.raw_solution !== nothing
        return sol.raw_solution(t)
    else
        error("No interpolatable solution available. Use saved times instead.")
    end
end

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, prob::PyrolysisProblem)
    t_start, t_end = prob.tspan
    n_cells = length(prob.mesh.cells)
    NC = n_components(prob.material)
    print(
        io,
        "PyrolysisProblem($(prob.material.name), $n_cells cells, $NC components, t∈[$t_start, $t_end])",
    )
end

function Base.show(io::IO, ::MIME"text/plain", prob::PyrolysisProblem)
    t_start, t_end = prob.tspan
    n_cells = length(prob.mesh.cells)
    NC = n_components(prob.material)
    NR = n_reactions(prob.material)

    println(io, "PyrolysisProblem")
    println(io, "  Material: $(prob.material.name)")
    println(io, "  Components: $NC")
    println(io, "  Reactions: $NR")
    println(io, "  Mesh: $n_cells cells")
    println(io, "  Time span: [$t_start, $t_end] s")
    println(io, "  Initial dt: $(prob.dt_initial) s")
    println(io, "  Output times: $(length(prob.output_times))")
    println(
        io,
        "  Thermal BCs: top=$(typeof(prob.bc_set.thermal.top)), bottom=$(typeof(prob.bc_set.thermal.bottom))",
    )
    println(
        io,
        "  Mass BCs: top=$(typeof(prob.bc_set.mass.top)), bottom=$(typeof(prob.bc_set.mass.bottom))",
    )
end

function Base.show(io::IO, sol::PyrolysisSolution)
    n_times = length(sol.t)
    print(io, "PyrolysisSolution($(sol.retcode), $n_times timesteps, $(sol.solve_time)s)")
end

function Base.show(io::IO, ::MIME"text/plain", sol::PyrolysisSolution)
    println(io, "PyrolysisSolution")
    println(io, "  Return code: $(sol.retcode)")
    println(io, "  Timesteps saved: $(length(sol.t))")
    if !isempty(sol.t)
        println(io, "  Time range: [$(sol.t[1]), $(sol.t[end])] s")
    end
    println(io, "  Solve time: $(round(sol.solve_time, digits=3)) s")
    if !isempty(sol.extras.gas_generation_rate)
        max_gen = maximum(sol.extras.gas_generation_rate)
        println(io, "  Peak gas generation: $(round(max_gen * 1000, digits=2)) g/(m²·s)")
    end
    if !isempty(sol.surface_T)
        max_T = maximum(sol.surface_T)
        println(io, "  Peak surface T: $(round(max_T, digits=1)) K")
    end
end
