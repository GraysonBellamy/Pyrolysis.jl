"""
PyrolysisProblem and PyrolysisSolution type definitions.

Reference: Design Document Section 6.7, 12
"""

# ============================================================================
# PyrolysisProblem Type
# ============================================================================

"""
	PyrolysisProblem{M,Mat,TBC,MBC,FT,Fξ,Fq,FTamb} <: AbstractPyrolysisProblem

Complete specification of a pyrolysis simulation problem.

# Type Parameters
- `M`: Mesh type
- `Mat`: Material type
- `BCS`: BoundaryConditionSet type (type-stable BC storage)
- `FT`: Initial temperature function type
- `Fξ`: Initial concentration functions type
- `Fq`: External flux function type
- `FTamb`: Ambient temperature function type

# Fields
- `mesh::M`: The computational mesh
- `material::Mat`: Material definition with components and reactions
- `bc_set::BCS`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial::FT`: Initial temperature function z -> T(z)
- `ξ_initial::Vector{Fξ}`: Initial concentration functions [z -> ξⱼ(z)]
- `external_flux::Fq`: External heat flux function t -> q"(t) [W/m²]
- `ambient_temperature::FTamb`: Ambient temperature function t -> T∞(t) [K]
- `tspan::Tuple{Float64,Float64}`: Simulation time span (t_start, t_end)
- `dt_initial::Float64`: Initial timestep [s]
- `output_times::Vector{Float64}`: Times at which to save output

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3),
	thermal_bottom = AdiabaticBC(),
	mass_top = ConvectiveMassBC(h_m=0.01),
	mass_bottom = ImpermeableBC()
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
struct PyrolysisProblem{M, Mat, BCS, FT, Fξ, Fq, FTamb} <: AbstractPyrolysisProblem
	mesh::M
	material::Mat
	bc_set::BCS
	T_initial::FT
	ξ_initial::Vector{Fξ}
	external_flux::Fq
	ambient_temperature::FTamb
	tspan::Tuple{Float64, Float64}
	dt_initial::Float64
	output_times::Vector{Float64}
end

# ============================================================================
# Constructors
# ============================================================================

"""
	PyrolysisProblem(; mesh, material, bc_set::BoundaryConditionSet,
					 T_initial=z->300.0, ξ_initial=nothing,
					 external_flux=t->0.0, ambient_temperature=t->300.0,
					 tspan, dt_initial=1e-4, output_times=Float64[])

Create a PyrolysisProblem with BoundaryConditionSet (preferred, type-stable).

# Arguments
- `mesh`: The computational mesh
- `material`: Material definition
- `bc_set`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial=z->300.0`: Initial temperature function
- `ξ_initial=nothing`: Initial concentration functions (defaults from material)
- `external_flux=t->0.0`: External heat flux function (fallback if no RadiativeFluxBC)
- `ambient_temperature=t->300.0`: Ambient temperature function
- `tspan`: Simulation time span (required)
- `dt_initial=1e-4`: Initial timestep
- `output_times=Float64[]`: Times to save output (auto-generated if empty)

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(Inf),
	mass_bottom = ImpermeableBC()
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
	external_flux = t -> 0.0,
	ambient_temperature = t -> 300.0,
	tspan::Tuple{Real, Real},
	dt_initial::Real = 1e-4,
	output_times::Vector{<:Real} = Float64[],
)
	NC = n_components(material)

	# Create default initial concentration functions if not provided
	if ξ_initial === nothing
		# Default: first solid component at full density, others at zero
		ξ_funcs = Vector{Function}(undef, NC)
		for j in 1:NC
			comp = material.components[j]
			if j == 1 && comp.phase == SOLID
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
		ξ_initial = [ξ_initial for _ in 1:NC]
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
		external_flux,
		ambient_temperature,
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
	thermal_top = HeatFluxBC(50e3),
	thermal_bottom = AdiabaticBC()
)
problem = PyrolysisProblem(mesh, pmma, bc_set, (0.0, 100.0))
```
"""
function PyrolysisProblem(mesh, material, bc_set::BoundaryConditionSet, tspan::Tuple; kwargs...)
	PyrolysisProblem(; mesh = mesh, material = material, bc_set = bc_set,
		tspan = tspan, kwargs...)
end

# ============================================================================
# PyrolysisSolution Type
# ============================================================================

"""
	PyrolysisSolution{P<:AbstractPyrolysisProblem} <: AbstractPyrolysisSolution

Container for pyrolysis simulation results.

# Fields
- `problem::P`: Reference to the original problem
- `times::Vector{Float64}`: Solution times
- `temperatures::Vector{Vector{Float64}}`: Temperature profiles at each time
- `concentrations::Vector{Matrix{Float64}}`: Concentration profiles [NC × n_cells]
- `gas_generation_rate::Vector{Float64}`: Volumetric gas generation rate at each time [kg/(m²·s)]
- `surface_temperature::Vector{Float64}`: Surface temperature at each time [K]
- `back_temperature::Vector{Float64}`: Back face temperature at each time [K]
- `total_mass::Vector{Float64}`: Total mass in domain at each time [kg/m²]
- `retcode::Symbol`: Return code from ODE solver
- `solve_time::Float64`: Wall-clock solve time [s]
- `raw_solution`: Raw DifferentialEquations.jl solution object

# Important Note
The `gas_generation_rate` represents the volumetric production of gases from chemical
reactions throughout the sample, NOT the actual mass flux leaving the domain.
The actual mass loss rate (what a cone calorimeter measures) is stored in 
`surface_mass_flux` and is computed using the exact boundary condition discretization.

# Boundary Fluxes (computed during solve for accuracy)
These are computed during solution extraction using the exact same discretization
as the solver, ensuring consistency and accuracy:
- `surface_mass_flux`: Actual MLR at top surface [kg/(m²·s)]
- `surface_heat_flux_net`: Net heat entering at surface [W/m²]
- `surface_heat_flux_absorbed`: Absorbed radiation [W/m²]
- `surface_heat_flux_reradiation`: Re-radiation losses [W/m²]
- `surface_heat_flux_convection`: Convective losses [W/m²]
- `surface_heat_flux_conduction`: Conductive flux into material [W/m²]
- `surface_emissivity`: Effective emissivity at surface [-]
- `surface_conductivity`: Effective thermal conductivity at surface [W/(m·K)]
- `back_heat_flux`: Heat flux at back face [W/m²]

# Accessors
Use interpolation to get values at arbitrary times:
```julia
T_profile = sol(t)  # Temperature profile at time t
T_surface = surface_temperature(sol, t)
gen_rate = gas_generation_rate(sol, t)
```
"""
mutable struct PyrolysisSolution{P <: AbstractPyrolysisProblem} <: AbstractPyrolysisSolution
	problem::P
	times::Vector{Float64}
	temperatures::Vector{Vector{Float64}}
	concentrations::Vector{Matrix{Float64}}
	gas_generation_rate::Vector{Float64}
	surface_temperature::Vector{Float64}
	back_temperature::Vector{Float64}
	total_mass::Vector{Float64}
	thickness_history::Vector{Float64}  # Thickness at each saved time (for ALE)

	# === Boundary fluxes (computed during solve for accuracy) ===
	# Surface (boundary_id=1) fluxes
	surface_mass_flux::Vector{Float64}            # Total mass leaving [kg/(m²·s)] - actual MLR
	surface_heat_flux_net::Vector{Float64}        # Net heat entering [W/m²]
	surface_heat_flux_absorbed::Vector{Float64}   # Absorbed radiation [W/m²]
	surface_heat_flux_reradiation::Vector{Float64}  # Re-radiation loss [W/m²]
	surface_heat_flux_convection::Vector{Float64}   # Convective loss [W/m²]
	surface_heat_flux_conduction::Vector{Float64}   # Conduction into material [W/m²]

	# Surface effective properties
	surface_emissivity::Vector{Float64}           # ε_eff at surface [-]
	surface_conductivity::Vector{Float64}         # k_eff at surface [W/(m·K)]

	# Back face (boundary_id=2) fluxes
	back_heat_flux::Vector{Float64}               # Heat flux at back [W/m²]

	retcode::Symbol
	solve_time::Float64
	raw_solution::Any  # Store raw ODE solution for interpolation
end

"""
	PyrolysisSolution(problem)

Create an empty solution container for a problem.
"""
function PyrolysisSolution(problem::PyrolysisProblem)
	PyrolysisSolution(
		problem,
		Float64[],
		Vector{Float64}[],
		Matrix{Float64}[],
		Float64[],
		Float64[],
		Float64[],
		Float64[],
		Float64[],  # thickness_history
		# Boundary fluxes
		Float64[],  # surface_mass_flux
		Float64[],  # surface_heat_flux_net
		Float64[],  # surface_heat_flux_absorbed
		Float64[],  # surface_heat_flux_reradiation
		Float64[],  # surface_heat_flux_convection
		Float64[],  # surface_heat_flux_conduction
		Float64[],  # surface_emissivity
		Float64[],  # surface_conductivity
		Float64[],  # back_heat_flux
		:Default,
		0.0,
		nothing,
	)
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
get_times(sol::PyrolysisSolution) = sol.times

"""
	get_temperatures(sol::PyrolysisSolution)

Get all temperature profiles.
"""
get_temperatures(sol::PyrolysisSolution) = sol.temperatures

"""
	get_concentrations(sol::PyrolysisSolution)

Get all concentration profiles.
"""
get_concentrations(sol::PyrolysisSolution) = sol.concentrations

"""
	get_gas_generation_rate(sol::PyrolysisSolution)

Get volumetric gas generation rate history.

Note: This is NOT the actual mass loss rate (surface flux). To compute the actual MLR,
calculate the surface flux: MLR = h_m * ξ_surface for convective mass BC.
"""
get_gas_generation_rate(sol::PyrolysisSolution) = sol.gas_generation_rate

"""
	n_timesteps(sol::PyrolysisSolution)

Get number of saved timesteps.
"""
n_timesteps(sol::PyrolysisSolution) = length(sol.times)

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
	print(io, "PyrolysisProblem($(prob.material.name), $n_cells cells, $NC components, t∈[$t_start, $t_end])")
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
	println(io, "  Thermal BCs: top=$(typeof(prob.bc_set.thermal_top)), bottom=$(typeof(prob.bc_set.thermal_bottom))")
	println(io, "  Mass BCs: top=$(typeof(prob.bc_set.mass_top)), bottom=$(typeof(prob.bc_set.mass_bottom))")
end

function Base.show(io::IO, sol::PyrolysisSolution)
	n_times = length(sol.times)
	print(io, "PyrolysisSolution($(sol.retcode), $n_times timesteps, $(sol.solve_time)s)")
end

function Base.show(io::IO, ::MIME"text/plain", sol::PyrolysisSolution)
	println(io, "PyrolysisSolution")
	println(io, "  Return code: $(sol.retcode)")
	println(io, "  Timesteps saved: $(length(sol.times))")
	if !isempty(sol.times)
		println(io, "  Time range: [$(sol.times[1]), $(sol.times[end])] s")
	end
	println(io, "  Solve time: $(round(sol.solve_time, digits=3)) s")
	if !isempty(sol.gas_generation_rate)
		max_gen = maximum(sol.gas_generation_rate)
		println(io, "  Peak gas generation: $(round(max_gen * 1000, digits=2)) g/(m²·s)")
	end
	if !isempty(sol.surface_temperature)
		max_T = maximum(sol.surface_temperature)
		println(io, "  Peak surface T: $(round(max_T, digits=1)) K")
	end
end
