"""
Comprehensive Profiling Script for Douglas Fir Pyrolysis Simulation

This script profiles the solve operation to identify:
1. Where time is being spent (hotspots)
2. Where allocations occur (memory pressure)
3. Comparison of analytical vs finite difference Jacobian
4. Opportunities for optimization

Uses simple, reliable timing methods (@time, @elapsed, @allocated) to avoid
GC issues with complex mutable structures. Uses Profile.jl for CPU profiling.

Uses a smaller mesh (fewer cells) and shorter time span for faster profiling.

Usage:
	julia --project=. examples/profile_douglas_fir.jl

For detailed profiling with VS Code:
	1. Run this script to see allocation/timing breakdown
	2. Use @profview for flame graph visualization
"""

using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Properties
using Pyrolysis.Properties: PropertyCache, update_properties!, effective_permeability,
	volume_fractions, effective_conductivity, effective_heat_capacity
using Pyrolysis.Physics
using Pyrolysis.Physics: ideal_gas_pressure, evaluate_reactions!
using Pyrolysis.Discretization
using Pyrolysis.Discretization: compute_rhs!, compute_rhs_darcy_fick!, RHSCache,
	update_darcy_properties!, boundary_face_ids
using Pyrolysis.TimeIntegration
using Pyrolysis.TimeIntegration: PyrolysisODEFunctionALE, StateVectorALE,
	create_state_vector_ale, unpack_state_ale!, JacobianCacheALE,
	compute_analytical_jacobian_ale!
using Pyrolysis.Geometry: GeometryCache
using Profile
using Printf
using InteractiveUtils
using LinearSolve: KLUFactorization

# ============================================================================
# Configuration - Smaller problem for profiling
# ============================================================================

# Profiling configuration
const PROFILE_N_CELLS = 40        # Smaller mesh for profiling (vs 240 in full example)
const PROFILE_T_FINAL = 50.0      # Short simulation (vs 1000s in full example)
const PROFILE_OUTPUT_POINTS = 51  # Fewer output points
const WARMUP_RUNS = 2             # Number of warmup iterations

println("=" ^ 80)
println("Douglas Fir Pyrolysis - Profiling Script")
println("=" ^ 80)
println()

# ============================================================================
# Material Definition (same as main example)
# ============================================================================

println("Setting up material and problem...")

# Temperature-dependent heat capacity functions
df_cp = CallableProperty(T -> -2131.0 + 9.119 * T)
df_int1_cp = CallableProperty(T -> -1983.0 + 8.7 * T)
df_int2_cp = CallableProperty(T -> -1419.0 + 7.0 * T)
df_int3_cp = CallableProperty(T -> 1.0 + 2.8 * T)
df_int4_cp = CallableProperty(T -> 176.0 + 2.3 * T)
char_cp = CallableProperty(T -> 347.2 + 1.784 * T)
moisture_cp = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2)
water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)

# State-dependent heat of sorption (simplified)
heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

# Gas transfer coefficients (matching 07_douglas_fir_cone.jl)
λ_wood = 1.0e-6     # m²/s - virgin wood (tight structure, tortuosity ~3-5)
λ_char = 1.0e-5     # m²/s - porous char (open structure, porosity ~0.7-0.8)
λ_gas = 1.0e-5      # m²/s - gas phase diffusion in porous media

# Permeability parameters for Darcy flow (matching 07_douglas_fir_cone.jl)
κ_virgin = 2.3e-15  # m² - matches ThermaKin DF transport (λ_TK=1e-5, λ_PJ=1e-6)
κ_char = 1.0e-14    # m² - matches ThermaKin char transport (λ_TK=5e-5, λ_PJ=1e-5)

# Create Douglas Fir material with component-based permeability (matching 07_douglas_fir_cone.jl)
const douglas_fir = Material(
	name = :DouglasFir,
	components = (
		SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood, κ = κ_virgin),
		SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood*1.5, κ = 2.7e-15),
		SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood*3.0, κ = 3.1e-15),
		SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood*6.0, κ = 3.6e-15),
		SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char*0.8, κ = 8.3e-15),
		SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char, κ = κ_char),
		LiquidComponent(:MOISTURE, ρ = 1e6, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0, κ = 0.0),
		GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas*0.85),
		GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas*1.0),
		GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas*0.6),
		GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas*1.3),
		GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas*1.4),
		GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas*1.6),
	),
	reactions = (
		Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
		Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
		Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
		Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
		Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
		Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
	),
	conductivity_mixing = WEIGHTED,
	β = 0.5,
	gas_transfer_mixing = WEIGHTED,
	β_λ = 0.5,
	permeability_mixing = WEIGHTED,  # Mix component permeabilities
	β_κ = 0.5,                        # 50/50 parallel/series mixing
)

# ============================================================================
# Helper function to create fresh problem
# ============================================================================

function create_profiling_problem(n_cells, t_final)
	sample_thickness = 12.0e-3   # 12 mm
	T_ambient = 300.0
	q_external = 50e3
	h_conv = 7.2
	moisture_content_db = 0.15
	ρ_DF_initial = 530.0
	ρ_moisture_initial = moisture_content_db * ρ_DF_initial

	# Create mesh
	mesh = create_uniform_mesh(
		sample_thickness, n_cells, Val(13);
		T_initial = T_ambient,
		ξ_initial = (ρ_DF_initial, 0.0, 0.0, 0.0, 0.0, 0.0, ρ_moisture_initial,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
	)

	# Boundary conditions (matching 07_douglas_fir_cone.jl)
	top_thermal_bc = RadiativeFluxBC(flux = q_external, absorptivity = 0.9) +
					 RadiativeBC(ε = 0.9, T_inf = 330.0) +
					 ConvectiveBC(h = h_conv, T_inf = 330.0)
	top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)
	bottom_thermal_bc = AdiabaticBC()
	bottom_mass_bc = ImpermeableBC()

	bc_set = BoundaryConditionSet(
		thermal_top = top_thermal_bc,
		thermal_bottom = bottom_thermal_bc,
		mass_top = top_mass_bc,
		mass_bottom = bottom_mass_bc,
	)

	# Initial condition functions
	T_initial_func = z -> T_ambient
	ξ_initial_funcs = [
		z -> ρ_DF_initial, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
		z -> ρ_moisture_initial, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
	]

	# Create problem
	return PyrolysisProblem(
		mesh = mesh,
		material = douglas_fir,
		bc_set = bc_set,
		T_initial = T_initial_func,
		ξ_initial = ξ_initial_funcs,
		external_flux = t -> 0.0,
		ambient_temperature = t -> T_ambient,
		tspan = (0.0, t_final),
		dt_initial = 1e-3,
		output_times = collect(range(0.0, t_final, length = 51)),
	)
end

# Create main problem
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)

println("  - Mesh: $(PROFILE_N_CELLS) cells")
println("  - DOFs: $(PROFILE_N_CELLS * 14) (N × (1 + NC))")
println("  - Time span: 0 to $(PROFILE_T_FINAL) s")
println()

# ============================================================================
# Section 0: Hot Path Allocation Diagnostics
# ============================================================================

println("=" ^ 80)
println("SECTION 0: Hot Path Allocation Diagnostics")
println("=" ^ 80)
println()
println("Testing allocations in every function on the hot path...")
println("Each function should show 0 bytes for optimal performance.")
println()

# --- Setup for allocation testing ---
const NC_TEST = 13  # Douglas Fir has 13 components
const NR_TEST = 6   # Douglas Fir has 6 reactions

# Create fresh mesh for testing
alloc_test_mesh = create_uniform_mesh(
	12.0e-3, PROFILE_N_CELLS, Val(NC_TEST);
	T_initial = 500.0,  # Hot temperature for active reactions
	ξ_initial = (400.0, 50.0, 50.0, 50.0, 30.0, 20.0, 50.0,  # Solids + moisture
		10.0, 5.0, 5.0, 5.0, 5.0, 5.0),  # Gases
)
n_cells_test = length(alloc_test_mesh.cells)
n_faces_test = length(alloc_test_mesh.faces)
n_nodes_test = length(alloc_test_mesh.nodes)

# Create caches
alloc_prop_cache = PropertyCache(n_cells_test, n_faces_test)
alloc_rhs_cache = RHSCache(n_cells_test, n_faces_test, NC_TEST, NR_TEST)

# Create BC set
alloc_bc_set = BoundaryConditionSet(
	thermal_top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.9) +
				  RadiativeBC(ε = 0.9, T_inf = 330) +
				  ConvectiveBC(h = 7.2, T_inf = 330),
	thermal_bottom = AdiabaticBC(),
	mass_top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
	mass_bottom = ImpermeableBC(),
)

# Pre-cache boundary face IDs
alloc_top_face_ids = boundary_face_ids(alloc_test_mesh, 1)
alloc_bottom_face_ids = boundary_face_ids(alloc_test_mesh, 2)

# Output arrays
alloc_dT = zeros(n_cells_test)
alloc_dξ = zeros(NC_TEST, n_cells_test)

# Get sample cell state for single-cell tests
test_state = alloc_test_mesh.cell_states[1]
test_ξ = test_state.ξ
test_T = test_state.T

# Work arrays for reaction rates
alloc_rates = zeros(NR_TEST)
alloc_heats = zeros(NR_TEST)

# --- Function to print allocation result ---
function print_alloc_result(name::String, allocs::Int; expected::Int = 0)
	status = allocs == expected ? "✓" : "✗"
	color = allocs == expected ? "\e[32m" : "\e[31m"
	reset = "\e[0m"
	@printf("  %-45s %10d bytes  %s%s%s\n", name, allocs, color, status, reset)
end

println("-" ^ 70)
println("1. Property Evaluation Functions (per-cell)")
println("-" ^ 70)

# Warmup
_ = volume_fractions(test_ξ, test_T, douglas_fir)
_ = effective_conductivity(test_ξ, test_T, douglas_fir)
_ = effective_heat_capacity(test_ξ, test_T, douglas_fir)
_ = effective_permeability(test_ξ, test_T, douglas_fir)
_ = ideal_gas_pressure(test_ξ, test_T, douglas_fir)
_ = Properties.effective_gas_viscosity(test_T)

# Test allocations
print_alloc_result("volume_fractions(ξ, T, material)",
	@allocated volume_fractions(test_ξ, test_T, douglas_fir))
print_alloc_result("effective_conductivity(ξ, T, material)",
	@allocated effective_conductivity(test_ξ, test_T, douglas_fir))
print_alloc_result("effective_heat_capacity(ξ, T, material)",
	@allocated effective_heat_capacity(test_ξ, test_T, douglas_fir))
print_alloc_result("effective_permeability(ξ, T, material)",
	@allocated effective_permeability(test_ξ, test_T, douglas_fir))
print_alloc_result("ideal_gas_pressure(ξ, T, material)",
	@allocated ideal_gas_pressure(test_ξ, test_T, douglas_fir))
print_alloc_result("effective_gas_viscosity(T)",
	@allocated Properties.effective_gas_viscosity(test_T))

println()
println("-" ^ 70)
println("2. Reaction Rate Computation")
println("-" ^ 70)

# Warmup
evaluate_reactions!(alloc_rates, alloc_heats, test_ξ, test_T, douglas_fir)

print_alloc_result("evaluate_reactions!(rates, heats, ξ, T, mat)",
	@allocated evaluate_reactions!(alloc_rates, alloc_heats, test_ξ, test_T, douglas_fir))

println()
println("-" ^ 70)
println("3. Property Cache Updates (full mesh)")
println("-" ^ 70)

# Warmup
update_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir)

print_alloc_result("update_properties!(cache, mesh, material)",
	@allocated update_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir))

# Test Darcy properties if available
if has_darcy_flow(douglas_fir)
	update_darcy_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir,
		alloc_bc_set, alloc_top_face_ids, alloc_bottom_face_ids)
	print_alloc_result("update_darcy_properties!(cache, mesh, material, ...)",
		@allocated update_darcy_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir,
			alloc_bc_set, alloc_top_face_ids, alloc_bottom_face_ids))
end

println()
println("-" ^ 70)
println("4. RHS Computation Functions")
println("-" ^ 70)

# Warmup - standard RHS
update_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir)
compute_rhs!(alloc_dT, alloc_dξ, alloc_test_mesh, alloc_prop_cache, alloc_rhs_cache,
	douglas_fir, alloc_bc_set, 0.0, alloc_top_face_ids, alloc_bottom_face_ids;
	I_radiation = 0.0, radiation_model = SURFACE_ABSORPTION)

print_alloc_result("compute_rhs!(dT, dξ, mesh, cache, ..., SURFACE_ABSORPTION)",
	@allocated compute_rhs!(alloc_dT, alloc_dξ, alloc_test_mesh, alloc_prop_cache, alloc_rhs_cache,
		douglas_fir, alloc_bc_set, 0.0, alloc_top_face_ids, alloc_bottom_face_ids;
		I_radiation = 0.0, radiation_model = SURFACE_ABSORPTION))

# Test Darcy-Fick RHS if Darcy flow is enabled
if has_darcy_flow(douglas_fir)
	update_darcy_properties!(alloc_prop_cache, alloc_test_mesh, douglas_fir,
		alloc_bc_set, alloc_top_face_ids, alloc_bottom_face_ids)
	compute_rhs_darcy_fick!(alloc_dT, alloc_dξ, alloc_test_mesh, alloc_prop_cache, alloc_rhs_cache,
		douglas_fir, alloc_bc_set, 0.0, alloc_top_face_ids, alloc_bottom_face_ids;
		I_radiation = 0.0, radiation_model = SURFACE_ABSORPTION)

	print_alloc_result("compute_rhs_darcy_fick!(dT, dξ, ...)",
		@allocated compute_rhs_darcy_fick!(alloc_dT, alloc_dξ, alloc_test_mesh, alloc_prop_cache,
			alloc_rhs_cache, douglas_fir, alloc_bc_set, 0.0,
			alloc_top_face_ids, alloc_bottom_face_ids;
			I_radiation = 0.0, radiation_model = SURFACE_ABSORPTION))
end

println()
println("-" ^ 70)
println("5. ALE ODE Function (full state vector RHS)")
println("-" ^ 70)

# Create ALE ODE function
alloc_ode_fn = PyrolysisODEFunctionALE(alloc_test_mesh, douglas_fir, alloc_bc_set;
	I_radiation = 0.0, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, min_density = 1.0)

# Create state vector
alloc_u0 = create_state_vector_ale(alloc_test_mesh, douglas_fir)
alloc_du = similar(alloc_u0)
n_dof_test = length(alloc_u0)

# Warmup ODE function call
alloc_ode_fn(alloc_du, alloc_u0, nothing, 0.0)

print_alloc_result("PyrolysisODEFunctionALE(du, u, p, t) [full RHS]",
	@allocated alloc_ode_fn(alloc_du, alloc_u0, nothing, 0.0))

println()
println("-" ^ 70)
println("6. Analytical Jacobian Computation")
println("-" ^ 70)

# Create Jacobian cache and matrix
alloc_jac_cache = JacobianCacheALE(n_cells_test, n_faces_test, n_nodes_test, NC_TEST, NR_TEST)
alloc_J = zeros(n_dof_test, n_dof_test)

# Warmup
compute_analytical_jacobian_ale!(alloc_J, alloc_u0, nothing, 0.0, alloc_ode_fn, alloc_jac_cache)

print_alloc_result("compute_analytical_jacobian_ale!(J, u, p, t, ode_fn, cache)",
	@allocated compute_analytical_jacobian_ale!(alloc_J, alloc_u0, nothing, 0.0, alloc_ode_fn, alloc_jac_cache))

println()
println("-" ^ 70)
println("7. Geometry Cache Updates")
println("-" ^ 70)

# Create geometry cache
alloc_geom_cache = GeometryCache(n_cells_test, n_faces_test)

# update_geometry_cache! requires the full state vector and StateVectorALE
# It's called internally by PyrolysisODEFunctionALE, so test it through the full ODE call
Pyrolysis.Geometry.update_geometry_cache!(alloc_geom_cache, alloc_u0, alloc_ode_fn.sv, alloc_test_mesh)

print_alloc_result("update_geometry_cache!(cache, u, sv, mesh)",
	@allocated Pyrolysis.Geometry.update_geometry_cache!(alloc_geom_cache, alloc_u0, alloc_ode_fn.sv, alloc_test_mesh))

println()
println("=" ^ 70)
println("Hot Path Allocation Summary")
println("=" ^ 70)

# Run all tests again and sum allocations
total_allocs = 0

# Property evaluations (run in a loop like real code would)
function test_property_loop(mesh, material)
	allocs = 0
	for state in mesh.cell_states
		allocs += @allocated volume_fractions(state.ξ, state.T, material)
		allocs += @allocated effective_conductivity(state.ξ, state.T, material)
		allocs += @allocated effective_heat_capacity(state.ξ, state.T, material)
	end
	return allocs
end

# Warmup
_ = test_property_loop(alloc_test_mesh, douglas_fir)
prop_loop_allocs = @allocated test_property_loop(alloc_test_mesh, douglas_fir)
total_allocs += prop_loop_allocs

# Full RHS
rhs_allocs = @allocated alloc_ode_fn(alloc_du, alloc_u0, nothing, 0.0)
total_allocs += rhs_allocs

# Jacobian
jac_allocs = @allocated compute_analytical_jacobian_ale!(alloc_J, alloc_u0, nothing, 0.0, alloc_ode_fn, alloc_jac_cache)
total_allocs += jac_allocs


println()
@printf("Property loop (%d cells):     %10d bytes\n", n_cells_test, prop_loop_allocs)
@printf("Full ODE RHS evaluation:       %10d bytes\n", rhs_allocs)
@printf("Analytical Jacobian:           %10d bytes\n", jac_allocs)
println("-" ^ 45)
@printf("Total hot path allocations:    %10d bytes\n", total_allocs)

if total_allocs == 0
	println("\n✓ All hot path functions are allocation-free!")
else
	println("\n⚠ WARNING: Hot path has allocations - investigate above functions marked with ✗")
end
println()

# ============================================================================
# Section 1: Warmup and Basic Timing
# ============================================================================

println("=" ^ 80)
println("SECTION 1: Basic Timing")
println("=" ^ 80)
println()

println("Warming up solver ($(WARMUP_RUNS) runs)...")
for i in 1:WARMUP_RUNS
	prob_warmup = create_profiling_problem(PROFILE_N_CELLS, 5.0)  # Short warmup
	_ = solve(prob_warmup; abstol = 1e-8, reltol = 1e-6, progress = false,
		use_ale = true, radiation_model = SURFACE_ABSORPTION,
		min_thickness = 0.5e-3, use_analytical_jacobian = true)
end
println("Warmup complete.")
println()

# Fresh problem for timing
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)

# Baseline timing with @time
println("Baseline timing:")
@time solution = solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)
println("  - Return code: $(solution.retcode)")
println("  - Solver wall time: $(round(solution.solve_time, digits=3)) s")

# ============================================================================
# Section 2: Allocation Analysis
# ============================================================================

println()
println("=" ^ 80)
println("SECTION 2: Allocation Analysis")
println("=" ^ 80)
println()

# Track allocations per call - use fresh problem
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)
println("Measuring allocations...")
allocs = @allocated solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)
println("  Total allocations: $(round(allocs / 1e6, digits=2)) MB")
println("  Per timestep estimate: $(round(allocs / 1e6 / 50 * PROFILE_T_FINAL, digits=2)) MB/s simulated")
println()

# ============================================================================
# Section 3: Jacobian Comparison
# ============================================================================

println("=" ^ 80)
println("SECTION 3: Jacobian Comparison")
println("=" ^ 80)
println()

println("3.1 Solve with Analytical Jacobian:")
println("-" ^ 40)
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)
time_aj = @elapsed sol_aj = solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)
memory_aj = (@allocated solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)) / 1e6
println("  Time: $(round(time_aj, digits=3)) s")
println("  Memory: $(round(memory_aj, digits=2)) MB")
println("  Return code: $(sol_aj.retcode)")
println()

println("3.2 Solve with Finite Difference Jacobian:")
println("-" ^ 40)
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)
time_fd = @elapsed sol_fd = solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = false)
memory_fd = (@allocated solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = false)) / 1e6
println("  Time: $(round(time_fd, digits=3)) s")
println("  Memory: $(round(memory_fd, digits=2)) MB")
println("  Return code: $(sol_fd.retcode)")
println()

println("Speedup from analytical Jacobian: $(round(time_fd / time_aj, digits=1))×")
println("Memory reduction: $(round(memory_fd / memory_aj, digits=1))×")
println()

# ============================================================================
# Section 4: Darcy Flow Impact Analysis
# ============================================================================

println("=" ^ 80)
println("SECTION 4: Darcy Flow Impact Analysis")
println("=" ^ 80)
println()

# Note: With component-based permeability, Darcy flow is automatically enabled
# when any solid component has κ > 0. The douglas_fir material has κ values set
# on all solid components, so Darcy flow is enabled by default.
#
# For this profiling section, we simply run with the default material which
# includes Darcy flow. To compare with diffusion-only, you would need to create
# a separate material with all κ = 0.

println("4.1 Darcy-Fick Transport (component-based κ):")
println("-" ^ 40)
println("  Material has per-component permeability:")
println("    Virgin wood: κ = 1e-15 m² (tight structure)")
println("    Intermediates: κ = 5e-15 to 5e-13 m² (gradual opening)")
println("    Char: κ = 1e-12 m² (porous)")
println()

problem_darcy = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)
time_darcy = @elapsed sol_darcy = solve(problem_darcy; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)
memory_darcy = (@allocated solve(problem_darcy; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)) / 1e6
println("  Time: $(round(time_darcy, digits=3)) s")
println("  Memory: $(round(memory_darcy, digits=2)) MB")
println("  Return code: $(sol_darcy.retcode)")
println("  Steps taken: $(length(sol_darcy.times))")
println()

println("  Note: With component-based permeability, effective κ evolves naturally")
println("  as virgin wood converts to char - no temperature-based blending needed.")
println()

# ============================================================================
# Section 5: Mesh Size Scaling
# ============================================================================

println("=" ^ 80)
println("SECTION 5: Mesh Size Scaling")
println("=" ^ 80)
println()

# Use FD Jacobian for scaling test (faster per iteration)
mesh_sizes = [20, 40, 60]
times_by_size = Float64[]
allocs_by_size = Float64[]

println("Running mesh scaling tests (using FD Jacobian for speed)...")
for n in mesh_sizes
	println("  Testing N=$n cells...")
	prob = create_profiling_problem(n, 10.0)  # 10s runs for scaling
	t = @elapsed solve(prob; abstol = 1e-8, reltol = 1e-6, progress = false,
		use_ale = true, radiation_model = SURFACE_ABSORPTION,
		min_thickness = 0.5e-3, use_analytical_jacobian = false)
	prob2 = create_profiling_problem(n, 10.0)
	a = @allocated solve(prob2; abstol = 1e-8, reltol = 1e-6, progress = false,
		use_ale = true, radiation_model = SURFACE_ABSORPTION,
		min_thickness = 0.5e-3, use_analytical_jacobian = false)
	push!(times_by_size, t)
	push!(allocs_by_size, a)
end

println()
println("Mesh Size Scaling (10s simulation with FD Jacobian):")
println("-" ^ 60)
@printf("  %8s  %12s  %12s  %12s\n", "N cells", "Time (s)", "Allocs (MB)", "Scaling")
for (i, n) in enumerate(mesh_sizes)
	scaling = i > 1 ? times_by_size[i] / times_by_size[1] : 1.0
	@printf("  %8d  %12.2f  %12.1f  %12.2f×\n", n, times_by_size[i], allocs_by_size[i]/1e6, scaling)
end
println()

# Expected scaling for O(N²) operations (Jacobian):
expected_scaling = [(n / mesh_sizes[1])^2 for n in mesh_sizes]
println("Expected O(N²) scaling: $(round.(expected_scaling, digits=2))")
actual_scaling = times_by_size ./ times_by_size[1]
println("Actual scaling:         $(round.(actual_scaling, digits=2))")
println()

# ============================================================================
# Section 6: CPU Profile
# ============================================================================

println("=" ^ 80)
println("SECTION 6: CPU Time Profile")
println("=" ^ 80)
println()

println("Collecting CPU profile...")
println("Running profiled solve...")
problem = create_profiling_problem(PROFILE_N_CELLS, PROFILE_T_FINAL)
Profile.clear()
@profile solve(problem; abstol = 1e-8, reltol = 1e-6, progress = false,
	use_ale = true, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, use_analytical_jacobian = true)
println()

println("Top functions by CPU time (showing top 30):")
println("-" ^ 60)
Profile.print(format = :flat, mincount = 50, sortedby = :count, noisefloor = 2.0, maxdepth = 100)
println()

# ============================================================================
# Section 7: Component Profiling (using fresh meshes)
# ============================================================================

println("=" ^ 80)
println("SECTION 7: Component Profiling")
println("=" ^ 80)
println()

# Import internal functions for profiling
using Pyrolysis: PropertyCache, update_properties!, RHSCache, compute_rhs!, resize_rhs_cache!
using Pyrolysis: compute_conductive_fluxes!, compute_gas_fluxes!
using Pyrolysis: divergence_1d!, divergence_species_1d!
using Pyrolysis: evaluate_reactions!, reaction_heat_source, species_source_terms!

println("Creating fresh mesh for component profiling...")

sample_thickness = 12.0e-3
T_ambient = 300.0
q_external = 50e3
ρ_DF_initial = 530.0
ρ_moisture_initial = 0.15 * ρ_DF_initial

# Create a FRESH mesh for profiling
mesh_profile = create_uniform_mesh(
	sample_thickness, PROFILE_N_CELLS, Val(13);
	T_initial = T_ambient,
	ξ_initial = (ρ_DF_initial, 0.0, 0.0, 0.0, 0.0, 0.0, ρ_moisture_initial,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
)

# Create property and RHS caches
NC = 13
NR = 6
n_cells = PROFILE_N_CELLS

prop_cache = PropertyCache(n_cells, NC)
rhs_cache = RHSCache(n_cells, n_cells + 1, NC, NR)

# Boundary conditions
top_thermal_bc = RadiativeFluxBC(flux = q_external, absorptivity = 0.9) +
				 RadiativeBC(ε = 0.9, T_inf = 330) +
				 ConvectiveBC(h = 7.2, T_inf = 330)
top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)
bottom_thermal_bc = AdiabaticBC()
bottom_mass_bc = ImpermeableBC()

bc_set = BoundaryConditionSet(
	thermal_top = top_thermal_bc,
	thermal_bottom = bottom_thermal_bc,
	mass_top = top_mass_bc,
	mass_bottom = bottom_mass_bc,
)

# Get face/cell IDs
top_face_ids = [length(mesh_profile.faces)]
bottom_face_ids = [1]

# Pre-allocate work arrays
dT = zeros(n_cells)
dξ = zeros(NC, n_cells)

println()
println("7.1 Timing Individual Components (manual timing to avoid GC issues):")
println("-" ^ 50)

# Number of iterations for accurate timing
const N_ITERS = 1000

# Property evaluation timing
println("\nProperty Evaluation:")
update_properties!(prop_cache, mesh_profile, douglas_fir)  # Warmup
alloc_prop = @allocated update_properties!(prop_cache, mesh_profile, douglas_fir)
t_prop = @elapsed for _ in 1:N_ITERS
	update_properties!(prop_cache, mesh_profile, douglas_fir)
end
t_prop_us = (t_prop / N_ITERS) * 1e6  # Convert to μs
println("  Time: $(round(t_prop_us, digits=2)) μs/call (avg of $N_ITERS)")
println("  Allocations: $(alloc_prop) bytes/call")

# RHS computation timing
println("\nFull RHS Computation:")
compute_rhs!(dT, dξ, mesh_profile, prop_cache, rhs_cache, douglas_fir, bc_set, 0.0,
	top_face_ids, bottom_face_ids;
	I_radiation = q_external*0.9, radiation_model = SURFACE_ABSORPTION)  # Warmup
alloc_rhs = @allocated compute_rhs!(dT, dξ, mesh_profile, prop_cache, rhs_cache, douglas_fir, bc_set, 0.0,
	top_face_ids, bottom_face_ids;
	I_radiation = q_external*0.9, radiation_model = SURFACE_ABSORPTION)
t_rhs = @elapsed for _ in 1:N_ITERS
	compute_rhs!(dT, dξ, mesh_profile, prop_cache, rhs_cache, douglas_fir, bc_set, 0.0,
		top_face_ids, bottom_face_ids;
		I_radiation = q_external*0.9, radiation_model = SURFACE_ABSORPTION)
end
t_rhs_us = (t_rhs / N_ITERS) * 1e6
println("  Time: $(round(t_rhs_us, digits=2)) μs/call (avg of $N_ITERS)")
println("  Allocations: $(alloc_rhs) bytes/call")

# Conductive fluxes timing
println("\nConductive Fluxes:")
compute_conductive_fluxes!(rhs_cache.q_cond, mesh_profile, prop_cache)  # Warmup
alloc_cond = @allocated compute_conductive_fluxes!(rhs_cache.q_cond, mesh_profile, prop_cache)
t_cond = @elapsed for _ in 1:N_ITERS
	compute_conductive_fluxes!(rhs_cache.q_cond, mesh_profile, prop_cache)
end
t_cond_us = (t_cond / N_ITERS) * 1e6
println("  Time: $(round(t_cond_us, digits=2)) μs/call (avg of $N_ITERS)")
println("  Allocations: $(alloc_cond) bytes/call")

# Gas fluxes timing
println("\nGas Transport Fluxes:")
compute_gas_fluxes!(rhs_cache.J_gas, mesh_profile, prop_cache, douglas_fir)  # Warmup
alloc_gas = @allocated compute_gas_fluxes!(rhs_cache.J_gas, mesh_profile, prop_cache, douglas_fir)
t_gas = @elapsed for _ in 1:N_ITERS
	compute_gas_fluxes!(rhs_cache.J_gas, mesh_profile, prop_cache, douglas_fir)
end
t_gas_us = (t_gas / N_ITERS) * 1e6
println("  Time: $(round(t_gas_us, digits=2)) μs/call (avg of $N_ITERS)")
println("  Allocations: $(alloc_gas) bytes/call")

# Compare with Darcy flow enabled
println("\nGas Transport with Darcy Flow:")
using Pyrolysis: compute_gas_fluxes_darcy!
# Add permeability to property cache for Darcy
prop_cache_darcy = PropertyCache(n_cells, NC)
update_properties!(prop_cache_darcy, mesh_profile, douglas_fir; perm_model = perm_model)
compute_gas_fluxes_darcy!(rhs_cache.J_gas, mesh_profile, prop_cache_darcy, douglas_fir)  # Warmup
alloc_gas_darcy = @allocated compute_gas_fluxes_darcy!(rhs_cache.J_gas, mesh_profile, prop_cache_darcy, douglas_fir)
t_gas_darcy = @elapsed for _ in 1:N_ITERS
	compute_gas_fluxes_darcy!(rhs_cache.J_gas, mesh_profile, prop_cache_darcy, douglas_fir)
end
t_gas_darcy_us = (t_gas_darcy / N_ITERS) * 1e6
println("  Time: $(round(t_gas_darcy_us, digits=2)) μs/call (avg of $N_ITERS)")
println("  Allocations: $(alloc_gas_darcy) bytes/call")
darcy_overhead_pct = ((t_gas_darcy_us - t_gas_us) / t_gas_us) * 100
println("  Overhead vs diffusion-only: $(round(darcy_overhead_pct, digits=1))%")

# Reaction rates timing (single cell) - use more iterations since it's fast
println("\nReaction Rates (single cell):")
rates = zeros(NR)
heats = zeros(NR)
state = mesh_profile.cell_states[1]
ξ = state.ξ
T = state.T
const N_ITERS_RXN = 100000
evaluate_reactions!(rates, heats, ξ, T, douglas_fir)  # Warmup
alloc_rxn = @allocated evaluate_reactions!(rates, heats, ξ, T, douglas_fir)
t_rxn = @elapsed for _ in 1:N_ITERS_RXN
	evaluate_reactions!(rates, heats, ξ, T, douglas_fir)
end
t_rxn_ns = (t_rxn / N_ITERS_RXN) * 1e9
println("  Time: $(round(t_rxn_ns, digits=2)) ns/call (avg of $N_ITERS_RXN)")
println("  Allocations: $(alloc_rxn) bytes/call")

# Store for summary (convert to consistent units)
prop_time = t_prop_us / 1e3  # ms
rhs_time = t_rhs_us / 1e3    # ms
cond_time = t_cond_us / 1e3  # ms
gas_time = t_gas_us / 1e3    # ms
rxn_time_us = t_rxn_ns / 1e3 # μs

prop_allocs = Int(alloc_prop)
rhs_allocs = Int(alloc_rhs)
cond_allocs = Int(alloc_cond)
gas_allocs = Int(alloc_gas)
rxn_allocs = Int(alloc_rxn)
println()

# ============================================================================
# Section 8: Type Stability Check
# ============================================================================

println()
println("=" ^ 80)
println("SECTION 8: Type Stability Analysis")
println("=" ^ 80)
println()

println("Checking key functions for type stability...")
println("(Any red/yellow output indicates type instability)")
println()

println("evaluate_reactions!:")
println("-" ^ 40)
@code_warntype evaluate_reactions!(rates, heats, ξ, T, douglas_fir)
println()

# ============================================================================
# Section 9: Summary and Recommendations
# ============================================================================

println("=" ^ 80)
println("PROFILING SUMMARY")
println("=" ^ 80)
println()

println("Component Timing Summary:")
println("-" ^ 50)
@printf("  %-25s %8.3f ms\n", "Property evaluation:", prop_time)
@printf("  %-25s %8.3f ms\n", "Full RHS computation:", rhs_time)
@printf("  %-25s %8.3f ms\n", "Conductive fluxes:", cond_time)
@printf("  %-25s %8.3f ms\n", "Gas transport fluxes:", gas_time)
@printf("  %-25s %8.3f μs\n", "Reaction rates (1 cell):", rxn_time_us)
println()

println("Allocation Summary:")
println("-" ^ 50)
@printf("  %-25s %8d bytes/call\n", "Property evaluation:", prop_allocs)
@printf("  %-25s %8d bytes/call\n", "Full RHS computation:", rhs_allocs)
@printf("  %-25s %8d bytes/call\n", "Conductive fluxes:", cond_allocs)
@printf("  %-25s %8d bytes/call\n", "Gas transport fluxes:", gas_allocs)
@printf("  %-25s %8d bytes/call\n", "Reaction rates:", rxn_allocs)
println()

println("Jacobian Analysis:")
println("-" ^ 50)
@printf("  Analytical Jacobian time: %8.2f s\n", time_aj)
@printf("  Finite Difference time:   %8.2f s\n", time_fd)
@printf("  Speedup:                  %8.1f×\n", time_fd / time_aj)
println()

# Identify issues
println("Potential Optimization Targets:")
println("-" ^ 50)

if prop_allocs > 0
	println("  ⚠ Property evaluation allocates ($(prop_allocs) allocs)")
	println("    → Consider pre-allocating interpolation buffers")
end

if rhs_allocs > 0
	println("  ⚠ RHS computation allocates ($(rhs_allocs) allocs)")
	println("    → Check for array slicing or temporary allocations")
end

if rxn_allocs > 0
	println("  ⚠ Reaction evaluation allocates ($(rxn_allocs) allocs)")
	println("    → Consider inlining or fusing operations")
end

# Scaling analysis
if length(actual_scaling) >= 3
	# Check if scaling is quadratic
	is_quadratic = actual_scaling[3] > 3.0 * actual_scaling[2]
	if is_quadratic
		println("  ℹ Scaling appears O(N²) - Jacobian operations dominate")
		println("    → Consider sparse Jacobian for large problems")
	else
		println("  ℹ Scaling appears sub-quadratic - RHS operations may dominate")
	end
end

println()
println("=" ^ 80)
println("Profiling complete!")
println("=" ^ 80)
println()
println("Next Steps:")
println("  1. Use Profile.print(format=:tree) for call graph")
println("  2. Use @profview (ProfileView.jl) for flame graph visualization")
println("  3. Focus on functions with >10% of CPU time")
println("  4. Eliminate allocations in hot loops")
println("  5. Check type stability with @code_warntype on suspect functions")
println()

# Save profile for external analysis
println("Saving profile data...")
open("profile_data.txt", "w") do io
	Profile.print(io, format = :flat, sortedby = :count, maxdepth = 100)
end
println("Profile saved to: profile_data.txt")
