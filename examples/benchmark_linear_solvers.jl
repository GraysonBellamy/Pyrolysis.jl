# Benchmark Linear Solvers: KLU vs UMFPACK
#
# This script compares the performance of KLU and UMFPACK linear solvers
# for the Douglas Fir pyrolysis simulation.
#
# Usage: julia --project examples/benchmark_linear_solvers.jl

using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using LinearSolve
using KLU
using OrdinaryDiffEq

println("="^60)
println("Linear Solver Benchmark: KLU vs UMFPACK")
println("="^60)

# ============================================================================
# Douglas Fir Material Definition
# ============================================================================

println("\nSetting up Douglas Fir material...")

# Heat capacities (temperature-dependent)
df_cp = CallableProperty(T -> 103.1 + 3.867*T)
df_int1_cp = CallableProperty(T -> 72.17 + 3.867*T)
df_int2_cp = CallableProperty(T -> 123.9 + 3.393*T)
df_int3_cp = CallableProperty(T -> 167.7 + 3.393*T)
df_int4_cp = CallableProperty(T -> -243.1 + 4.103*T)
char_cp = CallableProperty(T -> 193.7 + 2.273*T)
moisture_cp = CallableProperty(T -> 4190.0)
water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6*T + 2.0e-3*T^2)

# State-dependent heat of sorption
heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

# Gas transfer coefficients
λ_char = 5e-5
λ_wood = 1e-5
λ_gas = 1e-4

douglas_fir = Material(
	name = :DouglasFir,
	components = (
		SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, λ = λ_wood),
		SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, λ = λ_wood*1.2),
		SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, λ = λ_wood*1.5),
		SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, λ = λ_wood*2.0),
		SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, λ = λ_char*0.8),
		SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, λ = λ_char),
		LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9),
		GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
		GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
		GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas*0.8),
		GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
		GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
		GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas*1.2),
	),
	reactions = (
		Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
		Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
		Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
		Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
		Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
		Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
	),
	conductivity_mixing = WEIGHTED, β = 0.5,
	gas_transfer_mixing = WEIGHTED, β_λ = 0.5,
)

# ============================================================================
# Problem Setup
# ============================================================================

n_cells = 40
sample_thickness = 12.0e-3
mesh = create_uniform_mesh(sample_thickness, n_cells, Val(13))

# Boundary conditions
q_external = 50e3
T_ambient = 300.0
h_conv = 7.2
top_bc = HeatFluxBC(q_external) + RadiativeBC(ε = 0.9, T_inf = T_ambient) + ConvectiveBC(h = h_conv, T_inf = T_ambient)
bottom_bc = AdiabaticBC()
bc_set = BoundaryConditionSet(
	thermal_top = top_bc,
	thermal_bottom = bottom_bc,
	mass_top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
	mass_bottom = ImpermeableBC(),
)

# Initial conditions
ρ_DF_initial = 530.0
ρ_moisture_initial = 530.0 * 0.15
T_initial_func = z -> T_ambient

ξ_initial_funcs = [
	z -> ρ_DF_initial,     # 1. DF (virgin wood)
	z -> 0.0,              # 2. DF_int1
	z -> 0.0,              # 3. DF_int2
	z -> 0.0,              # 4. DF_int3
	z -> 0.0,              # 5. DF_int4
	z -> 0.0,              # 6. Char
	z -> ρ_moisture_initial, # 7. MOISTURE
	z -> 0.0,              # 8. DF_gas1
	z -> 0.0,              # 9. DF_gas2
	z -> 0.0,              # 10. DF_gas3
	z -> 0.0,              # 11. DF_gas4
	z -> 0.0,              # 12. DF_gas5
	z -> 0.0,              # 13. WATER_VAPOR
]

# Short simulation for benchmarking
benchmark_tspan = (0.0, 10.0)

problem = PyrolysisProblem(
	mesh = mesh,
	material = douglas_fir,
	bc_set = bc_set,
	T_initial = T_initial_func,
	ξ_initial = ξ_initial_funcs,
	tspan = benchmark_tspan,
)

println("Problem: $(n_cells) cells, 13 components, 6 reactions")
println("Time span: $(benchmark_tspan) seconds")
println("DOFs: $(n_cells * 14) (40 cells × 14 DOF/cell)")

# ============================================================================
# Benchmark Solvers
# ============================================================================

println("\n" * "="^60)
println("Running Benchmarks")
println("="^60)

# Test with default (UMFPACK)
println("\n--- UMFPACK (default) ---")
GC.gc()
t1 = @elapsed sol_umfpack = solve(problem,
	use_ale = true,
	min_thickness = 0.5e-3,
	use_analytical_jacobian = true,  # Enable sparse Jacobian
	solver = Rodas5P(),
)
println("Time: $(round(t1, digits=2))s")
println("Steps: $(length(sol_umfpack.times))")
println("Retcode: $(sol_umfpack.retcode)")

# Test with KLU
println("\n--- KLU ---")
GC.gc()
t2 = @elapsed sol_klu = solve(problem,
	use_ale = true,
	min_thickness = 0.5e-3,
	use_analytical_jacobian = true,  # Required for KLU (needs sparse matrix)
	solver = Rodas5P(linsolve = KLUFactorization()),
)
println("Time: $(round(t2, digits=2))s")
println("Steps: $(length(sol_klu.times))")
println("Retcode: $(sol_klu.retcode)")

# ============================================================================
# Results Summary
# ============================================================================

println("\n" * "="^60)
println("Results Summary")
println("="^60)
println("UMFPACK: $(round(t1, digits=2))s")
println("KLU:     $(round(t2, digits=2))s")

speedup = t1 / t2
if speedup > 1
	println("\nKLU is $(round(speedup, digits=2))× faster than UMFPACK")
else
	println("\nUMFPACK is $(round(1/speedup, digits=2))× faster than KLU")
end

# Verify solutions match
T_diff = maximum(abs.(sol_umfpack.temperatures[end] .- sol_klu.temperatures[end]))
println("\nMax temperature difference at t_final: $(round(T_diff, digits=4)) K")

if T_diff < 0.1
	println("✓ Solutions match (difference < 0.1 K)")
else
	println("⚠ Solutions differ significantly!")
end

println("\n" * "="^60)
println("Benchmark complete")
println("="^60)
