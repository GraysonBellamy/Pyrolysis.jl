#!/usr/bin/env julia
"""
Comprehensive Solver Comparison for Pyrolysis.jl

Tests all combinations of:
- ODE Solvers: Rodas5P, QNDF, FBDF, Rosenbrock23, KenCarp4
- Jacobian: Analytical vs Finite Difference
- Matrix: Sparse vs Dense
- Linear Solver: UMFPACK, KLU, (for sparse), LAPACK (for dense)

Outputs a summary table to help choose optimal settings.
"""

using Pkg
Pkg.activate(".")

using Pyrolysis
using OrdinaryDiffEq
using SparseArrays
using LinearAlgebra
using ADTypes
using Printf

# Import internals
using Pyrolysis.TimeIntegration: PyrolysisODEFunctionALE, StateVectorALE,
	JacobianCacheALE, compute_analytical_jacobian_ale!, jacobian_sparsity_ale,
	create_state_vector_ale
using Pyrolysis.Core: n_components, n_reactions

println("="^100)
println("COMPREHENSIVE SOLVER COMPARISON")
println("="^100)

# ============================================================================
# Create Douglas Fir problem
# ============================================================================

df_cp = CallableProperty(T -> -2131.0 + 9.119 * T)
df_int1_cp = CallableProperty(T -> -1983.0 + 8.7 * T)
df_int2_cp = CallableProperty(T -> -1419.0 + 7.0 * T)
df_int3_cp = CallableProperty(T -> 1.0 + 2.8 * T)
df_int4_cp = CallableProperty(T -> 176.0 + 2.3 * T)
char_cp = CallableProperty(T -> 347.2 + 1.784 * T)
moisture_cp = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2)
water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)
heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

λ_char = 5e-5
λ_wood = 1e-5
λ_gas = 1e-4
κ_virgin = 1e-15
κ_char = 1e-12

douglas_fir = Material(
	name = :DouglasFir,
	components = (
		SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood, κ = κ_virgin),
		SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood*1.2, κ = 5e-15),
		SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood*1.5, κ = 1e-14),
		SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood*2.0, κ = 1e-13),
		SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char*0.8, κ = 5e-13),
		SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char, κ = κ_char),
		LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0, κ = 0.0),
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
	conductivity_mixing = WEIGHTED,
	β = 0.5,
	gas_transfer_mixing = WEIGHTED,
	β_λ = 0.5,
	permeability_mixing = WEIGHTED,
	β_κ = 0.5,
)

n_cells = 400
sample_thickness = 12.0e-3
T_ambient = 300.0
q_external = 50e3

mesh = create_uniform_mesh(
	sample_thickness, n_cells, Val(13);
	T_initial = T_ambient,
	ξ_initial = (530.0, 0.0, 0.0, 0.0, 0.0, 0.0, 79.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
)

bc_set = BoundaryConditionSet(
	thermal_top = RadiativeFluxBC(flux = q_external) +
				  RadiativeBC(T_inf = 330) +
				  ConvectiveBC(h = 7.2, T_inf = 330),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
	mass_bottom = ImpermeableBC(),
)

material = douglas_fir
NC = n_components(material)
NR = n_reactions(material)

# Create ODE function
ode_fn = PyrolysisODEFunctionALE(
	mesh, material, bc_set;
	I_radiation = 50e3,
	radiation_model = Pyrolysis.Core.SURFACE_ABSORPTION,
	min_thickness = 1e-6,
)

u0 = create_state_vector_ale(mesh)
tspan = (0.0, 10.0)  # 10 seconds

println("Problem Setup:")
println("  Cells: $n_cells")
println("  Components: $NC")
println("  DOF: $(length(u0))")
println("  Simulation time: $(tspan[2]) s")
println()

# ============================================================================
# Prepare Jacobian infrastructure
# ============================================================================

n_faces = length(mesh.faces)
n_nodes = length(mesh.nodes)
jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
jac_prototype = jacobian_sparsity_ale(n_cells, NC)

function analytical_jac_ale!(J, u, p, t)
	compute_analytical_jacobian_ale!(J, u, p, t, ode_fn, jac_cache)
end

# Colorvec for sparse FD
using SparseDiffTools: matrix_colors
colorvec = matrix_colors(jac_prototype)

println("Jacobian Info:")
println("  Sparse: $(nnz(jac_prototype)) non-zeros ($(round(100*nnz(jac_prototype)/length(u0)^2, digits=2))% fill)")
println("  Dense: $(length(u0)^2) elements")
println("  Sparse memory: $(round(sizeof(jac_prototype.nzval)/1024, digits=1)) KB")
println("  Dense memory: $(round(length(u0)^2*8/1024, digits=1)) KB")
println()

# ============================================================================
# Define test configurations
# ============================================================================

using LinearSolve: KLUFactorization, UMFPACKFactorization, RFLUFactorization

struct SolverConfig
	name::String
	solver_fn::Function  # Returns configured solver
	ode_fn::Any          # ODEFunction (with or without jac)
	description::String
end

# Create different ODE functions
ode_fn_dense = ode_fn  # No jac_prototype = dense

ode_fn_sparse_fd = OrdinaryDiffEq.ODEFunction(ode_fn;
	jac_prototype = jac_prototype,
	colorvec = colorvec)

ode_fn_analytical_sparse = OrdinaryDiffEq.ODEFunction(ode_fn;
	jac = analytical_jac_ale!,
	jac_prototype = jac_prototype)

# Dense analytical Jacobian
jac_dense = zeros(length(u0), length(u0))
function analytical_jac_dense!(J, u, p, t)
	# Fill sparse then convert
	J_sparse = similar(jac_prototype)
	compute_analytical_jacobian_ale!(J_sparse, u, p, t, ode_fn, jac_cache)
	J .= Matrix(J_sparse)
end

ode_fn_analytical_dense = OrdinaryDiffEq.ODEFunction(ode_fn;
	jac = analytical_jac_dense!)

# Define all configurations to test
configs = SolverConfig[]

# ============================================================================
# Rosenbrock solvers (Rodas5P, Rosenbrock23)
# ============================================================================

# Dense FD (baseline - what we had before)
# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"Rodas5P + Dense FD",
# 	() -> Rodas5P(autodiff = AutoFiniteDiff()),
# 	ode_fn_dense,
# 	"Dense matrix, FD Jacobian, LAPACK LU",
# ))

# Sparse FD + UMFPACK
# push!(configs, SolverConfig(
# 	"Rodas5P + Sparse FD + UMFPACK",
# 	() -> Rodas5P(autodiff = AutoFiniteDiff()),
# 	ode_fn_sparse_fd,
# 	"Sparse matrix, colored FD Jacobian, UMFPACK",
# ))

# Sparse FD + KLU
# push!(configs, SolverConfig(
# 	"Rodas5P + Sparse FD + KLU",
# 	() -> Rodas5P(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization()),
# 	ode_fn_sparse_fd,
# 	"Sparse matrix, colored FD Jacobian, KLU",
# ))

# Analytical + UMFPACK
# push!(configs, SolverConfig(
# 	"Rodas5P + Analytical + UMFPACK",
# 	() -> Rodas5P(),
# 	ode_fn_analytical_sparse,
# 	"Sparse matrix, analytical Jacobian, UMFPACK",
# ))

# Analytical + KLU (expected best for Rosenbrock)
push!(configs, SolverConfig(
	"Rodas5P + Analytical + KLU",
	() -> Rodas5P(linsolve = KLUFactorization()),
	ode_fn_analytical_sparse,
	"Sparse matrix, analytical Jacobian, KLU",
))

# Dense Analytical
# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"Rodas5P + Analytical Dense",
# 	() -> Rodas5P(),
# 	ode_fn_analytical_dense,
# 	"Dense matrix, analytical Jacobian, LAPACK LU",
# ))

# ============================================================================
# BDF solvers (QNDF, FBDF)
# ============================================================================

# QNDF - Dense FD
# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"QNDF + Dense FD",
# 	() -> QNDF(autodiff = AutoFiniteDiff()),
# 	ode_fn_dense,
# 	"BDF, Dense FD, LAPACK",
# ))

# QNDF - Analytical + KLU
push!(configs, SolverConfig(
	"QNDF + Analytical + KLU",
	() -> QNDF(linsolve = KLUFactorization()),
	ode_fn_analytical_sparse,
	"BDF, Analytical sparse, KLU",
))

# FBDF - Dense FD
# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"FBDF + Dense FD",
# 	() -> FBDF(autodiff = AutoFiniteDiff()),
# 	ode_fn_dense,
# 	"Fixed-leading BDF, Dense FD, LAPACK",
# ))

# FBDF - Analytical + KLU
push!(configs, SolverConfig(
	"FBDF + Analytical + KLU",
	() -> FBDF(linsolve = KLUFactorization()),
	ode_fn_analytical_sparse,
	"Fixed-leading BDF, Analytical sparse, KLU",
))

# ============================================================================
# IMEX solvers (KenCarp4)
# ============================================================================

# KenCarp4 - Dense FD
# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"KenCarp4 + Dense FD",
# 	() -> KenCarp4(autodiff = AutoFiniteDiff()),
# 	ode_fn_dense,
# 	"IMEX ARK, Dense FD, LAPACK",
# ))

# KenCarp4 - Analytical + KLU
push!(configs, SolverConfig(
	"KenCarp4 + Analytical + KLU",
	() -> KenCarp4(linsolve = KLUFactorization()),
	ode_fn_analytical_sparse,
	"IMEX ARK, Analytical sparse, KLU",
))

# ============================================================================
# Rosenbrock23 (lower order, potentially faster per step)
# ============================================================================

# DISABLED FOR LARGE PROBLEM TEST
# push!(configs, SolverConfig(
# 	"Rosenbrock23 + Dense FD",
# 	() -> Rosenbrock23(autodiff = AutoFiniteDiff()),
# 	ode_fn_dense,
# 	"2nd order Rosenbrock, Dense FD",
# ))

push!(configs, SolverConfig(
	"Rosenbrock23 + Analytical + KLU",
	() -> Rosenbrock23(linsolve = KLUFactorization()),
	ode_fn_analytical_sparse,
	"2nd order Rosenbrock, Analytical sparse, KLU",
))

# ============================================================================
# Run benchmarks
# ============================================================================

struct BenchmarkResult
	name::String
	time::Float64
	alloc_mb::Float64
	nf::Int
	njacs::Int
	nw::Int
	naccept::Int
	nreject::Int
	success::Bool
	error_msg::String
end

results = BenchmarkResult[]

println("="^100)
println("RUNNING BENCHMARKS")
println("="^100)
println()

for (i, config) in enumerate(configs)
	println("[$i/$(length(configs))] $(config.name)")
	println("  $(config.description)")

	# Warmup
	try
		prob = ODEProblem(config.ode_fn, u0, (0.0, 0.1))
		sol = OrdinaryDiffEq.solve(prob, config.solver_fn();
			abstol = 1e-8, reltol = 1e-6, save_everystep = false, maxiters = 1000)
	catch e
		println("  Warmup failed: $(typeof(e))")
	end

	GC.gc()

	# Actual benchmark
	try
		prob = ODEProblem(config.ode_fn, u0, tspan)

		time_start = time()
		alloc = @allocated begin
			sol = OrdinaryDiffEq.solve(prob, config.solver_fn();
				abstol = 1e-8, reltol = 1e-6, save_everystep = false)
		end
		elapsed = time() - time_start

		if sol.retcode == ReturnCode.Success
			stats = sol.stats
			push!(results, BenchmarkResult(
				config.name,
				elapsed,
				alloc / 1e6,
				stats.nf,
				stats.njacs,
				stats.nw,
				stats.naccept,
				stats.nreject,
				true,
				"",
			))
			@printf("  ✓ Time: %.2fs, Alloc: %.1f MB, Steps: %d, Jacs: %d\n",
				elapsed, alloc/1e6, stats.naccept, stats.njacs)
		else
			push!(results, BenchmarkResult(
				config.name, elapsed, alloc/1e6, 0, 0, 0, 0, 0, false,
				"Solver returned: $(sol.retcode)",
			))
			println("  ✗ Failed: $(sol.retcode)")
		end
	catch e
		push!(results, BenchmarkResult(
			config.name, Inf, Inf, 0, 0, 0, 0, 0, false,
			"Exception: $(typeof(e))",
		))
		println("  ✗ Error: $(typeof(e))")
	end
	println()
end

# ============================================================================
# Results Summary
# ============================================================================

println("="^100)
println("RESULTS SUMMARY (sorted by time)")
println("="^100)
println()

# Sort by time
sorted_results = sort(filter(r -> r.success, results), by = r -> r.time)

# Find best for normalization
if !isempty(sorted_results)
	best_time = sorted_results[1].time
	best_alloc = minimum(r.alloc_mb for r in sorted_results)

	# Header
	@printf("%-40s %8s %8s %10s %8s %8s %8s\n",
		"Configuration", "Time(s)", "vs Best", "Alloc(MB)", "vs Best", "Steps", "Jacs")
	println("-"^100)

	for r in sorted_results
		time_ratio = r.time / best_time
		alloc_ratio = r.alloc_mb / best_alloc
		@printf("%-40s %8.2f %7.1fx %10.1f %7.1fx %8d %8d\n",
			r.name, r.time, time_ratio, r.alloc_mb, alloc_ratio, r.naccept, r.njacs)
	end

	println()
	println("-"^100)
	println("Failed configurations:")
	for r in filter(r -> !r.success, results)
		println("  $(r.name): $(r.error_msg)")
	end
end

# ============================================================================
# Detailed Analysis
# ============================================================================

println()
println("="^100)
println("ANALYSIS")
println("="^100)

if !isempty(sorted_results)
	best = sorted_results[1]
	println()
	println("WINNER: $(best.name)")
	println("  Time: $(round(best.time, digits=2)) s")
	println("  Allocations: $(round(best.alloc_mb, digits=1)) MB")
	println("  Steps: $(best.naccept)")
	println("  Jacobian evaluations: $(best.njacs)")
	println()

	# Group analysis
	println("By Jacobian Type (best in each category):")

	fd_dense = filter(r -> contains(r.name, "Dense FD"), sorted_results)
	fd_sparse = filter(r -> contains(r.name, "Sparse FD"), sorted_results)
	analytical = filter(r -> contains(r.name, "Analytical"), sorted_results)

	if !isempty(fd_dense)
		r = fd_dense[1]
		@printf("  Dense FD:     %-35s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end
	if !isempty(fd_sparse)
		r = fd_sparse[1]
		@printf("  Sparse FD:    %-35s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end
	if !isempty(analytical)
		r = analytical[1]
		@printf("  Analytical:   %-35s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end

	println()
	println("By Linear Solver (best in each category):")

	umfpack = filter(r -> contains(r.name, "UMFPACK"), sorted_results)
	klu = filter(r -> contains(r.name, "KLU"), sorted_results)
	lapack = filter(r -> !contains(r.name, "UMFPACK") && !contains(r.name, "KLU"), sorted_results)

	if !isempty(lapack)
		r = lapack[1]
		@printf("  LAPACK (dense): %-33s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end
	if !isempty(umfpack)
		r = umfpack[1]
		@printf("  UMFPACK:        %-33s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end
	if !isempty(klu)
		r = klu[1]
		@printf("  KLU:            %-33s  %.2fs, %.1f MB\n", r.name, r.time, r.alloc_mb)
	end

	println()
	println("By ODE Solver (best in each category):")

	for solver_name in ["Rodas5P", "QNDF", "FBDF", "KenCarp4", "Rosenbrock23"]
		matching = filter(r -> startswith(r.name, solver_name), sorted_results)
		if !isempty(matching)
			r = matching[1]
			@printf("  %-12s: %-33s  %.2fs, %.1f MB\n", solver_name, r.name, r.time, r.alloc_mb)
		end
	end
end

# ============================================================================
# Recommendations
# ============================================================================

println()
println("="^100)
println("RECOMMENDATIONS")
println("="^100)
println()

if !isempty(sorted_results)
	best = sorted_results[1]

	# Find best analytical
	analytical_results = filter(r -> contains(r.name, "Analytical"), sorted_results)
	best_analytical = isempty(analytical_results) ? nothing : analytical_results[1]

	# Find best for low memory
	low_mem = sort(sorted_results, by = r -> r.alloc_mb)[1]

	println("For FASTEST execution:")
	println("  $(best.name)")
	println()

	if best_analytical !== nothing
		println("For BEST analytical Jacobian:")
		println("  $(best_analytical.name)")
		println("  (Analytical is better when Jacobian is expensive to compute via FD)")
		println()
	end

	println("For LOWEST memory:")
	println("  $(low_mem.name)")
	println()

	println("General recommendations:")
	println("  1. Use KLU instead of UMFPACK for sparse linear solves")
	println("  2. Analytical Jacobian reduces RHS evaluations significantly")
	println("  3. Dense matrices work well for small problems (< 1000 DOF)")
	println("  4. For larger problems, sparse + KLU will scale better")
end

# ============================================================================
# Export results to CSV
# ============================================================================

using DelimitedFiles

csv_data = Any["Name" "Time_s" "Alloc_MB" "RHS_calls" "Jac_evals" "W_factors" "Steps" "Rejects" "Success"]
for r in results
	csv_data = vcat(csv_data, Any[r.name r.time r.alloc_mb r.nf r.njacs r.nw r.naccept r.nreject r.success])
end

csv_file = "solver_comparison_results.csv"
writedlm(csv_file, csv_data, ',')
println()
println("Results saved to: $csv_file")
