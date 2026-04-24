"""
Tests for ILU-GMRES Linear Solver Integration

Verifies that ILUGMRESFactorization works as a drop-in replacement for
KLUFactorization in ALE pyrolysis simulations with analytical Jacobian.
"""

using Test
using Pyrolysis
using OrdinaryDiffEq
using LinearSolve: KLUFactorization

@testset "ILU-GMRES Linear Solver" begin

	# ========================================================================
	# Helper: Create a simple ALE-compatible test problem
	# ========================================================================
	function create_test_problem(; n_cells = 20, t_final = 5.0)
		virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

		material = Material(
			name = :PMMA,
			components = (virgin, gas),
			reactions = (rxn,),
		)

		mesh = create_uniform_mesh(0.006, n_cells, Val(2);
			T_initial = 300.0,
			ξ_initial = (1190.0, 0.0))

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.86, T_inf = 300.0),
			thermal_bottom = AdiabaticBC(),
			mass_top = ConvectiveMassBC(h_m = 0.01, ξ_inf = 0.0),
			mass_bottom = ImpermeableBC(),
		)

		return PyrolysisProblem(
			mesh = mesh, material = material, bc_set = bc_set,
			T_initial = z -> 300.0,
			ξ_initial = [z -> 1190.0, z -> 0.0],
			tspan = (0.0, t_final),
			dt_initial = 1e-3,
			output_times = collect(range(0.0, t_final, length = 11)),
		)
	end

	# ========================================================================
	@testset "Basic solve with ILU-GMRES" begin
		problem = create_test_problem()
		sol = solve(problem;
			solver = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
			use_ale = true,
			use_analytical_jacobian = true,
			min_thickness = 0.1e-3,
			abstol = 1e-8, reltol = 1e-6,
			progress = false)

		@test sol.retcode == :Success
		@test length(sol.times) > 1
		@test all(isfinite, sol.surface_temperature)
		@test sol.surface_temperature[end] > 300.0  # Should heat up
	end

	# ========================================================================
	@testset "Solution matches KLU baseline" begin
		problem_klu = create_test_problem()
		sol_klu = solve(problem_klu;
			solver = KenCarp4(linsolve = KLUFactorization()),
			use_ale = true,
			use_analytical_jacobian = true,
			min_thickness = 0.1e-3,
			abstol = 1e-8, reltol = 1e-6,
			progress = false)

		problem_ilu = create_test_problem()
		sol_ilu = solve(problem_ilu;
			solver = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
			use_ale = true,
			use_analytical_jacobian = true,
			min_thickness = 0.1e-3,
			abstol = 1e-8, reltol = 1e-6,
			progress = false)

		@test sol_klu.retcode == :Success
		@test sol_ilu.retcode == :Success

		# Surface temperatures should match within 1 K
		T_diff = maximum(abs.(sol_klu.surface_temperature .- sol_ilu.surface_temperature))
		@test T_diff < 1.0
	end

	# ========================================================================
	@testset "Parameter variations" begin
		# τ = 0.01 (tighter ILU)
		problem1 = create_test_problem()
		sol1 = solve(problem1;
			solver = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.01)),
			use_ale = true,
			use_analytical_jacobian = true,
			min_thickness = 0.1e-3,
			progress = false)
		@test sol1.retcode == :Success

		# τ = 0.5 (sparser ILU)
		problem2 = create_test_problem()
		sol2 = solve(problem2;
			solver = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.5)),
			use_ale = true,
			use_analytical_jacobian = true,
			min_thickness = 0.1e-3,
			progress = false)
		@test sol2.retcode == :Success
	end

	# ========================================================================
	@testset "Constructor and display" begin
		alg = ILUGMRESFactorization()
		@test alg.τ == 0.1
		@test alg.gmres_rtol == 1e-8
		@test alg.restart == 30

		alg2 = ILUGMRESFactorization(τ = 0.05, gmres_rtol = 1e-6, restart = 50)
		@test alg2.τ == 0.05
		@test alg2.gmres_rtol == 1e-6
		@test alg2.restart == 50

		# Test show method
		buf = IOBuffer()
		show(buf, alg)
		s = String(take!(buf))
		@test occursin("ILUGMRESFactorization", s)
		@test occursin("0.1", s)
	end

end
