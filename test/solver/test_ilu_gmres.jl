"""
Tests for ILU-GMRES Linear Solver Integration

Verifies that ILUGMRESFactorization works as a drop-in replacement for
KLUFactorization in ALE pyrolysis simulations with analytical Jacobian.
"""

using Test
using Pyrolysis
using Pyrolysis.Solver:
    ILUGMRESFactorization, ILUBiCGSTABFactorization, KenCarp4, _krylov_solution, Krylov
using OrdinaryDiffEq
using LinearSolve: LinearSolve, KLUFactorization, LinearProblem
using SparseArrays: sparse
const ReturnCode = LinearSolve.SciMLBase.ReturnCode

@testset "ILU-GMRES Linear Solver" begin

    # ========================================================================
    # Helper: Create a simple ALE-compatible test problem
    # ========================================================================
    function create_test_problem(; n_cells = 20, t_final = 5.0)
        virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

        material = Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))

        mesh = create_uniform_mesh(
            0.006,
            n_cells,
            Val(2);
            T_initial = 300.0,
            ξ_initial = (1190.0, 0.0),
        )

        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(50000.0) + RadiativeBC(ε = 0.86, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = ConvectiveMassBC(h_m = 0.01, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        )

        return PyrolysisProblem(
            mesh = mesh,
            material = material,
            bc_set = bc_set,
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
        sol = solve(
            problem;
            integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
            use_ale = true,
            # The exact ALE Jacobian (DirectSolve assembles the full pattern).
            jacobian = Structured(strategy = DirectSolve()),
            min_thickness = 0.1e-3,
            abstol = 1e-8,
            reltol = 1e-6,
            progress = false,
        )

        @test sol.retcode == :Success
        @test length(sol.t) > 1
        @test all(isfinite, sol.surface_T)
        @test sol.surface_T[end] > 300.0  # Should heat up
    end

    # ========================================================================
    @testset "Solution matches KLU baseline" begin
        problem_klu = create_test_problem()
        sol_klu = solve(
            problem_klu;
            integrator = KenCarp4(linsolve = KLUFactorization()),
            use_ale = true,
            # The exact ALE Jacobian (DirectSolve assembles the full pattern).
            jacobian = Structured(strategy = DirectSolve()),
            min_thickness = 0.1e-3,
            abstol = 1e-8,
            reltol = 1e-6,
            progress = false,
        )

        problem_ilu = create_test_problem()
        sol_ilu = solve(
            problem_ilu;
            integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
            use_ale = true,
            # The exact ALE Jacobian (DirectSolve assembles the full pattern).
            jacobian = Structured(strategy = DirectSolve()),
            min_thickness = 0.1e-3,
            abstol = 1e-8,
            reltol = 1e-6,
            progress = false,
        )

        @test sol_klu.retcode == :Success
        @test sol_ilu.retcode == :Success

        # Surface temperatures should match within 1 K
        T_diff = maximum(abs.(sol_klu.surface_T .- sol_ilu.surface_T))
        @test T_diff < 1.0
    end

    # ========================================================================
    @testset "Parameter variations" begin
        # τ = 0.01 (tighter ILU)
        problem1 = create_test_problem()
        sol1 = solve(
            problem1;
            integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.01)),
            use_ale = true,
            # The exact ALE Jacobian (DirectSolve assembles the full pattern).
            jacobian = Structured(strategy = DirectSolve()),
            min_thickness = 0.1e-3,
            progress = false,
        )
        @test sol1.retcode == :Success

        # τ = 0.5 (sparser ILU)
        problem2 = create_test_problem()
        sol2 = solve(
            problem2;
            integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.5)),
            use_ale = true,
            # The exact ALE Jacobian (DirectSolve assembles the full pattern).
            jacobian = Structured(strategy = DirectSolve()),
            min_thickness = 0.1e-3,
            progress = false,
        )
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

    @testset "LinearSolve algorithm interface" begin
        @test isempty(LinearSolve.algorithm_interface_issues(ILUGMRESFactorization()))
        @test isempty(LinearSolve.algorithm_interface_issues(ILUBiCGSTABFactorization()))
    end

    @testset "an unconverged Krylov solve is a failure" begin
        A = sparse([4.0 1.0 0.0 0.5; 1.0 4.0 1.0 0.0; 0.0 1.0 4.0 1.0; 0.5 0.0 1.0 4.0])
        b = [1.0, 2.0, 3.0, 4.0]
        for alg in (ILUGMRESFactorization(), ILUBiCGSTABFactorization())
            cache = LinearSolve.init(LinearProblem(A, b), alg)
            sol = LinearSolve.solve!(cache)
            @test sol.retcode == ReturnCode.Success
            @test A * sol.u ≈ b
        end
        # One GMRES iteration without preconditioning cannot reach 1e-14 on a
        # 4×4 system: the solution reports the failure the Newton acts on.
        alg = ILUGMRESFactorization()
        cache = LinearSolve.init(LinearProblem(A, b), alg)
        ws = cache.cacheval.krylov_solver
        Krylov.gmres!(ws, A, b; itmax = 1, rtol = 1e-14, atol = 0.0)
        @test !Krylov.issolved(ws)
        sol = _krylov_solution(alg, cache, ws)
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.iters == 1
    end

end
