using Test
using SparseArrays
using Pyrolysis
using Pyrolysis: JacobianSpec, Structured, DirectSolve, ApproxSolve, compute_jacobian!
using Pyrolysis.Problem: chi_range

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "Jacobian Chi Coupling" begin
    @testset "Chi Columns Are Present For Lateral Law" begin
        problem = fixture_problem(n_cells = 4, lateral = true)
        def, ws, u, _ = fixture_state(problem; use_ale = true)
        u[chi_range(ws.layout)] .= range(0.1, 0.4, length = ws.layout.n_cells)

        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        S = spec.sparsity_phys .| spec.sparsity_geom
        J = SparseMatrixCSC{Float64,Int}(
            S.m,
            S.n,
            copy(S.colptr),
            copy(S.rowval),
            zeros(Float64, length(S.nzval)),
        )
        compute_jacobian!(J, spec, u, def, 0.0, ws)

        @test size(J) == (length(u), length(u))
        @test all(isfinite, nonzeros(J))
        @test any(!iszero, J[chi_range(ws.layout), chi_range(ws.layout)])
    end

    @testset "NoChi Sparsity Has No Chi Tail" begin
        problem = fixture_problem(n_cells = 4, lateral = false)
        _, ws, u, _ = fixture_state(problem; use_ale = true)
        @test isempty(chi_range(ws.layout))
    end
end

@testset "Reduced ALE Jacobian (Structured + ApproxSolve)" begin
    @testset "Reduced sparsity has nnz" begin
        problem = fixture_problem(n_cells = 5)
        def, ws, u, _ = fixture_state(problem; use_ale = true)
        spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
        @test size(spec.sparsity_phys, 1) == length(u)
        @test nnz(spec.sparsity_phys) > 0
    end

    @testset "Reduced Jacobian Evaluation" begin
        problem = fixture_problem(n_cells = 5)
        def, ws, u, _ = fixture_state(problem; use_ale = true)
        spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
        S = spec.sparsity_phys .| spec.sparsity_geom
        J = SparseMatrixCSC{Float64,Int}(
            S.m,
            S.n,
            copy(S.colptr),
            copy(S.rowval),
            zeros(Float64, length(S.nzval)),
        )
        compute_jacobian!(J, spec, u, def, 0.0, ws)
        @test size(J) == (length(u), length(u))
        @test all(isfinite, nonzeros(J))
    end

    @testset "Reduced Jacobian Solve" begin
        problem = fixture_problem(n_cells = 5, tspan = (0.0, 0.02))
        sol = solve(
            problem;
            use_ale = true,
            jacobian = Structured(strategy = ApproxSolve()),
            progress = false,
            validate = false,
        )
        @test sol.retcode == :Success
    end
end
