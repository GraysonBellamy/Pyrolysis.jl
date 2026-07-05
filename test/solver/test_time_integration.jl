using Test
using SparseArrays
using OrdinaryDiffEq
using Pyrolysis
using Pyrolysis: JacobianSpec, Structured, DirectSolve, compute_jacobian!
using Pyrolysis.Problem: temperature_range, concentration_range, z_range, chi_range

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "Unified Time Integration" begin
    @testset "StateLayout Blocks" begin
        problem = fixture_problem(n_cells = 4)

        def_e, ws_e, u_e, _ = fixture_state(problem)
        @test state_length(ws_e.layout) == 4 * 3
        @test temperature_range(ws_e.layout) == 1:4
        @test concentration_range(ws_e.layout, 2) == 9:12
        @test isempty(z_range(ws_e.layout))
        @test isempty(chi_range(ws_e.layout))
        @test length(u_e) == state_length(ws_e.layout)

        def_a, ws_a, u_a, _ = fixture_state(problem; use_ale = true)
        @test state_length(ws_a.layout) == 4 * 3 + 5
        @test z_range(ws_a.layout) == 13:17
        @test isempty(chi_range(ws_a.layout))
        @test get_node_positions(u_a, ws_a.layout) == [n.z for n in ws_a.mesh.nodes]
        @test def_e.mode !== def_a.mode
    end

    @testset "Residual Evaluation" begin
        problem = fixture_problem(n_cells = 5)

        def, ws, u, du = fixture_state(problem)
        residual!(du, u, def, 0.0, ws)
        @test all(isfinite, du)
        @test any(!iszero, du)

        def_ale, ws_ale, u_ale, du_ale = fixture_state(problem; use_ale = true)
        residual!(du_ale, u_ale, def_ale, 0.0, ws_ale)
        @test all(isfinite, du_ale)
        @test length(du_ale) == state_length(ws_ale.layout)
    end

    @testset "Structured Jacobians" begin
        problem = fixture_problem(n_cells = 4)
        def, ws, u, _ = fixture_state(problem)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        S = spec.sparsity_phys
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

        def_a, ws_a, u_a, _ = fixture_state(problem; use_ale = true)
        spec_a = JacobianSpec(Structured(strategy = DirectSolve()), def_a, ws_a)
        S_a = spec_a.sparsity_phys .| spec_a.sparsity_geom
        J_a = SparseMatrixCSC{Float64,Int}(
            S_a.m,
            S_a.n,
            copy(S_a.colptr),
            copy(S_a.rowval),
            zeros(Float64, length(S_a.nzval)),
        )
        compute_jacobian!(J_a, spec_a, u_a, def_a, 0.0, ws_a)
        @test size(J_a) == (length(u_a), length(u_a))
        @test all(isfinite, nonzeros(J_a))
    end
end
