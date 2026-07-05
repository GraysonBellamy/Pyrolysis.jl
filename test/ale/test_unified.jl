using Test
using Pyrolysis

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "ALE Unified Path" begin
    @testset "State Carries Mesh Positions" begin
        problem = fixture_problem(n_cells = 5)
        def, ws, u, du = fixture_state(problem; use_ale = true)

        @test length(get_node_positions(u, ws.layout)) == length(ws.mesh.nodes)
        @test get_thickness_from_state(u, ws.layout) ≈ 0.01

        residual!(du, u, def, 0.0, ws)
        zrng = Pyrolysis.Problem.z_range(ws.layout)
        @test all(isfinite, du[zrng])
    end

    @testset "ALE Solve (no Jacobian)" begin
        problem = fixture_problem(n_cells = 5, tspan = (0.0, 0.02))
        sol = solve(
            problem;
            use_ale = true,
            progress = false,
            validate = false,
            # Opt out of Jacobian assembly — integrator falls back to its
            # own colored FD.
            jacobian = nothing,
        )
        @test sol.retcode == :Success
        @test length(sol.t) >= 2
        @test all(isfinite, sol.extras.thickness_history)
    end

    @testset "ALE Structured Solve" begin
        problem = fixture_problem(n_cells = 5, tspan = (0.0, 0.02))
        sol = solve(
            problem;
            use_ale = true,
            progress = false,
            validate = false,
            jacobian = JacobianSpec(
                Structured(strategy = DirectSolve()),
                problem;
                use_ale = true,
            ),
        )
        @test sol.retcode == :Success
    end
end
