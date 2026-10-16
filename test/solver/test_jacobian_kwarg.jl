"""
The `jacobian` keyword of `solve`: a `Structured` backend is built for the
solve's own mesh copy and mode; nothing else is accepted.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian:
    JacobianSpec, Structured, ApproxSolve, DirectSolve, DenseAD_ForwardDiff

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "jacobian kwarg" begin
    problem = fixture_problem(n_cells = 6, tspan = (0.0, 0.01))

    @testset "bare backend builds the spec for the solve's mode" begin
        for use_ale in (false, true)
            sol = solve(
                problem;
                use_ale = use_ale,
                jacobian = Structured(strategy = use_ale ? ApproxSolve() : DirectSolve()),
                progress = false,
            )
            @test sol.retcode == :Success
        end
    end

    @testset "prebuilt specs and the DenseAD reference are rejected" begin
        def, ws, _ = initialize_state(problem; use_ale = true)
        spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
        @test_throws TypeError solve(problem; use_ale = true, jacobian = spec)
        @test_throws TypeError solve(problem; jacobian = DenseAD_ForwardDiff())
    end
end
