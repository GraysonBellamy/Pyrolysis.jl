"""
The `jacobian` keyword of `solve`: a bare backend is built for the solve's own
mode, a prebuilt spec is size-checked against the state vector, and
`state_length(spec)` reports the length a spec was built for.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian:
    JacobianSpec, Structured, ApproxSolve, DirectSolve, DenseAD_ForwardDiff
using Pyrolysis.Problem: state_length

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "jacobian kwarg" begin
    problem = fixture_problem(n_cells = 6, tspan = (0.0, 0.01))

    @testset "state_length(spec) for both cache families" begin
        def_e, ws_e, u_e = initialize_state(problem; use_ale = false)
        def_a, ws_a, u_a = initialize_state(problem; use_ale = true)
        @test state_length(JacobianSpec(Structured(), def_e, ws_e)) == length(u_e)
        @test state_length(
            JacobianSpec(Structured(strategy = ApproxSolve()), def_a, ws_a),
        ) == length(u_a)
        @test state_length(JacobianSpec(DenseAD_ForwardDiff(), def_a, ws_a)) == length(u_a)
        @test length(u_a) != length(u_e)
    end

    @testset "Eulerian-built spec fed to an ALE solve is rejected at entry" begin
        spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem)  # use_ale = false
        err = try
            solve(problem; use_ale = true, jacobian = spec, progress = false)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("use_ale", err.msg)
        @test occursin("$(state_length(spec))", err.msg)
    end

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

    @testset "matching prebuilt spec is accepted" begin
        spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
        sol = solve(problem; use_ale = true, jacobian = spec, progress = false)
        @test sol.retcode == :Success
    end
end
