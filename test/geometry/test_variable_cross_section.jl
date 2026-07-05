using Test
using Pyrolysis
using Pyrolysis.Problem: chi_range, chi_mode, WithChi

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "Variable Cross Section" begin
    @testset "WithChi Layout" begin
        problem = fixture_problem(n_cells = 5, lateral = true)
        def, ws, u, du = fixture_state(problem; use_ale = true)

        @test chi_mode(def.mode) isa WithChi
        @test length(chi_range(ws.layout)) == ws.layout.n_cells
        @test length(u) == 5 * (2 + 3) + 1

        u[chi_range(ws.layout)] .= 0.5
        A0 = ws.mesh.cross_section_area_initial
        A = cross_section_area_from_state(ws.mesh, problem.material, u, ws.layout)
        @test A ≈ A0 * 0.9

        residual!(du, u, def, 0.0, ws)
        @test all(isfinite, du)
        @test all(isfinite, du[chi_range(ws.layout)])
    end

    @testset "NoChi Layout Stays Compact" begin
        problem = fixture_problem(n_cells = 5, lateral = false)
        def, ws, u, _ = fixture_state(problem; use_ale = true)

        @test isempty(chi_range(ws.layout))
        @test length(u) == 5 * 3 + 6
        @test !(chi_mode(def.mode) isa WithChi)
    end

    @testset "Variable Cross Section Solve" begin
        problem = fixture_problem(n_cells = 5, lateral = true, tspan = (0.0, 0.02))
        sol = solve(
            problem;
            use_ale = true,
            progress = false,
            validate = false,
            # Opt out of Jacobian assembly — integrator falls back to FD.
            jacobian = nothing,
        )
        @test sol.retcode == :Success
        @test all(isfinite, sol.extras.cross_section_area_history)
    end
end
