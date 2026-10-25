using Test, Pyrolysis
using Pyrolysis.Materials: R_GAS

# Gas is a sealed, calorically passive tracer. Pore volume changes its pressure,
# but the condensed chemistry, heat equation and deformation have no φ term.
# This isolates the assumption used by chemistry/geometry-only fixtures.
function incidental_pore_problem(φ)
    mat=Material(
        name = :incidental_gas,
        components = (
            SolidComponent(:A; ρ = 500.0, φ, c = 1000.0, k = 1.0, λ = 0.0),
            SolidComponent(:B; ρ = 1000.0, φ, c = 1000.0, k = 1.0, λ = 0.0),
            GaseousComponent(:tracer; M = 0.03, c = 0.0, k = 0.03),
        ),
        reactions = (Reaction(:compact, 1=>2; A = 0.5, E = 0.0, h = 0.0),),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 3),
        material = mat,
        T_initial = z->500.0,
        ξ_initial = [z->500.0, z->0.0, z->φ*1e5*0.03/(R_GAS*500.0)],
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(100.0), bottom = AdiabaticBC()),
        ),
        tspan = (0.0, 3.0),
        output_times = collect(0.0:0.25:3.0),
    )
end
@testset "porosity independence of incidental gas ($moving, $energy)" for moving in
                                                                          (false, true),
    energy in (:temperature, :enthalpy)

    solutions=map((0.1, 0.3, 0.6)) do φ
        sol=solve(
            incidental_pore_problem(φ);
            use_ale = moving,
            energy_state = energy,
            reltol = 1e-9,
            abstol = 1e-11,
            pressure_tolerance = 0.001,
            dtmin = 1e-15,
        )
        @test completed(sol)
        @test all(get_porosity(sol) .> 0)
        @test sum(sol.m[1, :, end])≈5exp(-1.5) rtol=1e-7
        @test sol.m[3, :, end]≈sol.m[3, :, 1] rtol=1e-10
        sol
    end
    reference=solutions[2]
    for sol in solutions
        @test maximum(abs.(sol.T .- reference.T))<1e-6
        @test sol.m[1:2, :, :]≈reference.m[1:2, :, :] rtol=1e-7 atol=1e-10
        @test sol.z≈reference.z rtol=1e-7 atol=1e-12
    end
end
