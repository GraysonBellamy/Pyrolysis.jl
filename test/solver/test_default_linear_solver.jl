"""
The default Jacobian is `DirectSolve` on a fixed mesh and `ApproxSolve` under
ALE. The default linear solver follows the Jacobian pattern: UMFPACK when
`DirectSolve` materializes the ALE geometry block, KLU for every other sparse
Jacobian (the Eulerian Beer–Lambert suffix included). An explicit
`linear_solver` always wins.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: SURFACE_ABSORPTION, BEER_LAMBERT
using Pyrolysis.Solver: _resolve_integrator, _default_jacobian_backend
using LinearSolve: KLUFactorization, UMFPACKFactorization, SparspakFactorization

function _default_linsolve_problem()
    mat = Material(
        name = :default_linsolve_fixture,
        components = (
            SolidComponent(:s, ρ = 600.0, c = 1500.0, k = 0.15, α = 50.0, λ = 1e-6),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (Reaction(:r, 1 => 2, A = 1e10, E = 150e3, h = 300e3),),
    )
    mesh = create_uniform_mesh(0.005, 8, Val(2))
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 20e3) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, 0.5),
        output_times = [0.0, 0.5],
    )
end

_linsolve(sol) = sol.raw_solution.alg.linsolve

@testset "default Jacobian: DirectSolve on a fixed mesh, ApproxSolve under ALE" begin
    @test _default_jacobian_backend(false).strategy isa DirectSolve
    @test _default_jacobian_backend(true).strategy isa ApproxSolve
end

@testset "default linear solver follows the Jacobian pattern" begin
    @test _resolve_integrator(integrator = nothing, linear_solver = nothing).linsolve isa
          KLUFactorization
    @test _resolve_integrator(
        integrator = nothing,
        linear_solver = nothing,
        geometry_block = true,
    ).linsolve isa UMFPACKFactorization

    for (kw, expected) in (
        ((use_ale = false, radiation_model = SURFACE_ABSORPTION), KLUFactorization),
        ((use_ale = false, radiation_model = BEER_LAMBERT), KLUFactorization),
        ((use_ale = true, radiation_model = SURFACE_ABSORPTION), KLUFactorization),
        (
            (use_ale = true, radiation_model = SURFACE_ABSORPTION, handle_depletion = true),
            KLUFactorization,
        ),
        (
            (
                use_ale = true,
                radiation_model = SURFACE_ABSORPTION,
                jacobian = Structured(strategy = DirectSolve()),
            ),
            UMFPACKFactorization,
        ),
        (
            (
                use_ale = true,
                radiation_model = BEER_LAMBERT,
                jacobian = Structured(strategy = ApproxSolve()),
            ),
            KLUFactorization,
        ),
    )
        @testset "$(kw)" begin
            sol = solve(_default_linsolve_problem(); progress = false, kw...)
            @test sol.retcode == :Success
            @test _linsolve(sol) isa expected
        end
    end

    # An explicit choice is kept.
    sol = solve(
        _default_linsolve_problem();
        progress = false,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        linear_solver = SparspakFactorization(),
    )
    @test sol.retcode == :Success
    @test _linsolve(sol) isa SparspakFactorization
end
