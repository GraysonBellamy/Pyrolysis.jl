"""
Regression test: when the solve has no sparse Jacobian prototype (the
`handle_depletion` colored-FD path), an incompatible sparse linear solver
must be swapped for a dense LU **inside the user's algorithm**, with a
warning — pre-fix both fallback sites silently replaced the user's entire
integrator (type, linsolve, autodiff) with a literal
`KenCarp4(AutoFiniteDiff(), LUFactorization())`.
"""

using Test
using Pyrolysis
using LinearSolve: SparspakFactorization
using OrdinaryDiffEq: Rodas5P
using Pyrolysis.Solver: AutoFiniteDiff

function _fallback_problem()
    mat = Material(
        name = :fallback_fixture,
        components = (
            SolidComponent(:s, ρ = 600.0, c = 1500.0, k = 0.15),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )
    mesh = create_uniform_mesh(0.005, 5, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(5e3), bottom = AdiabaticBC()),
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

@testset "incompatible sparse linsolve warns, integrator survives" begin
    # handle_depletion forces jac_spec === nothing (colored-FD, no sparse
    # prototype) — the first fallback site. An explicit sparse linear_solver must
    # warn and run with dense LU inside the USER'S constructor (Rodas5P),
    # which the pre-fix code silently replaced with KenCarp4.
    sol =
        @test_logs (:warn, r"requires a sparse Jacobian prototype") match_mode = :any solve(
            _fallback_problem();
            use_ale = true,
            handle_depletion = true,
            integrator = Rodas5P,
            linear_solver = SparspakFactorization(),
        )
    @test sol.retcode == :Success

    # A fully constructed user algorithm with a sparse linsolve is REBUILT
    # (type + autodiff preserved), not replaced.
    alg = Rodas5P(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
    sol2 = @test_logs (:warn, r"Rebuilding Rodas5P") match_mode = :any solve(
        _fallback_problem();
        use_ale = true,
        handle_depletion = true,
        integrator = alg,
    )
    @test sol2.retcode == :Success

    # Default path under depletion: the resolver now picks the dense LU
    # directly — no warning, no silent integrator swap needed.
    sol3 = solve(_fallback_problem(); use_ale = true, handle_depletion = true)
    @test sol3.retcode == :Success
end
