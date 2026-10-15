"""
The integrator's type does not depend on the problem's. Two decks that differ
in component count, property types, reactions, boundary conditions, radiation
model and gas transport run, under the same solver options, on one integrator
type. The integrator, its nonlinear solver and the factorization therefore
compile once per solver configuration, and the precompile workloads cache
them for every problem.
"""

using Test
using Pyrolysis
using OrdinaryDiffEq: DiscreteCallback, derivative_discontinuity!

# Two constant-property components, a prescribed heat flux, Fickian transport.
function _plain_deck()
    mat = Material(
        name = :plain_deck,
        components = (
            SolidComponent(:resin, ρ = 1300.0, c = 1700.0, k = 0.20, ε = 0.9, λ = 2e-5),
            GaseousComponent(:vapor, M = 0.080, c = 1400.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (Reaction(:cure, 1 => 2, A = 1e12, E = 1.6e5, h = 300e3),),
    )
    bc = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(40e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(5e-3, 8, Val(2)),
        material = mat,
        bc_set = bc,
        tspan = (0.0, 5.0),
        output_times = [0.0, 5.0],
    )
end

# Three components with temperature-dependent closures, a charring reaction,
# a permeable solid (Darcy transport) and the cone boundary conditions.
function _rich_deck()
    k_virgin = CallableProperty(T -> 0.15 * (T / 300.0)^0.5)
    c_virgin = CallableProperty(T -> 1500.0 + 0.5 * (T - 300.0))
    mat = Material(
        name = :rich_deck,
        components = (
            SolidComponent(
                :virgin,
                ρ = 500.0,
                c = c_virgin,
                k = k_virgin,
                ε = 0.9,
                α = 5.0,
                κ = 1e-14,
                ρ_intrinsic = 1500.0,
                λ = 1e-5,
            ),
            SolidComponent(
                :char,
                ρ = 150.0,
                c = 1100.0,
                k = 0.08,
                ε = 0.95,
                α = 5.0,
                κ = 1e-12,
                ρ_intrinsic = 1500.0,
                λ = 1e-5,
            ),
            GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5),
        ),
        reactions = (
            Reaction(
                :pyrolysis,
                1 => (2, 3);
                A = 1e10,
                E = 150e3,
                h = 300e3,
                yields = (0.25, 0.75),
            ),
        ),
    )
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = RadiativeBC(ε = 0.9, T_inf = 300.0),
        ),
        mass = (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 8, Val(3)),
        material = mat,
        bc_set = bc,
        tspan = (0.0, 5.0),
        output_times = [0.0, 5.0],
    )
end

# A callback that records the integrator's type at the first step. Every solve
# gets a callback of the same type, so it does not tell the decks apart.
function _integrator_type_probe()
    seen = DataType[]
    probe = DiscreteCallback(
        (u, t, integrator) -> isempty(seen),
        integrator -> begin
            push!(seen, typeof(integrator))
            derivative_discontinuity!(integrator, false)
        end;
        save_positions = (false, false),
    )
    return probe, seen
end

function _solve_types(problem; kw...)
    probe, seen = _integrator_type_probe()
    sol = solve(problem; progress = false, callback = probe, kw...)
    @test sol.retcode == :Success
    @test length(seen) == 1
    return only(seen), typeof(sol.raw_solution)
end

@testset "integrator type is independent of the problem" begin
    configurations = [
        "fixed mesh" => (plain = (;), rich = (radiation_model = BEER_LAMBERT,)),
        "ALE" => (
            plain = (use_ale = true,),
            rich = (use_ale = true, radiation_model = BEER_LAMBERT),
        ),
        "ALE, parking" => (
            plain = (use_ale = true, handle_depletion = true),
            rich = (
                use_ale = true,
                radiation_model = BEER_LAMBERT,
                handle_depletion = true,
            ),
        ),
    ]
    for (name, kw) in configurations
        @testset "$name" begin
            integrator_plain, solution_plain = _solve_types(_plain_deck(); kw.plain...)
            integrator_rich, solution_rich = _solve_types(_rich_deck(); kw.rich...)
            @test integrator_plain === integrator_rich
            @test solution_plain === solution_rich
        end
    end
end
