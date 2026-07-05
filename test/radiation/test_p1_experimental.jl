using Test
using Pyrolysis
using Pyrolysis.Experimental:
    RadiationCache,
    P1Residual,
    compute_p1_residual!,
    default_radiation_bcs,
    get_radiation_intensity,
    P1_QUASI_STEADY
using Pyrolysis.Problem: state_length

include(joinpath(@__DIR__, "..", "helpers.jl"))

# An absorbing two-component material so the quasi-steady G actually couples
# to the state (the shared fixture material has no absorption coefficients).
function p1_material()
    solid =
        SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, γ = 1.0, α = 5.0)
    gas = GaseousComponent(:gas; M = 0.050, c = 1000.0, k = 0.02, λ = 1e-5, γ = 1.0)
    rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
    return Material(name = :p1_test, components = (solid, gas), reactions = (rxn,))
end

function p1_problem(top_thermal_bc; n_cells::Int = 6)
    mesh = create_uniform_mesh(
        0.01,
        n_cells,
        Val(2);
        T_initial = 300.0,
        ξ_initial = (1000.0, 0.0),
    )
    bcs = BoundaryConditionSet(
        thermal = (top = top_thermal_bc, bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(mesh, p1_material(), bcs, (0.0, 0.02))
    Pyrolysis.Problem.apply_initial_conditions!(problem)
    return problem
end

@testset "Experimental P1 Residual (unified layout)" begin
    @testset "G is not a state variable" begin
        # Regression: G lives in the RadiationCache, the state vector is
        # the ordinary unified layout, and the quasi-steady solve runs on
        # every residual evaluation — the field is never frozen.
        n = 6
        problem = p1_problem(ConvectiveBC(h = 10.0, T_inf = 300.0); n_cells = n)
        def, ws, u, du =
            fixture_state(problem; use_ale = true, radiation_model = P1_QUASI_STEADY)

        @test length(u) == state_length(ws.layout)  # T + ξ + z only, no G block

        rad_cache = RadiationCache(n)
        p1 = P1Residual(def, ws, rad_cache, default_radiation_bcs(2.5e4))

        p1(du, u, nothing, 0.0)
        @test all(isfinite, du)
        G1 = copy(rad_cache.G)
        @test all(>(0.0), G1)
        @test get_radiation_intensity(p1, 1) == rad_cache.G[1]

        # Heat the state: the next evaluation must re-solve G from the new
        # temperatures (the old state-vector G block stayed at its initial
        # value forever).
        u2 = copy(u)
        u2[1:n] .+= 300.0
        p1(du, u2, nothing, 0.0)
        @test any(i -> abs(rad_cache.G[i] - G1[i]) > 1e-3, 1:n)
    end

    @testset "Constructor requires the P1 radiation model" begin
        # Without radiation_model = P1_QUASI_STEADY in def.options the
        # core RHS would treat RadiativeFluxBC as a boundary heat flux while
        # the Marshak BC injects the same external flux — a double count.
        problem = p1_problem(ConvectiveBC(h = 10.0, T_inf = 300.0))
        def, ws, _, _ = fixture_state(problem; use_ale = true)  # NO_RADIATION
        rad_cache = RadiationCache(6)
        @test_throws ArgumentError P1Residual(
            def,
            ws,
            rad_cache,
            default_radiation_bcs(2.5e4),
        )
    end

    @testset "Boundary radiation BCs are inert under P1" begin
        # Regression: a cone-style stack RadiativeFluxBC + RadiativeBC
        # on top of the convective BC must contribute exactly nothing — the
        # external flux enters through the Marshak BC and emission/ambient
        # exchange is volumetric. du must be bit-identical with and without
        # the radiative BCs in the stack.
        bc_conv = ConvectiveBC(h = 10.0, T_inf = 300.0)
        bc_stack =
            RadiativeFluxBC(flux = 2.5e4) +
            RadiativeBC(T_inf = 300.0) +
            ConvectiveBC(h = 10.0, T_inf = 300.0)

        rad_bcs = default_radiation_bcs(2.5e4)

        results = map((bc_conv, bc_stack)) do top_bc
            problem = p1_problem(top_bc)
            def, ws, u, du = fixture_state(
                problem;
                use_ale = true,
                radiation_model = P1_QUASI_STEADY,
            )
            rad_cache = RadiationCache(6)
            p1 = P1Residual(def, ws, rad_cache, rad_bcs)
            p1(du, u, nothing, 0.0)
            du
        end

        @test results[1] == results[2]
    end

    @testset "P1 source enters the energy rows" begin
        # With a hot Marshak driving field the surface cells must heat faster
        # than the same problem with a cold (equilibrium) radiation field.
        problem = p1_problem(AdiabaticBC())
        def, ws, u, _ =
            fixture_state(problem; use_ale = true, radiation_model = P1_QUASI_STEADY)

        du_hot = similar(u)
        rad_hot = RadiationCache(6)
        P1Residual(def, ws, rad_hot, default_radiation_bcs(5e4))(du_hot, u, nothing, 0.0)

        du_cold = similar(u)
        rad_cold = RadiationCache(6)
        P1Residual(def, ws, rad_cold, default_radiation_bcs(0.0))(du_cold, u, nothing, 0.0)

        @test du_hot[6] > du_cold[6]  # top (surface) cell
        @test all(isfinite, du_hot) && all(isfinite, du_cold)
    end
end
