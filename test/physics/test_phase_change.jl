"""The phase change of a liquid and its vapour: moisture content, equilibrium
pressure, rate and heat as pure functions of the cell state."""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: activity, sorption_heat, R_GAS
using Pyrolysis.Physics: moisture_content, phase_change_rate, phase_change_heat
using Pyrolysis.Physics: pore_state
using Pyrolysis.Properties: skeletal_volume_fraction

const M_V = 0.018015       # molar mass of water vapour [kg/mol]
const L_REF = 2.4417e6     # latent heat of water at 298.15 K [J/kg]
const LAWS = (
    CAST_SPRUCE,
    CAST_BEECH,
    FreeLiquid(),
    ActivityFunction((X, T) -> 0.02 + 4X * exp(300 * (1 / 298.15 - 1 / T))),
)
const STATES =
    ((0.10, 350.0, 0.6), (0.03, 420.0, 0.7), (0.25, 300.0, 0.4), (0.06, 390.0, 0.65))

# The bulk vapour concentration in equilibrium with the liquid [kg/m³].
equilibrium_concentration(law, X, T, φ) =
    φ * M_V * equilibrium_vapour_pressure(law, WATER_SATURATION, X, T) / (R_GAS * T)

rate(law, K, X, T, φ, ξ_v) = phase_change_rate(law, WATER_SATURATION, M_V, K, X, T, φ, ξ_v)

function allocations(law, X, T, φ, ξ_v)
    rate(law, 1e3, X, T, φ, ξ_v)
    phase_change_heat(law, M_V, L_REF, X, T, 1.0e5, 2.0e5)
    return (@allocated rate(law, 1e3, X, T, φ, ξ_v)) +
           (@allocated phase_change_heat(law, M_V, L_REF, X, T, 1.0e5, 2.0e5)) +
           (@allocated moisture_content(ξ_v, 400.0))
end

@testset "Phase change: pure functions" begin
    @testset "moisture content" begin
        @test moisture_content(40.0, 400.0) == 40.0 / (400.0 + 1e-3)
        @test moisture_content(0.0, 400.0) == 0
        @test moisture_content(-4.0, 400.0) < 0
        # finite where the sorbing solids are gone: the liquid reads as free
        @test moisture_content(1.0, 0.0) == 1e3
        @test Pyrolysis.Physics.SORBENT_FLOOR == 1e-3
    end

    @testset "equilibrium vapour pressure" begin
        for law in LAWS, (X, T, φ) in STATES
            @test equilibrium_vapour_pressure(law, WATER_SATURATION, X, T) ==
                  activity(law, X, T) * WATER_SATURATION(T)
        end
        # free water boils at the saturation temperature to within X0/X
        p = equilibrium_vapour_pressure(FreeLiquid(), WATER_SATURATION, 0.25, 373.15)
        @test p / WATER_SATURATION(373.15) ≈ 0.25 / 0.2501 rtol = 1e-15
    end

    @testset "rate" begin
        for law in LAWS, (X, T, φ) in STATES, K in (1.0, 1e3, 1e4)
            ξ_eq = equilibrium_concentration(law, X, T, φ)
            @test ξ_eq > 0
            # zero at equilibrium
            @test abs(rate(law, K, X, T, φ, ξ_eq)) <= 1e-12 * K * ξ_eq
            # linear in the vapour, slope −K
            @test ForwardDiff.derivative(ξ -> rate(law, K, X, T, φ, ξ), 0.3ξ_eq) == -K
            @test rate(law, K, X, T, φ, 0.0) ≈ K * ξ_eq rtol = 1e-14
            @test rate(law, K, X, T, φ, 2ξ_eq) ≈ -K * ξ_eq rtol = 1e-14
            # evaporation below equilibrium, condensation above
            @test rate(law, K, X, T, φ, 0.5ξ_eq) > 0
            @test rate(law, K, X, T, φ, 1.5ξ_eq) < 0
            # a negative trial vapour mass raises the rate
            @test rate(law, K, X, T, φ, -1e-3) > rate(law, K, X, T, φ, 0.0) > 0
            # a negative trial water mass with no vapour draws vapour back
            @test rate(law, K, -0.01, T, φ, 0.0) < 0
        end
        # proportional to the pore fraction and the exchange rate at no vapour
        law, (X, T, φ) = CAST_SPRUCE, first(STATES)
        @test rate(law, 1e3, X, T, 2φ / 3, 0.0) ≈ 2 / 3 * rate(law, 1e3, X, T, φ, 0.0) rtol =
            1e-14
        @test rate(law, 2e3, X, T, φ, 0.0) ≈ 2 * rate(law, 1e3, X, T, φ, 0.0) rtol = 1e-14
        # the sample state by hand: K φ M a p_sat/(R T)
        a = activity(law, X, T)
        @test rate(law, 1e3, X, T, φ, 0.0) ≈
              1e3 * φ * M_V * a * WATER_SATURATION(T) / (R_GAS * T) rtol = 1e-14
    end

    @testset "heat" begin
        for law in LAWS, (X, T, φ) in STATES, (h_v, h_w) in ((0.0, 0.0), (1.0e5, 2.1e5))
            @test phase_change_heat(law, M_V, L_REF, X, T, h_v, h_w) ==
                  L_REF + h_v - h_w + sorption_heat(law, X, T, M_V)
        end
        # free water: Kirchhoff's latent heat alone
        @test phase_change_heat(FreeLiquid(), M_V, L_REF, 0.2, 373.15, 1.4e5, 3.15e5) ==
              L_REF + 1.4e5 - 3.15e5
        # the sorption heat adds to it below the cap
        @test phase_change_heat(CAST_SPRUCE, M_V, L_REF, 0.05, 300.0, 0.0, 0.0) >
              L_REF + 5e5
    end

    @testset "dual numbers, inference and allocations" begin
        J = Pyrolysis.ColoredJacobian.JacobianTag
        dX, dT, dξ = ForwardDiff.Dual{J}(0.1, 1.0),
        ForwardDiff.Dual{J}(350.0, 1.0),
        ForwardDiff.Dual{J}(0.05, 1.0)
        for law in LAWS
            g = ForwardDiff.gradient(
                v -> rate(law, 1e3, v[1], v[2], v[3], v[4]),
                [0.1, 350.0, 0.6, 0.05],
            )
            @test all(isfinite, g)
            h = ForwardDiff.gradient(
                v -> phase_change_heat(law, M_V, L_REF, v[1], v[2], v[3], v[4]),
                [0.1, 350.0, 1.0e5, 2.0e5],
            )
            @test all(isfinite, h)
            @test h[3] == 1 && h[4] == -1
            for (X, T, ξ) in
                ((0.1, 350.0, 0.05), (dX, 350.0, 0.05), (0.1, dT, 0.05), (0.1, 350.0, dξ))
                @inferred rate(law, 1e3, X, T, 0.6, ξ)
                @inferred phase_change_heat(law, M_V, L_REF, X, T, 1.0e5, 2.0e5)
                @test allocations(law, X, T, 0.6, ξ) == 0
            end
            @inferred moisture_content(dX, 400.0)
        end
    end
end

# Wood and char, both sorbing, with water, its vapour and air.
const WET = Material(
    name = :wet,
    components = (
        SolidComponent(
            :wood;
            ρ = 400.0,
            ρ_intrinsic = 1500.0,
            c = 1500.0,
            k = 0.12,
            λ = 5e-7,
        ),
        SolidComponent(
            :char;
            ρ = 100.0,
            ρ_intrinsic = 1800.0,
            c = 1100.0,
            k = 0.08,
            λ = 1e-6,
        ),
        LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6),
        GaseousComponent(:vapour; M = M_V, c = 1900.0, k = 0.025),
        GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026),
    ),
    phase_changes = (
        PhaseChange(
            :moisture,
            :water => :vapour;
            activity = CAST_SPRUCE,
            sorbents = (:wood, :char => 0.5),
            exchange_rate = 2e3,
        ),
    ),
)

@testset "Phase change: the resolved declaration and the start state" begin
    pc = only(WET.phase_changes)
    ξ = (300.0, 60.0, 30.0, 0.1, 0.9)
    ξ⁺ = (310.0, 50.0, 30.0, 0.1, 0.9)
    @testset "methods forward to the pure functions" begin
        # the sorbents weighted: S = 310 + 0.5·50 from ξ⁺, the water from ξ
        @test moisture_content(pc, ξ, ξ⁺) == moisture_content(30.0, 310.0 + 0.5 * 50.0)
        @test moisture_content(pc, (300.0, 60.0, -3.0, 0.1, 0.9), ξ⁺) < 0
        for (X, T, φ) in STATES, ξ_v in (0.0, 0.05, -1e-3)
            @test phase_change_rate(pc, X, T, φ, ξ_v) ===
                  phase_change_rate(CAST_SPRUCE, WATER_SATURATION, M_V, 2e3, X, T, φ, ξ_v)
            @test phase_change_heat(pc, X, T, 1.0e5, 2.0e5) ===
                  phase_change_heat(CAST_SPRUCE, M_V, L_REF, X, T, 1.0e5, 2.0e5)
            @test equilibrium_vapour_pressure(pc, X, T) ===
                  equilibrium_vapour_pressure(CAST_SPRUCE, WATER_SATURATION, X, T)
            @test equilibrium_vapour_concentration(pc, X, T, φ) ===
                  equilibrium_concentration(CAST_SPRUCE, X, T, φ)
            # the rate is K times the distance from the equilibrium concentration
            @test phase_change_rate(pc, X, T, φ, ξ_v) ≈
                  2e3 * (equilibrium_vapour_concentration(pc, X, T, φ) - ξ_v) rtol = 1e-15
        end
        J = Pyrolysis.ColoredJacobian.JacobianTag
        dξ = map(x -> ForwardDiff.Dual{J}(x, 1.0), ξ)
        @inferred moisture_content(pc, dξ, dξ)
        @test (@allocated moisture_content(pc, dξ, dξ)) == 0
        @inferred phase_change_rate(pc, dξ[3] / 400, 350.0, 0.6, dξ[4])
    end

    @testset "the start in equilibrium" begin
        for T in (280.0, 300.0, 350.0, 390.0), X in (0.02, 0.1, 0.25)
            ξc = (400.0, 50.0, X * 425.0, 0.0, 0.0)
            ξ_v = equilibrium_vapour_concentration(WET, ξc, T)
            p_eq = equilibrium_vapour_pressure(WET, ξc, T)
            φ = 1 - skeletal_volume_fraction(ξc, T, WET)
            Xc = moisture_content(pc, ξc, ξc)
            @test Xc ≈ X * 425.0 / (425.0 + 1e-3) rtol = 1e-15
            @test p_eq == equilibrium_vapour_pressure(pc, Xc, T)
            @test ξ_v == equilibrium_vapour_concentration(pc, Xc, T, φ)
            # the background gas at P − p_eq makes up the pressure
            P = max(101325.0, 2p_eq)
            ξ_air = background_gas_concentration(WET, :air, ξc, T; P = P - p_eq)
            state = pore_state((400.0, 50.0, ξc[3], ξ_v, ξ_air), T, WET)
            # the solver's pore state reads the vapour pressure back
            @test state.p[4] ≈ p_eq rtol = 4eps()
            @test state.P ≈ P rtol = 4eps()
            # the gas entries do not enter: the vapour first, the air second
            @test equilibrium_vapour_concentration(
                WET,
                (400.0, 50.0, ξc[3], ξ_v, ξ_air),
                T,
            ) == ξ_v
        end
        # by index or by name
        ξc = (400.0, 0.0, 40.0, 0.0, 0.0)
        @test equilibrium_vapour_concentration(WET, ξc, 300.0; phase_change = :moisture) ==
              equilibrium_vapour_concentration(WET, ξc, 300.0; phase_change = 1) ==
              equilibrium_vapour_concentration(WET, ξc, 300.0)
        @test equilibrium_vapour_pressure(WET, ξc, 300.0; phase_change = :moisture) ==
              equilibrium_vapour_pressure(WET, ξc, 300.0)
        # dry wood holds no vapour
        @test equilibrium_vapour_concentration(WET, (400.0, 0.0, 0.0, 0.0, 0.0), 300.0) == 0
        # refusals
        @test_throws ArgumentError equilibrium_vapour_concentration(
            WET,
            ξc,
            300.0;
            phase_change = 2,
        )
        @test_throws ArgumentError equilibrium_vapour_concentration(
            WET,
            ξc,
            300.0;
            phase_change = :drying,
        )
        @test_throws ArgumentError equilibrium_vapour_pressure(
            WET,
            ξc,
            300.0;
            phase_change = 0,
        )
        @test_throws ArgumentError equilibrium_vapour_concentration(WET, ξc, -1.0)
        @test_throws ArgumentError equilibrium_vapour_pressure(WET, ξc, NaN)
        @test_throws ArgumentError equilibrium_vapour_concentration(
            WET,
            (1500.0, 0.0, 0.0, 0.0, 0.0),
            300.0,
        )
        dry = Material(name = :dry, components = WET.components)
        @test_throws ArgumentError equilibrium_vapour_concentration(dry, ξc, 300.0)
    end
end
