"""The phase change of a liquid and its vapour: moisture content, equilibrium
pressure, rate and heat as pure functions of the cell state."""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: activity, sorption_heat, R_GAS
using Pyrolysis.Physics: moisture_content, phase_change_rate, phase_change_heat

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
