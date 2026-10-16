"""
Tests for Darcy-Fick infrastructure.

Tests cover:
- Ideal gas pressure computation
- Component-based permeability (effective_permeability)
- Gas viscosity (Sutherland's law)
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Problem: StateCache
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Properties
using Pyrolysis.Physics
using Pyrolysis.Physics: compute_gas_convective_source!
using Pyrolysis.Geometry
using Pyrolysis.Discretization

@testset "Darcy-Fick Infrastructure" begin

    # ========================================================================
    # Test Setup: Create test materials
    # ========================================================================

    # Two-component material: solid + gas (with component-based permeability)
    solid_only = Material(
        name = :solid_only,
        components = (SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-12),),
        reactions = (),
    )

    two_comp = Material(
        name = :two_comp,
        components = (
            SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-11),
            GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )

    three_comp = Material(
        name = :three_comp,
        components = (
            SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-11),
            SolidComponent(:char, ρ = 200.0, c = 1000.0, k = 0.1, κ = 1e-10),
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),  # N2-like
        ),
        reactions = (),
    )

    # ========================================================================
    @testset "Ideal Gas Pressure" begin
        # Pressure now uses P = N/φ where N = Σ ξⱼ R T / Mⱼ over gas components
        # and φ is the porosity (gas volume fraction).

        # ====================================================================
        @testset "Single gas component" begin
            # Pure gas → φ = 1, so P = ξ_gas R T / M (matches old formula)
            ξ_gas = 1.0     # kg/m³
            T = 400.0       # K
            M = 0.044       # kg/mol (CO2)

            ξ = (0.0, ξ_gas)  # (solid, gas)
            φ = porosity(ξ, T, two_comp)
            @test φ ≈ 1.0 rtol = 1e-10

            P_expected = ξ_gas * R_GAS * T / M  # φ = 1, no divisor effect
            P = ideal_gas_pressure(ξ, T, φ, two_comp)
            @test P ≈ P_expected rtol = 1e-10
        end

        # ====================================================================
        @testset "No gas components" begin
            # Solid-only material should have zero pressure
            ξ = (1000.0,)
            T = 400.0
            φ = porosity(ξ, T, solid_only)

            P = ideal_gas_pressure(ξ, T, φ, solid_only)
            @test P == 0.0
        end

        # ====================================================================
        @testset "Multiple components — porosity divisor" begin
            # Two solids + 1 gas, total solid vol < 1 so φ is meaningful.
            # Solid1 ρ = 1000 → vol = 500/1000 = 0.5
            # Solid2 ρ = 200  → vol = 50/200  = 0.25
            # φ = 1 - 0.5 - 0.25 = 0.25
            ξ_solid1 = 500.0
            ξ_solid2 = 50.0
            ξ_gas = 0.5
            T = 500.0
            M_gas = 0.028

            ξ = (ξ_solid1, ξ_solid2, ξ_gas)
            φ = porosity(ξ, T, three_comp)
            @test φ ≈ 0.25 rtol = 1e-10

            # P = (ξ_gas R T / M_gas) / φ
            N_expected = ξ_gas * R_GAS * T / M_gas
            P_expected = N_expected / φ
            P = ideal_gas_pressure(ξ, T, φ, three_comp)
            @test P ≈ P_expected rtol = 1e-10
        end

        # ====================================================================
        @testset "Zero gas concentration" begin
            # No gas → numerator is zero regardless of porosity
            ξ = (1000.0, 0.0)
            T = 400.0
            φ = porosity(ξ, T, two_comp)

            P = ideal_gas_pressure(ξ, T, φ, two_comp)
            @test P == 0.0
        end

        # ====================================================================
        @testset "High temperature scaling" begin
            # At fixed ξ, P scales linearly with T (φ depends only on ξ here)
            ξ = (0.0, 1.0)
            T1 = 300.0
            T2 = 600.0
            φ1 = porosity(ξ, T1, two_comp)
            φ2 = porosity(ξ, T2, two_comp)

            P1 = ideal_gas_pressure(ξ, T1, φ1, two_comp)
            P2 = ideal_gas_pressure(ξ, T2, φ2, two_comp)

            @test P2 / P1 ≈ T2 / T1 rtol = 1e-10
        end
    end

    # ========================================================================
    @testset "Pressure Derivatives (ForwardDiff)" begin
        # The analytic `porosity_derivatives` / `pressure_derivatives` family was
        # deleted (dead API; it ignored the porosity clamp and blew up as
        # N/φ_safe² at the 1e-12 floor). The live Jacobian path differentiates
        # `porosity` ∘ `ideal_gas_pressure` with ForwardDiff — verify it against
        # finite differences, and verify it stays finite AT the clamp.
        @testset "AD derivatives match finite difference" begin
            ξ = (500.0, 1.0)
            T = 400.0
            ε = 1e-6

            function _P(ξ_local, T_local)
                φ_l = porosity(ξ_local, T_local, two_comp)
                return ideal_gas_pressure(ξ_local, T_local, φ_l, two_comp)
            end

            g = ForwardDiff.gradient(x -> _P((x[1], x[2]), x[3]), [ξ[1], ξ[2], T])

            @test g[3] ≈ (_P(ξ, T + ε) - _P(ξ, T - ε)) / (2ε) rtol = 1e-6
            @test g[1] ≈ (_P((ξ[1] + ε, ξ[2]), T) - _P((ξ[1] - ε, ξ[2]), T)) / (2ε) rtol =
                1e-6
            @test g[2] ≈ (_P((ξ[1], ξ[2] + ε), T) - _P((ξ[1], ξ[2] - ε), T)) / (2ε) rtol =
                1e-6
        end

        @testset "AD stays finite at the porosity clamp" begin
            # Over-packed solid: φ clamps to 0 → P = N/max(φ, 1e-12). The old
            # analytic path formed N/φ_safe² ≈ 1e24·N here; the live AD path
            # must return finite values and finite (clamp-consistent, zero
            # through the φ channel) derivatives.
            ξ = (2000.0, 1.0)   # 2000/1000 = 2 > 1 → φ clamped to 0
            T = 400.0
            φ = porosity(ξ, T, two_comp)
            @test φ == 0.0

            function _P(ξ_local, T_local)
                φ_l = porosity(ξ_local, T_local, two_comp)
                return ideal_gas_pressure(ξ_local, T_local, φ_l, two_comp)
            end

            g = ForwardDiff.gradient(x -> _P((x[1], x[2]), x[3]), [ξ[1], ξ[2], T])
            @test all(isfinite, g)
            # Clamped porosity: ∂φ/∂ξ_solid is 0 through the clamp, so the
            # solid-column derivative carries no N/φ² blow-up.
            @test abs(g[1]) < 1e30
        end
    end

    @testset "Component-Based Permeability" begin
        # ====================================================================
        @testset "Effective Permeability - Single Component" begin
            # Material with single solid with permeability
            mat_single = Material(
                name = :single_perm,
                components = (
                    SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-12),
                ),
                reactions = (),
            )

            T = 500.0
            ξ = (1000.0,)

            κ_eff = effective_permeability(ξ, T, mat_single)
            @test κ_eff ≈ 1e-12 rtol = 0.01
        end

        # ====================================================================
        @testset "Effective Permeability - Two Solids Mixing" begin
            # Material with two solids of different permeability
            mat_mix = Material(
                name = :mix_perm,
                components = (
                    SolidComponent(:virgin, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-15),
                    SolidComponent(:char, ρ = 200.0, c = 800.0, k = 0.1, κ = 1e-12),
                    GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
                ),
                reactions = (),
                mixing = MaterialMixing(κ = MixingSpec(PARALLEL)),
            )

            T = 500.0

            # Pure virgin: κ should be ~1e-15
            ξ_virgin = (1000.0, 0.0, 0.0)
            κ_virgin = effective_permeability(ξ_virgin, T, mat_mix)
            @test κ_virgin ≈ 1e-15 rtol = 0.1

            # Pure char: κ should be ~1e-12
            ξ_char = (0.0, 200.0, 0.0)
            κ_char = effective_permeability(ξ_char, T, mat_mix)
            @test κ_char ≈ 1e-12 rtol = 0.1

            # 50/50 mixture (by volume): κ between extremes
            ξ_mix = (500.0, 100.0, 0.0)
            κ_mix = effective_permeability(ξ_mix, T, mat_mix)
            @test κ_mix > κ_virgin
            @test κ_mix < κ_char
        end

        # ====================================================================
        @testset "Solid-basis mixing with impermeable phases" begin
            # Regression: κ mixing operates on the RENORMALIZED SOLID
            # basis (like conductivity), and an impermeable phase with v > 0
            # zeroes the series (Wiener lower) bound.
            T = 400.0
            # v_solid1 = 0.5, v_solid2 = 0.2, v_gas = 0.3
            ξ = (500.0, 200.0, 0.4)

            mk(rule) = Material(
                name = :c3_uniform,
                components = (
                    SolidComponent(:s1, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-13),
                    SolidComponent(:s2, ρ = 1000.0, c = 1000.0, k = 0.1, κ = 1e-13),
                    GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
                ),
                reactions = (),
                mixing = MaterialMixing(κ = MixingSpec(rule)),
            )

            # Uniform solid κ₀ with 30% gas: both bounds must equal κ₀ on the
            # solid basis. (The old all-phase basis gave κ_series = 1.43e-13 >
            # κ₀ > κ_parallel = 7e-14 — the "lower bound" above the upper.)
            @test effective_permeability(ξ, T, mk(PARALLEL)) ≈ 1e-13 rtol = 1e-12
            @test effective_permeability(ξ, T, mk(SERIES)) ≈ 1e-13 rtol = 1e-12

            mk2(rule) = Material(
                name = :c3_impermeable,
                components = (
                    SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2),  # κ = 0 default
                    SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.1, κ = 1e-12),
                    GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
                ),
                reactions = (),
                mixing = MaterialMixing(κ = MixingSpec(rule)),
            )

            # Default-κ virgin (impermeable) + char + gas: the series bound is
            # EXACTLY zero (gas cannot cross a series-stacked impermeable
            # layer); parallel is the solid-basis volume-weighted mean; the
            # default WEIGHTED β = 0.5 is their average. (The old code gave
            # κ_series = 5e-12 — five times char's own κ through a
            # 71%-impermeable matrix.)
            @test effective_permeability(ξ, T, mk2(SERIES)) == 0.0
            @test effective_permeability(ξ, T, mk2(PARALLEL)) ≈ (0.2 / 0.7) * 1e-12 rtol =
                1e-12
            @test effective_permeability(ξ, T, mk2(WEIGHTED)) ≈ 0.5 * (0.2 / 0.7) * 1e-12 rtol =
                1e-12

            # ForwardDiff through effective_permeability is exact on both
            # branches: finite differences across the WEIGHTED rule.
            let mat = mk2(WEIGHTED)
                ∂ξ_ = ForwardDiff.gradient(
                    x -> effective_permeability(Tuple(x), T, mat),
                    [ξ...],
                )
                h = 1e-3
                for i = 1:3
                    ξp = ntuple(j -> j == i ? ξ[j] + h : ξ[j], 3)
                    ξm = ntuple(j -> j == i ? ξ[j] - h : ξ[j], 3)
                    fd =
                        (
                            effective_permeability(ξp, T, mat) -
                            effective_permeability(ξm, T, mat)
                        ) / (2h)
                    @test isapprox(∂ξ_[i], fd; rtol = 1e-5, atol = 1e-30)
                end
            end
        end

        # ====================================================================
        @testset "Gases Do Not Contribute to Permeability" begin
            mat = Material(
                name = :gas_test,
                components = (
                    SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-12),
                    GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
                ),
                reactions = (),
            )

            # Gas components should not have permeability set (implicitly κ = 0)
            # Check component κ values - only solids should have non-zero κ
            gas_component = mat.components[2]
            # GaseousComponent should have κ = 0 (or not set at all, defaulting to 0)
            # This is verified by checking the effective permeability with pure gas

            # Pure solid with permeability
            ξ_pure_solid = (1000.0, 0.0)
            κ_solid = effective_permeability(ξ_pure_solid, 500.0, mat)
            @test κ_solid == 1e-12

            # Check that gas contribution is minimal to total permeability
            # (gases are effectively treated as non-porous voids)
            ξ_with_gas = (1000.0, 10.0)
            κ_mixed = effective_permeability(ξ_with_gas, 500.0, mat)
            # Permeability can change but should still be driven by solid phase
            @test κ_mixed > 0  # Still has permeability from solids
        end

        # ====================================================================
        @testset "has_darcy_flow Detection" begin
            # Material with permeability
            mat_with_κ = Material(
                name = :HasDarcy,
                components = (
                    SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-12),
                ),
                reactions = (),
            )
            @test has_darcy_flow(mat_with_κ) == true

            # Material without permeability (default κ=0)
            mat_no_κ = Material(
                name = :NoDarcy,
                components = (
                    SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 0.0),
                ),
                reactions = (),
            )
            @test has_darcy_flow(mat_no_κ) == false
        end
    end

    # ========================================================================
    @testset "Gas Viscosity" begin
        # ====================================================================
        @testset "Reference condition" begin
            # At T_ref = 273.15 K, should return the canonical Sutherland
            # air reference μ_ref = 1.716e-5 Pa·s
            μ = SutherlandProperty()(273.15)
            @test μ ≈ 1.716e-5 rtol = 1e-10

            # Air at 300 K: 1.846e-5 Pa·s (the old default 1.81e-5 was the
            # ≈292.5 K value, a −2% bias on every Darcy velocity)
            μ_300 = SutherlandProperty()(300.0)
            @test μ_300 ≈ 1.846e-5 rtol = 1e-3
        end

        # ====================================================================
        @testset "Temperature dependence" begin
            # Viscosity should increase with temperature (gas property)
            μ_300 = SutherlandProperty()(300.0)
            μ_600 = SutherlandProperty()(600.0)
            μ_1000 = SutherlandProperty()(1000.0)

            @test μ_600 > μ_300
            @test μ_1000 > μ_600
        end

        # ====================================================================
        @testset "Sutherland's law verification" begin
            # Verify against known values for air. The 1% tolerance is tight
            # enough to reject the previous constants (μ_ref = 1.81e-5 @ 300 K,
            # S = 110.0), which deviate ≈2.2% at these points.
            # At 400K, air viscosity ≈ 2.29e-5 Pa·s
            μ_400 = SutherlandProperty()(400.0)
            @test μ_400 ≈ 2.29e-5 rtol = 0.01

            # At 500K, air viscosity ≈ 2.67e-5 Pa·s
            μ_500 = SutherlandProperty()(500.0)
            @test μ_500 ≈ 2.67e-5 rtol = 0.01
        end

        # ====================================================================
        @testset "Custom parameters" begin
            # Test with custom reference values
            μ_custom = SutherlandProperty(2.0e-5, 350.0, 120.0)(400.0)
            @test μ_custom > 0.0
            @test μ_custom != SutherlandProperty()(400.0)  # Different from default
        end

        # ====================================================================
        @testset "Mixture viscosity" begin
            ξ = (500.0, 1.0)
            T = 400.0

            μ_mix = mixture_gas_viscosity(ξ, T, two_comp)

            # Reads the material's own curve; `two_comp` leaves it at the
            # SutherlandProperty() default, so this is the bare air value.
            @test μ_mix == effective_gas_viscosity(T, two_comp)
            @test μ_mix ≈ SutherlandProperty()(T) rtol = 1e-10

            # TRIPWIRE — μ is a function of T alone. The production Darcy
            # Jacobian is the ANALYTICAL override in
            # `Jacobian/operators/darcy_transport.jl`, whose mobility chain rule
            # assumes dμ/dξ == 0 and takes the ξ-path through κ only. If a
            # composition-dependent rule (Wilke, option B) lands here, this
            # assert is the thing that must fail loudly: the rtol = 1e-12
            # analytical-vs-AD gate cannot catch the missing dμ/dξ terms while
            # every Darcy fixture carries a single gas component, because a
            # one-species mixture rule is the identity. When you intend to break
            # this, also add the `− mob·dμ_dξ[j]` terms to every `dvg_dξ_*_j`
            # expression in that override and land a multi-gas fixture first.
            # See docs/Gas_Viscosity_Recommendation_2026-09-21.md §1 and §3.
            for ξ_alt in ((500.0, 0.0), (500.0, 50.0), (10.0, 1.0), (0.0, 1.0))
                @test mixture_gas_viscosity(ξ_alt, T, two_comp) == μ_mix
            end
        end

        # ====================================================================
        @testset "Nonphysical trial-state hardening (T floor)" begin
            # Diverging Newton trial states can present nonphysical (even
            # negative) temperatures; T^1.5 must not throw a DomainError —
            # it killed an entire ApproxSolve ALE solve near burnout
            # (2026-06-11). The kernels clamp at T_GAS_PROPERTY_FLOOR.
            μ_neg = SutherlandProperty()(-43.0)
            @test isfinite(μ_neg) && μ_neg > 0.0
            # The floor binds: any T below it gives the floor's value.
            @test μ_neg == SutherlandProperty()(0.0)
            # The floor never binds on physical states (≥ ~200 K).
            @test SutherlandProperty()(300.0) ≈ 1.846e-5 rtol = 1e-3

            D_neg = chapman_enskog_diffusivity(-43.0, 0.06)
            @test isfinite(D_neg) && D_neg > 0.0
            @test D_neg == chapman_enskog_diffusivity(0.0, 0.06)
        end
    end

    # ========================================================================
    @testset "Chapman-Enskog Diffusivity" begin
        # ====================================================================
        @testset "Collision integral - Neufeld correlation" begin
            # Test against tabulated values
            # At T* = 1.0, Ω_D ≈ 1.44
            Ω_1 = collision_integral_diffusion(1.0)
            @test Ω_1 ≈ 1.44 rtol = 0.02

            # At T* = 2.0, Ω_D ≈ 1.075
            Ω_2 = collision_integral_diffusion(2.0)
            @test Ω_2 ≈ 1.075 rtol = 0.02

            # At T* = 5.0, Ω_D ≈ 0.84
            Ω_5 = collision_integral_diffusion(5.0)
            @test Ω_5 ≈ 0.84 rtol = 0.03

            # Ω_D should decrease with increasing T*
            @test Ω_2 < Ω_1
            @test Ω_5 < Ω_2
        end

        # ====================================================================
        @testset "Temperature dependence" begin
            M = 0.044  # CO2
            # Diffusivity should increase with temperature
            D_300 = chapman_enskog_diffusivity(300.0, M)
            D_500 = chapman_enskog_diffusivity(500.0, M)
            D_800 = chapman_enskog_diffusivity(800.0, M)

            @test D_500 > D_300
            @test D_800 > D_500

            # Roughly D ∝ T^1.5 to T^1.75
            ratio = D_500 / D_300
            expected_ratio_low = (500.0 / 300.0)^1.5  # ≈ 2.15
            expected_ratio_high = (500.0 / 300.0)^1.75  # ≈ 2.55
            @test ratio > expected_ratio_low * 0.8  # Within range
            @test ratio < expected_ratio_high * 1.2
        end

        # ====================================================================
        @testset "Pressure dependence" begin
            M = 0.028  # N2
            T = 500.0
            # D ∝ 1/P
            D_1atm = chapman_enskog_diffusivity(T, M; P = 101325.0)
            D_2atm = chapman_enskog_diffusivity(T, M; P = 202650.0)

            @test D_1atm / D_2atm ≈ 2.0 rtol = 0.01
        end

        # ====================================================================
        @testset "Molecular weight dependence" begin
            T = 500.0
            # Lighter molecules diffuse faster
            D_H2O = chapman_enskog_diffusivity(T, 0.018)   # Water
            D_N2 = chapman_enskog_diffusivity(T, 0.028)    # Nitrogen
            D_CO2 = chapman_enskog_diffusivity(T, 0.044)   # CO2
            D_MMA = chapman_enskog_diffusivity(T, 0.100)   # MMA

            @test D_H2O > D_N2
            @test D_N2 > D_CO2
            @test D_CO2 > D_MMA
        end

        # ====================================================================
        @testset "Known values for common gases" begin
            # Compare to literature values (at 1 atm)
            # Values are approximate due to Lennard-Jones parameter estimation

            # H2O in air at 373 K ≈ 2.2e-5 m²/s
            D_H2O_373 = chapman_enskog_diffusivity(373.0, 0.018)
            @test D_H2O_373 > 1e-5  # Order of magnitude check
            @test D_H2O_373 < 5e-5

            # CO2 in air at 298 K ≈ 1.6e-5 m²/s
            D_CO2_298 = chapman_enskog_diffusivity(298.0, 0.044)
            @test D_CO2_298 > 5e-6
            @test D_CO2_298 < 5e-5

            # N2 self-diffusion at 300 K ≈ 2.0e-5 m²/s
            D_N2_300 = chapman_enskog_diffusivity(300.0, 0.028)
            @test D_N2_300 > 1e-5
            @test D_N2_300 < 5e-5
        end

        # ====================================================================
        @testset "Custom Lennard-Jones parameters" begin
            T = 500.0
            M = 0.016  # CH4

            # With default parameters
            D_default = chapman_enskog_diffusivity(T, M)

            # With explicit CH4 Lennard-Jones parameters
            D_explicit = chapman_enskog_diffusivity(T, M; σ = 3.758, ε_k = 148.6)

            # Both should be positive and reasonable
            @test D_default > 0.0
            @test D_explicit > 0.0
            @test D_explicit > 1e-5
            @test D_explicit < 1e-4
        end

        # ====================================================================
        # ====================================================================
        @testset "Physical constants" begin
            @test CHAPMAN_ENSKOG_PREFACTOR > 0.0
            @test SIGMA_AIR ≈ 3.617 rtol = 0.01
            @test EPSILON_K_AIR ≈ 97.0 rtol = 0.01
            @test M_AIR ≈ 0.029 rtol = 0.01
        end
    end

    # ========================================================================
    # Darcy-Fick Gas Transport Tests
    # ========================================================================

    @testset "Darcy-Fick Single Face Flux" begin
        # Test gas_flux_darcy_fick function (volume fraction gradient formulation)

        @testset "Pure diffusion (v_g = 0)" begin
            # When v_g = 0, should equal gas_flux (diffusion only)
            λ_face = 0.01   # m²/s
            ξ_L = 0.5       # kg/m³
            ξ_R = 1.0       # kg/m³
            T_L = 400.0     # K
            T_R = 400.0     # K
            d_LR = 0.001    # m
            v_g = 0.0       # No advection

            J_df = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)
            J_diff = gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)

            @test J_df ≈ J_diff rtol = 1e-10
        end

        @testset "Pure advection (λ = 0)" begin
            # When λ = 0, should be pure upwind advection
            λ_face = 0.0
            ξ_L = 0.5
            ξ_R = 1.0
            T_L = 400.0
            T_R = 400.0
            d_LR = 0.001

            # Positive velocity: upwind from left
            v_g = 0.1  # m/s
            J_pos = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)
            @test J_pos ≈ v_g * ξ_L rtol = 1e-10

            # Negative velocity: upwind from right
            v_g = -0.1
            J_neg = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)
            @test J_neg ≈ v_g * ξ_R rtol = 1e-10
        end

        @testset "Advection carries the intrinsic density ξ/φ_face" begin
            # ξ is per BULK volume; the Darcy velocity is superficial, so the
            # advected scalar is ξ/φ. φ_face is the arithmetic face mean.
            λ_face = 0.0
            ξ_L = 0.5
            ξ_R = 1.0
            T_L = T_R = 400.0
            d_LR = 0.001
            v_g = 0.1
            φ_L, φ_R = 0.2, 0.4   # φ_face = 0.3
            J = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)
            @test J ≈ v_g * ξ_L / 0.3 rtol = 1e-12

            # Reversed flow upwinds the right state, same φ_face
            J_neg = gas_flux_darcy_fick(-v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)
            @test J_neg ≈ -v_g * ξ_R / 0.3 rtol = 1e-12

            # φ = 0.5 doubles the advective flux relative to the old ξ·v form
            J_half = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, 0.5, 0.5, λ_face, T_L, T_R, d_LR)
            @test J_half ≈ 2.0 * v_g * ξ_L rtol = 1e-12
        end

        @testset "Combined advection and diffusion" begin
            λ_face = 0.01
            ξ_L = 0.5
            ξ_R = 1.0
            T_L = 400.0
            T_R = 400.0
            d_LR = 0.001
            v_g = 0.05  # m/s

            J_df = gas_flux_darcy_fick(v_g, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)

            # Manual computation (volume fraction gradient)
            T_face = 0.5 * (T_L + T_R)
            J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
            J_adv = v_g * ξ_L  # v_g > 0 so upwind is left
            J_expected = J_adv + J_diff

            @test J_df ≈ J_expected rtol = 1e-10
        end

        @testset "Peclet number behavior" begin
            # High Pe: advection dominates
            # Low Pe: diffusion dominates
            λ_face = 0.01
            ξ_L = 0.5
            ξ_R = 1.0
            T_L = 400.0
            T_R = 400.0
            d_LR = 0.001

            # Compute characteristic diffusion velocity from volume fraction gradient
            T_face = 0.5 * (T_L + T_R)
            J_diff_ref = abs(λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR))
            v_diff = J_diff_ref / max(ξ_L, ξ_R)

            # High Peclet (advection >> diffusion)
            v_high = 10.0 * v_diff
            J_high = gas_flux_darcy_fick(v_high, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)
            J_adv_only = v_high * ξ_L
            # Should be close to pure advection (within 30%)
            @test abs(J_high - J_adv_only) / abs(J_high) < 0.3

            # Low Peclet (diffusion >> advection)
            v_low = 0.01 * v_diff
            J_low = gas_flux_darcy_fick(v_low, ξ_L, ξ_R, 1.0, 1.0, λ_face, T_L, T_R, d_LR)
            J_diff_only = gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)
            # Should be close to pure diffusion (within 30%)
            @test abs(J_low - J_diff_only) / abs(J_diff_only) < 0.3
        end
    end

    @testset "Darcy-Fick Full Mesh Flux" begin
        # Create mesh with 2 components (solid + gas)
        mesh = create_uniform_mesh(0.01, 10, Val(2))
        n_cells = length(mesh.cells)
        n_faces = length(mesh.faces)

        # Initialize uniform state
        T_uniform = 400.0
        ξ_solid = 500.0
        ξ_gas = 0.5
        for i = 1:n_cells
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] =
                Pyrolysis.Geometry.CellState{2}(T_uniform, (ξ_solid, ξ_gas), V)
        end

        # Create property cache
        cache = PropertyCache(n_cells, n_faces)

        # Set up cache values
        fill!(cache.ρ, ξ_solid + ξ_gas)
        fill!(cache.λ_face, 0.01)  # kg/(m·s)

        @testset "Zero velocity matches diffusion-only" begin
            v_gas = zeros(n_faces)

            J_df = zeros(2, n_faces)
            J_diff = zeros(2, n_faces)

            compute_gas_fluxes_darcy_fick!(
                J_df,
                StateCache(mesh, two_comp),
                mesh,
                cache,
                two_comp,
                v_gas,
                fill(1.0, n_cells),
            )
            compute_gas_fluxes!(J_diff, StateCache(mesh, two_comp), mesh, cache, two_comp)

            # Should be identical
            for f = 1:n_faces
                for j = 1:2
                    @test J_df[j, f] ≈ J_diff[j, f] atol = 1e-15
                end
            end
        end

        @testset "Non-zero velocity adds advection" begin
            # Create concentration gradient
            for i = 1:n_cells
                ξ_gas_i = 0.5 + 0.1 * (i - 1)  # Increasing gas concentration
                V = mesh.cell_volumes[i]
                mesh.cell_states[i] =
                    Pyrolysis.Geometry.CellState{2}(T_uniform, (ξ_solid, ξ_gas_i), V)
            end

            # Positive velocity (flow in +z direction)
            v_gas = fill(0.1, n_faces)
            v_gas[1] = 0.0   # Boundary
            v_gas[end] = 0.0  # Boundary

            J_df = zeros(2, n_faces)
            J_diff = zeros(2, n_faces)

            compute_gas_fluxes_darcy_fick!(
                J_df,
                StateCache(mesh, two_comp),
                mesh,
                cache,
                two_comp,
                v_gas,
                fill(1.0, n_cells),
            )
            compute_gas_fluxes!(J_diff, StateCache(mesh, two_comp), mesh, cache, two_comp)

            # Darcy-Fick should differ from diffusion-only due to advection
            interior_differ = false
            for f = 2:(n_faces-1)
                if abs(J_df[2, f] - J_diff[2, f]) > 1e-10
                    interior_differ = true
                    break
                end
            end
            @test interior_differ

            # For gas component (index 2), advection should add positive flux
            # since v_g > 0 and ξ_gas increases with z
            for f = 2:(n_faces-1)
                face = mesh.faces[f]
                if face.cell_ids[1] != 0 && face.cell_ids[2] != 0
                    # Advective flux = v_g * ξ_upwind > 0
                    # So J_df should be larger (more positive) than J_diff
                    @test J_df[2, f] > J_diff[2, f]
                end
            end
        end

        @testset "Solid component unaffected" begin
            v_gas = fill(0.1, n_faces)
            v_gas[1] = 0.0
            v_gas[end] = 0.0

            J_df = zeros(2, n_faces)

            compute_gas_fluxes_darcy_fick!(
                J_df,
                StateCache(mesh, two_comp),
                mesh,
                cache,
                two_comp,
                v_gas,
                fill(1.0, n_cells),
            )

            # Solid component (index 1) should have zero flux
            for f = 1:n_faces
                @test J_df[1, f] == 0.0
            end
        end
    end

end # main testset
# ===========================================================================
# Gas-coupled energy source
# ===========================================================================

@testset "Gas Convective Source" begin

    @testset "Gas Convective Source (Volumetric Form)" begin
        # Setup two-component material with constant gas cp = 1000 J/(kg·K)
        solid = Pyrolysis.Materials.SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2)
        gas = Pyrolysis.Materials.GaseousComponent(
            :gas,
            M = 0.029,
            c = 1000.0,
            k = 0.025,
            λ = 1e-5,
        )
        two_comp = Material(name = :test, components = (solid, gas), reactions = ())

        # Create mesh
        n_cells = 10
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(2))
        n_faces = length(mesh.faces)
        NC = 2

        # Initialize state: linear temperature gradient, uniform gas
        ξ_solid = 900.0
        ξ_gas = 10.0
        dTdx_exact = 20000.0  # 200 K over 0.01 m
        for i = 1:n_cells
            x_i = mesh.cell_centers[i]
            T_i = 300.0 + dTdx_exact * x_i  # Linear T(x) = 300 + 20000*x
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(T_i, (ξ_solid, ξ_gas), V)
        end

        # The boundary cells' gradient blends in the boundary-face slope
        # −q/k at the solved surface temperature. q_zero ⟹ T_s = T_cell
        # (an adiabatic surface); the Linear-T testset instead supplies the
        # conductive flux consistent with the manufactured profile.
        prop_cache = Pyrolysis.Properties.PropertyCache(n_cells, n_faces)
        q_zero = zeros(n_faces)
        f_bot = mesh.cells[1].face_ids[1]
        f_top = mesh.cells[n_cells].face_ids[2]

        @testset "Zero gas flux gives zero source" begin
            J_gas = zeros(NC, n_faces)
            S_conv = zeros(n_cells)
            state = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state, mesh, two_comp)
            compute_gas_convective_source!(
                S_conv,
                J_gas,
                q_zero,
                prop_cache.k,
                state,
                mesh,
                two_comp,
            )

            for i = 1:n_cells
                @test abs(S_conv[i]) < 1e-14
            end
        end

        @testset "Uniform temperature gives zero source" begin
            # Set uniform temperature
            for i = 1:n_cells
                V = mesh.cell_volumes[i]
                mesh.cell_states[i] =
                    Pyrolysis.Geometry.CellState{2}(400.0, (ξ_solid, ξ_gas), V)
            end

            J_gas = zeros(NC, n_faces)
            J_gas[2, :] .= 0.5  # Nonzero gas flux
            S_conv = zeros(n_cells)
            state = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state, mesh, two_comp)
            compute_gas_convective_source!(
                S_conv,
                J_gas,
                q_zero,
                prop_cache.k,
                state,
                mesh,
                two_comp,
            )

            for i = 1:n_cells
                @test abs(S_conv[i]) < 1e-12
            end

            # Restore temperature gradient
            for i = 1:n_cells
                x_i = mesh.cell_centers[i]
                T_i = 300.0 + dTdx_exact * x_i
                V = mesh.cell_volumes[i]
                mesh.cell_states[i] =
                    Pyrolysis.Geometry.CellState{2}(T_i, (ξ_solid, ξ_gas), V)
            end
        end

        @testset "Linear T with constant J gives correct S_conv" begin
            # S_conv = -c_g * J * dT/dx
            # Linear T: dT/dx = 200/(0.01) = 20000 K/m
            # With J = 0.5 kg/(m²·s), c_g = 1000 J/(kg·K):
            # S_conv = -1000 * 0.5 * 20000 = -10,000,000 W/m³
            J_gas = zeros(NC, n_faces)
            J_gas[2, :] .= 0.5  # Constant gas flux at all faces
            S_conv = zeros(n_cells)
            state = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state, mesh, two_comp)
            # Boundary fluxes consistent with the linear profile: q = −k·dT/dz
            # at both boundary faces (the +z convention), so the blended
            # boundary gradient sees the exact surface-side slope.
            q_lin = zeros(n_faces)
            q_lin[f_bot] = -prop_cache.k[1] * dTdx_exact
            q_lin[f_top] = -prop_cache.k[n_cells] * dTdx_exact
            compute_gas_convective_source!(
                S_conv,
                J_gas,
                q_lin,
                prop_cache.k,
                state,
                mesh,
                two_comp,
            )

            dTdx_expected = 200.0 / L  # = 20000 K/m
            cp_gas = 1000.0
            J_val = 0.5
            S_expected = -cp_gas * J_val * dTdx_expected

            # Interior cells should match exactly (central difference is exact for linear T)
            for i = 2:(n_cells-1)
                @test S_conv[i] ≈ S_expected rtol = 1e-10
            end

            # Boundary cells blend the interior slope with the exact
            # boundary-face slope −q/k — also exact for linear T.
            @test S_conv[1] ≈ S_expected rtol = 1e-10
            @test S_conv[n_cells] ≈ S_expected rtol = 1e-10
        end

        @testset "Source proportional to temperature gradient" begin
            J_gas = zeros(NC, n_faces)
            J_gas[2, :] .= 0.5
            S_conv1 = zeros(n_cells)
            state1 = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state1, mesh, two_comp)
            compute_gas_convective_source!(
                S_conv1,
                J_gas,
                q_zero,
                prop_cache.k,
                state1,
                mesh,
                two_comp,
            )

            # Double the temperature gradient
            dTdx_double = 2.0 * dTdx_exact
            for i = 1:n_cells
                x_i = mesh.cell_centers[i]
                T_i = 300.0 + dTdx_double * x_i
                V = mesh.cell_volumes[i]
                mesh.cell_states[i] =
                    Pyrolysis.Geometry.CellState{2}(T_i, (ξ_solid, ξ_gas), V)
            end

            S_conv2 = zeros(n_cells)
            state2 = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state2, mesh, two_comp)
            compute_gas_convective_source!(
                S_conv2,
                J_gas,
                q_zero,
                prop_cache.k,
                state2,
                mesh,
                two_comp,
            )

            # With doubled gradient, source should double
            for i = 1:n_cells
                if abs(S_conv1[i]) > 1e-10
                    @test S_conv2[i] / S_conv1[i] ≈ 2.0 rtol = 0.01
                end
            end

            # Restore original gradient
            for i = 1:n_cells
                x_i = mesh.cell_centers[i]
                T_i = 300.0 + dTdx_exact * x_i
                V = mesh.cell_volumes[i]
                mesh.cell_states[i] =
                    Pyrolysis.Geometry.CellState{2}(T_i, (ξ_solid, ξ_gas), V)
            end
        end

        @testset "Only gas components contribute" begin
            # Solid flux should not contribute to S_conv
            J_gas_solid_only = zeros(NC, n_faces)
            J_gas_solid_only[1, :] .= 0.5  # Flux in solid component
            S_conv = zeros(n_cells)
            state = StateCache(mesh, two_comp)
            Pyrolysis.Properties.update_properties!(prop_cache, state, mesh, two_comp)
            compute_gas_convective_source!(
                S_conv,
                J_gas_solid_only,
                q_zero,
                prop_cache.k,
                state,
                mesh,
                two_comp,
            )

            for i = 1:n_cells
                @test abs(S_conv[i]) < 1e-14
            end
        end
    end


end  # ALE Integration tests

# ============================================================================
# Pressure Boundary Conditions
# ============================================================================

@testset "Pressure Boundary Conditions" begin

    # Test materials. The Darcy-velocity assertions in
    # `Integration with update_darcy_properties!` need the solid component to
    # actually carry permeability — `κ = 0` would zero the effective
    # permeability and produce zero gas velocity regardless of pressure
    # gradient.
    two_comp = Material(
        name = :two_comp,
        components = (
            SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-11),
            GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )

    @testset "DirichletPressureBC Construction" begin
        # Default constructor
        bc_default = DirichletPressureBC()
        @test bc_default.P == 101325.0  # Standard atmospheric

        # Custom pressure
        bc_custom = DirichletPressureBC(50000.0)
        @test bc_custom.P == 50000.0
    end

    @testset "NeumannPressureBC Construction" begin
        # Default constructor (impermeable)
        bc_default = NeumannPressureBC()
        @test bc_default.∂P∂n == 0.0

        # Prescribed gradient
        bc_gradient = NeumannPressureBC(1000.0)
        @test bc_gradient.∂P∂n == 1000.0
    end

    @testset "ConvectivePressureBC Construction" begin
        bc = ConvectivePressureBC(1e-6, 101325.0)
        @test bc.h_P == 1e-6
        @test bc.P_ext == 101325.0
    end

    @testset "DirichletPressureBC Application" begin
        # Create mesh with gas at elevated pressure
        n_cells = 5
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(2))
        n_faces = length(mesh.faces)

        # Initialize with high internal pressure (gas at high T and concentration)
        T_high = 600.0
        ξ_solid = 800.0
        ξ_gas = 5.0  # kg/m³
        for i = 1:n_cells
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] =
                Pyrolysis.Geometry.CellState{2}(T_high, (ξ_solid, ξ_gas), V)
        end
        Pyrolysis.Geometry.update_geometry!(mesh)

        # Compute pressure and permeability fields
        P = zeros(n_cells)
        κ = fill(1e-12, n_cells)  # Constant permeability
        # Per-cell pore-gas viscosity, as `update_darcy_properties!` fills it.
        μ_cells = [
            mixture_gas_viscosity(mesh.cell_states[i].ξ, mesh.cell_states[i].T, two_comp) for i = 1:n_cells
        ]

        for i = 1:n_cells
            state = mesh.cell_states[i]
            φ_i = porosity(state.ξ, state.T, two_comp)
            P[i] = Physics.ideal_gas_pressure(state.ξ, state.T, φ_i, two_comp)
        end

        # P_internal should be higher than atmospheric
        P_atm = 101325.0
        @test all(P .> P_atm)

        # Initialize v_gas to zeros (boundary faces start at zero)
        v_gas = zeros(n_faces)

        # Apply DirichletPressureBC at top boundary (boundary_id=1, face at high z)
        # Build type-stable BoundaryConditionSet and pass boundary face lists
        bc_set = BoundaryConditionSet(
            mass = (top = DirichletPressureBC(P_atm), bottom = ImpermeableBC()),
        )
        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]
        apply_pressure_bc_to_darcy_velocity!(
            v_gas,
            StateCache(mesh, two_comp),
            mesh,
            P,
            κ,
            μ_cells,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        # Top boundary face should have positive velocity (outflow, gas escaping)
        # The top boundary face is the last face
        top_face_idx = n_faces  # Face at z = L
        @test v_gas[top_face_idx] > 0  # Outflow

        # Bottom boundary should remain zero (ImpermeableBC)
        bottom_face_idx = 1  # Face at z = 0
        @test v_gas[bottom_face_idx] == 0.0
    end

    @testset "NeumannPressureBC Application - Zero Gradient" begin
        n_cells = 5
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(2))
        n_faces = length(mesh.faces)

        # Initialize mesh
        for i = 1:n_cells
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(500.0, (900.0, 1.0), V)
        end
        Pyrolysis.Geometry.update_geometry!(mesh)

        P = fill(150000.0, n_cells)  # Elevated pressure
        κ = fill(1e-12, n_cells)
        μ_cells = [
            mixture_gas_viscosity(mesh.cell_states[i].ξ, mesh.cell_states[i].T, two_comp) for i = 1:n_cells
        ]
        v_gas = zeros(n_faces)

        # NeumannPressureBC with zero gradient = impermeable
        bc_set = BoundaryConditionSet(
            mass = (top = NeumannPressureBC(0.0), bottom = NeumannPressureBC(0.0)),
        )
        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]
        apply_pressure_bc_to_darcy_velocity!(
            v_gas,
            StateCache(mesh, two_comp),
            mesh,
            P,
            κ,
            μ_cells,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        # All boundary faces should have zero velocity
        @test v_gas[1] == 0.0
        @test v_gas[n_faces] == 0.0
    end

    @testset "ConvectivePressureBC Application" begin
        n_cells = 5
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(2))
        n_faces = length(mesh.faces)

        # Initialize with internal pressure higher than external
        T = 500.0
        for i = 1:n_cells
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(T, (800.0, 2.0), V)
        end
        Pyrolysis.Geometry.update_geometry!(mesh)

        P = zeros(n_cells)
        for i = 1:n_cells
            state = mesh.cell_states[i]
            φ_i = porosity(state.ξ, state.T, two_comp)
            P[i] = Physics.ideal_gas_pressure(state.ξ, state.T, φ_i, two_comp)
        end

        κ = fill(1e-12, n_cells)
        μ_cells = [
            mixture_gas_viscosity(mesh.cell_states[i].ξ, mesh.cell_states[i].T, two_comp) for i = 1:n_cells
        ]
        v_gas = zeros(n_faces)

        P_ext = 101325.0
        h_P = 1e-6  # Pressure transfer coefficient

        mass_bcs_top = ConvectivePressureBC(h_P, P_ext)

        bc_set = BoundaryConditionSet(mass = (top = mass_bcs_top, bottom = ImpermeableBC()))
        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]
        apply_pressure_bc_to_darcy_velocity!(
            v_gas,
            StateCache(mesh, two_comp),
            mesh,
            P,
            κ,
            μ_cells,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        # Top boundary: v = h_P × (P_int - P_ext)
        top_cell = n_cells
        expected_v = h_P * (P[top_cell] - P_ext)

        # Should be positive (outflow) since P_int > P_ext
        @test v_gas[n_faces] > 0
        @test v_gas[n_faces] ≈ expected_v rtol = 0.1
    end

    @testset "Integration with update_darcy_properties!" begin
        # Test that update_darcy_properties! correctly applies pressure BCs
        n_cells = 5
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(2))
        n_faces = length(mesh.faces)

        # Initialize with elevated gas concentration
        T = 600.0
        ξ_solid = 700.0
        ξ_gas = 3.0
        for i = 1:n_cells
            V = mesh.cell_volumes[i]
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(T, (ξ_solid, ξ_gas), V)
        end
        Pyrolysis.Geometry.update_geometry!(mesh)

        # Create property cache
        cache = PropertyCache(mesh)
        update_properties!(cache, StateCache(mesh, two_comp), mesh, two_comp)

        # Debug: Check that update_properties! worked
        if cache.ρ[1] == 0.0
            @warn "update_properties! failed: cache.ρ[1] = 0.0, expected ~1000"
            @warn "  T = $(mesh.cell_states[1].T), ξ = $(mesh.cell_states[1].ξ)"
            # Try direct computation
            ρ_direct = Properties.mixture_density(
                mesh.cell_states[1].ξ,
                mesh.cell_states[1].T,
                two_comp,
            )
            @warn "  Direct ρ computation: $ρ_direct"
        end

        # Update Darcy properties with BCs (component-based κ from two_comp)
        bc_set = BoundaryConditionSet(
            mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
        )
        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]
        update_darcy_properties!(
            cache,
            StateCache(mesh, two_comp),
            mesh,
            two_comp,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )

        # Debug: Check that update_darcy_properties! worked
        # (Darcy-only fields live under cache.darcy.*)
        if cache.darcy.P[1] == 0.0
            @warn "update_darcy_properties! failed: cache.darcy.P[1] = 0.0, expected ~340000"
            @warn "  cache.ρ = $(cache.ρ)"
            @warn "  cache.darcy.φ = $(cache.darcy.φ)"
            @warn "  cache.darcy.κ = $(cache.darcy.κ)"
        end

        # Internal pressure should be higher than atmospheric (gas at high T)
        @test all(cache.darcy.P .> 101325.0)

        # Top boundary should have outflow
        @test cache.darcy.v_gas[n_faces] > 0

        # Bottom boundary should be zero
        @test cache.darcy.v_gas[1] == 0.0
    end

end  # Pressure Boundary Conditions tests


# ==============================================================================
# Species-Specific Transport Tests
# ==============================================================================
@testset "Bruggeman Effective Medium Theory" begin

    @testset "Single component (trivial case)" begin
        # Single component should return that component's property
        v = (1.0,)
        k = (5.0,)
        k_eff = bruggeman_effective_property(v, k)
        @test k_eff ≈ 5.0 rtol = 1e-10
    end

    @testset "Two equal components" begin
        # 50/50 mix of equal properties should return that property
        v = (0.5, 0.5)
        k = (10.0, 10.0)
        k_eff = bruggeman_effective_property(v, k)
        @test k_eff ≈ 10.0 rtol = 1e-6
    end

    @testset "Bounds: Bruggeman between series and parallel" begin
        # Bruggeman should always lie between series and parallel bounds
        v = (0.7, 0.3)
        k = (1.0, 10.0)

        k_parallel = v[1] * k[1] + v[2] * k[2]  # = 3.7
        k_series = 1.0 / (v[1] / k[1] + v[2] / k[2])  # ≈ 1.43

        k_eff = bruggeman_effective_property(v, k)

        @test k_eff > k_series
        @test k_eff < k_parallel
        # Exact root of 2k² − 0.1k − 10 = 0
        @test k_eff ≈ 2.2612078 atol = 1e-4
    end

    @testset "Symmetry in component ordering" begin
        # Bruggeman is symmetric under relabeling of components
        # Swapping (v1,k1) with (v2,k2) should give same result
        v1 = (0.3, 0.7)
        k1 = (1.0, 10.0)
        v2 = (0.7, 0.3)
        k2 = (10.0, 1.0)

        k_eff1 = bruggeman_effective_property(v1, k1)
        k_eff2 = bruggeman_effective_property(v2, k2)

        # These SHOULD be equal (symmetric in component labeling)
        @test k_eff1 ≈ k_eff2 rtol = 1e-6

        # But NOT symmetric when only volume fractions swap (properties stay same)
        v3 = (0.7, 0.3)
        k3 = (1.0, 10.0)  # Same as k1
        k_eff3 = bruggeman_effective_property(v3, k3)

        # k_eff3 should differ from k_eff1 (different volume fractions for same properties)
        @test abs(k_eff3 - k_eff1) > 0.1
    end

    @testset "Three components" begin
        v = (0.5, 0.3, 0.2)
        k = (1.0, 5.0, 10.0)

        k_parallel = sum(v .* k)  # = 4.5
        # Series is complex for 3 components
        k_series = 1.0 / sum(v ./ k)  # ≈ 1.67

        k_eff = bruggeman_effective_property(v, k)

        @test k_eff > k_series
        @test k_eff < k_parallel
    end

    @testset "Zero volume fraction handled" begin
        # Component with zero volume should not contribute
        v = (0.0, 1.0)
        k = (100.0, 5.0)

        k_eff = bruggeman_effective_property(v, k)
        @test k_eff ≈ 5.0 rtol = 1e-6
    end

    @testset "High contrast ratio" begin
        # Char/virgin-like contrast: 10:1 ratio
        v = (0.6, 0.4)  # 60% high-k, 40% low-k
        k = (10.0, 1.0)

        k_eff = bruggeman_effective_property(v, k)

        # Should be dominated by high-k phase but reduced by low-k inclusions
        @test k_eff > 3.0  # Well above series
        @test k_eff < 6.0  # Well below parallel
    end

    @testset "Integration with effective_conductivity" begin
        # Test that BRUGGEMAN mixing rule works in effective_conductivity
        brugg_mat = Material(
            name = :brugg_test,
            components = (
                SolidComponent(:comp1, ρ = 1000.0, c = 1500.0, k = 1.0, κ = 0.0),
                SolidComponent(:comp2, ρ = 500.0, c = 1000.0, k = 10.0, κ = 0.0),
            ),
            reactions = (),
            mixing = MaterialMixing(
                k = MixingSpec(BRUGGEMAN, 0.5),  # β unused for Bruggeman
                λ = MixingSpec(PARALLEL, 0.5),
            ),
        )

        T = 400.0

        # 50/50 by mass → volume fractions depend on densities
        ξ1 = 500.0  # kg/m³
        ξ2 = 500.0  # kg/m³
        # v1 = ξ1/ρ1 / (ξ1/ρ1 + ξ2/ρ2) = 0.5 / (0.5 + 1.0) = 0.333
        # v2 = 0.667

        ξ = (ξ1, ξ2)
        k_eff = effective_conductivity(ξ, T, brugg_mat)

        # Parallel: 0.333 × 1.0 + 0.667 × 10.0 = 7.0
        # Series: 1.0 / (0.333/1.0 + 0.667/10.0) = 2.5
        # Bruggeman should be between
        @test k_eff > 2.5
        @test k_eff < 7.0
    end

    @testset "BRUGGEMAN rejected for gas transfer" begin
        # BRUGGEMAN is not a physically valid closure for gas-transfer
        # coefficients (or for permeability). It's now rejected at evaluation time.
        brugg_gas_mat = Material(
            name = :brugg_gas_test,
            components = (
                SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 0.0, λ = 1e-6),
                GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
            ),
            reactions = (),
            mixing = MaterialMixing(λ = MixingSpec(BRUGGEMAN)),
        )

        ξ = (800.0, 0.5)
        @test_throws ArgumentError effective_gas_transfer(ξ, 400.0, brugg_gas_mat)
    end

    @testset "Convergence verification" begin
        # Test that Bruggeman satisfies its defining equation
        v = (0.4, 0.6)
        k = (2.0, 8.0)

        k_eff = bruggeman_effective_property(v, k)

        # Verify: Σ vⱼ(kⱼ - k_eff)/(kⱼ + 2k_eff) ≈ 0
        residual = sum(v[i] * (k[i] - k_eff) / (k[i] + 2 * k_eff) for i = 1:2)
        @test abs(residual) < 1e-8
    end

    @testset "Percolation threshold with insulating phase" begin
        # An insulating phase contributes the exact k → 0 limit −v/2 to the
        # self-consistency equation. For conducting fraction p in insulator,
        # the analytical root is k_eff = k_c·(3p − 1)/2 for p ≥ 1/3, else 0.

        # Confirming case: v = (0.5, 0.5), k = (1, 0) → 0.25.
        # The pre-fix code skipped the insulating phase and returned 1.0.
        @test bruggeman_effective_property((0.5, 0.5), (1.0, 0.0)) ≈ 0.25 atol = 1e-8

        # General above-threshold point: p = 0.6, k_c = 2 → 2(0.8)/2 = 0.8
        @test bruggeman_effective_property((0.6, 0.4), (2.0, 0.0)) ≈ 0.8 atol = 1e-8

        # Below the v = 1/3 percolation threshold the conducting phase is
        # disconnected: k_eff = 0
        @test bruggeman_effective_property((0.3, 0.7), (1.0, 0.0)) == 0.0
        @test bruggeman_effective_property((0.1, 0.9), (5.0, 0.0)) == 0.0

        # At the threshold itself the root is 0
        @test bruggeman_effective_property((1 / 3, 2 / 3), (1.0, 0.0)) ≈ 0.0 atol = 1e-6

        # Three phases, one insulating: the returned root must satisfy the
        # full residual including the −v/2 insulator term
        v = (0.4, 0.3, 0.3)
        k = (1.0, 5.0, 0.0)
        k_eff = bruggeman_effective_property(v, k)
        @test k_eff > 0.0
        residual =
            v[1] * (k[1] - k_eff) / (k[1] + 2k_eff) +
            v[2] * (k[2] - k_eff) / (k[2] + 2k_eff) - 0.5 * v[3]
        @test abs(residual) < 1e-8
    end

end  # Bruggeman tests

@testset "Carman–Kozeny Permeability" begin
    # κ = ℓ²·φ³/(180·(1-φ)²): porosity-based permeability closure on the κ channel.

    ℓ = 50e-6  # 50 µm characteristic particle/pore diameter

    ck_mat = Material(
        name = :ck_test,
        components = (
            SolidComponent(
                :virgin,
                ρ = 400.0,
                c = 1500.0,
                k = 0.2,
                κ = 1e-14,
                ρ_intrinsic = 1500.0,
            ),
            GaseousComponent(:gas, M = 0.030, c = 1100.0, k = 0.05, λ = 1e-5),
        ),
        reactions = (),
        mixing = MaterialMixing(κ = MixingSpec(CARMAN_KOZENY; length_scale = ℓ)),
    )

    @testset "Construction validation" begin
        # CARMAN_KOZENY requires a positive, finite length_scale.
        @test_throws ArgumentError MixingSpec(CARMAN_KOZENY)
        @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = 0.0)
        @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = -1.0)
        @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = NaN)
        spec = MixingSpec(CARMAN_KOZENY; length_scale = ℓ)
        @test spec.rule == CARMAN_KOZENY
        @test spec.ℓ == ℓ
        # Other rules ignore length_scale and keep working as before.
        @test MixingSpec(WEIGHTED, 0.5).rule == WEIGHTED
    end

    @testset "Formula agreement across porosities" begin
        T = 400.0
        for ρs in (1400.0, 1000.0, 500.0, 100.0)
            ξ = (ρs, 0.05)
            φ = porosity(ξ, T, ck_mat)
            κ = effective_permeability(ξ, T, ck_mat)
            # Free-flow cap at ℓ²: the packed-bed closure
            # diverges as φ → 1; κ saturates at the physical free-flow scale.
            κ_manual = min(ℓ^2 * φ^3 / (180 * (1 - φ)^2), ℓ^2)
            @test isapprox(κ, κ_manual; rtol = 1e-12)
            @test κ ≥ 0.0
        end
    end

    @testset "Monotonic in porosity" begin
        # Lower solid loading → higher porosity → higher permeability.
        T = 400.0
        κ_dense = effective_permeability((1400.0, 0.05), T, ck_mat)
        κ_open = effective_permeability((300.0, 0.05), T, ck_mat)
        @test κ_open > κ_dense
    end

    # Probe AD the way the solver does: all components are perturbed into the
    # same Dual type (homogeneous NTuple), so the gas entry carries `one(s)`.
    _κ_of_solid(s) = effective_permeability((s, 0.05 * one(s)), 400.0, ck_mat)

    @testset "Free-flow cap at high porosity" begin
        T = 400.0
        # φ = 1 − 100/1500 ≈ 0.933: CK ratio φ³/(180(1−φ)²) ≈ 1.02 > 1 → capped.
        κ_high = effective_permeability((100.0, 0.05), T, ck_mat)
        @test κ_high == ℓ^2
        # Pure-gas cell (φ = 1 exactly): the old floor produced κ = ℓ²·1e30;
        # the cap keeps it at ℓ², and the AD derivative is finite (zero on the
        # capped branch — consistent with the residual's min()).
        κ_pure = effective_permeability((0.0, 0.05), T, ck_mat)
        @test κ_pure == ℓ^2
        g = ForwardDiff.derivative(_κ_of_solid, 0.0)
        @test isfinite(g)
        @test g == 0.0
    end

    @testset "AD-clean (ForwardDiff through effective_permeability)" begin
        g = ForwardDiff.derivative(_κ_of_solid, 800.0)
        @test isfinite(g)
        @test g < 0.0  # more solid → less porosity → less permeability
    end

    @testset "Cache-aware _permeability_from_v rejects CK (needs ξ)" begin
        # CK depends on porosity(ξ); the v-only entry point must say so clearly.
        v = (0.9, 0.1)
        @test_throws ArgumentError Pyrolysis.Properties._permeability_from_v(
            v,
            400.0,
            ck_mat,
        )
    end
end  # Carman–Kozeny tests
