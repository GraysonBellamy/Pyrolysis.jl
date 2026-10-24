"""
Pointwise gas closures of the shared layer.

- The ideal-gas pressure P = Σⱼ ξⱼRT/(Mⱼφ) over the gas components, at the
  pore fraction φ = 1 − Σ skeletal volumes that the state mapping floors
  smoothly.
- The permeability mixing rules on the renormalised solid basis, and the
  Carman–Kozeny closure of the pore fraction with its free-flow cap.
- The Chapman–Enskog binary diffusivity and its Neufeld collision integral.
- The volume-fraction-gradient gas flux J = −λ(ξ_R T_R − ξ_L T_L)/(T̄ d).
- The Bruggeman effective medium.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Properties:
    skeletal_volume_fraction,
    effective_permeability,
    effective_conductivity,
    effective_gas_transfer,
    bruggeman_effective_property,
    collision_integral_diffusion,
    chapman_enskog_diffusivity,
    CHAPMAN_ENSKOG_PREFACTOR,
    SIGMA_AIR,
    EPSILON_K_AIR,
    M_AIR,
    T_GAS_PROPERTY_FLOOR
using Pyrolysis.Physics: ideal_gas_pressure, gas_flux

# The permeability at the state's pore fraction 1 − Σ skeletal volumes, the
# porosity the state mapping forms; the per-component mixes do not read it.
_κ(ξ, T, material) =
    effective_permeability(ξ, T, material, 1 - skeletal_volume_fraction(ξ, T, material))

@testset "Gas closures" begin
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
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )

    # The pore fraction the state mapping floors smoothly.
    pore_fraction(ξ, T, mat) = 1 - skeletal_volume_fraction(ξ, T, mat)

    # ========================================================================
    @testset "ideal-gas pressure: P = Σ ξ R T/(M φ) over the gases" begin
        # Pure gas: φ = 1, so P = ξ_gas R T / M.
        ξ = (0.0, 1.0)
        T = 400.0
        φ = pore_fraction(ξ, T, two_comp)
        @test φ ≈ 1.0 rtol = 1e-10
        @test ideal_gas_pressure(ξ, T, φ, two_comp) ≈ 1.0 * R_GAS * T / 0.044 rtol = 1e-10

        # No gas component, or no gas: exactly zero.
        ξs = (1000.0,)
        @test ideal_gas_pressure(ξs, T, pore_fraction(ξs, T, solid_only), solid_only) == 0.0
        ξ0 = (1000.0, 0.0)
        @test ideal_gas_pressure(ξ0, T, pore_fraction(ξ0, T, two_comp), two_comp) == 0.0

        # At fixed ξ the pressure is linear in T (φ depends on ξ only here).
        P1 = ideal_gas_pressure(ξ, 300.0, pore_fraction(ξ, 300.0, two_comp), two_comp)
        P2 = ideal_gas_pressure(ξ, 600.0, pore_fraction(ξ, 600.0, two_comp), two_comp)
        @test P2 / P1 ≈ 2.0 rtol = 1e-10
    end

    @testset "the pore fraction sums the skeletal volumes and the pressure divides by it" begin
        # Two solids: 500/1000 + 50/200 = 0.75 of the cell, so φ = 0.25.
        ξ = (500.0, 50.0, 0.5)
        T = 500.0
        φ = pore_fraction(ξ, T, three_comp)
        @test φ ≈ 0.25 rtol = 1e-10
        @test ideal_gas_pressure(ξ, T, φ, three_comp) ≈ (0.5 * R_GAS * T / 0.028) / φ rtol =
            1e-10
    end

    @testset "pressure: ForwardDiff through the pore fraction equals central differences" begin
        ξ = (500.0, 1.0)
        T = 400.0
        h = 1e-6
        P(ξl, Tl) = ideal_gas_pressure(ξl, Tl, pore_fraction(ξl, Tl, two_comp), two_comp)
        g = ForwardDiff.gradient(x -> P((x[1], x[2]), x[3]), [ξ[1], ξ[2], T])
        @test g[3] ≈ (P(ξ, T + h) - P(ξ, T - h)) / (2h) rtol = 1e-6
        @test g[1] ≈ (P((ξ[1] + h, ξ[2]), T) - P((ξ[1] - h, ξ[2]), T)) / (2h) rtol = 1e-6
        @test g[2] ≈ (P((ξ[1], ξ[2] + h), T) - P((ξ[1], ξ[2] - h), T)) / (2h) rtol = 1e-6
    end

    # ========================================================================
    @testset "permeability: one solid, two solids in parallel, gases add nothing" begin
        single = Material(
            name = :single_perm,
            components = (
                SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-12),
            ),
            reactions = (),
        )
        @test _κ((1000.0,), 500.0, single) ≈ 1e-12 rtol = 0.01

        mix = Material(
            name = :mix_perm,
            components = (
                SolidComponent(:virgin, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-15),
                SolidComponent(:char, ρ = 200.0, c = 800.0, k = 0.1, κ = 1e-12),
                GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
            ),
            reactions = (),
            mixing = MaterialMixing(κ = MixingSpec(PARALLEL)),
        )
        κ_virgin = _κ((1000.0, 0.0, 0.0), 500.0, mix)
        κ_char = _κ((0.0, 200.0, 0.0), 500.0, mix)
        @test κ_virgin ≈ 1e-15 rtol = 0.1
        @test κ_char ≈ 1e-12 rtol = 0.1
        # Half the solid volume each: between the pure values.
        κ_mix = _κ((500.0, 100.0, 0.0), 500.0, mix)
        @test κ_virgin < κ_mix < κ_char

        # A gas component carries no permeability of its own.
        with_gas = Material(
            name = :gas_test,
            components = (
                SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, κ = 1e-12),
                GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
            ),
            reactions = (),
        )
        @test _κ((1000.0, 0.0), 500.0, with_gas) == 1e-12
        @test _κ((1000.0, 10.0), 500.0, with_gas) > 0
    end

    @testset "permeability mixes on the solid basis; an impermeable phase zeroes the series bound" begin
        # Volume fractions 0.5 and 0.2 of solid, 0.3 of gas.
        T = 400.0
        ξ = (500.0, 200.0, 0.4)
        uniform(rule) = Material(
            name = :c3_uniform,
            components = (
                SolidComponent(:s1, ρ = 1000.0, c = 1500.0, k = 0.2, κ = 1e-13),
                SolidComponent(:s2, ρ = 1000.0, c = 1000.0, k = 0.1, κ = 1e-13),
                GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
            ),
            reactions = (),
            mixing = MaterialMixing(κ = MixingSpec(rule)),
        )
        # A uniform solid κ₀ is κ₀ under both bounds, the gas notwithstanding.
        @test _κ(ξ, T, uniform(PARALLEL)) ≈ 1e-13 rtol = 1e-12
        @test _κ(ξ, T, uniform(SERIES)) ≈ 1e-13 rtol = 1e-12

        impermeable(rule) = Material(
            name = :c3_impermeable,
            components = (
                SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2),
                SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.1, κ = 1e-12),
                GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
            ),
            reactions = (),
            mixing = MaterialMixing(κ = MixingSpec(rule)),
        )
        @test _κ(ξ, T, impermeable(SERIES)) == 0.0
        @test _κ(ξ, T, impermeable(PARALLEL)) ≈ (0.2 / 0.7) * 1e-12 rtol = 1e-12
        @test _κ(ξ, T, impermeable(WEIGHTED)) ≈ 0.5 * (0.2 / 0.7) * 1e-12 rtol = 1e-12

        # ForwardDiff through the WEIGHTED rule equals central differences.
        mat = impermeable(WEIGHTED)
        g = ForwardDiff.gradient(x -> _κ(Tuple(x), T, mat), [ξ...])
        h = 1e-3
        for i = 1:3
            ξp = ntuple(j -> j == i ? ξ[j] + h : ξ[j], 3)
            ξm = ntuple(j -> j == i ? ξ[j] - h : ξ[j], 3)
            fd = (_κ(ξp, T, mat) - _κ(ξm, T, mat)) / (2h)
            @test isapprox(g[i], fd; rtol = 1e-5, atol = 1e-30)
        end
    end

    # ========================================================================
    @testset "Carman–Kozeny: κ = min(ℓ²φ³/(180(1 − φ)²), ℓ²)" begin
        ℓ = 50e-6
        ck = Material(
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

        @testset "the closure needs a positive finite length scale" begin
            @test_throws ArgumentError MixingSpec(CARMAN_KOZENY)
            @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = 0.0)
            @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = -1.0)
            @test_throws ArgumentError MixingSpec(CARMAN_KOZENY; length_scale = NaN)
            spec = MixingSpec(CARMAN_KOZENY; length_scale = ℓ)
            @test spec.rule == CARMAN_KOZENY
            @test spec.ℓ == ℓ
            # The other rules ignore the length scale.
            @test MixingSpec(WEIGHTED, 0.5).rule == WEIGHTED
        end

        # The permeability at the state's pore fraction, as the state mapping
        # evaluates it.
        κ_ck(ξ, T) = effective_permeability(ξ, T, ck, pore_fraction(ξ, T, ck))

        @testset "closed form across porosities, increasing with porosity" begin
            T = 400.0
            for ρs in (1400.0, 1000.0, 500.0, 100.0)
                φ = 1 - ρs / 1500.0
                κ = κ_ck((ρs, 0.05), T)
                @test isapprox(κ, min(ℓ^2 * φ^3 / (180 * (1 - φ)^2), ℓ^2); rtol = 1e-12)
                @test κ ≥ 0.0
            end
            @test κ_ck((300.0, 0.05), T) > κ_ck((1400.0, 0.05), T)
        end

        # The solver's AD perturbs every component into the same Dual type.
        κ_of_solid(s) = κ_ck((s, 0.05 * one(s)), 400.0)

        @testset "the free-flow cap binds at high porosity with zero slope" begin
            # φ ≈ 0.933: the packed-bed ratio exceeds 1, so κ = ℓ².
            @test κ_ck((100.0, 0.05), 400.0) == ℓ^2
            # A pure-gas cell (φ = 1) stays at the cap with a finite, zero slope.
            @test κ_ck((0.0, 0.05), 400.0) == ℓ^2
            g = ForwardDiff.derivative(κ_of_solid, 0.0)
            @test isfinite(g)
            @test g == 0.0
        end

        @testset "below the cap more solid lowers κ, with a finite AD slope" begin
            g = ForwardDiff.derivative(κ_of_solid, 800.0)
            @test isfinite(g)
            @test g < 0.0
        end

        @testset "the volume-fraction entry point refuses the closure" begin
            @test_throws ArgumentError Pyrolysis.Properties._permeability_from_v(
                (0.9, 0.1),
                400.0,
                ck,
            )
        end
    end

    # ========================================================================
    @testset "Chapman–Enskog diffusivity" begin
        @testset "Neufeld collision integral at T* = 1, 2, 5, falling" begin
            Ω_1 = collision_integral_diffusion(1.0)
            Ω_2 = collision_integral_diffusion(2.0)
            Ω_5 = collision_integral_diffusion(5.0)
            @test Ω_1 ≈ 1.44 rtol = 0.02
            @test Ω_2 ≈ 1.075 rtol = 0.02
            @test Ω_5 ≈ 0.84 rtol = 0.03
            @test Ω_5 < Ω_2 < Ω_1
        end

        @testset "rises with T between T^1.5 and T^1.75" begin
            M = 0.044
            D_300 = chapman_enskog_diffusivity(300.0, M)
            D_500 = chapman_enskog_diffusivity(500.0, M)
            D_800 = chapman_enskog_diffusivity(800.0, M)
            @test D_300 < D_500 < D_800
            ratio = D_500 / D_300
            @test ratio > (500.0 / 300.0)^1.5 * 0.8
            @test ratio < (500.0 / 300.0)^1.75 * 1.2
        end

        @testset "inverse in pressure, faster for lighter gases" begin
            D_1atm = chapman_enskog_diffusivity(500.0, 0.028; P = 101325.0)
            D_2atm = chapman_enskog_diffusivity(500.0, 0.028; P = 202650.0)
            @test D_1atm / D_2atm ≈ 2.0 rtol = 0.01
            D = [chapman_enskog_diffusivity(500.0, M) for M in (0.018, 0.028, 0.044, 0.100)]
            @test issorted(D; rev = true)
            @test allunique(D)
        end

        @testset "literature magnitudes in air at 1 atm" begin
            # H₂O at 373 K ≈ 2.2e-5, CO₂ at 298 K ≈ 1.6e-5, N₂ at 300 K ≈ 2.0e-5 m²/s.
            @test 1e-5 < chapman_enskog_diffusivity(373.0, 0.018) < 5e-5
            @test 5e-6 < chapman_enskog_diffusivity(298.0, 0.044) < 5e-5
            @test 1e-5 < chapman_enskog_diffusivity(300.0, 0.028) < 5e-5
        end

        @testset "custom Lennard-Jones parameters and the air constants" begin
            # Methane's σ and ε/k.
            @test chapman_enskog_diffusivity(500.0, 0.016) > 0.0
            D_ch4 = chapman_enskog_diffusivity(500.0, 0.016; σ = 3.758, ε_k = 148.6)
            @test 1e-5 < D_ch4 < 1e-4
            @test CHAPMAN_ENSKOG_PREFACTOR > 0.0
            @test SIGMA_AIR ≈ 3.617 rtol = 0.01
            @test EPSILON_K_AIR ≈ 97.0 rtol = 0.01
            @test M_AIR ≈ 0.029 rtol = 0.01
        end

        @testset "a nonphysical trial temperature evaluates at the floor" begin
            D_neg = chapman_enskog_diffusivity(-43.0, 0.06)
            @test isfinite(D_neg) && D_neg > 0.0
            @test D_neg == chapman_enskog_diffusivity(0.0, 0.06)
            @test D_neg == chapman_enskog_diffusivity(T_GAS_PROPERTY_FLOOR, 0.06)
        end
    end

    # ========================================================================
    @testset "gas flux: J = −λ(ξ_R T_R − ξ_L T_L)/(T̄ d)" begin
        @testset "equal concentrations at unequal temperatures flow hot to cold" begin
            J = gas_flux(2e-5, 10.0, 10.0, 400.0, 800.0, 0.001)
            @test J < 0.0
            @test J ≈ -2e-5 * (10.0 * 800.0 - 10.0 * 400.0) / (600.0 * 0.001) atol = 1e-12
            @test J ≈ -2e-5 * 6666.6667 rtol = 1e-4
        end

        @testset "at a uniform temperature it is Fick's law" begin
            J = gas_flux(2e-5, 5.0, 15.0, 500.0, 500.0, 0.001)
            @test J ≈ -2e-5 * (15.0 - 5.0) / 0.001 atol = 1e-12
            @test J ≈ -0.2 atol = 1e-12
        end

        @testset "the product form equals −λ[∂ξ + (ξ̄/T̄)∂T]" begin
            ξL, ξR, TL, TR, λ, d = 8.0, 12.0, 500.0, 700.0, 3e-5, 0.002
            T̄, ξ̄ = (TL + TR) / 2, (ξL + ξR) / 2
            J_expanded = -λ * ((ξR - ξL) / d + (ξ̄ / T̄) * (TR - TL) / d)
            @test gas_flux(λ, ξL, ξR, TL, TR, d) ≈ J_expanded atol = 1e-14
        end

        @testset "it is the partial-pressure form −(λM/(R T̄))∂P" begin
            ξL, ξR, TL, TR, λ, d, M = 10.0, 10.0, 400.0, 800.0, 2e-5, 0.001, 0.028
            ∂P = (ξR * R_GAS * TR / M - ξL * R_GAS * TL / M) / d
            @test gas_flux(λ, ξL, ξR, TL, TR, d) ≈ -(λ * M) / (R_GAS * (TL + TR) / 2) * ∂P atol =
                1e-12
        end
    end

    # ========================================================================
    @testset "Bruggeman effective medium" begin
        @testset "a single phase, two equal phases and an empty phase return the property" begin
            @test bruggeman_effective_property((1.0,), (5.0,)) ≈ 5.0 rtol = 1e-10
            @test bruggeman_effective_property((0.5, 0.5), (10.0, 10.0)) ≈ 10.0 rtol = 1e-6
            @test bruggeman_effective_property((0.0, 1.0), (100.0, 5.0)) ≈ 5.0 rtol = 1e-6
        end

        @testset "strictly between the series and parallel bounds" begin
            v, k = (0.7, 0.3), (1.0, 10.0)
            k_eff = bruggeman_effective_property(v, k)
            @test 1.0 / (v[1] / k[1] + v[2] / k[2]) < k_eff < v[1] * k[1] + v[2] * k[2]
            # The positive root of 2k² − 0.1k − 10 = 0.
            @test k_eff ≈ 2.2612078 atol = 1e-4

            v3, k3 = (0.5, 0.3, 0.2), (1.0, 5.0, 10.0)
            k3_eff = bruggeman_effective_property(v3, k3)
            @test 1.0 / sum(v3 ./ k3) < k3_eff < sum(v3 .* k3)

            # A 10:1 contrast: well inside the bounds.
            @test 3.0 < bruggeman_effective_property((0.6, 0.4), (10.0, 1.0)) < 6.0
        end

        @testset "symmetric under relabelling, not under swapping the fractions" begin
            a = bruggeman_effective_property((0.3, 0.7), (1.0, 10.0))
            @test bruggeman_effective_property((0.7, 0.3), (10.0, 1.0)) ≈ a rtol = 1e-6
            @test abs(bruggeman_effective_property((0.7, 0.3), (1.0, 10.0)) - a) > 0.1
        end

        @testset "the root satisfies Σ v(k − k_eff)/(k + 2k_eff) = 0" begin
            v, k = (0.4, 0.6), (2.0, 8.0)
            k_eff = bruggeman_effective_property(v, k)
            @test abs(sum(v[i] * (k[i] - k_eff) / (k[i] + 2k_eff) for i = 1:2)) < 1e-8
        end

        @testset "an insulating phase: k_c(3p − 1)/2 above p = 1/3, zero below" begin
            @test bruggeman_effective_property((0.5, 0.5), (1.0, 0.0)) ≈ 0.25 atol = 1e-8
            @test bruggeman_effective_property((0.6, 0.4), (2.0, 0.0)) ≈ 0.8 atol = 1e-8
            @test bruggeman_effective_property((0.3, 0.7), (1.0, 0.0)) == 0.0
            @test bruggeman_effective_property((0.1, 0.9), (5.0, 0.0)) == 0.0
            @test bruggeman_effective_property((1 / 3, 2 / 3), (1.0, 0.0)) ≈ 0.0 atol = 1e-6
            # Three phases, one insulating: the insulator enters as −v/2.
            v, k = (0.4, 0.3, 0.3), (1.0, 5.0, 0.0)
            k_eff = bruggeman_effective_property(v, k)
            @test k_eff > 0.0
            r =
                v[1] * (k[1] - k_eff) / (k[1] + 2k_eff) +
                v[2] * (k[2] - k_eff) / (k[2] + 2k_eff) - 0.5 * v[3]
            @test abs(r) < 1e-8
        end

        @testset "the conductivity rule lies between the bounds; the gas-transfer rule refuses it" begin
            brugg = Material(
                name = :brugg_test,
                components = (
                    SolidComponent(:comp1, ρ = 1000.0, c = 1500.0, k = 1.0, κ = 0.0),
                    SolidComponent(:comp2, ρ = 500.0, c = 1000.0, k = 10.0, κ = 0.0),
                ),
                reactions = (),
                mixing = MaterialMixing(
                    k = MixingSpec(BRUGGEMAN, 0.5),
                    λ = MixingSpec(PARALLEL, 0.5),
                ),
            )
            # Volume fractions 1/3 and 2/3: series 2.5, parallel 7.0.
            @test 2.5 < effective_conductivity((500.0, 500.0), 400.0, brugg) < 7.0

            brugg_gas = Material(
                name = :brugg_gas_test,
                components = (
                    SolidComponent(
                        :solid,
                        ρ = 1000.0,
                        c = 1500.0,
                        k = 0.2,
                        κ = 0.0,
                        λ = 1e-6,
                    ),
                    GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.03, λ = 1e-4),
                ),
                reactions = (),
                mixing = MaterialMixing(λ = MixingSpec(BRUGGEMAN)),
            )
            @test_throws ArgumentError effective_gas_transfer(
                (800.0, 0.5),
                400.0,
                brugg_gas,
            )
        end
    end
end
