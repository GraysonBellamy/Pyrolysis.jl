"""
Regression tests for the swelling-factor (γ) semantics.

- The swelling factor γ has ONE semantics everywhere — every component
  contributes v_j = γ_j·ξ_j/ρ_j to the dilation rate θ, the mixture
  density, and the conductivity/emissivity volume weights. Pre-fix the
  three kernels disagreed for partial-γ condensed phases (γ applied to
  every component in θ, to gases only in the density, and as a binary
  gas-exclusion predicate in the mixing weights).
- A fully-converted cell no longer becomes a perfect non-emitting
  insulator: k_eff falls back to the gas-mixture conductivity. Liquids
  default to a realistic emissivity (0.95) instead of 0.
- γ range validation at construction.
"""

using Test
using Pyrolysis
using Pyrolysis.Properties:
    mixture_density,
    solid_volume_fractions,
    mixing_volume_fractions,
    effective_conductivity,
    effective_emissivity
using Pyrolysis.Physics: volume_change_rate

function _gamma_material(; γ_liq = 0.0, ε_liq = nothing)
    liq_kwargs = ε_liq === nothing ? (;) : (; ε = ε_liq)
    Material(
        name = :gamma_test,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.30, ε = 0.90),
            LiquidComponent(:w, ρ = 1000.0, c = 4186.0, k = 0.60, γ = γ_liq, liq_kwargs...),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )
end

@testset "one γ semantics across θ, ρ_mix, and mixing weights" begin
    T = 300.0
    ξ = (500.0, 500.0, 0.0)

    # γ = 0 liquid: consistently volumeless.
    mat0 = _gamma_material(γ_liq = 0.0)
    # mixture density: only the solid's volume counts. Y_s = 0.5,
    # 1/ρ = 0.5/1000 → ρ = 2000 (pre-fix: condensed always contributed → 1000).
    @test mixture_density(ξ, T, mat0) ≈ 2000.0
    # mixing weights: all on the solid (pre-fix: split 50/50 with the liquid).
    @test solid_volume_fractions(ξ, T, mat0)[1] ≈ 1.0
    @test solid_volume_fractions(ξ, T, mat0)[2] == 0.0
    # θ: evaporating liquid adds no strain (this was already γ-weighted).
    @test volume_change_rate((0.0, -10.0, 10.0), T, mat0) == 0.0

    # γ = 1 liquid: full contribution in all three kernels.
    mat1 = _gamma_material(γ_liq = 1.0)
    @test mixture_density(ξ, T, mat1) ≈ 1000.0
    v1 = solid_volume_fractions(ξ, T, mat1)
    @test v1[1] ≈ 0.5 && v1[2] ≈ 0.5
    @test volume_change_rate((0.0, -10.0, 10.0), T, mat1) ≈ -10.0 / 1000.0

    # Partial γ: the SAME fraction shows up in all three kernels —
    # the defining unification.
    matp = _gamma_material(γ_liq = 0.5)
    # 1/ρ = 0.5/1000 + 0.5·0.5/1000 = 0.75e-3
    @test mixture_density(ξ, T, matp) ≈ 1000.0 / 0.75
    vp = solid_volume_fractions(ξ, T, matp)
    @test vp[1] ≈ 0.5 / 0.75 && vp[2] ≈ 0.25 / 0.75
    @test volume_change_rate((0.0, -10.0, 10.0), T, matp) ≈ 0.5 * (-10.0) / 1000.0

    # γ is validated to [0, 1] at construction.
    @test_throws ArgumentError SolidComponent(:bad, ρ = 1.0, c = 1.0, k = 1.0, γ = 1.5)
    @test_throws ArgumentError GaseousComponent(
        :bad,
        M = 0.05,
        c = 1.0,
        k = 1.0,
        λ = 1e-5,
        γ = -0.1,
    )
    @test_throws ArgumentError LiquidComponent(:bad, ρ = 1.0, c = 1.0, k = 1.0, γ = 2.0)
end

@testset "depleted-cell floor and liquid ε default" begin
    mat = _gamma_material()
    T = 600.0

    # Healthy cell: unchanged behavior (weights on the matrix).
    ξ_ok = (500.0, 0.0, 1.0)
    @test effective_conductivity(ξ_ok, T, mat) > 0.1
    @test effective_emissivity(ξ_ok, T, mat) ≈ 0.90

    # Fully converted cell (only pore gas left): pre-fix k_eff = ε_eff = 0 —
    # a perfect non-emitting insulator whose zero also wiped both adjacent
    # harmonic face means. Post-fix k_eff degrades to the gas conductivity.
    ξ_void = (0.0, 0.0, 1.0)
    @test mixing_volume_fractions(ξ_void, T, mat)[3] ≈ 1.0
    @test effective_conductivity(ξ_void, T, mat) ≈ 0.03 rtol = 1e-10
    # ε of a pure default-gas cell is genuinely 0 (transparent gas) — the
    # depletion callback owns this case; documented, not floored.
    @test effective_emissivity(ξ_void, T, mat) == 0.0

    # Liquid emissivity defaults realistically (water ≈ 0.96): a γ = 1
    # surface liquid must RAISE ε_eff toward water, not drag it to 0.
    liq = LiquidComponent(:w, ρ = 1000.0, c = 4186.0, k = 0.6)
    @test liq.emissivity(300.0) ≈ 0.95
    mat_wet = Material(
        name = :wet,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.30, ε = 0.90),
            LiquidComponent(:w, ρ = 1000.0, c = 4186.0, k = 0.60, γ = 1.0),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )
    ε_wet = effective_emissivity((900.0, 100.0, 0.0), 300.0, mat_wet)
    @test ε_wet > 0.90   # pre-fix (liquid ε = 0): 0.9·0.9 + 0.1·0 = 0.81
end
