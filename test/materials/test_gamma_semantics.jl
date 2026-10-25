"""Dilation and matrix weights use γ; skeletal occupancy is independent of γ."""

using Test
using Pyrolysis
using Pyrolysis.Properties:
    transport_matrix_state,
    skeletal_volume_fraction,
    solid_volume_fractions,
    mixing_volume_fractions,
    effective_conductivity,
    mixture_emissivity
using Pyrolysis.Physics: volume_change_rate

function _gamma_material(; γ_liq = 0.0, ε_liq = nothing)
    liq_kwargs = ε_liq === nothing ? (;) : (; ε = ε_liq)
    Material(
        name = :gamma_test,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.30, ε = 0.90, λ = 1e-5),
            LiquidComponent(
                :w,
                ρ = 1000.0,
                c = 4186.0,
                k = 0.60,
                γ = γ_liq,
                λ = 0.0,
                liq_kwargs...,
            ),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03),
        ),
        reactions = (),
    )
end

@testset "γ-weighted matrix and independent skeleton" begin
    T=300.0
    ξ=(500.0, 500.0, 0.0)
    for γ in (0.0, 0.5, 1.0)
        mat=_gamma_material(γ_liq = γ)
        q, w=transport_matrix_state(ξ, T, mat)
        @test q≈0.5*(1+γ)
        @test w[1]≈1/(1+γ)
        @test w[2]≈γ/(1+γ)
        @test w[3]==0
        @test all(isapprox.(solid_volume_fractions(ξ, T, mat), w))
        @test skeletal_volume_fraction(ξ, T, mat)==1.0
        @test volume_change_rate((0.0, -10.0, 10.0), T, mat)≈-γ*10/1000
    end
    # γ is validated to [0, 1] at construction.
    @test_throws ArgumentError SolidComponent(:bad, ρ = 1.0, c = 1.0, k = 1.0, γ = 1.5)
    @test_throws ArgumentError GaseousComponent(:bad, M = 0.05, c = 1.0, k = 1.0, γ = -0.1)
    @test_throws ArgumentError LiquidComponent(:bad, ρ = 1.0, c = 1.0, k = 1.0, γ = 2.0)
end

@testset "depleted-cell floor and liquid ε default" begin
    mat = _gamma_material()
    T = 600.0

    # Healthy cell: unchanged behavior (weights on the matrix).
    ξ_ok = (500.0, 0.0, 1.0)
    @test effective_conductivity(ξ_ok, T, mat) > 0.1
    @test mixture_emissivity(ξ_ok, T, mat) ≈ 0.90

    # Fully converted cell (only pore gas left): k_eff degrades to the gas
    # conductivity. A zero k_eff would make the cell a perfect insulator and
    # wipe both adjacent harmonic face means.
    ξ_void = (0.0, 0.0, 1.0)
    @test mixing_volume_fractions(ξ_void, T, mat)[3] ≈ 1.0
    @test effective_conductivity(ξ_void, T, mat) ≈ 0.03 rtol = 1e-10
    # ε of a pure default-gas cell is genuinely 0 (transparent gas) — the
    # depletion callback owns this case; documented, not floored.
    @test mixture_emissivity(ξ_void, T, mat) == 0.0

    # Liquid emissivity defaults realistically (water ≈ 0.96): a γ = 1
    # surface liquid must RAISE ε_eff toward water, not drag it to 0.
    liq = LiquidComponent(:w, ρ = 1000.0, c = 4186.0, k = 0.6)
    @test liq.emissivity(300.0) ≈ 0.95
    mat_wet = Material(
        name = :wet,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.30, ε = 0.90, λ = 1e-5),
            LiquidComponent(:w, ρ = 1000.0, c = 4186.0, k = 0.60, γ = 1.0, λ = 0.0),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03),
        ),
        reactions = (),
    )
    ε_wet = mixture_emissivity((900.0, 100.0, 0.0), 300.0, mat_wet)
    @test ε_wet > 0.90   # a liquid ε = 0 would give 0.9·0.9 + 0.1·0 = 0.81
    @test mixture_emissivity((900.0, 100.0, 0.0), 300.0, mat_wet) ≈ 0.9 * 0.9 + 0.1 * 0.95 rtol =
        1e-12
end
