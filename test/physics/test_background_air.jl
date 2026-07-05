"""
Regression tests for background pore air. A DarcyFick problem initialized
with all gas ξⱼ = 0 starts as a vacuum against the 101325 Pa boundary value,
so early product gas is spuriously advected inward down the
absolute-pressure slope instead of being resisted by ~1 atm of pore air.

The remedy is `background_gas_concentration`: the single-species inverse of
the pressure closure P = Σⱼ ξⱼRT/(Mⱼφ), used by the DarcyFick examples to
initialize an inert air component at ambient pressure. These tests pin the
helper's EOS round-trip and the equilibrium-hold property (uniform air at
1 atm ⇒ uniform pressure ⇒ zero Darcy driving force).
"""

using Test
using Pyrolysis
using Pyrolysis.Physics: ideal_gas_pressure
using Pyrolysis.Properties: porosity

function _air_fixture()
    virgin = SolidComponent(
        :virgin,
        ρ = 500.0,
        c = 2000.0,
        k = 0.15,
        κ = 1e-15,
        ρ_intrinsic = 1500.0,  # φ_virgin-only = 2/3
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 1500.0,
        k = 0.10,
        κ = 1e-12,
        ρ_intrinsic = 1950.0,
    )
    gas = GaseousComponent(:volatiles, M = 0.044, c = 1200.0, k = 0.03, λ = 2e-5)
    air = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
    rxn = Reaction(
        :pyro,
        1 => (2, 3),
        yields = (0.2, 0.8),
        A = 1e13,
        E = 150e3,
        h = 300e3,
        T_min = 400.0,
    )
    return Material(
        name = :background_air_fixture,
        components = (virgin, char, gas, air),
        reactions = (rxn,),
    )
end

@testset "Background pore air" begin
    mat = _air_fixture()
    T0 = 300.0
    ξ_condensed = (500.0, 0.0, 0.0, 0.0)

    @testset "EOS round-trip at ambient" begin
        ξ_air = background_gas_concentration(mat, :AIR, ξ_condensed, T0)
        @test ξ_air > 0.0

        # ξ_air = P·M·φ/(R·T) against hand numbers: φ = 1 − 500/1500 = 2/3
        φ0 = porosity(ξ_condensed, T0, mat)
        @test φ0 ≈ 2 / 3 rtol = 1e-12
        @test ξ_air ≈ 101325.0 * 0.02897 * φ0 / (8.31446261815324 * T0) rtol = 1e-12

        # Pressure closure reproduces exactly 1 atm with that concentration
        ξ_full = (500.0, 0.0, 0.0, ξ_air)
        φ = porosity(ξ_full, T0, mat)
        @test ideal_gas_pressure(ξ_full, T0, φ, mat) ≈ 101325.0 rtol = 1e-12

        # Index form matches name form; custom pressure scales linearly
        @test background_gas_concentration(mat, 4, ξ_condensed, T0) == ξ_air
        @test background_gas_concentration(mat, :AIR, ξ_condensed, T0; P = 2 * 101325.0) ≈
              2 * ξ_air rtol = 1e-12
    end

    @testset "equilibrium hold: uniform air ⇒ no Darcy driving force" begin
        # Uniform condensed composition + uniform air at 1 atm: the cell
        # pressure field is flat, so the interior Darcy driving force is zero
        # and the boundary face sees no pressure jump against a 101325 Pa
        # Dirichlet value. (The vacuum start had ΔP ≈ 1 atm at every face.)
        ξ_air = background_gas_concentration(mat, :AIR, ξ_condensed, T0)
        n = 8
        mesh = create_uniform_mesh(
            0.01,
            n,
            Val(4);
            T_initial = T0,
            ξ_initial = (500.0, 0.0, 0.0, ξ_air),
        )
        P = [
            let ξ_i = mesh.cell_states[i].ξ, T_i = mesh.cell_states[i].T
                ideal_gas_pressure(ξ_i, T_i, porosity(ξ_i, T_i, mat), mat)
            end for i = 1:n
        ]
        @test all(p -> isapprox(p, 101325.0; rtol = 1e-12), P)

        ∇P_max = maximum(
            abs(P[i+1] - P[i]) / (mesh.cell_centers[i+1] - mesh.cell_centers[i]) for
            i = 1:(n-1)
        )
        @test ∇P_max < 1e-6
    end

    @testset "argument validation" begin
        @test_throws ArgumentError background_gas_concentration(
            mat,
            :virgin,
            ξ_condensed,
            T0,
        )
        # unknown component name errors in component_index
        @test_throws Exception background_gas_concentration(mat, :nope, ξ_condensed, T0)
    end

    @testset "no-pore-space degeneracy is visible, not silent" begin
        # A material whose solids carry no ρ_intrinsic runs the legacy
        # "no-pores" model: φ = 0 at full bulk density, so the helper returns
        # exactly 0 — there is no pore volume to hold air. This documents why
        # the DarcyFick examples must set ρ_intrinsic.
        virgin0 = SolidComponent(:virgin, ρ = 500.0, c = 2000.0, k = 0.15, κ = 1e-15)
        gas0 = GaseousComponent(:volatiles, M = 0.044, c = 1200.0, k = 0.03, λ = 2e-5)
        air0 = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
        mat0 =
            Material(name = :no_pores, components = (virgin0, gas0, air0), reactions = ())
        @test background_gas_concentration(mat0, :AIR, (500.0, 0.0, 0.0), 300.0) == 0.0
    end
end
