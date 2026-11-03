"""Condensed components hold volume and mix by bulk occupancy, gases never do;
the transport matrix is the solids and the liquids that declare λ or κ; the
dry pyrolysis gases are those the solids give off."""

using Test
using Pyrolysis
using Pyrolysis.Materials: dry_gas_flags, participates_in_transport
using Pyrolysis.Properties:
    transport_matrix_state,
    skeletal_volume_fraction,
    volume_fractions,
    condensed_volume_fractions,
    mixing_volume_fractions,
    effective_conductivity,
    mixture_emissivity,
    mixture_absorptivity
using Pyrolysis.Physics: volume_change_rate
using Pyrolysis.Closures: LateralShrinkage

# A solid, a liquid with its own surface absorptivity, and a gas. The liquid
# declares λ on request.
function _participation_material(; liquid_λ = nothing, k_rule = WEIGHTED)
    Material(
        name = :participation,
        components = (
            SolidComponent(
                :s,
                ρ = 500.0,
                ρ_intrinsic = 1500.0,
                c = 1000.0,
                k = 0.30,
                ε = 0.90,
                λ = 1e-5,
            ),
            LiquidComponent(
                :w,
                ρ = 1000.0,
                c = 4186.0,
                k = 0.60,
                ε = 0.95,
                absorptivity = 0.5,
                λ = liquid_λ,
            ),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03),
        ),
        reactions = (),
        mixing = MaterialMixing(k = k_rule),
    )
end

@testset "condensed components hold their bulk volume, gases none" begin
    mat = _participation_material()
    T = 300.0
    # θ = Σ dξ/ρ over the solid and the liquid; a gas source adds nothing.
    @test volume_change_rate((-10.0, 0.0, 10.0), T, mat) == -10.0 / 500.0
    @test volume_change_rate((0.0, -10.0, 10.0), T, mat) == -10.0 / 1000.0
    @test volume_change_rate((-10.0, 5.0, 5.0), T, mat) ≈ -10.0 / 500.0 + 5.0 / 1000.0
    @test volume_change_rate((0.0, 0.0, 10.0), T, mat) == 0.0

    # The mixing weights are the same bulk occupancies, normalized.
    ξ = (250.0, 50.0, 1.0)  # v_s = 0.5, v_w = 0.05
    v = condensed_volume_fractions(ξ, T, mat)
    @test v[1] ≈ 0.5 / 0.55
    @test v[2] ≈ 0.05 / 0.55
    @test v[3] == 0.0
    @test mixing_volume_fractions(ξ, T, mat) == v
    # However much gas the pores hold, it carries no weight.
    @test condensed_volume_fractions((250.0, 50.0, 1e3), T, mat) == v

    # Skeletal occupancy reads the intrinsic densities.
    @test skeletal_volume_fraction(ξ, T, mat) ≈ 250.0 / 1500.0 + 50.0 / 1000.0
end

@testset "a liquid is weighted in k, ε and absorptivity" begin
    mat = _participation_material(k_rule = PARALLEL)
    T = 300.0
    ξ = (250.0, 50.0, 1.0)
    v = (0.5 / 0.55, 0.05 / 0.55)
    @test effective_conductivity(ξ, T, mat) ≈ v[1] * 0.30 + v[2] * 0.60 rtol = 1e-14
    @test mixture_emissivity(ξ, T, mat) ≈ v[1] * 0.90 + v[2] * 0.95 rtol = 1e-14
    # The solid follows Kirchhoff (a = ε = 0.9), the liquid sets a = 0.5.
    @test mixture_absorptivity(ξ, T, mat) ≈ v[1] * 0.90 + v[2] * 0.5 rtol = 1e-14

    # The weighted default rule mixes the same two condensed conductivities.
    mat_w = _participation_material()
    k_par = v[1] * 0.30 + v[2] * 0.60
    k_ser = 1 / (v[1] / 0.30 + v[2] / 0.60)
    @test effective_conductivity(ξ, T, mat_w) ≈ 0.5 * k_par + 0.5 * k_ser rtol = 1e-14
end

@testset "the depleted-cell fallback" begin
    mat = _participation_material()
    T = 600.0
    # A liquid alone is still a matrix: its volume share 0.05 is far above
    # the blend window, so it carries the whole weight.
    @test mixing_volume_fractions((0.0, 50.0, 1.0), T, mat) == (0.0, 1.0, 0.0)
    @test effective_conductivity((0.0, 50.0, 1.0), T, mat) ≈ 0.60 rtol = 1e-14

    # Halfway through the window (matrix fraction 5.5e-4) the weights are the
    # mean of the matrix fractions and the all-phase occupancy.
    ξ_mid = (0.0, 0.55, 1.0)
    w = mixing_volume_fractions(ξ_mid, T, mat)
    full = volume_fractions(ξ_mid, T, mat)
    @test collect(w) ≈ collect(0.5 .* (0.0, 1.0, 0.0) .+ 0.5 .* full) rtol = 1e-12

    # Only pore gas left: the all-phase occupancy, so k_eff is the gas's.
    ξ_void = (0.0, 0.0, 1.0)
    @test mixing_volume_fractions(ξ_void, T, mat)[3] ≈ 1.0
    @test effective_conductivity(ξ_void, T, mat) ≈ 0.03 rtol = 1e-10
    # A transparent gas has ε = 0; the depletion callback owns this case.
    @test mixture_emissivity(ξ_void, T, mat) == 0.0
end

@testset "transport participation" begin
    liquid(; kw...) = LiquidComponent(:l; ρ = 1000.0, c = 4186.0, k = 0.6, kw...)
    @test participates_in_transport(SolidComponent(:s, ρ = 500.0, c = 1000.0, k = 0.2))
    @test !participates_in_transport(GaseousComponent(:g, M = 0.05, c = 1.0, k = 0.03))
    # A liquid joins the matrix by declaring λ (even zero) or a permeability.
    @test !participates_in_transport(liquid())
    @test !participates_in_transport(liquid(κ = 0.0))
    @test participates_in_transport(liquid(λ = 0.0))
    @test participates_in_transport(liquid(λ = 2e-5))
    @test participates_in_transport(liquid(κ = 1e-14))

    T = 300.0
    ξ = (250.0, 50.0, 1.0)
    q, w = transport_matrix_state(ξ, T, _participation_material())
    @test q == 0.5
    @test w == (1.0, 0.0, 0.0)
    q, w = transport_matrix_state(ξ, T, _participation_material(liquid_λ = 2e-5))
    @test q ≈ 0.55
    @test collect(w) ≈ [0.5 / 0.55, 0.05 / 0.55, 0.0]

    # Under the lumped law a participant needs explicit λ: a solid cannot opt
    # out, and a liquid that declares only κ must supply it.
    gas = GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03)
    s = SolidComponent(:s, ρ = 500.0, c = 1000.0, k = 0.2, λ = 1e-5)
    bare = SolidComponent(:bare, ρ = 500.0, c = 1000.0, k = 0.2)
    @test_throws ArgumentError Material(name = :m, components = (s, bare, gas))
    @test_throws ArgumentError Material(name = :m, components = (s, liquid(κ = 1e-14), gas))
    @test Material(name = :m, components = (s, liquid(), gas)) isa Material
end

@testset "dry pyrolysis gases descend from the solids" begin
    wood = SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, λ = 1e-5)
    char = SolidComponent(:char, ρ = 150.0, c = 1100.0, k = 0.08, λ = 1e-5)
    moisture = LiquidComponent(:moisture, ρ = 1000.0, c = 4186.0, k = 0.6)
    volatiles = GaseousComponent(:volatiles, M = 0.05, c = 1500.0, k = 0.03)
    water = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025)
    vapour = GaseousComponent(:water_vapour, M = 0.018, c = 2000.0, k = 0.025)
    air = GaseousComponent(:AIR, M = 0.029, c = 1006.0, k = 0.026)
    pyrolysis =
        Reaction(:pyrolysis, 1 => (2, 4); A = 1.0, E = 0.0, h = 0.0, yields = (0.2, 0.8))
    dehydration =
        Reaction(:dehydration, 1 => (2, 5); A = 1.0, E = 0.0, h = 0.0, yields = (0.9, 0.1))
    drying = Reaction(:drying, 3 => 6; A = 1.0, E = 0.0, h = 0.0)
    law = χ -> 1.0 - 0.1χ
    wood_mat = Material(
        name = :lineage,
        components = (wood, char, moisture, volatiles, water, vapour, air),
        reactions = (pyrolysis, dehydration, drying),
        lateral_shrinkage_law = law,
    )
    # Water the wood gives off is dry gas whatever its name; the vapour of
    # the moisture and a gas no reaction makes are not.
    flags = (false, false, false, true, true, false, false)
    @test dry_gas_flags(wood_mat) == flags
    @test LateralShrinkage(wood_mat, 1.0, 1.0).dry_gas == flags

    # The lineage passes through a liquid product: polymer → melt → monomer,
    # listed in the order that needs a second pass.
    polymer = SolidComponent(:polymer, ρ = 1190.0, c = 1420.0, k = 0.2, λ = 1e-5)
    melt = LiquidComponent(:melt, ρ = 1100.0, c = 1800.0, k = 0.2)
    monomer = GaseousComponent(:monomer, M = 0.1, c = 1500.0, k = 0.02)
    melt_mat = Material(
        name = :melt,
        components = (polymer, melt, monomer, moisture, vapour),
        reactions = (
            Reaction(:boil, 2 => 3; A = 1.0, E = 0.0, h = 0.0),
            Reaction(:melting, 1 => 2; A = 1.0, E = 0.0, h = 0.0),
            Reaction(:drying, 4 => 5; A = 1.0, E = 0.0, h = 0.0),
        ),
    )
    @test dry_gas_flags(melt_mat) == (false, false, true, false, false)

    # A gas made by both lineages is refused, naming it and the reactions.
    mixed = Material(
        name = :mixed,
        components = (wood, char, moisture, volatiles, water),
        reactions = (dehydration, Reaction(:drying, 3 => 5; A = 1.0, E = 0.0, h = 0.0)),
        lateral_shrinkage_law = law,
    )
    err = try
        dry_gas_flags(mixed)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin(":H2O", err.msg)
    @test occursin(":dehydration", err.msg)
    @test occursin(":drying", err.msg)
    @test_throws ArgumentError LateralShrinkage(mixed, 1.0, 1.0)
end
