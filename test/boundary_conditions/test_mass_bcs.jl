"""
Mass boundary conditions of the shared layer: the pressure conditions'
defaults, the film law of `DiffusiveMassBC` and the ambient concentrations of
`ConvectiveMassBC`.
"""

using Test
using Pyrolysis
using Pyrolysis.BoundaryConditions: mass_bc_flux, ambient_concentration

@testset "Mass boundary conditions" begin
    @testset "pressure conditions: atmospheric and impermeable defaults, given values kept" begin
        @test DirichletPressureBC().P == 101325.0
        @test DirichletPressureBC(50000.0).P == 50000.0
        @test NeumannPressureBC().∂P∂n == 0.0
        @test NeumannPressureBC(1000.0).∂P∂n == 1000.0
        bc = ConvectivePressureBC(1e-6, 101325.0)
        @test bc.h_P == 1e-6
        @test bc.P_ext == 101325.0
    end

    @testset "the diffusive film flux reads λ, not λρ" begin
        # Without a film the face flux is −λξ/d for any carrier density, the
        # vanishing ones included; with a film coefficient the two
        # resistances add in series.
        λ, ξ, d = 1e-5, 2e-3, 1e-4
        for ρ in (1.0, 1e-20, 1e-30)
            @test mass_bc_flux(DiffusiveMassBC(), ξ, λ * ρ, ρ, d, 0.0, 1) ≈ -λ * ξ / d
        end
        h_m = 0.01
        @test mass_bc_flux(DiffusiveMassBC(h_m), ξ, λ * 1e-30, 1e-30, d, 0.0, 1) ≈
              -λ * h_m * ξ / (λ + h_m * d)
    end

    @testset "ambient concentrations: one for every species, one per species, or a function of time" begin
        scalar = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.3)
        @test ambient_concentration(scalar, 0.0, 1) == 0.3
        @test ambient_concentration(scalar, 0.0, 2) == 0.3

        per_species = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.27))
        @test ambient_concentration(per_species, 0.0, 1) == 0.0
        @test ambient_concentration(per_species, 0.0, 2) == 0.27

        timed = ConvectiveMassBC(h_m = 0.05, ξ_inf = t -> 0.1 * (1 + t))
        @test ambient_concentration(timed, 1.0, 2) ≈ 0.2

        # A negative ambient concentration is refused at construction.
        @test_throws ArgumentError ConvectiveMassBC(h_m = 0.05, ξ_inf = (-0.1, 0.0))
    end
end
