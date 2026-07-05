module WPKirchhoffFixture
include(joinpath(@__DIR__, "..", "..", "examples", "inverse_analysis", "wp_properties.jl"))
end

using Test
using ForwardDiff

@testset "White pine Kirchhoff heat lift" begin
    WP = WPKirchhoffFixture

    for i = 1:5
        h_lifted = WP.WP_KINETICS_KIRCHHOFF_ALL[i].h
        h_dsc = WP.WP_KINETICS[i].h
        @test h_lifted(WP.WP_T_PEAK_DSC[i]) ≈ h_dsc atol = 1e-8 rtol = 0
    end

    for i = 1:3
        h_lifted = WP.WP_KINETICS_KIRCHHOFF_123[i].h
        h_dsc = WP.WP_KINETICS[i].h
        @test h_lifted(WP.WP_T_PEAK_DSC[i]) ≈ h_dsc atol = 1e-8 rtol = 0
    end
    @test WP.WP_KINETICS_KIRCHHOFF_123[4].h == WP.WP_KINETICS[4].h
    @test WP.WP_KINETICS_KIRCHHOFF_123[5].h == WP.WP_KINETICS[5].h

    @test WP.WP_KINETICS_KIRCHHOFF_ALL[1].h(700.0) < WP.WP_KINETICS[1].h
    @test 0.0 < WP.WP_KINETICS_KIRCHHOFF_ALL[3].h(700.0) < WP.WP_KINETICS[3].h
    @test WP.WP_KINETICS_KIRCHHOFF_ALL[4].h(800.0) < WP.WP_KINETICS[4].h

    h3 = WP.WP_KINETICS_KIRCHHOFF_ALL[3].h
    dh3_dT = ForwardDiff.derivative(h3, 700.0)
    h3_dual = h3(ForwardDiff.Dual(700.0, 1.0))
    @test isfinite(h3(700.0))
    @test isfinite(dh3_dT)
    @test dh3_dT < 0.0
    @test ForwardDiff.value(h3_dual) ≈ h3(700.0)
end
