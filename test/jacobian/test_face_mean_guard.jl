"""
Regression tests: the residual and the structured-Jacobian operators must
share ONE face-mean guard band.

Pre-fix: the residual's `harmonic_mean_weighted` cut at 1e-20 (either side
below → face value 0) while the operators' `_harmonic_mean_with_derivs` cut
at 1e-30, so J ≠ residual whenever a conductivity fell in [1e-30, 1e-20).
Sharper: `DarcyTransportOperator` carried a private `_harmonic_λ_face` with a
one-sided fallback the residual's λ-face never had — disagreement whenever
exactly one side's λ was below the cut.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Properties: harmonic_mean_weighted
using Pyrolysis.Jacobian: _harmonic_mean_with_derivs

@testset "one face-mean guard band" begin
    d_L, d_R = 1e-3, 2e-3

    # Inside the old disagreement band [1e-30, 1e-20): both sides return 0.
    for k_small in (1e-25, 1e-21)
        @test harmonic_mean_weighted(k_small, 0.2, d_L, d_R) == 0.0
        k_face, dkL, dkR = _harmonic_mean_with_derivs(k_small, 0.2, d_L, d_R)
        @test k_face == 0.0 && dkL == 0.0 && dkR == 0.0
    end

    # Regular regime: identical values, and the derivatives match the
    # analytic harmonic-mean sensitivities.
    k_L, k_R = 0.15, 0.45
    k_face, dkL, dkR = _harmonic_mean_with_derivs(k_L, k_R, d_L, d_R)
    @test k_face ≈ harmonic_mean_weighted(k_L, k_R, d_L, d_R) rtol = 1e-14
    R = d_L / k_L + d_R / k_R
    @test dkL ≈ (d_L + d_R) * d_L / (k_L^2 * R^2) rtol = 1e-12
    @test dkR ≈ (d_L + d_R) * d_R / (k_R^2 * R^2) rtol = 1e-12

    # One-sided low: the residual zeroes the face; the operator's old private
    # helper returned the surviving side instead. The helper is gone — the
    # operator now calls harmonic_mean_weighted itself.
    @test harmonic_mean_weighted(0.0, 1e-5, d_L, d_R) == 0.0
    @test !isdefined(Pyrolysis.Jacobian, :_harmonic_λ_face)

    # The guard branch returns zero of the promoted argument type, as the
    # regular branch does.
    dual = ForwardDiff.Dual(1e-25, 1.0)
    @test harmonic_mean_weighted(dual, 0.2, d_L, d_R) isa typeof(dual)
    @test harmonic_mean_weighted(dual, 0.2, d_L, d_R) == 0
    @test harmonic_mean_weighted(0.0f0, 0.2f0, 1.0f-3, 2.0f-3) === 0.0f0
end
