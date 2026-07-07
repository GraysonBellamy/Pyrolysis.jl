# =============================================================================
# Self-tests of the exact-solution library
# =============================================================================
#
# Each exact solution is verified against an independent check before any
# solver output is compared to it (formula-provenance rule). These run as
# part of the verification suite.

@testset "Exact-solution library self-tests" begin
    @testset "semi-infinite convective erfc (C&J §2.7)" begin
        # Gpyro §3.4.1 parameters
        h, k = 20.0, 0.2
        α = 0.2 / (1000.0 * 1400.0)

        # 1. Stable erfcx form ≡ naive textbook form where the latter is
        #    representable (exponent small).
        for (x, t) in ((0.0, 1.0), (0.001, 10.0), (0.005, 100.0), (0.02, 180.0))
            η = x / (2 * sqrt(α * t))
            B = η + h * sqrt(α * t) / k
            naive = erfc(η) - exp(h * x / k + h^2 * α * t / k^2) * erfc(B)
            @test semiinf_convective_theta(x, t; h = h, k = k, α = α) ≈ naive rtol = 1e-12
        end

        # 2. h → ∞ recovers the Dirichlet solution θ = erfc(η).
        θ_dirichlet = semiinf_convective_theta(0.003, 60.0; h = 1e12, k = k, α = α)
        @test θ_dirichlet ≈ erfc(0.003 / (2 * sqrt(α * 60.0))) rtol = 1e-6

        # 3. Small-time surface expansion: θ_s = 1 − erfcx(a) ≈ 2a/√π,
        #    a = h√(αt)/k (relative error O(a·√π/2)).
        t_small = 1e-6
        a = h * sqrt(α * t_small) / k
        θ_s = semiinf_convective_theta(0.0, t_small; h = h, k = k, α = α)
        @test θ_s ≈ 2 * a / sqrt(π) rtol = 1e-3

        # 4. No overflow in the regime that kills the naive form
        #    (exponent h²αt/k² ~ 1e5): finite and bounded in [0, 1].
        θ_extreme = semiinf_convective_theta(1e-4, 100.0; h = 2e4, k = k, α = α)
        @test isfinite(θ_extreme) && 0.0 <= θ_extreme <= 1.0

        # 5. Monotonicity: θ decreases with depth, increases with time.
        θs = [semiinf_convective_theta(x, 60.0; h = h, k = k, α = α) for x = 0.0:0.002:0.02]
        @test issorted(θs, rev = true)
        θt = [
            semiinf_convective_theta(0.005, t; h = h, k = k, α = α) for t = 10.0:10.0:180.0
        ]
        @test issorted(θt)
    end
end
