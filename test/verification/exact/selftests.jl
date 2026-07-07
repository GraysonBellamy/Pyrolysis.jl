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

    @testset "slab convective eigenseries (C&J §3.11)" begin
        # 1. Eigencondition satisfied to bisection precision, roots in-branch.
        for Bi in (0.1, 1.0, 10.0, 100.0)
            λs = biot_eigenvalues(Bi, 60)
            @test all(abs(λ * sin(λ) - Bi * cos(λ)) < 1e-10 * max(1.0, λ * Bi) for λ in λs)
            @test all((n - 1) * π < λs[n] < (n - 1) * π + π / 2 for n in eachindex(λs))
        end

        # 2. Limits of λ₁: Bi → ∞ gives π/2 (Dirichlet); Bi → 0 gives √Bi (lumped).
        @test biot_eigenvalues(1e8, 1)[1] ≈ π / 2 rtol = 1e-6
        @test biot_eigenvalues(1e-4, 1)[1] ≈ sqrt(1e-4) rtol = 1e-3

        # 3. Initial condition recovered: θ/θᵢ → 1 as Fo → 0 (interior point).
        λs1 = biot_eigenvalues(1.0, 60)
        @test slab_convective_theta(0.5, 0.002, λs1) ≈ 1.0 atol = 1e-2

        # 4. Lumped limit: for Bi ≪ 1, θ ≈ exp(−Bi·Fo).
        λs_small = biot_eigenvalues(0.01, 60)
        @test slab_convective_theta(0.5, 10.0, λs_small) ≈ exp(-0.01 * 10.0) rtol = 1e-2

        # 5. Dirichlet limit: for Bi ≫ 1 matches the λₙ = (n−1/2)π series.
        λs_big = biot_eigenvalues(1e10, 60)
        θ_ref = sum(
            4 * sin(λ) / (2λ + sin(2λ)) * exp(-λ^2 * 0.1) * cos(λ * 0.3) for
            λ in [(n - 0.5) * π for n = 1:60]
        )
        @test slab_convective_theta(0.3, 0.1, λs_big) ≈ θ_ref rtol = 1e-6
    end

    @testset "isothermal kinetics closed forms" begin
        # 1. nth-order: dξ/dt = −kξⁿ verified by central finite difference.
        for n in (0.0, 2.0, 2 / 3)
            ξ0, k, t = 100.0, 0.5, 5.0
            δ = 1e-6
            dξdt =
                (
                    nth_order_xi(t + δ; ξ0 = ξ0, k = k, n = n) -
                    nth_order_xi(t - δ; ξ0 = ξ0, k = k, n = n)
                ) / (2δ)
            ξt = nth_order_xi(t; ξ0 = ξ0, k = k, n = n)
            @test dξdt ≈ -k * ξt^n rtol = 1e-6
        end

        # 2. Burnout: exact zero at t_b for n < 1, and m^(1/3) linear in t for
        #    the contracting-sphere order n = 2/3 (FDS NaHCO₃ marker).
        tb = nth_order_burnout_time(ξ0 = 100.0, k = 5.0, n = 0.0)
        @test tb == 20.0
        @test nth_order_xi(tb; ξ0 = 100.0, k = 5.0, n = 0.0) == 0.0
        tb3 = nth_order_burnout_time(ξ0 = 1.0, k = 0.3, n = 2 / 3)
        m13(t) = nth_order_xi(t; ξ0 = 1.0, k = 0.3, n = 2 / 3)^(1 / 3)
        @test m13(0.3 * tb3) ≈ 1 - 0.3 rtol = 1e-12
        @test m13(0.7 * tb3) ≈ 1 - 0.7 rtol = 1e-12

        # 3. Series chain: ODE residual by finite difference + conservation +
        #    B peaks at t* = ln(k₁/k₂)/(k₁ − k₂).
        ξ0, k1, k2 = 1000.0, 0.389, 0.262
        δ = 1e-6
        for t in (0.5, 2.0, 8.0)
            ξA, ξB, _ = series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)
            ξBp = series_chain_xi(t + δ; ξ0 = ξ0, k1 = k1, k2 = k2)[2]
            ξBm = series_chain_xi(t - δ; ξ0 = ξ0, k1 = k1, k2 = k2)[2]
            @test (ξBp - ξBm) / (2δ) ≈ k1 * ξA - k2 * ξB rtol = 1e-6
            @test sum(series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)) ≈ ξ0 rtol = 1e-14
        end
        tpeak = log(k1 / k2) / (k1 - k2)
        ξB_at(t) = series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)[2]
        @test ξB_at(tpeak) > ξB_at(tpeak - 0.01)
        @test ξB_at(tpeak) > ξB_at(tpeak + 0.01)
    end

    @testset "radiative cooling implicit solution" begin
        T∞, Ti = 300.0, 900.0

        # 1. G against composite-Simpson quadrature of ∫dT/(T⁴−T∞⁴) — the
        #    plan flags this antiderivative as sign-error-prone; this is the
        #    mandated independent check.
        for T_hi in (500.0, 700.0, 900.0)
            T_lo = 400.0
            n = 4000
            hstep = (T_hi - T_lo) / n
            f(T) = 1 / (T^4 - T∞^4)
            simpson =
                hstep / 3 * (
                    f(T_lo) +
                    f(T_hi) +
                    4 * sum(f(T_lo + hstep * i) for i = 1:2:(n-1)) +
                    2 * sum(f(T_lo + hstep * i) for i = 2:2:(n-2))
                )
            @test radiative_cooling_G(T_hi, T∞) - radiative_cooling_G(T_lo, T∞) ≈ simpson rtol =
                1e-10
        end

        # 2. Round trip: T(t(T)) = T through the bisection inverse.
        K = 5e-10
        for T in (850.0, 600.0, 400.0)
            t = (radiative_cooling_G(Ti, T∞) - radiative_cooling_G(T, T∞)) / K
            @test radiative_cooling_T(t; Ti = Ti, T∞ = T∞, K = K) ≈ T rtol = 1e-10
        end

        # 3. Early-time linearization: T ≈ Tᵢ − K(Tᵢ⁴−T∞⁴)t.
        t_small = 1e-4 / (K * Ti^3)
        T_lin = Ti - K * (Ti^4 - T∞^4) * t_small
        @test radiative_cooling_T(t_small; Ti = Ti, T∞ = T∞, K = K) ≈ T_lin rtol = 1e-3

        # 4. Steady balance root: residual is zero and root is bracketed
        #    sensibly (~732 K for these cone-like parameters).
        σ = 5.670374419e-8
        Ts = radiative_convective_equilibrium_Ts(
            q_abs = 20e3,
            ε = 0.85,
            h = 15.0,
            T∞ = 300.0,
            σ = σ,
        )
        @test abs(20e3 - 0.85 * σ * (Ts^4 - 300.0^4) - 15.0 * (Ts - 300.0)) < 1e-6
        @test 700.0 < Ts < 760.0
    end

    @testset "slab constant-flux series (C&J §3.8)" begin
        # 1. Fourier identity: φ(x̃, 0) = 0 for all x̃ (the series is exactly
        #    the cosine expansion of −(3x̃²−1)/6). Must be evaluated AT Fo = 0:
        #    for any Fo > 0 the face value physically rises like 2√(Fo/π).
        #    Residual is the 5000-term truncation tail, ~4e-5 at the face.
        for x̃ in (0.0, 0.25, 0.5, 0.75, 1.0)
            @test abs(slab_constflux_phi(x̃, 0.0)) < 1e-3
        end

        # 2. Front–back plateau: φ(1) − φ(0) → 1/2 for Fo ≫ 1.
        @test slab_constflux_phi(1.0, 2.0) - slab_constflux_phi(0.0, 2.0) ≈ 0.5 rtol = 1e-8

        # 3. Global energy: spatial mean ⟨φ⟩ = Fo exactly (trapezoid on a fine
        #    grid; the parabola and every cosine integrate to zero).
        for Fo in (0.05, 0.3, 1.0)
            xs = range(0.0, 1.0, length = 2001)
            φs = [slab_constflux_phi(x, Fo) for x in xs]
            mean_φ = (sum(φs) - 0.5 * (φs[1] + φs[end])) / (length(xs) - 1)
            @test mean_φ ≈ Fo atol = 1e-6
        end
    end
end
