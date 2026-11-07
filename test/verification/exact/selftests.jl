# =============================================================================
# Self-tests of the exact-solution library
# =============================================================================
#
# Each exact solution is verified against an independent check before any
# solver output is compared to it (formula-provenance rule). These run as
# part of the verification suite.

@testset "Exact-solution library self-tests" begin
    @testset "Gpyro-inspired independent references" begin
        # Robin uptake must equal the spatial integral of the concentration
        # and its time derivative must equal the imposed surface exchange.
        D, t = 2e-5, 150.0
        for b in (5.0, 500.0)
            response(x, time) = semiinf_convective_theta(x, time; h = b, k = 1.0, α = D)
            xmax, n = 16sqrt(D*t), 4000
            dx = xmax/n
            integral =
                dx/3*(
                    response(0.0, t) +
                    response(xmax, t) +
                    4sum(response(i*dx, t) for i = 1:2:(n-1)) +
                    2sum(response(i*dx, t) for i = 2:2:(n-2))
                )
            @test robin_diffusion_uptake(t; D, b) ≈ integral rtol = 1e-9
            δ = 1e-3
            derivative =
                (robin_diffusion_uptake(t+δ; D, b)-robin_diffusion_uptake(t-δ; D, b))/(2δ)
            @test derivative ≈ D*b*(1-response(0.0, t)) rtol = 1e-8
            # Independent one-sided derivative of the continuum surface BC.
            δx = 1e-7
            gradient = (-3response(0.0, t)+4response(δx, t)-response(2δx, t))/(2δx)
            @test -gradient ≈ b*(1-response(0.0, t)) rtol = 1e-7
        end
        @test robin_diffusion_uptake(0.0; D, b = 5.0) == 0.0

        # Average references are checked by quadrature of the point solution,
        # independently of the primitive used by the new error metric.
        for b in (6.5, 650.0),
            t in (50.0, 150.0, 400.0),
            width in (1/40, 1/160, 1/320),
            xlo in (0.0, 0.03, 0.3, 0.8)

            response(x) = semiinf_convective_theta(x, t; h = b, k = 1.0, α = D)
            n = 128
            dx = width/n
            quadrature =
                (
                    response(xlo) +
                    response(xlo+width) +
                    4sum(response(xlo+i*dx) for i = 1:2:(n-1)) +
                    2sum(response(xlo+i*dx) for i = 2:2:(n-2))
                )/(3n)
            average = robin_diffusion_cell_average(xlo, xlo+width, t; D, b)
            @test average ≈ quadrature rtol = 1e-9 atol = 2e-12
        end
        for b in (6.5, 650.0)
            width = 1/160
            integral = sum(
                robin_diffusion_cell_average(i*width, (i+1)*width, 400.0; D, b)*width
                for i = 0:159
            )
            @test integral ≈ robin_diffusion_uptake(400.0; D, b) rtol = 1e-11
        end
        @test robin_diffusion_cell_average(0.0, 0.01, 0.0; D, b = 5.0) == 0.0
        @test robin_diffusion_cell_average(0.0, 0.01, 50.0; D, b = 0.0) == 0.0

        # Darcy reference satisfies the constitutive law, with pressure and
        # gradient independently evaluated from the exact profile.
        Pout, flux, μ, T, κ, M = 1e5, 10.0, 2e-5, 300.0, 5e-10, 0.029
        p(x) = darcy_imposed_flux_pressure(x; Pout, flux, μ, T, κ, M)
        @test p(0.0) == Pout
        for x in (0.01, 0.05, 0.09)
            δ = 1e-6
            dpdx = (p(x+δ)-p(x-δ))/(2δ)
            @test κ*M/(μ*Pyrolysis.Materials.R_GAS*T)*p(x)*dpdx ≈ flux rtol = 1e-8
        end

        # Differentiate the temperature reference and check the independent
        # adiabatic energy equation, including the equal-capacity limit.
        T0, rate, ca = 1000.0, 0.2, 1200.0
        for cb in (ca, 1900.0, 1360.0), heat in (1e6, -1e5, 0.0)
            temp(t) = adiabatic_reaction_temperature(t; T0, rate, heat, ca, cb)
            @test temp(0.0) == T0
            for t in (0.5, 5.0, 20.0)
                δ = 1e-4
                derivative = (temp(t+δ)-temp(t-δ))/(2δ)
                a = exp(-rate*t)
                @test (ca*a+cb*(1-a))*derivative ≈ -heat*rate*a rtol = 1e-7 atol = 1e-6
                heat == 0 && @test temp(t) == T0
            end
        end
        for heat in (1e6, -1e5)
            equal = adiabatic_reaction_temperature(5.0; T0, rate, heat, ca, cb = ca)
            nearly = adiabatic_reaction_temperature(5.0; T0, rate, heat, ca, cb = ca+1e-8)
            @test nearly ≈ equal rtol = 1e-10
        end
    end

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
        # 1. nth-order: dξ/dt = −k ξ₀ (ξ/ξ₀)ⁿ verified by central finite
        #    difference, before the burnout of the orders below one.
        for n in (0.0, 2.0, 2 / 3)
            ξ0, k, t = 100.0, 0.5, 1.0
            δ = 1e-6
            dξdt =
                (
                    nth_order_xi(t + δ; ξ0 = ξ0, k = k, n = n) -
                    nth_order_xi(t - δ; ξ0 = ξ0, k = k, n = n)
                ) / (2δ)
            ξt = nth_order_xi(t; ξ0 = ξ0, k = k, n = n)
            @test dξdt ≈ -k * ξ0 * (ξt / ξ0)^n rtol = 1e-6
        end

        # 2. Burnout: exact zero at t_b for n < 1, and m^(1/3) linear in t for
        #    the contracting-sphere order n = 2/3 (FDS NaHCO₃ marker).
        tb = nth_order_burnout_time(k = 0.05, n = 0.0)
        @test tb == 20.0
        @test nth_order_xi(tb; ξ0 = 100.0, k = 0.05, n = 0.0) == 0.0
        tb3 = nth_order_burnout_time(k = 0.3, n = 2 / 3)
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

        # 4. Series scheme on conversions: the Bateman form at second-step order
        #    one, and the balance at another order by central finite difference.
        k1, k2 = 0.05, 0.1
        for t in (5.0, 20.0, 80.0)
            ξA, ξB, _ = series_chain_xi(t; ξ0 = 1.0, k1 = k1, k2 = k2)
            eA, eB = series_chain_conversion(t; k1 = k1, k2 = k2, n2 = 1.0)
            @test eA ≈ ξA rtol = 1e-13
            @test eB ≈ ξB rtol = 1e-13
            δ = 1e-3
            b₊ = series_chain_conversion(t + δ; k1 = k1, k2 = k2, n2 = 1.8)[2]
            b₋ = series_chain_conversion(t - δ; k1 = k1, k2 = k2, n2 = 1.8)[2]
            a, b = series_chain_conversion(t; k1 = k1, k2 = k2, n2 = 1.8)
            @test (b₊ - b₋) / 2δ ≈ k1 * a - k2 * b^1.8 rtol = 1e-6
        end
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

    @testset "Landau ablation traveling wave (QAM 8:81)" begin
        q, ρ, c, k = 2.8e6, 1000.0, 1000.0, 1.0
        α = k / (ρ * c)
        T0, Tp, ΔH = 300.0, 700.0, 1.0e6
        v = landau_vss(q = q, ρ = ρ, c = c, ΔH = ΔH, T0 = T0, Tp = Tp)

        # 1. Wave-frame ODE residual α T″ + v T′ = 0 by central differences.
        for ζ in (1e-4, 5e-4, 2e-3)
            dζ = 1e-7
            Tf(z) = landau_T(z; T0 = T0, Tp = Tp, v = v, α = α)
            T2 = (Tf(ζ + dζ) - 2 * Tf(ζ) + Tf(ζ - dζ)) / dζ^2
            T1 = (Tf(ζ + dζ) - Tf(ζ - dζ)) / (2 * dζ)
            scale = v^2 / α^2 * (Tp - T0) # |T″| scale at ζ = 0
            @test abs(α * T2 + v * T1) < 1e-5 * α * scale
        end

        # 2. Surface energy balance: q = k|dT/dζ|₀ + ρ v ΔH, with the gradient
        #    from a second-order one-sided difference (independent of the v_ss
        #    algebra; a first-order stencil has a (v/2α)·dζ ≈ 1e-5 rel. bias).
        dζ = 1e-8
        Tl(z) = landau_T(z; T0 = T0, Tp = Tp, v = v, α = α)
        grad0 = (-3 * Tl(0.0) + 4 * Tl(dζ) - Tl(2 * dζ)) / (2 * dζ)
        @test k * abs(grad0) + ρ * v * ΔH ≈ q rtol = 1e-6

        # 3. Boundary values and monotone decay to T₀.
        @test landau_T(0.0; T0 = T0, Tp = Tp, v = v, α = α) == Tp
        @test landau_T(20 * α / v; T0 = T0, Tp = Tp, v = v, α = α) ≈ T0 atol =
            1e-5 * (Tp - T0)
    end

    @testset "Neumann two-phase Stefan (C&J §11.2 / Gpyro §3.4.3)" begin
        p = stefan_gpyro_problem()
        λ = stefan_lambda(p)
        αs, αl = α_solid(p), α_liquid(p)

        # 1. Published eigenvalue: dissertation reports λ = 0.2317 (4 digits).
        @test abs(λ - 0.2317) < 5e-5
        # 2. Root satisfies the eigencondition to bisection precision.
        @test abs(stefan_eigencondition(λ, p)) < 1e-10

        # 3. Interface: both branches equal T_m at x = s(t).
        t = 150.0
        s = stefan_front(t; λ = λ, αs = αs)
        @test stefan_T(s - 1e-12, t, p, λ) ≈ p.T_m rtol = 1e-10
        @test stefan_T(s + 1e-12, t, p, λ) ≈ p.T_m rtol = 1e-10

        # 4. Stefan flux-jump condition k_s∂T_s/∂x − k_l∂T_l/∂x = ρΔH_m ds/dt
        #    via one-sided finite differences on each branch.
        dx = 1e-7
        g_s = (stefan_T(s, t, p, λ) - stefan_T(s - dx, t, p, λ)) / dx
        g_l = (stefan_T(s + dx, t, p, λ) - stefan_T(s, t, p, λ)) / dx
        dsdt = λ * sqrt(αs / t)
        @test p.k_s * g_s - p.k_l * g_l ≈ p.ρ * p.ΔH_m * dsdt rtol = 1e-5

        # 5. Both PDEs by finite differences (solid and liquid interior points).
        for (x0, αx) in ((0.4 * s, αs), (2.5 * s, αl))
            dxf, dtf = 1e-5, 1e-3
            Tt = (stefan_T(x0, t + dtf, p, λ) - stefan_T(x0, t - dtf, p, λ)) / (2 * dtf)
            Txx =
                (
                    stefan_T(x0 + dxf, t, p, λ) - 2 * stefan_T(x0, t, p, λ) +
                    stefan_T(x0 - dxf, t, p, λ)
                ) / dxf^2
            @test Tt ≈ αx * Txx rtol = 1e-4
        end

        # 6. Global energy balance by quadrature (independent check of the
        #    whole family, λ included): heat extracted at the wall
        #    Q_w(t) = 2 k_s (T_m − T∞) √t / (erf λ √(π αs))
        #    must equal latent release ρΔH_m·s(t) plus the sensible change
        #    ∫ ρ ĉ(x)(T₀ − T(x,t)) dx, where solidified material contributes
        #    c_l (T₀ → T_m) then c_s (T_m → T). Run on an asymmetric parameter
        #    set too, which exercises every term of the eigencondition.
        for p2 in (
            p,
            StefanProblem(
                k_s = 0.45,
                k_l = 0.20,
                ρ = 900.0,
                c_s = 2000.0,
                c_l = 1500.0,
                T_m = 433.0,
                ΔH_m = 80.0e3,
                T0 = 445.0,
                T∞ = 415.0,
            ),
        )
            λ2 = stefan_lambda(p2)
            αs2, αl2 = α_solid(p2), α_liquid(p2)
            t2 = 200.0
            s2 = stefan_front(t2; λ = λ2, αs = αs2)
            Q_wall = 2 * p2.k_s * (p2.T_m - p2.T∞) * sqrt(t2) / (erf(λ2) * sqrt(π * αs2))
            # sensible deficit by composite trapezoid, split at the front
            function sens_deficit(a, b, n, phase)
                xs = range(a, b, length = n + 1)
                function f(x)
                    T = stefan_T(x, t2, p2, λ2)
                    if phase === :solid
                        return p2.c_l * (p2.T0 - p2.T_m) + p2.c_s * (p2.T_m - T)
                    else
                        return p2.c_l * (p2.T0 - T)
                    end
                end
                vals = [f(x) for x in xs]
                return p2.ρ * (sum(vals) - 0.5 * (vals[1] + vals[end])) * (b - a) / n
            end
            x_far = s2 + 14 * sqrt(αl2 * t2)  # erfc tail < 1e-20 beyond
            E_sens =
                sens_deficit(0.0, s2, 40_000, :solid) +
                sens_deficit(s2, x_far, 200_000, :liquid)
            E_latent = p2.ρ * p2.ΔH_m * s2
            @test Q_wall ≈ E_latent + E_sens rtol = 1e-8
        end
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

# MMS driver self-tests (formula-provenance rule applied to the source
# generator itself): every ForwardDiff-assembled continuous-operator term is
# checked against an independently hand-derived closed form, then the driver
# is exercised end-to-end on a heat-only manufactured solution (spatial and
# temporal order).
@testset "MMS driver (exact/mms.jl)" begin
    @testset "conduction + storage source vs hand derivation" begin
        ρs, cs, ks = 1000.0, 1200.0, 0.4
        L, τ = 0.02, 5.0
        solid = SolidComponent(:s, ρ = ρs, c = cs, k = ks)
        mat = Material(name = :MMSSelfHC, components = (solid,), reactions = ())
        fields = MMSFields(
            (z, t) -> 300.0 + 50.0 * sin(π * z / L) * exp(-t / τ) + 20.0 * z / L,
            ((z, t) -> ρs,),
        )
        for z in (0.0, 0.004, 0.013, L), t in (0.0, 1.0, 4.0)
            dTdt = -(50.0 / τ) * sin(π * z / L) * exp(-t / τ)
            dTdz = 50.0 * (π / L) * cos(π * z / L) * exp(-t / τ) + 20.0 / L
            d2Tdz2 = -50.0 * (π / L)^2 * sin(π * z / L) * exp(-t / τ)
            @test mms_heat_flux(mat, fields, z, t) ≈ -ks * dTdz rtol = 1e-12
            S_hand = dTdt - ks * d2Tdz2 / (ρs * cs)
            @test mms_energy_source(mat, fields, z, t) ≈ S_hand rtol = 1e-12
        end
    end

    @testset "kinetics sources vs hand derivation (A→B, E = 0)" begin
        L = 0.02
        kA = 0.5          # rate constant [1/s]: A = kA, E = 0
        h_r = 1.0e5       # endothermic
        compA = SolidComponent(:A, ρ = 800.0, c = 1000.0, k = 0.3)
        compB = SolidComponent(:B, ρ = 600.0, c = 900.0, k = 0.2)
        rxn = Reaction(:r, 1 => 2; A = kA, E = 0.0, h = h_r)
        mat = Material(
            name = :MMSSelfRX,
            components = (compA, compB),
            reactions = (rxn,),
            depletion_limiter = DepletionLimiter(enabled = false),
        )
        τ = 3.0
        fields = MMSFields(
            (z, t) -> 500.0,
            ((z, t) -> 400.0 + 100.0 * (z / L) * exp(-t / τ), (z, t) -> 150.0),
        )
        for z in (0.0, 0.007, L), t in (0.0, 2.0)
            ξA = 400.0 + 100.0 * (z / L) * exp(-t / τ)
            dξA_dt = -(100.0 / τ) * (z / L) * exp(-t / τ)
            rate = kA * ξA
            @test mms_species_source(mat, fields, 1, z, t) ≈ dξA_dt + rate rtol = 1e-12
            @test mms_species_source(mat, fields, 2, z, t) ≈ -rate rtol = 1e-12
            # T* constant ⇒ q* ≡ 0; S_T = −Q_rxn/ρcp = +h·r/ρcp (h > 0 endo).
            ρcp = ξA * 1000.0 + 150.0 * 900.0
            @test mms_energy_source(mat, fields, z, t) ≈ h_r * rate / ρcp rtol = 1e-12
        end
    end

    @testset "gas flux, S_adv and species source (single gas)" begin
        L = 0.02
        λ = 2.0e-5
        cp_g = 1100.0
        solid = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3, λ = λ)
        gas = GaseousComponent(:g, M = 0.028, c = cp_g, k = 0.03)
        mat = Material(name = :MMSSelfGT, components = (solid, gas), reactions = ())
        # One condensed matrix gives Λ=λ; φ=1−700/800=1/8.
        τ = 4.0
        T_f(z, t) = 400.0 + 50.0 * z / L
        ξg_f(z, t) = 2.0 + sin(π * z / L) * exp(-t / τ)
        fields = MMSFields(T_f, ((z, t) -> 700.0, ξg_f))
        J_hand(z, t) = begin
            T = T_f(z, t)
            dTdz = 50.0 / L
            dξdz = (π / L) * cos(π * z / L) * exp(-t / τ)
            -(λ / (T*(1-700/800))) * (dξdz * T + ξg_f(z, t) * dTdz)
        end
        for z in (0.0, 0.006, 0.015, L), t in (0.0, 1.5)
            @test mms_gas_flux(mat, fields, 2, z, t) ≈ J_hand(z, t) rtol = 1e-12
            @test mms_S_adv(mat, fields, z, t) ≈ -cp_g * J_hand(z, t) * (50.0 / L) rtol =
                1e-12
            # Species source: ∂tξ + ∂zJ, with ∂zJ from a central difference of
            # the hand flux (independent of the ForwardDiff evaluation).
            δ = 1e-7 * L
            dJdz = (J_hand(z + δ, t) - J_hand(z - δ, t)) / (2δ)
            dξdt = -sin(π * z / L) * exp(-t / τ) / τ
            @test mms_species_source(mat, fields, 2, z, t) ≈ dξdt + dJdz rtol = 1e-6
        end
    end

    @testset "end-to-end heat-only MMS (spatial and temporal order)" begin
        ρs, cs, ks = 1000.0, 1200.0, 0.4
        L, τ = 0.02, 5.0
        fields = MMSFields(
            (z, t) -> 300.0 + 50.0 * sin(π * z / L) * exp(-t / τ) + 20.0 * z / L,
            ((z, t) -> ρs,),
        )
        function build(n)
            solid = SolidComponent(:s, ρ = ρs, c = cs, k = ks)
            mat = Material(name = :MMSSelfE2E, components = (solid,), reactions = ())
            mesh = uniform_mesh(L, n)
            bcs = mms_thermal_bcs(mat, fields, L)
            bc_set = BoundaryConditionSet(thermal = (top = bcs.top, bottom = bcs.bottom))
            return PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> fields.T(z, 0.0),
                ξ_initial = [z -> ρs],
                tspan = (0.0, 2.0),
                output_times = [0.0, 2.0],
            )
        end

        # Spatial order ≈ 2 (L2 of T at t_end; tight tolerances keep temporal
        # error subdominant).
        errs = [
            mms_L2_T(solve_mms(build(n), fields; abstol = 1e-11, reltol = 1e-9), fields) for n in (16, 32)
        ]
        order = log(errs[1] / errs[2]) / log(2.0)
        @test errs[2] < errs[1]
        @test order > 1.7

        # Temporal order of the default integrator (Rodas5P, observed 4.6;
        # finer steps reach roundoff): fixed mesh, fixed-dt sweep against a
        # small-dt reference on the SAME mesh so the spatial error cancels.
        ref = solve_mms(build(16), fields; dt = 0.0125)
        e_dt = [
            maximum(abs.(solve_mms(build(16), fields; dt = dt).T[:, end] - ref.T[:, end])) for dt in (0.2, 0.1)
        ]
        p_t = log(e_dt[1] / e_dt[2]) / log(2.0)
        @test p_t > 3.0
    end

    @testset "induced-motion machinery (θ*, w*, material-frame source)" begin
        # V-MMS-04 companion checks (formula-provenance rule): the package-
        # kinetics θ*, the Gauss–Legendre w* = ∫θ*, and the material-frame
        # condensed source are each verified against an independently
        # hand-derived closed form. Design mirrors the case: constant ρ_s,
        # the gas outside the strain, single first-order s → g, depletion limiter
        # numerically inactive (tanh(ξ_s/1.0) ≡ 1.0 in Float64 for
        # ξ_s ≳ 20 kg/m³ — here ξ_s ≥ 500).
        L = 0.01
        ρs = 800.0
        A_r, E_r = 5.0e4, 80.0e3
        solid = SolidComponent(:s, ρ = ρs, c = 1200.0, k = 0.2, λ = 2.0e-5)
        gas = GaseousComponent(:g, M = 0.028, c = 1000.0, k = 0.03)
        rxn = Reaction(:r, 1 => 2; A = A_r, E = E_r, h = 5.0e4)
        mat = Material(name = :MMSSelfIM, components = (solid, gas), reactions = (rxn,))

        τa, τb = 25.0, 30.0
        T_f(z, t) = 600.0 + 15.0 * sin(2π * t / 40.0)
        a_f(t) = 600.0 + 20.0 * cos(2π * t / τa)
        b_f(t) = 60.0 * exp(-t / τb)
        ξs_f(z, t) = a_f(t) + b_f(t) * sin(π * z / L)
        ξg_f(z, t) = 2.0 + 0.8 * sin(π * z / L) * exp(-t / 15.0)
        fields = MMSFields(T_f, (ξs_f, ξg_f))

        R = Pyrolysis.Materials.R_GAS
        k_hand(t) = A_r * exp(-E_r / (R * T_f(0.0, t)))
        rate_hand(z, t) = k_hand(t) * ξs_f(z, t)          # limiter ≡ 1 here
        θ_hand(z, t) = -rate_hand(z, t) / ρs               # the gas never enters θ
        # w* = ∫₀ᶻ θ* dζ, closed form for the sin profile:
        w_hand(z, t) =
            -(k_hand(t) / ρs) * (a_f(t) * z + b_f(t) * (L / π) * (1.0 - cos(π * z / L)))

        motion = MMSMotion((z, t) -> mms_mesh_velocity(mat, fields, z, t))
        for z in (0.0, 0.003, 0.007, L), t in (0.0, 6.0, 17.0)
            @test mms_theta(mat, fields, z, t) ≈ θ_hand(z, t) rtol = 1e-12
            # 32-node Gauss–Legendre vs the closed-form integral: machine
            # precision for this smooth integrand (requirement 4).
            if z > 0.0
                @test mms_mesh_velocity(mat, fields, z, t) ≈ w_hand(z, t) rtol = 1e-12
            else
                @test mms_mesh_velocity(mat, fields, z, t) == 0.0
            end

            # Material-frame solid source: S_s = ∂tξ* + w*·∂zξ* − ω̇* + ξ*·θ*
            # with ω̇_s = −rate.
            dξs_dt =
                -20.0 * (2π / τa) * sin(2π * t / τa) -
                (60.0 / τb) * exp(-t / τb) * sin(π * z / L)
            dξs_dz = b_f(t) * (π / L) * cos(π * z / L)
            S_hand =
                dξs_dt + w_hand(z, t) * dξs_dz + rate_hand(z, t) + ξs_f(z, t) * θ_hand(z, t)
            @test mms_species_source(mat, fields, 1, z, t; w = motion.w) ≈ S_hand rtol =
                1e-12
        end
    end

    @testset "Frank-Kamenetskii slab (FK limit + full-Arrhenius shooting)" begin
        # 1. Published criticality constants: δ_c = 0.878, θ_m,c ≈ 1.19.
        @test abs(fk_delta_crit() - 0.878) < 1e-3
        @test abs(fk_theta_max_crit() - 1.19) < 5e-3

        # 2. The FK profile satisfies θ″ + δe^θ = 0 (finite differences) and
        #    its wall BC θ(±1) = 0 to root precision.
        δ = 0.44
        θm = fk_theta_max(δ)
        prof(x̃) = fk_theta_profile(x̃; θm = θm, δ = δ)
        for x̃ in (0.0, 0.35, 0.8)
            dx = 1e-5
            θxx = (prof(x̃ + dx) - 2 * prof(x̃) + prof(x̃ - dx)) / dx^2
            @test θxx ≈ -δ * exp(prof(x̃)) rtol = 1e-5
        end
        @test abs(prof(1.0)) < 1e-12
        @test abs(prof(-1.0)) < 1e-12

        # 3. Small-δ limit: θ_m → δ/2 (source ≈ constant, parabolic profile).
        @test fk_theta_max(1e-4) ≈ 5e-5 rtol = 1e-3

        # 4. Supercritical δ must throw (no steady root).
        @test_throws ArgumentError fk_theta_max(0.9)

        # 5. Full-Arrhenius shooting: wall BC hit to bisection precision, and
        #    in the deep FK limit (E/RT₀ = 400) the center rise matches the
        #    FK closed form to O(RT₀/E).
        R = 8.314462618
        T0, k, L = 500.0, 0.2, 0.01
        for (E_over_RT0, rtol_fk) in ((400.0, 3e-3), (40.0, 5e-2))
            E = E_over_RT0 * R * T0
            εT = R * T0^2 / E
            S0 = δ * k * εT / L^2
            sh = fk_steady_shoot(S0 = S0, E = E, T0 = T0, k = k, L = L)
            @test abs(sh.T_at(L) - T0) < 1e-8 * T0
            θ_center = (sh.T_center - T0) / εT
            @test θ_center ≈ θm rtol = rtol_fk
            # symmetric interior point against the FK profile shape
            θ_half = (sh.T_at(0.5 * L) - T0) / εT
            @test θ_half ≈ prof(0.5) rtol = 2 * rtol_fk
        end
    end

    @testset "constant-heating-rate TGA reference" begin
        ξ0, A, E, β, T0 = 1000.0, 1.0e10, 130.0e3, 0.5, 300.0
        R = 8.314462618

        # 1. ODE residual by central finite difference: dξ/dt = −A e^(−E/RT) ξ.
        for t in (500.0, 600.0, 700.0)
            δt = 1e-3
            ξp = tga_xi(t + δt; ξ0 = ξ0, A = A, E = E, β = β, T0 = T0)
            ξm = tga_xi(t - δt; ξ0 = ξ0, A = A, E = E, β = β, T0 = T0)
            ξt = tga_xi(t; ξ0 = ξ0, A = A, E = E, β = β, T0 = T0)
            T = T0 + β * t
            @test (ξp - ξm) / (2 * δt) ≈ -A * exp(-E / (R * T)) * ξt rtol = 1e-6
        end

        # 2. Temperature integral vs the standard asymptotic expansion
        #    I ≈ (RT²/E)e^(−E/RT)(1 − 2u + 6u²), u = RT/E; next term is
        #    24u³ ≈ 4e-4 at u = 0.037 (T = 580 K).
        for T in (500.0, 580.0, 650.0)
            u = R * T / E
            asym = (R * T^2 / E) * exp(-E / (R * T)) * (1 - 2u + 6u^2)
            @test tga_temperature_integral(T; T0 = T0, E = E) ≈ asym rtol = 30 * u^3
        end

        # 3. Richardson: doubling the Simpson resolution changes nothing at
        #    the tolerance the cases rely on.
        I1 = tga_temperature_integral(650.0; T0 = T0, E = E, n = 20_000)
        I2 = tga_temperature_integral(650.0; T0 = T0, E = E, n = 40_000)
        @test I1 ≈ I2 rtol = 1e-10

        # 4. Monotone decay from ξ0 to burnout.
        @test tga_xi(0.0; ξ0 = ξ0, A = A, E = E, β = β, T0 = T0) == ξ0
        ξs = [tga_xi(t; ξ0 = ξ0, A = A, E = E, β = β, T0 = T0) for t = 0.0:50.0:800.0]
        @test issorted(ξs, rev = true)
        @test ξs[end] / ξ0 < 1e-40
    end
end

@testset "Gray-slab and two-stream radiation references" begin
    # A column of uneven cells: optical thicknesses around τ0, emissive powers
    # around 10 kW/m² (σT⁴ at about 650 K).
    column(N; τ0 = 0.3) = (
        [τ0 * (1 + 0.5 * sin(1.3i)) for i = 1:N],
        [1e4 * (1 + 0.3 * cos(0.7i)) for i = 1:N],
    )
    θ1, θ2 = deg2rad(20.0), deg2rad(73.0)

    @testset "gray slab: black faces against the jump sum" begin
        # The quadrature solution against the E₃ sum, which needs no angular
        # quadrature; τ0 = 3 and 60 cells reach the Rosseland regime.
        for (N, τ0) in ((1, 1e-3), (7, 0.05), (30, 0.4), (60, 3.0))
            τ, B = column(N; τ0)
            r = gray_slab_exact(τ, B; q = 2.5e4, B_amb = 450.0, B_b = 9e3)
            exact = gray_slab_black_faces(τ, B; q = 2.5e4, B_amb = 450.0, B_b = 9e3)
            @test maximum(abs, r.q .- exact) < 1e-12 * 2.5e4
        end
        # Limits of the closed form: a cold slab transmits 2qE₃(τ_L), an
        # isothermal one emits B(1 − 2E₃(τ_L)), with E₃(0) = 1/2.
        @test expint(3, 0.0) == 0.5
        for τL in (1e-4, 0.3, 10.0)
            τ = fill(τL / 5, 5)
            @test gray_slab_black_faces(τ, zeros(5); q = 1.0)[end] ≈ 2expint(3, τL) rtol =
                1e-14
            @test -gray_slab_black_faces(τ, fill(1e4, 5))[1] ≈ 1e4 * (1 - 2expint(3, τL)) rtol =
                1e-12
        end
    end

    @testset "gray slab: isothermal enclosures have zero source" begin
        faces = (
            (0.9, :specular, 1.0, :specular),
            (0.7, :diffuse, 0.4, :diffuse),
            (0.6, :specular, 0.0, :specular),
            (0.8, :diffuse, 0.0, :diffuse),
            (0.5, :specular, 0.3, :diffuse),
        )
        for τ0 in (1e-3, 0.3, 3.0), (ε, top, ε_b, back) in faces
            τ, _ = column(20; τ0)
            Bc = 1e4
            r = gray_slab_exact(
                τ,
                fill(Bc, 20);
                ε,
                top_reflection = top,
                B_amb = Bc,
                ε_b,
                B_b = Bc,
                back_reflection = back,
            )
            @test maximum(abs, r.S) < 1e-12 * Bc
            @test abs(r.q_out - r.q_enter) < 1e-12 * Bc
        end
    end

    @testset "gray slab: beams" begin
        τ, B = column(25; τ0 = 0.1)
        Z = zeros(25)
        depth = [0.0; cumsum(τ)]
        # Collimated and cone transmission through a cold, black-faced slab.
        r = gray_slab_exact(τ, Z; q = 1.0, heater = :collimated)
        @test maximum(abs, r.q .- exp.(-depth)) < 1e-15
        r = gray_slab_exact(τ, Z; q = 1.0, heater = (θ1, θ2))
        @test maximum(abs, r.q .- gray_slab_cone_transmission.(1.0, depth, θ1, θ2)) < 1e-12
        @test gray_slab_cone_transmission(1.0, 0.3, 0.0, π / 2) ≈ 2expint(3, 0.3) rtol =
            1e-14
        # A collimated beam meeting reflecting faces against a cone narrowed
        # to 0.1° about the normal, which reaches the radiosities through
        # the angular quadrature instead (4.4e-7 apart when pinned).
        for (ε, top, ε_b, back) in (
            (0.8, :specular, 0.0, :diffuse),
            (0.7, :diffuse, 0.5, :diffuse),
            (0.9, :diffuse, 0.0, :specular),
            (0.6, :specular, 0.0, :specular),
        )
            kw = (; ε, top_reflection = top, ε_b, back_reflection = back, q = 1.0)
            a = gray_slab_exact(τ, Z; kw..., heater = :collimated)
            b = gray_slab_exact(τ, Z; kw..., heater = (0.0, deg2rad(0.1)))
            @test sum(abs, a.S .- b.S) < 2e-6
        end
        # Conservation with every input at once.
        r = gray_slab_exact(
            τ,
            B;
            ε = 0.7,
            top_reflection = :diffuse,
            q = 2e4,
            heater = (θ1, θ2),
            B_amb = 500.0,
            ε_b = 0.5,
            B_b = 3e3,
        )
        @test abs(sum(r.S) - (r.q_enter - r.q_out - r.q_back)) < 1e-11 * 2e4
    end

    # Two-stream configurations: thick, thin and mixed cells; every back;
    # diffuse input on both faces; collimated and cone beams, reflected at the
    # back into q⁺ or specularly; the medium starting at the first cell or
    # deeper.
    τ_mixed = [1e-6, 0.01, 0.5, 5.0, 0.2, 50.0, 1e-3, 0.8, 2.0, 0.05, 0.3, 1.0]
    B_mixed = [1e4 * (1 + 0.3 * cos(0.7i)) for i = 1:12]
    cone = two_stream_cone_beams(3e4, θ1, θ2; a = 0.9)
    configurations = (
        (; ε = 1.0, ε_b = 1.0, Q_b = 2e3, beams = ()),
        (; ε = 0.9, ε_b = 0.0, beams = [TwoStreamBeam(1.0, 2e4)]),
        (;
            ε = 0.9,
            ε_b = 0.6,
            Q_b = 2e3,
            beams = [TwoStreamBeam(1.0, 2e4), TwoStreamBeam(0.6, 9e3)],
        ),
        (; ε = 0.3, ε_b = 0.6, Q_b = 2e3, beams = cone),
        (; ε = 0.96, ε_b = 0.0, beams = cone, beam_back = :specular),
        (; ε = 0.5, ε_b = 0.0, beams = [TwoStreamBeam(0.8, 1e4)], beam_back = :specular),
    )
    scale = 3e4

    @testset "two-stream: the rows hold, against a dense solve" begin
        for c in configurations, first in (1, 4)
            kw = (; c..., Q = 1.5e4, first)
            r = two_stream_solve(τ_mixed, B_mixed; kw...)
            res = two_stream_residual(r.qm, r.qp, r.Id, r.Iu, τ_mixed, B_mixed; kw...)
            @test maximum(abs, res.thermal) < 1e-12 * scale
            @test maximum(abs, res.down; init = 0.0) < 1e-12 * scale
            @test maximum(abs, res.up; init = 0.0) < 1e-12 * scale
            # The rows are affine in every flux: assembled as a matrix
            # column by column and solved by LU, independently of the sweep.
            N1, nb = length(τ_mixed) + 1, length(c.beams)
            n = 2N1 + 2N1 * nb
            unpack(u) = (
                u[1:N1],
                u[(N1+1):(2N1)],
                reshape(u[(2N1+1):(2N1+N1*nb)], N1, nb),
                reshape(u[(2N1+N1*nb+1):end], N1, nb),
            )
            rows(u) =
                let x = two_stream_residual(unpack(u)..., τ_mixed, B_mixed; kw...)
                    [vec(x.thermal); vec(x.down); vec(x.up)]
                end
            r0 = rows(zeros(n))
            M = reduce(hcat, (rows([k == j ? 1.0 : 0.0 for k = 1:n]) .- r0 for j = 1:n))
            dense = M \ (-r0)
            @test maximum(abs, dense .- [r.qm; r.qp; vec(r.Id); vec(r.Iu)]) < 1e-11 * scale
            # Energy: the sources add to what enters less what leaves.
            @test abs(sum(r.S) - (r.q_enter - r.q_out - r.q_back)) < 1e-11 * scale
            # Faces before the first cell of the medium hold nothing, and the
            # medium solves as if it began there.
            if first > 1
                @test all(iszero, r.qm[1:(first-1)]) && all(iszero, r.qp[1:(first-1)])
                t = two_stream_solve(
                    τ_mixed[first:end],
                    B_mixed[first:end];
                    kw...,
                    first = 1,
                )
                @test maximum(abs, r.qm[first:end] .- t.qm) < 1e-12 * scale
                @test maximum(abs, r.qp[first:end] .- t.qp) < 1e-12 * scale
            end
        end
    end

    @testset "two-stream: an open back reflects 1 − ε_b of each beam into q⁺" begin
        # Cells thin enough that the beams reach the back.
        c = configurations[3]
        cells = [1, 2, 3, 5, 7, 10]
        r = two_stream_solve(τ_mixed[cells], B_mixed[cells]; c..., Q = 1.5e4)
        arriving = sum(r.Id[end, :])
        @test arriving > 1e4
        @test r.qp[end] ≈ c.ε_b * c.Q_b + (1 - c.ε_b) * (r.qm[end] + arriving) rtol = 1e-14
        @test r.q_back ≈ c.ε_b * (r.qm[end] + arriving - c.Q_b) rtol = 1e-12
    end

    @testset "two-stream: isothermal equilibrium" begin
        Bc = 5.670374419e-8 * 650.0^4
        worst = 0.0
        for N in (1, 5, 40),
            τ0 in (1e-6, 1e-2, 1.0, 50.0),
            ε in (0.3, 0.96, 1.0),
            ε_b in (0.0, 0.6, 1.0)

            τ, _ = column(N; τ0)
            r = two_stream_solve(τ, fill(Bc, N); ε, Q = Bc, ε_b, Q_b = Bc)
            worst = max(worst, maximum(abs, r.S), abs(r.qm[1] - Bc), abs(r.qp[end] - Bc))
        end
        @test worst < 1e-14 * Bc
    end

    @testset "two-stream: closed forms" begin
        D = 1.5
        # An isothermal slab with diffuse inputs on both faces: the closed form
        # against the sweep, over uneven cells of the same total.
        for (ε, ε_b, Q, Q_b) in
            ((0.8, 0.0, 2e4, 0.0), (0.4, 0.6, 5e3, 2e4), (1.0, 1.0, 0.0, 0.0)),
            τL in (1e-3, 0.4, 6.0)

            τ = column(9; τ0 = 1.0)[1]
            τ .*= τL / sum(τ)
            Bc = 1e4
            r = two_stream_solve(τ, fill(Bc, 9); D, ε, Q, ε_b, Q_b)
            c = two_stream_isothermal_slab(; B = Bc, x = D * τL, ε, Q, ε_b, Q_b)
            @test maximum(
                abs,
                (
                    r.qm[1] - c.qm_e,
                    r.qp[1] - c.qp_e,
                    r.qm[end] - c.qm_b,
                    r.qp[end] - c.qp_b,
                ),
            ) < 1e-12 * 2e4
            @test abs(sum(r.S) - (c.intake - c.q_back)) < 1e-11 * 2e4
        end
        # The net intake of an isothermal slab before a mirror, facing Q:
        # ε(Q − B)(1 − t²)/(1 − (1 − ε)t²), t = e^{−Dτ_L}; zero at Q = B.
        c = two_stream_isothermal_slab(; B = 1e4, x = D * 0.1, ε = 0.8, Q = 2e4, ε_b = 0.0)
        t = exp(-D * 0.1)
        @test c.intake ≈ 0.8 * 1e4 * (1 - t^2) / (1 - 0.2t^2) rtol = 1e-14
        @test c.intake ≈ 2434.1 atol = 0.05
        c = two_stream_isothermal_slab(; B = 1e4, x = 0.7, ε = 0.8, Q = 1e4, ε_b = 0.0)
        @test c.intake == 0 && c.qm_e == 1e4 && c.qp_e == 1e4
        # Emission into 0 K: an open black back, a mirror, one material on
        # both faces; thin slabs keep their digits through expm1.
        for x in (1e-9, 1e-6, 1e-3, 0.3, 5.0), ε in (0.3, 0.9, 1.0)
            h, h2 = -expm1(-x), -expm1(-2x)
            t = exp(-x)
            B = 1e4
            @test two_stream_slab_emission(B, x, ε, 1.0) ≈ ε * B * h rtol = 1e-14
            @test two_stream_slab_emission(B, x, ε, 0.0) ≈ ε * B * h2 / (1 - (1 - ε) * t^2) rtol =
                1e-14
            @test two_stream_slab_emission(B, x, ε, ε) ≈ ε * B * h / (1 - (1 - ε) * t) rtol =
                1e-14
            for ε_b in (0.0, 0.5, 1.0)
                τ = fill(x / D / 4, 4)
                r = two_stream_solve(τ, fill(B, 4); D, ε, ε_b)
                @test abs(r.q_out - two_stream_slab_emission(B, x, ε, ε_b)) < 1e-12 * B
            end
        end
        # A linear emissive power, the boundary inputs set to the rows'
        # particular solution: every face carries g·coth(Dτ/2).
        for τc in (1e-3, 0.1, 1.0, 5.0)
            N, g, x = 40, 50.0, D * τc
            B = [2e4 - g * i for i = 1:N]
            Q = B[1] + g / (-expm1(-x))
            Q_b = B[N] - g / (-expm1(-x))
            r = two_stream_solve(fill(τc, N), B; D, ε = 1.0, Q, ε_b = 1.0, Q_b)
            @test maximum(abs, (r.qm .- r.qp) .- two_stream_linear_flux(g, x)) < 1e-12 * 2e4
        end
        @test two_stream_linear_flux(1.0, 0.7) ≈ coth(0.35) rtol = 1e-15
        # A beam through a slab with no emission, reflected at the back into
        # q⁺; with a mirror the absorbed flux is
        # q0[1 − εtT_c/(1 − (1 − ε)t²)], T_c = e^{−τ_L/μ}.
        for μ in (1.0, 0.6), ε in (0.3, 0.9), ε_b in (0.0, 0.6, 1.0), τL in (0.05, 1.0, 4.0)
            q0 = 2e4
            τ = column(10; τ0 = 1.0)[1]
            τ .*= τL / sum(τ)
            r = two_stream_solve(τ, zeros(10); D, ε, ε_b, beams = [TwoStreamBeam(μ, q0)])
            c = two_stream_beam_slab(q0, τL, ε, ε_b; D, μ)
            @test abs(r.q_out - c.q_out) < 1e-12 * q0
            @test abs(r.q_back - c.q_back) < 1e-12 * q0
            @test abs(sum(r.S) - c.absorbed) < 1e-12 * q0
            if ε_b == 0
                t, T_c = exp(-D * τL), exp(-τL / μ)
                @test c.absorbed ≈ q0 * (1 - ε * t * T_c / (1 - (1 - ε) * t^2)) rtol = 1e-14
            end
        end
    end

    @testset "two-stream: the sweep allocates nothing" begin
        w = TwoStreamWork(length(τ_mixed), length(cone))
        sweep(w) = @allocated two_stream_sweep!(
            w,
            τ_mixed,
            B_mixed;
            ε = 0.9,
            Q = 1.5e4,
            ε_b = 0.6,
            Q_b = 2e3,
            beams = cone,
            first = 3,
        )
        sweep(w)
        @test sweep(w) == 0
    end

    @testset "beams of both references agree" begin
        τ = column(120; τ0 = 0.15)[1]
        Z = zeros(120)
        # A collimated beam through black faces, and one reflected
        # specularly by a mirror and by the exposed face: identical.
        for (ε, ε_b, beam_back) in ((1.0, 1.0, :diffuse), (0.7, 0.0, :specular))
            a = two_stream_solve(
                τ,
                Z;
                ε,
                ε_b,
                beams = [TwoStreamBeam(1.0, ε * 1e4)],
                beam_back,
            )
            b = gray_slab_exact(τ, Z; ε, ε_b, q = 1e4, heater = :collimated)
            @test maximum(abs, a.S .- b.S) < 1e-12 * 1e4
        end
        # A cone carried by three directions against the exact cone: the
        # shares add up, and the source misplaced is a small fraction of the
        # flux (1.16e-4 of it when this tolerance was pinned; 2.1e-3 with two
        # directions, 7.3e-6 with four).
        @test sum(b.q0 for b in two_stream_cone_beams(1.0, θ1, θ2)) ≈ 1 rtol = 1e-15
        a = two_stream_solve(
            τ,
            Z;
            ε = 1.0,
            ε_b = 1.0,
            beams = two_stream_cone_beams(1e4, θ1, θ2),
        )
        b = gray_slab_exact(τ, Z; q = 1e4, heater = (θ1, θ2))
        @test sum(abs, a.S .- b.S) / 1e4 < 2e-4
    end

    # The continuous two-stream of the manufactured solutions (exact/mms.jl),
    # checked against its limits and against the discrete rows, an independent
    # code: the staircase of a smooth emissive power converges to it at second
    # order.
    @testset "continuous two-stream of the manufactured solutions" begin
        σ = Pyrolysis.Materials.STEFAN_BOLTZMANN
        D, L, κ, ε = 1.5, 0.01, 200.0, 0.8
        β = D * κ
        xs = range(0.0, L, length = 11)
        backs = ((0.0, 0.0), (0.6, 1.0))   # (ε_b, Q_b / B): a mirror, an open back

        # A medium facing its own temperature on both sides: q± = B, no source.
        B0 = σ * 700.0^4
        for (ε_b, r_b) in backs, x in xs
            qm, qp = mms_two_stream(_ -> B0, x, L, κ; ε, Q = B0, ε_b, Q_b = r_b * B0)
            @test abs(qm - B0) < 1e-12 * B0
            @test abs(qp - B0) < 1e-12 * B0
            @test abs(
                mms_radiative_source(_ -> B0, x, L, κ; ε, Q = B0, ε_b, Q_b = r_b * B0),
            ) < 1e-12 * β * B0
        end

        # B linear in depth with the face inputs of the particular solution
        # q∓ = B ∓ g/β: no source, and the net flux −2g/β, Rosseland's at D = 3/2.
        g, ε_b = 5e5, 0.6
        B(x) = 1e4 + g * x
        Q = (ε * B(0.0) - (2 - ε) * g / β) / ε
        Q_b = (B(L) + g / β - (1 - ε_b) * (B(L) - g / β)) / ε_b
        for x in xs
            qm, qp = mms_two_stream(B, x, L, κ; ε, Q, ε_b, Q_b)
            @test abs(qm - (B(x) - g / β)) < 1e-12 * B(L)
            @test abs(qp - (B(x) + g / β)) < 1e-12 * B(L)
            @test abs(mms_radiative_source(B, x, L, κ; ε, Q, ε_b, Q_b)) < 1e-12 * β * B(L)
        end
        @test -2g / β ≈ -(4 / (3κ)) * g rtol = 1e-15

        # A smooth emissive power: the discrete rows on staircases of N cells
        # (cell sources per unit volume at the cell centres) converge to the
        # continuous source at second order, for either back.
        Bs(x) = σ * (600.0 + 150.0 * sin(π * x / L) + 40.0 * x / L)^4
        Q = σ * 800.0^4
        for (ε_b, r_b) in backs
            errors = map((500, 1000, 2000)) do N
                Δ = L / N
                xc = [(i - 0.5) * Δ for i = 1:N]
                discrete = two_stream_solve(
                    fill(κ * Δ, N),
                    Bs.(xc);
                    D,
                    ε,
                    Q,
                    ε_b,
                    Q_b = r_b * Q,
                )
                exact = [
                    mms_radiative_source(Bs, x, L, κ; ε, Q, ε_b, Q_b = r_b * Q) for x in xc
                ]
                sqrt(sum(abs2, discrete.S ./ Δ .- exact) / N)
            end
            rates = log2.(errors[1:2] ./ errors[2:3])
            @test all(r -> abs(r - 2) < 0.01, rates)
            @test errors[3] < 1e-6 * β * Bs(L / 2)
        end

        # Positions are states on a moving mesh: dual numbers pass through,
        # and the derivative is the source's.
        S(x) = mms_radiative_source(Bs, x, L, κ; ε, Q, ε_b = 0.0)
        for x in (0.001, 0.005, 0.009)
            h = 1e-7
            @test ForwardDiff.derivative(S, x) ≈ (S(x + h) - S(x - h)) / 2h rtol = 1e-6
        end
    end
end
