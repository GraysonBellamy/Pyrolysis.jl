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

    @testset "in-depth absorption eigenseries (Lautenberger App. A)" begin
        δ, κ, k, h, q0, T0 = 0.01, 1000.0, 0.2, 20.0, 25.0e3, 300.0
        ρc = 1000.0 * 1400.0
        α = k / ρc
        ex = make_indepth_absorption_exact(
            δ = δ,
            κ = κ,
            k = k,
            h = h,
            α = α,
            q0 = q0,
            T0 = T0,
        )

        # 1. Initial condition reconstructed: series must cancel θ_ss at t=0⁺
        #    (60-term truncation; amplitude θ_ss(0) = q₀/h = 1250 K).
        for x in (0.0, 0.002, 0.005, 0.009)
            @test abs(ex.T_exact(x, 1e-12) - T0) < 1e-3 * (q0 / h)
        end

        # 2. Steady limit: series decayed ⇒ θ_ss; and steady front loss equals
        #    the deposited power q₀(1−e^(−κδ)) (identity of the derivation).
        #    Slowest mode decays like exp(−u₁²·αt/δ²) with u₁ ≈ 0.86 (Bi = 1),
        #    so αt/δ² = 60 gives e^(−44)·O(10³ K) ≪ machine precision.
        t_big = 60 * δ^2 / α
        @test ex.T_exact(0.004, t_big) ≈ T0 + ex.θ_ss(0.004) rtol = 1e-12
        @test h * ex.θ_ss(0.0) ≈ q0 * (1 - exp(-κ * δ)) rtol = 1e-12

        # 3. PDE residual by finite differences at an interior point and
        #    mid-transient time: θ_t − αθ_xx − q₀κe^(−κx)/ρc = 0.
        #    dx must be small: the exponential source gives θ_xxxx ~ q₀κ³/k
        #    ≈ 1.25e14 K/m⁴, so the central-difference truncation α·θxxxx·dx²/12
        #    only drops below 1e-4·source for dx ≲ 2e-6.
        x0, t0 = 0.004, 100.0
        dx, dt = 1e-6, 1e-2
        θt = (ex.T_exact(x0, t0 + dt) - ex.T_exact(x0, t0 - dt)) / (2dt)
        θxx =
            (ex.T_exact(x0 + dx, t0) - 2 * ex.T_exact(x0, t0) + ex.T_exact(x0 - dx, t0)) /
            dx^2
        source = q0 * κ * exp(-κ * x0) / ρc
        residual = θt - α * θxx - source
        @test abs(residual) < 3e-4 * source

        # 4. Front BC satisfied at finite t: k·θ_x(0) = h·θ(0). Second-order
        #    one-sided stencil (plain forward difference carries a θxx·dx/2
        #    bias ≈ 0.4% from the same strong surface curvature).
        t1 = 90.0
        dxb = 1e-6
        θx0 =
            (-3 * ex.T_exact(0.0, t1) + 4 * ex.T_exact(dxb, t1) - ex.T_exact(2 * dxb, t1)) /
            (2 * dxb)
        @test k * θx0 ≈ h * (ex.T_exact(0.0, t1) - T0) rtol = 1e-3
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
# is exercised end-to-end on a heat-only manufactured solution (spatial order
# ≈ 2, temporal order ≈ 4, analytic-vs-FD Jacobian agreement).
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

    @testset "gas flux, S_conv and species source (single gas)" begin
        L = 0.02
        λ = 2.0e-5
        cp_g = 1100.0
        solid = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3, λ = λ)
        gas = GaseousComponent(:g, M = 0.028, c = cp_g, k = 0.03, λ = λ)
        mat = Material(name = :MMSSelfGT, components = (solid, gas), reactions = ())
        # Equal component λ ⇒ λ_eff = λ exactly under any volume weighting, so
        # the hand flux needs no mixture algebra.
        τ = 4.0
        T_f(z, t) = 400.0 + 50.0 * z / L
        ξg_f(z, t) = 2.0 + sin(π * z / L) * exp(-t / τ)
        fields = MMSFields(T_f, ((z, t) -> 700.0, ξg_f))
        J_hand(z, t) = begin
            T = T_f(z, t)
            dTdz = 50.0 / L
            dξdz = (π / L) * cos(π * z / L) * exp(-t / τ)
            -(λ / T) * (dξdz * T + ξg_f(z, t) * dTdz)
        end
        for z in (0.0, 0.006, 0.015, L), t in (0.0, 1.5)
            @test mms_gas_flux(mat, fields, 2, z, t) ≈ J_hand(z, t) rtol = 1e-12
            @test mms_S_conv(mat, fields, z, t) ≈ -cp_g * J_hand(z, t) * (50.0 / L) rtol =
                1e-12
            # Species source: ∂tξ + ∂zJ, with ∂zJ from a central difference of
            # the hand flux (independent of the ForwardDiff evaluation).
            δ = 1e-7 * L
            dJdz = (J_hand(z + δ, t) - J_hand(z - δ, t)) / (2δ)
            dξdt = -sin(π * z / L) * exp(-t / τ) / τ
            @test mms_species_source(mat, fields, 2, z, t) ≈ dξdt + dJdz rtol = 1e-6
        end
    end

    @testset "end-to-end heat-only MMS (order + Jacobian agreement)" begin
        ρs, cs, ks = 1000.0, 1200.0, 0.4
        L, τ = 0.02, 5.0
        fields = MMSFields(
            (z, t) -> 300.0 + 50.0 * sin(π * z / L) * exp(-t / τ) + 20.0 * z / L,
            ((z, t) -> ρs,),
        )
        function build(n)
            solid = SolidComponent(:s, ρ = ρs, c = cs, k = ks)
            mat = Material(name = :MMSSelfE2E, components = (solid,), reactions = ())
            mesh = create_uniform_mesh(L, n, Val(1); T_initial = 300.0, ξ_initial = (ρs,))
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

        # Analytic Jacobian and colored-FD must agree (the injected source is
        # state-independent, so the Structured Jacobian stays exact).
        sol_an = solve_mms(build(16), fields; jacobian = :analytic)
        sol_fd = solve_mms(build(16), fields; jacobian = :fd)
        @test sol_an.retcode == :Success && sol_fd.retcode == :Success
        @test maximum(abs.(sol_an.T[:, end] - sol_fd.T[:, end])) < 1e-5

        # Temporal order ≈ 4 (KenCarp4): fixed mesh, fixed-dt sweep against a
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
        # γ_s = 1 / γ_g = 0, single first-order s → g, depletion limiter
        # numerically inactive (tanh(ξ_s/1.0) ≡ 1.0 in Float64 for
        # ξ_s ≳ 20 kg/m³ — here ξ_s ≥ 500).
        L = 0.01
        ρs = 800.0
        A_r, E_r = 5.0e4, 80.0e3
        solid = SolidComponent(:s, ρ = ρs, c = 1200.0, k = 0.2, λ = 2.0e-5, γ = 1.0)
        gas = GaseousComponent(:g, M = 0.028, c = 1000.0, k = 0.03, λ = 2.0e-5)  # γ = 0
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
        θ_hand(z, t) = -rate_hand(z, t) / ρs               # γ_g = 0 drops the gas term
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
