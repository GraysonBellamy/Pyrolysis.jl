# Depletion-limiter tests of the shared reaction rate.
#
# - The limiter is a multiplicative rate factor g(ξ/ξ_th) with no threshold
#   branch: zero at and below zero, the cubic (s − 2y)x³ + (3y − s)x² up to
#   the join (y = tanh 1, s = sech² 1 for the default join of 1), and tanh(x)
#   above, so value and slope are continuous everywhere. A base-applied
#   limiter `ξ_eff = ξ < ξ_th ? ξ·tanh(ξ/ξ_th) : ξ` would jump 23.84% in value
#   at ξ = ξ_th.
# - A base-applied limiter cannot gate order-0 reactants (ξ_eff^0 ≡ 1): the
#   rate would stay k(T) at ξ ≤ 0 and drive ξ negative without bound, and
#   orders 0 < n < 0.5 would keep a singular ∂r/∂ξ at depletion.

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Physics: reaction_rate

# The onset of `DepletionLimiter()` (join 1) and its slope, written out.
const _ONSET_Y = tanh(1.0)
const _ONSET_S = sech(1.0)^2
_onset(x) =
    x <= 0 ? 0.0 :
    x < 1 ? (_ONSET_S - 2_ONSET_Y) * x^3 + (3_ONSET_Y - _ONSET_S) * x^2 : tanh(x)
_onset_slope(x) = 3(_ONSET_S - 2_ONSET_Y) * x^2 + 2(3_ONSET_Y - _ONSET_S) * x

@testset "depletion limiter" begin
    T_kin = 700.0
    ξ_th = 1.0  # DepletionLimiter default threshold
    limiter = DepletionLimiter()
    no_limiter = DepletionLimiter(enabled = false)

    rxn1 = Reaction(:first_order, 1 => 2, A = 1e10, E = 150e3, h = 5e5)
    rxn0 = Reaction(:zero_order, 1 => 2, A = 1e10, E = 150e3, h = 5e5, n = 0.0)
    rxn03 = Reaction(:tga_fit, 1 => 2, A = 1e10, E = 150e3, h = 5e5, n = 0.3)

    # Arrhenius constant without touching R_GAS: n = 1 unlimited rate at ξ = 1.
    k = reaction_rate(rxn1, (1.0, 0.0), T_kin; limiter = no_limiter)
    @test k > 0

    rate1(ξ) = reaction_rate(rxn1, (ξ, zero(ξ)), T_kin; limiter)
    # The other orders against a reference of 1 kg/m³ with the conversion onset
    # at ξ_th, so the same onset reads on the same variable.
    on_conversion = DepletionLimiter(conversion_threshold = ξ_th)
    rate0(ξ) = reaction_rate(rxn0, (ξ, zero(ξ)), T_kin, 1.0; limiter = on_conversion)
    rate03(ξ) = reaction_rate(rxn03, (ξ, zero(ξ)), T_kin, 1.0; limiter = on_conversion)

    @testset "continuous and smooth across the threshold" begin
        # A base-applied limiter jumps by factor tanh(1) ≈ 0.762 crossing ξ_th.
        ε = 1e-9
        @test rate1(ξ_th * (1 + ε)) ≈ rate1(ξ_th * (1 - ε)) rtol = 1e-6
        d_lo = ForwardDiff.derivative(rate1, ξ_th * (1 - ε))
        d_hi = ForwardDiff.derivative(rate1, ξ_th * (1 + ε))
        @test d_lo ≈ d_hi rtol = 1e-6

        # The n = 1 limited rate is k·ξ·g(ξ/ξ_th): the cubic below the join,
        # k·ξ·tanh(ξ/ξ_th) above it.
        for ξ in (0.1, 0.5, 1.0, 3.0, 10.0)
            @test rate1(ξ) ≈ k * ξ * _onset(ξ / ξ_th) rtol = 1e-12
        end
    end

    @testset "bare Arrhenius rate recovered far above threshold" begin
        @test rate1(10 * ξ_th) ≈ k * 10 * ξ_th rtol = 1e-8
        # Saturated onset (ξ ≥ 20·ξ_th) is bit-identical to the unlimited rate.
        ξ_sat = 25 * ξ_th
        @test rate1(ξ_sat) == reaction_rate(rxn1, (ξ_sat, 0.0), T_kin; limiter = no_limiter)
        # AD through the saturated branch: dr/dξ = k exactly.
        @test ForwardDiff.derivative(rate1, ξ_sat) ≈ k rtol = 1e-12
    end

    @testset "zero-order reactions stop at depletion" begin
        # Without the multiplicative gate the rate stays k(T) at ξ = 0 and
        # beyond (unbounded negative ξ).
        @test rate0(0.0) == 0.0
        @test rate0(-1.0) == 0.0  # integrator overshoot consumes nothing
        # Near depletion the rate decays as k(3y − s)ξ², so the reactant
        # approaches zero without crossing it.
        @test rate0(1e-3) ≈ k * _onset(1e-3 / ξ_th) rtol = 1e-12
        # ∂r/∂ξ just above depletion is finite and vanishes with ξ.
        @test ForwardDiff.derivative(rate0, 1e-12) ≈ k * _onset_slope(1e-12) / ξ_th rtol =
            1e-12
        # Unlimited zero-order rate is k regardless of ξ — the hazard the
        # limiter exists to remove.
        @test reaction_rate(rxn0, (0.0, 0.0), T_kin, 1.0; limiter = no_limiter) ≈ k rtol =
            1e-12
    end

    @testset "bounded ∂r/∂ξ at depletion for n < 0.5" begin
        # A base-applied limiter gives ∂r/∂ξ ~ ξ^(2n−1) → ∞ for n < 0.5.
        # The multiplicative form gives r = k·b(ξ)ⁿ·g(ξ/ξ_th), so the
        # derivative vanishes at depletion instead of diverging.
        d4 = ForwardDiff.derivative(rate03, 1e-4)
        d8 = ForwardDiff.derivative(rate03, 1e-8)
        d12 = ForwardDiff.derivative(rate03, 1e-12)
        @test all(isfinite, (d4, d8, d12))
        @test d12 < d8 < d4
        # Small-ξ slope by the chain rule, with the finite base
        # b(ξ) = ξ + δ·log(1 + e^{−ξ/δ}) of slope 1/(1 + e^{−ξ/δ}), δ = 1e-8.
        n, δ, ξ = 0.3, on_conversion.conversion_base, 1e-12
        b = ξ + δ * log1p(exp(-ξ / δ))
        db = 1 / (1 + exp(-ξ / δ))
        @test d12 ≈
              k *
              (n * b^(n - 1) * db * _onset(ξ / ξ_th) + b^n * _onset_slope(ξ / ξ_th) / ξ_th) rtol =
            1e-9
    end

    @testset "temperature gate width" begin
        # Gate factor is 0.5·(1 + tanh((T − T_min)/gate_width)); default width
        # 1.0 K is the bare tanh(T − T_min) ramp bit-for-bit.
        T_min = 600.0
        rxn_default =
            Reaction(:gate_default, 1 => 2, A = 1e10, E = 150e3, h = 5e5, T_min = T_min)
        @test rxn_default.gate_width == 1.0
        r_ref(rxn, T) = reaction_rate(rxn, (1.0, 0.0), T; limiter = no_limiter)
        r_open(T) = reaction_rate(rxn1, (1.0, 0.0), T; limiter = no_limiter)
        for T in (T_min - 3.0, T_min, T_min + 3.0)
            @test r_ref(rxn_default, T) ≈ r_open(T) * 0.5 * (1 + tanh(T - T_min)) rtol =
                1e-14
        end

        # Narrow gate sharpens the ramp around T_min without moving it.
        w = 0.1
        rxn_narrow = Reaction(
            :gate_narrow,
            1 => 2,
            A = 1e10,
            E = 150e3,
            h = 5e5,
            T_min = T_min,
            gate_width = w,
        )
        for T in (T_min - 0.5, T_min - 0.05, T_min, T_min + 0.05, T_min + 0.5)
            @test r_ref(rxn_narrow, T) ≈ r_open(T) * 0.5 * (1 + tanh((T - T_min) / w)) rtol =
                1e-14
        end
        # Midpoint pinned at T_min for any width; 5 widths out the gate is open.
        @test r_ref(rxn_narrow, T_min) ≈ 0.5 * r_open(T_min) rtol = 1e-12
        @test r_ref(rxn_narrow, T_min + 5w) ≈ r_open(T_min + 5w) * 0.5 * (1 + tanh(5.0)) rtol =
            1e-12

        # T_max gate uses the same width (mirrored form).
        T_max = 900.0
        rxn_tmax = Reaction(
            :gate_tmax,
            1 => 2,
            A = 1e10,
            E = 150e3,
            h = 5e5,
            T_max = T_max,
            gate_width = w,
        )
        rxn_tmax_open = Reaction(:gate_tmax_open, 1 => 2, A = 1e10, E = 150e3, h = 5e5)
        for T in (T_max - 0.5, T_max, T_max + 0.5)
            @test r_ref(rxn_tmax, T) ≈
                  r_ref(rxn_tmax_open, T) * 0.5 * (1 + tanh((T_max - T) / w)) rtol = 1e-14
        end

        # AD through T stays finite and matches the analytic sech² slope at
        # the gate midpoint (the width threads through ForwardDiff untouched).
        rate_T(T) = reaction_rate(rxn_narrow, (1.0, 0.0), T; limiter = no_limiter)
        dT = ForwardDiff.derivative(rate_T, T_min)
        k_arr(T) = reaction_rate(rxn1, (1.0, 0.0), T; limiter = no_limiter)
        dk = ForwardDiff.derivative(k_arr, T_min)
        @test dT ≈ 0.5 * dk + k_arr(T_min) * 0.5 / w rtol = 1e-10

        # Construction rejects degenerate widths, numeric and symbolic forms.
        for bad in (0.0, -1.0, Inf, NaN)
            @test_throws ArgumentError Reaction(
                :bad_gate,
                1 => 2,
                A = 1e10,
                E = 150e3,
                h = 5e5,
                gate_width = bad,
            )
            @test_throws ArgumentError Reaction(
                :bad_gate,
                :virgin => :gas,
                A = 1e10,
                E = 150e3,
                h = 5e5,
                gate_width = bad,
            )
            @test_throws ArgumentError Reaction(
                :bad_gate,
                :virgin => (:char, :gas),
                A = 1e10,
                E = 150e3,
                h = 5e5,
                yields = (0.3, 0.7),
                gate_width = bad,
            )
        end

        # Symbolic constructors thread the width through Material resolution.
        mat = Material(
            name = :gate_thread,
            components = (
                SolidComponent(:virgin; λ = 0.0, ρ = 1000.0, c = 1000.0, k = 0.2),
                GaseousComponent(:gas; M = 0.05, c = 1000.0, k = 0.02),
            ),
            reactions = (
                Reaction(
                    :melt,
                    :virgin => :gas,
                    A = 1e10,
                    E = 150e3,
                    h = 5e5,
                    T_min = T_min,
                    gate_width = w,
                ),
            ),
        )
        @test mat.reactions[1].gate_width == w
    end

    @testset "construction rejects negative orders" begin
        @test_throws ArgumentError Reaction(
            :bad,
            1 => 2,
            A = 1e10,
            E = 150e3,
            h = 5e5,
            n = -0.5,
        )
        @test_throws ArgumentError Reaction(
            :bad,
            :virgin => :gas,
            A = 1e10,
            E = 150e3,
            h = 5e5,
            n = -1.0,
        )
        @test_throws ArgumentError Reaction(
            :bad,
            :virgin => (:char, :gas),
            A = 1e10,
            E = 150e3,
            h = 5e5,
            yields = (0.3, 0.7),
            n = -0.1,
        )
        # n = 0 and fractional orders are valid fits.
        @test rxn0 isa Reaction
        @test rxn03 isa Reaction
    end
end

@testset "conversion basis" begin
    T = 700.0
    limiter = DepletionLimiter()
    bare = DepletionLimiter(conversion_threshold = 1e-12)
    no_limiter = DepletionLimiter(enabled = false)
    rxn(n) = Reaction(:r, 1 => 2, A = 1e10, E = 150e3, h = 5e5, n = n)
    k = reaction_rate(rxn(1.0), (1.0, 0.0), T; limiter = no_limiter)
    x_th = limiter.conversion_threshold
    ξ_ref = 200.0
    rate(n, x; limiter = limiter) =
        reaction_rate(rxn(n), (x * ξ_ref, zero(x)), T, ξ_ref; limiter)

    @testset "the rate on the conversion" begin
        for n in (0.0, 0.5, 2.0), x in (0.9, 0.3, 0.01)
            @test rate(n, x) ≈
                  k * ξ_ref * x^n * Pyrolysis.Physics.one_sided_factor(x / x_th, 1.0)
        end
        # The reference cancels at order one: the same rate to the bit.
        for ξ in (-1.0, 0.0, 1e-12, 1e-3, 0.5, 1.0, 10.0, 1000.0), ξ_ref in (1e-3, 1.0, 1e3)
            a = reaction_rate(rxn(1.0), (ξ, 0.0), T; limiter)
            b = reaction_rate(rxn(1.0), (ξ, 0.0), T, ξ_ref; limiter)
            @test reinterpret(UInt64, a) == reinterpret(UInt64, b)
        end
        # At order one the reference is not read: zero is the kernel's slot.
        @test reaction_rate(rxn(1.0), (1.0, 0.0), T, 0.0) ==
              reaction_rate(rxn(1.0), (1.0, 0.0), T)
    end

    @testset "C¹ across the threshold and the join" begin
        for n in (0.0, 0.5, 2.0), x in (x_th, 3x_th)
            ε = 1e-9
            @test rate(n, x * (1 + ε)) ≈ rate(n, x * (1 - ε)) rtol = 1e-6
            d_lo = ForwardDiff.derivative(s -> rate(n, s), x * (1 - ε))
            d_hi = ForwardDiff.derivative(s -> rate(n, s), x * (1 + ε))
            @test d_lo ≈ d_hi rtol = 1e-6
        end
    end

    @testset "the bare rate above 20 x_th" begin
        for n in (0.0, 0.5, 2.0), x in (20x_th, 0.5, 1.0)
            @test rate(n, x) ≈ k * ξ_ref * x^n rtol = 1e-14
        end
    end

    @testset "both asymptotic ranges" begin
        y, s = tanh(1.0), sech(1.0)^2
        C = 3y - s
        δ = limiter.conversion_base
        @test C ≈ 1.8648 atol = 5e-5
        for n in (0.0, 0.5, 2.0)
            x = 1e-6
            @test rate(n, x) ≈ C * k * ξ_ref * x^(n + 2) / x_th^2 rtol = 1e-3
            x = 1e-12
            @test rate(n, x) ≈ C * k * ξ_ref * (δ * log(2))^n * x^2 / x_th^2 rtol = 1e-3
            d = ForwardDiff.derivative(s -> rate(n, s), 1e-12)
            @test 0 < d < 1e-5 * k * ξ_ref
        end
    end

    @testset "exactly zero at and below zero mass" begin
        for n in (0.0, 0.5, 2.0), x in (0.0, -0.0, -1e-6, -1.0)
            @test rate(n, x) == 0
            @test ForwardDiff.derivative(s -> rate(n, s), x) == 0
        end
    end

    @testset "the limiter disabled leaves the finite base" begin
        @test rate(2.0, 0.5; limiter = no_limiter) ≈ k * ξ_ref * 0.25
        @test rate(0.0, 0.5; limiter = no_limiter) ≈ k * ξ_ref
        @test rate(0.0, -1.0; limiter = no_limiter) ≈ k * ξ_ref
        @test isfinite(rate(0.5, -1.0; limiter = no_limiter))
        @test isfinite(ForwardDiff.derivative(s -> rate(0.5, s; limiter = no_limiter), 0.0))
    end

    @testset "the reference must be finite and positive" begin
        for ξ_ref in (0.0, -1.0, NaN, Inf)
            @test_throws ArgumentError reaction_rate(rxn(2.0), (1.0, 0.0), T, ξ_ref)
        end
    end

    @testset "an order other than one needs a reference" begin
        for n in (0.0, 0.5, 2.0)
            @test_throws "needs a finite positive reference" reaction_rate(
                rxn(n),
                (1.0, 0.0),
                T,
            )
        end
    end

    @testset "limiter scales" begin
        @test limiter.conversion_threshold == 1e-3
        @test limiter.conversion_base == 1e-8
        @test_throws ArgumentError DepletionLimiter(conversion_threshold = 0.0)
        @test_throws ArgumentError DepletionLimiter(conversion_base = -1.0)
        @test bare.conversion_threshold == 1e-12
    end

    @testset "basis" begin
        kw = (A = 1.0, E = 0.0, h = 0.0)
        @test isnan(rxn(2.0).basis)
        @test Reaction(:r, 1 => 2; kw..., n = 2.0, basis = 0.5).basis == 0.5
        @test Reaction(:r, 1 => 2; kw..., n = 2.0, basis = 1).basis == 1.0
        @test Reaction(:r, 1 => (2, 3); kw..., yields = (0.5, 0.5), n = 2.0, basis = 0.5).basis ==
              0.5
        for basis in (0.0, -0.5, 1.5, Inf)
            @test_throws "(0, 1]" Reaction(:r, 1 => 2; kw..., n = 2.0, basis)
            @test_throws "(0, 1]" Reaction(:r, :a => :b; kw..., n = 2.0, basis)
        end
        @test_throws "first order" Reaction(:r, 1 => 2; kw..., basis = 0.5)
        @test_throws "first order" Reaction(:r, :a => :b; kw..., basis = 0.5)
        @test_throws "first order" Reaction(
            :r,
            :a => (:b, :c);
            kw...,
            yields = (0.5, 0.5),
            basis = 0.5,
        )
        solid = SolidComponent(:a, ρ = 500.0, c = 1000.0, k = 0.2, λ = 1e-6)
        gas = GaseousComponent(:b, M = 0.03, c = 1000.0, k = 0.02)
        m = Material(
            name = :m,
            components = (solid, gas),
            reactions = (Reaction(:r, :a => :b; kw..., n = 2.0, basis = 0.5),),
        )
        @test only(m.reactions).basis == 0.5
    end
end
