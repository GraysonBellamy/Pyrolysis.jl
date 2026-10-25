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
    rate0(ξ) = reaction_rate(rxn0, (ξ, zero(ξ)), T_kin; limiter)
    rate03(ξ) = reaction_rate(rxn03, (ξ, zero(ξ)), T_kin; limiter)

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
        @test reaction_rate(rxn0, (0.0, 0.0), T_kin; limiter = no_limiter) ≈ k rtol = 1e-12
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
        n, δ, ξ = 0.3, limiter.base_scale, 1e-12
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
