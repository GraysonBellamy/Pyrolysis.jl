# =============================================================================
# Exact solutions: thermally-thin radiative cooling and the steady
# radiative–convective surface balance
# =============================================================================

"""
    radiative_cooling_G(T, T∞)

Antiderivative `G(T) = ∫ dT/(T⁴ − T∞⁴)` on the cooling branch `T > T∞ > 0`,
by partial fractions (re-derived; self-tested against quadrature):

    1/(T⁴−a⁴) = (1/2a²)·[1/(T²−a²) − 1/(T²+a²)]
    G(T) = (1/4a³)·ln((T−a)/(T+a)) − (1/2a³)·arctan(T/a),   a = T∞

The lumped radiative-cooling law `dT/dt = −K(T⁴ − T∞⁴)` then has the exact
implicit solution `t = (G(Tᵢ) − G(T))/K` (G is negative-valued and increasing
toward T∞... note dG/dT = 1/(T⁴−a⁴) > 0 on the branch).
"""
function radiative_cooling_G(T::Real, T∞::Real)
    T > T∞ > 0 || throw(DomainError((T, T∞), "cooling branch requires T > T∞ > 0"))
    a = T∞
    return log((T - a) / (T + a)) / (4 * a^3) - atan(T / a) / (2 * a^3)
end

"""
    radiative_cooling_T(t; Ti, T∞, K)

Temperature at time `t` for `dT/dt = −K(T⁴ − T∞⁴)`, `T(0) = Tᵢ`, inverted
from the implicit map by bisection on the monotone residual.
"""
function radiative_cooling_T(t::Real; Ti::Real, T∞::Real, K::Real)
    t <= 0 && return float(Ti)
    G_target = radiative_cooling_G(Ti, T∞) - K * t
    f(T) = radiative_cooling_G(T, T∞) - G_target
    return bisect_root(f, T∞ * (1 + 1e-12), float(Ti))
end

"""
    radiative_convective_equilibrium_Ts(; q_abs, ε, h, T∞, σ)

Steady surface temperature of an adiabatic-backed slab under constant
absorbed flux with radiative + convective losses — the algebraic root of

    q_abs = ε·σ·(Ts⁴ − T∞⁴) + h·(Ts − T∞)

bisected on [T∞, 5000] (LHS − RHS is strictly decreasing in Ts).
"""
function radiative_convective_equilibrium_Ts(;
    q_abs::Real,
    ε::Real,
    h::Real,
    T∞::Real,
    σ::Real,
)
    f(Ts) = q_abs - ε * σ * (Ts^4 - T∞^4) - h * (Ts - T∞)
    return bisect_root(f, float(T∞), 5000.0)
end
