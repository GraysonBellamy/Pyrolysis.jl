# =============================================================================
# Exact solutions: isothermal reaction kinetics
# =============================================================================
#
# The solver's rate law is r = A·e^(−E/RT)·ξ_ref·xⁿ on the conversion
# x = ξ/ξ_ref, with k = A·e^(−E/RT) in s⁻¹ for every order: a rate constant
# published for dx/dt = −k xⁿ is the solver's as it stands.

"""
    nth_order_xi(t; ξ0, k, n)

Isothermal nth-order decay of a reactant at the conversion `x = ξ/ξ₀`,
`dx/dt = −k xⁿ` (n ≠ 1), Khawam & Flanagan (2006):

    ξ(t) = ξ₀ [1 − (1−n)·k·t]^(1/(1−n))

For n < 1 there is a finite burnout time `t_b = 1/((1−n)k)`
([`nth_order_burnout_time`](@ref)); past it the solution is identically 0
(returned as 0 here).
"""
function nth_order_xi(t::Real; ξ0::Real, k::Real, n::Real)
    n == 1 && throw(ArgumentError("n = 1 is exponential decay; use exp(−kt) directly"))
    base = 1 - (1 - n) * k * t
    base <= 0 && return 0.0  # at/past burnout (only reachable for n < 1)
    return ξ0 * base^(1 / (1 - n))
end

"""
    nth_order_burnout_time(; k, n)

Exact burnout time `t_b = 1 / ((1−n)·k)` for n < 1.
"""
function nth_order_burnout_time(; k::Real, n::Real)
    n < 1 || throw(ArgumentError("finite burnout requires n < 1"))
    return 1 / ((1 - n) * k)
end

"""
    series_chain_xi(t; ξ0, k1, k2)

Isothermal two-step series chain A →(k₁) B →(k₂) C with unit yields
(Bateman): returns `(ξ_A, ξ_B, ξ_C)`.

    ξ_A = ξ₀ e^(−k₁t)
    ξ_B = ξ₀ k₁/(k₂−k₁) (e^(−k₁t) − e^(−k₂t))     (k₁ ≠ k₂)
    ξ_C = ξ₀ − ξ_A − ξ_B
"""
function series_chain_xi(t::Real; ξ0::Real, k1::Real, k2::Real)
    k1 == k2 && throw(ArgumentError("degenerate k₁ = k₂ form not implemented"))
    ξA = ξ0 * exp(-k1 * t)
    ξB = ξ0 * k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t))
    return (ξA, ξB, ξ0 - ξA - ξB)
end

"""
    series_chain_conversion(t; k1, k2, n1 = 1.0, n2, h = 1e-3)

Isothermal two-step series scheme on conversions, the Kinetics Neo convention
each intermediate is normalised to its own lineage:

    de_A/dt = −k₁ e_A^n₁,   de_B/dt = k₁ e_A^n₁ − k₂ e_B^n₂,

`e_A = m_A/m_{A,0}` and `e_B = m_B/(y m_{A,0})`, `y` the yield of B in the first
step, from `e_A = 1`, `e_B = 0`. Returns `(e_A, e_B)` by classical Runge–Kutta at
the step `h`, whose local error is `O((k h)⁵)`.
"""
function series_chain_conversion(
    t::Real;
    k1::Real,
    k2::Real,
    n1::Real = 1.0,
    n2::Real,
    h::Real = 1e-3,
)
    power(e, n) = e > 0 ? e^n : 0.0
    f(a, b) = (-k1 * power(a, n1), k1 * power(a, n1) - k2 * power(b, n2))
    steps = ceil(Int, t / h)
    dt = t / steps
    a, b = 1.0, 0.0
    for _ = 1:steps
        k₁ = f(a, b)
        k₂ = f(a + dt / 2 * k₁[1], b + dt / 2 * k₁[2])
        k₃ = f(a + dt / 2 * k₂[1], b + dt / 2 * k₂[2])
        k₄ = f(a + dt * k₃[1], b + dt * k₃[2])
        a += dt / 6 * (k₁[1] + 2k₂[1] + 2k₃[1] + k₄[1])
        b += dt / 6 * (k₁[2] + 2k₂[2] + 2k₃[2] + k₄[2])
    end
    return (a, b)
end
