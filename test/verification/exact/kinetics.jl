# =============================================================================
# Exact solutions: isothermal reaction kinetics
# =============================================================================
#
# UNIT CONVENTION (pitfall): the solver's rate law is r = A·e^(−E/RT)·ξⁿ with
# ξ the ABSOLUTE concentration [kg/m³] (src/Physics/kinetics.jl), not a
# normalized conversion. Published rate constants k_pub quoted for
# dm/dt = −k_pub·mⁿ with m = ξ/ξ₀ must be converted: A_ours = k_pub·ξ₀^(1−n).
# All functions below work directly in concentration units.

"""
    nth_order_xi(t; ξ0, k, n)

Isothermal nth-order decay `dξ/dt = −k ξⁿ` (n ≠ 1), Khawam & Flanagan (2006):

    ξ(t) = [ξ₀^(1−n) − (1−n)·k·t]^(1/(1−n))

For n < 1 there is a finite burnout time `t_b = ξ₀^(1−n)/((1−n)k)`
([`nth_order_burnout_time`](@ref)); past it the solution is identically 0
(returned as 0 here).
"""
function nth_order_xi(t::Real; ξ0::Real, k::Real, n::Real)
    n == 1 && throw(ArgumentError("n = 1 is exponential decay; use exp(−kt) directly"))
    base = ξ0^(1 - n) - (1 - n) * k * t
    base <= 0 && return 0.0  # at/past burnout (only reachable for n < 1)
    return base^(1 / (1 - n))
end

"""
    nth_order_burnout_time(; ξ0, k, n)

Exact burnout time `t_b = ξ₀^(1−n) / ((1−n)·k)` for n < 1.
"""
function nth_order_burnout_time(; ξ0::Real, k::Real, n::Real)
    n < 1 || throw(ArgumentError("finite burnout requires n < 1"))
    return ξ0^(1 - n) / ((1 - n) * k)
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
