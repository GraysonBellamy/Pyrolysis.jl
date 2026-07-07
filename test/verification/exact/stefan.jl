# =============================================================================
# Neumann two-phase Stefan problem (solidification front, exact solution)
# =============================================================================
#
# Carslaw & Jaeger §11.2; transcribed from Lautenberger dissertation (2007)
# §3.4.3 Eqs. 3.137–3.138 (Gpyro's flagship moving-front verification).
#
# Semi-infinite medium, liquid everywhere at T₀ > T_m at t = 0. The wall
# (depth x = 0) is dropped to T∞ < T_m; a solidification front s(t) propagates
# inward, releasing the latent heat ΔH_m. With αₛ, α_l the solid/liquid
# diffusivities:
#
#   s(t)      = 2 λ √(αₛ t)                                       (3.138a)
#   T_s(x,t)  = T∞ + (T_m − T∞) erf(x / 2√(αₛ t)) / erf(λ),  x ≤ s  (3.138b)
#   T_l(x,t)  = T₀ + (T_m − T₀) erfc(x / 2√(α_l t)) / erfc(λ√(αₛ/α_l)), x ≥ s (3.138c)
#
# and λ is the root of the transcendental eigencondition (3.138d):
#
#   exp(−λ²)/erf(λ)
#     − (√αₛ k_l (T₀ − T_m)) / (√α_l k_s (T_m − T∞))
#         · exp(−λ² αₛ/α_l) / erfc(λ √(αₛ/α_l))
#     = √π λ ΔH_m / (c_s (T_m − T∞))
#
# The classical Neumann solution assumes equal phase densities ρ (no flow).
#
# Self-tests (selftests.jl): λ = 0.2317 with the dissertation's
# polypropylene-like parameters; the Stefan flux-jump condition and both PDEs
# via finite differences; and a global energy balance by quadrature (wall heat
# extraction = latent release + sensible change), which pins the whole
# solution family — including λ — against an independent computation.

using SpecialFunctions: erf, erfc

"""
Parameter bundle for the Neumann two-phase Stefan problem. Fields are the
solid/liquid conductivities and diffusivities, common density, melting data,
and the boundary/initial temperatures (wall `T∞` < `T_m` < initial `T0`).
"""
Base.@kwdef struct StefanProblem
    k_s::Float64
    k_l::Float64
    ρ::Float64
    c_s::Float64
    c_l::Float64
    T_m::Float64
    ΔH_m::Float64
    T0::Float64      # initial liquid temperature (> T_m)
    T∞::Float64      # wall temperature (< T_m)
end

α_solid(p::StefanProblem) = p.k_s / (p.ρ * p.c_s)
α_liquid(p::StefanProblem) = p.k_l / (p.ρ * p.c_l)

"""
    stefan_eigencondition(λ, p::StefanProblem) -> residual

LHS − RHS of eigencondition (3.138d); the physical `λ` is its positive root.
"""
function stefan_eigencondition(λ::Real, p::StefanProblem)
    αs, αl = α_solid(p), α_liquid(p)
    ν = sqrt(αs / αl)
    lhs =
        exp(-λ^2) / erf(λ) -
        (sqrt(αs) * p.k_l * (p.T0 - p.T_m)) / (sqrt(αl) * p.k_s * (p.T_m - p.T∞)) *
        exp(-λ^2 * ν^2) / erfc(λ * ν)
    rhs = sqrt(π) * λ * p.ΔH_m / (p.c_s * (p.T_m - p.T∞))
    return lhs - rhs
end

"""
    stefan_lambda(p::StefanProblem; hi = 3.0) -> λ

Root of the eigencondition by bisection to ~1e-14. The LHS − RHS residual is
+∞ as λ → 0⁺ and decreases monotonically through the root, so the bracket
`[eps, hi]` is robust for any physical parameter set with T∞ < T_m < T0.
"""
function stefan_lambda(p::StefanProblem; hi::Float64 = 3.0)
    f(λ) = stefan_eigencondition(λ, p)
    lo = 1e-8
    f(lo) > 0 || throw(ArgumentError("eigencondition not positive at λ → 0⁺"))
    f(hi) < 0 || throw(ArgumentError("eigencondition has no sign change on [$lo, $hi]"))
    return bisect_root(f, lo, hi)
end

"""
    stefan_front(t; λ, αs) -> s(t) [m]

Front position `s(t) = 2λ√(αs t)`.
"""
stefan_front(t::Real; λ::Real, αs::Real) = 2.0 * λ * sqrt(αs * t)

"""
    stefan_T(x, t, p::StefanProblem, λ) -> T [K]

Exact temperature at depth `x ≥ 0` below the cold wall at time `t > 0`
(piecewise: solid branch for `x ≤ s(t)`, liquid branch beyond).
"""
function stefan_T(x::Real, t::Real, p::StefanProblem, λ::Real)
    αs, αl = α_solid(p), α_liquid(p)
    s = stefan_front(t; λ = λ, αs = αs)
    if x <= s
        return p.T∞ + (p.T_m - p.T∞) * erf(x / (2 * sqrt(αs * t))) / erf(λ)
    else
        ν = sqrt(αs / αl)
        return p.T0 + (p.T_m - p.T0) * erfc(x / (2 * sqrt(αl * t))) / erfc(λ * ν)
    end
end

"""
Dissertation §3.4.3 parameter set (polypropylene-like), for which
Eq. 3.138d gives λ = 0.2317.
"""
stefan_gpyro_problem() = StefanProblem(
    k_s = 0.20,
    k_l = 0.20,
    ρ = 900.0,
    c_s = 1700.0,
    c_l = 1700.0,
    T_m = 433.0,
    ΔH_m = 50.0e3,
    T0 = 438.0,
    T∞ = 428.0,
)
