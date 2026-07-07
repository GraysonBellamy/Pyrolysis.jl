# =============================================================================
# Exact solutions: semi-infinite solid conduction
# =============================================================================
#
# Requires `erfc` and `erfcx` (SpecialFunctions) in the including namespace.

"""
    semiinf_convective_theta(x, t; h, k, α)

Dimensionless temperature rise `θ = (T − T₀)/(T∞ − T₀)` for a semi-infinite
solid `x ≥ 0`, initially at `T₀`, exchanging heat by surface convection
`q″ = h(T∞ − T_s)` at `x = 0`.

Carslaw & Jaeger, *Conduction of Heat in Solids*, 2nd ed. (1959), §2.7:

    θ = erfc(η) − exp(hx/k + h²αt/k²)·erfc(η + h√(αt)/k),   η = x/(2√(αt))

Evaluated in a numerically stable form: with `B = η + h√(αt)/k` the exponent
satisfies `hx/k + h²αt/k² = B² − η²`, so

    θ = erfc(η) − exp(−η²)·erfcx(B)

which cannot overflow (the naive `exp(·)·erfc(·)` product does, for small `t`
or large `x` at large `h`).

A constant absorbed surface flux `q″_abs` superposed on the convective
exchange is equivalent to convection toward `T∞_eff = T∞ + q″_abs/h` (the
combined BC is linear in `T_s`), so this one function also covers the
Gpyro §3.4.1 configuration (25 kW/m² absorbed + h = 20 W/m²K loss).
"""
function semiinf_convective_theta(x::Real, t::Real; h::Real, k::Real, α::Real)
    t <= 0 && return 0.0
    s = sqrt(α * t)
    η = x / (2 * s)
    B = η + h * s / k
    return erfc(η) - exp(-η^2) * erfcx(B)
end

"""
    T_semiinf_convective(x, t; T0, Tinf, h, k, α)

Temperature at depth `x` below the exposed surface at time `t` for the
semi-infinite convective-BC problem (see [`semiinf_convective_theta`](@ref)).
"""
function T_semiinf_convective(
    x::Real,
    t::Real;
    T0::Real,
    Tinf::Real,
    h::Real,
    k::Real,
    α::Real,
)
    return T0 + (Tinf - T0) * semiinf_convective_theta(x, t; h = h, k = k, α = α)
end
