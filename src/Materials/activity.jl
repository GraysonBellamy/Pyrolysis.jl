"""
Activity laws of a liquid held in sorbing solids, and their heats of sorption.
"""

"""
	AbstractActivity

Supertype of the activity laws of a phase change. A law gives the activity
`a(X, T) = p_eq/p_sat(T)` of a liquid at moisture content `X` [kg per kg of
sorbing solid] and temperature `T` [K], through [`activity`](@ref), and its
heat of sorption, through [`sorption_heat`](@ref).

Every law is defined at every real `X`, negative trial values included, takes
the activity to zero with the liquid (`a = 0` at `X = 0`) and keeps it at or
below one.
"""
abstract type AbstractActivity end

"""
	activity(law::AbstractActivity, X, T) -> a

Activity of the liquid at moisture content `X` [kg/kg] and temperature `T`
[K]: the ratio of its equilibrium vapour pressure to the saturation pressure.
Signed with `X`. Generic in the number type.
"""
function activity end

"""
	sorption_heat(law::AbstractActivity, X, T, M) -> Q_s [J/kg]

Differential heat of sorption, the heat per kilogram of liquid released on
sorption above the latent heat of the free liquid: `(R T²/M) ∂ln a/∂T` at
constant `X`, with `M` [kg/mol] the molar mass of the vapour. Zero for a free
liquid. Generic in the number type.
"""
function sorption_heat end

# Defaults of the cap at one and of the factor that takes the activity to zero.
const ACTIVITY_CAP_WIDTH = 0.01
const ACTIVITY_VANISHING_SCALE = 1e-4

# Beyond this distance from the cap, in units of its width, the cap is the
# identity (below) or one (above) to round-off.
const ACTIVITY_CAP_CUTOFF = 40

function _check_activity_constants(name, δ, X0)
    isfinite(δ) && δ > 0 ||
        throw(ArgumentError("$name: the cap width δ must be finite and positive, got $δ"))
    isfinite(X0) && X0 > 0 || throw(
        ArgumentError("$name: the vanishing scale X0 must be finite and positive, got $X0"),
    )
    return nothing
end

# The factor X/(|X| + X0): zero at X = 0, signed with X, one to within X0/X.
@inline _vanishing_factor(X, X0) = X / (abs(X) + X0)

# The soft cap at one, smin_δ(x, 1) = 1 − δ ln(1 + e^{(1−x)/δ}): the argument
# itself where z = (1 − x)/δ > 40 and one where z < −40. The correction left out,
# δ ln(1 + e^{−|z|}), is below δ e^{−40} (4 × 10⁻²⁰ at δ = 0.01).
@inline function _soft_cap(x, δ)
    z = (1 - x) / δ
    z > ACTIVITY_CAP_CUTOFF && return x
    z < -ACTIVITY_CAP_CUTOFF && return one(x)
    return 1 - δ * log1p(exp(z))
end

# The heat of the capped activity from the heat q of the raw activity x: q times
# the logarithmic slope of the cap, x σ((1 − x)/δ)/smin_δ(x, 1), with the cap's
# two short cuts (q where the cap is the identity, zero where it is one). The
# product is formed without e^{(1−x)/δ}, which overflows far above the cap. The
# floor on the cap binds only for a cap wider than 1/40, whose join reaches x ≤ 0.
@inline function _capped_heat(q, x, δ)
    z = (1 - x) / δ
    z > ACTIVITY_CAP_CUTOFF && return q
    z < -ACTIVITY_CAP_CUTOFF && return zero(q)
    return q * x / (1 + exp(-z)) / max(_soft_cap(x, δ), 1e-12)
end

# ============================================================================
# CAST isotherm
# ============================================================================

"""
	CASTIsotherm(β₁, β₂, β₃, β₄, β₅, β₆; δ = 0.01, X0 = 1e-4)

The Comprehensive Analytical Sorption Thermodynamic isotherm of Dietenberger,
Glass and Boardman (2024), for water in wood and other cellulosic solids. The
differential Gibbs energy of sorption (their eq. 18), in kJ/mol with `u` the
moisture content [kg/kg] and `T` in K,

	Δg(u, T) = [u (β₁ e^{−β₂/T} + β₃) + β₄]^{−β₅} − β₆,

gives the activity `exp(−Δg/(R T))` (their eq. 27 without the fugacity
correction: the activity is read as `p_v/p_sat`) and the differential enthalpy
of sorption by the Gibbs–Helmholtz relation (their eq. 19),

	Δh = Δg + (Δg + β₆)/T · u β₁ β₂ β₅ e^{−β₂/T} / [u (β₁ e^{−β₂/T} + β₃) + β₄],

the heat of sorption being `Δh/M`. Two additions make it a law for every
state:

- the activity is multiplied by `X/(|X| + X0)`, which takes it to zero with the
  water and makes it negative below (the printed form keeps `a ≈ 0.002–0.024`
  at `u = 0` between 25 and 150 °C); `Δg` is evaluated at `u = max(X, 0)`, so
  the factor carries the sign and the activity is C¹ through `X = 0`;
- the result is capped at one by `1 − δ ln(1 + e^{(1 − a)/δ})`. Above the
  moisture content at which `Δg = 0` (their eq. 25: 29 % at 25 °C, 19 % at
  100 °C and 14 % at 150 °C for [`CAST_SPRUCE`](@ref)) the water is free and
  its vapour saturated. The heat of sorption carries the cap's logarithmic
  slope and falls to zero there.

The β are positive (the paper's statement of the model). The constructor
refuses a β that is not finite and positive, and a `δ` or `X0` that is not.
"""
struct CASTIsotherm <: AbstractActivity
    β::NTuple{6,Float64}
    δ::Float64
    X0::Float64
    function CASTIsotherm(
        β₁::Real,
        β₂::Real,
        β₃::Real,
        β₄::Real,
        β₅::Real,
        β₆::Real;
        δ::Real = ACTIVITY_CAP_WIDTH,
        X0::Real = ACTIVITY_VANISHING_SCALE,
    )
        β = Float64.((β₁, β₂, β₃, β₄, β₅, β₆))
        all(b -> isfinite(b) && b > 0, β) || throw(
            ArgumentError("CASTIsotherm: β₁ … β₆ must be finite and positive, got $β"),
        )
        _check_activity_constants(:CASTIsotherm, δ, X0)
        return new(β, Float64(δ), Float64(X0))
    end
end

"""
	CAST_SPRUCE

[`CASTIsotherm`](@ref) fitted to spruce by Dietenberger, Glass and Boardman
(2024, Table 3): `β = (304.2, 1661, 3.520, 0.3895, 2.886, 0.2006)`, from
Weichert's isotherms at 25–100 °C and Volbehr's heats of wetting. The
softwood set.
"""
const CAST_SPRUCE = CASTIsotherm(304.2, 1661.0, 3.520, 0.3895, 2.886, 0.2006)

"""
	CAST_BEECH

[`CASTIsotherm`](@ref) fitted to beech by Dietenberger, Glass and Boardman
(2024, Table 3): `β = (46.33, 777.6, 2.028, 0.2210, 1.725, 0.3664)`, from
Weichert's isotherms at 25–100 °C, Engelhardt's at 110–170 °C and the heats of
wetting of Hearmon and Burcham and of Nopens et al. The hardwood set.
"""
const CAST_BEECH = CASTIsotherm(46.33, 777.6, 2.028, 0.2210, 1.725, 0.3664)

# The moisture content max(X, 0) that enters the bracket, the Arrhenius factor
# e^{−β₂/T}, the bracket u (β₁ e^{−β₂/T} + β₃) + β₄ and Δg [J/mol].
@inline function _cast_terms(law::CASTIsotherm, X, T)
    β₁, β₂, β₃, β₄, β₅, β₆ = law.β
    u = max(X, zero(X))
    e = exp(-β₂ / T)
    B = u * (β₁ * e + β₃) + β₄
    return u, e, B, 1e3 * (B^(-β₅) - β₆)
end

@inline function _raw_activity(law::CASTIsotherm, X, T)
    _, _, _, g = _cast_terms(law, X, T)
    return exp(-g / (R_GAS * T)) * _vanishing_factor(X, law.X0)
end

@inline activity(law::CASTIsotherm, X, T) = _soft_cap(_raw_activity(law, X, T), law.δ)

@inline function sorption_heat(law::CASTIsotherm, X, T, M)
    β₁, β₂, _, _, β₅, β₆ = law.β
    u, e, B, g = _cast_terms(law, X, T)
    Δh = g + (g + 1e3 * β₆) / T * (u * β₁ * β₂ * β₅ * e / B)
    x = exp(-g / (R_GAS * T)) * _vanishing_factor(X, law.X0)
    return _capped_heat(Δh / M, x, law.δ)
end

# ============================================================================
# Free liquid
# ============================================================================

"""
	FreeLiquid(; X0 = 1e-4)

A liquid not held by the solid: its vapour is saturated whatever the moisture
content, down to the last of the liquid. The activity is

	a = X/(|X| + X0),

one to within `X0/X` (0.04 % at `X = 0.25`, 1 % at `X = 0.01`), never above
one, zero at `X = 0` and negative below, C¹ through zero. There is no cap and
no heat of sorption. The constructor refuses an `X0` that is not finite and
positive.
"""
struct FreeLiquid <: AbstractActivity
    X0::Float64
    function FreeLiquid(; X0::Real = ACTIVITY_VANISHING_SCALE)
        isfinite(X0) && X0 > 0 || throw(
            ArgumentError(
                "FreeLiquid: the vanishing scale X0 must be finite and positive, got $X0",
            ),
        )
        return new(Float64(X0))
    end
end

@inline activity(law::FreeLiquid, X, T) = _vanishing_factor(X, law.X0)

@inline sorption_heat(::FreeLiquid, X, T, M) = zero(X) * zero(T) * zero(M)

# ============================================================================
# Activity given as a function
# ============================================================================

"""
	ActivityFunction(f; δ = 0.01, X0 = 1e-4)

An activity law given as a function `f(X, T)` of the moisture content [kg/kg]
and the temperature [K], such as a measured isotherm. As for
[`CASTIsotherm`](@ref), `f` is evaluated at `max(X, 0)`, multiplied by
`X/(|X| + X0)` and capped at one by `1 − δ ln(1 + e^{(1 − a)/δ})`. Below
`X = 0` the activity is negative where `f(0, T) > 0` and zero where
`f(0, T) = 0`.

The heat of sorption is `(R T²/M) ∂ln f/∂T` by forward differentiation in `T`,
times the slope of the cap; where `f ≤ 0` it is zero. A law independent of
temperature has none, which is how an isothermal problem is posed.

`f` must accept dual numbers in both arguments. The constructor refuses a `δ`
or `X0` that is not finite and positive.

```julia
ActivityFunction((X, T) -> X / 0.2)   # linear in X, no heat of sorption
```
"""
struct ActivityFunction{F} <: AbstractActivity
    f::F
    δ::Float64
    X0::Float64
    function ActivityFunction(
        f::F;
        δ::Real = ACTIVITY_CAP_WIDTH,
        X0::Real = ACTIVITY_VANISHING_SCALE,
    ) where {F}
        _check_activity_constants(:ActivityFunction, δ, X0)
        return new{F}(f, Float64(δ), Float64(X0))
    end
end

@inline activity(law::ActivityFunction, X, T) =
    _soft_cap(law.f(max(X, zero(X)), T) * _vanishing_factor(X, law.X0), law.δ)

@inline function sorption_heat(law::ActivityFunction, X, T, M)
    u = max(X, zero(X))
    f = law.f(u, T)
    w = R_GAS * T^2 / M * ForwardDiff.derivative(t -> law.f(u, t), T)
    q = f > 0 ? w / f : zero(w) * zero(f)
    return _capped_heat(q, f * _vanishing_factor(X, law.X0), law.δ)
end
