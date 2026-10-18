# ============================================================================
# Blowing models: how outgassing through a face changes its convective
# coefficient. The model is a field of `ConvectiveBC` (default `NoBlowing()`),
# so the physics is declared where the surface is. Evaluation lives in
# `BoundaryConditions/blowing.jl`; this file holds the types, their constants and
# the gas-property data they need.
# ============================================================================

"""
	AbstractBlowingModel

A rule for the ratio `χ = h/h₀` of the convective coefficient with outflow to
the coefficient without it. Each `ConvectiveBC` carries one, and the surface
Newton applies it to the outward gas flux at that face at every iterate.
Concrete models: [`NoBlowing`](@ref), [`CouetteFilm`](@ref),
[`PropertyRatioFilm`](@ref), [`PoolingFilm`](@ref), [`FunctionBlowing`](@ref).
"""
abstract type AbstractBlowingModel end

# ----------------------------------------------------------------------------
# Gas property data
# ----------------------------------------------------------------------------

"""
	GasSpecies(name; M, c, k, μ = nothing, V_D = NaN, D = nothing, T_range = (200.0, 3000.0))

Transport and thermal properties of one gas, as the blowing models need them:
the free stream the face exchanges heat with (purge or ambient gas), or a
lumped description of the volatiles leaving the face.

# Arguments
- `name::Symbol`: identifier. A material component with the same name is
  treated as this gas (used to recognise purge gas venting from the pores).
- `M`: molar mass [kg/mol]
- `c`: heat capacity [J/(kg·K)] (number or property function of `T`)
- `k`: thermal conductivity [W/(m·K)] (number or property function of `T`)
- `μ = nothing`: dynamic viscosity [Pa·s]; `nothing` uses the Eucken estimate
  from `k`, `c` and `M` ([`eucken_viscosity`](@ref)). Only the mixing rule for
  the conductivity of a multi-species outflow reads it.
- `V_D = NaN`: Fuller diffusion volume [-] (Poling, Prausnitz and O'Connell
  Table 11-1), for the binary diffusivity into the free stream
  ([`fuller_diffusivity`](@ref)).
- `D = nothing`: binary diffusivity into the free stream [m²/s] (number or
  property function of `T`); overrides the Fuller estimate.
- `T_range = (200.0, 3000.0)`: temperature range [K] the property laws are
  valid on. The blowing models evaluate them at the surface-temperature
  iterate held inside this range, so a polynomial fit is never extrapolated
  by a wandering Newton iterate. `c` and `k` must be positive on it.

# Example
```julia
NITROGEN                       # shipped: NASA-7 cp, CFD-case k(T) and μ(T), V_D 18.5
mma = GasSpecies(:MMA; M = 0.1001, c = 1500.0, k = 0.03, V_D = 110.2)
```
"""
struct GasSpecies{Pc<:AbstractPropertyFunction,Pk<:AbstractPropertyFunction,Pμ,PD}
    name::Symbol
    molar_mass::Float64
    heat_capacity::Pc
    thermal_conductivity::Pk
    viscosity::Pμ
    diffusion_volume::Float64
    diffusivity::PD
    T_range::Tuple{Float64,Float64}
end

function GasSpecies(
    name::Symbol;
    M::Real,
    c,
    k,
    μ = nothing,
    V_D::Real = NaN,
    D = nothing,
    T_range::Tuple{<:Real,<:Real} = (200.0, 3000.0),
)
    isfinite(M) && M > 0 || throw(
        ArgumentError("GasSpecies :$name needs a positive molar mass M [kg/mol], got $M"),
    )
    isnan(V_D) ||
        V_D > 0 ||
        throw(
            ArgumentError("GasSpecies :$name diffusion volume must be positive, got $V_D"),
        )
    lo, hi = Float64(T_range[1]), Float64(T_range[2])
    isfinite(lo) && isfinite(hi) && 0 < lo < hi || throw(
        ArgumentError(
            "GasSpecies :$name T_range must be finite with 0 < T_lo < T_hi, got $T_range",
        ),
    )
    cp, kp = to_property(c), to_property(k)
    for (prop, what) in ((cp, "heat capacity"), (kp, "thermal conductivity"))
        for T in range(lo, hi, length = 61)
            v = prop(T)
            isfinite(v) && v > 0 || throw(
                ArgumentError(
                    "GasSpecies :$name: $what is $v at T = $T K; it must be positive on T_range = ($lo, $hi) K",
                ),
            )
        end
    end
    μp = μ === nothing ? nothing : to_property(μ)
    Dp = D === nothing ? nothing : to_property(D)
    return GasSpecies(name, Float64(M), cp, kp, μp, Float64(V_D), Dp, (lo, hi))
end

# The surface-temperature iterate held inside the species' property range.
@inline gas_property_temperature(gs::GasSpecies, T) =
    clamp(T, gs.T_range[1] + zero(T), gs.T_range[2] + zero(T))

"""
	nasa7_heat_capacity(a, M) -> PolynomialProperty
	nasa7_heat_capacity(a_low, a_high, M; T_common = 1000.0) -> CallableProperty

Heat capacity per unit mass [J/(kg·K)] from the NASA-7 polynomial coefficients
of `c_p/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴` (the first five of each set) and
the molar mass `M` [kg/mol]. With two sets, the low one applies below
`T_common` and the high one above it.
"""
function nasa7_heat_capacity(a::NTuple{5,<:Real}, M::Real)
    s = R_GAS / M
    return PolynomialProperty(ntuple(i -> Float64(a[i]) * s, 5))
end

function nasa7_heat_capacity(
    a_low::NTuple{5,<:Real},
    a_high::NTuple{5,<:Real},
    M::Real;
    T_common::Real = 1000.0,
)
    s = R_GAS / M
    lo = ntuple(i -> Float64(a_low[i]) * s, 5)
    hi = ntuple(i -> Float64(a_high[i]) * s, 5)
    Tc = Float64(T_common)
    return CallableProperty(T -> T < Tc ? evalpoly(T, lo) : evalpoly(T, hi))
end

"""
	eucken_viscosity(k, cp, M) -> μ

Eucken's relation between the conductivity and viscosity of a polyatomic gas,
`k = μ (c_p + 5R/(4M))`, solved for the viscosity [Pa·s]. `k` in W/(m·K),
`c_p` in J/(kg·K), `M` in kg/mol. Used where a species has no viscosity of
its own; the mixing rule for the outflow conductivity is the only consumer.
"""
@inline eucken_viscosity(k, cp, M) = k / (cp + 1.25 * R_GAS / M)

"""
	fuller_diffusivity(T, M_a, V_a, M_b, V_b; P = P_REF) -> D

Binary diffusion coefficient [m²/s] of gases `a` and `b` by the Fuller,
Schettler and Giddings correlation (Poling, Prausnitz and O'Connell, Eq.
11-4.4): `D = 0.00143 T^1.75 / (P √M_ab (V_a^⅓ + V_b^⅓)²)` in cm²/s with `P`
in bar, molar masses in g/mol and `M_ab = 2/(1/M_a + 1/M_b)`. Here `M` is in
kg/mol and `P` in Pa.
"""
@inline function fuller_diffusivity(T, M_a, V_a, M_b, V_b; P = P_REF)
    Ma, Mb = 1e3 * M_a, 1e3 * M_b
    M_ab = 2 / (1 / Ma + 1 / Mb)
    D_cm2 = 0.00143 * T^1.75 / ((P * 1e-5) * sqrt(M_ab) * (cbrt(V_a) + cbrt(V_b))^2)
    return 1e-4 * D_cm2
end

"""
	NITROGEN

Nitrogen as a free stream: NASA-7 heat capacity (Burcat's two sets, common
temperature 1000 K), conductivity and viscosity as 4th-order polynomials
fitted to NIST over 280–1250 K (the chamber CFD's nitrogen block), Fuller
diffusion volume 18.5. Property range 200–1250 K.
"""
const NITROGEN = GasSpecies(
    :N2;
    M = 0.028014,
    c = nasa7_heat_capacity(
        (3.53100528, -0.000123660988, -5.02999433e-07, 2.43530612e-09, -1.40881235e-12),
        (2.95257637, 1.39690040e-03, -4.92631603e-07, 7.86010195e-11, -4.60755204e-15),
        0.028014,
    ),
    k = PolynomialProperty((
        0.00025528103709,
        0.00010264086304,
        -6.7427534675e-08,
        3.9454196652e-11,
        -9.5578097215e-15,
    )),
    μ = PolynomialProperty((
        9.437486088e-07,
        6.9454771661e-08,
        -5.1584960289e-11,
        3.0001588211e-14,
        -7.2710794971e-18,
    )),
    V_D = 18.5,
    T_range = (200.0, 1250.0),
)

# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------

"""
	NoBlowing()

The coefficient is left alone whatever the outflow (`χ = 1`). The
`ConvectiveBC` default.
"""
struct NoBlowing <: AbstractBlowingModel end

"""
	CouetteFilm(λ = 0.5)

Film law with a fitted constant: `h = h₀ φ/(e^φ − 1)`,
`φ = 2λ ṁ″ c_p,g(T_s)/h₀`, with `ṁ″ c_p,g = Σ_g ṁ_g c_p,g(T_s)` the blowing
conductance of the gases leaving the face. `λ = 0.5` is the Couette film of
the classical theory (`φ = B = ṁ″ c_p,g/h₀`); a smaller `λ`, fitted to a
resolved flow for one gas, weakens the reduction at a given outflow and holds
only for that gas.
"""
struct CouetteFilm <: AbstractBlowingModel
    λ::Float64
    function CouetteFilm(λ::Real)
        isfinite(λ) && λ > 0 ||
            throw(ArgumentError("CouetteFilm: λ must be finite and > 0, got $λ"))
        return new(Float64(λ))
    end
end
CouetteFilm(; λ::Real = 0.5) = CouetteFilm(λ)

"""
	PropertyRatioFilm(free_stream; n_Γ = 0.6, volatiles = nothing)

Film law for any light gas, anchored on the Couette film of the free stream:

	h = h₀ φ/(e^φ − 1),  φ = (ṁ″ c_p,∞(T_s)/h₀) (Γ_∞/Γ_g)^n_Γ,  Γ = k/c_p

`c_p,∞`, `Γ_∞` are the free stream's at the surface temperature and `Γ_g` is
the leaving mixture's, its heat capacity weighted by mass fraction and its
conductivity by Wilke's rule over the species with outward flux. The gas
enters only through the property ratio, so no constant is refitted for a new
gas; `n_Γ = 0.6` is Kays's turbulent exponent on the heat-capacity ratio,
reproduced by the chamber CFD. With `n_Γ = 0` this is the Couette film in
the free-stream heat capacity.

`volatiles` replaces the per-species mixture by one `GasSpecies` (a surrogate
with its own `c_p`, `k`, `M`), which is how a model whose gas components are
lumped pseudo-species states the volatile properties it was built with.

The light-gas branch (`ψ = 1`) of [`PoolingFilm`](@ref).
"""
struct PropertyRatioFilm{G<:GasSpecies,V<:Union{Nothing,GasSpecies}} <: AbstractBlowingModel
    free_stream::G
    volatiles::V
    n_Γ::Float64
end

function PropertyRatioFilm(free_stream::GasSpecies; n_Γ::Real = 0.6, volatiles = nothing)
    isfinite(n_Γ) && n_Γ >= 0 ||
        throw(ArgumentError("PropertyRatioFilm: n_Γ must be finite and ≥ 0, got $n_Γ"))
    return PropertyRatioFilm(free_stream, volatiles, Float64(n_Γ))
end

"""
	PoolingFilm(free_stream; n_Γ = 0.6, n_D = 0.5, ψ0 = 0.59, Λc = 0.35, Λw = 0.10,
	            volatiles = nothing, diffusion_volume = (;))

The film law of [`PropertyRatioFilm`](@ref) with a pooling factor `ψ` for gas
that is denser than the free stream at the face:

	h = h₀ φ/(e^{φ/ψ} − 1)

`ψ` is the share of the no-outflow exchange the face keeps. It is 1 while the
plume lifts the gas off the face and falls to the floor `ψ0` when the gas
pools, switching over a logistic in the face buoyancy `Λ`:

	ψ = ψ0 + (1 − ψ0) / (1 + exp(−(Λ − Λc)/Λw))
	Λ = [1 − (M_face/T_s)(T_∞/M_∞)] / (1 − T_∞/T_s)

`M_face` is the molar mass of the gas at the face, from the face mass fractions
of the leaving species by the same film balance applied species by species:

	Y_i = ω_i [1 − exp(−ṁ″/(ψ g s_i^{n_D}))],  g = h₀/c̄_p,∞,  s_i = D_i,∞/α_∞

with `ω_i` the species' share of the outflow, `g` the gas-exchange conductance
without outflow (heat-transfer coefficient over the free-stream heat capacity
averaged between the free-stream and the surface temperature), `D_i,∞` the
binary diffusivity into the free stream, `α_∞` its thermal diffusivity and
`n_D = 1/2` the penetration-theory exponent. `ψ` and the face fractions are
solved together; where two states are self-consistent the mean of the two is
taken. A face no hotter than the free stream has no buoyancy to lose and keeps
`ψ = 1`.

The three switch constants belong to the apparatus flow, as `h₀` does; the
defaults are the CAPA chamber's (nitrogen purge, 185 SLPM). `diffusion_volume`
is a NamedTuple of Fuller diffusion volumes keyed by component name, for the
binary diffusivity of each gas component into the free stream
([`fuller_diffusivity`](@ref)); a component without an entry uses the
package's Chapman–Enskog estimate from its molar mass. Lumped `volatiles` use
their own `D`, else their diffusion volume, else the same estimate.
"""
struct PoolingFilm{G<:GasSpecies,V<:Union{Nothing,GasSpecies},N} <: AbstractBlowingModel
    free_stream::G
    volatiles::V
    n_Γ::Float64
    n_D::Float64
    ψ0::Float64
    Λc::Float64
    Λw::Float64
    diffusion_volume::NamedTuple{N,<:Tuple{Vararg{Float64}}}
end

function PoolingFilm(
    free_stream::GasSpecies;
    n_Γ::Real = 0.6,
    n_D::Real = 0.5,
    ψ0::Real = 0.59,
    Λc::Real = 0.35,
    Λw::Real = 0.10,
    volatiles = nothing,
    diffusion_volume = (;),
)
    isfinite(n_Γ) && n_Γ >= 0 ||
        throw(ArgumentError("PoolingFilm: n_Γ must be finite and ≥ 0, got $n_Γ"))
    isfinite(n_D) && n_D >= 0 ||
        throw(ArgumentError("PoolingFilm: n_D must be finite and ≥ 0, got $n_D"))
    isfinite(ψ0) && 0 < ψ0 < 1 ||
        throw(ArgumentError("PoolingFilm: ψ0 must lie in (0, 1), got $ψ0"))
    isfinite(Λc) || throw(ArgumentError("PoolingFilm: Λc must be finite, got $Λc"))
    isfinite(Λw) && Λw > 0 ||
        throw(ArgumentError("PoolingFilm: Λw must be finite and > 0, got $Λw"))
    all(v -> isfinite(v) && v > 0, values(diffusion_volume)) ||
        throw(ArgumentError("PoolingFilm: diffusion volumes must be finite and > 0"))
    Vnt = NamedTuple{keys(diffusion_volume)}(map(Float64, values(diffusion_volume)))
    return PoolingFilm(
        free_stream,
        volatiles,
        Float64(n_Γ),
        Float64(n_D),
        Float64(ψ0),
        Float64(Λc),
        Float64(Λw),
        Vnt,
    )
end

"""
	FunctionBlowing(f)

A user rule `χ = f(face, h, T_inf)`, where `face` is the
`BoundaryConditions.BlowingFace` at the current surface iterate (outward species
mass fluxes `face.ṁ`, surface temperature `face.T_s`, the `face.material`),
`h` the bare coefficient and `T_inf` the ambient temperature at this time.
`f` must return 1 at zero outflow, be smooth in its arguments and accept
`ForwardDiff.Dual` inputs (branch on `ForwardDiff.value` only).
"""
struct FunctionBlowing{F} <: AbstractBlowingModel
    f::F
end
