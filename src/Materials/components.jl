"""
Material component definitions.

Components represent distinct phases (solid, liquid, gas) within a material,
each with their own thermophysical properties.
"""

# ============================================================================
# Phase Enumeration
# ============================================================================

"""
	Phase

Enumeration of material phases.

- `SOLID`: Solid phase (does not transport)
- `LIQUID`: Liquid phase (may transport)
- `GAS`: Gas phase (transports through porous matrix)
"""
@enum Phase begin
    SOLID = 1
    LIQUID = 2
    GAS = 3
end

# ============================================================================
# Component Types — phase becomes a *type-level* property
# ============================================================================
#
# Three concrete subtypes of `AbstractComponent`. Phase is implicit in the
# type, so `is_gas`, `is_solid`, `is_liquid`, `phase`, `molar_mass`, and
# `has_permeability` all dispatch statically with no runtime field branching.
#
# Each subtype carries only its phase-relevant fields:
# - `SolidComponent`/`LiquidComponent`: matrix-forming, with `permeability`
# - `GaseousComponent`: transports through the matrix, with `molar_mass`,
#    no permeability (gases are not part of the matrix)
#
# Cross-phase property access (e.g. asking a `SolidComponent` for its molar
# mass) goes through the dispatched accessors `molar_mass(c)`, `phase(c)`, etc.
# rather than direct field access — those accessors return phase-appropriate
# defaults (0.0, etc.) for components that don't carry the field.

"""
	SolidComponent{Pρ,Pc,Pk,Pλ,Pκ,Pε,Pa,Pρi} <: AbstractComponent

A solid-phase component with temperature-dependent thermophysical properties
and Darcy-flow permeability.

`density` is the *bulk* (pure-phase) density used by mixing rules, mixture
density, and volume-change kernels. `intrinsic_density` is the *skeletal*
(cell-wall) density used only by the porosity formula
`φ = 1 − Σ ξⱼ/ρⱼ_intrinsic`. When not specified, `intrinsic_density` defaults
to `density`, recovering the implicit assumption that the pure
phase has no pores.

`emissivity` is a property function of temperature, evaluated at the cell
temperature by the mixing rule (`ε_eff = Σⱼ vⱼ εⱼ(T)`). `surface_absorptivity`
is an optional constant absorptivity for the *incident* (heater-spectrum)
radiation: `nothing` means Kirchhoff (`aⱼ = εⱼ(T)`), so a material whose
components all leave it unset absorbs at exactly `ε_eff`.
"""
struct SolidComponent{
    Pρ<:AbstractPropertyFunction,
    Pc<:AbstractPropertyFunction,
    Pk<:AbstractPropertyFunction,
    Pλ<:AbstractPropertyFunction,
    Pκ<:AbstractPropertyFunction,
    Pε<:AbstractPropertyFunction,
    Pa<:Union{Nothing,Float64},
    Pρi<:AbstractPropertyFunction,
} <: AbstractComponent
    name::Symbol
    density::Pρ
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    permeability::Pκ
    emissivity::Pε
    surface_absorptivity::Pa
    absorption_coeff::Float64
    swelling_factor::Float64
    intrinsic_density::Pρi
end

"""
	GaseousComponent{Pρ,Pc,Pk,Pλ,Pε} <: AbstractComponent

A gas-phase component. Density is typically the ideal-gas law. Gases flow
*through* the matrix, not part of it — there is no `permeability` field, and
no `surface_absorptivity` (a gas always follows Kirchhoff, `aⱼ = εⱼ(T)`).
"""
struct GaseousComponent{
    Pρ<:AbstractPropertyFunction,
    Pc<:AbstractPropertyFunction,
    Pk<:AbstractPropertyFunction,
    Pλ<:AbstractPropertyFunction,
    Pε<:AbstractPropertyFunction,
} <: AbstractComponent
    name::Symbol
    molar_mass::Float64           # kg/mol
    density::Pρ
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    emissivity::Pε
    absorption_coeff::Float64
    swelling_factor::Float64
end

"""
	LiquidComponent{Pρ,Pc,Pk,Pλ,Pκ,Pε,Pa,Pρi} <: AbstractComponent

A liquid-phase component (matrix-forming, with permeability).

See [`SolidComponent`](@ref) for the bulk vs. intrinsic density distinction
and for the `emissivity` / `surface_absorptivity` semantics.
"""
struct LiquidComponent{
    Pρ<:AbstractPropertyFunction,
    Pc<:AbstractPropertyFunction,
    Pk<:AbstractPropertyFunction,
    Pλ<:AbstractPropertyFunction,
    Pκ<:AbstractPropertyFunction,
    Pε<:AbstractPropertyFunction,
    Pa<:Union{Nothing,Float64},
    Pρi<:AbstractPropertyFunction,
} <: AbstractComponent
    name::Symbol
    density::Pρ
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    permeability::Pκ
    emissivity::Pε
    surface_absorptivity::Pa
    absorption_coeff::Float64
    swelling_factor::Float64
    intrinsic_density::Pρi
end

# ============================================================================
# Convenience Constructors
# ============================================================================

"""
	SolidComponent(name; ρ, c, k, ε=0.9, absorptivity=nothing, α=0.0, γ=1.0,
	               λ=0.0, κ=0.0, ρ_intrinsic=nothing, φ=nothing)

Create a solid component with simplified interface.

# Arguments
- `name::Symbol`: Component identifier
- `ρ`: Bulk (pure-phase) density [kg/m³] (number or property function)
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `ε=0.9`: Surface emissivity [-] (number or property function of `T`, same
  forms as `k`). Evaluated at the cell temperature in the mixing rule
  `ε_eff = Σⱼ vⱼ εⱼ(T)`, so the surface-emission and Beer-Lambert Jacobians
  carry `dε/dT`.
- `absorptivity=nothing`: Constant surface absorptivity [-] for *incident*
  radiation (a property of the heater spectrum, not of the surface
  temperature). `nothing` = Kirchhoff, `aⱼ = εⱼ(T)`. A `RadiativeFluxBC`
  without an explicit `absorptivity` absorbs at the mixture
  `a_eff = Σⱼ vⱼ aⱼ`, which equals `ε_eff` when no component sets this.
- `α=0.0`: Mass absorption coefficient [m²/kg] (Beer-Lambert in-depth
  absorption, unrelated to `absorptivity`)
- `γ=1.0`: Swelling factor [-] (default 1.0: solid contributes to volume change)
- `λ=0.0`: Gas transfer coefficient [m²/s] for mixture property calculation
- `κ=0.0`: Permeability [m²] for Darcy flow (0 = impermeable)
- `ρ_intrinsic=nothing`: Skeletal/intrinsic density of the cell-wall material
  [kg/m³]. Used **only** by the porosity formula `φ = 1 − Σ ξⱼ/ρⱼ_intrinsic`,
  not by mixing rules, mixture density, or volume-change kernels (those keep
  using bulk `ρ`). When omitted, defaults to `ρ` — i.e. the pure phase has no
  pores, the implicit assumption.
- `φ=nothing`: Alternative to `ρ_intrinsic`. Porosity of the pure phase
  (dimensionless, 0 ≤ φ < 1). When given with constant bulk `ρ`, sets
  `ρ_intrinsic = ρ / (1 − φ)`. Pass either `ρ_intrinsic` *or* `φ`, not both.

# Note on swelling_factor
γ is the fraction of this component's intrinsic volume that contributes to
the bulk (matrix) volume. It has **one** meaning everywhere: the dilation
rate θ = Σγⱼ(dξⱼ/dt)/ρⱼ, the mixture density 1/ρ = Σγⱼ·Yⱼ/ρⱼ, and the
conductivity/emissivity volume weights (`solid_volume_fractions`) all use
the same γⱼ·ξⱼ/ρⱼ weighting.
- γ=1.0: Full contribution (default for solids — when solid depletes, volume shrinks)
- γ=0.0: Consistently *volumeless* — no strain on conversion, no volume in
  the bulk density, **and no weight in conductivity/emissivity mixing**. Use
  with care: a γ=0 "rigid char" stops contributing its k and ε to the mixture.
- γ∈(0,1): Partial contribution (advanced use cases)

# Note on gas_transfer_coeff (λ)
Although solids don't transport themselves, they contribute to the effective gas
transfer coefficient of the mixture through volume-weighted mixing rules. Set to
a non-zero value (e.g., 2e-5 m²/s) to match ThermaKin conventions. With the
default `λ = 0`, Fickian gas transport vanishes in a cell of this solid that
holds no gas, and the gas field of a gas-producing material ends in a sharp
front on which the implicit Newton fails repeatedly; problem validation warns
in that configuration (`warn_degenerate_gas_transport`).

# Note on permeability (κ)
Permeability determines resistance to pressure-driven (Darcy) gas flow through
the solid matrix. Typical values:
- Virgin wood: ~10⁻¹⁵ m² (tight cellular structure)
- Char: ~10⁻¹² m² (porous with open channels)
- Impermeable: 0.0 m² (default)

# Example
```julia
virgin = SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86)
char = SolidComponent(:char, ρ=200, c=1000, k=0.10, ε=0.95, α=5.0, κ=1e-12)
# Char that doesn't contribute to volume change:
rigid_char = SolidComponent(:char, ρ=200, c=1000, k=0.10, γ=0.0)
# Emissivity falling with temperature, heater absorptivity fixed by the lamp spectrum:
wood = SolidComponent(:wood, ρ=500, c=1500, k=0.12, ε=(0.95, -1.8e-4), absorptivity=0.76)
```
"""
function SolidComponent(
    name::Symbol;
    ρ,
    c,
    k,
    ε = 0.9,
    absorptivity::Union{Nothing,Real} = nothing,
    α::Real = 0.0,
    γ::Real = 1.0,
    λ::Real = 0.0,
    κ::Real = 0.0,
    ρ_intrinsic = nothing,
    φ = nothing,
    validate_properties::Bool = true,
)
    _validate_swelling_factor(γ, name)
    bulk_ρ = to_property(ρ)
    intrinsic_ρ = _resolve_intrinsic_density(bulk_ρ, ρ, ρ_intrinsic, φ, name)
    comp = SolidComponent(
        name,
        bulk_ρ,
        to_property(c),
        to_property(k),
        to_property(λ),  # contributes to effective gas transfer coefficient
        to_property(κ),  # permeability for Darcy flow
        to_property(ε),
        _resolve_surface_absorptivity(absorptivity, name),
        Float64(α),
        Float64(γ),
        intrinsic_ρ,
    )
    validate_properties && _validate_condensed_properties(comp)
    return comp
end

# ============================================================================
# Property-law validation on the evaluation bracket
# ============================================================================
#
# Property laws are evaluated on [T_PROPERTY_MIN, T_PROPERTY_MAX] (their
# argument is clamped to it), so a law that is well behaved there is well
# behaved for every temperature the solver can hand it. The constructors
# sample each law on that bracket at 1 K steps: ρ, c and k must be positive
# (k ≤ 0 makes a cell a perfect insulator and stalls the integrator; ρ or
# c ≤ 0 inverts the energy equation), λ and κ non-negative, ε within [0, 1];
# a gas may have c = 0 (no sensible-enthalpy transport).
# `validate_properties = false` skips the check for laws built
# programmatically from search parameters, where an unphysical value should
# reach the objective function instead of throwing.

"""
	_validate_property_law(p, name, what, unit; lower, upper, strict_lower)

Sample the property law `p` at 1 K steps over `[T_PROPERTY_MIN, T_PROPERTY_MAX]`
and throw an `ArgumentError` naming the first temperature at which it leaves
`(lower, upper]` (or `[lower, upper]` when `strict_lower = false`) or is not
finite. `StateDependentProperty` laws need a composition and are not sampled.
"""
function _validate_property_law(
    p::AbstractPropertyFunction,
    name::Symbol,
    what::String,
    unit::String;
    lower::Float64 = 0.0,
    upper::Float64 = Inf,
    strict_lower::Bool = true,
)
    p isa StateDependentProperty && return nothing
    for T = T_PROPERTY_MIN:1.0:T_PROPERTY_MAX
        v = Float64(p(T))
        ok = isfinite(v) && (strict_lower ? v > lower : v >= lower) && v <= upper
        ok && continue
        bound =
            upper == Inf ? (strict_lower ? "positive" : "non-negative") :
            "within [$(lower), $(upper)]"
        throw(
            ArgumentError(
                "component :$name: $what is $v $unit at T = $T K; it must be $bound on " *
                "[$(T_PROPERTY_MIN), $(T_PROPERTY_MAX)] K, the bracket property laws are " *
                "evaluated on. Hold the fit constant outside its measured range " *
                "(`T -> law(clamp(T, T_lo, T_hi))`) or pass `validate_properties = false`.",
            ),
        )
    end
    return nothing
end

function _validate_condensed_properties(c::Union{SolidComponent,LiquidComponent})
    n = c.name
    _validate_property_law(c.density, n, "bulk density", "kg/m³")
    _validate_property_law(c.intrinsic_density, n, "intrinsic density", "kg/m³")
    _validate_property_law(c.heat_capacity, n, "heat capacity", "J/(kg·K)")
    _validate_property_law(c.thermal_conductivity, n, "thermal conductivity", "W/(m·K)")
    _validate_property_law(
        c.gas_transfer_coeff,
        n,
        "gas transfer coefficient",
        "m²/s";
        strict_lower = false,
    )
    _validate_property_law(c.permeability, n, "permeability", "m²"; strict_lower = false)
    _validate_property_law(
        c.emissivity,
        n,
        "emissivity",
        "";
        upper = 1.0,
        strict_lower = false,
    )
    return nothing
end

function _validate_gas_properties(c::GaseousComponent)
    n = c.name
    # A gas with c = 0 carries no sensible enthalpy (S_conv vanishes); the
    # matrix heat capacity floor keeps the energy equation
    # well posed, so zero is admissible here where it is not for a solid.
    _validate_property_law(
        c.heat_capacity,
        n,
        "heat capacity",
        "J/(kg·K)";
        strict_lower = false,
    )
    _validate_property_law(c.thermal_conductivity, n, "thermal conductivity", "W/(m·K)")
    _validate_property_law(
        c.gas_transfer_coeff,
        n,
        "gas transfer coefficient",
        "m²/s";
        strict_lower = false,
    )
    _validate_property_law(
        c.emissivity,
        n,
        "emissivity",
        "";
        upper = 1.0,
        strict_lower = false,
    )
    return nothing
end

# Surface absorptivity is a fraction of the incident flux; `nothing` keeps the
# Kirchhoff default (aⱼ = εⱼ(T)) at the type level so the mixing sum is static.
@inline _resolve_surface_absorptivity(::Nothing, ::Symbol) = nothing
@inline function _resolve_surface_absorptivity(a::Real, name::Symbol)
    0.0 <= a <= 1.0 || throw(
        ArgumentError(
            "component :$name has surface absorptivity $a outside [0, 1]. It is the " *
            "fraction of the incident radiation absorbed at the surface; leave it " *
            "unset (nothing) to use Kirchhoff's law a = ε(T).",
        ),
    )
    return Float64(a)
end

# γ is a volume fraction weight; values outside [0,1] have no physical
# meaning under the universal swelling semantics (θ, mixture density, and
# conductivity/emissivity weights all share it).
@inline function _validate_swelling_factor(γ::Real, name::Symbol)
    0.0 <= γ <= 1.0 || throw(
        ArgumentError(
            "component :$name has swelling factor γ = $γ outside [0, 1]. γ is " *
            "the fraction of the component's intrinsic volume contributing to " *
            "the bulk volume (θ, mixture density, and conductivity/emissivity " *
            "mixing all share this weighting).",
        ),
    )
    return nothing
end

"""
	GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0, γ=0.0)

Create a gaseous component with simplified interface.

For ideal gases, density is computed from the ideal gas law **at the
reference pressure**: ρ_g(T) = P_REF·M/(R·T). Property functions carry no
pressure argument, so in pressurized cells the gas volume fraction used in
λ/κ mixing (and the γ > 0 swelling-gas mixture density) is over-weighted by
~P/P_REF. The pressure closure P = N/φ and porosity are unaffected.

Gases do NOT have permeability - they flow THROUGH the solid matrix.

# Arguments
- `name::Symbol`: Component identifier
- `M`: Molar mass [kg/mol]
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `λ`: Gas transfer coefficient [m²/s] (number or property function)
- `ε=0.0`: Surface emissivity [-] (number or property function of `T`;
  typically 0 for gases)
- `α=0.0`: Mass absorption coefficient [m²/kg]
- `γ=0.0`: Swelling factor [-] (default 0.0: gas escapes, doesn't affect matrix)

# Note on swelling_factor
The swelling factor controls whether this component contributes to volume change:
- γ=0.0: No contribution (default for gases - gas escapes through pores)
- γ=1.0: Full contribution (for intumescent materials where gas creates new volume)

# Note on permeability
Gases do not have permeability because they flow through the solid matrix.
Only solid and liquid components contribute to the effective permeability.

# Example
```julia
MMA = GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
water_vapor = GaseousComponent(:H2O, M=0.018, c=2000, k=0.025, λ=2e-5)
# Intumescent gas that creates new volume:
intumescent_gas = GaseousComponent(:gas, M=0.040, c=1000, k=0.02, λ=1e-5, γ=1.0)
```
"""
function GaseousComponent(
    name::Symbol;
    M::Real,
    c,
    k,
    λ,
    ε = 0.0,
    α::Real = 0.0,
    γ::Real = 0.0,
    validate_properties::Bool = true,
)
    # A gas with M ≤ 0 would silently exert no pressure: the ideal-gas closure
    # P = Σ ξⱼRT/(Mⱼφ) filters on molar_mass > 0, while transport does not —
    # the species would advect/diffuse but never contribute to P.
    M > 0 || throw(
        ArgumentError(
            "GaseousComponent :$name requires a positive molar mass M [kg/mol], " *
            "got $M. A gas with M ≤ 0 is excluded from the pressure closure " *
            "P = Σ ξⱼRT/(Mⱼφ) but still transported — silently pressure-less.",
        ),
    )
    _validate_swelling_factor(γ, name)
    # Ideal gas density: ρ = PM/(RT), evaluated at the REFERENCE pressure.
    # Property functions are T-only, so every consumer of this density
    # (volume_fractions → λ/κ mixing; γ > 0 swelling-gas mixture density)
    # sees ρ_g(T, P_REF): in a pressurized cell the gas volume is
    # over-weighted by ~P/P_REF. The pressure closure
    # P = N/φ and the porosity formula are unaffected — they never use this
    # density. Acceptable at the ≤ atm-scale overpressures of validated
    # cases; revisit if property functions ever gain a pressure argument.
    M_float = Float64(M)
    ρ_func = CallableProperty(T -> P_REF * M_float / (R_GAS * T))

    comp = GaseousComponent(
        name,
        M_float,
        ρ_func,
        to_property(c),
        to_property(k),
        to_property(λ),
        to_property(ε),
        Float64(α),
        Float64(γ),
    )
    validate_properties && _validate_gas_properties(comp)
    return comp
end

"""
	LiquidComponent(name; ρ, c, k, λ=2e-5, κ=0.0, ε=0.95, absorptivity=nothing,
	                α=0.0, γ=0.0)

Create a liquid component with simplified interface.

# Arguments
- `name::Symbol`: Component identifier
- `ρ`: Density [kg/m³] (number or property function)
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `λ=2e-5`: Transport coefficient [m²/s] for mixture property calculation
- `κ=0.0`: Permeability [m²] for Darcy flow (0 = blocks flow, typical for liquids)
- `ε=0.95`: Surface emissivity [-] (number or property function of `T`).
  Condensed liquids are strong thermal emitters (water ≈ 0.96), so a low
  value on a γ > 0 surface liquid would *lower* ε_eff instead of raising it.
  Irrelevant for γ = 0 liquids, which carry no weight in emissivity mixing.
- `absorptivity=nothing`: Constant surface absorptivity [-] for incident
  radiation; `nothing` = Kirchhoff `aⱼ = εⱼ(T)`. See [`SolidComponent`](@ref).
- `α=0.0`: Mass absorption coefficient [m²/kg]
- `γ=0.0`: Swelling factor [-]

# Note on swelling_factor
The default γ = 0 makes a liquid **volumeless** under the universal swelling
semantics: it contributes no strain to θ, no volume to the bulk density, and
no weight in conductivity/emissivity mixing — ThermaKin's `L SWELLING: 0`
behavior for moisture bound in cell walls (this makes the historical
pseudo-density hack ρ ≈ 1e6 redundant, though harmless). Set γ = 1 for a
liquid that physically occupies bulk volume (e.g. a melt layer).

# Note on gas_transfer_coeff (λ)
Like solids, liquids contribute to the effective gas transfer coefficient of the
mixture through volume-weighted mixing rules. The default value (2e-5 m²/s)
matches ThermaKin conventions.

# Note on permeability (κ)
Liquids typically block gas flow (κ=0), as the liquid fills pore space.
Set to non-zero only for special cases like partially saturated porous media.
"""
function LiquidComponent(
    name::Symbol;
    ρ,
    c,
    k,
    λ = 2e-5,
    κ::Real = 0.0,
    ε = 0.95,
    absorptivity::Union{Nothing,Real} = nothing,
    α::Real = 0.0,
    γ::Real = 0.0,
    ρ_intrinsic = nothing,
    φ = nothing,
    validate_properties::Bool = true,
)
    _validate_swelling_factor(γ, name)
    bulk_ρ = to_property(ρ)
    intrinsic_ρ = _resolve_intrinsic_density(bulk_ρ, ρ, ρ_intrinsic, φ, name)
    comp = LiquidComponent(
        name,
        bulk_ρ,
        to_property(c),
        to_property(k),
        to_property(λ),
        to_property(κ),  # permeability for Darcy flow
        to_property(ε),
        _resolve_surface_absorptivity(absorptivity, name),
        Float64(α),
        Float64(γ),
        intrinsic_ρ,
    )
    validate_properties && _validate_condensed_properties(comp)
    return comp
end

# ============================================================================
# Intrinsic-density resolution
# ============================================================================

"""
	_resolve_intrinsic_density(bulk_ρ, raw_ρ, ρ_intrinsic, φ, name)

Pick the right intrinsic-density property function for a matrix-forming
component, given the user-supplied kwargs. Rules:

1. If both `ρ_intrinsic` and `φ` are given → error.
2. If `ρ_intrinsic` is given → wrap via `to_property`.
3. If `φ` is given → require `raw_ρ::Real` (porosity-from-bulk only well-defined
   for constant bulk density) and return `ConstantProperty(ρ / (1 − φ))`.
4. Else → fall back to the bulk density itself (no porosity, the historical
   default that keeps every existing example unchanged).
"""
@inline function _resolve_intrinsic_density(
    bulk_ρ::AbstractPropertyFunction,
    raw_ρ,
    ρ_intrinsic,
    φ,
    name::Symbol,
)
    if ρ_intrinsic !== nothing && φ !== nothing
        throw(ArgumentError("$name: pass either ρ_intrinsic or φ, not both"))
    elseif ρ_intrinsic !== nothing
        return to_property(ρ_intrinsic)
    elseif φ !== nothing
        raw_ρ isa Real || throw(
            ArgumentError(
                "$name: φ kwarg requires a constant `ρ`; for T-dependent bulk density, " *
                "supply ρ_intrinsic directly",
            ),
        )
        0.0 ≤ Float64(φ) < 1.0 ||
            throw(ArgumentError("$name: φ must satisfy 0 ≤ φ < 1, got $(φ)"))
        return ConstantProperty(Float64(raw_ρ) / (1.0 - Float64(φ)))
    else
        return bulk_ρ
    end
end

# ============================================================================
# Property Accessors
# ============================================================================

"""
	get_intrinsic_density(comp::AbstractComponent, T::Real) -> Real

Skeletal/intrinsic density at temperature T — the density of the cell-wall
material itself, ignoring pore volume. Used by the porosity formula
`φ = 1 − Σ ξⱼ/ρⱼ_intrinsic`.

For matrix-forming components (`SolidComponent`, `LiquidComponent`) this
returns the stored `intrinsic_density` field, which defaults to the bulk
`density` when the component was constructed without an explicit
`ρ_intrinsic` or `φ`. For `GaseousComponent` returns 0.0 (gases occupy the
pore space, they are not part of the matrix).
"""
@inline get_intrinsic_density(comp::SolidComponent, T::Real) = comp.intrinsic_density(T)
@inline get_intrinsic_density(comp::LiquidComponent, T::Real) = comp.intrinsic_density(T)
@inline get_intrinsic_density(::GaseousComponent, ::Real) = 0.0

"""
	get_permeability(comp::AbstractComponent, T::Real) -> Float64

Get component permeability at temperature T. `GaseousComponent` returns 0.0
(gases are not part of the matrix).
"""
@inline get_permeability(comp::SolidComponent, T::Real) = comp.permeability(T)
@inline get_permeability(comp::LiquidComponent, T::Real) = comp.permeability(T)
@inline get_permeability(::GaseousComponent, ::Real) = 0.0

"""
	get_emissivity(comp::AbstractComponent, T::Real) -> Real

Component surface emissivity at temperature `T`. AD-clean: a
`ForwardDiff.Dual` `T` returns a Dual carrying `dε/dT`.
"""
@inline get_emissivity(comp::AbstractComponent, T::Real) = comp.emissivity(T)

"""
	surface_absorptivity(comp::AbstractComponent) -> Union{Nothing,Float64}

The component's constant surface absorptivity for incident radiation, or
`nothing` when it follows Kirchhoff's law (`aⱼ = εⱼ(T)`). Always `nothing`
for a `GaseousComponent`.
"""
@inline surface_absorptivity(comp::SolidComponent) = comp.surface_absorptivity
@inline surface_absorptivity(comp::LiquidComponent) = comp.surface_absorptivity
@inline surface_absorptivity(::GaseousComponent) = nothing

"""
	get_absorptivity(comp::AbstractComponent, T::Real) -> Real

Component surface absorptivity for incident radiation at temperature `T`:
the constant `surface_absorptivity` when set, otherwise the Kirchhoff value
`get_emissivity(comp, T)`. The branch resolves at compile time (the field's
`Nothing`/`Float64` type is a struct parameter), so a material with no
absorptivity set mixes to exactly `ε_eff`.
"""
@inline get_absorptivity(comp::AbstractComponent, T::Real) =
    _absorptivity_or_kirchhoff(surface_absorptivity(comp), comp, T)
@inline _absorptivity_or_kirchhoff(::Nothing, comp, T) = get_emissivity(comp, T)
@inline _absorptivity_or_kirchhoff(a::Float64, comp, T) = a

"""
	has_permeability(comp::AbstractComponent) -> Bool

True for matrix-forming components (solids and liquids). False for gases.
"""
@inline has_permeability(::SolidComponent) = true
@inline has_permeability(::LiquidComponent) = true
@inline has_permeability(::GaseousComponent) = false

# ----- Phase predicates: static dispatch on type, no runtime branching -----

"""
	is_solid(comp::AbstractComponent) -> Bool

`true` for `SolidComponent`. Resolves at compile time via dispatch.
"""
@inline is_solid(::SolidComponent) = true
@inline is_solid(::AbstractComponent) = false

"""
	is_liquid(comp::AbstractComponent) -> Bool

`true` for `LiquidComponent`. Resolves at compile time via dispatch.
"""
@inline is_liquid(::LiquidComponent) = true
@inline is_liquid(::AbstractComponent) = false

"""
	is_gas(comp::AbstractComponent) -> Bool

`true` for `GaseousComponent`. Resolves at compile time via dispatch.
"""
@inline is_gas(::GaseousComponent) = true
@inline is_gas(::AbstractComponent) = false

"""
	phase(comp::AbstractComponent) -> Phase

Return the runtime `Phase` enum value for a component, derived from its type.
Use this when a `Phase` value is needed for printing, logging, or comparison;
use `is_gas`/`is_solid`/`is_liquid` for predicates in hot paths.
"""
@inline phase(::SolidComponent) = SOLID
@inline phase(::LiquidComponent) = LIQUID
@inline phase(::GaseousComponent) = GAS

"""
	molar_mass(comp::AbstractComponent) -> Float64

Molar mass [kg/mol]. Returns the stored value for gases; 0.0 for solids and
liquids (which do not carry a molar-mass field).
"""
@inline molar_mass(comp::GaseousComponent) = comp.molar_mass
@inline molar_mass(::AbstractComponent) = 0.0

# ============================================================================
# Radiation Property Queries
# ============================================================================

"""
	has_absorption_coefficients(components::Tuple) -> Bool

Returns true if any solid component has a non-zero absorption coefficient.

This is used to warn users when SURFACE_ABSORPTION radiation model is selected
but material components have α values set (which will be ignored).

# Arguments
- `components`: Tuple of Component objects from a Material

# Returns
- `true` if any component has absorption_coeff > 0
- `false` otherwise

# Example
```julia
if has_absorption_coefficients(material.components)
	@warn "Absorption coefficients will be ignored with SURFACE_ABSORPTION model"
end
```
"""
function has_absorption_coefficients(components::Tuple)
    for comp in components
        if comp.absorption_coeff > 0.0
            return true
        end
    end
    return false
end
