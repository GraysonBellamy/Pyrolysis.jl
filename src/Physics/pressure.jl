"""
Pressure field computation for gas phases.

Uses ideal gas law in a porous medium: P = Σⱼ (ξⱼ·R·T / Mⱼ) / φ for gas components,
where φ is the local porosity (gas volume fraction). The porosity divisor is
required; omitting it underestimates gas pressure at low porosity.
"""

using ..Materials: R_GAS, GAS, is_gas, molar_mass, component_index, n_phase_changes

const POROSITY_FLOOR = 1e-12

"""
    pore_state(ξ,T,material,raw = 1-skeletal_volume_fraction(ξ,T,material))

Central ideal-gas pore state from nonnegative property concentrations [kg/m³
bulk] and an admissible property temperature [K]. Raw porosity is one minus
condensed skeletal occupancy. The common C¹ divisor guard
is an identity for φ ≥ 2 POROSITY_FLOOR. Its use on trial states does not make
zero physical pore volume admissible. Vacuum has exactly zero pressures,
intrinsic densities and fractions, without an imposed background species.
"""
@inline function pore_state(
    ξ::NTuple{NC},
    T,
    material,
    raw = 1-skeletal_volume_fraction(ξ, T, material),
) where {NC}
    φ = smooth_floor(raw, POROSITY_FLOOR)
    ρ = map((c, x)->is_gas(c) ? x/φ : zero(x+T), material.components, ξ)
    p = map((c, r)->is_gas(c) ? r*R_GAS*T/molar_mass(c) : zero(r), material.components, ρ)
    P = sum(p)
    mass = sum(ρ)
    x = primal(P) == 0 ? map(zero, p) : map(v->v/P, p)
    Y = primal(mass) == 0 ? map(zero, ρ) : map(v->v/mass, ρ)
    return (; raw_porosity = raw, φ, ρ, p, P, x, Y)
end

"""Ideal-gas pressure [Pa] using the common guarded porosity divisor."""
@inline function ideal_gas_pressure(ξ::NTuple, T, φ, material)
    pore = smooth_floor(φ, POROSITY_FLOOR)
    return sum(
        map(
            (c, x)->is_gas(c) ? x*R_GAS*T/(molar_mass(c)*pore) : zero(x+T),
            material.components,
            ξ,
        ),
    )
end

"""
	background_gas_concentration(material, species, ξ, T; P = 101325.0) -> Float64

Bulk concentration [kg/m³] of gas component `species` whose partial pressure
equals `P` in the pore space of composition `ξ` at temperature `T` — the
single-species inverse of [`ideal_gas_pressure`](@ref):

	ξⱼ = P · Mⱼ · φ / (R · T),  φ = 1 − skeletal_volume_fraction(ξ, T)

Use this to initialize a background air component at ambient pressure in
gas-bearing problems: a pore space initialized to vacuum
(all gas ξⱼ = 0) against a 101325 Pa pressure boundary spuriously advects
early product gas *inward* down the absolute-pressure slope, instead of the
gas being resisted by ~1 atm of pore air.

`species` is a component index or name. The pore fraction is evaluated from
the condensed entries of `ξ` (gas entries carry no matrix volume), so passing
a solids-only initial composition is exact; a skeleton that leaves no pore
space is refused. Partial pressures of coexisting gas species
add — this function pressurizes `species` alone to `P`.

# Example
```julia
air = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026,)
# ... material includes `air`; ξ₀ is the condensed initial composition ...
ξ_air = background_gas_concentration(material, :AIR, ξ₀, 300.0)
```
"""
function background_gas_concentration(
    material::Material{NC},
    species::Symbol,
    ξ::NTuple{NC,<:Real},
    T::Real;
    P::Real = 101325.0,
) where {NC}
    return background_gas_concentration(
        material,
        component_index(material, species),
        ξ,
        T;
        P,
    )
end

function background_gas_concentration(
    material::Material{NC},
    j::Integer,
    ξ::NTuple{NC,<:Real},
    T::Real;
    P::Real = 101325.0,
) where {NC}
    c = material.components[j]
    is_gas(c) || throw(
        ArgumentError("background_gas_concentration: component $(c.name) is not gaseous"),
    )
    M = molar_mass(c)
    M > 0.0 || throw(
        ArgumentError(
            "background_gas_concentration: component $(c.name) has no molar mass",
        ),
    )
    isfinite(P) && P >= 0 ||
        throw(ArgumentError("background pressure must be finite and nonnegative"))
    isfinite(T) && T > 0 ||
        throw(ArgumentError("background temperature must be finite and positive"))
    φ = 1 - skeletal_volume_fraction(ξ, T, material)
    φ >= 2POROSITY_FLOOR ||
        throw(ArgumentError("background gas requires positive physical pore volume"))
    return P * M * φ / (R_GAS * T)
end

"""
	equilibrium_vapour_concentration(material, ξ, T; phase_change = 1) -> ξ_v [kg/m³]

Bulk concentration of the vapour in equilibrium with the liquid of the phase
change `phase_change` (an index into `material.phase_changes` or its name) in
the pore space of composition `ξ` at temperature `T`:

	ξ_v = φ M a(X, T) p_sat(T)/(R T),   φ = 1 − skeletal_volume_fraction(ξ, T),

with `X` the moisture content of `ξ` ([`moisture_content`](@ref)). Like
[`background_gas_concentration`](@ref), it reads the condensed entries of `ξ`
only, so a deck writes the vapour first and the background gas second, at the
pressure `P − p_eq` ([`equilibrium_vapour_pressure`](@ref)), for a pore gas in
equilibrium with the liquid at the start. A skeleton that leaves no pore space
is refused.

# Example
```julia
ξc = (ξ_wood, X₀ * ξ_wood, 0.0, 0.0)   # wood, water, vapour, air
ξ_v = equilibrium_vapour_concentration(material, ξc, 300.0)
p_v = equilibrium_vapour_pressure(material, ξc, 300.0)
ξ_air = background_gas_concentration(material, :AIR, ξc, 300.0; P = 101325.0 - p_v)
```
"""
function equilibrium_vapour_concentration(
    material::Material{NC},
    ξ::NTuple{NC,<:Real},
    T::Real;
    phase_change::Union{Integer,Symbol} = 1,
) where {NC}
    pc = _phase_change(material, phase_change)
    φ = _start_pore_fraction(material, ξ, T, "equilibrium_vapour_concentration")
    return equilibrium_vapour_concentration(pc, moisture_content(pc, ξ, ξ), T, φ)
end

"""
	equilibrium_vapour_pressure(material, ξ, T; phase_change = 1) -> p_eq [Pa]

Partial pressure of the vapour in equilibrium with the liquid of the phase
change `phase_change` (an index into `material.phase_changes` or its name) at
the composition `ξ` and temperature `T`, `a(X, T) p_sat(T)` with `X` the
moisture content of `ξ`; it reads the condensed entries of `ξ` only. The
partial pressure of the vapour of
[`equilibrium_vapour_concentration`](@ref).
"""
function equilibrium_vapour_pressure(
    material::Material{NC},
    ξ::NTuple{NC,<:Real},
    T::Real;
    phase_change::Union{Integer,Symbol} = 1,
) where {NC}
    pc = _phase_change(material, phase_change)
    (isfinite(T) && T > 0) || throw(
        ArgumentError(
            "equilibrium_vapour_pressure: the temperature must be finite and positive",
        ),
    )
    return equilibrium_vapour_pressure(pc, moisture_content(pc, ξ, ξ), T)
end

# The phase change of `material` named or indexed by `k`.
function _phase_change(material::Material, k::Integer)
    n = n_phase_changes(material)
    1 <= k <= n || throw(
        ArgumentError(
            "material $(material.name) has $n phase changes; phase_change = $k is " *
            "not one of them",
        ),
    )
    return material.phase_changes[k]
end

function _phase_change(material::Material, name::Symbol)
    k = findfirst(pc -> pc.name === name, material.phase_changes)
    k === nothing &&
        throw(ArgumentError("material $(material.name) has no phase change :$name"))
    return material.phase_changes[k]
end

# The pore fraction of the condensed entries of `ξ`, refused where it leaves
# no pore space.
function _start_pore_fraction(material, ξ, T, caller)
    (isfinite(T) && T > 0) ||
        throw(ArgumentError("$caller: the temperature must be finite and positive"))
    φ = 1 - skeletal_volume_fraction(ξ, T, material)
    φ >= 2POROSITY_FLOOR ||
        throw(ArgumentError("$caller: the composition leaves no pore volume"))
    return φ
end
