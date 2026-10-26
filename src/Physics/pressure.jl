"""
Pressure field computation for gas phases.

Uses ideal gas law in a porous medium: P = Σⱼ (ξⱼ·R·T / Mⱼ) / φ for gas components,
where φ is the local porosity (gas volume fraction). The porosity divisor is
required; omitting it underestimates gas pressure at low porosity.
"""

using ..Materials: R_GAS, GAS, is_gas, molar_mass, component_index

const POROSITY_FLOOR = 1e-12

"""
    pore_state(ξ,T,material)

Central ideal-gas pore state from nonnegative property concentrations [kg/m³
bulk] and an admissible property temperature [K]. Raw porosity is one minus
condensed skeletal occupancy, independent of γ. The common C¹ divisor guard
is an identity for φ ≥ 2 POROSITY_FLOOR. Its use on trial states does not make
zero physical pore volume admissible. Vacuum has exactly zero pressures,
intrinsic densities and fractions, without an imposed background species.
"""
@inline function pore_state(ξ::NTuple{NC}, T, material) where {NC}
    raw = 1-skeletal_volume_fraction(ξ, T, material)
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

	ξⱼ = P · Mⱼ · φ / (R · T),  φ = max(1 − skeletal_volume_fraction(ξ, T), 0)

Use this to initialize a background air component at ambient pressure in
gas-bearing problems: a pore space initialized to vacuum
(all gas ξⱼ = 0) against a 101325 Pa pressure boundary spuriously advects
early product gas *inward* down the absolute-pressure slope, instead of the
gas being resisted by ~1 atm of pore air.

`species` is a component index or name. The pore fraction is evaluated from
the condensed entries of `ξ` (gas entries carry no matrix volume), so passing
a solids-only initial composition is exact; a skeleton that fills the cell
leaves no pore space and no gas. Partial pressures of coexisting gas species
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
