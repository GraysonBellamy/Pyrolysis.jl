"""
Pressure field computation for gas phases.

Uses ideal gas law in a porous medium: P = Σⱼ (ξⱼ·R·T / Mⱼ) / φ for gas components,
where φ is the local porosity (gas volume fraction). The porosity divisor is
required; omitting it underestimates gas pressure at low porosity.
"""

using ..Materials: R_GAS, GAS, is_gas, molar_mass, component_index

# Floor on porosity used as a divisor — prevents 1/0 if a cell ever ends up
# fully solid.
const _PHI_FLOOR = 1e-12

# ============================================================================
# Ideal Gas Pressure Computation
# ============================================================================

# Sum of partial pressures Σⱼ ξⱼ R T / Mⱼ over gas components only.
@generated function _gas_partial_pressure_sum(ξ, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        N = zero(T)
    end
    for i = 1:NC
        push!(expr.args, quote
            c = comps[$i]
            if is_gas(c) && molar_mass(c) > 0.0
                N += ξ[$i] * R_GAS * T / molar_mass(c)
            end
        end)
    end
    push!(expr.args, :(return N))
    return expr
end

"""
	ideal_gas_pressure(ξ::NTuple{NC,Float64}, T::Float64, φ::Float64, material) -> Float64

Compute ideal gas pressure in a porous medium from mass concentrations of gas
components and local porosity.

	P = Σⱼ (ξⱼ × R × T / Mⱼ) / φ   for gaseous components j

# Arguments
- `ξ`: Mass concentrations of all components [kg/m³]
- `T`: Temperature [K]
- `φ`: Porosity (gas volume fraction) [-], expected in (0, 1]
- `material`: Material definition with components

# Returns
- `P`: Pressure [Pa]

# Notes
Only components with `phase == GAS` and `molar_mass > 0` contribute to the partial
pressure sum. Solid and liquid components are ignored. `φ` is floored at $(_PHI_FLOOR)
to keep the divisor well-defined if a cell becomes fully solid.
"""
@inline function ideal_gas_pressure(
    ξ::NTuple{NC,<:Real},
    T::Real,
    φ::Real,
    material,
) where {NC}
    N = _gas_partial_pressure_sum(ξ, T, material.components, Val(NC))
    return N / max(φ, _PHI_FLOOR)
end

"""
	background_gas_concentration(material, species, ξ, T; P = 101325.0) -> Float64

Bulk concentration [kg/m³] of gas component `species` whose partial pressure
equals `P` in the pore space of composition `ξ` at temperature `T` — the
single-species inverse of [`ideal_gas_pressure`](@ref):

	ξⱼ = P · Mⱼ · φ / (R · T),  φ = max(1 − skeletal_volume_fraction(ξ, T), 0)

Use this to initialize a background air component at ambient pressure in
DarcyFick-mode problems: a pore space initialized to vacuum
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
air = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
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
    φ = max(1 - skeletal_volume_fraction(ξ, T, material), 0.0)
    return P * M * φ / (R_GAS * T)
end
