"""
Pressure field computation for gas phases.

Uses ideal gas law in a porous medium: P = Σⱼ (ξⱼ·R·T / Mⱼ) / φ for gas components,
where φ is the local porosity (gas volume fraction). The porosity divisor is
required; omitting it underestimates gas pressure at low porosity.
"""

using ..Materials: R_GAS, GAS, is_gas, molar_mass, component_index

# Floor on porosity used as a divisor — prevents 1/0 if a cell ever ends up
# fully solid. Matches `_compute_porosity` clamping at 0.0.
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

	ξⱼ = P · Mⱼ · φ(ξ, T) / (R · T)

Use this to initialize a background air component at ambient pressure in
DarcyFick-mode problems: a pore space initialized to vacuum
(all gas ξⱼ = 0) against a 101325 Pa pressure boundary spuriously advects
early product gas *inward* down the absolute-pressure slope, instead of the
gas being resisted by ~1 atm of pore air.

`species` is a component index or name. Porosity is evaluated from the
condensed entries of `ξ` (gas entries carry no matrix volume), so passing a
solids-only initial composition is exact. Partial pressures of coexisting gas
species add — this function pressurizes `species` alone to `P`.

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
    φ = porosity(ξ, T, material)
    return P * M * φ / (R_GAS * T)
end

"""
	compute_pressure!(P, φ, state, mesh, material)

Compute pressure field from ideal gas law for all cells, reading T and ξ from
the typed `state::StateCache`. `φ` is a precomputed porosity vector (same
length as `state.T`).
"""
function compute_pressure!(
    P::AbstractVector,
    φ::AbstractVector,
    state,
    mesh,
    material::Material{NC},
) where {NC}
    @inbounds for i = 1:mesh.n_active_cells
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        P[i] = ideal_gas_pressure(ξ_i, T_i, φ[i], material)
    end
    return nothing
end

# ============================================================================
# Pressure Gradient Computation
# ============================================================================

"""
	compute_pressure_gradient!(∇P::AbstractVector, P::AbstractVector, state, mesh)

Compute pressure gradient at cell faces using central differencing.

∇P_f = (P_R - P_L) / d_LR

# Arguments
- `∇P`: Output gradient vector [Pa/m], length = n_faces
- `P`: Cell-centered pressure values [Pa]
- `state`: Live geometry provider (`state.cell_centers` is read)
- `mesh`: Mesh topology (`mesh.faces` is read)

# Notes
- Boundary faces (left_id=0 or right_id=0) are set to zero; boundary 
  conditions should be applied separately.
- Gradient is positive when pressure increases in +z direction.
"""
function compute_pressure_gradient!(∇P::AbstractVector, P::AbstractVector, state, mesh)
    @inbounds for f = 1:(mesh.n_active_cells+1)
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            ∇P[f] = zero(eltype(∇P))
        else
            d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])
            ∇P[f] = (P[right_id] - P[left_id]) / d_LR
        end
    end
    return nothing
end
