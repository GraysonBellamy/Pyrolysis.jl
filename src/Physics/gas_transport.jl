"""
The diffusive gas flux at a face.
Reference: ThermaKin2Ds Reference Manual, Section 1.2

Volume fraction gradient formulation:
	J_g = -ρ_gas × λ × ∂(ξ_g/ρ_gas)/∂x = -(λ/T) × ∂(ξT)/∂x

This is the volume fraction gradient formulation (equivalently partial pressure
gradient for ideal gas). It expands to J = -λ × [∂ξ/∂x + (ξ/T)×∂T/∂x].

At material interfaces, flux continuity requires the harmonic mean of λ,
which the gas-transport kernel forms at each face.
"""

"""
	gas_flux(λ_face::Float64, ξ_L::Float64, ξ_R::Float64,
			 T_L::Float64, T_R::Float64, d_LR::Float64) -> Float64

Compute gas mass flux at a face [kg/(m²·s)].

Uses the volume fraction gradient (partial pressure gradient) formulation:
	J = -(λ_face / T_face) × (ξ_R × T_R - ξ_L × T_L) / d_LR

where T_face = 0.5×(T_L + T_R) and λ_face is the harmonic mean of λ.
This is equivalent to J = -ρ_gas × λ × ∂(ξ/ρ_gas)/∂x, the volume fraction
gradient formulation.

# Validity (isobaric assumption)
`−(λ/T)∂(ξT)/∂z` equals Fick's law `−ρ_g·D·∇Y` (with `Y = ξ/(φρ_g)` and
`ρ_g = PM/RT`) **only when** `∇(φP) = 0`. The dropped term is proportional to
`−Y·∇(φP)/(φP)`. This is the inherited ThermaKin convention and is consistent
in the diffusion-only (Fickian) mode where the model carries no resolved
pressure field. In DarcyFick mode the code *does* resolve ∇P (and sharp char
fronts carry ∇φ), so this diffusive kernel is then a known approximation: the
pressure-driven part of the transport is carried by the Darcy advective term,
while the diffusive term keeps the isobaric form. Documented contradiction,
accepted by design; a full `∇(φP)` diffusive correction is future work.

# Arguments
- `λ_face`: Face transport coefficient (harmonic mean of λ) [m²/s]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `T_L, T_R`: Temperatures on left and right [K]
- `d_LR`: Distance between cell centers [m]
"""
@inline function gas_flux(
    λ_face::Real,
    ξ_L::Real,
    ξ_R::Real,
    T_L::Real,
    T_R::Real,
    d_LR::Real,
)
    # Volume fraction gradient: J = -(λ/T_face) × ∂(ξT)/∂x
    T_face = (T_L + T_R) / 2
    return -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
end
