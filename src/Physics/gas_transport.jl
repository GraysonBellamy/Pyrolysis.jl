"""
Gas transport (Darcy-Fick).
Reference: ThermaKin2Ds Reference Manual, Section 1.2

Volume fraction gradient formulation:
	J_g = -ρ_gas × λ × ∂(ξ_g/ρ_gas)/∂x = -(λ/T) × ∂(ξT)/∂x

This is the volume fraction gradient formulation (equivalently partial pressure
gradient for ideal gas). It expands to J = -λ × [∂ξ/∂x + (ξ/T)×∂T/∂x].

At material interfaces, flux continuity requires the harmonic mean of λ.
This is computed as λ_face in the property cache.
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

# ============================================================================
# Darcy-Fick Combined Advection-Diffusion Flux
# ============================================================================

"""
	gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR) -> flux

Compute gas mass flux at a face using full Darcy-Fick formulation [kg/(m²·s)].

Combines advective and diffusive transport:
	J = (ξ_upwind / φ_face) × v_g - (λ_face/T_face) × (ξ_R×T_R - ξ_L×T_L) / d_LR
		↑ advection                  ↑ diffusion (volume fraction gradient)

ξ is mass per **bulk** volume while the Darcy velocity `v_g = -(κ/μ)∇P` is
**superficial** (volumetric flux per total cross-section). Mass flux per bulk
area is intrinsic gas density × superficial velocity, J_adv = (ξ/φ)·v_g —
equivalently ξ · v_interstitial with v_i = v_g/φ. The face porosity is the
arithmetic mean φ_face = (φ_L + φ_R)/2, floored at `_PHI_FLOOR` consistently
with the ideal-gas pressure closure P = Σξ·R·T/(M·φ).

# Arguments
- `v_g`: superficial Darcy gas velocity at face [m/s] (positive = +z)
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `φ_L, φ_R`: Porosity (gas volume fraction) on left and right [-]
- `λ_face`: Face transport coefficient (harmonic mean of λ) [m²/s]
- `T_L, T_R`: Temperatures on left and right [K]
- `d_LR`: Distance between cell centers [m]

# Returns
Total flux [kg/(m²·s)], positive = flow in +z direction
"""
@inline function gas_flux_darcy_fick(
    v_g::Real,
    ξ_L::Real,
    ξ_R::Real,
    φ_L::Real,
    φ_R::Real,
    λ_face::Real,
    T_L::Real,
    T_R::Real,
    d_LR::Real,
)
    # Diffusive flux (volume fraction / partial pressure gradient)
    T_face = (T_L + T_R) / 2
    J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)

    # Advective flux: intrinsic density (ξ/φ_face) advected by the superficial
    # velocity, with first-order upwinding of ξ
    φ_face = max((φ_L + φ_R) / 2, _PHI_FLOOR)
    ξ_upwind = v_g >= 0 ? ξ_L : ξ_R
    J_adv = v_g * ξ_upwind / φ_face

    return J_adv + J_diff
end

"""
	compute_gas_fluxes!(J::Matrix{Float64}, mesh::Mesh1D{NC},
						cache::PropertyCache, material::Material{NC}) where NC

Compute gas mass fluxes at all faces for all gas components.

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache
- `material`: Material definition

Positive flux = flow in +z direction (toward bottom).

# Performance Note
This function uses compile-time unrolling via `compute_gas_fluxes_unrolled!`
to avoid type instability from heterogeneous tuple indexing.
"""
function compute_gas_fluxes!(
    J::AbstractMatrix,
    state,
    mesh,
    cache::PropertyCache,
    material::Material{NC},
) where {NC}
    n_faces_val = mesh.n_active_cells + 1

    fill!(J, zero(eltype(J)))

    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            continue
        end

        d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])
        λ_face = cache.λ_face[f]
        T_L = state.T[left_id]
        T_R = state.T[right_id]

        ξ_L = ntuple(j -> state.ξ[j, left_id], Val(NC))
        ξ_R = ntuple(j -> state.ξ[j, right_id], Val(NC))

        compute_gas_fluxes_unrolled!(J, f, ξ_L, ξ_R, λ_face, T_L, T_R, d_LR, material)
    end

    return nothing
end

"""
	compute_gas_fluxes_darcy_fick!(J::Matrix{Float64}, state, mesh,
								   cache::PropertyCache, material::Material{NC},
								   v_gas::Vector{Float64}, φ::Vector{Float64}) where NC

Compute gas mass fluxes at all faces using full Darcy-Fick formulation.

This combines advective transport (Darcy flow) with diffusive transport (volume fraction gradient):
	Jⱼ = (ξⱼ/φ_face)·v_g - (λ_face/T_face)·∂(ξⱼT)/∂x

ξ is per bulk volume and v_g is superficial, so the advected quantity is the
intrinsic gas density ξ/φ (see [`gas_flux_darcy_fick`](@ref)).

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `state`: Geometry/state cache (cell centers, face positions, T, ξ)
- `mesh`: Mesh with current state
- `cache`: Property cache (must have λ_face computed)
- `material`: Material definition
- `v_gas`: Darcy gas velocity at faces [m/s], **including** the boundary-face
  values written by `apply_pressure_bc_to_darcy_velocity!`
- `φ`: Porosity at cell centers [-] (e.g. `cache.darcy.φ`)

Positive flux = flow in +z direction (toward surface for our convention).

# Boundary faces
Boundary faces carry the **advective venting flux** from pressure BCs:
upwinding against a zero-concentration ambient ghost, only outflow advects,
J = v_b·ξ_int/φ_int (per gas species). The diffusive film-law part of the
boundary flux is the mass BC's domain and is **added** afterwards by
`apply_mass_bc_flux!`. With no pressure BC the boundary `v_gas` is 0 and the
boundary faces stay zero here.

# Notes
- When v_gas ≈ 0, this reduces to pure Fickian diffusion (same as compute_gas_fluxes!)
- Uses first-order upwinding for advective flux (stable, conservative)
"""
function compute_gas_fluxes_darcy_fick!(
    J::AbstractMatrix,
    state,
    mesh,
    cache::PropertyCache,
    material::Material{NC},
    v_gas::AbstractVector,
    φ::AbstractVector,
) where {NC}
    n_faces_val = mesh.n_active_cells + 1

    fill!(J, zero(eltype(J)))

    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 && right_id == 0
            continue
        end

        if left_id == 0 || right_id == 0
            # Boundary face: advective venting from the pressure BC's Darcy
            # velocity. The ambient ghost cell carries none of the tracked gas
            # species (ξ_ghost = 0), so upwinding admits only outflow:
            # interior-on-left vents at v_g ≥ 0 (top), interior-on-right at
            # v_g < 0 (bottom).
            v_g = v_gas[f]
            cell_id = left_id == 0 ? right_id : left_id
            outflow = left_id == 0 ? (v_g < 0) : (v_g >= 0)
            if outflow && v_g != 0
                v_over_φ = v_g / max(φ[cell_id], _PHI_FLOOR)
                ξ_int = ntuple(j -> state.ξ[j, cell_id], Val(NC))
                set_boundary_advective_flux_unrolled!(J, f, v_over_φ, ξ_int, material)
            end
            continue
        end

        d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])
        λ_face = cache.λ_face[f]
        T_L = state.T[left_id]
        T_R = state.T[right_id]
        v_g = v_gas[f]
        φ_L = φ[left_id]
        φ_R = φ[right_id]

        ξ_L = ntuple(j -> state.ξ[j, left_id], Val(NC))
        ξ_R = ntuple(j -> state.ξ[j, right_id], Val(NC))

        compute_gas_fluxes_darcy_fick_unrolled!(
            J,
            f,
            v_g,
            ξ_L,
            ξ_R,
            φ_L,
            φ_R,
            λ_face,
            T_L,
            T_R,
            d_LR,
            material,
        )
    end

    return nothing
end
