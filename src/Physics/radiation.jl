"""
Radiation transport models.

This module provides volumetric radiation source computation for different
radiation transport models: surface absorption and Beer-Lambert in-depth absorption.
"""

"""
	compute_beer_lambert_source!(q_rad, state, mesh, cache::PropertyCache, I_surface,
								  apply_kirchhoff_absorptivity::Bool)

Compute volumetric radiation absorption using the Beer-Lambert law.

Radiation enters at the live surface (face `mesh.n_active_cells + 1`) and
propagates toward the bottom (z=0):
	I(z) = α_surf × I_surface × exp(-∫_z^L α_eff(z') dz')
	q'''(z) = α_eff(z) × I(z)

Parked ghost cells above the live surface receive `q_rad = 0`.

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `state`: Live state cache (cell volumes, cross-section area)
- `mesh`: Mesh (live-cell count)
- `cache`: Property cache with absorption coefficients `α` and surface
  absorptivities `a`
- `I_surface`: Incident radiation intensity at the surface [W/m²]
- `apply_kirchhoff_absorptivity`: surface-transmissivity convention, see below

# Surface transmissivity α_surf
- `apply_kirchhoff_absorptivity = true` — the default path (the
  `RadiativeFluxBC` carries no explicit `absorptivity`): `α_surf = cache.a[end]`,
  the live surface cell's effective surface absorptivity
  `a_eff = Σⱼ vⱼ aⱼ`. With no component `surface_absorptivity` set this is
  Kirchhoff's law, `a_eff = ε_eff(T)`, so in-depth absorption stays
  consistent with the surface *emission* of the boundary Newton solve; the
  remaining (1 − a_eff)·I_surface is reflected at the surface and never
  enters the medium.
- `apply_kirchhoff_absorptivity = false` — the BC supplied an explicit
  `absorptivity`, which `extract_radiative_flux_from_bcs` has already folded
  into `I_surface`; `α_surf = 1` here so it is not applied twice.

# Physics
The Beer-Lambert law describes exponential attenuation of radiation through
an absorbing medium. The optical thickness τ = α × Δz determines how much
radiation is absorbed in each cell.

For each cell, the exact cell-averaged volumetric source is:
	q̄''' = I_in × (1 - exp(-τ)) / Δz

# Convention
- Surface is the live top cell `n_live = mesh.n_active_cells`, substrate is at z=0 (cell 1)
- Radiation propagates from cell `n_live` toward cell 1 (into depth)
"""
function compute_beer_lambert_source!(
    q_rad::AbstractVector,
    state,
    mesh,
    cache::PropertyCache,
    I_surface::Real,
    apply_kirchhoff_absorptivity::Bool,
)
    n_cells_val = mesh.n_active_cells
    fill!(q_rad, zero(eltype(q_rad)))

    if I_surface <= 0 || n_cells_val == 0
        return nothing
    end

    # Integrate from surface (z=L, cell N) into depth (toward z=0, cell 1).
    # The flux entering the medium is α_surf·I_surface (α_surf = a_eff of the
    # live surface cell — ε_eff under Kirchhoff — or 1 when the BC's explicit
    # absorptivity is already folded into I_surface upstream).
    I_current = oneunit(eltype(q_rad)) * I_surface
    if apply_kirchhoff_absorptivity
        I_current *= cache.a[n_cells_val]
    end

    A = state.A_section[]
    @inbounds for i = n_cells_val:-1:1
        α = cache.α[i]
        Δz = state.volume[i] / A

        if α > 1e-20
            τ = α * Δz
            absorption_factor = 1 - exp(-τ)
            q_rad[i] = I_current * absorption_factor / Δz
            I_current *= exp(-τ)
        else
            q_rad[i] = zero(eltype(q_rad))
        end
    end

    return nothing
end

"""
	compute_radiation_source!(q_rad, state, mesh, cache::PropertyCache,
							  I_surface, model::RadiationModel,
							  apply_kirchhoff_absorptivity::Bool)

Dispatch volumetric radiation source computation based on radiation model.

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `state`: Live state cache (cell volumes, cross-section area)
- `mesh`: Mesh (connectivity only; geometry comes from `state`)
- `cache`: Property cache with absorption coefficients (used by BEER_LAMBERT)
- `I_surface`: Incident radiation intensity at surface [W/m²]
- `model`: Radiation model enum value
- `apply_kirchhoff_absorptivity`: scale the entering flux by the surface
  cell's ε_eff (Kirchhoff default path); `false` when the BC's explicit
  absorptivity is already folded into `I_surface`. Only BEER_LAMBERT
  consumes it. See [`compute_beer_lambert_source!`](@ref).

# Behavior by Model
- `NO_RADIATION`: Sets q_rad to zero (no volumetric absorption)
- `SURFACE_ABSORPTION`: Sets q_rad to zero — the absorbed flux α·q_ext enters
  through the surface energy balance (boundary Newton solve), not the volume
- `BEER_LAMBERT`: Computes in-depth absorption via Beer-Lambert law (O(n))
- `P1_QUASI_STEADY`: Handled separately via radiation cache (q_rad set to zero here)
"""
function compute_radiation_source!(
    q_rad::AbstractVector,
    state,
    mesh,
    cache::PropertyCache,
    I_surface::Real,
    model::RadiationModel,
    apply_kirchhoff_absorptivity::Bool,
)
    if model == NO_RADIATION || model == P1_QUASI_STEADY || model == SURFACE_ABSORPTION
        fill!(q_rad, zero(eltype(q_rad)))
    elseif model == BEER_LAMBERT
        compute_beer_lambert_source!(
            q_rad,
            state,
            mesh,
            cache,
            I_surface,
            apply_kirchhoff_absorptivity,
        )
    else
        error("Unknown radiation model: $model")
    end
    return nothing
end
