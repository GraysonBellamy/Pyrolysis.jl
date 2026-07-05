"""
Radiation transport models.

This module provides volumetric radiation source computation for different
radiation transport models: surface absorption and Beer-Lambert in-depth absorption.
"""

"""
	compute_beer_lambert_source!(q_rad, state, cache::PropertyCache, I_surface,
								  apply_kirchhoff_absorptivity::Bool)

Compute volumetric radiation absorption using the Beer-Lambert law.

Radiation enters from the top surface (z=L) and propagates toward bottom (z=0):
	I(z) = α_surf × I_surface × exp(-∫_z^L α_eff(z') dz')
	q'''(z) = α_eff(z) × I(z)

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `state`: Live state cache (cell volumes, cross-section area)
- `cache`: Property cache with absorption coefficients `α` and surface
  emissivities `ε`
- `I_surface`: Incident radiation intensity at the surface [W/m²]
- `apply_kirchhoff_absorptivity`: surface-transmissivity convention, see below

# Surface transmissivity α_surf
- `apply_kirchhoff_absorptivity = true` — the default Kirchhoff path (the
  `RadiativeFluxBC` carries no explicit `absorptivity`): `α_surf = cache.ε[end]`,
  the live surface cell's effective emissivity. Kirchhoff's law (α = ε) makes
  in-depth absorption consistent with surface *emission*, which uses the same
  ε_eff in the boundary Newton solve; the remaining (1 − ε_eff)·I_surface is
  reflected at the surface and never enters the medium.
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
- Surface is at z=L (cell N), substrate is at z=0 (cell 1)
- Radiation propagates from cell N toward cell 1 (into depth)
"""
function compute_beer_lambert_source!(
    q_rad::AbstractVector,
    state,
    cache::PropertyCache,
    I_surface::Real,
    apply_kirchhoff_absorptivity::Bool,
)
    n_cells_val = length(state.T)

    if I_surface <= 0 || n_cells_val == 0
        fill!(q_rad, zero(eltype(q_rad)))
        return nothing
    end

    # Integrate from surface (z=L, cell N) into depth (toward z=0, cell 1).
    # The flux entering the medium is α_surf·I_surface (Kirchhoff α_surf = ε_eff
    # of the live surface cell, or 1 when the BC's explicit absorptivity is
    # already folded into I_surface upstream).
    I_current = oneunit(eltype(q_rad)) * I_surface
    if apply_kirchhoff_absorptivity
        I_current *= cache.ε[n_cells_val]
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
            cache,
            I_surface,
            apply_kirchhoff_absorptivity,
        )
    else
        error("Unknown radiation model: $model")
    end
    return nothing
end


# ============================================================================
# Utility Functions
# ============================================================================

"""
	total_absorbed_radiation(q_rad::Vector{Float64}, mesh::Mesh1D{NC}) where NC -> Float64

Compute total absorbed radiation [W/m²].

Integrates volumetric radiation source over the domain to verify energy conservation.
"""
function total_absorbed_radiation(q_rad::Vector{Float64}, mesh::Mesh1D{NC}) where {NC}
    Q_total = 0.0
    @inbounds for i in eachindex(q_rad)
        Q_total += q_rad[i] * cell_size(mesh, i)
    end
    return Q_total
end

"""
	transmitted_radiation(I_surface::Float64, mesh::Mesh1D{NC}, cache::PropertyCache) where NC -> Float64

Compute radiation transmitted through the material [W/m²].

Only applicable for Beer-Lambert model. Returns I_surface × exp(-τ_total) where
τ_total is the total optical depth through the domain.
"""
function transmitted_radiation(
    I_surface::Float64,
    mesh::Mesh1D{NC},
    cache::PropertyCache,
) where {NC}
    n_cells_val = length(mesh.cells)

    # Compute optical depth
    τ_total = 0.0
    @inbounds for i = 1:n_cells_val
        α = cache.α[i]
        Δz = cell_size(mesh, i)
        τ_total += α * Δz
    end

    return I_surface * exp(-τ_total)
end
