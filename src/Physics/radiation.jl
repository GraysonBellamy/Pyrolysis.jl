"""
Radiation transport models.

This module provides volumetric radiation source computation for different
radiation transport models: surface absorption and Beer-Lambert in-depth absorption.

Reference: Design Document Section 2.6, 4.3
"""

"""
	compute_surface_radiation_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
									   I_surface::Float64) where NC

Compute volumetric radiation source for surface-only absorption model.

All incident radiation is deposited as a volumetric source in the surface cell.
This is computationally efficient O(1) and appropriate for opaque materials.

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `mesh`: Mesh with geometry
- `I_surface`: Incident radiation intensity at surface [W/m²]

# Physics
For surface absorption, the volumetric source in the surface cell is:
	q'''_surface = I_surface / Δz_surface

All other cells have zero radiation source.

# Convention
- Surface is at z=L (cell N), substrate is at z=0 (cell 1)
- All radiation deposited in cell N

# Example
```julia
q_rad = zeros(n_cells)
compute_surface_radiation_source!(q_rad, mesh, 50e3)  # 50 kW/m² at surface
```
"""
function compute_surface_radiation_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
	I_surface::Float64) where NC
	n_cells = length(mesh.cells)

	# Zero all cells
	fill!(q_rad, 0.0)

	# Deposit all radiation in surface cell (cell n_cells is at z=L)
	if I_surface > 0.0 && n_cells > 0
		Δz_surface = cell_size(mesh, n_cells)
		@inbounds q_rad[n_cells] = I_surface / Δz_surface
	end

	return nothing
end

"""
	compute_beer_lambert_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
								  cache::PropertyCache, I_surface::Float64) where NC

Compute volumetric radiation absorption using Beer-Lambert law.

Radiation enters from the top surface (z=L) and propagates toward bottom (z=0):
	I(z) = I₀ × exp(-∫_L^z α_eff(z') dz')
	q'''(z) = α_eff(z) × I(z)

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `mesh`: Mesh with geometry
- `cache`: Property cache with absorption coefficients
- `I_surface`: Incident radiation intensity at surface [W/m²]

# Physics
The Beer-Lambert law describes exponential attenuation of radiation through
an absorbing medium. The optical thickness τ = α × Δz determines how much
radiation is absorbed in each cell.

For each cell, the exact cell-averaged volumetric source is:
	q̄''' = I_in × (1 - exp(-τ)) / Δz

# Convention
- Surface is at z=L (cell N), substrate is at z=0 (cell 1)
- Radiation propagates from cell N toward cell 1 (into depth)

# Example
```julia
# Semi-transparent PMMA with α = 1.47 m²/kg
cache.α .= α_pmma * ρ_pmma  # Convert to volumetric absorption coefficient
q_rad = zeros(n_cells)
compute_beer_lambert_source!(q_rad, mesh, cache, 50e3)
```
"""
function compute_beer_lambert_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
	cache::PropertyCache, I_surface::Float64) where NC
	n_cells_val = length(mesh.cells)

	if I_surface <= 0.0
		fill!(q_rad, 0.0)
		return nothing
	end

	# Integrate from surface (z=L, cell N) into depth (toward z=0, cell 1)
	I_current = I_surface

	@inbounds for i in n_cells_val:-1:1
		# Absorption coefficient for this cell
		α = cache.α[i]

		# Cell size
		Δz = cell_size(mesh, i)

		# Compute absorbed power density using exact cell-averaged Beer-Lambert
		if α > 1e-20
			# Optical thickness of this cell
			τ = α * Δz

			# Exact cell-averaged volumetric source:
			# q̄''' = (1/Δz) * ∫ α*I₀*exp(-α*z) dz = I_in * (1 - exp(-τ)) / Δz
			# This is the energy absorbed per unit volume
			absorption_factor = 1.0 - exp(-τ)
			q_rad[i] = I_current * absorption_factor / Δz

			# Update intensity for next cell (exact attenuation)
			I_current *= exp(-τ)
		else
			q_rad[i] = 0.0
		end
	end

	return nothing
end

"""
	compute_radiation_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
							  cache::PropertyCache, I_surface::Float64,
							  model::RadiationModel) where NC

Dispatch radiation source computation based on radiation model.

# Arguments
- `q_rad`: Output vector of volumetric heat source [W/m³]
- `mesh`: Mesh with geometry
- `cache`: Property cache with absorption coefficients (used by BEER_LAMBERT)
- `I_surface`: Incident radiation intensity at surface [W/m²]
- `model`: Radiation model enum value

# Behavior by Model
- `NO_RADIATION`: Sets q_rad to zero (no volumetric absorption)
- `SURFACE_ABSORPTION`: Deposits all flux in surface cell (O(1))
- `BEER_LAMBERT`: Computes in-depth absorption via Beer-Lambert law (O(n))
- `P1_QUASI_STEADY`: Handled separately via radiation cache (q_rad set to zero here)

# Example
```julia
# Surface absorption
compute_radiation_source!(q_rad, mesh, cache, 50e3, SURFACE_ABSORPTION)

# In-depth absorption
compute_radiation_source!(q_rad, mesh, cache, 50e3, BEER_LAMBERT)
```
"""
function compute_radiation_source!(q_rad::Vector{Float64}, mesh::Mesh1D{NC},
	cache::PropertyCache, I_surface::Float64,
	model::RadiationModel) where NC
	if model == NO_RADIATION || model == P1_QUASI_STEADY || model == SURFACE_ABSORPTION
		# NO_RADIATION: no volumetric source (radiation either ignored or handled at boundary)
		# SURFACE_ABSORPTION: radiation absorbed at boundary via thermal_bc_flux,
		#   NOT as volumetric source. This ensures boundary temperature solve gives
		#   physically correct T_surface (hotter than T_cell when net heating).
		# P1_QUASI_STEADY: handled via RadiationCache, not this function
		fill!(q_rad, 0.0)
	elseif model == BEER_LAMBERT
		compute_beer_lambert_source!(q_rad, mesh, cache, I_surface)
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
function total_absorbed_radiation(q_rad::Vector{Float64}, mesh::Mesh1D{NC}) where NC
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
function transmitted_radiation(I_surface::Float64, mesh::Mesh1D{NC}, cache::PropertyCache) where NC
	n_cells_val = length(mesh.cells)

	# Compute optical depth
	τ_total = 0.0
	@inbounds for i in 1:n_cells_val
		α = cache.α[i]
		Δz = cell_size(mesh, i)
		τ_total += α * Δz
	end

	return I_surface * exp(-τ_total)
end
