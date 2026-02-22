"""
P1 approximation for volumetric radiation transport.

Implements the first-order spherical harmonics (P1) method for solving
radiative transfer in optically thick to intermediate media. This is the
industry-standard approach for pyrolysis modeling where optical thickness
varies during decomposition.

The P1 method solves for incident radiation intensity G [W/m²]:
	∇·(D_rad ∇G) = α G - 4σεαT⁴

where:
- D_rad = 1/(3α) is the radiation diffusion coefficient
- α is absorption coefficient [m⁻¹]
- ε is emissivity [-]
- σ is Stefan-Boltzmann constant

Reference: 
- P1_Volumetric_Reradiation_Implementation.md
- Modest, M.F., Radiative Heat Transfer, 3rd Ed., 2013, Chapter 15
"""

# ============================================================================
# Core P1 Functions
# ============================================================================

"""
	compute_radiation_diffusion_coefficient(α::Float64, σ_s::Float64=0.0) -> Float64

Compute the P1 radiation diffusion coefficient.

	D_rad = 1 / (3(α + σ_s))

For pyrolysis materials, scattering is typically negligible (σ_s ≈ 0):
	D_rad = 1 / (3α)

# Arguments
- `α`: Absorption coefficient [m⁻¹]. Must be positive for finite D_rad.
- `σ_s`: Scattering coefficient [m⁻¹] (default: 0). Typically zero for pyrolysis.

# Returns
- `D_rad`: Radiation diffusion coefficient [m]

# Physical Interpretation
D_rad represents the "mean free path" for radiation diffusion. Large D_rad
(small α) means radiation travels far before being absorbed. Small D_rad
(large α) means radiation is absorbed locally.

# Examples
```julia
# Typical virgin polymer (α ≈ 1000-2000 m⁻¹)
D_rad = compute_radiation_diffusion_coefficient(1500.0)  # ≈ 2.2×10⁻⁴ m

# Char layer (α ≈ 5000-10000 m⁻¹)
D_rad = compute_radiation_diffusion_coefficient(8000.0)  # ≈ 4.2×10⁻⁵ m

# With scattering (rare for pyrolysis)
D_rad = compute_radiation_diffusion_coefficient(5000.0, 1000.0)
```
"""
@inline function compute_radiation_diffusion_coefficient(α::Float64, σ_s::Float64 = 0.0)
	# Guard against zero/negative absorption (transparent material)
	total_extinction = α + σ_s
	if total_extinction < 1e-10
		return 0.0  # Transparent material: no diffusion
	end
	return 1.0 / (3.0 * total_extinction)
end

"""
	compute_radiation_source_terms!(Q_emission::Vector{Float64},
									Q_absorption::Vector{Float64},
									mesh::Mesh1D{NC},
									rad_cache::RadiationCache,
									prop_cache::PropertyCache) where NC

Compute volumetric emission and absorption source terms for P1 equation.

# Emission Source (blackbody emission into radiation field):
	Q_emission[i] = 4σ × ε[i] × α[i] × T[i]⁴    [W/m³]

# Absorption (removal from radiation field):
	Q_absorption[i] = α[i] × G[i]                [W/m³]

The net radiation source for the energy equation is:
	Q_net = Q_absorption - Q_emission

Note: These are the sources for the G equation. For the T equation,
the radiation contribution comes through ∇·(D_rad ∇G).

# Arguments
- `Q_emission`: Output vector for emission sources [W/m³]
- `Q_absorption`: Output vector for absorption sources [W/m³]
- `mesh`: Mesh with current temperature state
- `rad_cache`: Radiation cache with G values
- `prop_cache`: Property cache with α, ε values
"""
function compute_radiation_source_terms!(Q_emission::Vector{Float64},
	Q_absorption::Vector{Float64},
	mesh,
	rad_cache::RadiationCache,
	prop_cache::PropertyCache)
	n_cells = length(mesh.cells)

	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T
		G = rad_cache.G[i]
		# Use volumetric α from rad_cache if populated, otherwise from prop_cache
		# Note: prop_cache.α is already volumetric (effective_absorption = Σⱼ ξⱼ αⱼ [m⁻¹])
		α = length(rad_cache.α) >= i && rad_cache.α[i] > 0 ?
			rad_cache.α[i] : prop_cache.α[i]

		# Emission: 4σαT⁴ (blackbody emission, NO emissivity factor)
		# Physical: Energy emitted per unit volume into radiation field
		# Note: In standard P1 approximation, volumetric emission is blackbody.
		# Emissivity only appears at boundaries (Marshak BC) where it represents
		# surface reflectivity/absorptivity. The interior is assumed to be a
		# perfect absorber/emitter (gray body with ε=1 locally).
		Q_emission[i] = 4.0 * STEFAN_BOLTZMANN * α * T^4

		# Absorption: αG
		# Physical: Energy absorbed from radiation field per unit volume
		Q_absorption[i] = α * G
	end

	return nothing
end

"""
	update_radiation_properties!(rad_cache::RadiationCache,
								 mesh::Mesh1D{NC},
								 prop_cache::PropertyCache) where NC

Update radiation-specific properties in the cache.

Computes:
1. Cell-centered D_rad[i] = 1/(3α[i])
2. Face-centered D_rad_face[f] = harmonic mean weighted by distance

The harmonic mean is used for D_rad at faces because radiation transport
follows resistance-in-series (like conduction). This ensures correct
flux at interfaces between materials with different optical properties.

# Arguments
- `rad_cache`: Radiation cache to update
- `mesh`: Mesh with geometry information
- `prop_cache`: Property cache with α values (must be updated first)
"""
function update_radiation_properties!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache)
	n_cells = length(mesh.cells)
	n_faces = length(mesh.faces)

	# Ensure cache is properly sized
	if length(rad_cache.D_rad) != n_cells
		resize_radiation_cache!(rad_cache, n_cells)
	end

	# Cell-centered D_rad and volumetric α
	@inbounds for i in 1:n_cells
		# prop_cache.α is already volumetric [m⁻¹] from effective_absorption = Σⱼ ξⱼ αⱼ
		# where ξⱼ [kg/m³] × αⱼ [m²/kg] = [m⁻¹]
		α_vol = prop_cache.α[i]

		# Store volumetric α in rad_cache
		rad_cache.α[i] = α_vol

		# D_rad = 1/(3α) for volumetric absorption coefficient
		rad_cache.D_rad[i] = compute_radiation_diffusion_coefficient(α_vol)
	end

	# Face-centered D_rad using distance-weighted harmonic mean
	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids
		face_pos = mesh.face_positions[f]

		if left_id == 0
			# Left boundary: use right cell value
			rad_cache.D_rad_face[f] = rad_cache.D_rad[right_id]
		elseif right_id == 0
			# Right boundary: use left cell value
			rad_cache.D_rad_face[f] = rad_cache.D_rad[left_id]
		else
			# Interior face: distance-weighted harmonic mean
			d_L = abs(face_pos - mesh.cell_centers[left_id])
			d_R = abs(mesh.cell_centers[right_id] - face_pos)

			rad_cache.D_rad_face[f] = harmonic_mean_weighted(
				rad_cache.D_rad[left_id],
				rad_cache.D_rad[right_id],
				d_L,
				d_R,
			)
		end
	end

	return nothing
end

# ============================================================================
# Radiation Flux Computation
# ============================================================================

"""
	compute_radiation_fluxes!(rad_cache::RadiationCache,
							  mesh::Mesh1D{NC},
							  prop_cache::PropertyCache,
							  radiation_bcs::RadiationBCSet,
							  t) where NC

Compute radiation diffusion fluxes at all faces using P1 approximation.

For interior faces:
	q_rad[f] = -D_rad_face × (G[R] - G[L]) / Δx

For boundary faces, the Marshak or other BC is applied.

# Arguments
- `rad_cache`: Radiation cache (updated in-place: q_rad)
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α, ε values
- `radiation_bcs`: RadiationBCSet with radiation boundary conditions
- `t`: Current simulation time [s] (supports dual numbers for AD)

# Notes
- Type-generic to support automatic differentiation
"""
function compute_radiation_fluxes!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	radiation_bcs::RadiationBCSet,
	t)
	n_faces = length(mesh.faces)

	# Interior faces (radiation diffusion)
	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			continue
		end

		D_rad_face = rad_cache.D_rad_face[f]
		Δx = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		rad_cache.q_rad[f] = -D_rad_face * (rad_cache.G[right_id] - rad_cache.G[left_id]) / Δx
	end

	_apply_all_radiation_bcs!(rad_cache, mesh, prop_cache, radiation_bcs, t)

	return nothing
end

# RadiationBCSet (type-stable) boundary application
@inline function _apply_all_radiation_bcs!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	radiation_bcs::RadiationBCSet,
	t)
	faces_top = get(mesh.boundary_faces, 1, nothing)
	if faces_top !== nothing && !isempty(faces_top)
		_apply_bc!(rad_cache, mesh, prop_cache, radiation_bcs.top, faces_top::Vector{Int}, t)
	end

	faces_bottom = get(mesh.boundary_faces, 2, nothing)
	if faces_bottom !== nothing && !isempty(faces_bottom)
		_apply_bc!(rad_cache, mesh, prop_cache, radiation_bcs.bottom, faces_bottom::Vector{Int}, t)
	end
	return nothing
end

@inline function _apply_bc_for_boundary!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	radiation_bcs::RadiationBCSet,
	::Val{B},
	t) where B
	face_ids_or_nothing = get(mesh.boundary_faces, B, nothing)
	face_ids_or_nothing === nothing && return nothing
	faces = face_ids_or_nothing::Vector{Int}
	isempty(faces) && return nothing

	bc = get_radiation_bc(radiation_bcs, Val(B))
	_apply_bc!(rad_cache, mesh, prop_cache, bc, faces, t)
	return nothing
end

# Dispatch helpers for RadiationBCSet (method-based, no isa chain)
@inline function _apply_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::MarshakRadiationBC,
	faces::Vector{Int},
	t)
	_apply_marshak_bc!(rad_cache, mesh, prop_cache, bc, faces, t)
	return nothing
end

@inline function _apply_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::ReflectiveRadiationBC,
	faces::Vector{Int},
	t)
	_apply_reflective_bc!(rad_cache, faces)
	return nothing
end

@inline function _apply_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::DirichletRadiationBC,
	faces::Vector{Int},
	t)
	_apply_dirichlet_bc!(rad_cache, faces)
	return nothing
end

function _apply_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::AbstractRadiationBC,
	faces::Vector{Int},
	t)
	_apply_generic_bc!(rad_cache, mesh, prop_cache, bc, faces, t)
	return nothing
end

# Helper to extract real value from Dual or Float64 for cache storage
# This is needed because RadiationCache uses Float64 arrays but may be called
# with Dual numbers during autodiff for time derivatives
@inline _value(x::Real) = ForwardDiff.value(x)
@inline _value(x::Float64) = x

# Type-specific implementations - these are called with concrete types
@inline function _apply_marshak_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::MarshakRadiationBC,
	faces::Vector{Int},
	t)
	@inbounds for f in faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids
		cell_id = left_id == 0 ? right_id : left_id

		T_surf = mesh.cell_states[cell_id].T
		G_surf = rad_cache.G[cell_id]
		α_surf = rad_cache.α[cell_id] > 0 ? rad_cache.α[cell_id] : prop_cache.α[cell_id]
		ε_surf = prop_cache.ε[cell_id]

		flux_into = radiation_bc_flux(bc, G_surf, T_surf, α_surf, ε_surf, t)
		# Extract value for storage in Float64 cache (needed for autodiff with Dual numbers)
		flux_value = _value(flux_into)
		if right_id == 0
			rad_cache.q_rad[f] = -flux_value
		else
			rad_cache.q_rad[f] = flux_value
		end
	end
	return nothing
end

@inline function _apply_reflective_bc!(rad_cache::RadiationCache, faces::Vector{Int})
	# Reflective BC: zero flux
	@inbounds for f in faces
		rad_cache.q_rad[f] = 0.0
	end
	return nothing
end

@inline function _apply_dirichlet_bc!(rad_cache::RadiationCache, faces::Vector{Int})
	# Dirichlet BC: flux determined by interior gradient (handled in discretization)
	# For now, set to zero as placeholder
	@inbounds for f in faces
		rad_cache.q_rad[f] = 0.0
	end
	return nothing
end

# Fallback for any other AbstractRadiationBC (may allocate but handles extensibility)
function _apply_generic_bc!(rad_cache::RadiationCache,
	mesh,
	prop_cache::PropertyCache,
	bc::AbstractRadiationBC,
	faces::Vector{Int},
	t)
	@inbounds for f in faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids
		cell_id = left_id == 0 ? right_id : left_id

		T_surf = mesh.cell_states[cell_id].T
		G_surf = rad_cache.G[cell_id]
		α_surf = rad_cache.α[cell_id] > 0 ? rad_cache.α[cell_id] : prop_cache.α[cell_id]
		ε_surf = prop_cache.ε[cell_id]

		flux_into = radiation_bc_flux(bc, G_surf, T_surf, α_surf, ε_surf, t)
		if right_id == 0
			rad_cache.q_rad[f] = -flux_into
		else
			rad_cache.q_rad[f] = flux_into
		end
	end
	return nothing
end

# ============================================================================
# Boundary Condition Flux Functions
# ============================================================================

"""
	radiation_bc_flux(bc::MarshakRadiationBC, G_surf, T_surf, α_surf, ε_surf, t)

Compute radiation flux at boundary with Marshak condition.

The standard P1 Marshak condition for a gray diffuse surface with emissivity ε:

	q = I_ext + (ε/2)(G_ext - G_surf)

where:
- I_ext is external incident radiation [W/m²]
- G_ext = 4σT_amb⁴ is the equilibrium incident radiation from ambient
- G_surf is the P1 incident radiation at the surface
- ε is surface emissivity (1-ε is reflected back into domain)

This formulation ensures:
1. At equilibrium (G = 4σT⁴, T_amb = T_surf, I_ext = 0): q = 0 ✓
2. External flux I_ext adds directly to the incoming radiation
3. Surface partially reflects (1-ε) of interior radiation back

# Arguments
- `bc`: Marshak boundary condition
- `G_surf`: Incident radiation at surface [W/m²]
- `T_surf`: Surface temperature [K] (not used in standard P1 Marshak)
- `α_surf`: Surface absorption coefficient [m⁻¹] (not used)
- `ε_surf`: Surface emissivity [-] (from prop_cache, may differ from BC)
- `t`: Current time [s] (supports dual numbers for AD)

# Returns
- Radiation flux [W/m²], positive = into domain

# Notes
- Type-generic to support automatic differentiation with Dual numbers
- Uses BC emissivity, not ε_surf from prop_cache, for consistency
"""
@inline function radiation_bc_flux(bc::MarshakRadiationBC,
	G_surf,
	T_surf,
	α_surf,
	ε_surf,
	t)
	# Evaluate BC parameters (preserves type for AD)
	I_ext = evaluate_param(bc.I_external, t)
	T_amb = evaluate_param(bc.T_ambient, t)
	ε = evaluate_param_T(bc.emissivity, T_surf)

	# External equilibrium radiation field from ambient temperature
	G_ext = 4.0 * STEFAN_BOLTZMANN * T_amb^4

	# Standard P1 Marshak BC flux (positive = into domain)
	# q = I_ext + (ε/2)(G_ext - G_surf)
	# Physical interpretation:
	# - I_ext: external source directly enters
	# - (ε/2)(G_ext - G_surf): net exchange with ambient, scaled by surface absorptivity
	flux = I_ext + 0.5 * ε * (G_ext - G_surf)

	return flux
end

"""
	radiation_bc_flux(::ReflectiveRadiationBC, args...)

Zero flux for reflective boundary (∂G/∂n = 0).
"""
@inline function radiation_bc_flux(::ReflectiveRadiationBC,
	G_surf,
	T_surf,
	α_surf,
	ε_surf,
	t)
	return zero(t)  # Preserve type for AD
end

"""
	radiation_bc_flux(bc::DirichletRadiationBC, args...)

For Dirichlet BCs, the flux is determined by the interior gradient.
This returns zero as a placeholder - actual implementation needs the interior G value.
"""
@inline function radiation_bc_flux(bc::DirichletRadiationBC,
	G_surf,
	T_surf,
	α_surf,
	ε_surf,
	t)
	# Dirichlet BC sets G directly; flux determined by interior gradient
	# This is handled specially in the discretization
	return zero(t)  # Preserve type for AD
end

# ============================================================================
# Optical Thickness Diagnostics
# ============================================================================

"""
	optical_thickness(mesh::Mesh1D{NC}, prop_cache::PropertyCache) where NC -> Float64

Compute the optical thickness τ = α × L for the domain.

This is a key parameter for assessing P1 validity:
- τ > 10: Optically thick (P1 excellent)
- 1 < τ < 10: Intermediate (P1 good)  
- 0.1 < τ < 1: Optically thin (P1 fair)
- τ < 0.1: Transparent (P1 poor)

# Arguments
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α values

# Returns
- τ: Optical thickness (dimensionless)

# Example
```julia
τ = optical_thickness(mesh, prop_cache)
if τ < 0.1
	@warn "Material is nearly transparent (τ = \$τ). P1 may be inaccurate."
end
```
"""
function optical_thickness(mesh, prop_cache::PropertyCache)
	n_cells = length(mesh.cells)

	# Compute total optical depth: τ = ∫ α dx
	τ_total = 0.0
	@inbounds for i in 1:n_cells
		α = prop_cache.α[i]
		Δx = cell_size(mesh, i)
		τ_total += α * Δx
	end

	return τ_total
end

"""
	optical_thickness_cell(mesh::Mesh1D{NC}, prop_cache::PropertyCache, i::Int) where NC -> Float64

Compute the optical thickness of a single cell: τ_i = α_i × Δx_i

This is useful for identifying which cells are optically thick vs thin.
"""
@inline function optical_thickness_cell(mesh, prop_cache::PropertyCache, i::Int)
	α = prop_cache.α[i]
	Δx = cell_size(mesh, i)
	return α * Δx
end

"""
	check_p1_validity(mesh, prop_cache::PropertyCache; 
					  warn::Bool=true) -> NamedTuple

Check if P1 approximation is valid for the current material state.

Returns a NamedTuple with:
- `τ_total`: Total optical thickness
- `τ_min`: Minimum cell optical thickness
- `τ_max`: Maximum cell optical thickness
- `regime`: String describing optical regime
- `valid`: Bool indicating if P1 is appropriate

# Arguments
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α values
- `warn`: If true, print warning for marginal cases

# Example
```julia
validity = check_p1_validity(mesh, prop_cache)
println("Optical regime: \$(validity.regime)")
println("P1 valid: \$(validity.valid)")
```
"""
function check_p1_validity(mesh, prop_cache::PropertyCache; warn::Bool = true)
	n_cells = length(mesh.cells)

	τ_total = optical_thickness(mesh, prop_cache)

	τ_min = Inf
	τ_max = 0.0
	@inbounds for i in 1:n_cells
		τ_i = optical_thickness_cell(mesh, prop_cache, i)
		τ_min = min(τ_min, τ_i)
		τ_max = max(τ_max, τ_i)
	end

	# Classify regime
	if τ_total >= 10
		regime = "optically thick (τ = $(round(τ_total, sigdigits=3)))"
		valid = true
	elseif τ_total > 1
		regime = "intermediate (τ = $(round(τ_total, sigdigits=3)))"
		valid = true
	elseif τ_total > 0.1
		regime = "optically thin (τ = $(round(τ_total, sigdigits=3)))"
		valid = true  # P1 still usable but less accurate
		if warn
			@warn "Material is optically thin (τ = $τ_total). P1 approximation may have reduced accuracy (15-30% error)."
		end
	else
		regime = "transparent (τ = $(round(τ_total, sigdigits=3)))"
		valid = false
		if warn
			@warn "Material is nearly transparent (τ = $τ_total). P1 approximation may be inaccurate (>30% error). Consider discrete ordinates method."
		end
	end

	return (τ_total = τ_total, τ_min = τ_min, τ_max = τ_max, regime = regime, valid = valid)
end

# ============================================================================
# Energy Conservation Diagnostics
# ============================================================================

"""
	radiation_energy_balance(rad_cache::RadiationCache, mesh) -> NamedTuple

Compute radiation energy balance for the domain.

Returns:
- `total_emission`: ∫ Q_emission dV [W/m²] (integrated over depth)
- `total_absorption`: ∫ Q_absorption dV [W/m²]
- `net_source`: total_absorption - total_emission [W/m²]
- `flux_top`: Radiation flux at top boundary [W/m²]
- `flux_bottom`: Radiation flux at bottom boundary [W/m²]
- `balance_error`: Should be close to zero for steady state
"""
function radiation_energy_balance(rad_cache::RadiationCache, mesh)
	n_cells = length(mesh.cells)

	total_emission = 0.0
	total_absorption = 0.0

	@inbounds for i in 1:n_cells
		Δx = cell_size(mesh, i)
		total_emission += rad_cache.Q_emission[i] * Δx
		total_absorption += rad_cache.Q_absorption[i] * Δx
	end

	net_source = total_absorption - total_emission

	# Boundary fluxes (top = last face, bottom = first face)
	flux_top = rad_cache.q_rad[end]
	flux_bottom = rad_cache.q_rad[1]

	# For steady state: flux_in - flux_out = net_source
	balance_error = (flux_top - flux_bottom) - net_source

	return (
		total_emission = total_emission,
		total_absorption = total_absorption,
		net_source = net_source,
		flux_top = flux_top,
		flux_bottom = flux_bottom,
		balance_error = balance_error,
	)
end

# ============================================================================
# P1 Approximation Validity Checking
# ============================================================================

# P1 validity thresholds based on optical thickness τ = α × L
const P1_TAU_MIN = 1.0      # Minimum for any validity
const P1_TAU_MARGINAL = 3.0 # Below this, P1 is approximate
const P1_TAU_GOOD = 5.0     # Above this, P1 is accurate

"""
	check_p1_validity(mesh::Mesh1D{NC}, material::Material{NC}, prop_cache::PropertyCache;
	                  warn::Bool=true, error_on_invalid::Bool=false) where NC -> Bool

Check if P1 approximation is valid for current mesh and material.

P1 requires optical thickness τ >> 1 (typically τ > 3-5) for accuracy.

The optical thickness is defined as:
	τ = α × L

where α is the absorption coefficient [m⁻¹] and L is the characteristic length [m].

# Validity Criteria
- τ < 1.0: Invalid - P1 not applicable
- 1.0 ≤ τ < 3.0: Marginal - P1 gives approximate results
- 3.0 ≤ τ < 5.0: Acceptable - P1 reasonably accurate
- τ ≥ 5.0: Good - P1 highly accurate

# Arguments
- `mesh`: Current mesh state
- `material`: Material with absorption coefficients
- `prop_cache`: Property cache with α values
- `warn`: If true, print warning for marginally valid cases
- `error_on_invalid`: If true, throw error when τ < 1

# Returns
`true` if P1 is valid (τ > τ_min for all cells), `false` otherwise

# Example
```julia
# Check validity before using P1
if !check_p1_validity(mesh, material, prop_cache)
    @warn "P1 radiation may be inaccurate - consider BEER_LAMBERT or SURFACE_ABSORPTION"
end
```
"""
function check_p1_validity(mesh::Mesh1D{NC}, material::Material{NC}, prop_cache::PropertyCache;
	warn::Bool = true,
	error_on_invalid::Bool = false) where NC

	τ_min = Inf
	invalid_cells = Int[]
	marginal_cells = Int[]
	acceptable_cells = Int[]

	@inbounds for i in eachindex(mesh.cells)
		if !mesh.active_cell_mask[i]
			continue
		end

		# Compute local optical thickness
		α = prop_cache.α[i]
		Δz = mesh.cell_volumes[i] / mesh.cross_section_area
		τ_local = α * Δz

		τ_min = min(τ_min, τ_local)

		if τ_local < P1_TAU_MIN
			push!(invalid_cells, i)
		elseif τ_local < P1_TAU_MARGINAL
			push!(marginal_cells, i)
		elseif τ_local < P1_TAU_GOOD
			push!(acceptable_cells, i)
		end
	end

	# Report findings
	if !isempty(invalid_cells) && error_on_invalid
		error("""
			P1 radiation invalid: $(length(invalid_cells)) cells have τ < $(P1_TAU_MIN)
			Minimum optical thickness: $τ_min
			Consider using BEER_LAMBERT or SURFACE_ABSORPTION model instead.
			""")
	end

	if warn
		if !isempty(invalid_cells)
			@warn "P1 radiation may be inaccurate" n_invalid = length(invalid_cells) τ_min
		elseif !isempty(marginal_cells)
			@warn "P1 radiation marginally valid" n_marginal = length(marginal_cells) τ_min
		elseif !isempty(acceptable_cells)
			# Only inform, not warn
			# @info "P1 radiation acceptable but not optimal" n_acceptable=length(acceptable_cells) τ_min
		end
	end

	return isempty(invalid_cells)
end

"""
	compute_optical_thickness(mesh::Mesh1D{NC}, prop_cache::PropertyCache) where NC -> Vector{Float64}

Compute optical thickness τ = α × Δz for each cell.

# Returns
Vector of optical thickness values [dimensionless] for each cell.

# Example
```julia
τ = compute_optical_thickness(mesh, prop_cache)
println("Minimum optical thickness: ", minimum(τ))
println("Maximum optical thickness: ", maximum(τ))
```
"""
function compute_optical_thickness(mesh::Mesh1D{NC}, prop_cache::PropertyCache) where NC
	n_cells = length(mesh.cells)
	τ = zeros(n_cells)

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			α = prop_cache.α[i]
			Δz = mesh.cell_volumes[i] / mesh.cross_section_area
			τ[i] = α * Δz
		end
	end

	return τ
end
