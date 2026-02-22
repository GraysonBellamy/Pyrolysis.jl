"""
Gas transport (Darcy-Fick).

Reference: Design Document Section 2.4, 4.4
Reference: ThermaKin2Ds Reference Manual, Section 1.2

ThermaKin formulation:
	J_g = -λ ρ ∂(ξ_g/ρ)/∂x = -λ ρ ∂Y_g/∂x

This is flux driven by mass fraction gradient (thermodynamically consistent
for porous media transport), not concentration gradient.

At material interfaces, flux continuity requires the harmonic mean of the (λ·ρ)
product, not separate interpolation of λ and ρ. This is computed as λρ_face
in the property cache.
"""

"""
	gas_flux(λρ_face::Float64, ξ_L::Float64, ξ_R::Float64, 
			 ρ_L::Float64, ρ_R::Float64, d_LR::Float64) -> Float64

Compute gas mass flux at a face [kg/(m²·s)].

Uses rigorous formulation with harmonic mean of λ·ρ product:
	J = -λρ_face (Y_R - Y_L) / d_LR

where Y = ξ/ρ is mass fraction and λρ_face is the harmonic mean of (λ·ρ).
This ensures physically correct flux at interfaces with sharp property jumps.

# Arguments
- `λρ_face`: Face transport coefficient (harmonic mean of λ·ρ) [kg/(m·s)]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `ρ_L, ρ_R`: Total mixture densities on left and right [kg/m³]
- `d_LR`: Distance between cell centers [m]
"""
@inline function gas_flux(λρ_face::Float64, ξ_L::Float64, ξ_R::Float64,
	ρ_L::Float64, ρ_R::Float64, d_LR::Float64)
	# Rigorous formulation: J = -λρ_face ∂Y/∂x
	# where λρ_face = harmonic_mean(λ_L·ρ_L, λ_R·ρ_R)
	Y_L = ξ_L / ρ_L
	Y_R = ξ_R / ρ_R
	return -λρ_face * (Y_R - Y_L) / d_LR
end

# ============================================================================
# Darcy-Fick Combined Advection-Diffusion Flux
# ============================================================================

"""
	gas_flux_darcy_fick(v_g::Float64, ξ_L::Float64, ξ_R::Float64,
						λρ_face::Float64, ρ_L::Float64, ρ_R::Float64, 
						d_LR::Float64) -> Float64

Compute gas mass flux at a face using full Darcy-Fick formulation [kg/(m²·s)].

Combines advective and diffusive transport:
	J = ξ_upwind × v_g - λρ_face × (Y_R - Y_L) / d_LR
		↑ advection      ↑ diffusion

# Arguments
- `v_g`: Darcy gas velocity at face [m/s] (positive = flow in +z direction)
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `λρ_face`: Face transport coefficient (harmonic mean of λ·ρ) [kg/(m·s)]
- `ρ_L, ρ_R`: Total mixture densities on left and right [kg/m³]
- `d_LR`: Distance between cell centers [m]

# Returns
Total flux [kg/(m²·s)], positive = flow in +z direction
"""
@inline function gas_flux_darcy_fick(v_g::Float64, ξ_L::Float64, ξ_R::Float64,
	λρ_face::Float64, ρ_L::Float64, ρ_R::Float64, d_LR::Float64)
	# Diffusive flux (same as gas_flux)
	Y_L = ξ_L / ρ_L
	Y_R = ξ_R / ρ_R
	J_diff = -λρ_face * (Y_R - Y_L) / d_LR

	# Advective flux with first-order upwinding
	ξ_upwind = v_g >= 0.0 ? ξ_L : ξ_R
	J_adv = v_g * ξ_upwind

	return J_adv + J_diff
end

"""
	gas_flux_advective(v_g::Float64, ξ_L::Float64, ξ_R::Float64) -> Float64

Compute purely advective gas flux at a face [kg/(m²·s)].

Uses first-order upwinding:
	J = ξ_upwind × v_g

# Arguments
- `v_g`: Gas velocity at face [m/s]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
"""
@inline function gas_flux_advective(v_g::Float64, ξ_L::Float64, ξ_R::Float64)
	ξ_upwind = v_g >= 0.0 ? ξ_L : ξ_R
	return v_g * ξ_upwind
end

"""
	compute_gas_fluxes!(J::Matrix{Float64}, mesh::Mesh1D{NC}, 
						cache::PropertyCache, material::Material{NC};
						use_species_scaling::Bool=false) where NC

Compute gas mass fluxes at all faces for all gas components.

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache
- `material`: Material definition
- `use_species_scaling`: If true, scale diffusive flux by √(M_ref/Mⱼ) for each species

Positive flux = flow in +z direction (toward bottom).

# Species-Specific Transport
When `use_species_scaling=true`, light gases (low molar mass) diffuse faster than
heavy gases (high molar mass). This follows kinetic theory where D ∝ 1/√M.

# Performance Note
This function uses compile-time unrolling via `compute_gas_fluxes_unrolled!`
to avoid type instability from heterogeneous tuple indexing.
"""
function compute_gas_fluxes!(J::Matrix{Float64}, mesh::Mesh1D{NC},
	cache::PropertyCache, material::Material{NC};
	use_species_scaling::Bool = false) where NC
	n_faces_val = length(mesh.faces)

	# Initialize fluxes to zero
	fill!(J, 0.0)

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: flux set by boundary conditions
			continue
		end

		# Interior face
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		λρ_face = cache.λρ_face[f]
		ρ_L = cache.ρ[left_id]
		ρ_R = cache.ρ[right_id]

		# Get concentration tuples for this face
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Use compile-time unrolled version for type stability
		if use_species_scaling
			compute_gas_fluxes_species_scaled_unrolled!(J, f, ξ_L, ξ_R, λρ_face, ρ_L, ρ_R, d_LR, material)
		else
			compute_gas_fluxes_unrolled!(J, f, ξ_L, ξ_R, λρ_face, ρ_L, ρ_R, d_LR, material)
		end
	end

	return nothing
end

"""
	compute_gas_fluxes_darcy_fick!(J::Matrix{Float64}, mesh::Mesh1D{NC}, 
								   cache::PropertyCache, material::Material{NC},
								   v_gas::Vector{Float64};
								   use_species_scaling::Bool=false) where NC

Compute gas mass fluxes at all faces using full Darcy-Fick formulation.

This combines advective transport (Darcy flow) with diffusive transport (Fick's law):
	Jⱼ = ξⱼ·v_g - λρ_face·∂Yⱼ/∂x

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache (must have λρ_face, ρ computed)
- `material`: Material definition
- `v_gas`: Darcy gas velocity at faces [m/s] (from compute_darcy_velocity!)
- `use_species_scaling`: If true, scale diffusive flux by √(M_ref/Mⱼ) for each species

Positive flux = flow in +z direction (toward surface for our convention).

# Species-Specific Transport
When `use_species_scaling=true`:
- Diffusive flux is scaled by √(M_ref/Mⱼ) - light gases diffuse faster
- Advective flux is NOT scaled - all species move with bulk gas velocity

# Notes
- When v_gas ≈ 0, this reduces to pure Fickian diffusion (same as compute_gas_fluxes!)
- Uses first-order upwinding for advective flux (stable, conservative)
- Boundary faces are skipped (handled by boundary conditions)
"""
function compute_gas_fluxes_darcy_fick!(J::Matrix{Float64}, mesh::Mesh1D{NC},
	cache::PropertyCache, material::Material{NC},
	v_gas::Vector{Float64};
	use_species_scaling::Bool = false) where NC
	n_faces_val = length(mesh.faces)

	# Initialize fluxes to zero
	fill!(J, 0.0)

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: flux set by boundary conditions
			continue
		end

		# Interior face
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		λρ_face = cache.λρ_face[f]
		ρ_L = cache.ρ[left_id]
		ρ_R = cache.ρ[right_id]
		v_g = v_gas[f]

		# Get concentration tuples for this face
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Use compile-time unrolled version for type stability
		if use_species_scaling
			compute_gas_fluxes_darcy_fick_species_scaled_unrolled!(J, f, v_g, ξ_L, ξ_R,
				λρ_face, ρ_L, ρ_R, d_LR, material)
		else
			compute_gas_fluxes_darcy_fick_unrolled!(J, f, v_g, ξ_L, ξ_R,
				λρ_face, ρ_L, ρ_R, d_LR, material)
		end
	end

	return nothing
end

"""
	surface_mass_flux(mesh::Mesh1D{NC}, cache::PropertyCache, material::Material{NC},
					  bc_set::BoundaryConditionSet, top_face_ids::Vector{Int},
					  t::Float64=0.0) where NC -> Float64

Compute total mass flux leaving the surface using type-stable BoundaryConditionSet.

This is the preferred signature for maximum type stability and zero allocations.
Uses pre-cached top_face_ids for the surface boundary.

# Arguments
- `mesh`: The mesh
- `cache`: Property cache with λ and ρ values
- `material`: Material definition (to identify gas components)
- `bc_set`: Type-stable BoundaryConditionSet containing mass BCs
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `t`: Current time [s] (for time-dependent BCs)

# Returns
Total mass flux [kg/(m²·s)] leaving the surface (positive = gas escaping).
"""
function surface_mass_flux(mesh::Mesh1D{NC}, cache::PropertyCache, material::Material{NC},
	bc_set::BoundaryConditionSet, top_face_ids::Vector{Int},
	t::Float64 = 0.0) where NC

	# Use first top face (all should give same result in 1D)
	if isempty(top_face_ids)
		return 0.0
	end

	f = top_face_ids[1]
	face = mesh.faces[f]
	left_id = face.cell_ids[1]  # Cell adjacent to surface

	if left_id == 0
		# Try other cell_id
		left_id = face.cell_ids[2]
		if left_id == 0
			return 0.0
		end
	end

	# Get surface cell state and properties
	ξ_surface = mesh.cell_states[left_id].ξ
	λ = cache.λ[left_id]
	ρ = cache.ρ[left_id]
	λρ = λ * ρ

	# Distance from cell center to surface
	face_pos = mesh.face_positions[f]
	cell_center = mesh.cell_centers[left_id]
	d = abs(face_pos - cell_center)

	# Get mass BC for top surface using type-stable accessor
	bc = get_mass_bc(bc_set, Val(1))

	# Compute total flux from all gas components
	J_total = 0.0
	for j in 1:NC
		if material.components[j].phase == GAS
			J_total += _compute_surface_flux(bc, ξ_surface[j], λρ, ρ, d, t)
		end
	end

	return J_total
end

# ============================================================================
# Surface Flux Dispatch Functions
# ============================================================================
# Dispatch-based flux calculation for different BC types (inlined here to avoid
# circular dependency with Discretization module).
#
# SIGN CONVENTION NOTE:
# These functions return POSITIVE flux for gas LEAVING the surface. This is the
# natural interpretation for ALE mesh motion: positive flux causes surface regression.
#
# This differs from mass_bc_flux() in boundary_fluxes.jl, which returns NEGATIVE
# flux for gas leaving (because the general BC convention is "positive into domain").
#
# The two conventions are reconciled in _apply_mass_bc_to_flux! (finite_volume.jl),
# which negates the flux at the top boundary: J[f] = -J_bc.
#
# Summary:
# - _compute_surface_flux(): Positive = gas escaping (for ALE regression)
# - mass_bc_flux():          Negative = gas escaping (for flux-into-domain convention)
# ============================================================================

# Surface flux formulation: At the surface boundary, gas escapes to environment.
# The flux is J = λρ (Y_surface - Y_env) / d, where Y_env ≈ 0 and Y = ξ/ρ.
# This gives J = λρ × (ξ/ρ) / d = λξ / d.
# So effectively, using λρ × Y is equivalent to the simpler λξ form.
@inline _compute_surface_flux(::Nothing, ξ, λρ, ρ, d, t) = λρ * (ξ / ρ) / d  # = λξ/d
@inline _compute_surface_flux(::ImpermeableBC, ξ, λρ, ρ, d, t) = 0.0
# Pressure BCs control Darcy velocity, not diffusive flux - return 0
# (The pressure BC contribution to mass flux is handled via Darcy velocity)
@inline _compute_surface_flux(::DirichletPressureBC, ξ, λρ, ρ, d, t) = 0.0
@inline _compute_surface_flux(::NeumannPressureBC, ξ, λρ, ρ, d, t) = 0.0
@inline _compute_surface_flux(::ConvectivePressureBC, ξ, λρ, ρ, d, t) = 0.0

@inline function _compute_surface_flux(bc::DiffusiveMassBC, ξ, λρ, ρ, d, t)
	h_m = bc.h_m
	# The flux uses mass fraction Y = ξ/ρ at the surface
	# J = λρ × Y / d = λρ × (ξ/ρ) / d = λξ/d for h_m=Inf
	# Extract effective λ for combined resistance calculation
	λ = λρ / max(ρ, 1e-20)
	if isinf(h_m)
		# Pure internal diffusion (no external resistance)
		# Positive value: gas leaving surface
		return λ * ξ / d
	else
		# Combined resistance: J = λ × h_m × ξ / (λ + h_m × d)
		# Positive value: gas leaving surface
		return (λ * h_m * ξ) / (λ + h_m * d)
	end
end

@inline function _compute_surface_flux(bc::ConvectiveMassBC, ξ, λρ, ρ, d, t)
	# External convection only: J = h_m × (ξ - ξ_inf)
	# Note: ξ_inf is typically 0 (ambient has no decomposition gas)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m
	ξ_inf = bc.ξ_inf isa Function ? bc.ξ_inf(t) : bc.ξ_inf
	# Positive value: gas leaving surface (when ξ > ξ_inf)
	return h_m * (ξ - ξ_inf)
end

# CombinedMassBC: sum contributions from all sub-BCs
# Use @generated to unroll at compile time and avoid allocation from heterogeneous tuple
@generated function _compute_surface_flux(bc::CombinedMassBC{BCs}, ξ, λρ, ρ, d, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(total += _compute_surface_flux(bc.bcs[$i], ξ, λρ, ρ, d, t)) for i in 1:N]
	return quote
		total = 0.0
		$(exprs...)
		return total
	end
end
