"""
Finite volume operators.

Reference: Design Document Section 8.2

This module provides:
- Divergence operator for 1D finite volumes
- RHS assembly for the complete system of equations
"""

using ..Core: R_GAS, n_components, n_reactions, DirichletThermalBC, RadiationBCSet
using ForwardDiff
using ..Geometry: cell_size, cell_volume, face_area, cell_center, face_position,
	boundary_face_ids, boundary_cell_id, face_area_value
using ..Properties: PropertyCache, update_properties!, effective_permeability
using ..Physics: compute_conductive_fluxes!, compute_radiation_source!,
	compute_gas_fluxes!, evaluate_reactions!, species_source_terms!,
	reaction_heat_source, SourceTermCache

# Helper to safely extract Float64 from Dual during AD passes when storing
# into Float64 caches. Preserves Float64 unchanged.
@inline _value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
@inline _value(x::Real) = x

# ============================================================================
# Divergence Operator
# ============================================================================

"""
	divergence_1d!(div::Vector{Float64}, fluxes::Vector{Float64}, mesh::Mesh1D)

Compute the divergence of a face-centered flux field.

For cell i with faces f_left and f_right:
	div[i] = (1/V_i) * (F[f_right] - F[f_left]) * A

Sign convention:
- Positive flux = flow in +z direction (from left to right / top to bottom)
- For a cell, right face flux is outflow (+), left face flux is inflow (-)

# Arguments
- `div`: Output vector of divergence values [units of flux / length]
- `fluxes`: Face-centered flux values [flux units]
- `mesh`: 1D mesh with geometry
"""
function divergence_1d!(div::Vector{Float64}, fluxes::Vector{Float64}, mesh)
	n_cells_val = length(mesh.cells)

	@inbounds for i in 1:n_cells_val
		cell = mesh.cells[i]
		f_left, f_right = cell.face_ids

		V_i = mesh.cell_volumes[i]
		A = mesh.face_areas[f_left]  # Same for both faces in 1D

		# Net flux into cell: inflow (left) - outflow (right)
		# For heat/mass conservation: dφ/dt = -∇·F = -(F_out - F_in)/V
		div[i] = (fluxes[f_right] - fluxes[f_left]) * A / V_i
	end

	return nothing
end

"""
	divergence_species_1d!(div::Matrix{Float64}, J::Matrix{Float64}, mesh::Mesh1D)

Compute divergence for species fluxes (multiple components).

# Arguments  
- `div`: Output matrix [NC × n_cells] of divergence values
- `J`: Face-centered species fluxes [NC × n_faces]
- `mesh`: 1D mesh
"""
function divergence_species_1d!(div::Matrix{Float64}, J::Matrix{Float64}, mesh)
	NC, n_cells_val = size(div)

	@inbounds for i in 1:n_cells_val
		cell = mesh.cells[i]
		f_left, f_right = cell.face_ids

		V_i = mesh.cell_volumes[i]
		A = mesh.face_areas[f_left]

		for j in 1:NC
			div[j, i] = (J[j, f_right] - J[j, f_left]) * A / V_i
		end
	end

	return nothing
end

# ============================================================================
# RHS Assembly for Fixed Mesh (Eulerian)
# ============================================================================

"""
	RHSCache

Pre-allocated workspace for RHS computation.
"""
mutable struct RHSCache
	# Work arrays for fluxes
	q_cond::Vector{Float64}      # Conductive heat fluxes at faces
	q_total::Vector{Float64}     # Total energy fluxes at faces
	q_conv::Vector{Float64}      # Convective heat fluxes at faces (advection form: J × cp × ΔT)
	J_gas::Matrix{Float64}       # Gas species fluxes [NC × n_faces]

	# Work arrays for sources
	Q_rxn::Vector{Float64}       # Reaction heat source per cell
	Q_rad::Vector{Float64}       # Radiation source per cell
	dξ_rxn::Matrix{Float64}      # Species source from reactions [NC × n_cells]

	# Per-cell reaction workspace
	rates::Vector{Float64}       # Reaction rates (NR)
	heats::Vector{Float64}       # Heats of reaction (NR)
	species_src::Vector{Float64} # Species sources (NC)

	# Divergence workspace
	div_q::Vector{Float64}       # Divergence of heat flux
	div_J::Matrix{Float64}       # Divergence of species flux [NC × n_cells]
end

"""
	RHSCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Create an RHSCache with pre-allocated arrays.
"""
function RHSCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)
	RHSCache(
		zeros(n_faces),           # q_cond
		zeros(n_faces),           # q_total
		zeros(n_faces),           # q_conv
		zeros(NC, n_faces),       # J_gas
		zeros(n_cells),           # Q_rxn
		zeros(n_cells),           # Q_rad
		zeros(NC, n_cells),       # dξ_rxn
		zeros(NR),                # rates
		zeros(NR),                # heats
		zeros(NC),                # species_src
		zeros(n_cells),           # div_q
		zeros(NC, n_cells),       # div_J
	)
end

"""
	RHSCache(mesh, material)

Create an RHSCache sized for the given mesh and material.
"""
function RHSCache(mesh, material)
	n_cells_val = length(mesh.cells)
	n_faces_val = length(mesh.faces)
	NC = n_components(material)
	NR = n_reactions(material)
	return RHSCache(n_cells_val, n_faces_val, NC, NR)
end

"""
	resize_rhs_cache!(cache::RHSCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Resize RHS cache arrays (e.g., after mesh refinement).
"""
function resize_rhs_cache!(cache::RHSCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)
	resize!(cache.q_cond, n_faces)
	resize!(cache.q_total, n_faces)
	resize!(cache.q_conv, n_faces)
	cache.J_gas = zeros(NC, n_faces)
	resize!(cache.Q_rxn, n_cells)
	resize!(cache.Q_rad, n_cells)
	cache.dξ_rxn = zeros(NC, n_cells)
	resize!(cache.rates, NR)
	resize!(cache.heats, NR)
	resize!(cache.species_src, NC)
	resize!(cache.div_q, n_cells)
	cache.div_J = zeros(NC, n_cells)
	return cache
end

"""
	compute_rhs!(dT::Vector{Float64}, dξ::Matrix{Float64},
				 mesh::Mesh1D{NC}, prop_cache::PropertyCache, rhs_cache::RHSCache,
				 material::Material{NC}, bc_set::BoundaryConditionSet, t,
				 top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
				 I_radiation::Float64=0.0) where NC

Compute the right-hand side of the governing equations using type-stable BCs.

This version uses pre-cached face lists and static dispatch for maximum
performance. All boundary condition handling is specialized at compile time.

# Arguments
- `dT`: Output vector of dT/dt [K/s]
- `dξ`: Output matrix of dξⱼ/dt [kg/(m³·s)], size [NC × n_cells]
- `mesh`: Mesh with current state (T, ξ in cell_states)
- `prop_cache`: Property cache (must be updated before calling)
- `rhs_cache`: RHS workspace cache
- `material`: Material definition
- `bc_set`: Type-stable BoundaryConditionSet with all boundary conditions
- `t`: Current time [s]
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `I_radiation`: Surface radiation intensity [W/m²] (default 0)

# Equations solved
Energy: ρcₚ dT/dt = -∇·q_cond + Q_rad + Q_rxn
Species: dξⱼ/dt = -∇·Jⱼ + Sⱼ (for gases with transport)
		 dξⱼ/dt = Sⱼ          (for solids, no transport)

The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations in the
inner loops.
"""
function compute_rhs!(dT::Vector{Float64}, dξ::Matrix{Float64},
	mesh::Mesh1D{NC}, prop_cache::PropertyCache, rhs_cache::RHSCache,
	material::Material{NC}, bc_set::BoundaryConditionSet, t,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
	I_radiation::Float64 = 0.0,
	radiation_model::RadiationModel = NO_RADIATION) where NC

	n_cells_val = length(mesh.cells)
	NR = n_reactions(material)

	# Ensure caches are properly sized
	if length(rhs_cache.q_cond) != length(mesh.faces)
		resize_rhs_cache!(rhs_cache, n_cells_val, length(mesh.faces), NC, NR)
	end

	# ========================================================================
	# Step 1: Compute conductive heat fluxes at interior faces
	# ========================================================================
	compute_conductive_fluxes!(rhs_cache.q_cond, mesh, prop_cache)

	# Copy to total flux (will add boundary contributions)
	copyto!(rhs_cache.q_total, rhs_cache.q_cond)

	# ========================================================================
	# Step 2: Apply thermal boundary conditions (static dispatch)
	# ========================================================================
	apply_thermal_bc_flux!(rhs_cache.q_total, mesh, bc_set, prop_cache, t,
		top_face_ids, bottom_face_ids, radiation_model)

	# ========================================================================
	# Step 3: Compute radiation source (model-dependent)
	# ========================================================================
	compute_radiation_source!(rhs_cache.Q_rad, mesh, prop_cache, I_radiation, radiation_model)

	# ========================================================================
	# Step 4: Compute reaction rates and source terms for each cell
	# ========================================================================
	@inbounds for i in 1:n_cells_val
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		# Evaluate reaction rates and heats
		evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)

		# Compute reaction heat source
		rhs_cache.Q_rxn[i] = reaction_heat_source(rhs_cache.rates, rhs_cache.heats)

		# Compute species source terms
		species_source_terms!(rhs_cache.species_src, rhs_cache.rates, material, NC)

		# Copy to output matrix
		for j in 1:NC
			rhs_cache.dξ_rxn[j, i] = rhs_cache.species_src[j]
		end
	end

	# ========================================================================
	# Step 5: Compute gas transport fluxes
	# ========================================================================
	compute_gas_fluxes!(rhs_cache.J_gas, mesh, prop_cache, material)

	# Apply mass boundary conditions (static dispatch)
	apply_mass_bc_flux!(rhs_cache.J_gas, mesh, prop_cache, bc_set, t,
		top_face_ids, bottom_face_ids, Val(NC))

	# ========================================================================
	# Step 5.5: Compute convective heat flux from gas mass fluxes
	# ========================================================================
	# Advection form: q_conv = Σⱼ Jⱼ × cp,j × (T_R - T_L)
	# This redistributes heat internally as gas flows through temperature gradients
	# Boundary faces have q_conv = 0 (no enthalpy flux at boundaries)
	Physics.compute_convective_heat_fluxes_from_mass_fluxes!(rhs_cache.q_conv,
		rhs_cache.J_gas, mesh, material)

	# Add convective flux to total heat flux
	n_faces_val = length(mesh.faces)
	@inbounds for f in 1:n_faces_val
		rhs_cache.q_total[f] += rhs_cache.q_conv[f]
	end

	# ========================================================================
	# Step 6: Compute divergences
	# ========================================================================
	divergence_1d!(rhs_cache.div_q, rhs_cache.q_total, mesh)
	divergence_species_1d!(rhs_cache.div_J, rhs_cache.J_gas, mesh)

	# ========================================================================
	# Step 7: Assemble final RHS
	# ========================================================================
	@inbounds for i in 1:n_cells_val
		# Energy equation: dT/dt = (1/ρcₚ) * (-∇·q + Q_rad + Q_rxn)
		ρcp = prop_cache.ρcp[i]
		dT[i] = (-rhs_cache.div_q[i] + rhs_cache.Q_rad[i] + rhs_cache.Q_rxn[i]) / ρcp

		# Species equations - use unrolled version for type stability
		apply_species_residual_unrolled!(dξ, i, rhs_cache.div_J, rhs_cache.dξ_rxn, material)
	end

	return nothing
end

# ============================================================================
# Helper functions for boundary condition application
# ============================================================================

using ..Core: BoundaryConditionSet, get_thermal_bc, get_mass_bc
using ..Core: RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY

"""
	apply_thermal_bc_flux!(q::Vector{Float64}, mesh, 
						   bc_set::BoundaryConditionSet, 
						   prop_cache, t,
						   top_face_ids::Vector{Int}, 
						   bottom_face_ids::Vector{Int},
						   radiation_model::RadiationModel = SURFACE_ABSORPTION)

Apply thermal boundary conditions using static dispatch (type-stable).

This version uses pre-cached face lists and BoundaryConditionSet for maximum
performance. Each boundary condition type is dispatched at compile time.

# Arguments
- `q`: Face flux array (modified in-place)
- `mesh`: The mesh
- `bc_set`: Type-stable boundary condition set
- `prop_cache`: Property cache
- `t`: Current time
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `radiation_model`: Radiation model (determines how RadiativeFluxBC handles absorption)
"""
function apply_thermal_bc_flux!(q::Vector{Float64}, mesh,
	bc_set::BoundaryConditionSet,
	prop_cache, t,
	top_face_ids::Vector{Int},
	bottom_face_ids::Vector{Int},
	radiation_model::RadiationModel = SURFACE_ABSORPTION)
	# Top boundary (boundary_id = 1, z = L)
	_apply_thermal_bc_static!(q, mesh, bc_set.thermal_top, prop_cache, t,
		top_face_ids, Val(:top), radiation_model)

	# Bottom boundary (boundary_id = 2, z = 0)
	_apply_thermal_bc_static!(q, mesh, bc_set.thermal_bottom, prop_cache, t,
		bottom_face_ids, Val(:bottom), radiation_model)

	return nothing
end

"""
Inner function for type-stable thermal BC application.

The `Val{boundary}` parameter enables compile-time dispatch on boundary type.
The `BC` type parameter enables compile-time dispatch on BC type.
The `radiation_model` determines how RadiativeFluxBC handles absorption.

For non-Dirichlet BCs, this now uses the boundary temperature solve to find
the surface temperature T_s that satisfies energy balance:
	q_BC(T_s) = k_eff × (T_cell - T_s) / Δz_half

This ensures proper energy conservation at boundaries, especially important
for radiative BCs where flux depends on T⁴.
"""
@inline function _apply_thermal_bc_static!(q::Vector{Float64}, mesh,
	bc::BC, prop_cache, t,
	face_ids::Vector{Int},
	::Val{boundary},
	radiation_model::RadiationModel = SURFACE_ABSORPTION) where {BC <: AbstractThermalBC, boundary}
	# Handle DirichletThermalBC specially (needs k and distance)
	if BC <: DirichletThermalBC
		@inbounds for f in face_ids
			face = mesh.faces[f]
			left_id, right_id = face.cell_ids
			cell_id = left_id == 0 ? right_id : left_id
			T_surface = mesh.cell_states[cell_id].T

			# Get prescribed boundary temperature
			T_bc = evaluate_param(bc.T, t)

			# Get cell thermal conductivity
			k_cell = prop_cache.k[cell_id]

			# Distance from cell center to boundary face (half cell width)
			dz = mesh.cell_volumes[cell_id] / mesh.cross_section_area
			dz_half = dz / 2.0

			# Conductive flux: q = -k * dT/dz
			q_bc = k_cell * (T_bc - T_surface) / dz_half

			# Apply sign convention based on boundary
			if boundary === :bottom
				q[f] = _value(q_bc)
			else  # :top
				q[f] = -_value(q_bc)
			end
		end
	else
		# ================================================================
		# Use boundary temperature solve for energy-conserving flux
		# ================================================================
		# This solves for T_surface such that q_BC(T_s) = k(T_cell - T_s)/Δz_half
		# ensuring perfect energy balance at the boundary.

		# Determine boundary_id from boundary symbol
		boundary_id = boundary === :top ? 1 : 2

		# Update boundary state (solves for T_surface satisfying energy balance)
		# Extract Float64 value from t - the solver doesn't need time derivatives
		t_val = _value(t)
		state = update_boundary_state!(mesh, boundary_id, bc, prop_cache, t_val, radiation_model)

		# Apply the conducted flux to all faces on this boundary
		# q_conducted is positive when heat flows into the domain
		@inbounds for f in face_ids
			if boundary === :bottom
				q[f] = _value(state.q_conducted)
			else  # :top
				q[f] = -_value(state.q_conducted)
			end
		end
	end
	return nothing
end

"""
	apply_mass_bc_flux!(J::Matrix{Float64}, mesh, prop_cache,
						bc_set::BoundaryConditionSet, t,
						top_face_ids::Vector{Int}, 
						bottom_face_ids::Vector{Int},
						::Val{NC}) where NC

Apply mass boundary conditions using static dispatch (type-stable).

# Arguments
- `J`: Species flux matrix [NC × n_faces] (modified in-place)
- `mesh`: The mesh
- `prop_cache`: Property cache (for λ, ρ)
- `bc_set`: Type-stable boundary condition set
- `t`: Current time
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `Val{NC}`: Number of components (for type stability)
"""
function apply_mass_bc_flux!(J::Matrix{Float64}, mesh, prop_cache,
	bc_set::BoundaryConditionSet, t,
	top_face_ids::Vector{Int},
	bottom_face_ids::Vector{Int},
	::Val{NC}) where NC
	# Top boundary (boundary_id = 1, z = L)
	_apply_mass_bc_static!(J, mesh, bc_set.mass_top, prop_cache, t,
		top_face_ids, Val(:top), Val(NC))

	# Bottom boundary (boundary_id = 2, z = 0)
	_apply_mass_bc_static!(J, mesh, bc_set.mass_bottom, prop_cache, t,
		bottom_face_ids, Val(:bottom), Val(NC))

	return nothing
end

"""
Inner function for type-stable mass BC application.
"""
@inline function _apply_mass_bc_static!(J::Matrix{Float64}, mesh,
	bc::BC, prop_cache, t,
	face_ids::Vector{Int},
	::Val{boundary},
	::Val{NC}) where {BC <: AbstractMassBC, boundary, NC}
	@inbounds for f in face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids
		cell_id = left_id == 0 ? right_id : left_id
		ξ_surface = mesh.cell_states[cell_id].ξ

		# Get λρ and ρ for diffusive BC (rigorous formulation)
		λ = prop_cache.λ[cell_id]
		ρ = prop_cache.ρ[cell_id]
		λρ = λ * ρ  # Combined transport coefficient
		face_pos = mesh.face_positions[f]
		cell_center = mesh.cell_centers[cell_id]
		d = abs(face_pos - cell_center)

		for j in 1:NC
			J_bc = mass_bc_flux(bc, ξ_surface[j], λρ, ρ, d, t)

			# Sign convention based on boundary
			if boundary === :bottom
				J[j, f] = _value(J_bc)   # Mass flux at bottom
			else  # :top
				J[j, f] = -_value(J_bc)  # Mass flux at top (gases escape here)
			end
		end
	end
	return nothing
end

# ============================================================================
# Simplified RHS for testing (no boundary conditions)
# ============================================================================

"""
	compute_rhs_interior!(dT::Vector{Float64}, dξ::Matrix{Float64},
						  mesh::Mesh1D{NC}, prop_cache::PropertyCache, rhs_cache::RHSCache,
						  material::Material{NC}) where NC

Compute RHS for interior cells only (no boundary condition handling).
Useful for testing and verification.

The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations.
"""
function compute_rhs_interior!(dT::Vector{Float64}, dξ::Matrix{Float64},
	mesh::Mesh1D{NC}, prop_cache::PropertyCache, rhs_cache::RHSCache,
	material::Material{NC}) where NC
	n_cells_val = length(mesh.cells)

	# Compute conductive fluxes
	compute_conductive_fluxes!(rhs_cache.q_cond, mesh, prop_cache)

	# Compute reaction sources for each cell
	@inbounds for i in 1:n_cells_val
		state = mesh.cell_states[i]
		evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, state.ξ, state.T, material)
		rhs_cache.Q_rxn[i] = reaction_heat_source(rhs_cache.rates, rhs_cache.heats)
		species_source_terms!(rhs_cache.species_src, rhs_cache.rates, material, NC)
		for j in 1:NC
			rhs_cache.dξ_rxn[j, i] = rhs_cache.species_src[j]
		end
	end

	# Compute gas fluxes
	compute_gas_fluxes!(rhs_cache.J_gas, mesh, prop_cache, material)

	# Compute divergences
	divergence_1d!(rhs_cache.div_q, rhs_cache.q_cond, mesh)
	divergence_species_1d!(rhs_cache.div_J, rhs_cache.J_gas, mesh)

	# Assemble RHS
	@inbounds for i in 1:n_cells_val
		ρcp = prop_cache.ρcp[i]
		dT[i] = (-rhs_cache.div_q[i] + rhs_cache.Q_rxn[i]) / ρcp

		for j in 1:NC
			if material.components[j].phase == GAS
				dξ[j, i] = -rhs_cache.div_J[j, i] + rhs_cache.dξ_rxn[j, i]
			else
				dξ[j, i] = rhs_cache.dξ_rxn[j, i]
			end
		end
	end

	return nothing
end

# ============================================================================
# Darcy-Fick RHS Computation
# ============================================================================

# NOTE: The old Dict-based and perm_model-based `update_darcy_properties!` have been removed
# in favor of the type-stable API using component-based permeability.
# The current API: `update_darcy_properties!(cache, mesh, material, bc_set, top_face_ids, bottom_face_ids)`
# computes effective permeability from component κ values via mixing rules.

"""
	update_darcy_properties!(cache::PropertyCache, mesh::Mesh1D{NC},
							  material::Material{NC},
							  bc_set::BoundaryConditionSet,
							  top_face_ids::Vector{Int},
							  bottom_face_ids::Vector{Int}) where NC

Update Darcy-related properties using component-based permeability.

This method computes effective permeability from component κ values via mixing rules,
rather than using a porosity-based permeability model. Use this when:
- Components have κ values specified (e.g., SolidComponent(:char, ..., κ=1e-12))
- No legacy permeability_model is provided

The effective permeability is computed as a volume-weighted mixture of
component permeabilities, similar to thermal conductivity.

# Arguments
- `cache`: PropertyCache (modified in place)
- `mesh`: Mesh with current cell states
- `material`: Material definition with component κ values and permeability_mixing
- `bc_set`: Type-stable BoundaryConditionSet containing mass BCs
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
"""
function update_darcy_properties!(cache::PropertyCache, mesh::Mesh1D{NC},
	material::Material{NC},
	bc_set::BoundaryConditionSet,
	top_face_ids::Vector{Int},
	bottom_face_ids::Vector{Int}) where NC
	n_cells = length(mesh.cells)
	n_faces = length(mesh.faces)

	# Step 1: Compute cell-centered pressure, porosity, permeability, viscosity
	@inbounds for i in 1:n_cells
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		# Pressure from ideal gas law
		cache.P[i] = Physics.ideal_gas_pressure(ξ, T, material)

		# Porosity (gas volume fraction)
		cache.φ[i] = _compute_porosity(ξ, T, material)

		# Permeability from component mixture (NEW: component-based)
		cache.κ[i] = effective_permeability(ξ, T, material)

		# Gas viscosity from Sutherland's law
		cache.μ[i] = Properties.effective_gas_viscosity(T)
	end

	# Step 2: Compute face-centered pressure gradient and Darcy velocity
	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: zero gradient (will be overwritten by BCs)
			cache.∇P[f] = 0.0
			cache.v_gas[f] = 0.0
			continue
		end

		# Distances from cell centers to face (for distance-weighted harmonic mean)
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R

		# Pressure gradient
		cache.∇P[f] = (cache.P[right_id] - cache.P[left_id]) / d_LR

		# Harmonic mean of mobility κ/μ at face
		κ_L = cache.κ[left_id]
		κ_R = cache.κ[right_id]
		μ_L = cache.μ[left_id]
		μ_R = cache.μ[right_id]

		# Mobility = κ/μ
		mob_L = μ_L > 1e-20 ? κ_L / μ_L : 0.0
		mob_R = μ_R > 1e-20 ? κ_R / μ_R : 0.0

		# Distance-weighted harmonic mean for series resistance (non-uniform meshes)
		# For uniform meshes (d_L = d_R), this reduces to standard harmonic mean
		if mob_L > 1e-30 && mob_R > 1e-30
			mob_face = (d_L + d_R) / (d_L / mob_L + d_R / mob_R)
		elseif mob_L > 1e-30
			mob_face = mob_L
		elseif mob_R > 1e-30
			mob_face = mob_R
		else
			mob_face = 0.0
		end

		# Darcy velocity: v_g = -κ/μ × ∇P
		cache.v_gas[f] = -mob_face * cache.∇P[f]
	end

	# Step 3: Apply pressure boundary conditions using type-stable bc_set
	apply_pressure_bc_to_darcy_velocity!(cache.v_gas, mesh, cache.P, cache.κ,
		bc_set, top_face_ids, bottom_face_ids)

	return cache
end

"""
	compute_rhs_darcy_fick!(dT::Vector{Float64}, dξ::Matrix{Float64},
							mesh::Mesh1D{NC}, prop_cache::PropertyCache,
							rhs_cache::RHSCache, material::Material{NC},
							bc_set::BoundaryConditionSet, t,
							top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
							I_radiation::Float64=0.0,
							radiation_model::RadiationModel=NO_RADIATION) where NC

Compute RHS with full Darcy-Fick gas transport using component-based permeability.

This version uses component κ values and mixing rules to compute effective
permeability, rather than a separate permeability model. Use this when:
- Components have κ values specified (e.g., SolidComponent(:char, ..., κ=1e-12))
- No legacy permeability_model is provided

# Arguments
- `dT`: Output vector of dT/dt [K/s]
- `dξ`: Output matrix of dξⱼ/dt [kg/(m³·s)], size [NC × n_cells]
- `mesh`: Mesh with current state
- `prop_cache`: Property cache (will be updated with Darcy properties)
- `rhs_cache`: RHS workspace cache (includes q_conv for convective heat fluxes)
- `material`: Material definition with component κ values
- `bc_set`: Type-stable BoundaryConditionSet
- `t`: Current time [s]
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `I_radiation`: Surface radiation intensity [W/m²]
- `radiation_model`: Radiation model type

# Notes
This is the preferred API for Darcy-Fick transport. The effective permeability
is computed from component values via mixing rules (PARALLEL, SERIES, WEIGHTED).
Convective heat fluxes are stored in rhs_cache.q_conv.
"""
function compute_rhs_darcy_fick!(dT::Vector{Float64}, dξ::Matrix{Float64},
	mesh::Mesh1D{NC}, prop_cache::PropertyCache,
	rhs_cache::RHSCache, material::Material{NC},
	bc_set::BoundaryConditionSet, t,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
	I_radiation::Float64 = 0.0,
	radiation_model::RadiationModel = NO_RADIATION) where NC

	n_cells_val = length(mesh.cells)
	n_faces_val = length(mesh.faces)
	NR = n_reactions(material)

	# Ensure caches are properly sized
	if length(rhs_cache.q_cond) != n_faces_val
		resize_rhs_cache!(rhs_cache, n_cells_val, n_faces_val, NC, NR)
	end

	# ========================================================================
	# Step 1: Update Darcy properties using component-based permeability
	# ========================================================================
	update_darcy_properties!(prop_cache, mesh, material,
		bc_set, top_face_ids, bottom_face_ids)

	# ========================================================================
	# Step 2: Compute conductive heat fluxes at interior faces
	# ========================================================================
	compute_conductive_fluxes!(rhs_cache.q_cond, mesh, prop_cache)

	# Copy to total flux (will add boundary and convective contributions)
	copyto!(rhs_cache.q_total, rhs_cache.q_cond)

	# ========================================================================
	# Step 3: Apply thermal boundary conditions (static dispatch)
	# ========================================================================
	apply_thermal_bc_flux!(rhs_cache.q_total, mesh, bc_set, prop_cache, t,
		top_face_ids, bottom_face_ids, radiation_model)

	# ========================================================================
	# Step 4: (Convective heat flux computed after gas fluxes in Step 8)
	# ========================================================================

	# ========================================================================
	# Step 5: Compute radiation source (model-dependent)
	# ========================================================================
	compute_radiation_source!(rhs_cache.Q_rad, mesh, prop_cache, I_radiation, radiation_model)

	# ========================================================================
	# Step 6: Compute reaction rates and source terms for each cell
	# ========================================================================
	@inbounds for i in 1:n_cells_val
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		# Evaluate reaction rates and heats
		evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)

		# Compute reaction heat source
		rhs_cache.Q_rxn[i] = reaction_heat_source(rhs_cache.rates, rhs_cache.heats)

		# Compute species source terms
		species_source_terms!(rhs_cache.species_src, rhs_cache.rates, material, NC)

		# Copy to output matrix
		for j in 1:NC
			rhs_cache.dξ_rxn[j, i] = rhs_cache.species_src[j]
		end
	end

	# ========================================================================
	# Step 7: Compute gas transport fluxes using Darcy-Fick
	# ========================================================================
	Physics.compute_gas_fluxes_darcy_fick!(rhs_cache.J_gas, mesh, prop_cache,
		material, prop_cache.v_gas)

	# Apply mass boundary conditions (static dispatch)
	apply_mass_bc_flux!(rhs_cache.J_gas, mesh, prop_cache, bc_set, t,
		top_face_ids, bottom_face_ids, Val(NC))

	# ========================================================================
	# Step 8: Compute convective heat flux from gas mass fluxes
	# ========================================================================
	Physics.compute_convective_heat_fluxes_from_mass_fluxes!(rhs_cache.q_conv, rhs_cache.J_gas, mesh, material)

	# Add convective flux to total heat flux
	@inbounds for f in 1:n_faces_val
		rhs_cache.q_total[f] += rhs_cache.q_conv[f]
	end

	# ========================================================================
	# Step 9: Compute divergences
	# ========================================================================
	divergence_1d!(rhs_cache.div_q, rhs_cache.q_total, mesh)
	divergence_species_1d!(rhs_cache.div_J, rhs_cache.J_gas, mesh)

	# ========================================================================
	# Step 10: Assemble final RHS
	# ========================================================================
	@inbounds for i in 1:n_cells_val
		ρcp = prop_cache.ρcp[i]

		dT[i] = (-rhs_cache.div_q[i] + rhs_cache.Q_rad[i] + rhs_cache.Q_rxn[i]) / ρcp

		# Species equations
		apply_species_residual_unrolled!(dξ, i, rhs_cache.div_J, rhs_cache.dξ_rxn, material)
	end

	return nothing
end

"""
Internal helper to compute porosity (gas volume fraction).
Uses dispatch-based unrolling to avoid allocation from heterogeneous tuple indexing.
"""
@inline function _compute_porosity(ξ::NTuple{NC, Float64}, T::Float64,
	material::Material{NC}) where NC
	return _compute_porosity_impl(ξ, T, material.components, Val(NC))
end

# Specialized versions for common component counts
@inline function _compute_porosity_impl(ξ::NTuple{1, Float64}, T::Float64,
	comps::Tuple{C1}, ::Val{1}) where {C1}
	solid_vol_frac = 0.0
	c1 = comps[1]
	if c1.phase == SOLID || c1.phase == LIQUID
		ρ1 = c1.density(T)
		if ρ1 > 0.0
			solid_vol_frac += ξ[1] / ρ1
		end
	end
	return clamp(1.0 - solid_vol_frac, 0.0, 1.0)
end

@inline function _compute_porosity_impl(ξ::NTuple{2, Float64}, T::Float64,
	comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	solid_vol_frac = 0.0
	c1 = comps[1]
	if c1.phase == SOLID || c1.phase == LIQUID
		ρ1 = c1.density(T)
		if ρ1 > 0.0
			solid_vol_frac += ξ[1] / ρ1
		end
	end
	c2 = comps[2]
	if c2.phase == SOLID || c2.phase == LIQUID
		ρ2 = c2.density(T)
		if ρ2 > 0.0
			solid_vol_frac += ξ[2] / ρ2
		end
	end
	return clamp(1.0 - solid_vol_frac, 0.0, 1.0)
end

@inline function _compute_porosity_impl(ξ::NTuple{3, Float64}, T::Float64,
	comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	solid_vol_frac = 0.0
	c1 = comps[1]
	if c1.phase == SOLID || c1.phase == LIQUID
		ρ1 = c1.density(T)
		if ρ1 > 0.0
			solid_vol_frac += ξ[1] / ρ1
		end
	end
	c2 = comps[2]
	if c2.phase == SOLID || c2.phase == LIQUID
		ρ2 = c2.density(T)
		if ρ2 > 0.0
			solid_vol_frac += ξ[2] / ρ2
		end
	end
	c3 = comps[3]
	if c3.phase == SOLID || c3.phase == LIQUID
		ρ3 = c3.density(T)
		if ρ3 > 0.0
			solid_vol_frac += ξ[3] / ρ3
		end
	end
	return clamp(1.0 - solid_vol_frac, 0.0, 1.0)
end

@inline function _compute_porosity_impl(ξ::NTuple{4, Float64}, T::Float64,
	comps::Tuple{C1, C2, C3, C4}, ::Val{4}) where {C1, C2, C3, C4}
	solid_vol_frac = 0.0
	c1 = comps[1]
	if c1.phase == SOLID || c1.phase == LIQUID
		ρ1 = c1.density(T)
		if ρ1 > 0.0
			solid_vol_frac += ξ[1] / ρ1
		end
	end
	c2 = comps[2]
	if c2.phase == SOLID || c2.phase == LIQUID
		ρ2 = c2.density(T)
		if ρ2 > 0.0
			solid_vol_frac += ξ[2] / ρ2
		end
	end
	c3 = comps[3]
	if c3.phase == SOLID || c3.phase == LIQUID
		ρ3 = c3.density(T)
		if ρ3 > 0.0
			solid_vol_frac += ξ[3] / ρ3
		end
	end
	c4 = comps[4]
	if c4.phase == SOLID || c4.phase == LIQUID
		ρ4 = c4.density(T)
		if ρ4 > 0.0
			solid_vol_frac += ξ[4] / ρ4
		end
	end
	return clamp(1.0 - solid_vol_frac, 0.0, 1.0)
end

# Generic fallback using generated function for NC > 4
@inline function _compute_porosity_impl(ξ::NTuple{NC, Float64}, T::Float64,
	comps::Tuple, ::Val{NC}) where NC
	return _compute_porosity_gen(ξ, T, comps)
end

@generated function _compute_porosity_gen(ξ::NTuple{NC, Float64}, T::Float64,
	comps::Tuple) where NC
	expr = quote
		solid_vol_frac = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			if c.phase == SOLID || c.phase == LIQUID
				ρ = c.density(T)
				if ρ > 0.0
					solid_vol_frac += ξ[$i] / ρ
				end
			end
		end)
	end
	push!(expr.args, :(return clamp(1.0 - solid_vol_frac, 0.0, 1.0)))
	return expr
end

# ============================================================================
# RHS computation with P1 radiation model
# ============================================================================

"""
	compute_rhs_with_p1_radiation!(dT::Vector{Float64}, dξ::Matrix{Float64},
								  mesh::Mesh1D{NC}, prop_cache::PropertyCache,
								  rhs_cache::RHSCache, material::Material{NC},
								  bc_set::BoundaryConditionSet, t::Float64,
								  top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
								  rad_cache::RadiationCache,
								  radiation_bcs::RadiationBCSet;
								  I_radiation::Float64=0.0,
								  use_quasi_steady::Bool=true) where NC

Compute RHS for energy and species equations with P1 volumetric radiation model.

This function extends `compute_rhs!` to include P1 radiative transport:
1. Computes standard RHS (conduction, reactions, gas transport)
2. Solves quasi-steady P1 equation if `use_quasi_steady=true`
3. Adds radiative flux divergence ∇·q_rad = αG - 4σεαT⁴ to energy equation

# Arguments
- `dT`: Output temperature rate of change [K/s]
- `dξ`: Output species concentration rate of change [kg/(m³·s)]
- `mesh`: 1D mesh with current state
- `prop_cache`: Cached material properties
- `rhs_cache`: Pre-allocated work arrays
- `material`: Material definition
- `bc_set`: Type-stable BoundaryConditionSet
- `t`: Current time [s]
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `rad_cache`: Radiation cache with G field and work arrays
- `radiation_bcs`: RadiationBCSet with radiation boundary conditions
- `I_radiation`: External radiation intensity [W/m²] (for Beer-Lambert, typically 0 with P1)
- `use_quasi_steady`: If true, solve quasi-steady P1; if false, G is already updated

# Notes
- The radiative flux divergence replaces/augments Beer-Lambert absorption
- For quasi-steady mode, G is solved iteratively at each timestep
- For transient mode, G should be integrated separately as part of the state
"""
function compute_rhs_with_p1_radiation!(dT::Vector{Float64}, dξ::Matrix{Float64},
	mesh::Mesh1D{NC}, prop_cache::PropertyCache,
	rhs_cache::RHSCache, material::Material{NC},
	bc_set::BoundaryConditionSet, t::Float64,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
	rad_cache::RadiationCache,
	radiation_bcs::RadiationBCSet;
	I_radiation::Float64 = 0.0,
	use_quasi_steady::Bool = true) where {NC}

	n_cells_val = mesh.n_active_cells

	# ========================================================================
	# Step 1: Solve quasi-steady P1 radiation equation (if enabled)
	# ========================================================================
	if use_quasi_steady
		# Update radiation properties in rad_cache from current temperature
		update_radiation_properties!(rad_cache, mesh, prop_cache)

		# Solve ∇·(D_rad∇G) = αG - 4σεαT⁴ with Newton iteration
		solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, t)
	end

	# ========================================================================
	# Step 2: Compute P1 radiative flux divergence
	# ========================================================================
	# ∇·q_rad = αG - 4σεαT⁴ (absorption - emission)
	# This is stored in rad_cache.div_q_rad
	compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

	# ========================================================================
	# Step 3: Compute standard RHS (conduction, reactions, transport)
	# ========================================================================
	# Use base compute_rhs! with BoundaryConditionSet (I_radiation=0 since P1 handles radiation)
	compute_rhs!(dT, dξ, mesh, prop_cache, rhs_cache, material,
		bc_set, t, top_face_ids, bottom_face_ids; I_radiation = I_radiation)

	# ========================================================================
	# Step 4: Add P1 radiation contribution to energy equation
	# ========================================================================
	# The compute_rhs! already added Beer-Lambert Q_rad (which should be 0)
	# Now add P1 radiation: dT += (-div_q_rad) / ρcₚ
	# Note: div_q_rad = αG - 4σεαT⁴ (positive = net absorption)
	# Energy equation: ∂T/∂t = -∇·q / ρcₚ, so dT -= div_q_rad / ρcₚ
	@inbounds for i in 1:n_cells_val
		ρcp = prop_cache.ρcp[i]
		# Subtract because div_q_rad is the divergence of radiative heat flux
		# and energy equation has -∇·q_total
		dT[i] -= rad_cache.div_q_rad[i] / ρcp
	end

	return nothing
end

"""
	update_radiation_properties!(rad_cache::RadiationCache, mesh::Mesh1D, prop_cache::PropertyCache)

Update radiation properties (α, D_rad) in the radiation cache from current state.

Called before solving the P1 equation to ensure properties match current temperature.
"""
function update_radiation_properties!(rad_cache::RadiationCache, mesh::Mesh1D, prop_cache::PropertyCache)
	n_cells = mesh.n_active_cells
	n_faces = n_cells + 1

	@inbounds for i in 1:n_cells
		# Get volumetric absorption coefficient from property cache
		# Note: prop_cache.α is already volumetric [m⁻¹] from effective_absorption = Σⱼ ξⱼ αⱼ
		# where ξⱼ [kg/m³] × αⱼ [m²/kg] = [m⁻¹]
		α_vol = prop_cache.α[i]

		rad_cache.α[i] = α_vol

		# Rosseland diffusion coefficient: D = 1/(3α)
		# Avoid division by zero for transparent regions
		if α_vol > 1e-10
			rad_cache.D_rad[i] = 1.0 / (3.0 * α_vol)
		else
			rad_cache.D_rad[i] = 1.0 / 3e-10  # Large but finite for transparent
		end
	end

	# Compute face-centered D_rad using harmonic mean
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

			# Harmonic mean weighted by distance
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

"""
	compute_radiation_flux_divergence!(rad_cache::RadiationCache, mesh::Mesh1D, prop_cache::PropertyCache)

Compute the divergence of radiative heat flux from P1 model.

∇·q_rad = α(G - 4σT⁴)

where:
- α is volumetric absorption coefficient [1/m]
- G is incident radiation [W/m²]
- σ is Stefan-Boltzmann constant
- T is temperature [K]

Note: For gray medium with emission, 4σT⁴ term represents blackbody emission.
"""
function compute_radiation_flux_divergence!(rad_cache::RadiationCache, mesh::Mesh1D, prop_cache::PropertyCache)
	n_cells = mesh.n_active_cells
	σ = 5.670374419e-8  # Stefan-Boltzmann constant

	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T
		G = rad_cache.G[i]
		α = rad_cache.α[i]

		# ∇·q_rad = α(G - 4σT⁴)
		# Positive when absorption > emission (net heating)
		rad_cache.div_q_rad[i] = α * (G - 4.0 * σ * T^4)
	end

	return nothing
end