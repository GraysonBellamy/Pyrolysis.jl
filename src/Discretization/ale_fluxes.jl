"""
ALE advection computation using upwind differencing.

Reference: Design Document Section 3.3, ALE_Energy_Reformulation.md

In the ALE (Arbitrary Lagrangian-Eulerian) formulation, the correct treatment
of mesh motion uses advection terms (w·∇φ) rather than conservative flux
divergence (∇·(φw)). This avoids spurious energy sources/sinks at boundaries.

The ALE equations are:
- Energy: ∂T/∂t = (physics terms) - w·∇T
- Species: ∂ξ/∂t = (physics terms) - w·∇ξ

The advection terms are discretized using first-order upwind differencing
based on the direction of mesh motion.
"""

using ..Core: GAS, SOLID

# ============================================================================
# ALE Convective Flux Computation
# ============================================================================

"""
	ale_convective_flux(φ_L::Float64, φ_R::Float64, v::Float64, w::Float64) -> Float64

Compute ALE convective flux at a face using first-order upwinding.

F = (v - w) * φ_upwind

# Arguments
- `φ_L, φ_R`: Field values in left and right cells
- `v`: Material velocity at face [m/s]
- `w`: Mesh velocity at face [m/s]

# Returns
Convective flux [φ units * m/s], positive = flow in +z direction
"""
@inline function ale_convective_flux(φ_L::Float64, φ_R::Float64,
	v::Float64, w::Float64)
	v_rel = v - w  # Relative velocity
	φ_face = v_rel >= 0.0 ? φ_L : φ_R
	return v_rel * φ_face
end

"""
	ale_convective_flux_limited(φ_LL::Float64, φ_L::Float64, φ_R::Float64, φ_RR::Float64,
								v::Float64, w::Float64, limiter::Limiter) -> Float64

Compute ALE convective flux with second-order TVD limiting.
"""
@inline function ale_convective_flux_limited(φ_LL::Float64, φ_L::Float64,
	φ_R::Float64, φ_RR::Float64,
	v::Float64, w::Float64,
	limiter::Limiter)
	v_rel = v - w
	φ_face = convective_face_value(φ_LL, φ_L, φ_R, φ_RR, v_rel, limiter)
	return v_rel * φ_face
end

# ============================================================================
# Gas Velocity Computation (Darcy flow)
# ============================================================================

"""
	compute_gas_velocity!(v_gas::Vector{Float64}, mesh::Mesh1D{NC}, cache::PropertyCache,
						  material, Δp::Vector{Float64}) where NC

Compute gas velocity at faces using Darcy's law.

v_g = -κ/μ * ∇p

For simplified pyrolysis models, gas velocity can be estimated from
mass generation and porosity.

Note: This uses λ (gas transfer coefficient) as a proxy for κ/μ, computed
as the harmonic mean at faces since the property cache now stores λρ_face.
"""
function compute_gas_velocity!(v_gas::Vector{Float64}, mesh::Mesh1D{NC}, cache,
	material, pressure_gradient::Vector{Float64}) where NC
	n_faces_val = length(mesh.faces)

	# Simplified: use gas transfer coefficient as effective permeability/viscosity
	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			v_gas[f] = 0.0
		else
			# Compute harmonic mean of λ at face (not stored in cache anymore)
			λ_L = cache.λ[left_id]
			λ_R = cache.λ[right_id]
			λ_sum = λ_L + λ_R
			λ_face = λ_sum > 1e-20 ? 2.0 * λ_L * λ_R / λ_sum : 0.0
			v_gas[f] = -λ_face * pressure_gradient[f]
		end
	end

	return nothing
end

"""
	estimate_gas_velocity_from_generation!(v_gas::Vector{Float64}, mesh::Mesh1D{NC}, 
										   mass_generation_rate::Vector{Float64},
										   ρ_gas::Float64) where NC

Estimate gas velocity from mass generation (simplified model).

Integrate mass generation from bottom to surface to get velocity.
"""
function estimate_gas_velocity_from_generation!(v_gas::Vector{Float64}, mesh::Mesh1D{NC},
	mass_generation_rate::Vector{Float64},
	ρ_gas::Float64) where NC
	n_faces_val = length(mesh.faces)
	n_cells_val = length(mesh.cells)

	# Integrate from bottom (last cell) toward surface
	# Mass flux at each face = integral of generation rate below
	cumulative_mass_flux = 0.0

	# Assume cells are ordered from surface (1) to bottom (n)
	# Gas flows toward surface (negative z direction in our convention)
	@inbounds for f in n_faces_val:-1:1
		face = mesh.faces[f]
		right_id = face.cell_ids[2]

		if right_id != 0
			# Add generation from this cell
			V = mesh.cell_volumes[right_id]
			A = mesh.face_areas[f]
			cumulative_mass_flux += mass_generation_rate[right_id] * V / A
		end

		# Velocity from mass flux
		if ρ_gas > 0.0
			v_gas[f] = -cumulative_mass_flux / ρ_gas  # Negative = toward surface
		else
			v_gas[f] = 0.0
		end
	end

	return nothing
end

# ============================================================================
# ALE Advection (Upwind Formulation) - CORRECTED FORMULATION
# ============================================================================
# 
# The correct ALE energy equation is:
#     ρcₚ (∂T/∂t + w·∇T) = -∇·q + Q
#
# NOT the conservative flux form:
#     ∂(ρcₚT)/∂t + ∇·(ρcₚTw) = ...  [WRONG]
#
# The advection term w·∇T is discretized using first-order upwind differencing.
# This avoids the spurious energy sources/sinks that occur with the flux-based
# formulation at fixed boundaries.
#
# Reference: ALE_Energy_Reformulation.md
# ============================================================================

"""
	compute_ale_advection_term!(advection::Vector{Float64},
								mesh::Mesh1D{NC},
								w_node::Vector{Float64}) where NC

Compute the ALE advection contribution to dT/dt using upwind differencing.

The ALE energy equation is:
	ρcₚ (∂T/∂t + w·∇T) = -∇·q + Q

This function computes (w·∇T)ᵢ for each cell using first-order upwind.

# Arguments
- `advection`: Output vector of advection terms [K/s] for each cell
- `mesh`: Mesh with current cell states (temperatures)
- `w_node`: Mesh velocities at nodes [m/s]

# Sign Convention
- w < 0: mesh contracting (nodes moving toward z=0)
- w > 0: mesh expanding (nodes moving away from z=0)
- The result is the advection term that should be SUBTRACTED from dT/dt

# Boundary Treatment
- Uses zero-gradient extrapolation (Neumann BC) at all boundaries
- This is physically appropriate for the advection/convection operator

# Algorithm
For cell i with mesh velocity wᵢ (average of bounding nodes):
- If wᵢ > 0: use backward difference (T_i - T_left) / Δz  [information from left]
- If wᵢ < 0: use forward difference (T_right - T_i) / Δz  [information from right]

This ensures numerical stability by using upwind differencing based on the
direction of mesh motion.
"""
function compute_ale_advection_term!(advection::Vector{Float64},
	mesh::Mesh1D{NC},
	w_node::Vector{Float64}) where NC

	n_cells = length(mesh.cells)

	@inbounds for i in 1:n_cells
		cell = mesh.cells[i]
		n_left, n_right = cell.node_ids

		# Cell-centered velocity (average of bounding nodes)
		w_cell = 0.5 * (w_node[n_left] + w_node[n_right])

		# Skip if velocity is essentially zero (no advection)
		if abs(w_cell) < 1e-20
			advection[i] = 0.0
			continue
		end

		T_i = mesh.cell_states[i].T
		Δz = mesh.cell_volumes[i] / mesh.cross_section_area

		# Get neighbor cell IDs via faces
		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]

		left_id = face_left.cell_ids[1]   # Cell to left (lower z)
		right_id = face_right.cell_ids[2] # Cell to right (higher z)

		# Get neighbor temperatures with appropriate boundary treatment
		T_left = _get_boundary_temperature(left_id, T_i, mesh)
		T_right = _get_boundary_temperature(right_id, T_i, mesh)

		# Upwind discretization: use gradient in direction mesh is moving
		if w_cell > 0
			# Mesh expanding toward +z: use backward difference
			# Information flows from left to right
			dT_dz = (T_i - T_left) / Δz
		else
			# Mesh contracting toward -z: use forward difference  
			# Information flows from right to left
			dT_dz = (T_right - T_i) / Δz
		end

		# Advection term: w · ∇T
		advection[i] = w_cell * dT_dz
	end

	return nothing
end

"""
Helper function to get temperature at boundary for ALE advection.
Uses zero-gradient extrapolation (cell temperature) for boundary cells.

The zero-gradient approach is physically appropriate for the advection operator.
"""
@inline function _get_boundary_temperature(neighbor_id::Int, T_cell::Float64,
	mesh::Mesh1D)

	# Interior cell - use neighbor temperature
	if neighbor_id != 0
		return mesh.cell_states[neighbor_id].T
	end

	# Boundary face - use zero-gradient (cell temperature)
	# This is physically appropriate for the convective/advective operator
	return T_cell
end

"""
	compute_ale_species_advection!(advection::Matrix{Float64},
								   mesh::Mesh1D{NC},
								   material::Material{NC},
								   w_node::Vector{Float64}) where NC

Compute the ALE advection contribution to dξ/dt using upwind differencing.

The ALE species equation is:
	∂ξⱼ/∂t + (vⱼ - w)·∇ξⱼ = -∇·Jⱼ + Sⱼ

Where (vⱼ - w) is the relative velocity of component j to the mesh.

# Arguments
- `advection`: Output matrix [NC × n_cells] of advection terms [kg/(m³·s)]
- `mesh`: Mesh with current cell states
- `material`: Material definition (to determine solid vs gas components)
- `w_node`: Mesh velocities at nodes [m/s]

# Physics
- For SOLID components: The material moves with the mesh (Lagrangian tracking),
  so v_solid = w and (v_solid - w) = 0. Therefore, NO advection term for solids.
- For LIQUID components: Treated same as solids (bound to solid matrix, moves with mesh).
- For GAS components: In the simplified model where bulk gas velocity v_gas ≈ 0,
  the relative velocity is (0 - w) = -w. The advection term is -w·∇ξ_gas.

This ensures mass-thickness consistency: solid mass is locked to the mesh,
so when the mesh shrinks due to mass loss, thickness tracks mass correctly.
"""
function compute_ale_species_advection!(advection::Matrix{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC},
	w_node::Vector{Float64}) where NC

	n_cells = length(mesh.cells)

	@inbounds for i in 1:n_cells
		cell = mesh.cells[i]
		n_left, n_right = cell.node_ids

		# Cell-centered velocity
		w_cell = 0.5 * (w_node[n_left] + w_node[n_right])

		# Get neighbor cell IDs
		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]

		left_id = face_left.cell_ids[1]
		right_id = face_right.cell_ids[2]

		Δz = mesh.cell_volumes[i] / mesh.cross_section_area

		# Skip if velocity is essentially zero
		if abs(w_cell) < 1e-20
			for j in 1:NC
				advection[j, i] = 0.0
			end
			continue
		end

		# Get cell state and neighbor states for upwind
		ξ_i = mesh.cell_states[i].ξ
		ξ_left_all = left_id == 0 ? ξ_i : mesh.cell_states[left_id].ξ
		ξ_right_all = right_id == 0 ? ξ_i : mesh.cell_states[right_id].ξ

		# Use unrolled version for type stability
		_compute_species_advection_unrolled!(advection, i, w_cell, Δz,
			ξ_i, ξ_left_all, ξ_right_all, material)
	end

	return nothing
end

# Dispatch-based unrolled species advection for type stability
@generated function _compute_species_advection_unrolled!(advection::Matrix{Float64},
	i::Int, w_cell::Float64, Δz::Float64,
	ξ_i::NTuple{NC}, ξ_left::NTuple{NC}, ξ_right::NTuple{NC},
	material::Material{NC}) where NC
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				# GAS components can move relative to mesh
				# With v_gas ≈ 0, relative velocity is -w
				ξ_i_j = ξ_i[$j]
				ξ_left_j = ξ_left[$j]
				ξ_right_j = ξ_right[$j]

				# Upwind discretization based on mesh velocity direction
				if w_cell > 0
					dξ_dz = (ξ_i_j - ξ_left_j) / Δz
				else
					dξ_dz = (ξ_right_j - ξ_i_j) / Δz
				end
				advection[$j, i] = w_cell * dξ_dz
			else
				# SOLID and LIQUID components move with the mesh (v = w)
				# So relative velocity (v - w) = 0 → no advection
				advection[$j, i] = 0.0
			end
		end)
	end
	return Expr(:block, exprs...)
end

# ============================================================================
# ALE Darcy-Fick Gas Fluxes
# ============================================================================
#
# In the ALE frame, gas transport uses the RELATIVE velocity:
#     v_rel = v_gas - w_mesh
#
# The combined Darcy-Fick flux in ALE frame is:
#     Jⱼ = ξⱼ × (v_gas - w) - λρ_face × ∂Yⱼ/∂x
#          ↑ advection (ALE)    ↑ diffusion
#
# This properly accounts for:
# - Gas moving through the porous medium (Darcy flow, v_gas)
# - Mesh motion (ALE, w)
# - Diffusive transport (Fick's law)
# ============================================================================

"""
	mesh_velocity_at_face(w_node::Vector{Float64}, mesh::Mesh1D{NC}, f::Int) -> Float64

Get mesh velocity at face f.

In 1D, the face is located at a single node, so we use the node velocity directly.
"""
@inline function mesh_velocity_at_face(w_node::Vector{Float64}, mesh::Mesh1D{NC}, f::Int) where NC
	node_id = mesh.faces[f].node_id
	return w_node[node_id]
end

"""
	ale_gas_flux_darcy_fick(v_g::Float64, w::Float64, ξ_L::Float64, ξ_R::Float64,
							λρ_face::Float64, ρ_L::Float64, ρ_R::Float64,
							d_LR::Float64) -> Float64

Compute gas flux at a face in ALE frame using Darcy-Fick formulation [kg/(m²·s)].

Combines advective transport with relative velocity and diffusive transport:
	J = ξ_upwind × (v_g - w) - λρ_face × (Y_R - Y_L) / d_LR
		↑ advection (ALE)       ↑ diffusion

# Arguments
- `v_g`: Darcy gas velocity at face [m/s]
- `w`: Mesh velocity at face [m/s]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `λρ_face`: Face transport coefficient [kg/(m·s)]
- `ρ_L, ρ_R`: Total mixture densities [kg/m³]
- `d_LR`: Distance between cell centers [m]

# Returns
Total flux [kg/(m²·s)], positive = flow in +z direction
"""
@inline function ale_gas_flux_darcy_fick(v_g::Float64, w::Float64,
	ξ_L::Float64, ξ_R::Float64,
	λρ_face::Float64, ρ_L::Float64, ρ_R::Float64,
	d_LR::Float64)
	# Relative velocity in ALE frame
	v_rel = v_g - w

	# Diffusive flux (same as standard Fick)
	Y_L = ξ_L / ρ_L
	Y_R = ξ_R / ρ_R
	J_diff = -λρ_face * (Y_R - Y_L) / d_LR

	# Advective flux with first-order upwinding based on relative velocity
	ξ_upwind = v_rel >= 0.0 ? ξ_L : ξ_R
	J_adv = v_rel * ξ_upwind

	return J_adv + J_diff
end

"""
	compute_ale_gas_fluxes_darcy_fick!(J::Matrix{Float64}, mesh::Mesh1D{NC},
									   cache::PropertyCache, material::Material{NC},
									   v_gas::Vector{Float64},
									   w_node::Vector{Float64}) where NC

Compute gas mass fluxes in ALE frame using full Darcy-Fick formulation.

This combines:
- Darcy advective transport (gas flow through porous medium)
- ALE mesh motion correction
- Fickian diffusive transport

The flux formula is:
	Jⱼ = ξⱼ × (v_gas - w) - λρ_face × ∂Yⱼ/∂x

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache (must have λρ_face, ρ computed)
- `material`: Material definition
- `v_gas`: Darcy gas velocity at faces [m/s]
- `w_node`: Mesh velocity at nodes [m/s]

# Notes
- When w_node = 0, reduces to standard Darcy-Fick
- When v_gas = 0, gives ALE advection due to mesh motion only
- Solid components have zero flux (transport handled by mesh motion)
"""
function compute_ale_gas_fluxes_darcy_fick!(J::Matrix{Float64}, mesh::Mesh1D{NC},
	cache::PropertyCache, material::Material{NC},
	v_gas::Vector{Float64},
	w_node::Vector{Float64}) where NC
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

		# Get mesh velocity at face
		w_face = mesh_velocity_at_face(w_node, mesh, f)

		# Interior face properties
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		λρ_face = cache.λρ_face[f]
		ρ_L = cache.ρ[left_id]
		ρ_R = cache.ρ[right_id]
		v_g = v_gas[f]

		# Get concentration tuples for this face
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Use compile-time unrolled version for type stability
		_compute_ale_gas_fluxes_darcy_fick_unrolled!(J, f, v_g, w_face,
			ξ_L, ξ_R, λρ_face, ρ_L, ρ_R, d_LR, material)
	end

	return nothing
end

"""
Compile-time unrolled ALE Darcy-Fick flux computation for type stability.
"""
@generated function _compute_ale_gas_fluxes_darcy_fick_unrolled!(J::Matrix, f::Int,
	v_g, w_face,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC},
	λρ_face, ρ_L, ρ_R, d_LR,
	material::Material{NC}) where NC

	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]

				# Relative velocity in ALE frame
				v_rel = v_g - w_face

				# Diffusive flux: J_diff = -λρ_face × ΔY / d
				J_diff = -λρ_face * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR

				# Advective flux with relative velocity (first-order upwind)
				ξ_upwind = v_rel >= 0.0 ? ξ_L_j : ξ_R_j
				J_adv = v_rel * ξ_upwind

				# Combined ALE Darcy-Fick flux
				J[$j, f] = J_adv + J_diff
			end
		end)
	end
	return Expr(:block, exprs...)
end

"""
	compute_ale_species_advection_darcy!(advection::Matrix{Float64},
										 mesh::Mesh1D{NC},
										 material::Material{NC},
										 v_gas::Vector{Float64},
										 w_node::Vector{Float64}) where NC

Compute species advection in ALE frame accounting for Darcy gas velocity.

This is an alternative formulation using advection terms (w·∇ξ) rather than
conservative fluxes. Use when you want to treat advection separately from
flux-based transport.

For gas components: advection = (v_gas - w) · ∇ξ
For solid components: advection = 0 (move with mesh)

# Arguments
- `advection`: Output matrix [NC × n_cells] of advection terms [kg/(m³·s)]
- `mesh`: Mesh with current cell states
- `material`: Material definition
- `v_gas`: Darcy gas velocity at faces [m/s]
- `w_node`: Mesh velocities at nodes [m/s]
"""
function compute_ale_species_advection_darcy!(advection::Matrix{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC},
	v_gas::Vector{Float64},
	w_node::Vector{Float64}) where NC

	n_cells = length(mesh.cells)

	@inbounds for i in 1:n_cells
		cell = mesh.cells[i]
		n_left, n_right = cell.node_ids

		# Cell-centered mesh velocity
		w_cell = 0.5 * (w_node[n_left] + w_node[n_right])

		# Get neighbor cell IDs
		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]

		left_id = face_left.cell_ids[1]
		right_id = face_right.cell_ids[2]

		# Average gas velocity at cell center (from face values)
		v_g_cell = 0.5 * (v_gas[f_left] + v_gas[f_right])

		# Relative velocity in ALE frame
		v_rel = v_g_cell - w_cell

		Δz = mesh.cell_volumes[i] / mesh.cross_section_area

		# Skip if relative velocity is essentially zero
		if abs(v_rel) < 1e-20
			for j in 1:NC
				advection[j, i] = 0.0
			end
			continue
		end

		# Get cell state and neighbor states for upwind
		ξ_i = mesh.cell_states[i].ξ
		ξ_left_all = left_id == 0 ? ξ_i : mesh.cell_states[left_id].ξ
		ξ_right_all = right_id == 0 ? ξ_i : mesh.cell_states[right_id].ξ

		# Use unrolled version for type stability
		_compute_ale_species_advection_darcy_unrolled!(advection, i, v_rel, Δz,
			ξ_i, ξ_left_all, ξ_right_all, material)
	end

	return nothing
end

@generated function _compute_ale_species_advection_darcy_unrolled!(advection::Matrix{Float64},
	i::Int, v_rel::Float64, Δz::Float64,
	ξ_i::NTuple{NC}, ξ_left::NTuple{NC}, ξ_right::NTuple{NC},
	material::Material{NC}) where NC
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_i_j = ξ_i[$j]
				ξ_left_j = ξ_left[$j]
				ξ_right_j = ξ_right[$j]

				# Upwind discretization based on relative velocity direction
				if v_rel > 0
					dξ_dz = (ξ_i_j - ξ_left_j) / Δz
				else
					dξ_dz = (ξ_right_j - ξ_i_j) / Δz
				end
				advection[$j, i] = v_rel * dξ_dz
			else
				# SOLID and LIQUID components move with the mesh
				advection[$j, i] = 0.0
			end
		end)
	end
	return Expr(:block, exprs...)
end
