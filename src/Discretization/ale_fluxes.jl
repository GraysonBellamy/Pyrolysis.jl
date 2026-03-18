"""
ALE advection computation using upwind differencing.

Reference: Design Document Section 3.2-3.3

In the ALE (Arbitrary Lagrangian-Eulerian) formulation, mesh motion is handled via
the ALE time derivative identity (chain rule for a moving sample point):

    ∂f/∂t|_χ = ∂f/∂t|_x + w·∇f

where ∂/∂t|_x is the Eulerian (fixed-point) rate and ∂/∂t|_χ is the ALE (mesh-point) rate.

The ALE equations are (Eulerian RHS + mesh correction):
- Energy: ∂T/∂t|_χ = (Eulerian RHS) + w·∇T      (T is intensive, no dilation)
- Gas species: ∂ξ/∂t|_χ = (Eulerian RHS) + w·∇ξ  (gas is spatially stationary)
- Mesh-bound species (solid/liquid): ∂ξ/∂t|_χ = S - ξ·θ  (dilation, handled separately)

The w·∇f advection terms are discretized using TVD-limited upwind differencing.
A 3-point stencil {i-2, i-1, i} (for w>0) is used with a slope limiter to achieve
second-order accuracy on smooth fields while preventing oscillations near
discontinuities. The limiter ψ(r) is applied to the ratio of consecutive gradients.
Falls back to first-order at boundaries where the upwind-upwind cell is unavailable.
Default limiter: VANLEER.
"""

using ..Core: GAS, SOLID
# Limiter, VANLEER, NONE, apply_limiter are available from flux_limiters.jl
# (included before this file in Discretization module)

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
as the harmonic mean at faces since the property cache now stores λ_face.
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
# ALE Advection (Upwind Formulation)
# ============================================================================
#
# ALE time derivative identity (chain rule):
#     ∂f/∂t|_χ = ∂f/∂t|_x + w·∇f
#
# For the non-conservative temperature equation:
#     ρcₚ ∂T/∂t|_x = -∇·q + Q     (Eulerian)
#     ∂T/∂t|_χ = (-∇·q + Q)/ρcₚ + w·∇T     (ALE = Eulerian + w·∇T)
#
# The w·∇T term is ADDED to the Eulerian RHS (not subtracted).
# It is discretized using first-order upwind differencing on the sign of w.
# ============================================================================

"""
	compute_ale_advection_term!(advection::Vector{Float64},
								mesh::Mesh1D{NC},
								w_node::Vector{Float64}) where NC

Compute the ALE advection term w·∇T using upwind differencing.

From the ALE identity: ∂T/∂t|_χ = ∂T/∂t|_x + w·∇T
The result w·∇T is ADDED to the Eulerian RHS to form the ALE rate.

# Arguments
- `advection`: Output vector of (w·∇T)ᵢ [K/s] for each cell
- `mesh`: Mesh with current cell states (temperatures)
- `w_node`: Mesh velocities at nodes [m/s]

# Sign Convention
- w < 0: mesh contracting (nodes moving toward z=0)
- w > 0: mesh expanding (nodes moving away from z=0)
- The result is ADDED to the Eulerian dT/dt to get the ALE dT/dt

# Boundary Treatment
- Uses zero-gradient extrapolation (Neumann BC) at all boundaries

# Algorithm
For cell i with mesh velocity wᵢ (average of bounding nodes):
- If wᵢ > 0: use backward difference (T_i - T_left) / Δz  [upwind from left]
- If wᵢ < 0: use forward difference (T_right - T_i) / Δz  [upwind from right]

Upwinding is based on the sign of w (the coefficient of ∇T in the added term).
"""
function compute_ale_advection_term!(advection::Vector{Float64},
	mesh::Mesh1D{NC},
	w_node::Vector{Float64};
	limiter::Limiter = VANLEER) where NC

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

		# Get neighbor cell IDs via faces
		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]

		left_id = face_left.cell_ids[1]   # Cell to left (lower z)
		right_id = face_right.cell_ids[2] # Cell to right (higher z)

		# Get neighbor temperatures with appropriate boundary treatment
		T_left = _get_boundary_temperature(left_id, T_i, mesh)
		T_right = _get_boundary_temperature(right_id, T_i, mesh)

		# TVD-limited upwind gradient (3-point stencil with slope limiter)
		# Falls back to first-order at boundaries where upwind-upwind unavailable
		if w_cell > 0
			# Upwind from left: stencil {left_left, left, i}
			left_left_id = _upwind_upwind_id(mesh, left_id, true)
			T_left_left = _get_boundary_temperature(left_left_id, T_left, mesh)
			d_i_L = _upwind_gradient_distance(mesh, i, left_id)
			d_L_LL = left_left_id != 0 ?
				abs(mesh.cell_centers[left_id] - mesh.cell_centers[left_left_id]) : d_i_L
			dT_dz = _limited_gradient(T_i, T_left, T_left_left,
				d_i_L, d_L_LL, left_left_id != 0, limiter, true)
		else
			# Upwind from right: stencil {right_right, right, i}
			right_right_id = _upwind_upwind_id(mesh, right_id, false)
			T_right_right = _get_boundary_temperature(right_right_id, T_right, mesh)
			d_i_R = _upwind_gradient_distance(mesh, i, right_id)
			d_R_RR = right_right_id != 0 ?
				abs(mesh.cell_centers[right_id] - mesh.cell_centers[right_right_id]) : d_i_R
			dT_dz = _limited_gradient(T_i, T_right, T_right_right,
				d_i_R, d_R_RR, right_right_id != 0, limiter, false)
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
Compute the correct upwind gradient distance: center-to-center distance to
the upwind neighbor cell.

On non-uniform meshes, this differs from the cell width. Using cell width
introduces a zeroth-order multiplicative bias in the gradient (the gradient
converges to (d/Δz)×f' instead of f' as h→0). Center-to-center distance
gives a consistent first-order upwind approximation.

At boundaries (neighbor_id == 0), the zero-gradient ghost cell is at the
face, so the distance is from cell center to face ≈ half cell width.
"""
@inline function _upwind_gradient_distance(mesh::Mesh1D, i::Int,
	neighbor_id::Int)
	if neighbor_id == 0
		# Boundary: distance from center to face (half cell width)
		return mesh.cell_volumes[i] / mesh.cross_section_area
	end
	return abs(mesh.cell_centers[i] - mesh.cell_centers[neighbor_id])
end

"""
Get the upwind-upwind neighbor cell ID for second-order TVD reconstruction.
Returns 0 if the upwind-upwind cell doesn't exist (boundary or upwind is boundary).
"""
@inline function _upwind_upwind_id(mesh::Mesh1D, upwind_id::Int, upwind_from_left::Bool)
	if upwind_id == 0
		return 0
	end
	cell = mesh.cells[upwind_id]
	if upwind_from_left
		f = cell.face_ids[1]
		return mesh.faces[f].cell_ids[1]
	else
		f = cell.face_ids[2]
		return mesh.faces[f].cell_ids[2]
	end
end

"""
Compute TVD-limited spatial gradient ∂T/∂z at a cell center using upwind
reconstruction. The gradient is always in the +z direction regardless of
which side is upwind.

For w > 0 (upwind from left): grad = (T_i - T_left) / d, positive when T increases with z
For w < 0 (upwind from right): grad = (T_right - T_i) / d, positive when T increases with z

The TVD correction uses a 3-point upwind stencil:
	grad = grad_1 + ψ(r)/2 × (grad_1 - grad_0)

where grad_1 is the first-order upwind gradient (in +z direction),
grad_0 is the same gradient one cell further upwind,
and r = grad_0/grad_1 is the slope ratio.

Arguments use directional convention:
- T_i, T_upwind, T_upup: values at current, upwind, and upwind-upwind cells
- upwind_from_left: if true, gradient = (T_i - T_upwind)/d; if false, = (T_upwind - T_i)/d

Falls back to first-order if upwind-upwind cell doesn't exist.
"""
@inline function _limited_gradient(T_i::Float64, T_upwind::Float64, T_upup::Float64,
	d_to_upwind::Float64, d_upwind_to_upup::Float64,
	has_upup::Bool, limiter::Limiter, upwind_from_left::Bool)

	# First-order gradient in +z direction
	if upwind_from_left
		# Upwind is to the left: grad = (T_i - T_left) / d
		grad_1 = (T_i - T_upwind) / d_to_upwind
	else
		# Upwind is to the right: grad = (T_right - T_i) / d
		grad_1 = (T_upwind - T_i) / d_to_upwind
	end

	if !has_upup || limiter == NONE
		return grad_1
	end

	# Previous gradient (one cell further upwind, same +z direction)
	if upwind_from_left
		grad_0 = (T_upwind - T_upup) / d_upwind_to_upup
	else
		grad_0 = (T_upup - T_upwind) / d_upwind_to_upup
	end

	if abs(grad_1) < 1e-20
		return grad_1
	end
	r = grad_0 / grad_1

	ψ = apply_limiter(r, limiter)

	return grad_1 + 0.5 * ψ * (grad_1 - grad_0)
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
	w_node::Vector{Float64};
	limiter::Limiter = VANLEER) where NC

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

		# Skip if velocity is essentially zero
		if abs(w_cell) < 1e-20
			for j in 1:NC
				advection[j, i] = 0.0
			end
			continue
		end

		# Compute upwind distances and upwind-upwind neighbor for TVD
		if w_cell > 0
			d_upwind = _upwind_gradient_distance(mesh, i, left_id)
			upup_id = _upwind_upwind_id(mesh, left_id, true)
			d_upup = upup_id != 0 ?
				abs(mesh.cell_centers[left_id] - mesh.cell_centers[upup_id]) : d_upwind
		else
			d_upwind = _upwind_gradient_distance(mesh, i, right_id)
			upup_id = _upwind_upwind_id(mesh, right_id, false)
			d_upup = upup_id != 0 ?
				abs(mesh.cell_centers[right_id] - mesh.cell_centers[upup_id]) : d_upwind
		end
		has_upup = upup_id != 0

		# Get cell state and neighbor states for upwind
		ξ_i = mesh.cell_states[i].ξ
		ξ_upwind = w_cell > 0 ?
			(left_id == 0 ? ξ_i : mesh.cell_states[left_id].ξ) :
			(right_id == 0 ? ξ_i : mesh.cell_states[right_id].ξ)
		ξ_upup = has_upup ? mesh.cell_states[upup_id].ξ : ξ_upwind

		# Use unrolled version for type stability
		upwind_from_left = w_cell > 0
		_compute_species_advection_limited_unrolled!(advection, i, w_cell,
			d_upwind, d_upup, has_upup, upwind_from_left,
			ξ_i, ξ_upwind, ξ_upup, material, limiter)
	end

	return nothing
end

# TVD-limited species advection (compile-time unrolled for type stability)
@generated function _compute_species_advection_limited_unrolled!(advection::Matrix{Float64},
	i::Int, w_cell::Float64, d_upwind::Float64, d_upup::Float64, has_upup::Bool,
	upwind_from_left::Bool,
	ξ_i::NTuple{NC}, ξ_upwind::NTuple{NC}, ξ_upup::NTuple{NC},
	material::Material{NC}, limiter::Limiter) where NC
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_i_j = ξ_i[$j]
				ξ_up_j = ξ_upwind[$j]
				ξ_upup_j = ξ_upup[$j]

				# TVD-limited gradient in +z direction
				dξ_dz = _limited_gradient(ξ_i_j, ξ_up_j, ξ_upup_j,
					d_upwind, d_upup, has_upup, limiter, upwind_from_left)
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
#     Jⱼ = ξⱼ × (v_gas - w) - (λ_face/T_face) × ∂(ξⱼT)/∂x
#          ↑ advection (ALE)    ↑ diffusion (volume fraction gradient)
#
# This properly accounts for:
# - Gas moving through the porous medium (Darcy flow, v_gas)
# - Mesh motion (ALE, w)
# - Diffusive transport (volume fraction gradient)
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
							λ_face::Float64, T_L::Float64, T_R::Float64,
							d_LR::Float64) -> Float64

Compute gas flux at a face in ALE frame using Darcy-Fick formulation [kg/(m²·s)].

Combines advective transport with relative velocity and diffusive transport:
	J = ξ_upwind × (v_g - w) - (λ_face/T_face) × (ξ_R×T_R - ξ_L×T_L) / d_LR
		↑ advection (ALE)       ↑ diffusion (volume fraction gradient)

# Arguments
- `v_g`: Darcy gas velocity at face [m/s]
- `w`: Mesh velocity at face [m/s]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `λ_face`: Face transport coefficient (harmonic mean of λ) [m²/s]
- `T_L, T_R`: Temperatures on left and right [K]
- `d_LR`: Distance between cell centers [m]

# Returns
Total flux [kg/(m²·s)], positive = flow in +z direction
"""
@inline function ale_gas_flux_darcy_fick(v_g::Float64, w::Float64,
	ξ_L::Float64, ξ_R::Float64,
	λ_face::Float64, T_L::Float64, T_R::Float64,
	d_LR::Float64)
	# Relative velocity in ALE frame
	v_rel = v_g - w

	# Diffusive flux (volume fraction / partial pressure gradient)
	T_face = 0.5 * (T_L + T_R)
	J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)

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
- Diffusive transport (volume fraction gradient)

The flux formula is:
	Jⱼ = ξⱼ × (v_gas - w) - (λ_face/T_face) × ∂(ξⱼT)/∂x

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache (must have λ_face computed)
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
		λ_face = cache.λ_face[f]
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		v_g = v_gas[f]

		# Get concentration tuples for this face
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Use compile-time unrolled version for type stability
		_compute_ale_gas_fluxes_darcy_fick_unrolled!(J, f, v_g, w_face,
			ξ_L, ξ_R, λ_face, T_L, T_R, d_LR, material)
	end

	return nothing
end

"""
Compile-time unrolled ALE Darcy-Fick flux computation for type stability.
"""
@generated function _compute_ale_gas_fluxes_darcy_fick_unrolled!(J::Matrix, f::Int,
	v_g, w_face,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC},
	λ_face, T_L, T_R, d_LR,
	material::Material{NC}) where NC

	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]

				# Relative velocity in ALE frame
				v_rel = v_g - w_face

				# Diffusive flux (volume fraction gradient)
				T_face = 0.5 * (T_L + T_R)
				J_diff = -λ_face * (ξ_R_j * T_R - ξ_L_j * T_L) / (T_face * d_LR)

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
	w_node::Vector{Float64};
	limiter::Limiter = VANLEER) where NC

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

		# Skip if relative velocity is essentially zero
		if abs(v_rel) < 1e-20
			for j in 1:NC
				advection[j, i] = 0.0
			end
			continue
		end

		# Compute upwind distances and upwind-upwind neighbor for TVD
		if v_rel > 0
			d_upwind = _upwind_gradient_distance(mesh, i, left_id)
			upup_id = _upwind_upwind_id(mesh, left_id, true)
			d_upup = upup_id != 0 ?
				abs(mesh.cell_centers[left_id] - mesh.cell_centers[upup_id]) : d_upwind
		else
			d_upwind = _upwind_gradient_distance(mesh, i, right_id)
			upup_id = _upwind_upwind_id(mesh, right_id, false)
			d_upup = upup_id != 0 ?
				abs(mesh.cell_centers[right_id] - mesh.cell_centers[upup_id]) : d_upwind
		end
		has_upup = upup_id != 0

		# Get cell state and neighbor states for upwind
		ξ_i = mesh.cell_states[i].ξ
		ξ_upwind = v_rel > 0 ?
			(left_id == 0 ? ξ_i : mesh.cell_states[left_id].ξ) :
			(right_id == 0 ? ξ_i : mesh.cell_states[right_id].ξ)
		ξ_upup = has_upup ? mesh.cell_states[upup_id].ξ : ξ_upwind

		# Reuse same TVD-limited unrolled helper as non-Darcy version
		upwind_from_left = v_rel > 0
		_compute_species_advection_limited_unrolled!(advection, i, v_rel,
			d_upwind, d_upup, has_upup, upwind_from_left,
			ξ_i, ξ_upwind, ξ_upup, material, limiter)
	end

	return nothing
end
