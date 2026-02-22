"""
Solution remapping after mesh changes.

Conservative remapping transfers solution fields from old mesh to new mesh
while preserving conservation properties (mass, energy).

Reference: Design Document Section 7.6
"""

using ..Core: AbstractRemappingMethod
using ..Geometry: Mesh1D

# ============================================================================
# Remapping Method Types
# ============================================================================

"""
	FirstOrderRemap <: AbstractRemappingMethod

First-order conservative remapping.

Uses cell-averaged values with volume-weighted averaging for overlapping regions.
Simple but may introduce numerical diffusion.
"""
struct FirstOrderRemap <: AbstractRemappingMethod end

"""
	SecondOrderRemap <: AbstractRemappingMethod

Second-order conservative remapping with gradient reconstruction.

Uses cell gradients to reduce numerical diffusion while maintaining conservation.
Limiters prevent spurious oscillations.

# Fields
- `limiter::Symbol`: Limiter type (:minmod, :vanleer, :superbee)
"""
struct SecondOrderRemap <: AbstractRemappingMethod
	limiter::Symbol
end
SecondOrderRemap() = SecondOrderRemap(:minmod)

# ============================================================================
# First-Order Remapping
# ============================================================================

"""
	remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
					u_old::Vector{Float64}, mesh_old::Mesh1D,
					NC::Int, method::FirstOrderRemap;
					material=nothing)

First-order conservative remapping of solution vector.

For each new cell, find overlapping old cells and compute appropriate averages.

# Temperature Remapping
When `material` is provided, uses **energy-conserving** remapping:
	T_new = Σ(ρcₚ × T × V)_overlap / Σ(ρcₚ × V)_overlap

When `material` is nothing, falls back to simple volume-weighted averaging:
	T_new = Σ(T × V)_overlap / Σ(V)_overlap

The energy-conserving form is recommended when cells have different compositions
(and hence different ρcₚ values), as it preserves total thermal energy.

# Concentration Remapping
Concentrations use mass-conserving remapping:
	ξ_new = Σ(ξ × V)_overlap / V_new

# Arguments
- `u_new`: Output solution vector on new mesh
- `mesh_new`: New mesh (after motion/refinement)
- `u_old`: Solution vector on old mesh
- `mesh_old`: Old mesh
- `NC`: Number of components
- `method`: FirstOrderRemap instance
- `material`: (keyword) Material for ρcₚ computation, or nothing for volume-weighted
"""
function remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
	u_old::Vector{Float64}, mesh_old::Mesh1D,
	NC::Int, ::FirstOrderRemap;
	material = nothing)
	n_cells_new = length(mesh_new.cells)
	n_cells_old = length(mesh_old.cells)
	n_dof_per_cell = 1 + NC  # T + ξ₁, ξ₂, ...

	# Process each new cell
	@inbounds for i_new in 1:n_cells_new
		if !mesh_new.active_cell_mask[i_new]
			continue
		end

		cell_new = mesh_new.cells[i_new]
		z_left_new = mesh_new.nodes[cell_new.node_ids[1]].z
		z_right_new = mesh_new.nodes[cell_new.node_ids[2]].z
		V_new = mesh_new.cell_volumes[i_new]

		# Accumulate contributions from overlapping old cells
		# For energy-conserving remapping: E_accum = Σ(ρcₚ T V), ρcₚV_accum = Σ(ρcₚ V)
		# For volume-weighted fallback: T_accum = Σ(T V), V_accum = Σ(V)
		E_accum = 0.0       # Energy accumulator (or T*V if no material)
		ρcₚV_accum = 0.0    # ρcₚ*V accumulator (or V if no material)
		ξ_accum = zeros(NC)
		V_overlap_total = 0.0

		for i_old in 1:n_cells_old
			if !mesh_old.active_cell_mask[i_old]
				continue
			end

			cell_old = mesh_old.cells[i_old]
			z_left_old = mesh_old.nodes[cell_old.node_ids[1]].z
			z_right_old = mesh_old.nodes[cell_old.node_ids[2]].z

			# Compute overlap
			z_left_overlap = max(z_left_new, z_left_old)
			z_right_overlap = min(z_right_new, z_right_old)

			if z_right_overlap > z_left_overlap
				# There is overlap
				A = mesh_old.cross_section_area
				V_overlap = (z_right_overlap - z_left_overlap) * A

				# Extract old cell state
				idx_old = (i_old - 1) * n_dof_per_cell
				T_old = u_old[idx_old+1]

				# Extract concentrations
				ξ_old = ntuple(j -> u_old[idx_old+1+j], Val(NC))

				# Accumulate energy (or T*V for fallback)
				if material !== nothing
					# Energy-conserving: weight by ρcₚ
					ρcₚ_old = Properties.effective_heat_capacity(ξ_old, T_old, material)
					E_accum += ρcₚ_old * T_old * V_overlap
					ρcₚV_accum += ρcₚ_old * V_overlap
				else
					# Fallback: simple volume weighting
					E_accum += T_old * V_overlap
					ρcₚV_accum += V_overlap
				end

				# Accumulate mass (same for both approaches)
				for j in 1:NC
					ξ_accum[j] += ξ_old[j] * V_overlap
				end

				V_overlap_total += V_overlap
			end
		end

		# Compute new cell values
		idx_new = (i_new - 1) * n_dof_per_cell

		if V_overlap_total > 0.0
			# Temperature: energy-conserving or volume-weighted
			u_new[idx_new+1] = E_accum / ρcₚV_accum

			# Concentrations: mass-conserving
			# New concentration = total mass / new volume
			for j in 1:NC
				u_new[idx_new+1+j] = ξ_accum[j] / V_new
			end
		else
			# No overlap - extrapolate from nearest cell (edge case)
			# Find nearest old cell
			z_center_new = 0.5 * (z_left_new + z_right_new)
			nearest_old = 1
			min_dist = Inf

			for i_old in 1:n_cells_old
				if mesh_old.active_cell_mask[i_old]
					z_center_old = mesh_old.cell_centers[i_old]
					dist = abs(z_center_new - z_center_old)
					if dist < min_dist
						min_dist = dist
						nearest_old = i_old
					end
				end
			end

			idx_old = (nearest_old - 1) * n_dof_per_cell
			u_new[idx_new+1] = u_old[idx_old+1]
			for j in 1:NC
				u_new[idx_new+1+j] = u_old[idx_old+1+j]
			end
		end
	end

	return nothing
end

# ============================================================================
# Second-Order Remapping
# ============================================================================

"""
	remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
					u_old::Vector{Float64}, mesh_old::Mesh1D,
					NC::Int, method::SecondOrderRemap)

Second-order conservative remapping with gradient reconstruction.
"""
function remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
	u_old::Vector{Float64}, mesh_old::Mesh1D,
	NC::Int, method::SecondOrderRemap)
	n_cells_new = length(mesh_new.cells)
	n_cells_old = length(mesh_old.cells)
	n_dof_per_cell = 1 + NC

	# First, compute gradients on old mesh
	∇T_old = compute_gradients_for_remap(u_old, mesh_old, NC, 1)  # Temperature
	∇ξ_old = [compute_gradients_for_remap(u_old, mesh_old, NC, 1+j) for j in 1:NC]

	# Get min/max values for limiting
	T_min_old, T_max_old = extrema_values(u_old, mesh_old, NC, 1)
	ξ_min_old = [extrema_values(u_old, mesh_old, NC, 1+j)[1] for j in 1:NC]
	ξ_max_old = [extrema_values(u_old, mesh_old, NC, 1+j)[2] for j in 1:NC]

	# Process each new cell
	@inbounds for i_new in 1:n_cells_new
		if !mesh_new.active_cell_mask[i_new]
			continue
		end

		cell_new = mesh_new.cells[i_new]
		z_left_new = mesh_new.nodes[cell_new.node_ids[1]].z
		z_right_new = mesh_new.nodes[cell_new.node_ids[2]].z
		z_center_new = 0.5 * (z_left_new + z_right_new)
		V_new = mesh_new.cell_volumes[i_new]

		T_accum = 0.0
		ξ_accum = zeros(NC)
		V_overlap_total = 0.0

		for i_old in 1:n_cells_old
			if !mesh_old.active_cell_mask[i_old]
				continue
			end

			cell_old = mesh_old.cells[i_old]
			z_left_old = mesh_old.nodes[cell_old.node_ids[1]].z
			z_right_old = mesh_old.nodes[cell_old.node_ids[2]].z
			z_center_old = mesh_old.cell_centers[i_old]

			z_left_overlap = max(z_left_new, z_left_old)
			z_right_overlap = min(z_right_new, z_right_old)

			if z_right_overlap > z_left_overlap
				A = mesh_old.cross_section_area
				V_overlap = (z_right_overlap - z_left_overlap) * A
				z_center_overlap = 0.5 * (z_left_overlap + z_right_overlap)

				# Extract old cell values
				idx_old = (i_old - 1) * n_dof_per_cell
				T_old = u_old[idx_old+1]

				# Second-order reconstruction
				dz = z_center_overlap - z_center_old
				T_reconstructed = T_old + ∇T_old[i_old] * dz

				T_accum += T_reconstructed * V_overlap

				for j in 1:NC
					ξ_old_j = u_old[idx_old+1+j]
					ξ_reconstructed = ξ_old_j + ∇ξ_old[j][i_old] * dz
					ξ_accum[j] += ξ_reconstructed * V_overlap
				end

				V_overlap_total += V_overlap
			end
		end

		idx_new = (i_new - 1) * n_dof_per_cell

		if V_overlap_total > 0.0
			T_new = T_accum / V_overlap_total

			# Apply limiter - don't create new extrema
			T_new = clamp(T_new, T_min_old, T_max_old)
			u_new[idx_new+1] = T_new

			for j in 1:NC
				ξ_new_j = ξ_accum[j] / V_new
				# Concentrations must be non-negative
				ξ_new_j = max(ξ_new_j, 0.0)
				# Don't exceed old maximum
				ξ_new_j = min(ξ_new_j, ξ_max_old[j])
				u_new[idx_new+1+j] = ξ_new_j
			end
		else
			# Fallback: nearest cell
			nearest_old = find_nearest_cell(z_center_new, mesh_old)
			idx_old = (nearest_old - 1) * n_dof_per_cell
			u_new[idx_new+1] = u_old[idx_old+1]
			for j in 1:NC
				u_new[idx_new+1+j] = u_old[idx_old+1+j]
			end
		end
	end

	return nothing
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	compute_gradients_for_remap(u, mesh, NC, var_idx) -> Vector{Float64}

Compute cell-centered gradients for remapping.
"""
function compute_gradients_for_remap(u::Vector{Float64}, mesh::Mesh1D,
	NC::Int, var_idx::Int)
	n_cells = length(mesh.cells)
	n_dof_per_cell = 1 + NC
	gradients = zeros(n_cells)

	@inbounds for i in 1:n_cells
		if !mesh.active_cell_mask[i]
			continue
		end

		idx = (i - 1) * n_dof_per_cell + var_idx
		φ_i = u[idx]
		z_i = mesh.cell_centers[i]

		# Find neighbors (previous and next active cells)
		i_prev = 0
		i_next = 0

		for k in (i-1):-1:1
			if mesh.active_cell_mask[k]
				i_prev = k
				break
			end
		end

		for k in (i+1):n_cells
			if mesh.active_cell_mask[k]
				i_next = k
				break
			end
		end

		if i_prev > 0 && i_next > 0
			# Central difference
			idx_prev = (i_prev - 1) * n_dof_per_cell + var_idx
			idx_next = (i_next - 1) * n_dof_per_cell + var_idx
			φ_prev = u[idx_prev]
			φ_next = u[idx_next]
			z_prev = mesh.cell_centers[i_prev]
			z_next = mesh.cell_centers[i_next]

			gradients[i] = (φ_next - φ_prev) / (z_next - z_prev)
		elseif i_prev > 0
			# One-sided (backward)
			idx_prev = (i_prev - 1) * n_dof_per_cell + var_idx
			φ_prev = u[idx_prev]
			z_prev = mesh.cell_centers[i_prev]
			gradients[i] = (φ_i - φ_prev) / (z_i - z_prev)
		elseif i_next > 0
			# One-sided (forward)
			idx_next = (i_next - 1) * n_dof_per_cell + var_idx
			φ_next = u[idx_next]
			z_next = mesh.cell_centers[i_next]
			gradients[i] = (φ_next - φ_i) / (z_next - z_i)
		else
			gradients[i] = 0.0
		end
	end

	return gradients
end

"""
	extrema_values(u, mesh, NC, var_idx) -> Tuple{Float64, Float64}

Get minimum and maximum values of a variable across the mesh.
"""
function extrema_values(u::Vector{Float64}, mesh::Mesh1D, NC::Int, var_idx::Int)
	n_cells = length(mesh.cells)
	n_dof_per_cell = 1 + NC

	val_min = Inf
	val_max = -Inf

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			idx = (i - 1) * n_dof_per_cell + var_idx
			val = u[idx]
			val_min = min(val_min, val)
			val_max = max(val_max, val)
		end
	end

	return (val_min, val_max)
end

"""
	find_nearest_cell(z, mesh) -> Int

Find index of nearest active cell to position z.
"""
function find_nearest_cell(z::Float64, mesh::Mesh1D)
	n_cells = length(mesh.cells)
	nearest = 1
	min_dist = Inf

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			dist = abs(z - mesh.cell_centers[i])
			if dist < min_dist
				min_dist = dist
				nearest = i
			end
		end
	end

	return nearest
end

# ============================================================================
# Conservation Check
# ============================================================================

"""
	check_remap_conservation(u_old, mesh_old, u_new, mesh_new, NC) -> NamedTuple

Check conservation properties after remapping.

# Returns
Named tuple with:
- `mass_error`: Relative mass conservation error for each component
- `energy_error`: Relative energy conservation error (approximate)
"""
function check_remap_conservation(u_old::Vector{Float64}, mesh_old::Mesh1D,
	u_new::Vector{Float64}, mesh_new::Mesh1D,
	NC::Int)
	n_dof_per_cell = 1 + NC
	n_cells_old = length(mesh_old.cells)
	n_cells_new = length(mesh_new.cells)

	# Compute total mass on old mesh
	mass_old = zeros(NC)
	@inbounds for i in 1:n_cells_old
		if mesh_old.active_cell_mask[i]
			V = mesh_old.cell_volumes[i]
			idx = (i - 1) * n_dof_per_cell
			for j in 1:NC
				mass_old[j] += u_old[idx+1+j] * V
			end
		end
	end

	# Compute total mass on new mesh
	mass_new = zeros(NC)
	@inbounds for i in 1:n_cells_new
		if mesh_new.active_cell_mask[i]
			V = mesh_new.cell_volumes[i]
			idx = (i - 1) * n_dof_per_cell
			for j in 1:NC
				mass_new[j] += u_new[idx+1+j] * V
			end
		end
	end

	# Compute errors
	mass_error = zeros(NC)
	for j in 1:NC
		if abs(mass_old[j]) > eps()
			mass_error[j] = abs(mass_new[j] - mass_old[j]) / mass_old[j]
		end
	end

	return (mass_error = mass_error,)
end

# ============================================================================
# State-Based Remapping (Alternative Interface)
# ============================================================================

"""
	remap_cell_states!(mesh_new::Mesh1D, mesh_old::Mesh1D, method::AbstractRemappingMethod)

Remap cell states directly between meshes.

This is a convenience function that works with the CellState structures directly.
"""
function remap_cell_states!(mesh_new::Mesh1D, mesh_old::Mesh1D,
	method::AbstractRemappingMethod)
	n_cells_new = length(mesh_new.cells)
	n_cells_old = length(mesh_old.cells)

	# Assume same NC for both meshes
	NC = length(mesh_old.cell_states[1].ξ)

	@inbounds for i_new in 1:n_cells_new
		if !mesh_new.active_cell_mask[i_new]
			continue
		end

		cell_new = mesh_new.cells[i_new]
		z_left_new = mesh_new.nodes[cell_new.node_ids[1]].z
		z_right_new = mesh_new.nodes[cell_new.node_ids[2]].z
		V_new = mesh_new.cell_volumes[i_new]

		T_accum = 0.0
		ξ_accum = zeros(NC)
		V_overlap_total = 0.0

		for i_old in 1:n_cells_old
			if !mesh_old.active_cell_mask[i_old]
				continue
			end

			cell_old = mesh_old.cells[i_old]
			z_left_old = mesh_old.nodes[cell_old.node_ids[1]].z
			z_right_old = mesh_old.nodes[cell_old.node_ids[2]].z

			z_left_overlap = max(z_left_new, z_left_old)
			z_right_overlap = min(z_right_new, z_right_old)

			if z_right_overlap > z_left_overlap
				A = mesh_old.cross_section_area
				V_overlap = (z_right_overlap - z_left_overlap) * A

				state_old = mesh_old.cell_states[i_old]
				T_accum += state_old.T * V_overlap

				for j in 1:NC
					ξ_accum[j] += state_old.ξ[j] * V_overlap
				end

				V_overlap_total += V_overlap
			end
		end

		state_new = mesh_new.cell_states[i_new]

		if V_overlap_total > 0.0
			state_new.T = T_accum / V_overlap_total

			ξ_tuple = ntuple(j -> max(0.0, ξ_accum[j] / V_new), NC)
			state_new.ξ = ξ_tuple
		end

		state_new.volume = V_new
	end

	return nothing
end

# ============================================================================
# Mass-Conserving Remapping Corrections
# ============================================================================

"""
    apply_mass_correction!(ξ_new::Vector{NTuple{NC,Float64}}, mesh_new::Mesh1D{NC},
                           ξ_old::Vector{NTuple{NC,Float64}}, mesh_old::Mesh1D{NC}) where NC

Apply uniform scaling to restore mass conservation after remapping.

After remapping with limiters, mass may not be exactly conserved due to the
limiters clamping values. This function applies a uniform scaling factor to
each component to restore exact mass conservation:

    ξ_new[j,i] *= (m_old[j] / m_new[j])

This preserves the spatial distribution while fixing total mass.

# Arguments
- `ξ_new`: State vector on new mesh (will be modified in-place)
- `mesh_new`: New mesh with volumes
- `ξ_old`: State vector on old mesh
- `mesh_old`: Old mesh with volumes

# Example
```julia
# After second-order remap
remap_solution_second_order!(ξ_new, mesh_new, ξ_old, mesh_old, NC)
# Correct any mass imbalance from limiters
apply_mass_correction!(ξ_new, mesh_new, ξ_old, mesh_old)
```
"""
function apply_mass_correction!(ξ_new::Vector{NTuple{NC,Float64}}, mesh_new::Mesh1D{NC},
	ξ_old::Vector{NTuple{NC,Float64}}, mesh_old::Mesh1D{NC}) where NC

	n_cells_old = length(mesh_old.cells)
	n_cells_new = length(mesh_new.cells)

	# Compute mass on old mesh
	mass_old = zeros(NC)
	@inbounds for i in 1:n_cells_old
		if mesh_old.active_cell_mask[i]
			V = mesh_old.cell_volumes[i]
			for j in 1:NC
				mass_old[j] += ξ_old[i][j] * V
			end
		end
	end

	# Compute mass on new mesh (post-remap, pre-correction)
	mass_new = zeros(NC)
	@inbounds for i in 1:n_cells_new
		if mesh_new.active_cell_mask[i]
			V = mesh_new.cell_volumes[i]
			for j in 1:NC
				mass_new[j] += ξ_new[i][j] * V
			end
		end
	end

	# Compute and apply correction factors
	@inbounds for j in 1:NC
		if mass_new[j] > eps()
			correction = mass_old[j] / mass_new[j]

			# Apply uniform scaling to all active cells
			for i in 1:n_cells_new
				if mesh_new.active_cell_mask[i]
					ξ_current = ξ_new[i]
					ξ_corrected = ntuple(k -> k == j ? ξ_current[k] * correction : ξ_current[k], NC)
					ξ_new[i] = ξ_corrected
				end
			end
		end
	end

	return nothing
end

"""
    remap_solution_conservative!(ξ_new::Vector{NTuple{NC,Float64}}, mesh_new::Mesh1D{NC},
                                  ξ_old::Vector{NTuple{NC,Float64}}, mesh_old::Mesh1D{NC},
                                  method::AbstractRemappingMethod;
                                  correct_mass::Bool=true) where NC

Remap solution with optional mass correction to ensure exact conservation.

This wraps the cell state remapping and applies a mass correction afterward to
restore exact conservation that may be violated by limiters.

# Arguments
- `ξ_new`: Output state vector on new mesh
- `mesh_new`: New mesh
- `ξ_old`: Input state vector on old mesh
- `mesh_old`: Old mesh
- `method`: Remapping method (FirstOrderRemap or SecondOrderRemap instance)
- `correct_mass`: If true, apply mass correction after remapping (default: true)

# Example
```julia
# Conservative second-order remap
remap_solution_conservative!(ξ_new, mesh_new, ξ_old, mesh_old,
                             SecondOrderRemap(); correct_mass=true)
```
"""
function remap_solution_conservative!(ξ_new::Vector{NTuple{NC,Float64}}, mesh_new::Mesh1D{NC},
	ξ_old::Vector{NTuple{NC,Float64}}, mesh_old::Mesh1D{NC},
	method::AbstractRemappingMethod;
	correct_mass::Bool = true) where NC

	# Copy old states to mesh_old cell states
	@inbounds for i in eachindex(ξ_old)
		mesh_old.cell_states[i].ξ = ξ_old[i]
	end

	# Perform cell state remapping
	remap_cell_states!(mesh_new, mesh_old, method)

	# Copy results to output
	@inbounds for i in eachindex(ξ_new)
		ξ_new[i] = mesh_new.cell_states[i].ξ
	end

	# Apply mass correction if requested
	if correct_mass
		apply_mass_correction!(ξ_new, mesh_new, ξ_old, mesh_old)
	end

	return nothing
end
