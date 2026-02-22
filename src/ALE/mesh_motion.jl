"""
Mesh motion algorithms.

Functions for moving mesh nodes and updating geometry.

Reference: Design Document Section 3.5, 7.3
"""

using ..Geometry: Mesh1D, Node1D, update_geometry!, cell_volume

# ============================================================================
# Basic Mesh Motion
# ============================================================================

"""
	move_mesh!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
			   use_gcl_volumes::Bool=true)

Move mesh nodes according to mesh velocity and update geometry.

# Arguments
- `mesh`: Mesh to modify
- `w`: Mesh velocity at each node [m/s]
- `dt`: Time step [s]
- `use_gcl_volumes`: If true, use GCL-compliant volume update (default: true)

# Notes
- Node positions are updated: z_new = z_old + dt * w
- Geometry cache is updated after motion
- If `use_gcl_volumes=true`, volumes are computed from GCL flux (exact conservation)
- If `use_gcl_volumes=false`, volumes are recomputed from geometry (may have roundoff)
"""
function move_mesh!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
	use_gcl_volumes::Bool = true)
	n_nodes = length(mesh.nodes)

	@assert length(w) == n_nodes "Mesh velocity vector length must match number of nodes"

	# Store old volumes for GCL update
	if use_gcl_volumes
		V_old = copy(mesh.cell_volumes)
	end

	# Update node positions
	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]
		z_new = node.z + dt * w[i]

		# Create new node with updated position (Node1D is immutable)
		mesh.nodes[i] = Node1D(
			node.id,
			z_new,
			node.z_ref,  # Reference position unchanged
			node.is_boundary,
			node.boundary_id,
		)
	end

	# Update mesh time
	mesh.mesh_time += dt

	# Update geometry
	if use_gcl_volumes
		# Update face positions and cell centers from geometry
		update_face_positions!(mesh)
		update_cell_centers!(mesh)

		# Update volumes using GCL formula (exact conservation)
		A = mesh.cross_section_area
		@inbounds for i in eachindex(mesh.cells)
			if !mesh.active_cell_mask[i]
				continue
			end

			cell = mesh.cells[i]
			node_left, node_right = cell.node_ids

			w_left = w[node_left]
			w_right = w[node_right]

			# GCL volume update
			V_new = V_old[i] + dt * (w_right - w_left) * A
			mesh.cell_volumes[i] = V_new
			mesh.cell_states[i].volume = V_new
		end
	else
		# Full geometry update from node positions
		update_geometry!(mesh)
	end

	# Store node velocities in mesh
	mesh.node_velocities .= w

	return nothing
end

"""
	update_face_positions!(mesh::Mesh1D)

Update face positions from current node positions.
"""
function update_face_positions!(mesh::Mesh1D)
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		mesh.face_positions[f] = mesh.nodes[face.node_id].z
	end
	return nothing
end

"""
	update_cell_centers!(mesh::Mesh1D)

Update cell center positions from current node positions.
"""
function update_cell_centers!(mesh::Mesh1D)
	@inbounds for i in eachindex(mesh.cells)
		cell = mesh.cells[i]
		z_left = mesh.nodes[cell.node_ids[1]].z
		z_right = mesh.nodes[cell.node_ids[2]].z
		mesh.cell_centers[i] = 0.5 * (z_left + z_right)
	end
	return nothing
end

# ============================================================================
# Mesh Quality Checks
# ============================================================================

"""
	check_mesh_quality(mesh::Mesh1D) -> NamedTuple

Check mesh quality metrics after motion.

# Returns
Named tuple with:
- `valid`: true if all cells have positive volume
- `min_volume`: minimum cell volume
- `max_aspect_ratio`: maximum ratio of adjacent cell sizes
- `tangled_cells`: indices of cells with non-positive volume
"""
function check_mesh_quality(mesh::Mesh1D)
	min_volume = Inf
	max_aspect_ratio = 0.0
	tangled_cells = Int[]

	prev_size = 0.0

	@inbounds for i in eachindex(mesh.cells)
		if !mesh.active_cell_mask[i]
			continue
		end

		V = mesh.cell_volumes[i]

		# Check for tangling (negative or zero volume)
		if V <= 0.0
			push!(tangled_cells, i)
		end

		min_volume = min(min_volume, V)

		# Cell size for aspect ratio
		cell = mesh.cells[i]
		z_left = mesh.nodes[cell.node_ids[1]].z
		z_right = mesh.nodes[cell.node_ids[2]].z
		cell_size = abs(z_right - z_left)

		# Aspect ratio with previous cell
		if prev_size > 0.0 && cell_size > 0.0
			ratio = max(cell_size / prev_size, prev_size / cell_size)
			max_aspect_ratio = max(max_aspect_ratio, ratio)
		end

		prev_size = cell_size
	end

	return (
		valid = isempty(tangled_cells) && min_volume > 0.0,
		min_volume = min_volume,
		max_aspect_ratio = max_aspect_ratio,
		tangled_cells = tangled_cells,
	)
end

"""
	mesh_is_valid(mesh::Mesh1D) -> Bool

Quick check if mesh has valid geometry (all positive volumes).
"""
function mesh_is_valid(mesh::Mesh1D)
	@inbounds for i in eachindex(mesh.cells)
		if mesh.active_cell_mask[i] && mesh.cell_volumes[i] <= 0.0
			return false
		end
	end
	return true
end

# ============================================================================
# Safe Mesh Motion (with quality check)
# ============================================================================

"""
	move_mesh_safe!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
					max_aspect_ratio::Float64=10.0) -> Bool

Move mesh with quality check. Returns false if motion would create invalid mesh.

# Arguments
- `mesh`: Mesh to modify
- `w`: Mesh velocity at each node [m/s]
- `dt`: Time step [s]
- `max_aspect_ratio`: Maximum allowed ratio of adjacent cell sizes

# Returns
`true` if mesh was moved successfully, `false` if motion was rejected
"""
function move_mesh_safe!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
	max_aspect_ratio::Float64 = 10.0)
	# Store original state for potential rollback
	original_nodes = copy(mesh.nodes)
	original_volumes = copy(mesh.cell_volumes)
	original_time = mesh.mesh_time

	# Attempt motion
	move_mesh!(mesh, w, dt)

	# Check quality
	quality = check_mesh_quality(mesh)

	if !quality.valid || quality.max_aspect_ratio > max_aspect_ratio
		# Rollback
		mesh.nodes .= original_nodes
		mesh.cell_volumes .= original_volumes
		mesh.mesh_time = original_time
		update_geometry!(mesh)
		return false
	end

	return true
end

# ============================================================================
# Utilities
# ============================================================================

"""
	mesh_displacement(mesh::Mesh1D) -> Vector{Float64}

Get displacement of each node from reference position.
"""
function mesh_displacement(mesh::Mesh1D)
	n_nodes = length(mesh.nodes)
	displacement = Vector{Float64}(undef, n_nodes)

	@inbounds for i in 1:n_nodes
		displacement[i] = mesh.nodes[i].z - mesh.nodes[i].z_ref
	end

	return displacement
end

"""
	max_mesh_displacement(mesh::Mesh1D) -> Float64

Get maximum absolute displacement from reference configuration.
"""
function max_mesh_displacement(mesh::Mesh1D)
	max_disp = 0.0
	@inbounds for node in mesh.nodes
		disp = abs(node.z - node.z_ref)
		max_disp = max(max_disp, disp)
	end
	return max_disp
end

"""
	reset_mesh_to_reference!(mesh::Mesh1D)

Reset mesh nodes to their reference positions.
"""
function reset_mesh_to_reference!(mesh::Mesh1D)
	@inbounds for i in eachindex(mesh.nodes)
		node = mesh.nodes[i]
		mesh.nodes[i] = Node1D(
			node.id,
			node.z_ref,  # Reset to reference
			node.z_ref,
			node.is_boundary,
			node.boundary_id,
		)
	end

	mesh.mesh_time = 0.0
	fill!(mesh.node_velocities, 0.0)
	update_geometry!(mesh)

	return nothing
end

# ============================================================================
# GCL-Enforced Mesh Motion
# ============================================================================

"""
    move_mesh_gcl!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
                   check_gcl::Bool=true, gcl_tol::Float64=1e-10,
                   warn_on_violation::Bool=true) -> Float64

Move mesh nodes and update volumes using GCL-compliant formula with optional verification.

This is a wrapper around `move_mesh!` that enforces GCL compliance and optionally
checks for violations.

# Arguments
- `mesh`: Mesh to modify
- `w`: Node velocities [m/s]
- `dt`: Time step [s]
- `check_gcl`: If true, compute and return GCL error
- `gcl_tol`: Tolerance for GCL warning
- `warn_on_violation`: If true, warn when GCL error exceeds tolerance

# Returns
The GCL error (should be ≈ machine precision if volumes are updated correctly).
Returns 0.0 if `check_gcl=false`.

# GCL Compliance
Uses the formula: V_new = V_old + dt × (w_right - w_left) × A
instead of geometric recomputation, ensuring exact GCL satisfaction.

# Example
```julia
gcl_err = move_mesh_gcl!(mesh, w, dt; gcl_tol=1e-10)
if gcl_err > 1e-8
    @warn "Large GCL error detected" error=gcl_err
end
```
"""
function move_mesh_gcl!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
	check_gcl::Bool = true, gcl_tol::Float64 = 1e-10,
	warn_on_violation::Bool = true)

	# Store old volumes for GCL check if requested
	V_old = check_gcl ? copy(mesh.cell_volumes) : Float64[]

	# Move mesh with GCL volume update
	move_mesh!(mesh, w, dt; use_gcl_volumes = true)

	# Check GCL compliance if requested
	gcl_err = 0.0
	if check_gcl
		# Import gcl_error_post_motion from parent module
		gcl_err = gcl_error_post_motion(mesh, V_old, w, dt)

		if warn_on_violation && gcl_err > gcl_tol
			@warn "GCL error exceeds tolerance" gcl_error = gcl_err tolerance = gcl_tol
		end
	end

	return gcl_err
end
