"""
Geometric quantity computations.

Functions for computing and updating cached geometric quantities.

Reference: Design Document Section 6.3
"""

# ============================================================================
# 1D Geometry Updates
# ============================================================================

"""
	update_geometry!(mesh::Mesh1D)

Recompute all cached geometric quantities from current node positions.

This should be called after any mesh motion.
"""
function update_geometry!(mesh::Mesh1D)
	n_cells_val = length(mesh.cells)
	n_faces_val = length(mesh.faces)
	A = mesh.cross_section_area

	# Update face positions and areas
	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		mesh.face_positions[f] = mesh.nodes[face.node_id].z
		mesh.face_areas[f] = A
	end

	# Update cell centers and volumes
	@inbounds for i in 1:n_cells_val
		cell = mesh.cells[i]
		z_left = mesh.nodes[cell.node_ids[1]].z
		z_right = mesh.nodes[cell.node_ids[2]].z

		mesh.cell_centers[i] = 0.5 * (z_left + z_right)
		mesh.cell_volumes[i] = abs(z_right - z_left) * A

		# Also update cell state volume
		mesh.cell_states[i].volume = mesh.cell_volumes[i]
	end

	return nothing
end

"""
	update_geometry_from_nodes!(mesh::Mesh1D)

Recompute all geometry (volumes, centers, face positions) from current node positions.

This is an alias for `update_geometry!` provided for clarity in ALE contexts
where geometry must be recomputed from the node positions stored in the state vector.

Used during ALE ODE RHS evaluation after unpacking state to ensure geometry
is consistent with the node positions.
"""
update_geometry_from_nodes!(mesh::Mesh1D) = update_geometry!(mesh)

"""
	cell_volume(mesh::Mesh1D, i::Int) -> Float64

Get the volume of the i-th cell.
"""
@inline cell_volume(mesh::Mesh1D, i::Int) = mesh.cell_volumes[i]

"""
	face_area(mesh::Mesh1D, f::Int) -> Float64

Get the area of the f-th face (cross-sectional area in 1D).
"""
@inline face_area(mesh::Mesh1D, f::Int) = mesh.face_areas[f]

"""
	cell_center(mesh::Mesh1D, i::Int) -> Float64

Get the center position of the i-th cell.
"""
@inline cell_center(mesh::Mesh1D, i::Int) = mesh.cell_centers[i]

"""
	face_position(mesh::Mesh1D, f::Int) -> Float64

Get the position of the f-th face.
"""
@inline face_position(mesh::Mesh1D, f::Int) = mesh.face_positions[f]

"""
	cell_size(mesh::Mesh1D, i::Int) -> Float64

Get the size (length) of the i-th cell.
"""
@inline function cell_size(mesh::Mesh1D, i::Int)
	cell = mesh.cells[i]
	z_left = mesh.nodes[cell.node_ids[1]].z
	z_right = mesh.nodes[cell.node_ids[2]].z
	return abs(z_right - z_left)
end

"""
	face_distance(mesh::Mesh1D, f::Int) -> Float64

Get the distance between cell centers across face f.
Returns the distance from left cell center to right cell center for interior faces,
or twice the distance from boundary to cell center for boundary faces.
"""
function face_distance(mesh::Mesh1D, f::Int)
	face = mesh.faces[f]
	left_id, right_id = face.cell_ids

	if left_id == 0
		# Left boundary: distance from boundary to right cell center
		return 2.0 * abs(mesh.cell_centers[right_id] - mesh.face_positions[f])
	elseif right_id == 0
		# Right boundary: distance from left cell center to boundary
		return 2.0 * abs(mesh.face_positions[f] - mesh.cell_centers[left_id])
	else
		# Interior face
		return abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
	end
end

"""
	total_volume(mesh::Mesh1D) -> Float64

Compute the total volume of the domain.
"""
function total_volume(mesh::Mesh1D)
	vol = 0.0
	@inbounds for i in 1:length(mesh.cell_volumes)
		if mesh.active_cell_mask[i]
			vol += mesh.cell_volumes[i]
		end
	end
	return vol
end

"""
	domain_bounds(mesh::Mesh1D) -> Tuple{Float64,Float64}

Get the domain bounds (z_min, z_max).
"""
function domain_bounds(mesh::Mesh1D)
	z_min = minimum(n.z for n in mesh.nodes)
	z_max = maximum(n.z for n in mesh.nodes)
	return (z_min, z_max)
end

# ============================================================================
# 2D Geometry Updates (Placeholder)
# ============================================================================

"""
	update_geometry!(mesh::Mesh2D)

Recompute all cached geometric quantities for a 2D mesh.
"""
function update_geometry!(mesh::Mesh2D)
	error("2D geometry update not yet implemented")
end

"""
	cell_volume(mesh::Mesh2D, i::Int) -> Float64

Get the volume of the i-th cell (includes 2πr factor for axisymmetric).
"""
@inline cell_volume(mesh::Mesh2D, i::Int) = mesh.cell_volumes[i]

"""
	face_area(mesh::Mesh2D, f::Int) -> Float64

Get the area of the f-th face (includes 2πr factor for axisymmetric).
"""
@inline face_area(mesh::Mesh2D, f::Int) = mesh.face_areas[f]

"""
	cell_center(mesh::Mesh2D, i::Int) -> Tuple{Float64,Float64}

Get the center position (r, z) of the i-th cell.
"""
@inline cell_center(mesh::Mesh2D, i::Int) = mesh.cell_centers[i]

# ============================================================================
# Mesh Quality Metrics
# ============================================================================

"""
	aspect_ratio(mesh::Mesh1D, i::Int) -> Float64

Compute aspect ratio of cell i relative to its neighbors.
Returns max(Δz_i/Δz_{i-1}, Δz_i/Δz_{i+1}).
"""
function aspect_ratio(mesh::Mesh1D, i::Int)
	n_cells_val = length(mesh.cells)
	Δz_i = cell_size(mesh, i)

	ratio = 1.0
	if i > 1
		Δz_prev = cell_size(mesh, i-1)
		ratio = max(ratio, Δz_i / Δz_prev, Δz_prev / Δz_i)
	end
	if i < n_cells_val
		Δz_next = cell_size(mesh, i+1)
		ratio = max(ratio, Δz_i / Δz_next, Δz_next / Δz_i)
	end

	return ratio
end

"""
	max_aspect_ratio(mesh::Mesh1D) -> Float64

Get the maximum aspect ratio across all cells.
"""
function max_aspect_ratio(mesh::Mesh1D)
	max_ar = 1.0
	for i in 1:length(mesh.cells)
		if mesh.active_cell_mask[i]
			max_ar = max(max_ar, aspect_ratio(mesh, i))
		end
	end
	return max_ar
end

"""
	min_cell_size(mesh::Mesh1D) -> Float64

Get the minimum cell size in the mesh.
"""
function min_cell_size(mesh::Mesh1D)
	min_size = Inf
	for i in 1:length(mesh.cells)
		if mesh.active_cell_mask[i]
			min_size = min(min_size, cell_size(mesh, i))
		end
	end
	return min_size
end

"""
	max_cell_size(mesh::Mesh1D) -> Float64

Get the maximum cell size in the mesh.
"""
function max_cell_size(mesh::Mesh1D)
	max_size = 0.0
	for i in 1:length(mesh.cells)
		if mesh.active_cell_mask[i]
			max_size = max(max_size, cell_size(mesh, i))
		end
	end
	return max_size
end

# ============================================================================
# Boundary Access Functions
# ============================================================================

"""
	boundary_cell_id(mesh::Mesh1D, face_id::Int) -> Int

Get the cell ID adjacent to a boundary face.
For a boundary face, returns the non-zero cell ID.
"""
@inline function boundary_cell_id(mesh::Mesh1D, face_id::Int)
	face = mesh.faces[face_id]
	left_id, right_id = face.cell_ids
	# One of these should be 0 for a boundary face
	return left_id == 0 ? right_id : left_id
end

"""
	face_area_value(mesh::Mesh1D, f::Int) -> Float64

Alias for face_area for consistency in boundary flux computation.
"""
@inline face_area_value(mesh::Mesh1D, f::Int) = mesh.face_areas[f]

"""
	boundary_cell_id(mesh::Mesh2D, face_id::Int) -> Int

Get the cell ID adjacent to a boundary face for 2D mesh.
"""
@inline function boundary_cell_id(mesh::Mesh2D, face_id::Int)
	face = mesh.faces[face_id]
	left_id, right_id = face.cell_ids
	return left_id == 0 ? right_id : left_id
end

"""
	face_area_value(mesh::Mesh2D, f::Int) -> Float64

Alias for face_area for consistency in boundary flux computation.
"""
@inline face_area_value(mesh::Mesh2D, f::Int) = mesh.face_areas[f]

# ============================================================================
# Geometry Cache for ALE ODE Evaluation
# ============================================================================

"""
	GeometryCache

Pre-allocated arrays for geometry computation during ALE ODE evaluation.

This struct holds work arrays that are updated from the node positions in the
state vector at each ODE RHS evaluation. Using pre-allocated arrays avoids
allocations during the time-critical RHS computation.

# Fields
- `volumes::Vector{Float64}`: Cell volumes (length n_cells)
- `centers::Vector{Float64}`: Cell center positions (length n_cells)
- `face_positions::Vector{Float64}`: Face positions (length n_faces)
- `face_distances::Vector{Float64}`: Distances between cell centers across faces (length n_faces)

# Usage
```julia
cache = GeometryCache(n_cells, n_faces)
update_geometry_cache!(cache, u, sv, mesh)
# Now cache.volumes, cache.centers, etc. contain current geometry
```
"""
struct GeometryCache
	volumes::Vector{Float64}
	centers::Vector{Float64}
	face_positions::Vector{Float64}
	face_distances::Vector{Float64}
end

"""
	GeometryCache(n_cells::Int, n_faces::Int)

Create a GeometryCache with pre-allocated arrays.
"""
function GeometryCache(n_cells::Int, n_faces::Int)
	GeometryCache(
		Vector{Float64}(undef, n_cells),
		Vector{Float64}(undef, n_cells),
		Vector{Float64}(undef, n_faces),
		Vector{Float64}(undef, n_faces),
	)
end

"""
	GeometryCache(mesh::Mesh1D)

Create a GeometryCache sized for the given mesh.
"""
GeometryCache(mesh::Mesh1D) = GeometryCache(length(mesh.cells), length(mesh.faces))

"""
	resize_geometry_cache!(cache::GeometryCache, n_cells::Int, n_faces::Int)

Resize a GeometryCache for a new mesh size.

This function is used during adaptive mesh operations (h-refinement, h-coarsening,
surface cell depletion) when the number of cells/faces changes.

# Arguments
- `cache::GeometryCache`: Cache to resize
- `n_cells::Int`: New number of cells
- `n_faces::Int`: New number of faces (typically n_cells + 1 for 1D)

# Notes
Arrays are resized in-place. Previous values are lost and should be
recomputed after calling this function.
"""
function resize_geometry_cache!(cache::GeometryCache, n_cells::Int, n_faces::Int)
	resize!(cache.volumes, n_cells)
	resize!(cache.centers, n_cells)
	resize!(cache.face_positions, n_faces)
	resize!(cache.face_distances, n_faces)
	return nothing
end

# ============================================================================
# Allocation-Free Geometry Computation from State Vector
# ============================================================================

# Import the state vector accessor functions we need
# These are defined in TimeIntegration.state_vector but we need them here
# We'll use a duck-typing approach - any sv with get_node_position will work

"""
	compute_geometry_from_state!(volumes, centers, face_positions, 
								 u, sv, A)

Compute geometry from node positions in state vector into pre-allocated arrays.

This is the allocation-free path used during ODE RHS evaluation. It computes
cell volumes, cell centers, and face positions directly from the node positions
stored in the state vector `u`.

# Arguments
- `volumes::Vector{Float64}`: Pre-allocated vector for cell volumes (length n_cells)
- `centers::Vector{Float64}`: Pre-allocated vector for cell centers (length n_cells)
- `face_positions::Vector{Float64}`: Pre-allocated vector for face positions (length n_nodes)
- `u::AbstractVector`: State vector containing node positions
- `sv`: StateVectorALE descriptor (must support `get_node_position`)
- `A::Float64`: Cross-sectional area

# Notes
This function does NOT modify the mesh. It computes geometry into work arrays
that can be used for flux calculations without touching the mesh struct.

For 1D meshes, face positions equal node positions (one face per node).
"""
function compute_geometry_from_state!(volumes::AbstractVector{Float64},
	centers::AbstractVector{Float64},
	face_positions::AbstractVector{Float64},
	u::AbstractVector,
	sv,
	A::Float64)
	n = sv.n_cells

	# Compute cell volumes and centers from node positions
	@inbounds for i in 1:n
		z_left = u[sv.z_offset+i]
		z_right = u[sv.z_offset+i+1]
		volumes[i] = abs(z_right - z_left) * A
		centers[i] = 0.5 * (z_left + z_right)
	end

	# Face positions = node positions in 1D
	@inbounds for i in 1:sv.n_nodes
		face_positions[i] = u[sv.z_offset+i]
	end

	return nothing
end

"""
	compute_face_distances_from_state!(distances, centers, face_positions, mesh)

Compute distances between cell centers across faces.

For interior faces: distance = |center_R - center_L|
For boundary faces: distance = 2 × |face_position - adjacent_center|

The factor of 2 for boundary faces accounts for the "ghost cell" approach where
the gradient at a boundary uses a reflected cell center.

# Arguments
- `distances::Vector{Float64}`: Pre-allocated output vector (length n_faces)
- `centers::Vector{Float64}`: Cell center positions (from compute_geometry_from_state!)
- `face_positions::Vector{Float64}`: Face positions (from compute_geometry_from_state!)
- `mesh::Mesh1D`: Mesh for topology information (cell_ids of faces)
"""
function compute_face_distances_from_state!(distances::AbstractVector{Float64},
	centers::AbstractVector{Float64},
	face_positions::AbstractVector{Float64},
	mesh::Mesh1D)
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0
			# Bottom boundary (z=0): distance from boundary to cell center
			distances[f] = 2.0 * abs(centers[right_id] - face_positions[f])
		elseif right_id == 0
			# Top boundary (z=L): distance from cell center to boundary
			distances[f] = 2.0 * abs(face_positions[f] - centers[left_id])
		else
			# Interior face: distance between adjacent cell centers
			distances[f] = abs(centers[right_id] - centers[left_id])
		end
	end

	return nothing
end

"""
	update_geometry_cache!(cache::GeometryCache, u::AbstractVector, sv, mesh::Mesh1D)

Update all geometry arrays in cache from current state vector.

This is the main entry point for geometry computation during ALE ODE RHS evaluation.
It computes volumes, centers, face positions, and face distances in a single call.

# Arguments
- `cache::GeometryCache`: Pre-allocated cache to update
- `u::AbstractVector`: ALE state vector containing node positions
- `sv`: StateVectorALE descriptor
- `mesh::Mesh1D`: Mesh for cross-sectional area and topology

# Example
```julia
# In ODE RHS evaluation:
update_geometry_cache!(geom_cache, u, sv, mesh)
# Now use geom_cache.volumes, geom_cache.centers, etc. for flux computation
```
"""
function update_geometry_cache!(cache::GeometryCache, u::AbstractVector,
	sv, mesh::Mesh1D)
	A = mesh.cross_section_area

	# Compute volumes, centers, face positions from node positions in state
	compute_geometry_from_state!(cache.volumes, cache.centers,
		cache.face_positions, u, sv, A)

	# Compute face distances using the computed centers and positions
	compute_face_distances_from_state!(cache.face_distances, cache.centers,
		cache.face_positions, mesh)

	return nothing
end

# ============================================================================
# GCL Verification Utilities
# ============================================================================

"""
	compute_gcl_volume_rate(w::AbstractVector{Float64}, mesh::Mesh1D) -> Vector{Float64}

Compute the volume change rate for each cell from node velocities using GCL formula.

For 1D: dVᵢ/dt = (w_{i+1} - wᵢ) × A

This should match the actual volume rate computed from dz/dt.

# Arguments
- `w::Vector{Float64}`: Node velocities (length n_nodes)
- `mesh::Mesh1D`: Mesh for cross-sectional area

# Returns
Vector of volume change rates (length n_cells)
"""
function compute_gcl_volume_rate(w::AbstractVector{Float64}, mesh::Mesh1D)
	n_cells = length(mesh.cells)
	A = mesh.cross_section_area

	dV_dt = Vector{Float64}(undef, n_cells)
	@inbounds for i in 1:n_cells
		# GCL: dV/dt = (w_right - w_left) * A
		dV_dt[i] = (w[i+1] - w[i]) * A
	end

	return dV_dt
end

"""
	verify_gcl_compliance(volumes_old, volumes_new, w, dt, mesh; atol=1e-12) -> Bool

Verify that volume changes satisfy the GCL to within tolerance.

Checks: (V_new - V_old) / dt ≈ (w_right - w_left) * A for each cell.

# Arguments
- `volumes_old::Vector{Float64}`: Cell volumes before motion
- `volumes_new::Vector{Float64}`: Cell volumes after motion
- `w::Vector{Float64}`: Node velocities used for motion
- `dt::Float64`: Time step
- `mesh::Mesh1D`: Mesh for cross-sectional area
- `atol::Float64`: Absolute tolerance for comparison

# Returns
`true` if GCL is satisfied to tolerance, `false` otherwise.
"""
function verify_gcl_compliance(volumes_old::AbstractVector{Float64},
	volumes_new::AbstractVector{Float64},
	w::AbstractVector{Float64},
	dt::Float64,
	mesh::Mesh1D;
	atol::Float64 = 1e-12)
	n_cells = length(mesh.cells)
	A = mesh.cross_section_area

	@inbounds for i in 1:n_cells
		dV_dt_actual = (volumes_new[i] - volumes_old[i]) / dt
		dV_dt_gcl = (w[i+1] - w[i]) * A

		if abs(dV_dt_actual - dV_dt_gcl) > atol
			return false
		end
	end

	return true
end
