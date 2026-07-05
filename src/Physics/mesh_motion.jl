"""
Mesh motion algorithms.

Functions for moving mesh nodes and updating geometry.
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
- `use_gcl_volumes`: If true, write cell volumes from the GCL flux increment
  ΔV = A·(w_R − w_L)·dt instead of recomputing from the moved node positions

# Notes
- Node positions are updated: z_new = z_old + dt * w
- Geometry cache is updated after motion
- The two volume updates agree to roundoff: in 1D with V = A·|Δz| the
  geometric recomputation already satisfies the discrete GCL as an algebraic
  identity (this is how the production RHS computes volumes from state every
  call). The flag only selects which equivalent form this utility writes.
"""
function move_mesh!(
    mesh::Mesh1D,
    w::Vector{Float64},
    dt::Float64;
    use_gcl_volumes::Bool = true,
)
    n_nodes = length(mesh.nodes)

    @assert length(w) == n_nodes "Mesh velocity vector length must match number of nodes"

    # Update node positions
    @inbounds for i = 1:n_nodes
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
        # Update face positions, face areas, and cell centers from geometry
        update_face_positions!(mesh)
        update_face_areas!(mesh)
        update_cell_centers!(mesh)

        # Update volumes analytically from the current global cross-section.
        # For variable cross-section, carrying V_old across an A(t) change would
        # preserve the old radial area. The geometric state is V = A(t) * L.
        A = mesh.cross_section_area
        @inbounds for i in eachindex(mesh.cells)
            if !mesh.active_cell_mask[i]
                continue
            end

            cell = mesh.cells[i]
            node_left, node_right = cell.node_ids
            L_new = abs(mesh.nodes[node_right].z - mesh.nodes[node_left].z)
            V_new = A * L_new
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
	update_face_areas!(mesh::Mesh1D)

Refresh all cached face areas from the current global cross-sectional area.
"""
function update_face_areas!(mesh::Mesh1D)
    A = mesh.cross_section_area
    @inbounds for f in eachindex(mesh.faces)
        mesh.face_areas[f] = A
    end
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

Cell widths are computed **signed** (`z_right − z_left` in node order): the
cached `cell_volumes` pass through `abs()`, so an inverted (tangled) cell would
otherwise look like a perfectly valid thin cell and rollback could never fire.
Adjacency for the aspect ratio is taken in spatial order (sorted by cell
center), not index order, which diverges after h-refinement.

# Returns
Named tuple with:
- `valid`: true if all cells have positive signed width
- `min_volume`: minimum signed cell volume
- `max_aspect_ratio`: maximum ratio of spatially adjacent cell sizes
- `tangled_cells`: indices of cells with non-positive signed width
"""
function check_mesh_quality(mesh::Mesh1D)
    min_volume = Inf
    max_aspect_ratio = 0.0
    tangled_cells = Int[]
    A = mesh.cross_section_area

    order = [i for i in eachindex(mesh.cells) if mesh.active_cell_mask[i]]
    sort!(order, by = i -> mesh.cell_centers[i])

    prev_size = 0.0

    @inbounds for i in order
        cell = mesh.cells[i]
        z_left = mesh.nodes[cell.node_ids[1]].z
        z_right = mesh.nodes[cell.node_ids[2]].z
        width = z_right - z_left  # signed: ≤ 0 ⟺ inverted or collapsed cell

        if width <= 0.0
            push!(tangled_cells, i)
        end

        min_volume = min(min_volume, width * A)

        # Aspect ratio with the spatially previous cell
        if prev_size > 0.0 && width > 0.0
            ratio = max(width / prev_size, prev_size / width)
            max_aspect_ratio = max(max_aspect_ratio, ratio)
        end

        prev_size = width
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

Quick check if mesh has valid geometry: every active cell has positive
**signed** width (z monotone in node order). The cached `cell_volumes` are
`abs()`-based and cannot detect inversion.
"""
function mesh_is_valid(mesh::Mesh1D)
    @inbounds for i in eachindex(mesh.cells)
        if mesh.active_cell_mask[i]
            cell = mesh.cells[i]
            width = mesh.nodes[cell.node_ids[2]].z - mesh.nodes[cell.node_ids[1]].z
            if width <= 0.0
                return false
            end
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
function move_mesh_safe!(
    mesh::Mesh1D,
    w::Vector{Float64},
    dt::Float64;
    max_aspect_ratio::Float64 = 10.0,
)
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

    @inbounds for i = 1:n_nodes
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
