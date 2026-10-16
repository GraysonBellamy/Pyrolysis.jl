"""
h-refinement (cell splitting).

h-refinement splits cells where error indicators are high. Each refined cell
is split into two children, with state interpolated conservatively.

## Conservation Property
When cell p splits into children a and b:
  V_p · φ_p = V_a · φ_a + V_b · φ_b

For uniform splitting with linear interpolation:
  φ_a = φ_p - (∇φ · Δx)/4
  φ_b = φ_p + (∇φ · Δx)/4
"""

using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, update_geometry!, n_cells

# ============================================================================
# Single Cell Refinement
# ============================================================================

"""
	refine_cell!(mesh::Mesh1D{NC}, cell_id::Int) -> Tuple{Int,Int}

Split a single cell into two children.

# Arguments
- `mesh`: The mesh to modify
- `cell_id`: Index of cell to refine

# Returns
- Tuple of (child_a_id, child_b_id)

# Algorithm
1. Create new node at cell center
2. Create two child cells (left and right)
3. Create new interior face between children
4. Update face connectivity
5. Interpolate state conservatively to children
6. Mark parent as inactive
"""
function refine_cell!(mesh::Mesh1D{NC}, cell_id::Int) where {NC}
    cell = mesh.cells[cell_id]

    # Check if already refined
    if cell.children != (0, 0)
        error("Cell $cell_id is already refined")
    end

    # Get parent geometry
    left_node_id = cell.node_ids[1]
    right_node_id = cell.node_ids[2]
    left_face_id = cell.face_ids[1]
    right_face_id = cell.face_ids[2]

    z_left = mesh.nodes[left_node_id].z
    z_right = mesh.nodes[right_node_id].z
    z_center = 0.5 * (z_left + z_right)
    z_ref_center = 0.5 * (mesh.nodes[left_node_id].z_ref + mesh.nodes[right_node_id].z_ref)

    # Parent state
    parent_state = mesh.cell_states[cell_id]
    V_parent = mesh.cell_volumes[cell_id]

    # ========================================================================
    # 1. Create new node at center
    # ========================================================================
    new_node_id = length(mesh.nodes) + 1
    new_node = Node1D(new_node_id, z_center, z_ref_center, false, 0)
    push!(mesh.nodes, new_node)
    push!(mesh.node_velocities, 0.0)

    # ========================================================================
    # 2. Create two child cells
    # ========================================================================
    new_level = cell.level + 1

    # Child A (left)
    child_a_id = length(mesh.cells) + 1
    # Child B (right)
    child_b_id = child_a_id + 1

    # New face between children
    new_face_id = length(mesh.faces) + 1

    # Child A: left boundary → center
    child_a = Cell1D(
        child_a_id,
        (left_node_id, new_node_id),
        (left_face_id, new_face_id),
        new_level,
        cell_id,
        (0, 0),
    )

    # Child B: center → right boundary
    child_b = Cell1D(
        child_b_id,
        (new_node_id, right_node_id),
        (new_face_id, right_face_id),
        new_level,
        cell_id,
        (0, 0),
    )

    push!(mesh.cells, child_a)
    push!(mesh.cells, child_b)

    # ========================================================================
    # 3. Create new interior face
    # ========================================================================
    new_face = Face1D(new_face_id, new_node_id, (child_a_id, child_b_id), 0)
    push!(mesh.faces, new_face)

    # ========================================================================
    # 4. Update existing face connectivity
    # ========================================================================
    # Left face: update right cell from parent to child_a
    left_face = mesh.faces[left_face_id]
    mesh.faces[left_face_id] = Face1D(
        left_face.id,
        left_face.node_id,
        (left_face.cell_ids[1], child_a_id),
        left_face.is_boundary,
        left_face.boundary_id,
    )

    # Right face: update left cell from parent to child_b
    right_face = mesh.faces[right_face_id]
    mesh.faces[right_face_id] = Face1D(
        right_face.id,
        right_face.node_id,
        (child_b_id, right_face.cell_ids[2]),
        right_face.is_boundary,
        right_face.boundary_id,
    )

    # ========================================================================
    # 5. Update parent cell with children info
    # ========================================================================
    mesh.cells[cell_id] = Cell1D(
        cell.id,
        cell.node_ids,
        cell.face_ids,
        cell.level,
        cell.parent,
        (child_a_id, child_b_id),
    )

    # ========================================================================
    # 6. Interpolate state conservatively to children
    # ========================================================================
    # Compute gradient for better interpolation
    T_grad = estimate_gradient(mesh, cell_id, parent_state.T, :T)
    ξ_grads = ntuple(j -> estimate_gradient(mesh, cell_id, parent_state.ξ[j], j), NC)

    Δz = z_center - mesh.cell_centers[cell_id]  # Distance from parent center to split

    # Child volumes (assuming uniform split)
    V_a = V_parent / 2
    V_b = V_parent / 2

    # Temperature: interpolate linearly
    T_a = parent_state.T - T_grad * abs(Δz) / 2
    T_b = parent_state.T + T_grad * abs(Δz) / 2

    # Concentrations: ensure conservation (extensive property)
    # Total mass must be preserved: V_p * ξ_p = V_a * ξ_a + V_b * ξ_b
    # With linear interpolation and equal volumes: ξ_a = ξ_b = ξ_p conserves mass
    ξ_a = parent_state.ξ
    ξ_b = parent_state.ξ

    # Create child states
    child_a_state = CellState{NC}(T_a, ξ_a, V_a, V_a)
    child_b_state = CellState{NC}(T_b, ξ_b, V_b, V_b)

    push!(mesh.cell_states, child_a_state)
    push!(mesh.cell_states, child_b_state)

    # ========================================================================
    # 7. Update geometry caches
    # ========================================================================
    push!(mesh.cell_centers, 0.5 * (z_left + z_center))  # child_a center
    push!(mesh.cell_centers, 0.5 * (z_center + z_right))  # child_b center
    push!(mesh.cell_volumes, V_a)
    push!(mesh.cell_volumes, V_b)
    push!(mesh.face_positions, z_center)
    push!(mesh.face_areas, mesh.cross_section_area)

    # Preserve initial dry-solid mass extensively across the split. With uniform
    # split (V_a = V_b = V_parent/2) the children inherit half each. A
    # non-uniform-split formula could refine this further.
    if !isempty(mesh.m_dry_solid_initial)
        m_parent = mesh.m_dry_solid_initial[cell_id]
        push!(mesh.m_dry_solid_initial, 0.5 * m_parent)
        push!(mesh.m_dry_solid_initial, 0.5 * m_parent)
    end

    # ========================================================================
    # 8. Update active cell mask
    # ========================================================================
    mesh.active_cell_mask[cell_id] = false
    push!(mesh.active_cell_mask, true)  # child_a
    push!(mesh.active_cell_mask, true)  # child_b

    mesh.n_active_cells += 1  # Net gain of 1 (removed 1, added 2)
    mesh.max_level = max(mesh.max_level, new_level)

    assert_amr_invariants(mesh)
    return (child_a_id, child_b_id)
end

"""
	estimate_gradient(mesh::Mesh1D, cell_id::Int, φ::Float64, var) -> Float64

Estimate gradient of a field at a cell for interpolation.
"""
function estimate_gradient(mesh::Mesh1D{NC}, cell_id::Int, φ::Float64, var) where {NC}
    cell = mesh.cells[cell_id]
    left_face = mesh.faces[cell.face_ids[1]]
    right_face = mesh.faces[cell.face_ids[2]]

    left_id = left_face.cell_ids[1]
    right_id = right_face.cell_ids[2]

    # Get neighbor values
    if left_id == 0 && right_id == 0
        return 0.0
    elseif left_id == 0
        φ_R = var == :T ? mesh.cell_states[right_id].T : mesh.cell_states[right_id].ξ[var]
        dz = mesh.cell_centers[right_id] - mesh.cell_centers[cell_id]
        return (φ_R - φ) / dz
    elseif right_id == 0
        φ_L = var == :T ? mesh.cell_states[left_id].T : mesh.cell_states[left_id].ξ[var]
        dz = mesh.cell_centers[cell_id] - mesh.cell_centers[left_id]
        return (φ - φ_L) / dz
    else
        φ_L = var == :T ? mesh.cell_states[left_id].T : mesh.cell_states[left_id].ξ[var]
        φ_R = var == :T ? mesh.cell_states[right_id].T : mesh.cell_states[right_id].ξ[var]
        # Opposite-distance-weighted slopes: second order on nonuniform
        # meshes
        d_L = mesh.cell_centers[cell_id] - mesh.cell_centers[left_id]
        d_R = mesh.cell_centers[right_id] - mesh.cell_centers[cell_id]
        return opposite_distance_gradient(φ_L, φ, φ_R, d_L, d_R)
    end
end

# ============================================================================
# Batch Refinement
# ============================================================================

"""
	mark_cells_for_refinement(mesh::Mesh1D, η::Vector{Float64}, 
							  θ_refine::Float64, h_min::Float64, max_level::Int)
							  -> Vector{Int}

Mark cells for refinement based on error indicators.

# Arguments
- `mesh`: The mesh
- `η`: Error indicator values
- `θ_refine`: Refinement threshold (refine if η > θ_refine * max(η))
- `h_min`: Minimum allowed cell size
- `max_level`: Maximum refinement level

# Returns
- Vector of cell indices to refine
"""
function mark_cells_for_refinement(
    mesh::Mesh1D,
    η::Vector{Float64},
    θ_refine::Float64,
    h_min::Float64,
    max_level::Int,
)
    marked = Int[]
    η_max = maximum(η[i] for i in eachindex(η) if mesh.active_cell_mask[i])

    if η_max ≤ 0
        return marked
    end

    threshold = θ_refine * η_max

    @inbounds for i in eachindex(η)
        if !mesh.active_cell_mask[i]
            continue
        end

        cell = mesh.cells[i]
        h = cell_size(mesh, i)

        # Check refinement criteria
        if η[i] > threshold && h > 2 * h_min && cell.level < max_level
            push!(marked, i)
        end
    end

    return marked
end

"""
	enforce_2to1_balance!(mesh::Mesh1D, marked::Vector{Int}, max_level::Int) -> Vector{Int}

Ensure 2:1 refinement level balance by marking additional cells.

Adjacent cells should differ by at most 1 refinement level.

# Returns
- Expanded vector of cells to refine (including balance cells)
"""
function enforce_2to1_balance!(mesh::Mesh1D, marked::Vector{Int}, max_level::Int)
    marked_set = Set(marked)
    changed = true

    while changed
        changed = false
        new_marks = Int[]

        for cell_id in marked_set
            cell = mesh.cells[cell_id]
            new_level = cell.level + 1

            # Check left neighbor
            left_face = mesh.faces[cell.face_ids[1]]
            left_id = left_face.cell_ids[1]
            if left_id > 0 && mesh.active_cell_mask[left_id]
                left_cell = mesh.cells[left_id]
                # If left neighbor would be 2+ levels coarser after refinement
                if left_cell.level < cell.level && left_cell.level < max_level
                    if !(left_id in marked_set)
                        push!(new_marks, left_id)
                        changed = true
                    end
                end
            end

            # Check right neighbor
            right_face = mesh.faces[cell.face_ids[2]]
            right_id = right_face.cell_ids[2]
            if right_id > 0 && mesh.active_cell_mask[right_id]
                right_cell = mesh.cells[right_id]
                if right_cell.level < cell.level && right_cell.level < max_level
                    if !(right_id in marked_set)
                        push!(new_marks, right_id)
                        changed = true
                    end
                end
            end
        end

        for id in new_marks
            push!(marked_set, id)
        end
    end

    return collect(marked_set)
end

"""
	perform_h_refinement!(mesh::Mesh1D{NC}, marked::Vector{Int}) where NC

Perform h-refinement on all marked cells.

# Arguments
- `mesh`: The mesh to modify
- `marked`: Vector of cell indices to refine

# Notes
- Cells are refined in order of decreasing index to avoid invalidating indices
- After refinement, call `update_geometry!(mesh)` to refresh caches
"""
function perform_h_refinement!(mesh::Mesh1D{NC}, marked::Vector{Int}) where {NC}
    if isempty(marked)
        return nothing
    end

    # Sort in decreasing order to avoid index invalidation issues
    # (actually, we append new cells, so order doesn't matter for indices)
    # But we still sort for deterministic behavior
    sorted_marks = sort(marked, rev = true)

    for cell_id in sorted_marks
        if mesh.active_cell_mask[cell_id]  # Double-check still active
            refine_cell!(mesh, cell_id)
        end
    end

    # Update geometry after all refinements
    update_geometry!(mesh)

    return nothing
end

# ============================================================================
# Conservation Verification
# ============================================================================

"""
	verify_refinement_conservation(mesh::Mesh1D{NC}, cell_id::Int) -> NamedTuple

Verify that refinement conserved mass and energy.

# Returns
- NamedTuple with fields: mass_error, energy_error
"""
function verify_refinement_conservation(mesh::Mesh1D{NC}, cell_id::Int) where {NC}
    cell = mesh.cells[cell_id]

    if cell.children == (0, 0)
        return (mass_error = 0.0, energy_error = 0.0)
    end

    child_a_id, child_b_id = cell.children
    parent_state = mesh.cell_states[cell_id]
    child_a_state = mesh.cell_states[child_a_id]
    child_b_state = mesh.cell_states[child_b_id]

    V_p = parent_state.volume_ref
    V_a = child_a_state.volume
    V_b = child_b_state.volume

    # Mass conservation for each component
    max_mass_error = 0.0
    for j = 1:NC
        mass_parent = V_p * parent_state.ξ[j]
        mass_children = V_a * child_a_state.ξ[j] + V_b * child_b_state.ξ[j]
        error = abs(mass_children - mass_parent) / max(mass_parent, 1e-10)
        max_mass_error = max(max_mass_error, error)
    end

    # Energy (thermal) conservation
    # Using ρcp*T*V as proxy (assuming constant properties)
    E_parent = V_p * parent_state.T
    E_children = V_a * child_a_state.T + V_b * child_b_state.T
    energy_error = abs(E_children - E_parent) / max(E_parent, 1e-10)

    return (mass_error = max_mass_error, energy_error = energy_error)
end
