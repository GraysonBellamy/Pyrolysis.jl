"""
h-coarsening (cell merging) — AMR sibling coarsening.

Merges sibling cells (created by a previous refinement) when error indicators
are low. This is the experimental h-AMR coarsening path. The adjacent-cell
merge machinery used by the core surface-depletion solve callback lives in
`Adaptivity/cell_merge.jl` and is imported here for the energy-conserving
temperature merge.
"""

using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, update_geometry!, n_cells
using ..Materials: AbstractMaterial, Material, SOLID
using ..Adaptivity: merge_temperatures_conservatively

# ============================================================================
# Single Cell Pair Coarsening
# ============================================================================

"""
	can_coarsen(mesh::Mesh1D, cell_id::Int) -> Bool

Check if a cell can be coarsened (merged with sibling).

A cell can be coarsened if:
1. It has a parent (is not a root cell)
2. It is a leaf (has no children)
3. Its sibling is also a leaf
"""
function can_coarsen(mesh::Mesh1D, cell_id::Int)
    cell = mesh.cells[cell_id]

    # Must have a parent
    if cell.parent == 0
        return false
    end

    # Must be a leaf
    if cell.children != (0, 0)
        return false
    end

    # Sibling must also be a leaf
    parent = mesh.cells[cell.parent]
    sibling_id = parent.children[1] == cell_id ? parent.children[2] : parent.children[1]
    sibling = mesh.cells[sibling_id]

    if sibling.children != (0, 0)
        return false
    end

    return true
end

"""
	get_sibling(mesh::Mesh1D, cell_id::Int) -> Int

Get the sibling cell ID (the other child of the same parent).
Returns 0 if no sibling exists.
"""
function get_sibling(mesh::Mesh1D, cell_id::Int)
    cell = mesh.cells[cell_id]

    if cell.parent == 0
        return 0
    end

    parent = mesh.cells[cell.parent]
    if parent.children[1] == cell_id
        return parent.children[2]
    else
        return parent.children[1]
    end
end

"""
	coarsen_cells!(mesh::Mesh1D{NC}, cell_a_id::Int, cell_b_id::Int; material=nothing) -> Int

Merge two sibling cells back into their parent.

# Arguments
- `mesh`: The mesh to modify
- `cell_a_id`: First child cell index
- `cell_b_id`: Second child cell index (sibling)
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- Parent cell index

# Algorithm
1. Verify cells are siblings and both are leaves
2. Compute conservative average of child states
3. Update parent state
4. Restore parent face connectivity
5. Mark children as inactive
6. Mark parent as active

# Energy Conservation Note
When `material` is provided, temperature is computed via energy conservation:
  E_total = V_a × ρcₚ_a × T_a + V_b × ρcₚ_b × T_b
  T_parent = E_total / (V_parent × ρcₚ_parent)

This is critical when merging cells with different compositions (e.g., solid vs gas),
as simple volume-weighted averaging can introduce significant temperature errors.

# Note
This does NOT remove the child cells or intermediate node from the mesh arrays.
They are simply marked inactive. This avoids index invalidation.
"""
function coarsen_cells!(
    mesh::Mesh1D{NC},
    cell_a_id::Int,
    cell_b_id::Int;
    material::Union{Nothing,Material{NC}} = nothing,
) where {NC}
    cell_a = mesh.cells[cell_a_id]
    cell_b = mesh.cells[cell_b_id]

    # Verify they are siblings
    if cell_a.parent != cell_b.parent || cell_a.parent == 0
        error("Cells $cell_a_id and $cell_b_id are not siblings")
    end

    # Verify both are leaves
    if cell_a.children != (0, 0) || cell_b.children != (0, 0)
        error("Cannot coarsen non-leaf cells")
    end

    parent_id = cell_a.parent
    parent = mesh.cells[parent_id]

    # Determine which is left (lower z) and right (higher z)
    if mesh.cell_centers[cell_a_id] < mesh.cell_centers[cell_b_id]
        left_child_id = cell_a_id
        right_child_id = cell_b_id
    else
        left_child_id = cell_b_id
        right_child_id = cell_a_id
    end

    left_child = mesh.cells[left_child_id]
    right_child = mesh.cells[right_child_id]

    # ========================================================================
    # 1. Get child states and volumes
    # ========================================================================
    state_a = mesh.cell_states[cell_a_id]
    state_b = mesh.cell_states[cell_b_id]

    V_a = mesh.cell_volumes[cell_a_id]
    V_b = mesh.cell_volumes[cell_b_id]
    V_parent = V_a + V_b

    # ========================================================================
    # 2. Compute conservative average
    # ========================================================================
    # Temperature: use energy-conserving merge if material is provided
    if material !== nothing
        T_parent = merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
    else
        # Fallback: volume-weighted average (less accurate but works without material)
        T_parent = (V_a * state_a.T + V_b * state_b.T) / V_parent
    end

    # Concentrations: conservative (mass-preserving)
    # ξ_parent = (V_a * ξ_a + V_b * ξ_b) / V_parent
    ξ_parent = ntuple(NC) do j
        (V_a * state_a.ξ[j] + V_b * state_b.ξ[j]) / V_parent
    end

    # ========================================================================
    # 3. Update parent state
    # ========================================================================
    mesh.cell_states[parent_id] =
        CellState{NC}(T_parent, ξ_parent, V_parent, mesh.cell_states[parent_id].volume_ref)
    mesh.cell_volumes[parent_id] = V_parent

    # Preserve initial dry-solid mass extensively across the coarsen.
    if !isempty(mesh.m_dry_solid_initial)
        mesh.m_dry_solid_initial[parent_id] =
            mesh.m_dry_solid_initial[cell_a_id] + mesh.m_dry_solid_initial[cell_b_id]
    end

    # ========================================================================
    # 4. Restore parent face connectivity
    # ========================================================================
    # The parent's faces should point to parent, not children

    # Left face of parent (was pointing to left child)
    left_face_id = parent.face_ids[1]
    left_face = mesh.faces[left_face_id]
    mesh.faces[left_face_id] = Face1D(
        left_face.id,
        left_face.node_id,
        (left_face.cell_ids[1], parent_id),  # Restore parent as right cell
        left_face.is_boundary,
        left_face.boundary_id,
    )

    # Right face of parent (was pointing to right child)
    right_face_id = parent.face_ids[2]
    right_face = mesh.faces[right_face_id]
    mesh.faces[right_face_id] = Face1D(
        right_face.id,
        right_face.node_id,
        (parent_id, right_face.cell_ids[2]),  # Restore parent as left cell
        right_face.is_boundary,
        right_face.boundary_id,
    )

    # ========================================================================
    # 5. Update parent cell to have no children
    # ========================================================================
    mesh.cells[parent_id] = Cell1D(
        parent.id,
        parent.node_ids,
        parent.face_ids,
        parent.level,
        parent.parent,
        (0, 0),  # No children anymore
    )

    # ========================================================================
    # 6. Update active cell mask
    # ========================================================================
    mesh.active_cell_mask[cell_a_id] = false
    mesh.active_cell_mask[cell_b_id] = false
    mesh.active_cell_mask[parent_id] = true

    mesh.n_active_cells -= 1  # Net loss of 1 (removed 2, added 1)

    # The intermediate face and node stay in the arrays, orphaned: compacting
    # them away would invalidate indices held throughout the code.

    assert_amr_invariants(mesh)
    return parent_id
end

# ============================================================================
# Batch Coarsening
# ============================================================================

"""
	mark_cells_for_coarsening(mesh::Mesh1D, η::Vector{Float64},
							  θ_coarsen::Float64, h_max::Float64)
							  -> Vector{Tuple{Int,Int}}

Mark sibling cell pairs for coarsening based on error indicators.

# Arguments
- `mesh`: The mesh
- `η`: Error indicator values
- `θ_coarsen`: Coarsening threshold (coarsen if η < θ_coarsen * max(η))
- `h_max`: Maximum allowed cell size

# Returns
- Vector of (cell_a, cell_b) tuples representing sibling pairs to coarsen
"""
function mark_cells_for_coarsening(
    mesh::Mesh1D,
    η::Vector{Float64},
    θ_coarsen::Float64,
    h_max::Float64,
)
    pairs = Tuple{Int,Int}[]
    processed = Set{Int}()

    η_max = maximum(η[i] for i in eachindex(η) if mesh.active_cell_mask[i]; init = 0.0)

    if η_max ≤ 0
        return pairs
    end

    threshold = θ_coarsen * η_max

    @inbounds for i in eachindex(η)
        if !mesh.active_cell_mask[i] || i in processed
            continue
        end

        if !can_coarsen(mesh, i)
            continue
        end

        cell = mesh.cells[i]
        sibling_id = get_sibling(mesh, i)

        if sibling_id == 0 || sibling_id in processed
            continue
        end

        # Check both cells have low error
        if η[i] < threshold && η[sibling_id] < threshold
            # Check resulting cell size
            h_combined = cell_size(mesh, i) + cell_size(mesh, sibling_id)
            if h_combined ≤ h_max
                push!(pairs, (i, sibling_id))
                push!(processed, i)
                push!(processed, sibling_id)
            end
        end
    end

    return pairs
end

"""
	perform_h_coarsening!(mesh::Mesh1D{NC}, pairs::Vector{Tuple{Int,Int}};
						  material=nothing) where NC

Perform h-coarsening on all marked sibling pairs.

# Arguments
- `mesh`: The mesh to modify
- `pairs`: Vector of (cell_a, cell_b) tuples to coarsen
- `material=nothing`: Material for the energy-conserving temperature merge.
  Without it every pair falls back to volume-weighted-T averaging, which does
  not conserve energy when the siblings have different compositions — pass the
  material whenever one is available.

# Notes
- After coarsening, call `update_geometry!(mesh)` to refresh caches
"""
function perform_h_coarsening!(
    mesh::Mesh1D{NC},
    pairs::Vector{Tuple{Int,Int}};
    material::Union{Nothing,Material{NC}} = nothing,
) where {NC}
    if isempty(pairs)
        return nothing
    end

    for (cell_a_id, cell_b_id) in pairs
        # Double-check both are still active (in case of complex scenarios)
        if mesh.active_cell_mask[cell_a_id] && mesh.active_cell_mask[cell_b_id]
            coarsen_cells!(mesh, cell_a_id, cell_b_id; material)
        end
    end

    # Update geometry after all coarsening
    update_geometry!(mesh)

    return nothing
end

# ============================================================================
# Conservation Verification  
# ============================================================================

"""
	verify_coarsening_conservation(V_a, state_a, V_b, state_b, V_p, state_p, NC) -> NamedTuple

Verify that coarsening conserved mass and energy.

# Arguments
- Pre-coarsening: V_a, state_a, V_b, state_b (children)
- Post-coarsening: V_p, state_p (parent)
- NC: Number of components

# Returns
- NamedTuple with fields: mass_error, energy_error
"""
function verify_coarsening_conservation(
    V_a::Float64,
    state_a::CellState{NC},
    V_b::Float64,
    state_b::CellState{NC},
    V_p::Float64,
    state_p::CellState{NC},
) where {NC}

    # Mass conservation for each component
    max_mass_error = 0.0
    for j = 1:NC
        mass_children = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
        mass_parent = V_p * state_p.ξ[j]
        error = abs(mass_parent - mass_children) / max(mass_children, 1e-10)
        max_mass_error = max(max_mass_error, error)
    end

    # Energy conservation (using T*V as proxy)
    E_children = V_a * state_a.T + V_b * state_b.T
    E_parent = V_p * state_p.T
    energy_error = abs(E_parent - E_children) / max(E_children, 1e-10)

    return (mass_error = max_mass_error, energy_error = energy_error)
end
