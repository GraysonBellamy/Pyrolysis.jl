"""
Cell merging for surface depletion (adjacent-cell merge).

Energy-conserving temperature merge plus the adjacent-cell / surface-cell merge
and mesh-compaction routines used by the core surface-depletion solve callback
(`Solver/callbacks.jl`). The AMR sibling-coarsening path lives in
`Experimental/h_coarsening.jl`.
"""

using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, update_geometry!, n_cells
using ..Materials: AbstractMaterial, Material, SOLID, T_REF, is_gas

# ============================================================================
# Energy-Conserving Temperature Merge
# ============================================================================

# 4-point Gauss–Legendre nodes/weights on [-1, 1]. Exact for polynomial cₚ(T)
# up to degree 7; merges are rare callback events, so the quadrature cost is
# irrelevant.
const _GL4_X =
    (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const _GL4_W =
    (0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385)

# ∫_{T_lo}^{T_hi} cₚ(T′) dT′ by 4-point Gauss–Legendre quadrature.
@inline function _integrate_cp(c_p, T_lo::Float64, T_hi::Float64)
    half = 0.5 * (T_hi - T_lo)
    mid = 0.5 * (T_hi + T_lo)
    s = 0.0
    @inbounds for k = 1:4
        s += _GL4_W[k] * Float64(c_p(mid + half * _GL4_X[k]))
    end
    return half * s
end

"""
	merge_temperatures_conservatively(V_a::Float64, state_a::CellState{NC},
									   V_b::Float64, state_b::CellState{NC},
									   material::Material{NC}) where NC -> Float64

Compute the merged temperature by conserving the condensed-phase sensible
enthalpy.

# Principle
The conserved functional is

```
H(T) = Σ_{j ∉ gas} m_j ∫_{T_REF}^{T} cₚ,j(T′) dT′ ,   m_j = V·ξ_j
```

and `T_merged` solves `H_merged(T_merged) = H_a(T_a) + H_b(T_b)`. Two choices
here are deliberate:

- **∫cₚ dT, not cₚ(T)·T.** For temperature-dependent cₚ the product form is
  not an enthalpy: with cₚ = a + bT it carries the b-term with coefficient 1
  instead of ½, a few-K systematic bias per merge when a hot surface cell
  merges into a cooler interior one. The integral form is exact for any cₚ(T)
  (4-point Gauss–Legendre per component). The datum `T_REF` cancels because
  each component's total mass is identical before and after the merge.
- **Gas components are excluded.** Formulation A stores thermal energy in the
  condensed matrix only (the gas is thermally quasi-steady), so gas heat
  capacity does not belong in the merge weighting.

The solve is a bracketed Newton iteration on `[min(T_a,T_b), max(T_a,T_b)]`;
`H` is strictly monotone (cₚ > 0) so the root is unique and always inside the
bracket.

# Arguments
- `V_a`: Volume of cell a [m³]
- `state_a`: Cell state for cell a (T, ξ)
- `V_b`: Volume of cell b [m³]
- `state_b`: Cell state for cell b (T, ξ)
- `material`: Material definition with component heat capacity functions

# Returns
Merged temperature [K] conserving condensed-phase sensible enthalpy.

# Edge case
If neither cell carries condensed heat capacity (both essentially pure gas),
the condensed energy balance is degenerate and the volume-weighted average is
returned.
"""
function merge_temperatures_conservatively(
    V_a::Float64,
    state_a::CellState{NC},
    V_b::Float64,
    state_b::CellState{NC},
    material::Material{NC},
) where {NC}
    T_a = state_a.T
    T_b = state_b.T
    V_total = V_a + V_b

    # Condensed-phase volumetric heat capacities (degenerate-case test and
    # Newton initial guess) and total condensed sensible enthalpy.
    ρcₚ_a = 0.0
    ρcₚ_b = 0.0
    H_total = 0.0
    @inbounds for j = 1:NC
        comp = material.components[j]
        is_gas(comp) && continue
        c_j = comp.heat_capacity
        ρcₚ_a += state_a.ξ[j] * Float64(c_j(T_a))
        ρcₚ_b += state_b.ξ[j] * Float64(c_j(T_b))
        H_total +=
            V_a * state_a.ξ[j] * _integrate_cp(c_j, T_REF, T_a) +
            V_b * state_b.ξ[j] * _integrate_cp(c_j, T_REF, T_b)
    end

    if ρcₚ_a < 1e-20 && ρcₚ_b < 1e-20
        # No condensed matter in either cell: the condensed energy balance
        # cannot determine T. Fall back to volume-weighted averaging.
        return (V_a * T_a + V_b * T_b) / V_total
    end

    T_lo, T_hi = minmax(T_a, T_b)
    T_hi - T_lo < 1e-12 && return T_a

    # Merged condensed enthalpy functional and its derivative. Masses are
    # merged exactly (ξ is volume-weighted), so m_j = V_a ξ_a,j + V_b ξ_b,j.
    H_merged = T -> begin
        s = 0.0
        @inbounds for j = 1:NC
            comp = material.components[j]
            is_gas(comp) && continue
            m_j = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
            s += m_j * _integrate_cp(comp.heat_capacity, T_REF, T)
        end
        s
    end
    Hprime_merged = T -> begin
        s = 0.0
        @inbounds for j = 1:NC
            comp = material.components[j]
            is_gas(comp) && continue
            m_j = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
            s += m_j * Float64(comp.heat_capacity(T))
        end
        s
    end

    # Bracketed Newton: H is strictly increasing, and the root lies in
    # [T_lo, T_hi] because H_merged(T_lo) ≤ H_total ≤ H_merged(T_hi).
    T = (V_a * ρcₚ_a * T_a + V_b * ρcₚ_b * T_b) / (V_a * ρcₚ_a + V_b * ρcₚ_b)
    T = clamp(T, T_lo, T_hi)
    lo, hi = T_lo, T_hi
    for _ = 1:60
        f = H_merged(T) - H_total
        if f > 0.0
            hi = T
        else
            lo = T
        end
        fp = Hprime_merged(T)
        T_new = fp > 1e-20 ? T - f / fp : 0.5 * (lo + hi)
        if !(lo < T_new < hi)
            T_new = 0.5 * (lo + hi)  # Newton left the bracket: bisect
        end
        if abs(T_new - T) < 1e-9
            return T_new
        end
        T = T_new
    end
    return T
end

"""
	merge_temperatures_conservatively_from_values(V_a::Float64, T_a::Float64, ξ_a::NTuple{NC,Float64},
												   V_b::Float64, T_b::Float64, ξ_b::NTuple{NC,Float64},
												   material::Material{NC}) where NC -> Float64

Alternate interface for energy-conserving temperature merge using raw values instead of CellState.
Useful when working directly with state vectors.
"""
function merge_temperatures_conservatively_from_values(
    V_a::Float64,
    T_a::Float64,
    ξ_a::NTuple{NC,Float64},
    V_b::Float64,
    T_b::Float64,
    ξ_b::NTuple{NC,Float64},
    material::Material{NC},
) where {NC}

    # Create temporary cell states
    state_a = CellState{NC}(T_a, ξ_a, V_a, V_a)
    state_b = CellState{NC}(T_b, ξ_b, V_b, V_b)

    return merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
end

# ============================================================================
# Adjacent Cell Merging (Extended h-coarsening for surface depletion)
# ============================================================================

"""
	find_shared_face(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int) -> Int

Find the face shared between two adjacent cells.

# Returns
- Face ID of shared face, or 0 if cells are not adjacent
"""
function find_shared_face(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int)
    cell_a = mesh.cells[cell_a_id]
    cell_b = mesh.cells[cell_b_id]

    # Check each face of cell_a
    for face_id in cell_a.face_ids
        face = mesh.faces[face_id]
        # Check if this face connects to cell_b
        if face.cell_ids[1] == cell_b_id || face.cell_ids[2] == cell_b_id
            return face_id
        end
    end

    return 0  # Not adjacent
end

"""
	are_cells_adjacent(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int) -> Bool

Check if two cells share a face.
"""
function are_cells_adjacent(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int)
    return find_shared_face(mesh, cell_a_id, cell_b_id) != 0
end

"""
	merge_adjacent_cells!(mesh::Mesh1D{NC}, cell_to_remove_id::Int, 
						  cell_to_keep_id::Int; material=nothing) where NC -> Int

Merge two adjacent cells, removing one and expanding the other.

This is an extended h-coarsening that works on ANY adjacent pair,
not just sibling cells from refinement. The key difference from
`coarsen_cells!` is that we don't require a parent relationship.

# Arguments
- `mesh`: The mesh to modify
- `cell_to_remove_id`: Cell to remove (typically depleted surface cell)
- `cell_to_keep_id`: Cell to expand (will absorb the removed cell's volume/mass)
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- ID of the merged cell (same as cell_to_keep_id)

# Algorithm
1. Verify cells are adjacent
2. Conservative merge: combine states with mass/energy conservation
3. Expand the kept cell's geometry to encompass both
4. Update face connectivity
5. Mark removed cell as inactive

# Physical Use Case
During pyrolysis, surface cells become depleted as solid decomposes.
When the surface cell is mostly gas, we merge it with the interior
neighbor to remove the singularity in dT/dt = Q/(ρcₚ).

# Conservation
- Mass: ξ_merged = (V_a * ξ_a + V_b * ξ_b) / V_merged
- Energy: When `material` is provided, uses ρcₚ-weighted energy conservation.
  Otherwise, uses simple volume-weighted temperature averaging.

# Note
Unlike `coarsen_cells!`, this function:
- Does NOT require a parent-child relationship
- Can work on any adjacent cell pair
- Is designed for removing depleted cells, not low-error coarsening
"""
function merge_adjacent_cells!(
    mesh::Mesh1D{NC},
    cell_to_remove_id::Int,
    cell_to_keep_id::Int;
    material::Union{Nothing,Material{NC}} = nothing,
) where {NC}

    # ========================================================================
    # 1. Validation
    # ========================================================================
    if !mesh.active_cell_mask[cell_to_remove_id]
        error("Cell $cell_to_remove_id is not active")
    end
    if !mesh.active_cell_mask[cell_to_keep_id]
        error("Cell $cell_to_keep_id is not active")
    end

    shared_face_id = find_shared_face(mesh, cell_to_remove_id, cell_to_keep_id)
    if shared_face_id == 0
        error("Cells $cell_to_remove_id and $cell_to_keep_id are not adjacent")
    end

    cell_remove = mesh.cells[cell_to_remove_id]
    cell_keep = mesh.cells[cell_to_keep_id]
    shared_face = mesh.faces[shared_face_id]

    # ========================================================================
    # 2. Determine spatial arrangement
    # ========================================================================
    # In 1D: left = lower z, right = higher z
    z_remove = mesh.cell_centers[cell_to_remove_id]
    z_keep = mesh.cell_centers[cell_to_keep_id]

    if z_remove > z_keep
        # Removing the right cell (surface cell)
        # Keep cell will expand to the right
        is_removing_right = true
    else
        # Removing the left cell (unusual but handle it)
        # Keep cell will expand to the left
        is_removing_right = false
    end

    # ========================================================================
    # 3. Conservative state merging
    # ========================================================================
    state_remove = mesh.cell_states[cell_to_remove_id]
    state_keep = mesh.cell_states[cell_to_keep_id]

    V_remove = mesh.cell_volumes[cell_to_remove_id]
    V_keep = mesh.cell_volumes[cell_to_keep_id]
    V_merged = V_remove + V_keep

    # Temperature: use energy-conserving merge if material is provided
    if material !== nothing
        T_merged = merge_temperatures_conservatively(
            V_remove,
            state_remove,
            V_keep,
            state_keep,
            material,
        )
    else
        # Fallback: volume-weighted average (less accurate but works without material)
        T_merged = (V_remove * state_remove.T + V_keep * state_keep.T) / V_merged
    end

    # Concentrations: conservative (mass-preserving)
    ξ_merged = ntuple(NC) do j
        (V_remove * state_remove.ξ[j] + V_keep * state_keep.ξ[j]) / V_merged
    end

    # Update the kept cell's state
    mesh.cell_states[cell_to_keep_id] = CellState{NC}(
        T_merged,
        ξ_merged,
        V_merged,
        mesh.cell_states[cell_to_keep_id].volume_ref +
        mesh.cell_states[cell_to_remove_id].volume_ref,
    )
    mesh.cell_volumes[cell_to_keep_id] = V_merged

    # Preserve initial dry-solid mass extensively across the merge.
    if !isempty(mesh.m_dry_solid_initial)
        mesh.m_dry_solid_initial[cell_to_keep_id] +=
            mesh.m_dry_solid_initial[cell_to_remove_id]
    end

    # ========================================================================
    # 4. Update node positions (expand kept cell)
    # ========================================================================
    if is_removing_right
        # The kept cell expands to take over the right node of removed cell
        new_right_node_id = cell_remove.node_ids[2]  # Right node of removed cell
        left_node_id = cell_keep.node_ids[1]  # Left node stays the same

        # Debug: check for node corruption
        if left_node_id == new_right_node_id
            @error "Node corruption in merge! is_removing_right=true" cell_to_remove_id cell_to_keep_id cell_remove.node_ids cell_keep.node_ids left_node_id new_right_node_id
        end

        # Update kept cell's node_ids
        mesh.cells[cell_to_keep_id] = Cell1D(
            cell_keep.id,
            (left_node_id, new_right_node_id),
            cell_keep.face_ids,  # Will update faces below
            cell_keep.level,
            cell_keep.parent,
            cell_keep.children,
        )
    else
        # The kept cell expands to take over the left node of removed cell
        new_left_node_id = cell_remove.node_ids[1]  # Left node of removed cell
        right_node_id = cell_keep.node_ids[2]  # Right node stays the same

        # Debug: check for node corruption
        if new_left_node_id == right_node_id
            @error "Node corruption in merge! is_removing_right=false" cell_to_remove_id cell_to_keep_id cell_remove.node_ids cell_keep.node_ids new_left_node_id right_node_id
        end

        # Update kept cell's node_ids
        mesh.cells[cell_to_keep_id] = Cell1D(
            cell_keep.id,
            (new_left_node_id, right_node_id),
            cell_keep.face_ids,  # Will update faces below
            cell_keep.level,
            cell_keep.parent,
            cell_keep.children,
        )
    end

    # ========================================================================
    # 5. Update face connectivity
    # ========================================================================
    if is_removing_right
        # Find the boundary face of the removed cell (its right face)
        # This becomes the new right face of the kept cell
        right_face_id_remove = cell_remove.face_ids[2]
        right_face = mesh.faces[right_face_id_remove]

        # Update the boundary face to point to the kept cell
        mesh.faces[right_face_id_remove] = Face1D(
            right_face.id,
            right_face.node_id,
            (cell_to_keep_id, right_face.cell_ids[2]),  # cell_ids[2] is boundary (0)
            right_face.is_boundary,
            right_face.boundary_id,
        )

        # Update kept cell's face_ids
        mesh.cells[cell_to_keep_id] = Cell1D(
            mesh.cells[cell_to_keep_id].id,
            mesh.cells[cell_to_keep_id].node_ids,
            (cell_keep.face_ids[1], right_face_id_remove),  # New right face
            mesh.cells[cell_to_keep_id].level,
            mesh.cells[cell_to_keep_id].parent,
            mesh.cells[cell_to_keep_id].children,
        )
    else
        # Find the boundary face of the removed cell (its left face)
        left_face_id_remove = cell_remove.face_ids[1]
        left_face = mesh.faces[left_face_id_remove]

        # Update the boundary face to point to the kept cell
        mesh.faces[left_face_id_remove] = Face1D(
            left_face.id,
            left_face.node_id,
            (left_face.cell_ids[1], cell_to_keep_id),  # cell_ids[1] is boundary (0)
            left_face.is_boundary,
            left_face.boundary_id,
        )

        # Update kept cell's face_ids
        mesh.cells[cell_to_keep_id] = Cell1D(
            mesh.cells[cell_to_keep_id].id,
            mesh.cells[cell_to_keep_id].node_ids,
            (left_face_id_remove, cell_keep.face_ids[2]),  # New left face
            mesh.cells[cell_to_keep_id].level,
            mesh.cells[cell_to_keep_id].parent,
            mesh.cells[cell_to_keep_id].children,
        )
    end

    # ========================================================================
    # 6. Boundary face dictionaries need no update: face connectivity was
    #    rewritten in step 5 and the boundary_faces dict keys/ids are stable.
    # ========================================================================

    # ========================================================================
    # 7. Mark removed cell as inactive
    # ========================================================================
    mesh.active_cell_mask[cell_to_remove_id] = false
    mesh.n_active_cells -= 1

    # ========================================================================
    # 8. Update geometry cache
    # ========================================================================
    # Update cell center of kept cell
    node_left = mesh.nodes[mesh.cells[cell_to_keep_id].node_ids[1]]
    node_right = mesh.nodes[mesh.cells[cell_to_keep_id].node_ids[2]]
    mesh.cell_centers[cell_to_keep_id] = (node_left.z + node_right.z) / 2

    assert_amr_invariants(mesh)
    return cell_to_keep_id
end

"""
	merge_surface_cell!(mesh::Mesh1D{NC}; material=nothing) where NC -> Int

Convenience function to merge the depleted surface cell with its interior neighbor.

# Arguments
- `mesh`: The mesh to modify
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- ID of the merged cell, or 0 if no merge was possible

# Usage
```julia
using Pyrolysis.Adaptivity

# When surface depletion is detected
if is_surface_depleted(mesh, material, 0.01)
	merged_id = merge_surface_cell!(mesh; material=material)
	if merged_id > 0
		# Update state vector etc.
	end
end
```
"""
function merge_surface_cell!(
    mesh::Mesh1D{NC};
    material::Union{Nothing,Material{NC}} = nothing,
) where {NC}
    # Import from surface_depletion.jl (same module)
    surface_id = find_surface_cell(mesh)
    if surface_id == 0
        return 0
    end

    neighbor_id = find_interior_neighbor(mesh, surface_id)
    if neighbor_id == 0
        return 0  # Surface cell has no interior neighbor (only 1 cell left)
    end

    return merge_adjacent_cells!(mesh, surface_id, neighbor_id; material = material)
end

"""
	compact_mesh_after_surface_merge!(mesh::Mesh1D{NC}) where NC

Remove the inactive surface cell and orphaned surface node/face from mesh arrays.

This function should be called AFTER `merge_surface_cell!` to actually remove
the inactive cell from the arrays. This is necessary because:
1. `merge_surface_cell!` only marks the surface cell inactive
2. The ODE state vector actually removes elements  
3. The mesh arrays must match the state vector length

# Algorithm
After surface merge, the mesh has:
- Surface cell marked inactive (last active becomes new surface)
- Surface node orphaned (not connected to any active cell)
- Interior face between merged cells orphaned

This function:
1. Identifies orphaned elements (cell, node, face)
2. Removes them from arrays
3. Renumbers all indices consistently
4. Updates connectivity
5. Updates geometry caches

# Returns
Nothing. Modifies mesh in place.

# Note
This is specifically for surface cell removal where the surface cell is always
the highest-indexed active cell. For general compaction with arbitrary inactive
cells, a more complex algorithm would be needed.
"""
function compact_mesh_after_surface_merge!(mesh::Mesh1D{NC}) where {NC}
    n_cells_old = length(mesh.cells)
    n_nodes_old = length(mesh.nodes)
    n_faces_old = length(mesh.faces)

    if n_cells_old < 2
        error("Cannot compact mesh with fewer than 2 cells")
    end

    # The surface cell (last one) should be inactive
    if mesh.active_cell_mask[n_cells_old]
        error("Last cell is active - did you call merge_surface_cell! first?")
    end

    # ========================================================================
    # Step 0: Identify what needs to be removed
    # ========================================================================
    removed_cell_id = n_cells_old
    removed_cell = mesh.cells[removed_cell_id]

    # The new surface cell is the one that absorbed the removed cell
    new_surface_cell_id = n_cells_old - 1

    # After merge_adjacent_cells!, the kept cell's node_ids span both old cells.
    # The orphaned node is the one between the two merged cells - this is the
    # node that belonged to both cells before merging, which is:
    #   - Right node of kept cell (new_surface_cell) BEFORE merge = left node of removed cell BEFORE merge
    #
    # From merge_adjacent_cells!, when merging surface cell (cell_remove) with its neighbor (cell_keep):
    #   cell_keep's new node_ids = (cell_keep.node_ids[1], cell_remove.node_ids[2])
    # So the orphaned node is cell_keep.node_ids[2] from BEFORE the merge, which equals
    # cell_remove.node_ids[1] from BEFORE the merge.
    #
    # BUT merge_adjacent_cells! already updated cell_keep's nodes, so we need to figure
    # out the orphaned node differently.

    # The orphaned node is the one NOT used by any active cell
    used_nodes = Set{Int}()
    for i = 1:n_cells_old
        if mesh.active_cell_mask[i]
            push!(used_nodes, mesh.cells[i].node_ids[1])
            push!(used_nodes, mesh.cells[i].node_ids[2])
        end
    end

    orphaned_node_id = 0
    for i = 1:n_nodes_old
        if !(i in used_nodes)
            orphaned_node_id = i
            break
        end
    end

    if orphaned_node_id == 0
        # This shouldn't happen - there should be exactly one orphaned node
        @warn "No orphaned node found after surface merge - this is unexpected"
        orphaned_node_id = n_nodes_old  # Fallback: assume surface node
    end

    # The orphaned face is the interior face that connected the two merged cells
    # This is the face that has the removed cell in its cell_ids
    orphaned_face_id = 0
    for i = 1:n_faces_old
        face = mesh.faces[i]
        if removed_cell_id in face.cell_ids && !face.is_boundary
            orphaned_face_id = i
            break
        end
    end

    # The boundary face needs to be kept but its cell_id updated
    # It's the right face of the removed surface cell
    boundary_face_id = removed_cell.face_ids[2]

    # ========================================================================
    # Step 1: Build index remapping
    # ========================================================================
    # We'll use a swap-and-pop approach: swap the removed element with the last,
    # then pop. This requires tracking how indices change.

    # Cell remapping: cell removed_cell_id is removed
    # Since removed_cell_id == n_cells_old, it's already last, just pop
    cell_remap = collect(1:n_cells_old)
    cell_remap[removed_cell_id] = 0  # Indicates removed

    # Node remapping
    node_remap = collect(1:n_nodes_old)
    if orphaned_node_id == n_nodes_old
        # Orphaned is last, just pop
        node_remap[orphaned_node_id] = 0
    else
        # Swap orphaned with last, then pop
        node_remap[orphaned_node_id] = 0  # Removed
        node_remap[n_nodes_old] = orphaned_node_id  # Last moves to orphaned's slot
    end

    # Face remapping
    face_remap = collect(1:n_faces_old)
    if orphaned_face_id == 0
        # No orphaned face found - this can happen if merge was at boundary
        @warn "No orphaned face found - skipping face removal"
    elseif orphaned_face_id == n_faces_old
        face_remap[orphaned_face_id] = 0
    else
        face_remap[orphaned_face_id] = 0
        face_remap[n_faces_old] = orphaned_face_id
    end

    # ========================================================================
    # Step 2: Update the boundary face to point to new surface cell
    # ========================================================================
    if boundary_face_id > 0 && boundary_face_id <= n_faces_old
        boundary_face = mesh.faces[boundary_face_id]
        mesh.faces[boundary_face_id] = Face1D(
            boundary_face.id,
            boundary_face.node_id,
            (new_surface_cell_id, 0),  # (left_cell=new_surface, right_cell=boundary)
            boundary_face.is_boundary,
            boundary_face.boundary_id,
        )
    end

    # ========================================================================
    # Step 3: Remove the cell (last cell)
    # ========================================================================
    pop!(mesh.cells)
    pop!(mesh.cell_states)
    pop!(mesh.cell_centers)
    pop!(mesh.cell_volumes)
    pop!(mesh.active_cell_mask)
    # χ denominator: parallel per-cell array, the removed (last) cell's slot is
    # stale after the merge summed its mass into the kept cell.
    if !isempty(mesh.m_dry_solid_initial)
        pop!(mesh.m_dry_solid_initial)
    end

    # ========================================================================
    # Step 4: Remove the orphaned face
    # ========================================================================
    if orphaned_face_id > 0
        if orphaned_face_id == n_faces_old
            pop!(mesh.faces)
            pop!(mesh.face_positions)
            pop!(mesh.face_areas)
        else
            # Swap with last then pop
            mesh.faces[orphaned_face_id] = mesh.faces[n_faces_old]
            mesh.face_positions[orphaned_face_id] = mesh.face_positions[n_faces_old]
            mesh.face_areas[orphaned_face_id] = mesh.face_areas[n_faces_old]

            # Update the moved face's ID
            moved_face = mesh.faces[orphaned_face_id]
            mesh.faces[orphaned_face_id] = Face1D(
                orphaned_face_id,
                moved_face.node_id,
                moved_face.cell_ids,
                moved_face.is_boundary,
                moved_face.boundary_id,
            )

            pop!(mesh.faces)
            pop!(mesh.face_positions)
            pop!(mesh.face_areas)
        end

        n_faces_new = length(mesh.faces)

        # Update cell face references
        for (i, cell) in enumerate(mesh.cells)
            new_face_ids = Tuple(
                f == n_faces_old ? orphaned_face_id : (f > n_faces_new ? f : f) for
                f in cell.face_ids
            )
            # Simpler: just remap faces
            new_face_ids = Tuple(face_remap[f] for f in cell.face_ids)
            # Make sure no face is 0 (removed)
            if any(f == 0 for f in new_face_ids)
                # The cell references the orphaned face - this shouldn't happen
                # after merge_adjacent_cells! properly updated things
                @warn "Cell $i references orphaned face"
            end
            if new_face_ids != cell.face_ids
                mesh.cells[i] = Cell1D(
                    cell.id,
                    cell.node_ids,
                    new_face_ids,
                    cell.level,
                    cell.parent,
                    cell.children,
                )
            end
        end

        # Update boundary_faces dict
        for (bid, fids) in mesh.boundary_faces
            mesh.boundary_faces[bid] = [face_remap[f] for f in fids if face_remap[f] != 0]
        end
    end

    # ========================================================================
    # Step 5: Remove the orphaned node
    # ========================================================================
    n_nodes_new = n_nodes_old - 1

    if orphaned_node_id == n_nodes_old
        # Orphaned node is last - just pop
        pop!(mesh.nodes)
        pop!(mesh.node_velocities)
    else
        # Swap with last then pop
        mesh.nodes[orphaned_node_id] = mesh.nodes[n_nodes_old]
        mesh.node_velocities[orphaned_node_id] = mesh.node_velocities[n_nodes_old]

        # Update the moved node's ID
        moved_node = mesh.nodes[orphaned_node_id]
        mesh.nodes[orphaned_node_id] = Node1D(
            orphaned_node_id,
            moved_node.z,
            moved_node.z_ref,
            moved_node.is_boundary,
            moved_node.boundary_id,
        )

        pop!(mesh.nodes)
        pop!(mesh.node_velocities)

        # Update cell node references (node n_nodes_old is now at orphaned_node_id)
        for (i, cell) in enumerate(mesh.cells)
            new_node_ids =
                Tuple(n == n_nodes_old ? orphaned_node_id : n for n in cell.node_ids)
            if new_node_ids != cell.node_ids
                mesh.cells[i] = Cell1D(
                    cell.id,
                    new_node_ids,
                    cell.face_ids,
                    cell.level,
                    cell.parent,
                    cell.children,
                )
            end
        end

        # Update face node references
        for (i, face) in enumerate(mesh.faces)
            if face.node_id == n_nodes_old
                mesh.faces[i] = Face1D(
                    face.id,
                    orphaned_node_id,
                    face.cell_ids,
                    face.is_boundary,
                    face.boundary_id,
                )
            end
        end
    end

    # ========================================================================
    # Step 6: Update counts and verify
    # ========================================================================
    mesh.n_active_cells = sum(mesh.active_cell_mask)

    # ========================================================================
    # Step 7: Update geometry cache from current node positions
    # ========================================================================
    update_geometry!(mesh)

    # ========================================================================
    # Step 8: Verify consistency
    # ========================================================================
    n_cells_final = length(mesh.cells)
    n_nodes_final = length(mesh.nodes)
    n_faces_final = length(mesh.faces)

    @assert n_nodes_final == n_cells_final + 1 "Node count mismatch: expected $(n_cells_final + 1) nodes, got $n_nodes_final"
    @assert n_faces_final == n_cells_final + 1 "Face count mismatch: expected $(n_cells_final + 1) faces, got $n_faces_final"
    @assert all(mesh.active_cell_mask) "Inactive cells remain after compaction"

    return nothing
end

"""
	compact_mesh_after_merge!(mesh::Mesh1D{NC}, removed_cell_id::Int) where NC

Remove an inactive cell and its orphaned node/face from mesh arrays.
Works for any cell position (surface, interior, or bottom-adjacent).

This is the generalized version of `compact_mesh_after_surface_merge!`.

# Algorithm
After merge_adjacent_cells!, the mesh has:
- One cell marked inactive (the removed cell)
- One orphaned node (not connected to any active cell)
- One orphaned interior face (was between the merged cells)

This function:
1. Identifies orphaned elements
2. Removes elements using SHIFT (O(N) per element, maintains ordering)
3. Renumbers all indices consistently (simple decrement for elements after removed)
4. Updates connectivity and geometry

Note: We use shift instead of swap-and-pop to maintain the invariant that
cell/node/face indices correspond to spatial ordering (lower index = lower z).
This is important for 1D mesh operations that rely on left/right semantics.

# Arguments
- `mesh`: Mesh to compact (modified in place)
- `removed_cell_id`: ID of the cell that was merged away (should be inactive)

# Returns
Nothing. Modifies mesh in place.
"""
function compact_mesh_after_merge!(mesh::Mesh1D{NC}, removed_cell_id::Int) where {NC}
    n_cells_old = length(mesh.cells)
    n_nodes_old = length(mesh.nodes)
    n_faces_old = length(mesh.faces)

    if n_cells_old < 2
        error("Cannot compact mesh with fewer than 2 cells")
    end

    if removed_cell_id < 1 || removed_cell_id > n_cells_old
        error("Invalid removed_cell_id: $removed_cell_id (mesh has $n_cells_old cells)")
    end

    if mesh.active_cell_mask[removed_cell_id]
        error(
            "Cell $removed_cell_id is still active - did you call merge_adjacent_cells! first?",
        )
    end

    # ========================================================================
    # Step 1: Identify orphaned node (not used by any active cell)
    # ========================================================================
    used_nodes = Set{Int}()
    for i = 1:n_cells_old
        if mesh.active_cell_mask[i]
            push!(used_nodes, mesh.cells[i].node_ids[1])
            push!(used_nodes, mesh.cells[i].node_ids[2])
        end
    end

    orphaned_node_id = 0
    for i = 1:n_nodes_old
        if !(i in used_nodes)
            orphaned_node_id = i
            break
        end
    end

    if orphaned_node_id == 0
        @warn "No orphaned node found after merge - this may indicate mesh corruption" n_cells_old n_nodes_old
    end

    # ========================================================================
    # Step 2: Identify orphaned interior face (connected to removed cell, not boundary)
    # ========================================================================
    orphaned_face_id = 0
    for i = 1:n_faces_old
        face = mesh.faces[i]
        if removed_cell_id in face.cell_ids && !face.is_boundary
            orphaned_face_id = i
            break
        end
    end

    # ========================================================================
    # Step 3: Build index remappings for SHIFT-based removal
    # With shift: elements after removed index get decremented by 1
    # ========================================================================

    # Cell remapping: old_id -> new_id (0 = removed)
    cell_remap = collect(1:n_cells_old)
    cell_remap[removed_cell_id] = 0  # Removed cell
    for i = (removed_cell_id+1):n_cells_old
        cell_remap[i] = i - 1  # Shift down by 1
    end

    # Node remapping
    node_remap = collect(1:n_nodes_old)
    if orphaned_node_id > 0
        node_remap[orphaned_node_id] = 0  # Removed node
        for i = (orphaned_node_id+1):n_nodes_old
            node_remap[i] = i - 1  # Shift down by 1
        end
    end

    # Face remapping
    face_remap = collect(1:n_faces_old)
    if orphaned_face_id > 0
        face_remap[orphaned_face_id] = 0  # Removed face
        for i = (orphaned_face_id+1):n_faces_old
            face_remap[i] = i - 1  # Shift down by 1
        end
    end

    # ========================================================================
    # Step 4: Remove cell using deleteat! (shift)
    # ========================================================================
    deleteat!(mesh.cells, removed_cell_id)
    deleteat!(mesh.cell_states, removed_cell_id)
    deleteat!(mesh.cell_centers, removed_cell_id)
    deleteat!(mesh.cell_volumes, removed_cell_id)
    deleteat!(mesh.active_cell_mask, removed_cell_id)
    # The χ denominator is a parallel per-cell array: its slot for the removed
    # cell is stale (the merge already added that mass to the kept cell), and
    # every consumer indexes it by compacted cell id.
    if !isempty(mesh.m_dry_solid_initial)
        deleteat!(mesh.m_dry_solid_initial, removed_cell_id)
    end

    # ========================================================================
    # Step 5: Remove orphaned face using deleteat! (shift)
    # ========================================================================
    if orphaned_face_id > 0
        deleteat!(mesh.faces, orphaned_face_id)
        deleteat!(mesh.face_positions, orphaned_face_id)
        deleteat!(mesh.face_areas, orphaned_face_id)
    end

    # ========================================================================
    # Step 6: Remove orphaned node using deleteat! (shift)
    # ========================================================================
    if orphaned_node_id > 0
        deleteat!(mesh.nodes, orphaned_node_id)
        deleteat!(mesh.node_velocities, orphaned_node_id)
    end

    # ========================================================================
    # Step 7: Reindex all cell references (nodes, faces, id)
    # ========================================================================
    for i in eachindex(mesh.cells)
        cell = mesh.cells[i]
        # Apply remap: if remap[x] > 0, use remap[x]; if remap[x] == 0, ERROR
        new_node_ids = Tuple(
            begin
                mapped = node_remap[n]
                if mapped == 0
                    error("Cell $i references removed node $n")
                end
                mapped
            end for n in cell.node_ids
        )

        new_face_ids = Tuple(
            begin
                mapped = face_remap[f]
                if mapped == 0
                    error("Cell $i references removed face $f")
                end
                mapped
            end for f in cell.face_ids
        )

        new_parent = cell.parent > 0 ? cell_remap[cell.parent] : 0
        new_children = (
            cell.children[1] > 0 ? cell_remap[cell.children[1]] : 0,
            cell.children[2] > 0 ? cell_remap[cell.children[2]] : 0,
        )

        mesh.cells[i] =
            Cell1D(i, new_node_ids, new_face_ids, cell.level, new_parent, new_children)
    end

    # ========================================================================
    # Step 8: Reindex all face references (nodes, cells, id)
    # ========================================================================
    for i in eachindex(mesh.faces)
        face = mesh.faces[i]
        new_node_id = node_remap[face.node_id]
        if new_node_id == 0
            error("Face $i references removed node $(face.node_id)")
        end

        # Handle cell_ids: 0 = boundary (keep as 0), otherwise apply remap
        new_cell_1 = face.cell_ids[1] == 0 ? 0 : cell_remap[face.cell_ids[1]]
        new_cell_2 = face.cell_ids[2] == 0 ? 0 : cell_remap[face.cell_ids[2]]

        mesh.faces[i] = Face1D(
            i,
            new_node_id,
            (new_cell_1, new_cell_2),
            face.is_boundary,
            face.boundary_id,
        )
    end

    # ========================================================================
    # Step 9: Reindex node IDs
    # ========================================================================
    for i in eachindex(mesh.nodes)
        node = mesh.nodes[i]
        mesh.nodes[i] = Node1D(i, node.z, node.z_ref, node.is_boundary, node.boundary_id)
    end

    # ========================================================================
    # Step 10: Update boundary_faces dict
    # ========================================================================
    for (bid, fids) in mesh.boundary_faces
        mesh.boundary_faces[bid] = [face_remap[f] for f in fids if face_remap[f] != 0]
    end

    # ========================================================================
    # Step 11: Update counts and geometry
    # ========================================================================
    mesh.n_active_cells = sum(mesh.active_cell_mask)
    update_geometry!(mesh)

    # ========================================================================
    # Step 12: Verify consistency
    # ========================================================================
    n_cells_final = length(mesh.cells)
    n_nodes_final = length(mesh.nodes)
    n_faces_final = length(mesh.faces)

    @assert n_nodes_final == n_cells_final + 1 "Node count mismatch: expected $(n_cells_final + 1) nodes, got $n_nodes_final"
    @assert n_faces_final == n_cells_final + 1 "Face count mismatch: expected $(n_cells_final + 1) faces, got $n_faces_final"
    @assert all(mesh.active_cell_mask) "Inactive cells remain after compaction"

    return nothing
end
