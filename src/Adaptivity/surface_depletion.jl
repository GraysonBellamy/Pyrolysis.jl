"""
Cell depletion utilities for mesh operations.

This module provides utility functions for finding cells and neighbors in
the mesh. These are used by the depletion callback in the Solver module.

## Usage

```julia
using Pyrolysis.Adaptivity

# Find surface and bottom cells
surface = find_surface_cell(mesh)
bottom = find_bottom_cell(mesh)

# Find neighbors
interior = find_interior_neighbor(mesh, cell_id)
exterior = find_exterior_neighbor(mesh, cell_id)
```
"""

using ..Geometry: Mesh1D, cell_size, cell_center, n_cells

# ============================================================================
# Surface Cell Detection
# ============================================================================

"""
	find_surface_cell(mesh::Mesh1D) -> Int

Find the ID of the current surface cell (highest z, active).

In the coordinate convention:
- z = 0 is the bottom (substrate)
- z = L is the top (exposed surface)

The surface cell is the active cell with the highest z center.
"""
function find_surface_cell(mesh::Mesh1D)
    max_z = -Inf
    surface_id = 0

    @inbounds for i in eachindex(mesh.cells)
        if mesh.active_cell_mask[i]
            z_center = mesh.cell_centers[i]
            if z_center > max_z
                max_z = z_center
                surface_id = i
            end
        end
    end

    return surface_id
end

"""
	find_bottom_cell(mesh::Mesh1D) -> Int

Find the ID of the bottom cell (lowest z, active).

This is the cell at z = 0 (substrate side).
"""
function find_bottom_cell(mesh::Mesh1D)
    min_z = Inf
    bottom_id = 0

    @inbounds for i in eachindex(mesh.cells)
        if mesh.active_cell_mask[i]
            z_center = mesh.cell_centers[i]
            if z_center < min_z
                min_z = z_center
                bottom_id = i
            end
        end
    end

    return bottom_id
end

"""
	find_interior_neighbor(mesh::Mesh1D, cell_id::Int) -> Int

Find the neighboring cell toward the interior (lower z).

For a surface cell at z = L, this returns the adjacent cell at lower z.
Returns 0 if no interior neighbor exists (cell is at bottom boundary).
"""
function find_interior_neighbor(mesh::Mesh1D, cell_id::Int)
    cell = mesh.cells[cell_id]

    # Left face (at lower z) connects to interior neighbor
    left_face_id = cell.face_ids[1]
    left_face = mesh.faces[left_face_id]

    # cell_ids[1] is the cell at lower z ("left" of face)
    neighbor_id = left_face.cell_ids[1]

    # neighbor_id == 0 means we're at the bottom boundary
    return neighbor_id
end

"""
	find_exterior_neighbor(mesh::Mesh1D, cell_id::Int) -> Int

Find the neighboring cell toward the exterior (higher z).

For an interior cell, this returns the adjacent cell at higher z.
Returns 0 if no exterior neighbor exists (cell is at top boundary).
"""
function find_exterior_neighbor(mesh::Mesh1D, cell_id::Int)
    cell = mesh.cells[cell_id]

    # Right face (at higher z) connects to exterior neighbor
    right_face_id = cell.face_ids[2]
    right_face = mesh.faces[right_face_id]

    # cell_ids[2] is the cell at higher z ("right" of face)
    neighbor_id = right_face.cell_ids[2]

    # neighbor_id == 0 means we're at the top boundary (surface)
    return neighbor_id
end

# ============================================================================
# Merge Neighbor Selection
# ============================================================================

"""
	find_merge_neighbor(mesh::Mesh1D, cell_id::Int) -> Int

Find the best neighbor to merge with for a thin cell.

Selection strategy:
1. If cell has only one neighbor (boundary cell): use that neighbor
2. If cell has two neighbors: prefer interior neighbor (toward substrate)

The interior neighbor is preferred for thermal stability since the substrate
side typically has more stable thermal conditions.

# Arguments
- `mesh`: Current mesh
- `cell_id`: ID of the thin cell

# Returns
Neighbor cell ID, or 0 if no valid neighbor exists
"""
function find_merge_neighbor(mesh::Mesh1D, cell_id::Int)
    interior = find_interior_neighbor(mesh, cell_id)
    exterior = find_exterior_neighbor(mesh, cell_id)

    # Boundary cases
    if interior == 0 && exterior == 0
        return 0  # No neighbors (shouldn't happen with valid mesh)
    elseif interior == 0
        return exterior  # At bottom boundary
    elseif exterior == 0
        return interior  # At surface
    end

    # Both neighbors exist - prefer interior for thermal stability
    return interior
end
