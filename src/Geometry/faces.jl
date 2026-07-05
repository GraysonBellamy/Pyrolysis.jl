"""
Face type definitions for mesh interfaces.

## Coordinate Convention
- z = 0 at bottom (substrate, boundary_id = 2)
- z = L at top (exposed surface, boundary_id = 1)
- "Left" (cell_ids[1]) is the cell at lower z, "Right" (cell_ids[2]) is at higher z
- At z=0 boundary: left_id = 0 (no cell below)
- At z=L boundary: right_id = 0 (no cell above)
"""

# ============================================================================
# 1D Face
# ============================================================================

"""
	Face1D <: AbstractFace

A face (interface) in a 1D mesh.

In 1D, a face is a single point between two cells or at a boundary.

# Fields
- `id::Int`: Unique face identifier
- `node_id::Int`: Node index at this face
- `cell_ids::Tuple{Int,Int}`: (left, right) cell indices, 0 = boundary
  - left_id = 0: face is at z=0 (bottom boundary)
  - right_id = 0: face is at z=L (top boundary)
- `is_boundary::Bool`: True if face is on domain boundary
- `boundary_id::Int`: Boundary identifier (0=interior, 1=top at z=L, 2=bottom at z=0)

# Coordinate Convention
- "Left" cell is at lower z (toward bottom/substrate)
- "Right" cell is at higher z (toward top/surface)
"""
struct Face1D <: AbstractFace
    id::Int
    node_id::Int                 # single node in 1D
    cell_ids::Tuple{Int,Int}     # (left, right), 0 = boundary
    is_boundary::Bool
    boundary_id::Int
end

# Convenience constructors
Face1D(id::Int, node_id::Int, cell_ids::Tuple{Int,Int}) =
    Face1D(id, node_id, cell_ids, cell_ids[1] == 0 || cell_ids[2] == 0, 0)

Face1D(id::Int, node_id::Int, cell_ids::Tuple{Int,Int}, boundary_id::Int) =
    Face1D(id, node_id, cell_ids, boundary_id != 0, boundary_id)
