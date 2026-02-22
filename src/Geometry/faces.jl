"""
Face type definitions for mesh interfaces.

Reference: Design Document Section 6.3

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
	cell_ids::Tuple{Int, Int}     # (left, right), 0 = boundary
	is_boundary::Bool
	boundary_id::Int
end

# Convenience constructors
Face1D(id::Int, node_id::Int, cell_ids::Tuple{Int, Int}) =
	Face1D(id, node_id, cell_ids, cell_ids[1] == 0 || cell_ids[2] == 0, 0)

Face1D(id::Int, node_id::Int, cell_ids::Tuple{Int, Int}, boundary_id::Int) =
	Face1D(id, node_id, cell_ids, boundary_id != 0, boundary_id)

# ============================================================================
# 2D Face
# ============================================================================

"""
	Face2D <: AbstractFace

A face (edge) in a 2D mesh.

# Fields
- `id::Int`: Unique face identifier
- `node_ids::Tuple{Int,Int}`: Two nodes defining the face
- `cell_ids::Tuple{Int,Int}`: (owner, neighbor) cell indices, 0 = boundary
- `normal::Tuple{Float64,Float64}`: Unit outward normal (from owner cell)
- `is_boundary::Bool`: True if face is on domain boundary
- `boundary_id::Int`: Boundary identifier
"""
struct Face2D <: AbstractFace
	id::Int
	node_ids::Tuple{Int, Int}     # two nodes defining face
	cell_ids::Tuple{Int, Int}     # (owner, neighbor), 0 = boundary
	normal::Tuple{Float64, Float64}  # unit outward normal (from owner)
	is_boundary::Bool
	boundary_id::Int
end

# ============================================================================
# Face Accessors
# ============================================================================

"""
	owner_cell(face::Face1D)

Get the "owner" (left) cell index.
"""
@inline owner_cell(face::Face1D) = face.cell_ids[1]

"""
	neighbor_cell(face::Face1D)

Get the "neighbor" (right) cell index.
"""
@inline neighbor_cell(face::Face1D) = face.cell_ids[2]

"""
	owner_cell(face::Face2D)

Get the owner cell index.
"""
@inline owner_cell(face::Face2D) = face.cell_ids[1]

"""
	neighbor_cell(face::Face2D)

Get the neighbor cell index.
"""
@inline neighbor_cell(face::Face2D) = face.cell_ids[2]

"""
	is_interior(face::Union{Face1D,Face2D})

Check if face is an interior face (not on boundary).
"""
@inline is_interior(face::Union{Face1D, Face2D}) = !face.is_boundary

"""
	left_cell(face::Face1D)

Get the left cell index (alias for owner_cell).
"""
@inline left_cell(face::Face1D) = face.cell_ids[1]

"""
	right_cell(face::Face1D)

Get the right cell index (alias for neighbor_cell).
"""
@inline right_cell(face::Face1D) = face.cell_ids[2]
