"""
Node type definitions for mesh vertices.

## Coordinate Convention
- z = 0 at bottom (substrate, boundary_id = 2)
- z = L at top (exposed surface, boundary_id = 1)
- For pyrolysis, heat enters from top, material shrinks toward bottom
"""

# ============================================================================
# 1D Node
# ============================================================================

"""
	Node1D <: AbstractNode

A node (vertex) in a 1D mesh.

# Fields
- `id::Int`: Unique node identifier
- `z::Float64`: Current position [m]
- `z_ref::Float64`: Reference (initial) position [m]
- `is_boundary::Bool`: True if node is on domain boundary
- `boundary_id::Int`: Boundary identifier (0=interior, 1=top at z=L, 2=bottom at z=0)

# Coordinate Convention
- z = 0: Bottom boundary (substrate, boundary_id = 2)
- z = L: Top boundary (exposed surface, boundary_id = 1)
"""
mutable struct Node1D <: AbstractNode
    id::Int
    z::Float64              # current position
    z_ref::Float64          # reference (initial) position
    is_boundary::Bool
    boundary_id::Int        # 0 = interior, 1 = top (z=L), 2 = bottom (z=0)
end

# Convenience constructor for interior nodes
Node1D(id::Int, z::Real) = Node1D(id, Float64(z), Float64(z), false, 0)

# Constructor for boundary nodes
Node1D(id::Int, z::Real, boundary_id::Int) =
    Node1D(id, Float64(z), Float64(z), boundary_id != 0, boundary_id)
