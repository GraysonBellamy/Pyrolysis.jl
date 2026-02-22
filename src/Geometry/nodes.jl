"""
Node type definitions for mesh vertices.

Reference: Design Document Section 6.3

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

# ============================================================================
# 2D Node (Axisymmetric)
# ============================================================================

"""
	Node2D <: AbstractNode

A node (vertex) in a 2D axisymmetric mesh.

# Fields
- `id::Int`: Unique node identifier
- `r::Float64`: Current radial position [m]
- `z::Float64`: Current axial position [m]
- `r_ref::Float64`: Reference radial position [m]
- `z_ref::Float64`: Reference axial position [m]
- `is_boundary::Bool`: True if node is on domain boundary
- `boundary_id::Int`: Boundary identifier (0=interior, 1=top, 2=bottom, 3=axis, 4=outer)
"""
struct Node2D <: AbstractNode
	id::Int
	r::Float64
	z::Float64
	r_ref::Float64
	z_ref::Float64
	is_boundary::Bool
	boundary_id::Int        # 0=interior, 1=top, 2=bottom, 3=axis, 4=outer
end

# Convenience constructor for interior nodes
Node2D(id::Int, r::Real, z::Real) =
	Node2D(id, Float64(r), Float64(z), Float64(r), Float64(z), false, 0)

# Constructor for boundary nodes
Node2D(id::Int, r::Real, z::Real, boundary_id::Int) =
	Node2D(id, Float64(r), Float64(z), Float64(r), Float64(z), boundary_id != 0, boundary_id)

# ============================================================================
# Node Accessors
# ============================================================================

"""
	position(node::Node1D)

Get the current position of a 1D node.
"""
@inline position(node::Node1D) = node.z

"""
	position(node::Node2D)

Get the current position of a 2D node as (r, z).
"""
@inline position(node::Node2D) = (node.r, node.z)

"""
	reference_position(node::Node1D)

Get the reference position of a 1D node.
"""
@inline reference_position(node::Node1D) = node.z_ref

"""
	reference_position(node::Node2D)

Get the reference position of a 2D node as (r, z).
"""
@inline reference_position(node::Node2D) = (node.r_ref, node.z_ref)

"""
	displacement(node::Node1D)

Get the displacement from reference position for a 1D node.
"""
@inline displacement(node::Node1D) = node.z - node.z_ref

"""
	displacement(node::Node2D)

Get the displacement from reference position for a 2D node.
"""
@inline displacement(node::Node2D) = (node.r - node.r_ref, node.z - node.z_ref)
