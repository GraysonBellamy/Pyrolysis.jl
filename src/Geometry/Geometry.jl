"""
	Geometry

Mesh and geometry infrastructure for the Pyrolysis.jl package.

This module provides:
- Node, cell, and face types for 1D and 2D meshes
- Mesh container types
- Mesh generation utilities
- Geometric quantity computations
"""
module Geometry

using ..Materials
using ..BoundaryConditions
using StaticArrays

# Include all subfiles
include("nodes.jl")
include("cells.jl")
include("faces.jl")
include("mesh.jl")
include("mesh_generation.jl")
include("geometric_quantities.jl")
# Pure topology / geometry split: Mesh1D carries mutable per-step state, while
# MeshTopology holds the immutable connectivity and MeshGeometry holds the
# mutable z/V/A buffers consumed by the residual / Jacobian.
include("topology.jl")
include("mesh_geometry.jl")

# ============================================================================
# Exports
# ============================================================================

# Node types
export Node1D

# Cell types
export Cell1D, CellState

# Face types
export Face1D

# Mesh types
export Mesh1D

# `MeshTopology` (immutable connectivity view) and `MeshGeometry` (mutable
# z/V/A buffers) are the user-facing types for the pure-residual layer.
# Builders / boundary-id translators (`mesh_topology`, `mesh_geometry`,
# `sync_geometry_from_mesh!`, `tag_to_boundary_id`, `boundary_id_to_tag`)
# are deliberately internal — field access on the structs is the
# canonical read path; the two exported accessors below dereference
# `Ref`-wrapped scalars and aren't pure field reads.
export MeshTopology, MeshGeometry
export resize_geometry!, cross_section_area, mesh_time

# Boundary state (for surface temperature solve)
export BoundaryState

# Mesh accessors
export n_nodes, n_cells, n_faces

# Mesh generation
export create_uniform_mesh, create_stretched_mesh

# Geometry operations
export update_geometry!
export cell_volume, face_area, cell_center
export total_volume, domain_bounds, cell_size
# Second-order cell-centered gradient on nonuniform stencils
export opposite_distance_gradient, opposite_distance_gradient_coeffs

# ALE geometry-from-state
export GeometryCache
export resize_geometry_cache!

end # module Geometry
