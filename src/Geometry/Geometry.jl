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
# Mesh1D carries the 1D mesh the solver uses. MeshTopology (immutable
# connectivity) and MeshGeometry (mutable z/V/A buffers) are the
# dimension-parameterised description for a multi-dimensional extension;
# the 1D residual and Jacobian do not read them.
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

# `MeshTopology` (immutable connectivity) and `MeshGeometry` (mutable z/V/A
# buffers). Builders / boundary-id translators (`mesh_topology`,
# `mesh_geometry`, `sync_geometry_from_mesh!`, `tag_to_boundary_id`,
# `boundary_id_to_tag`) are internal; the two exported accessors below
# dereference `Ref`-wrapped scalars.
export MeshTopology, MeshGeometry
export resize_geometry!, cross_section_area, mesh_time

# Boundary state (for surface temperature solve)
export BoundaryState

# Mesh accessors
export n_nodes, n_cells, n_faces

# Mesh generation
export create_uniform_mesh, create_stretched_mesh, create_mesh_from_nodes

# Geometry operations
export update_geometry!
export domain_bounds, cell_size
# Second-order cell-centered gradient on nonuniform stencils
export opposite_distance_gradient, opposite_distance_gradient_coeffs

end # module Geometry
