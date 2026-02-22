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

using ..Core
using StaticArrays

# Include all subfiles
include("nodes.jl")
include("cells.jl")
include("faces.jl")
include("mesh.jl")
include("mesh_generation.jl")
include("mesh_io.jl")
include("geometric_quantities.jl")

# ============================================================================
# Exports
# ============================================================================

# Node types
export Node1D, Node2D

# Cell types
export Cell1D, Cell2D, CellState

# Face types
export Face1D, Face2D

# Mesh types
export Mesh1D, Mesh2D

# Boundary state (for surface temperature solve)
export BoundaryState

# Mesh accessors
export n_nodes, n_cells, n_faces, n_active_cells
export get_node, get_cell, get_face, is_active

# Mesh generation
export create_uniform_mesh, create_stretched_mesh

# Geometry operations
export update_geometry!, update_geometry_from_nodes!
export cell_volume, face_area, cell_center
export total_volume, domain_bounds, boundary_face_ids, cell_size
export boundary_cell_id, face_area_value

# ALE geometry-from-state
export GeometryCache
export compute_geometry_from_state!
export compute_face_distances_from_state!
export update_geometry_cache!
export compute_gcl_volume_rate, verify_gcl_compliance
export resize_geometry_cache!

# I/O
export save_mesh, load_mesh

end # module Geometry
