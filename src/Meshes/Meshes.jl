"""
	Meshes

Face-based mesh description of the finite-volume core.

Faces carry an area, a unit normal, their two cells and a boundary tag; cells
carry a volume, their faces and their nodes; nodes carry positions, which are
the geometric state. Every metric is computed from the node positions by
[`geometry!`](@ref), the only writer of a [`MeshGeometry`](@ref), so the metric
element type follows the positions and automatic differentiation passes through
the geometry.

The 1D generators build planar, cylindrical and spherical meshes from node
positions over one topology; the coordinate system decides how face areas and
cell volumes follow from the positions.
"""
module Meshes

using StaticArrays

include("topology.jl")
include("generators.jl")
include("geometry.jl")

export CoordinateSystem, Planar, Cylindrical, Spherical
export geometric_dimension, face_area, shell_volume
export Face, Cell, MeshTopology, Mesh
export n_cells, n_faces, n_nodes, boundary_faces, boundary_tags, is_fixed
export line_mesh
export MeshGeometry, geometry!, cell_distance, center_to_face_distance

end # module Meshes
