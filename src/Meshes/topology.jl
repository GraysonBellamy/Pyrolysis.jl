# ============================================================================
# Coordinate systems
# ============================================================================

"""
	CoordinateSystem

How the face areas and cell volumes of a mesh follow from its node positions.
In 1D the position is the coordinate `x` along the column (planar) or the
radius `r` (cylindrical, spherical). The mesh's cross-section measure `A`
completes the metric: the face area in planar geometry, the axial length of
a cylinder, the solid angle of a sphere.
"""
abstract type CoordinateSystem end

"Planar coordinates: face area `A`, cell volume `A (x₊ − x₋)`."
struct Planar <: CoordinateSystem end

"Cylindrical radial coordinates: face area `2π r A`, cell volume `π (r₊² − r₋²) A`."
struct Cylindrical <: CoordinateSystem end

"Spherical radial coordinates: face area `A r²`, cell volume `A (r₊³ − r₋³) / 3`."
struct Spherical <: CoordinateSystem end

"""
	geometric_dimension(coordinates) -> Int

The exponent `d` of the coordinate system: face areas scale as `r^(d − 1)` and
every cell satisfies `V = (1/d) Σ_f A_f (x_f · n_f)`.
"""
geometric_dimension(::Planar) = 1
geometric_dimension(::Cylindrical) = 2
geometric_dimension(::Spherical) = 3

"""
	face_area(coordinates, r, A)

Area of the face at position `r` for the cross-section measure `A`.
"""
@inline face_area(::Planar, r, A) = A + zero(r)
@inline face_area(::Cylindrical, r, A) = 2π * r * A
@inline face_area(::Spherical, r, A) = A * r^2

"""
	shell_volume(coordinates, r₋, r₊, A)

Volume between the faces at `r₋` and `r₊` for the cross-section measure `A`.
Negative when the faces are inverted; the state mapping floors it.
"""
@inline shell_volume(::Planar, r0, r1, A) = A * (r1 - r0)
@inline shell_volume(::Cylindrical, r0, r1, A) = π * (r1^2 - r0^2) * A
@inline shell_volume(::Spherical, r0, r1, A) = A * (r1^3 - r0^3) / 3

# The cross-section measure a mesh gets when the caller names none: unit area,
# unit axial length, the full sphere.
default_cross_section(::Planar) = 1.0
default_cross_section(::Cylindrical) = 1.0
default_cross_section(::Spherical) = 4π

# ============================================================================
# Faces, cells, topology
# ============================================================================

"""
	Face{D,K}

A face of a `D`-dimensional mesh with `K` nodes (one in 1D).

`cells = (owner, neighbour)`; the neighbour is `0` on a boundary face. The
face normal points from the owner to the neighbour, so on a boundary face it
points out of the domain. `tag` is `:interior` or the boundary tag the
boundary conditions are keyed by (`:top`, `:bottom`, …).
"""
struct Face{D,K}
    nodes::NTuple{K,Int}
    cells::Tuple{Int,Int}
    tag::Symbol
end

"""
	Cell{D,F,K}

A cell of a `D`-dimensional mesh with `F` faces and `K` nodes (two of each in 1D).
"""
struct Cell{D,F,K}
    faces::NTuple{F,Int}
    nodes::NTuple{K,Int}
end

"""
	MeshTopology{D,C,FT,CT}

Immutable connectivity of a mesh: faces, cells, the node count, the boundary
tags (tag → face ids), the coordinate system and the cross-section measure
`A₀`.

`pinned` lists the nodes whose velocity is zero; a mesh with every node pinned
is fixed. `lagrangian` lists the nodes that always move with the solid (layer
interfaces); with the mesh velocity equal to the solid velocity every node does.
"""
struct MeshTopology{D,C<:CoordinateSystem,FT<:Face{D},CT<:Cell{D}}
    coordinates::C
    faces::Vector{FT}
    cells::Vector{CT}
    n_nodes::Int
    boundary_tags::Dict{Symbol,Vector{Int}}
    pinned::Vector{Int}
    lagrangian::Vector{Int}
    cross_section::Float64
end

@inline n_cells(t::MeshTopology) = length(t.cells)
@inline n_faces(t::MeshTopology) = length(t.faces)
@inline n_nodes(t::MeshTopology) = t.n_nodes

"""
	boundary_tags(topology) -> Vector{Symbol}

The boundary tags of the mesh, sorted.
"""
boundary_tags(t::MeshTopology) = sort!(collect(keys(t.boundary_tags)))

"""
	boundary_faces(topology, tag) -> Vector{Int}

The face ids carrying boundary tag `tag`.
"""
function boundary_faces(t::MeshTopology, tag::Symbol)
    haskey(t.boundary_tags, tag) ||
        throw(ArgumentError("mesh has no boundary tagged :$tag (has $(boundary_tags(t)))"))
    return t.boundary_tags[tag]
end

"""
	is_fixed(topology) -> Bool

Whether every node is pinned, so the geometry never changes.
"""
is_fixed(t::MeshTopology) = length(t.pinned) == t.n_nodes

"""
	Mesh{D,TopT}

A mesh: its topology and the reference positions of its nodes, the positions
at which the problem's initial conditions are evaluated.
"""
struct Mesh{D,TopT<:MeshTopology{D}}
    topology::TopT
    positions::Vector{SVector{D,Float64}}
end
