"""
MeshTopology — immutable connectivity.

A `MeshTopology{D, NC}` carries everything that does **not** change once the
mesh is built: nodes, cells, faces, boundary tags, and dimensional metadata.
Mutable per-step state (current node positions, cell volumes, face areas,
mesh velocities) lives in [`MeshGeometry`](@ref) inside the `Workspace`.

This split is the structural prerequisite for the unified variant surface
and the Jacobian backends (geometry derivatives need a clean
distinction between "what changes" and "what doesn't").

# Type Parameters
- `D`: Spatial dimension (1 today; 2/3 are the design target)
- `NC`: Number of components (carried for compile-time dispatch by consumers)

# Fields
- `nodes::Vector{NodeT}`: Node objects (Node1D for D=1)
- `cells::Vector{CellT}`: Cell objects
- `faces::Vector{FaceT}`: Face objects
- `boundary_tags::Dict{Symbol,Vector{Int}}`: tag => face IDs (e.g. `:top`, `:bottom`)
- `n_active_cells::Int`: Active (leaf) cell count for AMR
"""
struct MeshTopology{D,NC,NodeT,CellT,FaceT}
    nodes::Vector{NodeT}
    cells::Vector{CellT}
    faces::Vector{FaceT}
    boundary_tags::Dict{Symbol,Vector{Int}}
    n_active_cells::Int
end

"""
	MeshTopology{D,NC}(nodes, cells, faces, boundary_tags; n_active_cells=length(cells))

Direct constructor. `boundary_tags` is the user-facing tag-keyed map.
For 1D, the canonical tags are `:top` and `:bottom`.
"""
function MeshTopology{D,NC}(
    nodes::Vector{NodeT},
    cells::Vector{CellT},
    faces::Vector{FaceT},
    boundary_tags::Dict{Symbol,Vector{Int}};
    n_active_cells::Int = length(cells),
) where {D,NC,NodeT,CellT,FaceT}
    MeshTopology{D,NC,NodeT,CellT,FaceT}(nodes, cells, faces, boundary_tags, n_active_cells)
end

# ============================================================================
# Accessors
# ============================================================================

@inline n_nodes(t::MeshTopology) = length(t.nodes)
@inline n_cells(t::MeshTopology) = length(t.cells)
@inline n_faces(t::MeshTopology) = length(t.faces)

"""
	boundary_id_to_tag(id::Int) -> Symbol

Maps a boundary ID back to its tag. Returns `:interior` for unknown IDs.
"""
@inline function boundary_id_to_tag(id::Int)
    id == 1 ? :top : id == 2 ? :bottom : :interior
end

# ============================================================================
# Build a MeshTopology view from a Mesh1D
# ============================================================================

"""
	mesh_topology(mesh::Mesh1D{NC}) -> MeshTopology{1, NC}

Build an immutable `MeshTopology` view of an existing `Mesh1D`.
Boundary tags are derived from the integer `boundary_faces` map:
`1 => :top`, `2 => :bottom`.
"""
function mesh_topology(mesh::Mesh1D{NC}) where {NC}
    tags = Dict{Symbol,Vector{Int}}()
    for (id, faces) in mesh.boundary_faces
        tag = boundary_id_to_tag(id)
        tag === :interior && continue
        tags[tag] = faces
    end
    MeshTopology{1,NC}(
        mesh.nodes,
        mesh.cells,
        mesh.faces,
        tags;
        n_active_cells = mesh.n_active_cells,
    )
end
