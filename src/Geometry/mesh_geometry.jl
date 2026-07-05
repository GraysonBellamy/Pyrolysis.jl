"""
MeshGeometry — mutable per-step geometry.

A `MeshGeometry{D}` carries the geometric quantities that change as the mesh
moves (ALE) or is regenerated (AMR): node positions, cell volumes, face areas
and (for D ≥ 2) face normals. The cross-sectional area `A(t)` for the 1D
specialization lives here as a scalar.

Geometry instances live inside the `Workspace` so the immutable
[`MeshTopology`](@ref) does not have to be rebuilt every time `z` or `V`
changes.

# Type Parameters
- `D`: Spatial dimension (1 today; 2/3 are the design target)

# Fields (1D specialization layout)
- `node_positions::Vector{Float64}`: current `z` per node
- `node_positions_ref::Vector{Float64}`: reference (initial) `z` per node
- `cell_centers::Vector{Float64}`: cell center positions
- `cell_volumes::Vector{Float64}`: cell volumes
- `face_positions::Vector{Float64}`: face positions
- `face_areas::Vector{Float64}`: face cross-sectional areas
- `node_velocities::Vector{Float64}`: ALE mesh velocities
- `mesh_time::Base.RefValue{Float64}`: timestamp of the geometry snapshot
- `cross_section_area::Base.RefValue{Float64}`: current `A(t)` (scalar in 1D)
- `cross_section_area_initial::Float64`: calibration reference `A_0`

For `D ≥ 2`, `node_positions` and `face_normals` become `Vector{SVector{D,Float64}}`
and `cross_section_area*` becomes a per-face scalar.
"""
mutable struct MeshGeometry{D}
    node_positions::Vector{Float64}
    node_positions_ref::Vector{Float64}
    cell_centers::Vector{Float64}
    cell_volumes::Vector{Float64}
    face_positions::Vector{Float64}
    face_areas::Vector{Float64}
    node_velocities::Vector{Float64}
    mesh_time::Base.RefValue{Float64}
    cross_section_area::Base.RefValue{Float64}
    cross_section_area_initial::Float64
end

function MeshGeometry{D}(
    n_nodes::Int,
    n_cells::Int,
    n_faces::Int;
    cross_section_area::Real = 1.0,
) where {D}
    D == 1 || throw(ArgumentError("MeshGeometry: only D=1 is implemented (got D=$D)"))
    A0 = Float64(cross_section_area)
    MeshGeometry{D}(
        zeros(Float64, n_nodes),
        zeros(Float64, n_nodes),
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, n_faces),
        fill(A0, n_faces),
        zeros(Float64, n_nodes),
        Ref(0.0),
        Ref(A0),
        A0,
    )
end

# ============================================================================
# Accessors
# ============================================================================

@inline n_nodes(g::MeshGeometry) = length(g.node_positions)
@inline n_cells(g::MeshGeometry) = length(g.cell_volumes)
@inline n_faces(g::MeshGeometry) = length(g.face_positions)
@inline cross_section_area(g::MeshGeometry) = g.cross_section_area[]
@inline mesh_time(g::MeshGeometry) = g.mesh_time[]

"""
	resize_geometry!(g::MeshGeometry, n_nodes_new, n_cells_new, n_faces_new)

Resize all internal buffers in-place. Existing values are preserved up to
the new size; new entries are zero. Used after AMR or h-refinement.
"""
function resize_geometry!(
    g::MeshGeometry,
    n_nodes_new::Int,
    n_cells_new::Int,
    n_faces_new::Int,
)
    resize!(g.node_positions, n_nodes_new)
    resize!(g.node_positions_ref, n_nodes_new)
    resize!(g.cell_centers, n_cells_new)
    resize!(g.cell_volumes, n_cells_new)
    resize!(g.face_positions, n_faces_new)
    resize!(g.face_areas, n_faces_new)
    resize!(g.node_velocities, n_nodes_new)
    return g
end

# ============================================================================
# Adapter: Mesh1D -> MeshGeometry{1}
# ============================================================================

"""
	mesh_geometry(mesh::Mesh1D) -> MeshGeometry{1}

Build a `MeshGeometry{1}` view from an existing `Mesh1D` by **copying** the
mutable geometry fields. Does not alias the input — call
[`sync_geometry_from_mesh!`](@ref) to refresh after the underlying mesh
mutates.
"""
function mesh_geometry(mesh::Mesh1D)
    g = MeshGeometry{1}(
        length(mesh.nodes),
        length(mesh.cells),
        length(mesh.faces);
        cross_section_area = mesh.cross_section_area_initial,
    )
    sync_geometry_from_mesh!(g, mesh)
    return g
end

"""
	sync_geometry_from_mesh!(g::MeshGeometry, mesh::Mesh1D)

Refresh `g`'s buffers from the `Mesh1D`'s mutable fields. Called after
any operation that mutates `mesh` in place (ALE step, AMR refinement)
to keep the immutable `MeshGeometry` view consistent with the
underlying mutable mesh.
"""
function sync_geometry_from_mesh!(g::MeshGeometry, mesh::Mesh1D)
    if length(g.node_positions) != length(mesh.nodes) ||
       length(g.cell_volumes) != length(mesh.cells) ||
       length(g.face_areas) != length(mesh.faces)
        resize_geometry!(g, length(mesh.nodes), length(mesh.cells), length(mesh.faces))
    end
    @inbounds for i in eachindex(mesh.nodes)
        g.node_positions[i] = mesh.nodes[i].z
        g.node_positions_ref[i] = mesh.nodes[i].z_ref
    end
    @inbounds for i in eachindex(mesh.cells)
        g.cell_centers[i] = mesh.cell_centers[i]
        g.cell_volumes[i] = mesh.cell_volumes[i]
    end
    @inbounds for i in eachindex(mesh.faces)
        g.face_positions[i] = mesh.face_positions[i]
        g.face_areas[i] = mesh.face_areas[i]
    end
    if length(mesh.node_velocities) == length(g.node_velocities)
        copyto!(g.node_velocities, mesh.node_velocities)
    end
    g.mesh_time[] = mesh.mesh_time
    g.cross_section_area[] = mesh.cross_section_area
    return g
end
