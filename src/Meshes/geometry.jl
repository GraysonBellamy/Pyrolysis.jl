# ============================================================================
# Per-step metric buffers
# ============================================================================

"""
	MeshGeometry{D,T}

The metrics of a mesh at one set of node positions: node positions, face
areas, unit face normals and face centers, cell volumes and cell centers.
The element type `T` follows the positions, so a dual-number state carries
geometry sensitivities. [`geometry!`](@ref) is the only writer.

Volumes are the signed shell volumes between a cell's faces: an inverted
trial mesh yields a non-positive volume, which the state mapping floors for
evaluation and the admissibility check reads unfloored.
"""
struct MeshGeometry{D,T}
    node_position::Vector{SVector{D,T}}
    face_area::Vector{T}
    face_normal::Vector{SVector{D,T}}
    face_center::Vector{SVector{D,T}}
    cell_volume::Vector{T}
    cell_center::Vector{SVector{D,T}}
end

"""
	MeshGeometry{T}(topology) -> MeshGeometry{D,T}

Zero-filled metric buffers sized for `topology`.
"""
function MeshGeometry{T}(topo::MeshTopology{D}) where {T,D}
    z = zero(SVector{D,T})
    return MeshGeometry{D,T}(
        fill(z, n_nodes(topo)),
        zeros(T, n_faces(topo)),
        fill(z, n_faces(topo)),
        fill(z, n_faces(topo)),
        zeros(T, n_cells(topo)),
        fill(z, n_cells(topo)),
    )
end

"""
	geometry!(g::MeshGeometry{1}, x, A, topology::MeshTopology{1}) -> g

Compute every metric of a 1D mesh from the node positions `x` (one scalar per
node) and the cross-section measure `A`. Face areas and cell volumes follow
the coordinate system; cell centers are the midpoints of their two nodes;
face normals point from owner to neighbour, out of the domain on the boundary.
"""
function geometry!(
    g::MeshGeometry{1,T},
    x::AbstractVector,
    A,
    topo::MeshTopology{1},
) where {T}
    coords = topo.coordinates
    @inbounds for n = 1:n_nodes(topo)
        g.node_position[n] = SVector{1,T}(x[n])
    end
    @inbounds for f = 1:n_faces(topo)
        face = topo.faces[f]
        node = face.nodes[1]
        r = x[node]
        owner = face.cells[1]
        # The owner's upper node is this face when the owner lies below it.
        outward = topo.cells[owner].nodes[2] == node
        g.face_center[f] = SVector{1,T}(r)
        g.face_area[f] = face_area(coords, r, A)
        g.face_normal[f] = SVector{1,T}(outward ? one(T) : -one(T))
    end
    @inbounds for i = 1:n_cells(topo)
        n0, n1 = topo.cells[i].nodes
        r0, r1 = x[n0], x[n1]
        g.cell_volume[i] = shell_volume(coords, r0, r1, A)
        g.cell_center[i] = SVector{1,T}((r0 + r1) / 2)
    end
    return g
end

"""
	cell_distance(g, topology, f)

Distance between the centers of the two cells of interior face `f`, measured
along the face normal (positive on a valid mesh).
"""
@inline function cell_distance(g::MeshGeometry, topo::MeshTopology, f::Int)
    o, nb = topo.faces[f].cells
    return sum((g.cell_center[nb] - g.cell_center[o]) .* g.face_normal[f])
end

"""
	center_to_face_distance(g, topology, f, i)

Distance from the center of cell `i` to face `f`, measured along the normal
pointing out of cell `i` (positive on a valid mesh).
"""
@inline function center_to_face_distance(
    g::MeshGeometry,
    topo::MeshTopology,
    f::Int,
    i::Int,
)
    d = sum((g.face_center[f] - g.cell_center[i]) .* g.face_normal[f])
    return topo.faces[f].cells[1] == i ? d : -d
end
