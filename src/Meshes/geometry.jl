# ============================================================================
# Per-step metric buffers
# ============================================================================

"""
	MeshGeometry{D,T}

The metrics of a mesh at one set of node positions: node positions, face
areas, unit face normals, face centers and sweep vectors, cell volumes and
cell centers, and the cells' lumped corner volumes. The element type `T` follows the
positions, so a dual-number state carries geometry sensitivities.
[`geometry!`](@ref) is the only writer.

`face_sweep[q, f]` is `∫_f N_q n dA` for the `q`-th node of face `f`, with
`N_q` the node's shape function on the face: the rate at which the face
sweeps volume along its normal per unit velocity of that node, so a face
moving with node velocities `w` sweeps `Σ_q face_sweep[q, f] · w_q`, and the
sweep vectors of a face sum to its area vector `A_f n_f`.
`cell_node_volume[a, i]` is `∫_i N_a dV` for the `a`-th node of cell `i`: a
node's lumped volume is the sum over the cells around it, or over the live
ones.

Volumes are the signed shell volumes between a cell's faces: an inverted
trial mesh yields a non-positive volume, which the state mapping floors for
evaluation and the admissibility check reads unfloored.
"""
struct MeshGeometry{D,T}
    node_position::Vector{SVector{D,T}}
    face_area::Vector{T}
    face_normal::Vector{SVector{D,T}}
    face_center::Vector{SVector{D,T}}
    face_sweep::Matrix{SVector{D,T}}
    cell_volume::Vector{T}
    cell_center::Vector{SVector{D,T}}
    cell_node_volume::Matrix{T}
end

"""
	MeshGeometry{T}(topology) -> MeshGeometry{D,T}

Zero-filled metric buffers sized for `topology`.
"""
function MeshGeometry{T}(
    topo::MeshTopology{D,C,<:Face{D,K},<:Cell{D,F,KC}},
) where {T,D,C,K,F,KC}
    z = zero(SVector{D,T})
    return MeshGeometry{D,T}(
        fill(z, n_nodes(topo)),
        zeros(T, n_faces(topo)),
        fill(z, n_faces(topo)),
        fill(z, n_faces(topo)),
        fill(z, K, n_faces(topo)),
        zeros(T, n_cells(topo)),
        fill(z, n_cells(topo)),
        zeros(T, KC, n_cells(topo)),
    )
end

# The two-point Gauss–Legendre abscissae on [0, 1], each of weight 1/2: exact
# for polynomials of degree three in each variable, which covers every
# integrand of the multilinear cells below.
const _GAUSS = ((1 - 1 / sqrt(3)) / 2, (1 + 1 / sqrt(3)) / 2)

"""
	geometry!(g::MeshGeometry{1}, x, A, topology::MeshTopology{1}) -> g

Compute every metric of a 1D mesh from the node positions `x` (one scalar per
node) and the cross-section measure `A`. Face areas and cell volumes follow
the coordinate system; cell centers are the midpoints of their two nodes;
face normals point from owner to neighbour, out of the domain on the boundary.
A face's sweep vector is its area vector; a cell's corner volumes integrate
its two linear shape functions over the shell volume.
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
        g.face_sweep[1, f] = g.face_area[f] * g.face_normal[f]
    end
    @inbounds for i = 1:n_cells(topo)
        n0, n1 = topo.cells[i].nodes
        r0, r1 = x[n0], x[n1]
        g.cell_volume[i] = shell_volume(coords, r0, r1, A)
        g.cell_center[i] = SVector{1,T}((r0 + r1) / 2)
        # dV = a(r) dr with a the face area at r, of degree two at most.
        W0, W1 = zero(T), zero(T)
        for ξ in _GAUSS
            dV = (r1 - r0) * face_area(coords, r0 + ξ * (r1 - r0), A) / 2
            W0 += (1 - ξ) * dV
            W1 += ξ * dV
        end
        g.cell_node_volume[1, i] = W0
        g.cell_node_volume[2, i] = W1
    end
    return g
end

# ============================================================================
# Quadrilaterals and hexahedra
# ============================================================================

# Corner `a` of the reference cell [0, 1]^D in tensor order, the first axis
# fastest: the bits of `a − 1`.
@inline _corner(a::Int, k::Int) = ((a - 1) >> (k - 1)) & 1

# The multilinear shape functions of the 2^D corners at the reference point ξ
# and their gradients (one row per corner).
@inline function _shape(ξ::SVector{D}) where {D}
    K = 1 << D
    N = SVector{K}(
        ntuple(a -> prod(ntuple(k -> _corner(a, k) == 1 ? ξ[k] : 1 - ξ[k], Val(D))), K),
    )
    dN = SMatrix{K,D}(
        ntuple(K * D) do e
            a, l = (e - 1) % K + 1, (e - 1) ÷ K + 1
            s = _corner(a, l) == 1 ? 1 : -1
            s * prod(
                ntuple(k -> k == l ? 1.0 : (_corner(a, k) == 1 ? ξ[k] : 1 - ξ[k]), Val(D)),
            )
        end,
    )
    return N, dN
end

# The 2^D Gauss points of the reference cell, each of weight 1/2^D.
@inline _gauss_points(::Val{D}) where {D} =
    ntuple(q -> SVector{D}(ntuple(k -> _GAUSS[_corner(q, k)+1], Val(D))), 1 << D)

# The rotation that turns a 2D edge's tangent into its area vector: the edge's
# tangent crossed with the out-of-plane unit vector.
@inline _rot(v::SVector{2}) = SVector(v[2], -v[1])

# The sweep vectors ∫_f N_q n dA of a straight edge, each half its area vector.
@inline function _face_sweeps(p::NTuple{2,SVector{2,T}}) where {T}
    S = _rot(p[2] - p[1]) / 2
    return (S, S)
end

# The sweep vectors of a bilinear quadrilateral face with corners in tensor
# order: ∫∫ N_q (∂x/∂s × ∂x/∂t) ds dt by Gauss, exact for the degree-two
# integrand.
@inline function _face_sweeps(p::NTuple{4,SVector{3,T}}) where {T}
    z = zero(SVector{3,T})
    s = (z, z, z, z)
    for ξ in _gauss_points(Val(2))
        N, dN = _shape(ξ)
        xs = p[1] * dN[1, 1] + p[2] * dN[2, 1] + p[3] * dN[3, 1] + p[4] * dN[4, 1]
        xt = p[1] * dN[1, 2] + p[2] * dN[2, 2] + p[3] * dN[3, 2] + p[4] * dN[4, 2]
        S = cross(xs, xt) / 4
        s = ntuple(q -> s[q] + N[q] * S, Val(4))
    end
    return s
end

"""
	geometry!(g::MeshGeometry{D}, x, A, topology::MeshTopology{D}) -> g    (D = 2, 3)

Compute every metric of a mesh of quadrilaterals or hexahedra from the node
positions `x` (the flat vector of the `D` coordinates of every node in turn).
Each face's sweep vectors are exact integrals over its straight or bilinear
surface; its area vector is their sum, its area the vector's length, its
normal the vector over the area (owner to neighbour by the topology's node
order) and its center the normal-weighted mean of its nodes, the centroid of
a planar face. Cell volumes follow from the divergence identity over the
sweep vectors, `Vᵢ = (1/D) Σ_f σ_f Σ_q x_q · s_{q,f}`; cell centroids and
lumped corner volumes from Gauss quadrature of the multilinear map, exact for
these integrands. `A` is not read: 2D faces and cells are per unit depth.
"""
function geometry!(
    g::MeshGeometry{D,T},
    x::AbstractVector,
    A,
    topo::MeshTopology{D},
) where {D,T}
    @inbounds for n = 1:n_nodes(topo)
        g.node_position[n] = SVector{D,T}(ntuple(k -> x[D*(n-1)+k], Val(D)))
    end
    @inbounds for f = 1:n_faces(topo)
        p = map(n -> g.node_position[n], topo.faces[f].nodes)
        s = _face_sweeps(p)
        S = sum(s)
        a = norm(S)
        nf = S / a
        g.face_area[f] = a
        g.face_normal[f] = nf
        c = zero(SVector{D,T})
        for q in eachindex(s)
            g.face_sweep[q, f] = s[q]
            c += dot(nf, s[q]) * p[q]
        end
        g.face_center[f] = c / a
    end
    @inbounds for i = 1:n_cells(topo)
        cell = topo.cells[i]
        V = zero(T)
        for f in cell.faces
            σ = topo.faces[f].cells[1] == i ? 1 : -1
            nodes = topo.faces[f].nodes
            for q in eachindex(nodes)
                V += σ * dot(g.node_position[nodes[q]], g.face_sweep[q, f])
            end
        end
        V /= D
        g.cell_volume[i] = V
        p = map(n -> g.node_position[n], cell.nodes)
        m = zero(SVector{D,T})
        for a in eachindex(p)
            g.cell_node_volume[a, i] = zero(T)
        end
        for ξ in _gauss_points(Val(D))
            N, dN = _shape(ξ)
            J = sum(p[a] * dN[a, :]' for a in eachindex(p))
            dV = det(J) / (1 << D)
            m += dV * sum(N[a] * p[a] for a in eachindex(p))
            for a in eachindex(p)
                g.cell_node_volume[a, i] += N[a] * dV
            end
        end
        g.cell_center[i] = m / V
    end
    return g
end

"""
	cell_volume_quadrature(p) -> V

The volume `∫ det J` of the multilinear cell with corners `p` in tensor order,
by Gauss quadrature of the map's Jacobian: a path independent of the face
sweep vectors, for checking them.
"""
function cell_volume_quadrature(p::Tuple{SVector{D,T},Vararg{SVector{D,T}}}) where {D,T}
    V = zero(T)
    for ξ in _gauss_points(Val(D))
        _, dN = _shape(ξ)
        J = sum(p[a] * dN[a, :]' for a in eachindex(p))
        V += det(J) / (1 << D)
    end
    return V
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
