# ============================================================================
# 1D generators
# ============================================================================

"""
	line_mesh(x; coordinates = Planar(), cross_section = <default>,
	          bottom = :bottom, top = :top, fixed = false, lagrangian = Int[]) -> Mesh{1}

A 1D mesh with nodes at the strictly increasing positions `x` (a coordinate in
planar geometry, a non-negative radius in cylindrical and spherical geometry).

Cell `i` lies between nodes `i` and `i + 1`; face `f` sits at node `f`. The
face at `x[1]` is tagged `bottom`, the face at `x[end]` is tagged `top`, and
every other face is `:interior` with its normal pointing toward increasing
`x`. `cross_section` is the planar face area, the axial length of a cylinder
or the solid angle of a sphere (defaults 1, 1 and 4π).

The back node `x[1]` is pinned. With `fixed = true` every node is pinned and
the geometry never changes. `lagrangian` lists the nodes that always move with
the solid.
"""
function line_mesh(
    x::AbstractVector{<:Real};
    coordinates::CoordinateSystem = Planar(),
    cross_section::Real = default_cross_section(coordinates),
    bottom::Symbol = :bottom,
    top::Symbol = :top,
    fixed::Bool = false,
    lagrangian::AbstractVector{<:Integer} = Int[],
)
    n_nodes = length(x)
    n_nodes >= 2 || throw(ArgumentError("a line mesh needs at least two nodes"))
    all(isfinite, x) || throw(ArgumentError("node positions must be finite"))
    for n = 2:n_nodes
        x[n] > x[n-1] || throw(
            ArgumentError("node positions must increase strictly (x[$n] ≤ x[$(n-1)])"),
        )
    end
    if !(coordinates isa Planar) && x[1] < 0
        throw(ArgumentError("a radial mesh needs non-negative radii (x[1] = $(x[1]))"))
    end
    cross_section > 0 || throw(ArgumentError("cross_section must be positive"))
    bottom === top && throw(ArgumentError("the two boundary tags must differ"))
    (bottom === :interior || top === :interior) &&
        throw(ArgumentError(":interior is not a boundary tag"))

    n_cells = n_nodes - 1
    faces = Vector{Face{1,1}}(undef, n_nodes)
    faces[1] = Face{1,1}((1,), (1, 0), bottom)
    for f = 2:n_cells
        faces[f] = Face{1,1}((f,), (f - 1, f), :interior)
    end
    faces[n_nodes] = Face{1,1}((n_nodes,), (n_cells, 0), top)
    cells = [Cell{1,2,2}((i, i + 1), (i, i + 1)) for i = 1:n_cells]

    tags = Dict{Symbol,Vector{Int}}(bottom => [1], top => [n_nodes])
    pinned = fixed ? collect(1:n_nodes) : [1]
    lag = sort!(unique(Int.(lagrangian)))
    all(n -> 1 <= n <= n_nodes, lag) ||
        throw(ArgumentError("lagrangian node ids must lie in 1:$n_nodes"))

    topology = MeshTopology{1,typeof(coordinates),Face{1,1},Cell{1,2,2}}(
        coordinates,
        faces,
        cells,
        n_nodes,
        tags,
        pinned,
        lag,
        Float64(cross_section),
    )
    return Mesh(topology, [SVector{1,Float64}(xi) for xi in x])
end

# ============================================================================
# Quadrilaterals and hexahedra
# ============================================================================

# The local corners of side `s ∈ (0, 1)` of axis `k` of the reference cell, in
# the tensor order of the remaining axes, transposed where needed so that the
# side's area vector `∫ ∂x/∂s × ∂x/∂t` points out of the cell. The tensor
# order's area vector points along `(−1)^(k−1) e_k`; the outward normal is
# `(2s − 1) e_k`.
function _side_corners(D::Int, k::Int, s::Int)
    corners = [a for a = 1:(1<<D) if ((a - 1) >> (k - 1)) & 1 == s]
    if (2s - 1) != (isodd(k) ? 1 : -1)
        if D == 2
            reverse!(corners)
        else
            corners[2], corners[3] = corners[3], corners[2]
        end
    end
    return corners
end

"""
	cell_mesh(positions, cells; tag, pinned = (:bottom,), fixed = false,
	          lagrangian = Int[]) -> Mesh{D}

A planar mesh of quadrilaterals (`D = 2`) or hexahedra (`D = 3`) from the
node positions (`SVector{D}` each) and, per cell, its `2^D` node ids in the
tensor order of the reference cell `[0, 1]^D`, the first axis fastest (in 2D
`(0, 0), (1, 0), (0, 1), (1, 1)`). Every cell must be positively oriented.

Faces are the cells' sides, matched by their node sets. A side two cells
share is one interior face owned by the lower-numbered cell, its normal
pointing from the owner to the neighbour; a side no other cell shares is a
boundary face, its normal outward, tagged `tag(i, k, s)` for side
`s ∈ (0, 1)` of axis `k` of cell `i`. A face's nodes are in the tensor order
of its side's remaining axes, transposed where needed so that its area
vector points along its normal.

The nodes of the faces tagged with a tag in `pinned` are pinned; with
`fixed = true` every node is. `lagrangian` lists the nodes that always move
with the solid. The cross-section measure is 1: 2D faces and cells are per
unit depth.
"""
function cell_mesh(
    positions::AbstractVector{<:SVector{D,<:Real}},
    cells::AbstractVector{<:NTuple{K,Integer}};
    tag,
    pinned = (:bottom,),
    fixed::Bool = false,
    lagrangian::AbstractVector{<:Integer} = Int[],
) where {D,K}
    D in (2, 3) || throw(ArgumentError("cell meshes are two- or three-dimensional"))
    K == 1 << D || throw(ArgumentError("a $(D)D cell has $(1 << D) nodes, not $K"))
    nn = length(positions)
    all(p -> all(isfinite, p), positions) ||
        throw(ArgumentError("node positions must be finite"))
    isempty(cells) && throw(ArgumentError("a cell mesh needs at least one cell"))
    for (i, c) in enumerate(cells)
        all(n -> 1 <= n <= nn, c) ||
            throw(ArgumentError("cell $i names a node outside 1:$nn"))
        allunique(c) || throw(ArgumentError("cell $i repeats a node"))
        p = ntuple(a -> SVector{D,Float64}(positions[c[a]]), K)
        cell_volume_quadrature(p) > 0 ||
            throw(ArgumentError("cell $i is not positively oriented"))
    end

    KF = K ÷ 2
    sides = [_side_corners(D, k, s) for k = 1:D for s = 0:1]
    FT = Face{D,KF}
    faces = FT[]
    owner_side = Tuple{Int,Int,Int}[]
    lookup = Dict{NTuple{KF,Int},Int}()
    cell_faces = Vector{NTuple{2D,Int}}(undef, length(cells))
    for (i, c) in enumerate(cells)
        ids = ntuple(2D) do e
            k, s = (e - 1) ÷ 2 + 1, (e - 1) % 2
            nodes = ntuple(q -> Int(c[sides[e][q]]), KF)
            key = Tuple(sort!(collect(nodes)))
            f = get(lookup, key, 0)
            if f == 0
                push!(faces, FT(nodes, (i, 0), :interior))
                push!(owner_side, (i, k, s))
                f = lookup[key] = length(faces)
            else
                o, nb = faces[f].cells
                nb == 0 ||
                    throw(ArgumentError("a face is shared by more than two cells"))
                faces[f] = FT(faces[f].nodes, (o, i), :interior)
            end
            f
        end
        cell_faces[i] = ids
    end

    tags = Dict{Symbol,Vector{Int}}()
    for f in eachindex(faces)
        faces[f].cells[2] == 0 || continue
        t = tag(owner_side[f]...)
        t isa Symbol && t !== :interior ||
            throw(ArgumentError("boundary tags must be symbols other than :interior"))
        faces[f] = FT(faces[f].nodes, faces[f].cells, t)
        push!(get!(tags, t, Int[]), f)
    end
    for t in pinned
        haskey(tags, t) || throw(ArgumentError("cannot pin :$t, which tags no face"))
    end
    pin =
        fixed ? collect(1:nn) :
        sort!(unique([n for t in pinned for f in tags[t] for n in faces[f].nodes]))
    lag = sort!(unique(Int.(lagrangian)))
    all(n -> 1 <= n <= nn, lag) ||
        throw(ArgumentError("lagrangian node ids must lie in 1:$nn"))

    CT = Cell{D,2D,K}
    topology = MeshTopology{D,Planar,FT,CT}(
        Planar(),
        faces,
        [CT(cell_faces[i], ntuple(a -> Int(cells[i][a]), K)) for i in eachindex(cells)],
        nn,
        tags,
        pin,
        lag,
        1.0,
    )
    return Mesh(topology, [SVector{D,Float64}(p) for p in positions])
end

# The boundary tags of a box mesh, per axis and side: x, then y, then z, the
# last axis running from the back (`:bottom`) to the exposed face (`:top`).
_box_tags(::Val{2}) = ((:left, :right), (:bottom, :top))
_box_tags(::Val{3}) = ((:left, :right), (:front, :back), (:bottom, :top))

"""
	box_mesh(coordinates...; transform = identity, pinned = (:bottom,), fixed = false,
	         lagrangian = Int[]) -> Mesh{D}

A structured mesh of quadrilaterals (two coordinate vectors) or hexahedra
(three) on the tensor product of the strictly increasing node coordinates of
each axis, every node then moved by `transform(x::SVector) -> SVector`.
Nodes and cells are numbered with the first axis fastest. The sides are
tagged `:left` and `:right` along x, then `:bottom` and `:top` along the last
axis, and in 3D `:front` and `:back` along y. Keywords as [`cell_mesh`](@ref).
"""
function box_mesh(
    coordinates::Vararg{AbstractVector{<:Real},D};
    transform = identity,
    pinned = (:bottom,),
    fixed::Bool = false,
    lagrangian::AbstractVector{<:Integer} = Int[],
) where {D}
    D in (2, 3) || throw(ArgumentError("box meshes are two- or three-dimensional"))
    for (k, x) in enumerate(coordinates)
        length(x) >= 2 || throw(ArgumentError("axis $k needs at least two nodes"))
        all(n -> x[n] > x[n-1], 2:length(x)) ||
            throw(ArgumentError("the node coordinates of axis $k must increase strictly"))
    end
    dims = ntuple(k -> length(coordinates[k]), Val(D))
    ids = LinearIndices(dims)
    positions = [
        SVector{D,Float64}(
            transform(SVector{D,Float64}(ntuple(k -> coordinates[k][I[k]], Val(D)))),
        ) for I in CartesianIndices(dims)
    ][:]
    cells = [
        ntuple(
            a -> ids[CartesianIndex(ntuple(k -> I[k] + ((a - 1) >> (k - 1)) & 1, Val(D)))],
            1 << D,
        ) for I in CartesianIndices(dims .- 1)
    ][:]
    names = _box_tags(Val(D))
    return cell_mesh(
        positions,
        cells;
        tag = (i, k, s) -> names[k][s+1],
        pinned,
        fixed,
        lagrangian,
    )
end
