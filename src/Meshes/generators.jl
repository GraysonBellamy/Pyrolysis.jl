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
