# ============================================================================
# Balance-row closure
# ============================================================================

"""
	Closure

How the node velocities of a moving mesh are determined: [`closure_rows!`](@ref)
writes the residuals of their algebraic rows, [`consistent_velocities!`](@ref)
the velocities that satisfy them, and [`rank_check`](@ref) whether the rows
determine the velocities.
"""
abstract type Closure end

"""
	BalanceRowClosure{C,X} <: Closure

The node velocities of a line, or of a single column of cells, from the
discrete volume balance of every cell. The pinned node is at rest, `0 = w_p`; cell `k`'s balance sits on the
velocity row of its node farther from the pinned node,

	0 = Σ_f A_f ((w − c)·n)_f − (θ_k − θ_A) V_k,

in 1D `A_{k+1} (w − c)_{k+1} − A_k (w − c)_k − (θ_k − θ_A) V_k`, with `c` the
node corrections of the [`Correction`](@ref) (`ds.c`, zero under
[`NoCorrection`](@ref)), `θ_k` the strain rate of the cell's reaction mass
sources (`Physics.volume_change_rate` at the trial temperature) and `θ_A`
the lateral strain rate of a variable cross section: the solid velocity
`w − c` satisfies the volume balance. The rows are bidiagonal from the
pinned node and each reads one cell and its two nodes, besides what the
corrections of those nodes read, so the Jacobian stays banded where the
swept velocity would couple every node to every cell between it and the
pinned node. Since the cell volumes are computed from the same node
velocities' positions, `dVₖ/dt = θₖ Vₖ + sₖ` holds semi-discretely, with the
sweep `sₖ = Σ_f A_f (c·n)_f`.

Under a [`PrescribedCorrection`](@ref), a function `w(x, t)`, every velocity
row is `0 = w_n − w(x_n, t)` instead: a prescribed mesh velocity, for
verification, under which the mesh moves over a solid at rest. `lateral` is
the variable cross section's [`LateralShrinkage`](@ref), or `nothing`.

On a mesh of more than one dimension the closure holds a single column of
cells along the last axis in column mode: the lateral velocity components are
zero, the pinned nodes are at rest, cell `k`'s balance sits on the axial
velocity of one node of the face it adds to the column, and the face's other
nodes move with that node (`0 = w_n − w_leader`, the pairs in `ties`), so
the column's levels stay level and the rows are the 1D rows.

`pinned` lists the nodes at rest, `row_node[k]` the node whose velocity row
carries cell `k`'s balance and `order` the cells from the pinned nodes
outward.
"""
struct BalanceRowClosure{C<:Correction,X} <: Closure
    pinned::Vector{Int}
    row_node::Vector{Int}
    order::Vector{Int}
    ties::Vector{Tuple{Int,Int}}
    correction::C
    lateral::X
end

"""
	BalanceRowClosure(topology; correction = NoCorrection(), lateral = nothing)

The balance-row closure of a 1D mesh with one pinned node at an end of the
line.
"""
function BalanceRowClosure(
    topo::MeshTopology{1};
    correction::Correction = NoCorrection(),
    lateral = nothing,
)
    is_fixed(topo) && throw(ArgumentError("a fixed mesh has no node velocities to close"))
    length(topo.pinned) == 1 || throw(
        ArgumentError(
            "the balance rows need exactly one pinned node; the mesh pins " *
            "$(length(topo.pinned))",
        ),
    )
    p = only(topo.pinned)
    n = n_cells(topo)
    node_cells = [Int[] for _ = 1:n_nodes(topo)]
    for c = 1:n, nd in topo.cells[c].nodes
        push!(node_cells[nd], c)
    end
    length(node_cells[p]) == 1 ||
        throw(ArgumentError("the pinned node must lie at an end of the line"))
    row_node = zeros(Int, n)
    order = Int[]
    node, previous = p, 0
    while true
        next = filter(!=(previous), node_cells[node])
        isempty(next) && break
        c = only(next)
        nodes = topo.cells[c].nodes
        other = nodes[1] == node ? nodes[2] : nodes[1]
        row_node[c] = other
        push!(order, c)
        node, previous = other, c
    end
    length(order) == n ||
        throw(ArgumentError("the sweep from the pinned node does not reach every cell"))
    return BalanceRowClosure([p], row_node, order, Tuple{Int,Int}[], correction, lateral)
end

"""
	BalanceRowClosure(topology::MeshTopology{D}; correction = NoCorrection())    (D ≥ 2)

The balance-row closure of a single column of cells along the last axis,
swept from the pinned nodes: each cell must share a whole face with the
nodes already reached, and adds the nodes of the opposite face. The
correction is [`NoCorrection`](@ref), or a [`PrescribedCorrection`](@ref),
`w(x::SVector, t) -> SVector`, which prescribes every node's velocity on any
mesh.
"""
function BalanceRowClosure(
    topo::MeshTopology{D};
    correction::Correction = NoCorrection(),
    lateral = nothing,
) where {D}
    is_fixed(topo) && throw(ArgumentError("a fixed mesh has no node velocities to close"))
    correction isa Union{NoCorrection,PrescribedCorrection} || throw(
        ArgumentError("a column of more than one dimension takes no layer correction"),
    )
    lateral === nothing ||
        throw(ArgumentError("the variable cross section is a line's model"))
    n = n_cells(topo)
    correction isa PrescribedCorrection && return BalanceRowClosure(
        sort(topo.pinned),
        zeros(Int, n),
        Int[],
        Tuple{Int,Int}[],
        correction,
        lateral,
    )
    isempty(topo.pinned) && throw(ArgumentError("the column needs pinned nodes"))
    K = length(topo.cells[1].nodes)
    reached = falses(n_nodes(topo))
    reached[topo.pinned] .= true
    visited = falses(n)
    row_node = zeros(Int, n)
    order = Int[]
    ties = Tuple{Int,Int}[]
    while length(order) < n
        next = [
            i for i = 1:n if
            !visited[i] && count(nd -> reached[nd], topo.cells[i].nodes) == K ÷ 2
        ]
        length(next) == 1 || throw(
            ArgumentError(
                "the balance rows on more than one dimension need a single column of cells",
            ),
        )
        i = only(next)
        new = [nd for nd in topo.cells[i].nodes if !reached[nd]]
        visited[i] = true
        reached[new] .= true
        row_node[i] = new[1]
        append!(ties, ((nd, new[1]) for nd in new[2:end]))
        push!(order, i)
    end
    return BalanceRowClosure(sort(topo.pinned), row_node, order, ties, correction, lateral)
end

# The sweep coefficient σ A_f n_f of node `nd` in the balance of cell `k`: the
# area of the cell's face at that node, signed by its outward normal.
@inline function _coefficient(g::MeshGeometry, topo::MeshTopology{1}, k::Int, nd::Int)
    for f in topo.cells[k].faces
        topo.faces[f].nodes[1] == nd || continue
        σ = topo.faces[f].cells[1] == k ? 1 : -1
        return σ * g.face_area[f] * g.face_normal[f][1]
    end
    return zero(eltype(g.face_area))
end

# Σ_f A_f (w·n)_f over the faces of cell `k`: the rate of change of its volume.
@inline function _sweep_rate(g::MeshGeometry, topo::MeshTopology{1}, k::Int, w)
    s = zero(promote_type(eltype(g.face_area), eltype(w)))
    for f in topo.cells[k].faces
        σ = topo.faces[f].cells[1] == k ? 1 : -1
        s += σ * g.face_area[f] * g.face_normal[f][1] * w[topo.faces[f].nodes[1]]
    end
    return s
end

# Σ_f A_f ((w − c)·n)_f over the faces of cell `k`: the rate of change of its
# volume at the solid velocity.
@inline function _solid_sweep_rate(g::MeshGeometry, topo::MeshTopology{1}, k::Int, w, c)
    s = zero(promote_type(eltype(g.face_area), eltype(w), eltype(c)))
    for f in topo.cells[k].faces
        σ = topo.faces[f].cells[1] == k ? 1 : -1
        nd = topo.faces[f].nodes[1]
        s += σ * g.face_area[f] * g.face_normal[f][1] * (w[nd] - c[nd])
    end
    return s
end

"""
	correction!(ds, closure, setup)

Fill the node corrections `ds.c` at the derived state `ds`, once per residual
evaluation and before the kernels read them: nothing on a fixed mesh or
without a correction (`ds.c` stays zero), the node velocities under a
[`PrescribedCorrection`](@ref), and the layers' relaxed equidistribution
under a [`LayerCorrection`](@ref).
"""
@inline correction!(ds, ::Nothing, setup) = nothing
@inline correction!(ds, cl::BalanceRowClosure, setup) =
    _correction!(ds.c, cl.correction, ds, setup.reference)

# The velocity component `k` of node `nd` in a flat vector of `D` components per node.
@inline _component(w, ::Val{D}, nd::Int, k::Int) where {D} = w[D*(nd-1)+k]
@inline _node_vector(w, ::Val{D}, nd::Int) where {D} =
    SVector{D}(ntuple(k -> w[D*(nd-1)+k], Val(D)))

# Σ_f σ Σ_q s_{q,f} · w_q over the faces of cell `k` of a mesh of more than one
# dimension: the rate of change of its volume.
@inline function _sweep_rate(g::MeshGeometry, topo::MeshTopology{D}, k::Int, w) where {D}
    s = zero(promote_type(eltype(eltype(g.face_sweep)), eltype(w)))
    for f in topo.cells[k].faces
        σ = topo.faces[f].cells[1] == k ? 1 : -1
        nodes = topo.faces[f].nodes
        for q in eachindex(nodes)
            s += σ * dot(g.face_sweep[q, f], _node_vector(w, Val(D), nodes[q]))
        end
    end
    return s
end

# The coefficient of node `nd`'s component `axis` in the balance of cell `k`.
@inline function _coefficient(
    g::MeshGeometry{D},
    topo::MeshTopology{D},
    k::Int,
    nd::Int,
    axis::Int,
) where {D}
    a = zero(eltype(g.face_area))
    for f in topo.cells[k].faces
        σ = topo.faces[f].cells[1] == k ? 1 : -1
        nodes = topo.faces[f].nodes
        for q in eachindex(nodes)
            nodes[q] == nd && (a += σ * g.face_sweep[q, f][axis])
        end
    end
    return a
end

# The nodes cell `k` adds to the column: its row node and the nodes tied to it.
_new_nodes(cl::BalanceRowClosure, k::Int) =
    [cl.row_node[k]; [nd for (nd, leader) in cl.ties if leader == cl.row_node[k]]]

"""
	closure_rows!(c, closure, ds, topology, material, t)

Fill, from the kernels' reaction mass sources in `c.S_mass`, the cell strain
rates `c.θ` and, under a variable cross section, the history sources
`c.S_carried`; then the residuals `c.r_velocity` of the velocity rows. Runs
after the kernels and before any external source, so the strain rates are the
reactions' alone. A fixed mesh (`closure = nothing`) has none of these.
"""
@inline closure_rows!(c, ::Nothing, ds, topo, material, t) = nothing

function closure_rows!(
    c::Contributions{NC},
    cl::BalanceRowClosure,
    ds::DerivedState{NC},
    topo::MeshTopology,
    material::Material{NC},
    t,
) where {NC}
    @inbounds for i = 1:n_cells(topo)
        V = ds.V_eval[i]
        dξ = ntuple(j -> c.S_mass[j, i] / V, Val(NC))
        c.θ[i] = volume_change_rate(dξ, ds.T_trial[i], material)
    end
    _history_sources!(c, cl.lateral, topo)
    _velocity_rows!(c, cl, cl.correction, ds, topo, t)
    return nothing
end

@inline function _velocity_rows!(
    c,
    cl::BalanceRowClosure,
    ::NoCorrection,
    ds,
    topo::MeshTopology{1},
    t,
)
    g = ds.geometry
    θ_A = ds.θ_A[]
    for p in cl.pinned
        c.r_velocity[p] = ds.w[p]
    end
    @inbounds for k in cl.order
        c.r_velocity[cl.row_node[k]] =
            _sweep_rate(g, topo, k, ds.w) - (c.θ[k] - θ_A) * ds.V_eval[k]
    end
    return nothing
end

@inline function _velocity_rows!(c, cl::BalanceRowClosure, ::LayerCorrection, ds, topo, t)
    g = ds.geometry
    θ_A = ds.θ_A[]
    for p in cl.pinned
        c.r_velocity[p] = ds.w[p]
    end
    @inbounds for k in cl.order
        c.r_velocity[cl.row_node[k]] =
            _solid_sweep_rate(g, topo, k, ds.w, ds.c) - (c.θ[k] - θ_A) * ds.V_eval[k]
    end
    return nothing
end

@inline function _velocity_rows!(
    c,
    cl::BalanceRowClosure,
    corr::PrescribedCorrection,
    ds,
    topo::MeshTopology{1},
    t,
)
    @inbounds for nd in eachindex(ds.w)
        c.r_velocity[nd] = ds.w[nd] - corr.w(ds.x[nd], t)
    end
    return nothing
end

# A column of more than one dimension: the lateral components and the pinned
# nodes at rest, the ties, and the balances on the column axis.
function _velocity_rows!(
    c,
    cl::BalanceRowClosure,
    ::NoCorrection,
    ds,
    topo::MeshTopology{D},
    t,
) where {D}
    g = ds.geometry
    θ_A = ds.θ_A[]
    @inbounds for nd = 1:n_nodes(topo), k = 1:(D-1)
        c.r_velocity[D*(nd-1)+k] = ds.w[D*(nd-1)+k]
    end
    @inbounds for p in cl.pinned
        c.r_velocity[D*p] = ds.w[D*p]
    end
    @inbounds for (nd, leader) in cl.ties
        c.r_velocity[D*nd] = ds.w[D*nd] - ds.w[D*leader]
    end
    @inbounds for k in cl.order
        c.r_velocity[D*cl.row_node[k]] =
            _sweep_rate(g, topo, k, ds.w) - (c.θ[k] - θ_A) * ds.V_eval[k]
    end
    return nothing
end

function _velocity_rows!(
    c,
    cl::BalanceRowClosure,
    corr::PrescribedCorrection,
    ds,
    topo::MeshTopology{D},
    t,
) where {D}
    @inbounds for nd = 1:n_nodes(topo)
        w = corr.w(_node_vector(ds.x, Val(D), nd), t)
        for k = 1:D
            c.r_velocity[D*(nd-1)+k] = ds.w[D*(nd-1)+k] - w[k]
        end
    end
    return nothing
end

"""
	consistent_velocities!(w, closure, c, ds, topology, t) -> w

The node velocities that satisfy every velocity row at the derived state
`ds`, given the strain rates in `c.θ` (from [`closure_rows!`](@ref)), the
lateral rate in `ds.θ_A` and the node corrections in `ds.c`: one sweep of the
bidiagonal balance rows from the pinned node (the corrections read positions
and masses only, so one pass is exact), zero on the row node of a dead cell,
or the prescribed velocity at every node.
"""
function consistent_velocities!(
    w::AbstractVector,
    cl::BalanceRowClosure,
    c::Contributions,
    ds::DerivedState,
    topo::MeshTopology,
    t,
)
    _consistent_velocities!(w, cl, cl.correction, c, ds, topo, t)
    return w
end

function _consistent_velocities!(
    w,
    cl::BalanceRowClosure,
    ::NoCorrection,
    c,
    ds,
    topo::MeshTopology{1},
    t,
)
    g = ds.geometry
    θ_A = ds.θ_A[]
    w[cl.pinned] .= 0
    for k in cl.order
        row = cl.row_node[k]
        if !ds.active[k]
            w[row] = 0
            continue
        end
        nodes = topo.cells[k].nodes
        other = nodes[1] == row ? nodes[2] : nodes[1]
        rate = (c.θ[k] - θ_A) * ds.V_eval[k] - _coefficient(g, topo, k, other) * w[other]
        w[row] = rate / _coefficient(g, topo, k, row)
    end
    return w
end

function _consistent_velocities!(
    w,
    cl::BalanceRowClosure,
    ::LayerCorrection,
    c,
    ds,
    topo,
    t,
)
    g = ds.geometry
    θ_A = ds.θ_A[]
    w[cl.pinned] .= 0
    for k in cl.order
        row = cl.row_node[k]
        if !ds.active[k]
            w[row] = 0
            continue
        end
        nodes = topo.cells[k].nodes
        other = nodes[1] == row ? nodes[2] : nodes[1]
        rate =
            (c.θ[k] - θ_A) * ds.V_eval[k] -
            _coefficient(g, topo, k, other) * (w[other] - ds.c[other])
        w[row] = ds.c[row] + rate / _coefficient(g, topo, k, row)
    end
    return w
end

function _consistent_velocities!(
    w,
    cl::BalanceRowClosure,
    corr::PrescribedCorrection,
    c,
    ds,
    topo::MeshTopology{1},
    t,
)
    for nd in eachindex(w)
        w[nd] = corr.w(ds.x[nd], t)
    end
    return w
end

function _consistent_velocities!(
    w,
    cl::BalanceRowClosure,
    corr::PrescribedCorrection,
    c,
    ds,
    topo::MeshTopology{D},
    t,
) where {D}
    for nd = 1:n_nodes(topo)
        v = corr.w(_node_vector(ds.x, Val(D), nd), t)
        for k = 1:D
            w[D*(nd-1)+k] = v[k]
        end
    end
    return w
end

# A column of more than one dimension: the lateral components and the pinned
# nodes at rest; each cell's balance then fixes the common axial velocity of
# the nodes it adds, from those already reached.
function _consistent_velocities!(
    w,
    cl::BalanceRowClosure,
    ::NoCorrection,
    c,
    ds,
    topo::MeshTopology{D},
    t,
) where {D}
    g = ds.geometry
    θ_A = ds.θ_A[]
    fill!(w, 0)
    for k in cl.order
        new = _new_nodes(cl, k)
        ds.active[k] || continue
        known = zero(eltype(w))
        for nd in topo.cells[k].nodes
            nd in new || (known += _coefficient(g, topo, k, nd, D) * w[D*nd])
        end
        lead = sum(nd -> _coefficient(g, topo, k, nd, D), new)
        v = ((c.θ[k] - θ_A) * ds.V_eval[k] - known) / lead
        for nd in new
            w[D*nd] = v
        end
    end
    return w
end

"""
	rank_check(closure, topology, geometry) -> Bool

Whether the velocity rows determine every node velocity at the metrics
`geometry`: the balance rows do when every cell's coefficient on its row node
(the area of its face there) is nonzero; prescribed velocities always do.
"""
function rank_check(cl::BalanceRowClosure, topo::MeshTopology{1}, g::MeshGeometry)
    cl.correction isa PrescribedCorrection && return true
    return all(k -> _coefficient(g, topo, k, cl.row_node[k]) != 0, cl.order)
end

function rank_check(cl::BalanceRowClosure, topo::MeshTopology{D}, g::MeshGeometry) where {D}
    cl.correction isa PrescribedCorrection && return true
    return all(
        k -> sum(nd -> _coefficient(g, topo, k, nd, D), _new_nodes(cl, k)) != 0,
        cl.order,
    )
end
