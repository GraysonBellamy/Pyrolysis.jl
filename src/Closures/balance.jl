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
	BalanceRowClosure{O,X} <: Closure

The node velocities of a 1D mesh from the discrete volume balance of every
cell. The pinned node is at rest, `0 = w_p`; cell `k`'s balance sits on the
velocity row of its node farther from the pinned node,

	0 = Σ_f A_f (w·n)_f − (θ_k − θ_A) V_k,

in 1D `A_{k+1} w_{k+1} − A_k w_k − (θ_k − θ_A) V_k`, with `θ_k` the strain
rate of the cell's reaction mass sources (`Physics.volume_change_rate` at the
trial temperature) and `θ_A` the lateral strain rate of a variable cross
section. The rows are bidiagonal from the pinned node and each reads one
cell and its two nodes, so the Jacobian stays banded where the swept
velocity would couple every node to every cell between it and the pinned
node. Since the cell volumes are computed from the same node velocities'
positions, `dVₖ/dt = θₖ Vₖ` holds semi-discretely.

`override`, a function `w(x, t)`, replaces every velocity row by
`0 = w_n − w(x_n, t)`: a prescribed mesh velocity, for verification, under
which the mesh moves over a solid at rest. `lateral` is the variable cross
section's [`LateralShrinkage`](@ref), or `nothing`.

`row_node[k]` is the node whose velocity row carries cell `k`'s balance and
`order` the cells from the pinned node outward.
"""
struct BalanceRowClosure{O,X} <: Closure
    pinned::Int
    row_node::Vector{Int}
    order::Vector{Int}
    override::O
    lateral::X
end

"""
	BalanceRowClosure(topology; override = nothing, lateral = nothing)

The balance-row closure of a 1D mesh with one pinned node at an end of the
line.
"""
function BalanceRowClosure(topo::MeshTopology{1}; override = nothing, lateral = nothing)
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
    return BalanceRowClosure(p, row_node, order, override, lateral)
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
    topo::MeshTopology{1},
    material::Material{NC},
    t,
) where {NC}
    @inbounds for i = 1:n_cells(topo)
        V = ds.V_eval[i]
        dξ = ntuple(j -> c.S_mass[j, i] / V, Val(NC))
        c.θ[i] = volume_change_rate(dξ, ds.T_trial[i], material)
    end
    _history_sources!(c, cl.lateral, topo)
    _velocity_rows!(c, cl, cl.override, ds, topo, t)
    return nothing
end

@inline function _velocity_rows!(c, cl::BalanceRowClosure, ::Nothing, ds, topo, t)
    g = ds.geometry
    θ_A = ds.θ_A[]
    c.r_velocity[cl.pinned] = ds.w[cl.pinned]
    @inbounds for k in cl.order
        c.r_velocity[cl.row_node[k]] =
            _sweep_rate(g, topo, k, ds.w) - (c.θ[k] - θ_A) * ds.V_eval[k]
    end
    return nothing
end

@inline function _velocity_rows!(c, cl::BalanceRowClosure, override, ds, topo, t)
    @inbounds for nd in eachindex(ds.w)
        c.r_velocity[nd] = ds.w[nd] - override(ds.x[nd], t)
    end
    return nothing
end

"""
	consistent_velocities!(w, closure, c, ds, topology, t) -> w

The node velocities that satisfy every velocity row at the derived state
`ds`, given the strain rates in `c.θ` (from [`closure_rows!`](@ref)) and the
lateral rate in `ds.θ_A`: one sweep of the bidiagonal balance rows from the
pinned node, or the prescribed velocity at every node.
"""
function consistent_velocities!(
    w::AbstractVector,
    cl::BalanceRowClosure,
    c::Contributions,
    ds::DerivedState,
    topo::MeshTopology{1},
    t,
)
    _consistent_velocities!(w, cl, cl.override, c, ds, topo, t)
    return w
end

function _consistent_velocities!(w, cl::BalanceRowClosure, ::Nothing, c, ds, topo, t)
    g = ds.geometry
    θ_A = ds.θ_A[]
    w[cl.pinned] = 0
    for k in cl.order
        row = cl.row_node[k]
        nodes = topo.cells[k].nodes
        other = nodes[1] == row ? nodes[2] : nodes[1]
        rate = (c.θ[k] - θ_A) * ds.V_eval[k] - _coefficient(g, topo, k, other) * w[other]
        w[row] = rate / _coefficient(g, topo, k, row)
    end
    return w
end

function _consistent_velocities!(w, cl::BalanceRowClosure, override, c, ds, topo, t)
    for nd in eachindex(w)
        w[nd] = override(ds.x[nd], t)
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
    cl.override === nothing || return true
    return all(k -> _coefficient(g, topo, k, cl.row_node[k]) != 0, cl.order)
end
