# ============================================================================
# Potential closure
# ============================================================================

"""
	PotentialClosure{C} <: Closure

The node velocities of a moving mesh of any dimension from one multiplier `φᵢ`
per cell: the velocity of minimum discrete `L²` norm, `Σ_n W_n |w_n|²`, that
satisfies the volume balance of every live cell,

	0 = Σ_f σ_f Σ_q s_{q,f} · w_q − (θᵢ − θ_A) Vᵢ,

with `s_{q,f}` the sweep vectors of the faces (`MeshGeometry.face_sweep`) and
`θᵢ` the cell's reaction strain rate. Its stationarity condition is
`w = W⁻¹ Dᵀ φ` with `D` the sweep operator over the free velocity
components and `W` the lumped node volumes of the live cells: each node's
velocity is a local function of the multipliers of the live cells around
it, written before the kernels run, and each cell's balance is the
algebraic row of its multiplier. Both rows are local; the global system
`D W⁻¹ Dᵀ φ = (θ − θ_A) V` that the rows imply is solved only for the
consistent initial state ([`consistent_multipliers`](@ref)), never inside
the residual. The weight is the lumped volume: two unit squares side by side
under a uniform expansion with their lateral motion suppressed move
uniformly, where the unweighted form (`W = 1`) gives the shared node twice
the velocity of the others.

A node component is free unless its node is pinned or it is the fixed
component of a slip wall (`free`, `D` flags per node). `corners[n]` lists
the `(cell, corner)` pairs of node `n` and `sweeps[n]` the `(cell, face, q,
σ)` terms of its sweep vectors; a dead cell contributes to neither.
"""
struct PotentialClosure{C<:Correction} <: Closure
    free::BitVector
    corners::Vector{Vector{Tuple{Int,Int}}}
    sweeps::Vector{Vector{NTuple{4,Int}}}
    weighted::Bool
    correction::C
end

"""
	PotentialClosure(topology; slip = (), weighted = true)

The potential closure of a mesh: the pinned nodes at rest, and on the faces
of each `tag => k` pair of `slip` the velocity component `k` zero (a wall the
mesh slides along, on an axis-aligned face). `weighted = false` takes `W = 1`,
the unweighted form, for comparison only.
"""
function PotentialClosure(topo::MeshTopology{D}; slip = (), weighted::Bool = true) where {D}
    is_fixed(topo) && throw(ArgumentError("a fixed mesh has no node velocities to close"))
    nn = n_nodes(topo)
    free = trues(D * nn)
    for p in topo.pinned, k = 1:D
        free[D*(p-1)+k] = false
    end
    for (tag, k) in slip
        1 <= k <= D || throw(ArgumentError("a slip wall fixes a component in 1:$D"))
        for f in boundary_faces(topo, tag), nd in topo.faces[f].nodes
            free[D*(nd-1)+k] = false
        end
    end
    corners = [Tuple{Int,Int}[] for _ = 1:nn]
    sweeps = [NTuple{4,Int}[] for _ = 1:nn]
    for i = 1:n_cells(topo)
        for (a, nd) in enumerate(topo.cells[i].nodes)
            push!(corners[nd], (i, a))
        end
        for f in topo.cells[i].faces
            σ = topo.faces[f].cells[1] == i ? 1 : -1
            for (q, nd) in enumerate(topo.faces[f].nodes)
                push!(sweeps[nd], (i, f, q, σ))
            end
        end
    end
    return PotentialClosure(free, corners, sweeps, weighted, NoCorrection())
end

# The lumped volume of node `nd` over the live cells.
@inline function _node_weight(cl::PotentialClosure, g::MeshGeometry, active, nd::Int)
    W = zero(eltype(g.cell_volume))
    cl.weighted || return one(W)
    @inbounds for (i, a) in cl.corners[nd]
        active[i] && (W += g.cell_node_volume[a, i])
    end
    return W
end

# (Dᵀ φ) at node `nd`: Σ over the live cells around it of φᵢ Σ_f σ s_{q,f}.
@inline function _node_gradient(
    cl::PotentialClosure,
    g::MeshGeometry{D},
    active,
    φ,
    nd::Int,
) where {D}
    G = zero(SVector{D,promote_type(eltype(eltype(g.face_sweep)), eltype(φ))})
    @inbounds for (i, f, q, σ) in cl.sweeps[nd]
        active[i] && (G += (σ * φ[i]) * g.face_sweep[q, f])
    end
    return G
end

@inline correction!(ds, cl::PotentialClosure, setup) = _potential_velocities!(ds.w, cl, ds)

# The node velocities w = W⁻¹ Dᵀ φ over the free components, zero elsewhere.
function _potential_velocities!(
    w,
    cl::PotentialClosure,
    ds::DerivedState{NC,T,<:MeshGeometry{D}},
) where {NC,T,D}
    g = ds.geometry
    @inbounds for nd in eachindex(cl.corners)
        W = _node_weight(cl, g, ds.active, nd)
        G = _node_gradient(cl, g, ds.active, ds.multiplier, nd)
        for k = 1:D
            j = D * (nd - 1) + k
            w[j] = cl.free[j] && W > 0 ? G[k] / W : zero(eltype(w))
        end
    end
    return w
end

function closure_rows!(
    c::Contributions{NC},
    cl::PotentialClosure,
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
    multiplier_rows!(c, cl, ds, topo)
    return nothing
end

"""
	multiplier_rows!(c, closure::PotentialClosure, ds, topology)

The residuals `c.r_multiplier` of the cells' balances at the node velocities
`ds.w` and the strain rates `c.θ`.
"""
function multiplier_rows!(
    c::Contributions,
    cl::PotentialClosure,
    ds::DerivedState,
    topo::MeshTopology,
)
    g = ds.geometry
    θ_A = ds.θ_A[]
    @inbounds for i = 1:n_cells(topo)
        c.r_multiplier[i] = _sweep_rate(g, topo, i, ds.w) - (c.θ[i] - θ_A) * ds.V_eval[i]
    end
    return nothing
end

# The sweep operator over the live cells and the free components of the nodes
# they hold, with those components' weights: D (live cells × components) and W.
function _operator(
    cl::PotentialClosure,
    topo::MeshTopology{D},
    g::MeshGeometry,
    active,
) where {D}
    live = findall(active)
    row = zeros(Int, n_cells(topo))
    row[live] .= 1:length(live)
    # The free components of the nodes some live cell holds: a node only dead
    # cells hold has no balance to enter and no weight.
    held = [any(((i, a),) -> active[i], cl.corners[nd]) for nd in eachindex(cl.corners)]
    cols = [j for j in findall(cl.free) if held[(j-1)÷D+1]]
    col = zeros(Int, length(cl.free))
    col[cols] .= 1:length(cols)
    I_, J_, V_ = Int[], Int[], Float64[]
    for nd in eachindex(cl.sweeps), (i, f, q, σ) in cl.sweeps[nd], k = 1:D
        j = D * (nd - 1) + k
        (active[i] && col[j] != 0) || continue
        push!(I_, row[i])
        push!(J_, col[j])
        push!(V_, σ * g.face_sweep[q, f][k])
    end
    Dop = sparse(I_, J_, V_, length(live), length(cols))
    W = [_node_weight(cl, g, active, (j - 1) ÷ D + 1) for j in cols]
    return Dop, W, live
end

"""
	rank_check(closure::PotentialClosure, topology, geometry; active = all) -> Bool

Whether the balances of the live cells determine the multipliers at the
metrics `geometry`: the sweep operator over the free components has full row
rank, so `D W⁻¹ Dᵀ` is symmetric positive definite. A body with no free
component on some live cell's faces, or fully pinned, fails.
"""
function rank_check(
    cl::PotentialClosure,
    topo::MeshTopology,
    g::MeshGeometry;
    active::AbstractVector{Bool} = trues(n_cells(topo)),
)
    Dop, W, live = _operator(cl, topo, g, active)
    isempty(live) && return true
    all(>(0), W) || return false
    return rank(Matrix(Dop)) == length(live)
end

"""
	consistent_multipliers(closure::PotentialClosure, c, ds, topology) -> φ

The multipliers whose velocities satisfy every live cell's balance at the
derived state `ds`, given the strain rates in `c.θ` and the lateral rate in
`ds.θ_A`: one sparse symmetric positive-definite solve of
`D W⁻¹ Dᵀ φ = (θ − θ_A) V`, outside the residual, refused where the system
is singular (see [`rank_check`](@ref)). Dead cells get zero.
"""
function consistent_multipliers(
    cl::PotentialClosure,
    c::Contributions,
    ds::DerivedState,
    topo::MeshTopology,
)
    g = ds.geometry
    Dop, W, live = _operator(cl, topo, g, ds.active)
    b = [(c.θ[i] - ds.θ_A[]) * ds.V_eval[i] for i in live]
    A = Symmetric(Dop * Diagonal(1 ./ W) * Dop')
    F = cholesky(A; check = false)
    refused =
        ArgumentError("the potential closure's balances do not determine the multipliers")
    issuccess(F) || throw(refused)
    φ_live = F \ b
    # A singular system can factor through round-off pivots: the solve must
    # also satisfy it.
    norm(A * φ_live - b) <= 1e-10 * norm(b) || throw(refused)
    φ = zeros(n_cells(topo))
    φ[live] = φ_live
    return φ
end
