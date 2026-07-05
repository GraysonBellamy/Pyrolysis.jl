"""
The `GeometryOperator` — the ALE z-block AD pass exposed through
the operator interface.

Owns the geometry-block (z-column) Jacobian for ALE problems: ∂(physics
rows)/∂z and ∂(z rows)/∂z. ForwardDiff over the **z subvector only**: the
closure substitutes the perturbed z values into the full state vector
(T, ξ, χ slots come from the Float64 baseline `u`, untouched), then runs
the canonical pure `residual!` against a Dual-typed workspace. The result
is the geometry block: the n_total × n_nodes columns of `∂(residual)/∂z`.

Two-pass discipline: this operator runs **independently** of the cell-local
AD pass. Each pass writes Float64 entries into the shared `J::AbstractMatrix`;
the only Dual nesting is each pass's own Tag, so `eltype(J) == Float64`
after both passes.

Construction reuses [`Pyrolysis.Problem.build_ad_workspace`](@ref) — the
canonical Dual-workspace builder. Backends never reinvent it.

# Source

The math here lives alongside the operator interface. The `Structured`
backend dispatches via `contribute_to_jacobian!(GeometryOperator(), …)`.
"""

using ForwardDiff
using SparseArrays: SparseMatrixCSC
using ..Problem: ProblemDef, Workspace, StateLayout, state_length, build_ad_workspace
using ..Residual: residual!
using ..Materials: n_components

"""
	GeometryOperator <: Operator

Stateless tag. Per-problem cache lives on the backend's `geom_cache` field
(a [`GeometryADCache`](@ref)) and is built once at `JacobianSpec`
construction.

Note: this operator's stencil is global (z-column → all rows) by virtue of
running the full residual under AD; the conventional `(cells = …,)`
declaration doesn't apply. The default `stencil(::Operator) = (cells = 0:0,)`
fallback is left in place — sparsity for the z-columns is built directly
in `build_sparsity` from the geometry block's structural pattern, not from
this operator's `stencil` declaration.
"""
struct GeometryOperator <: Operator end

reads(::GeometryOperator) = ()
writes(::GeometryOperator) = ()

# ============================================================================
# Cache
# ============================================================================

"""
	GeometryADCache{ChunkSize, Tu, WSAD, CFG}

Type-stable per-problem cache for the [`GeometryOperator`]'s
`contribute_to_jacobian!`. Built once at `JacobianSpec` construction;
reused across every Jacobian-assembly call.

# Type parameters
- `ChunkSize`: ForwardDiff chunk size, capped at the z-block width.
- `Tu`: the `ForwardDiff.Dual{Tag, Float64, ChunkSize}` carried by `ws_ad`.
- `WSAD`: full `Workspace{...}` type with `Tu` as the scalar eltype.
- `CFG`: `ForwardDiff.JacobianConfig` parameterized to match the n_nodes
  input width and the chunk.

# Fields
- `ws_ad`: Dual-typed workspace built via `build_ad_workspace(def, ws, Tu)`
- `cfg`: ForwardDiff config (chunk + tag + Dual seeds)
- `u_dual`: full-length state vector with `Tu` eltype, refilled from the
  Float64 `u` on every call (T, ξ, χ slots) so the AD closure only needs
  to overwrite the z slots
- `f_cache`: Float64 vector of length `n_total` — the `y` argument to
  `ForwardDiff.jacobian!`. ForwardDiff uses it to extract Float64 output
  values from the Dual residual; the Dual du buffer itself is managed
  inside `cfg`.
- `z_seed`: Float64 vector of length `n_nodes` — the linearization point
  (the current value of `u[z_range]`) handed to `ForwardDiff.jacobian!`
- `J_dense`: scratch `n_total × n_nodes` `Matrix{Float64}` for the
  contiguous Jacobian block ForwardDiff writes into. The orchestrator
  scatters from this dense buffer into the user's `J` (sparse or dense).
- `z_off`: `n_cells * (1 + NC)` — base offset into `u` where the z block
  starts.
"""
struct GeometryADCache{ChunkSize,Tu,WSAD,CFG}
    ws_ad::WSAD
    cfg::CFG
    u_dual::Vector{Tu}
    f_cache::Vector{Float64}
    z_seed::Vector{Float64}
    J_dense::Matrix{Float64}
    z_off::Int
    n_total::Int
    n_nodes::Int
end

"""
	GeometryADCache(def::ProblemDef, ws::Workspace; chunk_size::Int = 8)

Build a [`GeometryADCache`](@ref) for the given problem. Allocates the
Dual-typed workspace, the `JacobianConfig`, and the dense scratch matrix
**eagerly** — call this once at `JacobianSpec` construction, never inside
the integrator hot path.

The effective chunk size is type-encoded (`ForwardDiff.Chunk{N}()`); the
returned cache's type captures it as the `ChunkSize` parameter so leaf
kernels see a fully concrete type.
"""
function GeometryADCache(def::ProblemDef, ws::Workspace; chunk_size::Int = 8)
    chunk_size >= 1 || throw(ArgumentError("chunk_size must be positive"))

    layout = ws.layout
    NC = n_components(layout)
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes
    n_total = state_length(layout)
    z_off = n_cells * (1 + NC)
    actual_chunk_size = min(chunk_size, n_nodes)

    # Build the JacobianConfig against the (du::Vector{Float64}, z::Vector{Float64})
    # pair to extract the chunk's Dual type without manually constructing a Tag.
    z_seed = zeros(Float64, n_nodes)
    f_cache_float = zeros(Float64, n_total)
    cfg = ForwardDiff.JacobianConfig(
        nothing,
        f_cache_float,
        z_seed,
        ForwardDiff.Chunk{actual_chunk_size}(),
    )
    Tu = eltype(cfg.duals[1])

    ws_ad = build_ad_workspace(def, ws, Tu)
    u_dual = zeros(Tu, n_total)
    J_dense = zeros(Float64, n_total, n_nodes)

    return GeometryADCache{actual_chunk_size,Tu,typeof(ws_ad),typeof(cfg)}(
        ws_ad,
        cfg,
        u_dual,
        f_cache_float,
        z_seed,
        J_dense,
        z_off,
        n_total,
        n_nodes,
    )
end

# ============================================================================
# Operator-pattern dispatcher
# ============================================================================

"""
	contribute_to_jacobian!(jac, ::GeometryOperator, u, ctx, t,
	                        cache::GeometryADCache, def::ProblemDef) -> nothing

Operator-pattern entry point used by the `Structured` backend. Thin wrapper
that delegates to [`apply_geometry_jacobian!`](@ref). Kept as a separate
method so the backend's call site can be a single uniform
`contribute_to_jacobian!` walk.
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::GeometryOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    cache::GeometryADCache,
    def::ProblemDef,
)
    apply_geometry_jacobian!(jac.A, cache, u, def, Float64(t))
    return nothing
end

# ============================================================================
# Assembler — public API consumed by the operator dispatcher above. The
# function remains separately callable so direct geometry-block assembly
# (e.g. for tests) can reuse the same AD plumbing.
# ============================================================================

"""
	apply_geometry_jacobian!(J, cache::GeometryADCache, u, def, t)

Assemble the geometry block of `J` (every (row, col) where col is a z
index) by ForwardDiff'ing the residual against the z subvector of `u`.

`J` may be any `AbstractMatrix{Float64}` — `Matrix` for the equivalence
tests, `SparseMatrixCSC{Float64,Int}` for the production hot path. The
function writes each `J[row, z_off + col]` directly; non-z columns of `J`
are **not** touched (the cell-local-AD or analytical kernel passes are
responsible for them).

For dense `J`, every entry in the (n_total × n_nodes) z-column block is
overwritten. For `SparseMatrixCSC`, only the structural entries already
present in the column are written — non-structural positions are
guaranteed zero by the sparsity contract (see [`build_sparsity`](@ref)),
so dropping them is correct.

Allocates only on first call (the `JacobianConfig`'s Dual seeds are inside
`cache`, not on the stack). Subsequent calls reuse `cache.u_dual`,
`cache.f_cache`, `cache.J_dense` in-place.
"""
function apply_geometry_jacobian!(
    J::AbstractMatrix{Float64},
    cache::GeometryADCache,
    u::AbstractVector{Float64},
    def::ProblemDef,
    t::Float64,
)
    _fill_geometry_jacobian_dense!(cache, u, def, t)

    J_dense = cache.J_dense
    z_off = cache.z_off
    n_total = cache.n_total
    n_nodes = cache.n_nodes

    size(J, 1) >= n_total && size(J, 2) >= z_off + n_nodes ||
        throw(ArgumentError("apply_geometry_jacobian!: J too small for layout"))

    # Dense scatter: copy every row of the contiguous block into J's z columns.
    @inbounds for col = 1:n_nodes
        J_col = z_off + col
        for row = 1:n_total
            J[row, J_col] = J_dense[row, col]
        end
    end

    return J
end

function apply_geometry_jacobian!(
    J::SparseMatrixCSC{Float64,Int},
    cache::GeometryADCache,
    u::AbstractVector{Float64},
    def::ProblemDef,
    t::Float64,
)
    _fill_geometry_jacobian_dense!(cache, u, def, t)

    J_dense = cache.J_dense
    z_off = cache.z_off
    n_total = cache.n_total
    n_nodes = cache.n_nodes

    size(J, 1) >= n_total && size(J, 2) >= z_off + n_nodes ||
        throw(ArgumentError("apply_geometry_jacobian!: J too small for layout"))

    # Sparse scatter — column-walk over the existing colptr/rowval, writing
    # only the structural entries. Non-structural positions in the dense
    # block are guaranteed zero by the sparsity contract; silently dropping
    # them is correct and `O(nnz)` rather than `O(n_total × n_nodes)`.
    colptr = J.colptr
    rowval = J.rowval
    nzval = J.nzval
    @inbounds for col = 1:n_nodes
        J_col = z_off + col
        r1 = colptr[J_col]
        r2 = colptr[J_col+1] - 1
        for idx = r1:r2
            row = rowval[idx]
            nzval[idx] = J_dense[row, col]
        end
    end

    return J
end

# ----------------------------------------------------------------------------
# Internal: run ForwardDiff on the residual w.r.t. the z subvector and leave
# the result in `cache.J_dense`. Shared by every `apply_geometry_jacobian!`
# scatter method above.
# ----------------------------------------------------------------------------
function _fill_geometry_jacobian_dense!(
    cache::GeometryADCache{ChunkSize,Tu},
    u::AbstractVector{Float64},
    def::ProblemDef,
    t::Float64,
) where {ChunkSize,Tu}
    ws_ad = cache.ws_ad
    u_dual = cache.u_dual
    f_cache = cache.f_cache
    z_seed = cache.z_seed
    J_dense = cache.J_dense
    z_off = cache.z_off
    n_total = cache.n_total
    n_nodes = cache.n_nodes

    length(u) == n_total || throw(
        ArgumentError(
            "apply_geometry_jacobian!: u has length $(length(u)), expected $n_total",
        ),
    )

    # Initialize u_dual with the current Float64 baseline (zero partials in
    # all slots). The closure overwrites the z slots with Dual-perturbed
    # values on each evaluation; the T, ξ, χ slots stay constant.
    @inbounds for i = 1:n_total
        u_dual[i] = Tu(u[i])
    end
    @inbounds for k = 1:n_nodes
        z_seed[k] = u[z_off+k]
    end

    # Closure: perturb z slots, then run the canonical pure residual on the
    # Dual-typed workspace. The `let` block makes the captures concrete so
    # Julia compiles a specialized method per call site.
    f_geom! =
        let u_dual = u_dual,
            z_off = z_off,
            n_nodes = n_nodes,
            def = def,
            t = t,
            ws_ad = ws_ad

            (du, z_in) -> begin
                @inbounds for k = 1:n_nodes
                    u_dual[z_off+k] = z_in[k]
                end
                residual!(du, u_dual, def, t, ws_ad)
                return du
            end
        end

    fill!(J_dense, 0.0)
    ForwardDiff.jacobian!(J_dense, f_geom!, f_cache, z_seed, cache.cfg)
    return J_dense
end
