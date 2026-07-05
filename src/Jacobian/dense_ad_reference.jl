"""
`DenseAD_ForwardDiff` backend.

Forward-mode automatic differentiation of the canonical `residual!` against
the entire state vector, written into a dense `Matrix{Float64}`.

This is the **`verify_jacobian` reference** — every other backend's output is
compared against `DenseAD_ForwardDiff` at canonical sample points (Step
B12/B13). It carries no sparsity assumption and no hand-coded math, so a
discrepancy between this and any other backend is unambiguously a bug in
the other backend.

Built on `Pyrolysis.Problem.build_ad_workspace`, which is the canonical
Dual-typed workspace constructor. The Dual eltype is captured from the
`ForwardDiff.JacobianConfig`'s seed at construction time so the workspace
type is fully concrete.

For ALE problems this backend differentiates the full state vector
(including z and χ) in one pass.
"""

using ForwardDiff

# ============================================================================
# Backend tag + cache
# ============================================================================

"""
	DenseAD_ForwardDiff <: JacobianBackend

Tag for the dense ForwardDiff backend. Carries no state.
"""
struct DenseAD_ForwardDiff <: JacobianBackend end

"""
	DenseADCache{ChunkSize, Tu, WSAD, CFG}

Per-problem cache. Holds the Dual-typed workspace, the `JacobianConfig`,
and a Float64 RHS scratch (`f_cache`) that ForwardDiff uses internally.

# Type parameters
- `ChunkSize`: ForwardDiff chunk size, capped at the state-vector width.
- `Tu`: `ForwardDiff.Dual{Tag, Float64, ChunkSize}` carried by `ws_ad`.
- `WSAD`: full `Workspace{...}` type with `Tu` as the scalar eltype.
- `CFG`: concrete `ForwardDiff.JacobianConfig` type (parameterized on
  the chunk and the closure tag).

# Fields
- `ws_ad`: Dual-typed workspace built via `build_ad_workspace(def, ws, Tu)`
- `cfg`: ForwardDiff JacobianConfig (chunk + tag + Dual seeds)
- `f_cache`: full-length Float64 RHS scratch buffer
"""
struct DenseADCache{ChunkSize,Tu,WSAD,CFG}
    ws_ad::WSAD
    cfg::CFG
    f_cache::Vector{Float64}
end

# ============================================================================
# Construction
# ============================================================================

function _build_jacobian_spec(
    ::DenseAD_ForwardDiff,
    def::ProblemDef,
    ws::Workspace;
    chunk_size::Int = 8,
    kwargs...,
)
    chunk_size >= 1 || throw(ArgumentError("chunk_size must be positive"))

    n = state_length(ws.layout)
    actual_chunk_size = min(chunk_size, n)
    x_template = zeros(Float64, n)
    f_cache = zeros(Float64, n)
    cfg = ForwardDiff.JacobianConfig(
        nothing,
        f_cache,
        x_template,
        ForwardDiff.Chunk{actual_chunk_size}(),
    )
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)
    cache =
        DenseADCache{actual_chunk_size,Tu,typeof(ws_ad),typeof(cfg)}(ws_ad, cfg, f_cache)
    return JacobianSpec(DenseAD_ForwardDiff(), nothing, nothing, cache, Int[], Int[])
end

# ============================================================================
# Compute
# ============================================================================

function _compute_jacobian!(
    J::AbstractMatrix{Float64},
    ::DenseAD_ForwardDiff,
    cache::DenseADCache,
    sparsity_phys,
    sparsity_geom,
    u,
    def,
    t,
    ws,
)
    # The closure dispatches `residual!` on the Dual-typed workspace so leaf
    # kernels see Duals end-to-end. `let` captures concrete bindings.
    f! = let def = def, ws_ad = cache.ws_ad, t = t
        (du, u_inner) -> residual!(du, u_inner, def, t, ws_ad)
    end
    ForwardDiff.jacobian!(J, f!, cache.f_cache, u, cache.cfg)
    return J
end
