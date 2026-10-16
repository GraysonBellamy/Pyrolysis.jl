"""
Workspace — mutable container that owns every per-step cache and work
array used by the unified residual.

The `Workspace` is the *only* place where mutable per-step state lives.
Material, boundary conditions, and the radiation model are read from
`def::ProblemDef` (the immutable side); the workspace does not duplicate
them.

Every cache is parameterized so that `Nothing` can appear as a **concrete
type parameter** rather than as a `Union{Nothing, T}`. Trait dispatch on
the type — not runtime `isnothing` checks — decides whether to touch a
cache.

# Type Parameters
- `NC`: number of components
- `Mode`: bundled `SimulationMode` (matches the parent `ProblemDef`)
- `Tu`: scalar eltype carried by the live state (`Float64` for the
        canonical workspace; `ForwardDiff.Dual{Tag, Float64, N}` for the
        AD workspace allocated by an AD Jacobian backend)
- `MeshT`:    Reference to the mutable `Mesh1D` carrying connectivity
              and per-step geometry. The mesh is **read-only** during
              residual / Jacobian assembly; writes happen only between
              integrator-accepted steps via `StepUpdateCallback`.
              Kernels only read connectivity (cell IDs, face IDs,
              neighbor lookups, boundary tags) off the mesh; live state
              and state-derived geometry live in `state_cache`.
- `GeomT`:    `MeshGeometry` (built from the mesh at construction, not updated
              afterwards; not read by the 1D residual)
- `LayoutT`:  `StateLayout` descriptor
- `RadCacheT`:  `Properties.RadiationCache` for P1; `Nothing` otherwise
- `StateCacheT`: `StateCache{NC, Tu, Mode}` — typed mirror of the live
                 state populated at every residual entry
- `DryIdxT`:  `NTuple` of the dry-pyrolysis-gas component indices
- `DarcyT`:   the `PropertyCache`'s Darcy sub-cache type (`Nothing` for
              Fickian transport)
"""
mutable struct Workspace{
    NC,
    Mode<:SimulationMode,
    Tu,
    MeshT,
    GeomT,
    LayoutT,
    RadCacheT,
    StateCacheT,
    DryIdxT,
    DarcyT,
}
    # ---- references to the working mesh / geometry / layout ----
    mesh::MeshT
    geometry::GeomT
    layout::LayoutT

    # ---- typed live-state mirror — single source of truth in residual ----
    state_cache::StateCacheT

    # ---- always-present caches (parameterized on Tu) ----
    # The `DarcyT` parameter pins the PropertyCache's concrete
    # type so `ws.prop_cache.darcy` is type-stable downstream.
    prop_cache::Properties.PropertyCache{Tu,DarcyT}
    rhs_cache::Discretization.RHSCache{Tu}

    # ---- mode-conditional caches ----
    rad_cache::RadCacheT       # RadiationCache (P1) or Nothing

    # ---- physics work arrays (always allocated to live cell count) ----
    dT_work::Vector{Tu}
    dξ_work::Matrix{Tu}

    # ---- ALE work arrays (empty for Eulerian) ----
    w_work::Vector{Tu}
    # Reaction-driven (physical) mesh velocity, saved before the `w_override`
    # hook replaces `w_work`. The energy equation advects T by the velocity of
    # the mesh RELATIVE to the mesh-bound condensed material,
    # w_rel = w_work − w_physical — identically zero in production (the mesh
    # rides the collapsing matrix), non-zero only under `w_override`. Only
    # written when the override is active.
    w_physical::Vector{Tu}
    θ_work::Vector{Tu}
    ale_advection::Vector{Tu}
    dξ_dt_work::Vector{Tu}

    # ---- pre-cached boundary / index info ----
    # dry_gas_indices stored as NTuple{N,Int} (parameter
    # `DryIdxT`) so per-residual iteration unrolls and avoids the runtime-
    # string-allocating @generated lookup at each call.
    dry_gas_indices::DryIdxT
    top_face_ids::Vector{Int}
    bottom_face_ids::Vector{Int}

    # ---- runtime-mutable settings ----
    I_radiation::Float64

    # ---- verification hook: prescribed ALE mesh velocity ----
    # When non-nothing (ALE mode only), the reaction-driven mesh velocity is
    # REPLACED node-by-node with `w_override(z_node, t)` after the physical
    # computation (θ_work keeps its reaction-driven values, so the dilation
    # correction stays physical). Used by the MMS verification suite to test
    # GCL / moving-mesh transport against manufactured solutions
    # (free-stream preservation, prescribed smooth motion). The caller is
    # responsible for keeping the bottom node pinned (w(z_bottom, t) = 0).
    # NOTE: the structured geometry Jacobian assumes reaction-driven motion —
    # use a finite-difference Jacobian when overriding (as the MMS driver does).
    w_override::Union{Nothing,Function}
end

# ---------------------------------------------------------------------------
# Trait-dispatched no-ops on `Nothing`-typed caches
# ---------------------------------------------------------------------------

"""
	update_radiation!(ws::Workspace, args...) -> Nothing

Refresh the radiation cache. The `Nothing` specialization is a no-op
(compiles away under specialization).
"""
@inline update_radiation!(
    ::Workspace{NC,Mode,Tu,MeshT,GeomT,LayoutT,Nothing},
    args...,
) where {NC,Mode,Tu,MeshT,GeomT,LayoutT} = nothing

# Internal: the face ids of a boundary of a Mesh1D. 1D meshes key boundaries
# by integer (1 = top, 2 = bottom); the user-facing tag→id translation lives in
# `Geometry.tag_to_boundary_id`.
_boundary_face_ids(mesh, boundary_id::Int) =
    collect(get(mesh.boundary_faces, boundary_id, Int[]))
