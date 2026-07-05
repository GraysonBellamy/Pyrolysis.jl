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
- `GeomT`:    `MeshGeometry` view
- `LayoutT`:  `StateLayout` descriptor
- `GeomCacheT`: `Geometry.GeometryCache` for ALE; `Nothing` for Eulerian
- `RadCacheT`:  `Properties.RadiationCache` for P1; `Nothing` otherwise
- `BCStateT`:   per-boundary-tag state container (NamedTuple)
- `StateCacheT`: `StateCache{NC, Tu, Mode}` — typed mirror of the live
                 state populated at every residual entry
"""
mutable struct Workspace{
    NC,
    Mode<:SimulationMode,
    Tu,
    MeshT,
    GeomT,
    LayoutT,
    GeomCacheT,
    RadCacheT,
    BCStateT,
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
    geom_cache::GeomCacheT     # GeometryCache (ALE) or Nothing
    rad_cache::RadCacheT       # RadiationCache (P1) or Nothing

    # ---- physics work arrays (always allocated to live cell count) ----
    dT_work::Vector{Tu}
    dξ_work::Matrix{Tu}

    # ---- ALE work arrays (empty for Eulerian) ----
    w_work::Vector{Tu}
    w_face::Vector{Tu}
    θ_work::Vector{Tu}
    ale_advection::Vector{Tu}
    ale_species_advection::Matrix{Tu}
    rates_work::Vector{Tu}
    heats_work::Vector{Tu}
    dξ_dt_work::Vector{Tu}
    S_conv_work::Vector{Tu}

    # ---- pre-cached boundary / index info ----
    # dry_gas_indices stored as NTuple{N,Int} (parameter
    # `DryIdxT`) so per-residual iteration unrolls and avoids the runtime-
    # string-allocating @generated lookup at each call.
    dry_gas_indices::DryIdxT
    top_face_ids::Vector{Int}
    bottom_face_ids::Vector{Int}
    top_cell_ids::Vector{Int}
    bottom_cell_ids::Vector{Int}

    # ---- runtime-mutable settings ----
    I_radiation::Float64
    ale_limiter::Discretization.Limiter
    min_thickness::Float64

    # ---- per-tag boundary state (NamedTuple keyed by tag) ----
    boundary_states::BCStateT
end

@inline n_components(::Workspace{NC}) where {NC} = NC
@inline simulation_mode(::Workspace{NC,Mode}) where {NC,Mode} = Mode.instance
@inline scalar_eltype(::Workspace{NC,Mode,Tu}) where {NC,Mode,Tu} = Tu

# ---------------------------------------------------------------------------
# Trait predicates used by residual / Jacobian
# ---------------------------------------------------------------------------

"""
	is_ale_workspace(ws::Workspace) -> Bool

Compile-time check: does the workspace's mode have ALE mesh motion?
"""
@inline is_ale_workspace(::Workspace{NC,Mode}) where {NC,Mode} =
    mesh_mode(Mode.instance) isa ALE

@inline is_eulerian_workspace(::Workspace{NC,Mode}) where {NC,Mode} =
    mesh_mode(Mode.instance) isa Eulerian

# ---------------------------------------------------------------------------
# Trait-dispatched no-ops on `Nothing`-typed caches
# ---------------------------------------------------------------------------

"""
	update_radiation!(ws::Workspace, args...) -> Nothing

Refresh the radiation cache. The `Nothing` specialization is a no-op
(compiles away under specialization).
"""
@inline update_radiation!(
    ::Workspace{NC,Mode,Tu,MeshT,GeomT,LayoutT,GeomCacheT,Nothing},
    args...,
) where {NC,Mode,Tu,MeshT,GeomT,LayoutT,GeomCacheT} = nothing

# ---------------------------------------------------------------------------
# Resize entry points
# ---------------------------------------------------------------------------

"""
	resize_workspace!(ws::Workspace, n_nodes_new, n_cells_new, n_faces_new)

Resize all internal buffers in-place. Single entry point for post-AMR or
post-h-refinement re-sizing; walks every owned buffer and delegates to its
own resize routine.
"""
function resize_workspace!(
    ws::Workspace{NC,Mode,Tu},
    n_nodes_new::Int,
    n_cells_new::Int,
    n_faces_new::Int,
) where {NC,Mode,Tu}
    Geometry.resize_geometry!(ws.geometry, n_nodes_new, n_cells_new, n_faces_new)

    # Layout descriptor is rebuilt with the new sizes.
    ws.layout = StateLayout(NC, n_cells_new, n_nodes_new, Mode.instance)

    # Typed state mirror — resize before the kernels touch it.
    resize_state_cache!(ws.state_cache, ws.layout, n_faces_new)

    # Always-present caches
    Properties.resize_cache!(ws.prop_cache, n_cells_new, n_faces_new)
    NR = length(ws.rhs_cache.rates)
    Discretization.resize_rhs_cache!(ws.rhs_cache, n_cells_new, n_faces_new, NC, NR)

    # Mode-conditional caches
    if ws.geom_cache !== nothing
        Geometry.resize_geometry_cache!(ws.geom_cache, n_cells_new, n_faces_new)
    end
    if ws.rad_cache !== nothing
        Properties.resize_radiation_cache!(ws.rad_cache, n_cells_new)
    end

    # Physics work arrays
    resize!(ws.dT_work, n_cells_new)
    ws.dξ_work = zeros(Tu, NC, n_cells_new)

    # ALE work arrays
    if mesh_mode(Mode.instance) isa ALE
        resize!(ws.w_work, n_nodes_new)
        resize!(ws.w_face, n_faces_new)
        resize!(ws.θ_work, n_cells_new)
        resize!(ws.ale_advection, n_cells_new)
        ws.ale_species_advection = zeros(Tu, NC, n_cells_new)
        resize!(ws.S_conv_work, n_cells_new)
    end

    # Refresh pre-cached boundary IDs from the mesh (mesh topology changed)
    top_face_ids, top_cell_ids = _boundary_face_cells(ws.mesh, 1)
    bottom_face_ids, bottom_cell_ids = _boundary_face_cells(ws.mesh, 2)
    ws.top_face_ids = top_face_ids
    ws.bottom_face_ids = bottom_face_ids
    ws.top_cell_ids = top_cell_ids
    ws.bottom_cell_ids = bottom_cell_ids

    return ws
end

# Internal: pull the (face_ids, cell_ids) for a boundary id off a Mesh1D.
# 1D meshes key boundaries by integer (1 = top, 2 = bottom); the user-facing
# tag→id translation lives in `Geometry.tag_to_boundary_id`.
function _boundary_face_cells(mesh, boundary_id::Int)
    face_ids = Int[]
    cell_ids = Int[]
    if haskey(mesh.boundary_faces, boundary_id)
        for f in mesh.boundary_faces[boundary_id]
            push!(face_ids, f)
            face = mesh.faces[f]
            left_id, right_id = face.cell_ids
            cell_id = left_id == 0 ? right_id : left_id
            push!(cell_ids, cell_id)
        end
    end
    return face_ids, cell_ids
end
