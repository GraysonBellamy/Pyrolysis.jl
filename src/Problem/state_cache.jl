"""
StateCache — typed, per-residual state mirror.

Holds T, ξ, z, χ and the geometry quantities (V_i, A_face_f, A_section)
that depend on those state variables. `populate!` reads `u` (whatever
scalar type it carries — `Float64`, `ForwardDiff.Dual`, …) and writes
into the cache without touching the working `Mesh1D`.

The cache is the **single source of truth** for state and state-derived
geometry inside the residual / Jacobian assembly. Leaf kernels read
T/ξ/z/χ/V/A_face from this cache, never from `mesh.cell_states[i]` or
`mesh.cell_volumes[i]`.

The mutable mesh is updated only by `StepUpdateCallback` between
integrator-accepted timesteps; during residual assembly it is read-only
(connectivity only).

# Type Parameters
- `NC`:   number of components
- `Tu`:   scalar eltype of the live state (`Float64` for the canonical
          workspace; `ForwardDiff.Dual{Tag, Float64, N}` for the AD
          workspace)
- `Mode`: bundled `SimulationMode` matching the parent `Workspace` /
          `ProblemDef`

# Fields
- `T::Vector{Tu}`           — temperatures, length `n_cells`
- `ξ::Matrix{Tu}`           — concentrations, `NC × n_cells`
- `z::Vector{Tu}`           — node positions, length `n_nodes`
                              (length 0 for non-ALE)
- `χ::Vector{Tu}`           — pyrolysis-progress, length `n_cells`
                              (length 0 for `NoChi`)
- `volume::Vector{Tu}`      — cell volumes, length `n_cells`
                              (mirrors `mesh.cell_volumes` for Eulerian;
                              recomputed from `z` and `A_section` for ALE)
- `A_face::Vector{Tu}`      — face cross-sectional areas, length `n_faces`
                              (mirrors `mesh.face_areas` for Eulerian)
- `cell_centers::Vector{Tu}`   — axial cell-center positions, length `n_cells`
                              (mirrors `mesh.cell_centers` for Eulerian;
                              recomputed from `z` for ALE)
- `face_positions::Vector{Tu}` — axial face positions, length `n_faces`
                              (mirrors `mesh.face_positions` for Eulerian;
                              equals `z` for ALE in 1D)
- `A_section::Base.RefValue{Tu}` — current column cross-sectional area
                              (`A_0 · law(χ̄)` for `WithChi`;
                              `mesh.cross_section_area_initial` otherwise)
"""
mutable struct StateCache{NC,Tu,Mode<:SimulationMode}
    T::Vector{Tu}
    ξ::Matrix{Tu}
    z::Vector{Tu}
    χ::Vector{Tu}
    volume::Vector{Tu}
    A_face::Vector{Tu}
    cell_centers::Vector{Tu}
    face_positions::Vector{Tu}
    A_section::Base.RefValue{Tu}
end

"""
	StateCache{NC, Tu}(layout::StateLayout, n_faces::Int, A_section_init)

Allocate a `StateCache` sized for the given layout and face count, with
the column cross-sectional area initialized to `A_section_init`.

For ALE layouts the `z` block has length `n_nodes`; for non-ALE the
block is length 0. Same for `χ` under `WithChi` / `NoChi`.

`volume` and `A_face` are always allocated (they are read by every
residual). For an Eulerian Float64 workspace, the caller is expected to
copy `mesh.cell_volumes` and `mesh.face_areas` into them once at
construction; for ALE they are recomputed each `update_geometry!` call.
"""
function StateCache{NC,Tu}(
    layout::StateLayout{NC,Mode},
    n_faces::Int,
    A_section_init::Real,
) where {NC,Tu,Mode}
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes
    z_len = mesh_mode(layout.mode) isa ALE ? n_nodes : 0
    χ_len = chi_mode(layout.mode) isa WithChi ? n_cells : 0
    StateCache{NC,Tu,Mode}(
        zeros(Tu, n_cells),
        zeros(Tu, NC, n_cells),
        zeros(Tu, z_len),
        zeros(Tu, χ_len),
        zeros(Tu, n_cells),
        zeros(Tu, n_faces),
        zeros(Tu, n_cells),
        zeros(Tu, n_faces),
        Ref{Tu}(Tu(A_section_init)),
    )
end

@inline n_components(::StateCache{NC}) where {NC} = NC

"""
	StateCache(mesh::Mesh1D{NC}; mode = SimulationMode())
	  -> StateCache{NC, Float64, typeof(mode)}

Build a Float64 `StateCache` whose blocks are populated from the live
`mesh.cell_states[i]` (T, ξ) and `mesh.nodes[i].z`, with geometry mirrored
from the mesh. Convenience constructor used by callers that want to
evaluate a kernel against a mesh-resident state outside the residual /
workspace pipeline (tests, post-processing, manual analyses).

The default `mode` is `SimulationMode()` (Eulerian/no-radiation/Fickian/no-χ);
override only if you need a non-default block layout (ALE / WithChi / …).
"""
function StateCache(mesh::Mesh1D{NC}; mode::SimulationMode = SimulationMode()) where {NC}
    n_cells = length(mesh.cells)
    n_faces = length(mesh.faces)
    n_nodes = length(mesh.nodes)
    layout = StateLayout(NC, n_cells, n_nodes, mode)
    cache = StateCache{NC,Float64}(layout, n_faces, mesh.cross_section_area)
    mirror_eulerian_geometry!(cache, mesh)

    # Mirror cell-centered T/ξ from the mesh.
    @inbounds for i = 1:n_cells
        cs = mesh.cell_states[i]
        cache.T[i] = cs.T
        for j = 1:NC
            cache.ξ[j, i] = cs.ξ[j]
        end
    end

    # z block (ALE only).
    if mesh_mode(mode) isa ALE
        @inbounds for i = 1:n_nodes
            cache.z[i] = mesh.nodes[i].z
        end
    end

    return cache
end

# ----------------------------------------------------------------------------
# Cross-section update
# ----------------------------------------------------------------------------

"""
	update_cross_section_area!(cache::StateCache, material, mesh,
	                            u::AbstractVector, layout, ::ChiMode)

Refresh `cache.A_section[]` from the χ block of the state vector under the
column-global lateral-shrinkage law.

For `NoChi` (or an identity-law material) this is a no-op — `A_section` was
seeded at workspace construction and never changes.

For `WithChi`, χ is mass-weighted by `mesh.m_dry_solid_initial` (Float64
constants) to form the column-average `χ̄`, and the new area is

    A = A_0 · law(χ̄)

`A_0` comes from `mesh.cross_section_area_initial`. Carries `eltype(u)` end
-to-end (no `Float64(real(...))` strip) so AD passes propagate sensitivities
through `A` into every kernel reading `cache.A_section`. Returns
`(A, χ̄, M)`.
"""
@inline update_cross_section_area!(
    cache::StateCache,
    material,
    mesh,
    u::AbstractVector,
    layout,
    ::NoChi,
) = (cache.A_section[], zero(eltype(u)), 0.0)

function update_cross_section_area!(
    cache::StateCache{NC,Tu},
    material,
    mesh,
    u::AbstractVector,
    layout::StateLayout{NC},
    ::WithChi,
) where {NC,Tu}
    # Identity law: A is constant; skip the χ̄ scan entirely.
    if material.lateral_shrinkage_law === nothing
        return (cache.A_section[], zero(Tu), 0.0)
    end

    χ̄, M = _chi_bar_from_state(u, layout, mesh)
    A = mesh.cross_section_area_initial * material.lateral_shrinkage_law(χ̄)
    ForwardDiff.value(A) > 0.0 || error(
        "lateral_shrinkage_law(χ̄ = $(ForwardDiff.value(χ̄))) gives a non-positive " *
        "cross-section area A = $(ForwardDiff.value(A)) m². The law must stay " *
        "strictly positive on χ ∈ [0, 1]; non-positive A makes every cell volume " *
        "non-positive.",
    )
    cache.A_section[] = A
    return (A, χ̄, M)
end

# Mass-weighted column average of χ from `u`, scalar-type-generic.
@inline function _chi_bar_from_state(
    u::AbstractVector,
    layout::StateLayout{NC},
    mesh,
) where {NC}
    m_init = mesh.m_dry_solid_initial
    (isempty(m_init)) && return (zero(eltype(u)), 0.0)

    rng = chi_range(layout)
    isempty(rng) && return (zero(eltype(u)), 0.0)
    χ_offset = first(rng) - 1
    n = layout.n_cells
    length(m_init) == n || throw(
        DimensionMismatch(
            "mesh.m_dry_solid_initial has $(length(m_init)) entries for " *
            "$n cells — the χ denominator vector must be compacted in " *
            "lockstep with the state vector",
        ),
    )

    num = zero(eltype(u))
    M = 0.0
    @inbounds for i = 1:mesh.n_active_cells
        m_i = m_init[i]
        num += u[χ_offset+i] * m_i
        M += m_i
    end
    M <= eps() && return (zero(eltype(u)), 0.0)
    return (num / M, M)
end

"""
	populate!(cache::StateCache, u::AbstractVector, layout::StateLayout) -> cache

Read the flat state vector `u` into the typed cache. Carries whatever
scalar type `u` has end-to-end — no projection to `Float64`.

The block layout is dictated by `layout`:
- T:  `u[1 : n]`
- ξ:  `u[n + (j-1)*n + i]` for `j ∈ 1:NC`, `i ∈ 1:n`
- z:  `u[n*(1+NC) + i]` for `i ∈ 1:n_nodes` (ALE only)
- χ:  trailing block of length `n` (WithChi only)

The cache, not the mesh, is the single source of truth during residual
assembly.
"""
@inline function populate!(
    cache::StateCache{NC,Tu,Mode},
    u::AbstractVector,
    layout::StateLayout{NC,Mode},
) where {NC,Tu,Mode}
    n = layout.n_cells

    @inbounds for i = 1:n
        cache.T[i] = u[i]
    end

    @inbounds for j = 1:NC, i = 1:n
        cache.ξ[j, i] = u[n+(j-1)*n+i]
    end

    if mesh_mode(layout.mode) isa ALE
        zoff = n * (1 + NC)
        @inbounds for i = 1:layout.n_nodes
            cache.z[i] = u[zoff+i]
        end
    end

    if chi_mode(layout.mode) isa WithChi
        χoff = first(chi_range(layout)) - 1
        @inbounds for i = 1:n
            cache.χ[i] = u[χoff+i]
        end
    end

    return cache
end

"""
	mirror_eulerian_geometry!(cache::StateCache, mesh) -> cache

Copy the immutable Eulerian geometry (`cell_volumes`, `face_areas`,
`cross_section_area`) from `mesh` into the cache. Used at workspace
construction so kernels reading from `cache.volume` / `cache.A_face` /
`cache.A_section` get the same values they would read directly from the
mesh.

For ALE workspaces this initial mirror is overwritten on the first
`update_geometry!(cache)` call.
"""
function mirror_eulerian_geometry!(cache::StateCache{NC,Tu}, mesh) where {NC,Tu}
    n_cells = length(cache.volume)
    n_faces = length(cache.A_face)
    @inbounds for i = 1:n_cells
        cache.volume[i] = Tu(mesh.cell_volumes[i])
        cache.cell_centers[i] = Tu(mesh.cell_centers[i])
    end
    @inbounds for f = 1:n_faces
        cache.A_face[f] = Tu(mesh.face_areas[f])
        cache.face_positions[f] = Tu(mesh.face_positions[f])
    end
    cache.A_section[] = Tu(mesh.cross_section_area)
    return cache
end

"""
	update_state_geometry!(cache::StateCache, mesh, ::Eulerian) -> cache

For Eulerian meshes the geometry is constant and was mirrored at workspace
construction; this is a no-op.

	update_state_geometry!(cache::StateCache, mesh, ::ALE) -> cache

For ALE, recompute face positions, face areas, cell centers, and cell
volumes from `cache.z` and `cache.A_section[]`. The math mirrors
`Geometry.update_geometry!(mesh::Mesh1D)` but writes into the typed cache,
preserving `eltype(cache.z) == Tu` end-to-end (so Dual-typed AD passes
carry geometry sensitivities through).
"""
@inline update_state_geometry!(cache::StateCache, mesh, ::Eulerian) = cache

function update_state_geometry!(cache::StateCache{NC,Tu}, mesh, ::ALE) where {NC,Tu}
    A = cache.A_section[]

    @inbounds for f = 1:length(cache.A_face)
        face = mesh.faces[f]
        cache.face_positions[f] = cache.z[face.node_id]
        cache.A_face[f] = A
    end

    @inbounds for i = 1:length(cache.volume)
        cell = mesh.cells[i]
        z_left = cache.z[cell.node_ids[1]]
        z_right = cache.z[cell.node_ids[2]]
        cache.cell_centers[i] = (z_left + z_right) / 2
        cache.volume[i] = abs(z_right - z_left) * A
    end

    return cache
end
