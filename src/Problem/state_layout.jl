"""
StateLayout — unified state-vector descriptor.

Describes the packing of the state vector `u` in terms of the bundled
[`SimulationMode`](@ref). Block layout (T, then ξ, then optional z, then
optional χ) is the single source of truth for every pack / unpack /
view operation in the package.

Block presence is **strictly** determined by `SimulationMode`:
- T and ξ blocks are always present.
- z block is present iff `mesh_mode(mode) isa ALE`.
- χ block is present iff `chi_mode(mode) isa WithChi`.

# Type Parameters
- `NC`: number of components
- `Mode <: SimulationMode`: bundled trait

# Block-major offsets (1-based, with `n = n_cells`)
- T:  `1 .. n`
- ξ:  `n + (j-1)*n + (1..n)` for `j ∈ 1:NC`
- z:  `n*(1+NC) + (1..n_nodes)` (ALE only)
- χ:  trailing block of length `n` (only if `WithChi`)
"""
struct StateLayout{NC,Mode<:SimulationMode}
    n_cells::Int
    n_nodes::Int
    mode::Mode
end

"""
	StateLayout(NC::Int, n_cells::Int, n_nodes::Int, mode::SimulationMode)

Build a `StateLayout` for the given component count, mesh size, and mode.
"""
@inline StateLayout(NC::Int, n_cells::Int, n_nodes::Int, mode::SimulationMode) =
    StateLayout{NC,typeof(mode)}(n_cells, n_nodes, mode)

@inline n_components(::StateLayout{NC}) where {NC} = NC

"""
	state_length(l::StateLayout) -> Int

Return the total length of the state vector encoded by this layout.
"""
@inline function state_length(l::StateLayout{NC}) where {NC}
    n = l.n_cells
    len = n * (1 + NC)
    mesh_mode(l.mode) isa ALE && (len += l.n_nodes)
    chi_mode(l.mode) isa WithChi && (len += n)
    return len
end

"""
	temperature_range(l::StateLayout) -> UnitRange

Index range for the temperature block.
"""
@inline temperature_range(l::StateLayout) = 1:l.n_cells

"""
	concentration_range(l::StateLayout, j::Int) -> UnitRange

Index range for component `j`'s concentration block.
"""
@inline function concentration_range(l::StateLayout{NC}, j::Int) where {NC}
    @boundscheck (1 <= j <= NC) || throw(BoundsError(l, j))
    n = l.n_cells
    start = n + (j - 1) * n + 1
    return start:(start+n-1)
end

"""
	z_range(l::StateLayout) -> UnitRange

Index range for the ALE node-position block. Returns an empty range for
non-ALE layouts.
"""
@inline function z_range(l::StateLayout{NC}) where {NC}
    mesh_mode(l.mode) isa ALE || return 1:0
    start = l.n_cells * (1 + NC) + 1
    return start:(start+l.n_nodes-1)
end

"""
	chi_range(l::StateLayout) -> UnitRange

Index range for the χ block. Empty unless `chi_mode(l.mode) isa WithChi`.
"""
@inline function chi_range(l::StateLayout{NC}) where {NC}
    chi_mode(l.mode) isa WithChi || return 1:0
    n = l.n_cells
    base = n * (1 + NC)
    mesh_mode(l.mode) isa ALE && (base += l.n_nodes)
    return (base+1):(base+n)
end

@inline z_offset(l::StateLayout) = isempty(z_range(l)) ? 0 : first(z_range(l)) - 1
@inline chi_offset(l::StateLayout) = isempty(chi_range(l)) ? 0 : first(chi_range(l)) - 1
@inline concentration_offset(l::StateLayout) = l.n_cells

# ============================================================================
# Compile-time block accessors via Val{:T} / Val{:ξ} / Val{:z} / Val{:χ}
# ============================================================================

"""
	block_view(u, l::StateLayout, ::Val{:T}) -> view of T block
	block_view(u, l::StateLayout, ::Val{:ξ}, j::Int) -> view of component j
	block_view(u, l::StateLayout, ::Val{:z}) -> view of node-position block (ALE)
	block_view(u, l::StateLayout, ::Val{:χ}) -> view of pyrolysis-progress block (WithChi)

Compile-time-dispatched view into the named block of `u`. Returns an
empty view when the block does not exist for the layout's mode.
"""
@inline block_view(u::AbstractVector, l::StateLayout, ::Val{:T}) =
    @view u[temperature_range(l)]

@inline function block_view(
    u::AbstractVector,
    l::StateLayout{NC},
    ::Val{:ξ},
    j::Int,
) where {NC}
    @view u[concentration_range(l, j)]
end

@inline block_view(u::AbstractVector, l::StateLayout, ::Val{:z}) = @view u[z_range(l)]

@inline block_view(u::AbstractVector, l::StateLayout, ::Val{:χ}) = @view u[chi_range(l)]

# ============================================================================
# Pack / unpack — single entry points dispatched on the layout's mode
# ============================================================================

"""
	pack_state!(u, mesh, layout::StateLayout) -> u

Pack the mesh's `cell_states` (and node positions when ALE; χ when WithChi)
into the flat ODE state vector `u`. Mode is read from the layout's type
parameter.
"""
function pack_state!(u::AbstractVector{Float64}, mesh, l::StateLayout{NC}) where {NC}
    n = l.n_cells

    @inbounds for i = 1:n
        state = mesh.cell_states[i]
        u[i] = state.T
        for j = 1:NC
            u[n+(j-1)*n+i] = state.ξ[j]
        end
    end

    if mesh_mode(l.mode) isa ALE
        zoff = l.n_cells * (1 + NC)
        @inbounds for i = 1:l.n_nodes
            u[zoff+i] = mesh.nodes[i].z
        end
    end

    if chi_mode(l.mode) isa WithChi
        # χ is reaction-driven and rebuilt at every RHS evaluation; pack zeros.
        χoff = first(chi_range(l)) - 1
        @inbounds for i = 1:n
            u[χoff+i] = 0.0
        end
    end

    return u
end

"""
	update_mesh_from_state!(mesh, u::AbstractVector{Float64}, layout::StateLayout) -> mesh

Mirror an integrator-accepted Float64 state vector into the mutable mesh.
Called by `StepUpdateCallback` between timesteps so AMR, depletion
checks, and post-processing read a consistent mesh.

This function is **not** invoked during residual / Jacobian assembly —
those are pure functions of `u` and the typed `state_cache`. The mesh is
read-only during assembly.
"""
function update_mesh_from_state!(
    mesh,
    u::AbstractVector{Float64},
    l::StateLayout{NC},
) where {NC}
    n = l.n_cells

    @inbounds for i = 1:n
        state = mesh.cell_states[i]
        state.T = u[i]
        state.ξ = ntuple(Val(NC)) do j
            u[n+(j-1)*n+i]
        end
    end

    if mesh_mode(l.mode) isa ALE
        zoff = l.n_cells * (1 + NC)
        @inbounds for i = 1:l.n_nodes
            mesh.nodes[i].z = u[zoff+i]
        end
        Geometry.update_geometry!(mesh)
    end

    return mesh
end

"""
	create_state_vector(layout::StateLayout, mesh) -> Vector{Float64}

Allocate a fresh state vector sized to `layout` and pack the current
mesh state into it. Used by `solve()` to build `u0`.
"""
function create_state_vector(l::StateLayout, mesh)
    u = Vector{Float64}(undef, state_length(l))
    pack_state!(u, mesh, l)
    return u
end

# ============================================================================
# Block accessors useful in extract / postprocess
# ============================================================================

"""
	get_temperature(u, layout, i) -> Real
"""
@inline function get_temperature(u::AbstractVector, l::StateLayout, i::Int)
    @boundscheck 1 <= i <= l.n_cells || throw(BoundsError(u, i))
    @inbounds u[i]
end

"""
	get_concentration(u, layout, i, j) -> Real
"""
@inline function get_concentration(
    u::AbstractVector,
    l::StateLayout{NC},
    i::Int,
    j::Int,
) where {NC}
    @boundscheck 1 <= i <= l.n_cells || throw(BoundsError())
    @boundscheck 1 <= j <= NC || throw(BoundsError())
    n = l.n_cells
    @inbounds u[n+(j-1)*n+i]
end

"""
	get_node_position(u, layout, i) -> Real

Only valid for ALE layouts.
"""
@inline function get_node_position(u::AbstractVector, l::StateLayout{NC}, i::Int) where {NC}
    @boundscheck 1 <= i <= l.n_nodes || throw(BoundsError())
    zoff = l.n_cells * (1 + NC)
    @inbounds u[zoff+i]
end

"""
	get_chi(u, layout, i) -> Real

Read χ_i. Returns zero(eltype(u)) when the layout has no χ block.
"""
@inline function get_chi(u::AbstractVector, l::StateLayout, i::Int)
    rng = chi_range(l)
    isempty(rng) && return zero(eltype(u))
    @boundscheck 1 <= i <= l.n_cells || throw(BoundsError())
    @inbounds u[first(rng)+i-1]
end

"""
	get_thickness_from_state(u, layout) -> Real

Compute current domain thickness `z_top - z_bottom` from the ALE state.
Returns 0 for Eulerian layouts.
"""
@inline function get_thickness_from_state(u::AbstractVector, l::StateLayout)
    mesh_mode(l.mode) isa ALE || return zero(eltype(u))
    z_bottom = get_node_position(u, l, 1)
    z_top = get_node_position(u, l, l.n_nodes)
    return z_top - z_bottom
end

@inline node_positions_view(u::AbstractVector, l::StateLayout) = @view u[z_range(l)]

"""
	get_node_positions(u, layout) -> view

View of node positions for ALE layouts.
"""
@inline get_node_positions(u::AbstractVector, l::StateLayout) = node_positions_view(u, l)

"""
	compute_volumes_from_state(u, layout, A) -> Vector{Float64}

Compute cell volumes from node positions in an ALE state vector.
"""
function compute_volumes_from_state(u::AbstractVector, l::StateLayout, A::Float64)
    n = l.n_cells
    volumes = Vector{Float64}(undef, n)
    @inbounds for i = 1:n
        z_left = get_node_position(u, l, i)
        z_right = get_node_position(u, l, i + 1)
        volumes[i] = abs(z_right - z_left) * A
    end
    return volumes
end
