"""
Solver callbacks.

ALE mesh motion is integrated directly into the unified workspace.
Thickness history is extracted from the state vector in extract_solution_ale!.
"""

using DiffEqCallbacks: DiscreteCallback, PeriodicCallback, ContinuousCallback
using ..Adaptivity:
    find_surface_cell,
    find_interior_neighbor,
    find_exterior_neighbor,
    merge_adjacent_cells!,
    compact_mesh_after_merge!,
    find_merge_neighbor,
    merge_temperatures_conservatively_from_values
using ..Problem:
    StateLayout,
    get_node_positions,
    get_temperature,
    get_concentration,
    get_node_position,
    get_chi,
    state_length,
    z_offset,
    chi_offset,
    chi_range,
    resize_workspace!,
    update_mesh_from_state!,
    mesh_mode,
    ALE
using ..Residual: cross_section_area_from_state
using ..Geometry: Mesh1D
using ..Materials: n_components, SOLID

# ============================================================================
# Cell Thickness Based Depletion Detection
# ============================================================================

"""
	is_cell_thin(u::AbstractVector, layout::StateLayout{NC}, 
				 cell_idx::Int, initial_thickness::Float64;
				 threshold::Float64=0.05) where NC -> Bool

Check if a cell has shrunk below a threshold fraction of its initial thickness.

This is the primary depletion detection mechanism for ALE simulations. When
a cell shrinks significantly, it should be merged with its neighbor to prevent:
1. Very small timesteps (CFL condition scales with Δx²)
2. Extreme aspect ratios between cells
3. Numerical instabilities from thin cells

# Arguments
- `u`: ODE state vector [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
- `layout`: StateLayout descriptor
- `cell_idx`: Index of the cell to check
- `initial_thickness`: Original thickness of **this** cell [m] — each cell is
  judged against its own reference (e.g. `state.initial_cell_thickness[cell_idx]`)
- `threshold`: Fraction of initial thickness below which cell is considered depleted

# Returns
- `true` if current cell thickness < threshold × initial_thickness

# Example
If threshold=0.05 (5%), a cell that started at 0.4mm and has shrunk to 0.02mm
(5% of original) will trigger merging.
"""
function is_cell_thin(
    u::AbstractVector,
    layout::StateLayout{NC},
    cell_idx::Int,
    initial_thickness::Float64;
    threshold::Float64 = 0.05,
) where {NC}
    # Get node positions from state vector
    z = get_node_positions(u, layout)

    # Compute current cell thickness
    # Cell i spans from node i to node i+1
    current_thickness = z[cell_idx+1] - z[cell_idx]

    # Check if cell has shrunk below threshold
    return current_thickness < threshold * initial_thickness
end

"""
	get_cell_thickness(u::AbstractVector, layout::StateLayout, cell_idx::Int) -> Float64

Get the current thickness of a specific cell from the state vector.
"""
function get_cell_thickness(u::AbstractVector, layout::StateLayout, cell_idx::Int)
    z = get_node_positions(u, layout)
    return z[cell_idx+1] - z[cell_idx]
end

"""
	find_thin_cells(u::AbstractVector, layout::StateLayout{NC},
					initial_thickness::AbstractVector{<:Real},
					threshold::Float64) where NC -> Vector{Int}

Find all cells thinner than `threshold` × their **own** initial thickness.

Each cell is compared against its own entry of the per-cell reference vector
`initial_thickness`. A single shared reference would misclassify on any
nonuniform mesh: with a max/min cell-size ratio above `1/threshold`, a
surface-cell reference either flags the fine cells at t = 0 or never flags
the coarse ones.

Returns the thin-cell IDs sorted by z-position in descending order (surface
cells first).

# Arguments
- `u`: ODE state vector [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
- `layout`: StateLayout descriptor
- `initial_thickness`: Initial thickness of each cell [m], one entry per cell
- `threshold`: Fraction of a cell's initial thickness below which it is thin

# Returns
Vector of cell IDs, sorted by z-position descending (surface first).

# Notes
Processing surface cells first (higher z) prevents "orphaning" interior cells
when multiple adjacent cells become thin simultaneously. This maintains mesh
connectivity during cascade scenarios.
"""
function find_thin_cells(
    u::AbstractVector,
    layout::StateLayout{NC},
    initial_thickness::AbstractVector{<:Real},
    threshold::Float64,
) where {NC}
    length(initial_thickness) == layout.n_cells || throw(
        DimensionMismatch(
            "initial_thickness has $(length(initial_thickness)) entries for " *
            "$(layout.n_cells) cells — the per-cell reference vector must be " *
            "compacted alongside the state vector after each merge",
        ),
    )
    thin = Int[]
    z = get_node_positions(u, layout)

    for i = 1:layout.n_cells
        thickness = z[i+1] - z[i]
        thickness > 0.0 || error(
            "Cell $i has non-positive thickness $(thickness) m: node positions " *
            "are no longer monotone in z (tangled mesh). This is a mesh-motion " *
            "failure, not depletion — refusing to merge the inverted cell away.",
        )
        if thickness < threshold * initial_thickness[i]
            push!(thin, i)
        end
    end

    # Sort by position descending (surface first)
    # Cell center approximated as average of node positions
    sort!(thin, by = i -> (z[i] + z[i+1]) / 2, rev = true)
    return thin
end

# ============================================================================
# ============================================================================
# Thickness Termination Callback (for ALE simulations)
# ============================================================================

"""
	create_thickness_termination_callback(sv_ref, min_thickness) -> ContinuousCallback

Create a callback that terminates the simulation when thickness falls below threshold.

This is used with ALE `StateLayout`s where node positions are part of
the state vector.

# Arguments
- `sv_ref`: Either a StateLayout directly, or a DepletionCallbackState (for shared reference)
- `min_thickness::Float64`: Minimum allowed thickness [m]

# Returns
A `ContinuousCallback` that terminates when thickness < min_thickness.

# Notes
When used with depletion handling, pass a DepletionCallbackState so that the
state vector descriptor is automatically updated when cells are merged.
"""
function create_thickness_termination_callback(sv_ref, min_thickness::Float64)
    # Condition function: returns 0 when thickness equals min_thickness
    condition = (u, t, integrator) -> begin
        # Get the current state vector descriptor (handles both direct layout and callback state)
        layout = _get_sv(sv_ref)
        z = get_node_positions(u, layout)
        thickness = z[end] - z[1]
        return thickness - min_thickness
    end

    # Affect function: terminate the integrator
    affect! = (integrator) -> terminate!(integrator)

    # Use ContinuousCallback for accurate event detection
    return ContinuousCallback(condition, affect!)
end

# ============================================================================
# Surface Cell Depletion Callback
# ============================================================================

"""
	DepletionCallbackState{NC}

Mutable state for tracking depletion callback progress.

This struct stores information needed across callback invocations,
particularly the current state vector descriptor which changes as cells are merged,
and the per-cell initial thicknesses for thickness-based depletion detection.

The reference vector is compacted in lockstep with the state vector: when a
thin cell is merged away, its entry is deleted and the kept cell retains its
own reference. (On a uniform mesh this reproduces the single-scalar behavior
exactly; on a stretched mesh it is what makes the thinness test per-cell.)

# Fields
- `layout::StateLayout{NC}`: Current state vector descriptor (updated after merges)
- `initial_cell_thickness::Vector{Float64}`: Initial thickness of each cell [m],
  one entry per currently-active cell
- `n_merges::Int`: Total number of cell merges performed
- `merge_times::Vector{Float64}`: Times at which merges occurred
- `merged_cells::Vector{Int}`: Original indices of cells that were merged (for diagnostics)
"""
mutable struct DepletionCallbackState{NC}
    layout::StateLayout{NC}
    initial_cell_thickness::Vector{Float64}
    n_merges::Int
    merge_times::Vector{Float64}
    merged_cells::Vector{Int}
end

function DepletionCallbackState(
    layout::StateLayout{NC},
    initial_cell_thickness::Vector{Float64},
) where {NC}
    length(initial_cell_thickness) == layout.n_cells || throw(
        DimensionMismatch(
            "initial_cell_thickness has $(length(initial_cell_thickness)) " *
            "entries for $(layout.n_cells) cells",
        ),
    )
    DepletionCallbackState{NC}(layout, initial_cell_thickness, 0, Float64[], Int[])
end

# Convenience constructor that computes per-cell initial thicknesses from the mesh
function DepletionCallbackState(layout::StateLayout{NC}, mesh::Mesh1D) where {NC}
    initial_thickness = mesh.cell_volumes ./ mesh.cross_section_area_initial
    DepletionCallbackState(layout, initial_thickness)
end

# Helper to extract state vector from either direct layout or DepletionCallbackState
_get_sv(layout::StateLayout) = layout
_get_sv(cb_state::DepletionCallbackState) = cb_state.layout

@inline function _set_temperature!(u, layout::StateLayout, i::Int, T)
    u[i] = T
    return u
end

@inline function _set_concentration!(
    u,
    layout::StateLayout{NC},
    i::Int,
    j::Int,
    ξ,
) where {NC}
    n = layout.n_cells
    u[n+(j-1)*n+i] = ξ
    return u
end

@inline function _set_chi!(u, layout::StateLayout, i::Int, χ)
    rng = chi_range(layout)
    isempty(rng) && return u
    u[first(rng)+i-1] = χ
    return u
end

function remove_cell_from_state_vector(
    u::AbstractVector{Float64},
    layout::StateLayout{NC},
    cell_id::Int,
    node_id::Int,
) where {NC}
    n_old = layout.n_cells
    n_nodes_old = layout.n_nodes
    n_new = n_old - 1
    n_nodes_new = n_nodes_old - 1
    n_new < 1 && error("Cannot remove cell: only 1 cell remaining")

    layout_new = StateLayout(NC, n_new, n_nodes_new, layout.mode)
    u_new = Vector{Float64}(undef, state_length(layout_new))

    dst = 0
    @inbounds for i = 1:n_old
        if i != cell_id
            dst += 1
            u_new[dst] = u[i]
        end
    end

    @inbounds for j = 1:NC
        dst = 0
        for i = 1:n_old
            if i != cell_id
                dst += 1
                src = n_old + (j - 1) * n_old + i
                dest = n_new + (j - 1) * n_new + dst
                u_new[dest] = u[src]
            end
        end
    end

    z_old = z_offset(layout)
    z_new = z_offset(layout_new)
    dst = 0
    @inbounds for i = 1:n_nodes_old
        if i != node_id
            dst += 1
            u_new[z_new+dst] = u[z_old+i]
        end
    end

    if !isempty(chi_range(layout))
        χ_old = chi_offset(layout)
        χ_new = chi_offset(layout_new)
        dst = 0
        @inbounds for i = 1:n_old
            if i != cell_id
                dst += 1
                u_new[χ_new+dst] = u[χ_old+i]
            end
        end
    end

    return u_new, layout_new
end

"""
	create_merged_state_vector(u, layout, cell_to_remove, cell_to_keep,
							   V_remove, V_keep;
							   material = nothing,
							   m_init_remove = 0.0, m_init_keep = 0.0)
		-> (u_new, layout_new)

Build the post-merge state vector: conservative ξ (volume-weighted bulk
concentrations), energy-conserving T (ρcₚ-weighted when `material` is given),
and — for WithChi layouts — the pyrolysis-progress merge

	χ_merged = (m_init_remove·χ_remove + m_init_keep·χ_keep) / (m_init_remove + m_init_keep)

weighted by the cells' *initial dry-solid masses* (`mesh.m_dry_solid_initial`),
because χᵢ = m_released,i / m_dry_init,i makes Σᵢ χᵢ·m_dry_init,i the extensive
invariant. Volume weighting would nearly discard the depleted cell's pyrolysis
record (V → 0 while its m_init is unchanged), stepping χ̄ down at each merge.
Callers must capture `m_init_remove` / `m_init_keep` **before** the mesh merge
sums and compacts `mesh.m_dry_solid_initial`; they are required (positive sum)
whenever the layout carries a χ block.
"""
function create_merged_state_vector(
    u::AbstractVector{Float64},
    layout::StateLayout{NC},
    cell_to_remove::Int,
    cell_to_keep::Int,
    V_remove::Float64,
    V_keep::Float64;
    material = nothing,
    m_init_remove::Float64 = 0.0,
    m_init_keep::Float64 = 0.0,
) where {NC}
    n_old = layout.n_cells
    n_old <= 1 && error("Cannot merge: only 1 cell remaining")
    V_total = V_remove + V_keep
    u_work = copy(u)

    T_remove = get_temperature(u_work, layout, cell_to_remove)
    T_keep = get_temperature(u_work, layout, cell_to_keep)
    ξ_remove = ntuple(Val(NC)) do j
        get_concentration(u_work, layout, cell_to_remove, j)
    end
    ξ_keep = ntuple(Val(NC)) do j
        get_concentration(u_work, layout, cell_to_keep, j)
    end

    ξ_merged = ntuple(Val(NC)) do j
        (V_remove * ξ_remove[j] + V_keep * ξ_keep[j]) / V_total
    end

    if material !== nothing
        # Canonical condensed-enthalpy-conserving merge (Adaptivity/cell_merge.jl).
        T_merged = merge_temperatures_conservatively_from_values(
            V_remove,
            T_remove,
            ξ_remove,
            V_keep,
            T_keep,
            ξ_keep,
            material,
        )
    else
        T_merged = (V_remove * T_remove + V_keep * T_keep) / V_total
    end

    _set_temperature!(u_work, layout, cell_to_keep, T_merged)
    for j = 1:NC
        _set_concentration!(u_work, layout, cell_to_keep, j, ξ_merged[j])
    end

    if !isempty(chi_range(layout))
        m_total = m_init_remove + m_init_keep
        m_total > 0.0 || throw(
            ArgumentError(
                "create_merged_state_vector: WithChi merge requires the cells' " *
                "initial dry-solid masses (m_init_remove, m_init_keep) captured " *
                "before the mesh merge",
            ),
        )
        χ_remove = get_chi(u_work, layout, cell_to_remove)
        χ_keep = get_chi(u_work, layout, cell_to_keep)
        χ_merged = (m_init_remove * χ_remove + m_init_keep * χ_keep) / m_total
        _set_chi!(u_work, layout, cell_to_keep, χ_merged)
    end

    node_to_remove = cell_to_remove > cell_to_keep ? cell_to_remove : cell_to_keep
    return remove_cell_from_state_vector(u_work, layout, cell_to_remove, node_to_remove)
end

"""
	create_depletion_callback(ws, mesh, material, layout;
							  threshold=0.05, min_cells=2) -> DiscreteCallback

Create a callback that detects and handles thin cells.

This callback checks at each timestep whether ANY cell has shrunk below a
threshold fraction of its initial thickness. When this occurs, the callback:
1. Finds all thin cells, sorted by z-position (surface first)
2. Merges the first thin cell with its best neighbor
3. Resizes the workspace's caches
4. Resizes the integrator's state vector using `resize!(integrator, ...)`
5. Remaps the solution to the new mesh

Only ONE merge is performed per callback invocation to allow the ODE solver
to adapt between merges.

This prevents numerical issues from extremely thin cells:
- Very small timesteps (CFL condition scales with Δx²)
- Extreme aspect ratios between cells
- Numerical instabilities

# Arguments
- `ws`: Mutable `Workspace` being solved
- `mesh::Mesh1D`: The mesh (will be modified during merges)
- `material`: Material definition (used for energy-conserving merge)
- `layout::StateLayout`: Initial state vector descriptor

# Keyword Arguments
- `threshold::Float64=0.05`: Thickness fraction threshold for depletion detection.
  Cell is considered depleted when current_thickness < threshold × its own
  initial thickness. Default 0.05 means cells are merged when shrunk to 5%
  of their original thickness.
- `min_cells::Int=2`: Minimum number of cells before terminating instead of merging.
  When only `min_cells` remain, simulation terminates rather than merging.

# Returns
A `DiscreteCallback` that handles cell depletion.

# Notes
- The callback uses `resize!(integrator, new_length)` to inform the ODE solver
  of the state vector size change
- Works with the unified ALE workspace and state layout
- The state vector descriptor in `ws.layout` is automatically updated
- All cells (surface and interior) are monitored for depletion

# Example
```julia
def = to_problem_def(problem; use_ale = true)
ws = build_workspace(problem, def)
u0 = create_state_vector(ws.layout, ws.mesh)
layout = ws.layout

# Merge cells when they shrink to 5% of original thickness
depletion_cb = create_depletion_callback(ws, mesh, material, layout;
										 threshold=0.05, min_cells=3)

prob = ODEProblem(ws, u0, tspan)
sol = solve(prob, Rodas5P(); callback=depletion_cb)
```
"""
function create_depletion_callback(
    ws,
    mesh::Mesh1D,
    material,
    sv_or_state;
    threshold::Float64 = 0.05,
    min_cells::Int = 2,
)

    # Accept either StateLayout or DepletionCallbackState
    # If just a StateLayout, create callback state with initial thickness from mesh
    if sv_or_state isa DepletionCallbackState
        cb_state = sv_or_state
    else
        cb_state = DepletionCallbackState(sv_or_state, mesh)
    end
    NC = typeof(cb_state).parameters[1]  # Get NC from DepletionCallbackState{NC}

    # Condition function: check if ANY cell has become too thin
    condition =
        (u, t, integrator) -> begin
            # Use the current state vector descriptor
            current_layout = cb_state.layout

            # Find ALL thin cells (not just surface). At n <= min_cells the
            # affect! terminates instead of merging, so the condition must
            # still fire — returning false here would leave cells free to
            # shrink indefinitely.
            thin_cells = find_thin_cells(
                u,
                current_layout,
                cb_state.initial_cell_thickness,
                threshold,
            )
            return !isempty(thin_cells)
        end

    # Affect function: merge the first thin cell (surface-first order)
    affect! =
        (integrator) -> begin
            current_layout = cb_state.layout
            n = current_layout.n_cells
            t = integrator.t

            # Check minimum cells
            if n <= min_cells
                @warn "At minimum cells ($min_cells), terminating simulation"
                terminate!(integrator)
                return
            end

            # Find thin cells (sorted surface-first)
            thin_cells = find_thin_cells(
                integrator.u,
                current_layout,
                cb_state.initial_cell_thickness,
                threshold,
            )

            if isempty(thin_cells)
                @warn "Condition triggered but no thin cell found"
                return
            end

            # Process only the FIRST thin cell (highest z / surface-most)
            cell_to_remove = thin_cells[1]

            # Find the best merge neighbor (prefers interior for thermal stability)
            cell_to_keep = find_merge_neighbor(mesh, cell_to_remove)
            if cell_to_keep == 0
                @warn "No valid merge neighbor for cell $cell_to_remove"
                return
            end

            # Get current cell thickness for logging
            current_thickness =
                get_cell_thickness(integrator.u, current_layout, cell_to_remove)
            thickness_percent =
                100 * current_thickness / cb_state.initial_cell_thickness[cell_to_remove]

            @info "Cell $cell_to_remove thin ($(round(thickness_percent, digits=1))% of initial) at t=$(round(t, digits=4))s, merging into cell $cell_to_keep, n_cells=$n"

            # Get cell volumes and initial dry-solid masses BEFORE mesh
            # modification: the merge sums m_dry_solid_initial into the kept
            # cell and the compaction reindexes it, so the χ-merge weights
            # must be captured here.
            V_remove = mesh.cell_volumes[cell_to_remove]
            V_keep = mesh.cell_volumes[cell_to_keep]
            has_m_init = !isempty(mesh.m_dry_solid_initial)
            m_init_remove = has_m_init ? mesh.m_dry_solid_initial[cell_to_remove] : 0.0
            m_init_keep = has_m_init ? mesh.m_dry_solid_initial[cell_to_keep] : 0.0

            # Store current state for conservative merge
            u_old = copy(integrator.u)

            # Merge cells in mesh (marks cell_to_remove as inactive)
            merged_id = merge_adjacent_cells!(
                mesh,
                cell_to_remove,
                cell_to_keep;
                material = material,
            )
            if merged_id == 0
                @warn "Failed to merge cell $cell_to_remove"
                return
            end

            # Compact the mesh arrays to remove the inactive cell
            # This is REQUIRED because the state vector actually removes elements,
            # so the mesh arrays must match in size
            compact_mesh_after_merge!(mesh, cell_to_remove)

            # Create new state vector with ENERGY-CONSERVING MERGE
            # Using material enables ρcₚ-weighted energy conservation instead of
            # simple volume-weighted temperature averaging
            u_new, layout_new = create_merged_state_vector(
                u_old,
                current_layout,
                cell_to_remove,
                cell_to_keep,
                V_remove,
                V_keep;
                material = material,
                m_init_remove = m_init_remove,
                m_init_keep = m_init_keep,
            )

            # Get dimensions from the new state vector descriptor (authoritative source)
            n_cells_new = layout_new.n_cells
            n_nodes_new = layout_new.n_nodes
            n_faces_new = n_cells_new + 1  # For 1D mesh
            new_state_length = length(u_new)

            # Update callback state. The per-cell reference vector is compacted
            # in lockstep with remove_cell_from_state_vector: the removed cell's
            # entry is dropped and the kept cell retains its own reference.
            cb_state.layout = layout_new
            deleteat!(cb_state.initial_cell_thickness, cell_to_remove)
            cb_state.n_merges += 1
            push!(cb_state.merge_times, t)
            push!(cb_state.merged_cells, cell_to_remove)

            # Resize workspace caches
            resize_workspace!(ws, n_nodes_new, n_cells_new, n_faces_new)

            # Resize the integrator state vector
            # This is the key DifferentialEquations.jl interface for size changes
            resize!(integrator, new_state_length)

            # Copy new state values into integrator
            integrator.u .= u_new

            # Force recomputation of derivatives
            u_modified!(integrator, true)

            @info "Merge complete: n_cells=$(n_cells_new), new state length=$(new_state_length)"
        end

    # Use DiscreteCallback - checked at end of each accepted step
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

"""
	get_depletion_callback_stats(cb::DiscreteCallback) -> NamedTuple

Extract statistics from a depletion callback.

This is a utility function to access the merge history from a depletion callback.
Note: This requires the callback to have been created with create_depletion_callback.

# Returns
NamedTuple with fields:
- `n_merges::Int`: Total number of cell merges
- `merge_times::Vector{Float64}`: Times at which merges occurred
- `merged_cells::Vector{Int}`: Original indices of merged cells
"""
function get_depletion_callback_stats(cb::DiscreteCallback)
    # Access the closure's captured state
    # This is somewhat fragile - depends on internal callback structure
    try
        # The affect! closure captures cb_state
        cb_state = cb.affect!.cb_state
        return (
            n_merges = cb_state.n_merges,
            merge_times = copy(cb_state.merge_times),
            merged_cells = copy(cb_state.merged_cells),
        )
    catch
        return (n_merges = -1, merge_times = Float64[], merged_cells = Int[])
    end
end

# ============================================================================
# StepUpdateCallback — mirrors the integrator's accepted state into the
# mutable mesh. The residual / Jacobian assembly is mesh-read-only; this
# callback runs at the end of each accepted step so AMR, depletion checks,
# and post-processing see a consistent mesh.
# ============================================================================

"""
	create_step_update_callback(ws, material) -> DiscreteCallback

Build a `DiscreteCallback` that fires after every accepted timestep and
mirrors the integrator's state vector `u` into the mutable mesh:

- Writes `u` cell temperatures and concentrations into `mesh.cell_states[i]`.
- Under ALE, also writes `u` node positions into `mesh.nodes[i].z`,
  refreshes `mesh.cross_section_area` from the χ block (when WithChi),
  and reruns `Geometry.update_geometry!` so derived geometry is in sync.

The callback dereferences `ws.layout` at every fire so AMR / depletion
callbacks that swap in a resized layout are picked up automatically.

The callback is unconditional (`condition` always true). With
`save_positions=(false, false)` it doesn't perturb the integrator's
output history.
"""
function create_step_update_callback(ws, material)
    mesh = ws.mesh

    condition = (u, t, integrator) -> true

    function affect!(integrator)
        u = integrator.u
        layout = ws.layout
        mesh.cross_section_area = cross_section_area_from_state(mesh, material, u, layout)
        update_mesh_from_state!(mesh, u, layout)
        return nothing
    end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end
