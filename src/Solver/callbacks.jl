"""
Solver callbacks.

ALE mesh motion is integrated directly into the unified workspace.
Thickness history is extracted from the state vector in extract_solution_ale!.
"""

using ..Adaptivity: find_merge_neighbor, park_cell!
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
    update_mesh_from_state!,
    mesh_mode,
    ALE
using ..Residual: cross_section_area_from_state
using ..Geometry: Mesh1D
using ..Materials: n_components, SOLID, is_gas
using ..Properties: MIXING_BLEND_LO, matrix_volume_fraction

# ============================================================================
# Cell Thickness Based Depletion Detection
# ============================================================================

"""
	is_cell_thin(u::AbstractVector, layout::StateLayout{NC},
				 cell_idx::Int, initial_thickness::Float64;
				 threshold::Float64=0.05) where NC -> Bool

Check if a cell has shrunk below a threshold fraction of its initial thickness.

This is the depletion trigger for ALE simulations: a cell that has shrunk to
`threshold` of its own initial thickness is parked (merged into a neighbour
and moved out of the live range) before its shrinking pitch drives the
implicit solve's conditioning (∝ 1/Δz²) beyond what the Newton can resolve.

# Arguments
- `u`: ODE state vector [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
- `layout`: StateLayout descriptor
- `cell_idx`: Index of the cell to check
- `initial_thickness`: Original thickness of **this** cell [m] — each cell is
  judged against its own reference
- `threshold`: Fraction of initial thickness below which cell is considered depleted

# Returns
- `true` if current cell thickness < threshold × initial_thickness
"""
function is_cell_thin(
    u::AbstractVector,
    layout::StateLayout{NC},
    cell_idx::Int,
    initial_thickness::Float64;
    threshold::Float64 = 0.05,
) where {NC}
    z = get_node_positions(u, layout)
    current_thickness = z[cell_idx+1] - z[cell_idx]
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
					threshold::Float64, n_live::Int) where NC -> Vector{Int}

Find the live cells (`1:n_live`) thinner than `threshold` × their **own**
initial thickness, sorted surface-first (descending z).

Each cell is compared against its own entry of the per-cell reference vector
`initial_thickness` (one entry per state slot, ghost entries unused). A
single shared reference would misclassify on any nonuniform mesh: with a
max/min cell-size ratio above `1/threshold`, a surface-cell reference either
flags the fine cells at t = 0 or never flags the coarse ones.

A live cell with non-positive thickness is a tangled mesh, not a depleted
cell, and raises an error instead of being parked.
"""
function find_thin_cells(
    u::AbstractVector,
    layout::StateLayout{NC},
    initial_thickness::AbstractVector{<:Real},
    threshold::Float64,
    n_live::Int,
) where {NC}
    length(initial_thickness) == layout.n_cells || throw(
        DimensionMismatch(
            "initial_thickness has $(length(initial_thickness)) entries for " *
            "$(layout.n_cells) state slots",
        ),
    )
    1 <= n_live <= layout.n_cells ||
        throw(ArgumentError("n_live = $(n_live) outside 1:$(layout.n_cells)"))
    thin = Int[]
    z = get_node_positions(u, layout)

    for i = 1:n_live
        thickness = z[i+1] - z[i]
        thickness > 0.0 || error(
            "Cell $i has non-positive thickness $(thickness) m: node positions " *
            "are no longer monotone in z (tangled mesh). This is a mesh-motion " *
            "failure, not depletion — refusing to park the inverted cell.",
        )
        if thickness < threshold * initial_thickness[i]
            push!(thin, i)
        end
    end

    # Surface first, so cascades of neighbouring thin cells park from the top.
    sort!(thin, by = i -> (z[i] + z[i+1]) / 2, rev = true)
    return thin
end

# ============================================================================
# Thickness Termination Callback (for ALE simulations)
# ============================================================================

"""
	create_thickness_termination_callback(sv_ref, min_thickness) -> ContinuousCallback

Create a callback that terminates the simulation when the live thickness
`z[n_live+1] − z[1]` falls below `min_thickness` [m].

`sv_ref` is either a `StateLayout` (every cell is live) or a
`DepletionCallbackState`, whose `n_live` tracks parked cells.
"""
function create_thickness_termination_callback(sv_ref, min_thickness::Float64)
    live_thickness = FunctionWrapper{Float64,Tuple{_StateVector}}(
        u -> begin
            layout, n_live = _layout_and_live(sv_ref)
            z = get_node_positions(u, layout)
            return z[n_live+1] - z[1]
        end,
    )
    condition = (u, t, integrator) -> live_thickness(u) - min_thickness
    affect! = (integrator) -> terminate!(integrator)
    return ContinuousCallback(condition, affect!)
end

_layout_and_live(layout::StateLayout) = (layout, layout.n_cells)

# ============================================================================
# Depletion Callback
# ============================================================================

"""
	DepletionCallbackState

Bookkeeping shared by the depletion callback, the thickness-termination
callback and solution extraction.

# Fields
- `layout::StateLayout`: the (fixed) state layout
- `n_live::Int`: number of live cells; cells `n_live+1:n` are parked ghosts
- `reference::Vector{Float64}`: per-slot depletion reference. Under ALE the
  cell's initial thickness [m]; on a fixed mesh the cell's initial condensed
  concentration Σ_{j∉gas} ξⱼ [kg/m³]. Shifts with the cells when an interior
  cell parks; ghost entries are unused.
- `park_times::Vector{Float64}`: time of each park event
- `parked_cells::Vector{Int}`: live index of the cell parked at each event
"""
mutable struct DepletionCallbackState{L<:StateLayout}
    layout::L
    n_live::Int
    reference::Vector{Float64}
    park_times::Vector{Float64}
    parked_cells::Vector{Int}
end

"""
	DepletionCallbackState(layout, mesh, material, u0)

Build the callback state for a fresh solve: every cell live, the reference
taken from the initial mesh (ALE: cell thickness) or the initial state
(fixed mesh: condensed concentration).
"""
function DepletionCallbackState(
    layout::StateLayout{NC},
    mesh::Mesh1D,
    material,
    u0::AbstractVector,
) where {NC}
    n = layout.n_cells
    length(mesh.cells) == n || throw(
        DimensionMismatch("mesh has $(length(mesh.cells)) cells for a $(n)-cell layout"),
    )
    reference = if mesh_mode(layout.mode) isa ALE
        mesh.cell_volumes ./ mesh.cross_section_area_initial
    else
        [_condensed_concentration(u0, layout, material, i) for i = 1:n]
    end
    return DepletionCallbackState(layout, n, reference, Float64[], Int[])
end

_layout_and_live(cb::DepletionCallbackState) = (cb.layout, cb.n_live)

"""
	n_live_cells(cb::DepletionCallbackState, t) -> Int

Number of live cells in the state saved at time `t`. Park events rewrite the
state after the saved values at the event time have been recorded, so a
state saved at exactly `t_park` still carries the pre-park layout.
"""
n_live_cells(cb::DepletionCallbackState, t::Real) =
    cb.layout.n_cells - count(<(t), cb.park_times)

# Σ_{j ∉ gas} ξⱼ of cell `i`: the condensed mass per unit volume, the fixed-mesh
# depletion measure.
@inline function _condensed_concentration(
    u::AbstractVector,
    layout::StateLayout{NC},
    material,
    i::Int,
) where {NC}
    s = 0.0
    @inbounds for j = 1:NC
        is_gas(material.components[j]) && continue
        s += get_concentration(u, layout, i, j)
    end
    return s
end

"""
	depletion_candidate(u, layout, mesh, material, cb::DepletionCallbackState,
	                    threshold) -> Int

The live cell to park next, or 0. ALE: the surface-most live cell thinner
than `threshold` × its initial thickness. Fixed mesh: the top live cell when
its condensed concentration has fallen below `threshold` × its initial value
(cells that started with no condensed mass never park).
"""
function depletion_candidate(
    u::AbstractVector,
    layout::StateLayout,
    mesh::Mesh1D,
    material,
    cb::DepletionCallbackState,
    threshold::Float64,
)
    n_live = cb.n_live
    if mesh_mode(layout.mode) isa ALE
        thin = find_thin_cells(u, layout, cb.reference, threshold, n_live)
        return isempty(thin) ? 0 : thin[1]
    end
    ref = cb.reference[n_live]
    ref > 0.0 || return 0
    return _condensed_concentration(u, layout, material, n_live) < threshold * ref ?
           n_live : 0
end

"""
	interior_void_cell(u, layout, material, cb::DepletionCallbackState) -> Int

Fixed mesh: the lowest live cell below the top live cell whose matrix volume
fraction ([`matrix_volume_fraction`](@ref)) has fallen below
`MIXING_BLEND_LO`, the property model's definition of a void, or 0. Cells
that started with no condensed mass are not tested.
"""
function interior_void_cell(
    u::AbstractVector,
    layout::StateLayout{NC},
    material,
    cb::DepletionCallbackState,
) where {NC}
    @inbounds for i = 1:(cb.n_live-1)
        cb.reference[i] > 0.0 || continue
        ξ = ntuple(j -> get_concentration(u, layout, i, j), Val(NC))
        φ_m = matrix_volume_fraction(ξ, get_temperature(u, layout, i), material)
        φ_m < MIXING_BLEND_LO && return i
    end
    return 0
end

"""
	park_cell_state!(u, layout, mesh, i_park, i_keep, χ_merged)

Write the post-park state vector after `park_cell!` has rearranged `mesh`:
T, ξ (and z under ALE) are read back from the mesh for every slot, so live
cells shift exactly as the mesh did and the ghost slot carries the merged
temperature, zero concentrations and the surface position. The χ block is
shifted the same way, with `χ_merged` in the keeper's slot and zero in the
ghost (the mesh holds no χ).
"""
function park_cell_state!(
    u::AbstractVector{Float64},
    layout::StateLayout{NC},
    mesh::Mesh1D,
    i_park::Int,
    i_keep::Int,
    χ_merged::Float64,
) where {NC}
    n = layout.n_cells
    @inbounds for i = 1:n
        st = mesh.cell_states[i]
        u[i] = st.T
        for j = 1:NC
            u[n+(j-1)*n+i] = st.ξ[j]
        end
    end
    ale = mesh_mode(layout.mode) isa ALE
    if ale
        zoff = z_offset(layout)
        @inbounds for i = 1:layout.n_nodes
            u[zoff+i] = mesh.nodes[i].z
        end
    end
    rng = chi_range(layout)
    if !isempty(rng)
        χoff = first(rng) - 1
        ghost = mesh.n_active_cells + 1
        if ale
            @inbounds for i = i_park:(ghost-1)
                u[χoff+i] = u[χoff+i+1]
            end
        end
        u[χoff+min(i_park, i_keep)] = χ_merged
        u[χoff+ghost] = 0.0
    end
    return u
end

"""
	park!(u, layout, mesh, material, cb::DepletionCallbackState, i_park, t) -> Int

Park live cell `i_park` at time `t`: choose the neighbour, merge
conservatively (`park_cell!`), rewrite `u`, and record the event. Returns the
new live-cell count.
"""
function park!(
    u::AbstractVector{Float64},
    layout::StateLayout,
    mesh::Mesh1D,
    material,
    cb::DepletionCallbackState,
    i_park::Int,
    t::Real,
)
    i_keep = find_merge_neighbor(mesh, i_park)
    i_keep > 0 || error("park!: cell $(i_park) has no live neighbour")
    ale = mesh_mode(layout.mode) isa ALE

    # χ merges by initial dry-solid mass (Σᵢ χᵢ·m_dry_init,i is the extensive
    # invariant); capture the weights before park_cell! sums them.
    χ_merged = 0.0
    if !isempty(chi_range(layout)) && !isempty(mesh.m_dry_solid_initial)
        m_p = mesh.m_dry_solid_initial[i_park]
        m_k = mesh.m_dry_solid_initial[i_keep]
        m_tot = m_p + m_k
        χ_merged =
            m_tot > 0.0 ?
            (m_p * get_chi(u, layout, i_park) + m_k * get_chi(u, layout, i_keep)) / m_tot :
            0.0
    end

    n_live = park_cell!(mesh, i_park, i_keep; material = material, ale = ale)
    park_cell_state!(u, layout, mesh, i_park, i_keep, χ_merged)

    if ale
        _shift_down!(cb.reference, i_park, n_live + 1)
    end
    cb.reference[n_live+1] = 0.0
    cb.n_live = n_live
    push!(cb.park_times, Float64(t))
    push!(cb.parked_cells, i_park)
    return n_live
end

@inline function _shift_down!(v::AbstractVector, from::Int, to::Int)
    @inbounds for i = from:(to-1)
        v[i] = v[i+1]
    end
    return v
end

"""
	DT_AFTER_PARK

Step size (s) proposed to the integrator after a park event.

The merged keeper cell is off the quasi-steady gas-concentration and
temperature profiles of its neighbours, so the residual at the parked state
is orders of magnitude larger than before the park and the integrator's
extrapolated predictor lands far from the solution. From this step the
Newton converges on the transient and the controller grows the step back
within a few steps. At the pre-park step size the Newton fails repeatedly
before the controller has cut the step, so the reset saves those failures at
every park.
"""
const DT_AFTER_PARK = 1e-6

# What a depletion check did; the callback then acts on the integrator.
@enum DepletionEvent NO_PARK CELL_PARKED STOP_SOLVE

"""
	create_depletion_callback(ws, mesh, material, cb_state::DepletionCallbackState;
							  threshold=0.05, min_cells=2) -> DiscreteCallback

Create the callback that parks depleted cells.

After every accepted step the condition scans the live cells
(`depletion_candidate`); when one qualifies the affect! parks it — at most
one cell per step, so the integrator adapts between events — or terminates
the solve when only `min_cells` live cells remain. The state vector keeps
its length; the park rewrites it in place, marks a derivative discontinuity
(`derivative_discontinuity!`) and proposes [`DT_AFTER_PARK`](@ref) for the
next step.

On a fixed mesh the condition also scans the cells below the top live cell
for a void ([`interior_void_cell`](@ref)). A void there means the material
has gasified from within. The model carries no heat across a void, so the
affect! warns and terminates the solve.

# Keyword Arguments
- `threshold::Float64=0.05`: ALE — fraction of a cell's initial thickness;
  fixed mesh — fraction of the top cell's initial condensed concentration.
- `min_cells::Int=2`: terminate instead of parking when `n_live ≤ min_cells`.
"""
function create_depletion_callback(
    ws,
    mesh::Mesh1D,
    material,
    cb_state::DepletionCallbackState;
    threshold::Float64 = 0.05,
    min_cells::Int = 2,
)
    layout = cb_state.layout
    fixed_mesh = !(mesh_mode(layout.mode) isa ALE)

    # The condition must keep firing at n_live ≤ min_cells: the affect!
    # terminates there, and returning false would let cells shrink without
    # bound.
    has_event = FunctionWrapper{Bool,Tuple{_StateVector}}(
        u ->
            depletion_candidate(u, layout, mesh, material, cb_state, threshold) != 0 ||
            fixed_mesh && interior_void_cell(u, layout, material, cb_state) != 0,
    )

    # Parks the candidate cell in `u` in place, or reports that the solve must stop.
    deplete! = FunctionWrapper{DepletionEvent,Tuple{_StateVector,Float64}}(
        (u, t) -> begin
            if fixed_mesh
                i_void = interior_void_cell(u, layout, material, cb_state)
                if i_void != 0
                    @warn "Cell $i_void, below the top live cell $(cb_state.n_live), is a void " *
                          "(matrix volume fraction < $MIXING_BLEND_LO) at t=$(round(t, digits=4)) s: " *
                          "the material has gasified from within. The model carries no heat " *
                          "across a void, so the fixed-mesh solve stops here. A regressing " *
                          "material belongs on the ALE mesh (use_ale = true)."
                    return STOP_SOLVE
                end
            end
            if cb_state.n_live <= min_cells
                @warn "At minimum cells ($min_cells), terminating simulation"
                return STOP_SOLVE
            end
            i_park =
                depletion_candidate(u, layout, mesh, material, cb_state, threshold)
            i_park == 0 && return NO_PARK
            n_live = park!(u, layout, mesh, material, cb_state, i_park, t)
            # The exposed boundary moved down one face; the workspace's
            # pre-cached boundary face id (shared with the Jacobian's AD
            # workspace) follows the mesh.
            ws.top_face_ids[1] = mesh.boundary_faces[1][1]
            @info "Parked cell $i_park at t=$(round(t, digits=4)) s; $(n_live) live cells"
            return CELL_PARKED
        end,
    )

    condition = (u, t, integrator) -> has_event(u)

    affect! = (integrator) -> begin
        # The park rewrites the integrator's own state, so it must receive
        # that vector, not a converted copy.
        event = deplete!(integrator.u::_StateVector, integrator.t)
        if event == STOP_SOLVE
            terminate!(integrator)
        elseif event == CELL_PARKED
            derivative_discontinuity!(integrator, true)
            set_proposed_dt!(integrator, DT_AFTER_PARK)
        end
        return nothing
    end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
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

The callback dereferences `ws.layout` at every fire so the experimental
h-AMR callbacks that swap in a resized layout are picked up automatically;
a park keeps the layout and only rewrites `u`.

The callback is unconditional (`condition` always true). With
`save_positions=(false, false)` it doesn't perturb the integrator's
output history.

The callback marks a derivative discontinuity after every accepted step, so
the integrator re-evaluates the derivative at the new state instead of
reusing the last Newton stage's, and treats its Jacobian as stale. The
residual does not read the mirrored fields, but on stiff fine meshes the
reused stage derivative and the longer Jacobian reuse cost accuracy far
beyond the one extra residual evaluation per step: on the finite-slab
heat-conduction case at 320 cells the L2 error grows 15× and the step count
2.6×.
"""
function create_step_update_callback(ws, material)
    mesh = ws.mesh

    sync_mesh! = FunctionWrapper{Nothing,Tuple{_StateVector}}(
        u -> begin
            layout = ws.layout
            mesh.cross_section_area =
                cross_section_area_from_state(mesh, material, u, layout)
            update_mesh_from_state!(mesh, u, layout)
        end,
    )

    condition = (u, t, integrator) -> true

    function affect!(integrator)
        sync_mesh!(integrator.u)
        derivative_discontinuity!(integrator, true)
        return nothing
    end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# ============================================================================
# Newton-failure step cap
# ============================================================================

"""
	FAILURE_MEMORY_GROWTH

Growth factor of the Newton-failure step cap: after a step whose Newton needed
retries, the proposed step stays at or below this multiple of that step
([`create_newton_failure_cap_callback`](@ref)).
"""
const FAILURE_MEMORY_GROWTH = 1.25

"""
	FAILURE_MEMORY_STEPS

Number of accepted steps over which the Newton-failure step cap holds
([`create_newton_failure_cap_callback`](@ref)).
"""
const FAILURE_MEMORY_STEPS = 20

"""
	create_newton_failure_cap_callback(; growth = FAILURE_MEMORY_GROWTH,
									   steps = FAILURE_MEMORY_STEPS) -> DiscreteCallback

Step cap that remembers Newton failures. After every accepted step the
callback reads the integrator's count of failed nonlinear solves. When the
step just accepted needed retries, the proposed next step is held at or below
`growth` times that step, and the cap stays in place for the following
`steps` accepted steps (a further failure restarts it from the new step).

The error controller regrows the step to its accuracy limit after a Newton
failure. On a strongly nonlinear front that limit lies beyond the range in
which the Newton converges, so without the cap the failure recurs every few
steps. A run without Newton failures never triggers the cap and keeps its
trajectory exactly.
"""
function create_newton_failure_cap_callback(;
    growth::Float64 = FAILURE_MEMORY_GROWTH,
    steps::Int = FAILURE_MEMORY_STEPS,
)
    n_failed = Ref(0)     # failed nonlinear solves up to the previous accepted step
    hold = Ref(0)         # accepted steps the cap still holds
    cap = Ref(Inf)

    affect! =
        integrator -> begin
            failures = integrator.stats.nnonlinconvfail
            if failures > n_failed[]
                n_failed[] = failures
                cap[] = growth * (integrator.t - integrator.tprev)
                hold[] = steps
            elseif hold[] > 0
                hold[] -= 1
                hold[] == 0 && (cap[] = Inf)
            end
            get_proposed_dt(integrator) > cap[] && set_proposed_dt!(integrator, cap[])
            derivative_discontinuity!(integrator, false)
            return nothing
        end

    return DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions = (false, false),
    )
end
