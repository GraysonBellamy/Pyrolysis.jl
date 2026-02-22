"""
Solver callbacks.

Reference: Design Document Section 12.4

ALE mesh motion is integrated directly into the ODE function via PyrolysisODEFunctionALE.
Thickness history is extracted from the state vector in extract_solution_ale!.
"""

using DiffEqCallbacks: DiscreteCallback, PeriodicCallback, ContinuousCallback
using ..Adaptivity: find_surface_cell, find_interior_neighbor, find_exterior_neighbor,
	merge_adjacent_cells!, compact_mesh_after_merge!, find_merge_neighbor
using ..TimeIntegration: StateVectorALE, get_node_positions, get_concentration_ale,
	remove_surface_cell_from_state_vector, resize_ode_function_ale!,
	create_merged_state_vector
using ..Geometry: Mesh1D
using ..Core: n_components, SOLID

# ============================================================================
# Cell Thickness Based Depletion Detection
# ============================================================================

"""
	is_cell_thin(u::AbstractVector, sv::StateVectorALE{NC}, 
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
- `sv`: StateVectorALE descriptor
- `cell_idx`: Index of the cell to check
- `initial_thickness`: Original thickness of the cell [m]
- `threshold`: Fraction of initial thickness below which cell is considered depleted

# Returns
- `true` if current cell thickness < threshold × initial_thickness

# Example
If threshold=0.05 (5%), a cell that started at 0.4mm and has shrunk to 0.02mm
(5% of original) will trigger merging.
"""
function is_cell_thin(u::AbstractVector, sv::StateVectorALE{NC},
	cell_idx::Int, initial_thickness::Float64;
	threshold::Float64 = 0.05) where NC
	# Get node positions from state vector
	z = get_node_positions(u, sv)

	# Compute current cell thickness
	# Cell i spans from node i to node i+1
	current_thickness = z[cell_idx+1] - z[cell_idx]

	# Check if cell has shrunk below threshold
	return current_thickness < threshold * initial_thickness
end

"""
	get_cell_thickness(u::AbstractVector, sv::StateVectorALE, cell_idx::Int) -> Float64

Get the current thickness of a specific cell from the state vector.
"""
function get_cell_thickness(u::AbstractVector, sv::StateVectorALE, cell_idx::Int)
	z = get_node_positions(u, sv)
	return z[cell_idx+1] - z[cell_idx]
end

"""
	get_surface_cell_thickness(u::AbstractVector, sv::StateVectorALE) -> Float64

Get the current thickness of the surface cell (last cell, at z=L).
"""
function get_surface_cell_thickness(u::AbstractVector, sv::StateVectorALE)
	return get_cell_thickness(u, sv, sv.n_cells)
end

"""
	find_thin_cells(u::AbstractVector, sv::StateVectorALE{NC},
					initial_thickness::Float64, threshold::Float64) where NC -> Vector{Int}

Find all cells thinner than threshold fraction of initial thickness.

This function identifies cells that have shrunk below the depletion threshold
and returns their IDs sorted by z-position in descending order (surface cells first).

# Arguments
- `u`: ODE state vector [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
- `sv`: StateVectorALE descriptor
- `initial_thickness`: Original thickness of each cell [m]
- `threshold`: Fraction of initial thickness below which cell is considered thin

# Returns
Vector of cell IDs, sorted by z-position descending (surface first).

# Notes
Processing surface cells first (higher z) prevents "orphaning" interior cells
when multiple adjacent cells become thin simultaneously. This maintains mesh
connectivity during cascade scenarios.
"""
function find_thin_cells(u::AbstractVector, sv::StateVectorALE{NC},
	initial_thickness::Float64, threshold::Float64) where NC
	thin = Int[]
	z = get_node_positions(u, sv)

	for i in 1:sv.n_cells
		thickness = z[i+1] - z[i]
		if thickness < threshold * initial_thickness
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

This is used with PyrolysisODEFunctionALE where node positions are part of
the state vector.

# Arguments
- `sv_ref`: Either a StateVectorALE directly, or a DepletionCallbackState (for shared reference)
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
		# Get the current state vector descriptor (handles both direct sv and callback state)
		sv = _get_sv(sv_ref)
		z = TimeIntegration.get_node_positions(u, sv)
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
and the initial cell thickness for thickness-based depletion detection.

# Fields
- `sv::StateVectorALE{NC}`: Current state vector descriptor (updated after merges)
- `initial_cell_thickness::Float64`: Initial thickness of each cell [m]
- `n_merges::Int`: Total number of cell merges performed
- `merge_times::Vector{Float64}`: Times at which merges occurred
- `merged_cells::Vector{Int}`: Original indices of cells that were merged (for diagnostics)
"""
mutable struct DepletionCallbackState{NC}
	sv::StateVectorALE{NC}
	initial_cell_thickness::Float64
	n_merges::Int
	merge_times::Vector{Float64}
	merged_cells::Vector{Int}
end

function DepletionCallbackState(sv::StateVectorALE{NC}, initial_cell_thickness::Float64) where NC
	DepletionCallbackState{NC}(sv, initial_cell_thickness, 0, Float64[], Int[])
end

# Convenience constructor that computes initial thickness from mesh
function DepletionCallbackState(sv::StateVectorALE{NC}, mesh::Mesh1D) where NC
	# For uniform mesh, all cells have same initial thickness
	# Use the surface cell (last cell) as reference
	n_cells = length(mesh.cells)
	initial_thickness = mesh.cell_volumes[n_cells] / mesh.cross_section_area
	DepletionCallbackState{NC}(sv, initial_thickness, 0, Float64[], Int[])
end

# Helper to extract state vector from either direct sv or DepletionCallbackState
_get_sv(sv::TimeIntegration.StateVectorALE) = sv
_get_sv(cb_state::DepletionCallbackState) = cb_state.sv

"""
	create_depletion_callback(ode_fn, mesh, material, sv;
							  threshold=0.05, min_cells=2) -> DiscreteCallback

Create a callback that detects and handles thin cells.

This callback checks at each timestep whether ANY cell has shrunk below a
threshold fraction of its initial thickness. When this occurs, the callback:
1. Finds all thin cells, sorted by z-position (surface first)
2. Merges the first thin cell with its best neighbor
3. Resizes the ODE function's caches
4. Resizes the integrator's state vector using `resize!(integrator, ...)`
5. Remaps the solution to the new mesh

Only ONE merge is performed per callback invocation to allow the ODE solver
to adapt between merges.

This prevents numerical issues from extremely thin cells:
- Very small timesteps (CFL condition scales with Δx²)
- Extreme aspect ratios between cells
- Numerical instabilities

# Arguments
- `ode_fn`: The PyrolysisODEFunctionALE or wrapper being solved
- `mesh::Mesh1D`: The mesh (will be modified during merges)
- `material`: Material definition (used for energy-conserving merge)
- `sv::StateVectorALE`: Initial state vector descriptor

# Keyword Arguments
- `threshold::Float64=0.05`: Thickness fraction threshold for depletion detection.
  Cell is considered depleted when current_thickness < threshold × initial_thickness.
  Default 0.05 means cells are merged when shrunk to 5% of original thickness.
- `min_cells::Int=2`: Minimum number of cells before terminating instead of merging.
  When only `min_cells` remain, simulation terminates rather than merging.

# Returns
A `DiscreteCallback` that handles cell depletion.

# Notes
- The callback uses `resize!(integrator, new_length)` to inform the ODE solver
  of the state vector size change
- Works with both `PyrolysisODEFunctionALE` and `PyrolysisODEFunctionALEWithRadiation`
- The state vector descriptor in `ode_fn.sv` is automatically updated
- All cells (surface and interior) are monitored for depletion

# Example
```julia
ode_fn = PyrolysisODEFunctionALE(mesh, material, thermal_bcs, mass_bcs)
u0 = create_state_vector_ale(mesh)
sv = StateVectorALE(mesh)

# Merge cells when they shrink to 5% of original thickness
depletion_cb = create_depletion_callback(ode_fn, mesh, material, sv;
										 threshold=0.05, min_cells=3)

prob = ODEProblem(ode_fn, u0, tspan)
sol = solve(prob, Rodas5P(); callback=depletion_cb)
```
"""
function create_depletion_callback(ode_fn, mesh::Mesh1D, material, sv_or_state;
	threshold::Float64 = 0.05,
	min_cells::Int = 2)

	# Accept either StateVectorALE or DepletionCallbackState
	# If just a StateVectorALE, create callback state with initial thickness from mesh
	if sv_or_state isa DepletionCallbackState
		cb_state = sv_or_state
	else
		cb_state = DepletionCallbackState(sv_or_state, mesh)
	end
	NC = typeof(cb_state).parameters[1]  # Get NC from DepletionCallbackState{NC}

	# Condition function: check if ANY cell has become too thin
	condition = (u, t, integrator) -> begin
		# Use the current state vector descriptor
		current_sv = cb_state.sv
		n = current_sv.n_cells

		if n <= min_cells
			return false  # Don't merge if at minimum cells
		end

		# Find ALL thin cells (not just surface)
		thin_cells = find_thin_cells(u, current_sv, cb_state.initial_cell_thickness, threshold)
		return !isempty(thin_cells)
	end

	# Affect function: merge the first thin cell (surface-first order)
	affect! = (integrator) -> begin
		current_sv = cb_state.sv
		n = current_sv.n_cells
		t = integrator.t

		# Check minimum cells
		if n <= min_cells
			@warn "At minimum cells ($min_cells), terminating simulation"
			terminate!(integrator)
			return
		end

		# Find thin cells (sorted surface-first)
		thin_cells = find_thin_cells(integrator.u, current_sv,
			cb_state.initial_cell_thickness, threshold)

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
		current_thickness = get_cell_thickness(integrator.u, current_sv, cell_to_remove)
		thickness_percent = 100 * current_thickness / cb_state.initial_cell_thickness

		@info "Cell $cell_to_remove thin ($(round(thickness_percent, digits=1))% of initial) at t=$(round(t, digits=4))s, merging into cell $cell_to_keep, n_cells=$n"

		# Get cell volumes BEFORE mesh modification
		# These are needed for conservative state merging
		V_remove = mesh.cell_volumes[cell_to_remove]
		V_keep = mesh.cell_volumes[cell_to_keep]

		# Store current state for conservative merge
		u_old = copy(integrator.u)

		# Merge cells in mesh (marks cell_to_remove as inactive)
		merged_id = merge_adjacent_cells!(mesh, cell_to_remove, cell_to_keep; material = material)
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
		u_new, sv_new = create_merged_state_vector(u_old, current_sv,
			cell_to_remove, cell_to_keep, V_remove, V_keep;
			material = material)

		# Get dimensions from the new state vector descriptor (authoritative source)
		n_cells_new = sv_new.n_cells
		n_nodes_new = sv_new.n_nodes
		n_faces_new = n_cells_new + 1  # For 1D mesh
		new_state_length = length(u_new)

		# Update callback state
		cb_state.sv = sv_new
		cb_state.n_merges += 1
		push!(cb_state.merge_times, t)
		push!(cb_state.merged_cells, cell_to_remove)

		# Resize ODE function caches
		resize_ode_function_ale!(ode_fn, n_cells_new, n_faces_new)

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
		return (n_merges = cb_state.n_merges,
			merge_times = copy(cb_state.merge_times),
			merged_cells = copy(cb_state.merged_cells))
	catch
		return (n_merges = -1, merge_times = Float64[], merged_cells = Int[])
	end
end
