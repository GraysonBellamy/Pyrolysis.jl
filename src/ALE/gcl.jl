"""
Geometric Conservation Law (GCL) compliance verification.

The GCL ensures that mesh motion conserves volume correctly:
	dVᵢ/dt = Σ_faces(w·n·A)_f

This is CRITICAL for conservation - violations create spurious sources/sinks.

Reference: Design Document Section 3.4
"""

using ..Geometry: Mesh1D, cell_volume, face_area

# ============================================================================
# GCL Error Computation
# ============================================================================

"""
	gcl_error(mesh::Mesh1D, w::Vector{Float64}, dt::Float64) -> Float64

Compute the maximum relative GCL error across all cells.

The GCL requires that the volume change equals the flux of mesh velocity:
	V_new - V_old = dt * Σ_faces(w·n·A)

# Arguments
- `mesh`: Current mesh (must have cell_volumes computed)
- `w`: Mesh velocity at each node [m/s]
- `dt`: Time step [s]

# Returns
Maximum relative GCL error (should be ≈ machine precision for correct implementation)
"""
function gcl_error(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)
	max_error = 0.0
	A = mesh.cross_section_area

	@inbounds for i in eachindex(mesh.cells)
		if !mesh.active_cell_mask[i]
			continue
		end

		cell = mesh.cells[i]
		V_i = mesh.cell_volumes[i]

		# Get node indices for this cell
		node_left, node_right = cell.node_ids

		# Mesh velocity at cell faces (at nodes in 1D)
		w_left = w[node_left]
		w_right = w[node_right]

		# GCL flux: net volume flux through faces
		# In 1D: positive w on right face adds volume, positive w on left subtracts
		# Volume flux = w * A (w is velocity in z direction, n points outward)
		gcl_flux = (w_right - w_left) * A

		# Expected volume change from GCL
		dV_gcl = dt * gcl_flux

		# Actual volume change would be:
		# dV_actual = V_new - V_old
		# For verification, we compute what V_new would be:
		z_left_new = mesh.nodes[node_left].z + dt * w_left
		z_right_new = mesh.nodes[node_right].z + dt * w_right
		V_new = abs(z_right_new - z_left_new) * A
		dV_actual = V_new - V_i

		# Relative error
		if abs(V_i) > eps()
			rel_error = abs(dV_actual - dV_gcl) / V_i
			max_error = max(max_error, rel_error)
		end
	end

	return max_error
end

"""
	gcl_error_post_motion(mesh::Mesh1D, V_old::Vector{Float64}, 
						  w::Vector{Float64}, dt::Float64) -> Float64

Compute GCL error after mesh motion has been performed.

Compare actual volume change to expected from GCL.

# Arguments
- `mesh`: Mesh after motion (with updated volumes)
- `V_old`: Cell volumes before motion
- `w`: Mesh velocity that was used [m/s]
- `dt`: Time step that was used [s]

# Returns
Maximum relative GCL error
"""
function gcl_error_post_motion(mesh::Mesh1D, V_old::Vector{Float64},
	w::Vector{Float64}, dt::Float64)
	max_error = 0.0
	A = mesh.cross_section_area

	@inbounds for i in eachindex(mesh.cells)
		if !mesh.active_cell_mask[i]
			continue
		end

		cell = mesh.cells[i]
		node_left, node_right = cell.node_ids

		# Mesh velocity at faces
		w_left = w[node_left]
		w_right = w[node_right]

		# Expected volume change from GCL
		dV_gcl = dt * (w_right - w_left) * A

		# Actual volume change
		dV_actual = mesh.cell_volumes[i] - V_old[i]

		# Relative error
		if abs(V_old[i]) > eps()
			rel_error = abs(dV_actual - dV_gcl) / V_old[i]
			max_error = max(max_error, rel_error)
		end
	end

	return max_error
end

# ============================================================================
# GCL Check (Boolean)
# ============================================================================

"""
	check_gcl(mesh::Mesh1D, w::Vector{Float64}, dt::Float64; 
			  tol::Float64=1e-12) -> Bool

Check if mesh motion satisfies the Geometric Conservation Law.

# Arguments
- `mesh`: Current mesh
- `w`: Mesh velocity at each node [m/s]
- `dt`: Time step [s]
- `tol`: Tolerance for GCL satisfaction (default: 1e-12)

# Returns
`true` if GCL is satisfied to within tolerance
"""
function check_gcl(mesh::Mesh1D, w::Vector{Float64}, dt::Float64;
	tol::Float64 = 1e-12)
	error = gcl_error(mesh, w, dt)
	return error < tol
end

# ============================================================================
# GCL-Compliant Volume Update
# ============================================================================

"""
	gcl_volume_update(V_old::Float64, w_left::Float64, w_right::Float64,
					  A::Float64, dt::Float64) -> Float64

Compute new cell volume using GCL-compliant formula.

Instead of computing volume from geometry (which can have roundoff errors),
update volume using the GCL flux directly:
	V_new = V_old + dt * (w_right - w_left) * A

This ensures exact GCL satisfaction.

# Arguments
- `V_old`: Previous volume [m³]
- `w_left`: Mesh velocity at left face [m/s]
- `w_right`: Mesh velocity at right face [m/s]
- `A`: Face area [m²]
- `dt`: Time step [s]

# Returns
New cell volume [m³]
"""
@inline function gcl_volume_update(V_old::Float64, w_left::Float64, w_right::Float64,
	A::Float64, dt::Float64)
	return V_old + dt * (w_right - w_left) * A
end

"""
	update_volumes_gcl!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)

Update cell volumes using GCL-compliant formula.

This should be used instead of recomputing volumes from geometry
when strict GCL compliance is required.
"""
function update_volumes_gcl!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)
	A = mesh.cross_section_area

	@inbounds for i in eachindex(mesh.cells)
		if !mesh.active_cell_mask[i]
			continue
		end

		cell = mesh.cells[i]
		node_left, node_right = cell.node_ids

		V_old = mesh.cell_volumes[i]
		w_left = w[node_left]
		w_right = w[node_right]

		# GCL-compliant volume update
		V_new = gcl_volume_update(V_old, w_left, w_right, A, dt)

		mesh.cell_volumes[i] = V_new
		mesh.cell_states[i].volume = V_new
	end

	return nothing
end

# ============================================================================
# GCL Tracking and Enforcement
# ============================================================================

"""
    GCLTracker

Tracks cumulative GCL errors throughout simulation.

# Fields
- `max_error::Float64`: Maximum GCL error observed
- `cumulative_error::Float64`: Sum of |GCL errors| × dt (integrated error)
- `error_count::Int`: Number of timesteps with error > tolerance
- `tolerance::Float64`: Tolerance threshold for error counting
- `history::Vector{Tuple{Float64,Float64}}`: (time, error) history
"""
mutable struct GCLTracker
    max_error::Float64
    cumulative_error::Float64
    error_count::Int
    tolerance::Float64
    history::Vector{Tuple{Float64,Float64}}
end

"""
    GCLTracker(; tolerance::Float64=1e-10)

Create a new GCL tracker with specified tolerance.

# Example
```julia
gcl_tracker = GCLTracker(tolerance=1e-10)
```
"""
GCLTracker(; tolerance::Float64 = 1e-10) =
    GCLTracker(0.0, 0.0, 0, tolerance, Tuple{Float64,Float64}[])

"""
    update!(tracker::GCLTracker, t::Float64, gcl_error::Float64, dt::Float64)

Update tracker with current GCL error.

# Arguments
- `tracker`: GCL tracker to update
- `t`: Current simulation time [s]
- `gcl_error`: Current GCL error (from gcl_error or gcl_error_post_motion)
- `dt`: Time step [s]
"""
function update!(tracker::GCLTracker, t::Float64, gcl_error::Float64, dt::Float64)
    tracker.max_error = max(tracker.max_error, gcl_error)
    tracker.cumulative_error += abs(gcl_error) * dt

    if gcl_error > tracker.tolerance
        tracker.error_count += 1
    end

    push!(tracker.history, (t, gcl_error))

    return nothing
end

"""
    report(tracker::GCLTracker) -> String

Generate a report of GCL compliance over the simulation.

# Returns
String containing formatted GCL compliance report.

# Example
```julia
println(report(gcl_tracker))
```
"""
function report(tracker::GCLTracker)
    return """
    GCL Compliance Report
    =====================
    Maximum error:     $(tracker.max_error)
    Cumulative error:  $(tracker.cumulative_error)
    Violations:        $(tracker.error_count) timesteps exceeded tolerance $(tracker.tolerance)
    Total timesteps:   $(length(tracker.history))
    """
end

"""
    reset!(tracker::GCLTracker)

Reset tracker to initial state.
"""
function reset!(tracker::GCLTracker)
    tracker.max_error = 0.0
    tracker.cumulative_error = 0.0
    tracker.error_count = 0
    empty!(tracker.history)
    return nothing
end
