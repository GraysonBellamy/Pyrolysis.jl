"""
Geometric Conservation Law (GCL) compliance verification.

The GCL ensures that mesh motion conserves volume correctly:
	dVᵢ/dt = Σ_faces(w·n·A)_f

This is CRITICAL for conservation - violations create spurious sources/sinks.
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
function check_gcl(mesh::Mesh1D, w::Vector{Float64}, dt::Float64; tol::Float64 = 1e-12)
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
@inline function gcl_volume_update(
    V_old::Float64,
    w_left::Float64,
    w_right::Float64,
    A::Float64,
    dt::Float64,
)
    return V_old + dt * (w_right - w_left) * A
end

"""
	update_volumes_from_geometry!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)

Update cell volumes from the analytic geometry `V = A · L` after advancing nodes
by `w · dt`. This is a *geometric* refresh of the volume cache, not a GCL flux
update; the previous name `update_volumes_gcl!` was misleading because the
function did not apply the GCL flux formula and silently used a static `A`
(ignoring time-varying cross-section). For the actual GCL diagnostic against the
flux formula, see `gcl_residual` below.
"""
function update_volumes_from_geometry!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)
    A = mesh.cross_section_area

    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        cell = mesh.cells[i]
        node_left, node_right = cell.node_ids

        z_left_new = mesh.nodes[node_left].z + dt * w[node_left]
        z_right_new = mesh.nodes[node_right].z + dt * w[node_right]
        L_new = abs(z_right_new - z_left_new)
        V_new = A * L_new

        mesh.cell_volumes[i] = V_new
        mesh.cell_states[i].volume = V_new
    end

    return nothing
end

"""
	gcl_residual(mesh::Mesh1D, V_old::Vector{Float64},
				 w::Vector{Float64}, dt::Float64;
				 A_old::Float64 = mesh.cross_section_area,
				 A_new::Float64 = mesh.cross_section_area) -> Float64

Discrete GCL residual against the post-motion mesh state, supporting time-varying
cross-section. The continuous GCL is `dV/dt = ∮(w·n)dA`; in 1D with cross-section
`A(t)` and faces moving at `w_left, w_right` it discretizes to:

	V_new - V_old ≈ dt · (w_right · A_new - w_left · A_new) + (A_new - A_old) · L_old

(the second term accounts for the area change at constant length over `dt`).
This helper returns the maximum relative residual across active cells. Pass
`A_old` and `A_new` explicitly when verifying a step that crossed an `A(t)` change.

For constant-A meshes (`A_old == A_new`), this collapses to the classical form
exercised by `gcl_error`.
"""
function gcl_residual(
    mesh::Mesh1D,
    V_old::Vector{Float64},
    w::Vector{Float64},
    dt::Float64;
    A_old::Float64 = mesh.cross_section_area,
    A_new::Float64 = mesh.cross_section_area,
)
    max_resid = 0.0

    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end
        cell = mesh.cells[i]
        node_left, node_right = cell.node_ids

        # Old cell length (volume / area is more robust than re-reading node z's).
        L_old = A_old > 0.0 ? V_old[i] / A_old : 0.0

        w_left = w[node_left]
        w_right = w[node_right]

        dV_expected = dt * (w_right - w_left) * A_new + (A_new - A_old) * L_old
        dV_actual = mesh.cell_volumes[i] - V_old[i]

        if abs(V_old[i]) > eps()
            rel_error = abs(dV_actual - dV_expected) / V_old[i]
            max_resid = max(max_resid, rel_error)
        end
    end

    return max_resid
end
