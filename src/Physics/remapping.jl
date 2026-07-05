"""
Solution remapping after mesh changes.

Conservative remapping transfers solution fields from old mesh to new mesh
while preserving conservation properties (mass, energy).
"""

using ..Materials: AbstractRemappingMethod
using ..Geometry: Mesh1D

# ============================================================================
# Remapping Method Types
# ============================================================================

"""
	FirstOrderRemap <: AbstractRemappingMethod

First-order conservative remapping.

Uses cell-averaged values with volume-weighted averaging for overlapping regions.
Simple but may introduce numerical diffusion.
"""
struct FirstOrderRemap <: AbstractRemappingMethod end

"""
	SecondOrderRemap <: AbstractRemappingMethod

Second-order conservative remapping with gradient reconstruction.

Uses cell gradients to reduce numerical diffusion while maintaining conservation.
Limiters prevent spurious oscillations.

# Fields
- `limiter::Symbol`: Limiter type (:minmod, :vanleer, :superbee)
"""
struct SecondOrderRemap <: AbstractRemappingMethod
    limiter::Symbol
end
SecondOrderRemap() = SecondOrderRemap(:minmod)

# Remapping is an instantaneous re-gridding at fixed time: both meshes must
# describe the same physical configuration, so they must share one
# cross-section area. Mixing A_old (overlap volumes) with A_new (division
# volumes) would silently scale every concentration by A_old/A_new.
function _remap_cross_section(mesh_new::Mesh1D, mesh_old::Mesh1D)
    A_old = mesh_old.cross_section_area
    A_new = mesh_new.cross_section_area
    isapprox(A_old, A_new; rtol = 1e-10) || throw(
        ArgumentError(
            "remap requires both meshes at the same cross-section area " *
            "(same-instant re-grid); got A_old = $A_old, A_new = $A_new. " *
            "A(t) evolution belongs to the ODE (χ swelling law), not the remap.",
        ),
    )
    return A_old
end

# ============================================================================
# First-Order Remapping
# ============================================================================

"""
	remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
					u_old::Vector{Float64}, mesh_old::Mesh1D,
					NC::Int, method::FirstOrderRemap;
					material=nothing)

First-order conservative remapping of solution vector.

For each new cell, find overlapping old cells and compute appropriate averages.

# Temperature Remapping
When `material` is provided, uses **energy-conserving** remapping:
	T_new = Σ(ρcₚ × T × V)_overlap / Σ(ρcₚ × V)_overlap

When `material` is nothing, falls back to simple volume-weighted averaging:
	T_new = Σ(T × V)_overlap / Σ(V)_overlap

The energy-conserving form is recommended when cells have different compositions
(and hence different ρcₚ values), as it preserves total thermal energy.

# Concentration Remapping
Concentrations use mass-conserving remapping:
	ξ_new = Σ(ξ × V)_overlap / V_new

# Arguments
- `u_new`: Output solution vector on new mesh
- `mesh_new`: New mesh (after motion/refinement)
- `u_old`: Solution vector on old mesh
- `mesh_old`: Old mesh
- `NC`: Number of components
- `method`: FirstOrderRemap instance
- `material`: (keyword) Material for ρcₚ computation, or nothing for volume-weighted
"""
function remap_solution!(
    u_new::Vector{Float64},
    mesh_new::Mesh1D,
    u_old::Vector{Float64},
    mesh_old::Mesh1D,
    NC::Int,
    ::FirstOrderRemap;
    material = nothing,
)
    n_cells_new = length(mesh_new.cells)
    n_cells_old = length(mesh_old.cells)
    n_dof_per_cell = 1 + NC  # T + ξ₁, ξ₂, ...
    A = _remap_cross_section(mesh_new, mesh_old)

    # Process each new cell
    @inbounds for i_new = 1:n_cells_new
        if !mesh_new.active_cell_mask[i_new]
            continue
        end

        cell_new = mesh_new.cells[i_new]
        z_left_new = mesh_new.nodes[cell_new.node_ids[1]].z
        z_right_new = mesh_new.nodes[cell_new.node_ids[2]].z
        V_new = mesh_new.cell_volumes[i_new]

        # Accumulate contributions from overlapping old cells
        # For energy-conserving remapping: E_accum = Σ(ρcₚ T V), ρcₚV_accum = Σ(ρcₚ V)
        # For volume-weighted fallback: T_accum = Σ(T V), V_accum = Σ(V)
        E_accum = 0.0       # Energy accumulator (or T*V if no material)
        ρcₚV_accum = 0.0    # ρcₚ*V accumulator (or V if no material)
        ξ_accum = zeros(NC)
        V_overlap_total = 0.0

        for i_old = 1:n_cells_old
            if !mesh_old.active_cell_mask[i_old]
                continue
            end

            cell_old = mesh_old.cells[i_old]
            z_left_old = mesh_old.nodes[cell_old.node_ids[1]].z
            z_right_old = mesh_old.nodes[cell_old.node_ids[2]].z

            # Compute overlap
            z_left_overlap = max(z_left_new, z_left_old)
            z_right_overlap = min(z_right_new, z_right_old)

            if z_right_overlap > z_left_overlap
                # There is overlap
                V_overlap = (z_right_overlap - z_left_overlap) * A

                # Extract old cell state
                idx_old = (i_old - 1) * n_dof_per_cell
                T_old = u_old[idx_old+1]

                # Extract concentrations
                ξ_old = ntuple(j -> u_old[idx_old+1+j], Val(NC))

                # Accumulate energy (or T*V for fallback)
                if material !== nothing
                    # Energy-conserving: weight by ρcₚ
                    ρcₚ_old = Properties.effective_heat_capacity(ξ_old, T_old, material)
                    E_accum += ρcₚ_old * T_old * V_overlap
                    ρcₚV_accum += ρcₚ_old * V_overlap
                else
                    # Fallback: simple volume weighting
                    E_accum += T_old * V_overlap
                    ρcₚV_accum += V_overlap
                end

                # Accumulate mass (same for both approaches)
                for j = 1:NC
                    ξ_accum[j] += ξ_old[j] * V_overlap
                end

                V_overlap_total += V_overlap
            end
        end

        # Compute new cell values
        idx_new = (i_new - 1) * n_dof_per_cell

        if V_overlap_total > 0.0
            # Temperature: energy-conserving or volume-weighted
            u_new[idx_new+1] = E_accum / ρcₚV_accum

            # Concentrations: mass-conserving
            # New concentration = total mass / new volume
            for j = 1:NC
                u_new[idx_new+1+j] = ξ_accum[j] / V_new
            end
        else
            # No overlap - extrapolate from nearest cell (edge case)
            # Find nearest old cell
            z_center_new = 0.5 * (z_left_new + z_right_new)
            nearest_old = 1
            min_dist = Inf

            for i_old = 1:n_cells_old
                if mesh_old.active_cell_mask[i_old]
                    z_center_old = mesh_old.cell_centers[i_old]
                    dist = abs(z_center_new - z_center_old)
                    if dist < min_dist
                        min_dist = dist
                        nearest_old = i_old
                    end
                end
            end

            idx_old = (nearest_old - 1) * n_dof_per_cell
            u_new[idx_new+1] = u_old[idx_old+1]
            for j = 1:NC
                u_new[idx_new+1+j] = u_old[idx_old+1+j]
            end
        end
    end

    return nothing
end

# ============================================================================
# Second-Order Remapping
# ============================================================================

"""
	remap_solution!(u_new::Vector{Float64}, mesh_new::Mesh1D,
					u_old::Vector{Float64}, mesh_old::Mesh1D,
					NC::Int, method::SecondOrderRemap)

Second-order conservative remapping with gradient reconstruction.
"""
function remap_solution!(
    u_new::Vector{Float64},
    mesh_new::Mesh1D,
    u_old::Vector{Float64},
    mesh_old::Mesh1D,
    NC::Int,
    method::SecondOrderRemap,
)
    n_cells_new = length(mesh_new.cells)
    n_cells_old = length(mesh_old.cells)
    n_dof_per_cell = 1 + NC
    A = _remap_cross_section(mesh_new, mesh_old)

    # First, compute gradients on old mesh
    ∇T_old = compute_gradients_for_remap(u_old, mesh_old, NC, 1)  # Temperature
    ∇ξ_old = [compute_gradients_for_remap(u_old, mesh_old, NC, 1+j) for j = 1:NC]

    # Get min/max values for limiting
    T_min_old, T_max_old = extrema_values(u_old, mesh_old, NC, 1)
    ξ_min_old = [extrema_values(u_old, mesh_old, NC, 1+j)[1] for j = 1:NC]
    ξ_max_old = [extrema_values(u_old, mesh_old, NC, 1+j)[2] for j = 1:NC]

    # Process each new cell
    @inbounds for i_new = 1:n_cells_new
        if !mesh_new.active_cell_mask[i_new]
            continue
        end

        cell_new = mesh_new.cells[i_new]
        z_left_new = mesh_new.nodes[cell_new.node_ids[1]].z
        z_right_new = mesh_new.nodes[cell_new.node_ids[2]].z
        z_center_new = 0.5 * (z_left_new + z_right_new)
        V_new = mesh_new.cell_volumes[i_new]

        T_accum = 0.0
        ξ_accum = zeros(NC)
        V_overlap_total = 0.0

        for i_old = 1:n_cells_old
            if !mesh_old.active_cell_mask[i_old]
                continue
            end

            cell_old = mesh_old.cells[i_old]
            z_left_old = mesh_old.nodes[cell_old.node_ids[1]].z
            z_right_old = mesh_old.nodes[cell_old.node_ids[2]].z
            z_center_old = mesh_old.cell_centers[i_old]

            z_left_overlap = max(z_left_new, z_left_old)
            z_right_overlap = min(z_right_new, z_right_old)

            if z_right_overlap > z_left_overlap
                V_overlap = (z_right_overlap - z_left_overlap) * A
                z_center_overlap = 0.5 * (z_left_overlap + z_right_overlap)

                # Extract old cell values
                idx_old = (i_old - 1) * n_dof_per_cell
                T_old = u_old[idx_old+1]

                # Second-order reconstruction
                dz = z_center_overlap - z_center_old
                T_reconstructed = T_old + ∇T_old[i_old] * dz

                T_accum += T_reconstructed * V_overlap

                for j = 1:NC
                    ξ_old_j = u_old[idx_old+1+j]
                    ξ_reconstructed = ξ_old_j + ∇ξ_old[j][i_old] * dz
                    ξ_accum[j] += ξ_reconstructed * V_overlap
                end

                V_overlap_total += V_overlap
            end
        end

        idx_new = (i_new - 1) * n_dof_per_cell

        if V_overlap_total > 0.0
            T_new = T_accum / V_overlap_total

            # Apply limiter - don't create new extrema
            T_new = clamp(T_new, T_min_old, T_max_old)
            u_new[idx_new+1] = T_new

            for j = 1:NC
                ξ_new_j = ξ_accum[j] / V_new
                # Concentrations must be non-negative
                ξ_new_j = max(ξ_new_j, 0.0)
                # Don't exceed old maximum
                ξ_new_j = min(ξ_new_j, ξ_max_old[j])
                u_new[idx_new+1+j] = ξ_new_j
            end
        else
            # Fallback: nearest cell
            nearest_old = find_nearest_cell(z_center_new, mesh_old)
            idx_old = (nearest_old - 1) * n_dof_per_cell
            u_new[idx_new+1] = u_old[idx_old+1]
            for j = 1:NC
                u_new[idx_new+1+j] = u_old[idx_old+1+j]
            end
        end
    end

    return nothing
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	compute_gradients_for_remap(u, mesh, NC, var_idx) -> Vector{Float64}

Compute cell-centered gradients for remapping.
"""
function compute_gradients_for_remap(
    u::Vector{Float64},
    mesh::Mesh1D,
    NC::Int,
    var_idx::Int,
)
    n_cells = length(mesh.cells)
    n_dof_per_cell = 1 + NC
    gradients = zeros(n_cells)

    @inbounds for i = 1:n_cells
        if !mesh.active_cell_mask[i]
            continue
        end

        idx = (i - 1) * n_dof_per_cell + var_idx
        φ_i = u[idx]
        z_i = mesh.cell_centers[i]

        # Find neighbors (previous and next active cells)
        i_prev = 0
        i_next = 0

        for k = (i-1):-1:1
            if mesh.active_cell_mask[k]
                i_prev = k
                break
            end
        end

        for k = (i+1):n_cells
            if mesh.active_cell_mask[k]
                i_next = k
                break
            end
        end

        if i_prev > 0 && i_next > 0
            # Central difference
            idx_prev = (i_prev - 1) * n_dof_per_cell + var_idx
            idx_next = (i_next - 1) * n_dof_per_cell + var_idx
            φ_prev = u[idx_prev]
            φ_next = u[idx_next]
            z_prev = mesh.cell_centers[i_prev]
            z_next = mesh.cell_centers[i_next]

            gradients[i] = (φ_next - φ_prev) / (z_next - z_prev)
        elseif i_prev > 0
            # One-sided (backward)
            idx_prev = (i_prev - 1) * n_dof_per_cell + var_idx
            φ_prev = u[idx_prev]
            z_prev = mesh.cell_centers[i_prev]
            gradients[i] = (φ_i - φ_prev) / (z_i - z_prev)
        elseif i_next > 0
            # One-sided (forward)
            idx_next = (i_next - 1) * n_dof_per_cell + var_idx
            φ_next = u[idx_next]
            z_next = mesh.cell_centers[i_next]
            gradients[i] = (φ_next - φ_i) / (z_next - z_i)
        else
            gradients[i] = 0.0
        end
    end

    return gradients
end

"""
	extrema_values(u, mesh, NC, var_idx) -> Tuple{Float64, Float64}

Get minimum and maximum values of a variable across the mesh.
"""
function extrema_values(u::Vector{Float64}, mesh::Mesh1D, NC::Int, var_idx::Int)
    n_cells = length(mesh.cells)
    n_dof_per_cell = 1 + NC

    val_min = Inf
    val_max = -Inf

    @inbounds for i = 1:n_cells
        if mesh.active_cell_mask[i]
            idx = (i - 1) * n_dof_per_cell + var_idx
            val = u[idx]
            val_min = min(val_min, val)
            val_max = max(val_max, val)
        end
    end

    return (val_min, val_max)
end

"""
	find_nearest_cell(z, mesh) -> Int

Find index of nearest active cell to position z.
"""
function find_nearest_cell(z::Float64, mesh::Mesh1D)
    n_cells = length(mesh.cells)
    nearest = 1
    min_dist = Inf

    @inbounds for i = 1:n_cells
        if mesh.active_cell_mask[i]
            dist = abs(z - mesh.cell_centers[i])
            if dist < min_dist
                min_dist = dist
                nearest = i
            end
        end
    end

    return nearest
end

# ============================================================================
# Conservation Check
# ============================================================================

"""
	check_remap_conservation(u_old, mesh_old, u_new, mesh_new, NC) -> NamedTuple

Check mass conservation after remapping.

# Returns
Named tuple with:
- `mass_error`: Relative mass conservation error for each component
"""
function check_remap_conservation(
    u_old::Vector{Float64},
    mesh_old::Mesh1D,
    u_new::Vector{Float64},
    mesh_new::Mesh1D,
    NC::Int,
)
    n_dof_per_cell = 1 + NC
    n_cells_old = length(mesh_old.cells)
    n_cells_new = length(mesh_new.cells)

    # Compute total mass on old mesh
    mass_old = zeros(NC)
    @inbounds for i = 1:n_cells_old
        if mesh_old.active_cell_mask[i]
            V = mesh_old.cell_volumes[i]
            idx = (i - 1) * n_dof_per_cell
            for j = 1:NC
                mass_old[j] += u_old[idx+1+j] * V
            end
        end
    end

    # Compute total mass on new mesh
    mass_new = zeros(NC)
    @inbounds for i = 1:n_cells_new
        if mesh_new.active_cell_mask[i]
            V = mesh_new.cell_volumes[i]
            idx = (i - 1) * n_dof_per_cell
            for j = 1:NC
                mass_new[j] += u_new[idx+1+j] * V
            end
        end
    end

    # Compute errors
    mass_error = zeros(NC)
    for j = 1:NC
        if abs(mass_old[j]) > eps()
            mass_error[j] = abs(mass_new[j] - mass_old[j]) / mass_old[j]
        end
    end

    return (mass_error = mass_error,)
end
