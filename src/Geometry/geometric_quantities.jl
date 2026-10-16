"""
Geometric quantity computations.

Functions for computing and updating cached geometric quantities.
"""

# ============================================================================
# 1D Geometry Updates
# ============================================================================

"""
	update_geometry!(mesh::Mesh1D)

Recompute all cached geometric quantities from current node positions.

This should be called after any mesh motion.
"""
function update_geometry!(mesh::Mesh1D)
    n_cells_val = length(mesh.cells)
    n_faces_val = length(mesh.faces)
    A = mesh.cross_section_area

    # Update face positions and areas
    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        mesh.face_positions[f] = mesh.nodes[face.node_id].z
        mesh.face_areas[f] = A
    end

    # Update cell centers and volumes
    @inbounds for i = 1:n_cells_val
        cell = mesh.cells[i]
        z_left = mesh.nodes[cell.node_ids[1]].z
        z_right = mesh.nodes[cell.node_ids[2]].z

        mesh.cell_centers[i] = 0.5 * (z_left + z_right)
        mesh.cell_volumes[i] = abs(z_right - z_left) * A

        # Also update cell state volume
        mesh.cell_states[i].volume = mesh.cell_volumes[i]
    end

    return nothing
end

"""
	cell_size(mesh::Mesh1D, i::Int) -> Float64

Get the size (length) of the i-th cell.
"""
@inline function cell_size(mesh::Mesh1D, i::Int)
    cell = mesh.cells[i]
    z_left = mesh.nodes[cell.node_ids[1]].z
    z_right = mesh.nodes[cell.node_ids[2]].z
    return abs(z_right - z_left)
end

"""
	domain_bounds(mesh::Mesh1D) -> Tuple{Float64,Float64}

Get the domain bounds (z_min, z_max).
"""
function domain_bounds(mesh::Mesh1D)
    z_min = minimum(n.z for n in mesh.nodes)
    z_max = maximum(n.z for n in mesh.nodes)
    return (z_min, z_max)
end

# ============================================================================
# Nonuniform-mesh cell-centered gradient
# ============================================================================

"""
	opposite_distance_gradient(φ_L, φ_C, φ_R, d_L, d_R)

Second-order cell-centered gradient on a nonuniform 3-point stencil.
`d_L` / `d_R` are the center-to-center distances to the left / right
neighbors.

The one-sided slopes `s_L = (φ_C − φ_L)/d_L` and `s_R = (φ_R − φ_C)/d_R` are
second-order estimates of φ′ at the face-midpoints `z_C − d_L/2` and
`z_C + d_R/2`; interpolating them linearly back to `z_C` weights each slope
by the *opposite* distance:

	φ′(z_C) = (d_R·s_L + d_L·s_R) / (d_L + d_R) + O(h²)

The naive wide stencil `(φ_R − φ_L)/(d_L + d_R)` weights each slope by its
*own* distance and carries a first-order error `½(d_R − d_L)·φ″` on any
nonuniform mesh; the unweighted mean `(s_L + s_R)/2` halves that error
constant but is still first-order. All three coincide on uniform meshes.
"""
@inline function opposite_distance_gradient(φ_L, φ_C, φ_R, d_L, d_R)
    s_L = (φ_C - φ_L) / d_L
    s_R = (φ_R - φ_C) / d_R
    return (d_R * s_L + d_L * s_R) / (d_L + d_R)
end

"""
	opposite_distance_gradient_coeffs(d_L, d_R) -> (c_L, c_C, c_R)

Stencil coefficients of [`opposite_distance_gradient`](@ref):
`φ′(z_C) = c_L·φ_L + c_C·φ_C + c_R·φ_R`, with `c_L + c_C + c_R = 0`.
Used by the hand-derived structured-Jacobian overrides, where the
coefficients are exactly the `∂(φ′)/∂φ_*` entries.
"""
@inline function opposite_distance_gradient_coeffs(d_L, d_R)
    inv_sum = 1.0 / (d_L + d_R)
    c_L = -d_R * inv_sum / d_L
    c_R = d_L * inv_sum / d_R
    return (c_L, -(c_L + c_R), c_R)
end
