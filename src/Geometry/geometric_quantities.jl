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
	cell_volume(mesh::Mesh1D, i::Int) -> Float64

Get the volume of the i-th cell.
"""
@inline cell_volume(mesh::Mesh1D, i::Int) = mesh.cell_volumes[i]

"""
	face_area(mesh::Mesh1D, f::Int) -> Float64

Get the area of the f-th face (cross-sectional area in 1D).
"""
@inline face_area(mesh::Mesh1D, f::Int) = mesh.face_areas[f]

"""
	cell_center(mesh::Mesh1D, i::Int) -> Float64

Get the center position of the i-th cell.
"""
@inline cell_center(mesh::Mesh1D, i::Int) = mesh.cell_centers[i]

"""
	face_position(mesh::Mesh1D, f::Int) -> Float64

Get the position of the f-th face.
"""
@inline face_position(mesh::Mesh1D, f::Int) = mesh.face_positions[f]

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
	total_volume(mesh::Mesh1D) -> Float64

Compute the total volume of the domain.
"""
function total_volume(mesh::Mesh1D)
    vol = 0.0
    @inbounds for i = 1:length(mesh.cell_volumes)
        if mesh.active_cell_mask[i]
            vol += mesh.cell_volumes[i]
        end
    end
    return vol
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
# Geometry Cache for ALE ODE Evaluation
# ============================================================================

"""
	GeometryCache

Pre-allocated arrays for geometry computation during ALE ODE evaluation.

This struct holds work arrays that are updated from the node positions in the
state vector at each ODE RHS evaluation. Using pre-allocated arrays avoids
allocations during the time-critical RHS computation.

# Fields
- `volumes::Vector{Float64}`: Cell volumes (length n_cells)
- `centers::Vector{Float64}`: Cell center positions (length n_cells)
- `face_positions::Vector{Float64}`: Face positions (length n_faces)
- `face_distances::Vector{Float64}`: Distances between cell centers across faces (length n_faces)

# Usage
```julia
cache = GeometryCache(n_cells, n_faces)
# Populate cache.volumes, cache.centers, etc. from the current state.
```
"""
struct GeometryCache
    volumes::Vector{Float64}
    centers::Vector{Float64}
    face_positions::Vector{Float64}
    face_distances::Vector{Float64}
end

"""
	GeometryCache(n_cells::Int, n_faces::Int)

Create a GeometryCache with pre-allocated arrays.
"""
function GeometryCache(n_cells::Int, n_faces::Int)
    GeometryCache(
        Vector{Float64}(undef, n_cells),
        Vector{Float64}(undef, n_cells),
        Vector{Float64}(undef, n_faces),
        Vector{Float64}(undef, n_faces),
    )
end

"""
	GeometryCache(mesh::Mesh1D)

Create a GeometryCache sized for the given mesh.
"""
GeometryCache(mesh::Mesh1D) = GeometryCache(length(mesh.cells), length(mesh.faces))

"""
	resize_geometry_cache!(cache::GeometryCache, n_cells::Int, n_faces::Int)

Resize a GeometryCache for a new mesh size.

This function is used during adaptive mesh operations (h-refinement, h-coarsening,
surface cell depletion) when the number of cells/faces changes.

# Arguments
- `cache::GeometryCache`: Cache to resize
- `n_cells::Int`: New number of cells
- `n_faces::Int`: New number of faces (typically n_cells + 1 for 1D)

# Notes
Arrays are resized in-place. Previous values are lost and should be
recomputed after calling this function.
"""
function resize_geometry_cache!(cache::GeometryCache, n_cells::Int, n_faces::Int)
    resize!(cache.volumes, n_cells)
    resize!(cache.centers, n_cells)
    resize!(cache.face_positions, n_faces)
    resize!(cache.face_distances, n_faces)
    return nothing
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
