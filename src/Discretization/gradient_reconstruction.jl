"""
Gradient reconstruction for finite volume methods.

This module provides:
- Cell-centered gradient computation
- Face value interpolation (central, upwind)
- Second-order reconstruction with limiters
"""

# ============================================================================
# Cell-Centered Gradients
# ============================================================================

"""
	compute_gradients!(∇φ::Vector{Float64}, φ::Vector{Float64}, mesh::Mesh1D)

Compute cell-centered gradients via finite differences.

For interior cells: opposite-distance-weighted mean of the one-sided slopes
(`opposite_distance_gradient`) — second order on nonuniform meshes;
reduces to (φ[i+1] - φ[i-1]) / (2Δz) on uniform meshes.
For boundary cells: one-sided difference.

# Arguments
- `∇φ`: Output vector of gradients [units of φ / length]
- `φ`: Cell-centered field values
- `mesh`: 1D mesh with geometry
"""
function compute_gradients!(∇φ::AbstractVector, φ::AbstractVector, state, mesh)
    n_cells_val = length(state.T)

    @inbounds for i = 1:n_cells_val
        cell = mesh.cells[i]
        f_left, f_right = cell.face_ids

        face_left = mesh.faces[f_left]
        face_right = mesh.faces[f_right]

        left_id = face_left.cell_ids[1]
        right_id = face_right.cell_ids[2]

        if left_id == 0 && right_id == 0
            ∇φ[i] = zero(eltype(∇φ))
        elseif left_id == 0
            Δz = state.cell_centers[right_id] - state.cell_centers[i]
            ∇φ[i] = (φ[right_id] - φ[i]) / Δz
        elseif right_id == 0
            Δz = state.cell_centers[i] - state.cell_centers[left_id]
            ∇φ[i] = (φ[i] - φ[left_id]) / Δz
        else
            d_L = state.cell_centers[i] - state.cell_centers[left_id]
            d_R = state.cell_centers[right_id] - state.cell_centers[i]
            ∇φ[i] = opposite_distance_gradient(φ[left_id], φ[i], φ[right_id], d_L, d_R)
        end
    end

    return nothing
end

"""
	compute_temperature_gradients!(∇T, state, mesh)

Compute temperature gradients from the typed state cache.
"""
function compute_temperature_gradients!(∇T::AbstractVector, state, mesh)
    n_cells_val = length(state.T)

    @inbounds for i = 1:n_cells_val
        cell = mesh.cells[i]
        f_left, f_right = cell.face_ids

        face_left = mesh.faces[f_left]
        face_right = mesh.faces[f_right]

        left_id = face_left.cell_ids[1]
        right_id = face_right.cell_ids[2]

        T_i = state.T[i]

        if left_id == 0 && right_id == 0
            ∇T[i] = zero(eltype(∇T))
        elseif left_id == 0
            T_R = state.T[right_id]
            Δz = state.cell_centers[right_id] - state.cell_centers[i]
            ∇T[i] = (T_R - T_i) / Δz
        elseif right_id == 0
            T_L = state.T[left_id]
            Δz = state.cell_centers[i] - state.cell_centers[left_id]
            ∇T[i] = (T_i - T_L) / Δz
        else
            T_L = state.T[left_id]
            T_R = state.T[right_id]
            d_L = state.cell_centers[i] - state.cell_centers[left_id]
            d_R = state.cell_centers[right_id] - state.cell_centers[i]
            ∇T[i] = opposite_distance_gradient(T_L, T_i, T_R, d_L, d_R)
        end
    end

    return nothing
end


# ============================================================================
# Face Value Interpolation
# ============================================================================

"""
	face_value_central(φ_L::Float64, φ_R::Float64) -> Float64

Central (linear) interpolation for face value.
Assumes uniform cell spacing or values at cell centers.
"""
@inline face_value_central(φ_L::Float64, φ_R::Float64) = 0.5 * (φ_L + φ_R)

"""
	face_value_upwind(φ_L::Float64, φ_R::Float64, v::Float64) -> Float64

First-order upwind interpolation for face value.

# Arguments
- `φ_L, φ_R`: Field values in left and right cells
- `v`: Velocity at face (positive = flow from L to R)

# Returns
Face value based on upwind cell.
"""
@inline function face_value_upwind(φ_L::Float64, φ_R::Float64, v::Float64)
    return v >= 0.0 ? φ_L : φ_R
end

"""
	face_value_weighted(φ_L::Float64, φ_R::Float64, w_L::Float64) -> Float64

Distance-weighted interpolation for face value.

# Arguments
- `φ_L, φ_R`: Field values in left and right cells
- `w_L`: Weight for left cell (0 to 1, typically d_face_R / d_LR)
"""
@inline function face_value_weighted(φ_L::Float64, φ_R::Float64, w_L::Float64)
    return w_L * φ_L + (1.0 - w_L) * φ_R
end

# ============================================================================
# Second-Order Reconstruction with Limiters
# ============================================================================

"""
	compute_slope_ratio(φ_LL::Float64, φ_L::Float64, φ_R::Float64) -> Float64

Compute the slope ratio r for TVD limiters.

r = (φ_L - φ_LL) / (φ_R - φ_L)

For the face between L and R cells, using the upwind cell L and its upwind neighbor LL.
"""
@inline function compute_slope_ratio(φ_LL::Float64, φ_L::Float64, φ_R::Float64)
    denom = φ_R - φ_L
    if abs(denom) < 1e-20
        return 1.0  # Uniform field, no limiting needed
    end
    return (φ_L - φ_LL) / denom
end

"""
	face_value_limited(φ_LL::Float64, φ_L::Float64, φ_R::Float64, 
					   limiter::Function) -> Float64

Second-order face value with TVD limiter.

φ_f = φ_L + 0.5 * ψ(r) * (φ_R - φ_L)

where r = (φ_L - φ_LL) / (φ_R - φ_L) and ψ is the limiter function.

# Arguments
- `φ_LL`: Value at upwind-upwind cell
- `φ_L`: Value at upwind cell
- `φ_R`: Value at downwind cell
- `limiter`: Limiter function (minmod, vanleer, or superbee)
"""
@inline function face_value_limited(
    φ_LL::Float64,
    φ_L::Float64,
    φ_R::Float64,
    limiter::Function,
)
    r = compute_slope_ratio(φ_LL, φ_L, φ_R)
    ψ = limiter(r)
    return φ_L + 0.5 * ψ * (φ_R - φ_L)
end

"""
	face_value_second_order(φ_L::Float64, φ_R::Float64, ∇φ_L::Float64, 
							d_Lf::Float64, v::Float64) -> Float64

Second-order face value using gradient reconstruction.

For positive velocity (flow from L to R):
	φ_f = φ_L + ∇φ_L * d_Lf

# Arguments
- `φ_L, φ_R`: Cell-centered values
- `∇φ_L`: Gradient in upwind cell
- `d_Lf`: Distance from upwind cell center to face
- `v`: Velocity at face
"""
@inline function face_value_second_order(
    φ_L::Float64,
    φ_R::Float64,
    ∇φ_L::Float64,
    ∇φ_R::Float64,
    d_Lf::Float64,
    d_Rf::Float64,
    v::Float64,
)
    if v >= 0.0
        return φ_L + ∇φ_L * d_Lf
    else
        return φ_R - ∇φ_R * d_Rf
    end
end

# ============================================================================
# Face Gradient Computation
# ============================================================================

"""
	face_gradient(φ_L::Float64, φ_R::Float64, d_LR::Float64) -> Float64

Compute gradient at a face using two-point difference.
"""
@inline face_gradient(φ_L::Float64, φ_R::Float64, d_LR::Float64) = (φ_R - φ_L) / d_LR

"""
	compute_face_gradients!(∇φ_face::Vector{Float64}, φ::Vector{Float64}, mesh::Mesh1D)

Compute gradients at all faces.
"""
function compute_face_gradients!(∇φ_face::AbstractVector, φ::AbstractVector, state, mesh)
    n_faces_val = length(state.A_face)

    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            ∇φ_face[f] = zero(eltype(∇φ_face))
        else
            d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])
            ∇φ_face[f] = face_gradient(φ[left_id], φ[right_id], d_LR)
        end
    end

    return nothing
end
