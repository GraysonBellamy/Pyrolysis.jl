"""
Heat conduction.
"""

"""
	conductive_flux(k_L, k_R, T_L, T_R, d_L, d_R) -> Real

Compute conductive heat flux at a face [W/m²] using distance-weighted harmonic mean.

	q = -k_f · (T_R − T_L) / (d_L + d_R)

where `k_f = (d_L + d_R) / (d_L/k_L + d_R/k_R)`. Scalar-type-generic.
"""
@inline function conductive_flux(k_L, k_R, T_L, T_R, d_L, d_R)
    k_face = harmonic_mean_weighted(k_L, k_R, d_L, d_R)
    return -k_face * (T_R - T_L) / (d_L + d_R)
end

"""
	conductive_flux(k_face, T_L, T_R, d_LR) -> Real

Compute conductive heat flux with pre-computed face conductivity.
"""
@inline conductive_flux(k_face, T_L, T_R, d_LR) = -k_face * (T_R - T_L) / d_LR

"""
	compute_conductive_fluxes!(q_cond, state, mesh, cache::PropertyCache)

Compute conductive heat fluxes at all faces, reading temperatures and
cell-center positions from the typed `state::StateCache` and face
connectivity from `mesh`.
"""
function compute_conductive_fluxes!(q_cond, state, mesh, cache::PropertyCache)
    n_faces_val = mesh.n_active_cells + 1

    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            q_cond[f] = zero(eltype(q_cond))
        else
            T_L = state.T[left_id]
            T_R = state.T[right_id]
            d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])

            q_cond[f] = conductive_flux(cache.k_face[f], T_L, T_R, d_LR)
        end
    end

    return nothing
end
