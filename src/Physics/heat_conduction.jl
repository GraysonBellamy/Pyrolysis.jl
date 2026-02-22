"""
Heat conduction.

Reference: Design Document Section 4.2
"""

"""
	conductive_flux(k_L::Float64, k_R::Float64, T_L::Float64, T_R::Float64,
					d_L::Float64, d_R::Float64) -> Float64

Compute conductive heat flux at a face [W/m²] using distance-weighted harmonic mean.

q = -k_f * (T_R - T_L) / (d_L + d_R)

where k_f = (d_L + d_R) / (d_L/k_L + d_R/k_R) is the distance-weighted harmonic mean.
This is the correct resistance-in-series interpolation for non-uniform meshes.
"""
@inline function conductive_flux(k_L::Float64, k_R::Float64,
	T_L::Float64, T_R::Float64,
	d_L::Float64, d_R::Float64)
	k_face = harmonic_mean_weighted(k_L, k_R, d_L, d_R)
	return -k_face * (T_R - T_L) / (d_L + d_R)
end

"""
	conductive_flux(k_face::Float64, T_L::Float64, T_R::Float64, d_LR::Float64) -> Float64

Compute conductive heat flux with pre-computed face conductivity.
"""
@inline function conductive_flux(k_face::Float64,
	T_L::Float64, T_R::Float64,
	d_LR::Float64)
	return -k_face * (T_R - T_L) / d_LR
end

"""
	compute_conductive_fluxes!(q_cond::Vector{Float64}, mesh::Mesh1D{NC}, cache::PropertyCache) where NC

Compute conductive heat fluxes at all faces.

# Arguments
- `q_cond`: Output vector of fluxes [W/m²], positive = flow in +z direction
- `mesh`: Mesh with current temperature state
- `cache`: Property cache with face conductivities
"""
function compute_conductive_fluxes!(q_cond::Vector{Float64}, mesh::Mesh1D{NC}, cache::PropertyCache) where NC
	n_faces_val = length(mesh.faces)

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: flux computed by boundary conditions
			q_cond[f] = 0.0
		else
			# Interior face
			T_L = mesh.cell_states[left_id].T
			T_R = mesh.cell_states[right_id].T
			d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])

			q_cond[f] = conductive_flux(cache.k_face[f], T_L, T_R, d_LR)
		end
	end

	return nothing
end
