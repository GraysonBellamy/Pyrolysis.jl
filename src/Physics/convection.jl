"""
Convective heat transfer (gas convection).

Reference: Design Document Section 4.5
"""

"""
	convective_heat_flux_advection(J::Float64, cp_gas::Float64, T_L::Float64, T_R::Float64) -> Float64

Compute convective heat flux using advection form [W/m²].

q_conv = J × cp × (T_R - T_L)

This redistributes heat internally without creating/destroying energy.
The advection form avoids double-counting energy already accounted for
in reaction enthalpies.

# Arguments
- `J`: Gas mass flux [kg/(m²·s)], positive = toward surface (higher z)
- `cp_gas`: Gas heat capacity [J/(kg·K)]
- `T_L`: Temperature of left (lower z) cell [K]
- `T_R`: Temperature of right (higher z) cell [K]

# Reference
ThermaKin equations 17-18
"""
@inline function convective_heat_flux_advection(J::Float64, cp_gas::Float64,
	T_L::Float64, T_R::Float64)
	return J * cp_gas * (T_R - T_L)
end

"""
	compute_convective_fluxes!(q_conv::Vector{Float64}, J::Matrix{Float64},
							   mesh::Mesh1D{NC}, material::Material{NC}) where NC

Compute convective heat fluxes at all faces.

# Arguments
- `q_conv`: Output vector of fluxes [W/m²]
- `J`: Gas mass fluxes J[j,f] [kg/(m²·s)]
- `mesh`: Mesh with current state
- `material`: Material definition
"""
function compute_convective_fluxes!(q_conv::Vector{Float64}, J::Matrix{Float64},
	mesh::Mesh1D{NC}, material::Material{NC}) where NC
	n_faces_val = length(mesh.faces)

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: advection form has no boundary contribution
			# Energy leaving system is handled by boundary conditions
			q_conv[f] = 0.0
			continue
		end

		# Interior face: advection form q = J × cp × (T_R - T_L)
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		T_avg = (T_L + T_R) * 0.5
		ΔT = T_R - T_L

		q = 0.0
		for j in 1:NC
			comp = material.components[j]
			if comp.phase == GAS
				J_j = J[j, f]
				cp_j = comp.heat_capacity(T_avg)
				q += convective_heat_flux_advection(J_j, cp_j, T_L, T_R)
			end
		end

		q_conv[f] = q
	end

	return nothing
end
# ============================================================================
# Darcy-Based Convective Heat Fluxes (Advection Form)
# ============================================================================
#
# Using advection form, the gas redistributes heat between cells:
#     q_conv = Σⱼ Jⱼ × cp,j × (T_R - T_L)
#
# This avoids double-counting energy already in reaction enthalpies.
# Boundary faces have q_conv = 0 (energy leaving system handled by BCs).
#
# Reference: ThermaKin equations 17-18
# ============================================================================

"""
	compute_convective_heat_fluxes_darcy!(q_conv::Vector{Float64},
										  mesh::Mesh1D{NC},
										  material::Material{NC},
										  v_gas::Vector{Float64},
										  w_node::Vector{Float64}) where NC

Compute convective heat fluxes at all faces using Darcy gas velocity in ALE frame.

This computes the enthalpy transport by advecting gas:
	q_conv = Σⱼ ξⱼ × (v_gas - w_mesh) × cp,j × T

# Arguments
- `q_conv`: Output vector of heat fluxes [W/m²]
- `mesh`: Mesh with current state
- `material`: Material definition
- `v_gas`: Darcy gas velocity at faces [m/s]
- `w_node`: Mesh velocity at nodes [m/s]

# Notes
- Only gas components contribute to convective transport
- Solid and liquid enthalpies are transported by mesh motion in ALE frame
- Uses first-order upwinding based on relative velocity direction
"""
function compute_convective_heat_fluxes_darcy!(q_conv::Vector{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC},
	v_gas::Vector{Float64},
	w_node::Vector{Float64}) where NC
	n_faces_val = length(mesh.faces)

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: handled by boundary conditions
			q_conv[f] = 0.0
			continue
		end

		# Get mesh velocity at face
		node_id = face.node_id
		w_face = w_node[node_id]

		# Darcy velocity at face
		v_g = v_gas[f]

		# Relative velocity in ALE frame
		v_rel = v_g - w_face

		# Cell temperatures and states
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Advection form: q = v_rel × ξ × cp × (T_R - T_L)
		T_avg = (T_L + T_R) * 0.5
		ΔT = T_R - T_L

		# Sum contributions from all gas components
		q = 0.0
		for j in 1:NC
			comp = material.components[j]
			if comp.phase == GAS
				# Upwind concentration (still needed for mass flux)
				ξ_upwind = v_rel >= 0.0 ? ξ_L[j] : ξ_R[j]

				# Gas heat capacity at average temperature
				cp_j = comp.heat_capacity(T_avg)

				# Convective heat flux: advection form redistributes heat
				q += v_rel * ξ_upwind * cp_j * ΔT
			end
		end

		q_conv[f] = q
	end

	return nothing
end

"""
	compute_convective_heat_fluxes_from_mass_fluxes!(q_conv::Vector{Float64},
													 J::Matrix{Float64},
													 mesh::Mesh1D{NC},
													 material::Material{NC}) where NC

Compute convective heat fluxes from gas mass fluxes.

This is an alternative to `compute_convective_heat_fluxes_darcy!` that uses
pre-computed mass fluxes (which may include diffusive contributions):
	q_conv = Σⱼ Jⱼ × cp,j × T_face

# Arguments
- `q_conv`: Output vector of heat fluxes [W/m²]
- `J`: Gas mass fluxes J[j,f] [kg/(m²·s)]
- `mesh`: Mesh with current state
- `material`: Material definition

# Notes
This formulation is more general as J may include both advective
and diffusive mass transport from Darcy-Fick.
Uses dispatch-based unrolling to avoid allocation from heterogeneous tuple indexing.

# Boundary Handling
At boundary faces, gas escaping carries enthalpy with it. The convective heat flux
is computed using the boundary cell temperature:
- Top boundary (right_id=0): gas escapes at T_surface, q_conv = Σⱼ Jⱼ × cₚ,ⱼ × T_surface
- Bottom boundary (left_id=0): typically impermeable, but if gas escapes, uses T_bottom

This ensures energy conservation: the heat carried by escaping gas is properly
subtracted from the surface cell's energy balance.
"""
function compute_convective_heat_fluxes_from_mass_fluxes!(q_conv::Vector{Float64},
	J::Matrix{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC}) where NC
	n_faces_val = length(mesh.faces)
	comps = material.components

	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: advection form has no boundary contribution
			# Energy leaving system is handled by boundary conditions.
			# Reaction enthalpies already account for gas creation energy.
			q_conv[f] = 0.0
		else
			# Interior face: advection form q = J × cp × (T_R - T_L)
			T_L = mesh.cell_states[left_id].T
			T_R = mesh.cell_states[right_id].T
			q_conv[f] = _sum_convective_heat_flux(J, f, T_L, T_R, comps, Val(NC))
		end
	end

	return nothing
end

# ============================================================================
# Unrolled helpers for interior convective heat flux summation
# ============================================================================

# Unrolled helpers for convective heat flux summation (advection form)
# q = J × cp(T_avg) × (T_R - T_L)

# Single component
@inline function _sum_convective_heat_flux(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::Tuple{C1}, ::Val{1}) where {C1}
	ΔT = T_R - T_L
	T_avg = (T_L + T_R) * 0.5
	q = 0.0
	c1 = comps[1]
	if c1.phase == GAS
		q += J[1, f] * c1.heat_capacity(T_avg) * ΔT
	end
	return q
end

# Two components
@inline function _sum_convective_heat_flux(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	ΔT = T_R - T_L
	T_avg = (T_L + T_R) * 0.5
	q = 0.0
	c1 = comps[1]
	if c1.phase == GAS
		q += J[1, f] * c1.heat_capacity(T_avg) * ΔT
	end
	c2 = comps[2]
	if c2.phase == GAS
		q += J[2, f] * c2.heat_capacity(T_avg) * ΔT
	end
	return q
end

# Three components
@inline function _sum_convective_heat_flux(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	ΔT = T_R - T_L
	T_avg = (T_L + T_R) * 0.5
	q = 0.0
	c1 = comps[1]
	if c1.phase == GAS
		q += J[1, f] * c1.heat_capacity(T_avg) * ΔT
	end
	c2 = comps[2]
	if c2.phase == GAS
		q += J[2, f] * c2.heat_capacity(T_avg) * ΔT
	end
	c3 = comps[3]
	if c3.phase == GAS
		q += J[3, f] * c3.heat_capacity(T_avg) * ΔT
	end
	return q
end

# Four components (common case)
@inline function _sum_convective_heat_flux(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::Tuple{C1, C2, C3, C4}, ::Val{4}) where {C1, C2, C3, C4}
	ΔT = T_R - T_L
	T_avg = (T_L + T_R) * 0.5
	q = 0.0
	c1 = comps[1]
	if c1.phase == GAS
		q += J[1, f] * c1.heat_capacity(T_avg) * ΔT
	end
	c2 = comps[2]
	if c2.phase == GAS
		q += J[2, f] * c2.heat_capacity(T_avg) * ΔT
	end
	c3 = comps[3]
	if c3.phase == GAS
		q += J[3, f] * c3.heat_capacity(T_avg) * ΔT
	end
	c4 = comps[4]
	if c4.phase == GAS
		q += J[4, f] * c4.heat_capacity(T_avg) * ΔT
	end
	return q
end

# Generic fallback using generated function for NC > 4
@inline function _sum_convective_heat_flux(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::Tuple, ::Val{NC}) where NC
	return _sum_convective_heat_flux_gen(J, f, T_L, T_R, comps)
end

@generated function _sum_convective_heat_flux_gen(J::Matrix{Float64}, f::Int, T_L::Float64, T_R::Float64,
	comps::NTuple{NC, Any}) where NC
	expr = quote
		ΔT = T_R - T_L
		T_avg = (T_L + T_R) * 0.5
		q = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			if c.phase == GAS
				q += J[$i, f] * c.heat_capacity(T_avg) * ΔT
			end
		end)
	end
	push!(expr.args, :(return q))
	return expr
end