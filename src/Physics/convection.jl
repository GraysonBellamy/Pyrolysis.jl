"""
Gas convective heat transport (volumetric source form).

The non-conservative temperature equation for a porous solid with flowing gas:

    ρcₚ ∂T/∂t = -Σ_g c_g J_g ∂T/∂x + ∇·(k∇T) + Q_rad + Q_rxn

The gas convective term is a cell-centered volumetric source:

    S_conv = -Σ_g c_g(T) × J_g × ∂T/∂x

This is derived from the conservative energy equation ∂(ρcₚT)/∂t + ∇·(Σ c_g J T) = ...
via product rule expansion. The T·c_g·∇·J piece is handled implicitly by the ODE
solver's simultaneous T-ξ evolution; the explicit c_g·J·∇T piece gives S_conv.

S_conv is a pure redistribution term — it moves energy between cells but creates none.
Reaction heats (Q_rxn) are unchanged by this formulation.

Reference: ThermaKin Eq 17, Bird/Stewart/Lightfoot Ch 11
"""

"""
	compute_gas_convective_source!(S_conv::Vector{Float64},
								   J_gas::Matrix{Float64},
								   mesh::Mesh1D{NC},
								   material::Material{NC}) where NC

Compute gas convective heat source at cell centers [W/m³].

S_conv,i = -Σ_g c_g(T_i) × J̄_g,i × (∂T/∂x)_i

where:
- J̄_g,i = average of gas mass fluxes at bounding faces
- (∂T/∂x)_i = average of face gradients (central difference)
- Boundary cells use one-sided gradients

# Arguments
- `S_conv`: Output vector of volumetric sources [W/m³], one per cell
- `J_gas`: Gas mass fluxes J[j,f] [kg/(m²·s)] at faces
- `mesh`: Mesh with current state
- `material`: Material definition
"""
function compute_gas_convective_source!(S_conv::Vector{Float64},
	J_gas::Matrix{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC}) where NC

	n_cells_val = length(mesh.cells)

	@inbounds for i in 1:n_cells_val
		T_i = mesh.cell_states[i].T

		# Bounding faces via cell face_ids
		cell = mesh.cells[i]
		f_L, f_R = cell.face_ids

		# Neighbor cell IDs (0 = boundary)
		left_id = mesh.faces[f_L].cell_ids[1]
		right_id = mesh.faces[f_R].cell_ids[2]

		# Compute ∂T/∂x at cell center (distances from cell_centers)
		if left_id == 0 && right_id == 0
			# Single-cell mesh: no gradient available
			S_conv[i] = 0.0
			continue
		elseif left_id == 0
			# Left boundary cell: one-sided (forward) gradient
			T_R = mesh.cell_states[right_id].T
			d_R = abs(mesh.cell_centers[right_id] - mesh.cell_centers[i])
			dTdx = (T_R - T_i) / d_R
		elseif right_id == 0
			# Right boundary cell: one-sided (backward) gradient
			T_L = mesh.cell_states[left_id].T
			d_L = abs(mesh.cell_centers[i] - mesh.cell_centers[left_id])
			dTdx = (T_i - T_L) / d_L
		else
			# Interior cell: average of left and right face gradients
			T_L = mesh.cell_states[left_id].T
			T_R = mesh.cell_states[right_id].T
			d_L = abs(mesh.cell_centers[i] - mesh.cell_centers[left_id])
			d_R = abs(mesh.cell_centers[right_id] - mesh.cell_centers[i])
			dTdx = ((T_i - T_L) / d_L + (T_R - T_i) / d_R) * 0.5
		end

		# S_conv = -Σ_g c_g(T_i) × J̄_g × ∂T/∂x
		# Use compile-time unrolled version for type stability
		S_conv[i] = _gas_convective_source_unrolled(T_i, dTdx, J_gas, f_L, f_R, material)
	end

	return nothing
end

"""
Compile-time unrolled gas convective source computation.
Avoids runtime indexing into heterogeneous component tuple.
"""
@generated function _gas_convective_source_unrolled(T_i::Float64, dTdx::Float64,
	J_gas::Matrix{Float64}, f_L::Int, f_R::Int,
	material::Material{NC}) where NC

	exprs = Expr[]
	push!(exprs, :(S = 0.0))
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				J_avg = (J_gas[$j, f_L] + J_gas[$j, f_R]) * 0.5
				cp_j = material.components[$j].heat_capacity(T_i)
				S -= cp_j * J_avg * dTdx
			end
		end)
	end
	push!(exprs, :(return S))
	return Expr(:block, exprs...)
end
