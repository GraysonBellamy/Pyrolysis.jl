# ============================================================================
# Reduced Analytical Jacobian for ALE Pyrolysis Solver
# ============================================================================
#
# Provides an approximate Jacobian that omits the dense ALE mesh velocity
# coupling (the lower-triangular C block, ~69% of non-zeros at 240 cells).
# The full monolithic RHS is still solved implicitly — Newton converges all
# terms to tolerance. Only the search direction is approximate, which means
# Newton may take slightly more iterations but each uses O(N) factorization
# instead of O(N^2.3).
#
# Contains:
#   - Physics derivatives: conduction, reactions, gas transport, BCs (Steps 4-7)
#   - ALE dilation correction: -ξ_j × θ (Steps 8a + 10)
#   - Zero for node position rows and z-column coupling
#
# Activated by solve(...; use_reduced_jacobian=true).
# Remove by deleting this file and its include/export lines.
# ============================================================================

"""
	jacobian_sparsity_ale_reduced(n_cells::Int, NC::Int;
		radiation_model::RadiationModel = NO_RADIATION)

Compute the reduced sparsity pattern for the ALE Jacobian.

Omits the dense lower-triangular ALE mesh velocity coupling (C block) and
geometry corrections (∂physics/∂z), keeping only:
- Physics-Physics block (T,ξ): block-tridiagonal (same as fixed-mesh)
- z diagonal: required for W = I/(γdt) - J identity term (values are zero)

Returns a SparseMatrixCSC{Float64, Int64} with the sparsity pattern.
"""
function jacobian_sparsity_ale_reduced(n_cells::Int, NC::Int;
	radiation_model::RadiationModel = NO_RADIATION)

	n_nodes = n_cells + 1
	n_physics = n_cells * (1 + NC)
	n_total = n_physics + n_nodes

	# Estimate non-zeros: physics block + z diagonal
	nnz_estimate = n_cells * 3 * (1 + NC)^2 + n_nodes
	if has_nonlocal_jacobian_coupling(radiation_model)
		nnz_estimate += n_cells * n_cells * NC ÷ 2
	end

	I = Int[]
	J = Int[]
	sizehint!(I, nnz_estimate)
	sizehint!(J, nnz_estimate)

	# ==========================================================================
	# Physics-Physics coupling (identical to fixed-mesh jacobian_sparsity)
	# ==========================================================================
	for i in 1:n_cells
		T_row = i

		# Temperature depends on neighbors' temperatures
		for neighbor in max(1, i-1):min(n_cells, i+1)
			push!(I, T_row)
			push!(J, neighbor)
		end

		# Temperature depends on neighbor concentrations
		for neighbor in max(1, i-1):min(n_cells, i+1)
			for j in 1:NC
				ξ_col = n_cells + (j-1)*n_cells + neighbor
				push!(I, T_row)
				push!(J, ξ_col)
			end
		end

		# Concentration rows for cell i
		for j in 1:NC
			ξ_row = n_cells + (j-1)*n_cells + i

			# Concentration depends on neighbor temperatures
			for neighbor in max(1, i-1):min(n_cells, i+1)
				push!(I, ξ_row)
				push!(J, neighbor)
			end

			# Concentration depends on neighbor concentrations
			for neighbor in max(1, i-1):min(n_cells, i+1)
				for k in 1:NC
					ξ_col = n_cells + (k-1)*n_cells + neighbor
					push!(I, ξ_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# ==========================================================================
	# Beer-Lambert radiation coupling (if active)
	# ==========================================================================
	if has_nonlocal_jacobian_coupling(radiation_model)
		for i in 1:n_cells
			T_row = i
			for k in (i+1):n_cells
				for j in 1:NC
					ξ_col = n_cells + (j-1)*n_cells + k
					push!(I, T_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# ==========================================================================
	# z-block diagonal (required for W = I/(γdt) - J)
	# The actual Jacobian values here are zero, but the sparsity entry must
	# exist so the solver can add the identity term to the diagonal.
	# ==========================================================================
	for k in 1:n_nodes
		z_idx = n_physics + k
		push!(I, z_idx)
		push!(J, z_idx)
	end

	return sparse(I, J, ones(length(I)), n_total, n_total)
end


"""
	compute_analytical_jacobian_ale_reduced!(J, u, p, t, f, jac_cache)

Compute the reduced analytical Jacobian for ALE simulations.

Computes physics derivatives (conduction, reactions, transport, BCs) and ALE
dilation correction, but omits the expensive ALE mesh velocity coupling:
- Node velocity derivatives (dense lower-triangular C block — 69% of nnz)
- ALE advection Jacobian
- Geometry corrections (∂physics/∂z)
- Node position equation Jacobian

The full monolithic RHS is still solved implicitly by Newton. This Jacobian
provides the search direction — Newton converges all terms to tolerance
regardless of the omitted entries.

This is Steps 0-7 + 8a + 10 from compute_analytical_jacobian_ale!, without
Steps 8b, 8c, 9, 11, 11b-11e, 12.
"""
function compute_analytical_jacobian_ale_reduced!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionALE{M, Mat, BCS, PC, RC, SV, GC},
	jac_cache::JacobianCacheALE{NC, NR}) where {NC, NR, M, Mat, BCS, PC, RC,
	SV <: StateVectorALE{NC}, GC}

	mesh = f.mesh
	material = f.material
	prop_cache = f.prop_cache
	sv = f.sv
	n = sv.n_cells
	base_cache = jac_cache.base

	# =========================================================================
	# Step 0: Initialize static caches (only on first call)
	# =========================================================================
	initialize!(base_cache, material)

	if base_cache.sparsity_map === nothing && J isa SparseMatrixCSC
		base_cache.sparsity_map = create_sparsity_map(J)
	end
	smap = base_cache.sparsity_map

	# =========================================================================
	# Step 1: Unpack state and update geometry/properties
	# =========================================================================
	unpack_state_ale!(mesh, u, sv)
	Properties.update_properties!(prop_cache, mesh, material)

	zero_jacobian!(J)

	# =========================================================================
	# Step 2: Compute current RHS values (needed for property derivative terms)
	# =========================================================================
	Discretization.compute_rhs!(f.dT_work, f.dξ_work, mesh, prop_cache, f.rhs_cache,
		material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
		I_radiation = f.I_radiation, radiation_model = f.radiation_model)

	# =========================================================================
	# Step 3: Compute θ (volumetric strain rate) for dilation correction
	# =========================================================================
	_compute_ale_mesh_velocity!(f.w_work, mesh, prop_cache, material,
		t, f.min_density;
		θ_out = f.θ_work,
		precomputed_dξ_rxn = f.rhs_cache.dξ_rxn,
		dξ_dt = f.dξ_dt_work)

	# =========================================================================
	# Step 4: Cell-local Jacobian terms (reactions, properties)
	# =========================================================================
	S_matrix = jac_cache.S_matrix
	fill!(S_matrix, 0.0)

	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		cache_all_component_properties!(base_cache, T, material)

		ρcp = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp

		evaluate_reactions!(f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, material)

		evaluate_reaction_derivatives!(base_cache.dr_dT, base_cache.dr_dξ,
			f.rhs_cache.rates, ξ, T, base_cache)

		dQ_rxn_dξ = base_cache.work_NC
		dQ_rxn_dT = heat_source_derivatives_cached!(dQ_rxn_dξ,
			base_cache.dr_dT, base_cache.dr_dξ,
			f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, base_cache, material)

		species_source_derivatives_cached!(@view(base_cache.dS_dT[:, i]), @view(base_cache.dS_dξ[:, :, i]),
			base_cache.dr_dT, base_cache.dr_dξ, base_cache, NC)

		for j in 1:NC
			S_matrix[j, i] = f.rhs_cache.dξ_rxn[j, i]
		end

		dρcp_dξ = base_cache.work_NC_2
		dρcp_dT = heat_capacity_derivatives_cached!(dρcp_dξ, ξ,
			base_cache.comp_heat_capacities, base_cache.comp_dcp_dT)
		dinv_ρcp_dT = inverse_heat_capacity_derivatives!(dρcp_dξ, ρcp, dρcp_dT)

		RHS_T = f.dT_work[i] * ρcp

		# Energy equation derivatives
		T_row = i
		fast_sparse_add!(J, smap, T_row, T_row, dinv_ρcp_dT * RHS_T)
		fast_sparse_add!(J, smap, T_row, T_row, inv_ρcp * dQ_rxn_dT)

		for j in 1:NC
			ξ_col = n + (j - 1) * n + i
			fast_sparse_add!(J, smap, T_row, ξ_col, dρcp_dξ[j] * RHS_T)
		end

		for j in 1:NC
			ξ_col = n + (j - 1) * n + i
			fast_sparse_add!(J, smap, T_row, ξ_col, inv_ρcp * dQ_rxn_dξ[j])
		end

		# Species equation derivatives
		for j in 1:NC
			ξ_row = n + (j - 1) * n + i
			fast_sparse_write!(J, smap, ξ_row, T_row, base_cache.dS_dT[j, i])

			for k in 1:NC
				ξ_col = n + (k - 1) * n + i
				fast_sparse_write!(J, smap, ξ_row, ξ_col, base_cache.dS_dξ[j, k, i])
			end
		end
	end

	# =========================================================================
	# Step 5: Conduction and gas transport Jacobian
	# =========================================================================
	precompute_all_cell_properties!(base_cache, mesh, material)

	compute_all_face_jacobians_fused!(J, smap, mesh, prop_cache, material,
		f.rhs_cache.J_gas, n, Val(NC), base_cache)

	compute_gas_convective_source_jacobian!(J, smap, mesh, prop_cache, material,
		f.rhs_cache.J_gas, n, Val(NC), base_cache)

	# =========================================================================
	# Step 5b: Darcy advective transport Jacobian (if active)
	# =========================================================================
	if has_darcy_flow(material)
		compute_darcy_transport_jacobian!(J, smap, mesh, prop_cache, material, n, Val(NC), jac_cache)
	end

	# =========================================================================
	# Step 6: Boundary condition Jacobian
	# =========================================================================
	compute_bc_jacobian!(J, smap, mesh, prop_cache, f.bc_set, n, NC, t,
		f.top_face_ids, f.bottom_face_ids, f.radiation_model)
	compute_mass_bc_jacobian!(J, smap, mesh, prop_cache, f.bc_set, material, n, NC, t,
		f.top_face_ids, f.bottom_face_ids, base_cache.is_gas, base_cache)

	# =========================================================================
	# Step 7: Radiation Jacobian (if active)
	# =========================================================================
	if f.I_radiation > 0.0 && f.radiation_model == BEER_LAMBERT
		compute_radiation_jacobian!(J, smap, mesh, prop_cache, f.I_radiation, n, NC,
			base_cache.comp_absorptions, base_cache.τ_cumulative)
	end

	# =========================================================================
	# Step 8a: Strain rate derivatives (needed for dilation correction)
	# =========================================================================
	compute_strain_rate_derivatives!(jac_cache.dθ_dT, jac_cache.dθ_dξ,
		base_cache.dS_dT, base_cache.dS_dξ, S_matrix, mesh, material, jac_cache)

	# =========================================================================
	# Step 10: ALE dilation correction Jacobian for mesh-bound species
	# =========================================================================
	_add_dilation_correction_jacobian!(J, smap, mesh, base_cache.is_gas, f.θ_work,
		jac_cache.dθ_dT, jac_cache.dθ_dξ, n, NC, sv.z_offset)

	# =========================================================================
	# OMITTED (approximate — these terms are small relative to physics):
	# - Step 8b: Node velocity derivatives (dense lower-triangular C block)
	# - Step 8c: Geometry derivatives
	# - Step 9:  ALE advection Jacobian
	# - Step 11: Geometry corrections for conduction/BCs/transport
	# - Step 12: Node position equation Jacobian (dz/dt = w)
	#
	# The z rows of J are all zero — the solver adds I/(γdt) to the diagonal,
	# making the z block trivially invertible.
	# =========================================================================

	return nothing
end
