"""
Boundary temperature solver for cell-centered finite volume method.

This module implements Newton iteration to solve for the surface temperature T_s
that satisfies energy balance at boundaries:

	q_BC(T_s) = q_cond(T_s, T_cell)

where:
- q_BC = q_ext - εσT_s⁴ - h(T_s - T_∞)  [applied heat flux, positive into domain]
- q_cond = k(T_s - T_cell) / Δz_half    [conduction into material, positive when T_s > T_cell]

This is necessary because in cell-centered FVM, we only know cell-center temperatures.
The surface temperature must be computed to correctly evaluate temperature-dependent
BC fluxes (especially radiation which scales as T⁴).

Reference: docs/boundary_temperature_enhancement.md
"""

using ..Core: STEFAN_BOLTZMANN, evaluate_param, evaluate_param_T
using ..Core: RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY

# ============================================================================
# Main Solver
# ============================================================================

"""
	solve_boundary_temperature(T_cell, k_eff, Δz_half, bc, t, ε_eff; kwargs...)

Solve for surface temperature satisfying energy balance using Newton iteration.

Finds T_s such that:
	q_BC(T_s) = k_eff × (T_s - T_cell) / Δz_half

Both sides are positive when heat flows INTO the material.

# Arguments
- `T_cell::Float64`: Temperature at adjacent cell center [K]
- `k_eff::Float64`: Effective thermal conductivity [W/(m·K)]
- `Δz_half::Float64`: Distance from boundary to cell center [m]
- `bc::AbstractThermalBC`: Thermal boundary condition
- `t::Float64`: Current time [s]
- `ε_eff::Float64`: Effective surface emissivity from material properties

# Keyword Arguments
- `T_guess::Float64`: Initial guess (default: T_cell)
- `tol::Float64`: Convergence tolerance (default: 1e-10)
- `max_iter::Int`: Maximum iterations (default: 50)
- `relaxation::Float64`: Under-relaxation factor (default: 1.0)

# Returns
- `BoundaryState`: Solved boundary state with T_surface, fluxes, convergence info
"""
function solve_boundary_temperature(
	T_cell::Float64,
	k_eff::Float64,
	Δz_half::Float64,
	bc::AbstractThermalBC,
	t::Float64,
	ε_eff::Float64;
	radiation_model::RadiationModel = SURFACE_ABSORPTION,
	T_guess::Float64 = T_cell,
	tol::Float64 = 1e-10,
	max_iter::Int = 50,
	relaxation::Float64 = 1.0,
)
	# Newton iteration
	T_s = T_guess
	converged = false
	iter = 0

	for i in 1:max_iter
		iter = i

		# Compute BC flux and its derivative at current T_s
		q_bc, dq_bc_dT = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

		# Conduction flux: q_cond = k(T_s - T_cell) / Δz_half
		# Positive when T_s > T_cell (heat flows from surface into material)
		q_cond = k_eff * (T_s - T_cell) / Δz_half
		dq_cond_dT = k_eff / Δz_half

		# Residual: f(T_s) = q_bc - q_cond = 0
		# At solution: heat applied = heat conducted
		f = q_bc - q_cond
		df = dq_bc_dT - dq_cond_dT

		# Check convergence (relative to flux magnitude)
		scale = max(abs(q_bc), abs(q_cond), 1.0)
		if abs(f) < tol * scale
			converged = true
			break
		end

		# Avoid division by zero (shouldn't happen with reasonable BCs)
		if abs(df) < 1e-30
			break
		end

		# Newton update with relaxation
		ΔT = -f / df
		T_s = T_s + relaxation * ΔT

		# Safety bounds to prevent numerical instability
		T_s = clamp(T_s, 200.0, 3000.0)
	end

	# Final flux evaluation for storage
	q_bc_final, _ = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)
	q_cond_final = k_eff * (T_s - T_cell) / Δz_half

	return BoundaryState(T_s, q_bc_final, q_cond_final, T_cell, k_eff, Δz_half,
		converged, iter)
end

"""
	solve_boundary_temperature!(state::BoundaryState, T_cell, k_eff, Δz_half, bc, t, ε_eff; kwargs...)

In-place version that updates an existing BoundaryState. Avoids allocation.
"""
@inline function solve_boundary_temperature!(
	state::BoundaryState,
	T_cell::Float64,
	k_eff::Float64,
	Δz_half::Float64,
	bc::AbstractThermalBC,
	t::Float64,
	ε_eff::Float64;
	radiation_model::RadiationModel = SURFACE_ABSORPTION,
	tol::Float64 = 1e-10,
	max_iter::Int = 50,
	relaxation::Float64 = 1.0,
)
	# Newton iteration - use previous T_surface as initial guess
	T_s = state.T_surface
	converged = false
	iter = 0

	for i in 1:max_iter
		iter = i

		# Compute BC flux and its derivative at current T_s
		q_bc, dq_bc_dT = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

		# Conduction flux: q_cond = k(T_s - T_cell) / Δz_half
		# Positive when T_s > T_cell (heat flows from surface into material)
		q_cond = k_eff * (T_s - T_cell) / Δz_half
		dq_cond_dT = k_eff / Δz_half

		# Residual: f(T_s) = q_bc - q_cond = 0
		# At solution: heat applied = heat conducted
		f = q_bc - q_cond
		df = dq_bc_dT - dq_cond_dT

		# Check convergence (relative to flux magnitude)
		scale = max(abs(q_bc), abs(q_cond), 1.0)
		if abs(f) < tol * scale
			converged = true
			break
		end

		# Avoid division by zero (shouldn't happen with reasonable BCs)
		if abs(df) < 1e-30
			break
		end

		# Newton update with relaxation
		ΔT = -f / df
		T_s = T_s + relaxation * ΔT

		# Safety bounds to prevent numerical instability
		T_s = clamp(T_s, 200.0, 3000.0)
	end

	# Final flux evaluation and update state in-place
	q_bc_final, _ = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)
	q_cond_final = k_eff * (T_s - T_cell) / Δz_half

	# Update fields in-place (no allocation)
	state.T_surface = T_s
	state.q_applied = q_bc_final
	state.q_conducted = q_cond_final
	state.T_cell = T_cell
	state.k_eff = k_eff
	state.Δz_half = Δz_half
	state.converged = converged
	state.iterations = iter

	return nothing
end

# ============================================================================
# BC Flux and Derivative Evaluation
# ============================================================================

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff)

Evaluate boundary condition flux and its temperature derivative.

# Arguments
- `bc::AbstractThermalBC`: Boundary condition
- `T_s::Float64`: Surface temperature [K]
- `t::Float64`: Current time [s]
- `ε_eff::Float64`: Effective emissivity from material properties

# Returns
- `(q, dq_dT)`: Flux [W/m²] and derivative [W/(m²·K)], positive into domain
"""
function evaluate_bc_flux_and_derivative end

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

Radiation-model aware variant used by the boundary temperature solve.

This is required to keep `RadiativeFluxBC` consistent with the main RHS, where
incident radiation is typically handled *volumetrically* (e.g. SURFACE_ABSORPTION)
and should therefore contribute **zero** boundary flux here.
"""
@inline function evaluate_bc_flux_and_derivative(
	bc::AbstractThermalBC,
	T_s::Float64,
	t::Float64,
	ε_eff::Float64,
	radiation_model::RadiationModel,
)
	return evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff)
end

# ----------------------------------------------------------------------------
# AdiabaticBC: q = 0
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(::AdiabaticBC, T_s::Float64, t::Float64, ε_eff::Float64)
	return 0.0, 0.0
end

# ----------------------------------------------------------------------------
# HeatFluxBC: q = prescribed flux (no T_s dependence)
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::HeatFluxBC{<:Real}, T_s::Float64, t::Float64, ε_eff::Float64)
	return Float64(bc.flux), 0.0
end

@inline function evaluate_bc_flux_and_derivative(bc::HeatFluxBC{<:Function}, T_s::Float64, t::Float64, ε_eff::Float64)
	return bc.flux(t), 0.0
end

# ----------------------------------------------------------------------------
# ConvectiveBC: q = h(T_inf - T_s)
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::ConvectiveBC, T_s::Float64, t::Float64, ε_eff::Float64)
	h = evaluate_param(bc.h, t)
	T_inf = evaluate_param(bc.T_inf, t)

	q = h * (T_inf - T_s)
	dq_dT = -h
	return q, dq_dT
end

# ----------------------------------------------------------------------------
# RadiativeBC: q = F × ε × σ × (T_inf⁴ - T_s⁴)
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::RadiativeBC, T_s::Float64, t::Float64, ε_eff::Float64)
	# Use effective emissivity if not specified
	ε = if bc.ε === nothing
		ε_eff
	else
		evaluate_param_T(bc.ε, T_s)
	end
	T_inf = evaluate_param(bc.T_inf, t)
	F = evaluate_param(bc.F, t)

	q = F * ε * STEFAN_BOLTZMANN * (T_inf^4 - T_s^4)
	dq_dT = -4 * F * ε * STEFAN_BOLTZMANN * T_s^3

	# Note: If ε depends on T, there's an additional term:
	# dq_dT += F * dε_dT * σ * (T_inf^4 - T_s^4)
	# This is usually small and we neglect it for simplicity

	return q, dq_dT
end

# ----------------------------------------------------------------------------
# RadiativeFluxBC: q = α × incident_flux (for NO_RADIATION mode)
# For in-depth radiation models, this is handled volumetrically
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::RadiativeFluxBC, T_s::Float64, t::Float64, ε_eff::Float64)
	# This BC represents incident radiation - absorptivity doesn't depend on T_s
	# If radiation is handled volumetrically, flux should be 0 here
	# For now, we return the absorbed flux (consistent with NO_RADIATION mode)
	flux = evaluate_param(bc.flux, t)
	α = if bc.absorptivity === nothing
		ε_eff  # Kirchhoff's law: α = ε
	else
		evaluate_param(bc.absorptivity, t)
	end
	return flux * α, 0.0
end

# Radiation-model aware handling for incident radiation.
# - NO_RADIATION/SURFACE_ABSORPTION: treat as boundary heat flux for correct T_surface
# - BEER_LAMBERT/P1: handled volumetrically, return 0 to avoid double-counting
@inline function evaluate_bc_flux_and_derivative(
	bc::RadiativeFluxBC,
	T_s::Float64,
	t::Float64,
	ε_eff::Float64,
	radiation_model::RadiationModel,
)
	if radiation_model == NO_RADIATION || radiation_model == SURFACE_ABSORPTION
		# Radiation absorbed at boundary - include in boundary balance
		# This ensures T_surface is physically correct (hotter than T_cell when net heating)
		flux = evaluate_param(bc.flux, t)
		α = if bc.absorptivity === nothing
			ε_eff
		else
			evaluate_param(bc.absorptivity, t)
		end
		return flux * α, 0.0
	end

	# Volumetric radiation models (BEER_LAMBERT, P1) handle incident
	# radiation separately, so we must not include it in the boundary balance.
	return 0.0, 0.0
end

# ----------------------------------------------------------------------------
# DirichletThermalBC: T = T_bc (special handling, not typical flux BC)
# For boundary temperature solve, we return the prescribed temperature directly
# since T_s = T_bc by definition
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::DirichletThermalBC, T_s::Float64, t::Float64, ε_eff::Float64)
	# Dirichlet BC prescribes temperature directly
	# The "flux" formulation: infinite conductance enforces T_s = T_bc
	# We handle this by setting large h with T_inf = T_bc
	T_bc = evaluate_param(bc.T, t)

	# Use a very large heat transfer coefficient to enforce T_s ≈ T_bc
	# But not so large that it causes numerical issues
	h_large = 1e8  # W/(m²·K) - reduced from 1e10 for stability
	q = h_large * (T_bc - T_s)
	dq_dT = -h_large

	return q, dq_dT
end

# ----------------------------------------------------------------------------
# CombinedThermalBC: Sum of component fluxes
# ----------------------------------------------------------------------------
# Use regular function instead of @generated to avoid excessive compilation time
@generated function evaluate_bc_flux_and_derivative(bc::CombinedThermalBC{BCs}, T_s::Float64, t::Float64, ε_eff::Float64) where {BCs}
	N = fieldcount(BCs)
	exprs = Expr[]
	push!(exprs, :(q_total = 0.0))
	push!(exprs, :(dq_total = 0.0))
	for i in 1:N
		push!(exprs, quote
			q_i, dq_i = evaluate_bc_flux_and_derivative(bc.bcs[$i], T_s, t, ε_eff)
			q_total += q_i
			dq_total += dq_i
		end)
	end
	push!(exprs, :(return q_total, dq_total))
	return Expr(:block, exprs...)
end

@generated function evaluate_bc_flux_and_derivative(
	bc::CombinedThermalBC{BCs},
	T_s::Float64,
	t::Float64,
	ε_eff::Float64,
	radiation_model::RadiationModel,
) where {BCs}
	N = fieldcount(BCs)
	exprs = Expr[]
	push!(exprs, :(q_total = 0.0))
	push!(exprs, :(dq_total = 0.0))
	for i in 1:N
		push!(exprs, quote
			q_i, dq_i = evaluate_bc_flux_and_derivative(bc.bcs[$i], T_s, t, ε_eff, radiation_model)
			q_total += q_i
			dq_total += dq_i
		end)
	end
	push!(exprs, :(return q_total, dq_total))
	return Expr(:block, exprs...)
end

# ----------------------------------------------------------------------------
# FlameBC: Deferred implementation (complex state-dependent BC)
# For now, return 0 - flame heat flux will be added separately
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(bc::FlameBC, T_s::Float64, t::Float64, ε_eff::Float64)
	# FlameBC depends on mass flux state, not just temperature
	# The flame contribution is handled separately in the RHS computation
	# Return 0 here to avoid double-counting
	return 0.0, 0.0
end

# ============================================================================
# Boundary State Update Function
# ============================================================================

"""
	update_boundary_state!(mesh, boundary_id, bc, cache, t)

Update the boundary state for a given boundary by solving for surface temperature.

This should be called at the beginning of each RHS evaluation to compute the
correct surface temperature and resulting fluxes.

# Arguments
- `mesh::Mesh1D`: The mesh containing boundary state storage
- `boundary_id::Int`: Boundary identifier (1 = top/surface, 2 = bottom)
- `bc::AbstractThermalBC`: Thermal boundary condition for this boundary
- `cache::PropertyCache`: Cached material properties (for k, ε)
- `t::Float64`: Current time [s]

# Returns
- `BoundaryState`: The updated boundary state

# Note
Uses the previous T_surface as initial guess for faster convergence.
"""
function update_boundary_state!(
	mesh::Mesh1D,
	boundary_id::Int,
	bc::AbstractThermalBC,
	cache::PropertyCache,
	t::Float64,
	radiation_model::RadiationModel = SURFACE_ABSORPTION,
)
	# Get the cell adjacent to this boundary
	# Convention: boundary_id=1 is top (z=L), boundary_id=2 is bottom (z=0)
	if boundary_id == 1  # Top surface (exposed, z = L)
		cell_id = mesh.n_active_cells  # Last cell (closest to top)
		face_id = mesh.cells[cell_id].face_ids[2]  # Right face of cell
	else  # boundary_id == 2, Bottom (substrate, z = 0)
		cell_id = 1  # First cell (closest to bottom)
		face_id = mesh.cells[cell_id].face_ids[1]  # Left face of cell
	end

	# Get cell temperature and effective conductivity
	T_cell = mesh.cell_states[cell_id].T
	k_eff = cache.k[cell_id]

	# Effective emissivity for radiation BCs
	ε_eff = cache.ε[cell_id]

	# Distance from boundary to cell center
	# For uniform mesh: Δz_half = cell_volume / face_area / 2
	# This gives half the cell thickness
	Δz_half = mesh.cell_volumes[cell_id] / mesh.face_areas[face_id] / 2

	# Solve for boundary temperature in-place (no allocation)
	state = mesh.boundary_states[boundary_id]
	solve_boundary_temperature!(state, T_cell, k_eff, Δz_half, bc, t, ε_eff;
		radiation_model = radiation_model)

	return state
end

# 2D version (placeholder - same logic but different cell/face indexing)
function update_boundary_state!(
	mesh::Mesh2D,
	boundary_id::Int,
	bc::AbstractThermalBC,
	cache::PropertyCache,
	t::Float64,
	radiation_model::RadiationModel = SURFACE_ABSORPTION,
)
	# TODO: Implement 2D boundary state update
	# For now, just return the existing state
	return mesh.boundary_states[boundary_id]
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
	check_boundary_energy_balance(mesh)

Verify that boundary energy balance is satisfied.

# Returns
Dict{Int, Float64}: Boundary ID to relative imbalance (should be ~machine precision).
"""
function check_boundary_energy_balance(mesh::AbstractMesh)
	results = Dict{Int, Float64}()

	for (boundary_id, state) in mesh.boundary_states
		if abs(state.q_applied) > 1e-10
			imbalance = abs(state.q_applied - state.q_conducted) / abs(state.q_applied)
		else
			imbalance = abs(state.q_applied - state.q_conducted)
		end
		results[boundary_id] = imbalance
	end

	return results
end

"""
	get_surface_temperature(mesh)

Get the top surface temperature from the solved boundary state.
"""
function get_surface_temperature(mesh::AbstractMesh)
	return mesh.boundary_states[1].T_surface
end

"""
	get_back_face_temperature(mesh)

Get the back face (bottom) temperature from the solved boundary state.
"""
function get_back_face_temperature(mesh::AbstractMesh)
	return mesh.boundary_states[2].T_surface
end
