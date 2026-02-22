"""
Radiation flux divergence computation for P1 radiation transport.

This module provides the discretization of radiation flux divergence
for integration with the energy equation, and the quasi-steady solver
for the G field.

The P1 energy equation coupling:
	ρcₚ ∂T/∂t = ∇·(k∇T) + ∇·(D_rad∇G) + Q_rxn + ...
						  ^^^^^^^^^^^^^^
						  This module computes this term

Reference: P1_Volumetric_Reradiation_Implementation.md Section 3.4
"""

using LinearAlgebra: norm
using SparseArrays: spzeros
using ..Core: RadiationBCSet, get_radiation_bc

# Type-stable accessor for radiation BCs
@inline _radiation_bc_at(rbc::RadiationBCSet, ::Val{1}) = get_radiation_bc(rbc, Val(1))
@inline _radiation_bc_at(rbc::RadiationBCSet, ::Val{2}) = get_radiation_bc(rbc, Val(2))
@inline _radiation_bc_at(rbc::RadiationBCSet, id::Int) = get_radiation_bc(rbc, id)

# ============================================================================
# Allocation-Free Tridiagonal Solver (Thomas Algorithm / TDMA)
# ============================================================================

"""
	thomas_solve!(dl, d, du, x, d_work)

Solve tridiagonal system Ax = b in-place using Thomas algorithm (TDMA).

The system has the form:
	[d₁  du₁  0   ...  0  ] [x₁]   [b₁]
	[dl₁ d₂  du₂  ...  0  ] [x₂]   [b₂]
	[0   dl₂  d₃  ...  0  ] [x₃] = [b₃]
	[...          ...     ] [..]   [..]
	[0   ...  dl_{n-1} dₙ ] [xₙ]   [bₙ]

# Arguments
- `dl::Vector{Float64}`: Lower diagonal (length n-1), unchanged on exit
- `d::Vector{Float64}`: Main diagonal (length n), unchanged on exit
- `du::Vector{Float64}`: Upper diagonal (length n-1), unchanged on exit
- `x::Vector{Float64}`: On entry: RHS vector b. On exit: solution x
- `d_work::Vector{Float64}`: Work array (length n) for modified diagonal

# Notes
- This is allocation-free when d_work is pre-allocated
- Uses partial pivoting implicitly (Thomas algorithm is stable for 
  diagonally dominant systems like the P1 radiation equation)
- The original diagonals dl, d, du are preserved
"""
@inline function thomas_solve!(dl::Vector{Float64}, d::Vector{Float64},
	du::Vector{Float64}, x::Vector{Float64},
	d_work::Vector{Float64})
	n = length(d)

	# Copy main diagonal to work array (preserve original for next iteration)
	@inbounds d_work[1] = d[1]

	# Forward elimination
	@inbounds for i in 2:n
		w = dl[i-1] / d_work[i-1]
		d_work[i] = d[i] - w * du[i-1]
		x[i] = x[i] - w * x[i-1]
	end

	# Back substitution
	@inbounds x[n] = x[n] / d_work[n]
	@inbounds for i in (n-1):-1:1
		x[i] = (x[i] - du[i] * x[i+1]) / d_work[i]
	end

	return nothing
end

# ============================================================================
# Radiation Flux Divergence for Energy Equation
# ============================================================================

"""
	apply_radiation_divergence!(dT::Vector{Float64}, mesh, rad_cache, prop_cache)

Add radiation flux divergence contribution to temperature RHS.

	dT[i] += (1/ρcₚ[i]) × (1/V[i]) × (q_rad[R] - q_rad[L]) × A

This implements the ∇·(D_rad∇G) term in the energy equation.

# Arguments
- `dT`: Temperature RHS vector [K/s] (modified in-place)
- `mesh`: Mesh with geometry
- `rad_cache`: Radiation cache with computed q_rad fluxes
- `prop_cache`: Property cache with ρcₚ values

# Notes
The radiation flux q_rad must be computed first using `compute_radiation_fluxes!`.
Sign convention: positive flux = into domain at boundaries, in +z direction for interior.

# Example
```julia
# First compute radiation fluxes
update_radiation_properties!(rad_cache, mesh, prop_cache)
compute_radiation_fluxes!(rad_cache, mesh, prop_cache, radiation_bcs, t)

# Then add to temperature RHS
apply_radiation_divergence!(dT, mesh, rad_cache, prop_cache)
```
"""
function apply_radiation_divergence!(dT::Vector{Float64},
	mesh,
	rad_cache,
	prop_cache)
	n_cells = length(mesh.cells)
	A = mesh.cross_section_area

	@inbounds for i in 1:n_cells
		# Face indices (1D: cell i has faces i and i+1)
		face_L = i
		face_R = i + 1

		# Cell volume and ρcₚ
		V = mesh.cell_volumes[i]
		ρcp = prop_cache.ρcp[i]

		# Radiation flux divergence: (q_out - q_in) / V × A
		# In 1D: flux at right face - flux at left face
		div_q_rad = (rad_cache.q_rad[face_R] - rad_cache.q_rad[face_L]) * A / V

		# Add to temperature RHS: dT/dt += ∇·(D_rad∇G) / ρcₚ
		# Note: We use +div_q_rad because q_rad is defined as heat flux INTO cell
		# from left face (positive in +z) and OUT of cell at right face
		# The divergence ∇·q = dq/dx represents net flux out, so we use -div
		# for conservation. But since q_rad is diffusive flux -D∇G, and
		# we want ∇·(D∇G), we get -∇·q_rad = -div_q_rad
		# 
		# Actually, let's be consistent with the conduction formulation:
		# q_cond = -k dT/dx (heat flows down gradient)
		# div_q = (q_R - q_L) / V × A
		# dT/dt = -div_q / ρcₚ  (minus because we want -∇·q)
		#
		# For radiation:
		# q_rad = -D_rad dG/dx (also flows down gradient)
		# div_q_rad = (q_rad_R - q_rad_L) / V × A
		# The energy equation has ∇·(D∇G) = -∇·q_rad = -div_q_rad
		#
		# So we ADD -div_q_rad to dT:
		dT[i] += -div_q_rad / ρcp
	end

	return nothing
end

# ============================================================================
# Transient Radiation Intensity RHS (for P1_TRANSIENT mode)
# ============================================================================

"""
	compute_radiation_intensity_rhs!(dG::Vector{Float64}, mesh, rad_cache, prop_cache)

Compute RHS for transient radiation intensity equation (dG/dt).

For the P1 transient formulation:
	dG/dt = (c/4σ) × ∇·(D_rad∇G) + α×G - 4σεαT⁴

where c is the speed of light (≈ 3×10⁸ m/s).

# Arguments
- `dG`: Output vector for dG/dt [W/(m²·s)]
- `mesh`: Mesh with temperature state
- `rad_cache`: Radiation cache with q_rad, Q_emission, Q_absorption
- `prop_cache`: Property cache with α values

# Notes
- This is only used for P1_TRANSIENT mode
- For quasi-steady (recommended), dG/dt = 0 is assumed and G is solved algebraically
- The radiation time scale τ_rad = 1/(c×α) ≈ 10⁻¹² s is extremely fast,
  so transient effects are rarely important for pyrolysis

# Example
```julia
if radiation_model == P1_TRANSIENT
	compute_radiation_intensity_rhs!(dG, mesh, rad_cache, prop_cache)
end
```
"""
function compute_radiation_intensity_rhs!(dG::Vector{Float64},
	mesh,
	rad_cache,
	prop_cache)
	n_cells = length(mesh.cells)
	A = mesh.cross_section_area
	c_light = 3.0e8  # Speed of light [m/s]
	σ = STEFAN_BOLTZMANN

	@inbounds for i in 1:n_cells
		# Face indices
		face_L = i
		face_R = i + 1

		# Cell geometry
		V = mesh.cell_volumes[i]

		# Radiation flux divergence
		div_q_rad = (rad_cache.q_rad[face_R] - rad_cache.q_rad[face_L]) * A / V

		# Get properties and source terms
		α = prop_cache.α[i]
		G = rad_cache.G[i]

		# Source terms were already computed
		Q_em = rad_cache.Q_emission[i]
		Q_abs = rad_cache.Q_absorption[i]

		# G equation: dG/dt = (c/4σ) × (-∇·q_rad) + α×G - 4σεαT⁴
		#                   = (c/4σ) × (-div_q_rad) + Q_abs - Q_em
		dG[i] = (c_light / (4.0 * σ)) * (-div_q_rad) + Q_abs - Q_em
	end

	return nothing
end

# ============================================================================
# Quasi-Steady Radiation Solver
# ============================================================================

"""
	solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, t;
								  tol=1e-10, maxiter=20, damping=1.0)

Solve the quasi-steady P1 radiation equation for G.

Since radiation equilibrates much faster than conduction/reaction time scales
(τ_rad ≈ 10⁻¹² s vs τ_cond ≈ 10² s), we can assume dG/dt = 0:

	∇·(D_rad∇G) = αG - 4σεαT⁴

This is a nonlinear algebraic equation (due to T⁴) solved using Newton-Raphson.

# Arguments
- `rad_cache`: Radiation cache (G is updated in-place)
- `mesh`: Mesh with temperature state
- `prop_cache`: Property cache (must be updated before calling)
- `radiation_bcs`: Dict of radiation boundary conditions
- `t`: Current time [s] (supports dual numbers for AD)
- `tol`: Convergence tolerance (default 1e-10)
- `maxiter`: Maximum iterations (default 20)
- `damping`: Damping factor for updates (default 1.0)

# Returns
- `converged::Bool`: True if converged within tolerance
- `iterations::Int`: Number of iterations taken

# Algorithm
Newton-Raphson iteration:
1. Assemble residual R[i] = ∇·(D_rad∇G) - αG + 4σεαT⁴
2. Assemble Jacobian J[i,j] = ∂R[i]/∂G[j]
3. Solve J × ΔG = -R
4. Update G = G + damping × ΔG
5. Check convergence: ||ΔG|| < tol

For 1D, the Jacobian is tridiagonal and can be solved efficiently.

# Notes
- Type-generic to support automatic differentiation

# Example
```julia
# First update radiation properties
update_radiation_properties!(rad_cache, mesh, prop_cache)

# Then solve for quasi-steady G
converged, iters = solve_quasi_steady_radiation!(
	rad_cache, mesh, prop_cache, radiation_bcs, t)

if !converged
	@warn "Quasi-steady radiation did not converge"
end
```
"""
function solve_quasi_steady_radiation!(rad_cache,
	mesh,
	prop_cache,
	radiation_bcs::RadiationBCSet,
	t;
	tol::Float64 = 1e-6,
	maxiter::Int = 20,
	damping::Float64 = 1.0,
	initialize::Bool = true)

	# Fast zero check without allocating a temporary Bool array
	@inline function _all_zero(vec)
		@inbounds for v in vec
			iszero(v) || return false
		end
		return true
	end
	n = length(mesh.cells)
	A = mesh.cross_section_area

	# Initialize G from temperature if requested and G appears uninitialized
	if initialize && _all_zero(rad_cache.G)
		initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
	end

	# Use pre-allocated work arrays from RadiationCache (avoid allocations in hot loop)
	R = rad_cache.newton_R        # Residual
	ΔG = rad_cache.newton_ΔG      # Update

	# Tridiagonal Jacobian storage (pre-allocated)
	dl = rad_cache.newton_dl   # Lower diagonal
	d = rad_cache.newton_d     # Main diagonal
	du = rad_cache.newton_du   # Upper diagonal

	# Compute reference scale for relative tolerance
	# Use max of current G or equilibrium estimate to avoid division issues
	G_scale = max(norm(rad_cache.G, Inf), 1.0)

	for iter in 1:maxiter
		# Update radiation fluxes with current G
		Physics.compute_radiation_fluxes!(rad_cache, mesh, prop_cache, radiation_bcs, t)

		# Compute source terms with current G
		Physics.compute_radiation_source_terms!(
			rad_cache.Q_emission,
			rad_cache.Q_absorption,
			mesh,
			rad_cache,
			prop_cache,
		)

		# Assemble residual: R[i] = ∇·(D_rad∇G) - αG + 4σεαT⁴
		# Note: The equation is  ∇·(D_rad∇G) = αG - 4σεαT⁴
		# So residual R = (αG - 4σεαT⁴) - ∇·(D_rad∇G)
		#               = Q_absorption - Q_emission - ∇·(D_rad∇G)
		@inbounds for i in 1:n
			face_L = i
			face_R = i + 1
			V = mesh.cell_volumes[i]

			# Divergence of radiation flux
			div_q_rad = (rad_cache.q_rad[face_R] - rad_cache.q_rad[face_L]) * A / V

			# Since q_rad = -D_rad dG/dx, div_q_rad = -∇·(D_rad∇G)
			# So ∇·(D_rad∇G) = -div_q_rad

			# Residual: R = αG - 4σεαT⁴ - ∇·(D_rad∇G)
			#            = Q_abs - Q_em - (-div_q_rad)
			#            = Q_abs - Q_em + div_q_rad
			R[i] = rad_cache.Q_absorption[i] - rad_cache.Q_emission[i] + div_q_rad
		end

		# Check convergence using relative tolerance
		# This accounts for the scale of G values which can be very large (10^4 to 10^6 W/m²)
		residual_norm = norm(R, Inf)
		if residual_norm < tol * G_scale
			return (converged = true, iterations = iter)
		end

		# Assemble Jacobian (tridiagonal)
		# J[i,j] = ∂R[i]/∂G[j]
		#
		# R[i] = α[i]×G[i] - Q_em[i] + div_q_rad[i]
		#      = α[i]×G[i] - Q_em[i] + (q_rad[R] - q_rad[L])/V × A
		#
		# q_rad[f] = -D_rad[f] × (G[R] - G[L]) / Δx[f]  for interior face
		#
		# ∂R[i]/∂G[i] = α[i] + (A/V) × [D_rad[R]/Δx[R] + D_rad[L]/Δx[L]]
		# ∂R[i]/∂G[i-1] = -(A/V) × D_rad[L]/Δx[L]
		# ∂R[i]/∂G[i+1] = -(A/V) × D_rad[R]/Δx[R]

		fill!(dl, 0.0)
		fill!(d, 0.0)
		fill!(du, 0.0)

		@inbounds for i in 1:n
			# Use volumetric absorption coefficient from rad_cache (already converted from mass-based)
			α = rad_cache.α[i]
			V = mesh.cell_volumes[i]
			face_L = i
			face_R = i + 1

			# Compute distance to neighbors for flux Jacobian
			if i > 1
				# Left interior face
				D_rad_L = rad_cache.D_rad_face[face_L]
				Δx_L = abs(mesh.cell_centers[i] - mesh.cell_centers[i-1])
				coeff_L = (A / V) * D_rad_L / Δx_L
			else
				coeff_L = 0.0
			end

			if i < n
				# Right interior face
				D_rad_R = rad_cache.D_rad_face[face_R]
				Δx_R = abs(mesh.cell_centers[i+1] - mesh.cell_centers[i])
				coeff_R = (A / V) * D_rad_R / Δx_R
			else
				coeff_R = 0.0
			end

			# Handle boundary faces with Marshak BC
			# For Marshak: flux_into = I_in - G/2 (positive = into domain)
			# ∂flux_into/∂G = -1/2
			# 
			# With our sign convention:
			#   Left boundary (left_id=0): q_rad = +flux_into (into domain = +z)
			#   Right boundary (right_id=0): q_rad = -flux_into (into domain = -z)
			#
			# For left boundary (cell 1):
			#   div_q_rad = (q_rad[R] - q_rad[L]) × A / V
			#   q_rad[L] = flux_into, so ∂q_rad[L]/∂G = -1/2
			#   ∂div_q_rad/∂G[1] from q_rad[L]: -(A/V) × ∂q_rad[L]/∂G[1] = -(A/V) × (-1/2) = +(A/V) × 0.5
			#
			# For right boundary (cell n):
			#   q_rad[R] = -flux_into, so ∂q_rad[R]/∂G = -(-1/2) = +1/2
			#   ∂div_q_rad/∂G[n] from q_rad[R]: +(A/V) × ∂q_rad[R]/∂G[n] = +(A/V) × (+1/2) = +(A/V) × 0.5

			# Left boundary (i == 1)
			if i == 1
				bc_left = _radiation_bc_at(radiation_bcs, 2)  # boundary_id=2 is bottom
				if bc_left !== nothing && bc_left isa MarshakRadiationBC
					# Left boundary: contributes +(A/V) × 0.5 to diagonal
					coeff_L = (A / V) * 0.5
				end
			end

			# Right boundary (i == n)
			if i == n
				bc_right = _radiation_bc_at(radiation_bcs, 1)  # boundary_id=1 is top
				if bc_right !== nothing && bc_right isa MarshakRadiationBC
					# Right boundary: contributes +(A/V) × 0.5 to diagonal  
					coeff_R = (A / V) * 0.5
				end
			end

			# Fill Jacobian
			# Main diagonal
			d[i] = α + coeff_L + coeff_R

			# Off-diagonals
			if i > 1
				dl[i-1] = -coeff_L  # Lower diagonal (∂R[i]/∂G[i-1])
			end
			if i < n
				du[i] = -coeff_R    # Upper diagonal (∂R[i]/∂G[i+1])
			end
		end

		# Solve tridiagonal system: J × ΔG = -R
		# Using allocation-free Thomas algorithm (TDMA)
		# Copy -R into ΔG (RHS for the solve)
		@inbounds for i in 1:n
			ΔG[i] = -R[i]
		end

		# Solve using Thomas algorithm with pre-allocated work array
		thomas_solve!(dl, d, du, ΔG, rad_cache.newton_d_work)

		# Check update norm for early convergence (before applying update)
		# This catches cases where the solution is converged but residual scaling differs
		update_norm = norm(ΔG, Inf)
		relative_update = update_norm / max(norm(rad_cache.G, Inf), 1.0)

		# Update G with damping
		@inbounds for i in 1:n
			rad_cache.G[i] += damping * ΔG[i]
			# Ensure G stays non-negative
			rad_cache.G[i] = max(rad_cache.G[i], 0.0)
		end

		# Converged if relative update is very small
		if relative_update < tol
			return (converged = true, iterations = iter)
		end
	end

	# Only warn if residual is actually large relative to G scale
	final_residual = norm(R, Inf)
	if final_residual > tol * G_scale * 100  # Only warn for significant non-convergence
		@warn "Quasi-steady radiation did not converge in $maxiter iterations (residual: $final_residual, G_scale: $G_scale)"
	end
	return (converged = false, iterations = maxiter)
end

# ============================================================================
# Initialization Helper
# ============================================================================

"""
	initialize_g_from_temperature!(rad_cache, mesh, prop_cache)

Initialize G field from local equilibrium with temperature.

For local thermodynamic equilibrium in P1: G = 4σT⁴ (blackbody)

The standard P1 approximation uses blackbody emission (no emissivity factor)
in the volumetric source term. At equilibrium: αG = 4σαT⁴ → G = 4σT⁴

# Arguments
- `rad_cache`: Radiation cache (G is updated in-place)
- `mesh`: Mesh with temperature state
- `prop_cache`: Property cache (not used, kept for API compatibility)
"""
function initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
	n_cells = length(mesh.cells)
	σ = STEFAN_BOLTZMANN

	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T

		# Local equilibrium: absorption = emission
		# αG = 4σαT⁴  →  G = 4σT⁴ (blackbody)
		rad_cache.G[i] = 4.0 * σ * T^4
	end

	return nothing
end
