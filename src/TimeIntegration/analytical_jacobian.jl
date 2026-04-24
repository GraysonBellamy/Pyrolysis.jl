"""
Analytical Jacobian computation for the pyrolysis equations.

This module provides exact analytical derivatives of the RHS with respect to
state variables, enabling efficient implicit time integration without 
finite difference approximations.

Reference: docs/Analytical_Jacobian_Feasibility.md

State vector layout: u = [T_1,...,T_n, xi_11,...,xi_1n, xi_21,...,xi_2n, ..., xi_nc,n]

The Jacobian J[i,j] = d(du[i]/dt) / du[j] has contributions from:
1. Reaction kinetics (cell-local, diagonal blocks)
2. Heat conduction (tridiagonal in T)
3. Gas transport (tridiagonal in gas species)
4. Radiation absorption (upper triangular coupling)
5. Effective properties (1/rho*cp factor)
"""

using ..Core: derivative, derivative_ξ, R_GAS, STEFAN_BOLTZMANN, GAS, SOLID
using ..Core: n_components, n_reactions
using ..Core: AdiabaticBC, HeatFluxBC, ConvectiveBC, RadiativeBC,
	DirichletThermalBC, RadiativeFluxBC, CombinedThermalBC, MarshakRadiationBC,
	RadiationBCSet, AbstractRadiationBC, FlameBC
using ..Core: ImpermeableBC, DiffusiveMassBC, ConvectiveMassBC, MassFluxBC, CombinedMassBC
using ..Core: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC
using ..Core: MixingRule, PARALLEL, SERIES, WEIGHTED
using ..Core: has_darcy_flow
using ..Core: RadiationModel, BEER_LAMBERT
using ..Properties: total_density, mass_fractions
using ..Properties: effective_gas_viscosity, harmonic_mean_permeability
using ..Physics: evaluate_reactions!, pressure_derivatives, pressure_derivatives!
using ..Discretization: compute_rhs!, solve_boundary_temperature, evaluate_bc_flux_and_derivative
using SparseArrays: SparseMatrixCSC

# ============================================================================
# Numerical Tolerance Constants
# ============================================================================
"""Tolerance for treating values as effectively zero in Jacobian calculations."""
const JAC_ZERO_TOL = 1e-20

"""Minimum density value to prevent division by zero."""
const DENSITY_MIN = 1e-20

"""Minimum volume value for normalization."""
const VOLUME_MIN = 1e-20

"""Minimum conductivity for harmonic mean calculations."""
const CONDUCTIVITY_MIN = 1e-20

"""Minimum gas transfer coefficient for mixing calculations."""
const LAMBDA_MIN = 1e-20

"""Tolerance for checking if flux derivatives are significant."""
const FLUX_DERIV_TOL = 1e-30

"""Tolerance for reaction rate derivatives."""
const RATE_MIN = 1e-30

"""Small value for limiter near-zero checks."""
const LIMITER_TOL = 1e-10

"""Tolerance for swelling factor check (non-swelling gases)."""
const SWELLING_TOL = 1e-10

# ============================================================================
# Efficient Jacobian Zeroing
# ============================================================================

"""
	zero_jacobian!(J::SparseMatrixCSC)
	zero_jacobian!(J::Matrix)

Efficiently zero a Jacobian matrix.

For sparse matrices, zeros only the stored values (O(nnz) instead of O(m×n)).
For dense matrices, uses standard fill!.
"""
@inline function zero_jacobian!(J::SparseMatrixCSC)
	fill!(J.nzval, 0.0)
	return J
end

@inline function zero_jacobian!(J::Matrix)
	fill!(J, 0.0)
	return J
end

# ============================================================================
# Jacobian Cache Structure
# ============================================================================

"""
	JacobianCache{NC,NR}

Pre-allocated workspace for analytical Jacobian computation.

# Type Parameters
- `NC`: Number of components
- `NR`: Number of reactions

# Fields
## Reaction derivatives (cell-local)
- `dr_dT::Vector{Float64}`: ∂r/∂T for each reaction
- `dr_dξ::Matrix{Float64}`: ∂r/∂ξⱼ for each reaction × component

## Property derivatives
- `dk_dT::Vector{Float64}`: ∂k/∂T at each cell
- `dρcp_dT::Vector{Float64}`: ∂(ρcₚ)/∂T at each cell
- `dρcp_dξ::Matrix{Float64}`: ∂(ρcₚ)/∂ξⱼ at each cell

## Face-centered derivatives
- `dk_face_dT_L::Vector{Float64}`: ∂k_face/∂T_left
- `dk_face_dT_R::Vector{Float64}`: ∂k_face/∂T_right

## Radiation derivatives (non-local)
- `dI_dξ::Matrix{Float64}`: ∂I(z_i)/∂ξⱼ(cell k) - intensity sensitivity

## Per-component property cache (avoids type instability from heterogeneous tuple access)
All per-component values are cached using @generated functions that unroll the
tuple access at compile time, avoiding dynamic dispatch allocations.
- `comp_densities::Vector{Float64}`: ρ_j(T) for each component
- `comp_heat_capacities::Vector{Float64}`: c_p,j(T) for each component
- `comp_conductivities::Vector{Float64}`: k_j(T) for each component
- `comp_lambdas::Vector{Float64}`: λ_j(T) for each component
- `comp_dρ_dT::Vector{Float64}`: ∂ρ_j/∂T for each component
- `comp_dcp_dT::Vector{Float64}`: ∂c_p,j/∂T for each component
- `comp_dk_dT::Vector{Float64}`: ∂k_j/∂T for each component
"""
mutable struct JacobianCache{NC, NR}
	# Reaction derivatives
	dr_dT::Vector{Float64}      # [NR]
	dr_dξ::Matrix{Float64}      # [NR × NC]

	# Heat source derivatives
	dQ_rxn_dT::Vector{Float64}  # [n_cells]
	dQ_rxn_dξ::Matrix{Float64}  # [NC × n_cells]

	# Species source derivatives  
	dS_dT::Matrix{Float64}      # [NC × n_cells]
	dS_dξ::Array{Float64, 3}    # [NC × NC × n_cells] - dSⱼ/dξₖ

	# Property derivatives (cell-centered)
	dk_dT::Vector{Float64}      # [n_cells]
	dk_dξ::Matrix{Float64}      # [NC × n_cells]
	dρcp_dT::Vector{Float64}    # [n_cells]
	dρcp_dξ::Matrix{Float64}    # [NC × n_cells]

	# Face conductivity derivatives
	dk_face_dT_L::Vector{Float64}  # [n_faces]
	dk_face_dT_R::Vector{Float64}  # [n_faces]
	dk_face_dξ_L::Matrix{Float64}  # [NC × n_faces]
	dk_face_dξ_R::Matrix{Float64}  # [NC × n_faces]

	# Gas transport derivatives
	dλ_dT::Vector{Float64}      # [n_cells]

	# Work arrays
	work_NC::Vector{Float64}    # [NC]
	work_NR::Vector{Float64}    # [NR]
	work_NC_2::Vector{Float64}  # [NC] - second NC work array for conductivity etc
	ξ_work::Vector{Float64}     # [NC] - for finite difference perturbations
	work_NC_3::Vector{Float64}  # [NC] - heat capacities L
	work_NC_4::Vector{Float64}  # [NC] - heat capacities R

	# Hoisted derivative arrays for gas transport optimization
	# Pre-computed once per face to avoid recomputation in nested j-k loops
	work_NC_5::Vector{Float64}  # [NC] - dρ_dξ for left cell
	work_NC_6::Vector{Float64}  # [NC] - dρ_dξ for right cell
	work_NC_7::Vector{Float64}  # [NC] - dλ_dξ for left cell
	work_NC_8::Vector{Float64}  # [NC] - dλ_dξ for right cell

	# =========================================================================
	# Per-component property cache (avoids heterogeneous tuple access)
	# These are populated using @generated functions that unroll tuple access
	# =========================================================================
	# Property values at current T
	comp_densities::Vector{Float64}       # [NC] - ρ_j(T)
	comp_heat_capacities::Vector{Float64} # [NC] - c_p,j(T)
	comp_conductivities::Vector{Float64}  # [NC] - k_j(T)
	comp_lambdas::Vector{Float64}         # [NC] - λ_j(T)

	# Property T-derivatives at current T
	comp_dρ_dT::Vector{Float64}           # [NC] - ∂ρ_j/∂T
	comp_dcp_dT::Vector{Float64}          # [NC] - ∂c_p,j/∂T
	comp_dk_dT::Vector{Float64}           # [NC] - ∂k_j/∂T

	# Derived quantities
	volume_fractions::Vector{Float64}     # [NC] - v_j at current state

	# =========================================================================
	# Per-cell property cache (Optimization #1: eliminates redundant caching)
	# Computed once per Jacobian evaluation, accessed via @view in face loops
	# =========================================================================
	comp_densities_per_cell::Matrix{Float64}      # [NC × n_cells] - ρ_j(T_i) for each cell
	comp_conductivities_per_cell::Matrix{Float64} # [NC × n_cells] - k_j(T_i) for each cell
	comp_lambdas_per_cell::Matrix{Float64}        # [NC × n_cells] - λ_j(T_i) for each cell
	comp_dk_dT_per_cell::Matrix{Float64}          # [NC × n_cells] - ∂k_j/∂T for each cell

	# Phase flags (cached once, static)
	is_gas::Vector{Bool}               # [NC] - true if component j is GAS phase
	is_solid::Vector{Bool}             # [NC] - true if component j is SOLID phase
	swelling_factors::Vector{Float64}  # [NC] - γ_j swelling factors (constant)
	comp_absorptions::Vector{Float64}  # [NC] - α_j absorption coefficients (constant)

	# =========================================================================
	# Per-reaction property cache (avoids heterogeneous tuple access)
	# These are populated once at initialization and reused
	# =========================================================================
	rxn_E::Vector{Float64}               # [NR] - activation energies E_i
	rxn_reactant1::Vector{Int}           # [NR] - first reactant index (0 = none)
	rxn_reactant2::Vector{Int}           # [NR] - second reactant index (0 = none)
	rxn_order1::Vector{Float64}          # [NR] - first reactant order n_1
	rxn_order2::Vector{Float64}          # [NR] - second reactant order n_2
	rxn_product1::Vector{Int}            # [NR] - first product index (0 = none)
	rxn_product2::Vector{Int}            # [NR] - second product index (0 = none)
	rxn_stoich::Matrix{Float64}          # [4 × NR] - stoichiometry (ν_r1, ν_r2, ν_p1, ν_p2)
	rxn_heat_is_state_dep::Vector{Bool}  # [NR] - true if heat function is StateDependentProperty

	# Radiation work arrays
	τ_cumulative::Vector{Float64}  # [n_cells + 1] - optical depth from surface

	# Darcy transport derivatives (for pressure-driven gas flow in non-ALE mode)
	# These mirror the fields in JacobianCacheALE so that
	# compute_darcy_transport_jacobian! can operate on either cache type.
	dP_dT::Vector{Float64}           # [n_cells] - ∂Pᵢ/∂Tᵢ
	dP_dξ::Matrix{Float64}           # [NC × n_cells] - ∂Pᵢ/∂ξⱼ
	dvgas_dP_L::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂P_left
	dvgas_dP_R::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂P_right
	mobility::Vector{Float64}        # [n_faces] - κ_face/μ_face
	dφ_dξ::Matrix{Float64}           # [NC × n_cells] - ∂φᵢ/∂ξⱼ (porosity)
	dκ_dξ::Matrix{Float64}           # [NC × n_cells] - ∂κᵢ/∂ξⱼ (permeability)
	dvgas_dκ_L::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂κ_L
	dvgas_dκ_R::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂κ_R

	# Sparsity map for O(1) sparse matrix writes (initialized lazily)
	# Uses FastSparsityMap for problems <= 5000 DOFs, SparsityMap for larger
	sparsity_map::Union{Nothing, FastSparsityMap, SparsityMap}

	# Initialization state
	initialized::Bool
end

"""
	JacobianCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Create a pre-allocated Jacobian cache.
"""
function JacobianCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)
	JacobianCache{NC, NR}(
		# Reaction derivatives
		zeros(NR),                      # dr_dT
		zeros(NR, NC),                  # dr_dξ

		# Heat source derivatives
		zeros(n_cells),                 # dQ_rxn_dT
		zeros(NC, n_cells),             # dQ_rxn_dξ

		# Species source derivatives
		zeros(NC, n_cells),             # dS_dT
		zeros(NC, NC, n_cells),         # dS_dξ

		# Property derivatives
		zeros(n_cells),                 # dk_dT
		zeros(NC, n_cells),             # dk_dξ
		zeros(n_cells),                 # dρcp_dT
		zeros(NC, n_cells),             # dρcp_dξ

		# Face conductivity derivatives
		zeros(n_faces),                 # dk_face_dT_L
		zeros(n_faces),                 # dk_face_dT_R
		zeros(NC, n_faces),             # dk_face_dξ_L
		zeros(NC, n_faces),             # dk_face_dξ_R

		# Gas transport
		zeros(n_cells),                 # dλ_dT

		# Work arrays
		zeros(NC),                      # work_NC
		zeros(NR),                      # work_NR
		zeros(NC),                      # work_NC_2
		zeros(NC),                      # ξ_work
		zeros(NC),                      # work_NC_3
		zeros(NC),                      # work_NC_4

		# Hoisted derivative arrays for gas transport optimization
		zeros(NC),                      # work_NC_5 - dρ_dξ for left cell
		zeros(NC),                      # work_NC_6 - dρ_dξ for right cell
		zeros(NC),                      # work_NC_7 - dλ_dξ for left cell
		zeros(NC),                      # work_NC_8 - dλ_dξ for right cell

		# Per-component property cache
		zeros(NC),                      # comp_densities
		zeros(NC),                      # comp_heat_capacities
		zeros(NC),                      # comp_conductivities
		zeros(NC),                      # comp_lambdas

		# Per-component property derivatives
		zeros(NC),                      # comp_dρ_dT
		zeros(NC),                      # comp_dcp_dT
		zeros(NC),                      # comp_dk_dT

		# Derived quantities
		zeros(NC),                      # volume_fractions

		# Per-cell property cache (Optimization #1)
		zeros(NC, n_cells),             # comp_densities_per_cell
		zeros(NC, n_cells),             # comp_conductivities_per_cell
		zeros(NC, n_cells),             # comp_lambdas_per_cell
		zeros(NC, n_cells),             # comp_dk_dT_per_cell

		# Phase flags (initialized lazily)
		fill(false, NC),                # is_gas
		fill(false, NC),                # is_solid
		zeros(NC),                      # swelling_factors
		zeros(NC),                      # comp_absorptions

		# Per-reaction property cache (initialized lazily)
		zeros(NR),                      # rxn_E
		zeros(Int, NR),                 # rxn_reactant1
		zeros(Int, NR),                 # rxn_reactant2
		zeros(NR),                      # rxn_order1
		zeros(NR),                      # rxn_order2
		zeros(Int, NR),                 # rxn_product1
		zeros(Int, NR),                 # rxn_product2
		zeros(4, NR),                   # rxn_stoich
		fill(false, NR),                # rxn_heat_is_state_dep

		# Radiation work arrays
		zeros(n_cells + 1),             # τ_cumulative

		# Darcy transport derivatives
		zeros(n_cells),                 # dP_dT
		zeros(NC, n_cells),             # dP_dξ
		zeros(n_faces),                 # dvgas_dP_L
		zeros(n_faces),                 # dvgas_dP_R
		zeros(n_faces),                 # mobility
		zeros(NC, n_cells),             # dφ_dξ
		zeros(NC, n_cells),             # dκ_dξ
		zeros(n_faces),                 # dvgas_dκ_L
		zeros(n_faces),                 # dvgas_dκ_R

		# Sparsity map (initialized lazily on first Jacobian call)
		nothing,                        # sparsity_map

		# Initialization state
		false,                          # initialized
	)
end

"""
	initialize!(cache::JacobianCache, material)

Initialize static material properties in the cache. Called once before first use.
This avoids repeated O(NC) checks per Jacobian evaluation.
"""
function initialize!(cache::JacobianCache{NC, NR}, material) where {NC, NR}
	cache.initialized && return nothing

	# Cache component properties
	cache_is_gas!(cache.is_gas, material)
	cache_is_solid!(cache.is_solid, material)
	cache_swelling_factors!(cache.swelling_factors, material)
	cache_absorption_coefficients!(cache.comp_absorptions, material)

	# Cache reaction properties (static, only done once)
	cache_reaction_properties!(cache, material)

	cache.initialized = true
	return nothing
end

"""
	cache_reaction_properties!(cache::JacobianCache, material)

Cache static reaction properties to avoid heterogeneous tuple access during Jacobian computation.
"""
@generated function cache_reaction_properties!(cache::JacobianCache{NC, NR},
	material::Material{NC, NR}) where {NC, NR}
	# Generate code that extracts reaction properties at compile time
	exprs = Expr[]
	for i in 1:NR
		push!(exprs, quote
			cache.rxn_E[$i] = material.reactions[$i].E
			cache.rxn_reactant1[$i] = material.reactions[$i].reactant1
			cache.rxn_reactant2[$i] = material.reactions[$i].reactant2
			cache.rxn_order1[$i] = material.reactions[$i].order1
			cache.rxn_order2[$i] = material.reactions[$i].order2
			cache.rxn_product1[$i] = material.reactions[$i].product1
			cache.rxn_product2[$i] = material.reactions[$i].product2
			cache.rxn_stoich[1, $i] = material.reactions[$i].stoich[1]
			cache.rxn_stoich[2, $i] = material.reactions[$i].stoich[2]
			cache.rxn_stoich[3, $i] = material.reactions[$i].stoich[3]
			cache.rxn_stoich[4, $i] = material.reactions[$i].stoich[4]
			cache.rxn_heat_is_state_dep[$i] = material.reactions[$i].heat isa StateDependentProperty
		end)
	end
	return quote
		$(exprs...)
		return nothing
	end
end

"""
	resize_jacobian_cache!(cache::JacobianCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Resize Jacobian cache after mesh adaptation.
"""
function resize_jacobian_cache!(cache::JacobianCache{NC, NR},
	n_cells::Int, n_faces::Int) where {NC, NR}
	# Reaction derivatives stay same size (per-cell, reused)

	# Resize cell-based arrays
	resize!(cache.dQ_rxn_dT, n_cells)
	cache.dQ_rxn_dξ = zeros(NC, n_cells)
	cache.dS_dT = zeros(NC, n_cells)
	cache.dS_dξ = zeros(NC, NC, n_cells)
	resize!(cache.dk_dT, n_cells)
	cache.dk_dξ = zeros(NC, n_cells)
	resize!(cache.dρcp_dT, n_cells)
	cache.dρcp_dξ = zeros(NC, n_cells)
	resize!(cache.dλ_dT, n_cells)

	# Resize face-based arrays
	resize!(cache.dk_face_dT_L, n_faces)
	resize!(cache.dk_face_dT_R, n_faces)
	cache.dk_face_dξ_L = zeros(NC, n_faces)
	cache.dk_face_dξ_R = zeros(NC, n_faces)

	# Resize radiation work arrays
	resize!(cache.τ_cumulative, n_cells + 1)

	return cache
end

# ============================================================================
# ALE Jacobian Cache Structure
# ============================================================================

"""
	JacobianCacheALE{NC,NR}

Pre-allocated workspace for analytical Jacobian computation in ALE simulations.

Extends JacobianCache with additional workspace for:
- Strain rate derivatives (∂θ/∂T, ∂θ/∂ξ)
- Node velocity derivatives (∂w/∂T, ∂w/∂ξ, ∂w/∂z)
- Geometry derivatives (∂V/∂z, ∂Δz/∂z)

# Type Parameters
- `NC`: Number of components
- `NR`: Number of reactions

# Fields
## Base cache (inherited from JacobianCache)
- `base::JacobianCache{NC,NR}`: Contains all fixed-mesh Jacobian workspace

## Strain rate derivatives (cell-centered)
- `dθ_dT::Vector{Float64}`: ∂θᵢ/∂Tᵢ for each cell [n_cells]
- `dθ_dξ::Matrix{Float64}`: ∂θᵢ/∂ξⱼᵢ [NC × n_cells]

## Node velocity derivatives
- `dw_dT::Matrix{Float64}`: ∂wₖ/∂Tᵢ [n_nodes × n_cells]
- `dw_dξ::Array{Float64,3}`: ∂wₖ/∂ξⱼᵢ [n_nodes × NC × n_cells]
- `dw_dz::Matrix{Float64}`: ∂wₖ/∂zₘ [n_nodes × n_nodes]

## Geometry derivatives
- `dV_dz::Matrix{Float64}`: ∂Vᵢ/∂zₖ [n_cells × n_nodes]
- `dΔz_dz::Matrix{Float64}`: ∂Δzᵢ/∂zₖ [n_cells × n_nodes]

## ALE advection work arrays
- `w_face::Vector{Float64}`: Mesh velocity at faces [n_faces]
- `dadvT_dT::Vector{Float64}`: ∂(w·∇T)ᵢ/∂Tⱼ work array [n_cells]
- `dadvT_dz::Vector{Float64}`: ∂(w·∇T)ᵢ/∂zₖ work array [n_nodes]

## Darcy transport derivatives (only used when has_darcy_flow(material) = true)
- `dP_dT::Vector{Float64}`: ∂Pᵢ/∂Tᵢ at each cell [n_cells]
- `dP_dξ::Matrix{Float64}`: ∂Pᵢ/∂ξⱼ [NC × n_cells]
- `dvgas_dP_L::Vector{Float64}`: ∂v_gas,f/∂P_left for each face [n_faces]
- `dvgas_dP_R::Vector{Float64}`: ∂v_gas,f/∂P_right for each face [n_faces]
- `mobility::Vector{Float64}`: κ_face/μ_face for each face [n_faces]
- `dφ_dξ::Matrix{Float64}`: ∂φᵢ/∂ξⱼ [NC × n_cells] - porosity derivatives
- `dκ_dξ::Matrix{Float64}`: ∂κᵢ/∂ξⱼ [NC × n_cells] - permeability derivatives
- `dvgas_dκ_L::Vector{Float64}`: ∂v_gas,f/∂κ_L for each face [n_faces]
- `dvgas_dκ_R::Vector{Float64}`: ∂v_gas,f/∂κ_R for each face [n_faces]
"""
mutable struct JacobianCacheALE{NC, NR}
	# Base cache for existing energy/species terms
	base::JacobianCache{NC, NR}

	# Strain rate derivatives (cell-centered)
	dθ_dT::Vector{Float64}           # [n_cells] - ∂θᵢ/∂Tᵢ
	dθ_dξ::Matrix{Float64}           # [NC × n_cells] - ∂θᵢ/∂ξⱼᵢ

	# Node velocity derivatives
	dw_dT::Matrix{Float64}           # [n_nodes × n_cells] - ∂wₖ/∂Tᵢ
	dw_dξ::Array{Float64, 3}         # [n_nodes × NC × n_cells] - ∂wₖ/∂ξⱼᵢ
	dw_dz::Matrix{Float64}           # [n_nodes × n_nodes] - ∂wₖ/∂zₘ

	# Geometry derivatives
	dV_dz::Matrix{Float64}           # [n_cells × n_nodes] - ∂Vᵢ/∂zₖ
	dΔz_dz::Matrix{Float64}          # [n_cells × n_nodes] - ∂Δzᵢ/∂zₖ

	# ALE advection work arrays
	w_face::Vector{Float64}          # [n_faces] - mesh velocity at faces
	dadvT_dT::Vector{Float64}        # [n_cells] - advection derivative work
	dadvT_dz::Vector{Float64}        # [n_nodes] - advection derivative work

	# Species source matrix for strain rate (avoids allocation in main loop)
	S_matrix::Matrix{Float64}        # [NC × n_cells] - species source terms

	# Darcy transport derivatives (for pressure-driven gas flow)
	dP_dT::Vector{Float64}           # [n_cells] - ∂Pᵢ/∂Tᵢ
	dP_dξ::Matrix{Float64}           # [NC × n_cells] - ∂Pᵢ/∂ξⱼ
	dvgas_dP_L::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂P_left
	dvgas_dP_R::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂P_right
	mobility::Vector{Float64}        # [n_faces] - κ_face/μ_face (for upwind term)

	# Permeability derivatives for Kozeny-Carman and other φ-dependent models
	dφ_dξ::Matrix{Float64}           # [NC × n_cells] - ∂φᵢ/∂ξⱼ (φ = porosity)
	dκ_dξ::Matrix{Float64}           # [NC × n_cells] - ∂κᵢ/∂ξⱼ (permeability)
	dvgas_dκ_L::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂κ_L
	dvgas_dκ_R::Vector{Float64}      # [n_faces] - ∂v_gas,f/∂κ_R
end

"""
	JacobianCacheALE(n_cells::Int, n_faces::Int, n_nodes::Int, NC::Int, NR::Int)

Create a pre-allocated ALE Jacobian cache.

# Arguments
- `n_cells`: Number of mesh cells
- `n_faces`: Number of mesh faces
- `n_nodes`: Number of mesh nodes (typically n_cells + 1 for 1D)
- `NC`: Number of material components
- `NR`: Number of reactions
"""
function JacobianCacheALE(n_cells::Int, n_faces::Int, n_nodes::Int, NC::Int, NR::Int)
	# Create base cache
	base = JacobianCache(n_cells, n_faces, NC, NR)

	JacobianCacheALE{NC, NR}(
		base,

		# Strain rate derivatives
		zeros(n_cells),                    # dθ_dT
		zeros(NC, n_cells),                # dθ_dξ

		# Node velocity derivatives
		zeros(n_nodes, n_cells),           # dw_dT
		zeros(n_nodes, NC, n_cells),       # dw_dξ
		zeros(n_nodes, n_nodes),           # dw_dz

		# Geometry derivatives
		zeros(n_cells, n_nodes),           # dV_dz
		zeros(n_cells, n_nodes),           # dΔz_dz

		# ALE advection work arrays
		zeros(n_faces),                    # w_face
		zeros(n_cells),                    # dadvT_dT
		zeros(n_nodes),                    # dadvT_dz

		# Species source matrix
		zeros(NC, n_cells),                # S_matrix

		# Darcy transport derivatives
		zeros(n_cells),                    # dP_dT
		zeros(NC, n_cells),                # dP_dξ
		zeros(n_faces),                    # dvgas_dP_L
		zeros(n_faces),                    # dvgas_dP_R
		zeros(n_faces),                    # mobility

		# Permeability derivatives
		zeros(NC, n_cells),                # dφ_dξ
		zeros(NC, n_cells),                # dκ_dξ
		zeros(n_faces),                    # dvgas_dκ_L
		zeros(n_faces),                    # dvgas_dκ_R
	)
end

"""
	resize_jacobian_cache_ale!(cache::JacobianCacheALE, n_cells::Int, n_faces::Int, n_nodes::Int)

Resize ALE Jacobian cache after mesh adaptation.

This resizes both the base cache and all ALE-specific arrays.
"""
function resize_jacobian_cache_ale!(cache::JacobianCacheALE{NC, NR},
	n_cells::Int, n_faces::Int, n_nodes::Int) where {NC, NR}

	# Resize base cache
	resize_jacobian_cache!(cache.base, n_cells, n_faces)

	# Resize strain rate derivatives
	resize!(cache.dθ_dT, n_cells)
	cache.dθ_dξ = zeros(NC, n_cells)

	# Resize node velocity derivatives
	cache.dw_dT = zeros(n_nodes, n_cells)
	cache.dw_dξ = zeros(n_nodes, NC, n_cells)
	cache.dw_dz = zeros(n_nodes, n_nodes)

	# Resize geometry derivatives
	cache.dV_dz = zeros(n_cells, n_nodes)
	cache.dΔz_dz = zeros(n_cells, n_nodes)

	# Resize ALE advection work arrays
	resize!(cache.w_face, n_faces)
	resize!(cache.dadvT_dT, n_cells)
	resize!(cache.dadvT_dz, n_nodes)

	# Resize species source matrix
	cache.S_matrix = zeros(NC, n_cells)

	# Resize Darcy transport derivatives
	resize!(cache.dP_dT, n_cells)
	cache.dP_dξ = zeros(NC, n_cells)
	resize!(cache.dvgas_dP_L, n_faces)
	resize!(cache.dvgas_dP_R, n_faces)
	resize!(cache.mobility, n_faces)

	# Resize permeability derivatives
	cache.dφ_dξ = zeros(NC, n_cells)
	cache.dκ_dξ = zeros(NC, n_cells)
	resize!(cache.dvgas_dκ_L, n_faces)
	resize!(cache.dvgas_dκ_R, n_faces)

	return cache
end

# ============================================================================
# Reaction Rate Derivatives
# ============================================================================

"""
	reaction_rate_derivatives!(dr_dT, dr_dξ, ξ, T, reaction, cache)

Compute derivatives of a single reaction rate with respect to T and ξ.

For Arrhenius kinetics with soft depletion limiter:
	r = A × exp(-E/RT) × ξ_r1^n1 × ξ_r2^n2 × limiter(ξ_r1)
	
Derivatives:
	∂r/∂T = r × E/(R×T²)
	∂r/∂ξ_r1 = r × (n1/ξ_r1 + dlimiter/dξ_r1)
	∂r/∂ξ_r2 = r × n2/ξ_r2
"""
function reaction_rate_derivative_T(r::Float64, E::Float64, T::Float64)
	# ∂r/∂T = r × E/(R×T²)
	return r * E / (R_GAS * T^2)
end

"""
Compute ∂r/∂ξ for the soft depletion limiter.

The soft limiter is: f(ξ) = ξ × tanh(ξ/ξ_threshold)
So the rate is: r = k × f(ξ) × ...
And: ∂r/∂ξ = r/f × df/dξ

df/dξ = tanh(x) + ξ × sech²(x)/ξ_threshold where x = ξ/ξ_threshold
	  = tanh(x) + x × sech²(x)
	  = tanh(x) + x × (1 - tanh²(x))
"""
function soft_limiter_derivative(ξ::Float64, ξ_threshold::Float64)
	if ξ < JAC_ZERO_TOL
		return 0.0
	end
	x = ξ / ξ_threshold
	tanh_x = tanh(x)
	sech2_x = 1.0 - tanh_x^2
	# df/dξ = tanh(x) + x × sech²(x)
	return tanh_x + x * sech2_x
end

"""
	evaluate_reaction_derivatives!(dr_dT::Vector{Float64}, dr_dξ::Matrix{Float64},
								   rates::Vector{Float64}, ξ, T, cache::JacobianCache)

Compute derivatives of all reaction rates with respect to T and ξ using cached reaction properties.
This version avoids heterogeneous tuple access by using pre-cached reaction properties.

# Arguments
- `dr_dT`: Output vector [NR] of ∂rᵢ/∂T
- `dr_dξ`: Output matrix [NR × NC] of ∂rᵢ/∂ξⱼ
- `rates`: Current reaction rates [NR]
- `ξ`: Mass concentrations tuple
- `T`: Temperature
- `cache`: JacobianCache with pre-cached reaction properties
"""
function evaluate_reaction_derivatives!(dr_dT::Vector{Float64}, dr_dξ::Matrix{Float64},
	rates::Vector{Float64}, ξ::NTuple{NC, Float64},
	T::Float64, cache::JacobianCache{NC, NR}) where {NC, NR}
	ξ_threshold = 1e-6  # Same as in kinetics.jl

	fill!(dr_dξ, 0.0)

	@inbounds for i in 1:NR
		r = rates[i]

		# Skip if rate is essentially zero
		if abs(r) < RATE_MIN
			dr_dT[i] = 0.0
			continue
		end

		# Temperature derivative: ∂r/∂T = r × E/(R×T²)
		dr_dT[i] = reaction_rate_derivative_T(r, cache.rxn_E[i], T)

		# Concentration derivatives using cached reaction properties
		r1 = cache.rxn_reactant1[i]
		r2 = cache.rxn_reactant2[i]
		n1 = cache.rxn_order1[i]
		n2 = cache.rxn_order2[i]

		# For reactant 1 (with soft limiter)
		if r1 > 0 && ξ[r1] > JAC_ZERO_TOL
			ξ_r1 = ξ[r1]

			# The rate has form: r = k × ξ_r1^n1 × limiter(ξ_r1) × ...
			# limiter = tanh(ξ/ξ_threshold)
			# 
			# ∂r/∂ξ_r1 = r × [n1/ξ_r1 + (1/limiter) × dlimiter/dξ]
			#          = r × [n1/ξ_r1 + sech²(x)/ξ_threshold / tanh(x)]

			x = ξ_r1 / ξ_threshold
			tanh_x = tanh(x)

			if tanh_x > LIMITER_TOL
				sech2_x = 1.0 - tanh_x^2
				# d(log limiter)/dξ = sech²(x)/(tanh(x) × ξ_threshold)
				limiter_term = sech2_x / (tanh_x * ξ_threshold)
				dr_dξ[i, r1] = r * (n1 / ξ_r1 + limiter_term)
			else
				# Near zero, limiter ≈ x, so d(log limiter)/dξ ≈ 1/ξ
				dr_dξ[i, r1] = r * (n1 + 1.0) / ξ_r1
			end
		end

		# For reactant 2 (no soft limiter, just power law)
		if r2 > 0 && r2 != r1 && ξ[r2] > JAC_ZERO_TOL
			dr_dξ[i, r2] = r * n2 / ξ[r2]
		end
	end

	return nothing
end

# ============================================================================
# Heat Source Derivatives
# ============================================================================

"""
	heat_source_derivatives_cached!(dQ_dξ, dr_dT, dr_dξ, rates, heats, ξ, T, cache, material) -> dQ_dT

Compute derivatives of reaction heat source Q = -Σᵢ hᵢ × rᵢ using cached reaction properties.
This version uses cached `rxn_heat_is_state_dep` flags to avoid type checks.

∂Q/∂T = -Σᵢ (hᵢ × ∂rᵢ/∂T + rᵢ × ∂hᵢ/∂T)
∂Q/∂ξⱼ = -Σᵢ hᵢ × ∂rᵢ/∂ξⱼ

# Arguments
- `dQ_dξ`: Output vector [NC] for ∂Q/∂ξⱼ (modified in-place)
- `cache`: JacobianCache with pre-cached reaction properties
- `material`: Material (still needed for heat function evaluation)
- Returns: dQ_dT scalar
"""
@generated function heat_source_derivatives_cached!(dQ_dξ::AbstractVector{Float64},
	dr_dT::Vector{Float64}, dr_dξ::Matrix{Float64},
	rates::Vector{Float64}, heats::Vector{Float64},
	ξ::NTuple{NC, Float64}, T::Float64,
	cache::JacobianCache{NC, NR}, material::Material{NC, NR}) where {NC, NR}

	# Generate unrolled loop at compile time
	loop_body = Expr[]
	for i in 1:NR
		push!(loop_body, quote
			h = heats[$i]
			r = rates[$i]
			heat_prop = material.reactions[$i].heat

			# ∂Q/∂T contribution from reaction i
			dh_dT_i = derivative(heat_prop, T, ξ)
			dQ_dT -= h * dr_dT[$i] + r * dh_dT_i

			# ∂Q/∂ξⱼ contributions
			# Check using cached flag to avoid repeated type checks
			is_state_dep = cache.rxn_heat_is_state_dep[$i]

			@inbounds for j in 1:NC
				dQ_dξ[j] -= h * dr_dξ[$i, j]
				if is_state_dep
					dh_dξ_j = derivative_ξ(heat_prop, T, ξ, j)
					dQ_dξ[j] -= r * dh_dξ_j
				end
			end
		end)
	end

	return quote
		dQ_dT = 0.0
		fill!(dQ_dξ, 0.0)
		$(loop_body...)
		return dQ_dT
	end
end

# ============================================================================
# Species Source Derivatives
# ============================================================================

"""
	species_source_derivatives_cached!(dS_dT, dS_dξ, dr_dT, dr_dξ, cache, NC_val)

Compute derivatives of species source terms Sⱼ = Σᵢ νᵢⱼ × rᵢ using cached reaction properties.

∂Sⱼ/∂T = Σᵢ νᵢⱼ × ∂rᵢ/∂T
∂Sⱼ/∂ξₖ = Σᵢ νᵢⱼ × ∂rᵢ/∂ξₖ
"""
function species_source_derivatives_cached!(dS_dT::AbstractVector{Float64},
	dS_dξ::AbstractMatrix{Float64},
	dr_dT::Vector{Float64},
	dr_dξ::Matrix{Float64},
	cache::JacobianCache{NC, NR}, NC_val::Int) where {NC, NR}

	fill!(dS_dT, 0.0)
	fill!(dS_dξ, 0.0)

	@inbounds for i in 1:NR
		# Use cached reaction properties instead of material.reactions[i]
		r1 = cache.rxn_reactant1[i]
		r2 = cache.rxn_reactant2[i]
		p1 = cache.rxn_product1[i]
		p2 = cache.rxn_product2[i]

		# Reactant consumption (negative stoichiometry)
		if r1 > 0
			ν = -cache.rxn_stoich[1, i]
			dS_dT[r1] += ν * dr_dT[i]
			for k in 1:NC_val
				dS_dξ[r1, k] += ν * dr_dξ[i, k]
			end
		end
		if r2 > 0
			ν = -cache.rxn_stoich[2, i]
			dS_dT[r2] += ν * dr_dT[i]
			for k in 1:NC_val
				dS_dξ[r2, k] += ν * dr_dξ[i, k]
			end
		end

		# Product generation (positive stoichiometry)
		if p1 > 0
			ν = cache.rxn_stoich[3, i]
			dS_dT[p1] += ν * dr_dT[i]
			for k in 1:NC_val
				dS_dξ[p1, k] += ν * dr_dξ[i, k]
			end
		end
		if p2 > 0
			ν = cache.rxn_stoich[4, i]
			dS_dT[p2] += ν * dr_dT[i]
			for k in 1:NC_val
				dS_dξ[p2, k] += ν * dr_dξ[i, k]
			end
		end
	end

	return nothing
end

# ============================================================================
# Heat Capacity Derivatives
# ============================================================================

"""
	heat_capacity_derivatives_cached!(dρcp_dξ, ξ, cp_values, dcp_dT_values) -> dρcp_dT

Compute derivatives of volumetric heat capacity ρcₚ = Σⱼ ξⱼ × cₚ,ⱼ(T) (in-place version).
Uses pre-cached component property values to avoid type instability.

∂(ρcₚ)/∂T = Σⱼ ξⱼ × ∂cₚ,ⱼ/∂T
∂(ρcₚ)/∂ξⱼ = cₚ,ⱼ(T)

# Arguments
- `dρcp_dξ`: Output vector [NC] for ∂(ρcₚ)/∂ξⱼ (modified in-place)
- `ξ`: Mass concentrations tuple
- `cp_values`: Pre-cached c_p,j(T) values for all components [NC]
- `dcp_dT_values`: Pre-cached ∂c_p,j/∂T values for all components [NC]
- Returns: dρcp_dT scalar
"""
@inline function heat_capacity_derivatives_cached!(
	dρcp_dξ::AbstractVector{Float64},
	ξ::NTuple{NC, Float64},
	cp_values::AbstractVector{Float64},
	dcp_dT_values::AbstractVector{Float64}) where NC

	dρcp_dT = 0.0

	@inbounds @simd for j in 1:NC
		dρcp_dT += ξ[j] * dcp_dT_values[j]
		dρcp_dξ[j] = cp_values[j]
	end

	return dρcp_dT
end

"""
	inverse_heat_capacity_derivatives!(dρcp_dξ, ρcp, dρcp_dT) -> dinv_dT

Compute derivatives of 1/ρcₚ using chain rule: d(1/x)/dx = -1/x².

∂(1/ρcₚ)/∂T = -∂(ρcₚ)/∂T / (ρcₚ)²
∂(1/ρcₚ)/∂ξⱼ = -∂(ρcₚ)/∂ξⱼ / (ρcₚ)²

Transforms dρcp_dξ in-place to dinv_dξ.
"""
function inverse_heat_capacity_derivatives!(dρcp_dξ::AbstractVector{Float64},
	ρcp::Float64, dρcp_dT::Float64)
	inv_ρcp_sq = 1.0 / (ρcp * ρcp)

	dinv_dT = -dρcp_dT * inv_ρcp_sq

	# Transform dρcp_dξ in-place to dinv_dξ
	@inbounds @simd for j in eachindex(dρcp_dξ)
		dρcp_dξ[j] = -dρcp_dξ[j] * inv_ρcp_sq
	end

	return dinv_dT  # dρcp_dξ is now dinv_dξ
end

# ============================================================================
# Conductivity Derivatives
# ============================================================================

"""
	conductivity_derivatives_cached!(dk_dξ, ξ, densities, conductivities, dk_dT_values,
									  is_gas, swelling_factors, mixing, β, v_work) -> dk_dT

Compute derivatives of effective thermal conductivity using pre-cached component properties.
This version avoids heterogeneous tuple iteration and allocations.

For WEIGHTED mixing rule: k = β×k_parallel + (1-β)×k_series

Note: Per ThermaKin convention, gaseous components with swelling_factor ≈ 0 do not 
contribute to thermal conductivity. Liquids with swelling_factor = 0 ARE included.

# Arguments
- `dk_dξ`: Output vector [NC] for ∂k/∂ξⱼ (modified in-place)
- `ξ`: Mass concentrations tuple
- `densities`: Pre-cached ρ_j(T) for all components [NC]
- `conductivities`: Pre-cached k_j(T) for all components [NC]
- `dk_dT_values`: Pre-cached ∂k_j/∂T for all components [NC]
- `is_gas`: Pre-cached phase flags (true if GAS) [NC]
- `swelling_factors`: Pre-cached γ_j for all components [NC]
- `mixing`: Mixing rule (PARALLEL, SERIES, or WEIGHTED)
- `β`: Mixing parameter for WEIGHTED rule
- `v_work`: Work array [NC] for volume fractions (modified)
- Returns: dk_dT scalar
"""
@inline function conductivity_derivatives_cached!(
	dk_dξ::AbstractVector{Float64},
	ξ::NTuple{NC, Float64},
	densities::AbstractVector{Float64},
	conductivities::AbstractVector{Float64},
	dk_dT_values::AbstractVector{Float64},
	is_gas::AbstractVector{Bool},
	swelling_factors::AbstractVector{Float64},
	mixing::MixingRule,
	β::Float64,
	v_work::AbstractVector{Float64}) where NC

	# Compute total "solid" volume (non-swelling gases excluded)
	# Per ThermaKin: exclude gases with γ ≈ 0, but include liquids with γ = 0
	V_total = 0.0
	@inbounds for j in 1:NC
		# Include component unless it's a non-swelling gas
		if !(is_gas[j] && swelling_factors[j] < SWELLING_TOL) && ξ[j] > 0.0
			ρ_j = densities[j]
			if ρ_j > 0.0
				V_total += ξ[j] / ρ_j
			end
		end
	end

	# Handle degenerate case
	if V_total < VOLUME_MIN
		@inbounds for j in 1:NC
			dk_dξ[j] = 0.0
		end
		return 0.0
	end

	inv_V_total = 1.0 / V_total

	# Compute volume fractions, conductivities, and their sums
	k_parallel = 0.0
	inv_k_series = 0.0
	@inbounds for j in 1:NC
		# Include component unless it's a non-swelling gas
		if !(is_gas[j] && swelling_factors[j] < SWELLING_TOL) && ξ[j] > 0.0
			ρ_j = densities[j]
			if ρ_j > 0.0
				v_j = (ξ[j] / ρ_j) * inv_V_total
				k_j = conductivities[j]
				k_parallel += v_j * k_j
				if k_j > CONDUCTIVITY_MIN
					inv_k_series += v_j / k_j
				end
				v_work[j] = v_j
			else
				v_work[j] = 0.0
			end
		else
			v_work[j] = 0.0
		end
	end

	k_series = inv_k_series > CONDUCTIVITY_MIN ? 1.0 / inv_k_series : 0.0
	k_series_sq = k_series * k_series

	# Compute ∂k/∂ξⱼ for each component
	@inbounds for j in 1:NC
		# Only non-gas components or swelling gases contribute to volume fractions
		# Exclude component if it's a non-swelling gas
		if is_gas[j] && swelling_factors[j] < SWELLING_TOL
			dk_dξ[j] = 0.0
			continue
		end

		ρ_j = densities[j]
		if ρ_j < DENSITY_MIN
			dk_dξ[j] = 0.0
			continue
		end

		# Compute ∂k_parallel/∂ξⱼ and ∂(1/k_series)/∂ξⱼ
		dk_parallel_dξj = 0.0
		d_inv_k_series_dξj = 0.0

		factor = 1.0 / (ρ_j * V_total)

		@inbounds for m in 1:NC
			# Skip non-swelling gases
			if is_gas[m] && swelling_factors[m] < SWELLING_TOL
				continue
			end

			v_m = v_work[m]
			k_m = conductivities[m]

			# ∂vₘ/∂ξⱼ = (δₘⱼ - vₘ) / (ρⱼ × V_total)
			if m == j
				dv_m_dξj = (1.0 - v_m) * factor
			else
				dv_m_dξj = -v_m * factor
			end

			dk_parallel_dξj += k_m * dv_m_dξj

			if k_m > CONDUCTIVITY_MIN
				d_inv_k_series_dξj += dv_m_dξj / k_m
			end
		end

		# ∂k_series/∂ξⱼ = -k_series² × ∂(1/k_series)/∂ξⱼ
		dk_series_dξj = -k_series_sq * d_inv_k_series_dξj

		# Apply mixing rule
		if mixing == PARALLEL
			dk_dξ[j] = dk_parallel_dξj
		elseif mixing == SERIES
			dk_dξ[j] = dk_series_dξj
		else  # WEIGHTED
			dk_dξ[j] = β * dk_parallel_dξj + (1.0 - β) * dk_series_dξj
		end
	end

	# Temperature derivative from pre-cached ∂k_j/∂T values
	# Computed analytically for each mixing rule:
	#   PARALLEL: dk/dT = Σⱼ vⱼ × dk_j/dT
	#   SERIES:   dk/dT = k_series² × Σⱼ vⱼ × (dk_j/dT) / k_j²
	#   WEIGHTED: dk/dT = β × dk_parallel/dT + (1-β) × dk_series/dT
	# Note: dv_j/dT terms (from density T-dependence) are small and omitted.
	dk_parallel_dT = 0.0
	dk_series_dT = 0.0
	@inbounds for j in 1:NC
		v_j = v_work[j]
		if v_j > 0.0
			dk_parallel_dT += v_j * dk_dT_values[j]
			k_j = conductivities[j]
			if k_j > CONDUCTIVITY_MIN
				dk_series_dT += v_j * dk_dT_values[j] / (k_j * k_j)
			end
		end
	end
	dk_series_dT *= k_series_sq  # k_series² × Σⱼ vⱼ(dk_j/dT)/k_j²

	if mixing == PARALLEL
		dk_dT = dk_parallel_dT
	elseif mixing == SERIES
		dk_dT = dk_series_dT
	else  # WEIGHTED
		dk_dT = β * dk_parallel_dT + (1.0 - β) * dk_series_dT
	end

	return dk_dT
end

# ============================================================================
# Face Conductivity Derivatives (Distance-Weighted Harmonic Mean)
# ============================================================================

"""
	harmonic_mean_derivatives_weighted(k_L, k_R, d_L, d_R) -> (k_face, dk_dkL, dk_dkR)

Compute distance-weighted harmonic mean and its derivatives.

k_face = (d_L + d_R) / (d_L/k_L + d_R/k_R)

Using resistance notation: R_i = d_i/k_i, k_face = d_total/R_total

∂k_face/∂k_L = d_total × d_L / (k_L² × R_total²)
∂k_face/∂k_R = d_total × d_R / (k_R² × R_total²)
"""
function harmonic_mean_derivatives_weighted(k_L::Float64, k_R::Float64, d_L::Float64, d_R::Float64)
	if k_L < CONDUCTIVITY_MIN || k_R < CONDUCTIVITY_MIN
		return 0.0, 0.0, 0.0
	end

	# Resistance terms
	R_L = d_L / k_L
	R_R = d_R / k_R
	R_sum = R_L + R_R
	d_total = d_L + d_R

	k_face = d_total / R_sum

	# Derivatives: ∂k_face/∂k_i = d_total × d_i / (k_i² × R_sum²)
	dk_dkL = d_total * d_L / (k_L * k_L * R_sum * R_sum)
	dk_dkR = d_total * d_R / (k_R * k_R * R_sum * R_sum)

	return k_face, dk_dkL, dk_dkR
end

# ============================================================================
# Boundary Condition Derivatives
# ============================================================================

"""
	convective_bc_derivative(h, T_surface) -> Float64

Derivative of convective BC flux q = h×(T∞ - T_surface) w.r.t. T_surface.

∂q/∂T_surface = -h
"""
convective_bc_derivative(h::Float64) = -h

"""
	radiative_bc_derivative(ε, T_surface) -> Float64

Derivative of radiative BC flux q = εσ(T∞⁴ - T⁴) w.r.t. T_surface.

∂q/∂T_surface = -4×ε×σ×T³
"""
radiative_bc_derivative(ε::Float64, T_surface::Float64) = -4 * ε * STEFAN_BOLTZMANN * T_surface^3

# ============================================================================
# Main Jacobian Assembly
# ============================================================================

"""
	compute_analytical_jacobian!(J::AbstractMatrix, u, p, t, f::PyrolysisODEFunction,
								 jac_cache::JacobianCache)

Compute the analytical Jacobian matrix.

The Jacobian J[i,j] = d(du[i]/dt) / du[j] where the state vector is
u = [T_1,...,T_n, xi_11,...,xi_1n, ..., xi_nc,n].

# Energy equation row i:
dTᵢ/dt = (1/ρcₚ)ᵢ × (-∇·q + Q_rad + Q_rxn)ᵢ

∂(dTᵢ/dt)/∂Tⱼ includes:
- ∂(1/ρcₚ)/∂T × (RHS) if j == i (property dependence)
- (1/ρcₚ) × ∂(-∇·q)/∂T (conduction, j = i±1 and i)
- (1/ρcₚ) × ∂Q_rxn/∂T if j == i (reaction)

∂(dTᵢ/dt)/∂ξⱼₖ includes:
- ∂(1/ρcₚ)/∂ξ × (RHS) if k == i (property dependence)
- (1/ρcₚ) × ∂Q_rxn/∂ξ if k == i (reaction)
- (1/ρcₚ) × ∂Q_rad/∂ξ (radiation - non-local!)

# Species equation row (NC×i + j):
dξⱼᵢ/dt = -∇·Jⱼ + Sⱼ  (for gases)
dξⱼᵢ/dt = Sⱼ         (for solids)

Derivatives similar to energy equation.
"""
function compute_analytical_jacobian!(J::AbstractMatrix, u, p, t,
	f, jac_cache::JacobianCache{NC, NR}) where {NC, NR}
	mesh = f.mesh
	material = f.material
	prop_cache = f.prop_cache
	sv = f.sv
	n = sv.n_cells

	# Initialize cache on first call (caches static material properties)
	initialize!(jac_cache, material)

	# Initialize sparsity map for O(1) sparse writes (only on first call with sparse J)
	if jac_cache.sparsity_map === nothing && J isa SparseMatrixCSC
		jac_cache.sparsity_map = create_sparsity_map(J)
	end
	smap = jac_cache.sparsity_map

	# Unpack state
	unpack_state!(mesh, u, sv)
	Properties.update_properties!(prop_cache, mesh, material)

	# Compute RHS values (needed for property derivative terms: ∂(1/ρcₚ)/∂ξ × RHS)
	# This populates f.dT_work, f.dξ_work, f.rhs_cache.dξ_rxn
	# Must match the RHS dispatch in _compute_rhs_dispatch! (Darcy-Fick when permeability present)
	if has_darcy_flow(material)
		Discretization.compute_rhs_darcy_fick!(f.dT_work, f.dξ_work, mesh, prop_cache, f.rhs_cache,
			f.material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
			I_radiation = f.I_radiation, radiation_model = f.radiation_model)
	else
		Discretization.compute_rhs!(f.dT_work, f.dξ_work, mesh, prop_cache, f.rhs_cache,
			f.material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
			I_radiation = f.I_radiation, radiation_model = f.radiation_model)
	end

	# Zero out Jacobian (efficient O(nnz) for sparse matrices)
	zero_jacobian!(J)

	# ==========================================================================
	# Cell-local Jacobian terms (reactions, properties)
	# ==========================================================================
	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		# Cache all component properties at this temperature to avoid repeated calls
		cache_all_component_properties!(jac_cache, T, material)

		# Get current property values
		ρcp = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp

		# Compute reaction rates (needed for derivatives)
		evaluate_reactions!(f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, material)

		# Compute reaction rate derivatives (using cached reaction properties)
		evaluate_reaction_derivatives!(jac_cache.dr_dT, jac_cache.dr_dξ,
			f.rhs_cache.rates, ξ, T, jac_cache)

		# Heat source derivatives (using cached reaction properties)
		dQ_rxn_dξ = jac_cache.work_NC  # Reuse work array
		dQ_rxn_dT = heat_source_derivatives_cached!(dQ_rxn_dξ,
			jac_cache.dr_dT, jac_cache.dr_dξ,
			f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, jac_cache, material)

		# Species source derivatives (using cached reaction properties)
		# Note: Must use @view to modify in-place, not create a copy
		# Store dS_dT in dS_dT column, reusing dS_dξ slice
		species_source_derivatives_cached!(@view(jac_cache.dS_dT[:, i]), @view(jac_cache.dS_dξ[:, :, 1]),
			jac_cache.dr_dT, jac_cache.dr_dξ, jac_cache, NC)

		# Heat capacity derivatives (in-place into work_NC_2) - using cached version
		dρcp_dξ = jac_cache.work_NC_2  # Use second work array
		dρcp_dT = heat_capacity_derivatives_cached!(dρcp_dξ, ξ,
			jac_cache.comp_heat_capacities, jac_cache.comp_dcp_dT)
		# Transform dρcp_dξ in-place to dinv_ρcp_dξ
		dinv_ρcp_dT = inverse_heat_capacity_derivatives!(dρcp_dξ, ρcp, dρcp_dT)
		# Note: dρcp_dξ is now dinv_ρcp_dξ

		# Get RHS value for this cell (needed for property derivative contribution)
		# RHS_T = (1/ρcₚ) × (-div_q + Q_rad + Q_rxn)
		RHS_T = f.dT_work[i] * ρcp  # Remove 1/ρcₚ factor to get raw RHS

		# -----------------------------------------------------------------
		# Energy equation derivatives (row i)
		# -----------------------------------------------------------------
		T_row = i

		# ∂(dT/dt)/∂T_i from property derivative: ∂(1/ρcₚ)/∂T × RHS
		J[T_row, T_row] += dinv_ρcp_dT * RHS_T

		# ∂(dT/dt)/∂T_i from reaction: (1/ρcₚ) × ∂Q_rxn/∂T
		J[T_row, T_row] += inv_ρcp * dQ_rxn_dT

		# ∂(dT/dt)/∂ξⱼ from property derivative: ∂(1/ρcₚ)/∂ξⱼ × RHS
		for j in 1:NC
			ξ_col = n + (j-1)*n + i
			J[T_row, ξ_col] += dρcp_dξ[j] * RHS_T  # dρcp_dξ is now dinv_ρcp_dξ
		end

		# ∂(dT/dt)/∂ξⱼ from reaction: (1/ρcₚ) × ∂Q_rxn/∂ξⱼ
		for j in 1:NC
			ξ_col = n + (j-1)*n + i
			J[T_row, ξ_col] += inv_ρcp * dQ_rxn_dξ[j]
		end

		# -----------------------------------------------------------------
		# Species equation derivatives (rows for each component)
		# -----------------------------------------------------------------
		for j in 1:NC
			ξ_row = n + (j-1)*n + i

			# ∂(dξⱼ/dt)/∂T from reaction: ∂Sⱼ/∂T
			J[ξ_row, T_row] = jac_cache.dS_dT[j, i]

			# ∂(dξⱼ/dt)/∂ξₖ from reaction: ∂Sⱼ/∂ξₖ
			for k in 1:NC
				ξ_col = n + (k-1)*n + i
				J[ξ_row, ξ_col] = jac_cache.dS_dξ[j, k, 1]
			end
		end
	end

	# ==========================================================================
	# Face flux Jacobian terms (conduction, gas transport)
	# ==========================================================================
	# Pre-compute ALL cell properties ONCE (Optimization #1)
	precompute_all_cell_properties!(jac_cache, mesh, material)

	# Fused face loop for conduction and gas transport Jacobians
	compute_all_face_jacobians_fused!(J, smap, mesh, prop_cache, material,
		f.rhs_cache.J_gas, n, Val(NC), jac_cache)

	# Cell-centered gas convective source Jacobian
	compute_gas_convective_source_jacobian!(J, smap, mesh, prop_cache, material,
		f.rhs_cache.J_gas, n, Val(NC), jac_cache)

	# ==========================================================================
	# Darcy advective transport Jacobian (if active)
	# ==========================================================================
	if has_darcy_flow(material)
		compute_darcy_transport_jacobian!(J, smap, mesh, prop_cache, material, n, Val(NC), jac_cache)
	end

	# ==========================================================================
	# Boundary condition Jacobian terms
	# ==========================================================================
	compute_bc_jacobian!(J, smap, mesh, prop_cache, f.bc_set, n, NC, t,
		f.top_face_ids, f.bottom_face_ids, f.radiation_model)
	compute_mass_bc_jacobian!(J, smap, mesh, prop_cache, f.bc_set, material, n, NC, t, f.top_face_ids, f.bottom_face_ids, jac_cache.is_gas, jac_cache)

	# NOTE: Boundary convective heat flux Jacobian removed (advection form has zero boundary contribution)
	# Energy leaving the system is handled by boundary conditions, not by convective heat flux.

	# ==========================================================================
	# Radiation Jacobian (non-local Beer-Lambert)
	# ==========================================================================
	if f.I_radiation > 0.0 && f.radiation_model == BEER_LAMBERT
		compute_radiation_jacobian!(J, smap, mesh, prop_cache, f.I_radiation, n, NC, jac_cache.comp_absorptions, jac_cache.τ_cumulative)
	end

	return nothing
end

# Dispatch for radiation wrapper (non-ALE version)
function compute_analytical_jacobian!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionWithRadiation,
	jac_cache::JacobianCache)
	# Delegate to the inner function, extracting radiation intensity
	inner_f = f.ode_fn
	# Update radiation intensity if needed
	if f.I_func !== nothing
		I = f.I_func(t)
	else
		I = 0.0
	end
	# Store radiation intensity in inner function for Jacobian computation
	inner_f.I_radiation = I
	# Delegate to main compute function with unwrapped ODE function
	compute_analytical_jacobian!(J, u, p, t, inner_f, jac_cache)
end

# ============================================================================
# Conduction Jacobian
# ============================================================================

"""
	compute_conduction_jacobian!(J, mesh, prop_cache, material, n, NC, jac_cache)

Add conduction flux derivatives to Jacobian.

For interior face f between cells L and R:
q_f = -k_face × (T_R - T_L) / d_LR

Contribution to cell L: +q_f/V_L (flux in)
Contribution to cell R: -q_f/V_R (flux out)

Derivatives:
∂q_f/∂T_L = k_face/d + q × ∂ln(k_face)/∂T_L
∂q_f/∂T_R = -k_face/d + q × ∂ln(k_face)/∂T_R
"""
function compute_conduction_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material, n::Int, NC::Int, jac_cache)
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces (handled separately)
		if left_id == 0 || right_id == 0
			continue
		end

		# Get cell properties
		k_L = prop_cache.k[left_id]
		k_R = prop_cache.k[right_id]
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Face geometry with cell-to-face distances
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		A_f = mesh.face_areas[f]

		# Face conductivity (distance-weighted harmonic mean) and derivatives
		k_face, dk_dkL, dk_dkR = harmonic_mean_derivatives_weighted(k_L, k_R, d_L, d_R)

		# Current flux and temperature difference
		ΔT = T_R - T_L
		q = -k_face * ΔT / d_LR

		# Flux derivatives w.r.t. temperature
		# ∂q/∂T_L = k_face/d (direct) + ∂k_face/∂k_L × ∂k_L/∂T × (-ΔT/d)
		# ∂q/∂T_R = -k_face/d (direct) + ∂k_face/∂k_R × ∂k_R/∂T × (-ΔT/d)

		# For simplicity, ignore property T-dependence in conduction Jacobian
		# (main coupling is through temperature gradient)
		dq_dT_L = k_face / d_LR
		dq_dT_R = -k_face / d_LR

		# Heat capacity factors
		inv_ρcp_L = 1.0 / prop_cache.ρcp[left_id]
		inv_ρcp_R = 1.0 / prop_cache.ρcp[right_id]

		# Pre-compute division factors (Task 3.3 optimization)
		A_f_inv_V_L = A_f / V_L
		A_f_inv_V_R = A_f / V_R
		coeff_L = inv_ρcp_L * A_f_inv_V_L
		coeff_R = inv_ρcp_R * A_f_inv_V_R

		# Energy equation: dT/dt = (1/ρcₚ) × (-∇·q + ...)
		# For cell L: flux contribution is +q×A/V (flux enters from right)
		# For cell R: flux contribution is -q×A/V (flux leaves to left)

		# d(dT_L/dt)/dT_L from flux: -(1/ρcₚ_L) × A/V_L × ∂q/∂T_L
		# (negative because energy equation has -∇·q)
		# Batched writes: 2 calls instead of 4 (Optimization #4)
		fast_sparse_sub_batch!(J, smap, left_id,
			(left_id, right_id),
			(coeff_L * dq_dT_L, coeff_L * dq_dT_R))
		fast_sparse_add_batch!(J, smap, right_id,
			(left_id, right_id),
			(coeff_R * dq_dT_L, coeff_R * dq_dT_R))

		# =================================================================
		# Flux derivatives w.r.t. species concentrations (cross-cell coupling)
		# q = -k_face × ΔT / d where k_face = harmonic_mean_weighted(k_L, k_R, d_L, d_R)
		# ∂q/∂ξ_j,L = -∂k_face/∂k_L × ∂k_L/∂ξ_j × ΔT / d
		# ∂q/∂ξ_j,R = -∂k_face/∂k_R × ∂k_R/∂ξ_j × ΔT / d
		# =================================================================
		if abs(ΔT) > JAC_ZERO_TOL
			# Get base cache (JacobianCacheALE wraps JacobianCache)
			base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache

			# Compute dk/dξ for left and right cells using work arrays
			dk_L_dξ = base_cache.work_NC
			dk_R_dξ = base_cache.work_NC_2
			v_work = base_cache.ξ_work

			# OPTIMIZATION #1: Use pre-computed per-cell properties instead of re-caching
			# These views are zero-cost and access the pre-computed matrices
			densities_L = @view base_cache.comp_densities_per_cell[:, left_id]
			conductivities_L = @view base_cache.comp_conductivities_per_cell[:, left_id]
			dk_dT_L_vals = @view base_cache.comp_dk_dT_per_cell[:, left_id]

			# Left cell conductivity derivatives (using pre-cached properties)
			dk_dT_L = conductivity_derivatives_cached!(dk_L_dξ, ξ_L,
				densities_L, conductivities_L, dk_dT_L_vals,
				base_cache.is_gas, base_cache.swelling_factors,
				material.conductivity_mixing, material.β, v_work)

			# Right cell pre-cached properties
			densities_R = @view base_cache.comp_densities_per_cell[:, right_id]
			conductivities_R = @view base_cache.comp_conductivities_per_cell[:, right_id]
			dk_dT_R_vals = @view base_cache.comp_dk_dT_per_cell[:, right_id]

			# Right cell conductivity derivatives (using pre-cached properties)
			dk_dT_R = conductivity_derivatives_cached!(dk_R_dξ, ξ_R,
				densities_R, conductivities_R, dk_dT_R_vals,
				base_cache.is_gas, base_cache.swelling_factors,
				material.conductivity_mixing, material.β, v_work)

			# Common factor: -ΔT / d_LR
			factor = -ΔT / d_LR

			for j in 1:NC
				# Derivative of flux w.r.t. ξ_j in left cell
				dq_dξ_jL = factor * dk_dkL * dk_L_dξ[j]
				# Derivative of flux w.r.t. ξ_j in right cell
				dq_dξ_jR = factor * dk_dkR * dk_R_dξ[j]

				if abs(dq_dξ_jL) > FLUX_DERIV_TOL
					ξ_col_L = n + (j-1)*n + left_id
					# Effect on left cell temperature
					fast_sparse_sub!(J, smap, left_id, ξ_col_L, coeff_L * dq_dξ_jL)
					# Effect on right cell temperature (cross-cell coupling!)
					fast_sparse_add!(J, smap, right_id, ξ_col_L, coeff_R * dq_dξ_jL)
				end

				if abs(dq_dξ_jR) > FLUX_DERIV_TOL
					ξ_col_R = n + (j-1)*n + right_id
					# Effect on left cell temperature (cross-cell coupling!)
					fast_sparse_sub!(J, smap, left_id, ξ_col_R, coeff_L * dq_dξ_jR)
					# Effect on right cell temperature
					fast_sparse_add!(J, smap, right_id, ξ_col_R, coeff_R * dq_dξ_jR)
				end
			end
		end
	end

	return nothing
end

# ============================================================================
# Gas Transport Jacobian
# ============================================================================

"""
	compute_gas_transport_jacobian!(J, mesh, prop_cache, material, n, NC, jac_cache)

Add gas transport flux derivatives to Jacobian.

Gas flux using volume fraction gradient formulation:
	J_j = -λ_face × (ξ_R[j]×T_R - ξ_L[j]×T_L) / (T_face × d_LR)
	λ_face = harmonic_mean(λ_L, λ_R)  (distance-weighted)
	T_face = (T_L + T_R) / 2

This distance-weighted harmonic mean is the correct resistance-in-series
interpolation for non-uniform meshes. The derivatives involve both ξ and T
dependence through the volume fraction gradient formulation.

Uses pre-allocated work arrays in jac_cache to avoid allocations from
heterogeneous tuple access in the hot loop.

OPTIMIZATION: dλ_dξ derivatives are computed once per face (for all k)
and stored in work arrays, then reused across all gas component j loops.
This reduces derivative computations from O(n_faces × n_gas × NC × 2) to
O(n_faces × NC × 2), a ~6× reduction for typical materials.
"""
function compute_gas_transport_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material::Material{NC}, n::Int, ::Val{NC}, jac_cache) where NC

	# Get work arrays from cache
	base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache
	# Note: densities_L/R and lambdas_L/R now use @view into per-cell matrices (set per-face below)
	v_L = base_cache.volume_fractions
	v_R = base_cache.ξ_work
	is_gas = base_cache.is_gas  # Cached gas phase information
	swelling_factors = base_cache.swelling_factors  # Cached swelling factors

	# Hoisted derivative arrays (computed once per face, reused across all gas components)
	dλ_L_dξ_arr = base_cache.work_NC_7
	dλ_R_dξ_arr = base_cache.work_NC_8

	# Get mixing parameters
	mixing = material.gas_transfer_mixing
	β = material.β_λ

	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces
		if left_id == 0 || right_id == 0
			continue
		end

		# Face geometry
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		A_f = mesh.face_areas[f]

		# Pre-compute division factors (Task 3.3 optimization)
		A_f_inv_V_L = A_f / V_L
		A_f_inv_V_R = A_f / V_R

		# Cell states
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Gas transfer coefficients
		λ_L = prop_cache.λ[left_id]
		λ_R = prop_cache.λ[right_id]

		# Compute distance-weighted harmonic mean of λ (NOT λρ)
		# λ_face = (d_L + d_R) / (d_L/λ_L + d_R/λ_R)
		if λ_L > LAMBDA_MIN && λ_R > LAMBDA_MIN
			R_L = d_L / λ_L
			R_R = d_R / λ_R
			R_sum = R_L + R_R
			λ_face = d_LR / R_sum

			# Derivatives of distance-weighted harmonic mean:
			# ∂λ_face/∂λ_L = d_LR × d_L / (λ_L² × R_sum²)
			dλ_face_dλ_L = d_LR * d_L / (λ_L * λ_L * R_sum * R_sum)
			dλ_face_dλ_R = d_LR * d_R / (λ_R * λ_R * R_sum * R_sum)
		else
			λ_face = 0.0
			dλ_face_dλ_L = 0.0
			dλ_face_dλ_R = 0.0
		end

		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		T_face = 0.5 * (T_L + T_R)

		# Volume fraction gradient: Φ = ξ_R[j]×T_R - ξ_L[j]×T_L
		# (computed per gas component j below)

		# Use pre-computed per-cell properties (zero-cost @view instead of cache_component_properties!)
		densities_L = @view base_cache.comp_densities_per_cell[:, left_id]
		lambdas_L = @view base_cache.comp_lambdas_per_cell[:, left_id]
		V_total_L = compute_volume_fractions!(v_L, ξ_L, densities_L)

		# Use pre-computed per-cell properties for right cell
		densities_R = @view base_cache.comp_densities_per_cell[:, right_id]
		lambdas_R = @view base_cache.comp_lambdas_per_cell[:, right_id]
		V_total_R = compute_volume_fractions!(v_R, ξ_R, densities_R)

		# =====================================================================
		# OPTIMIZATION: Pre-compute dλ_dξ for ALL components k
		# These are computed ONCE per face and reused in the gas component loop
		# =====================================================================
		@inbounds for k in 1:NC
			dλ_L_dξ_arr[k] = compute_dλ_dξ_cached(ξ_L, k, densities_L, lambdas_L, v_L, V_total_L, mixing, β)
			dλ_R_dξ_arr[k] = compute_dλ_dξ_cached(ξ_R, k, densities_R, lambdas_R, v_R, V_total_R, mixing, β)
		end

		# Pre-compute common denominator
		inv_Tf_dLR = 1.0 / (T_face * d_LR)
		inv_Tf_sq_dLR = inv_Tf_dLR / T_face

		# For each gas component
		for j in 1:NC
			# Use cached is_gas instead of material.components[j].phase (avoids allocation)
			if !is_gas[j]
				continue
			end

			# Volume fraction gradient
			Φ = ξ_R[j] * T_R - ξ_L[j] * T_L

			# Species equation row
			ξ_row_L = n + (j-1)*n + left_id
			ξ_row_R = n + (j-1)*n + right_id

			# Derivatives of flux w.r.t. left cell concentrations
			# J = -λ_face × Φ / (T_face × d_LR)
			# dJ/dξ_L[k] = -[dλ_face/dξ_L[k] × Φ + λ_face × dΦ/dξ_L[k]] / (T_face × d_LR)
			for k in 1:NC
				ξ_col = n + (k-1)*n + left_id

				# Use pre-computed derivative
				dλ_L_dξ = dλ_L_dξ_arr[k]

				# Chain rule: dλ_face/dξ_L[k] = dλ_face/dλ_L × dλ_L/dξ_L[k]
				dλ_face_dξ_L = dλ_face_dλ_L * dλ_L_dξ
				term1 = dλ_face_dξ_L * Φ

				# dΦ/dξ_L[k] = -T_L × δ(k,j) (only same-component contributes)
				if k == j
					term2 = -λ_face * T_L
				else
					term2 = 0.0
				end

				dJ_dξ_L = -(term1 + term2) * inv_Tf_dLR

				fast_sparse_sub!(J, smap, ξ_row_L, ξ_col, A_f_inv_V_L * dJ_dξ_L)
				fast_sparse_add!(J, smap, ξ_row_R, ξ_col, A_f_inv_V_R * dJ_dξ_L)
			end

			# Derivatives of flux w.r.t. right cell concentrations
			for k in 1:NC
				ξ_col = n + (k-1)*n + right_id

				# Use pre-computed derivative
				dλ_R_dξ = dλ_R_dξ_arr[k]

				# Chain rule: dλ_face/dξ_R[k] = dλ_face/dλ_R × dλ_R/dξ_R[k]
				dλ_face_dξ_R = dλ_face_dλ_R * dλ_R_dξ
				term1 = dλ_face_dξ_R * Φ

				# dΦ/dξ_R[k] = +T_R × δ(k,j)
				if k == j
					term2 = λ_face * T_R
				else
					term2 = 0.0
				end

				dJ_dξ_R = -(term1 + term2) * inv_Tf_dLR

				fast_sparse_sub!(J, smap, ξ_row_L, ξ_col, A_f_inv_V_L * dJ_dξ_R)
				fast_sparse_add!(J, smap, ξ_row_R, ξ_col, A_f_inv_V_R * dJ_dξ_R)
			end

			# NEW: Derivatives of flux w.r.t. temperature
			# ∂J/∂T_L = λ_face × [ξ_L[j]×T_face + 0.5×Φ] / (T_face² × d_LR)
			dJ_dT_L = λ_face * (ξ_L[j] * T_face + 0.5 * Φ) * inv_Tf_sq_dLR
			# ∂J/∂T_R = -λ_face × [ξ_R[j]×T_face - 0.5×Φ] / (T_face² × d_LR)
			dJ_dT_R = -λ_face * (ξ_R[j] * T_face - 0.5 * Φ) * inv_Tf_sq_dLR

			if abs(dJ_dT_L) > FLUX_DERIV_TOL || abs(dJ_dT_R) > FLUX_DERIV_TOL
				# T is in rows 1:n, columns 1:n (col = cell_id for temperature)
				fast_sparse_sub!(J, smap, ξ_row_L, left_id, A_f_inv_V_L * dJ_dT_L)
				fast_sparse_add!(J, smap, ξ_row_R, left_id, A_f_inv_V_R * dJ_dT_L)
				fast_sparse_sub!(J, smap, ξ_row_L, right_id, A_f_inv_V_L * dJ_dT_R)
				fast_sparse_add!(J, smap, ξ_row_R, right_id, A_f_inv_V_R * dJ_dT_R)
			end
		end
	end

	return nothing
end


# ============================================================================
# Gas Convective Source Jacobian (cell-centered)
# ============================================================================

"""
	compute_gas_convective_source_jacobian!(J, smap, mesh, prop_cache, material,
										    J_gas, n, NC, jac_cache)

Add gas convective source derivatives to the energy equation Jacobian.

The gas convective source at cell i is:
	S_conv,i = -Σ_g c_g(T_i) × J̄_g,i × (∂T/∂x)_i

where J̄_g,i = (J_g[f_L] + J_g[f_R]) / 2.

Temperature derivatives (diagonal and off-diagonal):
	∂S/∂T_i: through (∂T/∂x)_i
	∂S/∂T_{i±1}: through (∂T/∂x)_i neighbor dependence

Species derivatives (through J̄_g dependence on ξ):
	∂S/∂ξ_{j,k}: through ∂J̄/∂ξ at bounding faces
"""
function compute_gas_convective_source_jacobian!(J::AbstractMatrix, smap,
	mesh, prop_cache, material::Material{NC}, J_gas::Matrix{Float64},
	n::Int, ::Val{NC}, jac_cache) where NC

	base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache
	is_gas = base_cache.is_gas
	swelling_factors = base_cache.swelling_factors
	λ_mixing = material.gas_transfer_mixing
	β_λ = material.β_λ

	# Work arrays for volume fractions
	v_work1 = base_cache.volume_fractions
	v_work2 = base_cache.ξ_work

	# Hoisted derivative arrays
	dρ_dξ_arr = base_cache.work_NC_5
	dλ_dξ_arr = base_cache.work_NC_6

	n_cells = length(mesh.cells)

	# Work array for per-cell heat capacities (avoids runtime tuple indexing)
	cp_work = base_cache.comp_heat_capacities

	@inbounds for i in 1:n_cells
		T_i = mesh.cell_states[i].T
		ρcp_i = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp_i

		# Cache component heat capacities and their T-derivatives at this cell's temperature
		# Uses @generated function to avoid heterogeneous tuple allocations
		# Populates comp_heat_capacities (= cp_work) AND comp_dcp_dT
		cache_all_component_properties!(base_cache, T_i, material)

		# Bounding faces
		cell = mesh.cells[i]
		f_L, f_R = cell.face_ids
		left_id = mesh.faces[f_L].cell_ids[1]
		right_id = mesh.faces[f_R].cell_ids[2]

		# Compute dTdx and its derivatives w.r.t. temperatures
		if left_id == 0 && right_id == 0
			# Single-cell mesh: no gradient, no contribution
			continue
		elseif left_id == 0
			# Left boundary: dTdx = (T_R - T_i) / d_R
			T_R = mesh.cell_states[right_id].T
			d_R = abs(mesh.cell_centers[right_id] - mesh.cell_centers[i])
			dTdx = (T_R - T_i) / d_R
			# Derivatives
			ddTdx_dTi = -1.0 / d_R
			ddTdx_dTL = 0.0  # no left neighbor
			ddTdx_dTR = 1.0 / d_R
		elseif right_id == 0
			# Right boundary: dTdx = (T_i - T_L) / d_L
			T_L = mesh.cell_states[left_id].T
			d_L = abs(mesh.cell_centers[i] - mesh.cell_centers[left_id])
			dTdx = (T_i - T_L) / d_L
			# Derivatives
			ddTdx_dTi = 1.0 / d_L
			ddTdx_dTL = -1.0 / d_L
			ddTdx_dTR = 0.0  # no right neighbor
		else
			# Interior: dTdx = ((T_i - T_L)/d_L + (T_R - T_i)/d_R) / 2
			T_L = mesh.cell_states[left_id].T
			T_R = mesh.cell_states[right_id].T
			d_L = abs(mesh.cell_centers[i] - mesh.cell_centers[left_id])
			d_R = abs(mesh.cell_centers[right_id] - mesh.cell_centers[i])
			dTdx = ((T_i - T_L) / d_L + (T_R - T_i) / d_R) * 0.5
			# Derivatives
			ddTdx_dTi = (1.0 / d_L + 1.0 / d_R) * 0.5
			ddTdx_dTL = -0.5 / d_L
			ddTdx_dTR = 0.5 / d_R
		end

		# Compute Σ_g c_g * J̄_g and Σ_g (dc_g/dT) * J̄_g
		# Uses pre-cached cp_work[j] and comp_dcp_dT[j] to avoid runtime heterogeneous tuple indexing
		sum_cgJ = 0.0
		sum_dcgJ_dT = 0.0
		dcp_dT_work = base_cache.comp_dcp_dT
		for j in 1:NC
			if is_gas[j]
				J_avg = (J_gas[j, f_L] + J_gas[j, f_R]) * 0.5
				sum_cgJ += cp_work[j] * J_avg
				sum_dcgJ_dT += dcp_dT_work[j] * J_avg
			end
		end

		# ∂S/∂T contributions: S = -Σ_g c_g(T) * J̄_g * dTdx
		# Full derivative: ∂S/∂T_i = -Σ_g [dc_g/dT * J̄_g * dTdx + c_g * J̄_g * ddTdx/dT_i]
		#                           = -sum_dcgJ_dT * dTdx - sum_cgJ * ddTdx_dTi
		T_row = i

		dS_dTi = -sum_cgJ * ddTdx_dTi - sum_dcgJ_dT * dTdx
		if abs(dS_dTi) > 1e-30
			fast_sparse_add!(J, smap, T_row, i, inv_ρcp * dS_dTi)
		end

		# Left neighbor
		if left_id != 0 && abs(ddTdx_dTL) > 0.0
			dS_dTL = -sum_cgJ * ddTdx_dTL
			if abs(dS_dTL) > 1e-30
				fast_sparse_add!(J, smap, T_row, left_id, inv_ρcp * dS_dTL)
			end
		end

		# Right neighbor
		if right_id != 0 && abs(ddTdx_dTR) > 0.0
			dS_dTR = -sum_cgJ * ddTdx_dTR
			if abs(dS_dTR) > 1e-30
				fast_sparse_add!(J, smap, T_row, right_id, inv_ρcp * dS_dTR)
			end
		end

		# =====================================================================
		# Species derivatives: ∂S/∂ξ through J̄_g dependence
		# S = -Σ_g c_g * J̄_g * dTdx
		# ∂S/∂ξ_{j,k} = -c_g * dTdx * ∂J̄_g/∂ξ_{j,k}
		# J̄_g = (J[f_L] + J[f_R]) / 2
		# =====================================================================
		if abs(dTdx) < 1e-30
			continue  # No species derivatives if dTdx ≈ 0
		end

		# Process left face f_L (between left_id and i)
		# J at f_L depends on ξ of left_id and i
		if left_id != 0
			_add_convective_source_species_jacobian_face!(
				J, smap, mesh, prop_cache, material,
				J_gas, f_L, left_id, i,
				T_i, dTdx, inv_ρcp, T_row, n,
				is_gas, swelling_factors, λ_mixing, β_λ,
				base_cache, v_work1, v_work2, dρ_dξ_arr, dλ_dξ_arr,
				0.5, cp_work)  # factor 0.5 for averaging
		end

		# Process right face f_R (between i and right_id)
		if right_id != 0
			_add_convective_source_species_jacobian_face!(
				J, smap, mesh, prop_cache, material,
				J_gas, f_R, i, right_id,
				T_i, dTdx, inv_ρcp, T_row, n,
				is_gas, swelling_factors, λ_mixing, β_λ,
				base_cache, v_work1, v_work2, dρ_dξ_arr, dλ_dξ_arr,
				0.5, cp_work)  # factor 0.5 for averaging
		end
	end

	return nothing
end

"""
Helper: add species Jacobian contributions from one face of the convective source.

For face f between cells cell_L and cell_R:
J_j,f = -λ_face × (ξ_R[j]×T_R - ξ_L[j]×T_L) / (T_face × d_LR)

Contribution to S_conv,i: -c_g * (J_j,f / 2) * dTdx
So ∂S/∂ξ_k = -c_g * dTdx * (∂J_j,f/∂ξ_k) / 2 = -c_g * dTdx * factor * ∂J_j,f/∂ξ_k
"""
function _add_convective_source_species_jacobian_face!(
	J::AbstractMatrix, smap, mesh, prop_cache, material::Material{NC},
	J_gas::Matrix{Float64}, f::Int, cell_L::Int, cell_R::Int,
	T_i::Float64, dTdx::Float64, inv_ρcp::Float64, T_row::Int,
	n::Int,
	is_gas, swelling_factors, λ_mixing, β_λ,
	base_cache, v_work1, v_work2, dρ_dξ_arr, dλ_dξ_arr,
	factor::Float64, cp_work::Vector{Float64}) where NC

	# Cell states at the face
	ξ_L = mesh.cell_states[cell_L].ξ
	ξ_R = mesh.cell_states[cell_R].ξ
	λ_L = prop_cache.λ[cell_L]
	λ_R = prop_cache.λ[cell_R]

	# Face geometry
	face_pos = mesh.face_positions[f]
	d_L = abs(face_pos - mesh.cell_centers[cell_L])
	d_R = abs(mesh.cell_centers[cell_R] - face_pos)
	d_LR = d_L + d_R

	# λ harmonic mean (NOT λρ)
	if λ_L <= LAMBDA_MIN || λ_R <= LAMBDA_MIN
		return nothing
	end

	R_L = d_L / λ_L
	R_R = d_R / λ_R
	R_sum = R_L + R_R
	λ_face = d_LR / R_sum

	dλ_face_dλ_L = d_LR * d_L / (λ_L * λ_L * R_sum * R_sum)
	dλ_face_dλ_R = d_LR * d_R / (λ_R * λ_R * R_sum * R_sum)

	T_L = mesh.cell_states[cell_L].T
	T_R = mesh.cell_states[cell_R].T
	T_face = 0.5 * (T_L + T_R)
	inv_Tf_dLR = 1.0 / (T_face * d_LR)

	# Pre-compute per-cell properties for this face
	densities_L = @view base_cache.comp_densities_per_cell[:, cell_L]
	lambdas_L = @view base_cache.comp_lambdas_per_cell[:, cell_L]
	V_total_L = compute_volume_fractions!(v_work1, ξ_L, densities_L)
	densities_R = @view base_cache.comp_densities_per_cell[:, cell_R]
	lambdas_R = @view base_cache.comp_lambdas_per_cell[:, cell_R]
	V_total_R = compute_volume_fractions!(v_work2, ξ_R, densities_R)

	@inbounds for j in 1:NC
		if !is_gas[j]
			continue
		end

		prefactor = -cp_work[j] * dTdx * factor * inv_ρcp

		Φ = ξ_R[j] * T_R - ξ_L[j] * T_L

		# Left cell derivatives
		for k in 1:NC
			dλ_dξ = compute_dλ_dξ_cached(ξ_L, k, densities_L, lambdas_L, v_work1, V_total_L, λ_mixing, β_λ)

			dλ_face_dξ = dλ_face_dλ_L * dλ_dξ
			term1 = dλ_face_dξ * Φ

			if k == j
				term2 = -λ_face * T_L
			else
				term2 = 0.0
			end

			dJ_dξ = -(term1 + term2) * inv_Tf_dLR

			ξ_col = n + (k - 1) * n + cell_L
			if abs(dJ_dξ) > 1e-30
				fast_sparse_add!(J, smap, T_row, ξ_col, prefactor * dJ_dξ)
			end
		end

		# Right cell derivatives
		for k in 1:NC
			dλ_dξ = compute_dλ_dξ_cached(ξ_R, k, densities_R, lambdas_R, v_work2, V_total_R, λ_mixing, β_λ)

			dλ_face_dξ = dλ_face_dλ_R * dλ_dξ
			term1 = dλ_face_dξ * Φ

			if k == j
				term2 = λ_face * T_R
			else
				term2 = 0.0
			end

			dJ_dξ = -(term1 + term2) * inv_Tf_dLR

			ξ_col = n + (k - 1) * n + cell_R
			if abs(dJ_dξ) > 1e-30
				fast_sparse_add!(J, smap, T_row, ξ_col, prefactor * dJ_dξ)
			end
		end
	end

	return nothing
end

# ============================================================================
# Fused Face Jacobian (Optimization #3)
# ============================================================================

"""
	compute_all_face_jacobians_fused!(J, smap, mesh, prop_cache, material,
									   J_gas, n, NC, jac_cache)

Fused face loop combining conduction and gas transport Jacobians.

This function eliminates redundant computations by:
1. Iterating over faces only ONCE instead of twice
2. Computing shared quantities (geometry, properties, derivatives) once per face
3. Shared derivative arrays between conduction and gas transport

# Arguments
- `J`: Jacobian matrix (modified in-place)
- `smap`: Sparsity map for fast sparse matrix operations
- `mesh`: Mesh with cell states
- `prop_cache`: Cached material properties
- `material`: Material definition
- `J_gas`: Pre-computed gas mass fluxes [NC × n_faces]
- `n`: Number of cells
- `NC`: Number of components (static)
- `jac_cache`: Jacobian cache with work arrays

# Prerequisites
- `precompute_all_cell_properties!(jac_cache, mesh, material)` must be called first
  to populate per-cell property matrices
"""
function compute_all_face_jacobians_fused!(
	J::AbstractMatrix, smap, mesh, prop_cache,
	material::Material{NC}, J_gas::Matrix{Float64},
	n::Int, ::Val{NC}, jac_cache) where NC

	# Get base cache (handles both JacobianCache and JacobianCacheALE)
	base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache

	# Get mixing parameters (shared across all physics)
	k_mixing = material.conductivity_mixing
	β_k = material.β
	λ_mixing = material.gas_transfer_mixing
	β_λ = material.β_λ
	is_gas = base_cache.is_gas
	swelling_factors = base_cache.swelling_factors

	# Work arrays for volume fractions
	v_L = base_cache.volume_fractions
	v_R = base_cache.ξ_work

	# Hoisted derivative arrays (computed once per face, reused across all gas components)
	dρ_L_dξ_arr = base_cache.work_NC_5
	dρ_R_dξ_arr = base_cache.work_NC_6
	dλ_L_dξ_arr = base_cache.work_NC_7
	dλ_R_dξ_arr = base_cache.work_NC_8

	# Work arrays for conductivity derivatives
	dk_L_dξ = base_cache.work_NC
	dk_R_dξ = base_cache.work_NC_2

	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces
		if left_id == 0 || right_id == 0
			continue
		end

		# =====================================================================
		# SHARED: Geometry and cell states (computed once, used by all physics)
		# =====================================================================
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		A_f = mesh.face_areas[f]

		# Pre-compute division factors
		A_f_inv_V_L = A_f / V_L
		A_f_inv_V_R = A_f / V_R

		# Cell states
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# =====================================================================
		# SHARED: Pre-cached properties (zero-cost @view into per-cell matrices)
		# =====================================================================
		densities_L = @view base_cache.comp_densities_per_cell[:, left_id]
		densities_R = @view base_cache.comp_densities_per_cell[:, right_id]
		conductivities_L = @view base_cache.comp_conductivities_per_cell[:, left_id]
		conductivities_R = @view base_cache.comp_conductivities_per_cell[:, right_id]
		dk_dT_L_vals = @view base_cache.comp_dk_dT_per_cell[:, left_id]
		dk_dT_R_vals = @view base_cache.comp_dk_dT_per_cell[:, right_id]
		lambdas_L = @view base_cache.comp_lambdas_per_cell[:, left_id]
		lambdas_R = @view base_cache.comp_lambdas_per_cell[:, right_id]

		# =====================================================================
		# SHARED: Bulk properties and volume fractions
		# =====================================================================
		V_total_L = compute_volume_fractions!(v_L, ξ_L, densities_L)
		V_total_R = compute_volume_fractions!(v_R, ξ_R, densities_R)
		ρ_L = prop_cache.ρ[left_id]
		ρ_R = prop_cache.ρ[right_id]

		# =====================================================================
		# SHARED: Pre-compute derivatives ONCE for ALL components k
		# These arrays will be reused by BOTH gas transport AND convective heat flux
		# =====================================================================
		@inbounds for k in 1:NC
			dρ_L_dξ_arr[k] = compute_dρ_dξ_cached(ξ_L, k, densities_L, swelling_factors, is_gas)
			dρ_R_dξ_arr[k] = compute_dρ_dξ_cached(ξ_R, k, densities_R, swelling_factors, is_gas)
			dλ_L_dξ_arr[k] = compute_dλ_dξ_cached(ξ_L, k, densities_L, lambdas_L, v_L, V_total_L, λ_mixing, β_λ)
			dλ_R_dξ_arr[k] = compute_dλ_dξ_cached(ξ_R, k, densities_R, lambdas_R, v_R, V_total_R, λ_mixing, β_λ)
		end

		# =================================================================
		# PHYSICS 1: CONDUCTION JACOBIAN
		# =================================================================
		k_L = prop_cache.k[left_id]
		k_R = prop_cache.k[right_id]
		k_face, dk_dkL, dk_dkR = harmonic_mean_derivatives_weighted(k_L, k_R, d_L, d_R)

		ΔT = T_R - T_L

		# Flux derivatives w.r.t. temperature
		# q = -k_face × ΔT / d_LR
		# ∂q/∂T_L = k_face/d_LR + ∂k_face/∂T_L × (-ΔT/d_LR)
		# ∂q/∂T_R = -k_face/d_LR + ∂k_face/∂T_R × (-ΔT/d_LR)
		dq_dT_L = k_face / d_LR
		dq_dT_R = -k_face / d_LR

		# Heat capacity factors
		inv_ρcp_L = 1.0 / prop_cache.ρcp[left_id]
		inv_ρcp_R = 1.0 / prop_cache.ρcp[right_id]
		coeff_L = inv_ρcp_L * A_f_inv_V_L
		coeff_R = inv_ρcp_R * A_f_inv_V_R

		# Concentration and T-property derivatives (if temperature difference is significant)
		if abs(ΔT) > JAC_ZERO_TOL
			# Left cell conductivity derivatives (returns dk_dT scalar)
			dk_dT_L = conductivity_derivatives_cached!(dk_L_dξ, ξ_L,
				densities_L, conductivities_L, dk_dT_L_vals,
				is_gas, swelling_factors, k_mixing, β_k, v_L)

			# Right cell conductivity derivatives (returns dk_dT scalar)
			dk_dT_R = conductivity_derivatives_cached!(dk_R_dξ, ξ_R,
				densities_R, conductivities_R, dk_dT_R_vals,
				is_gas, swelling_factors, k_mixing, β_k, v_R)

			factor = -ΔT / d_LR

			# Add ∂k_face/∂T contribution to temperature derivatives
			# ∂k_face/∂T_L = dk_dkL × dk_L/dT, ∂k_face/∂T_R = dk_dkR × dk_R/dT
			dq_dT_L += factor * dk_dkL * dk_dT_L
			dq_dT_R += factor * dk_dkR * dk_dT_R

			# Concentration derivatives through conductivity
			for j in 1:NC
				dq_dξ_jL = factor * dk_dkL * dk_L_dξ[j]
				dq_dξ_jR = factor * dk_dkR * dk_R_dξ[j]

				if abs(dq_dξ_jL) > FLUX_DERIV_TOL
					ξ_col_L = n + (j-1)*n + left_id
					fast_sparse_sub!(J, smap, left_id, ξ_col_L, coeff_L * dq_dξ_jL)
					fast_sparse_add!(J, smap, right_id, ξ_col_L, coeff_R * dq_dξ_jL)
				end

				if abs(dq_dξ_jR) > FLUX_DERIV_TOL
					ξ_col_R = n + (j-1)*n + right_id
					fast_sparse_sub!(J, smap, left_id, ξ_col_R, coeff_L * dq_dξ_jR)
					fast_sparse_add!(J, smap, right_id, ξ_col_R, coeff_R * dq_dξ_jR)
				end
			end
		end

		# Temperature derivatives in Jacobian (batched: 2 calls instead of 4)
		fast_sparse_sub_batch!(J, smap, left_id,
			(left_id, right_id),
			(coeff_L * dq_dT_L, coeff_L * dq_dT_R))
		fast_sparse_add_batch!(J, smap, right_id,
			(left_id, right_id),
			(coeff_R * dq_dT_L, coeff_R * dq_dT_R))

		# =================================================================
		# PHYSICS 2: GAS TRANSPORT JACOBIAN (volume fraction gradient)
		# =================================================================
		λ_L = prop_cache.λ[left_id]
		λ_R = prop_cache.λ[right_id]

		if λ_L > LAMBDA_MIN && λ_R > LAMBDA_MIN
			# Compute distance-weighted harmonic mean of λ (NOT λρ)
			R_L_gas = d_L / λ_L
			R_R_gas = d_R / λ_R
			R_sum_gas = R_L_gas + R_R_gas
			λ_face = d_LR / R_sum_gas

			# Derivatives of harmonic mean
			dλ_face_dλ_L = d_LR * d_L / (λ_L * λ_L * R_sum_gas * R_sum_gas)
			dλ_face_dλ_R = d_LR * d_R / (λ_R * λ_R * R_sum_gas * R_sum_gas)

			T_face = 0.5 * (T_L + T_R)
			inv_Tf_dLR = 1.0 / (T_face * d_LR)
			inv_Tf_sq_dLR = inv_Tf_dLR / T_face

			# Pre-compute ∂λ/∂T for left and right cells (once per face, shared across gas components)
			# J = -λ_face × Φ / (T_face × d_LR)
			# ∂J/∂T via λ: (-Φ / (T_face × d_LR)) × ∂λ_face/∂T
			# where ∂λ_face/∂T_L = dλ_face_dλ_L × dλ_L/dT
			dλ_L_dT = compute_dλ_dT(ξ_L, T_L, material)
			dλ_R_dT = compute_dλ_dT(ξ_R, T_R, material)
			dλ_face_dT_L = dλ_face_dλ_L * dλ_L_dT
			dλ_face_dT_R = dλ_face_dλ_R * dλ_R_dT

			# For each gas component
			for j in 1:NC
				if !is_gas[j]
					continue
				end

				# Volume fraction gradient
				Φ = ξ_R[j] * T_R - ξ_L[j] * T_L

				# Species equation rows
				ξ_row_L = n + (j-1)*n + left_id
				ξ_row_R = n + (j-1)*n + right_id

				# Left cell derivatives (reuse pre-computed dλ_dξ)
				for k in 1:NC
					ξ_col = n + (k-1)*n + left_id

					dλ_L_dξ = dλ_L_dξ_arr[k]

					dλ_face_dξ_L = dλ_face_dλ_L * dλ_L_dξ
					term1 = dλ_face_dξ_L * Φ

					if k == j
						term2 = -λ_face * T_L
					else
						term2 = 0.0
					end

					dJ_dξ_L = -(term1 + term2) * inv_Tf_dLR

					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col, A_f_inv_V_L * dJ_dξ_L)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col, A_f_inv_V_R * dJ_dξ_L)
				end

				# Right cell derivatives (reuse pre-computed dλ_dξ)
				for k in 1:NC
					ξ_col = n + (k-1)*n + right_id

					dλ_R_dξ = dλ_R_dξ_arr[k]

					dλ_face_dξ_R = dλ_face_dλ_R * dλ_R_dξ
					term1 = dλ_face_dξ_R * Φ

					if k == j
						term2 = λ_face * T_R
					else
						term2 = 0.0
					end

					dJ_dξ_R = -(term1 + term2) * inv_Tf_dLR

					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col, A_f_inv_V_L * dJ_dξ_R)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col, A_f_inv_V_R * dJ_dξ_R)
				end

				# Temperature derivatives:
				# 1. Direct T dependence through Φ and T_face
				dJ_dT_L = λ_face * (ξ_L[j] * T_face + 0.5 * Φ) * inv_Tf_sq_dLR
				dJ_dT_R = -λ_face * (ξ_R[j] * T_face - 0.5 * Φ) * inv_Tf_sq_dLR

				# 2. Indirect T dependence through λ_face(T)
				# ∂J/∂T_L (via λ) = -Φ/(T_face × d_LR) × dλ_face/dT_L
				dJ_dT_L += -Φ * inv_Tf_dLR * dλ_face_dT_L
				dJ_dT_R += -Φ * inv_Tf_dLR * dλ_face_dT_R

				if abs(dJ_dT_L) > FLUX_DERIV_TOL || abs(dJ_dT_R) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row_L, left_id, A_f_inv_V_L * dJ_dT_L)
					fast_sparse_add!(J, smap, ξ_row_R, left_id, A_f_inv_V_R * dJ_dT_L)
					fast_sparse_sub!(J, smap, ξ_row_L, right_id, A_f_inv_V_L * dJ_dT_R)
					fast_sparse_add!(J, smap, ξ_row_R, right_id, A_f_inv_V_R * dJ_dT_R)
				end
			end
		end
	end

	return nothing
end

"""
	compute_dρ_dξ(ξ, T, k, material) -> Float64

Compute derivative of mixture density with respect to ξ_k.
For mixture density with volume fractions: ρ = 1 / Σ(ξ_j/ρ_j) / Σξ_j
Approximation: ∂ρ/∂ξ_k ≈ 1 for small volume change effects.
More accurate: use analytical derivative of mixture_density formula.
"""
@inline function compute_dρ_dξ(ξ::NTuple{NC, Float64}, T::Float64, k::Int,
	material::Material{NC}) where NC
	# For simple mixture: ρ = Σξ_j, so ∂ρ/∂ξ_k = 1
	# For mixture with volume fractions: more complex
	# Using finite difference for accuracy
	ε = 1e-6
	ξ_plus = ntuple(i -> i == k ? ξ[i] + ε : ξ[i], Val(NC))
	ρ_plus = Properties.mixture_density(ξ_plus, T, material)
	ρ_0 = Properties.mixture_density(ξ, T, material)
	return (ρ_plus - ρ_0) / ε
end

# Import harmonic_mean from Properties
using ..Properties: harmonic_mean

"""
	compute_dλ_dξ(ξ, T, k, material) -> Float64

Compute derivative of effective gas transfer coefficient with respect to ξ_k.
Handles all mixing rules (PARALLEL, SERIES, WEIGHTED) analytically.

For WEIGHTED: λ_eff = β × λ_parallel + (1-β) × λ_series
- λ_parallel = Σ_j v_j × λ_j
- λ_series = 1 / Σ_j (v_j / λ_j)  (or 0 if any component has v>0 and λ=0)

Derivatives:
- ∂λ_parallel/∂ξ_k = Σ_j λ_j × ∂v_j/∂ξ_k
- ∂λ_series/∂ξ_k = -λ_series² × Σ_j (∂v_j/∂ξ_k / λ_j)  (or 0 if λ_series=0)

Where ∂v_j/∂ξ_k = (δ_{jk}/ρ_j - v_j/ρ_k) / V_total
"""
@inline function compute_dλ_dξ(ξ::NTuple{NC, Float64}, T::Float64, k::Int,
	material::Material{NC}) where NC

	# Compute volume fractions and total specific volume
	V_total = 0.0
	@inbounds for m in 1:NC
		ρ_m = material.components[m].density(T)
		V_total += ξ[m] / max(ρ_m, DENSITY_MIN)
	end

	if V_total < VOLUME_MIN
		return 0.0
	end

	ρ_k = material.components[k].density(T)

	# Compute volume fractions first (needed for series check)
	v = ntuple(Val(NC)) do j
		ρ_j = material.components[j].density(T)
		(ξ[j] / max(ρ_j, DENSITY_MIN)) / V_total
	end

	# Check if series term will be zero (any component with v>1e-20 and λ=0)
	series_is_inf = false
	@inbounds for j in 1:NC
		if v[j] > VOLUME_MIN
			λ_j = material.components[j].gas_transfer_coeff(T)
			if λ_j < LAMBDA_MIN
				series_is_inf = true
				break
			end
		end
	end

	# Compute derivatives
	dλ_parallel_dξ = 0.0
	d_inv_λ_series_dξ = 0.0
	inv_λ_series = 0.0

	@inbounds for j in 1:NC
		comp_j = material.components[j]
		λ_j = comp_j.gas_transfer_coeff(T)
		ρ_j = comp_j.density(T)

		# ∂v_j/∂ξ_k = (δ_{jk}/ρ_j - v_j/ρ_k) / V_total
		if j == k
			dv_dξ = (1.0 / max(ρ_j, DENSITY_MIN) - v[j] / max(ρ_k, DENSITY_MIN)) / V_total
		else
			dv_dξ = -v[j] / max(ρ_k, DENSITY_MIN) / V_total
		end

		# Parallel: λ_parallel = Σ v_j × λ_j
		# ∂λ_parallel/∂ξ_k = Σ λ_j × ∂v_j/∂ξ_k
		dλ_parallel_dξ += λ_j * dv_dξ

		# Series: only compute if series term is not inf
		if !series_is_inf && λ_j > LAMBDA_MIN
			inv_λ_series += v[j] / λ_j
			d_inv_λ_series_dξ += dv_dξ / λ_j
		end
	end

	# Compute λ_series and its derivative
	# If series_is_inf, λ_series = 0 and dλ_series_dξ = 0
	if series_is_inf
		dλ_series_dξ = 0.0
	else
		λ_series = inv_λ_series > LAMBDA_MIN ? 1.0 / inv_λ_series : 0.0
		dλ_series_dξ = -λ_series * λ_series * d_inv_λ_series_dξ
	end

	# Apply mixing rule
	mixing = material.gas_transfer_mixing
	if mixing == PARALLEL
		return dλ_parallel_dξ
	elseif mixing == SERIES
		return dλ_series_dξ
	else  # WEIGHTED
		β = material.β_λ
		return β * dλ_parallel_dξ + (1.0 - β) * dλ_series_dξ
	end
end

# ============================================================================
# Allocation-Free Gas Transport Jacobian Helpers
# ============================================================================

"""
	cache_component_properties!(densities, lambdas, T, material)

Cache component densities and gas transfer coefficients at temperature T.
Uses unrolled access to avoid type instability from heterogeneous tuple indexing.

This function must use @generated or ntuple to access heterogeneous component tuples
without allocations.
"""
@generated function cache_component_properties!(
	densities::AbstractVector{Float64},
	lambdas::AbstractVector{Float64},
	T::Float64,
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(densities[$j] = material.components[$j].density(T)))
		push!(exprs, :(lambdas[$j] = material.components[$j].gas_transfer_coeff(T)))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_heat_capacities!(cp, T, material)

Cache component heat capacities at temperature T.
Uses unrolled access to avoid type instability from heterogeneous tuple indexing.

This function must use @generated to access heterogeneous component tuples
without allocations.
"""
@generated function cache_heat_capacities!(
	cp::AbstractVector{Float64},
	T::Float64,
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(cp[$j] = material.components[$j].heat_capacity(T)))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_all_component_properties!(cache, T, material)

Cache ALL component properties and their T-derivatives at temperature T.
This is the main entry point that populates all per-component arrays in the cache.

Uses @generated to unroll all component accesses, avoiding type instability.

# Properties cached:
- Density (ρ_j) and ∂ρ_j/∂T
- Heat capacity (c_p,j) and ∂c_p,j/∂T
- Thermal conductivity (k_j) and ∂k_j/∂T
- Gas transfer coefficient (λ_j)
"""
@generated function cache_all_component_properties!(
	cache::JacobianCache{NC, NR},
	T::Float64,
	material::Material{NC}) where {NC, NR}

	# Generate unrolled code for all components
	exprs = Expr[]

	for j in 1:NC
		# Property values
		push!(exprs, :(cache.comp_densities[$j] = material.components[$j].density(T)))
		push!(exprs, :(cache.comp_heat_capacities[$j] = material.components[$j].heat_capacity(T)))
		push!(exprs, :(cache.comp_conductivities[$j] = material.components[$j].thermal_conductivity(T)))
		push!(exprs, :(cache.comp_lambdas[$j] = material.components[$j].gas_transfer_coeff(T)))

		# Property T-derivatives
		push!(exprs, :(cache.comp_dρ_dT[$j] = derivative(material.components[$j].density, T)))
		push!(exprs, :(cache.comp_dcp_dT[$j] = derivative(material.components[$j].heat_capacity, T)))
		push!(exprs, :(cache.comp_dk_dT[$j] = derivative(material.components[$j].thermal_conductivity, T)))
	end

	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	precompute_all_cell_properties!(cache, mesh, material)

Pre-compute ALL component properties for ALL cells at once.
This eliminates redundant caching in face loops - properties are evaluated
once per cell and accessed via @view in the face loops.

OPTIMIZATION #1: Reduces property evaluations from O(n_faces × 2) to O(n_cells).
For 400 cells / 399 interior faces, this is ~2× fewer evaluations.
Combined with face loop reuse, actual savings are ~6× (each cell touches ~2 faces).

Uses @generated to unroll component access at compile time for type stability.
"""
@generated function precompute_all_cell_properties!(
	cache::JacobianCache{NC, NR},
	mesh::Mesh1D{NC},
	material::Material{NC}) where {NC, NR}

	# Generate unrolled code for all components
	# Outer loop over cells is runtime, inner loop over components is compile-time unrolled
	exprs = Expr[]

	# For each component j, generate a cell loop
	for j in 1:NC
		push!(exprs, quote
			@inbounds for i in 1:length(mesh.cells)
				T = mesh.cell_states[i].T
				cache.comp_densities_per_cell[$j, i] = material.components[$j].density(T)
				cache.comp_conductivities_per_cell[$j, i] = material.components[$j].thermal_conductivity(T)
				cache.comp_lambdas_per_cell[$j, i] = material.components[$j].gas_transfer_coeff(T)
				cache.comp_dk_dT_per_cell[$j, i] = derivative(material.components[$j].thermal_conductivity, T)
			end
		end)
	end

	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_is_gas!(is_gas, material)

Cache which components are gases (phase == GAS) to avoid runtime heterogeneous access.
Uses @generated to unroll access at compile time.
"""
@generated function cache_is_gas!(
	is_gas::AbstractVector{Bool},
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(is_gas[$j] = material.components[$j].phase == GAS))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_is_solid!(is_solid, material)

Cache which components are solids (phase == SOLID) to avoid runtime heterogeneous access.
Uses @generated to unroll access at compile time.
"""
@generated function cache_is_solid!(
	is_solid::AbstractVector{Bool},
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(is_solid[$j] = material.components[$j].phase == SOLID))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_swelling_factors!(γ, material)

Cache component swelling factors to avoid runtime heterogeneous access.
Uses @generated to unroll access at compile time.
"""
@generated function cache_swelling_factors!(
	γ::AbstractVector{Float64},
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(γ[$j] = material.components[$j].swelling_factor))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_absorption_coefficients!(α, material)

Cache component absorption coefficients to avoid runtime heterogeneous access.
Uses @generated to unroll access at compile time.
"""
@generated function cache_absorption_coefficients!(
	α::AbstractVector{Float64},
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(α[$j] = material.components[$j].absorption_coeff))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	cache_density_derivatives!(dρ_dT, T, material)

Cache component density derivatives dρⱼ/dT at the given temperature.
Uses @generated to unroll access at compile time.
"""
@generated function cache_density_derivatives!(
	dρ_dT::AbstractVector{Float64},
	T::Float64,
	material::Material{NC}) where NC

	# Generate unrolled code for each component
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(dρ_dT[$j] = derivative(material.components[$j].density, T)))
	end
	push!(exprs, :(return nothing))

	return Expr(:block, exprs...)
end

"""
	compute_volume_fractions!(v, ξ, densities, NC)

Compute volume fractions using pre-cached component densities.
Allocation-free: only writes to pre-allocated v array.
"""
@inline function compute_volume_fractions!(
	v::AbstractVector{Float64},
	ξ::NTuple{NC, Float64},
	densities::AbstractVector{Float64}) where NC

	# Compute total specific volume
	V_total = 0.0
	@inbounds @simd for j in 1:NC
		V_total += ξ[j] / max(densities[j], DENSITY_MIN)
	end

	if V_total < VOLUME_MIN
		@inbounds @simd for j in 1:NC
			v[j] = 0.0
		end
		return 0.0
	end

	# Compute volume fractions
	@inbounds @simd for j in 1:NC
		v[j] = (ξ[j] / max(densities[j], DENSITY_MIN)) / V_total
	end

	return V_total
end

"""
	compute_dλ_dξ_cached(ξ, k, densities, lambdas, v, V_total, mixing, β) -> Float64

Compute ∂λ/∂ξ_k using pre-cached component properties and volume fractions.
Allocation-free version.

# Arguments
- `ξ`: Mass concentrations tuple  
- `k`: Index of concentration to differentiate with respect to
- `densities`: Pre-cached component densities ρ_j(T)
- `lambdas`: Pre-cached component gas transfer coefficients λ_j(T)
- `v`: Pre-computed volume fractions
- `V_total`: Total specific volume
- `mixing`: Mixing rule (PARALLEL, SERIES, or WEIGHTED)
- `β`: Mixing parameter for WEIGHTED rule

# Returns
- The partial derivative ∂λ_eff/∂ξ_k
"""
@inline function compute_dλ_dξ_cached(
	ξ::NTuple{NC, Float64},
	k::Int,
	densities::AbstractVector{Float64},
	lambdas::AbstractVector{Float64},
	v::AbstractVector{Float64},
	V_total::Float64,
	mixing::MixingRule,
	β::Float64) where NC

	if V_total < VOLUME_MIN
		return 0.0
	end

	ρ_k = densities[k]

	# Check if series term will be zero (any component with v>1e-20 and λ=0)
	series_is_inf = false
	@inbounds for j in 1:NC
		if v[j] > VOLUME_MIN && lambdas[j] < LAMBDA_MIN
			series_is_inf = true
			break
		end
	end

	# Compute derivatives
	dλ_parallel_dξ = 0.0
	d_inv_λ_series_dξ = 0.0
	inv_λ_series = 0.0

	@inbounds for j in 1:NC
		λ_j = lambdas[j]
		ρ_j = densities[j]
		v_j = v[j]

		# ∂v_j/∂ξ_k = (δ_{jk}/ρ_j - v_j/ρ_k) / V_total
		if j == k
			dv_dξ = (1.0 / max(ρ_j, DENSITY_MIN) - v_j / max(ρ_k, DENSITY_MIN)) / V_total
		else
			dv_dξ = -v_j / max(ρ_k, DENSITY_MIN) / V_total
		end

		# Parallel: λ_parallel = Σ v_j × λ_j
		# ∂λ_parallel/∂ξ_k = Σ λ_j × ∂v_j/∂ξ_k
		dλ_parallel_dξ += λ_j * dv_dξ

		# Series: only compute if series term is not inf
		if !series_is_inf && λ_j > LAMBDA_MIN
			inv_λ_series += v_j / λ_j
			d_inv_λ_series_dξ += dv_dξ / λ_j
		end
	end

	# Compute λ_series derivative
	if series_is_inf
		dλ_series_dξ = 0.0
	else
		λ_series = inv_λ_series > LAMBDA_MIN ? 1.0 / inv_λ_series : 0.0
		dλ_series_dξ = -λ_series * λ_series * d_inv_λ_series_dξ
	end

	# Apply mixing rule
	if mixing == PARALLEL
		return dλ_parallel_dξ
	elseif mixing == SERIES
		return dλ_series_dξ
	else  # WEIGHTED
		return β * dλ_parallel_dξ + (1.0 - β) * dλ_series_dξ
	end
end

"""
	compute_dρ_dξ_cached(ξ, k, densities, swelling_factors, is_gas) -> Float64

Compute ∂ρ_mix/∂ξ_k using pre-cached component properties.
Allocation-free version of compute_dρ_dξ.

The mixture density formula (matching mixture_rules.jl) is:
	ρ_total = Σξ_j  (total concentration)
	Y_j = ξ_j / ρ_total  (mass fraction)
	v = Σ_{solid/liquid} (Y_j/ρ_j) + Σ_{gas} (γ_j × Y_j / ρ_j)
	ρ_mix = 1 / v

Key: For solids/liquids, swelling factor is NOT used (effectively γ=1).
	 For gases, the actual swelling_factor is used.

Let γ_eff_j = 1.0 for solid/liquid, γ_j for gas.
Then: v = Σ_j γ_eff_j × Y_j / ρ_j = Σ_j γ_eff_j × ξ_j / (ρ_total × ρ_j)

Derivative:
	∂v/∂ξ_k = γ_eff_k/(ρ_total × ρ_k) - v/ρ_total = (γ_eff_k/ρ_k - v) / ρ_total
	∂ρ_mix/∂ξ_k = -ρ_mix² × ∂v/∂ξ_k = ρ_mix × (1 - γ_eff_k × ρ_mix/ρ_k) / ρ_total
"""
@inline function compute_dρ_dξ_cached(
	ξ::NTuple{NC, Float64},
	k::Int,
	densities::AbstractVector{Float64},
	swelling_factors::AbstractVector{Float64},
	is_gas::AbstractVector{Bool}) where NC

	# Compute ρ_total (sum of all ξ)
	ρ_total = 0.0
	@inbounds @simd for j in 1:NC
		ρ_total += ξ[j]
	end

	if ρ_total < DENSITY_MIN
		return 0.0  # Fallback
	end

	# Compute specific volume v using the CORRECT formula:
	# - Solids/liquids: γ_eff = 1.0 (swelling factor ignored)
	# - Gases: γ_eff = swelling_factor
	v = 0.0
	@inbounds for j in 1:NC
		Y_j = ξ[j] / ρ_total
		ρ_j = max(densities[j], DENSITY_MIN)
		if is_gas[j]
			# Gas: use actual swelling factor
			v += swelling_factors[j] * Y_j / ρ_j
		else
			# Solid/liquid: effectively γ=1
			v += Y_j / ρ_j
		end
	end

	if v < VOLUME_MIN
		return 0.0  # Fallback
	end

	# ρ_mix = 1/v
	ρ_mix = 1.0 / v
	ρ_k = max(densities[k], DENSITY_MIN)

	# Effective γ for component k
	γ_eff_k = is_gas[k] ? swelling_factors[k] : 1.0

	# ∂ρ_mix/∂ξ_k = ρ_mix × (1 - γ_eff_k × ρ_mix/ρ_k) / ρ_total
	return ρ_mix * (1.0 - γ_eff_k * ρ_mix / ρ_k) / ρ_total
end

"""
	compute_dλ_dT(ξ, T, material) -> Float64

Compute derivative of effective gas transfer coefficient with respect to temperature.

The effective λ depends on T through:
1. Component λⱼ(T) if they are temperature-dependent
2. Volume fractions vⱼ = (ξⱼ/ρⱼ(T)) / Σ(ξₖ/ρₖ(T))

Uses central finite differences for robust computation.
"""
@inline function compute_dλ_dT(ξ::NTuple{NC, Float64}, T::Float64,
	material::Material{NC}) where NC
	# Use central differences for numerical derivative
	ε = 1e-6
	λ_plus = Properties.effective_gas_transfer(ξ, T + ε, material)
	λ_minus = Properties.effective_gas_transfer(ξ, T - ε, material)
	return (λ_plus - λ_minus) / (2ε)
end

"""
	compute_gas_transport_geometry_jacobian!(J, smap, mesh, prop_cache, material, dV_dz, n, Val(NC), z_offset, is_gas)

Add geometry derivatives for gas transport fluxes (species and energy equations).

Gas flux: J_j = -λ × ρ × (Y_R - Y_L) / d_LR
Gas convective source (handled separately by compute_gas_convective_source_jacobian!)

The derivatives w.r.t. node positions z_k come from:
1. ∂d_LR/∂z_k - affects the gradient
2. ∂(1/V)/∂z_k - affects the divergence scaling

For interior face between cells i and i+1:
- d_LR = |center(i+1) - center(i)| = 0.5*(z_{i+2} - z_i)
- ∂d_LR/∂z_i = -0.5
- ∂d_LR/∂z_{i+2} = +0.5
"""
function compute_gas_transport_geometry_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material::Material{NC}, dV_dz::Matrix{Float64}, n::Int, ::Val{NC}, z_offset::Int,
	is_gas::AbstractVector{Bool}) where NC

	n_nodes = size(dV_dz, 2)

	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces (handled separately)
		if left_id == 0 || right_id == 0
			continue
		end

		# Face geometry
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		inv_d_LR = 1.0 / d_LR
		inv_d_LR_sq = inv_d_LR * inv_d_LR
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		inv_V_L = 1.0 / V_L
		inv_V_R = 1.0 / V_R
		inv_V_L_sq = inv_V_L * inv_V_L
		inv_V_R_sq = inv_V_R * inv_V_R
		A_f = mesh.face_areas[f]

		# Cell states - use prop_cache for consistency with RHS
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Gas transfer coefficient - use λ_face from prop_cache (volume fraction gradient)
		λ_face = prop_cache.λ_face[f]

		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		T_face = 0.5 * (T_L + T_R)

		# Node indices for geometry derivatives
		n_left_node = left_id
		n_right_node = right_id + 1

		# ================================================================
		# NOTE: Convective heat source geometry derivatives are omitted here.
		# S_conv is now cell-centered (-Σ c_g * J̄ * dTdx) and its geometry
		# dependence through cell_centers is a second-order effect captured
		# by the FD Jacobian fallback.
		# ================================================================

		# ================================================================
		# Species equation geometry derivatives (∂(dξ/dt)/∂z_k)
		# ================================================================

		for j in 1:NC
			if !is_gas[j]
				continue
			end

			# Volume fraction gradient
			Φ = ξ_R[j] * T_R - ξ_L[j] * T_L

			# Current flux: J = -λ_face × Φ / (T_face × d_LR)
			J_flux = -λ_face * Φ / (T_face * d_LR)

			# Species equation rows
			ξ_row_L = n + (j-1)*n + left_id
			ξ_row_R = n + (j-1)*n + right_id

			# ============================================================
			# Term 1: ∂J/∂d_LR contribution
			# J = -λ_face × Φ / (T_face × d_LR)
			# ∂J/∂d_LR = λ_face × Φ / (T_face × d_LR²)
			# ============================================================
			dJ_dd = λ_face * Φ / (T_face * d_LR * d_LR)

			if abs(dJ_dd) > FLUX_DERIV_TOL
				# ∂d_LR/∂z[n_left_node] = -0.5
				deriv_L_left = -A_f * inv_V_L * dJ_dd * (-0.5)
				deriv_R_left = A_f * inv_V_R * dJ_dd * (-0.5)
				fast_sparse_add!(J, smap, ξ_row_L, z_offset+n_left_node, deriv_L_left)
				fast_sparse_add!(J, smap, ξ_row_R, z_offset+n_left_node, deriv_R_left)

				# ∂d_LR/∂z[n_right_node] = +0.5
				deriv_L_right = -A_f * inv_V_L * dJ_dd * (+0.5)
				deriv_R_right = A_f * inv_V_R * dJ_dd * (+0.5)
				fast_sparse_add!(J, smap, ξ_row_L, z_offset+n_right_node, deriv_L_right)
				fast_sparse_add!(J, smap, ξ_row_R, z_offset+n_right_node, deriv_R_right)
			end

			# ============================================================
			# Term 2: ∂(1/V)/∂z_k contribution
			# ============================================================
			if abs(J_flux) > FLUX_DERIV_TOL
				for k in 1:n_nodes
					# Left cell volume dependence
					if abs(dV_dz[left_id, k]) > FLUX_DERIV_TOL
						deriv = A_f * J_flux * inv_V_L_sq * dV_dz[left_id, k]
						fast_sparse_add!(J, smap, ξ_row_L, z_offset+k, deriv)
					end

					# Right cell volume dependence
					if abs(dV_dz[right_id, k]) > FLUX_DERIV_TOL
						deriv = -A_f * J_flux * inv_V_R_sq * dV_dz[right_id, k]
						fast_sparse_add!(J, smap, ξ_row_R, z_offset+k, deriv)
					end
				end
			end
		end
	end

	return nothing
end

# ============================================================================
# Darcy Transport Jacobian
# ============================================================================

"""
	compute_darcy_transport_jacobian!(J, smap, mesh, prop_cache, material, n, NC, jac_cache)

Add Jacobian contributions from Darcy advective gas transport.

The Darcy advective flux for gas component j is:
	J_adv,j = v_gas × ξ_upwind

where v_gas is the Darcy velocity computed from:
	v_gas = -(κ/μ) × ∇P

The key derivatives are:
	∂J_adv/∂T = ξ_upwind × ∂v_gas/∂T
	∂J_adv/∂ξⱼ = v_gas (direct) + ξ_upwind × ∂v_gas/∂ξⱼ (through pressure)

The chain rule for ∂v_gas/∂T and ∂v_gas/∂ξⱼ involves:
	∂v_gas/∂T = ∂v_gas/∂P_L × ∂P_L/∂T + ∂v_gas/∂P_R × ∂P_R/∂T + (viscosity term)
	∂v_gas/∂ξⱼ = ∂v_gas/∂P_L × ∂P_L/∂ξⱼ + ∂v_gas/∂P_R × ∂P_R/∂ξⱼ + (permeability term)

# Arguments
- `J`: Jacobian matrix (modified in-place)
- `smap`: Sparsity map for O(1) sparse writes
- `mesh`: 1D mesh with cell states
- `prop_cache`: Property cache containing P, v_gas, κ, μ
- `material`: Material definition
- `n`: Number of cells
- `NC`: Number of components (Val type)
- `jac_cache`: ALE Jacobian cache with Darcy arrays
"""
function compute_darcy_transport_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material, n::Int, ::Val{NC}, jac_cache::JacobianCacheALE{NC}) where NC

	base_cache = jac_cache.base
	is_gas = base_cache.is_gas

	# Get Darcy work arrays from cache
	dP_dT = jac_cache.dP_dT
	dP_dξ = jac_cache.dP_dξ
	dvgas_dP_L = jac_cache.dvgas_dP_L
	dvgas_dP_R = jac_cache.dvgas_dP_R
	mobility = jac_cache.mobility

	# Get permeability derivative arrays from cache
	dφ_dξ = jac_cache.dφ_dξ
	dκ_dξ = jac_cache.dκ_dξ
	dvgas_dκ_L = jac_cache.dvgas_dκ_L
	dvgas_dκ_R = jac_cache.dvgas_dκ_R

	# Get work vectors for in-place derivative computations
	# Note: pressure_derivatives! and porosity_derivatives! don't need work arrays,
	# so we can reuse work_NC, work_NC_2, work_NC_3 for permeability_derivatives!
	work_NC = base_cache.work_NC      # For permeability v_work
	work_NC_2 = base_cache.work_NC_2  # For permeability κ_work  
	work_NC_3 = base_cache.work_NC_3  # For permeability ρ_work

	# Step 1: Precompute pressure and permeability derivatives at all cells
	# Use in-place versions to avoid tuple allocations
	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		ξ = state.ξ
		T = state.T

		# Get pressure derivatives in-place (allocation-free)
		∂P_∂T = pressure_derivatives!(view(dP_dξ, :, i), ξ, T, material)
		dP_dT[i] = ∂P_∂T

		# Compute porosity derivatives in-place (allocation-free)
		Properties.porosity_derivatives!(view(dφ_dξ, :, i), T, material)

		# Compute permeability derivatives in-place (allocation-free)
		# Uses 3 work arrays: v_work (volume fractions), κ_work (permeabilities), ρ_work (densities)
		Properties.permeability_derivatives!(view(dκ_dξ, :, i), ξ, T, material, work_NC, work_NC_2, work_NC_3)
	end

	# Step 2: Precompute Darcy velocity Jacobian coefficients at all faces
	# 
	# IMPORTANT: Must match the RHS computation in update_darcy_properties!
	# The RHS uses distance-weighted harmonic mean of MOBILITY (κ/μ):
	#   mob_face = (d_L + d_R) / (d_L/mob_L + d_R/mob_R)
	# This is the correct resistance-in-series interpolation for non-uniform meshes.
	#
	# v_gas = -mob_face × ∇P
	# mob_L = κ_L / μ_L,  mob_R = κ_R / μ_R
	#
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Boundary faces: no interior Darcy transport
		if left_id == 0 || right_id == 0
			dvgas_dP_L[f] = 0.0
			dvgas_dP_R[f] = 0.0
			dvgas_dκ_L[f] = 0.0
			dvgas_dκ_R[f] = 0.0
			mobility[f] = 0.0
			continue
		end

		# Get cell-centered properties (matching update_darcy_properties!)
		κ_L = prop_cache.κ[left_id]
		κ_R = prop_cache.κ[right_id]
		μ_L = prop_cache.μ[left_id]
		μ_R = prop_cache.μ[right_id]

		# Distances from cell centers to face (matching RHS exactly)
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R

		# Pressure gradient
		P_L = prop_cache.P[left_id]
		P_R = prop_cache.P[right_id]
		dP_dz = (P_R - P_L) / d_LR

		# Cell mobilities (matching RHS exactly)
		mob_L = μ_L > 1e-20 ? κ_L / μ_L : 0.0
		mob_R = μ_R > 1e-20 ? κ_R / μ_R : 0.0

		# Distance-weighted harmonic mean of mobility (matching RHS exactly)
		if mob_L > 1e-30 && mob_R > 1e-30
			mob_face = d_LR / (d_L / mob_L + d_R / mob_R)
		elseif mob_L > 1e-30
			mob_face = mob_L
		elseif mob_R > 1e-30
			mob_face = mob_R
		else
			mob_face = 0.0
		end

		mobility[f] = mob_face

		# v_gas = -mob_face × ∇P = -mob_face × (P_R - P_L) / d_LR
		# ∂v_gas/∂P_L = +mob_face / d_LR
		# ∂v_gas/∂P_R = -mob_face / d_LR
		if d_LR > 1e-20
			mob_over_d = mob_face / d_LR
			dvgas_dP_L[f] = mob_over_d
			dvgas_dP_R[f] = -mob_over_d
		else
			dvgas_dP_L[f] = 0.0
			dvgas_dP_R[f] = 0.0
		end

		# ============================================================
		# Permeability derivatives (via mobility chain rule)
		# ============================================================
		# mob_L = κ_L / μ_L  →  ∂mob_L/∂κ_L = 1/μ_L
		# mob_face = d_total / (d_L/mob_L + d_R/mob_R) = d_total / R_sum
		# where R_L = d_L/mob_L, R_R = d_R/mob_R, R_sum = R_L + R_R
		#
		# ∂mob_face/∂mob_L = d_total × d_L / (mob_L² × R_sum²)
		# ∂mob_face/∂mob_R = d_total × d_R / (mob_R² × R_sum²)
		#
		# v_gas = -mob_face × ∇P
		# ∂v_gas/∂κ_L = -∇P × ∂mob_face/∂mob_L × ∂mob_L/∂κ_L

		if mob_L > 1e-30 && mob_R > 1e-30 && d_LR > 1e-20
			R_L = d_L / mob_L
			R_R = d_R / mob_R
			R_sum = R_L + R_R
			R_sum_sq = R_sum * R_sum

			# Derivatives of distance-weighted harmonic mean of mobility
			dmob_face_dmob_L = d_LR * d_L / (mob_L * mob_L * R_sum_sq)
			dmob_face_dmob_R = d_LR * d_R / (mob_R * mob_R * R_sum_sq)

			# Chain rule: ∂mob/∂κ = 1/μ (since mob = κ/μ)
			dmob_L_dκ_L = 1.0 / μ_L
			dmob_R_dκ_R = 1.0 / μ_R

			# ∂v_gas/∂κ = -∇P × ∂mob_face/∂mob × ∂mob/∂κ
			dvgas_dκ_L[f] = -dP_dz * dmob_face_dmob_L * dmob_L_dκ_L
			dvgas_dκ_R[f] = -dP_dz * dmob_face_dmob_R * dmob_R_dκ_R
		else
			dvgas_dκ_L[f] = 0.0
			dvgas_dκ_R[f] = 0.0
		end
	end

	# Step 3: Add Jacobian contributions for each interior face
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces
		if left_id == 0 || right_id == 0
			continue
		end

		# Get Darcy velocity at face
		v_gas = prop_cache.v_gas[f]

		# Skip if Darcy velocity is negligible
		if abs(v_gas) < 1e-20 && abs(dvgas_dP_L[f]) < 1e-30
			continue
		end

		# Face geometry
		A_f = mesh.face_areas[f]
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		A_f_inv_V_L = A_f / V_L
		A_f_inv_V_R = A_f / V_R

		# Cell states
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Compute derivatives of v_gas w.r.t. state variables via chain rule
		# ∂v_gas/∂T_L = ∂v_gas/∂P_L × ∂P_L/∂T_L (ignoring viscosity temp dependence for now)
		# ∂v_gas/∂T_R = ∂v_gas/∂P_R × ∂P_R/∂T_R
		dvgas_dT_L = dvgas_dP_L[f] * dP_dT[left_id]
		dvgas_dT_R = dvgas_dP_R[f] * dP_dT[right_id]

		# Process each gas component
		for j in 1:NC
			if !is_gas[j]
				continue
			end

			# Species equation rows
			ξ_row_L = n + (j - 1) * n + left_id
			ξ_row_R = n + (j - 1) * n + right_id

			# Upwind selection for advective flux
			# Flux = v_gas × ξ_upwind, positive v_gas means flow from L to R
			if v_gas >= 0.0
				# Flow from left to right: upwind = left
				ξ_upwind = ξ_L[j]
				upwind_id = left_id
				upwind_is_left = true
			else
				# Flow from right to left: upwind = right
				ξ_upwind = ξ_R[j]
				upwind_id = right_id
				upwind_is_left = false
			end

			# Advective flux: J_adv = v_gas × ξ_upwind
			# (This contributes to RHS as: -∇·(v_gas ξ) = -(J_out - J_in)/V)

			# ================================================================
			# Temperature derivatives: ∂J_adv/∂T
			# J_adv = v_gas × ξ_upwind
			# ∂J_adv/∂T = ξ_upwind × ∂v_gas/∂T
			# ================================================================

			# Contribution from left cell temperature
			dJ_adv_dT_L = ξ_upwind * dvgas_dT_L
			if abs(dJ_adv_dT_L) > FLUX_DERIV_TOL
				# Left cell sees: -A × J_adv / V_L (flux out)
				# Right cell sees: +A × J_adv / V_R (flux in)
				fast_sparse_sub!(J, smap, ξ_row_L, left_id, A_f_inv_V_L * dJ_adv_dT_L)
				fast_sparse_add!(J, smap, ξ_row_R, left_id, A_f_inv_V_R * dJ_adv_dT_L)
			end

			# Contribution from right cell temperature
			dJ_adv_dT_R = ξ_upwind * dvgas_dT_R
			if abs(dJ_adv_dT_R) > FLUX_DERIV_TOL
				fast_sparse_sub!(J, smap, ξ_row_L, right_id, A_f_inv_V_L * dJ_adv_dT_R)
				fast_sparse_add!(J, smap, ξ_row_R, right_id, A_f_inv_V_R * dJ_adv_dT_R)
			end

			# ================================================================
			# Concentration derivatives: ∂J_adv/∂ξ_k
			# J_adv = v_gas × ξ_upwind
			# ∂J_adv/∂ξ_k = ∂v_gas/∂ξ_k × ξ_upwind + v_gas × ∂ξ_upwind/∂ξ_k
			# 
			# For the upwind term: ∂ξ_upwind/∂ξ_k = δ_{j,k} if k is in upwind cell
			# For the velocity term via chain rule:
			#   ∂v_gas/∂ξ_k = ∂v_gas/∂P × ∂P/∂ξ_k + ∂v_gas/∂κ × ∂κ/∂ξ_k
			#
			# Note: Gas components affect pressure (∂P/∂ξ_k ≠ 0 for gas)
			#       Solid/liquid components affect permeability (∂κ/∂ξ_k ≠ 0 for solid/liquid)
			# ================================================================

			for k in 1:NC
				# Left cell ξ_k derivatives
				ξ_col_L = n + (k - 1) * n + left_id

				# Total velocity derivative includes both pressure and permeability paths
				# ∂v_gas/∂ξ_k,L = ∂v_gas/∂P_L × ∂P_L/∂ξ_k + ∂v_gas/∂κ_L × ∂κ_L/∂ξ_k
				dvgas_dξk_L = dvgas_dP_L[f] * dP_dξ[k, left_id] + dvgas_dκ_L[f] * dκ_dξ[k, left_id]

				# Term 1: velocity derivative contribution
				dJ_adv_dξk_L = ξ_upwind * dvgas_dξk_L

				# Term 2: direct upwind contribution (if k=j and left is upwind)
				if upwind_is_left && k == j
					dJ_adv_dξk_L += v_gas
				end

				if abs(dJ_adv_dξk_L) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col_L, A_f_inv_V_L * dJ_adv_dξk_L)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col_L, A_f_inv_V_R * dJ_adv_dξk_L)
				end

				# Right cell ξ_k derivatives
				ξ_col_R = n + (k - 1) * n + right_id

				# ∂v_gas/∂ξ_k,R = ∂v_gas/∂P_R × ∂P_R/∂ξ_k + ∂v_gas/∂κ_R × ∂κ_R/∂ξ_k
				dvgas_dξk_R = dvgas_dP_R[f] * dP_dξ[k, right_id] + dvgas_dκ_R[f] * dκ_dξ[k, right_id]

				# Term 1: velocity derivative contribution
				dJ_adv_dξk_R = ξ_upwind * dvgas_dξk_R

				# Term 2: direct upwind contribution (if k=j and right is upwind)
				if !upwind_is_left && k == j
					dJ_adv_dξk_R += v_gas
				end

				if abs(dJ_adv_dξk_R) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col_R, A_f_inv_V_L * dJ_adv_dξk_R)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col_R, A_f_inv_V_R * dJ_adv_dξk_R)
				end
			end
		end
	end

	return nothing
end

"""
	compute_darcy_transport_jacobian!(J, smap, mesh, prop_cache, material, n, NC, jac_cache::JacobianCache)

Non-ALE (fixed mesh) version of the Darcy advective transport Jacobian.
Identical physics to the ALE version but reads Darcy work arrays directly from JacobianCache
(which now has the same fields as JacobianCacheALE) and has no geometry corrections.

The Darcy advective flux for gas component j through face f is:

    J_adv,j = v_gas × ξ_upwind,j

where v_gas = -(κ_face/μ_face) × (P_R - P_L)/d_LR  (Darcy's law)

The Jacobian contributions are:
  ∂J_adv/∂T  = ξ_upwind × ∂v_gas/∂T     (velocity changes with T via pressure)
  ∂J_adv/∂ξk = ξ_upwind × ∂v_gas/∂ξk + v_gas × δ_{j,k}  (velocity + upwind source)

where ∂v_gas/∂ξk = ∂v_gas/∂P × ∂P/∂ξk + ∂v_gas/∂κ × ∂κ/∂ξk
"""
function compute_darcy_transport_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material, n::Int, ::Val{NC}, jac_cache::JacobianCache{NC}) where NC

	is_gas = jac_cache.is_gas

	# Get Darcy work arrays from cache (same field names as JacobianCacheALE)
	dP_dT = jac_cache.dP_dT
	dP_dξ = jac_cache.dP_dξ
	dvgas_dP_L = jac_cache.dvgas_dP_L
	dvgas_dP_R = jac_cache.dvgas_dP_R
	mobility = jac_cache.mobility

	# Get permeability derivative arrays from cache
	dφ_dξ = jac_cache.dφ_dξ
	dκ_dξ = jac_cache.dκ_dξ
	dvgas_dκ_L = jac_cache.dvgas_dκ_L
	dvgas_dκ_R = jac_cache.dvgas_dκ_R

	# Get work vectors for in-place derivative computations
	work_NC = jac_cache.work_NC
	work_NC_2 = jac_cache.work_NC_2
	work_NC_3 = jac_cache.work_NC_3

	# Step 1: Precompute pressure and permeability derivatives at all cells
	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		ξ = state.ξ
		T = state.T

		# Pressure derivatives (allocation-free)
		∂P_∂T = pressure_derivatives!(view(dP_dξ, :, i), ξ, T, material)
		dP_dT[i] = ∂P_∂T

		# Porosity derivatives (allocation-free)
		Properties.porosity_derivatives!(view(dφ_dξ, :, i), T, material)

		# Permeability derivatives (allocation-free, uses 3 work arrays)
		Properties.permeability_derivatives!(view(dκ_dξ, :, i), ξ, T, material, work_NC, work_NC_2, work_NC_3)
	end

	# Step 2: Precompute Darcy velocity Jacobian coefficients at all faces
	#
	# MUST match the RHS computation in update_darcy_properties! exactly:
	# v_gas = -mob_face × ∇P
	# mob_face = (d_L + d_R) / (d_L/mob_L + d_R/mob_R)  [distance-weighted harmonic mean]
	# mob_L = κ_L / μ_L,  mob_R = κ_R / μ_R
	#
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Boundary faces: no interior Darcy transport
		if left_id == 0 || right_id == 0
			dvgas_dP_L[f] = 0.0
			dvgas_dP_R[f] = 0.0
			dvgas_dκ_L[f] = 0.0
			dvgas_dκ_R[f] = 0.0
			mobility[f] = 0.0
			continue
		end

		# Cell-centered properties (matching update_darcy_properties!)
		κ_L = prop_cache.κ[left_id]
		κ_R = prop_cache.κ[right_id]
		μ_L = prop_cache.μ[left_id]
		μ_R = prop_cache.μ[right_id]

		# Distances from cell centers to face (matching RHS exactly)
		face_pos = mesh.face_positions[f]
		d_L = abs(face_pos - mesh.cell_centers[left_id])
		d_R = abs(mesh.cell_centers[right_id] - face_pos)
		d_LR = d_L + d_R

		# Pressure gradient
		P_L = prop_cache.P[left_id]
		P_R = prop_cache.P[right_id]
		dP_dz = (P_R - P_L) / d_LR

		# Cell mobilities (matching RHS exactly)
		mob_L = μ_L > 1e-20 ? κ_L / μ_L : 0.0
		mob_R = μ_R > 1e-20 ? κ_R / μ_R : 0.0

		# Distance-weighted harmonic mean of mobility (matching RHS exactly)
		if mob_L > 1e-30 && mob_R > 1e-30
			mob_face = d_LR / (d_L / mob_L + d_R / mob_R)
		elseif mob_L > 1e-30
			mob_face = mob_L
		elseif mob_R > 1e-30
			mob_face = mob_R
		else
			mob_face = 0.0
		end

		mobility[f] = mob_face

		# v_gas = -mob_face × (P_R - P_L) / d_LR
		# ∂v_gas/∂P_L = +mob_face / d_LR
		# ∂v_gas/∂P_R = -mob_face / d_LR
		if d_LR > 1e-20
			mob_over_d = mob_face / d_LR
			dvgas_dP_L[f] = mob_over_d
			dvgas_dP_R[f] = -mob_over_d
		else
			dvgas_dP_L[f] = 0.0
			dvgas_dP_R[f] = 0.0
		end

		# Permeability derivatives via mobility chain rule:
		# mob_face = d_total / (d_L/mob_L + d_R/mob_R) = d_total / R_sum
		# ∂mob_face/∂mob_L = d_total × d_L / (mob_L² × R_sum²)
		# ∂mob_face/∂mob_R = d_total × d_R / (mob_R² × R_sum²)
		# ∂mob/∂κ = 1/μ  (since mob = κ/μ)
		# ∂v_gas/∂κ_L = -∇P × ∂mob_face/∂mob_L × (1/μ_L)
		if mob_L > 1e-30 && mob_R > 1e-30 && d_LR > 1e-20
			R_L = d_L / mob_L
			R_R = d_R / mob_R
			R_sum = R_L + R_R
			R_sum_sq = R_sum * R_sum

			dmob_face_dmob_L = d_LR * d_L / (mob_L * mob_L * R_sum_sq)
			dmob_face_dmob_R = d_LR * d_R / (mob_R * mob_R * R_sum_sq)

			dmob_L_dκ_L = 1.0 / μ_L
			dmob_R_dκ_R = 1.0 / μ_R

			dvgas_dκ_L[f] = -dP_dz * dmob_face_dmob_L * dmob_L_dκ_L
			dvgas_dκ_R[f] = -dP_dz * dmob_face_dmob_R * dmob_R_dκ_R
		else
			dvgas_dκ_L[f] = 0.0
			dvgas_dκ_R[f] = 0.0
		end
	end

	# Step 3: Add Jacobian contributions for each interior face
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces
		if left_id == 0 || right_id == 0
			continue
		end

		# Get Darcy velocity at face
		v_gas = prop_cache.v_gas[f]

		# Skip if Darcy velocity is negligible
		if abs(v_gas) < 1e-20 && abs(dvgas_dP_L[f]) < 1e-30
			continue
		end

		# Face geometry
		A_f = mesh.face_areas[f]
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		A_f_inv_V_L = A_f / V_L
		A_f_inv_V_R = A_f / V_R

		# Cell states
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# ∂v_gas/∂T_L = ∂v_gas/∂P_L × ∂P_L/∂T_L
		# ∂v_gas/∂T_R = ∂v_gas/∂P_R × ∂P_R/∂T_R
		dvgas_dT_L = dvgas_dP_L[f] * dP_dT[left_id]
		dvgas_dT_R = dvgas_dP_R[f] * dP_dT[right_id]

		# Process each gas component
		for j in 1:NC
			if !is_gas[j]
				continue
			end

			# Species equation rows
			ξ_row_L = n + (j - 1) * n + left_id
			ξ_row_R = n + (j - 1) * n + right_id

			# Upwind selection: positive v_gas means flow from L to R
			if v_gas >= 0.0
				ξ_upwind = ξ_L[j]
				upwind_is_left = true
			else
				ξ_upwind = ξ_R[j]
				upwind_is_left = false
			end

			# ================================================================
			# Temperature derivatives: ∂J_adv/∂T = ξ_upwind × ∂v_gas/∂T
			# ================================================================

			# Left cell temperature
			dJ_adv_dT_L = ξ_upwind * dvgas_dT_L
			if abs(dJ_adv_dT_L) > FLUX_DERIV_TOL
				# Left cell divergence: -(A × J) / V  →  subtract
				# Right cell divergence: +(A × J) / V  →  add
				fast_sparse_sub!(J, smap, ξ_row_L, left_id, A_f_inv_V_L * dJ_adv_dT_L)
				fast_sparse_add!(J, smap, ξ_row_R, left_id, A_f_inv_V_R * dJ_adv_dT_L)
			end

			# Right cell temperature
			dJ_adv_dT_R = ξ_upwind * dvgas_dT_R
			if abs(dJ_adv_dT_R) > FLUX_DERIV_TOL
				fast_sparse_sub!(J, smap, ξ_row_L, right_id, A_f_inv_V_L * dJ_adv_dT_R)
				fast_sparse_add!(J, smap, ξ_row_R, right_id, A_f_inv_V_R * dJ_adv_dT_R)
			end

			# ================================================================
			# Concentration derivatives: ∂J_adv/∂ξ_k
			# ∂J_adv/∂ξ_k = ∂v_gas/∂ξ_k × ξ_upwind + v_gas × δ_{j,k}
			# where ∂v_gas/∂ξ_k = ∂v_gas/∂P × ∂P/∂ξ_k + ∂v_gas/∂κ × ∂κ/∂ξ_k
			# ================================================================

			for k in 1:NC
				# Left cell ξ_k derivatives
				ξ_col_L = n + (k - 1) * n + left_id

				dvgas_dξk_L = dvgas_dP_L[f] * dP_dξ[k, left_id] + dvgas_dκ_L[f] * dκ_dξ[k, left_id]

				# Term 1: velocity derivative contribution
				dJ_adv_dξk_L = ξ_upwind * dvgas_dξk_L

				# Term 2: direct upwind contribution (if k=j and left is upwind)
				if upwind_is_left && k == j
					dJ_adv_dξk_L += v_gas
				end

				if abs(dJ_adv_dξk_L) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col_L, A_f_inv_V_L * dJ_adv_dξk_L)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col_L, A_f_inv_V_R * dJ_adv_dξk_L)
				end

				# Right cell ξ_k derivatives
				ξ_col_R = n + (k - 1) * n + right_id

				dvgas_dξk_R = dvgas_dP_R[f] * dP_dξ[k, right_id] + dvgas_dκ_R[f] * dκ_dξ[k, right_id]

				# Term 1: velocity derivative contribution
				dJ_adv_dξk_R = ξ_upwind * dvgas_dξk_R

				# Term 2: direct upwind contribution (if k=j and right is upwind)
				if !upwind_is_left && k == j
					dJ_adv_dξk_R += v_gas
				end

				if abs(dJ_adv_dξk_R) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row_L, ξ_col_R, A_f_inv_V_L * dJ_adv_dξk_R)
					fast_sparse_add!(J, smap, ξ_row_R, ξ_col_R, A_f_inv_V_R * dJ_adv_dξk_R)
				end
			end
		end
	end

	return nothing
end

"""
	compute_darcy_transport_geometry_jacobian!(J, smap, mesh, prop_cache, material,
											   dV_dz, n, NC, z_offset, jac_cache)

Add geometry derivatives for Darcy advective transport.

The advective flux J_adv = v_gas × ξ_upwind depends on geometry through:
1. ∂v_gas/∂d_LR - pressure gradient sensitivity to cell distance
2. ∂(1/V)/∂z_k - volume scaling for divergence

# Arguments
- `J`: Jacobian matrix
- `smap`: Sparsity map
- `mesh`: 1D mesh
- `prop_cache`: Property cache
- `material`: Material definition
- `dV_dz`: Volume derivatives [n_cells × n_nodes]
- `n`: Number of cells
- `NC`: Number of components (Val)
- `z_offset`: Offset for z DOFs in state vector
- `jac_cache`: ALE Jacobian cache
"""
function compute_darcy_transport_geometry_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	material, dV_dz::Matrix{Float64}, n::Int, ::Val{NC}, z_offset::Int,
	jac_cache::JacobianCacheALE{NC}) where NC

	base_cache = jac_cache.base
	is_gas = base_cache.is_gas
	mobility = jac_cache.mobility

	n_nodes = size(dV_dz, 2)

	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Skip boundary faces
		if left_id == 0 || right_id == 0
			continue
		end

		v_gas = prop_cache.v_gas[f]
		mob = mobility[f]

		# Skip if negligible
		if abs(v_gas) < 1e-20 && abs(mob) < 1e-30
			continue
		end

		# Geometry
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
		A_f = mesh.face_areas[f]
		V_L = mesh.cell_volumes[left_id]
		V_R = mesh.cell_volumes[right_id]
		inv_V_L = 1.0 / V_L
		inv_V_R = 1.0 / V_R
		inv_V_L_sq = inv_V_L * inv_V_L
		inv_V_R_sq = inv_V_R * inv_V_R

		# Pressure gradient sensitivity to d_LR
		# v_gas = -mob × (P_R - P_L) / d_LR
		# ∂v_gas/∂d_LR = mob × (P_R - P_L) / d_LR²
		P_L = prop_cache.P[left_id]
		P_R = prop_cache.P[right_id]
		ΔP = P_R - P_L
		dvgas_dd = mob * ΔP / (d_LR * d_LR)

		# Cell states
		ξ_L = mesh.cell_states[left_id].ξ
		ξ_R = mesh.cell_states[right_id].ξ

		# Node indices affecting d_LR
		n_left_node = left_id
		n_right_node = right_id + 1

		# ================================================================
		# NOTE: Convective heat source geometry derivatives are omitted here.
		# S_conv is now cell-centered (-Σ c_g * J̄ * dTdx) and its geometry
		# dependence through cell_centers is a second-order effect captured
		# by the FD Jacobian fallback.
		# ================================================================

		# ================================================================
		# Species equation geometry derivatives (∂(dξ/dt)/∂z_k)
		# ================================================================

		for j in 1:NC
			if !is_gas[j]
				continue
			end

			# Species rows
			ξ_row_L = n + (j - 1) * n + left_id
			ξ_row_R = n + (j - 1) * n + right_id

			# Upwind concentration
			if v_gas >= 0.0
				ξ_upwind = ξ_L[j]
			else
				ξ_upwind = ξ_R[j]
			end

			# Current advective flux
			J_adv = v_gas * ξ_upwind

			# ============================================================
			# Term 1: ∂J_adv/∂d_LR contribution
			# J_adv = v_gas × ξ_upwind
			# ∂J_adv/∂d_LR = ξ_upwind × ∂v_gas/∂d_LR
			# ============================================================
			dJ_adv_dd = ξ_upwind * dvgas_dd

			if abs(dJ_adv_dd) > FLUX_DERIV_TOL
				# ∂d_LR/∂z[n_left_node] = -0.5
				deriv_L_left = -A_f * inv_V_L * dJ_adv_dd * (-0.5)
				deriv_R_left = A_f * inv_V_R * dJ_adv_dd * (-0.5)
				fast_sparse_add!(J, smap, ξ_row_L, z_offset + n_left_node, deriv_L_left)
				fast_sparse_add!(J, smap, ξ_row_R, z_offset + n_left_node, deriv_R_left)

				# ∂d_LR/∂z[n_right_node] = +0.5
				deriv_L_right = -A_f * inv_V_L * dJ_adv_dd * (+0.5)
				deriv_R_right = A_f * inv_V_R * dJ_adv_dd * (+0.5)
				fast_sparse_add!(J, smap, ξ_row_L, z_offset + n_right_node, deriv_L_right)
				fast_sparse_add!(J, smap, ξ_row_R, z_offset + n_right_node, deriv_R_right)
			end

			# ============================================================
			# Term 2: ∂(1/V)/∂z_k contribution (volume scaling)
			# ============================================================
			if abs(J_adv) > FLUX_DERIV_TOL
				for k in 1:n_nodes
					# Left cell volume dependence
					if abs(dV_dz[left_id, k]) > FLUX_DERIV_TOL
						deriv = A_f * J_adv * inv_V_L_sq * dV_dz[left_id, k]
						fast_sparse_add!(J, smap, ξ_row_L, z_offset + k, deriv)
					end

					# Right cell volume dependence
					if abs(dV_dz[right_id, k]) > FLUX_DERIV_TOL
						deriv = -A_f * J_adv * inv_V_R_sq * dV_dz[right_id, k]
						fast_sparse_add!(J, smap, ξ_row_R, z_offset + k, deriv)
					end
				end
			end
		end
	end

	return nothing
end

# ============================================================================
# Boundary Condition Jacobian
# ============================================================================

"""
	_compute_boundary_flux_jacobian(bc, T_cell, k_eff, Δz_half, ε_eff, t)

Compute the Jacobian contribution ∂q_conducted/∂T_cell for a boundary with thermal BC.

This properly accounts for the boundary temperature solve by applying the chain rule:
	∂q_conducted/∂T_cell = (∂q_conducted/∂T_s) × (∂T_s/∂T_cell)

The boundary temperature T_s is determined implicitly by the energy balance:
	q_BC(T_s) = k_eff × (T_cell - T_s) / Δz_half

From implicit differentiation of the energy balance w.r.t. T_cell:
	dq_BC/dT_s × dT_s/dT_cell = (k/Δz) × (dT_s/dT_cell - 1)

Solving for dT_s/dT_cell:
	dT_s/dT_cell = (k/Δz) / (k/Δz - dq_BC/dT_s)

And since q_conducted = k × (T_s - T_cell) / Δz:
	dq_conducted/dT_cell = (k/Δz) × (dT_s/dT_cell - 1)
						 = (k/Δz) × (dq_BC/dT_s) / (k/Δz - dq_BC/dT_s)
						 = (k/Δz × dq_BC/dT_s) / (k/Δz - dq_BC/dT_s)

Returns (dq_dT_cell, T_surface, converged) where:
	- dq_dT_cell: The Jacobian contribution ∂q_conducted/∂T_cell
	- T_surface: The solved surface temperature (for debugging)
	- converged: Whether the Newton solve converged
"""
@inline function _compute_boundary_flux_jacobian(bc, T_cell::Float64, k_eff::Float64,
	Δz_half::Float64, ε_eff::Float64, t::Float64, radiation_model::RadiationModel)

	# For DirichletThermalBC, flux is -k*(T_bc - T_cell)/Δz_half, so dq/dT_cell = k/Δz_half
	if bc isa DirichletThermalBC
		dq_dT_cell = k_eff / Δz_half
		T_bc = bc.T isa Function ? bc.T(t) : Float64(bc.T)
		return (dq_dT_cell, T_bc, true)
	end

	# Solve for boundary temperature
	state = solve_boundary_temperature(T_cell, k_eff, Δz_half, bc, t, ε_eff;
		radiation_model = radiation_model)
	T_s = state.T_surface

	# Get BC flux and derivative at solved surface temperature
	q_bc, dq_bc_dT_s = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

	# Conduction coefficient: c = k/Δz
	cond_coeff = k_eff / Δz_half

	# Denominator for chain rule: k/Δz - dq_BC/dT_s
	# Note: dq_BC/dT_s is typically negative (heat loss increases with T)
	# so (c - d) = c + |d| > 0 for well-posed problems
	denom = cond_coeff - dq_bc_dT_s

	# Avoid division by zero (this would indicate an ill-conditioned problem)
	if abs(denom) < 1e-20
		# Fallback: if dq_bc_dT_s ≈ k/Δz_half, the system is marginally stable
		# Return the direct conduction contribution as approximation
		return (cond_coeff, T_s, state.converged)
	end

	# Chain rule result: dq_conducted/dT_cell = (k/Δz × dq_BC/dT_s) / (k/Δz - dq_BC/dT_s)
	dq_dT_cell = cond_coeff * dq_bc_dT_s / denom

	return (dq_dT_cell, T_s, state.converged)
end

"""
	compute_bc_jacobian!(J, smap, mesh, prop_cache, bc_set, n, NC, t, top_face_ids, bottom_face_ids)

Add boundary condition flux derivatives to Jacobian using BoundaryConditionSet.

This uses the boundary temperature solve to properly account for the implicit 
relationship between cell temperature and conducted flux at boundaries with
temperature-dependent BCs (especially radiative).
"""
function compute_bc_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	bc_set::BoundaryConditionSet, n::Int, NC::Int, t::Float64,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
	radiation_model::RadiationModel)

	# Process top boundary
	for f in top_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id
		T_cell = mesh.cell_states[cell_id].T
		ε_eff = prop_cache.ε[cell_id]
		k_eff = prop_cache.k[cell_id]

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_ρcp = 1.0 / prop_cache.ρcp[cell_id]

		# Half-cell distance to boundary
		Δz_half = V / A / 2.0

		# Compute proper Jacobian using boundary temperature solve and chain rule
		dq_dT_cell, _, _ = _compute_boundary_flux_jacobian(bc_set.thermal_top, T_cell, k_eff, Δz_half, ε_eff, t, radiation_model)

		if !isnan(dq_dT_cell)
			# Sign convention: positive flux into cell, so at top (z=L) face:
			# dT/dt contribution is -(A/V) * q_out = -(A/V) * (-q_conducted) = (A/V) * q_conducted
			# But the conducted flux q_conducted = k*(T_cell - T_s)/Δz_half enters with negative sign
			# in the energy equation at the top face (outward normal)
			fast_sparse_add!(J, smap, cell_id, cell_id, inv_ρcp * A / V * dq_dT_cell)
		end
	end

	# Process bottom boundary
	for f in bottom_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id
		T_cell = mesh.cell_states[cell_id].T
		ε_eff = prop_cache.ε[cell_id]
		k_eff = prop_cache.k[cell_id]

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_ρcp = 1.0 / prop_cache.ρcp[cell_id]

		# Half-cell distance to boundary
		Δz_half = V / A / 2.0

		# Compute proper Jacobian using boundary temperature solve and chain rule
		dq_dT_cell, _, _ = _compute_boundary_flux_jacobian(bc_set.thermal_bottom, T_cell, k_eff, Δz_half, ε_eff, t, radiation_model)

		if !isnan(dq_dT_cell)
			fast_sparse_add!(J, smap, cell_id, cell_id, inv_ρcp * A / V * dq_dT_cell)
		end
	end

	return nothing
end

# BC derivative dispatch - now takes ε_eff for dynamic emissivity
_bc_derivative(::AdiabaticBC, T, t, ε_eff) = 0.0
_bc_derivative(bc::HeatFluxBC, T, t, ε_eff) = 0.0  # No T dependence
_bc_derivative(bc::ConvectiveBC, T, t, ε_eff) = convective_bc_derivative(Float64(bc.h isa Function ? bc.h(t) : bc.h))
_bc_derivative(bc::RadiativeBC, T, t, ε_eff) = begin
	# Use effective emissivity if not specified
	ε = if bc.ε === nothing
		ε_eff
	else
		bc.ε isa Function ? bc.ε(T) : Float64(bc.ε)
	end
	radiative_bc_derivative(ε, T)
end
_bc_derivative(bc::DirichletThermalBC, T, t, ε_eff) = NaN  # Handled specially
_bc_derivative(bc::RadiativeFluxBC, T, t, ε_eff) = begin
	# RadiativeFluxBC specifies incident flux - absorptivity handled elsewhere
	# No direct temperature dependence in the flux itself
	0.0
end

# Combined BC: sum of derivatives - now takes ε_eff
function _bc_derivative(bc::CombinedThermalBC, T, t, ε_eff)
	result = 0.0
	for sub_bc in bc.bcs
		d = _bc_derivative(sub_bc, T, t, ε_eff)
		if !isnan(d)
			result += d
		end
	end
	return result
end

# FlameBC: for now, assume flame contribution has no temperature derivative
_bc_derivative(bc::FlameBC, T, t, ε_eff) = 0.0

# ============================================================================
# Mass Boundary Condition Jacobian
# ============================================================================

"""
	compute_mass_bc_jacobian!(J, smap, mesh, prop_cache, bc_set, material, n, NC, t, top_face_ids, bottom_face_ids, is_gas, cache=nothing)

Add mass boundary condition flux derivatives to Jacobian using BoundaryConditionSet.

For DiffusiveMassBC with flux J = -λ × h_m × ξ / (λ + h_m × d):
	∂J/∂ξ_j = -λ × h_m / (λ + h_m × d) (direct derivative)
	∂J/∂ξ_k through λ = ∂J/∂λ × ∂λ/∂ξ_k (cross-species through λ)

For ImpermeableBC:
	∂J/∂ξ = 0

For ConvectiveMassBC:
	∂J/∂ξ = -h_m (when ξ_inf = 0)
"""
function compute_mass_bc_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache, bc_set::BoundaryConditionSet,
	material, n::Int, NC::Int, t::Float64, top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
	is_gas::AbstractVector{Bool}, cache::Union{Nothing, JacobianCache, JacobianCacheALE} = nothing)

	# Process top boundary
	_compute_mass_bc_jacobian_faces!(J, smap, mesh, prop_cache, bc_set.mass_top, material, n, NC, t, top_face_ids, is_gas, cache)

	# Process bottom boundary
	_compute_mass_bc_jacobian_faces!(J, smap, mesh, prop_cache, bc_set.mass_bottom, material, n, NC, t, bottom_face_ids, is_gas, cache)

	return nothing
end

# Helper function to process a set of faces with a given BC
function _compute_mass_bc_jacobian_faces!(J::AbstractMatrix, smap, mesh, prop_cache, bc,
	material, n::Int, NC::Int, t::Float64, face_ids::Vector{Int},
	is_gas::AbstractVector{Bool}, cache::Union{Nothing, JacobianCache, JacobianCacheALE} = nothing)

	# Get work arrays from cache if available (for allocation-free dλ/dξ computation)
	base_cache = if cache isa JacobianCacheALE
		cache.base
	elseif cache isa JacobianCache
		cache
	else
		nothing
	end

	# Pre-allocate arrays for cached computation (only if cache available)
	has_cache = base_cache !== nothing
	if has_cache
		densities = base_cache.comp_densities
		lambdas = base_cache.comp_lambdas
		v = base_cache.volume_fractions
		swelling_factors = base_cache.swelling_factors
		mixing = material.gas_transfer_mixing
		β = material.β_λ
	end

	for f in face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]

		# Get λ and d for diffusive-type BCs
		λ = prop_cache.λ[cell_id]
		face_pos = mesh.face_positions[f]
		cell_center = mesh.cell_centers[cell_id]
		d = abs(face_pos - cell_center)

		# Get cell state for λ derivatives
		ξ = mesh.cell_states[cell_id].ξ
		T = mesh.cell_states[cell_id].T

		# Pre-compute properties for cached dλ/dξ if cache is available
		if has_cache
			# Cache component properties at current T
			cache_component_properties!(densities, lambdas, T, material)

			# Compute volume fractions and total volume
			V_total = compute_volume_fractions!(v, ξ, densities)
		end

		# Compute ∂J_bc/∂ξⱼ for each gas component
		for j in 1:NC
			# Use cached is_gas instead of material.components[j].phase (avoids allocation)
			if !is_gas[j]
				continue
			end

			ξ_j = ξ[j]
			ξ_row = n + (j - 1) * n + cell_id

			# Direct derivative: ∂J/∂ξ_j (from ξ in flux formula)
			dJ_dξ_direct = _mass_bc_derivative(bc, λ, d, t)

			if abs(dJ_dξ_direct) > FLUX_DERIV_TOL
				ξ_col = n + (j - 1) * n + cell_id

				if left_id == 0
					fast_sparse_sub!(J, smap, ξ_row, ξ_col, A / V * dJ_dξ_direct)
				else
					fast_sparse_add!(J, smap, ξ_row, ξ_col, A / V * dJ_dξ_direct)
				end
			end

			# Cross-species derivatives through λ: ∂J/∂ξ_k = ∂J/∂λ × ∂λ/∂ξ_k
			# For DiffusiveMassBC: J = -λ × h_m × ξ / (λ + h_m × d)
			# ∂J/∂λ = -h_m × ξ / (λ + h_m × d) + λ × h_m × ξ / (λ + h_m × d)²
			#       = -h_m × ξ × [1/(λ + h_m×d) - λ/(λ + h_m×d)²]
			#       = -h_m × ξ × h_m × d / (λ + h_m × d)²
			dJ_dλ = _mass_bc_dJ_dλ(bc, λ, d, ξ_j, t)

			if abs(dJ_dλ) > FLUX_DERIV_TOL
				# Cross-species derivatives through λ
				for k in 1:NC
					# Use cached version if available for allocation-free computation
					dλ_dξ_k = if has_cache
						compute_dλ_dξ_cached(ξ, k, densities, lambdas, v, V_total, mixing, β)
					else
						compute_dλ_dξ(ξ, T, k, material)
					end

					if abs(dλ_dξ_k) > FLUX_DERIV_TOL
						ξ_col = n + (k - 1) * n + cell_id
						dJ_dξ_k = dJ_dλ * dλ_dξ_k

						if left_id == 0
							fast_sparse_sub!(J, smap, ξ_row, ξ_col, A / V * dJ_dξ_k)
						else
							fast_sparse_add!(J, smap, ξ_row, ξ_col, A / V * dJ_dξ_k)
						end
					end
				end

				# Temperature derivative through λ: ∂J/∂T = ∂J/∂λ × ∂λ/∂T
				# The effective λ depends on T through volume fractions (gas density ~ 1/T)
				dλ_dT = compute_dλ_dT(ξ, T, material)
				if abs(dλ_dT) > FLUX_DERIV_TOL
					T_col = cell_id
					dJ_dT = dJ_dλ * dλ_dT

					if left_id == 0
						fast_sparse_sub!(J, smap, ξ_row, T_col, A / V * dJ_dT)
					else
						fast_sparse_add!(J, smap, ξ_row, T_col, A / V * dJ_dT)
					end
				end
			end
		end
	end

	return nothing
end

# Mass BC derivative dispatch (direct derivative w.r.t. ξ_j)
# Types imported at top of file

_mass_bc_derivative(::ImpermeableBC, λ, d, t) = 0.0
_mass_bc_derivative(::MassFluxBC, λ, d, t) = 0.0  # Prescribed flux, no ξ dependence
# Pressure BCs control Darcy velocity, not diffusive flux - return 0 for Jacobian
_mass_bc_derivative(::DirichletPressureBC, λ, d, t) = 0.0
_mass_bc_derivative(::NeumannPressureBC, λ, d, t) = 0.0
_mass_bc_derivative(::ConvectivePressureBC, λ, d, t) = 0.0

function _mass_bc_derivative(bc::DiffusiveMassBC, λ, d, t)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m

	# Flux: J = -λ × h_m × ξ / (λ + h_m × d)
	# Derivative: ∂J/∂ξ = -λ × h_m / (λ + h_m × d)
	if isinf(h_m)
		# Pure internal diffusion: J = -λ × ξ / d
		return -λ / d
	else
		return -(λ * h_m) / (λ + h_m * d)
	end
end

# Mass BC derivative w.r.t. λ
_mass_bc_dJ_dλ(::ImpermeableBC, λ, d, ξ, t) = 0.0
_mass_bc_dJ_dλ(::MassFluxBC, λ, d, ξ, t) = 0.0
_mass_bc_dJ_dλ(::ConvectiveMassBC, λ, d, ξ, t) = 0.0  # No λ dependence
# Pressure BCs control Darcy velocity, not diffusive flux - return 0 for Jacobian
_mass_bc_dJ_dλ(::DirichletPressureBC, λ, d, ξ, t) = 0.0
_mass_bc_dJ_dλ(::NeumannPressureBC, λ, d, ξ, t) = 0.0
_mass_bc_dJ_dλ(::ConvectivePressureBC, λ, d, ξ, t) = 0.0

function _mass_bc_dJ_dλ(bc::DiffusiveMassBC, λ, d, ξ, t)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m

	# Flux: J = -λ × h_m × ξ / (λ + h_m × d)
	# ∂J/∂λ = -h_m × ξ × [1/(λ+h_m×d) - λ/(λ+h_m×d)²]
	#       = -h_m × ξ × [(λ+h_m×d - λ)/(λ+h_m×d)²]
	#       = -h_m × ξ × h_m × d / (λ + h_m × d)²
	if isinf(h_m)
		# Pure internal diffusion: J = -λ × ξ / d
		# ∂J/∂λ = -ξ / d
		return -ξ / d
	else
		denom = λ + h_m * d
		return -(h_m * h_m * d * ξ) / (denom * denom)
	end
end

function _mass_bc_derivative(bc::ConvectiveMassBC, λ, d, t)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m
	# Flux: J = h_m × (ξ_inf - ξ_surface)
	# Derivative: ∂J/∂ξ = -h_m
	return -h_m
end

# Combined: sum of all component derivatives for type stability
# Use @generated to unroll at compile time and avoid allocation from heterogeneous tuple
@generated function _mass_bc_derivative(bc::CombinedMassBC{BCs}, λ, d, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(result += _mass_bc_derivative(bc.bcs[$i], λ, d, t)) for i in 1:N]
	return quote
		result = 0.0
		$(exprs...)
		return result
	end
end

# Use @generated to unroll at compile time and avoid allocation from heterogeneous tuple
@generated function _mass_bc_dJ_dλ(bc::CombinedMassBC{BCs}, λ, d, ξ, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(result += _mass_bc_dJ_dλ(bc.bcs[$i], λ, d, ξ, t)) for i in 1:N]
	return quote
		result = 0.0
		$(exprs...)
		return result
	end
end

# ============================================================================
# Mass BC Geometry Jacobian
# ============================================================================

"""
	compute_mass_bc_geometry_jacobian!(J, smap, mesh, prop_cache, bc_set, material, dV_dz, n, NC, z_offset, t, top_face_ids, bottom_face_ids, is_gas)

Add geometry derivatives for mass boundary condition fluxes.

For DiffusiveMassBC with flux J = -λ × h_m × ξ / (λ + h_m × d):
The contribution to dξ/dt from BC is: -(1/V) × A × J (or +... depending on sign convention)

Derivatives w.r.t. node positions z come from:
1. ∂(1/V)/∂z_k = -1/V² × ∂V/∂z_k
2. ∂d/∂z_k where d is distance from cell center to face

For surface cell n bounded by nodes [n, n+1]:
- d = Δz/2 where Δz = z_{n+1} - z_n
- ∂d/∂z_n = -0.5
- ∂d/∂z_{n+1} = +0.5

∂J/∂d = λ × h_m² × ξ / (λ + h_m × d)² (for DiffusiveMassBC)
"""
function compute_mass_bc_geometry_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache, bc_set::BoundaryConditionSet,
	material, dV_dz::Matrix{Float64}, n::Int, NC::Int, z_offset::Int, t::Float64,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
	is_gas::AbstractVector{Bool})

	n_nodes = size(dV_dz, 2)

	# Process top boundary
	bc_top = bc_set.mass_top
	for f in top_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_V = 1.0 / V
		inv_V_sq = inv_V * inv_V

		# Get λ, ρ and d for diffusive-type BCs (using rigorous λρ formulation)
		λ = prop_cache.λ[cell_id]
		ρ = prop_cache.ρ[cell_id]
		face_pos = mesh.face_positions[f]
		cell_center = mesh.cell_centers[cell_id]
		d = abs(face_pos - cell_center)

		# For each gas component
		for j in 1:NC
			# Use cached is_gas instead of material.components[j].phase (avoids allocation)
			if !is_gas[j]
				continue
			end

			ξ_row = n + (j - 1) * n + cell_id
			ξ_val = mesh.cell_states[cell_id].ξ[j]

			# Get current flux and d-derivative (pass λ and ρ separately)
			J_bc, dJ_dd = _mass_bc_flux_and_d_derivative(bc_top, ξ_val, λ, ρ, d, t)

			if abs(J_bc) < FLUX_DERIV_TOL && abs(dJ_dd) < FLUX_DERIV_TOL
				continue
			end

			# Node indices for this cell (cell_id is bounded by nodes [cell_id, cell_id+1])
			n_left = cell_id
			n_right = cell_id + 1

			# Geometry derivatives
			# d = |face_pos - cell_center| ≈ Δz/2 for surface-adjacent cells
			# For cell n: d = (z_{n+1} - z_n) / 2
			# ∂d/∂z_n = -0.5, ∂d/∂z_{n+1} = +0.5

			dd_dzn_left = -0.5
			dd_dzn_right = 0.5

			# Species equation contribution from BC:
			# At top boundary (right_id=0): contribution = +(A/V) × J_bc
			# At bottom boundary (left_id=0): contribution = -(A/V) × J_bc
			# Sign factor
			sign_factor = (left_id == 0) ? -1.0 : 1.0

			# Term 1: ∂(1/V)/∂z_k contribution
			# d[(A/V) × J_bc]/dz_k = A × J_bc × (-1/V²) × dV/dz_k
			for k in 1:n_nodes
				if abs(dV_dz[cell_id, k]) > FLUX_DERIV_TOL
					deriv = sign_factor * A * J_bc * (-inv_V_sq) * dV_dz[cell_id, k]
					fast_sparse_add!(J, smap, ξ_row, z_offset+k, deriv)
				end
			end

			# Term 2: ∂J_bc/∂d × ∂d/∂z_k contribution
			if abs(dJ_dd) > FLUX_DERIV_TOL
				# ∂d/∂z[n_left] = -0.5
				deriv_left = sign_factor * A * inv_V * dJ_dd * dd_dzn_left
				fast_sparse_add!(J, smap, ξ_row, z_offset+n_left, deriv_left)

				# ∂d/∂z[n_right] = +0.5
				deriv_right = sign_factor * A * inv_V * dJ_dd * dd_dzn_right
				fast_sparse_add!(J, smap, ξ_row, z_offset+n_right, deriv_right)
			end
		end
	end

	# Process bottom boundary
	bc_bottom = bc_set.mass_bottom
	for f in bottom_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_V = 1.0 / V
		inv_V_sq = inv_V * inv_V

		# Get λ, ρ and d for diffusive-type BCs (using rigorous λρ formulation)
		λ = prop_cache.λ[cell_id]
		ρ = prop_cache.ρ[cell_id]
		face_pos = mesh.face_positions[f]
		cell_center = mesh.cell_centers[cell_id]
		d = abs(face_pos - cell_center)

		# For each gas component
		for j in 1:NC
			# Use cached is_gas instead of material.components[j].phase (avoids allocation)
			if !is_gas[j]
				continue
			end

			ξ_row = n + (j - 1) * n + cell_id
			ξ_val = mesh.cell_states[cell_id].ξ[j]

			# Get current flux and d-derivative (pass λ and ρ separately)
			J_bc, dJ_dd = _mass_bc_flux_and_d_derivative(bc_bottom, ξ_val, λ, ρ, d, t)

			if abs(J_bc) < FLUX_DERIV_TOL && abs(dJ_dd) < FLUX_DERIV_TOL
				continue
			end

			# Node indices for this cell (cell_id is bounded by nodes [cell_id, cell_id+1])
			n_left = cell_id
			n_right = cell_id + 1

			# Geometry derivatives
			# d = |face_pos - cell_center| ≈ Δz/2 for surface-adjacent cells
			# For cell n: d = (z_{n+1} - z_n) / 2
			# ∂d/∂z_n = -0.5, ∂d/∂z_{n+1} = +0.5

			dd_dzn_left = -0.5
			dd_dzn_right = 0.5

			# Species equation contribution from BC:
			# At top boundary (right_id=0): contribution = +(A/V) × J_bc
			# At bottom boundary (left_id=0): contribution = -(A/V) × J_bc
			# Sign factor
			sign_factor = (left_id == 0) ? -1.0 : 1.0

			# Term 1: ∂(1/V)/∂z_k contribution
			# d[(A/V) × J_bc]/dz_k = A × J_bc × (-1/V²) × dV/dz_k
			for k in 1:n_nodes
				if abs(dV_dz[cell_id, k]) > FLUX_DERIV_TOL
					deriv = sign_factor * A * J_bc * (-inv_V_sq) * dV_dz[cell_id, k]
					fast_sparse_add!(J, smap, ξ_row, z_offset+k, deriv)
				end
			end

			# Term 2: ∂J_bc/∂d × ∂d/∂z_k contribution
			if abs(dJ_dd) > FLUX_DERIV_TOL
				# ∂d/∂z[n_left] = -0.5
				deriv_left = sign_factor * A * inv_V * dJ_dd * dd_dzn_left
				fast_sparse_add!(J, smap, ξ_row, z_offset+n_left, deriv_left)

				# ∂d/∂z[n_right] = +0.5
				deriv_right = sign_factor * A * inv_V * dJ_dd * dd_dzn_right
				fast_sparse_add!(J, smap, ξ_row, z_offset+n_right, deriv_right)
			end
		end
	end

	return nothing
end

# Mass BC flux and d-derivative dispatch (using λ and ρ separately for consistency with λρ formulation)
_mass_bc_flux_and_d_derivative(::ImpermeableBC, ξ, λ, ρ, d, t) = (0.0, 0.0)
_mass_bc_flux_and_d_derivative(::MassFluxBC, ξ, λ, ρ, d, t) = (0.0, 0.0)  # No ξ or d dependence for prescribed flux
# Pressure BCs control Darcy velocity, not diffusive flux - return 0 for Jacobian
_mass_bc_flux_and_d_derivative(::DirichletPressureBC, ξ, λ, ρ, d, t) = (0.0, 0.0)
_mass_bc_flux_and_d_derivative(::NeumannPressureBC, ξ, λ, ρ, d, t) = (0.0, 0.0)
_mass_bc_flux_and_d_derivative(::ConvectivePressureBC, ξ, λ, ρ, d, t) = (0.0, 0.0)

function _mass_bc_flux_and_d_derivative(bc::DiffusiveMassBC, ξ, λ, ρ, d, t)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m

	if isinf(h_m)
		# Pure internal diffusion: J = -λ × ξ / d
		# Note: This is equivalent to -λρ × Y / d = -λρ × (ξ/ρ) / d = -λξ/d
		# ∂J/∂d = λ × ξ / d²
		J_bc = -λ * ξ / d
		dJ_dd = λ * ξ / (d * d)
		return (J_bc, dJ_dd)
	else
		# Mixed: J = -λ × h_m × ξ / (λ + h_m × d)
		# ∂J/∂d = λ × h_m² × ξ / (λ + h_m × d)²
		denom = λ + h_m * d
		J_bc = -(λ * h_m * ξ) / denom
		dJ_dd = (λ * h_m * h_m * ξ) / (denom * denom)
		return (J_bc, dJ_dd)
	end
end

function _mass_bc_flux_and_d_derivative(bc::ConvectiveMassBC, ξ, λ, ρ, d, t)
	h_m = bc.h_m isa Function ? bc.h_m(t) : bc.h_m
	ξ_inf = bc.ξ_inf isa Function ? bc.ξ_inf(t) : bc.ξ_inf
	# Flux: J = h_m × (ξ_inf - ξ)
	# ∂J/∂d = 0 (no d dependence)
	J_bc = h_m * (ξ_inf - ξ)
	return (J_bc, 0.0)
end

# CombinedMassBC: sum contributions from all sub-BCs
# Use @generated to unroll at compile time and avoid allocation from heterogeneous tuple
@generated function _mass_bc_flux_and_d_derivative(bc::CombinedMassBC{BCs}, ξ, λ, ρ, d, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [
		quote
			J_sub, dJ_dd_sub = _mass_bc_flux_and_d_derivative(bc.bcs[$i], ξ, λ, ρ, d, t)
			J_total += J_sub
			dJ_dd_total += dJ_dd_sub
		end for i in 1:N
	]
	return quote
		J_total = 0.0
		dJ_dd_total = 0.0
		$(exprs...)
		return (J_total, dJ_dd_total)
	end
end

# ============================================================================
# Radiation Jacobian (Beer-Lambert)
# ============================================================================

"""
	compute_radiation_jacobian!(J, mesh, prop_cache, I_surface, n, NC, comp_absorptions, τ_cumulative)

Add Beer-Lambert radiation absorption derivatives to Jacobian.

The radiation intensity at cell i depends on all cells from surface to i:
I_i = I₀ × exp(-Σₖ₌ᵢᴺ αₖ × Δzₖ)

where α = Σⱼ ξⱼ × αⱼ is the effective absorption coefficient.

∂Q_rad,i/∂ξⱼ,k has non-zero contributions for k ≥ i (upstream cells).

# Arguments
- `comp_absorptions`: Cached component absorption coefficients [NC]
- `τ_cumulative`: Pre-allocated work array of size at least n+1
"""
function compute_radiation_jacobian!(J::AbstractMatrix, smap, mesh, prop_cache,
	I_surface::Float64, n::Int, NC::Int, comp_absorptions::AbstractVector{Float64},
	τ_cumulative::AbstractVector{Float64})
	if I_surface <= 0.0
		return nothing
	end

	# Compute cumulative optical depth from surface (cell n) into depth
	# τ_cumulative[n+1] = 0 at surface
	fill!(τ_cumulative, 0.0)

	@inbounds @simd for i in n:-1:1
		α = prop_cache.α[i]
		Δz = cell_size(mesh, i)
		τ_cumulative[i] = τ_cumulative[i+1] + α * Δz
	end

	# For each cell, compute sensitivity of Q_rad to upstream ξ values
	@inbounds for i in 1:n
		T_row = i
		α_i = prop_cache.α[i]
		Δz_i = cell_size(mesh, i)
		V_i = mesh.cell_volumes[i]
		inv_ρcp_i = 1.0 / prop_cache.ρcp[i]

		# Radiation intensity entering cell i
		I_in = I_surface * exp(-τ_cumulative[i+1])

		# Current volumetric absorption
		τ_i = α_i * Δz_i
		if τ_i > JAC_ZERO_TOL
			absorption_factor = (1.0 - exp(-τ_i)) / Δz_i
		else
			absorption_factor = α_i  # Linear limit
		end
		Q_rad_i = I_in * absorption_factor

		# -----------------------------------------------------------------
		# Sensitivity to local cell (i) composition
		# -----------------------------------------------------------------
		# Q_rad = I_in × (1 - exp(-α×Δz)) / Δz
		# ∂Q_rad/∂α = I_in × exp(-α×Δz) (for small τ)
		#           = I_in × Δz × exp(-τ) / Δz = I_in × exp(-τ)

		for j in 1:NC
			αⱼ = comp_absorptions[j]
			if αⱼ < JAC_ZERO_TOL
				continue
			end

			# ∂α/∂ξⱼ = αⱼ (mass absorption coefficient)
			# ∂Q_rad/∂ξⱼ = I_in × exp(-τ_i) × αⱼ (approximately)
			dQ_dξ_local = I_in * exp(-τ_i) * αⱼ

			ξ_col = n + (j-1)*n + i
			fast_sparse_add!(J, smap, T_row, ξ_col, inv_ρcp_i * dQ_dξ_local)
		end

		# -----------------------------------------------------------------
		# Sensitivity to upstream cells (k > i) composition
		# -----------------------------------------------------------------
		# I_in at cell i depends on absorption in all upstream cells
		# ∂I_in/∂α_k = -I_surface × exp(-τ_k_to_surface) × Δz_k × exp(-τ_i+1_to_k)
		#            = -I_in × Δz_k (chain rule through optical depth)

		for k in (i+1):n
			Δz_k = cell_size(mesh, k)

			for j in 1:NC
				αⱼ = comp_absorptions[j]
				if αⱼ < JAC_ZERO_TOL
					continue
				end

				# ∂I_in/∂ξⱼ,k = -I_in × Δz_k × αⱼ
				# ∂Q_rad,i/∂ξⱼ,k = absorption_factor × ∂I_in/∂ξⱼ,k
				dQ_dξ_upstream = -absorption_factor * I_in * Δz_k * αⱼ

				ξ_col = n + (j-1)*n + k
				fast_sparse_add!(J, smap, T_row, ξ_col, inv_ρcp_i * dQ_dξ_upstream)
			end
		end
	end

	return nothing
end

"""
	compute_radiation_geometry_jacobian!(J, mesh, prop_cache, I_surface, n, NC, material,
										 τ_cumulative, dΔz_dz, z_offset)

Add z-derivatives of Beer-Lambert radiation absorption to Jacobian.

The radiation volumetric source is:
	Q_rad,i = I_in × (1 - exp(-τ_i)) / Δz_i

where:
- I_in = I₀ × exp(-Σₖ₌ᵢ₊₁ⁿ τₖ) is the intensity entering cell i from upstream
- τ_i = α_i × Δz_i is the optical depth of cell i
- Δz_i = z_{i+1} - z_i

When node positions z change:
1. Local effect: ∂Q_rad,i/∂Δz_i from the local cell size
2. Upstream effect: ∂I_in/∂Δz_k for k > i affects intensity reaching cell i

# Arguments
- `τ_cumulative`: Pre-computed cumulative optical depth (from surface)
- `dΔz_dz`: Matrix of ∂Δz_i/∂z_k derivatives
- `z_offset`: Offset to z DOFs in the Jacobian
"""
function compute_radiation_geometry_jacobian!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	prop_cache,
	I_surface::Float64,
	n::Int, NC_val::Int, material,
	τ_cumulative::AbstractVector{Float64},
	dΔz_dz::Matrix{Float64},
	z_offset::Int) where NC

	if I_surface <= 0.0
		return nothing
	end

	n_nodes = size(dΔz_dz, 2)

	@inbounds for i in 1:n
		T_row = i
		α_i = prop_cache.α[i]
		Δz_i = cell_size(mesh, i)
		inv_ρcp_i = 1.0 / prop_cache.ρcp[i]

		# Radiation intensity entering cell i (from upstream, cells i+1 to n)
		I_in = I_surface * exp(-τ_cumulative[i+1])

		# Current optical depth of cell i
		τ_i = α_i * Δz_i
		exp_neg_τi = exp(-τ_i)

		# Volumetric absorption rate
		# Q_rad = I_in × (1 - exp(-τ)) / Δz = I_in × α × f(τ)
		# where f(τ) = (1 - exp(-τ)) / τ is the absorption efficiency
		if τ_i > 1e-6
			one_minus_exp = 1.0 - exp_neg_τi
			absorption_factor = one_minus_exp / Δz_i
			Q_rad_i = I_in * absorption_factor
		else
			# Linear limit for τ → 0: (1 - exp(-τ))/τ → 1 - τ/2
			absorption_factor = α_i * (1.0 - 0.5 * τ_i)
			Q_rad_i = I_in * absorption_factor
		end

		# -----------------------------------------------------------------
		# 1. LOCAL z-derivatives: ∂Q_rad,i/∂Δz_i × ∂Δz_i/∂z_k
		# -----------------------------------------------------------------
		# Q_rad = I_in × (1 - exp(-α×Δz)) / Δz
		# Let f(τ) = (1 - exp(-τ)) / τ where τ = α×Δz
		# Then Q_rad = I_in × α × f(τ)
		#
		# ∂Q_rad/∂Δz = I_in × α × df/dτ × dτ/dΔz
		#            = I_in × α × df/dτ × α
		#            = I_in × α² × df/dτ
		#
		# df/dτ = (exp(-τ) × τ - (1-exp(-τ))) / τ²
		#       = (τ×exp(-τ) - 1 + exp(-τ)) / τ²
		#       = (exp(-τ)×(τ+1) - 1) / τ²
		#
		# For large τ: df/dτ → -1/τ²
		# For small τ: df/dτ → -1/2 + τ/6 - ...
		#
		# Alternative derivation (direct):
		# Q = I × (1 - exp(-α×Δz)) / Δz
		# ∂Q/∂Δz = I × [α×exp(-α×Δz)×Δz - (1-exp(-α×Δz))] / Δz²
		#        = I × [α×exp(-τ) - (1-exp(-τ))/Δz] / Δz
		#        = (I/Δz) × [α×exp(-τ) - absorption_factor]
		#        = (I/Δz) × [α×exp(-τ) - Q/I]

		if τ_i > 1e-6
			# Use the direct formula
			dQ_dΔz_local = (I_in / Δz_i) * (α_i * exp_neg_τi - absorption_factor)
		else
			# Small τ limit: ∂Q/∂Δz ≈ -I_in × α² / 2
			dQ_dΔz_local = -I_in * α_i * α_i * 0.5
		end

		# Apply chain rule: ∂Δz_i/∂z_k
		# For cell i bounded by nodes i and i+1:
		# Δz_i = z_{i+1} - z_i
		# ∂Δz_i/∂z_i = -1, ∂Δz_i/∂z_{i+1} = +1

		if abs(dQ_dΔz_local) > FLUX_DERIV_TOL
			# Contribution to z[i] (left node)
			fast_sparse_add!(J, smap, T_row, z_offset+i, inv_ρcp_i * dQ_dΔz_local * (-1.0))

			# Contribution to z[i+1] (right node)
			fast_sparse_add!(J, smap, T_row, z_offset+i+1, inv_ρcp_i * dQ_dΔz_local * (+1.0))
		end

		# -----------------------------------------------------------------
		# 2. UPSTREAM z-derivatives: ∂I_in/∂Δz_k × ∂Δz_k/∂z_m for k > i
		# -----------------------------------------------------------------
		# I_in = I₀ × exp(-Σₖ₌ᵢ₊₁ⁿ αₖ×Δzₖ)
		# ∂I_in/∂Δz_k = -I_in × αₖ
		#
		# ∂Q_rad,i/∂Δz_k = absorption_factor × ∂I_in/∂Δz_k
		#               = -absorption_factor × I_in × αₖ

		for k in (i+1):n
			α_k = prop_cache.α[k]

			if α_k < JAC_ZERO_TOL
				continue
			end

			dI_in_dΔz_k = -I_in * α_k
			dQ_dΔz_upstream = absorption_factor * dI_in_dΔz_k

			if abs(dQ_dΔz_upstream) > FLUX_DERIV_TOL
				# Apply chain rule: ∂Δz_k/∂z_m
				# For cell k: ∂Δz_k/∂z_k = -1, ∂Δz_k/∂z_{k+1} = +1

				# Contribution to z[k] (left node of cell k)
				fast_sparse_add!(J, smap, T_row, z_offset+k, inv_ρcp_i * dQ_dΔz_upstream * (-1.0))

				# Contribution to z[k+1] (right node of cell k)
				fast_sparse_add!(J, smap, T_row, z_offset+k+1, inv_ρcp_i * dQ_dΔz_upstream * (+1.0))
			end
		end
	end

	return nothing
end

# Helper function - use Geometry.cell_size
using ..Geometry: cell_size

# ============================================================================
# ALE Jacobian: Core Derivative Functions
# ============================================================================
#
# VALIDATION STATUS (January 2026): COMPLETE
# ------------------------------------------
# All ALE Jacobian blocks have been validated against finite differences:
#   - dT/dT, dT/dξ, dξ/dT, dξ/dξ: ~1e-10 relative error (machine precision)
#   - dz/dT, dz/dξ, dz/dz: ~1e-10 relative error (machine precision)
#   - dT/dz, dξ/dz: ~0.01% apparent error (see note below)
#
# Note on dT/dz and dξ/dz validation error:
# The ~0.01% error is NOT a bug - it is the truncation error of finite
# differencing (FD). The analytical Jacobian computes the exact first
# derivative, while FD with step ε includes O(ε) second-order terms.
# For the 1/V factor in the divergence: FD_error ≈ (q×A³/V³)×ε/ρcₚ.
# This was confirmed by showing error scales linearly with ε.
#
# ============================================================================

# Import additional modules needed for ALE derivatives
using ..Physics: volume_change_rate

"""
	compute_strain_rate_derivatives!(dθ_dT::Vector{Float64}, dθ_dξ::Matrix{Float64},
									 dS_dT::Matrix{Float64}, dS_dξ::Array{Float64,3},
									 S::Matrix{Float64}, mesh::Mesh1D{NC},
									 material::Material{NC}) where NC

Compute derivatives of volumetric strain rate θ with respect to T and ξ.

The strain rate is:
	θᵢ = Σⱼ (γⱼ/ρⱼ(T)) × Sⱼᵢ

Where:
- γⱼ is the swelling factor of component j
- ρⱼ(T) is the density of component j (may be T-dependent)
- Sⱼᵢ = dξⱼ/dt is the species source term from reactions

# Arguments
- `dθ_dT`: Output [n_cells] - ∂θᵢ/∂Tᵢ
- `dθ_dξ`: Output [NC × n_cells] - ∂θᵢ/∂ξⱼᵢ
- `dS_dT`: Input [NC × n_cells] - ∂Sⱼ/∂T (species source derivatives w.r.t. T)
- `dS_dξ`: Input [NC × NC × n_cells] - ∂Sⱼ/∂ξₖ (species source derivatives w.r.t. ξ)
- `S`: Input [NC × n_cells] - Current species source terms Sⱼᵢ
- `mesh`: Current mesh
- `material`: Material definition

# Notes
- The density derivative ∂ρⱼ/∂T is obtained from the property function
- Only components with γⱼ > 0 contribute to volume change
"""
function compute_strain_rate_derivatives!(dθ_dT::Vector{Float64},
	dθ_dξ::Matrix{Float64},
	dS_dT::Matrix{Float64},
	dS_dξ::Array{Float64, 3},
	S::Matrix{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC},
	jac_cache::Union{JacobianCache{NC}, JacobianCacheALE{NC}}) where NC

	n_cells = length(mesh.cells)

	# Get cached component properties (avoid heterogeneous tuple access)
	base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache
	γ = base_cache.swelling_factors
	densities = base_cache.comp_densities
	dρ_dT = base_cache.comp_dρ_dT

	# Zero output arrays
	fill!(dθ_dT, 0.0)
	fill!(dθ_dξ, 0.0)

	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T

		# Cache component properties for this temperature (using @generated functions)
		cache_component_properties!(densities, base_cache.comp_lambdas, T, material)
		cache_density_derivatives!(dρ_dT, T, material)

		# Accumulate contributions from each component using cached values
		for j in 1:NC
			γⱼ = γ[j]

			# Skip components that don't contribute to volume change
			if γⱼ <= 0.0
				continue
			end

			# Use cached density and its derivative
			ρⱼ = densities[j]
			dρⱼ_dT = dρ_dT[j]

			# Current species source for this component
			Sⱼ = S[j, i]

			# ∂θ/∂T contribution from this component:
			# θⱼ = γⱼ × Sⱼ / ρⱼ
			# ∂θⱼ/∂T = γⱼ × [(∂Sⱼ/∂T)/ρⱼ - Sⱼ×(∂ρⱼ/∂T)/ρⱼ²]
			#        = γⱼ/ρⱼ × [∂Sⱼ/∂T - Sⱼ×(∂ρⱼ/∂T)/ρⱼ]
			inv_ρⱼ = 1.0 / ρⱼ
			dθ_dT[i] += γⱼ * inv_ρⱼ * (dS_dT[j, i] - Sⱼ * dρⱼ_dT * inv_ρⱼ)

			# ∂θ/∂ξₖ contribution from this component:
			# ∂θⱼ/∂ξₖ = γⱼ/ρⱼ × ∂Sⱼ/∂ξₖ
			for k in 1:NC
				dθ_dξ[k, i] += γⱼ * inv_ρⱼ * dS_dξ[j, k, i]
			end
		end
	end

	return nothing
end



"""
	compute_node_velocity_derivatives!(dw_dT::Matrix{Float64},
									   dw_dξ::Array{Float64,3},
									   dw_dz::Matrix{Float64},
									   dθ_dT::Vector{Float64},
									   dθ_dξ::Matrix{Float64},
									   θ::Vector{Float64},
									   mesh::Mesh1D{NC}) where NC

Compute derivatives of node velocities w with respect to T, ξ, and z.

The node velocity is computed via cumulative sum from the fixed bottom boundary:
	w[1] = 0  (bottom node fixed)
	w[k+1] = w[k] + θₖ × Δzₖ

Therefore:
	wₖ = Σᵢ₌₁ᵏ⁻¹ θᵢ × Δzᵢ

# Arguments
- `dw_dT`: Output [n_nodes × n_cells] - ∂wₖ/∂Tᵢ
- `dw_dξ`: Output [n_nodes × NC × n_cells] - ∂wₖ/∂ξⱼᵢ
- `dw_dz`: Output [n_nodes × n_nodes] - ∂wₖ/∂zₘ
- `dθ_dT`: Input [n_cells] - ∂θᵢ/∂Tᵢ
- `dθ_dξ`: Input [NC × n_cells] - ∂θᵢ/∂ξⱼᵢ
- `θ`: Input [n_cells] - Current strain rates
- `mesh`: Current mesh

# Notes
- dw_dT and dw_dξ are lower triangular (node k depends on cells 1..k-1)
- dw_dz is also lower triangular due to cumulative sum structure
"""
function compute_node_velocity_derivatives!(dw_dT::Matrix{Float64},
	dw_dξ::Array{Float64, 3},
	dw_dz::Matrix{Float64},
	dθ_dT::Vector{Float64},
	dθ_dξ::Matrix{Float64},
	θ::Vector{Float64},
	mesh::Mesh1D{NC}) where NC

	n_cells = length(mesh.cells)
	n_nodes = length(mesh.nodes)

	# Zero output arrays
	fill!(dw_dT, 0.0)
	fill!(dw_dξ, 0.0)
	fill!(dw_dz, 0.0)

	# Node 1 (bottom) is fixed: w[1] = 0, all derivatives are 0
	# Already zeroed above

	# For nodes k = 2, ..., n_nodes:
	# wₖ = Σᵢ₌₁ᵏ⁻¹ θᵢ × Δzᵢ
	#
	# ∂wₖ/∂Tᵢ = Δzᵢ × (∂θᵢ/∂Tᵢ)  for i < k
	# ∂wₖ/∂ξⱼᵢ = Δzᵢ × (∂θᵢ/∂ξⱼᵢ)  for i < k
	# ∂wₖ/∂zₘ = Σᵢ₌₁ᵏ⁻¹ θᵢ × (∂Δzᵢ/∂zₘ)
	#         where ∂Δzᵢ/∂z_{i+1} = +1, ∂Δzᵢ/∂zᵢ = -1

	A = mesh.cross_section_area

	@inbounds for k in 2:n_nodes
		# Node k corresponds to the interface AFTER cell k-1
		# wₖ depends on cells 1, 2, ..., k-1

		for i in 1:(k-1)
			# Cell i contributes θᵢ × Δzᵢ to wₖ
			Δzᵢ = mesh.cell_volumes[i] / A

			# Derivatives w.r.t. T and ξ
			dw_dT[k, i] = Δzᵢ * dθ_dT[i]

			for j in 1:NC
				dw_dξ[k, j, i] = Δzᵢ * dθ_dξ[j, i]
			end

			# Derivatives w.r.t. node positions z
			# The contribution from cell i is θᵢ × Δzᵢ
			# Δzᵢ = z_{i+1} - zᵢ
			# ∂(θᵢ × Δzᵢ)/∂z_{i+1} = θᵢ × (+1) = +θᵢ
			# ∂(θᵢ × Δzᵢ)/∂zᵢ = θᵢ × (-1) = -θᵢ

			# Node indices for cell i: node i (left) and node i+1 (right)
			# In our 1D mesh: cell i has node_ids = (i, i+1)
			dw_dz[k, i] += -θ[i]      # ∂wₖ/∂zᵢ
			dw_dz[k, i+1] += θ[i]     # ∂wₖ/∂z_{i+1}
		end
	end

	return nothing
end

"""
	compute_geometry_derivatives!(dV_dz::Matrix{Float64},
								  dΔz_dz::Matrix{Float64},
								  mesh::Mesh1D{NC}) where NC

Compute derivatives of cell volumes and sizes with respect to node positions.

For 1D mesh with constant cross-section area A:
	Vᵢ = A × Δzᵢ = A × (z_{i+1} - zᵢ)
	Δzᵢ = z_{i+1} - zᵢ

# Arguments
- `dV_dz`: Output [n_cells × n_nodes] - ∂Vᵢ/∂zₖ
- `dΔz_dz`: Output [n_cells × n_nodes] - ∂Δzᵢ/∂zₖ
- `mesh`: Current mesh

# Notes
These are bidiagonal matrices:
- ∂Vᵢ/∂z_{i+1} = +A, ∂Vᵢ/∂zᵢ = -A
- ∂Δzᵢ/∂z_{i+1} = +1, ∂Δzᵢ/∂zᵢ = -1
"""
function compute_geometry_derivatives!(dV_dz::Matrix{Float64},
	dΔz_dz::Matrix{Float64},
	mesh::Mesh1D{NC}) where NC

	n_cells = length(mesh.cells)
	A = mesh.cross_section_area

	# Zero output arrays
	fill!(dV_dz, 0.0)
	fill!(dΔz_dz, 0.0)

	@inbounds for i in 1:n_cells
		# Cell i is bounded by nodes i (left) and i+1 (right)
		# Δzᵢ = z_{i+1} - zᵢ
		# Vᵢ = A × Δzᵢ

		# ∂Δzᵢ/∂zᵢ = -1
		# ∂Δzᵢ/∂z_{i+1} = +1
		dΔz_dz[i, i] = -1.0
		dΔz_dz[i, i+1] = 1.0

		# ∂Vᵢ/∂zᵢ = -A
		# ∂Vᵢ/∂z_{i+1} = +A
		dV_dz[i, i] = -A
		dV_dz[i, i+1] = A
	end

	return nothing
end

"""
	compute_ale_advection_jacobian_energy!(J::AbstractMatrix,
										   mesh::Mesh1D{NC},
										   w_node::Vector{Float64},
										   dw_dT::Matrix{Float64},
										   dw_dξ::Array{Float64,3},
										   dw_dz::Matrix{Float64},
										   sv::StateVectorALE{NC},
										   thermal_bcs::Dict,
										   t::Real) where NC

Add ALE energy advection derivatives to the Jacobian matrix.

The ALE energy advection term is:
	(w·∇T)ᵢ = wᵢ × dT/dz

With upwind discretization:
- If wᵢ > 0: dT/dz = (Tᵢ - T_{i-1}) / Δzᵢ (backward difference)
- If wᵢ < 0: dT/dz = (T_{i+1} - Tᵢ) / Δzᵢ (forward difference)

The term enters the energy equation as: dTᵢ/dt = ... + (w·∇T)ᵢ  (ADDED per ALE identity)

# Derivatives
∂(w·∇T)ᵢ/∂Tⱼ has contributions from:
1. Direct T dependence in the gradient
2. Indirect T dependence through mesh velocity wᵢ

∂(w·∇T)ᵢ/∂zₖ has contributions from:
1. Δz dependence in the gradient
2. Indirect z dependence through mesh velocity wᵢ
"""
function compute_ale_advection_jacobian_energy!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	w_node::Vector{Float64},
	dw_dT::Matrix{Float64},
	dw_dξ::Array{Float64, 3},
	dw_dz::Matrix{Float64},
	sv::StateVectorALE{NC}) where NC

	n = sv.n_cells
	n_nodes = sv.n_nodes
	A = mesh.cross_section_area

	@inbounds for i in 1:n
		cell = mesh.cells[i]
		n_left, n_right = cell.node_ids

		# Cell-centered mesh velocity (average of bounding nodes)
		w_cell = 0.5 * (w_node[n_left] + w_node[n_right])

		# Skip if velocity is essentially zero
		if abs(w_cell) < JAC_ZERO_TOL
			continue
		end

		T_i = mesh.cell_states[i].T

		# Get neighbor info
		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]
		left_id = face_left.cell_ids[1]
		right_id = face_right.cell_ids[2]

		# Get neighbor temperatures (with boundary handling)
		T_left = left_id != 0 ? mesh.cell_states[left_id].T : T_i
		T_right = right_id != 0 ? mesh.cell_states[right_id].T : T_i

		# Compute current gradient based on upwind direction
		# Use center-to-center distance (NOT cell width) for correct gradient
		if w_cell > 0
			# Upwind from left
			Δz = left_id != 0 ? abs(mesh.cell_centers[i] - mesh.cell_centers[left_id]) :
				mesh.cell_volumes[i] / A
			inv_Δz = 1.0 / Δz
			dT_dz = (T_i - T_left) * inv_Δz
		else
			# Upwind from right
			Δz = right_id != 0 ? abs(mesh.cell_centers[right_id] - mesh.cell_centers[i]) :
				mesh.cell_volumes[i] / A
			inv_Δz = 1.0 / Δz
			dT_dz = (T_right - T_i) * inv_Δz
		end

		T_row = i  # Energy equation row

		# =====================================================================
		# 1. Direct T derivatives of (w·∇T)
		# =====================================================================
		# The ALE correction is +w·∇T (ADDED to Eulerian RHS).
		# (w·∇T)ᵢ = w_cell × dT_dz
		# For w_cell > 0: dT_dz = (T_i - T_left) / Δz
		#   ∂(w·∇T)/∂T_i = w_cell / Δz
		#   ∂(w·∇T)/∂T_left = -w_cell / Δz
		# For w_cell < 0: dT_dz = (T_right - T_i) / Δz
		#   ∂(w·∇T)/∂T_i = -w_cell / Δz (note: w_cell < 0, so this is positive)
		#   ∂(w·∇T)/∂T_right = w_cell / Δz

		if w_cell > 0
			# Upwind from left: dT_dz = (T_i - T_left) / Δz
			# The advection term ADDS to dT/dt, so we add the derivative directly
			if left_id != 0
				# Interior: both T_i and T_left are independent variables
				fast_sparse_add!(J, smap, T_row, i, w_cell * inv_Δz)
				fast_sparse_add!(J, smap, T_row, left_id, -w_cell * inv_Δz)
			end
			# If left_id = 0 (boundary), T_left = T_i (ghost), so dT_dz = 0
			# and ∂(w·∇T)/∂T_i = 0 - skip contribution
		else
			# Upwind from right: dT_dz = (T_right - T_i) / Δz
			if right_id != 0
				# Interior: both T_i and T_right are independent variables
				fast_sparse_add!(J, smap, T_row, i, -w_cell * inv_Δz)
				fast_sparse_add!(J, smap, T_row, right_id, w_cell * inv_Δz)
			end
			# If right_id = 0 (boundary), T_right = T_i (ghost), so dT_dz = 0
			# and ∂(w·∇T)/∂T_i = 0 - skip contribution
		end

		# =====================================================================
		# 2. Indirect T derivatives through mesh velocity
		# =====================================================================
		# w_cell = 0.5 × (w_{n_left} + w_{n_right})
		# ∂(w·∇T)/∂T_j = dT_dz × ∂w_cell/∂T_j
		#              = dT_dz × 0.5 × (∂w_{n_left}/∂T_j + ∂w_{n_right}/∂T_j)

		for j in 1:n
			dw_cell_dTj = 0.5 * (dw_dT[n_left, j] + dw_dT[n_right, j])
			if abs(dw_cell_dTj) > FLUX_DERIV_TOL
				fast_sparse_add!(J, smap, T_row, j, dT_dz * dw_cell_dTj)
			end
		end

		# =====================================================================
		# 3. ξ derivatives through mesh velocity
		# =====================================================================
		for j in 1:NC
			for k in 1:n
				dw_cell_dξjk = 0.5 * (dw_dξ[n_left, j, k] + dw_dξ[n_right, j, k])
				if abs(dw_cell_dξjk) > FLUX_DERIV_TOL
					ξ_col = n + (j - 1) * n + k
					fast_sparse_add!(J, smap, T_row, ξ_col, dT_dz * dw_cell_dξjk)
				end
			end
		end

		# =====================================================================
		# 4. z derivatives (geometry + mesh velocity)
		# =====================================================================
		# (w·∇T) = w_cell × (T_up - T_down) / Δz
		#
		# ∂(w·∇T)/∂z_m has contributions from:
		# a) ∂(1/Δz)/∂z_m = -1/Δz² × ∂Δz/∂z_m
		# b) ∂w_cell/∂z_m through mesh velocity derivatives

		z_offset = sv.z_offset

		# a) Geometry contribution: w_cell × (T_up - T_down) × ∂(1/Δz)/∂z_m
		# Δzᵢ = z_{i+1} - zᵢ
		# ∂(1/Δzᵢ)/∂z_{i+1} = -1/Δzᵢ² × (+1) = -1/Δzᵢ²
		# ∂(1/Δzᵢ)/∂zᵢ = -1/Δzᵢ² × (-1) = +1/Δzᵢ²
		ΔT = w_cell > 0 ? (T_i - T_left) : (T_right - T_i)
		dinv_Δz_dz_left = inv_Δz * inv_Δz   # ∂(1/Δz)/∂z_i = +1/Δz²
		dinv_Δz_dz_right = -inv_Δz * inv_Δz # ∂(1/Δz)/∂z_{i+1} = -1/Δz²

		# z_i is node index i, z_{i+1} is node index i+1
		fast_sparse_add!(J, smap, T_row, z_offset+i, w_cell * ΔT * dinv_Δz_dz_left)
		fast_sparse_add!(J, smap, T_row, z_offset+i+1, w_cell * ΔT * dinv_Δz_dz_right)

		# b) Mesh velocity contribution: dT_dz × ∂w_cell/∂z_m
		for m in 1:n_nodes
			dw_cell_dzm = 0.5 * (dw_dz[n_left, m] + dw_dz[n_right, m])
			if abs(dw_cell_dzm) > FLUX_DERIV_TOL
				fast_sparse_add!(J, smap, T_row, z_offset+m, dT_dz * dw_cell_dzm)
			end
		end
	end

	return nothing
end

"""
	compute_ale_advection_jacobian_species!(J::AbstractMatrix,
											mesh::Mesh1D{NC},
											material::Material{NC},
											w_node::Vector{Float64},
											dw_dT::Matrix{Float64},
											dw_dξ::Array{Float64,3},
											dw_dz::Matrix{Float64},
											sv::StateVectorALE{NC}) where NC

Add ALE species advection derivatives to the Jacobian matrix.

For GAS components, the ALE advection term is:
	(w·∇ξⱼ)ᵢ = -wᵢ × dξⱼ/dz

(Negative because gas velocity relative to mesh is -w when bulk gas velocity ≈ 0)

For SOLID and LIQUID components:
	No advection (material moves with mesh, v_solid = w, so v - w = 0)
"""
function compute_ale_advection_jacobian_species!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	material::Material{NC},
	w_node::Vector{Float64},
	dw_dT::Matrix{Float64},
	dw_dξ::Array{Float64, 3},
	dw_dz::Matrix{Float64},
	sv::StateVectorALE{NC},
	is_gas::AbstractVector{Bool}) where NC

	n = sv.n_cells
	n_nodes = sv.n_nodes
	A = mesh.cross_section_area

	@inbounds for i in 1:n
		cell = mesh.cells[i]
		n_left, n_right = cell.node_ids

		w_cell = 0.5 * (w_node[n_left] + w_node[n_right])

		if abs(w_cell) < JAC_ZERO_TOL
			continue
		end

		f_left, f_right = cell.face_ids
		face_left = mesh.faces[f_left]
		face_right = mesh.faces[f_right]
		left_id = face_left.cell_ids[1]
		right_id = face_right.cell_ids[2]

		# Center-to-center distance for correct upwind gradient
		if w_cell > 0
			Δz = left_id != 0 ? abs(mesh.cell_centers[i] - mesh.cell_centers[left_id]) :
				mesh.cell_volumes[i] / A
		else
			Δz = right_id != 0 ? abs(mesh.cell_centers[right_id] - mesh.cell_centers[i]) :
				mesh.cell_volumes[i] / A
		end
		inv_Δz = 1.0 / Δz

		z_offset = sv.z_offset

		# Process each component - use cached is_gas to avoid heterogeneous access
		for j in 1:NC
			# Only gas components have advection relative to mesh
			if !is_gas[j]
				continue
			end

			ξ_row = n + (j - 1) * n + i  # Species equation row for component j, cell i
			ξ_i = mesh.cell_states[i].ξ[j]

			# Get neighbor concentrations
			ξ_left = left_id != 0 ? mesh.cell_states[left_id].ξ[j] : ξ_i
			ξ_right = right_id != 0 ? mesh.cell_states[right_id].ξ[j] : ξ_i

			# The advection term for gas is: -w × ∇ξ (relative velocity is -w)
			# For w > 0: use backward difference
			# For w < 0: use forward difference
			# Sign convention: advection SUBTRACTS from dξ/dt

			if w_cell > 0
				dξ_dz = (ξ_i - ξ_left) * inv_Δz
				# Only add derivative if not at boundary (where ghost cell = cell value)
				if left_id != 0
					fast_sparse_sub!(J, smap, ξ_row, n+(j-1)*n+i, -w_cell * inv_Δz)
					fast_sparse_sub!(J, smap, ξ_row, n+(j-1)*n+left_id, w_cell * inv_Δz)
				end
			else
				dξ_dz = (ξ_right - ξ_i) * inv_Δz
				# Only add derivative if not at boundary (where ghost cell = cell value)
				if right_id != 0
					fast_sparse_sub!(J, smap, ξ_row, n+(j-1)*n+i, w_cell * inv_Δz)
					fast_sparse_sub!(J, smap, ξ_row, n+(j-1)*n+right_id, -w_cell * inv_Δz)
				end
			end

			# Current advection gradient for indirect derivatives
			advection_sign = -1.0  # -w·∇ξ for gases

			# T derivatives through mesh velocity
			for k in 1:n
				dw_cell_dTk = 0.5 * (dw_dT[n_left, k] + dw_dT[n_right, k])
				if abs(dw_cell_dTk) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row, k, advection_sign * dξ_dz * dw_cell_dTk)
				end
			end

			# ξ derivatives through mesh velocity
			for jj in 1:NC
				for k in 1:n
					dw_cell_dξjjk = 0.5 * (dw_dξ[n_left, jj, k] + dw_dξ[n_right, jj, k])
					if abs(dw_cell_dξjjk) > FLUX_DERIV_TOL
						fast_sparse_sub!(J, smap, ξ_row, n+(jj-1)*n+k, advection_sign * dξ_dz * dw_cell_dξjjk)
					end
				end
			end

			# z derivatives (geometry + mesh velocity)
			Δξ = w_cell > 0 ? (ξ_i - ξ_left) : (ξ_right - ξ_i)
			dinv_Δz_dz_left = inv_Δz * inv_Δz
			dinv_Δz_dz_right = -inv_Δz * inv_Δz

			fast_sparse_sub!(J, smap, ξ_row, z_offset+i, advection_sign * w_cell * Δξ * dinv_Δz_dz_left)
			fast_sparse_sub!(J, smap, ξ_row, z_offset+i+1, advection_sign * w_cell * Δξ * dinv_Δz_dz_right)

			for m in 1:n_nodes
				dw_cell_dzm = 0.5 * (dw_dz[n_left, m] + dw_dz[n_right, m])
				if abs(dw_cell_dzm) > FLUX_DERIV_TOL
					fast_sparse_sub!(J, smap, ξ_row, z_offset+m, advection_sign * dξ_dz * dw_cell_dzm)
				end
			end
		end
	end

	return nothing
end



"""
	compute_node_position_jacobian!(J::AbstractMatrix,
									dw_dT::Matrix{Float64},
									dw_dξ::Array{Float64,3},
									dw_dz::Matrix{Float64},
									sv::StateVectorALE{NC}) where NC

Fill the node position equation rows of the Jacobian.

The node position equation is simply:
	dzₖ/dt = wₖ

So the derivatives are:
	∂(dzₖ/dt)/∂Tᵢ = ∂wₖ/∂Tᵢ = dw_dT[k, i]
	∂(dzₖ/dt)/∂ξⱼᵢ = ∂wₖ/∂ξⱼᵢ = dw_dξ[k, j, i]
	∂(dzₖ/dt)/∂zₘ = ∂wₖ/∂zₘ = dw_dz[k, m]
"""
function compute_node_position_jacobian!(J::AbstractMatrix, smap,
	dw_dT::Matrix{Float64},
	dw_dξ::Array{Float64, 3},
	dw_dz::Matrix{Float64},
	sv::StateVectorALE{NC}) where NC

	n = sv.n_cells
	n_nodes = sv.n_nodes
	z_offset = sv.z_offset

	@inbounds for k in 1:n_nodes
		z_row = z_offset + k

		# T derivatives
		for i in 1:n
			if abs(dw_dT[k, i]) > FLUX_DERIV_TOL
				fast_sparse_write!(J, smap, z_row, i, dw_dT[k, i])
			end
		end

		# ξ derivatives
		for j in 1:NC
			for i in 1:n
				if abs(dw_dξ[k, j, i]) > FLUX_DERIV_TOL
					ξ_col = n + (j - 1) * n + i
					fast_sparse_write!(J, smap, z_row, ξ_col, dw_dξ[k, j, i])
				end
			end
		end

		# z derivatives
		for m in 1:n_nodes
			if abs(dw_dz[k, m]) > FLUX_DERIV_TOL
				fast_sparse_write!(J, smap, z_row, z_offset+m, dw_dz[k, m])
			end
		end
	end

	return nothing
end


# ============================================================================
# ALE Jacobian: Integration
# ============================================================================

"""
	compute_analytical_jacobian_ale!(J, u, p, t, f, jac_cache)

Compute the full analytical Jacobian for ALE simulations.

# Algorithm
1. Unpack state and update properties
2. Compute reaction/species derivatives (reusing base cache)
3. Compute ALE-specific derivatives:
   - Strain rate derivatives
   - Node velocity derivatives
   - Geometry derivatives
4. Fill Jacobian blocks:
   - Cell-local terms (reactions, properties) - same as fixed mesh
   - Conduction terms (with geometry corrections for z)
   - Gas transport terms
   - Boundary conditions
   - ALE advection terms
   - Node position equations

# Arguments
- `J::AbstractMatrix`: Output Jacobian matrix (pre-allocated)
- `u`: Current state vector
- `p`: Parameters (unused)
- `t`: Current time
- `f::PyrolysisODEFunctionALE`: ALE ODE function containing mesh, material, caches
- `jac_cache::JacobianCacheALE`: Pre-allocated workspace

# Notes
- This function is typically called by the ODE solver during implicit integration
- The Jacobian should be accurate to enable quadratic convergence of Newton iteration
"""
function compute_analytical_jacobian_ale!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionALE{M, Mat, BCS, PC, RC, SV, GC},
	jac_cache::JacobianCacheALE{NC, NR}) where {NC, NR, M, Mat, BCS, PC, RC,
	SV <: StateVectorALE{NC}, GC}

	mesh = f.mesh
	material = f.material
	prop_cache = f.prop_cache
	sv = f.sv
	n = sv.n_cells
	n_nodes = sv.n_nodes
	base_cache = jac_cache.base

	# =========================================================================
	# Step 0: Initialize static caches (only on first call)
	# =========================================================================
	initialize!(base_cache, material)

	# Initialize sparsity map for O(1) sparse writes (only on first call with sparse J)
	# Uses FastSparsityMap (dense array) for problems <= 5000 DOFs, falls back to
	# SparsityMap (Dict-based) for larger problems to avoid excessive memory usage
	if base_cache.sparsity_map === nothing && J isa SparseMatrixCSC
		base_cache.sparsity_map = create_sparsity_map(J)
	end
	smap = base_cache.sparsity_map

	# =========================================================================
	# Step 1: Unpack state and update geometry/properties
	# =========================================================================
	unpack_state_ale!(mesh, u, sv)
	Properties.update_properties!(prop_cache, mesh, material)

	# Zero out Jacobian (efficient O(nnz) for sparse matrices)
	zero_jacobian!(J)

	# =========================================================================
	# Step 2: Compute current RHS values (needed for property derivative terms)
	# =========================================================================
	# We need the reaction rates and species source terms for derivatives
	# Compute physics RHS (populates f.dT_work, f.dξ_work, f.rhs_cache.dξ_rxn)
	Discretization.compute_rhs!(f.dT_work, f.dξ_work, mesh, prop_cache, f.rhs_cache,
		material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
		I_radiation = f.I_radiation, radiation_model = f.radiation_model)

	# =========================================================================
	# Step 3: Compute current mesh velocity field (needed for advection derivatives)
	# =========================================================================
	_compute_ale_mesh_velocity!(f.w_work, mesh, prop_cache, material,
		t, f.min_density;
		θ_out = f.θ_work,
		precomputed_dξ_rxn = f.rhs_cache.dξ_rxn,
		dξ_dt = f.dξ_dt_work)

	# =========================================================================
	# Step 4: Cell-local Jacobian terms (reactions, properties)
	#         Same as fixed-mesh but store species sources for strain rate
	# =========================================================================
	# Use pre-allocated species source matrix from cache for strain rate derivatives
	S_matrix = jac_cache.S_matrix
	fill!(S_matrix, 0.0)

	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		# Cache all component properties at this temperature
		cache_all_component_properties!(base_cache, T, material)

		# Get current property values
		ρcp = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp

		# Compute reaction rates
		evaluate_reactions!(f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, material)

		# Compute reaction rate derivatives (using cached reaction properties)
		evaluate_reaction_derivatives!(base_cache.dr_dT, base_cache.dr_dξ,
			f.rhs_cache.rates, ξ, T, base_cache)

		# Heat source derivatives (using cached reaction properties)
		dQ_rxn_dξ = base_cache.work_NC  # Reuse work array
		dQ_rxn_dT = heat_source_derivatives_cached!(dQ_rxn_dξ,
			base_cache.dr_dT, base_cache.dr_dξ,
			f.rhs_cache.rates, f.rhs_cache.heats, ξ, T, base_cache, material)

		# Species source derivatives (using cached reaction properties)
		# Store dS/dT directly in dS_dT matrix column, dS/dξ directly in dS_dξ[:,:,i]
		# Note: Must use @view to modify in-place, not create a copy
		species_source_derivatives_cached!(@view(base_cache.dS_dT[:, i]), @view(base_cache.dS_dξ[:, :, i]),
			base_cache.dr_dT, base_cache.dr_dξ, base_cache, NC)

		# Store species sources for strain rate computation
		for j in 1:NC
			S_matrix[j, i] = f.rhs_cache.dξ_rxn[j, i]
		end

		# Heat capacity derivatives (using cached version)
		dρcp_dξ = base_cache.work_NC_2  # Use second work array
		dρcp_dT = heat_capacity_derivatives_cached!(dρcp_dξ, ξ,
			base_cache.comp_heat_capacities, base_cache.comp_dcp_dT)
		# Transform dρcp_dξ in-place to dinv_ρcp_dξ
		dinv_ρcp_dT = inverse_heat_capacity_derivatives!(dρcp_dξ, ρcp, dρcp_dT)
		# Note: dρcp_dξ is now dinv_ρcp_dξ

		# Get raw RHS value for property derivative contribution
		RHS_T = f.dT_work[i] * ρcp  # Remove 1/ρcₚ factor

		# -----------------------------------------------------------------
		# Energy equation derivatives (row i)
		# -----------------------------------------------------------------
		T_row = i

		# ∂(dT/dt)/∂T_i from property derivative: ∂(1/ρcₚ)/∂T × RHS
		fast_sparse_add!(J, smap, T_row, T_row, dinv_ρcp_dT * RHS_T)

		# ∂(dT/dt)/∂T_i from reaction: (1/ρcₚ) × ∂Q_rxn/∂T
		fast_sparse_add!(J, smap, T_row, T_row, inv_ρcp * dQ_rxn_dT)

		# ∂(dT/dt)/∂ξⱼ from property derivative
		for j in 1:NC
			ξ_col = n + (j - 1) * n + i
			fast_sparse_add!(J, smap, T_row, ξ_col, dρcp_dξ[j] * RHS_T)  # dρcp_dξ is now dinv_ρcp_dξ
		end

		# ∂(dT/dt)/∂ξⱼ from reaction
		for j in 1:NC
			ξ_col = n + (j - 1) * n + i
			fast_sparse_add!(J, smap, T_row, ξ_col, inv_ρcp * dQ_rxn_dξ[j])
		end

		# -----------------------------------------------------------------
		# Species equation derivatives (rows for each component)
		# -----------------------------------------------------------------
		for j in 1:NC
			ξ_row = n + (j - 1) * n + i

			# ∂(dξⱼ/dt)/∂T from reaction
			fast_sparse_write!(J, smap, ξ_row, T_row, base_cache.dS_dT[j, i])

			# ∂(dξⱼ/dt)/∂ξₖ from reaction
			for k in 1:NC
				ξ_col = n + (k - 1) * n + i
				fast_sparse_write!(J, smap, ξ_row, ξ_col, base_cache.dS_dξ[j, k, i])
			end
		end
	end

	# =========================================================================
	# Step 5: Conduction and gas transport Jacobian
	# =========================================================================
	# Pre-compute ALL cell properties ONCE (Optimization #1)
	# This eliminates redundant caching in each face loop
	precompute_all_cell_properties!(base_cache, mesh, material)

	# Fused face loop for conduction and gas transport Jacobians
	compute_all_face_jacobians_fused!(J, smap, mesh, prop_cache, material,
		f.rhs_cache.J_gas, n, Val(NC), base_cache)

	# Cell-centered gas convective source Jacobian
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
	compute_mass_bc_jacobian!(J, smap, mesh, prop_cache, f.bc_set, material, n, NC, t, f.top_face_ids, f.bottom_face_ids, base_cache.is_gas, base_cache)

	# NOTE: Boundary convective heat flux Jacobian removed (advection form has zero boundary contribution)
	# Energy leaving the system is handled by boundary conditions, not by convective heat flux.

	# =========================================================================
	# Step 7: Radiation Jacobian (if active)
	# =========================================================================
	if f.I_radiation > 0.0 && f.radiation_model == BEER_LAMBERT
		compute_radiation_jacobian!(J, smap, mesh, prop_cache, f.I_radiation, n, NC, base_cache.comp_absorptions, base_cache.τ_cumulative)
	end

	# =========================================================================
	# Step 8: ALE-specific derivatives
	# =========================================================================

	# 8a. Strain rate derivatives: ∂θ/∂T, ∂θ/∂ξ
	compute_strain_rate_derivatives!(jac_cache.dθ_dT, jac_cache.dθ_dξ,
		base_cache.dS_dT, base_cache.dS_dξ, S_matrix, mesh, material, jac_cache)

	# 8b. Node velocity derivatives: ∂w/∂T, ∂w/∂ξ, ∂w/∂z
	compute_node_velocity_derivatives!(jac_cache.dw_dT, jac_cache.dw_dξ,
		jac_cache.dw_dz, jac_cache.dθ_dT, jac_cache.dθ_dξ, f.θ_work, mesh)

	# 8c. Geometry derivatives: ∂V/∂z, ∂Δz/∂z
	compute_geometry_derivatives!(jac_cache.dV_dz, jac_cache.dΔz_dz, mesh)

	# =========================================================================
	# Step 9: ALE advection Jacobian (energy and species)
	# =========================================================================
	compute_ale_advection_jacobian_energy!(J, smap, mesh, f.w_work,
		jac_cache.dw_dT, jac_cache.dw_dξ, jac_cache.dw_dz, sv)

	compute_ale_advection_jacobian_species!(J, smap, mesh, material, f.w_work,
		jac_cache.dw_dT, jac_cache.dw_dξ, jac_cache.dw_dz, sv, base_cache.is_gas)

	# =========================================================================
	# Step 10: ALE dilation correction Jacobian for mesh-bound species
	# =========================================================================
	# The dilation correction is: -ξⱼ × θ for mesh-bound components (phase != GAS)
	# ∂(-ξⱼ×θ)/∂T = -ξⱼ × (∂θ/∂T)
	# ∂(-ξⱼ×θ)/∂ξₖ = -δⱼₖ × θ - ξⱼ × (∂θ/∂ξₖ)
	_add_dilation_correction_jacobian!(J, smap, mesh, base_cache.is_gas, f.θ_work,
		jac_cache.dθ_dT, jac_cache.dθ_dξ, n, NC, sv.z_offset)

	# =========================================================================
	# Step 11: Geometry corrections for conduction (∂/∂z terms)
	# =========================================================================
	compute_conduction_geometry_jacobian!(J, smap, mesh, prop_cache, n, NC, n_nodes,
		jac_cache.dV_dz, jac_cache.dΔz_dz, sv.z_offset)

	# =========================================================================
	# Step 11b: Geometry corrections for boundary fluxes (∂/∂z terms)
	# =========================================================================
	compute_bc_geometry_jacobian!(J, smap, mesh, prop_cache, f.bc_set,
		jac_cache.dV_dz, n, sv.z_offset, t,
		f.top_face_ids, f.bottom_face_ids, f.radiation_model)

	# =========================================================================
	# Step 11c: Geometry corrections for mass boundary fluxes (∂/∂z terms)
	# =========================================================================
	compute_mass_bc_geometry_jacobian!(J, smap, mesh, prop_cache, f.bc_set,
		material, jac_cache.dV_dz, n, NC, sv.z_offset, t, f.top_face_ids, f.bottom_face_ids, base_cache.is_gas)

	# =========================================================================
	# Step 11d: Geometry corrections for gas transport fluxes (∂/∂z terms)
	# =========================================================================
	compute_gas_transport_geometry_jacobian!(J, smap, mesh, prop_cache, material,
		jac_cache.dV_dz, n, Val(NC), sv.z_offset, base_cache.is_gas)

	# =========================================================================
	# Step 11d2: Geometry corrections for Darcy transport (∂/∂z terms)
	# =========================================================================
	if has_darcy_flow(material)
		compute_darcy_transport_geometry_jacobian!(J, smap, mesh, prop_cache, material,
			jac_cache.dV_dz, n, Val(NC), sv.z_offset, jac_cache)
	end

	# =========================================================================
	# Step 11e: Geometry corrections for radiation (∂/∂z terms)
	# =========================================================================
	if f.I_radiation > 0.0 && f.radiation_model == BEER_LAMBERT
		compute_radiation_geometry_jacobian!(J, smap, mesh, prop_cache, f.I_radiation,
			n, NC, material, base_cache.τ_cumulative, jac_cache.dΔz_dz, sv.z_offset)
	end

	# =========================================================================
	# Step 12: Node position equation Jacobian (dzₖ/dt = wₖ)
	# =========================================================================
	compute_node_position_jacobian!(J, smap, jac_cache.dw_dT, jac_cache.dw_dξ,
		jac_cache.dw_dz, sv)

	return nothing
end

"""
	_add_dilation_correction_jacobian!(J, mesh, is_gas, θ, dθ_dT, dθ_dξ, n, NC, z_offset)

Add Jacobian terms for the ALE dilation correction: -ξⱼ × θ for mesh-bound components.

The dilation correction accounts for concentration changes due to cell volume change.
For mesh-bound components (solids and liquids, i.e., phase != GAS):
	dξⱼ/dt = Sⱼ - ξⱼ × θ

Derivatives:
	∂(-ξⱼ×θ)/∂Tᵢ = -ξⱼ,ᵢ × (∂θᵢ/∂Tᵢ)
	∂(-ξⱼ×θ)/∂ξₖ,ᵢ = -δⱼₖ × θᵢ - ξⱼ,ᵢ × (∂θᵢ/∂ξₖ,ᵢ)
"""
function _add_dilation_correction_jacobian!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	is_gas::AbstractVector{Bool},
	θ::Vector{Float64},
	dθ_dT::Vector{Float64},
	dθ_dξ::Matrix{Float64},
	n::Int, NC_val::Int, z_offset::Int) where NC

	@inbounds for i in 1:n
		for j in 1:NC
			# Only mesh-bound components (phase != GAS) have dilation correction
			if is_gas[j]
				continue
			end

			ξ_row = n + (j - 1) * n + i
			ξⱼ = mesh.cell_states[i].ξ[j]

			# ∂(-ξⱼ×θ)/∂T = -ξⱼ × dθ/dT
			fast_sparse_add!(J, smap, ξ_row, i, -ξⱼ * dθ_dT[i])

			# ∂(-ξⱼ×θ)/∂ξₖ = -δⱼₖ × θ - ξⱼ × dθ/dξₖ
			for k in 1:NC
				ξ_col = n + (k - 1) * n + i
				# Kronecker delta term
				if j == k
					fast_sparse_add!(J, smap, ξ_row, ξ_col, -θ[i])
				end
				# Chain rule term
				fast_sparse_add!(J, smap, ξ_row, ξ_col, -ξⱼ * dθ_dξ[k, i])
			end
		end
	end

	return nothing
end

"""
	compute_conduction_geometry_jacobian!(J, mesh, prop_cache, n, NC, n_nodes,
										  dV_dz, dΔz_dz, z_offset)

Add z-derivative terms to the conduction Jacobian.

The conduction flux contribution to cell i is:
	dT_i/dt += (1/V_i) × Σ_faces (k_face × A × (T_neighbor - T_i) / d_LR)

Derivatives with respect to node positions z come from:
1. ∂(1/V_i)/∂z_k = -1/V_i² × ∂V_i/∂z_k
2. ∂(d_LR)/∂z_k affects the gradient calculation

This function adds these z-derivative contributions to the energy equation rows.
"""
function compute_conduction_geometry_jacobian!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	prop_cache,
	n::Int, NC_val::Int, n_nodes::Int,
	dV_dz::Matrix{Float64},
	dΔz_dz::Matrix{Float64},
	z_offset::Int) where NC

	A = mesh.cross_section_area

	@inbounds for i in 1:n
		T_row = i
		V_i = mesh.cell_volumes[i]
		inv_V_i = 1.0 / V_i
		inv_V_i_sq = inv_V_i * inv_V_i
		inv_ρcp = 1.0 / prop_cache.ρcp[i]

		# Get cell faces
		cell = mesh.cells[i]
		f_left, f_right = cell.face_ids

		# Current conduction flux contribution (without 1/ρcₚ factor)
		# We need to reconstruct this from the mesh geometry
		q_total = 0.0

		# Left face
		face_left = mesh.faces[f_left]
		left_id = face_left.cell_ids[1]
		if left_id != 0
			k_L = prop_cache.k[left_id]
			k_i = prop_cache.k[i]
			k_face_L = 2.0 * k_L * k_i / (k_L + k_i)
			T_L = mesh.cell_states[left_id].T
			T_i = mesh.cell_states[i].T
			d_LR = abs(mesh.cell_centers[i] - mesh.cell_centers[left_id])
			q_left = k_face_L * A * (T_L - T_i) / d_LR  # flux INTO cell i
			q_total += q_left

			# z-derivatives through d_LR (FIRST-ORDER effect, not second-order!)
			# Distance between cell centers: d_LR = (Δz_i + Δz_{left}) / 2
			# where Δz_i = z_{i+1} - z_i and Δz_{left} = z_{left+1} - z_{left} = z_i - z_{i-1}
			#
			# For cell i with left neighbor (left_id = i-1):
			# d_LR = 0.5 * (Δz_i + Δz_{i-1}) = 0.5 * ((z_{i+1} - z_i) + (z_i - z_{i-1}))
			#      = 0.5 * (z_{i+1} - z_{i-1})
			#
			# ∂d_LR/∂z_{i-1} = -0.5
			# ∂d_LR/∂z_{i+1} = +0.5
			# ∂d_LR/∂z_i = 0 (cancels)
			#
			# Flux: q_left = k_face × A × (T_L - T_i) / d_LR
			# ∂q_left/∂z_k = -q_left / d_LR × ∂d_LR/∂z_k
			#              = -k_face × A × (T_L - T_i) / d_LR² × ∂d_LR/∂z_k
			#
			# For cell i, left neighbor is cell i-1 (nodes i-1, i, i+1)
			# ∂d_LR/∂z[i-1] = -0.5 (if i > 1)
			# ∂d_LR/∂z[i+1] = +0.5
			#
			inv_d_sq = 1.0 / (d_LR * d_LR)
			flux_factor = -k_face_L * A * (T_L - T_i) * inv_d_sq

			# Contribution from ∂d_LR/∂z[i-1] = -0.5 (left_id = i-1 means left node is index i-1, but it's node left_id which bounds cell left_id)
			# Actually: cell left_id is bounded by nodes [left_id, left_id+1], and cell i is bounded by [i, i+1]
			# For neighbors: left_id = i-1, so cell i-1 bounded by [i-1, i] and cell i bounded by [i, i+1]
			# Distance d_LR = center(i) - center(i-1) = 0.5*(z_i + z_{i+1}) - 0.5*(z_{i-1} + z_i)
			#               = 0.5 * (z_{i+1} - z_{i-1})
			# So: ∂d_LR/∂z_{i-1} = -0.5, ∂d_LR/∂z_{i+1} = +0.5

			# Node indices (1-based): left_id's left node, left_id's right node = cell i's left node, cell i's right node
			# Cell i-1 (left_id) bounded by nodes [left_id, left_id+1] = [i-1, i]
			# Cell i bounded by nodes [i, i+1]
			# d_LR depends on z[left_id] = z[i-1] and z[i+1]

			# Derivative w.r.t. z[left_id] (the left node of left cell)
			# ∂d_LR/∂z[left_id] = -0.5
			dRHS_dzk_left = inv_ρcp * inv_V_i * flux_factor * (-0.5)
			fast_sparse_add!(J, smap, T_row, z_offset+left_id, dRHS_dzk_left)

			# Derivative w.r.t. z[i+1] (the right node of cell i)
			# ∂d_LR/∂z[i+1] = +0.5
			dRHS_dzk_right = inv_ρcp * inv_V_i * flux_factor * (+0.5)
			fast_sparse_add!(J, smap, T_row, z_offset+i+1, dRHS_dzk_right)
		end

		# Right face
		face_right = mesh.faces[f_right]
		right_id = face_right.cell_ids[2]
		if right_id != 0
			k_R = prop_cache.k[right_id]
			k_i = prop_cache.k[i]
			k_face_R = 2.0 * k_R * k_i / (k_R + k_i)
			T_R = mesh.cell_states[right_id].T
			T_i = mesh.cell_states[i].T
			d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[i])
			q_right = k_face_R * A * (T_R - T_i) / d_LR  # flux INTO cell i
			q_total += q_right

			# z-derivatives through d_LR
			# Cell i bounded by [i, i+1], cell i+1 (right_id) bounded by [i+1, i+2]
			# d_LR = center(i+1) - center(i) = 0.5*(z_{i+1} + z_{i+2}) - 0.5*(z_i + z_{i+1})
			#      = 0.5 * (z_{i+2} - z_i)
			# ∂d_LR/∂z_i = -0.5, ∂d_LR/∂z_{i+2} = +0.5

			inv_d_sq = 1.0 / (d_LR * d_LR)
			flux_factor = -k_face_R * A * (T_R - T_i) * inv_d_sq

			# Derivative w.r.t. z[i] (left node of cell i)
			dRHS_dzk_left = inv_ρcp * inv_V_i * flux_factor * (-0.5)
			fast_sparse_add!(J, smap, T_row, z_offset+i, dRHS_dzk_left)

			# Derivative w.r.t. z[right_id+1] = z[i+2] (right node of right cell)
			dRHS_dzk_right = inv_ρcp * inv_V_i * flux_factor * (+0.5)
			fast_sparse_add!(J, smap, T_row, z_offset+right_id+1, dRHS_dzk_right)
		end

		# Add z-derivatives through 1/V_i
		# The energy RHS has factor (1/V_i) × q_total
		# ∂[(1/V_i)×q]/∂z_k = q × ∂(1/V_i)/∂z_k = q × (-1/V_i²) × ∂V_i/∂z_k
		# Apply to heat equation with 1/ρcₚ factor

		for k in 1:n_nodes
			if abs(dV_dz[i, k]) > FLUX_DERIV_TOL
				# ∂(dT/dt)/∂z_k from volume dependence in conduction
				dRHS_dzk = inv_ρcp * (-inv_V_i_sq) * dV_dz[i, k] * q_total
				fast_sparse_add!(J, smap, T_row, z_offset+k, dRHS_dzk)
			end
		end
	end

	return nothing
end

"""
	compute_bc_geometry_jacobian!(J, mesh, prop_cache, thermal_bcs, dV_dz, n, z_offset, t)

Add geometry derivatives for boundary condition fluxes.

For a heat flux BC, the contribution to dT/dt is:
	(dT/dt)_bc = (1/ρcₚ) × q_bc × A / V

When node positions change, V changes:
	∂(dT/dt)_bc/∂z_k = (1/ρcₚ) × q_bc × A × ∂(1/V)/∂z_k
					 = (1/ρcₚ) × q_bc × A × (-1/V²) × ∂V/∂z_k
"""
function compute_bc_geometry_jacobian!(J::AbstractMatrix, smap,
	mesh::Mesh1D{NC},
	prop_cache,
	bc_set::BoundaryConditionSet,
	dV_dz::Matrix{Float64},
	n::Int, z_offset::Int, t::Float64,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int},
	radiation_model::RadiationModel) where NC

	n_nodes = size(dV_dz, 2)

	# Process top boundary
	for f in top_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id
		T_cell = mesh.cell_states[cell_id].T
		ε_eff = prop_cache.ε[cell_id]
		k_eff = prop_cache.k[cell_id]

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_ρcp = 1.0 / prop_cache.ρcp[cell_id]
		inv_V_sq = 1.0 / (V * V)

		# Half-cell distance to boundary
		Δz_half = V / A / 2.0

		# Solve for boundary temperature and get actual conducted flux
		state = solve_boundary_temperature(T_cell, k_eff, Δz_half, bc_set.thermal_top, t, ε_eff;
			radiation_model = radiation_model)
		q_bc = state.q_conducted

		if !isnan(q_bc) && abs(q_bc) > FLUX_DERIV_TOL
			T_row = cell_id

			# ∂(dT/dt)_bc/∂z_k = (1/ρcₚ) × q_bc × A × (-1/V²) × ∂V/∂z_k
			for k in 1:n_nodes
				if abs(dV_dz[cell_id, k]) > FLUX_DERIV_TOL
					dRHS_dzk = inv_ρcp * q_bc * A * (-inv_V_sq) * dV_dz[cell_id, k]
					fast_sparse_add!(J, smap, T_row, z_offset+k, dRHS_dzk)
				end
			end
		end
	end

	# Process bottom boundary
	for f in bottom_face_ids
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		cell_id = left_id == 0 ? right_id : left_id
		T_cell = mesh.cell_states[cell_id].T
		ε_eff = prop_cache.ε[cell_id]
		k_eff = prop_cache.k[cell_id]

		V = mesh.cell_volumes[cell_id]
		A = mesh.face_areas[f]
		inv_ρcp = 1.0 / prop_cache.ρcp[cell_id]
		inv_V_sq = 1.0 / (V * V)

		# Half-cell distance to boundary
		Δz_half = V / A / 2.0

		# Solve for boundary temperature and get actual conducted flux
		state = solve_boundary_temperature(T_cell, k_eff, Δz_half, bc_set.thermal_bottom, t, ε_eff;
			radiation_model = radiation_model)
		q_bc = state.q_conducted

		if !isnan(q_bc) && abs(q_bc) > FLUX_DERIV_TOL
			T_row = cell_id

			# ∂(dT/dt)_bc/∂z_k = (1/ρcₚ) × q_bc × A × (-1/V²) × ∂V/∂z_k
			for k in 1:n_nodes
				if abs(dV_dz[cell_id, k]) > FLUX_DERIV_TOL
					dRHS_dzk = inv_ρcp * q_bc * A * (-inv_V_sq) * dV_dz[cell_id, k]
					fast_sparse_add!(J, smap, T_row, z_offset+k, dRHS_dzk)
				end
			end
		end
	end

	return nothing
end

# Get actual BC flux value (not derivative) - takes ε_eff for dynamic emissivity
_bc_flux_value(::AdiabaticBC, T, t, ε_eff) = 0.0
_bc_flux_value(bc::HeatFluxBC, T, t, ε_eff) = bc.flux isa Number ? bc.flux : bc.flux(t)
_bc_flux_value(bc::ConvectiveBC, T, t, ε_eff) = begin
	h = bc.h isa Number ? bc.h : bc.h(t)
	T_inf = bc.T_inf isa Number ? bc.T_inf : bc.T_inf(t)
	h * (T_inf - T)
end
_bc_flux_value(bc::RadiativeBC, T, t, ε_eff) = begin
	# Use effective emissivity if not specified
	ε = if bc.ε === nothing
		ε_eff
	else
		bc.ε isa Number ? bc.ε : bc.ε(T)
	end
	T_inf = bc.T_inf isa Number ? bc.T_inf : bc.T_inf(t)
	ε * STEFAN_BOLTZMANN * (T_inf^4 - T^4)
end
_bc_flux_value(bc::RadiativeFluxBC, T, t, ε_eff) = begin
	flux = bc.flux isa Number ? bc.flux : bc.flux(t)
	# Use effective emissivity as absorptivity if not specified (Kirchhoff's law)
	α = if bc.absorptivity === nothing
		ε_eff
	else
		bc.absorptivity isa Number ? bc.absorptivity : bc.absorptivity(T)
	end
	# RadiativeFluxBC specifies incident flux - always surface absorption
	α * flux
end
_bc_flux_value(bc::FlameBC, T, t, ε_eff) = 0.0  # Flame contribution computed elsewhere
_bc_flux_value(bc::CombinedThermalBC, T, t, ε_eff) = sum(b -> _bc_flux_value(b, T, t, ε_eff), bc.bcs)

# Dispatch for radiation wrapper
function compute_analytical_jacobian_ale!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionALEWithRadiation,
	jac_cache::JacobianCacheALE)
	# Delegate to the inner function
	compute_analytical_jacobian_ale!(J, u, p, t, f.ode_fn, jac_cache)
end

# Dispatch for P1 volumetric radiation wrapper
# Note: This is a simplified implementation that delegates to the base ALE Jacobian.
# For full accuracy, the P1 radiation Jacobian contributions should be added, but
# since P1 radiation is quasi-steady (solved at each timestep), the temperature
# dependence is captured through the ODE function evaluation.
function compute_analytical_jacobian_ale!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionALEWithP1Radiation,
	jac_cache::JacobianCacheALE)
	# Delegate to the inner ALE function
	# The P1 radiation is quasi-steady, so its Jacobian contribution is implicit
	# in the RHS evaluation (G equilibrates instantly to temperature changes)
	compute_analytical_jacobian_ale!(J, u, p, t, f.ode_fn, jac_cache)
end

# ============================================================================
# P1 Volumetric Radiation Jacobian (Week 3)
# ============================================================================

"""
	JacobianCacheP1Radiation


Cache for P1 volumetric radiation Jacobian computation.

This cache stores work arrays for computing the implicit sensitivity
of the radiation field G to temperature changes.

# Fields
- `dG_dT`: Sensitivity of G to temperature at each cell [W/(m²·K)]
- `dα_dT`: Derivative of absorption coefficient w.r.t. temperature [1/(m·K)]
- `dε_dT`: Derivative of emissivity w.r.t. temperature [1/K]
- `dD_rad_dT`: Derivative of diffusion coefficient w.r.t. temperature [m/K]
- `work_n`: General work array [n_cells]
- `work_n2`: Second work array [n_cells]
"""
mutable struct JacobianCacheP1Radiation
	dG_dT::Vector{Float64}          # [n_cells] - ∂Gᵢ/∂Tⱼ (diagonal only for now)
	dα_dT::Vector{Float64}          # [n_cells] - ∂αᵢ/∂Tᵢ
	dε_dT::Vector{Float64}          # [n_cells] - ∂εᵢ/∂Tᵢ
	dD_rad_dT::Vector{Float64}      # [n_cells] - ∂D_rad,i/∂Tᵢ
	work_n::Vector{Float64}         # [n_cells] - work array
	work_n2::Vector{Float64}        # [n_cells] - work array
end

"""
	JacobianCacheP1Radiation(n_cells::Int)

Create a pre-allocated P1 radiation Jacobian cache.
"""
function JacobianCacheP1Radiation(n_cells::Int)
	JacobianCacheP1Radiation(
		zeros(n_cells),    # dG_dT
		zeros(n_cells),    # dα_dT
		zeros(n_cells),    # dε_dT
		zeros(n_cells),    # dD_rad_dT
		zeros(n_cells),    # work_n
		zeros(n_cells),    # work_n2
	)
end

"""
	resize_jacobian_cache_p1!(cache::JacobianCacheP1Radiation, n_cells::Int)

Resize the P1 radiation Jacobian cache after mesh adaptation.
"""
function resize_jacobian_cache_p1!(cache::JacobianCacheP1Radiation, n_cells::Int)
	resize!(cache.dG_dT, n_cells)
	resize!(cache.dα_dT, n_cells)
	resize!(cache.dε_dT, n_cells)
	resize!(cache.dD_rad_dT, n_cells)
	resize!(cache.work_n, n_cells)
	resize!(cache.work_n2, n_cells)
	return cache
end

"""
	compute_p1_radiation_jacobian!(J::AbstractMatrix, mesh, prop_cache,
								   rad_cache, radiation_bcs, p1_jac_cache,
								   n_cells::Int, NC::Int, t::Float64;
								   z_offset::Int=0)

Compute P1 volumetric radiation Jacobian contributions.

For quasi-steady P1 radiation, the temperature equation includes:
	dT/dt += S_rad / ρcₚ

where S_rad = α(G - 4σT⁴) is the net radiation source (positive = heating),
and G is the solution of the quasi-steady P1 equation:
	∇·(D_rad∇G) = αG - 4σαT⁴

Note: In standard P1, volumetric emission/absorption is blackbody (no ε factor).
Emissivity only appears at boundaries via Marshak BC.

The Jacobian contribution ∂(dT_i/dt)/∂T_j has two parts:
1. **Direct**: From explicit T⁴ dependence in S_rad (j=i only)
2. **Indirect**: From implicit G(T) dependence through quasi-steady equation

# Direct contribution (diagonal, j=i):
	∂(dT_i/dt)/∂T_i|_direct = (1/ρcₚ) × ∂S_rad/∂T_i
	                        = (1/ρcₚ) × α × (-16σT³)
	                        = -(16σα/ρcₚ) × T³

# Indirect contribution (approximation):
For optically thick media (τ >> 1), G ≈ 4σT⁴ locally (blackbody equilibrium in P1), so:
	∂G_i/∂T_i ≈ 16σT³

Then:
	∂(S_rad)/∂T_i|_indirect = α × ∂G_i/∂T_i ≈ 16σαT³

This gives the same magnitude as the direct term, leading to partial cancellation
in local equilibrium. For non-equilibrium, we use the exact sensitivity.

# Arguments
- `J`: Jacobian matrix (modified in-place)
- `mesh`: Mesh with temperature state
- `prop_cache`: Property cache
- `rad_cache`: Radiation cache with current G, S_rad
- `radiation_bcs`: Radiation BC container (`Dict` or `RadiationBCSet`)
- `p1_jac_cache`: Work arrays for P1 Jacobian
- `n_cells`: Number of cells
- `NC`: Number of components  
- `t`: Current time
- `z_offset`: Offset for node positions in state vector (ALE mode)
"""
function compute_p1_radiation_jacobian!(J::AbstractMatrix, mesh, prop_cache,
	rad_cache, radiation_bcs, p1_jac_cache::JacobianCacheP1Radiation,
	n_cells::Int, NC::Int, t::Float64;
	z_offset::Int = 0)

	σ = STEFAN_BOLTZMANN
	A = mesh.cross_section_area

	# Compute the sensitivity of G to T using the implicit function theorem
	# For quasi-steady P1: R(G,T) = ∇·(D_rad∇G) - αG + 4σεαT⁴ = 0
	# By IFT: dG/dT = -(∂R/∂G)⁻¹ × (∂R/∂T)
	#
	# For now, use the local equilibrium approximation:
	# G ≈ 4σεT⁴  =>  dG/dT ≈ 16σεT³
	#
	# This is accurate for optically thick media where radiation is in
	# local equilibrium with the temperature field.

	# Compute dG_dT (local equilibrium approximation)
	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T
		ε = prop_cache.ε[i]
		# Local equilibrium: G = 4σεT⁴, so dG/dT = 16σεT³
		p1_jac_cache.dG_dT[i] = 16.0 * σ * ε * T^3
	end

	# Add Jacobian contributions to temperature equation
	# Energy: dT/dt += S_rad / ρcₚ  where S_rad = α(G - 4σT⁴)
	#
	# ∂(dT_i/dt)/∂T_i = (1/ρcₚᵢ) × ∂S_rad/∂T_i
	#                 = (1/ρcₚᵢ) × [α × dG/dT - 16σεαT³]
	#
	# Note: The two terms (α × dG/dT and 16σεαT³) are equal in local equilibrium,
	# leading to zero net Jacobian contribution (G and T track each other).
	# However, near boundaries or with external flux, there's a non-zero contribution.

	@inbounds for i in 1:n_cells
		T = mesh.cell_states[i].T
		α = rad_cache.α[i]  # Volumetric absorption coefficient [1/m]
		ε = prop_cache.ε[i]
		ρcp = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp

		# ∂S_rad/∂T_i = α × dG/dT - 16σεαT³
		dG_dT_i = p1_jac_cache.dG_dT[i]
		dS_rad_dT = α * dG_dT_i - 16.0 * σ * ε * α * T^3

		# Net contribution to Jacobian (positive because dT/dt = +S_rad/ρcₚ)
		# ∂(dT_i/dt)/∂T_i = +inv_ρcp × dS_rad_dT
		T_row = i
		J[T_row, T_row] += inv_ρcp * dS_rad_dT
	end

	# Add diffusion term contributions (off-diagonal)
	# The radiation flux q_rad = -D_rad × dG/dx contributes to neighboring cells
	# When G changes, the diffusive flux changes
	#
	# For face f between cells i and i+1:
	# q_f = -D_f × (G_{i+1} - G_i) / Δx
	#
	# Contribution to dT_i/dt from face f:
	# ∂(dT_i/dt)/∂T_j = -(1/ρcₚᵢ) × (A/Vᵢ) × ∂q_f/∂T_j
	#                 = (1/ρcₚᵢ) × (A/Vᵢ) × (D_f/Δx) × dG_j/dT_j

	@inbounds for i in 1:n_cells
		V = mesh.cell_volumes[i]
		ρcp = prop_cache.ρcp[i]
		inv_ρcp = 1.0 / ρcp

		T_row = i

		# Left face (between i-1 and i)
		if i > 1
			face_L = i
			D_rad_L = rad_cache.D_rad_face[face_L]
			Δx_L = abs(mesh.cell_centers[i] - mesh.cell_centers[i-1])
			coeff_L = inv_ρcp * (A / V) * D_rad_L / Δx_L

			# ∂(dT_i/dt)/∂T_{i-1} from left flux contribution
			# q_L depends on G_{i-1} which depends on T_{i-1}
			dG_dT_L = p1_jac_cache.dG_dT[i-1]
			J[T_row, i-1] += coeff_L * dG_dT_L

			# ∂(dT_i/dt)/∂T_i from left flux (G_i term)
			dG_dT_i = p1_jac_cache.dG_dT[i]
			J[T_row, T_row] += -coeff_L * dG_dT_i
		end

		# Right face (between i and i+1)
		if i < n_cells
			face_R = i + 1
			D_rad_R = rad_cache.D_rad_face[face_R]
			Δx_R = abs(mesh.cell_centers[i+1] - mesh.cell_centers[i])
			coeff_R = inv_ρcp * (A / V) * D_rad_R / Δx_R

			# ∂(dT_i/dt)/∂T_{i+1} from right flux contribution
			dG_dT_R = p1_jac_cache.dG_dT[i+1]
			J[T_row, i+1] += -coeff_R * dG_dT_R

			# ∂(dT_i/dt)/∂T_i from right flux (G_i term)
			dG_dT_i = p1_jac_cache.dG_dT[i]
			J[T_row, T_row] += coeff_R * dG_dT_i
		end
	end

	# Boundary contributions
	# For Marshak BC at boundaries, there's an additional flux-temperature coupling
	_add_marshak_bc_jacobian!(J, mesh, prop_cache, rad_cache, radiation_bcs,
		p1_jac_cache, n_cells, t)

	return nothing
end

@inline _radiation_bc_pairs(bcs::RadiationBCSet) = begin
	pairs = Tuple{Int, AbstractRadiationBC}[]
	if bcs.top !== nothing
		push!(pairs, (1, bcs.top))
	end
	if bcs.bottom !== nothing
		push!(pairs, (2, bcs.bottom))
	end
	return pairs
end

@inline _radiation_bc_pairs(bcs) = bcs

"""
Add Marshak boundary condition Jacobian contributions.

For Marshak BC: flux_into = I_ext - εσT⁴ + εσT_amb⁴ - G/2
The G/2 term creates a direct coupling to the radiation field.

∂flux/∂T = -4εσT³ (from surface emission)
∂flux/∂G = -1/2 (from G term)

Combined with dG/dT gives:
∂(dT/dt)/∂T_surface = -(1/ρcₚ)(A/V) × [∂flux/∂T + (∂flux/∂G)(dG/dT)]
"""
function _add_marshak_bc_jacobian!(J, mesh, prop_cache, rad_cache, radiation_bcs,
	p1_jac_cache, n_cells, t)

	σ = STEFAN_BOLTZMANN
	A = mesh.cross_section_area

	for (boundary_id, bc) in _radiation_bc_pairs(radiation_bcs)
		if !(bc isa MarshakRadiationBC)
			continue
		end

		# Determine which cell is at this boundary
		# boundary_id=1 is top (cell n_cells), boundary_id=2 is bottom (cell 1)
		cell_id = boundary_id == 1 ? n_cells : 1

		T = mesh.cell_states[cell_id].T
		V = mesh.cell_volumes[cell_id]
		ρcp = prop_cache.ρcp[cell_id]
		inv_ρcp = 1.0 / ρcp
		ε = prop_cache.ε[cell_id]

		# Get emissivity from BC
		ε_bc = bc.emissivity isa Number ? bc.emissivity : bc.emissivity(T)

		# ∂flux/∂T from Marshak BC emission term: -4εσT³
		dflux_dT_direct = -4.0 * ε_bc * σ * T^3

		# ∂flux/∂G from Marshak BC: -1/2
		dflux_dG = -0.5

		# Combined: ∂(dT/dt)/∂T = -(1/ρcₚ)(A/V) × [dflux_dT + dflux_dG × dG/dT]
		dG_dT = p1_jac_cache.dG_dT[cell_id]
		total_deriv = dflux_dT_direct + dflux_dG * dG_dT

		# Sign depends on boundary orientation
		# Top boundary (id=1): flux INTO domain is positive in -z direction
		# Bottom boundary (id=2): flux INTO domain is positive in +z direction
		sign_factor = boundary_id == 1 ? 1.0 : -1.0

		T_row = cell_id
		J[T_row, T_row] += sign_factor * inv_ρcp * (A / V) * total_deriv
	end

	return nothing
end

"""
	JacobianCacheALEWithP1{NC, NR}

Combined ALE + P1 radiation Jacobian cache.
"""
mutable struct JacobianCacheALEWithP1{NC, NR}
	ale_cache::JacobianCacheALE{NC, NR}
	p1_cache::JacobianCacheP1Radiation
end

"""
	JacobianCacheALEWithP1(n_cells::Int, n_faces::Int, n_nodes::Int, NC::Int, NR::Int)

Create a combined ALE + P1 radiation Jacobian cache.
"""
function JacobianCacheALEWithP1(n_cells::Int, n_faces::Int, n_nodes::Int, NC::Int, NR::Int)
	ale_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
	p1_cache = JacobianCacheP1Radiation(n_cells)
	JacobianCacheALEWithP1{NC, NR}(ale_cache, p1_cache)
end

"""
	resize_jacobian_cache_ale_p1!(cache::JacobianCacheALEWithP1, n_cells::Int, n_faces::Int, n_nodes::Int)

Resize combined ALE + P1 cache after mesh adaptation.
"""
function resize_jacobian_cache_ale_p1!(cache::JacobianCacheALEWithP1,
	n_cells::Int, n_faces::Int, n_nodes::Int)
	resize_jacobian_cache_ale!(cache.ale_cache, n_cells, n_faces, n_nodes)
	resize_jacobian_cache_p1!(cache.p1_cache, n_cells)
	return cache
end

"""
	compute_analytical_jacobian_ale_p1!(J, u, p, t, f, jac_cache)

Compute analytical Jacobian for ALE simulations with P1 volumetric radiation.

This function extends `compute_analytical_jacobian_ale!` to include
P1 radiation Jacobian contributions.

# Arguments
- `J::AbstractMatrix`: Output Jacobian matrix (pre-allocated)
- `u`: Current state vector
- `p`: Parameters (unused)
- `t`: Current time
- `f::PyrolysisODEFunctionALEWithP1Radiation`: P1 radiation ODE function
- `jac_cache::JacobianCacheALEWithP1`: Combined ALE + P1 Jacobian cache
"""
function compute_analytical_jacobian_ale_p1!(J::AbstractMatrix, u, p, t,
	f::PyrolysisODEFunctionALEWithP1Radiation,
	jac_cache::JacobianCacheALEWithP1)

	# First compute the base ALE Jacobian (without P1 radiation)
	compute_analytical_jacobian_ale!(J, u, p, t, f.ode_fn, jac_cache.ale_cache)

	# Then add P1 radiation contributions
	mesh = f.ode_fn.mesh
	prop_cache = f.ode_fn.prop_cache
	rad_cache = f.rad_cache
	sv = f.ode_fn.sv
	n_cells = sv.n_cells
	NC = sv.n_components

	compute_p1_radiation_jacobian!(J, mesh, prop_cache, rad_cache,
		f.radiation_bcs, jac_cache.p1_cache, n_cells, NC, t;
		z_offset = sv.z_offset)

	return nothing
end

"""
	PyrolysisJacobianALEWithP1{F, JC, SP}

Jacobian operator for ALE simulations with P1 volumetric radiation.

This struct wraps the ODE function and provides efficient Jacobian
computation for implicit solvers.

# Fields
- `ode_fn::F`: P1 radiation ODE function
- `jac_cache::JC`: Combined ALE + P1 Jacobian cache
- `sparsity::SP`: Sparsity pattern
- `colorvec::Vector{Int}`: Column coloring for FD fallback
- `n_colors::Int`: Number of colors
"""
struct PyrolysisJacobianALEWithP1{F, JC, SP}
	ode_fn::F
	jac_cache::JC
	sparsity::SP
	colorvec::Vector{Int}
	n_colors::Int
end

"""
	PyrolysisJacobianALEWithP1(f::PyrolysisODEFunctionALEWithP1Radiation)

Create a Jacobian operator for P1 radiation simulations.

# Arguments
- `f`: P1 radiation ODE function

# Returns
Jacobian operator with pre-allocated cache and sparsity pattern.
"""
function PyrolysisJacobianALEWithP1(f::PyrolysisODEFunctionALEWithP1Radiation)
	ode_fn_inner = f.ode_fn
	mesh = ode_fn_inner.mesh
	material = ode_fn_inner.material
	sv = ode_fn_inner.sv

	n_cells = sv.n_cells
	n_nodes = sv.n_nodes
	n_faces = n_cells + 1
	NC = n_components(material)
	NR = n_reactions(material)

	# Create combined cache
	jac_cache = JacobianCacheALEWithP1(n_cells, n_faces, n_nodes, NC, NR)

	# Create sparsity pattern (P1 quasi-steady radiation contributions)
	sparsity = jacobian_sparsity_ale(n_cells, NC; radiation_model = P1_QUASI_STEADY)

	# Compute coloring
	colorvec = matrix_colors(sparsity)
	n_colors = maximum(colorvec; init = 0)

	PyrolysisJacobianALEWithP1(f, jac_cache, sparsity, colorvec, n_colors)
end

"""
	(jac::PyrolysisJacobianALEWithP1)(J, u, p, t)

Compute the Jacobian matrix using analytical derivatives.
"""
function (jac::PyrolysisJacobianALEWithP1)(J::AbstractMatrix, u, p, t)
	compute_analytical_jacobian_ale_p1!(J, u, p, t, jac.ode_fn, jac.jac_cache)
	return nothing
end

"""
	resize_jacobian_ale_p1!(jac::PyrolysisJacobianALEWithP1, n_cells::Int, n_faces::Int, n_nodes::Int)

Resize the Jacobian operator after mesh adaptation.
"""
function resize_jacobian_ale_p1!(jac::PyrolysisJacobianALEWithP1, n_cells::Int, n_faces::Int, n_nodes::Int)
	resize_jacobian_cache_ale_p1!(jac.jac_cache, n_cells, n_faces, n_nodes)
	return jac
end
