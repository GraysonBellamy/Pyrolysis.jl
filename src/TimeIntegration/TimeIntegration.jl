"""
	TimeIntegration

Temporal integration for the Pyrolysis.jl package.

This module provides:
- State vector management for packing/unpacking mesh state to ODE vectors
- ODE RHS assembly interfacing with DifferentialEquations.jl
- Jacobian computation (sparse FD, dense FD, and AD options)
- CFL-based timestep control for stability
- ALE ODE function with integrated mesh motion

Reference: Design Document Sections 7.1, 7.2, 8.1, 8.3, 8.4
"""
module TimeIntegration

using ..Core
using ..Core: n_components, n_reactions, ImpermeableBC, AbstractMeshVelocityMethod,
	GAS, SOLID
using ..Geometry
using ..Properties
using ..Physics
using ..Discretization
using ..ALE
using OrdinaryDiffEq
using ForwardDiff
using SparseArrays

# Include subfiles
include("state_vector.jl")
include("ode_function.jl")
include("timestep_control.jl")
include("jacobian.jl")
include("analytical_jacobian.jl")
include("reduced_jacobian.jl")

# ============================================================================
# Exports
# ============================================================================

# State vector management
export StateVector
export state_vector_length
export create_state_vector
export pack_state!, unpack_state!
export get_temperature, get_concentration, get_concentrations
export set_temperature!, set_concentration!
export temperature_view, concentration_view
export validate_state_vector, copy_state!

# ALE state vector management
export StateVectorALE
export state_vector_length_ale
export create_state_vector_ale
export pack_state_ale!, unpack_state_ale!
export get_temperature_ale, get_concentration_ale, get_concentrations_ale
export get_node_position, get_node_positions
export set_temperature_ale!, set_concentration_ale!
export set_node_position!, set_node_velocity!
export temperature_view_ale, concentration_view_ale, node_positions_view
export validate_state_vector_ale
export get_thickness_from_state
export compute_volumes_from_state, compute_cell_centers_from_state
# State vector resizing for cell merging
export remove_surface_cell_from_state_vector, remove_cell_from_state_vector
export create_merged_state_vector, map_cell_indices_after_removal

# P1 Radiation state vector management
export StateVectorALERadiation
export state_vector_length_ale_radiation
export create_state_vector_ale_radiation
export pack_state_ale_radiation!, unpack_state_ale_radiation!
export get_temperature_ale_rad, get_concentration_ale_rad
export get_radiation_intensity, get_node_position_rad
export set_temperature_ale_rad!, set_concentration_ale_rad!
export set_radiation_intensity!, set_node_position_rad!
export temperature_view_ale_rad, radiation_view
export concentration_view_ale_rad, node_positions_view_rad
export validate_state_vector_ale_radiation

# ODE function
export PyrolysisODEFunction
export PyrolysisODEFunctionWithRadiation
export resize_ode_function!
export get_mesh, get_material, get_state_vector_descriptor
export set_radiation!

# ALE ODE function (integrated mesh motion)
export PyrolysisODEFunctionALE
export PyrolysisODEFunctionALEWithRadiation
export get_mesh_velocity
export resize_ode_function_ale!

# P1 Radiation ODE function (volumetric radiation)
export PyrolysisODEFunctionALEWithP1Radiation
export get_radiation_intensity

# Timestep control
export compute_timestep
export compute_reaction_timestep
export compute_stable_timestep
export TimestepController
export adjust_timestep!, reset_timestep!

# Jacobian
export jacobian_sparsity
export jacobian_sparsity_ale
export JacobianOperator
export matrix_colors
export optimal_fd_step, FD_RELATIVE_STEP, FD_ABSOLUTE_MIN
export compute_jacobian_fd!
export compute_jacobian_fd_sparse!
export PyrolysisJacobianAD
export PyrolysisJacobianAD_ALE
export make_ode_function_with_jac

# Analytical Jacobian
export JacobianCache
export resize_jacobian_cache!
export compute_analytical_jacobian!
export evaluate_reaction_derivatives!
export heat_capacity_derivatives
export conductivity_derivatives
export harmonic_mean_derivatives
export zero_jacobian!

# ALE Analytical Jacobian infrastructure
export JacobianCacheALE
export resize_jacobian_cache_ale!

# ALE Analytical Jacobian core derivatives
export compute_strain_rate_derivatives!
export compute_node_velocity_derivatives!
export compute_geometry_derivatives!
export compute_ale_advection_jacobian_energy!
export compute_ale_advection_jacobian_species!
export compute_node_position_jacobian!

# ALE Analytical Jacobian integration
export compute_analytical_jacobian_ale!
export compute_conduction_geometry_jacobian!

# Reduced ALE Jacobian (omits dense mesh velocity coupling for O(N) factorization)
export jacobian_sparsity_ale_reduced
export compute_analytical_jacobian_ale_reduced!

# Darcy transport Jacobian
export compute_darcy_transport_jacobian!
export compute_darcy_transport_geometry_jacobian!

# P1 Volumetric Radiation Jacobian (Week 3)
export JacobianCacheP1Radiation
export resize_jacobian_cache_p1!
export compute_p1_radiation_jacobian!
export JacobianCacheALEWithP1
export resize_jacobian_cache_ale_p1!
export compute_analytical_jacobian_ale_p1!
export PyrolysisJacobianALEWithP1
export resize_jacobian_ale_p1!

end