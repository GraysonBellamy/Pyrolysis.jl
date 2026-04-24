"""
	Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux using an Arbitrary Lagrangian-Eulerian (ALE) adaptive
meshing framework.

# Main Features
- Multi-component materials with temperature-dependent properties
- Arrhenius-type decomposition kinetics (single and multi-step)
- Gas generation and transport through porous char
- In-depth radiation absorption (Beer-Lambert)
- Material volume changes (shrinkage, swelling, intumescence)
- Surface regression and internal interface motion
- Adaptive mesh refinement (h- and r-adaptivity)

# Quick Start
```julia
using Pyrolysis

# Define a simple material
pmma = Material(
	name = :PMMA,
	components = (
		SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86),
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
	),
	reactions = (
		Reaction(:decomposition, :virgin => :MMA, A=8.5e12, E=188e3, h=-870e3),
	)
)

# Create mesh and solve
mesh = create_uniform_mesh(0.01, 100, 2)  # 1 cm, 100 cells, 2 components
problem = PyrolysisProblem(mesh, pmma, bcs, T_initial, tspan)
solution = solve(problem)
```
"""
module Pyrolysis

using LinearAlgebra
using SparseArrays
using StaticArrays
using ForwardDiff
using LinearSolve
using OrdinaryDiffEq


# ============================================================================
# Submodule includes
# ============================================================================

# Core types and constants (no dependencies on other submodules)
include("Core/Core.jl")
using .Core

# Geometry and mesh infrastructure
include("Geometry/Geometry.jl")
using .Geometry

# Material property evaluation
include("Properties/Properties.jl")
using .Properties

# Physical models
include("Physics/Physics.jl")
using .Physics

# Arbitrary Lagrangian-Eulerian mesh motion
include("ALE/ALE.jl")
using .ALE

# Adaptive mesh refinement
include("Adaptivity/Adaptivity.jl")
using .Adaptivity

# Spatial discretization (finite volume)
include("Discretization/Discretization.jl")
using .Discretization

# Temporal integration
include("TimeIntegration/TimeIntegration.jl")
using .TimeIntegration

# Problem definition and setup
include("Problem/Problem.jl")
using .Problem

# Solution algorithms
include("Solvers/Solvers.jl")
using .Solvers

# Post-processing and output
include("PostProcess/PostProcess.jl")
using .PostProcess

# ============================================================================
# Re-exports from submodules
# ============================================================================

# Core exports
export R_GAS, STEFAN_BOLTZMANN, T_REF, P_REF
export AbstractPropertyFunction, ConstantProperty, LinearProperty
export PolynomialProperty, ArrheniusProperty, TableProperty, CallableProperty
export StateDependentProperty  # For concentration-dependent heats of reaction
export Phase, SOLID, LIQUID, GAS
export AbstractComponent, Component, SolidComponent, GaseousComponent, LiquidComponent
export AbstractReaction, Reaction, ReactionSet
export AbstractMaterial, Material, MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN
export n_components, n_reactions, n_solids, n_gases
export AbstractBoundaryCondition, AbstractThermalBC, AbstractMassBC, AbstractMeshBC

# Radiation model exports
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
export Limiter, MINMOD, VANLEER, SUPERBEE, NONE, CENTRAL
export requires_absorption_coefficients, has_nonlocal_jacobian_coupling, requires_radiation_cache
export AbstractRadiationBC, MarshakRadiationBC, ReflectiveRadiationBC, DirichletRadiationBC
export RadiationBCSet, get_radiation_bc, default_radiation_bcs
export has_absorption_coefficients

# Unified parameter evaluation
export evaluate_param, evaluate_param_T

# Boundary condition exports
export AdiabaticBC, HeatFluxBC, ConvectiveBC, RadiativeBC, DirichletThermalBC
export CombinedThermalBC, FlameBC, RadiativeFluxBC
export ImpermeableBC, MassFluxBC, ConvectiveMassBC, DiffusiveMassBC
export DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC  # Pressure BCs for Darcy
export FixedMeshBC, FreeSurfaceBC, PrescribedVelocityBC, SlidingBC
export get_radiative_flux_bc, extract_radiative_flux_from_bcs, has_radiative_flux_bc
export BoundaryConditionSet, get_thermal_bc, get_mass_bc  # Type-stable BC container

# Geometry exports
export AbstractNode, AbstractCell, AbstractFace, AbstractMesh
export Node1D, Cell1D, Face1D, CellState, Mesh1D
export create_uniform_mesh, create_stretched_mesh, update_geometry!
export n_nodes, n_cells, n_faces, n_active_cells
export cell_volume, face_area, cell_center, cell_size
export total_volume, domain_bounds, boundary_face_ids
export boundary_cell_id, face_area_value

# Properties exports
export PropertyCache, update_properties!, resize_cache!, get_cell_properties
export mixture_density, volume_fractions, solid_volume_fractions, mass_fractions, total_density, porosity
export effective_heat_capacity, effective_specific_heat, effective_conductivity
export effective_emissivity, effective_absorption, effective_gas_transfer
export harmonic_mean_weighted, harmonic_mean_product_weighted
export RadiationCache, resize_radiation_cache!, initialize_radiation_intensity!
export radiation_energy_density, net_radiation_source, total_radiation_flux

# Physics exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export conductive_flux, compute_conductive_fluxes!
export compute_radiation_source!, compute_surface_radiation_source!, compute_beer_lambert_source!
export total_absorbed_radiation, transmitted_radiation
export gas_flux, compute_gas_fluxes!, surface_mass_flux
export compute_gas_convective_source!
export compute_energy_sources!, SourceTermCache, resize_source_cache!
export volume_change_rate, compute_cell_volume_rates!
export compute_radiation_diffusion_coefficient, compute_radiation_source_terms!
export compute_radiation_fluxes!, radiation_bc_flux
export optical_thickness, optical_thickness_cell, check_p1_validity, radiation_energy_balance

# ALE exports
export AbstractMeshVelocityMethod
export LaplacianSmoothing, SpringAnalogy, LinearElasticity, WeightedALE
export compute_mesh_velocity!, move_mesh!, check_gcl
export remap_solution!

# Adaptivity exports
export AbstractErrorIndicator, AbstractQualityMetric
export GradientIndicator, CurvatureIndicator, ReactionRateIndicator
export AMRController

# TimeIntegration exports
export StateVector, state_vector_length
export create_state_vector, pack_state!, unpack_state!
export PyrolysisODEFunction, PyrolysisODEFunctionWithRadiation
export compute_timestep, compute_stable_timestep, TimestepController
export jacobian_sparsity, jacobian_sparsity_ale, matrix_colors, make_ode_function_with_jac
export optimal_fd_step, FD_RELATIVE_STEP, FD_ABSOLUTE_MIN, zero_jacobian!

# ALE TimeIntegration exports
export StateVectorALE, state_vector_length_ale, create_state_vector_ale
export pack_state_ale!, unpack_state_ale!, temperature_view_ale
export PyrolysisODEFunctionALE, PyrolysisODEFunctionALEWithRadiation
export resize_ode_function_ale!, get_thickness_from_state

# P1 Radiation TimeIntegration exports (Week 2)
export StateVectorALERadiation, state_vector_length_ale_radiation
export create_state_vector_ale_radiation
export pack_state_ale_radiation!, unpack_state_ale_radiation!
export PyrolysisODEFunctionALEWithP1Radiation
export get_radiation_intensity, set_radiation_intensity!
export get_temperature_ale_rad, set_temperature_ale_rad!
export temperature_view_ale_rad, radiation_view
export validate_state_vector_ale_radiation

# P1 Radiation Jacobian exports (Week 3)
export JacobianCacheP1Radiation, resize_jacobian_cache_p1!
export compute_p1_radiation_jacobian!
export JacobianCacheALEWithP1, resize_jacobian_cache_ale_p1!
export compute_analytical_jacobian_ale_p1!
export PyrolysisJacobianALEWithP1, resize_jacobian_ale_p1!

# Discretization P1 radiation exports
export compute_rhs_with_p1_radiation!
export update_radiation_properties!, compute_radiation_flux_divergence!
export solve_quasi_steady_radiation!, initialize_g_from_temperature!
export apply_radiation_divergence!, compute_radiation_intensity_rhs!

# Problem exports
export PyrolysisProblem, PyrolysisSolution
export apply_initial_conditions!, set_initial_conditions!
export validate_problem, is_valid_problem
export step_profile, linear_profile, gaussian_profile, tanh_profile

# Solver exports
export solve
export ILUGMRESFactorization
export get_profile, get_cell_positions, interpolate_solution

# PostProcess exports
export mass_loss_rate, total_mass, surface_temperature

end # module Pyrolysis
