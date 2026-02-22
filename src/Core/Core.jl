"""
	Core

Fundamental types, constants, and abstractions for the Pyrolysis.jl package.

This module provides:
- Physical constants (R_GAS, STEFAN_BOLTZMANN, etc.)
- Temperature-dependent property functions
- Component and material type definitions
- Reaction type definitions
- Boundary condition abstractions
"""
module Core

using StaticArrays

# ============================================================================
# Abstract Type Hierarchy (must be defined BEFORE includes)
# ============================================================================

# Property functions
abstract type AbstractPropertyFunction end
abstract type AbstractThermalProperty <: AbstractPropertyFunction end
abstract type AbstractTransportProperty <: AbstractPropertyFunction end

# Material system
abstract type AbstractComponent end
abstract type AbstractReaction end
abstract type AbstractMaterial end

# Geometry (defined here, implemented in Geometry module)
abstract type AbstractNode end
abstract type AbstractCell end
abstract type AbstractFace end
abstract type AbstractMesh{DIM} end

# Boundary conditions
abstract type AbstractBoundaryCondition end
abstract type AbstractThermalBC <: AbstractBoundaryCondition end
abstract type AbstractMassBC <: AbstractBoundaryCondition end
abstract type AbstractMeshBC <: AbstractBoundaryCondition end
abstract type AbstractRadiationBC <: AbstractBoundaryCondition end

# ALE (defined here, implemented in ALE module)
abstract type AbstractMeshVelocityMethod end
abstract type AbstractRemappingMethod end
abstract type AbstractInterfaceTracker end

# Adaptivity (defined here, implemented in Adaptivity module)
abstract type AbstractErrorIndicator end
abstract type AbstractRefinementCriterion end
abstract type AbstractQualityMetric end

# Problem and solution
abstract type AbstractPyrolysisProblem end
abstract type AbstractPyrolysisSolution end

# ============================================================================
# Include all subfiles (AFTER abstract types are defined)
# ============================================================================
include("constants.jl")
include("property_functions.jl")
include("param_evaluation.jl")  # Unified parameter evaluation functions
include("components.jl")
include("reactions.jl")
include("materials.jl")
include("boundary_conditions.jl")
include("radiation_types.jl")
include("radiation_boundary.jl")

# ============================================================================
# Exports
# ============================================================================

# Constants
export R_GAS, STEFAN_BOLTZMANN, T_REF, P_REF

# Abstract types
export AbstractPropertyFunction, AbstractThermalProperty, AbstractTransportProperty
export AbstractComponent, AbstractReaction, AbstractMaterial
export AbstractNode, AbstractCell, AbstractFace, AbstractMesh
export AbstractBoundaryCondition, AbstractThermalBC, AbstractMassBC, AbstractMeshBC
export AbstractMeshVelocityMethod, AbstractRemappingMethod, AbstractInterfaceTracker
export AbstractErrorIndicator, AbstractRefinementCriterion, AbstractQualityMetric
export AbstractPyrolysisProblem, AbstractPyrolysisSolution

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty
export ArrheniusProperty, TableProperty, CallableProperty
export StateDependentProperty  # For state-dependent heats of reaction
export derivative, derivative_ξ  # Temperature and concentration derivatives of property functions

# Components
export Phase, SOLID, LIQUID, GAS
export Component, SolidComponent, GaseousComponent, LiquidComponent
export get_permeability, has_permeability  # Permeability accessors for components

# Reactions
export Reaction, ReactionSet

# Materials
export Material, MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN
export n_components, n_reactions, n_solids, n_gases
export is_gas_phase, is_solid_phase, get_swelling_factor, get_component_density
export foreach_component, foreach_gas_component, gas_component_indices
export compute_gas_fluxes_unrolled!, compute_gas_fluxes_darcy_fick_unrolled!
export compute_gas_fluxes_species_scaled_unrolled!, compute_gas_fluxes_darcy_fick_species_scaled_unrolled!
export apply_species_residual_unrolled!

# Boundary conditions
export AdiabaticBC, HeatFluxBC, ConvectiveBC, RadiativeBC, DirichletThermalBC
export CombinedThermalBC, FlameBC, RadiativeFluxBC
export ImpermeableBC, MassFluxBC, ConvectiveMassBC, DiffusiveMassBC, CombinedMassBC
export DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC  # Pressure BCs for Darcy flow
export FixedMeshBC, FreeSurfaceBC, PrescribedVelocityBC, SlidingBC
export get_radiative_flux_bc, extract_radiative_flux_from_bcs, has_radiative_flux_bc

# Type-stable boundary condition container
export BoundaryConditionSet, get_thermal_bc, get_mass_bc

# Gas transport (component-based permeability)
export has_darcy_flow

# Radiation types and boundary conditions
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
export requires_absorption_coefficients, has_nonlocal_jacobian_coupling, requires_radiation_cache
export AbstractRadiationBC, MarshakRadiationBC, ReflectiveRadiationBC, DirichletRadiationBC
export RadiationBCSet, get_radiation_bc
export default_radiation_bcs
export has_absorption_coefficients

# Unified parameter evaluation (replaces _evaluate_radiation_param*)
export evaluate_param, evaluate_param_T

end # module Core
