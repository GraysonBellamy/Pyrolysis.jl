"""
	Materials

Foundational types, constants, and abstractions for the Pyrolysis.jl package.

Holds the foundational layer — the abstract
type hierarchy, physical constants, temperature-dependent property functions,
the component/material/reaction type system, and the radiation-model enum
(`RadiationModel`, kept here for trait dispatch). The concrete boundary-condition
types were split out into the `BoundaryConditions` module, which depends on this
one for the abstract BC supertypes and property functions.

This module provides:
- Physical constants (R_GAS, STEFAN_BOLTZMANN, etc.)
- Temperature-dependent property functions
- Component and material type definitions
- Reaction type definitions
- The abstract boundary-condition / radiation-BC supertypes (concretes live in
  `BoundaryConditions`)
- The `RadiationModel` enum + predicates (trait dispatch)
"""
module Materials

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
abstract type AbstractRadiationBC <: AbstractBoundaryCondition end

# ALE (defined here, implemented in ALE module)
abstract type AbstractRemappingMethod end

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
include("components.jl")
include("reactions.jl")
include("materials.jl")
include("radiation_types.jl")

# ============================================================================
# Exports
# ============================================================================

# Constants
export R_GAS, STEFAN_BOLTZMANN, T_REF, P_REF, RHO_CP_FLOOR

# Abstract types
export AbstractPropertyFunction, AbstractThermalProperty, AbstractTransportProperty
export AbstractComponent, AbstractReaction, AbstractMaterial
export AbstractNode, AbstractCell, AbstractFace, AbstractMesh
export AbstractBoundaryCondition, AbstractThermalBC, AbstractMassBC
export AbstractRadiationBC
export AbstractRemappingMethod
export AbstractErrorIndicator, AbstractRefinementCriterion, AbstractQualityMetric
export AbstractPyrolysisProblem, AbstractPyrolysisSolution

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty
export ArrheniusProperty, TableProperty, CallableProperty
export StateDependentProperty  # For state-dependent heats of reaction
export derivative, derivative_ξ  # Temperature and concentration derivatives of property functions
export evaluate_param, evaluate_param_T  # Unified parameter evaluation (defined in property_functions.jl)

# Components
export Phase, SOLID, LIQUID, GAS
export SolidComponent, GaseousComponent, LiquidComponent
export phase, molar_mass, is_gas, is_solid, is_liquid  # Type-dispatch accessors
export delta_enthalpy                                    # Sensible-enthalpy difference (midpoint rule)
export get_permeability, has_permeability  # Permeability accessors for components
export get_density, get_intrinsic_density  # Bulk vs. skeletal density accessors

# Reactions
export Reaction, ReactionSet
export DepletionLimiter

# Materials
export Material, MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY
export MixingSpec, MaterialMixing
export n_components, n_reactions, n_solids, n_gases
export is_gas_phase, is_solid_phase, get_swelling_factor, get_component_density
export foreach_component,
    foreach_gas_component, gas_component_indices, dry_pyrolysis_gas_indices
export compute_gas_fluxes_unrolled!, compute_gas_fluxes_darcy_fick_unrolled!
export set_boundary_advective_flux_unrolled!
export apply_species_residual_unrolled!

# Gas transport (component-based permeability)
export has_darcy_flow

# Radiation model enum + predicates (trait dispatch; kept in the foundation)
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
export requires_absorption_coefficients,
    has_nonlocal_jacobian_coupling, requires_radiation_cache
export has_absorption_coefficients

end # module Materials
