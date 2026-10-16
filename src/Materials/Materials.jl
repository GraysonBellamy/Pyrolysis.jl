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

# ============================================================================
# Abstract Type Hierarchy (must be defined BEFORE includes)
# ============================================================================

# Property functions
abstract type AbstractPropertyFunction end

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

# Adaptivity (defined here, implemented in Adaptivity module)
abstract type AbstractErrorIndicator end
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
include("blowing_models.jl")

# ============================================================================
# Exports
# ============================================================================

# Constants
export R_GAS, STEFAN_BOLTZMANN, T_REF, P_REF, RHO_CP_FLOOR, T_PROPERTY_MIN, T_PROPERTY_MAX

# Abstract types
export AbstractPropertyFunction
export AbstractComponent, AbstractReaction, AbstractMaterial
export AbstractNode, AbstractCell, AbstractFace, AbstractMesh
export AbstractBoundaryCondition, AbstractThermalBC, AbstractMassBC
export AbstractRadiationBC
export AbstractErrorIndicator, AbstractQualityMetric
export AbstractPyrolysisProblem, AbstractPyrolysisSolution

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty
export ArrheniusProperty, TableProperty, CallableProperty
export SutherlandProperty  # Sutherland gas-viscosity correlation (Material.gas_viscosity default)
export StateDependentProperty  # For state-dependent heats of reaction
export evaluate_param, evaluate_param_T  # Unified parameter evaluation (defined in property_functions.jl)
export property_temperature  # Argument clamp shared by every temperature-dependent law

# Components
export Phase, SOLID, LIQUID, GAS
export SolidComponent, GaseousComponent, LiquidComponent
export phase, molar_mass, is_gas, is_solid, is_liquid  # Type-dispatch accessors
export get_permeability, has_permeability  # Permeability accessors for components
export get_intrinsic_density  # Skeletal density accessor
export get_emissivity, get_absorptivity, surface_absorptivity  # Radiative surface accessors

# Blowing models (evaluated in Discretization/blowing.jl)
export AbstractBlowingModel,
    NoBlowing, CouetteFilm, PropertyRatioFilm, PoolingFilm, FunctionBlowing
export GasSpecies, NITROGEN, nasa7_heat_capacity, eucken_viscosity, fuller_diffusivity
export gas_property_temperature

# Reactions
export Reaction
export DepletionLimiter

# Materials
export Material, MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY
export MixingSpec, MaterialMixing
export n_components, n_reactions
export dry_pyrolysis_gas_indices
export compute_gas_fluxes_unrolled!, compute_gas_fluxes_darcy_fick_unrolled!
export set_boundary_advective_flux_unrolled!
export apply_species_residual_unrolled!

# Gas transport (component-based permeability)
export has_darcy_flow

# Radiation model enum + predicates (trait dispatch; kept in the foundation)
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
export requires_radiation_cache
export has_absorption_coefficients

end # module Materials
