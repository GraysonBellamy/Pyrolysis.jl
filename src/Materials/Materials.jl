"""
	Materials

Foundational types, constants, and abstractions for the Pyrolysis.jl package.

Holds the foundational layer — the abstract
type hierarchy, physical constants, temperature-dependent property functions,
the component/material/reaction type system, and the radiation-model enum
(`RadiationModel`). The concrete boundary-condition
types live in the `BoundaryConditions` module, which depends on this one for the
abstract BC supertypes and property functions.

This module provides:
- Physical constants (R_GAS, STEFAN_BOLTZMANN, etc.)
- Temperature-dependent property functions
- Component and material type definitions
- Reaction type definitions
- The abstract boundary-condition supertypes (concretes live in
  `BoundaryConditions`)
- The `RadiationModel` enum
"""
module Materials

# ============================================================================
# Abstract Type Hierarchy (must be defined BEFORE includes)
# ============================================================================

# Property functions
abstract type AbstractPropertyFunction end

# Material system
abstract type AbstractComponent end

# Boundary conditions
abstract type AbstractThermalBC end
abstract type AbstractMassBC end

# ============================================================================
# Include all subfiles (AFTER abstract types are defined)
# ============================================================================
include("constants.jl")
include("guards.jl")
include("property_functions.jl")
include("components.jl")
include("reactions.jl")
include("gas_transport.jl")
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
export AbstractComponent
export AbstractThermalBC, AbstractMassBC

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty
export ArrheniusProperty, TableProperty, CallableProperty
export SutherlandProperty  # Sutherland gas-viscosity correlation (Material.gas_viscosity default)
export StateDependentProperty  # For state-dependent heats of reaction
export evaluate_param, evaluate_param_T  # Unified parameter evaluation (defined in property_functions.jl)
export property_temperature  # Argument clamp shared by every temperature-dependent law
export identically_zero  # Structural zero test of a property law

# Components
export Phase, SOLID, LIQUID, GAS
export SolidComponent, GaseousComponent, LiquidComponent
export phase, molar_mass, is_gas, is_solid, is_liquid  # Type-dispatch accessors
export get_permeability  # Permeability accessor for components
export participates_in_transport  # Gas-transport matrix membership
export get_intrinsic_density  # Skeletal density accessor
export get_emissivity, get_absorptivity, surface_absorptivity  # Radiative surface accessors

# Blowing models (evaluated in BoundaryConditions/blowing.jl)
export AbstractBlowingModel,
    NoBlowing, CouetteFilm, PropertyRatioFilm, PoolingFilm, FunctionBlowing
export GasSpecies, NITROGEN, nasa7_heat_capacity, eucken_viscosity, fuller_diffusivity
export gas_property_temperature

# Reactions
export Reaction, DepletionLimiter

# Materials
export Material, MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY
export MixingSpec, MaterialMixing
export n_components, n_reactions
export dry_gas_flags
export needs_lineage, lineage_graph, LineageGraph

# Gas transport (component-based permeability)
export gas_component_indices
export LumpedTransfer, DustyGas, LumpedEndpoint, DustyEndpoint, has_gas
export smooth_clamp, trial_temperature, positive_part, smooth_floor

# Radiation model enum
export RadiationModel, SURFACE_ABSORPTION, BEER_LAMBERT

end # module Materials
