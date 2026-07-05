"""
	Properties

Material property evaluation for the Pyrolysis.jl package.

This module provides:
- Mixture property computations
- Effective property calculations
- Property caching for performance
"""
module Properties

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ForwardDiff  # permeability_derivatives: forward-mode over effective_permeability

# Include subfiles
include("mixture_rules.jl")
include("effective_properties.jl")
include("property_cache.jl")
include("radiation_cache.jl")

# Exports
export mixture_density,
    volume_fractions, solid_volume_fractions, mass_fractions, total_density, porosity
export effective_heat_capacity, total_heat_capacity, effective_conductivity
export effective_emissivity, effective_absorption, effective_gas_transfer
export effective_permeability, permeability_derivatives  # Component-based permeability
export harmonic_mean_weighted
export effective_gas_viscosity, mixture_gas_viscosity
export collision_integral_diffusion, chapman_enskog_diffusivity
export CHAPMAN_ENSKOG_PREFACTOR, SIGMA_AIR, EPSILON_K_AIR, M_AIR
export bruggeman_effective_property
export PropertyCache, update_properties!, resize_cache!
export harmonic_mean_permeability

# Radiation cache (P1 approximation)
export RadiationCache, resize_radiation_cache!
export initialize_radiation_intensity!
export radiation_energy_density, net_radiation_source, total_radiation_flux

end # module Properties
