"""
	Properties

Material property evaluation for the Pyrolysis.jl package.

This module provides:
- Mixture property computations
- Effective property calculations
- Property caching for performance
"""
module Properties

using ..Core
using ..Geometry

# Include subfiles
include("mixture_rules.jl")
include("effective_properties.jl")
include("property_cache.jl")
include("radiation_cache.jl")

# Exports
export mixture_density, volume_fractions, solid_volume_fractions, mass_fractions, total_density, porosity
export porosity_derivatives, porosity_derivatives!
export effective_heat_capacity, effective_specific_heat, effective_conductivity
export effective_emissivity, effective_absorption, effective_gas_transfer
export effective_permeability, permeability_derivatives, permeability_derivatives!  # Component-based permeability
export harmonic_mean_weighted, harmonic_mean_product_weighted
export effective_gas_viscosity, mixture_gas_viscosity
export species_transport_scaling, M_REF_TRANSPORT
export collision_integral_diffusion, chapman_enskog_diffusivity, chapman_enskog_λ
export CHAPMAN_ENSKOG_PREFACTOR, SIGMA_AIR, EPSILON_K_AIR, M_AIR
export bruggeman_effective_property
export PropertyCache, update_properties!, resize_cache!, get_cell_properties
export harmonic_mean_permeability, harmonic_mean_permeability_weighted

# Radiation cache (P1 approximation)
export RadiationCache, resize_radiation_cache!
export initialize_radiation_intensity!
export radiation_energy_density, net_radiation_source, total_radiation_flux

end # module Properties
