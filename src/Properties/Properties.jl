"""
	Properties

Material property evaluation for the Pyrolysis.jl package: the mixture rules
(volume fractions, density, the skeletal volume) and the effective properties
(conductivity, surface emissivity and absorptivity, absorption, gas transfer,
permeability, the gas diffusivity) a cell's composition and temperature give.
"""
module Properties

using ..Materials

# Include subfiles
include("mixture_rules.jl")
include("effective_properties.jl")

# Exports
export mixture_density, volume_fractions, solid_volume_fractions
export skeletal_volume_fraction
export effective_conductivity
export effective_absorption
export mixture_emissivity, mixture_absorptivity
export effective_gas_transfer
export effective_permeability  # Component-based permeability
export guarded_harmonic_mean
export collision_integral_diffusion, chapman_enskog_diffusivity
export CHAPMAN_ENSKOG_PREFACTOR, SIGMA_AIR, EPSILON_K_AIR, M_AIR
export bruggeman_effective_property

end # module Properties
