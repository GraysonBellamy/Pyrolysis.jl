"""
	Physics

Physical models for pyrolysis simulation.

This module provides:
- Reaction kinetics
- Heat conduction
- Radiation transport
- Gas transport
- Source terms
"""
module Physics

using ..Core
using ..Geometry
using ..Properties
using StaticArrays
using ForwardDiff  # For _value() helper in radiation_p1.jl

# Include subfiles
include("kinetics.jl")
include("heat_conduction.jl")
include("radiation.jl")
include("radiation_p1.jl")
include("gas_transport.jl")
include("convection.jl")
include("source_terms.jl")
include("volume_change.jl")
include("pressure.jl")

# Exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export conductive_flux, compute_conductive_fluxes!
export compute_radiation_source!, compute_surface_radiation_source!, compute_beer_lambert_source!
export total_absorbed_radiation, transmitted_radiation
export gas_flux, compute_gas_fluxes!, surface_mass_flux
export gas_flux_darcy_fick, gas_flux_advective, compute_gas_fluxes_darcy_fick!
export compute_gas_convective_source!
export compute_energy_sources!, SourceTermCache, resize_source_cache!
export volume_change_rate, compute_cell_volume_rates!
export ideal_gas_pressure, compute_pressure!, compute_pressure_gradient!, pressure_derivatives, pressure_derivatives!

# P1 radiation transport exports
export compute_radiation_diffusion_coefficient
export compute_radiation_source_terms!, update_radiation_properties!
export compute_radiation_fluxes!, radiation_bc_flux
export optical_thickness, optical_thickness_cell, check_p1_validity
export radiation_energy_balance

end # module Physics
