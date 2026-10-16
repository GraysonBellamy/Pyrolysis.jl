"""
	Physics

Physical models for pyrolysis simulation.

This module provides:
- Reaction kinetics
- Heat conduction
- Radiation transport
- Gas transport
- Gas convective energy source
- Volume change and gas pressure
"""
module Physics

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties

# Include subfiles
include("kinetics.jl")
include("heat_conduction.jl")
include("radiation.jl")
include("gas_transport.jl")
include("convection.jl")
include("volume_change.jl")
include("pressure.jl")

# Exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export DepletionLimiter
export conductive_flux, compute_conductive_fluxes!
export compute_radiation_source!, compute_beer_lambert_source!
export gas_flux, compute_gas_fluxes!
export gas_flux_darcy_fick, compute_gas_fluxes_darcy_fick!
export compute_gas_convective_source!
export volume_change_rate
export ideal_gas_pressure
export background_gas_concentration

# P1 radiation transport (compute_radiation_diffusion_coefficient,
# compute_radiation_source_terms!, update_radiation_properties!,
# compute_radiation_fluxes!, radiation_bc_flux, optical_thickness,
# optical_thickness_cell, check_p1_validity, radiation_energy_balance)
# lives in Pyrolysis.Experimental.

end # module Physics
