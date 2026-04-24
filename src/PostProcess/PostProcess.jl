"""
	PostProcess

Post-processing and output for the Pyrolysis.jl package.

This module provides:
- Standard output computations
- Conservation checks
- Visualization recipes
- Data export
"""
module PostProcess

using ..Core
using ..Geometry
using ..Properties

# Include subfiles
include("conservation_checks.jl")

# Exports
export mass_loss_rate, total_mass, surface_temperature, back_face_temperature
export heat_release_rate
export MassConservationTracker, EnergyConservationTracker
export compute_domain_mass, compute_sensible_energy
export update_mass_tracker!, update_energy_tracker!
export check_mass_conservation, check_energy_conservation
export mass_conservation_report, energy_conservation_report

end # module PostProcess
