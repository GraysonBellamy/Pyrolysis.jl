"""
	Physics

Physical models for pyrolysis simulation.

This module provides:
- Reaction kinetics
- The diffusive gas flux at a face
- Volume change and gas pressure
"""
module Physics

using ..Materials
using ..BoundaryConditions
using ..Properties

# Include subfiles
include("kinetics.jl")
include("gas_transport.jl")
include("volume_change.jl")
include("pressure.jl")

# Exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export DepletionLimiter, one_sided_factor
export gas_flux
export volume_change_rate
export ideal_gas_pressure
export background_gas_concentration

end # module Physics
