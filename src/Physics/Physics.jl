"""
	Physics

Physical models for pyrolysis simulation.

This module provides:
- Reaction kinetics
- The diffusive gas flux at a face
- Volume change and gas pressure
- The phase change of a liquid and its vapour
"""
module Physics

using StaticArrays, ForwardDiff, LinearAlgebra
using ..Materials
using ..Properties

# Include subfiles
include("kinetics.jl")
include("constitutive_transport.jl")
include("interface_transport.jl")
include("volume_change.jl")
include("pressure.jl")
include("phase_change.jl")

# Exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export DepletionLimiter, one_sided_factor
export molar_wall_flux, pore_half_flux, interface_flux, lumped_interface_flux
export volume_change_rate
export pore_state, ideal_gas_pressure
export background_gas_concentration
export moisture_content, equilibrium_vapour_pressure, phase_change_rate, phase_change_heat

end # module Physics
