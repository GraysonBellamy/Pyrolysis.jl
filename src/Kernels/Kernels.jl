"""
	Kernels

Physics kernels of the finite-volume core. Each kernel reads the derived
state, the mesh topology and metrics, the material and its boundary
conditions, and adds face fluxes and cell sources to a
[`Contributions`](@ref) buffer; the assembler turns those into row
derivatives. Face fluxes are per unit area along the face normal; the
assembler scatters `±F_f A_f` into the two cells of each face.
"""
module Kernels

using StaticArrays
using ..Materials
using ..BoundaryConditions
using ..BoundaryConditions: BoundaryState, solve_boundary_temperature!
using ..Properties: guarded_harmonic_mean, mixture_density
using ..Physics: evaluate_reactions!, species_source_terms!, reaction_heat_source, gas_flux
using ..Meshes
using ..States

include("contributions.jl")
include("reactions.jl")
include("conduction.jl")
include("gas_transport.jl")
include("radiation.jl")
include("thermal_boundary.jl")

export Kernel, Contributions, reset!, contribute!, find_kernel
export ReactionKernel, ConductionKernel, normal_conductivity
export GasTransportKernel, outward_darcy_velocity, MOBILITY_GUARD
export BeerLambertKernel, incident_intensity, sweep_intensities!, transmitted_face
export BoundaryPatch, ThermalBoundaryKernel, boundary_energy_flux

end # module Kernels
