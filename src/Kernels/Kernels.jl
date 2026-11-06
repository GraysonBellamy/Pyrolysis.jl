"""
	Kernels

Physics kernels. Each kernel reads the derived
state, the mesh topology and metrics, the material and its boundary
conditions, and adds face fluxes and cell sources to a
[`Contributions`](@ref) buffer; the assembler turns those into row
derivatives. Face fluxes are per unit area along the face normal; the
assembler scatters `±F_f A_f` into the two cells of each face.
"""
module Kernels

using StaticArrays
using LinearAlgebra: I, dot, det, tr, issymmetric, isposdef
using ForwardDiff: value
using ..Materials
using ..BoundaryConditions
import ..BoundaryConditions
using ..BoundaryConditions: BoundaryState, solve_boundary_temperature!
using ..Properties:
    guarded_harmonic_mean, FACE_MEAN_GUARD, transport_coefficients, closed_transport_channel
using ..Physics:
    interface_flux,
    lumped_interface_flux,
    pore_half_flux,
    bounded_transport_root,
    implicit_transport_root,
    primal,
    TransportFailure
using ..Physics: evaluate_reactions!, species_source_terms!, reaction_heat_source
using ..Meshes
using ..States

include("contributions.jl")
include("reactions.jl")
include("conduction.jl")
include("gas_transport.jl")
include("radiation.jl")
include("radiative_transfer.jl")
include("thermal_boundary.jl")
include("gas_boundary.jl")
include("relative_flux.jl")

export Kernel, Contributions, reset!, contribute!, find_kernel
export ReactionKernel
export ConductionKernel, normal_conductivity, cell_gradient, gradient_support
export GasTransportKernel, face_peclet
export BeerLambertKernel, incident_intensity, sweep_intensities!, transmitted_face
export RadiativeTransferKernel, MirrorBack, OpenBack, consistent_radiation!
export radiative_transfer_rows, radiation_reach, incoming_radiation
export DIFFUSIVITY, INTERFACE_FLOOR
export exposed_face
export BoundaryPatch, ThermalBoundaryKernel, boundary_energy_flux, boundary_condition
export conducted_flux, boundary_conductivity
export RelativeFluxKernel, van_albada

end # module Kernels
