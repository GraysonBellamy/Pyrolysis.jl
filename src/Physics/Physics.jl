"""
	Physics

Physical models for pyrolysis simulation.

This module provides:
- Reaction kinetics
- Heat conduction
- Radiation transport
- Gas transport
- Source terms
- ALE mesh motion / GCL / conservative remapping
"""
module Physics

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties
using StaticArrays
using LinearAlgebra
using SparseArrays

# Include subfiles
include("kinetics.jl")
include("heat_conduction.jl")
include("radiation.jl")
include("gas_transport.jl")
include("convection.jl")
include("source_terms.jl")
include("volume_change.jl")
include("pressure.jl")
# ALE kernels (mesh motion with GCL compliance, Geometric Conservation Law
# verification, conservative remapping).
include("mesh_motion.jl")
include("gcl.jl")
include("remapping.jl")

# Exports
export reaction_rate, evaluate_reactions!, species_source_terms!, reaction_heat_source
export DepletionLimiter
export conductive_flux, compute_conductive_fluxes!
export compute_radiation_source!, compute_beer_lambert_source!
export total_absorbed_radiation, transmitted_radiation
export gas_flux, compute_gas_fluxes!, surface_mass_flux
export gas_flux_darcy_fick, gas_flux_advective, compute_gas_fluxes_darcy_fick!
export compute_gas_convective_source!
export compute_energy_sources!, SourceTermCache, resize_source_cache!
export volume_change_rate, compute_cell_volume_rates!
export ideal_gas_pressure, compute_pressure!, compute_pressure_gradient!
export background_gas_concentration

# P1 radiation transport (compute_radiation_diffusion_coefficient,
# compute_radiation_source_terms!, update_radiation_properties!,
# compute_radiation_fluxes!, radiation_bc_flux, optical_thickness,
# optical_thickness_cell, check_p1_validity, radiation_energy_balance)
# lives in Pyrolysis.Experimental.

# ALE mesh motion
export move_mesh!, move_mesh_safe!
export check_mesh_quality, mesh_is_valid
export mesh_displacement, max_mesh_displacement, reset_mesh_to_reference!
export update_face_positions!, update_cell_centers!

# ALE — Geometric Conservation Law
export check_gcl, gcl_error
export gcl_volume_update, update_volumes_from_geometry!, gcl_residual

# ALE — conservative remapping
export FirstOrderRemap, SecondOrderRemap
export remap_solution!
export check_remap_conservation

end # module Physics
