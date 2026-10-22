"""
	Outputs

The solution of the finite-volume core: the saved times and, at each, the
fields reconstructed from the state through the state mapping and the
kernels; and its accessors.
"""
module Outputs

using ..Meshes
using ..States
using ..States: state_length
using ..Materials: is_gas, T_REF
import ..Problem: get_times, get_temperatures, get_concentrations, get_gas_generation_rate
import ..Problem: n_timesteps, succeeded
using ..BoundaryConditions: compute_boundary_heat_fluxes, BlowingFace, DirichletThermalBC
using ..Kernels
using ..Assembly
using ..Deactivation

include("solution.jl")
include("accessors.jl")

export CoreSolution, CoreExtras, SolverStatistics
export get_times, get_temperatures, get_concentrations, get_gas_generation_rate
export n_timesteps, succeeded, completed
export get_profile, get_cell_positions, interpolate_solution

end # module Outputs
