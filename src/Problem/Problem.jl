"""
	Problem

Problem definition and setup for the Pyrolysis.jl package.

This module provides:
- PyrolysisProblem type
- Initial condition specification
- Problem validation
"""
module Problem

using ..Materials
using ..BoundaryConditions
using ..Materials: n_reactions
using ..BoundaryConditions: ImpermeableBC
import ..Materials: n_components  # Import to extend
using ..Geometry
# Qualified: Geometry exports names of its own mesh that Meshes shares.
import ..Meshes
using ..Properties
using ..Physics
using ..Discretization
using ForwardDiff
import ..Materials: verify_stoichiometry, mass_balance_error

# Include subfiles
include("pyrolysis_problem.jl")
include("initial_conditions.jl")
include("problem_validation.jl")

# Unified pure-residual + Workspace surface.
include("traits.jl")
include("state_layout.jl")
include("state_cache.jl")
include("problem_def.jl")
include("workspace.jl")
include("adapter.jl")

# Exports
# Types
export PyrolysisProblem, PyrolysisSolution

# Problem accessors
export get_mesh, get_material, tspan

# Solution accessors
export get_times, get_temperatures, get_concentrations, get_gas_generation_rate
export n_timesteps, succeeded

# Initial conditions
export apply_initial_conditions!, initial_fields
export step_profile, linear_profile, gaussian_profile, tanh_profile

# Validation
export validate_problem, is_valid_problem

# Pure-residual + Workspace surface.
#
# Trait subtype names (`Eulerian`, `ALE`, `Fickian`, ...) are deliberately
# **not** re-exported to keep the top-level surface small and avoid leaking
# short generic names (`ALE`, `Fickian`) into `using Pyrolysis`. Reach for
# trait structs as `Pyrolysis.Problem.ALE()`, `Pyrolysis.Problem.BeerLambert()`,
# etc., or import them explicitly: `using Pyrolysis.Problem: ALE, BeerLambert`.
export SimulationMode, MeshMode, RadiationMode, TransportMode, ChiMode
export mesh_mode, radiation_mode, transport_mode, chi_mode
export StateLayout,
    state_length, temperature_range, concentration_range, z_range, chi_range, block_view
export pack_state!, update_mesh_from_state!, create_state_vector
export StateCache,
    populate!, mirror_eulerian_geometry!, update_cross_section_area!, update_state_geometry!
export get_temperature,
    get_concentration,
    get_node_position,
    get_node_positions,
    get_chi,
    get_thickness_from_state,
    compute_volumes_from_state
export ProblemDef, get_option
export Workspace, update_radiation!
export initialize_state
export to_problem_def, build_workspace, build_ad_workspace

end # module Problem
