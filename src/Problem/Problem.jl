"""
	Problem

Problem definition and setup for the Pyrolysis.jl package.

This module provides:
- PyrolysisProblem type
- Initial condition specification
- Problem validation
"""
module Problem

using ..Core
using ..Core: n_reactions, ImpermeableBC, SOLID
import ..Core: n_components  # Import to extend
using ..Geometry
using ..Properties
using ..Physics
using ..ALE
using ..Adaptivity
using ..Discretization
using ..TimeIntegration

# Include subfiles
include("pyrolysis_problem.jl")
include("initial_conditions.jl")
include("problem_validation.jl")

# Exports
# Types
export PyrolysisProblem, PyrolysisSolution

# Problem accessors
export get_mesh, get_material, tspan

# Solution accessors
export get_times, get_temperatures, get_concentrations, get_gas_generation_rate
export n_timesteps, succeeded

# Initial conditions
export apply_initial_conditions!, set_initial_conditions!
export step_profile, linear_profile, gaussian_profile, tanh_profile
export verify_initial_conditions

# Validation
export validate_problem, is_valid_problem, check_consistency

end # module Problem
