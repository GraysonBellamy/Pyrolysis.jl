"""
	Problem

Problem definition for the Pyrolysis.jl package.

This module provides:
- PyrolysisProblem type
- Initial condition specification and the initial fields at the cell centres
- Problem validation
"""
module Problem

using ..Materials
using ..BoundaryConditions
using ..Materials: n_reactions
using ..BoundaryConditions: ImpermeableBC
import ..Materials: n_components  # Import to extend
import ..Meshes
using ..Properties
using ..Physics
import ..Materials: verify_stoichiometry, mass_balance_error

# Include subfiles
include("pyrolysis_problem.jl")
include("initial_conditions.jl")
include("problem_validation.jl")

# Exports
export PyrolysisProblem

# Initial conditions
export initial_fields
export step_profile, linear_profile, gaussian_profile, tanh_profile

# Validation
export validate_problem, is_valid_problem

end # module Problem
