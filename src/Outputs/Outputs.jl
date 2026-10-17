"""
	Outputs

The solution of the finite-volume core: the saved times and, at each, the
fields reconstructed from the state through the state mapping.
"""
module Outputs

using ..Meshes
using ..States

include("solution.jl")

export CoreSolution, SolverStatistics

end # module Outputs
