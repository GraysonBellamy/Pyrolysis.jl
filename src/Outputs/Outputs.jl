"""
	Outputs

The solution of the finite-volume core: the saved times and, at each, the
fields reconstructed from the state through the state mapping and the
kernels.
"""
module Outputs

using ..Meshes
using ..States
using ..States: state_length
using ..Materials: is_gas
using ..Kernels
using ..Assembly
using ..Deactivation

include("solution.jl")

export CoreSolution, SolverStatistics, cell_positions

end # module Outputs
