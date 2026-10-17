"""
	ColoredJacobian

The Jacobian of the finite-volume core: a sparsity pattern built from the
mesh connectivity, column coloring, and colored forward-mode automatic
differentiation of the residual, with dense forward-mode differentiation as
the reference that the pattern must contain and the values must match.
"""
module ColoredJacobian

using SparseArrays
using ForwardDiff
using ADTypes: AutoForwardDiff, AutoSparse, KnownJacobianSparsityDetector
using DifferentiationInterface: DifferentiationInterface, Constant, prepare_jacobian
using SparseMatrixColorings: GreedyColoringAlgorithm, column_groups, sparsity_pattern
using ..Meshes
using ..States
using ..Assembly

include("pattern.jl")
include("oracle.jl")
include("colored.jl")

export jacobian_pattern
export dense_jacobian, missing_entries, compare_entries
export ResidualFunction, ColoredJacobianCache, jac_prototype, jacobian!, same_structure

end # module ColoredJacobian
