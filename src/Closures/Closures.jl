"""
	Closures

How the node velocities of a moving mesh are determined. A closure turns the
cells' strain rates into algebraic rows on the node velocities, so every row
of the system stays local: in 1D the discrete volume balance of each cell,
swept from the pinned node, or under verification a prescribed mesh
velocity. A variable cross section adds its two algebraic rows, the area and
its lateral strain rate, which close the column's volume balance laterally.
"""
module Closures

using StaticArrays
using ForwardDiff
using ..Materials
using ..Physics: volume_change_rate
using ..Meshes
using ..States
using ..Kernels: Contributions

include("balance.jl")
include("lateral.jl")

export Closure, BalanceRowClosure, closure_rows!, consistent_velocities!, rank_check
export LateralShrinkage, LATERAL_CLAMP_WIDTH, chi_bar, lateral_rows, consistent_lateral

end # module Closures
