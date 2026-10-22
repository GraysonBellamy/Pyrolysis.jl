"""
	Closures

How the node velocities of a moving mesh are determined. A closure turns the
cells' strain rates into algebraic rows, so every row of the system stays
local: on a line or a single column the discrete volume balance of each
cell on the node velocities, swept from the pinned nodes; on any mesh the
balances on one multiplier per cell, the node velocities the multipliers'
local function of least lumped-volume norm; or under verification a
prescribed mesh velocity. A correction moves the mesh relative to the solid,
the node velocities then being the solid's plus the correction's. A variable
cross section adds its two algebraic rows, the area and its lateral strain
rate, which close the column's volume balance laterally.
"""
module Closures

using StaticArrays
using SparseArrays: sparse
using LinearAlgebra: dot, norm, rank, cholesky, issuccess, Symmetric, Diagonal
using ForwardDiff
using ..Materials
using ..Physics: volume_change_rate
using ..Meshes
using ..States
using ..Kernels: Contributions, DISTANCE_FLOOR_FRACTION

include("correction.jl")
include("balance.jl")
include("potential.jl")
include("lateral.jl")

export Closure, BalanceRowClosure, closure_rows!, consistent_velocities!, rank_check
export PotentialClosure, consistent_multipliers, multiplier_rows!
export Correction, NoCorrection, PrescribedCorrection, EquidistributionCorrection
export LayerCorrection, EQUIDISTRIBUTION_TIME, bind_correction, activation, correction!
export vanishing, vanishes
export LateralShrinkage, LATERAL_CLAMP_WIDTH, chi_bar, lateral_rows, consistent_lateral

end # module Closures
