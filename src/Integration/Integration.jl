"""
	Integration

Time integration of the finite-volume core: the adapter from a
`PyrolysisProblem`, the mass-matrix `ODEFunction` with the colored Jacobian
behind problem-independent function types, the integrator together with its
restart and per-step-flag policies, the admissibility check, and `solve`.
"""
module Integration

using LinearAlgebra
using SparseArrays
using FunctionWrappers: FunctionWrapper
using SciMLBase: NullParameters, NoInit, FullSpecialize, derivative_discontinuity!
using OrdinaryDiffEq: OrdinaryDiffEq, ODEFunction, ODEProblem, DiscreteCallback, CallbackSet
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve: KLUFactorization
using ADTypes: AutoFiniteDiff
using ..Materials
using ..BoundaryConditions: ImpermeableBC
using ..Meshes
using ..States
using ..Kernels
using ..Assembly
using ..ColoredJacobian
using ..Outputs
using ..Problem: PyrolysisProblem, apply_initial_conditions!, validate_problem

include("erased.jl")
include("policy.jl")
include("setup.jl")
include("solve.jl")

export ErasedRHS, ErasedJacobian
export IntegratorPolicy, default_policy, apply_observer_flag!
export CoreProblem, setup, inadmissible, solve

end # module Integration
