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
using StaticArrays: SMatrix
using FunctionWrappers: FunctionWrapper
using SciMLBase: NullParameters, NoInit, FullSpecialize, derivative_discontinuity!
using SciMLBase: terminate!
using OrdinaryDiffEq: OrdinaryDiffEq, ODEFunction, ODEProblem, DiscreteCallback, CallbackSet
using OrdinaryDiffEq: ContinuousCallback, Rodas5P, FBDF, ReturnCode, set_proposed_dt!
using SciMLBase: SciMLBase, get_proposed_dt, auto_dt_reset!
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve: LinearSolve, KLUFactorization
using ADTypes: AutoFiniteDiff
using ..Materials
using ..Properties: skeletal_volume_fraction
using ..BoundaryConditions: ImpermeableBC
using ..Meshes
using ..States
using ..Kernels
using ..Closures
using ..Assembly
using ..Deactivation
using ..ColoredJacobian
using ..Outputs
using ..Problem: PyrolysisProblem, apply_initial_conditions!, validate_problem

include("erased.jl")
include("policy.jl")
include("linear_solve.jl")
include("setup.jl")
include("checks.jl")
include("solve.jl")

export ErasedRHS, ErasedJacobian
export IntegratorPolicy, default_policy, apply_observer_flag!, DT_AFTER_EVENT
export StructuredKLU, default_linear_solver
export kencarp4_policy, rodas5p_policy, fbdf_policy
export DomainCheck, ConstraintMonitor, constraint_residual, MAX_CONSTRAINT_REJECTIONS
export FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS, NewtonCap, cap_step!
export CoreProblem,
    setup, consistent_state!, inadmissible, below_absolute_zero, solve, integrate
export CORRECTION_BAND

end # module Integration
