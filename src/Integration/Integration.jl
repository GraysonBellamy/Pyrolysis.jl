"""
	Integration

Time integration: the adapter from a
`PyrolysisProblem`, the mass-matrix `ODEFunction` with the colored Jacobian
behind problem-independent function types, the integrator together with its
restart and per-step-flag policies, the admissibility check, `solve`, its
read-only observers (the progress bar and the conservation recorder), and
parameter sweeps.
"""
module Integration
using ..Materials: gas_component_indices

using ..Properties: transport_domain, closed_transport_channel

using LinearAlgebra
using SparseArrays
using StaticArrays: SMatrix
using FunctionWrappers: FunctionWrapper
using SciMLBase: NullParameters, NoInit, FullSpecialize, derivative_discontinuity!
using SciMLBase: terminate!
using OrdinaryDiffEq: OrdinaryDiffEq, ODEFunction, ODEProblem, DiscreteCallback, CallbackSet
using OrdinaryDiffEq: ContinuousCallback, Rodas5P, FBDF, ReturnCode, set_proposed_dt!
using OrdinaryDiffEq: OrdinaryDiffEqCore
using SciMLBase: SciMLBase, get_proposed_dt, auto_dt_reset!
# `solve` is a method of the SciML ecosystem's common `solve`, so `using
# Pyrolysis` beside `using OrdinaryDiffEq` names one function.
import SciMLBase: solve
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve: LinearSolve, KLUFactorization
using ADTypes: AutoFiniteDiff
using ProgressMeter: ProgressMeter, Progress
using Distributed: pmap, nworkers
using ..Materials
import ..Physics
using ForwardDiff
using ..Properties: skeletal_volume_fraction
using ..BoundaryConditions: ImpermeableBC, MassFluxBC, validate_gas_boundary
using ..Meshes
using ..States
using ..Kernels
using ..Closures
using ..Assembly
using ..Deactivation
using ..ColoredJacobian
using ..Outputs
using ..Outputs: _extension
import ..Problem
using ..Problem: PyrolysisProblem, initial_fields, validate_problem

include("erased.jl")
include("policy.jl")
include("linear_solve.jl")
include("transport_support.jl")
include("setup.jl")
include("checks.jl")
include("pressure_error.jl")
include("observers.jl")
include("solve.jl")
include("parameter_sweep.jl")

export ErasedRHS, ErasedJacobian
export IntegratorPolicy, default_policy, apply_observer_flag!, DT_AFTER_EVENT
export StructuredKLU, default_linear_solver
export kencarp4_policy, rodas5p_policy, fbdf_policy
export DomainCheck, ConstraintMonitor, constraint_residual, MAX_CONSTRAINT_REJECTIONS
export FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS, NewtonCap, cap_step!
export CoreProblem,
    setup, consistent_state!, inadmissible, below_absolute_zero, solve, integrate
export parameter_sweep
export CORRECTION_BAND

end # module Integration
