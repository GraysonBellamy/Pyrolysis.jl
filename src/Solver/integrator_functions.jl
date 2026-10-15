# ============================================================================
# Functions the integrator calls, behind problem-independent types
# ============================================================================
#
# OrdinaryDiffEq compiles its integrator, nonlinear solver and factorizations
# for the types of the functions it is handed. A closure over a problem's
# definition and workspace carries the material, boundary-condition and mode
# types, so every new problem type would recompile that whole stack. `solve`
# hands the integrator the right-hand side, the Jacobian and the callbacks
# behind `FunctionWrapper`s with fixed signatures instead: the stack compiles
# once per solver configuration (and is cached in the package image by the
# precompile workloads), and a new problem type compiles only its own kernels.

const _StateVector = Vector{Float64}

const _RHSWrapper =
    FunctionWrapper{Nothing,Tuple{_StateVector,_StateVector,NullParameters,Float64}}
const _JacobianWrapper{M} =
    FunctionWrapper{Nothing,Tuple{M,_StateVector,NullParameters,Float64}}

"""
	ErasedRHS(rhs!)

The ODE right-hand side `rhs!(du, u, p, t)` with a problem-independent type.
Calls with `Float64` vectors and time go through a fixed-signature
`FunctionWrapper`. Any other argument types (dual numbers from a fully
constructed algorithm that differentiates the right-hand side itself) go to
`rhs!` by dynamic dispatch.
"""
struct ErasedRHS
    wrapped::_RHSWrapper
    generic::Any
    ErasedRHS(rhs!) = new(_RHSWrapper(rhs!), rhs!)
end

(f::ErasedRHS)(du::_StateVector, u::_StateVector, p::NullParameters, t::Float64) =
    f.wrapped(du, u, p, t)
(f::ErasedRHS)(du, u, p, t) = f.generic(du, u, p, t)

"""
	ErasedJacobian{M}(jac!)

The ODE Jacobian `jac!(J, u, p, t)` with a problem-independent type, for the
Jacobian matrix type `M` the integrator allocates from the prototype
(`SparseMatrixCSC{Float64,Int}`, or `Matrix{Float64}` without one). A call
with any other matrix type is a `MethodError`: the wrapper would otherwise
convert `J` to a copy and the integrator would never see the entries.
"""
struct ErasedJacobian{M<:AbstractMatrix{Float64}}
    wrapped::_JacobianWrapper{M}
    ErasedJacobian{M}(jac!) where {M<:AbstractMatrix{Float64}} =
        new{M}(_JacobianWrapper{M}(jac!))
end

(f::ErasedJacobian{M})(J::M, u::_StateVector, p::NullParameters, t::Float64) where {M} =
    f.wrapped(J, u, p, t)
