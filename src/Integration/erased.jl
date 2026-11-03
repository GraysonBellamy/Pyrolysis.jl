# ============================================================================
# Functions the integrator calls, behind problem-independent types
# ============================================================================
#
# OrdinaryDiffEq compiles its integrator, nonlinear solver and factorizations
# for the types of the functions it is handed. A closure over a problem's
# definition and workspace carries the material and boundary-condition types,
# so every new problem would recompile that whole stack. The right-hand side
# and the Jacobian therefore reach the integrator behind `FunctionWrapper`s
# with fixed signatures.

const _StateVector = Vector{Float64}

const _RHSWrapper =
    FunctionWrapper{Nothing,Tuple{_StateVector,_StateVector,NullParameters,Float64}}
const _JacobianWrapper{M} =
    FunctionWrapper{Nothing,Tuple{M,_StateVector,NullParameters,Float64}}

"""
	ErasedRHS(rhs!)

The right-hand side `rhs!(du, u, p, t)` with a problem-independent type. Calls
with `Float64` vectors and time go through a fixed-signature `FunctionWrapper`;
any other argument types go to `rhs!` by dynamic dispatch.
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

The Jacobian `jac!(J, u, p, t)` with a problem-independent type, for the
matrix type `M` the integrator allocates from the prototype. A call with any
other matrix type is a `MethodError`: a conversion would hand the integrator
a copy it never sees.
"""
struct ErasedJacobian{M<:AbstractMatrix{Float64}}
    wrapped::_JacobianWrapper{M}
    ErasedJacobian{M}(jac!) where {M<:AbstractMatrix{Float64}} =
        new{M}(_JacobianWrapper{M}(jac!))
end

(f::ErasedJacobian{M})(J::M, u::_StateVector, p::NullParameters, t::Float64) where {M} =
    f.wrapped(J, u, p, t)
