"""
The `Operator` interface of the `Structured` Jacobian backend.

Each term of the discretization (reactions, conduction, gas transport,
radiation, the boundary conditions, the ALE mesh velocity, χ and the ALE
geometry block) has an `<: Operator` tag type with two methods:

- `apply!` — the term's residual contribution, a thin wrapper over the
  residual kernel it differentiates;
- [`contribute_to_jacobian!`](@ref) — its Jacobian contribution, written
  into the sparse block `jac.A` (or pushed onto `jac.structured` as a
  [`StructuredCoupling`](@ref) when the block is non-local).

The backend (`_compute_jacobian!` in `structured_backend.jl`) calls the
operators in a fixed order. Each `contribute_to_jacobian!` is either
stencil-local ForwardDiff over the residual kernel or an analytical
derivative; the analytical ones are tested against a ForwardDiff reference
of the same kernel at rtol = 1e-12, and every mode is tested against
`DenseAD_ForwardDiff`.
"""

"""
	Operator

Abstract supertype of every Jacobian operator. Operators are stateless
tags: per-problem caches (AD workspaces, scratch buffers) live on the
backend cache and are passed to the operator's methods.
"""
abstract type Operator end

"""
	apply!(out_buffers..., op::Operator, args...) -> nothing

Operator's residual contribution. The buffers it writes and the arguments it
takes are operator-specific (conduction writes `dT`, gas transport writes
`dT` and `dξ`, the ALE mesh-velocity operator fills its workspace).
"""
function apply! end

"""
	contribute_to_jacobian!(jac::StructuredJacobian, op::Operator, u, t,
	                        ws, args...) -> nothing

Operator's contribution to the [`StructuredJacobian`](@ref) at state `u`
and time `t`. `ws` is the operator's workspace from the backend cache; the
remaining arguments are operator-specific. Local entries are added into
`jac.A`; non-local blocks are pushed onto `jac.structured`.
"""
function contribute_to_jacobian! end
