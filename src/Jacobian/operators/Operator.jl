"""
The `Operator` interface — the Layer-1 abstraction in the rebuilt
Jacobian architecture.

Each discretization operator (conduction, gas transport, ALE advection,
BC, reactions, …) ships as an `<: Operator` subtype that exposes both
its **residual contribution** (`apply!`) and its **Jacobian
contribution** (`contribute_to_jacobian!`). The default Jacobian path is
**face-local AD** over `apply!` — correct by construction. Hand-rolled
analytical derivatives are an *optimization* layer, not a parallel
codebase that has to be kept in sync with the residual.

Operators declare their dependencies on derived fields (mesh velocity,
T_s, properties) via [`reads`](@ref) / [`writes`](@ref) so the
[`ResidualContext`](@ref) can compose the chain rule unambiguously. The
specific cells / faces an operator touches are declared via
[`stencil`](@ref) so the AD workspace is sized once at problem build
time and reused per Newton step.

# Adding an operator

```julia
struct ConductionOperator <: Operator end

# Required:
apply!(dT, dξ, op::ConductionOperator, state, ctx, t) = ...

# Optional (no-op defaults shown):
reads(::ConductionOperator)  = (:k_eff,)         # consumes thermal conductivity
writes(::ConductionOperator) = ()                 # produces no derived field
stencil(::ConductionOperator) = (cells = -1:1,)  # face-centered tridiagonal

# Default contribute_to_jacobian! uses face-local ForwardDiff over
# apply!. Override only if a hand-rolled analytical derivative is
# benchmark-shown to win and equivalence-tested at rtol=1e-12.
```

The production backend now walks the per-discretization operator list
directly; whole-residual AD is retained only in `DenseAD_ForwardDiff` as a
verifier reference.
"""

"""
	Operator

Abstract supertype of every Jacobian operator. Concrete subtypes
implement the residual contribution (`apply!`) and either inherit the
default AD-based [`contribute_to_jacobian!`](@ref) or override it with
an analytical derivative.

Operators are stateless or near-stateless tags; per-problem caches
(AD workspaces, scratch buffers) belong on the *backend* cache, not on
the operator. This keeps the operator list cheap to construct and
re-walk per Newton step.
"""
abstract type Operator end

"""
	apply!(out_buffers..., op::Operator, state, ctx::ResidualContext, t) -> nothing

Operator's residual contribution. The exact `out_buffers...` shape is
operator-specific (e.g. conduction writes to `dT`, gas transport writes
to both `dT` and `dξ`). Concrete operators override this.

The default fallback errors. There is no sensible default contribution.
"""
function apply!(args...)
    error(
        "Pyrolysis.Jacobian.Operator: no `apply!` method for $(typeof(args[end-3])). " *
        "Each concrete operator must override `apply!`.",
    )
end

"""
	reads(op::Operator) -> Tuple{Vararg{Symbol}}

Derived fields the operator's `apply!` consumes. The framework uses this
to chain-rule via the field's provider when assembling
[`contribute_to_jacobian!`](@ref). Default: empty (operator depends only
on raw state).

Common derived fields:
- `:k_eff` — effective thermal conductivity (mixing rule over components)
- `:rho_cp` — volumetric heat capacity
- `:reaction_rate` — Arrhenius rates
- `:T_s` — boundary surface temperature (implicit Newton solve)
- `:w_cell`, `:w_node` — ALE mesh velocity (cumulative)
- `:theta_A` — column-global strain rate (WithChi only)
"""
reads(::Operator) = ()

"""
	writes(op::Operator) -> Tuple{Vararg{Symbol}}

Derived fields the operator *produces* for downstream consumers. The
[`ResidualContext`](@ref) registers each operator under its `writes`
list so consumers can find their providers. Default: empty.

Operators that write a derived field also typically register a
[`StructuredCoupling`](@ref) via [`contribute_to_structured!`](@ref)
so their consumers' Jacobian contributions can chain-rule through them
without densification.
"""
writes(::Operator) = ()

"""
	stencil(op::Operator) -> NamedTuple

Cell / face neighborhood the operator touches. Used at sparsity-build
time and AD-workspace-sizing time. The shape is operator-specific; the
framework only requires that the values are integer ranges:

```julia
stencil(::ConductionOperator)        = (cells = -1:1,)            # tridiagonal
stencil(::ALEAdvectionOperator)      = (cells = -2:2,)            # 4-cell upwind window
stencil(::ReactionSourceOperator)    = (cells = 0:0,)             # cell-local
stencil(::BCThermalOperator)         = (cells = -1:0, faces = 0:0)
```

Default: cell-local (`(cells = 0:0,)`), the most conservative shape.
Override when the operator pulls in neighbor cells or faces.
"""
stencil(::Operator) = (cells = 0:0,)

"""
	contribute_to_jacobian!(jac::StructuredJacobian, op::Operator,
	                        state, ctx::ResidualContext, t,
	                        ws::OperatorADWorkspace) -> nothing

Operator's contribution to the [`StructuredJacobian`](@ref).
**Default = face-local ForwardDiff over `apply!`** — runs ForwardDiff
on a stencil-local copy of the operator's apply, scatters the
result into `jac.A`, and composes through `reads(op)` via either:

1. Local AD (when the source field's stencil is small and the
   provider has *not* registered a `StructuredCoupling`), or
2. The provider's structured coupling (e.g. `PrefixSumCoupling` for
   the cumulative ALE mesh-velocity), inserted into `jac.structured`
   at assembly time so it isn't materialized into `A`.

Override only if a hand-rolled analytical derivative wins benchmark-
verified equivalence against the AD default at rtol=1e-12.
"""
function contribute_to_jacobian! end

"""
	contribute_to_structured!(jac::StructuredJacobian, op::Operator,
	                          state, ctx::ResidualContext, t) -> nothing

Operators that *produce* a structured coupling (today: only
[`ALEMeshVelocityOperator`](@ref) registers a [`PrefixSumCoupling`](@ref)
and optionally a [`Rank1Coupling`](@ref)) override this to push the
coupling onto `jac.structured`. Default: no-op.
"""
function contribute_to_structured!(jac, op::Operator, state, ctx, t)
    return nothing
end
