"""
Linear-solve strategies for [`StructuredJacobian`](@ref).

A strategy is the rule for inverting `(A + Σⱼ structuredⱼ + chain_rule_from_implicit)`.
Strategies are *orthogonal* to the assembling backend — at solve time, the
caller picks one of:

- [`DirectSolve`](@ref) — densify the structured terms and factor the
  whole sparse `A + dense(structured)` block with KLU/UMFPACK. Simple and
  exact; the verification anchor.
- [`StructuredSolve`](@ref) — sparse-direct on `A` plus fast `apply!`
  on each structured coupling, combined via preconditioned Richardson.
  O(N · n_structured) per iteration.
- [`ApproxSolve`](@ref) — *drop* the structured couplings; solve `A⁻¹ b`
  only. A Schur-style approximation. This strategy inherits the
  [`is_approximate`](@ref) trait, so the [`verify_jacobian_at`](@ref)
  support-complete check uses masked comparison when this strategy is
  in play.
- [`MatrixFreeSolve`](@ref) — solve with a GMRES iteration over the
  `StructuredJacobian` operator action, without densifying structured
  couplings.

Strategies are constructed parameterless today; future iterations may
carry per-strategy knobs (chunk size, Krylov tolerance, …) in fields.

This file ships only the type definitions. The concrete `ldiv!` methods
live alongside the [`StructuredJacobian`](@ref) implementation since they
need access to its internals.
"""

"""
	LinearSolveStrategy

Abstract supertype of every strategy for inverting a
[`StructuredJacobian`](@ref). Concrete strategies override
[`ldiv!`](@ref) on `(out, strat, jac, b)`.
"""
abstract type LinearSolveStrategy end

"""
	DirectSolve <: LinearSolveStrategy

Materialize `A + Σⱼ densify(structuredⱼ)` and factor with sparse-direct
(KLU / UMFPACK). Predictable, robust, the verification anchor for
correctness gates. Use this strategy first when proving a new backend
contributes correctly; switch to `StructuredSolve` once correctness is
established and perf becomes the goal.
"""
struct DirectSolve <: LinearSolveStrategy end

"""
	StructuredSolve <: LinearSolveStrategy

Sparse-direct on `A` plus fast `apply!` on each
[`StructuredCoupling`](@ref), combined via preconditioned Richardson.
Asymptotically `O(N · n_structured)` per iteration; the production
target for ALE problems at large N when exact dense assembly is too
expensive.
"""
struct StructuredSolve <: LinearSolveStrategy end

"""
	ApproxSolve <: LinearSolveStrategy

Drop every [`StructuredCoupling`](@ref) entirely; solve `A⁻¹ b` only.
A Schur-style approximation — Newton sees a Jacobian that omits the dense
ALE C block, takes a few extra iterations, but each linear solve is much
cheaper.

This is the **only** strategy where missing-entry verification is
permitted. The [`is_approximate`](@ref) trait extends to backends
parameterized over this strategy, opting them into masked comparison
in `verify_jacobian_at`.
"""
struct ApproxSolve <: LinearSolveStrategy end

"""
	MatrixFreeSolve <: LinearSolveStrategy

Use a GMRES iteration over the [`StructuredJacobian`](@ref) operator
action. The structured couplings are applied by their fast `apply!`
methods and are not densified into a concrete `A + C` matrix. Useful
for strategy-level experiments and for very large `N` where materializing
the coupled operator is itself expensive. Production callers should still
prefer `DirectSolve()` unless a problem-specific benchmark shows the
matrix-free path wins.
"""
struct MatrixFreeSolve <: LinearSolveStrategy end

"""
	sparsity_profile(strategy::LinearSolveStrategy) -> Symbol

Internal trait that maps a strategy to the sparsity / assembly profile the
`Structured` backend should build for it.

- `:full` — block-tridiagonal physics + Beer-Lambert suffix + χ block + ALE
  geometry block. Structured ALE prefix-sum couplings (z rows) are
  emitted; `GeometryADCache` is built and the `GeometryOperator` runs
  every step. Verifier uses the strict / support-complete comparison.
- `:ale_reduced` — block-tridiagonal physics + diagonal radiation +
  cell-local χ-row deps. Every contribution that would land on
  `jac.structured` is omitted at the producer:

    * z-row × (T, ξ, χ) prefix-sum and rank-1 blocks
      (`ALEMeshVelocityOperator`)
    * z-column geometry block (`S_geom` holds only the node-row diagonal,
      `GeometryADCache` not built)
    * χ-row × χ-col rank-1 (`ChiOperator`)
    * χ-row × z entries (no writer in this profile)
    * Beer-Lambert suffix `j > i` couplings (`RadiationOperator`)

  The masked-compare verifier (`is_approximate(::Structured{ApproxSolve})`)
  expects this set of omissions exactly.

`StructuredSolve` and `MatrixFreeSolve` stay on `:full` for now — their
`ldiv!` paths are not on OrdinaryDiffEq's hot path, so reducing their
assembly buys nothing measurable. Re-evaluate when those are wired into
the integrator's linear-solve protocol.

The trait is **internal**: users continue to write
`Structured(strategy = ApproxSolve())` and the backend reads the profile
during `JacobianSpec` construction.
"""
sparsity_profile(::DirectSolve) = :full
sparsity_profile(::ApproxSolve) = :ale_reduced
sparsity_profile(::StructuredSolve) = :full
sparsity_profile(::MatrixFreeSolve) = :full
