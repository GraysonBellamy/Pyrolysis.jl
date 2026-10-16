"""
Assembly strategies for the `Structured` Jacobian backend.

A strategy selects which couplings the assembled Jacobian carries:

- [`DirectSolve`](@ref) — the exact Jacobian: the structured couplings are
  written into the sparse matrix, which the integrator factors with
  KLU/UMFPACK. The verification anchor.
- [`ApproxSolve`](@ref) — *drop* the structured couplings and keep the local
  block only. A Schur-style approximation. This strategy inherits the
  [`is_approximate`](@ref) trait, so the [`verify_jacobian_at`](@ref)
  support-complete check uses masked comparison when it is in play.
"""

"""
	LinearSolveStrategy

Abstract supertype of the `Structured` backend's assembly strategies.
"""
abstract type LinearSolveStrategy end

"""
	DirectSolve <: LinearSolveStrategy

Materialize `A + Σⱼ densify(structuredⱼ)` for a sparse direct factorization
(KLU / UMFPACK). Predictable, robust, the verification anchor for
correctness gates.
"""
struct DirectSolve <: LinearSolveStrategy end

"""
	ApproxSolve <: LinearSolveStrategy

Drop every [`StructuredCoupling`](@ref) entirely and keep the local block `A`.
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

The trait is **internal**: users write `Structured(strategy = ApproxSolve())`
and the backend reads the profile during `JacobianSpec` construction.
"""
sparsity_profile(::DirectSolve) = :full
sparsity_profile(::ApproxSolve) = :ale_reduced
