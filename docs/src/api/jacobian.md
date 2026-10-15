```@meta
CurrentModule = Pyrolysis
```

# Jacobian backends and linear-solve strategies

The operator-based `Structured` Jacobian assembles per-physics operator
blocks and structured couplings; a linear-solve strategy then decides how the
implicit-integrator linear systems are solved. See the User Guide's
[solver-configuration chapter](../guide/user_guide/11_solver-jacobian-config.md)
and the Technical Reference's
[temporal-integration chapter](../guide/technical_reference/15_temporal-jacobian.md).

## Backends

```@docs
JacobianBackend
JacobianSpec
Structured
DenseAD_ForwardDiff
compute_jacobian!
```

## Linear-solve strategies

```@docs
LinearSolveStrategy
DirectSolve
StructuredSolve
ApproxSolve
MatrixFreeSolve
is_approximate
```

## Verification and benchmarking

```@docs
verify_jacobian_at
benchmark_jacobian_backends
```

## Operator architecture

```@docs
Operator
ResidualContext
StructuredJacobian
StructuredCoupling
PrefixSumCoupling
SuffixSumCoupling
Rank1Coupling
EmbeddedScaledSuffixSumCoupling
EmbeddedRank1Coupling
ImplicitTerm
```
