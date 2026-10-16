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
ApproxSolve
is_approximate
```

## Verification

```@docs
verify_jacobian_at
```

## Operator architecture

```@docs
Operator
StructuredJacobian
StructuredCoupling
PrefixSumCoupling
EmbeddedScaledSuffixSumCoupling
EmbeddedRank1Coupling
```
