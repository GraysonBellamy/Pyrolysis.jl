```@meta
CurrentModule = Pyrolysis
```

# Jacobian

The Jacobian is colored sparse forward-mode automatic differentiation of the
residual on a sparsity pattern built from the mesh connectivity; a dense
forward-mode Jacobian is the reference the pattern must contain and the
values must match. Reached through `Pyrolysis.ColoredJacobian` or
[`Pyrolysis.Internal`](internal.md).

```@docs
ColoredJacobian.jacobian_pattern
ColoredJacobian.ColoredJacobianCache
ColoredJacobian.jacobian!
ColoredJacobian.jac_prototype
ColoredJacobian.dense_jacobian
ColoredJacobian.compare_entries
ColoredJacobian.LateralRows
ColoredJacobian.dense_rows
ColoredJacobian.COLORING_ORDERS
ColoredJacobian.MAX_CHUNK
```
