```@meta
CurrentModule = Pyrolysis
```

# Internals

`Pyrolysis.Internal` is an umbrella module that re-exports the public
surface of every submodule (a name two submodules export as different
objects is left to its submodule):

```julia
using Pyrolysis
using Pyrolysis.Internal   # the residual, the colored Jacobian, the state layout, …
```

These names are **not** part of the stable public API and may change between
releases. The User Guide's
[import cheat-sheet](../guide/user_guide/16_api-reference.md) maps common
tasks to the right imports, and its
[custom-loop section](../guide/user_guide/11_solver-jacobian-config.md)
uses [`setup`](@ref) and the residual and Jacobian it builds.
