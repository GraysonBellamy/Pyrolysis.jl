```@meta
CurrentModule = Pyrolysis
```

# Internals

```@docs
Internal
```

The User Guide's
[import cheat-sheet](../guide/user_guide/16_api-reference.md) maps common
tasks to the right imports, and its
[custom-loop section](../guide/user_guide/11_solver-jacobian-config.md)
uses [`setup`](@ref) and the residual and Jacobian it builds.

## Submodule exports

```@autodocs
Modules = [Pyrolysis, Pyrolysis.Assembly, Pyrolysis.BoundaryConditions, Pyrolysis.Closures, Pyrolysis.ColoredJacobian, Pyrolysis.Deactivation, Pyrolysis.Integration, Pyrolysis.Kernels, Pyrolysis.Materials, Pyrolysis.Meshes, Pyrolysis.Outputs, Pyrolysis.Physics, Pyrolysis.Problem, Pyrolysis.Properties, Pyrolysis.States]
Private = false
Filter = t -> !(t in Main.API_DOCUMENTED)
```
