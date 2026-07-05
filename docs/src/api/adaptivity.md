```@meta
CurrentModule = Pyrolysis
```

# Adaptivity

Error indicators, monitor functions, r-refinement, mesh-quality metrics, and
the surface-depletion machinery. These names live in `Pyrolysis.Adaptivity`
(also re-exported by `Pyrolysis.Internal`); the production surface-depletion
path is enabled through [`solve`](@ref) keywords. See the User Guide's
[adaptive-meshing chapter](../guide/user_guide/10_adaptive-meshing.md) and
the Technical Reference's
[AMR chapter](../guide/technical_reference/14_adaptive-mesh-refinement.md).

```@autodocs
Modules = [Pyrolysis.Adaptivity]
```
