```@meta
CurrentModule = Pyrolysis
```

# Experimental

`Pyrolysis.Experimental` holds research-grade functionality whose API may
change between releases: the P1 (first-order spherical-harmonics) radiation
model (`P1_QUASI_STEADY` and its solver/boundary-condition machinery) and
experimental h-refinement/coarsening AMR. Prefer the stable top-level exports
when possible.

See the Technical Reference's
[thermal-radiation chapter](../guide/technical_reference/08_thermal-radiation.md)
(P1 formulation, Marshak boundary conditions, validity diagnostics) and
[AMR chapter](../guide/technical_reference/14_adaptive-mesh-refinement.md)
(h-refinement status).

```@autodocs
Modules = [Pyrolysis.Experimental]
```
