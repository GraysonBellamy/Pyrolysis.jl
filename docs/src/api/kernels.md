```@meta
CurrentModule = Pyrolysis
```

# Kernels and closures

Each kernel reads the derived state, the mesh topology and metrics, the
material and its boundary conditions, and adds face fluxes and cell sources
to a contributions buffer; the mesh closure sets the node velocities.
Reached through `Pyrolysis.Kernels` and `Pyrolysis.Closures` or
[`Pyrolysis.Internal`](internal.md).

## Kernels

```@docs
Kernels.ReactionKernel
Kernels.ConductionKernel
Kernels.GasTransportKernel
Kernels.BeerLambertKernel
Kernels.ThermalBoundaryKernel
Kernels.RelativeFluxKernel
Kernels.BoundaryPatch
Kernels.exposed_face
Kernels.van_albada
```

## Closures

```@docs
Closures.BalanceRowClosure
Closures.PotentialClosure
Closures.LayerCorrection
Closures.LateralShrinkage
Closures.Correction
Closures.PrescribedCorrection
Closures.bind_correction
Closures.residues
Closures.vanishes
Closures.consistent_multipliers
Closures.correction!
Closures.rank_check
```
