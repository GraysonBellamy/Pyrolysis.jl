```@meta
CurrentModule = Pyrolysis
```

# Radiation

The radiation model is selected per-solve via the `radiation_model` keyword
of [`solve`](@ref). The stable values of the `RadiationModel` enum are
`NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT`; the experimental
P1 model (`P1_QUASI_STEADY`) is exported from
[`Pyrolysis.Experimental`](experimental.md).

See the Technical Reference's
[thermal-radiation chapter](../guide/technical_reference/08_thermal-radiation.md)
for the model formulations and validity regimes.

```@docs
RadiationModel
```
