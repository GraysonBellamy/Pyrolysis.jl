```@meta
CurrentModule = Pyrolysis
```

# Radiation

The radiation model is selected per solve by the `radiation_model` keyword
of [`solve`](@ref): `NO_RADIATION`, `SURFACE_ABSORPTION` or `BEER_LAMBERT`.

See the Technical Reference's
[thermal-radiation chapter](../guide/technical_reference/08_thermal-radiation.md)
for the model formulations and validity regimes.

```@docs
RadiationModel
```
