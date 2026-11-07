```@meta
CurrentModule = Pyrolysis
```

# Radiation

The radiation model is selected per solve by the `radiation_model` keyword
of [`solve`](@ref): `SURFACE_ABSORPTION` (the default) or `RADIATIVE_TRANSFER`.
Under radiative transfer a [`RadiativeFluxBC`](@ref) arrives with the angular
distribution of its `distribution` keyword.

See the Technical Reference's
[thermal-radiation chapter](../guide/technical_reference/08_thermal-radiation.md)
for the model formulations and validity regimes.

```@docs
RadiationModel
```

## Angular distributions

```@docs
AngularDistribution
Diffuse
Collimated
Cone
```
