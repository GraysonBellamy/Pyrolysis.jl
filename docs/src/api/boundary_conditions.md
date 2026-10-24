```@meta
CurrentModule = Pyrolysis
```

# Boundary conditions

Thermal, mass, and pressure boundary conditions, collected per-face into a
[`BoundaryConditionSet`](@ref). Thermal and mass BCs compose additively with
`+`. See the User Guide's
[boundary-conditions chapter](../guide/user_guide/05_boundary-initial-conditions.md).

## Thermal

```@docs
AdiabaticBC
HeatFluxBC
ConvectiveBC
RadiativeBC
RadiativeFluxBC
DirichletThermalBC
CombinedThermalBC
```

## Mass (species)

```@docs
ImpermeableBC
MassFluxBC
ConvectiveMassBC
DiffusiveMassBC
CombinedMassBC
```

## Blowing

A [`ConvectiveBC`](@ref) reduces its film coefficient by a blowing model when
gas leaves through its face.

```@docs
AbstractBlowingModel
NoBlowing
CouetteFilm
PropertyRatioFilm
PoolingFilm
FunctionBlowing
GasSpecies
NITROGEN
nasa7_heat_capacity
eucken_viscosity
fuller_diffusivity
BlowingFace
blowing_state
blowing_ratio
blowing_factor
```

## Pressure (Darcy–Fick transport)

```@docs
DirichletPressureBC
NeumannPressureBC
ConvectivePressureBC
Physics.ideal_gas_pressure
```

## Container and helpers

```@docs
BoundaryConditionSet
background_gas_concentration
```
