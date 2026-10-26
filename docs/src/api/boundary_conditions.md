```@meta
CurrentModule = Pyrolysis
```

# Boundary conditions

Thermal, mass, and pressure boundary conditions, collected per-face into a
[`BoundaryConditionSet`](@ref). Thermal BCs compose additively with `+`; each mass condition supplies
one complete species and energy closure. See the User Guide's
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
GasReservoir
AmbientGasBC
MassTransferFilm
```

## Blowing

A [`ConvectiveBC`](@ref) evaluates its film coefficient with the selected
blowing closure. [`CouetteFilm`](@ref) supports signed exchange, including
ingress and species counterflow; the empirical closures have the support
restrictions stated in their docstrings.

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

## Pressure gradient and hydraulic resistance

```@docs
NeumannPressureBC
ConvectivePressureBC
Physics.ideal_gas_pressure
```

## Container and helpers

```@docs
BoundaryConditionSet
background_gas_concentration
```
