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

## Pressure (Darcy mode)

```@docs
DirichletPressureBC
NeumannPressureBC
ConvectivePressureBC
```

## Container and helpers

```@docs
BoundaryConditionSet
background_gas_concentration
```
