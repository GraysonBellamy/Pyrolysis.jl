```@meta
CurrentModule = Pyrolysis
```

# Materials

Components combine into a [`Material`](@ref) via mixing rules. See the User
Guide's [materials chapter](../guide/user_guide/03_defining-materials.md) for
the full walkthrough.

## Components

```@docs
Material
SolidComponent
LiquidComponent
GaseousComponent
```

## Phase queries

```@docs
phase
molar_mass
is_gas
is_solid
is_liquid
```

## Mixing rules

The `MixingRule` enum values are `PARALLEL`, `SERIES`, `WEIGHTED`,
`BRUGGEMAN`, and `CARMAN_KOZENY` (permeability only).

```@docs
MixingRule
MixingSpec
MaterialMixing
```

## Gas transport

```@docs
LumpedTransfer
DustyGas
LumpedEndpoint
DustyEndpoint
```
