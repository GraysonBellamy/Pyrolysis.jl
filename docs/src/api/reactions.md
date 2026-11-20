```@meta
CurrentModule = Pyrolysis
```

# Reactions and phase changes

Arrhenius decomposition steps and their assembly, and the two-way phase changes of
liquids into their vapours. See the User Guide's
[reactions chapter](../guide/user_guide/04_defining-reactions.md) and the
Technical Reference's
[kinetics chapter](../guide/technical_reference/06_reaction-kinetics.md).

```@docs
Reaction
DepletionLimiter
```

## Phase changes

A liquid held in sorbing solids evaporates and condenses toward sorption equilibrium
with its vapour. The activity laws and the saturation curve are on the
[Materials](materials.md) page.

```@docs
PhaseChange
equilibrium_vapour_pressure
equilibrium_vapour_concentration
Physics.moisture_content
Physics.phase_change_rate
Physics.phase_change_heat
Materials.n_phase_changes
```
