```@meta
CurrentModule = Pyrolysis
```

# Deactivation

A cell whose reactive inventory, the mass of every reaction's condensed
reactant, falls below a fraction of its reference is deactivated in an event:
its inventories go to the live cell beside it, on a moving mesh the live cells
above move down a slot, the outermost slot's rows are masked and the boundary
conditions move to the face of the outermost live cell. See the User Guide's
[burnout chapter](../guide/user_guide/10_burnout-correction.md).

```@docs
Deactivation.CellDeactivation
Deactivation.DeactivationEvent
Deactivation.deactivate!
Deactivation.below_threshold
Deactivation.VOID_VOLUME
Deactivation.reactive_mass
Deactivation.threshold_margin
Deactivation.SlotTables
Deactivation.permute_slots!
Deactivation.reset_cells!
Deactivation.replay_events!
```
