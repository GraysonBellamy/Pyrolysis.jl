```@meta
CurrentModule = Pyrolysis
```

# Deactivation

A cell whose condensed mass falls below a fraction of its initial
condensed mass is deactivated in an event: its inventories go to the live
cell beneath, its rows are masked and the boundary conditions move to the
face of the outermost live cell. See the User Guide's
[burnout chapter](../guide/user_guide/10_burnout-correction.md).

```@docs
Deactivation.CellDeactivation
Deactivation.DeactivationEvent
Deactivation.deactivate!
Deactivation.threshold_margin
Deactivation.ConservativeTransfer
```
