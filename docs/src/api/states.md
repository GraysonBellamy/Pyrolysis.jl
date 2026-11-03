```@meta
CurrentModule = Pyrolysis
```

# State and residual

The state the integrator carries (per-cell component masses and an energy
state; on a moving mesh the node positions and velocities; under
Beer–Lambert absorption the face intensities), the mapping from it to the
quantities the kernels read, and the residual as a pure function of the
state. Reached through `Pyrolysis.States` and `Pyrolysis.Assembly` or
[`Pyrolysis.Internal`](internal.md).

```@docs
States.StateLayout
States.StateSetup
States.DerivedState
States.derive!
States.TemperatureState
States.EnthalpyState
States.trial_temperature
States.block_order
States.mass_index
States.energy_index
States.is_differential
States.GasTransport
States.CaloricExtension
States.positive_part
States.energy_state_value
Assembly.residual!
```
