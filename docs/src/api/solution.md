```@meta
CurrentModule = Pyrolysis
```

# Solution

See the User Guide's
[solution chapter](../guide/user_guide/09_solution-plotting.md) for
post-processing.

```@docs
PyrolysisSolution
Outputs.SolutionExtras
SolverStatistics
```

## Accessors

```@docs
get_times
get_temperatures
get_concentrations
get_gas_generation_rate
n_timesteps
succeeded
completed
get_cell_positions
get_profile
interpolate_solution
```

## Conservation

```@docs
ConservationRecorder
conservation_audit
Outputs.ledger_rates
Outputs.live_inventories
```
