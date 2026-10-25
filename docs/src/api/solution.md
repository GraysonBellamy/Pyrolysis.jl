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
GasDiagnostics
```

## Accessors

```@docs
get_times
get_temperatures
get_concentrations
get_gas_generation_rate
get_pressure
get_partial_pressures
get_intrinsic_gas_densities
get_gas_composition
get_porosity
get_boundary_exchange
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
gas_species_audit
Outputs.ledger_rates
Outputs.live_inventories
```
