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
get_gas_composition
get_boundary_exchange
succeeded
completed
get_profile
interior_events
displaced_residue
```

## Conservation

```@docs
ConservationRecorder
conservation_audit
gas_species_audit
Outputs.ledger_rates
Outputs.live_inventories
```
