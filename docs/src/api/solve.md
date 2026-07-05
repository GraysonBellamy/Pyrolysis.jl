```@meta
CurrentModule = Pyrolysis
```

# Solving

The [`solve`](@ref) entry point drives the implicit time integration; its
keywords select the simulation modes (ALE, radiation model, Darcy transport),
tolerances, Jacobian backend, callbacks, and diagnostics. The User Guide's
[solver chapter](../guide/user_guide/08_running-solver.md) documents every
option with examples, and
[parameter sweeps](../guide/user_guide/12_parameter-sweeps.md) covers the
sweep harness.

```@docs
solve(::PyrolysisProblem)
parameter_sweep
```
