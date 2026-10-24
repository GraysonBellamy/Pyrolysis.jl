```@meta
CurrentModule = Pyrolysis
```

# Solving

[`solve`](@ref) integrates a problem in time; its keywords select the moving
mesh, the radiation model, deactivation, the energy state, tolerances, the
integrator policy and the linear solver, and observers. [`setup`](@ref) and
[`integrate`](@ref) split it for custom loops. The User Guide's
[solver chapter](../guide/user_guide/08_running-solver.md) documents every
option with examples, and
[parameter sweeps](../guide/user_guide/12_parameter-sweeps.md) covers the
sweep harness.

```@docs
solve(::PyrolysisProblem{<:Pyrolysis.Meshes.Mesh})
setup
integrate
CoreProblem
parameter_sweep
```

## Integrator policies and the linear solve

```@docs
default_policy
rodas5p_policy
fbdf_policy
kencarp4_policy
StructuredKLU
Integration.DomainCheck
Integration.inadmissible
Integration.consistent_state!
Integration.default_linear_solver
Integration.constraint_residual
Integration.bind_error_test!
Integration.below_absolute_zero
Integration.FAILURE_MEMORY_GROWTH
```

## The mesh correction

```@docs
NoCorrection
EquidistributionCorrection
ActivationTime
Integration.CORRECTION_BAND
```
