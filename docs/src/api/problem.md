```@meta
CurrentModule = Pyrolysis
```

# Problem

A [`PyrolysisProblem`](@ref) holds the mesh, the material, the boundary
conditions, the initial conditions as functions of position, the time span
and the output times. See the User Guide's
[problem-building chapter](../guide/user_guide/07_building-problem-modes.md).

```@docs
PyrolysisProblem
Problem.initial_fields
Problem.step_profile
Problem.linear_profile
Problem.gaussian_profile
Problem.tanh_profile
Problem.validate_problem
Problem.validate_darcy_initial_porosity
Problem.warn_degenerate_gas_transport
```
