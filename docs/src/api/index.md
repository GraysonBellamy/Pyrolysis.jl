# API Reference

`using Pyrolysis` exposes a deliberately curated top-level surface — the
types and functions needed to define materials, set boundary conditions,
build meshes and problems, solve them, and post-process results. Everything
else is reachable as `Pyrolysis.<Submodule>.name` or through the
[`Pyrolysis.Internal`](internal.md) umbrella module.

The User Guide's [API and parameter reference chapter](../guide/user_guide/16_api-reference.md)
complements these docstring pages with default-value tables and an import
cheat-sheet.

## Pages

| Page | Covers |
|:---|:---|
| [Materials](materials.md) | `Material`, components, phases, mixing rules |
| [Property functions](properties.md) | The eight property-function types |
| [Reactions](reactions.md) | `Reaction`, the depletion limiter |
| [Boundary conditions](boundary_conditions.md) | Thermal, mass, and pressure conditions; blowing; `BoundaryConditionSet` |
| [Radiation](radiation.md) | `RadiationModel` selection |
| [Meshes](meshes.md) | The mesh constructors, coordinate systems, topology and geometry |
| [Problem](problem.md) | `PyrolysisProblem` and the initial-condition profiles |
| [Solving](solve.md) | `solve`, `setup`, `integrate`, integrator policies, the mesh correction, `parameter_sweep` |
| [Solution](solution.md) | `PyrolysisSolution`, its accessors, the conservation recorder and audit |
| [State and residual](states.md) | The state layout, the state mapping, the residual |
| [Kernels and closures](kernels.md) | The physics kernels and the mesh closure |
| [Deactivation](deactivation.md) | Cell deactivation and its events |
| [Jacobian](jacobian.md) | The colored sparse Jacobian and its dense reference |
| [Extensions](extensions.md) | Optional-dependency extensions |
| [Internals](internal.md) | The `Internal` umbrella |
