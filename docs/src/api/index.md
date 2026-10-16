# API Reference

`using Pyrolysis` exposes a deliberately curated top-level surface — the
types and functions needed to define materials, set boundary conditions,
build problems, solve them, and post-process results. Everything else is
reachable as `Pyrolysis.<Submodule>.name` or through the
[`Pyrolysis.Internal`](internal.md) umbrella module, and research-grade
functionality lives in [`Pyrolysis.Experimental`](experimental.md).

The User Guide's [API and parameter reference chapter](../guide/user_guide/16_api-reference.md)
complements these docstring pages with default-value tables and an import
cheat-sheet.

## Pages

| Page | Covers |
|:---|:---|
| [Materials](materials.md) | `Material`, components, phases, mixing rules |
| [Property functions](properties.md) | The seven property-function types |
| [Reactions](reactions.md) | `Reaction`, depletion limiting |
| [Boundary conditions](boundary_conditions.md) | Thermal, mass, and pressure BCs; `BoundaryConditionSet` |
| [Radiation](radiation.md) | `RadiationModel` selection |
| [Geometry](geometry.md) | Mesh creation and mesh types |
| [Problem](problem.md) | `PyrolysisProblem`, `PyrolysisSolution` |
| [Solving](solve.md) | `solve` and `parameter_sweep` |
| [Residual layer](residual.md) | `ProblemDef`/`Workspace`, `residual!`, state accessors |
| [Jacobian](jacobian.md) | Backends, linear-solve strategies, verification |
| [Sensitivity](sensitivity.md) | Forward/adjoint sensitivity entry points |
| [Diagnostics](diagnostics.md) | Conservation trackers and reports |
| [Adaptivity](adaptivity.md) | Error indicators, monitors, mesh operations |
| [Experimental](experimental.md) | P1 radiation, h-AMR |
| [Extensions](extensions.md) | Optional-dependency extensions and their hooks |
| [Internals](internal.md) | The `Internal` umbrella and common internal helpers |
