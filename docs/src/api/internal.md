```@meta
CurrentModule = Pyrolysis
```

# Internals

`Pyrolysis.Internal` is an umbrella module that re-exports the public surface
of every submodule for advanced users:

```julia
using Pyrolysis
using Pyrolysis.Internal   # brings in Geometry accessors, Solver helpers, …
```

These names are **not** part of the stable public API and may change between
releases. The User Guide's
[import cheat-sheet](../guide/user_guide/16_api-reference.md) maps common
tasks to the right imports.

## Solution post-processing (`Pyrolysis.Solver`)

```@docs
Solver.get_profile
Solver.get_cell_positions
Solver.interpolate_solution
```

## Boundary temperatures (`Pyrolysis.Discretization`)

```@docs
Discretization.get_surface_temperature
Discretization.get_back_face_temperature
```

## Solver callbacks and iterative solvers (`Pyrolysis.Solver`)

```@docs
Solver.create_thickness_termination_callback
Solver.create_depletion_callback
Solver.ILUGMRESFactorization
Solver.ILUBiCGSTABFactorization
```
