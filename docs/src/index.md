# Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux. Pyrolysis.jl uses an Arbitrary Lagrangian–Eulerian (ALE)
finite-volume discretisation with adaptive meshing, Arrhenius kinetics,
in-depth radiation absorption, and pressure-driven gas transport (Darcy–Fick).

## Features

- Multi-component materials with temperature-dependent properties
- Single-step and multi-step Arrhenius decomposition kinetics
- Gas generation and transport through porous char (Fickian or Darcy–Fick)
- In-depth radiation absorption via Beer–Lambert (or surface absorption)
- Material volume changes — shrinkage, swelling, intumescence
- ALE mesh motion with surface regression and internal interface tracking
- Adaptive mesh refinement (h- and r-adaptivity)
- Operator-based `Structured` Jacobian with selectable solve strategies,
  a dense AD verifier, and pluggable linear solvers
- Runtime mass- and energy-conservation diagnostics
- Pluggable extensions for Plots recipes, SciML sensitivities, and
  symbolic-derivative verification

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/GraysonBellamy/Pyrolysis.jl")
```

Pyrolysis.jl requires Julia 1.10 or newer.

## Quick start

```julia
using Pyrolysis

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86),
        GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => :MMA;
                 A = 8.5e12, E = 188e3, h = -870e3),
    ),
)

mesh = create_uniform_mesh(0.01, 100, Val(2))
bcs  = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
)

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 100.0),
)

sol = solve(problem; progress = false, use_ale = true,
            radiation_model = SURFACE_ABSORPTION)

# Plot recipes (load `Plots` first)
using Plots
plot(sol)                     # surface temperature trace
plot(sol; kind = :heatmap)    # T(z, t) heatmap
plot(sol; kind = :mlr)        # mass-loss rate
plot(sol; kind = :residual)   # energy-conservation residual
```

The [quick-start tutorial](guide/user_guide/02_quickstart-cone.md) walks
through every line of a run like this one.

## How the documentation is organised

- **[User Guide](guide/user_guide/01_introduction-installation.md)** — how to
  *use* the package: defining materials, reactions, boundary conditions, and
  meshes; building and solving problems; post-processing, adaptivity,
  parameter sweeps, sensitivity analysis, and troubleshooting.
- **[Technical Reference](guide/technical_reference/01_introduction.md)** —
  the *model itself*: governing equations, constitutive closures,
  discretisation, ALE framework, Jacobian architecture, conservation
  diagnostics, and verification/validation — with systematic comparisons to
  ThermaKin2Ds, Gpyro, and the FDS solid-phase model.
- **[Examples](examples.md)** — end-to-end studies (cone calorimeter,
  Douglas fir, pressure-treated wood, inverse analysis, sensitivity) shipped
  in the repository.
- **[API Reference](api/index.md)** — docstrings for the exported surface,
  organised by area.
- **[Development](development.md)** — tests, quality tooling, benchmarks, and
  how these docs are built.
