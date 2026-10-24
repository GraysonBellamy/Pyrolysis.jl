# Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux. Pyrolysis.jl uses a finite-volume discretisation on a
mesh that moves with the solid (Arbitrary Lagrangian–Eulerian), with cell
deactivation at burnout, Arrhenius kinetics, in-depth radiation absorption,
and pressure-driven gas transport (Darcy–Fick).

## Features

- Multi-component materials with temperature-dependent properties
- Single-step and multi-step Arrhenius decomposition kinetics
- Gas generation and transport through porous char (Fickian or Darcy–Fick)
- In-depth radiation absorption via Beer–Lambert (or surface absorption)
- Material volume changes — shrinkage, swelling, intumescence
- Mesh motion with the solid, surface regression and internal interface
  tracking; component masses as the states, so mass is conserved exactly
- Surface burnout by cell deactivation in located events, with an
  equidistribution mesh correction that lets a layer vanish in one event
- Colored sparse automatic-differentiation Jacobian on a pattern built from
  the mesh connectivity, and a stiff Rosenbrock integrator by default
- A conservation recorder and an audit that closes mass and energy
- A Plots recipe extension

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
                 A = 8.5e12, E = 188e3, h = 870e3),
    ),
)

mesh = create_uniform_mesh(0.01, 100)
bcs  = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
)

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 100.0),
)

sol = solve(problem; use_ale = true, radiation_model = SURFACE_ABSORPTION)

# Plot recipes (load `Plots` first)
using Plots
plot(sol)                     # surface temperature trace
plot(sol; kind = :heatmap)    # T(z, t) heatmap
plot(sol; kind = :mlr)        # mass-loss rate
```

The [quick-start tutorial](guide/user_guide/02_quickstart-cone.md) walks
through every line of a run like this one.

## How the documentation is organised

- **[User Guide](guide/user_guide/01_introduction-installation.md)** — how to
  *use* the package: defining materials, reactions, boundary conditions, and
  meshes; building and solving problems; post-processing, burnout and the
  mesh correction, parameter sweeps, sensitivity analysis, and
  troubleshooting.
- **[Technical Reference](guide/technical_reference/01_introduction.md)** —
  the *model itself*: governing equations, constitutive closures,
  discretisation, the moving mesh, cell deactivation, time integration and the
  colored Jacobian, conservation ledgers, and verification/validation — with
  systematic comparisons to ThermaKin2Ds, Gpyro, and the FDS solid-phase model.
- **[Examples](examples.md)** — end-to-end studies (cone calorimeter,
  Douglas fir, pressure-treated wood, inverse analysis, sensitivity) shipped
  in the repository.
- **[API Reference](api/index.md)** — docstrings for the exported surface,
  organised by area.
- **[Development](development.md)** — tests, quality tooling, benchmarks, and
  how these docs are built.
