# Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux. Pyrolysis.jl uses a finite-volume discretisation on a
mesh that moves with the solid (Arbitrary Lagrangian–Eulerian), with cell
deactivation at burnout, Arrhenius kinetics, in-depth radiation absorption,
and lumped or dusty-gas transport through a shared pore state.

## Features

- Multi-component materials with temperature-dependent properties
- Single-step and multi-step Arrhenius decomposition kinetics
- Conservative species transport with lumped transfer or dusty-gas friction
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
    name = :porous_polymer,
    gas_transport = LumpedTransfer(endpoint=LumpedEndpoint(
        λ=1e-6, κ=0., max_skeletal_occupancy=.001)),
    components = (
        SolidComponent(:virgin, ρ = 1190, φ = .3, c = 1420, k = 0.21, ε = 0.86, λ = 1e-6),
        GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => :MMA;
                 A = 8.5e12, E = 188e3, h = 870e3),
    ),
)

mesh = create_uniform_mesh(0.01, 100)
bcs  = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
    mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 100.0),
)

# A mesh that moves with the solid; the cells the surface burns away deactivate
sol = solve(problem; use_ale = true, deactivation = true)

# Plot recipes (load `Plots` first)
using Plots
plot(sol)                     # surface temperature trace
plot(sol; kind = :heatmap)    # T(z, t) heatmap
plot(sol; kind = :mlr)        # mass-loss rate
```

The [quick-start tutorial](guide/user_guide/02_quickstart-cone.md) walks
through every line of a run like this one.

![Surface and back temperatures and actual boundary mass-loss rate from the full PMMA tutorial.](assets/figures/08_quickstart_compact.svg)

Actual results from the [full 1,600-second tutorial](guide/user_guide/02_quickstart-cone.md),
which uses a different material and boundary-condition deck from the short example above.
The full figure also compares gas generation with actual MLR and plots profiles at the current cell positions.

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
