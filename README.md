# Pyrolysis.jl

[![Documentation](https://img.shields.io/badge/docs-pyrolysis.graysonbellamy.dev-blue)](https://pyrolysis.graysonbellamy.dev)

A Julia package for simulating thermal decomposition (pyrolysis) of
materials under high heat flux. Uses a finite-volume discretisation on a mesh
that moves with the solid (Arbitrary Lagrangian–Eulerian), with cell
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
  the mesh connectivity, a stiff Rosenbrock integrator by default, and a
  conservation recorder and audit
- A Plots recipe extension

## Quick start

```julia
using Pyrolysis

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86, λ = 1e-6),
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
    mass    = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
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

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/GraysonBellamy/Pyrolysis.jl")
```

`Pyrolysis` requires Julia 1.10+.

## Public API

The top-level surface exports the problem-building and solver names.
Everything else is reachable as `Pyrolysis.<Submodule>.X` or via the
`Pyrolysis.Internal` umbrella module (`using Pyrolysis.Internal`) for
advanced users. Highlights:

| Group | Names |
|-------|-------|
| Materials | `Material`, `SolidComponent`, `GaseousComponent`, `LiquidComponent`, `Reaction`, `MaterialMixing`, `MixingSpec`, `PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN` |
| Properties | `ConstantProperty`, `LinearProperty`, `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`, `SutherlandProperty`, `StateDependentProperty` |
| Boundary conditions | `HeatFluxBC`, `ConvectiveBC`, `RadiativeBC`, `RadiativeFluxBC`, `CombinedThermalBC`, `DirichletThermalBC`, `AdiabaticBC`, `ImpermeableBC`, `MassFluxBC`, `ConvectiveMassBC`, `DiffusiveMassBC`, `BoundaryConditionSet` |
| Radiation | `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`, `BEER_LAMBERT` |
| Mesh | `create_uniform_mesh`, `create_stretched_mesh`, `create_mesh_from_nodes`, `Planar`, `Cylindrical`, `Spherical` |
| Problem & solve | `PyrolysisProblem`, `solve`, `setup`, `integrate`, `parameter_sweep` |
| Solver configuration | `default_policy`, `rodas5p_policy`, `fbdf_policy`, `kencarp4_policy`, `StructuredKLU`, `EquidistributionCorrection`, `ActivationTime`, `NoCorrection` |
| Solution | `PyrolysisSolution`, `get_times`, `get_temperatures`, `get_concentrations`, `get_profile`, `get_cell_positions`, `succeeded`, `completed`, `ConservationRecorder`, `conservation_audit` |

## Package extensions

Loaded automatically when the corresponding optional dependency is
present in the user's environment:

| Extension | Trigger | Provides |
|-----------|---------|----------|
| `PyrolysisPlotsExt` | `using RecipesBase` (or `Plots`) | Plot recipes for `PyrolysisSolution` |

## Parameter sweeps

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

@everywhere function modify(prob, p)
    m = prob.material
    material = Material(
        name = m.name,
        components = m.components,
        reactions = (Reaction(:decomp, :virgin => :MMA; A = p.A, E = p.E, h = 870e3),),
        mixing = m.mixing,
    )
    # A solve never changes the problem, so the points share its mesh.
    return PyrolysisProblem(
        mesh = prob.mesh,
        material = material,
        bc_set = prob.bc_set,
        T_initial = prob.T_initial,
        ξ_initial = prob.ξ_initial,
        tspan = prob.tspan,
        output_times = prob.output_times,
    )
end

@everywhere extract(sol) = (peak_T = maximum(sol.surface_T),
                            final_thickness = sol.extras.thickness_history[end])

results = parameter_sweep(problem, grid; modify, extract, use_ale = true)
```

## Examples

End-to-end studies under [`examples/`](examples/) cover:

- Cone calorimeter (PMMA, both Beer–Lambert and surface-absorption
  radiation models)
- Douglas fir cone calorimeter
- Pressure-treated wood with moisture
- Darcy–Fick gas-transport demonstration
- Inverse analysis of white-pine kinetics
- Sensitivity analysis (Latin Hypercube + Sobol indices)
- Literature comparison and validation

## Benchmarks

`benchmark/benchmarks.jl` is a [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl)
suite covering the canonical configurations (PMMA cone, white-pine,
charring polymer) for run-over-run regression detection. See
[`benchmark/README.md`](benchmark/README.md).

## Documentation

The full documentation — user guide, technical reference, examples, and API
reference — is hosted at
**[pyrolysis.graysonbellamy.dev](https://pyrolysis.graysonbellamy.dev)**.

The canonical sources live under [`docs/`](docs/):

- [`docs/User_Guide.md`](docs/User_Guide.md) and
  [`docs/Technical_Reference_Guide.md`](docs/Technical_Reference_Guide.md) —
  the hand-edited monolithic guides (browse per-chapter under
  [`docs/guide/`](docs/guide/))
- [`docs/make.jl`](docs/make.jl) — the Documenter.jl site build (see the
  site's Development page for build instructions)
