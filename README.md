# Pyrolysis.jl

[![Documentation](https://img.shields.io/badge/docs-pyrolysis.graysonbellamy.dev-blue)](https://pyrolysis.graysonbellamy.dev)

A Julia package for simulating thermal decomposition (pyrolysis) of
materials under high heat flux. Uses an Arbitrary Lagrangian-Eulerian
(ALE) finite-volume discretisation with adaptive meshing, Arrhenius
kinetics, in-depth radiation absorption, and pressure-driven gas
transport (Darcy–Fick).

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
- Pluggable extensions for Plots recipes, SciML sensitivities, and
  symbolic-derivative verification

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

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/GraysonBellamy/Pyrolysis.jl")
```

`Pyrolysis` requires Julia 1.10+.

## Public API

The top-level surface is intentionally curated relative to the pre-refactor
package, while keeping the primary problem-building and solver hooks easy to
reach. Everything else is reachable as `Pyrolysis.<Submodule>.X` or via the
`Pyrolysis.Internal` umbrella module (`using Pyrolysis.Internal`) for
advanced users. Highlights:

| Group | Names |
|-------|-------|
| Materials | `Material`, `SolidComponent`, `GaseousComponent`, `LiquidComponent`, `Reaction`, `ReactionSet`, `MaterialMixing`, `MixingSpec`, `PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN` |
| Properties | `ConstantProperty`, `LinearProperty`, `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`, `StateDependentProperty` |
| Boundary conditions | `HeatFluxBC`, `ConvectiveBC`, `RadiativeBC`, `RadiativeFluxBC`, `CombinedThermalBC`, `DirichletThermalBC`, `AdiabaticBC`, `ImpermeableBC`, `MassFluxBC`, `ConvectiveMassBC`, `DiffusiveMassBC`, `BoundaryConditionSet` |
| Radiation | `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`, `BEER_LAMBERT` |
| Geometry | `create_uniform_mesh`, `create_stretched_mesh` |
| Problem & solve | `PyrolysisProblem`, `PyrolysisSolution`, `solve`, `parameter_sweep` |
| Sensitivities | `forward_sensitivity`, `adjoint_sensitivity` (load `SciMLSensitivity`) |
| Verification | `verify_mixing_derivatives` (load `Symbolics`) |
| Jacobian backends | `JacobianBackend`, `JacobianSpec`, `Structured`, `LinearSolveStrategy`, `DirectSolve`, `StructuredSolve`, `ApproxSolve`, `MatrixFreeSolve`, `compute_jacobian!` |

## Package extensions

Loaded automatically when the corresponding optional dependency is
present in the user's environment:

| Extension | Trigger | Provides |
|-----------|---------|----------|
| `PyrolysisKLUExt` | `using KLU` | `klu_factorization` linear-solver hook |
| `PyrolysisAMGExt` | `using AlgebraicMultigrid` | `amg_preconditioner` |
| `PyrolysisPardisoExt` | `using Pardiso` | `pardiso_factorization` |
| `PyrolysisPlotsExt` | `using RecipesBase` (or `Plots`) | Plot recipes for `PyrolysisSolution` |
| `PyrolysisSciMLSensitivityExt` | `using SciMLSensitivity, Zygote` | Forward / adjoint sensitivity |
| `PyrolysisSymbolicsExt` | `using Symbolics` | Symbolic verification of mixing derivatives (test-only) |

## Parameter sweeps

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

mesh_factory() = create_uniform_mesh(0.01, 100, Val(2))

function modify(prob, p)
    material = Material(
        name = prob.material.name,
        components = prob.material.components,
        reactions = (Reaction(:decomp, :virgin => :MMA;
                              A = p.A, E = p.E, h = -870e3),),
        mixing = prob.material.mixing,
        lateral_shrinkage_law = prob.material.lateral_shrinkage_law,
    )

    PyrolysisProblem(
        mesh = mesh_factory(), # fresh mesh per run
        material = material,
        bc_set = prob.bc_set,
        T_initial = prob.T_initial,
        ξ_initial = prob.ξ_initial,
        tspan = prob.tspan,
        dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end

extract(sol) = (peak_T = maximum(sol.surface_T),
                final_thickness = sol.extras.thickness_history[end])

results = parameter_sweep(problem, grid; modify, extract)
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
- [`docs/archive/`](docs/archive/) — historical design memos retained
  for provenance
