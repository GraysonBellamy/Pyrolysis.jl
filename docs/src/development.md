# Development

This page covers working on Pyrolysis.jl itself. For using the package, start
with the [User Guide](guide/user_guide/01_introduction-installation.md).

## Repository layout

| Path | Contents |
|:---|:---|
| `src/` | The package, organised into submodules in load order: `Materials`, `Properties`, `BoundaryConditions`, `Physics` (the shared physics layers), `Meshes`, `States`, `Kernels`, `Closures`, `Assembly`, `Deactivation`, `ColoredJacobian` (the discretization), `Problem`, `Outputs`, `Integration`, plus the `Internal` umbrella. |
| `ext/` | Package extension (Plots recipes). |
| `test/` | Test suite mirroring the `src/` module layout (a directory per module, plus `jacobian/`, `multid/`, `type_analysis/` and `code_scan/`), the verification cases (`verification/`), and measurement scripts (`gates/`); `runtests.jl` walks the directories. |
| `benchmark/` | PkgBenchmark.jl regression suite over the canonical problems. |
| `examples/` | End-to-end studies with their own project environment. |
| `docs/` | The two canonical guide monoliths, the split chapters (`docs/guide/`), and the Documenter site (`docs/make.jl`). |
| `scripts/` | Maintenance scripts (formatting, guide splitting, profiling). |

## Running the tests

```julia
using Pkg
Pkg.test("Pyrolysis")
```

`test/runtests.jl` runs every `test_*.jl` file under the module directories,
each in a module of its own, the type analysis in a fresh process, and the
verification cases (`test/verification/runcases.jl`). The module suites test
both energy forms themselves; `PYROLYSIS_ENERGY_STATE=enthalpy` runs the
verification cases in the enthalpy form. `PYROLYSIS_TEST_GROUP` selects one
part (`jacobian`, `burnout`, `other`, `domains`, or `all`, the default) and
`PYROLYSIS_TEST_FAST=1` skips the verification cases' convergence spot
checks; the `tests.yml` workflow runs the parts side by side on every pull
request to `main`, and the whole suite weekly and on dispatch.

## Quality tooling

Formatting and static-quality checks run in isolated environments (the CI
`quality.yml` workflow runs the same commands):

```bash
julia --project=scripts scripts/format.jl --check   # JuliaFormatter (gating)
julia --project=test/quality test/quality/aqua.jl   # Aqua package hygiene
julia --project=test/quality test/quality/jet.jl    # JET analysis (advisory)
```

Run `scripts/format.jl` without `--check` to apply formatting.

## Benchmarks

`benchmark/benchmarks.jl` defines a PkgBenchmark.jl `SUITE` covering residual
evaluation, the colored Jacobian, and end-to-end solves on the canonical
configurations (PMMA cone, white pine, charring polymer). See
[`benchmark/README.md`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/benchmark/README.md)
for `benchmarkpkg`/`judge` usage.

## Working on the guides

The two monoliths — `docs/User_Guide.md` and
`docs/Technical_Reference_Guide.md` — are the canonical, hand-edited
documents. Two derived views are generated from them:

1. `docs/guide/**` (GitHub-browsable chapters) — regenerate after editing a
   monolith with:
   ```bash
   julia scripts/split_guides.jl
   ```
2. `docs/src/guide/**` (Documenter site sources) — generated automatically by
   the docs build; never edit these by hand. The generator
   (`docs/generate_guide.jl`) also converts `$`-delimited math to Documenter
   syntax and fails the build if anything is left unconverted.

## Building the documentation site

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The site is written to `docs/build/`. On pushes to `main`, the
`docs.yml` workflow builds the site and deploys it to Cloudflare Pages at
[pyrolysis.graysonbellamy.dev](https://pyrolysis.graysonbellamy.dev); pull
requests build the site for validation without deploying.
