# Development

This page covers working on Pyrolysis.jl itself. For using the package, start
with the [User Guide](guide/user_guide/01_introduction-installation.md).

## Repository layout

| Path | Contents |
|:---|:---|
| `src/` | The package, organised into submodules: `Materials`, `BoundaryConditions`, `Geometry`, `Properties`, `Physics`, `Adaptivity`, `Discretization`, `Problem`, `Residual`, `Jacobian`, `Experimental`, `Diagnostics`, `Solver`, plus the `Internal` umbrella. |
| `ext/` | Package extensions (Plots recipes, SciMLSensitivity). |
| `test/` | Test suite mirroring the `src/` module layout; `runtests.jl` walks the directories. |
| `benchmark/` | PkgBenchmark.jl regression suite over the canonical problems. |
| `examples/` | End-to-end studies with their own project environment. |
| `docs/` | The two canonical guide monoliths, the split chapters (`docs/guide/`), and the Documenter site (`docs/make.jl`). |
| `scripts/` | Maintenance scripts (formatting, guide splitting, profiling). |

## Running the tests

```julia
using Pkg
Pkg.test("Pyrolysis")
```

`test/runtests.jl` discovers the domain test directories automatically.
Performance invariants (allocation, type-stability, and method-count checks)
live under `test/perf/` and run as subprocesses; they are also exercised by
the `perf.yml` CI workflow on Julia 1.10 and the latest 1.x.

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
evaluation, Jacobian assembly, and end-to-end solves on the canonical
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
