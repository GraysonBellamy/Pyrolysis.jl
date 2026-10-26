# Pyrolysis.jl Benchmarks

PkgBenchmark.jl + BenchmarkTools.jl harness for run-over-run regression
detection.

## Layout

- [`benchmarks.jl`](benchmarks.jl) — entry point. Defines a top-level
  `SUITE::BenchmarkGroup` covering canonical problems (PMMA cone,
  white-pine-style charring, generic charring polymer) × the three hot
  operations (residual, colored Jacobian, end-to-end solve) on the problem
  `setup` builds.
- [`Project.toml`](Project.toml) — pinned `BenchmarkTools` /
  `PkgBenchmark` versions plus `Pyrolysis` itself.

## Running

From the package root:

```julia
using PkgBenchmark
results = benchmarkpkg("Pyrolysis")
export_markdown("benchmark/report.md", results)
```

## Regression detection

```julia
# Compare HEAD against main
using PkgBenchmark
judgment = judge("Pyrolysis", "main")
export_markdown("benchmark/judgment.md", judgment)
```

`judge` flags any benchmark whose mean time changes by more than the
configured threshold (default 5%). PRs that show meaningful regressions
should investigate the offending benchmark in isolation:

```julia
using BenchmarkTools, Pyrolysis
include("benchmark/benchmarks.jl")
run(SUITE["wood_moving_beer_lambert_darcy"]["jacobian"])
```

## CI

Wire `benchmarkpkg` into CI as a non-blocking job (the suite runs in
~2–5 minutes); the allocation and pattern tests under `test/jacobian/` remain
the hard correctness gate.
