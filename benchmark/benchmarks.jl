"""
PkgBenchmark.jl entry point for Pyrolysis.jl.

Defines a top-level `SUITE::BenchmarkGroup` covering three canonical
problems:

- PMMA cone (fixed mesh, surface absorption)
- White-pine-style charring (moving mesh, Beer–Lambert + lumped and Darcy)
- Generic charring polymer (moving mesh, lumped, no radiation)

For each problem it benchmarks the three hot operations on the problem
`setup` builds:

- One residual evaluation
- One colored Jacobian
- An end-to-end solve over a short time span

Run from the package root with:

```julia
using PkgBenchmark
results = benchmarkpkg("Pyrolysis")
export_markdown("benchmark/report.md", results)
```

For run-over-run regression detection, compare against a baseline:

```julia
judge = judge("Pyrolysis", "main")
```
"""

using BenchmarkTools
using Pyrolysis
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian: jac_prototype, jacobian!

const SUITE = BenchmarkGroup()

include(joinpath(@__DIR__, "decks.jl"))

for (name, config) in pairs(_CONFIGURATIONS)
    group = SUITE[string(name)] = BenchmarkGroup()
    cp = setup(_column(config.material(); config.darcy); config.options...)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    group["residual"] = @benchmarkable residual!($du, $u, $(cp.def), 0.0, $(cp.ws))
    group["jacobian"] = @benchmarkable jacobian!($J, $(cp.jacobian), $u, 0.0)
    problem = _column(config.material(); config.darcy, t_final = 5.0)
    group["solve"] = @benchmarkable solve($problem; $(config.options)...)
end
