"""
PkgBenchmark.jl entry point for Pyrolysis.jl.

Defines a top-level `SUITE::BenchmarkGroup` covering the canonical
problems from the Phase 9 plan:

- PMMA cone (Eulerian, surface absorption)
- White-pine-style charring (ALE, Beer-Lambert + Darcy-Fick)
- Generic charring polymer (ALE, Fickian)

For each problem we benchmark the three hot operations:

- Single residual evaluation
- Single Jacobian assembly
- End-to-end solve over a short tspan

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
using Pyrolysis: residual!
using Pyrolysis.Jacobian: JacobianSpec, Structured, DirectSolve, compute_jacobian!

# Reuse the canonical perf-test helpers — they already build every
# `SimulationMode` axis combination as `(def, ws, u, du)` tuples ready
# for the residual.
include(joinpath(@__DIR__, "..", "test", "perf", "perf_helpers.jl"))

const SUITE = BenchmarkGroup()

# ---------------------------------------------------------------------------
# Helper: register the three hot ops for a built (def, ws, u, du) tuple.
# ---------------------------------------------------------------------------
function _register_hot_ops!(group, label, builder)
    def, ws, u, du = builder()

    group["residual"] = @benchmarkable residual!($du, $u, $def, 0.0, $ws)

    # Build the structured Jacobian skeleton once; re-use across samples.
    jac_spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    J = Pyrolysis.Jacobian._alloc_target(jac_spec, length(u))
    compute_jacobian!(J, jac_spec, u, def, 0.0, ws)

    group["jacobian"] =
        @benchmarkable compute_jacobian!($J, $jac_spec, $u, $def, 0.0, $ws) setup=(fill!(
            $J,
            0.0,
        ))

    return nothing
end

# ---------------------------------------------------------------------------
# Helper: benchmark a short end-to-end solve.
# ---------------------------------------------------------------------------
function _bench_solve(material, bc_set, mesh; use_ale, radiation_model, t_final = 5.0)
    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        tspan = (0.0, t_final),
        dt_initial = 1e-3,
        output_times = [0.0, t_final / 2, t_final],
    )
    return @benchmarkable solve(
        $problem;
        progress = false,
        use_ale = $use_ale,
        radiation_model = $radiation_model,
    )
end

# ---------------------------------------------------------------------------
# PMMA cone — Eulerian, surface absorption
# ---------------------------------------------------------------------------
SUITE["pmma_eulerian_surface"] = BenchmarkGroup()
_register_hot_ops!(
    SUITE["pmma_eulerian_surface"],
    "PMMA Eulerian SA",
    build_eulerian_surface_absorption,
)

# ---------------------------------------------------------------------------
# White-pine-style charring — ALE, Beer-Lambert + Darcy-Fick
# ---------------------------------------------------------------------------
SUITE["wood_ale_beer_lambert_darcy"] = BenchmarkGroup()
_register_hot_ops!(
    SUITE["wood_ale_beer_lambert_darcy"],
    "Wood ALE BL+Darcy",
    build_ale_beer_lambert_darcy,
)

# ---------------------------------------------------------------------------
# Generic charring polymer — ALE, Fickian, no radiation
# ---------------------------------------------------------------------------
SUITE["polymer_ale_fickian"] = BenchmarkGroup()
_register_hot_ops!(
    SUITE["polymer_ale_fickian"],
    "Polymer ALE Fickian",
    build_ale_no_radiation,
)

# ---------------------------------------------------------------------------
# End-to-end solves (one per canonical configuration)
# ---------------------------------------------------------------------------
SUITE["solve"] = BenchmarkGroup()

let mesh =
        create_uniform_mesh(5e-3, 20, Val(2); T_initial = 300.0, ξ_initial = (1190.0, 0.0)),
    mat = Material(
        name = :PMMA_bench,
        components = (
            SolidComponent(:virgin, ρ = 1190.0, c = 1500.0, k = 0.21, ε = 0.96),
            GaseousComponent(:MMA, M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2, A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0),
        ),
    ),
    bc = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(25e3) + RadiativeBC(ε = 0.96, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )

    SUITE["solve"]["pmma_eulerian"] =
        _bench_solve(mat, bc, mesh; use_ale = false, radiation_model = SURFACE_ABSORPTION)
end
