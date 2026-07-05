"""
CI invariant — bounded-allocation Jacobian assembly (Structured/DirectSolve).

The structured-Jacobian hot path was de-allocated in 2026-06 (typed zeros in
the Properties helpers, mutable AD workspaces, `Val(NC)` closure captures,
concretely-typed `StructuredOpsCache` fields, stack-Dual ALE strain-rate
sweep): baseline_eulerian went from 167.5 KiB/eval at 18% GC to ~5 KiB/eval
at <1% GC, and the ALE configs from ~3 MiB/eval to ~100-280 KiB/eval.
A second pass removed the remaining per-eval churn (heterogeneous
component-tuple indexing in the ALE dilation block, copy-on-construct
structured couplings — they now alias workspace buffers refilled each
assembly, per-cell `ForwardDiff.gradient` in the Beer-Lambert α pass →
stack-Dual sweep, runtime-`NC` `ntuple`s in the advection kernels, and
per-boundary BC closures union-typed by branch-merged BC lookups →
function barriers): ~1.2 KiB/eval Eulerian, ~4 KiB/eval ALE configs.

The Eulerian configs now measure 0 B/eval; the ALE budgets are NOT zero
because each assembly pushes freshly-boxed coupling structs onto
`jac.structured` (the array payloads alias workspace buffers, but the
immutable struct wrappers are heap-allocated per push) — ~2.1 KiB/eval.

The point of the gate is to catch *regressions to per-cell allocation* —
those show up as tens-to-hundreds of KiB, far above these budgets.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian: JacobianSpec, Structured, DirectSolve, compute_jacobian!

include(joinpath(@__DIR__, "perf_helpers.jl"))

# Per-eval byte budgets (measured 0 B Eulerian, ~2.1 KiB ALE, ~2.2 KiB
# ALE+BeerLambert+Darcy on julia 1.12.6; generous absolute headroom for
# version-to-version codegen drift).
const JAC_ALLOC_BUDGETS = (
    baseline_eulerian = 2 * 1024,
    eulerian_surface_absorption = 2 * 1024,
    baseline_eulerian_darcy = 2 * 1024,
    ale_no_radiation = 6 * 1024,
    ale_beer_lambert_darcy = 6 * 1024,
)

function _check_jacobian_alloc(builder, label, budget)
    @testset "$label" begin
        def, ws, u, _ = builder()
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        J = Pyrolysis.Jacobian._alloc_target(spec, length(u))
        compute_jacobian!(J, spec, u, def, 0.0, ws)
        compute_jacobian!(J, spec, u, def, 0.0, ws)
        bytes = @allocated compute_jacobian!(J, spec, u, def, 0.0, ws)
        @test bytes <= budget
    end
end

@testset "bounded-allocation Jacobian" begin
    _check_jacobian_alloc(
        build_baseline_eulerian,
        "baseline {Eulerian, NoRadiation, Fickian, NoChi}",
        JAC_ALLOC_BUDGETS.baseline_eulerian,
    )
    _check_jacobian_alloc(
        build_eulerian_surface_absorption,
        "{Eulerian, SurfaceAbsorption, Fickian, NoChi}",
        JAC_ALLOC_BUDGETS.eulerian_surface_absorption,
    )
    _check_jacobian_alloc(
        build_baseline_eulerian_darcy,
        "{Eulerian, NoRadiation, DarcyFick, NoChi}",
        JAC_ALLOC_BUDGETS.baseline_eulerian_darcy,
    )
    _check_jacobian_alloc(
        build_ale_no_radiation,
        "{ALE, NoRadiation, Fickian, NoChi}",
        JAC_ALLOC_BUDGETS.ale_no_radiation,
    )
    _check_jacobian_alloc(
        build_ale_beer_lambert_darcy,
        "{ALE, BeerLambert, DarcyFick, NoChi}",
        JAC_ALLOC_BUDGETS.ale_beer_lambert_darcy,
    )
end
