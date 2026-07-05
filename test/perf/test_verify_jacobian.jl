"""
verify_jacobian CI matrix — `Structured` backend × every canonical
configuration × every strategy, all compared against
`DenseAD_ForwardDiff` (verifier reference) at `t = 0` on the structural
support of the active backend's sparsity pattern.

Threshold: `rtol = 1e-6`, `atol = 1e-10`.
"""

using Test
using LinearAlgebra
using SparseArrays
using Pyrolysis
using Pyrolysis: Jacobian
using Pyrolysis.Jacobian:
    JacobianSpec,
    JacobianBackend,
    compute_jacobian!,
    verify_jacobian_at,
    Structured,
    DirectSolve,
    ApproxSolve,
    DenseAD_ForwardDiff

include(joinpath(@__DIR__, "perf_helpers.jl"))

const _CONFIGS = (
    (name = "baseline_eulerian", builder = build_baseline_eulerian, is_ale = false),
    (name = "ale_no_radiation", builder = build_ale_no_radiation, is_ale = true),
    (
        name = "eulerian_surface_absorption",
        builder = build_eulerian_surface_absorption,
        is_ale = false,
    ),
    (
        name = "ale_beer_lambert_darcy",
        builder = build_ale_beer_lambert_darcy,
        is_ale = true,
    ),
)

@testset "verify_jacobian CI matrix" begin
    for cfg in _CONFIGS
        def, ws, u, _ = cfg.builder()

        # Build the verifier reference once per config.
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)

        @testset "$(cfg.name) × DenseAD_ForwardDiff (reference)" begin
            report = verify_jacobian_at(
                ref_spec,
                def,
                ws,
                u,
                0.0;
                reference_spec = ref_spec,
                rtol = 1e-12,
                atol = 1e-14,
            )
            @test report.passed
        end

        @testset "$(cfg.name) × Structured(DirectSolve)" begin
            spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
            report = verify_jacobian_at(
                spec,
                def,
                ws,
                u,
                0.0;
                reference_spec = ref_spec,
                rtol = 1e-6,
                atol = 1e-10,
            )
            @test report.passed
        end

        # ApproxSolve drops the structured couplings — masked compare via
        # the `is_approximate` trait. Only meaningful for ALE configs (the
        # strategy is a no-op on Eulerian where there are no structured
        # couplings to drop).
        if cfg.is_ale
            @testset "$(cfg.name) × Structured(ApproxSolve)" begin
                spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
                report = verify_jacobian_at(
                    spec,
                    def,
                    ws,
                    u,
                    0.0;
                    reference_spec = ref_spec,
                    rtol = 1e-6,
                    atol = 1e-10,
                )
                @test report.passed
            end
        end
    end
end
