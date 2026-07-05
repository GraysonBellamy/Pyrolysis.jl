"""
`Structured()` backend asserted regressions.

Two gates that depend on the operator-list assembly and limiter /
cumulative-coupling sparsity widening:

1. **Ramp-ALE regression** — turns a previously prose-only known gap
   into an asserted gate. The deleted legacy hand-coded analytical
   backend recorded
   `max_rel_dev ≈ 0.315` against `DenseAD_ForwardDiff` at the ramp-ALE
   probe state because it linearized simple upwind while the residual
   uses VANLEER. The current `Structured` backend captures the limiter
   chain rule via face-local AD over `_limited_gradient` plus the
   cumulative-`w` indirect chain rule via
   `EmbeddedScaledPrefixSumCoupling`. The strict support-complete check
   should pass: in-support entries agree with dense AD, and there are no
   missing outside-support entries.

2. **Vanleer-kink robustness** — the limiter `ψ(r) = (r + |r|)/(1 + |r|)`
   is non-smooth at `r = 0` (sign change of `grad_0/grad_1` across a
   local extremum). Face-local AD returns a one-sided derivative;
   Newton must still survive. Assert that the Jacobian at a
   temperature-extremum probe state is finite (no NaN / Inf) and that
   the in-support deviation against the dense AD reference is bounded
   (the AD reference has the same one-sided derivative, so they agree).

Both gates run against `Structured(strategy = DirectSolve())` because
that is the verification anchor (densifies structured couplings before
factor / compare).
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Pack a linear T ramp into u, and uniform ξ_vals across cells.
function _set_ramp_state!(u, n_cells, T_lo, T_hi, ξ_vals)
    @inbounds for i = 1:n_cells
        α = (i - 1) / max(1, n_cells - 1)
        u[i] = T_lo * (1 - α) + T_hi * α
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Pack a temperature-with-extremum profile: T monotonically rising, then
# a single dip at cell `i_extremum` (`T[i-1] > T[i] < T[i+1]`). Triggers
# vanleer's `r = 0` kink at the dip.
function _set_extremum_state!(u, n_cells, T_lo, T_hi, T_dip, i_extremum, ξ_vals)
    _set_ramp_state!(u, n_cells, T_lo, T_hi, ξ_vals)
    @inbounds u[i_extremum] = T_dip   # local minimum
    return u
end

@testset "Structured() asserted regressions" begin

    @testset "ramp-ALE limiter chain rule" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        # The probe state recorded when the legacy backend's limiter gap
        # was characterized — T_lo=600, T_hi=1100, ξ=(800.0, 1.5).
        _set_ramp_state!(u, n, 600.0, 1100.0, (800.0, 1.5))

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
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

        @test report.mode == :strict
        @test report.passed
        @test report.max_abs_dev_outside_support ≤ 1e-10

        # Every entry the widened sparsity declares should match the dense AD
        # reference to within rtol. This pins the limiter chain rule plus the
        # indirect cumulative-w coupling.
        ref_scale = max(maximum(abs, report.J_ref), 1e-12)
        in_rel = report.max_abs_dev_in_support / ref_scale
        @test in_rel ≤ 1e-6

        # Hard ceiling: the old hand-coded backend measured
        # `max_rel_dev = 0.315` at this state before it was deleted. The
        # Structured backend's in-support rel-dev must be at least 5 orders
        # of magnitude smaller. (Real value is ~5e-12 in development; the
        # 1e-7 ceiling is the regression guard.)
        @test in_rel ≤ 1e-7
    end

    @testset "Vanleer-kink probe — no NaN, AD agreement" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        # Ramp with a single local minimum at cell n÷2, ξ at moderate
        # pyrolysis-relevant values. The `_limited_gradient` ratio
        # `r = grad_0/grad_1` flips sign across cell `i_extremum`,
        # triggering vanleer's `r = 0` kink.
        _set_extremum_state!(u, n, 600.0, 1100.0, 700.0, n ÷ 2, (800.0, 1.5))

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)

        # Build the J target the same way `verify_jacobian_at` does.
        J = let S = spec.sparsity_phys .| spec.sparsity_geom
            SparseMatrixCSC{Float64,Int}(
                S.m,
                S.n,
                copy(S.colptr),
                copy(S.rowval),
                zeros(Float64, length(S.nzval)),
            )
        end
        Pyrolysis.Jacobian.compute_jacobian!(J, spec, u, def, 0.0, ws)

        # No NaN / Inf in the assembled Jacobian — the kink doesn't break
        # numerics.
        @test all(isfinite, J.nzval)

        # AD agreement on the in-support entries — vanleer's one-sided
        # derivative at `r = 0` is the same on both sides (operator's
        # face-local AD and the dense AD reference both walk the same
        # `_limited_gradient`), so they must agree.
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
        @test report.mode == :strict
        @test report.passed
        ref_scale = max(maximum(abs, report.J_ref), 1e-12)
        in_rel = report.max_abs_dev_in_support / ref_scale
        @test in_rel ≤ 1e-6
    end
end
