"""
`Structured()` backend asserted regressions.

Two gates on the ALE operator-list assembly against the dense AD reference:

1. **Ramp-ALE agreement** — a steep temperature ramp with active kinetics
   and mesh motion. The strict support-complete check must pass: every
   in-support entry agrees with dense AD and there is no missing
   outside-support entry. This pins the cumulative-`w` z-row couplings,
   the cell-local dilation rows and the geometry block together.

2. **Extremum-probe robustness** — the same ramp with a single local
   temperature minimum. Face-local AD through the upwinded Darcy–Fick and
   boundary kernels must still produce a finite Jacobian (no NaN / Inf)
   whose in-support deviation against dense AD stays bounded.

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
# a single dip at cell `i_extremum` (`T[i-1] > T[i] < T[i+1]`).
function _set_extremum_state!(u, n_cells, T_lo, T_hi, T_dip, i_extremum, ξ_vals)
    _set_ramp_state!(u, n_cells, T_lo, T_hi, ξ_vals)
    @inbounds u[i_extremum] = T_dip   # local minimum
    return u
end

@testset "Structured() asserted regressions" begin

    @testset "ramp-ALE strict agreement" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        # Probe state: T_lo=600, T_hi=1100, ξ=(800.0, 1.5).
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

        # Every entry the sparsity declares should match the dense AD
        # reference to within rtol. This pins the cumulative-w z-row
        # couplings and the dilation rows.
        ref_scale = max(maximum(abs, report.J_ref), 1e-12)
        in_rel = report.max_abs_dev_in_support / ref_scale
        @test in_rel ≤ 1e-6

        # Hard ceiling (real value ~5e-12; the 1e-7 ceiling is the
        # regression guard).
        @test in_rel ≤ 1e-7
    end

    @testset "extremum probe — no NaN, AD agreement" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        # Ramp with a single local minimum at cell n÷2, ξ at moderate
        # pyrolysis-relevant values.
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
