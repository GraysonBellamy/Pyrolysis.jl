"""
`Structured` backend correctness gate.

Asserts that `JacobianSpec(Structured(), …)` produces a Jacobian that
agrees with `DenseAD_ForwardDiff` under the **support-complete strict
check** at every canonical config × off-default state. This is the
gate that prevents the rebuilt architecture from landing as a path
with hidden missing entries.

The backend walks the per-discretization operator list and must agree with
the whole-residual `DenseAD_ForwardDiff` verifier on every canonical
configuration.

# Strategy
- `Structured(strategy = DirectSolve())` is the default and the
  correctness anchor; tested across every config × state.
- `Structured(strategy = ApproxSolve())` is tagged approximate, so
  `verify_jacobian_at` falls back to masked compare; the ApproxSolve
  variant is exercised on a single ALE config to confirm the trait
  plumbing is correct.
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian:
    verify_jacobian_at,
    JacobianSpec,
    DenseAD_ForwardDiff,
    Structured,
    DirectSolve,
    ApproxSolve,
    is_approximate

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Mirror the off-default state setter used elsewhere in the Jacobian tests.
function _set_state!(u, n_cells, T_val, ξ_vals)
    @inbounds for i = 1:n_cells
        u[i] = T_val
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

function _set_ramp_state!(u, n_cells, T_lo, T_hi, ξ_vals)
    @inbounds for i = 1:n_cells
        u[i] = T_lo + (T_hi - T_lo) * (i - 1) / max(n_cells - 1, 1)
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Canonical configurations × off-default states.
const _PROBES = (
    (
        name = "eulerian PMMA hot",
        builder = build_baseline_eulerian,
        T = 1000.0,
        ξ = (1100.0, 0.0),
    ),
    (
        name = "eulerian PMMA hot+gas",
        builder = build_baseline_eulerian,
        T = 1200.0,
        ξ = (400.0, 8.0),
    ),
    (
        name = "eulerian surface-absorption hot",
        builder = build_eulerian_surface_absorption,
        T = 1000.0,
        ξ = (700.0, 5.0),
    ),
    (
        name = "ale no-radiation hot",
        builder = build_ale_no_radiation,
        T = 1000.0,
        ξ = (1100.0, 0.0),
    ),
    (
        name = "ale Beer-Lambert+Darcy hot",
        builder = build_ale_beer_lambert_darcy,
        T = 1000.0,
        ξ = (400.0, 0.0),
    ),
    (
        name = "ale Beer-Lambert+Darcy hot+gas",
        builder = build_ale_beer_lambert_darcy,
        T = 1000.0,
        ξ = (400.0, 5.0),
    ),
)

@testset "Structured backend correctness gate" begin

    @testset "trait — is_approximate" begin
        # DirectSolve is exact; ApproxSolve is approximate.
        @test is_approximate(Structured(strategy = DirectSolve())) == false
        @test is_approximate(Structured(strategy = ApproxSolve())) == true
        # The default is DirectSolve.
        @test is_approximate(Structured()) == false
    end

    @testset "support-complete: $(p.name)" for p in _PROBES
        def, ws, u, _ = p.builder()
        n = ws.layout.n_cells
        _set_state!(u, n, p.T, p.ξ)

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(), def, ws)
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
    end

    @testset "support-complete: ALE ramp limiter geometry" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        _set_ramp_state!(u, n, 600.0, 1100.0, (700.0, 1.5))

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(), def, ws)
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
    end

    @testset "ApproxSolve trait routes to masked compare" begin
        # ALE config where the structured-coupling list can carry a
        # non-trivial PrefixSumCoupling. With an empty list, ApproxSolve
        # and DirectSolve produce the same matrix — but the trait flips
        # :strict to :masked, and that's what we're asserting here.
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (1100.0, 0.0))

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
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
        @test report.mode == :masked
        @test report.passed
    end

    # ---------------------------------------------------------------------
    # `:ale_reduced` producer-side gating regression coverage. Each round-
    # trip checks that under `Structured{ApproxSolve}`:
    #   1. Both `DirectSolve` and `ApproxSolve` JacobianSpecs build cleanly.
    #   2. The masked-compare verifier passes — i.e. for every entry that
    #      `ApproxSolve` emits a sparsity slot for, the value matches
    #      `DenseAD_ForwardDiff`. The pre-fix correctness leak (`(χ_row,
    #      z_col)` slots emitted but no operator wrote them under
    #      `:ale_reduced`) would surface here as a masked-compare failure.
    # ---------------------------------------------------------------------
    @testset "ApproxSolve verifier round-trip — ALE WithChi" begin
        # Small (5-cell) ALE WithChi configuration; exercises the
        # ApproxSolve verifier round-trip regression.
        solid = SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            κ = 1e-14,
            λ = 1e-5,
        )
        gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        mat = Material(
            name = :wood_chi,
            components = (solid, gas),
            reactions = (rxn,),
            lateral_shrinkage_law = χ -> 1.0 - 0.2 * χ,
        )
        mesh = create_uniform_mesh(0.01, 5, Val(2))
        bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        def, ws, u = initialize_state(
            problem;
            use_ale = true,
            radiation_model = Pyrolysis.Materials.SURFACE_ABSORPTION,
            radiation_intensity = 50e3,
        )
        _set_state!(u, ws.layout.n_cells, 1000.0, (400.0, 0.0))
        # Non-degenerate χ block: with χ (or the m_dry_solid_initial
        # denominator) left at zero this fixture sat at the χ̄ = 0 identity
        # point and the round-trip never exercised the shrinkage coupling
        # (2026-07-06 audit finding, fixed by the initialize_state migration).
        for (k, idx) in enumerate(Pyrolysis.Problem.chi_range(ws.layout))
            u[idx] = 0.2 + 0.05 * k
        end
        @test column_chi_bar_from_state(u, ws.layout, ws.mesh) > 0

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
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
        @test report.mode == :masked
        @test report.passed
    end

    @testset "ApproxSolve verifier round-trip — ALE BeerLambert" begin
        # Mirrors `build_ale_beer_lambert_darcy()` but with a verifier
        # round-trip under `ApproxSolve`. Catches anyone reintroducing the
        # suffix sparsity drop without the operator-side gate (or vice
        # versa).
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (400.0, 5.0))

        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
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
        @test report.mode == :masked
        @test report.passed
    end
end
