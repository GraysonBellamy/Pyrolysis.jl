"""
`RadiationOperator` correctness tests (BEER_LAMBERT branch).

Four gates:

1. **`apply!` no-op modes.** `NO_RADIATION` / `SURFACE_ABSORPTION` /
   `P1_QUASI_STEADY` leave the output unchanged.

2. **`apply!` BEER_LAMBERT.** The operator's Beer-Lambert apply!
   reproduces the existing `compute_beer_lambert_source!` from
   `Physics/radiation.jl` exactly — there is no math duplication; the
   operator wraps the shared physics function. Covered for both
   surface-transmissivity conventions (Kirchhoff ε_eff and explicit-α).

3. **`I_in` cell-edge buffer.** The workspace-aware apply!
   (`apply!(ws::BeerLambertWorkspace, op, …)`) populates
   `ws.q_rad` plus the upstream-cell incident-intensity buffer
   `ws.I_in`. This buffer is what the framework wiring multiplies
   the bare `SuffixSumCoupling` chain rule by (the per-row `−q_rad[i]`
   scaling). On the Kirchhoff path `I_in[n] = ε_eff·I_surface` — the
   flux that actually crosses the surface.

4. **Kirchhoff scaling.** The Kirchhoff path scales the
   whole profile by the surface cell's live ε_eff relative to the
   explicit-α (bare) path.

The bare `SuffixSumCoupling` builder is exercised in
`test_structured_couplings.jl` against a dense reference; its semantic
correctness for Beer-Lambert (the `j > i` non-local optical-depth
chain rule) is composed by the framework's `compose_into!` wiring and
tested there. The full structured-vs-DenseAD
agreement of the Kirchhoff ε_eff chain is in
`test_kirchhoff_eps_eff_jacobian.jl`.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
using Pyrolysis.Physics: compute_beer_lambert_source!
using Pyrolysis.Jacobian: RadiationOperator, build_beer_lambert_workspace, apply!
using Pyrolysis.Problem: populate!
using Pyrolysis.Properties: update_properties!

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# The perf fixture does not run the residual, so its property cache starts
# zeroed (α = ε = 0 would make every Beer-Lambert check trivially 0 ≈ 0).
# Populate the state mirror and evaluate properties first.
function _beer_lambert_fixture()
    def, ws, u, _ = build_ale_beer_lambert_darcy()
    populate!(ws.state_cache, u, ws.layout)
    update_properties!(ws.prop_cache, ws.state_cache, ws.mesh, def.material)
    return def, ws
end

@testset "RadiationOperator" begin

    @testset "no-op for non-Beer-Lambert modes" begin
        def, ws, _, _ = build_baseline_eulerian()
        material = def.material
        dT = ones(ws.layout.n_cells)
        dT_ref = copy(dT)

        # All three "no volumetric source" modes leave dT unchanged.
        for model in (NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY)
            for kirchhoff in (true, false)
                fill!(dT, 1.0)
                apply!(
                    dT,
                    RadiationOperator(),
                    ws.state_cache,
                    ws.mesh,
                    ws.prop_cache,
                    material,
                    50e3,
                    model,
                    kirchhoff,
                )
                @test dT == dT_ref
            end
        end
    end

    @testset "BEER_LAMBERT apply! matches compute_beer_lambert_source! (kirchhoff=$kirchhoff)" for kirchhoff in
                                                                                                   (
        true,
        false,
    )
        # Wood-like material has a non-zero absorption coefficient so
        # the Beer-Lambert path is exercised non-trivially.
        def, ws = _beer_lambert_fixture()
        material = def.material

        # Reference: existing physics path.
        q_ref = zeros(Float64, ws.layout.n_cells)
        compute_beer_lambert_source!(
            q_ref,
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            50e3,
            kirchhoff,
        )

        # Operator path: the bare `dT`-targeting overload.
        q_op = zeros(Float64, ws.layout.n_cells)
        apply!(
            q_op,
            RadiationOperator(),
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            material,
            50e3,
            BEER_LAMBERT,
            kirchhoff,
        )

        @test q_op ≈ q_ref rtol = 1e-12
    end

    @testset "Kirchhoff path scales the bare profile by surface ε_eff" begin
        def, ws = _beer_lambert_fixture()
        N = ws.layout.n_cells

        q_bare = zeros(Float64, N)
        compute_beer_lambert_source!(
            q_bare,
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            50e3,
            false,
        )
        q_kirchhoff = zeros(Float64, N)
        compute_beer_lambert_source!(
            q_kirchhoff,
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            50e3,
            true,
        )

        ε_surf = ws.prop_cache.ε[N]
        @test 0.0 < ε_surf < 1.0
        @test q_kirchhoff ≈ ε_surf .* q_bare rtol = 1e-13
    end

    @testset "BeerLambertWorkspace apply! populates q_rad and I_in" begin
        def, ws = _beer_lambert_fixture()
        material = def.material
        op = RadiationOperator()
        bl_ws = build_beer_lambert_workspace(material, ws.layout.n_cells)
        N = ws.layout.n_cells
        ε_surf = ws.prop_cache.ε[N]

        for kirchhoff in (true, false)
            # Run the workspace-aware apply!.
            apply!(
                bl_ws,
                op,
                ws.state_cache,
                ws.mesh,
                ws.prop_cache,
                material,
                50e3,
                BEER_LAMBERT,
                kirchhoff,
            )

            # 1) q_rad must match the reference exactly.
            q_ref = zeros(Float64, N)
            compute_beer_lambert_source!(
                q_ref,
                ws.state_cache,
                ws.mesh,
                ws.prop_cache,
                50e3,
                kirchhoff,
            )
            @test bl_ws.q_rad ≈ q_ref rtol = 1e-12

            # 2) I_in[i] = intensity entering cell i from above. I_in[N] is the
            # flux that actually crosses the surface: ε_eff·I on the Kirchhoff
            # path, the bare I otherwise.
            I_expected = kirchhoff ? ε_surf * 50e3 : 50e3
            @test bl_ws.I_in[N] ≈ I_expected rtol = 1e-12

            # 3) Per-cell consistency: q_rad[i]·Δz_i = I_in[i]·(1 − e^{−τ_i}).
            A_section = ws.state_cache.A_section[]
            for i = 1:N
                Δz = ws.state_cache.volume[i] / A_section
                τ = ws.prop_cache.α[i] * Δz
                @test bl_ws.q_rad[i] * Δz ≈ bl_ws.I_in[i] * (1 - exp(-τ)) atol = 1e-9
            end
        end

        # 4) Non-Beer-Lambert modes zero out both buffers.
        fill!(bl_ws.q_rad, 1.23)
        fill!(bl_ws.I_in, 4.56)
        apply!(
            bl_ws,
            op,
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            material,
            50e3,
            NO_RADIATION,
            true,
        )
        @test all(iszero, bl_ws.q_rad)
        @test all(iszero, bl_ws.I_in)
    end
end
