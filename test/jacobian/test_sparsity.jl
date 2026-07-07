using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian: build_sparsity, sparsity_profile
using Pyrolysis.Materials: BEER_LAMBERT, SURFACE_ABSORPTION

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# ---------------------------------------------------------------------------
# Builders for WithChi configurations used by the :ale_reduced gating tests.
# These mirror the inverse-analysis shape (lateral shrinkage law triggers
# `WithChi`; non-zero permeability triggers `DarcyFick`; SURFACE_ABSORPTION
# is the canonical radiation model for the cone-calorimeter shape).
# ---------------------------------------------------------------------------
function _build_ale_withchi_darcy_surface_absorption(; n_cells = PERF_N_CELLS)
    solid =
        SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, κ = 1e-14, λ = 1e-5)
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    mat = Material(
        name = :wood_chi,
        components = (solid, gas),
        reactions = (rxn,),
        lateral_shrinkage_law = χ -> 1.0 - 0.2 * χ,
    )
    mesh = create_uniform_mesh(0.01, n_cells, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    def, ws, u = initialize_state(
        problem;
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        radiation_intensity = 50e3,
    )
    return def, ws, u
end

@testset "sparsity source of truth" begin
    for (name, builder, is_ale) in (
        ("baseline_eulerian", build_baseline_eulerian, false),
        ("ale_no_radiation", build_ale_no_radiation, true),
        ("eulerian_surface_absorption", build_eulerian_surface_absorption, false),
        ("ale_beer_lambert_darcy", build_ale_beer_lambert_darcy, true),
    )
        @testset "$name" begin
            def, ws, u, _ = builder()
            S_phys_full, S_geom_full =
                build_sparsity(ws.layout, ws.layout.mode; profile = :full)
            n = length(u)
            @test size(S_phys_full) == (n, n)
            @test size(S_geom_full) == (n, n)
            @test eltype(S_phys_full) == Bool
            @test eltype(S_geom_full) == Bool
            @test nnz(S_phys_full) > 0
            @test (nnz(S_geom_full) > 0) == is_ale

            # DirectSolve gets the :full profile and matches the source-of-truth.
            spec_direct = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
            @test spec_direct.sparsity_phys == S_phys_full
            @test spec_direct.sparsity_geom == S_geom_full

            # ApproxSolve gets the :ale_reduced profile.
            spec_approx = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
            S_phys_red, S_geom_red =
                build_sparsity(ws.layout, ws.layout.mode; profile = :ale_reduced)
            @test spec_approx.sparsity_phys == S_phys_red
            @test spec_approx.sparsity_geom == S_geom_red

            if is_ale
                # ALE: reduced profile drops cumulative widening + entire geom block.
                @test nnz(spec_approx.sparsity_phys) < nnz(spec_direct.sparsity_phys)
                @test nnz(spec_approx.sparsity_geom) == 0
            else
                # Eulerian: nothing to reduce — both profiles produce identical
                # patterns (no ALE block exists to drop).
                @test spec_approx.sparsity_phys == spec_direct.sparsity_phys
                @test nnz(spec_approx.sparsity_geom) == nnz(spec_direct.sparsity_geom) == 0
            end
        end
    end

    @testset "sparsity_profile trait" begin
        @test sparsity_profile(DirectSolve()) === :full
        @test sparsity_profile(ApproxSolve()) === :ale_reduced
        @test sparsity_profile(StructuredSolve()) === :full
        @test sparsity_profile(MatrixFreeSolve()) === :full
    end
end

# ---------------------------------------------------------------------------
# `:ale_reduced` producer-side gating: assert specific sparsity entries are
# dropped that are populated only by `jac.structured` contributions which
# `_emit_jacobian!` itself drops under `:ale_reduced`. Mirrors the
# operator-side gates added in [`docs/jacobian_approxsolve_fix_plan.md`].
# ---------------------------------------------------------------------------
@testset ":ale_reduced specific-entry gating" begin
    @testset "ALE × WithChi × DarcyFick × SurfaceAbsorption" begin
        def, ws, u = _build_ale_withchi_darcy_surface_absorption()
        layout = ws.layout
        mode = layout.mode
        n_cells = layout.n_cells
        NC = Pyrolysis.Problem.n_components(layout)

        S_full, _ = build_sparsity(layout, mode; profile = :full)
        S_red, G_red = build_sparsity(layout, mode; profile = :ale_reduced)

        # Reduced strictly fewer nnz than full (drops geom block, χ-column
        # block, χ-row × z entries, ALE cumulative widening).
        @test nnz(S_red) < nnz(S_full)
        @test nnz(G_red) == 0

        χ_off = Pyrolysis.Problem.chi_offset(layout)
        z_off = Pyrolysis.Problem.z_offset(layout)
        xi_idx(j, i) = layout.n_cells + (j - 1) * layout.n_cells + i

        # Dense χ-column block: T/ξ/z/χ rows × χ cols all present in :full,
        # all absent in :ale_reduced.
        j_probe = 1
        χ_col = χ_off + j_probe
        @test S_full[1, χ_col]                                  # T row
        @test S_full[xi_idx(1, 1), χ_col]                       # ξ row
        @test S_full[z_off+2, χ_col]                          # z row (node 2)
        @test S_full[χ_off+1, χ_col]                          # χ row
        @test !S_red[1, χ_col]
        @test !S_red[xi_idx(1, 1), χ_col]
        @test !S_red[z_off+2, χ_col]
        @test !S_red[χ_off+1, χ_col]

        # χ-row × z entries: present in :full, absent in :ale_reduced (no
        # operator writes these under the reduced profile, so emitting the
        # slot strands a structurally-allocated zero).
        χ_row = χ_off + 1
        @test S_full[χ_row, z_off+1]
        @test S_full[χ_row, z_off+2]
        @test !S_red[χ_row, z_off+1]
        @test !S_red[χ_row, z_off+2]

        # Cell-local χ-row × {T_i, ξ_·,i} entries: present in BOTH profiles
        # (written cell-locally by `ChiOperator` regardless of profile).
        @test S_full[χ_row, 1]
        @test S_red[χ_row, 1]
        for j = 1:NC
            @test S_full[χ_row, xi_idx(j, 1)]
            @test S_red[χ_row, xi_idx(j, 1)]
        end
    end

    @testset "ALE × NoChi × Fickian × BeerLambert (suffix gating)" begin
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        layout = ws.layout
        mode = layout.mode
        n_cells = layout.n_cells
        NC = Pyrolysis.Problem.n_components(layout)

        S_full, _ = build_sparsity(layout, mode; profile = :full)
        S_red, _ = build_sparsity(layout, mode; profile = :ale_reduced)

        xi_idx(j, i) = layout.n_cells + (j - 1) * layout.n_cells + i

        # BeerLambert suffix: T_row[i] × ξ_{j, k} for k > i. To isolate the
        # suffix-only contribution, pick `k` outside the ALE advection
        # 5-cell limiter window `max(1,i-2):min(n,i+2)` — otherwise the
        # limiter emits the same (T_row, ξ_col) entry in both profiles.
        i = 1
        k = n_cells  # well outside i+2, BeerLambert is the only source
        @test k > i + 2
        for j = 1:NC
            @test S_full[i, xi_idx(j, k)]
            @test !S_red[i, xi_idx(j, k)]
        end

        # Diagonal BeerLambert entry T_row[i] × ξ_{j, i}: present in BOTH
        # profiles (written cell-locally by `RadiationOperator`'s diagonal
        # AD pass under any profile).
        for j = 1:NC
            @test S_full[i, xi_idx(j, i)]
            @test S_red[i, xi_idx(j, i)]
        end

        # Kirchhoff surface-transmissivity chain: every energy row
        # couples to the surface cell's TEMPERATURE through ε_eff(T_n, ξ_·,n)
        # in the injected intensity — T_row[i] × T_col[n] slots under :full
        # only (the suffix coupling that fills them is dropped under
        # :ale_reduced, like the rest of the suffix).
        @test S_full[1, n_cells]
        @test S_full[n_cells-2, n_cells]
        @test !S_red[1, n_cells]
    end
end
