using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian:
    build_sparsity, sparsity_profile, materializes_geometry_block, structural_zero_diagonal
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
                # ALE: the reduced profile drops the geometry block (keeping
                # only the structural diagonal of the node rows). The physics
                # pattern loses only the Beer–Lambert suffix when present;
                # the ALE dilation is cell-local and survives both profiles.
                if name == "ale_beer_lambert_darcy"
                    @test nnz(spec_approx.sparsity_phys) < nnz(spec_direct.sparsity_phys)
                else
                    @test spec_approx.sparsity_phys == spec_direct.sparsity_phys
                end
                z_rows = (Pyrolysis.Problem.z_offset(ws.layout)+1):n
                @test nnz(spec_approx.sparsity_geom) == length(z_rows)
                @test all(r -> spec_approx.sparsity_geom[r, r], z_rows)
            else
                # Eulerian: nothing to reduce — both profiles produce identical
                # patterns (no ALE block exists to drop).
                @test spec_approx.sparsity_phys == spec_direct.sparsity_phys
                @test nnz(spec_approx.sparsity_geom) == nnz(spec_direct.sparsity_geom) == 0
            end
        end
    end

    @testset "Beer–Lambert energy rows reach every node (ALE, :full)" begin
        for builder in
            (build_ale_beer_lambert_darcy_moderate, build_ale_beer_lambert_fickian_moderate)
            def, ws, u, _ = builder()
            layout = ws.layout
            n_cells = layout.n_cells
            z_off = Pyrolysis.Problem.z_offset(layout)
            z_top = z_off + layout.n_nodes
            _, G_full = build_sparsity(layout, layout.mode; profile = :full)
            _, G_red = build_sparsity(layout, layout.mode; profile = :ale_reduced)

            # The optical depth above every cell depends on the top node.
            @test all(i -> G_full[i, z_top], 1:n_cells)
            # Species rows keep the bounded reach 1..i+3.
            ξ_row = n_cells + 1
            @test G_full[ξ_row, z_off+4]
            @test !G_full[ξ_row, z_top]
            # The reduced profile keeps only the node rows' diagonal.
            @test nnz(G_red) == layout.n_nodes
            @test all(k -> G_red[z_off+k, z_off+k], 1:layout.n_nodes)
        end

        # Without Beer–Lambert the bottom cell's energy row stops at node 4.
        def, ws, u, _ = build_ale_no_radiation()
        z_off = Pyrolysis.Problem.z_offset(ws.layout)
        _, G = build_sparsity(ws.layout, ws.layout.mode; profile = :full)
        @test G[1, z_off+4]
        @test !G[1, z_off+5]
    end

    @testset "concurrent cache misses" begin
        # Threaded sweeps build their patterns at once; every caller must get
        # the one cached pattern for its key, and the cache must keep them all.
        _, ws = build_ale_no_radiation()[1:2]
        mode = ws.layout.mode
        NC = Pyrolysis.Problem.n_components(ws.layout)
        layouts = [Pyrolysis.Problem.StateLayout(NC, n, n + 1, mode) for n = 4:19]
        expected = [build_sparsity(l, mode; profile = :ale_reduced) for l in layouts]
        cache = Pyrolysis.Jacobian._SPARSITY_CACHE
        n_bad = 0
        for _ = 1:200
            lock(() -> empty!(cache), Pyrolysis.Jacobian._SPARSITY_CACHE_LOCK)
            tasks = [
                Threads.@spawn build_sparsity(l, mode; profile = :ale_reduced) for
                l in layouts
            ]
            got = fetch.(tasks)
            cached = [build_sparsity(l, mode; profile = :ale_reduced) for l in layouts]
            ok =
                got == expected &&
                all(got[k] === cached[k] for k in eachindex(layouts)) &&
                length(cache) == length(layouts)
            n_bad += !ok
        end
        @test n_bad == 0
    end

    @testset "sparsity_profile trait" begin
        @test sparsity_profile(DirectSolve()) === :full
        @test sparsity_profile(ApproxSolve()) === :ale_reduced
    end
end

# ---------------------------------------------------------------------------
# Every profile stores every row's diagonal, so the integrator's shifted
# W = J − M/(γΔt) has the pattern of J. Where the reduced profile has no
# writer (node rows, χ rows) the diagonal is a structural zero, so the extra
# slots leave J's values unchanged.
# ---------------------------------------------------------------------------
@testset "every row stores its diagonal" begin
    builders = (
        ("baseline_eulerian", build_baseline_eulerian),
        ("ale_no_radiation", build_ale_no_radiation),
        ("ale_beer_lambert_darcy_moderate", build_ale_beer_lambert_darcy_moderate),
        ("ale_withchi", _build_ale_withchi_darcy_surface_absorption),
    )
    for (name, builder) in builders, strategy in (DirectSolve(), ApproxSolve())
        @testset "$name, $(nameof(typeof(strategy)))" begin
            def, ws, u = builder()[1:3]
            spec = JacobianSpec(Structured(; strategy), def, ws)
            S = spec.sparsity_phys .| spec.sparsity_geom
            @test all(r -> S[r, r], 1:size(S, 1))
            layout = ws.layout
            profile = sparsity_profile(strategy)
            zero_rows = structural_zero_diagonal(layout, layout.mode, profile)
            if strategy isa ApproxSolve
                # The declared structural zeros are exactly the rows past the
                # T/ξ blocks (node rows, χ rows), and they stay zero.
                n_phys = layout.n_cells * (1 + Pyrolysis.Problem.n_components(layout))
                @test zero_rows == (n_phys+1):length(u)
                mixed_reacting_state!(u, layout.n_cells)
                J = Pyrolysis.Jacobian._alloc_target(spec, length(u))
                compute_jacobian!(J, spec, u, def, 0.0, ws)
                @test all(r -> J[r, r] == 0.0, zero_rows)
            else
                @test isempty(zero_rows)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# `materializes_geometry_block` (the default-factorization rule) must agree
# with the pattern: it holds exactly when the geometry pattern carries more
# than the node rows' structural diagonal.
# ---------------------------------------------------------------------------
function _build_eulerian_beer_lambert_darcy()
    mat = _wood_like_material(α = 0.5)
    mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    return _finalize(
        problem;
        use_ale = false,
        radiation_model = BEER_LAMBERT,
        radiation_intensity = 50e3,
    )
end

@testset "geometry-block rule matches the pattern" begin
    builders = (
        ("eulerian_surface_absorption", build_eulerian_surface_absorption),
        ("eulerian_beer_lambert_darcy", _build_eulerian_beer_lambert_darcy),
        ("ale_no_radiation", build_ale_no_radiation),
        ("ale_beer_lambert_darcy_moderate", build_ale_beer_lambert_darcy_moderate),
        ("ale_withchi", _build_ale_withchi_darcy_surface_absorption),
    )
    for (name, builder) in builders, profile in (:full, :ale_reduced)
        @testset "$name, $profile" begin
            def, ws = builder()[1:2]
            layout = ws.layout
            _, S_geom = build_sparsity(layout, layout.mode; profile)
            n_diagonal = length(Pyrolysis.Problem.z_range(layout))
            @test materializes_geometry_block(layout.mode, profile) ==
                  (nnz(S_geom) > n_diagonal)
        end
    end
end

# ---------------------------------------------------------------------------
# The reduced profile keeps every cell within the one-cell stencil of the
# block-tridiagonal physics (the widest local stencil), so KLU sees a band of
# width one cell.
# ---------------------------------------------------------------------------
function _max_cell_reach(S, layout)
    n = layout.n_cells
    z_off = n * (1 + Pyrolysis.Problem.n_components(layout))
    χ = Pyrolysis.Problem.chi_range(layout)
    chi_off = isempty(χ) ? typemax(Int) : first(χ) - 1
    cell(r) = r <= z_off ? mod1(r, n) : (r <= chi_off ? r - z_off : r - chi_off)
    return maximum(
        abs(cell(rowvals(S)[k]) - cell(c)) for c = 1:size(S, 2) for k in nzrange(S, c)
    )
end

@testset "reduced ALE pattern stays within the physics stencil" begin
    for builder in (
        build_ale_no_radiation,
        build_ale_beer_lambert_darcy_moderate,
        _build_ale_withchi_darcy_surface_absorption,
    )
        def, ws = builder()[1:2]
        layout = ws.layout
        S_phys, S_geom = build_sparsity(layout, layout.mode; profile = :ale_reduced)
        @test _max_cell_reach(S_phys .| S_geom, layout) == 1
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
        # block, χ-row × z entries); the geometry
        # pattern keeps only the node rows' structural diagonal.
        @test nnz(S_red) < nnz(S_full)
        @test nnz(G_red) == layout.n_nodes

        χ_off = Pyrolysis.Problem.chi_offset(layout)
        z_off = Pyrolysis.Problem.z_offset(layout)
        xi_idx(j, i) = layout.n_cells + (j - 1) * layout.n_cells + i

        # Dense χ-column block: T/ξ/z/χ rows × χ cols all present in :full,
        # all absent in :ale_reduced except the χ rows' structural diagonal.
        j_probe = 1
        χ_col = χ_off + j_probe
        @test S_full[1, χ_col]                                  # T row
        @test S_full[xi_idx(1, 1), χ_col]                       # ξ row
        @test S_full[z_off+2, χ_col]                          # z row (node 2)
        @test S_full[χ_off+2, χ_col]                          # χ row
        @test !S_red[1, χ_col]
        @test !S_red[xi_idx(1, 1), χ_col]
        @test !S_red[z_off+2, χ_col]
        @test !S_red[χ_off+2, χ_col]
        @test all(i -> S_red[χ_off+i, χ_off+i], 1:n_cells)

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

        # BeerLambert suffix: T_row[i] × ξ_{j, k} for k > i. Pick `k` outside
        # the block-tridiagonal physics stencil so the suffix is the only
        # source of the entry.
        i = 1
        k = n_cells
        @test k > i + 1
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
