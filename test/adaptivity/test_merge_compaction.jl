"""
Regression tests for the χ bookkeeping across depletion merges.

`compact_mesh_after_merge!` / `compact_mesh_after_surface_merge!` must
compact `mesh.m_dry_solid_initial` together with the other per-cell parallel
arrays. Before the fix, an interior merge left every cell i ≥ removed-index
reading the OLD cell i's χ denominator — the merged cell read the removed
cell's stale value (whose mass was also added to the keep cell: a double
count) and every cell above read its lower neighbor's.

The depletion χ merge must weight by the cells' initial dry-solid masses
(χ_i = m_released,i / m_dry_init,i ⇒ Σ χ_i·m_dry_init,i is the extensive
invariant). Volume weighting nearly discarded the depleted cell's pyrolysis
record (V → 0 while m_init is unchanged and its χ is the column maximum),
stepping χ̄ down at each merge and re-inflating A(t).
"""

using Test
using Pyrolysis
using Pyrolysis.Adaptivity:
    merge_adjacent_cells!,
    merge_surface_cell!,
    compact_mesh_after_merge!,
    compact_mesh_after_surface_merge!,
    assert_amr_invariants
using Pyrolysis.Solver: create_merged_state_vector
using Pyrolysis.Problem: chi_range

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "m_dry_solid_initial compaction" begin
    @testset "interior merge + deleteat! compaction" begin
        problem = fixture_problem(n_cells = 6)
        mesh = problem.mesh
        # Distinct values so index shifts are detectable.
        mesh.m_dry_solid_initial .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        total0 = sum(mesh.m_dry_solid_initial)

        # Interior merge: remove cell 3 into keep cell 2 (the production
        # path — find_merge_neighbor prefers the interior neighbor).
        merged_id = merge_adjacent_cells!(mesh, 3, 2; material = problem.material)
        @test merged_id == 2
        compact_mesh_after_merge!(mesh, 3)

        @test length(mesh.m_dry_solid_initial) == length(mesh.cells) == 5
        # Keep cell absorbed the removed cell's denominator extensively…
        @test mesh.m_dry_solid_initial[2] ≈ 2.0 + 3.0
        # …cells above shifted down by one (no stale double-counted slot)…
        @test mesh.m_dry_solid_initial[3:5] ≈ [4.0, 5.0, 6.0]
        @test mesh.m_dry_solid_initial[1] ≈ 1.0
        # …and the column total is conserved.
        @test sum(mesh.m_dry_solid_initial) ≈ total0
        assert_amr_invariants(mesh)
    end

    @testset "surface merge + swap-and-pop compaction" begin
        problem = fixture_problem(n_cells = 6)
        mesh = problem.mesh
        mesh.m_dry_solid_initial .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        total0 = sum(mesh.m_dry_solid_initial)

        merged_id = merge_surface_cell!(mesh; material = problem.material)
        @test merged_id == 5
        compact_mesh_after_surface_merge!(mesh)

        @test length(mesh.m_dry_solid_initial) == length(mesh.cells) == 5
        @test mesh.m_dry_solid_initial ≈ [1.0, 2.0, 3.0, 4.0, 11.0]
        @test sum(mesh.m_dry_solid_initial) ≈ total0
        assert_amr_invariants(mesh)
    end
end

@testset "χ merge is initial-dry-mass-weighted" begin
    problem = fixture_problem(n_cells = 4, lateral = true)
    _, ws, u, _ = fixture_state(problem; use_ale = true)
    layout = ws.layout
    χrng = chi_range(layout)
    @test !isempty(χrng)

    # Depleted surface cell 4 merging into interior cell 3: V_remove → 0
    # (that's why it merges) while its m_init is unchanged and its χ is the
    # column maximum.
    χ_remove, χ_keep = 0.9, 0.1
    u[χrng[4]] = χ_remove
    u[χrng[3]] = χ_keep
    V_remove, V_keep = 1e-9, 2.5e-4
    m_remove, m_keep = 1.0e-3, 1.0e-3

    u_new, layout_new = create_merged_state_vector(
        u,
        layout,
        4,
        3,
        V_remove,
        V_keep;
        material = problem.material,
        m_init_remove = m_remove,
        m_init_keep = m_keep,
    )
    χ_merged = u_new[chi_range(layout_new)[3]]

    # Extensive invariant: χ_m·(m_a + m_b) = χ_a·m_a + χ_b·m_b.
    @test χ_merged * (m_remove + m_keep) ≈ χ_remove * m_remove + χ_keep * m_keep
    @test χ_merged ≈ 0.5
    # Volume weighting would have nearly discarded the depleted cell's
    # pyrolysis record.
    χ_volume_weighted = (V_remove * χ_remove + V_keep * χ_keep) / (V_remove + V_keep)
    @test abs(χ_merged - χ_volume_weighted) > 0.39

    # WithChi merges must be given the pre-merge masses explicitly.
    @test_throws ArgumentError create_merged_state_vector(
        u,
        layout,
        4,
        3,
        V_remove,
        V_keep;
        material = problem.material,
    )
end
