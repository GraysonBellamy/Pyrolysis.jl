using Test
using Pyrolysis
using Pyrolysis.Solver:
    DepletionCallbackState,
    is_cell_thin,
    find_thin_cells,
    remove_cell_from_state_vector,
    create_merged_state_vector,
    create_thickness_termination_callback

include(joinpath(@__DIR__, "..", "helpers.jl"))

# Stretched-mesh variant of fixture_problem (per-cell depletion reference)
function fixture_stretched_state(; n_cells, stretch, stretch_from)
    material = fixture_material()
    mesh = create_stretched_mesh(
        0.01,
        n_cells,
        Val(2),
        stretch;
        ξ_initial = (1000.0, 0.0),
        stretch_from = stretch_from,
    )
    bcs = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(mesh, material, bcs, (0.0, 0.02))
    Pyrolysis.Problem.apply_initial_conditions!(problem)
    return fixture_state(problem; use_ale = true)
end

@testset "Surface Depletion Callbacks" begin
    @testset "Thin Cell Detection" begin
        problem = fixture_problem(n_cells = 4)
        _, ws, u, _ = fixture_state(problem; use_ale = true)
        state = DepletionCallbackState(ws.layout, ws.mesh)
        @test length(state.initial_cell_thickness) == ws.layout.n_cells
        @test !is_cell_thin(
            u,
            state.layout,
            1,
            state.initial_cell_thickness[1];
            threshold = 0.05,
        )

        zrng = Pyrolysis.Problem.z_range(ws.layout)
        u_thin = copy(u)
        u_thin[last(zrng)] = u_thin[last(zrng)-1] + 1e-8
        @test is_cell_thin(
            u_thin,
            state.layout,
            state.layout.n_cells,
            state.initial_cell_thickness[state.layout.n_cells];
            threshold = 0.05,
        )
    end

    @testset "Per-Cell Depletion Reference" begin
        # Geometric mesh with max/min cell-size ratio 1.35^11 ≈ 27 — above the
        # 1/threshold = 20 trigger where a single surface-cell reference breaks.
        n = 12
        stretch = 1.35

        # Coarse surface cell (fine cells at the bottom): a surface-cell scalar
        # reference flagged the fine bottom cells at t = 0 and merged away the
        # resolution the stretch provided. Per-cell references flag nothing.
        _, ws, u, _ =
            fixture_stretched_state(n_cells = n, stretch = stretch, stretch_from = :bottom)
        state = DepletionCallbackState(ws.layout, ws.mesh)
        zrng = Pyrolysis.Problem.z_range(ws.layout)
        @test state.initial_cell_thickness ≈ diff(u[zrng])
        @test isempty(find_thin_cells(u, ws.layout, state.initial_cell_thickness, 0.05))

        # Fine surface cell (default :top): a surface-cell scalar reference made
        # the criterion far too lax for coarse cells — a bottom cell at 4% of its
        # OWN initial thickness sat at ≈108% of the fine-surface reference and was
        # never flagged. Per-cell references flag it, and only it.
        _, ws, u, _ =
            fixture_stretched_state(n_cells = n, stretch = stretch, stretch_from = :top)
        state = DepletionCallbackState(ws.layout, ws.mesh)
        zrng = Pyrolysis.Problem.z_range(ws.layout)
        @test isempty(find_thin_cells(u, ws.layout, state.initial_cell_thickness, 0.05))

        Δz1 = state.initial_cell_thickness[1]
        u_thin = copy(u)
        u_thin[zrng[2]] = u_thin[zrng[1]] + 0.04 * Δz1
        @test find_thin_cells(u_thin, ws.layout, state.initial_cell_thickness, 0.05) == [1]

        u_ok = copy(u)
        u_ok[zrng[2]] = u_ok[zrng[1]] + 0.06 * Δz1
        @test isempty(find_thin_cells(u_ok, ws.layout, state.initial_cell_thickness, 0.05))

        # The reference vector must stay in lockstep with the layout: a stale
        # (uncompacted) vector after a merge is a hard error, not a silent skew.
        @test_throws DimensionMismatch find_thin_cells(
            u,
            ws.layout,
            state.initial_cell_thickness[1:(end-1)],
            0.05,
        )
        @test_throws DimensionMismatch DepletionCallbackState(ws.layout, ones(n - 1))
    end

    @testset "State Removal And Merge" begin
        problem = fixture_problem(n_cells = 4)
        _, ws, u, _ = fixture_state(problem; use_ale = true)

        u_removed, layout_removed = remove_cell_from_state_vector(u, ws.layout, 4, 5)
        @test layout_removed.n_cells == 3
        @test layout_removed.n_nodes == 4
        @test length(u_removed) == state_length(layout_removed)

        V = ws.mesh.cell_volumes
        u_merged, layout_merged = create_merged_state_vector(
            u,
            ws.layout,
            4,
            3,
            V[4],
            V[3];
            material = problem.material,
        )
        @test layout_merged.n_cells == 3
        @test length(u_merged) == state_length(layout_merged)
    end

    @testset "Termination Callback Construction" begin
        problem = fixture_problem(n_cells = 3)
        _, ws, _, _ = fixture_state(problem; use_ale = true)
        cb = create_thickness_termination_callback(ws.layout, 1e-6)
        @test cb !== nothing
    end
end
