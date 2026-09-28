using Test
using Pyrolysis
using Pyrolysis.Solver:
    DepletionCallbackState,
    DT_AFTER_PARK,
    is_cell_thin,
    find_thin_cells,
    create_thickness_termination_callback,
    create_depletion_callback,
    n_live_cells,
    derivative_discontinuity!,
    DiscreteCallback
using Pyrolysis.Problem: z_offset, z_range, get_node_positions

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
    return problem, fixture_state(problem; use_ale = true)...
end

@testset "Depletion detection" begin
    @testset "Thin Cell Detection" begin
        problem = fixture_problem(n_cells = 4)
        _, ws, u, _ = fixture_state(problem; use_ale = true)
        state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
        @test length(state.reference) == ws.layout.n_cells
        @test state.n_live == ws.layout.n_cells
        @test !is_cell_thin(u, ws.layout, 1, state.reference[1]; threshold = 0.05)

        zrng = z_range(ws.layout)
        u_thin = copy(u)
        u_thin[last(zrng)] = u_thin[last(zrng)-1] + 1e-8
        @test is_cell_thin(
            u_thin,
            ws.layout,
            ws.layout.n_cells,
            state.reference[ws.layout.n_cells];
            threshold = 0.05,
        )
    end

    @testset "Per-Cell Depletion Reference" begin
        # Geometric mesh with max/min cell-size ratio 1.35^11 ≈ 27 — above the
        # 1/threshold = 20 ratio at which a single surface-cell reference
        # misclassifies cells.
        n = 12
        stretch = 1.35

        # Coarse surface cell (fine cells at the bottom): a surface-cell scalar
        # reference would flag the fine bottom cells at t = 0. Per-cell
        # references flag nothing.
        problem, _, ws, u, _ =
            fixture_stretched_state(n_cells = n, stretch = stretch, stretch_from = :bottom)
        state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
        zrng = z_range(ws.layout)
        @test state.reference ≈ diff(u[zrng])
        @test isempty(find_thin_cells(u, ws.layout, state.reference, 0.05, n))

        # Fine surface cell (default :top): a bottom cell at 4% of its OWN
        # initial thickness sits at ≈108% of the fine-surface reference. The
        # per-cell reference flags it, and only it.
        problem, _, ws, u, _ =
            fixture_stretched_state(n_cells = n, stretch = stretch, stretch_from = :top)
        state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
        zrng = z_range(ws.layout)
        @test isempty(find_thin_cells(u, ws.layout, state.reference, 0.05, n))

        Δz1 = state.reference[1]
        u_thin = copy(u)
        u_thin[zrng[2]] = u_thin[zrng[1]] + 0.04 * Δz1
        @test find_thin_cells(u_thin, ws.layout, state.reference, 0.05, n) == [1]

        u_ok = copy(u)
        u_ok[zrng[2]] = u_ok[zrng[1]] + 0.06 * Δz1
        @test isempty(find_thin_cells(u_ok, ws.layout, state.reference, 0.05, n))

        # The reference vector has one entry per state slot.
        @test_throws DimensionMismatch find_thin_cells(
            u,
            ws.layout,
            state.reference[1:(end-1)],
            0.05,
            n,
        )
    end

    @testset "Termination callback tracks the live thickness" begin
        problem = fixture_problem(n_cells = 3)
        _, ws, u, _ = fixture_state(problem; use_ale = true)
        cb = create_thickness_termination_callback(ws.layout, 1e-6)
        @test cb !== nothing
        L = 0.01
        @test cb.condition(u, 0.0, nothing) ≈ L - 1e-6

        state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
        cb_live = create_thickness_termination_callback(state, 2e-3)
        @test cb_live.condition(u, 0.0, nothing) ≈ L - 2e-3
        # With one cell parked the live top node is node n_live + 1.
        state.n_live = 2
        z = get_node_positions(u, ws.layout)
        @test cb_live.condition(u, 0.0, nothing) ≈ (z[3] - z[1]) - 2e-3
    end

    @testset "park events give the live count of a saved state" begin
        problem = fixture_problem(n_cells = 5)
        _, ws, u, _ = fixture_state(problem; use_ale = true)
        state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
        push!(state.park_times, 1.0)
        push!(state.park_times, 2.0)
        @test n_live_cells(state, 0.5) == 5
        @test n_live_cells(state, 1.0) == 5     # saved before the event at t = 1
        @test n_live_cells(state, 1.5) == 4
        @test n_live_cells(state, 2.0) == 4
        @test n_live_cells(state, 3.0) == 3
    end

    @testset "the step after a park starts from DT_AFTER_PARK" begin
        # Isothermal zero-order burnout (the V-BO-04 configuration at a higher
        # rate): every cell reaches the threshold at once, so the first
        # cascade parks four of eight cells one step apart.
        ρ_s, T0, L0, r, E = 1000.0, 600.0, 10e-3, 100.0, 5.0e4
        A = r * exp(E / (8.314462618 * T0))
        mat = Material(
            name = :reset_fixture,
            components = (
                SolidComponent(:solid, ρ = ρ_s, c = 1000.0, k = 0.2),
                GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5),
            ),
            reactions = (Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = 0.0),),
        )
        mesh = create_uniform_mesh(L0, 8, Val(2); T_initial = T0, ξ_initial = (ρ_s, 0.0))
        bcs = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        )
        # Cascades at t_k = ln(2^k·20) ρ_s / r: 30.0 s (8 → 4 cells) and 36.9 s
        # (4 → 2); stop between them.
        t_end = 34.0
        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bcs,
            T_initial = z -> T0,
            ξ_initial = [z -> ρ_s, z -> 0.0],
            tspan = (0.0, t_end),
            output_times = [0.0, t_end],
        )
        steps = Tuple{Float64,Float64}[]           # (t, Δt) of every accepted step
        recorder = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                push!(steps, (integrator.t, integrator.t - integrator.tprev))
                derivative_discontinuity!(integrator, false)
            end;
            save_positions = (false, false),
        )
        sol = solve(
            problem;
            use_ale = true,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            callback = recorder,
            progress = false,
        )
        @test sol.retcode == :Success
        parks = sol.extras.park_times
        @test length(parks) == 4
        @test sol.extras.live_cells[end] == 4
        for t_park in parks
            i = findfirst(s -> s[1] > t_park, steps)
            @test i !== nothing
            @test steps[i][2] ≈ DT_AFTER_PARK rtol = 1e-6     # (t + Δt) − t roundoff
        end
        # Steps away from the cascade are far larger than the reset value.
        @test maximum(last, steps) > 1e3 * DT_AFTER_PARK
    end

    @testset "park times are recorded on a fixed mesh" begin
        # Isothermal zero-order burnout on a fixed mesh with the upper half at
        # half density: the upper three cells empty together near 4.75 s and
        # park from the top down; the lower half still holds 30 % of its solid
        # at the end.
        ρ_s, T0, r, E = 1000.0, 600.0, 100.0, 5.0e4
        A = r * exp(E / (8.314462618 * T0))
        mat = Material(
            name = :fixed_park_fixture,
            components = (
                SolidComponent(:solid, ρ = ρ_s, c = 1000.0, k = 0.2, λ = 1e-6),
                GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5),
            ),
            reactions = (Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = 0.0),),
        )
        L0 = 10e-3
        problem = PyrolysisProblem(
            mesh = create_uniform_mesh(
                L0,
                6,
                Val(2);
                T_initial = T0,
                ξ_initial = (ρ_s, 0.0),
            ),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (
                    top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0),
                    bottom = ImpermeableBC(),
                ),
            ),
            T_initial = z -> T0,
            ξ_initial = [z -> z > L0 / 2 ? ρ_s / 2 : ρ_s, z -> 0.0],
            tspan = (0.0, 7.0),
            output_times = [0.0, 7.0],
        )
        sol = solve(
            problem;
            use_ale = false,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            progress = false,
        )
        @test sol.retcode == :Success
        @test sol.extras.live_cells[end] == 3
        @test length(sol.extras.park_times) == 3
        @test issorted(sol.extras.park_times)
        @test all(t -> 0.0 < t < 7.0, sol.extras.park_times)
    end

    @testset "an interior void stops a fixed-mesh solve" begin
        # The same burnout with the lower half at half density: the lower
        # cells empty near 5.0 s beneath a top cell that still holds half its
        # solid, and the solve stops at the first void instead of parking.
        ρ_s, T0, r, E = 1000.0, 600.0, 100.0, 5.0e4
        A = r * exp(E / (8.314462618 * T0))
        mat = Material(
            name = :interior_void_fixture,
            components = (
                SolidComponent(:solid, ρ = ρ_s, c = 1000.0, k = 0.2, λ = 1e-6),
                GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5),
            ),
            reactions = (Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = 0.0),),
        )
        L0 = 10e-3
        problem = PyrolysisProblem(
            mesh = create_uniform_mesh(
                L0,
                6,
                Val(2);
                T_initial = T0,
                ξ_initial = (ρ_s, 0.0),
            ),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (
                    top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0),
                    bottom = ImpermeableBC(),
                ),
            ),
            T_initial = z -> T0,
            ξ_initial = [z -> z < L0 / 2 ? ρ_s / 2 : ρ_s, z -> 0.0],
            tspan = (0.0, 7.0),
            output_times = [0.0, 7.0],
        )
        # The stop falls between output times: read the last accepted step.
        n = 6
        last_step = Ref((t = NaN, ξ_bottom = NaN, ξ_top = NaN))
        recorder = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                u = integrator.u
                last_step[] = (t = integrator.t, ξ_bottom = u[n+1], ξ_top = u[2n])
                derivative_discontinuity!(integrator, false)
            end;
            save_positions = (false, false),
        )
        sol = @test_logs (:warn, r"gasified from within") match_mode = :any solve(
            problem;
            use_ale = false,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            callback = recorder,
            progress = false,
        )
        @test sol.retcode == :Terminated
        @test isempty(sol.extras.park_times)
        @test 4.9 < last_step[].t < 5.2
        @test last_step[].ξ_bottom < 0.1           # cell 1: φ_m = ξ/ρ_s < 1e-4
        @test last_step[].ξ_top > 0.4 * ρ_s        # the top cell is still solid
    end
end
