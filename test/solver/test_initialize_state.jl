"""
Unit tests for `initialize_state` — the canonical PyrolysisProblem →
(ProblemDef, Workspace, u0) setup path (ICs → `to_problem_def` →
`build_workspace` → `create_state_vector`).

Guards the IC-bypass trap that motivated it: `create_state_vector` packs the
mesh as-is, so a hand-rolled composition that skips
`apply_initial_conditions!(problem)` silently integrates from the mesh
constructor's placeholder state, and under a `lateral_shrinkage_law` leaves
`m_dry_solid_initial = 0` (degenerate χ̄) — a resolution-independent error
that still reports `retcode == :Success` (2026-07-06 audit finding).
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: chi_range, mesh_mode, get_option

@testset "initialize_state" begin
    function _init_state_material(; lateral::Bool = false)
        solid = SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, γ = 1.0)
        gas = GaseousComponent(:gas; M = 0.05, c = 1000.0, k = 0.02, λ = 1e-5, γ = 1.0)
        rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
        law = lateral ? (χ -> 1.0 - 0.2 * χ) : nothing
        return Material(
            name = :init_state_fixture,
            components = (solid, gas),
            reactions = (rxn,),
            lateral_shrinkage_law = law,
        )
    end

    # Callable ICs deliberately DIFFER from the mesh constructor's uniform
    # placeholder state so a pack-without-ICs bypass cannot pass.
    function _init_state_problem(; lateral::Bool = false)
        mesh = create_uniform_mesh(
            0.01,
            6,
            Val(2);
            T_initial = 300.0,
            ξ_initial = (1000.0, 0.0),
        )
        bcs =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()))
        return PyrolysisProblem(
            mesh = mesh,
            material = _init_state_material(lateral = lateral),
            bc_set = bcs,
            tspan = (0.0, 1.0),
            T_initial = z -> 300.0 + 5.0e4 * z,
            ξ_initial = [z -> 800.0 + 1.0e4 * z, z -> 0.0],
        )
    end

    @testset "callable ICs land in u0" begin
        problem = _init_state_problem()
        def, ws, u0 = initialize_state(problem)
        n = ws.layout.n_cells
        @test length(u0) == state_length(ws.layout)
        for i = 1:n
            z = problem.mesh.cell_centers[i]
            @test u0[i] ≈ 300.0 + 5.0e4 * z
            @test u0[n+i] ≈ 800.0 + 1.0e4 * z
            @test u0[2n+i] == 0.0
        end
    end

    @testset "kwargs forwarded" begin
        problem = _init_state_problem()
        def, ws, u0 = initialize_state(
            problem;
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            radiation_intensity = 42.0e3,
            min_thickness = 2.0e-6,
        )
        @test mesh_mode(def.mode) isa Pyrolysis.Problem.ALE
        @test get_option(def, :radiation_model, NO_RADIATION) == SURFACE_ABSORPTION
        @test get_option(def, :radiation_intensity, 0.0) == 42.0e3
        @test ws.I_radiation == 42.0e3
        @test ws.min_thickness == 2.0e-6
        # ALE layout carries the node block, packed from the live mesh.
        @test get_node_position(u0, ws.layout, 1) ≈ problem.mesh.nodes[1].z
    end

    @testset "WithChi material gets populated m_dry_solid_initial" begin
        problem = _init_state_problem(lateral = true)
        def, ws, u0 = initialize_state(problem; use_ale = true)
        mesh = problem.mesh
        @test length(mesh.m_dry_solid_initial) == length(mesh.cells)
        @test all(>(0.0), mesh.m_dry_solid_initial)
        # Denominator = V_i · Σ_{j ∈ solid} ξ_j(0) with the PROBLEM's callable
        # IC applied (not the mesh constructor's placeholder).
        for i = 1:length(mesh.cells)
            ξ_solid = 800.0 + 1.0e4 * mesh.cell_centers[i]
            @test mesh.m_dry_solid_initial[i] ≈ ξ_solid * mesh.cell_volumes[i]
        end
        @test !isempty(chi_range(ws.layout))
    end

    @testset "deleted two-arg apply_initial_conditions! stays deleted" begin
        problem = _init_state_problem()
        @test_throws MethodError Pyrolysis.Problem.apply_initial_conditions!(
            problem.mesh,
            problem,
        )
    end
end
