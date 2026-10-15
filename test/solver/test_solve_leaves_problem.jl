"""
`solve` leaves the caller's problem unchanged: the run integrates a copy with
its own mesh. Solving one problem twice therefore repeats the first solve bit
for bit, on ALE (where the first solve moves the nodes) and with depletion
(where it parks cells and moves the exposed boundary face), and a
`parameter_sweep` whose `modify` reuses the base mesh gives every point the
base geometry.
"""

using Test
using Pyrolysis

function _leaves_problem_material(; A = 8.5e12)
    return Material(
        name = :leaves_problem_pmma,
        components = (
            SolidComponent(:pmma, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, λ = 1e-6),
            GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (Reaction(:pyrolysis, 1 => 2, A = A, E = 188e3, h = 870e3),),
    )
end

function _leaves_problem_problem(material = _leaves_problem_material())
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 80e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.003, 20, Val(2)),
        material = material,
        bc_set = bc,
        tspan = (0.0, 90.0),
        dt_initial = 1e-3,
        output_times = collect(0.0:5.0:90.0),
    )
end

const _LEAVES_PROBLEM_KW = (
    radiation_model = SURFACE_ABSORPTION,
    use_ale = true,
    min_thickness = 1e-4,
    progress = false,
)

_leaves_problem_fingerprint(sol) = (
    sol.retcode,
    sol.raw_solution.stats.naccept,
    sol.raw_solution.stats.nf,
    sol.t,
    sol.T,
    sol.ξ,
    sol.z,
    sol.surface_T,
    sol.mass_loss_rate,
    sol.extras.park_times,
)

_leaves_problem_mesh(mesh) = (
    [node.z for node in mesh.nodes],
    copy(mesh.cell_volumes),
    mesh.n_active_cells,
    deepcopy(mesh.boundary_faces),
    mesh.cross_section_area,
)

@testset "solve leaves the problem unchanged" begin
    for (label, kw) in (
        ("ALE", _LEAVES_PROBLEM_KW),
        ("ALE with depletion", (; _LEAVES_PROBLEM_KW..., handle_depletion = true)),
    )
        @testset "$label" begin
            problem = _leaves_problem_problem()
            mesh_before = _leaves_problem_mesh(problem.mesh)
            sol1 = solve(problem; kw...)
            @test sol1.retcode == :Success
            # The run moved its own mesh, and parked cells under depletion.
            @test sol1.extras.thickness_history[end] < 0.95 * 0.003
            if haskey(kw, :handle_depletion)
                @test !isempty(sol1.extras.park_times)
                @test sol1.problem.mesh.n_active_cells < 20
            end
            @test isequal(_leaves_problem_mesh(problem.mesh), mesh_before)
            @test sol1.problem !== problem
            sol2 = solve(problem; kw...)
            @test isequal(
                _leaves_problem_fingerprint(sol2),
                _leaves_problem_fingerprint(sol1),
            )
        end
    end

    @testset "warm start from a parked solve is rejected" begin
        problem = _leaves_problem_problem()
        parked = solve(problem; _LEAVES_PROBLEM_KW..., handle_depletion = true)
        @test_throws ErrorException solve(
            problem;
            _LEAVES_PROBLEM_KW...,
            handle_depletion = true,
            warm_start = parked,
        )
    end

    @testset "parameter_sweep points start from the base geometry" begin
        base = _leaves_problem_problem()
        mesh_before = _leaves_problem_mesh(base.mesh)
        # `modify` reuses the base mesh, as the User Guide's sweep examples do.
        modify(prob, A) = PyrolysisProblem(
            mesh = prob.mesh,
            material = _leaves_problem_material(; A),
            bc_set = prob.bc_set,
            T_initial = prob.T_initial,
            ξ_initial = prob.ξ_initial,
            tspan = prob.tspan,
            dt_initial = prob.dt_initial,
            output_times = prob.output_times,
        )
        grid = [8.5e12, 2.0e13, 8.5e12]
        swept = parameter_sweep(
            base,
            grid;
            modify,
            extract = _leaves_problem_fingerprint,
            parallel = false,
            Base.structdiff(_LEAVES_PROBLEM_KW, (; progress = false))...,
        )
        @test isequal(_leaves_problem_mesh(base.mesh), mesh_before)
        for (A, result) in zip(grid, swept)
            fresh = solve(
                _leaves_problem_problem(_leaves_problem_material(; A));
                _LEAVES_PROBLEM_KW...,
            )
            @test isequal(result, _leaves_problem_fingerprint(fresh))
        end
        @test !isequal(swept[1], swept[2])
    end
end
