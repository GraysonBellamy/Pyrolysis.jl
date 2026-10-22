using Test
using Distributed
using Pyrolysis
using Pyrolysis.Integration

# The deck and the sweep's functions live in `Main` on every process, where
# the workers can deserialize the closures a problem carries.
const _SWEEP = quote
    using Pyrolysis
    _sweep_problem(q) = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 8),
        material = Material(;
            name = :inert,
            components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2),),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, 10.0),
        output_times = [0.0, 5.0, 10.0],
    )
    _sweep_modify(problem, q) = _sweep_problem(q)
    _sweep_extract(sol) = (T_s = sol.surface_T[end], retcode = sol.retcode)
end
Core.eval(Main, _SWEEP)

const _GRID = [1e3, 5e3, 1e4]

@testset "parameter sweep on the finite-volume core: serial" begin
    base = Main._sweep_problem(1e3)
    sweep = parameter_sweep(
        base,
        _GRID;
        modify = Main._sweep_modify,
        extract = Main._sweep_extract,
        parallel = false,
    )
    @test sweep ==
          [Main._sweep_extract(Integration.solve(Main._sweep_problem(q))) for q in _GRID]
    @test all(r -> r.retcode === :Success, sweep)
    @test issorted([r.T_s for r in sweep])
    # The solve keywords reach every solve; the default extract is the solution.
    whole = parameter_sweep(base, _GRID[1:1]; modify = Main._sweep_modify, reltol = 1e-8)
    @test only(whole) isa Pyrolysis.Outputs.CoreSolution
    @test only(whole).surface_T ≈ Integration.solve(base; reltol = 1e-8).surface_T
    @test_throws ArgumentError parameter_sweep(base, _GRID; parallel = false)
end

# An ablating solid heated at the top, its pre-exponential factor `A` the
# swept parameter: on a moving mesh its cells thin and deactivate.
function _ablator_material(A)
    return Material(
        name = :ablator,
        components = (
            SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 1.0, λ = 1e-6),
            GaseousComponent(:gas, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-4),
        ),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = 1.0e5, h = 1.0e6),),
    )
end

_ablator_problem(A) = PyrolysisProblem(
    mesh = create_uniform_mesh(3e-3, 8),
    material = _ablator_material(A),
    bc_set = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(2.8e6), bottom = AdiabaticBC()),
        mass = (top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 300.0 + 1e5 * z,
    ξ_initial = [z -> 1000.0, z -> 0.0],
    tspan = (0.0, 0.8),
    output_times = collect(0.0:0.1:0.8),
)

# A new problem on the base problem's own mesh.
_reuse_mesh(problem, A) = PyrolysisProblem(
    mesh = problem.mesh,
    material = _ablator_material(A),
    bc_set = problem.bc_set,
    T_initial = problem.T_initial,
    ξ_initial = problem.ξ_initial,
    tspan = problem.tspan,
    dt_initial = problem.dt_initial,
    output_times = problem.output_times,
)

_fields(sol) = (sol.t, sol.T, sol.m, sol.z, [e.t for e in sol.events])

@testset "parameter sweep: points that share the base mesh leave it unchanged" begin
    # The factor of a front receding at 2 mm/s from 700 K, and twice it.
    R, E, Tp, v = 8.314462618, 1.0e5, 700.0, 2.0e-3
    A = v^2 * E * (Tp - 300.0) * exp(E / (R * Tp)) / (1e-6 * R * Tp^2)
    grid = [A, 2A, A]
    kw = (use_ale = true, deactivation = true, abstol = 1e-7, reltol = 1e-5)
    base = _ablator_problem(A)
    before = deepcopy(base.mesh)
    swept = parameter_sweep(
        base,
        grid;
        modify = _reuse_mesh,
        extract = _fields,
        parallel = false,
        kw...,
    )
    @test isequal(base.mesh.positions, before.positions)
    topology, topology0 = base.mesh.topology, before.topology
    for name in fieldnames(typeof(topology))
        @test isequal(getfield(topology, name), getfield(topology0, name))
    end
    # Each point is the solve of its own problem on a fresh mesh, and the
    # solves moved their nodes and deactivated cells.
    for (a, result) in zip(grid, swept)
        @test isequal(result, _fields(Integration.solve(_ablator_problem(a); kw...)))
        @test !isempty(result[5])
    end
    @test !isequal(swept[1], swept[2])
end

@testset "parameter sweep on the finite-volume core: two workers" begin
    workers = addprocs(2; exeflags = "--project=$(Base.active_project())")
    try
        Distributed.remotecall_eval(Main, workers, _SWEEP)
        base = Main._sweep_problem(1e3)
        options = (; modify = Main._sweep_modify, extract = Main._sweep_extract)
        @test nworkers() == 2
        parallel = parameter_sweep(base, _GRID; options...)
        @test parallel == parameter_sweep(base, _GRID; options..., parallel = false)
    finally
        rmprocs(workers)
    end
end
