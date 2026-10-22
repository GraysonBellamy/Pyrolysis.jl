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
