isdefined(@__MODULE__, :SolverProbe) || include("solver_probe.jl")
isdefined(@__MODULE__, :RestartProbe) || include("restart_probe.jl")
using .RestartProbe

function write_callback_csv(name, rows)
    path = joinpath(@__DIR__, "results", name)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(string.(values(row)), ','))
        end
    end
    metric("CSV", path)
end

@testset "callback restart controls" begin
    controls = Dict(
        mode => callback_control(mode) for
        mode in (:clean, :dirty_then_clean, :clean_then_dirty)
    )
    clean = controls[:clean]
    @test SciMLBase.successful_retcode(clean.sol)
    @test maximum(clean.orders) >= 3
    for mode in (:dirty_then_clean, :clean_then_dirty)
        dirty = controls[mode]
        @test dirty.sol.retcode == SciMLBase.ReturnCode.MaxIters
        @test last(dirty.sol.t) < 1.0
        @test maximum(dirty.orders) == 1
        @test dirty.sol.stats.naccept > 10clean.sol.stats.naccept
    end
    @test controls[:dirty_then_clean].sol.stats.naccept ==
          controls[:clean_then_dirty].sol.stats.naccept
    write_callback_csv(
        "callback_restart_controls.csv",
        [
            (;
                mode,
                retcode = string(result.sol.retcode),
                final_time = last(result.sol.t),
                iteration_budget = 10000,
                accepted = result.sol.stats.naccept,
                rejected = result.sol.stats.nreject,
                rhs_calls = result.sol.stats.nf,
                jacobians = result.sol.stats.njacs,
                min_order = minimum(result.orders),
                max_order = maximum(result.orders),
                max_error = maximum(
                    abs(u[1]-exp(-t)) for (u, t) in zip(result.sol.u, result.sol.t)
                ),
            ) for (mode, result) in sort(collect(controls); by = first)
        ],
    )
    # A correctly flagged observer is numerically inert, including when
    # real events are present. SciML saves an extra identical left-limit
    # snapshot when an observer fires before an event (three records at
    # that time instead of two). Remove only exactly identical (t,u)
    # pairs; retain both distinct one-sided states and require bitwise
    # agreement of the resulting trajectories and step/RHS/Jacobian counters.
    for events in (0, 10)
        clean = run_restarts(8, :FBDF, :lifted, :colored, events)
        absent = run_restarts(8, :FBDF, :lifted, :colored, events; observer = :none)
        clean_history = unique(collect(zip(clean.sol.t, clean.sol.u)))
        absent_history = unique(collect(zip(absent.sol.t, absent.sol.u)))
        @test first.(clean_history) == first.(absent_history)
        @test last.(clean_history) == last.(absent_history)
        @test clean.sol.stats.naccept == absent.sol.stats.naccept
        @test clean.sol.stats.nreject == absent.sol.stats.nreject
        @test clean.sol.stats.nf == absent.sol.stats.nf
        @test clean.sol.stats.njacs == absent.sol.stats.njacs
    end
end
