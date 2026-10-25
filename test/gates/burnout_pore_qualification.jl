# Joint nominal PMMA pore storage, mesh and pressure-error qualification.
# The calorically passive gas makes thermal/kinetic independence an analytic
# invariant of this deck. Gas escape and deactivation venting are audited
# separately; a thermal regression pin does not qualify their transport law.
using Test, Pyrolysis
include(joinpath(@__DIR__, "../verification/casedefs.jl"))
load_verification_cases()

function qualify_burnout_pores(io)
    println(
        io,
        "energy,cells,porosity,pressure_tolerance,retcode,outcome,peak_generation,peak_temperature,final_thickness,gas_escape,gas_vented,gas_ledger_relative,accepted,rejected,minimum_raw_porosity",
    )
    @testset "burnout nominal pore volume" for energy in (:temperature, :enthalpy)
        values=Dict{Tuple{Int,Float64,Float64},Tuple{Float64,Float64}}()
        for n in (15, 30, 60), φ in (0.03, 0.01, 0.003), ptol in (0.001, 0.0001)
            p=_vbo_problem(n; porosity = φ)
            rec=ConservationRecorder()
            sol=solve(
                p;
                VBO_SOLVE_KWARGS...,
                energy_state = energy,
                pressure_tolerance = ptol,
                callback = rec,
            )
            @test completed(sol)
            a=gas_species_audit(sol, rec).species.MMA
            @test abs(a.relative)<1e-4
            peak=_vbo_measures(sol).peak/1000
            Tmax=maximum(filter(isfinite, vec(sol.T)))
            values[(n, φ, ptol)]=(peak, Tmax)
            row=(
                energy,
                n,
                φ,
                ptol,
                sol.retcode,
                sol.outcome,
                peak,
                Tmax,
                sol.extras.thickness_history[end],
                a.escape,
                a.vented,
                a.relative,
                sol.statistics.accepted,
                sol.statistics.rejected,
                minimum(filter(isfinite, vec(get_porosity(sol)))),
            )
            println(io, join(row, ','))
            flush(io)
            println(row)
            flush(stdout)
        end
        for n in (15, 30, 60)
            q=[values[(n, φ, ptol)] for φ in (0.03, 0.01, 0.003), ptol in (0.001, 0.0001)]
            @test maximum(first, q)-minimum(first, q)<0.001*maximum(first, q)
            @test maximum(last, q)-minimum(last, q)<0.1 # K
        end
    end
end
path=isempty(ARGS) ? "/tmp/pyrolysis-burnout-pore-qualification.csv" : ARGS[1]
open(path, "w") do io
    Base.invokelatest(qualify_burnout_pores, io)
end
