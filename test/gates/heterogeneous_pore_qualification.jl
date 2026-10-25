# Heterogeneous ALE transport does not preserve weighted matrix occupancy
# exactly at finite resolution. Qualify actual pore states and refinement,
# rather than interpreting nominal porosity as the evolving pore volume.
using Test, Pyrolysis
include(joinpath(@__DIR__, "../verification/casedefs.jl"))
load_verification_cases()

function qualify_heterogeneous_pores(io)
    println(
        io,
        "energy,cells,porosity,pressure_tolerance,peak_generation,peak_temperature,minimum_porosity,maximum_matrix_occupancy,gas_ledger_relative,accepted,rejected",
    )
    @testset "heterogeneous nominal pore volume" for energy in (:temperature, :enthalpy)
        values=Dict{Tuple{Int,Float64,Float64},NTuple{3,Float64}}()
        for n in (15, 30, 60), phi in (0.1, 0.15), ptol in (0.001, 0.0001)
            rec=ConservationRecorder()
            sol=solve(
                _vbo07_problem(n; porosity = phi);
                VBO_SOLVE_KWARGS...,
                energy_state = energy,
                pressure_tolerance = ptol,
                callback = rec,
            )
            @test completed(sol)
            raw=get_porosity(sol)
            minphi=minimum(filter(isfinite, vec(raw)))
            qmax=maximum(
                (1-raw[i, k])/(1-phi) for k in eachindex(sol.t), i in axes(raw, 1) if
                isfinite(raw[i, k]) && sol.ξ[3, i, k]==0
            )
            a=gas_species_audit(sol, rec).species.MMA
            @test abs(a.relative)<1e-4
            @test minphi>0
            @test qmax < (n==15 ? 1.10 : n==30 ? 1.04 : 1.02)
            peak=_vbo_measures(sol).peak/1000
            Tmax=maximum(filter(isfinite, vec(sol.T)))
            values[(n, phi, ptol)]=(peak, Tmax, qmax)
            row=(
                energy,
                n,
                phi,
                ptol,
                peak,
                Tmax,
                minphi,
                qmax,
                a.relative,
                sol.statistics.accepted,
                sol.statistics.rejected,
            )
            println(io, join(row, ','))
            flush(io)
            println(row)
            flush(stdout)
        end
        for n in (15, 30, 60)
            q=[values[(n, phi, ptol)] for phi in (0.1, 0.15), ptol in (0.001, 0.0001)]
            @test maximum(first, q)-minimum(first, q)<0.001*maximum(first, q)
            @test maximum(x->x[2], q)-minimum(x->x[2], q)<0.1
        end
        for phi in (0.1, 0.15), ptol in (0.001, 0.0001)
            @test values[(60, phi, ptol)][3]-1 < 0.6*(values[(30, phi, ptol)][3]-1)
            @test values[(30, phi, ptol)][3]-1 < 0.6*(values[(15, phi, ptol)][3]-1)
        end
    end
end
path=isempty(ARGS) ? "/tmp/pyrolysis-heterogeneous-pore-qualification.csv" : ARGS[1]
open(path, "w") do io
    Base.invokelatest(qualify_heterogeneous_pores, io)
end
