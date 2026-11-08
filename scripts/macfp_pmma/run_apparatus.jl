# Pyrolysis.jl on the MaCFP-3 NIST gasification apparatus tests: every deck
# with complete properties × every test, on the shared apparatus deck.
#
#   julia --project=scripts/macfp_pmma scripts/macfp_pmma/run_apparatus.jl [SET ...] [--inert]
#
# Writes results/pyrolysis/<set>/<variant>/apparatus_<test>.csv and
# results/pyrolysis/<set>/apparatus_runs.csv; with no SET, or with --inert, also
# the inert tests, results/pyrolysis/inert/<test>.csv.
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, Pyrolysis, Printf

const OUT = joinpath(@__DIR__, "results", "pyrolysis")

# The cone heater ("base"); for an in-depth deck also the diffuse and collimated
# heaters, which bound it.
function heaters(deck)
    h = Pair{String,Any}["base"=>APPARATUS_HEATER]
    deck["absorption"] == "in-depth" &&
        append!(h, ["diffuse" => Diffuse(), "collimated" => Collimated()])
    return h
end

function run_set(deck)
    id = deck["id"]
    mkpath(joinpath(OUT, id))
    io = open(joinpath(OUT, id, "apparatus_runs.csv"), "w")
    println(
        io,
        "set,variant,test,retcode,outcome,t_last,solve_time,accepted,rejected,rhs_evaluations",
    )
    for (vname, distribution) in heaters(deck), t in APPARATUS_TESTS
        try
            p, iif = apparatus_problem(deck, t; distribution)
            sol =
                solve(p; SOLVE_KWARGS..., radiation_model = MacfpPMMA.radiation_model(deck))
            write_series(
                joinpath(OUT, id, vname, "apparatus_$(t.test).csv"),
                apparatus_series(sol, iif),
            )
            s = sol.statistics
            println(
                io,
                @sprintf(
                    "%s,%s,%s,%s,%s,%.2f,%.3f,%d,%d,%d",
                    id,
                    vname,
                    t.test,
                    sol.retcode,
                    sol.outcome,
                    sol.t[end],
                    sol.solve_time,
                    s.accepted,
                    s.rejected,
                    s.rhs_evaluations
                )
            )
            flush(io)
            @printf(
                "%-14s %-10s %-12s %-10s t=%6.1f %6.2f s\n",
                id,
                vname,
                t.test,
                sol.retcode,
                sol.t[end],
                sol.solve_time
            )
            flush(stdout)
        catch err
            println(
                io,
                "$id,$vname,$(t.test),ERROR,$(replace(sprint(showerror, err), r"[\n,]" => " ")[1:min(end, 300)])",
            )
            flush(io)
            @warn "failed" id vname t.test exception = (err, catch_backtrace())
        end
    end
    close(io)
end

function run_inert()
    for t in INERT_TESTS
        sol = solve(inert_problem(t); radiation_model = SURFACE_ABSORPTION)
        write_series(joinpath(OUT, "inert", "$(t.test).csv"), inert_series(sol, t))
        @printf(
            "%-16s %-10s t=%6.1f %6.2f s\n",
            t.test,
            sol.retcode,
            sol.t[end],
            sol.solve_time
        )
        flush(stdout)
    end
end

function main(args)
    check_matl_db()
    sets = filter(a -> !startswith(a, "--"), args)
    decks = load_decks(; gasification_only = true)
    isempty(sets) || filter!(d -> d["id"] in sets, decks)
    foreach(run_set, decks)
    (isempty(sets) || "--inert" in args) && run_inert()
end

main(ARGS)
