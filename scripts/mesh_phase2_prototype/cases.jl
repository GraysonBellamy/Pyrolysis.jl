using Pyrolysis, Statistics
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO, "test", "verification", "casedefs.jl"))
ids = split(ARGS[1], ',')
cases = filter(
    c -> any(id -> c.id == id || startswith(c.id, id), ids),
    load_verification_cases(),
)
tag = VERIFICATION_CORE == "new" ? "new/$(VERIFICATION_ENERGY_STATE)" : "old"
nfail = 0
_stats(sol) =
    hasproperty(sol, :statistics) ? sol.statistics :
    (hasproperty(sol, :raw_solution) ? sol.raw_solution.stats : "")
_events(sol) =
    hasproperty(sol, :events) ? " events=$(length(sol.events)) outcome=$(sol.outcome)" :
    (
        hasproperty(sol, :extras) && hasproperty(sol.extras, :park_times) ?
        " parks=$(length(sol.extras.park_times))" : ""
    )
const NREP = get(ENV, "CASES_NREP", "3")
for case in cases
    try
        _case_solve(case, case.build(case.ci.n_cells); case.solve_kwargs...)
    catch err
        println(
            "[$(case.id) $tag] warm-up error: ",
            sprint(showerror, err)[1:min(end, 1200)],
        )
        global nfail += 1
        continue
    end
    times = Float64[]
    sol = nothing
    for r = 1:parse(Int, NREP)
        t = @elapsed s =
            _case_solve(case, case.build(case.ci.n_cells); case.solve_kwargs...)
        push!(times, t)
        sol = s
    end
    t = median(times)
    ok = sol.retcode == case.expected_retcode
    println(
        "[$(case.id) $tag] retcode=$(sol.retcode) wall_median=$(round(t, digits=3)) s of $(round.(times, digits=3)) (budget $(case.ci.budget_s)) t_end=$(sol.t[end])",
        _events(sol),
        " ",
        _stats(sol),
    )
    for q in Base.invokelatest(case.qois, sol)
        pass = q.error < q.tol
        ok &= pass
        println(
            "    ",
            pass ? "PASS" : "FAIL",
            "  ",
            q.name,
            ": value=",
            q.value,
            " error=",
            q.error,
            " tol=",
            q.tol,
        )
    end
    ok &= t < 5 * case.ci.budget_s
    global nfail += !ok
end
println(nfail == 0 ? "ALL PASS ($tag)" : "$nfail case(s) FAILED ($tag)")
