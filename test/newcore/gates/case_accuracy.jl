# Gate items 6 and 8 of the Phase 1 plan (§7.3): on every case the suite runs
# on the finite-volume core (`NEW_CORE_CASES`), in both energy forms, the
# largest scaled constraint residual of the accepted states and the achieved
# temperature error against the requested tolerance.
#
# Each case is solved as the suite solves it (its builder at the CI
# resolution, its solve keywords and its solve function), then again with
# both tolerances divided by 100, the reference. The achieved error is the
# largest temperature difference over the saved times and the cells live in
# both solves; under deactivation only the saved times at which both solves
# have the same live cells count (between two solves' event times a cell is
# merged in one and not in the other), and a case saved at every step (no
# output times) compares its final state. It is reported as a multiple of the
# requested relative tolerance times the reference's largest temperature, and
# is within its budget at ten times that or less (the budget of the leftover-2
# table). A reference that does not complete at a hundredth of the tolerances
# is taken at a tenth. Under deactivation a temperature just after an event
# follows the event's time, so those cases report the event-time shift too:
# the largest difference of the event times over the median interval between
# the reference's events, as a multiple of reltol, within budget at ten or
# less; their budget is judged on it. The constraint residual is the main
# solve's `statistics.largest_constraint_residual`; one or below passes.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/case_accuracy.jl
# Writes test/newcore/gates/measurements/case_accuracy.{csv,md}.

ENV["PYROLYSIS_CORE"] = "new"
using Pyrolysis
using Statistics: median

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const TIGHTEN = (100.0, 10.0)
const BUDGET = 10.0
# The tolerances of `Integration.solve`, which a case without a solve function
# and without its own tolerances is solved at.
const DEFAULT_TOL = (reltol = 1e-6, abstol = 1e-8)

function csv_row(io, vals...)
    println(io, join((v isa AbstractString ? "\"$v\"" : string(v) for v in vals), ","))
end

# The case runner reads the energy state from the environment when casedefs.jl
# is included, so each form is a separate module.
function load_cases(energy::String)
    ENV["PYROLYSIS_ENERGY_STATE"] = energy
    m = Module(Symbol("cases_", energy))
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
    cases = Base.invokelatest(m.load_verification_cases)
    return m, filter(c -> Base.invokelatest(m.runs_on_core, c), cases)
end

_live(sol) =
    hasproperty(sol, :extras) && hasproperty(sol.extras, :live_cells) ?
    sol.extras.live_cells : fill(0, length(sol.t))

# The largest temperature difference over the compared saved times and the
# cells live in both solves, the largest reference temperature, and the
# number of saved times compared.
function temperature_error(sol, ref)
    if size(sol.T) == size(ref.T)
        times = findall(sol.t .== ref.t .&& _live(sol) .== _live(ref))
    elseif sol.t[end] == ref.t[end] && size(sol.T, 1) == size(ref.T, 1)
        times = [0]
    else
        return NaN, NaN, 0
    end
    err, scale = 0.0, 0.0
    for k in times
        a = k == 0 ? view(sol.T, :, size(sol.T, 2)) : view(sol.T, :, k)
        b = k == 0 ? view(ref.T, :, size(ref.T, 2)) : view(ref.T, :, k)
        for i in eachindex(a, b)
            (isfinite(a[i]) && isfinite(b[i])) || continue
            err = max(err, abs(a[i] - b[i]))
            scale = max(scale, abs(b[i]))
        end
    end
    return err, scale, times == [0] ? 0 : length(times)
end

# The largest shift of the event times over the median interval between the
# reference's events (NaN without events or when the counts differ).
function event_shift(sol, ref)
    hasproperty(sol, :events) || return NaN
    a, b = [e.t for e in sol.events], [e.t for e in ref.events]
    (isempty(b) || length(a) != length(b)) && return NaN
    gaps = filter(>(0), diff(unique(b)))
    isempty(gaps) && return NaN
    return maximum(abs, a .- b) / median(gaps)
end

_statistic(sol, name) =
    hasproperty(sol, :statistics) && sol.statistics !== nothing ?
    getproperty(sol.statistics, name) : NaN

function survey(m, case, energy)
    kw = case.solve_kwargs
    if case.solve_fn !== nothing && !(haskey(kw, :reltol) && haskey(kw, :abstol))
        error("$(case.id): a case with a solve function needs explicit tolerances")
    end
    reltol = get(kw, :reltol, DEFAULT_TOL.reltol)
    abstol = get(kw, :abstol, DEFAULT_TOL.abstol)
    solve_case(; tol...) = Base.invokelatest(
        m._case_solve,
        case,
        Base.invokelatest(case.build, case.ci.n_cells);
        merge(kw, (; tol...))...,
    )
    sol = solve_case()
    wall = @elapsed sol = solve_case()
    ref, tightened = nothing, NaN
    for factor in TIGHTEN
        ref = solve_case(reltol = reltol / factor, abstol = abstol / factor)
        tightened = factor
        ref.retcode == case.expected_retcode && break
    end
    err, scale, compared = temperature_error(sol, ref)
    ratio = err / (reltol * scale)
    shift = event_shift(sol, ref) / reltol
    judged = isnan(shift) ? ratio : shift
    constraint = _statistic(sol, :largest_constraint_residual)
    return (;
        case = case.id,
        energy,
        retcode = sol.retcode,
        reference_retcode = ref.retcode,
        tightened,
        reltol,
        abstol,
        error = err,
        scale,
        compared,
        ratio,
        event_shift = shift,
        within_budget = sol.retcode == case.expected_retcode && judged <= BUDGET,
        constraint,
        constraint_passes = !(constraint > 1),
        steps = _statistic(sol, :accepted),
        wall,
    )
end

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
rows = []
for energy in ("temperature", "enthalpy")
    m, cases = load_cases(energy)
    for case in cases
        row = try
            survey(m, case, energy)
        catch err
            println("[$(case.id) $energy] error: ", sprint(showerror, err)[1:min(end, 400)])
            continue
        end
        push!(rows, row)
        println(join(("$k=$v" for (k, v) in pairs(row)), " "))
        flush(stdout)
    end
end
load_end = Sys.loadavg()

open(joinpath(OUT, "case_accuracy.csv"), "w") do io
    csv_row(io, string.(keys(first(rows)))...)
    for row in rows
        csv_row(io, values(row)...)
    end
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "case_accuracy.md"), "w") do io
    println(io, "# Gate items 6 and 8: constraint residuals and achieved error\n")
    println(
        io,
        "Produced by `test/newcore/gates/case_accuracy.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Every case of `NEW_CORE_CASES` at its CI resolution and solve ",
        "keywords, on the default integrator. Error: the largest temperature difference ",
        "[K] over the saved times and the cells live in both solves, against the same ",
        "case at both tolerances divided by $(Int(TIGHTEN[1])), or by $(Int(TIGHTEN[2])) ",
        "where that does not complete (tightened) (under deactivation the saved ",
        "times with the same live cells in both; a case saved at every step its final ",
        "state; compared: the saved times compared, zero for the final state); ratio: that error over ",
        "reltol times the reference's largest temperature, within budget at ",
        "$(Int(BUDGET)) or less; event shift: under deactivation the largest event-time ",
        "difference over the reference's median event interval, over reltol, on which ",
        "those cases' budget is judged. Constraint: the largest scaled constraint residual of ",
        "the accepted states (one or below passes; zero for a fixed-step solve or a ",
        "configuration without algebraic rows). Wall: one warm solve [s].\n",
    )
    println(
        io,
        "| case | form | retcode | reltol | tightened | error [K] | compared | ratio | event shift | within budget | constraint | steps | wall [s] |",
    )
    println(io, "|" * repeat("---|", 13))
    for r in rows
        println(
            io,
            "| $(r.case) | $(r.energy) | $(r.retcode) | $(r.reltol) | $(Int(r.tightened)) | ",
            "$(_g(r.error)) | $(r.compared) | $(_g(r.ratio)) | $(_g(r.event_shift)) | $(r.within_budget) | $(_g(r.constraint)) | ",
            "$(r.steps isa AbstractFloat && isnan(r.steps) ? "" : r.steps) | $(_g(r.wall)) |",
        )
    end
    n_budget = count(r -> r.within_budget, rows)
    n_constraint = count(r -> r.constraint_passes, rows)
    println(
        io,
        "\n$(n_budget) of $(length(rows)) solves within budget; $(n_constraint) of ",
        "$(length(rows)) with every constraint residual at or below one.",
    )
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
