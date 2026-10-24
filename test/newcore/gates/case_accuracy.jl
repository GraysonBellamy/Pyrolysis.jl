# Gate items 6 and 8 of the Phase 1 plan (§7.3) and item 7 of the Phase 2
# plan: on every case the suite runs on the finite-volume core
# (`NEW_CORE_CASES`), in both energy forms, the largest scaled constraint
# residual of the accepted states and the achieved temperature error against
# the requested tolerance, and on the deactivating cases the event timing with
# event location off and on.
#
# Each case is solved as the suite solves it (its builder at the CI
# resolution, its solve keywords and its solve function), then again with both
# tolerances divided by 100, the reference. The achieved error is the largest
# temperature difference over the saved times and the cells live in both
# solves; under deactivation only the saved times at which both solves have
# the same live cells count (between two solves' event times a cell is merged
# in one and not in the other), and a case saved at every step (no output
# times) compares its final state. When one solve stops before the other (a
# reference that does not complete), the saved times both reached are
# compared, and the last of them is reported. It is reported as a multiple of
# the requested relative tolerance times the reference's largest temperature,
# and is within its budget at ten times that or less (the budget of the
# leftover-2 table). A reference that does not complete at a hundredth of the
# tolerances is taken at a tenth. The constraint residual is the main solve's
# `statistics.largest_constraint_residual`; one or below passes.
#
# Under deactivation (`handle_depletion`) a temperature just after an event
# follows the event's time, so each such case is solved with event location
# off and on (`located_events`, its reference likewise) and judged on the
# event-time shift instead: the largest difference of the deactivation
# invocations' times between the solve and its reference, over the median
# interval between the reference's invocations, or, for a case with a single
# invocation, between its saved times (where its quantities are read). A case
# passes when the located shift is at most `LOCATED_BUDGET` of that interval
# and at most a `LOCATED_GAIN`-th of the shift without location.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/case_accuracy.jl [case ids...]
# With case ids, only those cases. Writes
# test/newcore/gates/measurements/case_accuracy.{csv,md}.

ENV["PYROLYSIS_CORE"] = "new"
using Pyrolysis
include(joinpath(@__DIR__, "event_shift.jl"))

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const TIGHTEN = (100.0, 10.0)
const BUDGET = 10.0
const LOCATED_BUDGET = 0.05
const LOCATED_GAIN = 3.0
const ONLY = ARGS
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

# The pairs of saved-time indices `(j, k)` with `sol.t[j] == ref.t[k]` and the
# same live cells in both solves.
function _common_times(sol, ref)
    live_s, live_r = _live(sol), _live(ref)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    pairs = Tuple{Int,Int}[]
    for (j, t) in enumerate(sol.t)
        k = get(index, t, 0)
        k > 0 && live_s[j] == live_r[k] && push!(pairs, (j, k))
    end
    return pairs
end

# The largest temperature difference over the compared saved times and the
# cells live in both solves, the largest reference temperature, the number of
# saved times compared (zero for the final state) and the last time compared.
function temperature_error(sol, ref)
    size(sol.T, 1) == size(ref.T, 1) || return NaN, NaN, 0, NaN
    if size(sol.T) == size(ref.T)
        pairs = [(k, k) for k in findall(sol.t .== ref.t .&& _live(sol) .== _live(ref))]
    elseif sol.t[end] == ref.t[end]
        pairs = [(length(sol.t), length(ref.t))]
    else
        pairs = filter(p -> sol.t[p[1]] > first(sol.t), _common_times(sol, ref))
    end
    isempty(pairs) && return NaN, NaN, 0, NaN
    err, scale = 0.0, 0.0
    for (j, k) in pairs
        a, b = view(sol.T, :, j), view(ref.T, :, k)
        for i in eachindex(a, b)
            (isfinite(a[i]) && isfinite(b[i])) || continue
            err = max(err, abs(a[i] - b[i]))
            scale = max(scale, abs(b[i]))
        end
    end
    final = size(sol.T) != size(ref.T) && sol.t[end] == ref.t[end]
    return err, scale, final ? 0 : length(pairs), sol.t[last(pairs)[1]]
end

_statistic(sol, name) =
    hasproperty(sol, :statistics) && sol.statistics !== nothing ?
    getproperty(sol.statistics, name) : NaN

# `located` is `nothing` for a case without deactivation, else whether its
# events are located.
function survey(m, case, energy; located = nothing)
    kw = case.solve_kwargs
    located === nothing || (kw = merge(kw, (; located_events = located)))
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
    err, scale, compared, compared_to = temperature_error(sol, ref)
    ratio = err / (reltol * scale)
    shift, interval = event_shift(sol, ref)
    constraint = _statistic(sol, :largest_constraint_residual)
    return (;
        case = case.id,
        energy,
        location = located === nothing ? "" : located ? "on" : "off",
        retcode = sol.retcode,
        reference_retcode = ref.retcode,
        tightened,
        reltol,
        abstol,
        error = err,
        scale,
        compared,
        compared_to,
        ratio,
        invocations = hasproperty(sol, :events) ? length(invocation_times(sol)) : 0,
        shift,
        interval,
        event_shift = shift / interval,
        completed = sol.retcode == case.expected_retcode,
        within_budget = sol.retcode == case.expected_retcode && ratio <= BUDGET,
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
        isempty(ONLY) || case.id in ONLY || continue
        located =
            get(case.solve_kwargs, :handle_depletion, false) ? (false, true) : (nothing,)
        for l in located
            row = try
                survey(m, case, energy; located = l)
            catch err
                println(
                    "[$(case.id) $energy] error: ",
                    sprint(showerror, err)[1:min(end, 400)],
                )
                continue
            end
            push!(rows, row)
            println(join(("$k=$v" for (k, v) in pairs(row)), " "))
            flush(stdout)
        end
    end
end
load_end = Sys.loadavg()

# Item 7, per deactivating case and form: the shifts with location off and
# on, as fractions of the median interval.
timing = []
for r in rows
    r.location == "on" || continue
    off = only(
        filter(o -> o.case == r.case && o.energy == r.energy && o.location == "off", rows),
    )
    gain = off.shift / r.shift
    push!(
        timing,
        (;
            r.case,
            r.energy,
            r.invocations,
            off_shift = off.shift,
            on_shift = r.shift,
            interval = r.interval,
            off_fraction = off.shift / r.interval,
            on_fraction = r.shift / r.interval,
            gain,
            passes = off.completed &&
                     r.completed &&
                     r.shift / r.interval <= LOCATED_BUDGET &&
                     gain >= LOCATED_GAIN,
        ),
    )
end
# A deactivating case's solves are judged on item 7.
rows = [
    r.location == "" ? r :
    merge(
        r,
        (;
            within_budget = any(
                t -> t.case == r.case && t.energy == r.energy && t.passes,
                timing,
            ),
        ),
    ) for r in rows
]

open(joinpath(OUT, "case_accuracy.csv"), "w") do io
    csv_row(io, string.(keys(first(rows)))...)
    for row in rows
        csv_row(io, values(row)...)
    end
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "case_accuracy.md"), "w") do io
    println(
        io,
        "# Gate items 6 and 8 (Phase 1) and 7 (Phase 2): constraint residuals, achieved error and event timing\n",
    )
    println(
        io,
        "Produced by `test/newcore/gates/case_accuracy.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. ",
        isempty(ONLY) ? "Every case of `NEW_CORE_CASES`" : "The cases $(join(ONLY, ", "))",
        ", each at its CI resolution and solve ",
        "keywords, on the default integrator. Error: the largest temperature difference ",
        "[K] over the saved times and the cells live in both solves, against the same ",
        "case at both tolerances divided by $(Int(TIGHTEN[1])), or by $(Int(TIGHTEN[2])) ",
        "where that does not complete (tightened) (under deactivation the saved ",
        "times with the same live cells in both; a case saved at every step its final ",
        "state; when one solve stops before the other, the saved times both reached; ",
        "compared: the saved times compared, zero for the final state; to: the last time ",
        "compared [s]); ratio: that error over ",
        "reltol times the reference's largest temperature, within budget at ",
        "$(Int(BUDGET)) or less. Under deactivation each case is solved with event ",
        "location off and on (location; its reference likewise) and judged on item 7 ",
        "below: invocations, the deactivation invocations; event shift, the largest ",
        "difference of the invocation times from the reference's over the median ",
        "interval between the reference's invocations (between its saved times for a ",
        "single invocation). Constraint: the largest scaled constraint residual of ",
        "the accepted states (one or below passes; zero for a fixed-step solve or a ",
        "configuration without algebraic rows). Wall: one warm solve [s].\n",
    )
    println(
        io,
        "| case | form | location | retcode | reltol | tightened | error [K] | compared | to [s] | ratio | invocations | event shift | within budget | constraint | steps | wall [s] |",
    )
    println(io, "|" * repeat("---|", 16))
    for r in rows
        println(
            io,
            "| $(r.case) | $(r.energy) | $(r.location) | $(r.retcode) | $(r.reltol) | $(Int(r.tightened)) | ",
            "$(_g(r.error)) | $(r.compared) | $(_g(r.compared_to)) | $(_g(r.ratio)) | $(r.invocations) | $(_g(r.event_shift)) | $(r.within_budget) | $(_g(r.constraint)) | ",
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
    isempty(timing) && return
    println(io, "\n## Item 7: event timing with location off and on\n")
    println(
        io,
        "Per deactivating case and form, the largest shift of the invocation times ",
        "between the solve and its reference [s] with location off and on, the median ",
        "interval between the reference's invocations, or between its saved times for ",
        "a single invocation [s], and the shifts as fractions of it, and gain, the shift off over the shift on. Passes when the located ",
        "shift is at most $(LOCATED_BUDGET) of the interval and the gain at least ",
        "$(LOCATED_GAIN), both solves completing as the case expects.\n",
    )
    println(
        io,
        "| case | form | invocations | shift off [s] | shift on [s] | interval [s] | off / interval | on / interval | gain | passes |",
    )
    println(io, "|" * repeat("---|", 10))
    for t in timing
        println(
            io,
            "| $(t.case) | $(t.energy) | $(t.invocations) | $(_g(t.off_shift)) | $(_g(t.on_shift)) | ",
            "$(_g(t.interval)) | $(_g(t.off_fraction)) | $(_g(t.on_fraction)) | $(_g(t.gain)) | $(t.passes) |",
        )
    end
    println(io, "\n$(count(t -> t.passes, timing)) of $(length(timing)) pass.")
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
