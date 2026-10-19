# Gate item 2 of the Phase 1 plan (§7.3): the audit decks within 1.5× of the
# old core's wall time at matched accuracy. Its times also complete the
# affordability decision of leftover 1 (§7.1).
#
# The decks are the `problem_for` decks of the Phase 0 harness (`eulerian`,
# `ale`, `ale_darcy`, 80 s, three outputs) at 25 and 100 cells. The two PMMA
# decks give the solid a gas-transfer coefficient of 1e-6 m²/s: with none,
# gas transport vanishes in the solid, the gas forms a sharp front, and no
# tolerance converges either core's solution (at 100 cells the new core at
# reltol 1e-10 and 1e-11 differs by 0.02 K, the old core at 1e-8 and 1e-10 by
# up to 0.16 K, both above their errors at the default tolerances), so an
# error against a reference measures nothing there. The original PMMA decks
# are timed at equal tolerance beside them. Each core runs
# its default integrator on a tolerance ladder (reltol from 1e-4 to 1e-9,
# abstol a hundredth of it): the old core `solve(problem; kw...)`, the new core
# `Integration.solve` in both energy forms. The achieved error is the largest
# temperature difference over the saved times and the cells against the same
# core's own reference (the same mesh at `REFERENCE`), so it is the time error
# alone and the two cores' spatial discretizations do not enter. The
# reference spread is the difference between that reference and one at a
# tenth of its tolerances; an error near it is not resolved.
#
# Matched accuracy: the target is the old core's error at its default
# tolerances (abstol 1e-8, reltol 1e-6). The new core's matched time is the
# warm time of its cheapest ladder solve whose error is at or below the target
# (an upper bound: the ladder's spacing is half a decade); the ratio is that
# over the old core's time at its default. The ratio at equal tolerance (both
# at the old default) is reported beside it. Timed are the solves these
# ratios read: every core at its default, and the new core's first two
# ladder solves at or below the target.
#
# Warm medians of `REPEATS` solves after one warm-up, one thread.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/audit_decks.jl
# Writes test/newcore/gates/measurements/audit_decks.{csv,md}.

ENV["PYROLYSIS_CORE"] = "old"
const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "jacobian_common.jl"))


const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REPEATS = 3
# (kind, cells, the solid's gas-transfer coefficient: `nothing` keeps the deck's)
const DECKS = (
    (:eulerian, 25, 1e-6),
    (:eulerian, 100, 1e-6),
    (:ale, 25, 1e-6),
    (:ale, 100, 1e-6),
    (:ale_darcy, 25, nothing),
    (:ale_darcy, 100, nothing),
)
const DEGENERATE_DECKS = ((:eulerian, 25), (:eulerian, 100), (:ale, 25), (:ale, 100))
const LADDER = (1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8)
const DEFAULT_RELTOL = 1e-6
const REFERENCE = 1e-11
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

# The median of `REPEATS` timed solves, after a warm-up the caller ran.
function timed(f)
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed f())
    end
    return median(times)
end

_steps(sol) =
    hasproperty(sol, :statistics) ? sol.statistics.accepted : sol.raw_solution.stats.naccept

time_error(sol, ref) = maximum(abs, sol.T .- ref.T)

# The ladder of one core on one deck; `target` (the old core's error at its
# default) selects the solves of the new core that are timed.
function core_rows(deck, run, solve_at; target = nothing)
    ref = solve_at(_tol(REFERENCE))
    finer = solve_at(_tol(REFERENCE / 10))
    spread = time_error(ref, finer)
    rows = []
    hits = 0
    for reltol in LADDER
        sol = solve_at(_tol(reltol))
        err = sol.retcode === :Success ? time_error(sol, ref) : NaN
        hit = target !== nothing && err <= target && hits < 2
        hits += hit
        t = reltol == DEFAULT_RELTOL || hit ? timed(() -> solve_at(_tol(reltol))) : NaN
        push!(
            rows,
            (;
                deck,
                run,
                reltol,
                retcode = sol.retcode,
                error = err,
                spread,
                steps = _steps(sol),
                time = t,
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(rows[end])), " "))
        flush(stdout)
    end
    return rows
end

_deck(kind, n, solid_λ) =
    solid_λ === nothing ? "$(kind) $(n)" : "$(kind) $(n), solid λ $(solid_λ)"

function deck_rows(kind, n, solid_λ)
    problem, kw = problem_for(kind, n; solid_λ)
    deck = _deck(kind, n, solid_λ)
    old_at(tol) = solve(problem; kw..., tol...)
    new_at(energy_state) =
        tol -> Pyrolysis.Integration.solve(
            problem;
            energy_state,
            use_ale = kw.use_ale,
            kw.radiation_model,
            tol...,
        )
    old = core_rows(deck, "old core", old_at)
    target = only(r.error for r in old if r.reltol == DEFAULT_RELTOL)
    return vcat(
        old,
        core_rows(deck, "new core, temperature form", new_at(:temperature); target),
        core_rows(deck, "new core, enthalpy form", new_at(:enthalpy); target),
    )
end

# The original PMMA decks at the default tolerances: warm times only.
function degenerate_rows(kind, n)
    problem, kw = problem_for(kind, n)
    runs = (
        ("old core", () -> solve(problem; kw...)),
        (
            "new core, temperature form",
            () -> Pyrolysis.Integration.solve(
                problem;
                use_ale = kw.use_ale,
                kw.radiation_model,
            ),
        ),
        (
            "new core, enthalpy form",
            () -> Pyrolysis.Integration.solve(
                problem;
                energy_state = :enthalpy,
                use_ale = kw.use_ale,
                kw.radiation_model,
            ),
        ),
    )
    rows = []
    for (run, f) in runs
        sol = f()
        push!(
            rows,
            (;
                deck = "$(kind) $(n)",
                run,
                retcode = sol.retcode,
                steps = _steps(sol),
                time = timed(f),
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(rows[end])), " "))
        flush(stdout)
    end
    return rows
end

# The cheapest ladder solve of `run` on `deck` whose error is at or below
# `target`.
function matched(rows, deck, run, target)
    hits = filter(
        r -> r.deck == deck && r.run == run && r.error <= target && isfinite(r.time),
        rows,
    )
    isempty(hits) && return nothing
    return argmin(r -> r.time, hits)
end

at_default(rows, deck, run) =
    only(filter(r -> r.deck == deck && r.run == run && r.reltol == DEFAULT_RELTOL, rows))

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
rows = []
for (kind, n, solid_λ) in DECKS
    append!(rows, deck_rows(kind, n, solid_λ))
end
degenerate = []
for (kind, n) in DEGENERATE_DECKS
    append!(degenerate, degenerate_rows(kind, n))
end
load_end = Sys.loadavg()

summary = []
for (kind, n, solid_λ) in DECKS
    deck = _deck(kind, n, solid_λ)
    old = at_default(rows, deck, "old core")
    for run in ("new core, temperature form", "new core, enthalpy form")
        new = at_default(rows, deck, run)
        m = matched(rows, deck, run, old.error)
        push!(
            summary,
            (;
                deck,
                run,
                target = old.error,
                spread = old.spread,
                old_time = old.time,
                equal_tolerance_ratio = new.time / old.time,
                new_error_at_default = new.error,
                matched_reltol = m === nothing ? NaN : m.reltol,
                matched_error = m === nothing ? NaN : m.error,
                matched_time = m === nothing ? NaN : m.time,
                matched_ratio = m === nothing ? NaN : m.time / old.time,
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(summary[end])), " "))
    end
end

function csv_out(path, rs)
    open(path, "w") do io
        csv_row(io, string.(keys(first(rs)))...)
        for r in rs
            csv_row(io, values(r)...)
        end
    end
end
csv_out(joinpath(OUT, "audit_decks.csv"), rows)
csv_out(joinpath(OUT, "audit_decks_summary.csv"), summary)
csv_out(joinpath(OUT, "audit_decks_degenerate.csv"), degenerate)

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "audit_decks.md"), "w") do io
    println(io, "# Gate item 2: the audit decks at matched accuracy\n")
    println(
        io,
        "Produced by `test/newcore/gates/audit_decks.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Each core's default integrator; error: the largest temperature ",
        "difference [K] over the saved times and cells against the same core at reltol ",
        "$(REFERENCE) (abstol a hundredth of reltol); spread: that reference against one ",
        "at a tenth of its tolerances. Target: the old core's error at its default ",
        "tolerances (reltol $(DEFAULT_RELTOL)). Matched: the new core's cheapest ladder ",
        "solve at or below the target. Warm medians of $(REPEATS) solves [s].\n",
    )
    println(
        io,
        "| deck | new core | target [K] | spread [K] | old [s] | equal-tolerance ratio | new error at default [K] | matched reltol | matched error [K] | matched [s] | matched ratio |",
    )
    println(io, "|" * repeat("---|", 11))
    for s in summary
        println(
            io,
            "| $(s.deck) | $(replace(s.run, "new core, " => "")) | $(_g(s.target)) | ",
            "$(_g(s.spread)) | $(_g(s.old_time)) | $(_g(s.equal_tolerance_ratio)) | ",
            "$(_g(s.new_error_at_default)) | $(_g(s.matched_reltol)) | ",
            "$(_g(s.matched_error)) | $(_g(s.matched_time)) | $(_g(s.matched_ratio)) |",
        )
    end
    println(
        io,
        "\n## The original PMMA decks (solid λ = 0) at the default tolerances\n\n",
        "No tolerance converges these decks (see the script's header); warm times only, ",
        "and the ratio to the old core's.\n",
    )
    println(io, "| deck | run | retcode | steps | warm [s] | ratio |")
    println(io, "|---|---|---|---|---|---|")
    for r in degenerate
        old = only(d for d in degenerate if d.deck == r.deck && d.run == "old core")
        println(
            io,
            "| $(r.deck) | $(r.run) | $(r.retcode) | $(r.steps) | $(_g(r.time)) | ",
            "$(_g(r.time / old.time)) |",
        )
    end
    println(io, "\n## Ladders\n")
    println(io, "| deck | run | reltol | retcode | error [K] | steps | warm [s] |")
    println(io, "|---|---|---|---|---|---|---|")
    for r in rows
        println(
            io,
            "| $(r.deck) | $(r.run) | $(r.reltol) | $(r.retcode) | $(_g(r.error)) | ",
            "$(r.steps) | $(_g(r.time)) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
