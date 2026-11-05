# Pyrolysis.jl's side of the cross-code cost comparison: for every case with
# matched FDS/Gpyro/ThermaKin runs, the warm solve time and the analytical error
# on exactly the curves and references the matched runs are scored on
# (test/verification/reference/matched/*.toml, scored by add_matched_comparisons!
# in scripts/verification/comparisons.jl), over the case's report resolution
# ladder, the comparison meshes (reference/comparison_targets.toml) and its
# cost-only rungs (`report.cost_n_cells`), each at the report's solve settings
# and at tolerances TOLERANCE_FACTORS times looser: the counterpart of the
# external codes' coarser steps. A looser level that fails is left out.
#
# Run from the root of the checkout being measured, pinned, one thread
# (scripts/timing/campaign.sh step `matched` does this):
#   julia --project=scripts/verification --startup-file=no --threads=1 \
#       scripts/timing/pyrolysis_matched.jl OUTDIR [case ids|all]
#
# Per rung and tolerance: one untimed warm-up solve for the case's first rung
# (compilation), then repeated solves with a full collection before each, 5
# when a solve takes under 1 s, 3 under 10 s, else 1. Writes
# OUTDIR/matched_times.csv (one row per timed solve) and
# OUTDIR/matched_errors.csv (one row per rung, tolerance, curve and reference),
# flushed after every case.

using Statistics: median
using TOML
using Printf

const ROOT = pwd()
include(joinpath(ROOT, "test", "verification", "casedefs.jl"))
include(joinpath(ROOT, "scripts", "verification", "report_harness.jl"))
include(joinpath(ROOT, "scripts", "verification", "comparisons.jl"))

const OUT = abspath(ARGS[1])
const WANTED = length(ARGS) >= 2 && ARGS[2] != "all" ? split(ARGS[2], ',') : nothing
mkpath(OUT)

repeats_for(seconds) = seconds < 1 ? 5 : seconds < 10 ? 3 : 1

# Multipliers on the report's reltol, abstol and pressure_tolerance.
const TOLERANCE_FACTORS = (1, 10, 100, 1000)

function loosened(skw, factor)
    factor == 1 && return skw
    haskey(skw, :reltol) && haskey(skw, :abstol) ||
        error("report solve settings name no reltol/abstol to loosen")
    keys_ = filter(k -> haskey(skw, k), (:reltol, :abstol, :pressure_tolerance))
    return merge(skw, NamedTuple{keys_}(Tuple(factor * skw[k] for k in keys_)))
end

"""Solve at `factor` and record it, or report a looser level as unavailable."""
function attempt(solve_once, case, label, factor)
    try
        times, bytes, sol = timed_solves(solve_once)
        sol.retcode == case.expected_retcode || error("retcode $(sol.retcode)")
        return times, bytes, sol
    catch err
        factor == 1 && rethrow()
        println(
            "$(case.id) $label tolerance ×$factor: unavailable: ",
            first(sprint(showerror, err), 200),
        )
        return nothing
    end
end

rungs(case) = sort(
    unique(
        vcat(
            [t["n_cells"] for t in TARGETS if t["case"] == case.id],
            collect(Int, get(case.report, :n_cells_ladder, [case.ci.n_cells])),
            collect(Int, get(case.report, :cost_n_cells, Int[])),
        ),
    ),
)

"""Timed solves of `solve_once()` after its warm-up; returns times, bytes and
the last solution."""
function timed_solves(solve_once)
    times, bytes = Float64[], Int[]
    local sol
    GC.gc()
    r = @timed solve_once()
    push!(times, r.time)
    push!(bytes, r.bytes)
    sol = r.value
    for _ = 2:repeats_for(r.time)
        GC.gc()
        r = @timed solve_once()
        push!(times, r.time)
        push!(bytes, r.bytes)
        sol = r.value
    end
    return times, bytes, sol
end

stat(sol, f) =
    try
        getproperty(sol.statistics, f)
    catch
        -1
    end

function errors_of(curves, entries)
    rows = Any[]
    for e in entries
        ci = findfirst(c -> string(c.name) == e["curve"], curves)
        ci === nothing && error("missing curve $(e["curve"])")
        s = findfirst(s -> !s.crosscode && s.label == e["reference"], curves[ci].series)
        s === nothing && error("missing reference $(e["reference"])")
        series = curves[ci].series[s]
        err = abs.(series.y .- series.exact.(series.x))
        push!(
            rows,
            (
                e["curve"],
                e["reference"],
                length(err),
                maximum(err),
                sqrt(sum(abs2, err) / length(err)),
            ),
        )
    end
    return rows
end

const TARGETS =
    TOML.parsefile(joinpath(REFERENCE_DIR, "comparison_targets.toml"))["targets"]
const CASES = Dict(c.id => c for c in load_verification_cases())

time_rows = Any[]
error_rows = Any[]

function flush_tables()
    write_csv(
        joinpath(OUT, "matched_times.csv"),
        [
            "case",
            "variant",
            "n_cells",
            "tolerance",
            "repeat",
            "seconds",
            "bytes",
            "accepted",
            "rejected",
            "rhs_evaluations",
            "jacobians",
        ],
        time_rows,
    )
    write_csv(
        joinpath(OUT, "matched_errors.csv"),
        [
            "case",
            "variant",
            "n_cells",
            "tolerance",
            "median_seconds",
            "curve",
            "reference",
            "samples",
            "Linf",
            "RMS",
        ],
        error_rows,
    )
end

function record!(id, variant, n, factor, times, bytes, sol, errs)
    for (i, t) in enumerate(times)
        push!(
            time_rows,
            (
                id,
                variant,
                n,
                factor,
                i,
                t,
                bytes[i],
                stat(sol, :accepted),
                stat(sol, :rejected),
                stat(sol, :rhs_evaluations),
                stat(sol, :jacobians),
            ),
        )
    end
    for e in errs
        push!(error_rows, (id, variant, n, factor, median(times), e...))
    end
    @printf(
        "%-9s %-3s n=%4d tol×%-4d %10.5f s (×%d)  steps %d\n",
        id,
        variant,
        n,
        factor,
        median(times),
        length(times),
        stat(sol, :accepted)
    )
    flush(stdout)
end

# V-CPL-03's report curves re-solve every path themselves, so its matched
# paths are timed and scored here, one problem per path.
function vcpl03!(case, entries)
    skw = report_solve_kwargs(case)
    ns = rungs(case)
    for e in entries
        perm = Symbol(e["reference"])
        build = n -> _vcpl03_problem(perm; n_cells = n, n_outputs = 33)
        _case_solve(case, build(ns[1]).problem; skw...)  # warm-up
        for factor in TOLERANCE_FACTORS, n in ns
            r = build(n)
            kw = loosened(skw, factor)
            run = attempt(
                () -> _case_solve(case, r.problem; kw...),
                case,
                "$perm n=$n",
                factor,
            )
            run === nothing && continue
            times, bytes, sol = run
            # The case's exact path for P1–P3 (no lateral law): its report curve.
            exact = t -> 1.0 + (r.exact_Vf - 1.0) * (1.0 - exp(-r.A_r * t))
            y = sol.extras.thickness_history ./ r.L
            err = abs.(y .- exact.(sol.t))
            record!(
                case.id,
                string(perm),
                n,
                factor,
                times,
                bytes,
                sol,
                [(
                    e["curve"],
                    e["reference"],
                    length(err),
                    maximum(err),
                    sqrt(sum(abs2, err) / length(err)),
                )],
            )
        end
    end
end

function generic!(case, entries)
    skw = report_solve_kwargs(case)
    ns = rungs(case)
    _case_solve(case, case.build(ns[1]); skw...)  # warm-up: compilation
    for factor in TOLERANCE_FACTORS, n in ns
        problem = case.build(n)
        kw = loosened(skw, factor)
        run = attempt(() -> _case_solve(case, problem; kw...), case, "n=$n", factor)
        run === nothing && continue
        times, bytes, sol = run
        curves = case.report.curves(sol)
        record!(case.id, "", n, factor, times, bytes, sol, errors_of(curves, entries))
    end
end

matched_ids = sort(unique(e["case"] for e in MATCHED_COMPARISONS))
family(id) = replace(id, r"(?<=\d)[a-z]+$" => "")
WANTED === nothing || filter!(id -> id in WANTED || family(id) in WANTED, matched_ids)
failures = String[]
for id in matched_ids
    case = CASES[id]
    entries = unique(
        e -> (e["curve"], e["reference"]),
        filter(e -> e["case"] == id, MATCHED_COMPARISONS),
    )
    try
        Base.invokelatest(startswith(id, "V-CPL-03") ? vcpl03! : generic!, case, entries)
    catch err
        push!(failures, id)
        @error "case failed" id exception = (err, catch_backtrace())
    end
    flush_tables()
end
isempty(failures) || println("FAILED: ", join(failures, ", "))
println("done: ", length(matched_ids) - length(failures), " cases")
