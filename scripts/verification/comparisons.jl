# Matched, native-input comparisons. A manifest identifies the exact case,
# curve, and analytical reference; it cannot silently borrow another variant.
using TOML

const MATCHED_LEVELS = (
    "base",
    "time_double",
    "time_quadruple",
    "time_half",
    "time_quarter",
    "space_half",
    "refined",
)
# Coarser steps than the matched one: a code may be unable to run them.
const MATCHED_COARSE = ("time_double", "time_quadruple")
# Each level's mesh and time step relative to the external run's matched
# (base) mesh h and step Δt.
const LEVEL_LABELS = Dict(
    "base" => "h, Δt",
    "time_double" => "h, 2Δt",
    "time_quadruple" => "h, 4Δt",
    "time_half" => "h, Δt/2",
    "time_quarter" => "h, Δt/4",
    "space_half" => "h/2, Δt/4",
    "refined" => "h/2, Δt/2",
)

const MATCHED_COMPARISONS = vcat(
    [
        TOML.parsefile(joinpath(REFERENCE_DIR, "matched", "$(solver).toml"))["series"]
        for solver in ("fds", "gpyro", "thermakin")
    ]...,
)

function add_matched_comparisons!(case_id, curves, accuracy_path)
    entries = filter(s -> s["case"] == case_id, MATCHED_COMPARISONS)
    isempty(entries) && return curves
    audit = Any[]
    for entry in entries
        ci = findfirst(c -> string(c.name) == entry["curve"], curves)
        ci === nothing && error("$(case_id): missing comparison curve $(entry["curve"])")
        curve = curves[ci]
        si = findfirst(s -> !s.crosscode && s.label == entry["reference"], curve.series)
        si === nothing &&
            error("$(case_id): missing analytical reference $(entry["reference"])")
        exact = curve.series[si].exact
        exact === nothing &&
            error("$(case_id): comparison requires an analytical reference")
        for level in MATCHED_LEVELS
            prefix = level == "base" ? "matched/" : "matched/$(level)/"
            if level in MATCHED_COARSE &&
               !isfile(joinpath(REFERENCE_DIR, prefix * entry["path"]))
                continue
            end
            indices, x, y = reference_columns(
                prefix * entry["path"],
                "series",
                "x",
                "y";
                header_rows = 1,
            )
            keep = findall(==(entry["index"]), indices)
            isempty(keep) && error("$(case_id): missing retained comparison samples")
            xx, yy = x[keep], y[keep]
            all(isfinite, xx) && all(isfinite, yy) || error("nonfinite comparison data")
            errors = abs.(yy .- exact.(xx))
            all(isfinite, errors) || error("$(case_id): nonfinite analytical reference")
            push!(
                audit,
                (
                    entry["label"],
                    curve.name,
                    entry["reference"],
                    level,
                    length(keep),
                    maximum(errors),
                    sqrt(sum(abs2, errors) / length(errors)),
                ),
            )
            if level == "base"
                push!(
                    curve.series,
                    CurveSeries(
                        label = entry["label"],
                        x = xx,
                        y = yy,
                        exact = exact,
                        crosscode = true,
                    ),
                )
            end
        end
    end
    write_csv(
        accuracy_path,
        ["code", "curve", "reference", "resolution", "samples", "Linf", "RMS"],
        audit,
    )
    return curves
end

has_matched_comparison(case_id) = any(s -> s["case"] == case_id, MATCHED_COMPARISONS)

# --- cost at matched accuracy -------------------------------------------------
# Times of every matched external run and of this code over its ladder, measured
# under scripts/timing/README.md and committed in reference/matched/timing/
# (scripts/timing/export_dataset.py). External errors are the analytical errors
# add_matched_comparisons! records; this code's are on the same curves and
# references, at its own output points.

const MATCHED_TIMING_DIR = joinpath(REFERENCE_DIR, "matched", "timing")
# Errors at or below this count as exact when matching accuracy: both runs are
# then at round-off for the quantities compared.
const ROUNDOFF_ERROR = 1e-12

# Rows of a CSV with a header line; fields may be double-quoted.
function _csv_rows(path)
    split_line(s) = begin
        out, buf, quoted = String[], IOBuffer(), false
        for c in s
            if c == '"'
                quoted = !quoted
            elseif c == ',' && !quoted
                push!(out, String(take!(buf)))
            else
                print(buf, c)
            end
        end
        push!(out, String(take!(buf)))
        out
    end
    lines = filter(!isempty, readlines(path))
    header = split_line(lines[1])
    return [Dict(zip(header, split_line(l))) for l in lines[2:end]]
end

const SOLVER_LABELS =
    Dict("fds" => "FDS 6.11.1", "gpyro" => "Gpyro V2", "thermakin" => "ThermaKin2Ds 7")
_solver_key(entry) = String(first(split(entry["path"], '/')))
_run_key(entry) = splitext(basename(entry["path"]))[1]

"""
    matched_costs(case_id, accuracy_path)

Per matched external run of `case_id` (solver and native deck): its
resolutions as (level, worst error, time) and this code's configurations as
(cells, tolerance multiplier, worst error, time), the worst error taken over
the curves and references that run is scored on; and the table rows: for each
external resolution, the external code's cheapest run and this code's cheapest
configuration at least as accurate on every one of them, and the ratio of the
two. `nothing` without a timing dataset or matched runs.
"""
function matched_costs(case_id, accuracy_path)
    entries = filter(s -> s["case"] == case_id, MATCHED_COMPARISONS)
    (isempty(entries) || !isdir(MATCHED_TIMING_DIR) || !isfile(accuracy_path)) &&
        return nothing
    acc = Dict(
        (r["code"], r["curve"], r["reference"], r["resolution"]) =>
            parse(Float64, r["Linf"]) for r in _csv_rows(accuracy_path)
    )
    ext = Dict(
        (r["solver"], r["case"], r["level"]) => parse(Float64, r["compared_s"]) for
        r in _csv_rows(joinpath(MATCHED_TIMING_DIR, "external.csv"))
    )
    ours = Dict{Tuple{Int,Int,String,String},Tuple{Float64,Float64,String}}()
    for r in _csv_rows(joinpath(MATCHED_TIMING_DIR, "pyrolysis.csv"))
        r["case"] == case_id || continue
        key = (
            parse(Int, r["n_cells"]),
            parse(Int, r["tolerance"]),
            r["curve"],
            r["reference"],
        )
        ours[key] = (parse(Float64, r["Linf"]), parse(Float64, r["median_s"]), r["variant"])
    end
    runs = Dict{Tuple{String,String},Vector{Dict{String,Any}}}()
    for e in entries
        push!(get!(runs, (_solver_key(e), _run_key(e)), Dict{String,Any}[]), e)
    end
    per_solver = Dict{String,Int}()
    for (s, _) in keys(runs)
        per_solver[s] = get(per_solver, s, 0) + 1
    end
    panels, rows = NamedTuple[], Any[]
    for key in sort(collect(keys(runs)))
        es = runs[key]
        solver, deck = key
        refs = unique([(e["curve"], e["reference"]) for e in es])
        name =
            per_solver[solver] > 1 ?
            SOLVER_LABELS[solver] * ", " * join(last.(refs), ", ") : SOLVER_LABELS[solver]
        external = NamedTuple[]
        for level in MATCHED_LEVELS
            haskey(ext, (solver, deck, level)) || continue
            errs =
                [get(acc, (e["label"], e["curve"], e["reference"], level), NaN) for e in es]
            any(isnan, errs) && continue
            push!(
                external,
                (; level, error = maximum(errs), time = ext[(solver, deck, level)]),
            )
        end
        ladder = NamedTuple[]
        for (n, tol) in sort(unique((k[1], k[2]) for k in keys(ours)))
            all(r -> haskey(ours, (n, tol, r...)), refs) || continue
            vals = [ours[(n, tol, r...)] for r in refs]
            # one solve per variant (V-CPL-03's paths are separate problems)
            time = sum(Dict(v[3] => v[2] for v in vals) |> values)
            push!(ladder, (; n, tol, error = maximum(first.(vals)), time))
        end
        (isempty(external) || isempty(ladder)) && continue
        push!(panels, (; solver, name, refs, external, ladder))
        # The coarsest configuration: fewest cells at the loosest tolerance.
        coarsest = (minimum(p.n for p in ladder), maximum(p.tol for p in ladder))
        for x in external
            target = max(x.error, ROUNDOFF_ERROR)
            # each code at its own cheapest run at least this accurate
            theirs = minimum(q.time for q in external if q.error <= target)
            ok = filter(p -> p.error <= target, ladder)
            if isempty(ok)
                push!(
                    rows,
                    (name, x.level, x.error, x.time, theirs, "not reached", NaN, ""),
                )
            else
                best = ok[argmin([p.time for p in ok])]
                bound = (best.n, best.tol) == coarsest
                push!(
                    rows,
                    (
                        name,
                        x.level,
                        x.error,
                        x.time,
                        theirs,
                        best.tol == 1 ? string(best.n) : "$(best.n), $(best.tol)× tol",
                        best.time,
                        (bound ? "≥ " : "") *
                        string(round(theirs / best.time; sigdigits = 3)),
                    ),
                )
            end
        end
    end
    isempty(panels) && return nothing
    return (; panels, rows)
end
