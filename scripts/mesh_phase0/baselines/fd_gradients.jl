# Finite-difference gradient smoothness of the V-BO-01 QoIs through park events.
#
# The V-BO-01 builder (test/verification/cases/V-BO-01.jl) exposes the PMMA
# absorption coefficient `pmma_alpha` (m² kg⁻¹, default 2870/ρ_PMMA). This
# script perturbs it by relative amounts ε ∈ {0, ±1e-3, ±2e-3, ±4e-3, ±8e-3} at
# the case's CI cell count, records the QoIs, the park sequence and the step
# statistics of every run, then forms central differences at
# h ∈ {1e-3, 2e-3, 4e-3, 8e-3} and a spread metric per QoI.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/fd_gradients.jl
# With `--analyze-only` the runs are read back from fd_gradients.csv and only
# the derivative CSV and the markdown are regenerated.
#
# Outputs (scripts/mesh_phase0/results/):
#   fd_gradients.csv             one row per run (written incrementally)
#   fd_gradient_derivatives.csv  central-difference estimates per QoI and h
#   fd_gradients.md              tables and findings

using Pyrolysis
using Logging
using Printf
using Statistics: mean
using OrdinaryDiffEq: DiscreteCallback, SciMLBase

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULTS = joinpath(ROOT, "scripts", "mesh_phase0", "results")
mkpath(RESULTS)
const RUNS_CSV = joinpath(RESULTS, "fd_gradients.csv")
const DERIV_CSV = joinpath(RESULTS, "fd_gradient_derivatives.csv")
const MD = joinpath(RESULTS, "fd_gradients.md")

include(joinpath(ROOT, "test", "verification", "casedefs.jl"))

const CASE = only(filter(c -> c.id == "V-BO-01", load_verification_cases()))
const N_PMMA = CASE.ci.n_cells
const P0 = 2870.0 / VBO_RHO_PMMA          # builder default of `pmma_alpha`
const EPS = [0.0, 1e-3, -1e-3, 2e-3, -2e-3, 4e-3, -4e-3, 8e-3, -8e-3]
const H = [1e-3, 2e-3, 4e-3, 8e-3]
const TIMEOUT_S = 600.0
const MAXITERS = 50_000

# ----------------------------------------------------------------------------
# Logger: captures the solver's "Parked cell i at t=… s" messages (the parked
# live index in event order), counts surface-Newton warnings, forwards warnings
# and errors, and drops the Info noise.
# ----------------------------------------------------------------------------
struct ParkCapture{L<:Logging.AbstractLogger} <: Logging.AbstractLogger
    inner::L
    cells::Vector{Int}
    times::Vector{Float64}
    surface_warnings::Base.RefValue{Int}
end
ParkCapture(inner) = ParkCapture(inner, Int[], Float64[], Ref(0))
Logging.min_enabled_level(::ParkCapture) = Logging.Info
Logging.shouldlog(::ParkCapture, args...) = true
Logging.catch_exceptions(::ParkCapture) = false
function Logging.handle_message(l::ParkCapture, level, message, args...; kwargs...)
    msg = string(message)
    m = match(r"^Parked cell (\d+) at t=([0-9.eE+-]+) s", msg)
    if m !== nothing
        push!(l.cells, parse(Int, m[1]))
        push!(l.times, parse(Float64, m[2]))
        return nothing
    end
    if level >= Logging.Warn
        occursin(r"^Surface", msg) && (l.surface_warnings[] += 1)
        return Logging.handle_message(l.inner, level, message, args...; kwargs...)
    end
    return nothing
end

# Same step recorder as the case's `solve_fn` (clears the derivative
# discontinuity flag after every accepted step) plus a wall-clock guard.
function solve_vbo(problem; limit_s = TIMEOUT_S, maxiters = MAXITERS)
    t0 = time()
    hit = Ref(false)
    recorder = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            if time() - t0 > limit_s
                hit[] = true
                SciMLBase.terminate!(integ)
            end
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
    cap = ParkCapture(Logging.current_logger())
    sol = Logging.with_logger(cap) do
        solve(problem; progress = false, callback = recorder, maxiters, CASE.solve_kwargs...)
    end
    return sol, time() - t0, hit[], cap
end

# Continuous companions of the grid-quantized QoIs: parabolic peak through the
# three output samples around the discrete maximum, and the linear crossing
# time of the live thickness through 12.02 mm.
function interp_qois(sol)
    mlr = sol.extras.gas_generation_rate .* 1000
    t = sol.t
    i = argmax(mlr)
    if 1 < i < length(mlr)
        y0, y1, y2 = mlr[i-1], mlr[i], mlr[i+1]
        d = y0 - 2y1 + y2
        δ = d == 0 ? 0.0 : 0.5 * (y0 - y2) / d
        t_pk = t[i] + δ * (t[i+1] - t[i])
        pk = y1 - 0.25 * (y0 - y2) * δ
    else
        t_pk, pk = t[i], mlr[i]
    end
    th = sol.extras.thickness_history
    L_bo = VBO_L_KW + 0.02e-3
    j = findfirst(<(L_bo), th)
    t_bo = NaN
    if j !== nothing && j > 1
        t_bo = t[j-1] + (th[j-1] - L_bo) / (th[j-1] - th[j]) * (t[j] - t[j-1])
    end
    return (; peak_mlr_interp = pk, t_peak_interp = t_pk, t_bo_interp = t_bo)
end

const QOI_KEYS = [
    "peak MLR" => :peak_mlr,
    "time of peak MLR" => :t_peak,
    "burnout time (live thickness < 12.02 mm)" => :t_burnout,
    "mass loss of the gasified PMMA layer" => :mass_loss,
    "complete burnout: live thickness at t_end is the substrate" => :L_end,
    "mass ledger closure" => :ledger,
]
const QOI_NAMES = [
    :peak_mlr,
    :t_peak,
    :t_burnout,
    :mass_loss,
    :L_end,
    :ledger,
    :peak_mlr_interp,
    :t_peak_interp,
    :t_bo_interp,
    :t_park_first,
    :t_park_last,
    :t_park_mean,
]

const COLUMNS = [
    "run",
    "eps",
    "p",
    "retcode",
    "timed_out",
    "wall_s",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nonliniter",
    "nonlinfail",
    "n_parks",
    "surface_warnings",
    string.(QOI_NAMES)...,
    "park_cells",
    "park_times",
]

function run_one(run, ε)
    p = P0 * (1 + ε)
    problem = _vbo_problem(N_PMMA; pmma_alpha = p)
    row = Dict{String,Any}("run" => run, "eps" => ε, "p" => p)
    for c in COLUMNS[4:end]
        row[c] = NaN
    end
    row["park_cells"] = ""
    row["park_times"] = ""
    try
        sol, wall, hit, cap = solve_vbo(problem)
        st = sol.raw_solution.stats
        row["retcode"] = string(sol.retcode)
        row["timed_out"] = hit
        row["wall_s"] = wall
        row["steps"] = st.naccept
        row["rejects"] = st.nreject
        row["nf"] = st.nf
        row["njacs"] = st.njacs
        row["nonliniter"] = st.nnonliniter
        row["nonlinfail"] = st.nnonlinconvfail
        pt = sol.extras.park_times
        row["n_parks"] = length(pt)
        row["surface_warnings"] = cap.surface_warnings[]
        VBO_SURFACE_WARNINGS[CASE.id] = cap.surface_warnings[]
        if sol.retcode == :Success
            for q in CASE.qois(sol)
                for (name, key) in QOI_KEYS
                    name == q.name && (row[string(key)] = q.value)
                end
            end
            iq = interp_qois(sol)
            for k in keys(iq)
                row[string(k)] = iq[k]
            end
        end
        if !isempty(pt)
            row["t_park_first"] = pt[1]
            row["t_park_last"] = pt[end]
            row["t_park_mean"] = mean(pt)
        end
        row["park_cells"] = join(string.(cap.cells), ";")
        row["park_times"] = join((@sprintf("%.6f", x) for x in pt), ";")
    catch err
        row["retcode"] = "Exception: " * sprint(showerror, err)
        row["timed_out"] = false
    end
    return row
end

csvcell(x) = x isa AbstractString ? "\"" * replace(x, "\"" => "'") * "\"" : string(x)

function write_row(io, row)
    println(io, join((csvcell(row[c]) for c in COLUMNS), ","))
    flush(io)
end

# ----------------------------------------------------------------------------
# Runs
# ----------------------------------------------------------------------------
function parse_csv_line(line)
    fields = String[]
    buf = IOBuffer()
    inq = false
    for ch in line
        if ch == '"'
            inq = !inq
        elseif ch == ',' && !inq
            push!(fields, String(take!(buf)))
        else
            write(buf, ch)
        end
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_runs(path)
    lines = readlines(path)
    header = parse_csv_line(lines[1])
    header == COLUMNS || error("unexpected header in $path")
    out = Dict{String,Any}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = parse_csv_line(line)
        r = Dict{String,Any}()
        for (c, x) in zip(COLUMNS, f)
            r[c] =
                c in ("retcode", "park_cells", "park_times") ? x :
                c == "timed_out" ? parse(Bool, x) :
                c in (
                    "run",
                    "steps",
                    "rejects",
                    "nf",
                    "njacs",
                    "nonliniter",
                    "nonlinfail",
                    "n_parks",
                    "surface_warnings",
                ) ? (x == "NaN" ? NaN : parse(Int, x)) : parse(Float64, x)
        end
        push!(out, r)
    end
    return out
end

t_script = time()
rows = Dict{String,Any}[]
if "--analyze-only" in ARGS
    rows = read_runs(RUNS_CSV)
    println("read ", length(rows), " runs from ", RUNS_CSV)
else
    println(
        "V-BO-01 FD gradients: n_pmma = $N_PMMA, p0 = $P0, timeout $TIMEOUT_S s, maxiters $MAXITERS",
    )
    flush(stdout)
    let
        r = run_one(0, 0.0)
        println(
            "warmup: ",
            r["retcode"],
            " ",
            @sprintf("%.1f", r["wall_s"]),
            " s, parks=",
            r["n_parks"],
        )
        flush(stdout)
    end
    open(RUNS_CSV, "w") do io
        println(io, join(COLUMNS, ","))
        flush(io)
        for (k, ε) in enumerate(EPS)
            r = run_one(k, ε)
            push!(rows, r)
            write_row(io, r)
            println(
                @sprintf(
                    "run %d eps=%+.0e: %s %.2f s steps=%s parks=%s peak=%s t_peak=%s t_bo=%s",
                    k,
                    ε,
                    r["retcode"],
                    r["wall_s"],
                    r["steps"],
                    r["n_parks"],
                    r["peak_mlr"],
                    r["t_peak"],
                    r["t_burnout"]
                ),
            )
            flush(stdout)
        end
    end
end

# ----------------------------------------------------------------------------
# Central differences and spread
# ----------------------------------------------------------------------------
byeps = Dict(r["eps"] => r for r in rows)
value(ε, q) = (r = byeps[ε]; v = r[string(q)]; v isa Number ? Float64(v) : NaN)

deriv = Dict{Symbol,Vector{Float64}}()
open(DERIV_CSV, "w") do io
    println(io, "qoi,h,Q_plus,Q_minus,dQ_dlnp,dQ_dp,dQ_dlnp_over_Q0")
    for q in QOI_NAMES
        ds = Float64[]
        q0 = value(0.0, q)
        for h in H
            qp, qm = value(h, q), value(-h, q)
            d = (qp - qm) / (2h)
            push!(ds, d)
            println(io, "$q,$h,$qp,$qm,$d,$(d / P0),$(q0 == 0 ? NaN : d / abs(q0))")
        end
        deriv[q] = ds
    end
end

function spread(ds)
    good = filter(isfinite, ds)
    isempty(good) && return NaN
    m = mean(good)
    s = maximum(good) - minimum(good)
    return m == 0 ? (s == 0 ? 0.0 : Inf) : s / abs(m)
end

# Largest logarithmic sensitivity |dQ/dln p| / |Q(0)| over the four estimates.
function rel_sens(q)
    q0 = value(0.0, q)
    q0 == 0 && return NaN
    good = filter(isfinite, deriv[q])
    return isempty(good) ? NaN : maximum(abs, good) / abs(q0)
end

const INSENSITIVE = 1e-6      # |D|/|Q0| below this: derivative at round-off
classify(q) =
    rel_sens(q) < INSENSITIVE ? "insensitive" :
    spread(deriv[q]) < 0.1 ? "smooth" : spread(deriv[q]) < 1 ? "moderate" : "erratic"

# ----------------------------------------------------------------------------
# Park-sequence comparison between neighbouring ε
# ----------------------------------------------------------------------------
sorted_eps = sort(EPS)
parkcells(ε) = (s = byeps[ε]["park_cells"]; isempty(s) ? Int[] : parse.(Int, split(s, ";")))
parktimes(ε) =
    (s = byeps[ε]["park_times"]; isempty(s) ? Float64[] : parse.(Float64, split(s, ";")))

park_rows = String[]
n_pairs = 0
n_same_count = 0
n_same_order = 0
max_shift_all = 0.0
for (k, ε) in enumerate(sorted_eps)
    c, t = parkcells(ε), parktimes(ε)
    same_count = same_order = "-"
    shift = "-"
    if k > 1
        global n_pairs += 1
        cp, tp = parkcells(sorted_eps[k-1]), parktimes(sorted_eps[k-1])
        sc = length(c) == length(cp)
        so = sc && c == cp
        global n_same_count += sc
        global n_same_order += so
        same_count = string(sc)
        same_order = sc ? string(so, " (", count(c .!= cp), " positions differ)") : "false"
        if sc && !isempty(t)
            ms = maximum(abs.(t .- tp))
            shift = @sprintf("%.3f", ms)
            global max_shift_all = max(max_shift_all, ms)
        end
    end
    push!(
        park_rows,
        @sprintf(
            "| %+.0e | %d | %s | %s | %s | %s | %s |",
            ε,
            length(t),
            isempty(t) ? "-" : @sprintf("%.3f", t[1]),
            isempty(t) ? "-" : @sprintf("%.3f", t[end]),
            same_count,
            same_order,
            shift
        ),
    )
end

monotone(xs) = all(diff(xs) .>= 0) || all(diff(xs) .<= 0)
first_parks = [value(ε, :t_park_first) for ε in sorted_eps]
last_parks = [value(ε, :t_park_last) for ε in sorted_eps]

# ----------------------------------------------------------------------------
# Markdown
# ----------------------------------------------------------------------------
fmtv(x) =
    !isfinite(x) ? string(x) :
    (abs(x) >= 1e4 || (0 < abs(x) < 1e-2) ? @sprintf("%.4e", x) : @sprintf("%.5g", x))

open(MD, "w") do io
    println(io, "# V-BO-01 finite-difference gradient smoothness through park events\n")
    println(
        io,
        "Generated by `scripts/mesh_phase0/baselines/fd_gradients.jl`. Raw runs in ",
        "`fd_gradients.csv`, central differences in `fd_gradient_derivatives.csv`.\n",
    )
    println(
        io,
        "Perturbed parameter: PMMA absorption coefficient `pmma_alpha` of the V-BO-01 ",
        "builder, p = p0 (1 + ε), p0 = ",
        @sprintf("%.6g", P0),
        " m² kg⁻¹ (2870 m⁻¹ / ρ_PMMA). ",
        "Case solved at its CI resolution (n_pmma = $N_PMMA, ",
        length(_vbo_problem(N_PMMA).mesh.cells),
        " cells) with the case's `solve_kwargs` and step recorder; wall-clock guard ",
        "$TIMEOUT_S s, maxiters $MAXITERS. Wall times are indicative (shared machine).\n",
    )
    println(
        io,
        "`t_peak` and `t_burnout` are the case's QoIs on the 0.5 s output grid ",
        "(quantized); `*_interp` are continuous companions (parabolic peak, linear ",
        "thickness crossing); `t_park_*` are statistics of `sol.extras.park_times`.\n",
    )
    println(io, "## Runs\n")
    println(
        io,
        "| ε | retcode | timed out | wall s | steps | rejects | nlfail | parks | peak MLR | t_peak | t_burnout | mass loss % | L_end m | ledger | peak interp | t_peak interp | t_bo interp | t_park first | t_park last |",
    )
    println(
        io,
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    )
    for ε in sorted_eps
        r = byeps[ε]
        println(
            io,
            @sprintf(
                "| %+.0e | %s | %s | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                ε,
                r["retcode"],
                r["timed_out"],
                r["wall_s"],
                r["steps"],
                r["rejects"],
                r["nonlinfail"],
                r["n_parks"],
                fmtv(value(ε, :peak_mlr)),
                fmtv(value(ε, :t_peak)),
                fmtv(value(ε, :t_burnout)),
                fmtv(value(ε, :mass_loss)),
                fmtv(value(ε, :L_end)),
                fmtv(value(ε, :ledger)),
                fmtv(value(ε, :peak_mlr_interp)),
                fmtv(value(ε, :t_peak_interp)),
                fmtv(value(ε, :t_bo_interp)),
                fmtv(value(ε, :t_park_first)),
                fmtv(value(ε, :t_park_last))
            )
        )
    end
    println(io, "\n## Central differences D(h) = dQ/dln p = (Q(+h) − Q(−h)) / 2h\n")
    println(
        io,
        "`|D|/Q0` is the largest logarithmic sensitivity over the four h; below ",
        INSENSITIVE,
        " the derivative is at round-off and the spread is meaningless (class `insensitive`). ",
        "Classes: smooth = spread < 0.1, moderate = spread < 1, erratic otherwise.\n",
    )
    println(
        io,
        "| QoI | h = 1e-3 | h = 2e-3 | h = 4e-3 | h = 8e-3 | mean | spread (max−min)/|mean| | |D|/Q0 | |D(1e-3)−D(2e-3)| | |D(2e-3)−D(4e-3)| | |D(4e-3)−D(8e-3)| | class |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for q in QOI_NAMES
        d = deriv[q]
        println(
            io,
            "| $q | ",
            join(fmtv.(d), " | "),
            " | ",
            fmtv(mean(filter(isfinite, d))),
            " | ",
            fmtv(spread(d)),
            " | ",
            fmtv(rel_sens(q)),
            " | ",
            fmtv(abs(d[1] - d[2])),
            " | ",
            fmtv(abs(d[2] - d[3])),
            " | ",
            fmtv(abs(d[3] - d[4])),
            " | ",
            classify(q),
            " |",
        )
    end
    println(io, "\n## Park sequence per ε (compared with the next-lower ε)\n")
    println(
        io,
        "`same cell order` compares the sequence of parked live-cell indices logged by the solver.\n",
    )
    println(
        io,
        "| ε | parks | first park s | last park s | same count | same cell order | max park-time shift s |",
    )
    println(io, "|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), park_rows)
    println(io, "\n## Findings\n")
    cls = Dict(q => classify(q) for q in QOI_NAMES)
    bycls(c) = [q for q in QOI_NAMES if cls[q] == c]
    describe(qs) =
        isempty(qs) ? "none" :
        join((string(q, " (", fmtv(spread(deriv[q])), ")") for q in qs), ", ")
    println(
        io,
        "1. Of the $(n_pairs) neighbouring-ε pairs, $(n_same_count) have the same park count (",
        join(sort(unique([byeps[ε]["n_parks"] for ε in sorted_eps])), ", "),
        " parks) and $(n_same_order) the same parked-cell order; ",
        "the largest shift of a matched park time between neighbours is ",
        @sprintf("%.3f", max_shift_all),
        " s, the first park time is ",
        monotone(first_parks) ? "monotone" : "not monotone",
        " in ε (",
        @sprintf("%.2f", minimum(first_parks)),
        "–",
        @sprintf("%.2f", maximum(first_parks)),
        " s) and the last park time is ",
        monotone(last_parks) ? "monotone" : "not monotone",
        " in ε (",
        @sprintf("%.2f", minimum(last_parks)),
        "–",
        @sprintf("%.2f", maximum(last_parks)),
        " s).",
    )
    println(
        io,
        "2. Spread of the four central differences: smooth (< 0.1) for ",
        describe(bycls("smooth")),
        "; moderate for ",
        describe(bycls("moderate")),
        "; erratic (≥ 1, estimates that change sign or magnitude with h) for ",
        describe(bycls("erratic")),
        "; insensitive (|D|/Q0 < ",
        INSENSITIVE,
        ", derivative at round-off) for ",
        join(string.(bycls("insensitive")), ", "),
        ".",
    )
end
println("wrote ", RUNS_CSV, ", ", DERIV_CSV, ", ", MD)
println(@sprintf("total script wall time %.1f s", time() - t_script))
