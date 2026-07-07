# =============================================================================
# Report-mode harness: ladder execution + artifact generation per case.
# =============================================================================
#
# Consumes `VerificationCase`s (test/verification/casedefs.jl) and writes the
# standard artifact set into `<outdir>/assets/<case id>/`:
#   figures  — overlay, pointwise error, convergence, sweeps (render.jl)
#   CSVs     — one per curve (long format), resolution ladder, QoIs at the
#              finest rung, ladder QoI errors, solver configuration
# Unlike assert mode, nothing here asserts: QoI errors are recorded against
# their pinned tolerances (within-tol yes/no) so coarse ladder rungs may
# legitimately exceed them.

using Printf
using Dates
import Pkg

# --- CSV helpers (hand-rolled: values are numbers + a few quoted strings; a
# --- CSV package dependency is not warranted) --------------------------------

_csv_cell(x::AbstractString) = '"' * replace(x, "\"" => "\"\"") * '"'
_csv_cell(x::Symbol) = _csv_cell(string(x))
_csv_cell(x::Bool) = x ? "true" : "false"
_csv_cell(x::Real) = string(x)
_csv_cell(::Nothing) = ""
_csv_cell(::Missing) = ""

function write_csv(path::String, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(map(_csv_cell, header), ","))
        for row in rows
            println(io, join(map(_csv_cell, collect(row)), ","))
        end
    end
    return path
end

# --- report-key resolution (fallbacks documented on `VerificationCase`) ------

function report_solve_kwargs(case::VerificationCase)
    conv = case.ci.convergence
    base =
        conv === nothing ? case.solve_kwargs : merge(case.solve_kwargs, conv.solve_kwargs)
    return merge(base, get(case.report, :solve_kwargs, (;)))
end

function report_error_norms_fn(case::VerificationCase)
    haskey(case.report, :error_norms) && return case.report.error_norms
    conv = case.ci.convergence
    conv === nothing && return nothing
    return sol -> (norm = conv.error_fn(sol),)
end

report_expected_order(case::VerificationCase) =
    haskey(case.report, :expected_order) ? Float64(case.report.expected_order) :
    case.ci.convergence === nothing ? nothing : case.ci.convergence.expected_order

"""
Other codes' wall times for a case, from `reference/timings.csv`
(case,code,label,res,wall_s,note — plain commas, notes use semicolons; `res`
is the resolution text shown next to the code's label on the figure).
Empty when the file or the case's rows are absent.
"""
function _reference_timings(case_id::AbstractString)
    path = joinpath(REFERENCE_DIR, "timings.csv")
    isfile(path) || return NamedTuple[]
    out = NamedTuple[]
    for ln in Iterators.drop(readlines(path), 1)
        isempty(strip(ln)) && continue
        f = split(ln, ','; limit = 6)
        length(f) == 6 || continue
        strip(f[1]) == case_id || continue
        push!(
            out,
            (
                code = String(strip(f[2])),
                label = String(strip(f[3])),
                res = String(strip(f[4])),
                wall = parse(Float64, f[5]),
                note = String(strip(f[6])),
            ),
        )
    end
    return out
end

# --- per-case execution -------------------------------------------------------

"""
Result of one report run of one case, consumed by the chapter generator.
`qois` and `norms` align with `ladder_ns`; `figures` are paths relative to the
report `outdir` (chapter markdown links use them verbatim).
"""
Base.@kwdef struct CaseReportResult
    case::VerificationCase
    ladder_ns::Vector{Int}
    ladder_wall::Vector{Float64}
    qois::Vector{Vector{QoIResult}}
    norms::Vector{NamedTuple}
    observed_order::Union{Nothing,Float64}
    figures::Vector{Pair{Symbol,String}}   # kind => relative path
    solve_kwargs::NamedTuple
end

_relfig(outdir, path) = relpath(path, outdir)

"""
Run `case` through its report ladder and write every artifact under
`<outdir>/assets/<id>/`. Returns a `CaseReportResult`, or `nothing` when the
case opts out of report mode (`report = (;)`).
"""
function run_case_report(case::VerificationCase, outdir::String)
    isempty(pairs(case.report)) && return nothing
    assets = joinpath(outdir, "assets", case.id)
    mkpath(assets)
    apath(suffix) = joinpath(assets, "$(case.id)_$(suffix)")

    ns = collect(Int, get(case.report, :n_cells_ladder, [case.ci.n_cells]))
    skw = report_solve_kwargs(case)
    normfn = report_error_norms_fn(case)

    # Warm-up solve, discarded: the first solve of a case pays Julia JIT
    # compilation for its concrete problem types, which would otherwise be
    # billed to the coarsest rung's wall time (and dwarf the solve-time
    # comparison). Errors are ignored here — the timed loop re-raises them.
    try
        _case_solve(case, case.build(ns[1]); skw...)
    catch
    end

    walls = Float64[]
    qois_per_n = Vector{QoIResult}[]
    norms_per_n = NamedTuple[]
    finest_sol = nothing
    for n in ns
        local sol
        wall = @elapsed sol = _case_solve(case, case.build(n); skw...)
        sol.retcode == :Success ||
            error("$(case.id): solve failed at n = $(n) (retcode $(sol.retcode))")
        push!(walls, wall)
        # QoI closures may assume the CI resolution (fit windows etc.); a
        # coarse rung that breaks them is recorded as empty, not fatal.
        qs = try
            case.qois(sol)
        catch err
            @warn "QoI evaluation failed on ladder rung" case.id n err
            QoIResult[]
        end
        push!(qois_per_n, qs)
        push!(norms_per_n, normfn === nothing ? (;) : normfn(sol))
        finest_sol = sol
        @info "report rung done" case.id n wall
    end

    figures = Pair{Symbol,String}[]

    # Curves: overlay + pointwise error figures, one long-format CSV per curve.
    # Curves carrying cross-code series additionally get `_vs_codes` variants:
    # the base figures stay ours-vs-exact only (so our error is discernible),
    # the variants add the other codes' data as a muted context layer.
    curves = haskey(case.report, :curves) ? case.report.curves(finest_sol) : ReportCurve[]
    for curve in curves
        p = render_overlay(curve, apath("$(curve.name).png"))
        push!(figures, curve.name => _relfig(outdir, p))
        pe = render_pointwise_error(curve, apath("$(curve.name)_error.png"))
        pe === nothing || push!(figures, Symbol(curve.name, :_error) => _relfig(outdir, pe))
        if any(s -> s.crosscode, curve.series)
            pc =
                render_overlay(curve, apath("$(curve.name)_vs_codes.png"); crosscode = true)
            push!(figures, Symbol(curve.name, :_vs_codes) => _relfig(outdir, pc))
            pec = render_pointwise_error(
                curve,
                apath("$(curve.name)_vs_codes_error.png");
                crosscode = true,
            )
            pec === nothing ||
                push!(figures, Symbol(curve.name, :_vs_codes_error) => _relfig(outdir, pec))
        end
        rows = Any[]
        for s in curve.series, i in eachindex(s.x)
            push!(
                rows,
                (s.label, s.x[i], s.y[i], s.exact === nothing ? nothing : s.exact(s.x[i])),
            )
        end
        write_csv(
            apath("$(curve.name).csv"),
            ["series", curve.xlabel, curve.ylabel, "exact"],
            rows,
        )
    end

    # Solve-time comparison, when other codes' reference runs exist for this
    # case (reference/timings.csv; provenance in reference/*/README.md).
    ref_times = _reference_timings(case.id)
    if !isempty(ref_times)
        p = render_solve_time(case.id, ns, walls, ref_times, apath("solve_time.png"))
        push!(figures, :solve_time => _relfig(outdir, p))
    end

    # Resolution ladder: CSV always; convergence figure when ≥ 2 usable rungs.
    observed = nothing
    if !isempty(norms_per_n) && !isempty(keys(norms_per_n[1]))
        hs = 1.0 ./ ns
        first_norm = keys(norms_per_n[1])[1]
        errs = Float64[getproperty(nt, first_norm) for nt in norms_per_n]
        if count(>(0), errs) >= 2
            keep = errs .> 0
            observed = fit_order(hs[keep], errs[keep])
        end
        norm_names = collect(keys(norms_per_n[1]))
        write_csv(
            apath("ladder.csv"),
            vcat(["n_cells", "h", "wall_s"], string.(norm_names)),
            [
                vcat(
                    Any[ns[i], hs[i], walls[i]],
                    Any[getproperty(norms_per_n[i], nm) for nm in norm_names],
                ) for i in eachindex(ns)
            ],
        )
        if length(ns) >= 2
            p = render_convergence(
                "$(case.id): spatial convergence",
                hs,
                norms_per_n,
                report_expected_order(case),
                apath("convergence.png"),
            )
            p === nothing || push!(figures, :convergence => _relfig(outdir, p))
        end
    end

    # Ladder QoI errors (long format) + finest-rung QoI table.
    write_csv(
        apath("ladder_qois.csv"),
        ["n_cells", "qoi", "value", "exact", "error", "tol", "within_tol"],
        [
            (ns[i], q.name, q.value, q.exact, q.error, q.tol, q.error < q.tol) for
            i in eachindex(ns) for q in qois_per_n[i]
        ],
    )
    write_csv(
        apath("qois.csv"),
        ["qoi", "value", "exact", "error", "tol", "within_tol", "comment"],
        [
            (q.name, q.value, q.exact, q.error, q.tol, q.error < q.tol, q.comment) for
            q in (isempty(qois_per_n) ? QoIResult[] : qois_per_n[end])
        ],
    )

    # Temporal sweep (fixed-dt time-order study, e.g. V-MMS-01). The error is
    # measured against a same-mesh reference solve at dt_min/8, NOT against
    # the exact solution — at these mesh resolutions the spatial error floor
    # dominates the exact-solution comparison and would read as order 0.
    if haskey(case.report, :temporal_sweep)
        ts = case.report.temporal_sweep
        dts = collect(Float64, ts.dts)
        # dt implies a fixed step: solve_mms sets adaptive=false itself
        ref = _case_solve(case, case.build(ts.n_cells); skw..., dt = minimum(dts) / 8)
        ref.retcode == :Success ||
            error("$(case.id): temporal-sweep reference solve failed")
        tnorms = NamedTuple[]
        for dt in dts
            sol = _case_solve(case, case.build(ts.n_cells); skw..., dt = dt)
            sol.retcode == :Success ||
                error("$(case.id): temporal-sweep solve failed at dt = $(dt)")
            e = sol.T[:, end] .- ref.T[:, end]
            push!(
                tnorms,
                (
                    L2_T_vs_ref = sqrt(sum(abs2, e) / length(e)),
                    Linf_T_vs_ref = maximum(abs, e),
                ),
            )
        end
        norm_names = collect(keys(tnorms[1]))
        write_csv(
            apath("temporal_sweep.csv"),
            vcat(["dt"], string.(norm_names)),
            [
                vcat(Any[dts[i]], Any[getproperty(tnorms[i], nm) for nm in norm_names]) for
                i in eachindex(dts)
            ],
        )
        p = render_convergence(
            "$(case.id): temporal convergence (n = $(ts.n_cells))",
            dts,
            tnorms,
            nothing,
            apath("convergence_dt.png");
            xlabel = "Δt (s)",
        )
        p === nothing || push!(figures, :convergence_dt => _relfig(outdir, p))
    end

    # Parameter sweep (e.g. V-CPL-02's activation-energy ladder).
    if haskey(case.report, :param_sweep)
        ps = case.report.param_sweep
        values = collect(Float64, ps.values)
        measures = NamedTuple[]
        for v in values
            sol = _case_solve(case, ps.build(case.ci.n_cells, v); skw...)
            sol.retcode == :Success ||
                error("$(case.id): param-sweep solve failed at $(ps.name) = $(v)")
            # measure may take the swept value too (when the reference
            # solution itself depends on it, e.g. V-CPL-05's β sweep).
            push!(
                measures,
                applicable(ps.measure, sol, v) ? ps.measure(sol, v) : ps.measure(sol),
            )
        end
        mnames = collect(keys(measures[1]))
        write_csv(
            apath("param_$(ps.name).csv"),
            vcat([String(ps.name)], string.(mnames)),
            [
                vcat(Any[values[i]], Any[getproperty(measures[i], nm) for nm in mnames]) for
                i in eachindex(values)
            ],
        )
        p = render_param_sweep(
            "$(case.id): $(ps.name) sweep (n = $(case.ci.n_cells))",
            String(ps.name),
            values,
            measures,
            apath("param_$(ps.name).png"),
        )
        p === nothing || push!(figures, Symbol(:param_, ps.name) => _relfig(outdir, p))
    end

    # Solver configuration, for reproducibility of the published artifacts.
    pyro_version = string(Pkg.TOML.parsefile(joinpath(ROOT, "Project.toml"))["version"])
    write_csv(
        apath("config.csv"),
        ["key", "value"],
        vcat(
            [
                ("integrator", haskey(skw, :alg) ? string(skw.alg) : "KenCarp4 (default)"),
                ("julia", string(VERSION)),
                ("Pyrolysis", pyro_version),
                ("n_cells_ladder", join(ns, " ")),
            ],
            [(string(k), repr(v)) for (k, v) in pairs(skw)],
        ),
    )

    return CaseReportResult(
        case = case,
        ladder_ns = ns,
        ladder_wall = walls,
        qois = qois_per_n,
        norms = norms_per_n,
        observed_order = observed,
        figures = figures,
        solve_kwargs = skw,
    )
end
