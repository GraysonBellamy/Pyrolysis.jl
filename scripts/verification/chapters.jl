# =============================================================================
# Guide chapter generator: one markdown chapter per verification case family.
# =============================================================================
#
# Writes `<outdir>/<base id>.md` (e.g. V-HC-02.md covering members a–d) plus
# `<outdir>/index.md` (methodology + machine-written dashboard), in the FDS
# Verification Guide shape: problem statement → figures overlaying expected vs
# computed → error/order statement. The problem statement is the case file's
# leading comment banner — single source, so the chapter can never drift from
# what the test actually runs. `docs/make.jl` copies `<outdir>` into
# `docs/src/guide/verification/` at build time.
#
# Validation (same hard-fail invariant as docs/generate_guide.jl): a stray `$`
# outside code would render as mangled prose under Documenter, so it aborts
# chapter generation.

using Printf
using Dates

"Family id: member ids differ only by a trailing letter (V-HC-02a → V-HC-02)."
base_case_id(id::String) = replace(id, r"(?<=\d)[a-z]+$" => "")

_fmt(x::Float64) = isnan(x) ? "—" : @sprintf("%.4g", x)
_md_cell(s::AbstractString) = replace(s, "|" => "\\|", "\n" => " ")

function _md_table(io, header::Vector{String}, rows)
    println(io, "| ", join(_md_cell.(header), " | "), " |")
    println(io, "|", repeat(" --- |", length(header)))
    for row in rows
        println(io, "| ", join(_md_cell.(string.(collect(row))), " | "), " |")
    end
    println(io)
end

"""
Leading comment banner of a case file, comment markers stripped, horizontal
rules and the redundant `V-XX-NN — title` line removed (the chapter supplies
its own heading).
"""
function case_banner(case_file::String)
    banner = String[]
    for line in readlines(case_file)
        startswith(line, "#") || break
        push!(banner, replace(line, r"^# ?" => ""))
    end
    filter!(l -> !occursin(r"^=+$", l), banner)
    while !isempty(banner) &&
          (occursin(r"^V-[A-Z]+-\d+", banner[1]) || isempty(strip(banner[1])))
        popfirst!(banner)
    end
    return strip(join(banner, "\n"))
end

function validate_chapter(label::String, text::String)
    in_fence = false
    for (ln, line) in enumerate(split(text, '\n'))
        if startswith(lstrip(line), "```")
            in_fence = !in_fence
            continue
        end
        in_fence && continue
        startswith(line, "    ") && continue          # indented code block
        stripped = replace(line, r"`[^`]*`" => "")     # inline code spans
        occursin('$', stripped) && error(
            "$label:$ln: stray \$ outside code — chapters must not emit " *
            "TeX-style math (use unicode or fenced blocks in the case banner)",
        )
    end
end

function _write_member!(io, r::CaseReportResult)
    case = r.case
    n_finest = r.ladder_ns[end]
    qois = isempty(r.qois) ? QoIResult[] : r.qois[end]

    println(io, "### Quantities of interest (n = $(n_finest))\n")
    if isempty(qois)
        println(io, "_QoI evaluation unavailable at this resolution._\n")
    else
        _md_table(
            io,
            ["QoI", "value", "exact", "error", "tolerance", "within tol", "provenance"],
            [
                (
                    q.name,
                    _fmt(q.value),
                    _fmt(q.exact),
                    _fmt(q.error),
                    _fmt(q.tol),
                    q.error < q.tol ? "yes" : "**no**",
                    q.comment,
                ) for q in qois
            ],
        )
    end

    # Figures: base results inline; cross-code material gets its own narrated
    # section below rather than being dumped in sequence.
    _is_vs_err(k) = endswith(string(k), "_vs_codes_error")
    _is_vs(k) = endswith(string(k), "_vs_codes") && !_is_vs_err(k)
    _is_cc(k) = _is_vs(k) || _is_vs_err(k) || k === :solve_time
    for (kind, path) in r.figures
        _is_cc(kind) && continue
        println(io, "![$(case.id) $(kind)]($(path))\n")
    end

    cc_err = [(k, p) for (k, p) in r.figures if _is_vs_err(k)]
    cc_over = [(k, p) for (k, p) in r.figures if _is_vs(k)]
    cc_time = [(k, p) for (k, p) in r.figures if k === :solve_time]
    if !(isempty(cc_err) && isempty(cc_time))
        println(io, "### Comparison with other codes\n")
        refs = _reference_timings(case.id)
        codes_txt =
            isempty(refs) ? "the reference codes" :
            join(("$(t.label) at $(t.res)" for t in refs), " and ")
        println(
            io,
            """
            The same case was solved with $(codes_txt); decks, outputs, and
            run provenance are committed under `test/verification/reference/`.
            Each code's error against the same exact solution is drawn below
            on a log scale, muted gray behind this solver's series — the
            signed linear-scale panel above shows where the error lives, this
            one compares magnitudes across codes.
            """,
        )
        for (kind, path) in cc_err
            println(io, "![$(case.id) $(kind)]($(path))\n")
        end
        if !isempty(cc_time)
            println(
                io,
                """
                Wall times at every ladder rung against the reference runs.
                Resolutions and simulated spans differ where noted (details in
                `reference/timings.csv`), so cross-code timings are indicative
                rather than a controlled benchmark; rungs at a matched
                resolution are directly comparable.
                """,
            )
            for (kind, path) in cc_time
                println(io, "![$(case.id) $(kind)]($(path))\n")
            end
        end
        if !isempty(cc_over)
            println(
                io,
                "Solution overlays including the other codes' points: ",
                join(("[$(k)]($(p))" for (k, p) in cc_over), ", "),
                ".\n",
            )
        end
    end

    if !isempty(r.norms) && !isempty(keys(r.norms[1]))
        println(io, "### Convergence\n")
        norm_names = collect(keys(r.norms[1]))
        _md_table(
            io,
            vcat(["n_cells", "h", "wall (s)"], string.(norm_names)),
            [
                vcat(
                    Any[r.ladder_ns[i], _fmt(1.0 / r.ladder_ns[i]), _fmt(r.ladder_wall[i])],
                    Any[_fmt(getproperty(r.norms[i], nm)) for nm in norm_names],
                ) for i in eachindex(r.ladder_ns)
            ],
        )
        expected = report_expected_order(case)
        if r.observed_order !== nothing
            print(io, "Observed order $( _fmt(r.observed_order)) ($(norm_names[1]))")
            expected === nothing || print(io, ", expected $(expected)")
            println(io, ".\n")
        end
    end

    println(io, "### Solver configuration\n")
    _md_table(
        io,
        ["setting", "value"],
        vcat(
            [(
                "integrator",
                haskey(r.solve_kwargs, :alg) ? string(r.solve_kwargs.alg) :
                "KenCarp4 (default)",
            )],
            [(string(k), repr(v)) for (k, v) in pairs(r.solve_kwargs)],
        ),
    )
end

"""
Write one chapter per case family plus `index.md`. `results` come from
`run_case_report` (non-reported cases already filtered out); `all_cases` is
the full built list, used to note report-mode opt-outs in the index.
"""
function write_chapters(
    results::Vector{CaseReportResult},
    all_cases::Vector{VerificationCase},
    outdir::String;
    write_index::Bool = true,
)
    families = Dict{String,Vector{CaseReportResult}}()
    order = String[]
    for r in results
        b = base_case_id(r.case.id)
        haskey(families, b) || (families[b] = CaseReportResult[]; push!(order, b))
        push!(families[b], r)
    end
    sort!(order)

    for base in order
        members = families[base]
        lead = members[1].case
        io = IOBuffer()
        println(io, "# $(base) — $(lead.title)\n")
        println(io, "**Tags:** ", join(("`$(t)`" for t in lead.tags), ", "), "\n")
        println(io, "**References:**\n")
        for ref in lead.references
            println(io, "- ", ref)
        end
        println(io)
        banner_file = joinpath(ROOT, "test", "verification", "cases", "$(base).jl")
        banner = isfile(banner_file) ? case_banner(banner_file) : ""
        if !isempty(banner)
            println(io, "## Problem statement\n")
            println(io, banner, "\n")
        end
        for r in members
            if length(members) > 1 || r.case.id != base
                println(io, "## $(r.case.id) — $(r.case.title)\n")
            end
            haskey(r.case.report, :note) && println(io, "_$(r.case.report.note)_\n")
            _write_member!(io, r)
        end
        text = String(take!(io))
        validate_chapter("$(base).md", text)
        write(joinpath(outdir, "$(base).md"), text)
    end

    # Dashboard index — only rewritten when this run covered every case; a
    # partial --cases run would otherwise drop the unselected rows.
    write_index || return nothing
    io = IOBuffer()
    println(
        io,
        """
        # Verification suite report

        Exact-solution verification of Pyrolysis.jl: every case solves a
        problem with a closed-form (or manufactured) reference solution and
        reports error norms, pinned tolerances with provenance, and observed
        convergence orders. Methodology, case selection, and the dual-mode
        (CI assert / report) architecture are described in the Technical
        Reference Guide, §18.1–§18.2; the assert-mode twin of
        this report runs in `Pkg.test()` from the same case definitions
        (`test/verification/`).

        Figures follow the FDS Verification Guide convention: lines are the
        exact/reference solution, markers are the solver. Tolerances are
        pinned at ~3× the observed error at pin time (the provenance column).

        Where a case has a published analog in another solver, it carries a
        "Comparison with other codes" section: FDS 6.11.0 and Gpyro 0.8200
        were run locally on the same problems (decks and outputs committed
        with provenance in `test/verification/reference/`), and their errors
        against the same exact solutions are overlaid — muted gray, log
        scale — together with a wall-time comparison at every resolution
        rung.

        Generated $(Dates.format(Dates.now(), "yyyy-mm-dd")) by
        `scripts/verification/run_reports.jl`. Regenerate with:
        `julia --project=scripts/verification scripts/verification/run_reports.jl`
        """,
    )
    rows = Any[]
    for base in order, r in families[base]
        qois = isempty(r.qois) ? QoIResult[] : r.qois[end]
        worst = isempty(qois) ? NaN : maximum(q.error / q.tol for q in qois)
        push!(
            rows,
            (
                "[$(r.case.id)]($(base).md)",
                r.case.title,
                join(string.(r.case.tags), ", "),
                isempty(qois) ? "—" :
                "$(count(q -> q.error < q.tol, qois))/$(length(qois))",
                _fmt(worst),
                r.observed_order === nothing ? "—" : _fmt(r.observed_order),
            ),
        )
    end
    _md_table(
        io,
        ["case", "title", "tags", "QoIs within tol", "worst error/tol", "observed order"],
        rows,
    )
    reported = Set(r.case.id for r in results)
    skipped = [c.id for c in all_cases if !(c.id in reported)]
    if !isempty(skipped)
        println(
            io,
            "Cases without report-mode output (no `report` metadata yet): ",
            join(skipped, ", "),
            ".",
        )
    end
    text = String(take!(io))
    validate_chapter("index.md", text)
    write(joinpath(outdir, "index.md"), text)
    return nothing
end
