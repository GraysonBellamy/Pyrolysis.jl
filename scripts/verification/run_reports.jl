# =============================================================================
# Verification-suite report mode.
# =============================================================================
#
# Runs verification cases through their report resolution ladders and writes
# the publishable artifact set (overlay/error/convergence figures, CSVs,
# generated guide chapters) to the output directory — the counterpart to the
# fast assert mode in test/verification/runcases.jl, built from the SAME case
# definitions.
#
# Usage:
#   julia --project=scripts/verification -e 'using Pkg; Pkg.instantiate()'   # once
#   julia --project=scripts/verification scripts/verification/run_reports.jl [options]
#
# Options (all --key=value; unknown keys are an error):
#   --cases=all              comma-separated case ids (members like V-HC-02a or
#                            whole families like V-HC-02), or "all"
#   --outdir=docs/guide/verification
#                            artifact root; chapters at its top level, figures
#                            and CSVs under assets/<case id>/. The default is
#                            the COMMITTED guide location — review the diff.
#   --chapters=true          also (re)generate <outdir>/*.md guide chapters
#
# Cases with `report = (;)` are skipped (listed in the index dashboard).
# Nothing here asserts: coarse ladder rungs may exceed the pinned CI
# tolerances by design; tables record within-tolerance yes/no instead.

const ROOT = dirname(dirname(@__DIR__))

include(joinpath(ROOT, "test", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "render.jl"))
include(joinpath(@__DIR__, "report_harness.jl"))
include(joinpath(@__DIR__, "chapters.jl"))

function parse_cli(args)
    opts = Dict{String,String}(
        "cases" => "all",
        "outdir" => joinpath(ROOT, "docs", "guide", "verification"),
        "chapters" => "true",
    )
    for a in args
        startswith(a, "--") || error("unrecognized argument $a (use --key=value)")
        kv = split(a[3:end], '='; limit = 2)
        length(kv) == 2 || error("missing value in $a (use --key=value)")
        haskey(opts, kv[1]) || error("unknown option --$(kv[1])")
        opts[kv[1]] = kv[2]
    end
    return opts
end

function select_cases(cases, spec::AbstractString)
    spec == "all" && return cases
    wanted = strip.(split(spec, ','))
    selected = [c for c in cases if c.id in wanted || base_case_id(c.id) in wanted]
    unmatched =
        [w for w in wanted if !any(c -> c.id == w || base_case_id(c.id) == w, cases)]
    isempty(unmatched) || error("unknown case id(s): $(join(unmatched, ", "))")
    return selected
end

opts = parse_cli(ARGS)
outdir = abspath(opts["outdir"])
mkpath(outdir)
set_theme!(report_theme())

all_cases = load_verification_cases()
selected = select_cases(all_cases, opts["cases"])
@info "report mode" n_selected = length(selected) n_total = length(all_cases) outdir

results = CaseReportResult[]
for case in selected
    if isempty(pairs(case.report))
        @info "skipping (empty report metadata)" case.id
        continue
    end
    @info "running case report" case.id
    # invokelatest: case closures were defined during load_verification_cases
    # (Julia 1.12 world-age semantics — see casedefs.jl).
    r = Base.invokelatest(run_case_report, case, outdir)
    if r !== nothing
        push!(results, r)
        # Persist completed reports even if a subsequent case fails. Passing
        # accumulated results retains every completed member of a case family.
        if opts["chapters"] == "true"
            Base.invokelatest(write_chapters, results, all_cases, outdir)
        end
    end
end

@info "report mode done" n_reported = length(results) outdir
