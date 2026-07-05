# JET.jl static type-inference checks (the Julia analogue of mypy).
#
# Run standalone:  julia --project=test/quality test/quality/jet.jl
#
# Policy (gating, owned-only): JET's *whole-package* analysis is dominated by
# dependency-internal noise we can't fix here — ProgressMeter version-skew
# FieldErrors, SciMLBase closure introspection (`numargs`/`isinplace` probing our
# solver closures for a field `r`), and Distributed/`asyncmap` `asyncrun` undefs.
# `target_modules=(Pyrolysis,)` restricts reports to code we own; we gate that set
# to zero. To inspect the full (noisy) report, drop `target_modules`.
using JET
using Pyrolysis
using Test

@testset "JET — type analysis (Pyrolysis-owned)" begin
    result =
        report_package(Pyrolysis; target_modules = (Pyrolysis,), toplevel_logger = nothing)
    reports = JET.get_reports(result)

    if isempty(reports)
        @info "JET: no Pyrolysis-owned reports."
    else
        @warn "JET: $(length(reports)) Pyrolysis-owned report(s) detected — see below"
        show(stdout, result)
        println(stdout)
    end

    @test isempty(reports)
end
