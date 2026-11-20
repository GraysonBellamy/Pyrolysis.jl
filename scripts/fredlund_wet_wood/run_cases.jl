# Solve cases under a variant and write their series to results/<variant>/<case>.csv.
#   julia --project=scripts/verification scripts/fredlund_wet_wood/run_cases.jl VARIANT [CASE ...]
# VARIANT is law[/tag] (variants.jl); the cases default to all three.
include(joinpath(@__DIR__, "FredlundWetWood.jl"))
using .FredlundWetWood
using .FredlundWetWood: PHASE_CHANGE_LAWS
using Printf
using TOML
include(joinpath(@__DIR__, "variants.jl"))

variant = ARGS[1]
cases = length(ARGS) > 1 ? ARGS[2:end] : collect(CASES)
tag = replace(variant, '/' => '_')
for case in cases
    t0 = time()
    s = try
        simulate(case; case_kwargs(variant, case)...)
    catch e
        @warn "$variant $case failed" exception = e
        continue
    end
    path = write_series(series_path(tag, case), s)
    @printf(
        "%-22s %-12s %s  %.0f s  mass %.1e energy %.1e -> %s\n",
        variant,
        case,
        s.retcode,
        time() - t0,
        s.audit.mass.relative,
        s.audit.energy.relative,
        path
    )
    flush(stdout)
end
