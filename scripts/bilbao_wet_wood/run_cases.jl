# Solve runs under variants and write their series to results/<variant>/<run>.csv.
#   julia --project=scripts/verification scripts/bilbao_wet_wood/run_cases.jl VARIANT [RUN ...]
# VARIANT is family/law[/tag] (variants.jl) or `raw:key=value,key=value` for a
# one-off; the runs default to all moist runs, plus the dry run under `dry`.
include(joinpath(@__DIR__, "BilbaoWetWood.jl"))
using .BilbaoWetWood
using Printf
using TOML
using Pyrolysis: SERIES, PARALLEL
include(joinpath(@__DIR__, "variants.jl"))

variant = ARGS[1]
runs = length(ARGS) > 1 ? ARGS[2:end] : sort([k for k in keys(RUNS) if RUNS[k].Y0 > 0])
kw = variant_kwargs(variant)
tag = replace(variant, '/' => '_')
for name in runs
    r = RUNS[name]
    t0 = time()
    s = try
        simulate(r; kw...)
    catch e
        @warn "$variant $name failed" exception = e
        continue
    end
    path = write_series(series_path(tag, name), s)
    @printf(
        "%-28s %-14s %s  %.0f s  mass %.1e energy %.1e -> %s\n",
        variant,
        name,
        s.retcode,
        time() - t0,
        s.audit.mass.relative,
        s.audit.energy.relative,
        path
    )
end
