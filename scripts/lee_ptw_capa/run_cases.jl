# Solve runs under a variant and write their series to results/<variant>/<run>.csv.
#   julia --project=scripts/verification scripts/lee_ptw_capa/run_cases.jl VARIANT [RUN ...]
# VARIANT is law[/tag] (variants.jl); the runs default to both, "63" and "29".
include(joinpath(@__DIR__, "LeePtwCapa.jl"))
using .LeePtwCapa
using Printf
include(joinpath(@__DIR__, "variants.jl"))

variant = ARGS[1]
runs = length(ARGS) > 1 ? ARGS[2:end] : ["63", "29"]
kw = variant_kwargs(variant)
tag = replace(variant, '/' => '_')
for name in runs
    t0 = time()
    s = try
        simulate(RUNS[name]; kw...)
    catch e
        @warn "$variant $name failed" exception = e
        continue
    end
    path = write_series(series_path(tag, name), s)
    @printf(
        "%-24s %-3s %s  %.0f s  mass %.1e energy %.1e -> %s\n",
        variant,
        name,
        s.retcode,
        time() - t0,
        s.audit.mass.relative,
        s.audit.energy.relative,
        path
    )
    flush(stdout)
end
