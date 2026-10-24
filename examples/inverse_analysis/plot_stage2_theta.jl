"""
Re-plot White Pine fit at a *given* Stage-2 θ — no optimization.

Forward-solves every fit/val scenario at the supplied 6-component conductivity
vector and writes the experiment-vs-simulation traces, reusing the exact
loaders and `plot_fits` from `inverse_white_pine.jl`.

USAGE
-----
    julia -t auto --project=examples examples/inverse_analysis/plot_stage2_theta.jl

Edit `THETA` below (or pass --theta=v1,v2,...,v6) to change the plotted point.
Stage-2 ordering: [k_virgin, k_int1, k_int2, k_int3, k_int4, k_char].
Optionally restrict the plotted fluxes with --fit=... / --val=... (same labels
as the driver); defaults to the registry roles (fit 20/30/40/50/60, val 70).
"""

# Reuse the driver's loaders + plot_fits. main() is guarded by a PROGRAM_FILE
# check, so including this does NOT kick off an optimization.
include("inverse_white_pine.jl")

# PYROLYSIS_SMOKE=1 (`SMOKE` and `smoke_dataset` come from inverse_white_pine.jl):
# every scenario over its first 20 s, the figure written beside the full run's
load_dataset(label) =
    SMOKE ? smoke_dataset(load_registry_dataset(label)) : load_registry_dataset(label)

const STAGE = 2

# Default to UMD-only fluxes (30/50 kW/m² fit, 70 kW/m² validation); the FSRI
# 20/40/60 cases are dropped. CLI --fit=/--val= still override these.
const DEFAULT_FIT_LABELS = ["30kW", "50kW"]
const DEFAULT_VAL_LABELS = ["70kW"]

# Fitted Stage-2 conductivities to plot.
const THETA_DEFAULT = [
    0.09902,  # k_virgin ± 0.00016
    0.1031,   # k_int1   ± 0.0025
    0.05003,  # k_int2   ± 0.0029
    0.1014,   # k_int3   ± 0.0012
    0.3337,   # k_int4   ± 0.0034
    0.1322,   # k_char   ± 0.0037
]

function _parse_theta(args)
    for a in args
        startswith(a, "--theta=") || continue
        return [parse(Float64, strip(s)) for s in split(split(a, "="; limit = 2)[2], ',')]
    end
    return THETA_DEFAULT
end

function plot_main()
    mkpath(RESULTS_DIR)
    θ = _parse_theta(ARGS)
    length(θ) == 6 || error("Stage-2 θ needs 6 values, got $(length(θ)): $θ")

    fit_labels = OPTS.fit_labels === nothing ? DEFAULT_FIT_LABELS : OPTS.fit_labels
    val_labels = OPTS.val_labels === nothing ? DEFAULT_VAL_LABELS : OPTS.val_labels
    @info "Re-plotting at supplied Stage-2 θ" θ fit = fit_labels val = val_labels

    fit_dataset = [load_dataset(l) for l in fit_labels]
    val_dataset = Any[]
    for l in val_labels
        try
            push!(val_dataset, load_dataset(l))
        catch err
            @warn "Validation data not loadable, skipping" label = l exception = err
        end
    end

    outpath = joinpath(
        RESULTS_DIR,
        SMOKE ? "fit_plots_stage2_manual_smoke.png" : "fit_plots_stage2_manual.png",
    )
    plot_fits(θ, fit_dataset, val_dataset, STAGE; outpath = outpath)
    @info "Done." outpath
end

plot_main()
