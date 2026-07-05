"""
White Pine — full diagnostic comparison plots at a given Stage-2 θ.

For each flux, forward-solves ONCE (reuse the optimizer's ApproxSolve Jacobian
with --approx-solve) and plots a column of diagnostics:

    row 1  back-surface T   [°C]   sim vs experiment (UMD CAPA, the fit target)
    row 2  front/surface T  [°C]   sim only (no UMD front-surface measurement;
                                   FSRI top-surface overlaid when present)
    row 3  mass             [g]    sim (total_mass·A0) vs experiment (UMD CAPA)
    row 4  MLR (gauge)      [g/m²/s] sim vs experiment (MaCFP detailed, ±SE band)
    row 5  thickness        [mm]   sim curve + measured post-test L_f marker
                                   (UMD CAPA replicate mean, wp_properties.jl)

Experimental curves are drawn only where data exists; otherwise the panel shows
the simulation alone. No re-fitting — θ is fixed.

USAGE
-----
  GKSwstype=100 julia -t 4 --project=examples \\
      examples/inverse_analysis/plot_diagnostics.jl \\
      --fit=30kW,50kW --val=70kW --approx-solve

Edit THETA_DEFAULT (or pass --theta=v1,...,v6). Stage-2 ordering:
[k_virgin, k_int1, k_int2, k_int3, k_int4, k_char].
"""

include("inverse_white_pine.jl")

const STAGE = 2

const THETA_DEFAULT = [
    0.09902,  # k_virgin ± 0.00016
    0.1031,   # k_int1   ± 0.0025
    0.05003,  # k_int2   ± 0.0029
    0.1014,   # k_int3   ± 0.0012
    0.3337,   # k_int4   ± 0.0034
    0.1322,   # k_char   ± 0.0037
]

# Default to UMD-only fluxes (30/50 kW/m² fit, 70 kW/m² validation); the FSRI
# 20/40/60 cases are dropped. CLI --fit=/--val= still override these.
const DEFAULT_FIT_LABELS = ["30kW", "50kW"]
const DEFAULT_VAL_LABELS = ["70kW"]

function _parse_theta(args)
    for a in args
        startswith(a, "--theta=") || continue
        return [parse(Float64, strip(s)) for s in split(split(a, "="; limit = 2)[2], ',')]
    end
    return THETA_DEFAULT
end

# ----------------------------------------------------------------------------
# Generic replicate-averaged experimental column loader.
# Tolerant variant of `load_back_temp_replicates`: rounds time to the nearest
# second and averages every replicate value landing in that bin. Returns
# `nothing` when no file matches or the column is absent everywhere.
# ----------------------------------------------------------------------------
function _load_col_reps(pattern::Regex, colpred; data_dir = DATA_DIR, unit = identity)
    files = sort([joinpath(data_dir, f) for f in readdir(data_dir) if occursin(pattern, f)])
    isempty(files) && return nothing
    acc = Dict{Int,Vector{Float64}}()
    nfiles = 0
    for f in files
        headers, data = _read_csv_with_header(f)
        t_col = findfirst(h -> startswith(h, "Time") || startswith(h, "Timestamp"), headers)
        y_col = findfirst(colpred, headers)
        (t_col === nothing || y_col === nothing) && continue
        nfiles += 1
        t_raw = data[:, t_col]
        y_raw = data[:, y_col]
        for i in eachindex(t_raw)
            (isfinite(t_raw[i]) && isfinite(y_raw[i])) || continue
            τ = round(Int, t_raw[i])
            push!(get!(acc, τ, Float64[]), unit(Float64(y_raw[i])))
        end
    end
    (nfiles == 0 || isempty(acc)) && return nothing
    ts = sort(collect(keys(acc)))
    return (
        t = Float64.(ts),
        y = [mean(acc[τ]) for τ in ts],
        y_std = [length(acc[τ]) > 1 ? std(acc[τ]) : 0.0 for τ in ts],
        n = nfiles,
    )
end

# Experimental mass [g] — UMD CAPA replicate files (same source as the back-T fit).
_exp_mass(flux) = _load_col_reps(
    Regex("^$(REPLICATE_PATTERN_PREFIX)_$(flux)_R\\d+\\.csv\$"),
    h -> startswith(h, "Mass"),
)

# Experimental gauge MLR + its SE [g/m²/s] — MaCFP White Pine detailed files.
# Value and SE are averaged on ONE shared time grid (only rows where BOTH are
# finite), so the ribbon is guaranteed length-aligned with the curve.
function _exp_mlr(flux; data_dir = DATA_DIR)
    pattern =
        Regex("^$(BC_DETAILED_PATTERN_PREFIX)_$(flux)_.*$(BC_DETAILED_PATTERN_SUFFIX)\$")
    files = sort([joinpath(data_dir, f) for f in readdir(data_dir) if occursin(pattern, f)])
    isempty(files) && return nothing
    acc_y = Dict{Int,Vector{Float64}}()
    acc_s = Dict{Int,Vector{Float64}}()
    for f in files
        headers, data = _read_csv_with_header(f)
        t_col = findfirst(h -> startswith(h, "Time"), headers)
        y_col = findfirst(h -> startswith(h, "MLR ("), headers)
        s_col = findfirst(h -> startswith(h, "MLR SE"), headers)
        (t_col === nothing || y_col === nothing || s_col === nothing) && continue
        for i in axes(data, 1)
            t, y, s = data[i, t_col], data[i, y_col], data[i, s_col]
            (isfinite(t) && isfinite(y) && isfinite(s)) || continue
            τ = round(Int, t)
            push!(get!(acc_y, τ, Float64[]), 1000 * Float64(y))
            push!(get!(acc_s, τ, Float64[]), 1000 * Float64(s))
        end
    end
    isempty(acc_y) && return nothing
    ts = sort(collect(keys(acc_y)))
    return (
        t = Float64.(ts),
        y = [mean(acc_y[τ]) for τ in ts],
        se = [mean(acc_s[τ]) for τ in ts],
        n = length(files),
    )
end

# Measured post-test final thickness L_f [mm] (UMD CAPA N₂ replicate means).
# Values are the canonical WP_LF_*_MM constants from wp_properties.jl; the
# tuple carries the replicate count for the legend.
const WP_LF_MEAS =
    Dict("30kW" => (WP_LF_30_MM, 3), "50kW" => (WP_LF_50_MM, 4), "70kW" => (WP_LF_70_MM, 3))

# ----------------------------------------------------------------------------
# One diagnostic column (5 stacked panels) for a single flux.
# ----------------------------------------------------------------------------
function flux_column(label_prefix, ds, θ; jacobian_strategy, panels)
    scenario, t_exp, T_back_exp, T_back_std, _ = ds
    flux = scenario.label
    q_kw = round(scenario.q_ext / 1000, digits = 1)
    head = "$(label_prefix): $(flux) ($(q_kw) kW/m²)"
    sol = build_and_solve(θ, scenario; stage = STAGE, jacobian_strategy = jacobian_strategy)
    sol === nothing && @warn "Forward solve returned nothing" flux

    # --- row 1: back-surface temperature (± replicate SD band) ---
    p_back = plot(
        t_exp,
        T_back_exp .- 273.15;
        ribbon = T_back_std,
        label = "experiment (±SD)",
        color = :black,
        lw = 2,
        fillalpha = 0.3,
        xlabel = "Time [s]",
        ylabel = "T_back [°C]",
        title = "$head — back T",
    )
    sol !== nothing && plot!(
        p_back,
        sol.t,
        sol.extras.back_temperature .- 273.15;
        label = "simulation",
        color = :red,
        lw = 2,
        linestyle = :dash,
    )
    panels[:back] = p_back

    # --- row 2: front / surface temperature (sim only for UMD) ---
    p_front = plot(xlabel = "Time [s]", ylabel = "T_front [°C]", title = "$head — front T")
    if startswith(flux, "FSRI") || startswith(scenario.label, "FSRI_")
        top = try
            load_fsri_top_surface_temp(replace(flux, "FSRI_" => ""); verbose = false)
        catch
            nothing
        end
        top !== nothing && scatter!(
            p_front,
            top.t,
            top.T .- 273.15;
            label = "experiment",
            color = :black,
            ms = 3,
            msw = 0,
        )
    end
    sol !== nothing && plot!(
        p_front,
        sol.t,
        sol.surface_T .- 273.15;
        label = "simulation",
        color = :red,
        lw = 2,
        linestyle = :dash,
    )
    panels[:front] = p_front

    # --- row 3: mass [g] ---
    p_mass = plot(xlabel = "Time [s]", ylabel = "Mass [g]", title = "$head — mass")
    em = _exp_mass(flux)
    em !== nothing &&
        plot!(p_mass, em.t, em.y; label = "experiment", color = :black, lw = 2)
    if sol !== nothing
        m_sim_g = sol.extras.total_mass .* scenario.A_0 .* 1000
        plot!(
            p_mass,
            sol.t,
            m_sim_g;
            label = "simulation",
            color = :red,
            lw = 2,
            linestyle = :dash,
        )
        @info "  mass m0 [g]" flux sim = round(m_sim_g[1], digits = 3) exp =
            (em === nothing ? NaN : round(em.y[1], digits = 3))
    end
    panels[:mass] = p_mass

    # --- row 4: gauge MLR [g/m²/s], x pinned to [0, t_final] ---
    p_mlr = plot(
        xlabel = "Time [s]",
        ylabel = "MLR [g/m²/s]",
        title = "$head — MLR",
        xlims = (0.0, scenario.t_final),
        ylims = (0.0, Inf),
    )
    emlr = _exp_mlr(flux)
    if emlr !== nothing
        plot!(
            p_mlr,
            emlr.t,
            emlr.y;
            ribbon = emlr.se,
            label = "experiment (±SE)",
            color = :black,
            lw = 2,
            fillalpha = 0.35,
        )
        @info "  MLR exp" flux peak_g_m2_s = round(maximum(emlr.y), digits = 2) max_SE =
            round(maximum(emlr.se), digits = 3)
    end
    sol !== nothing && plot!(
        p_mlr,
        sol.t,
        sol.extras.MLR_gauge .* 1000;
        label = "simulation",
        color = :red,
        lw = 2,
        linestyle = :dash,
    )
    panels[:mlr] = p_mlr

    # --- row 5: thickness [mm] — sim curve + measured post-test L_f marker ---
    p_thk =
        plot(xlabel = "Time [s]", ylabel = "Thickness [mm]", title = "$head — thickness")
    sol !== nothing && plot!(
        p_thk,
        sol.t,
        sol.extras.thickness_history .* 1000;
        label = "simulation",
        color = :red,
        lw = 2,
        linestyle = :dash,
    )
    if haskey(WP_LF_MEAS, flux)
        Lf_mm, n_rep = WP_LF_MEAS[flux]
        scatter!(
            p_thk,
            [scenario.t_final],
            [Lf_mm];
            label = "measured L_f (n=$n_rep)",
            color = :black,
            marker = :diamond,
            ms = 7,
            msw = 0,
        )
        sol !== nothing && @info "  thickness L_f [mm]" flux sim =
            round(sol.extras.thickness_history[end] * 1000, digits = 2) meas = Lf_mm
    end
    panels[:thk] = p_thk

    return nothing
end

function diag_main()
    θ = _parse_theta(ARGS)
    length(θ) == 6 || error("Stage-2 θ needs 6 values, got $(length(θ)): $θ")
    # ApproxSolve (~14% nnz, the :ale_reduced profile) is the default Jacobian
    # strategy here; the full DirectSolve path is not used.
    jac_strat = ApproxSolve()

    fit_labels = OPTS.fit_labels === nothing ? DEFAULT_FIT_LABELS : OPTS.fit_labels
    val_labels = OPTS.val_labels === nothing ? DEFAULT_VAL_LABELS : OPTS.val_labels
    @info "Diagnostic plots at Stage-2 θ" θ fit = fit_labels val = val_labels jacobian = "ApproxSolve"

    datasets = Tuple{String,Any}[]
    for l in fit_labels
        push!(datasets, ("Fit", load_registry_dataset(l)))
    end
    for l in val_labels
        try
            push!(datasets, ("Validation", load_registry_dataset(l)))
        catch err
            @warn "Validation data not loadable, skipping" label = l exception = err
        end
    end

    # Row order top→bottom; columns are fluxes left→right.
    rows = [:back, :front, :mass, :mlr, :thk]
    ncols = length(datasets)
    column_panels = Vector{Dict{Symbol,Any}}(undef, ncols)
    for (j, (prefix, ds)) in enumerate(datasets)
        panels = Dict{Symbol,Any}()
        flux_column(prefix, ds, θ; jacobian_strategy = jac_strat, panels = panels)
        column_panels[j] = panels
    end

    # Assemble row-major: (row r, col j) -> plots[(r-1)*ncols + j]
    plots = Plots.Plot[]
    for r in rows, j = 1:ncols
        push!(plots, column_panels[j][r])
    end
    fig = plot(
        plots...;
        layout = (length(rows), ncols),
        size = (430 * ncols, 260 * length(rows)),
        legend = :best,
        left_margin = 5Plots.mm,
    )
    outpath = joinpath(RESULTS_DIR, "diagnostics_stage2_manual.png")
    savefig(fig, outpath)
    @info "Diagnostic plots written" outpath
end

diag_main()
