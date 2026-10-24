"""
White Pine Inverse Analysis — Driver

Recovers component thermal conductivities by fitting CAPA back-surface
temperature time-series at multiple cone heat-flux levels.

USAGE
-----
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl                                       # default: stage 1, fit 20/30/40/50/60, val 70
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --stage=3
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --stage=1 --no-grid
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --prior-weight=0                      # data only, no HFM/LFA prior
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --prior-weight=1.0                    # max prior weight
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --bootstrap=200                       # opt-in bootstrap covariance
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --fit=30kW,50kW                       # UMD-only fit (drops all FSRI)
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --fit=20kW,40kW,60kW --val=          # FSRI-only fit, no validation
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --fit=20kW,30kW,40kW,50kW,60kW,70kW   # everything in fit (no val)
	julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --approx-solve                        # ApproxSolve Jacobian (faster, ~14% nnz)

Spread the heavy work across cores via `-t auto` (or `-t N`); the loss runs
all fit scenarios in parallel via `Threads.@spawn`.

`--fit=<csv>`  Comma-separated flux labels for the fit dataset.
                  Default: every flux marked `:fit` in `DATASET_REGISTRY`
                  (currently 20/30/40/50/60 kW). Pass an empty string
                  (`--fit=`) to error out — fit needs ≥1 dataset.
`--val=<csv>`  Comma-separated flux labels for the validation dataset.
                  Default: every flux marked `:val` (currently 70 kW). Pass
                  `--val=` to skip validation entirely.

`--prior-weight=W` (W ∈ [0, 1], default 0.5):
	0    → pure data fit (ignore HFM/LFA k(T) prior)
	0.5  → prior weighted equal to one flux scenario
	1.0  → prior weighted equal to the combined data term

Available datasets in examples/inverse_analysis/data/ (registry below):
	FSRI_Wood_Gasification_20kW_R*.csv (FSRI gasification, This-Study BCs, 3% MC)
	FSRI_Wood_CAPA_N2_40kW_R*.csv      (FSRI CAPA, This-Study BCs, 3% MC)
	FSRI_Wood_CAPA_N2_60kW_R*.csv      (FSRI CAPA, This-Study BCs, 3% MC)
	UMD_Wood_CAPA_N2_30kW_R*.csv       (UMD CAPA, MaCFP detailed BCs, dry)
	UMD_Wood_CAPA_N2_50kW_R*.csv       (UMD CAPA, MaCFP detailed BCs, dry)
	UMD_Wood_CAPA_N2_70kW_R*.csv       (UMD CAPA, MaCFP detailed BCs, dry)
	bcs_<flux>kW.csv                   (optional override on UMD BCs)

Outputs in examples/inverse_analysis/results/run_<YYYYMMDD_HHMMSS>_s<stage>_pw<weight>/:
	trace.csv         per-evaluation log (data + prior contributions)
	grid_scan.csv     residual surface (stage 1 only, if grid enabled)
	best.json         optimum θ + standard errors + correlation matrix
	fit_plots.png     experimental vs simulated traces

Each invocation writes to its own timestamped subdirectory so back-to-back
runs (e.g. comparing prior weights or stages) never overwrite each other.
"""

using Pkg
using Printf
using Statistics
using LinearAlgebra
using Dates
using Random

using Optim
using FiniteDiff
using Plots
using QuasiMonteCarlo

include("wp_objective.jl")

# ============================================================================
# CLI parsing
# ============================================================================

function parse_args(args)
    stage = 1
    do_grid = true
    n_starts = 3
    max_iter = 100   # NelderMead converges in fewer iterations than LBFGS+FD
    bootstrap = 0    # 0 = use Hessian; >0 = bootstrap with this many resamples
    prior_weight = 0.5   # weight on HFM/LFA k(T) prior; clamped to [0, 1]
    fit_labels = nothing   # nothing = use registry defaults
    val_labels = nothing
    use_approx_solve = false
    for a in args
        if startswith(a, "--stage=")
            stage = parse(Int, split(a, "=")[2])
        elseif a == "--no-grid"
            do_grid = false
        elseif startswith(a, "--starts=")
            n_starts = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--max-iter=")
            max_iter = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--bootstrap=")
            bootstrap = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--prior-weight=")
            prior_weight = parse(Float64, split(a, "=")[2])
            (0.0 <= prior_weight <= 1.0) ||
                error("--prior-weight must be in [0, 1], got $prior_weight")
        elseif startswith(a, "--fit=")
            rhs = split(a, "="; limit = 2)[2]
            fit_labels =
                isempty(rhs) ? String[] :
                String[strip(s) for s in split(rhs, ',') if !isempty(strip(s))]
        elseif startswith(a, "--val=")
            rhs = split(a, "="; limit = 2)[2]
            val_labels =
                isempty(rhs) ? String[] :
                String[strip(s) for s in split(rhs, ',') if !isempty(strip(s))]
        elseif a == "--approx-solve"
            use_approx_solve = true
        elseif a in ("-h", "--help")
            println(@doc("inverse_white_pine.jl"))
            exit(0)
        end
    end
    (;
        stage,
        do_grid,
        n_starts,
        max_iter,
        bootstrap,
        prior_weight,
        fit_labels,
        val_labels,
        use_approx_solve,
    )
end

const OPTS = parse_args(ARGS)

"Pick labels for a role from the registry, honoring CLI overrides."
function select_labels(role::Symbol, override::Union{Nothing,Vector{String}})
    override !== nothing && return override
    return String[p.first for p in DATASET_REGISTRY if p.second.default_role == role]
end

# ============================================================================
# Stage configuration
# ============================================================================

function stage_config(stage::Int)
    if stage == 1
        (;
            θ_init = [WP_K_VIRGIN_INIT, WP_K_CHAR_INIT],
            bounds = [WP_K_VIRGIN_BOUNDS, WP_K_CHAR_BOUNDS],
            names = ["k_virgin", "k_char"],
        )
    elseif stage == 2
        bounds = [
            WP_K_VIRGIN_BOUNDS,
            WP_K_INT_BOUNDS,
            WP_K_INT_BOUNDS,
            WP_K_INT_BOUNDS,
            WP_K_INT_BOUNDS,
            WP_K_CHAR_BOUNDS,
        ]
        (;
            θ_init = wp_stage2_init(),
            bounds = bounds,
            names = ["k_virgin", "k_int1", "k_int2", "k_int3", "k_int4", "k_char"],
        )
    elseif stage == 3
        (;
            θ_init = WP_STAGE3_INIT,
            bounds = WP_STAGE3_BOUNDS,
            names = ["k_v_lo", "dk_v", "k_c_lo", "dk_c"],
        )
    else
        error("unknown stage: $stage (must be 1, 2, or 3)")
    end
end

# ============================================================================
# Multi-start initialization
# ============================================================================

"""
	multi_start_initial_guesses(θ_init, bounds, n_starts; grid_seeds=nothing)

Build the multi-start list. First entry is always the nominal `θ_init`. If
`grid_seeds` is provided (a vector of θ vectors from grid-scan minima), those
are inserted next. The remaining slots are filled with Sobol-sequence points
scaled to `bounds`. Sobol gives much better box coverage than random LCG
draws for small `n_starts`.
"""
function multi_start_initial_guesses(
    θ_init,
    bounds,
    n_starts::Int;
    grid_seeds::Union{Nothing,AbstractVector} = nothing,
)

    starts = [collect(Float64.(θ_init))]
    if grid_seeds !== nothing
        for s in grid_seeds
            length(starts) < n_starts || break
            push!(starts, collect(Float64.(s)))
        end
    end
    n_remaining = n_starts - length(starts)
    if n_remaining > 0
        lo = Float64[b[1] for b in bounds]
        hi = Float64[b[2] for b in bounds]
        # QuasiMonteCarlo returns a (p × n) matrix
        pts = QuasiMonteCarlo.sample(n_remaining, lo, hi, SobolSample())
        for j = 1:n_remaining
            push!(starts, Vector{Float64}(pts[:, j]))
        end
    end
    return starts
end

"""
Read the Stage 1 grid scan output and return up to `k` lowest-loss θ vectors
to use as multi-start seeds. Returns an empty vector if the file is missing.
"""
function grid_scan_seeds(; k::Int = 2)
    path = joinpath(RUN_DIR[], "grid_scan.csv")
    isfile(path) || return Vector{Vector{Float64}}()
    rows = Tuple{Float64,Float64,Float64}[]
    open(path, "r") do io
        readline(io)  # header
        for line in eachline(io)
            parts = split(strip(line), ',')
            length(parts) == 3 || continue
            push!(
                rows,
                (
                    parse(Float64, parts[1]),
                    parse(Float64, parts[2]),
                    parse(Float64, parts[3]),
                ),
            )
        end
    end
    sort!(rows, by = r -> r[3])
    return [Float64[r[1], r[2]] for r in first(rows, min(k, length(rows)))]
end

# ============================================================================
# Optimization
# ============================================================================

function run_optimization(
    loss_fn,
    θ_init,
    bounds;
    max_iter::Int = 100,
    n_starts::Int = 3,
    verbose::Bool = true,
    grid_seeds::Union{Nothing,AbstractVector} = nothing,
)

    starts = multi_start_initial_guesses(θ_init, bounds, n_starts; grid_seeds = grid_seeds)
    @assert length(bounds) == length(θ_init) "bounds length must match parameter count"
    lower = [b[1] for b in bounds]
    upper = [b[2] for b in bounds]
    results = []
    for (i, θ0) in enumerate(starts)
        verbose && @printf(
            "\n--- Multi-start %d/%d, θ0 = %s ---\n",
            i,
            length(starts),
            _fmt_vec(θ0)
        )
        t0 = time()
        # NelderMead is derivative-free → no FD gradient evals, robust to the
        # loss-surface noise from solver tolerances + interpolation kinks.
        # Fminbox handles the box constraints; loss has no in-band penalty
        # anymore (only the FAILURE_PENALTY for solver failures).
        opts =
            Optim.Options(iterations = max_iter, show_trace = verbose, store_trace = false)
        res = optimize(loss_fn, lower, upper, θ0, Fminbox(NelderMead()), opts)
        push!(
            results,
            (
                θ = Optim.minimizer(res),
                loss = Optim.minimum(res),
                converged = Optim.converged(res),
                wall = time() - t0,
                start_idx = i,
            ),
        )
        verbose && @printf(
            "    finished in %.1fs, loss = %.6e, converged = %s\n",
            results[end].wall,
            results[end].loss,
            results[end].converged
        )
    end
    # Sort by loss ascending; first is best
    sort!(results, by = r -> r.loss)
    return results
end

# ============================================================================
# Parameter covariance — Hessian-based (default) and bootstrap (--bootstrap=N)
# ============================================================================
#
# Derivation (Hessian-based, σ-normalized weighted MSE loss):
#
#   Loss:         L(θ) = (1/N) · Σ_i (r_i(θ) / σ_i)²       with N total points
#   Gradient:     ∇L  = (2/N) · Jᵀ W r                      W = diag(1/σ_i²)
#   Hessian:     H_L ≈ (2/N) · Jᵀ W J                       (Gauss-Newton)
#   ⇒ JᵀWJ      = (N/2) · H_L
#
#   Standard weighted-NLS covariance:
#       Σ_θ̂ = (JᵀWJ)⁻¹ · ŝ²,   ŝ² = S_w / (N−p),  S_w = Σ(r_i/σ_i)² = N·L(θ̂)
#
#   ⇒  Σ_θ̂ = (2/N) · H_L⁻¹ · N·L(θ̂)/(N−p)
#          = (2 · L(θ̂) / (N−p)) · H_L⁻¹
#
# This assumes the σ_i correctly capture the noise. If they're underestimated,
# the resulting σ_θ is correspondingly underestimated; the χ²_red of the fit
# (≈ L(θ̂)) gives a sanity check (≈ 1 when σ is honest).

function estimate_covariance(
    loss_fn,
    θ_hat::AbstractVector,
    dataset;
    stage::Int,
    t_skip::Real = 1.0,
    fd_step::AbstractVector = fill(1e-3, length(θ_hat)),
    loss_at_theta::Union{Nothing,Real} = nothing,
)

    N = 0
    for (_, t_exp, _, _, _) in dataset
        N += count(>(t_skip), t_exp)
    end
    p = length(θ_hat)
    dof = max(N - p, 1)

    H = try
        FiniteDiff.finite_difference_hessian(
            loss_fn,
            θ_hat,
            Val(:hcentral);
            absstep = minimum(fd_step),
        )
    catch err
        @warn "Hessian FD failed" exception=err
        return nothing
    end

    L_hat = loss_at_theta === nothing ? loss_fn(θ_hat) : float(loss_at_theta)

    local Σ
    try
        Σ = (2.0 * L_hat / dof) * inv(H)
    catch err
        @warn "Hessian inversion failed (likely ill-conditioned at the bound)" exception=err
        return nothing
    end
    se = sqrt.(abs.(diag(Σ)))
    corr = Σ ./ (se * se')
    return (Σ = Σ, se = se, corr = corr, dof = dof, loss = L_hat, method = :hessian)
end

"""
	bootstrap_covariance(loss_fn, θ_hat, dataset; n_boot=200, stage, t_skip=1.0,
	                     refit_max_iter=50)

Residual bootstrap: at each draw, resample the per-time residuals (with
replacement), add them back to the simulated trajectory, and refit. Slow but
robust to loss-surface noise. Re-uses `loss_fn` so the same weighting and
priors apply.

Returns the same NamedTuple shape as `estimate_covariance`. `θ_hat` becomes
the median of the bootstrap distribution; `corr` is computed from the
bootstrap covariance.
"""
function bootstrap_covariance(
    loss_fn,
    θ_hat::AbstractVector,
    dataset,
    bounds;
    n_boot::Int = 200,
    stage::Int,
    t_skip::Real = 1.0,
    refit_max_iter::Int = 50,
    seed::Int = 42,
    verbose::Bool = true,
    jacobian_strategy::Union{Nothing,LinearSolveStrategy} = nothing,
)

    p = length(θ_hat)
    lower = [b[1] for b in bounds]
    upper = [b[2] for b in bounds]
    rng = MersenneTwister(seed)

    # Simulate the "best fit" trajectory once per scenario. The bootstrap
    # regenerates synthetic targets by adding resampled residuals to this.
    base_traj = []
    base_residuals = []
    for (scenario, t_exp, T_back_exp, T_back_std, n_reps) in dataset
        sol = build_and_solve(
            θ_hat,
            scenario;
            stage = stage,
            jacobian_strategy = jacobian_strategy,
        )
        sol === nothing && error("base solve failed at θ_hat — cannot bootstrap")
        mask = t_exp .> t_skip
        t_q = t_exp[mask]
        T_sim = interpolate_at(t_q, sol.t, sol.extras.back_temperature)
        r = T_sim .- T_back_exp[mask]
        push!(
            base_traj,
            (scenario, t_exp, T_back_exp, T_back_std, n_reps, t_q, T_sim, r, mask),
        )
        push!(base_residuals, r)
    end

    θ_samples = Matrix{Float64}(undef, p, n_boot)
    n_failed = 0
    for b = 1:n_boot
        # Build resampled "experimental" dataset
        boot_ds = []
        for (scenario, t_exp, T_back_exp, T_back_std, n_reps, _, T_sim, r, mask) in
            base_traj

            r_boot = r[rand(rng, 1:length(r), length(r))]
            T_back_boot = copy(T_back_exp)
            T_back_boot[mask] .= T_sim .- r_boot   # (T_sim - r_boot) so (sim - obs) = r_boot
            push!(boot_ds, (scenario, t_exp, T_back_boot, T_back_std, n_reps))
        end
        boot_loss =
            θ -> begin
                # Build a bound-respecting loss; same closure args as the user's
                # loss but with the resampled dataset.
                out = 0.0
                for (k, (scenario, t_exp, T_back, T_back_std, n_reps)) in enumerate(boot_ds)
                    rv = residual_vector(
                        θ,
                        scenario,
                        t_exp,
                        T_back;
                        stage = stage,
                        t_skip = t_skip,
                        jacobian_strategy = jacobian_strategy,
                    )
                    rv === nothing && return 1e6
                    m = t_exp .> t_skip
                    σ = max.(T_back_std[m], 0.5)
                    out += 0.5 * mean((rv ./ σ) .^ 2)
                end
                out
            end
        try
            res = optimize(
                boot_loss,
                lower,
                upper,
                copy(θ_hat),
                Fminbox(NelderMead()),
                Optim.Options(iterations = refit_max_iter, show_trace = false),
            )
            θ_samples[:, b] = Optim.minimizer(res)
        catch err
            @warn "bootstrap fit $b failed" exception=err
            θ_samples[:, b] .= θ_hat
            n_failed += 1
        end
        verbose &&
            (b % max(1, div(n_boot, 10)) == 0) &&
            @info "bootstrap $b/$n_boot done ($n_failed failed so far)"
    end

    θ_med = [median(θ_samples[i, :]) for i = 1:p]
    Σ = cov(θ_samples; dims = 2)
    se = sqrt.(abs.(diag(Σ)))
    corr = Σ ./ (se * se')
    return (
        Σ = Σ,
        se = se,
        corr = corr,
        dof = max(n_boot - n_failed, 1),
        loss = NaN,
        method = :bootstrap,
        samples = θ_samples,
        theta_median = θ_med,
    )
end

# ============================================================================
# Stage 1 grid scan (cheap, optional)
# ============================================================================

function grid_scan_stage1(loss_fn; n = 5)
    mkpath(RUN_DIR[])
    path = joinpath(RUN_DIR[], "grid_scan.csv")
    open(path, "w") do io
        println(io, "k_virgin,k_char,loss")
        klo, khi = WP_K_VIRGIN_BOUNDS
        clo, chi = WP_K_CHAR_BOUNDS
        for kv in range(klo, khi, length = n), kc in range(clo, chi, length = n)
            l = loss_fn([kv, kc])
            @printf(io, "%.4f,%.4f,%.6e\n", kv, kc, l)
        end
    end
    @info "Grid scan written" path
end

# ============================================================================
# Output / plotting
# ============================================================================

function plot_fits(
    θ_hat,
    fit_dataset,
    val_dataset,
    stage::Int;
    outpath::AbstractString,
    jacobian_strategy::Union{Nothing,LinearSolveStrategy} = nothing,
)
    plots = Plots.Plot[]

    for (label_prefix, ds) in (("Fit", fit_dataset), ("Validation", val_dataset))
        for (scenario, t_exp, T_back_exp, _, _) in ds
            sol = build_and_solve(
                θ_hat,
                scenario;
                stage = stage,
                jacobian_strategy = jacobian_strategy,
            )

            # FSRI scenarios also report a partial top-surface temperature.
            # Add a front-surface panel ABOVE the back-surface one so the
            # pair reads top→bottom (front, then back).
            if startswith(scenario.label, "FSRI_")
                flux_label = replace(scenario.label, "FSRI_" => "")
                top = try
                    load_fsri_top_surface_temp(flux_label; verbose = false)
                catch err
                    @warn "FSRI top-surface load failed for $(scenario.label)" exception=err
                    nothing
                end
                p_front = plot(
                    xlabel = "Time [s]",
                    ylabel = "T_front [°C]",
                    title = "$(label_prefix): $(scenario.label) ($(round(scenario.q_ext/1000, digits=1)) kW/m²) — front",
                )
                if top !== nothing
                    scatter!(
                        p_front,
                        top.t,
                        top.T .- 273.15;
                        label = "experiment ($(length(top.t)) pts)",
                        color = :black,
                        ms = 4,
                        msw = 0,
                    )
                end
                if sol !== nothing
                    plot!(
                        p_front,
                        sol.t,
                        sol.surface_T .- 273.15;
                        label = "simulation",
                        lw = 2,
                        color = :red,
                        linestyle = :dash,
                    )
                end
                push!(plots, p_front)
            end

            p_back = plot(
                t_exp,
                T_back_exp .- 273.15;
                label = "experiment",
                lw = 2,
                color = :black,
                xlabel = "Time [s]",
                ylabel = "T_back [°C]",
                title = "$(label_prefix): $(scenario.label) ($(round(scenario.q_ext/1000, digits=1)) kW/m²) — back",
            )
            if sol !== nothing
                plot!(
                    p_back,
                    sol.t,
                    sol.extras.back_temperature .- 273.15;
                    label = "simulation",
                    lw = 2,
                    color = :red,
                    linestyle = :dash,
                )
            end
            push!(plots, p_back)
        end
    end

    final = plot(plots...; layout = (length(plots), 1), size = (900, 350 * length(plots)))
    savefig(final, outpath)
    @info "Fit plots written" outpath
end

function write_best_json(path, θ_hat, names, loss_val, cov_info; stage)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "{")
        println(
            io,
            "  \"timestamp\": \"$(Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))\",",
        )
        println(io, "  \"stage\": $stage,")
        println(io, "  \"loss\": $(loss_val),")
        println(
            io,
            "  \"parameter_names\": [" * join(("\"$n\"" for n in names), ", ") * "],",
        )
        println(
            io,
            "  \"theta_hat\": [" * join((@sprintf("%.6e", x) for x in θ_hat), ", ") * "],",
        )
        if cov_info !== nothing
            println(
                io,
                "  \"std_error\": [" *
                join((@sprintf("%.6e", x) for x in cov_info.se), ", ") *
                "],",
            )
            println(io, "  \"residual_dof\": $(cov_info.dof),")
            println(io, "  \"covariance_method\": \"$(cov_info.method)\",")
            println(io, "  \"correlation\": [")
            for i = 1:size(cov_info.corr, 1)
                row =
                    "    [" *
                    join((@sprintf("%.4f", x) for x in cov_info.corr[i, :]), ", ") *
                    "]"
                println(io, row * (i < size(cov_info.corr, 1) ? "," : ""))
            end
            println(io, "  ]")
        else
            println(io, "  \"std_error\": null")
        end
        println(io, "}")
    end
    @info "Best parameters written" path
end

# ============================================================================
# Driver
# ============================================================================

function main()
    stage = OPTS.stage
    cfg = stage_config(stage)

    # Per-run output directory: results/run_<YYYYMMDD_HHMMSS>_s<stage>_pw<weight>/
    # Every output (trace.csv, grid_scan.csv, best.json, fit_plots.png) lands
    # here, so back-to-back runs with different settings never overwrite.
    run_stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    pw_tag = replace(@sprintf("%.2f", OPTS.prior_weight), "." => "p")
    run_dir = joinpath(RESULTS_DIR, "run_$(run_stamp)_s$(stage)_pw$(pw_tag)")
    set_run_dir!(run_dir)

    jac_strat = OPTS.use_approx_solve ? ApproxSolve() : nothing
    @info "Inverse analysis configuration" stage=stage names=cfg.names θ_init=cfg.θ_init bounds=cfg.bounds threads=Threads.nthreads() run_dir=run_dir jacobian=OPTS.use_approx_solve ?
                                                                                                                                                               "ApproxSolve" :
                                                                                                                                                               "DirectSolve"

    # --- Load datasets via registry ---
    # Registry-driven loader: each flux carries its lab tag (UMD vs FSRI gas
    # vs FSRI CAPA), which determines BC stack, cone calibration, sample
    # thickness, MC, and replicate-file prefix. CLI `--fit=...`/`--val=...`
    # override the registry defaults; pass `--val=` to skip validation.
    fit_labels_used = select_labels(:fit, OPTS.fit_labels)
    val_labels_used = select_labels(:val, OPTS.val_labels)
    isempty(fit_labels_used) && error(
        "--fit= produced an empty dataset; need ≥1 fit flux. Available: $(REGISTRY_LABELS)",
    )
    overlap = intersect(fit_labels_used, val_labels_used)
    isempty(overlap) || error("flux(es) appear in both --fit and --val: $overlap")

    @info "Loading scenarios..." umd_cal=HEAT_FLUX_CAL fsri_cal=FSRI_HEAT_FLUX_CAL fsri_mc_pct_dry=3.0 fit=fit_labels_used val=val_labels_used
    fit_dataset = []
    for label in fit_labels_used
        push!(fit_dataset, load_registry_dataset(label))
    end

    val_dataset = []
    for label in val_labels_used
        try
            push!(val_dataset, load_registry_dataset(label))
        catch err
            @warn "Validation data not loadable, skipping" label exception=err
        end
    end

    # --- Loss closure ---
    # Load LFA parallel priors (T ≤ 473 K) and include as k(T) Gaussian prior
    lfa_priors = try
        load_lfa_parallel_priors()
    catch e
        @warn "Failed to load LFA priors" exception=e
        nothing
    end
    hfm_priors = try
        load_hfm_priors()
    catch e
        @warn "Failed to load HFM priors" exception=e
        nothing
    end

    # Merge priors (HFM + LFA parallel) and filter nil/empty
    k_prior = Any[]
    if hfm_priors !== nothing && !isempty(hfm_priors)
        append!(k_prior, hfm_priors)
        @info "Loaded $(length(hfm_priors)) HFM prior points (T ≤ 473 K)"
    end
    if lfa_priors !== nothing && !isempty(lfa_priors)
        append!(k_prior, lfa_priors)
        @info "Loaded $(length(lfa_priors)) LFA prior points (parallel, T ≤ 473 K)"
    end

    # k_prior_weight: 0 disables the prior; 1 weights it equally with the
    # total data term; default 0.5 puts it on par with one of the two fluxes.
    if isempty(k_prior) || OPTS.prior_weight == 0.0
        reason = isempty(k_prior) ? "no HFM/LFA priors found" : "--prior-weight=0"
        @info "Running without k_prior ($reason)"
        loss_fit = make_loss(
            fit_dataset;
            stage = stage,
            verbose = true,
            jacobian_strategy = jac_strat,
        )
    else
        @info "Using $(length(k_prior)) combined HFM+LFA k(T) prior points" weight=OPTS.prior_weight
        loss_fit = make_loss(
            fit_dataset;
            stage = stage,
            verbose = true,
            k_prior = k_prior,
            k_prior_weight = OPTS.prior_weight,
            jacobian_strategy = jac_strat,
        )
    end

    # --- Optional grid scan (Stage 1 only) ---
    # Grid output is also used as multi-start seeds (top-k minima → starts 2..k+1).
    grid_seeds = nothing
    if stage == 1 && OPTS.do_grid
        @info "Running 5×5 grid scan..."
        grid_scan_stage1(loss_fit; n = 5)
        grid_seeds = grid_scan_seeds(; k = 2)
        isempty(grid_seeds) ||
            @info "Seeded multi-start with $(length(grid_seeds)) grid minima"
    end

    # --- Optimize ---
    @info "Starting optimization..."
    reset_eval_counter!()
    results = run_optimization(
        loss_fit,
        cfg.θ_init,
        cfg.bounds;
        max_iter = OPTS.max_iter,
        n_starts = OPTS.n_starts,
        verbose = true,
        grid_seeds = grid_seeds,
    )

    best = results[1]
    @info "Best result" loss=best.loss theta=best.θ start=best.start_idx

    # --- Covariance ---
    cov_info = if OPTS.bootstrap > 0
        @info "Estimating parameter covariance via residual bootstrap..." n_boot=OPTS.bootstrap
        bootstrap_covariance(
            loss_fit,
            best.θ,
            fit_dataset,
            cfg.bounds;
            n_boot = OPTS.bootstrap,
            stage = stage,
            verbose = true,
            jacobian_strategy = jac_strat,
        )
    else
        @info "Estimating parameter covariance via Hessian (FD central diff)..."
        estimate_covariance(
            loss_fit,
            best.θ,
            fit_dataset;
            stage = stage,
            loss_at_theta = best.loss,
        )
    end

    # --- Reporting ---
    @info "── Final estimates ──"
    for (i, name) in enumerate(cfg.names)
        if cov_info !== nothing
            @printf("  %-10s = %.4g  ± %.2g\n", name, best.θ[i], cov_info.se[i])
        else
            @printf("  %-10s = %.4g\n", name, best.θ[i])
        end
    end

    # --- Save outputs ---
    write_best_json(
        joinpath(RUN_DIR[], "best.json"),
        best.θ,
        cfg.names,
        best.loss,
        cov_info;
        stage = stage,
    )
    plot_fits(
        best.θ,
        fit_dataset,
        val_dataset,
        stage;
        outpath = joinpath(RUN_DIR[], "fit_plots.png"),
        jacobian_strategy = jac_strat,
    )

    @info "Done."
end

# ============================================================================
# Run when executed as a script
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
