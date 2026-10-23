"""
Data loading and loss function for the white pine inverse analysis.

Experimental data layout (in `data/`):

  Back-temperature replicates (one file per replicate):
      UMD_Wood_CAPA_N2_<flux>kW_R<n>.csv
      Columns: `Time (s)`, `Mass (g)`, `TC back 1 (K)`
      Multiple replicates per flux are averaged at each timestep, with NaN
      values skipped.

  Boundary-condition replicates (one file per replicate):
      MaCFP_White_Pine_CAPA_<flux>kW_<date>_R<n>_detailed.csv
      Columns include `Time (s)`, `Inner Wall Temperature (K)`,
      `Outer Wall Temperature (K)`. Multiple replicates per flux are averaged
      at each timestep. Mapping into the model BCs (per user):
          Inner Wall  →  T_bot_gas  AND  T_bot_surr
          Outer Wall  →  T_top_surr
      `T_top_gas` (top-surface gas temperature) is not in these files and
      defaults to constant ambient until provided by the user.

  Optional `bcs_<flux>kW.csv` override:
      Columns: `time` (required) plus any subset of
      `T_top_gas`, `T_top_surr`, `T_bot_gas`, `T_bot_surr`.
      Any column present overrides the corresponding default. Use this to
      drop in the top-gas time-series once measured.

`load_scenario(flux_label, q_ext)` returns `(scenario, t_exp, T_back_exp,
T_back_std, n_replicates_at_t)` where the `T_*` series are replicate-averaged
targets and per-timepoint replicate statistics.

`loss(θ, dataset; stage)` runs the forward model for each scenario in the
dataset, interpolates the simulated `T_back` onto the experimental grid,
computes the σ-normalized weighted SSE, and writes one row per evaluation
to `results/trace.csv`.
"""

using DelimitedFiles
using Statistics
using Printf
using Dates
using Logging

include("wp_forward.jl")

const DATA_DIR = joinpath(@__DIR__, "data")
const RESULTS_DIR = joinpath(@__DIR__, "results")

# Per-run output directory. Defaults to RESULTS_DIR; the driver overrides
# this with a timestamped subdirectory so concurrent / sequential runs don't
# overwrite each other.
const RUN_DIR = Ref{String}(RESULTS_DIR)
trace_path() = joinpath(RUN_DIR[], "trace.csv")
set_run_dir!(path::AbstractString) =
    (RUN_DIR[] = String(path); mkpath(RUN_DIR[]); RUN_DIR[])
const FAILURE_PENALTY = 1.0e6
# Per-timepoint σ for the residual weighting. We use the between-replicate
# population std with a 0.5 K floor (typical TC noise lower bound). Single-
# replicate or zero-spread points fall back to the floor so they don't blow
# the loss up.
const DATA_SIGMA_FLOOR = 0.5  # K

# Prior uncertainty handling for k(T) data
# - floor terms prevent over-confident weights when HFM/LFA disagree
# - scale terms allow tuning each dataset without editing the raw files
const K_PRIOR_SIGMA_FLOOR_ABS = 1.0e-3
const K_PRIOR_SIGMA_FLOOR_REL = Dict(:hfm => 0.03, :lfa => 0.05)
const K_PRIOR_SIGMA_SCALE = Dict(:hfm => 1.0, :lfa => 1.0)

const T_AMBIENT_DEFAULT = 300.0    # K — fallback when no BC data is provided
const REPLICATE_PATTERN_PREFIX = "UMD_Wood_CAPA_N2"
const FSRI_REPLICATE_PATTERN_PREFIX = "FSRI_Wood_Gasification"
const FSRI_CAPA_REPLICATE_PATTERN_PREFIX = "FSRI_Wood_CAPA_N2"
const BC_DETAILED_PATTERN_PREFIX = "MaCFP_White_Pine_CAPA"
const BC_DETAILED_PATTERN_SUFFIX = "_detailed.csv"

# ============================================================================
# Top-surface gas temperature analytic model (user-provided)
# ============================================================================
#
#   T_top_gas(t; q̇") = T1·exp(T2·t) + T3·exp(T4·t) + T_HFG
#
# Coefficients tabulated against NOMINAL cone heat flux (kW/m²). Linear
# interpolation in q̇" between rows; flat extrapolation outside [30, 80].
# 50 and 70 kW/m² rows are themselves linear interpolations of (40, 60) and
# (60, 80) respectively — exactly the values supplied in the user's spec.

const T_TOP_GAS_TABLE = (
    # (q_nom_kW, T1 [K], T2 [1/s],  T3 [K], T4 [1/s],   T_HFG [K])
    (30.0, 48.1, 3.11e-8, -26.6, -1.28e-2, 290.0),
    (40.0, 66.5, 9.80e-6, -35.0, -1.29e-2, 290.0),
    (50.0, 78.5, 1.085e-5, -44.85, -1.094e-2, 290.0),
    (60.0, 90.5, 1.19e-5, -54.7, -8.98e-3, 290.0),
    (70.0, 103.65, 2.605e-5, -62.15, -9.64e-3, 290.0),
    (80.0, 116.8, 4.02e-5, -69.6, -1.03e-2, 290.0),
)

"""
Linearly interpolate the top-gas model coefficients in nominal heat flux
[kW/m²]. Outside the tabulated range coefficients clamp to the endpoints.
Returns `(T1, T2, T3, T4, T_HFG)`.
"""
function _top_gas_coeffs(q_nom_kw::Real)
    tbl = T_TOP_GAS_TABLE
    q = float(q_nom_kw)
    q <= tbl[1][1] && return tbl[1][2:end]
    q >= tbl[end][1] && return tbl[end][2:end]
    for i = 1:(length(tbl)-1)
        qa, qb = tbl[i][1], tbl[i+1][1]
        if qa <= q <= qb
            f = (q - qa) / (qb - qa)
            return ntuple(j -> tbl[i][j+1] + f * (tbl[i+1][j+1] - tbl[i][j+1]), 5)
        end
    end
    error("unreachable: q=$q outside T_TOP_GAS_TABLE")
end

"""
	top_gas_temperature_model(q_nom_kw)

Returns a callable `t -> K` for the top-surface gas temperature at the given
nominal heat flux (kW/m²), per the analytic two-exponential model:

	T(t; q̇") = T1·exp(T2·t) + T3·exp(T4·t) + T_HFG

with coefficients interpolated from `T_TOP_GAS_TABLE`.
"""
function top_gas_temperature_model(q_nom_kw::Real)
    T1, T2, T3, T4, T_HFG = _top_gas_coeffs(q_nom_kw)
    return t -> T1 * exp(T2 * t) + T3 * exp(T4 * t) + T_HFG
end

"""
Parse a flux label like "30kW" / "50kW" / "70.5kW" → nominal kW/m² as Float.
Returns `nothing` if the label does not match.
"""
function _parse_nominal_kw(flux_label::AbstractString)
    m = match(r"^([0-9]+(?:\.[0-9]+)?)kW$", flux_label)
    m === nothing ? nothing : parse(Float64, m.captures[1])
end

# ============================================================================
# CSV helpers
# ============================================================================

function _read_csv_with_header(path::AbstractString)
    isfile(path) || error("missing file: $path")
    # Some sources (e.g. FSRI gasification CSVs) ship CRLF and have blank cells
    # where one channel didn't record a sample. readdlm trips on both. Strip
    # CRs first, then parse cell-by-cell so blanks/non-numerics become NaN —
    # downstream loaders already mask with `isfinite`.
    buf = IOBuffer(replace(read(path, String), "\r" => ""))
    raw, header = readdlm(buf, ','; header = true)
    data = Matrix{Float64}(undef, size(raw)...)
    @inbounds for j in axes(raw, 2), i in axes(raw, 1)
        x = raw[i, j]
        data[i, j] = if x isa Real
            Float64(x)
        else
            s = strip(string(x))
            isempty(s) ? NaN : (something(tryparse(Float64, s), NaN))
        end
    end
    return (vec(strip.(string.(header))), data)
end

"""
Detect Celsius vs Kelvin (max < 200 → assumed °C) and convert to K.
Pure NaN columns pass through unchanged.
"""
function _to_kelvin(T::AbstractVector{<:Real})
    finite_T = filter(isfinite, T)
    isempty(finite_T) && return Float64.(T)
    maximum(finite_T) < 200.0 ? T .+ 273.15 : Float64.(T)
end

"""
Linear interpolating callable with flat extrapolation outside [t[1], t[end]].
"""
function flat_extrap_interp(t::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    tt = collect(Float64.(t))
    yy = collect(Float64.(y))
    @assert length(tt) == length(yy) ≥ 2 "interpolator needs ≥2 points"
    @assert issorted(tt) "time vector must be sorted"
    t_first, t_last = tt[1], tt[end]
    y_first, y_last = yy[1], yy[end]
    function (x)
        x <= t_first && return y_first
        x >= t_last && return y_last
        lo, hi = 1, length(tt)
        while hi - lo > 1
            mid = (lo + hi) >> 1
            if tt[mid] <= x
                lo = mid
            else
                hi = mid
            end
        end
        f = (x - tt[lo]) / (tt[hi] - tt[lo])
        return yy[lo] + f * (yy[hi] - yy[lo])
    end
end

# ============================================================================
# Replicate-averaged back-temperature target
# ============================================================================

"""
	load_back_temp_replicates(flux_label; data_dir=DATA_DIR, verbose=true)

Discover all `UMD_Wood_CAPA_N2_<flux_label>_R*.csv` files in `data_dir`, parse
them, and produce a single averaged back-temperature curve. NaN values are
dropped per-replicate. At each time, the average is taken over whichever
replicates are still valid at that time.

Returns a NamedTuple `(t, T, T_std, n_replicates_at_t, files)` where:
  - `t`                  Vector{Float64}, integer-second grid
  - `T`                  mean back-surface temperature [K]
  - `T_std`              between-replicate std at each time [K] (0 if n=1)
  - `n_replicates_at_t`  count of valid replicates at each time
  - `files`              list of source filenames (basenames)
"""
function load_back_temp_replicates(
    flux_label::AbstractString;
    data_dir::AbstractString = DATA_DIR,
    verbose::Bool = true,
    replicate_prefix::AbstractString = REPLICATE_PATTERN_PREFIX,
)

    pattern = Regex("^$(replicate_prefix)_$(flux_label)_R\\d+\\.csv\$")
    files = sort([joinpath(data_dir, f) for f in readdir(data_dir) if occursin(pattern, f)])
    isempty(files) && error("No replicate files matching $pattern in $data_dir")

    # Read each replicate
    reps = []  # vector of (t::Vector{Float64}, T::Vector{Float64}, name::String)
    for f in files
        headers, data = _read_csv_with_header(f)
        # "Time" matches `Time (s)`; "Timestamp" (FSRI) also satisfies startswith.
        t_col = findfirst(h -> startswith(h, "Time"), headers)
        # Match UMD `TC back 1`, `T_back`, and FSRI `Back Surface Temperature`.
        T_col = findfirst(h -> occursin(r"TC.*back|T[_ ]back|Back Surface"i, h), headers)
        if t_col === nothing || T_col === nothing
            error(
                "$(basename(f)): expected `Time*` and back-temperature columns, got $headers",
            )
        end
        t_raw = data[:, t_col]
        T_raw = _to_kelvin(data[:, T_col])
        valid = isfinite.(T_raw) .& isfinite.(t_raw)
        t_v = Float64.(t_raw[valid])
        T_v = Float64.(T_raw[valid])
        # Per-file diagnostic
        n_total = length(t_raw)
        n_valid = sum(valid)
        if verbose && n_valid < n_total
            frac_valid = round(100 * n_valid / n_total, digits = 1)
            @info "  $(basename(f)): $n_valid / $n_total valid points ($frac_valid%)"
        end
        push!(reps, (t = t_v, T = T_v, name = basename(f)))
    end

    # Replicate times are 1 Hz integer seconds; round to Int keys so the dict
    # lookup is robust to float-equality surprises if a fractional-time CSV
    # ever shows up. Guardrail asserts the rounding error stays sub-millisecond.
    _round_keys(ts) = begin
        out = Vector{Int}(undef, length(ts))
        @inbounds for i in eachindex(ts)
            τi = round(Int, ts[i])
            abs(ts[i] - τi) < 1e-3 ||
                error("non-integer-second time $(ts[i]) in replicate")
            out[i] = τi
        end
        out
    end
    all_times = sort(unique(reduce(vcat, _round_keys(r.t) for r in reps)))
    rep_lookup = [Dict{Int,Float64}(zip(_round_keys(r.t), r.T)) for r in reps]

    t_out = Float64[]
    T_out = Float64[]
    std_out = Float64[]
    n_out = Int[]
    for τ in all_times
        vals = Float64[]
        for d in rep_lookup
            haskey(d, τ) && push!(vals, d[τ])
        end
        isempty(vals) && continue
        push!(t_out, float(τ))
        push!(T_out, mean(vals))
        push!(std_out, length(vals) > 1 ? std(vals) : 0.0)
        push!(n_out, length(vals))
    end

    if verbose
        max_n = maximum(n_out)
        full_n = sum(n_out .== max_n)
        @info "  averaged $(length(reps)) replicate(s); $(full_n)/$(length(t_out)) timesteps have all $(max_n) replicates"
    end

    return (
        t = t_out,
        T = T_out,
        T_std = std_out,
        n_replicates_at_t = n_out,
        files = [basename(f) for f in files],
    )
end

# ============================================================================
# Replicate-averaged BC time-series (Inner / Outer wall temperatures)
# ============================================================================

"""
	load_bc_wall_temps(flux_label; data_dir=DATA_DIR, verbose=true)

Discover all `MaCFP_White_Pine_CAPA_<flux_label>_*_detailed.csv` files in
`data_dir`, parse them, and produce replicate-averaged Inner Wall and Outer
Wall temperature curves.

Returns `nothing` if no matching files exist; otherwise a NamedTuple
`(t, T_inner, T_outer, T_inner_std, T_outer_std, n_replicates_at_t, files)`.
The two wall series are averaged on the same shared time grid (which is the
union of replicate time vectors, restricted to timesteps where at least one
replicate has both columns valid).
"""
function load_bc_wall_temps(
    flux_label::AbstractString;
    data_dir::AbstractString = DATA_DIR,
    verbose::Bool = true,
)

    pattern = Regex(
        "^$(BC_DETAILED_PATTERN_PREFIX)_$(flux_label)_.*$(BC_DETAILED_PATTERN_SUFFIX)\$",
    )
    files = sort([joinpath(data_dir, f) for f in readdir(data_dir) if occursin(pattern, f)])
    isempty(files) && return nothing

    reps = []  # vector of (t, T_inner, T_outer, name)
    for f in files
        headers, data = _read_csv_with_header(f)
        t_col = findfirst(h -> startswith(h, "Time"), headers)
        inner_col = findfirst(h -> occursin(r"Inner Wall"i, h), headers)
        outer_col = findfirst(h -> occursin(r"Outer Wall"i, h), headers)
        if t_col === nothing || inner_col === nothing || outer_col === nothing
            error(
                "$(basename(f)): expected `Time`, `Inner Wall`, `Outer Wall` columns, got $headers",
            )
        end
        t_raw = data[:, t_col]
        Ti = _to_kelvin(data[:, inner_col])
        To = _to_kelvin(data[:, outer_col])
        valid = isfinite.(t_raw) .& isfinite.(Ti) .& isfinite.(To)
        push!(
            reps,
            (
                t = Float64.(t_raw[valid]),
                T_inner = Float64.(Ti[valid]),
                T_outer = Float64.(To[valid]),
                name = basename(f),
            ),
        )
    end

    # Round to integer-second keys (BC files share the 1 Hz grid).
    _round_keys(ts) = begin
        out = Vector{Int}(undef, length(ts))
        @inbounds for i in eachindex(ts)
            τi = round(Int, ts[i])
            abs(ts[i] - τi) < 1e-3 ||
                error("non-integer-second time $(ts[i]) in BC replicate")
            out[i] = τi
        end
        out
    end
    all_times = sort(unique(reduce(vcat, _round_keys(r.t) for r in reps)))
    in_lookup = [Dict{Int,Float64}(zip(_round_keys(r.t), r.T_inner)) for r in reps]
    ou_lookup = [Dict{Int,Float64}(zip(_round_keys(r.t), r.T_outer)) for r in reps]

    t_out = Float64[]
    Ti_out = Float64[]
    To_out = Float64[]
    Ti_std = Float64[]
    To_std = Float64[]
    n_out = Int[]
    for τ in all_times
        i_vals = Float64[]
        o_vals = Float64[]
        for k in eachindex(reps)
            if haskey(in_lookup[k], τ)
                push!(i_vals, in_lookup[k][τ])
                push!(o_vals, ou_lookup[k][τ])
            end
        end
        isempty(i_vals) && continue
        push!(t_out, float(τ))
        push!(Ti_out, mean(i_vals))
        push!(To_out, mean(o_vals))
        push!(Ti_std, length(i_vals) > 1 ? std(i_vals) : 0.0)
        push!(To_std, length(o_vals) > 1 ? std(o_vals) : 0.0)
        push!(n_out, length(i_vals))
    end

    if verbose
        Ti_lo, Ti_hi =
            round(minimum(Ti_out), digits = 1), round(maximum(Ti_out), digits = 1)
        To_lo, To_hi =
            round(minimum(To_out), digits = 1), round(maximum(To_out), digits = 1)
        @info "  averaged $(length(reps)) BC replicate(s) for $flux_label" Inner_K=(
            Ti_lo,
            Ti_hi,
        ) Outer_K=(To_lo, To_hi)
    end

    return (
        t = t_out,
        T_inner = Ti_out,
        T_outer = To_out,
        T_inner_std = Ti_std,
        T_outer_std = To_std,
        n_replicates_at_t = n_out,
        files = [basename(f) for f in files],
    )
end

# ============================================================================
# BC time-series loader (detailed CSVs + analytic top-gas + optional override)
# ============================================================================
#
# Resolution order for each of the four BC series (later overrides earlier):
#   1. Constant `T_default` (last-resort fallback; warning emitted)
#   2. Analytic top-gas model `top_gas_temperature_model(q_nom_kw)` for
#      `T_top_gas` (only this one series). Requires nominal flux to be known.
#   3. `MaCFP_White_Pine_CAPA_<flux_label>_*_detailed.csv` replicate-averaged
#      (Inner Wall → T_bot_gas, T_bot_surr; Outer Wall → T_top_surr).
#   4. `bcs_<flux_label>.csv` override (any subset of the four columns).

function _load_bcs_or_default(
    flux_label::AbstractString;
    data_dir::AbstractString = DATA_DIR,
    T_default::Real = T_AMBIENT_DEFAULT,
    q_nominal_kw::Union{Nothing,Real} = nothing,
    verbose::Bool = true,
)

    f_const = let T = float(T_default)
        _ -> T
    end
    bcs = Dict{String,Function}(
        "T_top_gas" => f_const,
        "T_top_surr" => f_const,
        "T_bot_gas" => f_const,
        "T_bot_surr" => f_const,
    )
    sources = Dict{String,String}(
        "T_top_gas" => "constant $(T_default) K",
        "T_top_surr" => "constant $(T_default) K",
        "T_bot_gas" => "constant $(T_default) K",
        "T_bot_surr" => "constant $(T_default) K",
    )

    # (2) Analytic top-gas model from the user's two-exponential fit.
    q_kw = q_nominal_kw === nothing ? _parse_nominal_kw(flux_label) : float(q_nominal_kw)
    if q_kw !== nothing
        bcs["T_top_gas"] = top_gas_temperature_model(q_kw)
        sources["T_top_gas"] = "analytic model @ q̇\"=$(q_kw) kW/m²"
    end

    # (3) Detailed-file Inner/Outer wall temperatures (replicate-averaged)
    wall = try
        load_bc_wall_temps(flux_label; data_dir = data_dir, verbose = verbose)
    catch err
        verbose &&
            @warn "Failed to load detailed BC files for $flux_label" exception=err
        nothing
    end
    if wall !== nothing
        T_inner = flat_extrap_interp(wall.t, wall.T_inner)
        T_outer = flat_extrap_interp(wall.t, wall.T_outer)
        bcs["T_bot_gas"] = T_inner
        bcs["T_bot_surr"] = T_inner
        bcs["T_top_surr"] = T_outer
        tag = "Inner Wall avg ($(length(wall.files)) file(s))"
        sources["T_bot_gas"] = tag
        sources["T_bot_surr"] = tag
        sources["T_top_surr"] = "Outer Wall avg ($(length(wall.files)) file(s))"
    end

    # (4) Optional bcs_<flux>kW.csv override (any subset of columns)
    bc_path = joinpath(data_dir, "bcs_$(flux_label).csv")
    if isfile(bc_path)
        headers, data = _read_csv_with_header(bc_path)
        t_idx = findfirst(==("time"), headers)
        t_idx === nothing && error("$bc_path: missing column `time`, got $headers")
        t_bc = data[:, t_idx]
        @assert issorted(t_bc) "$bc_path: time column must be sorted"
        any_override = false
        for name in ("T_top_gas", "T_top_surr", "T_bot_gas", "T_bot_surr")
            i = findfirst(==(name), headers)
            if i !== nothing
                bcs[name] = flat_extrap_interp(t_bc, _to_kelvin(data[:, i]))
                sources[name] = "override from $(basename(bc_path))"
                any_override = true
            end
        end
        any_override ||
            @warn "$(basename(bc_path)) found but no recognized BC columns — ignored" headers=headers
    end

    # Diagnostic on the resolved source for each series
    if verbose
        for name in ("T_top_gas", "T_top_surr", "T_bot_gas", "T_bot_surr")
            if startswith(sources[name], "constant")
                @warn "  BC $name: $(sources[name]) — no time-series resolved"
            else
                @info "  BC $name: $(sources[name])"
            end
        end
    end

    return (
        T_top_gas = bcs["T_top_gas"],
        T_top_surr = bcs["T_top_surr"],
        T_bot_gas = bcs["T_bot_gas"],
        T_bot_surr = bcs["T_bot_surr"],
    )
end

# ============================================================================
# Scenario loader
# ============================================================================

"""
	load_scenario(flux_label, q_ext; sample_thickness=8.5e-3, n_cells=160,
	              h_top=10.0, h_bot=5.0, verbose=true)

`flux_label` is the heat-flux portion of the filename (e.g. "30kW", "50kW").
Discovers and averages all replicates, optionally loads `bcs_<flux_label>.csv`,
returns `(scenario, t_exp, T_back_exp, T_back_std, n_replicates_at_t)` where
the `T_*` series are replicate-averaged targets and per-timepoint replicate
statistics.
"""
function load_scenario(
    flux_label::AbstractString,
    q_ext::Real;
    sample_thickness::Real = 8.5e-3,
    n_cells::Int = 160,
    h_top::Real = 8.6,
    h_bot::Real = 4.0,
    ε_bottom::Real = 0.92,
    verbose::Bool = true,
    use_variable_cross_section::Bool = true,
)

    verbose && @info "Loading scenario: $flux_label ($(round(q_ext/1000, digits=1)) kW/m²)"

    rep = load_back_temp_replicates(flux_label; verbose = verbose)
    t_exp = rep.t
    T_back_exp = rep.T
    T_back_std = rep.T_std
    n_reps_at_t = rep.n_replicates_at_t

    bcs = _load_bcs_or_default(flux_label; verbose = verbose)

    t_final = float(t_exp[end])

    # Anchor the simulation IC to the first experimental back-surface point so
    # sim and experiment share the same starting temperature for each flux.
    T_initial = T_back_exp[1]
    verbose &&
        @info "  initial sample temperature set to first back-surface point" T_initial_K=round(
            T_initial,
            digits = 2,
        ) t0_s=t_exp[1]

    scenario = Scenario(
        q_ext = q_ext,
        T_top_gas = bcs.T_top_gas,
        T_top_surr = bcs.T_top_surr,
        T_bot_gas = bcs.T_bot_gas,
        T_bot_surr = bcs.T_bot_surr,
        h_top = h_top,
        h_bot = h_bot,
        ε_bottom = ε_bottom,
        sample_thickness = sample_thickness,
        n_cells = n_cells,
        t_final = t_final,
        label = flux_label,
        T_initial = T_initial,
        use_variable_cross_section = use_variable_cross_section,
    )
    return (scenario, t_exp, T_back_exp, T_back_std, n_reps_at_t)
end

# ============================================================================
# "This Study" (UMD chamber) BC anchors
# ============================================================================
#
# Wall temperatures: treated as constant in time (the measured "This Study"
# traces are flat-ish above ambient over the run window). Linearly interpolate
# the ΔT-above-ambient between the two anchor fluxes. Top-gas temperature is
# the same two-exponential form as the CAPA-II model with anchor-specific
# coefficients; T_HFG = T_initial of the scenario.
#
# Anchors: 25 and 60 kW/m². Source: UMD "This Study" chamber characterization —
# wall ΔT-above-ambient from chart readings, and the two-exponential top-gas
# coefficients from the Table 3.1 gas-temperature correlation fits.

const THIS_STUDY_ANCHORS_KW = (25.0, 60.0)

# (ΔT_outer, ΔT_inner) [K above ambient]
const THIS_STUDY_DT_OUTER = (12.0, 29.0)   # outer wall  → T_top_surr
const THIS_STUDY_DT_INNER = (0.5, 1.75)   # inner wall  → T_bot_gas, T_bot_surr

# Top-gas: (T1 [K], T2 [1/s], T3 [K], T4 [1/s])
const THIS_STUDY_TGAS_COEFFS = (
    (50.6, 2.18e-5, -40.39, -3.46e-2),   # 25 kW
    (110.7, 3.43e-5, -85.4, -3.25e-2),   # 60 kW
)

@inline function _this_study_interp_flux(q_kw::Real, ya::Real, yb::Real)
    qa, qb = THIS_STUDY_ANCHORS_KW
    f = (float(q_kw) - qa) / (qb - qa)
    ya + f * (yb - ya)
end

"""
	this_study_top_gas(q_nom_kw; T_HFG)

Return a callable `t -> K` for the top-surface gas temperature at the given
nominal flux, using the "This Study" two-exponential coefficients linearly
interp/extrapolated in flux between 25 and 60 kW/m². `T_HFG` is the constant
offset (= T_initial of the scenario).
"""
function this_study_top_gas(q_nom_kw::Real; T_HFG::Real)
    a, b = THIS_STUDY_TGAS_COEFFS
    T1 = _this_study_interp_flux(q_nom_kw, a[1], b[1])
    T2 = _this_study_interp_flux(q_nom_kw, a[2], b[2])
    T3 = _this_study_interp_flux(q_nom_kw, a[3], b[3])
    T4 = _this_study_interp_flux(q_nom_kw, a[4], b[4])
    THFG = float(T_HFG)
    return t -> T1 * exp(T2 * t) + T3 * exp(T4 * t) + THFG
end

"""
	this_study_wall_temps(q_nom_kw; T_amb)

Return `(T_outer, T_inner, ΔT_outer, ΔT_inner)` constants by linearly
interpolating the "This Study" wall ΔT-above-ambient between the 25 and 60
kW/m² anchors. Outer maps to `T_top_surr`; inner maps to both `T_bot_gas` and
`T_bot_surr`.
"""
function this_study_wall_temps(q_nom_kw::Real; T_amb::Real)
    ΔTo = _this_study_interp_flux(q_nom_kw, THIS_STUDY_DT_OUTER[1], THIS_STUDY_DT_OUTER[2])
    ΔTi = _this_study_interp_flux(q_nom_kw, THIS_STUDY_DT_INNER[1], THIS_STUDY_DT_INNER[2])
    (T_outer = T_amb + ΔTo, T_inner = T_amb + ΔTi, ΔT_outer = ΔTo, ΔT_inner = ΔTi)
end

# ============================================================================
# FSRI top-surface temperature loader (partial column; mostly NaN after t≈25 s)
# ============================================================================

"""
	load_fsri_top_surface_temp(flux_label; data_dir=DATA_DIR, verbose=true)

Read the `Top Surface Temperature [K]` column from FSRI replicate file(s)
and return the replicate-averaged finite-only series. The column is sparse
(only the first ~25 s typically have valid samples), so the returned vectors
contain only the timesteps where at least one replicate reports a finite
value. Returns `nothing` if no FSRI files exist for the requested flux.
"""
function load_fsri_top_surface_temp(
    flux_label::AbstractString;
    data_dir::AbstractString = DATA_DIR,
    verbose::Bool = true,
)

    pattern = Regex("^$(FSRI_REPLICATE_PATTERN_PREFIX)_$(flux_label)_R\\d+\\.csv\$")
    files = sort([joinpath(data_dir, f) for f in readdir(data_dir) if occursin(pattern, f)])
    isempty(files) && return nothing

    reps = Tuple{Vector{Int},Vector{Float64}}[]
    for f in files
        headers, data = _read_csv_with_header(f)
        t_col = findfirst(h -> startswith(h, "Time"), headers)
        T_col = findfirst(h -> occursin(r"Top Surface"i, h), headers)
        (t_col === nothing || T_col === nothing) && continue
        t_raw = data[:, t_col]
        T_raw = _to_kelvin(data[:, T_col])
        valid = isfinite.(t_raw) .& isfinite.(T_raw)
        t_v = [round(Int, x) for x in t_raw[valid]]
        T_v = Float64.(T_raw[valid])
        push!(reps, (t_v, T_v))
    end
    isempty(reps) && return nothing

    all_t = sort(unique(reduce(vcat, r[1] for r in reps)))
    lookups = [Dict(zip(r[1], r[2])) for r in reps]
    t_out = Float64[]
    T_out = Float64[]
    std_out = Float64[]
    n_out = Int[]
    for τ in all_t
        vals = Float64[]
        for d in lookups
            haskey(d, τ) && push!(vals, d[τ])
        end
        isempty(vals) && continue
        push!(t_out, float(τ))
        push!(T_out, mean(vals))
        push!(std_out, length(vals) > 1 ? std(vals) : 0.0)
        push!(n_out, length(vals))
    end
    verbose &&
        @info "  FSRI top-surface temp ($flux_label): $(length(t_out)) finite samples ($(round(t_out[1], digits=1))–$(round(t_out[end], digits=1)) s)"
    return (t = t_out, T = T_out, T_std = std_out, n_replicates_at_t = n_out)
end

# ============================================================================
# FSRI scenario loader (uses This-Study BCs; back-temp from FSRI replicates)
# ============================================================================

"""
	load_fsri_scenario(flux_label, q_ext; sample_thickness=8.0e-3, n_cells=160,
	                   h_top=8.6, h_bot=4.0, ε_bottom=0.92, verbose=true,
	                   c0_moisture=0.0)

Build a fit/validation scenario from FSRI replicate file(s) matching
`FSRI_Wood_Gasification_<flux_label>_R*.csv`. Boundary conditions come from
the "This Study" anchors (constant walls, two-exponential top gas, flux
linearly interp/extrapolated between 25 and 60 kW). T_initial is set to the
first back-surface point (same convention as `load_scenario`), and T_HFG of
the top-gas model is anchored to that initial temperature.

`c0_moisture` (kg/m³) is the initial moisture concentration injected into
slot 8 of the component vector. Pass `0.03 * WP_RHO_DRY_BASIS` for 3% MC.
"""
function load_fsri_scenario(
    flux_label::AbstractString,
    q_ext::Real;
    sample_thickness::Real = 8.0e-3,
    n_cells::Int = 160,
    h_top::Real = 8.6,
    h_bot::Real = 4.0,
    ε_bottom::Real = 0.92,
    verbose::Bool = true,
    c0_moisture::Real = 0.0,
    replicate_prefix::AbstractString = FSRI_REPLICATE_PATTERN_PREFIX,
    use_variable_cross_section::Bool = true,
)

    verbose &&
        @info "Loading FSRI scenario: $flux_label ($(round(q_ext/1000, digits=1)) kW/m²)" prefix=replicate_prefix

    rep = load_back_temp_replicates(
        flux_label;
        verbose = verbose,
        replicate_prefix = replicate_prefix,
    )
    t_exp = rep.t
    T_back_exp = rep.T
    T_back_std = rep.T_std
    n_reps_at_t = rep.n_replicates_at_t

    q_kw = _parse_nominal_kw(flux_label)
    q_kw === nothing && error("could not parse nominal kW from `$flux_label`")

    T_initial = T_back_exp[1]
    walls = this_study_wall_temps(q_kw; T_amb = T_initial)
    T_top_surr = let T = walls.T_outer
        _ -> T
    end
    T_bot = let T = walls.T_inner
        _ -> T
    end
    T_top_gas = this_study_top_gas(q_kw; T_HFG = T_initial)

    mc_pct = 100 * c0_moisture / WP_RHO_DRY_BASIS
    verbose &&
        @info "  This-Study BCs resolved" T_initial_K=round(T_initial, digits = 2) T_top_surr_K=round(
            walls.T_outer,
            digits = 2,
        ) T_bot_K=round(walls.T_inner, digits = 2) ΔT_outer_K=round(
            walls.ΔT_outer,
            digits = 2,
        ) ΔT_inner_K=round(walls.ΔT_inner, digits = 2) c0_moisture_kg_per_m3=round(
            float(c0_moisture),
            digits = 2,
        ) mc_pct_dry=round(mc_pct, digits = 2)

    t_final = float(t_exp[end])
    scenario = Scenario(
        q_ext = q_ext,
        T_top_gas = T_top_gas,
        T_top_surr = T_top_surr,
        T_bot_gas = T_bot,
        T_bot_surr = T_bot,
        h_top = h_top,
        h_bot = h_bot,
        ε_bottom = ε_bottom,
        sample_thickness = sample_thickness,
        n_cells = n_cells,
        t_final = t_final,
        label = "FSRI_$flux_label",
        T_initial = T_initial,
        c0_moisture = c0_moisture,
        use_variable_cross_section = use_variable_cross_section,
    )
    return (scenario, t_exp, T_back_exp, T_back_std, n_reps_at_t)
end

# ============================================================================
# Residual helpers
# ============================================================================

"""
	interpolate_at(t_query, t_grid, y_grid)

Linear interpolation of `y_grid(t_grid)` evaluated at every `t_query` point,
with flat extrapolation outside `[t_grid[1], t_grid[end]]`.
"""
function interpolate_at(
    t_query::AbstractVector,
    t_grid::AbstractVector,
    y_grid::AbstractVector,
)
    itp = flat_extrap_interp(t_grid, y_grid)
    return [itp(t) for t in t_query]
end

"""
	residual_vector(θ, scenario, t_exp, T_back_exp;
	                stage, t_skip=1.0)

Run forward, interpolate simulated `T_back` to experimental times, and return
`(T_back_sim - T_back_exp)` over `t > t_skip`. Returns `nothing` on failure.
"""
function residual_vector(
    θ::AbstractVector,
    scenario,
    t_exp,
    T_back_exp;
    stage::Int = 1,
    t_skip::Real = 1.0,
    jacobian_strategy::Union{Nothing,LinearSolveStrategy} = nothing,
)

    sol = build_and_solve(
        θ,
        scenario;
        stage = stage,
        log_level = Logging.Warn,
        jacobian_strategy = jacobian_strategy,
    )
    sol === nothing && return nothing
    # `PyrolysisSolution.retcode` is a `Symbol` (`:Success` on success).
    sol.retcode == :Success || return nothing
    any(!isfinite, sol.extras.back_temperature) && return nothing

    mask = t_exp .> t_skip
    t_q = t_exp[mask]
    T_obs = T_back_exp[mask]
    T_sim = interpolate_at(t_q, sol.t, sol.extras.back_temperature)
    return T_sim .- T_obs
end

# ============================================================================
# Loss function with logging
# ============================================================================

mutable struct EvalCounter
    n::Int
end
const EVAL_COUNTER = EvalCounter(0)

function _ensure_trace_header()
    mkpath(RUN_DIR[])
    tp = trace_path()
    if !isfile(tp)
        open(tp, "w") do io
            println(
                io,
                "iter,timestamp,stage,theta,loss,data_loss,prior_loss,wallclock_s,success",
            )
        end
    end
end

"""
	loss(θ, dataset; stage=1, t_skip=1.0, weights=nothing,
	     log_trace=true, verbose=false, k_prior=nothing,
	     k_prior_weight=0.5, sigma_floor=DATA_SIGMA_FLOOR)

`dataset` is a vector of `(scenario, t_exp, T_back_exp, T_back_std,
n_replicates_at_t)` tuples. Loss is σ-normalized weighted MSE summed across
scenarios. Per-timepoint σ uses between-replicate population std with
`sigma_floor` (K) as a lower bound. Failed solves contribute `FAILURE_PENALTY`
per scenario; box-bound enforcement is delegated to the optimizer (Fminbox).

Scenarios are evaluated in parallel via `Threads.@spawn` when more than one
thread is available; pin runs to `julia -t auto` (or `-t N`) for the speedup.

`k_prior_weight` defaults to 0.5 so the combined prior + data terms have
roughly equal weight when each flux carries `1/N_flux` (default 0.5 each for
two fluxes). Both contributions are written to the trace.
"""
function loss(
    θ::AbstractVector,
    dataset;
    stage::Int = 1,
    t_skip::Real = 1.0,
    weights::Union{Nothing,AbstractVector} = nothing,
    log_trace::Bool = true,
    verbose::Bool = false,
    sigma_floor::Real = DATA_SIGMA_FLOOR,
    # k_prior: vector of tuples (phase, T, k_meas, sigma_k[, source])
    k_prior::Union{Nothing,AbstractVector} = nothing,
    k_prior_weight::Real = 0.5,
    jacobian_strategy::Union{Nothing,LinearSolveStrategy} = nothing,
)

    tstart = time()
    _ensure_trace_header()
    EVAL_COUNTER.n += 1

    w = weights === nothing ? fill(1.0 / length(dataset), length(dataset)) : weights
    @assert length(w) == length(dataset) "weights length must match dataset"

    # Run scenarios in parallel when multi-threaded; serial otherwise so we
    # don't pay spawn overhead on a single-threaded run.
    residuals = if Threads.nthreads() > 1 && length(dataset) > 1
        tasks = [
            Threads.@spawn residual_vector(
                θ,
                ds[1],
                ds[2],
                ds[3];
                stage = stage,
                t_skip = t_skip,
                jacobian_strategy = jacobian_strategy,
            ) for ds in dataset
        ]
        fetch.(tasks)
    else
        [
            residual_vector(
                θ,
                ds[1],
                ds[2],
                ds[3];
                stage = stage,
                t_skip = t_skip,
                jacobian_strategy = jacobian_strategy,
            ) for ds in dataset
        ]
    end

    data_contrib = 0.0
    all_success = true
    for (k, (_, t_exp, _, T_back_std, n_reps)) in enumerate(dataset)
        r = residuals[k]
        if r === nothing
            data_contrib += w[k] * FAILURE_PENALTY
            all_success = false
            continue
        end
        mask = t_exp .> t_skip
        @assert length(T_back_std) == length(t_exp) == length(n_reps)
        std_q = T_back_std[mask]
        @assert length(r) == length(std_q)
        # Between-replicate population std with a floor.
        σ = max.(std_q, sigma_floor)
        data_contrib += w[k] * mean((r ./ σ) .^ 2)
    end

    # k(T) Gaussian prior penalty (HFM/LFA) — normalized to mean-z² so the
    # scale matches a data flux term.
    prior_contrib = 0.0
    if k_prior !== nothing && !isempty(k_prior)
        acc = 0.0
        for tup in k_prior
            phase, T_i, k_meas, sigma_k, source = _prior_tuple(tup)
            sigma_eff = _prior_sigma(k_meas, sigma_k, source)
            k_mod = model_k_at_T(θ, T_i, stage, phase; int_idx = 1)
            acc += ((k_mod - k_meas) / sigma_eff)^2
        end
        prior_contrib = acc / length(k_prior)
    end

    total = data_contrib + k_prior_weight * prior_contrib
    wall = time() - tstart

    verbose && @printf(
        "[eval %4d] stage=%d θ=%s loss=%.6e data=%.6e prior=%.6e wall=%.1fs\n",
        EVAL_COUNTER.n,
        stage,
        _fmt_vec(θ),
        total,
        data_contrib,
        prior_contrib,
        wall
    )
    log_trace && _append_trace(
        EVAL_COUNTER.n,
        stage,
        θ,
        total,
        data_contrib,
        k_prior_weight * prior_contrib,
        wall,
        all_success,
    )
    return total
end


"""
	load_lfa_parallel_priors(; path=joinpath(DATA_DIR, "lfa_summary.csv"), Tmax=473.0)

Read `lfa_summary.csv` and return a vector of tuples
	( :virgin, T_K, k_mean_W_per_mK, k_std_W_per_mK, :lfa )
filtered to `orientation == "parallel"` and T <= Tmax.
"""
function load_lfa_parallel_priors(;
    path = joinpath(DATA_DIR, "lfa_summary.csv"),
    Tmax::Real = 473.0,
)
    isfile(path) || (@warn "LFA file not found at $path"; return nothing)
    # Use readdlm directly so non-numeric columns (e.g. booleans) don't break parsing
    raw, header = readdlm(path, ','; header = true)
    headers = vec(strip.(string.(header)))
    @debug "LFA loader: found $(size(raw,1)) rows, headers=$headers"

    idx_orient = findfirst(==("orientation"), headers)
    idx_T = findfirst(==("T_K"), headers)
    idx_k = findfirst(==("k_mean_W_per_mK"), headers)
    idx_ks = findfirst(==("k_std_W_per_mK"), headers)
    if idx_orient === nothing || idx_T === nothing || idx_k === nothing
        @warn "LFA loader: missing required columns. orientation=$idx_orient T_K=$idx_T k_mean=$idx_k k_std=$idx_ks"
        return nothing
    end
    priors = Any[]
    for i = 1:size(raw, 1)
        orient = nothing
        T = nothing
        try
            orient = strip(string(raw[i, idx_orient]))
            T = parse(Float64, string(raw[i, idx_T]))
        catch e
            @debug "LFA row $i: parse error $e"
            continue
        end
        if orient !== nothing && T !== nothing && orient == "parallel" && T <= Tmax
            # parse k and k_std if possible, otherwise use 5% fallback
            k_val = try
                parse(Float64, string(raw[i, idx_k]))
            catch

                nothing
            end
            ks_val = try
                parse(Float64, string(raw[i, idx_ks]))
            catch

                nothing
            end
            if k_val === nothing
                @debug "LFA row $i: k_val is nothing"
                continue
            end
            k = k_val
            ks = ks_val === nothing ? max(1e-3, 0.05 * abs(k)) : ks_val
            push!(priors, (:virgin, T, k, ks, :lfa))
            @debug "LFA row $i (orient=$orient T=$T): added prior (k=$k, ks=$ks)"
        end
    end
    @debug "LFA loader: returning $(length(priors)) priors"
    return priors
end

"""
	load_hfm_priors(; path=joinpath(DATA_DIR, "hfm_dry_summary.csv"), Tmax=473.0)

Read `hfm_dry_summary.csv` and return a vector of tuples
	( :virgin, T_K, k_mean_W_per_mK, k_std_W_per_mK, :hfm )
for rows with T <= Tmax.
"""
function load_hfm_priors(;
    path = joinpath(DATA_DIR, "hfm_dry_summary.csv"),
    Tmax::Real = 473.0,
)
    isfile(path) || return nothing
    headers, data = _read_csv_with_header(path)
    idx_T = findfirst(==("T_K"), headers)
    idx_k = findfirst(==("k_mean_W_per_mK"), headers)
    idx_ks = findfirst(==("k_std_W_per_mK"), headers)
    priors = Any[]
    for i = 1:size(data, 1)
        T = float(data[i, idx_T])
        if T <= Tmax
            k = float(data[i, idx_k])
            ks = idx_ks === nothing ? max(1e-3, 0.05 * abs(k)) : float(data[i, idx_ks])
            push!(priors, (:virgin, T, k, ks, :hfm))
        end
    end
    return priors
end

# ----- k(T) prior helpers -----

@inline function _prior_tuple(tup)
    if length(tup) == 4
        phase, T_i, k_meas, sigma_k = tup
        return phase, T_i, k_meas, sigma_k, :unknown
    elseif length(tup) == 5
        phase, T_i, k_meas, sigma_k, source = tup
        return phase, T_i, k_meas, sigma_k, source
    else
        error("k_prior tuple must be (phase, T, k, sigma[, source])")
    end
end

@inline function _prior_sigma(k_meas::Real, sigma_meas::Real, source::Symbol)
    # Combine reported uncertainty with dataset floor to avoid over-weighting
    base = sigma_meas > 0 ? sigma_meas : 0.0
    scale = get(K_PRIOR_SIGMA_SCALE, source, 1.0)
    floor_rel = get(K_PRIOR_SIGMA_FLOOR_REL, source, 0.0)
    floor_abs = K_PRIOR_SIGMA_FLOOR_ABS
    sigma_floor = max(floor_abs, floor_rel * abs(k_meas))
    sigma_eff = sqrt((base * scale)^2 + sigma_floor^2)
    return sigma_eff > 0 ? sigma_eff : max(1e-6, 0.05 * abs(k_meas))
end

# ----- internals -----

_fmt_vec(θ) = "[" * join((@sprintf("%.4g", x) for x in θ), ", ") * "]"

function _append_trace(
    iter::Int,
    stage::Int,
    θ::AbstractVector,
    l::Real,
    data_l::Real,
    prior_l::Real,
    wall::Real,
    success::Bool,
)
    open(trace_path(), "a") do io
        theta_str = "[" * join((@sprintf("%.6e", x) for x in θ), ";") * "]"
        ts = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        @printf(
            io,
            "%d,%s,%d,%s,%.6e,%.6e,%.6e,%.3f,%s\n",
            iter,
            ts,
            stage,
            theta_str,
            l,
            data_l,
            prior_l,
            wall,
            success
        )
    end
end

# ============================================================================
# Convenience wrappers
# ============================================================================

function make_loss(dataset; kwargs...)
    θ -> loss(θ, dataset; kwargs...)
end

function reset_eval_counter!()
    EVAL_COUNTER.n = 0
    # With per-run RUN_DIR, the trace file should not pre-exist for a fresh
    # run. Guard anyway: archive any leftover trace alongside it so a re-run
    # of main() in the same Julia session doesn't append into a stale file.
    tp = trace_path()
    if isfile(tp)
        stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        mv(tp, joinpath(RUN_DIR[], "trace_$(stamp).csv"))
    end
end

# ============================================================================
# The scenarios
# ============================================================================

# Cone heater calibration: actual delivered flux is 97.3% of the nominal
# setpoint for the UMD CAPA cone. Apply to nominal labels (30/50/70 kW/m²)
# when constructing scenarios so the simulation sees the real radiative input.
const HEAT_FLUX_CAL = 0.973

# FSRI 20 kW gasification cone uses a 0.98 area-averaged multiplier per
# the FSRI calibration documentation.
const FSRI_HEAT_FLUX_CAL = 0.98

# Enable the variable-cross-section (radial shrinkage) law for the inverse fit.
# It is NOT free of the fitted trajectories: the ALE solver conserves cell volume
# V = A·Δz, and the density-driven volumetric strain rate θ_cell is A-independent,
# so the z-node velocity is the *residual* axial strain θ_L = θ_cell − θ_A (see
# src/Residual/residual.jl `_compute_ale_mesh_velocity!`). With the law off the
# radial channel θ_A ≡ 0 and the entire volume loss is forced into the thickness,
# making the sample ~13–17% thinner in z at end-of-test (A/A₀ ≈ 0.86 at χ̄≈0.65–
# 0.75) than the validated radial geometry — a shorter conduction path that biases
# back-surface temperatures and hence the fitted properties. The naive "q″·A and
# k·A·∂T/∂z both scale with A so it cancels" argument holds only at *fixed*
# geometry; A enters the geometry evolution, which is exactly where the law acts.
# The radial channel is calibrated/validated (validate_thickness.jl) against
# the experimental conversion χ̄=1−m_f/m_0; its
# fit residual sits at the ±0.013 D-measurement noise floor (see WP_SR2 block),
# so `true` is the more physical setting. Cost is per-eval only (χ-block
# tracking + a denser ∂physics/∂χⱼ Jacobian coupling).
const INVERSE_USE_VARIABLE_CROSS_SECTION = true

# ============================================================================
# Dataset registry — single source of truth for which fluxes exist and how
# each one is loaded. Each entry knows its lab (which sets BC stack, cone
# calibration, sample thickness, MC, replicate-file prefix) and a
# `default_role` that says whether the flux defaults to fit / val / off.
#
# CLI flags `--fit=<csv>` and `--val=<csv>` override the defaults — pass any
# subset of the labels below. Pass `--val=` (empty) to skip validation.
# ============================================================================

const FSRI_C0_MOISTURE = 0.03 * WP_RHO_DRY_BASIS   # 10.8 kg/m³
const FSRI_SAMPLE_THICK = 8.68e-3
const UMD_SAMPLE_THICK = 8.5e-3   # `load_scenario` default

# Order matters: this is the order in which defaults get assembled.
const DATASET_REGISTRY = [
    "20kW" => (lab = :fsri_gas, q_nom_kw = 20.0, default_role = :fit),
    "30kW" => (lab = :umd, q_nom_kw = 30.0, default_role = :fit),
    "40kW" => (lab = :fsri_capa, q_nom_kw = 40.0, default_role = :fit),
    "50kW" => (lab = :umd, q_nom_kw = 50.0, default_role = :fit),
    "60kW" => (lab = :fsri_capa, q_nom_kw = 60.0, default_role = :fit),
    "70kW" => (lab = :umd, q_nom_kw = 70.0, default_role = :val),
]

const REGISTRY_LABELS = String[p.first for p in DATASET_REGISTRY]
_registry_lookup(label) = first(filter(p -> p.first == label, DATASET_REGISTRY)).second

function load_registry_dataset(
    label::AbstractString;
    verbose::Bool = true,
    use_variable_cross_section::Bool = INVERSE_USE_VARIABLE_CROSS_SECTION,
)
    label in REGISTRY_LABELS ||
        error("unknown dataset $label; choose from $(REGISTRY_LABELS)")
    spec = _registry_lookup(label)
    if spec.lab == :umd
        return load_scenario(
            label,
            spec.q_nom_kw * 1e3 * HEAT_FLUX_CAL;
            sample_thickness = UMD_SAMPLE_THICK,
            verbose = verbose,
            use_variable_cross_section = use_variable_cross_section,
        )
    elseif spec.lab == :fsri_gas
        return load_fsri_scenario(
            label,
            spec.q_nom_kw * 1e3 * FSRI_HEAT_FLUX_CAL;
            c0_moisture = FSRI_C0_MOISTURE,
            sample_thickness = FSRI_SAMPLE_THICK,
            replicate_prefix = FSRI_REPLICATE_PATTERN_PREFIX,
            verbose = verbose,
            use_variable_cross_section = use_variable_cross_section,
        )
    elseif spec.lab == :fsri_capa
        return load_fsri_scenario(
            label,
            spec.q_nom_kw * 1e3 * FSRI_HEAT_FLUX_CAL;
            c0_moisture = FSRI_C0_MOISTURE,
            sample_thickness = FSRI_SAMPLE_THICK,
            replicate_prefix = FSRI_CAPA_REPLICATE_PATTERN_PREFIX,
            verbose = verbose,
            use_variable_cross_section = use_variable_cross_section,
        )
    else
        error("registry: unknown lab tag $(spec.lab) for $label")
    end
end
