# =============================================================================
# Rendering layer for the verification report runner (CairoMakie).
# =============================================================================
#
# All plotting lives here; the harness hands over the plain data structs from
# test/verification/casedefs.jl (`ReportCurve`/`CurveSeries`). Three figure
# kinds make up the standard artifact set:
#   1. overlay          — solver values as markers over densely-sampled exact
#                         reference lines (FDS Verification Guide style)
#   2. pointwise error  — y − exact(x) per series, shows WHERE the
#                         discretization struggles (norms alone hide this)
#   3. convergence      — log–log error norms vs h with fitted observed order
#                         and an expected-order reference slope
# Styling is fixed here (single theme, fixed size) so regenerated artifacts
# are deterministic and visually uniform across cases.

using CairoMakie

const FIG_SIZE = (880, 560)
const N_DENSE = 400          # samples for exact reference lines
const MAX_MARKERS = 60       # overlay markers per series (CSVs keep every point)

# Thin dense numeric series so markers don't bury the exact reference line;
# always keeps the first and last point.
function _marker_stride(v)
    step = max(1, length(v) ÷ MAX_MARKERS)
    idx = collect(1:step:length(v))
    idx[end] == length(v) || push!(idx, length(v))
    return idx
end

report_theme() = Theme(
    fontsize = 16,
    Axis = (
        xgridvisible = true,
        ygridvisible = true,
        titlealign = :left,
        subtitlefont = :italic,
    ),
)

_axis_scale(s::Symbol) = s === :log10 ? log10 : identity
_series_color(i) = Makie.wong_colors()[mod1(i, length(Makie.wong_colors()))]

_has_exact(curve::ReportCurve) = any(s -> s.exact !== nothing, curve.series)

# Cross-code (third-party) series are context, not focus: our series keep the
# full categorical palette; each other code gets one gray shade + one marker
# shape (identity by shape + legend, never color alone). Keyed by the label's
# leading word ("FDS", "Gpyro", …) in order of first appearance.
const CROSSCODE_GRAYS =
    [RGBf(0.40, 0.40, 0.40), RGBf(0.65, 0.65, 0.65), RGBf(0.52, 0.52, 0.52)]
const CROSSCODE_MARKERS = [:utriangle, :rect, :diamond, :cross]
_code_key(label::String) = first(split(label, (' ', ',')))

function _crosscode_style(series::Vector{CurveSeries})
    keys = String[]
    for s in series
        s.crosscode || continue
        k = _code_key(s.label)
        k in keys || push!(keys, k)
    end
    return Dict(
        k => (
            color = CROSSCODE_GRAYS[mod1(i, length(CROSSCODE_GRAYS))],
            marker = CROSSCODE_MARKERS[mod1(i, length(CROSSCODE_MARKERS))],
        ) for (i, k) in enumerate(keys)
    )
end

"""
Overlay: numeric markers over dense exact lines. Returns the saved path.

`crosscode = false` (default) draws only our series — the clean per-case
figure. `crosscode = true` additionally draws the third-party series
(`s.crosscode`), muted gray and *behind* ours so our data stays the focus.
"""
function render_overlay(curve::ReportCurve, path::String; crosscode::Bool = false)
    fig = Figure(size = FIG_SIZE)
    ours = [s for s in curve.series if !s.crosscode]
    codes = crosscode ? [s for s in curve.series if s.crosscode] : CurveSeries[]
    ax = Axis(
        fig[1, 1];
        title = curve.title,
        subtitle = if !_has_exact(curve)
            ""
        elseif isempty(codes)
            "lines: exact/reference — markers: solver"
        else
            "lines: exact — color: this solver — gray: other codes"
        end,
        xlabel = curve.xlabel,
        ylabel = curve.ylabel,
        xscale = _axis_scale(curve.xscale),
        yscale = _axis_scale(curve.yscale),
    )
    style = _crosscode_style(curve.series)
    for s in codes   # context layer first, so our markers draw on top
        st = style[_code_key(s.label)]
        idx = _marker_stride(s.x)
        scatter!(
            ax,
            s.x[idx],
            s.y[idx];
            color = st.color,
            marker = st.marker,
            markersize = 9,
            label = s.label,
        )
    end
    for (i, s) in enumerate(ours)
        color = _series_color(i)
        if s.exact !== nothing
            xd = collect(range(minimum(s.x), maximum(s.x); length = N_DENSE))
            lines!(ax, xd, Float64[s.exact(x) for x in xd]; color = color)
        end
        idx = _marker_stride(s.x)
        scatter!(
            ax,
            s.x[idx],
            s.y[idx];
            color = color,
            markersize = 8,
            strokewidth = 0.5,
            strokecolor = :black,
            label = s.label,
        )
    end
    (length(ours) + length(codes)) > 0 && Legend(fig[1, 2], ax)
    save(path, fig)
    return path
end

"""
Pointwise-error panel for every series carrying an `exact` callable:
e(x) = y − exact(x) at the numeric abscissa. Returns the saved path, or
`nothing` when no series has an exact reference.

`crosscode` as in `render_overlay`: false → our error only, linear *signed*
axis (structure and bias direction are the point); true → other codes' errors
added as muted dashed-gray context behind ours, on a log₁₀ |error| axis — the
cross-code comparison is about magnitudes spanning decades, which a linear
axis flattens (sign is still visible on the base variant). |e| is floored at
1e-12 so exact zeros survive the log.
"""
function render_pointwise_error(curve::ReportCurve, path::String; crosscode::Bool = false)
    ours = [s for s in curve.series if !s.crosscode && s.exact !== nothing]
    codes =
        crosscode ? [s for s in curve.series if s.crosscode && s.exact !== nothing] :
        CurveSeries[]
    isempty(ours) && isempty(codes) && return nothing
    err(s, idx) = begin
        e = s.y[idx] .- Float64[s.exact(x) for x in s.x[idx]]
        crosscode ? max.(abs.(e), 1e-12) : e
    end
    fig = Figure(size = FIG_SIZE)
    ax = Axis(
        fig[1, 1];
        title = "$(curve.title) — pointwise error",
        subtitle = isempty(codes) ? "" :
                   "color: this solver — gray: other codes — log |error|",
        xlabel = curve.xlabel,
        ylabel = crosscode ? "|solver − exact|" : "error (solver − exact)",
        xscale = _axis_scale(curve.xscale),
        yscale = crosscode ? log10 : identity,
    )
    crosscode || hlines!(ax, [0.0]; color = (:gray, 0.6), linestyle = :dash)
    style = _crosscode_style(curve.series)
    for s in codes   # context layer first, thinned so it can't bury ours
        st = style[_code_key(s.label)]
        idx = _marker_stride(s.x)
        scatterlines!(
            ax,
            s.x[idx],
            err(s, idx);
            color = st.color,
            marker = st.marker,
            linestyle = :dash,
            markersize = 7,
            label = s.label,
        )
    end
    for (i, s) in enumerate(ours)
        scatterlines!(
            ax,
            s.x,
            err(s, eachindex(s.x));
            color = _series_color(i),
            markersize = 7,
            label = s.label,
        )
    end
    Legend(fig[1, 2], ax)
    save(path, fig)
    return path
end

"""
Setup-and-simulation comparison: our warmed build-and-solve wall time at every
ladder rung vs the other codes' internal total elapsed wall time, as a log-x
lollipop (bars can't start at 0 on a log axis).
Ours colored, other codes muted gray — the reference runs differ in
resolution and simulated span (see reference/timings.csv notes), so this is
indicative, not a controlled benchmark; the subtitle says so.
"""
function render_solve_time(
    case_id::String,
    ns::Vector{Int},
    walls::Vector{Float64},
    refs,
    path::String,
)
    labels = vcat(
        ["this solver, warm (n = $(n))" for n in ns],
        String[r.label * " — " * r.res for r in refs],
    )
    values = vcat(walls, Float64[r.wall for r in refs])
    code_keys = unique(String[_code_key(r.label) for r in refs])
    colors = vcat(
        fill(_series_color(1), length(ns)),
        [
            CROSSCODE_GRAYS[mod1(
                findfirst(==(_code_key(r.label)), code_keys),
                length(CROSSCODE_GRAYS),
            )] for r in refs
        ],
    )
    N = length(values)
    ys = collect(Float64, N:-1:1)   # ours on top
    x0 = minimum(values) / 3
    fig = Figure(size = FIG_SIZE)
    ax = Axis(
        fig[1, 1];
        title = "$(case_id): setup and simulation wall time",
        subtitle = "warm Julia vs internal elapsed; indicative (see timing protocol)",
        xlabel = "elapsed wall time (s)",
        xscale = log10,
        yticks = (ys, labels),
    )
    linesegments!(
        ax,
        [Point2f(x0, y) => Point2f(v, y) for (v, y) in zip(values, ys)];
        color = colors,
        linewidth = 3,
    )
    scatter!(ax, values, ys; color = colors, markersize = 14)
    for (v, y) in zip(values, ys)
        text!(
            ax,
            v,
            y;
            text = " $(round(v; sigdigits = 3)) s",
            align = (:left, :center),
            fontsize = 11,
        )
    end
    xlims!(ax, x0, maximum(values) * 4)
    save(path, fig)
    return path
end

"Least-squares observed order: slope of ln(err) vs ln(h)."
function fit_order(hs::Vector{Float64}, errs::Vector{Float64})
    lx, ly = log.(hs), log.(errs)
    mx, my = mean(lx), mean(ly)
    return sum((lx .- mx) .* (ly .- my)) / sum(abs2, lx .- mx)
end

"""
Log–log convergence plot: one scatterline per norm in `norm_names` (columns of
`norms`, a Vector of NamedTuples, one per rung), fitted observed order in the
legend, and a dashed expected-order reference anchored at the finest rung of
the first norm. `hs` is the mesh-size proxy (1/n, or Δt for temporal sweeps).
Nonpositive errors (e.g. nodally exact cases) are dropped with a warning.
Returns the saved path, or `nothing` if nothing was plottable.
"""
function render_convergence(
    title::String,
    hs::Vector{Float64},
    norms::Vector{<:NamedTuple},
    expected_order::Union{Nothing,Float64},
    path::String;
    xlabel::String = "h ∝ 1/n_cells",
)
    isempty(norms) && return nothing
    norm_names = collect(keys(norms[1]))
    isempty(norm_names) && return nothing
    fig = Figure(size = FIG_SIZE)
    ax = Axis(
        fig[1, 1];
        title = title,
        xlabel = xlabel,
        ylabel = "error norm",
        xscale = log10,
        yscale = log10,
    )
    plotted = false
    for (i, name) in enumerate(norm_names)
        errs = Float64[getproperty(nt, name) for nt in norms]
        keep = errs .> 0
        if !all(keep)
            @warn "dropping nonpositive errors from convergence plot" title name
        end
        count(keep) >= 2 || continue
        p = fit_order(hs[keep], errs[keep])
        scatterlines!(
            ax,
            hs[keep],
            errs[keep];
            color = _series_color(i),
            markersize = 9,
            label = "$(name) (observed order $(round(p, digits = 2)))",
        )
        plotted = true
    end
    plotted || return nothing
    if expected_order !== nothing
        name = norm_names[1]
        errs = Float64[getproperty(nt, name) for nt in norms]
        anchor = findlast(>(0), errs)
        if anchor !== nothing
            href = sort(hs)
            eref = errs[anchor] .* (href ./ hs[anchor]) .^ expected_order
            lines!(
                ax,
                href,
                eref;
                color = :gray,
                linestyle = :dash,
                label = "expected order $(expected_order)",
            )
        end
    end
    Legend(fig[1, 2], ax)
    save(path, fig)
    return path
end

"""
Parameter-sweep plot (e.g. V-CPL-02's activation-energy ladder): each measure
returned by the sweep's `measure(sol)` as a scatterline over the swept values.
"""
function render_param_sweep(
    title::String,
    param_name::String,
    values::Vector{Float64},
    measures::Vector{<:NamedTuple},
    path::String,
)
    isempty(measures) && return nothing
    fig = Figure(size = FIG_SIZE)
    ax = Axis(fig[1, 1]; title = title, xlabel = param_name, ylabel = "measure")
    for (i, name) in enumerate(keys(measures[1]))
        scatterlines!(
            ax,
            values,
            Float64[getproperty(nt, name) for nt in measures];
            color = _series_color(i),
            markersize = 9,
            label = string(name),
        )
    end
    Legend(fig[1, 2], ax)
    save(path, fig)
    return path
end
