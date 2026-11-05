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

# Errors at or below this are drawn at it on the logarithmic axis (exact runs).
const COST_ERROR_FLOOR = 1e-14

# Each code keeps its colour and marker on every figure (fixed by identity).
const COST_STYLE = Dict(
    "this code" => (color = RGBf(0.0, 0.447, 0.698), marker = :circle),
    "fds" => (color = RGBf(0.902, 0.624, 0.0), marker = :utriangle),
    "gpyro" => (color = RGBf(0.0, 0.620, 0.451), marker = :rect),
    "thermakin" => (color = RGBf(0.800, 0.475, 0.655), marker = :diamond),
)

"""
Cost against accuracy for the matched external runs of a case
(`matched_costs`). External runs scored on the same curves and references
share one panel with this code's configurations, which are the same for all
of them; a case whose runs are scored on different references (V-CPL-03's
paths) gets a panel per reference set, all on the same axes. Filled markers
are each external code's matched mesh and step and this code's report
tolerances; open markers are the other runs. This code's line joins its
cheapest configuration at each accuracy. Lower left is better. A run every code solves exactly (all errors at the floor) is not drawn;
`nothing` when no run is left.
"""
function render_cost(case_id::String, panels, path::String)
    exact(p) = all(q -> q.error <= COST_ERROR_FLOOR, vcat(p.external, p.ladder))
    panels = filter(!exact, panels)
    isempty(panels) && return nothing
    keys_ = unique(p.refs for p in panels)
    groups = [filter(p -> p.refs == k, panels) for k in keys_]
    ncol = min(length(groups), 3)
    nrow = cld(length(groups), ncol)
    size = length(groups) == 1 ? FIG_SIZE : (380 * ncol + 120, 340 * nrow + 140)
    fig = Figure(; size)
    errors =
        [max(q.error, COST_ERROR_FLOOR) for p in panels for q in vcat(p.external, p.ladder)]
    times = [q.time for p in panels for q in vcat(p.external, p.ladder)]
    decades(v) =
        (a = floor(Int, log10(minimum(v))); (a, max(a + 1, ceil(Int, log10(maximum(v))))))
    xa, xb = decades(errors)
    ya, yb = decades(times)
    ticks(a, b) = (
        k = b - a > 7 ? 2 : 1;
        (10.0 .^ (a:k:b), [rich("10", superscript(string(e))) for e = a:k:b])
    )
    axes = Axis[]
    for (g, group) in enumerate(groups)
        row, col = cld(g, ncol), mod1(g, ncol)
        ax = Axis(
            fig[row, col];
            title = length(groups) > 1 ? join(last.(group[1].refs), ", ") : "",
            xscale = log10,
            yscale = log10,
            xticks = ticks(xa, xb),
            yticks = ticks(ya, yb),
        )
        col > 1 && hideydecorations!(ax; grid = false, ticks = false)
        row < nrow &&
            length(groups) - g >= ncol &&
            hidexdecorations!(ax; grid = false, ticks = false)
        ours = group[1].ladder
        st = COST_STYLE["this code"]
        e(q) = max(q.error, COST_ERROR_FLOOR)
        # Cheapest configuration at each accuracy: nothing at least as accurate is faster.
        front = sort(
            filter(q -> !any(p -> e(p) <= e(q) && p.time < q.time, ours), ours);
            by = e,
        )
        lines!(ax, e.(front), [q.time for q in front]; color = st.color, linewidth = 2)
        for q in ours
            report = q.tol == 1
            scatter!(
                ax,
                [e(q)],
                [q.time];
                marker = st.marker,
                markersize = 9,
                color = report ? st.color : (:white, 0.0),
                strokecolor = st.color,
                strokewidth = report ? 0 : 1.4,
            )
        end
        for p in group
            st = COST_STYLE[p.solver]
            for q in p.external
                base = q.level == "base"
                scatter!(
                    ax,
                    [max(q.error, COST_ERROR_FLOOR)],
                    [q.time];
                    marker = st.marker,
                    markersize = base ? 13 : 11,
                    color = base ? st.color : (:white, 0.0),
                    strokecolor = st.color,
                    strokewidth = base ? 0 : 1.6,
                )
            end
        end
        limits!(ax, 10.0^(xa - 0.2), 10.0^(xb + 0.2), 10.0^(ya - 0.2), 10.0^(yb + 0.2))
        push!(axes, ax)
    end
    Label(fig[1:nrow, 0], "solve time (s)"; rotation = pi / 2, tellheight = false)
    Label(
        fig[nrow+1, 1:ncol],
        "largest error over the scored references";
        tellwidth = false,
    )
    solvers = unique(p.solver for p in panels)
    labels =
        Dict("fds" => "FDS 6.11.1", "gpyro" => "Gpyro V2", "thermakin" => "ThermaKin2Ds 7")
    Legend(
        fig[nrow+2, 1:ncol],
        vcat(
            [[
                LineElement(color = COST_STYLE["this code"].color, linewidth = 2),
                MarkerElement(
                    color = COST_STYLE["this code"].color,
                    marker = :circle,
                    markersize = 9,
                ),
            ]],
            [
                MarkerElement(
                    color = COST_STYLE[s].color,
                    marker = COST_STYLE[s].marker,
                    markersize = 12,
                ) for s in solvers
            ],
        ),
        vcat(["Pyrolysis.jl"], [labels[s] for s in solvers]);
        orientation = :horizontal,
        framevisible = false,
        tellwidth = false,
    )
    colgap!(fig.layout, 12)
    rowgap!(fig.layout, 8)
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
