# The MaCFP PMMA validation pages: every code's runs of the MaCFP property sets
# against the institutions' submissions, against each other, and against the
# NIST gasification apparatus experiments.
#
#   julia --project=scripts/macfp_pmma scripts/macfp_pmma/report.jl
#
# Reads results/ (run_pyrolysis.jl, run_apparatus.jl, the external runners and,
# when present, a run_timing.py campaign in results/timing/) and the pinned
# matl-db checkout. Writes the overview docs/guide/validation/macfp_pmma.md, its
# detail pages in docs/guide/validation/macfp_pmma/ and the figures in
# docs/guide/validation/assets/macfp_pmma/. The same results give the same
# pages, byte for byte.
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, CairoMakie, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PAGE = joinpath(ROOT, "docs", "guide", "validation", "macfp_pmma.md")
const ASSETS = "assets/macfp_pmma"
const FIG = joinpath(dirname(PAGE), ASSETS)
const TIMING = joinpath(RESULTS, "timing")

const INK = colorant"#0b0b0b"
const GRID = colorant"#e4e3df"
const EXP = colorant"#9a9890"
const BLUE = colorant"#2a78d6"

# The runs on the page: (source, label, colour, line style).
const SOURCES = [
    ("pyrolysis/base", "Pyrolysis.jl", BLUE, :solid),
    ("fds/base", "FDS 6.11.1", colorant"#eb6834", (:dot, :dense)),
    ("gpyro/base", "Gpyro V2", colorant"#1baf7a", :dashdot),
    ("thermakin/base", "ThermaKin2Ds MAX", colorant"#e87ba4", (:dash, :loose)),
    ("thermakin/rand", "ThermaKin2Ds RAND", colorant"#6250d6", (:dot, :loose)),
]
# The variant that reproduces what a group evidently ran.
const AS_RAN = Dict(
    "BUW-FZJ_A" => "h5.33",
    "BUW-FZJ_B" => "h5.33",
    "DBI_1" => "tablesC",
    "Aalto_I" => "T303",
)
const AS_RAN_LABEL = Dict(
    "h5.33" => "Pyrolysis.jl, h = 5.33 W/(m²·K)",
    "tablesC" => "Pyrolysis.jl, tables read in °C",
    "T303" => "Pyrolysis.jl, 303.15 K start",
)
const VARIANTS = ["base", "collimated", "fds_closure", "fds_matched", "surface"]
const FDS_CELL = 5e-5

const CASES = [(q, L) for (q, L, _) in GASIFICATION_CASES]
casename(q, L) = "gasification_$(q)kW_$(L)mm"

# ---------------------------------------------------------------------------
# Numbers
# ---------------------------------------------------------------------------

function load_source(id, src, case)
    code, variant = split(src, "/")
    return load_run(code, id, variant, case)
end

rel(a, b) = (isfinite(a) && isfinite(b) && b != 0) ? (a - b) / b : NaN
span(v) = (w = filter(isfinite, v); "$(pct(minimum(w)))–$(pct(maximum(w)))")
howmany(k, n) = k == n ? "all $n" : "$k of $n"
medabs(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(abs.(w)))
maxabs(v) = (w = filter(isfinite, v); isempty(w) ? NaN : maximum(abs.(w)))
pct(x; d = 1) =
    isfinite(x) ? (d == 1 ? @sprintf("%.1f", 100x) : @sprintf("%.2f", 100x)) : "–"
function fmt(x, f = "%.2f")
    isfinite(x) || return "–"
    t = Printf.format(Printf.Format(f), x)
    return occursin(r"^-0\.?0*$", t) ? t[2:end] : t
end

# Metrics of every source of every deck: (set, case, source) => metrics.
function all_metrics(decks)
    M = Dict{Tuple{String,String,String},Any}()
    for deck in decks, (q, L) in CASES
        id, c = deck["id"], casename(q, L)
        sub = load_submission(deck, q, L)
        sub === nothing || (M[(id, c, "submission")] = metrics(sub))
        for code in ("pyrolysis", "fds", "gpyro", "thermakin")
            dir = joinpath(RESULTS, code, id)
            isdir(dir) || continue
            for v in sort(readdir(dir))
                s = load_run(code, id, v, c)
                s === nothing || (M[(id, c, "$code/$v")] = metrics(s))
            end
        end
    end
    return M
end

has(M, id, src) = any(haskey(M, (id, casename(q, L), src)) for (q, L) in CASES)
function deltas(M, id, a, b, field)
    return [
        haskey(M, (id, casename(q, L), a)) && haskey(M, (id, casename(q, L), b)) ?
        rel(
            getfield(M[(id, casename(q, L), a)], field),
            getfield(M[(id, casename(q, L), b)], field),
        ) : NaN for (q, L) in CASES
    ]
end
d50(M, id, a, b) = deltas(M, id, a, b, :t50)
dpk(M, id, a, b) = deltas(M, id, a, b, :peak_mlr)
cell(M, id, a, b) =
    has(M, id, a) && has(M, id, b) ?
    "$(pct(medabs(d50(M, id, a, b)))) ($(pct(medabs(dpk(M, id, a, b)))))" : "–"

# Series of `a` against `b` at b's saved times: the largest mass difference,
# and the RMS differences of the surface and back temperatures while both hold
# more than 10 % of their initial mass.
function series_diff(a, b; frac = 0.1)
    idx = [j for j in eachindex(b.time_s) if a.time_s[1] <= b.time_s[j] <= a.time_s[end]]
    am = [interp(a.time_s, a.mass_g, b.time_s[j]) for j in idx]
    dm = maximum(abs.(am .- b.mass_g[idx]))
    keep = [
        k for k in eachindex(idx) if
        am[k] > frac * a.mass_g[1] && b.mass_g[idx[k]] > frac * b.mass_g[1]
    ]
    function rms(col)
        v = filter(
            isfinite,
            [
                (
                    interp(a.time_s, getfield(a, col), b.time_s[idx[k]]) -
                    getfield(b, col)[idx[k]]
                )^2 for k in keep
            ],
        )
        return isempty(v) ? NaN : sqrt(sum(v) / length(v))
    end
    return dm, rms(:T_top_K), rms(:T_back_K)
end

# The match of a Pyrolysis.jl variant to another code's run over the six cases.
function match_stats(id, a, b)
    r = []
    for (q, L) in CASES
        x, y = load_source(id, a, casename(q, L)), load_source(id, b, casename(q, L))
        (x === nothing || y === nothing) && continue
        mx, my = metrics(x), metrics(y)
        push!(r, (rel(mx.t50, my.t50), rel(mx.peak_mlr, my.peak_mlr), series_diff(x, y)...))
    end
    isempty(r) && return nothing
    return (
        t50_max = maxabs(getindex.(r, 1)),
        t50_med = medabs(getindex.(r, 1)),
        peak_max = maxabs(getindex.(r, 2)),
        dm = maximum(getindex.(r, 3)),
        Ttop = maxabs(getindex.(r, 4)),
        Tback = maxabs(getindex.(r, 5)),
        n = length(r),
    )
end

kappas(deck) = [Float64(c["kappa"]) for c in deck["components"] if haskey(c, "kappa")]
fds_factor(κ) = log1p(2κ * FDS_CELL) / (1.5κ * FDS_CELL)

# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

function table(io, header, rows)
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", repeat("---|", length(header)))
    for r in rows
        println(io, "| ", join(r, " | "), " |")
    end
    println(io)
end

# A figure on the overview, or (`up`) on a detail page one directory down.
figure(io, file, alt; up = false) =
    println(io, "![$alt]($(up ? "../" : "")$ASSETS/$file)\n")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

function axis(pos; kw...)
    return Axis(
        pos;
        xgridcolor = GRID,
        ygridcolor = GRID,
        topspinevisible = false,
        rightspinevisible = false,
        kw...,
    )
end

function save_figure(name, fig)
    mkpath(FIG)
    save(joinpath(FIG, name), fig; px_per_unit = 1)
    return name
end

const T_DRAWN = (200.0, 1200.0)      # K, temperatures drawn in the gasification figures
const TOP_LABEL = Dict(
    :T_top_K => "T surface (K)",
    :T_avg1mm_K => "T, top 1 mm (K)",
    :T_3pt_K => "T, 0–1 mm (K)",
)

function gasification_figure(deck)
    id = deck["id"]
    fig = Figure(size = (1200, 800), fontsize = 11, backgroundcolor = :white)
    topcol = top_column(deck)
    rows = [
        (:mass_g, "mass (g)"),
        (:mlr_g_m2s, "MLR (g/(m²·s))"),
        (:top, TOP_LABEL[topcol]),
        (:T_back_K, "T back (K)"),
    ]
    sources = copy(SOURCES)
    haskey(AS_RAN, id) && push!(
        sources,
        ("pyrolysis/" * AS_RAN[id], AS_RAN_LABEL[AS_RAN[id]], BLUE, :dashdot),
    )
    used = String[]
    elems = Any[]
    for (j, (q, L)) in enumerate(CASES)
        c = casename(q, L)
        sub = load_submission(deck, q, L)
        for (i, (col, ylabel)) in enumerate(rows)
            ax = axis(
                fig[i+1, j];
                title = i == 1 ? "$q kW/m², $L mm" : "",
                xlabel = i == length(rows) ? "time (s)" : "",
                ylabel = j == 1 ? ylabel : "",
            )
            tmax = 0.0
            if sub !== nothing
                lines!(
                    ax,
                    sub.time_s,
                    getfield(sub, col == :top ? :T_top_K : col);
                    color = INK,
                    linewidth = 1.3,
                )
                tmax = max(tmax, burnout_time(sub))
                "submission" in used || (
                    push!(used, "submission");
                    push!(elems, LineElement(color = INK, linewidth = 2))
                )
            end
            for (src, label, colour, ls) in sources
                s = load_source(id, src, c)
                s === nothing && continue
                y =
                    col == :top ?
                    (startswith(src, "pyrolysis") ? getfield(s, topcol) : s.T_top_K) :
                    getfield(s, col)
                # ThermaKin2Ds' last outputs on UMD, where its integration diverges once
                # the deck's k turns negative, leave the physical range: not drawn.
                col in (:top, :T_back_K) &&
                    (y = [T_DRAWN[1] <= x <= T_DRAWN[2] ? x : NaN for x in y])
                lines!(ax, s.time_s, y; color = colour, linestyle = ls, linewidth = 1.6)
                tmax = max(tmax, burnout_time(s))
                label in used || (
                    push!(used, label);
                    push!(
                        elems,
                        LineElement(color = colour, linestyle = ls, linewidth = 2.5),
                    )
                )
            end
            tmax > 0 && xlims!(ax, 0, 1.1 * tmax)
        end
    end
    labels = [
        l == "submission" ?
        "$(deck["predictions"]) submission ($(deck["calibration_code"]))" : l for
        l in used
    ]
    haskey(deck, "predictions") || (labels = used)
    Legend(
        fig[1, 1:length(CASES)],
        elems,
        labels;
        orientation = :horizontal,
        framevisible = false,
        tellheight = true,
    )
    return save_figure("gasification_$(id).png", fig)
end

function parity_figure(M, ids)
    fig = Figure(size = (1000, 470), fontsize = 12, backgroundcolor = :white)
    used = Any[]
    for (j, (field, name)) in enumerate((
        (:t50, "time to 50 % mass loss (s)"),
        (:peak_mlr, "peak MLR (g/(m²·s))"),
    ))
        ax = axis(
            fig[2, j];
            title = name,
            xlabel = "submission",
            ylabel = j == 1 ? "code's run of the same deck" : "",
            xscale = log10,
            yscale = log10,
            aspect = 1,
        )
        lo, hi = Inf, -Inf
        for (k, (src, label, colour, _)) in enumerate(SOURCES)
            xs, ys = Float64[], Float64[]
            for id in ids, (q, L) in CASES
                a, b = get(M, (id, casename(q, L), src), nothing),
                get(M, (id, casename(q, L), "submission"), nothing)
                (a === nothing || b === nothing) && continue
                x, y = getfield(b, field), getfield(a, field)
                (isfinite(x) && isfinite(y) && x > 0 && y > 0) || continue
                push!(xs, x)
                push!(ys, y)
            end
            isempty(xs) && continue
            mk = (:circle, :rect, :diamond, :utriangle, :dtriangle)[k]
            scatter!(
                ax,
                xs,
                ys;
                color = colour,
                marker = mk,
                markersize = 8,
                strokecolor = :white,
                strokewidth = 0.5,
            )
            j == 1 && push!(
                used,
                (MarkerElement(color = colour, marker = mk, markersize = 10), label),
            )
            lo, hi = min(lo, minimum(xs), minimum(ys)), max(hi, maximum(xs), maximum(ys))
        end
        lines!(ax, [lo, hi], [lo, hi]; color = :gray50, linewidth = 1)
    end
    Legend(
        fig[1, 1:2],
        first.(used),
        last.(used);
        orientation = :horizontal,
        framevisible = false,
    )
    return save_figure("parity.png", fig)
end

# TGA submissions: Sandia files hold time, T [°C], T [K], mass-loss fraction,
# mass [mg]; the rest time, T [K], mass.
function tga_submission(deck, β)
    p = tga_prediction_path(deck, β)
    p === nothing && return nothing
    a = read_matl_csv(p)
    size(a, 1) > 1 || return nothing
    iT, im = size(a, 2) >= 5 ? (3, 5) : (2, 3)
    return (T = a[:, iT], m = a[:, im] ./ a[1, im])
end

# The temperature at which a TGA has lost the fraction `f` of its mass.
function tga_T(T, m, f)
    i = findfirst(x -> x <= 1 - f, m)
    (i === nothing || i == 1) && return NaN
    return T[i-1] + (T[i] - T[i-1]) * (m[i-1] - (1 - f)) / (m[i-1] - m[i])
end

function tga_figure(decks)
    n = length(decks)
    ncol = 5
    fig = Figure(
        size = (1200, 230 * cld(n, ncol) + 60),
        fontsize = 11,
        backgroundcolor = :white,
    )
    for (k, deck) in enumerate(decks)
        id = deck["id"]
        ax = axis(
            fig[cld(k, ncol)+1, mod1(k, ncol)];
            title = id,
            xlabel = "T (K)",
            ylabel = mod1(k, ncol) == 1 ? "m/m₀" : "",
        )
        for (β, ls) in zip(TGA_RATES, (:solid, :dash))
            s = load_run("pyrolysis", id, "base", "tga_$(β)K")
            s === nothing ||
                lines!(ax, s.T_K, s.mass_rel; color = BLUE, linestyle = ls, linewidth = 1.8)
            sub = tga_submission(deck, β)
            sub === nothing ||
                lines!(ax, sub.T, sub.m; color = INK, linestyle = ls, linewidth = 1.1)
        end
        xlims!(ax, 450, 800)
    end
    elems = [
        LineElement(color = BLUE, linewidth = 2.5),
        LineElement(color = INK, linewidth = 2.5),
        LineElement(color = :gray40, linestyle = :solid),
        LineElement(color = :gray40, linestyle = :dash),
    ]
    Legend(
        fig[1, 1:ncol],
        elems,
        ["Pyrolysis.jl", "submission", "10 K/min", "100 K/min"];
        orientation = :horizontal,
        framevisible = false,
    )
    return save_figure("tga.png", fig)
end

function apparatus_figure(deck)
    id = deck["id"]
    fig = Figure(size = (1200, 310), fontsize = 11, backgroundcolor = :white)
    for (j, (key, reps, col, ylabel, title)) in enumerate(APPARATUS_PANELS)
        ax = axis(fig[2, j]; title, xlabel = "time (s)", ylabel)
        tmax = 0.0
        for r in reps
            te, ye = apparatus_experiment(key, r)
            lines!(ax, te, ye; color = EXP, linewidth = 2.5)
            tmax = max(tmax, maximum(te))
        end
        for f in group_predictions(deck, key)
            a = read_matl_csv(f)
            size(a, 1) > 1 && lines!(ax, a[:, 1], a[:, 2]; color = INK, linewidth = 1.1)
        end
        for r in reps
            s = load_run(
                "pyrolysis",
                id,
                "base",
                "apparatus_" * replace(key, "MLR" => "Mass") * "_$r",
            )
            s === nothing ||
                lines!(ax, s.time_s, getfield(s, col); color = BLUE, linewidth = 1.3)
        end
        xlims!(ax, 0, 1.1 * tmax)
        col == :mass_g && ylims!(ax, -2, nothing)
    end
    group =
        haskey(deck, "validation") ? "$(basename(deck["validation"])) prediction" :
        "no group prediction"
    Legend(
        fig[1, 1:4],
        [
            LineElement(color = EXP, linewidth = 3),
            LineElement(color = INK, linewidth = 2),
            LineElement(color = BLUE, linewidth = 2.5),
        ],
        ["experiment (each repeat)", group, "Pyrolysis.jl, $id"];
        orientation = :horizontal,
        framevisible = false,
    )
    return save_figure("apparatus_$(id).png", fig)
end

function inert_figure()
    fig = Figure(size = (1100, 380), fontsize = 12, backgroundcolor = :white)
    for (j, t) in enumerate(INERT_TESTS)
        ax = axis(
            fig[2, j];
            title = t.test,
            xlabel = "time (s)",
            ylabel = j == 1 ? "T (K)" : "",
        )
        a = read_matl_csv(inert_experiment_path(t.test))
        s = load_run("pyrolysis", "inert", "", t.test; root = RESULTS)
        cols = s === nothing ? Symbol[] : collect(keys(s))[3:end]
        for (k, d) in enumerate(t.depths)
            lines!(ax, a[:, 1], a[:, 1+k]; color = EXP, linewidth = 2.5)
            s === nothing ||
                lines!(ax, s.time_s, getfield(s, cols[k]); color = BLUE, linewidth = 1.6)
            text!(
                ax,
                a[end, 1],
                a[end, 1+k];
                text = @sprintf(" %.2f mm", 1e3d),
                align = (:left, :center),
                fontsize = 11,
            )
        end
        xlims!(ax, 0, 1.18 * t.t_end)
    end
    Legend(
        fig[1, 1:2],
        [
            LineElement(color = EXP, linewidth = 3),
            LineElement(color = BLUE, linewidth = 2.5),
        ],
        ["experiment (mean of three)", "Pyrolysis.jl, shared apparatus deck"];
        orientation = :horizontal,
        framevisible = false,
    )
    return save_figure("inert.png", fig)
end

const COST_COLOUR = Dict(
    "Pyrolysis.jl" => BLUE,
    "fds" => colorant"#eb6834",
    "gpyro" => colorant"#1baf7a",
    "thermakin" => colorant"#e87ba4",
)
const COST_NAME = Dict(
    "Pyrolysis.jl" => "Pyrolysis.jl",
    "fds" => "FDS 6.11.1",
    "gpyro" => "Gpyro V2",
    "thermakin" => "ThermaKin2Ds",
)

function cost_figure(rows)
    keys_ = sort(unique((r[2], r[3]) for r in rows if r[1] != "Pyrolysis.jl"))
    n = length(keys_)
    ncol = 4
    SHAPE = Dict("10" => :circle, "20" => :rect, "40" => :diamond)
    marker(r) =
        r[1] == "Pyrolysis.jl" ? SHAPE[first(split(r[4], "/"))] :
        r[4] == "base" ? :circle : :rect
    filled(r) = r[1] != "Pyrolysis.jl" || endswith(r[4], "×10")
    fig = Figure(
        size = (1300, 290 * cld(n, ncol) + 90),
        fontsize = 12,
        backgroundcolor = :white,
    )
    for (k, key) in enumerate(keys_)
        ax = axis(
            fig[cld(k, ncol)+1, mod1(k, ncol)];
            title = "$(key[1]), $(replace(key[2], "gasification_" => ""))",
            xscale = log10,
            yscale = log10,
            limits = ((5e-4, 30), (0.05, 100)),
            xlabel = "error / tolerance",
            ylabel = mod1(k, ncol) == 1 ? "solve time (s)" : "",
        )
        vlines!(ax, [1.0]; color = :gray50, linestyle = :dash)
        for code in ("Pyrolysis.jl", "fds", "gpyro", "thermakin"),
            r in filter(r -> (r[2], r[3]) == key && r[1] == code && isfinite(r[6]), rows)

            c = COST_COLOUR[code]
            scatter!(
                ax,
                [max(r[6], 1e-3)],
                [r[5]];
                marker = marker(r),
                markersize = 10,
                color = filled(r) ? c : :white,
                strokecolor = c,
                strokewidth = 1.8,
            )
        end
    end
    Legend(
        fig[1, 1:ncol],
        [
            [
                MarkerElement(color = COST_COLOUR[c], marker = :circle, markersize = 10) for
                c in ("Pyrolysis.jl", "fds", "gpyro", "thermakin")
            ];
            [
                MarkerElement(color = :gray30, marker = SHAPE[s], markersize = 10) for
                s in ("10", "20", "40")
            ]
        ],
        [
            [COST_NAME[c] for c in ("Pyrolysis.jl", "fds", "gpyro", "thermakin")];
            [
                "10 cells/mm (external: base)",
                "20 cells/mm (external: refined)",
                "40 cells/mm",
            ]
        ];
        orientation = :horizontal,
        framevisible = false,
    )
    return save_figure("cost.png", fig)
end


# ---------------------------------------------------------------------------
# The pages
# ---------------------------------------------------------------------------

# The overview's detail pages, in navigation order.
const DETAIL = joinpath(dirname(PAGE), "macfp_pmma")
const DECKS_PAGE = "01_decks.md"
const GAS_PAGE = "02_gasification.md"
const HIST_PAGE = "03_histories.md"
const TGA_PAGE = "04_tga.md"
const APP_PAGE = "05_apparatus.md"
const COST_PAGE = "06_cost.md"
const OVERVIEW = "[MaCFP PMMA validation](../macfp_pmma.md)"

# A link from the overview to a detail page, or to one of its sections
# (Documenter's anchor: the heading with hyphens for spaces).
function link(file, text, section = "")
    anchor = isempty(section) ? "" : "#" * replace(section, " " => "-")
    return "[$text](macfp_pmma/$file$anchor)"
end

function write_page(path, io)
    mkpath(dirname(path))
    text = replace(String(take!(io)), r"\n{3,}" => "\n\n")
    write(path, rstrip(text) * "\n")
    return path
end

# A percentage as the pages print it, so that comparisons agree with the tables.
shown(x) = round(100x; digits = 1)

function main()
    check_matl_db()
    rm(FIG; force = true, recursive = true)
    rm(DETAIL; force = true, recursive = true)
    decks = load_decks()
    byid = Dict(d["id"] => d for d in decks)
    gas = filter(d -> !get(d, "kinetics_only", false), decks)
    withsub = filter(d -> haskey(d, "predictions"), gas)
    indepth = filter(d -> d["absorption"] == "in-depth", gas)
    M = all_metrics(gas)
    headline = filter(d -> d["id"] != "BUW-FZJ_A", withsub)

    # =========================================================================
    # Decks and scenario
    # =========================================================================
    io = IOBuffer()
    println(
        io,
        """
# Decks and scenario

The property sets, the scenario and the variants behind the $OVERVIEW.

Data come from [MaCFP/matl-db](https://github.com/MaCFP/matl-db) at commit
`$(MacfpPMMA.MATL_DB_COMMIT[1:8])`. The decks, the runners and the generator of these pages are
in `scripts/macfp_pmma/` (its `README.md` describes each choice); the translation of each
deck for FDS, Gpyro and ThermaKin2Ds is described beside their runners in
`scripts/macfp_pmma/external/`.

## Scenario

The MaCFP-2 gasification targets expose a 6 or 12 mm slab at 293 K to 10, 25 or
65 kW/m² with no convection; the slab re-radiates to 293 K and its back is perfectly
insulated. The flux is absorbed at the surface emissivity. Pyrolysis.jl runs each case
on a moving mesh of 20 cells/mm with deactivation, at the default tolerances. A deck
that absorbs at its surface runs under `SURFACE_ABSORPTION`. A deck with an absorption
coefficient κ runs under `RADIATIVE_TRANSFER`: every cell absorbs and emits, the
heater is `Diffuse()`, and the insulated back is a mirror. The TGA targets heat one
50 µm cell from 300 K at 10 and 100 K/min.

## Property sets
""",
    )
    rows = []
    for d in gas
        κ = kappas(d)
        push!(
            rows,
            [
                d["id"],
                string(d["year"]),
                d["calibration_code"],
                d["absorption"],
                d["absorption"] == "in-depth" ?
                join(unique(@sprintf("%.0f", k) for k in κ), ", ") : "–",
                string(length(d["reactions"])),
                haskey(d, "predictions") ? "yes" : "no",
            ],
        )
    end
    table(
        io,
        [
            "set",
            "year",
            "calibrated in",
            "absorption",
            "κ (1/m)",
            "reactions",
            "MaCFP-2 submission",
        ],
        rows,
    )
    kin = filter(d -> get(d, "kinetics_only", false), decks)
    println(
        io,
        "The sets ",
        join([d["id"] for d in kin], ", "),
        " give kinetics only and run the TGA targets alone.\n",
    )
    ran(code) = [d["id"] for d in gas if isdir(joinpath(RESULTS, code, d["id"], "base"))]
    println(
        io,
        """
FDS 6.11.1 runs $(length(ran("fds")) == length(gas) ? "all $(length(gas)) sets" : join(ran("fds"), ", ")),
Gpyro V2 runs $(join(ran("gpyro"), ", ")), and ThermaKin2Ds 7 runs $(join(ran("thermakin"), ", ")) in its
`MAX` mode (the absorbed radiation deposited in the element that absorbs most) and
$(join([d["id"] for d in gas if isdir(joinpath(RESULTS, "thermakin", d["id"], "rand"))], ", ")) also in its
`RAND` mode (in an element drawn from the absorption distribution).

## Variants

Each in-depth deck also runs variants that move one choice, and three variants run
what a group evidently ran:

| variant | what it changes |
|---|---|
| `collimated` | the heater as one beam along the normal, attenuated at κ rather than the diffuse field's 1.5κ |
| `fds_closure` | every component's κ × ln(1 + 2κΔz)/(1.5κΔz) at Δz = 50 µm: FDS's two-flux attenuation on its own cells |
| `fds_matched` | `fds_closure` with the one item FDS cannot express aligned with FDS: no pyrolyzate enthalpy in the solid (UMD, UMET_TK); UMET_GP's second reaction first order, as FDS runs it on a shrinking cell |
| `surface` | `SURFACE_ABSORPTION`: no radiation in depth |
| `h5.33` | BUW-FZJ_A/B: convection h = 5.33 W/(m²·K) to 293 K, which their submissions include |
| `tablesC` | DBI_1, DBI_calibrated: the kelvin-valued c and k tables read as °C, as DBI_1's FDS runs read them (DBI_calibrated's are not published) |
| `T303` | Aalto_I: start and surroundings at 303.15 K, where their runs start |
""",
    )
    write_page(joinpath(DETAIL, DECKS_PAGE), io)

    # =========================================================================
    # Gasification
    # =========================================================================
    io = IOBuffer()
    println(
        io,
        """
# Gasification

The MaCFP-2 gasification targets of the $OVERVIEW: each code's run of every deck
against the institution's submission, against the other codes' runs, and against FDS
with the attenuation matched. t₅₀ is the first saved time at or below half the initial
mass; the peak is the largest finite mass-loss rate. Unless a table says otherwise, a
cell is the median over the six cases of the relative difference in t₅₀, |Δt₅₀|, and in
brackets in the peak, |Δpeak| (%).

## Against the submissions

BUW-FZJ_A is reported on its own below.
""",
    )
    rows = []
    for d in headline
        id = d["id"]
        r = [id, d["calibration_code"]]
        for (src, _, _, _) in SOURCES
            push!(r, cell(M, id, src, "submission"))
        end
        push!(rows, r)
    end
    table(io, vcat(["set", "calibrated in"], [s[2] for s in SOURCES]), rows)
    println(
        io,
        """
Where a set was calibrated in a code, that code's run of the deck checks the deck: Gpyro
reproduces GIDAZE+ to $(pct(medabs(d50(M, "GIDAZE+", "gpyro/base", "submission")))) %, FDS reproduces NIST to
$(pct(medabs(d50(M, "NIST", "fds/base", "submission")))) % and Aalto_I to $(pct(medabs(d50(M, "Aalto_I", "fds/base", "submission")))) %, and ThermaKin2Ds reproduces UMD to
$(span(abs.(d50(M, "UMD", "thermakin/base", "submission")[1:2]))) % at 10 kW/m² in `MAX` mode and to
$(span(abs.(d50(M, "UMD", "thermakin/rand", "submission")[3:6]))) % at 25 and 65 kW/m² in `RAND` mode. Gpyro
reproduces UMET_GP to $(pct(medabs(d50(M, "UMET_GP", "gpyro/base", "submission")))) % with the density held at its 293 K value, as every
deck holds it, and to $(pct(medabs(d50(M, "UMET_GP", "gpyro/rhoT", "submission")))) % with UMET's temperature-dependent density.
ThermaKin2Ds does not reproduce UMET_TK ($(pct(medabs(d50(M, "UMET_TK", "thermakin/base", "submission")))) % in `MAX`,
$(pct(medabs(d50(M, "UMET_TK", "thermakin/rand", "submission")))) % in `RAND`, $(pct(medabs(d50(M, "UMET_TK", "thermakin/rhoT", "submission")))) % in `MAX` with UMET's ρ(T)).
""",
    )
    # The calibration code's own run: UMD's in ThermaKin2Ds is MAX at 10 kW/m² and
    # RAND at 25 and 65 kW/m², the reading that reproduces its submission; UMET_TK's
    # mode is not known.
    function own(id)
        id == "UMD" && return vcat(
            d50(M, id, "thermakin/base", "submission")[1:2],
            d50(M, id, "thermakin/rand", "submission")[3:6],
        )
        code = byid[id]["calibration_code"]
        code == "FDS" && return d50(M, id, "fds/base", "submission")
        code == "Gpyro" && return d50(M, id, "gpyro/base", "submission")
        return nothing
    end
    ind_sub = filter(d -> d["absorption"] == "in-depth", headline)
    ours = [
        (d["id"], medabs(d50(M, d["id"], "pyrolysis/base", "submission"))) for d in ind_sub
    ]
    known = filter(d -> own(d["id"]) !== nothing, ind_sub)
    closer = [
        d["id"] for d in known if
        medabs(d50(M, d["id"], "pyrolysis/base", "submission")) < medabs(own(d["id"]))
    ]
    farther = setdiff([d["id"] for d in known], closer)
    tk = [medabs(d50(M, "UMET_TK", "thermakin/$m", "submission")) for m in ("base", "rand")]
    ptk = medabs(d50(M, "UMET_TK", "pyrolysis/base", "submission"))
    println(
        io,
        """
On the in-depth sets Pyrolysis.jl's median |Δt₅₀| to the submission is
$(join(["$(i) $(pct(v))" for (i, v) in ours], ", ")) %. It is closer to the submission than the
calibration code's own run on $(isempty(closer) ? "none" : join(closer, ", ")) (UMD's ThermaKin2Ds run taken in `MAX` at
10 kW/m² and `RAND` above) and farther on $(isempty(farther) ? "none" : join(farther, ", ")). UMET_TK's ThermaKin2Ds mode is
not known; Pyrolysis.jl ($(pct(ptk)) %) is $(minimum(tk) <= ptk <= maximum(tk) ? "between" : "outside") `MAX` ($(pct(tk[1])) %) and `RAND`
($(pct(tk[2])) %) there.
""",
    )
    parity = parity_figure(M, [d["id"] for d in headline])
    figure(io, parity, "Parity of t50 and peak MLR against the submissions"; up = true)

    println(
        io,
        """
## Code to code

Pyrolysis.jl's run of each deck against each other code's run of the same deck. The
2023 sets NIST-StMU, DBI_calibrated and BUW-FZJ_C have no MaCFP-2 submission; FDS runs
them.
""",
    )
    rows = []
    for d in filter(d -> d["id"] != "BUW-FZJ_A", gas)
        id = d["id"]
        push!(
            rows,
            vcat(
                [id, d["absorption"]],
                [cell(M, id, "pyrolysis/base", s[1]) for s in SOURCES[2:end]],
            ),
        )
    end
    table(io, vcat(["set", "absorption"], [s[2] for s in SOURCES[2:end]]), rows)
    surf = filter(d -> d["absorption"] == "surface" && d["id"] != "BUW-FZJ_A", gas)
    sfds = [
        medabs(d50(M, d["id"], "pyrolysis/base", "fds/base")) for
        d in surf if has(M, d["id"], "fds/base")
    ]
    println(
        io,
        """
On the decks that absorb at the surface, Pyrolysis.jl and FDS solve the same problem and
agree to $(pct(minimum(sfds)))–$(pct(maximum(sfds))) % in t₅₀. On the in-depth decks the codes differ in how
the solid absorbs and emits the radiation (Technical Reference §8.9).

## Matched to FDS

FDS's in-depth model is a two-flux model with emission, like Pyrolysis.jl's, but its
sweep crosses a cell of width Δz with transmission 1/(1 + 2κΔz), where the diffuse
two-stream model transmits e^(−1.5κΔz). At an absorption coefficient of
ln(1 + 2κΔz)/(1.5Δz) the two agree cell by cell, in transmission and in emission
(`fds_closure`, taken at FDS's 50 µm cells). Where a deck holds an item FDS cannot
express, `fds_matched` also aligns it. FDS 6.11.1's runs of the same decks:
""",
    )
    rows = []
    matched = Dict{String,Any}()               # set => match_stats against FDS
    within = Dict{String,Tuple{Bool,Bool}}()   # t50, peak within FDS's own 50 → 25 µm change
    for d in indepth
        id = d["id"]
        src =
            has(M, id, "pyrolysis/fds_matched") ? "pyrolysis/fds_matched" :
            "pyrolysis/fds_closure"
        st = match_stats(id, src, "fds/base")
        st === nothing && continue
        matched[id] = st
        fdsref = match_stats(id, "fds/base_c2.5e-5", "fds/base")
        fdsref === nothing ||
            (within[id] = (st.t50_max <= fdsref.t50_max, st.peak_max <= fdsref.peak_max))
        f = fds_factor.(kappas(d))
        item =
            id in ("UMD", "UMET_TK") ? "no gas enthalpy" :
            id == "UMET_GP" ? "reaction 2 first order" : "–"
        push!(
            rows,
            [
                id,
                length(unique(f)) == 1 ? @sprintf("%.3f", f[1]) :
                join(unique(@sprintf("%.3f", x) for x in f), ", "),
                item,
                pct(st.t50_max; d = 2),
                pct(st.t50_med; d = 2),
                pct(st.peak_max; d = 2),
                "$(fmt(st.Ttop)) / $(fmt(st.Tback))",
                fmt(st.dm),
                fdsref === nothing ? "–" :
                "$(pct(fdsref.t50_max; d = 2)) / $(pct(fdsref.peak_max; d = 2))",
            ],
        )
    end
    table(
        io,
        [
            "set",
            "factor on κ",
            "aligned item",
            "max \\|Δt₅₀\\| (%)",
            "median \\|Δt₅₀\\| (%)",
            "max \\|Δpeak\\| (%)",
            "worst RMS ΔT top / back (K)",
            "max \\|Δm\\| (g)",
            "FDS 25 µm vs 50 µm: max \\|Δt₅₀\\| / \\|Δpeak\\| (%)",
        ],
        rows,
    )
    base_fds =
        Dict(d["id"] => match_stats(d["id"], "pyrolysis/base", "fds/base") for d in indepth)
    filter!(p -> p[2] !== nothing, base_fds)
    bf = collect(values(base_fds))
    println(
        io,
        """
Temperatures are compared at FDS's saved times while both runs hold more than 10 % of
their mass; masses are for the 10 cm × 10 cm sample. With the deck's κ as written and
the diffuse heater the same runs differ from FDS by $(pct(minimum(s.t50_med for s in bf)))–$(pct(maximum(s.t50_med for s in bf))) % in median t₅₀ and
up to $(pct(maximum(s.peak_max for s in bf))) % in peak. FDS also runs $(join(sort(collect(keys(within))), ", ")) on
25 µm cells. Matched, the largest |Δt₅₀| is within FDS's own change between its 50 and
25 µm runs on $(howmany(count(first, values(within)), length(within))) of these sets$(let w = sort([k for (k, v) in within if !v[1]]); isempty(w) ? "" : " (not on $(join(w, ", ")))" end),
and the largest |Δpeak| on $(howmany(count(last, values(within)), length(within)))$(let w = sort([k for (k, v) in within if !v[2]]); isempty(w) ? "" : " (not on $(join(w, ", ")))" end). The factor
belongs to FDS's cell size, so it changes with the cell FDS runs. BUW-FZJ_C is the
exception to a single factor: its residue's κ = 50 000 1/m takes a factor of $(@sprintf("%.2f", fds_factor(5e4))),
and a factor applied per component is not FDS's attenuation of the cell's mixture.

## Radiation in depth

Each Pyrolysis.jl variant of the in-depth decks against the submission and each code.
Bold: the variant closest to that column's run in t₅₀.
""",
    )
    refs = [("submission", "submission"); [(s[1], s[2]) for s in SOURCES[2:end]]]
    for d in indepth
        id = d["id"]
        cols = filter(r -> has(M, id, r[1]), refs)
        vs = filter(v -> has(M, id, "pyrolysis/$v"), VARIANTS)
        best = Dict(
            r[1] => argmin([
                let x = medabs(d50(M, id, "pyrolysis/$v", r[1]))
                    isfinite(x) ? x : Inf
                end for v in vs
            ]) for r in cols
        )
        rows = [
            [
                "`$v`";
                [
                    (best[r[1]] == k ? "**" : "") *
                    cell(M, id, "pyrolysis/$v", r[1]) *
                    (best[r[1]] == k ? "**" : "") for r in cols
                ]
            ] for (k, v) in enumerate(vs)
        ]
        println(io, "**$id** (calibrated in $(d["calibration_code"]))\n")
        table(io, ["variant"; [r[2] for r in cols]], rows)
    end
    g(id, v, r) = pct(medabs(d50(M, id, "pyrolysis/$v", r)))
    gv(id, v, r) = medabs(d50(M, id, "pyrolysis/$v", r))
    ids_with(r) = [d["id"] for d in indepth if has(M, d["id"], r)]
    # The in-depth variant closest to a run in median |Δt50|.
    closest(id, r, vs) = vs[argmin([
        let x = gv(id, v, r)
            isfinite(x) ? x : Inf
        end for v in vs
    ])]
    calib_sub = [
        d["id"] for d in indepth if
        d["calibration_code"] in ("Gpyro", "ThermaKin") && has(M, d["id"], "submission")
    ]
    heater_sub = [
        id for id in calib_sub if
        gv(id, "base", "submission") < gv(id, "collimated", "submission")
    ]
    heater_fds = [
        id for id in ids_with("fds/base") if
        gv(id, "base", "fds/base") < gv(id, "collimated", "fds/base")
    ]
    not_fds = setdiff(ids_with("fds/base"), heater_fds)
    gp = ids_with("gpyro/base")
    gp_close = [
        id for
        id in gp if closest(id, "gpyro/base", ["base", "collimated", "fds_closure"]) ==
        "fds_closure"
    ]
    tkr = ids_with("thermakin/base")
    tk_surf = [
        id for id in tkr if closest(
            id,
            "thermakin/base",
            filter(v -> has(M, id, "pyrolysis/$v"), VARIANTS),
        ) == "surface"
    ]
    rnd = ids_with("thermakin/rand")
    rv = [
        gv(id, v, "thermakin/rand") for id in rnd for
        v in ("base", "collimated", "fds_closure")
    ]
    surface_move = [gv(d["id"], "surface", "pyrolysis/base") for d in indepth]
    println(
        io,
        """
- **The heater.** The diffuse heater deposits at 1.5κ, the collimated one at κ. Diffuse
  is closer than collimated to the submission on $(length(heater_sub)) of the $(length(calib_sub)) sets
  calibrated in Gpyro or ThermaKin2Ds ($(join(["$id $(g(id, "base", "submission")) against $(g(id, "collimated", "submission")) %" for id in calib_sub], ", ")))
  and to FDS on $(length(heater_fds)) of $(length(ids_with("fds/base")))$(isempty(not_fds) ? "" : " (not on $(join(not_fds, ", ")))").
- **Gpyro** deposits the net surface flux along a Beer–Lambert march at κ and emits
  nothing from the depth (Technical Reference §8.9.2). Of the variants that absorb in
  depth, `fds_closure`, which deposits shallowest, is closest to Gpyro on
  $(length(gp_close)) of the $(length(gp)) sets Gpyro runs ($(join(["$id $(g(id, "fds_closure", "gpyro/base")) %" for id in gp], ", "))); no
  variant shares its net-flux deposition.
- **ThermaKin2Ds** deposits each step's absorbed radiation in one element. In `MAX` mode
  that is the surface element for a uniform κ, and `surface` is the closest variant on
  $(length(tk_surf)) of the $(length(tkr)) in-depth sets it runs ($(join(["$id $(g(id, "surface", "thermakin/base")) %" for id in tkr], ", "))). `RAND` draws the
  element from the absorption distribution and sits within $(pct(minimum(filter(isfinite, rv))))–$(pct(maximum(filter(isfinite, rv)))) % of
  the variants that absorb in depth.
- **UMET's runs** use a temperature-dependent density, which the decks hold at its
  293 K value; Gpyro with UMET's ρ(T) reproduces UMET_GP's submission (above).
- **No radiation in depth** (`surface`) moves t₅₀ by
  $(pct(minimum(surface_move)))–$(pct(maximum(surface_move))) % (median per set) from the diffuse base.

## Runs as the groups ran them

Three submissions depart from the MaCFP-2 scenario in a way their files show, and a
Pyrolysis.jl variant reproduces each; on NIST, a surface-absorbing set calibrated in
FDS, the two codes solve the same problem. |Δt₅₀| over the six cases (%):
""",
    )
    rows = []
    for (id, v, ref, what) in (
        ("BUW-FZJ_B", "h5.33", "submission", "convection h = 5.33 W/(m²·K)"),
        ("DBI_1", "tablesC", "submission", "c and k tables read as °C"),
        ("Aalto_I", "T303", "submission", "303.15 K start and surroundings"),
        ("NIST", "base", "fds/base", "the deck as written, against FDS"),
    )
        dv = d50(M, id, "pyrolysis/$v", ref)
        db = d50(M, id, "pyrolysis/base", ref)
        push!(
            rows,
            [
                id,
                "`$v`",
                what,
                ref == "submission" ? "submission" : "FDS 6.11.1",
                pct(medabs(dv); d = 2),
                pct(maxabs(dv); d = 2),
                v == "base" ? "–" : pct(medabs(db); d = 2),
            ],
        )
    end
    table(
        io,
        [
            "set",
            "variant",
            "what it reproduces",
            "against",
            "median",
            "max",
            "median, deck as written",
        ],
        rows,
    )

    da = byid["BUW-FZJ_A"]
    res = only(
        filter(
            c -> any(r -> get(r, "residue", "") == c["id"], da["reactions"]),
            da["components"],
        ),
    )
    vir = da["components"][1]
    println(
        io,
        """
## BUW-FZJ_A

BUW-FZJ_A absorbs at the surface, and its residue (ε = $(res["emissivity"]) against the virgin
material's $(vir["emissivity"])) collects in the surface cell as a film. The surface emissivity,
and with it the absorbed flux, then depends on how much of that film the surface cell
holds: in FDS on the size of the first cell, in Pyrolysis.jl on the cell size and on
how long the moving mesh keeps a burnt-out layer in the surface cell before it retires
it. The deck's mass-loss histories measure that discretization as much as the set, so
it is kept out of the medians above. Median / max |Δt₅₀| over the six cases (%):
""",
    )
    rows = []
    for (src, label, ref) in (
        ("pyrolysis/base", "Pyrolysis.jl, 20 cells/mm", "submission"),
        ("pyrolysis/base_c40", "Pyrolysis.jl, 40 cells/mm", "submission"),
        ("pyrolysis/h5.33", "Pyrolysis.jl, h = 5.33, 20 cells/mm", "submission"),
        ("pyrolysis/h5.33_c40", "Pyrolysis.jl, h = 5.33, 40 cells/mm", "submission"),
        ("fds/base", "FDS, 50 µm", "submission"),
        ("fds/base_c2.5e-5", "FDS, 25 µm", "submission"),
        ("gpyro/base", "Gpyro V2", "submission"),
        ("thermakin/base", "ThermaKin2Ds MAX", "submission"),
        ("pyrolysis/base_c40", "Pyrolysis.jl, 40 cells/mm", "pyrolysis/base"),
        ("fds/base_c2.5e-5", "FDS, 25 µm", "fds/base"),
        ("pyrolysis/base", "Pyrolysis.jl, 20 cells/mm", "fds/base"),
    )
        has(M, "BUW-FZJ_A", src) && has(M, "BUW-FZJ_A", ref) || continue
        v = d50(M, "BUW-FZJ_A", src, ref)
        against = Dict(
            "submission" => "submission",
            "pyrolysis/base" => "Pyrolysis.jl, 20 cells/mm",
            "fds/base" => "FDS, 50 µm",
        )[ref]
        push!(
            rows,
            [
                label,
                against,
                pct(medabs(v)),
                pct(maxabs(v)),
                pct(medabs(dpk(M, "BUW-FZJ_A", src, ref))),
            ],
        )
    end
    table(
        io,
        ["run", "against", "median \\|Δt₅₀\\|", "max \\|Δt₅₀\\|", "median \\|Δpeak\\|"],
        rows,
    )
    println(
        io,
        "The mass, mass-loss rate and temperature histories of every set are on ",
        "[Gasification histories]($HIST_PAGE).\n",
    )
    write_page(joinpath(DETAIL, GAS_PAGE), io)

    # =========================================================================
    # Gasification histories
    # =========================================================================
    io = IOBuffer()
    println(
        io,
        """
# Gasification histories

The MaCFP-2 gasification histories of the $OVERVIEW, per set: mass, mass-loss rate, the
submission's top temperature (the surface, or the mean over the top millimetre where the
group reported that) and the back temperature. The interval holding a deactivation event
is blanked in Pyrolysis.jl's MLR: the event vents the last 1 % of a layer's mass in one
step. ThermaKin2Ds' UMD runs end where its integration diverges, once the deck's
conductivity turns negative above 809.5 K; their last temperatures, outside
$(Int(T_DRAWN[1]))–$(Int(T_DRAWN[2])) K, are not drawn.
""",
    )
    for d in vcat(withsub, filter(d -> !haskey(d, "predictions"), gas))
        println(io, "## $(d["id"])\n")
        figure(io, gasification_figure(d), "Gasification histories, $(d["id"])"; up = true)
    end
    write_page(joinpath(DETAIL, HIST_PAGE), io)

    # =========================================================================
    # TGA
    # =========================================================================
    tgadecks = filter(d -> haskey(d, "tga_predictions"), decks)
    io = IOBuffer()
    println(
        io,
        """
# TGA

The MaCFP-2 TGA targets of the $OVERVIEW: the temperature at which 10, 50 and 90 % of
the mass is lost, Pyrolysis.jl minus the submission (K).
""",
    )
    rows = []
    tga_dev = Dict{String,Tuple{Float64,Int}}()    # set => largest |ΔT| (K) and its rate
    for d in tgadecks, β in TGA_RATES
        s = load_run("pyrolysis", d["id"], "base", "tga_$(β)K")
        sub = tga_submission(d, β)
        (s === nothing || sub === nothing) && continue
        dT = [tga_T(s.T_K, s.mass_rel, f) - tga_T(sub.T, sub.m, f) for f in (0.1, 0.5, 0.9)]
        m = maximum(abs, filter(isfinite, dT))
        m > get(tga_dev, d["id"], (-Inf, 0))[1] && (tga_dev[d["id"]] = (m, β))
        push!(rows, vcat([d["id"], string(β)], [fmt(x, "%.1f") for x in dT]))
    end
    table(io, ["set", "β (K/min)", "ΔT₁₀", "ΔT₅₀", "ΔT₉₀"], rows)
    tga_png = tga_figure(tgadecks)
    figure(io, tga_png, "TGA, every set"; up = true)
    write_page(joinpath(DETAIL, TGA_PAGE), io)

    # =========================================================================
    # Apparatus
    # =========================================================================
    θ = rad2deg(APPARATUS_HEATER.θ₂)
    io = IOBuffer()
    println(
        io,
        """
# Gasification apparatus

The MaCFP-3 experiments of the $OVERVIEW: every deck with complete properties runs the
eleven PMMA tests of the NIST gasification apparatus, scored against the measurements.

## Apparatus deck

The tests run on one shared deck: the sample at the thickness that holds its measured
mass at the deck's density on the 6.99 cm disc, bonded to 28.52 mm of Kaowool PM;
convection with h = $(MacfpPMMA.APPARATUS.h) W/(m²·K) to gas at $(MacfpPMMA.APPARATUS.T_gas) K and re-radiation to walls at
$(MacfpPMMA.APPARATUS.T_wall) K; adiabatic beneath the board. The heater is a cone of uniform intensity
from the normal to $(@sprintf("%.0f", θ))°: a coaxial disc fitted to the apparatus's radial flux profile
at 50 kW/m² is seen out to 33–62° within the profile's error bars, and no disc fits the
25 kW/m² profile within its bars. The Kaowool cell under the sample is opaque and absorbs
what the sample transmits; the board beneath it is transparent, because its measured
conductivity already carries its own radiation.

## Scores

The score is MaCFP's: the mean of the runs of the repeats, against the mean of the
measured repeats over the shortest repeat, as NRMS = RMS/(max − min) of the measured
mean; lower is better. For an in-depth deck, the range is that of the diffuse and
collimated heaters, which bound the cone.
""",
    )
    appdecks = filter(
        d -> isfile(
            joinpath(RESULTS, "pyrolysis", d["id"], "base", "apparatus_q50_Mass_R3.csv"),
        ),
        gas,
    )
    S = Dict{Tuple{String,String,Any},Float64}()     # (set, quantity, variant or nothing) => NRMS
    for d in appdecks,
        p in APPARATUS_PANELS,
        v in ("base", "diffuse", "collimated", nothing)

        S[(d["id"], p[1], v)] = apparatus_nrms(d, p; variant = v)
    end
    rows = []
    for d in appdecks
        id = d["id"]
        r = [id]
        for p in APPARATUS_PANELS
            v = S[(id, p[1], "base")]
            if d["absorption"] == "in-depth"
                b = [S[(id, p[1], h)] for h in ("diffuse", "collimated")]
                push!(
                    r,
                    "$(fmt(v, "%.3f")) ($(fmt(minimum(b), "%.3f"))–$(fmt(maximum(b), "%.3f")))",
                )
            else
                push!(r, fmt(v, "%.3f"))
            end
            push!(r, fmt(S[(id, p[1], nothing)], "%.3f"))
        end
        push!(rows, r)
    end
    qname = Dict(
        "q50_Mass" => "mass, 50",
        "q50_MLR" => "MLR, 50",
        "q50_Temp" => "T back, 50",
        "q25_Temp" => "T back, 25",
    )
    table(
        io,
        vcat(
            ["set"],
            vcat([["$(qname[p[1]]): Pyrolysis.jl", "group"] for p in APPARATUS_PANELS]...),
        ),
        rows,
    )
    # Per quantity: the sets with a group prediction, and those where Pyrolysis.jl's
    # score is as good or better.
    cmp = [
        (
            p[1],
            [d["id"] for d in appdecks if isfinite(S[(d["id"], p[1], nothing)])],
            [
                d["id"] for d in appdecks if isfinite(S[(d["id"], p[1], nothing)]) &&
                    S[(d["id"], p[1], "base")] <= S[(d["id"], p[1], nothing)]
            ],
        ) for p in APPARATUS_PANELS
    ]
    spread = maximum(
        abs(S[(d["id"], p[1], "diffuse")] - S[(d["id"], p[1], "collimated")]) for
        d in appdecks if d["absorption"] == "in-depth" for p in APPARATUS_PANELS
    )
    ins = load_run("pyrolysis", "inert", "", "Black-Insulation")
    ains = read_matl_csv(inert_experiment_path("Black-Insulation"))
    deep =
        ins === nothing ? NaN :
        ins.T_17_16mm_K[findfirst(==(1200.0), ins.time_s)] -
        ains[findfirst(==(1200.0), ains[:, 1]), 4]
    println(
        io,
        """
Against the groups' own predictions of the same tests, scored against the same measurements,
Pyrolysis.jl's score is as good or better on $(join(["$(length(c[3])) of $(length(c[2])) sets for $(qname[c[1]]) kW/m²" for c in cmp], ", ")).
The heater's distribution moves an in-depth deck's score by at most $(fmt(spread, "%.3f")). The
back temperatures are those of the sample's interface with the board, which the shared
deck's convection coefficient and the board's conductivity set as much as the
property set: on the inert Kaowool stack (below), the shared deck is
$(fmt(deep, "%.0f")) K above the measurement at 17.16 mm after 1200 s.

## Runs against the experiments
""",
    )
    for d in appdecks
        println(io, "### $(d["id"])\n")
        figure(
            io,
            apparatus_figure(d),
            "NIST gasification apparatus, $(d["id"])";
            up = true,
        )
    end
    println(
        io,
        """
## Inert tests

The apparatus's inert tests check the shared deck without a pyrolysing sample: five
Kaowool discs painted black (ε 0.95), and a copper disc painted black (ε 0.92) on
Kaowool, both at 50 kW/m². They test the shared convection coefficient and the
Kaowool's conductivity.
""",
    )
    rows = []
    for t in INERT_TESTS
        s = load_run("pyrolysis", "inert", "", t.test)
        s === nothing && continue
        a = read_matl_csv(inert_experiment_path(t.test))
        cols = collect(keys(s))[3:end]
        for (k, dpt) in enumerate(t.depths), tt in (t.t_end / 4, t.t_end / 2, t.t_end)
            i = findfirst(==(tt), s.time_s)
            j = findfirst(==(tt), a[:, 1])
            push!(
                rows,
                [
                    t.test,
                    @sprintf("%.2f", 1e3dpt),
                    @sprintf("%.0f", tt),
                    fmt(getfield(s, cols[k])[i], "%.1f"),
                    fmt(a[j, 1+k], "%.1f"),
                    fmt(getfield(s, cols[k])[i] - a[j, 1+k], "%+.1f"),
                ],
            )
        end
    end
    table(
        io,
        [
            "test",
            "depth (mm)",
            "time (s)",
            "Pyrolysis.jl (K)",
            "experiment (K)",
            "difference (K)",
        ],
        rows,
    )
    figure(io, inert_figure(), "Inert tests"; up = true)
    write_page(joinpath(DETAIL, APP_PAGE), io)

    # =========================================================================
    # Cost
    # =========================================================================
    timed = isfile(joinpath(TIMING, "pyrolysis_timing.csv"))
    cost_png, reached, pkeys, coarse, versus = "", [], [], [], String[]
    vs = Dict{String,NTuple{3,Int}}()     # code => (faster, shared cases, faster scaled)
    if timed
        crow = cost_rows(TIMING)
        io = IOBuffer()
        println(
            io,
            """
# Cost

The cost at matched accuracy of the $OVERVIEW. Each code is scored against its own finest
run of each case: Pyrolysis.jl against 80 cells/mm at reltol 1e-8, FDS against its 25 µm
run, Gpyro and ThermaKin2Ds against their half-cell, half-step runs. A run is within
tolerance when its peak MLR and t₅₀ are within 1 % and its surface and back temperature
histories within 2 K up to 90 % mass loss; the error is the largest of those four
ratios. Pyrolysis.jl runs a ladder of 10, 20 and 40 cells/mm at the default tolerances
and ten times them. Every time is the median of repeated solves on one pinned core of a
quiet machine (`scripts/timing/README.md`). FDS and ThermaKin2Ds step at a fixed step to
the case's end time; their times scaled to their own burnout are given beside them,
since Pyrolysis.jl and Gpyro stop at burnout.
""",
        )
        rows = []
        for key in sort(unique((r[2], r[3]) for r in crow)),
            code in ("Pyrolysis.jl", "fds", "gpyro", "thermakin")

            rs = filter(r -> (r[2], r[3]) == key && r[1] == code, crow)
            isempty(rs) && continue
            ok = filter(r -> r[6] <= 1, rs)
            best = isempty(ok) ? nothing : ok[argmin([r[5] for r in ok])]
            prac = filter(r -> r[4] in ("20/mm ×1", "base"), rs)
            ps = isempty(prac) ? NaN : span_scaled(code, key[1], key[2], "base", prac[1][5])
            push!(
                rows,
                [
                    key[1],
                    replace(key[2], "gasification_" => "", "apparatus_" => "apparatus "),
                    COST_NAME[code],
                    best === nothing ? "not reached" :
                    @sprintf("%.2f (%s)", best[5], best[4]),
                    isempty(prac) ? "–" : @sprintf("%.2f", prac[1][5]),
                    code in ("fds", "thermakin") && isfinite(ps) ? @sprintf("%.2f", ps) :
                    "–",
                ],
            )
        end
        table(
            io,
            [
                "set",
                "case",
                "code",
                "cheapest within tolerance (s)",
                "at 20 cells/mm or base (s)",
                "scaled to burnout (s)",
            ],
            rows,
        )
        ourrows(k) = filter(r -> (r[2], r[3]) == k && r[1] == "Pyrolysis.jl", crow)
        cheapest(k) = (
            ok = filter(r -> r[6] <= 1, ourrows(k));
            isempty(ok) ? NaN : minimum(r[5] for r in ok)
        )
        pkeys = sort(unique((r[2], r[3]) for r in crow if r[1] == "Pyrolysis.jl"))
        reached = filter(k -> isfinite(cheapest(k)), pkeys)
        coarse =
            filter(k -> any(r -> r[6] <= 1 && startswith(r[4], "10/"), ourrows(k)), reached)
        for code in ("fds", "gpyro", "thermakin")
            n, faster, scaled = 0, 0, 0
            for k in sort(unique((r[2], r[3]) for r in crow if r[1] == code))
                base =
                    filter(r -> (r[2], r[3]) == k && r[1] == code && r[4] == "base", crow)
                (isempty(base) || !isfinite(cheapest(k))) && continue
                n += 1
                faster += cheapest(k) < base[1][5]
                scaled += cheapest(k) < span_scaled(code, k[1], k[2], "base", base[1][5])
            end
            n == 0 && continue
            vs[code] = (faster, n, scaled)
            push!(
                versus,
                "$(COST_NAME[code])'s base run on $faster of $n" * (
                    code in ("fds", "thermakin") ?
                    " ($scaled of $n against its time scaled to burnout)" : ""
                ),
            )
        end
        println(
            io,
            """
Pyrolysis.jl reaches the tolerance on $(length(reached)) of its $(length(pkeys)) timed cases,
$(length(coarse)) of them at 10 cells/mm. On the cases where an external code is timed too,
its cheapest run within tolerance takes less time than $(join(versus, "; ")).
""",
        )
        cost_png = cost_figure(crow)
        figure(io, cost_png, "Solve time against error"; up = true)
        write_page(joinpath(DETAIL, COST_PAGE), io)
    end

    # =========================================================================
    # Overview
    # =========================================================================
    # Against the submissions: each headline set's deck as written, the
    # calibration code's own run of it, and the run as the group ran it.
    sub_rows = []
    as_written, own_runs = Dict{String,Float64}(), Dict{String,Float64}()
    for d in headline
        id = d["id"]
        p = medabs(d50(M, id, "pyrolysis/base", "submission"))
        as_written[id] = p
        o = own(id)
        o === nothing || (own_runs[id] = medabs(o))
        owncell =
            o !== nothing ? pct(medabs(o)) :
            "$(pct(tk[1])) / $(pct(tk[2])) (`MAX` / `RAND`)"
        asran =
            haskey(AS_RAN, id) ?
            pct(medabs(d50(M, id, "pyrolysis/$(AS_RAN[id])", "submission")); d = 2) : "–"
        push!(sub_rows, [id, d["calibration_code"], pct(p), owncell, asran])
    end
    nclose = count(id -> shown(as_written[id]) <= shown(own_runs[id]), keys(own_runs))
    nfar = count(id -> shown(as_written[id]) > shown(own_runs[id]), keys(own_runs))
    departed = [id for id in ("BUW-FZJ_B", "DBI_1") if haskey(as_written, id)]
    departed_to = [
        pct(medabs(d50(M, id, "pyrolysis/$(AS_RAN[id])", "submission")); d = 2) for
        id in departed
    ]
    # The two largest gaps are the scenario departures the bullet names.
    largest = sort(collect(keys(as_written)); by = id -> -as_written[id])[1:2]
    Set(largest) == Set(departed) ||
        error("the largest gaps to the submissions are $largest, not $departed")
    # Matched to FDS: every case of every in-depth set.
    ncases = sum(st.n for st in values(matched))
    worst_match = maximum(max(st.t50_max, st.peak_max) for st in values(matched))
    match_rows = [
        [
            d["id"],
            pct(base_fds[d["id"]].t50_med),
            pct(matched[d["id"]].t50_max; d = 2),
            pct(matched[d["id"]].peak_max; d = 2),
            "$(fmt(max(matched[d["id"]].Ttop, matched[d["id"]].Tback), "%.1f"))",
        ] for d in indepth if haskey(matched, d["id"]) && haskey(base_fds, d["id"])
    ]
    # TGA: the sets within 2 K, and the largest difference.
    tga_tol = 2.0
    tga_ok = count(v -> v[1] <= tga_tol, values(tga_dev))
    tga_worst = argmax(id -> tga_dev[id][1], collect(keys(tga_dev)))
    # Apparatus: per quantity, the medians over the sets the group also predicted.
    app_rows = [
        [
            "$(qname[c[1]]) kW/m²",
            "$(length(c[3])) of $(length(c[2]))",
            fmt(median([S[(id, c[1], "base")] for id in c[2]]), "%.3f"),
            fmt(median([S[(id, c[1], nothing)] for id in c[2]]), "%.3f"),
        ] for c in cmp
    ]
    # As good or better / sets with a group prediction, per kind of quantity.
    app_count(qs) = (
        sum(length(c[3]) for c in cmp if c[1] in qs),
        sum(length(c[2]) for c in cmp if c[1] in qs),
    )
    app_mass = app_count(("q50_Mass",))
    app_mlr = app_count(("q50_MLR",))
    app_T = app_count(("q50_Temp", "q25_Temp"))
    cost_bullet =
        timed && haskey(vs, "fds") ?
        "- **Cost.** Pyrolysis.jl reaches a 1 %, 2 K accuracy tolerance on $(length(reached)) of " *
        "$(length(pkeys)) timed\n  cases, and its cheapest such run takes less time than FDS's on " *
        "$(vs["fds"][1]) of $(vs["fds"][2]) shared cases." : ""

    io = IOBuffer()
    println(
        io,
        """
# MaCFP PMMA validation

The MaCFP condensed-phase working group's database holds property sets for one black
PMMA. Each set was calibrated by an institution in its own code (FDS, Gpyro or
ThermaKin2Ds) and comes with that institution's predictions: of the MaCFP-2 targets (TGA,
and gasification of 6 and 12 mm slabs at 10, 25 and 65 kW/m²) and, for most sets, of the
MaCFP-3 experiments in the NIST gasification apparatus. Pyrolysis.jl runs every set as
written, and FDS 6.11.1, Gpyro V2 and ThermaKin2Ds 7 run the same decks
($(link(DECKS_PAGE, "decks and scenario"))).

The submissions are calibrated models, not measurements: agreement with them shows how
faithfully a code reproduces another code's model. The apparatus experiments are the
measurements.

## Findings

- **Same deck, same answer.** On the decks that absorb at the surface, Pyrolysis.jl and
  FDS agree to $(pct(minimum(sfds)))–$(pct(maximum(sfds))) % in t₅₀. On the in-depth decks, with FDS's attenuation
  matched, all $ncases cases agree within $(pct(worst_match; d = 2)) % in t₅₀ and peak mass-loss rate.
- **The groups' predictions.** Pyrolysis.jl reproduces each submission about as well as
  the code the set was calibrated in: as close or closer on $nclose of $(length(own_runs)) sets. The two
  largest gaps are submissions that ran a different scenario; reproducing it brings
  $(join(["$id to $v %" for (id, v) in zip(departed, departed_to)], " and ")).
- **Where the codes differ**, the difference follows how each code deposits the
  radiation the solid absorbs.
- **TGA.** The kinetics reproduce the submitted curves within $(Int(tga_tol)) K on $tga_ok of
  $(length(tga_dev)) sets.
- **The apparatus experiments.** Scored against the measurements, Pyrolysis.jl is as
  good as or better than the groups' own predictions on mass-loss rate in $(app_mlr[1]) of $(app_mlr[2]) sets
  and on back temperature in $(app_T[1]) of $(app_T[2]), but on mass in $(app_mass[1]) of $(app_mass[2]). The shared
  apparatus deck limits the back-temperature scores.
$cost_bullet

## Against FDS on the same decks

On the decks that absorb at the surface, Pyrolysis.jl and FDS solve the same problem and
agree to $(pct(minimum(sfds)))–$(pct(maximum(sfds))) % in median t₅₀. On the in-depth decks both codes absorb
and emit in every cell, but FDS's two-flux sweep crosses a cell of width Δz with
transmission 1/(1 + 2κΔz) where Pyrolysis.jl's two-stream model transmits e^(−1.5κΔz).
Scaling κ so the two agree on FDS's 50 µm cells, and aligning the one item FDS cannot
express on UMD, UMET_GP and UMET_TK, leaves differences within $(pct(worst_match; d = 2)) % in every case.
In t₅₀ they are within FDS's own change between its 50 and 25 µm cells on
$(howmany(count(first, values(within)), length(within))) of the sets FDS also ran at 25 µm.
""",
    )
    table(
        io,
        [
            "set",
            "deck as written: median \\|Δt₅₀\\| (%)",
            "matched: max \\|Δt₅₀\\| (%)",
            "matched: max \\|Δpeak\\| (%)",
            "matched: worst RMS ΔT (K)",
        ],
        match_rows,
    )
    println(
        io,
        """
Per-code tables: $(link(GAS_PAGE, "code to code", "Code to code")) and
$(link(GAS_PAGE, "matched to FDS", "Matched to FDS")).

## Against the groups' predictions

The code a set was calibrated in, running the published deck, is the yardstick: it
reproduces the submission only as far as the deck and the scenario are what the group
ran. Median |Δt₅₀| over the six cases (%):
""",
    )
    table(
        io,
        [
            "set",
            "calibrated in",
            "Pyrolysis.jl",
            "calibration code",
            "Pyrolysis.jl, as the group ran",
        ],
        sub_rows,
    )
    println(
        io,
        """
UMD's ThermaKin2Ds run is read in `MAX` mode at 10 kW/m² and `RAND` above, which
reproduces its submission; UMET_TK's mode is not known. The last column reproduces what
three groups ran instead of the specified scenario: convection at h = 5.33 W/(m²·K)
(BUW-FZJ_B), c and k tables read as °C (DBI_1), and a 303.15 K start (Aalto_I).

BUW-FZJ_A is left out. Its low-emissivity residue collects as a film at the surface, and
how much of that film the surface cell holds, and so how much of the heater's flux the
sample absorbs, depends on the cell size in every code. Its results measure the
discretization as much as the property set ($(link(GAS_PAGE, "BUW-FZJ_A", "BUW-FZJ_A"))).
""",
    )
    figure(io, parity, "Parity of t50 and peak MLR against the submissions")
    println(
        io,
        """
Every code against every submission: $(link(GAS_PAGE, "gasification", "Against the submissions"));
the mass, mass-loss rate and temperature histories of each set:
$(link(HIST_PAGE, "gasification histories")).

## Where the codes differ

The differences to Gpyro and ThermaKin2Ds follow how each code deposits the radiation it
absorbs. Variants of each in-depth deck that move one choice at a time show it:

- Gpyro deposits the net surface flux along a Beer–Lambert march and emits nothing from
  the depth; the shallowest-depositing variant is closest to it on $(length(gp_close)) of $(length(gp)) sets.
- ThermaKin2Ds in `MAX` mode deposits the absorbed radiation in the surface element;
  surface absorption is the closest variant on $(length(tk_surf)) of $(length(tkr)) sets.
- A diffuse heater is closer than a collimated one to the submission on $(length(heater_sub)) of $(length(calib_sub))
  sets calibrated in Gpyro or ThermaKin2Ds, and to FDS on $(length(heater_fds)) of $(length(ids_with("fds/base"))).
- Absorbing at the surface instead of in depth moves t₅₀ by
  $(pct(minimum(surface_move)))–$(pct(maximum(surface_move))) % per set.

Every variant against every code: $(link(GAS_PAGE, "radiation in depth", "Radiation in depth")).

## TGA

The TGA targets test the kinetics alone. Pyrolysis.jl's temperatures at 10, 50 and
90 % mass loss are within $(Int(tga_tol)) K of the submission on $tga_ok of the $(length(tga_dev)) sets at both
heating rates; the largest difference is $(fmt(tga_dev[tga_worst][1], "%.1f")) K ($tga_worst, $(tga_dev[tga_worst][2]) K/min).
""",
    )
    figure(io, tga_png, "TGA, every set")
    println(
        io,
        """
Per set and rate: $(link(TGA_PAGE, "TGA")).

## Against the apparatus experiments

Every deck with complete properties runs the eleven PMMA tests of the NIST gasification
apparatus on one shared apparatus deck, scored against the measured repeats with MaCFP's
normalized RMS (lower is better). Against the groups' own predictions of the same tests,
over the sets with a group prediction:
""",
    )
    table(
        io,
        [
            "quantity",
            "Pyrolysis.jl as good or better",
            "median NRMS, Pyrolysis.jl",
            "median NRMS, groups",
        ],
        app_rows,
    )
    println(
        io,
        """
The back temperatures measure the shared deck as much as the property sets: on the
apparatus's inert Kaowool stack, with no sample, it runs $(fmt(deep, "%.0f")) K above the measurement at
17.16 mm after 1200 s.

Scores, histories and the inert tests: $(link(APP_PAGE, "gasification apparatus")).
""",
    )
    if timed
        println(
            io,
            """
## Cost

Each code is scored against its own finest run of each case, and a run is within
tolerance when its peak mass-loss rate and t₅₀ are within 1 % and its temperature
histories within 2 K. Pyrolysis.jl reaches the tolerance on $(length(reached)) of its $(length(pkeys)) timed
cases, $(length(coarse)) of them at 10 cells/mm. Its cheapest run within tolerance takes less time than
$(join(versus, "; ")).
""",
        )
        figure(io, cost_png, "Solve time against error")
        println(io, "Every case and code: $(link(COST_PAGE, "cost")).\n")
    end
    write_page(PAGE, io)

    total = sum(filesize(joinpath(FIG, f)) for f in readdir(FIG))
    @printf(
        "wrote %s, %d detail pages and %d figures (%.2f MB)\n",
        PAGE,
        length(readdir(DETAIL)),
        length(readdir(FIG)),
        total / 2^20
    )
end

main()
