# Figures and metrics for the MaCFP-2 PMMA targets: every code's run of each
# deck against the institution's own submission.
#
#   julia --project=scripts/macfp_pmma scripts/macfp_pmma/compare.jl
#
# Reads results/{pyrolysis,fds,gpyro,thermakin}/<set>/<variant>/*.csv and the
# matl-db submissions; writes results/figures/*.png and results/metrics.{csv,md}.
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, CairoMakie, Printf, Statistics

const RES = joinpath(@__DIR__, "results")
const FIG = joinpath(RES, "figures")

# Series identity: the code (fixed colour and line style), the submission in ink.
const INK = colorant"#0b0b0b"
const GRID = colorant"#e4e3df"
const SERIES = [
    # (code, variant, label, colour, linestyle)
    ("pyrolysis", "base", "Pyrolysis.jl", colorant"#2a78d6", :solid),
    (
        "pyrolysis",
        "collimated",
        "Pyrolysis.jl, collimated heater",
        colorant"#2a78d6",
        (:dash, :dense),
    ),
    (
        "pyrolysis",
        "fds_closure",
        "Pyrolysis.jl, FDS's cell attenuation",
        colorant"#8a4fd0",
        :dash,
    ),
    ("pyrolysis", "fds_matched", "Pyrolysis.jl, matched to FDS", colorant"#8a4fd0", :dot),
    (
        "pyrolysis",
        "surface",
        "Pyrolysis.jl, surface absorption",
        colorant"#2a78d6",
        (:dot, :loose),
    ),
    (
        "pyrolysis",
        "h5.33",
        "Pyrolysis.jl, h = 5.33 W/m²K (as their runs)",
        colorant"#2a78d6",
        :dashdot,
    ),
    (
        "pyrolysis",
        "tablesC",
        "Pyrolysis.jl, tables read in °C (as their FDS)",
        colorant"#2a78d6",
        :dashdot,
    ),
    (
        "pyrolysis",
        "T303",
        "Pyrolysis.jl, 303 K start (as their runs)",
        colorant"#2a78d6",
        :dashdot,
    ),
    ("fds", "base", "FDS 6.11.1", colorant"#eb6834", (:dot, :dense)),
    ("gpyro", "base", "Gpyro V2", colorant"#1baf7a", :dashdot),
    ("gpyro", "rhoT", "Gpyro V2, ρ(T) (as UMET ran)", colorant"#1baf7a", (:dot, :dense)),
    ("thermakin", "base", "ThermaKin2Ds (MAX)", colorant"#e87ba4", (:dash, :loose)),
    ("thermakin", "rand", "ThermaKin2Ds (RAND)", colorant"#6250d6", (:dot, :loose)),
    (
        "thermakin",
        "rhoT",
        "ThermaKin2Ds, ρ(T) (as UMET ran)",
        colorant"#e87ba4",
        (:dot, :dense),
    ),
]

const TOP_LABEL = Dict(
    :T_top_K => "T surface (K)",
    :T_avg1mm_K => "T mean of top 1 mm (K)",
    :T_3pt_K => "T mean of 0, 0.5, 1 mm (K)",
)

function gasification_figure(deck)
    id = deck["id"]
    cases = MacfpPMMA.GASIFICATION_CASES
    fig = Figure(size = (2600, 1500), fontsize = 15, backgroundcolor = :white)
    toplabel = TOP_LABEL[top_column(deck)]
    rows = [
        (:mass_g, "mass (g, 10×10 cm)"),
        (:mlr_g_m2s, "MLR (g/m²·s)"),
        (:top, toplabel),
        (:T_back_K, "T back (K)"),
    ]
    used = Dict{String,Any}()
    for (j, (q, L, t_end)) in enumerate(cases)
        case = "gasification_$(q)kW_$(L)mm"
        sub = load_submission(deck, q, L)
        for (i, (col, ylabel)) in enumerate(rows)
            ax = Axis(
                fig[i+1, j];
                title = i == 1 ? "$q kW/m², $L mm" : "",
                xlabel = i == length(rows) ? "time (s)" : "",
                ylabel = j == 1 ? ylabel : "",
                xgridcolor = GRID,
                ygridcolor = GRID,
                topspinevisible = false,
                rightspinevisible = false,
            )
            tmax = 0.0
            for (code, variant, label, colour, ls) in SERIES
                s = load_run(code, id, variant, case)
                s === nothing && continue
                c = col == :top ? top_column(deck) : col
                hasproperty(s, c) || continue      # external codes write the surface and point values only
                y = getfield(s, c)
                lines!(
                    ax,
                    s.time_s,
                    y;
                    color = colour,
                    linestyle = ls,
                    linewidth = 2,
                    label = label,
                )
                used[label] = (colour, ls)
                tmax = max(tmax, burnout_time(s))
            end
            if sub !== nothing
                c = col == :top ? :T_top_K : col
                lines!(
                    ax,
                    sub.time_s,
                    getfield(sub, c);
                    color = INK,
                    linewidth = 1.5,
                    label = "$(deck["predictions"]) submission ($(deck["calibration_code"]))",
                )
                used["$(deck["predictions"]) submission ($(deck["calibration_code"]))"] =
                    (INK, :solid)
                tmax = max(tmax, burnout_time(sub))
            end
            tmax > 0 && xlims!(ax, 0, 1.15 * tmax)
        end
    end
    labels = collect(keys(used))
    order = sortperm(
        labels;
        by = l -> occursin("submission", l) ? 0 : findfirst(s -> s[3] == l, SERIES),
    )
    labels = labels[order]
    elems = [
        LineElement(color = used[l][1], linestyle = used[l][2], linewidth = 3) for
        l in labels
    ]
    Legend(
        fig[1, 1:length(cases)],
        elems,
        labels;
        orientation = :horizontal,
        framevisible = false,
        tellheight = true,
    )
    Label(
        fig[0, 1:length(cases)],
        "MaCFP-2 PMMA gasification targets, property set $id (calibrated in $(deck["calibration_code"]))";
        fontsize = 22,
        font = :bold,
        halign = :left,
    )
    mkpath(FIG)
    save(joinpath(FIG, "gasification_$(id).png"), fig)
end

function tga_figure(decks)
    n = length(decks)
    ncol = 4
    fig = Figure(
        size = (1900, 420 * cld(n, ncol) + 80),
        fontsize = 15,
        backgroundcolor = :white,
    )
    for (k, deck) in enumerate(decks)
        id = deck["id"]
        ax = Axis(
            fig[cld(k, ncol)+1, mod1(k, ncol)];
            title = id,
            xlabel = "T (K)",
            ylabel = "m/m₀",
            xgridcolor = GRID,
            ygridcolor = GRID,
            topspinevisible = false,
            rightspinevisible = false,
        )
        for (β, ls) in zip(TGA_RATES, (:solid, :dash))
            s = load_run("pyrolysis", id, "base", "tga_$(β)K")
            s === nothing || lines!(
                ax,
                s.T_K,
                s.mass_rel;
                color = colorant"#2a78d6",
                linestyle = ls,
                linewidth = 2,
            )
            p = tga_prediction_path(deck, β)
            if p !== nothing
                a = read_matl_csv(p)
                # Sandia files: time, T [°C], T [K], mass-loss fraction, mass [mg]; the rest: time, T [K], mass.
                iT, im = size(a, 2) >= 5 ? (3, 5) : (2, 3)
                size(a, 1) > 1 && lines!(
                    ax,
                    a[:, iT],
                    a[:, im] ./ a[1, im];
                    color = INK,
                    linestyle = ls,
                    linewidth = 1.2,
                )
            end
        end
        xlims!(ax, 450, 800)
    end
    elems = [
        LineElement(color = colorant"#2a78d6", linewidth = 3),
        LineElement(color = INK, linewidth = 3),
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
    save(joinpath(FIG, "tga_all.png"), fig)
end

function metrics_table(decks)
    rows = String[]
    push!(rows, "set,case,source,peak_mlr,t_peak,t50,t90,m_end")
    for deck in decks, (q, L, _) in MacfpPMMA.GASIFICATION_CASES
        id = deck["id"]
        case = "gasification_$(q)kW_$(L)mm"
        sub = load_submission(deck, q, L)
        if sub !== nothing
            m = metrics(sub)
            push!(
                rows,
                @sprintf("%s,%s,submission,%.4g,%.4g,%.4g,%.4g,%.4g", id, case, m...)
            )
        end
        for (code, variant, label, _, _) in SERIES
            s = load_run(code, id, variant, case)
            s === nothing && continue
            m = metrics(s)
            push!(
                rows,
                @sprintf(
                    "%s,%s,%s/%s,%.4g,%.4g,%.4g,%.4g,%.4g",
                    id,
                    case,
                    code,
                    variant,
                    m...
                )
            )
        end
    end
    write(joinpath(RES, "metrics.csv"), join(rows, "\n") * "\n")
    return rows
end

# Parity of each code's metric against the submission's, all sets and cases.
function parity_figure(rows)
    data = [split(r, ",") for r in rows[2:end]]
    sub = Dict((d[1], d[2]) => d for d in data if d[3] == "submission")
    fig = Figure(size = (1800, 560), fontsize = 15, backgroundcolor = :white)
    quantities = [
        (4, "peak MLR (g/m²·s)"),
        (6, "time to 50 % mass loss (s)"),
        (7, "time to 90 % mass loss (s)"),
    ]
    used = Dict{String,Any}()
    for (j, (col, name)) in enumerate(quantities)
        ax = Axis(
            fig[2, j];
            title = name,
            xlabel = "submission",
            ylabel = j == 1 ? "code run of the same deck" : "",
            xgridcolor = GRID,
            ygridcolor = GRID,
            aspect = 1,
            topspinevisible = false,
            rightspinevisible = false,
        )
        lo, hi = Inf, -Inf
        groups = [
            (label, colour, [(code, variant)], :circle) for
            (code, variant, label, colour, _) in SERIES if variant == "base"
        ]
        push!(
            groups,
            (
                "Pyrolysis.jl, as the group ran (h5.33 / tablesC / T303)",
                colorant"#2a78d6",
                [("pyrolysis", "h5.33"), ("pyrolysis", "tablesC"), ("pyrolysis", "T303")],
                :utriangle,
            ),
        )
        for (label, colour, members, marker) in groups
            xs, ys = Float64[], Float64[]
            for d in data
                any(m -> d[3] == "$(m[1])/$(m[2])", members) || continue
                s = get(sub, (d[1], d[2]), nothing)
                s === nothing && continue
                x, y = parse(Float64, s[col]), parse(Float64, d[col])
                (isfinite(x) && isfinite(y)) || continue
                push!(xs, x)
                push!(ys, y)
            end
            isempty(xs) && continue
            scatter!(
                ax,
                xs,
                ys;
                color = marker == :circle ? colour : (:white, 0.0),
                marker,
                markersize = 11,
                strokecolor = marker == :circle ? :white : colour,
                strokewidth = marker == :circle ? 1 : 1.5,
            )
            used[label] = (colour, marker)
            lo, hi = min(lo, minimum(xs), minimum(ys)), max(hi, maximum(xs), maximum(ys))
        end
        isfinite(lo) && lines!(ax, [lo, hi], [lo, hi]; color = :gray50, linewidth = 1)
    end
    labels = collect(keys(used))
    Legend(
        fig[1, 1:3],
        [
            used[l][2] == :circle ?
            MarkerElement(color = used[l][1], marker = :circle, markersize = 12) :
            MarkerElement(
                color = (:white, 0.0),
                marker = used[l][2],
                markersize = 12,
                strokecolor = used[l][1],
                strokewidth = 1.5,
            ) for l in labels
        ],
        labels;
        orientation = :horizontal,
        framevisible = false,
    )
    save(joinpath(FIG, "parity.png"), fig)
end

# ---------------------------------------------------------------------------
# MaCFP-3 apparatus
# ---------------------------------------------------------------------------

const EXP = colorant"#9a9890"     # experiments: light ink, drawn first

function apparatus_figure(deck; variant = "base")
    id = deck["id"]
    fig = Figure(size = (2200, 560), fontsize = 15, backgroundcolor = :white)
    for (j, (key, reps, col, ylabel, title)) in enumerate(APPARATUS_PANELS)
        ax = Axis(
            fig[2, j];
            title,
            xlabel = "time (s)",
            ylabel,
            xgridcolor = GRID,
            ygridcolor = GRID,
            topspinevisible = false,
            rightspinevisible = false,
        )
        tmax = 0.0
        for r in reps
            te, ye = apparatus_experiment(key, r)
            lines!(ax, te, ye; color = EXP, linewidth = 3)
            tmax = max(tmax, maximum(te))
        end
        for f in group_predictions(deck, key)
            a = read_matl_csv(f)
            size(a, 1) > 1 && lines!(ax, a[:, 1], a[:, 2]; color = INK, linewidth = 1.2)
        end
        for r in reps
            test = replace(key, "MLR" => "Mass") * "_$r"
            s = load_run("pyrolysis", id, variant, "apparatus_$test")
            s === nothing && continue
            lines!(
                ax,
                s.time_s,
                getfield(s, col);
                color = colorant"#2a78d6",
                linewidth = 1.6,
            )
        end
        xlims!(ax, 0, 1.1 * tmax)
        col == :mass_g && ylims!(ax, -2, nothing)
    end
    labels = [
        "experiment (each repeat)",
        "$(get(deck, "validation", "—") |> basename) prediction",
        "Pyrolysis.jl (shared apparatus deck, cone heater)",
    ]
    elems = Any[
        LineElement(color = EXP, linewidth = 4),
        LineElement(color = INK, linewidth = 2),
        LineElement(color = colorant"#2a78d6", linewidth = 3),
    ]
    Legend(fig[1, 1:4], elems, labels; orientation = :horizontal, framevisible = false)
    Label(
        fig[0, 1:4],
        "MaCFP-3 NIST gasification apparatus, property set $id";
        fontsize = 20,
        font = :bold,
        halign = :left,
    )
    save(joinpath(FIG, "apparatus_$(id).png"), fig)
end

# The apparatus heaters Pyrolysis.jl runs: the cone, and the diffuse and
# collimated heaters that bound it (in-depth decks).
const APPARATUS_HEATERS = [
    ("base", "Pyrolysis.jl"),
    ("diffuse", "Pyrolysis.jl (diffuse heater)"),
    ("collimated", "Pyrolysis.jl (collimated heater)"),
]

function apparatus_scores(decks)
    rows = ["set,quantity,source,nrms"]
    for deck in decks, panel in APPARATUS_PANELS
        for (variant, source) in APPARATUS_HEATERS
            v = apparatus_nrms(deck, panel; variant)
            isnan(v) ||
                push!(rows, @sprintf("%s,%s,%s,%.4f", deck["id"], panel[1], source, v))
        end
        v = apparatus_nrms(deck, panel; variant = nothing)
        isnan(v) || push!(rows, @sprintf("%s,%s,group,%.4f", deck["id"], panel[1], v))
    end
    write(joinpath(RES, "apparatus_scores.csv"), join(rows, "\n") * "\n")
    return rows
end

function scores_figure(rows)
    data = [split(r, ",") for r in rows[2:end]]
    sets = unique(d[1] for d in data)
    fig = Figure(size = (2000, 520), fontsize = 15, backgroundcolor = :white)
    for (j, (key, _, _, _, title)) in enumerate(APPARATUS_PANELS)
        ax = Axis(
            fig[2, j];
            title,
            ylabel = j == 1 ? "NRMS vs experimental mean" : "",
            xticks = (1:length(sets), sets),
            xticklabelrotation = π / 3,
            xgridcolor = GRID,
            ygridcolor = GRID,
            topspinevisible = false,
            rightspinevisible = false,
        )
        for (src, colour, mk, dx) in (
            ("group", INK, :rect, 0.12),
            ("Pyrolysis.jl", colorant"#2a78d6", :circle, -0.12),
        )
            xs, ys = Float64[], Float64[]
            for (i, s) in enumerate(sets), d in data
                (d[1] == s && d[2] == key && d[3] == src) || continue
                v = parse(Float64, d[4])
                isfinite(v) && (push!(xs, i + dx); push!(ys, v))
            end
            scatter!(
                ax,
                xs,
                ys;
                color = colour,
                marker = mk,
                markersize = 12,
                strokecolor = :white,
                strokewidth = 1,
            )
        end
        ylims!(ax, 0, nothing)
    end
    Legend(
        fig[1, 1:4],
        [
            MarkerElement(color = colorant"#2a78d6", marker = :circle, markersize = 12),
            MarkerElement(color = INK, marker = :rect, markersize = 12),
        ],
        ["Pyrolysis.jl, shared apparatus deck", "the group's own MaCFP-3 prediction"];
        orientation = :horizontal,
        framevisible = false,
    )
    save(joinpath(FIG, "apparatus_scores.png"), fig)
end

function main(args = ARGS)
    check_matl_db()
    mkpath(FIG)
    all_decks = load_decks()
    gas_decks = filter(d -> !get(d, "kinetics_only", false), all_decks)
    if !("--apparatus-only" in args)
        foreach(gasification_figure, filter(d -> haskey(d, "predictions"), gas_decks))
        tga_figure(filter(d -> haskey(d, "tga_predictions"), all_decks))
        rows = metrics_table(filter(d -> haskey(d, "predictions"), gas_decks))
        parity_figure(rows)
    end
    have_apparatus = filter(
        d -> isfile(
            joinpath(RES, "pyrolysis", d["id"], "base", "apparatus_q50_Mass_R3.csv"),
        ),
        gas_decks,
    )
    if !isempty(have_apparatus)
        foreach(apparatus_figure, have_apparatus)
        scores_figure(apparatus_scores(have_apparatus))
    end
    println("figures in $FIG")
end

main()
