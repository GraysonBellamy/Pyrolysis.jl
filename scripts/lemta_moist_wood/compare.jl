# Score every variant run against the data: per-run metrics
# (results/metrics_<group>.csv), summary tables in Markdown
# (results/tables/*.md) and figures (results/figures/). The tables of the
# phase change (its summary on the long and on all runs, its runs and their
# cooling) are also printed.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/compare.jl
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using CairoMakie
using Printf
using Statistics

runs = Dict(r.name => r for r in load_runs())
tables = joinpath(RESULTS, "tables")
figdir = joinpath(RESULTS, "figures")
mkpath(tables)
mkpath(figdir)

variants_present() = sort([v for v in keys(VARIANTS) if isdir(joinpath(RESULTS, v))])
series_of(v) = Dict(
    splitext(f)[1] => read_series(joinpath(RESULTS, v, f)) for
    f in readdir(joinpath(RESULTS, v)) if
    endswith(f, ".csv") && !endswith(f, "_profiles.csv")
)
fmt(x; d = 3) = isnan(x) ? "–" : @sprintf("%.*f", d, x)
mm(d) = round(Int, 1000d)

# Thermocouple accuracy (IEC 60584 class quoted by Diarra et al.): ±2 K up to
# 100 °C, ±1.5 % of the reading in °C above.
tc_accuracy(T) = max(2.0, 0.015 * (T - 273.15))

"""Temperature errors of a run against a series over `window`: per depth (rms K, bias K, NRMS), all depths pooled."""
function temperature_errors(r, s; window = (0.0, Inf))
    t_hi = min(window[2], s.t[end], r.t[end])
    keep = findall(t -> window[1] <= t <= t_hi, r.t)
    per = Dict{Int,NamedTuple}()
    pooled = Float64[]
    within = 0
    for (j, d) in enumerate(r.depths)
        js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
        y = r.T[keep, j]
        m = [interp(s.t, view(s.T, :, js), t) for t in r.t[keep]]
        e = m .- y
        append!(pooled, e)
        within += count(abs.(e) .<= tc_accuracy.(y))
        per[mm(d)] = (
            rms = sqrt(mean(e .^ 2)),
            bias = mean(e),
            nrms = sqrt(mean(e .^ 2)) / (maximum(y) - minimum(y)),
        )
    end
    return (; per, rms = sqrt(mean(pooled .^ 2)), within = within / length(pooled))
end

# ---------------------------------------------------------------------------
# Dry runs
# ---------------------------------------------------------------------------

dry_variants = filter(v -> VARIANTS[v].runs == :dry, variants_present())
if !isempty(dry_variants)
    open(joinpath(RESULTS, "metrics_dry.csv"), "w") do io
        println(
            io,
            "variant,run,rms_K,within_accuracy,mean_nrms,bias_2mm,bias_10mm,bias_20mm,bias_30mm,bias_40mm",
        )
        for v in dry_variants, (name, s) in sort(collect(series_of(v)); by = first)
            r = runs[name]
            e = temperature_errors(r, s)
            b(d) = haskey(e.per, d) ? e.per[d].bias : NaN
            @printf(
                io,
                "%s,%s,%.3f,%.3f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                v,
                name,
                e.rms,
                e.within,
                mean(x.nrms for x in values(e.per)),
                b(2),
                b(10),
                b(20),
                b(30),
                b(40)
            )
        end
    end
    open(joinpath(tables, "dry.md"), "w") do io
        println(
            io,
            "| variant | run | rms (K) | within TC accuracy | mean NRMS | bias 2 / 10 / 20 / 30 / 40 mm (K) |",
        )
        println(io, "|---|---|---|---|---|---|")
        for v in dry_variants, (name, s) in sort(collect(series_of(v)); by = first)
            e = temperature_errors(runs[name], s)
            b(d) = haskey(e.per, d) ? @sprintf("%+.1f", e.per[d].bias) : "–"
            @printf(
                io,
                "| %s | %s | %.2f | %.0f %% | %.3f | %s / %s / %s / %s / %s |\n",
                v,
                name,
                e.rms,
                100e.within,
                mean(x.nrms for x in values(e.per)),
                b(2),
                b(10),
                b(20),
                b(30),
                b(40)
            )
        end
    end
    # per-depth bias profile at the end of heating, base deck
    if "dry" in dry_variants
        S = series_of("dry")
        local fig = Figure(size = (1200, 300 * cld(length(S), 2)))
        for (i, (name, s)) in enumerate(sort(collect(S); by = first))
            r = runs[name]
            ax = Axis(
                fig[cld(i, 2), mod1(i, 2)];
                title = name,
                xlabel = "time (h)",
                ylabel = "T (°C)",
            )
            for (j, d) in enumerate(r.depths)
                mm(d) in (2, 4, 6, 10, 14, 20, 30, 40) || continue
                js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
                c = Makie.wong_colors()[mod1(j, 7)]
                lines!(ax, r.t ./ 3600, r.T[:, j] .- 273.15; color = c)
                lines!(ax, s.t ./ 3600, s.T[:, js] .- 273.15; color = c, linestyle = :dash)
            end
        end
        save(joinpath(figdir, "dry_runs.png"), fig)
    end
end

# ---------------------------------------------------------------------------
# Moist runs
# ---------------------------------------------------------------------------

"""The moist-run metrics of one series against its run."""
function moist_metrics(r, s)
    t25 = min(9000.0, r.t[end], s.t[end])
    tend = min(r.t[end], s.t[end])
    e = temperature_errors(r, s)
    eh = temperature_errors(r, s; window = (0.0, t25))
    out = Dict{String,Float64}(
        "rms_K" => e.rms,
        "rms_heating_K" => eh.rms,
        "mean_nrms" => mean(x.nrms for x in values(e.per)),
        "loss_2p5h_exp" => mass_loss(r, t25),
        "loss_2p5h_mod" => mass_loss(s, t25),
        "loss_end_exp" => mass_loss(r, tend),
        "loss_end_mod" => mass_loss(s, tend),
        "water_g" => initial_water_g(r),
        "t25" => t25,
        "p_max_bar" => maximum(s.p_max) / 1e5,
        "pv_psat_max" => maximum(s.supersaturation),
    )
    # share of the water evaporated by 2.5 h that went where p_v > p_sat(T)
    k = findall(t -> t <= t25, s.t)
    tot = sum(s.evaporation[k])
    sup = sum(s.evaporation_supersaturated[k])
    out["evap_supersaturated_share"] = tot > 0 ? sup / tot : NaN
    for d in (6, 10, 14)
        j = findfirst(x -> mm(x) == d, r.depths)
        js = findfirst(x -> mm(x) == d, s.depths)
        if j === nothing
            out["cross_$(d)_exp"] =
                out["cross_$(d)_mod"] = out["ratio_$(d)_exp"] = out["ratio_$(d)_mod"] = NaN
            continue
        end
        out["cross_$(d)_exp"] = crossing_minutes(r.t, view(r.T, :, j))
        out["cross_$(d)_mod"] = crossing_minutes(s.t, view(s.T, :, js))
        out["ratio_$(d)_exp"] = rate_ratio(r.t, view(r.T, :, j))
        out["ratio_$(d)_mod"] = rate_ratio(s.t, view(s.T, :, js))
    end
    for (d, x) in e.per
        out["nrms_$(d)"] = x.nrms
        out["bias_$(d)"] = x.bias
    end
    # the phase change: the water evaporated and condensed over the record,
    # the vapour over its equilibrium pressure, the moisture content, the
    # steps and where it evaporates
    pc = s.phase_change
    if pc !== nothing
        out["evaporated_g"] = trapz(s.t, s.evaporation)
        out["condensed_g"] = trapz(s.t, pc.condensation)
        out["pv_peq_max"] = maximum(pc.pv_over_peq)
        out["X_max"] = maximum(pc.X_max)
        out["accepted"] = pc.accepted
        out["rejected"] = pc.rejected
        for (q, T) in zip((10, 50, 90), pc.T_evap)
            out["Tevap_p$(q)_C"] = T - 273.15
        end
    end
    return out
end

trapz(t, y) = sum(0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))

moist_variants = filter(v -> VARIANTS[v].runs == :moist, variants_present())
M = Dict{Tuple{String,String},Dict{String,Float64}}()
for v in moist_variants, (name, s) in series_of(v)
    M[(v, name)] = moist_metrics(runs[name], s)
end

if !isempty(moist_variants)
    keys_out = [
        "rms_K",
        "rms_heating_K",
        "mean_nrms",
        "water_g",
        "loss_2p5h_exp",
        "loss_2p5h_mod",
        "loss_end_exp",
        "loss_end_mod",
        "cross_10_exp",
        "cross_10_mod",
        "cross_14_exp",
        "cross_14_mod",
        "ratio_10_exp",
        "ratio_10_mod",
        "ratio_14_exp",
        "ratio_14_mod",
        "p_max_bar",
        "pv_psat_max",
        "evap_supersaturated_share",
        "pv_peq_max",
        "X_max",
        "evaporated_g",
        "condensed_g",
        "Tevap_p10_C",
        "Tevap_p50_C",
        "Tevap_p90_C",
        "accepted",
        "rejected",
    ]
    open(joinpath(RESULTS, "metrics_moist.csv"), "w") do io
        println(io, join(["variant", "run", keys_out...], ","))
        for ((v, name), m) in sort(collect(M); by = first)
            println(
                io,
                join(
                    [
                        v,
                        "\"$name\"",
                        (@sprintf("%.7g", get(m, k, NaN)) for k in keys_out)...,
                    ],
                    ",",
                ),
            )
        end
    end

    moist_runs = sort(unique(n for (_, n) in keys(M)))
    function law_tables(laws, tag)
        laws = filter(v -> v in moist_variants, laws)
        isempty(laws) && return
        open(joinpath(tables, "moist_runs$tag.md"), "w") do io
            println(
                io,
                "| run | water (g) | lost by 2.5 h (g): data / " *
                join(laws, " / ") *
                " | rms (K): " *
                join(laws, " / ") *
                " | 85→115 °C at 10 mm (min): data / " *
                join(laws, " / ") *
                " | rate ratio at 10 mm: data / " *
                join(laws, " / ") *
                " |",
            )
            println(io, "|---|---|---|---|---|---|")
            for name in moist_runs
                haskey(M, (laws[1], name)) || continue
                m1 = M[(laws[1], name)]
                g(l, k, d) = haskey(M, (l, name)) ? fmt(M[(l, name)][k]; d) : "–"
                println(
                    io,
                    "| $name | $(fmt(m1["water_g"]; d = 1)) | $(fmt(m1["loss_2p5h_exp"]; d = 1)) / " *
                    join([g(l, "loss_2p5h_mod", 1) for l in laws], " / ") *
                    " | " *
                    join([g(l, "rms_K", 1) for l in laws], " / ") *
                    " | $(fmt(m1["cross_10_exp"]; d = 0)) / " *
                    join([g(l, "cross_10_mod", 0) for l in laws], " / ") *
                    " | $(fmt(m1["ratio_10_exp"]; d = 2)) / " *
                    join([g(l, "ratio_10_mod", 2) for l in laws], " / ") *
                    " |",
                )
            end
        end
        # NRMS by depth, mean over runs
        depths = (2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40)
        open(joinpath(tables, "moist_nrms_by_depth$tag.md"), "w") do io
            println(io, "| law | " * join(["$d mm" for d in depths], " | ") * " |")
            println(io, "|---|" * repeat("---|", length(depths)))
            for l in laws
                cells = String[]
                for d in depths
                    x = [
                        M[(l, n)]["nrms_$d"] for n in moist_runs if
                        haskey(M, (l, n)) && haskey(M[(l, n)], "nrms_$d")
                    ]
                    push!(cells, isempty(x) ? "–" : @sprintf("%.3f", mean(x)))
                end
                println(io, "| $l | " * join(cells, " | ") * " |")
            end
        end
        # Mass change during cooling, runs with a cooling record
        open(joinpath(tables, "moist_cooling$tag.md"), "w") do io
            println(
                io,
                "| run | lost after the cone is withdrawn (g): data / " *
                join(laws, " / ") *
                " |",
            )
            println(io, "|---|---|")
            for name in moist_runs
                isfinite(runs[name].t_off) && haskey(M, (laws[1], name)) || continue
                c(m, w) = m["loss_end_$w"] - m["loss_2p5h_$w"]
                println(
                    io,
                    "| $name | $(fmt(c(M[(laws[1], name)], "exp"); d = 2)) / " *
                    join(
                        [
                            haskey(M, (l, name)) ? fmt(c(M[(l, name)], "mod"); d = 2) : "–"
                            for l in laws
                        ],
                        " / ",
                    ) *
                    " |",
                )
            end
        end
    end
    law_tables(["arrhenius_cal", "gated_cal", "chan_cal"], "")
    law_tables(["arrhenius", "gated", "chan"], "_published")
    law_tables(
        [
            "onset343_cal",
            "onset353_cal",
            "onset363_cal",
            "gated_cal",
            "onset373w10_cal",
            "onset383_cal",
        ],
        "_onset",
    )
    base = filter(
        v -> v in moist_variants,
        ["sorption_cal", "arrhenius_cal", "gated_cal", "chan_cal"],
    )

    # Summary per variant
    # The runs with a cooling record are the reused-sample tests of Diarra et al.
    # Table 1 (and two more at 7.3 %); the others are 2.3–2.6 h repeats.
    long_runs = [n for n in moist_runs if isfinite(runs[n].t_off)]
    subset = moist_runs
    summary(v, f) = [f(M[(v, n)]) for n in subset if haskey(M, (v, n))]
    function finite_mean(x)
        y = filter(isfinite, x)
        isempty(y) ? NaN : mean(y)
    end
    function finite_median(x)
        y = filter(isfinite, x)
        isempty(y) ? NaN : median(y)
    end
    for (sub, file) in
        ((moist_runs, "moist_summary.md"), (long_runs, "moist_summary_long.md"))
        global subset = sub
        open(joinpath(tables, file), "w") do io
            println(
                io,
                "| variant | runs | rms T (K), median | mean NRMS, median | mass lost at 2.5 h, model/data, mean (range) | 85→115 °C at 10 mm, model − data (min), mean | rate ratio 100/70 °C at 10 mm: model mean (data mean) | p_max (bar), max | p_v/p_sat, max | evaporation where p_v > p_sat |",
            )
            println(io, "|---|---|---|---|---|---|---|---|---|---|")
            for v in moist_variants
                n = length(summary(v, m -> 1))
                lr = summary(v, m -> m["loss_2p5h_mod"] / m["loss_2p5h_exp"])
                dc = summary(v, m -> m["cross_10_mod"] - m["cross_10_exp"])
                rm = summary(v, m -> m["ratio_10_mod"])
                re = summary(v, m -> m["ratio_10_exp"])
                @printf(
                    io,
                    "| %s | %d | %.1f | %.3f | %.2f (%.2f–%.2f) | %s | %s (%s) | %.2f | %.1f | %.0f %% |\n",
                    v,
                    n,
                    finite_median(summary(v, m -> m["rms_K"])),
                    finite_median(summary(v, m -> m["mean_nrms"])),
                    finite_mean(lr),
                    minimum(filter(isfinite, lr); init = Inf),
                    maximum(filter(isfinite, lr); init = -Inf),
                    fmt(finite_mean(dc); d = 1),
                    fmt(finite_mean(rm); d = 2),
                    fmt(finite_mean(re); d = 2),
                    maximum(summary(v, m -> m["p_max_bar"])),
                    maximum(summary(v, m -> m["pv_psat_max"])),
                    100 * finite_mean(summary(v, m -> m["evap_supersaturated_share"]))
                )
            end
        end
    end

    # The phase change: the summary of the model plan's columns on the long
    # runs and on all, against the gate; its runs and their cooling.
    pc_variants = filter(
        v ->
            startswith(v, "sorption") || v in ("gated_cal", "chan_cal", "arrhenius_cal"),
        moist_variants,
    )
    function pc_summary(sub, file, title)
        io = IOBuffer()
        println(
            io,
            "| variant | runs | rms (K), median | mass lost by 2.5 h, model/data, mean (range) | 85→115 °C at 10 mm, model − data (min), mean (runs) | rate ratio 100/70 °C at 10 mm, mean | p_max (bar), max | p_v/p_eq, max | X, max | T where it evaporates, p10 / p50 / p90 (°C), median | condensed / evaporated | accepted steps, total |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
        for v in pc_variants
            ms = [M[(v, n)] for n in sub if haskey(M, (v, n))]
            isempty(ms) && continue
            g(k) = [get(m, k, NaN) for m in ms]
            lr = g("loss_2p5h_mod") ./ g("loss_2p5h_exp")
            dc = filter(isfinite, g("cross_10_mod") .- g("cross_10_exp"))
            fl = filter(isfinite, lr)
            @printf(
                io,
                "| %s | %d | %.3f | %.4f (%.3f–%.3f) | %s (%d) | %s | %.3f | %s | %s | %s / %s / %s | %s | %s |\n",
                v,
                length(ms),
                finite_median(g("rms_K")),
                finite_mean(lr),
                minimum(fl; init = Inf),
                maximum(fl; init = -Inf),
                fmt(finite_mean(dc); d = 3),
                length(dc),
                fmt(finite_mean(g("ratio_10_mod")); d = 3),
                maximum(g("p_max_bar")),
                fmt(maximum(g("pv_peq_max")); d = 4),
                fmt(maximum(g("X_max")); d = 4),
                fmt(finite_median(g("Tevap_p10_C")); d = 1),
                fmt(finite_median(g("Tevap_p50_C")); d = 1),
                fmt(finite_median(g("Tevap_p90_C")); d = 1),
                fmt(sum(g("condensed_g")) / sum(g("evaporated_g")); d = 3),
                all(isfinite, g("accepted")) ? @sprintf("%d", sum(g("accepted"))) : "–"
            )
        end
        println(
            io,
            "\nData: rate ratio $(fmt(finite_mean([M[(pc_variants[1], n)]["ratio_10_exp"] for n in sub if haskey(M, (pc_variants[1], n))]); d = 3)).",
        )
        s = String(take!(io))
        write(joinpath(tables, file), s)
        println("\n$title\n", s)
    end
    if !isempty(pc_variants)
        pc_summary(long_runs, "phase_change_summary_long.md", "Phase change, long runs")
        pc_summary(moist_runs, "phase_change_summary.md", "Phase change, all moist runs")
        pc_shown = filter(
            v -> v in moist_variants,
            [
                "gated_cal",
                "sorption_cal_lambda5e-7",
                "sorption_cal",
                "sorption_cal_lambda5e-6",
                "sorption_cal_rh0",
                "sorption_cal_rh65",
                "sorption_cal_freewater",
            ],
        )
        law_tables(pc_shown, "_phase_change")
        for f in ("moist_runs_phase_change.md", "moist_cooling_phase_change.md")
            isfile(joinpath(tables, f)) || continue
            println("\n$f\n", read(joinpath(tables, f), String))
        end
        # each long run under the phase change: rms, mass ratio, transit, steps
        open(joinpath(tables, "phase_change_runs_long.md"), "w") do io
            sv = filter(v -> startswith(v, "sorption"), pc_variants)
            println(
                io,
                "| run | variant | rms (K) | mass ratio | 85→115 °C at 10 mm, model − data (min) | accepted steps |",
            )
            println(io, "|---|---|---|---|---|---|")
            for n in long_runs, v in sv
                haskey(M, (v, n)) || continue
                m = M[(v, n)]
                @printf(
                    io,
                    "| %s | %s | %.3f | %.4f | %s | %s |\n",
                    n,
                    v,
                    m["rms_K"],
                    m["loss_2p5h_mod"] / m["loss_2p5h_exp"],
                    fmt(m["cross_10_mod"] - m["cross_10_exp"]; d = 3),
                    haskey(m, "accepted") ? @sprintf("%d", m["accepted"]) : "–"
                )
            end
        end
    end

    # Figures: each moist run, data, the phase change and the three one-way laws
    for name in moist_runs
        r = runs[name]
        local S = Dict(
            l => read_series(series_path(l, name)) for
            l in base if isfile(series_path(l, name))
        )
        isempty(S) && continue
        local fig = Figure(size = (1300, 800))
        shown = [d for d in (2, 6, 10, 14, 20, 30, 40) if any(x -> mm(x) == d, r.depths)]
        styles = Dict(
            "sorption_cal" => :solid,
            "arrhenius_cal" => :dash,
            "gated_cal" => :dot,
            "chan_cal" => :dashdot,
        )
        for (i, d) in enumerate(shown)
            ax = Axis(
                fig[cld(i, 4), mod1(i, 4)];
                title = "$d mm",
                xlabel = "time (h)",
                ylabel = "T (°C)",
            )
            j = findfirst(x -> mm(x) == d, r.depths)
            lines!(ax, r.t ./ 3600, r.T[:, j] .- 273.15; color = :black, label = "data")
            for (l, s) in sort(collect(S); by = first)
                js = findfirst(x -> mm(x) == d, s.depths)
                lines!(
                    ax,
                    s.t ./ 3600,
                    s.T[:, js] .- 273.15;
                    linestyle = styles[l],
                    label = l,
                )
            end
            i == 1 && axislegend(ax; position = :rb, labelsize = 9)
        end
        ax = Axis(fig[2, 4]; title = "mass lost", xlabel = "time (h)", ylabel = "g")
        lines!(ax, r.t ./ 3600, r.mass[1] .- r.mass; color = :black)
        hlines!(ax, [initial_water_g(r)]; color = :gray, linestyle = :dot)
        for (l, s) in sort(collect(S); by = first)
            lines!(ax, s.t ./ 3600, s.mass[1] .- s.mass; linestyle = styles[l])
        end
        Label(fig[0, :], name; fontsize = 16)
        save(joinpath(figdir, "moist_$(replace(name, ',' => 'p', '%' => "pc")).png"), fig)
    end
end
println(
    "compare: $(length(dry_variants)) dry variants, $(length(moist_variants)) moist variants",
)
