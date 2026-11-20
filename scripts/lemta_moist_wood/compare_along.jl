# Score the along-grain variant runs against the data: per-run metrics
# (results/metrics_along.csv), tables (results/tables/along_*.md, also
# printed) and figures (results/figures/along_*.png).
#
# The mass is compared as the share of the initial water lost: the model block
# starts at the stated moisture content over the deck's density
# (initial_water_g), the measurement at MC/(1 + MC) × m₀ (measured_water_g),
# and the two disagree by 4–27 % (README.md). The runs without thermocouples
# are scored on the mass alone.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/compare_along.jl
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using CairoMakie
using Printf
using Statistics

runs = Dict(r.name => r for r in load_runs(; grain = :along))
tables = joinpath(RESULTS, "tables")
figdir = joinpath(RESULTS, "figures")
mkpath(tables)
mkpath(figdir)

fmt(x; d = 3) = isfinite(x) ? @sprintf("%.*f", d, x) : "–"
mm(d) = round(Int, 1000d)
"""The variant's along-grain series, by run name."""
series_of(v) = Dict(
    n => read_series(series_path(v, n)) for n in keys(runs) if isfile(series_path(v, n))
)
present(group) = sort([
    v for v in keys(VARIANTS) if
    VARIANTS[v].runs == group && isdir(joinpath(RESULTS, v)) && !isempty(series_of(v))
])
function report(name, text)
    write(joinpath(tables, name), text)
    println("\n$name\n", text)
end

# ---------------------------------------------------------------------------
# Dry runs
# ---------------------------------------------------------------------------

io = IOBuffer()
println(
    io,
    "| variant | run | rms (K) | within TC accuracy | bias 2 / 6 / 10 / 14 / 20 / 30 / 40 mm (K) |",
)
println(io, "|---|---|---|---|---|")
for v in present(:dry), (name, s) in sort(collect(series_of(v)); by = first)
    e = temperature_errors(runs[name], s)
    b(d) = haskey(e.per, d) ? @sprintf("%+.1f", e.per[d].bias) : "–"
    @printf(
        io,
        "| %s | %s | %.2f | %.0f %% | %s |\n",
        v,
        name,
        e.rms,
        100e.within,
        join(b.((2, 6, 10, 14, 20, 30, 40)), " / ")
    )
end
report("along_dry.md", String(take!(io)))

# ---------------------------------------------------------------------------
# Moist runs
# ---------------------------------------------------------------------------

"""The moist-run metrics of one series against its run."""
function moist_metrics(r, s)
    t25 = min(9000.0, r.t[end], s.t[end])
    tend = min(r.t_mass[end], s.t[end])
    e = temperature_errors(r, s)
    out = Dict{String,Float64}(
        "rms_K" => e.rms,
        "rms_heating_K" => temperature_errors(r, s; window = (0.0, t25)).rms,
        "water_exp_g" => measured_water_g(r),
        "water_mod_g" => s.water[1],
        "lost_2p5h_exp" => water_fraction_lost(r, t25),
        "lost_2p5h_mod" => water_fraction_lost(s, t25),
        "lost_end_exp" => water_fraction_lost(r, tend),
        "lost_end_mod" => water_fraction_lost(s, tend),
        "p_max_bar" => maximum(s.p_max) / 1e5,
    )
    # after the withdrawal: the change in the share lost (positive: still losing)
    out["after_off_exp"] = out["lost_end_exp"] - water_fraction_lost(r, r.t_off)
    out["after_off_mod"] = out["lost_end_mod"] - water_fraction_lost(s, r.t_off)
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
        out["bias_$(d)"] = x.bias
        out["rms_$(d)"] = x.rms
    end
    pc = s.phase_change
    if pc !== nothing
        out["evaporated_g"] = trapz(s.t, s.evaporation)
        out["condensed_g"] = trapz(s.t, pc.condensation)
        out["pv_peq_max"] = maximum(pc.pv_over_peq)
        out["X_max"] = maximum(pc.X_max)
        out["accepted"] = pc.accepted
        for (q, T) in zip((10, 50, 90), pc.T_evap)
            out["Tevap_p$(q)_C"] = T - 273.15
        end
    end
    return out
end

moist = present(:moist)
M = Dict{Tuple{String,String},Dict{String,Float64}}()
for v in moist, (name, s) in series_of(v)
    M[(v, name)] = moist_metrics(runs[name], s)
end
names = sort(unique(n for (_, n) in keys(M)))
keys_out = sort(unique(k for m in values(M) for k in keys(m)))
open(joinpath(RESULTS, "metrics_along.csv"), "w") do io
    println(io, join(["variant", "run", keys_out...], ","))
    for ((v, n), m) in sort(collect(M); by = first)
        println(
            io,
            join(
                [v, "\"$n\"", (@sprintf("%.7g", get(m, k, NaN)) for k in keys_out)...],
                ",",
            ),
        )
    end
end

finite(x) = filter(isfinite, x)
fmean(x) = isempty(finite(x)) ? NaN : mean(finite(x))

# Summary per variant: the runs with thermocouples for the temperatures, all
# for the mass
io = IOBuffer()
println(
    io,
    "| variant | runs | rms T (K), mean over the TC runs | water lost by 2.5 h, model/data, mean (range) | 85→115 °C at 10 mm, model − data (min), mean | rate ratio 100/70 °C at 10 mm: model (data) | after the withdrawal, share lost: model (data), mean | p_max (bar), max | X, max | condensed / evaporated |",
)
println(io, "|---|---|---|---|---|---|---|---|---|---|")
for v in moist
    ms = [M[(v, n)] for n in names if haskey(M, (v, n))]
    g(k) = [get(m, k, NaN) for m in ms]
    lr = finite(g("lost_2p5h_mod") ./ g("lost_2p5h_exp"))
    @printf(
        io,
        "| %s | %d | %s | %s (%s–%s) | %s | %s (%s) | %s (%s) | %.2f | %s | %s |\n",
        v,
        length(ms),
        fmt(fmean(g("rms_K")); d = 2),
        fmt(fmean(lr); d = 2),
        fmt(minimum(lr; init = Inf); d = 2),
        fmt(maximum(lr; init = -Inf); d = 2),
        fmt(fmean(g("cross_10_mod") .- g("cross_10_exp")); d = 1),
        fmt(fmean(g("ratio_10_mod")); d = 2),
        fmt(fmean(g("ratio_10_exp")); d = 2),
        fmt(fmean(g("after_off_mod")); d = 3),
        fmt(fmean(g("after_off_exp")); d = 3),
        maximum(g("p_max_bar")),
        fmt(maximum(finite(g("X_max")); init = -Inf); d = 3),
        fmt(sum(finite(g("condensed_g"))) / sum(finite(g("evaporated_g"))); d = 3)
    )
end
report("along_summary.md", String(take!(io)))

# Each run under the phase change and the gate
shown = filter(
    v -> v in moist,
    ["sorption_cal", "sorption", "gated_cal", "arrhenius_cal", "chan_cal"],
)
io = IOBuffer()
println(
    io,
    "| run | water (g): data / model | lost by 2.5 h: data / " *
    join(shown, " / ") *
    " | after the withdrawal: data / " *
    join(shown, " / ") *
    " | rms (K): " *
    join(shown, " / ") *
    " | 85→115 °C at 10 mm (min): data / " *
    join(shown, " / ") *
    " |",
)
println(io, "|---|---|---|---|---|---|")
for n in names
    haskey(M, (shown[1], n)) || continue
    m1 = M[(shown[1], n)]
    g(l, k, d) = haskey(M, (l, n)) ? fmt(M[(l, n)][k]; d) : "–"
    println(
        io,
        "| $n | $(fmt(m1["water_exp_g"]; d = 1)) / $(fmt(m1["water_mod_g"]; d = 1)) | " *
        "$(fmt(m1["lost_2p5h_exp"]; d = 2)) / " *
        join([g(l, "lost_2p5h_mod", 2) for l in shown], " / ") *
        " | $(fmt(m1["after_off_exp"]; d = 3)) / " *
        join([g(l, "after_off_mod", 3) for l in shown], " / ") *
        " | " *
        join([g(l, "rms_K", 1) for l in shown], " / ") *
        " | $(fmt(m1["cross_10_exp"]; d = 0)) / " *
        join([g(l, "cross_10_mod", 0) for l in shown], " / ") *
        " |",
    )
end
report("along_runs.md", String(take!(io)))

# Bias by depth over the heating and cooling record, mean over the TC runs
depths = (2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40)
io = IOBuffer()
println(io, "| variant | " * join(["$d mm" for d in depths], " | ") * " |")
println(io, "|---|" * repeat("---|", length(depths)))
for v in moist
    cells = [
        fmt(
            fmean([
                M[(v, n)]["bias_$d"] for
                n in names if haskey(M, (v, n)) && haskey(M[(v, n)], "bias_$d")
            ]);
            d = 1,
        ) for d in depths
    ]
    println(io, "| $v | " * join(cells, " | ") * " |")
end
report("along_bias_by_depth.md", String(take!(io)))

# Figures: each run, data against the phase change and the gate
styles = Dict("sorption_cal" => :solid, "gated_cal" => :dot, "sorption" => :dash)
for v in present(:dry)
    S = series_of(v)
    fig = Figure(size = (1200, 330 * length(S)))
    for (i, (name, s)) in enumerate(sort(collect(S); by = first))
        r = runs[name]
        ax = Axis(fig[i, 1]; title = "$v: $name", xlabel = "time (h)", ylabel = "T (°C)")
        for (j, d) in enumerate(r.depths)
            mm(d) in (2, 6, 10, 14, 20, 30, 40) || continue
            js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
            c = Makie.wong_colors()[mod1(j, 7)]
            lines!(ax, r.t ./ 3600, r.T[:, j] .- 273.15; color = c)
            lines!(ax, s.t ./ 3600, s.T[:, js] .- 273.15; color = c, linestyle = :dash)
        end
    end
    save(joinpath(figdir, "along_$v.png"), fig)
end
for name in names
    r = runs[name]
    S = Dict(
        l => read_series(series_path(l, name)) for
        l in keys(styles) if isfile(series_path(l, name))
    )
    isempty(S) && continue
    fig = Figure(size = (1300, 800))
    shown_d = [d for d in (2, 6, 10, 14, 20, 30, 40) if any(x -> mm(x) == d, r.depths)]
    for (i, d) in enumerate(shown_d)
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
            lines!(ax, s.t ./ 3600, s.T[:, js] .- 273.15; linestyle = styles[l], label = l)
        end
        i == 1 && axislegend(ax; position = :rb, labelsize = 9)
    end
    ax = Axis(fig[2, 4]; title = "share of the water lost", xlabel = "time (h)")
    lines!(
        ax,
        r.t_mass ./ 3600,
        (r.mass[1] .- r.mass) ./ measured_water_g(r);
        color = :black,
    )
    for (l, s) in sort(collect(S); by = first)
        lines!(ax, s.t ./ 3600, (s.mass[1] .- s.mass) ./ s.water[1]; linestyle = styles[l])
    end
    Label(fig[0, :], name; fontsize = 16)
    save(joinpath(figdir, "along_$(replace(name, ',' => 'p', '%' => "pc")).png"), fig)
end
println(
    "compare_along: $(length(present(:dry))) dry variants, $(length(moist)) moist variants",
)
