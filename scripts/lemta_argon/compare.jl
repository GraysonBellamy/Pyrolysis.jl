# Score every variant run against the data: per-test metrics
# (results/metrics_*.csv), summary tables in Markdown (results/tables/*.md) and
# figures (results/figures/compare_<species>.png, levers.png).
#
# - Mass: m/m₀ at 300, 600 and 698 s; the mean mass loss rate over an early
#   window (from the balance's settling, at least 30 s, to 300 s) and a late
#   one (300–698 s), from the change of m/m₀ and the test's m₀; NRMS of m/m₀
#   over the record after settling, over the measured loss 1 − min m/m₀.
# - Temperature: per thermocouple, on its usable record (cuts applied), the
#   rms, the bias (model − data) and NRMS = rms over the measured range; the
#   bias over the first 120 s; the readings at 300 and 600 s.
# - Char depth (I9CL2R inert rows at 12 min, the species and fluxes the
#   QRYAJ9 mass tests share): the model's char/virgin boundary at the
#   conversion thresholds CHAR_THRESHOLDS and its thickness lost, at 720 s.
#
#   julia --project=scripts/verification scripts/lemta_argon/compare.jl
include(joinpath(@__DIR__, "LemtaArgon.jl"))
using .LemtaArgon
include(joinpath(@__DIR__, "variants.jl"))
using CairoMakie
using Printf
using Statistics

const AREA = 0.01                # m², exposed face
const EARLY = 120.0              # s, the early-temperature window
const NEAR = (2, 4, 6)           # mm, depth groups of the summaries
const MID = (8, 10, 12, 14)
const DEEP = (16, 18, 20)
const BACK = (30, 40)
const COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]

tables = joinpath(RESULTS, "tables")
figdir = joinpath(RESULTS, "figures")
mkpath(tables)
mkpath(figdir)

mass = Dict(r.name => r for r in mass_runs())
temps = Dict(r.name => r for r in temperature_runs())

present = sort([v for v in keys(VARIANTS) if isdir(joinpath(RESULTS, v))])
"deck" in present || error("no deck results: run run_cases.jl deck first")
series = Dict(
    v => Dict(
        splitext(f)[1] => read_series(joinpath(RESULTS, v, f)) for
        f in readdir(joinpath(RESULTS, v)) if endswith(f, ".csv")
    ) for v in present
)
failures = Dict(
    v => [splitext(f)[1] for f in readdir(joinpath(RESULTS, v)) if endswith(f, ".failed")] for v in present
)

fmt(x; d = 3) = isfinite(x) ? @sprintf("%.*f", d, x) : "–"
mm(d) = round(Int, 1000d)
flux(r) = round(Int, r.flux)
nanmean(x) = (y = filter(isfinite, x); isempty(y) ? NaN : mean(y))

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

# The record's mean within ±5 s of t0 (the balance is noisy at 1 Hz).
data_at(t, y, t0) = (k = abs.(t .- t0) .<= 5; any(k) ? mean(y[k]) : NaN)

mlr_mean(m0, ma, mb, ta, tb) = m0 * (ma - mb) / (AREA * (tb - ta))   # g/(m² s)

function mass_metrics(r::MassRun, s)
    md(t0) = data_at(r.t, r.m, t0)
    mm_(t0) = interp(s.t, s.m, t0)
    t1 = max(r.t_settle, 30.0)
    keep = findall(t -> r.t_settle <= t <= min(698.0, s.t[end]), r.t)
    e = [interp(s.t, s.m, t) for t in r.t[keep]] .- r.m[keep]
    return (;
        m300 = (mm_(300.0), md(300.0)),
        m600 = (mm_(600.0), md(600.0)),
        m698 = (mm_(698.0), md(698.0)),
        mlr_early = (
            mlr_mean(r.m0, mm_(t1), mm_(300.0), t1, 300.0),
            mlr_mean(r.m0, md(t1), md(300.0), t1, 300.0),
        ),
        mlr_late = (
            mlr_mean(r.m0, mm_(300.0), mm_(698.0), 300.0, 698.0),
            mlr_mean(r.m0, md(300.0), md(698.0), 300.0, 698.0),
        ),
        nrms = sqrt(mean(e .^ 2)) / (1 - minimum(r.m[keep])),
    )
end

function temperature_metrics(r::TemperatureRun, s)
    out = Dict{Int,NamedTuple}()
    for (j, d) in enumerate(DEPTHS)
        live(r, j) || continue
        js = findfirst(x -> isapprox(x, d), s.depths)
        js === nothing && continue
        y = usable(r, j)
        keep = findall(k -> isfinite(y[k]) && r.t[k] <= s.t[end], eachindex(y))
        length(keep) < 10 && continue
        m = [interp(s.t, view(s.T, :, js), t) for t in r.t[keep]]
        e = m .- y[keep]
        early = findall(k -> r.t[keep[k]] <= EARLY, eachindex(keep))
        at(t0) = (
            interp(s.t, view(s.T, :, js), t0),
            (
                k = findall(i -> abs(r.t[i] - t0) <= 2 && isfinite(y[i]), eachindex(y));
                isempty(k) || t0 > r.t_end ? NaN : mean(y[k])
            ),
        )
        out[mm(d)] = (
            n = length(keep),
            rms = sqrt(mean(e .^ 2)),
            bias = mean(e),
            nrms = sqrt(mean(e .^ 2)) / (maximum(y[keep]) - minimum(y[keep])),
            bias_early = isempty(early) ? NaN : mean(e[early]),
            T300 = at(300.0),
            T600 = at(600.0),
        )
    end
    return out
end

group_bias(tm, group; field = :bias) =
    nanmean([getfield(tm[d], field) for d in group if haskey(tm, d)])

# ---------------------------------------------------------------------------
# Per-test metrics files
# ---------------------------------------------------------------------------

MM = Dict{Tuple{String,String},NamedTuple}()     # (variant, test) => mass metrics
TM = Dict{Tuple{String,String},Dict{Int,NamedTuple}}()
for v in present, (name, s) in series[v]
    if haskey(mass, name)
        MM[(v, name)] = mass_metrics(mass[name], s)
    elseif haskey(temps, name)
        TM[(v, name)] = temperature_metrics(temps[name], s)
    end
end

open(joinpath(RESULTS, "metrics_mass.csv"), "w") do io
    println(
        io,
        "variant,test,species,flux,rho,t_settle,m300_model,m300_data,m600_model,m600_data,m698_model,m698_data,mlr_early_model,mlr_early_data,mlr_late_model,mlr_late_data,nrms_m,retcode,window_retcodes,mass_relative,energy_relative",
    )
    for ((v, name), x) in sort(collect(MM); by = first)
        r, s = mass[name], series[v][name]
        @printf(
            io,
            "%s,%s,%s,%d,%g,%g,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.3f,%.3f,%.4f,%s,%s,%.2g,%.2g\n",
            v,
            name,
            r.species,
            flux(r),
            r.rho,
            r.t_settle,
            x.m300...,
            x.m600...,
            x.m698...,
            x.mlr_early...,
            x.mlr_late...,
            x.nrms,
            s.retcode,
            join(s.window_retcodes, ";"),
            s.mass_relative,
            s.energy_relative
        )
    end
end

open(joinpath(RESULTS, "metrics_temperature.csv"), "w") do io
    println(
        io,
        "variant,test,species,flux,rho,depth_mm,n,rms_K,bias_K,nrms,bias_first120s_K,T300_model_K,T300_data_K,T600_model_K,T600_data_K",
    )
    for ((v, name), tm) in sort(collect(TM); by = first), d in sort(collect(keys(tm)))
        r, x = temps[name], tm[d]
        @printf(
            io,
            "%s,%s,%s,%d,%g,%d,%d,%.2f,%.2f,%.4f,%.2f,%.1f,%.1f,%.1f,%.1f\n",
            v,
            name,
            r.species,
            flux(r),
            r.rho,
            d,
            x.n,
            x.rms,
            x.bias,
            x.nrms,
            x.bias_early,
            x.T300...,
            x.T600...
        )
    end
end

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

species_order = ["balsa", "poplar", "fir", "larch", "oak", "beech", "ipe"]
# The variant's test of a species at a flux (one per species and flux), or nothing.
function _test(runs, metrics, v, sp, q)
    ns =
        [n for (n, r) in runs if r.species == sp && flux(r) == q && haskey(metrics, (v, n))]
    return isempty(ns) ? nothing : only(ns)
end
mname(v, sp, q) = _test(mass, MM, v, sp, q)
tname(v, sp, q) = _test(temps, TM, v, sp, q)

# Solves: every series' return codes and every failure.
open(joinpath(tables, "solves.md"), "w") do io
    println(
        io,
        "| variant | tests | success | not success | failed (threw) | max mass / energy defect |",
    )
    println(io, "|---|---|---|---|---|---|")
    for v in present
        ss = collect(values(series[v]))
        ok = count(s -> s.retcode == "Success" && all(==("Success"), s.window_retcodes), ss)
        bad = [
            n for (n, s) in series[v] if
            !(s.retcode == "Success" && all(==("Success"), s.window_retcodes))
        ]
        @printf(
            io,
            "| %s | %d | %d | %s | %s | %.1e / %.1e |\n",
            v,
            length(ss) + length(failures[v]),
            ok,
            isempty(bad) ? "–" : join(bad, ", "),
            isempty(failures[v]) ? "–" : join(failures[v], ", "),
            maximum(s -> abs(s.mass_relative), ss; init = 0.0),
            maximum(s -> abs(s.energy_relative), ss; init = 0.0)
        )
    end
end

# The decks: mass by species and flux.
open(joinpath(tables, "deck_mass.md"), "w") do io
    println(
        io,
        "| species | flux | ρ | m/m₀ 300 s model / data | m/m₀ 698 s model / data | rate 30–300 s model / data | rate 300–698 s model / data | NRMS m/m₀ |",
    )
    println(io, "|---|---|---|---|---|---|---|---|")
    for sp in species_order, q in (25, 35, 45)
        n = mname("deck", sp, q)
        n === nothing && continue
        x = MM[("deck", n)]
        @printf(
            io,
            "| %s | %d | %g | %.3f / %.3f | %.3f / %.3f | %.2f / %.2f | %.2f / %.2f (%.0f %%) | %.2f |\n",
            sp,
            q,
            mass[n].rho,
            x.m300...,
            x.m698...,
            x.mlr_early...,
            x.mlr_late...,
            100x.mlr_late[1] / x.mlr_late[2],
            x.nrms
        )
    end
end

# The decks: temperature by test, per depth group.
open(joinpath(tables, "deck_temperature.md"), "w") do io
    println(
        io,
        "| species | flux | ρ | NRMS 2–20 mm | bias 2–6 | bias 8–14 | bias 16–20 | bias 30–40 | bias 2–4 mm, first 120 s | 10 mm at 600 s model / data |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    for sp in species_order, q in (35, 45)
        n = tname("deck", sp, q)
        n === nothing && continue
        tm = TM[("deck", n)]
        t10 = haskey(tm, 10) ? tm[10].T600 : (NaN, NaN)
        @printf(
            io,
            "| %s | %d | %g | %.3f | %s | %s | %s | %s | %s | %s / %s |\n",
            sp,
            q,
            temps[n].rho,
            nanmean([tm[d].nrms for d in (NEAR..., MID..., DEEP...) if haskey(tm, d)]),
            fmt(group_bias(tm, NEAR); d = 0),
            fmt(group_bias(tm, MID); d = 0),
            fmt(group_bias(tm, DEEP); d = 0),
            fmt(group_bias(tm, BACK); d = 0),
            fmt(group_bias(tm, (2, 4); field = :bias_early); d = 0),
            fmt(t10[1] - 273.15; d = 0),
            fmt(t10[2] - 273.15; d = 0)
        )
    end
end

# Levers: each variant against the decks on the same tests at 35–45 kW/m².
lever_rows = []
for v in present
    v == "deck" && continue
    ms = [n for (vv, n) in keys(MM) if vv == v && flux(mass[n]) >= 35]
    ts = [n for (vv, n) in keys(TM) if vv == v]
    dm = [
        MM[(v, n)].m698[1] - MM[("deck", n)].m698[1] for n in ms if haskey(MM, ("deck", n))
    ]
    rl(w, n) = MM[(w, n)].mlr_late[1] / MM[(w, n)].mlr_late[2]
    rate = [rl(v, n) for n in ms if haskey(MM, ("deck", n))]
    rate_deck = [rl("deck", n) for n in ms if haskey(MM, ("deck", n))]
    gb(w, g; field = :bias) =
        nanmean([group_bias(TM[(w, n)], g; field) for n in ts if haskey(TM, ("deck", n))])
    push!(
        lever_rows,
        (;
            v,
            label = VARIANTS[v].label,
            kind = VARIANTS[v].kind,
            nm = length(rate),
            nt = length(ts),
            dm = nanmean(dm),
            rate = nanmean(rate),
            rate_deck = nanmean(rate_deck),
            near = (gb(v, NEAR), gb("deck", NEAR)),
            mid = (gb(v, MID), gb("deck", MID)),
            deep = (gb(v, DEEP), gb("deck", DEEP)),
            early = (
                gb(v, (2, 4); field = :bias_early),
                gb("deck", (2, 4); field = :bias_early),
            ),
        ),
    )
end
open(joinpath(tables, "levers.md"), "w") do io
    println(
        io,
        "Each variant against the decks on the tests it applies to at 35 and 45 kW/m² (means over those tests; temperature biases model − data in K, pooled over each test's depth group; the deck's own value on the same tests in brackets).\n",
    )
    println(
        io,
        "| variant | mass tests | temperature tests | Δ m/m₀ at 698 s | late rate, model/data | bias 2–4 mm first 120 s | bias 2–6 mm | bias 8–14 mm | bias 16–20 mm |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for x in sort(lever_rows; by = x -> (x.kind, x.v))
        pair(p) = @sprintf("%s (%s)", fmt(p[1]; d = 0), fmt(p[2]; d = 0))
        @printf(
            io,
            "| %s | %d | %d | %s | %s (%s) | %s | %s | %s | %s |\n",
            x.label,
            x.nm,
            x.nt,
            fmt(x.dm),
            fmt(x.rate; d = 2),
            fmt(x.rate_deck; d = 2),
            pair(x.early),
            pair(x.near),
            pair(x.mid),
            pair(x.deep)
        )
    end
end

# Levers by species: each variant's own value (not its change) of three
# metrics at 35 and 45 kW/m², the mean of the two fluxes' tests; the decks first.
function by_species(v, sp, f)
    ns = [n for q in (35, 45) for n in (mname(v, sp, q), tname(v, sp, q)) if n !== nothing]
    return nanmean([f(n) for n in ns])
end
late_rate(v) =
    n -> haskey(MM, (v, n)) ? MM[(v, n)].mlr_late[1] / MM[(v, n)].mlr_late[2] : NaN
tbias(v, g; field = :bias) =
    n -> haskey(TM, (v, n)) ? group_bias(TM[(v, n)], g; field) : NaN
open(joinpath(tables, "levers_by_species.md"), "w") do io
    order = sort(
        [v for v in present if VARIANTS[v].kind !== :convergence];
        by = v -> (v != "deck", VARIANTS[v].kind, v),
    )
    for (title, metric, d) in (
        ("Late mass loss rate (300–698 s), model over data", late_rate, 2),
        ("Bias 8–14 mm (K), model − data", v -> tbias(v, MID), 0),
        ("Bias 16–20 mm (K)", v -> tbias(v, DEEP), 0),
        ("Bias 2–6 mm (K), usable record", v -> tbias(v, NEAR), 0),
        (
            "Bias 2–4 mm over the first 120 s (K)",
            v -> tbias(v, (2, 4); field = :bias_early),
            0,
        ),
    )
        println(io, "\n**$title**, mean of the 35 and 45 kW/m² tests\n")
        println(io, "| variant | ", join(species_order, " | "), " |")
        println(io, "|---|", repeat("---|", length(species_order)))
        for v in order
            vals = [by_species(v, sp, metric(v)) for sp in species_order]
            println(io, "| ", VARIANTS[v].label, " | ", join(fmt.(vals; d), " | "), " |")
        end
    end
end

# Convergence: the largest change from the decks at any thermocouple and time.
open(joinpath(tables, "convergence.md"), "w") do io
    println(
        io,
        "| variant | test | max |ΔT| at the thermocouples (K) | max |Δ m/m₀| | Δ char depth X 0.5 at the end (mm) |",
    )
    println(io, "|---|---|---|---|---|")
    for v in present
        VARIANTS[v].kind === :convergence || continue
        for (name, s) in sort(collect(series[v]); by = first)
            haskey(series["deck"], name) || continue
            b = series["deck"][name]
            dT = maximum(
                abs(interp(s.t, view(s.T, :, j), b.t[k]) - b.T[k, j]) for
                k in eachindex(b.t), j in eachindex(b.depths);
                init = 0.0,
            )
            dm = maximum(abs(interp(s.t, s.m, t) - b.m[k]) for (k, t) in enumerate(b.t))
            dc = 1000 * (s.char[end, 2] - b.char[end, 2])
            @printf(io, "| %s | %s | %.2f | %.1e | %+.2f |\n", v, name, dT, dm, dc)
        end
    end
end

# Char depth against I9CL2R at 720 s.
rows = char_depth_rows()
open(joinpath(RESULTS, "metrics_char.csv"), "w") do io
    println(
        io,
        "variant,test,species,flux,char_X20_mm,char_X50_mm,char_X80_mm,lost_mm,measured_char_mm,measured_lost_mm",
    )
    for v in present, (name, s) in sort(collect(series[v]); by = first)
        haskey(mass, name) || continue
        r = mass[name]
        key = (r.species, flux(r))
        haskey(rows, key) || continue
        k = findlast(<=(CHAR_EXPOSURE + 1e-6), s.t)
        @printf(
            io,
            "%s,%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%s,%s\n",
            v,
            name,
            r.species,
            flux(r),
            (1000 .* s.char[k, :])...,
            1000 * (s.thickness[1] - s.thickness[k]),
            join([@sprintf("%.2f", c) for (c, _) in rows[key]], ";"),
            join([@sprintf("%.2f", l) for (_, l) in rows[key]], ";")
        )
    end
end
open(joinpath(tables, "char_depth.md"), "w") do io
    println(
        io,
        "Char depth and thickness lost at 720 s: the decks (char/virgin boundary at conversion X = 0.5, with 0.2–0.8 in brackets) against the I9CL2R inert rows (each row; mm).\n",
    )
    println(
        io,
        "| species | flux | model char X 0.5 (0.2–0.8) | measured char | model lost | measured lost | model converted (char + lost) | measured converted |",
    )
    println(io, "|---|---|---|---|---|---|---|---|")
    for sp in species_order, q in (25, 35, 45)
        haskey(rows, (sp, q)) || continue
        n = mname("deck", sp, q)
        n === nothing && continue
        s = series["deck"][n]
        k = findlast(<=(CHAR_EXPOSURE + 1e-6), s.t)
        c = 1000 .* s.char[k, :]
        lost = 1000 * (s.thickness[1] - s.thickness[k])
        meas = rows[(sp, q)]
        @printf(
            io,
            "| %s | %d | %.1f (%.1f–%.1f) | %s | %.1f | %s | %.1f | %s |\n",
            sp,
            q,
            c[2],
            c[1],
            c[3],
            join([@sprintf("%.1f", x) for (x, _) in meas], ", "),
            lost,
            join([@sprintf("%.1f", x) for (_, x) in meas], ", "),
            c[2] + lost,
            join([@sprintf("%.1f", a + b) for (a, b) in meas], ", ")
        )
    end
end

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

for sp in species_order
    fig = Figure(size = (1300, 900))
    ax1 =
        Axis(fig[1, 1]; title = "$sp: m/m₀ (data solid, decks dashed)", xlabel = "time (s)")
    ax2 = Axis(fig[1, 2]; title = "mass loss rate (g/(m² s))", xlabel = "time (s)")
    for (i, q) in enumerate((25, 35, 45))
        n = mname("deck", sp, q)
        n === nothing && continue
        r, s = mass[n], series["deck"][n]
        lines!(ax1, r.t, r.m; color = COLORS[i], linewidth = 2, label = "$q kW/m²")
        lines!(ax1, s.t, s.m; color = COLORS[i], linewidth = 2, linestyle = :dash)
        lines!(ax2, r.t, r.mlr; color = (COLORS[i], 0.6), linewidth = 1)
        lines!(ax2, s.t, s.mlr; color = COLORS[i], linewidth = 2, linestyle = :dash)
    end
    axislegend(ax1; position = :lb)
    for (i, q) in enumerate((35, 45))
        n = tname("deck", sp, q)
        n === nothing && continue
        r, s = temps[n], series["deck"][n]
        ax = Axis(
            fig[2, i];
            title = "$sp $q kW/m², ρ $(r.rho)",
            xlabel = "time (s)",
            ylabel = "T (°C)",
        )
        for (k, d) in enumerate((0.002, 0.006, 0.010, 0.014, 0.020, 0.040))
            j = findfirst(x -> isapprox(x, d), DEPTHS)
            js = findfirst(x -> isapprox(x, d), s.depths)
            if live(r, j)
                lines!(
                    ax,
                    r.t,
                    usable(r, j) .- 273.15;
                    color = COLORS[k],
                    linewidth = 2,
                    label = "$(mm(d)) mm",
                )
            end
            js === nothing || lines!(
                ax,
                s.t,
                s.T[:, js] .- 273.15;
                color = COLORS[k],
                linewidth = 2,
                linestyle = :dash,
            )
        end
        axislegend(ax; position = :lt)
    end
    save(joinpath(figdir, "compare_$sp.png"), fig)
end

levers = filter(x -> x.kind in (:lever, :diagnostic), lever_rows)
if !isempty(levers)
    sort!(levers; by = x -> x.label)
    fig = Figure(size = (1300, 60 + 34length(levers)))
    labels = [x.label for x in levers]
    ys = 1:length(levers)
    for (i, (title, f)) in enumerate((
        ("Δ late mass rate, model/data (35–45 kW/m²)", x -> x.rate - x.rate_deck),
        ("Δ bias 8–14 mm (K)", x -> x.mid[1] - x.mid[2]),
        ("Δ bias 2–4 mm, first 120 s (K)", x -> x.early[1] - x.early[2]),
    ))
        ax = Axis(
            fig[1, i];
            title,
            yticks = (ys, i == 1 ? labels : fill("", length(ys))),
            yreversed = true,
        )
        vals = [f(x) for x in levers]
        colors = [x.kind === :diagnostic ? "#9a9a9a" : COLORS[1] for x in levers]
        barplot!(
            ax,
            ys,
            replace(vals, NaN => 0.0);
            direction = :x,
            color = colors,
            gap = 0.3,
        )
        vlines!(ax, [0.0]; color = :gray50, linewidth = 1)
    end
    save(joinpath(figdir, "levers.png"), fig)
end

println("wrote metrics_*.csv, tables/*.md and figures/compare_*.png under $RESULTS")
