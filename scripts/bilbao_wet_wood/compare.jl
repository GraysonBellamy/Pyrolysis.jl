# Score the solved variants against the digitised data: temperature rms by
# depth and the total conversion within each run's window, the plateau and
# crossing times at 9 and 19 mm, the 9 mm reading when the measured 1 mm
# record passes 440 and 530 K, the phase change's diagnostics, and figures.
#   julia --project=scripts/verification scripts/bilbao_wet_wood/compare.jl [FAMILY ...]
# Reads results/<family>_<law>[_tag]/<run>.csv; writes results/metrics.csv,
# results/tables/*.md and results/figures/*.png, and prints the tables.
include(joinpath(@__DIR__, "BilbaoWetWood.jl"))
using .BilbaoWetWood
using CairoMakie
using Printf
using Statistics

const LAWS_SHOWN = ("sorption", "gated", "chan", "arrhenius")
const LEVEL = 626.0   # K: the dry run's 19 mm reading at 24 min, the end of the calibration window
const MARKS = (440.0, 530.0)   # K: levels of the measured 1 mm record at which the 9 mm reading is compared
const P_ATM = 101325.0

"""
The scoring window of a run [s]: from 0 to the first time its measured 19 mm
record exceeds `LEVEL` (the thermal model is calibrated below it), or the end
of its record.
"""
function window(run)
    meas = temperatures(run)
    t, T = meas[19]
    k = findfirst(>(LEVEL), T)
    return (0.0, k === nothing ? t[end] : t[k])
end

"""Measured (t [s], T [K]) at `depth`, or empty."""
measured(run, depth) = get(temperatures(run), depth, (Float64[], Float64[]))

model_T(s, depth) = depth == 1 ? s.T1 : depth == 9 ? s.T9 : s.T19

function xt_rms(run, s; window)
    t, X = conversion(run)
    keep = findall(x -> window[1] <= x <= min(window[2], s.t[end]), t)
    isempty(keep) && return NaN
    return sqrt(mean((interp(s.t, s.X_T, t[i]) - X[i])^2 for i in keep))
end

# Minutes at which a series first reaches `level`, from a (t [s], T [K]) pair.
cross(t, T, level) = crossing_minutes(t, T, level)

trapz(t, y) = sum(0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))

# A key=value field of the series' `# retcode` line, as a number.
meta_value(s, key) =
    (m = match(Regex("$key=(\\S+)"), s.meta); m === nothing ? NaN : parse(Float64, m[1]))

function metrics(variant, run)
    path = series_path(variant, run)
    isfile(path) || return nothing
    s = read_series(path)
    w = window(run)
    rms = rms_by_depth(run, s; window = w)
    row = Dict{String,Any}(
        "variant" => variant,
        "run" => run,
        "window_min" => w[2] / 60,
        "rms_1mm" => get(rms, 1, NaN),
        "rms_9mm" => get(rms, 9, NaN),
        "rms_19mm" => get(rms, 19, NaN),
        "rms_XT" => xt_rms(run, s; window = w),
        "Tmin_19mm" => minimum(s.T19),
        "p_max_kPa" => (maximum(s.p_max) - P_ATM) / 1000,
        "mass_closure" => meta_value(s, "mass_closure"),
        "energy_closure" => meta_value(s, "energy_closure"),
        "meta" => s.meta,
    )
    for d in (9, 19), level in (373.0, 400.0, 440.0, 500.0)
        row[@sprintf("t%d_%dmm", round(Int, level), d)] = cross(s.t, model_T(s, d), level)
        tm, Tm = measured(run, d)
        row[@sprintf("t%d_%dmm_meas", round(Int, level), d)] =
            isempty(tm) ? NaN : cross(tm, Tm, level)
    end
    for d in (9, 19)
        row["plateau_$(d)mm"] = plateau_minutes(s.t, model_T(s, d))
        tm, Tm = measured(run, d)
        row["plateau_$(d)mm_meas"] = isempty(tm) ? NaN : plateau_minutes(tm, Tm)
    end
    # the 9 mm reading when the measured 1 mm record passes each mark
    t1, T1 = measured(run, 1)
    t9, T9 = measured(run, 9)
    for level in MARKS
        L = round(Int, level)
        tc = isempty(t1) ? NaN : 60 * cross(t1, T1, level)
        row["mark$(L)_min"] = tc / 60
        row["mark$(L)_9mm"] = isnan(tc) ? NaN : interp(s.t, s.T9, tc)
        row["mark$(L)_9mm_meas"] = isnan(tc) || isempty(t9) ? NaN : interp(t9, T9, tc)
    end
    pc = s.phase_change
    if pc === nothing
        for k in ("pv_peq_max", "X_max", "evaporated", "condensed", "water0")
            row[k] = NaN
        end
        row["accepted"] = NaN
    else
        row["pv_peq_max"] = maximum(s.pv_over_peq)
        row["X_max"] = maximum(s.X_max)
        row["evaporated"] = trapz(s.t, s.evaporation)
        row["condensed"] = trapz(s.t, s.condensation)
        row["water0"] = s.water[1]
        row["accepted"] = pc.accepted
    end
    row["condensed_over_water0"] = row["condensed"] / row["water0"]
    row["recondensed"] = row["condensed"] / row["evaporated"]
    return row
end

function variants_of(families)
    isdir(RESULTS) || return String[]
    v = filter(d -> isdir(joinpath(RESULTS, d)), readdir(RESULTS))
    return filter(d -> any(startswith(d, f * "_") for f in families), v)
end

function main(families)
    runs = ["dry_E", "y09_E", "y20_E", "y20_E_repeat", "y09_B"]
    rows = Dict{String,Any}[]
    for v in variants_of(families), run in runs
        m = metrics(v, run)
        m === nothing || push!(rows, m)
    end
    keys_ = [
        "variant",
        "run",
        "window_min",
        "rms_1mm",
        "rms_9mm",
        "rms_19mm",
        "rms_XT",
        "t373_9mm",
        "t373_9mm_meas",
        "t400_9mm",
        "t400_9mm_meas",
        "t500_9mm",
        "t500_9mm_meas",
        "t373_19mm",
        "t373_19mm_meas",
        "t400_19mm",
        "t400_19mm_meas",
        "t440_19mm",
        "t440_19mm_meas",
        "t500_19mm",
        "t500_19mm_meas",
        "Tmin_19mm",
        "t440_9mm",
        "t440_9mm_meas",
        "plateau_9mm",
        "plateau_9mm_meas",
        "plateau_19mm",
        "plateau_19mm_meas",
        "mark440_min",
        "mark440_9mm",
        "mark440_9mm_meas",
        "mark530_min",
        "mark530_9mm",
        "mark530_9mm_meas",
        "p_max_kPa",
        "pv_peq_max",
        "X_max",
        "water0",
        "evaporated",
        "condensed",
        "condensed_over_water0",
        "recondensed",
        "accepted",
        "mass_closure",
        "energy_closure",
    ]
    mkpath(RESULTS)
    open(joinpath(RESULTS, "metrics.csv"), "w") do io
        println(io, join(keys_, ","))
        for r in rows
            println(
                io,
                join(
                    [
                        k in ("variant", "run") ? r[k] : @sprintf("%.6g", r[k]) for
                        k in keys_
                    ],
                    ",",
                ),
            )
        end
    end
    mkpath(joinpath(RESULTS, "tables"))
    for fam in families
        fr = filter(r -> startswith(r["variant"], fam * "_"), rows)
        isempty(fr) && continue
        open(joinpath(RESULTS, "tables", "$fam.md"), "w") do io
            println(
                io,
                "| variant | run | window | rms 1 mm | rms 9 mm | rms 19 mm | rms X_T | 9 mm 400 K (meas) | 19 mm 400 K (meas) | 19 mm 500 K (meas) | min T19 |",
            )
            println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
            for r in sort(fr, by = r -> (r["run"], r["variant"]))
                @printf(
                    io,
                    "| %s | %s | %.0f min | %.1f | %.1f | %.1f | %.3f | %.1f (%.1f) | %.1f (%.1f) | %.1f (%.1f) | %.0f |\n",
                    r["variant"],
                    r["run"],
                    r["window_min"],
                    r["rms_1mm"],
                    r["rms_9mm"],
                    r["rms_19mm"],
                    r["rms_XT"],
                    r["t400_9mm"],
                    r["t400_9mm_meas"],
                    r["t400_19mm"],
                    r["t400_19mm_meas"],
                    r["t500_19mm"],
                    r["t500_19mm_meas"],
                    r["Tmin_19mm"]
                )
            end
        end
        phase_change_table(fam, fr, joinpath(RESULTS, "tables", "$(fam)_phase_change.md"))
        plateau_table(fam, fr, joinpath(RESULTS, "tables", "$(fam)_plateau.md"))
    end
    figures(families, runs)
    return rows
end

f1(x) = isnan(x) ? "–" : @sprintf("%.1f", x)
f2(x) = isnan(x) ? "–" : @sprintf("%.2f", x)

"""
The phase change's table, one row per run and variant of the family, the
measured row first: rms at 9 and 19 mm, the 9 mm reading at the two marks of
the measured 1 mm record, when 9 mm passes 373 K and 19 mm passes 373, 400
and 440 K, the largest gauge pore pressure, vapour over equilibrium pressure
and moisture content, the water condensed over the initial water and over
the water evaporated, and the accepted steps. Written to `path` and printed.
"""
function phase_change_table(fam, fr, path)
    io = IOBuffer()
    println(
        io,
        "| run | variant | rms 9 / 19 mm (K) | 9 mm at 440 / 530 K marks (K) | 9 mm passes 373 K (min) | 19 mm passes 373 / 400 / 440 K (min) | p_max (kPa gauge) | p_v/p_eq max | X_max | condensed / water0 | condensed / evaporated | accepted steps |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for run in sort(unique(r["run"] for r in fr))
        rr = sort(filter(r -> r["run"] == run, fr), by = r -> r["variant"])
        m = first(rr)
        @printf(
            io,
            "| %s | measured | — | %.0f / %.0f | %s | %s / %s / %s | — | — | — | — | — | — |\n",
            run,
            m["mark440_9mm_meas"],
            m["mark530_9mm_meas"],
            f1(m["t373_9mm_meas"]),
            f1(m["t373_19mm_meas"]),
            f1(m["t400_19mm_meas"]),
            f1(m["t440_19mm_meas"])
        )
        for r in rr
            @printf(
                io,
                "| %s | %s | %.2f / %.2f | %.2f / %.2f | %.3f | %.3f / %.3f / %.3f | %.1f | %s | %s | %s | %s | %s |\n",
                run,
                r["variant"],
                r["rms_9mm"],
                r["rms_19mm"],
                r["mark440_9mm"],
                r["mark530_9mm"],
                r["t373_9mm"],
                r["t373_19mm"],
                r["t400_19mm"],
                r["t440_19mm"],
                r["p_max_kPa"],
                isnan(r["pv_peq_max"]) ? "–" : @sprintf("%.4f", r["pv_peq_max"]),
                isnan(r["X_max"]) ? "–" : @sprintf("%.4f", r["X_max"]),
                isnan(r["condensed_over_water0"]) ? "–" :
                @sprintf("%.1f %%", 100r["condensed_over_water0"]),
                isnan(r["recondensed"]) ? "–" : @sprintf("%.1f %%", 100r["recondensed"]),
                isnan(r["accepted"]) ? "–" : @sprintf("%d", r["accepted"])
            )
        end
    end
    s = String(take!(io))
    write(path, s)
    println("\n$fam: phase change\n", s)
end

"""
Plateau and crossing times at 9 and 19 mm, one row per run and variant of the
family, the measured row first: the minutes from 365 to 400 K and the minutes
at which 400, 440 and 500 K are passed. Written to `path` and printed.
"""
function plateau_table(fam, fr, path)
    io = IOBuffer()
    println(
        io,
        "| run | variant | 9 mm 365→400 K (min) | 9 mm passes 400 / 440 / 500 K (min) | 19 mm 365→400 K (min) | 19 mm passes 400 / 440 / 500 K (min) | mass / energy closure |",
    )
    println(io, "|---|---|---|---|---|---|---|")
    for run in sort(unique(r["run"] for r in fr))
        rr = sort(filter(r -> r["run"] == run, fr), by = r -> r["variant"])
        m = first(rr)
        println(
            io,
            "| $run | measured | $(f1(m["plateau_9mm_meas"])) | $(f1(m["t400_9mm_meas"])) / $(f1(m["t440_9mm_meas"])) / $(f1(m["t500_9mm_meas"])) | $(f1(m["plateau_19mm_meas"])) | $(f1(m["t400_19mm_meas"])) / $(f1(m["t440_19mm_meas"])) / $(f1(m["t500_19mm_meas"])) | — |",
        )
        for r in rr
            @printf(
                io,
                "| %s | %s | %s | %s / %s / %s | %s | %s / %s / %s | %.1e / %.1e |\n",
                run,
                r["variant"],
                f1(r["plateau_9mm"]),
                f1(r["t400_9mm"]),
                f1(r["t440_9mm"]),
                f1(r["t500_9mm"]),
                f1(r["plateau_19mm"]),
                f1(r["t400_19mm"]),
                f1(r["t440_19mm"]),
                f1(r["t500_19mm"]),
                r["mass_closure"],
                r["energy_closure"]
            )
        end
    end
    s = String(take!(io))
    write(path, s)
    println("\n$fam: plateau and crossing times\n", s)
end

function figures(families, runs)
    figdir = joinpath(RESULTS, "figures")
    mkpath(figdir)
    colors = Dict(
        "sorption" => :seagreen,
        "arrhenius" => :firebrick,
        "gated" => :royalblue,
        "chan" => :darkorange,
    )
    for fam in families, run in runs
        fig = Figure(size = (1100, 420))
        ax1 =
            Axis(fig[1, 1], xlabel = "time [min]", ylabel = "T [K]", title = "$run — $fam")
        ax2 = Axis(
            fig[1, 2],
            xlabel = "time [min]",
            ylabel = "X_T",
            title = "total conversion",
        )
        w = window(run)
        vlines!(ax1, [w[2] / 60], color = :gray70, linestyle = :dot)
        for (d, mk) in ((1, :rect), (9, :diamond), (19, :utriangle))
            t, T = measured(run, d)
            isempty(t) || scatter!(
                ax1,
                t ./ 60,
                T,
                marker = mk,
                color = :black,
                markersize = 7,
                label = "$d mm",
            )
        end
        t, X = conversion(run)
        isempty(t) || scatter!(ax2, t ./ 60, X, color = :black, markersize = 6)
        any_line = false
        for law in LAWS_SHOWN
            p = series_path("$(fam)_$law", run)
            isfile(p) || continue
            s = read_series(p)
            any_line = true
            for (d, ls) in ((1, :dot), (9, :dash), (19, :solid))
                lines!(
                    ax1,
                    s.t ./ 60,
                    model_T(s, d),
                    color = colors[law],
                    linestyle = ls,
                    label = d == 19 ? law : nothing,
                )
            end
            lines!(ax2, s.t ./ 60, s.X_T, color = colors[law])
        end
        any_line || continue
        axislegend(ax1, position = :lt, merge = true, framevisible = false)
        save(joinpath(figdir, "$(fam)_$(run).png"), fig)
    end
end

families = isempty(ARGS) ? ["deck"] : ARGS
rows = main(families)
for r in sort(rows, by = r -> (r["run"], r["variant"]))
    @printf(
        "%-28s %-13s win %4.0f  rms 1/9/19 %5.1f %5.1f %5.1f  X_T %.3f  19mm@400 %5.1f (%5.1f)  9mm@400 %5.1f (%5.1f)\n",
        r["variant"],
        r["run"],
        r["window_min"],
        r["rms_1mm"],
        r["rms_9mm"],
        r["rms_19mm"],
        r["rms_XT"],
        r["t400_19mm"],
        r["t400_19mm_meas"],
        r["t400_9mm"],
        r["t400_9mm_meas"]
    )
end
