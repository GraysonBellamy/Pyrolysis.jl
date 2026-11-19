# Score the solved variants against the digitised data: temperature rms by
# depth and the total conversion within each run's window, the plateau and
# crossing times at 9 and 19 mm, and figures.
#   julia --project=scripts/verification scripts/bilbao_wet_wood/compare.jl [FAMILY ...]
# Reads results/<family>_<law>[_tag]/<run>.csv; writes results/metrics.csv,
# results/tables/*.md and results/figures/*.png.
include(joinpath(@__DIR__, "BilbaoWetWood.jl"))
using .BilbaoWetWood
using CairoMakie
using Printf
using Statistics

const LAWS_SHOWN = ("arrhenius", "gated", "chan")
const LEVEL = 626.0   # K: the dry run's 19 mm reading at 24 min, the end of the calibration window

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
        "meta" => s.meta,
    )
    for d in (9, 19), level in (373.0, 400.0, 500.0)
        row[@sprintf("t%d_%dmm", round(Int, level), d)] = cross(s.t, model_T(s, d), level)
        tm, Tm = measured(run, d)
        row[@sprintf("t%d_%dmm_meas", round(Int, level), d)] =
            isempty(tm) ? NaN : cross(tm, Tm, level)
    end
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
        "t500_19mm",
        "t500_19mm_meas",
        "Tmin_19mm",
    ]
    mkpath(RESULTS)
    open(joinpath(RESULTS, "metrics.csv"), "w") do io
        println(io, join(keys_, ","))
        for r in rows
            println(
                io,
                join(
                    [
                        k in ("variant", "run") ? r[k] : @sprintf("%.4g", r[k]) for
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
    end
    figures(families, runs)
    return rows
end

function figures(families, runs)
    figdir = joinpath(RESULTS, "figures")
    mkpath(figdir)
    colors = Dict("arrhenius" => :firebrick, "gated" => :royalblue, "chan" => :darkorange)
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
