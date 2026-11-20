# Score the solved variants against the digitised data: the bottom-face
# temperature (rms over the record and over the moisture window, the readings
# at marks, the crossing times and the time spent between 365 and 395 K), the
# mass loss rate, the mass lost, and the phase change's diagnostics; the
# paper's ThermaKin model (the red lines of Figs 6 and 7) is scored the same
# way. Figures of the main variants.
#   julia --project=scripts/verification scripts/lee_ptw_capa/compare.jl
# Reads results/<variant>/<run>.csv; writes results/metrics.csv,
# results/tables/<run>.md and results/figures/*.png, and prints the tables.
include(joinpath(@__DIR__, "LeePtwCapa.jl"))
using .LeePtwCapa
using CairoMakie
using Printf
using Statistics

const P_ATM = 101325.0
# the moisture window [s] (the bottom face's rise to the end of its shoulder)
const WINDOW = Dict("63" => 150.0, "29" => 300.0)
const MARKS = Dict("63" => (40, 60, 80, 100, 120), "29" => (100, 130, 160, 200, 250))
const LEVELS = (340.0, 360.0, 375.0, 400.0)
const PLATEAU = (365.0, 395.0)
const MLR_EARLY = Dict("63" => 30.0, "29" => 60.0)
const SHOWN = (
    "sorption",
    "thermakin",
    "sorption_lambda1e-7",
    "sorption_lambda2e-5",
    "sorption_kappa100",
    "freewater",
)

trapz(t, y) = sum(0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))
meta_value(s, key) =
    (m = match(Regex("$key=(\\S+)"), s.meta); m === nothing ? NaN : parse(Float64, m[1]))

"""rms of `f(t)` against (t, y) at the points `keep`."""
rms_at(f, t, y, keep) = isempty(keep) ? NaN : sqrt(mean((f(t[i]) - y[i])^2 for i in keep))

"""Seconds a (t, T) series spends between the levels, from first reaching the
lower to first reaching the upper."""
plateau_seconds(t, T) =
    crossing_seconds(t, T, PLATEAU[2]) - crossing_seconds(t, T, PLATEAU[1])

"""
The scores of a bottom-face history `Tf(t)` and mass loss rate `Mf(t)` of
`run` (functions of time [s]); `series` adds what a solve has.
"""
function score(run, Tf, Mf; series = nothing)
    D = bottom_temperature(run)
    M = mass_loss_rate(run)
    tend = RUNS[run].t_end
    keep = findall(i -> D.scored[i] && D.t[i] >= 5 && !isnan(D.T[i]), eachindex(D.t))
    early = filter(i -> D.t[i] <= WINDOW[run], keep)
    keepM = findall(i -> !isnan(M.MLR[i]) && M.t[i] >= 5, eachindex(M.t))
    tg = collect(0.0:1.0:tend)
    Tg = Tf.(tg)
    row = Dict{String,Any}(
        "rms_T" => rms_at(Tf, D.t, D.T, keep),
        "rms_T_window" => rms_at(Tf, D.t, D.T, early),
        "rms_MLR" => rms_at(Mf, M.t, M.MLR, keepM),
        "MLR_early" => mean(Mf.(0.0:1.0:MLR_EARLY[run])),
        "plateau_s" => plateau_seconds(tg, Tg),
    )
    for m in MARKS[run]
        row["T_at_$m"] = Tf(m)
    end
    for L in LEVELS
        row[@sprintf("t_%d", round(Int, L))] = crossing_seconds(tg, Tg, L)
    end
    if series !== nothing
        s = series
        row["mass_lost"] = s.mass[1] - s.mass[end]
        row["mass0"] = s.mass[1]
        row["p_max_kPa"] = (maximum(s.p_max) - P_ATM) / 1000
        row["X_back_max"] = maximum(s.X_back)
        row["X_max"] = maximum(s.X_max)
        row["T_bottom_min"] = minimum(s.T_bottom)
        row["evaporated"] = trapz(s.t, s.evaporation)
        row["condensed"] = trapz(s.t, s.condensation)
        row["mass_closure"] = meta_value(s, "mass_closure")
        row["energy_closure"] = meta_value(s, "energy_closure")
        row["retcode"] = (m = match(r"retcode=(\S+)", s.meta); m === nothing ? "" : m[1])
        pc = s.phase_change
        row["pv_peq_max"] = pc === nothing ? NaN : maximum(s.pv_over_peq)
        row["accepted"] = pc === nothing ? NaN : pc.accepted
    end
    return row
end

"""The measured curve scored as a model: its own crossing times and plateau."""
function measured_row(run)
    D = bottom_temperature(run)
    M = mass_loss_rate(run)
    ok = findall(i -> D.scored[i] && !isnan(D.T[i]), eachindex(D.t))
    Tf(x) = interp(D.t[ok], D.T[ok], x)
    okM = findall(i -> !isnan(M.MLR[i]), eachindex(M.t))
    Mf(x) = interp(M.t[okM], M.MLR[okM], x)
    return score(run, Tf, Mf)
end

function paper_row(run)
    D = bottom_temperature(run)
    M = mass_loss_rate(run)
    ok = findall(i -> !isnan(D.model[i]), eachindex(D.t))
    okM = findall(i -> !isnan(M.model[i]), eachindex(M.t))
    return score(
        run,
        x -> interp(D.t[ok], D.model[ok], x),
        x -> interp(M.t[okM], M.model[okM], x),
    )
end

function variants()
    isdir(RESULTS) || return String[]
    return sort(
        filter(
            d -> isdir(joinpath(RESULTS, d)) && d ∉ ("tables", "figures"),
            readdir(RESULTS),
        ),
    )
end

fmt(x; d = 1) =
    x isa Number && isnan(x) ? "–" : x isa AbstractString ? x : @sprintf("%.*f", d, x)

function table(io, run, rows)
    marks = MARKS[run]
    println(
        io,
        "| variant | rms T (K) | rms T, ≤ $(round(Int, WINDOW[run])) s | " *
        join(("T at $m s" for m in marks), " | ") *
        " | 340 / 360 / 375 / 400 K at (s) | 365→395 K (s) | rms MLR (g/m²s) | MLR 0–$(round(Int, MLR_EARLY[run])) s | mass lost (kg/m²) | p_max (kPa g) | X back max | p_v/p_eq max | steps | closures m / e |",
    )
    println(io, "|" * repeat("---|", 13 + length(marks)))
    for (name, r) in rows
        g(k; d = 1) = fmt(get(r, k, NaN); d)
        println(
            io,
            "| $name | $(g("rms_T")) | $(g("rms_T_window")) | " *
            join((g("T_at_$m") for m in marks), " | ") *
            " | $(g("t_340", d = 0)) / $(g("t_360", d = 0)) / $(g("t_375", d = 0)) / $(g("t_400", d = 0)) | $(g("plateau_s", d = 0)) | " *
            "$(fmt(1000 * get(r, "rms_MLR", NaN), d = 2)) | $(fmt(1000 * get(r, "MLR_early", NaN), d = 2)) | $(g("mass_lost", d = 3)) | " *
            "$(g("p_max_kPa", d = 0)) | $(g("X_back_max", d = 3)) | $(g("pv_peq_max", d = 4)) | $(g("accepted", d = 0)) | " *
            "$(haskey(r, "mass_closure") ? @sprintf("%.1e / %.1e", r["mass_closure"], r["energy_closure"]) : "–") |",
        )
    end
end

function figures(run, series)
    mkpath(joinpath(RESULTS, "figures"))
    D = bottom_temperature(run)
    M = mass_loss_rate(run)
    fig = Figure(size = (1100, 450))
    ax1 = Axis(fig[1, 1], xlabel = "t (s)", ylabel = "T bottom (K)", title = "$run kW/m²")
    ax2 = Axis(fig[1, 2], xlabel = "t (s)", ylabel = "MLR (kg m⁻² s⁻¹)")
    ok = findall(!isnan, D.T)
    band!(ax1, D.t[ok], D.lo[ok], D.hi[ok], color = (:gray, 0.4), label = "measured (band)")
    lines!(ax1, D.t, D.model, color = :red, linestyle = :dash, label = "L25 ThermaKin")
    okM = findall(!isnan, M.MLR)
    band!(ax2, M.t[okM], M.lo[okM], M.hi[okM], color = (:gray, 0.4))
    lines!(ax2, M.t, M.model, color = :red, linestyle = :dash)
    for (v, s) in series
        lines!(ax1, s.t, s.T_bottom, label = v)
        lines!(ax2, s.t, s.MLR)
    end
    run == "63" ? xlims!(ax1, 0, 250) : xlims!(ax1, 0, 600)
    ylims!(ax1, 295, run == "63" ? 600 : 560)
    axislegend(ax1, position = :lt, labelsize = 10)
    save(joinpath(RESULTS, "figures", "$(run).png"), fig)
end

function main()
    vs = variants()
    out = Dict{String,Any}[]
    mkpath(joinpath(RESULTS, "tables"))
    for run in ("63", "29")
        rows = Pair{String,Dict{String,Any}}[]
        push!(rows, "measured (Fig. 6, 7)" => measured_row(run))
        push!(rows, "L25 ThermaKin (Fig. 6, 7)" => paper_row(run))
        shown = Pair{String,Any}[]
        for v in vs
            path = series_path(v, run)
            isfile(path) || continue
            s = read_series(path)
            r = score(
                run,
                x -> interp(s.t, s.T_bottom, x),
                x -> interp(s.t, s.MLR, x);
                series = s,
            )
            push!(rows, v => r)
            push!(out, merge(r, Dict("variant" => v, "run" => run)))
            v in SHOWN && push!(shown, v => s)
        end
        open(joinpath(RESULTS, "tables", "$run.md"), "w") do io
            table(io, run, rows)
        end
        println("\n## $run kW/m²\n")
        table(stdout, run, rows)
        isempty(shown) || figures(run, shown)
    end
    keys_ = sort(unique(vcat(collect.(keys.(out))...)))
    open(joinpath(RESULTS, "metrics.csv"), "w") do io
        println(io, join(keys_, ","))
        for r in out
            println(
                io,
                join(
                    (x -> x isa Number ? @sprintf("%.6g", x) : string(x)).(
                        get.(Ref(r), keys_, NaN),
                    ),
                    ",",
                ),
            )
        end
    end
end

main()
