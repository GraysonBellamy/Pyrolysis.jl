# Score the solved variants against the data read from Fredlund's figures:
# temperatures at 20-70 mm (rms, the times 100, 150, 200 and 300 °C are
# reached, the time from 100 to 150 °C), the peak gauge pressures, and the
# moisture rise of the 1985 cube against its resistive record.
#   julia --project=scripts/verification scripts/fredlund_wet_wood/compare.jl
# Reads results/<variant>/<case>.csv; writes results/metrics.csv and
# results/tables/*.md, and prints the tables.
include(joinpath(@__DIR__, "FredlundWetWood.jl"))
using .FredlundWetWood
using Printf
using Statistics

const LEVELS = (373.15, 423.15, 473.15, 573.15)   # 100, 150, 200, 300 °C
const SCORED = Dict(
    "spruce14" => (20, 30, 40, 50, 70),
    "pine24" => (20, 30, 40, 50, 70),
    "spruce_dry" => (20, 30, 40, 50, 60),
)

meta_value(s, key) =
    (m = match(Regex("$key=(\\S+)"), s.meta); m === nothing ? NaN : parse(Float64, m[1]))
fmt(x; d = 1) = isnan(x) ? "–" : @sprintf("%.*f", d, x)
mins(x) = x / 60

function variants()
    isdir(RESULTS) || return String[]
    return sort(
        filter(d -> isdir(joinpath(RESULTS, d)) && d ∉ ("tables",), readdir(RESULTS)),
    )
end

"""Temperature scores of a solved `case` (series `s`) at each scored depth."""
function temperature_scores(case, s)
    meas = temperatures(case)
    out = Dict{Int,Any}()
    for dd in SCORED[case]
        haskey(meas, dd) || continue
        t, T = meas[dd]
        col = s[Symbol("T$dd")]
        keep = findall(<=(s.t[end]), t)
        rms =
            isempty(keep) ? NaN :
            sqrt(mean((interp(s.t, col, t[i]) - T[i])^2 for i in keep))
        bias = isempty(keep) ? NaN : mean(interp(s.t, col, t[i]) - T[i] for i in keep)
        out[dd] = (;
            rms,
            bias,
            model = [crossing(s.t, col, L) for L in LEVELS],
            meas = [crossing(t, T, L) for L in LEVELS],
        )
    end
    return out
end

"""The model's peak gauge pressure [Pa] and its time [s] at `depth`."""
function pressure_peak(s, depth)
    col = s[Symbol("P$depth")]
    k = argmax(col)
    return col[k], s.t[k]
end

function main()
    vs = variants()
    rows = String[]
    mkpath(joinpath(RESULTS, "tables"))
    io = IOBuffer()
    for case in ("spruce_dry", "spruce14", "pine24")
        println(io, "\n## $case: temperatures\n")
        println(
            io,
            "| variant | rms K (20 / 30 / 40 / 50 / $(last(SCORED[case])) mm) | bias K | " *
            "100 °C at (min) | 150 °C at | 200 °C at | 100→150 °C (min) | closures m / e | steps |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|")
        first_ = true
        for v in vs
            path = series_path(v, case)
            isfile(path) || continue
            s = read_series(path)
            sc = temperature_scores(case, s)
            ds = [dd for dd in SCORED[case] if haskey(sc, dd)]
            if first_
                m(j) = join((fmt(mins(sc[dd].meas[j])) for dd in ds), " / ")
                println(
                    io,
                    "| measured | – | – | $(m(1)) | $(m(2)) | $(m(3)) | " *
                    join((fmt(mins(sc[dd].meas[2] - sc[dd].meas[1])) for dd in ds), " / ") *
                    " | – | – |",
                )
                first_ = false
            end
            c(j) = join((fmt(mins(sc[dd].model[j])) for dd in ds), " / ")
            pc = s.phase_change
            println(
                io,
                "| $v | " *
                join((fmt(sc[dd].rms) for dd in ds), " / ") *
                " | " *
                join((fmt(sc[dd].bias) for dd in ds), " / ") *
                " | $(c(1)) | $(c(2)) | $(c(3)) | " *
                join((fmt(mins(sc[dd].model[2] - sc[dd].model[1])) for dd in ds), " / ") *
                " | " *
                @sprintf(
                    "%.1e / %.1e",
                    meta_value(s, "mass_closure"),
                    meta_value(s, "energy_closure")
                ) *
                " | $(pc === nothing ? "–" : pc.accepted) |",
            )
            for dd in ds
                push!(
                    rows,
                    join(
                        (
                            v,
                            case,
                            dd,
                            @sprintf("%.3f", sc[dd].rms),
                            @sprintf("%.3f", sc[dd].bias),
                            (@sprintf("%.2f", mins(x)) for x in sc[dd].model)...,
                            (@sprintf("%.2f", mins(x)) for x in sc[dd].meas)...,
                        ),
                        ",",
                    ),
                )
            end
        end
        case == "spruce_dry" && continue
        println(io, "\n## $case: peak gauge pressures (kPa) and their times (min)\n")
        P = pressures(case)
        println(
            io,
            "| variant | " *
            join(("$(p.depth) mm ($(p.status))" for p in P), " | ") *
            " | p_max | X_max |",
        )
        println(io, "|" * repeat("---|", length(P) + 3))
        println(
            io,
            "| measured | " *
            join(
                (
                    @sprintf(
                        "%.0f at %s",
                        p.p / 1e3,
                        isnan(p.t_from) ? "–" :
                        p.t_from == p.t_to ? fmt(mins(p.t_from)) :
                        fmt(mins(p.t_from)) * "–" * fmt(mins(p.t_to))
                    ) for p in P
                ),
                " | ",
            ) *
            " | – | $(fmt(d_moisture(case), d = 3)) |",
        )
        for v in vs
            path = series_path(v, case)
            isfile(path) || continue
            s = read_series(path)
            cells = String[]
            for p in P
                pk, tk = pressure_peak(s, p.depth)
                push!(cells, @sprintf("%.0f at %s", pk / 1e3, fmt(mins(tk))))
            end
            println(
                io,
                "| $v | " *
                join(cells, " | ") *
                " | $(fmt(maximum(s.p_max) / 1e3, d = 0)) | $(fmt(maximum(s.X_max), d = 3)) |",
            )
        end
    end
    println(
        io,
        "\n## pine14_1985: the moisture rise (model at the sensor depths; measured F85 Fig. 5.4)\n",
    )
    M = moisture_1985()
    X0 = deck()["cases"]["pine14_1985"]["moisture"]
    println(
        io,
        "| variant | X_peak/X0 (10 … 70 mm) | lead of the peak before local 100 °C (min) | T at the peak (°C) | width: X/X0 at 50, 60, 70 mm when the 40 mm peak passes |",
    )
    println(io, "|---|---|---|---|---|")
    println(
        io,
        "| measured | " *
        join((fmt(m.X_peak / X0, d = 2) for m in M), " / ") *
        " | " *
        join((fmt(mins(m.t_100C - m.t_peak)) for m in M), " / ") *
        " | " *
        join((fmt(m.T_at_peak - 273.15, d = 0) for m in M), " / ") *
        " | 1.14 / 1.09 / 1.07 |",
    )
    for v in vs
        path = series_path(v, "pine14_1985")
        isfile(path) || continue
        s = read_series(path)
        ratio, lead, Tp = String[], String[], String[]
        for m in M
            X = s[Symbol("X$(m.depth)")]
            T = s[Symbol("T$(m.depth)")]
            k = argmax(X)
            push!(ratio, fmt(X[k] / X0, d = 2))
            push!(lead, fmt(mins(crossing(s.t, T, 373.15) - s.t[k])))
            push!(Tp, fmt(T[k] - 273.15, d = 0))
        end
        k40 = argmax(s[:X40])
        width =
            join((fmt(s[Symbol("X$dd")][k40] / X0, d = 2) for dd in (50, 60, 70)), " / ")
        println(
            io,
            "| $v | " *
            join(ratio, " / ") *
            " | " *
            join(lead, " / ") *
            " | " *
            join(Tp, " / ") *
            " | $width |",
        )
    end
    text = String(take!(io))
    print(text)
    write(joinpath(RESULTS, "tables", "fredlund.md"), text)
    open(joinpath(RESULTS, "metrics.csv"), "w") do f
        println(
            f,
            "variant,case,depth_mm,rms_K,bias_K,t100_min,t150_min,t200_min,t300_min,t100_meas,t150_meas,t200_meas,t300_meas",
        )
        foreach(r -> println(f, r), rows)
    end
end

d_moisture(case) = deck()["cases"][case]["moisture"]

main()
