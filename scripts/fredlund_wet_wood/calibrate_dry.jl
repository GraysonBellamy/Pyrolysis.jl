# Fit factors on the solids' conductivity to the oven-dry spruce of F88
# Fig. 9.8 (its 20-60 mm readings, its 10 mm record as the boundary): `k` on
# the wood and the char, `kc` on the char alone.
#   julia --project=scripts/verification scripts/fredlund_wet_wood/calibrate_dry.jl K:KC [K:KC ...]
#       solves at each pair and appends its rms to results/calibration_scan.csv
#       (pairs can run as separate processes side by side);
#   julia --project=scripts/verification scripts/fredlund_wet_wood/calibrate_dry.jl pick
#       writes the pair of least rms in the scan to calibration.toml.
include(joinpath(@__DIR__, "FredlundWetWood.jl"))
using .FredlundWetWood
using Printf
using Statistics

const SCAN = joinpath(RESULTS, "calibration_scan.csv")

function rms_dry(k, kc)
    s = simulate(DRY; k_scale = k, kc_scale = kc)
    meas = temperatures(DRY)
    err = Float64[]
    for (dd, (t, T)) in meas
        col = s.cols[Symbol("T$dd")]
        append!(
            err,
            [interp(s.t, col, t[i]) - T[i] for i in eachindex(t) if t[i] <= s.t[end]],
        )
    end
    return sqrt(mean(abs2, err)), length(err), s.retcode
end

if get(ARGS, 1, "") == "pick"
    rows = [split(l, ',') for l in readlines(SCAN)[2:end]]
    best = rows[argmin([parse(Float64, r[3]) for r in rows])]
    open(joinpath(@__DIR__, "calibration.toml"), "w") do io
        println(
            io,
            "# Written by calibrate_dry.jl: the factors on the solids' conductivity of the",
        )
        println(
            io,
            "# least rms over the oven-dry spruce (F88 Fig. 9.8, 20-60 mm, the 10 mm record",
        )
        println(io, "# as the boundary) among the pairs in results/calibration_scan.csv.")
        println(io, "k_scale = $(best[1])\nkc_scale = $(best[2])\nrms_K = $(best[3])")
    end
    println(
        "k_scale = $(best[1]), kc_scale = $(best[2]) (rms $(best[3]) K) -> calibration.toml",
    )
else
    mkpath(RESULTS)
    for a in ARGS
        k, kc = parse.(Float64, split(a, ':'))
        r, n, rc = rms_dry(k, kc)
        @printf("k %.2f kc %.2f: rms %.1f K (%d points) %s\n", k, kc, r, n, rc)
        isfile(SCAN) || write(SCAN, "k_scale,kc_scale,rms_K,points,retcode\n")
        open(io -> @printf(io, "%.3f,%.3f,%.3f,%d,%s\n", k, kc, r, n, rc), SCAN, "a")
    end
end
