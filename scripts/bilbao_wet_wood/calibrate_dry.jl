# Fit a family's flux factor (or, with the 1 mm record as the boundary, the
# solids' conductivity factor) to the dry run's 9 and 19 mm records over 0–24 min, by golden-section
# search on the rms. Writes calibration.toml [family].
#   julia --project=scripts/verification scripts/bilbao_wet_wood/calibrate_dry.jl FAMILY
include(joinpath(@__DIR__, "BilbaoWetWood.jl"))
using .BilbaoWetWood
using Printf
using TOML
using Statistics
include(joinpath(@__DIR__, "variants.jl"))

const WINDOW = (0.0, 24 * 60.0)

function score(family; flux_factor, k_scale = 1.0)
    kw = FAMILIES[family]
    s = simulate(RUNS["dry_E"]; law = :gated, flux_factor, k_scale, kw...)
    rms = rms_by_depth("dry_E", s; window = WINDOW)
    v = sqrt((rms[9]^2 + rms[19]^2) / 2)
    @printf(
        "  %-9s f=%.4f k=%.3f: rms 9 mm %.1f, 19 mm %.1f (1 mm %.1f) -> %.2f K\n",
        family,
        flux_factor,
        k_scale,
        rms[9],
        rms[19],
        rms[1],
        v
    )
    return v, rms
end

function golden(g, a, b; tol = 0.01)
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    gc, gd = g(c), g(d)
    while b - a > tol
        if gc < gd
            b, d, gd = d, c, gc
            c = b - φ * (b - a)
            gc = g(c)
        else
            a, c, gc = c, d, gd
            d = a + φ * (b - a)
            gd = g(d)
        end
    end
    return gc < gd ? (c, gc) : (d, gd)
end

function fit(family)
    if FAMILIES[family].front === :tc1
        return (
            1.0,
            golden(
                x -> score(family; flux_factor = 1.0, k_scale = x)[1],
                0.8,
                3.0;
                tol = 0.01,
            )[1],
        )
    end
    return (golden(x -> score(family; flux_factor = x)[1], 1.4, 3.0)[1], 1.0)
end

family = ARGS[1]
best = fit(family)
v, rms = score(family; flux_factor = best[1], k_scale = best[2])
cal = isfile(CALIBRATION) ? TOML.parsefile(CALIBRATION) : Dict{String,Any}()
cal[family] = Dict(
    "flux_factor" => round(best[1], digits = 4),
    "k_scale" => round(best[2], digits = 3),
    "rms_9mm_K" => round(rms[9], digits = 1),
    "rms_19mm_K" => round(rms[19], digits = 1),
    "rms_1mm_K" => round(rms[1], digits = 1),
)
# re-read just before writing: other families may have finished meanwhile
cal = merge(
    isfile(CALIBRATION) ? TOML.parsefile(CALIBRATION) : Dict{String,Any}(),
    Dict(family => cal[family]),
)
open(CALIBRATION, "w") do io
    println(io, "# Flux factor (and wood k factor) per thermal family, fitted to the dry")
    println(io, "# run's 9 and 19 mm records over 0-24 min by calibrate_dry.jl.")
    TOML.print(io, cal; sorted = true)
end
@printf("%s: flux_factor %.4f, k_scale %.3f, rms %.2f K\n", family, best[1], best[2], v)
