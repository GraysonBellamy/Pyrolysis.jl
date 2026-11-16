# What the dry runs need: per species, the factor on the absorbed flux (α q)
# and on the conductivity across the grain that minimise the pooled rms of
# the temperatures over the heating period, on a grid
# (results/calibration_grid.csv, results/calibration.toml). spruce_test2_0%
# starts with the inside warmer than its faces and is left out.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/calibrate_dry.jl
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using Printf
using Statistics

const FQ = 0.95:0.025:1.25
const FK = 0.80:0.05:1.40
const LEFT_OUT = ("spruce_test2_0%_343",)

function pooled_rms(r, s)
    keep = findall(t -> t <= min(s.t[end], r.t[end]), r.t)
    e = Float64[]
    for (j, d) in enumerate(r.depths)
        js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
        append!(e, [interp(s.t, view(s.T, :, js), t) for t in r.t[keep]] .- r.T[keep, j])
    end
    return sqrt(mean(e .^ 2))
end

runs = filter(r -> r.mc == 0 && !(r.name in LEFT_OUT), load_runs())
best = Dict{String,Tuple{Float64,Float64,Float64,Float64}}()
mkpath(RESULTS)
open(joinpath(RESULTS, "calibration_grid.csv"), "w") do io
    println(io, "species,flux_factor,k_factor,rms_K")
    for species in sort(unique(r.species for r in runs))
        rs = filter(r -> r.species == species, runs)
        base = sqrt(mean(pooled_rms(r, simulate(r))^2 for r in rs))
        b = (Inf, NaN, NaN)
        for fq in FQ, fk in FK
            e = sqrt(
                mean(
                    pooled_rms(r, simulate(r; flux = 4.0e3 * fq, k_scale = fk))^2 for
                    r in rs
                ),
            )
            @printf(io, "%s,%.3f,%.3f,%.3f\n", species, fq, fk, e)
            e < b[1] && (b = (e, fq, fk))
        end
        best[species] = (b[2], b[3], b[1], base)
        @printf(
            "%s: base rms %.2f K; best flux ×%.3f, k ×%.2f, rms %.2f K\n",
            species,
            base,
            b[2],
            b[3],
            b[1]
        )
        flush(stdout)
    end
end
open(joinpath(RESULTS, "calibration.toml"), "w") do io
    println(
        io,
        "# factor on the absorbed flux and on k across the grain fitted to the dry runs",
    )
    for (sp, (fq, fk, e, base)) in sort(collect(best))
        println(
            io,
            "[$sp]\nflux_factor = $fq\nk_factor = $fk\nrms_K = $(round(e, digits = 3))\nrms_base_K = $(round(base, digits = 3))\n",
        )
    end
end
