# Score the model's moisture profiles against the slice profiles of Diarra et
# al. 2026, Fig. 11(a) (data/diarra_2026_fig11.csv): spruce at about 10 %,
# grain across the flux, 5 mm slices after 30, 60 and 120 min. The model run
# is the record spruce_11,42%_343 (11.4 % initially); both are compared as
# X/X0, the measured X0 the unexposed block's level (10.5 %) and the model's
# its initial moisture content, each model slice the mean of the cells whose
# centres lie within 2.5 mm of the measured slice centre. Prints, per variant,
# the slices at each time, the rms difference over the 22 points and the
# wettest slice; writes results/tables/profiles_diarra.md.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/profiles.jl [VARIANT ...]
#
# Reads results/<variant>/spruce_11,42%_343_profiles.csv (run_cases.jl writes
# it); without arguments, every variant that has it.
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using Printf
using Statistics

const RUN = "spruce_11,42%_343"
const DEEP_LEVEL = 10.5   # %, the measured profiles' unexposed level
const HALF_SLICE = 2.5    # mm

"""The digitised points: (time [h], depth [mm], moisture content [%])."""
function measured_points()
    lines = filter(
        l -> !startswith(l, "#"),
        readlines(joinpath(@__DIR__, "data", "diarra_2026_fig11.csv")),
    )
    return [Tuple(parse.(Float64, split(l, ','))) for l in lines[2:end]]
end

"""Mean moisture content of the model cells within half a slice of `depth`."""
slice_mean(profile, depth) =
    mean(X for (z, X) in profile if depth - HALF_SLICE <= z < depth + HALF_SLICE)

check_data()
r = load_run(RUN)
X0 = initial_water_g(r) / (1000 * r.rho * LemtaMoistWood.VOLUME)
points = measured_points()
times = sort(unique(p[1] for p in points))
variants =
    isempty(ARGS) ?
    sort([
        v for
        v in keys(VARIANTS) if VARIANTS[v].runs == :moist && isfile(profiles_path(v, RUN))
    ]) : ARGS

io = IOBuffer()
@printf(io, "Model X0 = %.5f (%s); measured X0 = %.1f %%.\n\n", X0, RUN, DEEP_LEVEL)
println(io, "| variant | rms difference in X/X0, 22 points | wettest slice, X/X0 |")
println(io, "|---|---|---|")
@printf(io, "| measured | — | %.2f |\n", maximum(p[3] for p in points) / DEEP_LEVEL)
details = IOBuffer()
println(details, "measured X/X0:")
for t in times
    println(
        details,
        @sprintf("  %.1f h: ", t) * join(
            (
                @sprintf("%4.1f mm %4.2f", p[2], p[3] / DEEP_LEVEL) for
                p in points if p[1] == t
            ),
            "  ",
        ),
    )
end
for v in variants
    path = profiles_path(v, RUN)
    isfile(path) || error("no profiles for $v: $path (run run_cases.jl $v)")
    prof = read_profiles(path)
    err = Float64[]
    peak = 0.0
    println(details, v)
    for t in times
        row = String[]
        for p in points
            p[1] == t || continue
            m = slice_mean(prof[t], p[2]) / X0
            push!(err, m - p[3] / DEEP_LEVEL)
            peak = max(peak, m)
            push!(row, @sprintf("%4.1f mm %4.2f", p[2], m))
        end
        println(details, @sprintf("  %.1f h: ", t) * join(row, "  "))
    end
    @printf(io, "| %s | %.3f | %.2f |\n", v, sqrt(mean(err .^ 2)), peak)
end
summary = String(take!(io))
mkpath(joinpath(RESULTS, "tables"))
write(joinpath(RESULTS, "tables", "profiles_diarra.md"), summary)
print(String(take!(details)))
println()
print(summary)
