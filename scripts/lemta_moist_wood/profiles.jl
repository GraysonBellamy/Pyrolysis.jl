# Score the model's moisture profiles against the slice profiles behind Diarra
# et al. 2026, Fig. 11: spruce blocks at about 10 % and 30 % moisture, the
# grain across or along the heat flux, cut into slices after 30, 60 and 120 min
# under the cone (the authors' values, csv/Spruce_MC_Profile_<10|30>__<Perp|Para>
# <time>.csv in the along-grain data). The sliced blocks are not the recorded
# ones; each set is compared with the record nearest to it, or with a run built
# from the deck where there is none:
#
#   10 % across   spruce_11,42%_343        (11.4 %)
#   30 % across   spruce_29,9%_343         (29.9 %)
#   10 % along    spruce_9%_para           (9 %)
#   30 % along    spruce_30%_para_slices   (30 %, deck only, 20 °C, slice_runs)
#
# Both sides are compared as X/X0: the measured X0 is the unexposed block's
# mean for the 30 % sets (8 slices, 29.1–30.7 %) and 10.5 % for the 10 % sets
# (the level the across-grain profiles end at; no unexposed block was sliced);
# the model's is its initial moisture content. A model slice is the mean of
# the cells whose centres lie within 2.5 mm of the measured slice centre (the
# paper's ≈ 5 mm slices). Prints, per set and variant, the slices at each time,
# the rms difference over the set's points and the wettest slice; writes
# results/tables/profiles_diarra.md.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/profiles.jl [VARIANT ...]
#
# Reads results/<variant>/<run>_profiles.csv (run_cases.jl writes them; the
# along-grain runs with --along); without arguments, every variant that has
# them all.
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using Printf
using Statistics

const HALF_SLICE = 2.5    # mm
const TEN_PERCENT_LEVEL = 10.5   # %, the 10 % sets' unexposed level
const TIMES = ("30min" => 0.5, "1h" => 1.0, "2h" => 2.0)
const SETS = (
    (label = "10 % across", mc = 10, sheet = "Perp", run = "spruce_11,42%_343"),
    (label = "30 % across", mc = 30, sheet = "Perp", run = "spruce_29,9%_343"),
    (label = "10 % along", mc = 10, sheet = "Para", run = "spruce_9%_para"),
    (label = "30 % along", mc = 30, sheet = "Para", run = "spruce_30%_para_slices"),
)

profile_csv(mc, sheet) =
    joinpath(lemta_along_data(), "csv", "Spruce_MC_Profile_$(mc)__$sheet.csv")

"""Rows (depth [mm], moisture content [%]) of one profile sheet."""
read_sheet(path) =
    [Tuple(parse.(Float64, split(l, ','))) for l in readlines(path)[2:end] if !isempty(l)]

"""The set's points (time [h], depth [mm], moisture content [%]) and its X0 [%]."""
function measured(set)
    points = [
        (h, z, x) for (tag, h) in TIMES for
        (z, x) in read_sheet(profile_csv(set.mc, "$(set.sheet) $tag"))
    ]
    X0 =
        set.mc == 30 ? mean(last.(read_sheet(profile_csv(30, "Unexposed")))) :
        TEN_PERCENT_LEVEL
    return points, X0
end

"""Mean moisture content of the model cells within half a slice of `depth`."""
slice_mean(profile, depth) =
    mean(X for (z, X) in profile if depth - HALF_SLICE <= z < depth + HALF_SLICE)

check_data(; grain = :across)
check_data(; grain = :along)
model_run(name) =
    endswith(name, "_slices") ? only(filter(r -> r.name == name, slice_runs())) :
    load_run(name)
variants =
    isempty(ARGS) ?
    sort([
        v for v in keys(VARIANTS) if
        VARIANTS[v].runs == :moist && all(s -> isfile(profiles_path(v, s.run)), SETS)
    ]) : ARGS

io = IOBuffer()
details = IOBuffer()
for set in SETS
    points, X0m = measured(set)
    r = model_run(set.run)
    X0 = initial_water_g(r) / (1000 * r.rho * LemtaMoistWood.VOLUME)
    @printf(
        io,
        "\n%s: %s, model X0 = %.4f; measured X0 = %.2f %%, %d points.\n\n",
        set.label,
        set.run,
        X0,
        X0m,
        length(points)
    )
    println(io, "| variant | rms difference in X/X0 | wettest slice, X/X0 |")
    println(io, "|---|---|---|")
    @printf(io, "| measured | — | %.2f |\n", maximum(p[3] for p in points) / X0m)
    println(details, "$(set.label), measured X/X0:")
    for (_, t) in TIMES
        println(
            details,
            @sprintf("  %.1f h: ", t) * join(
                (@sprintf("%4.1f mm %4.2f", p[2], p[3] / X0m) for p in points if p[1] == t),
                "  ",
            ),
        )
    end
    for v in variants
        path = profiles_path(v, set.run)
        isfile(path) || error("no profiles for $v: $path (run run_cases.jl $v)")
        prof = read_profiles(path)
        err = Float64[]
        peak = 0.0
        println(details, "  $v")
        for (_, t) in TIMES
            row = String[]
            for p in points
                p[1] == t || continue
                m = slice_mean(prof[t], p[2]) / X0
                push!(err, m - p[3] / X0m)
                peak = max(peak, m)
                push!(row, @sprintf("%4.1f mm %4.2f", p[2], m))
            end
            println(details, @sprintf("    %.1f h: ", t) * join(row, "  "))
        end
        @printf(io, "| %s | %.3f | %.2f |\n", v, sqrt(mean(err .^ 2)), peak)
    end
end
summary = String(take!(io))
mkpath(joinpath(RESULTS, "tables"))
write(joinpath(RESULTS, "tables", "profiles_diarra.md"), summary)
print(String(take!(details)))
print(summary)
