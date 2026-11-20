# Run variants (keys of LemtaMoistWood.VARIANTS) on their runs, or on the
# named ones, writing results/<variant>/<run>.csv and, for a moist run, its
# moisture profiles <run>_profiles.csv. Runs already written are skipped
# unless --force.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/run_cases.jl VARIANT... [RUN ...] [--force]
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using Printf

args = filter(a -> !startswith(a, "--"), ARGS)
force = "--force" in ARGS
variants = filter(a -> haskey(VARIANTS, a), args)
names = filter(a -> !haskey(VARIANTS, a), args)
isempty(variants) && error("give a variant: $(sort(collect(keys(VARIANTS))))")
check_data()
all_runs = load_runs()
unknown = setdiff(names, [r.name for r in all_runs])
isempty(unknown) || error("unknown variants or runs: $unknown")

for variant in variants
    v = VARIANTS[variant]
    runs = filter(r -> v.runs == :dry ? r.mc == 0 : r.mc > 0, all_runs)
    isempty(names) || filter!(r -> r.name in names, runs)
    for r in runs
        path = series_path(variant, r.name)
        if isfile(path) && !force
            println("$variant $(r.name): exists")
            continue
        end
        s = simulate(r; v.kwargs(r)...)
        write_series(path, s)
        s.profiles === nothing || write_profiles(profiles_path(variant, r.name), s)
        @printf(
            "%s %s: %s in %.1f s, mass lost %.2f g by %.0f s\n",
            variant,
            r.name,
            s.retcode,
            s.solve_time,
            mass_loss(s, min(9000.0, s.t[end])),
            min(9000.0, s.t[end])
        )
        flush(stdout)
    end
end
