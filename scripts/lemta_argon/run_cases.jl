# Run variants (variants.jl) on the tests they apply to, or on the named tests
# or species, writing results/<variant>/<test>.csv. Mass tests are solved to
# the char-depth tests' 12 min (CHAR_EXPOSURE). A test whose solve throws is
# written as results/<variant>/<test>.failed with the error; one that returns
# without success keeps its series with the return code. Tests already written
# are skipped unless --force; --part=i/n runs every n-th test from the i-th.
#
#   julia --project=scripts/verification scripts/lemta_argon/run_cases.jl VARIANT... [TEST|SPECIES ...] [--force] [--part=i/n]
include(joinpath(@__DIR__, "LemtaArgon.jl"))
using .LemtaArgon
include(joinpath(@__DIR__, "variants.jl"))
using Printf

args = filter(a -> !startswith(a, "--"), ARGS)
force = "--force" in ARGS
part = let p = findfirst(a -> startswith(a, "--part="), ARGS)
    p === nothing ? (1, 1) : Tuple(parse.(Int, split(ARGS[p][8:end], '/')))
end
variants = filter(a -> haskey(VARIANTS, a), args)
names = filter(a -> !haskey(VARIANTS, a), args)
isempty(variants) && error("give a variant: $(sort(collect(keys(VARIANTS))))")
check_data()
tests = [mass_runs(); temperature_runs()]
unknown = setdiff(names, [[r.name for r in tests]; collect(values(SPECIES))])
isempty(unknown) || error("unknown variants, tests or species: $unknown")

for variant in variants
    v = VARIANTS[variant]
    runs = filter(v.runs, tests)
    isempty(names) || filter!(r -> r.name in names || r.species in names, runs)
    runs = runs[part[1]:part[2]:end]
    for r in runs
        path = series_path(variant, r.name)
        failed = replace(path, r"\.csv$" => ".failed")
        if (isfile(path) || isfile(failed)) && !force
            println("$variant $(r.name): exists")
            continue
        end
        rm(failed; force = true)
        t_end = r isa MassRun ? max(r.t[end], CHAR_EXPOSURE) : LemtaArgon._t_end(r)
        wall = time()
        s = try
            simulate(r; t_end, v.kwargs(r)...)
        catch err
            mkpath(dirname(failed))
            open(failed, "w") do io
                showerror(io, err, catch_backtrace())
            end
            @printf("%s %s: FAILED (%s)\n", variant, r.name, sprint(showerror, err))
            flush(stdout)
            nothing
        end
        s === nothing && continue
        write_series(path, s)
        ok = s.retcode == "Success" && all(==("Success"), s.window_retcodes)
        @printf(
            "%s %s: %s (window passes %s) in %.0f s wall; m/m₀ %.3f at %.0f s; ledgers %.1e / %.1e\n",
            variant,
            r.name,
            ok ? "Success" : "NOT SUCCESS: " * s.retcode,
            join(s.window_retcodes, ","),
            time() - wall,
            s.m[end],
            s.t[end],
            s.audit.mass.relative,
            s.audit.energy.relative
        )
        flush(stdout)
    end
end
