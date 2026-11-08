# Cost and accuracy of Pyrolysis.jl on the full MaCFP cases.
#
#   taskset -c 1 julia --project=scripts/macfp_pmma --threads=1 scripts/macfp_pmma/time_pyrolysis.jl OUT.csv [SET ...]
#
# Run it under scripts/macfp_pmma/run_timing.py, which waits for a quiet
# machine and pins the process (scripts/timing/README.md).
#
# Every case runs a ladder of meshes and tolerances. Each rung is scored against
# the case's own converged reference (80 cells/mm, reltol 1e-8):
#   peak MLR and t50 (relative), the largest |ΔT| of the surface and back
#   histories up to the reference's burnout, and the largest mass difference
#   relative to the initial mass.
# Times: one discarded warm-up solve per process, then 5 / 3 / 1 repeats for
# solves under 1 s / 10 s / longer, each after GC; the median is reported.
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, Pyrolysis, Printf, Statistics

const LADDER = [(cells, factor) for cells in (10, 20, 40) for factor in (1, 10)]
const REFERENCE = (cells = 80, reltol = 1e-8, abstol = 1e-10)
const DEFAULT = (reltol = 1e-6, abstol = 1e-8, pressure_tolerance = 1.0)

repeats_for(seconds) = seconds < 1 ? 5 : seconds < 10 ? 3 : 1

function timed(f)
    GC.gc()
    first = @timed f()
    times = [first.time]
    for _ = 2:repeats_for(first.time)
        GC.gc()
        push!(times, (@timed f()).time)
    end
    return first.value, times
end

# The cases: every MaCFP-2 target, and two apparatus tests (q50 mass R4, q25 back T R2).
function cases(deck)
    out = []
    if haskey(deck, "predictions")
        for (q, L, t_end) in GASIFICATION_CASES
            push!(
                out,
                (
                    name = "gasification_$(q)kW_$(L)mm",
                    build = cells -> (
                        gasification_problem(deck, q, L; t_end, cells_per_mm = cells),
                        nothing,
                    ),
                ),
            )
        end
    end
    for t in APPARATUS_TESTS[[2, 10]]
        push!(
            out,
            (
                name = "apparatus_$(t.test)",
                build = cells -> apparatus_problem(deck, t; cells_per_mm = cells),
            ),
        )
    end
    return out
end

series(sol, iif) = iif === nothing ? gasification_series(sol) : apparatus_series(sol, iif)

burnout(s) = (b = crossing(s, 0.005); isnan(b) ? s.time_s[end] : b)
t50(s) = crossing(s, 0.5)

peak(s) = maximum(filter(isfinite, s.mlr_g_m2s))

# Histories are compared after the initial instant and up to the reference's
# 90 % mass loss. At t = 0 a step flux puts the surface value a half cell of
# conduction above the initial temperature (a transient of milliseconds), and
# the last seconds of a burn ramp the temperatures so steeply that a 0.1 s
# offset reads as kelvins.
function errors(s, ref)
    tb = (b = crossing(ref, 0.1); isnan(b) ? burnout(ref) : b)
    grid = filter(t -> 0 < t <= tb, ref.time_s)
    function maxdiff(col)
        d = filter(
            isfinite,
            [
                abs(
                    interp(s.time_s, getfield(s, col), x) -
                    interp(ref.time_s, getfield(ref, col), x),
                ) for x in grid
            ],
        )
        return isempty(d) ? NaN : maximum(d)
    end
    return (
        peak_mlr = abs(peak(s) / peak(ref) - 1),
        t50 = abs(t50(s) / t50(ref) - 1),
        T_top = maxdiff(:T_top_K),
        T_back = maxdiff(:T_back_K),
        mass = maxdiff(:mass_g) / ref.mass_g[1],
    )
end

# `--rescore`: recompute the error columns of an `--errors-only` CSV from the
# series it saved, without solving.
function rescore(out)
    lines = filter(!isempty, readlines(out))
    head = split(lines[1], ",")
    col(name) = findfirst(==(name), head)
    open(out, "w") do io
        println(io, lines[1])
        for l in lines[2:end]
            f = String.(split(l, ","))
            if f[col("retcode")] != "ERROR"
                dir = joinpath(dirname(out), "ladder", f[1])
                e = errors(
                    read_series(joinpath(dir, "$(f[2])_$(f[3])_$(f[4]).csv")),
                    read_series(joinpath(dir, "$(f[2])_reference.csv")),
                )
                f[col("err_peak_mlr"):col("err_mass_rel")] .=
                    [@sprintf("%.4g", x) for x in e]
            end
            println(io, join(f, ","))
        end
    end
end

# `--errors-only`: one solve per rung, no timing; every series saved beside OUT.
function main(args)
    "--rescore" in args && return rescore(only(filter(a -> !startswith(a, "--"), args)))
    errors_only = "--errors-only" in args
    args = filter(!=("--errors-only"), args)
    out = args[1]
    wanted = args[2:end]
    decks = load_decks(; gasification_only = true)
    isempty(wanted) || filter!(d -> d["id"] in wanted, decks)
    mkpath(dirname(out))
    io = open(out, "w")
    println(
        io,
        "set,case,cells_per_mm,tol_factor,median_s,min_s,max_s,repeats,retcode,steps,err_peak_mlr,err_t50,err_T_top_K,err_T_back_K,err_mass_rel",
    )
    warm = false
    for deck in decks, c in cases(deck)
        id = deck["id"]
        rm = MacfpPMMA.radiation_model(deck)
        p, iif = c.build(REFERENCE.cells)
        ref_sol = solve(
            p;
            SOLVE_KWARGS...,
            radiation_model = rm,
            reltol = REFERENCE.reltol,
            abstol = REFERENCE.abstol,
            pressure_tolerance = 1e-3,
        )
        ref = series(ref_sol, iif)
        if !warm      # the first solve of the process carries compilation
            pw, _ = c.build(10)
            solve(pw; SOLVE_KWARGS..., radiation_model = rm)
            warm = true
        end
        for (cells, factor) in LADDER
            p, iif = c.build(cells)
            kw = (;
                SOLVE_KWARGS...,
                radiation_model = rm,
                reltol = DEFAULT.reltol * factor,
                abstol = DEFAULT.abstol * factor,
                pressure_tolerance = DEFAULT.pressure_tolerance * factor,
            )
            sol, times = try
                errors_only ? (solve(p; kw...), [NaN]) : timed(() -> solve(p; kw...))
            catch err
                println(
                    io,
                    "$id,$(c.name),$cells,$factor,NaN,NaN,NaN,0,ERROR,0,NaN,NaN,NaN,NaN,NaN",
                )
                flush(io)
                continue
            end
            sr = series(sol, iif)
            e = errors(sr, ref)
            if errors_only
                write_series(
                    joinpath(
                        dirname(out),
                        "ladder",
                        id,
                        "$(c.name)_$(cells)_$(factor).csv",
                    ),
                    sr,
                )
                write_series(
                    joinpath(dirname(out), "ladder", id, "$(c.name)_reference.csv"),
                    ref,
                )
            end
            println(
                io,
                @sprintf(
                    "%s,%s,%d,%g,%.5g,%.5g,%.5g,%d,%s,%d,%.4g,%.4g,%.4g,%.4g,%.4g",
                    id,
                    c.name,
                    cells,
                    factor,
                    median(times),
                    minimum(times),
                    maximum(times),
                    length(times),
                    sol.retcode,
                    sol.statistics.rhs_evaluations,
                    e...
                )
            )
            flush(io)
            @printf(
                "%-12s %-28s %3d ×%-3g %8.3f s  peak %.2e t50 %.2e Ttop %.2f Tback %.2f\n",
                id,
                c.name,
                cells,
                factor,
                median(times),
                e.peak_mlr,
                e.t50,
                e.T_top,
                e.T_back
            )
        end
    end
    close(io)
end

main(ARGS)
