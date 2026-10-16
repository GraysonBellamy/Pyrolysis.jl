#!/usr/bin/env julia
# Deep-profiling harness for Pyrolysis.jl: runs the canonical configurations
# from test/perf/perf_helpers.jl under Profile (CPU flamegraphs) and
# Profile.Allocs (allocation sites), plus a wall-clock/GC summary, so
# slowdowns and allocations can be attributed to file:line.
#
# Usage:
#   julia --project=scripts/profile scripts/profile/profile.jl [options]
#
# Options (all --key=value):
#   --config=NAME    baseline_eulerian | ale_no_radiation |
#                    eulerian_surface_absorption | baseline_eulerian_darcy |
#                    ale_beer_lambert_darcy | all
#                    (default: ale_beer_lambert_darcy — the heaviest physics)
#   --mode=MODE      summary | cpu | allocs | all   (default: all)
#   --target=T       residual | jacobian | solve | all   (default: all)
#   --t-final=SEC    end-to-end solve horizon [s]   (default: 5.0)
#   --seconds=SEC    wall-clock budget per hot-op profile window (default: 3.0)
#   --alloc-rate=R   Profile.Allocs sample rate for the solve target
#                    (default: 0.01; residual/jacobian always use 1.0 over a
#                    record-budgeted window). Raise with care: every sampled
#                    allocation holds a full backtrace in memory, and an
#                    unbounded run has OOM-killed a 30 GiB machine.
#   --from-c=BOOL    include C/runtime frames (GC, libm) in CPU profiles
#                    (default: false — Julia-only frames attribute cleanest)
#   --strategy=S     direct | approx | structured | matrixfree — the
#                    Structured-backend linear-solve strategy (default:
#                    direct = full assembly, the verification anchor;
#                    approx = the reduced-ALE profile that drops the dense
#                    ALE couplings — the production target for ALE)
#
# Outputs:
#   - Terminal: timing/allocation/GC summary, integrator stats, top CPU frames,
#     top allocation sites and types.
#   - scripts/profile/output/<config>__<target>__{cpu,allocs}.pb.gz
#     Open interactively with:
#       julia --project=scripts/profile -e 'using PProf; PProf.refresh(file="<file>"); readline()'
#     or `pprof -http=:8080 <file>` if you have the pprof binary.
#
# Note: residual/jacobian hot ops are evaluated at the *initial* state (same
# convention as the zero-allocation gates in test/perf/). The end-to-end solve
# target exercises the realistic mixed-state path.

using Logging
using Profile
using Printf
using PProf
using Pyrolysis
using Pyrolysis: residual!
using Pyrolysis.Jacobian:
    JacobianSpec, Structured, DirectSolve, ApproxSolve, compute_jacobian!
using Pyrolysis.Materials: SURFACE_ABSORPTION, BEER_LAMBERT, NO_RADIATION

const ROOT = dirname(dirname(@__DIR__))
const OUTDIR = joinpath(@__DIR__, "output")

# Canonical (def, ws, u, du) builders + canonical materials.
include(joinpath(ROOT, "test", "perf", "perf_helpers.jl"))

const BUILDERS = Dict(
    "baseline_eulerian" => build_baseline_eulerian,
    "ale_no_radiation" => build_ale_no_radiation,
    "eulerian_surface_absorption" => build_eulerian_surface_absorption,
    "baseline_eulerian_darcy" => build_baseline_eulerian_darcy,
    "ale_beer_lambert_darcy" => build_ale_beer_lambert_darcy,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
function parse_cli(args)
    opts = Dict{String,String}(
        "config" => "ale_beer_lambert_darcy",
        "mode" => "all",
        "target" => "all",
        "t-final" => "5.0",
        "seconds" => "3.0",
        "alloc-rate" => "0.01",
        "from-c" => "false",
        "strategy" => "direct",
    )
    for a in args
        startswith(a, "--") || error("unrecognized argument $a (use --key=value)")
        kv = split(a[3:end], '='; limit = 2)
        length(kv) == 2 || error("missing value in $a (use --key=value)")
        haskey(opts, kv[1]) || error("unknown option --$(kv[1])")
        opts[kv[1]] = kv[2]
    end
    return opts
end

_expand(value, all_names) = value == "all" ? all_names : [value]

const STRATEGIES = Dict("direct" => DirectSolve, "approx" => ApproxSolve)

# ---------------------------------------------------------------------------
# End-to-end solve cases — mirror the canonical hot-op configs, but with the
# radiative flux delivered through a RadiativeFluxBC so `solve` picks it up
# (the hot-op builders inject intensity via `to_problem_def` instead).
# ---------------------------------------------------------------------------
function build_solve_case(name; t_final, strat_name = "direct")
    radiative_top() =
        RadiativeFluxBC(flux = 50e3) +
        RadiativeBC(T_inf = 300.0) +
        ConvectiveBC(h = 10.0, T_inf = 300.0)
    # Darcy needs somewhere for the pressure to vent — cone-style top mass BCs
    # (matches examples/07; the hot-op fixtures keep the impermeable default
    # because they never integrate).
    venting_mass_bcs() = (
        top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
        bottom = ImpermeableBC(),
    )
    # Solvable wood variant: the hot-op fixture (`_wood_like_material`) has no
    # ρ_intrinsic, so φ₀ = 0 and `solve` rejects it at validation for any
    # permeable material. Wood ρ_intrinsic ≈ 1500 kg/m³ gives φ₀ = 2/3. The
    # hot-op fixture itself must stay as-is — the perf gates pin it.
    function _wood_solvable()
        virgin = SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            α = 1000.0,
            κ = 1e-14,
            λ = 1e-5,
            ρ_intrinsic = 1500.0,
        )
        gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        return Material(name = :wood, components = (virgin, gas), reactions = (rxn,))
    end

    if name == "baseline_eulerian" || name == "ale_no_radiation"
        mat = _pmma_material()
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        radiation_model = NO_RADIATION
    elseif name == "eulerian_surface_absorption"
        mat = _pmma_material()
        bc = BoundaryConditionSet(thermal = (top = radiative_top(), bottom = AdiabaticBC()))
        radiation_model = SURFACE_ABSORPTION
    elseif name == "baseline_eulerian_darcy"
        mat = _wood_solvable()
        bc = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            mass = venting_mass_bcs(),
        )
        radiation_model = NO_RADIATION
    elseif name == "ale_beer_lambert_darcy"
        mat = _wood_solvable()
        bc = BoundaryConditionSet(
            thermal = (top = radiative_top(), bottom = AdiabaticBC()),
            mass = venting_mass_bcs(),
        )
        radiation_model = BEER_LAMBERT
    else
        error("unknown config $name")
    end

    use_ale = startswith(name, "ale_")
    mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, t_final),
        dt_initial = 1e-3,
        output_times = [0.0, t_final / 2, t_final],
    )
    kwargs = (progress = false, use_ale = use_ale, radiation_model = radiation_model)
    if strat_name != "direct"
        kwargs = (kwargs..., jacobian = Structured(strategy = STRATEGIES[strat_name]()))
    end
    return problem, kwargs
end

# ---------------------------------------------------------------------------
# Pretty-printing helpers
# ---------------------------------------------------------------------------
# Flush after every block: julia fully buffers redirected stdout, and a killed
# run (e.g. OOM) would otherwise lose the entire log.
banner(s) = (println("\n", "═"^78, "\n  ", s, "\n", "═"^78); flush(stdout))
section(s) = (println("\n  ── ", s, " ", "─"^max(0, 70 - length(s))); flush(stdout))

function fmt_bytes(b)
    b < 1024 && return @sprintf("%d B", b)
    b < 1024^2 && return @sprintf("%.1f KiB", b / 1024)
    b < 1024^3 && return @sprintf("%.1f MiB", b / 1024^2)
    return @sprintf("%.2f GiB", b / 1024^3)
end

fmt_time(t) =
    t < 1e-6 ? @sprintf("%.1f ns", t * 1e9) :
    t < 1e-3 ? @sprintf("%.2f μs", t * 1e6) :
    t < 1.0 ? @sprintf("%.2f ms", t * 1e3) : @sprintf("%.2f s", t)

_srcpath(file) = replace(String(file), ROOT * "/" => "")

# Fully-parametric type names can run to kilobytes — keep terminal lines sane.
_truncate(s, n = 140) = length(s) <= n ? s : first(s, n - 1) * "…"

# ---------------------------------------------------------------------------
# Summary mode — wall clock, allocation, and GC share per target
# ---------------------------------------------------------------------------

# Pick a repetition count so the timed window lasts ~`seconds` of wall clock.
function calibrate_reps(loop!; seconds)
    loop!(1)            # compile
    t1 = @elapsed loop!(1)
    return clamp(ceil(Int, seconds / max(t1, 1e-9)), 10, 5_000_000)
end

function summarize_hot_op(label, loop!, n)
    GC.gc()
    stats = @timed loop!(n)
    @printf(
        "  %-12s %12s/eval   %10s/eval   %8.1f%% GC   (%d evals, %s total)\n",
        label,
        fmt_time(stats.time / n),
        fmt_bytes(stats.bytes / n),
        100 * stats.gctime / max(stats.time, eps()),
        n,
        fmt_time(stats.time),
    )
end

function summarize_solve(config; t_final, strat_name)
    # Compilation is excluded by the warm-up solve in `run_config`.
    problem, kwargs = build_solve_case(config; t_final, strat_name)
    GC.gc()
    stats = @timed solve(problem; kwargs...)
    sol = stats.value
    @printf(
        "  %-12s %12s        %10s        %8.1f%% GC   (tspan 0–%.1f s, retcode %s)\n",
        "solve",
        fmt_time(stats.time),
        fmt_bytes(stats.bytes),
        100 * stats.gctime / max(stats.time, eps()),
        t_final,
        sol.retcode,
    )
    raw = sol.raw_solution
    if raw !== nothing && hasproperty(raw, :stats) && raw.stats !== nothing
        s = raw.stats
        @printf(
            "  %-12s nf=%d  njacs=%d  nsolve=%d  naccept=%d  nreject=%d\n",
            "  integrator",
            getfield_or(s, :nf, -1),
            getfield_or(s, :njacs, -1),
            getfield_or(s, :nsolve, -1),
            getfield_or(s, :naccept, -1),
            getfield_or(s, :nreject, -1),
        )
    end
    return nothing
end

getfield_or(x, name, default) = hasproperty(x, name) ? getproperty(x, name) : default

# Run `f`, reporting (not propagating) failures — a crashing solve case must
# not kill profiling of the remaining targets/configs.
function _guarded(f, what)
    try
        f()
        return true
    catch err
        msg = sprint(showerror, err; context = :limit => true)
        println("  ✗ $what failed: ", _truncate(first(split(msg, '\n')), 160))
        flush(stdout)
        return false
    end
end

# Allocation profiling holds a full backtrace per record and has OOM-killed a
# 30 GiB machine when left unbounded. Refuse to start without headroom.
function enough_memory(what; need_gib = 4.0)
    free = Sys.free_memory() / 2^30
    free >= need_gib && return true
    println(
        "  ⚠ skipping $what: only $(round(free, digits = 1)) GiB free ",
        "(needs ≥ $(need_gib) GiB headroom)",
    )
    flush(stdout)
    return false
end

# ---------------------------------------------------------------------------
# CPU mode — Profile.@profile → pprof file + terminal flat summary
# ---------------------------------------------------------------------------

# PProf on Julia ≥1.12 emits a known, harmless "Unexpected 0 in data" error
# for an occasional malformed sample-metadata block; the profile file is
# still written and loads fine. Filter just that message, pass the rest on.
struct DropPProfParseError <: Logging.AbstractLogger
    parent::Logging.AbstractLogger
end
Logging.min_enabled_level(l::DropPProfParseError) = Logging.min_enabled_level(l.parent)
Logging.shouldlog(l::DropPProfParseError, args...) = Logging.shouldlog(l.parent, args...)
Logging.catch_exceptions(l::DropPProfParseError) = Logging.catch_exceptions(l.parent)
function Logging.handle_message(l::DropPProfParseError, level, msg, args...; kwargs...)
    msg isa AbstractString && startswith(msg, "Unexpected 0 in data") && return nothing
    return Logging.handle_message(l.parent, level, msg, args...; kwargs...)
end

quiet_pprof(f) = Logging.with_logger(f, DropPProfParseError(Logging.current_logger()))

function profile_cpu(label, outfile, run_profiled!; delay, from_c)
    Profile.clear()
    Profile.init(n = 10_000_000, delay = delay)
    Profile.@profile run_profiled!()

    quiet_pprof() do
        pprof(; out = outfile, web = false, from_c = from_c)
    end
    println("  CPU profile written → ", _srcpath(outfile))

    # Terminal teaser: hottest frames by SELF time (the Overhead column).
    # Profile.print sorts ascending, so the hot frames are the last lines.
    # The idle profile-listener task contributes one all-overhead `poptask`
    # frame — drop it, it is harness noise, not workload.
    buf = IOBuffer()
    Profile.print(
        IOContext(buf, :displaysize => (1000, 200)),
        format = :flat,
        sortedby = :overhead,
        combine = true,
        mincount = 10,
    )
    lines = split(String(take!(buf)), '\n')
    header = lines[1:min(end, 2)]
    data = filter(lines[min(end, 3):end]) do l
        !isempty(l) && !occursin("poptask", l) && !occursin("Total snapshots", l)
    end
    keep = 25
    section("top CPU frames ($label, flat, by self time)")
    foreach(l -> println("  ", l), header)
    foreach(l -> println("  ", l), reverse(data)[1:min(end, keep)])
    length(data) > keep && println("  … (full detail in the pprof file)")
    flush(stdout)

    # Release the sample buffer eagerly — runs accumulate across configs.
    Profile.clear()
    GC.gc()
end

# ---------------------------------------------------------------------------
# Allocs mode — Profile.Allocs.@profile → pprof file + site/type aggregation
# ---------------------------------------------------------------------------
function profile_allocs(label, outfile, run_profiled!; sample_rate)
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate = sample_rate run_profiled!()
    results = Profile.Allocs.fetch()
    allocs = results.allocs

    section("allocations ($label, sample rate $sample_rate)")
    if isempty(allocs)
        println("  no allocations sampled ✓")
        return
    end

    # Stack resolution is expensive — cap how many records we attribute.
    cap = min(length(allocs), 100_000)
    cap < length(allocs) &&
        println("  ⚠ aggregating first $cap of $(length(allocs)) sampled allocations")

    bysite = Dict{String,Vector{Int}}()   # site → [count, bytes]
    bytype = Dict{String,Vector{Int}}()
    for a in view(allocs, 1:cap)
        site = "(outside Pyrolysis src)"
        for fr in a.stacktrace
            file = String(fr.file)
            if occursin(joinpath("Pyrolysis", "src"), file) ||
               startswith(file, joinpath(ROOT, "src"))
                site = string(fr.func, " @ ", _srcpath(file), ":", fr.line)
                break
            end
        end
        acc = get!(() -> [0, 0], bysite, site)
        acc[1] += 1
        acc[2] += a.size
        tacc = get!(() -> [0, 0], bytype, string(a.type))
        tacc[1] += 1
        tacc[2] += a.size
    end

    total_bytes = sum(a.size for a in view(allocs, 1:cap))
    println(
        "  $(length(allocs)) sampled allocations, ",
        fmt_bytes(total_bytes),
        " across $(length(bysite)) sites",
    )
    for (title, dict) in (("by Pyrolysis call site", bysite), ("by type", bytype))
        println("\n  top $title (bytes · count):")
        ranked = sort(collect(dict); by = kv -> -kv[2][2])
        for (key, (count, bytes)) in ranked[1:min(end, 12)]
            @printf("    %10s  %8d   %s\n", fmt_bytes(bytes), count, _truncate(key))
        end
    end

    # Hand PProf only the capped slice: every record carries a backtrace, and
    # converting millions of them is exactly what OOM-killed a prior run.
    capped = cap < length(allocs) ? Profile.Allocs.AllocResults(allocs[1:cap]) : results
    quiet_pprof() do
        PProf.Allocs.pprof(capped; out = outfile, web = false)
    end
    println("\n  allocation profile written → ", _srcpath(outfile))
    cap < length(allocs) && println("  (pprof file holds the same first $cap records)")
    flush(stdout)

    Profile.Allocs.clear()
    GC.gc()
end

# One profiled eval gives the records-per-eval rate, so the real window can be
# sized to stay under `target_records` total (bounds memory, keeps coverage).
function calibrate_alloc_evals(loop!; target_records = 200_000, max_evals = 100)
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate = 1.0 loop!(1)
    per_eval = length(Profile.Allocs.fetch().allocs)
    Profile.Allocs.clear()
    per_eval == 0 && return max_evals
    return clamp(target_records ÷ per_eval, 1, max_evals)
end

# ---------------------------------------------------------------------------
# Per-config driver
# ---------------------------------------------------------------------------
function run_config(
    config;
    modes,
    targets,
    t_final,
    seconds,
    alloc_rate,
    from_c,
    strat_name,
)
    # File stems keep the historical names for the default strategy.
    sfx = strat_name == "direct" ? "" : "__$(strat_name)"
    backend() = Structured(strategy = STRATEGIES[strat_name]())
    banner("config: $config")

    if "residual" in targets || "jacobian" in targets
        def, ws, u, du = BUILDERS[config]()

        jac_spec = JacobianSpec(backend(), def, ws)
        J = Pyrolysis.Jacobian._alloc_target(jac_spec, length(u))
        compute_jacobian!(J, jac_spec, u, def, 0.0, ws)

        residual_loop! = n -> begin
            for _ = 1:n
                residual!(du, u, def, 0.0, ws)
            end
        end
        jacobian_loop! = n -> begin
            for _ = 1:n
                fill!(J, 0.0)
                compute_jacobian!(J, jac_spec, u, def, 0.0, ws)
            end
        end

        hot_ops = Pair{String,Function}[]
        "residual" in targets && push!(hot_ops, "residual" => residual_loop!)
        "jacobian" in targets && push!(hot_ops, "jacobian" => jacobian_loop!)

        if "summary" in modes
            section("hot-op summary (per evaluation)")
            for (name, loop!) in hot_ops
                summarize_hot_op(name, loop!, calibrate_reps(loop!; seconds))
            end
        end
        for (name, loop!) in hot_ops
            n = calibrate_reps(loop!; seconds)
            if "cpu" in modes
                profile_cpu(
                    "$config/$name",
                    joinpath(OUTDIR, "$(config)__$(name)$(sfx)__cpu.pb.gz"),
                    () -> loop!(n);
                    delay = 1e-4,
                    from_c = from_c,
                )
            end
            if "allocs" in modes && enough_memory("$name allocation profile")
                # Hot ops should be allocation-free → sample everything, but
                # size the window from the measured records-per-eval rate.
                n_alloc = calibrate_alloc_evals(loop!)
                profile_allocs(
                    "$config/$name",
                    joinpath(OUTDIR, "$(config)__$(name)$(sfx)__allocs.pb.gz"),
                    () -> loop!(n_alloc);
                    sample_rate = 1.0,
                )
            end
        end
    end

    if "solve" in targets
        section("end-to-end solve (tspan 0–$(t_final) s)")
        # Warm-up: compiles the solve path and proves the case integrates at
        # all; skip the expensive profiled runs if it cannot.
        ok = _guarded("solve warm-up") do
            wp, wkw = build_solve_case(config; t_final = min(t_final, 0.5), strat_name)
            solve(wp; wkw...)
        end
        if ok && "summary" in modes
            _guarded(() -> summarize_solve(config; t_final, strat_name), "solve summary")
        end
        if ok && "cpu" in modes
            problem, kwargs = build_solve_case(config; t_final, strat_name)
            _guarded("solve CPU profile") do
                profile_cpu(
                    "$config/solve",
                    joinpath(OUTDIR, "$(config)__solve$(sfx)__cpu.pb.gz"),
                    () -> solve(problem; kwargs...);
                    delay = 1e-3,
                    from_c = from_c,
                )
            end
        end
        if ok && "allocs" in modes && enough_memory("solve allocation profile")
            problem, kwargs = build_solve_case(config; t_final, strat_name)
            _guarded("solve allocation profile") do
                profile_allocs(
                    "$config/solve",
                    joinpath(OUTDIR, "$(config)__solve$(sfx)__allocs.pb.gz"),
                    () -> solve(problem; kwargs...);
                    sample_rate = alloc_rate,
                )
            end
        end
    end
end

# ---------------------------------------------------------------------------
function main()
    opts = parse_cli(ARGS)
    configs = _expand(opts["config"], sort(collect(keys(BUILDERS))))
    modes = _expand(opts["mode"], ["summary", "cpu", "allocs"])
    targets = _expand(opts["target"], ["residual", "jacobian", "solve"])
    t_final = parse(Float64, opts["t-final"])
    seconds = parse(Float64, opts["seconds"])
    alloc_rate = parse(Float64, opts["alloc-rate"])
    from_c = parse(Bool, opts["from-c"])
    strat_name = opts["strategy"]

    for c in configs
        haskey(BUILDERS, c) ||
            error("unknown config $c — choose from $(join(keys(BUILDERS), ", "))")
    end
    haskey(STRATEGIES, strat_name) ||
        error("unknown strategy $strat_name — choose from $(join(keys(STRATEGIES), ", "))")

    mkpath(OUTDIR)
    println(
        "Pyrolysis.jl deep profiler — julia $(VERSION), $(Threads.nthreads()) thread(s)",
    )
    println(
        "configs: ",
        join(configs, ", "),
        " | modes: ",
        join(modes, ", "),
        " | targets: ",
        join(targets, ", "),
        " | strategy: ",
        strat_name,
    )

    for config in configs
        run_config(config; modes, targets, t_final, seconds, alloc_rate, from_c, strat_name)
        GC.gc()   # drop per-config workspaces/solutions before the next build
    end

    if "cpu" in modes || "allocs" in modes
        println("\nInspect any .pb.gz interactively:")
        println(
            "  julia --project=scripts/profile -e 'using PProf; " *
            "PProf.refresh(file=\"scripts/profile/output/<file>\"); readline()'",
        )
    end
end

main()
