# Gate item 8 of the Phase 2 plan (§7.3), as M4 measures it: the wall time
# of the burnout path against the old core's at matched accuracy, and the cost
# of a correction present but inactive.
#
# The burnout decks, both cores in this process: V-BO-01, V-BO-03 and V-BO-04
# at 30, 60 and 120 cells (PMMA cells on the first two). The finite-volume
# core solves each as the suite routes it, with the correction `setup`'s
# installation rule gives it: the equidistribution correction at its
# defaults, the PMMA–Kaowool interface Lagrangian, and on V-BO-04's uniform
# column a correction that stays exactly zero. The old core solves each with
# the case's keywords.
#
# Accuracy. A cell deactivates at the end of the step on which it crosses the
# threshold, so an event's time is set by that step, which the tolerance hardly
# moves, and the solution after an event follows its time: over the whole run
# the achieved error of either core does not fall with the tolerance. The
# error is therefore measured over the saved times before the first event of
# any solve of the deck (the window), each core against its own solve at
# reltol 1e-9 (abstol a hundredth of reltol throughout), as the RMS of the
# differences: of the temperatures over the cells on V-BO-01 and V-BO-03 [K],
# of the thickness relative to the reference's on V-BO-04, which stays at its
# initial temperature. The largest difference is reported beside it. The
# reference spread is that reference against one at reltol 1e-10. The new
# core's matched tolerance is the loosest of the ladder whose window RMS is at
# or below the old core's at the case's tolerance (the old default), and its
# ratio is its warm time there over the old default's (warm medians of three
# solves, of nine for the old default, whose time varies the most; a solve
# slower than fifteen times the old default is timed once and its core not
# taken to tighter tolerances).
# Beside it: the error over the whole run at the saved times with the same live
# cells, the largest event-time shift over the reference's median interval
# (both as `case_accuracy.jl` measures them), the L1 error of the mass loss
# rate history relative to the reference's integral, and the case's peak mass
# loss rate, burnout time and mass ledger.
#
# Where the time goes, on V-BO-01 at 60 cells: the sampling profiler over two
# warm solves at the matched tolerance, each sample given to the first of: the
# Jacobian (colored or structured), the domain test, the deactivation events,
# the factorization, the linear solves, the right-hand side, the formation of
# the iteration matrix, the per-step observers; what is left is the
# integrator's own work. Unit costs at the final state of a solve that ends at
# 200 s: one right-hand side and one Jacobian, with the state count and the
# Jacobian's stored entries.
#
# The inactive correction: the moving-mesh audit decks of Phase 1's gate item 2
# (`ale` with the solid's gas-transfer coefficient at 1e-6 m²/s and
# `ale_darcy`, at 25 and 100 cells) and V-CPL-02 and V-CPL-03 at their CI
# resolutions, on the finite-volume core without a correction and with an
# equidistribution correction at its defaults installed, at the case's
# tolerances (the audit decks at the old default, reltol 1e-6), beside the old
# core at the same tolerances: warm medians of three, the steps, and the
# largest correction velocity over the saved times (zero when the correction
# never activated).
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_timing.jl
# Writes test/newcore/gates/measurements/burnout_timing.{csv,md}.

ENV["PYROLYSIS_CORE"] = "new"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
include(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "scripts",
        "mesh_phase0",
        "baselines",
        "jacobian_common.jl",
    ),
)
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using Profile
using SparseArrays: nnz
using Statistics: median

const Ig = Pyrolysis.Integration
const Cl = Pyrolysis.Closures

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REPEATS = 3
const OLD_DEFAULT_REPEATS = 9
const LADDER = (1e-5, 3e-6, 1e-6, 3e-7, 1e-7)
const REFERENCE = 1e-9
const SPREAD = 1e-10
const UNIT_TIME = 200.0
const SLOW = 15.0
const CELLS = (30, 60, 120)
const ANATOMY = ("V-BO-01", 60)
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

const CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-03", "V-BO-04", "V-CPL-02", "V-CPL-03")
)

# ----------------------------------------------------------------------------
# The decks
# ----------------------------------------------------------------------------

deck_problem(deck, n; t_end = nothing) =
    deck == "V-BO-04" ?
    Base.invokelatest(_vbo04_setup, n; t_end = something(t_end, 480.0)).problem :
    Base.invokelatest(
        _vbo_problem,
        n;
        substrate_alpha = deck == "V-BO-03" ? 0.0 : 8.0e6 / VBO_RHO_KW,
        t_end = something(t_end, 500.0),
    )

# The case's tolerance: the old default.
default_reltol(deck) = CASES[deck].solve_kwargs.reltol
# The window's field: the cells' temperatures, or V-BO-04's thickness.
field(deck) = deck == "V-BO-04" ? :thickness : :T

function old_solve(deck, n, reltol; problem = deck_problem(deck, n))
    kw = CASES[deck].solve_kwargs
    return Logging.with_logger(Logging.NullLogger()) do
        solve(problem; progress = false, kw..., _tol(reltol)...)
    end
end

# The set-up problem of a deck on the finite-volume core, as the suite routes it.
function core_problem(deck, n, reltol; problem = deck_problem(deck, n))
    kw = CASES[deck].solve_kwargs
    rad = haskey(kw, :radiation_model) ? (; radiation_model = kw.radiation_model) : (;)
    return Ig.setup(
        problem;
        use_ale = kw.use_ale,
        deactivation = kw.handle_depletion,
        rad...,
        _tol(reltol)...,
    )
end

new_solve(deck, cp) =
    Ig.integrate(cp; min_thickness = CASES[deck].solve_kwargs.min_thickness)

# ----------------------------------------------------------------------------
# Error measures
# ----------------------------------------------------------------------------

_event_times(sol) =
    hasproperty(sol, :events) ? [e.t for e in sol.events] : collect(sol.extras.park_times)
_live(sol) = sol.extras.live_cells

# The differences at the saved times before `t1`: the temperatures over the
# cells, or the thickness relative to the reference's.
function _window_differences(sol, ref, t1, f)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    d = Float64[]
    for (j, t) in enumerate(sol.t)
        t < t1 || continue
        k = get(index, t, 0)
        k == 0 && continue
        if f === :thickness
            a, b = sol.extras.thickness_history[j], ref.extras.thickness_history[k]
            push!(d, (a - b) / b)
        else
            for i in axes(sol.T, 1)
                a, b = sol.T[i, j], ref.T[i, k]
                (isfinite(a) && isfinite(b)) && push!(d, a - b)
            end
        end
    end
    return d
end

function window_errors(sol, ref, t1, f)
    d = _window_differences(sol, ref, t1, f)
    return (
        rms = sqrt(sum(abs2, d; init = 0.0) / max(length(d), 1)),
        largest = maximum(abs, d; init = 0.0),
    )
end

# The largest difference at the saved times with the same live cells.
function run_error(sol, ref, f)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    ls, lr = _live(sol), _live(ref)
    err = 0.0
    for (j, t) in enumerate(sol.t)
        k = get(index, t, 0)
        (k == 0 || ls[j] != lr[k]) && continue
        if f === :thickness
            b = ref.extras.thickness_history[k]
            err = max(err, abs(sol.extras.thickness_history[j] - b) / b)
        else
            for i in axes(sol.T, 1)
                a, b = sol.T[i, j], ref.T[i, k]
                (isfinite(a) && isfinite(b)) && (err = max(err, abs(a - b)))
            end
        end
    end
    return err
end

function event_shift(sol, ref)
    a, b = _event_times(sol), _event_times(ref)
    (isempty(b) || length(a) != length(b)) && return NaN
    gaps = filter(>(0), diff(unique(b)))
    return isempty(gaps) ? NaN : maximum(abs, a .- b) / median(gaps)
end

function mlr_error(sol, ref)
    a, b = sol.extras.gas_generation_rate, ref.extras.gas_generation_rate
    n = min(length(a), length(b))
    t = sol.t[1:n]
    trap(f) = sum(0.5 * (f(j) + f(j + 1)) * (t[j+1] - t[j]) for j = 1:(n-1))
    return trap(j -> abs(a[j] - b[j])) / trap(j -> abs(b[j]))
end

# ----------------------------------------------------------------------------
# Cost anatomy
# ----------------------------------------------------------------------------

const PARTS = (
    :jacobian,
    :domain,
    :events,
    :factorization,
    :linsolve,
    :rhs,
    :wform,
    :observers,
    :integrator,
)
const _RHS_DIRS = (
    "/src/Assembly/",
    "/src/Kernels/",
    "/src/Closures/",
    "/src/States/",
    "/src/Residual/",
    "/src/Physics/",
    "/src/Discretization/",
    "/src/BoundaryConditions/",
    "/src/Properties/",
)

function _part(frame)
    f, file = string(frame.func), string(frame.file)
    (occursin("/src/ColoredJacobian/", file) || occursin("/src/Jacobian/", file)) &&
        return :jacobian
    occursin("Integration/checks.jl", file) && return :domain
    (occursin("/src/Deactivation/", file) || occursin("/src/Adaptivity/", file)) &&
        return :events
    (
        occursin("refactor", f) ||
        occursin("klu_l_factor", f) ||
        occursin("klu_l_analyze", f) ||
        occursin("klu_l_rgrowth", f) ||
        f in ("klu!", "lu!", "_load!", "_grow!", "_walk", "_union_matrix!")
    ) && return :factorization
    (occursin("klu_l_solve", f) || f == "ldiv!") && return :linsolve
    any(d -> occursin(d, file), _RHS_DIRS) && return :rhs
    f in ("jacobian2W!", "calc_W!", "calc_rosenbrock_differentiation!") && return :wform
    (
        f in ("_observe!", "cap_step!", "apply_observer_flag!") ||
        occursin("/src/Solver/callbacks.jl", file)
    ) && return :observers
    return :integrator
end

# The share of the profiler's samples over two warm solves in each part.
function anatomy(run)
    run()
    Profile.init(n = 20_000_000, delay = 0.0005)
    Profile.clear()
    for _ = 1:2
        @profile run()
    end
    data = Profile.fetch(include_meta = false)
    parts = Dict{UInt64,Vector{Symbol}}()
    counts = Dict(p => 0 for p in PARTS)
    seen = Set{Symbol}()
    for ip in data
        if ip == 0
            part = :integrator
            for p in PARTS
                p in seen && (part = p; break)
            end
            counts[part] += 1
            empty!(seen)
            continue
        end
        for p in get!(() -> map(_part, Profile.lookup(ip)), parts, ip)
            p === :integrator || push!(seen, p)
        end
    end
    total = max(sum(values(counts)), 1)
    return Dict(p => counts[p] / total for p in PARTS)
end

unit_time(f; samples = 101) = (f(); median([@elapsed(f()) for _ = 1:samples]))

# One right-hand side and one Jacobian of the integrator's functions at the
# final state of a solve that ends at `UNIT_TIME` (the cells deactivated by
# then are those of its layout).
function unit_costs(sol)
    raw = sol.raw_solution
    f = raw.prob.f
    u, t = copy(raw.u[end]), raw.t[end]
    du = similar(u)
    J = copy(f.jac_prototype)
    rhs = unit_time(() -> f(du, u, raw.prob.p, t))
    jac = unit_time(() -> f.jac(J, u, raw.prob.p, t))
    return (; rhs, jac, states = length(u), nnz = nnz(J))
end

# ----------------------------------------------------------------------------
# The measurements
# ----------------------------------------------------------------------------

# The warm median of `repeats` solves after one warm-up; a warm-up slower than
# `cap` [s] is not repeated and its own time is reported.
function timed(f; cap = Inf, repeats = REPEATS)
    sol = f()
    first = @elapsed sol = f()
    first > cap && return sol, first
    times = [first]
    for _ = 2:repeats
        GC.gc()
        push!(times, @elapsed sol = f())
    end
    return sol, median(times)
end

_counts(sol) =
    hasproperty(sol, :statistics) ?
    (
        steps = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        rhs = sol.statistics.rhs_evaluations,
        jacobians = sol.statistics.jacobians,
        nonlinear_failures = sol.statistics.nonlinear_failures,
    ) :
    (
        steps = sol.raw_solution.stats.naccept,
        rejected = sol.raw_solution.stats.nreject,
        rhs = sol.raw_solution.stats.nf,
        jacobians = sol.raw_solution.stats.njacs,
        nonlinear_failures = sol.raw_solution.stats.nnonlinconvfail,
    )

_measures(deck, sol) =
    deck == "V-BO-04" ? (; peak = NaN, t_bo = NaN, ledger = NaN) :
    Base.invokelatest(_vbo_measures, sol)

function deck_rows(deck, n)
    println("[$deck $n] references")
    flush(stdout)
    f = field(deck)
    default = default_reltol(deck)
    old_ref = old_solve(deck, n, REFERENCE)
    old_spread = old_solve(deck, n, SPREAD)
    new_ref = new_solve(deck, core_problem(deck, n, REFERENCE))
    new_spread = new_solve(deck, core_problem(deck, n, SPREAD))
    runs = Any[]
    for reltol in LADDER
        repeats = reltol == default ? OLD_DEFAULT_REPEATS : REPEATS
        sol, wall = timed(() -> old_solve(deck, n, reltol); repeats)
        push!(runs, (; core = "old", reltol, sol, wall, ref = old_ref))
    end
    # A core slower than `SLOW` times the old default is not taken to tighter
    # tolerances.
    cap = SLOW * only(r.wall for r in runs if r.reltol == default)
    for reltol in LADDER
        cp = core_problem(deck, n, reltol)
        sol, wall = timed(() -> new_solve(deck, cp); cap)
        push!(runs, (; core = "new", reltol, sol, wall, ref = new_ref))
        println("    new $(reltol): $(sol.retcode) $(round(wall, digits = 3)) s")
        flush(stdout)
        wall > cap && break
    end
    solves = vcat([r.sol for r in runs], [old_ref, old_spread, new_ref, new_spread])
    t1 = minimum(s -> minimum(_event_times(s); init = Inf), solves)
    spread = (
        old = window_errors(old_spread, old_ref, t1, f).rms,
        new = window_errors(new_spread, new_ref, t1, f).rms,
    )
    old_default = only(filter(r -> r.core == "old" && r.reltol == default, runs))
    target = window_errors(old_default.sol, old_ref, t1, f).rms
    rows = Any[]
    for r in runs
        m = _measures(deck, r.sol)
        w = window_errors(r.sol, r.ref, t1, f)
        push!(
            rows,
            (;
                deck,
                cells = n,
                core = r.core,
                reltol = r.reltol,
                retcode = string(r.sol.retcode),
                window_rms = w.rms,
                window_largest = w.largest,
                run_error = run_error(r.sol, r.ref, f),
                event_shift = event_shift(r.sol, r.ref),
                mlr_error = mlr_error(r.sol, r.ref),
                peak = m.peak,
                t_bo = m.t_bo,
                ledger = m.ledger,
                wall = r.wall,
                ratio = r.wall / old_default.wall,
                _counts(r.sol)...,
            ),
        )
    end
    candidates = filter(
        r -> r.core == "new" && r.window_rms <= target && r.retcode == "Success",
        rows,
    )
    matched = isempty(candidates) ? nothing : argmax(r -> r.reltol, candidates)
    old_row = only(filter(r -> r.core == "old" && r.reltol == default, rows))
    anatomy_rows = Any[]
    if (deck, n) == ANATOMY && matched !== nothing
        for (core, reltol) in (("old", default), ("new", matched.reltol))
            run =
                core == "old" ? () -> old_solve(deck, n, reltol) :
                let cp = core_problem(deck, n, reltol)

                    () -> new_solve(deck, cp)
                end
            shares = anatomy(run)
            short = deck_problem(deck, n; t_end = UNIT_TIME)
            unit_sol =
                core == "old" ? old_solve(deck, n, reltol; problem = short) :
                new_solve(deck, core_problem(deck, n, reltol; problem = short))
            push!(anatomy_rows, (; core, reltol, shares, unit = unit_costs(unit_sol)))
            println("    anatomy $(core) $(reltol) done")
            flush(stdout)
        end
    end
    return (;
        deck,
        cells = n,
        default,
        t1,
        spread,
        target,
        rows,
        old_row,
        matched,
        anatomy = anatomy_rows,
    )
end

# The inactive correction on the decks without burnout.
function inactive_rows()
    decks = Any[]
    for (kind, n, λ) in (
        (:ale, 25, 1e-6),
        (:ale, 100, 1e-6),
        (:ale_darcy, 25, nothing),
        (:ale_darcy, 100, nothing),
    )
        problem, kw = problem_for(kind, n; solid_λ = λ)
        push!(
            decks,
            (;
                name = "$(kind) $(n)",
                problem,
                old = (; use_ale = kw.use_ale, radiation_model = kw.radiation_model),
                new = (; use_ale = kw.use_ale, radiation_model = kw.radiation_model),
                tol = _tol(1e-6),
                min_thickness = 1e-6,
            ),
        )
    end
    for id in ("V-CPL-02", "V-CPL-03")
        case = CASES[id]
        kw = case.solve_kwargs
        old = (; (k => v for (k, v) in pairs(kw) if !(k in (:abstol, :reltol)))...)
        new = (; use_ale = kw.use_ale, deactivation = get(kw, :handle_depletion, false))
        push!(
            decks,
            (;
                name = id,
                problem = Base.invokelatest(case.build, case.ci.n_cells),
                old,
                new,
                tol = (abstol = kw.abstol, reltol = kw.reltol),
                min_thickness = kw.min_thickness,
            ),
        )
    end
    rows = Any[]
    for d in decks
        old_sol, old_wall = timed(() -> Logging.with_logger(Logging.NullLogger()) do
            solve(d.problem; progress = false, d.old..., d.tol...)
        end)
        results = Any[]
        for correction in (Cl.NoCorrection(), Cl.EquidistributionCorrection())
            cp = Ig.setup(d.problem; d.new..., correction, d.tol...)
            sol, wall = timed(() -> Ig.integrate(cp; min_thickness = d.min_thickness))
            c = filter(isfinite, sol.correction_velocity)
            push!(
                results,
                (;
                    sol,
                    wall,
                    largest = isempty(c) ? 0.0 : maximum(abs, c),
                    steps = sol.statistics.accepted,
                ),
            )
        end
        plain, corrected = results
        push!(
            rows,
            (;
                deck = d.name,
                reltol = d.tol.reltol,
                old_wall,
                old_steps = _counts(old_sol).steps,
                plain_wall = plain.wall,
                plain_steps = plain.steps,
                corrected_wall = corrected.wall,
                corrected_steps = corrected.steps,
                largest = corrected.largest,
                difference = maximum(
                    abs,
                    filter(isfinite, corrected.sol.T .- plain.sol.T);
                    init = 0.0,
                ),
            ),
        )
        println("    inactive $(d.name): $(rows[end])")
        flush(stdout)
    end
    return rows
end

load_start = Sys.loadavg()
started = time()
results = [deck_rows(deck, n) for deck in ("V-BO-01", "V-BO-03", "V-BO-04") for n in CELLS]
inactive = inactive_rows()
load_end = Sys.loadavg()

open(joinpath(OUT, "burnout_timing.csv"), "w") do io
    cols = keys(first(results[1].rows))
    println(io, join(cols, ","))
    for res in results, r in res.rows
        println(
            io,
            join((v isa AbstractString ? "\"$v\"" : string(v) for v in values(r)), ","),
        )
    end
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "burnout_timing.md"), "w") do io
    println(io, "# Gate item 8 at M4: the burnout path against the old core\n")
    println(
        io,
        "Produced by `test/newcore/gates/burnout_timing.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "The finite-volume core's default integrator (`Rodas5P`, heuristic restart); V-BO-01 and ",
        "V-BO-03 with the equidistribution correction at its defaults. Window RMS and largest: the ",
        "RMS and the largest of the differences over the saved times before the first event of any ",
        "solve of the deck, against the same core at reltol $(REFERENCE): of the temperatures over ",
        "the cells [K] on V-BO-01 and V-BO-03, of the thickness relative to the reference's on ",
        "V-BO-04; the new core is matched on the RMS. Run error: the largest difference over the ",
        "whole run at the saved times with the same live cells; shift: the largest event-time ",
        "difference over the reference's median interval; MLR: the L1 error of the mass loss rate ",
        "history over the reference's integral. Warm medians of $(REPEATS) solves, of ",
        "$(OLD_DEFAULT_REPEATS) for the old default (a solve slower than $(Int(SLOW)) times the old ",
        "default is timed once and its core not taken to tighter tolerances); ratio: over the old ",
        "core's time at the case's tolerance.\n",
    )
    println(io, "## Matched accuracy\n")
    println(
        io,
        "| deck | cells | window before [s] | spread old, new | target | old default [s] | new matched reltol | window RMS | window largest | new [s] | ratio | equal-tolerance ratio |",
    )
    println(io, "|" * repeat("---|", 12))
    for res in results
        m = res.matched
        k = findfirst(r -> r.core == "new" && r.reltol == res.default, res.rows)
        eq = k === nothing ? NaN : res.rows[k].ratio
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(_g(res.t1; d = 4)) | $(_g(res.spread.old)), ",
            "$(_g(res.spread.new)) | $(_g(res.target)) | $(_g(res.old_row.wall)) | ",
            m === nothing ? "none on the ladder | | | | | " :
            "$(m.reltol) | $(_g(m.window_rms)) | $(_g(m.window_largest)) | $(_g(m.wall)) | " *
            "$(_g(m.ratio)) | ",
            "$(_g(eq)) |",
        )
    end
    println(io, "\n## The ladders\n")
    println(
        io,
        "| deck | cells | core | reltol | retcode | window RMS | window largest | run error | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |",
    )
    println(io, "|" * repeat("---|", 20))
    for res in results, r in res.rows
        println(
            io,
            "| $(r.deck) | $(r.cells) | $(r.core) | $(r.reltol) | $(r.retcode) | ",
            "$(_g(r.window_rms)) | $(_g(r.window_largest)) | $(_g(r.run_error)) | ",
            "$(_g(r.event_shift)) | $(_g(r.mlr_error)) | $(_g(r.peak; d = 5)) | ",
            "$(_g(r.t_bo; d = 4)) | $(_g(r.ledger)) | $(_g(r.wall)) | $(_g(r.ratio)) | ",
            "$(r.steps) | $(r.rejected) | $(r.rhs) | $(r.jacobians) | $(r.nonlinear_failures) |",
        )
    end
    for res in results
        isempty(res.anatomy) && continue
        println(
            io,
            "\n## Where the time goes: $(res.deck) at $(res.cells) cells\n\nShares of the ",
            "profiler's samples at the old default and at the new core's matched tolerance, and ",
            "unit costs at $(UNIT_TIME) s [µs]:\n",
        )
        println(
            io,
            "| core | reltol | " *
            join(string.(PARTS), " | ") *
            " | RHS | Jacobian | states | nnz(J) |",
        )
        println(io, "|" * repeat("---|", 2 + length(PARTS) + 4))
        for a in res.anatomy
            u = a.unit
            println(
                io,
                "| $(a.core) | $(a.reltol) | " *
                join((_g(a.shares[p]; d = 2) for p in PARTS), " | ") *
                " | $(_g(1e6u.rhs)) | $(_g(1e6u.jac)) | $(u.states) | $(u.nnz) |",
            )
        end
    end
    println(io, "\n## The inactive correction\n")
    println(
        io,
        "The finite-volume core without a correction and with an equidistribution correction at ",
        "its defaults installed, at the same tolerance; the old core beside. Largest |c|: the ",
        "largest correction velocity over the saved times [m/s] (zero: never active); ",
        "difference: the largest temperature difference between the two new-core solves [K].\n",
    )
    println(
        io,
        "| deck | reltol | old [s] (steps) | new [s] (steps) | with the correction [s] (steps) | correction / none | with the correction / old | largest |c| | difference |",
    )
    println(io, "|" * repeat("---|", 9))
    for r in inactive
        println(
            io,
            "| $(r.deck) | $(r.reltol) | $(_g(r.old_wall)) ($(r.old_steps)) | ",
            "$(_g(r.plain_wall)) ($(r.plain_steps)) | $(_g(r.corrected_wall)) ",
            "($(r.corrected_steps)) | $(_g(r.corrected_wall / r.plain_wall)) | ",
            "$(_g(r.corrected_wall / r.old_wall)) | $(_g(r.largest)) | $(_g(r.difference)) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
