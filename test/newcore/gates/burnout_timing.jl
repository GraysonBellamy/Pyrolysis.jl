# Gate item 8 of the Phase 2 plan (§7.3) as M2 measures it (§5.7): the wall
# time of the burnout path against the old core's at matched accuracy, and
# where each solve's time goes.
#
# Decks, at 60 PMMA cells, both cores in this process:
#
#   V-BO-02, fixed mesh, surface absorption, sixty events, no correction;
#   V-BO-01, moving mesh, in-depth absorption, with the prototype's correction
#   on the finite-volume core (scripts/mesh_phase2_prototype/common.jl: τ of
#   0.2 s per cell, a band of ten thresholds), which ends the PMMA layer in one
#   invocation; the old core solves the case as it stands.
#
# The prototype's overrides stay in the process, so the script runs on its own.
#
# Accuracy. A cell deactivates at the end of the step on which it crosses the
# threshold, so an event's time is set by that step, which the tolerance hardly
# moves, and the temperatures after an event follow its time: over the whole
# run the achieved error of either core does not fall with the tolerance. The
# error is therefore the largest temperature difference over the cells and the
# saved times before the first event of any solve of the deck (the window),
# each core against its own solve at reltol 1e-9 (abstol a hundredth of
# reltol throughout); the reference spread is that reference against one at
# reltol 1e-10. A new-core configuration's matched tolerance is the loosest of
# the ladder whose window error is at or below the old core's at its default
# (reltol 1e-6), and its ratio is its warm time there over the old default's
# (warm medians of three solves; a configuration slower than fifteen times the
# old default is timed once and not taken to tighter tolerances).
# Beside it: the error over the whole run at the saved times with the same live
# cells, the largest event-time shift over the reference's median interval
# (both as `case_accuracy.jl` measures them), the L1 error of the mass loss
# rate history relative to the reference's integral, and the case's peak mass
# loss rate, burnout time and mass ledger.
#
# Where the time goes: the sampling profiler over two warm solves at the
# matched tolerance, each sample given to the first of: the Jacobian (colored
# or structured), the domain test, the deactivation events, the factorization,
# the linear solves, the right-hand side, the formation of the iteration
# matrix, the per-step observers; what is left is the integrator's own work.
# Unit costs at the final state of a solve that ends at 200 s: one right-hand side, one Jacobian,
# one KLU refactorization of `J − M/(0.1 s)` with a solve, and one solve.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_timing.jl
# Writes test/newcore/gates/measurements/burnout_timing.{csv,md}.

include(
    joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase2_prototype", "common.jl"),
)
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using Profile
using LinearSolve: LinearSolve, LinearProblem, KLUFactorization
using Pyrolysis.Integration: KenCarp4, AutoFiniteDiff

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const N_PMMA = 60
const REPEATS = 3
const DEFAULT_RELTOL = 1e-6
const LADDER = (1e-5, 3e-6, 1e-6, 3e-7, 1e-7)
const REFERENCE = 1e-9
const SPREAD = 1e-10
const UNIT_TIME = 200.0
const SLOW = 15.0
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

# The integrator configurations of the finite-volume core: the default and the
# Newton-based candidates with the restart of the old core after an event.
const CONFIGURATIONS = (
    ("Rodas5P (default)", ls -> Ig.default_policy(ls)),
    ("KenCarp4 raw (kencarp4_policy)", ls -> Ig.kencarp4_policy(ls)),
    (
        "KenCarp4 smoothed, forced restart, flag cleared",
        ls -> Ig.IntegratorPolicy(
            KenCarp4(autodiff = AutoFiniteDiff(), linsolve = ls),
            :forced_dt,
            :clear,
            Ig.DT_AFTER_EVENT,
        ),
    ),
)

# ----------------------------------------------------------------------------
# The decks
# ----------------------------------------------------------------------------

deck_problem(deck; t_end = 500.0) =
    deck == "V-BO-02" ?
    Base.invokelatest(
        _vbo_problem,
        N_PMMA;
        substrate_alpha = 0.0,
        pmma_alpha = 0.0,
        t_end,
    ) : Base.invokelatest(_vbo_problem, N_PMMA; t_end)
const OLD_KEYWORDS = (
    "V-BO-02" => (
        use_ale = false,
        handle_depletion = true,
        depletion_threshold = 0.05,
        min_cells = 2,
        radiation_model = SURFACE_ABSORPTION,
    ),
    "V-BO-01" => (
        use_ale = true,
        handle_depletion = true,
        depletion_threshold = 0.05,
        min_cells = 2,
        min_thickness = 1e-3,
        radiation_model = BEER_LAMBERT,
    ),
)

old_solve(deck, reltol; problem = deck_problem(deck)) =
    Logging.with_logger(Logging.NullLogger()) do
        solve(problem; progress = false, Dict(OLD_KEYWORDS)[deck]..., _tol(reltol)...)
    end

# The set-up problem of a deck on the finite-volume core.
function core_problem(deck, reltol; problem = deck_problem(deck))
    if deck == "V-BO-02"
        return Ig.setup(
            problem;
            use_ale = false,
            deactivation = true,
            radiation_model = SURFACE_ABSORPTION,
            time_scale = 1e-3,
            _tol(reltol)...,
        )
    end
    cp, _ = build_vbo(
        N_PMMA;
        builder = n -> problem,
        τ = 0.2,
        κ = 10.0,
        case_kw = (; radiation_model = BEER_LAMBERT),
        _tol(reltol)...,
    )
    return cp
end

function new_solve(cp, configuration)
    policy = configuration(Ig.default_linear_solver(cp.def.layout))
    return Ig.integrate(cp; integrator_policy = policy, min_thickness = 1e-3)
end

# ----------------------------------------------------------------------------
# Error measures
# ----------------------------------------------------------------------------

_event_times(sol) =
    hasproperty(sol, :events) ? [e.t for e in sol.events] : collect(sol.extras.park_times)
_live(sol) = sol.extras.live_cells

# The largest temperature difference over the cells at the saved times before `t1`.
function window_error(sol, ref, t1)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    err = 0.0
    for (j, t) in enumerate(sol.t)
        t < t1 || continue
        k = get(index, t, 0)
        k == 0 && continue
        for i in axes(sol.T, 1)
            a, b = sol.T[i, j], ref.T[i, k]
            (isfinite(a) && isfinite(b)) && (err = max(err, abs(a - b)))
        end
    end
    return err
end

# The largest temperature difference at the saved times with the same live cells.
function run_error(sol, ref)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    ls, lr = _live(sol), _live(ref)
    err = 0.0
    for (j, t) in enumerate(sol.t)
        k = get(index, t, 0)
        (k == 0 || ls[j] != lr[k]) && continue
        for i in axes(sol.T, 1)
            a, b = sol.T[i, j], ref.T[i, k]
            (isfinite(a) && isfinite(b)) && (err = max(err, abs(a - b)))
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
    (occursin("refactor", f) || occursin("klu_l_factor", f) || f in ("klu!", "lu!")) &&
        return :factorization
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

function _klu_costs(J, M)
    W = SparseMatrixCSC(J - M / 0.1)
    b = rand(size(W, 1))
    cache = LinearSolve.init(LinearProblem(W, b), KLUFactorization())
    LinearSolve.solve!(cache)
    factor = unit_time(() -> (cache.A = W; LinearSolve.solve!(cache)))
    solve_alone = unit_time(() -> (cache.b = b; LinearSolve.solve!(cache)))
    return factor, solve_alone, nnz(W)
end

# One right-hand side, Jacobian, refactorization with a solve and solve of the
# integrator's functions at the final state of a solve that ends at `UNIT_TIME`
# (the cells deactivated by then are those of its layout).
function unit_costs(sol)
    raw = sol.raw_solution
    f = raw.prob.f
    u, t = copy(raw.u[end]), raw.t[end]
    du = similar(u)
    J = copy(f.jac_prototype)
    rhs = unit_time(() -> f(du, u, raw.prob.p, t))
    jac = unit_time(() -> f.jac(J, u, raw.prob.p, t))
    M = f.mass_matrix isa UniformScaling ? sparse(1.0I, size(J)...) : sparse(f.mass_matrix)
    factor, solve_alone, nz = _klu_costs(J, M)
    return (; rhs, jac, factor, solve_alone, states = length(u), nnz = nz)
end

# ----------------------------------------------------------------------------
# The measurements
# ----------------------------------------------------------------------------

# The warm median of `REPEATS` solves after one warm-up; a warm-up slower than
# `cap` [s] is not repeated and its own time is reported.
function timed(f; cap = Inf)
    sol = f()
    first = @elapsed sol = f()
    first > cap && return sol, first
    times = [first]
    for _ = 2:REPEATS
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

function deck_rows(deck)
    println("[$deck] references")
    flush(stdout)
    old_ref = old_solve(deck, REFERENCE)
    old_spread = old_solve(deck, SPREAD)
    new_ref = new_solve(core_problem(deck, REFERENCE), CONFIGURATIONS[1][2])
    new_spread = new_solve(core_problem(deck, SPREAD), CONFIGURATIONS[1][2])
    runs = Any[]
    for reltol in LADDER
        sol, wall = timed(() -> old_solve(deck, reltol))
        push!(
            runs,
            (;
                core = "old",
                run = "old core default (KenCarp4)",
                reltol,
                sol,
                wall,
                ref = old_ref,
            ),
        )
    end
    # A configuration slower than `SLOW` times the old default is not taken to
    # tighter tolerances.
    cap = SLOW * only(r.wall for r in runs if r.reltol == DEFAULT_RELTOL)
    for (name, configuration) in CONFIGURATIONS
        for reltol in LADDER
            cp = core_problem(deck, reltol)
            sol, wall = timed(() -> new_solve(cp, configuration); cap)
            push!(runs, (; core = "new", run = name, reltol, sol, wall, ref = new_ref))
            println("    $(name) $(reltol): $(sol.retcode) $(round(wall, digits = 3)) s")
            flush(stdout)
            wall > cap && break
        end
    end
    solves = vcat([r.sol for r in runs], [old_ref, old_spread, new_ref, new_spread])
    t1 = minimum(s -> minimum(_event_times(s); init = Inf), solves)
    spread = (
        old = window_error(old_spread, old_ref, t1),
        new = window_error(new_spread, new_ref, t1),
    )
    old_default = only(filter(r -> r.core == "old" && r.reltol == DEFAULT_RELTOL, runs))
    target = window_error(old_default.sol, old_ref, t1)
    rows = Any[]
    for r in runs
        m = Base.invokelatest(_vbo_measures, r.sol)
        push!(
            rows,
            (;
                deck,
                core = r.core,
                run = r.run,
                reltol = r.reltol,
                retcode = string(r.sol.retcode),
                window_error = window_error(r.sol, r.ref, t1),
                run_error = run_error(r.sol, r.ref),
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
    matched = Any[only(filter(r -> r.core == "old" && r.reltol == DEFAULT_RELTOL, rows))]
    for name in first.(CONFIGURATIONS)
        candidates = filter(
            r -> r.run == name && r.window_error <= target && r.retcode == "Success",
            rows,
        )
        isempty(candidates) &&
            (push!(matched, (; deck, run = name, reltol = NaN)); continue)
        best = argmax(r -> r.reltol, candidates)
        push!(matched, best)
    end
    anatomy_rows = Any[]
    for m in matched
        isnan(m.reltol) && continue
        f =
            m.run == old_default.run ? () -> old_solve(deck, m.reltol) :
            let cp = core_problem(deck, m.reltol), c = Dict(CONFIGURATIONS)[m.run]

                () -> new_solve(cp, c)
            end
        shares = anatomy(f)
        short = deck_problem(deck; t_end = UNIT_TIME)
        unit_sol =
            m.run == old_default.run ? old_solve(deck, m.reltol; problem = short) :
            new_solve(
                core_problem(deck, m.reltol; problem = short),
                Dict(CONFIGURATIONS)[m.run],
            )
        push!(
            anatomy_rows,
            (; deck, run = m.run, reltol = m.reltol, shares, unit = unit_costs(unit_sol)),
        )
        println("    anatomy $(m.run) $(m.reltol) done")
        flush(stdout)
    end
    return (; deck, t1, spread, target, rows, matched, anatomy = anatomy_rows)
end

load_start = Sys.loadavg()
started = time()
results = [deck_rows(deck) for deck in ("V-BO-02", "V-BO-01")]
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
    println(io, "# Gate item 8 at M2: the burnout path against the old core\n")
    println(
        io,
        "Produced by `test/newcore/gates/burnout_timing.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "$(N_PMMA) PMMA cells; V-BO-01 on the finite-volume core with the prototype's correction ",
        "(τ 0.2 s per cell, band 10 ε_m). Window error: the largest temperature difference [K] over ",
        "the cells at the saved times before the first event of any solve of the deck, against the ",
        "same core at reltol $(REFERENCE); run error: over the whole run at the saved times with the ",
        "same live cells; shift: the largest event-time difference over the reference's median ",
        "interval; MLR: the L1 error of the mass loss rate history over the reference's integral. ",
        "Warm medians of $(REPEATS) solves (a solve slower than $(Int(SLOW)) times the old default ",
        "is timed once and its configuration not taken to tighter tolerances); ratio: over the ",
        "old core's default.\n",
    )
    for res in results
        println(io, "## $(res.deck)\n")
        println(
            io,
            "Window: the saved times before $(round(res.t1, digits = 1)) s. Reference spread ",
            "(reltol $(SPREAD) against $(REFERENCE)): old $(_g(res.spread.old)) K, new ",
            "$(_g(res.spread.new)) K. Target: the old core's window error at reltol ",
            "$(DEFAULT_RELTOL), $(_g(res.target)) K.\n",
        )
        println(
            io,
            "| run | reltol | retcode | window error [K] | run error [K] | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |",
        )
        println(io, "|" * repeat("---|", 17))
        for r in res.rows
            println(
                io,
                "| $(r.run) | $(r.reltol) | $(r.retcode) | $(_g(r.window_error)) | $(_g(r.run_error)) | ",
                "$(_g(r.event_shift)) | $(_g(r.mlr_error)) | $(_g(r.peak; d = 5)) | $(_g(r.t_bo; d = 4)) | ",
                "$(_g(r.ledger)) | $(_g(r.wall)) | $(_g(r.ratio)) | $(r.steps) | $(r.rejected) | $(r.rhs) | ",
                "$(r.jacobians) | $(r.nonlinear_failures) |",
            )
        end
        println(
            io,
            "\nMatched accuracy (the loosest ladder tolerance at or below the target):\n",
        )
        println(io, "| run | reltol | window error [K] | warm [s] | ratio |")
        println(io, "|---|---|---|---|---|")
        for m in res.matched
            isnan(m.reltol) ? println(io, "| $(m.run) | none on the ladder | | | |") :
            println(
                io,
                "| $(m.run) | $(m.reltol) | $(_g(m.window_error)) | $(_g(m.wall)) | $(_g(m.ratio)) |",
            )
        end
        println(
            io,
            "\nWhere the time goes at the matched tolerance (shares of the profiler's samples) and unit costs at $(UNIT_TIME) s [µs]:\n",
        )
        println(
            io,
            "| run | reltol | " *
            join(string.(PARTS), " | ") *
            " | RHS | Jacobian | refactor + solve | solve | states | nnz(W) |",
        )
        println(io, "|" * repeat("---|", 2 + length(PARTS) + 6))
        for a in res.anatomy
            u = a.unit
            println(
                io,
                "| $(a.run) | $(a.reltol) | " *
                join((_g(a.shares[p]; d = 2) for p in PARTS), " | ") *
                " | $(_g(1e6u.rhs)) | $(_g(1e6u.jac)) | $(_g(1e6u.factor)) | $(_g(1e6u.solve_alone)) | $(u.states) | $(u.nnz) |",
            )
        end
        println(io)
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
