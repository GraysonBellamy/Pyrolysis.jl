# Gate item 8 of the Phase 2 plan (§7.3): the wall time of the burnout path
# against the old core's at matched accuracy, and the cost of the correction
# where the installation rule installs it.
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
# taken to tighter tolerances). Item 8's verdict on that ratio: pass at 1.5
# or below, the owner's call above 1.5 and up to 2, a stop above 2 (§7.4).
# Beside it: the error over the whole run at the saved times with the same live
# cells, the largest invocation-time shift over the reference's median
# interval between invocations, or between its saved times for a single
# invocation (both as `case_accuracy.jl` measures them; the shift is
# `event_shift.jl`'s, shared with it), the L1 error of the mass loss
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
# The installed correction: the decks on which the installation rule installs
# it, the moving-mesh audit decks of Phase 1's gate item 2 with deactivation on
# (`ale` with the solid's gas-transfer coefficient at 1e-6 m²/s to 400 s and
# `ale_darcy` to 300 s, at 25 and 100 cells, default reltol 1e-6, the old
# core's depletion keywords at their defaults) and V-CPL-02 and V-CPL-02b at
# their CI resolutions and tolerances, measured as the burnout decks are: each
# core's ladder against its own reference, the window before the first event
# of any solve of the deck (the whole run where no solve has one), the new
# core matched on the window RMS. Beside it, the rule's solve at the default
# tolerance against the same deck with `NoCorrection()` (warm medians of
# three, the steps, the outcome and end time: without the correction a layer
# can end `:interior`) and the largest correction velocity over the saved
# times. A deck on which the rule installs nothing is bit-identical (item 5)
# and costs nothing.
#
# The cost split at the first activation, per matched row whose correction
# activates (the installed decks, V-BO-01 and V-BO-03): before `t_a`, the first
# saved time of the new core's matched solve with a nonzero correction
# velocity, the correction is exactly zero and the new core solves the
# problem without it on the wider pattern, so the ratio there is the
# integrator's and the per-step cost's; after it the gain `N²/τ` acts. Each
# core's accepted steps before and after `t_a` and its share of the time
# before it come from a read-only step recorder in a separate solve, the share
# applied to the timed median so the recorder's own cost stays out of the
# ratio; the gain's share of the new core's time is
# `(T_new,after − r_before · T_old,after) / T_new`. One diagnostic row,
# reported and not gated: `ale_darcy` at 100 cells with the layer time scaled
# by 16, whose per-node gain equals `ale_darcy` 25's at the default.
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
include(joinpath(@__DIR__, "event_shift.jl"))

using Logging
using OrdinaryDiffEq: DiscreteCallback
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

# Item 8 on a ratio at matched accuracy.
item8_verdict(ratio) = ratio <= 1.5 ? "pass" : ratio <= 2 ? "owner" : "stop"

const CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-03", "V-BO-04", "V-CPL-02", "V-CPL-02b")
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

# A deck of the tables: its name and cells, the problem, the old core's
# keywords and `setup`'s, the default relative tolerance (the old default),
# the window's field (`:T`, the cells' temperatures, or `:thickness`, the
# thickness relative to the reference's), `min_thickness`, and the problem
# ending at a given time for the unit costs (`nothing` where none are taken).
function vbo_deck(deck, n)
    kw = CASES[deck].solve_kwargs
    rad = haskey(kw, :radiation_model) ? (; radiation_model = kw.radiation_model) : (;)
    return (;
        name = deck,
        cells = n,
        problem = deck_problem(deck, n),
        old = kw,
        new = (; use_ale = kw.use_ale, deactivation = kw.handle_depletion, rad...),
        default = kw.reltol,
        field = deck == "V-BO-04" ? :thickness : :T,
        min_thickness = kw.min_thickness,
        until = t_end -> deck_problem(deck, n; t_end),
    )
end

# The decks on which the installation rule installs the correction.
function installed_decks()
    decks = Any[]
    for (kind, n, λ, t_end) in (
        (:ale, 25, 1e-6, 400.0),
        (:ale, 100, 1e-6, 400.0),
        (:ale_darcy, 25, nothing, 300.0),
        (:ale_darcy, 100, nothing, 300.0),
    )
        problem, kw =
            problem_for(kind, n; solid_λ = λ, t_end, outputs = round(Int, t_end) + 1)
        common = (; use_ale = kw.use_ale, radiation_model = kw.radiation_model)
        push!(
            decks,
            (;
                name = string(kind),
                cells = n,
                problem,
                old = (; common..., handle_depletion = true),
                new = (; common..., deactivation = true),
                default = 1e-6,
                field = :T,
                min_thickness = 1e-6,
                until = nothing,
            ),
        )
    end
    for id in ("V-CPL-02", "V-CPL-02b")
        case = CASES[id]
        kw = case.solve_kwargs
        push!(
            decks,
            (;
                name = id,
                cells = case.ci.n_cells,
                problem = Base.invokelatest(case.build, case.ci.n_cells),
                old = kw,
                new = (;
                    use_ale = kw.use_ale,
                    deactivation = get(kw, :handle_depletion, false),
                ),
                default = kw.reltol,
                field = :T,
                min_thickness = kw.min_thickness,
                until = nothing,
            ),
        )
    end
    return decks
end

old_solve(d, reltol; problem = d.problem, kw...) =
    Logging.with_logger(Logging.NullLogger()) do
        solve(problem; progress = false, d.old..., _tol(reltol)..., kw...)
    end

# The set-up problem of a deck on the finite-volume core, as the suite routes it.
core_problem(d, reltol; problem = d.problem, kw...) =
    Ig.setup(problem; d.new..., _tol(reltol)..., kw...)

new_solve(d, cp; kw...) = Ig.integrate(cp; min_thickness = d.min_thickness, kw...)

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

function mlr_error(sol, ref)
    a, b = sol.extras.gas_generation_rate, ref.extras.gas_generation_rate
    n = min(length(a), length(b))
    t = sol.t[1:n]
    trap(f) = sum(0.5 * (f(j) + f(j + 1)) * (t[j+1] - t[j]) for j = 1:(n-1))
    scale = trap(j -> abs(b[j]))
    return scale > 0 ? trap(j -> abs(a[j] - b[j])) / scale : NaN
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

_measures(d, sol) =
    d.name in ("V-BO-01", "V-BO-03") ? Base.invokelatest(_vbo_measures, sol) :
    (; peak = NaN, t_bo = NaN, ledger = NaN)

# ----------------------------------------------------------------------------
# The cost split at the first activation
# ----------------------------------------------------------------------------

# The first saved time with a nonzero correction velocity (NaN if none).
function first_activation(sol)
    c = sol.correction_velocity
    j = findfirst(j -> any(x -> isfinite(x) && !iszero(x), view(c, :, j)), axes(c, 2))
    return j === nothing ? NaN : sol.t[j]
end

# A read-only step recorder: the time and the wall clock at every accepted
# step, read by a condition that never fires.
function step_recorder()
    t, clock = Float64[], Float64[]
    callback = DiscreteCallback(
        (u, s, integrator) -> (push!(t, s); push!(clock, time()); false),
        integrator -> nothing;
        save_positions = (false, false),
    )
    return callback, t, clock
end

# The share of a solve's wall time up to `t_a` and its accepted steps before
# and after it, from the recorder in `run(callback)`, the second of two solves
# so that the first's compilation stays out.
function split_at(run, t_a)
    share, steps = NaN, (0, 0)
    for _ = 1:2
        callback, t, clock = step_recorder()
        start = time()
        run(callback)
        stop = time()
        k = count(<=(t_a), t)
        share = k == 0 ? 0.0 : (clock[k] - start) / (stop - start)
        steps = (k, length(t) - k)
    end
    return (; share, steps)
end

# The cost split of a matched solve at `t_a`: each core's time before and
# after it (the recorder's shares applied to the timed medians), the ratios
# before and after, and the gain's share of the new core's time, its time
# after `t_a` beyond the old core's at the ratio before.
function cost_split(d, matched, old_default, t_a)
    new = split_at(callback -> new_solve(d, matched.cp; callback), t_a)
    old = split_at(callback -> old_solve(d, d.default; callback), t_a)
    new_before, new_after = new.share * matched.wall, (1 - new.share) * matched.wall
    old_before, old_after = old.share * old_default.wall, (1 - old.share) * old_default.wall
    r_before = new_before / old_before
    return (;
        t_a,
        new_steps = new.steps,
        old_steps = old.steps,
        new_share = new.share,
        old_share = old.share,
        r_before,
        r_after = new_after / old_after,
        gain = (new_after - r_before * old_after) / matched.wall,
    )
end

# ----------------------------------------------------------------------------
# The decks' rows
# ----------------------------------------------------------------------------

function deck_rows(d)
    println("[$(d.name) $(d.cells)] references")
    flush(stdout)
    f = d.field
    old_ref = old_solve(d, REFERENCE)
    old_spread = old_solve(d, SPREAD)
    new_ref = new_solve(d, core_problem(d, REFERENCE))
    new_spread = new_solve(d, core_problem(d, SPREAD))
    runs = Any[]
    for reltol in LADDER
        repeats = reltol == d.default ? OLD_DEFAULT_REPEATS : REPEATS
        sol, wall = timed(() -> old_solve(d, reltol); repeats)
        push!(runs, (; core = "old", reltol, sol, wall, ref = old_ref, cp = nothing))
    end
    # A core slower than `SLOW` times the old default is not taken to tighter
    # tolerances.
    cap = SLOW * only(r.wall for r in runs if r.reltol == d.default)
    for reltol in LADDER
        cp = core_problem(d, reltol)
        sol, wall = timed(() -> new_solve(d, cp); cap)
        push!(runs, (; core = "new", reltol, sol, wall, ref = new_ref, cp))
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
    old_default = only(filter(r -> r.core == "old" && r.reltol == d.default, runs))
    target = window_errors(old_default.sol, old_ref, t1, f).rms
    rows = Any[]
    for r in runs
        m = _measures(d, r.sol)
        w = window_errors(r.sol, r.ref, t1, f)
        shift, interval = event_shift(r.sol, r.ref)
        push!(
            rows,
            (;
                deck = d.name,
                cells = d.cells,
                core = r.core,
                reltol = r.reltol,
                retcode = string(r.sol.retcode),
                t_last = r.sol.t[end],
                window_rms = w.rms,
                window_largest = w.largest,
                run_error = run_error(r.sol, r.ref, f),
                event_shift = shift / interval,
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
    old_row = only(filter(r -> r.core == "old" && r.reltol == d.default, rows))
    matched_run =
        matched === nothing ? nothing :
        only(filter(r -> r.core == "new" && r.reltol == matched.reltol, runs))
    anatomy_rows = Any[]
    if (d.name, d.cells) == ANATOMY && matched !== nothing
        for (core, reltol) in (("old", d.default), ("new", matched.reltol))
            run =
                core == "old" ? () -> old_solve(d, reltol) :
                let cp = core_problem(d, reltol)

                    () -> new_solve(d, cp)
                end
            shares = anatomy(run)
            short = d.until(UNIT_TIME)
            unit_sol =
                core == "old" ? old_solve(d, reltol; problem = short) :
                new_solve(d, core_problem(d, reltol; problem = short))
            push!(anatomy_rows, (; core, reltol, shares, unit = unit_costs(unit_sol)))
            println("    anatomy $(core) $(reltol) done")
            flush(stdout)
        end
    end
    t_a = matched_run === nothing ? NaN : first_activation(matched_run.sol)
    split = isfinite(t_a) ? cost_split(d, matched_run, old_default, t_a) : nothing
    # The rule's solve at the default tolerance, for the installed table.
    k = findfirst(r -> r.core == "new" && r.reltol == d.default, runs)
    rule =
        k === nothing ? nothing :
        let sol = runs[k].sol, c = filter(isfinite, sol.correction_velocity)
            (;
                runs[k].wall,
                steps = sol.statistics.accepted,
                sol.outcome,
                t_last = sol.t[end],
                sol.statistics.invocations,
                largest = isempty(c) ? 0.0 : maximum(abs, c),
            )
        end
    return (;
        deck = d.name,
        cells = d.cells,
        d.default,
        t1,
        spread,
        target,
        rows,
        old_row,
        matched,
        anatomy = anatomy_rows,
        t_a,
        split,
        rule,
    )
end

# The same deck with `NoCorrection()` at the default tolerance.
function without_correction(d)
    cp = core_problem(d, d.default; correction = Cl.NoCorrection())
    sol, wall = timed(() -> new_solve(d, cp))
    return (; wall, steps = sol.statistics.accepted, sol.outcome, t_last = sol.t[end])
end

# R41's diagnostic: `ale_darcy` at 100 cells with the layer time scaled by
# (100/25)², so its per-node gain `N²/τ` equals `ale_darcy` 25's at the
# default.
function equal_gain_row(decks)
    d = only(filter(x -> x.name == "ale_darcy" && x.cells == 100, decks))
    τ = (100 / 25)^2 * Cl.EQUIDISTRIBUTION_TIME
    cp = core_problem(d, d.default; correction = Cl.EquidistributionCorrection(; τ))
    sol, wall = timed(() -> new_solve(d, cp))
    return (;
        deck = d.name,
        cells = d.cells,
        τ,
        wall,
        steps = sol.statistics.accepted,
        sol.outcome,
        t_last = sol.t[end],
        sol.statistics.invocations,
    )
end

load_start = Sys.loadavg()
started = time()
results = [
    deck_rows(vbo_deck(deck, n)) for deck in ("V-BO-01", "V-BO-03", "V-BO-04") for
    n in CELLS
]
decks = installed_decks()
installed = map(decks) do d
    core_problem(d, d.default).def.closure.correction isa Cl.LayerCorrection ||
        error("the rule installs no correction on $(d.name) $(d.cells)")
    return (; deck_rows(d)..., none = without_correction(d))
end
equal_gain = equal_gain_row(decks)
load_end = Sys.loadavg()

open(joinpath(OUT, "burnout_timing.csv"), "w") do io
    cols = keys(first(results[1].rows))
    println(io, join(cols, ","))
    for res in vcat(results, installed), r in res.rows
        println(
            io,
            join((v isa AbstractString ? "\"$v\"" : string(v) for v in values(r)), ","),
        )
    end
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "burnout_timing.md"), "w") do io
    println(io, "# Gate item 8: the burnout path against the old core\n")
    println(
        io,
        "Produced by `test/newcore/gates/burnout_timing.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "The finite-volume core's default integrator (`Rodas5P`, heuristic restart); every deck ",
        "with the correction the installation rule gives it. Window RMS and largest: the ",
        "RMS and the largest of the differences over the saved times before the first event of any ",
        "solve of the deck (the whole run where none has one), against the same core at reltol ",
        "$(REFERENCE): of the thickness relative to the reference's on V-BO-04, of the temperatures ",
        "over the cells [K] on the other decks; the new core is matched on the RMS. Run error: the largest difference over the ",
        "whole run at the saved times with the same live cells; shift: the largest invocation-time ",
        "difference over the reference's median interval between invocations (between its saved ",
        "times for a single invocation); MLR: the L1 error of the mass loss rate ",
        "history over the reference's integral. Warm medians of $(REPEATS) solves, of ",
        "$(OLD_DEFAULT_REPEATS) for the old default (a solve slower than $(Int(SLOW)) times the old ",
        "default is timed once and its core not taken to tighter tolerances); ratio: over the old ",
        "core's time at the case's tolerance; item 8: pass at a matched ratio of 1.5 or below, the ",
        "owner's call above 1.5 and up to 2, stop above 2.\n",
    )
    println(io, "## Matched accuracy\n")
    println(
        io,
        "| deck | cells | window before [s] | spread old, new | target | old default [s] | new matched reltol | window RMS | window largest | new [s] | ratio | item 8 | equal-tolerance ratio |",
    )
    println(io, "|" * repeat("---|", 13))
    for res in results
        m = res.matched
        k = findfirst(r -> r.core == "new" && r.reltol == res.default, res.rows)
        eq = k === nothing ? NaN : res.rows[k].ratio
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(_g(res.t1; d = 4)) | $(_g(res.spread.old)), ",
            "$(_g(res.spread.new)) | $(_g(res.target)) | $(_g(res.old_row.wall)) | ",
            m === nothing ? "none on the ladder | | | | | unmatched | " :
            "$(m.reltol) | $(_g(m.window_rms)) | $(_g(m.window_largest)) | $(_g(m.wall)) | " *
            "$(_g(m.ratio)) | $(item8_verdict(m.ratio)) | ",
            "$(_g(eq)) |",
        )
    end
    println(io, "\n## The ladders\n")
    println(io, "The burnout decks, then the decks of the installed correction.\n")
    println(
        io,
        "| deck | cells | core | reltol | retcode | window RMS | window largest | run error | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |",
    )
    println(io, "|" * repeat("---|", 20))
    for res in vcat(results, installed), r in res.rows
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
    println(io, "\n## The installed correction\n")
    println(
        io,
        "The decks on which the installation rule installs the correction, measured as the ",
        "burnout decks are (window before [s]: the first event of any solve of the deck, `Inf` ",
        "where none has one; target: the old default's window RMS [K]). Beside the matched ",
        "ratio and its verdict: the equal-tolerance ratio; the rule's solve at the default ",
        "tolerance and the same deck with `NoCorrection()`, each with its warm time [s], ",
        "accepted steps, outcome and end time [s] (and the rule's invocations); largest |c|: the ",
        "rule's largest correction velocity over the saved times [m/s].\n",
    )
    println(
        io,
        "| deck | cells | default reltol | window before [s] | target | old default [s] (steps; end) | new matched reltol | window RMS | new [s] | ratio | item 8 | equal-tolerance ratio | rule [s] (steps; outcome, end; invocations) | none [s] (steps; outcome, end) | largest |c| |",
    )
    println(io, "|" * repeat("---|", 15))
    for res in installed
        m, o, a, b = res.matched, res.old_row, res.rule, res.none
        k = findfirst(r -> r.core == "new" && r.reltol == res.default, res.rows)
        eq = k === nothing ? NaN : res.rows[k].ratio
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(res.default) | $(_g(res.t1; d = 4)) | $(_g(res.target)) | ",
            "$(_g(o.wall)) ($(o.steps); $(_g(o.t_last; d = 4))) | ",
            m === nothing ? "none on the ladder | | | | unmatched | " :
            "$(m.reltol) | $(_g(m.window_rms)) | $(_g(m.wall)) | $(_g(m.ratio)) | " *
            "$(item8_verdict(m.ratio)) | ",
            "$(_g(eq)) | ",
            a === nothing ? " | " :
            "$(_g(a.wall)) ($(a.steps); $(a.outcome), $(_g(a.t_last; d = 4)); $(a.invocations)) | ",
            "$(_g(b.wall)) ($(b.steps); $(b.outcome), $(_g(b.t_last; d = 4))) | ",
            a === nothing ? " |" : "$(_g(a.largest)) |",
        )
    end
    println(io, "\n## The cost split at the first activation\n")
    println(
        io,
        "Per matched row whose correction activates: t_a, the first saved time of the new core's ",
        "matched solve with a nonzero correction velocity [s]; each core's accepted steps before ",
        "and after it and its share of the solve's time before it, from a read-only step ",
        "recorder in a separate solve (the second of two), applied to the timed medians above; ",
        "the ratio of the new core's time to the old default's before and after t_a; gain: ",
        "`(T_new,after − r_before · T_old,after) / T_new`, the share of the new core's time ",
        "beyond the old core's after t_a at the ratio before it. Before t_a the correction is ",
        "exactly zero, so the ratio there is the integrator's and the per-step cost's (R40); ",
        "after it the gain `N²/τ` acts (R41).\n",
    )
    println(
        io,
        "| deck | cells | new matched reltol | t_a [s] | old steps before, after | new steps before, after | old share before | new share before | ratio before | ratio after | matched ratio | gain |",
    )
    println(io, "|" * repeat("---|", 12))
    for res in vcat(results, installed)
        x = res.split
        x === nothing && continue
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(res.matched.reltol) | $(_g(x.t_a; d = 4)) | ",
            "$(x.old_steps[1]), $(x.old_steps[2]) | $(x.new_steps[1]), $(x.new_steps[2]) | ",
            "$(_g(x.old_share)) | $(_g(x.new_share)) | $(_g(x.r_before)) | $(_g(x.r_after)) | ",
            "$(_g(res.matched.ratio)) | $(_g(x.gain)) |",
        )
    end
    unsplit = [
        "$(res.deck) $(res.cells)" for
        res in vcat(results, installed) if res.split === nothing
    ]
    isempty(unsplit) || println(
        io,
        "\nWithout a split (no matched rung, or no activation in the matched solve): ",
        join(unsplit, ", "),
        ".",
    )
    e = equal_gain
    println(io, "\n## Equal gain (R41, reported)\n")
    println(
        io,
        "`$(e.deck)` at $(e.cells) cells with the layer time τ = $(_g(e.τ; d = 4)) s, whose ",
        "per-node gain `N²/τ` equals `$(e.deck)` 25's at the default, at the default tolerance: ",
        "$(_g(e.wall)) s, $(e.steps) steps, $(e.outcome) at $(_g(e.t_last; d = 4)) s, ",
        "$(e.invocations) invocations. Whether R41's steps are the `N²` itself reads against the ",
        "rule's rows of `$(e.deck)` 25 and 100 above.",
    )
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
