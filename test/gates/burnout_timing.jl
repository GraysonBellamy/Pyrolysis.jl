# The wall time and the accuracy of the burnout path on a tolerance ladder,
# and the cost of the correction where the installation rule installs it.
#
# The burnout decks: V-BO-01, V-BO-03 and V-BO-04 at 30, 60 and 120 cells
# (PMMA cells on the first two), each solved as the suite solves it, with the
# correction `setup`'s installation rule gives it: the equidistribution
# correction at its defaults, the PMMA–Kaowool interface Lagrangian, and on
# V-BO-04's uniform column a correction that stays exactly zero. Each deck
# runs the default integrator at reltol 1e-5, 3e-6, 1e-6, 3e-7 and 1e-7
# (abstol a hundredth of reltol throughout).
#
# Accuracy. A cell deactivates at the end of the step on which it crosses the
# threshold, so an event's time is set by that step, which the tolerance hardly
# moves, and the solution after an event follows its time: over the whole run
# the achieved error does not fall with the tolerance. The error is therefore
# measured over the saved times before the first event of any solve of the
# deck (the window), against the same deck at reltol 1e-9, as the RMS of the
# differences: of the temperatures over the cells on V-BO-01 and V-BO-03 [K],
# of the thickness relative to the reference's on V-BO-04, which stays at its
# initial temperature. The largest difference is reported beside it. The
# reference spread is that reference against one at reltol 1e-10. Beside
# them: the error over the whole run at the saved times with the same live
# cells, the largest invocation-time shift over the reference's median
# interval between invocations, or between its saved times for a single
# invocation (both as `case_accuracy.jl` measures them; the shift is
# `event_shift.jl`'s, shared with it), the L1 error of the mass loss rate
# history relative to the reference's integral, the case's peak mass loss
# rate and burnout time (its mass ledger is `burnout_ladder.jl`'s: the
# case's recorder would add to the timed solves' cost), and the warm wall
# time (median of three solves).
#
# Where the time goes, on V-BO-01 at 60 cells: the sampling profiler over two
# warm solves at the case's tolerance, each sample given to the first of: the
# Jacobian, the domain test, the deactivation events, the factorization, the
# linear solves, the right-hand side, the formation of the iteration matrix,
# the per-step observers; what is left is the integrator's own work. Unit
# costs at the final state of a solve that ends at 200 s: one right-hand side
# and one Jacobian, with the state count and the Jacobian's stored entries.
#
# The installed correction: the decks on which the installation rule installs
# it, the moving-mesh audit decks with deactivation on (`ale` with the solid's
# gas-transfer coefficient at 1e-6 m²/s to 400 s and `ale_darcy` to 300 s, at
# 25 and 100 cells, default reltol 1e-6) and V-CPL-02 and V-CPL-02b at their
# CI resolutions and tolerances, measured on the same ladder as the burnout
# decks (the window before the first event of any solve of the deck, the
# whole run where no solve has one). Beside it, the rule's solve at the
# default tolerance against the same deck with `NoCorrection()` (warm medians
# of three, the steps, the outcome and end time: without the correction a
# layer can end `:interior`) and the largest correction velocity over the
# saved times. A deck on which the rule installs nothing is bit-identical
# (`burnout_identity.jl`) and costs nothing.
#
# The cost split at the first activation, per deck whose correction
# activates at its default tolerance (the installed decks, V-BO-01 and
# V-BO-03): before `t_a`, the first saved time of that solve with a nonzero
# correction velocity, the correction is exactly zero and the solve is the
# problem without it on the wider pattern; after it the gain `N²/τ` acts. The
# accepted steps before and after `t_a` and the share of the time before it
# come from a read-only step recorder in a separate solve, the share applied
# to the timed median so the recorder's own cost stays out. One further row:
# `ale_darcy` at 100 cells with the layer time scaled by 16, whose per-node
# gain equals `ale_darcy` 25's at the default.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/gates/burnout_timing.jl
# Writes test/gates/measurements/burnout_timing.{csv,md}.

ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))
include(joinpath(@__DIR__, "event_shift.jl"))

using OrdinaryDiffEq: DiscreteCallback
using Profile
using SparseArrays: nnz
using Statistics: median

const Ig = Pyrolysis.Integration
const Cl = Pyrolysis.Closures

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REPEATS = 3
const LADDER = (1e-5, 3e-6, 1e-6, 3e-7, 1e-7)
const REFERENCE = 1e-9
const SPREAD = 1e-10
const UNIT_TIME = 200.0
const CELLS = (30, 60, 120)
const ANATOMY = ("V-BO-01", 60)
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

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

# A moving-mesh audit deck: a 10 mm column under 50 kW/m² with re-radiation
# and convection at the top and an adiabatic back, `n` cells, `t_end` [s]
# with `outputs` saved times. PMMA (`ale`) absorbs at its surface; the wood
# (`ale_darcy`) absorbs in depth and carries lumped and Darcy gas transport with a
# pressure condition at the top. Returns the problem and its solve keywords.
function problem_for(kind, n; t_end = 80.0, outputs = 3, solid_λ = nothing, flux = 50e3)
    wood = kind == :ale_darcy
    mat = if wood
        Material(
            name = :audit_wood,
            gas_transport = LumpedTransfer(
                endpoint = LumpedEndpoint(
                    λ = something(solid_λ, 1e-5),
                    κ = 1e-14,
                    max_skeletal_occupancy = 0.001,
                ),
            ),
            components = (
                SolidComponent(
                    :wood,
                    ρ = 500.0,
                    ρ_intrinsic = 1500.0,
                    c = 1500.0,
                    k = 0.15,
                    ε = 0.9,
                    α = 5.0,
                    κ = 1e-14,
                    λ = something(solid_λ, 1e-5),
                ),
                GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03),
            ),
            reactions = (Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = 300e3),),
        )
    else
        Material(
            name = :audit_pmma,
            gas_transport = LumpedTransfer(
                endpoint = LumpedEndpoint(
                    λ = something(solid_λ, 0.0),
                    κ = 0.0,
                    max_skeletal_occupancy = 0.001,
                ),
            ),
            components = (
                SolidComponent(
                    :pmma,
                    ρ = 1190.0,
                    φ = 0.01,
                    c = 1420.0,
                    k = 0.21,
                    ε = 0.86,
                    λ = something(solid_λ, 0.0),
                ),
                GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02),
            ),
            reactions = (Reaction(:pyrolysis, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
        )
    end
    thermal = (
        top = RadiativeFluxBC(flux = flux) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    )
    mass =
        wood ?
        (
            top = AmbientGasBC(
                GasReservoir(P = 101325.0, T = 300.0, mole_fractions = (gas = 1.0,)),
            ),
            bottom = ImpermeableBC(),
        ) : (top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)), bottom = ImpermeableBC())
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n),
        material = mat,
        bc_set = BoundaryConditionSet(; thermal, mass),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = outputs)),
    )
    return problem,
    (; use_ale = true, radiation_model = wood ? BEER_LAMBERT : SURFACE_ABSORPTION)
end

# A deck of the tables: its name and cells, the problem, `setup`'s keywords,
# the default relative tolerance, the window's field (`:T`, the cells'
# temperatures, or `:thickness`, the thickness relative to the reference's),
# `min_thickness`, and the problem ending at a given time for the unit costs
# (`nothing` where none are taken).
function vbo_deck(deck, n)
    kw = CASES[deck].solve_kwargs
    rad = haskey(kw, :radiation_model) ? (; radiation_model = kw.radiation_model) : (;)
    return (;
        name = deck,
        cells = n,
        problem = deck_problem(deck, n),
        keywords = (; use_ale = kw.use_ale, deactivation = kw.deactivation, rad...),
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
        push!(
            decks,
            (;
                name = string(kind),
                cells = n,
                problem,
                keywords = (; kw..., deactivation = true),
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
                keywords = (;
                    use_ale = kw.use_ale,
                    deactivation = get(kw, :deactivation, false),
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

# The set-up problem of a deck, as the suite solves it.
core_problem(d, reltol; problem = d.problem, kw...) =
    Ig.setup(problem; d.keywords..., _tol(reltol)..., kw...)

deck_solve(d, cp; kw...) = Ig.integrate(cp; min_thickness = d.min_thickness, kw...)

# ----------------------------------------------------------------------------
# Error measures
# ----------------------------------------------------------------------------

_event_times(sol) = [e.t for e in sol.events]
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
    "/src/Physics/",
    "/src/BoundaryConditions/",
    "/src/Properties/",
)

function _part(frame)
    f, file = string(frame.func), string(frame.file)
    occursin("/src/ColoredJacobian/", file) && return :jacobian
    occursin("Integration/checks.jl", file) && return :domain
    occursin("/src/Deactivation/", file) && return :events
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
    f in ("_observe!", "cap_step!", "apply_observer_flag!") && return :observers
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

# The warm median of `REPEATS` solves after one warm-up.
function timed(f)
    sol = f()
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed sol = f())
    end
    return sol, median(times)
end

_counts(sol) = (
    steps = sol.statistics.accepted,
    rejected = sol.statistics.rejected,
    rhs = sol.statistics.rhs_evaluations,
    jacobians = sol.statistics.jacobians,
    nonlinear_failures = sol.statistics.nonlinear_failures,
)

_measures(d, sol) =
    d.name in ("V-BO-01", "V-BO-03") ? Base.invokelatest(_vbo_measures, sol) :
    (; peak = NaN, t_bo = NaN)

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

# The cost split of the default tolerance's solve at `t_a`: its accepted steps
# before and after it, the share of its time before it, and its time before
# and after it (the share applied to the timed median `wall`).
function cost_split(d, cp, wall, t_a)
    s = split_at(callback -> deck_solve(d, cp; callback), t_a)
    return (;
        t_a,
        steps = s.steps,
        share = s.share,
        before = s.share * wall,
        after = (1 - s.share) * wall,
    )
end

# ----------------------------------------------------------------------------
# The decks' rows
# ----------------------------------------------------------------------------

function deck_rows(d)
    println("[$(d.name) $(d.cells)] references")
    flush(stdout)
    f = d.field
    ref = deck_solve(d, core_problem(d, REFERENCE))
    spread_solve = deck_solve(d, core_problem(d, SPREAD))
    runs = Any[]
    for reltol in LADDER
        cp = core_problem(d, reltol)
        sol, wall = timed(() -> deck_solve(d, cp))
        push!(runs, (; reltol, sol, wall, cp))
        println("    $(reltol): $(sol.retcode) $(round(wall, digits = 3)) s")
        flush(stdout)
    end
    solves = vcat([r.sol for r in runs], [ref, spread_solve])
    t1 = minimum(s -> minimum(_event_times(s); init = Inf), solves)
    spread = window_errors(spread_solve, ref, t1, f).rms
    rows = Any[]
    for r in runs
        m = _measures(d, r.sol)
        w = window_errors(r.sol, ref, t1, f)
        shift, interval = event_shift(r.sol, ref)
        push!(
            rows,
            (;
                deck = d.name,
                cells = d.cells,
                reltol = r.reltol,
                retcode = string(r.sol.retcode),
                t_last = r.sol.t[end],
                window_rms = w.rms,
                window_largest = w.largest,
                run_error = run_error(r.sol, ref, f),
                event_shift = shift / interval,
                mlr_error = mlr_error(r.sol, ref),
                peak = m.peak,
                t_bo = m.t_bo,
                wall = r.wall,
                _counts(r.sol)...,
            ),
        )
    end
    k = findfirst(r -> r.reltol == d.default, runs)
    default_run = k === nothing ? nothing : runs[k]
    anatomy_rows = Any[]
    if (d.name, d.cells) == ANATOMY
        cp = core_problem(d, d.default)
        shares = anatomy(() -> deck_solve(d, cp))
        short = d.until(UNIT_TIME)
        unit_sol = deck_solve(d, core_problem(d, d.default; problem = short))
        push!(anatomy_rows, (; reltol = d.default, shares, unit = unit_costs(unit_sol)))
        println("    anatomy done")
        flush(stdout)
    end
    t_a = default_run === nothing ? NaN : first_activation(default_run.sol)
    split = isfinite(t_a) ? cost_split(d, default_run.cp, default_run.wall, t_a) : nothing
    # The rule's solve at the default tolerance, for the installed table.
    rule =
        default_run === nothing ? nothing :
        let sol = default_run.sol, c = filter(isfinite, sol.correction_velocity)
            (;
                default_run.wall,
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
        rows,
        anatomy = anatomy_rows,
        t_a,
        split,
        rule,
    )
end

# The same deck with `NoCorrection()` at the default tolerance.
function without_correction(d)
    cp = core_problem(d, d.default; correction = Cl.NoCorrection())
    sol, wall = timed(() -> deck_solve(d, cp))
    return (; wall, steps = sol.statistics.accepted, sol.outcome, t_last = sol.t[end])
end

# `ale_darcy` at 100 cells with the layer time scaled by (100/25)², so its
# per-node gain `N²/τ` equals `ale_darcy` 25's at the default.
function equal_gain_row(decks)
    d = only(filter(x -> x.name == "ale_darcy" && x.cells == 100, decks))
    τ = (100 / 25)^2 * Cl.EQUIDISTRIBUTION_TIME
    cp = core_problem(d, d.default; correction = Cl.EquidistributionCorrection(; τ))
    sol, wall = timed(() -> deck_solve(d, cp))
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
    println(io, "# The burnout path: wall time and accuracy on a tolerance ladder\n")
    println(
        io,
        "Produced by `test/gates/burnout_timing.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "The default integrator (`Rodas5P`, heuristic restart); every deck ",
        "with the correction the installation rule gives it. Window RMS and largest: the ",
        "RMS and the largest of the differences over the saved times before the first event of any ",
        "solve of the deck (the whole run where none has one), against the same deck at reltol ",
        "$(REFERENCE): of the thickness relative to the reference's on V-BO-04, of the temperatures ",
        "over the cells [K] on the other decks. Run error: the largest difference over the ",
        "whole run at the saved times with the same live cells; shift: the largest invocation-time ",
        "difference over the reference's median interval between invocations (between its saved ",
        "times for a single invocation); MLR: the L1 error of the mass loss rate ",
        "history over the reference's integral. Warm medians of $(REPEATS) solves [s].\n",
    )
    println(io, "## The decks\n")
    println(
        io,
        "Per deck: its default tolerance, the window (before the first event of any solve, ",
        "`Inf` where none has one) and the reference spread (the reference against reltol ",
        "$(SPREAD), over the window).\n",
    )
    println(io, "| deck | cells | default reltol | window before [s] | spread |")
    println(io, "|---|---|---|---|---|")
    for res in vcat(results, installed)
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(res.default) | $(_g(res.t1; d = 4)) | $(_g(res.spread)) |",
        )
    end
    println(io, "\n## The ladders\n")
    println(io, "The burnout decks, then the decks of the installed correction.\n")
    println(
        io,
        "| deck | cells | reltol | retcode | window RMS | window largest | run error | shift | MLR | peak | t_bo [s] | warm [s] | steps | rejected | RHS | Jacobians | nonlinear failures |",
    )
    println(io, "|" * repeat("---|", 17))
    for res in vcat(results, installed), r in res.rows
        println(
            io,
            "| $(r.deck) | $(r.cells) | $(r.reltol) | $(r.retcode) | ",
            "$(_g(r.window_rms)) | $(_g(r.window_largest)) | $(_g(r.run_error)) | ",
            "$(_g(r.event_shift)) | $(_g(r.mlr_error)) | $(_g(r.peak; d = 5)) | ",
            "$(_g(r.t_bo; d = 4)) | $(_g(r.wall)) | ",
            "$(r.steps) | $(r.rejected) | $(r.rhs) | $(r.jacobians) | $(r.nonlinear_failures) |",
        )
    end
    for res in results
        isempty(res.anatomy) && continue
        println(
            io,
            "\n## Where the time goes: $(res.deck) at $(res.cells) cells\n\nShares of the ",
            "profiler's samples at the case's tolerance, and unit costs at $(UNIT_TIME) s [µs]:\n",
        )
        println(
            io,
            "| reltol | " *
            join(string.(PARTS), " | ") *
            " | RHS | Jacobian | states | nnz(J) |",
        )
        println(io, "|" * repeat("---|", 1 + length(PARTS) + 4))
        for a in res.anatomy
            u = a.unit
            println(
                io,
                "| $(a.reltol) | " *
                join((_g(a.shares[p]; d = 2) for p in PARTS), " | ") *
                " | $(_g(1e6u.rhs)) | $(_g(1e6u.jac)) | $(u.states) | $(u.nnz) |",
            )
        end
    end
    println(io, "\n## The installed correction\n")
    println(
        io,
        "The decks on which the installation rule installs the correction, at the default ",
        "tolerance: the rule's solve and the same deck with `NoCorrection()`, each with its warm ",
        "time [s], accepted steps, outcome and end time [s] (and the rule's invocations); largest ",
        "|c|: the rule's largest correction velocity over the saved times [m/s].\n",
    )
    println(
        io,
        "| deck | cells | default reltol | rule [s] (steps; outcome, end; invocations) | none [s] (steps; outcome, end) | largest |c| |",
    )
    println(io, "|" * repeat("---|", 6))
    for res in installed
        a, b = res.rule, res.none
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(res.default) | ",
            a === nothing ? " | " :
            "$(_g(a.wall)) ($(a.steps); $(a.outcome), $(_g(a.t_last; d = 4)); $(a.invocations)) | ",
            "$(_g(b.wall)) ($(b.steps); $(b.outcome), $(_g(b.t_last; d = 4))) | ",
            a === nothing ? " |" : "$(_g(a.largest)) |",
        )
    end
    println(io, "\n## The cost split at the first activation\n")
    println(
        io,
        "Per deck whose correction activates at its default tolerance: t_a, the first saved time ",
        "of that solve with a nonzero correction velocity [s]; the accepted steps before and after ",
        "it and the share of the solve's time before it, from a read-only step recorder in a ",
        "separate solve (the second of two), applied to the timed median: the time before and after ",
        "t_a [s]. Before t_a the correction is exactly zero; after it the gain `N²/τ` acts.\n",
    )
    println(
        io,
        "| deck | cells | reltol | t_a [s] | steps before, after | share before | before [s] | after [s] |",
    )
    println(io, "|" * repeat("---|", 8))
    for res in vcat(results, installed)
        x = res.split
        x === nothing && continue
        println(
            io,
            "| $(res.deck) | $(res.cells) | $(res.default) | $(_g(x.t_a; d = 4)) | ",
            "$(x.steps[1]), $(x.steps[2]) | $(_g(x.share)) | $(_g(x.before)) | $(_g(x.after)) |",
        )
    end
    unsplit = [
        "$(res.deck) $(res.cells)" for
        res in vcat(results, installed) if res.split === nothing
    ]
    isempty(unsplit) || println(
        io,
        "\nWithout a split (no activation at the default tolerance): ",
        join(unsplit, ", "),
        ".",
    )
    e = equal_gain
    println(io, "\n## Equal gain\n")
    println(
        io,
        "`$(e.deck)` at $(e.cells) cells with the layer time τ = $(_g(e.τ; d = 4)) s, whose ",
        "per-node gain `N²/τ` equals `$(e.deck)` 25's at the default, at the default tolerance: ",
        "$(_g(e.wall)) s, $(e.steps) steps, $(e.outcome) at $(_g(e.t_last; d = 4)) s, ",
        "$(e.invocations) invocations. Whether the steps follow the gain `N²/τ` reads against the ",
        "rule's rows of `$(e.deck)` 25 and 100 above.",
    )
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
