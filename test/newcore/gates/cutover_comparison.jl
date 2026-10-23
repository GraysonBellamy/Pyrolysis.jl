# The two solvers compared end to end before the old one is deleted: per
# deck, each core's accuracy against its own tight solve, their common
# physical quantities over the whole run against budgets, the wall time at
# matched accuracy with the integrators' counters, and each core's memory in
# a process of its own.
#
# The rows: the six audit decks (`problem_for`'s `eulerian`, `ale` and
# `ale_darcy` at 25 and 100 cells, the solid's gas-transfer coefficient
# 1e-6 m²/s on the two PMMA decks, 80 s, 81 outputs); V-BO-01 to V-BO-08 at
# their CI resolutions with their solve keywords; the Douglas fir deck of
# `examples/07_douglas_fir_cone.jl` over 600 s with 1 s outputs; the
# white-pine forward of `examples/inverse_analysis/` at the stage-2
# conductivities of `plot_stage2_theta.jl` on its two fit scenarios (30 and
# 50 kW/m²); the PMMA cone of `examples/06_cone_calorimeter.jl`. Each deck is
# built once per core: the cases from `casedefs.jl` loaded into two modules,
# one routed to each core, so each core's case QoIs are its own; the other
# decks from one description with each core's mesh constructor (the same
# nodes, cross section and initial-field functions). The old core runs its
# default integrator and linear solver with the deck's keywords, the
# finite-volume core its defaults in the temperature form, `deactivation`
# for `handle_depletion`.
#
# Accuracy. Each core is solved at reltol 1e-9 (the reference) and 1e-10
# (the spread), abstol at the deck's ratio to reltol; the old core at the
# deck's tolerance (its default), the finite-volume core on the ladder 1e-5
# to 1e-8. The window is the saved times before the first event of any solve
# of the deck, a stop before the end of the span counting as one (the whole
# run where none has one); its RMS is of the temperatures over the live cells
# against the same core's reference (of the thickness relative to the
# reference's on V-BO-04). The matched rung is the loosest of the ladder at or
# below the old default's window RMS whose own closure is within its budget
# (the case's ledger tolerance on the cases); the loosest on the window alone
# is reported beside it. Between the old default and the matched solve, over
# the saved times both hold: the temperature at the old core's live cell
# centres, interpolated linearly between the finite-volume core's live cell
# centres (RMS and largest; 1 K and 5 K); the mass loss (0.3 points); the
# burnout time (1 %) on the decks that burn out, the first saved time the
# thickness is within 0.02 mm of the larger of the two final thicknesses; the
# main peak of the gas-generation rate, the largest before the solve's last
# event, and its time (3 % each) where it is interior on both; the tail, the
# gas generated from the main peak to the end (10 %). A deck without events has a verdict on each. On a deck with events
# the temperatures before the first event have one; after the events a
# quantity has a common reference only where the case measures it against its
# own reference value and both cores pass that QoI, and is reported without a
# verdict otherwise. Each core's conservation closure is its own: the
# finite-volume core's from `ConservationRecorder` over the accepted steps
# (the case's per-step ledger on the cases), the old core's the total-mass
# drop against the trapezoid of ṁ″ on the saved times; on the cases it is
# held to the case's tolerance.
#
# Timing, both cores in one process: after a warm-up, the median of nine
# solves for the old default and of three for each of the finite-volume
# core's matched and default rungs (three for any configuration whose
# warm-up exceeds 15 s); the counters of the timed solves. Read only where
# every verdict of the row's accuracy passes: pass at a ratio of 1.5 or
# below, the owner's call to 2. Beside it, the finite-volume core's warm
# time (median and spread of five) and allocated bytes on V-BO-01 at 60
# cells, `ale` at 100 cells and the Douglas fir deck over 60 s, the numbers
# the deletion commit's recheck reads.
#
# Memory: each core on each row in a process of its own at the timing's
# tolerance, one warm-up solve, then the measured solve's allocated bytes and
# the process's peak resident set.
#
# The white-pine objective (the data term of `loss` at the fitted θ) on each
# core: the old core as fitted (`KenCarp4` with Sparspak, `build_and_solve`)
# and at its default, the finite-volume core at its default.
#
# The accuracy and memory parts do not depend on the load and run as
# parallel processes; the timing part runs alone, pinned to one core, at a
# one-minute load below 2. Run from the repository root, one thread each:
#   julia --project=ENV cutover_comparison.jl accuracy K DIR   (K = 1 to 4)
#   julia --project=ENV cutover_comparison.jl memory DIR [JOBS]
#   taskset -c 1 julia --project=ENV cutover_comparison.jl timing DIR
#   julia --project=ENV cutover_comparison.jl merge DIR
# `merge` writes test/newcore/gates/measurements/cutover_comparison.{csv,md}.
# With `--smoke` every mode runs a reduced matrix (the references at reltol
# 1e-7 and 3e-8, two rungs, one timed solve) and `merge` writes into DIR.

using Pyrolysis
using Logging
using Printf
using Serialization
using Statistics: median, mean

include(joinpath(@__DIR__, "stamp.jl"))

const Ig = Pyrolysis.Integration
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUT = joinpath(@__DIR__, "measurements")
const PARTS = 4
const BUDGET = (
    T_rms = 1.0,
    T_max = 5.0,
    mass_loss = 0.3,
    burnout = 0.01,
    peak = 0.03,
    t_peak = 0.03,
    tail = 0.10,
)
const BURNOUT_MARGIN = 0.02e-3
const SLOW_WARMUP = 15.0
const RECHECK_REPEATS = 5

# The full matrix and the smoke run's reduced one.
const FULL = (;
    reference = 1e-9,
    spread = 1e-10,
    ladder = (1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8),
    old_repeats = 9,
    new_repeats = 3,
    recheck = RECHECK_REPEATS,
)
const SMOKE = (;
    reference = 1e-7,
    spread = 3e-8,
    ladder = (1e-5, 1e-6),
    old_repeats = 1,
    new_repeats = 1,
    recheck = 1,
)
const CONFIG = Ref{Any}(FULL)

quiet(f) = Logging.with_logger(f, Logging.NullLogger())
latest(m, name, args...; kw...) =
    Base.invokelatest(Base.invokelatest(getglobal, m, name), args...; kw...)

# ----------------------------------------------------------------------------
# The decks
# ----------------------------------------------------------------------------

# `casedefs.jl` loaded into a module routed to `core` ("old" or "new").
const CASE_MODULES = Dict{String,Any}()
function cases_on(core)
    return get!(CASE_MODULES, core) do
        ENV["PYROLYSIS_CORE"] = core
        ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
        m = Module(Symbol("cases_", core))
        Core.eval(m, :(include(p) = Base.include($m, p)))
        Base.include(m, joinpath(REPO, "test", "verification", "casedefs.jl"))
        list = latest(m, :load_verification_cases)
        (; m, cases = Dict(c.id => c for c in list))
    end
end

# The audit decks' `problem_for`, loaded into a module of its own.
const AUDIT = Ref{Any}(nothing)
function audit_module()
    AUDIT[] === nothing || return AUDIT[]
    m = Module(:audit_decks)
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(
        m,
        joinpath(REPO, "scripts", "mesh_phase0", "baselines", "jacobian_common.jl"),
    )
    return AUDIT[] = m
end

# `problem` on the finite-volume core's mesh `mesh`.
native(problem, mesh) = PyrolysisProblem(
    mesh = mesh,
    material = problem.material,
    bc_set = problem.bc_set,
    T_initial = problem.T_initial,
    ξ_initial = problem.ξ_initial,
    tspan = problem.tspan,
    dt_initial = problem.dt_initial,
    output_times = problem.output_times,
)

# `problem` with its time span and saved times replaced.
function retimed(problem, t_end, outputs)
    return PyrolysisProblem(
        mesh = problem.mesh,
        material = problem.material,
        bc_set = problem.bc_set,
        T_initial = problem.T_initial,
        ξ_initial = problem.ξ_initial,
        tspan = (0.0, t_end),
        dt_initial = problem.dt_initial,
        output_times = collect(outputs),
    )
end

# An example script's statements from the line starting `first` to the line
# before the one starting `stop`, evaluated in a module of their own without
# their printing.
function example_module(name, file, first, stop)
    lines = readlines(joinpath(REPO, "examples", file))
    a = findfirst(startswith(first), lines)
    b = findfirst(startswith(stop), lines)
    (a === nothing || b === nothing) && error("$(file): no `$(first)` … `$(stop)` block")
    out = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") && continue
        push!(out, line)
    end
    m = Module(name)
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, Meta.parseall(join(out, "\n")))
    return m
end

const EXAMPLE_MODULES = Dict{Symbol,Any}()
douglas_fir() = get!(
    () -> example_module(
        :douglas_fir_deck,
        "07_douglas_fir_cone.jl",
        "const T_HFM_MAX",
        "solution = solve(",
    ),
    EXAMPLE_MODULES,
    :douglas_fir,
)
pmma_cone() = get!(
    () -> example_module(
        :pmma_cone_deck,
        "06_cone_calorimeter.jl",
        "pmma_cp = ",
        "solution = solve(",
    ),
    EXAMPLE_MODULES,
    :pmma_cone,
)

# The white-pine forward model and its data, and the fitted stage-2
# conductivities of `plot_stage2_theta.jl`.
function white_pine()
    return get!(EXAMPLE_MODULES, :white_pine) do
        m = Module(:white_pine)
        Core.eval(m, :(include(p) = Base.include($m, p)))
        Base.include(m, joinpath(REPO, "examples", "inverse_analysis", "wp_objective.jl"))
        file = joinpath(REPO, "examples", "inverse_analysis", "plot_stage2_theta.jl")
        theta = nothing
        for ex in Meta.parseall(read(file, String)).args
            ex isa Expr && ex.head === :const || continue
            a = ex.args[1]
            a isa Expr && a.head === :(=) && a.args[1] === :THETA_DEFAULT || continue
            theta = Core.eval(m, a.args[2])
        end
        theta === nothing && error("no THETA_DEFAULT in $(file)")
        dataset = Dict(
            l => quiet(() -> latest(m, :load_registry_dataset, l; verbose = false)) for
            l in ("30kW", "50kW")
        )
        (; m, theta = Float64.(theta), dataset)
    end
end
const WP_STAGE = 2

# The old core's keywords that the finite-volume core does not take, and the
# one it names differently.
const OLD_ONLY = (:depletion_threshold, :min_cells, :progress, :integrator, :jacobian)
function new_keywords(kw)
    nt = (; (k => v for (k, v) in pairs(kw) if !(k in OLD_ONLY))...)
    haskey(nt, :handle_depletion) || return nt
    return (;
        (k => v for (k, v) in pairs(nt) if k != :handle_depletion)...,
        deactivation = nt.handle_depletion,
    )
end
without_tolerances(kw) =
    (; (k => v for (k, v) in pairs(kw) if !(k in (:reltol, :abstol)))...)

# A row: its name; `build(core)`, the deck's problem on `core`; `keywords`,
# the old core's solve keywords without the tolerances; the deck's reltol and
# abstol; the window's field; whether the deck burns out; the case id or
# `nothing`; the white-pine scenario or `nothing`.
function row(
    name;
    build,
    keywords,
    reltol,
    abstol,
    field = :T,
    burnout = false,
    case = nothing,
    scenario = nothing,
)
    return (;
        name,
        build,
        keywords = without_tolerances(keywords),
        reltol,
        ratio = abstol / reltol,
        field,
        burnout,
        case,
        scenario,
    )
end

const AUDIT_DECKS = (
    (:eulerian, 25, 1e-6),
    (:eulerian, 100, 1e-6),
    (:ale, 25, 1e-6),
    (:ale, 100, 1e-6),
    (:ale_darcy, 25, nothing),
    (:ale_darcy, 100, nothing),
)
const CASE_IDS =
    ("V-BO-01", "V-BO-02", "V-BO-03", "V-BO-04", "V-BO-05", "V-BO-06", "V-BO-07", "V-BO-08")
const BURNS_OUT = ("V-BO-01", "V-BO-02", "V-BO-03", "V-BO-06", "V-BO-07", "V-BO-08")

function audit_row(kind, n, λ)
    old() =
        Base.invokelatest(audit_module().problem_for, kind, n; solid_λ = λ, outputs = 81)
    build(core) = begin
        problem, _ = old()
        core == "old" ? problem : native(problem, create_uniform_mesh(0.01, n))
    end
    _, kw = old()
    return row(
        audit_name(kind, n, λ);
        build,
        keywords = (; use_ale = kw.use_ale, radiation_model = kw.radiation_model),
        reltol = 1e-6,
        abstol = 1e-8,
    )
end

function case_row(id; core = "old")
    kw = cases_on(core).cases[id].solve_kwargs
    build(core) = begin
        c = cases_on(core).cases[id]
        Base.invokelatest(c.build, c.ci.n_cells)
    end
    return row(
        id;
        build,
        keywords = kw,
        reltol = kw.reltol,
        abstol = kw.abstol,
        field = id == "V-BO-04" ? :thickness : :T,
        burnout = id in BURNS_OUT,
        case = id,
    )
end

function douglas_fir_row(t_end)
    build(core) = begin
        m = douglas_fir()
        problem = Base.invokelatest(getglobal, m, :problem)
        mesh =
            core == "old" ? problem.mesh :
            create_uniform_mesh(
                Base.invokelatest(getglobal, m, :sample_thickness),
                Base.invokelatest(getglobal, m, :n_cells),
            )
        native(retimed(problem, t_end, 0.0:1.0:t_end), mesh)
    end
    return row(
        "Douglas fir, $(round(Int, t_end)) s";
        build,
        keywords = (;
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            dtmin = 1e-14,
            min_thickness = 0.5e-3,
            handle_depletion = false,
        ),
        reltol = 1e-7,
        abstol = 1e-9,
    )
end

function white_pine_row(label)
    build(core) = begin
        wp = white_pine()
        scenario = first(wp.dataset[label])
        problem = latest(wp.m, :build_problem, wp.theta, scenario; stage = WP_STAGE)
        core == "old" ? problem :
        native(
            problem,
            create_uniform_mesh(
                scenario.sample_thickness,
                scenario.n_cells;
                cross_section = get(
                    scenario,
                    :A_0,
                    Base.invokelatest(getglobal, wp.m, :WP_A0),
                ),
            ),
        )
    end
    kw = latest(white_pine().m, :solve_keywords)
    return row(
        "white pine, $(label)";
        build,
        keywords = (;
            (k => v for (k, v) in pairs(kw) if !(k in (:integrator, :jacobian)))...
        ),
        reltol = kw.reltol,
        abstol = kw.abstol,
        scenario = label,
    )
end

function pmma_cone_row()
    build(core) = begin
        m = pmma_cone()
        problem = Base.invokelatest(getglobal, m, :problem)
        core == "old" ? problem :
        native(
            problem,
            create_uniform_mesh(
                Base.invokelatest(getglobal, m, :sample_thickness),
                Base.invokelatest(getglobal, m, :n_cells),
            ),
        )
    end
    return row(
        "PMMA cone";
        build,
        keywords = (;
            use_ale = true,
            dtmin = 1e-14,
            min_thickness = 0.01e-3,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 1,
        ),
        reltol = 1e-6,
        abstol = 1e-8,
        burnout = true,
    )
end

# The rows by name, built when asked for: a process measuring one row loads
# that row's deck only.
audit_name(kind, n, λ) = λ === nothing ? "$(kind) $(n)" : "$(kind) $(n), solid λ $(λ)"
const ROWS = vcat(
    [audit_name(k...) => (core -> audit_row(k...)) for k in AUDIT_DECKS],
    [id => (core -> case_row(id; core)) for id in CASE_IDS],
    [
        "Douglas fir, 600 s" => (core -> douglas_fir_row(600.0)),
        "white pine, 30kW" => (core -> white_pine_row("30kW")),
        "white pine, 50kW" => (core -> white_pine_row("50kW")),
        "PMMA cone" => (core -> pmma_cone_row()),
    ],
)
const ROW_NAMES = first.(ROWS)
row_named(name; core = "old") = Dict(ROWS)[name](core)

# The parts of the accuracy measurement, by row name, balanced on the
# estimated cost (the white pine's two scenarios the longest).
part_of(name) =
    startswith(name, "white pine, 30") ? 1 :
    startswith(name, "white pine, 50") ? 2 : startswith(name, "V-BO") ? 4 : 3

# ----------------------------------------------------------------------------
# The solves
# ----------------------------------------------------------------------------

tolerances(r, reltol) = (; reltol, abstol = r.ratio * reltol)

# A plain solve of row `r` on `core` (timing, memory): the old core's
# `solve` with the deck's keywords, the finite-volume core's with its
# defaults; `callback` passed to either. Through `invokelatest`, since the
# decks' initial-field functions are defined after this script's world.
function plain_solve(r, core, problem, reltol; callback = nothing)
    tol = tolerances(r, reltol)
    cb = callback === nothing ? (;) : (; callback)
    return quiet() do
        core == "old" ?
        Base.invokelatest(solve, problem; progress = false, r.keywords..., tol..., cb...) :
        Base.invokelatest(
            Ig.solve,
            problem;
            energy_state = :temperature,
            new_keywords(r.keywords)...,
            tol...,
            cb...,
        )
    end
end

build(r, core) = Base.invokelatest(r.build, core)

# The accuracy part's solve: on a case through its own solve function, whose
# step records give the QoIs; on the finite-volume core otherwise with a
# conservation recorder. Returns the solution, its QoIs (`nothing` off the
# cases) and its closure.
function accuracy_solve(r, core, problem, reltol)
    tol = tolerances(r, reltol)
    if r.case !== nothing
        cm = cases_on(core)
        case = cm.cases[r.case]
        sol = quiet(
            () ->
                latest(cm.m, :_case_solve, case, problem; case.solve_kwargs..., tol...),
        )
        qois = Base.invokelatest(case.qois, sol)
        ledger = [q for q in qois if startswith(q.name, "mass ledger")]
        closure =
            isempty(ledger) ? old_closure(sol) :
            (; mass = only(ledger).value, energy = NaN, budget = only(ledger).tol)
        return sol, qois, closure
    end
    if core == "new"
        rec = Pyrolysis.ConservationRecorder()
        sol = plain_solve(r, core, problem, reltol; callback = rec)
        audit = Pyrolysis.conservation_audit(sol, rec)
        closure =
            (; mass = audit.mass.relative, energy = audit.energy.relative, budget = NaN)
        return sol, nothing, closure
    end
    sol = plain_solve(r, core, problem, reltol)
    return sol, nothing, old_closure(sol)
end

# The old core's closure: the total-mass drop against the trapezoid of ṁ″
# on the saved times, as the cases measure it, both per unit initial area
# (`MLR_gauge`; `mass_loss_rate` is per current area). A closure has a
# budget only on a case, its ledger QoI's tolerance.
function old_closure(sol)
    m, ṁ, t = sol.extras.total_mass, sol.extras.MLR_gauge, sol.t
    length(t) < 2 && return (; mass = NaN, energy = NaN, budget = NaN)
    gen = sum(0.5 * (ṁ[j] + ṁ[j+1]) * (t[j+1] - t[j]) for j = 1:(length(t)-1))
    drop = m[1] - m[end]
    return (; mass = drop == 0 ? NaN : abs(drop - gen) / drop, energy = NaN, budget = NaN)
end

# Whether a solve's own closure is within its budget (none off the cases).
closes(closure) = isnan(closure.budget) || abs(closure.mass) < closure.budget

is_new(sol) = sol isa Pyrolysis.Outputs.CoreSolution
event_times(sol) = is_new(sol) ? [e.t for e in sol.events] : collect(sol.extras.park_times)
outcome(sol) = is_new(sol) ? sol.outcome : :none
finished(sol) = sol.retcode in (:Success, :Terminated)

function counters(sol)
    if is_new(sol)
        s = sol.statistics
        return (;
            steps = s.accepted,
            rejected = s.rejected,
            rhs = s.rhs_evaluations,
            jacobians = s.jacobians,
            linear_solves = s.linear_solves,
            nonlinear_failures = s.nonlinear_failures,
            domain_rejections = s.domain_rejections,
            events = length(sol.events),
            invocations = s.invocations,
        )
    end
    s = sol.raw_solution.stats
    parks = sol.extras.park_times
    return (;
        steps = s.naccept,
        rejected = s.nreject,
        rhs = s.nf,
        jacobians = s.njacs,
        linear_solves = s.nsolve,
        nonlinear_failures = s.nnonlinconvfail,
        domain_rejections = -1,
        events = length(parks),
        invocations = length(unique(parks)),
    )
end

# ----------------------------------------------------------------------------
# The measures
# ----------------------------------------------------------------------------

# The differences at the saved times before `t1` against the reference:
# the live cells' temperatures, or the thickness relative to the reference's.
function window_errors(sol, ref, t1, field)
    index = Dict(t => k for (k, t) in enumerate(ref.t))
    d = Float64[]
    for (j, t) in enumerate(sol.t)
        t < t1 || continue
        k = get(index, t, 0)
        k == 0 && continue
        if field === :thickness
            a, b = sol.extras.thickness_history[j], ref.extras.thickness_history[k]
            push!(d, (a - b) / b)
        else
            for i in axes(sol.T, 1)
                a, b = sol.T[i, j], ref.T[i, k]
                (isfinite(a) && isfinite(b)) && push!(d, a - b)
            end
        end
    end
    return (;
        rms = isempty(d) ? NaN : sqrt(sum(abs2, d) / length(d)),
        largest = isempty(d) ? NaN : maximum(abs, d),
    )
end

# The live cells' centres [m] and temperatures at saved time `j`.
function live_profile(sol, j)
    T = view(sol.T, :, j)
    if is_new(sol)
        c = view(sol.cell_centers, :, j)
    elseif size(sol.z, 1) == length(T) + 1 && all(isfinite, view(sol.z, 1:2, j))
        z = view(sol.z, :, j)
        c = [0.5 * (z[i] + z[i+1]) for i in eachindex(T)]
    else
        c = sol.problem.mesh.cell_centers
    end
    keep = [i for i in eachindex(T) if isfinite(T[i]) && isfinite(c[i])]
    p = sortperm(c[keep])
    return c[keep][p], T[keep][p]
end

function interpolate(x, y, at)
    k = searchsortedlast(x, at)
    k == 0 && return y[1]
    k == length(x) && return y[end]
    w = (at - x[k]) / (x[k+1] - x[k])
    return (1 - w) * y[k] + w * y[k+1]
end

# The temperature difference at the old core's live cell centres inside the
# finite-volume core's live region, at the saved times before `t1` and over
# the whole run.
function temperature_difference(old, new, t1)
    index = Dict(t => k for (k, t) in enumerate(new.t))
    before, all = Float64[], Float64[]
    for (j, t) in enumerate(old.t)
        k = get(index, t, 0)
        k == 0 && continue
        xo, To = live_profile(old, j)
        xn, Tn = live_profile(new, k)
        length(xn) < 2 && continue
        for i in eachindex(xo)
            xn[1] <= xo[i] <= xn[end] || continue
            d = interpolate(xn, Tn, xo[i]) - To[i]
            push!(all, d)
            t < t1 && push!(before, d)
        end
    end
    rms(d) = isempty(d) ? NaN : sqrt(sum(abs2, d) / length(d))
    big(d) = isempty(d) ? NaN : maximum(abs, d)
    return (;
        rms = rms(all),
        largest = big(all),
        rms_window = rms(before),
        largest_window = big(before),
    )
end

trapezoid(t, y, lo, hi) = sum(
    (
        0.5 * (y[j] + y[j+1]) * (t[j+1] - t[j]) for
        j = 1:(length(t)-1) if t[j] >= lo && t[j+1] <= hi
    );
    init = 0.0,
)

# The physical quantities of one solution up to the saved time `t_c`: the
# mass loss [%], the main peak of the gas-generation rate [kg m⁻² s⁻¹] before
# the last event and its time, whether it is interior, the tail [kg/m²], the
# final thickness.
function physical(sol, t_c = sol.t[end])
    keep = findall(<=(t_c), sol.t)
    t = sol.t[keep]
    m = sol.extras.total_mass[keep]
    g = sol.extras.gas_generation_rate[keep]
    ev = filter(<=(t_c), event_times(sol))
    last_event = isempty(ev) ? Inf : maximum(ev)
    before = findall(<(last_event), t)
    isempty(before) && (before = collect(eachindex(t)))
    k = before[argmax(view(g, before))]
    th = sol.extras.thickness_history[keep]
    return (;
        mass_loss = 100 * (1 - m[end] / m[1]),
        peak = g[k],
        t_peak = t[k],
        interior = 1 < k < length(t),
        tail = trapezoid(t, g, t[k], t[end]),
        final_thickness = th[end],
        events = ev,
    )
end

function burnout_time(sol, level)
    th = sol.extras.thickness_history
    k = findfirst(x -> isfinite(x) && x <= level, th)
    return k === nothing ? NaN : sol.t[k]
end

# The case QoIs that measure a quantity against the case's own reference
# value, by the quantity's name.
const CASE_QOI = (
    mass_loss = "mass loss of the gasified",
    peak = "peak MLR",
    t_peak = "time of peak MLR",
    burnout = "burnout time",
)
function case_reference(qois_old, qois_new, quantity)
    (qois_old === nothing || qois_new === nothing) && return false
    prefix = getfield(CASE_QOI, quantity)
    a = [q for q in qois_old if startswith(q.name, prefix)]
    b = [q for q in qois_new if startswith(q.name, prefix)]
    (isempty(a) || isempty(b)) && return false
    return only(a).error < only(a).tol && only(b).error < only(b).tol
end

# One common quantity: both values, the difference (absolute in points for
# the mass loss, relative to the old core's value otherwise), its budget and
# the verdict (`:pass`, `:fail`, `:no_reference`, `:none` where it does not
# apply).
function quantity(name, a, b, budget; relative = true, verdict = true)
    (isfinite(a) && isfinite(b)) ||
        return (; name, old = a, new = b, difference = NaN, budget, verdict = :none)
    d = relative ? abs(b - a) / abs(a) : abs(b - a)
    return (;
        name,
        old = a,
        new = b,
        difference = d,
        budget,
        verdict = verdict ? (d <= budget ? :pass : :fail) : :no_reference,
    )
end

# Over the saved times both solutions hold, up to the last of them (the
# shorter run's end where one stops early).
function common_quantities(r, old, new, t1, qo, qn)
    events = isfinite(t1)
    t_c = maximum(intersect(old.t, new.t))
    po, pn = physical(old, t_c), physical(new, t_c)
    Td = temperature_difference(old, new, t1)
    after(name) = !events || case_reference(qo, qn, name)
    level = max(po.final_thickness, pn.final_thickness) + BURNOUT_MARGIN
    temperature(name, d, budget, judged) = (;
        name,
        old = NaN,
        new = NaN,
        difference = d,
        budget,
        verdict = !isfinite(d) ? :none :
                  !judged ? :no_reference : d <= budget ? :pass : :fail,
    )
    qs = Any[]
    if events
        push!(qs, temperature(:T_rms_window, Td.rms_window, BUDGET.T_rms, true))
        push!(qs, temperature(:T_max_window, Td.largest_window, BUDGET.T_max, true))
    end
    push!(qs, temperature(:T_rms, Td.rms, BUDGET.T_rms, !events))
    push!(qs, temperature(:T_max, Td.largest, BUDGET.T_max, !events))
    push!(
        qs,
        quantity(
            :mass_loss,
            po.mass_loss,
            pn.mass_loss,
            BUDGET.mass_loss;
            relative = false,
            verdict = after(:mass_loss),
        ),
    )
    if r.burnout
        push!(
            qs,
            quantity(
                :burnout,
                burnout_time(old, level),
                burnout_time(new, level),
                BUDGET.burnout;
                verdict = after(:burnout),
            ),
        )
    end
    if po.interior && pn.interior
        push!(qs, quantity(:peak, po.peak, pn.peak, BUDGET.peak; verdict = after(:peak)))
        push!(
            qs,
            quantity(
                :t_peak,
                po.t_peak,
                pn.t_peak,
                BUDGET.t_peak;
                verdict = after(:t_peak),
            ),
        )
        push!(qs, quantity(:tail, po.tail, pn.tail, BUDGET.tail; verdict = !events))
    end
    return qs, po, pn
end

# The white-pine objective's data term of `loss` from each scenario's
# solution: the σ-normalized mean square of the back-face temperature over
# t > 1 s, the replicates' spread floored, the scenarios weighted equally.
function wp_objective(solutions)
    wp = white_pine()
    floor = Base.invokelatest(getglobal, wp.m, :DATA_SIGMA_FLOOR)
    total, back = 0.0, Dict{String,Any}()
    for (label, sol) in solutions
        _, t_exp, T_exp, T_std, _ = wp.dataset[label]
        mask = t_exp .> 1.0
        sim = latest(wp.m, :interpolate_at, t_exp[mask], sol.t, sol.extras.back_temperature)
        r = sim .- T_exp[mask]
        σ = max.(T_std[mask], floor)
        total += mean((r ./ σ) .^ 2) / length(solutions)
        back[label] = sim
    end
    return total, back
end

# `loss` at θ over both scenarios, as the fit evaluates it (the old core's
# `build_and_solve`, `KenCarp4` with Sparspak), its trace in a scratch
# directory.
function wp_loss()
    wp = white_pine()
    Base.invokelatest(getglobal, wp.m, :RUN_DIR)[] = mktempdir()
    data = [wp.dataset[l] for l in ("30kW", "50kW")]
    return quiet(
        () -> latest(wp.m, :loss, wp.theta, data; stage = WP_STAGE, log_trace = false),
    )
end

# ----------------------------------------------------------------------------
# The accuracy part
# ----------------------------------------------------------------------------

function accuracy_row(r)
    cfg = CONFIG[]
    println("[$(r.name)] accuracy")
    flush(stdout)
    po, pn = build(r, "old"), build(r, "new")
    run(core, problem, reltol) = begin
        started = time()
        sol, qois, closure = accuracy_solve(r, core, problem, reltol)
        (; reltol, sol, qois, closure, wall = time() - started)
    end
    old_ref = run("old", po, cfg.reference)
    old_spread = run("old", po, cfg.spread)
    old_default = run("old", po, r.reltol)
    new_ref = run("new", pn, cfg.reference)
    new_spread = run("new", pn, cfg.spread)
    rungs = sort(unique(vcat(collect(cfg.ladder), r.reltol)); rev = true)
    ladder = [run("new", pn, x) for x in rungs]
    solves = vcat([old_ref, old_spread, old_default, new_ref, new_spread], ladder)
    # A solve that stops before the end of its span stops at an event of its
    # own (an interior void, the last cell, the thickness floor).
    t_span = po.tspan[2]
    stop(sol) = sol.t[end] < t_span * (1 - 1e-12) ? sol.t[end] : Inf
    t1 = minimum(s -> min(minimum(event_times(s.sol); init = Inf), stop(s.sol)), solves)
    window(s, ref) = window_errors(s.sol, ref.sol, t1, r.field)
    target = window(old_default, old_ref).rms
    rung_rows = map(ladder) do s
        w = window(s, new_ref)
        (;
            s.reltol,
            retcode = s.sol.retcode,
            outcome = outcome(s.sol),
            t_end = s.sol.t[end],
            window_rms = w.rms,
            window_largest = w.largest,
            s.wall,
            counters(s.sol)...,
        )
    end
    # Matched: the loosest rung on the window's target whose own closure is
    # within its budget; the loosest on the window alone is reported beside.
    on_target = [
        k for (k, x) in enumerate(rung_rows) if
        finished(ladder[k].sol) && x.window_rms <= target
    ]
    loosest(ks) = isempty(ks) ? nothing : ks[argmax([rung_rows[k].reltol for k in ks])]
    window_only = loosest(on_target)
    m = loosest([k for k in on_target if closes(ladder[k].closure)])
    matched = m === nothing ? nothing : ladder[m]
    common, phys_old, phys_new =
        matched === nothing ? (Any[], nothing, nothing) :
        common_quantities(
            r,
            old_default.sol,
            matched.sol,
            t1,
            old_default.qois,
            matched.qois,
        )
    if matched !== nothing
        for (name, c) in
            ((:closure_old, old_default.closure), (:closure_new, matched.closure))
            push!(
                common,
                (;
                    name,
                    old = NaN,
                    new = NaN,
                    difference = abs(c.mass),
                    budget = c.budget,
                    verdict = isnan(c.budget) || !isfinite(c.mass) ? :none :
                              closes(c) ? :pass : :fail,
                ),
            )
        end
    end
    qoi_summary(s) =
        s === nothing || s.qois === nothing ? nothing :
        [(; q.name, q.value, q.error, q.tol, pass = q.error < q.tol) for q in s.qois]
    w_old = window(old_default, old_ref)
    result = (;
        name = r.name,
        cells = size(old_default.sol.T, 1),
        reltol = r.reltol,
        ratio = r.ratio,
        t1,
        events = isfinite(t1),
        spread = (
            old = window(old_spread, old_ref).rms,
            new = window(new_spread, new_ref).rms,
        ),
        target,
        old_window_largest = w_old.largest,
        old = (;
            retcode = old_default.sol.retcode,
            t_end = old_default.sol.t[end],
            old_default.wall,
            closure = old_default.closure,
            qois = qoi_summary(old_default),
            counters(old_default.sol)...,
        ),
        ladder = rung_rows,
        matched_reltol = matched === nothing ? NaN : matched.reltol,
        window_reltol = window_only === nothing ? NaN : rung_rows[window_only].reltol,
        matched = matched === nothing ? nothing :
                  (;
            closure = matched.closure,
            qois = qoi_summary(matched),
            outcome = outcome(matched.sol),
        ),
        common,
        physical = (; old = phys_old, new = phys_new),
        references = (
            old = old_ref.wall,
            new = new_ref.wall,
            old_spread = old_spread.wall,
            new_spread = new_spread.wall,
        ),
        objective = nothing,
    )
    if r.scenario !== nothing
        result = (;
            result...,
            objective = (;
                label = r.scenario,
                old_default = old_default.sol,
                matched = matched === nothing ? nothing : matched.sol,
                new_default = ladder[findfirst(x -> x.reltol == r.reltol, ladder)].sol,
            ),
        )
    end
    println(
        "[$(r.name)] target $(_g(target)) matched $(result.matched_reltol) ",
        "verdicts $(join(("$(q.name)=$(q.verdict)" for q in common), " ")) ",
        "references $(round.(values(result.references); digits = 1)) s",
    )
    flush(stdout)
    return result
end

# The solutions are dropped before a part is written.
slim(x) = x.objective === nothing ? x : (; x..., objective = (; label = x.objective.label))

function accuracy_part(k, dir)
    load = Sys.loadavg()[1]
    started = time()
    rows = [row_named(n) for n in ROW_NAMES if part_of(n) == k]
    results = [accuracy_row(r) for r in rows]
    mkpath(dir)
    serialize(
        joinpath(dir, "accuracy$(k).jls"),
        (;
            results = map(slim, results),
            stamp = COMMIT_STAMP,
            load = (load, Sys.loadavg()[1]),
            elapsed = time() - started,
        ),
    )
    # The objective needs both scenarios' solutions, held by the parts that
    # solve them; each writes its own, with the scenario's solve as fitted,
    # and the merge combines them. The part of the 30 kW/m² scenario also
    # evaluates `loss` itself over both.
    for x in results
        x.objective === nothing && continue
        label = x.objective.label
        wp = white_pine()
        fitted = quiet(
            () -> latest(
                wp.m,
                :build_and_solve,
                wp.theta,
                first(wp.dataset[label]);
                stage = WP_STAGE,
            ),
        )
        serialize(
            joinpath(dir, "objective_$(label).jls"),
            (;
                label,
                old_default = series(x.objective.old_default),
                new_default = series(x.objective.new_default),
                matched = x.objective.matched === nothing ? nothing :
                          series(x.objective.matched),
                fitted = series(fitted),
                loss = label == "30kW" ? wp_loss() : NaN,
                stamp = COMMIT_STAMP,
            ),
        )
    end
    println("wrote accuracy part $(k) to $(dir); $(round(Int, time() - started)) s")
end

# What the objective reads of a solution.
series(sol) = (; t = sol.t, extras = (; back_temperature = sol.extras.back_temperature))

# Every part must have been measured at this commit; the smoke run, which
# spans the commits it tests, is only warned.
function check_stamp(stamp, what)
    stamp == COMMIT_STAMP && return nothing
    message = "$(what) was measured at $(stamp), not at $(COMMIT_STAMP)"
    CONFIG[] === SMOKE ? (@warn message) : error(message)
    return nothing
end

function read_parts(dir, pattern)
    files = sort(filter(f -> occursin(pattern, f), readdir(dir)))
    parts = [deserialize(joinpath(dir, f)) for f in files]
    foreach(p -> check_stamp(p.stamp, "$(dir): a part"), parts)
    return parts
end

function accuracy_results(dir)
    parts = read_parts(dir, r"^accuracy\d\.jls$")
    return vcat((p.results for p in parts)...), parts
end

function objective_from_parts(dir)
    parts = read_parts(dir, r"^objective_.*\.jls$")
    length(parts) == 2 || return nothing
    any(p -> p.matched === nothing, parts) && return nothing
    by(f) = Dict(p.label => getfield(p, f) for p in parts)
    old_default, _ = wp_objective(by(:old_default))
    new_default, nb = wp_objective(by(:new_default))
    new_matched, _ = wp_objective(by(:matched))
    fitted, fb = wp_objective(by(:fitted))
    back = Dict(
        l => (;
            rms = sqrt(mean(abs2, nb[l] .- fb[l])),
            largest = maximum(abs, nb[l] .- fb[l]),
        ) for l in keys(nb)
    )
    fitted_loss = only(p.loss for p in parts if p.label == "30kW")
    return (; fitted_loss, fitted, old_default, new_default, new_matched, back)
end

# ----------------------------------------------------------------------------
# The timing part
# ----------------------------------------------------------------------------

# The warm-up, then the median of `repeats` solves (three where the warm-up
# exceeds `SLOW_WARMUP`), their range, the allocated bytes of the last and
# the counters.
function timed(f, repeats)
    warm = @elapsed f()
    n = warm > SLOW_WARMUP ? min(repeats, 3) : repeats
    times, bytes, sol = Float64[], 0, nothing
    for _ = 1:n
        GC.gc()
        s = @timed f()
        push!(times, s.time)
        bytes, sol = s.bytes, s.value
    end
    return (;
        median = median(times),
        low = minimum(times),
        high = maximum(times),
        n,
        bytes,
        counters = counters(sol),
        retcode = sol.retcode,
    )
end

function timing_row(r, acc)
    cfg = CONFIG[]
    println("[$(r.name)] timing")
    flush(stdout)
    po, pn = build(r, "old"), build(r, "new")
    old = timed(() -> plain_solve(r, "old", po, r.reltol), cfg.old_repeats)
    default = timed(() -> plain_solve(r, "new", pn, r.reltol), cfg.new_repeats)
    matched =
        isfinite(acc.matched_reltol) ?
        (
            acc.matched_reltol == r.reltol ? default :
            timed(() -> plain_solve(r, "new", pn, acc.matched_reltol), cfg.new_repeats)
        ) : nothing
    return (; name = r.name, old, default, matched, matched_reltol = acc.matched_reltol)
end

const RECHECK = ("V-BO-01", "ale 100, solid λ 1.0e-6", "Douglas fir, 60 s")

function recheck_rows()
    cfg = CONFIG[]
    out = Any[]
    for name in RECHECK
        r = name == "Douglas fir, 60 s" ? douglas_fir_row(60.0) : row_named(name)
        pn = build(r, "new")
        push!(
            out,
            (; name, timed(() -> plain_solve(r, "new", pn, r.reltol), cfg.recheck)...),
        )
        println("[$(name)] recheck $(_g(out[end].median)) s")
        flush(stdout)
    end
    return out
end

function timing_part(dir)
    load = Sys.loadavg()
    started = time()
    acc, _ = accuracy_results(dir)
    results = Any[]
    for name in ROW_NAMES
        r = row_named(name)
        a = only(x for x in acc if x.name == r.name)
        push!(results, timing_row(r, a))
        x = results[end]
        println(
            "[$(r.name)] old $(_g(x.old.median)) s, new default $(_g(x.default.median)) s",
            x.matched === nothing ? ", unmatched" :
            ", matched $(x.matched_reltol) $(_g(x.matched.median)) s",
        )
        flush(stdout)
    end
    recheck = recheck_rows()
    serialize(
        joinpath(dir, "timing.jls"),
        (;
            results,
            recheck,
            stamp = COMMIT_STAMP,
            load = (load, Sys.loadavg()),
            elapsed = time() - started,
            cpus = affinity(),
        ),
    )
    println("wrote the timing part to $(dir); $(round(Int, time() - started)) s")
end

# The CPUs this process may run on, as `taskset` set them.
function affinity()
    try
        return strip(last(split(read(`taskset -cp $(getpid())`, String), ":")))
    catch
        return "unknown"
    end
end

# ----------------------------------------------------------------------------
# The memory mode
# ----------------------------------------------------------------------------

# One core on one row in this process: a warm-up, then the measured solve.
function memory_one(name, core, reltol, file)
    r = row_named(name; core)
    p = build(r, core)
    plain_solve(r, core, p, reltol)
    GC.gc()
    s = @timed plain_solve(r, core, p, reltol)
    serialize(
        file,
        (;
            name,
            core,
            reltol,
            bytes = s.bytes,
            rss = Sys.maxrss(),
            time = s.time,
            stamp = COMMIT_STAMP,
        ),
    )
end

function memory_part(dir, jobs)
    acc, _ = accuracy_results(dir)
    tasks = Any[]
    for x in acc, core in ("old", "new")
        reltol = core == "old" ? x.reltol : x.matched_reltol
        isfinite(reltol) || continue
        push!(tasks, (x.name, core, reltol))
    end
    mkpath(joinpath(dir, "memory"))
    julia = Base.julia_cmd()
    project = Base.active_project()
    queue = Channel{Any}(length(tasks))
    foreach(t -> put!(queue, t), tasks)
    close(queue)
    @sync for _ = 1:jobs
        @async for (name, core, reltol) in queue
            file = joinpath(
                dir,
                "memory",
                "$(replace(name, r"[^A-Za-z0-9]+" => "_"))_$(core).jls",
            )
            cmd = `$julia --project=$project --startup-file=no --threads=1 $(@__FILE__) memory-one $name $core $reltol $file`
            CONFIG[] === SMOKE && (cmd = `$cmd --smoke`)
            ok = success(pipeline(cmd; stdout = devnull, stderr = devnull))
            println("[$(name)] memory $(core): ", ok ? "done" : "FAILED")
            flush(stdout)
        end
    end
    results =
        [deserialize(joinpath(dir, "memory", f)) for f in readdir(joinpath(dir, "memory"))]
    serialize(joinpath(dir, "memory.jls"), (; results, stamp = COMMIT_STAMP))
end

# ----------------------------------------------------------------------------
# The tables
# ----------------------------------------------------------------------------

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
_v(v) =
    v === :pass ? "pass" :
    v === :fail ? "**fail**" : v === :no_reference ? "no common reference" : ""

accuracy_passes(x) =
    !isempty(x.common) && all(q -> q.verdict in (:pass, :none, :no_reference), x.common)

function performance_verdict(x, t)
    t.matched === nothing && return "unmatched"
    accuracy_passes(x) || return "not read (accuracy)"
    ratio = t.matched.median / t.old.median
    return ratio <= 1.5 ? "pass" : ratio <= 2 ? "owner" : "**over 2**"
end

function write_tables(dir_out, acc, acc_parts, timing, memory, objective)
    mkpath(dir_out)
    tim = Dict(t.name => t for t in timing.results)
    mem = Dict((m.name, m.core) => m for m in memory.results)
    open(joinpath(dir_out, "cutover_comparison.csv"), "w") do io
        println(io, "row,quantity,old,new,difference,budget,verdict")
        for x in acc, q in x.common
            println(
                io,
                "\"$(x.name)\",$(q.name),$(q.old),$(q.new),$(q.difference),$(q.budget),$(q.verdict)",
            )
        end
    end
    loads(ps) = join((round(p.load[1]; digits = 1) for p in ps), ", ")
    open(joinpath(dir_out, "cutover_comparison.md"), "w") do io
        println(io, "# The end-to-end comparison of the two cores\n")
        println(
            io,
            "Produced by `cutover_comparison.jl` at $(COMMIT_STAMP). Julia $(VERSION), one thread per process. ",
            "Accuracy: $(length(acc_parts)) processes at a one-minute load of $(loads(acc_parts)) at their starts. ",
            "Timing: one process on CPUs `$(timing.cpus)`, one-minute load $(round(timing.load[1][1]; digits = 1)) at the start and ",
            "$(round(timing.load[2][1]; digits = 1)) at the end, $(round(Int, timing.elapsed)) s. ",
            "Memory: one process per core and row. The script's header defines every row, measure and budget.\n",
        )
        println(io, "## The rows\n")
        println(
            io,
            "Window: the saved times before the first event of any solve of the row [s] (`Inf`: none). Target: the old ",
            "default's window RMS against its own reltol 1e-9 solve [K; relative thickness on V-BO-04]; matched: the ",
            "finite-volume core's loosest rung at or below it whose own closure is within its budget (on the window ",
            "alone: the loosest at or below it). Accuracy: every verdict of the common quantities (below). ",
            "Timing: warm medians [s]; ratio: matched over the old default.\n",
        )
        println(
            io,
            "| row | cells | reltol | window [s] | target | spreads old, new | matched reltol (on the window alone) | accuracy | old [s] | new matched [s] | ratio | equal-tolerance ratio | verdict |",
        )
        println(io, "|" * repeat("---|", 13))
        for x in acc
            t = tim[x.name]
            ratio = t.matched === nothing ? NaN : t.matched.median / t.old.median
            println(
                io,
                "| $(x.name) | $(x.cells) | $(x.reltol) | $(_g(x.t1; d = 4)) | $(_g(x.target)) | ",
                "$(_g(x.spread.old)), $(_g(x.spread.new)) | $(_g(x.matched_reltol)) ($(_g(x.window_reltol))) | ",
                isempty(x.common) ? "unmatched" :
                (accuracy_passes(x) ? "pass" : "**fail**"),
                " | ",
                "$(_g(t.old.median)) | $(t.matched === nothing ? "" : _g(t.matched.median)) | $(_g(ratio)) | ",
                "$(_g(t.default.median / t.old.median)) | $(performance_verdict(x, t)) |",
            )
        end
        println(io, "\n## The common quantities\n")
        println(
            io,
            "Between the old default and the finite-volume core's matched solve. Temperatures [K] at the old core's live ",
            "cell centres, over the window (the whole run on a row without events) and over the whole run; the mass loss ",
            "[%] and its difference in points; the burnout time [s], the main peak [g m⁻² s⁻¹] and its time [s], the tail ",
            "[g/m²], each with its difference relative to the old core's value; each core's closure against its case's ",
            "ledger tolerance.\n",
        )
        println(io, "| row | quantity | old | new | difference | budget | verdict |")
        println(io, "|---|---|---|---|---|---|---|")
        scale(q) = q.name in (:peak, :tail) ? 1000 : 1
        for x in acc, q in x.common
            println(
                io,
                "| $(x.name) | $(q.name) | $(_g(scale(q) * q.old; d = 5)) | $(_g(scale(q) * q.new; d = 5)) | ",
                "$(_g(q.difference)) | $(_g(q.budget)) | $(_v(q.verdict)) |",
            )
        end
        println(io, "\n## Closures and case QoIs\n")
        println(
            io,
            "Each core's own closure (relative): the finite-volume core's mass and energy from its conservation recorder ",
            "(mass from the case's ledger on the cases), the old core's mass from the saved-time trapezoid. QoIs: passing of ",
            "counted, the failing named.\n",
        )
        println(
            io,
            "| row | old: retcode, mass | old QoIs | new: outcome, mass, energy | new QoIs |",
        )
        println(io, "|---|---|---|---|---|")
        qs(q) =
            q === nothing ? "" :
            "$(count(z -> z.pass, q))/$(length(q))" *
            (all(z -> z.pass, q) ? "" : " ($(join((z.name for z in q if !z.pass), "; ")))")
        for x in acc
            m = x.matched
            println(
                io,
                "| $(x.name) | $(x.old.retcode), $(_g(x.old.closure.mass)) | $(qs(x.old.qois)) | ",
                m === nothing ? " | |" :
                "$(m.outcome), $(_g(m.closure.mass)), $(_g(m.closure.energy)) | $(qs(m.qois)) |",
            )
        end
        println(io, "\n## The counters\n")
        println(
            io,
            "Of the timed solves: accepted (rejected) steps, right-hand sides, Jacobians, linear solves, nonlinear ",
            "failures, domain rejections (the old core counts none), events and invocations; the old default and the ",
            "finite-volume core at its matched and default rungs.\n",
        )
        println(
            io,
            "| row | solve | reltol | warm [s] (range, n) | steps (rejected) | RHS | Jacobians | linear solves | nonlinear failures | domain rejections | events, invocations |",
        )
        println(io, "|" * repeat("---|", 11))
        for x in acc
            t = tim[x.name]
            for (label, s, rt) in (
                ("old default", t.old, x.reltol),
                ("matched", t.matched, x.matched_reltol),
                ("new default", t.default, x.reltol),
            )
                s === nothing && continue
                c = s.counters
                println(
                    io,
                    "| $(x.name) | $(label) | $(_g(rt)) | $(_g(s.median)) ($(_g(s.low))–$(_g(s.high)), $(s.n)) | ",
                    "$(c.steps) ($(c.rejected)) | $(c.rhs) | $(c.jacobians) | $(c.linear_solves) | $(c.nonlinear_failures) | ",
                    "$(c.domain_rejections < 0 ? "" : c.domain_rejections) | $(c.events), $(c.invocations) |",
                )
            end
        end
        println(io, "\n## Memory\n")
        println(
            io,
            "Each core on each row in a process of its own after one warm-up solve, at the timing's tolerance: the ",
            "measured solve's allocated bytes [MiB] and the process's peak resident set [MiB], and their ratios.\n",
        )
        println(
            io,
            "| row | old allocated | new allocated | ratio | old peak RSS | new peak RSS | ratio |",
        )
        println(io, "|---|---|---|---|---|---|---|")
        mib(b) = b / 2^20
        for x in acc
            o, n = get(mem, (x.name, "old"), nothing), get(mem, (x.name, "new"), nothing)
            (o === nothing || n === nothing) && continue
            println(
                io,
                "| $(x.name) | $(_g(mib(o.bytes))) | $(_g(mib(n.bytes))) | $(_g(n.bytes / o.bytes)) | ",
                "$(_g(mib(o.rss))) | $(_g(mib(n.rss))) | $(_g(n.rss / o.rss)) |",
            )
        end
        if objective !== nothing
            println(io, "\n## The white-pine objective\n")
            println(
                io,
                "The data term of `loss` at the fitted stage-2 conductivities over the 30 and 50 kW/m² scenarios: as ",
                "fitted (the old core's `KenCarp4` with Sparspak; `loss` itself gives $(_g(objective.fitted_loss; d = 6))), ",
                "the old core at its default, the finite-volume core at its default and matched rungs; the back-face ",
                "temperature of the finite-volume default against the fitted solve over the data times [K].\n",
            )
            println(
                io,
                "| as fitted | old default | new default | new matched | drift (new default − as fitted) |",
            )
            println(io, "|---|---|---|---|---|")
            println(
                io,
                "| $(_g(objective.fitted; d = 6)) | $(_g(objective.old_default; d = 6)) | $(_g(objective.new_default; d = 6)) | ",
                "$(_g(objective.new_matched; d = 6)) | $(_g(objective.new_default - objective.fitted; d = 3)) |",
            )
            for (l, b) in sort(collect(objective.back))
                println(
                    io,
                    "\nBack-face temperature, $(l): RMS $(_g(b.rms)) K, largest $(_g(b.largest)) K.",
                )
            end
        end
        println(io, "\n## The recheck configurations\n")
        println(
            io,
            "The finite-volume core at the deck's tolerance, warm median, range and count, and allocated bytes [MiB].\n",
        )
        println(
            io,
            "| configuration | warm [s] | range [s] | n | allocated [MiB] | steps |",
        )
        println(io, "|---|---|---|---|---|---|")
        for x in timing.recheck
            println(
                io,
                "| $(x.name) | $(_g(x.median)) | $(_g(x.low))–$(_g(x.high)) | $(x.n) | $(_g(mib(x.bytes))) | $(x.counters.steps) |",
            )
        end
    end
end

function merge_parts(dir; out = OUT)
    acc, acc_parts = accuracy_results(dir)
    timing = deserialize(joinpath(dir, "timing.jls"))
    check_stamp(timing.stamp, "the timing part")
    memory = deserialize(joinpath(dir, "memory.jls"))
    check_stamp(memory.stamp, "the memory part")
    acc = [
        only(x for x in acc if x.name == n) for n in ROW_NAMES if any(x -> x.name == n, acc)
    ]
    write_tables(out, acc, acc_parts, timing, memory, objective_from_parts(dir))
    println("wrote $(out)/cutover_comparison.{csv,md}")
end

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

function comparison_main(args)
    smoke = "--smoke" in args
    args = filter(!=("--smoke"), args)
    smoke && (CONFIG[] = SMOKE)
    mode = isempty(args) ? "" : args[1]
    if mode == "accuracy" && length(args) == 3
        accuracy_part(parse(Int, args[2]), args[3])
    elseif mode == "timing" && length(args) == 2
        timing_part(args[2])
    elseif mode == "memory" && length(args) in (2, 3)
        memory_part(args[2], length(args) == 3 ? parse(Int, args[3]) : 4)
    elseif mode == "memory-one" && length(args) == 5
        memory_one(args[2], args[3], parse(Float64, args[4]), args[5])
    elseif mode == "merge" && length(args) == 2
        merge_parts(args[2]; out = smoke ? args[2] : OUT)
    else
        error(
            "usage: cutover_comparison.jl accuracy K DIR | timing DIR | memory DIR [JOBS] | merge DIR, each with --smoke for the reduced matrix",
        )
    end
    return nothing
end

abspath(PROGRAM_FILE) == @__FILE__() && comparison_main(ARGS)
