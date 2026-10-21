# Gate items 1, 2, 4, 6 and 10 of the Phase 2 plan (§7.3) and R25's tail
# diagnostic: the burnout cases on the finite-volume core with the
# equidistribution correction, across grids, thresholds, `τ` and the monitor.
#
# The matrix, every run in both energy forms, the band at its default (0.1),
# the monitor without its mass term (β = 0) unless stated:
#
#   V-BO-01 at 30, 60 and 120 PMMA cells × ε_m ∈ {3e-3, 1e-2, 3e-2}, the
#   default τ;
#   V-BO-01 at 60 cells and ε_m = 1e-2 with τ = 360 and 3600 s (τ is the
#   layer's equidistribution time; the prototype's per-cell τ_p = τ/N² is
#   0.1, 0.2 and 1 s at 60 cells for τ = 360, 720 and 3600 s);
#   V-BO-03 at 30, 60 and 120 cells at the defaults, run past the case's
#   500 s to `VBO03_T_EXTENDED`, after its layer's invocation;
#   V-BO-04 at its 40 cells without and with the correction;
#   V-BO-06 (the graded mesh) at 60 cells and V-BO-08 (the thin layer) at
#   its 10 cells, each × the three thresholds; V-BO-08 at ε_m = 1e-2 also
#   without the correction, for its substrate before the first activation;
#   V-BO-07 (the mixed layer) at 60 cells with β = 0 and β = 1.
#
# Per run: the outcome; the deactivation invocations, each with its cells
# below ε_m, the cells the band took and its largest merged mass fraction
# (item 1: exactly one invocation per layer, the largest fraction at most
# 5 ε_m; a run with no invocation fails); every QoI of the case in assert mode against its fixed references,
# with the old core at the same resolution beside (item 2; its burst criteria
# are `burnout_burst.jl`'s); events, invocations and inadmissible results
# (item 4); the largest scaled constraint residual of the accepted states
# under C52's scale and under Phase 1's, which omits the sweep of the node
# corrections, and the accepted states Phase 1's would have rejected beside
# the domain test's rejections (item 6); the smallest live cell at the end of
# every step, before a deactivation event, over its initial thickness against
# 0.5 ε_m (item 10); the warm wall time (median of three solves without the
# recorders), the old core's at the same resolution beside. In the
# temperature form the plain solves take the instrumented solve's steps; in
# the enthalpy form a residual evaluation sets the warm start of the next
# caloric inversion, so the recorders move the instrumented solve off the
# plain one at round-off.
#
# The recorders read each accepted state before the deactivation event of its
# step: a continuous callback's condition sees the state at the end of the
# step before any discrete callback runs, and never crosses zero. On a step
# whose event is located they read the step's end before the event is
# located, a state after the event's time that the trajectory does not pass
# through: the thinnest cell there is thinner than at the event (item 10's
# reading is conservative), and the residuals are those of a proposal the
# domain test accepted rather than of a state on the trajectory (item 6).
#
# R25's tail diagnostic, for every V-BO-01 run: at the state that ends the
# step of the layer's invocation, before the merge, the PMMA mass `mᵢ` and
# the temperature `Tᵢ` of the merged cells and the receiver; after it the
# residue and the merged temperature `T_m`; the rate
# `Σ mᵢ k(Tᵢ) → m_res k(T̄) → m_res k(T_m)`, `T̄` the mass-weighted mean
# (lumping, the merged temperature); the mass gasified from the event to the end of the run
# against the residue; the largest mass loss rate after the event.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_ladder.jl [--item1]
# Writes test/newcore/gates/measurements/burnout_ladder.{csv,md}. With
# `--item1` nothing is solved: the item 1 verdicts are recomputed from the
# table's csv, which holds every quantity the criterion reads, and written
# into its item 1 column and the matching rows of the md.

ENV["PYROLYSIS_CORE"] = "new"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using Statistics: median
using OrdinaryDiffEq: ContinuousCallback, DiscreteCallback, CallbackSet

const Ig = Pyrolysis.Integration
const St = Pyrolysis.States
const As = Pyrolysis.Assembly
const Ke = Pyrolysis.Kernels
const Cl = Pyrolysis.Closures

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const FORMS = (:temperature, :enthalpy)
const REPEATS = 3
const TAU = Cl.EQUIDISTRIBUTION_TIME
const CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-03", "V-BO-04", "V-BO-06", "V-BO-07", "V-BO-08")
)
const THRESHOLDS = (3e-3, 1e-2, 3e-2)

# A run of the matrix: the case at `n` cells and threshold `ε_m`, with the
# equidistribution correction at `τ` and `monitor = (α, β)` or without it, to
# the case's end or to `t_end`.
run_row(
    case,
    n;
    ε_m = 1e-2,
    τ = TAU,
    corrected = true,
    monitor = (1.0, 0.0),
    t_end = nothing,
) = (; case, n, ε_m, τ, corrected, monitor, t_end)

const ROWS = [
    [run_row("V-BO-01", n; ε_m) for n in (30, 60, 120) for ε_m in THRESHOLDS]
    [run_row("V-BO-01", 60; τ) for τ in (360.0, 3600.0)]
    [run_row("V-BO-03", n; t_end = VBO03_T_EXTENDED) for n in (30, 60, 120)]
    [run_row("V-BO-04", 40; corrected) for corrected in (false, true)]
    [run_row("V-BO-06", 60; ε_m) for ε_m in THRESHOLDS]
    [run_row("V-BO-07", 60; monitor) for monitor in ((1.0, 0.0), (1.0, 1.0))]
    [run_row("V-BO-08", VBO08_N_PMMA; ε_m) for ε_m in THRESHOLDS]
    [run_row("V-BO-08", VBO08_N_PMMA; corrected = false)]
]

# ----------------------------------------------------------------------------
# The recorders
# ----------------------------------------------------------------------------

# What the recorders saw over a solve.
mutable struct Probe
    min_h::Float64
    min_ratio::Float64
    c52::Float64
    phase1::Float64
    phase1_above::Int
    mismatch::Float64
    steps::Int
    times::Vector{Float64}
    first_active::Float64
    stash::Vector{Float64}
    stash_t::Float64
    events::Int
    tails::Vector{Any}
end
Probe(n) = Probe(Inf, Inf, 0.0, 0.0, 0, 0.0, 0, Float64[], NaN, zeros(n), NaN, 0, Any[])

# The largest scaled residual of the algebraic rows at the derived state `ds`
# and contributions `c` of the residual in `check.du`, as
# `Integration.constraint_residual` computes it; `sweep = false` leaves the
# node corrections' sweep out of the balance rows' scale (Phase 1's scale).
function scaled_residual(check, ds, c; sweep::Bool)
    layout = check.def.layout
    du, tol, rtol = check.du, check.abstol, check.reltol
    mask = layout.mask
    worst = 0.0
    if St.is_moving(layout)
        g = ds.geometry
        topo = check.def.setup.topology
        for nd = 1:St.n_nodes(layout)
            row = St.velocity_index(layout, nd)
            mask[row] || continue
            k = check.balance_cell[nd]
            if k == 0
                S, a = abs(ds.w[nd]), tol[row]
            else
                S = (abs(c.θ[k]) + abs(ds.θ_A[])) * ds.V_eval[k]
                for f in topo.cells[k].faces
                    S += g.face_area[f] * abs(ds.w[topo.faces[f].nodes[1]])
                end
                if sweep
                    for f in topo.cells[k].faces
                        S += g.face_area[f] * abs(ds.c[topo.faces[f].nodes[1]])
                    end
                end
                a = tol[row] * g.face_area[nd]
            end
            worst = max(worst, abs(du[row]) / (rtol * S + a))
        end
    end
    if St.has_intensities(layout)
        for f = 1:St.n_nodes(layout)
            row = St.intensity_index(layout, f)
            mask[row] || continue
            worst = max(worst, abs(du[row]) / (rtol * check.intensity_scale + tol[row]))
        end
    end
    St.has_chi(layout) && error("the ladder's decks carry no conversion history")
    return worst
end

# One accepted state, before any discrete callback: the thinnest live cell
# and the scaled residuals.
function probe_step!(p::Probe, check, h0, u, t)
    layout = check.def.layout
    for i = 1:St.n_cells(layout)
        layout.mask[St.energy_index(layout, i)] || continue
        h = u[St.position_index(layout, i + 1)] - u[St.position_index(layout, i)]
        p.min_h = min(p.min_h, h)
        p.min_ratio = min(p.min_ratio, h / h0[i])
    end
    if !Ig.inadmissible(u, layout, check.abstol)
        As.residual!(check.du, u, check.def, t, check.ws)
        ds, c = As.buffers(check.ws, check.def, Float64)
        package = Ig.constraint_residual(check, ds, c, t)
        c52 = scaled_residual(check, ds, c; sweep = true)
        phase1 = scaled_residual(check, ds, c; sweep = false)
        p.mismatch = max(p.mismatch, abs(c52 - package) / max(package, eps()))
        p.c52 = max(p.c52, c52)
        p.phase1 = max(p.phase1, phase1)
        phase1 > 1 && (p.phase1_above += 1)
        isnan(p.first_active) && any(!iszero, ds.c) && (p.first_active = t)
    end
    p.steps += 1
    push!(p.times, t)
    copyto!(p.stash, u)
    p.stash_t = t
    return nothing
end

# The recorders of a set-up problem: the probe at the end of every step, and
# after the deactivation event of a step the states on both sides of it.
function instrument(cp, p::Probe)
    check = Ig.DomainCheck(cp)
    h0 = cp.def.setup.reference.thickness
    cells = cp.deactivation
    spy = ContinuousCallback(
        (u, t, integ) -> begin
            t == integ.t && probe_step!(p, check, h0, u, t)
            return 1.0
        end,
        integ -> nothing;
        save_positions = (false, false),
        interp_points = 0,
    )
    cells === nothing && return spy
    after = DiscreteCallback(
        (u, t, integ) -> length(cells.events) > p.events,
        integ -> begin
            push!(
                p.tails,
                (;
                    t = integ.t,
                    pre = copy(p.stash),
                    pre_t = p.stash_t,
                    post = copy(integ.u),
                    events = (p.events+1):length(cells.events),
                ),
            )
            p.events = length(cells.events)
        end;
        save_positions = (false, false),
    )
    return CallbackSet(spy, after)
end

# ----------------------------------------------------------------------------
# The solves
# ----------------------------------------------------------------------------

# The new core's keywords of a run beside the case's own. `setup` takes the
# layer interfaces as the Lagrangian nodes.
function run_keywords(row, form)
    correction =
        row.corrected ? Cl.EquidistributionCorrection(τ = row.τ, monitor = row.monitor) :
        Cl.NoCorrection()
    return (; energy_state = form, threshold = row.ε_m, correction)
end

# The deck of a run: the case's own, or V-BO-03's run on to `t_end`.
function deck(row)
    row.t_end === nothing && return Base.invokelatest(CASES[row.case].build, row.n)
    row.case == "V-BO-03" || error("only V-BO-03 runs past its case's end")
    return Base.invokelatest(_vbo_problem, row.n; substrate_alpha = 0.0, t_end = row.t_end)
end

# The case's measures, the mass ledger the case's recorder's; the variant
# decks read their substrate's thickness from the run, as their QoIs do.
function measures(case, sol)
    case == "V-BO-04" && return nothing
    case in ("V-BO-06", "V-BO-07", "V-BO-08") ||
        return Base.invokelatest(_vbo_measures, sol; id = case)
    l_layer = case == "V-BO-08" ? VBO08_L_PMMA : VBO_L_PMMA
    l_kw = first(sol.extras.thickness_history) - l_layer
    return Base.invokelatest(_vbo_measures, sol; l_kw, id = case)
end

# V-BO-08's substrate: the temperatures of the cells that start without PMMA
# at the saved times.
function substrate_history(sol)
    kw = findall(iszero, view(sol.ξ, 1, :, 1))
    return (; t = sol.t, T = sol.T[kw, :])
end

# The case's keywords in the new core's names.
function core_keywords(case)
    s = case.solve_kwargs
    rad = haskey(s, :radiation_model) ? (; radiation_model = s.radiation_model) : (;)
    return (;
        deactivation = s.handle_depletion,
        use_ale = s.use_ale,
        abstol = s.abstol,
        reltol = s.reltol,
        rad...,
    )
end

# The instrumented solve: through the case's own solve function where it has
# one, so its QoIs read what its recorder saw.
function instrumented_solve(row, form)
    case = CASES[row.case]
    problem = deck(row)
    kw = run_keywords(row, form)
    probe = Ref{Probe}()
    hook = cp -> (probe[] = Probe(length(cp.u0)); instrument(cp, probe[]))
    counter = _VBOWarnCounter(Logging.current_logger(), r"inadmissible", Ref(0))
    captured = Ref{Any}()
    sol = Logging.with_logger(counter) do
        if case.solve_fn === nothing
            cp = Ig.setup(problem; core_keywords(case)..., kw...)
            captured[] = cp
            Ig.integrate(
                cp;
                callback = hook(cp),
                min_thickness = case.solve_kwargs.min_thickness,
            )
        else
            Base.invokelatest(
                case.solve_fn,
                problem;
                case.solve_kwargs...,
                kw...,
                instrument = cp -> (captured[] = cp; hook(cp)),
            )
        end
    end
    return (; sol, cp = captured[], probe = probe[], inadmissible = counter.count[])
end

# The same solve without recorders, for its wall time.
function plain_solve(row, form)
    case = CASES[row.case]
    cp = Ig.setup(deck(row); core_keywords(case)..., run_keywords(row, form)...)
    return Ig.integrate(cp; min_thickness = case.solve_kwargs.min_thickness)
end

old_solve(row) = Logging.with_logger(Logging.NullLogger()) do
    solve(deck(row); progress = false, CASES[row.case].solve_kwargs...)
end

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

# ----------------------------------------------------------------------------
# The measures
# ----------------------------------------------------------------------------

# R25's decomposition at the layer's invocation (the first).
function tail_terms(cp, sol, tail, old)
    def, ws = cp.def, cp.ws
    layout, setup = def.layout, def.setup
    events = cp.deactivation.events[tail.events]
    receiver = first(events).receiver
    cells = vcat([e.cell for e in events], receiver)
    rx = only(setup.material.reactions)
    k(T) = rx.A * exp(-rx.E / (Pyrolysis.Materials.R_GAS * T))
    area = setup.reference.cross_section
    ds = ws.derived
    St.derive!(ds, tail.pre, layout, setup)
    m = [ds.m[1, i] for i in cells]
    T = [ds.T[i] for i in cells]
    T_receiver = ds.T[receiver]
    m_res = sum(m)
    T̄ = sum(m .* T) / m_res
    R0 = sum(m .* k.(T))
    Ra = m_res * k(T̄)
    St.derive!(ds, tail.post, layout, setup)
    T_m = ds.T[receiver]
    m_post = ds.m[1, receiver]
    Rc = m_post * k(T_m)
    # The mass gasified from the event on: the trapezoid from the rate just
    # after the merge to the first saved time after the event, then over the
    # saved times; the residue, which nothing else adds to.
    mlr = sol.extras.gas_generation_rate
    j = something(findfirst(>=(tail.t), sol.t), length(sol.t))
    burst =
        0.5 * (Rc / area + mlr[j]) * (sol.t[j] - tail.t) + sum(
            0.5 * (mlr[i] + mlr[i+1]) * (sol.t[i+1] - sol.t[i]) for i = j:(length(sol.t)-1);
            init = 0.0,
        )
    after = maximum(mlr[j:end])
    t_after = sol.t[j-1+argmax(mlr[j:end])]
    old_mlr = old.extras.gas_generation_rate
    jo = findfirst(>=(tail.t), old.t)
    old_after = jo === nothing ? NaN : maximum(old_mlr[jo:min(jo+10, end)])
    g = 1000 / area
    return (;
        t = tail.t,
        cells = length(events),
        m_res = 1000 * m_res / area,
        m_post = 1000 * m_post / area,
        T̄,
        T_receiver,
        T_m,
        weight = (T_m - T̄) / (T_receiver - T̄),
        R0 = g * R0,
        Ra = g * Ra,
        Rc = g * Rc,
        burst = 1000 * burst,
        after = 1000 * after,
        t_after,
        old_after = 1000 * old_after,
    )
end

# ----------------------------------------------------------------------------
# The runs
# ----------------------------------------------------------------------------

# Item 1: exactly one invocation per layer, the largest merged fraction of
# any of them at most 5 ε_m.
item1_verdict(invocations, largest, ε_m; layers = 1) =
    invocations == layers && largest <= 5 * ε_m ? "pass" : "fail"

# The tables' cells, shared by the runs and the item 1 recompute.
_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
_tau(r) = "$(Int(r.τ)) ($(_g(r.τ / r.n^2)))"
_run(r) =
    "$(r.case) | $(r.n) | $(r.form) | $(r.ε_m) | $(_tau(r)) | " *
    (r.corrected ? "β = $(_g(r.monitor[2]))" : "none")

function measure(row, form, old)
    case = CASES[row.case]
    r = instrumented_solve(row, form)
    sol, p = r.sol, r.probe
    qois = Base.invokelatest(case.qois, sol)
    plain, wall = timed(() -> plain_solve(row, form))
    st = sol.statistics
    inv = Base.invokelatest(_vbo_invocations, sol, row.ε_m)
    item1 = item1_verdict(length(inv), maximum(v.largest for v in inv; init = 0.0), row.ε_m)
    failing = [q.name for q in qois if !(q.error < q.tol)]
    h0 = minimum(r.cp.def.setup.reference.thickness)
    tail =
        row.case == "V-BO-01" && !isempty(p.tails) ?
        tail_terms(r.cp, sol, p.tails[1], old.sol) : nothing
    m = measures(row.case, sol)
    return (;
        row...,
        form,
        retcode = string(sol.retcode),
        outcome = sol.outcome,
        t_last = sol.t[end],
        invocations = inv,
        events = length(sol.events),
        item1,
        qois,
        failing,
        measures = m,
        inadmissible = r.inadmissible,
        min_h = p.min_h,
        min_ratio = p.min_ratio,
        item10 = p.min_ratio >= 0.5 * row.ε_m && sol.retcode != :DtLessThanMin ? "pass" :
                 "fail",
        h0,
        c52 = p.c52,
        phase1 = p.phase1,
        phase1_above = p.phase1_above,
        mismatch = p.mismatch,
        package_largest = st.largest_constraint_residual,
        domain = st.domain_rejections,
        steps = st.accepted,
        rejected = st.rejected,
        probed = p.steps,
        same_steps = plain.statistics.accepted == st.accepted &&
                     plain.statistics.rejected == st.rejected,
        wall,
        old_wall = old.wall,
        tail,
        substrate = row.case == "V-BO-08" ?
                    (; substrate_history(sol)..., p.times, p.first_active) : nothing,
    )
end

# A line of the csv split into its fields' ranges, a quoted field whole.
function _csv_ranges(line)
    ranges, start, quoted = UnitRange{Int}[], 1, false
    for (i, c) in pairs(line)
        if c == '"'
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(ranges, start:prevind(line, i))
            start = i + 1
        end
    end
    push!(ranges, start:lastindex(line))
    return ranges
end

_replace_range(line, r, new) =
    line[firstindex(line):prevind(line, first(r))] * new * line[nextind(line, last(r)):end]

# Recompute item 1 from the csv and write it into the csv's item 1 column and
# the item 1 cells of the md's outcome table, matched on the run's columns.
function recompute_item1()
    csv, md = joinpath(OUT, "burnout_ladder.csv"), joinpath(OUT, "burnout_ladder.md")
    lines = readlines(csv)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    verdicts = Dict{String,String}()
    for k = 2:length(lines)
        ranges = _csv_ranges(lines[k])
        field(name) = lines[k][ranges[col(name)]]
        num(name) = parse(Float64, field(name))
        verdict = item1_verdict(
            parse(Int, field("invocations")),
            num("largest_fraction"),
            num("eps_m"),
        )
        lines[k] = _replace_range(lines[k], ranges[col("item1")], verdict)
        row = (;
            case = field("case"),
            n = parse(Int, field("cells")),
            form = field("form"),
            ε_m = num("eps_m"),
            τ = num("tau"),
            corrected = field("corrected") == "true",
            monitor = (1.0, num("beta")),
        )
        verdicts[_run(row)] = verdict
    end
    open(io -> foreach(l -> println(io, l), lines), csv, "w")
    text = readlines(md)
    first_row = findfirst(startswith("| V-BO-"), text)
    last_row = findnext(!startswith("| "), text, first_row) - 1
    for k = first_row:last_row
        cells = split(text[k], " | ")
        key = join(cells[1:6], " | ")[3:end]
        cells[9] = verdicts[key]
        text[k] = join(cells, " | ")
    end
    stamp = findfirst(startswith("Produced by"), text)
    text[stamp] = replace(
        text[stamp],
        r"(at commit [^.]*)\." => SubstitutionString(
            "\\1; item 1 recomputed from this table's csv at $(COMMIT_STAMP).",
        );
        count = 1,
    )
    text .= replace.(
        text,
        "with at most one invocation per layer" => "with exactly one invocation per layer",
    )
    open(io -> foreach(l -> println(io, l), text), md, "w")
    return nothing
end

if "--item1" in ARGS
    recompute_item1()
    println("recomputed item 1 in $(OUT)")
    exit()
end

load_start = Sys.loadavg()
started = time()
old_runs = Dict{Tuple{String,Int},Any}()
for row in ROWS
    key = (row.case, row.n)
    haskey(old_runs, key) && continue
    sol, wall = timed(() -> old_solve(row))
    old_runs[key] = (;
        sol,
        wall,
        measures = measures(row.case, sol),
        parks = length(sol.extras.park_times),
    )
    println("old core $(key): $(sol.retcode) $(round(wall, digits = 3)) s")
    flush(stdout)
end

results = Any[]
for form in FORMS, row in ROWS
    res = measure(row, form, old_runs[(row.case, row.n)])
    push!(results, res)
    println(
        "$(row.case) n=$(row.n) ε_m=$(row.ε_m) τ=$(row.τ) ",
        "corrected=$(row.corrected) β=$(row.monitor[2]) $(form): $(res.retcode) $(res.outcome) ",
        "invocations=$(length(res.invocations)) item1=$(res.item1) ",
        "failing=$(res.failing) item10=$(res.item10) wall=$(round(res.wall, digits = 3)) s",
    )
    flush(stdout)
end
load_end = Sys.loadavg()

# V-BO-08's substrate before the first activation (§6.2): the corrected run at
# ε_m = 1e-2 against the run without the correction, at the saved times up to
# the end of the last accepted step before the first whose state has an active
# node (the step that crosses the activation moves every saved time its
# interpolant spans).
substrate_comparison = map(FORMS) do form
    pick(corrected) = only(
        r for r in results if r.case == "V-BO-08" &&
            r.form == form &&
            r.corrected == corrected &&
            r.ε_m == 1e-2
    )
    a, b = pick(true).substrate, pick(false).substrate
    t_inert = maximum(filter(<(a.first_active), a.times); init = 0.0)
    k = min(something(findlast(<=(t_inert), a.t), 1), length(b.t))
    d = maximum(abs, view(a.T, :, 1:k) .- view(b.T, :, 1:k))
    # The accepted steps up to the last compared saved time, and the first
    # at which the two solves part.
    ta, tb = filter(<=(a.t[k]), a.times), filter(<=(a.t[k]), b.times)
    parted = findfirst(j -> ta[j] != tb[j], 1:min(length(ta), length(tb)))
    (;
        form,
        t_active = a.first_active,
        t_inert,
        k,
        abs = d,
        rel = d / maximum(view(a.T, :, 1:k)),
        steps = (length(ta), length(tb)),
        parted = parted === nothing ? NaN : ta[parted],
    )
end

# ----------------------------------------------------------------------------
# The tables
# ----------------------------------------------------------------------------

_invocations(r) =
    isempty(r.invocations) ? "0" :
    join(
        (
            "$(v.cells) ($(v.below) below, $(v.band) band; $(_g(100 * v.largest)) %)" for
            v in r.invocations
        ),
        "; ",
    )

open(joinpath(OUT, "burnout_ladder.csv"), "w") do io
    println(
        io,
        "case,cells,form,eps_m,tau,corrected,beta,retcode,outcome,t_last,invocations,events,",
        "largest_fraction,item1,qois_failing,peak,t_peak,t_bo,mass_loss,ledger,",
        "inadmissible,min_h,min_ratio,item10,c52,phase1,phase1_above,domain,steps,rejected,wall,old_wall",
    )
    for r in results
        m = r.measures
        vals = (
            r.case,
            r.n,
            r.form,
            r.ε_m,
            r.τ,
            r.corrected,
            r.monitor[2],
            r.retcode,
            r.outcome,
            r.t_last,
            length(r.invocations),
            r.events,
            maximum(v.largest for v in r.invocations; init = 0.0),
            r.item1,
            "\"$(join(r.failing, "; "))\"",
            m === nothing ? NaN : m.peak,
            m === nothing ? NaN : m.t_peak,
            m === nothing ? NaN : m.t_bo,
            m === nothing ? NaN : m.mass_loss,
            m === nothing ? NaN : m.ledger,
            r.inadmissible,
            r.min_h,
            r.min_ratio,
            r.item10,
            r.c52,
            r.phase1,
            r.phase1_above,
            r.domain,
            r.steps,
            r.rejected,
            r.wall,
            r.old_wall,
        )
        println(io, join(string.(vals), ","))
    end
end

open(joinpath(OUT, "burnout_ladder.md"), "w") do io
    println(io, "# Gate items 1, 2, 4, 6 and 10: the burnout ladder\n")
    println(
        io,
        "Produced by `test/newcore/gates/burnout_ladder.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "The finite-volume core with the equidistribution correction (band $(Ig.CORRECTION_BAND), ",
        "activation 0.3/0.1, the monitor's mass term `β` per row), the layer interfaces Lagrangian; ",
        "τ [s] is the layer's ",
        "equidistribution time, with the prototype's per-cell `τ_p = τ/N²` [s] in parentheses. ",
        "The old core runs each case at the same resolution with its own keywords (its threshold ",
        "is a thickness fraction).\n",
    )
    println(io, "## Outcomes and QoIs (items 1, 2 and 4)\n")
    println(
        io,
        "Invocations: per invocation its cells (below ε_m, taken by the band; the largest merged ",
        "mass fraction). Item 1 passes with exactly one invocation per layer and the largest ",
        "fraction at most 5 ε_m. QoIs: the case's own, in assert mode against its fixed references; ",
        "the failing ones are named. Inadmissible: results of a deactivation the solve refused.\n",
    )
    println(
        io,
        "| case | cells | form | ε_m | τ (τ_p) | correction | retcode, outcome, end [s] | invocations | item 1 | QoIs | failing | peak, t_peak, t_bo (old core) | events | inadmissible |",
    )
    println(io, "|" * repeat("---|", 14))
    for r in results
        m = r.measures
        o = old_runs[(r.case, r.n)].measures
        pk =
            m === nothing ? "" :
            "$(_g(m.peak; d = 5)), $(_g(m.t_peak; d = 4)), $(_g(m.t_bo; d = 4)) " *
            "($(_g(o.peak; d = 5)), $(_g(o.t_peak; d = 4)), $(_g(o.t_bo; d = 4)))"
        passed = count(q -> q.error < q.tol, r.qois)
        println(
            io,
            "| $(_run(r)) | $(r.retcode), $(r.outcome), $(_g(r.t_last; d = 4)) | $(_invocations(r)) | ",
            "$(r.item1) | $(passed)/$(length(r.qois)) | $(join(r.failing, "; ")) | $(pk) | ",
            "$(r.events) | $(r.inadmissible) |",
        )
    end
    println(io, "\n## Cells, residuals and cost (items 6 and 10)\n")
    println(
        io,
        "Smallest cell: the thinnest live cell at the end of any step, before its deactivation ",
        "event [m], and the smallest ratio of a live cell to its initial thickness over ε_m ",
        "(item 10: at least 0.5, and no `DtLessThanMin`). On a step whose event is located, the ",
        "step's end is read before the event is located, a state after the event's time that the ",
        "trajectory does not pass through: the thinnest cell there is thinner than at the event ",
        "(item 10's reading is conservative), and the residuals are those of a proposal the domain ",
        "test accepted rather than of a state on the trajectory. C52 and Phase 1: the largest scaled ",
        "constraint residual of the accepted states under C52's scale and under Phase 1's, which ",
        "omits the sweep of the node corrections; Phase 1 > 1: the accepted states Phase 1's scale ",
        "would have rejected. Domain: the domain test's rejections (C52's scale, admissibility ",
        "included). Wall: the warm median of $(REPEATS) solves without the recorders [s], the old ",
        "core's at the same resolution in parentheses; same: whether those solves took the ",
        "instrumented solve's steps and rejections (in the enthalpy form the recorders' residual ",
        "evaluations move the warm start of the caloric inversion, so the two solves part at ",
        "round-off).\n",
    )
    println(
        io,
        "| case | cells | form | ε_m | τ (τ_p) | correction | smallest cell [m] | ratio / ε_m | item 10 | C52 | Phase 1 | Phase 1 > 1 | domain | steps (rejected) | wall (old core) | same |",
    )
    println(io, "|" * repeat("---|", 16))
    for r in results
        println(
            io,
            "| $(_run(r)) | $(_g(r.min_h)) | $(_g(r.min_ratio / r.ε_m)) | $(r.item10) | ",
            "$(_g(r.c52)) | $(_g(r.phase1)) | $(r.phase1_above) | $(r.domain) | ",
            "$(r.steps) ($(r.rejected)) | $(_g(r.wall)) ($(_g(r.old_wall))) | $(r.same_steps) |",
        )
    end
    worst = maximum(r -> r.mismatch, results)
    println(
        io,
        "\nThe script's C52 scale against `Integration.constraint_residual` at every probed ",
        "state: largest relative difference $(_g(worst)).\n",
    )
    println(io, "## R25: the end of the layer (V-BO-01)\n")
    println(
        io,
        "At the state that ends the step of the layer's invocation, before the merge, and after ",
        "it: the cells merged, the residue `m_res` [g/m²] (after the merge, in the receiver), the ",
        "mass-weighted temperature `T̄` of the merged PMMA [K], the receiver's temperature before ",
        "the merge and the merged temperature `T_m` [K], and the receiver's weight in it, ",
        "`(T_m − T̄)/(T_receiver − T̄)`; the PMMA gasification rate [g m⁻² s⁻¹] `R0 = Σ mᵢ k(Tᵢ)` ",
        "before the merge, `Ra = m_res k(T̄)` lumped and `Rc = m_res k(T_m)` after the merge; the ",
        "mass gasified from the event to the end of the run [g/m²] (from `Rc` at the ",
        "event to the first saved time after it, then the saved times); the largest mass loss rate at those saved times and its time, the old core's ",
        "largest over the ten saved times from the event's.\n",
    )
    println(
        io,
        "| cells | form | ε_m | τ (τ_p) | event [s] | merged | m_res (after) | T̄ | T receiver | T_m | weight | R0 | Ra | Rc | Rc/Ra | Rc/R0 | gasified after | largest after (t) | old core |",
    )
    println(io, "|" * repeat("---|", 19))
    for r in results
        r.tail === nothing && continue
        x = r.tail
        println(
            io,
            "| $(r.n) | $(r.form) | $(r.ε_m) | $(_tau(r)) | $(_g(x.t; d = 5)) | ",
            "$(x.cells) | $(_g(x.m_res)) ($(_g(x.m_post))) | $(_g(x.T̄; d = 5)) | ",
            "$(_g(x.T_receiver; d = 5)) | $(_g(x.T_m; d = 5)) | $(_g(x.weight)) | $(_g(x.R0)) | ",
            "$(_g(x.Ra)) | $(_g(x.Rc)) | $(_g(x.Rc / x.Ra)) | $(_g(x.Rc / x.R0)) | $(_g(x.burst)) | $(_g(x.after; d = 4)) ($(_g(x.t_after; d = 4))) | ",
            "$(_g(x.old_after; d = 4)) |",
        )
    end
    println(io, "\n## V-BO-08: the substrate before the first activation\n")
    println(
        io,
        "The temperatures of the cells that start without PMMA at the saved times up to the end ",
        "of the last accepted step before the first whose state has an active node (the step ",
        "that crosses the activation moves the saved times its interpolant spans): the ",
        "corrected run at ε_m = 1e-2 ",
        "against the run without the correction, the largest difference and that over the ",
        "largest temperature (§6.2: at most 1e-10 on the same accepted steps); the accepted ",
        "steps of the two solves up to the last saved time compared, and the time of the first ",
        "step at which they part (a changed sequence is reported with the difference, as item 5 ",
        "does). The run without the correction is this comparison's reference only: its layer ",
        "ends `:interior`, and it is not a run of items 1 and 2.\n",
    )
    println(
        io,
        "| form | first active step [s] | compared to [s] | saved times | largest difference [K] | relative | steps (without) | parted at [s] | verdict |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for x in substrate_comparison
        verdict = !isnan(x.parted) ? "steps differ" : x.rel <= 1e-10 ? "pass" : "fail"
        println(
            io,
            "| $(x.form) | $(_g(x.t_active; d = 5)) | $(_g(x.t_inert; d = 5)) | $(x.k) | ",
            "$(_g(x.abs)) | $(_g(x.rel)) | ",
            "$(x.steps[1]) ($(x.steps[2])) | $(_g(x.parted; d = 5)) | $(verdict) |",
        )
    end
    println(io, "\n## The old core\n")
    println(io, "| case | cells | retcode | parks | wall [s] |")
    println(io, "|---|---|---|---|---|")
    for ((id, n), o) in sort(collect(old_runs); by = first)
        println(io, "| $(id) | $(n) | $(o.sol.retcode) | $(o.parks) | $(_g(o.wall)) |")
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
