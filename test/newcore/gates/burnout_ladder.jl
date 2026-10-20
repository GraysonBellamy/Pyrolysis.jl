# Gate items 1, 2, 4, 6 and 10 of the Phase 2 plan (§7.3) and R25's tail
# diagnostic, as M4 measures them: the burnout cases on the finite-volume core
# with the equidistribution correction, across grids, thresholds, `τ` and the
# reactivity switch.
#
# The matrix, every run in both energy forms, the band at its default (0.1):
#
#   V-BO-01 at 30, 60 and 120 PMMA cells × ε_m ∈ {3e-3, 1e-2, 3e-2}, the
#   default τ;
#   V-BO-01 at 60 cells and ε_m = 1e-2 with τ = 360 and 3600 s (τ is the
#   layer's equidistribution time; the prototype's per-cell τ_p = τ/N² is
#   0.1, 0.2 and 1 s at 60 cells for τ = 360, 720 and 3600 s);
#   V-BO-01 at 30, 60 and 120 cells, ε_m = 1e-2, without the reactivity switch
#   (`switch_width = 0`);
#   V-BO-03 at 30, 60 and 120 cells at the defaults;
#   V-BO-04 at its 40 cells without and with the correction.
#
# Per run: the outcome; the deactivation invocations, each with its cells
# below ε_m, the cells the band took and its largest merged mass fraction
# (item 1: at most one invocation per layer, the largest fraction at most
# 5 ε_m); every QoI of the case in assert mode against its fixed references,
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
# step before any discrete callback runs, and never crosses zero.
#
# R25's tail diagnostic, for every V-BO-01 run: at the state that ends the
# step of the layer's invocation, before the merge, the PMMA mass `mᵢ`, the
# temperature `Tᵢ` and the switch `sᵢ` of the merged cells and the receiver;
# after it the residue and the merged temperature `T_m`; the rate
# `Σ mᵢ k(Tᵢ) sᵢ → m_res k(T̄) s̄ → m_res k(T̄) → m_res k(T_m) s_m`, `T̄` and
# `s̄` the mass-weighted means (lumping, the switch's release, the merged
# temperature); the mass gasified from the event to the end of the run
# against the residue; the largest mass loss rate after the event.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_ladder.jl
# Writes test/newcore/gates/measurements/burnout_ladder.{csv,md}.

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
const SWITCH = Ke.SWITCH_WIDTH
const CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-03", "V-BO-04")
)

const ROWS = [
    [
        (; case = "V-BO-01", n, ε_m, τ = TAU, switch = SWITCH, corrected = true) for
        n in (30, 60, 120) for ε_m in (3e-3, 1e-2, 3e-2)
    ]
    [
        (; case = "V-BO-01", n = 60, ε_m = 1e-2, τ, switch = SWITCH, corrected = true)
        for τ in (360.0, 3600.0)
    ]
    [
        (; case = "V-BO-01", n, ε_m = 1e-2, τ = TAU, switch = 0.0, corrected = true) for
        n in (30, 60, 120)
    ]
    [
        (; case = "V-BO-03", n, ε_m = 1e-2, τ = TAU, switch = SWITCH, corrected = true)
        for n in (30, 60, 120)
    ]
    [
        (; case = "V-BO-04", n = 40, ε_m = 1e-2, τ = TAU, switch = SWITCH, corrected)
        for corrected in (false, true)
    ]
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
    stash::Vector{Float64}
    stash_t::Float64
    events::Int
    tails::Vector{Any}
end
Probe(n) = Probe(Inf, Inf, 0.0, 0.0, 0, 0.0, 0, zeros(n), NaN, 0, Any[])

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
    end
    p.steps += 1
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

# The new core's keywords of a run beside the case's own.
function run_keywords(row, form, problem)
    corr =
        row.corrected ?
        (;
            correction = Cl.EquidistributionCorrection(τ = row.τ),
            lagrangian = row.case == "V-BO-04" ? Int[] :
                         [Base.invokelatest(_vbo_interface_node, problem)],
        ) : (;)
    return (; energy_state = form, threshold = row.ε_m, switch_width = row.switch, corr...)
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
    problem = Base.invokelatest(case.build, row.n)
    kw = run_keywords(row, form, problem)
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
    problem = Base.invokelatest(case.build, row.n)
    cp = Ig.setup(problem; core_keywords(case)..., run_keywords(row, form, problem)...)
    return Ig.integrate(cp; min_thickness = case.solve_kwargs.min_thickness)
end

old_solve(case, n) = Logging.with_logger(Logging.NullLogger()) do
    solve(Base.invokelatest(case.build, n); progress = false, case.solve_kwargs...)
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
    reaction = Ke.find_kernel(Ke.ReactionKernel, def.kernels)
    sw = reaction.switch
    s_of(ds, i) = sw === nothing ? 1.0 : Ke.reactivity(sw, ds, i)
    area = setup.reference.cross_section
    ds = ws.derived
    St.derive!(ds, tail.pre, layout, setup)
    m = [ds.m[1, i] for i in cells]
    T = [ds.T[i] for i in cells]
    s = [s_of(ds, i) for i in cells]
    T_receiver = ds.T[receiver]
    m_res = sum(m)
    T̄ = sum(m .* T) / m_res
    s̄ = sum(m .* s) / m_res
    R0 = sum(m .* k.(T) .* s)
    Ra = m_res * k(T̄) * s̄
    Rb = m_res * k(T̄)
    St.derive!(ds, tail.post, layout, setup)
    T_m = ds.T[receiver]
    m_post = ds.m[1, receiver]
    s_m = s_of(ds, receiver)
    Rc = m_post * k(T_m) * s_m
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
        s̄,
        T_receiver,
        T_m,
        weight = (T_m - T̄) / (T_receiver - T̄),
        s_m,
        R0 = g * R0,
        Ra = g * Ra,
        Rb = g * Rb,
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

function measure(row, form, old)
    case = CASES[row.case]
    r = instrumented_solve(row, form)
    sol, p = r.sol, r.probe
    qois = Base.invokelatest(case.qois, sol)
    plain, wall = timed(() -> plain_solve(row, form))
    st = sol.statistics
    inv = Base.invokelatest(_vbo_invocations, sol, row.ε_m)
    layers = 1
    item1 =
        length(inv) <= layers && all(v -> v.largest <= 5 * row.ε_m, inv) ? "pass" : "fail"
    failing = [q.name for q in qois if !(q.error < q.tol)]
    h0 = minimum(r.cp.def.setup.reference.thickness)
    tail =
        row.case == "V-BO-01" && !isempty(p.tails) ?
        tail_terms(r.cp, sol, p.tails[1], old.sol) : nothing
    m = row.case == "V-BO-04" ? nothing : Base.invokelatest(_vbo_measures, sol)
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
    )
end

load_start = Sys.loadavg()
started = time()
old_runs = Dict{Tuple{String,Int},Any}()
for row in ROWS
    key = (row.case, row.n)
    haskey(old_runs, key) && continue
    case = CASES[row.case]
    sol, wall = timed(() -> old_solve(case, row.n))
    old_runs[key] = (;
        sol,
        wall,
        measures = row.case == "V-BO-04" ? nothing : Base.invokelatest(_vbo_measures, sol),
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
        "$(row.case) n=$(row.n) ε_m=$(row.ε_m) τ=$(row.τ) switch=$(row.switch) ",
        "corrected=$(row.corrected) $(form): $(res.retcode) $(res.outcome) ",
        "invocations=$(length(res.invocations)) item1=$(res.item1) ",
        "failing=$(res.failing) item10=$(res.item10) wall=$(round(res.wall, digits = 3)) s",
    )
    flush(stdout)
end
load_end = Sys.loadavg()

# ----------------------------------------------------------------------------
# The tables
# ----------------------------------------------------------------------------

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
_tau(r) = "$(Int(r.τ)) ($(_g(r.τ / r.n^2)))"
_invocations(r) =
    isempty(r.invocations) ? "0" :
    join(
        (
            "$(v.cells) ($(v.below) below, $(v.band) band; $(_g(100 * v.largest)) %)" for
            v in r.invocations
        ),
        "; ",
    )
_run(r) =
    "$(r.case) | $(r.n) | $(r.form) | $(r.ε_m) | $(_tau(r)) | $(r.switch) | " *
    (r.corrected ? "yes" : "no")

open(joinpath(OUT, "burnout_ladder.csv"), "w") do io
    println(
        io,
        "case,cells,form,eps_m,tau,switch,corrected,retcode,outcome,t_last,invocations,events,",
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
            r.switch,
            r.corrected,
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
    println(io, "# Gate items 1, 2, 4, 6 and 10 at M4: the burnout ladder\n")
    println(
        io,
        "Produced by `test/newcore/gates/burnout_ladder.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) at the end. ",
        "The finite-volume core with the equidistribution correction (band $(Ig.CORRECTION_BAND), ",
        "activation 0.3/0.1, `β = 0`), the PMMA–Kaowool interface Lagrangian; τ [s] is the layer's ",
        "equidistribution time, with the prototype's per-cell `τ_p = τ/N²` [s] in parentheses; ",
        "switch: the reactivity switch's width in ε_m (0: none). The old core runs each case at the ",
        "same resolution with its own keywords (its threshold is a thickness fraction).\n",
    )
    println(io, "## Outcomes and QoIs (items 1, 2 and 4)\n")
    println(
        io,
        "Invocations: per invocation its cells (below ε_m, taken by the band; the largest merged ",
        "mass fraction). Item 1 passes with at most one invocation per layer and the largest ",
        "fraction at most 5 ε_m. QoIs: the case's own, in assert mode against its fixed references; ",
        "the failing ones are named. Inadmissible: results of a deactivation the solve refused.\n",
    )
    println(
        io,
        "| case | cells | form | ε_m | τ (τ_p) | switch | corrected | retcode, outcome, end [s] | invocations | item 1 | QoIs | failing | peak, t_peak, t_bo (old core) | events | inadmissible |",
    )
    println(io, "|" * repeat("---|", 15))
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
        "(item 10: at least 0.5, and no `DtLessThanMin`). C52 and Phase 1: the largest scaled ",
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
        "| case | cells | form | ε_m | τ (τ_p) | switch | corrected | smallest cell [m] | ratio / ε_m | item 10 | C52 | Phase 1 | Phase 1 > 1 | domain | steps (rejected) | wall (old core) | same |",
    )
    println(io, "|" * repeat("---|", 17))
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
        "mass-weighted temperature `T̄` and switch `s̄` of the merged PMMA [K], the receiver's ",
        "temperature before the merge and the merged temperature `T_m` [K] with its switch, and the ",
        "receiver's weight in it, `(T_m − T̄)/(T_receiver − T̄)`; the ",
        "PMMA gasification rate [g m⁻² s⁻¹] `R0 = Σ mᵢ k(Tᵢ) sᵢ` before the merge, `Ra = m_res k(T̄) s̄` ",
        "lumped, `Rb = m_res k(T̄)` with the switch released and `Rc = m_res k(T_m) s_m` after the ",
        "merge; the mass gasified from the event to the end of the run [g/m²] (from `Rc` at the ",
        "event to the first saved time after it, then the saved times); the largest mass loss rate at those saved times and its time, the old core's ",
        "largest over the ten saved times from the event's.\n",
    )
    println(
        io,
        "| cells | form | ε_m | τ (τ_p) | switch | event [s] | merged | m_res (after) | T̄ | s̄ | T receiver | T_m (s_m) | weight | R0 | Ra | Rb | Rc | Rb/Ra | Rc/Rb | Rc/R0 | gasified after | largest after (t) | old core |",
    )
    println(io, "|" * repeat("---|", 23))
    for r in results
        r.tail === nothing && continue
        x = r.tail
        println(
            io,
            "| $(r.n) | $(r.form) | $(r.ε_m) | $(_tau(r)) | $(r.switch) | $(_g(x.t; d = 5)) | ",
            "$(x.cells) | $(_g(x.m_res)) ($(_g(x.m_post))) | $(_g(x.T̄; d = 5)) | $(_g(x.s̄)) | ",
            "$(_g(x.T_receiver; d = 5)) | $(_g(x.T_m; d = 5)) ($(_g(x.s_m))) | $(_g(x.weight)) | ",
            "$(_g(x.R0)) | ",
            "$(_g(x.Ra)) | $(_g(x.Rb)) | $(_g(x.Rc)) | $(_g(x.Rb / x.Ra)) | $(_g(x.Rc / x.Rb)) | ",
            "$(_g(x.Rc / x.R0)) | $(_g(x.burst)) | $(_g(x.after; d = 4)) ($(_g(x.t_after; d = 4))) | ",
            "$(_g(x.old_after; d = 4)) |",
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
