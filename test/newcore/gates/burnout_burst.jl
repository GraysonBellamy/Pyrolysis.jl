# Gate item 2's burst criteria (Phase 2 plan §7.3): the end-of-layer burst of
# V-BO-01 judged by refinement on both cores, from the gas-generation rate at
# every accepted step and every post-event state.
#
# Runs: the finite-volume core through the case's solve at 30, 60, 120 and
# 240 PMMA cells, ε_m ∈ {1e-2, 3e-3, 3e-2}, in both energy forms; the old core
# at the same cells with the case's keywords. The criteria, over the window
# [t_pk, t_end] with t_pk the old core's 240-cell peak time (its peak QoI's):
#
#   (a) the largest per-step excess of the N-cell rate over the 2N-cell rate
#       after t_pk falls at least 1.5× per grid doubling and is no larger than
#       the old core's at the same N;
#   (b) the L1 of successive grids falls at least 1.5× per doubling, and the
#       L1 against the same family's 240-cell run is within 1.5× of the old
#       core's at each N;
#   (c) the L1 against the old core's 240-cell run falls with refinement.
#
# They gate the finite-volume core at ε_m of 1e-2 and 3e-3 in both forms;
# 3e-2, whose layer takes two invocations, is reported. A series is
# piecewise linear between its points, with a jump at every event (the state
# before the event and the state after it at the same time), and the
# integrals are exact for it. A located event moves the end of its step back
# to the crossing: the recorder then keeps the step's interpolated state at
# the crossing, before the event, and the state after it. In the enthalpy form
# the recorder's residual evaluations move the warm start of the caloric
# inversion, so its solves part from plain ones at round-off.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_burst.jl
# Writes test/newcore/gates/measurements/burnout_burst.{csv,md}.

ENV["PYROLYSIS_CORE"] = "new"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using OrdinaryDiffEq: ContinuousCallback, DiscreteCallback, CallbackSet, SciMLBase

const As = Pyrolysis.Assembly
const De = Pyrolysis.Deactivation

const OUT = joinpath(@__DIR__, "measurements")
const CASE =
    only(filter(c -> c.id == "V-BO-01", Base.invokelatest(load_verification_cases)))
const CELLS = (30, 60, 120, 240)
const THRESHOLDS = (1e-2, 3e-3, 3e-2)
const GATED = (1e-2, 3e-3)
const FORMS = (:temperature, :enthalpy)
const FALL = 1.5
const SELF = 1.5

# ----------------------------------------------------------------------------
# The series
# ----------------------------------------------------------------------------

# The gas-generation rate of a solve [kg m⁻² s⁻¹]: at the end of every
# accepted step (kind 0) and after every event (kind 1, at the event's time).
struct Series
    t::Vector{Float64}
    rate::Vector{Float64}
    kind::Vector{Int8}
end
Series() = Series(Float64[], Float64[], Int8[])

function record!(s::Series, t, rate, kind)
    push!(s.t, t)
    push!(s.rate, rate)
    push!(s.kind, kind)
    return nothing
end

# The series in time order, the state before an event ahead of the state after it.
function ordered(s::Series)
    p = sortperm(collect(zip(s.t, s.kind)))
    return Series(s.t[p], s.rate[p], s.kind[p])
end

# The finite-volume core's rate at `u`, as `extras.gas_generation_rate`: the
# gas the live cells' reactions make, over the initial cross section.
function core_rate(cp, du, u, t)
    def, ws = cp.def, cp.ws
    As.residual!(du, u, def, t, ws)
    ds, c = As.buffers(ws, def, Float64)
    gas = map(Pyrolysis.Materials.is_gas, def.setup.material.components)
    rate = 0.0
    for i in eachindex(ds.active)
        ds.active[i] || continue
        for j in eachindex(gas)
            gas[j] && (rate += c.S_mass[j, i])
        end
    end
    return rate / def.setup.reference.cross_section
end

# The recorder of a set-up problem. A continuous callback whose condition never
# crosses zero reads the state that ends each step before any event, and, when
# the threshold margin crosses inside the step, the step's interpolated state
# at the first crossing; a discrete callback, which runs after the step's
# events, keeps the step's end, or the state at the crossing when a located
# event moved the step's end back to it, and the state after any event.
mutable struct Pending
    on::Bool
    t::Float64
    tprev::Float64
    rate::Float64
    u::Vector{Float64}
    crossing::Bool
    t_cross::Float64
    rate_cross::Float64
end
Pending() = Pending(false, NaN, NaN, NaN, Float64[], false, NaN, NaN)

function recorder(cp, s::Series)
    du = similar(cp.u0)
    cells, layout = cp.deactivation, cp.def.layout
    margin(u) = De.threshold_margin(cells, layout, u)
    p = Pending()
    reader = ContinuousCallback(
        (u, t, integ) -> begin
            if t == integ.t && !(p.on && integ.tprev == p.tprev)
                p.on, p.t, p.tprev = true, t, integ.tprev
                p.rate = core_rate(cp, du, u, t)
                p.u = copy(u)
                p.crossing = false
                if margin(u) <= 0 && margin(integ.uprev) > 0
                    a, b = integ.tprev, t
                    for _ = 1:60
                        mid = (a + b) / 2
                        margin(integ(mid)) > 0 ? (a = mid) : (b = mid)
                    end
                    p.crossing, p.t_cross = true, b
                    p.rate_cross = core_rate(cp, du, integ(b), b)
                end
            end
            return 1.0
        end,
        integ -> nothing;
        save_positions = (false, false),
        interp_points = 0,
    )
    commit = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            if p.on && integ.t == p.t
                record!(s, p.t, p.rate, 0)
                integ.u != p.u &&
                    record!(s, integ.t, core_rate(cp, du, integ.u, integ.t), 1)
            elseif p.on && integ.t < p.t
                p.crossing && record!(s, p.t_cross, p.rate_cross, 0)
                record!(s, integ.t, core_rate(cp, du, integ.u, integ.t), 1)
            end
            p.on = false
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
    return CallbackSet(reader, commit)
end

function new_run(n, ε_m, form)
    problem = Base.invokelatest(CASE.build, n)
    s = Series()
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(
            CASE.solve_fn,
            problem;
            CASE.solve_kwargs...,
            energy_state = form,
            threshold = ε_m,
            instrument = cp -> recorder(cp, s),
        )
    end
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    return (;
        core = :new,
        n,
        ε_m,
        form,
        retcode = sol.retcode,
        steps = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        invocations = length(inv),
        last_event = isempty(sol.events) ? Inf : maximum(e.t for e in sol.events),
        series = ordered(s),
        Base.invokelatest(_vbo_measures, sol)[(:peak, :t_peak)]...,
    )
end

# The old core: the states that end its steps and those its steps start from
# when a park changed them, passed afterwards through its own extraction
# (a post-park state at its park's time plus 1e-9 s, so the park counts).
function old_run(n)
    problem = Base.invokelatest(CASE.build, n)
    ts, us, kinds = Float64[], Vector{Float64}[], Int8[]
    last_u, last_t = Ref(Float64[]), Ref(NaN)
    reader = ContinuousCallback(
        (u, t, integ) -> begin
            if t == integ.t && t != last_t[]
                if !isempty(last_u[]) && integ.tprev == last_t[] && integ.uprev != last_u[]
                    push!(ts, integ.tprev)
                    push!(us, copy(integ.uprev))
                    push!(kinds, Int8(1))
                end
                push!(ts, t)
                push!(us, copy(u))
                push!(kinds, Int8(0))
                last_u[], last_t[] = copy(u), t
            end
            return 1.0
        end,
        integ -> nothing;
        save_positions = (false, false),
        interp_points = 0,
    )
    sol = Logging.with_logger(Logging.NullLogger()) do
        solve(problem; progress = false, callback = reader, CASE.solve_kwargs...)
    end
    _, ws, _ = Logging.with_logger(Logging.NullLogger()) do
        Pyrolysis.initialize_state(
            Base.invokelatest(CASE.build, n);
            use_ale = true,
            radiation_model = BEER_LAMBERT,
        )
    end
    pt = copy(sol.extras.park_times)
    dep = Pyrolysis.Solver.DepletionCallbackState(ws.layout, 0, Float64[], pt, Int[])
    ex = Pyrolysis.Problem.PyrolysisSolution(Base.invokelatest(CASE.build, n))
    raw = (;
        t = [k == 1 ? t + 1e-9 : t for (t, k) in zip(ts, kinds)],
        u = us,
        retcode = :Success,
    )
    Logging.with_logger(Logging.NullLogger()) do
        Pyrolysis.Solver.extract_solution_ale!(
            ex,
            raw,
            ex.problem,
            ws.layout;
            depletion = dep,
            radiation_model = BEER_LAMBERT,
        )
    end
    st = sol.raw_solution.stats
    return (;
        core = :old,
        n,
        ε_m = NaN,
        form = :old,
        retcode = sol.retcode,
        steps = st.naccept,
        rejected = st.nreject,
        invocations = length(pt),
        last_event = isempty(pt) ? Inf : maximum(pt),
        series = ordered(Series(ts, copy(ex.extras.gas_generation_rate), kinds)),
        Base.invokelatest(_vbo_measures, sol)[(:peak, :t_peak)]...,
    )
end

# ----------------------------------------------------------------------------
# The measures of piecewise-linear series with jumps
# ----------------------------------------------------------------------------

# The segments between successive points of positive length; the two points
# of an event share a time and bound none.
struct Piecewise
    a::Vector{Float64}
    ya::Vector{Float64}
    b::Vector{Float64}
    yb::Vector{Float64}
end

function Piecewise(s::Series)
    keep = findall(isfinite, s.rate)
    t, y = s.t[keep], s.rate[keep]
    seg = [k for k = 1:(length(t)-1) if t[k+1] > t[k]]
    return Piecewise(t[seg], y[seg], t[seg .+ 1], y[seg .+ 1])
end

# The value at `x` of the segment starting at or before it (after a jump, the
# value after the event); `NaN` outside the series.
function value(f::Piecewise, x)
    i = searchsortedlast(f.a, x)
    (i == 0 || x > f.b[i]) && return NaN
    return f.ya[i] + (f.yb[i] - f.ya[i]) * (x - f.a[i]) / (f.b[i] - f.a[i])
end

breakpoints(f::Piecewise) = unique!(sort!(vcat(f.a, f.b)))

# ∫_lo^hi |f − g| dt, exact for piecewise-linear series: each interval between
# the union of their breakpoints is linear in both.
function l1(f::Piecewise, g::Piecewise, lo, hi)
    xs = filter(x -> lo < x < hi, vcat(breakpoints(f), breakpoints(g)))
    xs = unique!(sort!(vcat(lo, xs, hi)))
    total = 0.0
    for k = 1:(length(xs)-1)
        x0, x1 = xs[k], xs[k+1]
        h = x1 - x0
        δ = 1e-9h
        da = value(f, x0 + δ) - value(g, x0 + δ)
        db = value(f, x1 - δ) - value(g, x1 - δ)
        (isfinite(da) && isfinite(db)) || continue
        total +=
            da * db >= 0 ? 0.5 * (abs(da) + abs(db)) * h :
            0.5 * (da^2 + db^2) / (abs(da) + abs(db)) * h
    end
    return total
end

# The largest excess of `s` over `g` at the points of `s` after `lo`, within `g`.
function excess(s::Series, g::Piecewise, lo)
    hi = last(g.b)
    best = -Inf
    for (t, y) in zip(s.t, s.rate)
        (t > lo && t < hi && isfinite(y)) || continue
        best = max(best, y - value(g, t))
    end
    return best
end

# The largest per-step rate before the solve's last event, and its time (on
# the old core a park before its last can raise one above the main peak).
function main_peak(r)
    s = r.series
    k = findall(i -> s.t[i] < r.last_event && isfinite(s.rate[i]), eachindex(s.t))
    i = k[argmax(s.rate[k])]
    return s.rate[i], s.t[i]
end

ratios(x) = [x[k] / x[k+1] for k = 1:(length(x)-1)]

# The criteria's measures of one family (the runs of one core, threshold and
# form at the four resolutions) against the old core's runs.
function family(runs, old, t_pk, t_end)
    f = [Piecewise(r.series) for r in runs]
    fo = [Piecewise(r.series) for r in old]
    exc = [excess(runs[k].series, f[k+1], t_pk) for k = 1:(length(runs)-1)]
    exc_old = [excess(old[k].series, fo[k+1], t_pk) for k = 1:(length(old)-1)]
    succ = [l1(f[k], f[k+1], t_pk, t_end) for k = 1:(length(f)-1)]
    succ_old = [l1(fo[k], fo[k+1], t_pk, t_end) for k = 1:(length(fo)-1)]
    self = [l1(f[k], f[end], t_pk, t_end) for k = 1:(length(f)-1)]
    self_old = [l1(fo[k], fo[end], t_pk, t_end) for k = 1:(length(fo)-1)]
    cross = [l1(f[k], fo[end], t_pk, t_end) for k = 1:length(f)]
    a = all(>=(FALL), ratios(exc)) && all(exc .<= exc_old)
    b = all(>=(FALL), ratios(succ)) && all(self ./ self_old .<= SELF)
    c = all(<(0), diff(cross))
    return (;
        exc,
        exc_old,
        succ,
        succ_old,
        self,
        self_old,
        cross,
        a,
        b,
        c,
        pass = a && b && c,
    )
end

# ----------------------------------------------------------------------------
# The runs and the tables
# ----------------------------------------------------------------------------

function main()
    load_start = Sys.loadavg()
    started = time()
    old = Any[]
    for n in CELLS
        push!(old, old_run(n))
        r = old[end]
        println("old n=$(n): $(r.retcode) steps=$(r.steps) parks=$(r.invocations)")
        flush(stdout)
    end
    new = Dict{Any,Any}()
    for form in FORMS, ε_m in THRESHOLDS, n in CELLS
        r = new_run(n, ε_m, form)
        new[(form, ε_m, n)] = r
        println(
            "new $(form) ε_m=$(ε_m) n=$(n): $(r.retcode) steps=$(r.steps) invocations=$(r.invocations)",
        )
        flush(stdout)
    end
    load_end = Sys.loadavg()
    t_pk = old[end].t_peak
    t_end = Base.invokelatest(CASE.build, first(CELLS)).tspan[2]
    families = Any[]
    for form in FORMS, ε_m in THRESHOLDS
        m = family([new[(form, ε_m, n)] for n in CELLS], old, t_pk, t_end)
        push!(families, (; form, ε_m, gated = ε_m in GATED, m...))
    end
    write_tables(old, new, families, t_pk, t_end, load_start, load_end)
    gate = all(f -> f.pass, filter(f -> f.gated, families))
    println(
        "wrote $(OUT); the burst criteria ",
        gate ? "pass" : "fail",
        "; load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
    )
    return (; old, new, families)
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isfinite(x) ? string(round(x; sigdigits = d)) : string(x)) :
    string(x)
_list(x; scale = 1000) = join(_g.(scale .* x), ", ")
_yes(b) = b ? "pass" : "**fail**"

function write_tables(old, new, families, t_pk, t_end, load_start, load_end)
    mkpath(OUT)
    open(joinpath(OUT, "burnout_burst.csv"), "w") do io
        println(
            io,
            "core,form,eps_m,cells,retcode,steps,rejected,invocations,largest_before_last_event,t_largest,peak,t_peak,last_event",
        )
        for r in vcat(old, [new[k] for k in sort(collect(keys(new)); by = string)])
            pk, tp = main_peak(r)
            println(
                io,
                join(
                    (
                        r.core,
                        r.form,
                        r.ε_m,
                        r.n,
                        r.retcode,
                        r.steps,
                        r.rejected,
                        r.invocations,
                        1000pk,
                        tp,
                        r.peak,
                        r.t_peak,
                        r.last_event,
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(OUT, "burnout_burst.md"), "w") do io
        println(io, "# Gate item 2: the end-of-layer burst by refinement\n")
        println(
            io,
            "Produced by `test/newcore/gates/burnout_burst.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
            "$(Threads.nthreads()) thread; load average $(round.(load_start; digits = 1)) at the start ",
            "and $(round.(load_end; digits = 1)) at the end. V-BO-01, the finite-volume core through ",
            "the case's solve and the old core with the case's keywords, the gas-generation rate at every ",
            "accepted step and post-event state [g m⁻² s⁻¹]; window [$(_g(t_pk; d = 5)), $(_g(t_end; d = 4))] s ",
            "from the old core's 240-cell peak time. (a) The largest per-step excess of N cells over 2N ",
            "after the window's start, falling at least $(FALL)× per doubling and no larger than the old ",
            "core's; (b) the L1 of successive grids [g/m²], falling at least $(FALL)× per doubling, and ",
            "the L1 against the family's 240-cell run within $(SELF)× of the old core's; (c) the L1 ",
            "against the old core's 240-cell run falling with refinement. Gated: ε_m of ",
            "$(join(GATED, " and ")) in both forms.\n",
        )
        println(io, "## The criteria\n")
        println(
            io,
            "Each measure lists 30→60, 60→120 and 120→240 (the L1 against the old core's 240-cell run: ",
            "30, 60, 120 and 240), then its ratios per doubling.\n",
        )
        println(
            io,
            "| family | (a) excess of N over 2N; ratios | (a) | (b) successive L1; ratios | (b) self-L1 over the old core's | (b) | (c) L1 against the old core's 240 | (c) | gated |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|")
        f = first(families)
        println(
            io,
            "| old core | $(_list(f.exc_old)); $(_list(ratios(f.exc_old); scale = 1)) | | ",
            "$(_list(f.succ_old)); $(_list(ratios(f.succ_old); scale = 1)) | | | | | |",
        )
        for f in families
            println(
                io,
                "| $(f.form), ε_m = $(f.ε_m) | $(_list(f.exc)); $(_list(ratios(f.exc); scale = 1)) | $(_yes(f.a)) | ",
                "$(_list(f.succ)); $(_list(ratios(f.succ); scale = 1)) | $(_list(f.self ./ f.self_old; scale = 1)) | ",
                "$(_yes(f.b)) | $(_list(f.cross)) | $(_yes(f.c)) | $(f.gated ? "yes" : "reported") |",
            )
        end
        println(io, "\n## The runs\n")
        println(
            io,
            "| core | form | ε_m | cells | retcode | steps (rejected) | invocations or parks | largest per-step rate before the last event (t) | peak QoI (t) | last event [s] |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for r in vcat(
            old,
            [new[(form, ε_m, n)] for form in FORMS for ε_m in THRESHOLDS for n in CELLS],
        )
            pk, tp = main_peak(r)
            println(
                io,
                "| $(r.core) | $(r.form) | $(r.ε_m) | $(r.n) | $(r.retcode) | $(r.steps) ($(r.rejected)) | ",
                "$(r.invocations) | $(_g(1000pk; d = 5)) ($(_g(tp; d = 5))) | $(_g(r.peak; d = 5)) ($(_g(r.t_peak; d = 5))) | ",
                "$(_g(r.last_event; d = 6)) |",
            )
        end
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
