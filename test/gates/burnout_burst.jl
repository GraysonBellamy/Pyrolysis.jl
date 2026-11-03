# The end-of-layer burst judged by refinement, from the gas-generation rate
# at every accepted step and every post-event state.
#
# Runs: each deck through the case's solve at 30, 60, 120 and 240 PMMA cells,
# in both energy forms. V-BO-01 at ε_m ∈ {1e-2, 3e-3, 3e-2}; V-BO-03 run to
# `VBO03_T_EXTENDED`, past its layer's invocation, V-BO-06 and V-BO-07 at
# the default threshold, reported. The criteria, per family (a deck, form
# and threshold at the four resolutions) over the window [t_pk, t_end], with
# t_pk the peak time (its peak QoI's) of the deck's 240-cell run in the
# temperature form at ε_m = 1e-2:
#
#   (a) the largest per-step excess of the N-cell rate over the 2N-cell rate
#       after t_pk falls at least 1.5× per grid doubling;
#   (b) the L1 of successive grids falls at least 1.5× per doubling.
#
# The L1 of each resolution against the family's 240-cell run is reported
# beside. The criteria judge V-BO-01 at ε_m of 1e-2 and 3e-3 in both forms;
# its 3e-2, whose layer takes two invocations, and the other decks are
# reported. A series is piecewise linear between its points, with a jump at
# every event (the state before the event and the state after it at the same
# time), and the integrals are exact for it. A located event moves the end of
# its step back to the crossing: the recorder then keeps the step's
# interpolated state at the crossing, before the event, and the state after
# it. In the enthalpy form the recorder's residual evaluations move the warm
# start of the caloric inversion, so its solves part from plain ones at
# round-off.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/gates/burnout_burst.jl
# Writes test/gates/measurements/burnout_burst.{csv,md}.

ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using OrdinaryDiffEq: ContinuousCallback, DiscreteCallback, CallbackSet, SciMLBase

const As = Pyrolysis.Assembly
const De = Pyrolysis.Deactivation

const OUT = joinpath(@__DIR__, "measurements")
const CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-03", "V-BO-06", "V-BO-07")
)
const CELLS = (30, 60, 120, 240)
const JUDGED = (1e-2, 3e-3)
# The run whose peak time starts each deck's window: its form and threshold
# at the finest resolution.
const WINDOW_RUN = (:temperature, 1e-2)
# Each deck: its case, its builder at `n` PMMA cells and its thresholds.
const DECKS = [
    (; id = "V-BO-01", build = CASES["V-BO-01"].build, thresholds = (1e-2, 3e-3, 3e-2)),
    (;
        id = "V-BO-03",
        build = n -> _vbo_problem(n; substrate_alpha = 0.0, t_end = VBO03_T_EXTENDED),
        thresholds = (1e-2,),
    ),
    (; id = "V-BO-06", build = CASES["V-BO-06"].build, thresholds = (1e-2,)),
    (; id = "V-BO-07", build = CASES["V-BO-07"].build, thresholds = (1e-2,)),
]
judged(id, ε_m) = id == "V-BO-01" && ε_m in JUDGED
const FORMS = (:temperature, :enthalpy)
const FALL = 1.5

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

# The rate at `u`, as `extras.gas_generation_rate`: the gas the live cells'
# reactions make, over the initial cross section.
function rate_at(cp, du, u, t)
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
                p.rate = rate_at(cp, du, u, t)
                p.u = copy(u)
                p.crossing = false
                if margin(u) <= 0 && margin(integ.uprev) > 0
                    a, b = integ.tprev, t
                    for _ = 1:60
                        mid = (a + b) / 2
                        margin(integ(mid)) > 0 ? (a = mid) : (b = mid)
                    end
                    p.crossing, p.t_cross = true, b
                    p.rate_cross = rate_at(cp, du, integ(b), b)
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
                    record!(s, integ.t, rate_at(cp, du, integ.u, integ.t), 1)
            elseif p.on && integ.t < p.t
                p.crossing && record!(s, p.t_cross, p.rate_cross, 0)
                record!(s, integ.t, rate_at(cp, du, integ.u, integ.t), 1)
            end
            p.on = false
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
    return CallbackSet(reader, commit)
end

function burst_run(deck, n, ε_m, form)
    case = CASES[deck.id]
    problem = Base.invokelatest(deck.build, n)
    s = Series()
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(
            case.solve_fn,
            problem;
            case.solve_kwargs...,
            energy_state = form,
            threshold = ε_m,
            instrument = cp -> recorder(cp, s),
        )
    end
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    return (;
        id = deck.id,
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

# The largest per-step rate before the solve's last event, and its time.
function main_peak(r)
    s = r.series
    k = findall(i -> s.t[i] < r.last_event && isfinite(s.rate[i]), eachindex(s.t))
    i = k[argmax(s.rate[k])]
    return s.rate[i], s.t[i]
end

ratios(x) = [x[k] / x[k+1] for k = 1:(length(x)-1)]

# The criteria's measures of one family (the runs of one threshold and form
# at the four resolutions), and the L1 of each resolution against its finest.
function family(runs, t_pk, t_end)
    f = [Piecewise(r.series) for r in runs]
    exc = [excess(runs[k].series, f[k+1], t_pk) for k = 1:(length(runs)-1)]
    succ = [l1(f[k], f[k+1], t_pk, t_end) for k = 1:(length(f)-1)]
    self = [l1(f[k], f[end], t_pk, t_end) for k = 1:(length(f)-1)]
    a = all(>=(FALL), ratios(exc))
    b = all(>=(FALL), ratios(succ))
    return (; exc, succ, self, a, b, pass = a && b)
end

# ----------------------------------------------------------------------------
# The runs and the tables
# ----------------------------------------------------------------------------

function main()
    load_start = Sys.loadavg()
    started = time()
    runs = Dict{Any,Any}()
    windows = Dict{String,Any}()
    families = Any[]
    for deck in DECKS
        for form in FORMS, ε_m in deck.thresholds, n in CELLS
            r = burst_run(deck, n, ε_m, form)
            runs[(deck.id, form, ε_m, n)] = r
            println(
                "$(deck.id) $(form) ε_m=$(ε_m) n=$(n): $(r.retcode) steps=$(r.steps) invocations=$(r.invocations)",
            )
            flush(stdout)
        end
        t_pk = runs[(deck.id, WINDOW_RUN..., last(CELLS))].t_peak
        t_end = Base.invokelatest(deck.build, first(CELLS)).tspan[2]
        windows[deck.id] = (t_pk, t_end)
        for form in FORMS, ε_m in deck.thresholds
            m = family([runs[(deck.id, form, ε_m, n)] for n in CELLS], t_pk, t_end)
            push!(
                families,
                (; id = deck.id, form, ε_m, judged = judged(deck.id, ε_m), m...),
            )
        end
    end
    load_end = Sys.loadavg()
    write_tables(runs, families, windows, load_start, load_end)
    verdict = all(f -> f.pass, filter(f -> f.judged, families))
    println(
        "wrote $(OUT); the burst criteria ",
        verdict ? "pass" : "fail",
        "; load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
    )
    return (; runs, families)
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isfinite(x) ? string(round(x; sigdigits = d)) : string(x)) :
    string(x)
_list(x; scale = 1000) = join(_g.(scale .* x), ", ")
_yes(b) = b ? "pass" : "**fail**"

function write_tables(by_key, families, windows, load_start, load_end)
    mkpath(OUT)
    runs = [
        by_key[(deck.id, form, ε_m, n)] for deck in DECKS for form in FORMS for
        ε_m in deck.thresholds for n in CELLS
    ]
    open(joinpath(OUT, "burnout_burst.csv"), "w") do io
        println(
            io,
            "case,form,eps_m,cells,retcode,steps,rejected,invocations,largest_before_last_event,t_largest,peak,t_peak,last_event",
        )
        for r in runs
            pk, tp = main_peak(r)
            println(
                io,
                join(
                    (
                        r.id,
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
        println(io, "# The end-of-layer burst by refinement\n")
        println(
            io,
            "Produced by `test/gates/burnout_burst.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
            "$(Threads.nthreads()) thread; load average $(round.(load_start; digits = 1)) at the start ",
            "and $(round.(load_end; digits = 1)) at the end. Each deck through ",
            "the case's solve, the gas-generation rate at every ",
            "accepted step and post-event state [g m⁻² s⁻¹]; each deck's window from the peak time ",
            "of its 240-cell run in the $(WINDOW_RUN[1]) form at ε_m = $(WINDOW_RUN[2]) to its end (",
            join(
                (
                    "$(id) [$(_g(w[1]; d = 5)), $(_g(w[2]; d = 5))] s" for
                    (id, w) in sort(collect(windows))
                ),
                ", ",
            ),
            "; V-BO-03 is run to $(_g(VBO03_T_EXTENDED; d = 5)) s). (a) The largest per-step excess of N cells over 2N ",
            "after the window's start, falling at least $(FALL)× per doubling; (b) the L1 of successive ",
            "grids [g/m²], falling at least $(FALL)× per doubling. Reported beside: the L1 of each ",
            "resolution against the family's 240-cell run. Judged: V-BO-01 at ε_m of ",
            "$(join(JUDGED, " and ")) in both forms; the rest is reported.\n",
        )
        println(io, "## The criteria\n")
        println(
            io,
            "Each measure lists 30→60, 60→120 and 120→240 (the L1 against the family's 240-cell run: ",
            "30, 60 and 120), then its ratios per doubling.\n",
        )
        println(
            io,
            "| family | (a) excess of N over 2N; ratios | (a) | (b) successive L1; ratios | (b) | L1 against the family's 240; ratios | judged |",
        )
        println(io, "|---|---|---|---|---|---|---|")
        for f in families
            println(
                io,
                "| $(f.id), $(f.form), ε_m = $(f.ε_m) | $(_list(f.exc)); $(_list(ratios(f.exc); scale = 1)) | $(_yes(f.a)) | ",
                "$(_list(f.succ)); $(_list(ratios(f.succ); scale = 1)) | $(_yes(f.b)) | ",
                "$(_list(f.self)); $(_list(ratios(f.self); scale = 1)) | $(f.judged ? "yes" : "reported") |",
            )
        end
        println(io, "\n## The runs\n")
        println(
            io,
            "| case | form | ε_m | cells | retcode | steps (rejected) | invocations | largest per-step rate before the last event (t) | peak QoI (t) | last event [s] |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for r in runs
            pk, tp = main_peak(r)
            println(
                io,
                "| $(r.id) | $(r.form) | $(r.ε_m) | $(r.n) | $(r.retcode) | $(r.steps) ($(r.rejected)) | ",
                "$(r.invocations) | $(_g(1000pk; d = 5)) ($(_g(tp; d = 5))) | $(_g(r.peak; d = 5)) ($(_g(r.t_peak; d = 5))) | ",
                "$(_g(r.last_event; d = 6)) |",
            )
        end
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
