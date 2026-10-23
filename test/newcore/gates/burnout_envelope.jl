# W9 of the Phase 3 plan (§4) and its gate, §7 item 12(2): the burnout
# path's accuracy envelope under joint refinement, the resolutions,
# thresholds, bands, tolerances and activations inside which the hard
# criterion holds, with the deactivation-invocation count reported.
#
# The decks: V-BO-01, V-BO-06 (graded), V-BO-07 (mixed), V-BO-08 (thin) and
# the restricted deck at s = 1 (`tau_ladder.jl`'s `restricted_problem`: no
# re-radiation, an adiabatic back), at 30, 60, 120 and 240 PMMA cells
# (V-BO-08 at 10, 20, 40 and 80). The pairs (ε_m, band): (3e-3, 9e-3),
# (3e-3, 3e-2), (3e-3, 0.1), (1e-2, 3e-2) and (1e-2, 0.1), the bands 3 ε_m,
# 10 ε_m and 0.1 (at 1e-2 the second is the third). reltol 1e-6 and 1e-7,
# abstol the case's. The activation (r_hi, r_lo) of the correction (0.3, 0.1),
# the default, and (0.7, 0.3). The layer time is the default,
# `ActivationTime()`, throughout.
#
# The runs: the core matrix, every deck, resolution, pair, tolerance and
# activation in the temperature form; the same at each deck's CI resolution
# in the enthalpy form; two comparators at the default activation in the
# temperature form: the threshold-only removal (band = ε_m under the
# correction, at both thresholds, every deck and resolution) and the
# innermost receiver (`receiver = :innermost`, V-BO-07 at every resolution and
# V-BO-01 at 30 cells, every pair). A case deck solves through the case's own
# solve and its QoIs are the case's; the restricted deck's are the burnout,
# mass-loss and ledger QoIs, which need no reference (its adiabatic back
# reports no converged balance, so the surface QoI is left out). The
# invocation QoI is reported apart from the others.
#
# Per run: the QoIs, the failing named; the accepted-step peak, the largest
# per-step gas-generation rate over the run, and the main peak, the largest
# before the solve's last event (the peak QoI's, on the saved times); the
# integrated tail, the gas generated from the main peak's time to the end;
# the mass ledger over the accepted steps; the invocations, the cells each
# took and the largest merged fraction; the layer time τ and when it was set;
# the outcome, the steps and the wall time.
#
# Per family (a deck, pair, tolerance, activation and comparator at its four
# resolutions), the hard criterion: every run's QoIs (the invocation QoI
# aside); (a) the largest per-step excess of N cells over 2N after the finest
# run's main peak, and (b) the L1 of successive grids over that window, each
# falling at least 1.5× per doubling (an excess at or below zero counts as
# fallen); the accepted-step peak converging (successive differences not
# growing). The criterion is read over the four resolutions and, where it
# fails there, over the finest three: a family that passes from 60 cells (20
# on V-BO-08) is inside the envelope from that resolution.
#
# Run from the repository root, in one process:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_envelope.jl
# or in four, one per part, then merged:
#   julia ... burnout_envelope.jl part K DIR     (K = 1 to 4)
#   julia ... burnout_envelope.jl merge DIR
# Writes test/newcore/gates/measurements/burnout_envelope.{csv,md}.

include(joinpath(@__DIR__, "tau_ladder.jl"))

using Serialization

const DECK_NAMES = ("V-BO-01", "V-BO-06", "V-BO-07", "V-BO-08", "restricted")
const DECK_CELLS = Dict(
    "V-BO-01" => (30, 60, 120, 240),
    "V-BO-06" => (30, 60, 120, 240),
    "V-BO-07" => (30, 60, 120, 240),
    "V-BO-08" => (10, 20, 40, 80),
    "restricted" => (30, 60, 120, 240),
)
const PAIRS = ((3e-3, 9e-3), (3e-3, 3e-2), (3e-3, 0.1), (1e-2, 3e-2), (1e-2, 0.1))
const THRESHOLDS = (3e-3, 1e-2)
const RELTOLS = (1e-6, 1e-7)
const ACTIVATIONS = ((0.3, 0.1), (0.7, 0.3))
const DEFAULT_ACTIVATION = first(ACTIVATIONS)
const PARTS = (("V-BO-01",), ("restricted",), ("V-BO-07",), ("V-BO-06", "V-BO-08"))

# A deck: its builder at `n` PMMA cells, its solve, its QoIs and the key its
# solve records its step statistics under.
function deck(name)
    name == "restricted" && return (;
        name,
        build = n -> restricted_problem(n, 1.0),
        ci = 60,
        keywords = VBO_SOLVE_KWARGS,
        solve_fn = _vbo_solve_fn(name),
        qois = variant_qois(name; surface = false),
    )
    c = LADDER_CASES[name]
    return (;
        name,
        build = c.build,
        ci = c.ci.n_cells,
        keywords = c.solve_kwargs,
        solve_fn = c.solve_fn,
        qois = c.qois,
    )
end
const DECKS = Dict(name => deck(name) for name in DECK_NAMES)

# ----------------------------------------------------------------------------
# The matrix
# ----------------------------------------------------------------------------

# A run's key: (deck, cells, ε_m, band, reltol, activation, form, receiver).
# A family's: the key without the cells.
family_keys(name) = begin
    keys = Any[]
    for (ε_m, band) in PAIRS, reltol in RELTOLS, act in ACTIVATIONS
        push!(keys, (name, ε_m, band, reltol, act, :temperature, :beneath))
    end
    for ε_m in THRESHOLDS, reltol in RELTOLS
        push!(keys, (name, ε_m, ε_m, reltol, DEFAULT_ACTIVATION, :temperature, :beneath))
    end
    if name == "V-BO-07"
        for (ε_m, band) in PAIRS, reltol in RELTOLS
            push!(keys, (name, ε_m, band, reltol, DEFAULT_ACTIVATION, :temperature, :innermost))
        end
    end
    keys
end

function run_keys(name)
    keys = Any[]
    for f in family_keys(name), n in DECK_CELLS[name]
        push!(keys, (f[1], n, f[2:end]...))
    end
    ci = DECKS[name].ci
    for (ε_m, band) in PAIRS, reltol in RELTOLS, act in ACTIVATIONS
        push!(keys, (name, ci, ε_m, band, reltol, act, :enthalpy, :beneath))
    end
    if name == "V-BO-01"
        for (ε_m, band) in PAIRS, reltol in RELTOLS
            push!(
                keys,
                (name, 30, ε_m, band, reltol, DEFAULT_ACTIVATION, :temperature, :innermost),
            )
        end
    end
    return keys
end

# ----------------------------------------------------------------------------
# The runs
# ----------------------------------------------------------------------------

# ∫_lo^hi of a piecewise-linear series, exact.
function integral(f::Piecewise, lo, hi)
    total = 0.0
    for k in eachindex(f.a)
        a, b = max(f.a[k], lo), min(f.b[k], hi)
        b > a || continue
        total += 0.5 * (value(f, a) + value(f, b)) * (b - a)
    end
    return total
end

_isinvocation(q) = startswith(q.name, "deactivation invocations")

function envelope_run(key)
    name, n, ε_m, band, reltol, act, form, receiver = key
    d = DECKS[name]
    problem = Base.invokelatest(d.build, n)
    series = Series()
    bound = Ref{Any}(nothing)
    instrument = cp -> begin
        bound[] = cp
        recorder(cp, series)
    end
    started = time()
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(
            d.solve_fn,
            problem;
            d.keywords...,
            energy_state = form,
            threshold = ε_m,
            band,
            reltol,
            receiver,
            correction = Cl.EquidistributionCorrection(; activation = act),
            instrument,
        )
    end
    wall = time() - started
    corr = bound[].def.closure.correction
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    qois = Base.invokelatest(d.qois, sol)
    counted = filter(!_isinvocation, qois)
    failing = [q.name for q in counted if !(q.error < q.tol)]
    ordered_series = ordered(series)
    last_event = isempty(sol.events) ? Inf : maximum(e.t for e in sol.events)
    r = (;
        key,
        deck = name,
        n,
        ε_m,
        band,
        reltol,
        activation = act,
        form,
        receiver,
        τ = isempty(corr.time) ? NaN : first(corr.time),
        set_at = isempty(corr.set_at) ? NaN : first(corr.set_at),
        retcode = sol.retcode,
        outcome = sol.outcome,
        t_end = sol.t[end],
        invocations = inv,
        largest = isempty(inv) ? NaN : maximum(v.largest for v in inv),
        qois_total = length(counted),
        qois_failing = failing,
        qoi_pass = isempty(failing),
        ledger = VBO_STEP_STATS[name].ledger,
        steps = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        wall,
        series = ordered_series,
        last_event,
    )
    finite = findall(isfinite, ordered_series.rate)
    if length(finite) < 2
        return (; r..., peak = NaN, t_peak = NaN, main = NaN, t_main = NaN, tail = NaN)
    end
    k = finite[argmax(ordered_series.rate[finite])]
    before = filter(i -> ordered_series.t[i] < last_event, finite)
    main, t_main =
        isempty(before) ? (NaN, NaN) : begin
            j = before[argmax(ordered_series.rate[before])]
            ordered_series.rate[j], ordered_series.t[j]
        end
    tail = isnan(t_main) ? NaN : integral(Piecewise(ordered_series), t_main, sol.t[end])
    return (;
        r...,
        peak = 1000 * ordered_series.rate[k],
        t_peak = ordered_series.t[k],
        main = 1000main,
        t_main,
        tail = 1000tail,
    )
end

# ----------------------------------------------------------------------------
# The criteria
# ----------------------------------------------------------------------------

# The hard criterion on a family's runs, coarse to fine: (a), (b), the
# accepted-step peak's convergence and the QoIs; `NaN` measures and a failed
# criterion where a run gives no series.
function hard_criterion(runs)
    measures = (; exc = Float64[], succ = Float64[], peaks = [r.peak for r in runs])
    usable = all(r -> isfinite(r.main) && length(r.series.t) > 2, runs)
    qois = all(r -> r.qoi_pass, runs)
    usable || return (; measures..., a = false, b = false, peak = false, qois, hard = false)
    finest = runs[end]
    lo, hi = finest.t_main, finest.t_end
    f = [Piecewise(r.series) for r in runs]
    exc = 1000 .* [excess(runs[k].series, f[k+1], lo) for k = 1:(length(runs)-1)]
    succ = 1000 .* [l1(f[k], f[k+1], lo, hi) for k = 1:(length(f)-1)]
    steps = abs.(diff(measures.peaks))
    a = falls(exc)
    b = falls(succ)
    peak = all(steps[k+1] <= steps[k] for k = 1:(length(steps)-1))
    return (;
        exc,
        succ,
        peaks = measures.peaks,
        a,
        b,
        peak,
        qois,
        hard = a && b && peak && qois,
    )
end

# A family's runs and its criterion over the four resolutions and the finest
# three.
function envelope_family(fkey, runs)
    name = fkey[1]
    rs = [runs[(name, n, fkey[2:end]...)] for n in DECK_CELLS[name]]
    all4 = hard_criterion(rs)
    fine3 = all4.hard ? all4 : hard_criterion(rs[2:end])
    return (;
        key = fkey,
        deck = name,
        ε_m = fkey[2],
        band = fkey[3],
        reltol = fkey[4],
        activation = fkey[5],
        receiver = fkey[7],
        threshold_only = fkey[3] == fkey[2],
        all4,
        fine3,
        tails = [r.tail for r in rs],
        mains = [r.main for r in rs],
        invocations = [length(r.invocations) for r in rs],
        largest = [r.largest for r in rs],
        outcomes = [r.outcome for r in rs],
        failing = unique(vcat((r.qois_failing for r in rs)...)),
        from = all4.hard ? DECK_CELLS[name][1] : fine3.hard ? DECK_CELLS[name][2] : 0,
    )
end

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

function measure(names)
    runs = Dict{Any,Any}()
    for name in names, key in run_keys(name)
        r = envelope_run(key)
        runs[key] = r
        println(
            "$(key): τ=$(_g(r.τ)) $(r.retcode)/$(r.outcome) invocations=$(length(r.invocations)) ",
            "largest=$(_g(r.largest)) steps=$(r.steps) qois=$(r.qois_total - length(r.qois_failing))/",
            "$(r.qois_total) peak=$(_g(r.peak; d = 5)) main=$(_g(r.main; d = 5)) wall=$(_g(r.wall)) s",
        )
        flush(stdout)
    end
    families = [envelope_family(f, runs) for name in names for f in family_keys(name)]
    # The series are the families' only; the tables need the scalars.
    slim = Dict(k => Base.structdiff(r, NamedTuple{(:series,)}) for (k, r) in runs)
    return slim, families
end

function envelope_main(args)
    if length(args) == 3 && args[1] == "part"
        k = parse(Int, args[2])
        load_start, started = Sys.loadavg(), time()
        runs, families = measure(PARTS[k])
        mkpath(args[3])
        serialize(
            joinpath(args[3], "part$(k).jls"),
            (;
                runs,
                families,
                stamp = COMMIT_STAMP,
                load = (load_start, Sys.loadavg()),
                elapsed = time() - started,
            ),
        )
        println("wrote part $(k) to $(args[3]); $(round(Int, time() - started)) s")
        return nothing
    end
    if length(args) == 2 && args[1] == "merge"
        parts = [deserialize(joinpath(args[2], "part$(k).jls")) for k in eachindex(PARTS)]
        all(p -> p.stamp == COMMIT_STAMP, parts) || error(
            "the parts were measured at $(unique(p.stamp for p in parts)), not at $(COMMIT_STAMP)",
        )
        write_envelope(
            merge((p.runs for p in parts)...),
            vcat((p.families for p in parts)...),
            [p.load[1][1] for p in parts],
            [p.load[2][1] for p in parts],
            maximum(p.elapsed for p in parts),
        )
        println("wrote $(OUT)/burnout_envelope.{csv,md}")
        return nothing
    end
    isempty(args) || error("usage: burnout_envelope.jl [part K DIR | merge DIR]")
    load_start, started = Sys.loadavg(), time()
    runs, families = measure(DECK_NAMES)
    write_envelope(runs, families, [load_start[1]], [Sys.loadavg()[1]], time() - started)
    println("wrote $(OUT)/burnout_envelope.{csv,md}; $(round(Int, time() - started)) s")
    return nothing
end

# ----------------------------------------------------------------------------
# The tables
# ----------------------------------------------------------------------------

_act(a) = "$(a[1])/$(a[2])"
_pair(ε_m, band) = band == ε_m ? "$(ε_m), threshold only" : "$(ε_m), band $(band)"
_cells(f) = DECK_CELLS[f.deck]
# The QoIs by short names.
const _SHORT = (
    "complete burnout" => "burnout",
    "mass loss" => "mass loss",
    "burnout time" => "burnout time",
    "peak MLR" => "peak MLR",
    "time of peak" => "peak time",
    "mass ledger" => "ledger",
    "clamped surface" => "surface",
)
function _short(name)
    k = findfirst(x -> startswith(name, first(x)), _SHORT)
    return k === nothing ? name : last(_SHORT[k])
end
# What fails over the four resolutions: the criterion's parts and the QoIs.
function _failing(f)
    h = f.all4
    parts = [x for (x, ok) in (("(a)", h.a), ("(b)", h.b), ("peak", h.peak)) if !ok]
    return join(vcat(parts, unique(_short.(f.failing))), ", ")
end
_verdict(f) =
    (
        f.from == _cells(f)[1] ? "pass" :
        f.from == _cells(f)[2] ? "from $(f.from)" : "**fail**: $(_failing(f))"
    ) * _ended(f)
_ended(f) = any(==(:interior), f.outcomes) ? ", `:interior`" : ""
_invs(f) = join(("$(i)" for i in f.invocations), "/")
_load(l) = join(round.(l; digits = 1), ", ")

function _runs_in_order(runs)
    order = Any[]
    for name in DECK_NAMES, key in run_keys(name)
        haskey(runs, key) && push!(order, key)
    end
    return order
end

# The loads are each process's one-minute load average at its start and end.
function write_envelope(runs, families, load_start, load_end, elapsed)
    mkpath(OUT)
    order = _runs_in_order(runs)
    open(joinpath(OUT, "burnout_envelope.csv"), "w") do io
        println(
            io,
            "deck,cells,eps_m,band,reltol,r_hi,r_lo,form,receiver,tau,set_at,retcode,outcome,t_end,",
            "invocations,cells_merged,largest_fraction,qois_passing,qois_total,failing,ledger,",
            "peak,t_peak,main_peak,t_main,tail,steps,rejected,wall",
        )
        for k in order
            r = runs[k]
            println(
                io,
                join(
                    (
                        r.deck,
                        r.n,
                        r.ε_m,
                        r.band,
                        r.reltol,
                        r.activation[1],
                        r.activation[2],
                        r.form,
                        r.receiver,
                        r.τ,
                        r.set_at,
                        r.retcode,
                        r.outcome,
                        r.t_end,
                        length(r.invocations),
                        sum((v.cells for v in r.invocations); init = 0),
                        r.largest,
                        r.qois_total - length(r.qois_failing),
                        r.qois_total,
                        "\"$(join(r.qois_failing, "; "))\"",
                        r.ledger,
                        r.peak,
                        r.t_peak,
                        r.main,
                        r.t_main,
                        r.tail,
                        r.steps,
                        r.rejected,
                        r.wall,
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(OUT, "burnout_envelope.md"), "w") do io
        println(io, "# The burnout accuracy envelope (W9, §7 item 12(2))\n")
        println(
            io,
            "Produced by `test/newcore/gates/burnout_envelope.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
            "$(Threads.nthreads()) thread per process, $(length(load_start)) process(es); one-minute load ",
            "average $(_load(load_start)) at the start and $(_load(load_end)) at the end; the longest ",
            "process $(round(Int, elapsed)) s. ",
            "$(length(runs)) runs, $(length(families)) families. The finite-volume core with the ",
            "equidistribution correction at its default layer time (`ActivationTime()`, λ = ",
            "$(Cl.ACTIVATION_TIME_FACTOR)); the script's header defines every deck, run and measure.\n",
        )
        println(io, "## The envelope\n")
        println(
            io,
            "The hard criterion per family (a deck, threshold, band, tolerance and activation at its four ",
            "resolutions, the default receiver): pass over every resolution, \"from N\" over the finest ",
            "three, or fail with the failing QoIs named; the invocations at each resolution in brackets. ",
            "Every QoI within tolerance; (a) the largest per-step excess of N over 2N and (b) the successive ",
            "L1, each falling at least $(FALL)× per doubling; the accepted-step peak converging.\n",
        )
        cols = vcat([(ε, b) for (ε, b) in PAIRS], [(ε, ε) for ε in THRESHOLDS])
        println(
            io,
            "| deck | activation | reltol | ",
            join((_pair(c...) for c in cols), " | "),
            " |",
        )
        println(io, "|---|---|---|", repeat("---|", length(cols)))
        index = Dict(
            (f.deck, f.ε_m, f.band, f.reltol, f.activation, f.receiver) => f for
            f in families
        )
        for name in DECK_NAMES, act in ACTIVATIONS, reltol in RELTOLS
            any(f -> f.deck == name, families) || continue
            cells = map(cols) do (ε, b)
                f = get(index, (name, ε, b, reltol, act, :beneath), nothing)
                f === nothing ? "" : "$(_verdict(f)) [$(_invs(f))]"
            end
            println(io, "| $(name) | $(_act(act)) | $(reltol) | $(join(cells, " | ")) |")
        end
        println(io, "\n## The innermost receiver\n")
        println(
            io,
            "V-BO-07's families under each receiver, and V-BO-01's 30-cell runs (D142), at the default ",
            "activation: the verdict, the invocations, the largest merged fraction [%], the accepted-step ",
            "and main peaks [g m⁻² s⁻¹] and the tail [g/m²] (at the finest resolution for a family).\n",
        )
        println(
            io,
            "| deck | ε_m, band | reltol | receiver | verdict | invocations | largest | peak | main | tail |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for (ε, b) in PAIRS, reltol in RELTOLS, rec in (:beneath, :innermost)
            f = get(index, ("V-BO-07", ε, b, reltol, DEFAULT_ACTIVATION, rec), nothing)
            f === nothing && continue
            println(
                io,
                "| V-BO-07 | $(_pair(ε, b)) | $(reltol) | $(rec) | $(_verdict(f)) | $(_invs(f)) | ",
                "$(join(_g.(100 .* f.largest), ", ")) | $(_g(f.all4.peaks[end]; d = 5)) | ",
                "$(_g(f.mains[end]; d = 5)) | $(_g(f.tails[end]; d = 4)) |",
            )
        end
        for (ε, b) in PAIRS, reltol in RELTOLS, rec in (:beneath, :innermost)
            r = get(
                runs,
                ("V-BO-01", 30, ε, b, reltol, DEFAULT_ACTIVATION, :temperature, rec),
                nothing,
            )
            r === nothing && continue
            verdict =
                r.qoi_pass ? "QoIs pass" :
                "**fail**: $(join(unique(_short.(r.qois_failing)), ", "))"
            println(
                io,
                "| V-BO-01, 30 cells | $(_pair(ε, b)) | $(reltol) | $(rec) | $(verdict) | ",
                "$(length(r.invocations)) | $(_g(100r.largest)) | $(_g(r.peak; d = 5)) | ",
                "$(_g(r.main; d = 5)) | $(_g(r.tail; d = 4)) |",
            )
        end
        println(io, "\n## The energy forms at the CI resolution\n")
        println(
            io,
            "Each deck at its CI resolution in the enthalpy form against the temperature form's run with ",
            "the same parameters: the QoIs passing of counted and the failing named, the invocations, the ",
            "accepted-step peak [g m⁻² s⁻¹].\n",
        )
        println(
            io,
            "| deck | ε_m, band | reltol | activation | temperature: QoIs, invocations, peak | enthalpy: QoIs, invocations, peak |",
        )
        println(io, "|---|---|---|---|---|---|")
        for name in DECK_NAMES, (ε, b) in PAIRS, reltol in RELTOLS, act in ACTIVATIONS
            ci = DECKS[name].ci
            haskey(runs, (name, ci, ε, b, reltol, act, :enthalpy, :beneath)) || continue
            entry(form) = begin
                r = runs[(name, ci, ε, b, reltol, act, form, :beneath)]
                q = "$(r.qois_total - length(r.qois_failing))/$(r.qois_total)"
                fail =
                    isempty(r.qois_failing) ? "" :
                    " ($(join(unique(_short.(r.qois_failing)), ", ")))"
                "$(q)$(fail), $(length(r.invocations)), $(_g(r.peak; d = 5))"
            end
            println(
                io,
                "| $(name) | $(_pair(ε, b)) | $(reltol) | $(_act(act)) | $(entry(:temperature)) | ",
                "$(entry(:enthalpy)) |",
            )
        end
        println(io, "\n## The families\n")
        println(
            io,
            "(a) [g m⁻² s⁻¹] and (b) [g/m²] per doubling, then their ratios; the accepted-step peaks ",
            "[g m⁻² s⁻¹]; the tails [g/m²]; the largest merged fractions [%]; each part of the criterion ",
            "over the four resolutions, and the verdict.\n",
        )
        println(
            io,
            "| deck | ε_m, band | reltol | activation | receiver | (a); ratios | (b); ratios | peaks | tails | largest | (a) | (b) | peak | QoIs | verdict |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for f in families
            h = f.all4
            println(
                io,
                "| $(f.deck) | $(_pair(f.ε_m, f.band)) | $(f.reltol) | $(_act(f.activation)) | ",
                "$(f.receiver) | $(_list(h.exc; scale = 1)); $(_list(ratios(h.exc); scale = 1)) | ",
                "$(_list(h.succ; scale = 1)); $(_list(ratios(h.succ); scale = 1)) | ",
                "$(join(_g.(h.peaks; d = 5), ", ")) | $(join(_g.(f.tails; d = 4), ", ")) | ",
                "$(join(_g.(100 .* f.largest), ", ")) | $(_yes(h.a)) | $(_yes(h.b)) | ",
                "$(_yes(h.peak)) | $(_yes(h.qois)) | $(_verdict(f)) |",
            )
        end
        println(io, "\n## The runs\n")
        println(
            io,
            "τ: the layer's time [s] (set at t [s]). Invocations: each one's cells (largest merged ",
            "fraction). Peak: the accepted-step peak (t); main: the largest before the last event (t) ",
            "[g m⁻² s⁻¹]; tail [g/m²]; ledger: the mass ledger's closure. Wall: on a loaded machine.\n",
        )
        println(
            io,
            "| deck | cells | ε_m | band | reltol | activation | form | receiver | τ (set at) | retcode, outcome, end | invocations | QoIs | failing | ledger | peak (t) | main (t) | tail | steps (rejected) | wall [s] |",
        )
        println(
            io,
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        )
        for k in order
            r = runs[k]
            println(
                io,
                "| $(r.deck) | $(r.n) | $(r.ε_m) | $(r.band) | $(r.reltol) | $(_act(r.activation)) | ",
                "$(r.form) | $(r.receiver) | $(_g(r.τ; d = 4)) ($(_g(r.set_at; d = 5))) | ",
                "$(r.retcode), $(r.outcome), $(_g(r.t_end; d = 5)) | $(_inv(r)) | ",
                "$(r.qois_total - length(r.qois_failing))/$(r.qois_total) | $(join(unique(_short.(r.qois_failing)), ", ")) | ",
                "$(_g(r.ledger)) | $(_g(r.peak; d = 5)) ($(_g(r.t_peak; d = 5))) | ",
                "$(_g(r.main; d = 5)) ($(_g(r.t_main; d = 5))) | $(_g(r.tail; d = 4)) | ",
                "$(r.steps) ($(r.rejected)) | $(_g(r.wall)) |",
            )
        end
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && envelope_main(ARGS)
