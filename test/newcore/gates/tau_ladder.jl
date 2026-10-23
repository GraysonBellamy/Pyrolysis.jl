# W3 of the Phase 3 plan (§4) and its gate, §7 item 11: the equidistribution
# time of a vanishing layer across thickness, composition, heating, kinetics,
# rescaled time, a heating change after the layer's activation and a long
# layer, under three rules:
#
#   fixed — one time for every layer, EQUIDISTRIBUTION_TIME;
#   A     — that time scaled by the layer's initial thickness against
#           V-BO-01's 5.9 mm, where it was calibrated;
#   B     — `ActivationTime()`: the time set at the layer's first activation
#           from its thickness and thinning rate there;
#
# and, reported beside them, B with λ = 0.5 (`B(0.5)`), the factor under
# which the decks heated hardest keep one invocation.
#
# The configurations: V-BO-01, V-BO-06 (graded), V-BO-07 (mixed) and V-BO-08
# (thin); V-BO-01's deck at 25 and 100 kW/m² and at half and twice its
# pre-exponential factor; the restricted deck (V-BO-01's material and mesh
# with a convective loss for the re-radiation and an adiabatic back) at s = 1,
# 2 and 1/2, every rate
# times s and every time over s; V-BO-01's flux stepped to half and to one and
# a half times its value at 150 s, after the layer's first activation; and
# the long layer (a 10 mm Darcy–Fick wood column with deactivation, 300 s).
# The restricted deck replaces the re-radiation by a convective loss, since
# the emitted flux does not scale with the rates.
#
# The runs, in the temperature form at ε_m = 1e-2 and reltol 1e-6 unless
# stated: every configuration at its three resolutions (V-BO-08 at 10, 20 and
# 40 cells, the long layer at 25 and 100) under each rule; V-BO-01, V-BO-06
# and V-BO-08 at ε_m = 3e-3 and 3e-2 at their CI resolutions (Phase 2 item
# 1's rows); every configuration at reltol 1e-7 at its CI resolution, and the
# restricted deck at s = 1 at every resolution (the solve's own error, against
# which the rescaled pairs are read); V-BO-01, 06, 07 and 08 in the enthalpy
# form at their CI resolutions.
#
# Per run: the layer's time and, under B, when it was set; the first
# activation t_a (the first accepted step on which a layer cell is below
# r_hi of its thickness, interpolated linearly) and the activated life from
# it to the layer's first invocation, with T_sync = τ/π² over that life; the
# outcome; the invocations, each with its cells and largest merged fraction,
# and Phase 2 item 1 (one invocation for the layer, the largest fraction at
# most 5 ε_m); the case's QoIs in assert mode, the invocation QoI apart
# (accuracy is the hard criterion and the count a reported target, plan
# question 21; a variant deck's QoIs are the burnout, mass-loss, ledger and
# surface ones, which need no reference, the restricted deck's without the
# surface one, since its adiabatic back reports no converged balance); the
# main peak over the accepted
# steps; the steps, the Jacobians and the wall time of the instrumented solve.
#
# Per family (a configuration's resolutions under one rule), W9's hard
# criterion: every run's QoIs (the invocation QoI aside); (a) the largest
# per-step excess of N cells over 2N after the finest run's main peak, and
# (b) the L1 of successive grids over that window, each falling at least
# 1.5× per doubling (an excess at or below zero counts as fallen); the main
# peak converging (successive differences not growing).
#
# The rescaled pairs: the restricted deck at s = 2 and 1/2 against s = 1 at
# the same resolution and rule, in dimensionless time: the largest saved
# temperature and node-position differences and the event-time shift,
# against five times the s = 1 run's own differences between reltol 1e-6
# and 1e-7.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/tau_ladder.jl
# Writes test/newcore/gates/measurements/tau_ladder.{csv,md}.

include(joinpath(@__DIR__, "burnout_burst.jl"))

using Statistics: median

const Ig = Pyrolysis.Integration
const Cl = Pyrolysis.Closures

const LADDER_CASES = Dict(
    c.id => c for c in Base.invokelatest(load_verification_cases) if
    c.id in ("V-BO-01", "V-BO-06", "V-BO-07", "V-BO-08")
)
const RULES = (:fixed, :A, :B, :B05)
const H_REF = VBO_L_PMMA
const PAIR_BOUND = 5

_smoothstep(x) = x <= 0 ? 0.0 : x >= 1 ? 1.0 : x * x * (3 - 2x)

# ----------------------------------------------------------------------------
# The decks beyond the cases
# ----------------------------------------------------------------------------

# V-BO-01's material and mesh at `n` PMMA cells with every rate times `s`:
# the incident flux, the convective coefficient, the conductivities, the
# gas-transfer coefficients and the pre-exponential factor; no re-radiation,
# the top's convective coefficient 30 W/(m²·K) in its place (the re-radiation
# of a surface near 650 K, linearized), which keeps the bare substrate inside
# the property bracket after the layer is gone; an adiabatic back; the times
# over `s`.
function restricted_problem(n, s; t_end = 500.0)
    k_pmma(T) = T <= 394.0 ? 0.16 : 0.2156 * (T / 300.0)^-0.745
    c_pmma(T) = T <= 394.0 ? 1127.0793 * (T / 300.0)^1.9283 : 1768.9067 * (T / 300.0)^0.5552
    pmma = SolidComponent(
        :PMMA,
        ρ = VBO_RHO_PMMA,
        c = CallableProperty(c_pmma),
        k = CallableProperty(T -> s * k_pmma(T)),
        ε = 0.96,
        α = 2870.0 / VBO_RHO_PMMA,
        λ = s * 2e-5,
    )
    kaowool = SolidComponent(
        :KAOWOOL,
        ρ = VBO_RHO_KW,
        c = 1100.0,
        k = s * 0.058,
        ε = 1.0,
        λ = s * 2e-5,
        α = 8.0e6 / VBO_RHO_KW,
    )
    gas = GaseousComponent(:MMA, M = 0.044, c = 0.0, k = s * 0.02, λ = s * 2e-2)
    rx =
        Reaction(:pmma_gasification, 1 => 3; A = s * 1.35e11, E = 164e3, h = 8.7e5, n = 1.0)
    material = Material(
        name = :PMMA_restricted,
        components = (pmma, kaowool, gas),
        reactions = (rx,),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
    dz = VBO_L_PMMA / n
    n_kw = round(Int, VBO_L_KW / dz)
    z_if = n_kw * dz
    bc_set = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(; flux = s * 50e3, absorptivity = 0.96) +
                  ConvectiveBC(h = s * 30.0, T_inf = VBO_T_INF),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
    )
    te = t_end / s
    return PyrolysisProblem(;
        mesh = create_uniform_mesh((n + n_kw) * dz, n + n_kw),
        material,
        bc_set,
        T_initial = z -> VBO_T0,
        ξ_initial = [
            z -> (z > z_if ? VBO_RHO_PMMA : 0.0),
            z -> (z > z_if ? 0.0 : VBO_RHO_KW),
            z -> 0.0,
        ],
        tspan = (0.0, te),
        dt_initial = 1e-3 / s,
        output_times = collect(0.0:(0.5/s):te),
    )
end

# V-BO-01's deck with its incident flux taken to `f` times its value by a C¹
# ramp of `ramp` seconds from `t_step`.
stepped_problem(n, f; t_step = 150.0, ramp = 5.0, t_end = 800.0) = _vbo_problem(
    n;
    flux = t -> 50e3 * (1 + (f - 1) * _smoothstep((t - t_step) / ramp)),
    t_end,
)

# A 10 mm permeable wood column, Darcy–Fick transport and in-depth
# absorption, `n` cells.
function long_layer_problem(n; t_end = 300.0)
    wood = SolidComponent(
        :wood,
        ρ = 500.0,
        ρ_intrinsic = 1500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
        α = 5.0,
        κ = 1e-14,
        λ = 1e-5,
    )
    gas = GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    material = Material(
        name = :long_layer,
        components = (wood, gas),
        reactions = (Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = 300e3),),
    )
    thermal = (
        top = RadiativeFluxBC(flux = 50e3) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    )
    mass = (
        top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
        bottom = ImpermeableBC(),
    )
    return PyrolysisProblem(;
        mesh = create_uniform_mesh(0.01, n),
        material,
        bc_set = BoundaryConditionSet(; thermal, mass),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = round(Int, t_end) + 1)),
    )
end

# ----------------------------------------------------------------------------
# The configurations
# ----------------------------------------------------------------------------

# The QoIs of a variant of V-BO-01's deck: those that need no reference.
# `surface = false` leaves out the surface balances' QoI, which counts every
# thermal face whose balance did not converge: an adiabatic face reports none.
variant_qois(id; l_kw = VBO_L_KW, l_layer = VBO_L_PMMA, surface = true) =
    sol -> begin
        m = _vbo_measures(sol; l_kw, id)
        [
            _vbo_common_qois(
                m;
                id,
                mass_loss_exact = _vbo_mass_loss_exact(l_layer, l_kw),
                surface_warnings = surface,
            )...,
            _vbo_invocation_qoi(sol, id),
        ]
    end

case_config(id; cells, layer = VBO_L_PMMA, thresholds = ()) = (;
    name = id,
    build = LADDER_CASES[id].build,
    cells,
    ci = LADDER_CASES[id].ci.n_cells,
    layer,
    keywords = LADDER_CASES[id].solve_kwargs,
    solve_fn = LADDER_CASES[id].solve_fn,
    qois = LADDER_CASES[id].qois,
    thresholds,
    enthalpy = true,
    rescaled = nothing,
)

variant_config(name, build; cells = (30, 60, 120), rescaled = nothing) = (;
    name,
    build,
    cells,
    ci = 60,
    layer = VBO_L_PMMA,
    keywords = VBO_SOLVE_KWARGS,
    solve_fn = _vbo_solve_fn(name),
    qois = variant_qois(name; surface = rescaled === nothing),
    thresholds = (),
    enthalpy = false,
    rescaled,
)

const CONFIGS = [
    case_config("V-BO-01"; cells = (30, 60, 120), thresholds = (3e-3, 3e-2)),
    case_config("V-BO-06"; cells = (30, 60, 120), thresholds = (3e-3, 3e-2)),
    case_config("V-BO-07"; cells = (30, 60, 120)),
    case_config(
        "V-BO-08";
        cells = (10, 20, 40),
        layer = VBO08_L_PMMA,
        thresholds = (3e-3, 3e-2),
    ),
    variant_config("25 kW/m²", n -> _vbo_problem(n; flux = 25e3, t_end = 1500.0)),
    variant_config("100 kW/m²", n -> _vbo_problem(n; flux = 100e3)),
    variant_config("A/2", n -> _vbo_problem(n; A = 1.35e11 / 2)),
    variant_config("2A", n -> _vbo_problem(n; A = 2 * 1.35e11)),
    variant_config("restricted, s = 1", n -> restricted_problem(n, 1.0); rescaled = 1.0),
    variant_config("restricted, s = 2", n -> restricted_problem(n, 2.0); rescaled = 2.0),
    variant_config("restricted, s = 1/2", n -> restricted_problem(n, 0.5); rescaled = 0.5),
    variant_config("flux × 0.5 at 150 s", n -> stepped_problem(n, 0.5)),
    variant_config("flux × 1.5 at 150 s", n -> stepped_problem(n, 1.5)),
    (;
        name = "long layer",
        build = long_layer_problem,
        cells = (25, 100),
        ci = 100,
        layer = 0.01,
        keywords = (
            use_ale = true,
            handle_depletion = true,
            min_thickness = 1e-6,
            radiation_model = BEER_LAMBERT,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        solve_fn = _vbo_solve_fn("long layer"),
        qois = nothing,
        thresholds = (),
        enthalpy = false,
        rescaled = nothing,
    ),
]
const CONFIG = Dict(c.name => c for c in CONFIGS)

# The equidistribution time a rule gives the configuration `c`.
rule_time(rule, c) =
    rule === :fixed ? Cl.EQUIDISTRIBUTION_TIME :
    rule === :A ? Cl.EQUIDISTRIBUTION_TIME * c.layer / H_REF :
    rule === :B ? Cl.ActivationTime() : Cl.ActivationTime(λ = 0.5)

# ----------------------------------------------------------------------------
# The first activation
# ----------------------------------------------------------------------------

# The first accepted step on which a cell of a moving layer is below `r_hi` of
# its initial thickness, the crossing interpolated linearly between that
# step's start and end.
mutable struct Activation
    t_a::Float64
    last_t::Float64
    last_r::Float64
end
Activation() = Activation(NaN, NaN, NaN)

function activation_recorder(cp, act::Activation)
    layout = cp.def.layout
    corr = cp.def.closure.correction
    h0 = cp.def.setup.reference.thickness
    cells = [i for (lo, hi) in corr.layers for i = lo:(hi-1)]
    ratio(u) = minimum(
        (
            (
                u[Pyrolysis.States.position_index(layout, i + 1)] -
                u[Pyrolysis.States.position_index(layout, i)]
            ) / h0[i] for i in cells
        );
        init = Inf,
    )
    return DiscreteCallback(
        (u, t, integ) -> isnan(act.t_a),
        integ -> begin
            r = ratio(integ.u)
            if r < corr.r_hi
                act.t_a =
                    isnan(act.last_t) ? integ.t :
                    act.last_t +
                    (act.last_r - corr.r_hi) / (act.last_r - r) * (integ.t - act.last_t)
            end
            act.last_t, act.last_r = integ.t, r
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
end

# ----------------------------------------------------------------------------
# The runs
# ----------------------------------------------------------------------------

function ladder_run(c, rule, n; ε_m = 1e-2, reltol = 1e-6, form = :temperature)
    problem = Base.invokelatest(c.build, n)
    τ = rule_time(rule, c)
    series = Series()
    act = Activation()
    bound = Ref{Any}(nothing)
    instrument = cp -> begin
        bound[] = cp
        CallbackSet(recorder(cp, series), activation_recorder(cp, act))
    end
    started = time()
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(
            c.solve_fn,
            problem;
            c.keywords...,
            energy_state = form,
            threshold = ε_m,
            reltol,
            correction = Cl.EquidistributionCorrection(; τ),
            instrument,
        )
    end
    wall = time() - started
    corr = bound[].def.closure.correction
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    qois = c.qois === nothing ? [] : Base.invokelatest(c.qois, sol)
    counted = filter(q -> !startswith(q.name, "one invocation"), qois)
    failing = [q.name for q in counted if !(q.error < q.tol)]
    t_inv = isempty(inv) ? NaN : first(inv).t
    τ_layer = isempty(corr.time) ? NaN : first(corr.time)
    life = t_inv - act.t_a
    last_event = isempty(sol.events) ? Inf : maximum(e.t for e in sol.events)
    run = (;
        config = c.name,
        rule,
        n,
        ε_m,
        reltol,
        form,
        τ = τ_layer,
        set_at = isempty(corr.set_at) ? NaN : first(corr.set_at),
        t_a = act.t_a,
        life,
        sync = τ_layer / π^2 / life,
        retcode = sol.retcode,
        outcome = sol.outcome,
        t_end = sol.t[end],
        invocations = inv,
        item1 = length(inv) == 1 && first(inv).largest <= 5ε_m,
        qois_total = length(counted),
        qois_failing = failing,
        qoi_pass = isempty(failing),
        steps = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        jacobians = sol.statistics.jacobians,
        wall,
        series = ordered(series),
        last_event,
        t = copy(sol.t),
        T = copy(sol.T),
        z = copy(sol.z),
        events = [e.t for e in sol.events],
    )
    pk, t_pk = isempty(run.series.t) ? (NaN, NaN) : main_peak(run)
    return (; run..., peak = 1000pk, t_peak = t_pk)
end

# The runs of the ladder, in order.
function ladder_keys()
    keys = Any[]
    for c in CONFIGS, rule in RULES
        for n in c.cells
            push!(keys, (c.name, rule, n, 1e-2, 1e-6, :temperature))
        end
        for ε_m in c.thresholds
            push!(keys, (c.name, rule, c.ci, ε_m, 1e-6, :temperature))
        end
        for n in (c.rescaled == 1.0 ? c.cells : (c.ci,))
            push!(keys, (c.name, rule, n, 1e-2, 1e-7, :temperature))
        end
        c.enthalpy && push!(keys, (c.name, rule, c.ci, 1e-2, 1e-6, :enthalpy))
    end
    return keys
end

# ----------------------------------------------------------------------------
# The criteria
# ----------------------------------------------------------------------------

# Whether every successive value falls at least `FALL`× (a value at or below
# zero counts as fallen; one above zero after one at or below it does not).
falls(x) = all(x[k+1] <= 0 || (x[k] > 0 && x[k] / x[k+1] >= FALL) for k = 1:(length(x)-1))

# W9's hard criterion on a family: its runs at the configuration's
# resolutions, coarse to fine.
function refinement(runs)
    finest = runs[end]
    t_pk = finest.t_peak
    t_end = finest.t_end
    f = [Piecewise(r.series) for r in runs]
    exc = [excess(runs[k].series, f[k+1], t_pk) for k = 1:(length(runs)-1)]
    succ = [l1(f[k], f[k+1], t_pk, t_end) for k = 1:(length(f)-1)]
    peaks = [r.peak for r in runs]
    steps = abs.(diff(peaks))
    a = falls(exc)
    b = falls(succ)
    peak = all(steps[k+1] <= steps[k] for k = 1:(length(steps)-1))
    qois = all(r -> r.qoi_pass, runs)
    return (; exc, succ, peaks, a, b, peak, qois, hard = a && b && peak && qois)
end

# The restricted deck at `s` against s = 1, in dimensionless time.
function pair_differences(a, b, s)
    length(a.t) == length(b.t) || return (T = Inf, z = Inf, events = Inf)
    nan_same = isnan.(a.T) == isnan.(b.T) && length(a.events) == length(b.events)
    nan_same || return (T = Inf, z = Inf, events = Inf)
    return (
        T = maximum(abs, filter(isfinite, a.T .- b.T); init = 0.0),
        z = maximum(abs, filter(isfinite, a.z .- b.z); init = 0.0),
        events = maximum(abs, a.events .- b.events .* s; init = 0.0),
    )
end

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

function ladder_main()
    load_start = Sys.loadavg()
    started = time()
    runs = Dict{Any,Any}()
    for key in ladder_keys()
        name, rule, n, ε_m, reltol, form = key
        r = ladder_run(CONFIG[name], rule, n; ε_m, reltol, form)
        runs[key] = r
        println(
            "$(name) $(rule) n=$(n) ε_m=$(ε_m) reltol=$(reltol) $(form): τ=$(_g(r.τ)) ",
            "t_a=$(_g(r.t_a)) $(r.retcode)/$(r.outcome) invocations=$(length(r.invocations)) ",
            "steps=$(r.steps) qois=$(r.qois_total - length(r.qois_failing))/$(r.qois_total) ",
            "wall=$(_g(r.wall)) s",
        )
        flush(stdout)
    end
    families = Any[]
    for c in CONFIGS, rule in RULES
        length(c.cells) >= 3 || continue
        fam = [runs[(c.name, rule, n, 1e-2, 1e-6, :temperature)] for n in c.cells]
        push!(families, (; config = c.name, rule, refinement(fam)...))
    end
    pairs = Any[]
    base = CONFIG["restricted, s = 1"]
    for rule in RULES, n in base.cells
        ref = runs[(base.name, rule, n, 1e-2, 1e-6, :temperature)]
        self =
            pair_differences(ref, runs[(base.name, rule, n, 1e-2, 1e-7, :temperature)], 1.0)
        for s in (2.0, 0.5)
            name = s == 2 ? "restricted, s = 2" : "restricted, s = 1/2"
            d = pair_differences(ref, runs[(name, rule, n, 1e-2, 1e-6, :temperature)], s)
            pass =
                d.T <= PAIR_BOUND * self.T &&
                d.z <= PAIR_BOUND * self.z &&
                d.events <= PAIR_BOUND * self.events
            push!(pairs, (; rule, n, s, d, self, pass))
        end
    end
    load_end = Sys.loadavg()
    write_ladder(runs, families, pairs, load_start, load_end, time() - started)
    println(
        "wrote $(OUT)/tau_ladder.{csv,md}; load $(load_start) -> $(load_end); ",
        "$(round(Int, time() - started)) s",
    )
    return (; runs, families, pairs)
end

# ----------------------------------------------------------------------------
# The tables
# ----------------------------------------------------------------------------

_inv(r) =
    isempty(r.invocations) ? "0" :
    join(("$(v.cells) ($(_g(100v.largest)) %)" for v in r.invocations), "; ")
_rule(rule) = rule === :fixed ? "fixed" : rule === :B05 ? "B(0.5)" : string(rule)

# §7 item 11's reading of one rule.
function acceptance(rule, runs, families, pairs)
    item1_rows = [
        k for k in keys(runs) if
        k[2] === :fixed && k[5] == 1e-6 && runs[k].item1 && k[1] in keys(LADDER_CASES)
    ]
    held = all(runs[(k[1], rule, k[3:end]...)].item1 for k in item1_rows)
    thin = runs[("V-BO-08", rule, 10, 3e-2, 1e-6, :temperature)].item1
    rescaled = all(p -> p.pass, filter(p -> p.rule === rule, pairs))
    burnout = filter(f -> f.rule === rule, families)
    stepped = filter(f -> startswith(f.config, "flux ×"), burnout)
    long = runs[("long layer", rule, 100, 1e-2, 1e-6, :temperature)]
    return (;
        item1 = held,
        item1_rows = length(item1_rows),
        thin,
        rescaled,
        hard = all(f -> f.hard, burnout),
        hard_failing = [f.config for f in burnout if !f.hard],
        stepped = all(f -> f.hard, stepped),
        long_steps = long.steps,
        long_invocations = length(long.invocations),
        long_outcome = long.outcome,
    )
end

function write_ladder(runs, families, pairs, load_start, load_end, elapsed)
    mkpath(OUT)
    ordered_keys = ladder_keys()
    open(joinpath(OUT, "tau_ladder.csv"), "w") do io
        println(
            io,
            "config,rule,cells,eps_m,reltol,form,tau,set_at,t_a,life,sync,retcode,outcome,t_end,",
            "invocations,cells_merged,largest_fraction,item1,qois_passing,qois_total,peak,t_peak,",
            "steps,rejected,jacobians,wall",
        )
        for k in ordered_keys
            r = runs[k]
            println(
                io,
                join(
                    (
                        "\"$(r.config)\"",
                        r.rule,
                        r.n,
                        r.ε_m,
                        r.reltol,
                        r.form,
                        r.τ,
                        r.set_at,
                        r.t_a,
                        r.life,
                        r.sync,
                        r.retcode,
                        r.outcome,
                        r.t_end,
                        length(r.invocations),
                        sum((v.cells for v in r.invocations); init = 0),
                        isempty(r.invocations) ? NaN :
                        maximum(v.largest for v in r.invocations),
                        r.item1,
                        r.qois_total - length(r.qois_failing),
                        r.qois_total,
                        r.peak,
                        r.t_peak,
                        r.steps,
                        r.rejected,
                        r.jacobians,
                        r.wall,
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(OUT, "tau_ladder.md"), "w") do io
        println(io, "# The layer time's ladder (W3, §7 item 11)\n")
        println(
            io,
            "Produced by `test/newcore/gates/tau_ladder.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
            "$(Threads.nthreads()) thread; load average $(round.(load_start; digits = 1)) at the start ",
            "and $(round.(load_end; digits = 1)) at the end; $(round(Int, elapsed)) s. The finite-volume ",
            "core with the equidistribution correction (band 0.1, activation 0.3/0.1, monitor without its ",
            "mass term) under three rules for the layer's time: fixed, $(Cl.EQUIDISTRIBUTION_TIME) s; ",
            "A, that time times the layer's initial thickness over $(1000H_REF) mm; B, ",
            "`ActivationTime()` (λ = $(Cl.ACTIVATION_TIME_FACTOR)). The script's header defines every ",
            "configuration, run and measure.\n",
        )
        println(io, "## §7 item 11\n")
        println(
            io,
            "Per rule: Phase 2 item 1 on every case row the fixed time passes (one invocation, the ",
            "largest merged fraction at most 5 ε_m); V-BO-08 at ε_m = 3e-2; the rescaled pairs within ",
            "$(PAIR_BOUND)× the solve's own error; W9's hard criterion on every family, the flux-step ",
            "families among them; the long layer's steps at 100 cells (the band alone took 470 in ",
            "Phase 2's review).\n",
        )
        println(
            io,
            "| rule | item 1 where it held | V-BO-08 at 3e-2 | rescaled pairs | hard criterion (failing families) | flux steps | long layer, 100 cells: steps, invocations, outcome |",
        )
        println(io, "|---|---|---|---|---|---|---|")
        for rule in RULES
            a = acceptance(rule, runs, families, pairs)
            println(
                io,
                "| $(_rule(rule)) | $(_yes(a.item1)) ($(a.item1_rows) rows) | $(_yes(a.thin)) | ",
                "$(_yes(a.rescaled)) | $(_yes(a.hard))$(isempty(a.hard_failing) ? "" : " ($(join(a.hard_failing, "; ")))") | ",
                "$(_yes(a.stepped)) | $(a.long_steps), $(a.long_invocations), $(a.long_outcome) |",
            )
        end
        println(io, "\n## The families: W9's hard criterion\n")
        println(
            io,
            "(a) the largest per-step excess of N over 2N after the finest run's main peak ",
            "[g m⁻² s⁻¹], its ratio; (b) the successive L1 [g/m²], its ratio; both falling at least ",
            "$(FALL)× per doubling. The main peak [g m⁻² s⁻¹] at each resolution. QoIs: every run's, ",
            "the invocation QoI aside.\n",
        )
        println(
            io,
            "| configuration | rule | (a); ratio | (a) | (b); ratio | (b) | peaks | peak | QoIs | hard |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for f in families
            println(
                io,
                "| $(f.config) | $(_rule(f.rule)) | $(_list(f.exc)); $(_list(ratios(f.exc); scale = 1)) | ",
                "$(_yes(f.a)) | $(_list(f.succ)); $(_list(ratios(f.succ); scale = 1)) | $(_yes(f.b)) | ",
                "$(join(_g.(f.peaks; d = 5), ", ")) | $(_yes(f.peak)) | $(_yes(f.qois)) | $(_yes(f.hard)) |",
            )
        end
        println(io, "\n## The cost at the CI resolution\n")
        println(
            io,
            "Accepted steps (rejected) and Jacobians at each configuration's CI resolution, ε_m = 1e-2, ",
            "reltol 1e-6, the temperature form; the ratio of accepted steps to the fixed time's.\n",
        )
        println(
            io,
            "| configuration | cells | ",
            join(("$(_rule(r)): steps (rejected), Jacobians, ratio" for r in RULES), " | "),
            " |",
        )
        println(io, "|---|---|", repeat("---|", length(RULES)))
        for c in CONFIGS
            base = runs[(c.name, :fixed, c.ci, 1e-2, 1e-6, :temperature)]
            cells = join(
                (
                    begin
                        r = runs[(c.name, rule, c.ci, 1e-2, 1e-6, :temperature)]
                        "$(r.steps) ($(r.rejected)), $(r.jacobians), $(_g(r.steps / base.steps))"
                    end for rule in RULES
                ),
                " | ",
            )
            println(io, "| $(c.name) | $(c.ci) | $(cells) |")
        end
        println(io, "\n## The rescaled pairs\n")
        println(
            io,
            "The restricted deck at s against s = 1 in dimensionless time: the largest saved ",
            "temperature [K] and node-position [m] differences and the event-time shift [s], against ",
            "$(PAIR_BOUND)× the s = 1 run's own differences between reltol 1e-6 and 1e-7.\n",
        )
        println(
            io,
            "| rule | cells | s | ΔT | Δz | Δ events | own ΔT, Δz, Δ events | within |",
        )
        println(io, "|---|---|---|---|---|---|---|---|")
        for p in pairs
            println(
                io,
                "| $(_rule(p.rule)) | $(p.n) | $(p.s) | $(_g(p.d.T)) | $(_g(p.d.z)) | $(_g(p.d.events)) | ",
                "$(_g(p.self.T)), $(_g(p.self.z)), $(_g(p.self.events)) | $(_yes(p.pass)) |",
            )
        end
        println(io, "\n## The runs\n")
        println(
            io,
            "τ: the layer's time [s] (B: set at t [s]); t_a: the first activation [s]; life: from t_a ",
            "to the layer's first invocation [s]; T_sync/life: τ/π² over that life. Invocations: each ",
            "one's cells (largest merged fraction). QoIs: passing of counted, the failing named. Wall: ",
            "the instrumented solve's, on a loaded machine.\n",
        )
        println(
            io,
            "| configuration | rule | cells | ε_m | reltol | form | τ (set at) | t_a | life | T_sync/life | retcode, outcome, end | invocations | item 1 | QoIs | failing | peak (t) | steps (rejected) | Jacobians | wall [s] |",
        )
        println(
            io,
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        )
        for k in ordered_keys
            r = runs[k]
            set = isnan(r.set_at) || r.rule !== :B ? "" : " ($(_g(r.set_at; d = 5)))"
            println(
                io,
                "| $(r.config) | $(_rule(r.rule)) | $(r.n) | $(r.ε_m) | $(r.reltol) | $(r.form) | ",
                "$(_g(r.τ; d = 4))$(set) | $(_g(r.t_a; d = 5)) | $(_g(r.life; d = 4)) | $(_g(r.sync)) | ",
                "$(r.retcode), $(r.outcome), $(_g(r.t_end; d = 5)) | $(_inv(r)) | $(_yes(r.item1)) | ",
                "$(r.qois_total - length(r.qois_failing))/$(r.qois_total) | $(join(r.qois_failing, "; ")) | ",
                "$(_g(r.peak; d = 5)) ($(_g(r.t_peak; d = 5))) | $(r.steps) ($(r.rejected)) | ",
                "$(r.jacobians) | $(_g(r.wall)) |",
            )
        end
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && ladder_main()
