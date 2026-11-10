# ============================================================================
# Solve
# ============================================================================

"""
	solve(problem::PyrolysisProblem; energy_state = :temperature, abstol = 1e-8,
	      reltol = 1e-6, threshold = nothing, pressure_tolerance = 1.0,
	      maxiters = 1_000_000, dt = problem.dt_initial, adaptive = true,
	      dtmin = 1e-12, dtmax = Inf, saveat = problem.output_times,
	      callback = nothing, progress = false, use_ale = false,
	      min_thickness = 1e-6, radiation_model = SURFACE_ABSORPTION,
	      integrator_policy = default_policy(), sources = NoSources(),
	      deactivation = false, band = nothing, correction = nothing,
	      lagrangian = nothing) -> PyrolysisSolution

Solve `problem`: [`validate_problem`](@ref), [`setup`](@ref), [`integrate`](@ref).

The integrator is handed the mass-matrix system with the colored Jacobian, the
per-row absolute tolerances (`abstol` is the temperature tolerance [K]; the mass
and enthalpy rows take theirs from the cell's initial condensed mass and
`threshold`, whose default [`setup`](@ref) states), no initialization of its own,
and [`DomainCheck`](@ref) as its
domain test (admissibility, porosity and the constraint residuals of the
algebraic rows, on the proposals its error test accepts). `dt` is the first
step, and with `adaptive = false` the fixed step. `integrator_policy` holds the
method and its linear solver ([`rodas5p_policy`](@ref) and its siblings). A per-step
observer records the largest constraint residual of the accepted states, ends
the solve with `ReturnCode.Failure` when the domain test gave up, caps the step
after a nonlinear failure of a Newton-based policy
([`FAILURE_MEMORY_GROWTH`](@ref); a Rosenbrock method has none) and applies the
policy's flag rule. A discrete `callback`, or each discrete callback of a
`CallbackSet`, is a read-only observer: the policy's flag rule follows its
`affect!`, so it leaves the trajectory unchanged under every policy; a
continuous callback's event raises the flag as any event does. `callback` may
also be a [`ConservationRecorder`](@ref), or a tuple of recorders and
callbacks. `progress` draws a progress bar, a read-only observer as well.

`use_ale`, `correction`, `lagrangian`, `deactivation` and
`band` are those of [`setup`](@ref). On a moving mesh the solve ends with return code
`:Terminated` and outcome `:min_thickness` when the thickness falls below
`min_thickness` [m] (default 1e-6). With `deactivation`, cells are deactivated in events
([`deactivate!`](@ref)): an event is located inside the
step at the first time a live cell reaches its threshold (the first root of
[`threshold_margin`](@ref)) and applied at the interpolated state
there, the threshold and the band widened by `reltol` (`deactivate!`'s
`slack`): the cells that go are the runs the widened tests take, which hold the
crossing cell; at the end of every step the cells below the
threshold go in one event, the safety net for a crossing the location does not
see. The solve ends with `:Terminated` when the last cell goes or when a cell below the
threshold lies beneath live cells the pass cannot take it through, a cavity
(outcome `:interior`, with a warning naming the cell), and with `:Failure` when an event leaves
the state inadmissible (a merge's receiver held to its own mass tolerances,
as the domain test holds every later proposal) or no consistent node
velocities. The state a solve
ends on is always saved, whether or not `saveat` names its time.
"""
function solve(
    problem::PyrolysisProblem;
    energy_state::Symbol = :temperature,
    abstol::Real = 1e-8,
    reltol::Real = 1e-6,
    threshold::Union{Nothing,Real} = nothing,
    pressure_tolerance::Real = 1.0,
    maxiters::Int = 1_000_000,
    dt::Real = problem.dt_initial,
    adaptive::Bool = true,
    dtmin::Real = 1e-12,
    dtmax::Real = Inf,
    saveat = problem.output_times,
    callback = nothing,
    progress::Bool = false,
    use_ale::Bool = false,
    min_thickness::Real = 1e-6,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    integrator_policy::IntegratorPolicy = default_policy(),
    sources = NoSources(),
    deactivation::Bool = false,
    band::Union{Nothing,Real} = nothing,
    correction::Union{Nothing,Correction} = nothing,
    lagrangian::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    validate_problem(problem)
    started = time()
    cp = setup(
        problem;
        energy_state,
        reltol,
        abstol,
        threshold,
        pressure_tolerance,
        use_ale,
        radiation_model,
        time_scale = dt,
        sources,
        deactivation,
        band,
        correction,
        lagrangian,
    )
    return integrate(
        cp;
        maxiters,
        dt,
        adaptive,
        dtmin,
        dtmax,
        saveat,
        callback,
        min_thickness,
        integrator_policy,
        progress,
        started,
    )
end

"""
	integrate(cp::PreparedProblem; u0 = cp.u0, maxiters = 1_000_000,
	          dt = cp.problem.dt_initial, adaptive = true, dtmin = 1e-12,
	          dtmax = Inf, saveat = cp.problem.output_times, callback = nothing,
	          min_thickness = 1e-6, integrator_policy = default_policy(),
	          progress = false, started = time()) -> PyrolysisSolution

Integrate a set-up problem from the state `u0`, whose algebraic rows must be
consistent ([`consistent_state!`](@ref)), with the tolerances it was set up
with. The keywords are those of [`solve`](@ref); `started` is the time the
reported solve time counts from. An [`inadmissible`](@ref) `u0` is refused; a fixed-step solve,
which the domain test never sees, ends with `ReturnCode.Failure` at the
first accepted state that is inadmissible. A solve that fails ends on the
last state it accepted. Every call starts with every cell live and, in the inventory
form, the caloric inversion from the set-up temperatures, so repeated
integrations of one problem are bit-identical.
"""
function integrate(
    cp::PreparedProblem;
    u0::AbstractVector = cp.u0,
    maxiters::Int = 1_000_000,
    dt::Real = cp.problem.dt_initial,
    adaptive::Bool = true,
    dtmin::Real = 1e-12,
    dtmax::Real = Inf,
    saveat = cp.problem.output_times,
    callback = nothing,
    min_thickness::Real = 1e-6,
    integrator_policy::IntegratorPolicy = default_policy(),
    progress::Bool = false,
    started::Float64 = time(),
)
    min_thickness > 0 ||
        throw(ArgumentError("min_thickness must be positive, got $min_thickness"))
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    cells = cp.deactivation
    cells === nothing || reset_cells!(cells, layout)
    inadmissible(u0, layout, def.setup, cp.abstol) && throw(
        ArgumentError(
            "the initial state is inadmissible: a live component mass below its " *
            "tolerance, a live cell at or below absolute zero, or crossed nodes",
        ),
    )
    corr = _layer_correction(def)
    corr === nothing || reset_layer_times!(corr)
    copyto!(ws.T_warm, cp.T0)
    policy = integrator_policy

    check = DomainCheck(cp)
    rhs! = let check = check
        (du, u, p, t) -> begin
            right_hand_side!(du, check, u, t)
            return nothing
        end
    end
    evaluated! = let check = check
        (du, u, p, t) -> begin
            evaluate_right_hand_side!(du, check, u, t)
            return nothing
        end
    end
    J0 = jac_prototype(cj)
    jac! = let cj = cj
        (J, u, p, t) -> begin
            jacobian!(J, cj, u, t)
            return nothing
        end
    end
    M = all(is_differential(layout)) ? I : mass_matrix(layout)
    enthalpy = energy_state(layout) isa EnthalpyState
    f = ODEFunction{true,FullSpecialize}(
        ErasedRHS(rhs!);
        mass_matrix = M,
        jac = ErasedJacobian{typeof(J0)}(jac!),
        tgrad = ErasedRHS(_time_derivative(enthalpy ? evaluated! : rhs!, cp)),
        jac_prototype = J0,
    )
    ode = ODEProblem(f, Vector{Float64}(u0), cp.problem.tspan)

    # The domain test and every callback reach the problem through
    # fixed-signature wrappers, so the integrator's type does not depend on it.
    out_of_domain = _DomainWrapper(check)
    monitor = check.monitor
    cap = NewtonCap()
    violation = _ViolationWrapper(
        let layout = layout, setup = def.setup, abstol = cp.abstol
            (u, t) -> begin
                v=_violation(u, layout, setup, abstol)
                v==0 || return v
                ds=_derived(check, u, t)
                return _porosity_out(check, ds, setup.gas) ||
                       any(x->x<0, check.ws.contributions.surface_pressure_raw) ?
                       -1 : 0
            end
        end,
    )
    observer = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> _observe!(integrator, monitor, cap, policy, violation);
        save_positions = (false, false),
    )
    outcome = Ref(:none)
    thin = _thickness_termination(layout, Float64(min_thickness), outcome)
    events = _deactivation_callback(cp, check, policy)
    layer_times = _layer_time_callback(cp, check, policy)
    recs, callback = _split_recorders(callback)
    bound = [_BoundRecorder(r, cp, _crossings(cp, min_thickness)) for r in recs]
    # The last accepted state, after every callback of its step; then the
    # conservation recorders' commits and the step's saved states read again
    # (`_resave_step!`).
    accepted = AcceptedState(copy(ode.u0), Float64(cp.problem.tspan[1]))
    start = MethodDerivative(ode.u0, _algebraic_rows(def, cp.T0))
    deactivations = cells === nothing ? DeactivationEvent[] : cells.events
    recorder = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator ->
            _record!(integrator, accepted, policy, bound, start, deactivations);
        save_positions = (false, false),
    )
    observers = (
        map(b -> _reader(b, start), bound)...,
        (progress ? (_progress_observer(policy, cp.problem.tspan),) : ())...,
        (callback === nothing ? () : (_observing(callback, policy),))...,
    )
    callbacks =
        CallbackSet(observer, events..., layer_times..., thin..., observers..., recorder)

    pressure_norm=PressureErrorNorm(cp)
    integrator = SciMLBase.init(
        ode,
        policy.alg;
        abstol = cp.abstol,
        reltol = cp.reltol,
        maxiters,
        dtmin = Float64(dtmin),
        dtmax = Float64(dtmax),
        dt = Float64(dt),
        adaptive,
        saveat = saveat isa Number ? Float64(saveat) : collect(Float64, saveat),
        callback = callbacks,
        isoutofdomain = out_of_domain,
        internalnorm = ErasedPressureNorm(pressure_norm.callback),
        initializealg = NoInit(),
        save_end = true,
    )
    # The integrator takes a copy; an event's permutation reaches it through
    # this one array.
    integrator.opts.abstol = cp.abstol
    bind_pressure_norm!(pressure_norm, integrator)
    bind_error_test!(check, integrator)
    bind_tolerances!(integrator, cp.abstol)
    raw = SciMLBase.solve!(integrator)
    monitor.cold == 0 ||
        @warn "solve(): stopped at t = $(raw.t[end]): the last proposed state, at t = $(monitor.cold_t), takes cell $(monitor.cold) to or below absolute zero"
    SciMLBase.successful_retcode(raw.retcode) ||
        !(monitor.transport_t >= accepted.t) ||
        @warn "solve(): stopped at t = $(raw.t[end]): a proposal after the last accepted state, at t = $(monitor.transport_t), raised $(monitor.transport_message)"
    # A domain outcome the integrator ran into without the test giving up (the
    # step reached `dtmin` first) is the reason the solve stopped: name it.
    SciMLBase.successful_retcode(raw.retcode) ||
        monitor.failed ||
        monitor.outcome === :none ||
        @warn "solve(): stopped at t = $(raw.t[end]): the domain test rejected the proposals after the last accepted state, $(monitor.outcome): $(monitor.message)"
    # A solve its controller stopped saves its last proposal, which it
    # rejected, at the last accepted time: report the accepted state there.
    SciMLBase.successful_retcode(raw.retcode) ||
        raw.t[end] != accepted.t ||
        copyto!(raw.u[end], accepted.u)
    elapsed = time() - started
    s = raw.stats
    stats = SolverStatistics(
        accepted = s.naccept,
        rejected = s.nreject,
        rhs_evaluations = check.evaluations[],
        jacobians = s.njacs,
        linear_solves = s.nsolve,
        nonlinear_failures = s.nnonlinconvfail,
        domain_rejections = check.monitor.rejections,
        transport_rejections = check.monitor.transport,
        pressure_error_evaluations = pressure_norm.evaluations[],
        largest_constraint_residual = check.monitor.largest,
        invocations = cells === nothing || isempty(cells.events) ? 0 :
                      last(cells.events).invocation,
    )
    return PyrolysisSolution(
        cp.problem,
        raw.t,
        raw.u,
        def,
        ws,
        Symbol(raw.retcode),
        elapsed,
        stats,
        raw;
        deactivation = cells,
        outcome = monitor.outcome === :none ? outcome[] : monitor.outcome,
    )
end

# Whether the residual of `cp` can read the time: through its boundary
# conditions, its external sources or a prescribed mesh velocity.
function _reads_time(cp::PreparedProblem)
    closure = cp.def.closure
    return !(cp.def.sources isa NoSources) ||
           (closure isa BalanceRowClosure && closure.correction isa PrescribedCorrection) ||
           reads_time(cp.problem.bc_set)
end

# The step of the difference in time, relative to the time and at least itself.
const TIME_STEP = sqrt(eps(Float64))

# The signed step of the difference at `t`, inside the time span: forward, backward
# where the forward point lies beyond the span's end, and to the farther end of a
# span too short for either.
function _time_step(t, t_start, t_end)
    step = max(TIME_STEP * abs(t), TIME_STEP)
    t + step <= t_end && return step
    t - step >= t_start && return -step
    return t_end - t >= t - t_start ? t_end - t : t_start - t
end

# The derivative of `rhs!` in time, which a Rosenbrock method reads on every
# step. A residual that does not read the time has exactly zero; any other takes
# the difference of [`_time_step`](@ref). The caloric inversion's result depends by
# a few ulps on the temperatures it starts from, and that round-off over the step
# would read as a time derivative, so both residuals start from the same
# temperatures and leave them as they were. In the enthalpy form `rhs!` must
# evaluate: a residual the domain test remembers started from others.
function _time_derivative(rhs!, cp::PreparedProblem)
    _reads_time(cp) || return (dT, u, p, t) -> begin
        fill!(dT, 0.0)
        return nothing
    end
    ahead, here = zeros(length(cp.u0)), zeros(length(cp.u0))
    t_start, t_end = cp.problem.tspan
    warm_start, saved = cp.ws.T_warm, copy(cp.ws.T_warm)
    return (dT, u, p, t) -> begin
        h = _time_step(t, t_start, t_end)
        copyto!(saved, warm_start)
        rhs!(ahead, u, p, t + h)
        copyto!(warm_start, saved)
        rhs!(here, u, p, t)
        copyto!(warm_start, saved)
        @. dT = (ahead - here) / h
        return nothing
    end
end

"""
	FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS

After an accepted step whose nonlinear solves needed retries, the proposed
step stays at or below `FAILURE_MEMORY_GROWTH` times that step for the next
`FAILURE_MEMORY_STEPS` accepted steps.
"""
const FAILURE_MEMORY_GROWTH = 1.25
const FAILURE_MEMORY_STEPS = 20

# The Newton-failure step cap: failures seen so far, accepted steps the cap
# still holds, and the cap.
mutable struct NewtonCap
    failures::Int
    hold::Int
    cap::Float64
end

NewtonCap() = NewtonCap(0, 0, Inf)

const _DomainWrapper = FunctionWrapper{Bool,Tuple{_StateVector,NullParameters,Float64}}
const _ViolationWrapper = FunctionWrapper{Int,Tuple{_StateVector,Float64}}

# The last state the integrator accepted and its time, recorded by a
# read-only observer after every other callback of the step.
mutable struct AcceptedState
    u::Vector{Float64}
    t::Float64
end

function _record!(
    integrator,
    accepted::AcceptedState,
    policy,
    bound,
    start::MethodDerivative,
    deactivations::Vector{DeactivationEvent},
)
    copyto!(accepted.u, integrator.u)
    accepted.t = integrator.t
    for b in bound
        _commit_step!(b, integrator)
    end
    _resave_step!(integrator, start, length(deactivations))
    apply_observer_flag!(integrator, policy)
    return nothing
end

# The conditions of the events that can end a step before its end: the
# deactivation's threshold margin, the layer times' activation margin and the
# thickness stop's margin.
function _crossings(cp::PreparedProblem, min_thickness::Real)
    layout = cp.def.layout
    cells = cp.deactivation
    out = _ScalarWrapper[]
    cells === nothing || push!(out, _ScalarWrapper(u -> threshold_margin(cells, layout, u)))
    corr = _layer_correction(cp.def)
    if corr !== nothing
        reference = cp.def.setup.reference
        push!(out, _ScalarWrapper(u -> layer_time_margin(corr, layout, reference, u)))
    end
    if is_moving(layout)
        bottom = position_index(layout, 1)
        h = Float64(min_thickness)
        push!(
            out,
            _ScalarWrapper(
                u -> u[position_index(layout, _top_live_node(layout))] - u[bottom] - h,
            ),
        )
    end
    return out
end

# Zero for an admissible state; otherwise the first live cell at or below
# absolute zero, or -1 when the state is inadmissible otherwise.
function _violation(u, layout, setup, abstol)
    inadmissible(u, layout, setup, abstol) || return 0
    cold = below_absolute_zero(u, layout, setup)
    return cold == 0 ? -1 : cold
end

# The per-step observer: the largest constraint residual of the accepted
# states, the end of the solve once the domain test gave up or a fixed step
# left the physical domain, the step cap after a Newton failure, then the
# policy's flag rule. The domain test runs on adaptive steps only, so a
# fixed-step solve records no residual and its accepted states are held to
# admissibility here.
function _observe!(integrator, m::ConstraintMonitor, cap::NewtonCap, policy, violation)
    m.largest = max(m.largest, m.last)
    if !integrator.opts.adaptive
        v = violation(integrator.u, integrator.t)
        if v != 0
            m.failed = true
            m.message =
                "the fixed step accepted at t = $(integrator.t) is inadmissible" *
                (v > 0 ? ": cell $v is at or below absolute zero" : "")
        end
    end
    if m.failed
        @warn "solve(): $(m.message)"
        terminate!(integrator, ReturnCode.Failure)
        return nothing
    end
    cap_step!(integrator, cap)
    apply_observer_flag!(integrator, policy)
    return nothing
end

# A user callback as a read-only observer. The integrator raises the
# derivative-discontinuity flag before each discrete callback's `affect!` and
# folds every callback's own verdict whatever their order, so the policy's
# flag rule runs inside the user's `affect!`, after it.
_observing(callback, policy) = callback
_observing(callback::CallbackSet, policy) = CallbackSet(
    callback.continuous_callbacks...,
    map(c -> _observing(c, policy), callback.discrete_callbacks)...,
)
function _observing(callback::DiscreteCallback, policy)
    affect! = callback.affect!
    return DiscreteCallback(
        callback.condition,
        integrator -> begin
            affect!(integrator)
            apply_observer_flag!(integrator, policy)
            return nothing
        end;
        initialize = callback.initialize,
        finalize = callback.finalize,
        save_positions = Tuple(callback.save_positions),
        initializealg = callback.initializealg,
    )
end

"""
	cap_step!(integrator, cap::NewtonCap)

The Newton-failure step cap, applied after an accepted step: when the step's
nonlinear solves needed retries, the proposed step is held at or below
[`FAILURE_MEMORY_GROWTH`](@ref) times that step for the next
[`FAILURE_MEMORY_STEPS`](@ref FAILURE_MEMORY_GROWTH) accepted steps, then released.
"""
function cap_step!(integrator, cap::NewtonCap)
    failures = integrator.stats.nnonlinconvfail
    if failures > cap.failures
        cap.failures = failures
        cap.cap = FAILURE_MEMORY_GROWTH * (integrator.t - integrator.tprev)
        cap.hold = FAILURE_MEMORY_STEPS
    elseif cap.hold > 0
        cap.hold -= 1
        cap.hold == 0 && (cap.cap = Inf)
    end
    get_proposed_dt(integrator) > cap.cap && set_proposed_dt!(integrator, cap.cap)
    return nothing
end

# The solve's end when a moving mesh thins below `min_thickness`, with the
# outcome `:min_thickness`: the distance from the first node to the outermost
# live one.
function _thickness_termination(
    layout::StateLayout,
    min_thickness::Float64,
    outcome::Base.RefValue{Symbol},
)
    is_moving(layout) || return ()
    bottom = position_index(layout, 1)
    thickness = FunctionWrapper{Float64,Tuple{_StateVector}}(
        u -> u[position_index(layout, _top_live_node(layout))] - u[bottom],
    )
    return (
        ContinuousCallback(
            (u, t, integrator) -> thickness(u) - min_thickness,
            integrator -> begin
                outcome[] = :min_thickness
                terminate!(integrator)
            end;
            save_positions = (false, false),
        ),
    )
end

# The outermost node whose rows are live.
function _top_live_node(layout::StateLayout)
    for nd = n_nodes(layout):-1:1
        layout.mask[position_index(layout, nd)] && return nd
    end
    return 1
end

# The deactivation event: at the first root of the smallest threshold margin
# inside the step, or at the end of a step on which a live cell fell below the
# threshold (the safety net), one pass of `deactivate!` (at a root with the
# threshold widened by the relative tolerance, the mass rows' error scale
# there), then the consistent
# sweep of the algebraic rows, the admissibility of the result (a merge's
# receiver against its own mass tolerances, those the domain test applies to
# the proposals after the event), and the policy's restart step. A
# sweep that leaves a velocity row above round-off or an inadmissible result
# ends the solve with `ReturnCode.Failure`. The margin is tested for downward
# crossings, and location takes a cell at its first crossing, so a cell under
# a correction that would have recovered by the step's end goes there; only
# the end-of-step net, which sees the step's end alone, misses a dip and
# recovery inside one step.
function _deactivation_callback(cp::PreparedProblem, check::DomainCheck, policy)
    cells = cp.deactivation
    cells === nothing && return ()
    def, ws = cp.def, cp.ws
    layout = def.layout
    has_event = FunctionWrapper{Bool,Tuple{_StateVector}}(u -> any_below(cells, layout, u))
    # Deactivates in place with the threshold widened by `slack` and returns
    # `:transferred`, `:burned_out`, `:interior`, `:inconsistent` when the
    # sweep leaves a velocity row above round-off, or `:inadmissible` when the
    # result is not admissible. A void found after a merge ends on the sweep of
    # what was merged.
    deactivate = FunctionWrapper{Symbol,Tuple{_StateVector,Float64,Float64}}(
        (u, t, slack) -> begin
            made = length(cells.events)
            status = deactivate!(u, cells, def, t; slack)
            forget!(check)
            status === :transferred ||
                (status === :interior && length(cells.events) > made) ||
                return status
            _consistent_algebraic!(u, def, ws, t) === nothing || return :inconsistent
            inadmissible(u, layout, def.setup, check.abstol) && return :inadmissible
            _porosity_out(check, _derived(check, u, t), def.setup.gas) &&
                return :inadmissible
            return status
        end,
    )
    # Whether the solve has met an interior void (the located event and the
    # end-of-step net can both see one), and the cell below its threshold there.
    void_seen = FunctionWrapper{Bool,Tuple{}}(() -> cells.outcome[] === :interior)
    void_cell = FunctionWrapper{Int,Tuple{_StateVector,Float64}}(
        (u, slack) -> something(
            findfirst(i -> below_threshold(cells, layout, u, i; slack), 1:n_cells(layout)),
            0,
        ),
    )
    function affect!(integrator, slack::Float64)
        seen = void_seen()
        status = deactivate(integrator.u, integrator.t, slack)
        if status === :inadmissible
            @warn "solve(): the state after deactivating cell(s) at t = $(integrator.t) is inadmissible"
            terminate!(integrator, ReturnCode.Failure)
        elseif status === :inconsistent
            @warn "solve(): no consistent node velocities after deactivating cell(s) at t = $(integrator.t)"
            terminate!(integrator, ReturnCode.Failure)
        elseif status === :interior
            seen ||
                @warn "solve(): cell $(void_cell(integrator.u, slack)) is below its threshold beneath live cells it cannot leave through at t = $(integrator.t); the solve ends with outcome :interior"
            terminate!(integrator)
        elseif status !== :transferred
            terminate!(integrator)
        elseif policy.restart === :forced_dt
            set_proposed_dt!(integrator, policy.restart_dt)
        elseif policy.restart === :heuristic
            auto_dt_reset!(integrator)
            set_proposed_dt!(integrator, integrator.dt)
        end
        return nothing
    end
    safety_net = DiscreteCallback(
        (u, t, integrator) -> has_event(u),
        integrator -> affect!(integrator, 0.0);
        save_positions = (false, false),
    )
    margin = FunctionWrapper{Float64,Tuple{_StateVector}}(
        u -> threshold_margin(cells, layout, u),
    )
    reltol = cp.reltol
    return (
        ContinuousCallback(
            (u, t, integrator) -> margin(u),
            nothing;
            affect_neg! = integrator -> affect!(integrator, reltol),
            rootfind = SciMLBase.RightRootFind,
            save_positions = (false, false),
        ),
        safety_net,
    )
end

# The bound correction of a definition whose layers' equidistribution times
# are set during the solve ([`ActivationTime`](@ref)); `nothing` otherwise.
function _layer_correction(def::ResidualDefinition)
    closure = def.closure
    closure isa BalanceRowClosure || return nothing
    corr = closure.correction
    return estimates_layer_times(corr) ? corr : nothing
end

# The layer times' event, on every moving mesh: at the first root of the
# smallest margin of the layers still waiting for an estimate of their
# equidistribution time (to `r_hi` before the first, to `r_lo` before the
# second), each such layer that has crossed takes its estimate
# ([`estimate_layer_times!`](@ref), the margin widened by the relative
# tolerance as the deactivation's threshold is). Without such a layer the
# margin is infinite and the event never fires; the callback's types are the
# same either way, so the integrator's type does not depend on the correction.
# The state is not modified, and each estimate enters through a switch that is
# zero with zero slope where it is taken (the activation at `r_hi`, the second
# switch at `r_lo`), so the right-hand side does not jump and the node
# velocities stay consistent; the integrator re-evaluates its derivative there
# and the policy's flag rule applies as for an observer.
function _layer_time_callback(cp::PreparedProblem, check::DomainCheck, policy)
    layout = cp.def.layout
    is_moving(layout) || return ()
    corr = _layer_correction(cp.def)
    margin, estimate = if corr === nothing
        (
            FunctionWrapper{Float64,Tuple{_StateVector}}(u -> Inf),
            FunctionWrapper{Nothing,Tuple{_StateVector,Float64}}((u, t) -> nothing),
        )
    else
        reference = cp.def.setup.reference
        reltol = cp.reltol
        (
            FunctionWrapper{Float64,Tuple{_StateVector}}(
                u -> layer_time_margin(corr, layout, reference, u),
            ),
            FunctionWrapper{Nothing,Tuple{_StateVector,Float64}}(
                (u, t) -> begin
                    estimate_layer_times!(corr, layout, reference, u, t; slack = reltol)
                    forget!(check)
                    return nothing
                end,
            ),
        )
    end
    return (
        ContinuousCallback(
            (u, t, integrator) -> margin(u),
            nothing;
            affect_neg! = integrator -> begin
                estimate(integrator.u, integrator.t)
                apply_observer_flag!(integrator, policy)
                return nothing
            end,
            rootfind = SciMLBase.RightRootFind,
            save_positions = (false, false),
        ),
    )
end

# The derived state at `u`, from one residual evaluation into the check's buffer.
function _derived(check::DomainCheck, u, t)
    evaluate!(check, u, t)
    return check.ws.derived
end
