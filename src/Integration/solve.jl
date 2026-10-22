# ============================================================================
# Solve
# ============================================================================

const _SOLVE_KEYWORDS = (
    :energy_state,
    :abstol,
    :reltol,
    :threshold,
    :maxiters,
    :dt,
    :adaptive,
    :dtmin,
    :dtmax,
    :saveat,
    :callback,
    :progress,
    :validate,
    :use_ale,
    :w_override,
    :min_thickness,
    :radiation_model,
    :jacobian,
    :linear_solver,
    :integrator_policy,
    :sources,
    :deactivation,
    :band,
    :correction,
    :lagrangian,
    :located_events,
)

"""
	solve(problem::PyrolysisProblem; energy_state = :temperature, abstol = 1e-8,
	      reltol = 1e-6, threshold = 1e-2, maxiters = 1_000_000,
	      dt = problem.dt_initial, adaptive = true, dtmin = 1e-12, dtmax = Inf,
	      saveat = problem.output_times, callback = nothing, validate = true,
	      linear_solver = nothing, integrator_policy = nothing,
	      sources = NoSources(), deactivation = false, band = nothing,
	      correction = nothing, lagrangian = nothing,
	      located_events = true) -> CoreSolution

Solve `problem` on the finite-volume core.

The integrator is handed the mass-matrix system with the colored Jacobian, the
per-row absolute tolerances (`abstol` is the temperature tolerance [K]; the mass
and enthalpy rows take theirs from the cell's initial condensed mass and
`threshold`), no initialization of its own, and [`DomainCheck`](@ref) as its
domain test (admissibility, porosity and the constraint residuals of the
algebraic rows). `dt` is the first step, and with `adaptive = false` the fixed
step. `integrator_policy` defaults to [`default_policy`](@ref) with
`linear_solver` ([`default_linear_solver`](@ref) unless given). A per-step
observer records the largest constraint residual of the accepted states, ends
the solve with `ReturnCode.Failure` when the domain test gave up, caps the step
after a nonlinear failure of a Newton-based policy
([`FAILURE_MEMORY_GROWTH`](@ref); a Rosenbrock method has none) and applies the
policy's flag rule. A discrete `callback`, or each discrete callback of a
`CallbackSet`, is a read-only observer: the policy's flag rule follows its
`affect!`, so it leaves the trajectory unchanged under every policy; a
continuous callback's event raises the flag as any event does. `progress` is
accepted for the keyword set of the other solver and draws nothing; `jacobian`
is accepted and ignored, the colored Jacobian being the only one.

`use_ale`, `w_override`, `correction`, `lagrangian`, `deactivation` and
`band` are those of [`setup`](@ref). On a moving mesh the solve ends with return code
`:Terminated` when the thickness falls below `min_thickness` [m] (default
1e-6). With `deactivation`, cells are deactivated in events
([`deactivate!`](@ref)): with `located_events` an event is located inside the
step at the first time a live cell's condensed mass reaches its threshold (the
first root of [`threshold_margin`](@ref)) and applied at the interpolated state
there, the threshold and the band widened by `reltol` (`deactivate!`'s
`slack`): the cells that go are the contiguous run from the exposed side that
the widened tests take, which holds the crossing cell, and a cell within
`reltol` of its threshold inward of a cell above the band stays; at the end
of every step the cells below the
threshold go in one event, the only events without `located_events` and with
it the safety net for a crossing the location does not see. The
solve ends with `:Terminated` when the last cell goes or a cell below the
threshold lies inward of a live one, and with `:Failure` when an event leaves
the state inadmissible (a merge's receiver held to its own mass tolerances,
as the domain test holds every later proposal) or no consistent node
velocities. The state a solve
ends on is always saved, whether or not `saveat` names its time.
"""
function solve(problem::PyrolysisProblem{<:Mesh}; kwargs...)
    for k in keys(kwargs)
        k in _SOLVE_KEYWORDS || throw(ArgumentError("solve(): unknown keyword `$k`"))
    end
    return _solve(problem; kwargs...)
end

function _solve(
    problem::PyrolysisProblem{<:Mesh};
    energy_state::Symbol = :temperature,
    abstol::Real = 1e-8,
    reltol::Real = 1e-6,
    threshold::Real = 1e-2,
    maxiters::Int = 1_000_000,
    dt::Real = problem.dt_initial,
    adaptive::Bool = true,
    dtmin::Real = 1e-12,
    dtmax::Real = Inf,
    saveat = problem.output_times,
    callback = nothing,
    progress::Bool = false,
    validate::Bool = true,
    use_ale::Bool = false,
    w_override = nothing,
    min_thickness::Real = 1e-6,
    radiation_model = NO_RADIATION,
    jacobian = nothing,
    linear_solver = nothing,
    integrator_policy = nothing,
    sources = NoSources(),
    deactivation::Bool = false,
    band::Union{Nothing,Real} = nothing,
    correction::Union{Nothing,Correction} = nothing,
    lagrangian::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    located_events::Bool = true,
)
    jacobian === nothing ||
        @warn "solve(): `jacobian` is ignored by the finite-volume core" maxlog = 1
    validate && validate_problem(problem; strict = false)
    started = time()
    cp = setup(
        problem;
        energy_state,
        reltol,
        abstol,
        threshold,
        use_ale,
        radiation_model,
        w_override,
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
        linear_solver,
        integrator_policy,
        located_events,
        started,
    )
end

"""
	integrate(cp::CoreProblem; u0 = cp.u0, maxiters = 1_000_000,
	          dt = cp.problem.dt_initial, adaptive = true, dtmin = 1e-12,
	          dtmax = Inf, saveat = cp.problem.output_times, callback = nothing,
	          min_thickness = 1e-6, linear_solver = nothing,
	          integrator_policy = nothing, located_events = true,
	          max_constraint_rejections = MAX_CONSTRAINT_REJECTIONS,
	          started = time()) -> CoreSolution

Integrate a set-up problem from the state `u0`, whose algebraic rows must be
consistent ([`consistent_state!`](@ref)), with the tolerances it was set up
with. The keywords are those of [`solve`](@ref); `max_constraint_rejections`
bounds the domain test's run of consecutive constraint rejections
([`DomainCheck`](@ref)), and `started` is the time the reported solve time
counts from. An [`inadmissible`](@ref) `u0` is refused; a fixed-step solve,
which the domain test never sees, ends with `ReturnCode.Failure` at the
first accepted state that is inadmissible. A solve that fails ends on the
last state it accepted. Every call starts with every cell live and, in the inventory
form, the caloric inversion from the set-up temperatures, so repeated
integrations of one problem are bit-identical.
"""
function integrate(
    cp::CoreProblem;
    u0::AbstractVector = cp.u0,
    maxiters::Int = 1_000_000,
    dt::Real = cp.problem.dt_initial,
    adaptive::Bool = true,
    dtmin::Real = 1e-12,
    dtmax::Real = Inf,
    saveat = cp.problem.output_times,
    callback = nothing,
    min_thickness::Real = 1e-6,
    linear_solver = nothing,
    integrator_policy = nothing,
    located_events::Bool = true,
    max_constraint_rejections::Int = MAX_CONSTRAINT_REJECTIONS,
    started::Float64 = time(),
)
    min_thickness > 0 ||
        throw(ArgumentError("min_thickness must be positive, got $min_thickness"))
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    inadmissible(u0, layout, def.setup, cp.abstol) && throw(
        ArgumentError(
            "the initial state is inadmissible: a live component mass below its " *
            "tolerance, a live cell at or below absolute zero, or crossed nodes",
        ),
    )
    cells = cp.deactivation
    cells === nothing || reset_cells!(cells, layout)
    copyto!(ws.T_warm, cp.T0)
    policy =
        integrator_policy === nothing ?
        default_policy(
            linear_solver === nothing ? default_linear_solver() : linear_solver,
        ) : integrator_policy

    rhs! = let def = def, ws = ws
        (du, u, p, t) -> residual!(du, u, def, t, ws)
    end
    J0 = jac_prototype(cj)
    pattern = cj.pattern
    jac! = let cj = cj, pattern = pattern
        (J, u, p, t) -> begin
            same_structure(J, pattern) || error(
                "the integrator's Jacobian lost the colored pattern's structure",
            )
            jacobian!(J, cj, u, t)
            return nothing
        end
    end
    M = all(is_differential(layout)) ? I : mass_matrix(layout)
    f = ODEFunction{true,FullSpecialize}(
        ErasedRHS(rhs!);
        mass_matrix = M,
        jac = ErasedJacobian{typeof(J0)}(jac!),
        jac_prototype = J0,
    )
    ode = ODEProblem(f, Vector{Float64}(u0), cp.problem.tspan)

    # The domain test and every callback reach the problem through
    # fixed-signature wrappers, so the integrator's type does not depend on it.
    check = DomainCheck(cp; max_consecutive = max_constraint_rejections)
    out_of_domain = _DomainWrapper(check)
    monitor = check.monitor
    cap = NewtonCap()
    violation = _ViolationWrapper(let layout = layout, setup = def.setup, abstol = cp.abstol
        u -> _violation(u, layout, setup, abstol)
    end)
    observer = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> _observe!(integrator, monitor, cap, policy, violation);
        save_positions = (false, false),
    )
    # The last accepted state, after every callback of its step.
    accepted = AcceptedState(copy(ode.u0), Float64(cp.problem.tspan[1]))
    recorder = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> _record!(integrator, accepted, policy);
        save_positions = (false, false),
    )
    thin = _thickness_termination(layout, Float64(min_thickness))
    events = _deactivation_callback(cp, check, policy, located_events)
    callbacks =
        callback === nothing ? CallbackSet(observer, events..., thin..., recorder) :
        CallbackSet(observer, events..., thin..., _observing(callback, policy), recorder)

    raw = OrdinaryDiffEq.solve(
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
        initializealg = NoInit(),
        save_end = true,
    )
    monitor.cold == 0 ||
        @warn "solve(): stopped at t = $(raw.t[end]): the last proposed state, at t = $(monitor.cold_t), takes cell $(monitor.cold) to or below absolute zero"
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
        rhs_evaluations = s.nf,
        jacobians = s.njacs,
        linear_solves = s.nsolve,
        nonlinear_failures = s.nnonlinconvfail,
        domain_rejections = check.monitor.rejections,
        largest_constraint_residual = check.monitor.largest,
        invocations = cells === nothing || isempty(cells.events) ? 0 :
                      last(cells.events).invocation,
    )
    return CoreSolution(
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
    )
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
const _ViolationWrapper = FunctionWrapper{Int,Tuple{_StateVector}}

# The last state the integrator accepted and its time, recorded by a
# read-only observer after every other callback of the step.
mutable struct AcceptedState
    u::Vector{Float64}
    t::Float64
end

function _record!(integrator, accepted::AcceptedState, policy)
    copyto!(accepted.u, integrator.u)
    accepted.t = integrator.t
    apply_observer_flag!(integrator, policy)
    return nothing
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
        v = violation(integrator.u)
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
[`FAILURE_MEMORY_STEPS`](@ref) accepted steps, then released.
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

# The solve's end when a moving mesh thins below `min_thickness`: the distance
# from the first node to the outermost live one.
function _thickness_termination(layout::StateLayout, min_thickness::Float64)
    is_moving(layout) || return ()
    bottom = position_index(layout, 1)
    thickness = FunctionWrapper{Float64,Tuple{_StateVector}}(
        u -> u[position_index(layout, _top_live_node(layout))] - u[bottom],
    )
    return (
        ContinuousCallback(
            (u, t, integrator) -> thickness(u) - min_thickness,
            terminate!;
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

# The deactivation event: at the end of every step on which a live cell fell
# below the threshold, or, when `located`, at the first root of the smallest
# threshold margin inside the step (the end-of-step event then the safety
# net), one pass of `deactivate!` (at a root with the threshold widened by the
# relative tolerance, the mass rows' error scale there), then the consistent
# sweep of the algebraic rows, the admissibility of the result (a merge's
# receiver against its own mass tolerances, those the domain test applies to
# the proposals after the event), and the policy's restart step. A
# sweep that leaves a velocity row above round-off or an inadmissible result
# ends the solve with `ReturnCode.Failure`. The margin is tested for downward
# crossings, and location takes a cell at its first crossing, so a cell under
# a correction that would have recovered by the step's end goes there; only
# the end-of-step net, which sees the step's end alone, misses a dip and
# recovery inside one step.
function _deactivation_callback(cp::CoreProblem, check::DomainCheck, policy, located::Bool)
    cells = cp.deactivation
    cells === nothing && return ()
    def, ws = cp.def, cp.ws
    layout = def.layout
    has_event = FunctionWrapper{Bool,Tuple{_StateVector}}(u -> any_below(cells, layout, u))
    # Deactivates in place with the threshold widened by `slack` and returns
    # `:transferred`, `:burned_out`, `:interior`, `:inconsistent` when the
    # sweep leaves a velocity row above round-off, or `:inadmissible` when the
    # result is not admissible.
    deactivate = FunctionWrapper{Symbol,Tuple{_StateVector,Float64,Float64}}(
        (u, t, slack) -> begin
            status = deactivate!(u, cells, def, t; slack)
            status === :transferred || return status
            _consistent_intensities!(u, def, ws, t)
            _consistent_motion!(u, def, ws, t) === nothing || return :inconsistent
            inadmissible(u, layout, def.setup, check.abstol) && return :inadmissible
            _porosity_out(check, _derived(check, u, t), def.setup.gas) &&
                return :inadmissible
            return status
        end,
    )
    function affect!(integrator, slack::Float64)
        status = deactivate(integrator.u, integrator.t, slack)
        if status === :inadmissible
            @warn "solve(): the state after deactivating cell(s) at t = $(integrator.t) is inadmissible"
            terminate!(integrator, ReturnCode.Failure)
        elseif status === :inconsistent
            @warn "solve(): no consistent node velocities after deactivating cell(s) at t = $(integrator.t)"
            terminate!(integrator, ReturnCode.Failure)
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
    located || return (safety_net,)
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

# The derived state at `u`, from one residual evaluation into the check's buffer.
function _derived(check::DomainCheck, u, t)
    residual!(check.du, u, check.def, t, check.ws)
    return check.ws.derived
end
