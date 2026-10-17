module RestartProbe

using LinearAlgebra, SparseArrays, ForwardDiff
using OrdinaryDiffEq, OrdinaryDiffEqSDIRK
import LinearSolve, SciMLBase
import ..SolverProbe

export JumpLadder,
    event_state,
    structural_pattern,
    column_colors,
    ColorCache,
    colored_jac!,
    run_restarts,
    trajectory_error,
    jump_error,
    jacobian_cost,
    callback_control

# Same diffusion/recurrence equations as SolverProbe. For J events the jump is
# y+ = q*y- + (-1)^(k+1)*0.05*sin(n*pi*x), q = 0.6^(1/J).
# Thus total baseline contraction stays 0.6 rather than 0.6^J. Every event
# excites a bounded fast transient; this is not a production deactivation map.
# Event count changes both restart frequency and the amount of transient work.
mutable struct JumpLadder
    n::Int
    diffusion::Float64
    beta::Float64
    times::Vector{Float64}
    scales::Vector{Float64}
    amplitudes::Vector{Float64}
    slow_mode::Vector{Float64}
    fast_mode::Vector{Float64}
    segment::Int
    jumps::Int
    jac_derivative::Vector{Float64}
    jac_revision::Int
end

function JumpLadder(n, events; schedule = :uniform)
    schedule in (:uniform, :clustered) || error("unknown schedule")
    times = collect((1:events) ./ (events + 1))
    schedule == :clustered && (times .= 0.45 .+ 0.1 .* times)
    q = events == 0 ? 1.0 : 0.6^(1 / events)
    scales, amplitudes = ones(events + 1), zeros(events + 1)
    amplitudes[1] = 0.1
    for k = 1:events
        elapsed = times[k] - (k == 1 ? 0.0 : times[k-1])
        scales[k+1] = q * scales[k]
        amplitudes[k+1] = q * amplitudes[k] * exp(-200elapsed) + 0.05 * (-1)^(k+1)
    end
    x = (1:n) ./ (n + 1)
    return JumpLadder(
        n,
        0.25 * (n+1)^2,
        0.4,
        times,
        scales,
        amplitudes,
        sinpi.(x),
        sinpi.(n .* x),
        0,
        0,
        zeros(n),
        0,
    )
end

function event_state(p, t, segment)
    start = segment == 0 ? 0.0 : p.times[segment]
    slow = 0.2 .* p.slow_mode .* exp(-t)
    fast = p.amplitudes[segment+1] .* p.fast_mode .* exp(-200(t-start))
    scale = p.scales[segment+1]
    return scale .* (1 .+ slow) .+ fast, -scale .* slow .- 200 .* fast
end

# Explicitly choose the active segment in the RHS: the last stage at an event
# uses the left limit; only the event callback switches to the right limit.
SolverProbe.exact_state(p::JumpLadder, t) = event_state(p, t, p.segment)

# Connectivity pattern, independent of the state and of numerical zeros.
function structural_pattern(n, lifted)
    rows, cols = Int[], Int[]
    function add(i, j)
        push!(rows, i)
        push!(cols, j)
    end
    if lifted
        for i = 1:n
            for j = max(1, i-1):min(n, i+1)
                add(2i-1, 2j-1)
            end
            add(2i-1, 2i)
            add(2i, 2i-1)
            add(2i, 2i)
            i > 1 && add(2i, 2i-2)
        end
    else
        for i = 1:n, j = 1:min(n, i+1)
            add(i, j)
        end
    end
    m = lifted ? 2n : n
    return sparse(rows, cols, ones(length(rows)), m, m)
end

# A deterministic greedy coloring suffices for this prototype; no tracing or
# production color-count claim. One ForwardDiff directional pass per color.
function column_colors(pattern)
    occupied = [BitSet() for _ in axes(pattern, 1)]
    colors = zeros(Int, size(pattern, 2))
    for j in axes(pattern, 2)
        forbidden = BitSet()
        for k in nzrange(pattern, j)
            union!(forbidden, occupied[rowvals(pattern)[k]])
        end
        color = 1
        while color in forbidden
            color += 1
        end
        colors[j] = color
        for k in nzrange(pattern, j)
            push!(occupied[rowvals(pattern)[k]], color)
        end
    end
    return colors
end

struct ColorTag end
const Directional = ForwardDiff.Dual{ColorTag,Float64,1}
mutable struct ColorCache
    pattern::SparseMatrixCSC{Float64,Int}
    colors::Vector{Int}
    columns::Vector{Vector{Int}}
    input::Vector{Directional}
    output::Vector{Directional}
    calls::Int
    passes::Int
end
function ColorCache(n, lifted)
    pattern = structural_pattern(n, lifted)
    colors = column_colors(pattern)
    columns = [findall(==(c), colors) for c = 1:maximum(colors)]
    input = zeros(Directional, length(colors))
    return ColorCache(pattern, colors, columns, input, similar(input), 0, 0)
end

function colored_jac!(J, u, p, t, lifted, cache)
    fill!(J, 0)
    for (color, columns) in enumerate(cache.columns)
        for j in eachindex(u)
            cache.input[j] =
                ForwardDiff.Dual{ColorTag}(u[j], Float64(cache.colors[j] == color))
        end
        SolverProbe.ladder_rhs!(cache.output, cache.input, p, t, lifted)
        cache.passes += 1
        for j in columns, k in nzrange(cache.pattern, j)
            i = rowvals(cache.pattern)[k]
            J[i, j] = ForwardDiff.partials(cache.output[i])[1]
        end
    end
    # The exact bordered factorization also needs the recurrence derivative.
    y = lifted ? view(u, 1:2:length(u)) : u
    p.jac_derivative .= 2 .* y ./ p.n
    p.jac_revision += 1
    cache.calls += 1
    return nothing
end

function make_jump_problem(n, events, form, backend; schedule = :uniform)
    form in (:lifted, :bordered, :dense) || error("unknown form")
    backend in (:analytic, :colored) || error("unknown Jacobian backend")
    p = JumpLadder(n, events; schedule)
    lifted = form == :lifted
    y, _ = event_state(p, 0.0, 0)
    u = lifted ? vec(permutedims(hcat(y, SolverProbe.recurrence(y)))) : y
    colors = ColorCache(n, lifted)
    prototype = lifted ? copy(colors.pattern) : zeros(n, n)
    rhs!(du, u, p, t) = SolverProbe.ladder_rhs!(du, u, p, t, lifted)
    function jac!(J, u, p, t)
        if backend == :colored
            colored_jac!(J, u, p, t, lifted, colors)
        else
            SolverProbe.ladder_jac!(J, u, p, t, lifted)
        end
        return nothing
    end
    mass = lifted ? Diagonal(repeat([1.0, 0.0], n)) : I
    f = ODEFunction{true,SciMLBase.FullSpecialize}(
        rhs!;
        jac = jac!,
        jac_prototype = prototype,
        mass_matrix = mass,
    )
    return ODEProblem(f, u, (0.0, 1.0), p), colors
end

function run_restarts(
    n,
    method,
    form,
    backend,
    events;
    tolerance = 1e-6,
    schedule = :uniform,
    verify = false,
    observer = :clean,
    reset_dt = false,
)
    prob, colors = make_jump_problem(n, events, form, backend; schedule)
    p, lifted = prob.p, form == :lifted
    border = SolverProbe.BorderCache(n; verify)
    custom(A, b, u, parameters, fresh, Pl, Pr, val; kwargs...) =
        SolverProbe.border_solve!(border, A, b, u, p)
    linear =
        form == :bordered ? LinearSolve.LinearSolveFunction(custom) :
        lifted ? LinearSolve.UMFPACKFactorization() : LinearSolve.LHLFactorization()
    method in (:KenCarp4_raw, :FBDF, :Rodas5P) || error("unknown method")
    alg =
        method == :KenCarp4_raw ?
        KenCarp4(; linsolve = linear, concrete_jac = true, smooth_est = false) :
        method == :FBDF ? FBDF(; linsolve = linear, concrete_jac = true) :
        Rodas5P(; linsolve = linear, concrete_jac = true)
    order_samples = Int[]
    restart_orders = Int[]
    pending_restart = Ref(false)
    event_constraints = Float64[]
    event_times = Float64[]
    function record!(integrator)
        if method == :FBDF
            # cache.order is the BDF formula order; filter_order can differ in
            # versions with time filtering. This probe uses default FBDF.
            order = integrator.cache.order
            push!(order_samples, order)
            pending_restart[] && push!(restart_orders, order)
        end
        pending_restart[] = false
        derivative_discontinuity!(integrator, observer == :dirty)
        return nothing
    end
    function jump!(integrator)
        state = integrator.u
        y = lifted ? view(state, 1:2:2n) : state
        k = p.segment + 1
        q = p.scales[k+1] / p.scales[k]
        y .*= q
        y .+= (0.05 * (-1)^(k+1)) .* p.fast_mode
        lifted && (state[2:2:2n] .= SolverProbe.recurrence(y))
        p.segment = k
        p.jumps += 1
        push!(event_times, integrator.t)
        if lifted
            defect =
                maximum(abs(state[2i] - (i > 1 ? state[2i-2] : 0.0) - y[i]^2/n) for i = 1:n)
            push!(event_constraints, defect)
        end
        pending_restart[] = true
        derivative_discontinuity!(integrator, true)
        # Sensitivity only: 1e-6 in this nondimensional toy is not the
        # dimensional DT_AFTER_PARK=1e-6 s on the ~500 s production deck.
        reset_dt && set_proposed_dt!(integrator, 1e-6)
        return nothing
    end
    observe = DiscreteCallback((u, t, i) -> true, record!; save_positions = (false, false))
    jump = DiscreteCallback(
        (u, t, i) -> p.segment < length(p.times) && t == p.times[p.segment+1],
        jump!;
        save_positions = (true, true),
    )
    # Observe the accepted formula order before the jump changes history. A
    # clean observer explicitly clears its own flag and never changes state.
    callbacks = observer == :none ? CallbackSet(jump) : CallbackSet(observe, jump)
    sol = solve(
        prob,
        alg;
        callback = callbacks,
        tstops = p.times,
        reltol = tolerance,
        abstol = 0.01tolerance,
        initializealg = SciMLBase.NoInit(),
        internalnorm = (u, t) -> norm(u) / sqrt(length(u)),
        maxiters = 100000,
    )
    return (;
        sol,
        p,
        colors,
        border,
        order_samples,
        restart_orders,
        event_constraints,
        event_times,
    )
end

function trajectory_error(result, form)
    p, sol = result.p, result.sol
    worst = 0.0
    for (t, u) in zip(sol.t, sol.u)
        y = form == :lifted ? view(u, 1:2:length(u)) : u
        k = searchsortedlast(p.times, t)
        err = norm(y - event_state(p, t, k)[1], Inf)
        if k > 0 && t == p.times[k]
            err = min(err, norm(y - event_state(p, t, k-1)[1], Inf))
        end
        worst = max(worst, err)
    end
    return worst
end

# Separately require BOTH saved one-sided states and their actual jump map;
# otherwise the min(left,right) convention above could hide a missed jump.
function jump_error(result, form; map_only = false)
    worst = 0.0
    for (k, t) in enumerate(result.p.times)
        indices = findall(==(t), result.sol.t)
        length(indices) >= 2 || return Inf
        left, right = result.sol.u[first(indices)], result.sol.u[last(indices)]
        yl = form == :lifted ? left[1:2:end] : left
        yr = form == :lifted ? right[1:2:end] : right
        p = result.p
        q = p.scales[k+1] / p.scales[k]
        expected = q .* yl .+ (0.05 * (-1)^(k+1)) .* p.fast_mode
        worst = max(worst, norm(yr - expected, Inf))
        if !map_only
            worst = max(
                worst,
                norm(yl - event_state(p, t, k-1)[1], Inf),
                norm(yr - event_state(p, t, k)[1], Inf),
            )
        end
    end
    return worst
end

function jacobian_cost(n, form; repetitions = 200)
    prob, colors = make_jump_problem(n, 10, form, :colored)
    p, u, lifted = prob.p, prob.u0, form == :lifted
    J = copy(prob.f.jac_prototype)
    du = similar(u)
    # Timing small kernels in batches avoids timer resolution effects. The
    # median includes allocations, excludes compilation/preparation, and does
    # not insert a timer into each integrator Jacobian call.
    function batch(kind)
        for _ = 1:repetitions
            if kind == :rhs
                SolverProbe.ladder_rhs!(du, u, p, 0.0, lifted)
            elseif kind == :analytic
                SolverProbe.ladder_jac!(J, u, p, 0.0, lifted)
            else
                colored_jac!(J, u, p, 0.0, lifted, colors)
            end
        end
        return nothing
    end
    times = Float64[]
    for kind in (:rhs, :analytic, :colored)
        batch(kind)
        samples = [@elapsed(batch(kind)) for _ = 1:5]
        push!(times, sort(samples)[3] / repetitions)
    end
    return (;
        n,
        form,
        colors = length(colors.columns),
        rhs_seconds = times[1],
        analytic_seconds = times[2],
        colored_seconds = times[3],
        colored_rhs_equivalents = times[3] / times[1],
    )
end

# Exercise callback OR-folding independently of the ladder and record the
# actual FBDF order on a smooth scalar ODE. The last callback clearing its
# own flag must not erase a real discontinuity reported by an earlier one.
function callback_control(mode)
    orders = Int[]
    prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))
    function observe!(i)
        push!(orders, i.cache.order)
        derivative_discontinuity!(i, false)
    end
    clean = DiscreteCallback((u, t, i) -> true, observe!; save_positions = (false, false))
    dirty = DiscreteCallback(
        (u, t, i) -> true,
        i -> derivative_discontinuity!(i, true);
        save_positions = (false, false),
    )
    callback =
        mode == :clean ? CallbackSet(clean) :
        mode == :dirty_then_clean ? CallbackSet(dirty, clean) :
        mode == :clean_then_dirty ? CallbackSet(clean, dirty) : error("unknown control")
    # Bound the intentionally broken control. Reinitializing on every step
    # can exhaust the budget before t=1; that is a measured failure, not a
    # reason to grant this configuration a looser accuracy tolerance.
    sol = solve(
        prob,
        FBDF();
        callback,
        reltol = 1e-6,
        abstol = 1e-8,
        maxiters = 10000,
        verbose = DEVerbosity(SciMLLogging.None()),
    )
    return (; sol, orders)
end

end
