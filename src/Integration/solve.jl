# ============================================================================
# Admissibility
# ============================================================================

"""
	inadmissible(u, layout, abstol) -> Bool

Whether the state `u` lies outside the physical domain: a component mass more
negative than its absolute tolerance. Read from the raw state, never through a
guard. The integrator rejects such a step and halves it.
"""
function inadmissible(
    u::AbstractVector,
    layout::StateLayout{NC},
    abstol::Vector{Float64},
) where {NC}
    @inbounds for i = 1:n_cells(layout)
        layout.mask[energy_index(layout, i)] || continue
        for j = 1:NC
            k = mass_index(layout, i, j)
            u[k] < -abstol[k] && return true
        end
    end
    return false
end

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
    :radiation_model,
    :jacobian,
    :linear_solver,
    :integrator_policy,
    :sources,
)

"""
	solve(problem::PyrolysisProblem; energy_state = :temperature, abstol = 1e-8,
	      reltol = 1e-6, threshold = 1e-2, maxiters = 1_000_000,
	      dt = problem.dt_initial, adaptive = true, dtmin = 1e-12, dtmax = Inf,
	      saveat = problem.output_times, callback = nothing, validate = true,
	      linear_solver = nothing, integrator_policy = nothing,
	      sources = NoSources()) -> CoreSolution

Solve `problem` on the finite-volume core.

The integrator is handed the mass-matrix system with the colored Jacobian, the
per-row absolute tolerances (`abstol` is the temperature tolerance [K]; the mass
and enthalpy rows take theirs from the cell's initial condensed mass and
`threshold`), no initialization of its own, and the admissibility check as its
domain test. `dt` is the first step, and with `adaptive = false` the fixed
step. `integrator_policy` defaults to [`default_policy`](@ref) with
`linear_solver` (KLU unless given); a per-step observer applies the policy's
flag rule. `progress` is accepted for the keyword set of the other solver and
draws nothing; `jacobian` is accepted and ignored, the colored Jacobian being the
only one.
"""
function solve(problem::PyrolysisProblem; kwargs...)
    for k in keys(kwargs)
        k in _SOLVE_KEYWORDS || throw(ArgumentError("solve(): unknown keyword `$k`"))
    end
    return _solve(problem; kwargs...)
end

function _solve(
    problem::PyrolysisProblem;
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
    radiation_model = NO_RADIATION,
    jacobian = nothing,
    linear_solver = nothing,
    integrator_policy = nothing,
    sources = NoSources(),
)
    jacobian === nothing ||
        @warn "solve(): `jacobian` is ignored by the finite-volume core" maxlog = 1
    validate && validate_problem(problem; strict = false)
    start = time()
    cp = setup(
        problem;
        energy_state,
        reltol,
        abstol,
        threshold,
        use_ale,
        radiation_model,
        sources,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    policy =
        integrator_policy === nothing ?
        default_policy(linear_solver === nothing ? KLUFactorization() : linear_solver) :
        integrator_policy

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
    ode = ODEProblem(f, cp.u0, cp.problem.tspan)

    tol = cp.abstol
    out_of_domain = let layout = layout, tol = tol
        (u, p, t) -> inadmissible(u, layout, tol)
    end
    observer = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> apply_observer_flag!(integrator, policy);
        save_positions = (false, false),
    )
    callbacks =
        callback === nothing ? CallbackSet(observer) : CallbackSet(observer, callback)

    raw = OrdinaryDiffEq.solve(
        ode,
        policy.alg;
        abstol = tol,
        reltol = Float64(reltol),
        maxiters,
        dtmin = Float64(dtmin),
        dtmax = Float64(dtmax),
        dt = Float64(dt),
        adaptive,
        saveat = saveat isa Number ? Float64(saveat) : collect(Float64, saveat),
        callback = callbacks,
        isoutofdomain = out_of_domain,
        initializealg = NoInit(),
    )
    elapsed = time() - start
    s = raw.stats
    stats = SolverStatistics(
        accepted = s.naccept,
        rejected = s.nreject,
        rhs_evaluations = s.nf,
        jacobians = s.njacs,
        linear_solves = s.nsolve,
        nonlinear_failures = s.nnonlinconvfail,
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
        raw,
    )
end
