# ============================================================================
# Integrator policy
# ============================================================================

"""
	IntegratorPolicy{Alg}

An integrator together with the two policies that are part of its contract:

- `restart`: the step size after an event that modifies the state —
  `:solver_chosen` keeps the integrator's own proposal, `:forced_dt` proposes
  `restart_dt`;
- `observer_flag`: what a per-step observer does with the integrator's
  derivative-discontinuity flag, which the callback machinery raises before
  every callback — `:keep` leaves it raised (the integrator re-evaluates the
  derivative at the accepted state), `:clear` clears it (a read-only observer
  is then numerically inert).
"""
struct IntegratorPolicy{Alg}
    alg::Alg
    restart::Symbol
    observer_flag::Symbol
    restart_dt::Float64
    function IntegratorPolicy(
        alg::Alg,
        restart::Symbol,
        observer_flag::Symbol,
        restart_dt::Real,
    ) where {Alg}
        restart in (:solver_chosen, :forced_dt) ||
            throw(ArgumentError("restart must be :solver_chosen or :forced_dt"))
        observer_flag in (:keep, :clear) ||
            throw(ArgumentError("observer_flag must be :keep or :clear"))
        return new{Alg}(alg, restart, observer_flag, Float64(restart_dt))
    end
end

"""
	DT_AFTER_EVENT

The restart step [s] a `:forced_dt` policy proposes after an event that
modifies the state.
"""
const DT_AFTER_EVENT = 1e-6

"""
	kencarp4_policy(linear_solver)

`KenCarp4` with the raw (unsmoothed) error estimate, the observer flag kept and
the solver-chosen restart. The Jacobian is the colored one; the time
derivative of the right-hand side, where the method needs it, is a finite
difference.
"""
kencarp4_policy(linsolve) = IntegratorPolicy(
    KenCarp4(autodiff = AutoFiniteDiff(), linsolve = linsolve, smooth_est = false),
    :solver_chosen,
    :keep,
    DT_AFTER_EVENT,
)

"""
	rodas5p_policy(linear_solver; restart = :solver_chosen)

`Rodas5P`, a Rosenbrock method that evaluates the exact colored Jacobian at
every step and the time derivative of the right-hand side by finite
differences, with either restart. It rebuilds its Jacobian every step, so the
observer flag changes nothing; the policy clears it.
"""
rodas5p_policy(linsolve; restart::Symbol = :solver_chosen) = IntegratorPolicy(
    Rodas5P(autodiff = AutoFiniteDiff(), linsolve = linsolve),
    restart,
    :clear,
    DT_AFTER_EVENT,
)

"""
	fbdf_policy(linear_solver)

`FBDF`, the variable-order backward-differentiation method, with the observer
flag cleared (with it raised the method restarts at order one on every step)
and a forced small restart step after an event, where it restarts at order
one.
"""
fbdf_policy(linsolve) = IntegratorPolicy(
    FBDF(autodiff = AutoFiniteDiff(), linsolve = linsolve),
    :forced_dt,
    :clear,
    DT_AFTER_EVENT,
)

"""
	default_policy(linear_solver)

The default integrator: [`rodas5p_policy`](@ref) with the solver-chosen
restart.
"""
default_policy(linsolve) = rodas5p_policy(linsolve)

"""
	apply_observer_flag!(integrator, policy)

Leave or clear the derivative-discontinuity flag after a per-step observer,
as the policy says.
"""
@inline apply_observer_flag!(integrator, policy::IntegratorPolicy) =
    policy.observer_flag === :clear ? derivative_discontinuity!(integrator, false) : nothing
