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
	default_policy(linear_solver)

`KenCarp4` with the raw (unsmoothed) error estimate, the observer flag kept and
the solver-chosen restart.
"""
default_policy(linsolve) = IntegratorPolicy(
    KenCarp4(autodiff = AutoFiniteDiff(), linsolve = linsolve, smooth_est = false),
    :solver_chosen,
    :keep,
    1e-6,
)

"""
	apply_observer_flag!(integrator, policy)

Leave or clear the derivative-discontinuity flag after a per-step observer,
as the policy says.
"""
@inline apply_observer_flag!(integrator, policy::IntegratorPolicy) =
    policy.observer_flag === :clear ? derivative_discontinuity!(integrator, false) : nothing
