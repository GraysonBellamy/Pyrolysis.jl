"""
PyrolysisSciMLSensitivityExt — forward / adjoint sensitivities for
`PyrolysisProblem` via SciMLSensitivity.jl.

Loaded automatically when both `Pyrolysis` and `SciMLSensitivity` are
available in the user's environment. Attaches concrete methods to the
`Pyrolysis.forward_sensitivity` and `Pyrolysis.adjoint_sensitivity`
stubs.

The strategy is wrapper-based rather than residual-rewriting: `θ` flows
into a fresh `PyrolysisProblem` via the user-supplied `p_inject`
closure, the wrapper invokes `solve(modified_problem)` end-to-end, and
ForwardDiff (forward mode) or Zygote plus a SciMLSensitivity `sensealg`
(reverse mode) differentiates the chain. This avoids invasive refactor
of the residual signature while still giving users access to local
sensitivity information for the inverse-analysis workflow under
`examples/inverse_analysis/`.
"""
module PyrolysisSciMLSensitivityExt

using Pyrolysis
using Pyrolysis: forward_sensitivity, adjoint_sensitivity
using SciMLSensitivity
using ForwardDiff
using Zygote

# ---------------------------------------------------------------------------
# Forward-mode sensitivity (via ForwardDiff)
# ---------------------------------------------------------------------------
function Pyrolysis.forward_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function;
    output_fn::Function = identity,
    solve_kwargs...,
)

    # Wrap solve so ForwardDiff sees `θ -> output_fn(solve(p_inject(prob, θ)))`.
    wrapped = θ -> begin
        prob = p_inject(base_problem, θ)
        sol = solve(prob; progress = false, solve_kwargs...)
        return output_fn(sol)
    end

    output = wrapped(θ)
    ∂out_∂θ = ForwardDiff.jacobian(wrapped, θ)
    return (output, ∂out_∂θ)
end

# ---------------------------------------------------------------------------
# Reverse-mode (adjoint) sensitivity (via Zygote + InterpolatingAdjoint)
# ---------------------------------------------------------------------------
function Pyrolysis.adjoint_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function,
    loss_fn::Function;
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
    solve_kwargs...,
)

    wrapped =
        θ -> begin
            prob = p_inject(base_problem, θ)
            # `solve` forwards `sensealg` to OrdinaryDiffEq.solve; keeping it as a
            # first-class kwarg avoids relaxing solve's unknown-kwarg guard.
            sol = solve(prob; progress = false, sensealg = sensealg, solve_kwargs...)
            return loss_fn(sol)
        end

    loss = wrapped(θ)
    grad = first(Zygote.gradient(wrapped, θ))
    return (loss, grad)
end

end # module
