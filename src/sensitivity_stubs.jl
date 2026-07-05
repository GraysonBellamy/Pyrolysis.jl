"""
Sensitivity-analysis entry points (stubs).

These functions are method-stubs that error out by default. The
`PyrolysisSciMLSensitivityExt` extension (loaded automatically when
`SciMLSensitivity` is in the user's environment) attaches concrete
methods that wire the `PyrolysisProblem` into forward / adjoint
sensitivity workflows.

Calling without the matching `using SciMLSensitivity` raises a clear
`ErrorException` telling the user which package to load.
"""

"""
	forward_sensitivity(base_problem, θ, p_inject; output_fn = sol -> sol,
	                     solve_kwargs...) -> (output, ∂output_∂θ)

Compute the local forward sensitivity `∂output/∂θ` at the parameter
point `θ` for `base_problem`. `p_inject(base_problem, θ) -> problem'`
materialises the perturbed problem; `output_fn(sol) -> y` reduces a
solution to the quantity of interest.

Backed by `SciMLSensitivity.ForwardDiffSensitivity` — load
`SciMLSensitivity` to enable.
"""
function forward_sensitivity end

"""
	adjoint_sensitivity(base_problem, θ, p_inject, loss_fn;
	                    solve_kwargs...) -> (loss, ∂loss_∂θ)

Compute the adjoint gradient `∂loss/∂θ` at the parameter point `θ` for
`base_problem`. `p_inject(base_problem, θ)` materialises the perturbed
problem; `loss_fn(sol) -> ℝ` reduces a solution to a scalar loss.

Backed by `SciMLSensitivity` adjoint methods. Pass a `sensealg` through
`solve_kwargs...` when you need a specific reverse-mode algorithm. Load
`SciMLSensitivity` to enable.
"""
function adjoint_sensitivity end
