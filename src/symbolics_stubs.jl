"""
Symbolic-verification entry points (stubs).

Method-stubs that error out by default. The `PyrolysisSymbolicsExt`
extension (loaded automatically when `Symbolics` is in the user's
environment) attaches concrete methods that derive symbolic derivatives
of mixing rules via `Symbolics.jl`, then compare against the package's
AD path.

Intended for test-suite use only — not loaded by user code.
"""

"""
	verify_mixing_derivatives(material; rtol = 1e-8, atol = 1e-10) -> NamedTuple

For each mixing rule (`PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN` where
applicable) on the supplied `material`, derive the analytical
∂k_eff/∂ξ and ∂k_eff/∂T symbolically via `Symbolics.jl` and compare
against the package's AD path. Returns a `NamedTuple` of
`(rule_name => max_relative_error)`. Test code asserts each entry
satisfies `≤ rtol`.

Backed by `Symbolics.jl` — load `Symbolics` to enable.
"""
function verify_mixing_derivatives end
