"""
	PyrolysisKLUExt

Phase 5 Step C14 — KLU.jl linear-solver wiring.

Loaded automatically when a user does `using KLU` alongside `using Pyrolysis`.
The extension exists so the [`benchmark_jacobian_backends`](@ref) harness (and
out-of-tree consumers of `Pyrolysis.recommended_linear_solver`) can ask for a
KLU-backed `linsolve` without forcing core Pyrolysis to take a hard dependency
on KLU.

`KLUFactorization` is also reachable via `LinearSolve.KLUFactorization` (which
is a hard dep of Pyrolysis through LinearSolve.jl); this extension simply
exposes it under a stable Pyrolysis-level name and confirms that KLU itself is
loaded — useful for the `verify_jacobian` and benchmark CI configurations
where KLU's standalone build options matter.
"""
module PyrolysisKLUExt

using Pyrolysis
using LinearSolve
import KLU

"""
	Pyrolysis.klu_factorization()

Return a `KLUFactorization()` from `LinearSolve.jl`. Available only when
`using KLU` has been loaded; otherwise, calling this function errors with a
load hint.
"""
Pyrolysis.klu_factorization() = LinearSolve.KLUFactorization()

end # module
