"""
	PyrolysisPardisoExt

Phase 5 Step C14 — Pardiso.jl linear-solver wiring.

Loaded automatically when a user does `using Pardiso` alongside
`using Pyrolysis`. The extension exposes Pardiso's MKL-backed direct solver
via `Pyrolysis.pardiso_factorization()` so the benchmark harness (and
production users with an MKL license) can route the implicit linear solve
through Pardiso without making it a hard dependency of core Pyrolysis.
"""
module PyrolysisPardisoExt

using Pyrolysis
using LinearSolve: MKLPardisoFactorize
import Pardiso

"""
	Pyrolysis.pardiso_factorization()

Return an `MKLPardisoFactorize()` instance from `LinearSolve.jl`. Available
only when `using Pardiso` has been loaded; otherwise, calling this function
errors with a load hint.
"""
Pyrolysis.pardiso_factorization() = MKLPardisoFactorize()

end # module
