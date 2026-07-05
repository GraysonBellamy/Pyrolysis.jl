"""
	PyrolysisAMGExt

Phase 5 Step C14 — AlgebraicMultigrid.jl linear-solver wiring.

Loaded automatically when a user does `using AlgebraicMultigrid` alongside
`using Pyrolysis`. The extension exposes a smoothed-aggregation AMG
preconditioner constructor via `Pyrolysis.amg_preconditioner` for use as a
preconditioner inside Krylov-based `linsolve` choices (GMRES, BiCGStab) in
the benchmark harness and out-of-tree experiments.
"""
module PyrolysisAMGExt

using Pyrolysis
import AlgebraicMultigrid

"""
	Pyrolysis.amg_preconditioner(A)

Construct a smoothed-aggregation AMG preconditioner from
`AlgebraicMultigrid.jl` for the matrix `A`. Available only when
`using AlgebraicMultigrid` has been loaded; otherwise, calling this function
errors with a load hint.
"""
Pyrolysis.amg_preconditioner(A) =
    AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.smoothed_aggregation(A))

end # module
