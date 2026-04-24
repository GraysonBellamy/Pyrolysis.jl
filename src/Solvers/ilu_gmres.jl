"""
ILU-preconditioned GMRES linear solver for large ALE pyrolysis problems.

Drop-in replacement for `KLUFactorization()` that uses an incomplete LU
preconditioner with GMRES iteration instead of a full sparse LU factorization.

# Why this exists

The ALE mesh velocity formulation creates a dense O(N²) lower-triangular block
in the Jacobian (75% of all non-zeros at 240 cells).  This destroys KLU's sparse
factorization efficiency, causing O(N^2.3) scaling.  ILU(τ) builds an approximate
factorization in O(nnz) time by dropping small fill-in entries, then GMRES
converges in ~10-20 iterations using this preconditioner.

Benchmark results (50s Douglas Fir simulation):
  - 240 cells: ILU-GMRES 14.2s vs KLU 24.4s (1.72× faster)
  - 120 cells: roughly comparable (~4s each)
  -  40 cells: KLU faster (KLU overhead is minimal at small sizes)

# Usage

```julia
solver = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1))
solution = solve(problem; solver = solver, use_analytical_jacobian = true, use_ale = true)
```

# Parameters

- `τ = 0.1`: ILU drop tolerance.  Smaller values → denser ILU → better
  preconditioner but slower setup.  τ=0.1 is optimal for the ALE Jacobian
  structure (see `examples/linear_solver_benchmark.jl` §2 for sweep data).
- `gmres_rtol = 1e-8`: GMRES relative convergence tolerance.
- `restart = 30`: GMRES restart frequency (memory parameter).
"""

import LinearSolve
using LinearSolve: LinearCache, OperatorAssumptions
using IncompleteLU: ilu
# Krylov v0.10 renamed GmresSolver → GmresWorkspace
using Krylov
const KrylovGMRES = isdefined(Krylov, :GmresWorkspace) ? Krylov.GmresWorkspace : Krylov.GmresSolver
const krylov_gmres! = Krylov.gmres!
using SparseArrays: SparseMatrixCSC, getcolptr, rowvals, nonzeros, nnz
using LinearAlgebra: Diagonal, diag

# SciMLBase is re-exported through LinearSolve
const SciMLBase = LinearSolve.SciMLBase

"""
    ILUGMRESFactorization(; τ=0.1, gmres_rtol=1e-8, restart=30)

ILU-preconditioned GMRES as a drop-in replacement for `KLUFactorization()`.

Inherits `LinearSolve.AbstractSparseFactorization` so that OrdinaryDiffEq
materializes the concrete sparse W = I/(γdt) - J matrix (rather than wrapping
it in a matrix-free WOperator).  The `solve!` method then builds an ILU(τ)
preconditioner and solves with Krylov.jl's GMRES.
"""
Base.@kwdef struct ILUGMRESFactorization <: LinearSolve.AbstractSparseFactorization
	τ::Float64 = 0.1
	gmres_rtol::Float64 = 1e-8
	restart::Int = 30
end

function Base.show(io::IO, alg::ILUGMRESFactorization)
	print(io, "ILUGMRESFactorization(τ=", alg.τ, ", rtol=", alg.gmres_rtol,
		", restart=", alg.restart, ")")
end

"""Cache holding the ILU preconditioner and pre-allocated GMRES workspace."""
mutable struct ILUGMRESCacheVal
	P::Any                       # IncompleteLU.ILUFactorization or Diagonal fallback
	krylov_solver::Any           # Krylov GMRES workspace (GmresSolver or GmresWorkspace)
	τ::Float64
	gmres_rtol::Float64
	restart::Int
end

# --- LinearSolve interface methods ---

function LinearSolve.init_cacheval(alg::ILUGMRESFactorization,
		A::SparseMatrixCSC{Float64, Int}, b, u, Pl, Pr,
		maxiters::Int, abstol, reltol, verbose::Union{Bool, LinearSolve.LinearVerbosity},
		assumptions::OperatorAssumptions)
	P = ilu(A, τ = alg.τ)
	# Krylov v0.10: memory is keyword arg; v0.9: positional
	solver = isdefined(Krylov, :GmresWorkspace) ?
		KrylovGMRES(A, b; memory = alg.restart) :
		KrylovGMRES(A, b, alg.restart)
	return ILUGMRESCacheVal(P, solver, alg.τ, alg.gmres_rtol, alg.restart)
end

# Fallback for non-Float64 sparse or non-sparse matrices
function LinearSolve.init_cacheval(alg::ILUGMRESFactorization,
		A, b, u, Pl, Pr,
		maxiters::Int, abstol, reltol, verbose::Union{Bool, LinearSolve.LinearVerbosity},
		assumptions::OperatorAssumptions)
	return nothing
end

function SciMLBase.solve!(cache::LinearCache, alg::ILUGMRESFactorization; kwargs...)
	A = cache.A
	A = convert(AbstractMatrix, A)
	cv = cache.cacheval::ILUGMRESCacheVal

	if cache.isfresh
		# Jacobian changed — rebuild ILU preconditioner
		# Use the matrix directly (A is already a SparseMatrixCSC from jac_prototype)
		try
			cv.P = ilu(A, τ = cv.τ)
		catch
			# ILU can fail on singular/near-singular matrices during transients;
			# fall back to diagonal scaling which still helps GMRES conditioning
			cv.P = Diagonal(1.0 ./ max.(abs.(diag(A)), 1e-12))
		end
		cache.isfresh = false
	end

	# Solve (I/(γdt) - J) x = b  using pre-allocated GMRES with ILU preconditioner
	# gmres! reuses the GmresSolver workspace — no allocations per call
	# Note: memory/restart size is set at GmresSolver construction, not per-call
	gmres!(cv.krylov_solver, A, cache.b;
		M = cv.P, ldiv = true,
		rtol = cv.gmres_rtol,
		restart = true)
	cache.u .= cv.krylov_solver.x

	return SciMLBase.build_linear_solution(alg, cache.u, nothing, cache)
end

# Tell LinearSolve that this solver needs a concrete sparse matrix,
# not a matrix-free WOperator
LinearSolve.needs_concrete_A(::ILUGMRESFactorization) = true

# ============================================================================
# ILU-preconditioned BiCGSTAB
# ============================================================================

const KrylovBiCGSTAB = isdefined(Krylov, :BicgstabWorkspace) ? Krylov.BicgstabWorkspace : Krylov.BicgstabSolver

"""
    ILUBiCGSTABFactorization(; τ=0.1, rtol=1e-8)

ILU-preconditioned BiCGSTAB as a drop-in replacement for `SparspakFactorization()`.
Same approach as `ILUGMRESFactorization` but uses BiCGSTAB instead of GMRES.
BiCGSTAB has fixed memory cost (no restart parameter) and can be faster for
some matrix structures.
"""
Base.@kwdef struct ILUBiCGSTABFactorization <: LinearSolve.AbstractSparseFactorization
	τ::Float64 = 0.1
	rtol::Float64 = 1e-8
end

function Base.show(io::IO, alg::ILUBiCGSTABFactorization)
	print(io, "ILUBiCGSTABFactorization(τ=", alg.τ, ", rtol=", alg.rtol, ")")
end

mutable struct ILUBiCGSTABCacheVal
	P::Any
	krylov_solver::Any
	τ::Float64
	rtol::Float64
end

function LinearSolve.init_cacheval(alg::ILUBiCGSTABFactorization,
		A::SparseMatrixCSC{Float64, Int}, b, u, Pl, Pr,
		maxiters::Int, abstol, reltol, verbose::Union{Bool, LinearSolve.LinearVerbosity},
		assumptions::OperatorAssumptions)
	P = ilu(A, τ = alg.τ)
	solver = isdefined(Krylov, :BicgstabWorkspace) ?
		KrylovBiCGSTAB(A, b) : KrylovBiCGSTAB(A, b)
	return ILUBiCGSTABCacheVal(P, solver, alg.τ, alg.rtol)
end

function LinearSolve.init_cacheval(alg::ILUBiCGSTABFactorization,
		A, b, u, Pl, Pr,
		maxiters::Int, abstol, reltol, verbose::Union{Bool, LinearSolve.LinearVerbosity},
		assumptions::OperatorAssumptions)
	return nothing
end

function SciMLBase.solve!(cache::LinearCache, alg::ILUBiCGSTABFactorization; kwargs...)
	A = cache.A
	A = convert(AbstractMatrix, A)
	cv = cache.cacheval::ILUBiCGSTABCacheVal

	if cache.isfresh
		try
			cv.P = ilu(A, τ = cv.τ)
		catch
			cv.P = Diagonal(1.0 ./ max.(abs.(diag(A)), 1e-12))
		end
		cache.isfresh = false
	end

	Krylov.bicgstab!(cv.krylov_solver, A, cache.b;
		M = cv.P, ldiv = true,
		rtol = cv.rtol)
	cache.u .= cv.krylov_solver.x

	return SciMLBase.build_linear_solution(alg, cache.u, nothing, cache)
end

LinearSolve.needs_concrete_A(::ILUBiCGSTABFactorization) = true
