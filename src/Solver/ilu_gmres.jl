"""
ILU-preconditioned GMRES linear solver for large ALE pyrolysis problems.

Drop-in replacement for a sparse direct factorization that uses an incomplete LU
preconditioner with GMRES iteration instead of a full sparse LU factorization.

# Why this exists

ALE structured couplings can make direct sparse factorizations expensive on
large meshes. ILU(τ) builds an approximate factorization by dropping small
fill-in entries, then GMRES converges against that preconditioner. Treat this
as an experimental option: benchmark end-to-end solves against the default
direct factorization before using it for production studies.

# Usage

```julia
alg = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1))
solution = solve(problem; integrator = alg,
    jacobian = Structured(strategy = DirectSolve()),
    use_ale = true)
```

# Parameters

- `τ = 0.1`: ILU drop tolerance.  Smaller values → denser ILU → better
  preconditioner but slower setup.
- `gmres_rtol = 1e-8`: GMRES relative convergence tolerance.
- `restart = 30`: GMRES restart frequency (memory parameter).
"""

import LinearSolve
using LinearSolve: LinearCache, OperatorAssumptions
using IncompleteLU: ilu
using Krylov
using SparseArrays: SparseMatrixCSC, nnz
using LinearAlgebra: Diagonal, diag

# SciMLBase is re-exported through LinearSolve
const SciMLBase = LinearSolve.SciMLBase

# The Krylov iterate as a linear solution: `Success` when the iteration met its
# tolerance, `ConvergenceFailure` otherwise, which the implicit integrator's
# Newton treats as a divergent iteration (the step is retried smaller).
function _krylov_solution(alg, cache::LinearCache, workspace)
    cache.u .= workspace.x
    retcode =
        Krylov.issolved(workspace) ? SciMLBase.ReturnCode.Success :
        SciMLBase.ReturnCode.ConvergenceFailure
    return SciMLBase.build_linear_solution(
        alg,
        cache.u,
        nothing,
        cache;
        retcode = retcode,
        iters = workspace.stats.niter,
    )
end

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
    print(
        io,
        "ILUGMRESFactorization(τ=",
        alg.τ,
        ", rtol=",
        alg.gmres_rtol,
        ", restart=",
        alg.restart,
        ")",
    )
end

"""Cache holding the ILU preconditioner and pre-allocated GMRES workspace."""
mutable struct ILUGMRESCacheVal
    P::Any                       # IncompleteLU.ILUFactorization or Diagonal fallback
    krylov_solver::Any           # Krylov.GmresWorkspace
    τ::Float64
    gmres_rtol::Float64
end

# --- LinearSolve interface methods ---

function LinearSolve.init_cacheval(
    alg::ILUGMRESFactorization,
    A::SparseMatrixCSC{Float64,Int},
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{Bool,LinearSolve.LinearVerbosity},
    assumptions::OperatorAssumptions,
)
    P = ilu(A, τ = alg.τ)
    solver = Krylov.GmresWorkspace(A, b; memory = alg.restart)
    return ILUGMRESCacheVal(P, solver, alg.τ, alg.gmres_rtol)
end

# Fallback for non-Float64 sparse or non-sparse matrices
function LinearSolve.init_cacheval(
    alg::ILUGMRESFactorization,
    A,
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{Bool,LinearSolve.LinearVerbosity},
    assumptions::OperatorAssumptions,
)
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

    # Solve (I/(γdt) - J) x = b with the pre-allocated GMRES workspace and the
    # ILU preconditioner; the restart memory is fixed at workspace construction.
    gmres!(
        cv.krylov_solver,
        A,
        cache.b;
        M = cv.P,
        ldiv = true,
        rtol = cv.gmres_rtol,
        restart = true,
    )
    return _krylov_solution(alg, cache, cv.krylov_solver)
end

# Tell LinearSolve that this solver needs a concrete sparse matrix,
# not a matrix-free WOperator
LinearSolve.needs_concrete_A(::ILUGMRESFactorization) = true

# ============================================================================
# ILU-preconditioned BiCGSTAB
# ============================================================================

"""
    ILUBiCGSTABFactorization(; τ=0.1, rtol=1e-8)

ILU-preconditioned BiCGSTAB as a drop-in replacement for a sparse direct factorization.
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

function LinearSolve.init_cacheval(
    alg::ILUBiCGSTABFactorization,
    A::SparseMatrixCSC{Float64,Int},
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{Bool,LinearSolve.LinearVerbosity},
    assumptions::OperatorAssumptions,
)
    P = ilu(A, τ = alg.τ)
    solver = Krylov.BicgstabWorkspace(A, b)
    return ILUBiCGSTABCacheVal(P, solver, alg.τ, alg.rtol)
end

function LinearSolve.init_cacheval(
    alg::ILUBiCGSTABFactorization,
    A,
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{Bool,LinearSolve.LinearVerbosity},
    assumptions::OperatorAssumptions,
)
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

    Krylov.bicgstab!(cv.krylov_solver, A, cache.b; M = cv.P, ldiv = true, rtol = cv.rtol)
    return _krylov_solution(alg, cache, cv.krylov_solver)
end

LinearSolve.needs_concrete_A(::ILUBiCGSTABFactorization) = true
