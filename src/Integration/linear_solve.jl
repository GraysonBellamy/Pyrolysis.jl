# ============================================================================
# The linear solve of the integrator's W on a mass-matrix system
# ============================================================================

"""
	StructuredKLU(inner = KLUFactorization())

KLU for the integrator's `W = J − M/(γΔt)` of a mass-matrix system, with one
symbolic analysis for the whole solve. For a mass matrix that is not a
multiple of the identity, OrdinaryDiffEq forms `W` by a sparse broadcast that
stores only the entries that compute to nonzero, so `W`'s structure changes
whenever an entry of the Jacobian passes through zero, and KLU, which reuses
its analysis only on an unchanged structure, analyses afresh at every such
step. `W`'s stored entries are a subset of `jac_prototype`'s (the prototype
holds the diagonal), and the integrator hands `init` the all-ones `W` on the
prototype's whole structure: the cache keeps that matrix, scatters each fresh
`W` into it, zero elsewhere, and the inner KLU refactors it numerically on the
unchanged structure. An entry of `W` outside that structure is an error.
"""
struct StructuredKLU <: LinearSolve.AbstractSparseFactorization
    inner::KLUFactorization
end

StructuredKLU() = StructuredKLU(KLUFactorization())

# The prototype-structured matrix and the inner KLU cache that factors it.
struct StructuredKLUCache{C}
    full::SparseMatrixCSC{Float64,Int}
    inner::C
end

function LinearSolve.init_cacheval(
    alg::StructuredKLU,
    A::SparseMatrixCSC{Float64,Int},
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{Bool,LinearSolve.LinearVerbosity},
    assumptions::LinearSolve.OperatorAssumptions,
)
    full = copy(A)
    inner = LinearSolve.init(LinearSolve.LinearProblem(full, zeros(length(b))), alg.inner)
    return StructuredKLUCache(full, inner)
end

# Write `A`'s stored entries into `full`, whose structure holds them, and zero
# the rest; false when an entry of `A` lies outside that structure.
function _scatter!(full::SparseMatrixCSC, A::SparseMatrixCSC)
    fill!(nonzeros(full), 0.0)
    rows, vals = rowvals(full), nonzeros(full)
    @inbounds for col = 1:size(A, 2)
        k = full.colptr[col]
        last = full.colptr[col+1] - 1
        for j in nzrange(A, col)
            r = rowvals(A)[j]
            while k <= last && rows[k] < r
                k += 1
            end
            (k > last || rows[k] != r) && return false
            vals[k] = nonzeros(A)[j]
        end
    end
    return true
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::StructuredKLU; kwargs...)
    c = cache.cacheval
    if cache.isfresh
        _scatter!(c.full, cache.A) || error(
            "StructuredKLU: W has a stored entry outside the Jacobian prototype's structure",
        )
        c.inner.A = c.full
        cache.isfresh = false
    end
    c.inner.b = cache.b
    sol = SciMLBase.solve!(c.inner)
    copyto!(cache.u, sol.u)
    return SciMLBase.build_linear_solution(
        alg,
        cache.u,
        nothing,
        cache;
        retcode = sol.retcode,
    )
end

"""
	default_linear_solver(layout::StateLayout)

[`StructuredKLU`](@ref) on a layout with algebraic rows (a mass-matrix
system), `KLUFactorization()` otherwise.
"""
default_linear_solver(layout::StateLayout) =
    all(is_differential(layout)) ? KLUFactorization() : StructuredKLU()
