# ============================================================================
# The linear solve of the integrator's W
# ============================================================================

"""
	StructuredKLU(inner = KLUFactorization())

KLU for the integrator's `W = J − M/(γΔt)` on the entries of `W` that have
been nonzero at some step of the solve. The Jacobian's prototype holds every
entry the residual's connectivity allows, most of which are zero at any one
step, and KLU factors a stored zero like any other entry, so factoring `W` on
the prototype's structure carries the fill of entries that never matter. For
a mass matrix that is not a multiple of the identity OrdinaryDiffEq forms `W`
by a sparse broadcast that stores only the entries that compute to nonzero,
so `W`'s own structure changes whenever an entry passes through zero, and a
KLU on that structure would analyse afresh at every such step.

The cache keeps the union of the prototype's entries that have been nonzero
in this solve, with the diagonal, and copies each fresh `W` into a matrix of
that structure, where the inner KLU refactors it numerically. A nonzero
outside the union grows it, and the inner KLU analyses the grown structure
afresh; the union only grows, so the analyses stop once every entry the solve
reaches has appeared. `W`'s stored entries are a subset of the prototype's,
which the integrator hands `init` as an all-ones `W`; an entry of `W` outside
the prototype is an error.
"""
struct StructuredKLU <: LinearSolve.AbstractSparseFactorization
    inner::KLUFactorization
end

StructuredKLU() = StructuredKLU(KLUFactorization())

# The prototype's structure; which of its entries are in the union and where
# each lies in the union's matrix (zero outside it); that matrix; the inner
# KLU cache that factors it; and the number of times the union grew.
mutable struct StructuredKLUCache{C}
    const prototype::SparseMatrixCSC{Float64,Int}
    const in_union::BitVector
    const slot::Vector{Int}
    union::SparseMatrixCSC{Float64,Int}
    const inner::C
    growths::Int
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
    prototype = copy(A)
    in_union = falses(nnz(prototype))
    rows = rowvals(prototype)
    for col = 1:size(prototype, 2), k in nzrange(prototype, col)
        rows[k] == col && (in_union[k] = true)
    end
    slot = zeros(Int, nnz(prototype))
    union = _union_matrix!(slot, prototype, in_union)
    inner = LinearSolve.init(LinearSolve.LinearProblem(union, zeros(length(b))), alg.inner)
    return StructuredKLUCache(prototype, in_union, slot, union, inner, 0)
end

# The matrix on the prototype's entries in the union, zero-valued, with each
# entry's place in it written to `slot`.
function _union_matrix!(slot::Vector{Int}, prototype::SparseMatrixCSC, in_union::BitVector)
    n = size(prototype, 2)
    colptr = Vector{Int}(undef, n + 1)
    rowval = Vector{Int}(undef, count(in_union))
    rows = rowvals(prototype)
    colptr[1] = 1
    m = 0
    @inbounds for col = 1:n
        for k in nzrange(prototype, col)
            if in_union[k]
                m += 1
                rowval[m] = rows[k]
                slot[k] = m
            else
                slot[k] = 0
            end
        end
        colptr[col+1] = m + 1
    end
    return SparseMatrixCSC(size(prototype, 1), n, colptr, rowval, zeros(m))
end

# Walk `A`'s stored entries in the prototype's structure, calling `visit(k, v)`
# with each one's prototype entry `k` and value `v` until it returns false.
# `:outside` when an entry of `A` lies outside the prototype, `:stopped` when a
# visit returned false, `:done` otherwise.
function _walk(visit, prototype::SparseMatrixCSC, A::SparseMatrixCSC)
    rows, Arows, Avals = rowvals(prototype), rowvals(A), nonzeros(A)
    @inbounds for col = 1:size(A, 2)
        k = prototype.colptr[col]
        last = prototype.colptr[col+1] - 1
        for j in nzrange(A, col)
            r = Arows[j]
            while k <= last && rows[k] < r
                k += 1
            end
            (k > last || rows[k] != r) && return :outside
            visit(k, Avals[j]) || return :stopped
        end
    end
    return :done
end

# Copy `A` into the union's matrix; `:stopped` when a nonzero of `A` lies
# outside the union.
function _load!(c::StructuredKLUCache, A::SparseMatrixCSC)
    vals = nonzeros(c.union)
    fill!(vals, 0.0)
    slot = c.slot
    return _walk(c.prototype, A) do k, v
        s = slot[k]
        s > 0 && (@inbounds vals[s] = v; return true)
        return iszero(v)
    end
end

# Add `A`'s nonzeros to the union and rebuild its matrix.
function _grow!(c::StructuredKLUCache, A::SparseMatrixCSC)
    in_union = c.in_union
    _walk(c.prototype, A) do k, v
        iszero(v) || (@inbounds in_union[k] = true)
        return true
    end
    c.union = _union_matrix!(c.slot, c.prototype, in_union)
    c.growths += 1
    return c
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::StructuredKLU; kwargs...)
    c = cache.cacheval
    if cache.isfresh
        status = _load!(c, cache.A)
        if status === :stopped
            _grow!(c, cache.A)
            status = _load!(c, cache.A)
        end
        status === :done || error(
            "StructuredKLU: W has a stored entry outside the Jacobian prototype's structure",
        )
        c.inner.A = c.union
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
	default_linear_solver()

The linear solve of every layout: [`StructuredKLU`](@ref).
"""
default_linear_solver() = StructuredKLU()
