# ============================================================================
# The linear solve of the integrator's W
# ============================================================================

"""
	LumpedMass <: AbstractVector{Float64}

The diagonal of the integrator's mass matrix, one on the rows with a time
derivative and zero on the algebraic rows, with the place of the former's
diagonal entries in the Jacobian's stored values. The matrix is a `Diagonal` of
this vector, which forms `W = J − M/(γΔt)` in place on the Jacobian's own
structure: the sparse broadcast a general mass matrix takes would allocate a
matrix and search it on every step.
"""
struct LumpedMass <: AbstractVector{Float64}
    diagonal::Vector{Float64}
    shifted::Vector{Int}
end

Base.size(m::LumpedMass) = size(m.diagonal)
Base.IndexStyle(::Type{LumpedMass}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(m::LumpedMass, i::Int) = m.diagonal[i]

"""
	lumped_mass(layout, J) -> Diagonal

The integrator's mass matrix of `layout`, for the Jacobian prototype `J`.
"""
function lumped_mass(layout::StateLayout, J::SparseMatrixCSC{Float64,Int})
    diagonal = diag(mass_matrix(layout))
    rows = rowvals(J)
    shifted = [
        nzrange(J, i)[searchsortedfirst(view(rows, nzrange(J, i)), i)] for
        i in eachindex(diagonal) if !iszero(diagonal[i])
    ]
    return Diagonal(LumpedMass(diagonal, shifted))
end

function OrdinaryDiffEqDifferentiation.jacobian2W!(
    W::SparseMatrixCSC{Float64,Int},
    M::Diagonal{Float64,LumpedMass},
    dtgamma::Number,
    J::SparseMatrixCSC{Float64,Int},
)
    values = nonzeros(W)
    copyto!(values, nonzeros(J))
    shift = inv(dtgamma)
    @inbounds for k in M.diag.shifted
        values[k] -= shift
    end
    return nothing
end

"""
	StructuredKLU(inner = KLUFactorization(check_pattern = false))

KLU for the integrator's `W = J − M/(γΔt)` on the entries of `W` that have
been nonzero at some step of the solve. `W` has the structure of the Jacobian's
prototype, which holds every entry the residual's connectivity allows, most of
which are zero at any one step, and KLU factors a stored zero like any other
entry, so factoring `W` on that structure carries the fill of entries that
never matter.

The cache keeps the union of the prototype's entries that have been nonzero
in this solve, with the diagonal, and copies each fresh `W` into a matrix of
that structure, where the inner KLU refactors it numerically. A nonzero
outside the union grows it, and the inner KLU analyses the grown structure
afresh; the union only grows, so its size alone tells the inner KLU that the
structure changed, and the analyses stop once every entry the solve reaches has
appeared. A `W` on any other structure than the prototype's, which the
integrator hands `init` as an all-ones `W`, is an error.
"""
struct StructuredKLU <: LinearSolve.AbstractSparseFactorization
    inner::KLUFactorization
end

StructuredKLU() = StructuredKLU(KLUFactorization(check_pattern = false))

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

# Copy `A` into the union's matrix; false when a nonzero of `A` lies outside the
# union.
function _load!(c::StructuredKLUCache, A::SparseMatrixCSC)
    vals = nonzeros(c.union)
    from = nonzeros(A)
    slot = c.slot
    @inbounds for k in eachindex(from)
        s = slot[k]
        if s > 0
            vals[s] = from[k]
        elseif !iszero(from[k])
            return false
        end
    end
    return true
end

# Add `A`'s nonzeros to the union and rebuild its matrix.
function _grow!(c::StructuredKLUCache, A::SparseMatrixCSC)
    in_union = c.in_union
    from = nonzeros(A)
    @inbounds for k in eachindex(from)
        iszero(from[k]) || (in_union[k] = true)
    end
    c.union = _union_matrix!(c.slot, c.prototype, in_union)
    c.growths += 1
    return c
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::StructuredKLU; kwargs...)
    c = cache.cacheval
    if cache.isfresh
        A = cache.A
        same_structure(A, c.prototype) ||
            error("StructuredKLU: W is not on the Jacobian prototype's structure")
        if !_load!(c, A)
            _grow!(c, A)
            _load!(c, A)
        end
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

The default linear solve: [`StructuredKLU`](@ref).
"""
default_linear_solver() = StructuredKLU()
