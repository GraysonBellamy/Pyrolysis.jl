# ============================================================================
# The linear solve of the integrator's W
# ============================================================================

const KLU = LinearSolve.KLU

# The places of the diagonal entries of `mass`'s nonzero rows in `J`'s stored values.
function _locate_diagonal!(mass::LumpedMass, J::SparseMatrixCSC)
    rows = rowvals(J)
    for i in eachindex(mass.diagonal)
        iszero(mass.diagonal[i]) && continue
        range = nzrange(J, i)
        push!(mass.shifted, range[searchsortedfirst(view(rows, range), i)])
    end
    return mass
end

function OrdinaryDiffEqDifferentiation.jacobian2W!(
    W::SparseMatrixCSC{Float64,Int},
    M::Diagonal{Float64,LumpedMass},
    dtgamma::Number,
    J::SparseMatrixCSC{Float64,Int},
)
    mass = M.diag
    isempty(mass.shifted) && _locate_diagonal!(mass, J)
    values = nonzeros(W)
    copyto!(values, nonzeros(J))
    shift = inv(dtgamma)
    @inbounds for k in mass.shifted
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
outside the union grows it, and the grown structure is factored afresh, pivots
included, on the solve's one ordering: AMD of the prototype's pattern, taken
once, with KLU's block triangular form of each union. The ordering is symmetric
and the union holds the diagonal, so it suits every union, and nothing is ordered
again as the union grows. The union only grows, so its size alone tells the inner
KLU that the structure changed, and the growths stop once every entry the solve
reaches has appeared. A `W` on any other structure than the prototype's, which
the integrator hands `init` as an all-ones `W`, is an error.
"""
struct StructuredKLU <: LinearSolve.AbstractSparseFactorization
    inner::KLUFactorization
end

StructuredKLU() = StructuredKLU(KLUFactorization(check_pattern = false))

# The prototype's structure; which of its entries are in the union and where
# each lies in the union's matrix (zero outside it); the ordering of every
# factorization, zero-based; that matrix; the inner KLU cache that factors it;
# and the number of times the union grew.
mutable struct StructuredKLUCache{C}
    const prototype::SparseMatrixCSC{Float64,Int}
    const in_union::BitVector
    const slot::Vector{Int}
    const order::Vector{Int}
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
    union = _union_matrix!(slot, prototype, in_union, zeros(nnz(prototype)))
    # Which stored entries count is the union's to decide: a reduction of the stored
    # zeros inside the inner KLU would change the structure that a size alone stands
    # for, and with no pattern check a W would refactor on another W's structure.
    inner = LinearSolve.init(
        LinearSolve.LinearProblem(union, zeros(length(b))),
        alg.inner;
        assumptions = LinearSolve.OperatorAssumptions(
            true;
            nonstructural_zeros = LinearSolve.NonstructuralZeros.None,
        ),
    )
    return StructuredKLUCache(
        prototype,
        in_union,
        slot,
        _fill_reducing_order(prototype),
        union,
        inner,
        0,
    )
end

# AMD of the prototype's pattern, zero-based as KLU takes it.
function _fill_reducing_order(prototype::SparseMatrixCSC{Float64,Int})
    n = size(prototype, 2)
    order = Vector{Int}(undef, n)
    status = KLU.amd_l_order(
        n,
        prototype.colptr .- 1,
        rowvals(prototype) .- 1,
        order,
        C_NULL,
        C_NULL,
    )
    status >= 0 || error("StructuredKLU: AMD failed with status $status")
    return order
end

# The matrix on the prototype's entries in the union, holding `values`, the
# prototype's stored values, with each entry's place in it written to `slot`.
function _union_matrix!(
    slot::Vector{Int},
    prototype::SparseMatrixCSC,
    in_union::BitVector,
    values::Vector{Float64},
)
    n = size(prototype, 2)
    colptr = Vector{Int}(undef, n + 1)
    rowval = Vector{Int}(undef, count(in_union))
    nzval = Vector{Float64}(undef, length(rowval))
    rows = rowvals(prototype)
    colptr[1] = 1
    m = 0
    @inbounds for col = 1:n
        for k in nzrange(prototype, col)
            if in_union[k]
                m += 1
                rowval[m] = rows[k]
                nzval[m] = values[k]
                slot[k] = m
            else
                slot[k] = 0
            end
        end
        colptr[col+1] = m + 1
    end
    return SparseMatrixCSC(size(prototype, 1), n, colptr, rowval, nzval)
end

# Copy `A` into the union's matrix up to the first nonzero of `A` outside the
# union, whose place it returns; zero when there is none.
function _load!(c::StructuredKLUCache, A::SparseMatrixCSC)
    vals = nonzeros(c.union)
    from = nonzeros(A)
    slot = c.slot
    @inbounds for k in eachindex(from)
        s = slot[k]
        if s > 0
            vals[s] = from[k]
        elseif !iszero(from[k])
            return k
        end
    end
    return 0
end

# Add `A`'s nonzeros from place `first` on to the union, every nonzero before it
# being in the union already, and rebuild its matrix with `A`'s values.
function _grow!(c::StructuredKLUCache, A::SparseMatrixCSC, first::Int)
    in_union = c.in_union
    from = nonzeros(A)
    @inbounds for k = first:length(from)
        iszero(from[k]) || (in_union[k] = true)
    end
    c.union = _union_matrix!(c.slot, c.prototype, in_union, from)
    c.growths += 1
    return c
end

# Factor the grown union on the solve's ordering, choosing pivots afresh, and hand
# the factorization to the inner KLU, which refactors it from here on.
function _factor_union!(c::StructuredKLUCache)
    c.inner.A = c.union
    K = KLU.KLUFactorization(c.union)
    KLU.klu_analyze!(K, c.order, c.order)
    KLU.klu_factor!(K; check = false)
    c.inner.cacheval = K
    c.inner.isfresh = false
    return c
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::StructuredKLU; kwargs...)
    c = cache.cacheval
    if cache.isfresh
        A = cache.A
        same_structure(A, c.prototype) ||
            error("StructuredKLU: W is not on the Jacobian prototype's structure")
        outside = _load!(c, A)
        if iszero(outside)
            c.inner.A = c.union
        else
            _grow!(c, A, outside)
            _factor_union!(c)
        end
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
