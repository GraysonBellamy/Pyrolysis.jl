"""
The `StructuredJacobian`.

The Jacobian is assembled as a structured operator that carries:

- a local sparse block `A` (block-tridiagonal physics plus the ALE z and χ
  blocks);
- a list of `StructuredCoupling`s — non-local blocks kept in factored form
  (`EmbeddedNodePrefixSumCoupling` for the ALE cumulative mesh velocity,
  `EmbeddedRank1Coupling` for column-global modes, the Beer–Lambert suffix
  sum) and written into the sparse Jacobian only under the exact `:full`
  profile (see [`LinearSolveStrategy`](@ref)).
"""

using SparseArrays: SparseMatrixCSC

"""
	StructuredJacobian

Structured representation of the Jacobian.

# Fields
- `A::SparseMatrixCSC{Float64,Int}` — local sparse block. Block-
  tridiagonal physics for non-ALE, plus z-diagonal and χ block for
  ALE.
- `structured::Vector{StructuredCoupling}` — non-local blocks. Empty for
  Eulerian without Beer–Lambert; the node prefix sums for ALE, plus the
  rank-1 channels under WithChi.
- `smap::FastSparsityMap` — O(1) sparse-write helper for `A`.
- `n::Int` — full state dimension (rows == cols of `A`).

The operator is `A + Σⱼ structuredⱼ`. Assembly is performed by the backend
(`_compute_jacobian!` in [`structured_backend.jl`](structured_backend.jl)).
"""
struct StructuredJacobian
    A::SparseMatrixCSC{Float64,Int}
    structured::Vector{StructuredCoupling}
    smap::FastSparsityMap
    n::Int
end

function StructuredJacobian(A::SparseMatrixCSC{Float64,Int})
    n = size(A, 1)
    @assert size(A, 2) == n "StructuredJacobian: A must be square"
    return StructuredJacobian(A, StructuredCoupling[], FastSparsityMap(A), n)
end

"""
	zero!(jac::StructuredJacobian) -> jac

Reset the Jacobian's accumulator state: `A.nzval .= 0` and the structured
couplings are dropped (re-pushed on the next assembly). Called at the start
of each `_compute_jacobian!`.
"""
function zero!(jac::StructuredJacobian)
    fill!(jac.A.nzval, 0.0)
    empty!(jac.structured)
    return jac
end
