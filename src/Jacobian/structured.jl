"""
The `StructuredJacobian` — Layer 2 of the rebuilt architecture.

The Jacobian is *not* a sparse matrix; it is a structured operator that
carries:

- a local sparse block `A` (block-tridiagonal physics, plus any sparse
  augmentation from BC implicit terms);
- a list of `StructuredCoupling`s — typed fast-apply linear operators
  (`PrefixSumCoupling` for ALE cumulative mesh velocity, `Rank1Coupling`
  for column-global modes) that are *not* materialized into `A` unless
  the chosen [`LinearSolveStrategy`](@ref) is [`DirectSolve`](@ref);
- a list of `ImplicitTerm`s — IFT-resolved couplings folded into `A`
  as sparse cross-couplings when their sources/targets are sparse.

Backends produce a `StructuredJacobian`; strategies consume it. The two
axes (backend × strategy) are orthogonal — at solve time, the caller
picks the strategy independently of how the Jacobian was assembled.
"""

using SparseArrays: SparseMatrixCSC, sparse, spzeros
using LinearAlgebra: dot, mul!, norm, lu

"""
	ImplicitTerm

IFT-resolved coupling for a derived state that's the output of an
internal Newton solve (today: only `T_s` from
`solve_boundary_temperature!`). When the term's sources/targets are
small and bounded (BC `T_s` touches O(NC) boundary entries), the
contribution is folded into `StructuredJacobian.A` as auxiliary
sparse rows during assembly. When sources or targets are dense, the
implicit term is promoted to a [`StructuredCoupling`](@ref) instead.

# Fields
- `sources::Vector{Int}` — global state indices whose values determine
  the implicit value (e.g. boundary cell index for `T_s`).
- `rows::Vector{Int}` — global Jacobian rows the chain rule flows into.
- `dx_d_sources::Matrix{Float64}` — `∂x_implicit / ∂(state[sources])`
  evaluated at the converged Newton fixed point.
- `drhs_dx::Vector{Float64}` — `∂(rhs[rows]) / ∂x_implicit`.

No operator currently constructs an `ImplicitTerm`; the boundary `T_s`
derivative is taken by the AD-through-Newton path instead.
"""
struct ImplicitTerm
    sources::Vector{Int}
    rows::Vector{Int}
    dx_d_sources::Matrix{Float64}
    drhs_dx::Vector{Float64}
end

"""
	StructuredJacobian{Mode}

Structured representation of the Jacobian. The type parameter `Mode`
records the simulation mode (Eulerian / ALE / WithChi) so dispatch on
`Mode` selects between e.g. an empty `structured` list (Eulerian) and
the prefix-sum + rank-1 list (ALE WithChi).

# Fields
- `A::SparseMatrixCSC{Float64,Int}` — local sparse block. Block-
  tridiagonal physics for non-ALE, plus z-diagonal and χ block for
  ALE. May include the ALE limiter chain rule's enlarged stencil
  (within the 4-cell upwind window).
- `structured::Vector{StructuredCoupling}` — list of typed fast-apply
  linear operators. Empty for Eulerian; `[PrefixSumCoupling]` for ALE
  NoChi; `[PrefixSumCoupling, Rank1Coupling]` for ALE WithChi.
- `implicit::Vector{ImplicitTerm}` — IFT-resolved couplings. Empty in
  the AD-through-Newton path (the default).
- `smap::FastSparsityMap` — O(1) sparse-write helper for `A`.
- `n::Int` — full state dimension (rows == cols of `A` plus structured /
  implicit augmentation).

Apply:
```
(A + Σⱼ structuredⱼ + chain_rule_from_implicit) y = x
```

Assembly is performed by the backend (`_compute_jacobian!` in
[`structured_backend.jl`](structured_backend.jl)). Inversion
is the responsibility of [`LinearSolveStrategy`](@ref) — see
[`strategies.jl`](strategies.jl).
"""
struct StructuredJacobian{Mode}
    A::SparseMatrixCSC{Float64,Int}
    structured::Vector{StructuredCoupling}
    implicit::Vector{ImplicitTerm}
    smap::FastSparsityMap
    n::Int
end

# Mode is encoded as a Symbol (e.g. :eulerian, :ale_nochi, :ale_withchi).
# Constructed by the backend at `_build_jacobian_spec` time.
function StructuredJacobian(
    mode::Symbol,
    A::SparseMatrixCSC{Float64,Int},
    structured = StructuredCoupling[],
    implicit = ImplicitTerm[],
)
    n = size(A, 1)
    @assert size(A, 2) == n "StructuredJacobian: A must be square"
    smap = FastSparsityMap(A)
    return StructuredJacobian{mode}(A, structured, implicit, smap, n)
end

"""
	zero!(jac::StructuredJacobian) -> jac

Reset the Jacobian's accumulator state. `A.nzval .= 0`, structured
couplings are dropped (re-pushed on the next assembly), implicit
terms cleared. Called at the start of each `_compute_jacobian!`.
"""
function zero!(jac::StructuredJacobian)
    fill!(jac.A.nzval, 0.0)
    empty!(jac.structured)
    empty!(jac.implicit)
    return jac
end

# ============================================================================
# Linear-solve strategies (ldiv! methods on StructuredJacobian)
# ============================================================================

"""
	ldiv!(out, strat::LinearSolveStrategy, jac::StructuredJacobian, b) -> out

Apply `(StructuredJacobian)⁻¹ · b` into `out` according to `strat`'s
rule. Each strategy implements one method.

`DirectSolve` densifies the local sparse block plus all structured
couplings and factors that dense matrix. Other strategies provide
specialized approximations or iterative paths.
"""
function ldiv!(
    out::AbstractVector,
    strat::DirectSolve,
    jac::StructuredJacobian,
    b::AbstractVector,
)
    # Densify A + Σⱼ structuredⱼ. The production assembler keeps non-local
    # couplings factored until solve time; DirectSolve materializes them.
    M = Matrix(jac.A)
    for c in jac.structured
        densify!(M, c)
    end
    out .= M \ b
    return out
end

function _apply_structured_jacobian!(
    out::AbstractVector,
    jac::StructuredJacobian,
    x::AbstractVector,
)
    mul!(out, jac.A, x)
    for c in jac.structured
        apply!(out, c, x)
    end
    return out
end

"""
	ldiv!(out, ::StructuredSolve, jac::StructuredJacobian, b;
	      tol = 1e-12, max_iter = 50) -> out

Solve `(A + Σⱼ structuredⱼ) y = b` via Richardson iteration with `A⁻¹`
as the preconditioner. Each `structuredⱼ` is applied through its fast
`apply!` — no densification.

Algorithm:

```
F = lu(A)
y₀ = F \\ b
for iter = 1, 2, …
    r = b − A y − Σⱼ Cⱼ y      (residual of the full operator)
    if ‖r‖ ≤ tol·(‖b‖+1) break
    Δy = F \\ r
    y ← y + Δy
end
```

When the structured corrections are small relative to `A` (the typical
ALE case — `A` carries block-tridiagonal mass / conduction / reaction
terms whose diagonal magnitudes dominate the cumulative-coupling
entries), this converges in 2–5 iterations to machine precision. The
iteration falls back to the direct `A`-solve when `jac.structured` is
empty. If the residual has not met tolerance after `max_iter`
iterations, the solve throws an error instead of returning a stale
iterate.
"""
function ldiv!(
    out::AbstractVector,
    ::StructuredSolve,
    jac::StructuredJacobian,
    b::AbstractVector;
    tol::Real = 1e-12,
    max_iter::Integer = 50,
)
    # No structured corrections — the solve degenerates to A⁻¹ b.
    if isempty(jac.structured)
        out .= jac.A \ b
        return out
    end

    F = lu(jac.A)
    N = length(b)
    y = similar(b)
    r = similar(b)
    Δy = similar(b)
    scratch = similar(b)

    # Initial guess: y = A⁻¹ b.
    copyto!(y, b)
    y .= F \ y

    bnorm = norm(b) + 1.0
    threshold = tol * bnorm
    res_norm = Inf
    converged = false
    for _ = 1:max_iter
        # r ← A·y + Σⱼ Cⱼ·y
        _apply_structured_jacobian!(r, jac, y)
        # r ← b − r
        @inbounds for i = 1:N
            r[i] = b[i] - r[i]
        end

        res_norm = norm(r)
        if res_norm ≤ threshold
            converged = true
            break
        end

        # Δy ← A⁻¹ · r ; y ← y + Δy
        Δy .= F \ r
        @inbounds for i = 1:N
            y[i] += Δy[i]
        end
    end

    if !converged
        _apply_structured_jacobian!(r, jac, y)
        @inbounds for i = 1:N
            r[i] = b[i] - r[i]
        end
        res_norm = norm(r)
        converged = res_norm ≤ threshold
    end

    converged || error(
        "StructuredSolve failed to converge in $max_iter iterations " *
        "(residual norm $res_norm > tolerance $threshold)",
    )

    copyto!(out, y)
    return out
end

function ldiv!(
    out::AbstractVector,
    ::ApproxSolve,
    jac::StructuredJacobian,
    b::AbstractVector,
)
    # Drop structured couplings entirely; solve A only. A Schur-style
    # approximation: A alone is a complete (if approximate) operator.
    out .= jac.A \ b
    return out
end

function ldiv!(
    out::AbstractVector,
    ::MatrixFreeSolve,
    jac::StructuredJacobian,
    b::AbstractVector;
    tol::Real = 1e-10,
    max_iter::Integer = length(b),
)
    N = length(b)
    N == jac.n || throw(
        DimensionMismatch(
            "MatrixFreeSolve: rhs length $N does not match Jacobian size $(jac.n)",
        ),
    )
    length(out) == N || throw(
        DimensionMismatch(
            "MatrixFreeSolve: output length $(length(out)) does not match rhs length $N",
        ),
    )

    beta = norm(b)
    fill!(out, 0.0)
    beta == 0.0 && return out

    m = max(1, min(Int(max_iter), N))
    V = zeros(eltype(b), N, m + 1)
    H = zeros(eltype(b), m + 1, m)
    g = zeros(eltype(b), m + 1)
    w = similar(b)

    @inbounds for i = 1:N
        V[i, 1] = b[i] / beta
    end
    g[1] = beta

    best_res = Inf
    best_k = 0
    best_y = zeros(eltype(b), 0)
    threshold = tol * (beta + 1.0)

    for k = 1:m
        _apply_structured_jacobian!(w, jac, view(V, :, k))

        for j = 1:k
            h = dot(view(V, :, j), w)
            H[j, k] = h
            @inbounds for i = 1:N
                w[i] -= h * V[i, j]
            end
        end

        h_next = norm(w)
        H[k+1, k] = h_next
        if h_next > eps(Float64) * (norm(view(H, 1:k, k)) + 1.0)
            @inbounds for i = 1:N
                V[i, k+1] = w[i] / h_next
            end
        end

        Hk = view(H, 1:(k+1), 1:k)
        gk = view(g, 1:(k+1))
        yk = Hk \ gk
        res = norm(gk - Hk * yk)
        if res < best_res
            best_res = res
            best_k = k
            best_y = collect(yk)
        end
        if res ≤ threshold
            mul!(out, view(V, :, 1:k), yk)
            return out
        end

        h_next == 0.0 && break
    end

    best_k > 0 && mul!(out, view(V, :, 1:best_k), best_y)
    error(
        "MatrixFreeSolve failed to converge in $m iterations " *
        "(residual norm $best_res > tolerance $threshold)",
    )
end
