"""
Structured-coupling operators.

A `StructuredCoupling` is a typed linear operator that represents a
coupling block in the Jacobian via a *fast representation* — `O(N)`
apply / inverse — instead of as a dense matrix. Two kinds carry the
weight in this codebase:

- [`PrefixSumCoupling`](@ref) — represents `L · diag(v)` (lower-
  triangular all-ones × cell-local diagonal). Apply is an inclusive
  scan; inverse is a first-difference. This is the right shape for the
  ALE cumulative mesh-velocity coupling
  (`w_node(k) = ∫_0^{z_k} θ(z') dz'`), which is full-rank but has the
  prefix-sum factorization.
- [`SuffixSumCoupling`](@ref) — mirror image of `PrefixSumCoupling`.
  Represents `U · diag(v)` (upper-triangular all-ones × cell-local
  diagonal); apply is a backward inclusive scan, inverse is a
  forward first-difference. The right shape for non-local sources
  whose value at cell `i` depends on every cell *above* it — e.g.
  Beer-Lambert volumetric absorption, where `q_rad[i]` integrates
  over the optical depth from the surface down to cell `i`.
- [`Rank1Coupling`](@ref) — represents `u · vᵀ`. Apply is one inner
  product; inverse uses Sherman-Morrison internally. The right shape
  for genuinely low-rank channels like the column-global θ_A back-
  reaction in WithChi mode.

The [`StructuredJacobian`](@ref) carries a `Vector{StructuredCoupling}`
plus the local sparse block `A`. The full operator is
`A + Σⱼ structuredⱼ`, applied either by `DirectSolve` (densify each
coupling and factor the sum) or by `StructuredSolve` (apply each
coupling via its fast operator without densification).

The `apply!` / `apply_inverse!` / `densify!` methods are wired into
ALE assembly.
"""

"""
	StructuredCoupling

Abstract supertype of every fast-apply linear-operator block in a
[`StructuredJacobian`](@ref). Concrete subtypes override:

- `apply!(y, c, x)` — `y .+= c * x` in-place (additive accumulation).
- `apply_inverse!(y, c, x)` — `y = c⁻¹ x`. Optional; only required
  for `StructuredSolve`. Concretes that have no analytical inverse may
  leave this unimplemented (errors at solve time, not assembly time).
- `densify!(M, c)` — `M .+= dense(c)`. Used by `DirectSolve` to
  collapse the structured term into the sparse direct factorization.
- `Base.size(c)` — `(rows, cols)` of the operator's range / domain.

The coupling captures a derived field's chain-rule contribution; its
*sources* and *targets* (which rows/columns of the global Jacobian it
acts on) are encoded by how it is composed in the StructuredJacobian
assembly driver, not as fields on the coupling itself.
"""
abstract type StructuredCoupling end

# ============================================================================
# PrefixSumCoupling
# ============================================================================

"""
	PrefixSumCoupling{T} <: StructuredCoupling

Represents `L · diag(v)` where `L` is the lower-triangular all-ones
matrix and `v ∈ ℝ^N` is a per-cell scaling. Equivalently:
`(L · diag(v) · x)[k] = Σ_{i ≤ k-1} v[i] · x[i]` (exclusive prefix sum
of `v .* x`).

For ALE this captures the cumulative mesh-velocity coupling:
`v[i] = Δzᵢ · ∂θ_i/∂x_i`, and the apply produces `∂w_node/∂x` for
every node `k` from the cell-local strain-rate derivatives `∂θ_i/∂x_i`.

# Fields
- `v::Vector{T}` — per-cell scaling. Length `N` = number of cells.
- `_buf::Vector{T}` — internal scratch for in-place scans.

# Apply / inverse cost
- `apply!`, `apply_inverse!`: `O(N)` (one pass).
- `densify!`: `O(N²)` (writes the full lower-triangular block).

The apply, inverse, and densify paths are all implemented and validated
against dense references in
[`test_structured_couplings.jl`](../../test/jacobian/test_structured_couplings.jl).
"""
struct PrefixSumCoupling{T<:Real} <: StructuredCoupling
    v::Vector{T}
    _buf::Vector{T}
end

# Aliases `v` (does not copy) — the coupling is a per-eval view onto the
# caller's buffer; callers that need an independent snapshot copy themselves.
PrefixSumCoupling(v::Vector{T}) where {T<:Real} = PrefixSumCoupling{T}(v, similar(v))

Base.size(c::PrefixSumCoupling) = (length(c.v), length(c.v))

"""
	apply!(y, c::PrefixSumCoupling, x) -> y

Compute `y .+= L · diag(c.v) · x` in-place via inclusive scan. `O(N)`.
The accumulation is additive so this composes cleanly with sparse-block
contributions from `A`.
"""
function apply!(y::AbstractVector, c::PrefixSumCoupling, x::AbstractVector)
    N = length(c.v)
    @assert length(y) == N && length(x) == N
    acc = zero(eltype(y))
    @inbounds for k = 1:N
        # (L · diag(v) · x)[k] = Σ_{i ≤ k-1} v[i] · x[i] — exclusive sum.
        # Accumulate the previous cell's contribution before writing y[k].
        y[k] += acc
        acc += c.v[k] * x[k]
    end
    return y
end

"""
	apply_inverse!(y, c::PrefixSumCoupling, x) -> y

Solve `(L · diag(v)) · y = x` via first-difference. `O(N)`.

The matrix `M = L · diag(v)` is **singular**: its first row is zero
and its last column is zero. So `M y = x` requires `x[1] == 0` for
consistency, and `y[N]` is in the right-nullspace (any value satisfies
`M · y = x` once the other entries are pinned). This routine sets
`y[N] = 0` by convention; downstream solvers compose this inverse
with the rest of the system, where the rank-deficient direction is
pinned by the rest of `A`.

Errors if any `c.v[k] ≈ 0` — the first-difference recursion divides
by `v[k]`.
"""
function apply_inverse!(y::AbstractVector, c::PrefixSumCoupling, x::AbstractVector)
    N = length(c.v)
    @assert length(y) == N && length(x) == N
    # `apply!` writes x[k] = Σ_{j=1}^{k-1} v[j] · y[j]. So
    # x[k+1] - x[k] = v[k] · y[k], i.e. y[k] = (x[k+1] - x[k]) / v[k]
    # for k ∈ 1:N-1. y[N] is in the nullspace (last column of M is zero);
    # we set it to zero by convention.
    @inbounds begin
        for k = 1:(N-1)
            vk = c.v[k]
            vk == 0 && error("PrefixSumCoupling: zero scaling at index $k")
            y[k] = (x[k+1] - x[k]) / vk
        end
        y[N] = zero(eltype(y))
    end
    return y
end

"""
	densify!(M, c::PrefixSumCoupling) -> M

Write `M .+= L · diag(c.v)` into the dense `M`. `O(N²)`. Used by
[`DirectSolve`](@ref) when factoring the full `A + dense(structured)`.
"""
function densify!(M::AbstractMatrix, c::PrefixSumCoupling)
    N = length(c.v)
    @assert size(M) == (N, N)
    @inbounds for j = 1:N, i = 1:N
        # (L · diag(v))[i, j] = v[j] if j ≤ i-1 (exclusive lower triangle).
        if j <= i - 1
            M[i, j] += c.v[j]
        end
    end
    return M
end

# ============================================================================
# SuffixSumCoupling
# ============================================================================

"""
	SuffixSumCoupling{T} <: StructuredCoupling

Represents `U · diag(v)` where `U` is the strict upper-triangular all-
ones matrix (`U[i, j] = 1` for `j ≥ i+1`, else 0) and `v ∈ ℝ^N` is a
per-cell scaling. Equivalently:
`(U · diag(v) · x)[k] = Σ_{i ≥ k+1} v[i] · x[i]` (exclusive suffix sum
of `v .* x`).

This is the dual of [`PrefixSumCoupling`](@ref). The right shape for
non-local sources whose value at cell `i` depends on every cell with
index *greater* than `i`. Beer-Lambert volumetric absorption is the
motivating use: `q_rad[i]` depends on the optical depth integrated
from the surface (cell `N`) down to cell `i+1`, so cell `i`'s row has
nonzero entries from cell `i+1` through cell `N`.

# Fields
- `v::Vector{T}` — per-cell scaling. Length `N` = number of cells.
- `_buf::Vector{T}` — internal scratch (reserved for future buffered
  scans; the present implementation is allocation-free without it).

# Apply / inverse cost
- `apply!`, `apply_inverse!`: `O(N)` (one pass).
- `densify!`: `O(N²)` (writes the full upper-triangular block).

The matrix `M = U · diag(v)` is **singular**: its last row is zero
and its first column is zero. So `M y = x` requires `x[N] == 0` for
consistency, and `y[1]` is in the right-nullspace. `apply_inverse!`
sets `y[1] = 0` by convention; downstream solvers compose this
inverse with `A`, where the rank-deficient direction is pinned by
the rest of the system.
"""
struct SuffixSumCoupling{T<:Real} <: StructuredCoupling
    v::Vector{T}
    _buf::Vector{T}
end

# Aliases `v` (does not copy) — same contract as `PrefixSumCoupling`.
SuffixSumCoupling(v::Vector{T}) where {T<:Real} = SuffixSumCoupling{T}(v, similar(v))

Base.size(c::SuffixSumCoupling) = (length(c.v), length(c.v))

"""
	apply!(y, c::SuffixSumCoupling, x) -> y

Compute `y .+= U · diag(c.v) · x` in-place via backward inclusive scan.
`O(N)`. Accumulation is additive.
"""
function apply!(y::AbstractVector, c::SuffixSumCoupling, x::AbstractVector)
    N = length(c.v)
    @assert length(y) == N && length(x) == N
    acc = zero(eltype(y))
    @inbounds for k = N:-1:1
        # (U · diag(v) · x)[k] = Σ_{i ≥ k+1} v[i] · x[i] — exclusive
        # suffix sum. Accumulate the next cell's contribution before
        # writing y[k].
        y[k] += acc
        acc += c.v[k] * x[k]
    end
    return y
end

"""
	apply_inverse!(y, c::SuffixSumCoupling, x) -> y

Solve `(U · diag(v)) · y = x` via forward first-difference. `O(N)`.

Mirror of `PrefixSumCoupling.apply_inverse!`. The matrix is singular —
the last row is zero (consistency requires `x[N] == 0`) and the first
column is zero (`y[1]` is in the nullspace). This routine sets
`y[1] = 0` by convention; downstream solvers compose this inverse with
`A`, where the rank-deficient direction is pinned by the rest of the
system.

Errors if any `c.v[k] ≈ 0` for `k ∈ 2:N` — the first-difference
recursion divides by `v[k]`.
"""
function apply_inverse!(y::AbstractVector, c::SuffixSumCoupling, x::AbstractVector)
    N = length(c.v)
    @assert length(y) == N && length(x) == N
    # `apply!` writes x[k] = Σ_{j=k+1}^{N} v[j] · y[j]. So
    # x[k-1] - x[k] = v[k] · y[k], i.e. y[k] = (x[k-1] - x[k]) / v[k]
    # for k ∈ 2:N. y[1] is in the nullspace (first column of M is zero);
    # we set it to zero by convention.
    @inbounds begin
        y[1] = zero(eltype(y))
        for k = 2:N
            vk = c.v[k]
            vk == 0 && error("SuffixSumCoupling: zero scaling at index $k")
            y[k] = (x[k-1] - x[k]) / vk
        end
    end
    return y
end

"""
	densify!(M, c::SuffixSumCoupling) -> M

Write `M .+= U · diag(c.v)` into the dense `M`. `O(N²)`. Used by
[`DirectSolve`](@ref) when factoring the full `A + dense(structured)`.
"""
function densify!(M::AbstractMatrix, c::SuffixSumCoupling)
    N = length(c.v)
    @assert size(M) == (N, N)
    @inbounds for j = 1:N, i = 1:N
        # (U · diag(v))[i, j] = v[j] if j ≥ i+1 (exclusive upper triangle).
        if j >= i + 1
            M[i, j] += c.v[j]
        end
    end
    return M
end

# ============================================================================
# Rank1Coupling
# ============================================================================

"""
	Rank1Coupling{T} <: StructuredCoupling

Represents `u · vᵀ` — a dense rank-1 outer product. Captures genuinely
low-rank couplings (e.g. column-global θ_A back-reaction in WithChi
ALE mode).

# Fields
- `u::Vector{T}` — left factor (columns of `M = u * v'`).
- `v::Vector{T}` — right factor (rows of `M`).

# Apply / inverse cost
- `apply!`: `O(N)` (one inner product + one scaled add).
- `apply_inverse!`: not analytically invertible on its own (the rank-1
  matrix is singular); composes with the full system via Sherman-
  Morrison in [`StructuredJacobian`](@ref)'s solve path. Leaving
  unimplemented at the coupling level — the SMW combination lives at
  the StructuredSolve strategy.
- `densify!`: `O(N²)`.
"""
struct Rank1Coupling{T<:Real} <: StructuredCoupling
    u::Vector{T}
    v::Vector{T}
end

Base.size(c::Rank1Coupling) = (length(c.u), length(c.v))

function apply!(y::AbstractVector, c::Rank1Coupling, x::AbstractVector)
    @assert length(y) == length(c.u) && length(x) == length(c.v)
    dot_vx = zero(eltype(y))
    @inbounds for i in eachindex(c.v)
        dot_vx += c.v[i] * x[i]
    end
    @inbounds for i in eachindex(c.u)
        y[i] += c.u[i] * dot_vx
    end
    return y
end

function densify!(M::AbstractMatrix, c::Rank1Coupling)
    m, n = size(c)
    @assert size(M) == (m, n)
    @inbounds for j = 1:n, i = 1:m
        M[i, j] += c.u[i] * c.v[j]
    end
    return M
end

# ============================================================================
# Framework composition primitives — wrap bare provider-side primitives in
# globally-embedded couplings and push them onto a `StructuredJacobian`.
# Consumers call these instead of constructing the embedded types directly,
# preserving the provider/consumer separation.
# ============================================================================

# ============================================================================
# EmbeddedNodePrefixSumCoupling — ALE z-rows (dz_k/dt = w_node[k])
# ============================================================================

"""
	EmbeddedNodePrefixSumCoupling{T} <: StructuredCoupling

Embeds a *bare* per-cell prefix-sum coupling at **node rows** — the shape
of the ALE z-row equations. The z-row residual is

```
dz_k/dt = w_node[k] = Σ_{i ≤ k-1} (θ_i − θ_A) · Δzᵢ
```

so `∂(dz_k/dt)/∂x_j = 𝟙[j ≤ k-1] · v[j]` with `v[j] = Δzⱼ · ∂θⱼ/∂xⱼ`
(the `v` of the bare [`PrefixSumCoupling`](@ref) the mesh-velocity
provider builds). This block is the raw node-indexed prefix sum:
`n_nodes` rows × `n_cells` columns, unit row scaling. Node 1 is the fixed bottom node (`w_node[1] ≡ 0`); its
row receives no entries.

# What it represents

```
M[row_offset + k, col_offset + j] += (j ≤ k − 1) · v[j],   k ∈ 1:n_nodes
```

# Fields
- `v::Vector{T}` — per-cell scaling `Δzᵢ · ∂θᵢ/∂xᵢ`. Length `n_cells`.
- `row_offset::Int`, `col_offset::Int` — global offsets (row block starts
  at the layout's z offset).
- `n_nodes::Int` — number of node rows (`n_cells + 1` on a 1-D mesh).
"""
struct EmbeddedNodePrefixSumCoupling{T<:Real} <: StructuredCoupling
    v::Vector{T}
    row_offset::Int
    col_offset::Int
    n_nodes::Int
end

function EmbeddedNodePrefixSumCoupling(
    v::Vector{T},
    row_offset::Integer,
    col_offset::Integer,
    n_nodes::Integer,
    n_global::Integer,
) where {T<:Real}
    _ = n_global
    # Aliases `v` (no copy) — same per-eval rebuild contract as the other
    # embedded couplings.
    return EmbeddedNodePrefixSumCoupling{T}(
        v,
        Int(row_offset),
        Int(col_offset),
        Int(n_nodes),
    )
end

Base.size(c::EmbeddedNodePrefixSumCoupling) = (c.n_nodes, length(c.v))

function densify!(M::AbstractMatrix, c::EmbeddedNodePrefixSumCoupling)
    N = length(c.v)
    row_off = c.row_offset
    col_off = c.col_offset
    rmax = row_off + c.n_nodes
    cmax = col_off + N
    (rmax ≤ size(M, 1) && cmax ≤ size(M, 2)) ||
        throw(BoundsError(M, ((row_off+1):rmax, (col_off+1):cmax)))

    @inbounds for j = 1:N
        vj = c.v[j]
        vj == 0.0 && continue   # Common: ∂θ/∂x is zero outside the reactive zone.
        for k = (j+1):c.n_nodes
            M[row_off+k, col_off+j] += vj
        end
    end
    return M
end

function apply!(y::AbstractVector, c::EmbeddedNodePrefixSumCoupling, x::AbstractVector)
    N = length(c.v)
    row_off = c.row_offset
    col_off = c.col_offset
    (row_off + c.n_nodes ≤ length(y) && col_off + N ≤ length(x)) || throw(
        DimensionMismatch(
            "EmbeddedNodePrefixSumCoupling.apply!: y/x too small for " *
            "row_off=$row_off, col_off=$col_off, n_nodes=$(c.n_nodes), N=$N",
        ),
    )

    # Exclusive prefix scan: y[row_off + k] += Σ_{j ≤ k-1} v[j] · x[col_off + j].
    acc = zero(eltype(y))
    @inbounds for k = 1:c.n_nodes
        y[row_off+k] += acc
        k ≤ N && (acc += c.v[k] * x[col_off+k])
    end
    return y
end

function apply_inverse!(::AbstractVector, ::EmbeddedNodePrefixSumCoupling, ::AbstractVector)
    error(
        "EmbeddedNodePrefixSumCoupling.apply_inverse!: not defined — the " *
        "embedded coupling is rank-deficient on its own. StructuredSolve " *
        "inverts via A⁻¹ + apply! and does not need this method.",
    )
end

"""
	compose_into!(jac, prefix::PrefixSumCoupling, row_offset, col_offset,
	              n_nodes, n_global) -> jac

Wrap a bare [`PrefixSumCoupling`](@ref) in an
[`EmbeddedNodePrefixSumCoupling`](@ref) at node rows `row_offset .+
(1:n_nodes)` and source columns `col_offset .+ (1:length(prefix.v))`, then
push onto `jac.structured`. This is the z-row embedding of the cumulative
mesh-velocity coupling: `∂(dz_k/dt)/∂x_j = 𝟙[j ≤ k-1] · v[j]`.
"""
function compose_into!(
    jac,
    prefix::PrefixSumCoupling,
    row_offset::Integer,
    col_offset::Integer,
    n_nodes::Integer,
    n_global::Integer,
)
    push!(
        jac.structured,
        EmbeddedNodePrefixSumCoupling(
            prefix.v,
            Int(row_offset),
            Int(col_offset),
            Int(n_nodes),
            Int(n_global),
        ),
    )
    return jac
end

# ============================================================================
# EmbeddedScaledSuffixSumCoupling — Beer-Lambert volumetric absorption
# ============================================================================

"""
	EmbeddedScaledSuffixSumCoupling{T} <: StructuredCoupling

Embedded, per-row-scaled suffix sum for non-local sources that read every
cell **above** the current cell. Beer-Lambert volumetric
absorption is the canonical case: `q_rad[i]` integrates the optical depth
from the surface (cell `n`) down to cell `i`, so cell `i`'s residual row
depends on every cell `k ≥ i`.

# What it represents

```
M[row_offset + i, col_offset + j] =
    row_scale[i] · (j ≥ i) · v[j]
```

where `v[j]` carries the per-cell scaling (typically `∂α_j/∂x_j` × the
optical-depth Jacobian factor) and `row_scale[i]` carries the consumer's
per-cell coefficient (typically `−q_rad[i]` for the energy row).

# Fields
- `v::Vector{T}` — per-cell scaling.
- `row_scale::Vector{T}` — per-target-row consumer coefficient.
- `row_offset::Int`, `col_offset::Int` — global offsets.
- `n_cells::Int` — cached `length(v)`.
"""
struct EmbeddedScaledSuffixSumCoupling{T<:Real,RS<:AbstractVector{T}} <: StructuredCoupling
    v::Vector{T}
    row_scale::RS
    row_offset::Int
    col_offset::Int
    n_cells::Int
end

function EmbeddedScaledSuffixSumCoupling(
    v::Vector{T},
    row_scale::AbstractVector{T},
    row_offset::Integer,
    col_offset::Integer,
    n_cells::Integer,
    n_global::Integer,
) where {T<:Real}
    length(v) == n_cells || throw(
        DimensionMismatch(
            "EmbeddedScaledSuffixSumCoupling: v length $(length(v)) ≠ n_cells $(n_cells)",
        ),
    )
    length(row_scale) == n_cells || throw(
        DimensionMismatch(
            "EmbeddedScaledSuffixSumCoupling: row_scale length $(length(row_scale)) ≠ n_cells $(n_cells)",
        ),
    )
    _ = n_global
    # Aliases `v` and `row_scale` (no copies) — per-eval rebuild contract.
    return EmbeddedScaledSuffixSumCoupling{T,typeof(row_scale)}(
        v,
        row_scale,
        Int(row_offset),
        Int(col_offset),
        Int(n_cells),
    )
end

Base.size(c::EmbeddedScaledSuffixSumCoupling) = (c.n_cells, c.n_cells)

function densify!(M::AbstractMatrix, c::EmbeddedScaledSuffixSumCoupling)
    N = c.n_cells
    row_off = c.row_offset
    col_off = c.col_offset
    rmax = row_off + N
    cmax = col_off + N
    (rmax ≤ size(M, 1) && cmax ≤ size(M, 2)) ||
        throw(BoundsError(M, ((row_off+1):rmax, (col_off+1):cmax)))

    @inbounds for j = 1:N
        vj = c.v[j]
        vj == 0.0 && continue
        # Strict upper indicator: nonzero rows are `i < j` (j > i). The
        # diagonal `i == j` is not captured here; consumers add it
        # analytically.
        for i = 1:(j-1)
            M[row_off+i, col_off+j] += c.row_scale[i] * vj
        end
    end
    return M
end

function apply!(y::AbstractVector, c::EmbeddedScaledSuffixSumCoupling, x::AbstractVector)
    N = c.n_cells
    row_off = c.row_offset
    col_off = c.col_offset
    (row_off + N ≤ length(y) && col_off + N ≤ length(x)) || throw(
        DimensionMismatch(
            "EmbeddedScaledSuffixSumCoupling.apply!: y/x " *
            "too small for row_off=$row_off, col_off=$col_off, N=$N",
        ),
    )

    # Strict-suffix sum: S[i] = Σ_{j > i} v[j] · x[col_off + j]
    S = zeros(eltype(y), N)
    @inbounds begin
        acc = zero(eltype(y))
        for i = N:-1:1
            # Write before adding cell i's own contribution — exclusive scan.
            S[i] = acc
            acc += c.v[i] * x[col_off+i]
        end
        for i = 1:N
            y[row_off+i] += c.row_scale[i] * S[i]
        end
    end
    return y
end

function apply_inverse!(
    ::AbstractVector,
    ::EmbeddedScaledSuffixSumCoupling,
    ::AbstractVector,
)
    error(
        "EmbeddedScaledSuffixSumCoupling.apply_inverse!: not defined — " *
        "the embedded coupling is rank-deficient on its own.",
    )
end

"""
	compose_into!(jac, suffix::SuffixSumCoupling, row_scale,
	              row_offset, col_offset, n_cells, n_global) -> jac

Wrap a bare [`SuffixSumCoupling`](@ref) in an
[`EmbeddedScaledSuffixSumCoupling`](@ref) at the consumer's offsets and
per-cell `row_scale`, then push onto `jac.structured`.
"""
function compose_into!(
    jac,
    suffix::SuffixSumCoupling,
    row_scale::AbstractVector{<:Real},
    row_offset::Integer,
    col_offset::Integer,
    n_cells::Integer,
    n_global::Integer,
)
    rs = eltype(row_scale) === Float64 ? row_scale : collect(Float64, row_scale)
    push!(
        jac.structured,
        EmbeddedScaledSuffixSumCoupling(
            suffix.v,
            rs,
            Int(row_offset),
            Int(col_offset),
            Int(n_cells),
            Int(n_global),
        ),
    )
    return jac
end

# ============================================================================
# EmbeddedRank1Coupling — column-global χ back-reaction (WithChi)
# ============================================================================

"""
	EmbeddedRank1Coupling{T} <: StructuredCoupling

Embeds a bare [`Rank1Coupling`](@ref) `u · vᵀ` at consumer-specified
`(row_offset, col_offset)`. Captures genuinely low-rank channels — most
notably the column-global θ_A back-reaction in WithChi mode, where every
target row's coefficient depends on a single scalar functional of the χ
vector.

# What it represents

```
M[row_offset + i, col_offset + j] += u[i] · v[j]
```

# Fields
- `u::Vector{T}` — left factor (target rows). Length `n_rows`.
- `v::Vector{T}` — right factor (source columns). Length `n_cols`.
- `row_offset::Int`, `col_offset::Int` — global offsets.
"""
struct EmbeddedRank1Coupling{T<:Real} <: StructuredCoupling
    u::Vector{T}
    v::Vector{T}
    row_offset::Int
    col_offset::Int
end

function EmbeddedRank1Coupling(
    u::Vector{T},
    v::Vector{T},
    row_offset::Integer,
    col_offset::Integer,
    n_global::Integer,
) where {T<:Real}
    _ = n_global
    # Aliases `u` and `v` (no copies) — per-eval rebuild contract.
    return EmbeddedRank1Coupling{T}(u, v, Int(row_offset), Int(col_offset))
end

Base.size(c::EmbeddedRank1Coupling) = (length(c.u), length(c.v))

function densify!(M::AbstractMatrix, c::EmbeddedRank1Coupling)
    m = length(c.u)
    n = length(c.v)
    row_off = c.row_offset
    col_off = c.col_offset
    (row_off + m ≤ size(M, 1) && col_off + n ≤ size(M, 2)) ||
        throw(BoundsError(M, ((row_off+1):(row_off+m), (col_off+1):(col_off+n))))

    @inbounds for j = 1:n
        vj = c.v[j]
        vj == 0.0 && continue
        for i = 1:m
            ui = c.u[i]
            ui == 0.0 && continue
            M[row_off+i, col_off+j] += ui * vj
        end
    end
    return M
end

function apply!(y::AbstractVector, c::EmbeddedRank1Coupling, x::AbstractVector)
    m = length(c.u)
    n = length(c.v)
    row_off = c.row_offset
    col_off = c.col_offset
    (row_off + m ≤ length(y) && col_off + n ≤ length(x)) ||
        throw(DimensionMismatch("EmbeddedRank1Coupling.apply!"))

    dot_vx = zero(eltype(y))
    @inbounds for j = 1:n
        dot_vx += c.v[j] * x[col_off+j]
    end
    @inbounds for i = 1:m
        y[row_off+i] += c.u[i] * dot_vx
    end
    return y
end

function apply_inverse!(::AbstractVector, ::EmbeddedRank1Coupling, ::AbstractVector)
    error(
        "EmbeddedRank1Coupling.apply_inverse!: rank-1 inverse via " *
        "Sherman-Morrison is the StructuredSolve solver's responsibility.",
    )
end

"""
	compose_into!(jac, rank1::Rank1Coupling, row_offset, col_offset, n_global) -> jac

Wrap a bare [`Rank1Coupling`](@ref) in an [`EmbeddedRank1Coupling`](@ref)
at the consumer's `(row_offset, col_offset)`, then push onto
`jac.structured`. The wrapping accepts the bare rank-1's `u` and `v`
factors verbatim — sizes are read from `length(u) × length(v)`.
"""
function compose_into!(
    jac,
    rank1::Rank1Coupling,
    row_offset::Integer,
    col_offset::Integer,
    n_global::Integer,
)
    push!(
        jac.structured,
        EmbeddedRank1Coupling(
            rank1.u,
            rank1.v,
            Int(row_offset),
            Int(col_offset),
            Int(n_global),
        ),
    )
    return jac
end
