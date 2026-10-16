"""
Structured-coupling blocks.

A `StructuredCoupling` is a non-local Jacobian block kept in a factored form
during assembly and written into the sparse Jacobian (`densify!`) only by
the exact `:full` profile. Three shapes occur:

- [`EmbeddedNodePrefixSumCoupling`](@ref) — the ALE cumulative mesh-velocity
  coupling of every node row to every cell below it
  (`w_node(k) = ∫_0^{z_k} θ(z') dz'`), built from the per-cell factors of a
  [`PrefixSumCoupling`](@ref).
- [`EmbeddedScaledSuffixSumCoupling`](@ref) — sources whose value at cell `i`
  depends on every cell *above* it (Beer–Lambert volumetric absorption).
- [`EmbeddedRank1Coupling`](@ref) — genuinely low-rank channels such as the
  column-global θ_A back-reaction in WithChi mode.

The [`StructuredJacobian`](@ref) carries a `Vector{StructuredCoupling}` next to
its local sparse block `A`; the full operator is `A + Σⱼ structuredⱼ`.
"""

"""
	StructuredCoupling

Abstract supertype of every non-local block in a [`StructuredJacobian`](@ref).
Concrete subtypes implement `densify!(M, c)` — `M .+= dense(c)` — which the
`:full` assembly profile uses to write the block into the Jacobian. Where the
block acts in the global Jacobian is encoded in the coupling's offsets.
"""
abstract type StructuredCoupling end

# ============================================================================
# PrefixSumCoupling — per-cell factors of the ALE mesh-velocity coupling
# ============================================================================

"""
	PrefixSumCoupling{T}

Per-cell factors `v[i] = Δzᵢ · ∂θᵢ/∂xᵢ` of the ALE cumulative mesh-velocity
coupling `(L · diag(v) · x)[k] = Σ_{i ≤ k-1} v[i] · x[i]`. The mesh-velocity
operator fills one per source block; [`compose_into!`](@ref) embeds it at the
node rows as an [`EmbeddedNodePrefixSumCoupling`](@ref).
"""
struct PrefixSumCoupling{T<:Real}
    v::Vector{T}
end

"""
	compose_into!(jac, prefix::PrefixSumCoupling, row_offset, col_offset,
	              n_nodes, n_global) -> jac

Wrap a [`PrefixSumCoupling`](@ref) in an
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
(the `v` of the [`PrefixSumCoupling`](@ref) the mesh-velocity operator
fills). This block is the raw node-indexed prefix sum: `n_nodes` rows ×
`n_cells` columns, unit row scaling. Node 1 is the fixed bottom node
(`w_node[1] ≡ 0`); its row receives no entries.

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

# ============================================================================
# EmbeddedRank1Coupling — column-global χ back-reaction (WithChi)
# ============================================================================

"""
	EmbeddedRank1Coupling{T} <: StructuredCoupling

Embeds a rank-1 block `u · vᵀ` at consumer-specified
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
