"""
The `ReactionSourceOperator` — the per-cell reaction source extracted from the residual.

Represents the per-cell reaction source's contribution to the residual:

```
dT[i]    += Q_rxn(T_i, ξ_i) / ρcp_eff(T_i, ξ_i)
dξ[j, i] += dξ_rxn[j](T_i, ξ_i)
```

where `Q_rxn = -Σ heats[r] · rates[r]` and `dξ_rxn[j] = Σ ν[r,j] · rates[r]`.
The contribution is purely cell-local (`stencil = (cells = 0:0,)`); the chain
rule through `ρcp_eff(T, ξ)` is captured automatically by the cell-local AD
in the default `contribute_to_jacobian!`.

# Math correspondence with the existing residual

The existing `compute_rhs!` in [`finite_volume.jl`](../../Discretization/finite_volume.jl)
sums every per-cell source (Q_rxn, Q_rad, S_conv, -div_q) and divides
the lot by `ρcp` once. Splitting the operator means each operator divides its
own piece by `ρcp` and the sum still equals `(Σ sources) / ρcp` — exact-equal
to the existing residual at every state. The Jacobian-by-AD over each
operator's `apply!` captures the per-source chain rule independently; their
sum reconstructs the full `∂dT/∂state` row.

The existing residual path (`Pyrolysis.Problem.residual!`) keeps using the
cached `rhs_cache.Q_rxn` loop unchanged. This operator's `apply!` is
exercised through the `Structured` backend's Jacobian assembly only.
"""

using ForwardDiff
using ..Materials: Material, n_components, n_reactions
using ..Properties: effective_heat_capacity
using ..Physics: evaluate_reactions!, reaction_heat_source, species_source_terms!

"""
	ReactionSourceOperator <: Operator

Stateless tag. The per-problem AD workspace lives on the backend cache, not
on the operator. The workspace is plumbed through the `cache` argument
of `contribute_to_jacobian!`.
"""
struct ReactionSourceOperator <: Operator end

# Cell-local: no derived-field reads / writes at the framework level. Properties
# (ρcp_eff) are recomputed inside `apply_cell!` so the chain rule is captured
# by the local AD frontier; we deliberately do not declare `:ρcp` as a
# `:reads` field here — properties are inlined per operator; a dedicated
# `PropertyEvaluationOperator` is deferred.
reads(::ReactionSourceOperator) = ()
writes(::ReactionSourceOperator) = ()
stencil(::ReactionSourceOperator) = (cells = 0:0,)

"""
	ReactionScratch{Tu, NC, NR}

Per-cell scratch for `ReactionSourceOperator`. Holds the AD-typed `rates`,
`heats`, `species_src`, and the `(1 + NC)`-element output buffer. Owned by
the AD workspace; constructed once at problem-build time.
"""
struct ReactionScratch{Tu,NC,NR}
    rates::Vector{Tu}
    heats::Vector{Tu}
    species_src::Vector{Tu}
    out::Vector{Tu}
end

ReactionScratch{Tu}(NC::Int, NR::Int) where {Tu} = ReactionScratch{Tu,NC,NR}(
    zeros(Tu, NR),
    zeros(Tu, NR),
    zeros(Tu, NC),
    zeros(Tu, 1 + NC),
)

"""
	ReactionSourceADWorkspace{Tu, NC, NR, CFG}

AD workspace for `ReactionSourceOperator.contribute_to_jacobian!`. Cell-local,
so the AD chunk size is `1 + NC` (one Jacobian eval per cell).

# Fields
- `x_buf::Vector{Float64}` — `(1 + NC)`-length input buffer.
- `f_buf::Vector{Float64}` — `(1 + NC)`-length output buffer.
- `J_local::Matrix{Float64}` — `(1 + NC) × (1 + NC)` per-cell Jacobian.
- `cfg::CFG` — `ForwardDiff.JacobianConfig` sized for `(1 + NC)` inputs.
- `scratch::ReactionScratch{Tu, NC, NR}` — Dual-typed inner buffers.
"""
mutable struct ReactionSourceADWorkspace{Tu,NC,NR,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    scratch::ReactionScratch{Tu,NC,NR}
end

"""
	build_reaction_source_workspace(material::Material) -> ReactionSourceADWorkspace

Construct the AD workspace sized to `material`'s component count `NC` and
reaction count `NR`. Chunk size = `1 + NC` (cell-local). Called once at
`_build_jacobian_spec` time.
"""
function build_reaction_source_workspace(material::Material{NC,NR}) where {NC,NR}
    dim = 1 + NC
    x_buf = zeros(Float64, dim)
    f_buf = zeros(Float64, dim)
    J_local = zeros(Float64, dim, dim)
    cfg = ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{dim}())
    Tu = eltype(cfg.duals[1])
    scratch = ReactionScratch{Tu}(NC, NR)
    return ReactionSourceADWorkspace{Tu,NC,NR,typeof(cfg)}(
        x_buf,
        f_buf,
        J_local,
        cfg,
        scratch,
    )
end

"""
	apply_cell!(out, op::ReactionSourceOperator, T, ξ,
	            rates, heats, species_src, material) -> nothing

Per-cell reaction-source contribution. Writes into a length-`(1 + NC)` output
buffer:

- `out[1]     = Q_rxn(T, ξ) / ρcp_eff(T, ξ)` — energy contribution (per ρcp).
- `out[1 + j] = dξ_rxn[j](T, ξ)` for `j ∈ 1:NC` — species source.

Type-stable on whatever scalar type `T` and `ξ` carry — designed to be called
under `ForwardDiff.Dual`. Internal scratch (`rates`, `heats`, `species_src`)
is preallocated on a [`ReactionScratch`](@ref) carried by the AD workspace.
"""
function apply_cell!(
    out::AbstractVector,
    ::ReactionSourceOperator,
    T::Real,
    ξ::NTuple{NC,<:Real},
    rates::AbstractVector,
    heats::AbstractVector,
    species_src::AbstractVector,
    material::Material{NC,NR},
) where {NC,NR}

    fill!(rates, zero(eltype(rates)))
    fill!(heats, zero(eltype(heats)))
    evaluate_reactions!(rates, heats, ξ, T, material)
    Q_rxn = reaction_heat_source(rates, heats)
    ρcp = effective_heat_capacity(ξ, T, material)

    out[1] = Q_rxn / ρcp

    species_source_terms!(species_src, rates, material)
    @inbounds for j = 1:NC
        out[1+j] = species_src[j]
    end
    return nothing
end

"""
	apply!(dT, dξ, op::ReactionSourceOperator, state, t, material, scratch)

Operator-level residual contribution. Iterates over cells and accumulates the
per-cell `(Q_rxn / ρcp, dξ_rxn[:])` into `(dT, dξ)`.

Used by the standalone unit test to verify equivalence against the residual
slice. Not on the `Structured` backend's Jacobian hot path — that path goes
through `contribute_to_jacobian!`, which calls `apply_cell!` directly.
"""
function apply!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    op::ReactionSourceOperator,
    state,
    t,
    material::Material{NC,NR},
    scratch::ReactionScratch{Tu,NC,NR},
) where {Tu,NC,NR}

    n = length(state.T)
    @inbounds for i = 1:n
        T = state.T[i]
        ξ = ntuple(j -> state.ξ[j, i], Val(NC))
        apply_cell!(
            scratch.out,
            op,
            T,
            ξ,
            scratch.rates,
            scratch.heats,
            scratch.species_src,
            material,
        )
        dT[i] += scratch.out[1]
        for j = 1:NC
            dξ[j, i] += scratch.out[1+j]
        end
    end
    return nothing
end

"""
	contribute_to_jacobian!(jac, op::ReactionSourceOperator,
	                        u, ctx, t, ws, material, layout, mesh) -> nothing

Default cell-local AD path. For each live cell `i` (`1:mesh.n_active_cells`):

1. Pack `(T_i, ξ_{1,i}, …, ξ_{NC,i})` into `ws.x_buf`.
2. Run `ForwardDiff.jacobian!` over `apply_cell!` to fill `ws.J_local`
   (a `(1 + NC) × (1 + NC)` block).
3. Scatter `J_local` into the diagonal block of `jac.A` at the global
   indices for cell `i`.

The scatter uses `fast_sparse_add!` so off-pattern writes are silent — the
support-complete verifier flags a missing reaction-diagonal entry by
comparing against `DenseAD_ForwardDiff` reference.
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    op::ReactionSourceOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::ReactionSourceADWorkspace{Tu,NC,NR},
    material::Material{NC,NR},
    layout,
    mesh,
) where {Tu,NC,NR}

    n_cells = layout.n_cells
    scratch = ws.scratch
    f! = let op = op, scratch = scratch, material = material
        (out, x_inner) -> begin
            T = x_inner[1]
            ξ = ntuple(j -> x_inner[1+j], Val(NC))
            apply_cell!(
                out,
                op,
                T,
                ξ,
                scratch.rates,
                scratch.heats,
                scratch.species_src,
                material,
            )
        end
    end

    A = jac.A
    smap = jac.smap
    @inbounds for i = 1:mesh.n_active_cells
        # Pack stencil-local input from the global state vector.
        # Layout: T at u[i]; ξ_{j,i} at u[n_cells + (j-1)*n_cells + i].
        ws.x_buf[1] = u[i]
        for j = 1:NC
            ws.x_buf[1+j] = u[n_cells+(j-1)*n_cells+i]
        end

        ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

        # Scatter the (1 + NC) × (1 + NC) block. Row/col mapping mirrors the
        # packing above. Cell-local: the same i indexes both row and column
        # blocks (no inter-cell coupling for this operator).
        row_T = i
        col_T = i
        fast_sparse_add!(A, smap, row_T, col_T, ws.J_local[1, 1])
        for l = 1:NC
            col_ξ = n_cells + (l - 1) * n_cells + i
            fast_sparse_add!(A, smap, row_T, col_ξ, ws.J_local[1, 1+l])
        end
        for j = 1:NC
            row_ξ = n_cells + (j - 1) * n_cells + i
            fast_sparse_add!(A, smap, row_ξ, col_T, ws.J_local[1+j, 1])
            for l = 1:NC
                col_ξ = n_cells + (l - 1) * n_cells + i
                fast_sparse_add!(A, smap, row_ξ, col_ξ, ws.J_local[1+j, 1+l])
            end
        end
    end
    return nothing
end
