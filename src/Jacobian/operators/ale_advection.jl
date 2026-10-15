"""
The `ALEAdvectionOperator` — the conservative flux-form gas mesh-motion
transport extracted from the residual. **Gas species rows only**:

```
dξ[j, i] += (ξ_face,R·w_R − ξ_face,L·w_L)/Δzᵢ − ξ[j,i]·((w_R − w_L)/Δzᵢ + θ_A)
```

mirroring [`Pyrolysis.Discretization.compute_ale_species_advection!`](@ref)
(first-order donor face values; impermeable boundary faces contribute no
sweep flux; open boundary faces use the zero-gradient ghost). Solid / liquid
components have no advection (they are mesh-bound and move with the mesh).

The **temperature equation has no mesh-motion advection in production**:
condensed species are mesh-bound, so the correct energy correction advects T
by the mesh velocity *relative to the condensed material*,
`w_rel = w_applied − w_physical ≡ 0` when the mesh motion is reaction-driven
(tech-ref §11.2.2). A T-advection term fires only under the
`Workspace.w_override` verification hook — a path the structured Jacobian
does not support (see the note on `w_override`; MMS cases solve with
colored-FD Jacobians). This operator therefore assembles no T-row entries;
the sparsity pattern still over-emits them (harmless structural zeros, same
policy as the solid-component ξ rows).

# Two chain-rule pieces

The gas-transport Jacobian has *two* couplings, both captured:

1. **Direct chain rule** — the term is LINEAR in ξ at fixed `w`, so the
   stencil derivatives are analytic (no AD): `+w_R/Δz` to the right donor's
   column, `−w_L/Δz` to the left donor's column, and
   `−((w_R − w_L)/Δz + θ_A)` to the diagonal, all within the gas
   component's own ξ block (3-cell window, declared
   `stencil = (cells = -1:1,)`).

2. **Indirect chain rule** — `∂f_i/∂w_node[n]` for the two bounding nodes,
   each with its own scale:
   `∂f/∂w_R = (ξ_face,R − ξᵢ)/Δz`, `∂f/∂w_L = −(ξ_face,L − ξᵢ)/Δz`
   (impermeable boundary face: the flux part is absent, leaving the pure
   dilation `∓ξᵢ/Δz`; open boundary face: ghost = ξᵢ makes the scale 0).
   `w_node[k] = Σ_{i<k} (θ_pyrolysis,i − θ_A)·Δzᵢ` is the cumulative
   coupling from the mesh-velocity operator, captured via one
   [`EmbeddedScaledPrefixSumCoupling`](@ref) per (node side, source
   channel) — the `(n, n)` node-pair trick turns the coupling's
   `0.5·(𝟙[j≤n_left−1] + 𝟙[j≤n_right−1])` into a single-node `𝟙[j≤n−1]`.

# What this operator does

- `apply!` (wraps the residual kernel; no math duplication)
- `contribute_to_jacobian!`'s **direct** chain rule writes analytic
  entries into `jac.A`
- `contribute_to_jacobian!`'s **indirect** chain rule registers embedded
  prefix-sum couplings (and WithChi θ_A rank-1 channels) onto
  `jac.structured`

`Structured()` walks this operator directly. `DirectSolve` materializes
the embedded couplings, while structured strategies apply them in factored
form.

z-row contributions (through `∂Δz/∂zₘ` and `∂w/∂zₘ`) are owned by the
existing geometry block (`geometry_ad.jl`) and not duplicated here.
"""

using ..Materials: Material, n_components, n_reactions, is_gas
using ..Discretization: compute_ale_species_advection!

# ============================================================================
# Operator type — declared in `operator_declarations.jl` (the struct plus its
# `reads` / `writes` / `stencil` declarations live there so the type is visible
# to the topology sort in one place). This file provides the `apply!` and
# `contribute_to_jacobian!` implementations.
# ============================================================================

# ============================================================================
# Workspace
# ============================================================================

"""
	ALEAdvectionADWorkspace{NC}

Workspace for [`ALEAdvectionOperator`](@ref)`.contribute_to_jacobian!`.
The flux-form gas transport is linear in ξ at fixed `w`, so no ForwardDiff
machinery is needed — the workspace only caches mesh connectivity and the
per-node indirect-chain-rule scales.

# Fields
- `n_left::Vector{Int}`, `n_right::Vector{Int}` — bounding node ids per cell.
- `Δz_per_cell::Vector{Float64}` — cell heights.
- `scale_L::Matrix{Float64}`, `scale_R::Matrix{Float64}` — `(NC, n_cells)`
  per-gas-row scales `∂f_i/∂w_L`, `∂f_i/∂w_R` consumed (by reference, via
  row views) by the embedded prefix-sum couplings.
- `is_gas::Vector{Bool}` — component phase cache.
"""
mutable struct ALEAdvectionADWorkspace{NC}
    n_left::Vector{Int}
    n_right::Vector{Int}
    Δz_per_cell::Vector{Float64}
    scale_L::Matrix{Float64}
    scale_R::Matrix{Float64}
    is_gas::Vector{Bool}
end

"""
	build_ale_advection_workspace(material::Material{NC}, n_cells, n_nodes)
	    -> ALEAdvectionADWorkspace

Construct the workspace. Sized once at problem-build time, reused across
Newton steps.
"""
function build_ale_advection_workspace(
    material::Material{NC},
    n_cells::Integer,
    n_nodes::Integer,
) where {NC}
    is_gas_vec = [is_gas(material.components[j]) for j = 1:NC]
    return ALEAdvectionADWorkspace{NC}(
        zeros(Int, n_cells),
        zeros(Int, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, NC, n_cells),
        zeros(Float64, NC, n_cells),
        is_gas_vec,
    )
end

# ============================================================================
# apply! — wraps the residual kernel
# ============================================================================

"""
	apply!(dT, dξ, op::ALEAdvectionOperator, state, mesh, w_node, material;
	       θ_A = 0.0, top_gas_open = false, bottom_gas_open = false) -> nothing

Operator-level ALE gas-transport contribution. Calls
[`compute_ale_species_advection!`](@ref) (solid / liquid components are
no-ops). The results are *added* into `dξ` (not overwritten) so the
operator composes additively with the other operators. `dT` is accepted
for interface uniformity and **never touched** — the energy equation has
no mesh-motion advection in production (solids are mesh-bound; see the
module docstring).

`w_node` is supplied by the caller — typically populated by
[`ALEMeshVelocityOperator`](@ref)'s `apply!` running first in the
topo-sorted operator list.
"""
function apply!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    ::ALEAdvectionOperator,
    state,
    mesh,
    w_node::AbstractVector,
    material::Material{NC};
    θ_A::Real = 0.0,
    top_gas_open::Bool = false,
    bottom_gas_open::Bool = false,
) where {NC}

    n_cells = length(state.T)
    n_T = length(dT)
    n_T == n_cells || throw(
        DimensionMismatch(
            "ALEAdvectionOperator.apply!: dT length $(n_T) ≠ n_cells $(n_cells)",
        ),
    )
    size(dξ) == (NC, n_cells) || throw(
        DimensionMismatch(
            "ALEAdvectionOperator.apply!: dξ shape $(size(dξ)) ≠ ($NC, $n_cells)",
        ),
    )

    # Gas transport: zeros for non-gas rows, and the kernel only writes gas
    # rows, so we can scratch-allocate and then add.
    advection_ξ = similar(dξ)
    fill!(advection_ξ, zero(eltype(dξ)))
    compute_ale_species_advection!(
        advection_ξ,
        state,
        mesh,
        material,
        w_node;
        θ_A = θ_A,
        top_gas_open = top_gas_open,
        bottom_gas_open = bottom_gas_open,
    )

    @inbounds for i = 1:n_cells
        for j = 1:NC
            dξ[j, i] += advection_ξ[j, i]
        end
    end

    return nothing
end

# ============================================================================
# contribute_to_jacobian! — analytic direct entries + indirect couplings
# ============================================================================

"""
	contribute_to_jacobian!(jac, op::ALEAdvectionOperator,
	                        state, mesh, material, layout,
	                        w_node, ctx, ws;
	                        θ_A = 0.0, top_gas_open = false,
	                        bottom_gas_open = false,
	                        assembly_profile = :full) -> nothing

Compute and scatter the flux-form gas-transport Jacobian contribution.
No T-row entries are assembled — the energy equation has no mesh-motion
advection in production (see the module docstring).

# Direct chain rule (writes into `jac.A`)

The residual term for gas component j in cell i,

```
f = (ξ_face,R·w_R − ξ_face,L·w_L)/Δz − ξᵢ·((w_R − w_L)/Δz + θ_A)
```

is linear in ξ at fixed `w`, so the entries are analytic:
`+w_R/Δz` to the right-donor column, `−w_L/Δz` to the left-donor column
(folded into the diagonal when the donor is the cell itself, including the
open-boundary zero-gradient ghost), and `−((w_R − w_L)/Δz + θ_A)` to the
diagonal. Impermeable boundary faces contribute no flux entries.

# Indirect chain rule (composes against `ctx.derived_fields[:w_cell]`)

Reads `:w_cell` ← `ALEMeshVelocityOperator` (bare per-channel
[`PrefixSumCoupling`](@ref)s `Δzᵢ·∂θᵢ/∂xᵢ`, x ∈ {T, ξ₁..ξ_NC}). For each
gas row block, two [`compose_into!`](@ref) calls per source channel — one
per bounding node with its own scale (`ws.scale_L`, `ws.scale_R` row
views); passing the same node id for both `n_left`/`n_right` arguments
turns the embedded coupling's node-pair average into the single-node
indicator. WithChi additionally embeds the rank-1 θ_A channels with

```
u[i] = scale_R[j,i]·u_bare[n_right[i]] + scale_L[j,i]·u_bare[n_left[i]] − ξ[j,i]
```

(the last term is the explicit `−ξθ_A` dilation dependence).
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    op::ALEAdvectionOperator,
    state,
    mesh,
    material::Material{NC},
    layout,
    w_node::AbstractVector,
    ctx::ResidualContext,
    ws::ALEAdvectionADWorkspace{NC};
    θ_A::Real = 0.0,
    top_gas_open::Bool = false,
    bottom_gas_open::Bool = false,
    assembly_profile::Symbol = :full,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap
    # Global state length is the row dimension of the assembled sparse block.
    n_global = size(A, 1)

    A_section = Float64(state.A_section[])
    θ_A_val = Float64(θ_A)
    ξ = state.ξ

    # ---------------------------------------------------------------
    # Pre-pass: cache connectivity and Δz.
    # ---------------------------------------------------------------
    @inbounds for i = 1:n_cells
        cell = mesh.cells[i]
        n_left, n_right = cell.node_ids
        ws.n_left[i] = n_left
        ws.n_right[i] = n_right
        ws.Δz_per_cell[i] = state.volume[i] / A_section
    end

    # ---------------------------------------------------------------
    # Direct chain rule (analytic) + per-node indirect scales. Live cells
    # only; the scales of parked ghost slots are zero so the prefix-sum
    # couplings (which span every slot) contribute nothing there.
    # ---------------------------------------------------------------
    n_live = mesh.n_active_cells
    @inbounds for i = (n_live+1):n_cells
        for j = 1:NC
            ws.scale_L[j, i] = 0.0
            ws.scale_R[j, i] = 0.0
        end
    end
    @inbounds for i = 1:n_live
        w_L = Float64(w_node[ws.n_left[i]])
        w_R = Float64(w_node[ws.n_right[i]])
        Δz = ws.Δz_per_cell[i]

        cell = mesh.cells[i]
        f_left, f_right = cell.face_ids
        left_id = mesh.faces[f_left].cell_ids[1]
        right_id = mesh.faces[f_right].cell_ids[2]

        # Donor columns per face (0 = no flux entry: impermeable boundary).
        # Interior: donor is the side the face moves into. Open boundary:
        # zero-gradient ghost = the cell itself.
        donor_R = right_id != 0 ? (w_R > 0 ? right_id : i) : (top_gas_open ? i : 0)
        donor_L = left_id != 0 ? (w_L > 0 ? i : left_id) : (bottom_gas_open ? i : 0)

        for j = 1:NC
            if !ws.is_gas[j]
                ws.scale_L[j, i] = 0.0
                ws.scale_R[j, i] = 0.0
                continue
            end

            ξ_i = Float64(ξ[j, i])
            ξj_row = n_cells + (j - 1) * n_cells + i
            ξj_col_off = n_cells + (j - 1) * n_cells

            # Diagonal: dilation, plus any face whose donor is the cell.
            diag = -((w_R - w_L) / Δz + θ_A_val)
            if donor_R == i
                diag += w_R / Δz
            elseif donor_R != 0
                fast_sparse_add!(A, smap, ξj_row, ξj_col_off + donor_R, w_R / Δz)
            end
            if donor_L == i
                diag -= w_L / Δz
            elseif donor_L != 0
                fast_sparse_add!(A, smap, ξj_row, ξj_col_off + donor_L, -w_L / Δz)
            end
            fast_sparse_add!(A, smap, ξj_row, ξj_col_off + i, diag)

            # Indirect scales ∂f/∂w_node per bounding node.
            ξ_face_R = donor_R == 0 ? ξ_i : Float64(ξ[j, donor_R])  # ξ_i inert when 0
            ξ_face_L = donor_L == 0 ? ξ_i : Float64(ξ[j, donor_L])
            ws.scale_R[j, i] = donor_R == 0 ? -ξ_i / Δz : (ξ_face_R - ξ_i) / Δz
            ws.scale_L[j, i] = donor_L == 0 ? ξ_i / Δz : -(ξ_face_L - ξ_i) / Δz
        end
    end

    # ---------------------------------------------------------------
    # Indirect chain rule — compose ctx.derived_fields[:w_cell]
    # primitives with the per-node scales via the framework's
    # compose_into!. The provider (ALEMeshVelocityOperator) populates
    # the bare per-channel PrefixSumCoupling(Δzᵢ · ∂θ/∂x) entries
    # upstream.
    #
    # Skipped under `assembly_profile == :ale_reduced`: the reduced
    # sparsity drops both the cumulative-widening physics-block columns
    # and the entire z-column geometry block, so neither the densified
    # `compose_into!` writes nor the WithChi rank-1 χ-column entries
    # have anywhere to land. The direct entries above are unaffected.
    # ---------------------------------------------------------------
    if assembly_profile === :full
        w_cell_couplings = ctx.derived_fields[:w_cell]
        length(w_cell_couplings) == 1 + NC || error(
            "ALEAdvectionOperator: ctx.derived_fields[:w_cell] has " *
            "$(length(w_cell_couplings)) entries, expected 1 + NC = $(1 + NC). " *
            "ALEMeshVelocityOperator.contribute_to_structured! must run first.",
        )

        for j_species = 1:NC
            ws.is_gas[j_species] || continue
            ξj_row_off = n_cells + (j_species - 1) * n_cells
            scale_L = view(ws.scale_L, j_species, :)
            scale_R = view(ws.scale_R, j_species, :)

            for (ch, col_off) in enumerate_w_cell_channels(w_cell_couplings, n_cells)
                # Left node: factor 𝟙[j ≤ n_left − 1] with scale_L.
                compose_into!(
                    jac,
                    ch,
                    ws.n_left,
                    ws.n_left,
                    scale_L,
                    ξj_row_off,
                    col_off,
                    n_cells,
                    n_global,
                )
                # Right node: factor 𝟙[j ≤ n_right − 1] with scale_R.
                compose_into!(
                    jac,
                    ch,
                    ws.n_right,
                    ws.n_right,
                    scale_R,
                    ξj_row_off,
                    col_off,
                    n_cells,
                    n_global,
                )
            end
        end

        # ---------------------------------------------------------------
        # WithChi: rank-1 θ_A couplings. Reads `:w_node_chi` (χ channel) and
        # `:w_node_theta_A` (T + ξ channels) ← ALEMeshVelocityOperator. Every
        # channel is a bare Rank1Coupling sharing u[k] = -cumZ_node[k]; the v
        # factor is ∂θ_A/∂χⱼ, ∂θ_A/∂Tⱼ, or ∂θ_A/∂ξ_{k,j} respectively.
        #
        # Consumer-side embedding for ξⱼ-row × y-col (y ∈ {χ, T, ξ_k}):
        #   ∂f_i/∂yₘ = [scale_R·u_bare[n_right] + scale_L·u_bare[n_left]
        #               − ξ[j,i]]·v[m]
        # (the −ξ[j,i] term is the explicit −ξθ_A dilation dependence).
        # ---------------------------------------------------------------
        if haskey(ctx.derived_fields, :w_node_chi)
            w_node_chi = ctx.derived_fields[:w_node_chi][1]::Rank1Coupling
            u_bare = w_node_chi.u
            χ_off = chi_offset(layout)

            # Source channels: (v factor, global column offset) pairs. The χ
            # channel is always present; the θ_A → (T, ξ) back-reaction
            # channels are pushed by the same provider.
            channels = Tuple{Vector{Float64},Int}[(w_node_chi.v, χ_off)]
            if haskey(ctx.derived_fields, :w_node_theta_A)
                theta_A_chs = ctx.derived_fields[:w_node_theta_A]
                push!(channels, ((theta_A_chs[1]::Rank1Coupling).v, 0))   # T channel
                for k = 1:NC
                    ξk_col_off = n_cells + (k - 1) * n_cells
                    push!(channels, ((theta_A_chs[1+k]::Rank1Coupling).v, ξk_col_off))
                end
            end

            # Gas ξⱼ-row blocks
            for j_species = 1:NC
                ws.is_gas[j_species] || continue
                ξj_row_off = n_cells + (j_species - 1) * n_cells
                u_ξ = Vector{Float64}(undef, n_cells)
                @inbounds for i = 1:n_cells
                    u_ξ[i] =
                        ws.scale_R[j_species, i] * u_bare[ws.n_right[i]] +
                        ws.scale_L[j_species, i] * u_bare[ws.n_left[i]] -
                        Float64(ξ[j_species, i])
                end
                for (v_ch, col_off) in channels
                    push!(
                        jac.structured,
                        EmbeddedRank1Coupling(
                            u_ξ,
                            v_ch,
                            ξj_row_off,
                            col_off,
                            Int(n_global),
                        ),
                    )
                end
            end
        end
    end

    return nothing
end

# Iterate the provider's :w_cell channels as (coupling, global column offset)
# pairs: channel 1 is the T column block, channel 1+k the ξₖ block.
function enumerate_w_cell_channels(w_cell_couplings, n_cells::Int)
    return (
        (w_cell_couplings[c]::PrefixSumCoupling, c == 1 ? 0 : n_cells + (c - 2) * n_cells) for
        c = 1:length(w_cell_couplings)
    )
end
