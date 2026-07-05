"""
The `ALEAdvectionOperator` — the ALE advective contribution extracted from the residual.

Represents the per-cell ALE advective contribution to `dT` and (gas-only) `dξ`:

```
dT[i]    += w_cell · ∂T/∂z |_i              (T-row energy contribution)
dξ[j, i] += w_cell · ∂ξⱼ/∂z |_i   (gas j)   (ξ-row species contribution)
```

where `∂T/∂z` and `∂ξⱼ/∂z` are TVD-limited (default VANLEER) gradients
computed by [`Pyrolysis.Discretization._limited_gradient`](@ref). Solid /
liquid components have no advection (they move with the mesh).

# Two chain-rule pieces

The ALE advection Jacobian has *two* couplings, both of which are captured:

1. **Direct chain rule** — `∂(w·∇T)/∂(T_stencil, ξ_stencil)` with
   `w_cell` held fixed. Captured by face-local ForwardDiff over a 3-cell
   upwind window `(i, upwind, upwind-upwind)` (5-cell symmetric window
   when both upwind directions are exercised across the mesh, declared
   as `stencil = (cells = -2:2,)`). The limiter chain rule lives entirely
   here; AD picks it up automatically.

2. **Indirect chain rule** — `∂(w·∇T)/∂(T_j, ξ_{k,j}, χ_j)` for every
   upstream cell `j ≤ max(n_left, n_right) - 1`. The dependency flows
   through `w_cell = 0.5 · (w_node[n_left] + w_node[n_right])` and
   `w_node[k] = Σ_{i<k} (θ_pyrolysis,i − θ_A) · Δzᵢ`. This is the
   cumulative coupling from the mesh-velocity operator. It is captured via
   [`EmbeddedScaledPrefixSumCoupling`](@ref) operators registered onto
   `jac.structured` — one per (target row, source col) channel.

# What this operator does

- `apply!` (wraps existing residual functions; no math duplication)
- `contribute_to_jacobian!`'s **direct** chain rule via face-local AD
  scatters into `jac.A`
- `contribute_to_jacobian!`'s **indirect** chain rule registers
  embedded prefix-sum couplings onto `jac.structured`

`Structured()` walks this operator directly. `DirectSolve` materializes
the embedded couplings, while structured strategies apply them in factored
form.

z-row contributions from ALE advection (`dT_dz · ∂Δz/∂zₘ` and
`dT_dz · ∂w_cell/∂zₘ`) are owned by the existing geometry block
(`geometry_ad.jl`) and not duplicated here.
"""

using ForwardDiff
using ..Materials: Material, n_components, n_reactions, is_gas
using ..Discretization:
    Limiter,
    VANLEER,
    _limited_gradient,
    _get_boundary_temperature,
    _upwind_gradient_distance,
    _upwind_upwind_id

# `compute_ale_advection_term!` is exported by `..Discretization`; the species
# variant is also there.
using ..Discretization: compute_ale_advection_term!, compute_ale_species_advection!

# ============================================================================
# Operator type — declared in `operator_declarations.jl` (the struct plus its
# `reads` / `writes` / `stencil` declarations live there so the type is visible
# to the topology sort in one place). This file provides the `apply!` and
# `contribute_to_jacobian!` implementations.
# ============================================================================

# ============================================================================
# AD workspace
# ============================================================================

"""
	ALEAdvectionADWorkspace{NC, CFG_T, CFG_X}

AD workspace for [`ALEAdvectionOperator`](@ref)`.contribute_to_jacobian!`.
Holds two ForwardDiff configs:

- `cfg_T` — for the energy-row direct chain rule. Stencil is 3 cells
  (i, upwind, upwind-upwind) × `(1 + NC)` per cell = `3·(1+NC)` AD inputs.
  Output is scalar `(w·∇T)_i`.
- `cfg_x` — for the species-row direct chain rule. Same `3·(1+NC)` input
  size; output is scalar `(w·∇ξⱼ)_i` for one gas component at a time.

Cell-averaged mesh velocity `w_cell` is held fixed during these AD passes
— it enters as a Float64 closure capture, not as a Dual variable. The
indirect chain rule through `w_cell` is supplied separately by
[`EmbeddedScaledPrefixSumCoupling`](@ref) terms registered on
`jac.structured`.

# Fields
- `x_buf::Vector{Float64}` — `3·(1+NC)`-length packed stencil input.
- `f_buf::Vector{Float64}` — length-1 scalar output.
- `J_local::Matrix{Float64}` — `1 × 3·(1+NC)` per-cell Jacobian row.
- `cfg_T::CFG_T`, `cfg_x::CFG_X` — `ForwardDiff.JacobianConfig`s.

The two configs are identical at this signature; we keep them separate
to allow future variant chunk sizes (e.g., a smaller chunk if the
species closure has fewer active inputs in the `is_gas[j] == false`
case).
"""
mutable struct ALEAdvectionADWorkspace{NC,CFG_T,CFG_X}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg_T::CFG_T
    cfg_x::CFG_X

    # Dense buffers for the indirect-chain-rule prefix-sum.
    # `dθ_dT` and `dθ_dξ` come from the upstream `ALEMeshVelocityOperator`'s
    # workspace; we cache the per-Jacobian-step `n_left/n_right` mapping plus
    # precomputed `Δz` here.
    n_left::Vector{Int}
    n_right::Vector{Int}
    Δz_per_cell::Vector{Float64}

    # `dT_dz` per cell — populated during contribute_to_jacobian! and reused
    # by every embedded coupling that targets the T row.
    dT_dz::Vector{Float64}

    # `dξ_dz[j, i]` — populated for gas components only.
    dξ_dz::Matrix{Float64}

    # `is_gas[j]` cache.
    is_gas::Vector{Bool}
end

"""
	build_ale_advection_workspace(material::Material{NC}, n_cells, n_nodes)
	    -> ALEAdvectionADWorkspace

Construct the AD workspace. Sized once at problem-build time, reused
across Newton steps.
"""
function build_ale_advection_workspace(
    material::Material{NC},
    n_cells::Integer,
    n_nodes::Integer,
) where {NC}
    stencil_dim = 3 * (1 + NC)
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, 1)
    J_local = zeros(Float64, 1, stencil_dim)
    cfg_T =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    cfg_x =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())

    is_gas_vec = [is_gas(material.components[j]) for j = 1:NC]

    return ALEAdvectionADWorkspace{NC,typeof(cfg_T),typeof(cfg_x)}(
        x_buf,
        f_buf,
        J_local,
        cfg_T,
        cfg_x,
        zeros(Int, n_cells),
        zeros(Int, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, NC, n_cells),
        is_gas_vec,
    )
end

# ============================================================================
# apply! — wraps existing residual functions
# ============================================================================

"""
	apply!(dT, dξ, op::ALEAdvectionOperator, state, mesh, w_node, material;
	       limiter = VANLEER) -> nothing

Operator-level ALE advection contribution. Calls
[`compute_ale_advection_term!`](@ref) for the energy contribution and
[`compute_ale_species_advection!`](@ref) for the gas-species contribution
(solid / liquid components are no-ops). The results are *added* into
`dT` and `dξ` (not overwritten) so the operator composes additively
with the other operators.

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
    limiter::Limiter = VANLEER,
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

    # T-row contribution: scratch buffer for the per-cell advection_T values,
    # then add into dT.
    advection_T = similar(dT)
    fill!(advection_T, zero(eltype(dT)))
    compute_ale_advection_term!(advection_T, state, mesh, w_node; limiter = limiter)

    # Species advection: zeros for non-gas rows, and the function only writes
    # gas rows, so we can scratch-allocate and then add.
    advection_ξ = similar(dξ)
    fill!(advection_ξ, zero(eltype(dξ)))
    compute_ale_species_advection!(
        advection_ξ,
        state,
        mesh,
        material,
        w_node;
        limiter = limiter,
    )

    @inbounds for i = 1:n_cells
        dT[i] += advection_T[i]
        for j = 1:NC
            dξ[j, i] += advection_ξ[j, i]
        end
    end

    return nothing
end

# ============================================================================
# Direct-chain-rule helpers — face-local AD over `_limited_gradient`
# ============================================================================

# Look up the upwind / upwind-upwind cell ids and the relevant distances for
# cell `i`. The donor side follows the sign of the apparent transport velocity
# c = -w_cell (mirrors the residual's donor-side selection in
# `compute_ale_advection_term!` / `compute_ale_species_advection!`).
@inline function _ale_advection_geometry(state, mesh, i::Int, w_cell::Float64)
    cell = mesh.cells[i]
    f_left, f_right = cell.face_ids
    face_left = mesh.faces[f_left]
    face_right = mesh.faces[f_right]
    left_id = face_left.cell_ids[1]
    right_id = face_right.cell_ids[2]

    if w_cell > 0
        # c = -w < 0: donor on the right
        upwind_id = right_id
        upup_id = _upwind_upwind_id(mesh, right_id, false)
        d_to_upwind = _upwind_gradient_distance(state, i, upwind_id)
        d_upwind_to_upup =
            upup_id != 0 ?
            abs(state.cell_centers[upwind_id] - state.cell_centers[upup_id]) : d_to_upwind
        upwind_from_left = false
    else
        # c = -w > 0: donor on the left
        upwind_id = left_id
        upup_id = _upwind_upwind_id(mesh, left_id, true)
        d_to_upwind = _upwind_gradient_distance(state, i, upwind_id)
        d_upwind_to_upup =
            upup_id != 0 ?
            abs(state.cell_centers[upwind_id] - state.cell_centers[upup_id]) : d_to_upwind
        upwind_from_left = true
    end
    has_upup = upup_id != 0
    return (
        upwind_id = upwind_id,
        upup_id = upup_id,
        d_to_upwind = Float64(d_to_upwind),
        d_upwind_to_upup = Float64(d_upwind_to_upup),
        has_upup = has_upup,
        upwind_from_left = upwind_from_left,
    )
end

# Pack stencil-local (T, ξ) for cell i into a flat input vector.
# Layout: x[1] = T_i, x[2] = T_upwind, x[3] = T_upup,
#         x[4 .. 3+NC]      = ξ_{:,i}
#         x[4+NC .. 3+2NC]  = ξ_{:,upwind}
#         x[4+2NC .. 3+3NC] = ξ_{:,upup}.
# When the upwind / upup neighbor is the boundary (id == 0), self-fill with
# the cell's own values — the limiter falls back to first-order in those cases
# and the corresponding column derivative is correctly zero.
@inline function _pack_advection_stencil!(
    x::AbstractVector,
    state,
    i::Int,
    upwind_id::Int,
    upup_id::Int,
    NC::Int,
)
    # Resolve safe ids (self-fill for boundary cases).
    uw = upwind_id == 0 ? i : upwind_id
    uu = upup_id == 0 ? uw : upup_id
    @inbounds begin
        x[1] = state.T[i]
        x[2] = state.T[uw]
        x[3] = state.T[uu]
        base_i = 4
        base_u = 4 + NC
        base_uu = 4 + 2 * NC
        for j = 1:NC
            x[base_i-1+j] = state.ξ[j, i]
            x[base_u-1+j] = state.ξ[j, uw]
            x[base_uu-1+j] = state.ξ[j, uu]
        end
    end
    return nothing
end

# Unpack a Tu input vector into (T_i, T_upwind, T_upup, ξ_i, ξ_upwind, ξ_upup)
# tuples. Takes `Val(NC)` — a runtime-Int `ntuple` length makes the tuples
# unknown-length to inference, which heap-allocates them on every AD kernel
# evaluation.
@inline function _unpack_advection_stencil(x::AbstractVector, ::Val{NC}) where {NC}
    @inbounds T_i = x[1]
    @inbounds T_upwind = x[2]
    @inbounds T_upup = x[3]
    ξ_i = ntuple(j -> @inbounds(x[3+j]), Val(NC))
    ξ_upwind = ntuple(j -> @inbounds(x[3+NC+j]), Val(NC))
    ξ_upup = ntuple(j -> @inbounds(x[3+2*NC+j]), Val(NC))
    return T_i, T_upwind, T_upup, ξ_i, ξ_upwind, ξ_upup
end

# Cell-local closure for the energy advection: `(T·, ξ·) ↦ w_cell · dT/dz`.
function _ale_energy_advection_kernel(
    out::AbstractVector,
    x::AbstractVector,
    w_cell::Float64,
    d_to_upwind::Float64,
    d_upwind_to_upup::Float64,
    has_upup::Bool,
    upwind_from_left::Bool,
    limiter::Limiter,
    nc::Val,
)
    T_i, T_upwind, T_upup, _, _, _ = _unpack_advection_stencil(x, nc)
    dT_dz = _limited_gradient(
        T_i,
        T_upwind,
        T_upup,
        d_to_upwind,
        d_upwind_to_upup,
        has_upup,
        limiter,
        upwind_from_left,
    )
    out[1] = w_cell * dT_dz
    return nothing
end

# Cell-local closure for the species advection of one gas component
# `j_species`: `(T·, ξ·) ↦ +w_cell · dξⱼ/dz`. Mirrors the residual's
# `_compute_species_advection_limited_unrolled!`, which writes
# `advection[j, i] = w_cell · dξ_dz` (`ale_fluxes.jl`). The + sign is the
# physical one: expanding the Eulerian advective divergence gives
# ∂ξ/∂t|_χ = RHS − (v_g − w)·∇ξ − ξ∇·v_g, and with v_g ≈ 0 the correction
# ADDED to the RHS is +w·∇ξ — same sign as the energy term (guide eq 11.4).
function _ale_species_advection_kernel(
    out::AbstractVector,
    x::AbstractVector,
    w_cell::Float64,
    d_to_upwind::Float64,
    d_upwind_to_upup::Float64,
    has_upup::Bool,
    upwind_from_left::Bool,
    limiter::Limiter,
    nc::Val,
    j_species::Int,
)
    _, _, _, ξ_i, ξ_upwind, ξ_upup = _unpack_advection_stencil(x, nc)
    ξ_i_j = ξ_i[j_species]
    ξ_uw_j = ξ_upwind[j_species]
    ξ_uu_j = ξ_upup[j_species]
    dξ_dz = _limited_gradient(
        ξ_i_j,
        ξ_uw_j,
        ξ_uu_j,
        d_to_upwind,
        d_upwind_to_upup,
        has_upup,
        limiter,
        upwind_from_left,
    )
    out[1] = w_cell * dξ_dz
    return nothing
end

# ============================================================================
# contribute_to_jacobian! — direct chain rule + indirect coupling registration
# ============================================================================

"""
	contribute_to_jacobian!(jac, op::ALEAdvectionOperator,
	                        state, mesh, material, layout,
	                        w_node, ctx, ws;
	                        limiter = VANLEER) -> nothing

Compute and scatter the ALE advection Jacobian contribution.

# Direct chain rule (writes into `jac.A`)

Face-local ForwardDiff over the 3-cell upwind window `(i, upwind, upup)`:

1. Build the geometry slice for cell `i` (which side is upwind, distances,
   whether `upup` exists).
2. Pack the stencil state into `ws.x_buf`.
3. Run `ForwardDiff.jacobian!` on the energy closure
   `(T·, ξ·) → w_cell · _limited_gradient(...)` to fill `ws.J_local`.
4. Scatter into the T-row of `jac.A` at the columns `(i, upwind, upup) ×
   (T, ξ_:)`. Self-filled boundary slots are NOT scattered (those AD
   columns are identically zero by construction, and adding to the cell's
   own column would double-count).
5. Repeat steps 1–4 with the species closure for each gas component
   `j ∈ 1:NC` whose `is_gas[j]` is true; scatter into the `ξⱼ`-row.

# Indirect chain rule (composes against `ctx.derived_fields[:w_cell]`)

Reads `:w_cell` ← `ALEMeshVelocityOperator`. The upstream operator's
`contribute_to_structured!` must have populated
`ctx.derived_fields[:w_cell]` with bare [`PrefixSumCoupling`](@ref)
primitives — one per source channel:

```
ctx.derived_fields[:w_cell][1]      = PrefixSumCoupling(Δzᵢ · ∂θ_i/∂T_i)
ctx.derived_fields[:w_cell][1 + k]  = PrefixSumCoupling(Δzᵢ · ∂θ_i/∂ξ_{k,i})
```

This method:

1. Computes `dT_dz` and `dξ_dz` (gas-only) for every cell — the row-
   scaling factors.
2. Calls [`compose_into!`](@ref) once per (target row, source col)
   channel — the framework wraps each bare primitive in an
   [`EmbeddedScaledPrefixSumCoupling`](@ref) at the right global
   `(row_offset, col_offset)` with the consumer's `row_scale`:
   - T-row × T-col, T-row × ξₖ-col (k ∈ 1:NC)
   - For each gas j: ξⱼ-row × T-col, ξⱼ-row × ξₖ-col (k ∈ 1:NC)

The consumer never constructs the embedded type or reaches into the
provider's workspace.

# Defaults — boundary handling

Cells where `|w_cell| < 1e-20` skip the AD pass entirely (matches the
residual's early-return). The corresponding row of `jac.A` and any
embedded couplings are left untouched.
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
    limiter::Limiter = VANLEER,
    assembly_profile::Symbol = :full,
) where {NC}

    n_cells = layout.n_cells
    n_nodes = layout.n_nodes
    A = jac.A
    smap = jac.smap
    # Global state length is the row dimension of the assembled sparse block.
    # Avoids importing `state_length` here.
    n_global = size(A, 1)

    A_section = Float64(state.A_section[])

    # ---------------------------------------------------------------
    # Pre-pass: compute per-cell w_cell, dT_dz, dξ_dz, and remember
    # the n_left/n_right node ids + Δz for each cell. These are
    # reused in both the direct AD scatter and the indirect-coupling
    # row-scaling.
    # ---------------------------------------------------------------
    @inbounds for i = 1:n_cells
        cell = mesh.cells[i]
        n_left, n_right = cell.node_ids
        ws.n_left[i] = n_left
        ws.n_right[i] = n_right
        ws.Δz_per_cell[i] = state.volume[i] / A_section
    end

    # ---------------------------------------------------------------
    # Direct chain rule + dT_dz / dξ_dz row-scaling factor population
    # ---------------------------------------------------------------
    @inbounds for i = 1:n_cells
        w_cell = 0.5 * (w_node[ws.n_left[i]] + w_node[ws.n_right[i]])

        if abs(w_cell) < 1e-20
            ws.dT_dz[i] = 0.0
            for j = 1:NC
                ws.dξ_dz[j, i] = 0.0
            end
            continue
        end

        geom = _ale_advection_geometry(state, mesh, i, w_cell)
        _pack_advection_stencil!(ws.x_buf, state, i, geom.upwind_id, geom.upup_id, NC)

        # Resolve column indices for the scatter (self-filled boundary
        # slots map to 0 — skip those writes).
        col_i = i
        col_upwind = geom.upwind_id == 0 ? 0 : geom.upwind_id
        col_upup = geom.upup_id == 0 ? 0 : geom.upup_id

        # ----- T-row direct AD pass -----
        # Capture `Val(NC)` (not the runtime Int) so the kernel's stencil
        # unpack builds statically-sized tuples.
        f_T! = let w_cell = w_cell, geom = geom, limiter = limiter, nc = Val(NC)
            (out, x_inner) -> _ale_energy_advection_kernel(
                out,
                x_inner,
                w_cell,
                geom.d_to_upwind,
                geom.d_upwind_to_upup,
                geom.has_upup,
                geom.upwind_from_left,
                limiter,
                nc,
            )
        end
        ForwardDiff.jacobian!(ws.J_local, f_T!, ws.f_buf, ws.x_buf, ws.cfg_T)

        # Cache dT_dz (used as row-scaling for the indirect chain rule).
        # `ForwardDiff.jacobian!` extracts the function value into its `y`
        # argument (`ws.f_buf`), so the kernel value (w·∇T)_i is already
        # there — no re-run needed.
        ws.dT_dz[i] = ws.f_buf[1] / w_cell

        # Scatter J_local into J[T_row=i, T_col_*] and J[T_row=i, ξ_col_*].
        _scatter_advection_row!(
            A,
            smap,
            i,
            ws.J_local,
            col_i,
            col_upwind,
            col_upup,
            n_cells,
            NC,
        )

        # ----- Species-row direct AD pass (gas components only) -----
        for j_species = 1:NC
            if !ws.is_gas[j_species]
                ws.dξ_dz[j_species, i] = 0.0
                continue
            end

            f_x! =
                let w_cell = w_cell,
                    geom = geom,
                    limiter = limiter,
                    nc = Val(NC),
                    j_species = j_species

                    (out, x_inner) -> _ale_species_advection_kernel(
                        out,
                        x_inner,
                        w_cell,
                        geom.d_to_upwind,
                        geom.d_upwind_to_upup,
                        geom.has_upup,
                        geom.upwind_from_left,
                        limiter,
                        nc,
                        j_species,
                    )
                end
            ForwardDiff.jacobian!(ws.J_local, f_x!, ws.f_buf, ws.x_buf, ws.cfg_x)

            # Cache dξ_dz[j_species, i]. Kernel returns `+w·dξ/dz` and the
            # value is extracted into `ws.f_buf` by `jacobian!`; divide by
            # w_cell to recover `dξ/dz`.
            ws.dξ_dz[j_species, i] = ws.f_buf[1] / w_cell

            ξj_row = n_cells + (j_species - 1) * n_cells + i
            _scatter_advection_row!(
                A,
                smap,
                ξj_row,
                ws.J_local,
                col_i,
                col_upwind,
                col_upup,
                n_cells,
                NC,
            )
        end
    end

    # ---------------------------------------------------------------
    # Indirect chain rule — compose ctx.derived_fields[:w_cell]
    # primitives with this consumer's row-scaling factors via the
    # framework's compose_into!. The provider
    # (ALEMeshVelocityOperator) populates the bare per-channel
    # PrefixSumCoupling(Δzᵢ · ∂θ/∂x) entries upstream.
    #
    # Skipped under `assembly_profile == :ale_reduced`: the reduced
    # sparsity drops both the cumulative-widening physics-block columns
    # and the entire z-column geometry block, so neither the densified
    # `compose_into!` writes nor the WithChi rank-1 χ-column entries
    # have anywhere to land. The dT_dz / dξ_dz arrays computed above
    # are still used by the direct chain rule, so they're left alone.
    # ---------------------------------------------------------------
    if assembly_profile === :full
        w_cell_couplings = ctx.derived_fields[:w_cell]
        length(w_cell_couplings) == 1 + NC || error(
            "ALEAdvectionOperator: ctx.derived_fields[:w_cell] has " *
            "$(length(w_cell_couplings)) entries, expected 1 + NC = $(1 + NC). " *
            "ALEMeshVelocityOperator.contribute_to_structured! must run first.",
        )
        prefix_T = w_cell_couplings[1]::PrefixSumCoupling
        # prefix_ξ[k] = w_cell_couplings[1 + k]; resolved inline below.

        # T-row × T-col coupling: row scaling = dT_dz, col offset = 0.
        compose_into!(
            jac,
            prefix_T,
            ws.n_left,
            ws.n_right,
            ws.dT_dz,
            0,
            0,
            n_cells,
            n_global,
        )
        # T-row × ξₖ-col coupling: row scaling = dT_dz, col offset = k·n_cells.
        for k = 1:NC
            prefix_k = w_cell_couplings[1+k]::PrefixSumCoupling
            ξk_col_off = n_cells + (k - 1) * n_cells
            compose_into!(
                jac,
                prefix_k,
                ws.n_left,
                ws.n_right,
                ws.dT_dz,
                0,
                ξk_col_off,
                n_cells,
                n_global,
            )
        end

        # Gas ξⱼ-row × T-col coupling (and ξⱼ-row × ξₖ-col couplings).
        # Kernel writes `out = +w·dξ_dz` (matching the residual), so the indirect
        # chain rule is `∂(w·dξ_dz)/∂x_j = dξ_dz · ∂w_cell/∂x_j` — row-scaling
        # factor is `+dξ/dz`.
        for j_species = 1:NC
            ws.is_gas[j_species] || continue
            ξj_row_off = n_cells + (j_species - 1) * n_cells
            row_scale = view(ws.dξ_dz, j_species, :)

            # ξⱼ-row × T-col
            compose_into!(
                jac,
                prefix_T,
                ws.n_left,
                ws.n_right,
                row_scale,
                ξj_row_off,
                0,
                n_cells,
                n_global,
            )
            # ξⱼ-row × ξₖ-col (k ∈ 1:NC)
            for k = 1:NC
                prefix_k = w_cell_couplings[1+k]::PrefixSumCoupling
                ξk_col_off = n_cells + (k - 1) * n_cells
                compose_into!(
                    jac,
                    prefix_k,
                    ws.n_left,
                    ws.n_right,
                    row_scale,
                    ξj_row_off,
                    ξk_col_off,
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
        # Consumer-side embedding for T-row × y-col (y ∈ {χ, T, ξ_k}):
        #   ∂(w·∇T)_i/∂yⱼ = dT_dz[i] · 0.5·(u_bare[n_left[i]] + u_bare[n_right[i]]) · v[j]
        # i.e. EmbeddedRank1Coupling(u_T, v, T_row_off=0, y_col_off, n_global)
        # where u_T[i] = dT_dz[i] · 0.5·(u_bare[n_left[i]] + u_bare[n_right[i]]).
        # Same for gas-ξ rows with row_scale = dξ_dz[j_species, :].
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
                theta_A = ctx.derived_fields[:w_node_theta_A]
                push!(channels, ((theta_A[1]::Rank1Coupling).v, 0))   # T channel
                for k = 1:NC
                    ξk_col_off = n_cells + (k - 1) * n_cells
                    push!(channels, ((theta_A[1+k]::Rank1Coupling).v, ξk_col_off))
                end
            end

            # T-row block
            u_T = Vector{Float64}(undef, n_cells)
            @inbounds for i = 1:n_cells
                u_T[i] = ws.dT_dz[i] * 0.5 * (u_bare[ws.n_left[i]] + u_bare[ws.n_right[i]])
            end
            for (v_ch, col_off) in channels
                push!(
                    jac.structured,
                    EmbeddedRank1Coupling(u_T, v_ch, 0, col_off, Int(n_global)),
                )
            end

            # Gas ξⱼ-row blocks
            for j_species = 1:NC
                ws.is_gas[j_species] || continue
                ξj_row_off = n_cells + (j_species - 1) * n_cells
                u_ξ = Vector{Float64}(undef, n_cells)
                @inbounds for i = 1:n_cells
                    u_ξ[i] =
                        ws.dξ_dz[j_species, i] *
                        0.5 *
                        (u_bare[ws.n_left[i]] + u_bare[ws.n_right[i]])
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

# Helper — scatter the 1 × 3·(1+NC) row into J[row, cols] for the cell-
# centered + upwind + upup columns. `J_local` is interpreted as
# `[T_i, T_upwind, T_upup, ξ_{:,i}, ξ_{:,upwind}, ξ_{:,upup}]`.
#
# Boundary self-fill — `_pack_advection_stencil!` self-fills missing upwind
# / upup with the cell's own values (the residual's zero-gradient
# extrapolation). The residual treats the boundary "ghost" as a function
# of T_i (`T_left = T_i`), so the *true* global Jacobian column for T_i
# combines AD columns 1, 2, and 3 (when both upwind and upup are
# boundary), or columns 2 and 3 (when only upup is boundary). The
# corresponding global column for `T_upwind` or `T_upup` is identically
# zero — those cells don't exist as state variables.
#
# Without folding, AD's "as if independent" columns 2 and 3 leak into the
# scatter as if the boundary cells *were* state variables, creating
# spurious entries on `J[T_i, T_i]` for boundary cells.
@inline function _scatter_advection_row!(
    A,
    smap,
    row::Int,
    J_local::AbstractMatrix,
    col_i::Int,
    col_upwind::Int,
    col_upup::Int,
    n_cells::Int,
    NC::Int,
)
    @inbounds begin
        # T columns — fold self-filled boundary AD columns into col_i / col_upwind.
        T_to_i = J_local[1, 1]
        T_to_uw = 0.0
        T_to_uu = 0.0
        if col_upwind != 0
            T_to_uw = J_local[1, 2]
            if col_upup != 0
                T_to_uu = J_local[1, 3]
            else
                # upup is boundary; x[3] was self-filled with T_upwind →
                # AD col 3 belongs to T_upwind globally.
                T_to_uw += J_local[1, 3]
            end
        else
            # upwind is boundary → x[2] = x[3] = T_i → AD cols 2, 3 belong
            # to T_i.
            T_to_i += J_local[1, 2] + J_local[1, 3]
        end

        fast_sparse_add!(A, smap, row, col_i, T_to_i)
        col_upwind != 0 && fast_sparse_add!(A, smap, row, col_upwind, T_to_uw)
        col_upup != 0 && fast_sparse_add!(A, smap, row, col_upup, T_to_uu)

        # ξ columns — same fold logic per component.
        base_i = 4
        base_u = 4 + NC
        base_uu = 4 + 2 * NC
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            ξ_to_i = J_local[1, base_i-1+j]
            ξ_to_uw = 0.0
            ξ_to_uu = 0.0
            if col_upwind != 0
                ξ_to_uw = J_local[1, base_u-1+j]
                if col_upup != 0
                    ξ_to_uu = J_local[1, base_uu-1+j]
                else
                    ξ_to_uw += J_local[1, base_uu-1+j]
                end
            else
                ξ_to_i += J_local[1, base_u-1+j] + J_local[1, base_uu-1+j]
            end

            fast_sparse_add!(A, smap, row, ξ_off + col_i, ξ_to_i)
            col_upwind != 0 && fast_sparse_add!(A, smap, row, ξ_off + col_upwind, ξ_to_uw)
            col_upup != 0 && fast_sparse_add!(A, smap, row, ξ_off + col_upup, ξ_to_uu)
        end
    end
    return nothing
end
