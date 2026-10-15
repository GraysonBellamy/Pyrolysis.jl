"""
The `ConductionOperator` — the conductive-divergence contribution extracted from the residual.

Represents the per-cell conductive-divergence contribution to `dT`:

```
dT[i] += -[q_cond[f_right] - q_cond[f_left]] · A / V_i / ρcp_eff(T_i, ξ_i)
```

where each interior conductive flux is
`q_cond[f] = -k_face(ξ_L, T_L, ξ_R, T_R) · (T_R − T_L) / d_LR` and `k_face`
is the distance-weighted harmonic mean of the cell-centered effective
conductivities. Boundary faces contribute zero here (boundary fluxes are
the [`BCThermalOperator`](@ref)'s responsibility).

# Stencil

3-cell tridiagonal: `(cells = -1:1,)`. Cell `i`'s contribution depends on
`(T, ξ)` at cells `i-1, i, i+1` via the face-centered conductivities. For
boundary cells (`i = 1` or `i = n_cells`), the missing-neighbor face
contributes zero — the closure ignores the corresponding inputs and the AD
column is identically zero, which the operator's scatter respects (no write
to the missing entry).

# Math correspondence with the existing residual

`compute_rhs!` populates `q_cond` via `compute_conductive_fluxes!`, copies
into `q_total`, then BC application **overwrites** `q_total` at boundary
faces. The conduction-only piece is the divergence of the *un-overwritten*
`q_cond` (which is zero at boundary faces by construction).
`BCThermalOperator` adds the boundary-face flux's divergence to
`dT[boundary_cell]`.
"""

using ForwardDiff
using ..Materials: Material
using ..Properties: effective_conductivity, effective_heat_capacity, harmonic_mean_weighted

"""
	ConductionOperator <: Operator

Stateless tag. Per-problem AD workspace lives on the backend cache.
"""
struct ConductionOperator <: Operator end

reads(::ConductionOperator) = ()
writes(::ConductionOperator) = ()
stencil(::ConductionOperator) = (cells = -1:1,)

"""
	ConductionADWorkspace{NC, CFG}

Workspace for `ConductionOperator.contribute_to_jacobian!`. Carries both
the AD-default plumbing (`x_buf`, `f_buf`, `J_local`, `cfg`) and the
analytical override's per-cell scratch arrays (`dQ_dξ_L_buf` / `_C_buf`
/ `_R_buf`) so the per-cell ntuple-closure allocations are avoided.

# Fields
- `x_buf::Vector{Float64}` — `(3 · (1 + NC))`-length input buffer for
  the AD-default's stencil pack.
- `f_buf::Vector{Float64}` — length-1 output (`dT_cond[i]`).
- `J_local::Matrix{Float64}` — `1 × (3 · (1 + NC))` per-cell Jacobian row.
- `cfg::CFG` — `ForwardDiff.JacobianConfig` sized for the stencil dim.
- `dQ_dξ_L_buf, dQ_dξ_C_buf, dQ_dξ_R_buf::Vector{Float64}` — length-`NC`
  scratch for the analytical override's per-cell ξ-derivative
  accumulators.
"""
mutable struct ConductionADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    dQ_dξ_L_buf::Vector{Float64}
    dQ_dξ_C_buf::Vector{Float64}
    dQ_dξ_R_buf::Vector{Float64}
end

"""
	build_conduction_workspace(material::Material{NC}) -> ConductionADWorkspace

Construct the workspace. Chunk size = `3 · (1 + NC)`.
"""
function build_conduction_workspace(::Material{NC}) where {NC}
    stencil_dim = 3 * (1 + NC)
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, 1)
    J_local = zeros(Float64, 1, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return ConductionADWorkspace{NC,typeof(cfg)}(
        x_buf,
        f_buf,
        J_local,
        cfg,
        zeros(Float64, NC),
        zeros(Float64, NC),
        zeros(Float64, NC),
    )
end

# Pack stencil-local (T, ξ) for cell i into a flat input vector.
# Layout: x[1] = T_L, x[2] = T_C, x[3] = T_R,
#         x[4 .. 3+NC]      = ξ_{:,L},
#         x[4+NC .. 3+2*NC] = ξ_{:,C},
#         x[4+2*NC .. 3+3*NC] = ξ_{:,R}.
@inline function _pack_stencil!(
    x::AbstractVector,
    u::AbstractVector,
    n_cells::Int,
    NC::Int,
    i::Int,
)
    # Boundary cells: missing neighbor is filled with the cell's own values
    # so the closure (which conditionally drops the missing flux) sees stable
    # inputs. The AD will produce a zero column for the unused slot.
    iL = i == 1 ? i : i - 1
    iR = i == n_cells ? i : i + 1
    @inbounds begin
        x[1] = u[iL]
        x[2] = u[i]
        x[3] = u[iR]
        base_L = 4
        base_C = 4 + NC
        base_R = 4 + 2 * NC
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            x[base_L-1+j] = u[ξ_off+iL]
            x[base_C-1+j] = u[ξ_off+i]
            x[base_R-1+j] = u[ξ_off+iR]
        end
    end
    return nothing
end

# Unpack a Tu input vector into (T_L, T_C, T_R, ξ_L, ξ_C, ξ_R) tuples.
@inline function _unpack_stencil(x::AbstractVector, NC::Int)
    @inbounds T_L = x[1]
    @inbounds T_C = x[2]
    @inbounds T_R = x[3]
    ξ_L = ntuple(j -> @inbounds(x[3+j]), NC)
    ξ_C = ntuple(j -> @inbounds(x[3+NC+j]), NC)
    ξ_R = ntuple(j -> @inbounds(x[3+2*NC+j]), NC)
    return T_L, T_C, T_R, ξ_L, ξ_C, ξ_R
end

"""
	apply_cell!(out, op::ConductionOperator,
	            T_L, T_C, T_R, ξ_L, ξ_C, ξ_R,
	            geom, has_left, has_right, material) -> nothing

Per-cell conductive divergence contribution. Writes `out[1] = dT_cond[i]`.

`geom` is a tuple `(d_L_face, d_R_face, d_LC, d_CR, A_left, A_right, V_C)`
holding the face-distance and area constants. Boundary cases are gated by
`has_left` / `has_right` Booleans — when `false`, the corresponding flux
contribution is dropped.

Type-stable on whatever scalar `T` and `ξ` carry.
"""
function apply_cell!(
    out::AbstractVector,
    ::ConductionOperator,
    T_L::Real,
    T_C::Real,
    T_R::Real,
    ξ_L::NTuple{NC,<:Real},
    ξ_C::NTuple{NC,<:Real},
    ξ_R::NTuple{NC,<:Real},
    d_L_face::Real,
    d_R_face::Real,
    d_LC::Real,
    d_CR::Real,
    A_left::Real,
    A_right::Real,
    V_C::Real,
    has_left::Bool,
    has_right::Bool,
    material::Material{NC},
) where {NC}

    # Cell-local conductivity at C (always needed; both faces, if present,
    # share it).
    k_C = effective_conductivity(ξ_C, T_C, material)

    # Face flux at the left face (between cells L and C).
    q_left = if has_left
        k_L = effective_conductivity(ξ_L, T_L, material)
        k_face_L = harmonic_mean_weighted(k_L, k_C, d_L_face, d_LC - d_L_face)
        -k_face_L * (T_C - T_L) / d_LC
    else
        zero(T_C)
    end

    # Face flux at the right face (between cells C and R).
    q_right = if has_right
        k_R = effective_conductivity(ξ_R, T_R, material)
        k_face_R = harmonic_mean_weighted(k_C, k_R, d_R_face, d_CR - d_R_face)
        -k_face_R * (T_R - T_C) / d_CR
    else
        zero(T_C)
    end

    # Divergence × (-1) / ρcp_eff at C. Sign: +q_right - q_left, then negated
    # (so an outflow on the right *cools* the cell), divided by V_C / A and
    # scaled to dT/dt. Matches divergence_1d!: div[i] = (q_R - q_L) * A / V_i,
    # and the residual writes -div[i] / ρcp.
    # The face area carries the geometry; in 1D A_left ≈ A_right but for ALE
    # with χ-shrinkage they may differ slightly.
    div_q = (q_right * A_right - q_left * A_left) / V_C
    ρcp = effective_heat_capacity(ξ_C, T_C, material)
    out[1] = -div_q / ρcp
    return nothing
end

# Geometry slice for cell i — pulled once per cell at assembly time, held as
# Float64 constants in the AD closure.
@inline function _conduction_geometry(state, mesh, i::Int)
    cell = mesh.cells[i]
    f_left, f_right = cell.face_ids
    face_L = mesh.faces[f_left]
    face_R = mesh.faces[f_right]
    has_left = !(face_L.cell_ids[1] == 0 || face_L.cell_ids[2] == 0)
    has_right = !(face_R.cell_ids[1] == 0 || face_R.cell_ids[2] == 0)

    A_left = state.A_face[f_left]
    A_right = state.A_face[f_right]
    V_C = state.volume[i]

    # Distances. For boundary faces the values are not used (gated by
    # has_left/has_right) but must be finite to avoid NaN under AD.
    if has_left
        iL = face_L.cell_ids[1] == i ? face_L.cell_ids[2] : face_L.cell_ids[1]
        d_L_face = abs(state.face_positions[f_left] - state.cell_centers[iL])
        d_LC = abs(state.cell_centers[i] - state.cell_centers[iL])
    else
        d_L_face = 1.0
        d_LC = 1.0
    end
    if has_right
        iR = face_R.cell_ids[1] == i ? face_R.cell_ids[2] : face_R.cell_ids[1]
        d_R_face = abs(state.face_positions[f_right] - state.cell_centers[i])
        d_CR = abs(state.cell_centers[iR] - state.cell_centers[i])
    else
        d_R_face = 1.0
        d_CR = 1.0
    end
    return (
        d_L_face = Float64(d_L_face),
        d_R_face = Float64(d_R_face),
        d_LC = Float64(d_LC),
        d_CR = Float64(d_CR),
        A_left = Float64(A_left),
        A_right = Float64(A_right),
        V_C = Float64(V_C),
        has_left = has_left,
        has_right = has_right,
    )
end

"""
	apply!(dT, op::ConductionOperator, state, mesh, material) -> nothing

Operator-level conduction contribution. Iterates over cells and accumulates
the per-cell `-div_q_cond / ρcp` into `dT`. Used by the standalone unit
test for residual-slice equivalence.
"""
function apply!(
    dT::AbstractVector,
    ::ConductionOperator,
    state,
    mesh,
    material::Material{NC},
) where {NC}
    out = zeros(eltype(dT), 1)
    n = length(state.T)
    @inbounds for i = 1:n
        geom = _conduction_geometry(state, mesh, i)
        iL = i == 1 ? i : i - 1
        iR = i == n ? i : i + 1
        T_L = state.T[iL]
        T_C = state.T[i]
        T_R = state.T[iR]
        ξ_L = ntuple(j -> state.ξ[j, iL], Val(NC))
        ξ_C = ntuple(j -> state.ξ[j, i], Val(NC))
        ξ_R = ntuple(j -> state.ξ[j, iR], Val(NC))
        apply_cell!(
            out,
            ConductionOperator(),
            T_L,
            T_C,
            T_R,
            ξ_L,
            ξ_C,
            ξ_R,
            geom.d_L_face,
            geom.d_R_face,
            geom.d_LC,
            geom.d_CR,
            geom.A_left,
            geom.A_right,
            geom.V_C,
            geom.has_left,
            geom.has_right,
            material,
        )
        dT[i] += out[1]
    end
    return nothing
end

# ============================================================================
# Analytical override — hand-derived chain rule
#
# Computes the same Jacobian entries as `_contribute_to_jacobian_ad!` below
# (the AD-default), but without ForwardDiff over the full 3-cell stencil.
# Property derivatives `∂k/∂T, ∂k/∂ξ, ∂ρcp/∂T, ∂ρcp/∂ξ` are obtained by
# evaluating `effective_conductivity` / `effective_heat_capacity` once per
# cell under stack-allocated `Dual{Tag, Float64, NC+1}` (input dim = NC+1
# vs. 3·(NC+1) for the AD-default — so the dual arithmetic is much cheaper).
# Face-flux composition, divergence, and the `1/ρcp` chain rule are then
# done by hand.
#
# Equivalence with the AD-default is gated at rtol = 1e-12 in
# `test/jacobian/test_operator_conduction_analytical.jl`. Only floating-
# point summation order differs between the two paths; the math is the
# same.
# ============================================================================

# Tag for the per-cell property-derivative duals — keeps these duals from
# colliding with any outer AD that may be in flight (none today, but
# defensive).
struct _ConductionPropTag end

# Seed `(T, ξ_1, ..., ξ_NC)` as `Dual{Tag, Float64, NC+1}` with unit partials.
# `@generated` so the partial-tuple construction unrolls at compile time and
# the result is stack-allocated. Returns `(T_d, ξ_d::NTuple{NC, Dual})`.
# Generic over the tag so other operators (e.g. the ALE strain-rate pass)
# can run the same stack-Dual sweep under their own tag.
@generated function _stack_dual_seed(
    ::Type{Tag},
    T::Float64,
    ξ::NTuple{NC,Float64},
) where {Tag,NC}

    N = NC + 1   # static — `NC` is a type parameter
    T_partials = Expr(:tuple, [i == 1 ? :(1.0) : :(0.0) for i = 1:N]...)
    # ξ_d as a tuple of N expressions, each constructing a Dual with j-th seed.
    ξ_exprs = Expr[]
    for j = 1:NC
        jp = j + 1
        ξ_partials = Expr(:tuple, [i == jp ? :(1.0) : :(0.0) for i = 1:N]...)
        push!(
            ξ_exprs,
            quote
                ForwardDiff.Dual{Tag,Float64,$N}(
                    ξ[$j],
                    ForwardDiff.Partials{$N,Float64}($ξ_partials),
                )
            end,
        )
    end
    ξ_tuple = Expr(:tuple, ξ_exprs...)
    return quote
        T_d = ForwardDiff.Dual{Tag,Float64,$N}(
            T,
            ForwardDiff.Partials{$N,Float64}($T_partials),
        )
        ξ_d = $ξ_tuple
        return T_d, ξ_d
    end
end

@inline _cond_seed_inputs(T::Float64, ξ::NTuple{NC,Float64}) where {NC} =
    _stack_dual_seed(_ConductionPropTag, T, ξ)

# Extract value + partials of a Dual{Tag, Float64, NC+1} into
# `(value, ∂/∂T, NTuple{NC, Float64})` — the partials in `(T, ξ_1..ξ_NC)`
# order. `@generated` so the `ntuple` unrolls.
@generated function _cond_extract(
    d::ForwardDiff.Dual{_ConductionPropTag,Float64,N},
    ::Val{NC},
) where {N,NC}
    N == NC + 1 || error("_cond_extract: N=$N must equal NC+1=$(NC+1)")
    ξ_exprs = [:(p[$(j+1)]) for j = 1:NC]
    ξ_tuple = Expr(:tuple, ξ_exprs...)
    return quote
        v = ForwardDiff.value(d)
        p = ForwardDiff.partials(d)
        return v, p[1], $ξ_tuple
    end
end

"""
	_cond_props_at_cell(ξ::NTuple{NC, Float64}, T::Float64, material)
	    -> (k, dk_dT, dk_dξ::NTuple{NC, Float64},
	        ρcp, dρcp_dT, dρcp_dξ::NTuple{NC, Float64})

Evaluate `effective_conductivity` and `effective_heat_capacity` at
`(T, ξ)` and read off their derivatives in one Dual-mode pass.
Stack-allocated.
"""
@inline function _cond_props_at_cell(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _cond_seed_inputs(T, ξ)
    k_d = effective_conductivity(ξ_d, T_d, material)
    ρcp_d = effective_heat_capacity(ξ_d, T_d, material)
    k, dk_dT, dk_dξ = _cond_extract(k_d, Val(NC))
    ρcp, dρcp_dT, dρcp_dξ = _cond_extract(ρcp_d, Val(NC))
    return (k, dk_dT, dk_dξ, ρcp, dρcp_dT, dρcp_dξ)
end

"""
	_cond_k_with_derivs(ξ::NTuple{NC, Float64}, T::Float64, material)
	    -> (k, dk_dT, dk_dξ::NTuple{NC, Float64})

Variant for neighbor cells — only conductivity derivatives needed
(no ρcp). Cheaper because `effective_heat_capacity` is skipped.
"""
@inline function _cond_k_with_derivs(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _cond_seed_inputs(T, ξ)
    k_d = effective_conductivity(ξ_d, T_d, material)
    return _cond_extract(k_d, Val(NC))
end

"""
	_harmonic_mean_with_derivs(k_L, k_R, d_L, d_R)
	    -> (k_face, dk_dkL, dk_dkR)

Distance-weighted harmonic mean and its derivatives w.r.t. `k_L` and
`k_R`. This is the derivative companion of
`Properties.harmonic_mean_weighted` and must keep the **identical**
`1e-20` guard band: a different cutoff makes the Jacobian disagree with the
residual whenever a conductivity falls between the two thresholds.

`k_face = (d_L + d_R) / (d_L/k_L + d_R/k_R)`,
`∂k_face/∂k_i = (d_L + d_R) · d_i / (k_i² · R_sum²)` where `R_sum = d_L/k_L + d_R/k_R`.
Returns zeros if either k is below `1e-20` (the residual's face value is 0
there, with zero sensitivity to either side).
"""
@inline function _harmonic_mean_with_derivs(
    k_L::Float64,
    k_R::Float64,
    d_L::Float64,
    d_R::Float64,
)
    if k_L < 1e-20 || k_R < 1e-20
        return 0.0, 0.0, 0.0
    end
    R_L = d_L / k_L
    R_R = d_R / k_R
    R_sum = R_L + R_R
    d_total = d_L + d_R
    k_face = d_total / R_sum
    denom = (R_sum * R_sum)
    dk_dkL = d_total * d_L / (k_L * k_L * denom)
    dk_dkR = d_total * d_R / (k_R * k_R * denom)
    return k_face, dk_dkL, dk_dkR
end

"""
	_contribute_to_jacobian_analytical!(jac, op::ConductionOperator, u, ctx, t,
	                                     ws, state, mesh, material, layout)

Analytical override of `contribute_to_jacobian!` for
`ConductionOperator`.

Iterates cells; for each cell computes `(k_C, ρcp_C)` plus their
derivatives w.r.t. `(T_C, ξ_C)` once via stack-Dual; for each existing
neighbor computes `(k_N, dk_N/dT, dk_N/dξ)` via the same path. Then for
each face composes the harmonic-mean chain rule, the per-face
`q_face = -k_face · ΔT/d` derivatives, and finally divides through by
`ρcp` with the `−q · ∂ρcp/∂x / ρcp²` correction on the cell-center
columns.

Workspace is the same `ConductionADWorkspace{NC}` as the AD-default —
the only field used is `ws.x_buf` (as scratch for stencil packing during
debug); workspace fields could be slimmed in a follow-up.
"""
function _contribute_to_jacobian_analytical!(
    jac::StructuredJacobian,
    ::ConductionOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::ConductionADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap

    # Stable views of the per-cell scratch arrays — allocated once at
    # workspace-build time, reused per cell.
    dQ_dξ_L = ws.dQ_dξ_L_buf
    dQ_dξ_C = ws.dQ_dξ_C_buf
    dQ_dξ_R = ws.dQ_dξ_R_buf

    @inbounds for i = 1:mesh.n_active_cells
        geom = _conduction_geometry(state, mesh, i)
        has_left = geom.has_left
        has_right = geom.has_right

        iL = has_left ? (i - 1) : 0
        iR = has_right ? (i + 1) : 0

        # Cell-center state and properties (k, ρcp, plus their (T, ξ) derivs).
        T_C = u[i]
        ξ_C = ntuple(j -> u[n_cells+(j-1)*n_cells+i], Val(NC))
        (k_C, dk_C_dT, dk_C_dξ, ρcp, dρcp_dT, dρcp_dξ) =
            _cond_props_at_cell(ξ_C, T_C, material)
        inv_ρcp = 1.0 / ρcp
        inv_ρcp_sq = inv_ρcp * inv_ρcp

        # Reset per-cell ξ accumulators.
        for j = 1:NC
            dQ_dξ_L[j] = 0.0
            dQ_dξ_C[j] = 0.0
            dQ_dξ_R[j] = 0.0
        end
        dQ_dT_L = 0.0
        dQ_dT_C = 0.0
        dQ_dT_R = 0.0
        Q = 0.0

        # ----------- Left face (between iL and i) -----------
        if has_left
            T_L = u[iL]
            ξ_L = ntuple(j -> u[n_cells+(j-1)*n_cells+iL], Val(NC))
            (k_L, dk_L_dT, dk_L_dξ) = _cond_k_with_derivs(ξ_L, T_L, material)

            d1 = geom.d_L_face
            d2 = geom.d_LC - geom.d_L_face
            (k_face, dkf_dkL, dkf_dkR) = _harmonic_mean_with_derivs(k_L, k_C, d1, d2)

            ΔT_over_d = (T_C - T_L) / geom.d_LC
            q_left = -k_face * ΔT_over_d
            coef_AV = geom.A_left / geom.V_C
            Q += q_left * coef_AV

            inv_d = 1.0 / geom.d_LC
            dq_dT_L = (-k_face) * (-inv_d) - ΔT_over_d * (dkf_dkL * dk_L_dT)
            dQ_dT_L += dq_dT_L * coef_AV
            dq_dT_C = (-k_face) * (+inv_d) - ΔT_over_d * (dkf_dkR * dk_C_dT)
            dQ_dT_C += dq_dT_C * coef_AV

            factor = -ΔT_over_d * coef_AV
            for j = 1:NC
                dQ_dξ_L[j] += factor * dkf_dkL * dk_L_dξ[j]
                dQ_dξ_C[j] += factor * dkf_dkR * dk_C_dξ[j]
            end
        end

        # ----------- Right face (between i and iR) -----------
        if has_right
            T_R = u[iR]
            ξ_R = ntuple(j -> u[n_cells+(j-1)*n_cells+iR], Val(NC))
            (k_R, dk_R_dT, dk_R_dξ) = _cond_k_with_derivs(ξ_R, T_R, material)

            d1 = geom.d_R_face
            d2 = geom.d_CR - geom.d_R_face
            (k_face, dkf_dkL, dkf_dkR) = _harmonic_mean_with_derivs(k_C, k_R, d1, d2)

            ΔT_over_d = (T_R - T_C) / geom.d_CR
            q_right = -k_face * ΔT_over_d
            coef_AV = geom.A_right / geom.V_C
            Q -= q_right * coef_AV

            inv_d = 1.0 / geom.d_CR
            # Right-face is subtracted into Q (sign flipped on contributions):
            dq_dT_C = (-k_face) * (-inv_d) - ΔT_over_d * (dkf_dkL * dk_C_dT)
            dQ_dT_C += -(dq_dT_C * coef_AV)
            dq_dT_R_flux = (-k_face) * (+inv_d) - ΔT_over_d * (dkf_dkR * dk_R_dT)
            dQ_dT_R += -(dq_dT_R_flux * coef_AV)

            # ξ contribution from right face is also negated (Q_R term has −sign).
            factor = +ΔT_over_d * coef_AV   # = −(−ΔT_over_d) · coef_AV
            for j = 1:NC
                dQ_dξ_C[j] += factor * dkf_dkL * dk_C_dξ[j]
                dQ_dξ_R[j] += factor * dkf_dkR * dk_R_dξ[j]
            end
        end

        # Convert (Q, ∂Q/∂x) → (out, ∂out/∂x) via 1/ρcp chain rule and scatter.
        # Only T_C and ξ_C carry the −Q·∂ρcp/∂x/ρcp² correction.
        if has_left
            fast_sparse_add!(A, smap, i, iL, dQ_dT_L * inv_ρcp)
        end
        fast_sparse_add!(A, smap, i, i, dQ_dT_C * inv_ρcp - Q * dρcp_dT * inv_ρcp_sq)
        if has_right
            fast_sparse_add!(A, smap, i, iR, dQ_dT_R * inv_ρcp)
        end

        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            if has_left
                fast_sparse_add!(A, smap, i, ξ_off + iL, dQ_dξ_L[j] * inv_ρcp)
            end
            fast_sparse_add!(
                A,
                smap,
                i,
                ξ_off + i,
                dQ_dξ_C[j] * inv_ρcp - Q * dρcp_dξ[j] * inv_ρcp_sq,
            )
            if has_right
                fast_sparse_add!(A, smap, i, ξ_off + iR, dQ_dξ_R[j] * inv_ρcp)
            end
        end
    end
    return nothing
end

# ============================================================================
# AD-default path (preserved as `_contribute_to_jacobian_ad!` for the
# rtol=1e-12 equivalence test against the analytical override).
# ============================================================================

"""
	contribute_to_jacobian!(jac, op::ConductionOperator,
	                        u, ctx, t, ws, state, mesh, material, layout) -> nothing

Default per-cell AD path. For each cell `i`:

1. Pack `(T_{i-1}, T_i, T_{i+1}, ξ_{:,i-1}, ξ_{:,i}, ξ_{:,i+1})` into
   `ws.x_buf` (with self-fill for missing neighbors at boundaries).
2. Run `ForwardDiff.jacobian!` over `apply_cell!` to fill `ws.J_local`
   (a `1 × (3 · (1 + NC))` row vector).
3. Scatter `J_local` into `jac.A`'s row `i` (T row) at the columns for
   the 3-cell window — but only for cells that actually exist (the
   self-filled boundary slot is *not* scattered, since adding to the
   cell's own column twice is wrong).
"""
function _contribute_to_jacobian_ad!(
    jac::StructuredJacobian,
    ::ConductionOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::ConductionADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}

    n_cells = layout.n_cells
    n_live = mesh.n_active_cells
    A = jac.A
    smap = jac.smap

    @inbounds for i = 1:n_live
        geom = _conduction_geometry(state, mesh, i)

        # Capture geometry + material constants in the closure
        f! = let geom = geom, material = material, NC = NC
            (out, x_inner) -> begin
                T_L, T_C, T_R, ξ_L, ξ_C, ξ_R = _unpack_stencil(x_inner, NC)
                apply_cell!(
                    out,
                    ConductionOperator(),
                    T_L,
                    T_C,
                    T_R,
                    ξ_L,
                    ξ_C,
                    ξ_R,
                    geom.d_L_face,
                    geom.d_R_face,
                    geom.d_LC,
                    geom.d_CR,
                    geom.A_left,
                    geom.A_right,
                    geom.V_C,
                    geom.has_left,
                    geom.has_right,
                    material,
                )
            end
        end

        _pack_stencil!(ws.x_buf, u, n_cells, NC, i)
        ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

        # Scatter. J_local is 1 × stencil_dim; row index = i (T row).
        # Column mapping:
        #   x[1] → T_iL; x[2] → T_i; x[3] → T_iR
        #   x[3 + (k-1)*NC + j] → ξ_{j, i_k} for k ∈ {L=1, C=2, R=3}, j ∈ 1:NC
        # Skip self-filled boundary slots (iL == i for i==1, iR == i for i==n).
        iL = i == 1 ? 0 : i - 1   # 0 = no left neighbor; skip the L column
        iR = i == n_live ? 0 : i + 1

        # T columns
        if iL != 0
            fast_sparse_add!(A, smap, i, iL, ws.J_local[1, 1])
        end
        fast_sparse_add!(A, smap, i, i, ws.J_local[1, 2])
        if iR != 0
            fast_sparse_add!(A, smap, i, iR, ws.J_local[1, 3])
        end
        # ξ columns
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            if iL != 0
                fast_sparse_add!(A, smap, i, ξ_off + iL, ws.J_local[1, 3+j])
            end
            fast_sparse_add!(A, smap, i, ξ_off + i, ws.J_local[1, 3+NC+j])
            if iR != 0
                fast_sparse_add!(A, smap, i, ξ_off + iR, ws.J_local[1, 3+2*NC+j])
            end
        end
    end
    return nothing
end

# ============================================================================
# Public dispatch
# ============================================================================
#
# `contribute_to_jacobian!` is the public method consumed by the Structured
# backend's operator-list walk. It dispatches to the analytical override —
# the AD-default `_contribute_to_jacobian_ad!` is preserved alongside as the
# rtol=1e-12 equivalence reference.

function contribute_to_jacobian!(
    jac::StructuredJacobian,
    op::ConductionOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::ConductionADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}
    return _contribute_to_jacobian_analytical!(
        jac,
        op,
        u,
        ctx,
        t,
        ws,
        state,
        mesh,
        material,
        layout,
    )
end
