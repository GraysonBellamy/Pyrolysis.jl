"""
The `DarcyTransportOperator` — the Darcy-Fick gas-transport contribution extracted from the residual.

Mirrors [`GasTransportOperator`](@ref) but with the face flux formed from the
combined advective-diffusive Darcy-Fick rule. Per cell:

```
dT[i]    += S_conv(T, ξ) / ρcp_eff(T_i, ξ_i)
dξ[j, i] += -div_J_gas[j, i]            (only for j with `is_gas == true`)
```

with `J_gas[j, f] = gas_flux_darcy_fick(v_g_face, ξ_L[j], ξ_R[j], φ_L, φ_R,
λ_face, T_L, T_R, d_LR)` (the advected scalar is the intrinsic density
ξ/φ_face). The Darcy velocity `v_g_face` is reconstructed inside the
AD closure from cell-local `(P, κ, μ)`:

- `P_i = ideal_gas_pressure(ξ_i, T_i, φ_i, material)` with `φ_i =
  porosity(ξ_i, T_i, material)`,
- `κ_i = effective_permeability(ξ_i, T_i, material)`,
- `μ_i = mixture_gas_viscosity(ξ_i, T_i, material)` (the material's own curve;
  Sutherland air by default),
- `mob_face` = distance-weighted harmonic mean of `κ/μ`,
- `v_g_face = -mob_face · (P_R − P_L)/d_LR`.

Species-mass scaling is **deferred** to a later refinement.

# Stencil

3-cell tridiagonal: `(cells = -1:1,)`. Pressure is cell-local on `(T, ξ)`,
so the chain rule through `:gas_pressure` stays inside the existing 3-cell
support. `reads = (:gas_pressure,)` is declared as forward-compat metadata —
no operator currently registers as the provider, so the AD-default
`contribute_to_jacobian!` falls through to local stencil AD.

# Boundary handling

Boundary-face fluxes are zeroed (`has_left` / `has_right` gates), matching
[`GasTransportOperator`](@ref), and boundary cells contribute **no S_conv**
here — the whole boundary-cell S_conv energy row is owned by
[`BCConvectionOperator`](@ref). The boundary mass flux — the diffusive film
law AND the pressure-BC advective venting, on the species rows — belongs
to [`BCMassOperator`](@ref).

# Sharp-feature smoothness

- Upwind switch `v_g >= 0 ? ξ_L : ξ_R`: non-smooth at `v_g = 0`.
  ForwardDiff returns the right-derivative; documented and accepted (mirrors
  the vanleer kink in ALE advection).
- Mobility's `μ > 1e-20` and `mob > 1e-30` guards: one-sided derivative
  across the threshold; thresholds are well below physical floors so this
  never triggers in practice.
- Pressure formula floors `φ ≥ 1e-12` via `max(φ, _PHI_FLOOR)`. Test states
  must satisfy `Σⱼ ξⱼ_solid / ρⱼ_intrinsic ≤ 1` — otherwise φ saturates and
  the Jacobian disagrees with the residual at ~1e25.
"""

using ForwardDiff
using ..Materials: Material, is_gas
using ..Properties:
    effective_gas_transfer,
    effective_heat_capacity,
    effective_permeability,
    effective_gas_viscosity,
    mixture_gas_viscosity,
    porosity,
    harmonic_mean_weighted
using ..Physics: gas_flux_darcy_fick, ideal_gas_pressure

"""
	DarcyTransportOperator <: Operator

Stateless tag; per-problem AD workspace lives on the backend cache.
"""
struct DarcyTransportOperator <: Operator end

"""
	DarcyTransportADWorkspace{NC, CFG}

Workspace for `DarcyTransportOperator.contribute_to_jacobian!`. Carries the
AD-default's stencil-AD plumbing (`x_buf`, `f_buf`, `J_local`, `cfg`) plus
the analytical override's per-cell scratch arrays (`dS_dξ_L_buf`, `_C_buf`,
`_R_buf`) for the T-row's ξ-derivative accumulators. The scratch arrays
replace per-face `ntuple(j -> ..., Val(NC))` closures that would heap-allocate.
"""
mutable struct DarcyTransportADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    dS_dξ_L_buf::Vector{Float64}
    dS_dξ_C_buf::Vector{Float64}
    dS_dξ_R_buf::Vector{Float64}
end

function build_darcy_transport_workspace(::Material{NC}) where {NC}
    stencil_dim = 3 * (1 + NC)
    out_dim = 1 + NC
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, out_dim)
    J_local = zeros(Float64, out_dim, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return DarcyTransportADWorkspace{NC,typeof(cfg)}(
        x_buf,
        f_buf,
        J_local,
        cfg,
        zeros(Float64, NC),
        zeros(Float64, NC),
        zeros(Float64, NC),
    )
end

# Pack stencil-local (T, ξ) for cell i — identical layout to GasTransportOperator.
@inline function _dt_pack_stencil!(
    x::AbstractVector,
    u::AbstractVector,
    n_cells::Int,
    NC::Int,
    i::Int,
)
    iL = i == 1 ? i : i - 1
    iR = i == n_cells ? i : i + 1
    @inbounds begin
        x[1] = u[iL]
        x[2] = u[i]
        x[3] = u[iR]
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            x[3+j] = u[ξ_off+iL]
            x[3+NC+j] = u[ξ_off+i]
            x[3+2*NC+j] = u[ξ_off+iR]
        end
    end
    return nothing
end

@inline function _dt_unpack_stencil(x::AbstractVector, NC::Int)
    @inbounds T_L = x[1]
    @inbounds T_C = x[2]
    @inbounds T_R = x[3]
    ξ_L = ntuple(j -> @inbounds(x[3+j]), NC)
    ξ_C = ntuple(j -> @inbounds(x[3+NC+j]), NC)
    ξ_R = ntuple(j -> @inbounds(x[3+2*NC+j]), NC)
    return T_L, T_C, T_R, ξ_L, ξ_C, ξ_R
end

# Distance-weighted harmonic mean of two mobilities (matches
# `update_darcy_properties!`): an impermeable side is an infinite series
# resistance, so the face mobility is exactly 0 when either side is ≈ 0.
@inline function _dt_mobility_face(
    κ_L::Real,
    κ_R::Real,
    μ_L::Real,
    μ_R::Real,
    d_L::Real,
    d_R::Real,
)
    mob_L = μ_L > 1e-20 ? κ_L / μ_L : zero(promote_type(typeof(κ_L), typeof(μ_L)))
    mob_R = μ_R > 1e-20 ? κ_R / μ_R : zero(promote_type(typeof(κ_R), typeof(μ_R)))
    if mob_L > 1e-30 && mob_R > 1e-30
        return (d_L + d_R) / (d_L / mob_L + d_R / mob_R)
    else
        return zero(promote_type(typeof(mob_L), typeof(mob_R)))
    end
end

# Per-component face-flux for the Darcy-Fick path. Non-gas components zero.
@generated function _dt_face_fluxes(
    T_L::Real,
    T_C::Real,
    ξ_L::NTuple{NC,Tu},
    ξ_C::NTuple{NC,Tu},
    φ_L::Real,
    φ_C::Real,
    λ_face::Real,
    v_g::Real,
    d_LR::Real,
    material::Material{NC},
) where {NC,Tu}
    exprs = Expr[]
    for j = 1:NC
        push!(
            exprs,
            quote
                @inbounds if is_gas(material.components[$j])
                    gas_flux_darcy_fick(
                        v_g,
                        ξ_L[$j],
                        ξ_C[$j],
                        φ_L,
                        φ_C,
                        λ_face,
                        T_L,
                        T_C,
                        d_LR,
                    )
                else
                    zero(Tu)
                end
            end,
        )
    end
    return :(NTuple{NC,Tu}(($(exprs...),)))
end

"""
	apply_cell!(out, op::DarcyTransportOperator,
	            T_L, T_C, T_R, ξ_L, ξ_C, ξ_R, geom, has_left, has_right,
	            material) -> nothing

Per-cell Darcy-Fick contribution. Writes:
- `out[1]`     : `S_conv(T, ξ) / ρcp_eff(T_C, ξ_C)`.
- `out[1 + j]` : `-div_J_gas[j]` for gas component `j`; `0` otherwise.
"""
function apply_cell!(
    out::AbstractVector,
    ::DarcyTransportOperator,
    T_L::Real,
    T_C::Real,
    T_R::Real,
    ξ_L::NTuple{NC,Tu},
    ξ_C::NTuple{NC,Tu},
    ξ_R::NTuple{NC,Tu},
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
) where {NC,Tu}

    # Cell-local properties at C.
    λ_C = effective_gas_transfer(ξ_C, T_C, material)
    ρcp = effective_heat_capacity(ξ_C, T_C, material)

    # Pressure / permeability / viscosity at C — needed for both faces.
    φ_C = porosity(ξ_C, T_C, material)
    P_C = ideal_gas_pressure(ξ_C, T_C, φ_C, material)
    κ_C = effective_permeability(ξ_C, T_C, material)
    μ_C = mixture_gas_viscosity(ξ_C, T_C, material)

    # Left face flux per component.
    J_left = if has_left
        λ_L = effective_gas_transfer(ξ_L, T_L, material)
        λ_face_L = harmonic_mean_weighted(λ_L, λ_C, d_L_face, d_LC - d_L_face)
        φ_L = porosity(ξ_L, T_L, material)
        P_L = ideal_gas_pressure(ξ_L, T_L, φ_L, material)
        κ_L = effective_permeability(ξ_L, T_L, material)
        μ_L = mixture_gas_viscosity(ξ_L, T_L, material)
        mob_face = _dt_mobility_face(κ_L, κ_C, μ_L, μ_C, d_L_face, d_LC - d_L_face)
        v_g_L = -mob_face * (P_C - P_L) / d_LC
        _dt_face_fluxes(T_L, T_C, ξ_L, ξ_C, φ_L, φ_C, λ_face_L, v_g_L, d_LC, material)
    else
        ntuple(_ -> zero(Tu), Val(NC))
    end

    # Right face flux per component.
    J_right = if has_right
        λ_R = effective_gas_transfer(ξ_R, T_R, material)
        λ_face_R = harmonic_mean_weighted(λ_C, λ_R, d_R_face, d_CR - d_R_face)
        φ_R = porosity(ξ_R, T_R, material)
        P_R = ideal_gas_pressure(ξ_R, T_R, φ_R, material)
        κ_R = effective_permeability(ξ_R, T_R, material)
        μ_R = mixture_gas_viscosity(ξ_R, T_R, material)
        mob_face = _dt_mobility_face(κ_C, κ_R, μ_C, μ_R, d_R_face, d_CR - d_R_face)
        v_g_R = -mob_face * (P_R - P_C) / d_CR
        _dt_face_fluxes(T_C, T_R, ξ_C, ξ_R, φ_C, φ_R, λ_face_R, v_g_R, d_CR, material)
    else
        ntuple(_ -> zero(Tu), Val(NC))
    end

    # dT/dx at cell center: opposite-distance-weighted slopes (second order
    # on nonuniform meshes), mirroring the interior branch of
    # `Physics.compute_gas_convective_source!`. At boundary cells the WHOLE
    # S_conv energy row — boundary-face flux × T_s-aware gradient — is owned
    # by `BCConvectionOperator`, so this operator contributes no S_conv there
    # (dTdx = 0 zeroes it; the ξ-rows stay here).
    dTdx = if has_left && has_right
        opposite_distance_gradient(T_L, T_C, T_R, d_LC, d_CR)
    else
        zero(T_C)
    end

    S_conv = _dt_S_conv_unrolled(T_C, dTdx, J_left, J_right, material)

    @inbounds out[1] = S_conv / ρcp
    _dt_species_div_unrolled!(out, J_left, J_right, A_left, A_right, V_C, material)
    return nothing
end

@generated function _dt_S_conv_unrolled(
    T_C::Real,
    dTdx::Real,
    J_left::NTuple{NC,Tu},
    J_right::NTuple{NC,Tu},
    material::Material{NC},
) where {NC,Tu}
    exprs = Expr[]
    push!(exprs, :(S = zero(promote_type(typeof(T_C), typeof(dTdx), Tu))))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                J_avg = (J_left[$j] + J_right[$j]) / 2
                cp_j = material.components[$j].heat_capacity(T_C)
                S -= cp_j * J_avg * dTdx
            end
        end)
    end
    push!(exprs, :(return S))
    return Expr(:block, exprs...)
end

@generated function _dt_species_div_unrolled!(
    out::AbstractVector,
    J_left::NTuple{NC,Tu},
    J_right::NTuple{NC,Tu},
    A_left::Real,
    A_right::Real,
    V_C::Real,
    material::Material{NC},
) where {NC,Tu}
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                div_J = (J_right[$j] * A_right - J_left[$j] * A_left) / V_C
                out[1+$j] = -div_J
            else
                out[1+$j] = zero(Tu)
            end
        end)
    end
    push!(exprs, :(return nothing))
    return Expr(:block, exprs...)
end

@inline function _dt_geometry(state, mesh, i::Int)
    cell = mesh.cells[i]
    f_left, f_right = cell.face_ids
    face_L = mesh.faces[f_left]
    face_R = mesh.faces[f_right]
    has_left = !(face_L.cell_ids[1] == 0 || face_L.cell_ids[2] == 0)
    has_right = !(face_R.cell_ids[1] == 0 || face_R.cell_ids[2] == 0)

    A_left = state.A_face[f_left]
    A_right = state.A_face[f_right]
    V_C = state.volume[i]

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
	apply!(dT, dξ, op::DarcyTransportOperator, state, mesh, material) -> nothing

Operator-level loop accumulating per-cell `(S_conv / ρcp_eff, -div_J_gas[:])`
into `(dT, dξ)`. Boundary faces zero out (boundary mass flux is the
[`BCMassOperator`](@ref)'s domain).
"""
function apply!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    ::DarcyTransportOperator,
    state,
    mesh,
    material::Material{NC},
) where {NC}
    out = zeros(eltype(dT), 1 + NC)
    n = length(state.T)
    @inbounds for i = 1:n
        geom = _dt_geometry(state, mesh, i)
        iL = i == 1 ? i : i - 1
        iR = i == n ? i : i + 1
        Tu = eltype(state.T)
        ξ_L = ntuple(j -> state.ξ[j, iL], Val(NC))
        ξ_C = ntuple(j -> state.ξ[j, i], Val(NC))
        ξ_R = ntuple(j -> state.ξ[j, iR], Val(NC))
        apply_cell!(
            out,
            DarcyTransportOperator(),
            state.T[iL],
            state.T[i],
            state.T[iR],
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
        for j = 1:NC
            dξ[j, i] += out[1+j]
        end
    end
    return nothing
end

# ============================================================================
# Analytical override (DarcyTransport)
#
# Same Jacobian as the AD path below but with property derivatives extracted
# from one Dual-mode pass per cell (input dim = NC+1) — replacing the AD
# default's 3·(NC+1)-dim Dual chunk over the whole `apply_cell!`. Face fluxes
# follow the Darcy-Fick formula `J = v_g · ξ_upwind − λ_face · (ξ_R T_R − ξ_L T_L) /
# (T_face · d)`. The chain rule has three composed paths:
#   • λ-path:  λ_face is harmonic_mean(λ_L_cell, λ_R_cell) of cell-local λ.
#   • mobility-path:  mob_face is harmonic_mean(κ_L/μ_L, κ_R/μ_R) of cell-local
#     mobility, and v_g = -mob_face · ΔP / d.
#   • pressure-path:  P_L, P_R = ideal_gas_pressure(ξ, T, φ(ξ, T), …).
# `_harmonic_mean_with_derivs` (lifted in `operators/conduction.jl`) is reused
# for both the λ_face and mob_face derivatives. `μ` is the material's own
# `gas_viscosity` curve via `effective_gas_viscosity(T, material)` — `dμ/dT` from
# `ForwardDiff.derivative` on that scalar function, so any correlation works.
# ============================================================================

struct _DTPropTag end

@generated function _dt_seed_inputs(T::Float64, ξ::NTuple{NC,Float64}) where {NC}
    N = NC + 1
    T_p = Expr(:tuple, [i == 1 ? :(1.0) : :(0.0) for i = 1:N]...)
    ξ_exprs = Expr[]
    for j = 1:NC
        jp = j + 1
        p = Expr(:tuple, [i == jp ? :(1.0) : :(0.0) for i = 1:N]...)
        push!(
            ξ_exprs,
            quote
                ForwardDiff.Dual{_DTPropTag,Float64,$N}(
                    ξ[$j],
                    ForwardDiff.Partials{$N,Float64}($p),
                )
            end,
        )
    end
    return quote
        T_d = ForwardDiff.Dual{_DTPropTag,Float64,$N}(
            T,
            ForwardDiff.Partials{$N,Float64}($T_p),
        )
        ξ_d = $(Expr(:tuple, ξ_exprs...))
        return T_d, ξ_d
    end
end

@generated function _dt_extract(
    d::ForwardDiff.Dual{_DTPropTag,Float64,N},
    ::Val{NC},
) where {N,NC}
    N == NC + 1 || error("_dt_extract: N=$N must equal NC+1=$(NC+1)")
    ξ_t = Expr(:tuple, [:(p[$(j+1)]) for j = 1:NC]...)
    return quote
        v = ForwardDiff.value(d)
        p = ForwardDiff.partials(d)
        return v, p[1], $ξ_t
    end
end

# Cell-center properties: (λ, ρcp, P, κ, φ) and their (T, ξ)-derivatives in
# one Dual pass. P is computed via `ideal_gas_pressure(ξ, T, φ(ξ, T), …)` so
# the porosity chain rule is folded into `dP_dT, dP_dξ` automatically; φ and
# its derivatives are also extracted for the advective ξ/φ_face chain.
@inline function _dt_props_at_cell(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _dt_seed_inputs(T, ξ)
    λ_d = effective_gas_transfer(ξ_d, T_d, material)
    ρcp_d = effective_heat_capacity(ξ_d, T_d, material)
    φ_d = porosity(ξ_d, T_d, material)
    P_d = ideal_gas_pressure(ξ_d, T_d, φ_d, material)
    κ_d = effective_permeability(ξ_d, T_d, material)
    λ, dλ_dT, dλ_dξ = _dt_extract(λ_d, Val(NC))
    ρcp, dρcp_dT, dρcp_dξ = _dt_extract(ρcp_d, Val(NC))
    P, dP_dT, dP_dξ = _dt_extract(P_d, Val(NC))
    κ, dκ_dT, dκ_dξ = _dt_extract(κ_d, Val(NC))
    φ, dφ_dT, dφ_dξ = _dt_extract(φ_d, Val(NC))
    return (
        λ,
        dλ_dT,
        dλ_dξ,
        ρcp,
        dρcp_dT,
        dρcp_dξ,
        P,
        dP_dT,
        dP_dξ,
        κ,
        dκ_dT,
        dκ_dξ,
        φ,
        dφ_dT,
        dφ_dξ,
    )
end

# Neighbor-cell properties: (λ, P, κ, φ) only (ρcp not needed at neighbors).
@inline function _dt_neighbor_props(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _dt_seed_inputs(T, ξ)
    λ_d = effective_gas_transfer(ξ_d, T_d, material)
    φ_d = porosity(ξ_d, T_d, material)
    P_d = ideal_gas_pressure(ξ_d, T_d, φ_d, material)
    κ_d = effective_permeability(ξ_d, T_d, material)
    λ, dλ_dT, dλ_dξ = _dt_extract(λ_d, Val(NC))
    P, dP_dT, dP_dξ = _dt_extract(P_d, Val(NC))
    κ, dκ_dT, dκ_dξ = _dt_extract(κ_d, Val(NC))
    φ, dφ_dT, dφ_dξ = _dt_extract(φ_d, Val(NC))
    return (λ, dλ_dT, dλ_dξ, P, dP_dT, dP_dξ, κ, dκ_dT, dκ_dξ, φ, dφ_dT, dφ_dξ)
end

# Pore-gas viscosity + scalar derivative w.r.t. T, from the MATERIAL's own
# curve. Length-1 Dual via `ForwardDiff.derivative` so the derivative stays in
# lockstep with whatever correlation the material carries — no hand-coded
# constants here, and a user-supplied polynomial differentiates as correctly as
# the default Sutherland form.
#
# INVARIANT: μ is a function of T alone, so `dμ/dξ == 0` and the mobility path's
# ξ-partials below come only from κ. `mixture_gas_viscosity` ignores ξ by design
# (see its docstring). If that ever changes — the Wilke/option-B direction — this
# helper must also return `dμ_dξ::NTuple{NC}` and EVERY `dvg_dξ_*_j` expression
# below must gain a `− mob_*_cell * dμ_*_dξ[j]` term inside its `dmf` factor.
# Missing that leaves the Jacobian silently stale: the rtol = 1e-12
# analytical-vs-AD gate cannot see it while every Darcy fixture carries a single
# gas component, since a one-species mixture rule is the identity. See
# `docs/Gas_Viscosity_Recommendation_2026-09-21.md` §1.
@inline function _dt_viscosity_with_dT(T::Float64, material)
    μ = effective_gas_viscosity(T, material)
    dμ_dT = ForwardDiff.derivative(t -> effective_gas_viscosity(t, material), T)
    return μ, dμ_dT
end

# Direct partial derivatives of the Darcy-Fick face flux:
#   `J(λ_face, v_g, ξ_L, ξ_R, φ_L, φ_R, T_L, T_R, d) = v_g · ξ_upwind / φ_face
#     - λ_face · (ξ_R T_R − ξ_L T_L) / (T_face · d)`,
# with φ_face = max((φ_L + φ_R)/2, _PHI_FLOOR) — matches `gas_flux_darcy_fick`.
# Returns `(J, dJ/dλ, dJ/dv, dJ/dξ_L, dJ/dξ_R, dJ/dT_L, dJ/dT_R,
# dJ/dφ_L, dJ/dφ_R, adv_coef, upwind_is_left)` where
# `adv_coef = v_g/φ_face` is the advective direct partial w.r.t. the upwind ξ
# (added by the caller on the upwind side), `dJ/dφ_*` carries the ξ/φ_face²
# chain (zero when the floor binds — matching ForwardDiff through `max`), and
# `upwind_is_left = v_g >= 0` (matches the AD-default's right-derivative
# branch at v_g = 0).
@inline function _dt_J_face_partials(
    λ_face::Float64,
    v_g::Float64,
    ξ_L::Float64,
    ξ_R::Float64,
    φ_L::Float64,
    φ_R::Float64,
    T_L::Float64,
    T_R::Float64,
    d::Float64,
)
    T_face = 0.5 * (T_L + T_R)
    denom = T_face * d
    N_num = ξ_R * T_R - ξ_L * T_L
    J_diff = -λ_face * N_num / denom
    upwind_is_left = v_g >= 0.0
    ξ_upwind = upwind_is_left ? ξ_L : ξ_R
    φ_face_raw = 0.5 * (φ_L + φ_R)
    φ_face = max(φ_face_raw, 1e-12)  # _PHI_FLOOR, matches the pressure closure
    inv_φ_face = 1.0 / φ_face
    adv_coef = v_g * inv_φ_face
    J_adv = adv_coef * ξ_upwind
    J = J_adv + J_diff
    dJ_dλ = -N_num / denom
    dJ_dv = ξ_upwind * inv_φ_face
    # φ-chain: ∂J_adv/∂φ_face = −v_g·ξ_up/φ_face²; ∂φ_face/∂φ_L = ∂/∂φ_R = ½
    # (0 when the floor binds, like AD through max).
    dJ_dφ_common =
        φ_face_raw > 1e-12 ? -0.5 * v_g * ξ_upwind * inv_φ_face * inv_φ_face : 0.0
    # Diffusive direct partials (advective direct partials handled by caller via
    # the upwind branch).
    dJ_dξ_L = λ_face * T_L / denom
    dJ_dξ_R = -λ_face * T_R / denom
    dJ_dT_L = λ_face * ξ_L / denom - J_diff / (2 * T_face)
    dJ_dT_R = -λ_face * ξ_R / denom - J_diff / (2 * T_face)
    return J,
    dJ_dλ,
    dJ_dv,
    dJ_dξ_L,
    dJ_dξ_R,
    dJ_dT_L,
    dJ_dT_R,
    dJ_dφ_common,
    dJ_dφ_common,
    adv_coef,
    upwind_is_left
end

"""
	_contribute_to_jacobian_analytical!(jac, op::DarcyTransportOperator, …)

Analytical override of `contribute_to_jacobian!` for
`DarcyTransportOperator`. Iterates cells; for each cell evaluates
`(λ_C, ρcp, P_C, κ_C, μ_C, cp_g, all derivs)` once and assembles the
per-face Darcy-Fick flux chain rule plus the `S_conv/ρcp` and `div(J)`
divergences by hand.

Equivalence with the AD-default `_contribute_to_jacobian_ad!` is gated at
rtol = 1e-12 in `test/jacobian/test_operator_darcy_transport_analytical.jl`.
"""
function _contribute_to_jacobian_analytical!(
    jac::StructuredJacobian,
    ::DarcyTransportOperator,
    u::AbstractVector,
    t::Real,
    ws::DarcyTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap
    is_gas_t = _gt_is_gas_tuple(material.components, Val(NC))

    @inbounds for i = 1:mesh.n_active_cells
        geom = _dt_geometry(state, mesh, i)
        has_left = geom.has_left
        has_right = geom.has_right

        (!has_left && !has_right) && continue

        iL = has_left ? (i - 1) : 0
        iR = has_right ? (i + 1) : 0

        T_C = u[i]
        ξ_C = ntuple(j -> u[n_cells+(j-1)*n_cells+i], Val(NC))
        (
            λ_C,
            dλ_C_dT,
            dλ_C_dξ,
            ρcp,
            dρcp_dT,
            dρcp_dξ,
            P_C,
            dP_C_dT,
            dP_C_dξ,
            κ_C,
            dκ_C_dT,
            dκ_C_dξ,
            φ_C,
            dφ_C_dT,
            dφ_C_dξ,
        ) = _dt_props_at_cell(ξ_C, T_C, material)
        μ_C, dμ_C_dT = _dt_viscosity_with_dT(T_C, material)
        mob_C_cell = κ_C / μ_C
        inv_μ_C = 1.0 / μ_C
        inv_ρcp = 1.0 / ρcp
        inv_ρcp_sq = inv_ρcp * inv_ρcp

        (cp_t, dcp_dT_t) = _gt_cp_per_component(T_C, material.components, Val(NC))

        # dTdx and its (T_L, T_C, T_R) coefficients — same shape as
        # gas_transport: boundary cells have dTdx ≡ 0 here (their S_conv
        # energy row is owned whole by `BCConvectionOperator`), zeroing every
        # S_conv accumulator while the ξ-rows still assemble.
        if has_left && has_right
            T_L = u[iL]
            T_R = u[iR]
            dTdx = opposite_distance_gradient(T_L, T_C, T_R, geom.d_LC, geom.d_CR)
            (ddTdx_dT_L, ddTdx_dT_C, ddTdx_dT_R) =
                opposite_distance_gradient_coeffs(geom.d_LC, geom.d_CR)
        elseif has_right
            T_R = u[iR]
            T_L = T_C
            dTdx = 0.0
            ddTdx_dT_L = 0.0
            ddTdx_dT_C = 0.0
            ddTdx_dT_R = 0.0
        else
            T_L = u[iL]
            T_R = T_C
            dTdx = 0.0
            ddTdx_dT_L = 0.0
            ddTdx_dT_C = 0.0
            ddTdx_dT_R = 0.0
        end

        # T-row accumulators (per-cell scratch).
        dS_dT_L = 0.0
        dS_dT_C = 0.0
        dS_dT_R = 0.0
        dS_dξ_L_buf = ws.dS_dξ_L_buf
        dS_dξ_C_buf = ws.dS_dξ_C_buf
        dS_dξ_R_buf = ws.dS_dξ_R_buf
        for j = 1:NC
            dS_dξ_L_buf[j] = 0.0
            dS_dξ_C_buf[j] = 0.0
            dS_dξ_R_buf[j] = 0.0
        end
        ΣcgJavg = 0.0
        Σdcg_dT_Javg = 0.0

        # === Left face (between iL and i) ===
        if has_left
            ξ_L = ntuple(j -> u[n_cells+(j-1)*n_cells+iL], Val(NC))
            (
                λ_L,
                dλ_L_dT,
                dλ_L_dξ,
                P_L,
                dP_L_dT,
                dP_L_dξ,
                κ_L,
                dκ_L_dT,
                dκ_L_dξ,
                φ_L,
                dφ_L_dT,
                dφ_L_dξ,
            ) = _dt_neighbor_props(ξ_L, T_L, material)
            μ_L, dμ_L_dT = _dt_viscosity_with_dT(T_L, material)
            mob_L_cell = κ_L / μ_L
            inv_μ_L = 1.0 / μ_L

            d1 = geom.d_L_face
            d2 = geom.d_LC - geom.d_L_face

            (λ_face, dλf_dλL, dλf_dλC) = _harmonic_mean_with_derivs(λ_L, λ_C, d1, d2)
            (mob_face, dmf_dmL, dmf_dmC) =
                _harmonic_mean_with_derivs(mob_L_cell, mob_C_cell, d1, d2)

            inv_d = 1.0 / geom.d_LC
            dP = P_C - P_L
            v_g = -mob_face * dP * inv_d

            # v_g chain w.r.t. (T_L, T_C, ξ_L_j, ξ_C_j) via P-path and mob-path.
            # mob_L_cell = κ_L/μ_L, so dmob_L/dT = (dκ_L_dT - mob_L_cell · dμ_L_dT)/μ_L.
            dmob_L_dT_L = (dκ_L_dT - mob_L_cell * dμ_L_dT) * inv_μ_L
            dmob_C_dT_C = (dκ_C_dT - mob_C_cell * dμ_C_dT) * inv_μ_C
            coef_v_dPL = mob_face * inv_d
            coef_v_dPC = -mob_face * inv_d
            coef_v_dmf = -dP * inv_d

            dvg_dT_L = coef_v_dPL * dP_L_dT + coef_v_dmf * dmf_dmL * dmob_L_dT_L
            dvg_dT_C = coef_v_dPC * dP_C_dT + coef_v_dmf * dmf_dmC * dmob_C_dT_C

            coef_AV = geom.A_left / geom.V_C

            for g = 1:NC
                is_gas_t[g] || continue
                ξ_L_g = ξ_L[g]
                ξ_C_g = ξ_C[g]
                J,
                dJ_dλ,
                dJ_dv,
                dJ_dξL_diff,
                dJ_dξR_diff,
                dJ_dTL_diff,
                dJ_dTR_diff,
                dJ_dφL,
                dJ_dφC,
                adv_coef,
                upwind_is_left = _dt_J_face_partials(
                    λ_face,
                    v_g,
                    ξ_L_g,
                    ξ_C_g,
                    φ_L,
                    φ_C,
                    T_L,
                    T_C,
                    geom.d_LC,
                )

                cp_g = cp_t[g]
                dcp_g_dT = dcp_dT_t[g]

                # Compose face-flux derivatives w.r.t. stencil vars (λ-path,
                # mobility/pressure-path via v_g, and the φ_face-path of the
                # advected intrinsic density).
                dJ_dT_L_full =
                    dJ_dTL_diff +
                    dJ_dλ * dλf_dλL * dλ_L_dT +
                    dJ_dv * dvg_dT_L +
                    dJ_dφL * dφ_L_dT
                dJ_dT_C_full =
                    dJ_dTR_diff +
                    dJ_dλ * dλf_dλC * dλ_C_dT +
                    dJ_dv * dvg_dT_C +
                    dJ_dφC * dφ_C_dT

                # === ξ-row scatter for component g (left-face contributes
                # +J_L · A_L / V_C to dξ[g, i]).
                row_g = n_cells + (g - 1) * n_cells + i
                fast_sparse_add!(A, smap, row_g, iL, dJ_dT_L_full * coef_AV)
                fast_sparse_add!(A, smap, row_g, i, dJ_dT_C_full * coef_AV)
                for j = 1:NC
                    ξ_off = n_cells + (j - 1) * n_cells
                    # Pressure-path partials of v_g w.r.t. ξ_*_j (μ has no ξ
                    # dependence, so the mobility-path piece comes only from κ).
                    dvg_dξ_L_j =
                        coef_v_dPL * dP_L_dξ[j] +
                        coef_v_dmf * dmf_dmL * (inv_μ_L * dκ_L_dξ[j])
                    dvg_dξ_C_j =
                        coef_v_dPC * dP_C_dξ[j] +
                        coef_v_dmf * dmf_dmC * (inv_μ_C * dκ_C_dξ[j])

                    contrib_L =
                        dJ_dλ * dλf_dλL * dλ_L_dξ[j] +
                        dJ_dv * dvg_dξ_L_j +
                        dJ_dφL * dφ_L_dξ[j]
                    contrib_C =
                        dJ_dλ * dλf_dλC * dλ_C_dξ[j] +
                        dJ_dv * dvg_dξ_C_j +
                        dJ_dφC * dφ_C_dξ[j]
                    if j == g
                        contrib_L += dJ_dξL_diff
                        contrib_C += dJ_dξR_diff
                        # Advective direct upwind contribution: dJ_adv/dξ_{side, g}
                        # = v_g/φ_face on the upwind side, 0 elsewhere.
                        if upwind_is_left
                            contrib_L += adv_coef
                        else
                            contrib_C += adv_coef
                        end
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + iL, contrib_L * coef_AV)
                    fast_sparse_add!(A, smap, row_g, ξ_off + i, contrib_C * coef_AV)
                end

                # === T-row contribution from this face ===
                ΣcgJavg += cp_g * (J / 2)
                Σdcg_dT_Javg += dcp_g_dT * (J / 2)

                # Stencil-derivative factor of S_conv for this face.
                factor_T = -cp_g * 0.5 * dTdx

                dS_dT_L += factor_T * dJ_dT_L_full
                dS_dT_C += factor_T * dJ_dT_C_full

                for j = 1:NC
                    dvg_dξ_L_j =
                        coef_v_dPL * dP_L_dξ[j] +
                        coef_v_dmf * dmf_dmL * (inv_μ_L * dκ_L_dξ[j])
                    dvg_dξ_C_j =
                        coef_v_dPC * dP_C_dξ[j] +
                        coef_v_dmf * dmf_dmC * (inv_μ_C * dκ_C_dξ[j])
                    contrib_L =
                        dJ_dλ * dλf_dλL * dλ_L_dξ[j] +
                        dJ_dv * dvg_dξ_L_j +
                        dJ_dφL * dφ_L_dξ[j]
                    contrib_C =
                        dJ_dλ * dλf_dλC * dλ_C_dξ[j] +
                        dJ_dv * dvg_dξ_C_j +
                        dJ_dφC * dφ_C_dξ[j]
                    if j == g
                        contrib_L += dJ_dξL_diff
                        contrib_C += dJ_dξR_diff
                        if upwind_is_left
                            contrib_L += adv_coef
                        else
                            contrib_C += adv_coef
                        end
                    end
                    dS_dξ_L_buf[j] += factor_T * contrib_L
                    dS_dξ_C_buf[j] += factor_T * contrib_C
                end
            end
        end

        # === Right face (between i and iR) ===
        if has_right
            ξ_R = ntuple(j -> u[n_cells+(j-1)*n_cells+iR], Val(NC))
            (
                λ_R,
                dλ_R_dT,
                dλ_R_dξ,
                P_R,
                dP_R_dT,
                dP_R_dξ,
                κ_R,
                dκ_R_dT,
                dκ_R_dξ,
                φ_R,
                dφ_R_dT,
                dφ_R_dξ,
            ) = _dt_neighbor_props(ξ_R, T_R, material)
            μ_R, dμ_R_dT = _dt_viscosity_with_dT(T_R, material)
            mob_R_cell = κ_R / μ_R
            inv_μ_R = 1.0 / μ_R

            d1 = geom.d_R_face
            d2 = geom.d_CR - geom.d_R_face

            # face role: L-cell side = center (i); R-cell side = neighbor (iR).
            (λ_face, dλf_dλC, dλf_dλR) = _harmonic_mean_with_derivs(λ_C, λ_R, d1, d2)
            (mob_face, dmf_dmC, dmf_dmR) =
                _harmonic_mean_with_derivs(mob_C_cell, mob_R_cell, d1, d2)

            inv_d = 1.0 / geom.d_CR
            dP = P_R - P_C
            v_g = -mob_face * dP * inv_d

            dmob_C_dT_C = (dκ_C_dT - mob_C_cell * dμ_C_dT) * inv_μ_C
            dmob_R_dT_R = (dκ_R_dT - mob_R_cell * dμ_R_dT) * inv_μ_R
            coef_v_dPC = mob_face * inv_d
            coef_v_dPR = -mob_face * inv_d
            coef_v_dmf = -dP * inv_d

            dvg_dT_C = coef_v_dPC * dP_C_dT + coef_v_dmf * dmf_dmC * dmob_C_dT_C
            dvg_dT_R = coef_v_dPR * dP_R_dT + coef_v_dmf * dmf_dmR * dmob_R_dT_R

            coef_AV = geom.A_right / geom.V_C

            for g = 1:NC
                is_gas_t[g] || continue
                ξ_C_g = ξ_C[g]
                ξ_R_g = ξ_R[g]
                J,
                dJ_dλ,
                dJ_dv,
                dJ_dξL_diff,
                dJ_dξR_diff,
                dJ_dTL_diff,
                dJ_dTR_diff,
                dJ_dφC,
                dJ_dφR,
                adv_coef,
                upwind_is_left = _dt_J_face_partials(
                    λ_face,
                    v_g,
                    ξ_C_g,
                    ξ_R_g,
                    φ_C,
                    φ_R,
                    T_C,
                    T_R,
                    geom.d_CR,
                )

                cp_g = cp_t[g]
                dcp_g_dT = dcp_dT_t[g]

                dJ_dT_C_full =
                    dJ_dTL_diff +
                    dJ_dλ * dλf_dλC * dλ_C_dT +
                    dJ_dv * dvg_dT_C +
                    dJ_dφC * dφ_C_dT
                dJ_dT_R_full =
                    dJ_dTR_diff +
                    dJ_dλ * dλf_dλR * dλ_R_dT +
                    dJ_dv * dvg_dT_R +
                    dJ_dφR * dφ_R_dT

                # Right-face contributes -J_R · A_R / V_C to dξ[g, i].
                row_g = n_cells + (g - 1) * n_cells + i
                fast_sparse_add!(A, smap, row_g, i, -dJ_dT_C_full * coef_AV)
                fast_sparse_add!(A, smap, row_g, iR, -dJ_dT_R_full * coef_AV)
                for j = 1:NC
                    ξ_off = n_cells + (j - 1) * n_cells
                    dvg_dξ_C_j =
                        coef_v_dPC * dP_C_dξ[j] +
                        coef_v_dmf * dmf_dmC * (inv_μ_C * dκ_C_dξ[j])
                    dvg_dξ_R_j =
                        coef_v_dPR * dP_R_dξ[j] +
                        coef_v_dmf * dmf_dmR * (inv_μ_R * dκ_R_dξ[j])

                    contrib_C =
                        dJ_dλ * dλf_dλC * dλ_C_dξ[j] +
                        dJ_dv * dvg_dξ_C_j +
                        dJ_dφC * dφ_C_dξ[j]
                    contrib_R =
                        dJ_dλ * dλf_dλR * dλ_R_dξ[j] +
                        dJ_dv * dvg_dξ_R_j +
                        dJ_dφR * dφ_R_dξ[j]
                    if j == g
                        contrib_C += dJ_dξL_diff
                        contrib_R += dJ_dξR_diff
                        # For the right face, "left side" of the face is the
                        # center cell (i). upwind_is_left here means v_g >= 0,
                        # i.e., upwind = ξ_C.
                        if upwind_is_left
                            contrib_C += adv_coef
                        else
                            contrib_R += adv_coef
                        end
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + i, -contrib_C * coef_AV)
                    fast_sparse_add!(A, smap, row_g, ξ_off + iR, -contrib_R * coef_AV)
                end

                ΣcgJavg += cp_g * (J / 2)
                Σdcg_dT_Javg += dcp_g_dT * (J / 2)

                # Stencil-derivative factor of S_conv for this face.
                factor_T = -cp_g * 0.5 * dTdx

                dS_dT_C += factor_T * dJ_dT_C_full
                dS_dT_R += factor_T * dJ_dT_R_full

                for j = 1:NC
                    dvg_dξ_C_j =
                        coef_v_dPC * dP_C_dξ[j] +
                        coef_v_dmf * dmf_dmC * (inv_μ_C * dκ_C_dξ[j])
                    dvg_dξ_R_j =
                        coef_v_dPR * dP_R_dξ[j] +
                        coef_v_dmf * dmf_dmR * (inv_μ_R * dκ_R_dξ[j])
                    contrib_C =
                        dJ_dλ * dλf_dλC * dλ_C_dξ[j] +
                        dJ_dv * dvg_dξ_C_j +
                        dJ_dφC * dφ_C_dξ[j]
                    contrib_R =
                        dJ_dλ * dλf_dλR * dλ_R_dξ[j] +
                        dJ_dv * dvg_dξ_R_j +
                        dJ_dφR * dφ_R_dξ[j]
                    if j == g
                        contrib_C += dJ_dξL_diff
                        contrib_R += dJ_dξR_diff
                        if upwind_is_left
                            contrib_C += adv_coef
                        else
                            contrib_R += adv_coef
                        end
                    end
                    dS_dξ_C_buf[j] += factor_T * contrib_C
                    dS_dξ_R_buf[j] += factor_T * contrib_R
                end
            end
        end

        # === Compose T-row Jacobian ===
        S_conv = -ΣcgJavg * dTdx

        dS_dT_L += -ΣcgJavg * ddTdx_dT_L
        dS_dT_C += -ΣcgJavg * ddTdx_dT_C - Σdcg_dT_Javg * dTdx
        dS_dT_R += -ΣcgJavg * ddTdx_dT_R

        if has_left
            fast_sparse_add!(A, smap, i, iL, dS_dT_L * inv_ρcp)
        end
        fast_sparse_add!(A, smap, i, i, dS_dT_C * inv_ρcp - S_conv * dρcp_dT * inv_ρcp_sq)
        if has_right
            fast_sparse_add!(A, smap, i, iR, dS_dT_R * inv_ρcp)
        end
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            if has_left
                fast_sparse_add!(A, smap, i, ξ_off + iL, dS_dξ_L_buf[j] * inv_ρcp)
            end
            fast_sparse_add!(
                A,
                smap,
                i,
                ξ_off + i,
                dS_dξ_C_buf[j] * inv_ρcp - S_conv * dρcp_dξ[j] * inv_ρcp_sq,
            )
            if has_right
                fast_sparse_add!(A, smap, i, ξ_off + iR, dS_dξ_R_buf[j] * inv_ρcp)
            end
        end
    end
    return nothing
end

# Public dispatcher — the analytical override is the default. The AD path is
# available as `_contribute_to_jacobian_ad!` for the equivalence test in
# `test/jacobian/test_operator_darcy_transport_analytical.jl`.

function contribute_to_jacobian!(
    jac::StructuredJacobian,
    op::DarcyTransportOperator,
    u::AbstractVector,
    t::Real,
    ws::DarcyTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}
    return _contribute_to_jacobian_analytical!(
        jac,
        op,
        u,
        t,
        ws,
        state,
        mesh,
        material,
        layout,
    )
end

# ============================================================================
# AD-default path — preserved as `_contribute_to_jacobian_ad!` for the
# rtol=1e-12 equivalence reference.
# ============================================================================

"""
	contribute_to_jacobian!(jac, op::DarcyTransportOperator,
	                        u, t, ws, state, mesh, material, layout) -> nothing

Per-cell AD path. Each cell `i`:

1. Pack the 3-cell `(T, ξ)` stencil into `ws.x_buf`.
2. Run `ForwardDiff.jacobian!` over `apply_cell!` to produce a
   `(1 + NC) × (3 · (1 + NC))` block.
3. Scatter the T-row and gas-component ξ-rows into `jac.A`. Solid-component
   ξ-rows are identically zero.
"""
function _contribute_to_jacobian_ad!(
    jac::StructuredJacobian,
    ::DarcyTransportOperator,
    u::AbstractVector,
    t::Real,
    ws::DarcyTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap

    n_live = mesh.n_active_cells
    @inbounds for i = 1:n_live
        geom = _dt_geometry(state, mesh, i)

        f! = let geom = geom, material = material, NC = NC
            (out, x_inner) -> begin
                T_L, T_C, T_R, ξ_L, ξ_C, ξ_R = _dt_unpack_stencil(x_inner, NC)
                apply_cell!(
                    out,
                    DarcyTransportOperator(),
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

        _dt_pack_stencil!(ws.x_buf, u, n_cells, NC, i)
        ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

        iL = i == 1 ? 0 : i - 1
        iR = i == n_live ? 0 : i + 1

        function _scatter_row(block_row::Int, gr::Int)
            if iL != 0
                fast_sparse_add!(A, smap, gr, iL, ws.J_local[block_row, 1])
            end
            fast_sparse_add!(A, smap, gr, i, ws.J_local[block_row, 2])
            if iR != 0
                fast_sparse_add!(A, smap, gr, iR, ws.J_local[block_row, 3])
            end
            for j = 1:NC
                ξ_off = n_cells + (j - 1) * n_cells
                if iL != 0
                    fast_sparse_add!(A, smap, gr, ξ_off + iL, ws.J_local[block_row, 3+j])
                end
                fast_sparse_add!(A, smap, gr, ξ_off + i, ws.J_local[block_row, 3+NC+j])
                if iR != 0
                    fast_sparse_add!(
                        A,
                        smap,
                        gr,
                        ξ_off + iR,
                        ws.J_local[block_row, 3+2*NC+j],
                    )
                end
            end
        end

        _scatter_row(1, i)
        for j = 1:NC
            _scatter_row(1 + j, n_cells + (j - 1) * n_cells + i)
        end
    end
    return nothing
end
