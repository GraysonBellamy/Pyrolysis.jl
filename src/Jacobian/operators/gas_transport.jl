"""
The `GasTransportOperator` — the Fickian gas-transport contribution extracted from the residual.

Represents the Fickian gas-transport contribution to the residual:

```
dT[i]    += S_conv(T, ξ) / ρcp_eff(T_i, ξ_i)
dξ[j, i] += -div_J_gas[j, i]            (only for j with `is_gas == true`)
```

where `J_gas[j, f] = gas_flux(λ_face, ξ_L[j], ξ_R[j], T_L, T_R, d_LR)` per
gas component and `S_conv = −Σ_{g} cp_g(T_i) · J̄_g,i · ∂T/∂x|_i`. The
optional generation-sink source `S_gen = −Σ_g Δh_g(T_C, T_ref) · (∇·J_g)`
(with `Δh_g = cp_g(T_m)·(T_C − T_ref)`, `T_m = (T_C+T_ref)/2`; the true ∫cp
midpoint rule, matching `Physics.compute_gas_generation_sink!`) is included via
the `T_ref` keyword: the sentinel `T_ref = NaN` (default) disables the term, and
a finite `T_ref` enables it (resolved upstream from the BC set in the
`Structured` backend when `energy_form === :with_generation_sink`).

# Stencil

3-cell tridiagonal: `(cells = -1:1,)`. Boundary handling mirrors
[`ConductionOperator`](@ref): `has_left` / `has_right` gates drop the missing-
neighbor flux contribution. At boundary cells this operator contributes **no
S_conv** — the whole boundary-cell S_conv energy row (interior + BC face flux
× the T_s-aware blended gradient) is owned by [`BCConvectionOperator`](@ref).
The boundary-face flux contributions to `dξ_gas` and to the energy row's
`S_gen` slice are the [`BCMassOperator`](@ref)'s responsibility.

# Math correspondence with the existing residual

`compute_rhs!` computes `J_gas` via `compute_gas_fluxes!`, overwrites the
boundary entries via `apply_mass_bc_flux!`, then takes the species
divergence into `dξ` and computes `S_conv` from cell-center T gradients +
face-averaged J_gas. `GasTransportOperator.apply!` reproduces the *interior-
only* slice (boundary faces zeroed) plus, at interior cells, the
`S_conv / ρcp_eff` energy contribution. Together with `BCMassOperator` and
`BCConvectionOperator`, the three cover the same residual contribution as
the existing path.

# Simplifications

- No species-mass scaling (`use_species_scaling = false`).
- Interior face fluxes only (boundary handled by the BC operators).
"""

using ForwardDiff
using ..Materials: Material, is_gas, delta_enthalpy
using ..Properties: effective_gas_transfer, effective_heat_capacity, harmonic_mean_weighted
using ..Physics: gas_flux

"""
	GasTransportOperator <: Operator

Stateless tag; per-problem AD workspace lives on the backend cache.
"""
struct GasTransportOperator <: Operator end

reads(::GasTransportOperator) = ()
writes(::GasTransportOperator) = ()
stencil(::GasTransportOperator) = (cells = -1:1,)

"""
	GasTransportADWorkspace{NC, CFG}

Workspace for `GasTransportOperator.contribute_to_jacobian!`. Carries
the AD-default's stencil-AD plumbing (`x_buf`, `f_buf`, `J_local`,
`cfg`) and the analytical override's per-cell scratch arrays for the
T-row's ξ-derivative accumulators (`dS_dξ_L_buf`, `_C_buf`, `_R_buf`).
The scratch arrays replace per-face `ntuple(j -> ..., Val(NC))`
closures, which heap-allocate.

# Fields
- AD plumbing: see field docs above.
- `dS_dξ_L_buf, dS_dξ_C_buf, dS_dξ_R_buf::Vector{Float64}` — length-`NC`
  scratch for the T-row's ξ-derivative accumulators.
"""
mutable struct GasTransportADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    dS_dξ_L_buf::Vector{Float64}
    dS_dξ_C_buf::Vector{Float64}
    dS_dξ_R_buf::Vector{Float64}
end

function build_gas_transport_workspace(::Material{NC}) where {NC}
    stencil_dim = 3 * (1 + NC)
    out_dim = 1 + NC
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, out_dim)
    J_local = zeros(Float64, out_dim, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return GasTransportADWorkspace{NC,typeof(cfg)}(
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
# Identical layout to ConductionOperator's _pack_stencil!.
@inline function _gt_pack_stencil!(
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

@inline function _gt_unpack_stencil(x::AbstractVector, NC::Int)
    @inbounds T_L = x[1]
    @inbounds T_C = x[2]
    @inbounds T_R = x[3]
    ξ_L = ntuple(j -> @inbounds(x[3+j]), NC)
    ξ_C = ntuple(j -> @inbounds(x[3+NC+j]), NC)
    ξ_R = ntuple(j -> @inbounds(x[3+2*NC+j]), NC)
    return T_L, T_C, T_R, ξ_L, ξ_C, ξ_R
end

# Compute the gas-flux row for one face (between L and C, say).
# Returns NTuple{NC, Tu} of fluxes; non-gas components are 0.
@generated function _gt_face_fluxes(
    T_L::Real,
    T_C::Real,
    ξ_L::NTuple{NC,Tu},
    ξ_C::NTuple{NC,Tu},
    λ_face::Real,
    d_LR::Real,
    material::Material{NC},
) where {NC,Tu}
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                gas_flux(λ_face, ξ_L[$j], ξ_C[$j], T_L, T_C, d_LR)
            else
                zero(Tu)
            end
        end)
    end
    return :(NTuple{NC,Tu}(($(exprs...),)))
end

"""
	apply_cell!(out, op::GasTransportOperator,
	            T_L, T_C, T_R, ξ_L, ξ_C, ξ_R, geom, has_left, has_right,
	            material) -> nothing

Per-cell gas-transport contribution. Writes:
- `out[1]`     : `S_conv(T, ξ) / ρcp_eff(T_C, ξ_C)`.
- `out[1 + j]` : `-div_J_gas[j]` for gas component `j`; `0` otherwise.

`geom` is a NamedTuple `(d_L_face, d_R_face, d_LC, d_CR, A_left, A_right, V_C)`
of Float64 constants. `has_left` / `has_right` gate the boundary-face flux
(missing neighbor → zero face flux from this operator).
"""
function apply_cell!(
    out::AbstractVector,
    ::GasTransportOperator,
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
    T_ref::Real = NaN,
) where {NC,Tu}

    # Cell-local properties at C.
    λ_C = effective_gas_transfer(ξ_C, T_C, material)
    ρcp = effective_heat_capacity(ξ_C, T_C, material)

    # Left face flux per component (NTuple{NC, Tu}; non-gas components zero).
    J_left = if has_left
        λ_L = effective_gas_transfer(ξ_L, T_L, material)
        λ_face_L = harmonic_mean_weighted(λ_L, λ_C, d_L_face, d_LC - d_L_face)
        _gt_face_fluxes(T_L, T_C, ξ_L, ξ_C, λ_face_L, d_LC, material)
    else
        ntuple(_ -> zero(Tu), Val(NC))
    end

    # Right face flux per component.
    J_right = if has_right
        λ_R = effective_gas_transfer(ξ_R, T_R, material)
        λ_face_R = harmonic_mean_weighted(λ_C, λ_R, d_R_face, d_CR - d_R_face)
        _gt_face_fluxes(T_C, T_R, ξ_C, ξ_R, λ_face_R, d_CR, material)
    else
        ntuple(_ -> zero(Tu), Val(NC))
    end

    # dT/dx at cell center: opposite-distance-weighted slopes (second order
    # on nonuniform meshes), mirroring the interior branch of
    # `Physics.compute_gas_convective_source!`. At boundary cells the WHOLE
    # S_conv energy row — boundary-face flux × T_s-aware gradient — is owned
    # by `BCConvectionOperator`, so this operator contributes no S_conv there
    # (dTdx = 0 zeroes it; the ξ-rows and the linear S_gen slice stay here).
    dTdx = if has_left && has_right
        opposite_distance_gradient(T_L, T_C, T_R, d_LC, d_CR)
    else
        zero(T_C)
    end

    # S_conv = -Σ_{g: is_gas} cp_g(T_C) · J_avg[g] · dTdx
    S_conv = _gt_S_conv_unrolled(T_C, dTdx, J_left, J_right, material)

    # S_gen = -Σ_g Δh_g(T_C, T_ref) · div_J[g]
    #         (only when energy_form === :with_generation_sink; sentinel
    #          T_ref = NaN disables the term).
    S_gen = if isnan(T_ref)
        zero(promote_type(typeof(T_C), Tu))
    else
        _gt_S_gen_unrolled(T_C, J_left, J_right, A_left, A_right, V_C, T_ref, material)
    end

    # Output assembly
    @inbounds out[1] = (S_conv + S_gen) / ρcp
    # Species rows: -div_J for gases, 0 for solids.
    _gt_species_div_unrolled!(out, J_left, J_right, A_left, A_right, V_C, material)
    return nothing
end

# S_gen,C = -Σ_{g: is_gas} Δh_g(T_C, T_ref) · div_J[g]
# (matches Physics.compute_gas_generation_sink! — true ∫cp via delta_enthalpy)
@generated function _gt_S_gen_unrolled(
    T_C::Real,
    J_left::NTuple{NC,Tu},
    J_right::NTuple{NC,Tu},
    A_left::Real,
    A_right::Real,
    V_C::Real,
    T_ref::Real,
    material::Material{NC},
) where {NC,Tu}
    exprs = Expr[]
    push!(exprs, :(s = zero(promote_type(typeof(T_C), Tu))))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                div_J = (J_right[$j] * A_right - J_left[$j] * A_left) / V_C
                cmp = material.components[$j]
                s += delta_enthalpy(cmp, T_C, T_ref) * div_J
            end
        end)
    end
    push!(exprs, :(return -s))
    return Expr(:block, exprs...)
end

@generated function _gt_S_conv_unrolled(
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

@generated function _gt_species_div_unrolled!(
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

# Geometry slice for cell i — Float64 constants captured in the AD closure.
@inline function _gt_geometry(state, mesh, i::Int)
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
	apply!(dT, dξ, op::GasTransportOperator, state, mesh, material) -> nothing

Operator-level loop accumulating per-cell `(S_conv / ρcp_eff, -div_J_gas[:])`
into `(dT, dξ)`.
"""
function apply!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    ::GasTransportOperator,
    state,
    mesh,
    material::Material{NC},
) where {NC}
    out = zeros(eltype(dT), 1 + NC)
    n = length(state.T)
    @inbounds for i = 1:n
        geom = _gt_geometry(state, mesh, i)
        iL = i == 1 ? i : i - 1
        iR = i == n ? i : i + 1
        Tu = eltype(state.T)
        ξ_L = ntuple(j -> state.ξ[j, iL], Val(NC))
        ξ_C = ntuple(j -> state.ξ[j, i], Val(NC))
        ξ_R = ntuple(j -> state.ξ[j, iR], Val(NC))
        apply_cell!(
            out,
            GasTransportOperator(),
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
# Analytical override (GasTransport)
#
# Same Jacobian as the AD path below but with property derivatives extracted
# from one Dual-mode pass per cell (input dim = NC+1) — replacing the AD
# default's 3·(NC+1)-dim Dual chunk over the whole `apply_cell!`. The face
# fluxes' chain rule is hand-derived: the direct (T, ξ) terms come from the
# closed-form `gas_flux` formula, and the indirect path through
# `λ_face = harmonic_mean(λ_L, λ_C, …)` composes cell-`λ` derivatives via
# the harmonic-mean formula (lifted as `_harmonic_mean_with_derivs` in
# `operators/conduction.jl`; reused here).
# ============================================================================

struct _GTPropTag end

@generated function _gt_seed_inputs(T::Float64, ξ::NTuple{NC,Float64}) where {NC}
    N = NC + 1
    T_p = Expr(:tuple, [i == 1 ? :(1.0) : :(0.0) for i = 1:N]...)
    ξ_exprs = Expr[]
    for j = 1:NC
        jp = j + 1
        p = Expr(:tuple, [i == jp ? :(1.0) : :(0.0) for i = 1:N]...)
        push!(
            ξ_exprs,
            quote
                ForwardDiff.Dual{_GTPropTag,Float64,$N}(
                    ξ[$j],
                    ForwardDiff.Partials{$N,Float64}($p),
                )
            end,
        )
    end
    return quote
        T_d = ForwardDiff.Dual{_GTPropTag,Float64,$N}(
            T,
            ForwardDiff.Partials{$N,Float64}($T_p),
        )
        ξ_d = $(Expr(:tuple, ξ_exprs...))
        return T_d, ξ_d
    end
end

@generated function _gt_extract(
    d::ForwardDiff.Dual{_GTPropTag,Float64,N},
    ::Val{NC},
) where {N,NC}
    N == NC + 1 || error("_gt_extract: N=$N must equal NC+1=$(NC+1)")
    ξ_t = Expr(:tuple, [:(p[$(j+1)]) for j = 1:NC]...)
    return quote
        v = ForwardDiff.value(d)
        p = ForwardDiff.partials(d)
        return v, p[1], $ξ_t
    end
end

# Cell-center properties: (λ, ρcp) and their (T, ξ)-derivatives in one Dual
# pass.
@inline function _gt_props_at_cell(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _gt_seed_inputs(T, ξ)
    λ_d = effective_gas_transfer(ξ_d, T_d, material)
    ρcp_d = effective_heat_capacity(ξ_d, T_d, material)
    λ, dλ_dT, dλ_dξ = _gt_extract(λ_d, Val(NC))
    ρcp, dρcp_dT, dρcp_dξ = _gt_extract(ρcp_d, Val(NC))
    return (λ, dλ_dT, dλ_dξ, ρcp, dρcp_dT, dρcp_dξ)
end

# Neighbor-cell λ only.
@inline function _gt_λ_with_derivs(
    ξ::NTuple{NC,Float64},
    T::Float64,
    material::Material{NC},
) where {NC}
    T_d, ξ_d = _gt_seed_inputs(T, ξ)
    λ_d = effective_gas_transfer(ξ_d, T_d, material)
    return _gt_extract(λ_d, Val(NC))
end

# `cp_g(T_C)` and `dcp_g/dT(T_C)` for every component, returned as parallel
# NTuples. Each component's `heat_capacity` is a closure on T; we get its
# scalar derivative via a length-1 Dual seed. `@generated` to unroll the
# per-component component-tuple traversal.
@generated function _gt_cp_per_component(T::Float64, comps::Tuple, ::Val{NC}) where {NC}
    cp_exprs = Expr[]
    dcp_exprs = Expr[]
    for j = 1:NC
        push!(cp_exprs, :(comps[$j].heat_capacity(T)))
        push!(dcp_exprs, :(ForwardDiff.derivative(comps[$j].heat_capacity, T)))
    end
    cp_t = Expr(:tuple, cp_exprs...)
    dcp_t = Expr(:tuple, dcp_exprs...)
    return quote
        cp = $cp_t
        dcp_dT = $dcp_t
        return cp, dcp_dT
    end
end

# `is_gas` per component as a compile-time NTuple{NC, Bool}.
@generated function _gt_is_gas_tuple(comps::Tuple, ::Val{NC}) where {NC}
    exprs = [:(is_gas(comps[$j])) for j = 1:NC]
    return :($(Expr(:tuple, exprs...)))
end

# Reuse the harmonic-mean derivative helper from operators/conduction.jl.
# Both files are included into the same `Pyrolysis.Jacobian` module, so
# `_harmonic_mean_with_derivs` is in scope here.

# Direct partial derivatives of `gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d)`
# w.r.t. its scalar arguments. Returns a 5-tuple
# `(∂J/∂λ, ∂J/∂ξ_L, ∂J/∂ξ_R, ∂J/∂T_L, ∂J/∂T_R)`. Math comes from
# `J = -λ · (ξ_R T_R − ξ_L T_L) / (T_face · d)` with `T_face = (T_L + T_R)/2`.
@inline function _gt_J_direct_partials(
    λ::Float64,
    ξ_L::Float64,
    ξ_R::Float64,
    T_L::Float64,
    T_R::Float64,
    d::Float64,
)
    T_face = 0.5 * (T_L + T_R)
    denom = T_face * d
    N_num = ξ_R * T_R - ξ_L * T_L
    J = -λ * N_num / denom
    dJ_dλ = N_num == 0.0 ? 0.0 : -N_num / denom
    dJ_dξ_L = λ * T_L / denom
    dJ_dξ_R = -λ * T_R / denom
    # ∂J/∂T_L = λ·ξ_L/(T_face·d) − J/(2·T_face)
    dJ_dT_L = λ * ξ_L / denom - J / (2 * T_face)
    dJ_dT_R = -λ * ξ_R / denom - J / (2 * T_face)
    return J, dJ_dλ, dJ_dξ_L, dJ_dξ_R, dJ_dT_L, dJ_dT_R
end

"""
	_contribute_to_jacobian_analytical!(jac, op::GasTransportOperator, …)

Analytical override of `contribute_to_jacobian!` for
`GasTransportOperator`. Iterates cells; for each cell evaluates
`(λ_C, ρcp, cp_g, all derivs)` once and assembles the per-face flux
chain rule plus the `S_conv/ρcp` and `div(J)` divergences by hand.

Equivalence with the AD-default `_contribute_to_jacobian_ad!` is gated
at rtol = 1e-12 in
`test/jacobian/test_operator_gas_transport_analytical.jl`.
"""
function _contribute_to_jacobian_analytical!(
    jac::StructuredJacobian,
    ::GasTransportOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::GasTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    T_ref::Real = NaN,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap
    is_gas_t = _gt_is_gas_tuple(material.components, Val(NC))
    is_gen = !isnan(T_ref)

    @inbounds for i = 1:mesh.n_active_cells
        geom = _gt_geometry(state, mesh, i)
        has_left = geom.has_left
        has_right = geom.has_right

        # Single-cell mesh edge case — contributes nothing.
        (!has_left && !has_right) && continue

        iL = has_left ? (i - 1) : 0
        iR = has_right ? (i + 1) : 0

        T_C = u[i]
        ξ_C = ntuple(j -> u[n_cells+(j-1)*n_cells+i], Val(NC))
        (λ_C, dλ_C_dT, dλ_C_dξ, ρcp, dρcp_dT, dρcp_dξ) =
            _gt_props_at_cell(ξ_C, T_C, material)
        inv_ρcp = 1.0 / ρcp
        inv_ρcp_sq = inv_ρcp * inv_ρcp

        (cp_t, dcp_dT_t) = _gt_cp_per_component(T_C, material.components, Val(NC))
        # S_gen uses the midpoint-rule sensible enthalpy Δh_g = cp_g(T_m)·(T_C−T_ref)
        # with T_m = (T_C+T_ref)/2 (matching `delta_enthalpy`), so its cp factors are
        # evaluated at T_m — distinct from S_conv, which uses cp_g(T_C). The explicit
        # ∂Δh_g/∂T_C cp-chain is ½·cp_g′(T_m). When not in generation-sink mode the
        # midpoint tuples are unused (T_ref = NaN sentinel).
        (cpm_t, dcpm_dT_t) =
            is_gen ?
            _gt_cp_per_component(0.5 * (T_C + T_ref), material.components, Val(NC)) :
            (cp_t, dcp_dT_t)

        # dTdx and its (T_L, T_C, T_R) coefficients — function of geometry only.
        # Matches `apply_cell!`: interior cells use the opposite-distance
        # gradient; boundary cells have dTdx ≡ 0 here (their S_conv energy row
        # is owned whole by `BCConvectionOperator`), which zeroes every S_conv
        # accumulator below while the ξ-rows and S_gen pieces still assemble.
        if has_left && has_right
            T_L = u[iL]
            T_R = u[iR]
            dTdx = opposite_distance_gradient(T_L, T_C, T_R, geom.d_LC, geom.d_CR)
            (ddTdx_dT_L, ddTdx_dT_C, ddTdx_dT_R) =
                opposite_distance_gradient_coeffs(geom.d_LC, geom.d_CR)
        elseif has_right
            T_R = u[iR]
            T_L = T_C   # not used; placeholder
            dTdx = 0.0
            ddTdx_dT_L = 0.0
            ddTdx_dT_C = 0.0
            ddTdx_dT_R = 0.0
        else  # has_left only
            T_L = u[iL]
            T_R = T_C   # placeholder
            dTdx = 0.0
            ddTdx_dT_L = 0.0
            ddTdx_dT_C = 0.0
            ddTdx_dT_R = 0.0
        end

        # Σ_g cp_g · J_avg_g, used for the T-row's S_conv = -Σ · dTdx.
        # We accumulate this and its derivatives in stencil order.
        # The (∂S/∂x) decomposition has three parts:
        #   ∂S/∂x = -Σ_g [(∂cp_g/∂x) · J_avg · dTdx
        #                  + cp_g · (∂J_avg/∂x) · dTdx
        #                  + cp_g · J_avg · (∂dTdx/∂x)]
        # `∂cp_g/∂x` is non-zero only for x = T_C (and equals dcp_g/dT(T_C)).
        # `∂dTdx/∂x` is non-zero only for x ∈ {T_L, T_C, T_R} (geom coeffs).

        # Accumulators for ∂S/∂x and div(J)·... contributions.
        # T-row accumulators (entries that go into row `i` of A).
        dS_dT_L = 0.0
        dS_dT_C = 0.0
        dS_dT_R = 0.0
        # Reset workspace ξ-deriv accumulators in place.
        dS_dξ_L_buf = ws.dS_dξ_L_buf
        dS_dξ_C_buf = ws.dS_dξ_C_buf
        dS_dξ_R_buf = ws.dS_dξ_R_buf
        for j = 1:NC
            dS_dξ_L_buf[j] = 0.0
            dS_dξ_C_buf[j] = 0.0
            dS_dξ_R_buf[j] = 0.0
        end
        # Σ_g cp_g · J_avg_g — for the dTdx and ρcp chain-rule terms.
        ΣcgJavg = 0.0
        Σdcg_dT_Javg = 0.0
        # S_gen accumulators (used only when is_gen).
        # ΣcgDivJ = Σ_g cp_g(T_m) · div_J_g (assembled face-by-face from
        # div_J_g = (J_R · A_R − J_L · A_L)/V_C). Used for the T_C explicit
        # `(T_C − T_ref)` derivative and the `S_gen` value (= -dT_gen·ΣcgDivJ
        # = -Σ Δh_g·div_J_g).
        # ΣdcgDT_DivJ = Σ_g ½·(dcp_g/dT)(T_m) · div_J_g — the explicit cp-chain of
        # Δh_g w.r.t. T_C (∂Δh_g/∂T_C = cp_g(T_m) + ½(T_C−T_ref)cp_g′(T_m)).
        ΣcgDivJ = 0.0
        ΣdcgDT_DivJ = 0.0
        dT_gen = is_gen ? (T_C - T_ref) : 0.0

        # ξ-row accumulators (per gas component g) — hold ∂(out[1+g])/∂x
        # for stencil x. We scatter immediately per face / component to
        # avoid keeping NC × stencil scratch arrays.

        # --- Process both faces. We compute J_g and its derivatives at each
        # face, accumulate the T-row pieces, and immediately scatter the
        # ξ-row pieces into A.

        # Closure factory replaced by inline branches to keep the code
        # stack-allocated.

        # === Left face (between iL and i, when has_left) ===
        # face role: L-cell side = neighbor (iL); R-cell side = center (i).
        if has_left
            ξ_L = ntuple(j -> u[n_cells+(j-1)*n_cells+iL], Val(NC))
            (λ_L, dλ_L_dT, dλ_L_dξ) = _gt_λ_with_derivs(ξ_L, T_L, material)
            d1 = geom.d_L_face
            d2 = geom.d_LC - geom.d_L_face
            (λ_face, dλf_dλL, dλf_dλC) = _harmonic_mean_with_derivs(λ_L, λ_C, d1, d2)

            # Per-gas processing.
            for g = 1:NC
                is_gas_t[g] || continue
                ξ_L_g = ξ_L[g]
                ξ_C_g = ξ_C[g]
                J, dJ_dλ, dJ_dξL, dJ_dξR, dJ_dTL, dJ_dTR =
                    _gt_J_direct_partials(λ_face, ξ_L_g, ξ_C_g, T_L, T_C, geom.d_LC)

                cp_g = cp_t[g]
                dcp_g_dT = dcp_dT_t[g]

                # Compose face-flux derivatives w.r.t. stencil vars.
                #   ∂J/∂T_L_full = dJ_dTL + dJ_dλ · dλf_dλL · dλ_L_dT
                #   ∂J/∂T_C_full = dJ_dTR + dJ_dλ · dλf_dλC · dλ_C_dT
                #   ∂J/∂ξ_L_g_direct = dJ_dξL                (g-only direct)
                #   ∂J/∂ξ_L_j_via_λ  = dJ_dλ · dλf_dλL · dλ_L_dξ[j]
                #   ∂J/∂ξ_C_g_direct = dJ_dξR
                #   ∂J/∂ξ_C_j_via_λ  = dJ_dλ · dλf_dλC · dλ_C_dξ[j]
                dJ_dT_L_full = dJ_dTL + dJ_dλ * dλf_dλL * dλ_L_dT
                dJ_dT_C_full = dJ_dTR + dJ_dλ * dλf_dλC * dλ_C_dT
                # Note: T_R doesn't enter J_left; ξ_R doesn't enter J_left.

                # === ξ-row scatter for component g (the AD-default writes
                # `out[1+g] = -(J_R·A_R - J_L·A_L)/V_C`; here we add the
                # left-face contribution: `+J_L · A_L / V_C`).
                row_g = n_cells + (g - 1) * n_cells + i
                coef_AV_L = geom.A_left / geom.V_C
                # T_L
                fast_sparse_add!(A, smap, row_g, iL, dJ_dT_L_full * coef_AV_L)
                # T_C
                fast_sparse_add!(A, smap, row_g, i, dJ_dT_C_full * coef_AV_L)
                # ξ columns (j from 1:NC). Direct path only on j==g; via-λ on all j.
                for j = 1:NC
                    ξ_off = n_cells + (j - 1) * n_cells
                    # left ξ_L_j column
                    contrib_L = dJ_dλ * dλf_dλL * dλ_L_dξ[j]
                    if j == g
                        contrib_L += dJ_dξL
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + iL, contrib_L * coef_AV_L)
                    # center ξ_C_j column
                    contrib_C = dJ_dλ * dλf_dλC * dλ_C_dξ[j]
                    if j == g
                        contrib_C += dJ_dξR
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + i, contrib_C * coef_AV_L)
                end

                # === T-row contribution from this face ===
                # J_avg_g picks up J_left/2; dJ/dx pieces multiply 1/2 too.
                ΣcgJavg += cp_g * (J / 2)
                Σdcg_dT_Javg += dcp_g_dT * (J / 2)

                # S_gen contribution: L-face contributes -J_L · A_L / V_C to
                # div_J_g. Stencil-derivative factor for this face (S_conv +
                # S_gen combined): cp_g · (-0.5·dTdx + dT_gen·coef_AV_L).
                factor_T = -cp_g * 0.5 * dTdx
                if is_gen
                    # S_gen uses Δh_g = cp_g(T_m)·(T_C−T_ref): cp at T_m here,
                    # and the explicit cp-chain factor is ½·cp_g′(T_m).
                    cp_g_m = cpm_t[g]
                    dcp_g_m = dcpm_dT_t[g]
                    factor_T += dT_gen * cp_g_m * coef_AV_L
                    ΣcgDivJ += cp_g_m * (-J * coef_AV_L)
                    ΣdcgDT_DivJ += 0.5 * dcp_g_m * (-J * coef_AV_L)
                end

                # ∂S/∂T_L and ∂S/∂T_C contributions (no cp_g(T_L)/cp_g(T_R)
                # chain — cp_g depends only on T_C).
                dS_dT_L += factor_T * dJ_dT_L_full
                dS_dT_C += factor_T * dJ_dT_C_full

                # ∂S/∂ξ_L_j and ∂S/∂ξ_C_j contributions for this face.
                for j = 1:NC
                    contrib_L = dJ_dλ * dλf_dλL * dλ_L_dξ[j]
                    if j == g
                        contrib_L += dJ_dξL
                    end
                    dS_dξ_L_buf[j] += factor_T * contrib_L

                    contrib_C = dJ_dλ * dλf_dλC * dλ_C_dξ[j]
                    if j == g
                        contrib_C += dJ_dξR
                    end
                    dS_dξ_C_buf[j] += factor_T * contrib_C
                end
            end
        end

        # === Right face (between i and iR, when has_right) ===
        if has_right
            ξ_R = ntuple(j -> u[n_cells+(j-1)*n_cells+iR], Val(NC))
            (λ_R, dλ_R_dT, dλ_R_dξ) = _gt_λ_with_derivs(ξ_R, T_R, material)
            d1 = geom.d_R_face
            d2 = geom.d_CR - geom.d_R_face
            (λ_face, dλf_dλL, dλf_dλC) = _harmonic_mean_with_derivs(λ_C, λ_R, d1, d2)
            # face role: L-cell side = center (i); R-cell side = neighbor (iR).
            # So `dλf_dλL` corresponds to ∂λ_face/∂λ_C, `dλf_dλC` to ∂λ_face/∂λ_R.

            for g = 1:NC
                is_gas_t[g] || continue
                ξ_C_g = ξ_C[g]
                ξ_R_g = ξ_R[g]
                J, dJ_dλ, dJ_dξL, dJ_dξR, dJ_dTL, dJ_dTR =
                    _gt_J_direct_partials(λ_face, ξ_C_g, ξ_R_g, T_C, T_R, geom.d_CR)

                cp_g = cp_t[g]
                dcp_g_dT = dcp_dT_t[g]

                # Compose face-flux derivatives w.r.t. stencil vars.
                dJ_dT_C_full = dJ_dTL + dJ_dλ * dλf_dλL * dλ_C_dT
                dJ_dT_R_full = dJ_dTR + dJ_dλ * dλf_dλC * dλ_R_dT
                # T_L doesn't enter J_right.

                row_g = n_cells + (g - 1) * n_cells + i
                coef_AV_R = geom.A_right / geom.V_C
                # Right-face contributes -J_R·A_R/V_C to out[1+g].
                # T_C
                fast_sparse_add!(A, smap, row_g, i, -dJ_dT_C_full * coef_AV_R)
                # T_R
                fast_sparse_add!(A, smap, row_g, iR, -dJ_dT_R_full * coef_AV_R)
                for j = 1:NC
                    ξ_off = n_cells + (j - 1) * n_cells
                    contrib_C = dJ_dλ * dλf_dλL * dλ_C_dξ[j]
                    if j == g
                        contrib_C += dJ_dξL
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + i, -contrib_C * coef_AV_R)
                    contrib_R = dJ_dλ * dλf_dλC * dλ_R_dξ[j]
                    if j == g
                        contrib_R += dJ_dξR
                    end
                    fast_sparse_add!(A, smap, row_g, ξ_off + iR, -contrib_R * coef_AV_R)
                end

                ΣcgJavg += cp_g * (J / 2)
                Σdcg_dT_Javg += dcp_g_dT * (J / 2)

                # S_gen contribution: R-face contributes +J_R · A_R / V_C to
                # div_J_g. Stencil-derivative factor for this face (S_conv +
                # S_gen combined): cp_g · (-0.5·dTdx − dT_gen·coef_AV_R).
                factor_T = -cp_g * 0.5 * dTdx
                if is_gen
                    # S_gen uses Δh_g = cp_g(T_m)·(T_C−T_ref): cp at T_m here,
                    # and the explicit cp-chain factor is ½·cp_g′(T_m).
                    cp_g_m = cpm_t[g]
                    dcp_g_m = dcpm_dT_t[g]
                    factor_T -= dT_gen * cp_g_m * coef_AV_R
                    ΣcgDivJ += cp_g_m * (J * coef_AV_R)
                    ΣdcgDT_DivJ += 0.5 * dcp_g_m * (J * coef_AV_R)
                end

                dS_dT_C += factor_T * dJ_dT_C_full
                dS_dT_R += factor_T * dJ_dT_R_full

                for j = 1:NC
                    contrib_C = dJ_dλ * dλf_dλL * dλ_C_dξ[j]
                    if j == g
                        contrib_C += dJ_dξL
                    end
                    dS_dξ_C_buf[j] += factor_T * contrib_C

                    contrib_R = dJ_dλ * dλf_dλC * dλ_R_dξ[j]
                    if j == g
                        contrib_R += dJ_dξR
                    end
                    dS_dξ_R_buf[j] += factor_T * contrib_R
                end
            end
        end

        # === Compose T-row Jacobian ===
        # S_conv = -ΣcgJavg · dTdx; S_gen = -dT_gen · ΣcgDivJ.
        S_conv = -ΣcgJavg * dTdx
        S_gen = is_gen ? -dT_gen * ΣcgDivJ : 0.0
        S_total = S_conv + S_gen
        out_T = S_total * inv_ρcp

        # ∂S/∂x = (already accumulated in dS_dT_*, dS_dξ_*_buf for the J-chain
        #         and the dTdx-coefficient pieces); plus:
        #  - the cp_g(T_C) chain rule: -Σ_g (dcp_g/dT)·J_avg·dTdx contributes
        #    only to ∂S/∂T_C.
        #  - the dTdx chain rule via geometry coefficients on T_L/T_C/T_R:
        #    -ΣcgJavg · ∂dTdx/∂x_T.

        dS_dT_L += -ΣcgJavg * ddTdx_dT_L
        dS_dT_C += -ΣcgJavg * ddTdx_dT_C - Σdcg_dT_Javg * dTdx
        dS_dT_R += -ΣcgJavg * ddTdx_dT_R

        # T_C explicit S_gen terms: ∂Δh_g/∂T_C = cp_g(T_m) + ½(T_C−T_ref)cp_g′(T_m).
        # First piece via ΣcgDivJ (= Σ cp_g(T_m)·div_J), second via ΣdcgDT_DivJ
        # (= Σ ½·cp_g′(T_m)·div_J), scaled by dT_gen = (T_C−T_ref).
        if is_gen
            dS_dT_C += -ΣcgDivJ - dT_gen * ΣdcgDT_DivJ
        end

        # ∂out[1]/∂x = ∂S/∂x · inv_ρcp − S_total · ∂ρcp/∂x · inv_ρcp²
        # Only T_C and ξ_C carry the −S·∂ρcp correction.
        if has_left
            fast_sparse_add!(A, smap, i, iL, dS_dT_L * inv_ρcp)
        end
        fast_sparse_add!(A, smap, i, i, dS_dT_C * inv_ρcp - S_total * dρcp_dT * inv_ρcp_sq)
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
                dS_dξ_C_buf[j] * inv_ρcp - S_total * dρcp_dξ[j] * inv_ρcp_sq,
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
# `test/jacobian/test_operator_gas_transport_analytical.jl`.

function contribute_to_jacobian!(
    jac::StructuredJacobian,
    op::GasTransportOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::GasTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    T_ref::Real = NaN,
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
        T_ref,
    )
end

# ============================================================================
# AD-default path — preserved as `_contribute_to_jacobian_ad!` for the
# rtol=1e-12 equivalence reference.
# ============================================================================

"""
	contribute_to_jacobian!(jac, op::GasTransportOperator,
	                        u, ctx, t, ws, state, mesh, material, layout) -> nothing

Per-cell AD path. Each cell `i`:

1. Pack the 3-cell `(T, ξ)` stencil into `ws.x_buf`.
2. Run `ForwardDiff.jacobian!` over `apply_cell!` to produce a
   `(1 + NC) × (3 · (1 + NC))` block.
3. Scatter the T-row and gas-component ξ-rows into `jac.A`. Solid-component
   ξ-rows are identically zero (the closure outputs zero for them) — those
   scatter-adds are no-ops.
"""
function _contribute_to_jacobian_ad!(
    jac::StructuredJacobian,
    ::GasTransportOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::GasTransportADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    T_ref::Real = NaN,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap
    out_dim = 1 + NC

    n_live = mesh.n_active_cells
    @inbounds for i = 1:n_live
        geom = _gt_geometry(state, mesh, i)

        f! = let geom = geom, material = material, NC = NC, T_ref = T_ref
            (out, x_inner) -> begin
                T_L, T_C, T_R, ξ_L, ξ_C, ξ_R = _gt_unpack_stencil(x_inner, NC)
                apply_cell!(
                    out,
                    GasTransportOperator(),
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
                    T_ref,
                )
            end
        end

        _gt_pack_stencil!(ws.x_buf, u, n_cells, NC, i)
        ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

        # Scatter. Output rows: 1 → T_i; 1+j → ξ_{j, i}.
        # Input cols (stencil): 1,2,3 → T_{iL,i,iR}; 3+(k-1)*NC+j → ξ_{j, i_k}.
        iL = i == 1 ? 0 : i - 1
        iR = i == n_live ? 0 : i + 1

        # Helper: write block_row's contributions into global row `gr`.
        # Skips boundary self-fill columns.
        function _scatter_row(block_row::Int, gr::Int)
            # T columns
            if iL != 0
                fast_sparse_add!(A, smap, gr, iL, ws.J_local[block_row, 1])
            end
            fast_sparse_add!(A, smap, gr, i, ws.J_local[block_row, 2])
            if iR != 0
                fast_sparse_add!(A, smap, gr, iR, ws.J_local[block_row, 3])
            end
            # ξ columns
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

        # T row
        _scatter_row(1, i)
        # ξ rows (gas + solid; solid rows are zero so adds are no-ops in value
        # but still O(1) operations — acceptable).
        for j = 1:NC
            _scatter_row(1 + j, n_cells + (j - 1) * n_cells + i)
        end
    end
    return nothing
end
