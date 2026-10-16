"""
The `ChiOperator` — pyrolysis-progress (χ) row Jacobian.

Owns the χ-row residual and its Jacobian for the variable cross-section
(WithChi) family of materials. The residual reads:

```
dχᵢ/dt = S_dry,i · Vᵢ / m_init,i
```

where `S_dry,i = Σⱼ ∈ dry-pyrolysis-gas Sⱼ,i(T_i, ξ_·,i)` and `Vᵢ = A · Δzᵢ`.
For identity-law materials (no `lateral_shrinkage_law`), `dχᵢ/dt ≡ 0`; the
operator skips both `apply!` and the Jacobian writes.

# Jacobian decomposition

- Cell-local in `(T_i, ξ_·,i)`: captured by ForwardDiff over
  `(T, ξ_·) ↦ S_dry(T, ξ_·) · V_i / m_init,i` with `V_i` and `m_init,i`
  fixed Float64 closures.
- z columns: `V_i = A · (z_{i+1} - z_i)` makes `dχᵢ/dt` depend on `z_i` and
  `z_{i+1}`; the [`GeometryOperator`](@ref)'s AD pass over the z subvector
  writes those entries (exact `:full` profile only).
- χ-column rank-1: `∂(dχᵢ/dt)/∂χⱼ = (S_dry,i · L_i / m_init,i) · ∂A/∂χⱼ`,
  factored as `u[i] · v[j]`. Pushed onto `jac.structured` as an
  [`EmbeddedRank1Coupling`](@ref) at `(χ_off, χ_off, n_global)`.

`ws.dA_dχ` must have been populated by
[`compute_theta_A_derivatives!`](@ref) on the upstream
[`ALEMeshVelocityOperator`](@ref) workspace before this operator runs.
"""

using ForwardDiff
using ..Materials: Material, n_components, n_reactions, dry_pyrolysis_gas_indices
using ..Physics: evaluate_reactions!, species_source_terms!

"""
	ChiOperator <: Operator

Stateless tag. Populates the χ-row block of the Jacobian for variable
cross-section (WithChi) materials.
"""
struct ChiOperator <: Operator end

"""
	ChiADWorkspace{Tu, NC, NR, CFG}

Per-problem AD workspace for [`ChiOperator`](@ref). Cell-local with chunk
size `1 + NC` (T + ξ_·).

# Fields
- `x_buf::Vector{Float64}` — `(1 + NC)`-length input buffer.
- `f_buf::Vector{Float64}` — length-1 scalar output buffer.
- `J_local::Matrix{Float64}` — `1 × (1 + NC)` per-cell Jacobian row.
- `cfg::CFG` — `ForwardDiff.JacobianConfig` for the closure.
- `rates`, `heats`, `dξ_dt::Vector{Tu}` — Dual-typed inner scratch.
"""
mutable struct ChiADWorkspace{Tu,NC,NR,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
    rates::Vector{Tu}
    heats::Vector{Tu}
    dξ_dt::Vector{Tu}
end

"""
	build_chi_workspace(material::Material) -> ChiADWorkspace
"""
function build_chi_workspace(material::Material{NC,NR}) where {NC,NR}
    dim = 1 + NC
    x_buf = zeros(Float64, dim)
    f_buf = zeros(Float64, 1)
    J_local = zeros(Float64, 1, dim)
    cfg = ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{dim}())
    Tu = eltype(cfg.duals[1])
    return ChiADWorkspace{Tu,NC,NR,typeof(cfg)}(
        x_buf,
        f_buf,
        J_local,
        cfg,
        zeros(Tu, NR),
        zeros(Tu, NR),
        zeros(Tu, NC),
    )
end

"""
	apply!(du, op::ChiOperator, state, mesh, material, layout) -> nothing

Operator-level χ-row residual contribution. Adds
`du[χ_off + i] += S_dry · V_i / m_init` for each cell `i`. No-op for
identity-law materials.
"""
function apply!(
    du::AbstractVector,
    ::ChiOperator,
    state,
    mesh,
    material::Material{NC,NR},
    layout,
) where {NC,NR}
    material.lateral_shrinkage_law === nothing && return nothing
    idx_dry = dry_pyrolysis_gas_indices(material)
    χ_off = chi_offset(layout)
    n_cells = layout.n_cells
    m_init = mesh.m_dry_solid_initial
    rates = zeros(NR)
    heats = zeros(NR)
    dξ = zeros(NC)
    @inbounds for i = 1:n_cells
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        fill!(rates, 0.0)
        fill!(heats, 0.0)
        evaluate_reactions!(rates, heats, ξ_i, T_i, material)
        species_source_terms!(dξ, rates, material)
        S_dry = 0.0
        for j in idx_dry
            S_dry += dξ[j]
        end
        denom = max(m_init[i], eps())
        du[χ_off+i] += S_dry * state.volume[i] / denom
    end
    return nothing
end

# Cell-local closure: (T, ξ_·) ↦ S_dry(T, ξ_·) · V_i / m_init.
function _chi_residual_kernel!(
    out::AbstractVector,
    x::AbstractVector,
    V_over_m::Float64,
    rates::AbstractVector,
    heats::AbstractVector,
    dξ_dt::AbstractVector,
    material::Material{NC,NR},
    idx_dry,
) where {NC,NR}
    T = x[1]
    ξ = ntuple(j -> x[1+j], Val(NC))
    fill!(rates, zero(eltype(rates)))
    fill!(heats, zero(eltype(heats)))
    evaluate_reactions!(rates, heats, ξ, T, material)
    species_source_terms!(dξ_dt, rates, material)
    S_dry = zero(eltype(rates))
    for j in idx_dry
        S_dry += dξ_dt[j]
    end
    out[1] = S_dry * V_over_m
    return nothing
end

"""
	contribute_to_jacobian!(jac, op::ChiOperator, u, t, ws,
	                        state, mesh, material, layout, ale_mesh_vel_ws)
	    -> nothing

Populate the χ-row block of `jac.A` plus push the χ-row × χ-col rank-1
onto `jac.structured`.

Cell-local AD captures `∂(dχᵢ/dt)/∂{T_i, ξ_·,i}` in one ForwardDiff pass
per cell (chunk = 1 + NC). The z-column entries belong to the
[`GeometryOperator`](@ref).
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::ChiOperator,
    u::AbstractVector,
    t::Real,
    ws::ChiADWorkspace{Tu,NC,NR},
    state,
    mesh,
    material::Material{NC,NR},
    layout,
    ale_mesh_vel_ws;
    assembly_profile::Symbol = :full,
) where {Tu,NC,NR}

    material.lateral_shrinkage_law === nothing && return nothing

    idx_dry = dry_pyrolysis_gas_indices(material)
    A_section = Float64(state.A_section[])
    n_cells = layout.n_cells
    χ_off = chi_offset(layout)
    m_init = mesh.m_dry_solid_initial
    A_jac = jac.A
    smap = jac.smap
    is_full = assembly_profile === :full

    # u_chi drives the χ-row × χ-col rank-1; only allocate / build it under
    # :full, since the rank-1 push is gated below. Ghost slots carry zero.
    u_chi = is_full ? zeros(Float64, n_cells) : nothing

    @inbounds for i = 1:mesh.n_active_cells
        denom = max(m_init[i], eps())
        m_inv = 1.0 / denom
        V_i = state.volume[i]
        V_over_m = V_i * m_inv

        # Pack inputs: (T_i, ξ_·,i).
        ws.x_buf[1] = state.T[i]
        for j = 1:NC
            ws.x_buf[1+j] = state.ξ[j, i]
        end

        f! =
            let V_over_m = V_over_m,
                rates = ws.rates,
                heats = ws.heats,
                dξ_dt = ws.dξ_dt,
                material = material,
                idx_dry = idx_dry

                (out, x_inner) -> _chi_residual_kernel!(
                    out,
                    x_inner,
                    V_over_m,
                    rates,
                    heats,
                    dξ_dt,
                    material,
                    idx_dry,
                )
            end
        ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

        χ_row = χ_off + i
        fast_sparse_add!(A_jac, smap, χ_row, i, ws.J_local[1, 1])
        for k = 1:NC
            ξk_col = n_cells + (k - 1) * n_cells + i
            fast_sparse_add!(A_jac, smap, χ_row, ξk_col, ws.J_local[1, 1+k])
        end

        # z-column entries (ALE only) are owned by `GeometryOperator`'s AD
        # pass over the z subvector; ChiOperator does NOT write them.
        # u_chi[i] = S_dry,i · L_i / m_init,i drives the χ-row × χ-col rank-1.
        if is_full
            S_dry_f = _eval_S_dry(
                material,
                idx_dry,
                state.T[i],
                ntuple(j -> state.ξ[j, i], Val(NC)),
            )
            L_i = V_i / A_section
            u_chi[i] = S_dry_f * L_i * m_inv
        end
    end

    # Push the χ-row × χ-col rank-1 onto jac.structured. Skipped under
    # :ale_reduced — the EmbeddedRank1Coupling would be dropped by
    # `_emit_jacobian!` anyway, and the dense χ-column sparsity slots that
    # would receive its contribution are not emitted under that profile.
    if is_full
        # Aliases the provider's buffer — the coupling only reads it, and it
        # stays valid until the next assembly rebuilds both.
        v_chi = ale_mesh_vel_ws.dA_dχ
        push!(
            jac.structured,
            EmbeddedRank1Coupling(u_chi, v_chi, χ_off, χ_off, Int(size(A_jac, 1))),
        )
    end

    return nothing
end

# Float64 evaluator for S_dry,i — separate from the Dual-typed AD path.
function _eval_S_dry(
    material::Material{NC,NR},
    idx_dry,
    T::Real,
    ξ::NTuple{NC,<:Real},
) where {NC,NR}
    rates = zeros(Float64, NR)
    heats = zeros(Float64, NR)
    dξ = zeros(Float64, NC)
    evaluate_reactions!(rates, heats, ξ, T, material)
    species_source_terms!(dξ, rates, material)
    S = 0.0
    @inbounds for j in idx_dry
        S += dξ[j]
    end
    return S
end
