"""
The `ALEMeshVelocityOperator`.

Provides the ALE cumulative mesh-velocity field
`w_node(k) = ∫_0^{z_k} θ(z') dz' = Σ_{i<k} (θ_pyrolysis,i − θ_A) · Δzᵢ`,
the per-cell strain rate `θ`, and the column-global `θ_A` (zero in NoChi
mode), together with their cell-local derivatives.

The operator's `apply!` populates these fields by wrapping the existing
Dual-aware physics functions ([`_compute_ale_mesh_velocity!`](@ref) and
[`_compute_theta_A`](@ref)) — no math is duplicated.

Every component rides the mesh, so the T and ξ rows carry no mesh-motion
coupling through `w`; the only mesh-motion term in the species rows is the
cell-local dilation `−ξⱼθᵢ`, which the backend assembles from this
operator's cached `∂θ/∂T`, `∂θ/∂ξ`. The cumulative velocity couples only
the z rows (`dz_k/dt = w_node[k]`), which this operator embeds itself in
`contribute_to_structured!` as [`EmbeddedNodePrefixSumCoupling`](@ref)
blocks plus, under WithChi, the θ_A rank-1 channels.

The per-cell factors `Δzᵢ · ∂θ_i/∂x_i` of those blocks live in the
workspace's [`PrefixSumCoupling`](@ref)s, refilled every assembly.
"""

using ForwardDiff
using ..Materials: Material, n_components, n_reactions
using ..Properties: PropertyCache
using ..Physics: volume_change_rate, evaluate_reactions!, species_source_terms!
using ..Problem: NoChi, WithChi
using ..Residual: _compute_ale_mesh_velocity!, _compute_theta_A

"""
	ALEMeshVelocityOperator <: Operator

Stateless tag. Every contribution is a structured coupling, pushed by
`contribute_to_structured!`: the z-row × (T, ξ) prefix-sum blocks and, under
WithChi, the θ_A rank-1 channels. The backend assembles the cell-local ξ-row
dilation `−ξⱼθᵢ` from this operator's cached strain-rate derivatives.
"""
struct ALEMeshVelocityOperator <: Operator end

"""
	ALEMeshVelocityWorkspace{Tu, NC, NR}

Per-problem workspace for [`ALEMeshVelocityOperator`](@ref). Holds the
derived-field output buffers (`w_node`, `theta`, `theta_A`) plus the
per-cell derivative buffers (`dθ_dT`, `dθ_dξ`, `dθA_dχ`) that the
structured-coupling builders consume.

# Fields
- `w_node::Vector{Float64}` — length `n_nodes`. Output of `apply!`.
- `theta::Vector{Float64}`  — length `n_cells`. Per-cell strain rate
  `θ_pyrolysis = Σⱼ γⱼ/ρⱼ · Sⱼ`.
- `theta_A::Base.RefValue{Float64}` — column-global radial channel.
  Zero in `NoChi`; populated by `_compute_theta_A` in `WithChi`.
- `dθ_dT::Vector{Float64}` — `∂θ_i/∂T_i`, per cell.
- `dθ_dξ::Matrix{Float64}` — `∂θ_i/∂ξ_{j,i}`, `(NC, n_cells)`.
- `dθA_dχ::Vector{Float64}` — `∂θ_A/∂χ_j`, length `n_cells`. WithChi
  only; zeros for NoChi.
- `dθA_dT::Vector{Float64}` — `∂θ_A/∂T_j = G(χ̄)·V_j·∂S_dry,j/∂T_j/M`,
  length `n_cells`. WithChi only; zeros for NoChi.
- `dθA_dξ::Matrix{Float64}` — `∂θ_A/∂ξ_{k,j}`, `(NC, n_cells)`. WithChi
  only; zeros for NoChi.
- internal scratch (`rates`, `heats`, `dξ_dt`) for the apply!.
"""
# Tag for the strain-rate stack-Dual sweep — keeps these duals from
# colliding with any outer AD that may be in flight.
struct _ALEStrainTag end

mutable struct ALEMeshVelocityWorkspace{Tu,NC,NR,D}
    w_node::Vector{Float64}
    theta::Vector{Float64}
    theta_A::Base.RefValue{Float64}

    dθ_dT::Vector{Float64}
    dθ_dξ::Matrix{Float64}
    dθA_dχ::Vector{Float64}
    dA_dχ::Vector{Float64}
    dθA_dT::Vector{Float64}
    dθA_dξ::Matrix{Float64}

    # Cumulative Δz from bottom to each node (length n_nodes). Populated by
    # `apply!` and reused by the consumer-side z-row / advection χ-coupling.
    cumZ_node::Vector{Float64}

    rates::Vector{Tu}
    heats::Vector{Tu}
    dξ_dt::Vector{Tu}

    # Dual-typed scratch for `compute_strain_rate_derivatives!`'s stack-Dual
    # sweep: `D = Dual{_ALEStrainTag, Float64, NC + 1}`. Preallocated here so
    # the per-cell derivative pass is allocation-free.
    rates_d::Vector{D}
    heats_d::Vector{D}
    dξ_dt_d::Vector{D}

    # Persistent bare prefix-sum couplings for `contribute_to_structured!`:
    # entry 1 = T channel, entries 1+j = ξ_j channels. The `.v` payloads are
    # refilled in place every Jacobian assembly; downstream embedded
    # couplings alias the same vectors, so no per-eval copies are needed.
    prefix_couplings::Vector{PrefixSumCoupling{Float64}}
end

"""
	build_ale_mesh_velocity_workspace(material, n_cells, n_nodes) -> ALEMeshVelocityWorkspace

Allocate the workspace once at problem-build time. Buffers are sized
from the material (`NC`, `NR`) and the mesh (`n_cells`, `n_nodes`).
"""
function build_ale_mesh_velocity_workspace(
    material::Material{NC,NR},
    n_cells::Integer,
    n_nodes::Integer,
) where {NC,NR}
    D = ForwardDiff.Dual{_ALEStrainTag,Float64,NC + 1}
    prefix_couplings = [PrefixSumCoupling(zeros(Float64, n_cells)) for _ = 1:(1+NC)]
    return ALEMeshVelocityWorkspace{Float64,NC,NR,D}(
        zeros(Float64, n_nodes),
        zeros(Float64, n_cells),
        Ref(0.0),
        zeros(Float64, n_cells),
        zeros(Float64, NC, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, n_cells),
        zeros(Float64, NC, n_cells),
        zeros(Float64, n_nodes),
        zeros(Float64, NR),
        zeros(Float64, NR),
        zeros(Float64, NC),
        zeros(D, NR),
        zeros(D, NR),
        zeros(D, NC),
        prefix_couplings,
    )
end

"""
	apply!(ws::ALEMeshVelocityWorkspace, op::ALEMeshVelocityOperator,
	       state, mesh, material, t; chi_mode = NoChi(),
	       χ̄ = nothing, M_ξ = nothing,
	       layout = nothing, dry_gas_indices = nothing) -> nothing

Populate the derived-field buffers on `ws` from `state`. Wraps the
existing [`_compute_ale_mesh_velocity!`](@ref) path so the math is
exactly the same one the residual sees — there is no duplicate
implementation of the cumulative integral.

Computes:
1. `θ_A` (zero in NoChi; column-global in WithChi).
2. `w_node[k]` (cumulative from bottom) and `θ[i]` (per-cell strain).
3. `cumZ_node[k]` (cumulative Δz from the bottom).

The species source terms are evaluated per cell from the state.

# Required keyword arguments by mode

- **NoChi:** none — `θ_A ≡ 0`, the operator falls through.
- **WithChi:** `χ̄`, `M_ξ`, `layout`, `dry_gas_indices`.
"""
function apply!(
    ws::ALEMeshVelocityWorkspace{Tu,NC,NR},
    ::ALEMeshVelocityOperator,
    state,
    mesh,
    material::Material{NC,NR},
    t::Real;
    chi_mode = NoChi(),
    χ̄::Union{Nothing,Real} = nothing,
    M_ξ::Union{Nothing,Real} = nothing,
    layout = nothing,
    dry_gas_indices = nothing,
) where {Tu,NC,NR}

    # 1) Column-global θ_A (zero for NoChi), computing S_dry per cell on the
    # fly via `ws.rates / ws.heats / ws.dξ_dt`.
    θ_A_value = if chi_mode isa WithChi
        χ̄ === nothing && error("ALEMeshVelocityOperator(WithChi): χ̄ required.")
        M_ξ === nothing && error("ALEMeshVelocityOperator(WithChi): M_ξ required.")
        layout === nothing &&
            error("ALEMeshVelocityOperator(WithChi): layout required.")
        dry_gas_indices === nothing &&
            error("ALEMeshVelocityOperator(WithChi): dry_gas_indices required.")
        _compute_theta_A_inplace(
            material,
            state,
            χ̄,
            M_ξ,
            layout,
            dry_gas_indices,
            ws.rates,
            ws.heats,
            ws.dξ_dt,
        )
    else
        0.0
    end
    ws.theta_A[] = Float64(θ_A_value)

    # 2) Cumulative mesh velocity + per-cell θ. The shared physics
    # function writes both ws.w_node (cumulative from bottom) and
    # ws.theta (per-cell strain), recomputing the species source terms via
    # ws.rates / ws.heats / ws.dξ_dt.
    _compute_ale_mesh_velocity!(
        ws.w_node,
        state,
        mesh,
        material;
        θ_out = ws.theta,
        rates = ws.rates,
        heats = ws.heats,
        dξ_dt = ws.dξ_dt,
        θ_A = ws.theta_A[],
    )

    # 3) Cumulative Δz from bottom to each node — reused by both the rank-1
    # χ builder (u_node[k] = -cumZ[k]) and the consumer-side z-row × χ-col
    # composition. cumZ[1] = 0; cumZ[k] = Σ_{i<k} Δzᵢ.
    A_section = Float64(state.A_section[])
    ws.cumZ_node[1] = 0.0
    @inbounds for k = 2:length(ws.cumZ_node)
        Δz_prev = state.volume[k-1] / A_section
        ws.cumZ_node[k] = ws.cumZ_node[k-1] + Δz_prev
    end

    return nothing
end

# In-place θ_A computation when no precomputed dξ_rxn matrix is available.
# Mirrors `Pyrolysis.Problem._compute_theta_A` but evaluates the per-cell
# species source on the fly using caller-supplied scratch (`rates_buf`,
# `heats_buf`, `dξ_buf`) so it's allocation-free.
function _compute_theta_A_inplace(
    material::Material{NC,NR},
    state,
    χ̄::Real,
    M::Real,
    l,
    idx_dry,
    rates_buf::AbstractVector,
    heats_buf::AbstractVector,
    dξ_buf::AbstractVector,
) where {NC,NR}
    law = material.lateral_shrinkage_law
    law === nothing && return zero(typeof(χ̄))
    M <= eps() && return zero(typeof(χ̄))
    n = l.n_cells
    n == 0 && return zero(typeof(χ̄))

    dχ̄_dt = zero(eltype(state.T))
    @inbounds for i = 1:n
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        fill!(rates_buf, zero(eltype(rates_buf)))
        fill!(heats_buf, zero(eltype(heats_buf)))
        evaluate_reactions!(rates_buf, heats_buf, ξ_i, T_i, material)
        species_source_terms!(dξ_buf, rates_buf, material)
        S_dry = zero(eltype(dξ_buf))
        for j in idx_dry
            S_dry += dξ_buf[j]
        end
        dχ̄_dt += S_dry * state.volume[i]
    end
    dχ̄_dt /= M

    A_ratio = law(χ̄)
    A_ratio > eps() || return zero(typeof(χ̄))
    dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)
    return dlaw_dχ̄ * dχ̄_dt / A_ratio
end

# ============================================================================
# Per-cell strain-rate derivatives (∂θ/∂T, ∂θ/∂ξ) via cell-local AD
# ============================================================================

"""
	compute_strain_rate_derivatives!(ws::ALEMeshVelocityWorkspace,
	                                 state, mesh, material) -> nothing

Populate `ws.dθ_dT` and `ws.dθ_dξ` for the live cells from the current
`state` using cell-local ForwardDiff over `θ_cell(T, ξ) = γ ⨁ S(T, ξ) / ρ(T)`;
parked ghost slots are zero. The
existing analytical kernel
[`Pyrolysis.Jacobian.compute_strain_rate_derivatives!`](@ref) computes
the same quantity by hand-rolled chain rule; using AD here keeps the
operator drift-free and keeps the analytical kernel reachable as a
perf override.

The closure is type-stable on the `Real` superparameter so ForwardDiff
Duals propagate cleanly. Per-cell cost is O(NC + 1) — same scaling as
the existing analytical kernel.
"""
function compute_strain_rate_derivatives!(
    ws::ALEMeshVelocityWorkspace{Tu,NC,NR},
    state,
    mesh,
    material::Material{NC,NR},
) where {Tu,NC,NR}

    n_cells = mesh.n_active_cells
    rates_d = ws.rates_d
    heats_d = ws.heats_d
    dξ_dt_d = ws.dξ_dt_d
    @inbounds for i = (n_cells+1):length(ws.dθ_dT)
        ws.dθ_dT[i] = 0.0
        for j = 1:NC
            ws.dθ_dξ[j, i] = 0.0
        end
    end
    @inbounds for i = 1:n_cells
        # Seed (T_i, ξ_{1,i}, …, ξ_{NC,i}) as stack Duals with unit partials
        # (same pattern as the conduction operator) and run the
        # θ_cell(T, ξ) chain once; the gradient falls out of the partials.
        # Scratch buffers are Dual-typed and preallocated on the workspace,
        # so the whole sweep is allocation-free.
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        T_d, ξ_d = _stack_dual_seed(_ALEStrainTag, state.T[i], ξ_i)
        evaluate_reactions!(rates_d, heats_d, ξ_d, T_d, material)
        species_source_terms!(dξ_dt_d, rates_d, material)
        θ_d = volume_change_rate(dξ_dt_d, T_d, material)

        p = ForwardDiff.partials(θ_d)
        ws.dθ_dT[i] = p[1]
        for j = 1:NC
            ws.dθ_dξ[j, i] = p[1+j]
        end
    end
    return nothing
end

# ============================================================================
# Column-global χ-coupling scalars (∂A/∂χⱼ, ∂θ_A/∂χⱼ) — WithChi only
# ============================================================================

"""
	compute_theta_A_derivatives!(ws, state, mesh, material, layout, dry_gas_indices)
	    -> θ_A_value

Populate `ws.dA_dχ[j]`, `ws.dθA_dχ[j]`, `ws.dθA_dT[j]`, and
`ws.dθA_dξ[k, j]` for `j = 1:n_cells` via the column-global
`(χ̄, T, ξ) → θ_A` chain:

```
M       = Σⱼ m_init,ⱼ
χ̄       = (Σⱼ χⱼ · m_init,ⱼ) / M
A       = A_0 · law(χ̄)
sum_S_V = Σᵢ S_dry,i(T_i, ξ_i) · V_i
dχ̄/dt   = sum_S_V / M
G(χ̄)    = law'(χ̄) / law(χ̄)
G'(χ̄)   = (law''·law − law'²) / law²

∂A/∂χⱼ      = A_0 · law'(χ̄) · m_init,ⱼ / M
∂θ_A/∂χⱼ    = G' · (m_init,ⱼ/M) · dχ̄/dt
            + G · (∂A/∂χⱼ / M) · sum_S_V/A
∂θ_A/∂Tⱼ    = G · V_j · ∂S_dry,j/∂Tⱼ / M
∂θ_A/∂ξ_{k,j} = G · V_j · ∂S_dry,j/∂ξ_{k,j} / M
```

The `(T, ξ)` channels are the θ_A back-reaction that every `w_node[k]`
(and hence every z-row and every ALE-advection row) carries under
WithChi; the provider embeds them as rank-1 couplings.

Returns the current `θ_A` value so the caller can avoid a recompute. For
identity-law / zero-mass / collapsed-A_ratio cases, all arrays are zeroed
and the function returns `0.0`. Computes per-cell dry-gas reaction sources
(value + cell-local gradient) on the fly.
"""
function compute_theta_A_derivatives!(
    ws::ALEMeshVelocityWorkspace{Tu,NC,NR},
    state,
    mesh,
    material::Material{NC,NR},
    layout,
    dry_gas_indices,
) where {Tu,NC,NR}

    dA_dχ = ws.dA_dχ
    dθA_dχ = ws.dθA_dχ
    dθA_dT = ws.dθA_dT
    dθA_dξ = ws.dθA_dξ
    fill!(dA_dχ, 0.0)
    fill!(dθA_dχ, 0.0)
    fill!(dθA_dT, 0.0)
    fill!(dθA_dξ, 0.0)

    law = material.lateral_shrinkage_law
    law === nothing && return 0.0

    m_init = mesh.m_dry_solid_initial
    n_cells = mesh.n_active_cells
    (n_cells == 0 || length(m_init) < n_cells) && return 0.0

    M = 0.0
    χ̄ = 0.0
    @inbounds for j = 1:n_cells
        m_j = m_init[j]
        χ̄ += state.χ[j] * m_j
        M += m_j
    end
    M <= eps() && return 0.0
    inv_M = 1.0 / M
    χ̄ *= inv_M

    A = Float64(state.A_section[])
    A_0 = mesh.cross_section_area_initial

    # sum_S_V = Σᵢ S_dry,i · V_i. Per cell, evaluate S_dry and its cell-local
    # gradient ∂S_dry,i/∂(T_i, ξ_{·,i}) in one ForwardDiff pass; the gradient
    # feeds the θ_A → (T, ξ) back-reaction channels (dθA_dT / dθA_dξ below).
    sum_S_V = 0.0
    @inbounds for i = 1:n_cells
        x = zeros(1 + NC)
        x[1] = state.T[i]
        for j = 1:NC
            x[1+j] = state.ξ[j, i]
        end

        f = let material = material, dry_gas_indices = dry_gas_indices
            (x_inner) -> begin
                T_d = x_inner[1]
                ξ_d = ntuple(j -> x_inner[1+j], Val(NC))
                rates_d = zeros(eltype(x_inner), NR)
                heats_d = zeros(eltype(x_inner), NR)
                dξ_dt_d = zeros(eltype(x_inner), NC)
                evaluate_reactions!(rates_d, heats_d, ξ_d, T_d, material)
                species_source_terms!(dξ_dt_d, rates_d, material)
                S_dry = zero(eltype(x_inner))
                for j in dry_gas_indices
                    S_dry += dξ_dt_d[j]
                end
                S_dry
            end
        end

        grad = ForwardDiff.gradient(f, x)
        V_i = Float64(state.volume[i])
        sum_S_V += f(x) * V_i
        # Stash V_i · ∂S_dry,i/∂x — scaled by G/M once G is known below.
        dθA_dT[i] = V_i * grad[1]
        for j = 1:NC
            dθA_dξ[j, i] = V_i * grad[1+j]
        end
    end
    dχ̄_dt = sum_S_V * inv_M

    # law(χ̄), law'(χ̄), law''(χ̄) via Dual-of-Dual (single eval).
    dual_inner = ForwardDiff.Dual(χ̄, 1.0)
    dual_outer = ForwardDiff.Dual(dual_inner, ForwardDiff.Dual(1.0, 0.0))
    law_dd = law(dual_outer)
    law_val = ForwardDiff.value(ForwardDiff.value(law_dd))
    law_p1 = ForwardDiff.value(ForwardDiff.partials(law_dd, 1))
    law_p2 = ForwardDiff.partials(ForwardDiff.partials(law_dd, 1), 1)

    if !(law_val > eps())
        fill!(dθA_dT, 0.0)
        fill!(dθA_dξ, 0.0)
        return 0.0
    end

    G = law_p1 / law_val
    G_p1 = (law_p2 * law_val - law_p1 * law_p1) / (law_val * law_val)

    sum_S_L_over_M = sum_S_V / (A * M)
    coef_TA = G * inv_M
    @inbounds for j = 1:n_cells
        w_j = m_init[j] * inv_M
        dA_dχ[j] = A_0 * law_p1 * w_j
        dθA_dχ[j] = G_p1 * w_j * dχ̄_dt + G * dA_dχ[j] * sum_S_L_over_M
        dθA_dT[j] *= coef_TA
        for k = 1:NC
            dθA_dξ[k, j] *= coef_TA
        end
    end

    return G * dχ̄_dt
end

# ============================================================================
# Framework hooks
# ============================================================================

"""
	contribute_to_structured!(jac, op::ALEMeshVelocityOperator, state, t,
	                          mesh_vel_ws, layout, chi_mode_tag, n_global;
	                          assembly_profile = :full) -> nothing

The mesh-velocity provider's structured-coupling contribution. Two
responsibilities:

1. **z-row × (T, ξ) block.** The z-row equations `dz_k/dt = w_node[k] =
   Σ_{i<k}(θ_i − θ_A)·Δzᵢ` are the consumers of the cumulative coupling:
   `∂(dz_k/dt)/∂x_j = 𝟙[j ≤ k-1]·Δzⱼ·∂θⱼ/∂xⱼ`. The bare prefix-sum
   primitives (`v = Δzᵢ · ∂θ_i/∂x_i` per source channel, entry 1 = T,
   entry 1+j = ξ_j) are embedded at the z-row offsets as
   [`EmbeddedNodePrefixSumCoupling`](@ref) blocks (rows = nodes,
   lower-triangular). The geometry AD block covers z *columns* only and
   never builds these rows.

2. **WithChi θ_A back-reaction.** When `chi_mode_tag isa WithChi`, θ_A
   couples every `w_node[k]` to the full χ vector *and*, through
   `dχ̄/dt = Σᵢ S_dry,i(T_i, ξ_i)·V_i/M`, to every cell's `(T, ξ)`:

   ```
   ∂w_node[k]/∂y_j = −cumZ_node[k] · ∂θ_A/∂y_j,   y ∈ {χ, T, ξ_1..ξ_NC}
   ```

   Each channel lands as a z-row [`EmbeddedRank1Coupling`](@ref).

`compute_strain_rate_derivatives!` (and, for WithChi,
`compute_theta_A_derivatives!`) must have populated the workspace for the
current state before this method is called — the operator-list dispatcher
in `structured_backend.jl` arranges that.
"""
function contribute_to_structured!(
    jac::StructuredJacobian,
    ::ALEMeshVelocityOperator,
    state,
    ::Real,
    mesh_vel_ws::ALEMeshVelocityWorkspace,
    layout,
    chi_mode_tag,
    n_global::Integer;
    assembly_profile::Symbol = :full,
)
    # Under :ale_reduced no entry on `jac.structured` is materialized
    # (everything this method pushes onto `jac.structured` is dropped by
    # `_emit_jacobian!`). Skip the producer entirely.
    assembly_profile === :full || return nothing

    NC = size(mesh_vel_ws.dθ_dξ, 1)
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes
    z_off = z_offset(layout)

    # Refill the persistent bare prefix-sum couplings in place (entry 1 = T
    # channel, 1+j = ξ_j): `v[i] = Δzᵢ · ∂θᵢ/∂xᵢ`, written into
    # workspace-owned buffers so the hot path allocates nothing. The strain-rate derivatives
    # are zero in parked ghost slots, so the prefix sums stop at the live
    # surface and every ghost node row reproduces the surface node's row.
    couplings = mesh_vel_ws.prefix_couplings
    A_section = state.A_section[]
    v_T_pc = mesh_vel_ws.prefix_couplings[1].v
    @inbounds for i = 1:n_cells
        Δzᵢ = state.volume[i] / A_section
        v_T_pc[i] = Δzᵢ * mesh_vel_ws.dθ_dT[i]
    end
    for j = 1:NC
        v_j = mesh_vel_ws.prefix_couplings[1+j].v
        @inbounds for i = 1:n_cells
            Δzᵢ = state.volume[i] / A_section
            v_j[i] = Δzᵢ * mesh_vel_ws.dθ_dξ[j, i]
        end
    end

    # z-row × (T, ξ) block: embed the bare prefix sums at the z-row offsets.
    # Row block = nodes 1..n_nodes (node 1 receives nothing — w_node[1] ≡ 0).
    compose_into!(jac, couplings[1], z_off, 0, n_nodes, Int(n_global))
    for k = 1:NC
        ξk_col_off = n_cells + (k - 1) * n_cells
        compose_into!(jac, couplings[1+k], z_off, ξk_col_off, n_nodes, Int(n_global))
    end

    if chi_mode_tag isa WithChi
        # Bare rank-1 u factor shared by every θ_A channel:
        # u[k] = -cumZ_node[k] (so node 1's row gets exactly zero).
        u_bare = Vector{Float64}(undef, n_nodes)
        @inbounds for k = 1:n_nodes
            u_bare[k] = -mesh_vel_ws.cumZ_node[k]
        end

        # χ channel: z-row × χ-col. Aliases the workspace buffer (couplings
        # are read-only and rebuilt every eval).
        v_chi = mesh_vel_ws.dθA_dχ
        χ_off = chi_offset(layout)
        push!(
            jac.structured,
            EmbeddedRank1Coupling(u_bare, v_chi, z_off, χ_off, Int(n_global)),
        )

        # (T, ξ) channels — the θ_A back-reaction, embedded at the z-rows.
        v_T = mesh_vel_ws.dθA_dT
        push!(jac.structured, EmbeddedRank1Coupling(u_bare, v_T, z_off, 0, Int(n_global)))
        for k = 1:NC
            v_ξk = Vector{Float64}(undef, n_cells)
            @inbounds for j = 1:n_cells
                v_ξk[j] = mesh_vel_ws.dθA_dξ[k, j]
            end
            ξk_col_off = n_cells + (k - 1) * n_cells
            push!(
                jac.structured,
                EmbeddedRank1Coupling(u_bare, v_ξk, z_off, ξk_col_off, Int(n_global)),
            )
        end
    end

    return nothing
end
