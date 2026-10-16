"""
Unified ODE right-hand side for the pyrolysis package.

```
function residual!(du, u, def::ProblemDef, t, ws::Workspace) -> nothing
```

This is the **single** RHS for all non-P1 simulations. The trait dispatch
chain is:

	residual!(du, u, def, t, ws)
		→ _residual!(du, u, def, t, ws, mesh_mode(def.mode))
			::Eulerian → Eulerian path
			::ALE      → ALE path

Inside each path, the gas transport law follows the property cache (a Darcy
sub-cache exists only in `DarcyFick` mode) and the χ block is included iff
`chi_mode(def.mode) isa WithChi`.

P1 radiation has its own residual in `Pyrolysis.Experimental.RadiationP1`
and is **not** routed through this entry point.

# Contract
- `def` is **immutable** — safe to share across threads / sweep workers.
- `ws` is **mutable** — one per worker.
- `t` is a plain `Real`; the kernel must accept whatever scalar type
  `du`/`u` carry (`Float64`, `ForwardDiff.Dual`, …).
"""
@inline function residual!(
    du::AbstractVector,
    u::AbstractVector,
    def::ProblemDef,
    t::Real,
    ws::Workspace,
)
    _residual!(du, u, def, t, ws, mesh_mode(def.mode))
    return nothing
end

# ============================================================================
# Eulerian residual
# ============================================================================

function _residual!(du, u, def::ProblemDef, t, ws::Workspace, ::Eulerian)
    mesh = ws.mesh
    material = def.material
    bc_set = def.bcs
    layout = ws.layout
    n = layout.n_cells
    NC = n_components(layout)
    radiation_model = get_option(def, :radiation_model, NO_RADIATION)
    radiation_kirchhoff = get_option(def, :radiation_kirchhoff, true)

    # 1) Populate the typed state cache from u — single source of truth
    #    for state during residual assembly. Eulerian geometry is constant
    #    (mirrored at workspace construction), so no geometry update needed.
    populate!(ws.state_cache, u, layout)

    # 2) Update material properties at current state
    Properties.update_properties!(ws.prop_cache, ws.state_cache, mesh, material)

    # 3) Physics RHS
    _compute_physics_rhs!(
        ws,
        mesh,
        material,
        bc_set,
        t,
        radiation_model,
        radiation_kirchhoff,
    )

    # 4) Pack the live rows; parked ghost rows are identically zero.
    _pack_cell_rows!(du, ws.dT_work, ws.dξ_work, n, mesh.n_active_cells, NC)

    return nothing
end

# Live cells 1:n_live take their T and ξ rates; ghost slots n_live+1:n are
# inert and get exactly zero, so W = I/(γΔt) − J is the identity there.
@inline function _pack_cell_rows!(du, dT, dξ, n::Int, n_live::Int, NC::Int)
    @inbounds for i = 1:n_live
        du[i] = dT[i]
        for j = 1:NC
            du[n+(j-1)*n+i] = dξ[j, i]
        end
    end
    z = zero(eltype(du))
    @inbounds for i = (n_live+1):n
        du[i] = z
        for j = 1:NC
            du[n+(j-1)*n+i] = z
        end
    end
    return nothing
end

# ============================================================================
# ALE residual
# ============================================================================

function _residual!(du, u, def::ProblemDef, t, ws::Workspace, ::ALE)
    mesh = ws.mesh
    material = def.material
    bc_set = def.bcs
    layout = ws.layout
    n = layout.n_cells
    n_nodes = layout.n_nodes
    NC = n_components(layout)
    radiation_model = get_option(def, :radiation_model, NO_RADIATION)
    radiation_kirchhoff = get_option(def, :radiation_kirchhoff, true)

    # 1) Populate the typed state cache from u (T, ξ, z, χ).
    populate!(ws.state_cache, u, layout)

    # 2) Update χ-driven cross-section area (Tu-clean) and the rest of the
    #    state-derived geometry (volume, A_face, cell_centers, face_positions).
    #    These all carry eltype(u) end-to-end so AD passes through.
    _, χ̄_cur, M_cur = update_cross_section_area!(
        ws.state_cache,
        material,
        mesh,
        u,
        layout,
        chi_mode(def.mode),
    )
    update_state_geometry!(ws.state_cache, mesh, mesh_mode(def.mode))

    # 3) Update properties
    Properties.update_properties!(ws.prop_cache, ws.state_cache, mesh, material)

    # 4) Physics RHS — must run BEFORE mesh velocity so `dξ_rxn` is populated.
    _compute_physics_rhs!(
        ws,
        mesh,
        material,
        bc_set,
        t,
        radiation_model,
        radiation_kirchhoff,
    )

    # 5) Mesh velocity. θ_A_global picks up the column-global radial channel
    #    when WithChi; with NoChi it is exactly 0.
    θ_A_global = _compute_theta_A_dispatched(
        material,
        ws.state_cache,
        χ̄_cur,
        M_cur,
        ws.rhs_cache.dξ_rxn,
        ws.layout,
        ws.dry_gas_indices,
        chi_mode(def.mode),
    )
    _compute_ale_mesh_velocity!(
        ws.w_work,
        ws.state_cache,
        mesh,
        material;
        θ_out = ws.θ_work,
        precomputed_dξ_rxn = ws.rhs_cache.dξ_rxn,
        dξ_dt = ws.dξ_dt_work,
        θ_A = θ_A_global,
    )

    # Verification hook: replace the reaction-driven mesh velocity with a
    # prescribed w(z_node, t) — θ_work above keeps its physical values so the
    # dilation correction is unaffected. See `Workspace.w_override`. Local
    # binding so the Union{Nothing,Function} field narrows across the branch
    # (mutable-struct field loads don't flow through `!==` checks).
    w_ovr = ws.w_override
    if w_ovr !== nothing
        copyto!(ws.w_physical, ws.w_work)
        _apply_mesh_velocity_override!(ws.w_work, u, layout, w_ovr, t)
    end

    # 6) Mesh-relative temperature advection. Every component rides the
    #    skeleton — the mesh IS the skeleton, and the Darcy–Fick flux in the
    #    step-4 RHS is relative to it — so the whole thermal mass ρcp (gas
    #    included) moves with the mesh and T is advected only by the mesh
    #    velocity RELATIVE to the material, w_rel = w_applied − w_physical
    #    (tech-ref §11.2.2). In production w_rel ≡ 0 and nothing fires; the
    #    kernel runs only under the `w_override` verification hook.
    if w_ovr !== nothing
        @inbounds for k in eachindex(ws.w_physical)
            ws.w_physical[k] = ws.w_work[k] - ws.w_physical[k]  # now w_rel
        end
        Discretization.compute_ale_advection_term!(
            ws.ale_advection,
            ws.state_cache,
            mesh,
            ws.w_physical,
        )
        _apply_ale_T_advection!(ws.dT_work, ws.ale_advection, mesh)
    end

    # 7) Dilation. For every component d(ξⱼV)/dt = ω̇ⱼV − [JⱼA]_faces on the
    #    moving cell, so the concentration rate carries −ξⱼθᵢ with the cell's
    #    own strain rate (tech-ref §11.2.3). Gas is included: it is carried
    #    with the pores it occupies and leaves only through Jⱼ.
    _apply_ale_dilation_correction!(ws.dξ_work, ws.θ_work, ws.state_cache, mesh)

    # 8) Pack du. Ghost T/ξ rows are zero; every node — ghost nodes included,
    #    which ride the live surface — takes its mesh velocity.
    _pack_cell_rows!(du, ws.dT_work, ws.dξ_work, n, mesh.n_active_cells, NC)
    zoff = n * (1 + NC)
    @inbounds for i = 1:n_nodes
        du[zoff+i] = ws.w_work[i]
    end

    # 9) χ block (only when WithChi; with NoChi the chi_range is empty)
    _pack_chi_derivative_dispatched!(
        du,
        layout,
        ws.state_cache,
        mesh,
        material,
        ws.rhs_cache.dξ_rxn,
        ws.dry_gas_indices,
        chi_mode(def.mode),
    )

    return nothing
end

# ============================================================================
# Physics RHS
# ============================================================================

# The transport law (Fickian or Darcy–Fick) follows the property cache's Darcy
# sub-cache, which `build_workspace` sizes from the simulation mode.
@inline function _compute_physics_rhs!(
    ws,
    mesh,
    material,
    bc_set,
    t,
    radiation_model,
    radiation_kirchhoff::Bool,
)
    Discretization.compute_rhs!(
        ws.dT_work,
        ws.dξ_work,
        ws.state_cache,
        mesh,
        ws.prop_cache,
        ws.rhs_cache,
        material,
        bc_set,
        t,
        ws.top_face_ids,
        ws.bottom_face_ids;
        I_radiation = ws.I_radiation,
        radiation_model = radiation_model,
        radiation_kirchhoff = radiation_kirchhoff,
    )
end

# ============================================================================
# χ-mode dispatch helpers
# ============================================================================


@inline _compute_theta_A_dispatched(
    material,
    state,
    χ̄,
    M,
    dξ_rxn,
    layout,
    idx_dry,
    ::NoChi,
) = 0.0

@inline _compute_theta_A_dispatched(
    material,
    state,
    χ̄,
    M,
    dξ_rxn,
    layout,
    idx_dry,
    ::WithChi,
) = _compute_theta_A(material, state, χ̄, M, dξ_rxn, layout, idx_dry)

@inline _pack_chi_derivative_dispatched!(
    du,
    layout,
    state,
    mesh,
    material,
    dξ_rxn,
    idx_dry,
    ::NoChi,
) = nothing

@inline _pack_chi_derivative_dispatched!(
    du,
    layout,
    state,
    mesh,
    material,
    dξ_rxn,
    idx_dry,
    ::WithChi,
) = _pack_chi_derivative!(du, layout, state, mesh, material, dξ_rxn, idx_dry)

# ============================================================================
# ALE helpers — layout-keyed
# ============================================================================

@inline function _lateral_area_ratio(material, χ̄)
    law = material.lateral_shrinkage_law
    return law === nothing ? one(χ̄) : law(χ̄)
end

"""
	_compute_ale_mesh_velocity!(w, state, mesh, material; ...)

Cumulative mesh velocity from per-cell volume change of every cell from
the bottom node up. With the default identity lateral-shrinkage law,
the column-global radial channel `θ_A ≡ 0` and behavior reduces to the
constant-cross-section formulation. State (T, ξ, geometry) is read from
the typed `state::StateCache`.
"""
function _compute_ale_mesh_velocity!(
    w::AbstractVector,
    state,
    mesh,
    material::Material;
    θ_out::Union{Nothing,AbstractVector} = nothing,
    precomputed_dξ_rxn::Union{Nothing,AbstractMatrix} = nothing,
    rates::Union{Nothing,AbstractVector} = nothing,
    heats::Union{Nothing,AbstractVector} = nothing,
    dξ_dt::Union{Nothing,AbstractVector} = nothing,
    θ_A::Real = 0.0,
)
    NC = n_components(material)
    NR = n_reactions(material)
    n_cells = mesh.n_active_cells

    use_precomputed = precomputed_dξ_rxn !== nothing
    if !use_precomputed
        Tu = eltype(state.T)
        rates === nothing && (rates = zeros(Tu, NR))
        heats === nothing && (heats = zeros(Tu, NR))
    end
    dξ_dt === nothing && (dξ_dt = zeros(eltype(state.T), NC))

    w[1] = zero(eltype(w))
    cumulative_velocity = zero(eltype(w))

    @inbounds for cell_idx = 1:n_cells
        T_cell = state.T[cell_idx]
        Δz_cell = state.volume[cell_idx] / state.A_section[]

        if use_precomputed
            for j = 1:NC
                dξ_dt[j] = precomputed_dξ_rxn[j, cell_idx]
            end
        else
            ξ_cell = ntuple(j -> state.ξ[j, cell_idx], Val(NC))
            fill!(rates, zero(eltype(rates)))
            fill!(heats, zero(eltype(heats)))
            Physics.evaluate_reactions!(rates, heats, ξ_cell, T_cell, material)
            Physics.species_source_terms!(dξ_dt, rates, material)
        end

        θ_cell = Physics.volume_change_rate(dξ_dt, T_cell, material)
        if θ_out !== nothing
            θ_out[cell_idx] = θ_cell
        end

        # Axial residual after column-global radial channel removed.
        θ_L = θ_cell - θ_A
        cumulative_velocity += θ_L * Δz_cell

        right_node = mesh.cells[cell_idx].node_ids[2]
        w[right_node] = cumulative_velocity
    end

    # Parked ghost nodes sit on the live surface and move with it.
    @inbounds for k = (n_cells+2):length(w)
        w[k] = cumulative_velocity
    end

    return nothing
end

"""
	_apply_mesh_velocity_override!(w, u, layout, override, t)

Replace every node velocity with `override(z_node, t)`, reading node
positions from the live state's z-block. Verification-only path (see
`Workspace.w_override`); called after the physical velocity computation so
`θ_work` retains its reaction-driven values.
"""
function _apply_mesh_velocity_override!(
    w::AbstractVector,
    u::AbstractVector,
    layout,
    override::Function,
    t::Real,
)
    zoff = layout.n_cells * (1 + n_components(layout))
    @inbounds for i in eachindex(w)
        w[i] = override(u[zoff+i], t)
    end
    return nothing
end

function _apply_ale_T_advection!(dT::AbstractVector, ale_advection::AbstractVector, mesh)
    @inbounds for i = 1:mesh.n_active_cells
        dT[i] += ale_advection[i]
    end
    return nothing
end

"""
	_apply_ale_dilation_correction!(dξ, θ, state, mesh)

Add the mesh-motion dilation `−ξⱼθᵢ` to every component's rate: the cell
volume changes at the skeleton strain rate θᵢ and every component, gas
included, is carried with it, so its concentration rises as the cell
shrinks and falls as it swells.
"""
function _apply_ale_dilation_correction!(dξ::AbstractMatrix, θ::AbstractVector, state, mesh)
    NC = size(dξ, 1)
    @inbounds for i = 1:mesh.n_active_cells
        θ_i = θ[i]
        for j = 1:NC
            dξ[j, i] -= state.ξ[j, i] * θ_i
        end
    end
    return nothing
end

"""
	_pack_chi_derivative!(du, layout, mesh, material, dξ_rxn, idx_dry)

Write `du[χ_offset + i] = (Σⱼ ∈ dry-pyrolysis-gas Sⱼ,i) · V_i / m_dry_solid_initial,i`
for each cell i.

Identity-law materials skip the dry-gas summation and zero the χ block —
χ is undefined when the cross-section is fixed.
"""
@inline function _pack_chi_derivative!(
    du::AbstractVector,
    l::StateLayout,
    state,
    mesh,
    material,
    dξ_rxn::AbstractMatrix,
    idx_dry,
)
    rng = chi_range(l)
    isempty(rng) && return nothing
    χ_offset = first(rng) - 1
    n = l.n_cells

    if material.lateral_shrinkage_law === nothing
        @inbounds for i = 1:n
            du[χ_offset+i] = zero(eltype(du))
        end
        return nothing
    end

    m_init = mesh.m_dry_solid_initial
    n_live = mesh.n_active_cells

    @inbounds for i = 1:n_live
        S_dry = zero(eltype(dξ_rxn))
        for j in idx_dry
            S_dry += dξ_rxn[j, i]
        end
        V_i = state.volume[i]
        denom = max(m_init[i], eps())
        du[χ_offset+i] = S_dry * V_i / denom
    end
    @inbounds for i = (n_live+1):n
        du[χ_offset+i] = zero(eltype(du))
    end
    return nothing
end

"""
	_compute_theta_A(material, state, χ̄, M, dξ_rxn, layout, idx_dry) -> Float64

Column-global radial volumetric strain rate

	θ_A = (1/A) · dA/dt = law'(χ̄) · dχ̄/dt / law(χ̄)

implied by the lateral shrinkage law `A(t)/A_0 = law(χ̄)` and the current
state. Returns `0.0` for identity-law materials. Parked ghost cells have
zero volume (ALE) and zero reaction source, so summing over every state slot
is the live-cell sum.
"""
@inline function _compute_theta_A(
    material,
    state,
    χ̄::Real,
    M::Real,
    dξ_rxn::Union{Nothing,AbstractMatrix},
    l::StateLayout,
    idx_dry,
)
    law = material.lateral_shrinkage_law
    law === nothing && return zero(typeof(χ̄))
    dξ_rxn === nothing && return zero(typeof(χ̄))
    M <= eps() && return zero(typeof(χ̄))

    n = l.n_cells
    n == 0 && return zero(typeof(χ̄))

    dχ̄_dt = zero(eltype(dξ_rxn))
    @inbounds for i = 1:n
        S_dry = zero(eltype(dξ_rxn))
        for j in idx_dry
            S_dry += dξ_rxn[j, i]
        end
        dχ̄_dt += S_dry * state.volume[i]
    end
    dχ̄_dt /= M

    A_ratio = law(χ̄)
    A_ratio > eps() || return zero(typeof(χ̄))

    dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)
    return dlaw_dχ̄ * dχ̄_dt / A_ratio
end

"""
	column_chi_bar_from_state(u, layout, mesh) -> Real

Mass-weighted column-averaged χ̄. AD-Dual safe (numerator carries duals,
denominator from the immutable `m_dry_solid_initial`). Thin wrapper over the
single implementation `_chi_bar_from_state`, which also returns the column
dry-solid mass used by `update_cross_section_area!`.
"""
@inline column_chi_bar_from_state(u::AbstractVector, l::StateLayout, mesh) =
    _chi_bar_from_state(u, l, mesh)[1]

"""
	cross_section_area_from_state(mesh, material, u, layout) -> Real

Current column cross-sectional area implied by the χ block. Returns
`mesh.cross_section_area_initial` for identity-law materials.
"""
function cross_section_area_from_state(mesh, material, u::AbstractVector, l::StateLayout)
    χ̄ = column_chi_bar_from_state(u, l, mesh)
    area_ratio = _lateral_area_ratio(material, χ̄)
    return mesh.cross_section_area_initial * area_ratio
end
