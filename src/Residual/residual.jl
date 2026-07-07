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

Inside each path, transport is dispatched via `transport_mode(def.mode)`
(`Fickian` vs `DarcyFick`) and the χ block is included iff
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
    energy_form = get_option(def, :energy_form, :standard)
    use_transpiration_bc = get_option(def, :use_transpiration_bc, false)
    use_blowing_correction = get_option(def, :use_blowing_correction, false)

    # 1) Populate the typed state cache from u — single source of truth
    #    for state during residual assembly. Eulerian geometry is constant
    #    (mirrored at workspace construction), so no geometry update needed.
    populate!(ws.state_cache, u, layout)

    # 2) Update material properties at current state
    Properties.update_properties!(ws.prop_cache, ws.state_cache, mesh, material)

    # 3) Physics RHS — dispatch on transport
    _compute_physics_rhs!(
        ws,
        mesh,
        material,
        bc_set,
        t,
        radiation_model,
        radiation_kirchhoff,
        energy_form,
        use_transpiration_bc,
        use_blowing_correction,
        transport_mode(def.mode),
    )

    # 4) Pack
    @inbounds for i = 1:n
        du[i] = ws.dT_work[i]
        for j = 1:NC
            du[n+(j-1)*n+i] = ws.dξ_work[j, i]
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
    energy_form = get_option(def, :energy_form, :standard)
    use_transpiration_bc = get_option(def, :use_transpiration_bc, false)
    use_blowing_correction = get_option(def, :use_blowing_correction, false)

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
        energy_form,
        use_transpiration_bc,
        use_blowing_correction,
        transport_mode(def.mode),
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
        ws.prop_cache,
        material,
        t,
        θ_out = ws.θ_work,
        precomputed_dξ_rxn = ws.rhs_cache.dξ_rxn,
        rates = ws.rates_work,
        heats = ws.heats_work,
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

    _interpolate_node_to_face!(ws.w_face, ws.w_work, mesh)

    # 6) ALE mesh-motion advection (the Darcy gas flux is already in the
    #    lab-frame RHS from step 4).
    #
    #    Energy: condensed species are mesh-bound (v_solid ≡ w, tech-ref
    #    §11.1.2), so the matrix-carried thermal mass (ρcp_eff is matrix-only)
    #    is advected relative to the mesh only by w_rel = w_applied −
    #    w_physical. In production the mesh velocity IS the material velocity
    #    — w_rel ≡ 0 and the energy kernel does not fire. It fires only under
    #    the `w_override` verification hook, where w_applied ≠ w_physical.
    #
    #    Gas species are NOT mesh-bound: their mesh-motion transport is the
    #    conservative flux form driven by the full applied mesh velocity
    #    (solids get no advection term; see tech-ref §11.2.3 — the same
    #    relative-velocity logic). Boundary sweep fluxes follow the mass BCs:
    #    an impermeable face carries its gas along (compression), an open
    #    face exchanges at the zero-gradient ghost value.
    if w_ovr !== nothing
        @inbounds for k in eachindex(ws.w_physical)
            ws.w_physical[k] = ws.w_work[k] - ws.w_physical[k]  # now w_rel
        end
        Discretization.compute_ale_advection_term!(
            ws.ale_advection,
            ws.state_cache,
            mesh,
            ws.w_physical;
            limiter = ws.ale_limiter,
        )
    end
    Discretization.compute_ale_species_advection!(
        ws.ale_species_advection,
        ws.state_cache,
        mesh,
        material,
        ws.w_work;
        θ_A = θ_A_global,
        top_gas_open = !(bc_set.mass.top isa ImpermeableBC),
        bottom_gas_open = !(bc_set.mass.bottom isa ImpermeableBC),
    )

    # 7) Apply ALE corrections (T advection only when it was computed)
    if w_ovr !== nothing
        _apply_ale_T_advection!(ws.dT_work, ws.ale_advection)
    end
    _apply_ale_species_advection!(ws.dξ_work, ws.ale_species_advection, NC)

    # 8) Dilation correction for mesh-bound (solid/liquid) species. Gas is
    #    NOT included here: its volume-rate term is part of the conservative
    #    flux form in step 6.
    _apply_ale_dilation_correction!(ws.dξ_work, ws.θ_work, ws.state_cache, material)

    # 9) Pack du
    @inbounds for i = 1:n
        du[i] = ws.dT_work[i]
    end
    @inbounds for i = 1:n
        for j = 1:NC
            du[n+(j-1)*n+i] = ws.dξ_work[j, i]
        end
    end
    zoff = n * (1 + NC)
    @inbounds for i = 1:n_nodes
        du[zoff+i] = ws.w_work[i]
    end

    # 10) χ block (only when WithChi; with NoChi the chi_range is empty)
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
# Transport-mode dispatch
# ============================================================================

@inline function _compute_physics_rhs!(
    ws,
    mesh,
    material,
    bc_set,
    t,
    radiation_model,
    radiation_kirchhoff::Bool,
    energy_form::Symbol,
    use_transpiration_bc::Bool,
    use_blowing_correction::Bool,
    ::Fickian,
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
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
    )
end

@inline function _compute_physics_rhs!(
    ws,
    mesh,
    material,
    bc_set,
    t,
    radiation_model,
    radiation_kirchhoff::Bool,
    energy_form::Symbol,
    use_transpiration_bc::Bool,
    use_blowing_correction::Bool,
    ::DarcyFick,
)
    Discretization.compute_rhs_darcy_fick!(
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
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
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
	_compute_ale_mesh_velocity!(w, state, mesh, prop_cache, material, t; ...)

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
    prop_cache,
    material::Material,
    t::Real;
    θ_out::Union{Nothing,AbstractVector} = nothing,
    precomputed_dξ_rxn::Union{Nothing,AbstractMatrix} = nothing,
    rates::Union{Nothing,AbstractVector} = nothing,
    heats::Union{Nothing,AbstractVector} = nothing,
    dξ_dt::Union{Nothing,AbstractVector} = nothing,
    θ_A::Real = 0.0,
)
    NC = n_components(material)
    NR = n_reactions(material)
    n_cells = length(state.T)

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
            Physics.species_source_terms!(dξ_dt, rates, material, NC)
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

function _interpolate_node_to_face!(w_face::AbstractVector, w_nodes::AbstractVector, mesh)
    @inbounds for f in eachindex(mesh.faces)
        face = mesh.faces[f]
        w_face[f] = w_nodes[face.node_id]
    end
    return nothing
end

function _apply_ale_T_advection!(dT::AbstractVector, ale_advection::AbstractVector)
    @inbounds for i in eachindex(dT)
        dT[i] += ale_advection[i]
    end
    return nothing
end

function _apply_ale_species_advection!(
    dξ::AbstractMatrix,
    ale_species_advection::AbstractMatrix,
    NC::Int,
)
    n_cells = size(dξ, 2)
    @inbounds for i = 1:n_cells
        for j = 1:NC
            dξ[j, i] += ale_species_advection[j, i]
        end
    end
    return nothing
end

function _apply_ale_dilation_correction!(
    dξ::AbstractMatrix,
    θ::AbstractVector,
    state,
    material::Material{NC},
) where {NC}
    n_cells = length(state.T)
    @inbounds for i = 1:n_cells
        θ_i = θ[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        _apply_dilation_correction_unrolled!(dξ, i, ξ_i, θ_i, material)
    end
    return nothing
end

@generated function _apply_dilation_correction_unrolled!(
    dξ::AbstractMatrix,
    i::Int,
    ξ_i::NTuple{NC,<:Real},
    θ_i::Real,
    material,
) where {NC}
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if !is_gas(material.components[$j])
                dξ[$j, i] -= ξ_i[$j] * θ_i
            end
        end)
    end
    return Expr(:block, exprs...)
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

    @inbounds for i = 1:n
        S_dry = zero(eltype(dξ_rxn))
        for j in idx_dry
            S_dry += dξ_rxn[j, i]
        end
        V_i = state.volume[i]
        denom = max(m_init[i], eps())
        du[χ_offset+i] = S_dry * V_i / denom
    end
    return nothing
end

"""
	_compute_theta_A(material, mesh, χ̄, M, dξ_rxn, layout, idx_dry) -> Float64

Column-global radial volumetric strain rate

	θ_A = (1/A) · dA/dt = law'(χ̄) · dχ̄/dt / law(χ̄)

implied by the lateral shrinkage law `A(t)/A_0 = law(χ̄)` and the current
state. Returns `0.0` for identity-law materials.
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

    # m_dry_solid_initial is a Float64 vector on mesh — but we need it
    # only for the length comparison and via the caller's M. The volume
    # weights come from state_cache.
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
