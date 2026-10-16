"""
The `Structured` backend: operator-by-operator Jacobian assembly.

`Structured(strategy = DirectSolve())` produces a `JacobianSpec`
backed by:

- a [`StructuredJacobian`](@ref) (sparse `A` + structured couplings);
- per-operator AD workspaces, walked in a fixed order;
- the chosen [`LinearSolveStrategy`](@ref), which fixes the assembly profile.

# Assembly path

Builds per-operator AD workspaces and dispatches each operator's
`contribute_to_jacobian!` directly. The operator-list path covers the
core residual matrix:

- `{Eulerian, ALE}` × `{Fickian, DarcyFick}` ×
  `{NoRadiation, SurfaceAbsorption, BeerLambert}` ×
  `{NoChi, WithChi}`, with `WithChi` under ALE only.

# Backend × strategy is orthogonal

`Structured` is the *backend* (the assembling rule). `DirectSolve` /
`ApproxSolve` is the *strategy*: `DirectSolve` includes the structured
couplings, `ApproxSolve` deliberately drops them.

# Approximation trait

`is_approximate(::Structured{ApproxSolve})` is `true` because the
strategy *drops* the structured-coupling contribution at solve time —
which is what makes `verify_jacobian_at` use masked compare rather
than strict. Every other strategy is exact.
"""

using SparseArrays: SparseMatrixCSC
using ..Materials:
    RadiationModel,
    NO_RADIATION,
    SURFACE_ABSORPTION,
    BEER_LAMBERT,
    n_components,
    is_gas,
    dry_pyrolysis_gas_indices
using ..Problem:
    ProblemDef,
    Workspace,
    StateLayout,
    mesh_mode,
    chi_mode,
    transport_mode,
    get_option,
    populate!,
    update_cross_section_area!,
    update_state_geometry!,
    Fickian,
    DarcyFick
using ..Residual: column_chi_bar_from_state
using ..Properties: update_properties!
# `ALE` the trait (from Problem) collides with `ALE` the submodule, which is
# already in scope at the Jacobian module level. Import the trait under an
# explicit alias so `isa` dispatches against the type, not the module.
import ..Problem as PyProblem

# ============================================================================
# Backend tag
# ============================================================================

"""
	Structured{S<:LinearSolveStrategy} <: JacobianBackend

Backend that assembles a [`StructuredJacobian`](@ref) by walking an
operator list. The strategy parameter `S` fixes the assembly profile
(`sparsity_profile`) and is reported by `is_approximate`.

# Constructor
`Structured(; strategy = DirectSolve())` — the standard form. The
strategy is a constructed instance, not a type, so future variants
can carry per-strategy knobs in fields without touching this signature.
"""
struct Structured{S<:LinearSolveStrategy} <: JacobianBackend
    strategy::S
end

Structured(; strategy::LinearSolveStrategy = DirectSolve()) = Structured(strategy)

# Approximation trait — only `ApproxSolve` strategy is approximate.
is_approximate(::Structured{<:LinearSolveStrategy}) = false
is_approximate(::Structured{ApproxSolve}) = true

# ============================================================================
# Cache types
# ============================================================================

# Every field is concretely typed (note `DefT`/`JacT` for the parametric
# `ProblemDef`/`StructuredJacobian`): an abstractly-typed field here turns
# every operator call in `_compute_jacobian!` into a dynamic dispatch that
# boxes all of its arguments (about 5 KiB per Jacobian evaluation).
"""
	StructuredOpsCache

Operator-list cache. Carries
per-operator AD workspaces, mode flags, and references to `def` /
material / mesh / layout that the operator dispatchers thread into each
`contribute_to_jacobian!` call.

Fields with type `Union{Nothing, …}` correspond to operators that are
mode-conditional (e.g. `darcy_ws` is non-`nothing` only when transport
mode is `DarcyFick`). The dispatcher walks every field and skips the
`nothing` slots; the resulting union dispatch is overshadowed by the
per-operator AD cost.

# Fields
- `def::ProblemDef`, `material`, `bc_set`, `layout`, `mesh` — references
  used by every operator dispatcher.
- `radiation_model::RadiationModel` — extracted from `def.options`.
- `radiation_kirchhoff::Bool` — Beer-Lambert surface-transmissivity
  convention (apply the live surface-cell ε_eff at injection), extracted
  from `def.options`; must match the residual's
  `compute_beer_lambert_source!` call.
- `is_ale::Bool` — toggle for the ALE-only stages (mesh-velocity
  workspace populate, dilation rows, geometry AD).
- Per-operator AD workspaces (`reaction_ws`, `conduction_ws`,
  `gas_transport_ws`, `darcy_ws`, `bc_thermal_ws`, `bc_mass_ws`,
  `bc_convection_ws`, `ale_mesh_vel_ws`).
- `geom_cache::Union{Nothing, GeometryADCache}` — for the ALE z-block
  AD pass that runs after the operator walk. `nothing` under
  `assembly_profile == :ale_reduced` (the geometry block is dropped).
- `assembly_profile::Symbol` — the assembly profile this cache was built
  for (`:full` or `:ale_reduced`). Read by per-step assembly to gate
  structured-coupling emission and the geometry AD pass.
- `jac::StructuredJacobian`.
"""
struct StructuredOpsCache{
    Mat,
    BCS,
    MeshT,
    LayoutT,
    WSReact,
    WSCond,
    WSGT,
    WSD,
    WSBT,
    WSBM,
    WSBC,
    WSAMV,
    WSChi,
    WSBL,
    GC,
    DefT<:ProblemDef,
    JacT<:StructuredJacobian,
}
    def::DefT
    material::Mat
    bc_set::BCS
    layout::LayoutT
    mesh::MeshT
    radiation_model::RadiationModel
    radiation_intensity::Float64
    radiation_kirchhoff::Bool
    is_ale::Bool
    assembly_profile::Symbol

    reaction_ws::WSReact
    conduction_ws::WSCond
    gas_transport_ws::WSGT
    darcy_ws::WSD
    bc_thermal_ws::WSBT
    bc_mass_ws::WSBM
    bc_convection_ws::WSBC
    ale_mesh_vel_ws::WSAMV
    chi_ws::WSChi
    beer_lambert_ws::WSBL
    geom_cache::GC

    jac::JacT
end

# ============================================================================
# Construction
# ============================================================================

function _build_jacobian_spec(
    backend::Structured,
    def::ProblemDef,
    ws::Workspace;
    chunk_size::Int = 8,
    kwargs...,
)
    chunk_size >= 1 || throw(ArgumentError("chunk_size must be positive"))

    profile = sparsity_profile(backend.strategy)
    S_phys, S_geom = build_sparsity(ws.layout, def.mode; profile = profile)
    S_union = S_phys .| S_geom
    A = SparseMatrixCSC{Float64,Int}(
        S_union.m,
        S_union.n,
        copy(S_union.colptr),
        copy(S_union.rowval),
        zeros(Float64, length(S_union.nzval)),
    )
    jac = StructuredJacobian(A)

    cache =
        _build_ops_cache(def, ws, jac; chunk_size = chunk_size, assembly_profile = profile)
    return JacobianSpec(backend, S_phys, S_geom, cache)
end

function _build_ops_cache(
    def::ProblemDef,
    ws::Workspace,
    jac;
    chunk_size::Int,
    assembly_profile::Symbol = :full,
)
    material = def.material
    layout = ws.layout
    mesh = ws.mesh
    bc_set = def.bcs
    is_ale = mesh_mode(def.mode) isa PyProblem.ALE
    is_darcy = transport_mode(def.mode) isa DarcyFick
    is_with_chi = chi_mode(def.mode) isa PyProblem.WithChi
    rad = get_option(def, :radiation_model, NO_RADIATION)
    is_beer_lambert = rad == BEER_LAMBERT
    radiation_intensity = Float64(get_option(def, :radiation_intensity, 0.0))
    radiation_kirchhoff = get_option(def, :radiation_kirchhoff, true)
    is_reduced = assembly_profile === :ale_reduced

    # Required workspaces — every supported mode needs these.
    reaction_ws = build_reaction_source_workspace(material)
    conduction_ws = build_conduction_workspace(material)
    bc_thermal_ws = build_bc_thermal_workspace(material)
    bc_mass_ws = build_bc_mass_workspace(material)
    bc_convection_ws = build_bc_convection_workspace(material)

    # Transport — exactly one of these is non-nothing.
    gas_transport_ws = is_darcy ? nothing : build_gas_transport_workspace(material)
    darcy_ws = is_darcy ? build_darcy_transport_workspace(material) : nothing

    # ALE workspaces. `ChiOperator` reads `dA_dχ` from the mesh-velocity
    # workspace; WithChi exists only under ALE (`to_problem_def`).
    ale_mesh_vel_ws =
        is_ale ?
        build_ale_mesh_velocity_workspace(material, layout.n_cells, layout.n_nodes) :
        nothing
    chi_ws = is_with_chi ? build_chi_workspace(material) : nothing
    beer_lambert_ws =
        is_beer_lambert ? build_beer_lambert_workspace(material, layout.n_cells) : nothing
    # Geometry AD cache: built only when the geometry block participates in
    # the assembly. Under `:ale_reduced` the z-column couplings are dropped
    # from sparsity entirely, so building the cache would be pure waste.
    geom_cache =
        (is_ale && !is_reduced) ? GeometryADCache(def, ws; chunk_size = chunk_size) :
        nothing

    return StructuredOpsCache{
        typeof(material),
        typeof(bc_set),
        typeof(mesh),
        typeof(layout),
        typeof(reaction_ws),
        typeof(conduction_ws),
        typeof(gas_transport_ws),
        typeof(darcy_ws),
        typeof(bc_thermal_ws),
        typeof(bc_mass_ws),
        typeof(bc_convection_ws),
        typeof(ale_mesh_vel_ws),
        typeof(chi_ws),
        typeof(beer_lambert_ws),
        typeof(geom_cache),
        typeof(def),
        typeof(jac),
    }(
        def,
        material,
        bc_set,
        layout,
        mesh,
        rad,
        radiation_intensity,
        radiation_kirchhoff,
        is_ale,
        assembly_profile,
        reaction_ws,
        conduction_ws,
        gas_transport_ws,
        darcy_ws,
        bc_thermal_ws,
        bc_mass_ws,
        bc_convection_ws,
        ale_mesh_vel_ws,
        chi_ws,
        beer_lambert_ws,
        geom_cache,
        jac,
    )
end

# ============================================================================
# Compute
# ============================================================================

function _compute_jacobian!(
    J::SparseMatrixCSC{Float64,Int},
    b::Structured,
    cache::StructuredOpsCache,
    sparsity_phys,
    sparsity_geom,
    u,
    def,
    t,
    ws,
)
    zero!(cache.jac)

    # 1. Populate the typed state cache from `u` (Float64), refresh ALE
    #    geometry, then properties. Mirrors the order in `Problem._residual!`.
    populate!(ws.state_cache, u, ws.layout)
    if cache.is_ale
        update_cross_section_area!(
            ws.state_cache,
            cache.material,
            cache.mesh,
            u,
            ws.layout,
            chi_mode(def.mode),
        )
        update_state_geometry!(ws.state_cache, cache.mesh, mesh_mode(def.mode))
    end
    update_properties!(ws.prop_cache, ws.state_cache, cache.mesh, cache.material)

    # 2a. ALE prep: populate the mesh-velocity workspace (w_node, theta,
    #     cumZ) and cell-local strain-rate derivatives (∂θ/∂T, ∂θ/∂ξ). The
    #     dilation rows read the derivatives; the provider-side
    #     `contribute_to_structured!` emits the z-row × (T, ξ) prefix-sum
    #     blocks onto `jac.structured` and (for WithChi) the θ_A rank-1
    #     z-row embeddings.
    if cache.is_ale
        chi_mode_tag = chi_mode(def.mode)
        if chi_mode_tag isa PyProblem.WithChi
            χ̄ = column_chi_bar_from_state(u, cache.layout, cache.mesh)
            M = sum(cache.mesh.m_dry_solid_initial)
            idx_dry = dry_pyrolysis_gas_indices(cache.material)
            apply!(
                cache.ale_mesh_vel_ws,
                ALEMeshVelocityOperator(),
                ws.state_cache,
                cache.mesh,
                cache.material,
                t;
                chi_mode = chi_mode_tag,
                χ̄ = χ̄,
                M_ξ = M,
                layout = cache.layout,
                dry_gas_indices = idx_dry,
            )
            compute_theta_A_derivatives!(
                cache.ale_mesh_vel_ws,
                ws.state_cache,
                cache.mesh,
                cache.material,
                cache.layout,
                idx_dry,
            )
        else
            apply!(
                cache.ale_mesh_vel_ws,
                ALEMeshVelocityOperator(),
                ws.state_cache,
                cache.mesh,
                cache.material,
                t,
            )
        end
        compute_strain_rate_derivatives!(
            cache.ale_mesh_vel_ws,
            ws.state_cache,
            cache.mesh,
            cache.material,
        )
        contribute_to_structured!(
            cache.jac,
            ALEMeshVelocityOperator(),
            ws.state_cache,
            t,
            cache.ale_mesh_vel_ws,
            cache.layout,
            chi_mode_tag,
            size(cache.jac.A, 1);
            assembly_profile = cache.assembly_profile,
        )
    end

    # 3. Walk the operators in fixed order. Each branch threads only the
    #    args its operator needs.
    contribute_to_jacobian!(
        cache.jac,
        ReactionSourceOperator(),
        u,
        t,
        cache.reaction_ws,
        cache.material,
        cache.layout,
        cache.mesh,
    )

    contribute_to_jacobian!(
        cache.jac,
        ConductionOperator(),
        u,
        t,
        cache.conduction_ws,
        ws.state_cache,
        cache.mesh,
        cache.material,
        cache.layout,
    )

    if cache.darcy_ws !== nothing
        contribute_to_jacobian!(
            cache.jac,
            DarcyTransportOperator(),
            u,
            t,
            cache.darcy_ws,
            ws.state_cache,
            cache.mesh,
            cache.material,
            cache.layout,
        )
    else
        contribute_to_jacobian!(
            cache.jac,
            GasTransportOperator(),
            u,
            t,
            cache.gas_transport_ws,
            ws.state_cache,
            cache.mesh,
            cache.material,
            cache.layout,
        )
    end

    # Radiation — no-op for NO_RADIATION / SURFACE_ABSORPTION. For
    # BEER_LAMBERT we populate the workspace's q_rad / I_in / α-derivatives
    # first, then dispatch to the workspace overload.
    if cache.beer_lambert_ws !== nothing
        apply!(
            cache.beer_lambert_ws,
            RadiationOperator(),
            ws.state_cache,
            cache.mesh,
            ws.prop_cache,
            cache.material,
            cache.radiation_intensity,
            cache.radiation_model,
            cache.radiation_kirchhoff,
        )
        compute_alpha_derivatives!(
            cache.beer_lambert_ws,
            ws.state_cache,
            cache.mesh,
            ws.prop_cache,
            cache.material,
        )
        contribute_to_jacobian!(
            cache.jac,
            RadiationOperator(),
            u,
            t,
            cache.beer_lambert_ws,
            ws.state_cache,
            cache.mesh,
            ws.prop_cache,
            cache.material,
            cache.layout,
            cache.radiation_intensity,
            cache.radiation_model,
            cache.radiation_kirchhoff;
            assembly_profile = cache.assembly_profile,
        )
    else
        contribute_to_jacobian!(
            cache.jac,
            RadiationOperator(),
            u,
            t,
            nothing,
            cache.radiation_model,
        )
    end

    contribute_to_jacobian!(
        cache.jac,
        BCThermalOperator(),
        u,
        t,
        cache.bc_thermal_ws,
        ws.state_cache,
        cache.mesh,
        cache.material,
        cache.layout,
        cache.bc_set,
        cache.radiation_model,
        cache.darcy_ws !== nothing,  # blowing's ṁ″ includes Darcy venting
    )

    contribute_to_jacobian!(
        cache.jac,
        BCMassOperator(),
        u,
        t,
        cache.bc_mass_ws,
        ws.state_cache,
        cache.mesh,
        cache.material,
        cache.layout,
        cache.bc_set,
        cache.darcy_ws !== nothing,  # venting term exists only on the Darcy path
    )

    # Boundary-cell S_conv energy row — single owner of the J̄ × ∂T/∂x(T_s)
    # product spanning the thermal-BC Newton, the mass-BC face flux, and the
    # mode-appropriate interior face flux.
    contribute_to_jacobian!(
        cache.jac,
        BCConvectionOperator(),
        u,
        t,
        cache.bc_convection_ws,
        ws.state_cache,
        cache.mesh,
        cache.material,
        cache.layout,
        cache.bc_set,
        cache.darcy_ws !== nothing,  # interior-face kernel + venting follow the mode
        cache.radiation_model,
    )

    if cache.is_ale
        # ALE dilation for every ξ row. Mirrors
        # `_apply_ale_dilation_correction!` in `Residual.residual!` —
        # `dξ[j, i] -= ξ[j, i] · θ[i]` for every component. The cell-local
        # Jacobian:
        #   ∂(-ξ[j,i]·θ[i])/∂T[i]    = -ξ[j,i] · ∂θ[i]/∂T[i]
        #   ∂(-ξ[j,i]·θ[i])/∂ξ[k,i]  = -ξ[j,i] · ∂θ[i]/∂ξ[k,i]   (k ≠ j)
        #   ∂(-ξ[j,i]·θ[i])/∂ξ[j,i]  = -θ[i] - ξ[j,i] · ∂θ[i]/∂ξ[j,i]
        # Owned here (rather than on `ALEMeshVelocityOperator`) because the
        # correction acts on the ξ row block using the upstream operator's
        # cached strain-rate derivatives.
        _apply_ale_dilation_jacobian!(
            cache.jac,
            cache.ale_mesh_vel_ws,
            ws.state_cache,
            cache.mesh,
            cache.material,
            cache.layout,
        )
    end

    # 4. Geometry block (z columns) — ALE only. AD runs the full residual
    #    against the z subvector. Writes directly into `cache.jac.A` since
    #    the sparsity pattern includes the geometry columns.
    if cache.geom_cache !== nothing
        contribute_to_jacobian!(cache.jac, GeometryOperator(), u, t, cache.geom_cache, def)
    end

    # 5. χ-row block (variable cross-section) — WithChi only. Writes the
    #    χ-row × {T_i, ξ_·,i, z_i, z_{i+1}} cell-local entries to `jac.A`,
    #    plus pushes the χ-row × χ-col rank-1 onto `jac.structured`.
    if cache.chi_ws !== nothing
        contribute_to_jacobian!(
            cache.jac,
            ChiOperator(),
            u,
            t,
            cache.chi_ws,
            ws.state_cache,
            cache.mesh,
            cache.material,
            cache.layout,
            cache.ale_mesh_vel_ws;
            assembly_profile = cache.assembly_profile,
        )
    end

    _emit_jacobian!(J, cache.jac, cache.assembly_profile)
    return J
end

# ALE dilation Jacobian. Mirrors `_apply_ale_dilation_correction!` from
# `Residual/residual.jl`. Writes cell-local entries on the ξ rows of every
# component.
function _apply_ale_dilation_jacobian!(
    jac::StructuredJacobian,
    ale_ws,
    state,
    mesh,
    material,
    layout,
)
    A = jac.A
    smap = jac.smap
    n_cells = layout.n_cells
    NC = n_components(material)
    θ = ale_ws.theta
    dθ_dT = ale_ws.dθ_dT
    dθ_dξ = ale_ws.dθ_dξ

    @inbounds for i = 1:mesh.n_active_cells
        θ_i = θ[i]
        for j = 1:NC
            ξ_ji = state.ξ[j, i]
            ξj_row = n_cells + (j - 1) * n_cells + i

            # T column: ∂(-ξ[j,i]·θ[i])/∂T[i] = -ξ[j,i] · ∂θ[i]/∂T[i]
            fast_sparse_add!(A, smap, ξj_row, i, -ξ_ji * dθ_dT[i])

            # ξ columns
            for k = 1:NC
                ξk_col = n_cells + (k - 1) * n_cells + i
                val = -ξ_ji * dθ_dξ[k, i]
                if k == j
                    val -= θ_i
                end
                fast_sparse_add!(A, smap, ξj_row, ξk_col, val)
            end
        end
    end
    return nothing
end

# Materialize `cache.jac` (sparse `A` + structured couplings) into the
# caller's `J`. Under `assembly_profile == :full` the structured-coupling
# contributions are densified into `J` so downstream consumers (linear
# solvers, `verify_jacobian_at`) see the full operator. Under
# `:ale_reduced` the couplings are *dropped* — that's the profile's
# defining behavior, mirrored on the assembly side by the producers
# skipping their `compose_into!` calls (so `jac.structured` is typically
# empty already), and the masked-compare verifier accounts for it via the
# `is_approximate` trait on `Structured{ApproxSolve}`.
function _emit_jacobian!(
    J::SparseMatrixCSC{Float64,Int},
    jac::StructuredJacobian,
    assembly_profile::Symbol,
)
    copyto!(J.nzval, jac.A.nzval)
    if assembly_profile === :full
        for c in jac.structured
            densify!(J, c)
        end
    end
    return J
end
