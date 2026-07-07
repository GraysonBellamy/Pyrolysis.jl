"""
Adapter — bridge between the user-facing `PyrolysisProblem` and the
unified `ProblemDef` + `Workspace` pair that the residual + Jacobian
operate on.

The canonical front door is [`initialize_state`](@ref), which applies the
problem's initial conditions and returns the ready `(def, ws, u0)` triple:

```julia
def, ws, u0 = initialize_state(problem; use_ale=true, radiation_model=BEER_LAMBERT)
du = similar(u0)
residual!(du, u0, def, 0.0, ws)
```

The underlying primitives stay public as the advanced layer:

- [`to_problem_def`](@ref) — produces the immutable `ProblemDef`.
- [`build_workspace`](@ref) — produces a fresh mutable `Workspace` with
  every cache allocated up front.

Composing them by hand requires `apply_initial_conditions!(problem)`
first — `create_state_vector` packs the mesh as-is (see
`initialize_state`'s docstring for the full contract).
"""

using ..Geometry: Mesh1D, MeshTopology, MeshGeometry, mesh_topology, mesh_geometry
using ..Geometry: GeometryCache
using ..Materials:
    NO_RADIATION,
    SURFACE_ABSORPTION,
    BEER_LAMBERT,
    P1_QUASI_STEADY,
    has_darcy_flow,
    n_components,
    n_reactions
using ..Materials: dry_pyrolysis_gas_indices
using ..Properties: PropertyCache
using ..Discretization: RHSCache, Limiter, VANLEER

# Internal: derive a SimulationMode from a problem and the runtime switches
# `to_problem_def` accepts. Inlined into `to_problem_def`; not part of the
# public surface.
@inline function _simulation_mode_from(
    problem::PyrolysisProblem,
    use_ale::Bool,
    radiation_model,
)
    mesh_t = use_ale ? ALE() : Eulerian()
    rad_t =
        radiation_model === SURFACE_ABSORPTION ? SurfaceAbsorption() :
        radiation_model === BEER_LAMBERT ? BeerLambert() : NoRadiation()
    trans_t = has_darcy_flow(problem.material) ? DarcyFick() : Fickian()
    # χ trait follows the material: identity (no law) → NoChi; otherwise WithChi.
    chi_t = problem.material.lateral_shrinkage_law === nothing ? NoChi() : WithChi()
    return SimulationMode(;
        mesh = mesh_t,
        radiation = rad_t,
        transport = trans_t,
        chi = chi_t,
    )
end

"""
	to_problem_def(problem::PyrolysisProblem; use_ale=false,
	                radiation_model=NO_RADIATION,
	                radiation_intensity=0.0,
	                radiation_kirchhoff=true,
	                energy_form=:standard,
	                use_transpiration_bc=false,
	                use_blowing_correction=false) -> ProblemDef

Lift a user-facing `PyrolysisProblem` into the immutable `ProblemDef`.
The returned object shares the underlying material/BC/initial-condition
*data* with the input but no longer carries any mutable mesh state.

`radiation_kirchhoff` records the surface-transmissivity convention for the
Beer–Lambert injection: `true` (the default) means `radiation_intensity` is
the bare incident flux and the live surface-cell ε_eff is applied at
injection (Kirchhoff α = ε); `false` means the `RadiativeFluxBC`'s explicit
absorptivity is already folded into the intensity. `solve()` derives it from
the top thermal BC.

`use_blowing_correction` (off by default) rescales every `ConvectiveBC` in the
surface energy balance by the Couette-film blowing factor `h → h·B/(e^B − 1)`,
`B = ṁ″·c_p,g/h`, built from the outward gas mass flux at each boundary face.
This is the momentum-side companion of the transpiration enthalpy term and is
fully compatible with the default S_conv-only convention (it modifies the
convective *coefficient*, not the gas enthalpy accounting). At peak burning
(B ≳ 1) it suppresses convective cooling substantially — e.g. B = 4 gives
h_eff/h = 0.075. See `Discretization.blowing_factor`.

!!! warning "Calibrated examples"
    Enabling blowing alone on a case calibrated *without* it (and without
    flame heat-flux feedback) will raise T_s and MLR; blowing and flame
    feedback partially cancel and should be introduced together before any
    recalibration pass.

!!! warning "energy_form = :with_generation_sink requires T_ref-datum heats"
    The gas-generation enthalpy sink is consistent only when every reaction
    heat `h` is referenced to the `T_ref` (≈298.15 K) ambient datum. With
    DSC-style heats quoted at the local reaction temperature (the usual
    interpretation of calibrated / ThermaKin / Gpyro heats), the gas
    products' sensible enthalpy is already inside `h(T)` and the sink
    subtracts it a second time (~0.7–1.2 MJ per kg of volatiles at 900 K).
    Under DSC-style heats the default `:standard` is the exact form. A
    one-time warning is emitted when the sink is selected. See
    `compute_gas_generation_sink!`.
"""
function to_problem_def(
    problem::PyrolysisProblem;
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    radiation_intensity::Real = 0.0,
    radiation_kirchhoff::Bool = true,
    energy_form::Symbol = :standard,
    use_transpiration_bc::Bool = false,
    use_blowing_correction::Bool = false,
)
    # Validate the energy_form switch once at construction. Acceptable values:
    #   :standard               — default Formulation A (drops T·cp·∇·J),
    #                             matches FDS/Gpyro/ThermaKin convention.
    #   :with_generation_sink  — adds the −Σ Δh_g·∇·J_g cell source (the
    #                             gas-generation enthalpy sink: the sensible
    #                             enthalpy carried away by gas at the local
    #                             temperature where it is generated; this is the
    #                             gas-enthalpy-divergence completion term, distinct
    #                             from compression work).
    # DATUM REQUIREMENT: :with_generation_sink is consistent ONLY when every
    # reaction heat h is referenced to the T_ref (≈298.15 K) ambient datum
    # (formation-style h = Σν_j·h_j(T_ref)). With DSC-style heats quoted at the
    # local reaction temperature — the usual interpretation of calibrated /
    # ThermaKin / Gpyro heats — the gas products' sensible enthalpy is already
    # inside h(T), and adding S_gen subtracts it a SECOND time (~0.7–1.2 MJ per
    # kg of volatiles at 900 K). Under DSC-style heats :standard is the exact
    # form. See `compute_gas_generation_sink!` and Guide §6.6.3/§7.5.4.
    # Note: enabling :with_generation_sink is a *partial* move toward
    # conservative form; full conservative form would also need gases back in
    # ρcp_eff (which Formulation A excludes).
    if !(energy_form === :standard || energy_form === :with_generation_sink)
        throw(
            ArgumentError(
                "energy_form must be :standard or :with_generation_sink, got $(repr(energy_form))",
            ),
        )
    end
    if energy_form === :with_generation_sink
        @warn(
            "energy_form = :with_generation_sink requires every reaction heat h " *
            "to be referenced to the T_ref (≈298.15 K) ambient datum. With " *
            "DSC-style heats quoted at the reaction temperature (the usual " *
            "calibrated/ThermaKin/Gpyro interpretation), the gas sensible " *
            "enthalpy is already inside h and this option double-counts it " *
            "(~0.7–1.2 MJ per kg of volatiles at 900 K — comparable to the heat " *
            "of pyrolysis itself). Use the default :standard unless your heats " *
            "are T_ref-datum.",
            maxlog = 1
        )
    end

    # Guard against the gas-convective-enthalpy DOUBLE-COUNT. The volumetric
    # advection source S_conv = -Σ_g cp_g J_g ∂T/∂x is ALWAYS active (both
    # energy_form values include it). The surface transpiration term
    # Σ_g ṁ_g·Δh_g(T_s,T_∞), folded into the boundary energy balance, overlaps
    # S_conv by exactly ∫S_conv — applying both produces a spurious ~5-10 kW/m²
    # sink that under-predicts surface temperature and mass-loss rate. The
    # supported convention (ThermaKin Eq.17 / FDS / Gpyro) is S_conv only with a
    # plain convective+radiative BC. The transpiration path is retained for
    # closed-ledger experiments only and must never be combined with S_conv.
    if use_transpiration_bc
        throw(
            ArgumentError(
                "use_transpiration_bc=true double-counts gas convective enthalpy: " *
                "the volumetric advection source S_conv is always active and the " *
                "surface transpiration enthalpy term overlaps it by ∫S_conv " *
                "(a spurious ~5-10 kW/m² sink). Use the supported convention " *
                "use_transpiration_bc=false (S_conv only, ThermaKin Eq.17). The " *
                "transpiration BC is retained for closed-ledger experiments only " *
                "and must not be combined with the volumetric advection source. " *
                "For mass-flux coupling to the convective BC, use the supported " *
                "use_blowing_correction=true instead (Couette-film h-rescaling; " *
                "no enthalpy double-count).",
            ),
        )
    end
    mode = _simulation_mode_from(problem, use_ale, radiation_model)
    topology = mesh_topology(problem.mesh)
    initial = InitialConditions(problem.T_initial, problem.ξ_initial)
    options = (
        radiation_model = radiation_model,
        radiation_intensity = Float64(radiation_intensity),
        radiation_kirchhoff = radiation_kirchhoff,
        dt_initial = problem.dt_initial,
        output_times = problem.output_times,
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
    )
    return ProblemDef(
        material = problem.material,
        topology = topology,
        bcs = problem.bc_set,
        initial = initial,
        tspan = problem.tspan,
        mode = mode,
        options = options,
    )
end

"""
	build_workspace(problem::PyrolysisProblem, def::ProblemDef;
	                  min_thickness=1e-6,
	                  ale_limiter=VANLEER, scalar_eltype::Type=Float64) -> Workspace

Build a `Workspace` that owns every per-step cache and work array used
by the unified residual. Each cache is sized to the current mesh and
component count.

`scalar_eltype` is the live state's eltype (`Tu` parameter) — `Float64`
for the canonical workspace `solve()` allocates; `ForwardDiff.Dual{...}`
for the AD workspace allocated lazily by an AD Jacobian backend.

The `mesh` carried by the workspace is the same mutable `Mesh1D`
referenced by `problem`. The mesh is **read-only** during residual /
Jacobian assembly; live state lives in `ws.state_cache` and is written
back to the mesh only between accepted timesteps via `StepUpdateCallback`.
"""
function build_workspace(
    problem::PyrolysisProblem,
    def::ProblemDef;
    min_thickness::Real = 1e-6,
    ale_limiter::Limiter = VANLEER,
    scalar_eltype::Type{Tu} = Float64,
) where {Tu}
    mode = def.mode
    NC = n_components(problem.material)
    NR = n_reactions(problem.material)
    mesh = problem.mesh
    geom = mesh_geometry(mesh)

    n_cells = length(mesh.cells)
    n_faces = length(mesh.faces)
    n_nodes = length(mesh.nodes)

    layout = StateLayout(NC, n_cells, n_nodes, mode)

    # ------------------------------------------------------------------
    # Typed live-state mirror — populated at every residual entry.
    # Initialized with the mesh's geometry so kernels see consistent
    # values from the first call (before populate! has run).
    # ------------------------------------------------------------------
    state_cache = StateCache{NC,Tu}(layout, n_faces, mesh.cross_section_area)
    mirror_eulerian_geometry!(state_cache, mesh)

    # ------------------------------------------------------------------
    # Always-present caches (parameterized on Tu)
    # Skip Darcy sub-cache for Fickian-only transport.
    # ------------------------------------------------------------------
    with_darcy = transport_mode(mode) isa DarcyFick
    prop_cache = PropertyCache{Tu}(n_cells, n_faces; with_darcy = with_darcy)
    rhs_cache = RHSCache{Tu}(n_cells, n_faces, NC, NR)

    # ------------------------------------------------------------------
    # Mode-conditional caches — Nothing concrete typed (no Union)
    # ------------------------------------------------------------------
    geom_cache = mesh_mode(mode) isa ALE ? GeometryCache(n_cells, n_faces) : nothing
    rad_cache = nothing  # P1 lives in Pyrolysis.Experimental

    # ------------------------------------------------------------------
    # Physics work arrays (always allocated)
    # ------------------------------------------------------------------
    dT_work = zeros(Tu, n_cells)
    dξ_work = zeros(Tu, NC, n_cells)

    # ------------------------------------------------------------------
    # ALE work arrays — sized for ALE; zero-length for Eulerian so the
    # fields exist (concrete type) but cost nothing.
    # ------------------------------------------------------------------
    if mesh_mode(mode) isa ALE
        w_work = zeros(Tu, n_nodes)
        w_physical = zeros(Tu, n_nodes)
        w_face = zeros(Tu, n_faces)
        θ_work = zeros(Tu, n_cells)
        ale_advection = zeros(Tu, n_cells)
        ale_species_advection = zeros(Tu, NC, n_cells)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = zeros(Tu, n_cells)
    else
        w_work = Tu[]
        w_physical = Tu[]
        w_face = Tu[]
        θ_work = Tu[]
        ale_advection = Tu[]
        ale_species_advection = zeros(Tu, NC, 0)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = Tu[]
    end

    # ------------------------------------------------------------------
    # Pre-cached boundary IDs and dry-pyrolysis gas indices.
    # `dry_pyrolysis_gas_indices` is @generated but its body still
    # allocates per-call (Symbol→String for the water-vapor name filter),
    # so we evaluate once at construction and cache the resulting NTuple.
    # ------------------------------------------------------------------
    dry_gas_indices = dry_pyrolysis_gas_indices(problem.material)
    top_face_ids, top_cell_ids = _boundary_face_cells(mesh, 1)
    bottom_face_ids, bottom_cell_ids = _boundary_face_cells(mesh, 2)

    # ------------------------------------------------------------------
    # Boundary states — tag-keyed view of the mesh's per-id Dict.
    # ------------------------------------------------------------------
    bs_dict = mesh.boundary_states
    bs_top = get(bs_dict, 1, BoundaryState())
    bs_bottom = get(bs_dict, 2, BoundaryState())
    bs_named = (top = bs_top, bottom = bs_bottom)

    # ------------------------------------------------------------------
    # Initial radiation intensity from the problem options.
    # ------------------------------------------------------------------
    I_radiation = Float64(get_option(def, :radiation_intensity, 0.0))

    return Workspace{
        NC,
        typeof(mode),
        Tu,
        typeof(mesh),
        typeof(geom),
        typeof(layout),
        typeof(geom_cache),
        typeof(rad_cache),
        typeof(bs_named),
        typeof(state_cache),
        typeof(dry_gas_indices),
        typeof(prop_cache.darcy),
    }(
        mesh,
        geom,
        layout,
        state_cache,
        prop_cache,
        rhs_cache,
        geom_cache,
        rad_cache,
        dT_work,
        dξ_work,
        w_work,
        w_physical,
        w_face,
        θ_work,
        ale_advection,
        ale_species_advection,
        rates_work,
        heats_work,
        dξ_dt_work,
        S_conv_work,
        dry_gas_indices,
        top_face_ids,
        bottom_face_ids,
        top_cell_ids,
        bottom_cell_ids,
        I_radiation,
        ale_limiter,
        Float64(min_thickness),
        bs_named,
        nothing,  # w_override — verification hook, off by default
    )
end

"""
	initialize_state(problem::PyrolysisProblem;
	                 min_thickness=1e-6,
	                 use_ale=false,
	                 radiation_model=NO_RADIATION,
	                 radiation_intensity=0.0,
	                 radiation_kirchhoff=true,
	                 energy_form=:standard,
	                 use_transpiration_bc=false,
	                 use_blowing_correction=false) -> (def, ws, u0)

Canonical setup path from a `PyrolysisProblem` to a ready-to-integrate
`(ProblemDef, Workspace, u0)` triple. Mirrors `solve()`'s order exactly:

1. [`apply_initial_conditions!`](@ref) — write the problem's T/ξ callables
   into `mesh.cell_states` and populate `mesh.m_dry_solid_initial` (the χ
   progress denominator used under a `lateral_shrinkage_law`).
2. [`to_problem_def`](@ref) — lift to the immutable `ProblemDef` (all
   keyword switches forwarded).
3. [`build_workspace`](@ref) — allocate every residual/Jacobian cache
   (`min_thickness` forwarded).
4. [`create_state_vector`](@ref) — pack the freshly initialized mesh into
   `u0`.

`solve()` itself calls this function; use it whenever you need the raw
residual/Jacobian composition, e.g.

```julia
def, ws, u0 = initialize_state(problem; use_ale = true)
du = similar(u0)
residual!(du, u0, def, 0.0, ws)
```

# Advanced manual composition

The low-level primitives stay public for the genuinely advanced paths —
AD workspaces via [`build_ad_workspace`](@ref), mid-simulation repacking
via [`pack_state!`](@ref), warm starts. If you compose them by hand, the
IC contract is on you: `create_state_vector` packs the mesh **as-is**, so
`apply_initial_conditions!(problem)` must run first or `u0` silently
carries the mesh constructor's placeholder state instead of the problem's
initial conditions (and `m_dry_solid_initial` stays zero, degenerating the
χ block under `WithChi`).
"""
function initialize_state(
    problem::PyrolysisProblem;
    min_thickness::Real = 1e-6,
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    radiation_intensity::Real = 0.0,
    radiation_kirchhoff::Bool = true,
    energy_form::Symbol = :standard,
    use_transpiration_bc::Bool = false,
    use_blowing_correction::Bool = false,
)
    apply_initial_conditions!(problem)
    def = to_problem_def(
        problem;
        use_ale = use_ale,
        radiation_model = radiation_model,
        radiation_intensity = radiation_intensity,
        radiation_kirchhoff = radiation_kirchhoff,
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
    )
    ws = build_workspace(problem, def; min_thickness = min_thickness)
    u0 = create_state_vector(ws.layout, ws.mesh)
    return def, ws, u0
end

"""
	build_ad_workspace(def::ProblemDef, ws::Workspace, ::Type{Tu}) -> Workspace

Build a parallel `Workspace{NC, Mode, Tu}` whose state cache, property
cache, RHS cache, and physics work arrays carry scalar eltype `Tu`
(typically `ForwardDiff.Dual{Tag, Float64, N}`), sharing the canonical
workspace's mesh, layout, boundary IDs, and runtime settings.

Used by AD Jacobian backends to keep a Dual-typed workspace alongside the
Float64 one — the residual then dispatches on `eltype(ws.dT_work) == Tu`
end-to-end without per-call reallocation.

The mesh is **shared** (same identity), reflecting the read-only mesh
contract during residual / Jacobian assembly.
"""
function build_ad_workspace(
    def::ProblemDef,
    ws::Workspace{NC,Mode,Float64},
    ::Type{Tu},
) where {NC,Mode,Tu}
    mode = def.mode
    material = def.material
    NR = n_reactions(material)
    mesh = ws.mesh
    geom = ws.geometry
    layout = ws.layout

    n_cells = layout.n_cells
    n_faces = length(mesh.faces)
    n_nodes = layout.n_nodes

    state_cache = StateCache{NC,Tu}(layout, n_faces, mesh.cross_section_area)
    mirror_eulerian_geometry!(state_cache, mesh)

    # Mirror the parent workspace's Darcy decision so the AD prop_cache has
    # the same shape.
    with_darcy = ws.prop_cache.darcy !== nothing
    prop_cache = PropertyCache{Tu}(n_cells, n_faces; with_darcy = with_darcy)
    rhs_cache = RHSCache{Tu}(n_cells, n_faces, NC, NR)

    geom_cache = ws.geom_cache
    rad_cache = ws.rad_cache

    dT_work = zeros(Tu, n_cells)
    dξ_work = zeros(Tu, NC, n_cells)

    if mesh_mode(mode) isa ALE
        w_work = zeros(Tu, n_nodes)
        w_physical = zeros(Tu, n_nodes)
        w_face = zeros(Tu, n_faces)
        θ_work = zeros(Tu, n_cells)
        ale_advection = zeros(Tu, n_cells)
        ale_species_advection = zeros(Tu, NC, n_cells)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = zeros(Tu, n_cells)
    else
        w_work = Tu[]
        w_physical = Tu[]
        w_face = Tu[]
        θ_work = Tu[]
        ale_advection = Tu[]
        ale_species_advection = zeros(Tu, NC, 0)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = Tu[]
    end

    return Workspace{
        NC,
        typeof(mode),
        Tu,
        typeof(mesh),
        typeof(geom),
        typeof(layout),
        typeof(geom_cache),
        typeof(rad_cache),
        typeof(ws.boundary_states),
        typeof(state_cache),
        typeof(ws.dry_gas_indices),
        typeof(prop_cache.darcy),
    }(
        mesh,
        geom,
        layout,
        state_cache,
        prop_cache,
        rhs_cache,
        geom_cache,
        rad_cache,
        dT_work,
        dξ_work,
        w_work,
        w_physical,
        w_face,
        θ_work,
        ale_advection,
        ale_species_advection,
        rates_work,
        heats_work,
        dξ_dt_work,
        S_conv_work,
        ws.dry_gas_indices,
        copy(ws.top_face_ids),
        copy(ws.bottom_face_ids),
        copy(ws.top_cell_ids),
        copy(ws.bottom_cell_ids),
        ws.I_radiation,
        ws.ale_limiter,
        ws.min_thickness,
        ws.boundary_states,
        ws.w_override,
    )
end
