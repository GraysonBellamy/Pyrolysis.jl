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
    # Only the ALE residual evolves χ and the cross-section it drives.
    law = problem.material.lateral_shrinkage_law
    if law !== nothing && !use_ale
        throw(
            ArgumentError(
                "material $(problem.material.name) has a lateral_shrinkage_law, which " *
                "needs the ALE mesh: on the fixed mesh neither the pyrolysis progress " *
                "χ nor the cross-section area evolves. Solve with use_ale = true, or " *
                "build the material without a lateral_shrinkage_law.",
            ),
        )
    end
    chi_t = law === nothing ? NoChi() : WithChi()
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
	                radiation_kirchhoff=true) -> ProblemDef

Lift a user-facing `PyrolysisProblem` into the immutable `ProblemDef`.
The returned object shares the underlying material/BC/initial-condition
*data* with the input but no longer carries any mutable mesh state.

`radiation_kirchhoff` records the surface-transmissivity convention for the
Beer–Lambert injection: `true` (the default) means `radiation_intensity` is
the bare incident flux and the live surface-cell ε_eff is applied at
injection (Kirchhoff α = ε); `false` means the `RadiativeFluxBC`'s explicit
absorptivity is already folded into the intensity. `solve()` derives it from
the top thermal BC.

The effect of outgassing on convective exchange is a property of each
`ConvectiveBC` (its `blowing` model), not a solve option: the surface balance
feeds every convective member the outward gas mass flux at its face and the
model returns the coefficient ratio. This modifies the convective
*coefficient*, not the gas enthalpy accounting, which stays with the
volumetric source `S_conv`. See `Discretization.blowing_state`.

!!! warning "Calibrated examples"
    Enabling blowing alone on a case calibrated *without* it (and without
    flame heat-flux feedback) will raise T_s and MLR; blowing and flame
    feedback partially cancel and should be introduced together before any
    recalibration pass.
"""
function to_problem_def(
    problem::PyrolysisProblem;
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    radiation_intensity::Real = 0.0,
    radiation_kirchhoff::Bool = true,
)
    mode = _simulation_mode_from(problem, use_ale, radiation_model)
    topology = mesh_topology(problem.mesh)
    initial = InitialConditions{Any,Any}(problem.T_initial, problem.ξ_initial)
    options = (
        radiation_model = radiation_model,
        radiation_intensity = Float64(radiation_intensity),
        radiation_kirchhoff = radiation_kirchhoff,
        dt_initial = problem.dt_initial,
        output_times = problem.output_times,
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
        θ_work = zeros(Tu, n_cells)
        ale_advection = zeros(Tu, n_cells)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = zeros(Tu, n_cells)
    else
        w_work = Tu[]
        w_physical = Tu[]
        θ_work = Tu[]
        ale_advection = Tu[]
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
        θ_work,
        ale_advection,
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
	                 radiation_kirchhoff=true) -> (def, ws, u0)

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
)
    apply_initial_conditions!(problem)
    def = to_problem_def(
        problem;
        use_ale = use_ale,
        radiation_model = radiation_model,
        radiation_intensity = radiation_intensity,
        radiation_kirchhoff = radiation_kirchhoff,
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
        θ_work = zeros(Tu, n_cells)
        ale_advection = zeros(Tu, n_cells)
        rates_work = zeros(Tu, NR)
        heats_work = zeros(Tu, NR)
        dξ_dt_work = zeros(Tu, NC)
        S_conv_work = zeros(Tu, n_cells)
    else
        w_work = Tu[]
        w_physical = Tu[]
        θ_work = Tu[]
        ale_advection = Tu[]
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
        θ_work,
        ale_advection,
        rates_work,
        heats_work,
        dξ_dt_work,
        S_conv_work,
        ws.dry_gas_indices,
        # Shared with the canonical workspace: a park event moves the exposed
        # boundary face, and both workspaces must see the new face id.
        ws.top_face_ids,
        ws.bottom_face_ids,
        ws.top_cell_ids,
        ws.bottom_cell_ids,
        ws.I_radiation,
        ws.ale_limiter,
        ws.min_thickness,
        ws.boundary_states,
        ws.w_override,
    )
end
