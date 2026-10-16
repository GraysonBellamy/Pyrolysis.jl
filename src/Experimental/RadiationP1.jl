"""
Pyrolysis.Experimental.RadiationP1 — P1 quasi-steady volumetric radiation.

The P1 model augments the unified residual with a per-cell incident-radiation
field G. G is **not** part of the ODE state vector: it is an algebraic
(quasi-steady) field, re-solved from the current (T, ξ, z) state at every
residual evaluation and stored in the `RadiationCache`. The state vector is
therefore the ordinary unified `StateLayout` (T, ξ, z, χ) — there is no
P1-specific layout.

This module ships:

- `compute_p1_residual!(du, u, def, t, ws, p1)` — the unified residual
  evaluated with `radiation_model = P1_QUASI_STEADY` (so the boundary
  `RadiativeFluxBC`/`RadiativeBC` contributions are zeroed and the incident
  energy is counted exactly once, through the Marshak BC), plus the P1
  volumetric source S_rad = α(G − 4σT⁴) on the energy rows.
- `P1Residual` — the callable wrapper consumed by `OrdinaryDiffEq`.
- `get_radiation_intensity(p::P1Residual, i)` — accessor for the
  quasi-steady G field of the most recent residual evaluation.

P1 is **Experimental**: research-grade, Float64-only (the G Newton solve and
the radiation cache are Float64; state-vector `Dual`s are not supported), and
out of scope for the unified trait dispatch. New code paths should use
`SURFACE_ABSORPTION` or `BEER_LAMBERT`.
"""

using ..Materials: NO_RADIATION, P1_QUASI_STEADY
using ..Problem:
    ProblemDef, Workspace, get_option, chi_mode, WithChi, update_mesh_from_state!
using ..Properties: RadiationCache, resize_radiation_cache!
using ..Residual: residual!

# ============================================================================
# P1 residual — unified residual + quasi-steady G solve + S_rad source
# ============================================================================

@inline function _require_p1_radiation_model(def::ProblemDef)
    model = get_option(def, :radiation_model, NO_RADIATION)
    model === P1_QUASI_STEADY || throw(
        ArgumentError(
            "the P1 residual requires def.options.radiation_model == " *
            "P1_QUASI_STEADY (got $model). Build the def with " *
            "to_problem_def(problem; radiation_model = P1_QUASI_STEADY): the " *
            "core RHS then zeroes the boundary RadiativeFluxBC/RadiativeBC " *
            "contributions, so the incident flux and the surface emission are " *
            "counted exactly once — through the Marshak BC and the volumetric " *
            "4σαT⁴ emission.",
        ),
    )
    return nothing
end

"""
	compute_p1_residual!(du, u, def, t, ws, p1) -> nothing

Unified residual augmented with P1 quasi-steady volumetric radiation.

`u` is an ordinary unified-layout state vector (T, ξ, z when ALE, χ when
WithChi) — G is not a state variable. `p1` provides `rad_cache` and
`radiation_bcs` (a `P1Residual` or any object with those fields).

Steps:
1. Evaluate the unified `residual!`. `def` must carry
   `radiation_model = P1_QUASI_STEADY` (enforced), which zeroes the boundary
   `RadiativeFluxBC`/`RadiativeBC` fluxes and the core volumetric radiation
   source — P1 owns all radiation energy.
2. Mirror `u` into the `Mesh1D` (the P1 kernels are Float64-only and
   read cell temperatures and live geometry from the mesh).
3. Update radiation properties (α, D_rad) and solve the quasi-steady P1
   system for G. This happens on **every** evaluation — G is an algebraic
   field, never frozen.
4. Add S_rad/ρcₚ = α(G − 4σT⁴)/ρcₚ to the energy rows of `du`.
"""
function compute_p1_residual!(
    du::AbstractVector,
    u::AbstractVector,
    def::ProblemDef,
    t::Real,
    ws::Workspace,
    p1,
)
    _require_p1_radiation_model(def)

    mesh = ws.mesh
    rad_cache = p1.rad_cache
    radiation_bcs = p1.radiation_bcs
    n = ws.layout.n_cells

    # 1) Unified residual (conduction, reactions, transport, ALE, χ). Leaves
    #    ws.state_cache / ws.prop_cache consistent with u.
    residual!(du, u, def, t, ws)

    # 2) Mesh sync for the Float64-only P1 kernels. residual! is
    #    mesh-read-only, so the mesh must be refreshed here. The χ-driven
    #    cross-section area must land before the geometry update so cell
    #    volumes are rebuilt against the current A(t).
    if chi_mode(def.mode) isa WithChi
        mesh.cross_section_area = Float64(ws.state_cache.A_section[])
    end
    update_mesh_from_state!(mesh, u, ws.layout)

    # 3) Quasi-steady G solve at the current state.
    update_radiation_properties!(rad_cache, mesh, ws.prop_cache)
    solve_quasi_steady_radiation!(rad_cache, mesh, ws.prop_cache, radiation_bcs, t)

    # 4) Net volumetric radiation source on the energy rows.
    compute_radiation_flux_divergence!(rad_cache, mesh, ws.prop_cache)
    @inbounds for i = 1:n
        du[i] += rad_cache.S_rad[i] / ws.prop_cache.ρcp[i]
    end

    return nothing
end

# ============================================================================
# Callable wrapper for OrdinaryDiffEq
# ============================================================================

"""
	P1Residual{D, W, RC, RBC}

Callable struct holding everything needed to evaluate the P1 residual: the
immutable `def`, the mutable `ws`, the radiation cache (which owns the
quasi-steady G field), and the radiation BCs.

The state vector is the unified layout from `ws.layout` — build `def`, `ws`,
and `u0` with `initialize_state(problem; radiation_model = P1_QUASI_STEADY)`
(which applies the problem's initial conditions before packing; a manual
`create_state_vector` call packs the mesh as-is). Read the G field via
`p.rad_cache.G` or [`get_radiation_intensity`](@ref).

The constructor requires `def.options.radiation_model == P1_QUASI_STEADY`;
this is what keeps the incident radiation single-counted between the Marshak
BC and the thermal boundary conditions.

Construct explicitly in `Pyrolysis.Experimental`; the core `solve` path keeps
P1 out of the unified residual.
"""
struct P1Residual{D,W,RC,RBC}
    def::D
    ws::W
    rad_cache::RC
    radiation_bcs::RBC

    function P1Residual(
        def::ProblemDef,
        ws::Workspace,
        rad_cache::RadiationCache,
        radiation_bcs::RadiationBCSet,
    )
        _require_p1_radiation_model(def)
        return new{typeof(def),typeof(ws),typeof(rad_cache),typeof(radiation_bcs)}(
            def,
            ws,
            rad_cache,
            radiation_bcs,
        )
    end
end

(p::P1Residual)(du, u, _p, t) = compute_p1_residual!(du, u, p.def, t, p.ws, p)

function (p::P1Residual)(u, _p, t)
    du = similar(u)
    p(du, u, _p, t)
    return du
end

"""
	get_radiation_intensity(p::P1Residual, i::Int) -> Float64

Quasi-steady incident radiation G in cell `i` [W/m²], from the most recent
residual evaluation. G is an algebraic field stored in the radiation cache —
it is not part of the ODE state vector.
"""
@inline get_radiation_intensity(p::P1Residual, i::Int) = p.rad_cache.G[i]

# ============================================================================
# Exports
# ============================================================================

export compute_p1_residual!, P1Residual, get_radiation_intensity
