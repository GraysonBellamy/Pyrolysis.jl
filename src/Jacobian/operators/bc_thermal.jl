"""
The `BCThermalOperator` — the boundary-condition thermal contribution extracted from the residual.

Boundary-face thermal flux contribution to `dT[boundary_cell]`. For non-
Dirichlet BCs the surface temperature `T_s` is the solution of the surface
energy balance:

```
q_BC(T_s) = k_eff·(T_s − T_cell)/Δz_half + Σ_g ṁ_g · Δh_g(T_s, T_∞)
```

The transpiration term is present only in the opt-in path (off by default);
the default S_conv-only convention zeroes `ṁ_g` here.

The existing `solve_boundary_temperature!` is **already Dual-aware** (see
[`boundary_temperature.jl`](../Discretization/boundary_temperature.jl)),
so the AD path through Newton is just "call the solver with Dual inputs and
ride the derivative through `q_conducted`."

# Stencil

Cell-local on the boundary cell: `(cells = 0:0,)`. Geometry (`Δz_half`,
`A_face`, `V`) is captured as Float64 constants in the AD closure. The
outward mass flux fed to the surface Newton's *transpiration* argument is
**hard-zeroed**, mirroring the residual: `to_problem_def` rejects
`use_transpiration_bc = true` at construction, and `apply_thermal_bc_flux!`
zeroes `J_boundary_outward` whenever transpiration is off — so the
residual's surface energy balance never sees ṁ_g in the enthalpy term, and
neither must this operator. (Rebuilding ṁ_g from the mass BC here — as an
earlier revision did — makes the Jacobian differentiate a *different*
surface balance than the residual applies whenever a permeable mass BC is
configured.) The boundary mass flux's effect on the energy row flows
through the volumetric `S_conv`/`S_gen` instead, owned by
[`BCMassOperator`](@ref).

# Blowing correction

With `use_blowing = true` the residual's surface Newton rescales every
`ConvectiveBC` by `blowing_factor(ṁ″·c_p,g/h)` built from the **true**
boundary mass flux (`J_gas[:, f]`). This operator then rebuilds that flux
from the stencil state via the shared [`_bc_boundary_face_fluxes`](@ref)
kernel (film law + Darcy venting when `darcy_active`, the same kernel the
residual's mass-BC pass applies) and feeds it to the Newton's `J_blowing`,
so the `∂(h·χ)/∂(T_b, ξ_b)` chain lands exactly. The transpiration argument
stays hard-zeroed — blowing and transpiration are independent terms. Under
ALE the face distance `d_to_face` is a Float64 constant in the closure,
consistent with this operator's existing frozen-geometry convention (the
omitted z-channel degrades only Newton convergence, never the residual).

# Sign convention

For both boundaries, `dT[bdy_cell] += q_conducted · A_face / V_cell / ρcp_eff`.
The conducted flux is

```
q_conducted = k_eff · (T_s − T_cell) / Δz_half
```

with `T_s` from the Newton solve (or `T_bc` for Dirichlet). The signs cancel
out the orientation handling that `_apply_thermal_bc_static!` does on
`q_total[f]` — see derivation in the test reference.

# Notes

- An optional analytical IFT override is deferred.
- `writes = (:T_s,)` is forward-compat metadata; no operator currently
  reads `:T_s`.
"""

using ForwardDiff
using ..Materials: Material, T_REF
using ..Materials: evaluate_param
using ..BoundaryConditions: ambient_temperature
using ..Materials: AbstractThermalBC, RadiationModel, SURFACE_ABSORPTION
using ..BoundaryConditions:
    DirichletThermalBC, BoundaryConditionSet, get_thermal_bc, get_mass_bc
using ..Geometry: BoundaryState
using ..Properties: effective_conductivity, effective_heat_capacity, effective_emissivity
using ..Discretization: solve_boundary_temperature!

"""
	BCThermalOperator <: Operator

Stateless tag.
"""
struct BCThermalOperator <: Operator end

reads(::BCThermalOperator) = ()
writes(::BCThermalOperator) = (:T_s,)
stencil(::BCThermalOperator) = (cells = 0:0,)

"""
	BCThermalADWorkspace{NC, CFG}

Per-boundary AD workspace. Stencil dim = `1 + NC` (T + ξ of the boundary
cell); output dim = 1 (the boundary cell's `dT` contribution).
"""
mutable struct BCThermalADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
end

function build_bc_thermal_workspace(::Material{NC}) where {NC}
    stencil_dim = 1 + NC
    out_dim = 1
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, out_dim)
    J_local = zeros(Float64, out_dim, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return BCThermalADWorkspace{NC,typeof(cfg)}(x_buf, f_buf, J_local, cfg)
end

# ----------------------------------------------------------------------------
# Boundary-cell-local energy contribution (under Duals from contribute_to_jacobian!)
# ----------------------------------------------------------------------------

# Per-boundary signed-flux scalar that goes into `dT[cell_id] += contribution`.
# Generic on T, ξ scalar type — accepts ForwardDiff.Dual.
function _bc_thermal_contribution(
    T_cell::Real,
    ξ_cell::NTuple{NC,Tu},
    thermal_bc::AbstractThermalBC,
    material::Material{NC},
    A_face::Real,
    V_cell::Real,
    Δz_half::Real,
    t::Real,
    T_ref::Real,
    radiation_model::RadiationModel,
    mass_bc = nothing,
    d_to_face::Real = Δz_half,
    tag::Val = Val(:top),
    darcy_active::Bool = false,
    use_blowing::Bool = false,
) where {NC,Tu}

    # Cell-local properties (Dual-typed when ξ_cell, T_cell are Dual).
    k_eff = effective_conductivity(ξ_cell, T_cell, material)
    ρcp = effective_heat_capacity(ξ_cell, T_cell, material)
    ε_eff = effective_emissivity(ξ_cell, T_cell, material)

    # Outward gas flux at this boundary face — hard-zeroed, mirroring the
    # residual: transpiration is rejected at construction
    # (`to_problem_def(use_transpiration_bc = true)` throws) and
    # `apply_thermal_bc_flux!` zeroes `J_boundary_outward` when it is off,
    # so the residual's surface Newton never sees ṁ_g in the enthalpy term.
    # The Jacobian must differentiate the same surface balance the residual
    # applies.
    J_boundary_outward = ntuple(_ -> zero(Tu), Val(NC))

    # Blowing reads the TRUE boundary mass flux: rebuild it from the stencil
    # state with the same kernel the residual's mass-BC pass applies.
    # `_bc_boundary_face_fluxes` returns F positive INTO the domain, so the
    # outward flux is −F at either face.
    J_blowing = if use_blowing
        F = _bc_boundary_face_fluxes(
            T_cell,
            ξ_cell,
            mass_bc,
            material,
            d_to_face,
            t,
            tag,
            darcy_active,
        )
        ntuple(j -> -F[j], Val(NC))
    else
        J_boundary_outward
    end

    # Solve for T_s (Dirichlet branch is closed-form; non-Dirichlet uses
    # the Newton solver, which is already Dual-aware).
    q_conducted = if thermal_bc isa DirichletThermalBC
        T_bc = evaluate_param(thermal_bc.T, t)
        k_eff * (T_bc - T_cell) / Δz_half
    else
        # Use a fresh BoundaryState (Float64 storage; the solver strips Duals
        # for storage but the returned tuple carries Tu).
        state = BoundaryState(
            Float64(ForwardDiff.value(T_cell)),
            0.0,
            0.0,
            0.0,
            Float64(ForwardDiff.value(T_cell)),
            Float64(ForwardDiff.value(k_eff)),
            Float64(ForwardDiff.value(Δz_half)),
            false,
            0,
        )
        result = solve_boundary_temperature!(
            state,
            T_cell,
            k_eff,
            Δz_half,
            thermal_bc,
            material,
            J_boundary_outward,
            T_ref,
            t,
            ε_eff;
            radiation_model = radiation_model,
            use_blowing = use_blowing,
            J_blowing = J_blowing,
        )
        result.q_conducted
    end

    # `dT[bdy] += q_conducted · A_face / V_cell / ρcp_eff` (sign derivation in
    # this file's docstring — the q-orientation cancels with divergence sign
    # at both top and bottom).
    return q_conducted * A_face / V_cell / ρcp
end

# Dispatcher: extract per-boundary scalars and call the typed kernel.
@inline function _bc_thermal_geometry(state, mesh, boundary_id::Int)
    if boundary_id == 1
        cell_id = mesh.n_active_cells
        face_id = mesh.cells[cell_id].face_ids[2]
    else
        cell_id = 1
        face_id = mesh.cells[cell_id].face_ids[1]
    end
    A_face = state.A_face[face_id]
    V_cell = state.volume[cell_id]
    Δz_half = V_cell / A_face / 2
    face_pos = state.face_positions[face_id]
    cell_ctr = state.cell_centers[cell_id]
    d_to_face = abs(face_pos - cell_ctr)
    return (
        cell_id = cell_id,
        face_id = face_id,
        A_face = Float64(A_face),
        V_cell = Float64(V_cell),
        Δz_half = Float64(Δz_half),
        d_to_face = Float64(d_to_face),
    )
end

"""
	apply!(dT, ::BCThermalOperator, state, mesh, material, bc_set, t,
	       radiation_model, darcy_active = false, use_blowing = false) -> nothing

Accumulate `dT[boundary_cell] += q_conducted · A_face / V_cell / ρcp_eff` for
both top (boundary_id = 1) and bottom (boundary_id = 2) boundaries.
"""
function apply!(
    dT::AbstractVector,
    ::BCThermalOperator,
    state,
    mesh,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    t::Real,
    radiation_model::RadiationModel,
    darcy_active::Bool = false,
    use_blowing::Bool = false,
) where {NC}

    for (boundary_id, tag) in ((1, Val(:top)), (2, Val(:bottom)))
        geom = _bc_thermal_geometry(state, mesh, boundary_id)
        thermal_bc =
            boundary_id == 1 ? get_thermal_bc(bc_set, Val(:top)) :
            get_thermal_bc(bc_set, Val(:bottom))
        mass_bc =
            boundary_id == 1 ? get_mass_bc(bc_set, Val(:top)) :
            get_mass_bc(bc_set, Val(:bottom))
        T_ref = ambient_temperature(thermal_bc, t, Float64(T_REF))

        T_cell = state.T[geom.cell_id]
        ξ_cell = ntuple(j -> state.ξ[j, geom.cell_id], Val(NC))

        dT_contrib = _bc_thermal_contribution(
            T_cell,
            ξ_cell,
            thermal_bc,
            material,
            geom.A_face,
            geom.V_cell,
            geom.Δz_half,
            t,
            T_ref,
            radiation_model,
            mass_bc,
            geom.d_to_face,
            tag,
            darcy_active,
            use_blowing,
        )

        dT[geom.cell_id] += dT_contrib
    end
    return nothing
end

# Pack/unpack helpers for the AD closure: input vector is `[T_cell, ξ_1, ..., ξ_NC]`.
@inline function _bct_pack!(
    x::AbstractVector,
    u::AbstractVector,
    n_cells::Int,
    NC::Int,
    cell_id::Int,
)
    @inbounds begin
        x[1] = u[cell_id]
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            x[1+j] = u[ξ_off+cell_id]
        end
    end
    return nothing
end

# Takes `Val(NC)` (not a runtime `Int`) so the `ntuple` constant-folds: the
# AD closure captures only the `Val`, keeping the unpacked ξ tuple
# stack-allocated inside the ForwardDiff sweep.
@inline function _bct_unpack(x::AbstractVector, ::Val{NC}) where {NC}
    @inbounds T_cell = x[1]
    ξ_cell = ntuple(j -> @inbounds(x[1+j]), Val(NC))
    return T_cell, ξ_cell
end

"""
	contribute_to_jacobian!(jac, op::BCThermalOperator, u, ctx, t, ws, state, mesh,
	                        material, layout, bc_set, radiation_model,
	                        darcy_active = false, use_blowing = false) -> nothing

Per-boundary AD path. Each boundary face:

1. Pack `(T_cell, ξ_cell)` of the boundary cell into `ws.x_buf`.
2. `ForwardDiff.jacobian!` over the `_bc_thermal_contribution` closure.
3. Scatter the resulting `1 × (1 + NC)` row into `jac.A` at row =
   T-row[cell_id], cols = (T-col, ξ-cols of cell_id).
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::BCThermalOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::BCThermalADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    bc_set::BoundaryConditionSet,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    darcy_active::Bool = false,
    use_blowing::Bool = false,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap

    # Two statically-typed calls instead of a loop: a branch-merged
    # `top-or-bottom` BC lookup is a Union, which heap-allocates the AD
    # closure and dynamic-dispatches `jacobian!`. With per-branch lookups
    # the barrier specializes on each boundary's concrete BC types.
    geom_top = _bc_thermal_geometry(state, mesh, 1)
    _bct_assemble_boundary!(
        A,
        smap,
        ws,
        u,
        n_cells,
        geom_top,
        get_thermal_bc(bc_set, Val(:top)),
        get_mass_bc(bc_set, Val(:top)),
        Val(:top),
        material,
        t,
        radiation_model,
        darcy_active,
        use_blowing,
    )
    geom_bot = _bc_thermal_geometry(state, mesh, 2)
    _bct_assemble_boundary!(
        A,
        smap,
        ws,
        u,
        n_cells,
        geom_bot,
        get_thermal_bc(bc_set, Val(:bottom)),
        get_mass_bc(bc_set, Val(:bottom)),
        Val(:bottom),
        material,
        t,
        radiation_model,
        darcy_active,
        use_blowing,
    )
    return nothing
end

function _bct_assemble_boundary!(
    A,
    smap,
    ws::BCThermalADWorkspace{NC},
    u::AbstractVector,
    n_cells::Int,
    geom,
    thermal_bc,
    mass_bc,
    tag::Val,
    material::Material{NC},
    t::Real,
    radiation_model::RadiationModel,
    darcy_active::Bool,
    use_blowing::Bool,
) where {NC}
    T_ref = ambient_temperature(thermal_bc, t, Float64(T_REF))

    f! =
        let geom = geom,
            thermal_bc = thermal_bc,
            mass_bc = mass_bc,
            material = material,
            t = t,
            T_ref = T_ref,
            radiation_model = radiation_model,
            val_nc = Val(NC),
            tag = tag,
            darcy_active = darcy_active,
            use_blowing = use_blowing

            (out, x_inner) -> begin
                T_cell, ξ_cell = _bct_unpack(x_inner, val_nc)
                out[1] = _bc_thermal_contribution(
                    T_cell,
                    ξ_cell,
                    thermal_bc,
                    material,
                    geom.A_face,
                    geom.V_cell,
                    geom.Δz_half,
                    t,
                    T_ref,
                    radiation_model,
                    mass_bc,
                    geom.d_to_face,
                    tag,
                    darcy_active,
                    use_blowing,
                )
                return nothing
            end
        end

    _bct_pack!(ws.x_buf, u, n_cells, NC, geom.cell_id)
    ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

    # Scatter row 1 (T-row of boundary cell) — input col 1 = T_cell, cols 2..1+NC = ξ_cell[j].
    gr = geom.cell_id
    @inbounds fast_sparse_add!(A, smap, gr, geom.cell_id, ws.J_local[1, 1])
    @inbounds for j = 1:NC
        ξ_off = n_cells + (j - 1) * n_cells
        fast_sparse_add!(A, smap, gr, ξ_off + geom.cell_id, ws.J_local[1, 1+j])
    end
    return nothing
end
