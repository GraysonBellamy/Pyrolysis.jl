"""
The `BCMassOperator` — the boundary-condition mass contribution extracted from the residual.

Owns the species rows of the boundary-face gas mass flux. The transport
operators ([`GasTransportOperator`](@ref) / [`DarcyTransportOperator`](@ref))
zero their boundary faces; this operator supplies what they drop on those
rows. The boundary cell's **`S_conv` energy row is NOT owned here** — its
gradient runs through the surface-Newton `T_s` and multiplies the full
face-averaged `J̄_g`, so the whole row has a single owner,
[`BCConvectionOperator`](@ref) (which reuses this file's
[`_bc_boundary_face_fluxes`](@ref) kernel for the boundary-face flux).

Species rows (per gas component `g`, at the boundary cell):

```
dξ[g, bdy] += F_g · A_face / V_cell
F_g = mass_bc_flux(mass_bc, ξ_bdy[g], λρ, ρ, d, t, g)     (film law)
    − v_out · ξ_bdy[g] / φ_bdy                            (Darcy mode only)
```

`F_g` is the boundary-face mass flux per area, **positive into the
domain** (negative when gas escapes). The advective venting term mirrors
the pressure-BC boundary flux that `compute_gas_fluxes_darcy_fick!`
writes: `v_out` is the outward boundary Darcy speed from
[`boundary_darcy_velocity`](@ref) (clamped at 0 — inflow advects only the
zero-concentration ambient) and ξ/φ is the intrinsic gas density.

# Stencil

Cell-local on the boundary cell: the mass-BC closure and the venting chain
(P, κ, μ, φ) depend only on `(T_bdy, ξ_bdy)`.
"""

using ForwardDiff
using ..Materials: Material, is_gas
using ..BoundaryConditions: BoundaryConditionSet
using ..BoundaryConditions: get_mass_bc
using ..Properties:
    effective_gas_transfer,
    mixture_density,
    porosity,
    effective_permeability,
    effective_gas_viscosity,
    mixture_gas_viscosity
using ..Physics: ideal_gas_pressure
using ..Discretization: mass_bc_flux, boundary_darcy_velocity

"""
	BCMassOperator <: Operator

Stateless tag.
"""
struct BCMassOperator <: Operator end

"""
	BCMassADWorkspace{NC, CFG}

Per-boundary AD workspace. Stencil dim = `1 + NC` (`T_bdy`, `ξ_bdy[1:NC]`);
output dim = `NC` (one entry per ξ-component; non-gas entries are identically
zero).
"""
mutable struct BCMassADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
end

function build_bc_mass_workspace(::Material{NC}) where {NC}
    stencil_dim = 1 + NC
    out_dim = NC
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, out_dim)
    J_local = zeros(Float64, out_dim, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return BCMassADWorkspace{NC,typeof(cfg)}(x_buf, f_buf, J_local, cfg)
end

"""
	_bc_boundary_face_fluxes(T_cell, ξ_cell, mass_bc, material,
	                         d_to_face, t, ::Val{tag}, darcy_active)
	    -> NTuple{NC}

Per-component boundary-face gas mass flux per area, **positive into the
domain**: the diffusive film law plus (when `darcy_active`) the pressure-BC
advective venting that mirrors the boundary-face branch of
`compute_gas_fluxes_darcy_fick!`. Non-gas components are zero. Shared by
[`BCMassOperator`](@ref) (species rows) and
[`BCConvectionOperator`](@ref) (the `J̄_g` face average of the boundary
cell's S_conv energy row); generic over Dual scalar types.
"""
function _bc_boundary_face_fluxes(
    T_cell::Real,
    ξ_cell::NTuple{NC,Tu},
    mass_bc,
    material::Material{NC},
    d_to_face::Real,
    t::Real,
    ::Val{tag},
    darcy_active::Bool,
) where {NC,Tu,tag}

    λ_eff = effective_gas_transfer(ξ_cell, T_cell, material)
    ρ_eff = mixture_density(ξ_cell, T_cell, material)
    λρ = λ_eff * ρ_eff

    Tv = promote_type(typeof(T_cell), Tu)
    F = ntuple(Val(NC)) do j
        if is_gas(material.components[j])
            Tv(mass_bc_flux(mass_bc, ξ_cell[j], λρ, ρ_eff, d_to_face, t, j))
        else
            zero(Tv)
        end
    end

    if darcy_active
        # Advective venting: outward Darcy speed advecting the intrinsic
        # density ξ/φ out of the boundary cell (inflow advects only the
        # zero-concentration ambient, hence the one-sided clamp).
        φ = porosity(ξ_cell, T_cell, material)
        P_int = ideal_gas_pressure(ξ_cell, T_cell, φ, material)
        κ_int = effective_permeability(ξ_cell, T_cell, material)
        μ_int = mixture_gas_viscosity(ξ_cell, T_cell, material)
        v_z = boundary_darcy_velocity(mass_bc, P_int, κ_int, μ_int, d_to_face, Val(tag))
        v_out = tag === :top ? v_z : -v_z
        if v_out > 0
            coef = v_out / max(φ, 1e-12)
            # Return a fresh tuple instead of rebinding `F` — a self-referential
            # rebind would box the variable (type-unstable, and inference can't
            # prove the boxed read defined).
            return ntuple(Val(NC)) do j
                if is_gas(material.components[j])
                    F[j] - coef * ξ_cell[j]
                else
                    F[j]
                end
            end
        end
    end

    return F
end

# Generic per-cell contribution: species rows in `out[j]`. `tag` is the
# boundary orientation (:top/:bottom); `darcy_active` enables the pressure-BC
# advective venting term, which exists only on the Darcy-Fick transport path
# (the Fickian residual never consumes a boundary Darcy velocity).
function _bc_mass_contribution!(
    out::AbstractVector,
    T_cell::Real,
    ξ_cell::NTuple{NC,Tu},
    mass_bc,
    material::Material{NC},
    A_face::Real,
    V_cell::Real,
    d_to_face::Real,
    t::Real,
    tag::Val,
    darcy_active::Bool,
) where {NC,Tu}

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

    # Species rows: scale the per-area flux to the cell update. The boundary
    # cell's energy row is owned whole by BCConvectionOperator.
    AV = A_face / V_cell
    @inbounds for j = 1:NC
        out[j] = F[j] * AV
    end
    return nothing
end

@inline function _bcm_geometry(state, mesh, boundary_id::Int)
    n = mesh.n_active_cells
    if boundary_id == 1
        cell_id = n
        face_id = mesh.cells[cell_id].face_ids[2]
    else
        cell_id = 1
        face_id = mesh.cells[cell_id].face_ids[1]
    end
    A_face = state.A_face[face_id]
    V_cell = state.volume[cell_id]
    face_pos = state.face_positions[face_id]
    cell_ctr = state.cell_centers[cell_id]
    d_to_face = abs(face_pos - cell_ctr)
    return (
        cell_id = cell_id,
        face_id = face_id,
        A_face = Float64(A_face),
        V_cell = Float64(V_cell),
        d_to_face = Float64(d_to_face),
    )
end

"""
	apply!(dξ, ::BCMassOperator, state, mesh, material, bc_set, t,
	       darcy_active) -> nothing

Accumulate the boundary mass-flux contribution — the diffusive film law plus,
when `darcy_active`, the pressure-BC advective venting — into
`dξ[g, bdy_cell]` for both boundaries.
"""
function apply!(
    dξ::AbstractMatrix,
    ::BCMassOperator,
    state,
    mesh,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    t::Real,
    darcy_active::Bool,
) where {NC}

    out = zeros(eltype(dξ), NC)
    for boundary_id in (1, 2)
        geom = _bcm_geometry(state, mesh, boundary_id)
        mass_bc =
            boundary_id == 1 ? get_mass_bc(bc_set, Val(:top)) :
            get_mass_bc(bc_set, Val(:bottom))
        tag = boundary_id == 1 ? Val(:top) : Val(:bottom)

        T_cell = state.T[geom.cell_id]
        ξ_cell = ntuple(j -> state.ξ[j, geom.cell_id], Val(NC))

        fill!(out, zero(eltype(dξ)))
        _bc_mass_contribution!(
            out,
            T_cell,
            ξ_cell,
            mass_bc,
            material,
            geom.A_face,
            geom.V_cell,
            geom.d_to_face,
            t,
            tag,
            darcy_active,
        )
        @inbounds for j = 1:NC
            dξ[j, geom.cell_id] += out[j]
        end
    end
    return nothing
end

@inline function _bcm_pack!(
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

# `Val(NC)` so the `ntuple` constant-folds inside the AD closure (see
# `_bct_unpack`).
@inline function _bcm_unpack(x::AbstractVector, ::Val{NC}) where {NC}
    @inbounds T_cell = x[1]
    ξ_cell = ntuple(j -> @inbounds(x[1+j]), Val(NC))
    return T_cell, ξ_cell
end

"""
	contribute_to_jacobian!(jac, op::BCMassOperator, u, t, ws, state, mesh,
	                        material, layout, bc_set, darcy_active) -> nothing

Per-boundary AD path. Each boundary face writes an `NC × (1 + NC)` block of
ξ-rows over the boundary cell's `[T-col, ξ-cols]`. When `darcy_active`, the
block includes the pressure-BC advective venting chain (through P, κ, μ, φ
of the boundary cell).
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::BCMassOperator,
    u::AbstractVector,
    t::Real,
    ws::BCMassADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    bc_set::BoundaryConditionSet,
    darcy_active::Bool,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap

    # Per-branch BC lookups + barrier (see `_bct_assemble_boundary!`): keeps
    # the AD closure concretely typed per boundary.
    geom_top = _bcm_geometry(state, mesh, 1)
    _bcm_assemble_boundary!(
        A,
        smap,
        ws,
        u,
        n_cells,
        geom_top,
        get_mass_bc(bc_set, Val(:top)),
        Val(:top),
        material,
        t,
        darcy_active,
    )
    geom_bot = _bcm_geometry(state, mesh, 2)
    _bcm_assemble_boundary!(
        A,
        smap,
        ws,
        u,
        n_cells,
        geom_bot,
        get_mass_bc(bc_set, Val(:bottom)),
        Val(:bottom),
        material,
        t,
        darcy_active,
    )
    return nothing
end

function _bcm_assemble_boundary!(
    A,
    smap,
    ws::BCMassADWorkspace{NC},
    u::AbstractVector,
    n_cells::Int,
    geom,
    mass_bc,
    tag::Val,
    material::Material{NC},
    t::Real,
    darcy_active::Bool,
) where {NC}
    f! =
        let geom = geom,
            mass_bc = mass_bc,
            material = material,
            t = t,
            val_nc = Val(NC),
            tag = tag,
            darcy_active = darcy_active

            (out, x_inner) -> begin
                T_cell, ξ_cell = _bcm_unpack(x_inner, val_nc)
                _bc_mass_contribution!(
                    out,
                    T_cell,
                    ξ_cell,
                    mass_bc,
                    material,
                    geom.A_face,
                    geom.V_cell,
                    geom.d_to_face,
                    t,
                    tag,
                    darcy_active,
                )
                return nothing
            end
        end

    _bcm_pack!(ws.x_buf, u, n_cells, NC, geom.cell_id)
    ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

    # Scatter the NC ξ-rows: cols 1 = T_cell, 2..1+NC = ξ_cell[j].
    @inbounds for r = 1:NC
        gr = n_cells + (r - 1) * n_cells + geom.cell_id
        fast_sparse_add!(A, smap, gr, geom.cell_id, ws.J_local[r, 1])
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            fast_sparse_add!(A, smap, gr, ξ_off + geom.cell_id, ws.J_local[r, 1+j])
        end
    end
    return nothing
end
