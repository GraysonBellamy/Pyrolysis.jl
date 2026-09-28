"""
The `BCConvectionOperator` — single owner of the boundary cell's `S_conv` energy row.

The residual's boundary-cell gas convective source is

```
S_conv,b = −Σ_g cp_g(T_b) · J̄_g,b · (∂T/∂x)_b
J̄_g,b   = (J_g[f_bc] + J_g[f_int]) / 2
(∂T/∂x)_b = (d_nb·s_face + Δz_half·s_int) / (d_nb + Δz_half)
```

where `s_face = ±(T_s − T_b)/Δz_half` is the boundary-face slope at the
surface temperature solved by the thermal-BC Newton (the residual recovers
it AD-cleanly as `−q_total[f]/k_b`; here it is formed directly from the
Dual-aware Newton — algebraically identical because the residual's
`q_total[f] = ∓k_b(T_s − T_b)/Δz_half` makes the `k_b` cancel exactly) and
`s_int` is the interior one-sided slope. The same sign-consistent
expression holds at both boundaries with `s_face`/`s_int` written in the
+z orientation.

This row is a **product** of pieces owned by three different operators —
the interior-face flux kernel (transport operators), the boundary-face BC
flux (mass-BC kernel), and the `T_s` implicit chain (thermal-BC Newton) —
so a split assembly would need `J̄ × ∂(∂T/∂x)/∂x` cross-products spanning
owners. Instead this operator owns the whole row: the transport operators
zero their boundary-cell `S_conv` (they keep the boundary cell's ξ-rows
and the linear `S_gen` interior slice), and [`BCMassOperator`](@ref) keeps
only the linear `ΔS_gen` boundary-face slice on the energy row.

# Kernels

- Boundary-face flux: [`_bc_boundary_face_fluxes`](@ref) (shared with
  `BCMassOperator`) — film law plus pressure-BC venting when
  `darcy_active`, oriented to +z with `σ = ∓1`.
- Interior-face flux: mode-appropriate — Fickian [`gas_flux`](@ref) with
  the distance-weighted harmonic `λ_face`, or Darcy-Fick
  [`gas_flux_darcy_fick`](@ref) with the cell-local `(P, κ, μ, φ)`
  reconstruction and harmonic face mobility, byte-matching the transport
  operators' kernels.
- `T_s`: the Dual-aware `solve_boundary_temperature!` Newton with the
  *transpiration* mass flux hard-zeroed (the residual's surface balance
  never sees ṁ_g in the enthalpy term — transpiration is rejected at
  construction). With `use_blowing` the Newton's `J_blowing` is fed the
  true outward flux −F from the shared boundary-face kernel, so the
  convective `h·blowing_factor(ṁ″c_p,g/h)` chain is differentiated
  exactly. Dirichlet BCs take the closed form.

# Stencil

`(cells = -1:1,)`: the energy row of each boundary cell, columns
`(T, ξ)` of the boundary cell and its single interior neighbor — all
inside the physics tridiagonal. A single-cell mesh contributes nothing
(the residual sets `S_conv = 0` on isolated cells).

This file relies on sibling-operator imports (all operator files share
the `Pyrolysis.Jacobian` module): `gas_flux`/`harmonic_mean_weighted`
(gas_transport.jl), `gas_flux_darcy_fick`/`porosity`/`ideal_gas_pressure`/
`effective_permeability`/`mixture_gas_viscosity`/`_dt_mobility_face`
(darcy_transport.jl), `_bc_boundary_face_fluxes`/`get_mass_bc`
(bc_mass.jl), and the thermal-BC Newton plumbing (bc_thermal.jl).
"""

"""
	BCConvectionOperator <: Operator

Stateless tag.
"""
struct BCConvectionOperator <: Operator end

reads(::BCConvectionOperator) = (:T_s,)
writes(::BCConvectionOperator) = ()
stencil(::BCConvectionOperator) = (cells = -1:1,)

"""
	BCConvectionADWorkspace{NC, CFG}

Per-boundary AD workspace. Stencil dim = `2·(1 + NC)`
(`T_bdy, ξ_bdy[1:NC], T_nb, ξ_nb[1:NC]`); output dim = 1 (the boundary
cell's `S_conv/ρcp_eff` energy contribution).
"""
mutable struct BCConvectionADWorkspace{NC,CFG}
    x_buf::Vector{Float64}
    f_buf::Vector{Float64}
    J_local::Matrix{Float64}
    cfg::CFG
end

function build_bc_convection_workspace(::Material{NC}) where {NC}
    stencil_dim = 2 * (1 + NC)
    out_dim = 1
    x_buf = zeros(Float64, stencil_dim)
    f_buf = zeros(Float64, out_dim)
    J_local = zeros(Float64, out_dim, stencil_dim)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x_buf, ForwardDiff.Chunk{stencil_dim}())
    return BCConvectionADWorkspace{NC,typeof(cfg)}(x_buf, f_buf, J_local, cfg)
end

# Per-boundary geometry slice (Float64 constants captured in the AD closure).
# `d1_int`/`d2_int` are the interior face's harmonic-weighting distances in
# +z orientation (d1 = lower-z cell to face, d2 = face to upper-z cell),
# matching `_gt_geometry`/`_dt_geometry`.
@inline function _bcc_geometry(state, mesh, boundary_id::Int)
    n = mesh.n_active_cells
    if boundary_id == 1
        cell_id = n
        face_bc = mesh.cells[cell_id].face_ids[2]
        face_int = mesh.cells[cell_id].face_ids[1]
        nb_id = cell_id - 1
    else
        cell_id = 1
        face_bc = mesh.cells[cell_id].face_ids[1]
        face_int = mesh.cells[cell_id].face_ids[2]
        nb_id = 2
    end
    has_nb = 1 <= nb_id <= n && nb_id != cell_id
    A_bc = state.A_face[face_bc]
    V_cell = state.volume[cell_id]
    Δz_half = V_cell / A_bc / 2
    d_to_face = abs(state.face_positions[face_bc] - state.cell_centers[cell_id])
    if has_nb
        d_nb = abs(state.cell_centers[nb_id] - state.cell_centers[cell_id])
        lower_id = boundary_id == 1 ? nb_id : cell_id
        d1_int = abs(state.face_positions[face_int] - state.cell_centers[lower_id])
        d2_int = d_nb - d1_int
    else
        d_nb = 1.0
        d1_int = 0.5
        d2_int = 0.5
    end
    return (
        cell_id = cell_id,
        nb_id = has_nb ? nb_id : 0,
        has_nb = has_nb,
        Δz_half = Float64(Δz_half),
        d_to_face = Float64(d_to_face),
        d_nb = Float64(d_nb),
        d1_int = Float64(d1_int),
        d2_int = Float64(d2_int),
    )
end

# Interior-face per-component flux in +z orientation. `(T_lo, ξ_lo)` is the
# lower-z side of the face, `(T_hi, ξ_hi)` the upper-z side. Mirrors the
# transport operators' face kernels (`_gt_face_fluxes` / `_dt_face_fluxes`
# reconstruction) so the residual's interior flux is reproduced exactly.
function _bcc_interior_face_fluxes(
    T_lo::Real,
    T_hi::Real,
    ξ_lo::NTuple{NC,Tu},
    ξ_hi::NTuple{NC,Tu},
    d1::Real,
    d2::Real,
    d_LR::Real,
    material::Material{NC},
    darcy_active::Bool,
) where {NC,Tu}
    λ_lo = effective_gas_transfer(ξ_lo, T_lo, material)
    λ_hi = effective_gas_transfer(ξ_hi, T_hi, material)
    λ_face = harmonic_mean_weighted(λ_lo, λ_hi, d1, d2)

    if darcy_active
        φ_lo = porosity(ξ_lo, T_lo, material)
        φ_hi = porosity(ξ_hi, T_hi, material)
        P_lo = ideal_gas_pressure(ξ_lo, T_lo, φ_lo, material)
        P_hi = ideal_gas_pressure(ξ_hi, T_hi, φ_hi, material)
        κ_lo = effective_permeability(ξ_lo, T_lo, material)
        κ_hi = effective_permeability(ξ_hi, T_hi, material)
        μ_lo = mixture_gas_viscosity(ξ_lo, T_lo, material)
        μ_hi = mixture_gas_viscosity(ξ_hi, T_hi, material)
        mob_face = _dt_mobility_face(κ_lo, κ_hi, μ_lo, μ_hi, d1, d2)
        v_g = -mob_face * (P_hi - P_lo) / d_LR
        return _dt_face_fluxes(
            T_lo,
            T_hi,
            ξ_lo,
            ξ_hi,
            φ_lo,
            φ_hi,
            λ_face,
            v_g,
            d_LR,
            material,
        )
    else
        return _gt_face_fluxes(T_lo, T_hi, ξ_lo, ξ_hi, λ_face, d_LR, material)
    end
end

# Boundary-cell S_conv/ρcp contribution. Generic over Dual scalar types; the
# whole chain — surface Newton, both face-flux kernels, the blended gradient —
# is differentiated in one piece so the `J̄ × ∂T_s` cross-products are exact.
function _bc_convection_contribution(
    T_b::Real,
    ξ_b::NTuple{NC,Tu},
    T_nb::Real,
    ξ_nb::NTuple{NC,Tu},
    thermal_bc::AbstractThermalBC,
    mass_bc,
    material::Material{NC},
    Δz_half::Real,
    d_to_face::Real,
    d_nb::Real,
    d1_int::Real,
    d2_int::Real,
    t::Real,
    T_ref::Real,
    ::Val{tag},
    darcy_active::Bool,
    radiation_model::RadiationModel,
    use_blowing::Bool = false,
) where {NC,Tu,tag}

    # Cell-local properties at the boundary cell (Dual-typed under AD).
    k_eff = effective_conductivity(ξ_b, T_b, material)
    ρcp = effective_heat_capacity(ξ_b, T_b, material)
    ε_eff = effective_emissivity(ξ_b, T_b, material)
    a_eff = effective_absorptivity(ξ_b, T_b, material)

    # Boundary-face flux in +z orientation: F is positive into the domain,
    # so J_z = σ·F with σ = −1 at the top face, +1 at the bottom face.
    # Computed BEFORE the surface Newton because the blowing correction
    # rebuilds B = ṁ″·c_p,g/h from the outward flux −F.
    F = _bc_boundary_face_fluxes(
        T_b,
        ξ_b,
        mass_bc,
        material,
        d_to_face,
        t,
        Val(tag),
        darcy_active,
    )
    σ = tag === :top ? -1.0 : 1.0
    J_bc = ntuple(j -> σ * F[j], Val(NC))

    # Surface temperature — mirrors `_bc_thermal_contribution`: the
    # transpiration mass flux is hard-zeroed (the residual's surface balance
    # never sees ṁ_g in the enthalpy term), while blowing reads the true
    # outward flux −F; Dirichlet closed-form, Dual-aware Newton otherwise.
    T_s = if thermal_bc isa DirichletThermalBC
        evaluate_param(thermal_bc.T, t)
    else
        bstate = BoundaryState(
            Float64(ForwardDiff.value(T_b)),
            0.0,
            0.0,
            0.0,
            Float64(ForwardDiff.value(T_b)),
            Float64(ForwardDiff.value(k_eff)),
            Float64(ForwardDiff.value(Δz_half)),
            false,
            0,
        )
        J_zero = ntuple(_ -> zero(Tu), Val(NC))
        J_blowing = use_blowing ? ntuple(j -> -F[j], Val(NC)) : J_zero
        result = solve_boundary_temperature!(
            bstate,
            T_b,
            k_eff,
            Δz_half,
            thermal_bc,
            material,
            J_zero,
            T_ref,
            t,
            ε_eff;
            a_eff = a_eff,
            radiation_model = radiation_model,
            use_blowing = use_blowing,
            J_blowing = J_blowing,
        )
        result.T_s
    end

    # Blended cell-center gradient in +z orientation: the surface is the
    # pseudo-neighbor at distance Δz_half. Identical at both boundaries once
    # the slopes are written in +z:
    #   top:    s_face = (T_s − T_b)/Δz_half,  s_int = (T_b − T_nb)/d_nb
    #   bottom: s_face = (T_b − T_s)/Δz_half,  s_int = (T_nb − T_b)/d_nb
    s_face = tag === :top ? (T_s - T_b) / Δz_half : (T_b - T_s) / Δz_half
    s_int = tag === :top ? (T_b - T_nb) / d_nb : (T_nb - T_b) / d_nb
    dTdx = (d_nb * s_face + Δz_half * s_int) / (d_nb + Δz_half)

    # Interior-face flux in +z orientation.
    J_int = if tag === :top
        _bcc_interior_face_fluxes(
            T_nb,
            T_b,
            ξ_nb,
            ξ_b,
            d1_int,
            d2_int,
            d_nb,
            material,
            darcy_active,
        )
    else
        _bcc_interior_face_fluxes(
            T_b,
            T_nb,
            ξ_b,
            ξ_nb,
            d1_int,
            d2_int,
            d_nb,
            material,
            darcy_active,
        )
    end

    # S_conv = −Σ_g cp_g(T_b) · (J_bc + J_int)/2 · dTdx (face order irrelevant
    # to the average; reuses gas_transport's unrolled kernel).
    S_conv = _gt_S_conv_unrolled(T_b, dTdx, J_bc, J_int, material)
    return S_conv / ρcp
end

"""
	apply!(dT, ::BCConvectionOperator, state, mesh, material, bc_set, t,
	       darcy_active, radiation_model, use_blowing = false) -> nothing

Accumulate the full boundary-cell `S_conv/ρcp_eff` contribution into
`dT[bdy_cell]` for both boundaries. No-op on a single-cell mesh.
"""
function apply!(
    dT::AbstractVector,
    ::BCConvectionOperator,
    state,
    mesh,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    t::Real,
    darcy_active::Bool,
    radiation_model::RadiationModel,
    use_blowing::Bool = false,
) where {NC}

    for boundary_id in (1, 2)
        geom = _bcc_geometry(state, mesh, boundary_id)
        geom.has_nb || continue
        thermal_bc =
            boundary_id == 1 ? get_thermal_bc(bc_set, Val(:top)) :
            get_thermal_bc(bc_set, Val(:bottom))
        mass_bc =
            boundary_id == 1 ? get_mass_bc(bc_set, Val(:top)) :
            get_mass_bc(bc_set, Val(:bottom))
        tag = boundary_id == 1 ? Val(:top) : Val(:bottom)
        T_ref = ambient_temperature(thermal_bc, t, Float64(T_REF))

        T_b = state.T[geom.cell_id]
        T_nb = state.T[geom.nb_id]
        ξ_b = ntuple(j -> state.ξ[j, geom.cell_id], Val(NC))
        ξ_nb = ntuple(j -> state.ξ[j, geom.nb_id], Val(NC))

        dT[geom.cell_id] += _bc_convection_contribution(
            T_b,
            ξ_b,
            T_nb,
            ξ_nb,
            thermal_bc,
            mass_bc,
            material,
            geom.Δz_half,
            geom.d_to_face,
            geom.d_nb,
            geom.d1_int,
            geom.d2_int,
            t,
            T_ref,
            tag,
            darcy_active,
            radiation_model,
            use_blowing,
        )
    end
    return nothing
end

# Pack/unpack: input vector is `[T_b, ξ_b…, T_nb, ξ_nb…]`.
@inline function _bcc_pack!(
    x::AbstractVector,
    u::AbstractVector,
    n_cells::Int,
    NC::Int,
    cell_id::Int,
    nb_id::Int,
)
    @inbounds begin
        x[1] = u[cell_id]
        x[2+NC] = u[nb_id]
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            x[1+j] = u[ξ_off+cell_id]
            x[2+NC+j] = u[ξ_off+nb_id]
        end
    end
    return nothing
end

# `Val(NC)` so the `ntuple`s constant-fold inside the AD closure (see
# `_bct_unpack`).
@inline function _bcc_unpack(x::AbstractVector, ::Val{NC}) where {NC}
    @inbounds T_b = x[1]
    ξ_b = ntuple(j -> @inbounds(x[1+j]), Val(NC))
    @inbounds T_nb = x[2+NC]
    ξ_nb = ntuple(j -> @inbounds(x[2+NC+j]), Val(NC))
    return T_b, ξ_b, T_nb, ξ_nb
end

"""
	contribute_to_jacobian!(jac, op::BCConvectionOperator, u, ctx, t, ws,
	                        state, mesh, material, layout, bc_set,
	                        darcy_active, radiation_model) -> nothing

Per-boundary AD path. Each boundary writes a `1 × 2(1 + NC)` energy-row
block at T-row[bdy] over the `(T, ξ)` columns of the boundary cell and
its interior neighbor. AD runs through the surface Newton, both face-flux
kernels, and the blended gradient in one piece, so every cross-product
(`J̄ × ∂T_s`, `∂J̄ × dTdx`) lands exactly.
"""
function contribute_to_jacobian!(
    jac::StructuredJacobian,
    ::BCConvectionOperator,
    u::AbstractVector,
    ctx::ResidualContext,
    t::Real,
    ws::BCConvectionADWorkspace{NC},
    state,
    mesh,
    material::Material{NC},
    layout,
    bc_set::BoundaryConditionSet,
    darcy_active::Bool,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    use_blowing::Bool = false,
) where {NC}

    n_cells = layout.n_cells
    A = jac.A
    smap = jac.smap

    # Per-branch BC lookups + barrier (see `_bct_assemble_boundary!`): keeps
    # the AD closure concretely typed per boundary.
    geom_top = _bcc_geometry(state, mesh, 1)
    geom_top.has_nb && _bcc_assemble_boundary!(
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
        darcy_active,
        radiation_model,
        use_blowing,
    )
    geom_bot = _bcc_geometry(state, mesh, 2)
    geom_bot.has_nb && _bcc_assemble_boundary!(
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
        darcy_active,
        radiation_model,
        use_blowing,
    )
    return nothing
end

function _bcc_assemble_boundary!(
    A,
    smap,
    ws::BCConvectionADWorkspace{NC},
    u::AbstractVector,
    n_cells::Int,
    geom,
    thermal_bc,
    mass_bc,
    tag::Val,
    material::Material{NC},
    t::Real,
    darcy_active::Bool,
    radiation_model::RadiationModel,
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
            val_nc = Val(NC),
            tag = tag,
            darcy_active = darcy_active,
            radiation_model = radiation_model,
            use_blowing = use_blowing

            (out, x_inner) -> begin
                T_b, ξ_b, T_nb, ξ_nb = _bcc_unpack(x_inner, val_nc)
                out[1] = _bc_convection_contribution(
                    T_b,
                    ξ_b,
                    T_nb,
                    ξ_nb,
                    thermal_bc,
                    mass_bc,
                    material,
                    geom.Δz_half,
                    geom.d_to_face,
                    geom.d_nb,
                    geom.d1_int,
                    geom.d2_int,
                    t,
                    T_ref,
                    tag,
                    darcy_active,
                    radiation_model,
                    use_blowing,
                )
                return nothing
            end
        end

    _bcc_pack!(ws.x_buf, u, n_cells, NC, geom.cell_id, geom.nb_id)
    ForwardDiff.jacobian!(ws.J_local, f!, ws.f_buf, ws.x_buf, ws.cfg)

    # Scatter the single energy row: cols 1 = T_b, 2..1+NC = ξ_b[j],
    # 2+NC = T_nb, 3+NC..2+2NC = ξ_nb[j].
    gr = geom.cell_id
    @inbounds begin
        fast_sparse_add!(A, smap, gr, geom.cell_id, ws.J_local[1, 1])
        fast_sparse_add!(A, smap, gr, geom.nb_id, ws.J_local[1, 2+NC])
        for j = 1:NC
            ξ_off = n_cells + (j - 1) * n_cells
            fast_sparse_add!(A, smap, gr, ξ_off + geom.cell_id, ws.J_local[1, 1+j])
            fast_sparse_add!(A, smap, gr, ξ_off + geom.nb_id, ws.J_local[1, 2+NC+j])
        end
    end
    return nothing
end
