"""
Gas convective energy source of the temperature equation

    ρc_p ∂T/∂t = −∇·q_cond + S_conv + Q_rad + Q_rxn,    S_conv = −Σ_g cp_g J_g ∂T/∂z

`ρc_p` sums every component, the pore gas included (see
`Properties.effective_heat_capacity`), and `J_g` is the gas mass flux
relative to the skeleton. `S_conv` is the advection piece of the gas-enthalpy
divergence: expanding `∂_z(Σ_g h_g J_g)` against the species balances leaves
`Σ_g cp_g J_g ∂_z T` once the reaction enthalpies are collected into the heats
of reaction. With those heats quoted at the local reaction temperature
(DSC-style, the convention of calibrated, ThermaKin and Gpyro heats) the
equation above is exact and has no further gas term.

The surface energy balance carries no gas-enthalpy term: the gas crosses the
surface at `T_s`, and its heating from formation to `T_s` is `S_conv`
integrated along its path. A surface term `Σ_g ṁ_g Δh_g(T_s, T_∞)` would count
that heating twice.

Reference: ThermaKin Eq. 17, Bird/Stewart/Lightfoot Ch. 11.
"""

"""
	compute_gas_convective_source!(S_conv, J_gas, q_faces, k_cell,
								   state, mesh, material)

Compute gas convective heat source at cell centers [W/m³].

S_conv,i = -Σ_g c_g(T_i) × J̄_g,i × (∂T/∂x)_i

where:
- J̄_g,i = average of gas mass fluxes at bounding faces
- (∂T/∂x)_i = opposite-distance-weighted mean of the one-sided slopes
  (`opposite_distance_gradient`) — second order on nonuniform meshes
- Boundary cells blend the interior one-sided slope with the **boundary-face
  slope at the solved surface temperature**, recovered AD-cleanly from the
  applied thermal-BC flux: `s_face = −q_faces[f]/k_cell[i]`. Since the
  thermal BC writes `q[f] = ∓k_i·(T_s − T_i)/Δz_half` (with `T_s` from the
  surface Newton), the ratio is exactly the `±(T_s − T_i)/Δz_half` face slope
  with full AD sensitivities — `boundary_states.T_surface` is Float64-stripped
  and must not be read here. The surface acts as a pseudo-neighbor at distance
  `Δz_half = V/(A_face·2)` in the same opposite-distance weighting interior
  cells use, so the scheme stays second order up to the boundary.
  (A single-cell mesh has no usable stencil; `S_conv = 0` there.)

**Must be called after the thermal BC pass** (`apply_thermal_bc_flux!`) has
filled the boundary entries of `q_faces`.

# Arguments
- `S_conv`: Output vector of volumetric sources [W/m³], one per cell
- `J_gas`: Gas mass fluxes J[j,f] [kg/(m²·s)] at faces
- `q_faces`: Total heat flux at faces [W/m²] with thermal-BC boundary entries
- `k_cell`: Cell-centered effective conductivities [W/(m·K)] — the same values
  the surface Newton used (`prop_cache.k`)
- `state`: Geometry/state cache
- `mesh`: Mesh with current state
- `material`: Material definition
"""
function compute_gas_convective_source!(
    S_conv::AbstractVector,
    J_gas::AbstractMatrix,
    q_faces::AbstractVector,
    k_cell::AbstractVector,
    state,
    mesh,
    material::Material{NC},
) where {NC}

    n_cells_val = mesh.n_active_cells

    @inbounds for i = 1:n_cells_val
        T_i = state.T[i]

        # Bounding faces via cell face_ids
        cell = mesh.cells[i]
        f_L, f_R = cell.face_ids

        # Neighbor cell IDs (0 = boundary)
        left_id = mesh.faces[f_L].cell_ids[1]
        right_id = mesh.faces[f_R].cell_ids[2]

        if left_id == 0 && right_id == 0
            S_conv[i] = zero(eltype(S_conv))
            continue
        elseif left_id == 0
            # Bottom boundary cell: surface is the pseudo-neighbor on the left
            # at distance Δz_half, with slope −q/k (= (T_i − T_s)/Δz_half at
            # the solved T_s).
            s_face = -q_faces[f_L] / k_cell[i]
            dz_half = state.volume[i] / state.A_face[f_L] / 2
            T_R = state.T[right_id]
            d_R = abs(state.cell_centers[right_id] - state.cell_centers[i])
            s_int = (T_R - T_i) / d_R
            dTdx = (d_R * s_face + dz_half * s_int) / (dz_half + d_R)
        elseif right_id == 0
            # Top boundary cell: surface is the pseudo-neighbor on the right.
            s_face = -q_faces[f_R] / k_cell[i]
            dz_half = state.volume[i] / state.A_face[f_R] / 2
            T_L = state.T[left_id]
            d_L = abs(state.cell_centers[i] - state.cell_centers[left_id])
            s_int = (T_i - T_L) / d_L
            dTdx = (dz_half * s_int + d_L * s_face) / (d_L + dz_half)
        else
            T_L = state.T[left_id]
            T_R = state.T[right_id]
            d_L = abs(state.cell_centers[i] - state.cell_centers[left_id])
            d_R = abs(state.cell_centers[right_id] - state.cell_centers[i])
            dTdx = opposite_distance_gradient(T_L, T_i, T_R, d_L, d_R)
        end

        S_conv[i] = _gas_convective_source_unrolled(T_i, dTdx, J_gas, f_L, f_R, material)
    end

    return nothing
end

"""
Compile-time unrolled gas convective source computation.
Avoids runtime indexing into heterogeneous component tuple.
"""
@generated function _gas_convective_source_unrolled(
    T_i::Real,
    dTdx::Real,
    J_gas::AbstractMatrix,
    f_L::Int,
    f_R::Int,
    material::Material{NC},
) where {NC}

    exprs = Expr[]
    push!(exprs, :(S = zero(promote_type(typeof(T_i), typeof(dTdx), eltype(J_gas)))))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                J_avg = (J_gas[$j, f_L] + J_gas[$j, f_R]) / 2
                cp_j = material.components[$j].heat_capacity(T_i)
                S -= cp_j * J_avg * dTdx
            end
        end)
    end
    push!(exprs, :(return S))
    return Expr(:block, exprs...)
end
