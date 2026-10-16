"""
Gas-coupled energy sources for the non-conservative (Formulation A)
temperature equation:

    ρcₚ_eff ∂T/∂t = −∇·q_cond − Σ_g cp_g J_g ∂T/∂x − Σ_g Δh_g(T,T_ref) (∇·J_g) + Q_rad + Q_rxn

(`ρcₚ_eff` excludes gas storage — see `Properties.effective_heat_capacity`.)

The two gas-coupled volumetric sources implemented here:

1. **`S_conv = −Σ_g cp_g(T) · J̄_g · ∂T/∂x`** — the **advection** piece of the
   conservative-form gas-enthalpy divergence
   `∇·(Σ_g cp_g (T−T_ref) J_g) = Σ_g cp_g J_g · ∇T + (T−T_ref) Σ_g cp_g ∇·J_g`.
   Always present.

2. **`S_gen = −Σ_g Δh_g(T, T_ref) · (∇·J_g)`** — the gas-generation enthalpy
   sink (the divergence-completion piece; historically mislabeled "compression
   work"). Using `∇·J_g = ṁ'''_g` this is `−Σ_g Δh_g(T) · ṁ'''_g`: the sensible
   enthalpy carried away by gas *at the local temperature where it is generated*,
   measured from the `T_ref` datum. Dropped in `:standard` and added back as a
   cell-centered source when `energy_form = :with_generation_sink`. `T_ref` is
   the surface ambient temperature.
   **Datum requirement:** `S_gen` is consistent only with `T_ref`-datum reaction
   heats — see `compute_gas_generation_sink!`.

**Convention (S_conv only).** Gas enthalpy transport is carried **entirely**
by the volumetric advection term `S_conv`, with a plain convective+radiative
surface BC and no surface gas-enthalpy term (the ThermaKin Eq.17 / FDS / Gpyro
convention). Note the algebraic identity (constant cp):

```
∫(S_conv + S_gen) dz = −G,    G ≡ Σ_g ṁ_g · Δh_g(T_s, T_ref) > 0   (surface outflow)
```

so `∫S_conv = −G − ∫S_gen`: `S_conv` already carries all but `|∫S_gen|` of the
outflow enthalpy `G`. A surface term applying the full `G` on top of `S_conv`
would double-count the advective enthalpy by `|∫S_conv|`.

Reference: ThermaKin Eq 17 (advection term only), Bird/Stewart/Lightfoot Ch 11.
The approximation costs are visible at runtime via the diagnostics callback.
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

# ============================================================================
# Gas-Generation Enthalpy Sink (opt-in via `energy_form = :with_generation_sink`)
# (distinct from compression work)
# ============================================================================

"""
	compute_gas_generation_sink!(S_gen::Vector{Float64},
	                              div_J::Matrix{Float64},
	                              mesh::Mesh1D{NC},
	                              material::Material{NC},
	                              T_ref::Float64) where NC

Cell-centered gas-generation enthalpy sink for the optional `:with_generation_sink`
energy form (the sensible enthalpy carried away by gas at the local temperature
where it is generated; distinct from compression work):

```
S_gen,i = −Σ_{g: is_gas} Δh_g(T_i, T_ref) · (∇·J_g)_i      [W/m³]
```

where `Δh_g(T_i, T_ref) = ∫_{T_ref}^{T_i} cp_g dτ` is the true sensible-enthalpy
integral (midpoint rule, via `delta_enthalpy`). Using the integral rather than the
single-point product `(T_i − T_ref)·cp_g(T_i)` keeps the term an exact enthalpy
difference; `T_ref` is the ambient temperature from the BC.

**DATUM REQUIREMENT — read before enabling.** `S_gen` is consistent ONLY when
every reaction heat `h_r` is referenced to the `T_ref` (≈298.15 K) ambient
datum, formation-style: `h_r = Σ_j ν_j·h_j(T_ref)`. If `h_r` is instead the
heat of reaction *at the local reaction temperature* — the standard DSC /
inverse-calibration interpretation, and the usual meaning of ThermaKin- or
Gpyro-calibrated heats — then the gas products' sensible enthalpy
`Δh_g(T, T_ref)` is already inside `h_r(T)`, and the exact temperature
equation contains conduction + `S_conv` + `Q_rxn` only. Enabling this option
then subtracts `Σ_g Δh_g·ṁ'''_g` a SECOND time: 0.66–1.20 MJ per kg of
volatiles at 900 K for `cp_g` = 1.1–2.0 kJ/(kg·K) — comparable to typical
heats of pyrolysis. With DSC-style heats the default `:standard` form is
already the consistent configuration; leave this term off.
`to_problem_def` warns once when the option is selected.

**Caveat**: even with `T_ref`-datum heats, enabling this option does NOT
recover full conservative form on its own. The full conservative-equivalent
assembly also requires gases in `ρcp_eff` (currently excluded by
Formulation A; see `Properties.effective_heat_capacity`) and the
condensed-product sensible-generation term `−Σ_{S/L} Δh_j·ω̇_j` (not
implemented).
"""
function compute_gas_generation_sink!(
    S_gen::AbstractVector,
    div_J::AbstractMatrix,
    state,
    mesh,
    material::Material{NC},
    T_ref::Real,
) where {NC}

    n_cells_val = mesh.n_active_cells

    @inbounds for i = 1:n_cells_val
        T_i = state.T[i]
        S_gen[i] = -_gas_generation_sink_unrolled(T_i, div_J, i, material, T_ref)
    end

    return nothing
end

# Compile-time unrolled inner sum: Σ_{g: is_gas} Δh_g(T_i, T_ref) · (∇·J_g)_i.
# Uses the true sensible-enthalpy integral `delta_enthalpy` (midpoint rule).
@generated function _gas_generation_sink_unrolled(
    T_i::Real,
    div_J::AbstractMatrix,
    i::Int,
    material::Material{NC},
    T_ref::Real,
) where {NC}

    exprs = Expr[]
    push!(exprs, :(s = zero(promote_type(typeof(T_i), eltype(div_J), typeof(T_ref)))))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                cmp = material.components[$j]
                s += delta_enthalpy(cmp, T_i, T_ref) * div_J[$j, i]
            end
        end)
    end
    push!(exprs, :(return s))
    return Expr(:block, exprs...)
end
