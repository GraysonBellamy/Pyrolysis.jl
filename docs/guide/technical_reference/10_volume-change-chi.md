# 10. Volume Change, Swelling, Intumescence, and Lateral Shrinkage

This chapter derives the complete volume-change model of Pyrolysis.jl: how reaction-driven mass exchange between condensed and gaseous phases produces a per-cell volumetric strain rate $\theta$, how that strain is converted into axial mesh motion (the ALE mesh velocity), and how a column-global *lateral* shrinkage / swelling / intumescence law deforms the cross-sectional area $A(t)$ through the pyrolysis-progress state variable $\chi$. We give the continuous balance, the exact discrete update rules as implemented, the dilation correction that keeps mesh-bound solids Lagrangian, and a careful comparison to the deformation treatments in ThermaKin2Ds, Gpyro, and FDS. Every equation here is grounded in the source: `src/Physics/volume_change.jl`, `src/Residual/residual.jl`, `src/Problem/state_cache.jl`, `src/Problem/initial_conditions.jl`, and `src/Materials/materials.jl`.

The axial mesh-motion machinery (cumulative integration of $\theta_L$, the Geometric Conservation Law, the ALE advection operators, and conservative remapping) is developed in §11, *The ALE Framework*; here we develop only the *constitutive* content — the strain rates and the cross-section law that feed §11 — and the species dilation correction that closes the species balance under mesh motion. The state-vector layout that carries the $\chi$ block is detailed in §13, *Finite-Volume Discretization*, and §15, *Temporal Integration & Structured Jacobian*.

**Index convention for this chapter.** Following the binding resolution in §2 (overload I1), $i$ is the **cell index** ($1\le i\le n$), $j$ is the **component/species index** ($1\le j\le N_C$), and $r$ is the **reaction index**. The symbol $\chi$ is *always* pyrolysis progress (overload A1); the cross-section ratio is $A/A_0=\mathrm{law}(\bar\chi)$ and is never written as $\chi$. The symbol $\theta$ is *always* a volumetric strain rate (overload A2); stoichiometry is $\nu_{r,j}$.

---

## 10.1 The volumetric strain rate $\theta$

### 10.1.1 Continuous derivation from the mixture continuity equation

Consider a material control volume $V$ whose condensed-phase matrix (solids and liquids) forms a skeleton through which gas can percolate. The matrix occupies a volume that is the sum of the partial volumes of each matrix-forming component. If component $j$ has bulk (pure-phase) density $\rho_j$ and mass concentration $\xi_j$ (kg of $j$ per m³ of bulk material; the **primary species state**, see §2), then the mass of $j$ in the cell is $m_j=\xi_j V$, and the partial (pure-phase) volume it occupies is $m_j/\rho_j=\xi_j V/\rho_j$.

The total matrix volume is the sum of these partial volumes weighted by how much each phase contributes to the *skeleton* that the mesh must track. Introduce the dimensionless **swelling factor** $\gamma_j\in[0,1]$ to encode that weighting:

$$
V \;=\; \sum_{j=1}^{N_C} \gamma_j \,\frac{m_j}{\rho_j}
\;=\; V\sum_{j=1}^{N_C} \gamma_j \,\frac{\xi_j}{\rho_j}.
\tag{10.1}
$$

Here $\gamma_j=1$ means component $j$ fully claims pure-phase volume in the skeleton (a normal solid: when it depletes, the matrix shrinks), $\gamma_j=0$ means it claims none (a normal gas: it escapes through the pores and does not set the matrix dimensions), and $\gamma_j\in(0,1)$ is a partial claim (an intumescent gas that is partly trapped and inflates the char). Equation (10.1) is the constitutive *definition* of the matrix volume in terms of composition; it is exact for an additive-volume mixture.

Differentiate (10.1) with respect to time at a fixed material point, holding $\rho_j$ as a (slowly varying) property and the $\gamma_j$ as constants:

$$
\frac{dV}{dt} \;=\; V\sum_{j=1}^{N_C}\gamma_j\,\frac{1}{\rho_j}\frac{d\xi_j}{dt}
\;+\; V\sum_{j=1}^{N_C}\gamma_j\,\xi_j\,\frac{d}{dt}\!\left(\frac{1}{\rho_j}\right)
\;+\; \frac{dV}{dt}\sum_{j=1}^{N_C}\gamma_j\frac{\xi_j}{\rho_j}.
$$

Using (10.1), $\sum_j\gamma_j\xi_j/\rho_j=1$, so the last term reproduces $dV/dt$. The standard reduction — the one the code adopts — is to *neglect the density-rate term* $\xi_j\,d(1/\rho_j)/dt$ (thermal expansion of the pure phases is small compared with reaction-driven mass exchange) and to read off the **relative volume change rate** — the volumetric strain rate — directly from the mass-source term:

$$
\boxed{\;\theta \;\equiv\; \frac{1}{V}\frac{dV}{dt}
\;=\; \sum_{j=1}^{N_C}\gamma_j\,\frac{1}{\rho_j(T)}\,\frac{d\xi_j}{dt}\;}
\tag{10.2}
$$

with SI units of $\mathrm{s^{-1}}$. The interpretation is transparent: a unit decrease in the mass concentration of solid $j$ ($d\xi_j/dt<0$) frees a pure-phase volume of $1/\rho_j$ per unit mass, contracting the matrix ($\theta<0$); a gas that does not contribute ($\gamma_g=0$) is invisible to $\theta$ regardless of how fast it is produced.

> **Implementation note.** $\theta$ is computed by `volume_change_rate(dξ, T, material)` / the dispatch helpers `_volume_change_rate`, in `src/Physics/volume_change.jl` (the inline $N_C\in\{1,2,3\}$ overloads at lines 15–63 and the `@generated` form for $N_C>3$ at lines 66–90). The exact loop is
> ```julia
> θ += γⱼ * dξ[j] / ρⱼ        # ρⱼ = comps[j].density(T)
> ```
> evaluated **only when `γⱼ > 0`** (an early-exit guard: components with $\gamma_j=0$ are skipped, so no density evaluation is wasted and the term is exactly zero). There is **no minus sign** in the code: the sign of $\theta$ is carried entirely by the sign of $d\xi_j/dt$, which is negative for a depleting reactant (the species source term packs reactants with coefficient $-1$; see §6.5). The density $\rho_j(T)$ is the *bulk* density property evaluated at the cell temperature — not the intrinsic/skeletal density $\rho_{i,j}$ (which is reserved for the porosity formula, see §5.3 and §2 overload A7).

### 10.1.2 What $d\xi_j/dt$ is in (10.2)

A subtle and load-bearing point: the $d\xi_j/dt$ fed into (10.2) is the **reaction source term only**, not the full species RHS (which also includes gas transport and ALE advection). This is correct because $\theta$ measures the rate at which mass is *converted between phases by chemistry* — the only process that changes the additive-volume balance (10.1). Transport of an already-gaseous species through the pores does not change the matrix volume, and the ALE advection / dilation terms are bookkeeping for the moving frame, not new mass exchange. Including them would double-count the volume change that ALE motion is itself representing.

> **Implementation note.** In the ALE residual (`_residual!(..., ::ALE)`, `src/Residual/residual.jl:91–222`) the reaction source matrix `ws.rhs_cache.dξ_rxn` is populated by the physics RHS at step 4 (line 123) and is the array passed as `precomputed_dξ_rxn` to `_compute_ale_mesh_velocity!` (line 156). Inside that function (lines 389–401), `dξ_dt[j] = precomputed_dξ_rxn[j, cell_idx]`, and then `θ_cell = Physics.volume_change_rate(dξ_dt, T_cell, material)`. So $\theta$ is built from `dξ_rxn`, the pure reaction source — exactly the per-cell strain rate of (10.2).

### 10.1.3 Absolute cell volume rate

The absolute rate of change of a cell's volume follows immediately:

$$
\frac{dV_i}{dt} \;=\; V_i\,\theta_i.
\tag{10.3}
$$

> **Implementation note.** `compute_cell_volume_rates!(dV, state, dξ_matrix, material)` (`src/Physics/volume_change.jl:126–149`) reads $V_i$ and $T_i$ from the typed `StateCache`, takes $d\xi/dt$ from a precomputed matrix, evaluates $\theta_i$ via `volume_change_rate`, and writes `dV[i] = V * θ`. This batch helper is provided for diagnostics; the production ALE path uses $\theta_i$ directly in the mesh-velocity integration (§10.6) rather than materializing $dV_i/dt$.

---

## 10.2 Swelling-factor semantics

The swelling factor $\gamma_j$ is a per-component scalar field on every `AbstractComponent`. Its defaults encode the standard physics:

| Phase | Default $\gamma$ | Meaning | Effect on $\theta$ |
|---|---|---|---|
| Solid | $1.0$ | Skeleton-forming; depletion frees pure-phase volume | $\theta<0$ when solid is consumed |
| Liquid | $0.0$ | Fills pore space; default treats it as volume-neutral | none unless overridden |
| Gas | $0.0$ | Escapes through pores; does not set matrix dimensions | none |

> **Implementation note.** Defaults are set in the convenience constructors of `src/Materials/components.jl`: `SolidComponent(...; γ::Real = 1.0, ...)`, `GaseousComponent(...; γ::Real = 0.0, ...)`, and `LiquidComponent(...; γ::Real = 0.0, ...)`. The field is `swelling_factor::Float64` on each struct, validated to $\gamma\in[0,1]$ at construction. All three γ consumers — the θ kernels, the mixture density (5.6), and the mixing-rule volume fractions (5.4) — use the identical `γ > 0` participation test and the identical $\gamma_j\,\xi_j/\rho_j$ weighting.

**Modeling guidance.**

- **Non-charring shrinkage (e.g. melting/charring polymer to gas):** keep the virgin solid at $\gamma=1$ and the gas product at $\gamma=0$. As the solid converts to gas, $\theta<0$ and the column recedes — the classic shrinking-slab picture.
- **Charring with a rigid char:** give the char $\gamma=0$ to "freeze" its volume contribution (`rigid_char = SolidComponent(:char, ...; γ=0.0)`, the documented example at `components.jl:182`); then only the mass converted out of the cell drives axial motion, and the char does not re-inflate or collapse with the mesh. (This is a modeling choice, not a default; the default char is $\gamma=1$.)
- **Intumescence (gas-producing expansion):** assign the trapped/blowing gas a swelling factor $\gamma_g\in(0,1]$. A positive $d\xi_g/dt$ (gas being produced) with $\gamma_g>0$ contributes $\theta>0$, i.e. the matrix *inflates*. Combined with a swelling lateral law (§10.5) this reproduces the order-of-magnitude expansion of intumescent coatings. Axial intumescence (through-thickness swell) is captured by $\theta>0$ in (10.2); radial intumescence (lateral swell of the column) is captured by $\mathrm{law}(\bar\chi)>1$ in (10.8).

Because $\gamma_j$ is a free per-component parameter, the model spans the full range from pure shrinkage ($\gamma=1$ solid, $\gamma=0$ gas) through volume-neutral ($\gamma=0$ everywhere) to intumescent ($\gamma_g>0$).

---

## 10.3 The pyrolysis-progress state $\chi$

### 10.3.1 Definition and governing ODE

To make the *lateral* (cross-section) deformation depend on how far decomposition has advanced, Pyrolysis.jl introduces a per-cell **pyrolysis progress** $\chi_i$ (dimensionless): the cumulative mass of *dry* pyrolysis gas released in cell $i$, normalized by the cell's initial dry-solid mass. Its governing ODE is

$$
\boxed{\;\frac{d\chi_i}{dt}
\;=\; \frac{S_{\mathrm{dry},i}\,V_i}{m_{\mathrm{dry},i}}\;},
\qquad
S_{\mathrm{dry},i} \;=\; \sum_{j\in\mathrm{dry\ gas}}\left.\frac{d\xi_j}{dt}\right|_{\mathrm{rxn}},
\tag{10.4}
$$

where $S_{\mathrm{dry},i}$ (units $\mathrm{kg\,m^{-3}\,s^{-1}}$) is the **dry-gas production rate** — the sum of reaction source terms over the components flagged as dry pyrolysis gas — $V_i$ is the current cell volume, and $m_{\mathrm{dry},i}$ is the constant initial dry-solid mass in the cell. The product $S_{\mathrm{dry},i}V_i$ is the dry-gas mass-release *rate* (kg/s); integrating (10.4) gives

$$
\chi_i(t) \;=\; \frac{1}{m_{\mathrm{dry},i}}\int_0^t S_{\mathrm{dry},i}(t')\,V_i(t')\,dt'
\;=\; \frac{m_{\mathrm{released,dry},i}(t)}{m_{\mathrm{dry},i}},
\tag{10.5}
$$

i.e. exactly the fractional dry-gas mass release. By construction $\chi_i$ starts at $0$ and, for a fully decomposing solid that converts entirely to dry gas, tends to $1$; with char yield it saturates below $1$.

> **Implementation note.** Computed in `_pack_chi_derivative!` (`src/Residual/residual.jl:484–517`). The exact write is
> ```julia
> S_dry = Σ_{j ∈ idx_dry} dξ_rxn[j, i]
> denom = max(m_init[i], eps())
> du[χ_offset + i] = S_dry * V_i / denom
> ```
> with `m_init = mesh.m_dry_solid_initial`. The `max(..., eps())` guard prevents division by zero in cells that began with no solid. Note that $V_i$ here is the *current* `state.volume[i]`, so $\chi$ accumulates the true mass release even as the cell volume changes with swelling/shrinkage — consistent with the integral form (10.5).

### 10.3.2 The initial dry-solid mass $m_{\mathrm{dry},i}$

The denominator $m_{\mathrm{dry},i}$ is the **solid** mass present in cell $i$ at $t=0$:

$$
m_{\mathrm{dry},i} \;=\; V_i(0)\sum_{j\in\mathrm{solid}}\xi_{j,i}(0).
\tag{10.6}
$$

> **Implementation note.** `populate_m_dry_solid_initial!(mesh, material)` (`src/Problem/initial_conditions.jl:256–272`) sums `mesh.cell_states[i].ξ[j]` over components for which `is_solid(material.components[j])` is true, multiplied by `mesh.cell_volumes[i]`. Crucially the sum is over **solids only** — liquids (e.g. bound moisture) are *not* counted in the denominator, even though they are matrix-forming. This is deliberate: $\chi$ is a *dry-decomposition* progress variable, and water mass is excluded from both numerator (dry-gas only) and denominator (solids only). The array is stored as `Float64` on the immutable mesh and is constant for the whole solve.

### 10.3.3 Dry-gas identification and the water-vapor caveat

The index set $\{j\in\mathrm{dry\ gas}\}$ collects all GAS-phase components *except* water vapor / moisture-derived gases. Phase is a compile-time (type-level) property, but the water-vapor exclusion is done by **name matching**:

> **Implementation note.** `dry_pyrolysis_gas_indices(material)` (`src/Materials/materials.jl:625–652`) recurses over the component tuple at compile time (`_collect_dry_pyrolysis_gas_indices`) and includes index $j$ iff `is_gas(comp) && !_is_water_vapor_name(comp.name)`. The name filter `_is_water_vapor_name` (line 645) is a **case-insensitive substring match**: a name is treated as water vapor if its string contains `"H2O"`, or (lowercased) contains `"water"`, or equals `"h2o"`. The resulting tuple is cached once at workspace construction (`ws.dry_gas_indices`) to avoid a per-residual `Symbol→String` allocation. **Caveat:** because identification is by substring, a dry-gas species accidentally named (say) `:water_gas` would be misclassified as moisture and excluded from $\chi$. Users modeling moisture must name the water-vapor component so it matches (e.g. `:water_vapor` or `:H2O`), and must avoid those substrings in genuine dry-gas names.

### 10.3.4 Why $\chi$ is reconstructed, not freely integrated

Although $\chi$ occupies a trailing block of the state vector $\mathbf{u}$ (the "$\chi$ block"; see §10.9 and §15.1) and is advanced by the ODE integrator, its derivative $d\chi_i/dt$ is *recomputed from the reaction sources on every residual evaluation* by `_pack_chi_derivative!` rather than being a free dynamical degree of freedom with independent physics. Because the RHS for $\chi$ is a pure function of the (reaction part of the) other states and geometry, $\chi$ stays exactly synchronized with cumulative dry-gas production regardless of the integrator's internal representation. For identity-law materials (no lateral shrinkage), the $\chi$ block is absent from the state vector entirely (the `chi_range` is empty), and `_pack_chi_derivative!` returns immediately (line 494) — so $\chi$ carries zero memory/compute cost when it is not needed. When a law is present but identity-valued, the block is zeroed (lines 498–502) because $\chi$ is undefined when the cross-section is fixed.

---

## 10.4 Mass-weighted column average $\bar\chi$

The lateral-shrinkage law acts on a single, *column-global* progress scalar, the **mass-weighted column average**:

$$
\boxed{\;\bar\chi \;=\; \frac{\displaystyle\sum_{i=1}^{n} m_{\mathrm{dry},i}\,\chi_i}
{\displaystyle\sum_{i=1}^{n} m_{\mathrm{dry},i}}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i=1}^{n} m_{\mathrm{dry},i}\,\chi_i\;},
\qquad
M_{\mathrm{dry}}=\sum_{i=1}^{n} m_{\mathrm{dry},i}.
\tag{10.7}
$$

Weighting by the *initial* dry-solid mass (rather than by volume, or uniformly) means $\bar\chi$ is the fraction of the column's *total decomposable solid mass* that has been released as dry gas — a physically meaningful, conserved-denominator average. Because $M_{\mathrm{dry}}$ is fixed at setup, $\bar\chi\in[0,1]$ moves monotonically as the column decomposes.

> **Implementation note.** `_chi_bar_from_state(u, layout, mesh)` (`src/Problem/state_cache.jl:195–217`) reads the $\chi$ block from `u` and the weights from `mesh.m_dry_solid_initial`, returning the tuple $(\bar\chi, M_{\mathrm{dry}})$. The numerator `num += u[χ_offset+i] * m_i` accumulates in the scalar type of `u` (which is `ForwardDiff.Dual` in a Jacobian/sensitivity pass), while the denominator `M` accumulates in `Float64`. The division `num / M` therefore returns a Dual when `u` is Dual — the average is **AD-safe**, carrying sensitivities of $\bar\chi$ with respect to every state through which $\chi_i$ depends. Guards: an empty `m_dry_solid_initial` or empty `chi_range` returns $(0,0)$, and $M_{\mathrm{dry}}\le\varepsilon$ returns $(0,0)$ (line 215). The public wrapper is `column_chi_bar_from_state(u, layout, mesh)` (`src/Residual/residual.jl:574–575`).

---

## 10.5 The lateral-shrinkage law $A(t)=A_0\,\mathrm{law}(\bar\chi)$

### 10.5.1 Definition

The column cross-sectional area evolves as a prescribed function of the global progress:

$$
\boxed{\;A(t) \;=\; A_0\,\mathrm{law}(\bar\chi)\;},
\tag{10.8}
$$

where $A_0$ is the calibration-reference area fixed at setup and $\mathrm{law}:\,[0,1]\to(0,\infty)$ is a user-supplied smooth function (or `nothing`, meaning identity). The area $A(t)$ is uniform along $z$ (the column has one cross-section at any instant) but time-varying.

> **Implementation note.** Stored as `material.lateral_shrinkage_law::LSL` (`src/Materials/materials.jl:171`), a concrete `Function` subtype or `Nothing`. The constructor validates that it is a `Function` or `nothing` (line 186) and the public `Material(...; lateral_shrinkage_law = nothing)` keyword defaults to identity. $A_0$ is `mesh.cross_section_area_initial` (`src/Geometry/mesh.jl:109`), set once at mesh construction. The update is performed by `update_cross_section_area!(cache, material, mesh, u, layout, ::WithChi)` (`src/Problem/state_cache.jl:175–192`):
> ```julia
> χ̄, M = _chi_bar_from_state(u, layout, mesh)
> A = mesh.cross_section_area_initial * material.lateral_shrinkage_law(χ̄)
> cache.A_section[] = A
> ```
> The result is stored in `cache.A_section[]` carrying `eltype(u)` (Dual-safe), and the function returns $(A,\bar\chi,M_{\mathrm{dry}})$ for downstream use. For identity-law materials it short-circuits and returns the unchanged $A_0$ (lines 184–186); for `NoChi` mode it is a no-op returning the current area (lines 166–173). Note the code uses the field name `cross_section_area_initial` for $A_0$; the mutable `cross_section_area` holds the current $A(t)$.

### 10.5.2 Canonical laws

| Law | $\mathrm{law}(\bar\chi)$ | Behavior |
|---|---|---|
| Identity (default) | $1$ | Constant cross-section $A=A_0$ (standard 1D pyrolysis) |
| Linear shrink | $1-k\bar\chi$, $0<k\le 1$ | Area contracts linearly with extent of decomposition |
| Linear swell | $1+k\bar\chi$, $k>0$ | Area expands (intumescence) |
| Smooth (sigmoid) | $1-\tfrac{k}{1+e^{-s(\bar\chi-\bar\chi_0)}}$ | Onset-localized transition |

For example, the verification fixture uses $\mathrm{law}(\bar\chi)=1-0.2\,\bar\chi$, so at $\bar\chi=0.5$ the area is $A=0.9\,A_0$ (`test/geometry/test_variable_cross_section.jl:9–19`, which asserts `A ≈ A0 * 0.9`).

### 10.5.3 First-principles status

The lateral law is **empirical / phenomenological** — a calibration handle, not a derived constitutive relation. There is at present no first-principles model (e.g. a percolation- or homogenization-based predictor) in the code that selects the functional form of $\mathrm{law}(\bar\chi)$ from microstructure; the user supplies it. This is by design: a 1D solver cannot resolve the lateral (transverse) dynamics that actually set $A(t)$, so the law condenses that unresolved physics into a single $\bar\chi$-dependence. Its primary uses are (i) matching cone-calorimeter samples that visibly shrink or intumesce so the per-unit-*original-area* mass-loss-rate gauge is correct (§10.8.2), and (ii) order-of-magnitude intumescent-thickness studies. Section 10.12 lists the open questions this raises.

---

## 10.6 The radial strain rate $\theta_A$ and axial/radial decomposition

### 10.6.1 Definition and chain-rule derivation

A time-varying cross-section contributes its own volumetric strain. Define the **column-global radial volumetric strain rate**

$$
\theta_A \;\equiv\; \frac{1}{A}\frac{dA}{dt}.
\tag{10.9}
$$

Differentiating (10.8) and using the chain rule (treating $A_0$ as constant):

$$
\frac{dA}{dt} \;=\; A_0\,\mathrm{law}'(\bar\chi)\,\frac{d\bar\chi}{dt},
\qquad\Longrightarrow\qquad
\boxed{\;\theta_A \;=\; \frac{A_0\,\mathrm{law}'(\bar\chi)\,\dot{\bar\chi}}{A_0\,\mathrm{law}(\bar\chi)}
\;=\; \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\,\frac{d\bar\chi}{dt}\;}.
\tag{10.10}
$$

The total column-average progress rate $\dot{\bar\chi}$ is obtained by mass-weighting the per-cell rates. Using (10.4) and (10.7),

$$
\frac{d\bar\chi}{dt}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i=1}^{n} m_{\mathrm{dry},i}\,\frac{d\chi_i}{dt}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i=1}^{n} m_{\mathrm{dry},i}\,\frac{S_{\mathrm{dry},i}V_i}{m_{\mathrm{dry},i}}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i=1}^{n} S_{\mathrm{dry},i}\,V_i.
\tag{10.11}
$$

The $m_{\mathrm{dry},i}$ cancels per cell, leaving a clean volume-weighted sum of the dry-gas source over the column, divided by the total initial dry-solid mass. This is exactly the discrete form used.

> **Implementation note.** `_compute_theta_A(material, state, χ̄, M, dξ_rxn, layout, idx_dry)` (`src/Residual/residual.jl:529–564`) computes
> ```julia
> dχ̄_dt = Σ_i ( Σ_{j ∈ idx_dry} dξ_rxn[j,i] ) * state.volume[i]   # = Σ S_dry,i V_i
> dχ̄_dt /= M                                                       # M = M_dry
> A_ratio = law(χ̄)
> dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)                        # law'(χ̄) by AD
> return dlaw_dχ̄ * dχ̄_dt / A_ratio
> ```
> matching (10.10)–(10.11) term for term. The derivative $\mathrm{law}'(\bar\chi)$ is obtained at runtime by `ForwardDiff.derivative` (line 562), so the user supplies only $\mathrm{law}$, never its derivative, and any smooth law works. **Guards** return $0$ when: the law is identity (`law === nothing`, line 539); the reaction source is unavailable (`dξ_rxn === nothing`, line 540); $M_{\mathrm{dry}}\le\varepsilon$ (line 541); the cell count is zero (line 547); or the area ratio $\mathrm{law}(\bar\chi)\le\varepsilon$ (line 560) — the last prevents the $1/\mathrm{law}(\bar\chi)$ in (10.10) from blowing up as the cross-section vanishes. The mode dispatch `_compute_theta_A_dispatched(...)` returns `0.0` for `NoChi` (line 302) and calls `_compute_theta_A` for `WithChi` (line 313).

### 10.6.2 Axial residual strain and double-counting avoidance

The per-cell strain rate $\theta_i$ from (10.2) is the *total* volumetric strain of the cell — it includes both the through-thickness (axial) deformation that the mesh tracks by moving nodes along $z$ **and** the lateral (radial) deformation that the cross-section law tracks by changing $A$. If we fed the full $\theta_i$ into the axial mesh-velocity integration *and* also changed $A(t)$, the radial volume change would be counted twice. The resolution is to subtract the column-global radial channel and feed only the **axial residual strain**

$$
\boxed{\;\theta_{L,i} \;=\; \theta_i - \theta_A\;}
\tag{10.12}
$$

into the mesh-velocity integral. Here $\theta_A$ is the same for every cell (it is a column-global quantity), and $\theta_i$ is cell-local. The cumulative axial mesh velocity is then the bottom-up integral of $\theta_{L,i}$ over cell thicknesses:

$$
w_1 = 0,
\qquad
w_{i+1} \;=\; w_i + \theta_{L,i}\,\Delta z_i
\;=\; w_i + (\theta_i-\theta_A)\,\Delta z_i,
\qquad
\Delta z_i = \frac{V_i}{A}.
\tag{10.13}
$$

Equivalently $w_{i+1}=\sum_{k=1}^{i}\theta_{L,k}\Delta z_k$, with the bottom node pinned at $w_1=0$ (the substrate is immobile; material recedes downward and the surface moves). The full development of (10.13) — its discrete Geometric Conservation Law, the node-position update $z_i^{\text{new}}=z_i^{\text{old}}+\Delta t\,w_i$, and the generalized GCL with the area-change term — is the subject of §11, *The ALE Framework*. Here we only establish the *constitutive split* (10.12) and why $\theta_A$ must be removed.

> **Implementation note.** `_compute_ale_mesh_velocity!` (`src/Residual/residual.jl:355–415`) initializes `w[1] = 0` and `cumulative_velocity = 0` (lines 382–383), then for each cell computes `Δz_cell = state.volume[cell_idx] / state.A_section[]` (line 387, i.e. $\Delta z_i=V_i/A$), evaluates `θ_cell` via `volume_change_rate` (line 401), forms `θ_L = θ_cell - θ_A` (line 407), accumulates `cumulative_velocity += θ_L * Δz_cell` (line 408), and assigns the result to the cell's right node `w[right_node] = cumulative_velocity` (lines 410–411). The radial channel `θ_A` is passed in as the keyword `θ_A` (default `0.0`), so for `NoChi` it is identically zero and (10.13) collapses to the constant-cross-section integral $w_{i+1}=w_i+\theta_i\Delta z_i$. The per-cell $\theta_i$ is also exported through the optional `θ_out` array (line 403), which is later consumed by the dilation correction (§10.7).

### 10.6.3 Sequencing of the area update and the strain integration

The two channels are computed in a specific order each residual evaluation, because $\theta_A$ depends on $\bar\chi$ (through $\mathrm{law}$ and $\mathrm{law}'$) and on the reaction sources, while the axial integration depends on $\theta_A$ and on the *current* geometry (volumes $V_i$ and area $A$). The ALE residual therefore does, in sequence (`src/Residual/residual.jl:104–161`):

1. `populate!` the state cache from $\mathbf{u}$ (T, ξ, z, χ).
2. `update_cross_section_area!` → returns $(A,\bar\chi,M_{\mathrm{dry}})$ and stores $A$ in `cache.A_section[]`.
3. `update_state_geometry!(..., ALE())` → recompute $V_i$, face areas, cell centers from $z$ and the **new** $A$ (so the $\Delta z_i$ in step 6 use the updated cross-section).
4. Properties update; then physics RHS (populates `dξ_rxn`).
5. `_compute_theta_A_dispatched(...)` → $\theta_A$ from (10.10)–(10.11).
6. `_compute_ale_mesh_velocity!(...; θ_A)` → integrate (10.13).

This ordering guarantees the radial and axial channels are mutually consistent at the single instant of the residual evaluation: the same $\bar\chi$ sets both $A$ (hence $\Delta z_i$) and $\theta_A$.

---

## 10.7 The ALE dilation correction for mesh-bound species

### 10.7.1 The problem

In the moving (ALE) frame, the mesh velocity $w$ already represents the volume change of each cell: when a cell's volume shrinks because its solid is being consumed, the nodes move so that $\Delta z_i$ decreases. But the species concentration $\xi_j$ (mass per unit *bulk* volume) is being advanced by an RHS that contains a *transport/advection* contribution written in flux form for the moving frame. For components that are **bound to the mesh** (solids and liquids, whose material velocity equals the mesh velocity), the advection flux must produce *no net change* in their cell-integrated mass when there is no reaction — they simply ride along with the shrinking cell. The discrete operator, however, leaves a residual term proportional to the cell's dilation, which would spuriously over- or under-deplete the solid. The **dilation correction** removes it.

### 10.7.2 Derivation

For a mesh-bound (Lagrangian) component, the total mass in a cell, $m_j=\xi_j V$, changes only by reaction: $dm_j/dt = S_j^{\mathrm{rxn}}V$. Expanding,

$$
\frac{d m_j}{dt} = \frac{d(\xi_j V)}{dt} = V\frac{d\xi_j}{dt} + \xi_j\frac{dV}{dt}
= V\frac{d\xi_j}{dt} + \xi_j V\,\theta,
$$

so the *concentration* rate that keeps cell mass on its reaction budget is

$$
\frac{d\xi_j}{dt} = \frac{1}{V}\frac{dm_j}{dt} - \xi_j\,\theta
= S_j^{\mathrm{rxn}} - \xi_j\,\theta.
\tag{10.14}
$$

The first term ($S_j^{\mathrm{rxn}}$) is the chemical source already in the species RHS; the second term, $-\xi_j\theta$, is the **dilation correction**. It says: when a cell shrinks ($\theta<0$), the concentration of a mesh-bound species must *rise* (mass conserved in less volume), so we *subtract* a negative number, increasing $d\xi_j/dt$; when a cell swells ($\theta>0$), the concentration falls. This is precisely the geometric concentration change that the moving frame must account for in addition to the advective flux form already applied (see §11.2 for how the advective flux and this term combine):

$$
\boxed{\;\left.\frac{d\xi_j}{dt}\right|_{\mathrm{dilation}} = -\,\xi_j\,\theta_i,
\qquad j\in\{\text{solids, liquids}\}\;}
\tag{10.15}
$$

> **Implementation note.** `_apply_ale_dilation_correction!(dξ, θ, state, material)` (`src/Residual/residual.jl:442–455`) loops over cells, reads the per-cell strain `θ_i = θ[i]` (the `θ_out` array filled during the mesh-velocity computation, §10.6.2), and calls the compile-time-unrolled `_apply_dilation_correction_unrolled!` (lines 457–473), which executes
> ```julia
> if !is_gas(material.components[j])
>     dξ[j, i] -= ξ_i[j] * θ_i
> end
> ```
> for each component. The correction is **applied to every non-gas component** (solids *and* liquids — the guard is `!is_gas`), and **skipped for gases** — the gas volume-rate term lives inside the conservative flux-form transport instead (§11.2.3, §13.8.3), so the pairing telescopes. The sign is exactly the $-\xi_j\theta_i$ of (10.15). It is applied *after* the ALE mesh-motion terms (step 7 of the residual) and *after* the reaction sources, so it sits on top of the already-assembled species RHS.

### 10.7.3 Why gases are excluded

Gases are **not** bound to the mesh: they percolate through the pores with their own velocity, and in the simplified upwind ALE formulation the gas frame is taken stationary relative to the mesh (the gas advection is handled in the lab frame by the Darcy–Fick flux plus the $+w\,\partial_z\xi$ ALE correction, §11.2/§13.8. A gas does not "carry mass with the shrinking skeleton"; its concentration responds to transport and reaction, not to the matrix dilation. Applying (10.15) to a gas would incorrectly couple the escaping gas to the solid's volume change and over-deplete it. Hence the `!is_gas` guard. This is consistent with the swelling-factor semantics of §10.2: a $\gamma=0$ gas contributes nothing to $\theta$ in the first place, so its concentration should not be driven by $\theta$ either.

---

## 10.8 Geometry coupling and observable consequences

### 10.8.1 Cell volume from state

In ALE mode every cell volume is recomputed from the nodal positions and the current cross-section on each residual evaluation:

$$
V_i \;=\; \lvert z_{i+1}-z_i\rvert\,A,
\qquad A \equiv A(t)=A_0\,\mathrm{law}(\bar\chi).
\tag{10.16}
$$

> **Implementation note.** `update_state_geometry!(cache, mesh, ::ALE)` (`src/Problem/state_cache.jl:336–354`) sets `cache.volume[i] = abs(z_right - z_left) * A` with `A = cache.A_section[]`, and face areas `cache.A_face[f] = A`, cell centers `(z_left + z_right)/2`. All quantities preserve `eltype(cache.z)` so AD passes carry geometry sensitivities (e.g. $\partial V_i/\partial\chi$ via $A$) through into every kernel that reads volumes or face areas. This is the link by which a perturbation in $\chi$ propagates to the energy and species balances: a larger $A$ enlarges every $V_i$ and every face area, rescaling divergences and the cell-integrated masses. The absolute value $\lvert\cdot\rvert$ guards against transient node inversion during a Newton/AD probe.

### 10.8.2 Mass-loss-rate gauge for variable cross-section

A measurable consequence of $A(t)$ is the mass-loss-rate (MLR) reported to the user. The surface gas flux $J_{\mathrm{top}}$ computed by the boundary solve is a *per-current-area* flux $\dot m''_{\mathrm{local}}$ (kg·m⁻²·s⁻¹ through the present surface $A(t)$). The total mass-loss rate is $\dot m_{\mathrm{total}}=\dot m''_{\mathrm{local}}\,A(t)$. For cone-calorimeter comparison, the convention is to report mass loss per unit *original* (gauge) area $A_0$:

$$
\dot m''_{\mathrm{gauge}} \;=\; \frac{\dot m_{\mathrm{total}}}{A_0}
\;=\; \dot m''_{\mathrm{local}}\,\frac{A(t)}{A_0}
\;=\; \dot m''_{\mathrm{local}}\,\mathrm{law}(\bar\chi).
\tag{10.17}
$$

> **Implementation note.** In `src/Solver/solve.jl` (ALE post-processing) the solver stores `solution.mass_loss_rate[out_idx] = total_mass_flux + bottom_mass_flux` — the local per-current-area flux leaving the slab through **both** faces (the bottom term is the film/venting flux of a permeable back-face mass BC, and is zero for the default impermeable holder) — and the gauge flavor `solution.extras.MLR_gauge[out_idx] = (total_mass_flux + bottom_mass_flux) * A_t / mesh.cross_section_area_initial`, exactly (10.17). The top-face-only flux remains separately available as `solution.extras.surface_mass_flux`. It also records `cross_section_area_history[out_idx] = A_t`. For identity-law / `NoChi` runs $A(t)\equiv A_0$ and the two flavors coincide. This gauge normalization is what makes a shrinking or intumescing sample's predicted MLR directly comparable to a cone measurement, which always reports per nominal specimen area. The energy-balance and conservation diagnostics that use the trapezoidal time-integral of $A$ for ALE efflux are described in §16, *Conservation Diagnostics*.

---

## 10.9 State-vector layout for variable cross-section

In the fully-featured ALE + `WithChi` mode the flat state vector is block-major:

$$
\mathbf{u} = \big[\underbrace{T_1,\dots,T_n}_{n},\;
\underbrace{\xi_{1,1},\dots,\xi_{1,n},\dots,\xi_{N_C,n}}_{N_C n},\;
\underbrace{z_1,\dots,z_{n_{\mathrm{nodes}}}}_{n_{\mathrm{nodes}}},\;
\underbrace{\chi_1,\dots,\chi_n}_{n}\big],
\tag{10.18}
$$

with total length $N_{\mathrm{state}} = n(1+N_C) + \delta_{\mathrm{ALE}}\,n_{\mathrm{nodes}} + \delta_{\mathrm{WithChi}}\,n$. The $z$ block is present only in ALE mode and the $\chi$ block only in `WithChi` mode. For example, with $N_C=2$ components, $n=5$ cells (hence $n_{\mathrm{nodes}}=6$), the ALE+WithChi length is $5\cdot 3 + 6 + 5 = 26$, which the test `test/geometry/test_variable_cross_section.jl:14` asserts as `5*(2+3)+1`; the same problem with `lateral=false` (NoChi) has length $5\cdot 3 + 6 = 21$ asserted as `5*3+6` (line 31). The block accessors `temperature_range`, `concentration_range`, `z_range`, `chi_range` (`src/Problem/state_layout.jl`) return empty ranges for absent blocks, and the residual's `_pack_chi_derivative_dispatched!` is a no-op when `chi_range` is empty. The `WithChi`/`NoChi` axis is one of the four bundled `SimulationMode` traits; the detailed layout descriptor `StateLayout` and its role in the structured Jacobian are covered in §13 and §15.1.

---

## 10.10 Worked sign check

To verify the sign conventions end to end, consider a single charring step `virgin → ν_char·char + ν_gas·gas` in one cell, with $N_C=3$ ordered (virgin, char, gas), swelling factors $\gamma=(1,1,0)$, and reaction rate $r>0$ (kg·m⁻³·s⁻¹). The reaction species sources are (from §6.5: reactant $-1$, products $+\nu$):

$$
\left.\frac{d\xi_{\mathrm{virgin}}}{dt}\right|_{\mathrm{rxn}}=-r,
\quad
\left.\frac{d\xi_{\mathrm{char}}}{dt}\right|_{\mathrm{rxn}}=+\nu_{\mathrm{char}}\,r,
\quad
\left.\frac{d\xi_{\mathrm{gas}}}{dt}\right|_{\mathrm{rxn}}=+\nu_{\mathrm{gas}}\,r,
$$

with $\nu_{\mathrm{char}}+\nu_{\mathrm{gas}}=1$ (mass conservation, the build-time invariant of §6.5). Then from (10.2):

$$
\theta = \frac{1}{\rho_{\mathrm{virgin}}}(-r) + \frac{1}{\rho_{\mathrm{char}}}(+\nu_{\mathrm{char}}r) + 0
= r\!\left(\frac{\nu_{\mathrm{char}}}{\rho_{\mathrm{char}}} - \frac{1}{\rho_{\mathrm{virgin}}}\right).
$$

Char is typically far less dense than virgin ($\rho_{\mathrm{char}}\ll\rho_{\mathrm{virgin}}$) but is produced in smaller mass ($\nu_{\mathrm{char}}<1$). The sign of $\theta$ is set by the competition $\nu_{\mathrm{char}}/\rho_{\mathrm{char}}$ vs. $1/\rho_{\mathrm{virgin}}$: for a high-char-yield, very-low-density char the matrix can *swell* ($\theta>0$, the char occupies more pure-phase volume than the virgin it replaced); for a low-char-yield, dense char it *shrinks* ($\theta<0$). The gas, with $\gamma=0$, never contributes — it escapes through the pores. With $\theta<0$ and a mesh-bound char, the dilation correction (10.15) gives $-\xi_{\mathrm{char}}\theta>0$, *raising* the char concentration as the cell contracts, conserving char mass. This is the correct, physically consistent behavior the code produces. Through (10.13), $\theta<0$ in every cell drives a downward cumulative mesh velocity ($w<0$ for nodes above the bottom), so the surface node recedes — material burns away from the top while the substrate stays pinned.

---

## 10.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

All three reference solid-phase pyrolysis codes confront the same physics — a condensed phase that loses mass and changes volume — but each encodes the deformation differently. We map their notation to ours on first use.

### 10.11.1 ThermaKin2Ds — the single-reference-point deformation source

ThermaKin2Ds (Stoliarov & Lyon) advances mass concentrations on a *fixed* spatial grid and represents volume change through a **deformation source term** rather than a moving mesh. Its species balance in 1D Cartesian is, in ThermaKin notation,

$$
\frac{\partial \xi_i}{\partial t}
= \sum_{\alpha\beta} \theta_i^{\alpha\beta}\,r_{\alpha\beta}
\;-\; \xi_i\,\frac{1}{\bar\rho}\frac{\partial \bar\rho}{\partial t}
\;-\; \frac{\partial J_i}{\partial x},
$$

where ThermaKin's $\theta_i^{\alpha\beta}$ is a **stoichiometric coefficient** (mapped to our $\nu_{r,j}$; recall the binding resolution A2 that reserves $\theta$ for strain rate), $r_{\alpha\beta}$ a reaction rate (our $r_r$), and $\bar\rho$ the mixture density. The middle term $-\xi_i\,\bar\rho^{-1}\,\partial_t\bar\rho$ is ThermaKin's deformation/dilation source: as the local mixture density changes, the fixed grid must "concentrate" or "dilute" the surviving species. Integrating $\bar\rho^{-1}\partial_x\rho$ (the single-reference-point construction) yields the displacement of material relative to the fixed grid.

The Pyrolysis.jl equivalent is split into two pieces that together reproduce ThermaKin's source: (i) the **mesh moves** with velocity (10.13), built from $\theta_i$ of (10.2), so that the grid follows the deformation instead of staying fixed; and (ii) the **dilation correction** $-\xi_j\theta_i$ of (10.15) provides the concentration rescaling for mesh-bound species. Pyrolysis.jl's $\theta=\sum_j\gamma_j\rho_j^{-1}\,d\xi_j/dt$ is the additive-volume analog of ThermaKin's $-\bar\rho^{-1}\partial_t\bar\rho$: both express that the matrix volume tracks the per-phase pure-volume balance, and indeed $\bar\rho^{-1}\partial_t\bar\rho = \sum_j \rho_j^{-1}\partial_t\xi_j$ for an additive-volume mixture (the $\gamma_j$ generalize this to selective swelling). The key structural difference is **ALE vs. fixed-grid-with-source**: ThermaKin keeps $x$ fixed and carries the deformation algebraically; Pyrolysis.jl moves the nodes and recomputes geometry. The two are mathematically equivalent in the continuum limit, but the ALE formulation keeps the cell-centered states co-located with the material they describe. ThermaKin treats the cross-section as fixed in its core 1D form; Pyrolysis.jl's lateral-shrinkage law (§10.5) is an *added* capability with no direct ThermaKin counterpart.

### 10.11.2 Gpyro — solid-fraction / mechanism-selection swelling

Gpyro (Lautenberger & Fernandez-Pello) solves control-volume conservation of $\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, gas mass, and energy on a variable-$\Delta z$ grid. Volume change is handled by allowing each control volume's thickness $\Delta z$ to change so that the **mass-weighted bulk density** stays consistent with the prescribed solid-fraction behavior of each species. Gpyro lets the user pick, per condensed species, how its volume responds: a species can be **non-shrinking** (its contribution to $\Delta z$ is held), **shrinking** (its volume is removed as it is consumed), or **swelling/intumescent** (it adds volume). In Gpyro's notation the bulk density $\bar\rho$ and per-species solid fraction control the $\Delta z$ update; conversion $\alpha_{\mathrm{conv}}$ tracks reaction extent.

The Pyrolysis.jl mapping is direct: Gpyro's per-species shrink/non-shrink/swell selection is our **swelling factor** $\gamma_j$ — $\gamma_j=1$ reproduces a shrinking species, $\gamma_j=0$ a non-shrinking one, and $\gamma_j>0$ on a gas reproduces an intumescent (volume-adding) species. Gpyro's variable $\Delta z$ is our ALE $\Delta z_i=V_i/A$ driven by the integrated $\theta_{L,i}$. Gpyro's conversion variable $\alpha_{\mathrm{conv}}$ is the close analog of our pyrolysis progress $\chi$, but Pyrolysis.jl uses $\chi$ for a *different* purpose — driving the *lateral* (cross-sectional) law (10.8) — whereas Gpyro's swelling acts on the axial $\Delta z$. In effect, Pyrolysis.jl has two channels (axial via $\theta$, radial via $\mathrm{law}(\bar\chi)$) where Gpyro's 1D model has only the axial one; the radial channel is Pyrolysis.jl's surrogate for the transverse swelling Gpyro would need a multi-D model to resolve.

### 10.11.3 FDS — the max/sum density-ratio swell/shrink rule

The FDS solid-phase pyrolysis model (McGrattan et al., *FDS Technical Reference Guide*, Solid Phase chapter) updates the layer thickness using a **density-ratio** rule. Each material component $\alpha$ has a reference (initial) density; as reactions change the component densities $\rho_{s,\alpha}$, FDS adjusts the cell thickness $\delta$ by the ratio of current to initial total density. FDS provides two limiting behaviors via a user flag:

- **Shrinking / "sum" behavior:** the layer thickness scales with the *sum* of component density ratios, so the cell shrinks in proportion to total mass lost (the standard charring-solid recession).
- **Swelling / "max" behavior:** the thickness scales with the *maximum* density ratio (the most-expanding component dominates), capturing intumescence where one product inflates the layer.

In FDS notation (their $\rho_{s,\alpha}$ for solid component densities, $X_\alpha$ for volume fractions), the thickness factor is built from $\sum_\alpha \rho_{s,\alpha}/\rho_{s,\alpha}^{0}$ (shrink) or a max-based variant (swell). Mapping to ours: FDS's $X_\alpha$ is our volume fraction $v_j$, and FDS's $\kappa_s$ (solid absorption) maps to our $\alpha$ (overload A5). The Pyrolysis.jl axial volume change (10.2) is the differential, additive-volume analog of FDS's *sum* rule: $\theta=\sum_j\gamma_j\rho_j^{-1}d\xi_j/dt$ is exactly a density-weighted sum of per-component mass-change rates, and integrating it reproduces a thickness that scales with total mass lost. The FDS *max* (swelling) behavior has no exact Pyrolysis.jl counterpart in the axial channel — Pyrolysis.jl represents swell additively (a $\gamma_g>0$ gas adds a positive term to the *sum*), not by a max selection. Radial/lateral swell is instead handled by the explicit $\mathrm{law}(\bar\chi)$ (§10.5), which FDS does not have (FDS keeps a fixed surface area and changes only thickness). A practical implication: FDS's max-rule can produce abrupt, single-component-dominated swelling, whereas Pyrolysis.jl's additive $\theta$ and smooth $\mathrm{law}$ give continuously differentiable deformation — important for the AD-based Jacobian and sensitivity machinery (§15, §17).

### 10.11.4 Summary of mappings

| Concept | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| Axial deformation mechanism | ALE mesh motion, $w=\int\theta_L\,dz$ | Fixed grid + deformation source $-\xi_i\bar\rho^{-1}\partial_t\bar\rho$ | Variable $\Delta z$ control volumes | Density-ratio thickness update |
| Per-species volume response | $\gamma_j\in[0,1]$ | implicit via $\bar\rho$ | shrink/non-shrink/swell flag | shrink (sum) / swell (max) flag |
| Strain rate | $\theta=\sum_j\gamma_j\rho_j^{-1}d\xi_j/dt$ | $\bar\rho^{-1}\partial_t\bar\rho$ (sign-equiv.) | bulk-density consistency | $\sum_\alpha\rho_{s,\alpha}/\rho_{s,\alpha}^0$ |
| Lateral (cross-section) law | $A=A_0\,\mathrm{law}(\bar\chi)$ | none (fixed $A$) | none in 1D | none (fixed area) |
| Progress / conversion variable | $\chi$ (dry-gas, for $\mathrm{law}$) | reaction extent | $\alpha_{\mathrm{conv}}$ | density ratio |
| Concentration rescaling under shrink | dilation correction $-\xi_j\theta$ | deformation source term | $\Delta z$ rebalance | thickness rebalance |

---

## 10.12 Limitations and open issues

The following limitations are documented in the source and KB notes and should guide use and future work:

1. **Density time-averaging in $\theta$.** Equation (10.2) uses $\rho_j(T)$ evaluated at the *current* cell temperature, and the density-rate term $\xi_j\,d(1/\rho_j)/dt$ is neglected (§10.1.1). For materials with strongly temperature-dependent pure-phase densities and steep thermal gradients this is an approximation; whether an implicit or time-averaged density would materially change $\theta$ has not been quantified.

2. **Empirical lateral law.** $\mathrm{law}(\bar\chi)$ is user-specified with no first-principles closure (§10.5.3). A percolation- or homogenization-based predictor of the functional form would remove a calibration degree of freedom, especially for intumescent coatings whose expansion factor can exceed an order of magnitude. The law must also be *smooth*: `_compute_theta_A` differentiates it by `ForwardDiff`, so a non-smooth (e.g. piecewise-linear) law gives undefined $\theta_A$ at the kinks and degrades the Jacobian/sensitivity quality (see §17.7).

3. **Vanishing cross-section.** The radial strain rate (10.10) carries a $1/\mathrm{law}(\bar\chi)$ that diverges as the cross-section collapses; the code guards with `law(χ̄) > eps()` (returning $\theta_A=0$ below the floor), but there is no *physical* lower bound on shrinkage — a law allowed to reach zero is unphysical and the guard merely prevents a numerical blow-up. Users should design $\mathrm{law}$ so it stays bounded away from zero.

4. **$\bar\chi$ weighting choice.** The column average (10.7) is mass-weighted by *initial* dry-solid mass. Strongly non-uniform decomposition (e.g. a thin, fully-charred surface over an intact interior) is reduced to a single scalar $\bar\chi$, so the lateral law cannot represent depth-dependent cross-section change. A different weighting (or a depth-resolved area) would be required for that regime.

5. **Single-reactant kinetics.** Volume change inherits the single-reactant restriction of the kinetics (§6.5, §6.10): `species_source_terms!` enforces $N_{\mathrm{Rin}}=1$, so $\theta$ and $\chi$ are built from single-reactant sources only. Bimolecular or multi-reactant stoichiometry (e.g. char oxidation consuming both char and O₂) would require extending the source-term machinery before $\theta$/$\chi$ could account for it.

6. **Interaction with adaptive mesh refinement.** The $\chi$ block, the per-cell strain $\theta_i$, and $m_{\mathrm{dry},i}$ must all be carried through any AMR remap (surface-depletion merge, h-refinement). The energy-conserving merge and the mass-exact remap are detailed in §14, *Adaptive Mesh Refinement*; the consistency of $\chi$ accumulation across a merge (summing $m_{\mathrm{dry}}$ denominators) is a subtle invariant that must be preserved, and the interaction of the depletion callback with the pluggable Jacobian (which falls back to colored finite differences) is an open coupling discussed in §14 and §15.

7. **Substring-based water-vapor exclusion.** The dry-gas index set relies on a name substring match (§10.3.3), which is fragile to naming. A type-level or explicit-flag mechanism for "this gas is moisture-derived" would be more robust than the current `"H2O"`/`"water"` heuristic.

8. **Axial/radial split is sequential, not unified.** As noted in the ALE KB, the axial (mesh-motion) and radial (cross-section) channels are computed in sequence with $\theta_A$ subtracted from each cell's $\theta_i$; a fully unified space-time ALE that handles both simultaneously (and a higher-order-in-time mesh update than the forward-Euler node advance) remains future work (see §11.11).

9. **Mesh motion is solid motion (solid-Lagrangian ALE).** The volume-change pipeline of this chapter is what *defines* the mesh velocity: $w$ is the prefix sum of the reaction-driven strain rate $\theta$ (§11.3), and condensed phases ride the mesh with no discrete advection operator (§10.7, §11.1.2). Mesh motion decoupled from the solid material velocity — remeshing, smoothing, prescribed motion, AMR-in-ALE — is not representable without adding a condensed-phase advection operator; see §11.11 item 7 for the full statement of this capability boundary and its consequence for MMS verification (the spatially-uniform-solid rule under `w_override`).
