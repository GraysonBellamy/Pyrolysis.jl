# 10. Volume Change, Swelling, Intumescence, and Lateral Shrinkage

This chapter derives the volume-change model of Pyrolysis.jl: how reaction-driven mass exchange between condensed and gaseous phases produces a per-cell volumetric strain rate $\theta$, how that strain enters the balance rows that fix the node velocities of a moving mesh, and how a column-global *lateral* shrinkage / swelling / intumescence law deforms the cross-sectional area $A(t)$ through the column's conversion $\bar\chi$, which a per-cell *conversion-history pair* carries. We give the continuous balance, the discrete rows as implemented, why integrating the component masses makes the dilation of every concentration implicit, and a comparison to the deformation treatments in ThermaKin2Ds, Gpyro, and FDS. Every equation here is grounded in the source: `src/Physics/volume_change.jl` (the strain rate), `src/Closures/balance.jl` (the balance rows), `src/Closures/lateral.jl` (the variable cross section), `src/Kernels/relative_flux.jl` (the history pair on a mesh that moves relative to the solid), `src/Assembly/residual.jl` (the row assembly), `src/Integration/setup.jl` (the initial state and the validation), and `src/Materials/materials.jl` (the dry-gas set).

The mesh-motion machinery (the node positions as differential unknowns, the balance rows that fix the node velocities, the geometric volume from the node positions, the equidistribution correction and the relative fluxes) is developed in Chapter 11, *The ALE Framework*; here we develop the *constitutive* content — the strain rates and the cross-section law that feed Chapter 11 — and the bookkeeping of the history pair. The state layout that carries the history pair and the cross-section rows is detailed in §15.1, and the Jacobian of the two dense cross-section rows in §15.9.

**Index convention for this chapter.** Following the binding resolution in §2 (overload I1), $i$ is the **cell index** ($1\le i\le n$), $j$ is the **component/species index** ($1\le j\le N_C$), and $r$ is the **reaction index**. The symbol $\chi$ is *always* pyrolysis progress (overload A1); the cross-section ratio is $A/A_0=\mathrm{law}(\bar\chi)$ and is never written as $\chi$. The symbol $\theta$ is *always* a volumetric strain rate (overload A2); stoichiometry is $\nu_{r,j}$. The second member of the conversion-history pair is written $m_{\mathrm{dry},i}$ (the code calls it `M`, a letter the nomenclature reserves for molar mass).

---

## 10.1 The volumetric strain rate $\theta$

### 10.1.1 Continuous derivation from the mixture continuity equation

Consider a material control volume $V$ whose condensed-phase matrix (solids and liquids) forms a skeleton through which gas can percolate. The matrix occupies a volume that is the sum of the partial volumes of each matrix-forming component. If component $j$ has bulk (pure-phase) density $\rho_j$ and mass concentration $\xi_j$ (kg of $j$ per m³ of bulk material, see §2), then the mass of $j$ in the cell is $m_j=\xi_j V$, and the partial (pure-phase) volume it occupies is $m_j/\rho_j=\xi_j V/\rho_j$.

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

> **Implementation note.** $\theta$ is computed by `volume_change_rate(dξ, T, material)` in `src/Physics/volume_change.jl` (inline methods for $N_C\in\{1,2,3\}$ and a `@generated` method for $N_C>3$, so the loop over the heterogeneous component tuple allocates nothing). The loop is
> ```julia
> θ += γⱼ * dξ[j] / ρⱼ        # ρⱼ = comps[j].density(T)
> ```
> evaluated **only when `γⱼ > 0`**: components with $\gamma_j=0$ are skipped, so no density is evaluated and the term is exactly zero. There is **no minus sign** in the code: the sign of $\theta$ is carried entirely by the sign of $d\xi_j/dt$, which is negative for a depleting reactant (the species source consumes the reactant with coefficient $-1$; see §6.5). The density $\rho_j(T)$ is the *bulk* density property evaluated at the cell's trial temperature $T_{\mathrm{trial}}$ (the cell temperature smoothly clamped into the property bracket $[200, 3000]$ K, the identity on every physical state) — not the intrinsic/skeletal density $\rho_{i,j}$, which is reserved for the porosity (§5.3 and §2, overload A7).

### 10.1.2 What $d\xi_j/dt$ is in (10.2)

A subtle and load-bearing point: the $d\xi_j/dt$ fed into (10.2) is the **reaction source only**, $S_{j}=\sum_r\nu_{r,j}r_r$, not the full rate of change of the concentration (which also includes gas transport, the fluxes relative to a moving mesh, and the change of the cell volume itself). This is correct because $\theta$ measures the rate at which mass is *converted between phases by chemistry* — the only process that changes the additive-volume balance (10.1). Transport of an already-gaseous species through the pores does not change the matrix volume, and the concentration change that the volume change itself produces is the consequence of $\theta$, not a cause of it (§10.7). Including them would double-count the volume change that the mesh motion is itself representing.

> **Implementation note.** `closure_rows!` (`src/Closures/balance.jl`) runs after every kernel and before the external sources are added (`add_sources!`), and reads the cell mass sources `c.S_mass` [kg/s], which only the reaction kernel writes (the transport kernels write face fluxes). For each cell it forms `dξ = S_mass[:, i] / V_eval[i]`, with $\tilde V_i$ the smoothly floored cell volume, and sets `c.θ[i] = volume_change_rate(dξ, ds.T_trial[i], material)`. So $\theta_i$ is the strain rate of the reactions alone, exactly (10.2) with $d\xi_j/dt = S_{j,i}$; a manufactured mass source of a verification run does not move the mesh.

### 10.1.3 Absolute cell volume rate

The absolute rate of change of a cell's volume follows immediately:

$$
\frac{dV_i}{dt} \;=\; V_i\,\theta_i.
\tag{10.3}
$$

> **Implementation note.** No routine materializes $dV_i/dt$. The cell volume is computed from the node positions (and, under a lateral law, the area unknown) at every residual evaluation (§11.4.1), and the balance rows (§10.6.2, §11.3.2) move the nodes so that (10.3) holds semi-discretely; when the mesh moves relative to the solid under a correction, $dV_i/dt = \theta_i V_i + s_i$ with $s_i=\sum_f s_{i,f}A_f(c\cdot n)_f$ the volume the node corrections sweep (§11.7).

---

## 10.2 Swelling-factor semantics

The swelling factor $\gamma_j$ is a per-component scalar field on every `AbstractComponent`. Its defaults encode the standard physics:

| Phase | Default $\gamma$ | Meaning | Effect on $\theta$ |
|---|---|---|---|
| Solid | $1.0$ | Skeleton-forming; depletion frees pure-phase volume | $\theta<0$ when solid is consumed |
| Liquid | $0.0$ | Fills pore space; default treats it as volume-neutral | none unless overridden |
| Gas | $0.0$ | Escapes through pores; does not set matrix dimensions | none |

> **Implementation note.** Defaults are set in the convenience constructors of `src/Materials/components.jl`: `SolidComponent(...; γ::Real = 1.0, ...)`, `GaseousComponent(...; γ::Real = 0.0, ...)`, and `LiquidComponent(...; γ::Real = 0.0, ...)`. The field is `swelling_factor::Float64` on each struct, validated to $\gamma\in[0,1]$ at construction (an `ArgumentError` names the component otherwise). The factor has one meaning everywhere: the strain rate (10.2), the condensed transport weights and the thermal mixing volume fractions that weight the conductivity and the surface emissivity (§5.1.3) all use the same $\gamma_j\,\xi_j/\rho_j$ weighting.

**Modeling guidance.**

- **Non-charring shrinkage (e.g. a polymer converting to gas):** keep the virgin solid at $\gamma=1$ and the gas product at $\gamma=0$. As the solid converts to gas, $\theta<0$ and the column recedes — the classic shrinking-slab picture.
- **Zero-swelling char:** `γ=0` removes the char's direct strain and its matrix property weight. Its skeletal occupancy, mass and heat storage remain. Such a phase must have zero transport coefficients; a state with no γ-weighted matrix still needs a declared endpoint whose skeletal-occupancy domain includes that state. This choice does not make residual char an unobstructed void.

- **Intumescence (gas-producing expansion):** assign the trapped/blowing gas a swelling factor $\gamma_g\in(0,1]$. A positive $d\xi_g/dt$ (gas being produced) with $\gamma_g>0$ contributes $\theta>0$, i.e. the matrix *inflates*. Combined with a swelling lateral law (§10.5) this reproduces the order-of-magnitude expansion of intumescent coatings. Axial intumescence (through-thickness swell) is captured by $\theta>0$ in (10.2); radial intumescence (lateral swell of the column) is captured by $\mathrm{law}(\bar\chi)>1$ in (10.8). A factor $\gamma>1$ is not representable.

Because $\gamma_j$ is a free per-component parameter, the model spans the full range from pure shrinkage ($\gamma=1$ solid, $\gamma=0$ gas) through volume-neutral ($\gamma=0$ everywhere) to intumescent ($\gamma_g>0$).

---

## 10.3 The pyrolysis progress and the conversion-history pair

### 10.3.1 Definition and governing equations

To make the *lateral* (cross-section) deformation depend on how far decomposition has advanced, Pyrolysis.jl defines a per-cell **pyrolysis progress** $\chi_i$ (dimensionless): the cumulative mass of *dry* pyrolysis gas produced by the reactions of cell $i$, normalized by the cell's initial dry-solid mass,

$$
\chi_i \;=\; \frac{G_i}{m_{\mathrm{dry},i}},
\qquad
\boxed{\;\frac{dG_i}{dt} \;=\; S_{\mathrm{dry},i}\,V_i\;},
\qquad
S_{\mathrm{dry},i} \;=\; \sum_{j\in\mathrm{dry\ gas}}S_{j,i},
\tag{10.4}
$$

where $G_i$ [kg] is the cumulative dry gas produced in cell $i$, $S_{\mathrm{dry},i}$ (units $\mathrm{kg\,m^{-3}\,s^{-1}}$) is the **dry-gas production rate** — the sum of the reaction sources over the components flagged as dry pyrolysis gas — $V_i$ is the current cell volume, and $m_{\mathrm{dry},i}$ is the cell's initial dry-solid mass (§10.3.2). The product $S_{\mathrm{dry},i}V_i$ is the dry-gas mass-production *rate* (kg/s); integrating (10.4) gives

$$
\chi_i(t) \;=\; \frac{1}{m_{\mathrm{dry},i}}\int_0^t S_{\mathrm{dry},i}(t')\,V_i(t')\,dt'
\;=\; \frac{m_{\mathrm{produced,dry},i}(t)}{m_{\mathrm{dry},i}},
\tag{10.5}
$$

i.e. exactly the fractional dry-gas mass production. By construction $\chi_i$ starts at $0$ and, for a fully decomposing solid that converts entirely to dry gas, tends to $1$; with char yield it saturates below $1$.

> **Implementation note.** The pair $(G_i, m_{\mathrm{dry},i})$ — the *conversion-history pair* — is carried per cell as two differential unknowns of the state, the cell block's carried field `:history` (`history_index(layout, i, 1)` for $G_i$, `history_index(layout, i, 2)` for $m_{\mathrm{dry},i}$; `src/States/layout.jl`). Their sources are written by `closure_rows!` through `_history_sources!` (`src/Closures/lateral.jl`) from the reaction mass sources:
> ```julia
> c.S_carried[1, i] = Σ_{j ∈ dry gas} c.S_mass[j, i]   # dG/dt source [kg/s]
> c.S_carried[2, i] = 0                                 # m_dry has no source
> ```
> The row assembly adds the face fluxes of the pair (§10.3.4) and writes the two rows (`_history_rows!`, `src/Assembly/residual.jl`). The cell volume enters through the reaction kernel's cell source (rate × $\tilde V_i$), so $G_i$ accumulates the true mass production while the cell changes volume, consistent with the integral form (10.5). $\chi_i$ itself is never formed: no per-cell division by $m_{\mathrm{dry},i}$ takes place, and only the column ratio of §10.4 enters the model. A material without a lateral shrinkage law carries no pair (the layout has no `:history` field), so the pair costs nothing when it is not needed.

### 10.3.2 The initial dry-solid mass $m_{\mathrm{dry},i}$

The second member of the pair is the **solid** mass present in cell $i$ at $t=0$:

$$
m_{\mathrm{dry},i} \;=\; \sum_{j\in\mathrm{solid}} m_{j,i}(0) \;=\; V_i(0)\sum_{j\in\mathrm{solid}}\xi_{j,i}(0).
\tag{10.6}
$$

> **Implementation note.** `setup` (`src/Integration/setup.jl`) evaluates the initial fields at the cell centres (§12.8), forms the initial inventories $m_{j,i}(0)=\xi_{j,i}(0)V_i(0)$, sums them over the components for which `is_solid` is true, and writes the result into the second history row of the initial state; $G_i$ starts at zero. The sum is over **solids only** — liquids (e.g. bound moisture) are *not* counted in the denominator, even though they are matrix-forming. This is deliberate: $\chi$ is a *dry-decomposition* progress variable, and water mass is excluded from both numerator (dry gas only) and denominator (solids only). The row of $m_{\mathrm{dry},i}$ has no source; it changes only through the relative fluxes, when the mesh moves relative to the solid (§10.3.4), and at a deactivation event, which adds the dead cell's pair to the receiver's (Chapter 14). Both conserve the column total $M_{\mathrm{dry}}=\sum_i m_{\mathrm{dry},i}$ over the live cells, which is therefore fixed for the whole solve (until the last live cell goes with no receiver, taking its pair with it).

### 10.3.3 Dry-gas identification and the water-vapor caveat

The index set $\{j\in\mathrm{dry\ gas}\}$ collects all GAS-phase components *except* water vapor / moisture-derived gases. Phase is a compile-time (type-level) property, but the water-vapor exclusion is done by **name matching**:

> **Implementation note.** `dry_pyrolysis_gas_indices(material)` (`src/Materials/materials.jl`) recurses over the component tuple at compile time (`_collect_dry_pyrolysis_gas_indices`) and includes index $j$ iff `is_gas(comp) && !_is_water_vapor_name(comp.name)`. The name filter `_is_water_vapor_name` is a **substring match**: a name is treated as water vapor if it contains `"H2O"` (case-sensitive), if its lowercase form contains `"water"`, or if its lowercase form is exactly `"h2o"`. `setup` turns the result into the flag tuple `dry_gas::NTuple{NC,Bool}` of the `LateralShrinkage` it builds, once. **Caveat:** because identification is by substring, a dry-gas species accidentally named (say) `:water_gas` would be misclassified as moisture and excluded from $\chi$. Users modeling moisture must name the water-vapor component so it matches (e.g. `:water_vapor` or `:H2O`), and must avoid those substrings in genuine dry-gas names.

### 10.3.4 The history pair on a moving mesh and at deactivation

The pair is *extensive*: it belongs to the solid of a cell, and wherever that solid goes the pair goes with it.

- **A mesh that moves relative to the solid.** Under the equidistribution correction the mesh velocity is $w = v_s + c$ and every inventory crosses the faces at the relative velocity $-c$ (§11.6). The pair moves with the solid: through each face its flux is the solid mass flux $F_{\mathrm{solid}}=\sum_{j\in\mathrm{solid}}F_j$ times the per-solid-mass value of the upwind reconstruction,
  $$
  F_{G,f} = F_{\mathrm{solid},f}\,\Big(\frac{G}{m_s}\Big)_f,
  \qquad
  F_{m_{\mathrm{dry}},f} = F_{\mathrm{solid},f}\,\Big(\frac{m_{\mathrm{dry}}}{m_s}\Big)_f,
  $$
  with $m_s$ a cell's solid mass smoothly floored at its mass scale (`RelativeFluxKernel`, `src/Kernels/relative_flux.jl`; §11.6.4). The face values are second-order van Albada reconstructions; unlike the concentrations they are not inventories and carry no positivity bound. The correction vanishes at the end nodes (§11.7.5), so the boundary faces carry no pair and the interior fluxes telescope: the column sums are unchanged.
- **Deactivation.** When a cell deactivates, its pair is added to the receiver's and zeroed (§14.3.3). The column sums $\sum_i G_i$ and $\sum_i m_{\mathrm{dry},i}$ over the live cells, hence $\bar\chi$ and the area, are unchanged by the event.

---

## 10.4 Mass-weighted column average $\bar\chi$

The lateral-shrinkage law acts on a single, *column-global* progress scalar, the **mass-weighted column average**, the column's conversion:

$$
\boxed{\;\bar\chi \;=\; \frac{\displaystyle\sum_{i} m_{\mathrm{dry},i}\,\chi_i}
{\displaystyle\sum_{i} m_{\mathrm{dry},i}}
\;=\; \frac{\displaystyle\sum_{i} G_i}{\displaystyle\sum_{i} m_{\mathrm{dry},i}}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i} G_i\;},
\qquad
M_{\mathrm{dry}}=\sum_{i} m_{\mathrm{dry},i},
\tag{10.7}
$$

the sums running over the live cells. Weighting by the *initial* dry-solid mass (rather than by volume, or uniformly) means $\bar\chi$ is the fraction of the column's *total decomposable solid mass* that has been produced as dry gas — a physically meaningful, conserved-denominator average. Because $M_{\mathrm{dry}}$ is fixed (§10.3.2), $\bar\chi$ moves monotonically as the column decomposes, and lies in $[0,1]$ when the dry gas comes from the solid.

> **Implementation note.** `chi_bar(lateral, ΣG, ΣM) = ΣG / smooth_floor(ΣM, M_floor)` (`src/Closures/lateral.jl`), with the floor $10^{-8}M_{\mathrm{dry}}(0)$ (the identity once the denominator is twice the floor), so the ratio stays finite with a finite derivative on any trial state. The sums are formed by `_lateral_rows!` (`src/Assembly/residual.jl`) over the cells the layout marks live, in the element type of the state, so under automatic differentiation they carry the derivatives of every history unknown. Before it reaches the law the conversion is smoothly clamped,
> $\bar\chi_c = \mathrm{smooth\_clamp}\big(\bar\chi,\,-10^{-3},\,1+10^{-3},\,10^{-3}\big)$,
> the identity on $[0,1]$, constant outside $[-10^{-3}, 1+10^{-3}]$ and $C^1$ in the joins (`LATERAL_CLAMP_WIDTH`): every physical conversion passes unchanged, and the law's argument is bounded on trial states.

---

## 10.5 The lateral-shrinkage law $A(t)=A_0\,\mathrm{law}(\bar\chi)$

### 10.5.1 Definition

The column cross-sectional area evolves as a prescribed function of the column's conversion:

$$
\boxed{\;A(t) \;=\; A_0\,\mathrm{law}(\bar\chi_c)\;},
\tag{10.8}
$$

where $A_0$ is the cross-section measure of the mesh, $\bar\chi_c$ the clamped conversion of §10.4, and $\mathrm{law}$ a user-supplied smooth function with $\mathrm{law}(0)=1$, so that $A(0)=A_0$. The area $A(t)$ is uniform along $z$ (the column has one cross-section at any instant) but time-varying. In the discrete system the area is an **algebraic unknown** of the state, and (10.8) is its row, the *area row*,

$$
0 \;=\; A - A_0\,\mathrm{law}(\bar\chi_c),
$$

which holds exactly at the consistent start and to the integrator's tolerance at every accepted state (§10.6.3).

> **Implementation note.** The law is the material's `lateral_shrinkage_law` (`Material(...; lateral_shrinkage_law = nothing)`, `src/Materials/materials.jl`), a `Function` of one argument or `nothing`, which means a constant cross section and no cross-section rows at all. $A_0$ is the cross-section measure of the mesh (the `cross_section` keyword of `create_uniform_mesh` and its siblings, 1 m² by default on a planar line). `setup` builds `LateralShrinkage(material, A₀, M₀)` (`src/Closures/lateral.jl`) with $M_0=M_{\mathrm{dry}}(0)$, which refuses with an `ArgumentError`:
>
> - a law that is not positive and finite at every one of 201 equally spaced points of $[-10^{-3},\,1+10^{-3}]$, the clamp band of $\bar\chi_c$;
> - a law with $|\mathrm{law}(0)-1|>10^{-12}$;
> - a column without dry solid ($M_0 = 0$), whose conversion is undefined.
>
> `setup` also refuses a lateral law without `use_ale = true` ("material … has a lateral_shrinkage_law, which needs a moving mesh: solve with use_ale = true"), on cylindrical or spherical coordinates ("a lateral shrinkage law needs planar coordinates"), and on a mesh of more than one dimension. The area row's residual is the first component of `lateral_rows(lateral, s)`, a function of six column sums (§10.6.1).

### 10.5.2 Canonical laws

| Law | $\mathrm{law}(\bar\chi)$ | Behavior |
|---|---|---|
| Identity (default, `nothing`) | $1$ | Constant cross-section $A=A_0$ (standard 1D pyrolysis) |
| Linear shrink | $1-k\bar\chi$, $0<k<1$ | Area contracts linearly with extent of decomposition |
| Linear swell | $1+k\bar\chi$, $k>0$ | Area expands (intumescence) |
| Smooth (sigmoid) | $1-k\big[\sigma(s(\bar\chi-\bar\chi_0))-\sigma(-s\bar\chi_0)\big]$, $\sigma(x)=1/(1+e^{-x})$ | Onset-localized transition; the offset makes $\mathrm{law}(0)=1$ |

The linear shrink law must keep $k<1/(1+10^{-3})$ to stay positive on the clamp band. The unit test of the cross-section rows (`test/closures/test_variable_cross_section.jl`) uses $\mathrm{law}(\bar\chi)=1-0.2\,\bar\chi$: with every cell's $G_i$ set to half its $m_{\mathrm{dry},i}$ ($\bar\chi=0.5$) the area row is not satisfied at $A_0$ and is satisfied to $10^{-14}A_0$ at $A=0.9\,A_0$. The verification case V-CPL-03 checks the law against an exact endpoint: with the chemistry virgin $\to$ 0.5 char + 0.5 dry gas, a char of half the virgin density (so $V_f=V_0$) and $\mathrm{law}(\bar\chi)=1-0.5\,\bar\chi$, complete conversion gives $\bar\chi_f=1/2$, $A_f/A_0=0.75$ and $L_f/L_0=(V_f/V_0)/(A_f/A_0)=4/3$, which the case asserts to $3\times10^{-6}$.

A complete solve with the linear shrink law, a 1 cm column of a solid that converts entirely to dry gas, on a moving mesh:

```julia
using Pyrolysis

solid = SolidComponent(:solid; ρ = 1000.0, φ=.3, c = 1000.0, k = 0.2, ε = 0.9, λ = 1e-6)
gas = GaseousComponent(:gas; M = 0.050, c = 1000.0, k = 0.02,)
shrink = Material(
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6,κ=0.,max_skeletal_occupancy=.001)),

    name = :shrink,
    components = (solid, gas),
    reactions = (Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
    lateral_shrinkage_law = χ -> 1.0 - 0.2χ,
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, 10),
    material = shrink,
    bc_set = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(3e4), bottom = AdiabaticBC()),
        mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
    ),
    ξ_initial = [z -> 1000.0, z -> 0.0],
    tspan = (0.0, 30.0),
    output_times = [0.0, 10.0, 20.0, 30.0],
)
sol = solve(problem; use_ale = true)

A = sol.extras.cross_section_area_history            # A(t) [m²]
χ̄ = (1 .- A ./ A[1]) ./ 0.2                          # the conversion the law inverts to
ratio = sol.extras.MLR_gauge ./ sol.mass_loss_rate   # A/A₀ after t = 0 (§10.8.2)
```

It returns `:Success` with $A = 1,\ 0.99288,\ 0.97408,\ 0.95451$ m² at $t = 0, 10, 20, 30$ s ($\bar\chi = 0,\ 0.0356,\ 0.1296,\ 0.2275$) while the thickness falls from 10 mm to 9.713, 8.936 and 8.094 mm. The same problem solved without `use_ale = true` is refused with the `ArgumentError` quoted in §10.5.1.

### 10.5.3 First-principles status

The lateral law is **empirical / phenomenological** — a calibration handle, not a derived constitutive relation. There is no first-principles model (e.g. a percolation- or homogenization-based predictor) in the code that selects the functional form of $\mathrm{law}(\bar\chi)$ from microstructure; the user supplies it. This is by design: a 1D solver cannot resolve the lateral (transverse) dynamics that actually set $A(t)$, so the law condenses that unresolved physics into a single $\bar\chi$-dependence. Its primary uses are (i) matching cone-calorimeter samples that visibly shrink or intumesce so the per-unit-*original-area* mass-loss-rate gauge is correct (§10.8.2), and (ii) order-of-magnitude intumescent-thickness studies. Section 10.12 lists the open questions this raises.

---

## 10.6 The lateral strain rate $\theta_A$ and the axial/lateral decomposition

![Initial, axially shrunken, laterally shrunken and swollen one-dimensional samples, showing volume as area times thickness and two MLR normalizations.](../../assets/figures/11_shrinkage.svg)

**Figure 11. Axial versus lateral volume change (schematic).** A planar cell has $V=A\Delta z$ and total strain rate $\theta=\theta_L+\theta_A$. The axial balance uses $\theta-\theta_A$, so lateral volume change is counted once. The drawn height stands for cross-sectional area, not a resolved lateral mesh. For the same total outward mass rate, `mass_loss_rate` divides by the current exposed area and `extras.MLR_gauge` by its initial area. Shrinking area makes the former larger than the latter.

### 10.6.1 Definition, chain-rule derivation, and the lateral-rate row

A time-varying cross-section contributes its own volumetric strain. Define the **column-global lateral (radial) volumetric strain rate**

$$
\theta_A \;\equiv\; \frac{1}{A}\frac{dA}{dt}.
\tag{10.9}
$$

Differentiating (10.8) and using the chain rule (with $A_0$ constant, and on $[0,1]$ where the clamp is the identity):

$$
\frac{dA}{dt} \;=\; A_0\,\mathrm{law}'(\bar\chi)\,\frac{d\bar\chi}{dt},
\qquad\Longrightarrow\qquad
\boxed{\;\theta_A \;=\; \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\,\frac{d\bar\chi}{dt}
\;=\; \frac{d}{dt}\log\mathrm{law}(\bar\chi)\;}.
\tag{10.10}
$$

The column progress rate follows from the history rows. With $\bar\chi=\sum G_i/\sum m_{\mathrm{dry},i}$ and $\sum_i\dot m_{\mathrm{dry},i}=0$ (§10.3.2), (10.4) gives

$$
\frac{d\bar\chi}{dt}
\;=\; \frac{\sum_i \dot G_i - \bar\chi\sum_i\dot m_{\mathrm{dry},i}}{M_{\mathrm{dry}}}
\;=\; \frac{1}{M_{\mathrm{dry}}}\sum_{i} S_{\mathrm{dry},i}\,V_i ,
\tag{10.11}
$$

a volume-weighted sum of the dry-gas source over the column divided by the total initial dry-solid mass (in per-cell form, $\sum_i m_{\mathrm{dry},i}\dot\chi_i/M_{\mathrm{dry}}$ with the $m_{\mathrm{dry},i}$ cancelling cell by cell).

$\theta_A$ is the second algebraic unknown of the cross section. Its row, the *lateral-rate row*, is the logarithmic rate of the area row along the column sums of the history rates:

$$
0 \;=\; \theta_A - \frac{d}{d\tau}\,\log\mathrm{law}\!\left(\bar\chi_c\!\left(\frac{\sum_i G_i+\tau\sum_i\dot G_i}{\lfloor\sum_i m_{\mathrm{dry},i}+\tau\sum_i\dot m_{\mathrm{dry},i}\rfloor}\right)\right)\bigg|_{\tau=0},
$$

with $\lfloor\cdot\rfloor$ the smooth floor of §10.4. On $[0,1]$ and above the floor this is (10.10)–(10.11); everywhere it is the exact time derivative of $\log$ of the area row's right-hand side, so on the solution, where both rows hold, the area's own rate is exactly $dA/dt=\theta_A A$, the definition (10.9).

> **Implementation note.** `_lateral_rate` (`src/Closures/lateral.jl`) takes this derivative by forward-mode automatic differentiation, `ForwardDiff.derivative(τ -> log(law(χ̄_c(ΣG + τΣĠ, ΣM + τΣṀ))), 0.0)`, so the user supplies only $\mathrm{law}$, never its derivative, and any differentiable law works. The rates $\sum\dot G$ and $\sum\dot m_{\mathrm{dry}}$ are the assembled history rows (sources plus face fluxes) of the same residual evaluation: `_lateral_rows!` runs after `_history_rows!` in `assemble_rows!` and reads them from `du`. No guard switches $\theta_A$ off: the validation of §10.5.1 keeps $\mathrm{law}$ positive on the whole clamp band, so $1/\mathrm{law}$ is bounded. Both rows depend on the state only through the six column sums $\boldsymbol s=(\sum G,\sum m_{\mathrm{dry}},\sum\dot G,\sum\dot m_{\mathrm{dry}},\theta_A,A)$ over the live cells (`lateral_rows`); they read every live cell's pair, the lateral one also its rates, and are the system's two dense rows, assembled through those sums (§15.9).

### 10.6.2 Axial residual strain and the balance rows

The per-cell strain rate $\theta_i$ from (10.2) is the *total* volumetric strain of the cell — it includes both the through-thickness (axial) deformation that the mesh tracks by moving nodes along $z$ **and** the lateral deformation that the cross-section law tracks by changing $A$. If the full $\theta_i$ drove the node velocities *and* $A(t)$ also changed, the lateral volume change would be counted twice. The resolution is to subtract the column-global lateral channel and drive the nodes with the **axial residual strain**

$$
\boxed{\;\theta_{L,i} \;=\; \theta_i - \theta_A\;}.
\tag{10.12}
$$

Here $\theta_A$ is the same for every cell (it is a column-global quantity), and $\theta_i$ is cell-local. The node velocities are algebraic unknowns held by the **balance rows** (`BalanceRowClosure`, §11.3.2): the back node is pinned, and each cell's volume balance at the solid velocity $v_s = w-c$ sits on the velocity row of its outer node,

$$
0 = w_1,
\qquad
0 = A_{i+1}\,(w-c)_{i+1} - A_i\,(w-c)_i - (\theta_i-\theta_A)\,V_i
\quad (i=1,\dots,n).
\tag{10.13}
$$

On a planar line every face has the area $A$ of the cross section, and the rows reduce to the recursion $v_{s,1}=0$, $v_{s,i+1}=v_{s,i}+(\theta_i-\theta_A)\,\Delta z_i$ with $\Delta z_i=V_i/A$: the solid velocity at node $i+1$ is the prefix sum $\sum_{k\le i}\theta_{L,k}\Delta z_k$, the substrate is immobile, and the surface moves at the total axial volume-change rate of the column. That the subtraction is exactly right follows from the volume of a planar cell, $V_i = A\,(z_{i+1}-z_i)$: with $c=0$,

$$
\frac{dV_i}{dt} = \frac{dA}{dt}\,\Delta z_i + A\,(w_{i+1}-w_i) = \theta_A V_i + (\theta_i-\theta_A)\,V_i = \theta_i V_i ,
$$

which is (10.3): the area carries the lateral part of the strain and the node motion the rest. The full development of (10.13) — the node positions as state ($dz_n/dt = w_n$), the geometric volume recomputed from them, the correction $c$ — is the subject of Chapter 11.

> **Implementation note.** `closure_rows!` fills the strain rates `c.θ` (§10.1.2) and then the velocity rows (`_velocity_rows!`, `src/Closures/balance.jl`), reading $\theta_A$ from the state (`ds.θ_A`, zero without a lateral law, when (10.13) is the constant-cross-section balance), the floored volumes and the face areas of the current geometry. Each row reads one cell and its two nodes (and, under a correction, what the corrections of those nodes read), so the Jacobian stays banded; the prefix sum is formed only at the consistent start (§10.6.3). The row node of a dead cell carries an identity row (Chapter 14).

### 10.6.3 The rows of one DAE, and the consistent start

The area, the lateral rate and the node velocities are algebraic unknowns; the masses, the energy states, the node positions and the history pairs are differential ones. Every row is evaluated in one residual pass at one state: the state mapping reads $A$ from the state and computes the cell volumes and face areas with it; the kernels run; `closure_rows!` computes the strain rates, the history sources and the velocity rows with the $\theta_A$ of the state; the row assembly writes the history rows and then the two cross-section rows from the column sums. No channel is computed from another in sequence: the integrator advances all rows together as one semi-explicit DAE with the mass matrix zero on the algebraic rows (the default `Rodas5P`, a Rosenbrock method, linearizes them in each stage with the colored Jacobian; Chapter 15). The DAE is of index one: given the differential unknowns, the cross-section rows are solved explicitly by $A=A_0\,\mathrm{law}(\bar\chi_c)$ and the lateral rate, and the balance rows by one sweep from the pinned node.

The integrator starts from a state whose algebraic rows hold. `consistent_state!` (through `_consistent_motion!`, `src/Integration/setup.jl`) makes one residual pass for the strain and history rates at the initial state, sets $A$ and $\theta_A$ from the column sums (`consistent_lateral`), makes a second pass at the consistent lateral rate, and sweeps the balance rows once from the pinned node (`consistent_velocities!`); a velocity row left above $10^{-10}$ of its scale is an error at setup. A deactivation event repeats these steps after its transfer (§14.3.5). On every proposal its error test accepts, the integrator's domain test (§15.5.4) rejects a non-positive area and requires the residuals of the area row and the lateral-rate row to lie within $r_{\mathrm{tol}}|A|$ and $r_{\mathrm{tol}}(|\theta_A|+\bar{|\theta|})$ plus their absolute tolerances, with $\bar{|\theta|}$ the volume-weighted mean of $|\theta_i|$.

---

## 10.7 Concentrations on a moving mesh

### 10.7.1 The problem

On a moving mesh the volume of a cell follows its strain: the nodes move so that a cell whose solid is being consumed shrinks. The concentration $\xi_j$ (mass per unit *bulk* volume) must respond to that volume change even when no mass crosses the cell's faces: when a cell shrinks, every component it carries is concentrated. Every component rides the mesh — the condensed components because they *are* the skeleton, the gas because its transport flux is defined relative to the skeleton (§3.2.3, §11.1.2), and every component crosses the faces at the relative velocity $-c$ when the mesh moves relative to the solid (§11.6) — so each component's cell mass changes only by reaction and by the fluxes through the faces. A formulation that integrated the concentrations would need a *dilation term* in every concentration row to produce the response to the volume change. Pyrolysis.jl integrates the masses instead, and the dilation becomes an identity.

### 10.7.2 Masses as states: the dilation is implicit

The integrated unknowns are the cell masses $m_{j,i}$, whose rows carry the reaction source and the face fluxes and nothing else:

$$
\frac{d m_{j,i}}{dt} \;=\; S_{j,i}\,V_i \;-\; \sum_{f\in\partial i} s_{i,f}\,A_f\,F_{j,f},
\tag{10.14}
$$

with $F_{j,f}\equiv 0$ for a condensed component on a mesh that moves with the solid. The state mapping recovers the concentration from the current volume, $\xi_{j,i}=m_{j,i}/V_i$. Differentiating,

$$
\frac{d\xi_{j,i}}{dt} = \frac{1}{V_i}\frac{dm_{j,i}}{dt} - \frac{m_{j,i}}{V_i^2}\frac{dV_i}{dt}
= S_{j,i} - \frac{1}{V_i}\sum_{f\in\partial i}s_{i,f}A_fF_{j,f} - \xi_{j,i}\,\theta_i ,
$$

using $dV_i/dt=\theta_iV_i$ (10.3), so that

$$
\boxed{\;\left.\frac{d\xi_j}{dt}\right|_{\mathrm{dilation}} = -\,\xi_j\,\theta_i,
\qquad \text{every component } j\;}
\tag{10.15}
$$

holds as an identity of the discrete system, not as a term the residual adds. It says: when a cell shrinks ($\theta<0$) the concentration of every component it carries rises (mass conserved in less volume); when a cell swells ($\theta>0$) the concentration falls. For the pore gas this is the compression of the gas in the shrinking pores, which the pressure closure $P = N/\varphi$ (§9.4) turns into a Darcy response.

> **Implementation note.** The row assembly writes the mass rows as `du[mass_index(layout, i, j)] = c.S_mass[j, i] + c.M_flux[j, i]` (`assemble_rows!`, `src/Assembly/residual.jl`): the reaction source plus the scattered face fluxes, with no geometric term. The state mapping (`derive!`, `src/States/derived.jl`) forms `ξ[j, i] = m / V_eval[i]` from the current geometry, with $\tilde V_i$ the volume smoothly floored at $10^{-8}V_i(0)$, and every kernel reads the concentrations from there. The derivation is repeated, in the ALE setting, in §11.3.4.

### 10.7.3 Why masses, not concentrations

Two properties follow from carrying the masses.

- **Exact mass conservation.** Every face flux enters its owner and its neighbour with opposite signs, and each reaction conserves mass, so the column mass changes only through the boundary fluxes, to round-off, whatever the mesh does (Chapter 16). There is no dilation term whose discretization could drift from the volume change.
- **The dilation follows whatever changes the volume.** It lives in $V_i$, which the geometry computes from the node positions and the area: the reactions' strain through the balance rows, the lateral law through $A$, a correction, or a prescribed mesh velocity all concentrate and dilute the inventories consistently. A uniform concentration field stays uniform when a prescribed velocity moves the mesh over a solid at rest (free-stream preservation; exact on the planar mesh in `test/kernels/test_relative_flux.jl`, and the verification case V-MMS-02).

---

## 10.8 Geometry coupling and observable consequences

### 10.8.1 Cell volume from state

On a moving mesh every cell volume is recomputed from the node positions and the current cross-section on each residual evaluation:

$$
V_i \;=\; (z_{i+1}-z_i)\,A,
\qquad A \equiv A(t)=A_0\,\mathrm{law}(\bar\chi_c).
\tag{10.16}
$$

> **Implementation note.** `geometry!` (`src/Meshes/geometry.jl`), called by the state mapping, sets the cell volumes, the face areas $A_f=A$ and the cell centres $(z_i+z_{i+1})/2$ from the node positions and the area unknown. All quantities carry the element type of the state, so automatic differentiation propagates the geometry sensitivities (e.g. $\partial V_i/\partial A$ and $\partial A_f/\partial A$) into every kernel that reads volumes or face areas. This is the link by which the area couples to the energy and species balances: a larger $A$ enlarges every $V_i$ and every face area, rescaling the face fluxes, the cell-integrated sources and every concentration $m/V$. The volume is signed: an inverted trial mesh gives a non-positive volume, which the state mapping floors smoothly for $m/V$ and the domain test reads unfloored, rejecting crossed nodes and a non-positive area (§11.8, §15.5.4).

### 10.8.2 Mass-loss-rate gauge for variable cross-section

A measurable consequence of $A(t)$ is the mass-loss rate (MLR) reported to the user. The outward gas flux of the exposed face is a *per-current-area* flux $\dot m''_{\mathrm{local}}$ (kg·m⁻²·s⁻¹ through the present surface $A(t)$). The total mass-loss rate is $\dot m_{\mathrm{total}}=\dot m''_{\mathrm{local}}\,A(t)$. For cone-calorimeter comparison, the convention is to report mass loss per unit *original* (gauge) area $A_0$:

$$
\dot m''_{\mathrm{gauge}} \;=\; \frac{\dot m_{\mathrm{total}}}{A_0}
\;=\; \dot m''_{\mathrm{local}}\,\frac{A(t)}{A_0}
\;=\; \dot m''_{\mathrm{local}}\,\mathrm{law}(\bar\chi_c).
\tag{10.17}
$$

> **Implementation note.** The solution assembly (`src/Outputs/solution.jl`) sums, at every saved time, the outward mass flux through every live boundary face times its current area, $\dot m_{\mathrm{total}}=\sum_f A_f\sum_jF_{j,f}$ [kg/s] — the top face and, when a permeable back-face condition vents, the bottom face too — and stores `sol.mass_loss_rate` $=\dot m_{\mathrm{total}}$ over the current area of the `:top` faces and `sol.extras.MLR_gauge` $=\dot m_{\mathrm{total}}$ over their initial area, exactly (10.17). The top-face flux alone is `sol.extras.surface_mass_flux`, and the area is `sol.extras.cross_section_area_history`. For a constant cross section the two flavors coincide. In the solve of §10.5.2 the ratio `MLR_gauge ./ mass_loss_rate` equals $A/A_0$ to round-off (a relative $10^{-16}$) at $t=10, 20, 30$ s (both rates are zero at $t=0$, where the ratio is `NaN`). This gauge normalization is what makes a shrinking or intumescing sample's predicted MLR directly comparable to a cone measurement, which always reports per nominal specimen area. The conservation ledgers use the current face areas (Chapter 16).

---

## 10.9 State-vector layout for variable cross-section

A lateral law adds three things to the state of a moving mesh: the history pair in every cell block, and the area and the lateral rate in the global block. On a line the state is a sequence of interleaved node and cell blocks followed by the global block (§15.1):

$$
\mathbf{u} = \big[\,\underbrace{z_1,\,w_1}_{\text{node }1},\;
\underbrace{m_{1,1},\dots,m_{N_C,1},\,T_1,\,G_1,\,m_{\mathrm{dry},1}}_{\text{cell }1},\;
\underbrace{z_2,\,w_2}_{\text{node }2},\;\dots,\;
\underbrace{\cdots}_{\text{cell }n},\;
\underbrace{z_{n+1},\,w_{n+1}}_{\text{node }n+1},\;
\underbrace{A,\,\theta_A}_{\text{global}}\,\big],
\tag{10.18}
$$

with the enthalpy inventory $H_i$ in place of $T_i$ under `energy_state = :enthalpy` and a face intensity appended to every node block under `BEER_LAMBERT`. The total length is $N_{\mathrm{state}} = n(N_C+3) + 2\,n_{\mathrm{nodes}} + 2$ (plus $n_{\mathrm{nodes}}$ intensities under `BEER_LAMBERT`). For $N_C=2$ and $n=5$ (so $n_{\mathrm{nodes}}=6$) this is $5\cdot5+6\cdot2+2=39$, which `test/closures/test_variable_cross_section.jl` asserts as `5 * (1 + 2 + 2) + 6 * 2 + 2`; the same problem without a lateral law has $5\cdot3+6\cdot2=27$, asserted as `5 * (1 + 2) + 6 * 2`. The solve of §10.5.2 ($N_C=2$, $n=10$) has 74 unknowns. A fixed mesh carries none of these rows: `setup` refuses a lateral law without `use_ale = true`.

The history rows are differential (ones on the mass matrix), the area and lateral-rate rows algebraic (zeros). Their absolute tolerances are the cell's mass tolerance $r_{\mathrm{tol}}\,\varepsilon_m\,m^s_i(0)$ for both members of the pair, $r_{\mathrm{tol}}A_0$ for the area, and for the lateral rate the smallest node-velocity tolerance divided by the initial thickness (`abstol_vector`, `src/States/tolerances.jl`). The rows are reached through `history_index(layout, i, k)`, `area_index(layout)` and `lateral_rate_index(layout)`, and their presence through `has_history(layout)` and `has_chi(layout)` (`src/States/layout.jl`). The colored Jacobian writes the two dense cross-section rows after its colored pass, through the six column sums (`LateralRows`, §15.9).

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

Char is typically far less dense than virgin ($\rho_{\mathrm{char}}\ll\rho_{\mathrm{virgin}}$) but is produced in smaller mass ($\nu_{\mathrm{char}}<1$). The sign of $\theta$ is set by the competition $\nu_{\mathrm{char}}/\rho_{\mathrm{char}}$ vs. $1/\rho_{\mathrm{virgin}}$: for a high-char-yield, very-low-density char the matrix can *swell* ($\theta>0$, the char occupies more pure-phase volume than the virgin it replaced); for a low-char-yield, dense char it *shrinks* ($\theta<0$). The gas, with $\gamma=0$, never contributes — it escapes through the pores. With $\theta<0$ the char mass $m_{\mathrm{char}}$, a state, is unchanged by the contraction while $V$ falls, so the char concentration $\xi_{\mathrm{char}}=m_{\mathrm{char}}/V$ rises at the rate $-\xi_{\mathrm{char}}\theta>0$ of (10.15), conserving char mass. Through the balance rows (10.13), $\theta<0$ in every cell gives a downward solid velocity ($w<0$ at every node above the bottom), so the surface node recedes — material burns away from the top while the substrate stays pinned.

**An exact endpoint.** Because the strain rate is linear in the reaction rate, the volume change of a cell under a single reaction $A\to\sum_j y_jB_j$ is linear in the converted mass and path-independent: with $r\,V\,dt=-dm_A$,

$$
dV = \Big(\frac{\gamma_A}{\rho_A}-\sum_j\frac{y_j\gamma_j}{\rho_j}\Big)\,dm_A ,
\qquad
\frac{V_f}{V_0} = \sum_j y_j\gamma_j\frac{\rho_A}{\rho_j}
\quad(\gamma_A=1,\ \xi_{A}(0)=\rho_A),
$$

whatever the temperature history. The verification case V-CPL-03 asserts this endpoint for seven permutations (a swelling product of half the reactant density, $V_f/V_0=2$; a shrinking one, $1/2$; a neutral one, $1$; char plus gas, $0.3$; mass loss at constant volume, $1$; a partial swelling factor, $1/2$; and the lateral law of §10.5.2, $A_f/A_0=0.75$), each to $3\times10^{-6}$, together with total-mass conservation.

---

## 10.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

### 10.11.1 ThermaKin2Ds — deformation and remeshing

TK Eq. 15 gives, in one dimension,

$$
\partial_t\xi_j=\sum_r\theta_r^j r_r-\partial_xJ_j+
\partial_x\left[\xi_j\int_0^x\rho^{-1}\partial_t\rho\,dx'\right].
$$

The full divergence contains both transport and dilation; it cannot be
replaced by $-\xi_j\rho^{-1}\partial_t\rho$. The reference plane is
stationary. TK §3 describes changing element thicknesses and conservative
splitting/combining, not a permanently fixed mesh. Its density law includes
solid and liquid volumes and a mixture-dependent gas swelling factor
(Eqs. 3–4). Our reaction-only $\sum_j\gamma_j\omega_j/\rho_j$ closure
does not include all density-change effects and is not a general equivalent.
TK restricts deformation to its primary coordinate; our empirical scalar
lateral-area law is a separate constitutive choice.

### 10.11.2 Gpyro — reaction yields and changing control volumes

GP Eqs. 22–24 define a **reaction-specific** solid mass yield

$$\mathrm{SF}_k=1+\left(\rho_{B_k}/\rho_{A_k}-1\right)\chi_k.$$

Together with bulk component densities and the mass balances, this determines
the changing $\Delta z$ (GP Eqs. 3, 56–57). It is not a per-species
shrink/non-shrink/swell flag. $\chi_k$ is not the conversion variable
$\alpha_A$ and neither is our dry-gas progress $\chi$.

For constant densities, fixed area, no additional volume contributions and
condensed $\gamma_j=1$, our volume change for $A\to y_BB+\text{gas}$
is proportional to $y_B/\rho_B-1/\rho_A$. Matching $y_B=\mathrm{SF}_k$
can match this restricted Gpyro volume balance. There is no universal
$\chi_k\leftrightarrow\gamma_j$ conversion. Our lateral law is not the
analogue of Gpyro's axial shrinkage law.

### 10.11.3 FDS — implemented density-sum rule and constraints

In FDS 6.11.0 `Source/wall.f90`, `REGRID_FACTOR` is accumulated as the sum
of updated component densities divided by their adjusted reference densities.
A remaining component with `ALLOW_SHRINKING=.FALSE.` bounds that factor below
by one; `ALLOW_SWELLING=.FALSE.` bounds it above by one. Component densities
are divided by the resulting factor to conserve mass, coordinates are
updated, and layers can be remeshed.

The same release's solid-phase manual still prints a max/sum expression.
That expression does **not** describe this source path: in particular, the
source does not select the most-expanding component as the swelling factor.
This comparison follows the [implemented regrid calculation](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90#L2388-L2455).
The V-CPL-03 endpoint examples are restricted checks of matched volume
balances, not proof that every FDS material maps to our $\gamma_j$.

### 10.11.4 Summary of mappings

| Model | Volume-change input | Numerical geometry |
|---|---|---|
| Pyrolysis.jl | reaction-volume sum with component $\gamma_j$; optional lateral law | coupled node positions/velocities and scalar area |
| ThermaKin2Ds | mixture density with gas swelling factor, TK Eqs. 3–4 | changing elements, splitting and combining |
| Gpyro 1D | reaction yield, component bulk densities and mass conservation | changing control-volume thicknesses |
| FDS 6.11.0 | density sum constrained by material shrink/swell flags | thickness and coordinate update, remeshing |

Primary source definitions are listed in §1.6.5. In particular, equality of
endpoint thicknesses does not establish equality of transient mesh dynamics.

---

## 10.12 Limitations and open issues

The following limitations are documented in the source and should guide use and future work:

1. **Density in $\theta$.** Equation (10.2) uses $\rho_j(T)$ evaluated at the *current* cell (trial) temperature, and the density-rate term $\xi_j\,d(1/\rho_j)/dt$ is neglected (§10.1.1). For materials with strongly temperature-dependent pure-phase densities and steep thermal gradients this is an approximation; whether retaining the term would materially change $\theta$ has not been quantified.

2. **Empirical lateral law, which must be smooth.** $\mathrm{law}(\bar\chi)$ is user-specified with no first-principles closure (§10.5.3). A percolation- or homogenization-based predictor of the functional form would remove a calibration degree of freedom, especially for intumescent coatings whose expansion factor can exceed an order of magnitude. The law must also be *smooth*: the lateral-rate row differentiates $\log\mathrm{law}$ by `ForwardDiff`, and the colored Jacobian differentiates that row again, so an exact Jacobian needs a twice-differentiable law. A piecewise-linear or table law gives a discontinuous $\theta_A$ at its kinks, degrades the Jacobian there, and makes any sensitivity taken across a kink ill-defined (§17.7.2).

3. **Vanishing cross-section.** The lateral rate (10.10) carries a $1/\mathrm{law}(\bar\chi)$ that would diverge as the cross-section collapses. `setup` refuses a law that is not positive and finite on the clamp band, and the domain test rejects a non-positive area, but there is no *physical* lower bound on shrinkage: a law that comes close to zero on $[0,1]$ makes $\theta_A$, and the stiffness of the cross-section rows, large. Users should design $\mathrm{law}$ so it stays bounded away from zero.

4. **$\bar\chi$ weighting choice.** The column average (10.7) is mass-weighted by *initial* dry-solid mass. Strongly non-uniform decomposition (e.g. a thin, fully-charred surface over an intact interior) is reduced to a single scalar $\bar\chi$, so the lateral law cannot represent depth-dependent cross-section change. A different weighting (or a depth-resolved area) would be required for that regime.

5. **Single-reactant kinetics.** Volume change inherits the single-reactant restriction of the kinetics (§6.5, §6.10): a reaction with more than one reactant is refused by problem validation and by source evaluation (`_apply_reaction_source!`, `src/Physics/kinetics.jl`), so $\theta$ and the history pair are built from single-reactant sources only. Bimolecular or multi-reactant stoichiometry (e.g. char oxidation consuming both char and O₂) would require extending the source-term machinery before $\theta$ and $G_i$ could account for it.

6. **Interaction with deactivation.** A deactivated cell's history pair goes to the receiver (§10.3.4, §14.3.3), which keeps $\bar\chi$ and $A$ continuous across the event; the conversion is a column quantity, so how the pair of a burned-out surface cell is attributed among the live cells does not matter to the area. The last live cell, lost with no receiver, takes its pair with it.

7. **Substring-based water-vapor exclusion.** The dry-gas index set relies on a name substring match (§10.3.3), which is fragile to naming. A type-level or explicit-flag mechanism for "this gas is moisture-derived" would be more robust than the `"H2O"`/`"water"` heuristic.

8. **One lateral strain for the whole column.** The axial (node-motion) and lateral (cross-section) channels are rows of one DAE, solved together, but the lateral channel is a single column-global $\theta_A$ subtracted from every cell's $\theta_i$. A cross-section that varies with depth, or lateral motion that differs between layers, would need a multi-dimensional model; the variable cross section is a line's (planar coordinates, one-dimensional mesh).

9. **The mesh follows the solid unless corrected.** The volume-change pipeline of this chapter defines the *solid* velocity $v_s$ through the balance rows; the mesh velocity is $w=v_s+c$, and it differs from the solid velocity only under the equidistribution correction (installed for burnout) or a prescribed mesh velocity for verification, with every inventory and the history pair then crossing the faces at the relative velocity $-c$ (§11.6, §11.7). A prescribed mesh velocity moves the mesh over a solid at rest, so its material must not react (§11.7.8).
