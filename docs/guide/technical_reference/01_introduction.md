# 1. Introduction, Scope, and Modeling Philosophy

This chapter establishes what Pyrolysis.jl is, the class of physical problems it solves, and the modeling decisions that govern every later chapter. Pyrolysis.jl is a one-dimensional (through-thickness) condensed-phase pyrolysis solver: it integrates the coupled heat-conduction, multi-species reaction-transport, and (optionally) gas-pressure equations that describe a solid heated to thermal decomposition. We position the code relative to the three established reference solvers it is most often compared to — ThermaKin2Ds, Gpyro, and the FDS solid-phase pyrolysis model — and we state the binding conventions (coordinate orientation, flux sign, divergence sign, heat-of-reaction sign) and the central "Formulation A" energy approximation that the remainder of the Technical Reference Guide derives and exploits. Everything asserted here is grounded in the actual implementation; where a textbook form and the code differ, the code is authoritative and the deviation is flagged.

---

## 1.1 Purpose and intended audience

Pyrolysis.jl simulates the thermal decomposition (pyrolysis) of condensed-phase materials subjected to an external heat flux — the canonical fire-science and materials-degradation problem of a solid that heats up, undergoes one or more Arrhenius-rate reactions, releases gaseous volatiles, and (in the charring/intumescent cases) changes its conductivity, porosity, and physical dimensions as it does so. The numerical core is a cell-centered finite-volume discretization on a possibly moving (Arbitrary Lagrangian–Eulerian, ALE) mesh, advanced in time by stiff implicit ODE integrators, with a structured (operator-assembled, automatic-differentiation-verified) Jacobian.

This document — the **Technical Reference Guide** — derives and justifies all of the mathematics and physics: the continuum governing equations (§3, Governing Equations of the Continuum Model), the constitutive closures (§4, Materials & Property Functions; §5, Effective Properties & Mixing Rules; §6, Reaction Kinetics; §8, Thermal Radiation), the discretization (§7, Heat Conduction, Convection & Energy Assembly; §13, Finite-Volume Discretization), the moving-mesh framework (§10, Volume Change & Lateral Shrinkage; §11, The ALE Framework), the boundary treatment (§12, Boundary & Initial Conditions), the temporal integration and Jacobian machinery (§15, Temporal Integration & Structured Jacobian), and the verification, conservation, and sensitivity apparatus (§16–§18). It is a *reference*, not a tutorial: it states the equations as the code implements them, with file-and-function attribution, and explains why each modeling choice was made. The companion **User Guide** (Part II of the documentation) teaches workflow — how to build materials, set boundary conditions, choose simulation modes, and run the solver.

The reader is assumed to be comfortable with the heat equation, the Arrhenius law, finite-volume / finite-difference discretization, and the basic vocabulary of stiff ODE integration. No prior knowledge of ThermaKin2Ds, Gpyro, or FDS is required, but readers familiar with any of them will find explicit notation-mapping and contrast subsections throughout.

### 1.1.1 How to read the two guides

The Technical Reference and the User Guide share a single, authoritative symbol table — the **unified nomenclature** reproduced verbatim in §2 (Nomenclature). Every chapter conforms to it. Three conventions in that table are *load-bearing* and are restated below in §1.5 because misreading any of them inverts a sign in a governing equation:

1. the spatial-coordinate orientation and the meaning of "top" and "bottom",
2. the positive-flux direction and the outflow-positive divergence operator, and
3. the heat-of-reaction sign convention.

A typical reading path for a new user is: §1 (this chapter) → §3 (the PDE system) → the User Guide quick-start, returning to §4–§9 as constitutive questions arise, and to §13–§15 when numerical behavior needs explaining. A developer extending the code should additionally read §11 (ALE) and §15 (Jacobian) in full.

---

## 1.2 Scope: the class of problems modeled

Pyrolysis.jl targets the **one-dimensional, through-thickness** decomposition of a planar sample. The spatial domain is a column of material of thickness $L$ and cross-sectional area $A$; the only resolved spatial coordinate is the axial coordinate $z$ (see §1.5.1). The unresolved lateral plane carries the area $A$, which may itself vary in time to represent lateral (radial) shrinkage or swelling (§10).

The state of the material is described by:

- the **temperature field** $T(z,t)$ [K], one scalar per finite-volume cell, and
- the **mass-concentration field** $\xi_j(z,t)$ [kg·m⁻³] for each chemical component $j \in \{1,\dots,N_C\}$.

The choice of mass concentration $\xi_j$ (rather than mass fraction $Y_j = \xi_j/\rho$) as the *primary* species state variable is deliberate: it makes the species conservation law a plain transport equation with a reaction source, with no embedded density derivative, and it is the form the reaction kinetics naturally produce (§6). Auxiliary states are introduced only when the corresponding physics is enabled — node positions $z_i$ for a moving mesh (ALE; §11) and the cumulative pyrolysis-progress variable $\chi_i$ for variable cross-section (§10).

### 1.2.1 Intended applications

The solver is designed for, and validated against, the following experiment and material classes (the shipped `examples/` directory exercises each):

- **Cone calorimetry.** A sample is irradiated at a controlled flux (commonly 25–75 kW·m⁻²) while surface temperature, mass-loss rate (MLR), and a heat-release-rate proxy are recorded. The exposed top face takes a radiative-plus-convective boundary condition; the back face is typically adiabatic and impermeable. Examples `06_cone_calorimeter.jl` (PMMA, 50 kW·m⁻²) and `07_douglas_fir_cone.jl` (charring wood with ALE and Darcy–Fick transport) are the reference configurations.
- **Idealized thermal-analysis experiments.** Thermogravimetric analysis (TGA), differential scanning calorimetry (DSC), and microscale combustion calorimetry (MCC) can be idealized as a thin, thermally-thin sample with a prescribed heating rate. In the limit of a single well-mixed cell, the species equations reduce to the 0-D Arrhenius mass-loss ODE that TGA inversion assumes; this idealization underpins the single-reaction verification case in §18 (Verification & Validation).
- **Charring, non-charring, and intumescent solids.** Non-charring polymers (PMMA) leave essentially no residue; charring solids (wood) form an insulating char layer with conductivity and porosity that differ sharply from the virgin material; intumescent coatings swell by orders of magnitude. The component/reaction/property-function machinery (§4) plus the swelling factor $\gamma_j$ and lateral-shrinkage law (§10) span all three.
- **Pressure-treated and wet wood / moisture transport.** Hygroscopic materials carry liquid water that evaporates endothermically and migrates as vapor under pressure gradients. The Darcy–Fick transport mode (§9, Gas-Phase Mass Transport & Pressure) and the `StateDependentProperty` heat-of-sorption closure (§4) handle these; examples `08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl`, and `moisture_hos_comparison.jl` are the reference cases.

### 1.2.2 What is **not** modeled

The scope boundary is as important as the scope itself. The following are explicitly out of scope in the current core:

- **Multi-dimensional geometry.** The model is strictly 1-D in $z$. The cross-section $A$ is a scalar gauge, not a resolved 2-D/3-D field; lateral shrinkage is a column-global area scaling driven by an empirical law (§10), not a resolved lateral transport. A genuine 2-D/3-D generalization of the ALE framework is listed as an open issue (§11).
- **Full two-temperature gas.** The gas phase is *not* given its own energy equation with a separate gas temperature. Gas is assumed to be in local thermal equilibrium with the matrix at the cell temperature $T$, and its sensible-heat *storage* is neglected (the quasi-steady-gas approximation, §1.4.1). Gas sensible-enthalpy *transport* is retained, carried by the advective source $S_{\text{conv}}$ (§7).
- **Multi-reactant reactions.** Each `Reaction` supports exactly one reactant ($N_{\text{in}}=1$ enforced at construction; see `Reaction` / `src/Materials/reactions.jl`). The reactant is consumed with unit mass stoichiometry and products are explicit yields summing to 1.0 (§6). Bimolecular or competitive oxygen-sensitive schemes (which ThermaKin and Gpyro support natively) must be approximated as sequences of single-reactant steps. Extending the reactant tuple to a 3-ary `(idx, order, mass_stoich)` form is an enumerated open question.
- **Pressure- or oxygen-dependent kinetics, and distributed activation energy.** The Arrhenius pre-exponential $A_i$ is temperature-independent, the rate has no explicit pressure or oxygen factor, and there is no built-in distributed-activation-energy (DAEM) tooling (§6).

These exclusions are not accidental: they are the modeling simplifications that make the 1-D solver fast, type-stable, and differentiable end-to-end, and they are shared (to varying degrees) by the reference codes discussed in §1.6.

---

## 1.3 Design philosophy: a generalized model with optional physics

The organizing principle of Pyrolysis.jl is a **single generalized model with opt-in physics layers**. Rather than maintaining separate solvers for "simple conduction," "charring wood," and "pressure-driven char," the code expresses one general residual and selects which physical terms are active through a compile-time trait, the `SimulationMode`. This keeps the simplest cases cheap (no unused machinery is evaluated) while allowing the full coupled problem to be assembled from the same kernel.

### 1.3.1 The `SimulationMode` trait and its four axes

The `SimulationMode{Mesh,Rad,Trans,Chi}` type (file `src/Problem/traits.jl`) bundles four orthogonal physics switches into one type parameter, so that the residual specializes at compile time without runtime branching and without method explosion:

| Axis | Trait type | Options | Physics toggled |
|---|---|---|---|
| Mesh motion | `MeshMode` | `Eulerian` / `ALE` | fixed mesh vs. moving (Lagrangian-solid) mesh; node positions $z_i$ become state in ALE |
| Radiation | `RadiationMode` | `NoRadiation` / `SurfaceAbsorption` / `BeerLambert` (P1 has no trait subtype; it is quarantined in `Experimental`) | how incident radiation deposits energy in depth |
| Gas transport | `TransportMode` | `Fickian` / `DarcyFick` | diffusion only vs. pressure-driven Darcy flow plus diffusion |
| Cross-section | `ChiMode` | `NoChi` / `WithChi` | constant area $A=A_0$ vs. variable $A(t)=A_0\,\mathrm{law}(\bar\chi)$ |

The compile-time accessors (`mesh_mode`, `radiation_mode`, `transport_mode`, `chi_mode`) extract the type parameters with zero runtime cost:

```julia
@inline mesh_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = M()
```

The unified entry point `residual!(du, u, def, t, ws)` dispatches first on `mesh_mode(def.mode)` into the Eulerian or ALE path, then on `transport_mode` into the Fickian or Darcy–Fick path (file `src/Residual/residual.jl`). The four axes combine into the small number of residual paths described in §15.1; modes that are not selected contribute *no* terms and allocate *no* work arrays (ALE work buffers are zero-length for Eulerian problems; the Darcy sub-cache exists only for `DarcyFick`; see `build_workspace`, `src/Problem/workspace.jl`).

> **Implementation note.** The trait is derived from user-level flags in `to_problem_def` (`src/Problem/adapter.jl`): `use_ale` selects `Eulerian`/`ALE`; `radiation_model` selects the radiation axis; the material's Darcy/χ traits (presence of permeability, presence of a `lateral_shrinkage_law`) select the transport and cross-section axes. The trait subtypes (`Eulerian`, `ALE`, …) are intentionally **not** exported at the top level to keep the public surface clean; they are reached via `Pyrolysis.Problem.Eulerian()` or explicit import.

### 1.3.2 Opt-in physics layers

The same opt-in philosophy governs the constitutive layers. Each of the following is off (or set to its cheapest closure) by default and is enabled only when the problem requires it:

- **Radiation in depth.** `NO_RADIATION` performs no volumetric absorption; `SURFACE_ABSORPTION` routes the absorbed flux through the surface energy balance (boundary Newton; no volumetric deposit) at $O(1)$ cost and preserves the tridiagonal Jacobian; `BEER_LAMBERT` resolves exponential in-depth absorption at the cost of a non-local (upper-triangular) Jacobian coupling; `P1_QUASI_STEADY` solves a diffusion equation for the radiation field. These are mutually exclusive and are queried by traits such as `requires_absorption_coefficients`, `has_nonlocal_jacobian_coupling`, and `requires_radiation_cache` (§8).
- **Darcy gas flow.** With `Fickian` transport, gas migrates only by (volume-fraction-gradient) diffusion. Selecting `DarcyFick` adds pressure-driven bulk flow $v_g=-(\kappa/\mu)\nabla P$ with the ideal-gas pressure closure (§9), and allocates the Darcy sub-cache.
- **ALE moving mesh.** A fixed Eulerian mesh suffices for thermally-thin or non-receding problems. Enabling `ALE` makes the mesh track the receding/expanding surface, with the node positions integrated as part of the coupled state and the Geometric Conservation Law enforced to machine precision (§11).
- **Variable cross-section.** `WithChi` enables lateral shrinkage/swelling through the pyrolysis-progress variable $\chi$ and an empirical area law (§10).
- **Adaptivity.** Surface-cell depletion merging is the production adaptivity path; r-refinement (node relocation) and the Experimental h-refinement are additional opt-in layers (§14, Adaptive Mesh Refinement).

The cost of an optional layer is therefore paid only when it is used, and the modeling assumptions of the default (cheapest) configuration are exactly those of the reference codes (§1.6).

### 1.3.3 Implementation goals stated up front

Three implementation properties are treated as first-class requirements throughout the code, and they drive many modeling micro-decisions documented in later chapters:

1. **Type stability.** Every hot-path function is written so the Julia compiler can infer concrete types, using heterogeneous tuples for the component/reaction lists and `@generated` functions to unroll loops over $N_C$ components and $N_R$ reactions (e.g. `_sum_heat_capacity`, `_gas_convective_source_unrolled`, `dry_pyrolysis_gas_indices`). The number of components $N_C$ and reactions $N_R$ are static type parameters of `Material{NC,NR,…}` (§4).
2. **Allocation-free assembly.** The residual evaluates with zero heap allocation. All scratch storage lives in pre-allocated caches (`RHSCache`, `SourceTermCache`, `PropertyCache`, `StateCache`) owned by the `Workspace`; the mutable mesh is read-only during assembly, and the typed `StateCache` is the single source of truth for state and state-derived geometry (§15.2). A dedicated test (`test/perf/test_zero_allocation_rhs.jl`) enforces this.
3. **AD-compatibility (differentiability).** The entire residual — including geometry updates, the boundary-temperature Newton solve, and property evaluations — propagates `ForwardDiff.Dual` numbers end-to-end. The `StateCache{NC,Tu,Mode}` carries a scalar eltype `Tu` that is `Float64` for ordinary evaluation and `ForwardDiff.Dual` for Jacobian/sensitivity passes. This requires every non-smooth construct to be replaced by a smooth surrogate: hard temperature gates become $\tfrac12(1+\tanh(\cdot))$ ramps, the depletion cutoff becomes a smooth multiplicative $\tanh(\xi/\text{threshold})$ rate factor, and the sensible-enthalpy integral is taken by the AD-clean midpoint rule (§4, §6). Differentiability is what makes the structured analytic Jacobian (§15) and forward/adjoint sensitivity analysis (§17, Sensitivity Analysis) possible.

These three goals are not merely engineering niceties: AD-compatibility in particular *constrains the physics formulation*, because every modeling term must be expressible as a smooth function of the state.

---

## 1.4 Formulation A: the canonical energy convention

The single most important modeling decision in Pyrolysis.jl is the form of the energy equation. The code solves the **non-conservative (Formulation A) temperature equation** with **gas sensible-heat storage excluded** from the volumetric heat capacity. This subsection states the equation, derives why the gas-storage term is dropped, and explains why this matches the reference codes.

### 1.4.1 The energy equation and the quasi-steady-gas approximation

The continuous energy equation solved is (verified against `compute_rhs!`, `src/Discretization/finite_volume.jl`, and the governing-equation docstring in `src/Properties/effective_properties.jl`):

$$
\rho c_{p}^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}}
\;+\; S_{\text{conv}}
\;+\; S_{\text{gen}}
\;+\; Q_{\text{rad}}
\;+\; Q_{\text{rxn}}.
\tag{1.1}
$$

The crucial feature is the left-hand side. The **effective (matrix-only) volumetric heat capacity** is

$$
\rho c_{p}^{\text{eff}}(T,\boldsymbol\xi)
= \sum_{j\in\{\text{solid,\,liquid}\}} \xi_j\,c_{p,j}(T),
\tag{1.2}
$$

i.e. the sum runs over **solid and liquid components only**; gaseous components are excluded. This is verified directly in the generated function `_sum_heat_capacity` (`effective_properties.jl`), which guards each term with `if !is_gas(comps[i])`. The total (matrix + gas) heat capacity $\rho c_p^{\text{total}}=\sum_{j}\xi_j c_{p,j}(T)$ is computed *only* for diagnostics (`total_heat_capacity`), never for the time-derivative term.

**Derivation of the exclusion.** Start from the full mixture energy balance, in which every phase stores sensible heat:

$$
\Big(\sum_{j}\xi_j c_{p,j}\Big)\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}} + (\text{gas enthalpy transport}) + Q_{\text{rad}} + Q_{\text{rxn}}.
$$

Split the storage sum into matrix and gas parts and move the gas storage to the right:

$$
\Big(\sum_{\text{S,L}}\xi_j c_{p,j}\Big)\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}} + (\cdots)
- \underbrace{\Big(\sum_{g}\xi_g c_{p,g}\Big)\frac{\partial T}{\partial t}}_{\text{gas sensible storage}}.
$$

The **quasi-steady-gas approximation** drops the underbraced term. Its justification is twofold:

1. *Timescale separation.* The gas residence time in the pore network (set by advection out of the heated layer) is much shorter than the matrix thermal-diffusion timescale. Over a matrix conduction time, the in-pore gas inventory turns over many times, so the gas is effectively in steady state at the local matrix temperature and its transient storage is negligible.
2. *Magnitude.* The in-pore gas mass concentration $\xi_g$ is tiny compared with the solid: at 1 atm the gas-storage term is on the order of **0.2 %** of the matrix term, with the ratio scaling linearly in pressure (the docstring in `effective_properties.jl` records this estimate). The approximation therefore costs $\sim$0.2 % at atmospheric pressure and degrades only at elevated pore pressure.

The dropped storage does not silently disappear: the gas *sensible-enthalpy transport* is retained through the advective source $S_{\text{conv}}$ (§1.4.3, §7), and the diagnostics machinery (§16, Conservation Diagnostics) reports the residual cost of the approximation as the **physical-form energy residual**, computed against $\rho c_p^{\text{total}}$, so the user can quantify the error for any given run.

### 1.4.2 Why "non-conservative"

Equation (1.1) is *non-conservative* because the temperature time derivative $\partial T/\partial t$, not a conserved quantity, appears on the left. The conserved quantity is the volumetric enthalpy $\sum_j \xi_j h_j(T)$ — **not** $\rho c_p^{\text{eff}}\,T$, which is not an enthalpy density ($h_j(T) \ne c_{p,j}\,T$ in general). A strictly conservative form would write $\partial_t\!\sum_j \xi_j h_j(T)$, whose product-rule expansion is $\rho c_p\,\partial_t T + \sum_j h_j(T)\,\partial_t\xi_j$: the companion term tracks the enthalpy carried by composition change. Substituting the species equations shows that most of this term is not lost — its reaction part becomes the reaction heat $Q_{\text{rxn}}$ (exactly, when the $h_r$ datum is the local reaction temperature) and its gas-transport part becomes the gas-enthalpy divergence carried by $S_{\text{conv}}$ (and, optionally, $S_{\text{gen}}$); the full derivation is in §3.3.2 (see also §6.6). Pyrolysis.jl uses the non-conservative form for three reasons: it is the form for which $T$ is the natural ODE unknown; it matches the reference codes (§1.6); and the conservation diagnostics (§16) audit the resulting bookkeeping so the discrete energy ledger still closes to integrator tolerance. What the temperature form gives up is exact *discrete* enthalpy conservation; that gap is one of the known limitations enumerated in §16.

### 1.4.3 The source terms and the $S_{\text{conv}}$ / transpiration mutual exclusivity

The right-hand-side sources in (1.1) are (full derivations in §7, §8, §6):

- $-\nabla\!\cdot\mathbf{q}_{\text{cond}}$, the conductive flux divergence, with Fourier's law $\mathbf{q}_{\text{cond}}=-k\,\nabla T$ and a distance-weighted harmonic-mean face conductivity (§1.4.4).
- $S_{\text{conv}} = -\sum_g c_{p,g}(T)\,\overline{J}_{g}\,\partial T/\partial z$, the **gas advective source**, carrying the sensible heat convected by the through-flow of volatiles. This term is **always active** in the standard energy form and is what compensates for the dropped gas storage.
- $S_{\text{gen}} = -\sum_g \Delta h_g(T,T_{\text{ref}})\,(\nabla\!\cdot\mathbf{J}_g)$, the **optional gas-generation enthalpy sink**, included only when `energy_form = :with_generation_sink`. It is a *partial* move toward conservative form and is consistent **only with $T_{\text{ref}}$-datum reaction heats** — with DSC-style heats quoted at the reaction temperature it double-counts the gas sensible enthalpy (§6.6.3, §7.5.4).
- $Q_{\text{rad}}$, the volumetric radiation absorption (model-dependent; §8).
- $Q_{\text{rxn}} = -\sum_r h_r r_r$, the reaction heat source (sign convention in §1.5.3).

A subtle but critical convention follows. The two source terms $S_{\text{conv}}$ and $S_{\text{gen}}$ together are the volumetric decomposition of the conservative-form sensible-enthalpy divergence:

$$
\nabla\!\cdot\!\Big(\sum_g \Delta h_g(T,T_{\text{ref}})\,\mathbf{J}_g\Big)
= \underbrace{\sum_g c_{p,g}\mathbf{J}_g\!\cdot\!\nabla T}_{-S_{\text{conv}}}
+ \underbrace{\sum_g \Delta h_g(T,T_{\text{ref}})\,\nabla\!\cdot\mathbf{J}_g}_{-S_{\text{gen}}},
\tag{1.3}
$$

an identity that is exact even for temperature-dependent $c_{p,g}$, since $\nabla\Delta h_g = c_{p,g}(T)\,\nabla T$. Integrated over the column, the two terms telescope to the surface transpiration outflow:

$$
\int (S_{\text{conv}} + S_{\text{gen}})\,dz = -\sum_g \dot m_g\,\Delta h_g(T_s,T_{\text{ref}}).
$$

Because $S_{\text{conv}}$ is *always* present, applying a **surface transpiration boundary term** $\sum_g \dot m_g \Delta h_g$ on top of it would double-count gas convective enthalpy by exactly $\int S_{\text{conv}}$, producing a spurious $\sim$5–10 kW·m⁻² sink that under-predicts surface temperature and MLR. For this reason the transpiration BC is **rejected at problem construction**: `to_problem_def` throws an `ArgumentError` if `use_transpiration_bc=true` (file `src/Problem/adapter.jl`). The supported convention is $S_{\text{conv}}$-only with a plain convective + radiative surface balance, which is exactly the ThermaKin Eq. 17 / FDS / Gpyro convention.
### 1.4.4 The discrete energy update

Discretizing (1.1) over a finite-volume cell $i$ gives the per-cell ODE that the integrator advances (verified against `compute_rhs!`, Step 8, `finite_volume.jl`):

$$
\frac{dT_i}{dt}
= \frac{1}{\rho c_{p,i}^{\text{eff}}}
\Big(
-(\nabla\!\cdot\mathbf{q})_i
+ S_{\text{conv},i}
+ S_{\text{gen},i}
+ Q_{\text{rad},i}
+ Q_{\text{rxn},i}
\Big),
\tag{1.4}
$$

with the conductive flux divergence assembled from face fluxes

$$
(\nabla\!\cdot\mathbf{q})_i = \frac{(q_{f_R}-q_{f_L})\,A}{V_i},
\qquad
q_f = -k_f\,\frac{T_R-T_L}{d_L+d_R},
\qquad
k_f = \frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}.
\tag{1.5}
$$

The face conductivity $k_f$ is the **distance-weighted harmonic mean** (resistance-in-series) of the two adjacent cell-center conductivities. This is derived from flux continuity across the interface (the total thermal resistance $R = d_L/(k_L A) + d_R/(k_R A)$ inverts to the harmonic mean) and reduces to the textbook $2k_Lk_R/(k_L+k_R)$ on a uniform mesh ($d_L=d_R$); on stretched meshes it improves near-boundary flux accuracy by roughly 1–5 % (§7). The same distance-weighted harmonic mean is applied to permeability and gas-transfer face values (§5, §9).

> **Implementation note.** Conductive face fluxes are computed in `conductive_flux` / `compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`); boundary faces are explicitly zeroed there and overwritten by the boundary-temperature solver. The matrix-only heat capacity is `effective_heat_capacity` (`src/Properties/effective_properties.jl`); the gas advective source is `compute_gas_convective_source!` (`src/Physics/convection.jl`); the divergence operator is `divergence_1d!` (`src/Discretization/finite_volume.jl`). The full eight-step RHS assembly sequence — interior conductive fluxes → reactions → gas transport → thermal/mass BCs → radiation → $S_{\text{conv}}$ → divergences/$S_{\text{gen}}$ → final RHS — is documented in §7.

The species equation, by contrast, is solved in conservative form (the conserved quantity $\xi_j$ is itself the unknown), and transport applies only to gas components (§3):

$$
\frac{d\xi_{j,i}}{dt}
= -\frac{(J_{j,f_R}-J_{j,f_L})\,A}{V_i}
+ S_{j,i},
\qquad
S_{j,i}=\sum_r \nu_{r,j}\,r_{r,i},
\tag{1.6}
$$

where the flux term is present for gaseous $j$ and zeroed for solids/liquids (`apply_species_residual_unrolled!`, §4).

---

## 1.5 Units and coordinate conventions

All quantities are SI (see §2 for the complete table). Three conventions are restated here because they fix signs in the governing equations.

### 1.5.1 Coordinate orientation

The single spatial coordinate is the **axial / through-thickness coordinate $z$** [m]:

- $z=0$ is the **bottom / substrate**, with boundary id **2** and tag `:bottom`;
- $z=L$ is the **top / exposed surface**, with boundary id **1** and tag `:top`.

Heat enters from the **top** (the irradiated face); a charring/receding material shrinks **downward** toward the substrate. The lateral plane carries the cross-sectional area $A$. The tag→id bridge (`:top→1`, `:bottom→2`) is resolved at compile time via `Val{tag}` lookup (§12). A user who needs to index raw state vectors should remember this orientation: temperature and concentration blocks are ordered by cell index $i=1,\dots,n$ from bottom to top.

### 1.5.2 Flux and divergence sign conventions

- **Positive flux = transport in the $+z$ direction** (bottom → top). Thus a positive conductive flux $q_f>0$ carries heat upward, and Fourier's law $q_f=-k_f(T_R-T_L)/(d_L+d_R)$ gives $q_f>0$ when $T_L>T_R$ (heat flows from hot to cold, opposite the gradient).
- **The discrete divergence is outflow-positive**: $(\nabla\!\cdot F)_i = (F_{f_R}-F_{f_L})A/V_i$ is positive when the quantity flows *out* of cell $i$. Consequently the energy RHS (1.4) carries $-(\nabla\!\cdot\mathbf{q})_i$, so a net conductive outflow ($>0$) cools the cell, as physically required.

### 1.5.3 Heat-of-reaction sign convention

The internal (storage) convention is **$h>0$ endothermic (cools)** and **$h<0$ exothermic (heats)**, with $h$ measured per kilogram of the (single) reactant [J·kg⁻¹]. The volumetric reaction heat source carries the minus sign:

$$
Q_{\text{rxn}} = -\sum_r h_r\,r_r,
\tag{1.7}
$$

verified verbatim in `reaction_heat_source` (`src/Physics/kinetics.jl`, `Q -= heats[i]*rates[i]`). Hence an endothermic reaction ($h>0$, $r>0$) gives $Q_{\text{rxn}}<0$ (cooling), and an exothermic reaction ($h<0$) gives $Q_{\text{rxn}}>0$ (heating).

> **Overload caution (binding).** ThermaKin and Gpyro *publish* the opposite sign — in their manuals $h>0$ denotes an exothermic reaction. When mapping literature heats of reaction into Pyrolysis.jl inputs, flip the sign. This is overload resolution H1 of the nomenclature (§2); the minus sign in (1.7) is what reconciles the two storage conventions.

Two further nomenclature overloads that recur in this guide and are worth flagging now (full list in §2): the symbol $A$ denotes **cross-sectional area** in the physics chapters, while the Arrhenius pre-exponential is always written subscripted, $A_i$ (or $A_{\alpha\beta}$); and $\chi$ always denotes **pyrolysis progress**, never a cross-section ratio (the area ratio is $A/A_0=\mathrm{law}(\bar\chi)$).

### 1.5.4 Fixed physical constants

The constants used throughout are defined in `src/Materials/constants.jl` and match the nomenclature table exactly:

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | $R_g$ | 8.31446261815324 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | $\sigma$ | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference temperature | $T_{\text{ref}}$ | 298.15 | K |
| Reference pressure | $P_{\text{ref}}$ | 101325 | Pa |
| Reference molar mass (air, Chapman–Enskog estimate) | $M_{\text{ref}}$ | 0.029 | kg·mol⁻¹ |

---

## 1.6 Relation to ThermaKin2Ds, Gpyro, and the FDS solid-phase model

Pyrolysis.jl deliberately adopts the established modeling conventions of three reference solvers so that its results are directly comparable and its inputs largely portable. This subsection maps notation and contrasts conventions; detailed equation-by-equation comparisons appear in the relevant later chapters (§3, §6, §7, §9, §15).

### 1.6.1 ThermaKin2Ds

ThermaKin2Ds is a 1-D/quasi-2-D condensed-phase pyrolysis solver widely used for inverse parameter estimation from cone/TGA data. The correspondence:

- **Energy convention.** ThermaKin's energy equation (its Eq. 17) carries gas sensible-enthalpy transport through a convective term $c_g(\mathbf{J}_g\!\cdot\!\nabla T)$ — this is exactly the Pyrolysis.jl $S_{\text{conv}}$ term, and adopting it (rather than a surface transpiration BC) is why transpiration is forbidden alongside $S_{\text{conv}}$ (§1.4.3). The one *deviation*: ThermaKin's storage term sums $c_p$ over **all** components including gas ($\sum_j \xi_j c_{p,j}$), whereas Pyrolysis.jl excludes gas storage (Formulation A, §1.4.1). Pyrolysis.jl thus matches ThermaKin's *transport* convention but takes the FDS/Gpyro *storage* approximation. The cost is the $\sim$0.2 %/atm physical-form residual.
- **Stoichiometry.** ThermaKin writes stoichiometric coefficients as $\theta_i^j$; Pyrolysis.jl uses $\nu_{i,j}$ (negative for reactants, positive for products). The symbol $\theta$ is reserved here for the volumetric strain rate, so $\theta_i^j$ appears only in ThermaKin-comparison passages and is mapped immediately to $\nu_{i,j}$ (overload A2, §2).
- **Kinetics.** ThermaKin supports bimolecular rates ($\xi_{\text{COMP1}}\xi_{\text{COMP2}}$); Pyrolysis.jl is single-reactant only (§6).
- **Boundary fluxes and mesh.** ThermaKin uses piecewise external/flame BCs with reflected incident flux, and a dynamically adjusted mesh; Pyrolysis.jl's ALE framework (§11) and boundary-temperature Newton solve (§12) are the analogs.

### 1.6.2 Gpyro

Gpyro is a generalized 1-D (with limited multi-D) condensed-phase pyrolysis model built on Patankar-style control-volume discretization.

- **Conservation form.** Gpyro writes control-volume conservation of $\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, gas mass, and energy. Its bulk density $\bar\rho$ and volume fractions $X_\alpha$ map to Pyrolysis.jl's mixture density $\rho$ and volume fractions $v_j$; its pre-exponential $Z$ maps to $A_i$; its permeability $K$ maps to $\kappa$; its conversion variable $\alpha_{\text{conv}}$ is distinct from the absorption coefficient $\alpha$ (overloads A3, A4, A5, §2).
- **Energy storage.** Like FDS, Gpyro carries the matrix heat capacity (gas storage excluded) — Pyrolysis.jl follows this, hence "FDS/Gpyro convention" in (1.2).
- **Gas momentum.** Gpyro uses a Darcian mass flux $\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ with a pressure-evolution equation; Pyrolysis.jl's `DarcyFick` mode uses $v_g=-(\kappa/\mu)\nabla P$ with an ideal-gas pressure closure (§9), omitting the buoyancy term (negligible in a heated thin sample).
- **Discretization and solve.** Gpyro uses Patankar's Peclet-based scheme selection, harmonic interface interpolation, and TDMA within a fully-implicit backward-Euler linearization. Pyrolysis.jl shares the harmonic interface interpolation (§1.4.4) but uses a general stiff ODE integrator with a structured analytic Jacobian instead of hand-linearized backward Euler (§13, §15).
- **Mixing rules.** Gpyro averages properties by volume fraction with power-law $T$-dependence; Pyrolysis.jl offers parallel, series, weighted, Bruggeman, and Carman–Kozeny rules per property (§5).

### 1.6.3 The FDS solid-phase pyrolysis model

FDS (Fire Dynamics Simulator) embeds a 1-D solid-phase pyrolysis model behind each wall/boundary cell of its CFD gas solver.

- **Energy and storage.** FDS solves a 1-D Fourier conduction equation with matrix-only heat capacity and Arrhenius source terms — the same Formulation A storage convention Pyrolysis.jl adopts (1.2). The FDS boundary condition $-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$ (convective + radiative) corresponds to Pyrolysis.jl's surface energy balance (§12).
- **Kinetics.** FDS rates can include density, oxygen, and temperature powers ($\rho^n O_2^n T^n$) and provide a recipe for recovering $A_i,E_i$ from a TGA reference temperature and peak rate; Pyrolysis.jl supports the same TGA-recovery recipe (§6, User Guide UG-4) but without the oxygen/density power factors.
- **Swelling/shrinkage.** FDS uses max/sum density-ratio rules for swelling and shrinkage; Pyrolysis.jl uses the swelling factor $\gamma_j$ and the $\chi$-driven lateral-shrinkage law (§10).
- **Radiation.** FDS uses a two-flux (Schuster–Schwarzschild) in-depth radiation treatment; Pyrolysis.jl offers surface absorption, Beer–Lambert, and an Experimental P1 diffusion model with diffusion coefficient $D_{\text{rad}}=1/(3(\alpha+\sigma_s))$ (§8).
- **Energy accounting.** The Pyrolysis.jl conservation diagnostics (§16) are the analog of FDS's energy-conservation accounting (its Appendix P), additionally quantifying the quasi-steady-gas approximation cost.

### 1.6.4 Summary contrast

| Aspect | ThermaKin2Ds | Gpyro | FDS solid phase | **Pyrolysis.jl** |
|---|---|---|---|---|
| Dimensionality | 1-D / quasi-2-D | 1-D (limited multi-D) | 1-D per face | **1-D** |
| Gas sensible storage in $\rho c_p$ | included | excluded | excluded | **excluded (Formulation A)** |
| Gas enthalpy transport | $c_g\mathbf J_g\!\cdot\!\nabla T$ | conservative div | surface BC | **$S_{\text{conv}}$ (= ThermaKin)** |
| Heat-of-reaction sign | $h>0$ exo | $h>0$ exo | — | **$h>0$ endo**, $Q=-\sum h r$ |
| Reactants per reaction | bimolecular | nth-order conv. | $\rho^nO_2^nT^n$ | **single reactant** |
| Gas momentum | Boyle-law Darcy | Darcy + buoyancy | — | **Darcy (no buoyancy)** |
| Time integration | Crank–Nicolson | implicit BE + TDMA | implicit | **stiff ODE + structured Jacobian** |
| Differentiable end-to-end | no | no | no | **yes (AD)** |

The last row is the distinguishing capability: Pyrolysis.jl is built to be automatically differentiable, enabling exact Jacobians and gradient-based sensitivity/inverse analysis (§15, §17) that the reference codes obtain only by finite differences.

---

## 1.7 Document roadmap

The Technical Reference Guide proceeds from continuum theory to discretization to numerics to verification:

- **§2 Nomenclature** — the verbatim symbol table and the three sign conventions; binding overload resolutions.
- **§3 Governing Equations of the Continuum Model** — the full coupled PDE system, the gas-vs-solid species dispatch, the Formulation A energy equation, the ideal-gas pressure closure, and side-by-side comparison tables with the reference codes.
- **§4 Materials & Property Functions** and **§5 Effective Properties & Mixing Rules** — the constitutive layer: components, phases, property-function types and their analytic derivatives, mixture density and porosity, and the conductivity/transfer/permeability mixing rules.
- **§6 Reaction Kinetics** — the Arrhenius rate law, smooth temperature gates, the depletion limiter, stoichiometry, and the heat-of-reaction source.
- **§7 Heat Conduction, Convection & Energy Assembly** — Fourier conduction, the gas advective and generation-sink sources, and the complete RHS assembly sequence.
- **§8 Thermal Radiation** — the four radiation models, their volumetric sources, and the P1 boundary conditions.
- **§9 Gas-Phase Mass Transport & Pressure** and **§10 Volume Change & Lateral Shrinkage** — diffusive and Darcy fluxes, the ideal-gas pressure, the strain-rate / $\chi$ machinery.
- **§11 ALE Framework**, **§12 Boundary & Initial Conditions**, **§13 Finite-Volume Discretization**, **§14 Adaptive Mesh Refinement** — the discretization and moving-mesh apparatus.
- **§15 Temporal Integration & Structured Jacobian**, **§16 Conservation Diagnostics**, **§17 Sensitivity Analysis**, **§18 Verification & Validation** — the numerical engine, the audits, and the V&V program.

Part II (the User Guide) mirrors this structure at the workflow level, beginning with installation and a quick-start cone-calorimeter run.

---

## 1.8 Limitations and open issues

The current model carries several known limitations, recorded here so that users size their expectations correctly (each is revisited in the relevant chapter):

1. **Dimensionality and lateral transport.** The model is strictly 1-D; lateral shrinkage is a column-global area scaling driven by an *empirical* law $A/A_0=\mathrm{law}(\bar\chi)$, not a resolved 2-D field. The physical justification of the shrinkage law and its interaction with multi-dimensional geometry is an open question (§10).
2. **Gas-storage approximation has no automatic high-pressure switch.** Formulation A drops gas sensible storage; the cost ($\sim$0.2 %/atm) is quantified by diagnostics but the code does not automatically restore gas storage at elevated pore pressure (§7, §16 open issues).
3. **Single-reactant kinetics only.** Multi-reactant (bimolecular, competitive) reactions are unsupported, as are temperature-dependent activation energy, pressure dependence, and distributed-activation-energy tooling (§6).
4. **Composition-change enthalpy not conserved discretely.** The temperature form does not evolve the composition-change enthalpy term $\sum_j h_j(T)\,\partial_t\xi_j$ of a conserved enthalpy density (§3.3.2); its reaction and gas-transport parts reappear as $Q_{\text{rxn}}$ and $S_{\text{conv}}$/$S_{\text{gen}}$, while the condensed-product sensible-generation part is omitted. The discrete energy ledger nonetheless closes to integrator tolerance, but the continuous closure for strongly composition-varying problems is approximate (§16).
5. **Asymmetric ambient + permeable back face.** With `:with_generation_sink`, a single scalar $T_{\text{ref}}$ is applied to every cell; if the two boundaries have different ambient temperatures *and* the back face is permeable, the energy ledger fails to close at the bottom by exactly $\sum_g\dot m_{g,\text{bottom}}c_{p,g}(T_{\infty,\text{top}}-T_{\infty,\text{bottom}})$ (harmless for the usual impermeable back face; §7, §12).
6. **Transport closures are first-order / heuristic in places.** Gas advection uses first-order upwind, diffusion applies one mixture transfer coefficient to every gas species (no per-species discrimination, with only heuristic Lennard-Jones parameter estimates available for sizing it), and viscosity mixing omits Wilke's rule (§5, §9).
7. **Experimental layers.** P1 radiation and h-AMR refinement/coarsening live in `Pyrolysis.Experimental`, use separate code paths (P1 has its own residual), and are explicitly unstable; the surface-depletion-callback adaptivity path is incompatible with the pluggable structured Jacobian and falls back to colored finite differences (§8, §14, §15).
8. **Imported calibrations must respect the S_conv-only convention.** Pyrolysis.jl never applies a surface transpiration term on top of $S_{\text{conv}}$ (§1.4.3); parameters calibrated against a solver that does will predict higher $T_s$ and MLR here and should be re-validated (§18).

These limitations bound the model's domain of applicability; within that domain — 1-D, single-reactant, atmospheric-to-moderate-pressure condensed-phase pyrolysis — the formulation is rigorous, differentiable, and conservation-audited.
