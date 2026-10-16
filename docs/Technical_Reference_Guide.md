# Pyrolysis.jl — Technical Reference Guide

A complete mathematical and physical description of the Pyrolysis.jl condensed-phase pyrolysis solver: governing equations with derivations and justifications, constitutive models, the finite-volume and ALE discretizations, the structured-Jacobian time integration, and the verification & validation basis. Companion to the **User Guide**.

_This document is the canonical source; the per-chapter views in `docs/guide/technical_reference` are regenerated from it with `julia scripts/split_guides.jl`. Each chapter is grounded in the Pyrolysis.jl source and was adversarially verified against the implementation._

---

## Contents

- [1. Introduction, Scope, and Modeling Philosophy](#1-introduction-scope-and-modeling-philosophy)
- [2. Nomenclature](#2-nomenclature)
- [3. Governing Equations of the Continuum Model](#3-governing-equations-of-the-continuum-model)
- [4. Materials, Components, and Property Functions](#4-materials-components-and-property-functions)
- [5. Effective (Mixture) Properties and Mixing Rules](#5-effective-mixture-properties-and-mixing-rules)
- [6. Reaction Kinetics and Stoichiometry](#6-reaction-kinetics-and-stoichiometry)
- [7. Heat Conduction, Gas Convection, and Energy Assembly](#7-heat-conduction-gas-convection-and-energy-assembly)
- [8. Thermal Radiation](#8-thermal-radiation)
- [9. Gas-Phase Mass Transport and Pressure](#9-gas-phase-mass-transport-and-pressure)
- [10. Volume Change, Swelling, Intumescence, and Lateral Shrinkage](#10-volume-change-swelling-intumescence-and-lateral-shrinkage)
- [11. The Arbitrary Lagrangian–Eulerian (ALE) Framework](#11-the-arbitrary-lagrangianeulerian-ale-framework)
- [12. Boundary and Initial Conditions](#12-boundary-and-initial-conditions)
- [13. Finite-Volume Spatial Discretization](#13-finite-volume-spatial-discretization)
- [14. Adaptive Mesh Refinement and Surface Depletion](#14-adaptive-mesh-refinement-and-surface-depletion)
- [15. Temporal Integration, DAE/Residual Formulation, and the Structured Jacobian](#15-temporal-integration-daeresidual-formulation-and-the-structured-jacobian)
- [16. Conservation Diagnostics](#16-conservation-diagnostics)
- [17. Sensitivity Analysis (Theory)](#17-sensitivity-analysis-theory)
- [18. Verification and Validation](#18-verification-and-validation)

---

# 1. Introduction, Scope, and Modeling Philosophy

This chapter establishes what Pyrolysis.jl is, the class of physical problems it solves, and the modeling decisions that govern every later chapter. Pyrolysis.jl is a one-dimensional (through-thickness) condensed-phase pyrolysis solver: it integrates the coupled heat-conduction, multi-species reaction-transport, and (optionally) gas-pressure equations that describe a solid heated to thermal decomposition. We position the code relative to the three established reference solvers it is most often compared to — ThermaKin2Ds, Gpyro, and the FDS solid-phase pyrolysis model — and we state the binding conventions (coordinate orientation, flux sign, divergence sign, heat-of-reaction sign) and the energy equation that the remainder of the Technical Reference Guide derives and exploits. Everything asserted here is grounded in the actual implementation; where a textbook form and the code differ, the code is authoritative and the deviation is flagged.

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
- **Full two-temperature gas.** The gas phase is *not* given its own energy equation with a separate gas temperature. Gas is assumed to be in local thermal equilibrium with the matrix at the cell temperature $T$: it stores sensible heat like every other phase, and its flux relative to the skeleton carries sensible enthalpy through the advective source $S_{\text{conv}}$ (§1.4.1, §7).
- **Multi-reactant reactions.** Each `Reaction` supports exactly one reactant ($N_{\text{in}}=1$; a reaction with more is rejected when its species source is applied, `_apply_reaction_source!` in `src/Physics/kinetics.jl`). The reactant is consumed with unit mass stoichiometry and products are explicit yields summing to 1.0 (§6). Bimolecular or competitive oxygen-sensitive schemes (which ThermaKin and Gpyro support natively) must be approximated as sequences of single-reactant steps. Extending the reactant tuple to a 3-ary `(idx, order, mass_stoich)` form is an open item (§1.8.1).
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

- **Radiation in depth.** `NO_RADIATION` performs no volumetric absorption; `SURFACE_ABSORPTION` routes the absorbed flux through the surface energy balance (boundary Newton; no volumetric deposit) at $O(1)$ cost and preserves the tridiagonal Jacobian; `BEER_LAMBERT` resolves exponential in-depth absorption at the cost of a non-local (upper-triangular) Jacobian coupling; `P1_QUASI_STEADY` solves a diffusion equation for the radiation field, lives in `Pyrolysis.Experimental`, and is rejected by `solve`. These are mutually exclusive (§8).
- **Darcy gas flow.** With `Fickian` transport, gas migrates only by (volume-fraction-gradient) diffusion. Selecting `DarcyFick` adds pressure-driven bulk flow $v_g=-(\kappa/\mu)\nabla P$ with the ideal-gas pressure closure (§9), and allocates the Darcy sub-cache.
- **ALE moving mesh.** A fixed Eulerian mesh suffices for thermally-thin or non-receding problems. Enabling `ALE` makes the mesh track the receding/expanding surface, with the node positions integrated as part of the coupled state and the Geometric Conservation Law enforced to machine precision (§11).
- **Variable cross-section.** `WithChi` enables lateral shrinkage/swelling through the pyrolysis-progress variable $\chi$ and an empirical area law (§10).
- **Adaptivity.** Depleted-cell parking is the only adaptivity `solve` uses. Error indicators, r-refinement (node relocation), and mesh-quality metrics are building blocks that `solve` does not call; h-refinement lives in `Pyrolysis.Experimental`, and its meshes cannot be passed to `solve` (§14, Adaptive Mesh Refinement).

The cost of an optional layer is therefore paid only when it is used, and the modeling assumptions of the default (cheapest) configuration are exactly those of the reference codes (§1.6).

### 1.3.3 Implementation goals stated up front

Three implementation properties are treated as first-class requirements throughout the code, and they drive many modeling micro-decisions documented in later chapters:

1. **Type stability.** Every hot-path function is written so the Julia compiler can infer concrete types, using heterogeneous tuples for the component/reaction lists and `@generated` functions to unroll loops over $N_C$ components and $N_R$ reactions (e.g. `_sum_heat_capacity`, `_gas_convective_source_unrolled`, `dry_pyrolysis_gas_indices`). The number of components $N_C$ and reactions $N_R$ are static type parameters of `Material{NC,NR,…}` (§4).
2. **Allocation-free assembly.** The residual evaluates with zero heap allocation. All scratch storage lives in pre-allocated caches (`RHSCache`, `PropertyCache`, `StateCache`) owned by the `Workspace`; the mutable mesh is read-only during assembly, and the typed `StateCache` is the single source of truth for state and state-derived geometry (§15.2). A dedicated test (`test/perf/test_zero_allocation_rhs.jl`) enforces this.
3. **AD-compatibility (differentiability).** The entire residual — including geometry updates, the boundary-temperature Newton solve, and property evaluations — propagates `ForwardDiff.Dual` numbers end-to-end. The `StateCache{NC,Tu,Mode}` carries a scalar eltype `Tu` that is `Float64` for ordinary evaluation and `ForwardDiff.Dual` for Jacobian/sensitivity passes. This requires every non-smooth construct to be replaced by a smooth surrogate: hard temperature gates become $\tfrac12(1+\tanh(\cdot))$ ramps, the depletion cutoff becomes a smooth multiplicative $\tanh(\xi/\text{threshold})$ rate factor, and every property law is evaluated on a temperature clamped to its validity bracket (§4.4, §6). Differentiability is what makes the structured analytic Jacobian (§15) and forward/adjoint sensitivity analysis (§17, Sensitivity Analysis) possible.

These three goals are not merely engineering niceties: AD-compatibility in particular *constrains the physics formulation*, because every modeling term must be expressible as a smooth function of the state.

---

## 1.4 The energy equation

The single most important modeling decision in Pyrolysis.jl is the form of the energy equation. The code solves the **non-conservative temperature equation** of a one-temperature porous medium: the pore gas is in local thermal equilibrium with the matrix, **every component stores sensible heat**, and the gas flux relative to the skeleton carries sensible enthalpy through a volumetric advective source. This subsection states the equation, derives it from the mixture enthalpy balance, and explains how it maps onto the reference codes.

### 1.4.1 The energy equation and its derivation

The continuous energy equation solved is (verified against `compute_rhs!`, `src/Discretization/finite_volume.jl`, and the docstring of `effective_heat_capacity` in `src/Properties/effective_properties.jl`):

$$
\rho c_{p}^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}}
\;+\; S_{\text{conv}}
\;+\; Q_{\text{rad}}
\;+\; Q_{\text{rxn}}.
\tag{1.1}
$$

The storage coefficient is the **mixture volumetric heat capacity**

$$
\rho c_{p}^{\text{eff}}(T,\boldsymbol\xi)
= \sum_{j=1}^{N_C} \xi_j\,c_{p,j}(T),
\tag{1.2}
$$

summed over **every** component, gas included. This is the generated function `_sum_heat_capacity` behind `effective_heat_capacity` (`effective_properties.jl`); the same value drives the time-derivative term, the structured Jacobian, the park-and-mask temperature merge, and the energy ledger.

**Derivation.** Let $\xi_j$ be the bulk mass concentration of component $j$, $v_s$ the velocity of the skeleton (every condensed component moves with it), and $\mathbf J_g = \xi_g(\mathbf v_g - \mathbf v_s)$ the gas mass flux *relative to the skeleton* — the quantity Darcy's law and the volume-fraction-gradient diffusion law return (§3.4, §3.5, §11.1.2). With $h_j = h_j^0 + \Delta h_j(T)$ and $\Delta h_j = \int c_{p,j}\,dT$, the enthalpy balance per bulk volume is

$$
\partial_t\Big(\sum_j \xi_j h_j\Big) + \nabla\!\cdot\Big(\mathbf v_s \sum_j \xi_j h_j\Big) + \nabla\!\cdot\Big(\sum_g h_g \mathbf J_g\Big)
= -\nabla\!\cdot\mathbf q_{\text{cond}} + Q_{\text{rad}} .
$$

Expanding each $\xi_j h_j$ with the product rule and substituting the species balance $\partial_t\xi_j + \nabla\!\cdot(\mathbf v_s\xi_j) = \dot\omega_j - \nabla\!\cdot\mathbf J_j$ ($\mathbf J_j = 0$ for condensed $j$) gives

$$
\Big(\sum_j \xi_j c_{p,j}\Big)\frac{D_s T}{Dt}
+ \sum_j h_j\dot\omega_j
- \sum_g h_g\nabla\!\cdot\mathbf J_g + \nabla\!\cdot\Big(\sum_g h_g\mathbf J_g\Big)
= -\nabla\!\cdot\mathbf q_{\text{cond}} + Q_{\text{rad}},
\qquad \frac{D_s}{Dt} \equiv \partial_t + \mathbf v_s\!\cdot\!\nabla .
$$

The two gas-flux pieces combine to $\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T = -S_{\text{conv}}$, and the reaction piece is $\sum_j h_j\dot\omega_j = \sum_r r_r\sum_j\nu_{r,j}h_j(T) = \sum_r h_r(T)\,r_r = -Q_{\text{rxn}}$ with $h_r(T)$ the heat of reaction at the local temperature (the calibrated, DSC-style datum; §6.6.3). Hence

$$
\Big(\sum_j \xi_j c_{p,j}\Big)\frac{D_s T}{Dt}
= -\nabla\!\cdot\mathbf q_{\text{cond}} + S_{\text{conv}} + Q_{\text{rad}} + Q_{\text{rxn}},
\tag{1.3}
$$

which is ThermaKin Eq. 17 term for term: the all-component storage sum, the convective term $c_g(\mathbf J_g\!\cdot\!\nabla T)$, and the deformation term, which is the $\mathbf v_s\!\cdot\!\nabla T$ part of $D_sT/Dt$. On a fixed mesh $D_s/Dt = \partial_t$; in ALE the mesh moves with the skeleton, so $D_sT/Dt$ is the mesh-frame rate and (1.1) holds with no advection term of any kind (§11.2.2).

**Why the gas is stored.** The gas concentration $\xi_g$ is an integrated state, so the term $\xi_g c_{p,g}\,\partial_t T$ is available at no cost and dropping it would be an approximation with no benefit. Gpyro drops it in its default thermal-equilibrium form (the quasi-steady-gas approximation, restored by its `FULL_QSG` option); FDS has no in-depth gas, so the term is identically zero there. On white-pine cone decks the term is at most about half a percent of $\rho c_p$ and shifts the response by a uniform time lag of about a tenth of a percent — small, but a thousand times above the solver's noise floor and therefore worth keeping exact.

### 1.4.2 Why "non-conservative"

Equation (1.1) is *non-conservative* because the temperature time derivative $\partial T/\partial t$, not a conserved quantity, appears on the left. The conserved quantity is the volumetric enthalpy $\sum_j \xi_j h_j(T)$ — **not** $\rho c_p^{\text{eff}}\,T$, which is not an enthalpy density ($h_j(T) \ne c_{p,j}\,T$ in general). A strictly conservative form would write $\partial_t\!\sum_j \xi_j h_j(T)$, whose product-rule expansion is $\rho c_p\,\partial_t T + \sum_j h_j(T)\,\partial_t\xi_j$: the companion term tracks the enthalpy carried by composition change. The derivation above shows where it goes — its reaction part becomes the reaction heat $Q_{\text{rxn}}$ (exactly, when the $h_r$ datum is the local reaction temperature) and its gas-transport part becomes $S_{\text{conv}}$ (see also §3.3.2 and §6.6). Pyrolysis.jl uses the non-conservative form for three reasons: it is the form for which $T$ is the natural ODE unknown; it matches the reference codes (§1.6); and the conservation diagnostics (§16) audit the resulting bookkeeping; the energy ledger closes to integrator tolerance only at constant composition and cell volume (§16.7). What the temperature form gives up is exact *discrete* enthalpy conservation; that gap is one of the known limitations enumerated in §16.

### 1.4.3 The source terms and the surface gas enthalpy

The right-hand-side sources in (1.1) are (full derivations in §7, §8, §6):

- $-\nabla\!\cdot\mathbf{q}_{\text{cond}}$, the conductive flux divergence, with Fourier's law $\mathbf{q}_{\text{cond}}=-k\,\nabla T$ and a distance-weighted harmonic-mean face conductivity (§1.4.4).
- $S_{\text{conv}} = -\sum_g c_{p,g}(T)\,\overline{J}_{g}\,\partial T/\partial z$, the **gas advective source**, carrying the sensible heat convected by the through-flow of volatiles relative to the skeleton. This term is always active; it is the only place the gas flux enters the energy equation.
- $Q_{\text{rad}}$, the volumetric radiation absorption (model-dependent; §8).
- $Q_{\text{rxn}} = -\sum_r h_r r_r$, the reaction heat source (sign convention in §1.5.3).

The surface itself carries no gas-enthalpy term. Under local thermal equilibrium the gas crosses the exposed face at $T_s$ on both sides, so its enthalpy flux is continuous there and the surface balance reduces to conducted flux equals applied flux, $q_{\text{BC}}(T_s) = q_{\text{cond}}(T_s)$. The sensible heat the escaping gas carries away is already accounted for in depth: its enthalpy at formation sits in $h_r(T)$, and its heating from the formation temperature to $T_s$ is $S_{\text{conv}}$ integrated along its path. A surface term $\sum_g \dot m_g \Delta h_g$ applied on top of $S_{\text{conv}}$ would count that heating a second time, producing a spurious $\sim$5–10 kW·m⁻² sink that under-predicts surface temperature and MLR. The surface balance is therefore a plain convective + radiative balance — the ThermaKin Eq. 17 convention, which Gpyro's thermal-equilibrium form shares. The only boundary use of the gas flux is the blowing model of a convective BC, which reads the outward gas flux to reduce the heat-transfer coefficient (§12.3.3).

### 1.4.4 The discrete energy update

Discretizing (1.1) over a finite-volume cell $i$ gives the per-cell ODE that the integrator advances (verified against `compute_rhs!`, Step 8, `finite_volume.jl`):

$$
\frac{dT_i}{dt}
= \frac{1}{\rho c_{p,i}^{\text{eff}}}
\Big(
-(\nabla\!\cdot\mathbf{q})_i
+ S_{\text{conv},i}
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

> **Implementation note.** Conductive face fluxes are computed in `conductive_flux` / `compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`); boundary faces are explicitly zeroed there and overwritten by the boundary-temperature solver. The mixture heat capacity is `effective_heat_capacity` (`src/Properties/effective_properties.jl`); the gas advective source is `compute_gas_convective_source!` (`src/Physics/convection.jl`); the divergence operator is `divergence_1d!` (`src/Discretization/finite_volume.jl`). The full eight-step RHS assembly sequence — interior conductive fluxes → reactions → gas transport with the mass BCs → thermal BCs → radiation → $S_{\text{conv}}$ → divergences → final RHS — is documented in §7.

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

The constants used throughout are defined in `src/Materials/constants.jl` (the air molar mass of the Chapman–Enskog estimate is `M_AIR` in `src/Properties/effective_properties.jl`) and match the nomenclature table exactly:

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

- **Energy convention.** ThermaKin's energy equation (its Eq. 17) is the same equation as (1.3): the storage term sums $c_p$ over **all** components including gas, the convective term $c_g(\mathbf{J}_g\!\cdot\!\nabla T)$ is exactly the Pyrolysis.jl $S_{\text{conv}}$ term, and its deformation term is the skeleton-following derivative that the ALE mesh frame absorbs (§11). Because the convective term carries the gas enthalpy, the surface balance has no gas-enthalpy term (§1.4.3). ThermaKin's Eq. 15 applies the same deformation term to every component, gas included — the frame choice Pyrolysis.jl makes in §11.1.2.
- **Stoichiometry.** ThermaKin writes stoichiometric coefficients as $\theta_i^j$; Pyrolysis.jl uses $\nu_{i,j}$ (negative for reactants, positive for products). The symbol $\theta$ is reserved here for the volumetric strain rate, so $\theta_i^j$ appears only in ThermaKin-comparison passages and is mapped immediately to $\nu_{i,j}$ (overload A2, §2).
- **Kinetics.** ThermaKin supports bimolecular rates ($\xi_{\text{COMP1}}\xi_{\text{COMP2}}$); Pyrolysis.jl is single-reactant only (§6).
- **Boundary fluxes and mesh.** ThermaKin uses piecewise external/flame BCs with reflected incident flux, and a dynamically adjusted mesh; Pyrolysis.jl's ALE framework (§11) and boundary-temperature Newton solve (§12) are the analogs.

### 1.6.2 Gpyro

Gpyro is a generalized 1-D (with limited multi-D) condensed-phase pyrolysis model built on Patankar-style control-volume discretization.

- **Conservation form.** Gpyro writes control-volume conservation of $\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, gas mass, and energy. Its bulk density $\bar\rho$ and volume fractions $X_\alpha$ map to Pyrolysis.jl's mixture density $\rho$ and volume fractions $v_j$; its pre-exponential $Z$ maps to $A_i$; its permeability $K$ maps to $\kappa$; its conversion variable $\alpha_{\text{conv}}$ is distinct from the absorption coefficient $\alpha$ (overloads A3, A4, A5, §2).
- **Energy storage and gas frame.** Gpyro's thermal-equilibrium form drops the gas storage $\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default (its Eq. 65) and restores it with `FULL_QSG`; Pyrolysis.jl always keeps it (1.2). Gpyro writes its gas balance (Eq. 58) on control volumes that deform with the solid, with only the Darcy flux crossing their faces — the same skeleton-bound frame as the ALE mesh (§11.1.2).
- **Gas momentum.** Gpyro uses a Darcian mass flux $\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ with a pressure-evolution equation; Pyrolysis.jl's `DarcyFick` mode uses $v_g=-(\kappa/\mu)\nabla P$ with an ideal-gas pressure closure (§9), omitting the buoyancy term (negligible in a heated thin sample).
- **Discretization and solve.** Gpyro uses Patankar's Peclet-based scheme selection, harmonic interface interpolation, and TDMA within a fully-implicit backward-Euler linearization. Pyrolysis.jl shares the harmonic interface interpolation (§1.4.4) but uses a general stiff ODE integrator with a structured analytic Jacobian instead of hand-linearized backward Euler (§13, §15).
- **Mixing rules.** Gpyro averages properties by volume fraction with power-law $T$-dependence; Pyrolysis.jl offers parallel, series, weighted, Bruggeman, and Carman–Kozeny rules per property (§5).

### 1.6.3 The FDS solid-phase pyrolysis model

FDS (Fire Dynamics Simulator) embeds a 1-D solid-phase pyrolysis model behind each wall/boundary cell of its CFD gas solver.

- **Energy and storage.** FDS solves a 1-D Fourier conduction equation with the condensed-phase heat capacity and Arrhenius source terms; it carries no in-depth gas, so the gas term of (1.2) is identically zero there. The FDS boundary condition $-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$ (convective + radiative) corresponds to Pyrolysis.jl's surface energy balance (§12).
- **Kinetics.** FDS rates can include density, oxygen, and temperature powers ($\rho^n O_2^n T^n$) and provide a recipe for recovering $A_i,E_i$ from a TGA reference temperature and peak rate; Pyrolysis.jl supports the same TGA-recovery recipe (§6, User Guide UG-4) but without the oxygen/density power factors.
- **Swelling/shrinkage.** FDS uses max/sum density-ratio rules for swelling and shrinkage; Pyrolysis.jl uses the swelling factor $\gamma_j$ and the $\chi$-driven lateral-shrinkage law (§10).
- **Radiation.** FDS uses a two-flux (Schuster–Schwarzschild) in-depth radiation treatment; Pyrolysis.jl offers surface absorption, Beer–Lambert, and an Experimental P1 diffusion model with diffusion coefficient $D_{\text{rad}}=1/(3(\alpha+\sigma_s))$ (§8).
- **Energy accounting.** The Pyrolysis.jl conservation diagnostics (§16) are the analog of FDS's energy-conservation accounting (its Appendix P).

### 1.6.4 Summary contrast

| Aspect | ThermaKin2Ds | Gpyro | FDS solid phase | **Pyrolysis.jl** |
|---|---|---|---|---|
| Dimensionality | 1-D / quasi-2-D | 1-D (limited multi-D) | 1-D per face | **1-D** |
| Gas sensible storage in $\rho c_p$ | included | dropped by default (`FULL_QSG` restores) | no in-depth gas | **included** |
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
- **§3 Governing Equations of the Continuum Model** — the full coupled PDE system, the gas-vs-solid species dispatch, the energy equation, the ideal-gas pressure closure, and side-by-side comparison tables with the reference codes.
- **§4 Materials & Property Functions** and **§5 Effective Properties & Mixing Rules** — the constitutive layer: components, phases, property-function types and their analytic derivatives, mixture density and porosity, and the conductivity/transfer/permeability mixing rules.
- **§6 Reaction Kinetics** — the Arrhenius rate law, smooth temperature gates, the depletion limiter, stoichiometry, and the heat-of-reaction source.
- **§7 Heat Conduction, Convection & Energy Assembly** — Fourier conduction, the gas advective source, and the complete RHS assembly sequence.
- **§8 Thermal Radiation** — the four radiation models, their volumetric sources, and the P1 boundary conditions.
- **§9 Gas-Phase Mass Transport & Pressure** and **§10 Volume Change & Lateral Shrinkage** — diffusive and Darcy fluxes, the ideal-gas pressure, the strain-rate / $\chi$ machinery.
- **§11 ALE Framework**, **§12 Boundary & Initial Conditions**, **§13 Finite-Volume Discretization**, **§14 Adaptive Mesh Refinement** — the discretization and moving-mesh apparatus.
- **§15 Temporal Integration & Structured Jacobian**, **§16 Conservation Diagnostics**, **§17 Sensitivity Analysis**, **§18 Verification & Validation** — the numerical engine, the audits, and the V&V program.

Part II (the User Guide) mirrors this structure at the workflow level, beginning with installation and a quick-start cone-calorimeter run.

---

## 1.8 Limitations and open issues

The current model carries several known limitations, recorded here so that users size their expectations correctly (each is revisited in the relevant chapter):

1. **Dimensionality and lateral transport.** The model is strictly 1-D; lateral shrinkage is a column-global area scaling driven by an *empirical* law $A/A_0=\mathrm{law}(\bar\chi)$, not a resolved 2-D field. The physical justification of the shrinkage law and its interaction with multi-dimensional geometry is an open question (§10).
2. **One temperature.** The matrix and the pore gas share one temperature (local thermal equilibrium); there is no two-temperature gas energy equation, so regimes with very fast, hot gas throughput are outside the model (§3.9).
3. **Single-reactant kinetics only.** Multi-reactant (bimolecular, competitive) reactions are unsupported, as are temperature-dependent activation energy, pressure dependence, and distributed-activation-energy tooling (§6).
4. **Composition-change enthalpy not conserved discretely.** The temperature form does not evolve the composition-change enthalpy term $\sum_j h_j(T)\,\partial_t\xi_j$ of a conserved enthalpy density (§3.3.2); its reaction and gas-transport parts reappear as $Q_{\text{rxn}}$ and $S_{\text{conv}}$, exactly when the heats of reaction are quoted at the reaction temperature (§3.3.2). The discrete energy ledger closes to integrator tolerance only for constant composition and cell volume, and the continuous closure for strongly composition-varying problems is approximate (§16.7).
5. **Transport closures are first-order / heuristic in places.** Gas advection uses first-order upwind, diffusion applies one mixture transfer coefficient to every gas species (no per-species discrimination, with only heuristic Lennard-Jones parameter estimates available for sizing it), and the pore-gas viscosity is one curve per material with no composition-dependent (Wilke) mixing (§5, §9).
6. **Experimental layers.** P1 radiation and h-refinement/coarsening live in `Pyrolysis.Experimental`: `solve` rejects `P1_QUASI_STEADY` (the P1 kernels are called directly), and an h-refined mesh cannot be passed to `solve`. Depletion handling (§14.3) parks cells without changing the state size, so it keeps the structured Jacobian; parameter sensitivities are valid only between park events (§17.7.3).
7. **Imported calibrations must respect the $S_{\text{conv}}$-only convention.** Pyrolysis.jl carries the gas enthalpy in $S_{\text{conv}}$ and applies no surface gas-enthalpy term (§1.4.3); parameters calibrated against a solver that applies both will predict higher $T_s$ and MLR here and should be re-validated (§18).

These limitations bound the model's domain of applicability; within that domain — 1-D, single-reactant, atmospheric-to-moderate-pressure condensed-phase pyrolysis — the formulation is rigorous, differentiable, and conservation-audited.

### 1.8.1 Known gaps and missing physics

The physics below is absent from the model. It is ordered by its expected effect on the current validation cases; items already listed above are cross-referenced rather than repeated.

**Likely to bias the current validation cases**

1. **Flame heat-flux feedback.** There is no ignition-triggered flame heat-flux boundary condition. After ignition a cone-calorimeter flame returns roughly 20–30 kW·m⁻² (PMMA) or 10–15 kW·m⁻² (wood) to the surface, 40–60 % on top of a 50 kW·m⁻² exposure, so parameters calibrated without it absorb that net-flux bias. The natural form is a flame-flux BC switched on by a critical-MLR or critical-$T_s$ callback. Flame feedback and the blowing reduction of the convective coefficient (the `blowing` model of `ConvectiveBC`, §12.3.3) act in opposite directions and partially cancel: enabling a blowing model alone on a cone example calibrated without either raises $T_s$ and MLR. Both belong in place before the cone examples are recalibrated.
2. **Multi-reactant kinetics** (item 3 above). A reaction with more than one reactant is rejected when its species source is applied (`_apply_reaction_source!`, `src/Physics/kinetics.jl`). This blocks char oxidation, recondensation, and autocatalysis (Prout–Tompkins fits cannot be represented). The extension point is a 3-ary reactant tuple `(idx, order, mass_stoich)`; it touches the reaction type, the kinetics, the volumetric strain rate, and the reaction-source Jacobian operator.

**Relevant to near-term applications**

3. **Oxidative pyrolysis and char oxidation with O₂ ingress** (about 32.8 MJ per kg of carbon). It dominates at sub-critical fluxes, after flameout, and in smoldering, and it controls the late-time MLR tail and the residue mass of the Douglas-fir case. It needs multi-reactant kinetics (item 2) and a per-species ambient concentration, which `ConvectiveMassBC` already accepts (`ξ_inf` per species). The upwind Darcy outflow suppresses ingress at high MLR.
4. **Equilibrium (pressure-dependent) moisture evaporation and recondensation.** Drying is a pressure-independent Arrhenius reaction; no saturation pressure $P_{\text{sat}}(T)$ appears anywhere. Vapor recondensation ahead of the drying front (which raises the local moisture content by 30–50 %) is absent, and at the 1.5–2 atm pore pressure of pressure-treated wood the boiling point shifts to 112–121 °C, which Arrhenius drying cannot represent. It needs item 2 and a smooth $r \propto (P_v - P_{\text{sat}}(T))$ closure; the pressure is already in the property cache.
5. **Backing layer and contact resistance at the back face.** There is no contact-resistance or substrate model. An adiabatic back face over-predicts late-time back-half temperatures (the missing loss is about 0.7–3 kW·m⁻² for 25 mm of wood over 10–20 min); a lumped Robin condition on the bottom face (radiative + convective, as in example `06_cone_calorimeter.jl`) is the available approximation. A contact conductance plus a semi-infinite-substrate $h(t)$ would close the gap.
6. **Kinetics expressiveness.** Normalized-conversion ($\alpha$) tracking, a $T^n$ pre-factor, and DAEM-style stacking are absent. Any real order $n \ge 0$ per reactant ($\xi^n$) is supported and equals $\xi_0^n(1-\alpha)^n$ up to a rescaling of $A_i$, so the gap is $\alpha$-normalization, not the conversion order.
7. **Transport-closure accuracy** (item 5 above): no composition-dependent (Wilke) viscosity mixing, air-default Lennard-Jones parameters for the Chapman–Enskog estimate, and first-order upwind Darcy advection with no TVD limiter; and in `DarcyFick` mode the diffusive flux keeps the isobaric driving potential $-(\lambda/T)\,\partial(\xi T)/\partial z$, dropping the term proportional to $\nabla(\varphi P)$ that a resolved pressure field and sharp char fronts carry (`gas_flux`, `src/Physics/gas_transport.jl`).

**Niche or future**

- Spectral (two-band) radiation and in-depth emission from glowing char (the regime for which §8 points to P1).
- Substrate reflection; scattering ($\sigma_s$ is plumbed only in the Experimental P1 model and never exercised); view-factor edge fall-off.
- Knudsen/Klinkenberg slip, $\kappa_{\text{eff}} = \kappa(1 + b/P)$: $b \approx 28$ kPa for 1 µm pores but about 280 kPa for 0.1 µm, up to about 4× $\kappa$ at 1 atm in sub-µm char.
- Stefan–Maxwell multicomponent diffusion (first relevant with O₂ ingress against the pyrolyzate efflux) and a tortuosity factor ($\phi/\tau$, a factor of 2–10) on the diffusivity and the gas-transfer coefficient.
- Liquid free-water capillarity (relevant only above the fiber-saturation point, about 30 % moisture content) and bound-water diffusion ($\sqrt{Dt} \approx 0.25$ mm over a cone test, negligible at fire heating rates).
- Mechanics (char cracking, delamination, spalling), anisotropic conductivity (about 2× along the wood grain), melting and dripping, and 2-D/3-D edge burning (about a 10–20 % effect on late-time MLR in cone tests; item 1 above).


---

# 2. Nomenclature

This is the **single, authoritative symbol table** for the Pyrolysis.jl Technical
Reference Guide and User Guide. Every writer **MUST** conform to these symbols,
meanings, and SI units verbatim. Where a symbol is used differently in different
source subsystems, the **conflict has been resolved here** and the resolution is
binding. Conflicts and overloads are called out explicitly in the
"Overloaded / reserved symbols" section at the end.

Conventions:
- The spatial coordinate is **z** (axial, through-thickness). `z = 0` is the
  **bottom / substrate** (boundary id 2, tag `:bottom`); `z = L` is the **top /
  exposed surface** (boundary id 1, tag `:top`). Heat enters from the top; the
  material shrinks downward. The transverse/lateral plane carries the
  cross-sectional area `A`.
- Sign convention for fluxes: **positive flux = transport in the +z direction**
  (bottom → top). Divergence `(F_R − F_L)·A/V` is positive when the quantity
  flows *out* of a cell.
- Heat-of-reaction sign convention: **h > 0 endothermic (cools)**,
  **h < 0 exothermic (heats)**; the volumetric heat source is `Q_rxn = −Σ h_r r_r`,
  so endothermic reactions give `Q_rxn < 0`. (This is the storage convention used
  internally; note ThermaKin/Gpyro publish `h > 0 = exothermic` — see overload
  note H1.)
- "Bulk" density `ρ_j` is the pure-phase density used by mixing rules and mixture
  density; "intrinsic/skeletal" density `ρ_{i,j}` (cell-wall density) is used
  ONLY in the porosity formula.
- Mass concentration `ξ_j` (kg/m³) is the **primary species state variable**, not
  mass fraction. `ξ_j = Y_j ρ`.

---

## 1. Latin letters

| Symbol | Meaning | SI units |
|---|---|---|
| `A` | Cross-sectional area of the 1D column (may vary in time via lateral-shrinkage law); also face area in the divergence operator | m² |
| `A_0` | Initial / calibration-reference cross-sectional area (fixed at setup) | m² |
| `A_{αβ}`, `A_i` | Arrhenius pre-exponential factor for reaction *i* (or *αβ*) | s⁻¹ (1st-order); m³·kg⁻¹·s⁻¹ (bimolecular) |
| `B` | Spalding blowing number `ṁ''c_{p,g}/h` of the convective blowing models (§12.3.3) | – |
| `c_p`, `c_{p,j}` | Specific heat capacity (of component *j*) | J·kg⁻¹·K⁻¹ |
| `c` | Speed of light (radiation energy density diagnostic only) | m·s⁻¹ |
| `D_{AB}`, `D_j` | Binary (Chapman–Enskog) diffusion coefficient | m²·s⁻¹ |
| `D_rad` | P1 radiation diffusion coefficient, `1/(3(α+σ_s))` | m |
| `d`, `d_L`, `d_R`, `d_{LR}` | Distance from cell center to face / between adjacent cell centers | m |
| `E`, `E_{αβ}`, `E_i` | Activation energy for reaction *i* | J·mol⁻¹ |
| `E_sensible` | Sensible energy of the domain, `Σ ρc_p^eff T V` over every component | J |
| `e` | Euler's number / base of natural log (literal) | – |
| `F` | Radiative view factor | – |
| `F_L`, `F_R` | Generic face flux (left/right) in the divergence operator | (flux units)·m⁻² |
| `f` | Generic field, or Newton-residual function `f(T_s)` | varies |
| `G` | Incident radiation intensity (P1 primary variable) | W·m⁻² |
| `G_ext` | Ambient radiation field at boundary, `4σT_∞⁴` | W·m⁻² |
| `g` | Gravitational acceleration (Darcy buoyancy, comparison only) | m·s⁻² |
| `H(·)` | Heaviside step (implemented as smooth tanh ramp) | – |
| `h` (kinetics) | Heat of reaction, per kg of first reactant | J·kg⁻¹ |
| `h_conv` | Convective heat-transfer coefficient | W·m⁻²·K⁻¹ |
| `h_m` | External mass-transfer (film) coefficient | m·s⁻¹ |
| `h_P` | Pressure-transfer coefficient (convective pressure BC) | m·Pa⁻¹·s⁻¹ |
| `Δh`, `Δh_g` | Sensible enthalpy difference `∫ c_p dT` | J·kg⁻¹ |
| `I(z)`, `I_in` | Radiation intensity along path / at cell entrance (Beer–Lambert) | W·m⁻² |
| `I_surface`, `I_ext` | Incident radiation intensity at domain surface / external source | W·m⁻² |
| `J_j` | Mass flux of (gas) component *j* at a face | kg·m⁻²·s⁻¹ |
| `J_j^diff`, `J_j^adv` | Diffusive / advective part of the gas flux | kg·m⁻²·s⁻¹ |
| `K`, `κ` | Permeability — **use `κ`**; `K` reserved for Gpyro-comparison text only | m² |
| `k`, `k_j` | Thermal conductivity (of component *j*) | W·m⁻¹·K⁻¹ |
| `k_eff`, `k_f` | Effective (mixture) / face-averaged thermal conductivity | W·m⁻¹·K⁻¹ |
| `k_B` | Boltzmann constant (reduced temperature `T* = k_B T/ε`) | J·K⁻¹ |
| `L` | Domain thickness (total material length); also characteristic length in convection correlations | m |
| `M`, `M_j` | Molar mass (of component *j*) | kg·mol⁻¹ |
| `M_ref` | Reference molar mass in the Chapman–Enskog combining rules (= 0.029, air; §5.10) | kg·mol⁻¹ |
| `M_dry` | Total initial dry-solid mass in the column | kg |
| `m_j` | Total mass of component *j* in domain | kg |
| `m_{dry,i}` | Initial dry-solid mass in cell *i* (χ denominator) | kg |
| `ṁ_g` | Outward gas mass flux at surface | kg·m⁻²·s⁻¹ |
| `N` | Partial-pressure sum (numerator of ideal-gas pressure) | Pa |
| `N_C` (`NC`) | Number of chemical components (type parameter) | – |
| `N_R` (`NR`) | Number of reactions (type parameter) | – |
| `n` | Number of finite-volume cells | – |
| `n_nodes`, `n_faces` | Number of mesh nodes / faces (= n+1 in 1D) | – |
| `n_{i,j}`, `n_s` | Reaction order of reaction *i* w.r.t. reactant *j* | – |
| `Nu` | Nusselt number (convection BC) | – |
| `P` | Gas pressure in the porous medium | Pa |
| `P_ref` | Reference pressure (= 101325) | Pa |
| `Q_rxn` | Volumetric reaction heat source, `−Σ h_r r_r` | W·m⁻³ |
| `Q_rad` | Volumetric radiation absorption source | W·m⁻³ |
| `Q_em`, `Q_abs` | P1 volumetric emission `4σαT⁴` / absorption `αG` | W·m⁻³ |
| `q`, `q_cond` | (Conductive) heat flux | W·m⁻² |
| `q_BC` | Boundary-condition heat flux (positive into domain) | W·m⁻² |
| `q'''`, `q_rad'''` | Volumetric radiation source (Beer–Lambert) | W·m⁻³ |
| `R_g` (`R`, `R_gas`) | Universal gas constant (= 8.314462618) | J·mol⁻¹·K⁻¹ |
| `r`, `r_i` | Reaction rate (per unit volume) of reaction *i* | kg·m⁻³·s⁻¹ |
| `r` (limiters) | Slope ratio argument of a flux limiter `ψ(r)` | – |
| `S` | Sutherland constant (= 110.4 K, air) in the viscosity law | K |
| `S_j`, `S_{j,i}` | Species source term for component *j* (in cell *i*) | kg·m⁻³·s⁻¹ |
| `S_conv` | Volumetric gas-advection (convection) energy source | W·m⁻³ |
| `S_rad` | Net volumetric P1 radiation source, `α(G − 4σT⁴)` | W·m⁻³ |
| `S_dry,i` | Dry-pyrolysis-gas source term in cell *i* | kg·m⁻³·s⁻¹ |
| `T` | Absolute temperature | K |
| `T_s`, `T_w` | Surface (boundary face) temperature | K |
| `T_cell` | Temperature at the cell center adjacent to a boundary | K |
| `T_∞` | Ambient / surroundings temperature | K |
| `T_ref` | Reference temperature for enthalpy datum (= 298.15) | K |
| `T_face`, `T_f` | Arithmetic-mean temperature at a face | K |
| `T_min,i`, `T_max,i` | Lower / upper temperature gate for reaction *i* | K |
| `T*` | Reduced temperature `k_B T / ε_AB` (collision integral) | – |
| `t` | Time | s |
| `u` | Flat ODE/DAE state vector | mixed |
| `du` | Time derivative of state vector | mixed |
| `V`, `V_i` | Cell volume (of cell *i*) | m³ |
| `v_g` | Bulk (Darcy) gas velocity at a face | m·s⁻¹ |
| `v_j` | Volume fraction of component *j* | – |
| `v_j^solid` | Solid volume fraction (excludes non-swelling gases) | – |
| `W` | Implicit-RK "W" matrix, `I/(γΔt) − J` | – |
| `w`, `w_i` | ALE mesh node velocity (at node *i*) | m·s⁻¹ |
| `w_cell` | Cell-centered ALE mesh velocity | m·s⁻¹ |
| `X_α` | Volume fraction (FDS/Gpyro comparison notation; internally `v_j`) | – |
| `x`, `y` | Transverse/lateral coordinates (multi-D comparison only) | m |
| `Y_j` | Mass fraction of component *j*, `ξ_j/ρ` | – |
| `Z` | Pre-exponential factor (Gpyro/FDS comparison notation; internally `A`) | varies |
| `z` | Axial / through-thickness coordinate (primary spatial variable) | m |
| `z_i` | Position of mesh node *i* (ALE state) | m |
| `z_i^c` | Cell-center position of cell *i* | m |
| `Δz`, `Δz_i` | Cell thickness (of cell *i*) | m |

---

## 2. Greek letters

| Symbol | Meaning | SI units |
|---|---|---|
| `α` (`α_eff`) | Volumetric absorption (extinction) coefficient | m⁻¹ |
| `α_j` | Mass-basis absorption coefficient of component *j* | m²·kg⁻¹ |
| `α_conv` | Conversion / extent variable (Gpyro kinetics comparison) | – |
| `β` | Weighting parameter of the WEIGHTED mixing rule, `k = βk_∥ + (1−β)k_series` | – |
| `γ` (`γ_j`) | Swelling factor of component *j* (1 = solids, 0 = non-swelling gas) | – |
| `γ_RK` | Implicit-RK stage coefficient (in `W = I/(γΔt) − J`) | – |
| `Δ` | Forward-difference / change operator (prefix) | varies |
| `δ_{ij}` | Kronecker delta | – |
| `δT`, `δξ_j` | Finite-difference perturbation step | K; kg·m⁻³ |
| `ε`, `ε_j` | Emissivity (of component *j*) | – |
| `ε_eff` | Effective (composition-weighted) emissivity | – |
| `ε_AB` | Lennard-Jones well depth (collision integral) | J |
| `ζ` | Dummy integration variable along the optical path | m |
| `η_i` | Refinement error-indicator value for cell *i* | varies |
| `θ`, `θ_i` | Volumetric strain rate (dilation) of cell *i*, `(1/V) dV/dt` | s⁻¹ |
| `θ_A` | Column-global **radial** volumetric strain rate, `(1/A) dA/dt` | s⁻¹ |
| `θ_L` | **Axial** volumetric strain rate (`θ − θ_A`) | s⁻¹ |
| `θ_i^j` | Stoichiometric coefficient (ThermaKin comparison notation; internally `ν_{i,j}`) | – |
| `κ` (`κ_j`, `κ_face`) | Permeability (Darcy) of component *j* / at a face | m² |
| `κ_s` | Solid absorption coefficient (FDS comparison notation; internally `α`) | m⁻¹ |
| `λ` (`λ_j`, `λ_face`, `λ_eff`) | Gas transfer (diffusion) coefficient | m²·s⁻¹ |
| `ℓ` | Characteristic particle/pore diameter (Carman–Kozeny) | m |
| `μ`, `μ_face` | Gas dynamic (Sutherland) viscosity | Pa·s |
| `μ_ref` | Reference viscosity (= 1.716×10⁻⁵ at 273.15 K, air) | Pa·s |
| `ν` | Kinematic viscosity (Gpyro Darcy comparison) | m²·s⁻¹ |
| `ν_{i,j}` | Stoichiometric coefficient (yield) of component *j* in reaction *i* (negative reactant, positive product) | – |
| `ξ_j`, `ξ_{j,i}` | Mass concentration of component *j* (in cell *i*) — **primary species state** | kg·m⁻³ |
| `tanh(ξ/ξ_th)` | Depletion rate factor (per reactant) | – |
| `ρ` | Mixture / bulk density (`Σ ξ_j` for total; mixture-density formula for effective) | kg·m⁻³ |
| `ρ_j` | Bulk (pure-phase) density of component *j* | kg·m⁻³ |
| `ρ_{i,j}` | Intrinsic / skeletal density of component *j* (porosity formula only) | kg·m⁻³ |
| `ρ_g` | Gas-phase density (ideal-gas law) | kg·m⁻³ |
| `ρc_p^eff` | Effective (mixture) volumetric heat capacity, every component | J·m⁻³·K⁻¹ |
| `σ` | Stefan–Boltzmann constant (= 5.670374419×10⁻⁸) | W·m⁻²·K⁻⁴ |
| `σ_s` | Scattering coefficient (≈ 0 for pyrolysis) | m⁻¹ |
| `σ_AB` | Lennard-Jones collision diameter | Å |
| `σ_mon` | Gaussian width of the interface monitor (r-refinement) | m |
| `τ`, `τ_i` | Optical thickness (total / cell-local), `α Δz` | – |
| `τ_ILU` | ILU drop tolerance | – |
| `φ` | Porosity (gas void fraction) | – |
| `χ`, `χ_i` | **Pyrolysis progress** in cell *i* (cumulative dry-gas mass released per unit initial dry-solid mass) | – |
| `χ̄` | Mass-weighted column-average pyrolysis progress | – |
| `ψ(r)` | Flux-limiter function | – |
| `Ω_D` | Diffusion collision integral (Neufeld) | – |
| `ω` (`ω(z)`) | Monitor function (r-refinement equidistribution) | varies |
| `ω̇'''` | Volumetric reaction rate (Gpyro comparison notation; internally `r`) | kg·m⁻³·s⁻¹ |

---

## 3. Subscripts

| Subscript | Meaning |
|---|---|
| `j` | Component / species index (1…N_C) |
| `i` | Cell index (1…n); also reaction index in kinetics context (see overload I1) |
| `r`, `αβ` | Reaction index |
| `g` | Gas phase / gaseous component |
| `s` | Solid / condensed phase |
| `ℓ`, `liq` | Liquid phase |
| `L`, `R` | Left / right of a face |
| `f` | Face-centered quantity |
| `face` | Face-averaged property |
| `eff` | Effective (mixture) property |
| `cell` | Cell-centered quantity |
| `s` (boundary) | Surface (boundary face) — e.g. `T_s`, `ξ_{j,s}` |
| `∞` | Ambient / surroundings |
| `ref` | Reference value |
| `in` | Cell-entrance (radiation) |
| `ext` | External source |
| `min`, `max` | Lower / upper bound |
| `dry` | Dry-pyrolysis-gas (excludes water vapor) |
| `init`, `0` | Initial value |
| `top`, `bot` | Top (z=L) / bottom (z=0) boundary |
| `rxn` | Reaction-source contribution |
| `cond`, `conv`, `rad` | Conduction / convection / radiation contribution |
| `∥`, `series` | Parallel / series mixing limit |
| `total` | All-component total (e.g. `ρ_total = Σ_j ξ_j`) |

## 4. Superscripts and accents

| Mark | Meaning |
|---|---|
| `°`, `^o`, `_init` | Value at previous time step / initial value |
| `+`, `−` (radiation) | Forward / backward hemispheric intensity (two-flux) |
| `solid` | Restricted to solid/liquid (non-swelling) components |
| `diff`, `adv` | Diffusive / advective flux component |
| `^c` | Cell-center (e.g. `z_i^c`) |
| `‾` (overbar) | Column / spatial average (e.g. `χ̄`, `J̄_g`) or Gpyro-style mixture average (`ρ̄`) |
| `·` (overdot) | Time rate / flux (e.g. `ṁ`, `ω̇`) |
| `'`, `''`, `'''` | Per-length / per-area / per-volume (FDS/Gpyro comparison notation) |

## 5. Operators

| Operator | Meaning |
|---|---|
| `∂/∂t` | Partial time derivative (Eulerian, fixed point) |
| `∂/∂t\|_χ` | ALE material-point (moving-mesh) time derivative |
| `d/dt` | Total time derivative |
| `∇`, `∇·` | Gradient / divergence (1D: `∂/∂z`, `(F_R−F_L)A/V`) |
| `∇²` | Laplacian / second derivative (curvature indicator) |
| `Σ`, `∫` | Sum / integral |
| `clamp(x,a,b)` | Clamp `x` to `[a,b]` |
| `tanh` | Hyperbolic tangent (smooth gates and depletion limiter) |
| `‖·‖_∞` | Infinity norm |
| `diag(v)` | Diagonal matrix from vector `v` |
| `J`, `A`, `C_j` | Jacobian / its sparse block / structured coupling *j* |

---

## 6. Physical constants (fixed values)

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | `R_g` | 8.314462618 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | `σ` | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference pressure | `P_ref` | 101325 | Pa |
| Reference temperature | `T_ref` | 298.15 | K |
| Reference molar mass (air) | `M_ref` | 0.029 | kg·mol⁻¹ |
| Reference viscosity (air, 273.15 K) | `μ_ref` | 1.716×10⁻⁵ | Pa·s |
| Sutherland constant (air) | `S` | 110.4 | K |
| Kozeny constant (spheres) | — | 180 | – |

---

## 7. Overloaded / reserved symbols — resolutions (BINDING)

These symbols are overloaded somewhere in the source subsystems. The resolution
below is authoritative; writers must use the disambiguating subscript/context.

- **χ (chi) — A1 [pyrolysis progress vs. cross-section ratio].** The
  `volume_change`, `problem_residual`, and `geometry` subsystems define `χ` as
  the **pyrolysis progress** (dimensionless, cumulative dry-gas mass released per
  unit initial dry-solid mass); it is easily misread as the cross-section area
  ratio. **Resolution:** `χ` is *always* pyrolysis progress;
  the cross-section ratio is `A/A_0 = law(χ̄)` and is never written as `χ`.

- **θ (theta) — A2 [strain rate vs. stoichiometry].** Internally `θ` is the
  **volumetric strain rate** (s⁻¹), with `θ_A` radial and `θ_L` axial. ThermaKin
  uses `θ_i^j` for stoichiometric coefficients. **Resolution:** use `ν_{i,j}` for
  stoichiometry everywhere; reserve `θ` for strain rate. `θ_i^j` appears only in
  the ThermaKin-comparison passage and is immediately mapped to `ν_{i,j}`.

- **A — A3 [area vs. pre-exponential factor].** `A` is the **cross-sectional
  area**; the Arrhenius pre-exponential is `A_i` / `A_{αβ}` (always subscripted by
  reaction). In Jacobian/linear-algebra context `A` also denotes the sparse
  Jacobian block. **Resolution:** bare `A` = area in physics chapters; `A_i` =
  pre-exponential (kinetics); `A` (linear algebra) only inside the Jacobian
  chapter where the meaning is unambiguous and stated. The FDS/Gpyro
  pre-exponential `Z` is mapped to `A_i` on first use.

- **κ vs. K — A4 [permeability].** Use **`κ`** for permeability throughout.
  `K` (Gpyro permeability) appears only in comparison text and is mapped to `κ`.

- **α — A5 [absorption vs. conversion vs. heat-transfer].** `α` (or `α_eff`) =
  **volumetric absorption coefficient** (m⁻¹); `α_j` = mass-basis absorption
  coefficient (m²·kg⁻¹). Conversion (Gpyro) is `α_conv`; FDS `κ_s` maps to `α`.
  The mass-transfer "absorptivity" of a surface is also `α` but always paired
  with a BC context and Kirchhoff `α = ε_eff`.

- **λ — A6 [gas transfer vs. anything else].** `λ` is reserved for the **gas
  transfer (diffusion) coefficient** (m²·s⁻¹). Do not use `λ` for thermal
  conductivity (that is `k`) or eigenvalues.

- **ρ subscripts — A7 [bulk vs. intrinsic].** `ρ_j` = bulk density (mixing,
  mixture density); `ρ_{i,j}` = intrinsic/skeletal density (porosity only). Never
  write a bare `ρ_i` for skeletal density without the explicit "intrinsic"
  qualifier; prefer `ρ_{i,j}`.

- **h — H1 [heat of reaction; sign].** `h` (kinetics) = heat of reaction per kg
  of first reactant, with **h > 0 endothermic** internally; ThermaKin/Gpyro
  publish `h > 0 = exothermic`. The conversion to the volumetric source carries
  the minus sign: `Q_rxn = −Σ h_r r_r`. The convective coefficient is `h_conv`,
  the mass coefficient `h_m`, the pressure coefficient `h_P`, all subscripted.

- **i — I1 [cell index vs. reaction index].** `i` is the **cell index** in
  discretization/geometry/ALE; in the kinetics chapter `i` is the **reaction
  index** (with `j` the component index). Each chapter states its convention in
  its first paragraph. Where both appear together, use `i` = cell, `r` =
  reaction.

- **φ — A8 [porosity only].** `φ` is **porosity**; do not reuse for view factor
  (that is `F`) or any potential. `φ` denotes the intrinsic-density form
  `1 − Σ ξ_j/ρ_{i,j}`; the all-phase gas volume fraction `Σ v_j^gas` is a
  different quantity and is not written `φ`.

- **S — S1 [source vs. Sutherland constant].** `S` with a subscript
  (`S_conv`, `S_rad`, `S_j`) is a source term; bare `S` is the
  Sutherland constant in the viscosity law only.


---

# 3. Governing Equations of the Continuum Model

This chapter states the full coupled system of partial differential equations that
Pyrolysis.jl solves *before* it is discretized. We introduce the primary state
variables (cell temperature $T$ and component mass concentrations $\xi_j$) and the
auxiliary states carried only in certain simulation modes (node positions $z_i$ for
the Arbitrary Lagrangian–Eulerian mesh, and the pyrolysis-progress scalar $\chi_i$
for variable cross-section), then derive each conservation law: condensed-phase
mass, per-species transport, the non-conservative energy equation,
the ideal-gas pressure closure, and Darcy momentum. The chapter is deliberately a
*statement and justification* of the continuum model; the discrete operators that
realize it are derived in §7 (Heat Conduction, Convection & Energy Assembly), §9
(Gas-Phase Mass Transport & Pressure), §11 (ALE Framework), and §13 (Finite-Volume
Discretization). Every equation here is verified against the implementing source and,
where the textbook form and the code differ, the chapter documents what the **code**
actually does.

All symbols, units, and sign conventions are those fixed in §2 (Nomenclature). The
three binding sign conventions are repeated here because the entire chapter rests on
them:

1. **Flux positivity.** A flux is positive when it transports its quantity in the
   $+z$ direction, i.e. from the bottom/substrate ($z=0$, boundary id 2,
   `:bottom`) toward the top/exposed surface ($z=L$, boundary id 1, `:top`).
2. **Divergence positivity.** The discrete divergence
   $(\nabla\cdot\boldsymbol F)_i = (F_R - F_L)\,A/V_i$ is positive when the
   quantity flows *out* of cell $i$ on net. In the energy and species balances the
   divergence enters with a leading minus sign, so net outflow lowers the stored
   quantity.
3. **Heat-of-reaction sign.** Internally $h_r > 0$ is **endothermic** (the reaction
   cools the matrix) and $h_r < 0$ is **exothermic**. The volumetric reaction heat
   source therefore carries an explicit minus sign,
   $Q_{\mathrm{rxn}} = -\sum_r h_r r_r$, so an endothermic reaction
   ($h_r>0$, $r_r>0$) gives $Q_{\mathrm{rxn}}<0$. (ThermaKin and Gpyro publish the
   opposite convention $h>0=$ exothermic; the mapping is given in §3.7 and in
   overload note H1 of §2.)

---

## 3.1 The continuum model and its state variables

### 3.1.1 Domain and idealization

Pyrolysis.jl models a one-dimensional column of reacting condensed-phase material
of thickness $L$ and (possibly time-varying) cross-sectional area $A(t)$. The column
exchanges heat and mass with its surroundings through two boundaries: the **top**
($z=L$, exposed to the heat source) and the **bottom** ($z=0$, the substrate). The
physical picture is a slab in a cone calorimeter, a TGA/DSC/MCC sample idealized as a
thin film, or a charring/intumescing solid. The transverse plane is not resolved; all
lateral effects (lateral shrinkage, swelling of the cross-section) are folded into the
scalar area $A(t)$ via the law of §3.6 and §10.

The condensed phase is a **multi-component, multi-reaction mixture**. Each component
$j \in \{1,\dots,N_C\}$ is one of three phases — solid, liquid, or gas — encoded at
the type level (`SolidComponent`, `LiquidComponent`, `GaseousComponent`) so that
phase-specific behaviour dispatches at compile time (§4). Solids and liquids form the
load-bearing *matrix*; gases are decomposition products that percolate through the
matrix pores and ultimately escape. Reactions convert a single reactant component into
one or more product components (§6).

### 3.1.2 Primary and auxiliary state variables

The continuum state at every point $z$ and time $t$ is

$$
\boxed{\; T(z,t) \quad\text{and}\quad \xi_j(z,t),\; j=1,\dots,N_C. \;}
$$

- $T$ [K] is the **single** (condensed-phase + in-pore gas) temperature. Pyrolysis.jl
  is a *one-temperature* model: the gas and the matrix are assumed in local thermal
  equilibrium at each point, justified by the fast heat exchange between the pore
  gas and the pore walls relative to the matrix thermal-diffusion time. Every
  component, gas included, stores sensible heat at that temperature (§3.3).
- $\xi_j$ [kg·m⁻³] is the **mass concentration** of component $j$ — mass of $j$ per
  unit *bulk* volume of the column. This is the canonical species state variable, *not*
  the mass fraction $Y_j = \xi_j/\rho$. Choosing $\xi_j$ avoids ever forming $\rho$ at
  intermediate points and degrades gracefully as a component vanishes.

In the discrete model the state is cell-centered: $T_i$ and $\xi_{j,i}$ for cells
$i=1,\dots,n$. Two **auxiliary state blocks** are appended only in the modes that need
them (the bundled `SimulationMode{Mesh,Rad,Trans,Chi}` trait selects which; see §15):

- **Node positions** $z_i$ [m], $i=1,\dots,n_{\mathrm{nodes}}$, present only in ALE
  mode (`Mesh = ALE`). The mesh becomes part of the coupled ODE system; its evolution is
  governed by a mesh-velocity field $w$ (§3.6, §11).
- **Pyrolysis progress** $\chi_i$ [–], present only in variable-cross-section mode
  (`Chi = WithChi`). Each $\chi_i$ accumulates dry-pyrolysis-gas mass released per unit
  initial dry-solid mass and drives the lateral-shrinkage law (§3.6, §10).

The flat state vector is packed **block-major**,

$$
\boldsymbol u =
\bigl[\, \underbrace{T_1,\dots,T_n}_{T\text{-block}},\;
\underbrace{\xi_{1,1},\dots,\xi_{1,n},\dots,\xi_{N_C,n}}_{\xi\text{-blocks}},\;
\underbrace{z_1,\dots,z_{n_{\mathrm{nodes}}}}_{z\text{-block (ALE)}},\;
\underbrace{\chi_1,\dots,\chi_n}_{\chi\text{-block (WithChi)}} \,\bigr],
$$

with total length
$N_{\mathrm{state}} = n(1+N_C) + \delta_{\mathrm{ALE}}\,n_{\mathrm{nodes}} + \delta_{\mathrm{WithChi}}\,n$,
where $\delta_X = 1$ if the mode includes block $X$ and $0$ otherwise. The layout
indicator functions (`temperature_range`, `concentration_range`, `z_range`,
`chi_range`) and the pack/unpack routines are detailed in §13 and §15.

> **Implementation notes.** The state layout and block ranges are defined by
> `StateLayout` in `src/Problem/state_layout.jl`; the typed mirror used during
> residual assembly is `StateCache` in `src/Problem/state_cache.jl`; the bundled
> simulation-mode trait is in `src/Problem/traits.jl`. The unified right-hand side
> is `residual!` in `src/Residual/residual.jl`, which dispatches on `mesh_mode`.

### 3.1.3 Overview of the coupled system

Collecting the laws derived below, the continuum model is the coupled system

$$
\begin{aligned}
&\text{Energy:}\quad
\rho c_p^{\mathrm{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\cdot\mathbf q_{\mathrm{cond}} + S_{\mathrm{conv}}
+ Q_{\mathrm{rad}} + Q_{\mathrm{rxn}}, \\[4pt]
&\text{Gas species:}\quad
\frac{\partial \xi_j}{\partial t}
= -\nabla\cdot\mathbf J_j + \left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}, \\[4pt]
&\text{Condensed (solid/liquid) species:}\quad
\frac{\partial \xi_j}{\partial t}
= \left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}},
\end{aligned}
$$

closed by the constitutive laws inventoried in §3.8 (Fourier conduction, Arrhenius
kinetics, mixture properties, the ideal-gas pressure $P = N(\xi,T)/\varphi$, Darcy
momentum $v_g = -(\kappa/\mu)\nabla P$, and the volume-change/ALE relations). The
energy equation is **non-conservative** — it is written for $T$, not for an enthalpy
density — with every component, gas included, in the storage coefficient; the
choice is justified in §3.3. Each term and its sign are derived in the sections
that follow.

---

## 3.2 Condensed-phase mass and total density

### 3.2.1 Mixture (total) density

The bulk **total density** of the mixture is the sum of all component mass
concentrations,

$$
\boxed{\;\rho = \sum_{j=1}^{N_C} \xi_j. \;}
$$

This is exact by the definition of $\xi_j$ as mass-per-bulk-volume and requires no
property evaluation. Mass fractions $Y_j = \xi_j/\rho$ are formed from it wherever
they are needed.

The **mixture density used by the swelling/specific-volume bookkeeping**,
`mixture_density`, is a *different* quantity: it is the density a parcel would have if
its components occupied their pure-phase volumes, every component contributing volume
through its swelling factor $\gamma_j$ (the universal γ weighting, shared with the
dilation rate θ and the conductivity/emissivity volume fractions, §5):

$$
\frac{1}{\rho_{\mathrm{mix}}}
= \sum_{j} \gamma_j\,\frac{Y_j}{\rho_j},
$$

where $\rho_j(T)$ is the **bulk pure-phase density** of component $j$ (overload note A7
of §2: $\rho_j$ is the bulk density, distinct from the skeletal/intrinsic density
$\rho_{i,j}$ that appears only in porosity). For canonical materials
($\gamma_j = 1$ condensed, $\gamma_j = 0$ gas) the gas terms drop out and the
mixture occupies only the solid/liquid volume; for intumescent gases
($\gamma_j=1$) the escaping gas creates new matrix volume; a $\gamma_j = 0$
condensed phase (the `LiquidComponent` default) is volumeless. The
distinction between $\rho=\sum_j\xi_j$ (total) and $\rho_{\mathrm{mix}}$ (specific-volume
form) matters because the two enter different parts of the model: $\rho=\sum\xi_j$ is
the conserved bulk inventory; $\rho_{\mathrm{mix}}$ and the per-component specific
volumes feed the volume-change rate of §3.6. The full derivation of the swelling
mixture density and its derivatives is in §5.

> **Implementation notes.** `mixture_density` is in `src/Properties/mixture_rules.jl`;
> the $\rho_{\mathrm{mix}}$ swelling form is `_specific_vol_with_swelling_direct` (same
> file).

### 3.2.2 No separate bulk-continuity PDE

Pyrolysis.jl does **not** solve an independent continuity equation for $\rho$. Total
mass is a *derived* diagnostic: summing the per-species equations of §3.2.3 over $j$
gives the evolution of $\rho$ automatically. Because every reaction is mass-balanced
at build time (the single reactant is consumed with stoichiometric mass $1$ and the
product yields sum to $1$, §6), the reaction contribution telescopes,

$$
\sum_{j=1}^{N_C} \left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}
= \sum_r \Bigl(\underbrace{-1}_{\text{reactant}} + \underbrace{\textstyle\sum_{\text{products}}\nu}_{=\,1}\Bigr) r_r
= 0,
$$

so reactions move mass *between* components without creating or destroying it, and the
only net change in $\rho$ comes from gas transport out of the cell. This build-time
invariant is enforced by the `Reaction` constructor (`_validate_stoichiometry`,
relative tolerance $10^{-10}$); `mass_balance_error` and `verify_stoichiometry` report
it (§6).

### 3.2.3 Phase-resolved mass: solids/liquids are immobile, only gases transport

The per-species conservation law is the foundation of the model. In integral
(finite-volume) form over a cell $i$,

$$
\frac{d}{dt}\!\int_{V_i}\!\xi_j\,dV
= -\oint_{\partial V_i}\!\mathbf J_j\cdot\mathbf n\,dA
+ \int_{V_i}\! \left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}\!dV ,
$$

where $\mathbf J_j$ is the **mass flux** of component $j$ *relative to the skeleton*.
The decisive modeling choice is *which components transport*. In Pyrolysis.jl the
matrix is rigid relative to itself: solids and liquids do **not** flux between cells
— they are consumed/produced in place by reactions and (in ALE) carried by the mesh,
never by molecular transport. Only gas components diffuse and advect through the
pore network, and they do so *relative to the skeleton*: Darcy's law and the
volume-fraction-gradient diffusion law both return the gas velocity relative to the
pore walls, so the laboratory-frame gas flux is $\xi_g\mathbf v_s + \mathbf J_g$ and the
pore gas is carried with the pores it occupies (§11.1.2). Hence the **gas-vs-solid
dispatch**:

$$
\boxed{\;
\frac{\partial \xi_j}{\partial t} =
\begin{cases}
-\,\nabla\cdot \mathbf J_j + \displaystyle\left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}},
& \text{if } j \text{ is gas},\\[10pt]
\displaystyle\left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}},
& \text{if } j \text{ is solid or liquid.}
\end{cases}
\;}
$$

The continuous gas flux $\mathbf J_j$ is the Fickian or Darcy–Fick flux closed in §3.4,
§3.5 and detailed in §9; the reaction source
$\left.\partial\xi_j/\partial t\right|_{\mathrm{rxn}} = \sum_r \nu_{r,j}\,r_r$ is the
stoichiometry-weighted sum of reaction rates ($\nu_{r,j}<0$ reactant, $>0$ product;
§6). The discrete realization writes, for each cell $i$,

$$
\frac{d\xi_{j,i}}{dt} =
\begin{cases}
-\dfrac{(J_{j,f_R} - J_{j,f_L})\,A}{V_i} + S_{j,i}, & j \text{ gas},\\[10pt]
S_{j,i}, & j \text{ solid/liquid},
\end{cases}
\qquad S_{j,i} = \sum_r \nu_{r,j}\,r_{r,i}.
$$

This is exactly the branch taken by the generated kernel
`apply_species_residual_unrolled!`: for a gas component it assigns
`dξ = -div_J + dξ_rxn`; for any non-gas component it assigns `dξ = dξ_rxn` only — the
divergence is never even formed for the matrix phases.

> **Implementation notes.** The compile-time gas/solid branch is
> `apply_species_residual_unrolled!` in `src/Materials/materials.jl`.
> Species sources $S_{j,i}=\sum_r\nu_{r,j}r_r$ come from `species_source_terms!` in
> `src/Physics/kinetics.jl`; the per-reaction application enforces the single-reactant
> mass-balance convention in `_apply_reaction_source!` (same file). The species
> divergence $(J_R-J_L)A/V_i$ is `divergence_species_1d!` in
> `src/Discretization/finite_volume.jl`.

---

## 3.3 The energy equation (non-conservative)

### 3.3.1 Statement

Pyrolysis.jl solves the **non-conservative** energy equation written directly for the
temperature $T$, with the time-storage term carried by the **mixture** volumetric
heat capacity:

$$
\boxed{\;
\rho c_p^{\mathrm{eff}}(T,\xi)\,\frac{\partial T}{\partial t}
= -\,\nabla\cdot \mathbf q_{\mathrm{cond}}
+ S_{\mathrm{conv}}
+ Q_{\mathrm{rad}}
+ Q_{\mathrm{rxn}} .
\;}
$$

The terms, with sign conventions, are:

| Term | Meaning | Continuous form | Sign / units |
|---|---|---|---|
| $\rho c_p^{\mathrm{eff}}\,\partial_t T$ | sensible-heat storage, every component | $\bigl(\sum_{j}\xi_j c_{p,j}(T)\bigr)\partial_t T$ | J·m⁻³·s⁻¹ |
| $-\nabla\cdot\mathbf q_{\mathrm{cond}}$ | Fourier conduction | $\mathbf q_{\mathrm{cond}} = -k_{\mathrm{eff}}\nabla T$ | net outflow cools |
| $S_{\mathrm{conv}}$ | gas advection of sensible enthalpy | $-\sum_g c_{p,g}(T)\,\mathbf J_g\!\cdot\!\nabla T$ | W·m⁻³ |
| $Q_{\mathrm{rad}}$ | volumetric radiation absorption | model-dependent | W·m⁻³ |
| $Q_{\mathrm{rxn}}$ | reaction heat | $-\sum_r h_r r_r$ | W·m⁻³ |

The discrete cell update is

$$
\frac{dT_i}{dt}
= \frac{1}{\rho c_{p,i}^{\mathrm{eff}}}
\Bigl(
-(\nabla\cdot\mathbf q)_i + S_{\mathrm{conv},i}
+ Q_{\mathrm{rad},i} + Q_{\mathrm{rxn},i}
\Bigr),
$$

which is precisely the final assembly loop of `compute_rhs!` (the divergence enters
with a minus sign, the three volumetric sources add, the whole is divided by the
mixture $\rho c_p$). The detailed discretization of each term — the distance-weighted
harmonic-mean face conductivity, the second-order cell-centered temperature gradient in
$S_{\mathrm{conv}}$, the eight-step assembly order — belongs to §7; here we justify the
*continuous* form.

> **Implementation notes.** The final RHS assembly
> $dT_i/dt = (-\,\mathrm{div\_q} + S_{\mathrm{conv}} + Q_{\mathrm{rad}} + Q_{\mathrm{rxn}})/\rho c_p$
> is in `compute_rhs!`, `src/Discretization/finite_volume.jl` (Step 8).
> The mixture heat capacity is `effective_heat_capacity` and the conductive flux is
> `conductive_flux` / `compute_conductive_fluxes!`.

### 3.3.2 The mixture heat capacity and where the gas enthalpy goes

The storage coefficient is the **mixture** (solid + liquid + gas) volumetric heat
capacity,

$$
\rho c_p^{\mathrm{eff}}(T,\xi)
= \sum_{j=1}^{N_C} \xi_j\, c_{p,j}(T).
$$

The pore gas is at the matrix temperature (one-temperature model, §3.1.2) and its
concentration $\xi_g$ is an integrated state, so its sensible storage
$\xi_g c_{p,g}\,\partial_t T$ is available at no cost and is kept. Dropping it — the
*quasi-steady-gas approximation* of Gpyro's default thermal-equilibrium form and
of FDS (which carries no in-depth gas at all) — would be an approximation with no
benefit here. On white-pine cone decks the term is at most about half a percent of
$\rho c_p$ and produces a uniform time lag of about a tenth of a percent.

Where the rest of the gas enthalpy goes follows from the derivation of §1.4.1,
summarized here. Write the enthalpy of component $j$ as $h_j = h_j^0 + \Delta h_j(T)$
with $\Delta h_j = \int c_{p,j}\,dT$, let every condensed component move at the
skeleton velocity $\mathbf v_s$, and let $\mathbf J_g = \xi_g(\mathbf v_g - \mathbf v_s)$
be the gas flux relative to the skeleton (§3.2.3). The enthalpy balance per bulk
volume,

$$
\partial_t\Big(\sum_j \xi_j h_j\Big) + \nabla\!\cdot\Big(\mathbf v_s \sum_j \xi_j h_j\Big)
+ \nabla\!\cdot\Big(\sum_g h_g \mathbf J_g\Big)
= -\nabla\!\cdot\mathbf q_{\mathrm{cond}} + Q_{\mathrm{rad}},
$$

expanded with the product rule and the species balances
$\partial_t\xi_j + \nabla\!\cdot(\mathbf v_s\xi_j) = \dot\omega_j - \nabla\!\cdot\mathbf J_j$
(with $\mathbf J_j = 0$ for condensed $j$), becomes

$$
\Big(\sum_j \xi_j c_{p,j}\Big)\Big(\partial_t T + \mathbf v_s\!\cdot\!\nabla T\Big)
= -\nabla\!\cdot\mathbf q_{\mathrm{cond}}
\underbrace{-\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T}_{S_{\mathrm{conv}}}
\underbrace{-\sum_r h_r(T)\,r_r}_{Q_{\mathrm{rxn}}}
+ Q_{\mathrm{rad}},
$$

where $h_r(T) = \sum_j \nu_{r,j} h_j(T)$ is the heat of reaction at the local
temperature. Three things follow:

1. **$S_{\mathrm{conv}}$ is exact**, with $\mathbf J_g$ the skeleton-relative flux
   the transport closures return. It is the only place the gas flux enters the
   energy equation.
2. **No further term exists when the heats of reaction are quoted at the reaction
   temperature** — the DSC-style, inverse-calibrated convention (§6.6.3). The
   composition-change enthalpy $\sum_j h_j\,\partial_t\xi_j$ of the conservative form
   has been absorbed exactly into $Q_{\mathrm{rxn}}$ and $S_{\mathrm{conv}}$. If the
   heats are instead quoted at a datum $T_{\mathrm{ref}}$, Kirchhoff's law
   $h_r(T) = h_r^0 + \sum_j\nu_{r,j}\Delta h_j(T,T_{\mathrm{ref}})$ adds
   $-\sum_j\Delta h_j(T,T_{\mathrm{ref}})\,\dot\omega_j$ over *every* component, a
   term the energy equation does not carry. Heats of reaction are therefore supplied
   at the reaction temperature (§6.6.3).
3. **The skeleton-following derivative $\partial_t T + \mathbf v_s\!\cdot\!\nabla T$ is
   the mesh-frame rate in ALE**, where the mesh moves with the skeleton (§3.6, §11).
   The temperature equation therefore carries no advection term of any kind in
   production, and the statement above holds unchanged on a fixed mesh with
   $\mathbf v_s = 0$.

This is ThermaKin Eq. 17 term for term (all-component storage, the convective term
$c_g\mathbf J_g\!\cdot\!\nabla T$, and the deformation term as the $\mathbf v_s\!\cdot\!\nabla T$
piece); Gpyro's `FULL_QSG` form is the same equation (§3.7).

> **Implementation notes.** `effective_heat_capacity` (every component) is in
> `src/Properties/effective_properties.jl`; the generated summation
> `_sum_heat_capacity` has no phase branch. The same value feeds the residual, the
> structured-Jacobian operators (through Dual evaluation), the park-and-mask
> temperature merge (`merge_temperatures_conservatively`,
> `src/Adaptivity/cell_merge.jl`) and the energy ledger (§16).

### 3.3.3 Conduction, radiation, and reaction-heat terms

**Conduction.** $\mathbf q_{\mathrm{cond}} = -k_{\mathrm{eff}}(T,\xi)\nabla T$ is
Fourier's law with the composition- and temperature-dependent effective conductivity
$k_{\mathrm{eff}}$ (mixing rules in §5). The minus sign in $-\nabla\cdot\mathbf q_{\mathrm{cond}}$
makes net conductive *outflow* cool the cell. At a face the flux is discretized with a
distance-weighted harmonic-mean conductivity (resistance-in-series), reducing to the
ordinary harmonic mean on a uniform mesh; the full derivation is in §7.

**Radiation.** $Q_{\mathrm{rad}}$ is the volumetric radiation absorption and is
entirely model-dependent. For `NO_RADIATION` it is zero; for `SURFACE_ABSORPTION` it is **also zero** —
the absorbed flux $\alpha\,q_{\mathrm{ext}}$ enters through the *surface energy
balance* (boundary Newton solve, §12.3.5), not as a volumetric deposit;
for `BEER_LAMBERT` it follows the in-depth attenuation $q''' = \alpha I(z)$ with
$I(z) = I_{\mathrm{surface}}\exp(-\int\alpha\,dz')$; for the experimental `P1_QUASI_STEADY` model (which `solve`
rejects; §8.12) it is the net
P1 source $\alpha(G-4\sigma T^4)$ computed from a separate quasi-steady radiation field.
The four models, their mutual exclusivity, and their numerics are the subject of §8.

**Reaction heat.** With the internal sign convention $h_r>0$ endothermic,

$$
Q_{\mathrm{rxn}} = -\sum_{r=1}^{N_R} h_r\, r_r,
$$

so an exothermic reaction ($h_r<0$, $r_r>0$) gives $Q_{\mathrm{rxn}}>0$ (heats the
cell) and an endothermic reaction ($h_r>0$) gives $Q_{\mathrm{rxn}}<0$ (cools it). The
units check: $h_r$ [J·kg⁻¹ of first reactant] $\times$ $r_r$ [kg·m⁻³·s⁻¹] gives
W·m⁻³. The implementing loop literally accumulates `Q -= heats[i] * rates[i]`. The heat
of reaction $h_r$ may itself be temperature- or state-dependent (e.g. a
moisture-sorption heat $h(T,\xi)$ via `StateDependentProperty`); the rate law $r_r$ is
the smoothly-gated Arrhenius form of §6.

> **Implementation notes.** $Q_{\mathrm{rxn}} = -\sum_r h_r r_r$ is
> `reaction_heat_source` in `src/Physics/kinetics.jl`, with the negative
> sign documented in-source ("Note the negative sign!"). Reaction rates and heats are
> evaluated together by `evaluate_reactions!` (same file).

### 3.3.4 The advection source $S_{\mathrm{conv}}$ and the surface gas enthalpy

**$S_{\mathrm{conv}}$ (always active).** The gas-advection source carries the sensible
enthalpy convected by the gas flux,

$$
S_{\mathrm{conv}} = -\sum_{g} c_{p,g}(T)\,\mathbf J_g\cdot\nabla T,
\qquad\text{discretely}\qquad
S_{\mathrm{conv},i} = -\sum_{g} c_{p,g}(T_i)\,\bar J_{g,i}\left(\frac{\partial T}{\partial z}\right)_i,
$$

where $\bar J_{g,i} = \tfrac12 (J_{g,f_L} + J_{g,f_R})$ is the cell-center average of the
two face fluxes and $(\partial T/\partial z)_i$ is the second-order opposite-distance
gradient (boundary cells blend in the boundary-face slope at the solved $T_s$; §7.4.2).
It is the only carrier of gas sensible enthalpy in the temperature equation.

**No surface gas-enthalpy term.** The conservative-form gas-enthalpy divergence splits
by the product rule as

$$
\nabla\cdot\Big(\sum_g \Delta h_g\,\mathbf J_g\Big)
= \sum_g c_{p,g}\,\mathbf J_g\cdot\nabla T + \sum_g \Delta h_g\,\nabla\cdot\mathbf J_g ,
$$

exactly even for temperature-dependent $c_{p,g}$, since $\nabla\Delta h_g = c_{p,g}(T)\,\nabla T$.
The first piece is $-S_{\mathrm{conv}}$; the second cancels exactly against the
transport part of the gas composition-change term ($\partial_t\xi_g = \dot\omega_g -
\nabla\cdot\mathbf J_g$, §3.3.2), whose reaction part is $Q_{\mathrm{rxn}}$ when the
heats of reaction are quoted at the reaction temperature. The gas enthalpy is
therefore carried entirely *volumetrically*, and the surface boundary condition is a
plain convective + radiative balance, $q_{\mathrm{BC}}(T_s) = q_{\mathrm{cond}}(T_s)$,
with **no** gas-enthalpy term. A surface term $\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$
would apply the *entire* outflow enthalpy on top of $\int S_{\mathrm{conv}}\,dz$,
double-counting roughly $5$–$10$ kW·m⁻² as a spurious sink. The gas flux reaches the
surface balance only through the blowing model of a `ConvectiveBC`, which reduces the
convective coefficient (§12.3.3). The surface energy balance is detailed in §12.

> **Implementation notes.** $S_{\mathrm{conv}}$ is `compute_gas_convective_source!`
> in `src/Physics/convection.jl`.

### 3.3.5 Why temperature, not enthalpy: trade-offs of the non-conservative form

Solving for $T$ (rather than an enthalpy density $\rho h$) keeps the primary unknown
directly observable, makes property evaluation and the surface Newton solve natural,
and is the FDS/Gpyro/ThermaKin tradition for solid-phase pyrolysis. The price is that
the *discrete* energy balance is only approximately conservative: it does not evolve the
composition-change enthalpy term $\sum_j h_j(T)\,\partial_t\xi_j$ as part of a conserved
enthalpy density (§3.3.2) — its reaction and gas-transport parts are re-expressed as
$Q_{\text{rxn}}$ and $S_{\text{conv}}$ — and (in ALE) the
volume-change contribution to storage is likewise outside the ledger. These
discrete-closure gaps are bounded and are made visible by the conservation diagnostics
(§16).

---

## 3.4 Gas-phase mass conservation and the pressure closure

### 3.4.1 The gas species balance and its flux

For each gas component the balance is the gas branch of §3.2.3,
$\partial_t \xi_j = -\nabla\cdot\mathbf J_j + \left.\partial_t\xi_j\right|_{\mathrm{rxn}}$.
The flux $\mathbf J_j$ closes by one of two constitutive models depending on the
transport mode (§3.5). In the diffusion-only (`Fickian`) mode the flux is the
volume-fraction-gradient Fickian form,

$$
J_j = -\,\lambda_{\mathrm{face}}\,
\frac{\xi_R T_R - \xi_L T_L}{T_{\mathrm{face}}\,d_{LR}},
\qquad T_{\mathrm{face}} = \tfrac12(T_L+T_R),
$$

with $\lambda_{\mathrm{face}}$ the distance-weighted harmonic mean of the cell gas-transfer
coefficients. The reason the gradient is of $\xi T$, not of $\xi$, is the ideal-gas
density dependence: for a gas $\rho_g\propto P/T$, so the mass-fraction gradient that
actually drives Fickian diffusion expands by the product rule as

$$
\xi_R T_R - \xi_L T_L = T_{\mathrm{face}}(\xi_R-\xi_L) + \xi_{\mathrm{face}}(T_R-T_L),
$$

cleanly separating ordinary (composition) diffusion from thermal (Soret) diffusion. The
derivation, the partial-pressure-gradient equivalence, and the boundary fluxes
are all in §9;
here the essential point for the governing system is simply that $\mathbf J_j$ is a
known constitutive function of the state.

> **Implementation notes.** The Fickian gas flux is `gas_flux` /
> `compute_gas_fluxes!` in `src/Physics/gas_transport.jl`. The harmonic-mean
> $\lambda_{\mathrm{face}}$ is `harmonic_mean_weighted` in
> `src/Properties/effective_properties.jl`.

### 3.4.2 Algebraic gas-pressure closure: the ideal-gas law in a porous medium

The gas pressure is not an independent evolved variable; it is closed algebraically
— an **exact equation-of-state evaluation**, not a quasi-steady approximation: the
gas storage term ∂ξ_g/∂t is fully retained in the species ODEs, and the summed
species equations reproduce total gas continuity with no spurious source. The closure is built
from the gas concentrations, temperature, and porosity by the ideal-gas law applied to
the **pore volume**. Summing partial pressures over gas components and dividing by the
porosity $\varphi$ (the gas occupies only the void fraction of the bulk volume),

$$
\boxed{\;
P = \frac{N(\xi,T)}{\varphi},
\qquad
N(\xi,T) = \sum_{j:\,\mathrm{is\_gas}(j)} \frac{\xi_j R_g T}{M_j} .
\;}
$$

Only components with `phase == GAS` and $M_j>0$ contribute to $N$; solids and liquids
are skipped. The division by $\varphi$ is physically essential: the same gas mass
concentration in a nearly-solid cell occupies a much smaller pore volume and therefore
exerts a higher pressure (confined gas). To keep the divisor well-defined as a cell
approaches fully solid, $\varphi$ is floored at $10^{-12}$, i.e. the code evaluates
$P = N/\max(\varphi,10^{-12})$.

The porosity is the **canonical intrinsic-density form** (overload note A8 of §2),

$$
\varphi = \operatorname{clamp}\!\Bigl(1 - \sum_{j\in\{\text{solid},\text{liquid}\}} \frac{\xi_j}{\rho_{i,j}},\;0,\;1\Bigr),
$$

where $\rho_{i,j}$ is the **intrinsic (skeletal, cell-wall) density** — the density of
the condensed material itself, ignoring pore volume — and gases do not enter the sum.
When a component is constructed without an explicit intrinsic density, $\rho_{i,j}$
falls back to its bulk density, and the component contributes no pore volume. (The
all-phase gas volume fraction $\sum_j v_j^{\text{gas}}$ of §5.1.2 is a different
quantity; the intrinsic-density form above is the porosity the property updates and the
pressure closure use.) The porosity derivatives the Jacobian
needs ($\partial\varphi/\partial\xi_j = -1/\rho_{i,j}$ for solid/liquid, $0$ for gas, in
the clamp-inactive interior) are obtained by forward-mode AD through `porosity`
itself (§5.3.2, §15).

Pressure matters only when Darcy flow is active (`DarcyFick` mode). In the diffusion-only
mode the pressure field is never assembled; the gas flux is purely Fickian and the
pressure closure is inert. The pressure derivatives $\partial P/\partial T$ and
$\partial P/\partial\xi_k$ (gas direct term plus the porosity-coupling term) that the
Jacobian needs are obtained by forward-mode differentiation of the closure itself —
see §9.6.

> **Implementation notes.** $P = N/\max(\varphi,10^{-12})$ is `ideal_gas_pressure`
> in `src/Physics/pressure.jl` (porosity floor `_PHI_FLOOR = 1e-12`; gas-only
> partial-pressure sum `_gas_partial_pressure_sum`). The
> canonical porosity is `porosity` in `src/Properties/mixture_rules.jl`; its
> derivatives are taken by AD through the same function (§5.3.2).

---

## 3.5 Gas-phase momentum: Darcy's law

### 3.5.1 Darcy velocity and the Darcy–Fick flux

When pressure-driven flow is significant (the material has finite permeability and the
problem is built in `DarcyFick` transport mode), the bulk gas velocity is given by
**Darcy's law**,

$$
\boxed{\; v_g = -\frac{\kappa}{\mu}\,\nabla P. \;}
$$

At a face the production path forms the per-cell **mobility** $\kappa/\mu$ (with
$\mu$ the material's pore-gas viscosity at the cell temperature, Sutherland's law for
air by default) and combines the two sides
with the distance-weighted harmonic mean (§9.5.2); if either side's mobility is
$\le 10^{-30}$ the face mobility is exactly $0$ (an impermeable cell is an
infinite series resistance). The pressure gradient is a central difference
$\nabla P = (P_R-P_L)/d_{LR}$:

$$
v_g = -\,\mathrm{mob}_{\mathrm{face}}\,\frac{P_R-P_L}{d_{LR}} .
$$

Gases flow *through* the matrix, so only solid/liquid components carry
permeability; `get_permeability` returns $0$ for a `GaseousComponent`.

The full gas flux in Darcy–Fick mode superposes advection and diffusion,

$$
J_j^{\mathrm{DF}} = \underbrace{\frac{\xi_j^{\mathrm{upwind}}}{\varphi_{\mathrm{face}}}\, v_g}_{\text{advection}}
\;+\; \underbrace{\Bigl(-\lambda_{\mathrm{face}}\frac{\xi_R T_R - \xi_L T_L}{T_{\mathrm{face}}\,d_{LR}}\Bigr)}_{\text{diffusion}},
\qquad
\xi_j^{\mathrm{upwind}} =
\begin{cases}\xi_{L,j}, & v_g\ge 0,\\ \xi_{R,j}, & v_g<0,\end{cases}
$$

a first-order upwind selection for the advected concentration; the $1/\varphi_{\mathrm{face}}$
converts the bulk concentration to the intrinsic pore density that the superficial
Darcy velocity transports (§9.4). All gas species are
carried equally by the bulk velocity $v_g$, and they share a single diffusive
transfer coefficient — there is no per-species molecular-weight scaling (§9.2).
The permeability and viscosity closures and the velocity's Jacobian treatment
are detailed in §9.

> **Implementation notes.** The production Darcy velocity (distance-weighted
> harmonic mobility) is computed in `update_darcy_properties!`
> (`src/Discretization/finite_volume.jl`); the boundary-face velocity comes from
> the shared kernel `boundary_darcy_velocity`
> (`src/Discretization/darcy_velocity.jl`). The combined flux is
> `gas_flux_darcy_fick` / `compute_gas_fluxes_darcy_fick!` in
> `src/Physics/gas_transport.jl` (advection of the intrinsic density ξ/φ_face).
> The viscosity is `effective_gas_viscosity(T, material)`, the material's
> `gas_viscosity` property with a positive floor (§5.9).

### 3.5.2 Where momentum is not solved

There is no separate momentum PDE, no inertial or viscous-stress terms in the gas, and
no two-temperature treatment. Darcy's law is an algebraic constitutive closure for the
bulk velocity given the (algebraically closed) pressure field of §3.4. This is the
standard porous-medium-pyrolysis assumption: pore Reynolds numbers are tiny and the
flow is creeping. The pressure boundary conditions (Dirichlet/Neumann/Convective) enter
through the Darcy velocity at the boundary faces, not through a pressure evolution
equation (§12).

---

## 3.6 Deformation: from ThermaKin's single-reference-point source to ALE mesh motion

### 3.6.1 The volumetric strain rate

As reactions consume condensed mass and (optionally) swell the matrix, the material
deforms. The local **volumetric strain rate** (relative rate of cell-volume change) is

$$
\boxed{\;
\theta_i = \frac{1}{V_i}\frac{dV_i}{dt}
= \sum_{j=1}^{N_C} \gamma_j\,\frac{1}{\rho_j(T_i)}\,\left.\frac{d\xi_j}{dt}\right|_i ,
\;}
$$

derived from the mixture continuity relation: when component $j$ (with bulk density
$\rho_j$) changes its mass concentration at rate $d\xi_j/dt$, it frees or claims volume
at rate $(1/\rho_j)\,d\xi_j/dt$, weighted by the swelling factor $\gamma_j$ that decides
whether the component participates in the matrix skeleton. Solids default to $\gamma=1$
(depletion shrinks the matrix); gases default to $\gamma=0$ (they escape without
changing matrix volume); intumescent gases may take $\gamma=1$ (gas creates new
volume); partial values $\gamma\in(0,1)$ are allowed. The implementing kernel
`volume_change_rate` sums exactly $\sum_j \gamma_j\,\dot\xi_j/\rho_j(T)$, skipping any
component with $\gamma_j\le 0$, and uses the *reaction* species source $d\xi_j/dt$ as the
driver. The absolute rate is $dV_i/dt = V_i\theta_i$.

> **Implementation notes.** $\theta = \sum_j \gamma_j\,\dot\xi_j/\rho_j(T)$ is
> `volume_change_rate` (`_volume_change_rate`) in `src/Physics/volume_change.jl`
> (with `if γ > 0` guards). It is the swelling-weighted per-component form, not
> $-(\sum_j\dot\xi_j)/\rho$.

### 3.6.2 The ThermaKin/Gpyro single-reference-point deformation source

In ThermaKin and Gpyro the same physics is expressed as a *deformation source* added to
the transport equations on a fixed (Eulerian) grid that is then re-mapped to a single
reference point. Schematically, ThermaKin introduces a convective-like term
$\int (1/\rho)(\partial\rho/\partial t)\,dx$ that displaces material relative to a
reference coordinate as the density changes; Gpyro tracks a per-cell solid fraction /
conversion and rescales control-volume thicknesses. Both encode the *same* local strain
$\theta$ but realize the resulting mesh deformation by re-referencing the grid.

### 3.6.3 Pyrolysis.jl's recasting via ALE mesh motion

Pyrolysis.jl instead makes the strain drive an **Arbitrary Lagrangian–Eulerian** mesh:
the nodes move with a velocity field $w$ chosen so that each cell's geometric volume
tracks its strain. Integrating the axial strain cumulatively from the pinned bottom
node ($w_1=0$) upward,

$$
w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i,
\qquad \Delta z_i = \frac{V_i}{A},
$$

where $\theta_A$ is the **column-global radial strain rate** subtracted to avoid
double-counting volume change that is instead expressed through the cross-section (see
below). The node positions then advance as $dz_i/dt = w_i$, and the geometric volume is
recomputed as $V_i = A(t)\,|z_{i+1}-z_i|$, which satisfies the discrete Geometric
Conservation Law $V_i^{\text{new}}-V_i^{\text{old}} = \Delta t\,(w_R-w_L)A$ to machine
precision. Because the mesh moves with the skeleton and every component — the pore
gas included, since its transport flux is relative to the skeleton — rides with it,
the mesh frame introduces **no advection term** in any equation; the only mesh-motion
term is the **dilation correction** $\Delta(d\xi_j/dt) = -\xi_j\theta_i$ applied to
every component, which keeps each cell's component masses on their reaction-and-flux
budget as its volume changes. The mesh kinematics, the dilation correction, and the
geometric conservation law are derived in full in §11 (The ALE Framework) and §10 (Volume Change &
Lateral Shrinkage). This chapter records only that the deformation physics is *the same* as
ThermaKin/Gpyro's; the implementation choice is a moving mesh rather than a
reference-point source.

### 3.6.4 Variable cross-section and the pyrolysis-progress state $\chi$

Lateral (radial) shrinkage or swelling is folded into the scalar cross-sectional area

$$
A(t) = A_0\cdot \mathrm{law}(\bar\chi),
$$

driven by the **mass-weighted column-average pyrolysis progress**

$$
\bar\chi = \frac{\sum_i m_{\mathrm{dry},i}\,\chi_i}{\sum_i m_{\mathrm{dry},i}},
\qquad
\frac{d\chi_i}{dt} = \frac{S_{\mathrm{dry},i}\,V_i}{m_{\mathrm{dry},i,\mathrm{init}}},
\quad
S_{\mathrm{dry},i} = \sum_{j\in\text{dry gas}} \left.\frac{d\xi_j}{dt}\right|_{\mathrm{rxn}}.
$$

Here $\chi_i$ is the auxiliary state of §3.1.2 (cumulative dry-pyrolysis-gas mass
released per unit initial dry-solid mass in the cell), and $m_{\mathrm{dry},i,\mathrm{init}}$
is the constant initial dry-solid mass populated at setup. Water vapour is excluded from
the dry-gas sum (it is sorption/desorption, not pyrolysis). The column-global radial
strain rate that feeds the axial mesh velocity above follows by the chain rule,

$$
\theta_A = \frac{1}{A}\frac{dA}{dt}
= \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\,\frac{d\bar\chi}{dt},
$$

with $\mathrm{law}'$ obtained by automatic differentiation and guards against
$\mathrm{law}(\bar\chi)\to 0$ (vanishing cross-section). The lateral-shrinkage law, its
identity/linear-shrink/swell variants, and the first-principles-vs-empirical discussion
are in §10.

> **Implementation notes.** The cumulative mesh velocity
> $w_{i+1}=w_i+(\theta_i-\theta_A)\Delta z_i$ is in `_compute_ale_mesh_velocity!`,
> `src/Residual/residual.jl` (with `w[1]=0`). The χ dynamics
> $d\chi_i/dt = S_{\mathrm{dry},i}V_i/m_{\mathrm{dry},i}$ is `_pack_chi_derivative!`;
> $\theta_A$ is `_compute_theta_A` (same file). The area update
> $A=A_0\,\mathrm{law}(\bar\chi)$ is `update_cross_section_area!` in
> `src/Problem/state_cache.jl`.

---

## 3.7 Comparison to ThermaKin2Ds, Gpyro, and FDS

Pyrolysis.jl is, by design, a synthesis of the three established condensed-phase
pyrolysis solvers. This section maps their governing forms onto ours and contrasts the
modeling choices. On first use we map each code's notation to the unified symbols of §2:
ThermaKin's stoichiometric coefficient $\theta_i^j \to \nu_{r,j}$; Gpyro/FDS
pre-exponential $Z \to A_r$, permeability $K\to\kappa$, kinematic viscosity
$\nu\to\mu/\rho_g$, volumetric reaction rate $\dot\omega''' \to r$, volume fraction
$X_\alpha\to v_j$; FDS solid absorption $\kappa_s\to\alpha$. **Heat-of-reaction sign:**
ThermaKin and Gpyro publish $h>0=$ exothermic; internally Pyrolysis.jl uses
$h>0=$ endothermic and carries the minus sign in $Q_{\mathrm{rxn}}=-\sum_r h_r r_r$, so a
ThermaKin/Gpyro exothermic $h_{\text{pub}}>0$ corresponds to our $h_r<0$.

### 3.7.1 Species / mass conservation

| Quantity | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS solid phase |
|---|---|---|---|---|
| Species state | mass concentration $\xi_j$ [kg·m⁻³] | mass concentration | $\bar\rho Y_i$ (mass per bulk volume) | $\rho_{s,\alpha}$ component densities |
| Solid transport | none (in-place; ALE carries) | none (in-place) | none; control-volume $\Delta z$ rescaled | none; layer thicknesses adjusted |
| Gas transport | $-\nabla\cdot\mathbf J_j$ (Fickian or Darcy–Fick) | $-\nabla\cdot\mathbf J_j$, $J=-\rho_g\lambda\,\partial(\xi_g/\rho_g)/\partial x$ | Darcian $\dot m'' = -(K/\nu)(\partial P/\partial z - g\rho_g)$ + pressure-evolution eq. | lumped gas release |
| Reaction source | $\sum_r\nu_{r,j}r_r$, single reactant, $\sum_{\text{prod}}\nu=1$ | $\sum_i\theta_i^j r_i$ (bimolecular $\xi_{\mathrm{C1}}\xi_{\mathrm{C2}}$ allowed) | $n$-th order in conversion, O₂ sensitivity | $\rho^n\,Y_{\mathrm O_2}^{n}$ |

All four share the central idea that condensed components are immobile and only gas
fluxes. ThermaKin's gas flux $J=-\rho_g\lambda\,\partial(\xi_g/\rho_g)/\partial x$ is the
mass-fraction-gradient form; Pyrolysis.jl's $\xi T$-gradient form (§3.4) is the same
physics written so the $\rho_g\propto P/T$ density dependence appears explicitly,
splitting ordinary and thermal diffusion. Gpyro additionally evolves the pressure field
through a gas-mass-balance PDE and uses Darcian mass flux; Pyrolysis.jl closes pressure
*algebraically* ($P=N/\varphi$) and uses Darcy's law only for the velocity. ThermaKin
supports bimolecular reactions; Pyrolysis.jl currently enforces a single reactant (§6).

### 3.7.2 Energy equation

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| Primary unknown | temperature $T$ | $T$ | $T$ | $T$ |
| Storage coefficient | all phases $\sum_j\xi_j c_{p,j}$ | all phases $\sum_j\xi_j c_{p,j}$ (Eq. 17) | matrix by default; all phases with `FULL_QSG` | matrix (no in-depth gas) |
| Gas sensible transport | $S_{\mathrm{conv}}=-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$ | convective term $c_g(\mathbf J_g\!\cdot\!\nabla T)$ | $\dot m'' c_{pg}\,\partial T/\partial z$ (thermal-equilibrium form) | boundary outflow |
| Conduction | Fourier, harmonic-mean face $k$ | Fourier | Fourier, Patankar harmonic interface | 1D Fourier |
| Reaction heat | $-\sum_r h_r r_r$ ($h>0$ endo) | $\sum h_i r_i$ ($h>0$ exo) | $\sum(\Delta H)\dot\omega'''$ | $\sum h_r\dot\omega'''$ |

Pyrolysis.jl's energy equation is ThermaKin Eq. 17 term for term: the all-phase
storage coefficient, the convective term $c_g(\mathbf J_g\cdot\nabla T)$ (our
$S_{\mathrm{conv}}$), and the deformation term, which the ALE mesh frame absorbs
(§3.3.2, §11.2.2). Gpyro's thermal-equilibrium form drops the gas storage
$\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default and restores it with
`FULL_QSG`, which makes its Eq. 65 the same equation; FDS carries no in-depth gas, so
the gas storage is identically zero there.
The one-dimensional Cartesian energy equation is otherwise identical in all four
codes; ThermaKin additionally supports a cylindrical
(radial) Laplacian $\frac1r\partial_r(rk\partial_r T)$, which Pyrolysis.jl approximates
only through the scalar lateral-shrinkage area $A(t)$ rather than resolving a radial
coordinate.

### 3.7.3 Control-volume conservation forms

Gpyro writes its discretization as control-volume balances of $\bar\rho\,\Delta z$ (bulk
mass per area), $\bar\rho Y_i\,\Delta z$ (per-species mass), gas mass, and energy, with
$\Delta z$ itself a state (variable control-volume thickness); its gas balance (Eq. 58)
is written on those deforming volumes with only the Darcy flux $\dot m''$ crossing
their faces. Pyrolysis.jl's ALE mesh motion (§3.6.3) is the moving-mesh analog of the
same frame: instead of rescaling $\Delta z$ inside fixed control volumes, the nodes
move so that $V_i = A|z_{i+1}-z_i|$ tracks the strain, every component including the
pore gas rides the cell, and the discrete GCL guarantees the volume bookkeeping is
exact (§11). ThermaKin's Eqs. 15 and 17 apply one deformation term to every
component and to the energy equation — again the same frame. FDS adjusts solid-layer
thicknesses according to a max/sum density-ratio swell/shrink rule (§10). All three
strategies enforce the same physical volume change; they differ in *bookkeeping*.

### 3.7.4 Reaction kinetics

FDS recovers Arrhenius parameters from a TGA reference temperature and peak rate
($A$, $E$ from $T_{\mathrm{peak}}$ and $r_{\mathrm{peak}}$) and uses rate laws of the form
$\rho^n Y_{\mathrm O_2}^{n} T^{n}$; Gpyro uses $n$-th order in conversion $\alpha_{\mathrm{conv}}$
with optional oxygen sensitivity $g(Y_{\mathrm O_2})$ and critical-temperature Heaviside
gates; ThermaKin permits bimolecular products $\xi_{\mathrm{C1}}\xi_{\mathrm{C2}}$.
Pyrolysis.jl uses the mass-concentration Arrhenius law $r = A\,e^{-E/R_gT}\prod_j\xi_j^{n_{r,j}}$
with smooth `tanh` temperature gates and a `tanh` depletion limiter, restricted to a
single reactant (§6). These are *kinetics* details; for the governing system the only
relevant fact is that $r_r$ and $\nu_{r,j}$ are known functions of $(\xi,T)$.

---

## 3.8 Closure inventory: which constitutive law closes each equation

The PDE system of §3.1.3 is incomplete until every flux and source is expressed in terms
of the state $(T,\xi)$ (and the auxiliary $z$, $\chi$). The table below inventories the
closures and points to the chapter where each is derived. This is the contract that makes
the residual a *pure function* of the state.

| Governing equation / term | Closure | Constitutive law | Chapter |
|---|---|---|---|
| Total density $\rho$ | exact | $\rho=\sum_j\xi_j$ | §3.2, §5 |
| Storage $\rho c_p^{\mathrm{eff}}$ | property + mixing | $\sum_j\xi_j c_{p,j}(T)$, every component | §4, §5 |
| Conductive flux $\mathbf q_{\mathrm{cond}}$ | Fourier + mixing | $-k_{\mathrm{eff}}(T,\xi)\nabla T$, harmonic-mean face $k$ | §5, §7 |
| Reaction rate $r_r$ | Arrhenius | $A_r e^{-E_r/R_gT}\prod_j\xi_j^{n_{r,j}}$, tanh gates/limiter | §6 |
| Species source $S_{j}$ | stoichiometry | $\sum_r\nu_{r,j}r_r$ | §6 |
| Reaction heat $Q_{\mathrm{rxn}}$ | enthalpy of reaction | $-\sum_r h_r(T,\xi)\,r_r$ | §6 |
| Radiation source $Q_{\mathrm{rad}}$ | radiation model | surface / Beer–Lambert / P1 | §8 |
| Gas advection $S_{\mathrm{conv}}$ | gas flux + $c_{p,g}$ | $-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$ | §7, §9 |
| Gas diffusive flux $\mathbf J_j^{\mathrm{diff}}$ | Fick (volume-fraction gradient) | $-\lambda_{\mathrm{face}}(\xi_R T_R-\xi_L T_L)/(T_{\mathrm{face}}d_{LR})$ | §9 |
| Gas advective flux $\mathbf J_j^{\mathrm{adv}}$ | upwind intrinsic density $\times$ Darcy | $(\xi_j^{\mathrm{upwind}}/\varphi_{\mathrm{face}})\,v_g$ | §9 |
| Pressure $P$ | ideal gas in pores | $N(\xi,T)/\varphi$ | §9 |
| Porosity $\varphi$ | intrinsic density | $\operatorname{clamp}(1-\sum_{S/L}\xi_j/\rho_{i,j},0,1)$ | §5 |
| Gas velocity $v_g$ | Darcy | $-(\kappa/\mu)\nabla P$, distance-weighted harmonic mobility $\kappa/\mu$, material $\mu(T)$ | §9 |
| Effective $k_{\mathrm{eff}},\lambda_{\mathrm{eff}},\kappa_{\mathrm{eff}},\varepsilon_{\mathrm{eff}},\alpha_{\mathrm{eff}}$ | mixing rules | parallel/series/weighted/Bruggeman/Carman–Kozeny | §5 |
| Volume-change rate $\theta$ | swelling continuity | $\sum_j\gamma_j\dot\xi_j/\rho_j(T)$ | §10 |
| Mesh velocity $w$ | cumulative strain | $w_{i+1}=w_i+(\theta_i-\theta_A)\Delta z_i$ | §11 |
| Cross-section $A(t)$ | lateral-shrinkage law | $A_0\,\mathrm{law}(\bar\chi)$ | §10 |
| Boundary fluxes $q_{\mathrm{BC}}, J_{j,s}, v_{g,s}$ | thermal/mass/pressure BCs | Newton surface balance, film resistance, Darcy-from-gradient | §12 |
| Initial state | initial conditions | $T_{\mathrm{initial}}(z)$, $\xi_{j,\mathrm{initial}}(z)$ | §12 |

With every entry of this table substituted, the right-hand side of the system in §3.1.3
is a fully determined function $\mathrm{rhs}(\boldsymbol u, t)$. The unified
`residual!` evaluates exactly this function, dispatching on the bundled
`SimulationMode` to assemble only the active terms (Eulerian vs ALE; Fickian vs
Darcy–Fick; with/without radiation; with/without $\chi$). The discretization of the
operators (finite-volume divergence, face averaging, flux limiters), the temporal
integration (stiff implicit Runge–Kutta with the structured $W = I/(\gamma\Delta t)-J$
matrix), and the Jacobian assembly are the subjects of §13–§15.

---

## 3.9 Limitations and open issues

The continuum model rests on several approximations that are exact only in limiting
regimes. The most important, with the section that quantifies or discusses each:

- **One temperature.** Matrix and in-pore gas are assumed in local thermal equilibrium;
  there is no two-temperature gas energy equation. This breaks for very fast, hot gas
  throughput, where the pore gas does not equilibrate with the walls within a cell.
- **Algebraic (not evolved) pressure.** Pressure is closed instantaneously by
  $P=N/\varphi$ (§3.4.2) with a $\varphi$ floor of $10^{-12}$; there is no transient
  gas-mass-balance pressure PDE as in Gpyro. Acoustic and fast compressibility transients
  are not represented.
- **First-order upwind advection** for the Darcy–Fick advective flux (§3.5.1); no TVD or
  MUSCL limiter on the gas advection. Over many steps, upwinding can push
  concentrations slightly non-physical if not damped by reaction/heat coupling.
- **Single reactant per reaction.** Multi-reactant (bimolecular) kinetics, supported by
  ThermaKin, is not implemented; the per-reactant tuple would need extending to carry a
  mass stoichiometry (§6). There is also no temperature-dependent $A$ or pressure
  dependence, and no distributed-activation-energy tooling.
- **Scalar lateral shrinkage.** Radial deformation is captured only through the single
  scalar $A(t)=A_0\,\mathrm{law}(\bar\chi)$ (§3.6.4); there is no resolved transverse
  coordinate, so genuinely 2D/axisymmetric effects (ThermaKin's cylindrical Laplacian) are
  not represented. The functional form of $\mathrm{law}(\chi)$ is user-supplied and
  largely empirical; a first-principles closure (e.g. percolation-based) is an open
  question (§10). Guarding against $\mathrm{law}(\bar\chi)\to 0$ (cross-section collapse) is
  ad hoc.
- **Discrete energy-conservation gaps.** The non-conservative temperature form does not
  carry the composition-change enthalpy $\sum_j h_j(T)\,\partial_t\xi_j$ as a conserved
  quantity (§3.3.2) and (in ALE) leaves volume-change storage outside the discrete energy
  ledger; these gaps are bounded and made visible by diagnostics (§16) but are not closed
  exactly.

These limitations are inherited by every downstream chapter; where a later chapter
relaxes or quantifies one of them (e.g. §16 for the energy/mass ledgers, §18 for
regime-validity checks), the cross-reference is given there.


---

# 4. Materials, Components, and Property Functions

This chapter specifies the *constitutive layer* of Pyrolysis.jl: how a pyrolyzing
substance is decomposed into discrete phase-typed **components**, how each
component carries its thermophysical **property functions**, how those functions
are evaluated and differentiated, and how components and reactions are assembled
into the static `Material` object that every physics kernel reads. The governing
PDE system (Chapter 3, *Governing Equations*) closes only once these laws are
supplied; the homogenization of per-component properties into mixture (effective)
properties is the subject of Chapter 5, *Effective Properties and Mixing Rules*,
to which this chapter forwards repeatedly. Throughout, the spatial coordinate is
$z$ (axial, through-thickness, $z=0$ substrate / bottom, $z=L$ exposed surface),
the primary species state variable is the **mass concentration** $\xi_j$ [kg·m⁻³],
and the heat-of-reaction sign convention is $h>0$ endothermic with the volumetric
source $Q_\mathrm{rxn} = -\sum_r h_r r_r$ (overload note H1 of Chapter 2). In the
kinetics passages of this chapter the index $i$ denotes a **reaction** and $j$ a
**component** (overload note I1); elsewhere $i$ is a cell index.

---

## 4.1 The phase model: solids, liquids, and gases

### 4.1.1 Why phase is a *type*, not a field

A pyrolyzing material is a heterogeneous mixture of condensed-phase species
(virgin polymer, char, moisture) and gaseous species (pyrolyzate, water vapor,
oxidizer). The three phases obey structurally different conservation laws
(Chapter 3): condensed species are **Lagrangian** — they do not transport
relative to the matrix, so their concentration changes only through reaction —
whereas gaseous species are **Eulerian** — they flow through the porous matrix
and therefore carry a flux divergence in addition to a reaction source. The
gas-vs-condensed dispatch appears in the species residual:

$$
\frac{\partial \xi_j}{\partial t} =
\begin{cases}
-\,\nabla\!\cdot\mathbf{J}_j \;+\; \left.\dfrac{\partial \xi_j}{\partial t}\right|_\mathrm{rxn}
   & j \text{ is gas},\\[2.2ex]
\left.\dfrac{\partial \xi_j}{\partial t}\right|_\mathrm{rxn}
   & j \text{ solid or liquid.}
\end{cases}
\tag{4.1}
$$

Because this branch is taken inside the innermost residual loop, over every cell
and every component, it must not incur a runtime type check. Pyrolysis.jl
encodes the phase **at the type level**: there are three concrete subtypes of
`AbstractComponent`, and the phase predicate is a method that dispatches on the
singleton type, so the compiler eliminates the branch entirely at code
generation. The runtime `Phase` enum (`SOLID=1`, `LIQUID=2`, `GAS=3`) exists only
for printing, logging, and comparison; it is *not* consulted in hot paths.

```julia
@enum Phase begin SOLID = 1; LIQUID = 2; GAS = 3 end
```

> **Implementation note.** Phase types and predicates are defined in
> `SolidComponent`, `LiquidComponent`, `GaseousComponent`, and `is_solid` /
> `is_liquid` / `is_gas` / `phase` in `src/Materials/components.jl`. The
> compile-time-unrolled application of Eq. (4.1) is
> `apply_species_residual_unrolled!`, `src/Materials/materials.jl`, which emits
> one statically dispatched `is_gas(material.components[j])` branch per component
> — folded to a constant at codegen.

### 4.1.2 The three component types

**`SolidComponent`** is the matrix-forming condensed phase. It carries seven
property-function fields, an optional constant surface absorptivity, and two
scalar (`Float64`) fields:

```julia
struct SolidComponent{Pρ,Pc,Pk,Pλ,Pκ,Pε,Pa,Pρi} <: AbstractComponent
    name::Symbol
    density::Pρ                # bulk density ρ_j           [kg/m³]
    heat_capacity::Pc          # c_{p,j}(T)                 [J/(kg·K)]
    thermal_conductivity::Pk   # k_j(T)                     [W/(m·K)]
    gas_transfer_coeff::Pλ     # λ_j(T)                     [m²/s]
    permeability::Pκ           # κ_j(T)                     [m²]
    emissivity::Pε             # ε_j(T)                     [-]
    surface_absorptivity::Pa   # a_j, Nothing = Kirchhoff   [-]
    absorption_coeff::Float64  # α_j (mass basis)           [m²/kg]
    swelling_factor::Float64   # γ_j                        [-]
    intrinsic_density::Pρi     # skeletal density ρ_{i,j}   [kg/m³]
end
```

The seven property fields are each parameterized by their own type (`Pρ`, `Pc`, …),
so a component built with a constant conductivity and a polynomial heat capacity
is a *distinct concrete type* from one built with two constants. This preserves
type stability and allocation-free evaluation through the heterogeneous component
tuple. The emissivity `ε_j(T)` is a property function like `k_j(T)` (the `ε`
keyword accepts the same forms) and is evaluated at the cell temperature by the
mixing rule (§5.8). `surface_absorptivity` is the optional constant absorptivity
`a_j` for *incident* (heater-spectrum) radiation, `Pa <: Union{Nothing,Float64}`:
`Nothing` selects Kirchhoff's law `a_j = ε_j(T)` at the type level, so the
absorptivity mixture (§5.8.2) is a static sum. The two scalar fields (`α_j`, `γ_j`)
are stored as plain `Float64` because they are not given temperature dependence
in the current model.

**`LiquidComponent`** has the *identical* field layout to `SolidComponent`
(bulk and intrinsic density, the five property functions, the optional surface
absorptivity, the two scalars). It
is the matrix-forming condensed phase used for moisture / liquid water in
hygroscopic materials. Its convenience constructor differs only in defaults
(see §4.1.5): $\lambda$ defaults to $2\times10^{-5}$ m²·s⁻¹, $\varepsilon$ to
0.95 (condensed liquids are strong emitters — water ≈ 0.96),
$\gamma$ to 0 (volumeless bound moisture, ThermaKin `L SWELLING: 0`), and
$\kappa$ to 0 (a liquid normally fills pore space and blocks gas flow).

**`GaseousComponent`** represents species that flow *through* the matrix:

```julia
struct GaseousComponent{Pρ,Pc,Pk,Pλ,Pε} <: AbstractComponent
    name::Symbol
    molar_mass::Float64        # M_j                        [kg/mol]
    density::Pρ                # ρ_g(T), default ideal gas  [kg/m³]
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    emissivity::Pε             # ε_j(T)
    absorption_coeff::Float64
    swelling_factor::Float64
end
```

A gas has **no permeability field, no intrinsic density, and no surface
absorptivity** (`get_absorptivity` always returns its `ε_j(T)`) — it is not part
of the matrix skeleton. It carries instead a `molar_mass` $M_j$, used by the ideal-
gas density law (§4.3), the ideal-gas pressure closure (Chapter 9), and the
Chapman–Enskog diffusivity estimate (§5.10).

### 4.1.3 Compile-time phase accessors

All cross-phase property access goes through dispatched accessors that return
phase-appropriate defaults, so a kernel may ask any component for any quantity
without a runtime field-existence test:

| Accessor | Solid / Liquid | Gas |
|---|---|---|
| `is_solid(c)` | `true` / `false` | `false` |
| `is_liquid(c)` | `false` / `true` | `false` |
| `is_gas(c)` | `false` | `true` |
| `phase(c)` | `SOLID` / `LIQUID` | `GAS` |
| `molar_mass(c)` | `0.0` | `c.molar_mass` |
| `get_permeability(c, T)` | `c.permeability(T)` | `0.0` |
| `has_permeability(c)` | `true` | `false` |
| `get_intrinsic_density(c, T)` | `c.intrinsic_density(T)` | `0.0` |

Each is `@inline` and resolves at compile time on the singleton type; none branch
on a stored field. The gas's `get_permeability` and `get_intrinsic_density`
return exactly `0.0`, which is the algebraically correct value for the matrix
mixing rules (a gas contributes nothing to permeability or to the solid
skeleton).

> **Implementation note.** Accessors `get_intrinsic_density`, `get_permeability`,
> `has_permeability`, `get_emissivity`, `get_absorptivity`,
> `is_solid/is_liquid/is_gas`, `phase`, `molar_mass`, all in
> `src/Materials/components.jl`. Kernels read the remaining property fields
> (`density`, `heat_capacity`, …) directly.

### 4.1.4 Bulk versus intrinsic (skeletal) density

The single most important constitutive distinction at the component level is
between two densities (overload note A7 of Chapter 2):

- **Bulk (pure-phase) density** $\rho_j$ [kg·m⁻³] — the mass of pure phase $j$ per
  unit of the *bulk volume it would occupy in the mixture*. This is the density
  used by every mixing rule, by the mixture-density formula, and by the volume-
  change kernels (Chapter 10). Field `density`.

- **Intrinsic / skeletal density** $\rho_{i,j}$ [kg·m⁻³] — the density of the
  cell-wall material itself, ignoring any pore volume internal to the component.
  This is used by **one formula only**, the porosity (Chapter 5, §5.3):

  $$
  \phi = \mathrm{clamp}\!\left(1 - \sum_{j\in\{S,L\}} \frac{\xi_j}{\rho_{i,j}},\; 0,\; 1\right).
  \tag{4.2}
  $$

  Field `intrinsic_density`.

When a component is constructed without an explicit intrinsic density, the field
defaults to the bulk density, $\rho_{i,j}=\rho_j$: the pure phase is assumed to
have no internal pores. The resolution logic enforces that a
caller pass *either* `ρ_intrinsic` *or* a pure-phase porosity $\phi$ (from which
$\rho_{i,j}=\rho_j/(1-\phi)$ is computed, requiring a constant bulk $\rho_j$),
never both:

```
1. both ρ_intrinsic and φ given        → ArgumentError
2. ρ_intrinsic given                    → to_property(ρ_intrinsic)
3. φ given (requires ρ::Real, 0≤φ<1)    → ConstantProperty(ρ / (1 − φ))
4. neither                              → fall back to bulk ρ
```

> **Implementation note.** `_resolve_intrinsic_density`, `src/Materials/components.jl`.
> The `φ` path requires `raw_ρ isa Real` and `0 ≤ φ < 1`, throwing an
> `ArgumentError` otherwise. Gas components do not store an intrinsic density;
> `get_intrinsic_density(::GaseousComponent, ::Real) = 0.0`, which is consistent
> with Eq. (4.2) summing only over solid/liquid.

### 4.1.5 The swelling factor

The swelling factor $\gamma_j\in[0,1]$ controls whether changes in component $j$'s
concentration change the *bulk volume* of the cell. Physically it answers: when
this component is created or destroyed, does the material expand/shrink, or does
the mass simply enter/leave the existing pore space?

- $\gamma_j = 1$ (default for solids): full contribution. Depletion of the solid
  shrinks the bulk volume; this is the charring/recession behavior.
- $\gamma_j = 0$ (default for gases): no contribution. Gas fills existing pores
  and escapes without expanding the bulk.
- $\gamma_j \in (0,1)$: partial contribution, the intumescent case, where evolved
  gas creates new bulk volume.

$\gamma_j$ enters two places: the mixture-density formula (§4.2.2, where a gas with
$\gamma_j=0$ contributes zero specific volume) and the volumetric strain rate
$\theta = \sum_j \gamma_j \rho_j^{-1}\,\mathrm{d}\xi_j/\mathrm{d}t$ that drives the
ALE mesh motion (Chapter 10). It also weights the **matrix volume fractions**
used for conductivity/emissivity (Chapter 5, §5.1): every component carries
weight $\gamma_j\,\xi_j/\rho_j$ there, so a $\gamma_j = 0$ component (the
canonical gas) carries no weight in the matrix radiative surface.

> **Why the cap $\gamma_j \le 1$ is a model decision, not a missing feature.**
> Construction rejects $\gamma_j > 1$ (`_validate_swelling_factor`) because
> the *same* $\gamma$ weight feeds three mixing rules with volume-fraction
> semantics — the strain kernel, the mixture density, and the $k/\varepsilon$
> matrix volume fractions — and a weight above 1 corrupts all three (volume
> fractions would not sum consistently). Expansion beyond the reactant's
> volume is already expressible the physical way: give the product a **low
> bulk density** ($V_f/V_0 = \sum_j y_j\gamma_j\rho_A/\rho_j$ — the V-CPL-03
> P1 verification permutation demonstrates ×2 growth), with $\gamma \in
> (0,1]$ tuning how much of that specific volume the matrix keeps. A
> *decoupled* intumescent expansion factor (strain contribution independent
> of the mixing weight) would be a separate modeling design; no current use
> case demands it.

---

## 4.2 Component property fields and their roles

### 4.2.1 Inventory

Each property field of a component, its symbol, its SI units, and where it is
consumed downstream:

| Field | Symbol | Units | Phase | Primary consumer |
|---|---|---|---|---|
| `density` | $\rho_j$ | kg·m⁻³ | all | mixing rules, $\rho$ (§4.2.2), volume change (Ch. 10) |
| `intrinsic_density` | $\rho_{i,j}$ | kg·m⁻³ | S, L | porosity $\phi$ (§4.2.3, Eq. 4.2) |
| `heat_capacity` | $c_{p,j}(T)$ | J·kg⁻¹·K⁻¹ | all | $\rho c_p^\mathrm{eff}$, advection source, park-merge enthalpy (Ch. 5, 7, 14) |
| `thermal_conductivity` | $k_j(T)$ | W·m⁻¹·K⁻¹ | all | $k_\mathrm{eff}$ (Ch. 5), Fourier flux (Ch. 7) |
| `gas_transfer_coeff` | $\lambda_j(T)$ | m²·s⁻¹ | all | $\lambda_\mathrm{eff}$ (Ch. 5), gas diffusion (Ch. 9) |
| `permeability` | $\kappa_j(T)$ | m² | S, L | $\kappa_\mathrm{eff}$ (Ch. 5), Darcy velocity (Ch. 9) |
| `emissivity` | $\varepsilon_j(T)$ | – | all | $\varepsilon_\mathrm{eff}$ (Ch. 5), radiative BC (Ch. 12) |
| `surface_absorptivity` | $a_j$ (or `nothing` ⇒ $\varepsilon_j(T)$) | – | S, L | $a_\mathrm{eff}$ (§5.8.2), incident-flux BC (Ch. 12), Beer–Lambert surface transmissivity (Ch. 8) |
| `absorption_coeff` | $\alpha_j$ | m²·kg⁻¹ | all | $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ (Ch. 5, 8) |
| `swelling_factor` | $\gamma_j$ | – | all | $\rho$, strain rate $\theta$ (Ch. 10) |
| `molar_mass` | $M_j$ | kg·mol⁻¹ | G | $\rho_g$ (§4.3), $P$ (Ch. 9), Chapman–Enskog estimate (§5.10) |

Note that $\alpha_j$ is a **mass-basis** absorption coefficient [m²·kg⁻¹]; it is
multiplied by concentration to give the volumetric (extinction) absorption
coefficient $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ [m⁻¹] used in Beer–Lambert
and P1 radiation (Chapter 8). A solid usually carries a nonzero $\lambda_j$
in practice (e.g. $2\times10^{-5}$ m²·s⁻¹; the `SolidComponent` default is 0, and
`validate_problem` warns when a gas-producing Fickian material leaves it there,
§5.6): although solids do not themselves
diffuse, their $\lambda_j$ contributes to the *mixture* gas-transfer coefficient
through the volume-weighted mixing rule, matching the ThermaKin convention.

### 4.2.2 Mixture density with per-component swelling

The mixture density combines per-component bulk densities through specific
volumes (volume per unit mass), every component weighted by its own swelling
factor (the universal γ semantics shared with θ and the mixing-rule volume
fractions, §5.1.3):

$$
\frac{1}{\rho}
= \sum_{j} \gamma_j\,\frac{Y_j}{\rho_j},
\qquad Y_j = \frac{\xi_j}{\rho_\mathrm{tot}},\quad
\rho_\mathrm{tot}=\sum_k \xi_k.
\tag{4.3}
$$

Here $Y_j$ is the mass fraction. A non-swelling gas ($\gamma_j=0$) contributes
**zero** specific volume — its mass is "free" in the pore space and does not
change the bulk volume — whereas a swelling gas ($\gamma_j=1$) contributes its
full specific volume like a solid. This is the discrete realization of the
swelling semantics of §4.1.5.

> **Important — correction to a common gloss.** The swelling weight in Eq. (4.3)
> is applied **per component**, with each component's *own* `swelling_factor` —
> not an averaged $\gamma_\mathrm{avg}$ multiplying the gas sum, and not gases
> only. The source evaluates `γ_j · Y_j/ρ_j` term by term for every component
> (a γ = 0 term contributes nothing).

> **Implementation note.** `mixture_density` and the inline, allocation-free
> kernel `_specific_vol_with_swelling_direct`,
> `src/Properties/mixture_rules.jl`. The
> inline form computes $Y_j=\xi_j\,/\rho_\mathrm{tot}$ on the fly to avoid
> allocating an intermediate mass-fraction tuple. Degeneracy guards: if
> $\rho_\mathrm{tot}<10^{-20}$ or the accumulated specific volume
> $<10^{-20}$, the function returns `0.0`.

### 4.2.3 Porosity and its single use

Porosity $\phi$ (Eq. 4.2) is the void fraction of the bulk *not* occupied by the
solid/liquid skeletal frame. It is computed from **intrinsic** densities and is
the canonical definition $\phi = \mathrm{clamp}(1-\sum_{j\in\{S,L\}}\xi_j/\rho_{i,j},0,1)$.
Its derivatives are exactly (in the unclamped interior; at an active clamp the
derivative is 0, and the AD path used in production returns exactly that)

$$
\frac{\partial \phi}{\partial \xi_j} =
\begin{cases}
-\,1/\rho_{i,j} & j\in\{S,L\},\\
0 & j\in G.
\end{cases}
\tag{4.4}
$$

A gas does not enter $\phi$ because it occupies the void that the solid frame
leaves; adding gas mass does not displace skeletal material. Porosity feeds (i)
the ideal-gas pressure closure $P=N(\xi,T)/\phi$ (Chapter 9) and (ii) the
Carman–Kozeny permeability closure (Chapter 5, §5.7). The all-phase gas volume
fraction $\sum_{j\in G} v_j$ is a different quantity; the **intrinsic-density form
Eq. (4.2) is canonical** and is what `porosity` computes and `update_darcy_properties!` uses
(overload note A8 of Chapter 2).

> **Implementation note.** Canonical porosity: `porosity` →
> `_porosity_void_frac`, in `src/Properties/mixture_rules.jl` (derivatives are
> taken by forward-mode AD through `porosity` itself; §5.3.2). In
> `_porosity_void_frac` a component
> contributes $\xi_j/\rho_{i,j}$ only when `is_solid(c) || is_liquid(c)` *and*
> $\rho_{i,j}>0$. Full development is in Chapter 5.

---

## 4.3 Ideal-gas density and the gaseous-component default

By default a `GaseousComponent` constructed via the convenience constructor is
given a density obeying the ideal-gas law at the **reference pressure**:

$$
\rho_g(T) = \frac{P_\mathrm{ref}\,M_j}{R_g\,T},
\tag{4.5}
$$

with $P_\mathrm{ref}=101325$ Pa, $R_g=8.31446261815324$ J·mol⁻¹·K⁻¹, and $M_j$
the gas molar mass. The constructor wraps this in a `CallableProperty`:

```julia
ρ_func = CallableProperty(T -> P_REF * M_float / (R_GAS * T))
```

Two modeling points deserve emphasis:

1. **The density used for mixing/volume is at $P_\mathrm{ref}$, not the local
   pressure $P$.** The component's `density(T)` is consulted for volume fractions,
   specific volumes, and the mixture density. The *local* pressure-dependent gas
   density (which matters for the pore-pressure closure $P=N/\phi$) is handled
   separately by the ideal-gas formula in the Darcy machinery (Chapter 9), not by
   the component's `density` field. The component-level density at $P_\mathrm{ref}$
   is what enters the **specific volume** in mixing rules. In the
   swelling-weighted bases (mixture density, conductivity/emissivity weights)
   the gas's contribution is usually gated to zero by $\gamma_j=0$ (Eq. 4.3);
   in the **all-phase** basis used for λ mixing (§5.6) the gas volume is
   present regardless of γ and is over-weighted by ~$P/P_\mathrm{ref}$ in a
   pressurized cell — acceptable at the ≤ atm-scale overpressures of the
   validated cases. The pressure closure $P=N/\phi$ and the porosity are
   unaffected (they never use this density).

2. **The user may override the gas density** by passing any property function
   (constant, polynomial, table, callable) — Eq. (4.5) is merely the default
   produced when one constructs a gas with `GaseousComponent(:name; M=…, c=…, k=…,
   λ=…)`.

> **Implementation note.** `GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0,
> γ=0.0)` constructor, `src/Materials/components.jl`. Defaults for a gas are
> $\varepsilon=0$, $\alpha=0$, $\gamma=0$ (a gas escapes without swelling and is
> not part of the radiating surface).

---

## 4.4 Property-function types

A *property function* is a callable object $p$ such that $p(T)$ (or, for
state-dependent properties, $p(T,\xi)$) returns the property value at the given
state. All eight concrete types subtype `AbstractPropertyFunction` and implement
both the `(p)(T)` and `(p)(T, ξ)` call signatures (the latter ignoring $\xi$ for
all but `StateDependentProperty`). This uniformity lets a single component field
hold any property law, and lets the kinetics evaluator pass $\xi$ to a heat-of-
reaction property without knowing whether it is composition-dependent.

**The property bracket.** Every temperature-dependent law evaluates on the
clamped argument $\tilde T = \operatorname{clamp}(T, T_{\min}, T_{\max})$ with
$[T_{\min}, T_{\max}] = [200, 3000]$ K (`T_PROPERTY_MIN`, `T_PROPERTY_MAX`,
`src/Materials/constants.jl`; the same bracket clamps the surface-temperature
Newton iterate, §12). Inside the bracket the clamp is the identity with unit
derivative, so nothing changes for any physical state; outside it the law sees
the nearest bracket end. A Newton trial state with a nonphysical (even
negative) temperature therefore produces a failed Newton step, not a
`DomainError` out of a fractional power or a fitted polynomial's far
extrapolation. No derivative is hand-coded: the Jacobian and sensitivity passes
evaluate each law on `ForwardDiff.Dual` arguments (§4.5), so the derivatives
quoted below are what automatic differentiation of the law returns — and, through
the clamp, zero outside the bracket.

The component constructors sample every law on the bracket at 1 K steps and
reject a component whose bulk or intrinsic density, heat capacity or thermal
conductivity is not positive there (a gas may have $c_p = 0$), whose gas
transfer coefficient or permeability is negative, or whose emissivity leaves
$[0, 1]$. The error names the first offending temperature: a linear fit
$k = 0.3 - 4.2\times10^{-4}\,T$ is rejected at 715 K because it turns a cell
into a perfect insulator once it reaches that temperature. Hold a fit constant
outside its measured range (`T -> law(clamp(T, T_lo, T_hi))`) or pass
`validate_properties = false` for laws built programmatically from search
parameters, where an unphysical value should reach the objective function
instead of throwing.

### 4.4.1 `ConstantProperty`

$$
p(T) = \texttt{value}, \qquad \frac{\mathrm{d}p}{\mathrm{d}T}=0.
$$

The workhorse for properties given a single number. It has a *fast path* in the
parameter-evaluation interface (§4.7) that bypasses the dynamic-dispatch hop.

### 4.4.2 `LinearProperty`

$$
p(T) = a + bT, \qquad \frac{\mathrm{d}p}{\mathrm{d}T}=b.
$$

Fields `a` (intercept at $T=0$) and `b` (slope). Constructed automatically from a
two-tuple via `to_property((a,b))`.

### 4.4.3 `PolynomialProperty`

$$
p(T) = \sum_{m=0}^{N-1} c_m\, T^{m} = c_0 + c_1 T + c_2 T^2 + \cdots,
$$

evaluated with `evalpoly` (Horner's scheme) for numerical stability and to avoid
recomputing powers. Forward-mode AD through Horner's scheme returns the exact
derivative

$$
\frac{\mathrm{d}p}{\mathrm{d}T} = \sum_{m=1}^{N-1} m\,c_m\,T^{m-1}
= c_1 + 2c_2 T + 3c_3 T^2 + \cdots .
\tag{4.6}
$$

### 4.4.4 `ArrheniusProperty`

$$
p(T) = A\,\exp\!\left(-\frac{E}{R_g\,T}\right),
\qquad
\frac{\mathrm{d}p}{\mathrm{d}T} = p(T)\,\frac{E}{R_g\,T^2}.
\tag{4.7}
$$

Fields `A` (pre-exponential), `E` (activation energy [J·mol⁻¹]), and `T_ref` (a
reference temperature, defaulting to $T_\mathrm{ref}=298.15$ K, carried for
normalization purposes but not used in the value itself). The law shares its
functional form with the reaction-rate constant (§4.5.3, Chapter 6).

### 4.4.5 `TableProperty`

Tabulated $(T_m, v_m)$ pairs with **linear interpolation** and **flat
extrapolation**. The constructor validates that the temperatures are strictly
increasing (`error` otherwise). For a query $T_q$:

$$
p(T_q) =
\begin{cases}
v_1 & T_q \le T_1,\\[1ex]
v_m + \dfrac{T_q - T_m}{T_{m+1}-T_m}\,(v_{m+1}-v_m)
   & T_m \le T_q < T_{m+1},\\[2ex]
v_N & T_q \ge T_N.
\end{cases}
\tag{4.8}
$$

The interval search is **linear for $N\le 8$** (fewer comparisons than a binary
search on small tables) and **binary for $N>8$** (logarithmic). Under AD the
derivative is the slope of the segment containing $T_q$ (the right-hand segment at
a knot) and zero in the flat extrapolation; it jumps at every knot.

### 4.4.6 `CallableProperty`

Wraps an arbitrary callable $f(T)$. Its value is $p(T)=f(\tilde T)$ on the
clamped argument of the property bracket; its derivative is whatever AD
propagates through $f$, so $f$ must accept `ForwardDiff.Dual` arguments (no
`Float64` annotations or conversions inside). This is the type produced by
`to_property(f)` for a function `f` and is how the ideal-gas density default
(Eq. 4.5) is stored.

### 4.4.7 `SutherlandProperty`

Sutherland's law for a gas transport property, carrying its own three constants:

$$
p(T) = \mu_{\text{ref}}
\left(\frac{T}{T_{\text{ref}}}\right)^{3/2}
\frac{T_{\text{ref}}+S}{T+S}.
$$

`SutherlandProperty()` is the canonical **air** triple
($\mu_{\text{ref}}=1.716\times10^{-5}$ Pa·s, $T_{\text{ref}}=273.15$ K,
$S=110.4$ K) and is the default `Material.gas_viscosity` (§5.9), so a material
that says nothing about its pore gas uses air. Its derivative, as AD returns it, is

$$
\frac{\mathrm{d}p}{\mathrm{d}T} = p(T)\left(\frac{3}{2T} - \frac{1}{T+S}\right),
$$

(the AD-cleanliness of the type is tested in
`test/materials/test_gas_viscosity_property.jl`).

Like every property law it evaluates on the clamped argument $\tilde T$ of the
property bracket. The positivity floor `MU_FLOOR` is applied once, centrally, in
`effective_gas_viscosity(T, material)` (§5.9), so that every curve a material may
carry is protected rather than just this one.

### 4.4.8 `StateDependentProperty`

A property that depends on **both** temperature and the full concentration tuple,
$p(T,\xi)=f(T,\xi)$. Calling it with $T$ alone raises an error (state is
mandatory). Its canonical use is a moisture-dependent heat of sorption: the heat
released/absorbed when water desorbs depends on the current moisture content, so

$$
h(T,\xi) = f\!\left(T,\; \mathrm{MC}(\xi)\right),
$$

where $\mathrm{MC}$ is, e.g., the moisture content reconstructed from the
concentration tuple. The kinetics evaluator always passes $\xi$ to the heat
property (§4.5, Chapter 6), so a `StateDependentProperty` heat is evaluated with
the live composition.

Its derivatives in $T$ and in each $\xi_j$ are propagated by AD through $f$
(§4.5), so $f$ must accept a `Dual` temperature and a tuple of `Dual`
concentrations.

### 4.4.9 `to_property` conversions

The convenience converter normalizes user inputs to property functions, so that
component constructors and reaction constructors accept either property objects
or bare values:

| Input | Result |
|---|---|
| an `AbstractPropertyFunction` | passes through unchanged |
| a `Real` number | `ConstantProperty(Float64(x))` |
| a `Tuple{Real,Real}` $(a,b)$ | `LinearProperty(a, b)` |
| a `Function` | `CallableProperty(f)` |

> **Implementation note.** All property types, their call methods,
> `linear_interpolate`, and `to_property` are in
> `src/Materials/property_functions.jl`.

---

## 4.5 Property derivatives and their role in the Jacobian

Implicit time integration (Chapter 15) requires the Jacobian
$\partial(\text{RHS})/\partial u$. Because effective properties depend on $T$ and
$\xi$, and the RHS depends on effective properties, the per-component derivatives
$\mathrm{d}p_j/\mathrm{d}T$ and $\partial p_j/\partial \xi_k$ propagate through the
chain rule into the structured Jacobian (Chapter 5, §5.15 for the mixing-rule
chain; Chapter 15 for assembly). The property types carry no derivative methods.
Each Jacobian operator either runs ForwardDiff over a stencil-local kernel that
evaluates the property laws on `Dual` arguments, or implements an analytical
override tested against a ForwardDiff reference (Chapter 15). Either way the
property derivative is the derivative of the law as written.

### 4.5.1 Temperature derivative $\mathrm{d}p/\mathrm{d}T$

What forward-mode AD returns for each type, inside the property bracket:

| Type | $\mathrm{d}p/\mathrm{d}T$ |
|---|---|
| `ConstantProperty` | $0$ |
| `LinearProperty` | $b$ |
| `PolynomialProperty` | $\sum_{m\ge 1} m c_m T^{m-1}$ (Eq. 4.6) |
| `ArrheniusProperty` | $p(T)\,E/(R_g T^2)$ (Eq. 4.7) |
| `TableProperty` | slope of the active segment; $0$ in the flat extrapolation |
| `CallableProperty` | $f'(T)$, as propagated through `f` |
| `SutherlandProperty` | $p(T)\,(3/(2T) - 1/(T+S))$ |
| `StateDependentProperty` | $\partial f/\partial T$ at fixed $\xi$ |

These values are exact derivatives of the evaluated law, with no truncation
error. A `TableProperty` has a discontinuous derivative at every knot, which the
Newton iteration sees as a kink; a smooth law (polynomial, Arrhenius, Sutherland)
avoids it.

### 4.5.2 Concentration derivative $\partial p/\partial \xi_j$

For every property type *except* `StateDependentProperty`, the two-argument call
`p(T, ξ)` ignores $\xi$, so $\partial p/\partial \xi_j = 0$ exactly. For a
`StateDependentProperty` the derivative is $\partial f/\partial \xi_j$, propagated by
AD through `f(T, ξ)` with the concentration tuple carrying `Dual` entries. The
kinetics evaluator passes the live $\xi$ to every heat-of-reaction property, so a
moisture-dependent heat of sorption contributes its composition derivative to the
reaction-source block of the Jacobian.

> **Implementation note.** The call methods `(p)(T)` and `(p)(T, ξ)` are in
> `src/Materials/property_functions.jl`; the per-cell reaction-source Jacobian that
> differentiates the heats is `ReactionSourceOperator`
> (`src/Jacobian/operators/reaction_source.jl`).

### 4.5.3 The reaction-rate connection

The Arrhenius derivative (Eq. 4.7) is identical in form to the
temperature-derivative of the reaction-rate constant
$k(T)=A\exp(-E/(R_g T))$, since $\mathrm{d}k/\mathrm{d}T = k\,E/(R_g T^2)$; the
reaction-source Jacobian obtains it by ForwardDiff over the per-cell kernel
(Chapter 15). The full reaction-rate derivative including the depletion limiter and
temperature gates is developed in Chapter 6.

---

## 4.6 Material assembly

### 4.6.1 The `Material` type

A `Material` is the static, fully type-specialized aggregate that every physics
kernel reads. It is **immutable** and **heterogeneously typed** so the compiler
can specialize every loop over its components and reactions:

```julia
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ} <: AbstractMaterial
    name::Symbol
    components::C          # heterogeneous component tuple
    reactions::RT          # heterogeneous reaction tuple
    mixing::NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}
    lateral_shrinkage_law::LSL   # A(t)/A₀ = law(χ̄), or nothing
    depletion_limiter::DepletionLimiter   # tanh(ξ/threshold) rate roll-off (§6.3)
    gas_viscosity::Pμ      # pore-gas μ(T), default SutherlandProperty() (§5.9)
end
```

The type parameters carry compile-time information:

- `NC` — number of components, a static integer driving every `@generated`
  unroll (`Val(NC)`).
- `NR` — number of reactions, likewise static.
- `C`, `RT` — the *exact* concrete tuple types of the components and reactions.
  Because these are full tuple types (not `Vector{AbstractComponent}`), iterating
  over them is type-stable and allocation-free; the cost is that each distinct
  material *shape* is a distinct compiled type.
- `LSL` — the lateral-shrinkage-law type, a concrete `Function` subtype or
  `Nothing`. Parameterizing it here lets calls through
  `material.lateral_shrinkage_law` specialize statically rather than dispatch
  through a `Union{Function,Nothing}` (Chapter 10). The inner constructor
  enforces `LSL === Nothing || LSL <: Function`.
- `Pμ` — the type of the pore-gas viscosity property function, parameterized for
  the same reason so `material.gas_viscosity(T)` stays monomorphic. The inner
  constructor probes the curve at 200–1800 K and rejects one that is not positive
  and finite there (§5.9).

The `mixing` field bundles three `MixingSpec`s — one each for thermal
conductivity ($k$), gas transfer ($\lambda$), and permeability ($\kappa$) — in a
fixed-order `NamedTuple` so that `material.mixing.k` compiles to a direct field
read. Each `MixingSpec` carries `(rule, β, ℓ)`: the mixing rule, the WEIGHTED
blend weight $\beta$ (default 0.5), and the Carman–Kozeny characteristic length
$\ell$ (NaN unless the rule is `CARMAN_KOZENY`). The default channel rules are

$$
k:\ \text{WEIGHTED},\ \beta=0.5;\qquad
\lambda:\ \text{PARALLEL};\qquad
\kappa:\ \text{WEIGHTED},\ \beta=0.5,
$$

built by `MaterialMixing()`. The `MixingSpec` constructor validates
$0\le\beta\le1$ and, for `CARMAN_KOZENY`, requires a finite positive
`length_scale`, failing fast otherwise. The semantics of each rule (PARALLEL,
SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY), the BRUGGEMAN rejection for $\lambda$
and $\kappa$, and the Newton solve for BRUGGEMAN, are all developed in Chapter 5.

The `depletion_limiter` field carries the material's own `DepletionLimiter`
(default `threshold = 1.0` kg·m⁻³, `enabled = true`), set via the constructor's
`depletion_limiter` keyword and read by `evaluate_reactions!` for every
reaction-rate evaluation (§6.3).

> **Implementation note.** `Material`, its inner/outer constructors,
> `_construct_material`, `MixingRule`, `MixingSpec`, `MaterialMixing`, in
> `src/Materials/materials.jl`. Counts via `n_components`, `n_reactions`.

### 4.6.2 `PendingReaction` and symbolic resolution

Reactions may be written against component **symbols** rather than integer
indices, which is far less error-prone for the user:

```julia
Reaction(:decomposition, :virgin => :MMA;          A=8.5e12, E=188e3, h=870e3)
Reaction(:charring,      :virgin => (:char, :gas); A=1e15,  E=2e5, h=5e5, yields=(0.2,0.8))
```

A symbolic constructor produces a `PendingReaction`, which stores the symbol
payloads plus the kinetic parameters. When `Material(...)` is called, it builds a
`name → index` map from the components tuple (erroring on a duplicate component
name) and resolves every `PendingReaction` into a numeric `Reaction`. The
resolution path is taken **only if** at least one reaction is pending — the common
all-numeric case skips the `Dict` allocation entirely:

```
needs_resolution = any(rxn isa PendingReaction for rxn in reactions)
needs_resolution || return reactions          # fast path, no Dict
```

Numeric reactions (built against integer indices directly) bypass resolution.

> **Implementation note.** `PendingReaction`, the symbolic `Reaction` constructors,
> `resolve_reaction`, and `_resolve_reactions` (called from the `Material`
> constructor) in `src/Materials/reactions.jl` and `src/Materials/materials.jl`.

### 4.6.3 The single-reactant mass-balance convention

The `Reaction` type is n-ary in both reactants and products in principle, but the
current model **supports only single-reactant reactions** ($N_\mathrm{in}=1$). The
mass-balance convention (developed fully in Chapter 6) is:

- The lone reactant is consumed with implicit mass stoichiometry $-1$ per unit
  reaction extent.
- Each product $j$ is produced with an explicit mass yield $\nu_j > 0$.
- Mass conservation requires $\sum_{\text{products}} \nu_j = 1$.

The `Reaction` inner constructor enforces this at build time — every
construction path, including the bare positional constructor, runs the
validators. Reaction orders must be $n \ge 0$ (`_validate_reaction_order`);
every product yield must be $\ge 0$ (`_validate_yield`, enforced
unconditionally — a negative yield consumes a species the rate law does not
see and drives its concentration negative without bound; note $(1.5, -0.5)$
sums to 1 and would pass a sum-only check); and for a single-reactant reaction
the product yields must sum to 1 (`_validate_stoichiometry`, relative tolerance
$\max(10^{-10}\cdot(1+\Sigma\nu), 10^{-14})$), raising an `error` on violation.
Passing `validate_mass=false` disables only the sum check (not recommended);
the non-negativity and order checks always run.

> **Implementation note.** `_validate_stoichiometry`, `verify_stoichiometry`
> (tol $10^{-10}$, returns `false` + `@warn` on failure), and `mass_balance_error`
> in `src/Materials/reactions.jl`.

### 4.6.4 Dry-pyrolysis-gas identification

For variable-cross-section (lateral shrinkage, Chapter 10) the pyrolysis-progress
state $\chi$ counts only **dry** pyrolysis gas — water vapor is excluded because
moisture evaporation should not drive radial shrinkage. The set of dry-gas
component indices is computed at compile time by a recursive tuple walk that
includes a component index iff it is a gas **and** its name is not a water-vapor
name. The water-vapor test is a substring/equality heuristic on the component
name (`"H2O"`, `"water"`/`"Water"`, or exactly `"h2o"` case-insensitively):

```
is_gas(c) && !_is_water_vapor_name(c.name)  ⇒ include index
```

This name-based identification is a documented heuristic (open issue §4.9); a
gas whose name does not contain a recognized water token is treated as dry
pyrolyzate regardless of its actual chemistry.

> **Implementation note.** `dry_pyrolysis_gas_indices` →
> `_collect_dry_pyrolysis_gas_indices` and the heuristic `_is_water_vapor_name`,
> `src/Materials/materials.jl`. The result is a statically known tuple of indices,
> constant-folded by the compiler.

### 4.6.5 Type-stable iteration over heterogeneous tuples

The defining performance pattern of the Materials subsystem is the use of
`@generated` functions and recursive tuple decomposition to **unroll** every
loop over the (heterogeneous) component and reaction tuples at compile time. A
naive `for c in material.components` over a heterogeneous tuple would box each
element and allocate; instead, e.g., the mixture heat capacity is

```julia
@generated function _sum_heat_capacity(ξ, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        result = zero(_prop_type(ξ, T))
    end
    for i = 1:NC
        push!(expr.args, :(result += ξ[$i] * comps[$i].heat_capacity(T)))
    end
    push!(expr.args, :(return result))
    return expr
end
```

which emits $N_C$ statically dispatched terms with no runtime loop. The same idiom
appears in `evaluate_reactions!`, the species-source and gas-flux kernels, the
mixing-rule sums (Chapter 5), and the porosity sums. The net effect is that the
entire property-evaluation loop in `update_properties!` is **allocation-free**
once the cache is sized (Chapter 5, §5.13): all composition-dependent dispatch is
generated, all intermediate tuples are statically sized, and there are no
dictionary lookups in hot paths.

---

## 4.7 Time- and temperature-parameter evaluation

Several kernels accept "flexible" parameters that may be a constant, a function
of time, or a property function — for example a boundary heat flux $q(t)$, a
ramped convective coefficient, or a temperature-dependent property given as a
property object. Two thin, AD-preserving dispatch helpers normalize these.

**Time-dependent parameters** (`evaluate_param(param, t)`):

$$
\texttt{evaluate\_param}(p, t) =
\begin{cases}
\texttt{convert}(\mathrm{typeof}(t),\,p) & p::\text{Real},\\
p(t) & p::\text{Function}.
\end{cases}
$$

The `convert(typeof(t), p)` for a constant is the crucial detail: when $t$ is a
`ForwardDiff.Dual` (during a sensitivity solve, Chapter 17), a bare `Float64`
constant is promoted to the Dual type so the computation stays in the AD number
type end to end. Returning the raw `Float64` would break the Dual propagation.

**Temperature-dependent parameters** (`evaluate_param_T(param, T)`):

$$
\texttt{evaluate\_param\_T}(p, T) =
\begin{cases}
\texttt{convert}(\mathrm{typeof}(T),\,p) & p::\text{Real},\\
p(T) & p::\text{Function},\\
p(T) & p::\text{AbstractPropertyFunction}.
\end{cases}
$$

with the same Dual-preserving `convert` for constants. A dedicated fast-path
method short-circuits the most common case:

```julia
@inline evaluate_param_T(param::ConstantProperty, T) = param.value
```

This bypasses the `param(T)` dynamic-dispatch hop and reads the value field
directly — measurable in tight loops where the same constant property is queried
per cell per step.

> **Implementation note.** `evaluate_param`, `evaluate_param_T`, and the
> `ConstantProperty` fast path in `src/Materials/property_functions.jl`.

---

## 4.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

Pyrolysis.jl's material model is a deliberate synthesis of the three reference
condensed-phase pyrolysis codes. We map their notation to ours on first use.

### 4.8.1 Component / species representation

**ThermaKin2Ds** represents a material as a set of *components* with mass
concentrations as the primary state — the convention Pyrolysis.jl adopts wholesale
($\xi_j$ in kg·m⁻³ is the species state, not mass fraction). ThermaKin distinguishes
condensed from gaseous components and applies the *solid-volume-fraction*
convention to conductivity and emissivity — exactly Pyrolysis.jl's
`solid_volume_fractions` (Chapter 5, §5.1), which excludes non-swelling gases
(ThermaKin: "when $\gamma=0$, gaseous components do not contribute to the thermal
conductivity"). ThermaKin's stoichiometric coefficient $\theta_i^{\,j}$ maps to our
$\nu_{i,j}$ (overload note A2).

**Gpyro** uses a control-volume formulation with $\bar\rho\,\Delta z$ (mass per
unit area) and mass fractions $Y_i$ of *condensed-phase species*; its gas species
are tracked separately. Gpyro's volume-fraction property averaging $X_\alpha$
maps to our $v_j$, and its power-law temperature dependence of properties is
representable here as a `PolynomialProperty` or `CallableProperty`. Gpyro's
solid-fraction $\mathrm{SF}_k$ / conversion $\chi_k$ mechanism for selecting
non-charring / charring / intumescent / shrinkage behavior corresponds to our
swelling-factor $\gamma_j$ plus the lateral-shrinkage law (Chapter 10).

**FDS** (solid phase) uses a layered material model with per-material density,
$c_p$, $k$, and emissivity, plus Arrhenius reaction parameters; its
pre-exponential $Z$ maps to our $A_i$ (overload note A3) and its solid absorption
coefficient $\kappa_s$ maps to our volumetric $\alpha$ (overload note A5). FDS's
material density-ratio swell/shrink rule corresponds to our $\gamma_j$-weighted
volume change (Chapter 10).

### 4.8.2 Heat capacity and the gas-storage convention

ThermaKin's energy equation (its Eq. 17) sums the volumetric heat capacity
over **all** components, gas included:
$\rho c_p^{\,\text{ThermaKin}} = \sum_{j\in\text{all}} \xi_j c_{p,j}$. Pyrolysis.jl
uses the same sum in the temperature-storage term (Chapter 3):

$$
\rho c_p^\mathrm{eff} = \sum_{j=1}^{N_C} \xi_j\, c_{p,j}(T),
\tag{4.9}
$$

because the pore gas shares the matrix temperature and its concentration is an
integrated state, so its sensible storage is available at no cost (Chapter 5,
§5.4 and Chapter 7). Gpyro's default thermal-equilibrium form drops the gas
term (restored by its `FULL_QSG` option) and FDS has no in-depth gas; both are
approximations of (4.9) that Pyrolysis.jl does not make. The gas *flux*
relative to the skeleton enters the energy equation only through the advective
source $S_\mathrm{conv}$.

> **Implementation note.** `effective_heat_capacity` → `_sum_heat_capacity`
> (every component, no phase branch), `src/Properties/effective_properties.jl`.

### 4.8.3 Reaction kinetics representation

ThermaKin supports bimolecular rates $\propto \xi_{\text{COMP1}}\xi_{\text{COMP2}}$;
Gpyro uses an nth-order-in-conversion rate with optional O₂ sensitivity
$g(Y_{O_2})$ and a critical-temperature Heaviside; FDS uses
$\rho^n\cdot Y_{O_2}^{n}\cdot$ a temperature factor with kinetic-parameter recovery
from a TGA reference temperature and peak rate. Pyrolysis.jl restricts to
**single-reactant** mass-action kinetics in mass-concentration basis with smooth
tanh temperature gates and a smooth depletion limiter (Chapter 6). The
representation of heats of reaction as property functions — including the
state-dependent moisture-sorption heat via `StateDependentProperty` — generalizes
all three codes' scalar heats. The full kinetic comparison is in Chapter 6.

### 4.8.4 Property functions and temperature dependence

ThermaKin and Gpyro both permit temperature-dependent properties (tabulated or
power-law). Pyrolysis.jl's property-function abstraction subsumes both: a Gpyro
power law $p(T)=p_0(T/T_0)^m$ is a `CallableProperty`; a ThermaKin table is a
`TableProperty`; FDS ramp-style piecewise data is a `TableProperty` with linear
interpolation and flat extrapolation. The novel addition is the uniform
$(T,\xi)$ call signature, which allows composition-dependent properties (heat of
sorption) without a special-case code path.

---

## 4.9 Limitations and open issues

The following are documented limitations of the constitutive layer; several are
revisited where the relevant physics is developed.

1. **Single-reactant reactions only.** $N_\mathrm{in}=1$ is enforced. Multi-step
   schemes must be modeled as *sequences* of single-reactant reactions.
   Multi-reactant (e.g. genuine bimolecular char oxidation) would require
   extending the reactant tuple to 3-ary `(idx, order, mass_stoich)` and a
   revised mass-balance convention (Chapter 6).

2. **Intrinsic-density fallback may be unphysical.** When $\rho_{i,j}$ is not
   supplied it defaults to the bulk $\rho_j$, asserting the pure phase has no
   internal pores. For genuinely microporous solids this under-predicts porosity
   and over-predicts the matrix volume.

3. **Name-based dry-gas identification.** Water vapor is excluded from the
   pyrolysis-progress source by a substring match on the component name
   (`_is_water_vapor_name`). A water species named outside the recognized tokens
   would be (incorrectly) counted as dry pyrolyzate, and a non-water species named
   with a water token would be (incorrectly) excluded.

4. **Kinked and AD-dependent property laws.** A `TableProperty` is differentiated
   segment by segment, so its derivative jumps at every knot and the Newton
   iteration sees a kink there. `CallableProperty` and `StateDependentProperty`
   are differentiated by propagating `ForwardDiff.Dual` numbers through the user's
   function, which must therefore be generic in its argument types (§4.5).

5. **Heterogeneous-tuple compilation cost.** Preserving the exact component and
   reaction tuple types yields allocation-free, type-stable hot paths, but at the
   cost of recompiling the residual for each distinct material *shape*. Materials
   that differ only in numerical parameters share a compiled type; materials that
   differ in the *types* of their property functions do not.

6. **No temperature-dependent pre-exponential or pressure-dependent kinetics.**
   The Arrhenius property and `Reaction` carry constant $A$ and $E$ (no Kissinger/
   Ozawa $T$-dependent activation energy, no distributed-activation-energy
   tooling, no explicit pressure dependence). See Chapter 6.

7. **Heat-of-reaction sign mapping to reference codes.** Internally $h>0$ is
   endothermic ($Q_\mathrm{rxn}=-\sum_r h_r r_r$), whereas ThermaKin/Gpyro publish
   $h>0$ as exothermic. Users porting parameters from those codes must flip the
   sign (Chapter 2, overload note H1; Chapter 6).


---

# 5. Effective (Mixture) Properties and Mixing Rules

This chapter derives, in full, how Pyrolysis.jl collapses the per-component
thermophysical and transport properties of a multi-phase pyrolyzing solid into
the *single* effective property each governing equation needs. It covers the two
averaging bases (mass fractions and volume fractions, and the two flavours of the
latter), the swelling-aware mixture density, the intrinsic-density porosity, the
mixture effective heat capacity, the four
conductivity mixing rules (parallel / series / weighted / Bruggeman) together
with the Bruggeman Newton solve, the effective gas-transfer and permeability
closures (component mixing versus Carman–Kozeny), the radiative properties
(emissivity and absorption), the gas viscosity (Sutherland) and Chapman–Enskog
diffusivity, the absence of per-species transport scaling, the distance-weighted
harmonic face averaging, the property cache, and the mixing derivatives the
Jacobian consumes. Every equation, constant, guard, threshold, and sign below has
been confirmed line-by-line against `src/Properties/mixture_rules.jl`,
`src/Properties/effective_properties.jl`, `src/Properties/property_cache.jl`,
`src/Materials/components.jl`, and `src/Materials/materials.jl`; where the
implementation departs from the textbook form (most importantly which basis each
property uses, the harmonic-mean face average, and which mixing rules are
rejected for which channel), **the code is authoritative** and the deviation is
called out explicitly.

Throughout this chapter `j` is the **component (species) index** running over the
`N_C` components of a `Material`, and `i` is the **cell index** running over the
`n` finite volumes (overload resolution I1, §2). The primary species state is the
**mass concentration** `ξ_j` (kg·m⁻³), not the mass fraction (§3, Governing
Equations). Two density symbols are used with care: the **bulk** (pure-phase)
density `ρ_j` drives all mixing rules and the mixture density, while the
**intrinsic** (skeletal, cell-wall) density `ρ_{i,j}` appears *only* in the
porosity formula (overload resolution A7, §2). The three binding sign conventions
of §2 are assumed without restatement.

---

## 5.1 Averaging bases: mass fractions and volume fractions

A homogenized property must be a weighted combination of its component values,
and the *correct weight* depends on the physics. Pyrolysis.jl uses three distinct
weighting bases, each computed from the cell's mass-concentration tuple
`ξ = (ξ_1,…,ξ_{N_C})` and (for the volume bases) the temperature `T`.

### 5.1.1 Total density and mass fractions

The **total (bulk) density** of a cell is simply the sum of mass concentrations,

$$
\rho_{\text{tot}}(\boldsymbol\xi) \;=\; \sum_{j=1}^{N_C} \xi_j ,
\tag{5.1}
$$

This is the density in *bulk-density space* — the mass per unit bulk volume — before any
swelling correction (§5.2). It is the natural denominator for **per-unit-mass**
properties, because such a property averaged over the mixture must weight each
component by its share of the mass.

The **mass fraction** of component `j` is

$$
Y_j \;=\; \frac{\xi_j}{\rho_{\text{tot}}}, \qquad \sum_j Y_j = 1 ,
\tag{5.2}
$$

Mass fractions are the conceptual weights for the swelling mixture density
(§5.2), where `Y_j` is folded in inline (with an empty-cell guard
`ρ_tot < 1e-20`) so that no mass-fraction tuple is built; the heat-capacity
storage term weights by the raw concentrations `ξ_j` directly (§5.4). They are
*not* the right weights for transport properties, which act over volumes.

### 5.1.2 All-phase volume fractions

A transport property such as thermal conductivity is governed by *how much
volume* each phase occupies along the heat-flow path, not by how much mass it
carries. The **bulk volume fraction** of component `j` is

$$
v_j \;=\; \frac{\xi_j/\rho_j(T)}{\displaystyle\sum_{k=1}^{N_C}\xi_k/\rho_k(T)},
\qquad \sum_j v_j = 1 ,
\tag{5.3}
$$

where `ρ_j(T)` is the component's **bulk** density evaluated through its property
function (constant, polynomial, table, or callable — §4). The quantity
`ξ_j/ρ_j` is the volume occupied by component `j` per unit bulk volume, so (5.3)
normalizes those volumes to a partition of unity. This is `volume_fractions(ξ, T,
material)` (`mixture_rules.jl`).

The function builds the denominator `v_total = Σ_k ξ_k/ρ_k` via
`_sum_specific_vol` and, if `v_total < 1e-20` (a vanished cell), returns the zero
tuple. Each per-component term is itself guarded `ξ_k > 0 && ρ_k > 0`, so a
component with zero concentration or a (non-physical) zero density contributes
nothing. The all-phase basis is used for **gas transfer** (§5.6), where every
phase — including the gas that occupies pore space — participates in the
volumetric averaging; component-based permeability (§5.7.1) uses the
swelling-weighted basis (5.4).

### 5.1.3 Swelling-weighted (matrix) volume fractions and why surface properties use them

Thermal conductivity and surface emissivity are determined by the *condensed
matrix*, not by gas transiently occupying the pores. Pyrolysis.jl therefore
defines a **swelling-weighted volume fraction** in which every component
contributes the γ-fraction of its intrinsic volume:

$$
v_j^{\text{solid}} \;=\;
\dfrac{\gamma_j\,\xi_j/\rho_j(T)}{\displaystyle\sum_{k}\gamma_k\,\xi_k/\rho_k(T)},
\tag{5.4}
$$

so the fractions sum to one over the contributing (γ > 0) components. This is
`solid_volume_fractions(ξ, T, material)` (`mixture_rules.jl`). The swelling
factor here has **one universal meaning** — the fraction of a component's
intrinsic volume contributing to the bulk volume — shared with the dilation
rate $\theta = \sum_j\gamma_j(\mathrm{d}\xi_j/\mathrm{d}t)/\rho_j$ (§3.6,
§10) and the swelling mixture density (5.6).

The physical justification is the **ThermaKin convention**: "when γ = 0,
gaseous components do not contribute to the thermal conductivity"
(ThermaKin2Ds Reference, Eqs. 8–10), reproduced verbatim in the source
docstring. A non-swelling gas (`γ = 0`, the default for `GaseousComponent`,
§4) fills *existing* pore space and escapes; it neither forms a continuous
solid conduction path nor presents a radiating surface, so it carries zero
weight. A gas given `γ > 0` (an intumescent foaming agent that *creates* new
volume) is retained with weight γ.

> **Subtlety (γ = 0 condensed phases).** Because the weighting is universal,
> a *condensed* component with `γ = 0` (the `LiquidComponent` default, and
> the "rigid char" pattern) is consistently **volumeless**: zero weight here,
> zero strain in θ, and zero volume in the mixture density. Such a component
> contributes its mass and heat capacity but not its k or ε.

> **Void-cell blend.** Conductivity and emissivity consumers go through
> `mixing_volume_fractions`, which returns the swelling-weighted matrix
> fractions while the matrix volume fraction $\varphi_m = \sum_j \gamma_j\xi_j/\rho_j$
> is above $10^{-3}$, the all-phase occupancy (5.3) below $10^{-4}$, and a cubic
> smoothstep blend of the two in between. A fully converted cell thus degrades
> to the gas-mixture conductivity instead of 0 (which would zero both adjacent
> harmonic face means and, via Kirchhoff, its absorption), and because the
> matrix fractions are normalised by $\varphi_m$ their derivatives grow as
> $1/\varphi_m^2$ — the smooth blend keeps $\partial k_{\rm eff}/\partial\xi$
> bounded while a cell empties on a fixed mesh, which parks only its top live
> cell (§14.3.2).

The three bases and their consumers are summarized below; this mapping is the
single most error-prone part of effective-property bookkeeping, so it is stated
once authoritatively from the code in `update_properties!`
(`property_cache.jl`):

| Basis | Symbol | Built by | Consumed by |
|---|---|---|---|
| Mass fraction | `Y_j` | inline in `_specific_vol_with_swelling_direct` | swelling mixture density (§5.2; the heat-capacity storage term weights by raw `ξ_j`, §5.4) |
| All-phase volume fraction | `v_j` | `volume_fractions` | `effective_gas_transfer` |
| Swelling-weighted volume fraction | `v_j^solid` | `solid_volume_fractions` (via `mixing_volume_fractions`, which adds the depleted-cell all-phase fallback) | `effective_conductivity`, `effective_emissivity`, `effective_absorptivity`; `effective_permeability` reads `solid_volume_fractions` directly |

Note in particular that **gas transfer uses the all-phase basis `v_j`** while
**conductivity and emissivity use the solid basis `v_j^solid`** — this is what
`update_properties!` does (`cache.λ[i] = _gas_transfer_from_v(v_full, …)` versus
`cache.k[i] = _conductivity_from_v(v_solid, …)`).

---

## 5.2 Mixture density with swelling

The **mixture (bulk) density** `ρ` answers a different question than
`ρ_tot`: it is the mass per unit *bulk volume* once each phase is allowed to
contribute, or not contribute, to that bulk volume according to its swelling
factor `γ_j`. The constitutive relation is a swelling-weighted harmonic
(specific-volume) average with the **universal** γ weighting (every
component, not gases only — the same weighting as θ and (5.4)),

$$
\frac{1}{\rho}
\;=\; \sum_{j} \gamma_j\,\frac{Y_j}{\rho_j(T)},
\tag{5.6}
$$

with `Y_j = ξ_j/ρ_tot` the mass fraction. This is `mixture_density(ξ, T,
material)` (`mixture_rules.jl`).

The semantics of `γ_j` follow directly from (5.6):

- **`γ_j = 1` (default for solids, `SolidComponent`).** The
  component contributes its full specific volume. When such a component depletes,
  its specific volume disappears and `1/ρ` decreases, i.e. the bulk *shrinks*.
  This is the ordinary charring/recession behaviour.
- **`γ_j = 0` (default for gases, `GaseousComponent`, **and for liquids**,
  `LiquidComponent`).** The term drops out of (5.6) entirely. Gas generated by
  a reaction fills pre-existing pore space and escapes through the matrix
  without inflating the bulk; a γ = 0 liquid is moisture bound in cell walls
  (ThermaKin `L SWELLING: 0`) — consistently volumeless across (5.4), (5.6),
  and θ. Set `γ = 1` for a liquid that physically occupies bulk volume (a
  melt layer).
- **`γ_j ∈ (0,1)`.** Partial contribution — the model for an intumescent
  blowing agent that *does* create new bulk volume as it is produced. Setting a
  gas's `γ_j` near 1 makes the generated gas expand the solid.

The volumetric strain rate `θ` that converts (5.6) into actual cell-volume change
(and hence into mesh motion under ALE) is derived in §10 (Volume Change and
Lateral Shrinkage) and §11 (The ALE Framework); here we only fix the
constitutive density.

> **Implementation note.** `mixture_density` computes the mass fractions
> *inline* via `_specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps, Val(NC))`
> (`mixture_rules.jl`), folding `Y_j = ξ_j · (1/ρ_tot)` into the
> specific-volume accumulation so that no intermediate `Y` tuple is allocated.
> Two empty-cell guards bracket the computation: `ρ_tot < 1e-20 ⇒ 0`, and
> `v_specific < 1e-20 ⇒ 0` (`mixture_rules.jl`). The swelling factor is
> applied per component (one universal γ semantics, §5.1.3), so a mixture may
> contain swelling and non-swelling components of any phase simultaneously.

---

## 5.3 Porosity

### 5.3.1 The canonical intrinsic-density form

**Porosity** `φ` is the bulk-volume fraction *not* occupied by the solid/liquid
skeleton. It is the only place the **intrinsic (skeletal) density** `ρ_{i,j}` is
used. The canonical definition is

$$
\phi \;=\; \operatorname{clamp}\!\left(
1 - \sum_{j\in\{\text{solid,liquid}\}} \frac{\xi_j}{\rho_{i,j}(T)},
\;0,\;1\right),
\tag{5.7}
$$

implemented as `porosity(ξ, T, material)` →
`_porosity_void_frac` (`mixture_rules.jl`). The sum runs over **solid and
liquid** components only — the generated code branches on
`is_solid(c) || is_liquid(c)`, and each term is taken only if
`ρ_intr > 0`. Gases are *excluded from the sum* because they occupy the pore
space being measured, not the skeleton. The result is clamped to `[0,1]`; a
clamp to 0 occurs only in the (numerical-error) case where the computed solid
volume exceeds the bulk volume.

The distinction between bulk and intrinsic density is essential here. The
*intrinsic* density is the density of the cell-wall material itself, ignoring any
micro-pores inside the pure phase; the *bulk* density is the apparent density of
that phase as it sits in the composite. For a porous solid these differ, and
using the intrinsic density in (5.7) correctly recovers the void fraction. When a
component is constructed **without** an explicit intrinsic density,
`get_intrinsic_density` falls back to the bulk density (`components.jl`),
which makes the summand `ξ_j/ρ_j`: the pure phase is taken to have no internal
pores. For a gas, `get_intrinsic_density(::GaseousComponent,
::Real) = 0.0`, consistent with its exclusion from the sum.

Porosity feeds two downstream models:

1. the **ideal-gas pressure closure** `P = N(ξ,T)/φ`, with
   `N = Σ_{g} ξ_g R_g T / M_g` the partial-pressure numerator (§9, Gas-Phase
   Mass Transport and Pressure); and
2. the **Carman–Kozeny permeability** closure (§5.7).

### 5.3.2 Porosity derivatives

Differentiating (5.7) term by term, and recognizing that the clamp is inactive in
the physical interior `0 < φ < 1` (where derivatives are needed), gives

$$
\frac{\partial \phi}{\partial \xi_j}
=
\begin{cases}
-\dfrac{1}{\rho_{i,j}(T)}, & j\in\{\text{solid, liquid}\}\ \text{and}\ \rho_{i,j}>0,\\[8pt]
0, & j\in\text{gas (or }\rho_{i,j}\le 0).
\end{cases}
\tag{5.8}
$$

The sign is negative for the skeleton: adding solid mass reduces the void
fraction. Gases contribute zero — they fill pore space already accounted for by
the solid complement, so changing a gas concentration does not change `φ`.
These derivatives are not hand-coded: the Jacobian path that needs them
(the `DarcyTransportOperator`'s Carman–Kozeny permeability and pressure
chains, §5.15, §15.8.3) differentiates
`porosity` itself with forward-mode `Dual`s, which is clamp-consistent by
construction — at the `clamp(·, 0, 1)` bounds the propagated derivative is 0,
exactly matching the residual's clamped value.

---

## 5.4 Effective volumetric heat capacity

### 5.4.1 Mixture heat capacity

The temperature-storage term of the energy equation (§3, §7) uses the
**mixture** volumetric heat capacity,

$$
(\rho c_p)_{\text{eff}}
\;=\; \sum_{j=1}^{N_C} \xi_j\, c_{p,j}(T) ,
\tag{5.9}
$$

implemented as `effective_heat_capacity(ξ, T, material)` →
`_sum_heat_capacity` (`effective_properties.jl`). The generated helper adds
`ξ[j] * comps[j].heat_capacity(T)` for every component with no phase branch,
so **gas sensible storage is part of the left-hand side of the energy
equation**. The units check: `[kg·m⁻³] · [J·kg⁻¹·K⁻¹] = [J·m⁻³·K⁻¹]`.

The result is floored at `RHO_CP_FLOOR = 10³` J·m⁻³·K⁻¹
(`src/Materials/constants.jl`), the thermal mass of pore gas at ambient
(ρ_g c_{p,g} ≈ 1.2 × 1000; Gpyro floors the same quantity). Because every
component is summed, ρc_p reaches the floor only in a cell emptied of all
mass — a burned-out cell whose pore gas has also left; without it the 1/ρc_p
factor in dT/dt would blow up there. Virgin solids carry ρc_p ~ 10⁵–10⁶
J·m⁻³·K⁻¹, so the floor never binds in a cell that holds material. The floor
is applied *inside* `effective_heat_capacity`, so every consumer — the
residual, the property cache and the
structured-Jacobian operators that extract ∂ρc_p/∂(T, ξ) through Duals —
sees the identical clamped value and the consistent (zero) derivative at the
floor.

### 5.4.2 Why the gas is stored

The one-temperature model puts the pore gas at the matrix temperature (§3.1.2),
so the gas stores sensible heat `ξ_g c_{p,g} ∂T/∂t` exactly like every other
component. Its concentration `ξ_g` is an integrated state, so the term costs
nothing; dropping it — Gpyro's default *quasi-steady-gas approximation*, which
its `FULL_QSG` option reverses, and the implicit convention of FDS, which has
no in-depth gas — would be an approximation with no benefit. The in-pore gas
heat capacity is roughly `φ ρ_g c_{p,g}` against the matrix `(1−φ) ρ_s c_{p,s}`;
at 1 atm and pyrolysis temperatures the ratio is a few tenths of a percent and
it scales with pressure. On white-pine cone decks the term is at most about
half a percent of ρc_p and shifts the response by a uniform time lag of about
a tenth of a percent — small, but a thousand times above the solver's noise
floor and therefore worth keeping exact.

The gas *flux* relative to the skeleton enters the energy equation only
through the advective source `S_conv = −Σ_g c_{p,g}(T) J̄_g (∂T/∂z)` (§7). This is
the ThermaKin Eq. 17 convention term for term (§3.3.2).

No separate mass-averaged specific heat is computed or cached; consumers that
need `c_p` per unit mass divide `(ρc_p)_eff` by `ρ_tot`.

---

## 5.5 Thermal conductivity mixing

The effective thermal conductivity is computed from the swelling-weighted
volume fractions `v_j^solid` (5.4) and the per-component conductivities
`k_j(T)` under one of four mixing rules selected by `material.mixing.k.rule`:
`effective_conductivity(ξ, T, material)` builds
`v = mixing_volume_fractions(…)` (the (5.4) weights, with the depleted-cell
all-phase fallback of §5.1.3) then dispatches to `_conductivity_from_v`
(`effective_properties.jl`). The default rule is `WEIGHTED` with
`β = 0.5` (`materials.jl`).

### 5.5.1 Parallel and series bounds

The two classical bounds correspond to the geometric extremes of phase
arrangement relative to the heat-flow direction. With

$$
k_\parallel = \sum_j v_j^{\text{solid}}\, k_j(T),
\qquad
\frac{1}{k_{\text{series}}} = \sum_j \frac{v_j^{\text{solid}}}{k_j(T)},
\tag{5.13}
$$

the **PARALLEL** rule returns `k = k_∥` and the **SERIES** rule returns
`k = k_series`. Physically:

- **Parallel (upper bound, Voigt/Wiener):** all phases span the slab in the
  heat-flow direction and conduct independently in parallel; thermal
  conductances add. Appropriate when every phase forms a continuous path along
  the gradient.
- **Series (lower bound, Reuss/Wiener):** phases are stacked perpendicular to
  the gradient so heat passes through each in turn; thermal *resistances* add,
  giving a harmonic mean. Appropriate for stratified layering.

The true effective conductivity of any two-phase mixture lies between these
Hashin–Shtrikman-bracketed bounds.

> **Implementation note.** Both sums are produced in one pass by the generated
> `_conductivity_sums(v, T, comps, Val(NC))`
> (`effective_properties.jl`). The series accumulation is guarded:
> `inv_k_series += v[i]/k` is taken only when `v[i] > 1e-20 && k > 1e-20`, so a
> phase with zero volume or (pathologically) zero conductivity does not blow up
> the harmonic sum. After the loop
> `k_series = inv_k_series > 1e-20 ? 1/inv_k_series : 0`
> (`effective_properties.jl`).

### 5.5.2 Weighted blend

The **WEIGHTED** rule linearly interpolates between the bounds with a single
user-supplied parameter `β ∈ [0,1]`:

$$
k_{\text{eff}} \;=\; \beta\, k_\parallel + (1-\beta)\, k_{\text{series}} .
\tag{5.14}
$$

`β = 1` recovers the parallel upper bound, `β = 0` the series lower bound, and the
default `β = 0.5` (`MixingSpec` default, `materials.jl`) is the arithmetic mean
of the two bounds. The blend is the package's default conductivity model
(`materials.jl`) because, lacking microstructural information, the mean of the
rigorous bounds is the least-biased single estimate. The `β` value is validated to
lie in `[0,1]` at `MixingSpec` construction (`materials.jl`).

### 5.5.3 Bruggeman effective-medium theory

For random heterogeneous microstructures — and especially for the self-similar
pore structure of forming char — the **BRUGGEMAN** symmetric effective-medium
closure is more faithful than a fixed blend. It treats each phase as a spherical
inclusion embedded in the *effective* medium and demands self-consistency: the
volume-averaged polarization vanishes,

$$
\sum_j v_j^{\text{solid}}\,
\frac{k_j - k_{\text{eff}}}{k_j + 2 k_{\text{eff}}} \;=\; 0 .
\tag{5.15}
$$

Unlike the parallel/series/weighted forms, (5.15) is *implicit* in `k_eff` and
can capture percolation-like behaviour (the effective conductivity dropping
sharply when a continuous phase fractures) and genuinely non-linear composition
dependence. The factor `+2 k_eff` is the spherical-inclusion depolarization
factor for an isotropic 3-D medium; the symmetric form (no distinguished host
phase) reflects an isotropic, statistically self-similar microstructure.

#### Newton–Raphson solve

Equation (5.15) is solved by Newton's method in
`bruggeman_effective_property(v, k; tol, maxiter, k_init)`
(`effective_properties.jl`). Writing the residual
`f(k_{\text{eff}}) = \sum_j v_j (k_j - k_{\text{eff}})/(k_j + 2 k_{\text{eff}})`,
the exact derivative is

$$
\frac{df}{dk_{\text{eff}}}
= \sum_j v_j\,\frac{-3\,k_j}{(k_j + 2 k_{\text{eff}})^2} ,
\tag{5.16}
$$

derived in the source comment line-by-line (`effective_properties.jl`):
the numerator of `d/dk_eff [v(k−k_eff)/(k+2k_eff)]` collapses
`−(k+2k_eff) − 2(k−k_eff)` to `−3k`. The Newton update is
`k_eff ← k_eff − f/df`, with the following exact guards and conventions taken from
the implementation:

- **Degenerate volume.** If `Σ v_j < 1e-20`, return 0
  (`effective_properties.jl`).
- **Single active phase.** Active phases are those with `v_i > 1e-20`; if exactly
  one is active, return its (volume-weighted) parallel mean `k_parallel`
  immediately (`effective_properties.jl`). This is exact, since a
  single-phase Bruggeman medium is that phase.
- **Warm start.** The Newton seed is `k_init` when a positive `k_init` is supplied
  (`k_init > 1e-20`), otherwise the parallel mean `k_parallel`
  (`effective_properties.jl`). In production `k_init` is the *previous time
  step's converged* `k_eff` for the cell (passed by `update_properties!` as
  `cache.k[i]`, `property_cache.jl`), so on smooth time-stepping the solve
  needs only 1–3 iterations.
- **Per-term contributions.** A phase enters the residual when `v_i > 1e-20`.
  A conducting phase (`k_i > 1e-20`) contributes the standard term; an
  **insulating** phase (`k_i ≤ 1e-20`) contributes its exact k → 0 limit
  `−v_i/2` to `f` (zero `df` contribution). The insulating term is what
  produces the percolation threshold: a conducting phase dispersed in an
  insulator yields `k_eff = 0` until its volume fraction exceeds 1/3
  (e.g. `v = (0.5, 0.5)`, `k = (1, 0)` → `k_eff = 0.25`).
- **Below-threshold collapse.** Below the percolation threshold the root is
  `k_eff = 0` and `f < 0` for every `k_eff > 0`, so the step-capped Newton
  decays `k_eff` geometrically; once `k_eff < 1e-30` the solve returns 0
  immediately instead of burning `maxiter`.
- **Singular-derivative bail-out.** If `|df| < 1e-20` the iteration breaks
  (it cannot take a Newton step; the problem is locally ill-posed)
  (`effective_properties.jl`).
- **Step limiter (positivity).** The raw step `dk = −f/df` is clipped so the
  update never reduces `k_eff` by more than 90 %: if
  `k_eff + dk < 0.1·k_eff`, set `dk = −0.9·k_eff`
  (`effective_properties.jl`). This prevents Newton from overshooting
  into a negative conductivity.
- **Convergence and return.** Iterate up to `maxiter = 100`, stopping when
  `|f| < tol = 1e-10`; the return is `max(k_eff, 0.0)`
  (`effective_properties.jl`).

> **Implementation note.** The component conductivities are gathered into a tuple
> by the generated `_conductivity_tuple(T, comps, Val(NC))`
> (`effective_properties.jl`), so the Bruggeman path is allocation-free.
> The dispatch in `_conductivity_from_v` checks `mixing == BRUGGEMAN` first and
> only otherwise computes the parallel/series sums
> (`effective_properties.jl`).

---

## 5.6 Effective gas-transfer coefficient

The diffusive gas flux (§9) uses an effective Fickian transfer coefficient
`λ_eff` [m²·s⁻¹], homogenized from the per-component `λ_j(T)`. Crucially — and
unlike conductivity — gas transfer is averaged over the **all-phase** volume
fractions `v_j` (5.3), *not* the solid basis:

$$
\lambda_{\text{eff}}(T,\boldsymbol\xi)
\;=\; \text{MixingRule}\bigl(v_j,\ \lambda_j(T)\bigr),
\tag{5.17}
$$

`effective_gas_transfer(ξ, T, material)` builds `v = volume_fractions(…)` and
dispatches to `_gas_transfer_from_v` (`effective_properties.jl`). This is
confirmed in `update_properties!` where the cached gas-transfer property is
`cache.λ[i] = _gas_transfer_from_v(v_full, T_i, material)`
(`property_cache.jl`). The supported rules mirror conductivity:

$$
\lambda_\parallel = \sum_j v_j\,\lambda_j,\qquad
\frac{1}{\lambda_{\text{series}}}=\sum_j\frac{v_j}{\lambda_j},\qquad
\lambda_{\text{eff}}=\beta\lambda_\parallel+(1-\beta)\lambda_{\text{series}}.
\tag{5.18}
$$

The default rule for `λ` is **PARALLEL** with `β = 0.5` (`materials.jl`).

**Degenerate transport.** `SolidComponent` defaults to $\lambda_j = 0$. A
condensed component with no gas-transfer coefficient contributes nothing to
(5.18), so a cell of that component holding no gas has $\lambda_{\text{eff}}
\approx v_{\rm gas}\lambda_{\rm gas} \to 0$ under PARALLEL (and exactly zero
under SERIES), and the face harmonic mean (5.29) vanishes when either side does.
Gas transport then stops where the gas fraction stops: the gas field of a
gas-producing material ends in a sharp front, and a Newton iterate that drives
the front cell's gas concentration below zero closes both of its faces and
flattens its residual row. The implicit solve fails repeatedly on that front.
On the Landau ablation case V-CPL-02 a solid coefficient of
$10^{-6}$ m²/s halves the step count and removes more than half of the
Newton failures, while the measured regression rate moves by $10^{-5}$
relative.
`validate_problem`, which `solve` calls, therefore warns when a material with a
gas-producing reaction and Fickian transport (no permeable component) has a
condensed component with $\lambda_j = 0$
(`warn_degenerate_gas_transport`, `src/Problem/problem_validation.jl`).

### 5.6.1 Bruggeman is rejected for gas transfer

Applying the effective-medium closure (5.15) to a diffusivity is **rejected**:
`_gas_transfer_from_v` throws an `ArgumentError` if `mixing == BRUGGEMAN`
(`effective_properties.jl`), with the message that the EMT closure "has
no physical justification for diffusivities." The reasoning is that a Fickian
transfer coefficient is already an average over a heterogeneous pore network
(molecular and Knudsen regimes), so wrapping it in a *second* effective-medium
self-consistency relation derived for a thermodynamic-like coefficient is not
physically supported. (The `effective_gas_transfer` docstring
(`effective_properties.jl`) explicitly states that `BRUGGEMAN` is **not**
supported for gas transfer and lists only `PARALLEL`/`SERIES`/`WEIGHTED`,
matching the running code, which is authoritative.)

> **Implementation note.** The series guard for `λ` is more elaborate than for
> `k`: in `_gas_transfer_sums` (`effective_properties.jl`), once any
> active phase has `λ ≤ 1e-20` the running `inv_λ_series` is set to `Inf` and the
> rest of the sum is short-circuited, which drives `λ_series → 0`. This is the
> correct series limit: a single impermeable-to-diffusion layer blocks the
> series path entirely.

---

## 5.7 Effective permeability

Permeability `κ` [m²] is the resistance-to-flow property in Darcy's law
`v_g = −(κ/μ)∇P` (§9). Only **solid and liquid** components carry a permeability;
gases flow *through* the matrix and `get_permeability(::GaseousComponent, T) =
0.0` (`components.jl`). Pyrolysis.jl offers two **mutually exclusive**
closures, selected by `material.mixing.κ.rule`
(`effective_permeability(ξ, T, material)`, `effective_properties.jl`).
The default rule for `κ` is **WEIGHTED** with `β = 0.5` (`materials.jl`).

### 5.7.1 Component mixing

When the rule is `PARALLEL`, `SERIES`, or `WEIGHTED`, the effective permeability
is a volume-fraction mixture of the component permeabilities, exactly as for
conductivity: on the **renormalized swelling-weighted basis** `v_j`
(`solid_volume_fractions` (5.4) — γ = 0 components carry no weight and the
surviving weights sum to 1):

$$
\kappa_\parallel = \sum_j v_j\,\kappa_j,\qquad
\frac{1}{\kappa_{\text{series}}}=\sum_j\frac{v_j}{\kappa_j},\qquad
\kappa_{\text{eff}}=\beta\kappa_\parallel+(1-\beta)\kappa_{\text{series}},
\tag{5.19}
$$

`_permeability_from_v` (`effective_properties.jl`), with the parallel/series
sums from `_permeability_sums`. An **impermeable phase** (`κ ≈ 0`, the
`SolidComponent` default) with `v_j > 0` makes the series (Wiener lower) bound
**exactly zero** — gas cannot cross a series-stacked impermeable layer; the
accumulator signals this with `1/κ_series = ∞`, the same convention as the
gas-transfer channel. The series guard floor is tighter than for
conductivity — `κ > 1e-30` rather than `1e-20` — because permeabilities
themselves are tiny (virgin wood ~10⁻¹⁵ m², char ~10⁻¹² m²). Representative
component values from the source docstring: virgin wood `~10⁻¹⁵ m²` (tight
cellular structure), char `~10⁻¹² m²` (open channels), impermeable `0.0 m²`.

### 5.7.2 Carman–Kozeny porosity closure

When the rule is `CARMAN_KOZENY`, permeability is *not* a per-component mix but a
porosity-based closure:

$$
\kappa \;=\; \min\!\left(\frac{\ell^2\,\phi^3}{180\,(1-\phi)^2},\;\ell^2\right),
\tag{5.20}
$$

`_carman_kozeny_permeability(φ, ℓ)` (`effective_properties.jl`), where
`ℓ` is the characteristic particle/pore diameter [m] supplied through the
`MixingSpec` `length_scale` and `φ` is the *canonical* porosity (5.7). The
constant **180** is the Kozeny constant for a packed bed of spheres. The
`φ³/(1−φ)²` dependence captures the steep growth of permeability as porosity
rises during virgin → char conversion: as `φ → 0` (fully solid), `κ → 0`
(impermeable). As `φ → 1` the packed-bed closure loses validity and diverges,
so `κ` is capped at the physical free-flow scale `ℓ²` (an open channel of width
`ℓ` has `κ ~ ℓ²/32`; `ℓ²` is a generous upper bound). The cap binds for
`φ ≳ 0.93`, where `180(1−φ)² < φ³`; a tiny `1e-30` denominator guard remains
solely to avoid `Inf`/`NaN` at exactly `φ = 1` before the `min` resolves
(`effective_properties.jl`).

Because Carman–Kozeny needs the porosity (which depends on `ξ` through the
intrinsic densities), it must be called through `effective_permeability(ξ, T,
material)` — calling the cache-aware `_permeability_from_v` with this rule throws
an `ArgumentError` explaining that concentrations are required
(`effective_properties.jl`). `MixingSpec(CARMAN_KOZENY; length_scale=ℓ)`
fails fast at construction if `ℓ` is not a positive, finite number
(`materials.jl`).

The closure derivative needed for the Jacobian (§5.15) is not hand-coded: the
`DarcyTransportOperator` (§15) differentiates `effective_permeability` itself by
forward-mode AD, so the Carman–Kozeny chain
(through the porosity's composition dependence) is picked up automatically. On
the capped branch the `min` gives the value-consistent zero `dκ/dφ` — value and
derivative saturate at the same point.

### 5.7.3 Bruggeman is rejected for permeability

As with gas transfer, `BRUGGEMAN` is rejected for permeability:
`_permeability_from_v` throws an `ArgumentError` if `mixing == BRUGGEMAN`
("EMT is only defined for scalar transport coefficients … use
PARALLEL/SERIES/WEIGHTED, or CARMAN_KOZENY", `effective_properties.jl`).
Permeability is a flow-resistance of the matrix geometry, not a polarizable
scalar property to which the symmetric EMT applies. Carman–Kozeny is the
physically motivated porosity-based alternative for the percolation-like
virgin → char transition.

### 5.7.4 Face-averaged permeability

For the Darcy velocity at a face, the resistance-in-series combination is
applied to the **mobility** κ/μ rather than to κ alone: the production path
(`update_darcy_properties!`, src/Discretization/finite_volume.jl) forms each
cell's mobility mobᵢ = κᵢ/μ(Tᵢ) and combines the two sides with the
distance-weighted harmonic mean (§5.12),

    mob_face = (d_L + d_R) / (d_L/mob_L + d_R/mob_R),

with mob_face = 0 exactly when either side's mobility is below `1e-30` — an
impermeable cell is an infinite series resistance, so gas cannot cross into it.
The Jacobian operator's `_dt_mobility_face` applies the identical rule.

---

## 5.8 Radiative properties: emissivity and absorption

### 5.8.1 Effective emissivity

The surface emissivity is a **solid**-volume-fraction average of the component
emissivities, each evaluated at the cell temperature, clamped to the physical
range:

$$
\varepsilon_{\text{eff}}(T,\xi) \;=\; \operatorname{clamp}\!\left(
\sum_j v_j^{\text{solid}}\,\varepsilon_j(T),\ 0,\ 1\right),
\tag{5.22}
$$

`effective_emissivity(ξ, T, material)` builds the solid basis and calls
`_emissivity_from_v(v, T, material)`, which evaluates `get_emissivity(comp, T)`
per component and applies `clamp(ε_eff, 0, 1)` (`effective_properties.jl`). The
solid basis is used for the same physical reason as conductivity (§5.1.3): the
radiating surface is the solid/liquid material, not trapped non-swelling gas.
Component emissivities are property functions (`ε` accepts the same forms as
`k`: number, linear tuple, callable, property function) with defaults of `0.9`
for `SolidComponent`, `0.95` for `LiquidComponent` (condensed liquids are strong
emitters — water ≈ 0.96), and `0.0` for `GaseousComponent`
(`src/Materials/components.jl`). Because the sum is evaluated generically over
the scalar type, a `ForwardDiff.Dual` temperature returns
$\partial\varepsilon_{\text{eff}}/\partial T$, which the surface-emission
(§12) and Beer–Lambert (§8) Jacobians pick up without operator-level
bookkeeping. The effective emissivity is the surface property used by the
re-radiation boundary condition (§12) and as the gray-body factor in the
volumetric P1 emission source of the experimental P1 model (§8), which `solve`
does not accept.

### 5.8.2 Effective surface absorptivity

The absorptivity for *incident* (heater-spectrum) radiation is mixed with the
same weights,

$$
a_{\text{eff}}(T,\xi) \;=\; \operatorname{clamp}\!\left(
\sum_j v_j^{\text{solid}}\,a_j,\ 0,\ 1\right),\qquad
a_j = \begin{cases}
\texttt{surface\_absorptivity}_j & \text{if set}\\
\varepsilon_j(T) & \text{otherwise (Kirchhoff)}
\end{cases}
\tag{5.22b}
$$

via `effective_absorptivity(ξ, T, material)` → `_absorptivity_from_v` →
`get_absorptivity(comp, T)`. The Kirchhoff branch resolves at compile time
(`Nothing` vs `Float64` is a type parameter of the component), and for a
material in which no component sets a surface absorptivity the two generated
sums perform the same floating-point operations in the same order, so
$a_{\text{eff}} \equiv \varepsilon_{\text{eff}}$ bit-for-bit — a material that sets no absorptivity
absorbs exactly as Kirchhoff's law predicts. A component with an explicit
$a_j$ decouples the absorbed heater flux from the temperature-dependent
emission: an emissivity that falls with temperature (integrating-sphere data for
virgin wood, for instance) does not drag the absorbed cone flux down with it.
$a_{\text{eff}}$ is consumed by a `RadiativeFluxBC` with
`absorptivity = nothing` (§12.2.6) and as the Beer–Lambert surface
transmissivity (§8.3); surface emission always uses
$\varepsilon_{\text{eff}}$.

### 5.8.3 Effective absorption coefficient

The volumetric (extinction) absorption coefficient used by the in-depth radiation
models (Beer–Lambert, and the experimental P1 model, §8) is a **concentration**-weighted sum of the
component *mass* absorption coefficients `α_j` [m²·kg⁻¹]:

$$
\alpha_{\text{eff}} \;=\; \sum_j \xi_j\, \alpha_j \quad [\text{m}^{-1}],
\tag{5.23}
$$

`effective_absorption(ξ, material)` → `_sum_absorption`
(`effective_properties.jl`). The dimensional bookkeeping is the key
point: a *mass* absorption coefficient `α_j` [m²·kg⁻¹] times a mass concentration
`ξ_j` [kg·m⁻³] yields a *volumetric* absorption coefficient [m⁻¹], which is what
the Beer–Lambert optical thickness `τ = α Δz` and the P1 diffusion coefficient
`D_rad = 1/(3α)` require. Because the weight is the raw concentration `ξ_j` (not a
normalized fraction), `α_eff` correctly scales down as material is consumed — a
thinning, depleting layer becomes optically thinner. Component absorption
coefficients default to `0.0` (no absorption) and carry no temperature
dependence.

---

## 5.9 Gas viscosity (Sutherland's law)

The dynamic viscosity of the escaping gas, needed for Darcy's law (§9), is the
material's `gas_viscosity` property. Its default is **Sutherland's law** for air,

$$
\mu(T) \;=\; \mu_{\text{ref}}
\left(\frac{T}{T_{\text{ref}}}\right)^{3/2}
\frac{T_{\text{ref}}+S}{T+S},
\tag{5.24}
$$

carried by `SutherlandProperty()` (§4.4.7) with the **canonical Sutherland air
constants** (e.g. White, *Viscous Fluid Flow*):

| Parameter | Default | Source |
|---|---|---|
| `μ_ref` | `1.716 × 10⁻⁵` Pa·s | `SutherlandProperty()`, `property_functions.jl` |
| `T_ref` | `273.15` K | `SutherlandProperty()`, `property_functions.jl` |
| `S` (Sutherland constant) | `110.4` K | `SutherlandProperty()`, `property_functions.jl` |

These give `μ(300 K) = 1.846 × 10⁻⁵` Pa·s and `μ(600 K) = 3.016 × 10⁻⁵` Pa·s.

> **Caution (reference temperature).** The viscosity reference temperature here
> is **`T_ref = 273.15 K`**, the temperature at which `μ_ref` is quoted for air —
> this is a *local default of the viscosity correlation* and is distinct from the
> global enthalpy datum `T_REF = 298.15 K` (§2). Do not conflate the two;
> `SutherlandProperty()` carries `T_ref = 273.15` in its own field.

Sutherland's law is the Chapman–Enskog rigid-sphere result for a single gas and
captures the *increase* of gas viscosity with temperature (`μ ∝ T^{1/2}` in the
high-`T` limit), the opposite of liquid behaviour. It is valid roughly over
`200–2000 K` for air-like gases.

### Material-supplied viscosity

The air constants are only the **default**. Each `Material` carries its own
pore-gas viscosity curve in the `gas_viscosity` field, and the Darcy path reads
it through `effective_gas_viscosity(T, material)`, which evaluates
`material.gas_viscosity(T)` on the property bracket and floors the result at
`MU_FLOOR`. A material that says nothing
carries `SutherlandProperty()` — the air triple above — so this reduces to
(5.24).

The field accepts anything `to_property` understands, so a fitted correlation
drops in directly:

```julia
# Cantera mixture-averaged fit for wood pyrolyzates, ascending coefficients,
# valid 280–1200 K.
Material(name = :WhitePine, components = (...), reactions = (...),
         gas_viscosity = PolynomialProperty(mu_coeffs))
```

**Why the guards matter.** μ enters the model *only* as the Darcy mobility
κ/μ (§9), never on its own. An error in μ is therefore exactly the reciprocal
error in permeability — and a μ that reaches zero or changes sign does not
degrade the answer, it reverses or blows up the gas velocity. Hence:

- the curve is a property law, so it evaluates on the property bracket
  $[200, 3000]$ K (§4.4): diverging Newton trial states cannot kill the solve
  on `T^{3/2}` or evaluate a fitted polynomial absurdly far outside its range
  (the Chapman–Enskog diffusivity formula, which is not a property law, keeps
  its own 100 K floor `T_GAS_PROPERTY_FLOOR`);
- the result is floored at `MU_FLOOR` (1e-12 Pa·s, seven orders below any real
  gas) so an extrapolated fit cannot go non-positive;
- `Material`'s constructor probes the curve over 200–1800 K and refuses one that
  is already non-positive or non-finite there, in the spirit of `MixingSpec`'s
  `CARMAN_KOZENY`/`length_scale` check.

Because the bracket belongs to every property law and the floor is applied in
`effective_gas_viscosity(T, material)`, every curve a material can carry is
protected, not just the built-in Sutherland form.

**Accuracy.** Air is within ~2 % for an N₂-swept pore space, but overstates the
viscosity of wood pyrolyzates by roughly 1.7× at 300 K falling to 1.3× at
1000 K — a bias in one direction everywhere, growing as the pore space fills
with volatiles.

**Composition dependence.** `mixture_gas_viscosity(ξ, T, material)` takes `ξ`
and every Darcy call site passes it, but the function returns the
material's own curve regardless: the model carries **one curve per material**, an average over the
composition that material is expected to produce. It cannot resolve the pore gas
shifting from background air to volatiles as pyrolysis proceeds. A
composition-dependent rule (Wilke) would be introduced by changing that one
function — but see §5.16 item 1 for why that is not a local change.

---

## 5.10 Chapman–Enskog binary diffusivity

For materials where an explicit per-component `λ_j` is not supplied, the package
can compute a first-principles binary diffusion coefficient from
**Chapman–Enskog kinetic theory** (Bird, Stewart & Lightfoot). For a species `A`
diffusing through a background gas `B` (taken as air),

$$
D_{AB} \;=\; 5.95\times10^{-4}\,
\frac{T^{3/2}\sqrt{\,1/M_A + 1/M_B\,}}{P\,\sigma_{AB}^2\,\Omega_D(T^*)}
\quad [\text{m}^2\text{s}^{-1}],
\tag{5.25}
$$

with `T` in K, `M_A, M_B` in **kg·mol⁻¹**, `P` in **Pa**, `σ_{AB}` in **Å**, and
`D_{AB}` in **m²·s⁻¹**.

`chapman_enskog_diffusivity(T, M; P, σ, ε_k)`
(`effective_properties.jl`). The prefactor and unit handling deserve
care, because the *textbook* Chapman–Enskog constant is `1.8583 × 10⁻³` in
**cgs-mixed** units (`T` in K, `M` in g·mol⁻¹, `P` in atm, `σ` in Å, `D` in
cm²·s⁻¹), stored as `CHAPMAN_ENSKOG_PREFACTOR = 1.8583e-3`
(`effective_properties.jl`). The implementation converts inputs and output to
SI explicitly (`effective_properties.jl`), and **every** conversion
contributes to the lumped prefactor — not just the `D` output conversion:

- `M`: kg·mol⁻¹ → g·mol⁻¹ (`×1000`). This appears *inside* `√(1/M)`, so it
  introduces a factor `1/√1000` in the diffusivity.
- `P`: Pa → atm (`/101325`). Because `P` is in the denominator, using SI `P` in
  Pa introduces a factor `×101325`.
- `D`: cm²·s⁻¹ → m²·s⁻¹ (`×1e-4`).

The net SI prefactor is therefore

$$
1.8583\times10^{-3}\;\times\;10^{-4}\;\times\;\frac{101325}{\sqrt{1000}}
\;\approx\; 5.95\times10^{-4},
$$

which is the constant quoted in the SI form of (5.25) (`M` in kg·mol⁻¹, `P` in
Pa). Folding *only* the `×1e-4` output conversion into the cgs prefactor would
give `1.8583 × 10⁻⁷`, but that value is valid **only** when the inputs are kept
in their cgs-mixed units (`M` in g·mol⁻¹, `P` in atm); pairing `1.8583 × 10⁻⁷`
with the SI symbols defined for (5.25) is dimensionally inconsistent and
under-predicts `D` by a factor `101325/√1000 ≈ 3.2 × 10³`. The physical scalings
are `D ∝ T^{3/2}` (through the explicit `T^{1.5}` and weak `Ω_D` dependence),
`D ∝ 1/P`, and `D ∝ √(1/M_A + 1/M_B)` (lighter species diffuse faster).

### 5.10.1 Combining rules and the collision integral

The binary interaction parameters use the standard Lennard-Jones combining rules
(`effective_properties.jl`):

$$
\sigma_{AB} = \tfrac12(\sigma_A + \sigma_B), \qquad
\varepsilon_{AB} = \sqrt{\varepsilon_A\,\varepsilon_B},
$$

and the **reduced temperature** is `T* = T/(ε_AB/k_B)`
(`effective_properties.jl`). The **diffusion collision integral** `Ω_D(T*)`
is evaluated with the Neufeld et al. (1972) correlation,

$$
\Omega_D(T^*) =
\frac{A}{(T^*)^B} + \frac{C}{e^{D T^*}} + \frac{E}{e^{F T^*}} + \frac{G}{e^{H T^*}},
\tag{5.26}
$$

`collision_integral_diffusion(T*)` (`effective_properties.jl`), with the
coefficients verified against the source (`effective_properties.jl`):

| | | | |
|---|---|---|---|
| `A = 1.06036` | `B = 0.15610` | `C = 0.19300` | `D = 0.47635` |
| `E = 1.03587` | `F = 1.52996` | `G = 1.76474` | `H = 3.89411` |

The fit is valid for `0.3 ≤ T* ≤ 400` with <0.1 % error.

> **Note (Neufeld G).** The coefficient `G` is **`1.76474`**, as in the code
> and the original Neufeld, Janzen & Aziz (1972) correlation.

### 5.10.2 Heuristic Lennard-Jones defaults and their unvalidated status

When the user does not supply `σ` or `ε/k_B` for a species, the code estimates
them from the molar mass via rough correlations anchored to small organics
(`effective_properties.jl`):

$$
\sigma_{\text{species}} \approx 3.5 + 0.5\log_{10}\!\Bigl(\max(M/0.016,\,1)\Bigr)\ [\text{Å}],
\qquad
\frac{\varepsilon}{k_B} \approx 150 + 100\,\frac{M}{0.044}\ [\text{K}],
\tag{5.27}
$$

with the air background fixed at `σ_air = 3.617 Å`, `ε/k_B|_air = 97 K`,
`M_air = 0.029 kg·mol⁻¹` (`effective_properties.jl`). These are **rough
estimates** suitable for organic vapours; they have **not** been validated
against tabulated LJ parameters for the specific multi-component gas mixtures that
arise in pyrolysis, and precise values should be supplied when available (§5.16).

### 5.10.3 Conversion to a gas-transfer coefficient

The model's diffusive flux is written in terms of a transfer coefficient `λ`
[m²·s⁻¹] with the convention `J = −λ ∂ξ/∂x` (the gas-density factor is folded
into `λ`), so for gas-in-gas diffusion `λ ≈ D` is the consistent
identification: assign the `chapman_enskog_diffusivity` value (optionally
porosity/tortuosity-corrected, and optionally wrapped in a temperature-dependent
property function) as the gaseous component's `λ`. There is no separate
conversion routine — the identification is one-to-one (§9 derives the precise
volume-fraction-gradient flux that `λ` multiplies).

---

## 5.11 No per-species transport scaling

The effective gas-transfer coefficient `λ_eff` of §5.6 is a property of the
**mixture**, and the gas flux kernels apply the same face value to every gas
species (§9.2). There is no kinetic-theory per-species scaling of the diffusive
flux: a light gas and a heavy gas in the same simulation diffuse identically.
Species-faithful magnitudes enter only through the values the user assigns to
the component `λ_j(T)` properties — for which the Chapman–Enskog estimate of
§5.10 is the first-principles starting point. (The nomenclature symbol `M_ref`
appears only inside the Chapman–Enskog combining rules of §5.10.)

---

## 5.12 Face averaging: distance-weighted harmonic mean

All cell-centered transport properties (`k`, `λ`, `κ`) must be transferred to
**faces** to form the finite-volume fluxes (§13). For a conserved flux crossing a
material interface, the correct face value is the one that preserves flux
continuity, which makes the property behave as a **resistance in series**. With
distances `d_L, d_R` from the two cell centers to the face, the
**distance-weighted harmonic mean** is

$$
\varphi_{\text{face}}
\;=\; \frac{d_L + d_R}{\,d_L/\varphi_L + d_R/\varphi_R\,},
\tag{5.29}
$$

`harmonic_mean_weighted(a, b, d_a, d_b)`
(`effective_properties.jl`). The derivation: the conductance of each
half-cell along the flux path is `φ/d`, the two resistances `d_L/φ_L` and
`d_R/φ_R` add, and the effective face property over the total span `d_L + d_R` is
(5.29). On a uniform mesh `d_L = d_R` and (5.29) reduces to the standard harmonic
mean `2 φ_L φ_R/(φ_L + φ_R)`. On a stretched mesh (e.g. boundary refinement) the
distance weighting corrects the averaging and improves boundary fluxes by ~1–5 %.

> **Implementation note.** `update_properties!` applies (5.29) to fill `k_face`
> and `λ_face` for every interior face; at a boundary face (where one
> `cell_id == 0`) it simply copies the adjacent interior cell's value
> (`property_cache.jl`). The distances are computed from the live state
> positions `d_L = |z_face − z^c_L|`, `d_R = |z^c_R − z_face|`, so the averaging
> is correct under ALE mesh motion. A small-value guard returns 0 if either input
> is below `1e-20` (`effective_properties.jl`), preventing `Inf`/`NaN`
> from a vanished cell. The harmonic mean is applied to `λ` itself, **not** to
> the product `λρ` that a mass-fraction-gradient discretization
> `J = −λρ ∂Y/∂x` would average: the volume-fraction-gradient flux used by the
> gas transport (§9) carries the density dependence in the `ξT` driving
> potential and averages `λ` alone.

The same resistance-in-series averaging is used for permeability at faces
(§5.7.4) and for the P1 radiation diffusion coefficient `D_rad` at faces (§8).

---

## 5.13 The property cache

Effective properties are not recomputed inside every inner flux loop; they are
evaluated **once per cell per residual** and stored in a `PropertyCache`, then
read by the conduction, transport, radiation, and boundary kernels.

### 5.13.1 Layout

`PropertyCache{Tu, DarcyT}` (`property_cache.jl`) is parameterized on the
live scalar type `Tu` so the entire cache holds `ForwardDiff.Dual` values
unchanged during automatic differentiation (§17), and on `DarcyT` so the
optional Darcy sub-cache is either present (`DarcyCache{Tu}`) or `Nothing`. The
per-cell arrays (length `n`) are

| Field | Symbol | Units | Built from |
|---|---|---|---|
| `ρ` | mixture density | kg·m⁻³ | `mixture_density` (5.6) |
| `ρcp` | `(ρc_p)_eff` | J·m⁻³·K⁻¹ | `effective_heat_capacity` (5.9) |
| `k` | `k_eff` | W·m⁻¹·K⁻¹ | `_conductivity_from_v` (§5.5) |
| `λ` | `λ_eff` | m²·s⁻¹ | `_gas_transfer_from_v` (§5.6) |
| `ε` | `ε_eff` | – | `_emissivity_from_v` (5.22) |
| `a` | `a_eff` | – | `_absorptivity_from_v` (5.22b) |
| `α` | `α_eff` | m⁻¹ | `effective_absorption` (5.23) |

and the per-face arrays (length `n+1`) are `k_face` and `λ_face` (5.29). The
**Darcy sub-cache** `DarcyCache{Tu}` (`property_cache.jl`) is allocated only
for pressure-driven modes (`with_darcy = true`) and adds per-cell `P, φ, κ, μ` and
per-face `v_gas, ∇P` — keeping the four extra per-cell and two per-face arrays out
of memory for the common Fickian-only configuration. Note `ρcp` caches the
mixture `(ρc_p)_eff` of (5.9) — the storage term the energy
RHS divides by (§7).

### 5.13.2 The update kernel

`update_properties!(cache, state, mesh, material)`
(`property_cache.jl`) is the single entry point. It (i) resizes the cache
if its length differs from the state's; (ii) loops over the live cells
`1:n_active_cells` (parked ghost slots, §14.3, are left untouched), building `v_full =
volume_fractions(ξ_i, T_i, …)` and `v_solid = mixing_volume_fractions(ξ_i, T_i,
…)` (the swelling-weighted basis with the depleted-cell fallback, §5.1.3)
once each and threading them into the effective-property helpers; and (iii)
loops over faces applying the distance-weighted harmonic mean. The cell loop is
exactly:

```julia
cache.ρ[i]   = mixture_density(ξ_i, T_i, material)
cache.ρcp[i] = effective_heat_capacity(ξ_i, T_i, material)
cache.k[i]   = _conductivity_from_v(v_solid, T_i, material; k_init = k_prev)
cache.λ[i]   = _gas_transfer_from_v(v_full, T_i, material)
cache.ε[i]   = _emissivity_from_v(v_solid, T_i, material)
cache.a[i]   = _absorptivity_from_v(v_solid, T_i, material)
cache.α[i]   = effective_absorption(ξ_i, material)
```

This is the authoritative statement of which basis each property uses:
**conductivity, emissivity, and surface absorptivity from `v_solid`, gas
transfer from `v_full`**
(§5.1.3). The previous step's conductivity `k_prev = cache.k[i]` is passed as the
Bruggeman warm-start `k_init` (§5.5.3). After sizing, the loop is **allocation
free**: all composition-dependent dispatch is generated (compile-time unrolled),
all intermediate tuples are statically sized, and the Bruggeman Newton iteration
reuses scalars.

`resize_cache!` grows all arrays (and the Darcy sub-cache, if
present) in place after a mesh change.

---

## 5.14 Comparison to ThermaKin2Ds, Gpyro, and FDS

The mixing-rule choices map cleanly onto the established condensed-phase pyrolysis
codes. Mapping their notation onto ours on first use: ThermaKin's stoichiometric
`θ_i^j` is our `ν_{i,j}`, Gpyro's volume fraction `X_α` is our `v_j`, Gpyro's
mixture average `\barρ` corresponds to our swelling mixture density `ρ`, and the
FDS/Gpyro pre-exponential `Z` is our `A_i` (§2).

**Effective heat capacity / gas storage.** ThermaKin (Eq. 17) keeps the gas
sensible-storage term — it sums `Σ_j ξ_j c_j` over *all* components, gas included
— and Pyrolysis.jl does the same in (5.9). Gpyro's default thermal-equilibrium
form drops gas storage from `ρc_p` (its `FULL_QSG` option restores it) and FDS
carries no in-depth gas; both are approximations of (5.9) that Pyrolysis.jl
does not make (§5.4.2).

**Conductivity / emissivity basis.** ThermaKin determines surface properties from
the *solid* matrix: "when γ = 0, gaseous components do not contribute to the
thermal conductivity" (Eqs. 8–10). Pyrolysis.jl reproduces this with the
**swelling-weighted volume-fraction basis** (5.4) for both `k_eff` and `ε_eff`;
the `γ_j·ξ_j/ρ_j` weighting reduces exactly to ThermaKin's `γ = 0` rule for
canonical gases (zero weight) and condensed phases (full weight).

**Conductivity blend.** ThermaKin blends the series and parallel limits with a
weighting parameter — our `WEIGHTED` rule (5.14) with `β` is the direct analogue
(`k = β k_∥ + (1−β) k_series`), default `β = 0.5`. Gpyro homogenizes by
volume-fraction averaging (its `X_α`-weighted properties) with power-law
temperature dependence built into the component properties; our `PARALLEL` rule
(5.13, the volume-weighted arithmetic mean) is the volume-fraction average, and
arbitrary `k_j(T)` (including power laws) is supported through the property
functions (§4). The **BRUGGEMAN** effective-medium option (5.15) is an addition
beyond the volume-average/blend menu of ThermaKin and Gpyro, motivated by the
self-similar char microstructure.

**Permeability.** Gpyro uses a porosity-dependent permeability for its Darcian
gas transport; our **Carman–Kozeny** closure (5.20) is the canonical
porosity-based form (Kozeny constant 180), and is offered alongside component
mixing as a mutually exclusive alternative.

**Radiation properties.** FDS works with a solid absorption coefficient `κ_s`
(our `α`, §2) and a surface emissivity/absorptivity tied by Kirchhoff's law; our
`α_eff = Σ ξ_j α_j` (5.23) and `ε_eff` (5.22) provide exactly these, with the
`α = ε_eff` Kirchhoff identity used at the boundary (§12).

---

## 5.15 Mixing derivatives for the Jacobian

The structured Jacobian (§15) needs the sensitivity of each effective property to
the state `(T, ξ)`. No mixing rule carries a hand-coded derivative: the Jacobian
operators evaluate the mixing functions on `ForwardDiff.Dual` arguments, so every
derivative follows the closure as implemented.

### 5.15.1 Permeability and porosity derivatives

The `DarcyTransportOperator` (§15) differentiates `effective_permeability` and
`porosity` per cell with a single stack-allocated `Dual` pass over `(T, ξ)`
(`NC+1` partials) and feeds the results into its analytical Darcy-flux Jacobian.
The mixing rule therefore has a single source of truth: the derivatives track
every closure — the renormalized-solid-basis PARALLEL/SERIES/WEIGHTED mixing
including the impermeable-phase series hard-zero (where the derivative is
exactly 0 on the zero branch), the CARMAN_KOZENY porosity chain, and any
user-supplied `κ_j(T)` (so `∂κ/∂T` is exact for temperature-dependent κ).
BRUGGEMAN is rejected for permeability (§5.7.3).

### 5.15.2 Conductivity, gas-transfer, and other derivatives

For `k_eff`, `λ_eff`, `ε_eff`, `a_eff`, and `α_eff`, the temperature dependence
enters through the component laws `k_j(T)`, `λ_j(T)`, `ε_j(T)` (§4.5) and the
composition dependence through the volume fractions (5.3)/(5.4) (for `α_eff`,
through ξ directly). The conduction, gas-transport and Beer–Lambert operators
evaluate these mixing functions once per cell on `Dual` arguments over `(T, ξ)`
and compose the face fluxes or the transmission chain by hand in their
analytical Jacobians; the reaction-source and boundary operators run ForwardDiff
over stencil-local kernels (§15). The Bruggeman `k_eff` is differentiated by AD
through its Newton solve (started from `k_parallel` in the Jacobian, warm-started
from the cached `k` in `update_properties!`),
which is differentiable because every operation in the iteration is smooth (the
step limiter is inactive and the per-term branch selection is locally constant
at a converged interior solution).
This AD-cleanliness — together with the smooth `tanh` kinetics gates (§6) and the
clamped property laws (§4.4) — is what makes the end-to-end sensitivities of §17
well-defined.

**Warm-start caveat.** The Bruggeman iteration's convergence test is on the
*value* only, and in an AD pass the Dual-typed cache means the warm start
`k_init` carries the *previous* step's partials. This is benign in the normal
path: at least one Dual-propagating Newton update always executes before the
tolerance break, and by the implicit-function property a single update
annihilates stale partials to O(tol) — the returned derivatives are those of
the implicit closure, not of the stale seed. The one failure mode is the
pre-update `|df| < 1e-20` degenerate-input break, which can return before any
update has run and thus pass the stale partials through; it requires
degenerate inputs (all contrasts ≈ 0) that do not arise from valid materials.

---

## 5.16 Limitations and open issues

The following modelling limitations of the effective-property subsystem are
recorded from the source; several are flagged here for the
verification and validation discussion of §18.

1. **Single-value mixture gas viscosity.** `mixture_gas_viscosity` ignores
   composition and returns the material's single viscosity curve. The
   *material* dependence is resolved — a material supplies the blend correlation
   for the gas it produces — but the *composition* dependence within a material
   is not: one curve averages over the pore gas shifting from background air to
   volatiles as pyrolysis proceeds. A Wilke mixing rule over per-component
   μⱼ(T) would resolve that transition.

   Two things make this more than a one-function change, and both are worth
   knowing before attempting it:

   - The production Darcy Jacobian is the **analytical override** in
     `Jacobian/operators/darcy_transport.jl`, not the AD path (which is retained
     only as the rtol = 1e-12 equivalence reference). Its mobility chain rule
     assumes `dμ/dξ = 0` and takes the ξ-path through κ alone; every
     `dvg_dξ_*_j` expression would need a `− mob·dμ_dξ[j]` term.
   - That equivalence gate **cannot catch the omission** while every Darcy test
     fixture carries a single gas component, because a one-species mixture rule
     is the identity. A multi-gas Darcy fixture has to land first.

   `test/discretization/test_darcy.jl` asserts ξ-invariance of
   `mixture_gas_viscosity` as a tripwire, so the day this changes the failure is
   loud and points here.

2. **Unvalidated Chapman–Enskog LJ defaults.** The heuristic `σ(M)` and
   `ε/k_B(M)` correlations (5.27) are rough estimates anchored to small organics
   and have not been validated against literature LJ parameters for specific
   pyrolysis gas species. Supply explicit `σ`, `ε_k` where known.

3. **Permeability microstructure evolution.** Real materials may exhibit
   permeability changes with thermal stress and cracking beyond what the
   `κ_j(T)` property functions and porosity-based closures can represent;
   the Darcy Jacobian differentiates whatever closure is supplied
   (including `κ_j(T)` and the Carman–Kozeny porosity chain), but the closure
   itself remains the modelling limitation.

4. **Intrinsic-density temperature dependence.** Intrinsic (skeletal) density
   shares the same functional form as bulk density. In reality the cell-wall
   material may have a different thermal-expansion behaviour, which would shift
   the porosity (5.7) and hence the pressure and Carman–Kozeny permeability.

5. **Bruggeman near percolation.** The Newton iteration (5.15)–(5.16) is
   safeguarded by the 90 % step limiter, the insulating-phase exact limit, and
   the below-threshold collapse return (§5.5.3). If a continuous phase
   approaches zero volume while its conductivity diverges relative to the
   others, convergence can still stall; the step limiter is heuristic, not a
   percolation-aware bound.

6. **Carman–Kozeny microstructure assumption.** The `φ³/(1−φ)²` closure with
   Kozeny constant 180 is empirical for a packed bed of spheres. Real porous char
   may have a fractal or anisotropic pore structure for which a different closure
   (or a measured `κ(φ)` curve) would be more accurate.

The matching subsystem discussions are: the energy assembly that consumes
`(ρc_p)_eff` and the face conductivities — §7 (Heat Conduction, Convection &
Energy Assembly); the gas fluxes that consume `λ_eff`, `λ_j`, `κ_eff`, `μ`, and
`φ` — §9 (Gas-Phase Mass Transport & Pressure); the radiation models that consume
`α_eff` and `ε_eff` — §8 (Thermal Radiation); the swelling/volume-change use of
`ρ`, `γ_j` — §10 (Volume Change and Lateral Shrinkage); and the structured
assembly of every mixing derivative — §15 (Temporal Integration and Structured
Jacobian).


---

# 6. Reaction Kinetics and Stoichiometry

This chapter develops the complete condensed-phase reaction model of Pyrolysis.jl:
the Arrhenius rate law on a mass-concentration basis, the smooth temperature gates
and depletion limiter that replace the hard step functions of classical pyrolysis
codes, the stoichiometric assembly of species sources, and the heat-of-reaction
energy source with its (deliberately storage-oriented) sign convention. Every
governing equation, default value, and sign is derived from first principles and
then matched line-for-line to the implementation in `src/Physics/kinetics.jl` and
`src/Materials/reactions.jl`. We close with a careful contrast against ThermaKin2Ds,
Gpyro, and FDS, and a statement of the model's present limitations.

Throughout this chapter the index `i` denotes the **reaction index** and `j` the
**component (species) index** — the kinetics convention of overload note I1 in
§2, Nomenclature. The cell index, where it appears in the integration discussion,
is written explicitly. Mass concentration `ξ_j` (kg·m⁻³) is the primary species
state variable; stoichiometric coefficients are `ν_{i,j}` (overload note A2); the
Arrhenius pre-exponential is `A_i` (overload note A3); and the heat of reaction
`h` uses the storage sign convention `h > 0` endothermic (overload note H1).

---

## 6.1 The Arrhenius Reaction Rate Law

### 6.1.1 Continuous form

A single condensed-phase decomposition reaction $i$ proceeds at a volumetric rate
$r_i$ (kg·m⁻³·s⁻¹) given by an Arrhenius rate constant multiplied by a power-law
dependence on the mass concentrations of its reactants and gated to a temperature
window:

$$
r_i(T,\boldsymbol{\xi}) \;=\; A_i \,
\exp\!\left(-\frac{E_i}{R_g\,T}\right)
\;\prod_{j\,\in\,\mathrm{reactants}(i)} \xi_j^{\,n_{i,j}}
\,\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr)
\;\cdot\; H(T - T_{\min,i})\,H(T_{\max,i} - T).
\tag{6.1}
$$

Here

- $A_i$ is the pre-exponential factor (s⁻¹ for a first-order reaction; more
  generally with units that make $r_i$ come out in kg·m⁻³·s⁻¹ — see §6.4),
- $E_i$ is the activation energy (J·mol⁻¹),
- $R_g = 8.31446261815324$ J·mol⁻¹·K⁻¹ is the universal gas constant
  (`R_GAS`, src/Materials/constants.jl),
- $n_{i,j}$ is the reaction order of reaction $i$ with respect to reactant $j$
  ($n_{i,j} \ge 0$, enforced at construction),
- $\tanh(\xi_j/\xi_{\mathrm{th}})$ is the smooth **depletion rate factor**
  (§6.3), and
- $H(\cdot)$ is the Heaviside step, implemented as a **smooth tanh ramp** (§6.2).

The factorisation of $r_i$ into a temperature-only Arrhenius constant
$k_i(T) \equiv A_i \exp(-E_i/R_g T)$ and a composition-only product
$\prod_j \xi_j^{n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ is exact and is exploited
directly by the code:

```julia
k = rxn.A * exp(-rxn.E / (R_GAS * T))
```

The rate constant carries the classical Arrhenius physical interpretation: $A_i$
is the collision/attempt frequency (modulated by an entropic steric factor), and
$\exp(-E_i/R_g T)$ is the Boltzmann fraction of molecular configurations that
surmount the activation barrier $E_i$. Both $A_i$ and $E_i$ are stored as plain
`Float64` scalars on the `Reaction` struct (src/Materials/reactions.jl), so
they are temperature-independent constants — see §6.10 for the consequences.

**Mass-concentration (not molar) basis.** The composition dependence in (6.1) is
written in terms of the **mass concentration** $\xi_j$ (kg of component $j$ per m³
of mixture), not molar concentration. This is the decisive modelling choice that
distinguishes a pyrolysis solver from a gas-phase combustion code. In a
decomposing solid the molar mass of the "reactant" (e.g. virgin wood, a polymer)
is ill-defined — the condensed phase is a continuum of macromolecules, not a
well-mixed ideal gas of identifiable molecules. Working in $\xi_j$ sidesteps the
need for a molar mass of the solid, makes $r_i$ directly the mass conversion rate,
and lets the stoichiometric coefficients $\nu_{i,j}$ be pure **mass** yields that
sum to unity (§6.5). The price is that $A_i$ inherits whatever units make the
product dimensionally consistent (§6.4); this is the universal convention of the
condensed-phase pyrolysis literature (ThermaKin, Gpyro, FDS all do the same — see
§6.9).

### 6.1.2 Discrete (per-cell) form

The finite-volume discretisation (§13, Finite-Volume Discretization) collocates
all state at cell centres. Within `compute_rhs!` the reaction model is evaluated
once per cell from the cell-centred state $(T, \xi_1,\dots,\xi_{N_C})$
(src/Discretization/finite_volume.jl):

```julia
T = state.T[i]
ξ = ntuple(j -> state.ξ[j, i], Val(NC))
evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)
```

so the continuous law (6.1) becomes simply $r_i^{\,\mathrm{cell}} =
r_i(T_{\mathrm{cell}}, \boldsymbol{\xi}_{\mathrm{cell}})$ — there is **no spatial
operator inside the reaction term**. Reactions are local (point-wise in space);
all spatial coupling enters through transport (conduction §7, gas diffusion/Darcy
§9) and through the ALE mesh motion (§11). This is the cell-as-batch-reactor view:
each control volume behaves, kinetically, like an isolated TGA sample, identical
to the Gpyro and FDS interpretation (§6.9).

The single-reaction kernel is `reaction_rate` (src/Physics/kinetics.jl). It
is fully scalar-type-generic (`<:Real` arguments), so `ForwardDiff.Dual` numbers
propagate cleanly through it for Jacobian and sensitivity work (§15, §17).

> **Implementation note.** Single-reaction rate: `reaction_rate`,
> src/Physics/kinetics.jl. All-reaction evaluation: `evaluate_reactions!`
> dispatching to the `@generated` kernel `_evaluate_reactions_gen!`,
> src/Physics/kinetics.jl.

---

## 6.2 Smooth Temperature Gates

### 6.2.1 Why hard steps are unacceptable

The temperature window $[T_{\min,i}, T_{\max,i}]$ serves two purposes (following
ThermaKin's original motivation): it spares the solver from evaluating a reaction
where it is kinetically negligible, and it lets a reaction be used to represent a
thermodynamically-controlled phase transition that is sharply switched on at a
threshold. The textbook expression of this is the product of two Heaviside steps,
$H(T - T_{\min,i})\,H(T_{\max,i} - T)$, equal to 1 inside the window and 0 outside.

A hard step is, however, **non-differentiable** at the threshold and has a zero
derivative everywhere else. Both properties are poison for the implicit machinery
of Pyrolysis.jl:

1. The stiff time integrators solve nonlinear systems with Newton's method, whose
   convergence depends on a meaningful Jacobian $\partial r_i/\partial T$. At a
   hard step the Jacobian is a Dirac delta (numerically: undefined), and away from
   it the step contributes nothing — Newton sees a flat rate that suddenly jumps.
2. Automatic differentiation (`ForwardDiff`) of a discontinuous function returns a
   derivative of zero almost everywhere and is simply wrong at the breakpoint,
   corrupting both the structured Jacobian (§15) and forward-mode sensitivities
   (§17).

### 6.2.2 The tanh ramp

Pyrolysis.jl replaces each Heaviside by a hyperbolic-tangent ramp. The lower gate
is applied only when $T_{\min,i}$ is finite and strictly positive; the upper gate
only when $T_{\max,i}$ is finite:

$$
g_{\min,i}(T) =
\begin{cases}
1, & T_{\min,i}\ \text{not finite or}\ \le 0,\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh\!\big((T - T_{\min,i})/w_{\text{gate},i}\big)\bigr), & \text{otherwise},
\end{cases}
\qquad
g_{\max,i}(T) =
\begin{cases}
1, & T_{\max,i}\ \text{not finite},\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh\!\big((T_{\max,i} - T)/w_{\text{gate},i}\big)\bigr), & \text{otherwise},
\end{cases}
\tag{6.2}
$$

and the gated rate is $r_i \leftarrow r_i\, g_{\min,i}(T)\, g_{\max,i}(T)$, with
$w_{\text{gate},i}$ the per-reaction **gate width** (`Reaction(...;
gate_width)`, default $1.0$ K, validated finite and $>0$). The implementation
is verbatim:

```julia
if isfinite(rxn.T_min) && rxn.T_min > 0.0
    r *= 0.5 * (1 + tanh((T - rxn.T_min) / rxn.gate_width))
end
if isfinite(rxn.T_max)
    r *= 0.5 * (1 + tanh((rxn.T_max - T) / rxn.gate_width))
end
```

(src/Physics/kinetics.jl). The default gates are $T_{\min}=0$ and
$T_{\max}=\infty$ (src/Materials/reactions.jl), so an
unconstrained reaction incurs no gating cost at all — the `T_min > 0.0` guard
short-circuits the lower ramp and `isfinite(Inf) == false` short-circuits the
upper ramp.

The ramp $\tfrac12(1 + \tanh(\Delta T/w_{\text{gate}}))$ rises from $\approx 0$
to $\approx 1$ over a width set by $w_{\text{gate}}$. At the default
$w_{\text{gate}} = 1$ K the transition has a characteristic half-width of
$\mathcal O(1\text{–}3)$ K and reaches $\tfrac12(1+\tanh 5)\approx 0.99995$
only about 5 K inside the window; equivalently, the gate is
$\approx 4.5\times10^{-5}$ about 5 K **outside** the window. All widths scale
linearly with $w_{\text{gate}}$.

### 6.2.3 The ~10·w_gate leakage caveat

The smoothing is not free: the reaction is "weakly active" within roughly
$10\,w_{\text{gate}}$ of each boundary, where a hard step would have it
identically off. At the lower edge
$\tfrac12(1+\tanh((T - T_{\min})/w_{\text{gate}}))$ takes the value $0.5$
exactly at $T = T_{\min}$, $\approx 0.018$ at $T_{\min} - 2\,w_{\text{gate}}$,
and $\approx 4.5\times10^{-5}$ at $T_{\min}-5\,w_{\text{gate}}$. For
decomposition kinetics, where $T_{\min}$ is a numerical convenience rather
than a sharp physical threshold, this leakage is immaterial because the
Arrhenius factor is itself vanishingly small that far below the reaction's
effective onset. For a reaction used as a sharp phase-transition switch
(melting fronts, T-gated latent heat as in V-CPL-01b), narrow the gate
(`gate_width = 0.1` K makes the smear a tenth of the default) — every
Jacobian coupling differentiates through `reaction_rate` with ForwardDiff, so
the width threads through the structured Jacobian automatically, at the cost
of a stiffer gate-adjacent Jacobian entry. On a driving ΔT of only a few
kelvin the default 1 K ramp visibly shifts the effective switching point;
either narrow the gate or offset $T_{\min}$ to recentre the ramp.

> **Implementation note.** Temperature gates: `reaction_rate`,
> src/Physics/kinetics.jl. `gate_width` field, validation, and the default
> `T_min = 0.0`, `T_max = Inf` in the convenience constructors,
> src/Materials/reactions.jl.

---

## 6.3 The Depletion Limiter

### 6.3.1 The stiffness problem near zero concentration

As a reactant is consumed, $\xi_j \to 0$. For a first-order reaction ($n_{i,j}=1$)
the bare rate $r_i = k_i(T)\,\xi_j$ approaches zero linearly, which is benign. But
the **Jacobian entry** $\partial r_i/\partial \xi_j = k_i(T)$ does not vanish: at
high temperature $k_i(T)$ can be enormous ($A_i$ ranges over $10^8$–$10^{20}$ s⁻¹
in typical fits). A control volume holding a numerically tiny,
possibly slightly-negative-from-roundoff residue of reactant therefore presents
the implicit solver with a stiff, fast-relaxing mode whose eigenvalue is
$-k_i(T)$ even though essentially no mass remains to react. This manifests as
time-step throttling, Newton non-convergence, and — if $\xi_j$ overshoots below
zero and $n_{i,j}$ is non-integer — complex or `NaN` rates from $\xi_j^{n_{i,j}}$.

Sub-first-order fits make this strictly worse: for $0 < n_{i,j} < 1$ the bare
$\partial r_i/\partial\xi_j \sim \xi_j^{\,n-1}$ **diverges** at depletion, and for
a zero-order fit ($n_{i,j}=0$, common in TGA analysis) the bare rate does not
vanish at all — left unguarded it keeps consuming reactant at $k_i(T)$ straight
through $\xi_j = 0$ and integrates the concentration negative without bound.

### 6.3.2 The tanh rate factor

Pyrolysis.jl rolls the rate off smoothly as a reactant depletes. Each reactant
contributes a multiplicative **depletion factor to the rate** (not to the
power-law base):

$$
r_i \;=\; k_i(T)\prod_{j\,\in\,\mathrm{reactants}(i)}
\xi_j^{\,n_{i,j}}\;\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr),
\tag{6.3}
$$

with a hard floor at zero applied first: $\xi_j \leftarrow \max(\xi_j, 0)$. The
implementation:

```julia
ξ_raw = max(ξ[idx], zero(eltype(ξ)))
conc_factor *= _fast_pow(ξ_raw, order)
if limiter.enabled
    conc_factor *= _depletion_factor(ξ_raw, limiter.threshold)
end
```

(src/Physics/kinetics.jl). The `max(ξ, 0)` clip guarantees a non-negative
base for the power even if the integrator transiently produces $\xi_j < 0$, and
it pins the depletion factor at $\tanh 0 = 0$ there — an overshooting cell
consumes nothing. `_depletion_factor` short-circuits to `one(ξ)` above
$20\,\xi_{\mathrm{th}}$, which is **bit-identical** to evaluating the `tanh`
(Float64 `tanh(x)` rounds to exactly `1.0` for $x \gtrsim 19.1$), so the common
path — bulk concentrations orders of magnitude above the threshold — pays no
transcendental (src/Physics/kinetics.jl).

Applying the factor to the rate rather than substituting an effective
concentration into the power-law base matters for low orders: a base-applied
limiter leaves a zero-order reactant unguarded ($\xi_{\mathrm{eff}}^{\,0}\equiv 1$)
and, for $0 < n < 0.5$, behaves as $(\xi^2/\xi_{\mathrm{th}})^n$ near zero with a
still-singular $\partial r/\partial\xi \sim \xi^{\,2n-1}$. The multiplicative
form removes the singular eigenvalue for **every** order $n \ge 0$.

### 6.3.3 Behaviour and default values

The limiter is configured by the `DepletionLimiter` struct
(defined in src/Materials/reactions.jl; re-exported from `Pyrolysis.Physics`),
whose defaults are

| field | default | meaning |
|---|---|---|
| `threshold` ($\xi_{\mathrm{th}}$) | $1.0$ kg·m⁻³ | rate roll-off scale |
| `enabled` | `true` | set `false` to recover the bare Arrhenius rate |

Each `Material` carries its own limiter in `material.depletion_limiter`, set
with the `Material` constructor's `depletion_limiter` keyword (default
`DepletionLimiter()`); the solver's kinetics path `evaluate_reactions!` reads it
for every rate evaluation. `DEFAULT_DEPLETION_LIMITER = DepletionLimiter()`
(src/Physics/kinetics.jl) backs only the standalone `reaction_rate` entry point
when no `limiter` keyword is passed.

The factor $\tanh(\xi/\xi_{\mathrm{th}})$ has the key properties:

- **Smooth everywhere.** There is no threshold branch: the rate is $C^\infty$ in
  $\xi > 0$ (the only non-smooth point left is the `max(ξ, 0)` clip exactly at
  zero). Neither the rate nor $\partial r/\partial\xi$ jumps anywhere.
- $\xi \gg \xi_{\mathrm{th}}$: the bare rate is recovered double-exponentially —
  the relative deficit is $2e^{-2\xi/\xi_{\mathrm{th}}}$, i.e. $\sim\!3.7\%$ at
  $2\,\xi_{\mathrm{th}}$, $4\times10^{-9}$ at $10\,\xi_{\mathrm{th}}$, and exactly
  zero in Float64 beyond $\approx 19.1\,\xi_{\mathrm{th}}$.
- $\xi \to 0$: the rate vanishes as $k\,\xi^{\,n+1}/\xi_{\mathrm{th}}$ with
  $\partial r/\partial\xi \to 0$ for $n > 0$ and $\to k/\xi_{\mathrm{th}}$
  (finite) for $n = 0$ — the stiff eigenvalue is removed for every order
  $n \ge 0$. For $n = 1$ the limited rate is exactly
  $k\,\xi\tanh(\xi/\xi_{\mathrm{th}})$.

The robustness is bought with a rate **bias near depletion**: the suppression is
$\sim\!24\%$ at $\xi = \xi_{\mathrm{th}}$ ($\tanh 1 \approx 0.762$) and grows as
the reactant runs out — $90\%$ at $0.1\,\xi_{\mathrm{th}}$. For first-order
kinetics the limited rate behaves as $k\,\xi^2/\xi_{\mathrm{th}}$ near zero,
which converts the exponential burnout tail into an algebraic
$\xi \approx \xi_{\mathrm{th}}/(kt)$ tail. Choose $\xi_{\mathrm{th}}$ well below
the concentrations that matter for the quantities of interest.

Reaction orders are validated at construction: $n < 0$ is rejected
(`_validate_reaction_order`, src/Materials/reactions.jl), since a negative-order
rate is singular at depletion and no limiter can repair it.

For validation studies that need the unbiased Arrhenius rate, the limiter is
disabled per call:

```julia
r = reaction_rate(rxn, ξ, T; limiter = DepletionLimiter(enabled = false))
```

or for an entire solve, per material:

```julia
Material(name = …, components = …, reactions = …,
         depletion_limiter = DepletionLimiter(enabled = false))
```

Disabling removes the only depletion guard: zero-order reactions then consume
reactant through $\xi = 0$ without bound, and orders $0 < n < 1$ have a singular
$\partial r/\partial\xi$ at depletion.

> **Implementation note.** Limiter struct: `DepletionLimiter`,
> src/Materials/reactions.jl (carried per material in
> `material.depletion_limiter`; re-exported from Pyrolysis.Physics).
> Per-reactant application inside the concentration product loop:
> src/Physics/kinetics.jl.

---

## 6.4 Reaction Order and the Concentration Product

The composition factor is the product over all reactants of the clipped
concentration raised to its per-reactant order, times the per-reactant
depletion factor of §6.3:

$$
\Phi_i(\boldsymbol{\xi}) = \prod_{j\,\in\,\mathrm{reactants}(i)}
\xi_j^{\,n_{i,j}}\,\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr).
\tag{6.4}
$$

Each reactant is stored as an `(component_index, reaction_order)` pair in the
`reactants::NTuple{NRin, Tuple{Int,Float64}}` field
(src/Materials/reactions.jl), so the order $n_{i,j}$ is a **`Float64`** and may
be non-integer (e.g. $n=0.5$, $n=1.5$, common in TGA fits and in nucleation/
diffusion reaction models). Orders must satisfy $n_{i,j} \ge 0$ — negative
orders are rejected at construction (§6.3.3). The tuple is iterated with the
loop fully unrolled by the compiler because `NRin` is a static type parameter.

The exponentiation routes through a fast-path helper that special-cases the two
most common orders:

```julia
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)
```

(src/Physics/kinetics.jl). For $n=1$ it returns $x$ directly (no `pow` call);
for $n=0$ it returns `one(x)` — a zero-order reactant contributes no
concentration dependence to the power law, but its depletion factor
$\tanh(\xi/\xi_{\mathrm{th}})$ still gates the rate to zero as it runs out
(§6.3.2); otherwise it falls back to the general
`x^n`. Because the `max(ξ,0)` clip (§6.3) guarantees a non-negative base,
$x^n$ is always real-valued even for fractional $n$. The accumulator is
initialised at the promoted scalar type so Duals carry through:

```julia
conc_factor = one(promote_type(eltype(ξ), typeof(T)))
```

(src/Physics/kinetics.jl).

**Units of the pre-exponential.** Dimensional consistency of (6.1) forces

$$
[A_i] = \text{kg}\,\text{m}^{-3}\,\text{s}^{-1}
\Big/ \Big(\text{kg}\,\text{m}^{-3}\Big)^{\sum_j n_{i,j}}
= \bigl(\text{kg}\,\text{m}^{-3}\bigr)^{1-\sum_j n_{i,j}}\,\text{s}^{-1}.
$$

For a first-order single reactant ($\sum_j n_{i,j}=1$) this reduces to s⁻¹; for a
second-order/bimolecular case it becomes m³·kg⁻¹·s⁻¹. The struct docstring records
this convention ("[s⁻¹] or [m³/(kg·s)]", src/Materials/reactions.jl) but the
code does **not** check it — the user is responsible for supplying $A_i$ in units
consistent with the chosen orders.

---

## 6.5 Stoichiometry and Species Source Terms

### 6.5.1 The species source

Summing the per-reaction mass conversion over all reactions gives the volumetric
chemical source for component $j$ (kg·m⁻³·s⁻¹):

$$
S_j(T,\boldsymbol{\xi}) \;=\;
\left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}
\;=\; \sum_{i=1}^{N_R} \nu_{i,j}\, r_i(T,\boldsymbol{\xi}),
\tag{6.5}
$$

where $\nu_{i,j}$ is the stoichiometric (mass) coefficient: **negative** when $j$
is a reactant of reaction $i$, **positive** when $j$ is a product. This $S_j$ is
the chemical contribution to the species conservation equation of §3, Governing
Equations; for gases it adds to the transport divergence, for solids/liquids it is
the entire right-hand side (§6.5.4).

### 6.5.2 The single-reactant mass convention

The current model supports exactly one reactant per reaction. Per unit reaction
extent, that reactant is consumed with **unit mass stoichiometry** (an implicit
$\nu = -1$), and the products are produced with their explicit per-product mass
yields. The kernel that applies one reaction to the source vector is:

```julia
@inline function _apply_reaction_source!(src, rxn::Reaction{NRin}, r) where {NRin}
    NRin == 1 || error("Multi-reactant reaction $(rxn.name) (NRin=$NRin) ...")
    @inbounds src[rxn.reactants[1][1]] -= r          # reactant: ν = -1
    @inbounds for (idx, ν) in rxn.products
        src[idx] += ν * r                            # products: ν = yield ≥ 0
    end
    return nothing
end
```

(src/Physics/kinetics.jl). The `NRin == 1` guard **errors at runtime** if
a multi-reactant reaction is ever applied — multi-reactant support would require
extending the reactant tuple to a 3-ary `(idx, order, mass_stoich)` form so that a
mass stoichiometry (not just a rate order) could be attached to each reactant
(src/Materials/reactions.jl). Note the asymmetry baked into the model: a
reactant carries a rate **order** $n_{i,j}$ (how its concentration scales the
rate) but its **mass** stoichiometry is fixed at $1$; a product carries a mass
**yield** but no rate order (products never appear in $\Phi_i$).

For a single-reactant reaction, then, the coefficients are

$$
\nu_{i,\text{reactant}} = -1, \qquad
\nu_{i,\text{product }p} = y_{i,p} \ge 0,
$$

and (6.5) becomes, for reaction $i$ acting on its reactant $a$ and products $p$:

$$
S_a \mathrel{-}= r_i, \qquad S_p \mathrel{+}= y_{i,p}\, r_i.
$$

### 6.5.3 Mass conservation as a build-time invariant

Total mass is conserved by reaction $i$ if and only if the consumed reactant mass
equals the produced product mass, i.e.

$$
\sum_{p\,\in\,\mathrm{products}(i)} y_{i,p} = 1.
\tag{6.6}
$$

This is enforced **at construction time**, not at run time — in the `Reaction`
inner constructor itself, so no construction path bypasses it.
`_validate_stoichiometry` (src/Materials/reactions.jl) sums the product yields
and requires $\sum_p y_{i,p} = 1$ against the implicit unit reactant mass of
the single-reactant convention. Each yield is additionally required to be
non-negative (`_validate_yield`, enforced even with `validate_mass = false` —
a pair like $(1.5, -0.5)$ sums to 1 but would silently drive the
negative-yield species' concentration negative without bound), and each
reaction order must satisfy $n \ge 0$ (`_validate_reaction_order`). The
tolerance is relative to the total coefficient mass with an absolute floor:

$$
\text{tol} = \max\!\bigl(10^{-10}\,(m_{\mathrm{react}} + m_{\mathrm{prod}}),\ 10^{-14}\bigr),
$$

(src/Materials/reactions.jl). A violation raises an error at construction unless
the user passes `validate_mass = false`, which disables only the sum check.

Worked example (the canonical charring reaction):

```julia
Reaction(:charring, 1 => (2, 3); A=1e15, E=200e3, h=500e3, yields=(0.2, 0.8))
```

assigns $\nu_{i,1} = -1$ (virgin consumed), $\nu_{i,2} = 0.2$ (char produced),
$\nu_{i,3} = 0.8$ (gas produced); the check $0.2 + 0.8 = 1.0$ passes. The species
sources at rate $r_i$ are then $S_1 = -r_i$, $S_2 = +0.2\,r_i$, $S_3 = +0.8\,r_i$,
whose sum is exactly $0$ — mass is conserved instantaneously.

Because conservation is guaranteed by the build-time check and the implicit
$\nu=-1$ convention, **no run-time projection or Lagrange multiplier is used**; the
species sum $\sum_j S_j$ is identically zero up to floating-point roundoff. The
validation utilities of §6.8 exist purely as diagnostics.

### 6.5.4 Gas-vs-solid dispatch

The species source $S_j$ is the same expression for every phase; the phase
distinction lives in the *conservation equation* into which $S_j$ is inserted, not
in the kinetics. At final RHS assembly (src/Discretization/finite_volume.jl)
the code calls `apply_species_residual_unrolled!`, which assembles

$$
\frac{d\xi_j}{dt} =
\begin{cases}
-(\nabla\!\cdot\! \mathbf J_j)_{\mathrm{cell}} + S_j, & j\ \text{gaseous},\\[4pt]
S_j, & j\ \text{solid or liquid},
\end{cases}
\tag{6.7}
$$

so solids/liquids have **no transport term** — they are Lagrangian, fixed to the
material point (in the ALE frame, §11) — while gases carry the divergence of their
mass flux $\mathbf J_j$ (§9). The dispatch is resolved at compile time from the
component's static phase type (`is_gas`, src/Materials/components.jl), with no
run-time branching. The full energy/species RHS assembly sequence is documented in
§7, Heat Conduction, Convection & Energy Assembly.

> **Implementation note.** Species source: `species_source_terms!` and the
> `@generated` kernel `_species_source_terms_gen!`, src/Physics/kinetics.jl;
> per-reaction application `_apply_reaction_source!` (same file).
> Build-time balance: `_validate_stoichiometry`, src/Materials/reactions.jl.

---

## 6.6 Heat of Reaction

### 6.6.1 The reaction heat source and its sign convention

The volumetric heat liberated or absorbed by chemistry is

$$
Q_{\mathrm{rxn}}(T,\boldsymbol{\xi}) \;=\; -\sum_{i=1}^{N_R} h_i(T,\boldsymbol{\xi})\, r_i(T,\boldsymbol{\xi}),
\tag{6.8}
$$

(W·m⁻³), implemented as

```julia
@inline function reaction_heat_source(rates, heats)
    Tu = promote_type(eltype(rates), eltype(heats))
    Q = zero(Tu)
    @inbounds for i in eachindex(rates, heats)
        Q -= heats[i] * rates[i]          # Note the negative sign!
    end
    return Q
end
```

(src/Physics/kinetics.jl). The single minus sign is the entire content of
the sign convention, and it is worth stating precisely because it is the most
error-prone point in the whole model.

The field `rxn.heat` stores $h_i$ as an **enthalpy-of-reaction**: the heat
*absorbed* per kilogram of (first) reactant consumed. The storage convention is
(overload note H1):

$$
\boxed{\;h_i > 0 \Rightarrow \text{endothermic (absorbs heat, cools)};\quad
h_i < 0 \Rightarrow \text{exothermic (releases heat, heats)}.\;}
$$

The minus sign in (6.8) converts this enthalpy convention into the **heat-source**
convention required by the energy equation, where a positive source heats the
material:

- **Endothermic** ($h_i > 0$, $r_i > 0$): $Q_{\mathrm{rxn}} = -h_i r_i < 0$ — the
  cell is cooled, as a decomposition that consumes heat should.
- **Exothermic** ($h_i < 0$, $r_i > 0$): $Q_{\mathrm{rxn}} = -h_i r_i > 0$ — the
  cell is heated, as char oxidation should.

The source comment states this explicitly (src/Physics/kinetics.jl). In
the energy equation $Q_{\mathrm{rxn}}$ enters additively on the right-hand side
together with the conduction divergence, gas-advection source, and radiation
source, and the whole is divided by the mixture volumetric heat capacity
$\rho c_p^{\mathrm{eff}}$ (every component; see §3, §5, §7):

$$
\frac{dT_{\mathrm{cell}}}{dt} =
\frac{-\,(\nabla\!\cdot\! q_{\mathrm{cond}}) + S_{\mathrm{conv}}
+ Q_{\mathrm{rad}} + Q_{\mathrm{rxn}}}{\rho c_p^{\mathrm{eff}}}
\tag{6.9}
$$

(`compute_rhs!`, Step 8, src/Discretization/finite_volume.jl).

**Units check.** $[h_i] = $ J·kg⁻¹ (per kg of first reactant) and $[r_i] = $
kg·m⁻³·s⁻¹, so $[h_i r_i] = $ J·m⁻³·s⁻¹ $=$ W·m⁻³, the volumetric power density
demanded by (6.9). The "per kg of first reactant" basis is consistent with the
single-reactant mass convention of §6.5: the reactant is consumed at mass rate
$r_i$, and $h_i$ is the enthalpy per unit of that mass.

### 6.6.2 Temperature- and state-dependent heat

$h_i$ is not restricted to a constant. The `heat::Ph` field is typed by any
`AbstractPropertyFunction` (src/Materials/reactions.jl), and in
`evaluate_reactions!` the heat is evaluated with **both** $T$ and the full
concentration tuple $\boldsymbol{\xi}$:

```julia
heats[$i] = rxn.heat(T, ξ)   # Pass ξ for state-dependent heats
```

(src/Physics/kinetics.jl). The two-argument call dispatches as follows
(src/Materials/property_functions.jl):

- `ConstantProperty(h)` $\to h$ (the `(T,ξ)` method ignores both arguments);
- `LinearProperty(a,b)` $\to a + bT$ (affine in temperature, $\xi$ ignored);
- `PolynomialProperty` $\to \mathrm{evalpoly}(T,\text{coeffs})$;
- `StateDependentProperty(f)` $\to f(T, \boldsymbol{\xi})$ — the only type that
  actually consumes $\boldsymbol{\xi}$.

A constant heat is created by passing a number, an affine heat by passing a
`(a,b)` tuple, through `to_property` (src/Materials/property_functions.jl):

```julia
Reaction(:decomp, 1 => 2; A=8.5e12, E=188e3, h=870e3)          # ConstantProperty
Reaction(:decomp, 1 => 2; A=8.5e12, E=188e3, h=(870e3, 100.0)) # LinearProperty
```

(Positive `h` = endothermic, the storage convention of §6.6.1; PMMA-style
depolymerization absorbs ≈870 kJ·kg⁻¹.)

The principal use of `StateDependentProperty` is the **heat of moisture sorption**
in hygroscopic materials, where the differential heat depends on the current
moisture content $MC$. A representative closure (the form used in the moisture
examples) is

$$
h_{\mathrm{sorp}}(T,\boldsymbol{\xi}) =
\Bigl(2275 + 860\,e^{-0.120\,MC}\Bigr)\times 10^{3}\ \text{J·kg}^{-1},
\qquad MC = 100\,\frac{\xi_{\mathrm{moist}}}{\sum_{k}\xi_k},
$$

encoded as a closure over the moisture index (src/Materials/property_functions.jl).
Because the closure is passed $\boldsymbol{\xi}$ at evaluation, the
sorption heat self-consistently rises as the wood dries — a behaviour impossible
to represent with a temperature-only property. Calling a `StateDependentProperty`
with only $T$ raises an error (src/Materials/property_functions.jl),
guaranteeing the state argument is never silently dropped.

Whether the reaction heat is "separate from" the sensible-enthalpy transport of
the released gases is precisely the **datum question** of §6.6.3 — the partition
between "chemical" and "sensible" enthalpy depends on the temperature at which
$h_i$ is referenced.

> **Implementation note.** Heat source: `reaction_heat_source`,
> src/Physics/kinetics.jl. Heat evaluation with state:
> `_evaluate_reactions_gen!` (same file). Property dispatch:
> src/Materials/property_functions.jl.

### 6.6.3 The enthalpy datum: heats at reaction temperature vs the $T_{\text{ref}}$ datum

A heat of reaction is only defined relative to a **datum**: $h_i$ may be quoted

- **at the local reaction temperature** — $h_i(T) = \sum_j \nu_{ji}\,h_j(T)$,
  the heat a DSC measures at the decomposition temperature, and the usual
  meaning of inverse-calibrated (ThermaKin/Gpyro-style) effective heats; or
- **at the fixed ambient datum** $T_{\text{ref}} = 298.15$ K —
  $h_i = \sum_j \nu_{ji}\,h_j(T_{\text{ref}})$, the formation-enthalpy style.

The energy equation (§3.3) uses the first. With heats at the local reaction
temperature, subtracting the species equations from the total-enthalpy balance
gives the **exact** temperature equation
$\rho c_p\,\partial_t T = -\nabla\cdot q + S_{\text{conv}} + Q_{\text{rad}} +
Q_{\text{rxn}}$ — no further term exists, because the gas products' sensible
enthalpy $\Delta h_g(T, T_{\text{ref}}) = \int_{T_{\text{ref}}}^{T} c_{p,g}\,\mathrm{d}T'$ is already *inside* $h_i(T)$.

A heat quoted at the $T_{\text{ref}}$ datum is converted by Kirchhoff's law,
$h_i(T) = h_i(T_{\text{ref}}) + \sum_j \nu_{ji}\,\Delta h_j(T, T_{\text{ref}})$,
with the sum running over the reactant ($\nu = -1$) and the products; the
converted heat can be supplied as a temperature-dependent heat property (§6.6.2).
The conversion is the integral of the heat-capacity change across the reaction,
$\int_{T_{\text{ref}}}^{T}\bigl(\sum_p y_{i,p}c_{p,p} - c_{p,\text{reactant}}\bigr)\,\mathrm{d}T'$.
The gas products' share alone, $\Delta h_g(900\,\text{K}, 298.15\,\text{K}) =
0.66$–$1.20$ MJ·kg⁻¹ for $c_{p,g} = 1.1$–$2.0$ kJ·kg⁻¹·K⁻¹, is largely offset by
the reactant's, but a heat-capacity change of 0.5 kJ·kg⁻¹·K⁻¹ still shifts the heat
by 0.3 MJ·kg⁻¹ at 900 K, comparable to typical wood/polymer pyrolysis endotherms
(0.1–0.9 MJ·kg⁻¹). A $T_{\text{ref}}$-datum heat used unconverted mis-states the
reaction endotherm by that integral.

---

## 6.7 Allocation-Free Evaluation

The reaction model sits in the innermost loop of every residual and Jacobian
evaluation, so it is engineered to allocate nothing and to specialise completely
on the (compile-time-known) reaction and component counts $N_R$ and $N_C$. Two
`@generated` kernels and one inline reduction achieve this.

**Rate and heat evaluation.** `evaluate_reactions!` forwards to
`_evaluate_reactions_gen!` (src/Physics/kinetics.jl), which the compiler
expands into a flat block — one explicit `rates[i] = reaction_rate(...)` and
`heats[i] = rxn.heat(T, ξ)` statement per reaction — with no run-time loop over
the heterogeneous `reactions` tuple. Heterogeneity (each `Reaction` may have a
different `NRin`, `NRout`, and heat property type) is therefore handled by static
dispatch with no `Union` boxing and no allocation.

**Species sources.** `species_source_terms!(sources, rates, material)` zeroes the
output and calls `_species_source_terms_gen!`, again unrolled to one
`_apply_reaction_source!` per reaction. It is the single method, called from
`compute_rhs!` and from the Jacobian operators that need the reaction sources.

**Heat sum.** `reaction_heat_source` is a single `@inline` reduction over the
pre-filled `rates`/`heats` vectors (src/Physics/kinetics.jl), promoting the
accumulator type so Duals carry through.

All three use the promoted scalar type so that, in a Dual-typed (AD) pass, the
caches hold Duals and derivatives propagate end-to-end. The per-cell scratch
(`rates`, `heats`, `dξ_rxn`, `Q_rxn`) lives in the RHS cache and is sized by
$N_R$, $N_C$, and the cell count (`RHSCache`,
src/Discretization/finite_volume.jl), so the hot loop reuses storage
across cells and time steps.

> **Implementation note.** Generated kernels: `_evaluate_reactions_gen!` and
> `_species_source_terms_gen!` (src/Physics/kinetics.jl). Compile-time constants $N_C, N_R$ are the
> `Material{NC,NR}` type parameters (src/Materials/materials.jl).

---

## 6.8 Validation Utilities

Stoichiometry errors are caught when a reaction is built (§6.5.3); two small
diagnostics check a reaction afterwards.

**`verify_stoichiometry(rxn; tol=1e-10)`** computes the mass-balance error of a
single reaction via `mass_balance_error`:

$$
\varepsilon_{\mathrm{mb}}(i) = \Bigl|\,1 - \textstyle\sum_p y_{i,p}\,\Bigr|.
$$

For a single-reactant reaction it compares the implicit unit reactant mass against
the product-yield sum; for $N_{\mathrm{in}} > 1$ it returns `NaN` (unsupported).
The function warns and returns `false` if $\varepsilon_{\mathrm{mb}} > $ `tol`.
`validate_problem`, which `solve` calls unless `validate = false`, checks that the
material has at least one component and that every reactant and product index
lies in $1..N_C$ (numeric reactions are not index-checked when the `Material` is
built); with `strict = true` it also runs `verify_stoichiometry` over every
reaction of the material.

Symbolic reactions are resolved to integer indices by `resolve_reaction` against
the component name→index map when the `Material` is built, so a reaction that
names a non-existent component errors with a clear message before any solve.

> **Implementation note.** `verify_stoichiometry`, `mass_balance_error`,
> `resolve_reaction`: src/Materials/reactions.jl. `validate_problem`:
> src/Problem/problem_validation.jl.

---

## 6.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

All three reference codes solve essentially the same condensed-phase chemistry —
mass-concentration Arrhenius kinetics with mass-yield stoichiometry — and
Pyrolysis.jl deliberately reproduces their physics. The differences are in
generality, in the gating mechanism, and most consequentially in the **sign and
basis of the heat of reaction**. We map each code's notation to ours on first use.

### 6.9.1 ThermaKin2Ds

ThermaKin (Stoliarov & Lyon, *ThermaKin2Ds Reference v2*, Eq. 5) writes a reaction
as a generalised two-reactant, two-product step

$$
\theta_1\,\text{COMP1} + \theta_2\,\text{COMP2} \to
\theta_3\,\text{COMP3} + \theta_4\,\text{COMP4} + h,
$$

with stoichiometric coefficients $\theta_i^{\,j}$ (their notation) which map to our
$\nu_{i,j}$ (overload note A2 — $\theta$ is reserved for strain rate in our
nomenclature). The rate is **explicitly bimolecular** (ThermaKin Eq. 6):

$$
r = A\,\exp\!\Bigl(-\frac{E}{RT}\Bigr)\,\xi_{\mathrm{COMP1}}\,\xi_{\mathrm{COMP2}},
\qquad \xi = mf\cdot\rho,
$$

with the convention that "in the absence of a second reactant, $\xi_{\mathrm{COMP2}}$
is set to 1." Pyrolysis.jl differs in two ways. First, our rate is **n-ary in the
order but single-reactant in the mass**: we allow an arbitrary order $n_{i,j}$ on
the single supported reactant (a strict generalisation of ThermaKin's
unit-exponent factor for COMP1), but we do not yet support a genuine second
mass-bearing reactant (§6.10). ThermaKin's COMP2 mechanism is precisely what our
multi-reactant extension would need. Second, ThermaKin's temperature limit is a
**hard cutoff** ("if temperature decreases below the lower limit or increases above
the upper limit, the rate of reaction is set to 0"), whereas we use a smooth tanh
ramp (§6.2) for AD/Newton compatibility. ThermaKin's species/heat assembly
($\partial\xi_j/\partial t = \sum_i \theta_i^{\,j} r_i + \text{transport}$, their
Eqs. 15–16, and heat $= h\,r$) is identical in form to our (6.5) and (6.8).

The pivotal difference is the **heat-of-reaction sign**. ThermaKin states
explicitly: "Following a convention frequently adopted in pyrolysis studies, a
positive $h$ is used to represent an **exothermic** reaction." This is the
**opposite** of our storage convention ($h>0$ endothermic). The two are reconciled
by the minus sign in our $Q_{\mathrm{rxn}} = -\sum_i h_i r_i$: a ThermaKin input
of $h_{\mathrm{TK}} = +870$ kJ·kg⁻¹ (their exothermic) corresponds to our
$h = -870$ kJ·kg⁻¹, and both yield a positive (heating) source. A user porting a
ThermaKin material **must flip the sign of every heat of reaction.** (ThermaKin
also notes its reaction parameters are not constrained by mass or energy
conservation — the user owns that; Pyrolysis.jl instead enforces mass balance at
construction, §6.5.3.)

### 6.9.2 Gpyro

Gpyro (Lautenberger, *Gpyro Technical Reference*) tracks a **conversion** variable
$\alpha_{A_k}$ for condensed-phase reactant $A_k$ and writes the $n^{\text{th}}$
-order heterogeneous rate (their Eq. 39, with $Z_k \to A_i$, $E_k\to E_i$,
$n_k\to n_{i,j}$, and the Gpyro volumetric rate $\dot\omega''' \to r$):

$$
r_k = f(\alpha_{A_k})\,
\frac{(\bar\rho\,Y_{A_k}\Delta z)_\Sigma}{\Delta z}\,
Z_k\,e^{-E_k/RT}\,g(Y_{\mathrm{O_2}})\,H(T - T_{\mathrm{crit},k}),
\qquad f(\alpha_{A_k}) = (1-\alpha_{A_k})^{n_k}.
$$

Two Gpyro features have no analogue in Pyrolysis.jl: the oxygen-sensitivity factor
$g(Y_{\mathrm{O_2}})$ (their Eq. 40 — our kinetics have no $Y_{\mathrm{O_2}}$
dependence, §6.10), and the rich library of specialised reaction models
$f(\alpha)$ (nucleation, diffusion, Avrami, etc., their Table 1 — we offer only the
$n^{\text{th}}$-order power law of (6.4)). Gpyro's critical temperature
$H(T-T_{\mathrm{crit},k})$ is again a **hard** Heaviside, where we use a tanh ramp.

The conversion formalism is worth dwelling on because it reveals when the two codes
agree. Gpyro normalises by the cumulative species mass
$(\bar\rho Y_{A_k}\Delta z)_\Sigma$ (the initial plus all subsequently formed mass
of $A_k$), so that each cell "behaves as if it was a TGA sample." When $n_k=1$ and
$n_{\mathrm{O_2},k}=0$, Gpyro's rate collapses (their Eq. 41) to

$$
r_k = \bar\rho\,Y_{A_k}\,Z_k\,e^{-E_k/RT} = \xi_{A_k}\,A_i\,e^{-E_i/RT},
$$

since $\bar\rho Y_{A_k} = \xi_{A_k}$ is exactly our mass concentration. **For
first-order, oxygen-insensitive kinetics with no specialised model, our (6.1) and
Gpyro's rate are identical** — the conversion bookkeeping reduces to plain mass
concentration. Gpyro's stoichiometry (solid fraction $\mathrm{SF}_k$ and the
gas-yield split, their Eqs. 23–24) is a parameterisation of the same mass-balanced
yields we store directly as $\nu_{i,j}$; Gpyro's $\chi_k$ swelling parameter
belongs to the volume-change model (§10), not the kinetics.

On the **heat of reaction**, Gpyro (Eq. 43) splits the enthalpy into a solid part
$\Delta H_{\mathrm{sol},k}$ and a volatilisation part $\Delta H_{\mathrm{vol},k}$
and forms $\dot Q_{s,k}''' = -\dot\omega_{fB,k}'''\Delta H_{\mathrm{sol},k} -
\dot\omega_{fg,k}'''\Delta H_{\mathrm{vol},k}$. Setting
$\Delta H_{\mathrm{sol},k}=\Delta H_{\mathrm{vol},k}=\Delta H_k$ recovers (their
Eq. 45) the per-reactant form $\dot Q_{s,k}''' = -\dot\omega_{dA_k}'''\Delta H_k$,
which is **exactly our** (6.8) with $\Delta H_k \equiv h_i$ and
$\dot\omega_{dA_k}''' \equiv r_i$. Gpyro's sign convention agrees with ours:
"a positive value of $\Delta H_{\mathrm{vol},k}$ denotes an **endothermic**
reaction, and a negative value … an exothermic reaction." Our single heat $h_i$ is
thus Gpyro's per-reactant-consumed heat (the $\Delta H_{\mathrm{sol}} =
\Delta H_{\mathrm{vol}}$ special case), which Gpyro itself documents as the way to
specify a heat "per kg of reactant consumed" — precisely our per-first-reactant
basis (§6.6.1). Gpyro additionally separates sensible-enthalpy transport into the
$\sum(\dot\omega_f - \dot\omega_d)h_i$ term of its gas energy equation; in
Pyrolysis.jl that bookkeeping is the gas-advection source $S_{\text{conv}}$ of
§7 together with the reaction-temperature datum of $h_i$ (§6.6.3).

### 6.9.3 FDS

The FDS solid-phase pyrolysis model (FDS *Technical Reference*, Eqs. 7.34–7.35)
maps most closely term-for-term. With $\alpha,\beta$ the FDS material/reaction
indices and $\rho_{s,\alpha}\to\xi_j$, $A_{\alpha\beta}\to A_i$,
$E_{\alpha\beta}\to E_i$, FDS writes

$$
r_{\alpha\beta} = \rho_{s,\alpha}^{\,n_{s,\alpha\beta}}\,
A_{\alpha\beta}\,e^{-E_{\alpha\beta}/RT_s}\,
[X_{\mathrm{O_2}}(x)]^{\,n_{\mathrm{O_2},\alpha\beta}}\,
T_s^{\,n_{t,\alpha\beta}},
$$

a product of a reactant-concentration power, the Arrhenius factor, an oxygen
power, **and a temperature power** $T_s^{n_t}$. Pyrolysis.jl reproduces the first
two factors (with the same mass-concentration basis — FDS notes $A_{\alpha\beta}$
has units 1/s for $n_s=1$, exactly our convention) and omits the oxygen and
explicit-temperature powers (§6.10). FDS's species production
$S_\alpha = \sum v_{\alpha,\alpha'\beta} r_{\alpha'\beta}$ (Eq. 7.37) and gas
production $\dot m_\gamma''' = \sum v_{\gamma,\alpha\beta} r_{\alpha\beta}$
(Eq. 7.38) with yields $v$ summing to 1 are our (6.5)–(6.6). FDS's chemical heat
source $\dot q_{s,c}''' = -\sum r_{\alpha\beta} H_{r,\alpha\beta}$ (Eq. 7.42) is
our (6.8) with $H_{r,\alpha\beta}\to h_i$, and FDS uses the **same** sign as Gpyro
and Pyrolysis.jl ($H_r > 0$ endothermic). FDS interprets each in-depth cell as a
virtual TGA sample, matching our cell-as-batch-reactor view.

As Gpyro notes, FDS Version 5's kinetics are the special case of Gpyro's Eq. 42
with $n_{\mathrm{O_2}}=0$, conventional reaction order, and no volume change — and
therefore also the special case of our (6.1) for first/$n$-th-order
oxygen-insensitive kinetics.

**Kinetic-parameter recovery from TGA.** A practical FDS capability worth
reproducing in user workflows is the closed-form estimate of $A_i, E_i$ from a TGA
peak. Given a reference (peak) temperature $T_p$, the normalised peak mass-loss
rate $r_p/Y_s(0)$ (s⁻¹), a constant heating rate $\dot T$ (K·s⁻¹), and residue
yield $v_s$, FDS uses (FDS *User Guide*, Eqs. 9.4–9.5; the Lyon–Stoliarov
relations)

$$
E_i = \frac{e\,r_p}{Y_s(0)}\,\frac{R\,T_p^2}{\dot T},
\qquad
A_i = \frac{e\,r_p}{Y_s(0)}\,e^{E_i/(R\,T_p)},
\qquad
\frac{r_p}{Y_s(0)} = \frac{2\dot T}{\Delta T}(1 - v_s),
$$

with $e$ Euler's number and $\Delta T$ the (triangular) pyrolysis-range width
(default 80 K), $\dot T$ default 5 K·min⁻¹. These produce a first-order
$(A_i, E_i)$ pair that, fed into our `Reaction(...)` constructor, reproduces the
TGA peak — a recommended starting point for cone/TGA calibration (see §18,
Verification & Validation, and the User Guide).

### 6.9.4 Summary of correspondences

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| Rate basis | mass conc. $\xi_j$ | mass conc. $\xi$ | $\bar\rho Y_i$ (= $\xi$) | $\rho_{s,\alpha}$ (= $\xi$) |
| Reactant factor | $\xi_j^{n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ (1 reactant) | $\xi_{\mathrm{C1}}\xi_{\mathrm{C2}}$ (bimolec.) | $(1-\alpha)^{n_k}$ + models | $\rho_s^{n_s}$ |
| Extra factors | — | — | $g(Y_{\mathrm{O_2}})$ | $X_{\mathrm{O_2}}^{n_{\mathrm{O_2}}} T_s^{n_t}$ |
| Temp. gate | smooth tanh | hard cutoff | hard $H(T-T_{\mathrm{crit}})$ | (none / power) |
| Stoich. yields | $\nu_{i,j}$, $\sum y=1$ | $\theta_i^j$ | $\mathrm{SF}_k$ / yields | $v$, $\sum v=1$ |
| Heat sign ($h>0$) | **endothermic** | **exothermic** | endothermic | endothermic |
| Heat source | $-\sum h_i r_i$ | $+h\,r$ | $-\sum \Delta H_k r_k$ | $-\sum H_r r$ |

The lone sign reversal (ThermaKin's $h>0$ exothermic) is the single most important
porting hazard; Gpyro and FDS agree with us.

---

## 6.10 Limitations and Open Issues

The kinetics model is deliberately the simplest one that covers the target
applications (cone calorimetry, TGA/DSC/MCC idealisation, charring/non-charring/
intumescent solids, moisture). Its boundaries, with the source-of-truth guards:

1. **Single reactant only.** `_apply_reaction_source!` (src/Physics/kinetics.jl)
   errors for $N_{\mathrm{in}} > 1$, and `mass_balance_error` returns `NaN`
   for multi-reactant reactions (src/Materials/reactions.jl). True bimolecular
   chemistry (ThermaKin's COMP1·COMP2, char oxidation consuming a gaseous O₂
   reactant) is not representable. Extending it requires a 3-ary reactant tuple
   `(idx, order, mass_stoich)` so a mass stoichiometry can be attached to each
   reactant; the outer struct shape already anticipates this
   (src/Materials/reactions.jl). Multi-step schemes are instead built as
   **sequences** of single-reactant reactions.

2. **No oxygen or pressure dependence.** The rate (6.1) depends only on $T$ and
   condensed-phase mass concentrations. There is no $g(Y_{\mathrm{O_2}})$ (Gpyro)
   or $X_{\mathrm{O_2}}^{n_{\mathrm{O_2}}}$ (FDS) factor and no partial-pressure
   term, so oxidative pyrolysis and pressure-enhanced decomposition cannot be
   modelled directly. Pressure effects enter only through the separate gas
   transport/permeability block (§9), never the rate constant.

3. **Constant pre-exponential and activation energy.** $A_i$ and $E_i$ are scalar
   `Float64` fields (src/Materials/reactions.jl). There is no
   temperature-dependent $A(T)$, no explicit $T^{n_t}$ power (FDS's fourth factor),
   no kinetic compensation rule $A(E)=A_0 e^{\beta E}$, and no isokinetic shift. A
   temperature-dependent rate constant would require wrapping the entire constant
   in a `CallableProperty` and providing a custom evaluation path — not part of the
   public API.

4. **No distributed-activation-energy (DAEM) tooling.** A DAEM model can be
   approximated by stacking several parallel `Reaction` instances with different
   $(A_i, E_i)$, but there is no built-in quadrature over an activation-energy
   distribution or fitting helper for it.

5. **Depletion-limiter bias.** The limiter (6.3) under-predicts the rate near
   and below $\xi_{\mathrm{th}} = 1$ kg·m⁻³ (~24% at the threshold, growing as
   the reactant runs out) and converts exponential burnout tails into algebraic
   ones (§6.3.3). For materials whose physically meaningful residue
   concentrations are comparable to the default threshold, the threshold should
   be lowered (at some cost in stiffness) or the limiter disabled for
   validation.

These limitations are consistent across the kinetics and materials chapters; where a richer model is required (bimolecular oxidation, oxygen-sensitive
char burnout, DAEM), the present architecture's chief virtue is that the
allocation-free, AD-clean, compile-time-specialised kernels of §6.7 provide a
clean foundation onto which the additional factors can be grafted without
disturbing the surrounding finite-volume and Jacobian machinery (§13, §15).


---

# 7. Heat Conduction, Gas Convection, and Energy Assembly

This chapter derives the condensed-phase energy equation as Pyrolysis.jl solves it: the
non-conservative temperature equation, in which the volumetric thermal storage sums the
heat capacity of every component, pore gas included, while the sensible heat carried by the
gas flux relative to the skeleton is transported by an explicit volumetric advection source.
We develop each term from its continuous form to its finite-volume discretization — Fourier
conduction with distance-weighted harmonic-mean face conductivity, the discrete divergence
operator, the gas convective source $S_{\text{conv}}$, the reaction-heat source
$Q_{\text{rxn}}$, and the radiation source $Q_{\text{rad}}$ — and we trace the exact
eight-step assembly sequence in `compute_rhs!` that produces $dT_i/dt$ and
$d\xi_{j,i}/dt$. Particular care is taken with the three sign conventions of §2
(Nomenclature), with the precise set of terms the code includes and excludes, and with why
the surface energy balance carries no gas-enthalpy term alongside $S_{\text{conv}}$. Where
useful, the Pyrolysis.jl forms are mapped to and
contrasted with ThermaKin2Ds, Gpyro, and FDS.

Index conventions for this chapter follow §2 overload note I1: **$i$ is the cell index**
($1\le i\le n$) and **$r$ is the reaction index**; the component/species index is $j$
($1\le j\le N_C$). The spatial coordinate is $z$ with $z=0$ the bottom/substrate (boundary
id 2) and $z=L$ the top/exposed surface (boundary id 1); heat enters from the top. Positive
flux is transport in the $+z$ direction. The symbol $x$ appears in a few places (e.g.
$\partial T/\partial x$) only because the implementation uses it as a synonym for the
through-thickness coordinate $z$; the two are identical in this 1D model.

> **Frame of this chapter.** The Eulerian energy balance developed here has
> no solid-convection term — it assumes the condensed phase carrying
> $\rho c_p^{\text{eff}}$ is **stationary in the lab frame**. On a fixed
> (Eulerian) mesh that is exact. Under ALE the matrix *moves* (it subsides
> with the reaction-driven mesh), and the correct time derivative of the
> matrix-carried temperature is the *mesh-frame* one, with no additional
> advection term — see §11.2.2 for the relative-velocity argument.

---

## 7.1 Fourier Conduction and the Face Heat Flux

### 7.1.1 Continuous law

Conductive transport of sensible heat in the condensed phase obeys Fourier's law,

$$
\mathbf{q}_{\text{cond}} = -k(T,\xi)\,\nabla T ,
\qquad
q_{\text{cond}}(z) = -k(T,\xi)\,\frac{\partial T}{\partial z}
\quad [\text{W}\,\text{m}^{-2}],
$$

with the minus sign expressing the second law: heat flows down the temperature gradient, in
the direction of decreasing $T$. The effective conductivity $k=k_{\text{eff}}(T,\xi)$ is the
composition- and temperature-dependent mixture conductivity built from the component
conductivities $k_j(T)$ by the chosen mixing rule (PARALLEL / SERIES / WEIGHTED / BRUGGEMAN)
on the **solid** volume-fraction basis; its construction is the subject of §5 (Effective
Properties & Mixing Rules), and only its *face value* enters here.

### 7.1.2 Discrete face flux: distance-weighted harmonic mean

In the cell-centered finite-volume layout (§13, Finite-Volume Discretization) the unknown
$T$ is stored at cell centers $z_i^c$, and fluxes are required at the faces that separate
adjacent cells. For an interior face $f$ separating a left cell $L$ and a right cell $R$,
with cell-center temperatures $T_L,T_R$ and center-to-face distances $d_L,d_R$, the
conductive flux is

$$
\boxed{\;
q_f \;=\; -\,k_f\,\frac{T_R-T_L}{d_L+d_R}
\;}
\tag{7.1}
$$

where the **face conductivity** $k_f$ is the **distance-weighted harmonic mean** of the two
cell-center conductivities:

$$
\boxed{\;
k_f \;=\; \frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
\;}
\tag{7.2}
$$

The denominator $d_L+d_R$ is exactly the center-to-center distance $d_{LR}=|z_R^c-z_L^c|$, so
on a mesh where each cell center sits midway between its faces, $d_L+d_R=d_{LR}$ and (7.1)
reduces to $q_f=-k_f(T_R-T_L)/d_{LR}$.

**Derivation of the harmonic mean from flux continuity.** Treat the half-cell on either side
of the face as a thermal resistance. With cross-sectional area $A$ the resistance of the left
half is $R_L=d_L/(k_L A)$ and of the right half $R_R=d_R/(k_R A)$. At steady state the same
flux $q_f$ passes through both half-cells in series, so the temperature drop from $T_L$ to
$T_R$ across the combined resistance is

$$
T_L-T_R \;=\; q_f\,A\,(R_L+R_R)
\;=\; q_f\left(\frac{d_L}{k_L}+\frac{d_R}{k_R}\right).
$$

Solving for $q_f$ and writing it as an effective conductivity $k_f$ acting over the total
distance $d_L+d_R$,

$$
q_f \;=\; -\,\frac{T_R-T_L}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
\;\equiv\; -\,k_f\,\frac{T_R-T_L}{d_L+d_R},
\qquad
k_f=\frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}},
$$

which is (7.1)–(7.2). The harmonic (resistance-in-series) mean is the *only* face averaging
that makes the discrete flux exact for a piecewise-constant conductivity field (a sharp
material interface at the face): an arithmetic mean would over-predict the flux through the
more-resistive layer. This is the same series-resistance averaging applied to the
gas-transfer coefficient $\lambda$ and the permeability $\kappa$ (§9, Gas-Phase Mass
Transport).

**Uniform-mesh limit.** When the cell centers are equidistant from the face, $d_L=d_R=d$,
(7.2) collapses to the standard (unweighted) harmonic mean,

$$
k_f \;\xrightarrow{\,d_L=d_R\,}\; \frac{2d}{\dfrac{d}{k_L}+\dfrac{d}{k_R}}
\;=\; \frac{2k_Lk_R}{k_L+k_R}.
$$

On **stretched** meshes — for example, boundary-refined grids used to resolve a steep
near-surface thermal gradient — $d_L\ne d_R$, and the distance weighting corrects the
averaging. Empirically the distance-weighted form recovers $1$–$5\%$ accuracy in the
near-boundary fluxes relative to the unweighted harmonic mean on such meshes, at no extra
cost.

**Small-value guard.** The implementation of (7.2) returns $0$ if either $k_L<10^{-20}$ or
$k_R<10^{-20}\;\text{W}\,\text{m}^{-1}\text{K}^{-1}$, preventing a division blow-up when a
cell has been driven to a vanishing conductivity (e.g. a near-fully-depleted cell). This is a
numerical safety net; physically realistic conductivities are many orders above the
threshold.

> **Implementation note.** The single-face flux is `conductive_flux(k_face, T_L, T_R, d_LR)`
> in `src/Physics/heat_conduction.jl`; it consumes the face conductivity $k_f$ that
> `update_properties!` pre-computes with `harmonic_mean_weighted`
> (`src/Properties/effective_properties.jl`), where the `<1e-20` guard lives.

### 7.1.3 Interior flux pass and the boundary-face zeroing convention

All interior face fluxes are computed in one pass before any boundary condition is applied,
by `compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`). For each face $f$ the
routine reads the left/right cell ids from the mesh connectivity; for an **interior** face it
reads $T_L,T_R$ from the state, the center-to-center distance $d_{LR}=|z_R^c-z_L^c|$ from the
(possibly ALE-moving) cell-center array, and the cached face conductivity $k_f$, then writes
$q_f$ via (7.1). For a **boundary** face — one whose left or right cell id is $0$ — it writes

$$
q_f \;=\; 0 .
\tag{7.3}
$$

The zeroing in (7.3) is deliberate, not a skip: the boundary face flux is **overwritten**
later (Step 4 of §7.6) by the surface energy-balance solver, which computes the true
conductive flux into the boundary cell from the solved surface temperature $T_s$ (see §12,
Boundary & Initial Conditions). Zeroing first guarantees that if the BC step did not touch a
face, that face contributes nothing rather than a stale value, and it keeps the interior
conduction operator and the boundary operator cleanly separated for Jacobian assembly (the
`ConductionOperator` fills only interior tridiagonal entries; the `BCThermalOperator` fills
the boundary self-coupling — see §15, Temporal Integration & Structured Jacobian).

The face conductivities $k_f$ themselves are pre-computed once per RHS evaluation in
`update_properties!` and stored in the `PropertyCache` field `k_face`, so the flux pass does
no mixing-rule work. The property cache is parameterized on the live scalar type $T_u$
(`Float64` in normal operation, `ForwardDiff.Dual` during automatic differentiation), so the
entire conduction computation carries forward-mode sensitivities without re-allocation.

---

## 7.2 The Discrete Divergence Operator

The conductive flux enters the energy balance through its divergence. Integrating
$\nabla\!\cdot\mathbf{q}_{\text{cond}}$ over cell $i$ and applying the divergence theorem
gives the net outflow through the two bounding faces; dividing by the cell volume $V_i$ yields
the cell-averaged divergence

$$
\boxed{\;
(\nabla\!\cdot\mathbf{q})_i \;=\; \frac{1}{V_i}\bigl(q_{f_R}-q_{f_L}\bigr)\,A
\;}
\tag{7.4}
$$

where $q_{f_L},q_{f_R}$ are the fluxes at the left and right faces of cell $i$, $A$ is the
face area (read from the state as `A_face[f_left]`, which may vary in time in ALE mode), and
$V_i$ is the cell volume.

**Sign convention (outflow-positive).** With the global rule "positive flux = transport in
$+z$," the right face flux $q_{f_R}$ is an *outflow* (it carries heat out of cell $i$ in the
$+z$ direction) and the left face flux $q_{f_L}$ is an *inflow*. Hence the $(q_{f_R}-q_{f_L})$
ordering makes $(\nabla\!\cdot\mathbf{q})_i>0$ when there is a *net loss* of heat from the
cell. Because the energy equation carries $-(\nabla\!\cdot\mathbf{q})_i$ (see §7.3), a positive
divergence (net outflow) correctly *reduces* $dT_i/dt$: a cell losing more heat than it gains
cools.

The same operator applies component-wise to the gas species fluxes,

$$
(\nabla\!\cdot\mathbf{J}_j)_i \;=\; \frac{1}{V_i}\bigl(J_{j,f_R}-J_{j,f_L}\bigr)\,A ,
\qquad j=1,\dots,N_C ,
\tag{7.5}
$$

with units $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ when $J$ is in $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$ — i.e.
$(\nabla\!\cdot\mathbf{J}_j)_i$ is the volumetric net mass-efflux rate of species $j$.

> **Implementation note.** `divergence_1d!` (scalar heat flux) and `divergence_species_1d!`
> (matrix species flux, `[N_C × n_cells]`) are in `src/Discretization/finite_volume.jl`. Both
> use the *current* state volumes and face areas, so in ALE mode they automatically pick up
> the geometric volume change; in Eulerian mode the volume is constant. The species variant
> uses the identical formula component-wise inside an inner `@simd` loop.

A constant flux field gives zero divergence everywhere, and a telescoping sum
$\sum_i (\nabla\!\cdot\mathbf{F})_i V_i = F_{\text{boundary}}^{R}-F_{\text{boundary}}^{L}$
recovers the global balance — the discrete conservation property exploited by the diagnostics
ledgers (§16, Conservation Diagnostics).

---

## 7.3 Effective Heat Capacity and the Energy Equation

### 7.3.1 The non-conservative temperature equation

Pyrolysis.jl solves the energy balance in **non-conservative (temperature-form)** shape:

$$
\boxed{\;
\rho c_{p}^{\text{eff}}(T,\xi)\,\frac{\partial T}{\partial t}
\;=\;
-\,\nabla\!\cdot\mathbf{q}_{\text{cond}}
\;+\; S_{\text{conv}}
\;+\; Q_{\text{rad}}
\;+\; Q_{\text{rxn}}
\;}
\tag{7.6}
$$

with, term by term,

| term | expression | role | § |
|---|---|---|---|
| $\mathbf{q}_{\text{cond}}$ | $-k\,\nabla T$ | Fourier conduction | 7.1 |
| $S_{\text{conv}}$ | $-\sum_g c_{p,g}(T)\,\mathbf{J}_g\!\cdot\!\nabla T$ | gas sensible-enthalpy advection | 7.4 |
| $Q_{\text{rad}}$ | model-dependent | volumetric radiation absorption | 7.5, §8 |
| $Q_{\text{rxn}}$ | $-\sum_r h_r r_r$ | reaction heat | 7.5 |

The subscript $g$ runs over the **gas** components only (the predicate `is_gas`). The storage
coefficient $\rho c_p^{\text{eff}}$ is the **mixture** volumetric heat capacity, defined
next.

### 7.3.2 Mixture heat capacity: every component stores sensible heat

The volumetric thermal storage sums every component, gas included:

$$
\boxed{\;
\rho c_p^{\text{eff}}(T,\xi)
\;=\!\!\sum_{j=1}^{N_C}\!\! \xi_j\,c_{p,j}(T)
\;}
\qquad [\text{J}\,\text{m}^{-3}\text{K}^{-1}]
\tag{7.7}
$$

The pore gas is at the matrix temperature (one-temperature model, §3.1.2) and its
concentration $\xi_g$ is an integrated state, so its storage term costs nothing and is
kept. Dropping it — the *quasi-steady-gas approximation* of Gpyro's default
thermal-equilibrium form (restored there by `FULL_QSG`) and the implicit convention of
FDS, which has no in-depth gas — would be an approximation with no benefit. At
$1\,\text{atm}$ the in-pore gas heat capacity is a few tenths of a percent of the matrix
heat capacity (the gas-to-solid density ratio is $\mathcal{O}(10^{-3})$ and the specific
heats are comparable), scaling with pressure; on white-pine cone decks it is at most
about half a percent and shifts the response by a uniform time lag of about a tenth
of a percent.

The gas enters the energy equation in two places and two places only: its storage in
(7.7) and its *flux relative to the skeleton* through the advective source
$S_{\text{conv}}$ (§7.4). The derivation of §3.3.2 shows these two terms, plus the
reaction heat with heats quoted at the reaction temperature, to be the complete
equation — ThermaKin Eq. 17 term for term. In ALE the whole $\rho c_p^{\text{eff}}$
rides the mesh with the skeleton, so there is no mesh-motion term in the energy
equation (§11.2.2).

> **Implementation note.** $\rho c_p^{\text{eff}}$ is `effective_heat_capacity` in
> `src/Properties/effective_properties.jl`, a `@generated` unrolled sum over the
> component tuple with no phase branch. The per-cell value used by the assembly is
> cached as `prop_cache.ρcp[i]`; the same function feeds the structured-Jacobian
> operators (through Dual evaluation) and the energy ledger; the park-and-mask
> temperature merge (`merge_temperatures_conservatively`, `src/Adaptivity/cell_merge.jl`)
> sums the same every-component heat capacities directly.

### 7.3.3 The discrete temperature update

Substituting the discrete divergence (7.4) into (7.6) and dividing by the cached mixture heat
capacity gives the cell temperature time derivative used by the integrator:

$$
\boxed{\;
\frac{dT_i}{dt}
\;=\;
\frac{1}{\rho c_{p,i}^{\text{eff}}}
\Bigl(
-\,(\nabla\!\cdot\mathbf{q})_i
\;+\; S_{\text{conv},i}
\;+\; Q_{\text{rad},i}
\;+\; Q_{\text{rxn},i}
\Bigr)
\;}
\tag{7.9}
$$

Equivalently, expanding the divergence,

$$
\frac{dT_i}{dt}
=
\frac{1}{V_i\,\rho c_{p,i}^{\text{eff}}}
\Bigl[
-\bigl(q_{f_R}-q_{f_L}\bigr)A
+ V_i\bigl(S_{\text{conv},i}+Q_{\text{rad},i}+Q_{\text{rxn},i}\bigr)
\Bigr].
$$

The three volumetric source terms are *cell-centered* (units $\text{W}\,\text{m}^{-3}$); only
the conduction term arrives as a flux divergence. The exact code is

```julia
ρcp = prop_cache.ρcp[i]
dT[i] =
    (
        -rhs_cache.div_q[i] +
        rhs_cache.S_conv[i] +
        rhs_cache.Q_rad[i] +
        rhs_cache.Q_rxn[i]
    ) / ρcp
```

so the storage coefficient enters at exactly one place — the divisor `ρcp`, the
mixture value (7.7).

### 7.3.4 Species equations (companion to the energy equation)

The same RHS pass assembles the species derivatives. With the gas-vs-solid dispatch (transport
for gases only — see §3, Governing Equations, and §9),

$$
\frac{d\xi_{j,i}}{dt}
=
\begin{cases}
-\,(\nabla\!\cdot\mathbf{J}_j)_i + S_{j,\text{rxn},i}, & \text{gas } j,\\[4pt]
S_{j,\text{rxn},i}, & \text{solid/liquid } j,
\end{cases}
\tag{7.10}
$$

where $S_{j,\text{rxn},i}=\sum_r \nu_{j,r}\,r_r$ is the stoichiometric reaction source (§6,
Reaction Kinetics). Solids and liquids have no transport flux ($\mathbf{J}_j\equiv 0$); they
move only through reactions (and, in ALE, with the mesh — §10, §11). This dispatch is
performed by the `@generated` `apply_species_residual_unrolled!` (`src/Materials/materials.jl`),
which emits `dξ[j,i] = -div_J[j,i] + dξ_rxn[j,i]` for gas components and
`dξ[j,i] = dξ_rxn[j,i]` otherwise, with the gas/solid branch resolved at compile time.

---

## 7.4 The Gas Convective (Advection) Source $S_{\text{conv}}$

### 7.4.1 Origin: the advection piece of the gas-enthalpy divergence

In a fully conservative energy equation, the sensible enthalpy carried by the gas appears as
the divergence of an enthalpy flux,

$$
\nabla\!\cdot\!\Bigl(\sum_g \Delta h_g(T,T_{\text{ref}})\,\mathbf{J}_g\Bigr),
\qquad
\Delta h_g(T,T_{\text{ref}}) = \int_{T_{\text{ref}}}^{T} c_{p,g}(\tau)\,d\tau .
$$

Expanding this divergence with the product rule (exact even for temperature-dependent
$c_{p,g}$, since $\nabla\Delta h_g = c_{p,g}(T)\,\nabla T$),

$$
\nabla\!\cdot\!\Bigl(\sum_g \Delta h_g(T,T_{\text{ref}})\,\mathbf{J}_g\Bigr)
=
\underbrace{\sum_g c_{p,g}\,\mathbf{J}_g\!\cdot\!\nabla T}_{\text{advection}}
\;+\;
\underbrace{\sum_g \Delta h_g(T,T_{\text{ref}})\,\nabla\!\cdot\!\mathbf{J}_g}_{\text{generation/divergence}} .
\tag{7.11}
$$

The energy equation keeps the **advection** piece as a volumetric source on the right-hand side of
the temperature equation, with a leading minus sign (because the conservative term sits on the
storage/transport side; moving it to the source side flips the sign):

$$
\boxed{\;
S_{\text{conv}} \;=\; -\sum_g c_{p,g}(T)\,\mathbf{J}_g\!\cdot\!\nabla T
\;}
\qquad [\text{W}\,\text{m}^{-3}]
\tag{7.12}
$$

The second (generation) piece of (7.11) needs no term of its own: with heats of reaction
quoted at the reaction temperature it is already part of $Q_{\text{rxn}}$ (§3.3.2, §6.6.3).
Physically, $S_{\text{conv}}$ is the rate at which flowing gas deposits or removes
sensible heat as it crosses a temperature gradient: gas moving from hot toward cold heats the
cell it passes through, and vice versa. It carries the enthalpy of the gas *flux* relative to
the skeleton; the gas *storage* sits in $\rho c_p^{\text{eff}}$ (§7.3.2).

### 7.4.2 Discretization

The cell-centered discrete form is

$$
\boxed{\;
S_{\text{conv},i}
\;=\;
-\sum_{g:\,\text{is\_gas}} c_{p,g}(T_i)\,\overline{J}_{g,i}\,
\left(\frac{\partial T}{\partial z}\right)_i
\;}
\tag{7.13}
$$

with two cell-centered ingredients:

**(a) Average gas flux $\overline{J}_{g,i}$.** The flux that "passes through" cell $i$ is the
arithmetic mean of the fluxes at its two bounding faces,

$$
\overline{J}_{g,i} \;=\; \tfrac12\bigl(J_{g,f_L}+J_{g,f_R}\bigr),
\tag{7.14}
$$

evaluated component-wise from the gas-flux matrix populated in Step 3 of the assembly.

**(b) Cell-centered temperature gradient $(\partial T/\partial z)_i$.** For an interior cell
the two one-sided differences are combined with **opposite-distance weights**
(`opposite_distance_gradient`, src/Geometry/geometric_quantities.jl),

$$
\left(\frac{\partial T}{\partial z}\right)_i
=
\frac{d_R\,\dfrac{T_i-T_L}{d_L} + d_L\,\dfrac{T_R-T_i}{d_R}}{d_L+d_R},
\tag{7.15}
$$

where $T_L,T_R$ are the neighbor cell-center temperatures and $d_L,d_R$ the respective
center-to-center distances. Each one-sided slope is a second-order estimate of $T'$ at its
segment midpoint ($z_i - d_L/2$, resp. $z_i + d_R/2$); weighting by the *opposite* distance
interpolates them linearly back to $z_i$, cancelling the leading $T''$ truncation term, so
(7.15) is second-order on arbitrary nonuniform meshes.

A **boundary cell** applies the same interpolation with the *surface* as the missing
neighbor: the solved surface temperature $T_s$ (Step 4 of the assembly) lives on the
boundary face at distance $\Delta z_{1/2} = V_i/(2A_f)$ from the center, giving the
face-side slope $s_{\rm face} = \pm(T_s - T_i)/\Delta z_{1/2}$ alongside the interior
one-sided slope $s_{\rm int}$, and

$$
\left(\frac{\partial T}{\partial z}\right)_i
= \frac{d_{\rm nb}\,s_{\rm face} + \Delta z_{1/2}\,s_{\rm int}}{d_{\rm nb}+\Delta z_{1/2}}
\tag{7.15b}
$$

— second order up to the boundary (exact on quadratics), and sensitive to the
surface-side gradient, where the outflowing gas is hottest, that a plain interior
one-sided difference misses. The implementation recovers the face slope AD-cleanly as
$s_{\rm face} = -q_{\rm total}[f]/k_i$, which is *algebraically* identical because the
thermal-BC pass writes the boundary flux as exactly
$q_{\rm total}[f] = \mp k_i (T_s - T_i)/\Delta z_{1/2}$ at the solved $T_s$ (the $k_i$
cancels). This is why $S_{\rm conv}$ is computed *after* the thermal BC pass. A fully
isolated cell (single-cell mesh, which does not occur in a well-formed problem) gives
$S_{\text{conv},i}=0$.

On a uniform mesh ($d_L=d_R=\Delta z$), (7.15) reduces to the standard central difference
$(T_R-T_L)/(2\Delta z)$. The unweighted mean of one-sided slopes and the wide-stencil form
$(T_R-T_L)/(z_R^c-z_L^c)$ both carry a first-order $\mathcal{O}(|d_R-d_L|)\,T''$ error on
stretched meshes and are not used.

> **Implementation note.** `compute_gas_convective_source!` in `src/Physics/convection.jl`
> builds $\overline{J}_{g,i}$ and the gradient per cell; the inner sum over gas components is
> the `@generated` `_gas_convective_source_unrolled`, which emits, for each gas $j$,
> `S -= cp_j * J_avg * dTdx` with `cp_j = comp.heat_capacity(T_i)` and `J_avg` from (7.14).
> Solid/liquid components are skipped by the compile-time `is_gas` branch. The accumulator is
> initialized to `zero(promote_type(typeof(T_i), typeof(dTdx), eltype(J_gas)))` so it is
> type-stable under AD.

### 7.4.3 The surface carries no gas-enthalpy term

$S_{\text{conv}}$ is always active and is the sole carrier of gas sensible-enthalpy
transport. It is also the reason the surface boundary condition is a *plain*
convective + radiative balance, $q_{\text{BC}}(T_s) = q_{\text{cond}}(T_s)$, with no
gas-enthalpy term: the gas energy leaving the domain has already been accounted for
volumetrically as it advected up to the surface.

Integrating (7.11) over the column shows the bookkeeping. With heats of reaction at the
reaction temperature, the generation piece of (7.11)
is inside $Q_{\text{rxn}}$, and the remaining advection piece, $\int_0^L S_{\text{conv}}\,dz$,
is the heating of the gas from its formation temperature to $T_s$ along its path. Under local
thermal equilibrium the gas crosses the surface at $T_s$ on both sides, so its enthalpy flux
is continuous there and the surface balance is conducted flux equals applied flux (§1.4.3).
A surface term $\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$ added to that balance would apply
the outflow a second time, on top of $\int S_{\text{conv}}\,dz$ — a spurious sink of order
$5$–$10\,\text{kW}\,\text{m}^{-2}$ that under-predicts the surface temperature and the
mass-loss rate. The convention is ThermaKin Eq. 17; Gpyro's thermal-equilibrium form uses the
same term. Parameter sets calibrated against a solver that applies a surface gas-enthalpy term
*in addition to* a volumetric advection source predict lower $T_s$ and MLR there than here;
re-validate imported calibrations (§18).

**The mass-flux/convection coupling is the blowing model of each `ConvectiveBC`.** The
*momentum-side* effect of the gas efflux — boundary-layer thickening that suppresses the
convective exchange coefficient, and the pooling of a dense gas on the face — is independent
of the enthalpy bookkeeping above and is the `blowing` model a `ConvectiveBC` carries
(§12.3.3): the coefficient becomes $h\,\chi$ with $\chi$ from the model and the outward gas
flux at the face. Because it modifies only the coefficient, it composes exactly with the
$S_{\text{conv}}$-only convention; no double-count arises. This is the only use of the gas
flux in the surface balance.

> **Implementation note.** The surface Newton solve is in
> `src/Discretization/boundary_temperature.jl`; it takes the outward gas flux at the face as
> its only gas-flux argument and passes it to the blowing models.

---

## 7.5 Reaction-Heat and Radiation Sources

### 7.5.1 Reaction heat $Q_{\text{rxn}}$ and its sign

The volumetric reaction-heat source is

$$
\boxed{\;
Q_{\text{rxn},i} \;=\; -\sum_r h_r\,r_r
\;}
\qquad [\text{W}\,\text{m}^{-3}]
\tag{7.16}
$$

where $r_r$ is the rate of reaction $r$ (units $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$; §6)
and $h_r$ is its heat of reaction per unit mass of the first reactant (units
$\text{J}\,\text{kg}^{-1}$). The **internal storage sign convention** (§2, overload note H1)
is

- $h_r > 0$ **endothermic** (the reaction *absorbs* heat). Then $Q_{\text{rxn}}<0$: the cell
  **cools**.
- $h_r < 0$ **exothermic** (the reaction *releases* heat). Then $Q_{\text{rxn}}>0$: the cell
  **heats**.

The leading minus sign in (7.16) converts the enthalpy-of-reaction convention into the
heat-source convention. This sign is the most error-prone in the whole assembly, so it is
worth stating precisely against the implementation: the accumulator computes
`Q -= heats[i] * rates[i]` — i.e. $Q=-\sum_r h_r r_r$ — with rates $r_r\ge 0$ for forward
reactions, so an exothermic reaction ($h_r<0$) produces $Q>0$.

> **Caution on external conventions.** ThermaKin and Gpyro publish $h>0=$ *exothermic*, the
> opposite of the internal storage convention here. The mapping is absorbed by the minus sign
> in (7.16): with their published $h$ one would write $Q_{\text{rxn}}=+\sum h r$.
> Pyrolysis.jl stores $h$ with $h>0=$ endothermic and forms $Q_{\text{rxn}}=-\sum h r$, so a
> user porting a ThermaKin heat of reaction must flip its sign (see §6).

> **Implementation note.** `reaction_heat_source(rates, heats)` in `src/Physics/kinetics.jl`
> performs the negated sum. The per-cell reaction rates and heats are evaluated by
> `evaluate_reactions!` (§6) at the cell's local $(T_i,\xi_i)$, and the stoichiometric species
> sources by `species_source_terms!`; `compute_rhs!` evaluates them inline, cell by cell.

### 7.5.2 Radiation source $Q_{\text{rad}}$ (dispatch only)

The volumetric radiation absorption source is computed by `compute_radiation_source!` and
dispatched on the active radiation model:

- **`NO_RADIATION`**: $Q_{\text{rad},i}=0$.
- **`SURFACE_ABSORPTION`**: $Q_{\text{rad},i}=0$ everywhere — the absorbed flux
  $\alpha\,q_{\mathrm{ext}}$ is applied through the surface energy balance
  (boundary Newton, §12.3.5), not as an $I/\Delta z$ volumetric deposit (opaque
  surface, better physics: the absorbed power raises $T_s$ directly).
- **`BEER_LAMBERT`**: in-depth absorption, $q'''(z)=\alpha(z)I(z)$ with
  $I(z)=I_{\text{surface}}\exp(-\!\int\!\alpha\,dz')$, assembled cell-exact via a backward
  sweep from the surface.
- **`P1_QUASI_STEADY`**: $Q_{\text{rad},i}=0$ here; P1 lives in `Pyrolysis.Experimental` with its own residual wrapper, and `solve` rejects it (§8.12).

The full derivation of each model — Beer–Lambert cell-averaged form, P1 diffusion, Marshak
boundary conditions, the energy-equation coupling, and the regime diagnostics — is the subject
of **§8 (Thermal Radiation)**. For the purposes of the energy assembly here, $Q_{\text{rad},i}$
is simply a cell-centered source term ($\text{W}\,\text{m}^{-3}$) added in (7.9). Note that the
boundary-incident radiative flux can also be applied through the *surface energy balance* (the
`RadiativeFluxBC` / boundary Newton solve, §12); which path a given model takes is part of §8's
coupling discussion.

---

## 7.6 The `compute_rhs!` Eight-Step Assembly Sequence

The complete energy + species RHS is assembled in `compute_rhs!`
(`src/Discretization/finite_volume.jl`) in a **fixed eight-step order**. The ordering is not
arbitrary: several steps consume results produced by earlier steps, and getting the order
wrong would break the surface energy balance, the boundary-cell $S_{\text{conv}}$, or the heat
divergence. One function serves
both transport laws: the law follows the property cache's Darcy sub-cache (`nothing` ⇔
Fickian), which selects the Darcy property refresh in Step 1 and the gas-flux form in Step 3
(§9). The sequence is:

**Step 1 — Transport properties and interior conductive fluxes.** In Darcy–Fick mode
`update_darcy_properties!` refreshes the cell pressure, porosity, permeability, and viscosity
and the face Darcy velocity (§9); in Fickian mode this is a no-op. `compute_conductive_fluxes!`
then fills $q_f$ at all interior faces via (7.1) and zeroes boundary faces (7.3). The result is
copied into the working flux vector `q_total` (`copyto!`), which will be modified at the
boundaries in Step 4.

**Step 2 — Reaction rates, heat source, and species sources.** For each cell,
`evaluate_reactions!` gives the rates $r_r$ and heats $h_r$; `reaction_heat_source` forms
$Q_{\text{rxn},i}$ (7.16); `species_source_terms!` forms $S_{j,\text{rxn},i}$, which is stored
in `dξ_rxn`. Step 8 reads both; in the ALE residual the strain rate and the mesh velocity
are built from `dξ_rxn` after `compute_rhs!` returns (§10.1.2, §11.3).

**Step 3 — Gas transport fluxes (interior + boundary).** `compute_gas_fluxes!` (Fickian
volume-fraction-gradient form) or `compute_gas_fluxes_darcy_fick!` (Darcy–Fick, with the
pressure-BC venting flux on the boundary faces) fills the gas-flux matrix `J_gas` (§9), then
`apply_mass_bc_flux!` adds the mass boundary conditions at the top and bottom faces. *This
must complete before the thermal-BC solve* (Step 4): the blowing model of each
`ConvectiveBC` reads the outward gas flux at its face from `J_gas` (§7.4.3, §12.3.3).

**Step 4 — Thermal boundary conditions.** `apply_thermal_bc_flux!` solves, at each
non-Dirichlet boundary, the surface energy balance for the surface temperature $T_s$ by Newton
iteration, and **overwrites** the boundary entries of `q_total` with the resulting conductive
flux $k_{\text{eff}}(T_s-T_{\text{cell}})/\Delta z_{\text{half}}$ (with the boundary sign
convention: $+$ at the bottom, $-$ at the top, so the divergence sees the correct
inflow/outflow). Dirichlet BCs set $T_s$ directly and evaluate the flux without iteration. The
full surface-balance derivation — the Newton residual $f(T_s)=q_{\text{BC}}(T_s)-
k_{\text{eff}}(T_s-T_{\text{cell}})/\Delta z_{\text{half}}$, its derivative, the
$[200,3000]\,\text{K}$ clamp, and the convergence criterion — is in §12.

**Step 5 — Radiation source.** `compute_radiation_source!` fills $Q_{\text{rad}}$ per the
active model (§7.5.2, §8).

**Step 6 — Gas convective source.** `compute_gas_convective_source!` fills $S_{\text{conv}}$
(7.13) from `J_gas`, the cell-centered temperature gradients, and — at the two boundary
cells — the boundary entries of `q_total` written in Step 4, from which the boundary-face
slope at the solved $T_s$ is recovered as $-q_{\rm total}[f]/k_i$ (7.15b). (It uses the
full face-averaged gas fluxes, including the boundary-applied values, now in `J_gas`.)

**Step 7 — Divergences.** `divergence_1d!` forms $(\nabla\!\cdot\mathbf{q})_i$ from
`q_total` (now including the boundary fluxes from Step 4), and `divergence_species_1d!` forms
$(\nabla\!\cdot\mathbf{J}_j)_i$.

**Step 8 — Final RHS assembly.** For each cell, the temperature derivative (7.9) is formed by
combining $-(\nabla\!\cdot\mathbf{q})_i$, $S_{\text{conv},i}$,
$Q_{\text{rad},i}$, $Q_{\text{rxn},i}$ and dividing by `prop_cache.ρcp[i]`; the species
derivative (7.10) is formed by `apply_species_residual_unrolled!`.

### 7.6.1 The ordering dependency graph

The dependencies that *fix* the order are:

1. **Gas fluxes (Step 3) → thermal BC solve (Step 4):** the blowing models of the
   convective BCs need the outward gas flux at the boundary.
2. **Thermal BC solve (Step 4) → gas convective source (Step 6) and heat divergence
   (Step 7):** both read the boundary fluxes that Step 4 wrote into `q_total` — the
   divergence directly, and $S_{\text{conv}}$ through the boundary cells' face slope
   $-q_{\rm total}[f]/k_i$ at the solved $T_s$ (7.15b).
3. **Reaction sources (Step 2) → everything downstream:** $Q_{\text{rxn}}$ and the species
   reaction sources are needed in Step 8; in the *ALE* residual (§11, §15) the reaction species
   sources also feed the strain rate $\theta_i$ and hence the cumulative mesh velocity $w$ — so
   reaction rates must be available *before* the mesh velocity is computed.

The full-residual entry point `residual!`, its two trait-dispatched paths (Eulerian and ALE),
and the ALE-specific extensions to this sequence (strain rate, mesh velocity, dilation
correction, $\chi$ dynamics) are described in §15
(Temporal Integration & Structured Jacobian) and §11 (ALE Framework). This chapter's
`compute_rhs!` is the Eulerian core that those paths build on.

---

## 7.7 Caches and Type Parameterization for AD

### 7.7.1 `RHSCache{Tu}`

All RHS intermediates are stored in a pre-allocated `RHSCache{Tu}`
(`src/Discretization/finite_volume.jl`), eliminating per-call allocations (verified by the
zero-allocation RHS test). Its fields:

| field | shape | meaning | units |
|---|---|---|---|
| `q_cond` | $n_{\text{faces}}$ | conductive face fluxes | $\text{W}\,\text{m}^{-2}$ |
| `q_total` | $n_{\text{faces}}$ | total energy face fluxes ($q_{\text{cond}}$ + BC) | $\text{W}\,\text{m}^{-2}$ |
| `S_conv` | $n_{\text{cells}}$ | gas convective source | $\text{W}\,\text{m}^{-3}$ |
| `J_gas` | $N_C\times n_{\text{faces}}$ | gas species face fluxes | $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$ |
| `Q_rxn` | $n_{\text{cells}}$ | reaction heat source | $\text{W}\,\text{m}^{-3}$ |
| `Q_rad` | $n_{\text{cells}}$ | radiation source | $\text{W}\,\text{m}^{-3}$ |
| `dξ_rxn` | $N_C\times n_{\text{cells}}$ | species reaction sources | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |
| `rates`, `heats` | $N_R$ | per-cell reaction workspace | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$, $\text{J}\,\text{kg}^{-1}$ |
| `species_src` | $N_C$ | per-cell species-source workspace | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |
| `div_q` | $n_{\text{cells}}$ | $\nabla\!\cdot\mathbf{q}$ | $\text{W}\,\text{m}^{-3}$ |
| `div_J` | $N_C\times n_{\text{cells}}$ | $\nabla\!\cdot\mathbf{J}_j$ | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |

The cache keeps a fixed size for the whole solve: depleted-cell parking (§14.3) leaves
ghost slots in place instead of changing the cell/face counts, and `compute_rhs!` loops
over the live cells `1:n_active_cells` only. It is constructed sized to a mesh and material
via `RHSCache(mesh, material)`, defaulting to `Tu = Float64`.

### 7.7.2 The scalar type parameter `Tu`

`RHSCache{Tu}`, `PropertyCache{Tu,DarcyT}`, and the typed state cache are all parameterized on
a single live scalar type `Tu`. In ordinary integration `Tu = Float64`. During forward-mode
automatic differentiation — used by the structured-Jacobian operators' stencil-local
ForwardDiff kernels, by the dense reference Jacobian of the verifier, and end-to-end for
sensitivity analysis (§17) —
`Tu = ForwardDiff.Dual`. Because every flux, source, property, and geometric quantity is held
at type `Tu`, the entire RHS assembly propagates dual sensitivities without a single
re-allocation or type conversion. The compile-time component unrolling (`@generated` functions
for the heat-capacity sums, $S_{\text{conv}}$, and the species residual) is
what makes this type stability achievable over the heterogeneous component tuple: the gas/solid
branches are resolved at code generation, so there is no runtime tuple indexing (which would
allocate) and the compiler specializes each arithmetic line to the concrete `Tu`.

> **Implementation note.** The production assembly uses `RHSCache` exclusively. The
> `PropertyCache` (§5) supplies `ρ`, `ρcp`,
> `k`, `λ`, `ε`, `a`, `α` per cell and `k_face`, `λ_face` per face (plus an optional Darcy sub-cache
> of `P`, `φ`, `κ`, `μ` per cell and `v_gas`, `∇P` per face); it must be refreshed by `update_properties!` before `compute_rhs!` is
> called.

---

## 7.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to the unified nomenclature on first use and contrast the
energy assembly. (The reference codes are condensed-phase pyrolysis / fire solvers:
ThermaKin2Ds by Stoliarov et al.; Gpyro by Lautenberger & Fernandez-Pello; the FDS solid-phase
model by McGrattan et al.)

### 7.8.1 The energy-equation skeleton

All four codes solve a 1D (through-thickness) Fourier-conduction energy balance with
volumetric reaction heat and a gas-enthalpy coupling. In our symbols the common skeleton is

$$
(\rho c_p)\,\frac{\partial T}{\partial t}
=
\frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right)
+ Q_{\text{rxn}} + Q_{\text{rad}} + (\text{gas enthalpy term}),
$$

and the codes differ mainly in (i) whether gas heat capacity is in the storage term, (ii) how
the gas enthalpy coupling is written, and (iii) the radiation closure.

### 7.8.2 ThermaKin2Ds

ThermaKin writes the through-thickness energy equation with the all-component storage
$\sum_j \xi_j c_{p,j}\,\partial T/\partial t$, the gas-advection term
$c_g\,(\mathbf{J}_g\!\cdot\!\partial T/\partial x)$ (its Eq. 17) — *exactly* the
$S_{\text{conv}}$ of (7.12) — and a deformation term that is the skeleton-following part
of the time derivative, which the ALE mesh frame absorbs (§11.2.2). Pyrolysis.jl's energy
equation is therefore ThermaKin Eq. 17 term for term. ThermaKin maps to our symbols as:
ThermaKin's stoichiometric $\theta_i^{\,j}\to \nu_{j,r}$ (the strain rate $\theta$ is
reserved here for dilation; §2 note A2); ThermaKin's bimolecular concentration product
$\xi_{\text{COMP1}}\xi_{\text{COMP2}}$ has no analog because Pyrolysis.jl enforces a single
reactant (§6); ThermaKin's gas flux $J=-\rho_g\lambda\,\partial(\xi_g/\rho_g)/\partial x$ is
the volume-fraction-gradient form implemented here (§9).

### 7.8.3 Gpyro

Gpyro's energy equation (Lautenberger) uses a control-volume conservative formulation on
$\bar\rho\,\Delta z$, with the bulk (averaged) density $\bar\rho$ (our $\rho$),
volume-fraction ($X_\alpha\to v_j$) property averaging, and power-law temperature-dependent
properties. Its thermal-equilibrium form carries the gas through the solid-to-gas heat
transfer $\dot Q_{s-g}''' \approx \dot m'' c_{pg}\,\partial T/\partial z$ (its Eq. 65) — our
$S_{\text{conv}}$ with the Darcy mass flux — and drops the gas storage
$\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default; its `FULL_QSG` option restores it,
which makes the equation (7.6) with (7.7). The mapping: Gpyro's pre-exponential
$Z\to A_r$; its conversion variable $\alpha_{\text{conv}}$ has no single analog
(Pyrolysis.jl tracks mass concentrations $\xi_j$ directly, with $\chi$ as a separate
progress variable for cross-section, §10); its permeability $K\to\kappa$; its mass flux
$\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ is the Darcy term of §9 (Pyrolysis.jl
omits the gravity/buoyancy term $g\rho_g$ as negligible in the thin condensed phase).
Numerically Gpyro uses Patankar control-volume discretization with Peclet-based scheme
selection, harmonic interface interpolation, and a TDMA solve; Pyrolysis.jl uses the same
harmonic interface interpolation (7.2) but a method-of-lines + stiff-ODE integrator with a
structured Jacobian (§15) rather than a hand-rolled implicit TDMA sweep.

### 7.8.4 FDS

The FDS solid-phase model solves a 1D Fourier conduction equation with the surface balance
$-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$, Arrhenius reaction kinetics
($\rho^n O_2^n T^n$ rate form), and no in-depth gas: volatiles leave instantaneously, so
there is no gas storage and no in-depth gas-enthalpy coupling; the gas enthalpy is closed
at the surface. FDS maps as: pre-exponential $Z\to A_r$; solid absorption coefficient
$\kappa_s\to\alpha$; the surface energy balance with $\dot q_c''$ (convective) and
$\dot q_r''$ (radiative) is the boundary Newton solve of §12. The radiation closure differs
(FDS uses a two-flux Schuster–Schwarzschild model in-depth; Pyrolysis.jl offers
SURFACE_ABSORPTION / BEER_LAMBERT / P1, §8). Pyrolysis.jl, following ThermaKin, resolves
the gas in depth and distributes its enthalpy transport volumetrically as $S_{\text{conv}}$
rather than as a surface outflow term, and the surface balance carries no gas-enthalpy term
on top of it (§7.4.3).

### 7.8.5 Summary table

| feature | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| storage $\rho c_p$ | all components (7.7) | all components (incl. gas) | matrix by default; all with `FULL_QSG` | matrix (no in-depth gas) |
| gas enthalpy coupling | $S_{\text{conv}}$ volumetric | $S_{\text{conv}}$ volumetric (Eq. 17) | $\dot m'' c_{pg}\,\partial T/\partial z$ (thermal equilibrium) | surface outflow |
| face conductivity | distance-weighted harmonic mean | harmonic | harmonic (Patankar) | — |
| conduction discretization | cell-centered FV, MOL + stiff ODE | Crank–Nicolson (x) | Patankar CV + TDMA | implicit FD |
| radiation | none / surface / Beer–Lambert / P1 | Beer–Lambert + re-radiation | in-depth exp absorption | two-flux S–S |
| heat-of-reaction sign | $h>0$ endothermic, $Q_{\text{rxn}}=-\sum hr$ | $h>0$ exothermic (flip) | $h>0$ exothermic (flip) | per-reaction |

---

## 7.9 Limitations and Open Issues

The energy assembly inherits a small number of known limitations, drawn from the source notes
and discussed where relevant above:

1. **One temperature.** The pore gas is stored at the matrix temperature (§7.3.2); there is
   no two-temperature gas energy equation, so regimes where the gas does not equilibrate
   with the pore walls within a cell are outside the model (§3.9).

2. **One global radiation model.** $Q_{\text{rad}}$ is dispatched on a single model for the
   whole domain; per-cell or temperature-dependent model switching is not supported (§8).

3. **Discrete energy-balance closure for composition-varying mixtures.** The energy ledger
   (§16) audits bookkeeping closure to integrator tolerance but does not carry the
   composition-change enthalpy $\sum_j h_j(T)\,\partial_t\xi_j$ as a conserved quantity
   (§3.3.2) or the ALE volume-change work, so the *physical* balance (as opposed to the
   discrete-bookkeeping balance) carries those as part of the quantified approximation cost
   rather than as an exactly-closed term.

These do not affect the correctness of the discrete update (7.9)–(7.10); they bound the
*modeling* fidelity of the energy equation and are exposed, not hidden, through the
conservation diagnostics of §16.


---

# 8. Thermal Radiation

This chapter derives and documents every thermal-radiation model implemented in
Pyrolysis.jl and the way each couples into the energy equation. The solver offers
four mutually exclusive volumetric-radiation treatments — `NO_RADIATION`,
`SURFACE_ABSORPTION`, `BEER_LAMBERT`, and a quasi-steady P1 (first-order
spherical-harmonics) approximation — each with a distinct cost, Jacobian
structure, and physical regime of validity. We give the continuous radiative
transfer equation (RTE), its reduction to the Beer–Lambert and P1 forms, the exact
cell-averaged discretizations, the Marshak/Reflective/Dirichlet boundary closures,
the allocation-free Newton/Thomas solver for the P1 field, and the non-local
(upper-triangular) Jacobian that Beer–Lambert induces. Every governing form,
discrete update, constant, sign, factor, guard, and limiter below has been verified
against the implementation; where a docstring or textbook form diverges from
the code, **the code is authoritative** and the divergence is
called out explicitly.

Throughout, the spatial coordinate is $z$ (axial, through-thickness), with $z=0$
the bottom/substrate (boundary id 2, `:bottom`) and $z=L$ the exposed/top surface
(boundary id 1, `:top`). Heat enters from the top. The flux sign convention is
**positive flux = transport in the $+z$ direction**; the discrete divergence
$(F_R-F_L)A/V$ is positive when the quantity flows *out* of a cell. The cell index
is $i$ ($1\le i\le n$), the component index is $j$ ($1\le j\le N_C$). The
Stefan–Boltzmann constant is $\sigma = 5.670374419\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}$
(`STEFAN_BOLTZMANN`, src/Materials/constants.jl).

---

## 8.1 The Radiation Model Family and its Traits

### 8.1.1 The four models

Volumetric radiation in Pyrolysis.jl is selected at `solve` time by a single
enumeration value passed as `radiation_model`. The enum is defined in
src/Materials/radiation_types.jl:

```julia
@enum RadiationModel begin
    NO_RADIATION = 0
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
    P1_QUASI_STEADY = 3
end
```

The models are **mutually exclusive** — exactly one is active for the whole domain
in a given run — and they differ in where the incident external radiation is
deposited and how the material's own re-emission is treated:

1. **`NO_RADIATION`.** No volumetric radiation source ($Q_{\text{rad}}=0$).
   External radiative *boundary* conditions (`RadiativeBC` surface re-radiation,
   `RadiativeFluxBC` incident flux) still act at the surface through the boundary
   energy balance (§8.8, and §12, Boundary & Initial Conditions). Use for
   optically transparent media ($\tau \ll 0.1$) or when in-depth absorption is
   negligible.

2. **`SURFACE_ABSORPTION`.** The incident external radiation is absorbed
   entirely at the exposed surface. The absorbed flux is **not** deposited as a
   volumetric source; it enters the energy balance as a *boundary heat flux*
   through the surface-temperature Newton solve (§8.2, §8.8). Component absorption
   coefficients $\alpha_j$ are ignored. Cost $O(1)$; tridiagonal Jacobian. This is
   the recommended choice for opaque solids with first-cell optical thickness
   $\tau\gg 1$.

3. **`BEER_LAMBERT`.** In-depth exponential absorption following the Beer–Lambert
   law, $I(z)=I_{\text{surface}}\exp\!\big(-\int_L^z\alpha_{\text{eff}}\,\mathrm{d}z'\big)$,
   with volumetric source $q'''=\alpha_{\text{eff}}I$. Radiation penetrates over the
   absorption length $1/\alpha_{\text{eff}}$. Cost $O(n)$; the intensity at cell $i$
   depends on all cells above it, producing a **non-local, upper-triangular**
   Jacobian coupling (§8.10). Appropriate for semi-transparent polymers and foams
   ($0.1\lesssim\tau\lesssim 10$).

4. **`P1_QUASI_STEADY`.** The P1 diffusion approximation, solving an auxiliary
   incident-radiation field $G$ from a quasi-steady diffusion equation, including
   *volumetric re-emission* $4\sigma\alpha T^4$. Cost $\approx O(10n)$ (a small
   Newton loop per step); the $G$ equation is tridiagonal. Recommended when
   internal re-radiation matters ($T\gtrsim 800$ K, optically thick char). P1 lives
   in the `Experimental` module and runs through a *separate residual* that
   re-solves the quasi-steady $G$ field into its radiation cache at every
   evaluation — $G$ is not part of the ODE state. `solve` rejects this model;
   the P1 kernels are called directly (§8.12).

### 8.1.2 Model traits

The enum carries one trait query: `requires_radiation_cache(model)`
(src/Materials/radiation_types.jl) is true for `P1_QUASI_STEADY` only, where a
`RadiationCache` for $G$, $D_{\mathrm{rad}}$, and Newton work arrays is allocated.
The other model-dependent decisions are made where they arise. The Jacobian
sparsity builder widens the pattern to the upper-triangular Beer–Lambert block when
the simulation mode's radiation axis is `BeerLambert` (§8.10, §15.11), and `solve`
checks the material against the model (`_validate_radiation_model`,
src/Solver/solve.jl): it warns when `SURFACE_ABSORPTION` meets non-zero
$\alpha_j$ (ignored) and when `BEER_LAMBERT` meets none (no absorption).

The performance/structure summary (radiation_types.jl):

| Model | Loop cost | Jacobian pattern | Extra memory |
|---|---|---|---|
| `NO_RADIATION` | $O(1)$ | tridiagonal | baseline |
| `SURFACE_ABSORPTION` | $O(1)$ | tridiagonal | baseline |
| `BEER_LAMBERT` | $O(n)$ | upper-triangular | $+O(n^2)$ (dense fill) |
| `P1_QUASI_STEADY` | $O(10n)$ | tridiagonal | $+O(n)$ |

> **Implementation note.** The dispatch `compute_radiation_source!`
> (`q_rad, state, mesh, cache, I_surface, model, apply_kirchhoff_absorptivity`,
> src/Physics/radiation.jl)
> *zeros* `q_rad` for `NO_RADIATION`, `P1_QUASI_STEADY`, **and**
> `SURFACE_ABSORPTION`, and only calls `compute_beer_lambert_source!` for
> `BEER_LAMBERT`. P1 is added on a separate residual path (§8.12). The boundary
> energy balance accounts for the surface-absorbed radiation in `NO_RADIATION` and
> `SURFACE_ABSORPTION` modes.

---

## 8.2 Surface Absorption

### 8.2.1 Surface absorption is a boundary flux, not a volumetric source

In the assembled RHS, the `SURFACE_ABSORPTION` model contributes *zero*
volumetric source — `compute_radiation_source!` zeroes `q_rad`
(src/Physics/radiation.jl) — and the incident radiation is folded into the
surface energy balance through the `RadiativeFluxBC` dispatch
(boundary_temperature.jl):

$$
q_{\text{BC}}(T_s) = \Phi_{\text{incident}}\,\alpha
\quad\text{for }\texttt{NO\_RADIATION}\text{ or }\texttt{SURFACE\_ABSORPTION},
$$

where $\Phi_{\text{incident}}=I_{\text{ext}}(t)$ is the incident intensity
[W·m⁻²] and $\alpha$ is the surface absorptivity. The default absorptivity is the
material's effective surface absorptivity
$\alpha=a_{\text{eff}}=\sum_j v_j^{\text{solid}}a_j$ (§5.8.2, clamped to
$[0,1]$), where $a_j$ is the component's `absorptivity` when set and otherwise
its emissivity $\varepsilon_j(T)$ — so with no component absorptivity set this
is Kirchhoff's law $\alpha=\varepsilon_{\text{eff}}$; an explicit BC
`absorptivity` overrides this. The temperature derivative
of this contribution is $\mathrm{d}q_{\text{BC}}/\mathrm{d}T_s=0$ (absorbed
incident flux is independent of $T_s$). Surface re-radiation loss is supplied
separately by a `RadiativeBC` term, $q=F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$,
whose $T_s$-dependence does enter the Newton derivative (§12).

The textbook alternative — depositing the absorbed flux as a volumetric source
$I/\Delta z_N$ in the top cell — is deliberately **not** implemented. The
boundary-flux route is the better construction: it lets the surface-temperature
Newton solve resolve the (possibly large) jump between $T_s$ and the cell-center
temperature $T_N$ consistently with convection and re-radiation; for an opaque
material with a well-resolved first cell the two constructions would be
numerically close anyway.

### 8.2.2 Assumptions and validity

Surface absorption presumes the material is optically thick over the first cell,
$\tau_N=\alpha_{\text{eff}}\Delta z_N\gg 1$, so that essentially no radiation
penetrates beyond cell $N$. It is $O(1)$, adds no Jacobian non-locality (the
tridiagonal conduction stencil is preserved), and is the recommended choice for
charring solids and fast cone-calorimeter analysis.

---

## 8.3 The Beer–Lambert In-Depth Model

### 8.3.1 Continuous form

Beer–Lambert follows from the 1D radiative transfer equation for a collimated beam
along $z$ with absorption coefficient $\alpha_{\text{eff}}(z)$ and negligible
scattering/emission of the *incident* beam:

$$
\frac{\partial I}{\partial z} = -\alpha_{\text{eff}}\, I .
$$

With the beam entering at the top surface $z=L$ and travelling into depth (toward
$z=0$), integration gives exponential attenuation:

$$
I(z) = I_{\text{surface}}\,\exp\!\left(-\int_{L}^{z}\alpha_{\text{eff}}(z')\,\mathrm{d}z'\right).
$$

The local energy deposited per unit volume is the absorbed power density,

$$
q'''(z) = \alpha_{\text{eff}}(z)\,I(z)\qquad[\mathrm{W\,m^{-3}}].
$$

The effective absorption coefficient is the composition-weighted volumetric
extinction,

$$
\alpha_{\text{eff}} = \sum_j \xi_j\,\alpha_j,
$$

where $\xi_j$ [kg·m⁻³] is the mass concentration and $\alpha_j$ [m²·kg⁻¹] is the
per-component mass absorption coefficient (`effective_absorption`,
src/Properties/effective_properties.jl). The product $\xi_j\alpha_j$ has
units m⁻¹, so $\alpha_{\text{eff}}$ is a *volumetric* coefficient; it is stored in
the property cache as `cache.α[i]` (see §5, Effective Properties & Mixing Rules).
The optical thickness of the whole layer is
$\tau = \int_0^L \alpha_{\text{eff}}\,\mathrm{d}z$.

### 8.3.2 The exact cell-averaged discretization

A naive evaluation $q'''_i \approx \alpha_i I(z_i^c)$ at the cell center is only
first-order accurate and can over- or under-deposit energy depending on
$\tau_i$. Pyrolysis.jl instead uses the **exact cell-integrated** source: integrate
$q'''=\alpha I$ over a cell of constant $\alpha_i$ and entrance intensity
$I_{\text{in}}$,

$$
\int_{\text{cell }i} \alpha_i\,I\,\mathrm{d}z
= \int_0^{\Delta z_i}\!\alpha_i\, I_{\text{in}}\, e^{-\alpha_i s}\,\mathrm{d}s
= I_{\text{in}}\big(1 - e^{-\alpha_i \Delta z_i}\big),
$$

so the *cell-averaged* volumetric source is the absorbed power divided by the cell
thickness:

$$
\boxed{\;q^{\text{rad}}_i = I_{\text{in}}\,\frac{1 - e^{-\tau_i}}{\Delta z_i},
\qquad \tau_i \equiv \alpha_i\,\Delta z_i\;}
$$

This is exact for piecewise-constant $\alpha$ and conserves energy by
construction: the power absorbed in cell $i$ is $q^{\text{rad}}_i\,\Delta z_i =
I_{\text{in}}(1-e^{-\tau_i})$, exactly the drop in beam power across the cell.

The flux that actually crosses the surface is $\alpha_{\text{surf}}\,
I_{\text{surface}}$, where $\alpha_{\text{surf}}$ is the surface
transmissivity: on the default **Kirchhoff path** (the `RadiativeFluxBC`
carries no explicit `absorptivity`) it is the *live* surface cell's effective
surface absorptivity $a_{\text{eff}}(T_N,\xi_{\cdot,N})$ (`cache.a[N]`,
§5.8.2) — identical to $\varepsilon_{\text{eff}}(T_N,\xi_{\cdot,N})$ unless a
component sets its own `absorptivity`, so by default absorption and surface
emission use the same coefficient (Kirchhoff's law $\alpha=\varepsilon$,
detailed balance) — and the remaining $(1-a_{\text{eff}})I_{\text{surface}}$
is reflected and never enters the medium. With an **explicit BC absorptivity**, $\alpha(t)$ is folded
into the extracted intensity upstream and the kernel does not rescale
(`apply_kirchhoff_absorptivity = false`). The entrance intensity is then the
cumulative attenuation from the surface down to the top face of cell $i$:

$$
I_{\text{in}}(i) = \alpha_{\text{surf}}\,I_{\text{surface}}\prod_{j=i+1}^{N} e^{-\tau_j}
= \alpha_{\text{surf}}\,I_{\text{surface}}\,\exp\!\left(-\sum_{j>i}\alpha_j\Delta z_j\right).
$$

The transmitted (un-absorbed) fraction leaving the bottom is
$\alpha_{\text{surf}}\,I_{\text{surface}}\,e^{-\tau}$; it leaves through the back face (§8.3.3).

### 8.3.3 Backward sweep and guards

Because $I_{\text{in}}(i)$ depends only on cells *above* $i$ (larger index, nearer
the surface), the implementation sweeps from the live top cell
$n_{\rm live} = $ `mesh.n_active_cells` down to cell 1, carrying a running
intensity `I_current`; parked ghost slots above the live surface (§14.3) get
$q^{\text{rad}} = 0$. Forward iteration would require not-yet-computed
information; the backward sweep is the natural causal order.

```julia
n_cells_val = mesh.n_active_cells
fill!(q_rad, zero(eltype(q_rad)))
I_current = oneunit(eltype(q_rad)) * I_surface
if apply_kirchhoff_absorptivity
    I_current *= cache.a[n_cells_val]   # α_surf = a_eff of the live surface cell
end
A = state.A_section[]
@inbounds for i = n_cells_val:-1:1
    α  = cache.α[i]
    Δz = state.volume[i] / A
    if α > 1e-20
        τ = α * Δz
        absorption_factor = 1 - exp(-τ)
        q_rad[i] = I_current * absorption_factor / Δz
        I_current *= exp(-τ)
    else
        q_rad[i] = zero(eltype(q_rad))
    end
end
```
(`compute_beer_lambert_source!`, src/Physics/radiation.jl.)

Four details are load-bearing:

- **Cell thickness is geometric.** $\Delta z_i = V_i / A$ from `state.volume[i]`
  and the live cross-section `state.A_section[]`, so the source remains exact under
  ALE mesh motion and variable cross-section (the lateral-shrinkage area cancels
  between $V_i$ and $A$; see §10, Volume Change & Lateral Shrinkage, and §11, The
  ALE Framework).
- **Transparency guard.** Cells with $\alpha\le 10^{-20}$ get $q^{\text{rad}}_i=0$
  and do not attenuate the beam, preventing `exp(0)`-degenerate work and spurious
  absorption in radiation-free regions.
- **Zero-flux short-circuit.** If $I_{\text{surface}}\le 0$, the routine returns
  after zeroing the array.
- **Transmitted flux.** Whatever intensity is left after cell 1,
  $I_{\text{surface}}\,e^{-\tau_{\text{tot}}}$, leaves through the back face: the
  sweep deposits nothing there and the back-face boundary condition describes
  only the surface it models. That is the right answer for a slab whose back
  face is open and a modelling error for one that rests on an opaque holder
  with a transparent absorption coefficient, so the extractor reports the
  transmitted flux per saved time (`solution.extras.radiation_transmitted`)
  and a solve that loses more than 1 % of the entering flux at any saved time
  warns once (`_warn_transmitted_radiation`, src/Solver/solve.jl).

The source is **positive = heating**: it is added to the energy RHS as
$+q^{\text{rad}}_i/\rho c_p^{\text{eff}}$ (§8.8). Energy conservation holds in the
strong sense $\sum_i q^{\text{rad}}_i\,V_i = \alpha_{\text{surf}}\,
I_{\text{surface}}\,(1-e^{-\tau_{\text{tot}}})\,A \le \alpha_{\text{surf}}\,
I_{\text{surface}}\,A$, with equality only if the layer is fully opaque
($e^{-\tau}\to 0$) — consistent with `SURFACE_ABSORPTION`, which applies
$\varepsilon_{\text{eff}}\,I_{\text{surface}}$ at the face (switching models
therefore does not change the absorbed power in the opaque limit).

### 8.3.4 Where $I_{\text{surface}}$ comes from

The surface intensity fed to Beer–Lambert is the incident flux carried by the
top-boundary `RadiativeFluxBC`. `extract_radiative_flux_from_bcs`
(src/BoundaryConditions/boundary_conditions.jl) returns
`(I, has_explicit_absorptivity)`: with an explicit BC `absorptivity` it returns
$I_{\text{ext}}\alpha(t)$ and `true`; otherwise the bare incident value and
`false` — the live surface-cell $\varepsilon_{\text{eff}}$ is then applied at
injection inside `compute_beer_lambert_source!` (the Kirchhoff path). `solve()`
derives the static `radiation_kirchhoff` flag from the BC once and threads it
through `to_problem_def` to both the residual and the structured Jacobian, so
the two conventions can never mix.

Beer–Lambert models **top-incident** radiation only: the sweep starts at the
exposed surface, and the per-BC face flux of a `RadiativeFluxBC` is zeroed at
*every* boundary under this model. A `RadiativeFluxBC` on the **bottom**
boundary therefore cannot be honored, and `solve()` rejects the combination
with an `ArgumentError` rather than silently dropping the back-side heating
(src/Solver/solve.jl); use `SURFACE_ABSORPTION` (face deposit at either
boundary) or a `HeatFluxBC` for back-side radiant heating.

The solver feeds this flux to the ODE by re-sampling it at **every** RHS
evaluation — `ws.I_radiation = radiation_flux(t)` — whenever the active model
consumes it (`SURFACE_ABSORPTION` or `BEER_LAMBERT`; src/Solver/solve.jl). The
Jacobian closure refreshes `ws.I_radiation` the same way at its own evaluation
time, and the radiation operator reads it from there, so the Beer–Lambert
Jacobian sees the flux at the time it is evaluated. There is no constant-flux
detection: one user-function call per RHS is
negligible against the residual cost, and unconditional re-sampling guarantees
that a ramped or modulated cone flux can never be frozen at its $t=0$ value.

> **Implementation note.** Snapshotting a time-varying incident flux at $t=0$
> would silently freeze a cone-heater ramp — a silent error. The unconditional
> per-RHS re-sampling forecloses this failure mode by construction; a regression
> with $q_0(1+\sin 2\pi t)$ guards it (§18).

### 8.3.5 Worked numbers

For a 10 mm semi-transparent slab with $\alpha_{\text{eff}}=100\ \mathrm{m^{-1}}$,
uniform $\Delta z=1\ \mathrm{mm}$, $I_{\text{surface}}=50\ \mathrm{kW\,m^{-2}}$,
$\varepsilon_{\text{eff}}=0.9$ (Kirchhoff path): the flux entering the medium
is $0.9\times50 = 45\ \mathrm{kW\,m^{-2}}$. With $\tau_i=0.1$ the top cell
absorbs $q_1=45\,000\,(1-e^{-0.1})/0.001 \approx 4.28\ \mathrm{MW\,m^{-3}}$;
the beam entering the second cell is $45\,000\,e^{-0.1}=40.7\ \mathrm{kW\,m^{-2}}$,
giving $q_2\approx 3.88\ \mathrm{MW\,m^{-3}}$, and so on. The total optical
thickness is $\tau=\alpha L=1$, so $e^{-\tau}\approx 0.37$ of the entering flux
is transmitted through the slab — Beer–Lambert is the right tool in this regime,
whereas surface absorption would (incorrectly) deposit everything at the surface.

---

## 8.4 The P1 (First-Order Spherical-Harmonics) Approximation

### 8.4.1 From the RTE to a diffusion equation for $G$

For an absorbing–emitting (and, in principle, isotropically scattering) gray
medium, expanding the directional intensity in spherical harmonics and truncating
at first order (P1) yields a diffusion equation for the **incident radiation**
$G=\int_{4\pi} I\,\mathrm{d}\Omega$ [W·m⁻²]:

$$
\boxed{\;\nabla\!\cdot\!\big(D_{\mathrm{rad}}\,\nabla G\big) = \alpha\,G - 4\sigma\alpha T^4\;}
$$

with the P1 radiation diffusion coefficient

$$
D_{\mathrm{rad}} = \frac{1}{3(\alpha+\sigma_s)}\;\xrightarrow{\sigma_s\to0}\;\frac{1}{3\alpha}
\qquad[\mathrm{m}],
$$

where $\alpha$ is the volumetric absorption coefficient, $\sigma_s$ is the
scattering coefficient (taken $\approx 0$ for pyrolysis), and $\sigma$ is the
Stefan–Boltzmann constant. The left side is radiative *diffusion* (transport by
repeated absorption/re-emission); the right side balances removal from the field
by absorption ($\alpha G$) against blackbody emission ($4\sigma\alpha T^4$).

> **Implementation note.** `compute_radiation_diffusion_coefficient`
> (src/Experimental/radiation_p1.jl) returns
> $1/\big(3\max(\alpha+\sigma_s,\,10^{-10})\big)$ with $\sigma_s=0$ by default —
> the transparent limit is clamped to the large finite value
> $1/(3\cdot 10^{-10})$ rather than zeroed, because $D=0$ would delete the $G$
> row of the quasi-steady Newton system wherever $\alpha=0$ (singular
> tridiagonal matrix). This is the **single** canonical guard: the one
> `update_radiation_properties!` method routes every $D_{\mathrm{rad}}$ through
> it (see open issue 8 in §8.13).

### 8.4.2 The net volumetric source and its sign

Integrating the RTE over solid angle gives the radiative-flux divergence
$\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}} = \alpha(4\sigma T^4 - G)$. The energy
equation carries $-\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}}$, so the **net
volumetric radiation source** that heats the matrix is

$$
\boxed{\;S_{\mathrm{rad}} = \alpha\,(G - 4\sigma T^4)\;}
$$

with the convention

- $S_{\mathrm{rad}}>0$: net absorption ($G>4\sigma T^4$) — the material **heats**;
- $S_{\mathrm{rad}}<0$: net emission ($T^4$ exceeds $G$) — the material **cools**;
- $S_{\mathrm{rad}}=0$: local radiative equilibrium.

This is what the active P1 residual uses: `compute_radiation_flux_divergence!`
(src/Experimental/radiation_fluxes.jl) computes the **algebraic** source
$S_{\mathrm{rad}}[i] = \alpha_i(G_i - 4\sigma T_i^4)$ cell-by-cell and stores it in
`rad_cache.S_rad`, and the residual adds $\mathrm{d}T_i/\mathrm{d}t \mathrel{+}=
S_{\mathrm{rad}}[i]/\rho c_p^{\text{eff}}_i$ (`compute_p1_residual!`,
RadiationP1.jl). This algebraic identity is equivalent to $-\nabla\!\cdot\!
\mathbf{q}_{\mathrm{rad}}$ *only when the discrete $G$ field exactly satisfies the
P1 equation*; the Newton solve below drives it there.

> **Algebraic, not flux, divergence.** The flux form
> $S=(q_{\mathrm{rad}}[R]-q_{\mathrm{rad}}[L])A/V$ built from the diffusive face fluxes
> $q_{\mathrm{rad}}=-D_{\mathrm{rad}}\nabla G$ would add
> $\mathrm{d}T_i/\mathrm{d}t \mathrel{+}= -S/\rho c_p$ (the minus sign because
> $\nabla\!\cdot\!(D\nabla G) = -\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}}$). The P1
> residual (`compute_p1_residual!`, RadiationP1.jl) uses the *algebraic*
> `compute_radiation_flux_divergence!` instead; the two agree at the converged $G$.

### 8.4.3 The quasi-steady assumption

Radiation equilibrates orders of magnitude faster than conduction or reaction
($\tau_{\mathrm{rad}}\sim 10^{-12}$ s vs. $\tau_{\mathrm{cond}}\sim 10^{2}$ s), so
the field is treated as instantaneous:

$$
\frac{\partial G}{\partial t}\approx 0 .
$$

Operationally, $G$ is **not** a state variable at all: it lives in the
`RadiationCache` and is re-solved as a nonlinear algebraic system by the Newton
iteration of §8.7 at *every* residual evaluation (§8.12). The temperature
non-linearity $T^4$ makes the $G$ system nonlinear even though it is linear in
$G$ alone at fixed $T$.

---

## 8.5 P1 Volumetric Source Terms

The two pieces of the right-hand side of the P1 equation are computed by
`compute_radiation_source_terms!` (src/Experimental/radiation_p1.jl):

**Emission** [W·m⁻³]:

$$
Q_{\mathrm{em}}[i] = 4\sigma\,\alpha_i\,T_i^4 .
$$

**Absorption** [W·m⁻³]:

$$
Q_{\mathrm{abs}}[i] = \alpha_i\,G_i .
$$

The net source used in the energy equation is $S_{\mathrm{rad}} = Q_{\mathrm{abs}} -
Q_{\mathrm{em}} = \alpha(G - 4\sigma T^4)$ (§8.4.2), consistent with
`net_radiation_source` (src/Properties/radiation_cache.jl).

> **No emissivity factor in the volumetric emission — the code is authoritative.**
> Several references (and even the docstrings of `compute_radiation_source_terms!`
> and the `RadiationCache.Q_emission` field) write emission as $4\sigma
> \varepsilon\alpha T^4$. The executed code (radiation_p1.jl) uses **blackbody**
> emission $4\sigma\alpha T^4$ with **no** $\varepsilon$ factor. The physical
> justification is that in the P1 approximation the interior is a gray body in local
> thermodynamic equilibrium ($\varepsilon=1$ locally): a control volume that absorbs
> $\alpha G$ must, at equilibrium $G=4\sigma T^4$, emit exactly $4\sigma\alpha T^4$,
> giving $S_{\mathrm{rad}}=0$. Surface emissivity $\varepsilon<1$ enters **only** at
> the boundaries through the Marshak condition (§8.6), where it represents partial
> reflection of the diffuse field, not reduced volumetric emission. Writing
> $\varepsilon$ in the interior emission would break the equilibrium identity.

The $\alpha$ used here is read from `rad_cache.α` when populated (and positive),
otherwise from `prop_cache.α` (radiation_p1.jl); both hold the same
volumetric $\alpha_{\text{eff}}=\sum_j\xi_j\alpha_j$. Caching $\alpha$ on the
radiation cache avoids recomputing it in the source, diffusion-coefficient, and
Jacobian kernels.

---

## 8.6 P1 Boundary Conditions

P1 boundary conditions act on the radiation field $G$, not on $T$, and are stored
in a type-stable `RadiationBCSet{B1,B2}` with `top` and `bottom` slots
(src/Experimental/radiation_boundary.jl), accessed by
`get_radiation_bc(set, Val(1))` (top) / `Val(2)` (bottom). The **default** set is
a Marshak condition at the top and a reflective condition at the bottom
(`default_radiation_bcs`, radiation_boundary.jl).

### 8.6.1 Marshak condition (exposed surface)

The Marshak (mixed) condition closes the P1 half-range fluxes at a gray,
diffuse surface of emissivity $\varepsilon$. With $q_{\text{in}} = G/4 + q/2$
and $q_{\text{out}} = G/4 - q/2$ (q = net flux along the inward normal), the
inward half-range flux just inside the surface equals the transmitted external
radiosity plus the reflected outgoing half-range flux,

$$
\frac{G}{4} + \frac{q}{2}
= \varepsilon\big(I_{\text{ext}} + \sigma T_{\text{amb}}^4\big)
+ (1-\varepsilon)\Big(\frac{G}{4} - \frac{q}{2}\Big).
$$

Solving for the net flux into the domain gives the relation implemented by
`radiation_bc_flux` (src/Experimental/radiation_p1.jl):

$$
\boxed{\;q_{\text{BC}} = c\,\big(G_{\text{ext}} + 4 I_{\text{ext}} -
G_{\text{surf}}\big),\qquad c = \frac{\varepsilon}{2(2-\varepsilon)}\;}
$$

with

- $I_{\text{ext}}=I_{\text{external}}(t)$ — external incident radiation [W·m⁻²]
  (cone heater, radiant panel, flame). It carries the half-range weight
  $4c = 2\varepsilon/(2-\varepsilon)$ — **the same** weight the ambient term
  carries ($c\,G_{\text{ext}} = c\cdot 4\sigma T_{\text{amb}}^4 \equiv
  4c\cdot\sigma T_{\text{amb}}^4$): 2.0 at $\varepsilon=1$, 1.64 at
  $\varepsilon=0.9$. A bare weight of 1 on $I_{\text{ext}}$ under-injects the
  external flux by that factor (for a cold optically thick absorber at
  $\varepsilon = 0.9$ the self-consistent absorbed fraction is
  $4c/(1+\sqrt{3}c) \approx 0.96$ with the correct weight versus
  $1/(1+\sqrt{3}c) \approx 0.59$ with weight 1);
- $G_{\text{ext}} = 4\sigma T_{\text{amb}}^4$ — equilibrium incident field of the
  surroundings at the ambient temperature $T_{\text{amb}}(t)$;
- $G_{\text{surf}}$ — the P1 incident field at the surface;
- $\varepsilon=\varepsilon(T_s)$ — surface emissivity (constant, $\varepsilon(T)$,
  or BC-supplied).

```julia
c     = _marshak_exchange_coefficient(bc, T_surf)        # ε/(2(2−ε))
flux  = c * (_marshak_driving_field(bc, t) - G_surf)     # driving = 4σT_amb⁴ + 4I_ext
```

This formulation satisfies the required limits: at equilibrium
($G=4\sigma T^4$, $T_{\text{amb}}=T_{\text{surf}}$, $I_{\text{ext}}=0$) the flux
vanishes; the external and ambient sources enter with consistent half-range
weights; and a fraction $(1-\varepsilon)$ of the interior radiation is reflected
back into the domain.

**Half-cell-corrected application.** The solver applies the relation at the
boundary-cell *center*, eliminating the face value $G_{\text{face}}$ through
the discrete diffusive flux $q = D_{\mathrm{rad}}(G_{\text{face}} - G_N)/\delta$
($\delta$ = center-to-face distance). The result is the same relation with the
series-resistance coefficient

$$
q_{\text{BC}} = c_{\text{eff}}\,\big(G_{\text{ext}} + 4I_{\text{ext}} - G_N\big),
\qquad c_{\text{eff}} = \frac{1}{1/c + \delta/D_{\mathrm{rad}}},
$$

(`_marshak_effective_coefficient`), which removes the $O(\delta)$ error of a
cell-centered evaluation. The flux application (`_apply_marshak_bc!`) and the
quasi-steady Newton Jacobian use the identical $c_{\text{eff}}$.
$\varepsilon(T)$ and $\alpha$ are evaluated at the boundary-cell temperature —
the surface-Newton $T_s$ is not available inside the radiation solve.

> **The $\varepsilon/(2(2-\varepsilon))$ coefficient is the correct Marshak form**
> (Modest, *Radiative Heat Transfer*, 3rd ed., Eq. 16.48), **not** $\tfrac{1}{2}\varepsilon$.
> A common simplification, $0.5\varepsilon$, silently *underestimates* radiative
> exchange for $\varepsilon<1$ (the implemented coefficient is
> `_marshak_exchange_coefficient`, radiation_p1.jl). For $\varepsilon=1$ both
> reduce to $\tfrac12$; the discrepancy grows as $\varepsilon$ falls (e.g. for
> $\varepsilon=0.5$ the correct coefficient is $0.5/(2\cdot1.5)=1/6$ vs. $0.25$).

The flux is applied at the boundary face with the sign-orientation convention of
the divergence operator (radiation_p1.jl): at a **right** boundary
(`right_id==0`, the top) the stored face flux is $-q_{\text{BC}}$; at a **left**
boundary (the bottom) it is $+q_{\text{BC}}$. Because the `RadiationCache.q_rad`
array is `Float64`, the Marshak flux is stored via `ForwardDiff.value` to strip any
dual part during the Float64-only P1 solve (radiation_p1.jl).

The Marshak BC is constructed as `MarshakRadiationBC(; I_external, T_ambient=300.0,
emissivity=0.9)` (radiation_boundary.jl) with validation that
$I_{\text{external}}\ge 0$, $T_{\text{ambient}}>0$, and $\varepsilon\in[0,1]$. Any of
the three fields may be a constant or a callable of $t$ (or $T$ for emissivity).

### 8.6.2 Reflective condition (insulated / symmetry surface)

`ReflectiveRadiationBC` (radiation_boundary.jl) imposes a zero normal gradient,

$$
\frac{\partial G}{\partial n} = 0,
$$

realized by setting the boundary face flux to zero (`_apply_reflective_bc!`,
radiation_p1.jl). It is the natural choice for an insulated back face, a
symmetry plane, or a reflective substrate. It is the default bottom BC.

### 8.6.3 Dirichlet condition (prescribed $G$)

`DirichletRadiationBC(; G)` (radiation_boundary.jl) prescribes the
boundary intensity directly. The current implementation is a placeholder: the
boundary face flux is set to zero (`_apply_dirichlet_bc!`, radiation_p1.jl;
`radiation_bc_flux` returns `zero(t)`), so $G$ is effectively pinned only through
the interior gradient. It is used mainly for validation against analytic solutions
and is rarely used in practice; Marshak is preferred because it couples the
external source to surface emission physically.

---

## 8.7 P1 Numerics: Diffusion Faces, Fluxes, and the Newton/Thomas Solve

### 8.7.1 Face diffusion coefficient (harmonic mean)

Interior face values of $D_{\mathrm{rad}}$ use a **distance-weighted harmonic
mean** (`update_radiation_properties!`, radiation_p1.jl):

$$
D_{\mathrm{rad,face}} = \frac{d_L + d_R}{\dfrac{d_L}{D_{\mathrm{rad},L}} + \dfrac{d_R}{D_{\mathrm{rad},R}}},
$$

where $d_L,d_R$ are the distances from the left/right cell centers to the face. The
harmonic mean enforces flux continuity across the interface (resistance in series),
exactly as for thermal conductivity (§7, Heat Conduction) and the transport
coefficients (§5). At a boundary face the value of the adjacent interior cell is
used.

### 8.7.2 Interior diffusive flux

The radiative flux at an interior face is the discrete diffusive flux
(radiation_p1.jl):

$$
q_{\mathrm{rad}}[f] = -D_{\mathrm{rad,face}}\,\frac{G_R - G_L}{\Delta x},
\qquad \Delta x = |z_R^c - z_L^c|,
$$

with the standard sign convention (positive flux toward $+z$). Boundary faces are
overwritten by the BC handlers of §8.6 in `_apply_all_radiation_bcs!`
(radiation_p1.jl).

### 8.7.3 The quasi-steady residual

`solve_quasi_steady_radiation!` (src/Experimental/radiation_fluxes.jl)
solves the discrete P1 equation for $G$ by Newton–Raphson. The per-cell residual
combines the algebraic source terms with the discrete flux divergence
(radiation_fluxes.jl):

$$
R_i = Q_{\mathrm{abs}}[i] - Q_{\mathrm{em}}[i] + S_{\mathrm{rad,div}}[i],
\qquad
S_{\mathrm{rad,div}}[i] = \big(q_{\mathrm{rad}}[R] - q_{\mathrm{rad}}[L]\big)\frac{A}{V_i}.
$$

This is the discrete form of $R = \alpha G - 4\sigma\alpha T^4 -
\nabla\!\cdot\!(D_{\mathrm{rad}}\nabla G)$: since $q_{\mathrm{rad}}=-D_{\mathrm{rad}}
\nabla G$, the flux divergence $S_{\mathrm{rad,div}} = -\nabla\!\cdot\!(D_{\mathrm{rad}}
\nabla G)$, and $R = (Q_{\mathrm{abs}}-Q_{\mathrm{em}}) + S_{\mathrm{rad,div}}$.

### 8.7.4 The tridiagonal Jacobian

Differentiating $R_i$ with respect to the unknowns $G$ gives a tridiagonal Jacobian
(radiation_fluxes.jl). With
$\mathrm{coeff}_L = \tfrac{A}{V_i}\tfrac{D_{\mathrm{rad},L}}{\Delta x_L}$ and
$\mathrm{coeff}_R = \tfrac{A}{V_i}\tfrac{D_{\mathrm{rad},R}}{\Delta x_R}$ for
interior faces,

$$
J_{ii} = \alpha_i + \mathrm{coeff}_L + \mathrm{coeff}_R,\qquad
J_{i,i-1} = -\mathrm{coeff}_L,\qquad
J_{i,i+1} = -\mathrm{coeff}_R .
$$

At a Marshak boundary the boundary flux contributes $\partial q_{\text{BC}}/\partial
G_N$ to the diagonal. The applied flux is the half-cell-corrected form
$q_{\text{BC}} = c_{\text{eff}}(G_{\text{ext}} + 4I_{\text{ext}} - G_N)$
(§8.6.1), so the assembled Jacobian uses the **same** effective coefficient as
the residual,

$$
c_{\text{eff}} = \frac{1}{1/c_M + \delta/D_{\mathrm{rad}}},\qquad
c_M(T)=\frac{\varepsilon(T)}{2(2-\varepsilon(T))},
$$

with $\delta$ the boundary-cell center-to-face distance, evaluated from the
boundary condition's emissivity law at the boundary-cell temperature
(`_marshak_effective_coefficient`, radiation_p1.jl):

$$
J_{NN}\mathrel{+}= \frac{A}{V_N}\,c_{\text{eff}}\quad\text{(top Marshak)},
\qquad
J_{11}\mathrel{+}= \frac{A}{V_1}\,c_{\text{eff}}\quad\text{(bottom Marshak)}.
$$

This keeps the Newton matrix exactly consistent with the Marshak residual for
non-blackbody boundaries. Temperature derivatives of emissivity are not included
because this Jacobian is with respect to the radiation unknowns $G$, while the
thermal state is held fixed during the quasi-steady P1 solve.

### 8.7.5 Allocation-free Thomas (TDMA) solve

The tridiagonal system $J\,\Delta G = -R$ is solved by the Thomas algorithm
(`thomas_solve!`, radiation_fluxes.jl) using pre-allocated diagonals
(`newton_dl`, `newton_d`, `newton_du`), a residual buffer `newton_R`, an update
buffer `newton_ΔG`, and a work array `newton_d_work` that preserves the original
diagonal across iterations (all on the `RadiationCache`). The solve is therefore
allocation-free and $O(n)$ per Newton step. Thomas is stable without pivoting for
the diagonally dominant P1 system.

### 8.7.6 The Newton loop, scaling, and clipping

Per iteration (radiation_fluxes.jl):

1. recompute the radiation fluxes from the current $G$ (`compute_radiation_fluxes!`);
2. recompute $Q_{\mathrm{em}}$, $Q_{\mathrm{abs}}$ (`compute_radiation_source_terms!`);
3. assemble $R$;
4. test convergence on the residual (below);
5. assemble the tridiagonal Jacobian;
6. solve $J\,\Delta G=-R$ via Thomas;
7. update $G_i \mathrel{+}= \text{damping}\cdot\Delta G_i$ and clip $G_i \leftarrow \max(G_i,0)$;
8. test convergence on the update.

**Defaults (the code is authoritative):** `tol = 1e-6`, `maxiter = 20`,
`damping = 1.0`, `initialize = true` (radiation_fluxes.jl). (A docstring
elsewhere says `tol=1e-10`; the executed signature uses `1e-6`.)

**Relative convergence.** $G$ is large ($10^4$–$10^6\ \mathrm{W\,m^{-2}}$), so an
absolute tolerance would be unphysical. The residual test scales by
$G_{\text{scale}}=\max(\|G\|_\infty, 1)$,

$$
\|R\|_\infty < \mathrm{tol}\cdot G_{\text{scale}}\ \Rightarrow\ \text{converged},
$$

and a second test on the relative update,
$\|\Delta G\|_\infty / \max(\|G\|_\infty,1) < \mathrm{tol}$, catches a
machine-zero plateau when $G$ is small. A warning is emitted only if the final
residual exceeds $100\,\mathrm{tol}\,G_{\text{scale}}$ (significant
non-convergence), avoiding noise on benign near-convergence.

**Non-negativity.** $G\ge0$ is enforced by `G[i] = max(G[i], 0.0)` after each
update, since the incident radiation cannot be negative.

### 8.7.7 Initialization

When $G$ is all-zero on entry (cold start) the field is initialized from local
thermodynamic equilibrium (`initialize_g_from_temperature!(rad_cache, mesh)`,
radiation_fluxes.jl):

$$
G_i = 4\sigma T_i^4 ,
$$

which makes $S_{\mathrm{rad}}=\alpha(G-4\sigma T^4)$ vanish initially — a good first
guess that avoids cold-start oscillation. (This same blackbody equilibrium is
exposed by `initialize_radiation_intensity!`,
src/Properties/radiation_cache.jl.)

---

## 8.8 Energy-Equation Coupling for Each Model

The condensed-phase energy equation (§3 Governing Equations and §7)
carries a volumetric radiation source $Q_{\text{rad}}$ on its right-hand side:

$$
\rho c_p^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\!\mathbf{q}_{\text{cond}} + S_{\text{conv}}
+ Q_{\text{rad}} + Q_{\text{rxn}} .
$$

The four models populate $Q_{\text{rad}}$ — or bypass it for the surface route — as
follows.

- **`NO_RADIATION`.** $Q_{\text{rad}}=0$. Surface re-radiation/incident handled only
  by the boundary balance (a `RadiativeBC`/`RadiativeFluxBC`).

- **`SURFACE_ABSORPTION`.** $Q_{\text{rad}}=0$ in the dispatched RHS
  (`compute_radiation_source!` zeroes it). The absorbed incident flux enters through
  the surface-temperature Newton solve as the boundary flux $q_{\text{BC}}=I_{\text{ext}}
  \alpha$ (boundary_temperature.jl); see §8.2 and §12.

- **`BEER_LAMBERT`.** $Q_{\text{rad},i}=q^{\text{rad}}_i$ from §8.3, added to the
  RHS in step 5 of `compute_rhs!` (src/Discretization/finite_volume.jl) and
  divided by $\rho c_p^{\text{eff}}$ in the final assembly (finite_volume.jl):

  $$
  \Big(\frac{\mathrm{d}T_i}{\mathrm{d}t}\Big)_{\text{rad}}
  = \frac{q^{\text{rad}}_i}{\rho c_p^{\text{eff}}_i}.
  $$

  The `RadiativeFluxBC` boundary flux is set to **zero** in Beer–Lambert mode
  (boundary_temperature.jl) so that the incident energy is not double-counted
  (it is already deposited volumetrically). The structured-operator wrapper
  `RadiationOperator` (src/Jacobian/operators/radiation.jl) routes the same physics
  through `compute_beer_lambert_source!`, treating `dT` as the `q_rad` buffer before
  the $\rho c_p$ divide.

- **`P1_QUASI_STEADY`.** Handled on the separate P1 residual
  (`compute_p1_residual!`, src/Experimental/RadiationP1.jl): the unified
  `residual!` runs first with `radiation_model = P1_QUASI_STEADY` from
  `def.options` (the `P1Residual` constructor enforces this), which zeroes the
  core volumetric source, the boundary `RadiativeFluxBC` **and** the boundary
  `RadiativeBC`; then the quasi-steady $G$ is solved and
  $S_{\mathrm{rad},i}=\alpha_i(G_i-4\sigma T_i^4)$ is added as
  $\mathrm{d}T_i/\mathrm{d}t\mathrel{+}=S_{\mathrm{rad},i}/\rho c_p^{\text{eff}}_i$.
  All radiation energy is owned by P1: the external flux enters through the
  Marshak BC (§8.6.1), emission and ambient exchange through the volumetric
  $4\sigma\alpha T^4$ term plus the Marshak exchange.

**The radiation-model BC dispatch is the linchpin.** The model-dependent
behavior of `RadiativeFluxBC` (boundary_temperature.jl, boundary_fluxes.jl) is
what keeps the incident energy single-counted across all four models:
surface-deposited for `NO_RADIATION`/`SURFACE_ABSORPTION`, volumetric (and
hence boundary-zero) for `BEER_LAMBERT`/`P1_QUASI_STEADY`. `RadiativeBC` is
additionally zeroed under `P1_QUASI_STEADY` only — under P1 the medium emits
volumetrically and exchanges with the ambient through the Marshak BC, so a
surface re-radiation law on top would double-count both channels; under
`BEER_LAMBERT` (which deposits only the *incident* flux volumetrically) the
surface emission stays at the boundary and `RadiativeBC` remains active. The
same dispatch threads through `CombinedThermalBC`, so a cone-calorimeter BC
such as
`RadiativeFluxBC(50e3) + RadiativeBC(ε=0.9,T_inf=300) + ConvectiveBC(h=10,T_inf=300)`
behaves correctly under any radiation model — under P1 the two radiative
members contribute exactly zero and only the convective member survives at
the boundary.

---

## 8.9 Regime Diagnostics

Three diagnostics quantify whether a chosen model is appropriate.

**Optical thickness.** `optical_thickness(mesh, prop_cache)`
(radiation_p1.jl) returns $\tau=\sum_i\alpha_i\Delta z_i=\int_0^L\alpha\,
\mathrm{d}z$, and `optical_thickness_cell` gives the per-cell
$\tau_i=\alpha_i\Delta z_i$ (a resolution diagnostic, not a validity
criterion). The standard validity bands (radiation_types.jl) are:

| $\tau=\alpha L$ | regime | recommendation |
|---|---|---|
| $\tau>10$ | optically thick | P1 excellent; surface absorption acceptable |
| $1<\tau<10$ | intermediate | P1 good; Beer–Lambert recommended |
| $0.1<\tau<1$ | optically thin | Beer–Lambert recommended |
| $\tau<0.1$ | transparent | `NO_RADIATION` |

**P1 validity.** `check_p1_validity(mesh, prop_cache)` (radiation_p1.jl)
classifies the regime by the **domain** optical thickness and warns when P1 is
marginal (15–30 % error for $0.1<\tau<1$) or invalid ($\tau<0.1$, $>30$ % error —
consider discrete ordinates). Validity is a property of the domain $\tau$
only: a per-cell criterion ($\tau_i \ge 1$ in every cell) would be
mesh-dependent — refining the mesh shrinks every $\tau_i$ with no change in the
physics — so no per-cell validity test is applied.

**Radiation energy balance.** `radiation_energy_balance`
(radiation_p1.jl) integrates $\int Q_{\mathrm{em}}\,\mathrm{d}V$ and
$\int Q_{\mathrm{abs}}\,\mathrm{d}V$, reports the net source and the top/bottom
boundary fluxes, and forms the integral residual of the quasi-steady G
equation, $(q_{\text{top}}-q_{\text{bottom}}) + (\text{net source})$, which
vanishes at a converged quasi-steady state (the boundary-flux difference and
the net source *cancel* — they do not equal each other: `q_rad` stores the
physical $+z$ flux at every face, and $\int(Q_{\mathrm{abs}}-Q_{\mathrm{em}})
\,\mathrm{d}V + \oint \mathbf{q}\cdot\hat{\mathbf{n}}\,\mathrm{d}A = 0$). A
speed-of-light radiation energy density $E_{\text{rad}}=G/c$
($c=2.998\times10^8\ \mathrm{m\,s^{-1}}$) is also available as a diagnostic
(radiation_cache.jl).

---

## 8.10 The Beer–Lambert Non-Local Jacobian

### 8.10.1 Why it is upper-triangular

Beer–Lambert makes $q^{\text{rad}}_i$ depend on the absorption coefficient of every
cell **above** it, through the cumulative attenuation $I_{\text{in}}(i)=
I_{\text{surface}}\exp(-\sum_{j>i}\alpha_j\Delta z_j)$. Differentiating with respect
to a state variable $x_j$ (either $T_j$ or $\xi_{k,j}$) gives a structurally
**upper-triangular** coupling (src/Jacobian/operators/radiation.jl):

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_j} =
\begin{cases}
\big(I_{\text{in}}(i) - q^{\text{rad}}_i\,\Delta z_i\big)\,\dfrac{\partial \alpha_i}{\partial x_i}, & j=i,\\[2.0ex]
-\,q^{\text{rad}}_i\,\Delta z_j\,\dfrac{\partial \alpha_j}{\partial x_j}, & j>i,\\[1.0ex]
0, & j<i.
\end{cases}
$$

The **diagonal** ($j=i$) term has a clean identity. From
$q^{\text{rad}}_i=I_{\text{in}}(1-e^{-\alpha_i\Delta z_i})/\Delta z_i$,

$$
\frac{\partial q^{\text{rad}}_i}{\partial\alpha_i}
= I_{\text{in}}\,e^{-\alpha_i\Delta z_i}
= I_{\text{in}} - q^{\text{rad}}_i\,\Delta z_i ,
$$

so the diagonal coefficient is exactly $\big(I_{\text{in}}-q^{\text{rad}}_i\Delta z_i\big)\partial\alpha_i/\partial x_i$.
(The form $I_{\text{in}}(1-(1-e^{-\tau})/\tau)$ is sometimes quoted; it is *not*
equal to this and is incorrect — the code form above is the verified derivative.)

The **off-diagonal** ($j>i$) term comes from
$\partial I_{\text{in}}(i)/\partial\alpha_j = -I_{\text{in}}(i)\,\Delta z_j$ for
$j>i$, giving $\partial q^{\text{rad}}_i/\partial x_j = -q^{\text{rad}}_i\,\Delta z_j\,
\partial\alpha_j/\partial x_j$.

**The Kirchhoff surface-transmissivity chain.** On the default Kirchhoff path
(§8.3.2) the injected intensity is
$\varepsilon_{\text{eff}}(T_n,\xi_{\cdot,n})\,I_{\text{surface}}$, so every
$q^{\text{rad}}_i \propto \varepsilon_{\text{eff}}$ and each energy row gains a
dependence on the **surface cell's** state:

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_n} \mathrel{+}=
q^{\text{rad}}_i\,\frac{1}{\varepsilon_{\text{eff}}}
\frac{\partial\varepsilon_{\text{eff}}}{\partial x_n}\qquad(i<n).
$$

In the structured factorization this folds into the same suffix coupling via
$v_n \mathrel{-}= (\partial\varepsilon_{\text{eff}}/\partial x_n)/\varepsilon_{\text{eff}}$;
the surface cell's own diagonal closure recomputes
$\varepsilon_{\text{eff}}(T,\xi)$ live so its
$\partial\varepsilon_{\text{eff}}/\partial(T_n,\xi_{\cdot,n})$ enters exactly;
and the sparsity pattern gains $T$-row($i$) × $T$-column($n$) entries for $i<n$
(the $\xi_{\cdot,n}$ columns already exist in the suffix block;
`_emit_radiation_pattern!`, src/Jacobian/sparsity.jl). With an explicit BC
absorptivity (`radiation_kirchhoff = false`) the scaling is folded into
$I_{\text{surface}}$ upstream and none of these terms appear. The chain is
verified entry-by-entry against the dense AD reference (§18.3).

### 8.10.2 Structured factorization (suffix-sum coupling)

The $j>i$ block factors as $\mathrm{Diag}(-q^{\text{rad}})\cdot U\cdot
\mathrm{diag}(v)$ where $v_j=\Delta z_j\,\partial\alpha_j/\partial x_j -
\delta_{j,n}(\partial\varepsilon_{\text{eff}}/\partial x_n)/\varepsilon_{\text{eff}}$
(the last term on the Kirchhoff path only) and $U$ is
the strict upper-triangular all-ones matrix — i.e. a **suffix scan**
$\sum_{k>i}v_k$. The operator pushes it as an `EmbeddedScaledSuffixSumCoupling`
per channel ($T$ and each $\xi_j$) with per-row scale
$-q^{\text{rad}}_i/\rho c_p^{\text{eff}}_i$ (radiation.jl); the coupling stores
$O(n)$ data and is written into the sparse Jacobian's upper-triangular block by
`densify!` under the exact `:full` profile (§15.9.2). The diagonal term — including
the $\rho c_p$ chain rule $-q^{\text{rad}}/(\rho c_p)^2\,\partial\rho c_p/\partial x_i$
— is captured by a cell-local ForwardDiff over the closure
$(T,\xi)\mapsto I_{\text{in}}(1-e^{-\alpha\Delta z})/(\Delta z\,\rho c_p^{\text{eff}})$
with $I_{\text{in}}$ held fixed (radiation.jl). For the standard linear
mixing $\alpha=\sum_j\xi_j\alpha_j$ the derivatives reduce to $\partial\alpha/\partial
T=0$ and $\partial\alpha/\partial\xi_j=\alpha_j$, so the suffix coupling is cheap;
the AD path keeps it generic for future $T$-dependent absorption laws. (See §15,
Temporal Integration & the Structured Jacobian, for how the embedded couplings are
assembled.)

**ALE geometry columns.** Under ALE the deposit also depends on the live geometry
($\Delta z_i = V_i/A$) through the per-cell factor $(1-e^{-\alpha_i\Delta z_i})/\Delta z_i$
and the transmission suffix $\exp(-\sum_{j>i}\alpha_j\Delta z_j)$. The radiation operator
fills the $T$ and $\xi$ columns only; the $z$ columns are written by the
`GeometryOperator`'s AD pass over the full residual under the exact `:full` profile
(§15.11.3).

> **Solver impact.** Under `BEER_LAMBERT` the sparsity pattern is upper-triangular rather
> than tridiagonal; the other three models keep the tridiagonal pattern. Under the
> `:ale_reduced` assembly profile (`ApproxSolve`) the suffix couplings and every $z$ column
> are dropped, trading Jacobian inexactness for speed (§15.10).

---

## 8.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to ours on first use. All three reference
solvers are condensed-phase pyrolysis models; their radiation treatments span the
same idioms Pyrolysis.jl supports, plus a two-flux variant in FDS.

### 8.11.1 ThermaKin2Ds

ThermaKin (Stoliarov & Lyon) uses **in-depth exponential (Beer–Lambert)
absorption** of the external flux together with explicit surface **re-radiation**.
In ThermaKin notation the external radiant flux is $\dot q''_{ex}$ (our
$I_{\text{ext}}$), the absorption coefficient is per-component and forms a
volumetric coefficient (our $\alpha_{\text{eff}}=\sum_j\xi_j\alpha_j$), and the
in-depth source is $\dot q''_{ex}\,\alpha\,e^{-\int\alpha\,\mathrm{d}x}$ — identical
in form to our $q'''=\alpha I$ with $I=I_{\text{surface}}e^{-\int\alpha\,\mathrm{d}z}$
(our `BEER_LAMBERT`). ThermaKin additionally models re-radiation $I_{rr}=
\varepsilon\sigma T_s^4$ at the surface (often with reflection of a fraction of
$\dot q''_{ex}$), which Pyrolysis.jl supplies via a top `RadiativeBC` re-radiation
term $F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$ in the boundary balance
(§12). ThermaKin also offers a radiative-conductivity augmentation $k\propto T^3$
to mimic radiative diffusion in char — the physical content of our P1 diffusion
term $\nabla\!\cdot\!(D_{\mathrm{rad}}\nabla G)$, which Pyrolysis.jl models
explicitly through `P1_QUASI_STEADY` rather than by inflating $k$. ThermaKin
publishes heats of reaction with the opposite sign convention ($h>0$ exothermic);
this does not touch the radiation terms (overload note H1, §2, Nomenclature).

### 8.11.2 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) deposits the external radiant flux
**in-depth** with an exponential absorption profile and treats the surface energy
balance with combined radiative–convective exchange. In Gpyro notation the
pre-exponential of any property is $Z$ (our $A_i$), the volume fraction is $X_\alpha$
(our $v_j$), and the in-depth absorbed flux follows $\dot q''_e\,\kappa\,e^{-\kappa
x}$ — again our `BEER_LAMBERT` form with Gpyro's $\kappa$ (the absorption
coefficient, distinct from the permeability $K$ in §9) mapping to our volumetric
$\alpha_{\text{eff}}$. Gpyro's effective radiative conductivity in optically thick
char (a temperature-cubed conductivity contribution) plays the same role as our P1
diffusion term; Pyrolysis.jl keeps the two mechanisms separate (a genuine $k(T,\xi)$
mixture conductivity in §5/§7 and an explicit $G$ field in P1). Gpyro's surface
emissivity enters the boundary balance just as our $\varepsilon_{\text{eff}}=
\sum_j v_j^{\text{solid}}\varepsilon_j$ does.

### 8.11.3 FDS

The FDS solid-phase model uses a **two-flux (Schuster–Schwarzschild)** in-depth
radiation model for the condensed phase: forward and backward hemispheric
intensities $I^+,I^-$ are propagated with absorption coefficient $\kappa_s$ (our
$\alpha$), and the radiative source in the 1D Fourier conduction equation is the
divergence of the net radiative flux $I^+-I^-$. Pyrolysis.jl's `BEER_LAMBERT`
corresponds to the single-flux limit (only the inward-going beam, no internal
re-emission), while `P1_QUASI_STEADY` is the diffusion-limit analog of the two-flux
model (P1 and two-flux/Schuster–Schwarzschild are the two classical first-order
closures of the RTE; both reduce to a single second-order ODE for a combined field).
FDS applies surface absorptivity/emissivity through Kirchhoff's law at the gas–solid
interface — the same default ($\alpha=\varepsilon_{\text{eff}}$) Pyrolysis.jl uses
for `RadiativeFluxBC`. FDS's gray-band absorption maps to our gray (wavelength-
independent) $\alpha_{\text{eff}}$. The principal difference is dimensionality and
the explicit two-stream bookkeeping: FDS resolves $I^\pm$ separately, whereas
Pyrolysis.jl resolves either the inward beam ($I$, Beer–Lambert) or the isotropic
incident field ($G$, P1).

### 8.11.4 Summary table

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS (solid) |
|---|---|---|---|---|
| In-depth absorption | Beer–Lambert $q'''=\alpha I$ | Beer–Lambert | exp. in-depth | two-flux $I^\pm$ |
| Re-radiation diffusion | explicit P1 ($G$) | radiative $k\propto T^3$ | radiative $k$ | two-flux divergence |
| Surface absorptivity | $\varepsilon_{\text{eff}}$ (Kirchhoff) | $\varepsilon$ + reflection | $\varepsilon$ | Kirchhoff |
| Spectral model | gray | gray | gray | gray-band |
| Surface re-radiation | `RadiativeBC` $\varepsilon\sigma(T_\infty^4-T_s^4)$ | $\varepsilon\sigma T_s^4$ | radiative–convective | $\varepsilon\sigma T_s^4$ |

---

## 8.12 Experimental Status of P1: the Separate Residual

P1 is **Experimental**: it is research-grade and out of scope for the unified
trait-dispatched residual that drives `SURFACE_ABSORPTION`/`BEER_LAMBERT` runs
(`RadiationP1.jl` header). `solve` rejects `radiation_model = P1_QUASI_STEADY`;
the P1 kernels and the residual wrapper below are called directly. Three structural
consequences follow.

**No dedicated state layout.** $G$ is an algebraic (quasi-steady) field, not a
state variable: it lives in the `RadiationCache` and is re-solved from the
current $(T,\xi,z)$ state at **every** residual evaluation. The ODE state
vector is the ordinary unified `StateLayout` (T, ξ, z when ALE, χ when
WithChi); build `def`, `ws`, and `u0` with
`initialize_state(problem; radiation_model = P1_QUASI_STEADY)` (which applies
the problem's initial conditions before packing — a bare
`create_state_vector` call packs the mesh as-is, §15.2) and read
the field with `p.rad_cache.G` or `get_radiation_intensity(p::P1Residual, i)`
— always the field of the most recent residual evaluation, never a frozen
copy.

**A thin dedicated residual.** `compute_p1_residual!`
(src/Experimental/RadiationP1.jl) delegates to the unified `residual!` —
conduction, reactions, transport, ALE, χ all come from the production code
path — and then (1) mirrors `u` into the `Mesh1D`, which the
Float64-only P1 kernels read, (2) updates $\alpha$/$D_{\mathrm{rad}}$ and solves
the quasi-steady $G$ system, and (3) adds
$S_{\mathrm{rad}}/\rho c_p = \alpha(G-4\sigma T^4)/\rho c_p$ to the energy
rows. The `P1Residual` constructor **requires**
`def.options.radiation_model == P1_QUASI_STEADY` (build the def with
`to_problem_def(problem; radiation_model = P1_QUASI_STEADY)`): the unified
residual then runs with zero core volumetric radiation source and zero
boundary `RadiativeFluxBC`/`RadiativeBC` flux, so the incident energy and the
surface emission are counted exactly once — through the Marshak BC and the
volumetric $4\sigma\alpha T^4$ emission (§8.8).

**Float64-only.** The `RadiationCache` arrays are `Float64`, and the P1 Newton
solve runs in Float64; AD `Dual` parts of the time variable are stripped via
`ForwardDiff.value` when storing boundary fluxes (radiation_p1.jl), and
state-vector `Dual`s are not supported. Consequently P1 sensitivities through
the $G$ field are not propagated by forward-mode AD the way the core models'
Jacobians are; P1 stays here until the trait system gains a first-class
`P1Radiation` mode (§15, §17 Sensitivity Analysis).

---

## 8.13 Limitations and Open Issues

The radiation subsystem makes several deliberate simplifications; the following are
the open issues most likely to matter in practice.

1. **Scattering is not exercised.** $\sigma_s$ is plumbed through
   $D_{\mathrm{rad}}=1/(3(\alpha+\sigma_s))$ but is taken $=0$ everywhere for
   pyrolysis. Strongly scattering chars or foams would need $\sigma_s>0$ and, for
   anisotropy, a higher-order closure (P3 or discrete ordinates).

2. **Gray (wavelength-independent) properties.** $\alpha_j$ and $\varepsilon_j$ are
   gray; real polymers and chars are spectrally selective. A gray-to-band extension
   would require per-band absorption coefficients and a banded $G$ field.

3. **P1 is quasi-steady only.** $\partial G/\partial t=0$ is built in: $G$ is
   an algebraic field in the `RadiationCache`, re-solved at every residual
   evaluation, and is not part of the ODE state at all. A fully transient P1
   would integrate $G$ as a genuine state variable with its own time scale —
   rarely needed given $\tau_{\mathrm{rad}}\ll\tau_{\mathrm{cond}}$, but
   unavailable today.

4. **The surface energy-balance Newton Jacobian neglects $\mathrm{d}\varepsilon/\mathrm{d}T$**
   (§8.8, §12.2). The Marshak coefficient itself is consistent — the half-cell
   effective coefficient $c_{\text{eff}} = 1/(1/c_M + \delta/D_{\mathrm{rad}})$,
   $c_M = \varepsilon/(2(2-\varepsilon))$, appears identically in the residual
   and in the Jacobian (§8.7.4). The remaining approximation is that,
   where a state-dependent emissivity is used, the temperature dependence of
   $\varepsilon$ is dropped from the surface-temperature derivative. For the
   near-constant emissivities used in practice this is negligible; a strongly
   $T$-dependent $\varepsilon$ can slow Newton convergence.

5. **No interior emissivity factor by design** (§8.5). This is physically correct
   for a gray body in LTE, but it means surface emissivity influences the interior
   field *only* through the Marshak boundary; users expecting a volumetric
   $\varepsilon\sigma T^4$ emission will not find one.

6. **P1 is AD-opaque (Float64-only).** Sensitivities through the $G$ field are not
   propagated by the forward-mode AD path used elsewhere; P1 remains Experimental and
   outside the unified trait dispatch (§8.12, §15, §17).

7. **`SURFACE_ABSORPTION` assumes an opaque first cell** (§8.2). Absorbed
   incident radiation enters through the *boundary balance*; the model is
   accurate only when the first cell is optically thick ($\tau_N\gg1$). If the
   first cell is coarse or only marginally opaque, refine the surface cell or
   switch to `BEER_LAMBERT`.

8. **Transparent-limit clamp on $D_{\mathrm{rad}}$.** Below
   $10^{-10}\ \mathrm{m^{-1}}$ extinction,
   `compute_radiation_diffusion_coefficient` — the single canonical guard, used
   by the one `update_radiation_properties!` method — clamps the diffusion
   coefficient to the large finite value $1/(3\cdot 10^{-10})$ rather than
   letting it diverge (or zeroing it, which would make the quasi-steady Newton
   system singular wherever $\alpha = 0$). A transparent cell thus behaves as a
   near-perfect conductor of $G$. P1 is not recommended for $\tau<0.1$
   regardless.


---

# 9. Gas-Phase Mass Transport and Pressure

This chapter derives the constitutive and discrete laws by which gaseous
decomposition products move through the condensed-phase matrix in Pyrolysis.jl.
Two transport mechanisms are treated: ordinary (Fickian) diffusion driven by a
volume-fraction / partial-pressure gradient, and Darcy advection driven by the
gas-phase pressure field in a porous medium; together they form the *Darcy–Fick*
flux. We establish the ideal-gas pressure closure in pore space, the Darcy
velocity and its analytic pressure derivatives, the surface mass-flux boundary
treatments, the gas-only finite-volume divergence that closes the species
equations, and the analytic Jacobian contributions of every flux term. Every
equation in this chapter has been checked against the implementing source and,
where the textbook form and the code differ, the chapter documents what the code
actually does.

Throughout, the spatial coordinate is the through-thickness axis $z$, with
$z=0$ the substrate (bottom, boundary id 2) and $z=L$ the exposed surface (top,
boundary id 1). The face-flux sign convention is **positive = transport in the
$+z$ direction** (bottom $\to$ top), and the discrete divergence
$(F_R-F_L)A/V$ is positive when the quantity flows *out* of a cell. In this
chapter the index $i$ denotes a **cell**, $j$ a **component/species**, and $f$ a
**face**; $L,R$ denote the left/right cells bounding a face. The primary species
state is the mass concentration $\xi_j$ [kg·m⁻³], not the mass fraction
$Y_j=\xi_j/\rho$.

---

## 9.1 The volume-fraction-gradient diffusive flux

### 9.1.1 Continuous form and physical origin

Gaseous decomposition products in a hot porous solid diffuse down their
concentration gradient. The natural driving potential for an ideal gas is not the
mass concentration $\xi_j$ but the *mass fraction* $Y_j=\xi_j/\rho_g$, because
Fick's first law for a dilute species in a carrier of density $\rho_g$ is

$$
J_j = -\rho_g\,\lambda\,\frac{\partial Y_j}{\partial z}
    = -\rho_g\,\lambda\,\frac{\partial}{\partial z}\!\left(\frac{\xi_j}{\rho_g}\right),
$$

where $\lambda$ [m²·s⁻¹] is the effective gas-transfer (diffusion) coefficient of
the porous medium (a tortuosity- and porosity-corrected molecular diffusivity;
see §5 *Effective Properties and Mixing Rules*, and §9.2 below). For an ideal gas
in the pore space, the gas density at fixed pressure scales as

$$
\rho_g \propto \frac{P}{T},
$$

so along an isobar $\rho_g \propto 1/T$ and the volume fraction (equivalently, the
mass-fraction-like ratio) of species $j$ obeys $\xi_j/\rho_g \propto \xi_j T$. Up
to the pressure scale, then,

$$
J_j = -\rho_g\,\lambda\,\frac{\partial}{\partial z}\!\left(\frac{\xi_j}{\rho_g}\right)
    = -\frac{\lambda}{T}\,\frac{\partial(\xi_j T)}{\partial z}.
$$

This is the **volume-fraction-gradient formulation** documented at the head of
`src/Physics/gas_transport.jl`. Expanding the product derivative,

$$
\boxed{\;
J_j = -\lambda\left(\frac{\partial \xi_j}{\partial z}
       + \frac{\xi_j}{T}\,\frac{\partial T}{\partial z}\right)
\;}
$$

exposes the two physical contributions:

- **Ordinary diffusion** $-\lambda\,\partial\xi_j/\partial z$: transport down the
  composition gradient.
- **Thermal coupling** $-\lambda(\xi_j/T)\,\partial T/\partial z$: a Soret-like
  term arising *not* from genuine thermodiffusion but from the temperature
  dependence of the ideal-gas density. A hotter cell holds the same mass
  concentration in a more dilute gas (lower $\rho_g$), so equal $\xi_j$ across a
  temperature gradient still corresponds to a volume-fraction gradient and hence a
  flux.

Casting the driving potential in terms of $\xi_j T$ rather than $Y_j$ avoids ever
computing the gas density $\rho_g$ at the face, which would require evaluating $P$
and $T$ at an intermediate point. The formulation is well behaved in the
low-density (high-temperature, low-$P$) regime where $\rho_g\to0$ would otherwise
appear in a denominator.

### 9.1.2 Discrete face flux

Discretizing $-\frac{\lambda}{T}\,\partial(\xi T)/\partial z$ at an interior face
$f$ between cells $L$ and $R$, with the temperature in the prefactor taken at the
arithmetic-mean face value $T_{\text{face}}=\tfrac12(T_L+T_R)$, gives the
implemented kernel:

$$
\boxed{\;
J_j = -\,\lambda_{\text{face}}\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\;d_{LR}}
\;}
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

where $d_{LR}=|z_R^c - z_L^c|$ is the distance between the two cell centers and
$\lambda_{\text{face}}$ is the distance-weighted harmonic mean of the cell
transfer coefficients (§9.1.4). This is exactly

```julia
@inline function gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)
    T_face = (T_L + T_R) / 2
    return -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
end
```

*Implementation note.* Implemented in `gas_flux`, `src/Physics/gas_transport.jl`,
and assembled face-by-face by `compute_gas_fluxes!` (with the
compile-time-unrolled inner helper `compute_gas_fluxes_unrolled!`).

### 9.1.3 The exact product-rule identity

The composition / thermal split of §9.1.1 is not an approximation of the discrete
kernel — it is algebraically identical to it when arithmetic-mean face values are
used. The identity

$$
\xi_R T_R - \xi_L T_L
   = T_{\text{face}}\,(\xi_R - \xi_L) + \xi_{\text{face}}\,(T_R - T_L),
\qquad \xi_{\text{face}}=\tfrac12(\xi_L+\xi_R),
$$

is proved by expanding both terms on the right:

$$
\begin{aligned}
T_{\text{face}}(\xi_R-\xi_L) + \xi_{\text{face}}(T_R-T_L)
&= \tfrac12(T_L+T_R)(\xi_R-\xi_L) + \tfrac12(\xi_L+\xi_R)(T_R-T_L) \\
&= \tfrac12\big(2\,\xi_R T_R - 2\,\xi_L T_L\big)
 = \xi_R T_R - \xi_L T_L .
\end{aligned}
$$

Dividing by $T_{\text{face}}\,d_{LR}$ and multiplying by $-\lambda_{\text{face}}$
recovers the discrete composition + thermal split

$$
J_j = -\lambda_{\text{face}}\!\left(\frac{\xi_R-\xi_L}{d_{LR}}
       + \frac{\xi_{\text{face}}}{T_{\text{face}}}\frac{T_R-T_L}{d_{LR}}\right).
$$

This equivalence is verified to machine precision (`atol = 1e-14`) in
`test/discretization/test_gas_flux.jl`.

### 9.1.4 Why the harmonic mean is applied to $\lambda$, not $\lambda\rho$

At a material interface, *flux continuity* — not property continuity — is the
physical constraint. Treating the two half-cells as transport resistances in
series, the face transfer coefficient is the distance-weighted harmonic mean

$$
\boxed{\;
\lambda_{\text{face}} = \frac{d_L + d_R}{\dfrac{d_L}{\lambda_L} + \dfrac{d_R}{\lambda_R}}
\;}
$$

where $d_L,d_R$ are the distances from the two cell centers to the face. On a
uniform mesh ($d_L=d_R$) this reduces to the standard harmonic mean
$2\lambda_L\lambda_R/(\lambda_L+\lambda_R)$. The distance weighting makes the
scheme correct on stretched meshes (e.g. surface-refined grids); see §13.7.1 for
the general resistance-in-series argument and §7.1 for the identical treatment of
thermal conductivity. A small-value guard returns $0$ if either $\lambda_L$ or
$\lambda_R$ is below $10^{-20}$.

A subtle but important modeling choice: the code takes the harmonic mean of
$\lambda$ **alone**, not of the product $\lambda\rho_g$ that would appear in a
true mass-fraction-gradient discretization $J=-\lambda\rho_g\,\partial Y/\partial z$.
This is consistent with — indeed *required by* — the volume-fraction formulation,
in which $\rho_g$ never appears explicitly: the density dependence is folded into
the $\xi T$ driving potential, which uses *cell-center* concentrations $\xi_L,\xi_R$
and the arithmetic-mean face temperature, with no face density. No
$\lambda\rho$-averaging helper exists in the package; the distinction is
exercised in `test/discretization/test_gas_flux.jl`.

*Implementation note.* `λ_face` is precomputed once per step in the
`PropertyCache` via `harmonic_mean_weighted`,
`src/Properties/effective_properties.jl`; per-cell $\lambda_i$
comes from `effective_gas_transfer` (§9.2 and §5).

### 9.1.5 Partial-pressure equivalence

Because the formulation is anchored to the ideal-gas relation, the same flux can
be written as a partial-pressure gradient. With $P_j = \xi_j R_g T / M_j$ the
partial pressure of species $j$,

$$
J_j = -\frac{\lambda_{\text{face}}\,M_j}{R_g\,T_{\text{face}}}\,\frac{\partial P_j}{\partial z}.
$$

To see the equivalence, note $\xi_j T = P_j M_j / R_g$, so
$\xi_R T_R - \xi_L T_L = (M_j/R_g)(P_{j,R}-P_{j,L})$, and substituting into the
discrete kernel gives the boxed pressure form directly. This is verified to
`atol = 1e-12` in `test/discretization/test_gas_flux.jl`. The
partial-pressure reading reinforces the physical picture: gas diffuses down its
own partial-pressure gradient, with mobility $\lambda M_j/(R_g T)$.

---

## 9.2 One transfer coefficient for all species

Every gas component diffuses with the **same** effective transfer coefficient
$\lambda_{\text{eff}}$: the face coefficient $\lambda_{\text{face}}$ entering
the kernel of §9.1.2 is a single mixture property (§5.6), shared by all species
at that face. The model therefore does not discriminate between light and heavy
gases within one simulation — a light product (H₂O) and a heavy one (MMA)
diffuse at the same rate through the same matrix. The per-component `λ_j(T)`
property controls how each *condensed/gas phase* contributes to the mixture
coefficient as composition evolves, not how individual gas species move
relative to one another.

Where a species-faithful magnitude is wanted, the Chapman–Enskog helper
`chapman_enskog_diffusivity(T, M; P, σ, ε_k)` (§5.10) estimates the binary
diffusivity of a given vapor in air; the user can fold that estimate (with a
porosity/tortuosity correction) into the component `λ` values. Genuine
multicomponent (Stefan–Maxwell) diffusion, or even a per-species kinetic-theory
scaling of the diffusive flux, is not implemented (§9.12).

---

## 9.3 The Darcy–Fick combined flux

When the material has nonzero permeability, pressure-driven (Darcy) advection
augments diffusion. The total face flux is the sum of an advective and a diffusive
part:

$$
\boxed{\;
J_j^{\text{DF}} = \underbrace{\frac{\xi_j^{\text{upwind}}}{\varphi_{\text{face}}}\,v_g}_{J_j^{\text{adv}}}
   \;+\; \underbrace{\Big(\!-\lambda_{\text{face}}\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d_{LR}}\Big)}_{J_j^{\text{diff}}}
\;}
$$

where $v_g$ is the Darcy velocity at the face (§9.5, positive in $+z$). The
$1/\varphi_{\text{face}}$ is forced by the bookkeeping: $\xi_j$ is mass per
**bulk** volume (the pressure closure $P=\sum\xi_jR_gT/(M_j\varphi)$ divides by
$\varphi$ to recover the intrinsic pore concentration) while the Darcy velocity
is **superficial** (volumetric flux per total cross-section), so the mass flux
per bulk area is intrinsic density × superficial velocity,
$J^{\text{adv}}=(\xi/\varphi)\,v_g$ — equivalently $\xi\,v_i$ with the
interstitial velocity $v_i=v_g/\varphi$. The face porosity is the arithmetic
mean $\varphi_{\text{face}}=\tfrac12(\varphi_L+\varphi_R)$, floored at
$10^{-12}$ consistently with the pressure closure. The discrete kernel is

```julia
@inline function gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)
    T_face = (T_L + T_R) / 2
    J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
    φ_face = max((φ_L + φ_R) / 2, _PHI_FLOOR)
    ξ_upwind = v_g >= 0 ? ξ_L : ξ_R
    J_adv = v_g * ξ_upwind / φ_face
    return J_adv + J_diff
end
```

### 9.3.1 First-order upwind selection

The advected concentration is taken from the **upwind** cell — the cell from which
the flow originates:

$$
\xi_j^{\text{upwind}} = \begin{cases}
\xi_{j,L}, & v_g \ge 0 \quad (\text{flow toward }+z), \\
\xi_{j,R}, & v_g < 0 \quad (\text{flow toward }-z).
\end{cases}
$$

This first-order upwind rule is unconditionally stable and conservative (the same
face value enters the two cells sharing the face with opposite divergence signs).
It is total-variation-diminishing but only first-order accurate; no TVD/MUSCL
slope limiting is applied to the Eulerian Darcy–Fick advective term. (Higher-order
limited advection exists only for the mesh-relative temperature advection of the
`w_override` verification hook; see §11.6,
*The ALE Framework*, and §13.4–13.5 for the limiter library.) A practical
consequence: when upwinding selects a cell with zero concentration adjacent to a
cell with positive concentration, the advective contribution can momentarily
produce a "wrong-side" zero flux; this is conservative but, over many steps,
relies on reaction and heat coupling to keep states physical.

### 9.3.2 Advection carries all species equally

The Darcy–Fick flux uses the bulk velocity $v_g$ in the advective term for every
component, and the diffusive term likewise uses one shared $\lambda_{\text{face}}$
for all species (§9.2). Physically, bulk pore flow is a mixture motion; the model
applies no molecular-weight discrimination in either term (§9.2, §9.12).

In the pure-Fickian limit $v_g\to0$, the Darcy–Fick flux reduces exactly to the
diffusion-only flux of §9.1 — the two assembly paths agree, which is checked in
the gas-flux test suite.

*Implementation note.* `gas_flux_darcy_fick`,
`src/Physics/gas_transport.jl`; assembled by
`compute_gas_fluxes_darcy_fick!`.

---

## 9.4 Ideal-gas pressure in a porous medium

The pressure that drives Darcy flow is the gas-phase pressure inside the pore
space. Summing partial pressures of the gaseous components and dividing by the
porosity (the fraction of bulk volume actually occupied by gas) gives the
porous-medium ideal-gas law:

$$
\boxed{\;
P = \frac{1}{\varphi}\sum_{j:\,\text{is\_gas}(j)} \frac{\xi_j R_g T}{M_j}
  \equiv \frac{N(\xi,T)}{\varphi}
\;},
\qquad N \equiv \sum_{j:\,\text{gas}} \frac{\xi_j R_g T}{M_j}.
$$

The numerator $N$ is the sum of partial pressures *referred to the bulk volume*
(mass per bulk volume times $R_g T/M_j$). Dividing by the porosity $\varphi$
converts to the pressure *in the gas phase*, which occupies only the fraction
$\varphi$ of the bulk: the same gas mass confined to a smaller void volume exerts
a proportionally higher pressure. Without the $1/\varphi$ divisor the pressure
would be underestimated in nearly-solid cells.

Two implementation details are load-bearing:

- **Only true gases contribute to $N$.** A component contributes iff
  `is_gas(j)` is true *and* `molar_mass(j) > 0`; solids and liquids are skipped.
  This is a compile-time-eliminated branch in the generated sum
  `_gas_partial_pressure_sum`.
- **Porosity floor.** The divisor is clamped from below at
  $\varphi_{\min}=10^{-12}$ via `max(φ, _PHI_FLOOR)`, guaranteeing a finite
  pressure even if a cell becomes fully solid ($\varphi\to0$). This matches the
  clamping in the porosity routine (porosity is itself clamped to $[0,1]$).

The porosity $\varphi$ is computed from *intrinsic* (skeletal) component densities
$\rho_{i,j}$, not bulk densities:
$\varphi = \mathrm{clamp}\!\big(1-\sum_{j\in\text{S/L}}\xi_j/\rho_{i,j},\,0,\,1\big)$
(see §5.3). It is the same porosity used by the Carman–Kozeny permeability closure
(§9.5.2), so the pressure and the permeability share a single consistent void
fraction.

*Implementation note.* `ideal_gas_pressure`, `src/Physics/pressure.jl`; the cell
field is assembled in `update_darcy_properties!` (`src/Discretization/finite_volume.jl`).
Constant $R_g = 8.314462618$ J·mol⁻¹·K⁻¹.

---

## 9.5 Pressure gradient and Darcy velocity

### 9.5.1 Pressure gradient

At an interior face the pressure gradient is a plain central difference of
cell-centered pressures:

$$
(\nabla P)_f = \frac{P_R - P_L}{d_{LR}}.
$$

Boundary faces (one neighbor id $=0$) are assigned a zero gradient here; the
physical pressure boundary conditions are applied separately to the Darcy velocity
(§9.5.3, §12.5). The sign convention is the usual one: $(\nabla P)_f>0$ when
pressure increases toward $+z$.

*Implementation note.* The face gradient is formed in `update_darcy_properties!`
(`src/Discretization/finite_volume.jl`), together with the Darcy velocity.

### 9.5.2 Darcy velocity

Darcy's law gives the superficial (bulk) gas velocity at a face from the pressure
gradient, the permeability, and the gas viscosity:

$$
\boxed{\;
v_g = -\frac{\kappa_{\text{face}}}{\mu_{\text{face}}}\,(\nabla P)_f
    = -\frac{\kappa_{\text{face}}}{\mu_{\text{face}}}\,\frac{P_R - P_L}{d_{LR}}
\;}
$$

with $v_g>0$ for flow toward $+z$ (the surface). The face properties combine the
two cells consistently with their physical roles:

- **Mobility — distance-weighted harmonic mean.** Permeability is a flow
  *resistance* in series across the interface. The production path
  (`update_darcy_properties!`) forms the per-cell mobility
  $\text{mob}=\kappa/\mu$ and combines the two sides with the
  distance-weighted harmonic mean

  $$
  \text{mob}_{\text{face}} = \frac{d_L+d_R}{d_L/\text{mob}_L + d_R/\text{mob}_R},
  $$

  which reduces to the unweighted harmonic mean $2ab/(a+b)$ on uniform meshes.
  The harmonic mean is conservative: it tends toward the *smaller* of the two
  mobilities, so a single tight (low-$\kappa$) cell correctly throttles the
  interface flow. If **either** side's mobility is below $10^{-30}$ the face
  mobility is **exactly 0** — an impermeable cell is an infinite series
  resistance, and gas cannot cross into it. The Jacobian operator's
  `_dt_mobility_face` mirrors this rule.

- **Viscosity — per-cell temperature.** Viscosity is an intensive,
  smoothly varying gas property, so it is evaluated per cell at the cell
  temperature, $\mu_i=\mu_g(T_i)$, and enters the face through the per-cell
  mobilities $\kappa_i/\mu_i$ that the distance-weighted harmonic mean
  combines. $\mu_g$ is the material's `gas_viscosity` (§5.9); the default is
  Sutherland's law for air,

  $$
  \mu_g(T) = \mu_{\text{ref}}\left(\frac{T}{T_{\text{ref}}^{\mu}}\right)^{3/2}
             \frac{T_{\text{ref}}^{\mu}+S}{T+S},
  $$

  with the canonical air constants $\mu_{\text{ref}}=1.716\times10^{-5}$ Pa·s,
  $S=110.4$ K. **Caveat on the reference temperature:** the Sutherland reference
  in the code is $T_{\text{ref}}^{\mu}=273.15$ K (the temperature at which
  $\mu_{\text{ref}}$ is quoted for air), which is *distinct* from the enthalpy
  datum $T_{\text{ref}}=298.15$ K used elsewhere in the model. The two should not
  be conflated.

- **Guards.** A cell with $\kappa_i\le0$ or $\mu_i\le10^{-20}$ Pa·s has zero
  mobility, which zeroes the mobility of both adjacent faces (a fully solid or
  numerically degenerate cell carries no Darcy flow).

*Implementation note.* The face mobility and velocity are assembled in
`update_darcy_properties!` (`src/Discretization/finite_volume.jl`); the viscosity is
`effective_gas_viscosity(T, material)`, `src/Properties/effective_properties.jl`.

The effective permeability $\kappa$ entering the cell-center field is either a
component mixing rule (PARALLEL/SERIES/WEIGHTED on the renormalized solid
basis, with the series bound exactly 0 when any condensed phase is
impermeable; §5.7.1) or the porosity-based **Carman–Kozeny** closure

$$
\kappa = \min\!\left(\frac{\ell^2\,\varphi^3}{180\,(1-\varphi)^2},\;\ell^2\right),
$$

with $\ell$ the characteristic pore/particle diameter and $180$ the Kozeny
constant for a packed bed of spheres; as $\varphi\to1$ the packed-bed form
diverges and $\kappa$ is capped at the free-flow scale $\ell^2$ (binding for
$\varphi\gtrsim0.93$), with a $10^{-30}$ denominator guard solely against
$\mathrm{Inf}$ at exactly $\varphi=1$. Carman–Kozeny is
evaluated on the *canonical* porosity $\varphi(\xi,T)$ — the same porosity the
pressure formula uses — not on bulk-density volume fractions (see §5.7). The two
closures (component mixing vs. Carman–Kozeny) are mutually exclusive, and
BRUGGEMAN EMT is rejected for permeability (no physical basis for a transport
resistance; an `ArgumentError` is thrown). As virgin material converts to porous
char, $\varphi$ rises and $\kappa$ increases, easing gas escape.

### 9.5.3 Pressure boundary conditions and the Darcy venting flux

Pressure BCs do not enter the diffusive surface flux; they set the boundary-face
Darcy velocity, which drives the **advective venting flux into the species
equations**. The chain is: `boundary_darcy_velocity` (one shared scalar kernel
per BC type) → `apply_pressure_bc_to_darcy_velocity!` writes $v_b$ at the
boundary faces → `compute_gas_fluxes_darcy_fick!` writes the boundary face flux

$$
J_j = v_b\,\frac{\xi_j^{\text{int}}}{\varphi_{\text{int}}}
\quad\text{(outflow only)},
$$

upwinding against a zero-concentration ambient ghost: only outflow advects the
interior gas, inflow advects nothing (the ambient carries none of the tracked
species). The diffusive film-law part of the boundary flux is then **added** by
`apply_mass_bc_flux!`. Venting therefore depends on permeability and the
pressure gradient — the quantity Darcy–Fick mode exists to model.

With $d$ the cell-center-to-face distance,
$\kappa_{\text{int}},\mu_{\text{int}}$ the adjacent interior cell's
properties, and the tag-aware sign (kernel methods in
`boundary_darcy_velocity`):

- **Dirichlet** ($P=P_{\text{bc}}$):
  $\;v_b = -\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,(\nabla P)$, with
  $(\nabla P)=(P_{\text{bc}}-P_{\text{int}})/d$ at the top and
  $(P_{\text{int}}-P_{\text{bc}})/d$ at the bottom (default
  $P_{\text{bc}}=101325$ Pa).
- **Neumann** ($\partial P/\partial n$ given along the **outward normal**):
  $\;v_b = -\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\partial P/\partial n$
  at the top, and $+\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\partial
  P/\partial n$ at the bottom ($\hat n=-\hat z$ there, so
  $\partial P/\partial z=-\partial P/\partial n$);
  zero gradient $\Rightarrow$ impermeable.
- **Convective (Robin)**: $\;v_b = \pm\,h_P\,(P_{\text{int}}-P_{\text{ext}})$,
  with the $+$ sign at the top boundary and $-$ at the bottom (so positive
  outflow is consistent with the $+z$ convention), and $h_P$ [m·Pa⁻¹·s⁻¹] the
  pressure-transfer coefficient — gated on $\kappa_{\text{int}}>0$ (an
  impermeable surface cell cannot vent regardless of $h_P$). As
  $h_P\to\infty$ this approaches Dirichlet; as $h_P\to0$ it approaches a
  zero-gradient (impermeable) Neumann condition.

These are detailed in §12.5; the shared kernel is `boundary_darcy_velocity`
(`src/Discretization/darcy_velocity.jl`), used identically by the residual, the
`BCMassOperator` Jacobian, and the saved-output MLR diagnostics. Each pressure
BC also re-applies the $\kappa>0,\ \mu>10^{-20}$ guards, and Dirichlet
additionally zeroes the velocity if $d<10^{-15}$ m.

### 9.5.4 Initializing the pore gas: background air

Darcy–Fick problems must **not** start with an evacuated pore space. With all
gas $\xi_j = 0$ the closure gives $P \approx 0$ inside the solid against a
101325 Pa boundary value, so early product gas is spuriously advected *inward*
down the absolute-pressure slope instead of being resisted by ~1 atm of pore
air (peer codes initialize $P = P_\infty$ with a background species or solve
gauge pressure). Two ingredients make the initial state physical:

1. **A real pore space.** The condensed components must declare an intrinsic
   (skeletal) density via `ρ_intrinsic` (or the `φ` porosity kwarg) on
   `SolidComponent`/`LiquidComponent`. Without it, `get_intrinsic_density`
   falls back to the bulk density — a phase with no pores of its own — and a cell
   at full bulk density has $\varphi = 0$: there is no pore volume to hold
   air, and any generated gas hits the $\varphi$ floor ($10^{-12}$) in
   $P = N/\varphi$, producing absurd pressures.
2. **An inert background air component** initialized at ambient pressure:

   ```julia
   air  = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
   ξ_air = background_gas_concentration(material, :AIR, ξ_condensed₀, T₀)
   # → use ξ_air as the :AIR entry of the initial composition
   ```

   `background_gas_concentration` (src/Physics/pressure.jl) is the
   single-species inverse of the pressure closure,
   $\xi_j = P\,M_j\,\varphi(\xi,T)/(R\,T)$, so the initial cell pressure is
   exactly `P` (101325 Pa by default).

The shipped DarcyFick examples (`07_douglas_fir_cone.jl`,
`08_darcy_fick_demo.jl`, `08_pressure_treated_wood_cone.jl`) follow this
pattern. Caveat: the diffusive film BC drives every species toward a
zero-concentration ambient, so the background air slowly leaks out over an
internal-diffusion timescale (~minutes for cone-scale slabs); the
initialization matters most early in heating, before significant gas
generation. Advective air ingress from the ambient is not
modeled (boundary Darcy inflow upwinds against a zero-concentration ghost,
§9.5.3).

---

## 9.6 Pressure derivatives in the Jacobian

The structured-operator Jacobian (§15) requires the sensitivity of
$P=N(\xi,T)/\varphi(\xi,T)$ to the state. In the clamp-inactive interior
($0<\varphi<1$, $\varphi$ above its $10^{-12}$ floor) the quotient rule gives

$$
\frac{\partial P}{\partial T} = \frac{1}{\varphi}\,\frac{\partial N}{\partial T},\qquad
\frac{\partial P}{\partial \xi_k} = \frac{1}{\varphi}\frac{\partial N}{\partial \xi_k}
- \frac{N}{\varphi^{2}}\frac{\partial\varphi}{\partial\xi_k},
$$

with $\partial N/\partial\xi_k = R_gT/M_k$ for a gas and $0$ otherwise, and
$\partial\varphi/\partial\xi_k = -1/\rho_{i,k}$ for a condensed phase and $0$ for
a gas. Two physical readings follow: a **gas** component raises $P$ directly
through the partial-pressure sum; a **solid/liquid** component raises $P$
indirectly — adding condensed mass shrinks the void fraction, confining the same
gas to a smaller volume.

These derivatives are not hand-coded. The `DarcyTransportOperator`
(src/Jacobian/operators/darcy_transport.jl) seeds the stencil state as
`ForwardDiff.Dual`s and evaluates the *same* functions the residual calls —
`porosity` → `ideal_gas_pressure` → `effective_permeability` →
`effective_gas_viscosity` — extracting
$(P,\partial P/\partial T,\partial P/\partial\xi)$, $(\kappa,\dots)$,
$(\varphi,\dots)$ and $(\mu,d\mu/dT)$ from the Dual partials. Differentiating
through the executed code makes the Jacobian consistent with every guard by
construction: at the porosity floor, AD through `max(φ, 1e-12)` returns a zero
$\varphi$-partial, exactly matching the residual's clamped value (a hand-written
$N/\varphi^2$ would blow up as $10^{24}N$ there). The chain also picks up
temperature dependence the interior formulas above elide — intrinsic densities
$\rho_{i,j}(T)$ inside $\varphi$, and user-supplied $\kappa_j(T)$ — with no
separate code path to keep in sync.

---

## 9.7 The Darcy property cache and the velocity Jacobian

Darcy quantities are recomputed **once per residual evaluation** and stored in
the conditional `DarcyCache` (per-cell $P,\varphi,\kappa,\mu$; per-face
$v_g,\nabla P$), which is allocated only when Darcy flow is enabled (§5.13).
`update_darcy_properties!` (src/Discretization/finite_volume.jl) fills the
per-cell fields from the state ($\varphi$ → $P$ → $\kappa$ → $\mu$), forms each
interior face's mobility as the distance-weighted harmonic mean of the per-cell
$\kappa_i/\mu_i$ (exactly 0 if either side is $\le10^{-30}$; §9.5.2), writes
$v_g[f]=-\text{mob}_f\,(\nabla P)_f$, and finishes by applying the pressure BCs
to the boundary faces (§9.5.3).

Holding the mobility fixed, the velocity's pressure sensitivity is exact and
constant in $P$:

$$
\frac{\partial v_g[f]}{\partial P_L} = +\frac{\text{mob}_f}{d_{LR}},\qquad
\frac{\partial v_g[f]}{\partial P_R} = -\frac{\text{mob}_f}{d_{LR}}.
$$

The structured Jacobian does not consume a cached coefficient table for this:
the `DarcyTransportOperator` (§15.8.3) reconstructs $v_g$ from the stencil's
cell-local $(P,\kappa,\mu)$ inside its own closure (`_dt_mobility_face`,
src/Jacobian/operators/darcy_transport.jl — the same face rule as the residual)
and obtains both the direct $\partial v_g/\partial P$ and the full chain to
$(T,\xi)$ — through $P$, $\varphi$, $\kappa(T,\xi)$, and $\mu(T)$ — from the
forward-mode pass of §9.6, with an analytical override gated against the AD path
at `rtol = 1e-12`.

---

## 9.8 Surface diffusive mass flux

At the exposed surface gas escapes to the environment. The surface diffusive flux
is obtained from the internal-diffusion-only form $J=-\lambda\,\partial\xi/\partial z$
with the environment assumed to hold *zero* concentration of decomposition
products, so the boundary gradient is approximated as $(\xi-0)/d=\xi/d$ with $d$
the distance from the adjacent cell center to the boundary face.

### 9.8.1 Infinite vs. finite external resistance

- **Pure internal diffusion** (no external film resistance, $h_m\to\infty$):

  $$
  J_j^{\text{surface}} = \frac{\lambda\,\xi_j}{d}.
  $$

  The implementation carries the surface coefficient as $\lambda\rho$ and writes
  $J=\lambda\rho\,(\xi/\rho)/d=\lambda\xi/d$; the two forms are identical because
  $\rho$ cancels (consistent with the volume-fraction convention of §9.1.4). The
  effective $\lambda$ is recovered by exact division $\lambda=\lambda\rho/\rho$
  (the kernel returns zero flux for $\rho \le 0$; no floor is applied — the
  callers pass $\lambda\rho = \lambda\cdot\rho$, so the division is exact for any
  positive density).

- **Combined internal + external resistance** (finite film coefficient
  $h_m$ [m·s⁻¹]):

  $$
  \boxed{\;
  J_j^{\text{surface}} = \frac{\lambda\,h_m\,\xi_j}{\lambda + h_m\,d}
  \;}
  $$

  This is two resistances in series. Writing the flux as
  $J=\xi_j/(R_{\text{int}}+R_{\text{ext}})$ with internal resistance
  $R_{\text{int}}=d/\lambda$ and external resistance $R_{\text{ext}}=1/h_m$,

  $$
  J_j = \frac{\xi_j}{\dfrac{d}{\lambda}+\dfrac{1}{h_m}}
      = \frac{\lambda h_m \xi_j}{h_m d + \lambda},
  $$

  recovering the boxed form. As $h_m\to\infty$, $R_{\text{ext}}\to0$ and
  $J\to\lambda\xi_j/d$ (pure internal limit); as $h_m\to0$, $J\to h_m\xi_j$
  (external-controlled).

A purely **convective** mass BC (external transfer only) gives
$J_j^{\text{surface}}=h_m(\xi_j-\xi_{j,\infty})$, with the ambient concentration
$\xi_{j,\infty}$ typically zero for decomposition products venting to atmosphere.
Pressure and impermeable BCs return zero surface *diffusive* flux (a pressure BC
acts through the Darcy velocity, §9.5.3).

### 9.8.2 Two sign conventions, reconciled

The diagnostic `compute_surface_mass_flux` (`Discretization/boundary_fluxes.jl`)
returns **positive flux for gas escaping** — the natural sense for a surface efflux.
This is the *opposite* of the general BC convention used by `mass_bc_flux`
(`Discretization/boundary_fluxes.jl`), which returns **negative for gas escaping**
because its convention is "positive = into the domain."

The two are reconciled when the flux is written into the face array. Because the
finite-volume divergence is $(J_R-J_L)A/V$ and the top face is the *right* face of
the surface cell, the BC flux is negated at the top boundary and applied directly
at the bottom:

$$
J[f_{\text{top}}] = -J_{\text{bc}}, \qquad J[f_{\text{bottom}}] = +J_{\text{bc}}.
$$

The sign flip at the top is **not** a bug; it is the consistent realization of the
signed divergence, so that gas escaping at $z=L$ produces a positive flux in the
$+z$ direction and a positive (outflow) contribution to the surface cell's species
balance. The full sign-convention summary is:

| Quantity | Convention |
|---|---|
| Interior face flux $J_j[f]$ | $+$ = transport in $+z$ |
| `compute_surface_mass_flux` | $+$ = gas escaping |
| `mass_bc_flux` (BC) | $-$ = gas escaping |
| Stored top-face flux | $J[f]=-J_{\text{bc}}$ |
| Stored bottom-face flux | $J[f]=+J_{\text{bc}}$ |

*Implementation note.* `mass_bc_flux(bc, ξ_surface, λρ, ρ, d, t, j)` (with the
sign-convention block on the `DiffusiveMassBC` method) and
`compute_surface_mass_flux`, `src/Discretization/boundary_fluxes.jl`; face application in
`apply_mass_bc_flux!`, `src/Discretization/finite_volume.jl`. See §12.4 for the
full mass-BC catalogue (Impermeable, MassFlux, Convective, Diffusive, Combined).

---

## 9.9 Species conservation and the gas-only divergence

The species mass-conservation PDE for a gaseous component $j$ is

$$
\frac{\partial \xi_j}{\partial t} = -\nabla\!\cdot J_j + S_{j}^{\text{rxn}},
$$

while a solid or liquid component has no transport flux and evolves by reaction
(and, in ALE, mesh dilation) only:

$$
\frac{\partial \xi_j}{\partial t} = S_{j}^{\text{rxn}}
\qquad (j \text{ solid/liquid}),
$$

(see §3.2 for the gas-vs-solid dispatch and §11.2 for the ALE corrections). The
reaction source $S_j^{\text{rxn}}=\sum_r \nu_{j,r}\,r_r$ is treated in §6.

### 9.9.1 Discrete finite-volume divergence

Integrating the flux divergence over cell $i$ and applying the divergence theorem
in 1D yields the face-difference operator

$$
\boxed{\;
\left.\frac{d\xi_j}{dt}\right|_{\text{transport},\,i}
   = -\frac{A}{V_i}\Big(J_j^{R} - J_j^{L}\Big)
\;}
$$

where $A$ is the cross-sectional area and $V_i$ the cell volume. The discrete
operator uses **one** area for both faces (`divergence_species_1d!` reads
`state.A_face[f_left]`): all faces share the global scalar $A_{\text{section}}$
(also under `WithChi`, where $A(t)$ is column-global), so this is exactly
conservative — the operator does not support per-face areas. The
multi-component form is assembled by `divergence_species_1d!` over all components,
but only gas components feed the species residual — the dispatch
`apply_species_residual_unrolled!` writes

$$
\frac{d\xi_{j,i}}{dt} = -\,\text{div}_J[j,i] + S_{j,i}^{\text{rxn}}
   \quad(j\text{ gas}),
\qquad
\frac{d\xi_{j,i}}{dt} = S_{j,i}^{\text{rxn}}
   \quad(j\text{ solid/liquid}).
$$

This gas-vs-solid branch is a compile-time-eliminated `is_gas` test (zero runtime
cost). Solid components have $J_j\equiv0$ at all faces, so even if their divergence
were formed it would vanish; the explicit dispatch makes the intent — and the
mass-thickness locking of solids to the mesh (§11) — unambiguous.

The gas flux that feeds this divergence is the diffusion-only flux (Eulerian
Fickian mode) or the Darcy–Fick flux (when permeability is present), selected by
the simulation-mode trait (§15.1). In ALE mode the same skeleton-relative flux feeds
the divergence unchanged: every component rides the mesh, so the moving frame adds
only the dilation of §11.3.4 and no advective correction.

### 9.9.2 Coupling to the energy equation

Gas transport feeds back on temperature through one volumetric source of the energy
equation, the advective source (the advection piece of the gas-enthalpy divergence,
derived in detail in §7.4):

$$
S_{\text{conv},i} = -\sum_{j:\,\text{gas}} c_{p,j}(T_i)\,\bar{J}_{j,i}\,
    \left(\frac{\partial T}{\partial z}\right)_i,
\qquad \bar{J}_{j,i}=\tfrac12\big(J_j^{L}+J_j^{R}\big),
$$

with $(\partial T/\partial z)_i$ the opposite-distance-weighted gradient in the
interior; boundary cells blend the interior one-sided slope with the
boundary-face slope at the solved $T_s$ (§7.4.2).

It appears in the assembled temperature RHS
$\rho c_p^{\text{eff}}\,dT_i/dt = -(\nabla\!\cdot q)_i + S_{\text{conv},i}
+ Q_{\text{rad},i} + Q_{\text{rxn},i}$, whose storage coefficient sums every
component including the gas (§5.4, §7.3); the gas *flux* enters only through
$S_{\text{conv}}$. The surface energy balance therefore carries no gas-enthalpy
term: one on top of $S_{\text{conv}}$ would count the advected enthalpy a second
time (§7.4.3, §12.2).

*Implementation note.* `divergence_species_1d!`,
`src/Discretization/finite_volume.jl`; `apply_species_residual_unrolled!`,
`src/Materials/materials.jl`; `compute_gas_convective_source!`,
`src/Physics/convection.jl`. The 8-step `compute_rhs!` ordering (transport
properties and interior conductive fluxes $\to$ reactions $\to$ gas transport
$\to$ thermal BCs $\to$ radiation $\to$ $S_{\text{conv}}$ $\to$ divergences
$\to$ final RHS) is described in §7.6.

---

## 9.10 Direct gas-flux Jacobian derivatives

The `GasTransportOperator` contributes the analytic Jacobian of the gas-transport
terms over a 3-cell tridiagonal stencil ($i-1,i,i+1$). The core building block is
the set of direct partial derivatives of the single-face diffusive flux

$$
J = -\lambda\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d},
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

with respect to each scalar argument. Writing $N_{\text{num}}=\xi_R T_R-\xi_L T_L$
and $\text{denom}=T_{\text{face}}\,d$:

$$
\frac{\partial J}{\partial \lambda} = -\frac{N_{\text{num}}}{\text{denom}},
\qquad
\frac{\partial J}{\partial \xi_L} = \frac{\lambda T_L}{\text{denom}},
\qquad
\frac{\partial J}{\partial \xi_R} = -\frac{\lambda T_R}{\text{denom}},
$$

$$
\frac{\partial J}{\partial T_L} = \frac{\lambda \xi_L}{\text{denom}} - \frac{J}{2\,T_{\text{face}}},
\qquad
\frac{\partial J}{\partial T_R} = -\frac{\lambda \xi_R}{\text{denom}} - \frac{J}{2\,T_{\text{face}}}.
$$

The temperature derivatives have two pieces: a direct effect through the numerator
($\pm\lambda\xi/\text{denom}$) and an indirect effect through the
$T_{\text{face}}$ prefactor in the denominator
($-J/(2T_{\text{face}})$, since $\partial T_{\text{face}}/\partial T_{L,R}=\tfrac12$
and $\partial(1/T_{\text{face}})/\partial T_{\text{face}}=-1/T_{\text{face}}^2$).
A guard sets $\partial J/\partial\lambda=0$ when $N_{\text{num}}=0$ exactly (a flat
$\xi T$ field has no diffusive flux and no $\lambda$-sensitivity).

Beyond these direct partials, the operator composes two additional chains:

- **Mixing-rule chain (via $\lambda_{\text{face}}$).** Since
  $\lambda_{\text{face}}$ is a harmonic mean of the cell coefficients, and each
  $\lambda_{L/R}$ is itself a composition-dependent mixing rule, the flux
  sensitivity to a concentration also includes

  $$
  \frac{\partial J}{\partial \xi_{L/R,k}}
  \;\supset\;
  \frac{\partial J}{\partial \lambda}\cdot
  \frac{\partial \lambda_{\text{face}}}{\partial \lambda_{L/R}}\cdot
  \frac{\partial \lambda_{L/R}}{\partial \xi_{L/R,k}},
  $$

  with the last two factors supplied by the harmonic-mean derivative helper and
  the effective-gas-transfer mixing derivatives (§5.6, §5.15).

- **Convective-source chain (via $S_{\text{conv}}$).** The energy-equation
  coupling contributes terms such as

  $$
  \frac{\partial S_{\text{conv}}}{\partial T_C} \;\supset\;
  \left(\frac{\partial c_p}{\partial T}\Big|_{T_C}\right)\bar{J}\,
  \left(\frac{\partial T}{\partial z}\right),
  $$

  i.e. the explicit $c_p$-chain through the temperature dependence of the gas heat
  capacity, plus the dependence of $\bar J$ and of the gradient on the stencil
  temperatures. The species rows likewise receive $-\partial(\text{div}_J)/\partial(\cdot)$
  via the divergence $\text{div}_J=(J_R A_R-J_L A_L)/V_C$, written into the species
  block as $\text{out}[1+j]=-\text{div}_J$ for gas components and $0$ for
  solids/liquids.

The operator's `contribute_to_jacobian!` is the analytical override
(`_contribute_to_jacobian_analytical!`); a ForwardDiff reference over the full
`apply_cell!` closure (`_contribute_to_jacobian_ad!`) is kept for validation and is
gated against it at `rtol = 1e-12` in
`test/jacobian/test_operator_gas_transport_analytical.jl`. Boundary-face
contributions are zeroed in the operator — the mass BCs and pressure BCs are
handled by `BCMassOperator`, and the boundary cell's $S_{\rm conv}$ energy row
by `BCConvectionOperator` (§15.8.8) — keeping
the stencil pure and avoiding double-counting.

*Implementation note.* `_gt_J_direct_partials`, `_gt_species_div_unrolled!`, and
the analytic operator `_contribute_to_jacobian_analytical!`, all in
`src/Jacobian/operators/gas_transport.jl`. The Darcy-velocity pressure
derivatives that close the advective term's Jacobian are in §9.7. See §15.7–15.9
for how these operator contributions assemble into the structured Jacobian.

---

## 9.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The Pyrolysis.jl gas-transport model is a faithful descendant of the ThermaKin2Ds
and Gpyro condensed-phase formulations, with the FDS solid-phase model as the
minimal-transport reference. The notation map below is used on first contact;
internally all chapters use the unified symbols of §2.

| External | Pyrolysis.jl | Meaning |
|---|---|---|
| ThermaKin $\theta_i^{\,j}$ | $\nu_{i,j}$ | stoichiometric yield |
| Gpyro $K$ | $\kappa$ | permeability |
| Gpyro $\bar\rho$, $X_\alpha$ | $\rho$, $v_j$ | mixture density, volume fraction |
| Gpyro $\nu$ (Darcy) | $\mu/\rho_g$ | kinematic viscosity |
| FDS $Z$ | $A_i$ | Arrhenius pre-exponential |

### ThermaKin2Ds

ThermaKin writes the gas mass flux as a volume-fraction-gradient diffusion law,
$J = -\rho_g\,\lambda\,\partial(\xi_g/\rho_g)/\partial x$, which — for an ideal gas
with $\rho_g\propto P/T$ — is identical to the Pyrolysis.jl form
$J=-(\lambda/T)\,\partial(\xi T)/\partial x$ (§9.1; the header of
`src/Physics/gas_transport.jl` cites ThermaKin §1.2 explicitly). The two codes
therefore share the same diffusive driving potential and the same
composition + thermal split. Pyrolysis.jl additionally offers the Darcy advective
augmentation (§9.3) for permeable/pressure-driven problems; in the rigid,
non-swelling limit ($\gamma=0$ for gases) ThermaKin's Boyle's-law Darcy treatment
maps onto the porous-medium ideal-gas pressure (§9.4) with $\kappa$ controlling the
flow resistance. The energy-equation convention is the same as ThermaKin's:
Eq. 17 sums $\sum_j \xi_j c_j$ over *all* components including gas and carries
the gas flux through $c_g(\mathbf J_g\cdot\nabla T)$, exactly $\rho c_p^{\text{eff}}$
(§5.4) and $S_{\text{conv}}$ (§9.9.2). ThermaKin's Eq. 15 applies its deformation
term to every component, gas included — the skeleton-bound frame of §11.1.2.

### Gpyro

Gpyro solves pressure-driven gas transport with a Darcian mass flux
$\dot m'' = -(K/\nu)(\partial P/\partial z - g\rho_g)$ together with a
pressure-evolution equation derived from gas continuity and the ideal-gas law.
Pyrolysis.jl uses the same Darcy constitutive law,
$v_g=-(\kappa/\mu)\nabla P$ (§9.5.2), with the distance-weighted harmonic mean
of the per-cell mobilities $\kappa/\mu$ (the material's viscosity curve at the
cell temperature); the gravity/buoyancy term $g\rho_g$ is omitted (negligible
for thin condensed-phase samples in the through-thickness direction). The pressure
is closed algebraically by the porous-medium ideal-gas law $P=N/\varphi$ (§9.4)
rather than time-integrated as a separate state. This is **not** an additional
approximation relative to Gpyro's pressure ODE: gas storage $\partial\xi_g/\partial t$
is fully retained in the species ODEs, and $P=N/\varphi$ is an exact EOS
evaluation of that evolving state — summing the species equations reproduces gas
continuity with no spurious source. Gpyro's control-volume conservation of
$\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, and gas mass — the gas balance written
on control volumes that deform with the solid, with only $\dot m''$ crossing their
faces — corresponds to the cell-integrated species divergence of §9.9.1 on the
moving (ALE) mesh of §11, where every component rides the mesh; the harmonic
interface interpolation of Gpyro's Patankar scheme is the same resistance-in-series
face averaging used here (§9.1.4, §13.11).

### FDS

The FDS solid-phase model treats the condensed phase with 1D Fourier conduction
and Arrhenius kinetics but, in its standard form, does **not** resolve in-depth
gaseous mass transport — volatiles are assumed to escape instantaneously to the
gas phase. The Pyrolysis.jl `NO_RADIATION`/Fickian configuration with negligible
permeability reduces toward this minimal-transport limit; the explicit
diffusive/Darcy machinery of this chapter is the generalization FDS lacks, and
with no in-depth gas FDS has neither gas storage nor an in-depth $S_{\text{conv}}$.

---

## 9.12 Limitations and open issues

- **First-order upwind advection only.** The Eulerian Darcy–Fick advective term
  uses unlimited first-order upwinding (§9.3.1). It is stable and conservative but
  diffusive; sharp composition fronts are smeared. Higher-order TVD limiting is
  implemented only for the mesh-relative temperature advection of the
  verification hook (§11.6) and is not exposed for the Eulerian gas flux.

- **No species discrimination in diffusion.** All gas species share one
  effective transfer coefficient $\lambda_{\text{eff}}$ (§9.2): light and heavy
  vapors diffuse identically within a simulation. This is a simplified
  surrogate for a full multicomponent (Stefan–Maxwell) diffusion model; the
  Chapman–Enskog Lennard-Jones parameter estimates available for sizing
  $\lambda$ (§5.10) are heuristic and unvalidated for pyrolysis vapors.

- **No Wilke mixing for viscosity.** The Darcy viscosity is one curve per material
  (Sutherland's law for air by default; §5.9, §9.5.2), independent of the pore-gas
  composition; a Wilke mixing rule would be required to capture mixtures of very
  different molecular weights (e.g. H₂O + CO₂ at high temperature; §5.16). The
  Sutherland reference temperature (273.15 K) differs from the enthalpy datum
  (298.15 K) and must not be conflated.

- **Zero-exterior-concentration surface assumption.** Surface diffusive fluxes
  assume the environment holds no decomposition products ($\xi_\infty=0$ implicit;
  §9.8.1). Re-entry or accumulation of products in the boundary layer is not
  modeled.

- **One-temperature gas energy.** The pore gas is stored at the matrix
  temperature and its flux enthalpy is carried by $S_{\text{conv}}$; there is no
  two-temperature gas energy equation, so regimes where the gas does not
  equilibrate with the pore walls within a cell are outside the model. The
  pressure itself is *not* approximated: $P=N/\varphi$ is an exact
  equation-of-state evaluation with gas mass storage fully retained in the
  species ODEs (§3.4.2).

- **Negative-advective-flux robustness.** When upwinding selects a zero-concentration
  donor cell, the advective term can transiently produce a "wrong-side" flux
  (§9.3.1); physicality over long runs relies on reaction/heat coupling rather
  than an explicit positivity-preserving limiter.


---

# 10. Volume Change, Swelling, Intumescence, and Lateral Shrinkage

This chapter derives the complete volume-change model of Pyrolysis.jl: how reaction-driven mass exchange between condensed and gaseous phases produces a per-cell volumetric strain rate $\theta$, how that strain is converted into axial mesh motion (the ALE mesh velocity), and how a column-global *lateral* shrinkage / swelling / intumescence law deforms the cross-sectional area $A(t)$ through the pyrolysis-progress state variable $\chi$. We give the continuous balance, the exact discrete update rules as implemented, the dilation correction that keeps every component's cell mass on its reaction-and-flux budget as the mesh moves, and a careful comparison to the deformation treatments in ThermaKin2Ds, Gpyro, and FDS. Every equation here is grounded in the source: `src/Physics/volume_change.jl`, `src/Residual/residual.jl`, `src/Problem/state_cache.jl`, `src/Problem/initial_conditions.jl`, and `src/Materials/materials.jl`.

The axial mesh-motion machinery (cumulative integration of $\theta_L$, the geometric volume from the node positions, and the verification-only mesh-relative temperature advection) is developed in §11, *The ALE Framework*; here we develop only the *constitutive* content — the strain rates and the cross-section law that feed §11 — and the species dilation correction that closes the species balance under mesh motion. The state-vector layout that carries the $\chi$ block is detailed in §13, *Finite-Volume Discretization*, and §15, *Temporal Integration & Structured Jacobian*.

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

> **Implementation note.** $\theta$ is computed by `volume_change_rate(dξ, T, material)` / the dispatch helpers `_volume_change_rate`, in `src/Physics/volume_change.jl` (inline overloads for $N_C\in\{1,2,3\}$ and a `@generated` form for $N_C>3$). The exact loop is
> ```julia
> θ += γⱼ * dξ[j] / ρⱼ        # ρⱼ = comps[j].density(T)
> ```
> evaluated **only when `γⱼ > 0`** (an early-exit guard: components with $\gamma_j=0$ are skipped, so no density evaluation is wasted and the term is exactly zero). There is **no minus sign** in the code: the sign of $\theta$ is carried entirely by the sign of $d\xi_j/dt$, which is negative for a depleting reactant (the species source term packs reactants with coefficient $-1$; see §6.5). The density $\rho_j(T)$ is the *bulk* density property evaluated at the cell temperature — not the intrinsic/skeletal density $\rho_{i,j}$ (which is reserved for the porosity formula, see §5.3 and §2 overload A7).

### 10.1.2 What $d\xi_j/dt$ is in (10.2)

A subtle and load-bearing point: the $d\xi_j/dt$ fed into (10.2) is the **reaction source term only**, not the full species RHS (which also includes gas transport and the ALE dilation). This is correct because $\theta$ measures the rate at which mass is *converted between phases by chemistry* — the only process that changes the additive-volume balance (10.1). Transport of an already-gaseous species through the pores does not change the matrix volume, and the ALE dilation term is bookkeeping for the moving frame, not new mass exchange. Including them would double-count the volume change that ALE motion is itself representing.

> **Implementation note.** In the ALE residual (`_residual!(..., ::ALE)`, `src/Residual/residual.jl`) the reaction source matrix `ws.rhs_cache.dξ_rxn` is populated by the physics RHS at step 4 and is passed as the positional argument `dξ_rxn` of `_compute_ale_mesh_velocity!`. Inside that function, `dξ_dt[j] = dξ_rxn[j, cell_idx]`, and then `θ_cell = Physics.volume_change_rate(dξ_dt, T_cell, material)`. So $\theta$ is built from `dξ_rxn`, the pure reaction source — exactly the per-cell strain rate of (10.2).

### 10.1.3 Absolute cell volume rate

The absolute rate of change of a cell's volume follows immediately:

$$
\frac{dV_i}{dt} \;=\; V_i\,\theta_i.
\tag{10.3}
$$

> **Implementation note.** No routine materializes $dV_i/dt$: the ALE path evaluates $\theta_i$ with `volume_change_rate` (`src/Physics/volume_change.jl`) inside the mesh-velocity integration (§10.6) and recomputes $V_i$ from the node positions (§11.4.1).

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
- **Charring with a rigid char:** give the char $\gamma=0$ to "freeze" its volume contribution (`rigid_char = SolidComponent(:char, ...; γ=0.0)`, the documented example at `components.jl`); then only the mass converted out of the cell drives axial motion, and the char does not re-inflate or collapse with the mesh. (This is a modeling choice, not a default; the default char is $\gamma=1$.)
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

> **Implementation note.** Computed in `_pack_chi_derivative!` (`src/Residual/residual.jl`). The exact write is
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

> **Implementation note.** `populate_m_dry_solid_initial!(mesh, material)` (`src/Problem/initial_conditions.jl`) sums `mesh.cell_states[i].ξ[j]` over components for which `is_solid(material.components[j])` is true, multiplied by `mesh.cell_volumes[i]`. Crucially the sum is over **solids only** — liquids (e.g. bound moisture) are *not* counted in the denominator, even though they are matrix-forming. This is deliberate: $\chi$ is a *dry-decomposition* progress variable, and water mass is excluded from both numerator (dry-gas only) and denominator (solids only). The array is stored as `Float64` on the mesh and changes only at a park event, which adds the parked cell's entry to the keeper's (§14.3.4), so the column total $M_{\mathrm{dry}}$ is fixed for the whole solve.

### 10.3.3 Dry-gas identification and the water-vapor caveat

The index set $\{j\in\mathrm{dry\ gas}\}$ collects all GAS-phase components *except* water vapor / moisture-derived gases. Phase is a compile-time (type-level) property, but the water-vapor exclusion is done by **name matching**:

> **Implementation note.** `dry_pyrolysis_gas_indices(material)` (`src/Materials/materials.jl`) recurses over the component tuple at compile time (`_collect_dry_pyrolysis_gas_indices`) and includes index $j$ iff `is_gas(comp) && !_is_water_vapor_name(comp.name)`. The name filter `_is_water_vapor_name` is a **case-insensitive substring match**: a name is treated as water vapor if its string contains `"H2O"`, or (lowercased) contains `"water"`, or equals `"h2o"`. The resulting tuple is cached once at workspace construction (`ws.dry_gas_indices`) to avoid a per-residual `Symbol→String` allocation. **Caveat:** because identification is by substring, a dry-gas species accidentally named (say) `:water_gas` would be misclassified as moisture and excluded from $\chi$. Users modeling moisture must name the water-vapor component so it matches (e.g. `:water_vapor` or `:H2O`), and must avoid those substrings in genuine dry-gas names.

### 10.3.4 Why $\chi$ is reconstructed, not freely integrated

Although $\chi$ occupies a trailing block of the state vector $\mathbf{u}$ (the "$\chi$ block"; see §10.9 and §15.1) and is advanced by the ODE integrator, its derivative $d\chi_i/dt$ is *recomputed from the reaction sources on every residual evaluation* by `_pack_chi_derivative!` rather than being a free dynamical degree of freedom with independent physics. Because the RHS for $\chi$ is a pure function of the (reaction part of the) other states and geometry, $\chi$ stays exactly synchronized with cumulative dry-gas production regardless of the integrator's internal representation. For identity-law materials (no lateral shrinkage), the $\chi$ block is absent from the state vector entirely (the `chi_range` is empty), and `_pack_chi_derivative!` returns immediately — so $\chi$ carries zero memory/compute cost when it is not needed.

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

> **Implementation note.** `_chi_bar_from_state(u, layout, mesh)` (`src/Problem/state_cache.jl`) reads the $\chi$ block from `u` and the weights from `mesh.m_dry_solid_initial`, returning the tuple $(\bar\chi, M_{\mathrm{dry}})$. The numerator `num += u[χ_offset+i] * m_i` accumulates in the scalar type of `u` (which is `ForwardDiff.Dual` in a Jacobian/sensitivity pass), while the denominator `M` accumulates in `Float64`. The division `num / M` therefore returns a Dual when `u` is Dual — the average is **AD-safe**, carrying sensitivities of $\bar\chi$ with respect to every state through which $\chi_i$ depends. Guards: an empty `m_dry_solid_initial` or empty `chi_range` returns $(0,0)$, and $M_{\mathrm{dry}}\le\varepsilon$ returns $(0,0)$. The public wrapper is `column_chi_bar_from_state(u, layout, mesh)` (`src/Residual/residual.jl`).

---

## 10.5 The lateral-shrinkage law $A(t)=A_0\,\mathrm{law}(\bar\chi)$

### 10.5.1 Definition

The column cross-sectional area evolves as a prescribed function of the global progress:

$$
\boxed{\;A(t) \;=\; A_0\,\mathrm{law}(\bar\chi)\;},
\tag{10.8}
$$

where $A_0$ is the calibration-reference area fixed at setup and $\mathrm{law}:\,[0,1]\to(0,\infty)$ is a user-supplied smooth function (or `nothing`, meaning identity). The area $A(t)$ is uniform along $z$ (the column has one cross-section at any instant) but time-varying.

> **Implementation note.** Stored as `material.lateral_shrinkage_law::LSL` (`src/Materials/materials.jl`), a concrete `Function` subtype or `Nothing`. The constructor validates that it is a `Function` or `nothing` and the public `Material(...; lateral_shrinkage_law = nothing)` keyword defaults to identity. $A_0$ is `mesh.cross_section_area_initial` (`src/Geometry/mesh.jl`), set once at mesh construction. The update is performed by `update_cross_section_area!(cache, material, mesh, u, layout, ::WithChi)` (`src/Problem/state_cache.jl`):
> ```julia
> χ̄, M = _chi_bar_from_state(u, layout, mesh)
> A = mesh.cross_section_area_initial * material.lateral_shrinkage_law(χ̄)
> cache.A_section[] = A
> ```
> The result is stored in `cache.A_section[]` carrying `eltype(u)` (Dual-safe), and the function returns $(A,\bar\chi,M_{\mathrm{dry}})$ for downstream use. For identity-law materials it short-circuits and returns the unchanged $A_0$; for `NoChi` mode it is a no-op returning the current area. Note the code uses the field name `cross_section_area_initial` for $A_0$; the mutable `cross_section_area` holds the current $A(t)$.

### 10.5.2 Canonical laws

| Law | $\mathrm{law}(\bar\chi)$ | Behavior |
|---|---|---|
| Identity (default) | $1$ | Constant cross-section $A=A_0$ (standard 1D pyrolysis) |
| Linear shrink | $1-k\bar\chi$, $0<k\le 1$ | Area contracts linearly with extent of decomposition |
| Linear swell | $1+k\bar\chi$, $k>0$ | Area expands (intumescence) |
| Smooth (sigmoid) | $1-\tfrac{k}{1+e^{-s(\bar\chi-\bar\chi_0)}}$ | Onset-localized transition |

For example, the verification fixture uses $\mathrm{law}(\bar\chi)=1-0.2\,\bar\chi$, so at $\bar\chi=0.5$ the area is $A=0.9\,A_0$ (`test/geometry/test_variable_cross_section.jl`, which asserts `A ≈ A0 * 0.9`).

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

> **Implementation note.** `_compute_theta_A(material, state, χ̄, M, dξ_rxn, layout, idx_dry)` (`src/Residual/residual.jl`) computes
> ```julia
> dχ̄_dt = Σ_i ( Σ_{j ∈ idx_dry} dξ_rxn[j,i] ) * state.volume[i]   # = Σ S_dry,i V_i
> dχ̄_dt /= M                                                       # M = M_dry
> A_ratio = law(χ̄)
> dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)                        # law'(χ̄) by AD
> return dlaw_dχ̄ * dχ̄_dt / A_ratio
> ```
> matching (10.10)–(10.11) term for term. The derivative $\mathrm{law}'(\bar\chi)$ is obtained at runtime by `ForwardDiff.derivative`, so the user supplies only $\mathrm{law}$, never its derivative, and any smooth law works. **Guards** return $0$ when: the law is identity (`law === nothing`); the reaction source is unavailable (`dξ_rxn === nothing`); $M_{\mathrm{dry}}\le\varepsilon$; the cell count is zero; or the area ratio $\mathrm{law}(\bar\chi)\le\varepsilon$ — the last prevents the $1/\mathrm{law}(\bar\chi)$ in (10.10) from blowing up as the cross-section vanishes. The mode dispatch `_compute_theta_A_dispatched(...)` returns `0.0` for `NoChi` and calls `_compute_theta_A` for `WithChi`.

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

Equivalently $w_{i+1}=\sum_{k=1}^{i}\theta_{L,k}\Delta z_k$, with the bottom node pinned at $w_1=0$ (the substrate is immobile; material recedes downward and the surface moves). The full development of (10.13) — the node positions as state ($dz_i/dt = w_i$), the geometric volume $V_i = A|z_{i+1}-z_i|$ recomputed from them, and the role of the area-change term — is the subject of §11, *The ALE Framework*. Here we only establish the *constitutive split* (10.12) and why $\theta_A$ must be removed.

> **Implementation note.** `_compute_ale_mesh_velocity!(w, state, mesh, material, dξ_rxn; θ_out, dξ_dt, θ_A)` (`src/Residual/residual.jl`) takes the per-cell reaction species sources `dξ_rxn` (`NC × n_cells`) as a positional argument, initializes `w[1] = 0`, then for each live cell computes `Δz_cell = state.volume[cell_idx] / state.A_section[]` (i.e. $\Delta z_i=V_i/A$), evaluates `θ_cell` via `volume_change_rate`, forms `θ_L = θ_cell - θ_A`, accumulates `cumulative_velocity += θ_L * Δz_cell`, and assigns the result to the cell's right node; parked ghost nodes take the surface velocity. The radial channel `θ_A` is passed in as the keyword `θ_A` (default `0.0`), so for `NoChi` it is identically zero and (10.13) collapses to the constant-cross-section integral $w_{i+1}=w_i+\theta_i\Delta z_i$. The per-cell $\theta_i$ is also exported through the optional `θ_out` array, which is later consumed by the dilation correction (§10.7). The residual and the `ALEMeshVelocityOperator` of the Jacobian call this one function (and `_compute_theta_A`); the operator first re-evaluates the per-cell reaction sources into its own `dξ_rxn` buffer from the same state.

### 10.6.3 Sequencing of the area update and the strain integration

The two channels are computed in a specific order each residual evaluation, because $\theta_A$ depends on $\bar\chi$ (through $\mathrm{law}$ and $\mathrm{law}'$) and on the reaction sources, while the axial integration depends on $\theta_A$ and on the *current* geometry (volumes $V_i$ and area $A$). The ALE residual therefore does, in sequence (`src/Residual/residual.jl`):

1. `populate!` the state cache from $\mathbf{u}$ (T, ξ, z, χ).
2. `update_cross_section_area!` → returns $(A,\bar\chi,M_{\mathrm{dry}})$ and stores $A$ in `cache.A_section[]`.
3. `update_state_geometry!(..., ALE())` → recompute $V_i$, face areas, cell centers from $z$ and the **new** $A$ (so the $\Delta z_i$ in step 6 use the updated cross-section).
4. Properties update; then physics RHS (populates `dξ_rxn`).
5. `_compute_theta_A_dispatched(...)` → $\theta_A$ from (10.10)–(10.11).
6. `_compute_ale_mesh_velocity!(...; θ_A)` → integrate (10.13).

This ordering guarantees the radial and axial channels are mutually consistent at the single instant of the residual evaluation: the same $\bar\chi$ sets both $A$ (hence $\Delta z_i$) and $\theta_A$.

---

## 10.7 The ALE dilation correction

### 10.7.1 The problem

In the moving (ALE) frame, the mesh velocity $w$ already represents the volume change of each cell: when a cell's volume shrinks because its solid is being consumed, the nodes move so that $\Delta z_i$ decreases. But the species concentration $\xi_j$ (mass per unit *bulk* volume) is advanced by a right-hand side written for a fixed volume. Every component rides the mesh — the condensed components because they *are* the skeleton, the gas because the transport flux $\mathbf J_g$ the closures return is defined relative to the skeleton (§3.2.3, §11.1.2) — so each component's cell mass must change only by reaction and by $\mathbf J_j$ through the faces. The concentration, however, must also respond to the changing volume; the **dilation correction** supplies that response.

### 10.7.2 Derivation

For any component the mass in a cell, $m_j=\xi_j V$, changes by reaction and by the skeleton-relative flux through the faces: $dm_j/dt = S_j^{\mathrm{rxn}}V - [J_jA]_{\text{faces}}$, with $J_j \equiv 0$ for condensed components. Expanding,

$$
\frac{d m_j}{dt} = \frac{d(\xi_j V)}{dt} = V\frac{d\xi_j}{dt} + \xi_j\frac{dV}{dt}
= V\frac{d\xi_j}{dt} + \xi_j V\,\theta,
$$

so the *concentration* rate that keeps cell mass on its budget is

$$
\frac{d\xi_j}{dt} = \frac{1}{V}\frac{dm_j}{dt} - \xi_j\,\theta
= S_j^{\mathrm{rxn}} - \frac{[J_jA]_{\text{faces}}}{V} - \xi_j\,\theta.
\tag{10.14}
$$

The first two terms are the reaction source and (for gases) the flux divergence already in the species RHS; the last term, $-\xi_j\theta$, is the **dilation correction**. It says: when a cell shrinks ($\theta<0$), the concentration of every component it carries must *rise* (mass conserved in less volume), so we *subtract* a negative number, increasing $d\xi_j/dt$; when a cell swells ($\theta>0$), the concentration falls. For the pore gas this is the compression of the gas in the shrinking pores, which the pressure closure $P = N/\varphi$ (§9.4) then turns into a Darcy response:

$$
\boxed{\;\left.\frac{d\xi_j}{dt}\right|_{\mathrm{dilation}} = -\,\xi_j\,\theta_i,
\qquad \text{every component } j\;}
\tag{10.15}
$$

> **Implementation note.** `_apply_ale_dilation_correction!(dξ, θ, state, mesh)` (`src/Residual/residual.jl`) loops over cells and components, reads the per-cell strain `θ_i = θ[i]` (the `θ_out` array filled during the mesh-velocity computation, §10.6.2), and executes `dξ[j, i] -= state.ξ[j, i] * θ_i` for every component. The sign is exactly the $-\xi_j\theta_i$ of (10.15). It is applied *after* the reaction sources and the flux divergences, so it sits on top of the already-assembled species RHS, and it is the only mesh-motion term in the species equations.

### 10.7.3 Why the strain rate, not the mesh velocity

The volume rate of a cell is $\dot V/V = (w_R - w_L)/\Delta z + \theta_A$, which equals $\theta_i$ in production by construction of the mesh velocity (§10.6.2, §11.3). The correction is written with the cell's own strain rate $\theta_i$ rather than with the node velocities because $\theta_i$ is local to the cell — the structured Jacobian then carries the dilation as a cell-local block (§15.8.4) — and because under the verification hook `Workspace.w_override` the mesh is driven by a prescribed velocity while $\theta_i$ keeps its reaction-driven value, so a reaction-free prescribed compression leaves every concentration unchanged (free-stream preservation, V-MMS-02).

---

## 10.8 Geometry coupling and observable consequences

### 10.8.1 Cell volume from state

In ALE mode every cell volume is recomputed from the nodal positions and the current cross-section on each residual evaluation:

$$
V_i \;=\; \lvert z_{i+1}-z_i\rvert\,A,
\qquad A \equiv A(t)=A_0\,\mathrm{law}(\bar\chi).
\tag{10.16}
$$

> **Implementation note.** `update_state_geometry!(cache, mesh, ::ALE)` (`src/Problem/state_cache.jl`) sets `cache.volume[i] = abs(z_right - z_left) * A` with `A = cache.A_section[]`, and face areas `cache.A_face[f] = A`, cell centers `(z_left + z_right)/2`. All quantities preserve `eltype(cache.z)` so AD passes carry geometry sensitivities (e.g. $\partial V_i/\partial\chi$ via $A$) through into every kernel that reads volumes or face areas. This is the link by which a perturbation in $\chi$ propagates to the energy and species balances: a larger $A$ enlarges every $V_i$ and every face area, rescaling divergences and the cell-integrated masses. The absolute value $\lvert\cdot\rvert$ guards against transient node inversion during a Newton/AD probe.

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

with total length $N_{\mathrm{state}} = n(1+N_C) + \delta_{\mathrm{ALE}}\,n_{\mathrm{nodes}} + \delta_{\mathrm{WithChi}}\,n$. The $z$ block is present only in ALE mode and the $\chi$ block only in `WithChi` mode, which requires ALE: the fixed-mesh residual evolves neither $\chi$ nor $A$, so `to_problem_def` rejects a lateral shrinkage law with `use_ale = false`. For example, with $N_C=2$ components, $n=5$ cells (hence $n_{\mathrm{nodes}}=6$), the ALE+WithChi length is $5\cdot 3 + 6 + 5 = 26$, which the test `test/geometry/test_variable_cross_section.jl` asserts as `5*(2+3)+1`; the same problem with `lateral=false` (NoChi) has length $5\cdot 3 + 6 = 21$ asserted as `5*3+6`. The block accessors `temperature_range`, `concentration_range`, `z_range`, `chi_range` (`src/Problem/state_layout.jl`) return empty ranges for absent blocks, and the residual's `_pack_chi_derivative_dispatched!` is a no-op when `chi_range` is empty. The `WithChi`/`NoChi` axis is one of the four bundled `SimulationMode` traits; the detailed layout descriptor `StateLayout` and its role in the structured Jacobian are covered in §13 and §15.1.

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

The following limitations are documented in the source and should guide use and future work:

1. **Density time-averaging in $\theta$.** Equation (10.2) uses $\rho_j(T)$ evaluated at the *current* cell temperature, and the density-rate term $\xi_j\,d(1/\rho_j)/dt$ is neglected (§10.1.1). For materials with strongly temperature-dependent pure-phase densities and steep thermal gradients this is an approximation; whether an implicit or time-averaged density would materially change $\theta$ has not been quantified.

2. **Empirical lateral law.** $\mathrm{law}(\bar\chi)$ is user-specified with no first-principles closure (§10.5.3). A percolation- or homogenization-based predictor of the functional form would remove a calibration degree of freedom, especially for intumescent coatings whose expansion factor can exceed an order of magnitude. The law must also be *smooth*: `_compute_theta_A` differentiates it by `ForwardDiff`, so a non-smooth (e.g. piecewise-linear) law gives undefined $\theta_A$ at the kinks and degrades the Jacobian/sensitivity quality (see §17.7).

3. **Vanishing cross-section.** The radial strain rate (10.10) carries a $1/\mathrm{law}(\bar\chi)$ that diverges as the cross-section collapses; the code guards with `law(χ̄) > eps()` (returning $\theta_A=0$ below the floor), but there is no *physical* lower bound on shrinkage — a law allowed to reach zero is unphysical and the guard merely prevents a numerical blow-up. Users should design $\mathrm{law}$ so it stays bounded away from zero.

4. **$\bar\chi$ weighting choice.** The column average (10.7) is mass-weighted by *initial* dry-solid mass. Strongly non-uniform decomposition (e.g. a thin, fully-charred surface over an intact interior) is reduced to a single scalar $\bar\chi$, so the lateral law cannot represent depth-dependent cross-section change. A different weighting (or a depth-resolved area) would be required for that regime.

5. **Single-reactant kinetics.** Volume change inherits the single-reactant restriction of the kinetics (§6.5, §6.10): `species_source_terms!` enforces $N_{\mathrm{Rin}}=1$, so $\theta$ and $\chi$ are built from single-reactant sources only. Bimolecular or multi-reactant stoichiometry (e.g. char oxidation consuming both char and O₂) would require extending the source-term machinery before $\theta$/$\chi$ could account for it.

6. **Interaction with depletion parking.** The $\chi$ block, the per-cell strain $\theta_i$, and $m_{\mathrm{dry},i}$ must all be carried through a park event. The energy-conserving merge is detailed in §14.3; the consistency of $\chi$ accumulation across a park (summing $m_{\mathrm{dry}}$ denominators and weighting the merged $\chi$ by them) is a subtle invariant that the park primitive preserves (§14.3.4).

7. **Substring-based water-vapor exclusion.** The dry-gas index set relies on a name substring match (§10.3.3), which is fragile to naming. A type-level or explicit-flag mechanism for "this gas is moisture-derived" would be more robust than the current `"H2O"`/`"water"` heuristic.

8. **Axial/radial split is sequential, not unified.** The axial (mesh-motion) and radial (cross-section) channels are computed in sequence with $\theta_A$ subtracted from each cell's $\theta_i$; a fully unified space-time ALE that handles both simultaneously remains future work (see §11.9).

9. **Mesh motion is solid motion (solid-Lagrangian ALE).** The volume-change pipeline of this chapter is what *defines* the mesh velocity: $w$ is the prefix sum of the reaction-driven strain rate $\theta$ (§11.3), and every component rides the mesh with no discrete advection operator (§10.7, §11.1.2). Mesh motion decoupled from the solid material velocity — remeshing, smoothing, prescribed motion, refinement under ALE — is not representable without adding a species advection operator; see §11.9 item 5 for the full statement of this capability boundary and its consequence for MMS verification (the spatially-uniform-solid rule under `w_override`).


---

# 11. The Arbitrary Lagrangian–Eulerian (ALE) Framework

Pyrolysis is, intrinsically, a *moving-domain* problem: the condensed phase recedes
as it gasifies, intumescent chars swell, and the position of the exposed surface is
itself an unknown. This chapter develops the Arbitrary Lagrangian–Eulerian (ALE)
formulation that Pyrolysis.jl uses to track that motion without solving a moving-
boundary PDE. We derive the ALE time-derivative identity from the kinematics of a
moving sample point, specialize it to the energy and species equations (every
component rides the skeleton, so the mesh frame adds only the dilation), construct the
cumulative mesh velocity from the per-cell volumetric strain rate, show that the
discrete Geometric Conservation Law (GCL) holds identically because every cell volume
is recomputed from the node positions in the state (including with a time-varying
cross-section), and document the TVD-limited upwind discretization of the
mesh-relative temperature advection used by the verification hook. Throughout, the
governing relations are tied line-by-line to the implementation in
`src/Physics/volume_change.jl`, `src/Problem/state_cache.jl`,
`src/Discretization/ale_fluxes.jl`, and `src/Residual/residual.jl`, and every sign
convention is verified against the code.

The coordinate convention is fixed by the nomenclature (Chapter 2): $z$ is the
through-thickness axial coordinate, $z=0$ is the bottom/substrate (boundary id 2,
`:bottom`), $z=L$ is the top/exposed surface (boundary id 1, `:top`), heat enters
from the top, and the material shrinks downward. Cell index $i$ increases from
bottom to top; the *left* face/node of a cell is at lower $z$, the *right* at higher
$z$. The flux sign convention is **positive flux = transport in the $+z$
direction**. In this chapter $i$ denotes the **cell index** and $r$ the reaction
index (overload note I1).

---

## 11.1 Motivation: tracking a receding/expanding surface

### 11.1.1 The moving-boundary problem and three classical responses

A pyrolyzing slab loses mass through gasification reactions of the form
$\text{solid} \to \nu\,\text{char} + (1-\nu)\,\text{gas}$. The condensed material
that supplies the gas physically disappears, so the geometric thickness of the
column,
$$
L(t) = z_{\text{top}}(t) - z_{\text{bottom}} ,
$$
decreases (or, for intumescent systems, increases) in time. Three numerical
strategies exist:

1. **Fixed Eulerian grid with a sharp moving front.** Keep the mesh stationary and
   track the surface as an interface crossing cells. This requires either a
   level-set / VOF reconstruction or front-tracking, both of which inject
   discontinuous coefficients into the energy and species operators. In 1D with a
   single receding surface this is overkill and produces stair-step mass-loss
   artifacts.

2. **Pure Lagrangian grid.** Pin every mesh node to a material particle. Mass is
   then trivially conserved per cell, but as the solid is consumed the cells near the
   surface collapse, the time step crashes against the CFL limit of the smallest
   cell, and the mesh tangles ($V_i \le 0$). For deep burns this is fatal.

3. **Arbitrary Lagrangian–Eulerian (ALE).** Let the mesh move with a velocity
   $w(z,t)$ that is *neither* the material velocity (pure Lagrangian) *nor* zero
   (pure Eulerian), but is chosen to keep the mesh well-conditioned while still
   tracking the receding surface. Pyrolysis.jl adopts this approach.

The defining freedom of ALE is that the mesh velocity $w$ is a *modeling choice*,
constrained only by boundary conditions and the requirement that the mesh follow the
domain boundaries. Pyrolysis.jl makes the natural choice for a 1D consumed-solid
problem: **the mesh velocity is the cumulative integral of the local volumetric
strain rate from the pinned substrate upward**, so that each node moves exactly fast
enough to absorb the local volume change of the material below it (§11.3). The
result is a mesh that shrinks with the material — its nodes ride on the deforming
solid skeleton — yet is rebuilt analytically each step so it never tangles unless
the physical domain itself collapses.

### 11.1.2 Every component rides the skeleton

The frame of each component follows from the frame of its transport closure:

- **Solid and liquid components** are bound to the deforming matrix skeleton. Their
  material velocity equals the mesh velocity, $\mathbf v_{\text{solid}} = \mathbf w$.
  In the mesh frame they do not advect; they are **Lagrangian**. Their mass per cell
  is locked to the mesh, so when the mesh shrinks their *concentration* $\xi_j$ rises
  even with no reaction — the "dilation" effect handled by a dedicated correction
  (§11.3.4).

- **Gas components** percolate through the porous skeleton. Darcy's law is a drag
  balance between the gas and the pore walls, and the volume-fraction-gradient
  diffusion law is likewise written relative to the pore walls, so the flux the
  closures return, $\mathbf J_g = \xi_g(\mathbf v_{\text{gas}} - \mathbf v_{\text{solid}})$,
  is the gas flux **relative to the skeleton** (§3.2.3, §9). The laboratory-frame gas
  flux is $\xi_g\mathbf v_{\text{solid}} + \mathbf J_g$: the pore gas is carried with
  the pores it occupies and slips through them only by $\mathbf J_g$. In the mesh
  frame the gas is therefore **also Lagrangian** apart from $\mathbf J_g$, which the
  Eulerian flux divergence already handles. It receives the same dilation
  correction as the condensed components and no advection term.

This is the frame ThermaKin writes (its Eq. 15 applies one deformation term to
every component, gas included) and the frame Gpyro writes (its gas balance, Eq. 58,
lives on control volumes that deform with the solid, with only $\dot m''$ crossing
their faces). Treating the gas as fixed in the laboratory instead would count the
skeleton's motion twice: once inside the frame of $\mathbf J_g$ and once as a sweep
of the gas by the moving mesh.

> **Formulation-level fact: this is a *solid-Lagrangian* ALE.** The identity
> $\mathbf v_{\text{solid}} = \mathbf w$ is not merely an assumption about the
> physics — it is built into the discretization. The mesh velocity is
> *computed from* the solid volume-change rate (the $\theta$ prefix sum,
> §11.3), and the residual gives **no component** a discrete advection
> operator: the only mesh-motion term in the species equations is the
> dilation $-\xi_j\theta_i$ applied to every component. The temperature
> equation follows the same logic: the thermal mass $\rho c_p^{\text{eff}}$
> (every component) rides the mesh, so **production ALE has no advection
> term in any equation** (§11.2.2). The "arbitrary" in ALE is therefore
> constrained: any mesh motion **decoupled from the solid material
> velocity** — remeshing, mesh smoothing, prescribed node motion,
> refinement inside ALE — is *not representable* without first adding an advection
> operator for the slip velocity $w - v_{\text{solid}}$. The verification
> hook `Workspace.w_override` (which replaces the reaction-driven $w$ with a
> prescribed $w(z,t)$ for MMS runs, §18.2) inherits exactly this constraint:
> under prescribed non-material motion a spatially varying lab-frame-static
> field has a $w\,\partial_z\xi$ term with no discrete operator to represent
> it, so the manufactured-solution driver assembles every species source in
> the material-following frame (see test/verification/exact/mms.jl and
> §11.9 item 5). The temperature field, by contrast, *does* get its slip
> advection under the override — $T$ is advected by
> $w_{\text{applied}} - w_{\text{physical}}$, which is non-zero exactly and
> only there (§11.2.2).

---

## 11.2 The ALE time-derivative identity

### 11.2.1 Kinematics of a moving sample point

Let $f(z,t)$ be any field (temperature, concentration). Introduce a moving sample
point whose position $z=Z(\chi,t)$ is labeled by a material-like coordinate $\chi$
that is constant along the mesh trajectory, with mesh velocity
$$
w \;=\; \left.\frac{\partial Z}{\partial t}\right|_{\chi}.
$$
The value of $f$ seen by an observer riding on that point is $\hat f(\chi,t) =
f(Z(\chi,t),t)$. Differentiating with the chain rule,
$$
\left.\frac{\partial f}{\partial t}\right|_{\chi}
\;=\;
\left.\frac{\partial f}{\partial t}\right|_{z}
\;+\;
\left.\frac{\partial f}{\partial z}\right|_{t}
\left.\frac{\partial Z}{\partial t}\right|_{\chi},
$$
which is the **ALE time-derivative identity**
$$
\boxed{\;
\left.\frac{\partial f}{\partial t}\right|_{\chi}
=
\left.\frac{\partial f}{\partial t}\right|_{x}
+ \mathbf w\cdot\nabla f
\;}
\tag{11.1}
$$
(`ale_fluxes.jl`). Here $\partial/\partial t|_{\chi}$ is the **ALE
(mesh-point)** rate — the rate at which $f$ changes for a node carried by the mesh —
and $\partial/\partial t|_{x}$ is the **Eulerian (fixed-point)** rate that the
underlying conservation laws of Chapter 3 are written in. The advection correction
$\mathbf w\cdot\nabla f$ converts between the two. In 1D, $\mathbf w\cdot\nabla f =
w\,\partial f/\partial z$.

The integrator advances the *mesh-frame* state — the node positions $z_i$ and the
cell-centered $(T_i,\xi_{j,i})$ that ride with the mesh — so its right-hand side is
$\partial/\partial t|_{\chi}$. The physics modules, however, compute the Eulerian
right-hand side $\partial/\partial t|_{x}$ (conduction divergence, reaction sources,
gas-flux divergence). Identity (11.1) tells us exactly what to add:
$$
\underbrace{\left.\frac{\partial f}{\partial t}\right|_{\chi}}_{\text{what the solver integrates}}
=
\underbrace{\left.\frac{\partial f}{\partial t}\right|_{x}}_{\text{Eulerian RHS from physics}}
+ \underbrace{w\,\frac{\partial f}{\partial z}}_{\text{ALE advection correction}}.
\tag{11.2}
$$

### 11.2.2 Specialization to the energy equation — the relative-velocity form

The Eulerian energy balance (Chapter 7) reads
$$
\rho c_p^{\text{eff}}\left.\frac{\partial T}{\partial t}\right|_{x}
= -\nabla\cdot\mathbf q + Q,
$$
with $\mathbf q$ the conductive flux and $Q$ the lumped volumetric source
($Q_{\text{rxn}}+Q_{\text{rad}}+S_{\text{conv}}+\dots$; $\rho c_p^{\text{eff}}$ is the
mixture heat capacity, every component included — see §7 and §11.2.4).
Crucially, that balance — like every fixed-frame heat equation with no
solid-convection term — is written for a medium that is **stationary in the lab
frame**. Here it is not: the thermal mass $\rho c_p^{\text{eff}}$ is carried by the
skeleton together with the pore gas it holds, and the skeleton *moves* — at
exactly the mesh velocity, $\mathbf v_{\text{solid}} = \mathbf w$ (§11.1.2). The
correct lab-frame balance
for the moving matrix therefore carries a solid-convection term,
$\rho c_p^{\text{eff}}\,(\partial_t T|_x + \mathbf v_{\text{solid}}\cdot\nabla T)
= -\nabla\cdot\mathbf q + Q$, and converting to the mesh frame with (11.1) leaves
the **relative-velocity** correction:
$$
\boxed{\;
\left.\frac{\partial T}{\partial t}\right|_{\chi}
= \underbrace{\frac{1}{\rho c_p^{\text{eff}}}\!\left(-\nabla\cdot\mathbf q + Q\right)}_{\text{Eulerian RHS}}
+\; \left(\mathbf w - \mathbf v_{\text{solid}}\right)\cdot\nabla T
\;=\; \text{Eulerian RHS} \quad (\text{production: } \mathbf v_{\text{solid}} = \mathbf w)
\;}
\tag{11.3}
$$
This is the exact energy-equation counterpart of the species logic in §11.2.3:
the advecting velocity is always the velocity of the mesh *relative to the
material that carries the field's capacity*, and for the mesh-bound matrix that
relative velocity vanishes identically. **Production ALE therefore applies no
temperature advection at all** — the mesh-frame rate of the matrix-carried
temperature *is* the material rate. (Temperature is intensive, so there is also
no dilation contribution: a static temperature field is unchanged by uniform
stretching of the mesh.)

> **Why the distinction matters.** Applying the *full* chain-rule term
> $+w\,\partial T/\partial z$ from (11.2) here would implicitly treat the
> solid as lab-stationary while §11.2.3 simultaneously treats it as
> mesh-bound — two mutually inconsistent statements. The full term advects
> the temperature field relative to the very material carrying it: a
> spurious energy sink $\int \rho c_p\, w\,\partial_z T\,dz \approx \dot
> m''\,c_s\,(T_s - \bar T_{\text{form}})$ — a few percent of the effective
> heat of gasification in every ALE run with volume change and thermal
> gradients, and a systematic low-MLR bias against regressing-mesh codes
> (FDS/Gpyro/ThermaKin advect enthalpy with the material). The Landau
> ablation case V-CPL-02 verifies the traveling-wave identity that exposes
> such a term, and V-CPL-06b pins the production T-advection at exactly
> zero.

**The `w_override` exception.** Under the MMS verification hook
(`Workspace.w_override`, §18.2) the applied mesh velocity is decoupled from the
material velocity, so the relative velocity
$w_{\text{rel}} = w_{\text{applied}} - w_{\text{physical}}$ is non-zero and the
advection kernel runs with exactly that velocity (`residual.jl`, ALE step 6;
the physical, reaction-driven velocity is saved in `Workspace.w_physical`
before the override is applied). Free-stream and prescribed-motion MMS cases
have $w_{\text{physical}} = 0$, so they exercise the kernel with the full
override velocity — the regression net for the advection kernel itself.

**Sign of the energy advection (override path).** With $w_{\text{rel}}>0$ (mesh
sliding in $+z$ relative to the material) and $\partial T/\partial z>0$, the
term $w_{\text{rel}}\,\partial T/\partial z>0$: the node moves into warmer
material, so the temperature it carries increases — the physically correct
convective interpretation of the chain rule.

### 11.2.3 Specialization to the species equations

For component $j$ with laboratory-frame velocity $\mathbf v_j$ and skeleton-relative
flux $\mathbf J_j = \xi_j(\mathbf v_j - \mathbf v_{\text{solid}})$, the Eulerian species
balance is
$$
\left.\frac{\partial \xi_j}{\partial t}\right|_{x}
+ \nabla\cdot(\mathbf v_{\text{solid}}\, \xi_j)
= -\nabla\cdot\mathbf J_j + S_j ,
$$
with $\mathbf J_j \equiv 0$ for condensed components and $\mathbf J_g$ the Darcy–Fick
flux for gases (§11.1.2). The general ALE form follows from (11.1): substituting
$\partial\xi/\partial t|_x = -\nabla\cdot(\mathbf v_{\text{solid}}\xi_j) - \nabla\cdot\mathbf J_j + S_j$
and expanding $\nabla\cdot(\mathbf v_{\text{solid}}\xi_j) = \mathbf v_{\text{solid}}\cdot\nabla\xi_j + \xi_j\nabla\cdot\mathbf v_{\text{solid}}$ gives
$$
\left.\frac{\partial \xi_j}{\partial t}\right|_{\chi}
= \text{(Eulerian RHS)} - (\mathbf v_{\text{solid}} - \mathbf w)\cdot\nabla\xi_j
\;-\; \xi_j\,\theta,
\tag{11.4}
$$
where the last term, $-\xi_j\theta$, is the **volumetric dilation** of the skeleton,
with $\theta$ the **full** volumetric strain rate $\theta_{\text{cell}}$ — axial
$\partial w/\partial z$ *plus* the lateral channel $\theta_A$. Identifying $\theta$
with $\nabla\cdot\mathbf w$ is correct only without lateral shrinkage: in 1D,
$\partial w/\partial z = \theta - \theta_A \ne \theta$ when $\theta_A \ne 0$, and
applying the full $\theta_{\text{cell}}$ is the unique choice that makes
$d(\xi V)/dt = SV - [JA]_{\text{faces}}$ exact (§11.5.2).

In production the mesh moves with the skeleton, $\mathbf w = \mathbf v_{\text{solid}}$,
so the slip term vanishes for **every** component and (11.4) reduces to
$$
\boxed{\;
\left.\frac{\partial \xi_j}{\partial t}\right|_{\chi}
= S_j - \nabla\cdot\mathbf J_j - \xi_j\,\theta
\qquad (\text{every component; } \mathbf J_j = 0 \text{ for solid/liquid}).
\;}
\tag{11.5}
$$
The species equations therefore carry **no mesh-motion advection**: the Eulerian
RHS supplies $S_j - \nabla\cdot\mathbf J_j$ unchanged (the skeleton-relative flux is
the right flux in the mesh frame) and the residual adds the dilation (11.5) to
every component (§11.3.4). Integrated over a cell this is
$d(\xi_jV)/dt = S_jV - [J_jA]_{\text{faces}}$: the mass of every component in a
cell changes only by reaction and by $\mathbf J_j$ through its faces, so the
column total changes only through the boundary faces, exactly as on a fixed
mesh. In a sealed collapsing column the pore gas is compressed cell by cell
($\xi_g \uparrow$ as $V \downarrow$, V-CPL-03 P4), and under Darcy–Fick transport
the pressure closure turns that compression into a venting response.

**Why not a face-sweep flux for the gas.** Treating the gas as fixed in the
laboratory while the mesh sweeps past it would add $\nabla\cdot(\xi_g\mathbf w)$ to
(11.5). With $\mathbf J_g$ defined relative to the skeleton that term counts the
skeleton's motion a second time: in a uniformly collapsing sealed column it
leaves every interior cell's gas untouched and piles the whole compression into
the top cell, and at an open regressing surface it releases pore gas at
$\xi_g\,w$ that the diffusive and Darcy fluxes already account for.

### 11.2.4 Why the gas is in $\rho c_p^{\text{eff}}$ (cross-reference)

The advection term in (11.3) is the mesh-frame transport of *sensible* enthalpy,
and it vanishes in production because the thermal mass rides the mesh. That
thermal mass is the mixture heat capacity $\rho c_p^{\text{eff}} = \sum_j\xi_jc_{p,j}$
over every component (Chapter 7): the pore gas is carried with the pores, so its
sensible heat moves with the mesh like the matrix's and needs no separate
advection term either. The complementary gas enthalpy carried by the *flux
relative to the skeleton* is accounted for by the volumetric source
$S_{\text{conv}} = -\sum_g c_{p,g}\,\bar J_g\,\partial T/\partial z$ (Chapter 7). This is
the ThermaKin Eq. 17 convention term for term, including its deformation term,
which is the $\mathbf v_{\text{solid}}\cdot\nabla T$ piece the mesh frame absorbs. The
surface balance carries no gas-enthalpy term, which would count the gas heating a
second time (§7.4.3, Chapter 12).

---

## 11.3 Cumulative mesh velocity

The ALE freedom is fixed by *choosing* $w(z,t)$. Pyrolysis.jl chooses it so that the
mesh exactly absorbs the local volume change of the material, with the substrate
pinned. This section derives the cumulative-integral formula, its discrete form, and
the dilation correction that the mesh motion induces in the species equations.

### 11.3.1 The volumetric strain rate $\theta$

For a fixed-mass control volume of the deforming matrix, the relative volume change
rate (volumetric strain rate, dilation) is obtained from the mixture continuity
relation by summing the volume freed/created per unit mass change of each
matrix-contributing component:
$$
\boxed{\;
\theta \;\equiv\; \frac{1}{V}\frac{dV}{dt}
\;=\; \sum_{j=1}^{N_C} \gamma_j\,\frac{1}{\rho_j}\,\frac{d\xi_j}{dt}
\;}
\tag{11.6}
$$
(`volume_change.jl`). Here:

- $\gamma_j$ is the **swelling factor** of component $j$ (component field
  `swelling_factor`): $\gamma=1$ for solids by default (when a solid depletes, the
  matrix loses that volume), $\gamma=0$ for gases by default (escaping gas does not
  change the skeleton volume), and $\gamma\in[0,1]$ for intumescent/partial
  contributors. Only components with $\gamma_j>0$ enter the sum (the implementation
  early-exits on $\gamma_j\le 0$, `volume_change.jl`).
- $\rho_j(T)$ is the **bulk** (pure-phase) density of component $j$, evaluated at the
  cell temperature; the temperature dependence is honored through the property
  function `comp.density(T)` (`volume_change.jl`).
- $d\xi_j/dt$ is the **species reaction source term** $S_j$ — the *reaction* part of
  the species RHS, supplied as the pre-computed `dξ_rxn` matrix
  (`residual.jl`). It is negative for a consumed solid and positive for a
  produced char/gas.

**Physical reading and sign.** When a solid of density $\rho_j$ is consumed
($d\xi_j/dt<0$, $\gamma_j=1$), the term $\gamma_j(1/\rho_j)(d\xi_j/dt)<0$, so
$\theta<0$: the cell shrinks. A net char-producing or swelling reaction with the
appropriate $\gamma$ can yield $\theta>0$ (expansion/intumescence). Note that
(11.6) is the *swelling-factor-weighted* form; the loose statement "$\theta =
-(\sum_j d\xi_j/dt)/\rho_{\text{total}}$" sometimes used as shorthand
is only correct in the special case of a single solid with
$\gamma=1$ and constant $\rho$ — the code always uses (11.6).

> **Implementation note.** Relative strain rate: `volume_change_rate` /
> `_volume_change_rate` (dispatch-unrolled for $N_C\in\{1,2,3\}$, `@generated` for
> $N_C>3$), `src/Physics/volume_change.jl`.

### 11.3.2 Continuous-to-discrete: integrating strain into mesh velocity

The mesh velocity must reproduce the volume change of every cell below a given node.
In the continuous picture the velocity at height $z$ is the integral of the
*axial* strain rate from the pinned bottom:
$$
w(z,t) = \int_0^{z}\big(\theta(z',t) - \theta_A(t)\big)\,dz',
\qquad w(0,t)=0,
\tag{11.7}
$$
where $\theta_A$ is the **column-global radial channel** subtracted to avoid
double-counting any cross-section change (§11.3.3); for constant cross-section
$\theta_A\equiv 0$ and $w$ integrates the full local strain. The axial residual
strain is $\theta_L = \theta - \theta_A$ (`_compute_ale_mesh_velocity!`,
`residual.jl`). Discretizing the integral cell-by-cell with the bottom node
pinned gives the **cumulative mesh-velocity recursion**
$$
\boxed{\;
w_1 = 0,\qquad
w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i
\;}
\tag{11.8}
$$
with cell thickness $\Delta z_i = V_i/A_{\text{section}}$
(`residual.jl`). Equivalently, $w$ at the right node of cell $i$ is the
exclusive prefix sum
$$
w_{i+1} = \sum_{k=1}^{i}(\theta_k - \theta_A)\,\Delta z_k
= \sum_{k\le i}\theta_{L,k}\,\Delta z_k .
\tag{11.9}
$$

The bottom-up sweep (loop `cell_idx = 1:n_cells`, `residual.jl`) starts at
`w[1]=0` (substrate pinned, `residual.jl`) and accumulates
`cumulative_velocity += θ_L * Δz_cell`, assigning the running total to the cell's
right node (`residual.jl`). Because the strain rates of all cells below a
node are summed, the top node moves at $w_{n+1} = \sum_k \theta_{L,k}\Delta z_k$,
i.e. the surface velocity equals the total axial volume-change rate of the column,
exactly tracking the receding (or expanding) surface.

**Substrate boundary condition.** Pinning $w_1=0$ enforces zero normal mesh velocity
at the impermeable substrate: the bottom node never moves, all motion is taken up
above it, and the material "stacks up" toward the surface as it is consumed.
This is the natural choice
for a slab on a non-receding backing.

### 11.3.3 The column-global radial channel $\theta_A$ (WithChi)

When the material carries a *lateral shrinkage law* (variable cross-section,
`WithChi` mode, Chapter 10), the column cross-section evolves as $A(t) =
A_0\,\mathrm{law}(\bar\chi)$ with $\bar\chi$ the mass-weighted column-average
pyrolysis progress. Part of the total volume change is then radial (cross-section)
and must be removed from the *axial* mesh-velocity integral so it is not counted
twice. The column-global radial volumetric strain rate is, by the chain rule,
$$
\boxed{\;
\theta_A = \frac{1}{A}\frac{dA}{dt}
= \frac{1}{\mathrm{law}(\bar\chi)}\,\mathrm{law}'(\bar\chi)\,\frac{d\bar\chi}{dt}
\;}
\tag{11.10}
$$
(`residual.jl`, `_compute_theta_A`). The column-average progress rate is the
volume-weighted dry-gas production normalized by the total initial dry-solid mass,
$$
\frac{d\bar\chi}{dt} = \frac{1}{M_{\text{dry}}}\sum_{i=1}^{n}
\Big(\sum_{j\in\text{dry gas}} S_{j,i}\Big)\,V_i ,
\qquad
M_{\text{dry}} = \sum_i m_{\text{dry},i},
\tag{11.11}
$$
computed at `residual.jl`, with $\mathrm{law}'(\bar\chi)$ obtained by
automatic differentiation `ForwardDiff.derivative(law, χ̄)` (`residual.jl`).
Guards: $\theta_A=0$ for identity-law materials, for missing `dξ_rxn`, for
$M_{\text{dry}}\le\epsilon$, and for $\mathrm{law}(\bar\chi)\le\epsilon$ (preventing
division by a vanishing cross-section, `residual.jl`). For the default
constant-cross-section material the dispatch returns exactly `0.0`
(`residual.jl`, `NoChi`), and (11.8) reduces to the pure cumulative-strain
integral. The interaction of $\theta_A$ with $\chi$-block dynamics and the lateral
law is developed fully in Chapter 10; here it enters only as the scalar subtracted
in (11.8).

### 11.3.4 Cell-centered velocity and the dilation correction

The verification-only temperature advection (§11.2.2, §11.6) needs a *cell-centered*
mesh velocity. It is the arithmetic mean of the bounding node velocities:
$$
\boxed{\;
w_{\text{cell},i} = \tfrac12\big(w_{\text{left}} + w_{\text{right}}\big)
= \tfrac12\big(w_i + w_{i+1}\big)
\;}
\tag{11.12}
$$
(`compute_ale_advection_term!`, `ale_fluxes.jl`).

The **dilation correction** implements the mesh-motion term of (11.5) for every
component. For each cell $i$ and each component $j$,
$$
\boxed{\;
\left.\frac{d\xi_j}{dt}\right|_{\text{dilation}} = -\,\xi_{j,i}\,\theta_i
\qquad (\text{every } j)
\;}
\tag{11.13}
$$
(`_apply_ale_dilation_correction!`, `src/Residual/residual.jl`). Physically: in the
mesh-following frame, a shrinking cell ($\theta_i<0$) must *concentrate* the mass it
carries, so subtracting $\xi_{j,i}\theta_i$ (a positive amount when $\theta_i<0$)
raises $\dot\xi_j$; conversely an expanding cell dilutes it. Without this term every
species would be spuriously depleted by the mesh motion. The pore gas is included
because it is carried with the pores it occupies (§11.1.2); its compression under
collapse is what the pressure closure sees.

The combined ALE species update assembled into the residual is therefore
$$
\left.\frac{d\xi_{j,i}}{dt}\right|_{\chi}
=
\underbrace{-\frac{1}{V_i}(J_{j,R}-J_{j,L})A + S_{j,i}}_{\text{Eulerian RHS}}
\;-\;
\underbrace{\xi_{j,i}\,\theta_i}_{\text{dilation (every component)}} ,
\tag{11.14}
$$
and the energy update is (11.3) with the divergence and sources from Chapter 7.

> **Implementation note.** The full ALE residual assembly order is: populate
> typed state cache → update $A(t)$ from $\bar\chi$ → recompute state geometry from
> node positions → properties → Eulerian physics RHS (`dξ_rxn` cached) →
> $\theta_A$ → cumulative mesh velocity → (override only) mesh-relative T
> advection → dilation correction → pack `du` (incl. $du[z_{\text{off}}+i]=w_i$)
> → pack $\chi$ derivatives. `_residual!(…, ::ALE)`, `src/Residual/residual.jl`.

---

## 11.4 Node-position update and geometric volume recompute

In ALE mode the node positions $z_i$ are part of the ODE state vector (Chapter 13;
`StateLayout`), packed in the $z$-block at offset $n(1+N_C)$.
Their evolution equation is simply the mesh velocity,
$$
\boxed{\;
\frac{dz_i}{dt} = w_i
\;}
\tag{11.15}
$$
written into the residual as `du[zoff+i] = ws.w_work[i]` (`residual.jl`).
The time integrator (Chapter 15) advances $z_i$ together with $(T,\xi,\chi)$ as one
coupled stiff system; the structured Jacobian carries the $z$-*column* couplings
in the geometry AD block and the $z$-*row* dependence on $(T,\xi)$ as embedded
prefix-sum couplings (plus the WithChi $\theta_A$ rank-1 channels) emitted by
the `ALEMeshVelocityOperator` (Chapter 15, §15.8.4 and §15.11.3).

### 11.4.1 Geometric volume recompute

After the nodes are at their new positions, all cell volumes are recomputed
**analytically from geometry**, not carried forward by a flux formula:
$$
\boxed{\;
V_i^{\text{new}} = A(t)\,\big|z_{\text{right}} - z_{\text{left}}\big|
= A(t)\,|z_{i+1}-z_i|
\;}
\tag{11.16}
$$
with the current (possibly $\chi$-updated) cross-section $A(t)$. Cell centers are the
node midpoints,
$$
z_i^c = \tfrac12\big(z_i + z_{i+1}\big),
\tag{11.17}
$$
and in 1D all face areas equal the global cross-section, $A_f = A(t)$. During
residual evaluation this geometry is rebuilt directly from the state-vector node
positions, type-stably (carrying `eltype(u)` so AD propagates through geometry), by
`update_state_geometry!(…, ::ALE)` (`src/Problem/state_cache.jl`), after
`update_cross_section_area!` has set $A(t)$. This keeps the residual a pure function
of $u$ that never mutates the mesh; no separate volume is evolved, and the residual
never reads the mesh's mirrored `cell_volumes`.

**Why geometric recompute rather than a flux update?** Recomputing $V = A(t)\,L$
from current node positions is robust to a simultaneous change in the cross-section
$A(t)$ (variable-area mode): carrying $V_{\text{old}}$ across an $A(t)$ change via a
pure axial-flux formula would freeze the old radial area into the volume. The
geometric path bakes the current $A(t)$ into every volume automatically. As shown
in §11.5, this path satisfies the GCL as an algebraic identity.

---

## 11.5 The Geometric Conservation Law (GCL)

### 11.5.1 Statement and discrete form

The Geometric Conservation Law is the requirement that a *uniform* field be
preserved exactly under mesh motion — equivalently, that the discrete rate of
change of each cell volume equal the net flux of mesh velocity through its faces. In
continuous form,
$$
\frac{dV_i}{dt} = \oint_{\partial\Omega_i}\mathbf w\cdot\mathbf n\,dA .
$$
In 1D, cell $i$ has a left face (outward normal $-\hat z$, velocity $w_{\text{left}}$)
and a right face (outward normal $+\hat z$, velocity $w_{\text{right}}$), both of
area $A$, so
$$
\frac{dV_i}{dt} = A\,w_{\text{right}} - A\,w_{\text{left}} = A\,(w_{\text{right}} - w_{\text{left}})
\tag{11.18}
$$
for a constant cross-section. Why this matters: a GCL violation $\delta V$ acts as a
*spurious mass and energy source*. If the discrete volume does not change by exactly
the swept volume, the cell appears to gain or lose conserved quantity from nothing,
and even tiny violations accumulate over thousands of time steps until they can
dominate the genuine reaction-driven mass loss. Exact GCL is therefore
non-negotiable for long-time accuracy.

### 11.5.2 Exactness of the geometric path

Pyrolysis.jl never evolves a cell volume. The node positions are state (11.15) and
every residual evaluation recomputes $V_i = A(t)\,|z_{i+1}-z_i|$ from them (11.16);
every volume the residual uses (the divergences, the face distances, the
$\chi$ source) is this one, so there is no second copy that could drift.
Differentiating (11.16) along the semi-discrete trajectory, with $dz_i/dt = w_i$,
$$
\frac{dV_i}{dt} = A\,(w_{i+1} - w_i) + |z_{i+1}-z_i|\,\frac{dA}{dt},
$$
which is (11.18) plus the radial area change at constant length — the GCL for a
time-varying cross-section holds by construction, with no separate correction. With
the cumulative mesh velocity (11.8), $A(w_{i+1}-w_i) = A\,(\theta_i-\theta_A)\,\Delta z_i
= V_i(\theta_i - \theta_A)$, and with $\theta_A = (1/A)\,dA/dt$ (11.10) the area term is
$V_i\theta_A$, so
$$
\boxed{\;
\frac{dV_i}{dt} = V_i\,\theta_i
\;}
$$
exactly. Combined with the dilation correction (11.13), every component then obeys
$d(\xi_{j,i}V_i)/dt = \big(S_{j,i} - (\nabla\!\cdot\mathbf J_j)_i\big)V_i$: a uniform
field with no sources and no fluxes is preserved exactly, and a cell's mass changes
only by reaction and by the skeleton-relative flux through its faces. After time
discretization the integrator advances $z$ (hence $V$) and $\xi$ as one coupled
system, so the identity holds to the integrator's local error with no additional
geometric error term. The induced-motion MMS case V-MMS-04, whose mesh velocity is
reaction-driven, verifies it; the free-stream and prescribed-motion cases V-MMS-02
and V-MMS-03 (§18.2) drive the mesh through `Workspace.w_override`, where
$\theta_i$ keeps its reaction-driven value, so they check the geometric path (a
uniform field stays uniform, smooth fields converge at the design order) rather
than this identity.

---

## 11.6 ALE advection discretization

One kernel implements a mesh-motion advection, and only under the verification
hook: **energy** (`compute_ale_advection_term!`), the pointwise correction
$w_{\text{rel}}\,\partial T/\partial z$ with a donor-side-upwinded, TVD-limited,
three-point stencil — derived in this subsection. Since the relative-velocity
form of §11.2.2 makes $w_{\text{rel}} \equiv 0$ in production, this kernel fires
**only under `Workspace.w_override`** (MMS verification); the derivation below
applies with $w = w_{\text{rel}}$. The species equations have no advection kernel
(§11.2.3).

### 11.6.1 Upwinding on the sign of the apparent transport velocity $c=-w$

The advection correction is *added* to the Eulerian RHS, so the mesh-frame
equation reads
$$
\frac{\partial f}{\partial \tau} - w\,\frac{\partial f}{\partial z} = \text{RHS},
$$
i.e. in standard form $\partial f/\partial\tau + c\,\partial f/\partial z = \text{RHS}$
with apparent transport velocity $\boxed{c = -w_{\text{cell}}}$: the lab-frame
field drifts at $-w$ relative to the moving mesh. The stable donor (upwind) side
follows the sign of $c$, **not** the sign of $w$ — a von Neumann analysis of the
$w$-sign choice gives $\mathrm{Re}\,\lambda = (|w|/\Delta z)(1-\cos k\Delta z) > 0$
(growth), equivalent to an injected anti-diffusion coefficient $-|w|\Delta z/2$.

For each cell $i$, the cell-centered mesh velocity $w_{\text{cell}}$ (11.12) is
computed; if $|w_{\text{cell}}|<10^{-20}$ the advection is set to zero. Otherwise:

- **$w_{\text{cell}}>0$** ($c<0$): the donor is the **right** neighbor, forward
  difference,
  $$
  \left.\frac{\partial f}{\partial z}\right|_i^{(1)} = \frac{f_{\text{right}} - f_i}{d_{i\rightarrow\text{right}}},
  $$
  with the donor-donor cell taken further right.

- **$w_{\text{cell}}<0$** ($c>0$, the standard shrinking case): the donor is the
  **left** neighbor, backward difference,
  $$
  \left.\frac{\partial f}{\partial z}\right|_i^{(1)} = \frac{f_i - f_{\text{left}}}{d_{i\leftarrow\text{left}}},
  $$
  donor-donor further left. In particular the receding-surface cell takes its
  gradient from the interior — where $\nabla T$ is steepest — rather than from
  the zero-gradient boundary ghost.

Both forms are written so the gradient is expressed in the $+z$ direction (positive
when $f$ increases with $z$), and the advection term is then
$$
\text{advection}[i] = w_{\text{cell}}\,\left.\frac{\partial f}{\partial z}\right|_i
\tag{11.19}
$$
(`compute_ale_advection_term!` for $T$ — override-only, with
$w = w_{\text{rel}}$; `src/Discretization/ale_fluxes.jl`).

### 11.6.2 TVD-limited three-point gradient

To recover second-order accuracy on smooth fields while suppressing oscillations
near fronts (e.g. the sharp pyrolysis temperature front), the first-order
donor-side gradient is corrected with a slope limiter. With $\text{grad}_1$ the
first-order donor-side gradient (toward the donor cell) and $\text{grad}_0$ the
next gradient one cell further to the donor side, both in the $+z$ convention,
the limited gradient is
$$
\boxed{\;
\left.\frac{\partial f}{\partial z}\right|_i
= \text{grad}_1 + \beta\,\psi(r)\,\big(\text{grad}_1 - \text{grad}_0\big),
\qquad
r = \frac{\text{grad}_1}{\text{grad}_0},
\qquad
\beta = \frac{d_0}{d_0+d_1}
\;}
\tag{11.20}
$$
(`_limited_gradient`, `src/Discretization/ale_fluxes.jl`), where $d_0$ is the
center-to-center distance of the near slope segment (cell $i$ to the donor) and
$d_1$ that of the far segment (donor to donor-donor). The blend
$\beta = d_0/(d_0+d_1)$ is the distance-correct extrapolation of the near
slope (which sits $d_0/2$ to the donor side of the cell center) back to the
cell center: $\text{grad}_1 - \text{grad}_0 \approx \tfrac12 (d_0+d_1) f''$,
so the correction adds exactly the missing $\tfrac{d_0}{2}f''$. On a uniform
mesh $\beta = \tfrac12$. The ratio is
**near/far**: because every implemented limiter satisfies
$0 \le \psi(r) \le \min(2r, 2)$ (with $\psi(r)=0$ for $r\le 0$ and $\psi(1)=1$),
the correction stays bounded by the one-sided slopes on monotone data and
vanishes at extrema, while $r\to 1$ on smooth fields recovers the second-order
one-sided extrapolation $\text{grad}_1 + \beta(\text{grad}_1-\text{grad}_0)$. The
inverse ratio $\text{grad}_0/\text{grad}_1$ would scale the correction with the
*far* slope: $\Phi(r) = 1 + \tfrac12\psi(r)(1-r) \to -\infty$ as $r\to\infty$,
producing sign-flipped, hull-violating gradients exactly at the leading edge of
steep fronts (for van Leer the sign flips beyond $r = 1+\sqrt2$).
The limiter $\psi(r)$ is one of (Chapter 13; `flux_limiters.jl`):
$$
\psi_{\text{minmod}}(r)=\max(0,\min(1,r)),\quad
\psi_{\text{vanleer}}(r)=\frac{r+|r|}{1+|r|},\quad
\psi_{\text{superbee}}(r)=\max\!\big(0,\min(2r,1),\min(r,2)\big).
$$
The default for ALE advection is **Van Leer** (the `limiter::Limiter = VANLEER`
keyword of `compute_ale_advection_term!`). Guards: if there is no donor-donor cell, or the limiter is
`NONE`, or $|\text{grad}_0|<10^{-20}$ (flat far slope — the ratio is undefined),
the scheme falls back to first-order $\text{grad}_1$.

### 11.6.3 Boundary treatment: zero-gradient ghost and the boundary gradient distance

At a domain boundary the "neighbor" cell does not exist (id $=0$). Two choices make
the stencil well-posed:

1. **Zero-gradient ghost value.** A boundary ghost takes the boundary cell's own
   value, $f_{\text{ghost}} = f_{\text{cell}}$ (`_get_boundary_temperature`,
   `ale_fluxes.jl`). This is a Neumann (no-flux)
   condition for scalar advection — physically appropriate, since there is no
   advective transport across the closed domain edge.

2. **Boundary gradient distance.** When the upwind neighbor is a boundary, the
   distance returned by `_upwind_gradient_distance` is the full cell thickness,
   $$
   d = \frac{V_i}{A_{\text{section}}} = \Delta z_i,
   $$
   (`ale_fluxes.jl`). This is the distance across the upwind face: the
   **center-to-center** spacing $|z_i^c - z_{\text{neighbor}}^c|$ for
   an interior neighbor, and the center-to-ghost-cell spacing
   $\Delta z_i = 2\times(\text{center-to-face})$ for a boundary. Reusing the genuine
   center-to-center distance (rather than a naive cell width) preserves first-order
   consistency on non-uniform meshes — a mismatched width would introduce a
   zeroth-order $d/\Delta z$ bias that does not vanish as $h\to 0$. At a boundary the
   value of $d$ is in any case **immaterial**: the closed edge uses a zero-gradient
   ghost (the upwind value is set equal to the cell's own value), so the first-order
   gradient is identically zero there and the TVD branch is disabled. Consequently
   $d$ divides a zero numerator in the residual, and no residual or Jacobian entry
   depends on it (the structured Jacobian carries no advection operator, §15.8.5;
   override runs use finite-difference Jacobians).

> **Implementation note.** Energy advection (override-only):
> `compute_ale_advection_term!` with the shared limited gradient
> `_limited_gradient`, in `src/Discretization/ale_fluxes.jl`. The species
> dilation is `_apply_ale_dilation_correction!` in `src/Residual/residual.jl`
> (§11.3.4).

---

## 11.7 Mesh safety

ALE motion can in principle tangle the mesh ($V_i\le 0$) if a node overtakes its
neighbor, e.g. under an aggressive time step near complete burnout. Within the
implicit solve the integrator's own step-size control (Chapter 15) is the first line
of defense. With `handle_depletion = true`, depletion parking (Chapter 14) removes a cell that has thinned below a
fraction of its initial thickness before it collapses; its detector
(`find_thin_cells`) throws on a non-positive signed thickness, reporting an inverted
cell as a mesh-motion failure rather than as depletion. The `min_thickness` keyword
of `solve` (default
$10^{-6}$ m) stops the solve through a termination callback when the live thickness
$z_{n_{\text{live}}+1} - z_1$ falls below it (`create_thickness_termination_callback`,
`src/Solver/callbacks.jl`). The mesh-quality metrics of §14.5 are diagnostics that
`solve` does not call.

---

## 11.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

The dynamic-mesh strategies of the three reference pyrolysis solvers map onto the
Pyrolysis.jl ALE framework as follows. On first use we map their notation to ours:
ThermaKin's stoichiometric coefficient $\theta_i^j$ is our $\nu_{i,j}$ (overload
note A2); Gpyro's pre-exponential $Z$ is our $A_i$, its permeability $K$ is our
$\kappa$, its volume fraction $X_\alpha$ is our $v_j$, and its conversion $\alpha$
is our $\alpha_{\text{conv}}$; FDS's absorption $\kappa_s$ is our $\alpha$.

**ThermaKin2Ds — dynamic mesh adjustment.** ThermaKin2Ds (Stoliarov and co-workers)
solves the condensed-phase equations on a mesh that contracts/expands as material is
consumed or swells, with element thicknesses adjusted at each step to follow the
local volume change — conceptually the same as Pyrolysis.jl's cumulative
strain-driven mesh velocity (11.8). Two differences are worth noting. First,
Pyrolysis.jl makes the mesh velocity an explicit *state variable* $z_i$ integrated by
the same stiff implicit solver as $(T,\xi)$, so the mesh deformation is fully coupled
in the Jacobian (Chapter 15) rather than applied as an operator-split mesh-adjustment
step. Second, the volume change in Pyrolysis.jl is governed by the per-component
swelling factor $\gamma_j$ and bulk density $\rho_j$ through (11.6), giving a single
unified rule for shrinkage, swelling, and intumescence, whereas ThermaKin handles
the cross-sectional (radial) channel separately. Pyrolysis.jl mirrors that
separation via $\theta_A$ (11.10): the *axial* mesh motion is the cumulative strain
(11.8), the *radial* cross-section is updated independently through the lateral
shrinkage law $A=A_0\,\mathrm{law}(\bar\chi)$ (Chapter 10), and the two are kept from
double-counting by subtracting $\theta_A$ from each local rate. The
heat-of-reaction sign convention also differs at the publication level: ThermaKin
publishes $h>0$ as exothermic, whereas Pyrolysis.jl's internal storage uses $h>0$
endothermic with $Q_{\text{rxn}}=-\sum_r h_r r_r$ (overload note H1) — a labeling
difference, not a physics difference.

**Gpyro — variable-$\Delta z$ fixed control volumes.** Gpyro (Lautenberger &
Fernandez-Pello) uses a *fixed* (Eulerian) set of control volumes whose thicknesses
$\Delta z$ can be specified non-uniformly, and accounts for shrinkage/swelling
through a volume-fraction bookkeeping rather than by physically moving the grid.
Pyrolysis.jl can reproduce the Eulerian limit (mesh mode `Eulerian`, $w\equiv 0$,
no $z$-block in the state vector), in which case the residual collapses to the
standard fixed-mesh divergence form (Chapter 13) and the ALE dilation/GCL
machinery is entirely inactive. Gpyro's gas balance is written on its deforming
control volumes with only the Darcy flux crossing their faces — the same
skeleton-bound frame as §11.1.2. The ALE mode adds, on top of
Gpyro-style physics, the option of *physically* tracking the surface with a moving
grid, which keeps in-depth resolution concentrated where it is needed without the
non-uniform-$\Delta z$ tuning Gpyro relies on. Both codes share the same
distance-weighted harmonic-mean face conductivity and volume-fraction/partial-
pressure gradient gas-transport formulation (Chapters 7, 9), so on a constant cross-section with $w=0$ the two
discretizations are directly comparable. Gpyro's lateral (in-plane) area changes are
captured by Pyrolysis.jl's $\mathrm{law}(\bar\chi)$ and $\theta_A$ path.

**FDS — fixed solid grid, no domain motion.** The FDS pyrolysis model (McDermott,
McGrattan, et al.) integrates the condensed-phase mass and energy equations on a
*fixed* one-dimensional solid grid and represents shrinkage/swelling through layer
thickness updates and a cell-shrinking algorithm rather than a continuously moving
ALE mesh. FDS carries no in-depth gas, so its energy convention is the
condensed-phase limit of Pyrolysis.jl's (Chapter 7); it shares the in-depth
Beer–Lambert / surface radiation absorption options (Chapter 8). The principal contrast is again the treatment of the moving boundary:
FDS adjusts discrete layer thicknesses, while Pyrolysis.jl's ALE evolves node
positions as continuous state with an exactly enforced GCL. Where FDS uses a fixed
grid, Pyrolysis.jl's `Eulerian` mode is the corresponding regime; where the surface
recession must be tracked with sub-cell fidelity, the ALE mode supplies it.

In all three comparisons the GCL discussion of §11.5 has no direct counterpart in
the fixed-grid codes (a fixed grid trivially satisfies a degenerate GCL with
$w\equiv 0$); it is specific to the continuously moving ALE mesh and is what
guarantees that the moving grid introduces *no* spurious mass or energy.

---

## 11.9 Limitations and open issues

The principal limitations and open issues are:

1. **Dimensionality.** The framework is implemented only for 1D ($D=1$). The mesh
   data structures `MeshTopology{D,…}` and `MeshGeometry{D}` are parameterized by
   dimension but are placeholders for a multi-dimensional extension — the 1D
   residual reads `Mesh1D` and the `StateCache`, not them. A 2D/3D extension would
   require vector node positions, per-face normals and variable areas, the full
   divergence-theorem GCL, and surface-integral fluxes. The current code assumes a
   spatially uniform (though time-varying) cross-section.

2. **Separation of axial and radial motion.** Variable cross-section is handled by
   splitting axial ALE motion (with $\theta_A$ subtracted) from a separate radial
   cross-section update via $\mathrm{law}(\bar\chi)$, requiring careful sequencing
   (compute axial velocity, then update area). A unified space-time ALE treatment of
   both channels simultaneously is not implemented.

3. **Density evaluation in the strain rate.** The volume-change rate (11.6) uses
   $\rho_j(T)$ at the instantaneous cell temperature; whether an implicit or
   time-averaged density would improve consistency near steep thermal transients is
   an open modeling question, as is the extension of the volume-change accounting to
   multi-reactant stoichiometry.

4. **Jacobian profile trade-offs.** The structured Jacobian (Chapter 15) represents
   the cumulative mesh velocity as embedded prefix-sum couplings, the column-global
   $\theta_A$ back-reaction as `EmbeddedRank1Coupling`s, and the $z$-columns as a
   dedicated geometry AD block. Under the `:ale_reduced` profile (used by `ApproxSolve`) every structured
   coupling (the prefix sums and the $\theta_A$ rank-1 channels, $\chi$ column
   included) and the geometry block are *dropped*, trading exact
   coupling for speed; the WithChi$\times$Darcy interaction is exercised mainly with
   Fickian transport. These are solver-accuracy/cost choices rather than physical
   limitations of the ALE formulation itself.

5. **Solid-Lagrangian capability boundary: no general (non-material) mesh
   motion.** As stated in §11.1.2, the formulation identifies mesh velocity
   with solid material velocity and gives no component a discrete advection
   operator. Consequently mesh motion decoupled from the solid — remeshing,
   mesh smoothing/redistribution, prescribed node motion, refinement performed
   *within* the ALE frame — is outside the representable class: a field that
   is static in the laboratory frame but non-uniform in $z$ cannot be held
   stationary under such motion, because its $w\,\partial_z\xi$ contribution
   has no operator. Supporting general ALE would require adding a species
   advection operator on the slip velocity $w - v_{\text{solid}}$ (with its
   own limiter and Jacobian couplings). The verification hook
   `Workspace.w_override` (prescribed $w(z,t)$ for MMS runs, §18.2)
   deliberately does *not* lift this boundary — the manufactured-solution
   driver assembles every species source in the material-following frame so
   the missing operator is never exercised, and the structured geometry
   Jacobian assumes reaction-driven motion, so override runs use
   finite-difference Jacobians.

---

### Cross-references

- The volumetric strain rate $\theta$, swelling factor $\gamma_j$, lateral shrinkage
  law, and the $\chi$-block dynamics are developed in detail in **Chapter 10, Volume
  Change & Lateral Shrinkage**.
- The Eulerian energy assembly, $\rho c_p^{\text{eff}}$ (every component),
  $S_{\text{conv}}$, and harmonic-mean face conductivity are in **Chapter 7, Heat
  Conduction, Convection & Energy Assembly**.
- The volume-fraction/partial-pressure gas flux and Darcy velocity are in **Chapter
  9, Gas-Phase Mass Transport & Pressure**.
- Flux limiters and the finite-volume divergence operator are in **Chapter 13,
  Finite-Volume Discretization**; depletion parking, which changes the live cell
  count under ALE, is in **Chapter 14, Adaptive Mesh Refinement and Surface
  Depletion**.
- The structured Jacobian couplings for ALE mesh velocity, dilation, and the
  geometry block, and the time integration of the coupled $(T,\xi,z,\chi)$ system
  are in **Chapter 15, Temporal Integration & Structured Jacobian**.
- The mass and energy ledgers are in **Chapter 16, Conservation Diagnostics**.


---

# 12. Boundary and Initial Conditions

This chapter specifies how Pyrolysis.jl closes the coupled energy, species, and
gas-pressure equations at the two physical surfaces of the 1D column ($z=0$
substrate / `:bottom`, $z=L$ exposed surface / `:top`) and how it seeds the
initial state. It derives the surface energy balance and the Newton iteration
that solves for the true surface temperature $T_s$ in a cell-centered
finite-volume scheme; documents every concrete thermal, mass, pressure, and
mesh-motion boundary condition with its exact flux law, temperature derivative,
and sign convention as implemented; explains the type-stable
`BoundaryConditionSet` container and the tag$\to$id bridge; and treats the
initial-condition closures and their auto-defaults. Throughout, the equations are
checked term-by-term against the source so that the math and the code agree
exactly — including the boundary-face sign flips, the Kirchhoff emissivity
coupling, the neglected $d\varepsilon/dT$ term, and the absence of a surface
gas-enthalpy term alongside the volumetric advection source $S_{\text{conv}}$.

This chapter conforms to the unified nomenclature of §2, Nomenclature, and uses
the three binding sign conventions: positive flux is transport in the $+z$
direction (bottom$\to$top); the discrete divergence $(F_R-F_L)A/V$ is positive
when the quantity leaves a cell; and for heats of reaction $h>0$ is endothermic
with $Q_{\text{rxn}}=-\sum_r h_r r_r$. In this chapter the index $i$ denotes the
cell index and $j$ the component index, consistent with the discretization
chapters.

---

## 12.1 Overview and design philosophy

### 12.1.1 What a boundary condition must supply

The continuum model of §3, Governing Equations, is a system of conservation laws
for the condensed-phase temperature $T(z,t)$ and the per-component mass
concentrations $\xi_j(z,t)$ (kg·m⁻³), with the gas-phase pressure $P$ as an
algebraic closure in Darcy mode. In flux (divergence) form the spatial operators
are

$$
\rho c_p^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}} + S_{\text{conv}}
   + Q_{\text{rad}} + Q_{\text{rxn}},
\qquad
\frac{\partial \xi_j}{\partial t}
= -\nabla\!\cdot\mathbf{J}_j + S_{j,\text{rxn}} .
$$

A cell-centered finite-volume discretization (see §13, Finite-Volume
Discretization) integrates each equation over a cell $i$ and replaces the
divergence by a difference of face fluxes,

$$
(\nabla\!\cdot \mathbf{F})_i = \frac{(F_{f_R}-F_{f_L})\,A}{V_i}.
$$

For an interior face the flux $F_f$ is a constitutive law evaluated from the two
adjacent cells (Fourier, Fick, Darcy–Fick). For a **boundary face** there is no
neighboring cell, so the flux must instead be supplied by a boundary condition.
The boundary-condition subsystem therefore has exactly one job: for each boundary
face, given the interior state and the time $t$, produce the heat flux
$q_{\text{BC}}$, the species fluxes $J_j$, and (in Darcy mode) the gas velocity
$v_g$ that enter the divergence operator with the correct sign.

There is one subtlety unique to the **thermal** boundary: because temperature is
stored at cell centers and radiative and convective fluxes depend nonlinearly on
the *surface* temperature $T_s$ (which is not a state variable), the code does not
simply evaluate $q_{\text{BC}}$ at the adjacent cell-center temperature. Doing so
produces an $\sim$20% energy imbalance for $T^4$ radiation. Instead it solves a
small surface energy balance for the true $T_s$ (§12.2). The species and pressure
boundaries have no such nonlinearity and are evaluated directly.

### 12.1.2 Type hierarchy and dispatch

Boundary conditions are organized by physical kind, each rooted in an abstract
supertype declared in the `Materials` module:

- `AbstractThermalBC` — energy/heat boundary conditions;
- `AbstractMassBC` — species and pressure-driven (Darcy) boundary conditions;
- `AbstractRadiationBC` — P1 radiative-transfer boundary conditions (in the
  `Experimental` module; see §8, Thermal Radiation).

Concrete types live in the `BoundaryConditions` module
(`src/BoundaryConditions/boundary_conditions.jl`); the thermal flux laws with their
$T_s$ derivatives (`evaluate_bc_flux_and_derivative`) and the surface-temperature
Newton solve live in `src/Discretization/boundary_temperature.jl`, and the mass-BC
flux laws (`mass_bc_flux`) and the diagnostic heat-flux decomposition
(`compute_boundary_heat_fluxes`) in `src/Discretization/boundary_fluxes.jl`. The whole set is stored in
a type-stable, tag-keyed `BoundaryConditionSet` (§12.7) so that the residual
specializes statically with zero dispatch overhead.

> **Implementation note.** Concrete BC types and the `BoundaryConditionSet`
> container are defined in `src/BoundaryConditions/boundary_conditions.jl`; flux
> laws in `src/Discretization/boundary_fluxes.jl` (mass BCs) and, for thermal BCs
> with their $T_s$ derivatives, `src/Discretization/boundary_temperature.jl`; BC application into the flux
> arrays in `src/Discretization/finite_volume.jl` (`apply_thermal_bc_flux!`,
> `apply_mass_bc_flux!`); pressure BCs in `src/Discretization/darcy_velocity.jl`
> (`apply_pressure_bc_to_darcy_velocity!`).

### 12.1.3 The boundary-face sign convention

Two boundaries exist in 1D, with the **binding** id/tag/position mapping
(`src/Geometry/topology.jl`, `src/Geometry/mesh_generation.jl`):

| Tag | `boundary_id` | Position | Face id |
|---|---|---|---|
| `:top` | 1 | $z=L$, exposed surface | last node ($n+1$) |
| `:bottom` | 2 | $z=0$, substrate | first node ($1$) |

`tag_to_boundary_id(:top)=1`, `tag_to_boundary_id(:bottom)=2`,
`boundary_id_to_tag(1)=:top`.

Because the discrete divergence is $(F_{f_R}-F_{f_L})A/V_i$ with the convention
that a positive face flux points in $+z$, the BC application must orient the
computed *into-the-domain* flux correctly for each boundary. A flux
$q_{\text{cond}}$ that carries heat **into** the slab is in the $+z$ direction at
the bottom face (left face of cell 1) but in the $-z$ direction at the top face
(right face of cell $n$). Hence the implementation stores

$$
q[f_{\text{bottom}}] = +q_{\text{conducted}}, \qquad
q[f_{\text{top}}] = -q_{\text{conducted}},
$$

and identically for mass fluxes (`apply_thermal_bc_flux!`, `apply_mass_bc_flux!`,
`finite_volume.jl`). This is not a bug; it is the consistent realization of the signed
divergence. We return to it in §12.2.5 and §12.4.4.

---

## 12.2 The surface energy balance and the Newton solve

The surface temperature $T_s$ is the keystone of the thermal boundary. It is the
temperature of the exposed face itself, distinct from the adjacent cell-center
temperature $T_{\text{cell}}$, and it is determined by an instantaneous balance
between the externally applied flux and the heat conducted into the slab.

### 12.2.1 Continuous statement

Consider the exposed face. An energy balance on an infinitesimally thin control
volume at the surface (zero thermal mass) equates the net flux applied by the
environment to the flux conducted into the material:

$$
q_{\text{BC}}(T_s) \;=\; q_{\text{cond}}(T_s, T_{\text{cell}}),
$$

with, for the most general thermal BC,

$$
q_{\text{BC}}(T_s) = q_{\text{ext}}\;-\;\varepsilon_{\text{eff}}\sigma T_s^4
\;-\;h_{\text{conv}}\,(T_s-T_\infty)
$$

(absorbed external flux, re-radiation loss, convective loss; all signed positive
into the domain) and the conducted flux written as a one-sided Fourier difference
from the surface to the first cell center,

$$
q_{\text{cond}}(T_s,T_{\text{cell}})
= k_{\text{eff}}\,\frac{T_s - T_{\text{cell}}}{\Delta z_{1/2}},
$$

where $\Delta z_{1/2}$ is the distance from the boundary face to the first cell
center and $k_{\text{eff}}$ is the effective conductivity of the boundary cell
(§5, Effective Properties & Mixing Rules). Both sides are positive when heat
flows into the material. The gas enthalpy carried by the flux relative to the
skeleton is transported by the volumetric source $S_{\text{conv}}$ (§7), not by
a surface term; under local thermal equilibrium the gas crosses the surface at
$T_s$ on both sides, so the balance has no gas-enthalpy term (§1.4.3, §7.4.3). The
outward gas flux enters the balance only through the blowing model of a
`ConvectiveBC`, which scales the convective coefficient (§12.3.3).

The conductive distance is exactly half the boundary cell's thickness, computed
from state-dependent geometry (AD- and ALE-clean) as

$$
\Delta z_{1/2}
= \frac{1}{2}\,\frac{V_{\text{cell}}}{A_{\text{face}}}
\qquad\text{(`boundary\_temperature.jl`)}.
$$

Using the half-width (rather than the full cell width) places the surface at the
physical face and the datum $T_{\text{cell}}$ at the cell center, which is the
correct one-sided gradient for a cell-centered scheme.

### 12.2.2 Residual and Newton iteration

Define the surface residual

$$
f(T_s) \;=\; q_{\text{BC}}(T_s)\;-\;q_{\text{cond}}(T_s).
$$

The root $f(T_s)=0$ is found by Newton iteration

$$
T_s^{(k+1)} = T_s^{(k)} - \frac{f(T_s^{(k)})}{f'(T_s^{(k)})},
\qquad
f'(T_s) = \frac{dq_{\text{BC}}}{dT_s} - \frac{k_{\text{eff}}}{\Delta z_{1/2}}
$$

(`solve_boundary_temperature!`, `boundary_temperature.jl`).
The initial guess is the adjacent cell temperature clamped to the property
bracket, $T_s^{(0)} = \min\bigl(\max(T_{\text{cell}}, 200\ \text{K}),
3000\ \text{K}\bigr)$, so the solve is a pure function of its arguments. The
per-boundary `BoundaryState` receives the converged $T_s$ and fluxes for
diagnostics and is never read back; a non-finite final iterate is not stored.
A residual evaluation therefore never depends on the ones before it: a trial
state of a diverging outer Newton, with a cell temperature far outside the
bracket, fails that trial and leaves every later evaluation untouched. The
conduction derivative
$dq_{\text{cond}}/dT_s = k_{\text{eff}}/\Delta z_{1/2}$ is exact and always
present; the BC derivatives are given per-type in §12.2.4.

**Convergence.** The residual is compared against a relative tolerance scaled by
the dominant physical flux,

$$
|f(T_s)| < \max\!\Bigl(\texttt{tol}\cdot
\max\!\bigl(|q_{\text{BC}}|,\,|q_{\text{cond}}|,\,1\bigr),\ f_{\text{floor}}\Bigr),
\qquad \texttt{tol}=10^{-10},
$$

with a maximum of 50 iterations. The floor
$f_{\text{floor}} = 32\,\epsilon_{\text{mach}}\max(|T_s|,1)\,\bigl(k_{\text{eff}}/\Delta z_{1/2}+|dq_{\text{BC}}/dT_s|\bigr)$
is the flux change from a few ulps of $T_s$: the conducted flux is a difference of
two temperatures, so a face with near-zero net flux cannot resolve its residual
below it. Scaling by
the dominant flux makes the criterion meaningful whether the surface flux is
$\sim 10^4$ W·m⁻² (cone) or near zero (near-adiabatic). **Nothing aborts on
non-convergence**: on max-iterations or a vanishing Jacobian
($|f'|<10^{-30}$ — e.g. an adiabatic boundary with no temperature dependence
anywhere) the loop emits a throttled `@warn` and falls through to the final
flux evaluation at the unconverged iterate; the ODE integrator's step-error
control is the backstop that rejects a step built on a bad surface state.

**Clamping.** Each Newton iterate is clamped to the physical window
$[200,\,3000]$ K to prevent unphysical
excursions during the transient. The clamp constructs the clamped value with
**zero ForwardDiff partials** (ForwardDiff `convert` semantics) — derivatives
do *not* ride along through a clamped iterate. This self-heals in the normal
path because a subsequent un-clamped Newton update restores the implicit-
function derivatives; only when the *final* iterate is pinned at 200/3000 K
(root outside the window) are the returned $T_s$/$q$ sensitivities zero —
which is then the mathematically correct clamp derivative, and the solver
warns (throttled) with the AD consequence spelled out.

For the physical BCs, convergence is fast and robust: convection contributes a
constant negative slope $-h_{\text{conv}}$ (one Newton step is exact), and
radiation contributes a monotone cubic slope $-4F\varepsilon_{\text{eff}}\sigma
T_s^3$ that Newton handles in 2–5 iterations. Combined BCs sum monotone
contributions, so the residual is monotone decreasing in $T_s$ and the root is
unique.

### 12.2.3 Automatic differentiation through the solve

The solver is written to carry whatever scalar type its inputs carry. It promotes
to the widest type among $(T_{\text{cell}}, k_{\text{eff}}, \Delta z_{1/2},
t, \varepsilon_{\text{eff}}, a_{\text{eff}}, \dot m_g)$; under sensitivity analysis (§17, Sensitivity
Analysis) this is a `ForwardDiff.Dual`. The Newton *convergence test* uses
`ForwardDiff.value(f)` (the primal residual), so the iteration count is identical
with or without AD, while the returned tuple $(T_s,\,q_{\text{applied}},\,
q_{\text{conducted}})$ is fully Dual-typed and carries sensitivities into the
Jacobian. One guard is needed: on the very first iteration under AD the
perturbed residual may pass the primal tolerance while the dual part has not yet
propagated, so convergence is suppressed for `dual_mode && i==1`. The diagnostic
`BoundaryState` stores plain `Float64` (it is bookkeeping, never differentiated).

### 12.2.4 Per-type flux and derivative

The Newton kernel calls `evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff,
radiation_model, face)` (the `BlowingFace` at the iterate, read only by convective
BCs, §12.3.3), which returns $(q_{\text{BC}}, dq_{\text{BC}}/dT_s)$ for each
concrete thermal BC. The exact pairs (verified against
`boundary_temperature.jl`) are collected below; the
$\varepsilon_{\text{eff}}$ argument is the composition-weighted effective surface
emissivity from the property cache and the keyword `a_eff` (default
$\varepsilon_{\text{eff}}$) the effective surface absorptivity for incident
radiation (§12.2.6).

| BC | $q_{\text{BC}}(T_s)$ | $dq_{\text{BC}}/dT_s$ |
|---|---|---|
| `AdiabaticBC` | $0$ | $0$ |
| `HeatFluxBC` | $q(t)$ (const or $q(t)$) | $0$ |
| `ConvectiveBC` | $h_{\text{conv}}(t)\,(T_\infty(t)-T_s)$ | $-h_{\text{conv}}(t)$ |
| `RadiativeBC` | $F\,\varepsilon\,\sigma\,(T_\infty^4-T_s^4)$ | $-4F\,\varepsilon\,\sigma\,T_s^3$ |
| `RadiativeBC` (under `P1_QUASI_STEADY`) | $0$ | $0$ |
| `RadiativeFluxBC` (surface modes) | $\alpha\,q_{\text{ext}}(t)$ | $0$ |
| `RadiativeFluxBC` (volumetric modes) | $0$ | $0$ |
| `DirichletThermalBC` (penalty) | $h_{\text{large}}\,(T_{\text{bc}}-T_s)$ | $-h_{\text{large}}$ |
| `CombinedThermalBC` | $\sum_m q_{\text{BC},m}$ | $\sum_m dq_{\text{BC},m}/dT_s$ |

Each BC is derived and motivated in §12.3.

**Neglected emissivity slope.** For a temperature-dependent emissivity
$\varepsilon=\varepsilon(T_s)$ the exact radiative derivative would be

$$
\frac{dq_{\text{BC}}}{dT_s}
= -4F\varepsilon(T_s)\sigma T_s^3
  \;+\; F\,\frac{d\varepsilon}{dT_s}\,\sigma\,(T_\infty^4 - T_s^4),
$$

but the code keeps only the first term. The
omitted term is part of the Newton *Jacobian* of the iteration, not the residual
itself — the residual still uses the true $\varepsilon(T_s)$ — so the converged
$T_s$ is unbiased; only the convergence rate degrades from quadratic to
super-linear when $d\varepsilon/dT$ is large. This is documented as a known
limitation (§12.12).

### 12.2.5 Application into the flux array and the sign flip

After the Newton solve, `_apply_thermal_bc_static!` writes the **conducted** flux
$q_{\text{conducted}}$ (not the applied flux) into the energy face-flux array,
with the boundary-orientation sign

$$
q[f] =
\begin{cases}
+q_{\text{conducted}}, & \text{bottom (}z=0\text{)},\\[2pt]
-q_{\text{conducted}}, & \text{top (}z=L\text{)},
\end{cases}
$$

(`finite_volume.jl`). Writing the conducted flux (rather than the applied
flux) is what makes the scheme energy-conserving: the divergence operator then
sees the flux that actually enters the first interior cell, and the surface energy
balance $q_{\text{BC}}=q_{\text{cond}}$ guarantees that the energy delivered at
the surface equals the energy conducted inward. The diagnostic check
`check_boundary_energy_balance` returns the relative imbalance
$|q_{\text{applied}}-q_{\text{conducted}}|/|q_{\text{applied}}|$ (the absolute
difference when $|q_{\text{applied}}|\le10^{-10}$)
per boundary, which is $\approx$ machine precision when Newton has converged
(`boundary_temperature.jl`).

Recall (§7, §13) that the interior conductive fluxes are explicitly **zeroed at
boundary faces** in `compute_conductive_fluxes!` before BC application, so the
overwrite here is the sole source of the boundary face flux; there is no
double-counting.

### 12.2.6 Effective surface emissivity and absorptivity (composition-weighted)

When a `RadiativeBC` is constructed with $\varepsilon=\texttt{nothing}$, the
re-radiation uses the boundary cell's effective emissivity
$\varepsilon_{\text{eff}}$, computed in the property cache as the volume-weighted
mean over the solid/liquid components with each component emissivity evaluated
at the cell temperature,

$$
\varepsilon_{\text{eff}} = \sum_j v_j\,\varepsilon_j(T_{\text{cell}}),
$$

clamped to $[0,1]$ (see §5.8.1). When a `RadiativeFluxBC` is constructed with
$\alpha=\texttt{nothing}$, the absorbed incident flux uses the boundary cell's
effective surface absorptivity

$$
a_{\text{eff}} = \sum_j v_j\,a_j,\qquad
a_j = \texttt{surface\_absorptivity}_j \text{ if set, else } \varepsilon_j(T_{\text{cell}}),
$$

(§5.8.2). Two physical principles are encoded here. First, **Kirchhoff's law**
$\alpha=\varepsilon$ at thermal equilibrium is the default: with no component
absorptivity set, $a_{\text{eff}}\equiv\varepsilon_{\text{eff}}$ and the same
quantity serves as both emissivity (re-radiation) and absorptivity (absorbed
incident flux). A component-level `absorptivity` relaxes this where the heater
spectrum differs from the surface's own emission spectrum, keeping the two
mixtures on the same volume weights. Second, surface emissivity is a
*geometric* property of the exposed mixture, so a volume-weighted (not
mass-weighted) average is appropriate and it **automatically tracks composition
changes** during pyrolysis (e.g. polymer $\to$ char as $v_j$ shifts). The
boundary solver reads these as `cache.ε[cell_id]` and `cache.a[cell_id]` and
threads them into `solve_boundary_temperature!` as the positional
$\varepsilon_{\text{eff}}$ and the keyword `a_eff` (`boundary_temperature.jl`);
the Jacobian operators recompute both on the Dual-seeded $(T,\xi)$, so
$\partial\varepsilon_{\text{eff}}/\partial T_{\text{cell}}$ and
$\partial a_{\text{eff}}/\partial T_{\text{cell}}$ from a temperature-dependent
component emissivity enter the surface energy row exactly.

### 12.2.7 `BoundaryState` storage

The per-boundary diagnostic record (`src/Geometry/mesh.jl`) holds the
solved surface temperature, the applied flux, the conducted flux, the inputs
$(T_{\text{cell}}, k_{\text{eff}}, \Delta z_{1/2})$, convergence metadata
(`converged`, `iterations`), and the blowing diagnostics at the converged $T_s$
(`blowing_mcp`, `blowing_ratio`, and the pooling-model fields `pooling_factor`,
`face_fraction`, `face_buoyancy`; §12.3.3). The fields are plain `Float64`
(`converged` a `Bool`, `iterations` an `Int`);
they are diagnostics and are never differentiated. The live (possibly Dual-typed)
quantities are returned from the solve and flow into the residual.

---

## 12.3 Thermal boundary conditions

This section derives each concrete thermal BC, its physical meaning, and how it
plugs into §12.2. All flux laws are signed **positive into the domain**.

### 12.3.1 `AdiabaticBC`

Perfect insulation: $q_{\text{BC}}=0$, $dq_{\text{BC}}/dT_s=0$. The surface energy
balance reduces to $q_{\text{cond}}=0$, hence $T_s=T_{\text{cell}}$ — the surface
floats to the adjacent cell temperature with zero conduction across the face.
Typical for the substrate face in a cone-calorimeter idealization with a
well-insulated back. (`boundary_fluxes.jl`,
`boundary_temperature.jl`.)

### 12.3.2 `HeatFluxBC`

A prescribed flux $q_{\text{BC}}=q(t)$, constant (`HeatFluxBC{<:Real}`) or a
function of time (`HeatFluxBC{<:Function}`), independent of $T_s$ so
$dq_{\text{BC}}/dT_s=0$. The Newton solve still runs to find the $T_s$ that
balances this against conduction. Sign: positive $q$ heats the material.
(`boundary_fluxes.jl`, `boundary_temperature.jl`.)

> **$T_s$ is a half-cell extrapolation.** For a pure prescribed-flux boundary
> the surface Newton reduces to the linear extrapolation
> $T_s = T_{\text{cell}} + q''\,\Delta z_{1/2}/k_{\text{eff}}$, which assumes
> the *entire* prescribed flux crosses the boundary half-cell by conduction.
> When a volumetric consumer sits inside that half-cell (a reaction zone under
> MW-scale flux, in-depth absorption), the reported $T_s$ is biased by up to
> tens of K at coarse resolution — a real but $n$-shrinking effect. It enters
> the Landau traveling-wave identity of V-CPL-02 through
> $v_{ss}(\text{measured}\ T_s)$: with the reaction zone under-resolved
> ($\Delta z \approx 2\delta_r$) the identity misses by ≈1.3 %, collapsing to
> ≈5×10⁻⁴ once $\Delta z \lesssim \delta_r$. (A *mesh-independent* residual
> of up to about 1 % has a different source, the wave-development transient of
> a uniform-$T_0$ start; the case measures after at least six development time
> constants to exclude it.) Interpret surface-temperature outputs under `HeatFluxBC`
> accordingly, and refine the near-surface mesh when $T_s$ itself is a target
> quantity.

### 12.3.3 `ConvectiveBC` (Newton's law of cooling)

$$
q_{\text{BC}} = h_{\text{conv}}(t)\,\bigl(T_\infty(t) - T_s\bigr),
\qquad
\frac{dq_{\text{BC}}}{dT_s} = -h_{\text{conv}}(t).
$$

Both $h_{\text{conv}}$ and $T_\infty$ may be constants or functions of $t$
(evaluated through `evaluate_param`). The flux is positive (heating) when
$T_\infty>T_s$ and negative (cooling) otherwise, so for a hot surface in a cool
ambient this is a convective loss. Because the derivative is the constant
$-h_{\text{conv}}$, the contribution to the Newton residual is linear and exactly
resolved in a single step. (`boundary_temperature.jl`.)

**Blowing.** Every `ConvectiveBC` carries a blowing model (`bc.blowing`, an
`AbstractBlowingModel`, default `NoBlowing()`), and the surface balance applies
it at every Newton iterate to the face state: the **signed outward** gas mass
flux of every component, read from the same `J_gas[:, f]` the mass-flux pass
populates (film-law diffusion plus, in Darcy–Fick mode, the pressure-BC venting
flux), the surface-temperature iterate and the material (`BlowingFace`,
`Discretization/blowing.jl`). The convective flux is $h\,\chi\,(T_\infty - T_s)$
with $\chi = h_{\text{eff}}/h$ from the model (`blowing_state`). It modifies
only the convective *coefficient*, not the gas enthalpy accounting, so it
composes exactly with the $S_{\text{conv}}$-only convention (§7.4.3) — no
double-count is possible. Zero outflow gives $\chi = 1$ exactly for every
shipped model, so the back face is unaffected automatically.

*Film law with a fitted constant* (`CouetteFilm(λ)`):

$$
\chi \;=\; \frac{\phi}{e^{\phi}-1},
\qquad
\phi \;=\; \frac{2\lambda\,\dot m''\,c_{p,g}}{h},
\qquad
\dot m''\,c_{p,g} \;=\; \sum_{g}\dot m_g\,c_{p,g}(T_s).
$$

$\lambda = 1/2$ is the classical Couette-film (Spalding) factor, $\phi = B =
\dot m'' c_{p,g}/h$; a smaller $\lambda$ fitted to a resolved flow weakens the
reduction at a given outflow and holds for that gas only. Limits in $B$:
$B = 0$ recovers the bare $h$ exactly; $B = 4$ (peak-PMMA outgassing against
canonical cone $h$) gives $h_{\text{eff}}/h = 0.075$; $B < 0$ (net inflow) is
the smooth suction enhancement.

*Property-ratio film for any light gas* (`PropertyRatioFilm(free_stream; n_Γ)`),
anchored on the Couette film of the free stream:

$$
\phi \;=\; \frac{\dot m''\,c_{p,\infty}(T_s)}{h}\left(\frac{\Gamma_\infty}{\Gamma_g}\right)^{n_\Gamma},
\qquad \Gamma = \frac{k}{c_p},
$$

with $\Gamma_\infty$ of the free stream and $\Gamma_g$ of the mixture leaving the
face (heat capacity weighted by mass fraction, conductivity by Wilke's rule
over the components with outward flux, each species' viscosity from Eucken's
relation when the component carries none). $n_\Gamma = 0.6$ is Kays's turbulent
exponent on the heat-capacity ratio, reproduced by the chamber CFD; with
$n_\Gamma = 0$ this is the Couette film in the free-stream heat capacity. Under
net inflow the arriving gas is the free stream and the ratio is 1.

*Pooling* (`PoolingFilm(free_stream; n_Γ, n_D, ψ0, Λc, Λw)`): the same
$\phi$ with a pooling factor $\psi$, the share of the no-outflow exchange the
face keeps,

$$
\chi \;=\; \frac{\phi}{e^{\phi/\psi}-1} \;=\; \psi\,f(\phi/\psi),
\qquad
\psi \;=\; \psi_0 + \frac{1-\psi_0}{1+\exp[-(\Lambda-\Lambda_c)/\Lambda_w]},
\qquad
\Lambda \;=\; \frac{1-(M_{\text{face}}/T_s)(T_\infty/M_\infty)}{1-T_\infty/T_s},
$$

where $f$ is the film factor and $\Lambda$ the face buoyancy relative to the
hot free stream, from the face molar mass $M_{\text{face}}$ of the species
mass fractions

$$
Y_i \;=\; \omega_i\left[1-\exp\!\left(-\frac{\dot m''}{\psi\,g\,s_i^{\,n_D}}\right)\right],
\qquad g = \frac{h}{\bar c_{p,\infty}},
\qquad s_i = \frac{D_{i,\infty}}{\alpha_\infty},
$$

with $\omega_i$ the species' share of the outflow, $g$ the gas-exchange
conductance without outflow (the film analogy: $h$ over the free-stream heat
capacity averaged between $T_\infty$ and $T_s$, a 5-point Gauss–Legendre mean),
$D_{i,\infty}$ the binary diffusivity into the free stream (Fuller, from the
component's diffusion volume on the model, else the Chapman–Enskog estimate
from its molar mass; a lumped `volatiles` gas states its own) and $n_D = 1/2$
the penetration-theory exponent. Components whose name matches the free
stream count in $\dot m''$ but not in the face composition. $\psi$ and the
$Y_i$ are solved together: $G(\psi) = \Psi(\Lambda(Y(\psi))) - \psi$ has
$G(\psi_0) > 0$ and $G(1) < 0$, so a 16-interval scan of $[\psi_0, 1]$ on the
primal values brackets every root, each refined by safeguarded Newton; where
two states are self-consistent the mean of the lowest and highest is taken
(the resolved flow alternates between them). One Newton step in the working
scalar type at the primal root then gives the implicit-function derivative,
so `ForwardDiff` partials ride through the solve. A face no hotter than the
free stream keeps $\psi = 1$. The three switch constants belong to the
apparatus flow, as $h$ does.

*User rule* (`FunctionBlowing(f)`): $\chi = f(\text{face}, h, T_\infty)$.

The factor is evaluated by the numerically stable `blowing_factor` (Taylor
branch at $|\phi|<10^{-5}$, overflow-safe $\phi\,e^{-\phi}/(1-e^{-\phi})$ form for
large $\phi$). The Newton derivative uses the frozen-$\chi$ form $-h\,\chi$ —
first-order consistent, so convergence is marginally sub-quadratic and the
converged root is exact. Each convective member of a `CombinedThermalBC` applies its own model
against its own $h$ and $T_\infty$; all other BC types are untouched. At the
converged $T_s$ the `BoundaryState` records the blowing conductance
`blowing_mcp`, the applied ratio `blowing_ratio` and, for the pooling model,
`pooling_factor`, `face_fraction` and `face_buoyancy` (the first convective
member's), and the solution extras carry the ratio and the pooling fields per
output time (`surface_blowing_ratio`, `surface_pooling_factor`,
`surface_face_fraction`, `surface_face_buoyancy`). The structured
Jacobian's `BCThermalOperator` and `BCConvectionOperator` rebuild $\dot m''$
from the stencil state via the shared boundary-face kernel, so the
$\partial(h\chi)/\partial(T_b,\xi_b)$ chain is exact (support-complete strict
verification passes with the film and pooling models in Fickian, Darcy–Fick
and ALE modes). Giving a convective BC a model on a case calibrated with the
bare coefficient will raise $T_s$ and MLR; introduce it together with flame
heat-flux feedback before any recalibration pass (the two partially cancel;
§1.8.1).
(`Materials/blowing_models.jl` types and gas data,
`Discretization/blowing.jl` evaluation.)

### 12.3.4 `RadiativeBC` (gray surface exchange with surroundings)

$$
q_{\text{BC}} = F\,\varepsilon\,\sigma\,\bigl(T_\infty^4 - T_s^4\bigr),
\qquad
\frac{dq_{\text{BC}}}{dT_s} = -4\,F\,\varepsilon\,\sigma\,T_s^3,
$$

with the Stefan–Boltzmann constant $\sigma=5.670374419\times10^{-8}$ W·m⁻²·K⁻⁴
(`Materials/constants.jl`, exported as `STEFAN_BOLTZMANN`). The fields are: the
**view factor** $F\in[0,1]$ (fraction of the surface's hemisphere that exchanges
with surroundings at $T_\infty$; default $1$), the **emissivity** $\varepsilon$
(either a supplied constant/function or, when `nothing`, the composition-weighted
$\varepsilon_{\text{eff}}$ of §12.2.6), and the surroundings temperature
$T_\infty$. This term is the surface's own re-radiation to and absorption from a
gray environment; it is the dominant cooling mechanism for a hot pyrolyzing
surface. The cubic temperature derivative is the source of the $\sim$20% error
that would result from naively evaluating radiation at the cell-center
temperature — and the reason the Newton solve exists. (`boundary_temperature.jl`.)

The neglected $d\varepsilon/dT$ term (§12.2.4) applies here when $\varepsilon$ is
temperature-dependent.

### 12.3.5 `RadiativeFluxBC` (external incident radiation; model-dependent)

This BC represents incident radiation from an *external* source (cone heater,
furnace, radiant panel) — not the surface's own emission. The absorbed fraction is
$\alpha\,q_{\text{ext}}$, where $q_{\text{ext}}=$ `flux` (constant or $q(t)$) and
the absorptivity $\alpha$ is either supplied or, when `nothing`, defaults to
the mixture surface absorptivity $a_{\text{eff}}$, which equals $\varepsilon_{\text{eff}}$
by Kirchhoff's law unless a component sets its own `absorptivity` (§12.2.6).

Its behavior is **dispatched on the global `radiation_model`** passed to
`solve()` (§8, Thermal Radiation). Note the defaults diverge: `solve()`
defaults to `NO_RADIATION`, while several internal helpers
(`apply_thermal_bc_flux!`, `update_boundary_state!`, `solve_boundary_temperature!`
and the solution extractors) default to
`SURFACE_ABSORPTION` — a trap for direct callers of the internals; always
pass `radiation_model` explicitly when bypassing `solve()`. The crucial
design point is that incident radiation must be deposited exactly once:

- **`NO_RADIATION` and `SURFACE_ABSORPTION`**: the radiation is
  treated as a *surface* heat flux. The Newton solve receives
  $$q_{\text{BC}} = \alpha\,q_{\text{ext}}(t), \qquad dq_{\text{BC}}/dT_s = 0,$$
  so the absorbed flux enters the surface energy balance and correctly raises
  $T_s$ above $T_{\text{cell}}$ when there is net heating
  (`boundary_temperature.jl`).
- **`BEER_LAMBERT` and `P1_QUASI_STEADY`**: the radiation is deposited
  *volumetrically* by `compute_radiation_source!` (Beer–Lambert) or by the P1
  residual path. To avoid double-counting, the boundary flux is forced to zero,
  $$q_{\text{BC}} = 0, \qquad dq_{\text{BC}}/dT_s = 0$$
  (`boundary_temperature.jl`, `boundary_fluxes.jl`).

A note on consistency. The *Newton-solve* evaluator
(`evaluate_bc_flux_and_derivative`) applies $\alpha q_{\text{ext}}$ as a boundary
flux for `SURFACE_ABSORPTION`, while `compute_radiation_source!` deposits *zero*
volumetric radiation for `SURFACE_ABSORPTION`. These are consistent:
`SURFACE_ABSORPTION` absorbs at the boundary face only, so the energy enters via
the boundary flux and the volumetric source is zero. Only `BEER_LAMBERT` (and the
experimental P1 path) routes energy through the volume, and only then is the
boundary flux suppressed. The diagnostic decomposition
`compute_boundary_heat_fluxes` follows the *same* dispatch: its `net` is the
flux the Newton balance actually applies — $\alpha q_{\text{ext}}$ is inside
`net` under `NO_RADIATION`/`SURFACE_ABSORPTION` and absent under
`BEER_LAMBERT`/P1, where the volumetric deposit is reported separately as
`extras.volumetric_absorption` (§12.9.1).

> **Implementation note.** The incident intensity is extracted at compile time
> from the (possibly combined) BC by `extract_radiative_flux_from_bcs(bc_set,
> Val(tag), t)`, which returns `(I, has_explicit_absorptivity)` — the flux with
> an explicit absorptivity already folded in, or the bare incident flux whose
> $a_{\text{eff}}$ scaling happens at injection
> (`boundary_conditions.jl`, §8.3.4). The intensity is passed as `I_radiation`
> to `compute_rhs!` together with the `radiation_kirchhoff` flag.

### 12.3.6 `DirichletThermalBC` (prescribed surface temperature)

A prescribed surface temperature $T=T_{\text{bc}}(t)$ is enforced by one of two
mechanisms, selected by how the BC appears, rather than a hard algebraic
constraint:

1. **In the Newton solve** (`boundary_temperature.jl`): a penalty "flux"
   with a very large conductance,
   $$q_{\text{BC}} = h_{\text{large}}\,(T_{\text{bc}} - T_s),\qquad
   dq_{\text{BC}}/dT_s = -h_{\text{large}}, \qquad h_{\text{large}}=10^{8}
   \;\text{W·m}^{-2}\text{·K}^{-1}.$$
   The large slope drives $T_s\to T_{\text{bc}}$ to high accuracy in one step. The
   value $10^8$ is large enough to pin $T_s$ to $T_{\text{bc}}$ within Newton
   tolerance while keeping the iteration well conditioned (a still-larger penalty
   would degrade conditioning with no accuracy benefit).
2. **Standalone, in the RHS assembly** (`finite_volume.jl`): a bare `DirichletThermalBC`
   is special-cased and bypasses the Newton path entirely. The surface temperature is the
   prescribed value, $T_s = T_{\text{bc}}(t)$, and the conducted flux is written directly,
   $$q_{\text{bc}} = k_{\text{cell}}\,\frac{T_{\text{bc}} - T_{\text{cell}}}{\Delta z_{1/2}},
   \qquad
   q[f]= \begin{cases}+q_{\text{bc}},&\text{bottom},\\ -q_{\text{bc}},&\text{top},\end{cases}$$
   with $\Delta z_{1/2}=\tfrac12\,V_{\text{cell}}/A_{\text{face}}$ — the boundary face's own
   area, identical to the Newton path's distance, so the two mechanisms agree under a
   variable cross-section. The bypass also updates the diagnostic `BoundaryState`
   ($T_{\text{surface}}=T_{\text{bc}}$, $q_{\text{applied}}=q_{\text{conducted}}=q_{\text{bc}}$,
   `converged = true`), so surface-temperature and flux diagnostics report the prescribed
   state. The penalty mechanism of item 1 is used only when the Dirichlet BC is *nested
   inside* a `CombinedThermalBC` (the composite goes through the Newton path; the two
   mechanisms agree to $\sim10^{-4}$ K).

The penalty approach is acknowledged as less clean than a direct constraint (a
genuine algebraic constraint would require a DAE index reduction); it is robust in
practice because the prescribed temperature pins the boundary tightly.

### 12.3.7 `CombinedThermalBC` (additive superposition)

Real surfaces experience several mechanisms at once. The `+` operator builds an
additive composite, e.g. the canonical cone-calorimeter top BC

```julia
bc = RadiativeFluxBC(flux = 50e3)            # absorbed incident cone flux
   + RadiativeBC(ε = nothing, T_inf = 300)   # surface re-radiation loss
   + ConvectiveBC(h = 10.0, T_inf = 300)     # convective loss
```

Fluxes and derivatives sum termwise,

$$
q_{\text{BC}} = \sum_m q_{\text{BC},m}(T_s,t),\qquad
\frac{dq_{\text{BC}}}{dT_s} = \sum_m \frac{dq_{\text{BC},m}}{dT_s}.
$$

The sum is built by a `@generated` function that unrolls the tuple at compile time
for type stability (`boundary_temperature.jl`). Because the radiation model can affect a
`RadiativeFluxBC` member, the combined evaluator threads `radiation_model` through
to every child. The `+` operator overloads handle every nesting case (BC+BC,
combined+BC, BC+combined, combined+combined) so associativity holds
(`boundary_conditions.jl`).

---

## 12.4 Mass (species) boundary conditions

The species boundaries supply the gas mass flux $J_j$ at each boundary face.
Unlike the thermal boundary, there is no surface-temperature nonlinearity, so the
flux is evaluated directly from the adjacent cell state. The internal diffusive
flux model is the volume-fraction (partial-pressure) gradient form of §9,
Gas-Phase Mass Transport & Pressure; here it is closed at the surface by an
external resistance.

### 12.4.1 `ImpermeableBC` (default)

$J_j = 0$ for all components — no gas crosses the boundary. This is the default
for both faces (`boundary_fluxes.jl`) and is the physically correct choice for
the substrate (sealed back face) and any non-venting surface. It is also the
implicit closure for solid/liquid components everywhere, since only gas components
are transported.

### 12.4.2 `MassFluxBC`

A prescribed species flux $J_j = J(t)$ (constant or function), positive into the
domain (`boundary_fluxes.jl`). Used to impose a known gas influx/efflux,
typically zero or negative (gas escaping).

### 12.4.3 `ConvectiveMassBC` (external film transfer)

$$
J_j = h_m(t)\,\bigl(\xi_{j,\infty} - \xi_{j,s}\bigr),
$$

the mass-transfer analog of Newton's cooling law, with the film coefficient $h_m$
(m·s⁻¹, velocity-based) and the ambient concentration $\xi_{j,\infty}$ (typically
$0$ for pyrolysis gas escaping to atmosphere) (`boundary_fluxes.jl`). The
coefficient may be estimated by the heat–mass analogy,

$$
h_m \approx \frac{h_{\text{conv}}}{\rho_{\text{air}}\,c_{p,\text{air}}},
$$

e.g. $h_{\text{conv}}\approx10$ W·m⁻²·K⁻¹ gives $h_m\approx0.008$ m·s⁻¹. With
$\xi_{j,\infty}=0$ and $\xi_{j,s}>0$, $J_j<0$ (escaping), consistent with the
flux-into-domain sign.

The ambient $\xi_{j,\infty}$ accepts three forms, resolved per gas species $j$ by
`ambient_concentration(bc, t, j)`: a `Real` (the same ambient for every gas species —
typically 0, the atmosphere carrying no pyrolyzate), a `Function` of $t$ (time-varying
scalar, broadcast to every gas species), or an `NTuple{NC}` of per-component ambient
concentrations indexed by the material's component order (condensed-phase entries are
ignored; entries are validated non-negative at construction and the tuple length against
$N_C$ at problem validation). The per-species form is what represents ambient air ingress
without injecting pyrolyzate: a scalar $\xi_\infty>0$ would feed that concentration of
*every* gas species into the domain.

### 12.4.4 `DiffusiveMassBC` (combined internal/external resistance)

This BC couples the **internal** Fickian diffusion resistance from the surface
cell to the face with an **external** film resistance, in series. Working in the
volume-fraction formulation $J=-\lambda\,\partial\xi/\partial x$ with the surface
gradient approximated as $(\xi_s-0)/d$ (zero exterior concentration), the two
resistances are

$$
R_{\text{internal}} = \frac{d}{\lambda},
\qquad
R_{\text{external}} = \frac{1}{h_m},
$$

where $d=\Delta z_{1/2}$ is the cell-center-to-face distance and $\lambda$ is the
effective gas-transfer coefficient. The combined flux is

$$
J_j = -\frac{\xi_{j,s}}{R_{\text{internal}}+R_{\text{external}}}
    = -\frac{\lambda\,h_m\,\xi_{j,s}}{\lambda + h_m\,d}.
$$

The minus sign places escaping gas at $J_j<0$ (flux-into-domain convention). Two
limits follow:

- **No external resistance** ($h_m\to\infty$, the default `DiffusiveMassBC()` with
  `h_m = Inf`): pure internal diffusion,
  $$J_j = -\frac{\lambda\,\xi_{j,s}}{d};$$
- **Impermeable** ($h_m\to0$): $J_j\to0$.

> **Implementation note.** The flux law is in `mass_bc_flux(bc::DiffusiveMassBC,
> ξ_surface, λρ, ρ, d, t, j)` (`boundary_fluxes.jl`) — every mass-BC flux method carries
> the trailing species index `j` so per-species ambients resolve correctly. The
> discretization passes the *combined transport coefficient* $\lambda\rho$ (consistent
> with the rigorous gas-transport closure) and the density $\rho$; for $\rho \le 0$ the
> flux is zero, and otherwise the function recovers $\lambda = \lambda\rho/\rho$ by exact
> division (no floor — the callers construct $\lambda\rho$ as the product, so division is
> exact for any positive density), because $h_m$ is defined per concentration, not per
> mass fraction. The surface-cell quantities $(\xi_{j,s},\,\lambda,\,\rho,\,d)$ are read
> in `_apply_mass_bc_static!` (`finite_volume.jl`).

**Sign reconciliation.** Three sign conventions coexist and must be reconciled
(`boundary_fluxes.jl`):

| Quantity | Convention |
|---|---|
| `mass_bc_flux(...)` | **negative** = gas escaping (flux-into-domain) |
| `compute_surface_mass_flux(...)` (diagnostic) | **positive** = gas escaping |
| `J[face]` stored in the flux array | **positive** = flux in $+z$ |

The discretization bridges these by orienting per boundary
(`_apply_mass_bc_static!`, `finite_volume.jl`): at the bottom
$J[j,f]\mathrel{+}=J_{\text{bc}}$, at the top $J[j,f]\mathrel{-}=J_{\text{bc}}$.
With $J_{\text{bc}}<0$ for escaping gas, the top-face stored flux is positive —
gas leaving in $+z$ at the exposed surface — exactly what the divergence
operator needs. The application **accumulates** rather than overwrites: in
Darcy–Fick mode the boundary face already carries the pressure-BC advective
venting flux (§9.5.3), and the film-law flux adds to it; in Fickian mode the
boundary faces arrive zero so the accumulation is equivalent to assignment.

### 12.4.5 `CombinedMassBC`

Like the thermal composite, $J_j=\sum_m J_{j,m}$, built by a `@generated` unrolled
sum (`boundary_fluxes.jl`), e.g. `DiffusiveMassBC(Inf) +
ConvectiveMassBC(h_m=0.01, ξ_inf=0.0)`. Pressure BCs (§12.5) participate in the
composite by contributing zero *diffusive* flux (they act on the Darcy velocity,
not the diffusive flux), so they can be combined freely with `DiffusiveMassBC`
(`boundary_fluxes.jl`).

---

## 12.5 Pressure boundary conditions (Darcy flow)

In Darcy–Fick transport mode (§9), gas advection is driven by the Darcy velocity
$v_g = -(\kappa/\mu)\nabla P$, and the pressure boundaries set $v_g$ at the
boundary faces. These BCs are subtypes of `AbstractMassBC` so they live in the
mass slot of the `BoundaryConditionSet`; they do **not** contribute a diffusive
flux (returning $0$ from `mass_bc_flux`, §12.4.5) — they act through
`apply_pressure_bc_to_darcy_velocity!` (`darcy_velocity.jl`), which runs *after*
the interior Darcy velocities are computed and sets the boundary face velocity
$v_b$. That velocity then enters the species equations as the boundary advective
venting flux $J_j=v_b\,\xi_j^{\text{int}}/\varphi_{\text{int}}$ (outflow
only; §9.5.3), written by `compute_gas_fluxes_darcy_fick!` before the diffusive
mass BC is added.

### 12.5.1 `DirichletPressureBC` (prescribed pressure)

Prescribes the boundary pressure $P_{\text{bc}}$ (default $101325$ Pa,
atmospheric). The boundary face velocity follows from Darcy's law with a one-sided
pressure gradient between the boundary and the interior cell:

$$
v_g = -\frac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\nabla P,
\qquad
\nabla P =
\begin{cases}
(P_{\text{bc}} - P_{\text{int}})/d, & \text{top},\\[2pt]
(P_{\text{int}} - P_{\text{bc}})/d, & \text{bottom},
\end{cases}
$$

with $d$ the cell-center-to-face distance, $\mu_{\text{int}}$ the material's
pore-gas viscosity at the cell temperature (§5.9), and the velocity clamped to zero if
$\kappa_{\text{int}}\le0$ or $\mu_{\text{int}}<10^{-20}$ Pa·s
(`darcy_velocity.jl`). The orientation of the gradient ensures the correct
outflow direction at each face: at the top a higher interior pressure
($P_{\text{int}}>P_{\text{bc}}$) yields $v_g>0$ (outflow in $+z$); at the bottom
the same drives $v_g<0$ (outflow in $-z$).

### 12.5.2 `NeumannPressureBC` (prescribed gradient)

Prescribes the pressure gradient $\partial P/\partial n$ along the **outward
normal** (default $0$ = impermeable). For zero gradient the boundary velocity is
set to zero; for a nonzero prescribed gradient,
$v_g=-(\kappa_{\text{int}}/\mu_{\text{int}})\,\partial P/\partial n$ at the
top and $v_g=+(\kappa_{\text{int}}/\mu_{\text{int}})\,\partial P/\partial n$
at the bottom (where $\hat n=-\hat z$, so $\partial P/\partial z=-\partial
P/\partial n$), with the same $\kappa,\mu$ guards (`boundary_darcy_velocity`,
`darcy_velocity.jl`). The zero-gradient default makes a `NeumannPressureBC()`
equivalent to an impermeable wall for Darcy flow.

### 12.5.3 `ConvectivePressureBC` (Robin / surface resistance)

A Robin-type law for a surface with finite resistance to gas escape (char layer,
protective coating, blowing resistance):

$$
v_g\cdot\hat n = h_P\,\bigl(P_{\text{int}} - P_{\text{ext}}\bigr),
$$

with $h_P$ the pressure-transfer coefficient (m·Pa⁻¹·s⁻¹) and $P_{\text{ext}}$ the
exterior pressure, gated on $\kappa_{\text{int}}>0$ — an impermeable surface
cell cannot vent regardless of $h_P$. The implementation orients the outward
normal per boundary (`boundary_darcy_velocity`, `darcy_velocity.jl`):

$$
v_g[f] =
\begin{cases}
+h_P\,(P_{\text{int}}-P_{\text{ext}}), & \text{top},\\[2pt]
-h_P\,(P_{\text{int}}-P_{\text{ext}}), & \text{bottom}.
\end{cases}
$$

The limits bracket the other two pressure BCs: $h_P\to\infty$ pins the boundary
pressure (Dirichlet) and $h_P\to0$ seals it (zero-gradient Neumann / impermeable).

> **Implementation note.** All three pressure BCs are scalar methods of the
> shared kernel `boundary_darcy_velocity(bc, P_int, κ_int, μ_int, d, Val(tag))`
> (`darcy_velocity.jl`); a `@generated` overload sums the contributions of a
> `CombinedMassBC`, and the `AbstractMassBC` fallback contributes $0$ (e.g. a
> `DiffusiveMassBC` in the mass slot does not set a Darcy velocity). The
> face loop `_apply_pressure_bc_faces!` writes the kernel value per boundary
> face; the orchestration entry `apply_pressure_bc_to_darcy_velocity!` applies
> top then bottom. The same kernel is evaluated (with `Dual` arguments) inside
> the `BCMassOperator` Jacobian and (with the saved state) in the MLR
> diagnostics, so the three paths cannot drift apart.

---

## 12.6 Mesh-motion boundary conditions (ALE)

In ALE mode (§11, The ALE Framework) the node positions $z_i$ are state variables
and the mesh-velocity boundary treatment is **fixed by the formulation, not
user-configurable**: the bottom node is pinned ($w_1 = 0$, hard-coded in the
residual's mesh-velocity integration) and the velocity at every other node
follows cumulatively from the species dilation $\theta_i$ and the column-global
radial strain $\theta_A$, $w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i$
(§11.3, §15.3). The surface node is therefore automatically material-following
(Lagrangian) — it recedes or expands with the integrated dilation.

There are no user-facing mesh-motion BC types: `BoundaryConditionSet` carries
only thermal and mass slots, and the ALE bottom anchor ($w_1 = 0$) is
hard-coded in the residual's mesh-velocity assembly.
See §11 for the full ALE mesh-velocity derivation and the geometric
conservation law.

---

## 12.7 The `BoundaryConditionSet` container

### 12.7.1 Structure and tag keying

All boundaries are bundled into a single type-stable container
(`boundary_conditions.jl`):

```julia
struct BoundaryConditionSet{TT<:NamedTuple, MT<:NamedTuple}
    thermal::TT   # NamedTuple{tags, Tuple{<:AbstractThermalBC,...}}
    mass::MT      # NamedTuple{tags, Tuple{<:AbstractMassBC,...}}
end
```

Each boundary tag maps to exactly one thermal BC and one mass BC. The 1D
convention uses `:top` ($z=L$) and `:bottom` ($z=0$); the field shape is designed
so that a future 2D/3D constructor can add tags without changing the runtime
contract. Storing the BCs in a typed `NamedTuple` (rather than, say, a `Dict`)
means a symbol lookup compiles to a direct field read and the residual specializes
statically.

### 12.7.2 Constructor, defaults, and validation

The keyword constructor (`boundary_conditions.jl`) merges the
user-supplied partial NamedTuples with the 1D defaults

$$
\texttt{thermal} = (\,\texttt{top}=\texttt{AdiabaticBC()},\ \texttt{bottom}=\texttt{AdiabaticBC()}\,),
\qquad
\texttt{mass} = (\,\texttt{top}=\texttt{ImpermeableBC()},\ \texttt{bottom}=\texttt{ImpermeableBC()}\,),
$$

so a caller may specify only the slots it cares about. After merging it
**validates** that every thermal payload is an `AbstractThermalBC` and every mass
payload is an `AbstractMassBC`, throwing an `ArgumentError` if a BC was placed in
the wrong slot. The default — adiabatic and
impermeable on both faces — is the inert closure: a column with no defined BCs
neither gains energy nor loses mass.

### 12.7.3 Compile-time accessors

Access is through `Val{tag}` dispatch:

```julia
@inline get_thermal_bc(bc_set, ::Val{tag}) where {tag} = bc_set.thermal[tag]
@inline get_mass_bc(bc_set,    ::Val{tag}) where {tag} = bc_set.mass[tag]
```

At a call site such as `get_thermal_bc(bc_set, Val(:top))`, the `Val(:top)`
resolves at compile time to a direct NamedTuple field read with zero overhead
(`boundary_conditions.jl`). The RHS assembly uses these accessors so that
each boundary's BC type is known statically and the per-type flux law is inlined.
The tag$\to$id bridge `tag_to_boundary_id` (§12.1.3) connects the user-facing
symbol API to the integer-keyed `Mesh1D.boundary_faces` / `boundary_states`
dictionaries.

### 12.7.4 Constructor parameter validation

Individual BC constructors optionally validate their physical parameters at build
time (`validate_bc_params`, `boundary_conditions.jl`): $h\ge0$, $T>0$ K,
$\varepsilon\in[0,1]$, $F\in[0,1]$, $\alpha\in[0,1]$, $h_m\ge0$, $P>0$.
Function-valued parameters cannot be checked at construction and are skipped
(`_get_scalar_value` returns `nothing` for functions and property functions).

---

## 12.8 Initial conditions

### 12.8.1 The state seed

The initial state is specified by a temperature closure $T_{\text{initial}}(z)$
and per-component concentration closures $\xi_{j,\text{initial}}(z)$, evaluated at
each cell center $z_i^c$ (`initial_conditions.jl`):

$$
T_i(0) = T_{\text{initial}}(z_i^c), \qquad
\xi_{j,i}(0) = \xi_{j,\text{initial}}(z_i^c), \quad j=1,\dots,N_C .
$$

The closures are evaluated and stored in the mutable mesh's cell states, then
packed into the flat state vector (§15, Temporal Integration & Structured
Jacobian). Cell volumes are not altered by initialization.

The single problem-level entry point is `apply_initial_conditions!(problem)`,
which performs the cell-state write above **and** populates the dry-solid
denominator of §12.8.3 (the populator reads the freshly written cell states,
so the two steps are ordered). `solve()` and the exported composition helper
`initialize_state(problem; …) -> (def, ws, u0)` both call it before packing
`u0`; the mesh-level utility `apply_initial_conditions!(mesh, material; T, ξ)`
writes only uniform cell states and deliberately does **not** touch
`m_dry_solid_initial`. Composing the pure-residual layer by hand without the
problem-level IC step packs the mesh constructor's placeholder state instead
of the problem's initial conditions (§15.2) — a resolution-independent error
with a perfectly plausible-looking solution.

### 12.8.2 Auto-defaults

If the user supplies no initial conditions, the `PyrolysisProblem` constructor
synthesizes them (`pyrolysis_problem.jl`):

- $T_{\text{initial}}(z)=300$ K uniformly;
- $\xi_{1,\text{initial}}(z)=\rho_1(300\,\text{K})$ if the first component is a
  solid (the virgin solid at its reference density), and $\xi_{j,\text{initial}}=0$
  for all other components.

This "first solid at reference density, everything else zero" default models a
fresh, uniform virgin slab at room temperature — the standard cone/TGA starting
point. The uniform-IC variant `apply_initial_conditions!(mesh, material; T, ξ)`
follows the same rule when $\xi=\texttt{nothing}$, with the first solid's density
evaluated at the given $T$ rather than at 300 K (`initial_conditions.jl`). Profile helpers (`step_profile`,
`linear_profile`, `gaussian_profile`, `tanh_profile`) are provided for non-uniform
seeds (e.g. a pre-heated layer).

### 12.8.3 The dry-solid mass denominator $m_{\text{dry},i}$

A second, BC-adjacent quantity is fixed at initialization: the per-cell initial
dry-solid mass, which is the constant denominator of the pyrolysis-progress state
$\chi_i = m_{\text{released,dry},i}/m_{\text{dry},i}$ used by the variable
cross-section law (§10, Volume Change & Lateral Shrinkage). It is populated by
`populate_m_dry_solid_initial!` (`initial_conditions.jl`), wired into
`apply_initial_conditions!(problem)` after the cell-state write:

$$
m_{\text{dry},i} = V_i(0)\,\sum_{j\in\text{solid}} \xi_{j,i}(0),
$$

summing only the solid components (the `is_solid` predicate — not liquids, not
gases) at the initial state. This is computed once and stored on the mutable mesh
as `Float64`; under AD the $\chi$ numerator is Dual-typed but this denominator
stays `Float64` (see §15, §17).

---

## 12.9 Boundary-flux diagnostics

### 12.9.1 Flux decomposition

For auditing and post-processing, `compute_boundary_heat_fluxes(bc, T_s, t,
ε_eff, radiation_model[, face])` decomposes the boundary flux into named
components (`boundary_fluxes.jl`). The optional trailing `face` is the
`BlowingFace` at the converged $T_s$ (outward species fluxes, surface
temperature, material); when supplied, the convective component is reported
with the same ratio $\chi$ its blowing model gave the surface Newton
(§12.3.3), keeping the decomposition consistent with `net`. Components:

- `absorbed` — incident radiation absorbed *at the boundary face* [W·m⁻²];
- `reradiation` — re-radiation loss (positive = leaving);
- `convection` — convective loss (positive = leaving);
- `prescribed` — prescribed flux;
- `net` — the flux the surface energy balance actually applies at the face
  (positive into the domain).

For `ConvectiveBC` and `RadiativeBC`, the loss components are reported with the
*opposite* sign to the into-domain convention (so a hot surface re-radiating shows
a positive `reradiation`). Per component, `net` matches
`evaluate_bc_flux_and_derivative` — the Newton-applied flux — for every
radiation model, `DirichletThermalBC` excepted (its decomposition is all zero; for a
bare Dirichlet BC the saved extras report the solved $q_{\text{conducted}}$ as the
net flux), so the ledger identity
$\texttt{net} = q_{\text{cond}}$ holds at convergence
(§12.9.2, §16.6). For combined BCs the components **and** `net` are summed over
children; `net` is never recomputed from the components. Consequences of the
model dispatch:

- `NO_RADIATION` / `SURFACE_ABSORPTION`: a `RadiativeFluxBC` contributes
  `absorbed = net = α·q_ext` (α the explicit BC absorptivity or a_eff) — the
  flux is applied at the face.
- `BEER_LAMBERT` / P1: a `RadiativeFluxBC` contributes zero `absorbed` and zero
  `net` — the radiation enters volumetrically (§8.3) and is reported separately
  as `extras.volumetric_absorption` $= \alpha_{\text{surf}}\,q_{\text{ext}}\,
  (1-e^{-\tau_{\text{tot}}})$, the exact telescoped sum of the Beer–Lambert
  deposit on the saved state.

### 12.9.2 Energy-balance check

`check_boundary_energy_balance(mesh)` returns, per boundary, the relative
imbalance $|q_{\text{applied}}-q_{\text{conducted}}|/|q_{\text{applied}}|$ (the
absolute difference when $|q_{\text{applied}}|\le10^{-10}$), which should be
$\approx$ machine precision ($\sim10^{-15}$) when the
Newton solve has converged (`boundary_temperature.jl`). This is a direct
continuous self-check that the surface energy balance $q_{\text{BC}}=q_{\text{cond}}$
holds. The accessors `get_surface_temperature(mesh)` and
`get_back_face_temperature(mesh)` return the solved $T_s$ at the top (id 1) and
bottom (id 2) boundaries (`boundary_temperature.jl`).

---

## 12.10 Worked example: the cone-calorimeter top boundary

Assembling the canonical top boundary makes the pieces concrete. A 50 kW·m⁻² cone
with surface re-radiation and natural convection to a 300 K ambient, on a
sealed-but-venting slab, is

```julia
bc_set = BoundaryConditionSet(
    thermal = (
        top    = RadiativeFluxBC(flux = 50e3)            # absorbed incident, α = ε_eff
               + RadiativeBC(ε = nothing, T_inf = 300)   # re-radiation loss
               + ConvectiveBC(h = 10.0, T_inf = 300),    # convective loss
        bottom = AdiabaticBC(),                          # insulated back
    ),
    mass = (
        top    = DiffusiveMassBC(),                      # gas vents, internal diffusion limit
        bottom = ImpermeableBC(),                        # sealed back
    ),
)
```

With `radiation_model = SURFACE_ABSORPTION`, the surface energy
balance solved by Newton at the top is

$$
\underbrace{\alpha\,(50{,}000)}_{\text{absorbed}}
\;-\;\underbrace{\varepsilon_{\text{eff}}\,\sigma\,(T_s^4 - 300^4)}_{\text{re-radiation}}
\;-\;\underbrace{10\,(T_s - 300)}_{\text{convection}}
\;=\;
k_{\text{eff}}\,\frac{T_s - T_{\text{cell}}}{\Delta z_{1/2}},
$$

with $\alpha=\varepsilon_{\text{eff}}$ (Kirchhoff, both `nothing`). As the surface
chars, $\varepsilon_{\text{eff}}$ tracks the composition automatically, raising
both the absorbed and re-radiated terms. The Newton derivative is
$f'(T_s) = -4\varepsilon_{\text{eff}}\sigma T_s^3 - 10 - k_{\text{eff}}/\Delta
z_{1/2}$ (the absorbed term and convection-ambient are $T_s$-independent), which is
strictly negative — guaranteeing a unique root reached in a few iterations. The
conducted flux $q_{\text{conducted}}$ is written to the top face as
$q[f_{\text{top}}]=-q_{\text{conducted}}$, and the vented gas leaves the top with
$\dot m_g$ accumulated by `DiffusiveMassBC` and oriented to $+z$ via the top-face
sign flip. If `radiation_model = BEER_LAMBERT` is selected instead, the
`RadiativeFluxBC` contributes zero to this balance and the 50 kW·m⁻² enters the
energy equation as a volumetric Beer–Lambert source $q'''=\alpha I$ (§8).

---

## 12.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The three reference solvers express the same surface physics in their own
notation; the subsections below map their symbols to ours on first use and
contrast the modeling choices. (ThermaKin/Gpyro publish $h>0=$ exothermic; our
internal convention is $h>0=$ endothermic with $Q_{\text{rxn}}=-\sum h_r r_r$, but
this sign issue is confined to the volumetric reaction heat, not the surface BCs.)

### 12.11.1 ThermaKin2Ds

ThermaKin2Ds applies piecewise external/flame boundary conditions and supports an
external radiant source with surface re-radiation. The exposed-face energy
balance, in ThermaKin's form, is

$$
-k\frac{\partial T}{\partial x}\Big|_{\text{surf}}
= \alpha\,\dot q''_{\text{ex}}
  - \varepsilon\sigma\bigl(T_s^4 - T_\infty^4\bigr)
  - h_{\text{conv}}\bigl(T_s - T_\infty\bigr),
$$

with $\dot q''_{\text{ex}}$ the incident external flux, $\alpha$ the absorptivity
and $\varepsilon$ the emissivity (Kirchhoff $\alpha=\varepsilon$). Pyrolysis.jl's
default cone BC `RadiativeFluxBC + RadiativeBC + ConvectiveBC` is exactly this
balance, with our $q_{\text{ext}}=\dot q''_{\text{ex}}$ and our
$\varepsilon_{\text{eff}}$ playing the role of ThermaKin's surface $\varepsilon$.
ThermaKin also tracks a reflected component $I_{\text{ex}}^0$; our
composition-weighted $\varepsilon_{\text{eff}}$ supplies
$\alpha=\varepsilon_{\text{eff}}$ and the reflected fraction $(1-\alpha)$ is
implicitly discarded (opaque-surface assumption), consistent with the
SURFACE_ABSORPTION model. The shared idea — that both the absorbed external flux
and the surface re-radiation depend on the *true* surface temperature — is
realized by ThermaKin through its Crank–Nicolson surface node and by Pyrolysis.jl
through the explicit Newton surface-balance solve of §12.2. ThermaKin carries gas
enthalpy via a convective term $c_g(J_g\cdot\partial T/\partial x)$ (its Eq. 17),
which is precisely our $S_{\text{conv}}$ and the reason the surface balance carries
no gas-enthalpy term (§7.4.3).

### 12.11.2 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) writes the exposed-surface balance as a
radiative–convective balance with an impermeable back face by default. In Gpyro
notation,

$$
-\bar k\frac{\partial T}{\partial z}\Big|_{\text{surf}}
= \dot q''_e - \varepsilon\sigma\bigl(T_s^4 - T_\infty^4\bigr)
  - h_c\bigl(T_s - T_\infty\bigr),
$$

and gas escape at the surface uses a heat–mass transfer analogy for the film
coefficient, $h_m \sim h_c/(\rho c_p)$ — the same estimate we recommend for
`ConvectiveMassBC` (§12.4.3). Gpyro's control-volume formulation places the
surface temperature at the face via its harmonic interface interpolation; our
Newton solve achieves the equivalent by solving the one-sided Fourier balance for
$T_s$. Gpyro permits a non-zero back-face pressure/mass boundary in its
pressure-evolution model; our `DirichletPressureBC` / `NeumannPressureBC` /
`ConvectivePressureBC` (§12.5) are the analogous Darcy boundary closures, with
Gpyro's permeability $K$ mapping to our $\kappa$. Gpyro's default impermeable back
face corresponds to our default `ImpermeableBC` + `NeumannPressureBC()` (zero
gradient).

### 12.11.3 FDS

The FDS solid-phase model imposes the 1D surface boundary condition

$$
-k\frac{\partial T_s}{\partial x}\Big|_{x=0}
= \dot q''_c + \dot q''_r,
$$

i.e. the conducted flux equals the sum of the convective and net radiative surface
fluxes, with $\dot q''_r$ including absorbed incident and emitted re-radiation.
This is the same content as our $q_{\text{cond}}=q_{\text{BC}}$ with
$q_{\text{BC}}=q_{\text{ext}}-\varepsilon\sigma T_s^4-h_{\text{conv}}(T_s-
T_\infty)$. FDS additionally performs a *surface oxygen* Newton iteration for
oxidative pyrolysis (coupling surface $O_2$ to the gas phase); Pyrolysis.jl is a
condensed-phase-only solver with no gas-phase $O_2$ feedback, so our surface
Newton iterates on $T_s$ alone. FDS recovers Arrhenius parameters from TGA
reference temperature/rate (a kinetics, not a BC, concern; see §6, §18). The FDS
in-depth radiation (a two-flux Schuster–Schwarzschild model) corresponds to our
BEER_LAMBERT / P1 volumetric handling, under which the `RadiativeFluxBC` boundary
flux is suppressed to avoid double-counting (§12.3.5), just as FDS deposits the
in-depth absorbed fraction volumetrically rather than at the face.

---

## 12.12 Limitations and open issues

1. **Temperature-dependent emissivity derivative.** For a `RadiativeBC` whose
   own `ε` is a function of $T_s$, the radiative Newton derivative omits the
   $F\,(d\varepsilon/dT)\,\sigma\,(T_\infty^4-T_s^4)$ term (§12.2.4). The
   converged $T_s$ is unbiased, but convergence can become sub-quadratic for
   strongly temperature-sensitive $\varepsilon(T_s)$. This does not affect
   component-level $\varepsilon_j(T)$: the mixture $\varepsilon_{\text{eff}}$
   is evaluated at the *cell* temperature and is a constant of the surface
   Newton iteration, so that derivative is exact.

2. **Dirichlet thermal BC via penalty.** `DirichletThermalBC` uses a large
   conductance $h_{\text{large}}=10^8$ in the Newton solve and a special-cased
   direct flux in the RHS (§12.3.6). A genuine algebraic constraint (DAE) would be
   cleaner and would remove the residual conditioning sensitivity of the penalty.

3. **View-factor meaning in 1D.** The `RadiativeBC` view factor $F$ is a scalar
   surrogate for a 3D geometric quantity; in a strictly 1D model it is best read as
   a tuning factor for partial enclosure rather than a rigorously computed exchange
   factor.

4. **Mass-transfer coefficient estimation.** The film coefficient $h_m$ for
   `ConvectiveMassBC` / `DiffusiveMassBC` is supplied by the user, typically via
   the heat–mass analogy $h_m\approx h_{\text{conv}}/(\rho_{\text{air}}
   c_{p,\text{air}})$; there is no built-in correlation, so its accuracy depends on
   the user's estimate of $h_{\text{conv}}$ and the flow regime.

5. **Pressure-BC scope.** The three pressure BCs act only under Darcy–Fick
   transport, i.e. when a condensed component is permeable (`has_darcy_flow`); on a
   Fickian material they are accepted and contribute nothing. Their venting is
   outflow-only: the ambient ghost carries none of the tracked gas species, so an
   exterior pressure above the interior pressure draws no gas into the slab
   (§9.5.3). The residual, the Jacobian and the saved-output diagnostics evaluate
   the one kernel `boundary_darcy_velocity` (§12.5).


---

# 13. Finite-Volume Spatial Discretization

This chapter derives the cell-centered finite-volume (FVM) spatial discretization that converts the continuous coupled PDE system of §3 (Governing Equations) into the spatially discrete semi-discrete system of ordinary differential / differential-algebraic equations integrated in §15 (Temporal Integration & Structured Jacobian). It documents the 1D mesh layout (nodes, cells, faces), the geometric quantities that enter every flux, the divergence operator and its sign conventions, the conductive / gas-diffusive / Darcy–Fick face-value formulas, the TVD flux limiters and gradient-reconstruction machinery, the moving-mesh (ALE) advection discretization, mesh generation (uniform and geometric-stretch), and the block-major state-vector packing with allocation-free geometry-from-state. Every equation, sign, factor, guard, and limiter below has been confirmed against the implementation; where the code deviates from the textbook form (or from the design-document prose), the code is reported as ground truth and the deviation is called out.

The unifying design principle is that **all transport is expressed as fluxes evaluated at cell faces**, and the time derivative of every conserved or quasi-conserved cell quantity is the discrete divergence of those face fluxes plus volumetric sources. This guarantees discrete conservation up to the boundary-flux treatment (see §16, Conservation Diagnostics) and produces a sparse, banded coupling structure that the structured Jacobian of §15 exploits.

Throughout this chapter the index `i` is the **cell index** (overload resolution I1: `i` = cell in discretization; the reaction index is `r`), `j` is the **component index** (`1…N_C`), and `f` is the **face index**. The spatial coordinate is `z` with `z = 0` at the bottom/substrate (boundary id 2, tag `:bottom`) and `z = L` at the top/exposed surface (boundary id 1, tag `:top`); heat enters from the top and the material recedes downward. The three sign conventions of §2 (Nomenclature) are in force: positive flux is transport in `+z`; the divergence `(F_R − F_L)A/V` is positive for net outflow; and the heat-of-reaction convention `h > 0` endothermic, `Q_rxn = −Σ_r h_r r_r`.

---

## 13.1 Cell-Centered Finite-Volume Layout in 1D

### 13.1.1 Topology: nodes, cells, faces

Pyrolysis.jl stores its primary unknowns — temperature `T` and the mass concentrations `ξ_j` — at **cell centers**, not at faces or nodes. This is the classic cell-centered (collocated) finite-volume arrangement. The 1D mesh consists of three topological entity types whose counts are rigidly linked:

$$
n_{\text{nodes}} = n_{\text{faces}} = n + 1, \qquad n \equiv n_{\text{cells}}.
$$

- A **node** (`Node1D`) is a mesh vertex at position `z_i`. Node 1 sits at `z = 0` (bottom, `boundary_id = 2`), node `n+1` sits at `z = L` (top, `boundary_id = 1`); the rest are interior (`boundary_id = 0`). In ALE mode the node positions become part of the ODE state (§13.10).
- A **cell** (`Cell1D`) is the finite volume. Cell `i` spans the two adjacent nodes `(i, i+1)` and is bounded by the two faces `(i, i+1)`. Cells carry the unknowns through their `CellState{NC}` (`T`, the `NTuple{NC}` of `ξ_j`, and a diagnostic `volume`). Cells are immutable in connectivity; depletion parking (§14.3) changes only the live cell count `n_active_cells`, never the arrays.
- A **face** (`Face1D`) is, in 1D, a single point co-located with a node. Face `f` carries its `cell_ids = (left, right)` with `0` marking a domain boundary: the bottom face (id 1, at `z = 0`) has `cell_ids = (0, 1)`; the top face (id `n+1`, at `z = L`) has `cell_ids = (n, 0)`; interior face `i` has `cell_ids = (i−1, i)`.

"Left" always denotes the lower-`z` (toward substrate) side of a face and "Right" the higher-`z` (toward surface) side. Cell indices increase from bottom to top, so node `i` is the left endpoint of cell `i` and the right endpoint of cell `i−1`; face `i` co-locates with node `i`.

> **Implementation note.** The data structures are defined in `Node1D` (src/Geometry/nodes.jl), `Cell1D`/`CellState` (src/Geometry/cells.jl), `Face1D` (src/Geometry/faces.jl), and the container `Mesh1D{NC}` (src/Geometry/mesh.jl). The connectivity view `MeshTopology` and the per-step `MeshGeometry` (src/Geometry/topology.jl, src/Geometry/mesh_geometry.jl) are placeholders for a multi-dimensional extension: the 1D residual reads the `Mesh1D` connectivity and the state-derived geometry held in the typed `StateCache` (§13.10.2), which is what makes pure (mesh-non-mutating) residual and Jacobian evaluation possible.

### 13.1.2 Boundary identifier convention

The internal representation uses integer boundary ids; the user-facing API uses symbols. The bridge is

$$
\texttt{tag\_to\_boundary\_id}: \;\; :\texttt{top} \mapsto 1, \quad :\texttt{bottom} \mapsto 2, \quad \text{otherwise} \mapsto 0,
$$

with the inverse `boundary_id_to_tag` (src/Geometry/topology.jl). The boundary-face map `mesh.boundary_faces` stores `1 => [n+1]` (top) and `2 => [1]` (bottom). This convention propagates into the boundary-temperature solver, the mass/pressure/mesh boundary conditions of §12, and the diagnostics of §16.

### 13.1.3 What a "cell-centered" scheme means for accuracy

Because unknowns live at cell centers but fluxes are needed at faces, every flux requires an *interpolation/reconstruction* from cell values to faces. For diffusion this is a face-property average plus a two-point gradient (§13.7); for advection it is an upwind or limited reconstruction (§13.4, §13.5). On a uniform mesh the central two-point gradient at a face is second-order; the cell-centered gradient (§13.6) is second-order in the interior and degrades to first-order one-sided at boundaries. The harmonic-mean face conductivity (§13.7) preserves flux continuity exactly across material interfaces, which is the physically correct behavior for a layered/heterogeneous medium and the main reason the cell-centered scheme is preferred over a vertex-centered one here.

---

## 13.2 Geometric Quantities

All discrete operators are built from a small set of geometric quantities derived from node positions and the (possibly time-varying) cross-sectional area `A(t)`.

### 13.2.1 Cell volume and cell center

For cell `i` spanning nodes `{z_i, z_{i+1}}`:

$$
V_i = |z_{i+1} - z_i|\, A, \qquad z_i^c = \tfrac{1}{2}\,(z_i + z_{i+1}).
$$

The absolute value guards against transiently inverted nodes during ALE motion; the cross-sectional area `A` is the **current** column area (it equals `A_0` unless a lateral-shrinkage law is active — see §10, Volume Change). All faces in 1D carry the same area `A_f = A`.

> **Implementation note.** `update_geometry!` (src/Geometry/geometric_quantities.jl) writes the mesh-resident `cell_centers[i]`, `cell_volumes[i] = abs(z_right − z_left) * A`, and `face_areas[f] = A`; it is used between integrator steps. During residual/Jacobian assembly the same quantities are rebuilt from the state by `update_state_geometry!` (src/Problem/state_cache.jl; §13.10).

### 13.2.2 Interior face distance

Diffusive fluxes use the **center-to-center distance** across a face, not the cell width. For an interior face `f` with neighbors `L` and `R`:

$$
d_{LR} = |z_R^c - z_L^c|.
$$

This is the denominator of every two-point gradient in §13.7 and is the correct choice for second-order accuracy on non-uniform meshes (using cell width instead would introduce a zeroth-order multiplicative bias).

### 13.2.3 Boundary distances

At a boundary face the "other side" cell is absent. Every boundary flux is formed over the **center-to-face** distance of the adjacent cell,

$$
d_f^{\text{bnd}} = |z_f - z^c_{\text{adj}}| = \tfrac12\,\frac{V_{\text{adj}}}{A},
$$

the one-sided gradient from the surface to the first cell center: the surface energy balance (§12.2.1), the film-law mass BCs (§12.4), and the pressure BCs (§12.5) all use it. The interior conductive pass zeroes the boundary faces (§7.1.3), so no other boundary distance enters the diffusive fluxes.

### 13.2.4 The ALE advection boundary distance

The ALE advection operator (§13.8) needs an *upwind* distance, computed by `_upwind_gradient_distance` (src/Discretization/ale_fluxes.jl). For an interior upwind neighbor it returns the center-to-center distance, consistent with §13.2.2. For a boundary (`neighbor_id == 0`) it returns

$$
d^{\text{bnd}}_{\text{upwind}} = \frac{V_i}{A} = \Delta z_i,
$$

i.e. the **full cell width**: the spacing from the cell center to a mirror ghost cell reflected across the boundary face, twice the center-to-face distance of §13.2.3 — the same footing as the interior center-to-center distance. At a boundary the value is in any case inert: the zero-gradient ghost makes the first-order gradient identically zero in that direction, so no residual or Jacobian entry depends on it (§11.6.3).

> **Implementation note.** The ALE state-vector geometry lives in the typed `StateCache` (`volume`, `cell_centers`, `face_positions`, `A_face`, `A_section`), rebuilt from the node positions by `update_state_geometry!` at every residual evaluation (§13.10.2). Face distances are not stored; each kernel forms the center-to-center or center-to-face distance it needs from `cell_centers` and `face_positions`.

---

## 13.3 The Divergence Operator

### 13.3.1 Discrete divergence of a scalar flux field

The control-volume integral of a divergence reduces, via the divergence theorem, to a sum of outward face fluxes. In 1D with a left face (outward normal `−ẑ`) and a right face (outward normal `+ẑ`):

$$
(\nabla\!\cdot\!\mathbf F)_i \;=\; \frac{1}{V_i}\oint_{\partial V_i}\mathbf F\cdot \hat{\mathbf n}\,\mathrm dS
\;\xrightarrow{\text{1D}}\;
\boxed{\;\text{div}[i] = \frac{(F_{f_R} - F_{f_L})\,A}{V_i}\;}
$$

where `F_{f_L}`, `F_{f_R}` are the face fluxes at the left and right faces of cell `i`, `A` is the (common) face area, and `V_i` is the cell volume. Because positive flux is `+z` transport, the right-face flux is outflow `(+)` and the left-face flux is inflow `(−)`; a positive divergence therefore means **net efflux from the cell**, which in the energy equation enters with a leading minus sign so that net efflux cools the cell.

> **Implementation note.** `divergence_1d!` (src/Discretization/finite_volume.jl) reads `V_i = state.volume[i]` and `A = state.A_face[f_left]` and computes `(fluxes[f_right] − fluxes[f_left]) * A / V_i`. Both `V_i` and `A` carry the live scalar type, so in ALE/variable-cross-section mode the divergence is differentiable through the geometry (essential for the Jacobian and for the sensitivities of §17).

### 13.3.2 Multi-component (species) divergence

For the gas mass fluxes, the same operator is applied component-wise to the `NC × n_faces` flux matrix `J`, producing the `NC × n_cells` divergence matrix:

$$
\text{div}_j[i] = \frac{(J_{j,f_R} - J_{j,f_L})\,A}{V_i}, \qquad j = 1,\dots,N_C.
$$

> **Implementation note.** `divergence_species_1d!` (src/Discretization/finite_volume.jl) loops over cells and components with the identical formula. Solid/liquid components carry zero face flux (they are mesh-bound, §13.7.4), so their divergence is identically zero and the species residual reduces to reactions only (§13.3.4).

### 13.3.3 Discrete conservation property

Summing the scalar divergence over all cells telescopes:

$$
\sum_{i=1}^{n} V_i\,(\nabla\!\cdot\!\mathbf F)_i = A\sum_{i=1}^{n}(F_{f_R(i)} - F_{f_L(i)}) = A\,(F_{\text{top face}} - F_{\text{bottom face}}),
$$

so the integral of the divergence equals the net boundary flux exactly — the discrete analog of Gauss's theorem. This identity is the foundation of the conservation diagnostics in §16 and is exercised by the unit tests (constant flux → zero divergence; linear flux → constant divergence; integral closure to boundary fluxes; test/discretization/test_discretization.jl, test/discretization/test_conservation.jl).

### 13.3.4 Sign assembly in the final residual

The cell-wise temperature derivative assembled at the end of `compute_rhs!` is

$$
\frac{dT_i}{dt} = \frac{1}{\rho c_{p,\text{eff},i}}\Bigl[\,-(\nabla\!\cdot\!\mathbf q)_i + S_{\text{conv},i} + Q_{\text{rad},i} + Q_{\text{rxn},i}\,\Bigr],
$$

and the species derivative (per cell, per component) is

$$
\frac{d\xi_{j,i}}{dt} =
\begin{cases}
-(\nabla\!\cdot\!\mathbf J_j)_i + S_{j,\text{rxn},i}, & \text{component } j \text{ is gas},\\[4pt]
S_{j,\text{rxn},i}, & \text{component } j \text{ is solid/liquid}.
\end{cases}
$$

The leading `1/ρc_{p,\text{eff}}` uses the mixture volumetric heat capacity (every component; see §7, Heat Conduction & Energy Assembly, and §13.7.5). The gas/solid dispatch is exactly the `is_gas`-branch in `apply_species_residual_unrolled!` (src/Materials/materials.jl): gas → `-div_J + dξ_rxn`, solid/liquid → `dξ_rxn`. The full eight-step assembly order is given in §7.6; this chapter is concerned only with the spatial operators that produce `div_q`, `div_J`, and the face fluxes feeding them.

---

## 13.4 Flux Limiters (TVD Schemes)

Advective transport with a centered reconstruction is second-order but oscillatory near steep fronts; a pure upwind reconstruction is monotone but only first-order and excessively diffusive. Total-Variation-Diminishing (TVD) flux limiters blend the two: they recover second order in smooth regions and degrade gracefully to first order near extrema, provably not increasing the total variation `Σ|φ_{i+1}−φ_i|`. Pyrolysis.jl uses TVD limiting only in the mesh-relative temperature advection of the ALE verification hook (§13.8).

### 13.4.1 The limiter enumeration

A small `@enum Limiter` selects the scheme (src/Discretization/flux_limiters.jl):

| Enum | `ψ(r)` | Character |
|---|---|---|
| `MINMOD` | `max(0, min(1, r))` | Most diffusive; strictly TVD |
| `VANLEER` | `(r + \|r\|)/(1 + \|r\|)` | Balanced; **default in ALE** |
| `SUPERBEE` | `max(0, min(2r, 1), min(r, 2))` | Least diffusive; compressive |
| `NONE` | `0` | First-order upwind (no limiting) |
| `CENTRAL` | `1` | Second-order central (no limiting) |

The dispatcher `apply_limiter(r, limiter)` returns the corresponding `ψ(r)`.

### 13.4.2 Limiter functions and the TVD region

The slope ratio `r` measures the smoothness of the field across the upwind direction (defined in §13.5). Each limiter is a function `ψ(r)` of that ratio:

$$
\psi_{\text{minmod}}(r) = \max(0,\min(1,r)), \qquad
\psi_{\text{vanLeer}}(r) = \frac{r + |r|}{1 + |r|}, \qquad
\psi_{\text{superbee}}(r) = \max\bigl(0,\min(2r,1),\min(r,2)\bigr).
$$

All three satisfy `ψ(r) = 0` for `r ≤ 0` (extremum → first-order upwind), `ψ(1) = 1` (smooth → second order), and lie within Sweby's TVD region `0 ≤ ψ(r) ≤ min(2r, 2)`. Minmod hugs the lower edge of the region (maximal dissipation, maximal robustness); superbee hugs the upper edge (minimal dissipation, mild steepening/compression that can clip smooth extrema); van Leer is the smooth symmetric interior choice — it is the analytic harmonic mean of `r` and `1`, is differentiable for `r ≠ 0` (which matters for the AD-based Jacobian and sensitivities of §15/§17), and is the reason it is the package default.

> **Implementation note.** `minmod`, `vanleer`, `superbee`, `no_limiter`, `central_limiter`, and the `apply_limiter` dispatcher are all `@inline` scalar functions (src/Discretization/flux_limiters.jl). The hard `max`/`min` in minmod/superbee make those two non-differentiable at the kink points; van Leer is preferred when smooth derivatives are required end-to-end.

---

## 13.5 Limited Reconstruction: What Is Live

The TVD machinery in `src/Discretization/flux_limiters.jl` consists of the
limiter functions $\psi(r)$ (minmod, van Leer, superbee; §13.8.2) and the
`Limiter` enum / `apply_limiter` dispatcher. Their consumer is the mesh-relative
temperature advection's cell-centered limited gradient `_limited_gradient`
(§11.6.2, §13.8.2), which reconstructs $\partial f/\partial z$ donor-side with the
distance-correct blend $\beta = d_0/(d_0+d_1)$ and runs only under the
`w_override` verification hook.

There is no limited *face-value* reconstruction in production: the Eulerian
gas-advection face value is first-order upwind on the Darcy velocity
(`gas_flux_darcy_fick`, §9; a TVD upgrade is a documented accuracy-debt item),
and conduction/diffusion use the exact two-point face gradients of §13.7.

---

## 13.6 Gradient Reconstruction

### 13.6.1 Cell-centered gradient on non-uniform meshes

The gas-advection source `S_conv` of §7.4 and the error indicators and $r$-refinement monitors of §14 require a cell-centered gradient. The interior formula combines the two one-sided center-to-center slopes with **opposite-distance weights** (`opposite_distance_gradient`, src/Geometry/geometric_quantities.jl):

$$
\nabla\varphi_i = \frac{d_R\,s_L + d_L\,s_R}{d_L+d_R},
\qquad
s_L = \frac{\varphi_i - \varphi_{i-1}}{d_L},\quad
s_R = \frac{\varphi_{i+1} - \varphi_i}{d_R},
$$

with $d_L = z_i^c - z_{i-1}^c$ and $d_R = z_{i+1}^c - z_i^c$. Each one-sided slope estimates $\varphi'$ at its segment midpoint; the opposite-distance weighting interpolates them back to $z_i^c$, cancelling the leading $\varphi''$ truncation term — second order on arbitrary nonuniform meshes. On a uniform mesh this reduces to the familiar central difference $(\varphi_{i+1}-\varphi_{i-1})/(2\Delta z)$. The wide-stencil form $(\varphi_{i+1}-\varphi_{i-1})/(z_{i+1}^c-z_{i-1}^c)$ — which weights each slope by its *own* distance — carries a first-order $\tfrac12(d_R-d_L)\varphi''$ error on stretched meshes and is not used.

### 13.6.2 One-sided boundary gradients

At domain edges, one neighbor is missing, so the gradient degrades to a one-sided difference:

$$
\nabla\varphi_i =
\begin{cases}
\dfrac{\varphi_{i+1} - \varphi_i}{z_{i+1}^c - z_i^c}, & \text{bottom cell (left\_id = 0)},\\[10pt]
\dfrac{\varphi_i - \varphi_{i-1}}{z_i^c - z_{i-1}^c}, & \text{top cell (right\_id = 0)},\\[10pt]
0, & \text{isolated cell (both ids = 0; degenerate)}.
\end{cases}
$$

> **Implementation note.** `compute_gradient_1d(φ, mesh, i)` (src/Adaptivity/error_indicators.jl) implements these branches on the mesh-resident cell centers for the error indicators and the $r$-refinement monitors (src/Adaptivity/error_indicators.jl, r_refinement.jl). The interior `opposite_distance_gradient` also serves the `S_conv` temperature gradient (src/Physics/convection.jl), which reads the state's `cell_centers` so the gradient carries AD sensitivities through the geometry and upgrades the boundary cells to the surface-slope blend of §7.4.2 (the thermal BC supplies a solved $T_s$ there; a generic field has no boundary value, so the one-sided fallback above is the right one for `compute_gradient_1d`).

---

## 13.7 Conductive, Gas-Diffusive, and Darcy–Fick Face Values

This section states the discrete face fluxes as they appear in the FVM. The physical derivations are in §7 (conduction, convection) and §9 (gas transport); here the focus is on the *discrete* form, the face averaging, and the signs.

### 13.7.1 Conductive heat flux

Fourier's law at a face, discretized with a distance-weighted harmonic-mean conductivity:

$$
q_f = -k_f\,\frac{T_R - T_L}{d_L + d_R}, \qquad
k_f = \frac{d_L + d_R}{\dfrac{d_L}{k_L} + \dfrac{d_R}{k_R}}.
$$

The harmonic mean is the **resistance-in-series** average: across an interface the thermal resistances `d_L/k_L` and `d_R/k_R` add, and continuity of heat flux then fixes `k_f` as above. On a uniform mesh `d_L = d_R` and this collapses to the standard harmonic mean `2k_Lk_R/(k_L+k_R)`; on stretched meshes the distance weighting improves near-boundary flux accuracy by roughly 1–5%. The minus sign is Fourier's law: `T_R > T_L` gives `q_f < 0` (heat flows toward `−z`, down the gradient).

> **Implementation note.** `conductive_flux(k_face, T_L, T_R, d_LR) = -k_face*(T_R−T_L)/d_LR` (src/Physics/heat_conduction.jl) consumes the face conductivity `update_properties!` pre-computes with `harmonic_mean_weighted`. The assembly `compute_conductive_fluxes!` loops all faces using the precomputed `cache.k_face[f]` and `d_LR = |z_R^c − z_L^c|`, and **explicitly zeros boundary faces** (`left_id == 0 || right_id == 0 → q_cond[f] = 0`). The zero is then overwritten by the surface-energy-balance Newton solve (§12); zeroing first prevents double-counting (the interior operator must not see a boundary flux that the BC solver will also supply).

### 13.7.2 Gas-diffusive (volume-fraction-gradient) flux

The gas diffusive flux is written in **volume-fraction (equivalently partial-pressure) gradient** form. Starting from `J_g = −ρ_g λ ∂(ξ_g/ρ_g)/∂x` with the ideal-gas relation `ρ_g ∝ P/T`, this collapses to `J = −(λ/T)·∂(ξT)/∂x`, whose discrete face form is

$$
\boxed{\;J_f = -\lambda_f\,\frac{\xi_R T_R - \xi_L T_L}{T_f\,d_{LR}}\;}, \qquad T_f = \tfrac12(T_L + T_R).
$$

The product `ξT` is proportional to partial pressure (`P_j = ξ_j R_g T / M_j`), so this is diffusion down the partial-pressure gradient. The algebraic identity

$$
\xi_R T_R - \xi_L T_L = T_f\,(\xi_R - \xi_L) + \xi_f\,(T_R - T_L)
$$

shows the flux contains both an ordinary-diffusion term (`∝ ∂ξ/∂x`) and a thermal-diffusion / Soret-like term (`∝ ξ/T·∂T/∂x`); the single combined expression captures both without ever forming `ρ_g` at the face.

Crucially, the **harmonic mean is taken of `λ` alone, not of `λρ`**. The density dependence is carried implicitly by the mass concentrations `ξ` in the numerator, so averaging `λ` (a transport resistance) is the correct interface treatment:

$$
\lambda_f = \frac{d_L + d_R}{\dfrac{d_L}{\lambda_L} + \dfrac{d_R}{\lambda_R}}.
$$

> **Implementation note.** `gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)` (src/Physics/gas_transport.jl) returns exactly `-λ_face*(ξ_R*T_R − ξ_L*T_L)/(T_face*d_LR)` with `T_face = (T_L+T_R)/2`. The face coefficient `λ_face` is the distance-weighted harmonic mean from `harmonic_mean_weighted` (src/Properties/effective_properties.jl). All gas species share the same `λ_face`; there is no per-species scaling of the diffusive term (§9.2).

### 13.7.3 Darcy–Fick combined advection–diffusion flux

When pressure-driven Darcy flow is active, the face flux adds a first-order-upwinded advective term to the diffusive term:

$$
J_f^{\text{DF}} = \underbrace{\frac{\xi^{\text{upwind}}}{\varphi_f}\,v_g}_{\text{advection}} \;+\; \underbrace{\Bigl(-\lambda_f\,\frac{\xi_R T_R - \xi_L T_L}{T_f\,d_{LR}}\Bigr)}_{\text{diffusion}},
\qquad
\xi^{\text{upwind}} =
\begin{cases}
\xi_L, & v_g \ge 0,\\
\xi_R, & v_g < 0,
\end{cases}
$$

with $\varphi_f=\max(\tfrac12(\varphi_L+\varphi_R),10^{-12})$ the face porosity: $\xi$ is per bulk volume and $v_g$ is superficial, so the advected scalar is the intrinsic density $\xi/\varphi$ (§9.3). The upwinding chooses the cell *from which* the bulk gas flows, which is the monotone first-order choice. The bulk Darcy velocity at the face is

$$
v_g = -\text{mob}_f\,\frac{P_R - P_L}{d_{LR}},\qquad
\text{mob}_f = \frac{d_L+d_R}{d_L\mu_L/\kappa_L + d_R\mu_R/\kappa_R},
$$

the **distance-weighted harmonic mean of the cell mobilities** $\kappa/\mu$, with `mob_f = 0` exactly when either side's mobility is below $10^{-30}$ (an impermeable side is an infinite series resistance) and `μ > 10^{-20}` guards per cell. The ideal-gas pressure closure `P = (Σ_g ξ_j R_g T/M_j)/φ` (porosity floor `10^{-12}`) and the pressure gradient are derived in §9; here they enter only as the velocity driving the advective face flux. All gas species share the same `v_g` (advection carries every species equally) and the same `λ_face` (§9.2).

> **Implementation note.** `gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)` (src/Physics/gas_transport.jl) computes `J_diff` exactly as above, sets `ξ_upwind = v_g ≥ 0 ? ξ_L : ξ_R`, and returns `J_adv + J_diff` with `J_adv = v_g·ξ_upwind/φ_face`. The face velocity and mobility come from `update_darcy_properties!` (src/Discretization/finite_volume.jl). `compute_rhs!` takes this path when the property cache carries a Darcy sub-cache (`nothing` ⇔ Fickian): it refreshes pressure and velocity first, calls the Darcy–Fick flux instead of the diffusion-only flux, and writes the boundary-face advective venting flux from the pressure BCs (§9.5.3).

### 13.7.4 Solid/liquid components carry no transport flux

Solid and liquid components are mesh-bound (Lagrangian): they are not transported by diffusion or Darcy flow. Their face fluxes are identically zero, so `divergence_species_1d!` produces zero for them and the species residual reduces to reactions only (§13.3.4). This is what locks solid mass-thickness `ρ_s·Δz` to the mesh and makes the ALE motion (§13.8) the mechanism by which solids are advected. The gas fluxes of §13.7.2–§13.7.3 are defined relative to that skeleton, so under ALE the gas rides the mesh too and the same divergence is the whole of its transport (§11.1.2).

### 13.7.5 The energy-storage convention

The conductive divergence feeds the energy equation with the mixture heat capacity `ρc_{p,\text{eff}} = Σ_j ξ_j c_{p,j}(T)` in the denominator — every component, gas included, stores sensible heat. The enthalpy carried by the gas flux relative to the skeleton is the volumetric advection source `S_conv`, always on. This is the ThermaKin Eq. 17 convention. The discretization layer simply assembles `div_q`, `S_conv`, `Q_rad`, `Q_rxn` and divides by `prop_cache.ρcp[i]`. See §7.3–§7.4 for the full derivation, and §13.3.4 for the assembled formula.

---

## 13.8 ALE Flux Discretization (Recap)

This section summarizes how the moving mesh modifies the spatial operators; the full ALE theory (mesh-velocity construction, GCL) is in §11, The ALE Framework.

### 13.8.1 The ALE advection correction

The arbitrary Lagrangian–Eulerian time-derivative identity for a moving sample point is

$$
\frac{\partial f}{\partial t}\bigg|_{\chi} = \frac{\partial f}{\partial t}\bigg|_{x} + w\,\frac{\partial f}{\partial z},
$$

so a field whose capacity is carried by material moving at velocity $v$ picks up the advection term $(w - v)\cdot\nabla f$ in the mesh frame. Every component rides the skeleton and the mesh moves with it (§11.1.2, §11.2.2), so in production $w = v$ and no advection term is assembled. The term is computed only under the `w_override` verification hook, with $w$ replaced by the relative velocity $w_{\text{rel}} = w_{\text{applied}} - w_{\text{physical}}$, cell by cell, and **added** to the assembled Eulerian RHS:

$$
\frac{dT_i}{dt}\bigg|_{\text{ALE}} = \frac{dT_i}{dt}\bigg|_{\text{Euler}} + w_{\text{cell},i}\,\Bigl(\frac{\partial T}{\partial z}\Bigr)_i,
\qquad
w_{\text{cell},i} = \tfrac12\,(w_{n_L} + w_{n_R}).
$$

The cell velocity is the average of its two node velocities; node velocities are built by bottom-up integration of the per-cell dilation (`w_1 = 0` pinned at the substrate; §11.3).

### 13.8.2 Upwind, TVD-limited gradient for the advection term

Because the correction `w·∇φ` is *added* to the Eulerian RHS, the mesh-frame equation is `∂φ/∂τ + c·∂φ/∂z = RHS` with apparent transport velocity `c = −w` (the lab-frame field drifts at `−w` relative to the mesh). The advected gradient `∂φ/∂z` is therefore reconstructed **donor-side with respect to `c = −w`** — the right neighbor for `w > 0`, the left for `w < 0` — with an optional TVD limiter (default `VANLEER`). The first-order donor-side gradient, always written in the `+z` sense, is

$$
\text{grad}_1 =
\begin{cases}
(T_{\text{donor}} - T_i)/d, & w_{\text{cell}} > 0 \;(c<0,\ \text{donor on the right}),\\[4pt]
(T_i - T_{\text{donor}})/d, & w_{\text{cell}} < 0 \;(c>0,\ \text{donor on the left}),
\end{cases}
$$

and the second-order TVD correction, using the next-donor-side gradient `grad_0`, the **near/far** ratio `r = grad_1/grad_0` (see §11.6.2 for why the inverse ratio is unbounded), and the distance-correct blend `β = d₀/(d₀+d₁)` (½ on uniform meshes; §11.6.2), is

$$
\text{grad}_{\text{TVD}} = \text{grad}_1 + \beta\,\psi(r)\,(\text{grad}_1 - \text{grad}_0).
$$

The correction is dropped (first order only) when there is no donor-donor cell, when `limiter == NONE`, or when `|grad_0| < 10^{-20}` (flat far slope — ratio undefined). Boundary neighbors use zero-gradient (ghost = cell value) extrapolation; in the standard shrinking case (`w < 0`) the surface cell's donor is the interior neighbor, so the receding surface receives a real gradient rather than the boundary ghost's zero. A negligible mesh velocity `|w_cell| < 10^{-20}` short-circuits the advection to zero.

> **Implementation note.** `compute_ale_advection_term!` (energy) and `_limited_gradient` (src/Discretization/ale_fluxes.jl) implement these formulas; the assignment is `advection[i] = w_cell * dT_dz` and this is added to `dT` in `_apply_ale_T_advection!` (src/Residual/residual.jl) — **only under `Workspace.w_override`**, with `w = w_rel = w_applied − w_physical` (§11.2.2): production ALE has no energy advection because the matrix carrying the thermal mass moves with the mesh.

### 13.8.3 Species: dilation only

Every component rides the mesh (§11.1.2), so the species equations carry no mesh-motion advection at all. The Eulerian species RHS (`-div_J + dξ_rxn` for gases, `dξ_rxn` for solids/liquids) is used unchanged, and the residual adds the dilation `−ξ_{j,i}·θ_i` to every component (§11.2.3, §11.3.4): per cell, `d(ξ_j V)/dt = S_j V − [J_j A]_faces`, so each component's mass changes only by reaction and by its skeleton-relative flux through the faces, and the column total changes only through the boundary faces. In a sealed collapsing column the pore gas compresses cell by cell (V-CPL-03 P4 pins exact conservation); at an open face the gas leaves only through the mass BC's diffusive or Darcy flux.

> **Implementation note.** `_apply_ale_dilation_correction!` (src/Residual/residual.jl): `dξ[j, i] -= state.ξ[j, i] * θ[i]` for every component, with `θ` the cell's reaction-driven strain rate. Its Jacobian is the cell-local block of §15.8.4.

---

## 13.9 Mesh Generation

### 13.9.1 Uniform mesh

The simplest mesh divides `[0, L]` into `n` equal cells:

$$
\Delta z = \frac{L}{n}, \qquad z_i = (i-1)\,\Delta z \;\; (i = 1,\dots,n+1), \qquad V_i = \Delta z\,A.
$$

Node 1 (`z = 0`) is the bottom boundary (`boundary_id = 2`); node `n+1` (`z = L`) is the top (`boundary_id = 1`). Faces, cells, the active-cell mask (all `true`), the per-cell dry-solid-mass denominator `m_dry_solid_initial` (zeroed, populated later at problem setup), and the boundary-face map (`1 => [n+1]`, `2 => [1]`) are all initialized, and `update_geometry!` computes the cached centers/volumes/areas.

> **Implementation note.** `create_uniform_mesh(L, n_cells, Val(NC); T_initial, ξ_initial, cross_section_area)` (src/Geometry/mesh_generation.jl).

### 13.9.2 Geometric-stretch mesh

To resolve the steep thermal and reaction fronts near the heated surface without paying for fine resolution everywhere, a geometric stretch concentrates cells near one boundary. With a stretch ratio `s` (`Δz_{i+1}/Δz_i = s`), the cell widths form a geometric series whose sum is `L`:

$$
L = \sum_{i=1}^{n}\Delta z_i = \Delta z_0\,\frac{s^{\,n} - 1}{s - 1}
\;\;\Longrightarrow\;\;
\Delta z_0 = L\,\frac{s - 1}{s^{\,n} - 1}, \qquad \Delta z_i = \Delta z_0\,s^{\,i-1}.
$$

Node positions are accumulated `z_{i} = z_{i-1} + Δz_{i-1}` and then **renormalized** by `L/z_{n+1}` so that `z_{n+1} = L` exactly (removing floating-point drift in the partial sums). For `|s − 1| < 10^{-10}` the code falls back to the uniform `range(0, L, n+1)`.

The `stretch_from` keyword names the boundary that gets the **finest** cells. The series is always built fine-end-first from `z = 0` (smallest cell `Δz_0` at the bottom, cells growing toward `z = L`); with the default `stretch_from = :top` the positions are then reflected so the finest cells sit at the exposed surface `z = L`:

$$
\tilde z_i = L - z_{n+2-i}.
$$

`stretch_from = :bottom` keeps the unreflected series (finest cells at the substrate). The constructor requires `stretch ≥ 1` and `stretch_from ∈ (:top, :bottom)`; orientation is selected only by `stretch_from`, never by a sub-unity ratio.

> **Implementation note.** `create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from=:top, ...)` (src/Geometry/mesh_generation.jl). The geometric series uses `Δz_i = Δz_0 * stretch^(i-2)` in the node-accumulation loop (so the first interior gap is `Δz_0 * s^0 = Δz_0`), followed by the `positions .*= L/positions[end]` renormalization, and the `:top` reflection `positions = L .- reverse(positions)`. For convergence studies on stretched meshes see §18 (Verification & Validation), e.g. `convergence_douglas_fir`.

---

## 13.10 State-Vector Packing and Geometry-from-State

### 13.10.1 Block-major layout `[T | ξ | z | χ]`

The flat ODE/DAE state vector `u` is packed **block-major**: all temperatures first, then each component's concentration block, then (in ALE) the node positions, then (in WithChi) the pyrolysis-progress block. With `n = n_cells`:

$$
\underbrace{T_1,\dots,T_n}_{\text{T block: }1..n}\;\bigl|\;
\underbrace{\xi_{1,1},\dots,\xi_{1,n}}_{\xi_1}\,,\dots,\,\underbrace{\xi_{N_C,1},\dots,\xi_{N_C,n}}_{\xi_{N_C}}\;\bigl|\;
\underbrace{z_1,\dots,z_{n+1}}_{\text{z block (ALE)}}\;\bigl|\;
\underbrace{\chi_1,\dots,\chi_n}_{\text{$\chi$ block (WithChi)}}.
$$

The exact ranges (1-based) are:

$$
T:\;1..n, \qquad
\xi_j:\;n + (j-1)n + (1..n), \qquad
z:\;n(1+N_C) + (1..n_{\text{nodes}}), \qquad
\chi:\;\text{trailing }n.
$$

Block presence is determined **strictly** by the `SimulationMode` trait: `T` and `ξ` always present; `z` iff `mesh_mode(mode) isa ALE`; `χ` iff `chi_mode(mode) isa WithChi`. The total length is `n(1+N_C)` plus `n_nodes` (ALE) plus `n` (WithChi).

> **Implementation note.** `StateLayout{NC,Mode}` (src/Problem/state_layout.jl) is the single source of truth for every pack/unpack/view. Accessors: `temperature_range`, `concentration_range(l, j)`, `z_range`, `chi_range`, `z_offset`, `chi_offset`, `state_length`, plus the compile-time `block_view(u, l, Val{:T}|Val{:ξ}|Val{:z}|Val{:χ})`. `pack_state!` writes `mesh.cell_states` and (ALE) node positions into `u`; the χ block is initialized to zero because χ is reconstructed each residual rather than integrated as raw state (§10.3).

### 13.10.2 Geometry-from-state during assembly

In ALE mode the geometry is a function of the `z` block, so it must be recomputed from `u` *without mutating the mesh* (the residual and Jacobian must be pure functions of `u`). The geometry kernel reads node positions directly from `u` at the `z_offset`:

$$
V_i = |u_{z_{\text{off}}+i+1} - u_{z_{\text{off}}+i}|\,A, \quad
z_i^c = \tfrac12\,(u_{z_{\text{off}}+i} + u_{z_{\text{off}}+i+1}), \quad
z_f = u_{z_{\text{off}}+f},
$$

with face distances formed by each kernel from these centers and positions (§13.2). Because `A` and the `u`-derived positions carry the live scalar type, the entire geometry — and every flux and divergence built on it — is differentiable end-to-end.

> **Implementation note.** `update_state_geometry!(cache, mesh, ::ALE)` (src/Problem/state_cache.jl) recomputes `face_positions`, `A_face`, `cell_centers`, and `volume` from `cache.z` and `cache.A_section[]` once per RHS evaluation in ALE mode (the Eulerian method is a no-op: that geometry is mirrored once at workspace construction). Between *accepted* integrator steps, `update_mesh_from_state!` (src/Problem/state_layout.jl) mirrors the state back into `mesh.nodes[i].z` for the depletion checks and post-processing. The `StateCache` holds `T`, `ξ`, `z`, `χ`, the geometry quantities `volume`/`A_face`/`A_section`, `cell_centers`, and `face_positions`, all parameterized on the scalar type `Tu`.

### 13.10.3 RHS cache

The face fluxes, sources, and divergences are stored in a preallocated `RHSCache{Tu}` (src/Discretization/finite_volume.jl) carrying `q_cond`, `q_total`, `S_conv`, `J_gas` (`NC × n_faces`), `Q_rxn`, `Q_rad`, `dξ_rxn` (`NC × n_cells`), the reaction workspaces `rates`/`heats`/`species_src`, and the divergences `div_q`/`div_J`. The type parameter `Tu` is `Float64` for ordinary evaluation and `ForwardDiff.Dual` for AD-based Jacobian/sensitivity passes, so no reallocation occurs when switching modes. The cache keeps its size for the whole solve: depletion parking (§14.3) leaves ghost slots in place and the kernels loop over the live cells only.

---

## 13.11 Comparison to Gpyro and ThermaKin2Ds Discretizations

We map the reference codes' notation onto the unified nomenclature of §2 on first use, then contrast the discretization choices.

### 13.11.1 Gpyro (Lautenberger & Fernandez-Pello)

Gpyro is a 1D (optionally pseudo-2D) generalized condensed-phase pyrolysis solver using a **Patankar-style finite-volume** discretization (control volumes about cell-centered nodes). Mapping its symbols to ours: Gpyro's pre-exponential `Z → A_i`, permeability `K → κ`, kinematic viscosity `ν` (in `ṁ'' = −(K/ν)(∂P/∂z − gρ_g)`), conversion `α_conv`, and the overbarred mixture density `ρ̄`.

- **Conservation form vs. non-conservative.** Gpyro writes integral control-volume balances of `ρ̄Δz`, `ρ̄Y_iΔz`, gas mass, and energy on volumes that deform with the solid. Pyrolysis.jl uses the same cell-centered FVM divergence operator (§13.3) on the ALE mesh, with every component riding the mesh, but solves the **non-conservative** temperature equation, storing every component in `ρc_p` and carrying the gas-flux enthalpy via `S_conv` (§13.7.5). Both telescope to the same boundary-flux balance (§13.3.3).
- **Interface interpolation.** Patankar's hallmark is the **harmonic interface conductivity** for diffusion coefficients, justified by flux continuity across a control-volume face — identical in spirit to the harmonic-mean `k_f`, `λ_f`, `κ_f` of §13.7.1–§13.7.3. Pyrolysis.jl additionally distance-weights the harmonic mean for non-uniform meshes.
- **Convection scheme selection.** Patankar/Gpyro selects the face scheme by the cell **Péclet number** (central below ≈2, upwind/power-law above), a built-in switching. Pyrolysis.jl does not auto-switch on Péclet: diffusion is always the volume-fraction-gradient flux and gas advection is first-order upwind (Darcy–Fick; under ALE the gas rides the mesh and has no mesh-motion advection, §13.8.3). Choosing resolution to keep the cell Péclet number moderate is left to the user (see §18.7).
- **Linear solver.** Gpyro is fully implicit (backward Euler) with Patankar linearization solved by **TDMA (Thomas algorithm)** for the tridiagonal system. Pyrolysis.jl delegates time integration to stiff OrdinaryDiffEq integrators with a structured/sparse Jacobian (§15); a TDMA appears internally only in the P1 radiation solve (§8.7).
- **Pressure/Darcy.** Gpyro's Darcian gas flux `ṁ'' = −(K/ν)(∂P/∂z − gρ_g)` includes buoyancy; Pyrolysis.jl's `v_g = −(κ/μ)∇P` (§13.7.3) omits the gravity term (1D through-thickness, gravity negligible) and uses dynamic viscosity `μ` rather than kinematic `ν`.

### 13.11.2 ThermaKin2Ds (Stoliarov et al.)

ThermaKin is a 1D (with 2D extensions) finite-difference pyrolysis solver. Mapping its symbols: stoichiometric coefficient `θ_i^j → ν_{i,j}` (overload A2; ThermaKin's `θ` is stoichiometry, ours is strain rate), and its published heat-of-reaction sign is `h > 0 = exothermic` (overload H1; internally we use `h > 0` endothermic with `Q_rxn = −Σh_r r_r`).

- **Gas transport.** Pyrolysis.jl's volume-fraction-gradient flux `J = −(λ/T)∂(ξT)/∂x` (§13.7.2) is taken directly from ThermaKin's Section 1.2 gas-transport formulation `J = −ρ_g λ ∂(ξ_g/ρ_g)/∂x`. For a rigid non-swelling matrix (`γ = 0`) ThermaKin recovers a Boyle's-law Darcy form; Pyrolysis.jl provides the explicit Darcy–Fick path (§13.7.3) for that regime.
- **Convective energy term.** ThermaKin's Eq. 17 carries gas sensible enthalpy as `c_g (J_g·∂T/∂x)`, which is exactly the `S_conv` source of §13.7.5 (the S_conv-only convention), and sums every component in the storage term — both as here.
- **Time discretization.** ThermaKin uses Crank–Nicolson in the primary (through-thickness, `x`) direction and explicit updates in the transverse `(y, z)` directions; the spatial operators are second-order finite differences on a structured grid. Pyrolysis.jl is finite-**volume** (flux-based, conservative-divergence) rather than finite-difference, and offloads temporal integration to adaptive stiff solvers (§15).
- **Moving mesh.** ThermaKin performs dynamic mesh adjustment around a single reference point as the material deforms; Pyrolysis.jl recasts the same deformation through the ALE framework (§11) with cumulative mesh velocity built from per-cell dilation, and (optionally) a separate lateral-shrinkage `χ` channel (§10).

The shared lineage (harmonic interface averaging, S_conv gas-enthalpy convention, volume-fraction-gradient gas flux) is intentional — Pyrolysis.jl is designed to reproduce ThermaKin/Gpyro physics while adding type-stable, allocation-free, AD-compatible assembly and a structured Jacobian.

---

## 13.12 Limitations and Open Issues

The following discretization-specific issues are open (drawn from the source annotations and the verification against the implementation):

1. **ALE advection limiter is fixed.** `VANLEER` is the default `limiter` keyword of `compute_ale_advection_term!` (the mesh-relative temperature advection, override-only), and the residual calls it with that default; neither the workspace nor `solve` exposes another choice (minmod/superbee/none), which would aid verification studies.

2. **First-order Darcy–Fick advection only.** The pressure-driven advective face flux uses first-order upwinding (§13.7.3); no TVD/MUSCL option exists for Darcy advection (the TVD machinery is wired only into the verification-hook ALE temperature advection). Sharp pressure-driven composition fronts are therefore smeared at first order unless the mesh is refined.

3. **Multi-dimensional generalization.** All geometry, the divergence operator, the GCL, and the face structures are specialized to 1D (one face per node, scalar face areas, all faces sharing `A`). The `MeshTopology{D,...}` and `MeshGeometry{D}` types are parameterized by dimension `D` but are placeholders the 1D residual does not read; 2D/3D would require vector node positions, per-face normals and variable areas, and a full divergence-theorem GCL.

For the temporal integration and Jacobian that consume the operators of this chapter, see §15 (Temporal Integration & Structured Jacobian); for depletion parking and the refinement building blocks, see §14 (Adaptive Mesh Refinement and Surface Depletion); for the conservation audits that exploit the telescoping divergence of §13.3.3, see §16 (Conservation Diagnostics).


---

# 14. Adaptive Mesh Refinement and Surface Depletion

This chapter develops the adaptive-mesh-refinement (AMR) machinery of Pyrolysis.jl
and, in particular, the *surface-depletion* mechanism that is the only AMR path
exercised in the production solver. We treat the local error indicators that flag
where the mesh is under-resolved, the equidistribution theory behind
$r$-refinement (node relocation), the energy-conserving cell-merge algebra that
removes the $\rho c_p\!\to\!0$ singularity of a depleting surface cell, the
experimental $h$-refinement/$h$-coarsening (cell split/merge) topology operators,
the mesh-quality metrics that keep the discretization well-conditioned, the
active-cell bookkeeping that all AMR consumers must respect, and the way the
depletion callback parks a burned-out cell without changing the size of the
state. Every discrete rule, threshold, guard, and sign convention below is grounded
in the implementing function so that the mathematics is traceable to the code;
where the code makes an approximation, a guard, or a choice that differs from the
idealized form (including a temperature-interpolation term that vanishes
identically), we say so and trust the code.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$ (axial, through-thickness). **$z=0$ is the
  bottom/substrate** (boundary id 2, tag `:bottom`); **$z=L$ is the top/exposed
  surface** (boundary id 1, tag `:top`). Heat enters from the top; the material
  recedes downward. Cell indices increase from bottom to top, so the "left" cell
  of a face is at lower $z$ and the "right" cell is at higher $z$.
- The cell index is $i$ ($i=1,\dots,n$); $j$ is a component index
  ($j=1,\dots,N_C$); $r$ a reaction index; $f$ a face index. This follows
  overload resolution I1 of §2.
- Following the adaptivity subsystem we write the **cell size** $h_i\equiv\Delta z_i$
  for the through-thickness width of cell $i$, to match the AMR literature; this
  is the same quantity called $\Delta z_i$ elsewhere in the guide. We reserve
  $\eta_i$ for an error indicator value and $\omega(z)$ for a monitor function, as
  in §2.
- $\phi$ in this chapter is a **generic scalar field** (temperature $T$ or a
  concentration $\xi_j$) being tracked by an indicator or a monitor, *not* the
  porosity of §5. The porosity does not appear in this chapter.

The subsystem splits along a maturity line. The **core** `Adaptivity` module
(`src/Adaptivity/`) holds the depleted-cell park (`park_cell!`) with its
energy-conserving temperature merge and neighbor finders — the only adaptivity
`solve` uses — and building blocks that `solve` does not call: the error
indicators, $r$-refinement, the mesh-quality metrics, and the invariant checker.
The **experimental** $h$-AMR cluster (`src/Experimental/`) holds cell splitting
(`refine_cell!`), sibling coarsening (`coarsen_cells!`), the `AMRController`
orchestrator, and the combined `hr_adapt!` driver. It is not wired into `solve`,
and a mesh it refines cannot be passed to `solve`: refinement leaves inactive
parent cells among the active ones, and the solver's kernels assume the active
cells come first (`src/Experimental/Experimental.jl`). This chapter covers both,
flagging the status of each throughout.

---

## 14.1 Error indicators

An *error indicator* $\eta_i\ge 0$ assigns to each active cell a scalar estimate
of the local discretization error, used to decide where to refine (large $\eta$)
or coarsen (small $\eta$). Pyrolysis.jl implements four physically motivated
indicators plus a composite, all subtypes of the abstract type
`AbstractErrorIndicator` (declared in `src/Materials/Materials.jl`, so that the
type is visible to both the `Adaptivity` and `Experimental` modules). All
indicator math lives in `src/Adaptivity/error_indicators.jl`.

The unifying design is **dimensional consistency through the cell size**: an
$O(h^p)$ truncation-error term is mimicked by multiplying a discrete derivative of
the field by the appropriate power of $h_i$. This makes $\eta_i$ behave like a
local error contribution, so that thresholding $\eta_i$ against the global maximum
$\max_j\eta_j$ adapts the mesh toward equidistributed error.

### 14.1.1 Gradient indicator

The gradient indicator targets steep solution fronts (thermal waves, reaction
zones). For a tracked field $\phi$,

$$
\eta_i^{\mathrm{grad}} = w\, h_i \,\bigl|\nabla\phi\bigr|_i ,
$$

where $w$ is a user weight (default $1$), $h_i$ is the cell size, and
$|\nabla\phi|_i$ is the magnitude of the discrete first derivative at cell $i$.
The factor $h_i$ makes $\eta^{\mathrm{grad}}$ scale like the first-order term in a
Taylor expansion: a fixed jump across a cell produces a fixed indicator regardless
of how the cell is sized, so the indicator flags *resolution relative to the
feature*, not absolute slope.

The gradient itself (`compute_gradient_1d`) is the
opposite-distance-weighted mean of the two one-sided slopes
(`opposite_distance_gradient`, §13.6.1) — second order on non-uniform
meshes — and degrades gracefully at boundaries:

$$
(\nabla\phi)_i =
\begin{cases}
\dfrac{d_R\,\frac{\phi_i-\phi_L}{d_L} + d_L\,\frac{\phi_R-\phi_i}{d_R}}{d_L+d_R}, & \text{both neighbors present},\\[2.2ex]
\dfrac{\phi_R-\phi_i}{z^c_R-z^c_i}, & \text{left is a boundary (one-sided right)},\\[2.2ex]
\dfrac{\phi_i-\phi_L}{z^c_i-z^c_L}, & \text{right is a boundary (one-sided left)},\\[2.2ex]
0, & \text{a lone cell (no neighbors)} .
\end{cases}
$$

Here $L,R$ index the cells immediately below/above cell $i$, found through the
cell's left/right faces, $z^c$ are cell-center positions, and
$d_L = z^c_i - z^c_L$, $d_R = z^c_R - z^c_i$ are the center-to-center
distances. A boundary cell takes the one-sided difference toward its only
neighbor, so it never reads a missing value, and a zero center-to-center distance
returns a zero gradient.

### 14.1.2 Curvature indicator

The curvature indicator is *more* sensitive to sharp features than the gradient
indicator, because it carries two powers of the cell size:

$$
\eta_i^{\mathrm{curv}} = w\, h_i^2 \,\bigl|\nabla^2\phi\bigr|_i .
$$

The discrete second derivative (`compute_curvature_1d`) uses the
three-point **non-uniform mesh** stencil:

$$
(\nabla^2\phi)_i =
\frac{2}{h_L+h_R}
\left[
  \frac{\phi_R-\phi_i}{h_R}
  -
  \frac{\phi_i-\phi_L}{h_L}
\right],
\qquad
h_L=z^c_i-z^c_L,\quad h_R=z^c_R-z^c_i .
$$

For $h_L=h_R=h$, this reduces to the standard centered formula
$(\phi_R-2\phi_i+\phi_L)/h^2$. On graded meshes it preserves zero curvature for
linear fields, avoiding false refinement from mesh stretching alone.

Because three collinear points are required, the curvature is returned as
**exactly zero at any cell that lacks an interior neighbor on either side**
(boundary cells); the indicator therefore never drives refinement of the first or
last cell. This is a deliberate restriction: a one-sided second difference at a
boundary is ill-conditioned and would falsely flag the surface cell, which is
already handled by the depletion mechanism (§14.3).

> **Implementation note.** Gradient and curvature indicators are
> `compute_indicator(::GradientIndicator, …)` /
> `compute_indicator(::CurvatureIndicator, …)`, with the derivative kernels
> `compute_gradient_1d` and `compute_curvature_1d`, all in
> `src/Adaptivity/error_indicators.jl`. Both skip cells with
> `active_cell_mask[i] == false`, returning $\eta_i=0$ there. The variable to
> track is the `variable::Union{Symbol,Int}` field — `:T` for temperature, or an
> integer component index `j` to track $\xi_{j,\cdot}$; any other symbol throws
> (component names are not resolvable without a material in scope).

### 14.1.3 Reaction-rate indicator

Conduction- and diffusion-based indicators cannot see a reaction front that has
not yet imprinted a steep temperature or concentration gradient. The
reaction-rate indicator refines wherever chemistry is locally fast relative to the
available reactant:

$$
\eta_i^{\mathrm{react}} = w\, h_i \,\max_j \frac{|r_j|_i}{\max(\xi_{j,i},\,\varepsilon)},
\qquad \varepsilon = 10^{-10}\ \mathrm{kg\,m^{-3}} .
$$

The normalization by $\xi_{j,i}$ converts an absolute rate (kg m$^{-3}$ s$^{-1}$)
into an inverse time scale (s$^{-1}$), i.e. a *fractional* consumption rate; the
$\max_j$ takes the worst-case over the rate matrix rows, so a single fast species
suffices to flag the cell. The floor $\varepsilon$ prevents division by zero in a
depleted cell.

This indicator is special in that it **requires** the per-cell, per-species
reaction rate matrix `rates`, passed as a keyword:
`compute_indicator(ind, mesh, T, ξ; rates = r)` with `r` sized like ξ for the
*current* mesh. Calling it without rates throws an `ArgumentError`, and a `rates` matrix whose column count does not match
the mesh throws `DimensionMismatch`. Because `adapt_mesh!` recomputes
indicators *after* topology changes — when any precomputed rates are stale —
this indicator cannot run inside the controller cycle and is standalone-only;
the default controller uses a gradient + curvature composite instead.
Supplying live rates inside the time-stepping loop would mean recomputing the
Arrhenius source at the indicator-evaluation point (see §6, Reaction
Kinetics) — this coupling gap is an open issue (§14.10).

> **Implementation note.** `compute_indicator(::ReactionRateIndicator, mesh, T, ξ;
> rates=nothing)` in `src/Adaptivity/error_indicators.jl`.

### 14.1.4 Interface-proximity indicator

The previous three indicators are *solution-driven*. The interface-proximity
indicator is purely *geometric*: it forces refinement in a neighborhood of a
tracked interface position $z_{\mathrm{int}}$ regardless of the current solution,
which is useful for a moving pyrolysis front that the solution has not yet
resolved:

$$
\eta_i^{\mathrm{iface}} = \frac{w}{1 + d_i/h_{\mathrm{target}}},
\qquad
d_i = |z^c_i - z_{\mathrm{int}}| ,
$$

with $h_{\mathrm{target}}$ the desired cell size near the interface. The indicator
is a Lorentzian bump: it equals $w$ exactly at the interface ($d_i=0$), falls to
$w/2$ one target-size away, and decays as $w\,h_{\mathrm{target}}/d_i$ far from it.
Because $\eta^{\mathrm{iface}}$ is bounded and smooth, it composes cleanly with the
solution-driven indicators.

> **Implementation note.** `compute_indicator(::InterfaceProximityIndicator, …)`
> in `src/Adaptivity/error_indicators.jl`. The struct carries
> `interface_position`, `h_target`, and `weight`; no derivative is computed.

### 14.1.5 Composite indicator

Physical problems usually need several indicators at once (e.g. a thermal-gradient
indicator together with a curvature indicator). `CompositeIndicator` evaluates
each component indicator $\eta_k(i)$, **normalizes each by its own maximum
over active cells**, $\hat\eta_k(i) = \eta_k(i)/\max_{j\in\text{active}}
\eta_k(j)$, and combines the normalized values per-cell by one of three rules:

$$
\eta_i^{\max} = \max_k \hat\eta_k(i),
\qquad
\eta_i^{\ell_2} = \sqrt{\sum_k \hat\eta_k(i)^2},
\qquad
\eta_i^{\Sigma} = \sum_k \hat\eta_k(i) .
$$

The normalization is what makes the combination meaningful: component
indicators carry incommensurable dimensions (a temperature-gradient indicator
is in K, a reaction-rate indicator in s⁻¹), so combining raw values would let
whichever indicator has the largest unit-system numbers dominate. The $\max$
rule (the default) is the most conservative: a cell is flagged if *any*
component demands it. The $\ell_2$ and $\sum$ rules blend contributions,
useful when several weak indicators should reinforce one another. The
per-component weights $w_k$ are baked into each $\eta_k$ *before*
normalization, so they tilt the within-indicator shape, while the relative
scale across indicators is set by the normalization.

> **Implementation note.** `compute_indicator(::CompositeIndicator, …)` in
> `src/Adaptivity/error_indicators.jl` dispatches on the
> `combination::Symbol` field (`:max`, `:l2`, `:sum`) and errors on any other
> value. A convenience overload `compute_indicator(ind, mesh)` pulls $T$ and $\xi$ from the mesh cell states via `extract_fields` and forwards
> to the field-taking method; this dense extraction happens on every call and is
> not cached.

### 14.1.6 Indicator interface and the active-cell contract

Every indicator computes over `n_cells(mesh)` cells but **writes $\eta_i=0$ for any
inactive cell** (`active_cell_mask[i] == false`). This is the first instance of a
pervasive AMR contract: because $h$-refinement appends children to `mesh.cells`
without ever compacting (§14.6), `length(mesh.cells)` can exceed
`n_active_cells`, and *every* consumer must guard on the mask rather than iterate a
raw index range. The indicators honor this; downstream marking routines
(§14.4.2) compute $\max\eta$ only over masked cells.

---

## 14.2 $r$-refinement: node relocation by equidistribution

$r$-refinement ("relocation refinement") moves the interior nodes of a *fixed
topology* mesh to concentrate resolution where it is needed, without changing the
number of cells. Unlike $h$-refinement it adds no degrees of freedom — it is a
pure remeshing of the existing cells. `solve` does not call it: the ALE mesh moves
with the solid (§11.1.2), and node relocation decoupled from the solid is not
representable there (§11.9 item 5). The governing principle is **equidistribution** of a
user-chosen *monitor function* $\omega(z)>0$. All $r$-refinement code is in
`src/Adaptivity/r_refinement.jl`.

### 14.2.1 The equidistribution principle

Given a monitor $\omega(z)$ that is large where mesh resolution is desired, the
equidistribution principle places node $z_n^\ast$ so that each subinterval carries
an equal share of the total "monitor mass":

$$
\int_{0}^{z_n^\ast}\!\omega(z)\,dz
= \frac{n}{N+1}\int_{0}^{L}\!\omega(z)\,dz,
\qquad n = 1,\dots,N,
$$

where $N$ is the number of interior nodes (so there are $N+1$ cells). Equivalently,
each cell $[z_n^\ast,z_{n+1}^\ast]$ contains the same $\int\omega\,dz$. Where
$\omega$ is large the cells must be small to keep the integral constant — exactly
the desired refinement. The continuous monitor integral is discretized as a sum
over cells,

$$
\Omega_{\mathrm{total}} = \int_0^L\omega\,dz \approx \sum_i \omega_i\,h_i ,
$$

with $\omega_i$ the monitor evaluated at cell center $z^c_i$ and $h_i$ the cell
size.

### 14.2.2 Monitor functions

**Gradient monitor.** The most common monitor is built from the solution
gradient, regularized so it never vanishes:

$$
\omega^{\mathrm{grad}}(z) = \sqrt{1 + \alpha\,|\nabla\phi|^2}\, .
$$

The "$1+$" floor guarantees $\omega\ge 1$ everywhere (a uniform mesh in the
gradient-free limit), while $\alpha$ tunes the sensitivity to the gradient. This
is the classical *arc-length monitor* of moving-mesh PDE methods: in 1D,
$\int\sqrt{1+\alpha\,\phi_z^2}\,dz$ is (for $\alpha=1$) the arc length of the
solution curve, so equidistributing $\omega^{\mathrm{grad}}$ places nodes at equal
arc-length intervals along the solution — naturally clustering them on steep
fronts. The gradient is the same distance-weighted center difference used by the
indicators (`compute_gradient_1d`, `src/Adaptivity/error_indicators.jl`), with
one-sided forms at boundaries.

**Interface monitor.** A Gaussian bump concentrates nodes near a prescribed
interface $z_{\mathrm{int}}$ irrespective of the solution:

$$
\omega^{\mathrm{iface}}(z) = \omega_{\mathrm{base}}
+ \omega_{\mathrm{peak}}\,\exp\!\left(-\frac{(z-z_{\mathrm{int}})^2}{2\sigma^2}\right).
$$

The defaults are $\sigma=10^{-3}\,$m, $\omega_{\mathrm{base}}=1$,
$\omega_{\mathrm{peak}}=10$, so the monitor is ten-fold larger at the interface
than in the far field, with a Gaussian width $\sigma$. Here $\sigma$ is the
interface-monitor width $\sigma_{\mathrm{mon}}$ of §2.

**Composite monitor.** `CompositeMonitor` combines several monitors by `:max`,
`:sum`, or `:product`. The product rule is multiplicative reinforcement (both
monitors must be large); the max rule is the conservative union.

> **Implementation note.** `evaluate_monitor` (`src/Adaptivity/r_refinement.jl`)
> dispatches on `GradientMonitor`, `InterfaceMonitor`, and
> `CompositeMonitor`; monitors are subtypes of the abstract `MonitorFunction`.
> Like the indicators, monitor evaluation skips inactive cells. Note that monitor
> functions are problem-dependent and carry arbitrary units; the user is
> responsible for keeping the disparity in $\omega$ bounded, because a large
> dynamic range in $\omega$ makes the equidistribution solve ill-conditioned.

### 14.2.3 Computing target node positions

`compute_target_positions(mesh, ω)` inverts the discrete
equidistribution relation by piecewise-linear interpolation of the cumulative
monitor integral. The algorithm:

1. **Sort active cells by $z$**, building the list
   $\{(i, z^c_i, \omega_i, h_i)\}$ ordered by center position. This step is what
   makes $r$-refinement robust to the non-monotone cell ordering that AMR can
   leave behind.
2. **Form the cumulative monitor mass.** With $C_0=0$ and
   $C_k = C_{k-1} + \omega_k h_k$ for $k=1,\dots,M$ over the $M$ sorted cells,
   $C_M = \Omega_{\mathrm{total}}$. A companion array $z_k$ holds running cell
   right-boundary estimates, with $z_0 = z_{\min}$ and $z_M = z_{\max}$.
3. **Invert for each interior node.** For interior node $n=1,\dots,N$ (with
   $N$ the number of *active* interior nodes — nodes referenced by an active
   cell, minus the two boundaries; $h$-coarsening orphans are excluded via
   `_active_node_mask`), set the target cumulative value
   $\Omega_n^\ast = \tfrac{n}{N+1}\,\Omega_{\mathrm{total}}$, locate the cumulative
   interval $[C_{k},C_{k+1}]$ containing it by `searchsortedlast`, and
   linearly interpolate within it:

$$
z_n^\ast = z_k + \frac{\Omega_n^\ast - C_k}{C_{k+1}-C_k}\,(z_{k+1}-z_k).
$$

   If a cumulative interval is degenerate ($C_{k+1}=C_k$, a zero-monitor cell) the
   fraction defaults to $\tfrac12$.

The routine returns the empty vector if there are no interior nodes or if
$\Omega_{\mathrm{total}}\le 0$ (which then disables relocation for that step).

The inversion is *node-count* based ($\tfrac{n}{N+1}\Omega_{\mathrm{total}}$), not
*cell-count* based, which is the correct choice for placing the $N$ interior nodes.

### 14.2.4 Relaxed relocation with quality constraints

`relocate_nodes!` drives the interior nodes toward their targets
by damped fixed-point iteration. Per iteration it:

1. Evaluates the monitor and computes target positions.
2. For each *active* interior node $n$ (sorted by position; nodes orphaned by
   $h$-coarsening are skipped via `_active_node_mask`) takes a **relaxed step**

$$
z_n^{\mathrm{new}} = z_n^{\mathrm{old}} + \beta\,(z_n^\ast - z_n^{\mathrm{old}}),
\qquad \beta\in(0,1],
$$

   with default relaxation $\beta=0.5$. Relaxation prevents node overshoot and
   tangling when the monitor changes sharply between iterations.

3. **Enforces a minimum cell size** by clamping the new position against its
   already-updated neighbors:

$$
z_{n-1}^{\mathrm{new}} + h_{\min} \;\le\; z_n^{\mathrm{new}} \;\le\; z_{n+1}^{\mathrm{old}} - h_{\min},
$$

   with the domain boundaries $z_{\min},z_{\max}$ standing in for absent neighbors
   at the ends. The default $h_{\min}=10^{-8}\,$m. This clamp can prevent a node
   from reaching its target if the equidistribution would demand a cell thinner
   than $h_{\min}$.

4. Calls `update_geometry!(mesh)` to refresh centers, volumes, and face positions.

5. **Checks convergence** on the maximum *relative* node movement,

$$
\max_n \frac{|z_n^{\mathrm{new}} - z_n^{\mathrm{old}}|}{L} < \mathrm{tol},
$$

   default $\mathrm{tol}=10^{-3}$. The function returns the iteration count at
   convergence, or `max_iter` (default $10$) if not converged.

**Boundary nodes never move** — only interior nodes are relocated. The boundary
positions are physical domain limits and are fixed by construction.

A critical caveat: `relocate_nodes!` relocates nodes but does **not remap the
solution onto the moved mesh**. It reads the cell states once at entry
(`extract_fields`) to evaluate the monitor and then moves nodes; the cell
*state* (the $T$ and $\xi$ stored per cell) is left attached to its resized
cell, so a cell's mass $\xi_j V_i$ changes with its volume. The package has no
conservative remap; `relocate_nodes!` is a *geometry* operation, conservation of
the transported fields is the caller's responsibility, and `solve` never calls
it.

---

## 14.3 Depletion handling: parked cells (the production path)

The single AMR capability wired into `solve` is **depletion handling**, enabled by
`solve(…; handle_depletion = true)` (default `false`). As
pyrolysis consumes the solid, an ALE cell (§11) thins toward zero thickness; on a
fixed mesh a cell empties instead. A vanishing or emptied cell is numerically
catastrophic for two reasons. First, its matrix heat capacity
$\rho c_p^{\mathrm{eff}} V_i \to 0$ as the solid gasifies, so the energy equation
$\rho c_p^{\mathrm{eff}}\,\dot T_i = (\text{fluxes})$ becomes singular. Second, the
transport rows of a thin cell scale as $1/\Delta z_i^2$, so the condition number
of the implicit stage matrix $W = I/(\gamma\Delta t) - J$ grows as
$\Delta z_i^{-2}$: at $\Delta z \approx 10^{-7}$ m, $\operatorname{cond}(W)$
exceeds $10^{13}$ for any $\Delta t \gtrsim 10^{-3}$ s and Newton's linear solves
carry no usable digits, so the integrator halves $\Delta t$ until the cell's own
time scale is resolved. No Jacobian accuracy and no regularisation of the
properties changes this; the cell must leave the nonlinear system.

Depletion handling **parks** the cell: it merges the cell's mass and sensible
enthalpy into a neighbour, then moves the freed state slot above the live
surface as an inert ghost. The state vector, the mesh arrays, the
`StateLayout`, the sparsity pattern and the Jacobian caches keep their size for
the whole solve; every kernel runs over the live range, and the boundary
conditions act at the live surface.

### 14.3.1 Live range and ghost slots

With $n$ the (fixed) number of state slots and $n_{\rm live}$ =
`mesh.n_active_cells`:

- live cells are `1:n_live`, live nodes and faces `1:n_live+1`; face
  `n_live+1` is the exposed boundary face (`Face1D` with `cell_ids =
  (n_live, 0)`, `boundary_id = 1`; `mesh.boundary_faces[1] == [n_live+1]`) and
  cell `n_live` the boundary cell;
- ghost slots are `n_live+1:n` (cells) and `n_live+2:n+1` (nodes);
  `active_cell_mask[i] == (i ≤ n_live)`;
- a ghost carries $\xi_j = 0$ for every component, the temperature it was
  parked with, $\chi = 0$ and `m_dry_solid_initial = 0`. Under ALE its nodes
  sit on the live surface (zero thickness) and move with it,
  $\dot z_k = w_{n_{\rm live}+1}$ for $k > n_{\rm live}+1$; on a fixed mesh
  they keep their geometry;
- every residual kernel and every Jacobian operator loops over
  `1:mesh.n_active_cells` (cells) or `1:mesh.n_active_cells+1` (faces) and
  detects the top boundary through the face connectivity, so no kernel reads a
  ghost. The pack step writes exactly zero into the ghost $T$, $\xi$ and $\chi$
  rows. Their Jacobian rows are therefore zero and $W$ is block-diagonal with
  $I/(\gamma\Delta t)$ in the ghost block — perfectly conditioned — while the
  error estimate of a ghost component is identically zero;
- the live thickness is $z_{n_{\rm live}+1} - z_1$; every ledger, mass and
  energy sum runs over the live cells.

The correctness statement is an **equivalence**: the residual (and the
Jacobian) of the live rows of a parked state equal, to roundoff, those of the
same physical configuration built as an $n_{\rm live}$-cell problem
(`test/residual/test_park_equivalence.jl`, `test/jacobian/test_park_jacobian.jl`,
for the Eulerian, ALE, Darcy, Beer–Lambert and variable-cross-section modes).

### 14.3.2 Trigger

Under ALE a live cell is depleted when its current thickness falls below a
fixed fraction of its *initial* thickness,

$$
\Delta z_i = z_{i+1}-z_i \;<\; \tau_{\mathrm{dep}}\,\Delta z_{i,\mathrm{init}},
$$

with the default `depletion_threshold` $\tau_{\mathrm{dep}} = 0.05$. Node
positions are read from the ODE state, so detection sees the live integrator
state. The reference thicknesses are stored per state slot in
`DepletionCallbackState.reference` (computed once as $V_{i,\mathrm{init}}/A_0$),
so each cell is judged against its *own* initial thickness — the only detector
that is correct on a stretched mesh, where a single shared reference would, for
any cell-size ratio above $1/\tau_{\mathrm{dep}} = 20$, flag the fine cells at
$t = 0$ or never flag the coarse ones. When an interior cell parks, the
reference entries shift with the cells. `find_thin_cells` returns the thin live
cells surface first (descending $z$) and throws on a non-positive signed
thickness: an inverted cell is a mesh-motion failure, not depletion.

The threshold is not a lever for Newton failures. Raising $\tau_{\rm dep}$
merges the hottest cell into a colder one earlier, which moves the surface
temperature: on the Landau case V-CPL-02 (2 mm/s on 62.5 µm cells), going
from 0.05 to 0.3 shifts the Landau identity residual from $7\times10^{-4}$ to
$2\times10^{-3}$. The Newton failures on a fast ablation front have two
other causes, degenerate gas transport (§5.6) and a step controller that
regrows the step beyond the Newton's range (§15.5.3), and are removed there.
The default stays 0.05.

On a fixed mesh cells cannot thin; they empty. The trigger is the **top live
cell's condensed concentration**,
$\sum_{j\notin\mathrm{gas}}\xi_{j,n_{\rm live}} < \tau_{\mathrm{dep}}\,
\sum_{j\notin\mathrm{gas}}\xi_{j,n_{\rm live}}^{\mathrm{init}}$, with the same
threshold (a cell that starts with no condensed mass never parks). Only the
top live cell parks, so a fixed mesh burns a layer top down: the surface cell
empties and parks, and the cell beneath it becomes the surface. That requires
the surface cell to be the hottest, as under surface absorption. On the PMMA
burnout of V-BO-02 all 60 cells park one by one from 63 to 452 s, and the
maximum live-cell temperature stays at 915 K.

Under in-depth absorption the temperature maximum of a regressing
non-charring polymer lies below the surface, so the surface cell reacts
slowest. Under ALE it shrinks and parks, and the surface advances into the
hotter material. On a fixed mesh the surface cell keeps its volume and loses
concentration: it persists as a low-density skin while the cells beneath it
empty, and the layer gasifies from within. A void cell conducts with the gas
conductivity and absorbs nothing, and nothing in the model carries heat across
a void by radiation. The flux transmitted through the void column is deposited
in the first condensed cell below it, which can return that heat upward only
by conduction through the gas, and the temperatures run away. Continued past its first void, the PMMA
deck heats the top substrate cell past 1000 K within about 1.5 s and past
3000 K, the end of the property bracket, about 60 s later. From the first
interior void on, the model's equations have no solution with valid
properties.

A fixed-mesh solve therefore stops at its first interior void. After every
accepted step the depletion callback tests each live cell below the top live
cell for a matrix volume fraction
$\varphi_m = \sum_j \gamma_j\xi_j/\rho_j < 10^{-4}$ (`MIXING_BLEND_LO`), the
property model's own definition of a void (§5.1); a cell that starts with no
condensed mass is not tested. On a hit the callback warns that the material
has gasified from within, that the model carries no heat across a void and
that a regressing material belongs on the ALE mesh, and terminates the
integration (retcode `:Terminated`). The in-depth PMMA deck stops at
337.0 s, when the lowest PMMA cell becomes a void, with the maximum live-cell
temperature at 912 K (V-BO-05). A fixed
mesh with depletion is the porous-residue model (a char that keeps its
volume, burning top down); a material that regresses belongs on the ALE
mesh.

> **Implementation note.** `is_cell_thin`, `find_thin_cells`,
> `depletion_candidate` and `interior_void_cell` in
> `src/Solver/callbacks.jl`; `matrix_volume_fraction` in
> `src/Properties/mixture_rules.jl`. The callback *condition* fires whenever a
> candidate or (fixed mesh) an interior void exists; the *affect* terminates
> the integration on a void, parks when `n_live > min_cells` (default
> `min_cells = 2`) and **terminates the integration** otherwise.

### 14.3.3 Energy-conserving temperature merge

When cells $a,b$ of volumes $V_a,V_b$ are merged into one, the concentrations
combine by **volume-weighted averaging**, which conserves species mass exactly
because mass is the extensive product $\xi_j V$:

$$
\xi_{\mathrm{merged},j} = \frac{V_a\,\xi_{a,j} + V_b\,\xi_{b,j}}{V_{\mathrm{keep}}},
$$

with $V_{\mathrm{keep}} = V_a + V_b$ under ALE (the keeper spans both cells) and
$V_{\mathrm{keep}} = V_b$ on a fixed mesh (the keeper's volume is fixed; the
parked cell's mass raises its concentration). Either way
$\xi_{\mathrm{merged},j}V_{\mathrm{keep}} = \xi_{a,j}V_a + \xi_{b,j}V_b$.

Temperature is **intensive**, and a naive volume-weighted average of $T$ does
*not* conserve internal energy when the two cells have different heat capacities —
the very situation in depletion, where a hot, dense solid cell is merged with a
cool, mostly-converted (low-capacity) one. Pyrolysis.jl therefore conserves the
**sensible enthalpy** of every component

$$
H(T) \;=\; \sum_{j=1}^{N_C} m_j \int_{T_{\mathrm{ref}}}^{T}
c_{p,j}(T')\,\mathrm{d}T',
\qquad m_j = V\,\xi_j ,
$$

and `T_merged` solves $H_{\mathrm{merged}}(T_{\mathrm{merged}}) = H_a(T_a) +
H_b(T_b)$. Two choices are deliberate:

- **$\int c_p\,\mathrm{d}T$, not $c_p(T)\cdot T$.** For temperature-dependent
  $c_p$ the product form is not an enthalpy: with $c_p = a + bT$ it carries the
  $b$-term with coefficient 1 instead of $\tfrac12$, a few-K systematic bias
  per merge (always the same sign when a hot surface cell merges into a cooler
  interior one). The integral is evaluated per component by 4-point
  Gauss–Legendre quadrature — exact for polynomial $c_p$ up to degree 7, and
  park events are rare so the cost is irrelevant. The datum
  $T_{\mathrm{ref}}$ cancels because each component's total mass is identical
  before and after the merge.
- **Every component is included.** The merged enthalpy sums the gas with the
  condensed components, the same $\rho c_p$ the temperature equation stores
  (§5, §7).

The solve is a bracketed Newton iteration on
$[\min(T_a,T_b),\,\max(T_a,T_b)]$: $H$ is strictly monotone ($c_p > 0$), and
$H_{\mathrm{merged}}(T_{\min}) \le H_{\mathrm{total}} \le
H_{\mathrm{merged}}(T_{\max})$ guarantees the root lies inside the bracket, so
the iteration cannot fail. If neither cell carries any heat capacity (both
empty) the energy balance is degenerate and the routine falls back to the
volume-weighted average.

The accuracy stakes: merging a $1\,$m$^3$ condensed-rich cell at $800\,$K with
$\rho c_p = 2000\,$kJ m$^{-3}$K$^{-1}$ and a $1\,$m$^3$
mostly-converted cell at $400\,$K with
$\rho c_p = 500\,$kJ m$^{-3}$K$^{-1}$ gives, by the wrong volume-weighted rule,
$T=(800+400)/2 = 600\,$K, whereas conserving the enthalpy (constant
$c_p$, where the integral form reduces to the familiar
$\rho c_p$-weighted mean) gives

$$
T_{\mathrm{merged}} = \frac{2000\cdot 800 + 500\cdot 400}{2000 + 500} = 720\ \mathrm{K},
$$

a $120\,$K (20 % of the temperature difference) error avoided.

> **Implementation note.** `merge_temperatures_conservatively(V_a, state_a, V_b,
> state_b, material)` in `src/Adaptivity/cell_merge.jl`. The component heat
> capacity is `material.components[j].heat_capacity`, a callable $c_{p,j}(T)$.
> This is the **single** canonical merge: the park primitive `park_cell!` and
> the experimental `coarsen_cells!` both call it.

### 14.3.4 The park event

`park_cell!(mesh, i_park, i_keep; material, ale)` (`src/Adaptivity/cell_park.jl`)
performs the mesh side of a park; `park!` (`src/Solver/callbacks.jl`) drives it
and rewrites the state vector. The sequence:

1. **Neighbour selection.** `find_merge_neighbor(mesh, i_park)`
   (`src/Adaptivity/surface_depletion.jl`): a boundary cell takes its only
   neighbour; an interior cell takes the **interior** (lower-$z$) one. The
   substrate side has the more slowly varying thermal state, so absorbing the
   thin cell into the cooler, larger interior cell produces the less abrupt
   merged temperature. Neighbours are found through the cell's faces
   (`find_interior_neighbor` / `find_exterior_neighbor`; either is 0 at a
   domain boundary, which is how the live top cell reports no cell above it).
2. **Conservative merge** of the two `CellState`s per §14.3.3, with
   $V_{\mathrm{keep}}$ as defined there; `m_dry_solid_initial` is summed into the
   keeper so the per-cell pyrolysis-progress denominator stays consistent
   (§10). The freed slot becomes the ghost state ($\xi = 0$, $T = T_{\rm merged}$,
   zero volume under ALE, its own volume on a fixed mesh, `m_dry_solid_initial = 0`).
3. **ALE slot permutation.** Cells `i_park+1:n_live` move down one slot and
   the nodes above the shared node `max(i_park, i_keep)` move down one slot
   (per-cell and per-node *values*: `cell_states`, `cell_volumes`,
   `m_dry_solid_initial`, node positions and velocities). The 1-D connectivity
   is regular — cell $i$ always owns nodes and faces $(i, i+1)$ — so no cell,
   face or node reference changes; the freed cell slot `n_live` (the last live
   index) becomes the ghost and the freed node slot lands on the surface. On a
   fixed mesh only the top live cell parks and nothing shifts.
4. **Boundary flags.** Face `n_live + 1` (`n_live` the new count) becomes the exposed
   boundary face, `cell_ids = (n_live, 0)`, `boundary_id = 1`; the previous top
   face becomes an inert ghost face; the matching node flags follow;
   `mesh.boundary_faces[1]` is updated; `n_active_cells` is decremented and
   `active_cell_mask` cleared for the ghost; the mesh geometry is recomputed
   (`update_geometry!`) and the AMR invariants (§14.6) asserted.
5. **State vector.** $T$, $\xi$ and (ALE) $z$ are read back from the mesh for
   every slot — the mesh mirrors the accepted state when the callback runs, so
   live cells shift exactly as the mesh did and the ghost carries the merged
   temperature, zero concentrations and the surface position. The $\chi$ block
   is shifted the same way with
   $\chi_{\rm merged} = (m_a\chi_a + m_b\chi_b)/(m_a+m_b)$ weighted by the
   **initial dry-solid masses** (captured before step 2 sums them), because
   $\chi_i = m_{\mathrm{released},i}/m_{\mathrm{dry,init},i}$ makes
   $\sum_i \chi_i\,m_{\mathrm{dry,init},i}$ the extensive invariant; the ghost
   gets $\chi = 0$. The workspace's cached top boundary face id follows the
   mesh (it is shared with the Jacobian's AD workspace, so one update serves
   both), and `derivative_discontinuity!(integrator, true)` forces a derivative
   recomputation.

Mass, sensible enthalpy, $\chi$-weighted dry mass and (under ALE) live volume are
conserved to roundoff at every event (`test/adaptivity/test_cell_park.jl`).

### 14.3.5 The depletion callback

The live machinery is a `DiscreteCallback` checked at the end of every
accepted integrator step (`create_depletion_callback`,
`src/Solver/callbacks.jl`). Its condition is `depletion_candidate(...) != 0`,
or on a fixed mesh an interior void (`interior_void_cell(...) != 0`); its
`affect!` terminates the integration with a warning on an interior void
(§14.3.2), and otherwise parks **exactly one cell per step**, surface first,
so the integrator adapts its step between events, or terminates the
integration when only `min_cells` live cells remain. It records every event in the shared
`DepletionCallbackState` (`park_times`, `parked_cells`, `n_live`), which the
thickness-termination callback (live thickness) and solution extraction also
read.

After the park the `affect!` flags the state as modified and proposes
`DT_AFTER_PARK` ($10^{-6}$ s) for the next step. The merged keeper cell
conserves mass and enthalpy but sits off the quasi-steady gas-concentration and
temperature profiles of its neighbours, so the residual at the parked state is
orders of magnitude above that of the accepted pre-park state, and the
integrator's extrapolated stage predictor, built from the pre-park derivative,
starts far from the post-park solution. At the pre-park step size the Newton
fails until the controller has cut the step far enough; from the proposed step
it converges at once, and the controller grows the step back within a few
steps. The reset is a cost setting, not a condition for completing the run: on
the 60-park PMMA burnout (V-BO-01) the run takes 858 steps and 8 Newton
failures with it and 1686 steps and 649 failures without it, with the same
mass-loss history. A cascade (several cells reaching the threshold together)
proceeds one park per $10^{-6}$-s step.

The integrator drains `saveat` with the pre-modification interpolant before
any `DiscreteCallback.affect!` runs, so a state saved at exactly the park time
carries the pre-park layout: `n_live_cells(cb, t) = n − #{t_park < t}`.
Extraction NaN-pads the ghost slots of `T`, `ξ` and `z`, reports the live
thickness $z_{n_{\rm live}+1}-z_1$, sums mass and gas generation over the live
cells, and exposes `solution.extras.live_cells` and
`solution.extras.park_times`. No dense output is lost and `saveat`
interpolation is unaffected because the state never changes size.

### 14.3.6 The Jacobian through burnout

Because the state keeps its size, the structured Jacobian (§15.6) stays valid
through every park:

- the sparsity pattern is built once for $n$ slots and is a superset of every
  live pattern;
- every operator assembles the live range only, so ghost rows and columns of
  the sparse block are zero;
- the per-cell payload vectors of the structured couplings (the prefix sums of
  the ALE mesh velocity, the Beer–Lambert suffix sums, the $\chi$ rank-1) are
  zero above $n_{\rm live}$, so the sums stop at the live surface and every
  ghost node row reproduces the surface node row, exactly as the residual's
  $\dot z_{\rm ghost} = w_{n_{\rm live}+1}$ requires;
- the geometry AD pass (§15.7) differentiates the masked residual and is
  consistent by construction.

`verify_jacobian` against the dense ForwardDiff reference passes strictly on
parked states for `Structured(strategy = DirectSolve())`, and the parked
`ApproxSolve` assembly equals that of the compacted problem on the live block.

The ALE default Jacobian, `Structured(strategy = ApproxSolve())` (§15.10.2),
stays the default through parking. The `:full` ALE pattern carries the
prefix-sum couplings of every node to every cell below it, which makes the
factorized $W$ effectively dense and the sparse direct solve 60–90× more
expensive per step than the reduced pattern; parking keeps every live cell
above $\tau_{\rm dep}$ of its initial thickness, and the two strategies follow
the same trajectory. The post-park transient is not a
Jacobian effect: `DirectSolve` restarted from $5\times10^{-2}$ s after every
park fails its Newton more than 500 times on the 60-park PMMA burnout, while
the reduced strategy with the reset of §14.3.5 needs 858 steps and 8
failures at a hundredth of the wall time.

The same holds on fast ablation. On a surface regressing at mm/s through a
reaction front one cell wide (the Landau case V-CPL-02, $v = 2$ mm/s,
$\Delta z = 62.5$ µm) the dropped mesh-velocity couplings are larger, but the
reduced modified Newton still contracts by a factor of 0.03 to 0.18 per
iteration at the controller's steps of 2.5 to 5 ms, and the two strategies
follow the same trajectory: without the step cap of §15.5.3, 1088 steps and
894 Newton failures for the reduced strategy against 1064 and 823 for
`DirectSolve` with UMFPACK, which takes several times the wall time. Those
Newton failures are not a Jacobian effect either; a Newton that re-evaluates
the exact Jacobian at every iterate stalls on the same states. They come from
degenerate gas transport (§5.6) and from the step controller (§15.5.3).

---

## 14.4 $h$-refinement and $h$-coarsening (experimental)

The $h$-adaptivity cluster changes mesh *topology* by splitting a cell into two
children ($h$-refinement) or merging two refinement siblings back into their
parent ($h$-coarsening). It lives in `src/Experimental/` and is **not invoked by
`solve`**; it is exercised by tests and the experimental `AMRController`. A mesh
it refines cannot be passed to `solve`: refinement leaves inactive parent cells
among the active ones, while the solver's kernels assume the active cells are the
contiguous prefix `1:n_active_cells`; it also changes the number of degrees of
freedom, which the fixed-length state vector of `solve` cannot follow.

### 14.4.1 Single-cell refinement

`refine_cell!(mesh, cell_id)` (`src/Experimental/h_refinement.jl`)
splits a leaf cell $p$ into a left child $a$ and a right child $b$:

1. **New center node.** Place a node at the geometric midpoint
   $z_{\mathrm{center}} = \tfrac12(z_{\mathrm{left}} + z_{\mathrm{right}})$, with
   its reference coordinate the midpoint of the parent's node references.
2. **Two children.** Child $a$ spans `(left_node, center_node)` and inherits the
   parent's left face; child $b$ spans `(center_node, right_node)` and inherits
   the right face; a new interior face joins them. Both children carry
   `level = parent.level + 1` and `parent = cell_id`.
3. **Repoint the parent's boundary faces** to the children, and record the
   children on the parent (`children = (a,b)`), which marks the parent a non-leaf.
4. **Interpolate state to the children** (see below).
5. **Append geometry caches** (centers, half volumes, the new face position/area)
   and split the dry-solid denominator evenly,
   $m_{\mathrm{dry},a}=m_{\mathrm{dry},b}=\tfrac12 m_{\mathrm{dry},p}$.
6. **Update the active mask**: mark the parent inactive, append two active
   children, increment `n_active_cells` by one (net $-1+2$), and raise
   `max_level`.
7. **Assert AMR invariants.**

**Conservation of mass** is exact. Concentrations are *not* gradient-interpolated;
the children simply inherit the parent value, $\xi_a = \xi_b = \xi_p$, with equal
volumes $V_a = V_b = V_p/2$. Then

$$
V_a\,\xi_a + V_b\,\xi_b = \tfrac12 V_p\,\xi_p + \tfrac12 V_p\,\xi_p = V_p\,\xi_p ,
$$

so each species mass is conserved exactly.

**Temperature interpolation** is where the *code differs from its docstring*, a
discrepancy that must be recorded for "no mistakes." The docstring advertises a
linearly interpolated split,
$\phi_{a,b} = \phi_p \mp (\nabla\phi\cdot\Delta x)/4$. The code computes a gradient
$T_{\mathrm{grad}}$ and forms

$$
T_a = T_p - T_{\mathrm{grad}}\,\frac{|\Delta z|}{2},
\qquad
T_b = T_p + T_{\mathrm{grad}}\,\frac{|\Delta z|}{2},
\qquad
\Delta z = z_{\mathrm{center}} - z^c_{p} .
$$

But $z_{\mathrm{center}} = \tfrac12(z_{\mathrm{left}}+z_{\mathrm{right}})$ is the
*same point* as the parent cell center $z^c_p = \tfrac12(z_{\mathrm{left}} +
z_{\mathrm{right}})$. Hence **$\Delta z = 0$ identically**, the correction
vanishes, and in the current implementation

$$
T_a = T_b = T_p .
$$

The split is thus a piecewise-constant (zeroth-order) injection of the parent
temperature into both children. This is exactly energy-conserving for the
$T\!\cdot\!V$ proxy ($V_a T_a + V_b T_b = V_p T_p$), so `verify_refinement_conservation`
reports zero error, but it does not realize the advertised second-order
reconstruction. With variable heat capacity the constant-$T$ split is also
energy-conserving in the true $\rho c_p T V$ sense because $\rho c_p$ and $T$ are
identical in both children. The gradient routine `estimate_gradient` is correct
and ready; only the geometric $\Delta z$ argument nullifies it.

> **Implementation note.** `refine_cell!` and `estimate_gradient` in
> `src/Experimental/h_refinement.jl`. The conservation check
> `verify_refinement_conservation` compares parent vs. summed
> children for mass (per species, against `volume_ref`) and a $T\!\cdot\!V$ energy
> proxy.

### 14.4.2 Marking and 2:1 balance

`mark_cells_for_refinement(mesh, η, θ_refine, h_min, max_level)`
flags a cell when **all three** criteria hold:

$$
\eta_i > \theta_{\mathrm{refine}}\,\max_j\eta_j,
\qquad
h_i > 2\,h_{\min},
\qquad
\mathrm{level}_i < \mathrm{max\_level}.
$$

The first is the relative-error threshold (default $\theta_{\mathrm{refine}}=0.7$);
the second ensures the children will still exceed $h_{\min}$; the third caps the
refinement depth. The maximum $\eta$ is taken only over active cells, and the
routine returns empty if $\max\eta\le 0$.

`enforce_2to1_balance!(mesh, marked, max_level)` augments the
marked set so that, after refinement, adjacent cells differ by at most one
refinement level. It is a **fixed-point iteration**: for each currently-marked
cell it inspects both neighbors; if a neighbor is coarser
($\mathrm{level}_{\mathrm{nb}} < \mathrm{level}_i$) and below `max_level`, it marks
the neighbor; the loop repeats until no new cell is marked. The 2:1 balance keeps
level transitions gradual, which bounds the interpolation error and the
conditioning of the discretization across refinement jumps. Two known
limitations: the balance is enforced only on *newly* marked cells (an
already-coarse cell that becomes unbalanced after a complex refinement pattern is
not retroactively un-marked), and the fixed-point approach is conservative and can
over-refine.

`perform_h_refinement!(mesh, marked)` refines the marked cells in
decreasing index order (for determinism) and calls `update_geometry!` once at the
end.

### 14.4.3 Sibling coarsening

`coarsen_cells!(mesh, a, b; material)` (`src/Experimental/h_coarsening.jl`) merges two *siblings* (same parent, both leaves) back into their parent.
Its preconditions distinguish it sharply from the depletion park's adjacent-cell merge:
it requires a parent–child relationship (`can_coarsen` checks parent $\ne 0$, both
leaves), whereas the park (`park_cell!`, §14.3.4) works on any
adjacent live pair. The merge
algebra is identical to §14.3.3 — energy-conserving $T$ (or volume-weighted
without material), volume-weighted $\xi$, summed $m_{\mathrm{dry}}$ — but the
*topology* operation restores the parent: the parent's faces are repointed back to
the parent, the parent is marked active with `children = (0,0)`, the two children
are marked inactive, and `n_active_cells` decrements by one.

Crucially, `coarsen_cells!` **does not remove the orphaned intermediate node and
interior face from the arrays** — it only flips the mask. This is the deliberate
"mark-inactive" strategy (§14.6): removing elements would invalidate indices
throughout the code, so the experimental $h$-AMR path leaves orphans in place,
trading memory for index stability. (Contrast the depletion park, which keeps the
arrays' length and shifts the live values down so the live cells stay a contiguous
prefix, §14.3.4.)

`mark_cells_for_coarsening(mesh, η, θ_coarsen, h_max)` marks a
sibling pair when **both** siblings have low error and the merged cell is not too
large:

$$
\eta_i < \theta_{\mathrm{coarsen}}\,\max_j\eta_j,
\qquad
\eta_{\mathrm{sib}} < \theta_{\mathrm{coarsen}}\,\max_j\eta_j,
\qquad
h_i + h_{\mathrm{sib}} \le h_{\max},
$$

with default $\theta_{\mathrm{coarsen}}=0.1$. `perform_h_coarsening!` iterates the
pairs (re-checking both are still active) and refreshes geometry.

### 14.4.4 The AMR controller and combined $hr$-adaptation

`AMRController{EI,QM}` (`src/Experimental/amr_controller.jl`) bundles an error
indicator, a quality metric, the refine/coarsen thresholds, the size and level
bounds, and a check cadence. Its defaults are $\theta_{\mathrm{refine}}=0.7$,
$\theta_{\mathrm{coarsen}}=0.1$, $h_{\min}=10^{-6}\,$m, $h_{\max}=10^{-2}\,$m,
`max_level=4`, `check_every=10`. `adapt_mesh!(mesh, controller)` runs the full cycle:

1. Compute indicators on the current state.
2. Mark for refinement, enforce 2:1 balance, perform $h$-refinement.
3. *Recompute* indicators on the refined mesh.
4. Mark sibling pairs for coarsening, perform $h$-coarsening.
5. Check the quality metric; if it fails, apply Laplacian smoothing
   (`enforce_quality!`, §14.5) and report whether all violations cleared.

It returns `(n_refined, n_coarsened, quality_ok)`. `should_adapt`/`maybe_adapt!`
gate the cycle to fire every `check_every` steps; `maybe_adapt!` (like
`adapt_mesh!`) accepts a `material` keyword that is forwarded to the coarsening pass for the energy-conserving
temperature merge (without it, sibling merges fall back to volume-weighted-T
averaging). `create_default_controller(; h_target)` builds a controller tuned
for pyrolysis: a `CompositeIndicator` of a temperature gradient ($w=1$) and a
temperature curvature ($w=0.5$) combined by `:max` (a `ReactionRateIndicator`
cannot run inside the cycle — see §14.1.3), a `CompositeQualityMetric` of
aspect-ratio ($\le 2$), smoothness ($\le 0.5$), and minimum-size
($h_{\mathrm{target}}/10$), with $h_{\min}=h_{\mathrm{target}}/10$ and
$h_{\max}=100\,h_{\mathrm{target}}$.

`hr_adapt!(mesh; …)` (`src/Experimental/hr_adapt.jl`) runs $h$-adaptivity *first*
(to capture sharp features by topology change) then $r$-adaptivity *second* (to
smooth and optimize node positions via `relocate_nodes!`). The ordering reflects
the design principle that topology changes should precede node relocation so that
relocation operates on the final cell count.

---

## 14.5 Mesh-quality metrics

Quality metrics guard against discretizations that, while topologically valid, are
numerically dangerous: extreme size jumps degrade interpolation accuracy and
matrix conditioning; non-positive volumes signal tangling. `solve` does not call
them; they serve the experimental `AMRController` (§14.4.4) and diagnostics. All
metrics are subtypes of `AbstractQualityMetric` (`src/Materials/Materials.jl`) and live in
`src/Adaptivity/mesh_quality.jl`. Each metric supports `check_quality` (a boolean)
and `quality_violations` (the offending cell indices).

**Aspect-ratio metric.** Bounds the ratio of adjacent active cell sizes:

$$
\frac{\max(h_i,h_j)}{\min(h_i,h_j)} \le \mathrm{max\_ratio},
\qquad \text{for adjacent } i,j,
$$

default $\mathrm{max\_ratio}=2$. Steep size transitions amplify the leading
truncation-error term across a face and worsen the conditioning of the flux
stencil, so a ratio near unity is desirable. The check sorts active cells by
position and compares each consecutive pair.

**Smoothness metric.** Bounds the *relative* size jump:

$$
\frac{|h_j - h_i|}{\max(h_i,h_j)} \le \mathrm{max\_size\_jump},
$$

default $0.5$. This is closely related to the aspect ratio but measures the jump
as a fraction of the larger cell; it is the more natural bound when cell sizes vary
smoothly. (Note the violation test inside `quality_violations` for smoothness uses
$|h_j-h_i|/h_i$ and $|h_j-h_i|/h_j$ per cell — a per-cell normalization — while the
summary statistic `compute_smoothness` reports $|h_j-h_i|/h$ per neighbor.)

**Jacobian metric.** Ensures positive cell volumes (no inverted/tangled cells):

$$
V_i > \mathrm{min\_jacobian}, \qquad \text{default } \mathrm{min\_jacobian}=0 .
$$

The name follows finite-element usage where a non-positive element Jacobian marks
an inverted element; in 1D it is simply the cell volume.

**Minimum-cell-size metric.** Prevents pathologically small cells:

$$
h_i \ge h_{\min}, \qquad \text{default } h_{\min}=10^{-8}\ \mathrm{m} .
$$

**Composite quality metric.** `CompositeQualityMetric` passes only if *all*
sub-metrics pass; its violation set is the union.

**Enforcement by Laplacian smoothing.** `enforce_quality!(mesh, metric;
max_iter=10)` attempts to repair violations by relaxed Laplacian
smoothing of the interior nodes of violating cells:

$$
z_{\mathrm{new}} = \tfrac12\bigl(z_{\mathrm{old}} + \bar z_{\mathrm{neigh}}\bigr),
$$

where $\bar z_{\mathrm{neigh}}$ averages the positions of the nodes of all active
cells sharing the node. Boundary nodes are never moved. The routine iterates until
no violations remain or `max_iter` is reached, refreshing geometry each pass, and
returns the residual violation count. This is a *local, incremental* fix: it can
fail to converge for clustered or severe violations (e.g. after several large
refinements), in which case the AMR controller reports `quality_ok = false`.

**Summary statistics.** `mesh_quality_summary(mesh)` returns the active-cell count,
min/max cell sizes, max aspect ratio, max smoothness deviation, and whether all
volumes are positive — a one-call diagnostic of mesh health.

> **Implementation note.** All metric types, `check_quality`/`quality_violations`,
> `compute_aspect_ratios`/`compute_smoothness`, `enforce_quality!`, and
> `mesh_quality_summary` in `src/Adaptivity/mesh_quality.jl`; all of them mask out
> inactive cells.

---

## 14.6 AMR bookkeeping and invariants

The AMR data model is split for performance (see §13):
**immutable connectivity** (`Cell1D`, `Face1D` are immutable structs; the topology
graph does not change during a residual evaluation) and **mutable geometry**
(positions, volumes, areas, the active mask). Topology changes — the park inside
`solve`, refinement and coarsening outside it — happen only between residual
evaluations; the residual never mutates the mesh.

The central bookkeeping object is the **active-cell mask** `active_cell_mask::
BitVector`, with `n_active_cells` its population count. The hard rule — call it the
*active-cell contract* — is that **every AMR consumer must iterate via the mask (or
the masked accessors), never over a raw `eachindex(mesh.cells)` range**. The reason
is that the experimental $h$-refinement *appends* children to `mesh.cells` and
*marks-inactive* on coarsening without compacting, so `length(mesh.cells)` can
exceed `n_active_cells`. Iterating the raw range would touch stale parent/child
cells — the failure mode the contract (and the masked accessors) exists to
prevent.

`assert_amr_invariants(mesh)` (`src/Adaptivity/Adaptivity.jl`) is the
cheap consistency check called after every refine/coarsen/merge step. It asserts
that all five parallel per-cell arrays — `active_cell_mask`, `cell_states`,
`cell_volumes`, `cell_centers`, and the $\chi$ denominator `m_dry_solid_initial`
(empty allowed only for hand-built meshes that never sized it) — have length
equal to `length(mesh.cells)`, and
that `n_active_cells == count(active_cell_mask)`. It does not check connectivity
correctness; it catches the most common drift, a counter that disagrees with the
mask or a parallel array of the wrong length.

Two distinct slot strategies coexist, by design:

| Operation | Strategy | Arrays | Index validity |
|---|---|---|---|
| $h$-refinement / $h$-coarsening (experimental) | **mark-inactive** | grow / left in place | preserved (no shift) |
| depletion (production) | **park**: shift live values down, ghost at the end | fixed length | live range renumbered (shift), connectivity regular |

The production depletion path shifts because the live cells must stay a
contiguous prefix `1:n_active_cells` in spatial order (§14.3.1); the arrays keep
their length and the 1-D connectivity never changes. The experimental $h$-AMR
path leaves orphans (so its arrays grow monotonically across refinement cycles)
precisely to avoid invalidating indices held elsewhere. A periodic
full-compaction pass for the
$h$-AMR path is an open design question (§14.10).

No workspace buffer changes size on a park: the caches are allocated once for
$n$ slots, and the Jacobian's sparsity pattern is a superset of every live
pattern (§14.3.6).

---

## 14.7 Interaction with the solver

The relationship between depletion parking and the time integrator is mediated
entirely by callbacks (§15) and the read-only-mesh discipline of the residual (§3,
§13); the $h$- and $r$-adaptivity of §14.2 and §14.4 has no such link. Three
facts complete the picture:

1. **Mesh sync is callback-driven.** The residual and Jacobian assembly treat the
   mutable mesh as read-only; the `StepUpdateCallback` mirrors the integrator's
   accepted state (including ALE node positions and the
   $\chi$-driven cross-section area) back into the mesh at the end of each step,
   so that the depletion check, the thickness-termination check, and
   post-processing all see a consistent mesh. The depletion callback runs *after*
   the step-update callback in the `CallbackSet` (see §15).

2. **Termination guards the singular limit.** Alongside depletion, under ALE, the
   thickness-termination callback (`create_thickness_termination_callback`, a
   `ContinuousCallback` on the live thickness
   $z_{n_{\mathrm{live}}+1}-z_1 - \mathrm{min\_thickness}$, default
   $\mathrm{min\_thickness}=10^{-6}\,$m) stops the integration when the whole
   sample has receded to near-zero thickness, before $\rho c_p^{\mathrm{eff}}V$
   becomes singular globally. Depletion handles *per-cell* thinning; termination
   handles *whole-domain* exhaustion. Depletion also self-terminates when only
   `min_cells` remain.

3. **Depletion keeps the Jacobian.** A park rewrites the state in place and
   leaves its size unchanged (§14.3.6), so the structured Jacobian, its sparse
   linear solver and `verify_jacobian` all run through burnout. A park event is
   nonetheless a discrete rewrite of the state: parameter sensitivities (§17)
   are smooth between events and only piecewise smooth across them.

---

## 14.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

The reference codes for condensed-phase pyrolysis take markedly different stances
on mesh adaptivity, because each makes a different choice about how to represent a
receding surface. Mapping their notation to ours on first use:

- **ThermaKin2Ds** (Stoliarov–Lyon) uses a *single moving reference point* and a
  dynamic-mesh adjustment to represent volume change: the deformation source
  $\int (1/\rho)(\partial\rho/\partial t)\,dx$ drives a redistribution of the mesh
  rather than a topology change. There is no error-indicator-driven AMR. The
  Pyrolysis.jl equivalent of ThermaKin's volume-change-driven motion is the ALE
  mesh velocity of §11; what Pyrolysis.jl adds on top is the *surface-depletion
  merge* (§14.3), which removes a cell once it thins past
  $\tau_{\mathrm{dep}}$ — an operation ThermaKin handles by its reference-point
  rescaling rather than by deleting cells. Both share the goal of keeping the
  discretization well-conditioned as the sample recedes; Pyrolysis.jl's merge is
  the discrete, mass/energy-conserving realization of that goal.

- **Gpyro** (Lautenberger–Fernandez-Pello) uses a variable-$\Delta z$
  control-volume mesh: the cell thicknesses themselves are state-like and shrink as
  the solid is consumed, with the front-most control volume allowed to vanish.
  Gpyro's grid stretching toward the heated face is the static analogue of
  Pyrolysis.jl's geometrically stretched mesh (§13.9) plus the interface-proximity
  indicator (§14.1.4) and the interface monitor (§14.2.2), which dynamically pull
  resolution toward the front. Gpyro's "cell disappears" behavior corresponds to
  Pyrolysis.jl's depleted-cell park; the energy-conserving temperature combine
  (§14.3.3) is the careful version of folding a vanishing control volume's energy
  into its neighbor.

- **FDS** (McGrattan et al.) solves the 1D solid-phase heat equation on a *fixed*
  through-thickness grid with a node-collapsing rule as the surface regresses (the
  cell nearest the surface is removed once its mass drops below a threshold), and
  swelling/shrinking is applied through a density-ratio rule on the cell sizes. The
  FDS surface-node removal is conceptually identical to Pyrolysis.jl's
  surface-depletion merge, with the difference that Pyrolysis.jl folds the removed
  cell's sensible enthalpy into its interior neighbor by the enthalpy-conserving
  merge of §14.3.3, whereas FDS's regression bookkeeping is tied to its fixed-grid
  layout.
  Neither FDS nor Gpyro carries an indicator-driven $h$-refinement of the kind in
  Pyrolysis.jl's experimental cluster (§14.4).

The common thread: all four codes must reconcile a fixed (or slowly varying)
numerical grid with a physically receding boundary, and all four ultimately
*remove* resolution at the surface as material is consumed. Pyrolysis.jl's
contribution is to make that removal a clean, conservation-checked,
energy-conserving cell merge with a documented thinness criterion, decoupled from
the rest of the physics by the callback architecture.

> Note: the reference comparison here is drawn from the general published
> descriptions of these codes, so the mapping above is at the level of modeling
> strategy rather than equation-by-equation correspondence. The Pyrolysis.jl side
> of every statement is grounded in the source as cited throughout this chapter.

---

## 14.9 Worked summary of the key relations

For reference, the binding discrete relations of this chapter, all verified against
the source:

**Error indicators** (`src/Adaptivity/error_indicators.jl`):

$$
\eta_i^{\mathrm{grad}} = w\,h_i|\nabla\phi|_i,
\quad
\eta_i^{\mathrm{curv}} = w\,h_i^2|\nabla^2\phi|_i,
\quad
\eta_i^{\mathrm{react}} = w\,h_i\max_j\frac{|r_j|_i}{\max(\xi_{j,i},10^{-10})},
\quad
\eta_i^{\mathrm{iface}} = \frac{w}{1+d_i/h_{\mathrm{target}}} .
$$

**Composite combination**:
$\eta^{\max}=\max_k\hat\eta_k,\ \eta^{\ell_2}=\sqrt{\textstyle\sum_k\hat\eta_k^2},\ \eta^{\Sigma}=\sum_k\hat\eta_k$, with $\hat\eta_k=\eta_k/\max_{\text{active}}\eta_k$.

**Equidistribution** (`src/Adaptivity/r_refinement.jl`):

$$
\int_0^{z_n^\ast}\!\omega\,dz = \frac{n}{N+1}\int_0^L\!\omega\,dz,
\qquad
\omega^{\mathrm{grad}}=\sqrt{1+\alpha|\nabla\phi|^2},
\qquad
\omega^{\mathrm{iface}}=\omega_{\mathrm{base}}+\omega_{\mathrm{peak}}e^{-(z-z_{\mathrm{int}})^2/2\sigma^2}.
$$

**Relaxed node step**:
$z_n^{\mathrm{new}} = z_n^{\mathrm{old}} + \beta(z_n^\ast - z_n^{\mathrm{old}})$,
clamped to $[z_{n-1}^{\mathrm{new}}+h_{\min},\,z_{n+1}^{\mathrm{old}}-h_{\min}]$,
$\beta=0.5$, $\mathrm{tol}=10^{-3}$.

**Thinness criterion** (`src/Solver/callbacks.jl`):
$\Delta z_i < \tau_{\mathrm{dep}}\,\Delta z_{i,\mathrm{init}}$, $\tau_{\mathrm{dep}}=0.05$.

**Conservative merge** (`src/Adaptivity/cell_merge.jl`):

$$
\xi_{\mathrm{merged},j} = \frac{V_a\xi_{a,j}+V_b\xi_{b,j}}{V_{\mathrm{keep}}},
\qquad
\sum_{j} m_j\!\int_{T_{\mathrm{ref}}}^{T_{\mathrm{merged}}}\! c_{p,j}\,dT'
= H_a(T_a) + H_b(T_b),
\qquad
m_{\mathrm{dry},\mathrm{keep}}\mathrel{+}=m_{\mathrm{dry},\mathrm{remove}},
$$

with $V_{\mathrm{keep}} = V_a+V_b$ under ALE and $V_b$ on a fixed mesh, $m_j = V\,\xi_j$ the mass of component $j$, the sum running over every component ($T_{\mathrm{merged}}$ by
bracketed Newton; reduces to the $\rho c_p$-weighted mean for constant $c_p$;
§14.3.3), and $\chi_{\mathrm{merged}} = (m_a\chi_a+m_b\chi_b)/(m_a+m_b)$ with $m$
the initial dry-solid mass (§14.3.4).

**Refinement marking** (`src/Experimental/h_refinement.jl`):
$\eta_i>\theta_{\mathrm{refine}}\max\eta \ \wedge\ h_i>2h_{\min}\ \wedge\ \mathrm{level}_i<\mathrm{max\_level}$,
$\theta_{\mathrm{refine}}=0.7$.

**Coarsening marking** (`src/Experimental/h_coarsening.jl`):
both siblings $\eta<\theta_{\mathrm{coarsen}}\max\eta \ \wedge\ h_i+h_{\mathrm{sib}}\le h_{\max}$,
$\theta_{\mathrm{coarsen}}=0.1$.

**Quality bounds** (`src/Adaptivity/mesh_quality.jl`):
$\max(h_i,h_j)/\min(h_i,h_j)\le 2$ (aspect),
$|h_j-h_i|/\max(h_i,h_j)\le 0.5$ (smoothness),
$V_i>0$ (Jacobian),
$h_i\ge 10^{-8}$ (min size).

**AMR invariant** (`src/Adaptivity/Adaptivity.jl`):
$\mathrm{len}(\text{mask})=\mathrm{len}(\text{states})=\mathrm{len}(V)=\mathrm{len}(\text{centers})=\mathrm{len}(\text{cells})$
and $\mathrm{n\_active\_cells}=\mathrm{count}(\text{mask})$.

---

## 14.10 Limitations and open issues

The adaptivity subsystem is functional in the depletion-parking form `solve` uses
and experimental in its full $h$-/$r$-adaptive form, which `solve` does not use. The following are the principal
open issues, several of which are discrepancies or design questions surfaced by
the code itself.

1. **$h$-refinement temperature split is zeroth-order, not the advertised
   first-order.** Because the split point coincides with the parent cell center,
   $\Delta z = 0$ and the gradient correction in `refine_cell!` vanishes
   identically: both children inherit $T_a=T_b=T_p$ (§14.4.1). This is
   energy-conserving but loses the second-order reconstruction the docstring
   describes. A correct fix would interpolate to the *child* centers
   ($z^c_p \pm \Delta z_p/4$) rather than to the split point.

2. **Mesh-compaction strategy for $h$-AMR.** The experimental refine/coarsen path
   marks inactive without compacting, so `mesh.cells` grows monotonically across
   adaptation cycles. Whether to add a periodic full-compaction pass (reclaiming
   memory at the cost of a global reindex) or to accept the growth for
   long-running simulations is unresolved (§14.6).

3. **Reaction-rate indicator coupling.** The reaction-rate indicator requires
   the per-cell rate matrix (`rates = r`) and therefore cannot run inside the
   `adapt_mesh!` cycle, where indicators are recomputed after topology changes
   and any precomputed rates are stale (it throws rather than running silently).
   Integrating it with the implicit time-stepping — sampling the rates at the
   adaptation point — is not yet done (§14.1.3).

4. **Sensitivities across park events.** A park is a discrete rewrite of the
   state, so trajectories are only piecewise smooth in the parameters across an
   event (§17.7.3); adjoint sensitivity through a park event is unsupported.

5. **$r$-refinement defaults are untuned, and relocation does not remap.**
   The convergence tolerance ($10^{-3}$) and relaxation ($0.5$) defaults have not
   been validated on realistic steep-front pyrolysis problems, and
   `relocate_nodes!` moves nodes without conservatively remapping $T,\xi$ onto the
   new layout — that remap is the caller's responsibility (§14.2.4). A non-smooth
   relocation or a missing remap would also break the AD-cleanliness on which
   forward sensitivity relies (§17).

6. **2:1 balance is one-sided and Laplacian smoothing is fragile.** Balance is
   enforced only on newly marked cells and may leave residual imbalances for
   complex patterns (§14.4.2); the Laplacian quality enforcer is a local
   incremental fix that can fail to converge for clustered or severe violations
   (§14.5).

---

# 15. Temporal Integration, DAE/Residual Formulation, and the Structured Jacobian

This chapter describes how the spatially-discretized governing equations of §13
(Finite-Volume Spatial Discretization) are assembled into a single flat
ODE/DAE right-hand side, integrated in time by stiff implicit Runge–Kutta
methods, and linearized through a structured-operator Jacobian. We derive the
block-major state-vector layout and the two dispatched residual paths; the
unified `residual!` entry point and the immutable-definition / mutable-workspace
split that makes it pure and thread-safe; the full ALE residual assembly
sequence; the link to OrdinaryDiffEq's implicit-RK $W$-matrix; the
three-layer structured-Jacobian architecture (operators → couplings →
linear-solve strategies); the exact per-operator Jacobian formulas with their
sign conventions; the O(N) prefix/suffix/rank-1 couplings that represent the
non-local terms; the linear-solve strategies and iterative solvers; the
two-pattern sparsity construction; and the verification
machinery that anchors every backend to a dense AD reference. Every discrete
rule is grounded in the implementing function so the mathematics is traceable to
the code; where the code makes an approximation, a guard, or a sign choice, we
state it and trust the code over any idealized form.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$ (axial, through-thickness). **$z=0$ is the
  bottom/substrate** (boundary id 2, tag `:bottom`); **$z=L$ is the top/exposed
  surface** (boundary id 1, tag `:top`). Cell indices increase from bottom to
  top.
- **Flux sign convention:** positive flux is transport in the $+z$ direction.
- **Divergence convention:** the discrete divergence $(F_R-F_L)A/V$ is positive
  when the quantity flows *out* of the cell; the energy RHS carries
  $-(\nabla\cdot\mathbf{q})$.
- **Heat-of-reaction sign convention:** $h>0$ endothermic (cools), $h<0$
  exothermic (heats); the volumetric source is $Q_{\rm rxn}=-\sum_r h_r r_r$
  (overload resolution H1 of §2).
- In this chapter $i$ is the **cell index** ($1,\dots,n$), $f$ a face index, $j$
  a component index, $r$ a reaction index, $k$ either a node index or a generic
  loop index where stated. In the linear-algebra sections, $\mathsf{A}$ denotes
  the sparse Jacobian block (not the cross-sectional area $A$); the disambiguation
  follows overload resolution A3 of §2 and is stated where both appear.

The state vector is denoted $\boldsymbol{u}$; its time derivative is
$\boldsymbol{du}=d\boldsymbol{u}/dt$. The residual function (in the
method-of-lines sense, an *explicit* RHS, not an implicit DAE residual
$F(\boldsymbol{u},\boldsymbol{du},t)=0$) is `residual!`, which writes
$\boldsymbol{du}=\mathrm{RHS}(\boldsymbol{u},t)$.

---

## 15.1 State layout and the simulation-mode trait

### 15.1.1 Block-major packing

All time-dependent unknowns are packed into a single flat vector
$\boldsymbol{u}\in\mathbb{R}^{N_{\rm state}}$ in **block-major** order, with the
block set determined by the active `SimulationMode`. The descriptor is
`StateLayout`:

$$
\boldsymbol{u} =
\big[\,\underbrace{T_1,\dots,T_n}_{T\ \text{block}}\ \big|\
\underbrace{\xi_{1,1},\dots,\xi_{1,n}}_{\xi_1\ \text{block}},\dots,
\underbrace{\xi_{N_C,1},\dots,\xi_{N_C,n}}_{\xi_{N_C}\ \text{block}}\ \big|\
\underbrace{z_1,\dots,z_{n_{\rm nodes}}}_{z\ \text{block (ALE)}}\ \big|\
\underbrace{\chi_1,\dots,\chi_n}_{\chi\ \text{block (WithChi)}}\,\big].
$$

The explicit index map (used verbatim by every pack/unpack routine and by the
Jacobian's `FastSparsityMap`) is

$$
\begin{aligned}
\text{T index}        &: \quad i, &&i=1,\dots,n,\\[2pt]
\text{$\xi_{j,i}$ index} &: \quad n + (j-1)\,n + i, &&j=1,\dots,N_C,\ i=1,\dots,n,\\[2pt]
\text{$z_k$ index}     &: \quad n(1+N_C) + k, &&k=1,\dots,n_{\rm nodes}\ \text{(ALE only)},\\[2pt]
\text{$\chi_i$ index}  &: \quad n(1+N_C) + \delta_{\rm ALE}\,n_{\rm nodes} + i, &&i=1,\dots,n\ \text{(WithChi only)},
\end{aligned}
$$

and the total length is

$$
N_{\rm state} = n\,(1+N_C) + \delta_{\rm ALE}\,n_{\rm nodes} + \delta_{\rm WithChi}\,n,
$$

with $\delta_X=1$ if the mode includes block $X$ and $0$ otherwise. In 1D
$n_{\rm nodes}=n+1$. The $\chi$ block is **appended after** the $z$ block, not
interleaved; its offset arithmetic is therefore dense and is the single source
of subtle indexing bugs the code guards against by always going through the
range accessors `temperature_range`, `concentration_range(l,j)`, `z_range(l)`,
and `chi_range(l)`.

> **Implementation note.** The layout, accessors, and the
> `pack_state!`/`update_mesh_from_state!`/`create_state_vector` operations live
> in `StateLayout`, src/Problem/state_layout.jl. The same packing is mirrored in
> the Jacobian sparsity module's index helpers `_xi_idx` and `_z_idx`,
> src/Jacobian/sparsity.jl, and in solution extraction,
> src/Solver/solve.jl. The block order $[T\mid\xi\mid z\mid\chi]$ is binding: any
> code that reads or writes raw state must respect it.
>
> **IC contract.** `create_state_vector` packs the mesh **as-is** — it does not
> consult the problem's initial-condition closures (which live on the
> `PyrolysisProblem`, not the mesh) and does not populate
> `mesh.m_dry_solid_initial`. The canonical setup path is the exported
> `initialize_state(problem; …) -> (def, ws, u0)`
> (src/Problem/adapter.jl), which mirrors `solve()` exactly:
> `apply_initial_conditions!(problem)` → `to_problem_def` → `build_workspace` →
> `create_state_vector`. Hand-rolled compositions that skip the IC step
> integrate from the mesh constructor's placeholder state — a
> resolution-independent error that still returns `retcode == :Success` —
> and, under `WithChi`, leave the χ denominator at zero (§12.8.3).

The block-major (rather than cell-major / interleaved) choice is deliberate. It
makes each physical field a contiguous slice, which (i) lets `block_view` expose
a strided sub-vector for diagnostics and callbacks without
copying, (ii) keeps the conduction/transport stencils as small offsets within a
block, and (iii) makes the ALE geometry coupling a clean rectangular sub-block
(the $z$-columns), which the two-pattern sparsity construction of §15.11
exploits.

The experimental P1 radiation path (quarantined to
`Pyrolysis.Experimental.RadiationP1`) uses this *same* layout: the
incident-radiation field $G$ is an algebraic (quasi-steady) field stored in
the `RadiationCache` and re-solved at every residual evaluation, never a
state block (see §8.12). The layout above is therefore the canonical one for
every mode, production and experimental alike.

### 15.1.2 The bundled `SimulationMode` trait

Rather than carry four orthogonal trait axes and suffer combinatorial method
explosion, Pyrolysis.jl bundles all four into a single type parameter:

$$
\texttt{SimulationMode}\{\,\underbrace{\text{Mesh}}_{\rm Eulerian/ALE},\ \underbrace{\text{Rad}}_{\rm NoRadiation/SurfaceAbsorption/BeerLambert},\ \underbrace{\text{Trans}}_{\rm Fickian/DarcyFick},\ \underbrace{\text{Chi}}_{\rm NoChi/WithChi}\,\}.
$$

The four axes are:

| Axis | Type | Values | Effect |
|---|---|---|---|
| Mesh motion | `MeshMode` | `Eulerian`, `ALE` | $z$ block present iff ALE; geometry constant iff Eulerian |
| Radiation transport | `RadiationMode` | `NoRadiation`, `SurfaceAbsorption`, `BeerLambert` | volumetric source $Q_{\rm rad}$ and Jacobian non-locality |
| Gas transport | `TransportMode` | `Fickian`, `DarcyFick` | diffusion-only vs. pressure-driven Darcy + diffusion |
| Cross-section | `ChiMode` | `NoChi`, `WithChi` | $\chi$ block present iff `WithChi`; $A(t)=A_0\,\mathrm{law}(\bar\chi)$ |

Accessors extract the type parameters at compile time:

```julia
@inline mesh_mode(::SimulationMode{M,R,T,C})      where {M,R,T,C} = M()
@inline radiation_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = R()
@inline transport_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = T()
@inline chi_mode(::SimulationMode{M,R,T,C})       where {M,R,T,C} = C()
```

Because the mode is a type parameter, `mesh_mode(def.mode) isa ALE` and friends
resolve to compile-time constants, and the residual's dispatch chain
specializes fully — no runtime branching on the mode. The residual entry
dispatches on the mesh mode only, giving **two residual paths**, Eulerian and
ALE. Inside each, the transport law follows the property cache (a Darcy sub-cache
exists only in `DarcyFick` mode, and `compute_rhs!` selects the Darcy–Fick flux
from it), the $\chi$ block is included iff `WithChi`, and the radiation model
specializes the source.

> **Implementation note.** `SimulationMode`, the `MeshMode`/`RadiationMode`/
> `TransportMode`/`ChiMode` abstract hierarchies, and the accessors are in
> src/Problem/traits.jl. The mode is *derived*, not user-set directly:
> `to_problem_def(problem; use_ale, radiation_model, radiation_intensity,
> radiation_kirchhoff)` (src/Problem/adapter.jl) assembles it from the `use_ale`
> flag, the radiation
> model, and the material's Darcy/$\chi$ traits (`has_darcy_flow`, the presence
> of a `lateral_shrinkage_law`). The `RadiationModel` enum value
> `P1_QUASI_STEADY` is rejected at solve setup (it lives in `Experimental`), so
> the production `RadiationMode` axis spans only the first three values.

---

## 15.2 The problem definition / workspace separation

The single most important structural decision in the temporal layer is that the
residual is a **pure function** of state and time, parameterized by an
*immutable* problem definition and a *mutable* scratch workspace:

```julia
residual!(du, u, def::ProblemDef, t, ws::Workspace)  ->  nothing
```

- `def::ProblemDef` is **immutable** and shareable across threads / sweep
  workers. It carries the mesh connectivity (`MeshTopology`, a placeholder the
  1D residual does not read) and the *unchanging* specification: the material,
  the BC set, the `SimulationMode`, and a NamedTuple of runtime options
  (`radiation_model`, `radiation_intensity`, `radiation_kirchhoff`).
- `ws::Workspace` is **mutable** and owned by exactly one worker. It holds the
  typed state mirror (`StateCache`), property cache (`PropertyCache`), RHS work
  arrays (`RHSCache`), and (in ALE mode) the ALE work arrays. It also caches the
  dry-gas component indices and the boundary face IDs.
- `t` is a plain `Real`; the kernel accepts whatever scalar type $\boldsymbol{u}$
  and $\boldsymbol{du}$ carry — `Float64` for ordinary evaluation, or
  `ForwardDiff.Dual` for AD Jacobian assembly and forward sensitivity (§17).

This split is what makes parameter sweeps (§ User Guide UG-12) and the
AD Jacobians possible without re-allocating problem metadata: every worker shares
one `def` and owns one `ws`; the AD Jacobian builds a *parallel* `Dual`-typed
workspace (`build_ad_workspace`) that shares the mesh, layout, and boundary
indices but carries independent `Dual`-typed caches.

A defining invariant is that **the mutable `Mesh1D` is read-only during residual
and Jacobian assembly** — connectivity only. All state and all state-derived
geometry (cell centers, volumes, face areas, cross-section) live in the typed
`StateCache`, which is populated once at residual entry by `populate!` and
carries `eltype(u)` end-to-end. The mutable mesh is re-synchronized from the
accepted integrator state *only between steps*, by the `StepUpdateCallback`
(§15.5). This guarantees that the many residual evaluations inside one implicit
step (Newton iterations, Jacobian evaluations) all see a
consistent geometry, and that AD sensitivities propagate cleanly through
geometry because the geometry is a function of `cache.z` and `cache.A_section`,
both `Dual`-typed in the AD pass.

> **Implementation note.** `ProblemDef` and the runtime-option NamedTuple:
> src/Problem/problem_def.jl. `Workspace`: src/Problem/workspace.jl;
> `build_workspace` and `build_ad_workspace`: src/Problem/adapter.jl. `StateCache`, `populate!`,
> `update_state_geometry!`, `update_cross_section_area!`:
> src/Problem/state_cache.jl. The user-facing `PyrolysisProblem` and the
> `to_problem_def` adapter: src/Problem/pyrolysis_problem.jl and
> src/Problem/adapter.jl. The canonical `(def, ws, u0)` front door is
> `initialize_state(problem; …)` (same file) — it applies the problem's
> initial conditions before packing and is the path `solve()` itself takes;
> `to_problem_def`/`build_workspace`/`create_state_vector` are public as
> the advanced layer (AD workspaces, mid-simulation repacking, warm starts)
> with the IC step on the caller (§15.1.1).

---

## 15.3 The unified residual and its dispatch

The top-level entry is a one-line dispatcher on the mesh mode:

```julia
@inline function residual!(du, u, def::ProblemDef, t, ws::Workspace)
    _residual!(du, u, def, t, ws, mesh_mode(def.mode))
    return nothing
end
```

### 15.3.1 Eulerian path

For a fixed mesh the geometry is constant (mirrored at workspace construction),
so the residual is a pure function of $(\boldsymbol{u},t)$ in four steps:

1. `populate!(ws.state_cache, u, layout)` — unpack $T,\xi$ from $\boldsymbol{u}$
   into the typed cache (single source of truth).
2. `Properties.update_properties!` — evaluate all state-dependent properties
   ($k_{\rm eff}$, $\rho c_p^{\rm eff}$, $\lambda$, $\varepsilon$, $\alpha$, and
   face-averaged $k_{\rm face}$, $\lambda_{\rm face}$).
3. `_compute_physics_rhs!` — calls `Discretization.compute_rhs!`, which picks
   the transport law from the property cache's Darcy sub-cache (`nothing` ⇔
   Fickian), refreshes that sub-cache ($P,\varphi,\kappa,\mu$ and the Darcy
   velocity, `update_darcy_properties!`) in its first step under `DarcyFick`, and
   populates `ws.dT_work[i]`, `ws.dξ_work[j,i]`, and
   `ws.rhs_cache.dξ_rxn[j,i]`.
4. Pack $\boldsymbol{du}$: $du[i]=dT_i/dt$, $du[n+(j-1)n+i]=d\xi_{j,i}/dt$.

The discrete updates assembled in step 3 are exactly those derived in §7 and
§13 (repeated here for completeness):

$$
\frac{dT_i}{dt} = \frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i + S_{{\rm conv},i} + Q_{{\rm rad},i} + Q_{{\rm rxn},i}\Big],
$$

$$
\frac{d\xi_{j,i}}{dt} = -(\nabla\cdot\mathbf{J}_j)_i + S_{j,{\rm rxn},i}\quad(\text{gas}),\qquad
\frac{d\xi_{j,i}}{dt} = S_{j,{\rm rxn},i}\quad(\text{solid/liquid}),
$$

with the face fluxes, harmonic-mean conductivity, smooth-gated Arrhenius rates,
and the radiation model all as specified in their respective chapters. The
reaction heat source carries the storage sign convention,
$Q_{{\rm rxn},i}=-\sum_r h_r\,r_{r}(T_i,\xi_i)$ (endothermic $h>0$ gives
$Q_{\rm rxn}<0$); the conduction divergence enters with a leading minus so that
outflow cools the cell.

> **Implementation note.** Eulerian path: `_residual!(…, ::Eulerian)`,
> src/Residual/residual.jl. The 8-step `compute_rhs!` assembly is in
> src/Discretization/finite_volume.jl (see §7.6 for the strict step ordering).

### 15.3.2 Order dependencies inside the physics RHS

The 8-step `compute_rhs!` sequence (transport properties and conductive fluxes →
reactions → gas transport → thermal BCs → radiation → $S_{\rm conv}$ →
divergences → final RHS) is *ordered* for correctness, not convenience:

- **Gas fluxes before thermal BCs.** The mass BCs (step 3) must finish before
  the boundary-temperature Newton solve (step 4), because the blowing model of
  each convective BC reads the outward gas mass flux $\dot{m}_g$ at the face.
- **Thermal BCs before $S_{\rm conv}$.** The boundary cells' temperature
  gradient in $S_{\rm conv}$ (step 6) blends in the boundary-face slope at the
  solved $T_s$, read from the flux that step 4 wrote.
- **Reaction rates before mesh velocity (ALE).** The ALE mesh velocity is built
  from the per-cell strain rate $\theta_i$, which is a function of the species
  reaction sources `dξ_rxn`; hence the physics RHS (which fills `dξ_rxn`) runs
  before the mesh-velocity integration.

---

## 15.4 The full ALE residual assembly

When `mesh_mode(def.mode) isa ALE`, mesh motion is part of the coupled system,
and the residual additionally evolves the node positions $z_k$ and (under
`WithChi`) the pyrolysis-progress $\chi_i$. The assembly is the following
sequence; we present it as the canonical pipeline and give the governing
equation for each step. Cross-references: the ALE kinematics and GCL are in §11;
the volume-change / $\chi$ model is in §10.

1. **Populate state cache.** Unpack $T,\xi,z,\chi$ from $\boldsymbol{u}$ into the
   typed cache.

2. **Update $\chi$-driven cross-section** (`update_cross_section_area!`,
   dispatched on `chi_mode`). Under `WithChi`, compute the mass-weighted
   column-average progress and set the column area:
   $$
   \bar\chi = \frac{\sum_i \chi_i\,m_{{\rm dry},i}}{\sum_i m_{{\rm dry},i}},\qquad
   A = A_0\,\mathrm{law}(\bar\chi),
   $$
   storing $A$ in `cache.A_section[]` carrying `eltype(u)` (AD-safe). The
   numerator accumulates `Dual` types; the denominator $M_{\rm dry}$ is a
   `Float64` constant stored on the mesh. Under `NoChi` this is a no-op
   returning $A_0$.

3. **Update state-derived geometry** (`update_state_geometry!(…, ALE())`).
   Recompute face positions $z_{f}=z_{{\rm node}(f)}$, face areas $A_f=A$, cell
   centers $z_i^c=\tfrac12(z_L+z_R)$, and cell volumes $V_i=|z_R-z_L|\,A$, all
   from `cache.z` and `cache.A_section`, preserving the scalar type.

4. **Update properties.** As in the Eulerian path, but now on the moving
   geometry.

5. **Physics RHS.** Same `compute_rhs!` as Eulerian; importantly it also fills
   `ws.rhs_cache.dξ_rxn`, which the ALE steps reuse.

6. **Column-global radial strain and cumulative mesh velocity.** First compute
   $\theta_A$ (the column-global radial channel; $0$ under `NoChi`), then
   integrate the axial strain bottom-up to obtain the node velocities. The
   per-cell volumetric strain rate is the volume-change law of §10,
   $$
   \theta_i = \frac{1}{V_i}\frac{dV_i}{dt} = \sum_{j}\gamma_j\,\frac{1}{\rho_j}\,\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn},
   $$
   (only $\gamma_j>0$ components contribute; gases with $\gamma=0$ escape without
   changing matrix volume). The axial part removes the radial channel,
   $\theta_{L,i}=\theta_i-\theta_A$, and the cumulative integration with pinned
   bottom node $w_1=0$ gives
   $$
   w_{k+1} = w_k + (\theta_i-\theta_A)\,\Delta z_i,\qquad \Delta z_i = \frac{V_i}{A},
   $$
   i.e. $w_k=\sum_{i<k}\theta_{L,i}\,\Delta z_i$ (an exclusive cumulative sum;
   this is exactly the structure of the prefix-sum couplings in the Jacobian,
   §15.9). The radial channel under `WithChi` is, by the chain rule,
   $$
   \theta_A = \frac{1}{A}\frac{dA}{dt} = \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\frac{d\bar\chi}{dt},\qquad
   \frac{d\bar\chi}{dt} = \frac{1}{M_{\rm dry}}\sum_i S_{{\rm dry},i}\,V_i,
   $$
   with $\mathrm{law}'(\bar\chi)$ obtained by `ForwardDiff.derivative(law, χ̄)`
   and $S_{{\rm dry},i}=\sum_{j\in\text{dry-gas}}(d\xi_{j,i}/dt)|_{\rm rxn}$. The
   guard `law(χ̄) > ε` prevents the radial channel from diverging as the
   cross-section collapses.

7. **Mesh-relative temperature advection** (`compute_ale_advection_term!`,
   override-only). Every component rides the mesh, so the thermal mass
   $\rho c_p^{\rm eff}$ moves with it and $T$ is advected only by the mesh
   velocity relative to the material, $w_{\rm rel} = w_{\rm applied} −
   w_{\rm physical}$ (§11.2.2): in production $w_{\rm rel} \equiv 0$ and nothing
   fires; under `Workspace.w_override` the TVD-limited kernel adds
   $w_{\rm rel}\,\partial T/\partial z$ to $dT$.

8. **Dilation correction** for every component
   (`_apply_ale_dilation_correction!`):
   $$
   \frac{d\xi_{j,i}}{dt}\bigg|_{\rm dilation} \mathrel{-}= \xi_{j,i}\,\theta_i\quad(\text{every } j),
   $$
   which keeps each component's cell mass on its reaction-and-flux budget as
   the cell volume changes: $d(\xi_jV)/dt = S_jV − [J_jA]_{\rm faces}$. The gas
   is included — it is carried with the pores it occupies, and its flux
   $J_g$ is relative to the skeleton (§11.2.3).

9. **Pack the $T,\xi,z$ blocks.** $du[i]=dT_i/dt$ (with the ALE correction),
    $du[n+(j-1)n+i]=d\xi_{j,i}/dt$ (with dilation), and
    $du[n(1+N_C)+k]=w_k$ — i.e. $dz_k/dt=w_k$ (Euler-forward node motion
    realized by the integrator).

10. **Pack the $\chi$ block** (`_pack_chi_derivative_dispatched!`, `WithChi`
    only). The progress variable obeys
    $$
    \frac{d\chi_i}{dt} = \frac{S_{{\rm dry},i}\,V_i}{m_{{\rm dry},i}},\qquad
    S_{{\rm dry},i}=\sum_{j\in\text{dry-gas}}\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn},
    $$
    with a `max(m_init, eps())` guard against division by zero. Note $\chi$ is
    *reconstructed* from `dξ_rxn` each residual rather than carried as an
    independent integral — it is still a genuine ODE state (the integrator
    advances it), but its derivative is recomputed from first principles every
    evaluation, keeping it exactly in sync with cumulative dry-gas production.
    Under identity-law materials the $\chi$ derivatives are zeroed.

The complete ALE discrete updates therefore read

$$
\frac{dz_k}{dt}=w_k,
$$
$$
\frac{dT_i}{dt}=\frac{-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}}{\rho c_p^{\rm eff}{}_{,i}}
\;\Bigl[\,+ w_{{\rm rel},i}\Big(\frac{\partial T}{\partial z}\Big)_i \text{ under \texttt{w\_override} only; } w_{\rm rel}\equiv 0 \text{ in production}\Bigr],
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i}
- \underbrace{\xi_{j,i}\,\theta_i}_{\text{dilation (every component)}}.
$$

> **Implementation note.** ALE path: `_residual!(…, ::ALE)`,
> `_compute_ale_mesh_velocity!`, `_compute_theta_A` / `_compute_theta_A_dispatched`,
> `_apply_ale_T_advection!`, `_apply_ale_dilation_correction!`, and
> `_pack_chi_derivative!` are all in src/Residual/residual.jl. The strain-rate
> law `volume_change_rate` is in src/Physics/volume_change.jl. The ALE flux
> kernels are in src/Discretization/ale_fluxes.jl.

---

## 15.5 Time integration

### 15.5.1 The method-of-lines ODE and the integrator

The packed residual turns the PDE system into the autonomous-in-form ODE
$\boldsymbol{du}/dt=\mathrm{RHS}(\boldsymbol{u},t)$, which `solve` hands to
OrdinaryDiffEq.jl. Pyrolysis kinetics are *stiff* (Arrhenius source terms span
many decades over a narrow temperature window; conduction adds the usual
parabolic stiffness $\sim\Delta t/\Delta z^2$), so the default integrator is an
L-stable singly-diagonally-implicit Runge–Kutta (SDIRK) method:

$$
\texttt{KenCarp4}\big(\text{autodiff}=\text{AutoFiniteDiff()}\big).
$$

KenCarp4 is a 4th-order ESDIRK pair with an embedded error estimator for
adaptive step-size control; its singly-diagonal structure means every implicit
stage shares the same Newton iteration matrix, which we exploit below.

**Default linear solver.** `_resolve_integrator` picks the sparse direct
factorization from the Jacobian's pattern: `UMFPACKFactorization` when the
Jacobian materializes the ALE geometry block (`DirectSolve` under ALE,
`materializes_geometry_block`, §15.11.2), and `KLUFactorization` for every
other sparse Jacobian. The rule depends on the mode and the strategy, never on
the mesh size. It rests on end-to-end KenCarp4 runs of 100-cell PMMA and
charring-wood decks (Eulerian and ALE, surface and Beer–Lambert absorption,
with and without a lateral shrinkage law) and of V-BO-01:

- On the local patterns (every `ApproxSolve` Jacobian, Eulerian `DirectSolve`
  with surface absorption) KLU is 2.6–4.4× faster than Sparspak and 1.1–1.4×
  faster than UMFPACK, allocates 8–83 MB per solve against 19–534 MB for
  Sparspak, and takes the same steps to within 1%.
- On the exact ALE Jacobian UMFPACK is 1.5–1.9× faster than Sparspak and
  2.6–2.7× faster than KLU, with identical step and Newton counts.
- On the Eulerian Beer–Lambert pattern of `DirectSolve` KLU is 1.7× faster
  than Sparspak. UMFPACK's threshold pivoting (relative pivot tolerance 0.1)
  produced a factor with growth near $10^{50}$ on a $W$ of condition number
  about 30; the resulting Newton steps failed, and the run took 50 times the
  steps. That pattern therefore uses KLU.
- All three factorizations are backward stable on these systems (normwise
  backward error below $3\times10^{-15}$). Sparspak does not pivot: on the
  Darcy decks its componentwise backward error reaches 1 in individual rows,
  where KLU's and UMFPACK's stay below $10^{-13}$. On the exact ALE Jacobian the
  relative residuals $\lVert Wx-b\rVert/\lVert b\rVert$ reach
  $6\times10^{-11}$ for UMFPACK, $2\times10^{-11}$ for KLU and
  $3\times10^{-9}$ for Sparspak, without effect on the Newton iteration.

Sparspak, KLU, UMFPACK and the ILU-Krylov solvers of §15.12 are interchangeable
through `linear_solver`.

> **Implementation note.** Integrator resolution: `_resolve_integrator`,
> src/Solver/solve.jl. A user may pass a constructor (`integrator=Rodas5P`) plus
> a `linear_solver=`, or a fully-built algorithm. The RHS closure re-samples
> the incident radiation at every evaluation (`ws.I_radiation = radiation_flux(t)`)
> whenever the active model consumes it (`SURFACE_ABSORPTION` / `BEER_LAMBERT`)
> — constant and time-varying fluxes alike (§8.3.4). The Jacobian closure
> refreshes `ws.I_radiation` at its own $t$ the same way, so a time-varying flux
> enters the Beer–Lambert Jacobian at the time the Jacobian is evaluated.

**Problem-independent integrator types.** OrdinaryDiffEq compiles its stepper,
nonlinear solver and factorization for the types of the functions it is handed.
The RHS and Jacobian closures capture the `ProblemDef` and workspace, whose
types carry the material, boundary-condition and mode types. `solve` therefore
hands them to the integrator behind fixed-signature `FunctionWrapper`s:
`ErasedRHS` for `(du, u, p, t)` with `Vector{Float64}` state and `Float64`
time, and `ErasedJacobian{M}` for the Jacobian matrix type `M`. The callbacks
do the same (§15.5.3). The integrator's type then depends only on the solver
configuration: the mesh mode, depletion handling, the installed optional
callbacks, the integrator and its factorization. A problem with new material or
boundary-condition types compiles only the Pyrolysis kernels (initialization,
residual, Jacobian, solution extraction). The precompile workloads
(`src/precompile.jl`) cache the integrator for the default configurations. An
algorithm that calls the RHS with other argument types, such as dual numbers
for a ForwardDiff time derivative, reaches the closure by dynamic dispatch. A
call of the Jacobian with a different matrix type is a `MethodError`, because
a converted copy would silently drop the entries.

> **Implementation note.** src/Solver/integrator_functions.jl (`ErasedRHS`,
> `ErasedJacobian`); the `ODEFunction` is built in `solve`, src/Solver/solve.jl.
> `solve` also converts the tolerances, step bounds and `saveat` to `Float64`,
> because the solve keywords enter the integrator's type.
> test/solver/test_integrator_type.jl checks that two decks differing in every
> problem type share one integrator type.

### 15.5.2 The implicit-RK $W$ matrix and warm start

Each implicit SDIRK stage solves a nonlinear algebraic system for the stage
increment by Newton's method. The Newton iteration matrix — OrdinaryDiffEq's
"$W$" matrix — is

$$
\boxed{\,W = \frac{1}{\gamma_{\rm RK}\,\Delta t}\,I - J\,,}\qquad J = \frac{\partial\,\mathrm{RHS}}{\partial\boldsymbol{u}},
$$

where $\gamma_{\rm RK}$ is the (common, by the singly-diagonal property) stage
coefficient and $J$ is the residual Jacobian. (Equivalently $W$ may be scaled by
$\gamma_{\rm RK}\Delta t$ to $W=I-\gamma_{\rm RK}\Delta t\,J$; the scaling is an
internal detail of the integrator. The factor that matters for our purposes is
that $W$ is a *diagonal shift of $-J$*, so the diagonal of $J$ adds to the
$1/(\gamma_{\rm RK}\Delta t)$ scale.) Each Newton step solves
$W\,\delta\boldsymbol{u}=-r$ for the stage residual $r$ with the sparse
factorization of §15.5.1; the strategy (`DirectSolve` / `ApproxSolve`, §15.10)
decides which entries of $J$ are assembled.

The Jacobian $J$ is supplied to the integrator through an `ODEFunction` carrying
both `jac = jac!` and a `jac_prototype` sparsity pattern (the union of the
physics and geometry patterns of §15.11). `solve` builds the spec for its own
`(def, ws)`:

```julia
jac_spec = JacobianSpec(backend, def, ws)
jac! = (J, u, p, t) -> compute_jacobian!(J, jac_spec, u, def, t, ws)
jac_prototype = Float64.(jac_spec.sparsity_phys .| jac_spec.sparsity_geom)
ode_fn = ODEFunction{true,FullSpecialize}(ErasedRHS(rhs!);
    jac = ErasedJacobian{typeof(jac_prototype)}(jac!),
    jac_prototype = jac_prototype)
```

(with the incident flux refreshed at the Jacobian's own $t$ inside `jac!`,
§15.5.1).

**Warm start.** A prior solution's final state can seed a new solve's initial
condition (`warm_start=prev_solution`), provided the component count, mesh size,
ALE mode, and $\chi$ mode all agree (enforced implicitly by a state-vector
length check) and the prior solve parked no cells. The integrator's `tspan[1]` still defines the start time, so
warm-starting treats the prior final state as fresh initial data — useful for
continuation studies and for sensitivity sweeps, with the caveat (§17.7) that
warm-start artifacts can perturb finite-difference sensitivities.

> **Implementation note.** ODE construction, Jacobian wiring, and warm-start
> application: src/Solver/solve.jl (`solve`, `_apply_warm_start!`). The
> per-step linear-solve protocol is OrdinaryDiffEq's; Pyrolysis supplies only
> $J$ and the strategy.

### 15.5.3 Output-time handling and callbacks

Snapshots are saved at `saveat=problem.output_times` (11–101 evenly spaced
samples by default). Callbacks are assembled in a fixed order and combined into
a `CallbackSet`:

| Order | Callback | Type | Role |
|---|---|---|---|
| 1 | `StepUpdateCallback` | discrete (every accepted step) | mirror integrator $\boldsymbol{u}\to$ mutable mesh; sync geometry & cross-section |
| 2 | thickness termination | continuous (ALE) | terminate when the live thickness $z_{n_{\rm live}+1}-z_1<$ `min_thickness` |
| 3 | `DepletionCallback` | discrete (opt-in) | park cells below their depletion threshold; the state keeps its size |
| 4 | Newton-failure cap | discrete (every accepted step; Newton integrators) | after a step whose Newton needed retries, hold the proposed step at or below 1.25× that step for 20 steps |
| 5 | progress bar | discrete (every accepted step) | progress bar |
| 6 | `DiagnosticsCallback` | discrete (every `diagnostics_stride`) | log a sensible-energy snapshot row (ledger fields NaN; §16.8.3) |
| 7 | Jacobian verification | discrete (at sample times) | compare active Jacobian vs. dense AD reference |
| 8 | user callback | — | custom |

Every callback of rows 1–7 keeps its problem-typed work behind a fixed-signature
`FunctionWrapper` on the state vector and the time: the mesh sync, the
live-thickness evaluation, the depletion check and park, the diagnostics row and
the Jacobian comparison. Only the actions on the integrator itself stay in the
callback: `terminate!`, `set_proposed_dt!`, `derivative_discontinuity!` and the
reads of its step and Newton-failure counters. The `CallbackSet`'s type
therefore depends on which callbacks are installed, not on the problem (§15.5.1).

The **`StepUpdateCallback`** is the linchpin of the read-only-mesh design: after
each accepted step it sets `mesh.cross_section_area = cross_section_area_from_state(…)`
and calls `update_mesh_from_state!(mesh, u, layout)`, so all between-step
consumers (depletion detection, diagnostics, the surface BC solve in
post-processing) see a mesh consistent with the integrator state.

The **`DepletionCallback`** (§14.3) detects cells whose current thickness has
fallen below `depletion_threshold × Δz_init,i` (default $5\%$ of each cell's
own initial thickness; on a fixed mesh, the top live cell whose condensed
concentration has fallen below the same fraction of its initial value), merges
it into a neighbour with the sensible-enthalpy-conserving temperature merge of
§14.3.3,

$$
\sum_{j} m_j\!\int_{T_{\rm ref}}^{T_{\rm merged}}\! c_{p,j}\,dT'
= H_a(T_a) + H_b(T_b)
\ \text{(bracketed Newton; every component, gas included)},\qquad
\xi_{{\rm merged},j} = \frac{V_a\,\xi_{a,j}+V_b\,\xi_{b,j}}{V_{\rm keep}},
$$

and parks the freed slot above the live surface
(`derivative_discontinuity!(integrator, true)`). The state keeps its size, so
the Jacobian's sparsity pattern and caches stay valid and `saveat`
interpolation is unaffected (§14.3.5–§14.3.6).

The **Newton-failure cap** (`create_newton_failure_cap_callback`) is a step cap with a memory of Newton failures. After
every accepted step it reads the integrator's count of failed nonlinear solves.
When the step just accepted needed retries, the proposed next step is held at or
below $g\,\Delta t_{\rm acc}$ with $g = 1.25$ (`FAILURE_MEMORY_GROWTH`) for the
next $M = 20$ accepted steps (`FAILURE_MEMORY_STEPS`); a further failure
restarts the memory from the new step. The error controller alone regrows the
step to its accuracy limit after a failure, and on a strongly nonlinear front
that limit lies beyond the range in which the stage Newton converges, so the
integrator settles into a fail–halve–regrow rhythm. On the Landau ablation case
V-CPL-02 the Newton fails near $\Delta t = 5$ ms and converges at 2.5 ms while
the controller keeps proposing 5 ms: the cap takes the case from 1087 steps and
892 Newton failures to 1032 steps and 78 failures. A run without Newton
failures never triggers the cap and keeps its trajectory exactly. It is
installed for every integrator that solves its stages by Newton (an
algorithm with an `nlsolve` field).

> **Implementation note.** Callbacks: src/Solver/callbacks.jl
> (`create_step_update_callback`, `create_thickness_termination_callback`,
> `create_depletion_callback`, `park!`, `create_newton_failure_cap_callback`)
> and src/Solver/solve.jl (assembly, diagnostics, verify-Jacobian). The
> cross-keyword checks are in `_validate_solve_kwargs`, src/Solver/solve.jl.

---

## 15.6 Jacobian backends

The production backend is **`Structured{S}`**, parameterized by a
`LinearSolveStrategy` $S$ (`DirectSolve` or `ApproxSolve`). It assembles the
Jacobian *operator by operator*, exploiting the block-tridiagonal physics
structure and representing the non-local terms as O(N) couplings until they are
written into the sparse matrix. A second backend, **`DenseAD_ForwardDiff`**,
runs `ForwardDiff` over the *entire* residual vector and returns a dense
Jacobian. It is correct by construction (it differentiates the same `residual!`
the integrator uses) but $O(N^2)$ in storage and slow; it is the
**verification reference** (§15.13), not a solve backend.

`solve` takes the bare backend: `jacobian = Structured(strategy =
DirectSolve())` or `Structured(strategy = ApproxSolve())`. With `jacobian =
nothing` (the default) it picks `DirectSolve` on a fixed mesh — exact assembly
plus a sparse direct solve — and `ApproxSolve` under ALE, which drops the
structured couplings (§15.10.2). `solve` builds the spec for its own mesh copy,
mode, and options; the `jacobian` keyword is typed `Union{Structured,Nothing}`,
so passing a prebuilt spec or `DenseAD_ForwardDiff()` raises a `TypeError`.

`JacobianSpec` bundles the backend with everything the integrator needs:

```julia
struct JacobianSpec{B<:JacobianBackend, S, C}
    backend::B
    sparsity_phys::S      # SparseMatrixCSC{Bool,Int}, or nothing (dense reference)
    sparsity_geom::S      # geometry (z-column) pattern; nothing for Eulerian
    cache::C              # backend-private workspace (not shareable)
end
```

It is built by `JacobianSpec(backend, def, ws)` — the constructor `solve` calls,
and the one the verifier and other advanced callers use with a `(def, ws)` pair
from `initialize_state`. A single entry point
`compute_jacobian!(J, spec, u, def, t, ws)` dispatches to
`_compute_jacobian!(J, spec.backend, …)`.

> **Implementation note.** Backend contract and `JacobianSpec`:
> src/Jacobian/interface.jl. Dense reference: src/Jacobian/dense_ad_reference.jl.
> Structured backend assembly: src/Jacobian/structured_backend.jl. The
> `is_approximate(backend)` trait (default `false`; `true` for
> `Structured{ApproxSolve}`) selects strict vs. masked verification (§15.13).

---

## 15.7 The structured-operator architecture

The structured Jacobian is organized into three orthogonal layers.

### 15.7.1 Layer 1 — Operators

Each physics block is an `Operator` subtype: a lightweight stateless tag with two
methods, its residual contribution `apply!` (a thin wrapper over the residual
kernel it differentiates) and its Jacobian contribution
`contribute_to_jacobian!(jac, op, u, t, ws, args...)`, where `ws` is the
operator's workspace from the backend cache and the remaining arguments are
operator-specific. There is no default fallback: each operator implements its
own contribution, as `ForwardDiff` over a stencil-local kernel (the
`GeometryOperator` differentiates the full residual with respect to the $z$
block), or as an analytical override (conduction, gas transport, Darcy transport)
tested against a ForwardDiff reference of the same kernel at `rtol = 1e-12`. Local entries are added into the
sparse block `jac.A`; non-local blocks are pushed onto `jac.structured` as
`StructuredCoupling`s. The ALE mesh-velocity operator emits its non-local
blocks through a separate `contribute_to_structured!`; the $\chi$ operator pushes
its rank-1 block from its own `contribute_to_jacobian!`.

The backend walks the operators in a fixed order:

```
1. (ALE)      ALEMeshVelocityOperator      apply! + strain-rate (and θ_A) derivatives,
                                           then contribute_to_structured! (z-row prefix sums, θ_A rank-1)
2. ReactionSourceOperator
3. ConductionOperator
4. GasTransportOperator | DarcyTransportOperator
5. RadiationOperator
6. BCThermalOperator
7. BCMassOperator
8. BCConvectionOperator                    (boundary-cell S_conv energy row)
9. (ALE)      dilation rows                (cell-local, from the cached ∂θ)
10.(ALE & :full) GeometryOperator          (full-residual AD over z-columns)
11.(WithChi)  ChiOperator                  (cell-local rows; χ-column rank-1)
```

The order is fixed so that the data each step reads exists: the state cache,
the ALE geometry, and the properties are refreshed first (mirroring the
residual), the mesh-velocity workspace and its strain-rate derivatives are filled
before the dilation rows read them, and the $\chi$ operator reads the
cross-section derivatives the mesh-velocity workspace holds. The radiation,
mesh-velocity and $\chi$ operators build their embedded couplings at the global
row/column offsets of the blocks they fill; `compose_into!` embeds the prefix-sum
factor vectors at the node rows.

> **Implementation note.** Operator interface: src/Jacobian/operators/Operator.jl;
> per-operator implementations under src/Jacobian/operators/; the walk is
> `_compute_jacobian!` in src/Jacobian/structured_backend.jl. Per-operator
> `ForwardDiff` chunk sizes: reaction `1+N_C`; conduction $3(1+N_C)$; gas
> transport $3(1+N_C)$; BC thermal $1+N_C$; BC mass $1+N_C$; BC convection
> $2(1+N_C)$; $\chi$ $1+N_C$. The ALE dilation derivatives come from the
> mesh-velocity operator's strain-rate sweep.

### 15.7.2 Layer 2 — `StructuredJacobian`

The assembled Jacobian is *not* a monolithic sparse matrix. It is

```julia
struct StructuredJacobian
    A::SparseMatrixCSC{Float64,Int}        # local block: block-tridiag physics + sparse augments
    structured::Vector{StructuredCoupling} # non-local couplings in factored form
    smap::FastSparsityMap                  # O(1) sparse-write helper
    n::Int                                 # full state dimension
end
```

The split is the architectural heart of the design: the sparse block
$\mathsf{A}$ carries everything *local* (block-tridiagonal mass / conduction /
reaction / transport plus the cell-local self-terms), and the `structured`
vector carries everything *non-local* (the cumulative ALE mesh velocity, the
Beer–Lambert optical-depth suffix, the column-global $\theta_A$ and $\chi$
back-reactions) in factored form. The assembly writes them into the sparse
Jacobian with `densify!` only under the exact `:full` profile.

`FastSparsityMap` is an O(1) sparse-write helper: a dense lookup table
`indices[row + (col-1)*n_rows] = nz_idx` (zero meaning "not in pattern"). It
trades quadratic memory ($\sim2.5\,\mathrm{MB}$ at 560 DOFs, $\sim90\,\mathrm{MB}$ at
3360 DOFs) for O(1) scattered additions via `fast_sparse_add!`.

### 15.7.3 Layer 3 — `LinearSolveStrategy`

A strategy selects which couplings the assembled Jacobian carries; it is
orthogonal to the assembling backend. `DirectSolve` writes every structured
coupling into the sparse matrix (the exact Jacobian), and `ApproxSolve` drops
them (the reduced pattern). The integrator then factors the resulting sparse $W$
with KLU or UMFPACK (§15.5.1). The two strategies are detailed in §15.10.

---

## 15.8 Per-operator Jacobian formulas

We give the exact discrete derivative each operator contributes, with the sign
conventions verified against the code. Throughout, the residual writes
$dT_i\mathrel{+}=(\text{source})/\rho c_p^{\rm eff}$ and the conduction
divergence enters with a leading minus, so signs are tracked *per operator* and
accumulated additively.

### 15.8.1 Reaction source

The cell-local contribution is

$$
\frac{dT_i}{dt}\bigg|_{\rm rxn}=\frac{Q_{\rm rxn}(T_i,\xi_i)}{\rho c_p^{\rm eff}(T_i,\xi_i)},\qquad
\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn}=S_{j,{\rm rxn}}(T_i,\xi_i),
$$
$$
Q_{\rm rxn}=-\sum_r h_r\,r_r,\qquad S_{j,{\rm rxn}}=\sum_r \nu_{r,j}\,r_r.
$$

The Jacobian is purely cell-local — the $(1+N_C)\times(1+N_C)$ block coupling
$(T_i,\xi_{1,i},\dots,\xi_{N_C,i})$ to itself — obtained by `ForwardDiff` over
`apply_cell!` (which evaluates the reactions, the heat source *including the
$1/\rho c_p$ division*, and the species sources in place). Note the energy row
differentiates the *quotient* $Q_{\rm rxn}/\rho c_p^{\rm eff}$, so it captures
both the explicit $T,\xi$ dependence of the rates *and* the $\xi$-dependence of
the heat capacity. The smooth tanh temperature gates and the multiplicative
$\tanh(\xi/\text{threshold})$ depletion rate factor (§6) make these
derivatives well-defined and free of the spurious spikes a hard Heaviside or
$\min$ would produce — this is the reason the kinetics are AD-clean.

> **Implementation note.** `ReactionSourceOperator`, `apply_cell!`,
> src/Jacobian/operators/reaction_source.jl. Chunk size $1+N_C$.

### 15.8.2 Conduction (tridiagonal divergence)

The energy row depends on the three cells $\{i-1,i,i+1\}$ through the divergence
of the harmonic-mean face flux. With the distance-weighted harmonic-mean
conductivity (§7.1)

$$
q_f=-k_f\frac{T_R-T_L}{d_L+d_R},\qquad k_f=\frac{d_L+d_R}{d_L/k_L+d_R/k_R},
$$

the operator computes (verified in src/Jacobian/operators/conduction.jl)

$$
q_{\rm right}=-k_{f,R}\frac{T_R-T_C}{d_{CR}},\quad
q_{\rm left}=-k_{f,L}\frac{T_C-T_L}{d_{LC}},
$$
$$
\frac{dT_i}{dt}\bigg|_{\rm cond}=-\frac{(q_{\rm right}A_{\rm right}-q_{\rm left}A_{\rm left})/V_C}{\rho c_p^{\rm eff}(T_C,\xi_C)}.
$$

The leading minus realizes "outflow on the right cools the cell," matching the
discrete divergence $\mathrm{div}[i]=(q_R-q_L)A/V_i$ and the residual's
$-\mathrm{div}[i]/\rho c_p$. In ALE with $\chi$-shrinkage, $A_{\rm left}$ and
$A_{\rm right}$ can differ slightly (the face areas carry the geometry), so the
code keeps them distinct rather than factoring out a single $A$. The Jacobian is
the tridiagonal block, assembled by a hand-coded analytical override
(`_contribute_to_jacobian_analytical!`) that is tested against a `ForwardDiff`
reference over the 3-cell stencil (`_contribute_to_jacobian_ad!`, chunk size
$3(1+N_C)$) at `rtol = 1e-12`; each of the three cells contributes
$(T,\xi_1,\dots,\xi_{N_C})$ inputs because $k_f$ and $\rho c_p^{\rm eff}$ both
depend on composition.

**Boundary-cell self-fill.** At a boundary cell (cell 1 has no left neighbor,
cell $n$ no right), the missing neighbor's $(T,\xi)$ in the AD stencil is filled
with the cell's *own* value, and the missing face flux is zero (the BC operator
owns the boundary flux). AD over a zero face produces structural-zero columns,
which the sparse scatter respects. The genuine boundary flux contribution is the
`BCThermalOperator`'s responsibility, computed through the surface-temperature
Newton solve (§12.2).

> **Implementation note.** `ConductionOperator`, `apply_cell!`, `_pack_stencil!`,
> src/Jacobian/operators/conduction.jl.

### 15.8.3 Gas transport / Darcy

The Fickian gas-transport operator differentiates the volume-fraction-gradient
flux (§9.1) over a $\le\pm1$-cell stencil (chunk $3(1+N_C)$), giving the direct
derivatives $\partial J/\partial\lambda_f$, $\partial J/\partial\xi_L$,
$\partial J/\partial\xi_R$ and their $T$-coupling through $T_{\rm face}$. The
Darcy operator (`DarcyTransportOperator`) additionally carries the
pressure-driven advection $(\xi_{\rm upwind}/\varphi_{\rm face})\,v_g$: it
reconstructs the face velocity inside its own stencil closure from the
cell-local state — $P_i$ from the ideal-gas closure, $\kappa_i$ from the mixing
rule, $\mu_i$ from the material's viscosity curve, the face mobility by the residual's
distance-weighted harmonic rule (`_dt_mobility_face`, zero when either side
$\le10^{-30}$) — so the chain $\partial J/\partial\xi$ and
$\partial J/\partial T$ through $P$, $\varphi$, $\kappa$, $\mu$ is
differentiated through the same formulas the residual executes (including the
$\xi/\varphi_{\rm face}^2$ porosity chain, zero when the $\varphi$ floor binds).
First-order upwinding keeps the advective stencil local; the diffusive part uses
harmonic-mean $\lambda_f$ as in the residual. Boundary faces are zeroed here —
the boundary venting flux is `BCMassOperator`'s row, and the boundary cell's
whole $S_{\rm conv}$ energy row (whose gradient runs through the surface-Newton
$T_s$) is `BCConvectionOperator`'s (§15.8.8); the transport operators keep the
boundary cell's ξ-rows but contribute no boundary-cell $S_{\rm conv}$. Each
operator's `contribute_to_jacobian!` is a hand-coded analytical override, gated
against a ForwardDiff reference of the same kernel (`_contribute_to_jacobian_ad!`)
at `rtol = 1e-12`.

> **Implementation note.** `GasTransportOperator`,
> src/Jacobian/operators/gas_transport.jl; `DarcyTransportOperator`,
> src/Jacobian/operators/darcy_transport.jl (`_dt_props` Dual pass,
> `_dt_mobility_face`, `_dt_J_face_partials`, analytical override
> `_contribute_to_jacobian_analytical!`).

### 15.8.4 ALE mesh velocity and the dilation correction

The mesh velocity operator does **not** AD a single residual; it builds the
analytical chain-rule pieces and emits a structured coupling. Two derivative
sweeps populate workspace buffers:

- `compute_strain_rate_derivatives!` — cell-local `ForwardDiff` over
  $\theta(T,\xi)=\big(\sum_j\gamma_j S_j(T,\xi)\big)/\rho(T)$, producing
  $d\theta/dT_i$ and $d\theta/d\xi_{j,i}$.
- `compute_theta_A_derivatives!` (WithChi) — analytical $d\theta_A/d\chi_j$ and
  $dA/d\chi_j$.

The cumulative velocity $w_k=\sum_{i<k}\theta_{L,i}\Delta z_i$ has the
*exclusive prefix-sum* structure: $\partial w_k/\partial x_i=\Delta z_i\,
\partial\theta_i/\partial x_i$ for $i<k$ (and $0$ for $i\ge k$). The per-cell
scaling $v[i]=\Delta z_i\,\partial\theta_i/\partial x_i$ is held in a
`PrefixSumCoupling{Float64}` (a plain carrier of the per-cell factors) per source
channel ($T$ and each $\xi_j$), which `compose_into!` embeds at the z-row offsets
as `EmbeddedNodePrefixSumCoupling` blocks. The z-row equations $dz_k/dt = w_k$
are themselves consumers of the cumulative coupling, with $\partial(dz_k/dt)/\partial x_j =
\mathbb{1}[j\le k-1]\,v[j]$. The geometry AD pass (§15.11.3) covers z
*columns* only and plays no part in these rows; the z-row × $(T,\xi)$ block
is owned entirely by these embedded couplings.

Under `WithChi`, the column-global $\theta_A$ enters every
$w_k = \sum_{i<k}(\theta_i-\theta_A)\Delta z_i$, so each source channel
$y\in\{\chi,\,T,\,\xi_1..\xi_{N_C}\}$ contributes a genuinely rank-1 block
$\partial w_k/\partial y_j = -\text{cumZ}_{\rm node}[k]\cdot
\partial\theta_A/\partial y_j$ with

$$
\frac{\partial\theta_A}{\partial\chi_j} = G'\,\frac{m_{{\rm init},j}}{M}\,
\frac{d\bar\chi}{dt} + G\,\frac{\partial A/\partial\chi_j}{M}\,
\frac{\sum_i S_{{\rm dry},i}V_i}{A},\qquad
\frac{\partial\theta_A}{\partial T_j} = \frac{G\,V_j}{M}\,
\frac{\partial S_{{\rm dry},j}}{\partial T_j},\qquad
\frac{\partial\theta_A}{\partial \xi_{k,j}} = \frac{G\,V_j}{M}\,
\frac{\partial S_{{\rm dry},j}}{\partial \xi_{k,j}},
$$

where $G=\mathrm{law}'(\bar\chi)/\mathrm{law}(\bar\chi)$. The rank-1s (shared
left factor $u[k]=-\text{cumZ}_{\rm node}[k]$) are embedded directly at the
z-rows as `EmbeddedRank1Coupling` blocks.

The **dilation correction** $d\xi_{j,i}\mathrel{-}=\xi_{j,i}\theta_i$ adds the
cell-local Jacobian (every component $j$)

$$
\frac{\partial}{\partial T_i}:\ -\xi_{j,i}\frac{\partial\theta_i}{\partial T_i},\qquad
\frac{\partial}{\partial\xi_{k,i}}:\ \begin{cases}
-\theta_i-\xi_{j,i}\dfrac{\partial\theta_i}{\partial\xi_{j,i}} & k=j,\\[6pt]
-\xi_{j,i}\dfrac{\partial\theta_i}{\partial\xi_{k,i}} & k\ne j,
\end{cases}
$$

reusing the cached $d\theta/dT$, $d\theta/d\xi$ from the strain-rate sweep
(`_apply_ale_dilation_jacobian!`).

> **Implementation note.** `ALEMeshVelocityOperator`,
> `compute_strain_rate_derivatives!`, `compute_theta_A_derivatives!`,
> src/Jacobian/operators/ale_mesh_velocity.jl; dilation Jacobian in
> src/Jacobian/structured_backend.jl. A documented sharp edge: the
> $d\theta/dT$, $d\theta/d\xi$ buffers persist in the workspace and are reused by
> both the prefix-sum builder and the dilation Jacobian; assembly is atomic per
> Newton step, so the reuse is safe in practice.

### 15.8.5 No ALE advection operator

Every component rides the mesh and the thermal mass rides with it (§11.1.2,
§11.2.2), so neither the ξ rows nor the T row carry a mesh-motion coupling
through $w$ in production: the only mesh-motion term in the species rows is the
cell-local dilation of §15.8.4, and the mesh-relative temperature kernel runs
only under the `w_override` verification hook, which the structured Jacobian
does not support (override runs use finite-difference Jacobians). The cumulative
mesh velocity therefore couples only the z rows, through the prefix-sum blocks
the mesh-velocity operator embeds itself. The block-tridiagonal physics
pattern is identical under the `:full` and `:ale_reduced` profiles; the profiles
differ in the geometry block and the structured couplings (§15.11.2).

### 15.8.6 Variable cross-section ($\chi$)

Under `WithChi` (and only for non-identity laws; identity-law materials have
$d\chi/dt\equiv0$), the $\chi$ row is

$$
\frac{d\chi_i}{dt}=\frac{S_{{\rm dry},i}(T_i,\xi_i)\,V_i}{m_{{\rm init},i}}.
$$

Its Jacobian has three parts:

- **Cell-local $(T,\xi)$**: `ForwardDiff` over $S_{\rm dry}\,V/m_{\rm init}$ gives
  $\partial(d\chi_i/dt)/\partial T_i$ and $\partial(d\chi_i/dt)/\partial\xi_{j,i}$.
- **$z$-columns (ALE)**: because $V_i=A|z_{i+1}-z_i|$ is linear in the two
  bounding node positions, the entries are
  $\partial(d\chi_i/dt)/\partial z_i=-A\,S_{\rm dry}/m_{\rm init}$ and
  $\partial(d\chi_i/dt)/\partial z_{i+1}=+A\,S_{\rm dry}/m_{\rm init}$. The
  `ChiOperator` does not write them: the `GeometryOperator`'s AD pass over the
  $z$ subvector fills them under `:full` (§15.11.3), and `:ale_reduced` drops them.
- **$\chi$-columns (rank-1)**: the area $A=A_0\,\mathrm{law}(\bar\chi)$ depends on
  *every* $\chi_j$ through $\bar\chi$, giving a dense column coupling factored as
  an `EmbeddedRank1Coupling` with
  $u[i]=S_{{\rm dry},i}\,\Delta z_i/m_{{\rm init},i}$ ($\Delta z_i=V_i/A$) and
  $v[j]=dA/d\chi_j$, pushed only under `:full`.

> **Implementation note.** `ChiOperator`, `apply!`, `_chi_residual_kernel!`,
> src/Jacobian/operators/chi.jl. Chunk size $1+N_C$. The $\chi$ block is appended
> after $(T,\xi,z)$; the offset is $\chi_{\rm off}=n(1+N_C)+\delta_{\rm ALE}n_{\rm nodes}$.

### 15.8.7 Beer–Lambert non-local radiation

In `BeerLambert` mode the absorbed flux in cell $i$ depends on the optical depth
of *all cells above it*,
$q_{\rm rad}[i]=I_{\rm in}\,(1-e^{-\tau_i})/\Delta z_i$ with
$\tau_i=\alpha_i\Delta z_i$ and $I_{\rm in}$ the cumulative-attenuated intensity
from the top. Differentiating (verified in src/Jacobian/operators/radiation.jl):

$$
\frac{\partial q_{\rm rad}[i]}{\partial x_j}=
\begin{cases}
I_{\rm in}\,e^{-\tau_i}\,\dfrac{\partial\alpha_i}{\partial x_j}=\big(I_{\rm in}-q_{\rm rad}[i]\,\Delta z_i\big)\dfrac{\partial\alpha_i}{\partial x_j} & j=i,\\[10pt]
-\,q_{\rm rad}[i]\,\Delta z_j\,\dfrac{\partial\alpha_j}{\partial x_j} & j>i,\\[6pt]
0 & j<i.
\end{cases}
$$

The energy row $T[i]$ therefore couples to the $\xi$-columns of every cell
$j\ge i$ — an **upper-triangular** (strictly, a *suffix*) pattern. The diagonal
$j=i$ term is added by the operator through a cell-local `ForwardDiff` pass (with
$I_{\rm in}$ held fixed); the strict-suffix $j>i$
couplings are captured as an `EmbeddedScaledSuffixSumCoupling`. The split
($j=i$ cell-local, $j>i$ structured) avoids double-counting the diagonal.

On the Kirchhoff default path (`radiation_kirchhoff = true`, §8.3.2) the
injected intensity is $a_{\rm eff}(T_n,\xi_{\cdot,n})\,I_{\rm surface}$ (the
effective surface absorptivity of §5.8.2, equal to $\varepsilon_{\rm eff}$ unless a
component sets its own `absorptivity`), so every energy row $i<n$ additionally
couples to the **surface cell**: the suffix vector carries
$v_n \mathrel{-}= (\partial a_{\rm eff}/\partial x_n)/a_{\rm eff}$,
the surface cell's diagonal closure recomputes $a_{\rm eff}(T,\xi)$
live, and the sparsity pattern gains $T$-row($i$) × $T$-column($n$) slots for
$i<n$ (`_emit_radiation_pattern!`). With an explicit BC absorptivity the factor
is folded into $I_{\rm surface}$ upstream and these terms vanish.

Under ALE the operator fills **T and ξ columns only**. The dependence of the
deposit on the live geometry ($\Delta z_i = V_i/A$ from the z block) lands in the
z columns, which the `GeometryOperator`'s AD pass over the full residual writes
under the exact `:full` profile (§15.11.3); the `:ale_reduced` profile drops
every z column (src/Jacobian/operators/radiation.jl).

> **Implementation note.** `RadiationOperator`,
> src/Jacobian/operators/radiation.jl. Active only in `BeerLambert` mode;
> `SurfaceAbsorption`/`NoRadiation` leave a tridiagonal pattern.

### 15.8.8 Boundary conditions

`BCThermalOperator` overwrites the boundary-cell conductive flux with the
BC-specified surface flux (convective + radiative + fixed-$T$), obtained from the
surface-temperature Newton solve $f(T_s)=q_{\rm BC}(T_s)-k_{\rm eff}(T_s-T_{\rm cell})/\Delta z_{1/2}=0$
(§12.2). Its Jacobian is taken by AD through that Newton solve, producing the
boundary-cell row's dependence on $(T_{\rm cell},\xi_{\rm cell})$. The outward
gas mass flux fed to the Newton solve is rebuilt from the stencil state by the
shared boundary-face kernel `_bc_boundary_face_fluxes` (film law plus, in Darcy
mode, the venting flux) — the same flux the residual passes to the blowing models
— so the $(T_b,\xi_b)$ derivatives of a blowing model's coefficient $h\chi$ ($\chi$
the blowing factor of §12.3.3, not the pyrolysis progress) are exact
(§12.3.3). Chunk size $1+N_C$.

`BCMassOperator` owns the species rows of the boundary-face gas flux $F_g$
(film law plus, in Darcy mode, the pressure-BC advective venting),
$d\xi_{g,\rm bdy} \mathrel{+}= F_g A_f/V$, because the residual's
$\nabla\!\cdot\!J_{\rm gas}$ is formed *after* `apply_mass_bc_flux!` fills the
boundary faces, while the transport operators zero those faces. Everything is
cell-local on the boundary cell; its block is $N_C\times(1+N_C)$ and the chunk
size $1+N_C$. The shared per-component flux kernel is
`_bc_boundary_face_fluxes` (film law + venting, positive into the domain).

`BCConvectionOperator` is the **single owner of the boundary cell's
$S_{\rm conv}$ energy row**. At a boundary cell the residual's gradient
blends the interior one-sided slope with the boundary-face slope at the
solved surface temperature (§7.4.2), so the row is the product

$$
S_{\rm conv,b} = -\sum_g c_{p,g}(T_b)\,
\underbrace{\tfrac12\bigl(J_{g,f_{\rm bc}} + J_{g,f_{\rm int}}\bigr)}_{\bar J_g}
\,\underbrace{\frac{d_{\rm nb}\,s_{\rm face}(T_s) + \Delta z_{1/2}\,s_{\rm int}}
{d_{\rm nb}+\Delta z_{1/2}}}_{(\partial T/\partial z)_b},
$$

whose factors are owned by three different operators — the BC face flux
(mass-BC kernel), the interior face flux (transport kernel), and the $T_s$
implicit chain (thermal-BC Newton). A split assembly would need
$\bar J \times \partial(\partial T/\partial z)/\partial x$ cross-products
spanning owners, so the whole row lives here instead: the operator re-runs
the Dual-aware surface Newton (fed the outward flux from the shared
boundary-face kernel, so each convective member's blowing model is differentiated
exactly; $s_{\rm face} = \pm(T_s-T_b)/\Delta z_{1/2}$ is algebraically
identical to the residual's $-q_{\rm total}[f]/k_b$ because the applied flux
is exactly $\mp k_b(T_s-T_b)/\Delta z_{1/2}$), evaluates both face-flux
kernels mode-appropriately (Fickian `gas_flux` or Darcy-Fick
`gas_flux_darcy_fick` for the interior face; `_bc_boundary_face_fluxes` for
the BC face), and ADs the whole product in one piece. The transport
operators zero their boundary-cell $S_{\rm conv}$ (keeping the boundary
cell's ξ-rows). Stencil: the
boundary cell and its interior neighbor's $(T,\xi)$ — inside the physics
tridiagonal; chunk size $2(1+N_C)$
(src/Jacobian/operators/bc_convection.jl).

---

## 15.9 Structured couplings

The non-local terms are represented by typed coupling objects, subtypes of
`StructuredCoupling`, each carrying O(N) data — a per-cell factor vector and its
global row/column offsets — and implementing `densify!(M, c)`, the addition
$M \mathrel{+}= C$ that the `:full` profile uses to write the block into the
sparse Jacobian. All three shapes are *embedded*: they know where they sit in
the full state vector.

### 15.9.1 `EmbeddedNodePrefixSumCoupling` (exclusive prefix sum at the node rows)

The cumulative ALE mesh velocity: node $k$'s velocity depends on the strain rates
of all *lower* cells,

$$
M[\text{row}_k,\text{col}_j] \mathrel{+}= \mathbb{1}[j\le k-1]\,v[j],
\qquad v[j]=\Delta z_j\,\frac{\partial\theta_j}{\partial x_j},
$$

i.e. $C=L\,\mathrm{diag}(v)$ with $L$ strictly lower-triangular all-ones
(exclusive prefix sum). Node 1, the pinned bottom node, gets no entries. The
per-cell factors arrive in a `PrefixSumCoupling` — a plain carrier of the vector
$v$, not itself a `StructuredCoupling` — which `compose_into!` embeds at the node
rows.

### 15.9.2 `EmbeddedScaledSuffixSumCoupling` (strict suffix sum)

The mirror image, with a per-row scale:

$$
M[\text{row}_i,\text{col}_j] \mathrel{+}= \text{row\_scale}[i]\;\mathbb{1}[j\ge i+1]\;v[j],
$$

i.e. $C=\mathrm{diag}(\text{row\_scale})\,U\,\mathrm{diag}(v)$ with $U$ strictly
upper-triangular all-ones. This is the Beer–Lambert non-local absorption: cell
$i$'s absorbed flux depends on all cells *above* it ($j>i$), with
$\text{row\_scale}[i]=-q_{\rm rad}[i]/\rho c_{p,i}^{\rm eff}$ for the energy row
(§15.8.7). The diagonal $j=i$ is not part of the coupling; the radiation operator
adds it into $\mathsf{A}$ through a cell-local `ForwardDiff` pass.

### 15.9.3 `EmbeddedRank1Coupling`

A rank-1 block $C=u\,v^{\mathsf T}$ at given offsets:

$$
M[\text{row}_i,\text{col}_j] \mathrel{+}= u[i]\,v[j].
$$

This is the column-global $\theta_A$ / $\chi$ back-reaction in `WithChi`: the
cross-section $A=A_0\,\mathrm{law}(\bar\chi)$ depends on *all* $\chi_j$, so the
affected rows get a dense column dependence captured as a single outer product
(the $\theta_A$ channels at the z rows, §15.8.4, and the $\chi$-column block of
the $\chi$ rows, §15.8.6).

> **Implementation note.** `PrefixSumCoupling`, the three embedded coupling types,
> their `densify!`, and `compose_into!` are in src/Jacobian/structured_couplings.jl.
> The couplings alias their factor vectors (no copies) and are rebuilt at every
> Jacobian evaluation; `zero!` empties `jac.structured` at the start of each
> assembly.

---

## 15.10 Linear-solve strategies

A strategy decides which couplings the assembled Jacobian carries. Both operate on
the same assembled `StructuredJacobian`; the integrator then factors the sparse
$W$ with the factorization of §15.5.1.

### 15.10.1 `DirectSolve` (anchor)

Write every structured coupling into the sparse matrix:
$J=\mathsf{A}+\sum_j\mathrm{densify}(C_j)$ under the `:full` sparsity profile.
The linear solve then sees the exact Jacobian. It is the verification anchor and
the default on a fixed mesh. Under ALE the prefix-sum and geometry blocks couple
every node to every cell below it, which makes the factorized $W$ effectively
dense — the reason for the other strategy.

### 15.10.2 `ApproxSolve` (Schur-style, drop couplings)

Keep the local block only, $J=\mathsf{A}$, under the `:ale_reduced` sparsity
profile: every structured coupling is omitted at its producer and the geometry
block is not assembled. Newton then sees a Jacobian that omits the dense ALE
block, while each linear solve stays local and cheap. `ApproxSolve` is the
**only** strategy carrying the `is_approximate` trait, which switches the
verifier to *masked* comparison (in-support only) since the dropped couplings are
*expected* omissions.

It is the default under ALE. With $W = I-\gamma\Delta t\,J$ the modified
Newton iteration that factors the reduced matrix has the error operator
$E = I - W_{\rm red}^{-1}W_{\rm full}$, whose spectral radius is its
asymptotic contraction per iteration. It was estimated at every accepted step
of KenCarp4 runs, at $\gamma$ times the larger of the step taken and the step
proposed, by power iteration in the integrator's weighted norm and checked
against dense eigenvalues, on plain-ALE decks without parking at 100 cells:
PMMA at 25, 50 and 100 kW/m² (solid $\lambda = 10^{-6}$, and 0 at 50 kW/m²)
run to 5 % of its thickness; a charring wood (Darcy, Beer–Lambert) whose char
shrinks the cells to 0.83 of their volume under a lateral shrinkage law; an
intumescent coating that swells 3.6× in volume under a widening lateral law;
the Landau ablation deck of V-CPL-02; and 80-s cone decks of PMMA and charring
wood. It never exceeded 0.19, and it peaked before any step collapse, so the
iteration contracts by at least a factor of five per iteration. $E$ is
strongly non-normal: because every $z$ column is dropped, its one-iteration
weighted norm reaches $10^2$ on the cone decks and far more once cells
collapse, with no visible cost in steps or Newton failures.

The two strategies follow the same trajectories:

- The decks that run to their end take the same steps to within two (54/54,
  979/979, 82/83, 80/80, 337/336, 286/288 reduced/exact) and agree to
  $3\times10^{-4}$ K in surface temperature and $5\times10^{-6}$ of the peak
  in gas generation. PMMA at 25 kW/m² reaches 5 % of its thickness at the same
  time, 1617.5 s, in 21,036 and 22,607 steps.
- Without parking, surface cells burn down to $10^{-10}$ to $10^{-13}$ m. On
  the faster burns the step then collapses, and both strategies stall alike:
  at equal step counts they reach the same simulated time. PMMA at 50 kW/m²
  with solid $\lambda = 0$ ends in `DtLessThanMin` at 225 s (reduced) and
  232.5 s (exact), the ablation deck at 0.25 s with both, and PMMA at 50 and
  100 kW/m² had not reached 5 % of its thickness after an hour of either. In
  that regime the surface mass flux changes sign between outputs with either
  Jacobian, while the gas generation and the thickness stay smooth and agree
  to $10^{-7}$. Parking (§14.3) is what carries a run through burnout; the
  exact Jacobian does not.

Per unit of work the reduced strategy with KLU is 8–17× faster than
`DirectSolve` with UMFPACK on these decks (median CPU time over five
interleaved rounds: 0.057–0.121 of it where the decks run to their end, and
0.107–0.120 for the stalling decks over an equal budget of 4000 integrator
iterations), and allocates 20–8000 times less.

| Strategy | Jacobian assembled | Exact? | Sparsity profile | Verification |
|---|---|---|---|---|
| `DirectSolve` | $\mathsf{A}+\sum C_j$ densified | yes | `:full` | strict |
| `ApproxSolve` | $\mathsf{A}$ only | no | `:ale_reduced` | masked |

> **Implementation note.** Strategy types: src/Jacobian/strategies.jl;
> `sparsity_profile(::Strategy)` maps `DirectSolve` → `:full` and `ApproxSolve` →
> `:ale_reduced`. The emit step `_emit_jacobian!` (src/Jacobian/structured_backend.jl)
> copies $\mathsf{A}$ into the integrator's matrix and densifies the couplings only
> under `:full`.

---

## 15.11 Sparsity construction

### 15.11.1 Two patterns: physics and geometry

The ALE geometry coupling motivates a **two-pattern** model:

- **$S_{\rm phys}$** — every coupling with $T$/$\xi$/$\chi$ columns:
  block-tridiagonal physics ($T\!\leftrightarrow\!T$, $T\!\leftrightarrow\!\xi$,
  $\xi\!\leftrightarrow\!\xi$ at each cell and its neighbors), the Beer–Lambert
  suffix (row $i$ depends on $\xi$-columns $j\ge i$, plus — on the Kirchhoff
  default path — the surface temperature column $T_{\rm col}(n)$ for every row
  $i<n$, the $a_{\rm eff}(T_n,\xi_{\cdot,n})$ transmissivity chain of
  §15.8.7), and the $\chi$-row block
  (cell-local plus the dense $\chi$-column). The ALE dilation is cell-local and
  lies inside the block-tridiagonal physics.
- **$S_{\rm geom}$** — every coupling with $z$ columns (ALE only; empty for
  Eulerian): physics$\to z$ (cells read nodes; reach $\le i+3$ without $\chi$,
  full $1{:}n_{\rm nodes}$ with `WithChi` via the $\theta_A$ path, and full for
  the energy rows under Beer–Lambert, whose optical depth
  $\tau_i=\sum_{k>i}\alpha_k\Delta z_k$ depends on every node above the cell),
  $z\to$physics (node $k$ depends on cells $i<k$, lower-triangular), and $z\to z$
  (lower-triangular including the diagonal, the cumulative nature of mesh
  velocity).

Four emitter functions compose these patterns. `build_sparsity(layout, mode;
profile)` builds $S_{\rm phys}$ from the first three, in this order, and
$S_{\rm geom}$ from the fourth:

1. `_emit_physics_pattern!` — block-tridiagonal;
2. `_emit_radiation_pattern!` — Beer–Lambert suffix (BeerLambert);
3. `_emit_chi_pattern!` — $\chi$-row and dense $\chi$-column (WithChi);
4. `_emit_radiation_pattern!` — Beer–Lambert suffix (BeerLambert).

The integrator's `jac_prototype` is the union $S_{\rm phys}\cup S_{\rm geom}$.
Every row of the union stores its diagonal in every profile, so the shifted
matrix $W$ of §15.5.2 has exactly the pattern of $J$; where no operator writes
the diagonal (the node rows and the $\chi$ rows under `:ale_reduced`) it holds
a structural zero. Patterns are memoized by the fingerprint
$(\mathrm{typeof}(\text{mode}),N_C,n_{\rm cells},n_{\rm nodes},\text{profile})$ in
`_SPARSITY_CACHE`.

### 15.11.2 Profiles

- `:full` (default) — all entries assembled; structured couplings emitted; the
  geometry AD pass runs.
- `:ale_reduced` (`ApproxSolve` only) — drops every structured coupling (the
  z-row prefix-sum blocks, the $\theta_A$ rank-1
  channels, the Beer–Lambert suffix, and the $\chi$-column rank-1) and the
  geometry block, of which only the node rows' structural diagonal remains.
  The masked verifier expects exactly this set of omissions, and counts the
  stored but never-written diagonals of the node and $\chi$ rows
  (`structural_zero_diagonal`) among them.

Only `:full` couples cells an arbitrary distance apart, through the couplings
it materializes: the ALE geometry block and prefix sums, the Beer–Lambert
suffix, and the $\chi$-column rank-1. Each adds $O(N^2)$ entries; every other
pattern couples a cell only to its immediate neighbours.
`materializes_geometry_block(mode, profile)` (`:full` under ALE) selects the
default factorization (§15.5.1).

### 15.11.3 The geometry operator (z-column AD)

Because $z$ couples through many indirect paths (mesh velocity, dilation,
$\chi$-dependent area), the $z$-columns are filled by a *separate* full-residual
`ForwardDiff` pass with respect to the $z$ subvector (`GeometryOperator`,
`GeometryADCache`), run *after* the operator walk and skipped under
`:ale_reduced`. This decouples the $z$-column computation from the per-operator
framework, at the cost that its entries must be consistent with what the
operators scattered into the $T/\xi/\chi$ rows independently — consistency that
verification (§15.13) checks.

The geometry pass covers **columns only**: it says nothing about the z-*rows*'
dependence on $(T,\xi,\chi)$. Those entries are structured couplings emitted by
`ALEMeshVelocityOperator` (§15.8.4) — the `EmbeddedNodePrefixSumCoupling`
prefix-sum blocks and, under `WithChi`, the $\theta_A$ rank-1 channels.

> **Implementation note.** `build_sparsity` and the four emitters:
> src/Jacobian/sparsity.jl. `GeometryOperator`: src/Jacobian/operators/geometry.jl.
> `FastSparsityMap`: src/Jacobian/fast_sparsity_map.jl.

---

## 15.12 ILU-preconditioned Krylov solvers

For very large ALE problems where the direct factorization of $W$ dominates, two
iterative `LinearSolve.jl` factorizations are provided. They operate on the
*concrete* sparse $W=I/(\gamma_{\rm RK}\Delta t)-J$ matrix the integrator forms,
not on the `StructuredJacobian` (these are integrator-level linear solvers,
orthogonal to the §15.10 strategies):

- **`ILUGMRESFactorization(τ=0.1, gmres_rtol=1e-8, restart=30)`** — build an
  ILU($\tau$) preconditioner $P=\mathrm{ilu}(W,\tau)$ and solve $Wx=b$ by
  preconditioned GMRES. If ILU fails on a near-singular transient, it falls back
  to diagonal scaling $P=\mathrm{Diagonal}(1/\max(|{\rm diag}(W)|,10^{-12}))$.
- **`ILUBiCGSTABFactorization(τ=0.1, rtol=1e-8)`** — same preconditioner, BiCGSTAB
  iteration (fixed memory, no restart).

The drop tolerance default $\tau=0.1$ trades fill for preconditioner quality;
inexact solves can increase Newton rejections on Arrhenius transients, so the
iterative path is recommended only when the direct factorization is the
bottleneck.

**Linear-solver family.** Every `solve` attaches a sparse `jac_prototype`, which
pairs with a sparse direct default chosen from its pattern (`KLUFactorization` or
`UMFPACKFactorization`, §15.5.1). An explicit `linear_solver` — one of the
iterative solvers above, `SparspakFactorization()`, or any other LinearSolve
backend, such as `MKLPardisoFactorize()` (with `using LinearSolve` and the Pardiso
/ MKL packages loaded) — always wins.

> **Implementation note.** `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`,
> and the Krylov integration: src/Solver/ilu_gmres.jl.

---

## 15.13 Verification

Every backend is anchored to `DenseAD_ForwardDiff` — `ForwardDiff` over the whole
residual, which is correct by definition because it differentiates the same
`residual!` the integrator integrates. `verify_jacobian_at(spec, def, ws, u, t;
reference=DenseAD_ForwardDiff(), rtol, atol)` compares the active backend's
Jacobian against this reference at a state point and returns
`(passed, max_abs_dev, max_rel_dev, worst_index, …, J_active, J_ref)`.

Two comparison modes, selected by the `is_approximate` trait:

- **Strict / support-complete** (exact backends, `is_approximate==false`):
  compares the *full* Jacobian under a **per-entry** tolerance with a
  row-magnitude floor,

  $$
  |\Delta J_{ij}| \le \text{atol} + \text{rtol}\cdot
  \max\big(|J_{\rm ref}[i,j]|,\ \max_j |J_{\rm ref}[i,:]|\big).
  $$

  Pyrolysis Jacobians span ~10 orders of magnitude across rows (Arrhenius
  ξ-rows ~$10^6$ vs ALE z-rows ~$10^{-2}$), so a global max-$|J_{\rm ref}|$
  normalization would let a wholly missing small-magnitude block hide behind
  the matrix maximum, with pass/fail depending on problem size rather than
  correctness. Under the per-entry rule, a missing entry contributes to its
  own row's reference scale and cannot hide, while honest roundoff at true
  zeros of a stiff row stays within $\text{rtol}\cdot\text{row max}$. Any missing
  entry the reference has nonzero (beyond its per-entry tolerance) is a
  failure.
- **Masked** (`is_approximate==true`, i.e. `Structured{ApproxSolve}`): compares
  only inside the declared support, under the coarse global-scale rule
  ($\max|\Delta J| \le \text{atol}$ or
  $\max|\Delta J|/\max|J_{\rm ref}| \le \text{rtol}$); entries the strategy
  *intentionally* omits (the dropped ALE/radiation/$\chi$ couplings of
  `:ale_reduced`) are informational, not failures. The in-support values of an
  approximate backend legitimately deviate by the dropped coupling magnitude,
  so masked mode is a sanity check, not an exactness gate.

Verification can run as a setup-time check or as a callback (`verify_jacobian = (rtol=, atol=, sample_times=)`)
firing at chosen `sample_times` during a live solve. The callback *warns* on
deviation rather than aborting — a single transient outlier should not kill an
otherwise-good production run; the user reads `max_abs_dev`/`max_rel_dev` to
judge severity.

Depletion handling (§14.3) keeps the state size and the Jacobian's pattern and
caches, so `verify_jacobian` runs through park events; the parked states of every
mode are part of the verifier's own matrix
(`test/jacobian/test_park_jacobian.jl`).

> **Implementation note.** `verify_jacobian_at`, strict/masked compare:
> src/Jacobian/verify.jl; `is_approximate`: src/Jacobian/interface.jl;
> the verify callback: src/Solver/solve.jl.

---

## 15.14 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to ours on first use.

### 15.14.1 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) solves the condensed phase by a
**fully-implicit backward-Euler** time discretization with **Patankar-style
linearization** of the source terms and a **tridiagonal (TDMA / Thomas)** linear
solve, sweeping the coupled energy/species/gas/pressure equations to convergence
each step. Mapping: Gpyro's volume fractions $X_\alpha$ are our $v_j$; its
pre-exponential $Z$ is our $A_i$; its permeability $K$ is our $\kappa$; its
volumetric reaction rate $\dot\omega'''$ is our $r$; its mixture average
$\bar\rho$ corresponds to our $\rho=\sum_j\xi_j$.

Contrasts with Pyrolysis.jl:

- **Time scheme.** Gpyro is first-order backward-Euler with an outer
  source-linearization sweep; Pyrolysis.jl uses an adaptive higher-order SDIRK
  (KenCarp4) with a true Newton solve per stage. Both are L-stable and implicit;
  Pyrolysis.jl gains temporal order and step-size adaptivity, at the cost of a
  more elaborate $W=I/(\gamma_{\rm RK}\Delta t)-J$ inversion than Gpyro's single
  TDMA pass.
- **Jacobian.** Gpyro never forms a global Jacobian; its Patankar linearization
  produces a tridiagonal coefficient matrix per equation, solved by TDMA. The
  coupling between equations is handled by the outer iteration. Pyrolysis.jl
  forms the *full coupled* Jacobian as a structured operator (block-tridiagonal
  $\mathsf{A}$ plus O(N) couplings), which captures cross-equation coupling
  directly inside the Newton step. The block-tridiagonal core of our
  $\mathsf{A}$ is the multi-field analog of Gpyro's per-equation tridiagonal
  matrix; our prefix/suffix couplings have no Gpyro analog because Gpyro's
  fixed-grid control volumes do not carry a cumulative-mesh-velocity term.
- **Mesh.** Gpyro uses fixed (or prescribed variable-$\Delta z$) control volumes
  with a solid-fraction / $\chi_k$ mechanism for swell/shrink; Pyrolysis.jl moves
  the mesh by ALE and represents the cumulative volume-change kinematics as a
  prefix-sum Jacobian coupling (§15.9).

### 15.14.2 ThermaKin2Ds

ThermaKin2Ds (Stoliarov et al.) solves the 1D (and 2D) condensed phase with a
**Crank–Nicolson** scheme in the primary (through-thickness, here $x$/our $z$)
direction and explicit treatment in the transverse directions, on a **dynamically
adjusted mesh**. Mapping: ThermaKin's stoichiometric coefficients $\theta_i^j$
are our $\nu_{i,j}$ (we reserve $\theta$ for strain rate, overload resolution
A2); its component concentrations are our $\xi_j$.

Contrasts:

- **Time scheme.** Crank–Nicolson is second-order and A-stable (not L-stable),
  trapezoidal in time; KenCarp4 is fourth-order, L-stable, with embedded error
  control. L-stability matters for the stiff Arrhenius transients at ignition,
  where Crank–Nicolson can ring; the L-stable SDIRK damps the fast modes
  cleanly. ThermaKin's fixed-$\theta$ trapezoidal rule contrasts with our
  adaptive SDIRK $\gamma_{\rm RK}$.
- **Energy convention.** ThermaKin (Eq. 17) carries gas enthalpy by the
  convective term $c_g(\mathbf{J}_g\cdot\partial T/\partial x)$ — exactly our
  $S_{\rm conv}$ — with every component in the storage term, which is the
  energy-equation convention Pyrolysis.jl adopts. The sign convention
  for heat of reaction differs in *publication* ($h>0$ exothermic in
  ThermaKin/Gpyro vs. $h>0$ endothermic internally here), but the volumetric
  source $Q_{\rm rxn}=-\sum h_r r_r$ is consistent once the convention is mapped
  (overload H1).
- **Mesh.** ThermaKin's dynamic mesh adjustment is the closest analog to our ALE
  with depletion parking; both track a receding surface. Pyrolysis.jl's distinction is the
  *Jacobian-aware* treatment: the moving-mesh kinematics enter the Newton matrix
  through the geometry-AD $z$-columns and the prefix-sum couplings, so the
  implicit step "knows" the mesh is moving.

### 15.14.3 FDS

The FDS (McDermott et al.) solid-phase pyrolysis model integrates a 1D Fourier
conduction equation with Arrhenius source terms; the surface boundary couples to
the gas-phase CFD. FDS uses a node-based 1D solid grid and an implicit/explicit
split appropriate to its multiphysics coupling. Pyrolysis.jl is a *standalone*
condensed-phase solver: it does not solve a gas-phase flow field, but it provides
a far richer in-solid model (multi-component effective properties, four radiation
models, Darcy gas transport, ALE with variable cross-section). FDS carries no
in-depth gas, so its energy equation is the condensed-phase limit of ours; it
shares the Arrhenius kinetic form, and its
two-flux (Schuster–Schwarzschild) radiation maps onto our Beer–Lambert / P1
choices (§8). On the *temporal/Jacobian* axis specifically, FDS does not expose a
structured-operator Jacobian of the kind described here — the relevant comparison
is methodological (implicit conduction + explicit/implicit coupling) rather than
a direct mapping of Jacobian internals.

---

## 15.15 Worked summary of the discrete update

Collecting the pieces, the production residual writes, for each cell $i$ and
component $j$:

$$
\frac{dT_i}{dt}=\frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}\Big],
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i}-\underbrace{\xi_{j,i}\theta_i}_{\text{ALE, every component}},
$$
$$
\frac{dz_k}{dt}=w_k=\sum_{i<k}(\theta_i-\theta_A)\Delta z_i\quad(\text{ALE}),\qquad
\frac{d\chi_i}{dt}=\frac{S_{{\rm dry},i}V_i}{m_{{\rm dry},i}}\quad(\text{WithChi}).
$$

The Jacobian of this RHS, $J=\partial\,\mathrm{RHS}/\partial\boldsymbol{u}$, enters
the implicit-RK $W=I/(\gamma_{\rm RK}\Delta t)-J$ and is assembled as the
structured operator $J=\mathsf{A}+\sum_j C_j$, where $\mathsf{A}$ is the
block-tridiagonal local block (reaction, conduction, transport, BC, dilation,
$\chi$-local, $z$-columns) and the $C_j$ are the prefix-sum (mesh velocity),
suffix-sum (Beer–Lambert), and rank-1 ($\theta_A$/$\chi$) couplings, written into
the sparse matrix under `DirectSolve` and dropped under `ApproxSolve`. The linear
solve $W\delta\boldsymbol{u}=-r$ is a sparse factorization of that matrix, and the
whole construction is verified against the dense AD reference.

---

## 15.16 Limitations and open issues

- **Geometry operator under `:ale_reduced`.** The $z$-columns are dropped
  entirely under `ApproxSolve`. Whether *any* residual $z$-coupling survives the
  profile (i.e. whether the profile is guaranteed to leave exactly zero $z$-rows)
  is asserted by the masked verifier rather than proven; the design intent is
  zero, and verification is the safeguard.
- **Chunk-size tuning.** Per-operator `ForwardDiff` chunk sizes are hardcoded
  (e.g. $3(1+N_C)$ for conduction). They are the informationally-necessary
  stencil widths, not adaptively tuned for cache behavior; no auto-tuning is
  present.
- **WithChi × Darcy.** `WithChi` is tested mainly with Fickian transport; the
  interaction of the $\chi$ rank-1 coupling with the full Darcy operator set is
  not exhaustively exercised.
- **Adjoint sensitivity across park events.** A park is a discrete rewrite of
  the state (§14.3.4); adjoint sensitivity through such an event is
  unsupported, and forward sensitivities are valid only between events.
- **ILU drop tolerance.** The default $\tau=0.1$ for the ILU-Krylov solvers has
  not been swept for optimality across pyrolysis regimes; inexact solves can
  raise Newton rejections on Arrhenius transients.
- **Warm-start artifacts.** Reusing a prior final state as fresh initial data can
  perturb finite-difference sensitivities and continuation studies; its effect is
  unquantified (§17.7).
- **Diagnostics cadence.** The diagnostics snapshots are stride-based (every
  `diagnostics_stride` steps), not interpolated to output times, and record the
  sensible energy only; the ledger fields are NaN (§16.8.3).


---

# 16. Conservation Diagnostics

This chapter derives the runtime conservation diagnostics that Pyrolysis.jl uses
to audit its own numerics. We develop a per-component **mass ledger** and a
**sensible-energy ledger** that closes to integrator tolerance by construction
for constant composition and volume. For every accumulator we give the
continuous balance law, its exact discrete counterpart as the code assembles it,
the precise sign conventions, and the conditions under which the residual is
expected to vanish. The chapter closes by stating exactly which balance terms
are *not* yet accumulated, why the energy residual cannot close for time-varying
composition without them, and how the production diagnostics callback operates
in snapshot mode.

All symbols follow the unified nomenclature of §2 (Nomenclature). In this chapter
`i` is the **cell index** and `r` the **reaction index**; `j` is the **component
index** (1…$N_C$). The three binding sign conventions are in force throughout:
flux is positive in the $+z$ direction (bottom→top); the discrete divergence
$(F_R-F_L)A/V$ is positive for *outflow*; and the heat-of-reaction convention is
$h>0$ endothermic, $h<0$ exothermic, with the volumetric heat source
$Q_\mathrm{rxn}=-\sum_r h_r r_r$ (see §6, Reaction Kinetics, and §7, Heat
Conduction).

The implementing module is `Diagnostics`
(`src/Diagnostics/Diagnostics.jl`), comprising the conservation-check kernels
(`src/Diagnostics/conservation_checks.jl`) and the callback infrastructure
(`src/Diagnostics/callback.jl`).

---

## 16.1 Purpose and design philosophy

A pyrolysis solver couples a stiff reaction network, in-depth radiation, gas
transport, a moving (ALE) mesh, and surface-cell depletion, all advanced by an
implicit time integrator. With this many interacting mechanisms, a silent error
in flux or source assembly can masquerade as physics. The diagnostics subsystem
exists to make conservation **observable at runtime** rather than only in offline
post-processing, so that a broken term announces itself as a non-closing ledger.

The question the energy ledger answers is one of bookkeeping: do the
accumulated boundary, reaction, and gas-coupling terms add up to the change in
stored energy that the RHS actually evolves? It is designed to close to
integrator tolerance *by construction* — it is assembled from exactly the terms
the RHS computes. A failure is a bug in flux/source assembly or in the callback
plumbing, not a physics error.

The key design subtlety, stated plainly in the source documentation of
`EnergyConservationTracker`, is that this ledger audits *bookkeeping, not the
physical model*: because it is built from the terms the RHS actually assembles,
**it will close even if the RHS itself double-counts a mechanism.** Catching that
class of error is the job of the independent closure tests (the sealed-column and
gas-enthalpy verification cases of §18) and of the formulation itself (the surface
balance carries no gas-enthalpy term alongside $S_\mathrm{conv}$, §7.4.3).

> **Implementation note.** The two trackers and the callback live in module
> `Diagnostics`. They are pure functions of the *current mesh state*
> (`mesh.cell_states`, `mesh.cell_volumes`, `mesh.active_cell_mask`) plus a small
> set of mutable accumulators; they never read the raw state vector directly. The
> solver synchronizes the mesh from the integrator's `u` before each diagnostic
> snapshot via `update_mesh_from_state!` (`src/Solver/solve.jl`).

---

## 16.2 The continuous balance laws being audited

Before discretizing, we state the integral conservation laws that the ledgers
approximate. Let $\Omega(t)=[0,L(t)]\times A(t)$ be the (possibly moving,
variable-area) domain.

### 16.2.1 Mass

Integrating the component conservation law (§3.2) over the domain and applying
the divergence theorem, for each component $j$:

$$
\frac{d}{dt}\int_\Omega \xi_j \,dV
= -\int_{\partial\Omega} J_j\cdot \hat{\mathbf n}\,dA
\;+\; \int_\Omega \Big(\tfrac{\partial \xi_j}{\partial t}\Big)_\mathrm{rxn}\,dV .
$$

The reaction term integrates to zero *summed over all components* because the
single-reactant stoichiometry conserves mass exactly ($\sum_j \nu_{r,j}=0$, with
the reactant coefficient $-1$ balancing product yields summing to $1$; see §6.5).
Therefore the **total** mass changes only through the boundary efflux. For a
component-by-component statement we keep the reaction term explicit but rely on
build-time stoichiometric closure (§6.5.3) to guarantee $\sum_j S_j=0$. Under ALE
every component rides the mesh (§11.1.2), so $J_j$ — the flux relative to the
skeleton — is the only transport channel on a moving mesh too. With the two
boundaries (top $z=L$, bottom $z=0$) and the sign convention that $J_j>0$ is
$+z$ transport, the *outward* normal flux is $+J_j$ at the top and $-J_j$ at the
bottom. The diagnostics layer adopts the **into-domain** convention for its
boundary inputs (positive $J$ = inflow), under which the net efflux is
$-(J_{\mathrm{top},j}+J_{\mathrm{bot},j})A$ — see §16.3.

### 16.2.2 Energy

Multiplying the temperature equation (§7.3) by $\rho c_p^\mathrm{eff}\,V$ and
integrating, the sensible energy the equation evolves obeys

$$
\frac{d}{dt}\int_\Omega \rho c_p^\mathrm{eff}\,T\,dV
\;\approx\; \int_{\partial\Omega} q_\mathrm{cond}\,dA
+ \int_\Omega Q_\mathrm{rxn}\,dV
+ \int_\Omega S_\mathrm{conv}\,dV ,
$$

with $\rho c_p^\mathrm{eff}=\sum_j \xi_j c_{p,j}(T)$ the mixture heat capacity over
**every** component, gas included (§7.3.2), and $q_\mathrm{cond}$ the *conducted*
boundary flux. The surface balance has no gas-enthalpy term, so the conducted
flux equals the applied boundary flux: the gas crosses the surface at $T_s$ on
both sides (§1.4.3). The volumetric radiation source $\int_\Omega Q_\mathrm{rad}\,dV$ of
the temperature equation is left out (it is zero except under `BEER_LAMBERT`;
§16.11). This is the law the energy ledger
discretizes. The $\approx$ flags two terms missing from the storage
time-derivative — the composition-change term and the ALE volume-change term —
which are the subject of §16.7.

---

## 16.3 Mass conservation ledger

### 16.3.1 Domain mass

The total mass of component $j$ in the domain at the current state is the
volume-weighted sum of cell concentrations over the **active** cells only:

$$
m_j(t) = \sum_{i\in\mathrm{active}} \xi_{j,i}\,V_i .
$$

This is a direct first-order (midpoint) quadrature of $\int_\Omega \xi_j\,dV$:
$\xi_{j,i}$ is the cell-centered mass concentration [kg·m⁻³] and $V_i$ the cell
volume [m³] (`mesh.cell_volumes[i]`, which is time-dependent in ALE mode and
carries the variable cross-section). The active-cell mask
(`mesh.active_cell_mask[i]`) excludes parked cells (see §14, Adaptive Mesh
Refinement), so the ledger always reflects only the live solution.

> **Implementation note.** `compute_domain_mass`,
> `src/Diagnostics/conservation_checks.jl`. The loop iterates
> `eachindex(mesh.cells)` and `continue`s on inactive cells; it accumulates into
> a `MVector{NC,Float64}` and returns an `NTuple{NC,Float64}`. The mask test must
> be performed cell-by-cell — never `eachindex` over a compacted array — which is
> the AMR invariant noted in §14.

### 16.3.2 Boundary efflux accumulation

The boundary mass flux update integrates the surface fluxes over one substep.
With the diagnostics **into-domain** sign convention ($J>0$ = inflow), the mass
that *leaves* during a substep of size $\Delta t$ is

$$
\Delta m_{\mathrm{out},j} = -\big(J_{\mathrm{top},j} + J_{\mathrm{bot},j}\big)\,A\,\Delta t ,
$$

accumulated into a running `cumulative_efflux`:

$$
m_{\mathrm{efflux},j}(t) = \int_0^t -\big(J_{\mathrm{top},j}+J_{\mathrm{bot},j}\big)\,A\,d\tau
\;\;\approx\;\; \sum_{\text{substeps}} \Delta m_{\mathrm{out},j}.
$$

A positive flux (inflow) contributes *negatively* to the efflux (mass stays in
the domain); a negative flux (the typical case — gas escaping at the top,
$J_{\mathrm{top},j}<0$) contributes *positively*. The accumulator therefore
counts up the total mass that has left. The canonical configuration is
$J_{\mathrm{top},j}<0$ (gas escaping) and $J_{\mathrm{bot},j}=0$ (impermeable
back face), giving a monotonically increasing efflux for the gaseous products.

> **Implementation note.** `update_mass_tracker!`,
> `conservation_checks.jl`. `J_boundary` is an `NC×2` matrix with column 1
> = top (boundary id 1), column 2 = bottom (boundary id 2); the code computes
> `mass_out = -(J_top + J_bottom) * A * dt` per component and adds it to
> `cumulative_efflux`. It also advances `tracker.last_time += dt`.

**Variable cross-section (ALE).** For a moving mesh the area $A(t)$ changes within
the substep. To keep the boundary contribution consistent with the trapezoidal
geometry update used elsewhere (§11, ALE Framework), the caller passes the
**trapezoidal-in-area average** over the substep,

$$
A = \tfrac12\big(A_t^\mathrm{pre} + A_t^\mathrm{post}\big),
$$

so that $\Delta m_{\mathrm{out},j}$ uses the substep-mean area. With the identity
lateral-shrinkage law (`NoChi`) or a non-ALE (Eulerian) run, $A$ is constant and
equals `mesh.cross_section_area_initial`.

### 16.3.3 The balance check and tolerance

Mass conservation is the statement that the current mass equals the initial mass
minus everything that has escaped:

$$
m_j(t) \;\overset{?}{=}\; m_{\mathrm{initial},j} - m_{\mathrm{efflux},j}(t)
\quad\Longleftrightarrow\quad
m_j(t) + m_{\mathrm{efflux},j}(t) = m_{\mathrm{initial},j}.
$$

The check defines the **expected mass** $m_{\mathrm{expected},j}=m_{\mathrm{initial},j}-m_{\mathrm{efflux},j}$
and a per-component relative error

$$
\epsilon_{m,j} =
\begin{cases}
\dfrac{\big|m_j(t) - m_{\mathrm{expected},j}\big|}{\,|m_{\mathrm{initial},j}|\,},
& |m_{\mathrm{initial},j}| > \texttt{eps()},\\[2.2ex]
\big|m_j(t) - m_{\mathrm{expected},j}\big|, & \text{otherwise (absolute fallback).}
\end{cases}
$$

The **zero-initial-mass fallback** is essential: gaseous products typically start
with $m_{\mathrm{initial},j}\approx 0$, so a relative error would divide by zero.
Switching to absolute error correctly flags even a tiny *spurious* mass when
there was no baseline to normalize against. The verdict is

$$
\text{conserved} \iff \max_j \epsilon_{m,j} < \texttt{tol},\qquad \texttt{tol}=10^{-10}.
$$

This is a tight, *bookkeeping-level* tolerance: it should hold to roundoff
because both sides are computed by exactly consistent quadratures. The check
returns the full per-component breakdown (`current_mass`, `expected_mass`,
`relative_error`, `max_error`, `conserved`).

**Scope of the $10^{-10}$ tolerance.** The roundoff-level closure holds wherever
the discrete ledger telescopes exactly: every component on a fixed (Eulerian)
mesh, and every component under ALE, where the only mesh-motion term is the
dilation $-\xi_j\theta_i$ (which conserves each cell's mass $\xi_jV$) and the
transport flux $J_j$ at the boundary faces is the whole of the efflux (§11.2.3).
A sealed collapsing column conserves each component to integrator tolerance
(V-CPL-03 P4).

> **Implementation note.** `check_mass_conservation` and
> `mass_conservation_report`, `conservation_checks.jl`; the latter renders the
> per-component table and an overall PASSED/FAILED line. The state is held in
> `MassConservationTracker{NC}`: `initial_mass`, `cumulative_efflux`,
> `last_time`.

### 16.3.4 Sources of (small) mass-ledger discrepancy

Even when the implementation is correct the ledger carries avoidable error
sources. (i) **Offline vs. online integration.** If the efflux is reconstructed
offline by trapezoidal integration of saved surface fluxes (rather than
accumulated online per substep), the integration error of the post-processing
quadrature appears in the residual. (ii) **Substep granularity.** The online
accumulator uses the substep-mean area; on a rapidly shrinking ALE mesh the
trapezoidal-area approximation is second-order in $\Delta t$. (iii) **Active-mask
transitions.** When a surface cell is parked (§14.3), its mass is merged into
its neighbour exactly; the ledger picks up the post-merge masses, so a correctly
conserving merge leaves no residual but a buggy one would surface here.

---

## 16.4 Energy storage

The stored quantity is the sensible energy the temperature equation evolves,

$$
E(t) = \sum_{i\in\mathrm{active}} \rho c_p^\mathrm{eff}(\xi_i, T_i)\, T_i\, V_i,
\qquad
\rho c_p^\mathrm{eff}(\xi,T)=\sum_{j=1}^{N_C}\xi_j\,c_{p,j}(T),
$$

with the mixture heat capacity over every component (`@generated` summation
kernel in `effective_properties.jl`, no phase branch). The sum respects the
active-cell mask. The product $\rho c_p\cdot T$ uses the **absolute** temperature
$T$ (not $T-T_\mathrm{ref}$): the diagnostics measure a sensible-energy *content*
relative to $0\,\mathrm{K}$, which is acceptable because the ledger only ever uses
*differences* $E(t)-E(0)$ in which the datum cancels (provided $\rho c_p$ is held
fixed; see the composition-change caveat in §16.7).

> **Implementation note.** `compute_sensible_energy`,
> `conservation_checks.jl`, calling `Properties.effective_heat_capacity`.

---

## 16.5 The energy ledger

### 16.5.1 Continuous-to-discrete construction

The ledger is assembled from *exactly the terms the RHS computes*. Integrating
the energy balance of §16.2.2 in time and discretizing each term by the same
quadrature the RHS uses gives

$$
E(t) = E(0)
+ \int_0^t \Big[\, q_\mathrm{cond,boundary}\,A
+ Q_\mathrm{rxn}^\mathrm{tot}
+ S_\mathrm{conv}^\mathrm{tot} \,\Big]\,d\tau,
$$

where each integrand is the *volume-integral of the corresponding RHS field*:

- $q_\mathrm{cond,boundary}$ — net **conducted** boundary heat flux [W·m⁻²],
  positive into the slab. This is the flux the Newton boundary solver drives
  into the first/last cell (§12.2); the surface balance has no gas-enthalpy
  term, so it equals the applied BC flux.
- $Q_\mathrm{rxn}^\mathrm{tot}=\sum_i Q_{\mathrm{rxn},i}V_i = \int_\Omega Q_\mathrm{rxn}\,dV$ [W],
  the total reaction heat ($Q_\mathrm{rxn}=-\sum_r h_r r_r$; §6.6).
- $S_\mathrm{conv}^\mathrm{tot}=\int_\Omega S_\mathrm{conv}\,dV$ [W], the total gas
  advective enthalpy sink; **typically negative** (heat carried out by gas; §7.4).

### 16.5.2 Per-substep accumulation

The energy tracker accumulates the time-integrals as Euler products over each
substep $\Delta t$:

$$
\begin{aligned}
\texttt{cum\_boundary\_conduction} &\mathrel{+}= q_\mathrm{cond,boundary}\,A\,\Delta t,\\
\texttt{cum\_reaction\_heat}       &\mathrel{+}= Q_\mathrm{rxn}^\mathrm{tot}\,\Delta t,\\
\texttt{cum\_advection\_sink}      &\mathrel{+}= S_\mathrm{conv}^\mathrm{tot}\,\Delta t .
\end{aligned}
$$

The boundary term carries the area factor $A$ and has flux units [W·m⁻²], so
$q\cdot A\cdot\Delta t$ has units of joules; the source terms are already
volume-integrals [W], so $\cdot\,\Delta t$ gives joules.

> **Implementation note.** `update_energy_tracker!`, `conservation_checks.jl`.
> The accumulators live in `EnergyConservationTracker`: `initial_energy`,
> `cumulative_boundary_conduction`, `cumulative_reaction_heat`,
> `cumulative_advection_sink`, `last_time`.

### 16.5.3 The residual

At a check point the *expected* energy is the initial energy plus the
accumulated sources,

$$
E^\mathrm{exp}(t) = E(0)
+ \texttt{cum\_boundary\_conduction}
+ \texttt{cum\_reaction\_heat}
+ \texttt{cum\_advection\_sink},
$$

and the **residual** is the gap between the measured and expected stored energy:

$$
\boxed{\,r(t) = E(t) - E^\mathrm{exp}(t)\,}
$$

with relative magnitude $r/|E(0)|$ (absolute if $|E(0)|\le\texttt{eps()}$) and
the verdict $r^\mathrm{rel}<\texttt{tol}$, default $\texttt{tol}=10^{-6}$.

> **Implementation note.** `check_energy_conservation`, `conservation_checks.jl`,
> returning `residual`, `relative`, `conserved`, `current_energy`.

### 16.5.4 Why it closes — and what it does *not* prove

This residual closes to integrator tolerance *by construction* whenever the
heat capacity and volume are constant over the trajectory and no radiation is
absorbed in depth, because both
$E(t)$ and its expectation are assembled from the same RHS fields under exactly
consistent quadrature; the only error is the integrator's local truncation error
in advancing $T$, which the implicit solver controls to `reltol`. A non-closing
residual therefore points to a **bookkeeping bug**: a flux computed but not
accumulated, a callback firing at the wrong time, a sign error in a source term,
or a missing factor of $A$ or $V$.

The crucial limitation, repeated from §16.1: the ledger **cannot detect a
double-count**. If the RHS were (incorrectly) to add a surface gas-enthalpy term
on top of $S_\mathrm{conv}$, the accumulator would faithfully record both, the
expectation would include both, and the residual would still close. Catching that
class of error is the job of the independent verification cases (§18). The ledger is also subject to the
composition-change and ALE volume-change omissions of §16.7, which prevent
*exact* closure whenever the mixture heat capacity or cell volume vary in time.

---

## 16.6 Boundary-flux semantics

The energy ledger accumulates the **conducted flux** $q_\mathrm{cond,boundary}$:
the heat actually conducted into the first/last interior cell, i.e.
$k_\mathrm{eff}(T_s-T_\mathrm{cell})/\Delta z_{1/2}$ evaluated at the converged
surface temperature $T_s$ from the Newton boundary solve (§12.2). This is what the
energy equation receives. The surface energy balance the Newton solver enforces
is (§12.2)

$$
q_\mathrm{BC,applied}(T_s) = q_\mathrm{cond}(T_s),
$$

with no gas-enthalpy term (§1.4.3, §7.4.3), so the conducted flux *is* the
applied flux: the externally
driven flux that the boundary condition imposes (for a `RadiativeFluxBC` the
absorbed incident radiation less the re-radiation and convective exchange; for a
pure `HeatFluxBC` the prescribed $q(t)$). The solver exposes these as
`extras.surface_heat_flux_net` ($q_\mathrm{BC,applied}$ — the applied face flux in
every radiation regime; under `BEER_LAMBERT` it excludes the incident radiation,
which enters volumetrically), `extras.surface_heat_flux_conduction`
($q_\mathrm{cond}$), plus `extras.volumetric_absorption` (the Beer–Lambert in-depth
deposit $\alpha_{\text{surf}}q_{\text{ext}}(1-e^{-\tau_{\text{tot}}})$; zero for other
models) and `extras.back_heat_flux` for the bottom boundary
(`src/Problem/pyrolysis_problem.jl`).

---

## 16.7 What the energy ledger does *not* yet accumulate

The residual of §16.5 closes exactly only under *constant heat capacity and
constant cell volume*. Pyrolysis violates both: the mixture composition $\xi$
evolves (changing $\rho c_p^\mathrm{eff}$ even at fixed $T$), and in ALE mode the
cell volumes $V_i$ change. Differentiating $E=\sum_i \rho c_{p,i}\,T_i\,V_i$ in
time gives three terms:

$$
\frac{dE}{dt}
= \sum_i \rho c_{p,i}\,V_i\,\frac{dT_i}{dt}
\;+\; \sum_i T_i\,V_i\,\frac{d(\rho c_{p,i})}{dt}
\;+\; \sum_i \rho c_{p,i}\,T_i\,\frac{dV_i}{dt} .
$$

Only the **first** term is what the source accumulators reconstruct (it is the
RHS times $\rho c_p V$). The temperature equation evolves $T$, *not* $\rho c_p\,T$,
so the ledger is **missing**:

1. **Composition-change proxy term** $\displaystyle\int T\,\frac{\partial(\rho c_p^\mathrm{eff})}{\partial t}\,dV\,dt$
   — the change of the tracked $\rho c_p\,T$ proxy density at constant temperature as
   reactions convert high-$c_p$ virgin solid into low-$c_p$ char and release gas,
   and as the gas concentration changes by transport and dilation. (This term
   belongs to the ledger's proxy $E$, not to a physical enthalpy density — the
   conserved quantity of the continuous formulation is $\sum_j\xi_j h_j(T)$,
   §3.3.2.) It is generally nonzero whenever reactions are active.
2. **ALE volume-change term** $\displaystyle\int \rho c_p^\mathrm{eff}\,T\,\frac{\partial V}{\partial t}\,dt$
   — the stored energy carried by a shrinking/expanding cell as its volume
   changes under mesh motion (§11). For the dilation part this cancels the
   corresponding piece of term 1 exactly (dilation conserves $\xi_jV$).

Consequently the residual as currently coded **cannot close to integrator
tolerance during active pyrolysis with time-varying composition or a moving
mesh**. Closing it would require per-substep, per-cell accumulation of the
*incremental* stored energy $\rho c_{p,i}(t)\,\Delta T_i\,V_i(t)$ (which absorbs
the composition and volume effects into the measured $\Delta(\rho c_p T)$), i.e.
wiring `record_substep!` into `compute_rhs!` so the tracker sees every cell's
state at every substep. That plumbing **does not exist in the package**; the
verification case V-CPL-06 carries its own per-step ledger with the hand-off and
dilation terms made explicit and closes the ALE energy budget to $5\times10^{-5}$
of the boundary input (§18.4).

This is the most important caveat for any user reading a runtime residual: a
small *but nonzero* residual during active reaction is **expected**, not a bug,
and reflects the unaccumulated composition/volume terms — not an error in the
flux/source assembly.

---

## 16.8 Callback infrastructure and runtime logging

### 16.8.1 DiagnosticsState and DiagnosticsRow

The running state of both trackers, plus the emitted log, is held in a
single mutable container:

```julia
mutable struct DiagnosticsState{NC}
    mass::MassConservationTracker{NC}
    energy::EnergyConservationTracker
    stride::Int
    step_counter::Int
    last_t::Float64
    log::Vector{DiagnosticsRow}
end
```

It is constructed with `DiagnosticsState(mesh, material; stride=10)`, which seeds
each tracker from the current mesh state — establishing the $t=0$ baseline
($E(0)$, $m_{\mathrm{initial},j}$). The default `stride` of 10 means one log row
is emitted per ten accepted substeps, trading output volume for temporal
resolution.

Each log row is a concretely typed `NamedTuple`,

$$
\texttt{DiagnosticsRow} = \big(t,\ E_\mathrm{sensible},\ r_\mathrm{energy},\ r_\mathrm{energy}^\mathrm{rel},\ \epsilon_{m}^{\max}\big),
$$

with all-`Float64` fields (`callback.jl`). Using a concrete `NamedTuple` (rather
than `Dict{Symbol,Any}`) keeps the log type-stable and means adding a new
diagnostic is a field addition, not a stringly-typed key.

### 16.8.2 record_substep! — the full cumulative path

`record_substep!` (`callback.jl`) is the entry point for the *full* cumulative
ledger. Per call it:

1. computes the substep size $\Delta t = t - \texttt{state.last\_t}$ and updates
   `last_t`;
2. **guards** against $\Delta t\le0$ (a callback fired at a duplicate time, e.g.
   the initial point) by returning early, before `step_counter` is incremented — a
   duplicate time neither pollutes the accumulators nor advances the stride;
3. increments `step_counter` and calls `update_mass_tracker!` and `update_energy_tracker!` with the substep's
   conducted boundary flux, the two volume-integrated source totals
   ($Q_\mathrm{rxn}$, $S_\mathrm{conv}$), the boundary mass-flux matrix, and the
   area;
4. **only on stride boundaries** (`step_counter % stride == 0`) calls
   `check_energy_conservation` and `check_mass_conservation` against the current
   mesh and pushes a fully populated `DiagnosticsRow` (with real residuals).

This path requires the caller to supply, every substep, the conducted boundary
flux, the two volume-integrated source totals, the boundary mass-flux matrix,
and the area — i.e. it must be threaded through the residual assembly. As noted
in §16.7, that threading is not yet wired, so this path is exercised by tests
rather than by `solve()`.

### 16.8.3 snapshot_diagnostics_row — the production path

What `solve()` *actually* uses is the **snapshot** path.
`snapshot_diagnostics_row(t, mesh, material)` (`callback.jl`) computes only the
*instantaneous* sensible energy $E(t)$ and fills the cumulative fields with
`NaN`:

$$
r_\mathrm{energy}=r_\mathrm{energy}^\mathrm{rel}=\epsilon_m^{\max}=\texttt{NaN}.
$$

The `NaN` sentinel is a deliberate signal that those ledger residuals are *not
available in snapshot mode* (they require per-substep accumulation). The snapshot
path is cheap — one reduction over the mesh — and gives the user the quantity
that *is* well-defined point-by-point: the energy content.

### 16.8.4 DiagnosticsCallback and solver integration

`DiagnosticsCallback(record!)` (src/Diagnostics/callback.jl) wraps a user
closure in a `SciMLBase.DiscreteCallback` whose condition is `true` (fire on
every accepted step) and whose `affect!` runs `record!(u, t)` on the
integrator's state and time, with `save_positions = (false, false)` so the
callback does not force extra solution saves. `record!` sits behind a
fixed-signature `FunctionWrapper`, so the callback's type does not depend on the
closure's (§15.5.3).

In `solve()`, enabling diagnostics (`diagnostics=true`, default `false`;
`diagnostics_stride` default `10`, validated $\ge1$) installs a snapshot-mode
closure (`solve`, src/Solver/solve.jl). The closure maintains its own
stride counter; on stride boundaries it (i) synchronizes the mesh from the
integrator's state with `update_mesh_from_state!(mesh, u, layout)` so
the trackers read the current cell states, then (ii) pushes a
`snapshot_diagnostics_row` to a local `Vector{DiagnosticsRow}`. After the solve,
the rows are appended to `solution.extras.diagnostics_log` (the `extras`
schema, `_extras` in src/Problem/pyrolysis_problem.jl, reserves
`diagnostics_log`). Reading `sol.extras.diagnostics_log` therefore yields the
energy time series, with `NaN` in the cumulative-residual columns under the
current production wiring.

---

## 16.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

The conservation-diagnostics design adds an explicit, runtime ledger audit that
is not a standard feature of any of the three reference codes.

**Energy accounting.** **FDS** performs an internal solid-phase energy budget and
exposes global energy-conservation accounting (the FDS Technical Reference
Guide's energy-budget appendix, Appendix P analog), reporting cumulative
boundary, radiative, and chemical contributions against stored enthalpy.
Pyrolysis.jl's energy ledger is the direct analog at the cell level: both
reconstruct stored energy from accumulated boundary and source fluxes and report
the closure residual.

**Mass accounting.** **Gpyro** advances control-volume conservation of
$\bar\rho\,\Delta z$ and $\bar\rho Y_i\,\Delta z$ with a fully-implicit
backward-Euler/Patankar scheme; conservation is enforced by the discretization
itself, with closure verified through the solid-density bookkeeping. Pyrolysis.jl
instead conserves species mass by construction (stoichiometry validated when each
reaction is built, §6.5.3) and *audits* it with the global mass ledger.

**Reference-point deformation vs. ALE.** ThermaKin handles condensed-phase
deformation through a single-reference-point density-integral source term, and
Gpyro through a variable-$\Delta z$ control volume; Pyrolysis.jl uses ALE mesh
motion (§10, §11). This is why the Pyrolysis.jl ledgers must carry a
time-dependent area $A(t)$ (trapezoidal average over a substep) and why the
unaccumulated ALE volume-change term (§16.7) is a closure caveat unique to the
ALE formulation that the reference codes, with their different deformation
treatments, do not incur in the same form.

---

## 16.10 Worked summary of the ledger equations

For quick reference, the complete set of accumulated and checked quantities:

| Ledger | Stored quantity | Expected | Residual | Tolerance |
|---|---|---|---|---|
| Mass (per comp. $j$) | $m_j=\sum_{i\in\text{act}}\xi_{j,i}V_i$ | $m_{\mathrm{init},j}-m_{\mathrm{efflux},j}$ | $m_j-m_{\mathrm{init},j}+m_{\mathrm{efflux},j}$ | $10^{-10}$ (rel.; abs. if $m_\mathrm{init}\approx0$) |
| Energy | $E=\sum \rho c_p^\mathrm{eff}T_iV_i$ | $E(0)+\!\int\!(q_\mathrm{cond}A+Q_\mathrm{rxn}+S_\mathrm{conv})$ | $E-E^\mathrm{exp}$ | $10^{-6}$ (rel.) |

Sign-convention reminders carried through every row: $J>0$ = inflow (so efflux
$=-J A\Delta t$); divergence positive = outflow; $Q_\mathrm{rxn}=-\sum_r h_r r_r$
with $h>0$ endothermic; $S_\mathrm{conv}$ typically negative.

---

## 16.11 Limitations and open issues

The following are the known gaps in the conservation diagnostics. They are
caveats on *interpretation*, not correctness bugs in the kernels.

1. **Composition-varying closure (primary limitation).** As derived in §16.7,
   the residual omits the composition-change proxy term
   $\int T\,\partial(\rho c_p^\mathrm{eff})/\partial t\,dV\,dt$ and the ALE
   volume-change term $\int \rho c_p^\mathrm{eff}\,T\,\partial V/\partial t\,dt$.
   A residual that closes to integrator tolerance during active pyrolysis
   requires per-substep, per-cell accumulation of $\rho c_{p,i}\,\Delta T_i\,V_i$,
   i.e. wiring `record_substep!` into `compute_rhs!`. This is **not implemented**;
   in production, `solve()` runs the snapshot path with `NaN` residual columns.

2. **Snapshot incompleteness.** Because `solve()` uses
   `snapshot_diagnostics_row`, the production `diagnostics_log` carries only
   $E_\mathrm{sensible}$; the residual and the mass-error column are `NaN`. The
   full cumulative ledger is available only by manually threading
   `record_substep!`.

3. **Boundary-flux decomposition.** The energy ledger treats the boundary flux
   as a single scalar per boundary. Decomposing it into separate radiation,
   convection, and conduction channels would let each mechanism be audited
   independently — the solver already exposes the *components* in `extras`
   (`surface_heat_flux_absorbed`, `_reradiation`, `_convection`, `_conduction`),
   but the tracker consumes only the aggregate.

4. **No radiation term.** The energy ledger accumulates the conducted boundary
   flux, $Q_\mathrm{rxn}$, and $S_\mathrm{conv}$ only; a volumetric Beer–Lambert
   deposit $Q_\mathrm{rad}$ is not among its accumulators, so under `BEER_LAMBERT`
   the residual also carries the absorbed in-depth radiation.

---

*Cross-references.* The governing equations these ledgers audit are stated in §3
(Governing Equations); the energy-source term $S_\mathrm{conv}$, the mixture heat
capacity, and the absence of a surface gas-enthalpy term in §7 (Heat Conduction,
Convection & Energy Assembly); the
heat-of-reaction sign and stoichiometric closure in §6 (Reaction Kinetics); the
surface energy balance and boundary fluxes in §12 (Boundary & Initial
Conditions); the time-dependent area and volume from §10 (Volume Change) and §11
(ALE Framework); the active-cell mask from §14 (Adaptive Mesh Refinement and
Surface Depletion); and
the callback/integration machinery from §15 (Temporal Integration). The
diagnostics are themselves a verification tool used in §18 (Verification &
Validation).


---

# 17. Sensitivity Analysis (Theory)

This chapter develops the mathematics of parametric sensitivity for the
Pyrolysis.jl solver: how the package computes the derivative of a simulation
output with respect to model parameters, why those derivatives are
well-defined, and how they support inverse parameter estimation, uncertainty
quantification, and gradient-based optimization. We treat both *local*
sensitivity (forward-mode algorithmic differentiation through the entire
`solve` chain, and reverse-mode adjoint differentiation of a scalar loss) and
*global*, variance-based sensitivity (Sobol indices computed through the
parameter-sweep harness). Throughout, the central enabling fact is that the
unified residual of §15 is written to be smooth and scalar-type-generic: the
state cache, the state-derived geometry, the boundary-temperature Newton solve,
and every property evaluation carry the scalar element type `Tu` end to end, so
a `ForwardDiff.Dual` (or a Zygote pullback) propagates cleanly from a parameter
all the way to a quantity of interest. The implementation lives in the optional
`PyrolysisSciMLSensitivityExt` extension behind two stub functions; the
mathematics below is what those wrappers actually differentiate.

Conventions used in this chapter (binding, from §2, Nomenclature):

- $\boldsymbol{\theta}\in\mathbb{R}^{p}$ is the vector of **parameters** with
  respect to which sensitivities are taken (e.g. $A_r$, $E_r$, $h_r$, a density
  $\rho_j$, a conductivity $k_j$, or a boundary-condition coefficient). $p$ is
  the parameter count.
- $\mathbf{u}(t)\in\mathbb{R}^{N_{\text{state}}}$ is the flat ODE state vector
  packed block-major as $[\,T \mid \xi_1\dots\xi_{N_C}\mid z\mid\chi\,]$ (§15.1).
- $\mathbf{y}=\mathcal{O}(\mathbf{u},\boldsymbol{\theta})\in\mathbb{R}^{m}$ is a
  vector **output** (a quantity of interest, "QoI") extracted from the
  solution; $\ell(\boldsymbol{\theta})\in\mathbb{R}$ is a **scalar loss**.
- The index $r$ is a reaction index, $j$ a component index, $i$ a cell index
  (overload resolution I1 of §2).
- $\sigma$ denotes the Stefan–Boltzmann constant where it appears in physics,
  not a standard deviation; variances in the global-sensitivity section are
  written $\mathbb{V}[\cdot]$ to avoid the clash.

---

## 17.1 Motivation: why sensitivities

A pyrolysis simulation maps a finite set of material and boundary parameters
$\boldsymbol{\theta}$ — Arrhenius pre-exponentials $A_r$ and activation energies
$E_r$, heats of reaction $h_r$, component densities $\rho_j$, conductivities
$k_j$, emissivities $\varepsilon_j$, incident-flux histories, convective
coefficients $h_{\text{conv}}$ — to time-resolved fields and to derived scalar
observables: surface temperature $T_s(t)$, mass-loss rate (MLR), peak heat
release rate, time to ignition, residual char fraction. The **sensitivity** of
an output to a parameter,

$$
S_{kq} \;=\; \frac{\partial y_k}{\partial \theta_q},
$$

is the workhorse quantity for three distinct tasks.

**(1) Inverse parameter estimation (calibration).** Cone-calorimeter or TGA
data fix observables $\mathbf{y}^{\text{data}}$; one seeks $\boldsymbol{\theta}$
minimizing a misfit
$\ell(\boldsymbol{\theta})=\tfrac12\lVert\mathbf{y}(\boldsymbol{\theta})-\mathbf{y}^{\text{data}}\rVert_W^2$.
Gradient-based optimizers (Gauss–Newton, L-BFGS, Levenberg–Marquardt) require
$\nabla_{\boldsymbol{\theta}}\ell$, which is exactly an adjoint or a
Jacobian-vector contraction of $S_{kq}$. The shipped inverse-analysis workflow
(`examples/inverse_analysis/inverse_white_pine.jl`, with `wp_objective.jl` and
`wp_forward.jl`) minimizes such a misfit for White-Pine kinetics and properties
against multi-flux cone data, derivative-free (`Fminbox(NelderMead())`), with a
`FiniteDiff` Hessian at the optimum for the standard errors.

**(2) Uncertainty quantification (UQ).** Given parameter uncertainties
(covariances from instrument calibration and per-curve regression), the
first-order propagated output covariance is
$\boldsymbol{\Sigma}_y \approx S\,\boldsymbol{\Sigma}_\theta\,S^{\mathsf T}$,
the *delta method*. Variance-based (Sobol) decomposition (§17.5) goes beyond
this linearization to apportion output variance to parameter groups.

**(3) Optimization and design.** Maximizing fire-retardant performance,
selecting char-forming additives, or designing a TGA heating program are all
gradient-amenable once $\nabla_{\boldsymbol{\theta}}\ell$ is available.

Pyrolysis.jl provides two complementary *local* derivative engines and a
*global* harness:

| Mode | Engine | Cost scaling | Best when |
|---|---|---|---|
| Forward (local) | `forward_sensitivity` / ForwardDiff | $\mathcal{O}(p)$ solves' worth | small $p$ ($\lesssim 20$), many outputs |
| Adjoint (local) | `adjoint_sensitivity` / Zygote + SciMLSensitivity | $\mathcal{O}(1)$ in $p$ per scalar loss | large $p$, scalar loss (inverse problems) |
| Global (variance) | parameter-sweep + surrogate (Sobol) | $\mathcal{O}(N_{\text{samples}})$ solves | apportioning variance, ranking groups |

Both local engines are *wrapper-based*: they differentiate the existing
`solve(problem)` end to end rather than rewriting the residual. This is the
single most important design decision in this chapter and is treated in §17.4.

---

## 17.2 Forward-mode local sensitivity (ForwardDiff)

### 17.2.1 The differentiated map

`forward_sensitivity` computes the full Jacobian of a vector output with respect
to the parameters by forward-mode algorithmic differentiation of the composite
map

$$
\boldsymbol{\theta}
\;\xmapsto{\;p_{\text{inject}}\;}\;
\text{problem}'
\;\xmapsto{\;\texttt{solve}\;}\;
\text{solution}
\;\xmapsto{\;\mathcal{O}\;}\;
\mathbf{y}.
$$

The extension wraps these three stages into a single closure and hands it to
`ForwardDiff.jacobian`:

```julia
wrapped = θ -> output_fn(solve(p_inject(base_problem, θ); progress=false, solve_kwargs...))
output   = wrapped(θ)
∂out_∂θ  = ForwardDiff.jacobian(wrapped, θ)
return (output, ∂out_∂θ)
```

Mathematically the returned matrix is

$$
\big[\partial\mathbf{y}/\partial\boldsymbol{\theta}\big]_{kq}
=\frac{\partial\, \mathcal{O}_k\big(\texttt{solve}(p_{\text{inject}}(\text{base},\boldsymbol{\theta}))\big)}{\partial \theta_q},
\qquad k=1\dots m,\;\; q=1\dots p ,
$$

evaluated at the supplied base point $\boldsymbol{\theta}$. The signature is
`forward_sensitivity(base_problem, θ, p_inject; output_fn=identity, solve_kwargs...)`;
`p_inject` is a positional closure and `output_fn` defaults to `identity`, which
returns the `PyrolysisSolution` itself; `ForwardDiff.jacobian` needs an array, so
pass an `output_fn` that returns a vector (e.g. `sol -> sol.surface_T`).

> **Implementation note.** Forward-mode sensitivity is implemented in
> `Pyrolysis.forward_sensitivity`, `ext/PyrolysisSciMLSensitivityExt.jl`, attached
> to the stub `forward_sensitivity` declared in `src/sensitivity_stubs.jl`. The
> stub has no methods until `using SciMLSensitivity, Zygote` loads the extension;
> calling it before raises a `MethodError`. ForwardDiff is the only differentiation engine used by the
> *forward* path; SciMLSensitivity is a load trigger, not an active component
> here.

### 17.2.2 Dual propagation end to end

Forward mode works because the package carries the scalar element type `Tu`
of the state through every cache and every kernel. `ForwardDiff.jacobian`
seeds $\boldsymbol{\theta}$ with $p$-component dual numbers
$\theta_q + \sum_{s} \delta_{qs}\,\epsilon_s$ (one partial $\epsilon_s$ per
parameter; in practice ForwardDiff chunks the $p$ partials). Each dual
satisfies the algebra
$f(\theta+\dot\theta\,\epsilon) = f(\theta) + f'(\theta)\,\dot\theta\,\epsilon$,
so propagating duals through the solve simultaneously evaluates the output and
its directional derivatives. For the chain to be exact, no stage may drop the
dual part or silently coerce to `Float64`. The package satisfies this:

- **State cache.** `StateCache{NC,Tu,Mode}` is parameterized on `Tu`; for AD the
  workspace is built with `Tu = ForwardDiff.Dual{Tag,Float64,N}` (§15,
  `build_ad_workspace`). The column cross-section area is stored as
  `A_section::Base.RefValue{Tu}` so that the $\chi$-driven area carries
  sensitivities. (See `src/Problem/state_cache.jl`.)
- **State-derived geometry.** In ALE mode, cell volumes
  $V_i=\lvert z_{i+1}-z_i\rvert\,A_{\text{section}}$, cell centers, and face
  positions are recomputed from the dual-typed node positions in
  `update_state_geometry!(cache, mesh, ::ALE)`; the mutable mesh is read-only
  (connectivity only). Geometry sensitivities therefore propagate (§11, The ALE
  Framework).
- **Property evaluation.** Property functions are callable on `Real`, hence on
  duals, and carry no derivative methods (§4.5): ForwardDiff differentiates each
  law's own arithmetic, so even `TableProperty` (slope of the active segment) and
  `CallableProperty` are differentiated exactly with respect to $T$ and to any
  parameter inside them.
- **Boundary-temperature Newton solve.** The surface-temperature solve
  $f(T_s)=0$ (§12.2) detects dual inputs and returns a dual $T_s$; the implicit
  function theorem is realized automatically because Newton's iterate
  $T_s\leftarrow T_s - f/f'$ is itself differentiable arithmetic once converged.
  (See `solve_boundary_temperature`, `src/Discretization/boundary_temperature.jl`.)
- **Kinetics.** `reaction_rate` is written with
  `one(promote_type(eltype(ξ), typeof(T)))` and `_fast_pow` so that the
  concentration product, tanh gates, and depletion limiter all stay in the
  promoted dual type (`src/Physics/kinetics.jl`).

The time integrator (`KenCarp4` by default, with the structured Jacobian of §15
as its Newton Jacobian and the default sparse factorization for the linear solve) does **not** need to be
dual-aware: ForwardDiff differentiates the *output of the converged solve*, not
the integrator's internal Jacobian. OrdinaryDiffEq's `solve` is itself a
differentiable program, and SciMLSensitivity supplies the dual-compatible
solve path so that `ForwardDiff.Dual` parameters flow through the time loop.

### 17.2.3 What the derivative actually is

For an output $y=\mathcal{O}(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta})$
evaluated at the final state, the total derivative obeys the chain rule

$$
\frac{\mathrm{d}y}{\mathrm{d}\theta_q}
=\frac{\partial\mathcal{O}}{\partial\boldsymbol{\theta}}\bigg|_q
+\frac{\partial\mathcal{O}}{\partial\mathbf{u}}\cdot
\frac{\partial\mathbf{u}(t_{\text{end}})}{\partial\theta_q}.
$$

The second factor is the *state sensitivity*
$\mathbf{s}_q(t)\equiv\partial\mathbf{u}(t)/\partial\theta_q$, which (formally)
satisfies the **forward sensitivity ODE** obtained by differentiating the
residual $\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$:

$$
\dot{\mathbf{s}}_q
=\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\,\mathbf{s}_q
+\frac{\partial\mathbf{f}}{\partial\theta_q},
\qquad
\mathbf{s}_q(0)=\frac{\partial\mathbf{u}_0}{\partial\theta_q}.
\tag{17.1}
$$

Here $\partial\mathbf{f}/\partial\mathbf{u}$ is precisely the structured
Jacobian of §15. Pyrolysis.jl does **not** integrate (17.1) explicitly; the
forward-mode wrapper instead lets ForwardDiff carry the dual partials through
the existing time loop, which is algebraically equivalent to solving (17.1) by
the *same* discretization the integrator uses for the primal — i.e. it inherits
the integrator's order and step-size control automatically, including the
contribution of $\partial\mathbf{u}_0/\partial\theta_q$ when a parameter
appears in the initial condition (e.g. an initial density). Equation (17.1) is
the conceptual object; the dual-propagation implementation is the realized one.

### 17.2.4 Cost and when to use forward mode

Forward mode evaluates $p$ directional derivatives in lockstep with the primal.
Its cost is roughly that of $1+p$ primal solves (ForwardDiff chunks reduce the
constant but not the linear-in-$p$ scaling). It is therefore the right choice
when:

- the parameter count is small ($p\lesssim 20$ — a handful of reactions and
  property coefficients), and
- many outputs are wanted at once (the *full* $m\times p$ Jacobian comes for
  free; reverse mode would need $m$ passes for $m$ scalar outputs).

For a typical single-reaction calibration with $\boldsymbol{\theta}=(A,E,h)$ and
outputs $\{T_s(t_k),\text{MLR}(t_k)\}$, forward mode delivers the entire
sensitivity matrix in three dual-augmented solves' worth of work — ideal for
Gauss–Newton, whose normal equations need $S$ itself, not just
$S^{\mathsf T}\mathbf{r}$.

---

## 17.3 Adjoint-mode local sensitivity (Zygote + SciMLSensitivity)

### 17.3.1 The differentiated map and the gradient

When the quantity of interest is a *scalar* loss $\ell(\boldsymbol{\theta})$
(the inverse-problem case), reverse mode computes the entire gradient
$\nabla_{\boldsymbol{\theta}}\ell$ at a cost essentially independent of $p$.
`adjoint_sensitivity` wraps the loss-producing solve and differentiates it with
Zygote:

```julia
wrapped = θ -> loss_fn(solve(p_inject(base_problem, θ);
                             progress=false, sensealg=sensealg, solve_kwargs...))
loss = wrapped(θ)
grad = first(Zygote.gradient(wrapped, θ))
return (loss, grad)
```

The returned gradient is
$\big[\nabla_{\boldsymbol{\theta}}\ell\big]_q=\partial\ell/\partial\theta_q$.
The signature is
`adjoint_sensitivity(base_problem, θ, p_inject, loss_fn; sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), solve_kwargs...)`.

> **Implementation note.** Adjoint sensitivity is
> `Pyrolysis.adjoint_sensitivity`, `ext/PyrolysisSciMLSensitivityExt.jl`. The `sensealg` is **forwarded to `solve` as a first-class
> keyword** — `solve` accepts `sensealg` precisely so SciMLSensitivity can
> intercept the ODE solve and augment it with adjoint dynamics without tripping
> `solve`'s unknown-keyword guard (`src/Solver/solve.jl`). Zygote provides the
> reverse-mode pullback over the *non-ODE* parts of the chain (the
> `p_inject`/`loss_fn` glue); SciMLSensitivity provides the differentiable ODE
> solve in between.

### 17.3.2 The continuous adjoint and the augmented ODE

The mathematics SciMLSensitivity realizes is the *continuous adjoint*. Let the
loss be an integral (and/or endpoint) functional of the trajectory,

$$
\ell=\int_{0}^{t_{\text{end}}} g\big(\mathbf{u}(t),\boldsymbol{\theta},t\big)\,\mathrm{d}t
\;+\;\Phi\big(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta}\big),
$$

subject to $\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$.
Introducing the adjoint (costate) $\boldsymbol{\lambda}(t)\in\mathbb{R}^{N_{\text{state}}}$,
stationarity of the Lagrangian gives the **backward adjoint ODE**

$$
\dot{\boldsymbol{\lambda}}
=-\left(\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\right)^{\!\mathsf T}\boldsymbol{\lambda}
-\left(\frac{\partial g}{\partial\mathbf{u}}\right)^{\!\mathsf T},
\qquad
\boldsymbol{\lambda}(t_{\text{end}})=\left(\frac{\partial\Phi}{\partial\mathbf{u}}\right)^{\!\mathsf T}\bigg|_{t_{\text{end}}},
\tag{17.2}
$$

integrated from $t_{\text{end}}$ back to $0$, after which the gradient is the
quadrature

$$
\frac{\mathrm{d}\ell}{\mathrm{d}\theta_q}
=\frac{\partial\Phi}{\partial\theta_q}
+\int_{0}^{t_{\text{end}}}
\left(
\boldsymbol{\lambda}^{\mathsf T}\frac{\partial\mathbf{f}}{\partial\theta_q}
+\frac{\partial g}{\partial\theta_q}
\right)\mathrm{d}t .
\tag{17.3}
$$

The key structural facts:

- (17.2) is a single ODE of the same dimension as the state, **independent of
  $p$**. One backward solve yields the full gradient via the $p$ cheap
  quadratures (17.3). This is why reverse mode dominates when $p$ is large.
- The operator $(\partial\mathbf{f}/\partial\mathbf{u})^{\mathsf T}$ is the
  transpose of the structured Jacobian of §15. SciMLSensitivity never forms it
  densely; with `autojacvec=ReverseDiffVJP()` it evaluates *vector–Jacobian
  products* $\boldsymbol{\lambda}^{\mathsf T}(\partial\mathbf{f}/\partial\mathbf{u})$
  by reverse-mode AD of the residual.
- The default `InterpolatingAdjoint` stores (and interpolates) the forward
  trajectory so the backward pass can re-evaluate $\mathbf{f}$ and its VJPs at
  the correct states; this trades memory for the avoidance of a backward
  re-solve of the primal.

### 17.3.3 sensealg selection

The `sensealg` keyword chooses the adjoint algorithm. The default is
`InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`. The practical menu (all
from SciMLSensitivity) and their trade-offs:

| `sensealg` | Memory | Robustness | Notes |
|---|---|---|---|
| `InterpolatingAdjoint` (default) | high (stores forward sol) | high | recompute-free backward pass |
| `QuadratureAdjoint` | high | high | separate quadrature for (17.3); fast for few params |
| `BacksolveAdjoint` | low | lower (can drift for stiff Arrhenius) | re-solves primal backward; risky here |
| `GaussAdjoint` | moderate | high | checkpointed quadrature |
| `ForwardDiffSensitivity` | — | exact | reverts to forward mode; small $p$ only |

For pyrolysis the reactions are stiff and the trajectory can vary by orders of
magnitude across a decomposition front, so `BacksolveAdjoint` is generally
discouraged (its backward re-solve can leave the stiff manifold); the
interpolating or Gauss variants are preferred. The `autojacvec` choice
(`ReverseDiffVJP`, `ZygoteVJP`, `EnzymeVJP`, or `TrackerVJP`) sets how the VJP
$\boldsymbol{\lambda}^{\mathsf T}\partial\mathbf{f}/\partial\mathbf{u}$ is taken;
`ReverseDiffVJP` is robust for the scalar-heavy residual kernels here.

### 17.3.4 Cost and when to use adjoint mode

Reverse mode costs $\mathcal{O}(1)$ scalar-loss gradients independent of $p$
(one forward solve + one backward adjoint solve), at the price of higher memory
and the need for a fully reverse-differentiable solve path. Use it when:

- the loss is a single scalar (sum-of-squares misfit, integrated MLR error),
  and
- the parameter count is moderate-to-large (many property coefficients across
  several reactions, distributed-activation surrogates, spatially-varying
  initial fields).

The dividing line between §17.2 and §17.3 is essentially $p$ versus $m$: forward
mode for "few parameters, many outputs," adjoint mode for "many parameters, one
loss." For the in-between regime the wrapper lets a user simply swap
`forward_sensitivity` ↔ `adjoint_sensitivity` with the *same* `p_inject`.

---

## 17.4 The `p_inject` wrapper pattern and parameter handles

### 17.4.1 Why wrap rather than rewrite

Both engines avoid touching the residual signature. The residual
$\texttt{residual!}(du,u,\text{def},t,ws)$ takes an *immutable* `ProblemDef` and
a *mutable* `Workspace` (§15.2); it has no explicit parameter argument. Rather
than thread a parameter vector $\boldsymbol{\theta}$ through `def` and every
operator, sensitivity is obtained by differentiating the *construction* of the
problem followed by the solve. The user supplies a closure

$$
p_{\text{inject}}:\ (\text{base\_problem},\boldsymbol{\theta})\ \longmapsto\ \text{problem}',
$$

that materializes a *fresh* `PyrolysisProblem` with $\boldsymbol{\theta}$ baked
into the relevant fields, and the engine differentiates
$\boldsymbol{\theta}\mapsto\texttt{solve}(p_{\text{inject}}(\dots))$.

This has three consequences worth stating precisely:

1. **No invasive refactor.** The residual, the operators, the Jacobian backends,
   and the state layout are untouched. The cost is that AD must re-run *problem
   construction* (cheap: building a `Material` and `BoundaryConditionSet`) on
   every evaluation; the solve dominates anyway.
2. **Any constructible parameter is differentiable.** Because $p_{\text{inject}}$
   builds the problem, *anything* that flows into a `Material`,
   `BoundaryConditionSet`, mesh, or initial condition is a valid handle — there
   is no fixed parameter registry.
3. **The base problem must not be mutated.** `p_inject` returns a *new*
   problem; mutating `base_problem` in place would break the purity that lets
   sweeps and AD share immutable state (the same contract the parameter-sweep
   harness enforces, §17.5).

### 17.4.2 Parameter handles

Common $\boldsymbol{\theta}$ entries and how `p_inject` injects them:

- **Kinetic parameters** $A_r$, $E_r$. Rebuild the reaction(s) with new
  pre-exponential/activation energy and assemble a new `Material`. The rate
  $r=A\,e^{-E/(R_gT)}\prod_j\xi_j^{n_{r,j}}$ is analytic in both $A$ (linear) and
  $E$ (through the exponential), so $\partial r/\partial A = r/A$ and
  $\partial r/\partial E = -r/(R_gT)$ are exact under AD. Because $A$ and $E$ are
  strongly correlated along the Arrhenius "compensation" direction, the
  sensitivity matrix is typically ill-conditioned in $(A,E)$ — see §17.7.
- **Heat of reaction** $h_r$ (or a `LinearProperty` $a+bT$, or a
  `StateDependentProperty` for moisture sorption). It enters only the energy
  source $Q_{\text{rxn}}=-\sum_r h_r r_r$ linearly, so
  $\partial Q_{\text{rxn}}/\partial h_r = -r_r$ exactly. With the sign
  convention $h>0$ endothermic, increasing $h_r$ *cools* the cell, lowering
  $T_s$ — the sensitivity $\partial T_s/\partial h_r<0$.
- **Densities / conductivities / emissivities** $\rho_j,k_j,\varepsilon_j$.
  These feed the mixture rules of §5 and the boundary energy balance; their
  derivatives propagate through the effective-property closures (ForwardDiff
  through every mixing rule, Bruggeman's Newton solve included; §5.15). A
  bulk-density change also shifts the initial dry-solid mass $m_{\text{dry},i}$
  (and hence $\chi$, §10), so density is a particularly broad-reaching handle.
- **Boundary-condition parameters** — incident flux level, convective
  coefficient $h_{\text{conv}}$, ambient $T_\infty$, view factor $F$,
  emissivity in $q_{\text{rad}}=F\varepsilon\sigma(T_\infty^4-T_s^4)$. These
  enter the surface energy balance solved for $T_s$; their sensitivities pass
  through the implicit Newton solve as in §17.2.2.

### 17.4.3 Worked structure of a forward-sensitivity call

The canonical usage that the extension differentiates is:

```julia
using Pyrolysis, SciMLSensitivity, Zygote   # loads the extension

θ0 = [A0, E0, h0]                            # base parameter point

function p_inject(prob, θ)
    rxn  = Reaction(:decomp, 1 => 2; A = θ[1], E = θ[2], h = θ[3])
    m    = prob.material
    mat  = Material(; name = m.name, components = m.components, reactions = (rxn,),
                    mixing = m.mixing, lateral_shrinkage_law = m.lateral_shrinkage_law,
                    depletion_limiter = m.depletion_limiter, gas_viscosity = m.gas_viscosity)
    return PyrolysisProblem(; mesh = deepcopy(prob.mesh), material = mat,
                            bc_set = prob.bc_set, T_initial = prob.T_initial,
                            ξ_initial = prob.ξ_initial, tspan = prob.tspan,
                            dt_initial = prob.dt_initial, output_times = prob.output_times)
end

output_fn(sol) = sol.surface_T               # vector QoI over saved times

(y, S) = forward_sensitivity(base_problem, θ0, p_inject; output_fn = output_fn)
# y :: surface-T trace; S[k,q] = ∂(T_s at t_k)/∂θ_q
```

`S` is the $m\times 3$ sensitivity matrix used directly in the Gauss–Newton
normal equations $S^{\mathsf T}WS\,\Delta\boldsymbol{\theta}
= -S^{\mathsf T}W\,\mathbf{r}$ for the residual
$\mathbf{r}=\mathbf{y}-\mathbf{y}^{\text{data}}$.

---

## 17.5 Global / variance-based sensitivity (Sobol)

### 17.5.1 What global sensitivity adds

Local sensitivities $S_{kq}$ are first-order Taylor coefficients at one point.
When parameter uncertainties are large or the response is nonlinear, one wants
*variance-based* (global) indices that hold over the whole input distribution.
Treat the parameters as random $\boldsymbol{\Theta}$ with a joint distribution
and a scalar QoI $Y=\mathcal{Q}(\boldsymbol{\Theta})$. The Sobol/ANOVA
decomposition apportions the output variance $\mathbb{V}[Y]$ to subsets of
inputs. For an input (or group) $\theta_q$ the **first-order** and **total**
indices are

$$
S_q^{\,1}=\frac{\mathbb{V}_{\theta_q}\!\big(\mathbb{E}_{\sim q}[Y\mid\theta_q]\big)}{\mathbb{V}[Y]},
\qquad
S_q^{\,T}=1-\frac{\mathbb{V}_{\sim q}\!\big(\mathbb{E}_{\theta_q}[Y\mid\boldsymbol{\theta}_{\sim q}]\big)}{\mathbb{V}[Y]},
\tag{17.4}
$$

where $\sim q$ denotes "all parameters except $q$." $S_q^1$ measures the variance
explained by $\theta_q$ acting alone; $S_q^T$ adds all interactions involving
$\theta_q$. They bracket the importance: $S_q^1\le S_q^T$, with equality iff
$\theta_q$ has no interactions.

### 17.5.2 Design, sampling, and the surrogate

Pyrolysis.jl does not bundle a Sobol estimator into the core solver; instead the
global workflow is built on the **parameter-sweep harness** (§17.5.3) plus a
*surrogate* fitted to the swept solutions. The shipped reference implementation
is `examples/sensitivity_analysis/sensitivity_sobol.jl`, which for the
Douglas-fir cone problem:

1. Draws correlated parameter samples (16 uncertain parameters in 7
   *independent* groups — calibration biases, low/high-T virgin $c_p$, char
   $c_p$, heat of sorption, heats of reaction). Correlations are *within-group
   only*, so the standard Sobol formulas (17.4) apply at the group level without
   correlated-input extensions.
2. Transforms the correlated draws to independent standard normals (Cholesky
   inverse for Gaussian groups; Rosenblatt sequential conditioning for
   Student-$t$ groups).
3. Solves the pyrolysis problem at each sample (the expensive step — one full
   cone solve per sample) and reduces each solution to scalar QoIs.
4. Fits an adaptive-degree (2–4) **Polynomial Chaos Expansion (PCE)** surrogate
   to each QoI, selecting degree by leave-one-out $Q^2$ and using a sparse LASSO
   fit when the oversampling ratio is low.
5. Computes first-order and total *group* Sobol indices **analytically from the
   PCE coefficients** — because the PCE basis is orthonormal under the input
   measure, each Sobol index is a normalized sum of squared coefficients, so no
   second Monte-Carlo loop is needed.
6. Quantifies index uncertainty with a BCa bootstrap and assesses convergence on
   nested subsamples.

The PCE surrogate is the practical enabler: it converts an
$\mathcal{O}(N_{\text{samples}})$ *solve* budget into closed-form indices,
avoiding the $\mathcal{O}((p{+}2)N)$ solves of a raw Saltelli/Sobol pick-freeze
estimator. The degree-selection rationale is documented in
`examples/sensitivity_analysis/pce_degree_selection.md`; `compare_degrees.jl`
sweeps the degree explicitly.

### 17.5.3 Relation to the parameter-sweep harness

The variance design is executed through `parameter_sweep(problem, grid; modify, extract)`
(§ User Guide UG-12). Its contract dovetails with the sampling above:

$$
\texttt{modify}(\text{base\_problem},\,\mathbf{p}_s)\ \longmapsto\ \text{problem}_s,
\qquad
\texttt{extract}(\text{solution}_s)\ \longmapsto\ \mathbf{y}_s,
$$

for each sample point $\mathbf{p}_s$ in the design `grid`. `modify` plays the
same role as `p_inject` in §17.4 — it must build a *fresh* problem and must not
mutate the base. The harness dispatches the (embarrassingly parallel) solves
over `Distributed.pmap` when workers are available, giving each worker its own
mutable `Workspace` while sharing the immutable `ProblemDef`:

```julia
work = point -> extract(solve(modify(problem, point); progress=false, solve_kwargs...))
results = (nworkers() > 1) ? pmap(work, pts) : map(work, pts)
```

(`src/Solver/parameter_sweep.jl`.) The `extract` closure should return small
NamedTuples of scalar QoIs (peak $T_s$, integrated MLR, time-to-peak) to
minimize worker→head transfer; those scalars are exactly the $Y$ samples the PCE
surrogate consumes. Thus the *same* solve and the *same* immutable-problem
discipline underpin local sensitivity (§17.2–17.4), parameter sweeps, and global
Sobol analysis — one harness, three derivative products.

### 17.5.4 Local-vs-global consistency

The two views are linked by the *delta method*: near a base point, the
first-order propagated variance is
$\mathbb{V}[Y]\approx\sum_{q,q'}S_q S_{q'}\,\mathrm{Cov}(\theta_q,\theta_{q'})$
with $S_q=\partial Y/\partial\theta_q$ from §17.2/§17.3. When the response is
linear in the inputs over their uncertainty range, the Sobol $S_q^1$ equals the
normalized squared local sensitivity
$S_q^2\,\mathbb{V}[\theta_q]/\mathbb{V}[Y]$ and $S_q^T=S_q^1$ (no interactions).
Discrepancies between the local and global rankings are therefore *diagnostics
of nonlinearity/interaction*, not contradictions — a useful cross-check when
both are available.

---

## 17.6 Interaction with AD-clean design choices

Sensitivities are only meaningful if the map $\boldsymbol{\theta}\mapsto\mathbf{y}$
is (at least) continuously differentiable. Several modeling choices in
Pyrolysis.jl exist *specifically* to make this true, and they are precisely the
smooth replacements for the classical non-smooth pyrolysis closures. Each is
verified against the source below.

**(1) Smooth temperature gates.** The hard Heaviside window
$H(T-T_{\min})H(T_{\max}-T)$ would make $\partial r/\partial T$ a delta function
at the gate edges and destroy both the integrator's Newton Jacobian and any AD
derivative. The code replaces it with the product of tanh ramps

$$
r \;\mathrel{*}=\; \tfrac12\big(1+\tanh((T-T_{\min})/w_{\text{gate}})\big)\,
                   \tfrac12\big(1+\tanh((T_{\max}-T)/w_{\text{gate}})\big),
$$

applied only when the bounds are finite/positive (`src/Physics/kinetics.jl`), with
the per-reaction gate width $w_{\text{gate}}$ (default 1 K, §6.2.2).
This is $C^\infty$ in $T$, so all temperature sensitivities of the rate are
finite. The cost is a leakage of reactivity of about $10\,w_{\text{gate}}$ outside
the nominal window (the ramp is not unit-step), discussed in §6.2; for sensitivity
purposes the benefit — globally finite $\partial r/\partial T$ — outweighs it.

**(2) Depletion limiter.** Near zero reactant the bare Arrhenius rate is stiff
and its parameter sensitivities blow up. Each reactant multiplies the rate by
the smooth factor

$$
\tanh(\xi/\xi_{\text{thr}}),
$$

with default threshold $\xi_{\text{thr}}=1.0\ \text{kg·m}^{-3}$
(`src/Physics/kinetics.jl`; reactant clamped via $\max(\xi,0)$ first).
The rate vanishes at depletion with a bounded $\partial r/\partial\xi$ for
every reaction order $n \ge 0$, keeping the reaction-rate Jacobian — and hence
the parameter sensitivities of MLR near burnout — bounded and well-defined.
The factor is $C^\infty$ in $\xi > 0$ with no threshold branch; the only
non-smooth point left is the $\max(\xi,0)$ clip exactly at $\xi=0$ (§17.7).

**(3) Smooth flux limiters and lateral-shrinkage derivative.** The TVD limiters
(Van Leer is the ALE default) are continuous, and the $\chi$-driven radial
strain rate uses an *automatic-differentiation* derivative of the
lateral-shrinkage law,
$\theta_A=\dfrac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\dfrac{\mathrm{d}\bar\chi}{\mathrm{d}t}$
with `dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)` and a guard
`A_ratio > eps()` (`src/Residual/residual.jl`). Because $\mathrm{law}'$
is taken by AD, *any differentiable* law contributes correct geometry
sensitivities — but a non-smooth law breaks this (see §17.7).

**(4) Harmonic-mean face averaging and tanh-based smoothing elsewhere.** The
distance-weighted harmonic mean used for face $k_f,\lambda_f,\kappa_f$ (§5.12)
and the smooth (rather than `min`/`max`) constructions in the gates keep the
discrete operators differentiable. The few genuinely non-smooth operations that
remain — `clamp` on porosity/temperature, `max` in the depletion clamp, and the
state rewrite at park events — are localized and are exactly the pitfalls
catalogued next.

> **Implementation note.** The smoothness-critical kernels are
> `reaction_rate` (`src/Physics/kinetics.jl`), the limiter library
> (`src/Discretization/flux_limiters.jl`), and `_compute_theta_A`
> (`src/Residual/residual.jl`). Each was written with `tanh`/AD
> derivatives *in place of* the textbook step/cutoff/finite-difference form so
> that the chain $\boldsymbol{\theta}\mapsto\mathbf{y}$ is differentiable.

---

## 17.7 Pitfalls

Even with an AD-clean residual, sensitivity computations can return wrong or
meaningless numbers if certain solver features are active or the parameter point
sits on a non-smooth manifold. The pitfalls below are specific to Pyrolysis.jl.

### 17.7.1 Warm-start artifacts

`solve(...; warm_start=prev_solution)` reuses the prior solve's *final* state as
the new initial condition (`_apply_warm_start!`, `src/Solver/solve.jl`). Inside
a sensitivity loop this introduces a **hidden parameter dependence that AD does
not see**: the warm-start state was produced at a *different* $\boldsymbol{\theta}$,
so $\partial\mathbf{u}_0/\partial\theta_q$ is silently set to zero even though
the true initial condition (had the perturbed problem been integrated from the
real start) would depend on $\boldsymbol{\theta}$. Differentiating a
warm-started solve therefore yields a derivative of the *wrong map*. **Rule:** do not warm-start
solves that are inside `forward_sensitivity`/`adjoint_sensitivity` or a
sensitivity sweep; integrate each sample from its true initial condition.

### 17.7.2 Non-smooth lateral-shrinkage laws

The radial strain rate uses $\mathrm{law}'(\bar\chi)$ via AD (§17.6, item 3). Identity
(`NoChi`) and the smooth linear laws $1\mp k\bar\chi$ are differentiable, but a
piecewise or table-defined law with kinks makes $\mathrm{law}'$ discontinuous at
the breakpoints. A sensitivity evaluated *at* a breakpoint is then ill-defined
(left/right derivatives differ), and AD returns whichever branch the
floating-point comparison happens to take. The guard `A_ratio > eps()` further
*zeroes* $\theta_A$ — and thus its sensitivity — when the cross-section nearly
vanishes, a hard switch that is correct physically but non-differentiable.
**Rule:** use $C^1$ shrinkage laws for sensitivity studies, and avoid evaluating
at $\bar\chi$ values where the law has a corner or where $\mathrm{law}(\bar\chi)\to 0$.

### 17.7.3 Park-event discontinuities

The depletion path (§14.3) parks a depleted cell: it merges the cell into a
neighbour and moves the freed slot above the live surface, rewriting the state
vector in place. The state keeps its size and the structured Jacobian keeps
working, but the rewrite is a **discrete event**: a perturbation of
$\boldsymbol{\theta}$ that shifts a park *time* produces a jump in the
trajectory, so $\mathbf{u}(t)$ is only piecewise smooth in the parameters. Forward
sensitivities are valid on horizons without park events; adjoint sensitivity,
which needs a differentiable solve path, cannot be run through a park event.

### 17.7.4 Arrhenius $(A,E)$ ill-conditioning

The pre-exponential and activation energy are nearly redundant along the
compensation direction $\delta\ln A \approx \delta E/(R_g T_{\text{peak}})$: a
joint perturbation that holds the rate fixed at the dominant reaction
temperature leaves the output almost unchanged. The sensitivity matrix $S$ is
consequently near-rank-deficient in the $(A,E)$ block, so Gauss–Newton normal
equations $S^{\mathsf T}WS$ are ill-conditioned. This is a property of the
*physics*, not a bug; mitigate by reparameterizing (e.g. fit $E$ and the rate at
a reference temperature instead of $A$ directly), by regularization
(Tikhonov/priors), or by using observables at *multiple* heating rates to break
the degeneracy.

### 17.7.5 Integrator-tolerance noise in finite-difference checks

When validating AD sensitivities against finite differences, remember the primal
solve is only accurate to the integrator tolerances (`abstol=1e-8`,
`reltol=1e-6` by default). A finite-difference step $\delta\theta$ smaller than
the solve noise yields garbage; a step too large incurs truncation error. AD
(forward or adjoint) avoids this entirely — the dual/adjoint derivative is
*exact* for the discretized solve at the chosen tolerances — which is precisely
why the package differentiates the solve rather than finite-differencing it. If a
finite-difference cross-check is needed, tighten `abstol`/`reltol` first.

### 17.7.6 Non-differentiable guards in the kinetics and clamps

The reactant clamp $\max(\xi,0)$, the temperature clamp to $[200,3000]\ \text{K}$
in the boundary Newton solve, and the porosity `clamp(\cdot,0,1)` are all
non-smooth at their boundaries. In normal operation states stay in the interior
and these are inactive, but a sensitivity evaluated exactly *at* a clamp boundary
(e.g. a fully depleted reactant, $\xi=0$) inherits a one-sided or zero
derivative. This is rarely an issue for QoIs integrated over time but can surface
in pointwise late-time MLR sensitivities near burnout.

---

## 17.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

Pyrolysis.jl's sensitivity capability is structurally different from its
condensed-phase peers, which expose parameter sensitivity only
*externally* (through driver scripts and finite differences) rather than through
differentiable internals.

**ThermaKin2Ds.** ThermaKin couples a forward pyrolysis solver to inverse
analysis via external optimization (e.g. Stochastic Hill-Climbing / Shuffled
Complex Evolution over the kinetic and property parameters), evaluating the
forward model as a black box and forming sensitivities by repeated forward runs
(finite differences). Mapping notation: ThermaKin's stoichiometric coefficients
$\theta_i^j$ are our yields $\nu_{r,j}$, and its heats of reaction carry the
*opposite* sign ($h>0$ exothermic) which our $Q_{\text{rxn}}=-\sum_r h_r r_r$
converts to (overload H1, §2). Pyrolysis.jl's advantage is that the *same*
parameter handles ($A_r,E_r,h_r,\rho_j,k_j$) are differentiated *analytically
through the solve* by AD, so gradients are exact at the integrator tolerance and
cost $\mathcal{O}(p)$ (forward) or $\mathcal{O}(1)$ (adjoint) rather than $p{+}1$
black-box forward runs per gradient.

**Gpyro.** Gpyro's calibration (Lautenberger & Fernandez-Pello) is built on a
*genetic algorithm* (GA) that minimizes the mismatch between simulated and
measured MLR/temperature over a population of parameter sets — derivative-free by
design. Gpyro maps the pre-exponential $Z$ → our $A_r$ and the conversion-based
nth-order kinetics to our mass-concentration form (§6.9). A GA explores a global
landscape without gradients (robust to the $(A,E)$ degeneracy and to non-smooth
responses) but needs many forward solves and gives no local sensitivity matrix.
Pyrolysis.jl is complementary: AD gradients drive fast *local* refinement
(Gauss–Newton/L-BFGS) once a GA or sweep has found a basin, and the
variance-based Sobol harness (§17.5) supplies the global parameter ranking that
a GA does not.

**FDS (solid phase).** The FDS solid-phase pyrolysis model recovers kinetic
parameters from a reference TGA temperature and peak rate (the
$A,E$-from-$T_{\text{peak}}$ recipe, §6.9) rather than offering a built-in
sensitivity/adjoint capability; FDS's broader uncertainty treatment is the
verification-and-validation accounting of its guides, not parametric
differentiation of the solid solver. Pyrolysis.jl differs in providing
first-class, AD-exact local sensitivities and a surrogate-based global
decomposition as part of the solver API.

In short: where the reference codes treat the forward model as a black box and
obtain sensitivities by repeated evaluation (finite difference) or avoid them
entirely (GA), Pyrolysis.jl exposes the forward map as a differentiable program
and offers both exact local gradients (forward and adjoint) and variance-based
global indices, all sharing one immutable-problem / mutable-workspace harness.

---

## 17.9 Limitations and open issues

1. **Adjoint sensitivity across park events (open).** A park is a discrete
   rewrite of the state, so the trajectory is only piecewise smooth across
   events and adjoint sensitivity cannot run through a depleting
   (`handle_depletion=true`) horizon that contains one. Reconciling the event
   with a differentiable solve path is an open design problem.

2. **Warm-start in sensitivity loops (open).** Whether — and how — reusing a
   prior solution's final state as a fresh initial condition perturbs sensitivity
   estimates is unquantified. The current guidance is to avoid
   warm-starting inside differentiated solves (§17.7.1).

3. **No native forward-sensitivity ODE / adjoint.** The wrapper approach
   re-runs problem construction and relies on SciMLSensitivity's generic
   machinery rather than a hand-coded sensitivity/adjoint exploiting the
   structured Jacobian and its $\mathcal{O}(N)$ couplings (§15.9). A native
   adjoint that transposes the structured operators directly would be faster and
   lower-memory.

4. **Non-smoothness residue.** A handful of guards (`max`/`clamp`, the
   `A_ratio > eps()` switch, hard temperature/
   porosity clamps) remain non-differentiable at their boundaries (§17.7.2,
   §17.7.6). They are inactive in the interior of normal operation but can
   contaminate sensitivities evaluated exactly on the boundary; a fully smooth
   reformulation of these guards is not yet in place.

5. **Parameter handles are unregistered.** Because injection is by user closure
   (`p_inject`/`modify`), there is no enforced catalogue of "supported"
   parameters or automatic correctness check that the injected field actually
   flows differentiably into the residual; responsibility for a clean,
   non-mutating, $C^1$ injection rests with the user (§17.4.1).

6. **Global-sensitivity tooling lives in examples, not the core.** Sobol/PCE is
   implemented in `examples/sensitivity_analysis/`, not as a packaged function;
   the design and degree-selection choices (PCE order, BCa bootstrap,
   group-correlation handling) are reference recipes rather than a stable API.

Cross-references: the structured Jacobian whose transpose underlies the adjoint
is §15 (Temporal Integration & Structured Jacobian); the smooth kinetics gates
and depletion limiter are §6 (Reaction Kinetics); the mixture-property
derivatives are §5 (Effective Properties); the lateral-shrinkage law and $\chi$
dynamics are §10 (Volume Change & Lateral Shrinkage); the ALE geometry whose
dual-typed recompute makes geometry sensitivities possible is §11 (The ALE
Framework); and the conservation diagnostics used to sanity-check
sensitivity-driven calibrations are §16 (Conservation Diagnostics). Validation
of sensitivity-driven fits against experiment and reference codes is §18
(Verification & Validation).


---

# 18. Verification and Validation

This chapter sets out the verification and validation (V&V) methodology for
Pyrolysis.jl. Following Roache's now-standard distinction, *verification* asks
whether the discrete equations are solved correctly (a mathematics question,
internal to the code) and *validation* asks whether the right equations are
solved (a physics question, answered against experiment and trusted reference
codes). We catalogue the analytical and manufactured-solution checks, the
Jacobian-against-reference machinery, the runtime conservation ledgers used as
continuous self-checks, the spatial/temporal/tolerance convergence studies, the
experimental and code-to-code validation cases, the regime-validity diagnostics
that tell the user when a modeling approximation is being stressed, and the known
discrepancies that still require attention. Every tolerance, sign convention, and
discrete formula stated here has been traced to the implementing source; where a
formula appears it is the one the code evaluates, not an idealized textbook form.

All symbols conform to the unified nomenclature of §2; the three binding sign
conventions are reused throughout: (i) positive flux is transport in the $+z$
direction (bottom→top); (ii) the discrete divergence $(\mathrm F_R-\mathrm
F_L)\,A/V_i$ is positive for net outflow; (iii) heat of reaction $h>0$ is
endothermic, with the volumetric source $Q_\mathrm{rxn}=-\sum_r h_r r_r$.

---

## 18.1 The V&V philosophy

### 18.1.1 Verification versus validation

A pyrolysis solver is a chain of approximations: a continuum PDE system (§3) is
closed by constitutive laws (§4–§10), projected onto a finite-volume mesh (§13),
advanced by a stiff implicit integrator with a structured Jacobian (§15), and
read out through post-processing (§16). Errors can enter at every link, and they
have qualitatively different remedies. We therefore separate the two questions
explicitly.

- **Verification** — *Are we solving the equations right?* This is purely
  internal. It compares the code's discrete output against (a) closed-form
  solutions of the *same* discrete or continuous equations, (b) the
  finite-difference / automatic-differentiation reference of its own Jacobian,
  and (c) the conservation laws the discretization is designed to satisfy. A
  verification failure is a *bug* (a wrong sign, a missing factor, a stale cache,
  an incorrect limiter), not a modeling deficiency. Verification can, in
  principle, be made exact: a correctly implemented conservative scheme conserves
  mass to round-off, a correctly assembled Jacobian matches its AD reference to
  AD precision, and cell volumes recomputed from the node positions satisfy the
  geometric conservation law as an identity.

- **Validation** — *Are we solving the right equations?* This is external. It
  compares the code's predictions against measured cone-calorimeter mass-loss
  and surface-temperature histories, against thermogravimetric (TGA) curves, and
  against the published output of established codes (ThermaKin2Ds, Gpyro, the FDS
  solid phase). A validation discrepancy can be a modeling choice
  (one-temperature gas, gray-body radiation, single-reactant kinetics), a
  parameter-calibration issue, or a genuine missing physics. Validation is never
  exact: it is bounded by measurement uncertainty and by the irreducible
  approximations of a 1D condensed-phase model.

The two are kept methodologically distinct because conflating them is the most
common way to draw a wrong conclusion. A model that agrees beautifully with a
TGA curve after the activation energy has been tuned to it has been *calibrated*,
not validated; conversely, a perfectly verified code can still disagree with
experiment because its physics is incomplete. Throughout this chapter we keep the
ledger checks (verification) separate from the experimental comparisons
(validation), and we report the regime diagnostics (§18.7) that quantify how hard
a given run leans on an approximation.

### 18.1.2 The verification hierarchy in Pyrolysis.jl

Pyrolysis.jl layers its verification from cheap, fundamental checks to experimental
validation; of these, the Jacobian comparison can also be armed during a
production solve (§18.3.4):

1. **Unit-level analytical checks** (§18.2): isolated kernels — conduction,
   single-reaction kinetics, Beer–Lambert attenuation, the GCL, the P1 sign —
   compared against hand-derivable solutions.
2. **Jacobian verification** (§18.3): the production structured Jacobian compared
   against a whole-residual ForwardDiff reference, either as a one-shot call or
   as an in-solve callback.
3. **Conservation ledgers** (§18.4): the mass and energy trackers of §16, driven
   by the tests and by the verification cases that carry their own per-step
   ledgers; in a production solve `diagnostics = true` logs sensible-energy
   snapshots (§16.8.3).
4. **Convergence studies** (§18.5): spatial, temporal, and integrator-tolerance
   refinement, confirming that the discretization error shrinks at the designed
   rate.
5. **Code- and experiment-level validation** (§18.6): cone calorimeter, Douglas
   fir, pressure-treated wood, moisture transport, against measurements and
   reference-code output.

The first three run in the test suite (the verification suite of §18.1.3 also
checks convergence orders), while the convergence studies of §18.5 and the validation cases are scripted
examples shipped in `examples/`.

### 18.1.3 The verification suite: test mode and report mode

The exact-solution and manufactured-solution cases (catalogued below; the MMS
driver is described in §18.2) form a
verification suite that runs in two modes from a single case definition. The
whole suite has to stay cheap enough for `Pkg.test()`, while publishable
evidence — fine meshes, resolution ladders, figures — is far too expensive for it
and needs plotting dependencies the test environment should not carry.

- **Test (assert) mode** runs inside `Pkg.test()`. `test/runtests.jl` picks up the
  `.jl` files directly under `test/verification/` (non-recursively), and
  `runcases.jl` among them is the suite's entry point: it runs every case at its CI resolution, asserts each quantity of interest (QoI)
  against a pinned tolerance, runs a short spatial-convergence check where the
  case defines one, and asserts the wall time (after a discarded warm-up solve)
  against five times the case's budget. It plots nothing and writes no files.
  Setting `PYROLYSIS_VERIFICATION_PRINT_QOIS=1` prints every observed QoI and
  convergence rate, which is how tolerances are re-pinned.
- **Report mode** runs in the isolated environment `scripts/verification/` (its
  own `Project.toml` with CairoMakie, so the test manifest stays free of plotting
  dependencies; the same pattern as `scripts/profile/` and `test/quality/`):

  ```
  julia --project=scripts/verification scripts/verification/run_reports.jl [--cases=… --outdir=… --chapters=…]
  ```

  `run_reports.jl` runs each case over its report resolution ladder with its
  report solver settings and writes the artifact set under
  `<outdir>/assets/<case id>/`: overlay figures (solver values as markers over
  densely sampled exact curves), pointwise-error figures that show *where* the
  discretization struggles, a log–log convergence figure with the fitted observed
  order and an expected-order reference slope, and CSVs of every curve, the
  resolution ladder, the QoIs at the finest rung, the ladder QoI errors, and the
  solver configuration (`report_harness.jl`, with all plotting in `render.jl`).
  `chapters.jl` then writes one Markdown chapter per case family — problem
  statement, figures, error and order tables — plus an index dashboard. Report
  mode does not assert the QoIs (a failed solve still aborts the run): coarse ladder rungs may exceed the pinned tolerances,
  and the tables record within-tolerance yes/no instead. The default output
  directory is the committed `docs/guide/verification/`; the documentation build
  copies it into the published site, and the manual `verification-report`
  workflow regenerates it and opens a pull request with the diff.

**The shared case definition.** `test/verification/casedefs.jl` declares the
types both modes consume and depends on nothing from `Test`, so the report runner
can include it from its isolated environment. Each file in
`test/verification/cases/` defines a `VerificationCase` with

- identity and provenance: `id`, `title`, physics `tags`, `references`;
- `build(n_cells) -> PyrolysisProblem`, optional `solve_kwargs`, a `solve_fn`
  when the case solves through its own path (the MMS cases use the wrapped-RHS
  `solve_mms`), and `expected_retcode` (`:Success` by default, `:Terminated`
  when a callback stopping the solve is the answer);
- `qois(sol) -> Vector{QoIResult}`: each QoI's value, exact value, error, and
  tolerance, with a `comment` recording the tolerance provenance (tolerances are
  pinned at about three times the error observed when pinned);
- `ci::CIConfig`: the CI resolution, the wall-time budget, and an optional
  `ConvergenceSpec` (resolutions, error functional, expected and minimum observed
  order);
- `report::NamedTuple`: the report-mode keys — the resolution ladder
  `n_cells_ladder`, `curves(sol) -> Vector{ReportCurve}` (plain data:
  `CurveSeries` samples, each with an optional `exact` callable on the same
  abscissa from which the pointwise error is derived), `error_norms`,
  `expected_order`, report `solve_kwargs`, and optional temporal or parameter
  sweeps. A case with `report = (;)` is skipped in report mode.

The case file's leading comment banner is the problem statement the report
chapter prints, so the published chapter cannot drift from what the test runs.
Reference solutions are computed in the run at high precision (series summed to
machine tolerance, transcendental roots by bisection), never hard-coded; they
live in the shared library `test/verification/exact/` (series, roots, the MMS
driver), which carries its own self-tests.

**Case catalogue.** Each case family has one file in `test/verification/cases/`;
lettered members such as V-HC-02a–d share the family file and its report chapter.
Most families compare against an exact or manufactured solution; the burnout
family and V-CPL-06 check budgets and a pinned reference run instead.

| Family | Cases | Reference solution | What it verifies |
|---|---|---|---|
| Heat conduction | V-HC-01 – V-HC-05 | Carslaw–Jaeger erfc and eigenfunction series, resistance network, Kirchhoff transform | the convective-BC surface Newton (V-HC-01; Biot sweep V-HC-02), the constant-flux slab to its back-face plateau (V-HC-03), the harmonic-mean face conductivity at a layer jump (V-HC-04), $k(T)$ (V-HC-05) |
| Kinetics and surface radiation | V-RX-01 – V-RX-04 | Bateman chain, Khawam–Flanagan isothermal integrals, parallel-reaction splits, implicit radiative-cooling law | multi-step coupling, reaction order and the depletion limiter, two-product yields and heat superposition, `RadiativeBC` through the surface Newton |
| In-depth radiation | V-RAD-01 | Lautenberger eigenfunction series | Beer–Lambert absorption with conduction and a convective face, including the sub-surface temperature peak |
| Gas transport | V-GT-01, V-GT-02 | Gaussian plane source; steady compressible Darcy $P^2$ profile | the Fickian operator and mass conservation; the pressure closure, Darcy mobility and outflow BC |
| Coupled physics | V-CPL-01 – V-CPL-06 | Neumann two-phase Stefan solution, Landau ablation wave, exact shrink/swell endpoints, Frank-Kamenetskii slab, 0D TGA curve, discrete gas-enthalpy ledger | a latent-heat front, ALE surface regression, volume change and the lateral-shrinkage law, the reaction–conduction bifurcation, the BC/volume coupling of kinetics, $S_\mathrm{conv}$ (§18.4.3) |
| Manufactured solutions | V-MMS-01 – V-MMS-04 | the MMS driver of §18.2 | design order of the coupled Eulerian residual in space and time, free-stream preservation under prescribed motion, ALE order under prescribed and reaction-driven motion |
| Burnout | V-BO-01 – V-BO-05 | exact exponential thickness (V-BO-04); otherwise mass budgets and a pinned reference run | parked-cell depletion through complete burnout on ALE and fixed meshes (§14.3), the transmitted-flux warning, the interior-void stop (`expected_retcode = :Terminated`, V-BO-05) |

Kernel-level analytic checks are in `test/verification/test_analytical.jl`, which
`Pkg.test()` runs beside the suite: semi-infinite conduction, the Beer–Lambert
profile and energy split, isothermal and constant-heating-rate kinetics, and the
mass and energy ledgers. Two conventions matter when a case transcribes
literature parameters: $h>0$ is endothermic here, the opposite of ThermaKin,
Gpyro and FDS; and most exact-kinetics cases disable the depletion limiter
(`DepletionLimiter(enabled = false)`) and open the temperature gates (§18.2.2),
while V-RX-02d keeps the limiter on as a guard.

**Cross-code overlays.** Where a case has a published analog in another code,
that code's numerical output is committed under
`test/verification/reference/<code>/` with a provenance README (FDS and Gpyro runs
on the same problems; MaCFP-2 PMMA community submissions, which no case reads) and
read only inside a case's
report `curves`. A series marked `crosscode = true` is drawn as a muted context
layer in a separate `_vs_codes` variant of each affected figure, so the base
figure keeps the solver's own error on its own scale; for exact-solution cases the
other codes' *errors* against the same exact solution are compared.
`reference/timings.csv` feeds a wall-time comparison at every ladder rung, and a
case's ladder includes a rung matched to the reference run's resolution where
needed.

---

## 18.2 Method of manufactured solutions and analytical checks

The most rigorous verification tool is the *method of manufactured solutions*
(MMS): one chooses an analytic field, substitutes it into the governing PDE to
obtain the exact forcing, supplies that forcing to the discrete solver, and
confirms that the numerical solution converges to the chosen field at the
designed order. The verification suite ships a coupled MMS driver
(`test/verification/exact/mms.jl`: `MMSFields`, a ForwardDiff exact-source
generator, manufactured-BC builders, and the wrapped-RHS `solve_mms` path)
exercised by the V-MMS cases in `test/verification/cases/` — coupled Eulerian
MMS, ALE free-stream/GCL preservation, ALE MMS under prescribed mesh motion via
the `Workspace.w_override` hook, and ALE MMS under reaction-driven motion
(V-MMS-04, `MMSMotion`). Two constraints of that driver
follow from the formulation: the structured geometry Jacobian assumes
reaction-driven motion, so `w_override` runs use finite-difference Jacobians;
and because the scheme is a **solid-Lagrangian ALE** (mesh velocity ≡ solid
material velocity; no component carries a discrete advection operator —
§11.1.2, §11.9 item 5), every manufactured species source on a moving mesh
is assembled in the material-following frame, with the prescribed or the
reaction-driven mesh velocity as its $w$.

> **What MMS cannot see: wrong-PDE defects.** MMS verifies that the
> discretization solves the code's *own* target PDE at the designed order —
> the manufactured source is built from that same PDE, so any term the
> formulation wrongly includes (or omits) is self-consistently compensated
> in the source and invisible at every order and every resolution. A
> spurious ALE advection term, for example, would pass every MMS case
> untouched. Only checks with an *independent physical identity* the code
> does not get to define — an exact solution (the Landau traveling wave of
> V-CPL-02), a conservation ledger (V-CPL-03, V-CPL-06) — can catch a
> wrong-PDE defect. This is why the suite pairs MMS (discretization order)
> with exact-solution and conservation anchors (formulation correctness):
> the two families have disjoint blind spots.

In addition, Pyrolysis.jl supports the *analytical-solution* special case of
MMS for each isolated physics kernel, in which a parameter regime is chosen
so that an exact solution of the governing equation is already known and no
manufactured forcing is needed. These are the building blocks; each is
independently verifiable and each anchors a specific term of `compute_rhs!`.

### 18.2.1 Pure conduction

With all gas transport, reaction, and radiation terms disabled, the energy
equation of §7 collapses to the linear heat equation

$$
\rho c_p^\mathrm{eff}\,\frac{\partial T}{\partial t}
   = -\nabla\!\cdot\mathbf q_\mathrm{cond},
\qquad
\mathbf q_\mathrm{cond} = -k\,\nabla T .
$$

For constant $\rho c_p^\mathrm{eff}$ and $k$ on a fixed (Eulerian) mesh this is
$\partial_t T = \alpha_\mathrm{th}\,\partial_{zz}T$ with thermal diffusivity
$\alpha_\mathrm{th}=k/\rho c_p^\mathrm{eff}$, which admits a wealth of
closed-form solutions: the error-function semi-infinite response to a step
surface temperature, the Fourier-series response to a constant surface flux, and
the decaying-cosine eigenmode

$$
T(z,t) = T_\infty + (T_0-T_\infty)\,
   \cos\!\Big(\tfrac{\pi z}{2L}\Big)\,
   \exp\!\Big(-\alpha_\mathrm{th}\,\tfrac{\pi^2}{4L^2}\,t\Big)
$$

on $[0,L]$ with an insulated bottom ($\partial_z T|_{z=0}=0$, an `AdiabaticBC`)
and a fixed top ($T(L,t)=T_\infty$, a `DirichletThermalBC`). The single decay
constant $\lambda=\alpha_\mathrm{th}\pi^2/4L^2$ makes this an exacting check on
the whole conduction pathway: the harmonic-mean face conductivity, the divergence
operator, the diffusivity assembly, and the implicit time integration must all be
correct for the numerical decay rate to match.

The discrete face flux verified here is (§7.1, `conductive_flux`,
`src/Physics/heat_conduction.jl`)

$$
q_f = -k_f\,\frac{T_R-T_L}{d_L+d_R},
\qquad
k_f = \frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
$$

(distance-weighted harmonic mean). On a *uniform* mesh, $d_L=d_R=\tfrac12\Delta
z$ and $k_f$ reduces to the ordinary harmonic mean $2k_Lk_R/(k_L+k_R)$, and for
constant $k$ the scheme is the textbook second-order central difference; the
numerical decay rate then matches $\lambda$ to $O(\Delta z^2)+O(\Delta t^p)$ where
$p$ is the integrator order. The harmonic mean (rather than arithmetic) is
verified by a *two-layer* slab with a conductivity jump: the steady-state
solution has a continuous flux and a kinked temperature profile, and only the
harmonic (resistance-in-series) face value reproduces the analytic interface
temperature.

**Implementation note.** Interior conductive fluxes are assembled by
`compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`), which zeroes
boundary faces (the BC solver overwrites them); the divergence
$(\nabla\!\cdot\mathbf q)_i=(q_{f_R}-q_{f_L})A/V_i$ is `divergence_1d!`
(`src/Discretization/finite_volume.jl`). The RHS assembles
$dT_i/dt = -(\nabla\!\cdot\mathbf q)_i/\rho c_p^\mathrm{eff}_{,i}$ in step 8 of
`compute_rhs!`.

### 18.2.2 Single-reaction TGA (0D-idealized cell)

A thermogravimetric experiment heats a small, thermally-thin sample at a constant
rate $\beta = dT/dt$ and records the residual mass fraction. In the limit of a
single, well-mixed cell with a uniform temperature ramp, all transport vanishes
and the species balance for a one-step first-order decomposition $A_1\!\to\!
\text{products}$ becomes the 0D ODE

$$
\frac{d\xi_1}{dt}
   = -A_1\exp\!\Big(-\frac{E_1}{R_g T}\Big)\,\xi_1 ,
\qquad T(t)=T_0+\beta t .
$$

Changing the independent variable to temperature and integrating gives the
classic non-isothermal conversion law

$$
\frac{\xi_1(T)}{\xi_1(0)}
   = \exp\!\left[
       -\frac{A_1}{\beta}\int_{T_0}^{T}
         \exp\!\Big(-\frac{E_1}{R_g\,T'}\Big)\,dT'
     \right],
$$

whose integral is the Arrhenius temperature integral (approximated by the
Coats–Redfern or Senum–Yang series, or evaluated numerically). The conversion
$\chi_\mathrm{conv}=1-\xi_1/\xi_1(0)$ traces the familiar sigmoid; its inflection
(maximum reaction rate) sits at the temperature where
$d^2\chi_\mathrm{conv}/dT^2=0$, giving the reference relation between the peak
temperature $T_p$ and the kinetic pair $(A_1,E_1)$,

$$
\frac{E_1\,\beta}{R_g\,T_p^2}
  = A_1\exp\!\Big(-\frac{E_1}{R_g\,T_p}\Big),
$$

which is exactly the FDS-style recipe for recovering $(A_1,E_1)$ from a measured
$T_p$ and peak rate (see §6.9 and §18.6.2). Reproducing $T_p$, the peak rate, and
the residual fraction in a single-cell Pyrolysis run verifies the Arrhenius
kernel, the species stoichiometry, and the depletion limiter together.

Two implementation subtleties must be accounted for when this check is made
*quantitatively* exact rather than merely qualitative:

- **Mass-concentration basis.** The rate law is evaluated on mass concentration
  $\xi_j$ (kg m$^{-3}$), not mole fraction or normalized conversion:
  $r_i=A_i\exp(-E_i/R_gT)\prod_j\xi_j^{\,n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ (Eq. 1 of §6,
  `reaction_rate`, `src/Physics/kinetics.jl`). For a first-order single
  reactant this is linear in $\xi_1$ and the analytic solution above holds with
  $\xi_1(0)$ the *initial concentration*.
- **Smooth gates and the depletion limiter perturb the tails.** The hard
  Heaviside gates $H(T-T_\mathrm{min})H(T_\mathrm{max}-T)$ are replaced by tanh
  ramps $\tfrac12(1+\tanh(T-T_\mathrm{min}))\cdot\tfrac12(1+\tanh(T_\mathrm{max}
  -T))$ (§6.2, `kinetics.jl`), which leak $\sim\!10$ K outside the nominal
  window; and the depletion limiter multiplies the rate by
  $\tanh(\xi_1/\text{threshold})$ (default threshold $1.0$ kg m$^{-3}$, §6.3,
  `kinetics.jl`). To compare against the bare-Arrhenius analytic solution one
  must either set $T_\mathrm{min},T_\mathrm{max}$ far outside the active range
  (so the gates read $\approx 1$) and keep $\xi_1$ well above the limiter
  threshold throughout, or disable the limiter via
  `DepletionLimiter(enabled=false)`. With those settings the numerical curve
  matches the analytic conversion to integrator tolerance; the residual gap *is*
  the smoothing error and is itself a quantitative verification of the gate
  width.

### 18.2.3 Beer–Lambert transmitted fraction

The in-depth radiation model (§8.3) integrates the Beer–Lambert law
$\partial I/\partial z=-\alpha_\mathrm{eff}I$, giving
$I(z)=I_\mathrm{surface}\exp\!\big(-\!\int_z^L\alpha_\mathrm{eff}\,dz'\big)$ and
the exact cell-averaged volumetric source

$$
q'''_i = I_\mathrm{in}\,\frac{1-\exp(-\tau_i)}{\Delta z_i},
\qquad
\tau_i = \alpha_i\,\Delta z_i,
$$

assembled by `compute_beer_lambert_source!` in a top-to-bottom sweep
($i=n\!:\!-1\!:\!1$) so that absorbed flux upstream correctly attenuates the
intensity reaching downstream cells (`src/Physics/radiation.jl`). Two exact
energy-balance identities verify this kernel:

- **Total absorbed flux.** Summing the volumetric source times cell thickness,

  $$
  \sum_i q'''_i\,\Delta z_i
     = \alpha_\mathrm{surf}\,I_\mathrm{surface}\big(1-e^{-\tau_\mathrm{total}}\big),
  \qquad
  \tau_\mathrm{total}=\sum_i\alpha_i\Delta z_i,
  $$

  where $\alpha_\mathrm{surf}$ is the surface transmissivity of §8.3.2
  ($\varepsilon_{\text{eff}}$ on the Kirchhoff path; 1 when the BC's explicit
  absorptivity is already folded into the intensity), is a telescoping sum that
  is *exact at the discrete level* — independent of cell count — because each
  cell deposits $I_\mathrm{in}(1-e^{-\tau_i})$ and passes
  $I_\mathrm{in}e^{-\tau_i}$ to the next. The solution reports it as
  `extras.volumetric_absorption`, and the sum of the deposited source over the
  cells must equal this analytic value.

- **Transmitted fraction.** The flux escaping the back face is
  $\alpha_\mathrm{surf}\,I_\mathrm{surface}\,e^{-\tau_\mathrm{total}}$, reported as
  `extras.radiation_transmitted` (§8.3.2). Conservation requires
  $\text{absorbed}+\text{transmitted}=\alpha_\mathrm{surf}\,I_\mathrm{surface}$
  to round-off — the remaining $(1-\alpha_\mathrm{surf})\,I_\mathrm{surface}$ is
  reflected at the surface and never enters the medium — which is the strongest
  single check on the Beer–Lambert path.

The transparent-cell guard ($\alpha>10^{-20}$, else $q'''_i=0$,
`radiation.jl`) is verified by the $\alpha\to0$ limit, where $1-e^{-\tau}\to
0$ and all flux is transmitted. The opposite limit $\tau_1\gg1$ in the surface
cell deposits the entire entering flux
$\alpha_\mathrm{surf}\,I_\mathrm{surface}$ in the top cell — consistent with
`SURFACE_ABSORPTION`, which applies the same
$\varepsilon_{\text{eff}}\,I_\mathrm{surface}$ at the face (§8.2), providing a
cross-check between the two non-P1 radiation models.

### 18.2.4 Geometric Conservation Law to machine precision

The single most decisive ALE verification is the Geometric Conservation Law
(GCL, §11.5). A moving-mesh scheme that does not reproduce the swept cell volume
exactly injects *phantom* mass and energy sources proportional to the volume
error; over hundreds of stiff steps these can swamp the genuine reaction-driven
mass loss. The continuous 1D GCL is $dV_i/dt = A\,(w_R-w_L)$ (plus the
area-change term $|z_R-z_L|\,dA/dt$ under a lateral shrinkage law). Pyrolysis.jl
satisfies it *geometrically*: the node positions are state ($dz_i/dt=w_i$) and
every residual evaluation recomputes $V_i=A(t)\,|z_R-z_L|$ from them, so the
volume the scheme uses is always the swept one and no separately evolved volume
can drift (§11.5.2). With the cumulative mesh velocity this gives
$dV_i/dt=V_i\theta_i$ identically, and with the dilation correction a uniform
field with no sources and no fluxes stays uniform.

The verification case V-MMS-02 asserts exactly that: a spatially uniform state
with no reactions and zero-flux boundaries, on a mesh compressed to 70 % of its
thickness through the `Workspace.w_override` hook with the uniform-stretch
velocity $w^*(z,t)=(\dot s/s)\,z$, must stay uniform to machine precision. The
node trajectories $z_i(t)=z_i(0)\,s(t)$ are known exactly, so the top-cell-center
position doubles as a check on the prescribed-motion path. V-MMS-03 then measures
the observed order of a smooth manufactured solution under prescribed smooth
mesh motion. A drift in V-MMS-02 is a moving-mesh transport defect that the
smooth-field order study could hide inside its truncation error.

### 18.2.5 P1 sign test

The P1 quasi-steady radiation model (§8.4, Experimental) solves
$\nabla\!\cdot(D_\mathrm{rad}\nabla G)=\alpha G-4\sigma\alpha T^4$ for the
incident field $G$ and feeds the net source

$$
S_\mathrm{rad} = \alpha\,(G-4\sigma T^4)
$$

into the energy equation. The sign is the part most easily gotten wrong, and it
is verified directly: $S_\mathrm{rad}>0$ when $G>4\sigma T^4$ (the local field
exceeds the blackbody level → net absorption → heating), $S_\mathrm{rad}<0$ when
$T^4$ dominates (net emission → cooling), and $S_\mathrm{rad}=0$ at local
equilibrium. Two corollary checks pin the boundary and the volumetric coupling:

- **Marshak equilibrium.** With no external source ($I_\mathrm{ext}=0$) and the
  surface in equilibrium with its surroundings ($T_s=T_\infty$, so
  $G_\mathrm{surf}=G_\mathrm{ext}=4\sigma T_\infty^4$), the Marshak boundary flux
  $q_\mathrm{BC}=\frac{\varepsilon}{2(2-\varepsilon)}
  (G_\mathrm{ext}+4I_\mathrm{ext}-G_\mathrm{surf})$ must vanish identically. This
  verifies the exact $\frac{\varepsilon}{2(2-\varepsilon)}$ coefficient (Modest
  Eq. 16.48) rather than the incorrect simplification $0.5\varepsilon$. A second
  regression pins the **external-flux weight**: the $I_\mathrm{ext}$ response is
  $2\varepsilon/(2-\varepsilon)\cdot I_\mathrm{ext}$ — the same half-range weight
  the ambient term carries — and a cold optically-thick absorber at
  $\varepsilon=0.9$ absorbs the continuum-P1 fraction $4c/(1+\sqrt3\,c)\approx
  0.96$ of the incident flux (a weight-1 injection would absorb only
  $\approx 0.59$).
- **No volumetric emissivity factor.** The volumetric emission is the blackbody
  value $4\sigma\alpha T^4$ with *no* $\varepsilon$ factor (gray-body local
  equilibrium, $\varepsilon=1$ in the interior); $\varepsilon<1$ enters *only* at
  the boundary through the Marshak term. A test that compares interior emission
  against $4\sigma\alpha T^4$ verifies this convention.

The Newton solve for $G$ uses an allocation-free Thomas (TDMA) factorization with
relative convergence $\|R\|_\infty<\text{tol}\cdot\|G\|_\infty$ (G fields are
$10^4$–$10^6$ W m$^{-2}$, so absolute tolerance would be unphysical) and a
non-negativity clip $G\ge0$; these are verified against analytic optically-thick
and optically-thin limits in `test/radiation/test_p1_validation.jl`. Because P1
lives in `Pyrolysis.Experimental` and uses a *separate* residual
(`RadiationP1`), the production `solve` does not route
`radiation_model=P1_QUASI_STEADY` through the unified residual path; it rejects
the setting outright at keyword validation with an error directing the user to
the `Pyrolysis.Experimental.RadiationP1` primitives, which are invoked directly
(`_validate_solve_kwargs`, src/Solver/solve.jl). P1's verification is therefore confined to
the Experimental test set (§8.12).

---

## 18.3 Jacobian verification

### 18.3.1 Why and against what

The implicit integrators used for these stiff systems (default `KenCarp4`, an
SDIRK method) solve a Newton system at every stage with the iteration matrix
$W=I/(\gamma_\mathrm{RK}\Delta t)-J$ (§15.5). A wrong Jacobian $J$ does not change
the *fixed point* of the Newton iteration — so a converged step is still a
correct solution of the discrete ODE — but it degrades or destroys convergence,
inflating step rejections and, at worst, stalling the solve on an Arrhenius
transient. The production Jacobian is the *structured* operator assembly of §15
(a sparse block $A$ plus O($N$) prefix/suffix/rank-1 couplings); its analytical
and AD-per-operator pieces are optimizations of a quantity that must equal the
true Jacobian of the residual. The reference against which "true" is defined is
`DenseAD_ForwardDiff`: a single ForwardDiff pass over the *entire* residual
vector, which is correct by construction (it differentiates exactly the function
the integrator evaluates) and slow (dense, full-width). Verification compares the
two.

### 18.3.2 `verify_jacobian_at`

The standalone comparator is

```julia
verify_jacobian_at(spec_active, def, ws, u, t;
    reference = DenseAD_ForwardDiff(), rtol = 1e-6, atol = 1e-10) -> NamedTuple
```

(`src/Jacobian/verify.jl`). It assembles the active backend's Jacobian and the
reference Jacobian at the same state $(u,t)$ and reports the worst absolute and
relative deviation and the $(\text{row},\text{col})$ where it occurs.

### 18.3.3 Strict versus masked comparison

The comparison mode is selected by the `is_approximate` trait of the active
backend (§15.6, §15.10):

- **Strict / support-complete** (exact backends, `is_approximate == false`):
  every entry of the dense reference is compared against the active backend's
  dense view, with the backend's view zero outside its declared support, under
  the **per-entry** criterion

  $$
  |\Delta J_{ij}| \;\le\; \texttt{atol} + \texttt{rtol}\cdot
  \max\big(|J^{\rm ref}_{ij}|,\ \max\nolimits_j |J^{\rm ref}_{i,:}|\big).
  $$

  The row-magnitude floor keeps honest roundoff at structural zeros of a
  stiff row from failing (a $10^{-9}$ residue in a $10^6$-scale Arrhenius row
  passes), while a wholly missing block in a *small*-magnitude row — the ALE
  z-rows are ~10 orders below the ξ-rows — contributes to its own row's
  reference scale and fails regardless of problem size. (A global
  max-$|J^{\rm ref}|$ normalization would make pass/fail scale-dependent —
  small-magnitude rows such as the ALE z-rows could lose an entire block
  behind the Arrhenius rows' magnitude — which is why the criterion is
  per-entry with a row floor.) A reference non-zero that falls outside the backend's
  support surfaces as a deviation — a *missing-entry* failure. The verifier
  separately reports `max_abs_dev_in_support` (genuine disagreement on a
  claimed entry) and `max_abs_dev_outside_support` (magnitude of reference
  entries the support omits). The exact backend, `Structured(DirectSolve())`,
  must pass this regime.

- **Masked** (approximate backends, `is_approximate == true`): entries outside
  the active backend's support are zeroed in the difference, so only in-support
  disagreement is asserted, under the coarse global-scale criterion
  $\texttt{max\_abs\_dev}\le\texttt{atol} \lor
  \texttt{max\_abs\_dev}/\max_{ij}|J^{\rm ref}_{ij}|\le\texttt{rtol}$. This is
  the regime for `Structured(ApproxSolve())`, which *deliberately* drops the
  structured ALE couplings (the cumulative mesh-velocity prefix sum, the z-row
  prefix blocks, the Beer–Lambert suffix, the $\theta_A$ rank-1 back-reaction,
  the geometry $z$-columns) under the `:ale_reduced` sparsity profile. Its
  in-support values legitimately deviate by the dropped coupling magnitudes,
  so masked comparison is a sanity check
  rather than an exactness gate; the dropped couplings are by design (they
  cost accuracy in the Newton step, not in the converged solution).

### 18.3.4 In-solve verification callback

The same comparator can be installed as a `DiscreteCallback` that fires when the
integrator crosses requested sample times:

```julia
solve(problem; verify_jacobian = (rtol=1e-6, atol=1e-10, sample_times=[0.0, 1.0, 10.0]))
```

(`src/Solver/solve.jl`). The callback compares the active backend against the
`DenseAD_ForwardDiff` reference at each sample time and `@warn`s on a failure but
*does not throw*: a single transient outlier (e.g. an extreme Arrhenius spike
where the analytical and AD paths differ at the last AD digit) should not abort
an otherwise-good production run. This lets a long calibration sweep be run with
verification armed at negligible extra cost — the reference Jacobian is built
only at the sample times.

### 18.3.5 Parked states

Depletion handling (§14.3) parks cells without changing the state size, so the
structured Jacobian is the one the integrator uses through burnout and the
verifier covers it: `test/jacobian/test_park_jacobian.jl` builds the specs on
the fresh problem, parks surface, interior and bottom cells in every mode
(Eulerian, ALE, Darcy, Beer–Lambert, variable cross-section), and checks
`Structured(strategy = DirectSolve())` strictly against the dense reference, the
zero ghost rows, the ghost node rows equal to the surface node row, and the
live block of the parked Jacobian against the compacted problem's — for
`ApproxSolve` too, since that strategy drops couplings by design and is not
compared against the dense reference on a reacting state.

---

## 18.4 Conservation-based verification

The conservation ledgers of §16 serve as verification where the tests and the
verification cases drive them per step (`record_substep!`, or a case's own
per-step ledger callback); a production `solve(problem; diagnostics=true)` logs
only sensible-energy snapshots, with the ledger columns NaN (§16.8.3). The checks
divide cleanly into those that must close by construction (so a failure is a bug)
and those that quantify an approximation (so a large value is a regime warning,
not a bug).

### 18.4.1 Mass ledger

The per-component mass tracker (`MassConservationTracker`,
`src/Diagnostics/conservation_checks.jl`) maintains
$m_j(t)=\sum_{i\in\mathrm{active}}\xi_j^{(i)}V_i$ (masked: only active cells,
§16.3.1) and accumulates boundary efflux
$\Delta m_{\mathrm{out},j}=-(J_{\mathrm{top},j}+J_{\mathrm{bot},j})\,A\,\Delta t$
with the convention that positive $J$ is *inflow*, so positive flux *reduces*
efflux. The balance check is

$$
\max_j
  \frac{\big|m_{\mathrm{current},j}
          -(m_{\mathrm{initial},j}-m_{\mathrm{efflux},j})\big|}
       {|m_{\mathrm{initial},j}|}
  < 10^{-10}
$$

(`check_mass_conservation`, default tol $10^{-10}$). When a component starts at
zero or near-zero mass (e.g. char before any reaction), the check switches from
relative to *absolute* error to avoid division by a vanishing baseline. For ALE
runs with a time-varying cross-section the area in the efflux term is the
trapezoidal substep average $A=\tfrac12(A^\mathrm{pre}+A^\mathrm{post})$, which
keeps the contribution symmetric and consistent with the geometric volume update.
A mass-balance failure indicates a flux-assembly or callback bug.

### 18.4.2 Energy ledger

Energy verification uses the sensible-energy ledger of §16.4–16.5. The stored
quantity is $E(t)=\sum_i\rho c_p^\mathrm{eff}\,T_i\,V_i$ with the mixture
$\rho c_p^\mathrm{eff}$ (every component, gas included, §16.4), balanced against
the *exact terms the RHS assembles*:

$$
E(t) = E(0)
  + \int_0^t\!\Big[
      q_\mathrm{cond,boundary}\,A
      + Q_\mathrm{rxn}
      + S_\mathrm{conv}
    \Big]\,d\tau .
$$

The residual $r=E(t)-E(0)-\sum(\text{sources})$ should close to integrator
tolerance ($\sim10^{-6}$ relative, `check_energy_conservation` default tol
$10^{-6}$) because it is built from exactly the terms `compute_rhs!` evaluates. A
failure here is a bookkeeping bug. The one caveat (§16.7, §18.8) is that this
ledger *cannot* close exactly when composition varies in time, because it omits
the $\int T\,\partial(\rho c_p)/\partial t\,dV\,dt$ heat-capacity-change term of its
$\rho c_p\,T$ proxy density (§16.7) and the ALE volume-change term; the residual
then carries a small, composition-rate-dependent offset that is a known
incompleteness of the *audit*, not of the *solver*.

**The boundary flux.** The ledger uses $q_\mathrm{cond,boundary}$ (the conducted
flux from the surface Newton solve). The surface balance has no gas-enthalpy
term, so it equals the applied boundary flux $q_\mathrm{BC,applied}$
(the externally driven flux, which may include radiation/convection; §16.6).
Crucially, the discrete ledger
*audits bookkeeping, not double-counting*: if the RHS itself double-counts a
mechanism, the discrete ledger will still close, because it is built from the same
(buggy) terms. Detecting a double-count requires an *independent* check — an
exact solution or a ledger built from outside the RHS (§18.4.3).

### 18.4.3 Independent closure of the gas-enthalpy coupling

The gas enthalpy is carried entirely by the always-on volumetric advection source

$$
S_{\mathrm{conv},i}
  = -\sum_g c_{p,g}(T_i)\,\bar J_{g,i}\,
      \Big(\frac{\partial T}{\partial z}\Big)_i
$$

(§7.4, `compute_gas_convective_source!`), and the surface balance carries no
gas-enthalpy term (§7.4.3). A surface term $\sum_g\dot m_g\,\Delta h_g(T_s,T_\infty)$
added on top of $S_\mathrm{conv}$ would count the advected enthalpy twice, yet
the energy ledger of §18.4.2, built from the same RHS terms, would record both
and still close. The coupling is therefore pinned by a ledger built independently
of the RHS: the verification case V-CPL-06 (a char-layer column whose gas forms in
a cool band and vents through a hot layer, $h=0$ so only the sensible ledger acts)
checks

- **V-CPL-06a (Eulerian):** the discrete summation-by-parts identity
  $\int S_\mathrm{conv}\,dV\,dt = -c_g\big(\int J_\mathrm{top}T_{f,\mathrm{top}}A\,dt
  - \int\sum_i(\mathrm{gen}_i-\dot\xi_{g,i})\,T_i V_i\,dt\big)$ — the energy
  $S_\mathrm{conv}$ removes from the field equals the energy the species ledger
  says the gas picks up between formation and exit — and the physical reheat
  reference $-c_g\sum\mathrm{gen}\,(T_\mathrm{exit}-T_\mathrm{form})\,V$ to
  quasi-steady accuracy;
- **V-CPL-06b (ALE):** the same column with the fuel band collapsing; the ALE
  energy budget closes with $E=\sum\rho c_p T V$ over every component, and the
  mesh-relative $T$-advection accumulator is exactly zero (§11.2.2).

Its ledger runs in a per-accepted-step callback around the pure `residual!` with
trapezoidal time quadrature (a rectangle rule leaves an $O(\Delta t\,\dot X)$
residual that is pure quadrature). V-CPL-06a shows that $S_\mathrm{conv}$ alone
removes the whole gas reheat, so a surface gas-enthalpy term would remove it a
second time; a mesh-velocity $T$-advection term would appear in V-CPL-06b as a
sink of about 66 kJ per kg of volatiles. The absence of the surface term is
asserted directly by `test/boundary_conditions/test_boundary_temperature.jl`
(testset "Outward gas flux enters only through blowing"): under a prescribed-flux
BC, the solved surface temperature and conducted flux are identical with and
without a gas flux across the face.

---

## 18.5 Convergence studies

A scheme can be conservative and have a correct Jacobian yet still be *under-
resolved*. Convergence studies confirm that the discretization error shrinks at
the designed rate as the mesh, time step, or integrator tolerance is refined, and
they let the user pick a resolution that is converged for the quantity of
interest (QoI) at acceptable cost.

### 18.5.1 Spatial convergence (uniform and stretched)

The cell-centered FVM of §13 is second-order accurate in space for smooth
solutions: the harmonic-mean face conductivity and the central-difference
gradient are each $O(\Delta z^2)$ on a uniform mesh, and the distance-weighted
generalizations preserve formal second order on smooth, slowly-stretched meshes.
The error of a QoI $f$ should therefore scale as $f(\Delta z)-f_\mathrm{exact}
=C\,\Delta z^p+\cdots$ with $p\approx2$. With an exact solution unavailable, the
finest mesh is used as the reference (Richardson-style self-convergence): the
shipped `examples/convergence_douglas_fir.jl` solves a Douglas-fir cone case at a
sequence of cell counts (finest first = reference), then reports per-QoI relative
errors

$$
\Delta f(n) = \frac{|f(n)-f_\mathrm{ref}|}{|f_\mathrm{ref}|}\times100\%
$$

for scalar QoIs — peak mass-loss rate, time-to-ignition, final mass, peak surface
temperature — and *curve-level* $L_2$ relative errors

$$
\frac{\lVert f(n)-f_\mathrm{ref}\rVert_2}{\lVert f_\mathrm{ref}\rVert_2}\times100\%
$$

for the full MLR$(t)$, $T_s(t)$, and back-face $T(t)$ histories. Plotting
$\Delta f$ against cell count on log–log axes recovers the slope $-p$; a slope
near $-2$ confirms second-order spatial convergence, while a shallower slope
flags a first-order contamination (e.g. the first-order upwind in the Darcy–Fick
advection of §9.3, or the surface-cell radiation deposition).

Stretched meshes (`create_stretched_mesh`, geometric ratio, §13.9) are tested
separately because the boundary-refinement they provide near the heated surface
is where the steep thermal and reaction fronts live; the distance-weighted
harmonic mean is what preserves accuracy there (a 1–5% flux-accuracy gain near
boundaries versus the naive arithmetic mean, §7.1). A convergence study on a
stretched mesh verifies that the non-uniform-mesh corrections in the gradient and
face-value reconstructions (§13.6) are correct: an *incorrect* non-uniform
treatment shows up as a loss of order on the stretched run that is absent on the
uniform run.

### 18.5.2 Temporal convergence

The default integrator `KenCarp4` is a 4th-order SDIRK method; an integrator-order
study fixes the mesh, sweeps the maximum step (via `dtmax` or by tightening
tolerances), and confirms the temporal error scales at the method's design order
for smooth intervals. In practice the *Arrhenius transients* — the sharp onset of
decomposition as a reaction gate opens — dominate the temporal error and force
the adaptive controller to small steps; the tanh-smoothed gates and the
tanh depletion limiter (§6.2–6.3) are essential here precisely because they keep
the RHS and its Jacobian differentiable, so the error estimator and the Newton
solver behave. A temporal study that holds the mesh fixed and tightens `dtmax`
isolates the time-integration error from the spatial error.

### 18.5.3 Integrator-tolerance study

Because the step size is chosen adaptively to meet `abstol`/`reltol` (defaults
$10^{-8}$ and $10^{-6}$, `src/Solver/solve.jl`), the practical temporal-error
control knob is the tolerance pair. A tolerance study sweeps `reltol` (and
`abstol` in proportion) and tracks a QoI to its tolerance-converged plateau; the
convergence example uses tightened tolerances (`abstol=1e-9`, `reltol=1e-7`) for
its reference runs so that the *spatial* error dominates the comparison and is not
masked by temporal error. The energy-ledger residual (§18.4.2) provides an
independent confirmation that the tolerance is tight enough: it closes to roughly
`reltol`, so a residual that fails to shrink as `reltol` is tightened
indicates a bookkeeping issue rather than under-resolution. The three knobs —
$\Delta z$, integrator order/`dtmax`, and tolerance — are studied separately so
that the dominant error source is unambiguous; only when all three are converged
is a quoted QoI a property of the *model* rather than of the *discretization*.

---

## 18.6 Validation against experiments and reference codes

Validation compares converged Pyrolysis.jl predictions against measured data and
against the published output of established condensed-phase pyrolysis codes. The
notation of the reference codes is mapped to ours on first use; the comparisons
draw on the governing-form, kinetics, transport, and radiation contrasts
documented throughout §3–§9.

### 18.6.1 Cone calorimeter: HRR/MLR/$T_s$

The cone calorimeter is the canonical bench-scale validation experiment: a sample
is exposed to a constant radiant flux (typically 25–75 kW m$^{-2}$) and the
mass-loss rate (MLR), surface temperature, and oxygen-consumption heat-release
rate (HRR) are recorded. Pyrolysis.jl, being a condensed-phase solver, predicts
MLR and $T_s$ directly and an HRR *proxy* via the effective heat of combustion
times the gas-generation rate; the gas-phase flame is not modeled. The shipped
`examples/06_cone_calorimeter.jl` (PMMA, 25 kW m$^{-2}$, top convection +
re-radiation, back-face re-radiation and convection with an impermeable back,
single-step Arrhenius, temperature-dependent properties) is the reference workflow; its surface-
absorption variant `06_cone_calorimeter_surface_absorption.jl` exercises the
opaque-material radiation path.

The validation QoIs are the peak MLR and its timing, the residual (char) mass
fraction, the time-to-ignition surrogate (a surface-temperature or
gas-generation threshold), and the $T_s(t)$ history. Two implementation points
materially affect quantitative agreement and are themselves verified:

- **Surface temperature in cell-centered FV.** Temperature lives at cell
  *centers*, but a $T^4$ radiation BC must be evaluated at the *surface*. Using
  the boundary-cell-center temperature directly produces an $\sim$20% energy
  imbalance for strong radiative BCs; `solve_boundary_temperature` instead Newton-
  iterates for the true $T_s$ satisfying $q_\mathrm{applied}(T_s)=
  q_\mathrm{conducted}(T_s,T_\mathrm{cell})$ and stores it for read-out
  (`get_surface_temperature`, §12.2). Validation against measured $T_s$ requires
  this solved value, not the cell-center temperature.
- **MLR gauge normalization.** For a fixed cross-section MLR is per unit area; for
  a shrinking variable cross-section (WithChi) the *gauge* MLR is normalized to
  the initial area, $\mathrm{MLR}_\mathrm{gauge}=\mathrm{MLR}_\mathrm{local}\cdot
  A(t)/A_0$ (§10, `src/Solver/solve.jl`), matching the cone gauge convention so
  the prediction is comparable to a measurement reported per *initial* sample
  area.

### 18.6.2 Douglas fir and the convergence case

`examples/07_douglas_fir_cone.jl` is the most fully-featured validation case:
multi-step wood decomposition kinetics, surface-absorption radiation (wood is
opaque), Darcy–Fick gas transport (species molar masses enter the pressure
closure; all gases share one diffusive transfer coefficient, §5.11), ALE
tracking of the shrinking sample, and the `ApproxSolve` reduced Jacobian for a
$\sim$8× speedup over the full structured Jacobian. The companion
`convergence_douglas_fir.jl` (§18.5.1) establishes the mesh resolution at which
the validation QoIs are converged, so that any residual disagreement with
experiment is attributable to physics or parameters rather than to
under-resolution. Validation here checks the MLR$(t)$ shape (the characteristic
two-peak signature of hemicellulose/cellulose and lignin decomposition), the
char yield, and the back-face temperature rise (a sensitive probe of the in-depth
conduction and the effective conductivity mixing rule, §5.5).

The FDS-style kinetic recovery of §18.2.2 is the practical link between TGA data
and the cone model: $(A_1,E_1)$ are recovered from a measured TGA peak temperature
and peak rate via the inflection relation, then carried into the cone run. This is
a *calibration* step (the kinetics are fit to TGA), so the cone comparison that
follows is a genuine *validation* only insofar as the cone data were not used in
the fit — a distinction §18.1.1 insists on.

### 18.6.3 Pressure-treated wood and moisture transport

`examples/08_pressure_treated_wood_cone.jl` (and `ptw_comparison.jl`) validate the
Darcy-driven transport path: a permeable charring solid develops an internal
pressure $P=\sum_g\xi_jR_gT/M_j/\varphi$ (§9.4), and the Darcy velocity
$v_g=-(\kappa_\mathrm{face}/\mu_\mathrm{face})\nabla P$ (§9.5) advects the
pyrolysis gases outward. The validation QoI is the MLR augmentation and the
internal-pressure history relative to a Fickian-only run; the Darcy–Fick demo
`08_darcy_fick_demo.jl` isolates the pressure-driven transport for a controlled
check. Moisture transport (`moisture_hos_comparison.jl`) validates the
state-dependent heat-of-sorption path: the evaporation "reaction" uses a
`StateDependentProperty` heat $h(T,\xi)$ that depends on the local moisture
concentration (§6.6), and the moisture front, the evaporation-plateau in $T(t)$
(the constant-temperature stall while water boils off), and the dried-mass profile
are the comparison targets.

### 18.6.4 Comparison to ThermaKin2Ds, Gpyro, and FDS

Code-to-code comparison (`examples/literature_comparison/`) is run alongside
experimental validation because it isolates *modeling-convention* differences
from *numerical* ones. The three reference codes differ from Pyrolysis.jl in
documented, specific ways; mapping their notation to ours:

- **ThermaKin2Ds.** Maps its stoichiometric coefficients $\theta_i^{\,j}$ to our
  $\nu_{i,j}$ and its component concentrations to our $\xi_j$. ThermaKin's energy
  equation (its Eq. 17) is Pyrolysis.jl's energy equation term for term: the
  storage term $\sum_j\xi_j c_{p,j}$ summed over *all* components including gas,
  the convective term, and a deformation term that the ALE mesh frame absorbs
  (§3.3.2). ThermaKin's convective energy term
  $c_g(J_g\cdot\partial T/\partial x)$ is identical in form to our $S_\mathrm{conv}$
  (the central-difference gradient matches ThermaKin's discretization, §7.4), so
  the convective coupling is a clean cross-check.
  ThermaKin supports *bimolecular* reactions ($\xi_\mathrm{COMP1}\,
  \xi_\mathrm{COMP2}$); Pyrolysis.jl enforces single-reactant kinetics
  ($N_\mathrm{Rin}=1$, §6.5), so multi-reactant ThermaKin schemes must be
  re-expressed as reaction sequences before comparison.

- **Gpyro.** Maps its pre-exponential $Z$ to our $A_i$, its permeability $K$ to
  our $\kappa$, its volume fractions $X_\alpha$ to our $v_j$, and its volumetric
  rate $\dot\omega'''$ to our $r$. Gpyro is fully conservative (control-volume
  conservation of $\bar\rho\Delta z$, $\bar\rho Y_i\Delta z$, gas mass, and
  energy) and uses volume-fraction property averaging with power-law
  temperature-dependence; Pyrolysis.jl's WEIGHTED/PARALLEL/SERIES/BRUGGEMAN
  mixing rules (§5.5) span and generalize Gpyro's averaging. Gpyro's Darcian mass
  flux $\dot m''=-(K/\nu)(\partial P/\partial z-g\rho_g)$ matches our Darcy form
  (§9.5) up to the buoyancy term $g\rho_g$, which Pyrolysis.jl omits (negligible
  for thin condensed-phase slabs). Gpyro selects swelling/shrinkage mechanisms by
  solid-fraction state; Pyrolysis.jl's lateral-shrinkage law $A(t)=A_0\cdot
  \mathrm{law}(\bar\chi)$ (§10.5) is the 1D analog. Gpyro's solver is a
  fully-implicit backward-Euler / Patankar linearization with TDMA; Pyrolysis.jl
  uses higher-order SDIRK with a structured Jacobian — so code-to-code agreement
  on a steady or slowly-varying case validates the *physics* independent of the
  time integrator.

- **FDS solid phase.** Maps its absorption $\kappa_s$ to our $\alpha$ and its
  pre-exponential $Z$ to our $A_i$. FDS solves 1D Fourier conduction with
  Arrhenius kinetics ($\rho^n O_2^n T^n$ rate forms) and a two-flux
  (Schuster–Schwarzschild) in-depth radiation model; Pyrolysis.jl's Beer–Lambert
  and P1 models bracket the FDS two-flux treatment. FDS carries no in-depth gas,
  so its energy equation is the condensed-phase limit of ours; the FDS
  TGA-reference kinetic recovery (peak temperature and
  rate → $A,E$) is the same recipe used in §18.2.2 and §18.6.2. FDS's surface
  energy balance $-k\,\partial T/\partial x=q_c''+q_r''$ is the multi-mechanism
  analog of our `solve_boundary_temperature` Newton balance (§12.2).

In all three comparisons, agreement *after* matching conventions (gas treatment,
radiation model, kinetic form) is the verification that the shared physics is
implemented consistently; *residual* disagreement after convention-matching is
attributable to the documented modeling differences, not to error.

---

## 18.7 Regime-validity checks

A model is only valid in the regime where its approximations hold. Pyrolysis.jl
exposes runtime diagnostics that tell the user when a run is leaving the
validated regime — these are validation aids, not verification checks, because
they assess the *applicability of the equations*, not the *correctness of the
solve*.

### 18.7.1 Optical-thickness and P1 validity

The optical thickness $\tau_\mathrm{total}=\int\alpha\,dz=\sum_i\alpha_i\Delta z_i$
selects the appropriate radiation model and is computed by `optical_thickness`
(Experimental). The validated regimes (§8): $\tau>10$ optically thick (P1
excellent, surface absorption acceptable); $1<\tau<10$ intermediate (P1 good,
Beer–Lambert recommended); $0.1<\tau<1$ optically thin (Beer–Lambert); $\tau<0.1$
transparent (use `NO_RADIATION`). `check_p1_validity` flags a run whose optical
thickness places it outside the P1-accurate band, and `radiation_energy_balance`
confirms the radiation source integrates to the absorbed-minus-emitted flux. The
absorbed-plus-transmitted identity of §18.2.3 is the underlying conservation
check; the validity query layers a *regime* judgement on top of it.

### 18.7.2 Pore-gas thermal equilibrium

The energy equation stores the pore gas at the matrix temperature (§7.3). That
one-temperature assumption holds when the gas equilibrates with the pore walls
within a cell, which is the case for pore-scale heat-transfer lengths far below
the cell size — every validated regime here (1 atm, pore sizes well below the
mesh spacing, gas residence times of milliseconds). A run with coarse open pores
or very fast, hot throughput, where the gas does not equilibrate with the
walls, is outside the model: there is no two-temperature gas option (§3.9), and
the prediction should be treated with caution.

### 18.7.3 Péclet number and upwind diffusion

The Darcy–Fick species transport uses first-order upwinding for the advective
part (§9.3). First-order upwind introduces numerical diffusion proportional to the
cell Péclet number $\mathrm{Pe}=|v_g|\Delta z/\lambda$; when $\mathrm{Pe}\gg1$
the scheme is robust but the artificial diffusion can smear the gas-composition
front. The convergence study (§18.5.1) exposes this as a first-order contribution
to the spatial error of pressure-driven cases. The practical regime check is to
confirm that the cell Péclet number is moderate over the resolved front — if it is
large, the user should refine the mesh near the front (a stretched mesh) until the
upwind diffusion is below the physical diffusion. The conduction analog is benign
(the conductive flux is central-differenced and unconditionally second-order), so
this regime concern is specific to the advective gas transport.

### 18.7.4 Depletion-threshold sensitivity

Two thresholds govern near-zero-concentration behavior, and a validated run should
be insensitive to both within their physically-negligible range. The *depletion
limiter* threshold (default $1.0$ kg m$^{-3}$, min $0.01$, §6.3) sets where the
tanh rate-suppression begins; the *surface-depletion* threshold (default 5% of
each cell's own initial thickness, `depletion_threshold=0.05`, §14.3,
`src/Solver/solve.jl`) sets when a thinned (ALE) or emptied (fixed mesh) cell is
parked. A sensitivity check varies
each threshold by an order of magnitude and confirms the QoIs (peak MLR, char
yield) move by less than the tolerance; a strong dependence means the threshold is
acting in a physically-active region rather than only suppressing
numerical stiffness near true depletion, and the value should be lowered. These
are AD-relevant: the limiter is tanh-smooth (so sensitivities are well-defined,
§17.6), but a park is a discrete rewrite of the state (whose size stays fixed)
and is *not* differentiable — a caveat for sensitivity studies (§17.7, §18.8.3).

---

## 18.8 Known discrepancies and re-validation needs

Honest V&V records the open items. The following are documented incompletenesses
that affect, or require re-checking of, validated results.

### 18.8.1 Surface gas-enthalpy convention and calibrated examples

The energy balance follows the $S_\mathrm{conv}$-only convention (ThermaKin
Eq. 17 / FDS): $S_\mathrm{conv}$ carries the escaping-gas enthalpy and the
surface balance has no gas-enthalpy term (§7.4.3), since a surface term would
count that enthalpy a second time. A balance that adds a surface term alongside
$S_\mathrm{conv}$ artificially depresses the surface temperature, so parameter
sets calibrated in a solver that does so read a *higher* $T_s$ and mass-loss rate
here (§1.8, item 7). The fitted kinetic and property parameters of the bundled
calibrated examples — the inverse-analysis cases and the cone runs — *need
re-validation* against data, together with the flame heat-flux feedback and blowing that a
recalibration should include (§1.8.1). This is a validation-state issue, not a
verification failure.

### 18.8.2 Composition-varying discrete energy residual does not close exactly

The energy ledger (§18.4.2) is missing the
$\int T\,\partial(\rho c_p)/\partial t\,dV\,dt$ heat-capacity-change term of its
$\rho c_p\,T$ proxy density (heat capacity changing at constant temperature as
composition reacts — §16.7) and an ALE volume-change term. As a result it cannot close to integrator tolerance when
composition varies in time; the residual carries a small composition-rate-
dependent offset. A ledger that closes exactly would require per-substep
accumulation of $\rho c_p(t_i)\,\Delta T_i\,V(t_i)$ per cell, which is not
implemented. This is a limitation of the *audit*, not of the solver: the solver's
energy update is correct; the bookkeeping formula used to *check* it is
incomplete. Until it is extended, a non-zero residual during active
reaction should be interpreted against this known offset, not as a bug.

### 18.8.3 Sensitivities across park events

Depletion handling is a sequence of discrete state rewrites (§14.3.4, §17.7).
The Jacobian and the verifier run through them (§18.3.5), but smooth
sensitivity analysis does not: forward sensitivities are valid only between
park events and adjoint sensitivity through an event is unsupported. Runs that
need smooth sensitivities must stay on horizons without park events.

### 18.8.4 Heuristic and unvalidated closures

Several constitutive closures are heuristic and carry their own validation debt:
the Lennard-Jones parameter estimates for Chapman–Enskog diffusivity
($\sigma\approx3.5+0.5\log_{10}(M/0.016)$,
$\varepsilon/k\approx150+100(M/0.044)$, §5.10) are explicitly flagged unvalidated;
the lateral-shrinkage laws (§10.5) are empirical rather than first-principles; and
the single global radiation model is a simplification. Validation cases that depend
sensitively on these closures should report that dependence.

### 18.8.5 Stale-flux and Marshak-coefficient regressions

Two silent-failure modes carry standing regression guards. The incident
radiation flux is re-sampled at *every* RHS evaluation with no constancy
detection (§8.3.4), so a time-varying flux can never be frozen at its initial
value; and the Marshak BC uses the exact
$\frac{\varepsilon}{2(2-\varepsilon)}$ coefficient rather than
$0.5\varepsilon$, with the external flux at its half-range weight
$2\varepsilon/(2-\varepsilon)$ rather than 1 (§18.2.5). Both would be silent
errors if regressed, so the equilibrium P1 sign test, the cold-absorber
absorbed-fraction test, and a modulated-flux comparison
($q_0(1+\sin 2\pi t)$) are kept as permanent checks.

---

## 18.9 Reproducibility

Verifiable results must be *reproducible*, which requires a controlled software
environment and an automated quality gate.

### 18.9.1 Quality tooling

Three static/dynamic quality tools run in isolated environments (so the
package's own Manifest is not pinned to their dependencies):

- **JuliaFormatter** — enforces a single canonical source format, so diffs are
  semantic and review is reliable.
- **Aqua.jl** — audits package hygiene: undefined exports, method ambiguities,
  stale dependencies, unbound type parameters, piracy. A clean Aqua run is a
  precondition for trusting that the public API surface (UG-16) is consistent.
- **JET.jl** — static type-stability and dispatch analysis. Because the hot paths
  rely on `@generated` unrolling over heterogeneous component tuples and on
  end-to-end `Tu`-typed (Float64 or `ForwardDiff.Dual`) caches for AD-cleanliness
  (§4.6.5, §7.7), a JET regression that introduces a type instability or a dynamic
  dispatch would both slow the solve and threaten the AD/sensitivity paths.

These run via dedicated project environments (`--project=scripts`,
`--project=test/quality`); the main test manifest is intentionally separate
because the quality tools' pins are fragile.

### 18.9.2 Test layout

The test suite mirrors the source module tree (`test/` mirrors `src/`), and
`runtests.jl` auto-walks the directories so a new test file is picked up without
manual registration. The full suite is run via `Pkg.test()`. The verification
machinery of this chapter has direct test coverage:

- Jacobian verification: `test/perf/test_verify_jacobian.jl` drives
  `verify_jacobian_at` over the canonical-configuration × backend matrix;
  per-operator tests (`test/jacobian/test_operator_*.jl`) check each operator's
  Jacobian row against AD.
- Conservation: the tracker tests exercise the mass and energy ledgers, and the
  verification case V-CPL-06 pins the gas-enthalpy coupling, and
  `test/boundary_conditions/test_boundary_temperature.jl` asserts that the surface
  balance carries no gas-enthalpy term (§18.4.3).
- The verification suite: `test/verification/runcases.jl` runs every case of
  `test/verification/cases/` in assert mode (§18.1.3).
- Radiation: `test/radiation/test_p1_sign.jl` (sign), `test_p1_validation.jl`
  (analytic limits), `test/jacobian/test_operator_radiation.jl` (Beer–Lambert
  upper-triangular structure).
- Allocation discipline: `test/perf/test_zero_allocation_rhs.jl` asserts the RHS
  assembly allocates zero bytes — a reproducibility-relevant performance
  invariant, since allocations in the hot path would couple results to GC timing.

### 18.9.3 Determinism and environment

Reproducible numbers require a pinned `Manifest.toml` (the integrator,
linear-solver, and AD package versions all affect the last digits) and a fixed
choice of integrator, tolerances, Jacobian backend, and linear solver. The
defaults — `KenCarp4` with `AutoFiniteDiff`, the pattern-chosen factorization
(`KLUFactorization` or `UMFPACKFactorization`), `abstol=1e-8`/`reltol=1e-6`,
`Structured(DirectSolve())` (`ApproxSolve` under ALE) — are the reference
configuration; deviations (e.g. `ApproxSolve` on a fixed mesh, `ILUGMRESFactorization`
for large ALE runs)
trade a controlled amount of accuracy for performance and should be re-checked
against the reference configuration with `verify_jacobian` armed when first
adopted on a new problem.

---

## 18.10 Limitations and open issues

The V&V coverage above has gaps that a user should weigh:

1. **MMS coverage is bounded by the solid-Lagrangian formulation.** The
   coupled MMS driver (§18.2, `test/verification/exact/mms.jl`) closes the
   design-order convergence gap for the fully coupled residual, but mesh
   motion in the MMS cases is either reaction-driven or prescribed via the
   verification-only `w_override` hook; under prescribed (non-material)
   motion the species sources are assembled in the material-following frame
   and the Jacobian falls back to finite differences (§11.9 item 5). General
   ALE (a species advection operator for mesh motion decoupled from the
   solid) remains out of scope.

2. **Convergence reference is self-consistent, not exact.** The convergence
   studies use the finest mesh as reference (Richardson self-convergence); they
   confirm *grid convergence* and the *rate*, but not the absolute error against
   a true solution. For QoIs with a known analytic value (the conduction
   eigenmode, the Beer–Lambert transmitted fraction) the absolute error can be
   bounded; for the full nonlinear cone case it cannot.

3. **Audit ledgers are incomplete for varying composition and ALE volume change**
   (§18.8.2). The discrete energy residual does not close exactly during active
   reaction; the missing heat-capacity-change and volume-change terms must be
   accounted for when interpreting it.

4. **Depletion is non-differentiable** (§18.8.3): the structured Jacobian and its
   verifier run through park events (§18.3.5), but sensitivities are valid only
   between them, so the most strongly-shrinking cases (deep char, ablation) have
   the narrowest sensitivity horizons.

5. **P1 radiation is Experimental and separately verified** (§18.2.5, §8.12): it
   has its own residual outside the production `solve` path, so its V&V does not
   ride the unified test matrix and must be exercised through the Experimental
   tests.

6. **Validation is parameter-entangled.** Several cone/inverse cases were
   calibrated against the data they are compared to, and those calibrations need
   re-fitting (§18.8.1). Re-fitting against held-out data
   — keeping calibration and validation data disjoint — is the outstanding
   validation task, and is the discipline §18.1.1 requires for a comparison to
   count as validation rather than calibration.

7. **Heuristic closures carry validation debt** (§18.8.4): the LJ diffusivity
   estimates and the empirical shrinkage laws are unvalidated or simplified, and cases that lean on them inherit that
   uncertainty.


---

