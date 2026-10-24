# 1. Introduction, Scope, and Modeling Philosophy

This chapter establishes what Pyrolysis.jl is, the class of physical problems it solves, and the modeling decisions that govern every later chapter. Pyrolysis.jl is a one-dimensional (through-thickness) condensed-phase pyrolysis solver: it integrates the coupled heat-conduction, multi-species reaction-transport, and (optionally) gas-pressure equations that describe a solid heated to thermal decomposition. We position the code relative to the three established reference solvers it is most often compared to — ThermaKin2Ds, Gpyro, and the FDS solid-phase pyrolysis model — and we state the binding conventions (coordinate orientation, flux sign, divergence sign, heat-of-reaction sign) and the energy equation that the remainder of the Technical Reference Guide derives and exploits. Everything asserted here is grounded in the actual implementation; where a textbook form and the code differ, the code is authoritative and the deviation is flagged.

---

## 1.1 Purpose and intended audience

Pyrolysis.jl simulates the thermal decomposition (pyrolysis) of condensed-phase materials subjected to an external heat flux — the canonical fire-science and materials-degradation problem of a solid that heats up, undergoes one or more Arrhenius-rate reactions, releases gaseous volatiles, and (in the charring/intumescent cases) changes its conductivity, porosity, and physical dimensions as it does so. The numerical core is a cell-centered finite-volume discretization on a face-based mesh that may move with the solid (Arbitrary Lagrangian–Eulerian, ALE; the mesh velocity is the solid velocity plus a correction, $w = v_s + c$). The semi-discrete equations form a mass-matrix differential-algebraic system $\mathsf M\,d\mathbf u/dt = \mathbf f(\mathbf u, t)$, advanced by the stiff Rosenbrock integrator `Rodas5P` with a Jacobian obtained by colored sparse forward-mode automatic differentiation of the residual on a sparsity pattern built from the mesh connectivity.

This document — the **Technical Reference Guide** — derives and justifies all of the mathematics and physics: the continuum governing equations (§3, Governing Equations of the Continuum Model), the constitutive closures (§4, Materials & Property Functions; §5, Effective Properties & Mixing Rules; §6, Reaction Kinetics; §8, Thermal Radiation), the discretization (§7, Heat Conduction, Gas Convection & Energy Assembly; §13, Finite-Volume Discretization), the moving-mesh framework (§10, Volume Change & Lateral Shrinkage; §11, The ALE Framework), the boundary treatment (§12, Boundary & Initial Conditions), burnout (§14, Cell Deactivation and Burnout), the temporal integration and Jacobian machinery (§15, Temporal Integration, the DAE, and the Colored Jacobian), and the conservation, sensitivity, and verification apparatus (§16 Conservation Ledgers and the Audit, §17, §18). It is a *reference*, not a tutorial: it states the equations as the code implements them, with file-and-function attribution, and explains why each modeling choice was made. The companion **User Guide** (Part II of the documentation) teaches workflow — how to build materials, set boundary conditions, choose the physics options, and run the solver.

The reader is assumed to be comfortable with the heat equation, the Arrhenius law, finite-volume / finite-difference discretization, and the basic vocabulary of stiff ODE integration. No prior knowledge of ThermaKin2Ds, Gpyro, or FDS is required, but readers familiar with any of them will find explicit notation-mapping and contrast subsections throughout.

### 1.1.1 How to read the two guides

The Technical Reference and the User Guide share a single, authoritative symbol table — the **unified nomenclature** reproduced verbatim in §2 (Nomenclature). Every chapter conforms to it. Three conventions in that table are *load-bearing* and are restated below in §1.5 because misreading any of them inverts a sign in a governing equation:

1. the spatial-coordinate orientation and the meaning of "top" and "bottom",
2. the positive-flux direction and the outflow-positive divergence operator, and
3. the heat-of-reaction sign convention.

A typical reading path for a first-time reader is: §1 (this chapter) → §3 (the PDE system) → the User Guide quick-start, returning to §4–§9 as constitutive questions arise, and to §13–§16 when numerical behavior needs explaining. A developer extending the code should additionally read §11 (ALE) and §15 (the DAE and the colored Jacobian) in full.

---

## 1.2 Scope: the class of problems modeled

Pyrolysis.jl targets the **one-dimensional, through-thickness** decomposition of a planar sample. The spatial domain is a column of material of thickness $L$ and cross-sectional area $A$; the only resolved spatial coordinate is the axial coordinate $z$ (see §1.5.1). The unresolved lateral plane carries the area $A$, which may itself vary in time to represent lateral (radial) shrinkage or swelling (§10). A line may also be read in cylindrical or spherical coordinates, $z$ then being the radius (the mesh's `coordinates`, §13).

The state of the material is described, in every finite-volume cell $i$, by:

- the **component masses** $m_{j,i}$ [kg] of every chemical component $j \in \{1,\dots,N_C\}$ — the integrated unknowns; the **mass concentrations** $\xi_{j,i} = m_{j,i}/V_i$ [kg·m⁻³] that every constitutive law reads are derived from them — and
- one **energy state**: the cell temperature $T_i$ [K] (the default), or the cell's sensible-enthalpy inventory $H_i = \sum_j m_{j,i}\,h_j(T_i)$ [J] (`energy_state = :enthalpy`), from which the temperature is recovered (§1.4.2).

The choice of mass concentration $\xi_j$ (rather than mass fraction $Y_j = \xi_j/\rho$) as the constitutive species variable is deliberate: it makes the species conservation law a plain transport equation with a reaction source, with no embedded density derivative, and it is the form the reaction kinetics naturally produce (§6). Integrating the cell masses themselves makes the discrete mass balance exact: a cell's mass changes only by its reaction sources and its face fluxes, and a moving cell's change of volume needs no correction term (§3.2, §11.3.4). Further unknowns are introduced only when the corresponding physics is enabled — on a moving mesh (`use_ale = true`) the node positions $z_n$ (differential) and node velocities $w_n$ (algebraic, held by the cells' volume balances; §11.3); under `BEER_LAMBERT` the face intensities $I_f$ (algebraic; §8); under a lateral-shrinkage law the conversion-history pair $(G_i, m_{\mathrm{dry},i})$ of every cell and the cross-section area $A$ with its lateral strain rate $\theta_A$ (algebraic; §10).

### 1.2.1 Intended applications

The solver is designed for, and validated against, the following experiment and material classes (the shipped `examples/` directory exercises each):

- **Cone calorimetry.** A sample is irradiated at a controlled flux (commonly 25–75 kW·m⁻²) while surface temperature, mass-loss rate (MLR), and a heat-release-rate proxy are recorded. The exposed top face takes a radiative-plus-convective boundary condition; the back face is typically adiabatic or a lumped loss, and impermeable. Examples `06_cone_calorimeter.jl` (5.8 mm of PMMA under 25 kW·m⁻² on a moving mesh, with cell deactivation) and `07_douglas_fir_cone.jl` (12 mm of charring Douglas fir under 75 kW·m⁻² on a moving mesh, with Darcy–Fick transport and surface absorption) are the reference configurations.
- **Idealized thermal-analysis experiments.** Thermogravimetric analysis (TGA), differential scanning calorimetry (DSC), and microscale combustion calorimetry (MCC) can be idealized as a thin, thermally-thin sample with a prescribed heating rate. In the limit of a single well-mixed cell, the species equations reduce to the 0-D Arrhenius mass-loss ODE that TGA inversion assumes; this idealization underpins the single-reaction verification case in §18 (Verification & Validation).
- **Charring, non-charring, and intumescent solids.** Non-charring polymers (PMMA) leave essentially no residue; charring solids (wood) form an insulating char layer with conductivity and porosity that differ sharply from the virgin material; intumescent coatings swell by orders of magnitude. The component/reaction/property-function machinery (§4) plus the swelling factor $\gamma_j$ and lateral-shrinkage law (§10) span all three.
- **Pressure-treated and wet wood / moisture transport.** Hygroscopic materials carry liquid water that evaporates endothermically and migrates as vapor under pressure gradients. Darcy–Fick transport (§9, Gas-Phase Mass Transport & Pressure; selected when a condensed component is permeable) and the `StateDependentProperty` heat-of-sorption closure (§4) handle these; examples `08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl`, and `moisture_hos_comparison.jl` are the reference cases.

### 1.2.2 What is **not** modeled

The scope boundary is as important as the scope itself. The following are explicitly out of scope:

- **Multi-dimensional geometry.** Full solves are one-dimensional in $z$. The cross-section $A$ is a scalar gauge, not a resolved 2-D/3-D field; lateral shrinkage is a column-global area scaling driven by an empirical law (§10), not a resolved lateral transport. Meshes of two and three dimensions (`box_mesh`, `cell_mesh`) exist at the residual level — a potential closure of the mesh motion, conductivity tensors and element death are implemented and tested there (§11.3.5, §14.5) — but `setup` refuses gas transport, in-depth radiation, the variable cross section, the equidistribution correction and cell deactivation on them, and they are not a supported path for full solves.
- **Full two-temperature gas.** The gas phase is *not* given its own energy equation with a separate gas temperature. Gas is assumed to be in local thermal equilibrium with the matrix at the cell temperature $T$: it stores sensible heat like every other phase, and its flux relative to the skeleton carries its sensible enthalpy across the cell faces (§1.4.3, §7).
- **Multi-reactant reactions.** Each `Reaction` supports exactly one reactant ($N_{\text{in}}=1$; a reaction with more is rejected by problem validation; `_apply_reaction_source!` in `src/Physics/kinetics.jl` also guards source evaluation). The reactant is consumed with unit mass stoichiometry and products are explicit yields summing to 1.0 (§6). Bimolecular or competitive oxygen-sensitive schemes (which ThermaKin and Gpyro support natively) must be approximated as sequences of single-reactant steps. Extending the reactant tuple to a 3-ary `(idx, order, mass_stoich)` form is an open item (§1.8.1).
- **Pressure- or oxygen-dependent kinetics, and distributed activation energy.** The Arrhenius pre-exponential $A_i$ is temperature-independent, the rate has no explicit pressure or oxygen factor, and there is no built-in distributed-activation-energy (DAEM) tooling (§6).

These exclusions are not accidental: they are the modeling simplifications that make the 1-D solver fast, type-stable, and differentiable with respect to its state, and they are shared (to varying degrees) by the reference codes discussed in §1.6.

---

## 1.3 Design philosophy: a generalized model with optional physics

The organizing principle of Pyrolysis.jl is a **single generalized model with opt-in physics layers**. Rather than maintaining separate solvers for "simple conduction," "charring wood," and "pressure-driven char," the code expresses one general residual and selects which physical terms are active through the *types* of the definition it reads. `setup` turns a `PyrolysisProblem` and the solve keywords into a `CoreDefinition` whose state layout, kernel tuple, gas-transport model and mesh closure are type parameters, so the residual specializes at compile time without runtime branching and without method explosion: there is exactly one residual method (`test/jacobian/test_allocations.jl` asserts it), and a physics option that is off contributes no kernel, no rows and no work.

### 1.3.1 What `setup` decides

`setup(problem; ...)` (file `src/Integration/setup.jl`; `solve` calls it) makes each of the following choices once, from the problem and the keywords:

| Choice | Decided by | What it changes |
|---|---|---|
| Energy state | `energy_state` (`:temperature`, the default, or `:enthalpy`) | the energy row of every cell carries $T_i$ or $H_i$ (§1.4.2, §3.3.5) |
| Mesh motion | `use_ale` | node positions (differential) and node velocities (algebraic, the balance rows of `BalanceRowClosure`) join the state; without it every node is pinned |
| Mesh correction | `correction`; by default a rule | `EquidistributionCorrection` where the mesh moves, cells are deactivated and some cell, converted completely, would fall below the deactivation threshold; the inventories then cross the faces at the relative velocity (`RelativeFluxKernel`, §11.7) |
| Radiation | `radiation_model` | under `BEER_LAMBERT` the face intensities (algebraic) and `BeerLambertKernel`; otherwise the incident flux enters the surface energy balance (§8) |
| Gas transport | the material and the mass boundary conditions | none (no gas component, or nothing can move it), Fickian (some component carries a gas-transfer coefficient or a boundary condition passes mass), Darcy–Fick (a condensed component is permeable); `GasTransportKernel` (§9) |
| Cross section | the material's `lateral_shrinkage_law` | the conversion-history pair per cell and the global rows of $A$ and $\theta_A$ (`LateralShrinkage`, §10); needs `use_ale = true` |
| Burnout | `deactivation`, `threshold`, `band` | `CellDeactivation`: an event that removes cells whose condensed mass falls below the threshold (§14) |

The kernels run in a fixed order, each present only when its physics is: `ReactionKernel`, `GasTransportKernel`, `BeerLambertKernel`, `ConductionKernel`, `ThermalBoundaryKernel`, `RelativeFluxKernel` (`_kernels`, `src/Integration/setup.jl`). The order matters where one kernel reads another's output: the surface energy balance of the thermal boundary kernel reads the outward gas fluxes of the gas transport for its blowing models and for the enthalpy those fluxes carry (§7.6, §15.3).

The layout's type parameters are read by `@inline` accessors with zero runtime cost:

```julia
@inline is_moving(::StateLayout{NC,E,Mov}) where {NC,E,Mov} = Mov
```

The entry point `residual!(du, u, def::CoreDefinition, t, ws::CoreWorkspace)` (file `src/Assembly/residual.jl`, reached through `Pyrolysis.Internal`) runs the same sequence for every configuration — the state mapping `derive!`, the node corrections, the kernels in order, the mesh closure's rows, the external sources, the row assembly (§15.3). Its buffers are built once per problem, and once per dual-number type the Jacobian evaluates the residual at (`CoreWorkspace`, `src/Assembly/definition.jl`); a configuration without a moving mesh carries no position or velocity rows, one without Darcy–Fick transport evaluates no porosity, pressure or mobility.

> **Implementation note.** `_gas_transport` chooses the transport model, `_installed_correction` the default correction, `_kernels` the kernel tuple, and `_require_supported` refuses a problem that needs physics its mesh does not support (`src/Integration/setup.jl`). The layout, closure and kernel types are not exported at the top level; they are reached through `Pyrolysis.Internal` or their submodules (`Pyrolysis.States`, `Pyrolysis.Closures`, `Pyrolysis.Kernels`).

### 1.3.2 Opt-in physics layers

The same opt-in philosophy governs the constitutive layers. Each of the following is off (or set to its cheapest closure) by default and is enabled only when the problem requires it:

- **Radiation in depth.** `NO_RADIATION` and `SURFACE_ABSORPTION` both route the incident flux of a `RadiativeFluxBC` through the surface energy balance, with no volumetric deposit, at $O(1)$ cost (the two give identical residuals); `BEER_LAMBERT` resolves exponential in-depth absorption through face intensities that are algebraic unknowns of the state, which keeps the Jacobian banded. The models are mutually exclusive (§8).
- **Darcy gas flow.** With Fickian transport, gas migrates only by (volume-fraction-gradient) diffusion. When a condensed component carries a permeability, `setup` selects Darcy–Fick transport, which adds pressure-driven bulk flow $v_g=-(\kappa/\mu)\nabla P$ with the ideal-gas pressure closure (§9); the state mapping then also evaluates every cell's porosity, pressure and mobility.
- **ALE moving mesh.** A fixed mesh suffices for thermally-thin or non-receding problems. With `use_ale = true` the mesh moves with the solid: the node positions are differential unknowns, the node velocities algebraic unknowns held by each cell's volume balance, and the cell volumes follow from the positions, so the Geometric Conservation Law holds by construction (§11). Where cells must burn out, an equidistribution correction moves the mesh relative to the solid so that a vanishing layer's cells thin together (§11.7).
- **Variable cross-section.** A material's `lateral_shrinkage_law` enables lateral shrinkage/swelling through the column conversion $\bar\chi$ and an empirical area law (§10).
- **Burnout.** With `deactivation = true`, a cell whose condensed mass falls below the fraction `threshold` ($\varepsilon_m$, default $10^{-2}$) of its initial condensed mass is deactivated, together with the near-threshold cells outward of it (the band), in an event located inside the step; its inventories go to the live cell beneath, and the state keeps its size (§14, User Guide §10). There is no mesh refinement.

The cost of an optional layer is therefore paid only when it is used, and the modeling assumptions of the default (cheapest) configuration are exactly those of the reference codes (§1.6).

### 1.3.3 Implementation goals stated up front

Three implementation properties are treated as first-class requirements throughout the code, and they drive many modeling micro-decisions documented in later chapters:

1. **Type stability.** Every hot-path function is written so the Julia compiler can infer concrete types, using heterogeneous tuples for the component/reaction lists and `@generated` or recursive tuple functions to unroll loops over $N_C$ components and $N_R$ reactions (e.g. `_gas_partial_pressure_sum`, `_heat_capacity_sum`, `dry_pyrolysis_gas_indices`). The number of components $N_C$ and reactions $N_R$ are static type parameters of `Material{NC,NR,…}` (§4), and every physics option is a type parameter of the definition the residual reads.
2. **Allocation-free assembly.** The residual and the colored Jacobian evaluate with zero heap allocation per call. All scratch storage lives in the derived state and the kernel contributions of the `CoreWorkspace`, with twins for each dual-number element type built on first use; the state mapping `derive!` is the only reader of the state vector and the row assembly the only writer of its derivative (§15.2). A dedicated test (`test/jacobian/test_allocations.jl`) asserts zero bytes for one residual and one Jacobian evaluation on every configuration it builds, in both energy forms.
3. **AD-compatibility (differentiability).** The entire residual — including the geometry from the node positions, the surface-temperature Newton solve, the caloric inversion of the enthalpy form, and the property evaluations — propagates `ForwardDiff.Dual` numbers, and the colored Jacobian is its exact derivative. This requires every non-smooth construct to take the form of a C¹ surrogate or by a guard that is the identity on physical states: the kernels evaluate the property laws at the *trial temperature*, the cell temperature smoothly clamped into the property bracket [200, 3000] K with 1 K joins; masses enter the property laws through a one-sided positive part; the floors on the heat capacity, the conductivity, the porosity and the viscosity are smooth (§13); hard temperature gates become $\tfrac12(1+\tanh(\cdot))$ ramps; and the depletion limiter multiplies each reactant's rate by a C¹ factor that is exactly zero at and below zero concentration (§6.3). The one deliberate exception is an upwind side, which is read from the value of its velocity, so an upwind flux has a kink where its velocity changes sign (§9, §11.6). Differentiability is what makes the exact Jacobian (§15) possible, and the same smoothness keeps the finite-difference sensitivities of §17 well-defined.

These three goals are not merely engineering niceties: AD-compatibility in particular *constrains the physics formulation*, because every modeling term must be expressible as a smooth function of the state.

---

## 1.4 The energy equation

The single most important modeling decision in Pyrolysis.jl is the form of the energy equation. The code solves the energy balance of a one-temperature porous medium: the pore gas is in local thermal equilibrium with the matrix, **every component stores sensible heat**, and the gas flux relative to the skeleton carries its sensible enthalpy across the cell faces. The conserved quantity is the sensible-enthalpy content of the material, and the integrator carries either the cell temperature or the cell's sensible-enthalpy inventory (§1.4.2). This subsection states the balance, derives the temperature equation from it, and explains how it maps onto the reference codes.

### 1.4.1 The energy equation and its derivation

Let $\xi_j$ be the bulk mass concentration of component $j$, $\mathbf v_s$ the velocity of the skeleton (every condensed component moves with it), and $\mathbf J_g = \xi_g(\mathbf v_g - \mathbf v_s)$ the gas mass flux *relative to the skeleton* — the quantity Darcy's law and the volume-fraction-gradient diffusion law return (§3.4, §3.5, §11.1.2). Every component's enthalpy is written on one datum as $h_j^0 + h_j(T)$, with $h_j^0$ its formation enthalpy at $T_{\text{ref}}$ and $h_j(T)=\int_{T_{\text{ref}}}^{T} c_{p,j}\,dT'$ its sensible enthalpy. The species balance is $\partial_t\xi_j + \nabla\!\cdot(\mathbf v_s\xi_j) = \omega_j - \nabla\!\cdot\mathbf J_j$, with $\mathbf J_j = 0$ for condensed $j$ and $\omega_j = \sum_r \nu_{r,j}\,r_r$ the net production by reaction.

**The conservative balance.** The total enthalpy has no reaction source: per bulk volume,
$\partial_t\sum_j \xi_j(h_j^0+h_j) + \nabla\!\cdot\big(\mathbf v_s\sum_j\xi_j(h_j^0+h_j)\big) + \nabla\!\cdot\sum_g (h_g^0+h_g)\mathbf J_g = -\nabla\!\cdot\mathbf q_{\text{cond}} + Q_{\text{rad}}$.
Subtracting $h_j^0$ times each species balance leaves the formation term $-\sum_j h_j^0\,\omega_j$ on the right. With the heat of reaction at the local temperature, $h_r(T) = \sum_j \nu_{r,j}\big(h_j^0 + h_j(T)\big)$ — the calibrated, DSC-style datum (§6.6.3) — that term is $-\sum_r r_r\big(h_r(T) - \sum_j\nu_{r,j}h_j(T)\big) = Q_{\text{rxn}} + \sum_j \omega_j h_j(T)$, so the sensible enthalpy obeys

$$
\partial_t\Big(\sum_j \xi_j h_j\Big) + \nabla\!\cdot\Big(\mathbf v_s \sum_j \xi_j h_j\Big) + \nabla\!\cdot\Big(\sum_g h_g \mathbf J_g\Big)
= -\nabla\!\cdot\mathbf q_{\text{cond}} + Q_{\text{rad}} + Q_{\text{rxn}} + \sum_j \omega_j\,h_j ,
\tag{1.1}
$$

with $Q_{\text{rxn}} = -\sum_r h_r r_r$ the reaction heat (§1.5.3). The last term converts the heat of reaction, quoted at the local temperature, to the sensible-enthalpy datum. Equation (1.1) is what the cell rows of the code integrate (§1.4.4).

**The temperature equation.** The storage coefficient is the **mixture volumetric heat capacity**

$$
\rho c_{p}^{\text{eff}}(T,\boldsymbol\xi)
= \sum_{j=1}^{N_C} \xi_j\,c_{p,j}(T),
\tag{1.2}
$$

summed over **every** component, gas included (discretely, the cell's heat-capacity sum $C_i = \sum_j m_{j,i}c_{p,j}(T_i)$, §1.4.4). Expanding each $\xi_j h_j$ in (1.1) with the product rule, substituting the species balances, and using $\nabla h_g = c_{p,g}\nabla T$, the composition terms cancel: $\sum_j h_j\omega_j$ appears on both sides, and $-\sum_g h_g\nabla\!\cdot\mathbf J_g$ cancels the second half of $\nabla\!\cdot(\sum_g h_g\mathbf J_g) = \sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T + \sum_g h_g\nabla\!\cdot\mathbf J_g$. Hence

$$
\Big(\sum_j \xi_j c_{p,j}\Big)\frac{D_s T}{Dt}
= -\nabla\!\cdot\mathbf q_{\text{cond}} - \sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T + Q_{\text{rad}} + Q_{\text{rxn}},
\qquad \frac{D_s}{Dt} \equiv \partial_t + \mathbf v_s\!\cdot\!\nabla ,
\tag{1.3}
$$

which is ThermaKin Eq. 17 term for term: the all-component storage sum, the convective term $c_g(\mathbf J_g\!\cdot\!\nabla T)$, and the deformation term, which is the $\mathbf v_s\!\cdot\!\nabla T$ part of $D_sT/Dt$. On a fixed mesh $D_s/Dt = \partial_t$; on a mesh that moves with the skeleton $D_sT/Dt$ is the mesh-frame rate, and neither (1.1) nor (1.3) carries an advection term of any kind. Where the equidistribution correction moves the mesh relative to the skeleton, every inventory crosses the faces at the relative velocity $-\mathbf c$ (§11.2.2, §11.6).

**Why the gas is stored.** The gas concentration $\xi_g$ is an integrated state, so the term $\xi_g c_{p,g}\,\partial_t T$ is available at no cost and dropping it would be an approximation with no benefit. Gpyro drops it in its default thermal-equilibrium form (the quasi-steady-gas approximation, restored by its `FULL_QSG` option); FDS has no in-depth gas, so the term is identically zero there. Its size is $\varphi\rho_g c_{p,g}$ against the matrix heat capacity: for a char of porosity 0.8 and bulk heat capacity $1.5\times10^5$ J·m⁻³·K⁻¹ holding pore gas at 1 atm and 800 K ($\rho_g \approx 0.44$ kg·m⁻³, $c_{p,g}\approx 1100$ J·kg⁻¹·K⁻¹), about a quarter of a percent of $\rho c_p$ — small, but exact at no cost.

### 1.4.2 The two energy states

The integrator carries one energy unknown per cell, chosen by `energy_state`. Both rows are the discrete form of (1.1) on a cell (§1.4.4); they differ in which variable is integrated:

- **Temperature form** (`energy_state = :temperature`, the default): the cell temperature $T_i$, with the row $C_i\,dT_i/dt = E_i - \sum_j h_j(T_i)\,\dot M_{j,i}$, the cell integral of (1.3). Energy is conserved to the integrator's order.
- **Inventory form** (`energy_state = :enthalpy`): the sensible-enthalpy inventory $H_i = \sum_j m_{j,i}h_j(T_i)$, with the row $dH_i/dt = E_i + \sum_j \omega_{j,i}V_i\,h_j(T_i)$, the cell integral of (1.1). The inventory then changes by exactly the integral of its rate (the inventory less the integrated rate is a linear invariant of the semi-discrete system, which every Runge–Kutta step conserves to round-off); the temperature is recovered from $H_i$ and the masses by a regularized inversion of the caloric relation (§3.3.5).

The two rows are one equation: applying the product rule to $H_i$ with the mass rows $dm_{j,i}/dt = \omega_{j,i}V_i + \dot M_{j,i}$ turns either into the other. The conserved density is $\sum_j\xi_j h_j(T)$, not $\rho c_p^{\text{eff}}\,T$, which is not an enthalpy density ($h_j(T) \ne c_{p,j}\,T$ in general); both forms conserve it, and the energy ledger closes exactly in both, to round-off in the inventory form and at the integrator's order in the temperature form (§16.5.2). The temperature form is the default because the temperature is the observable that the property laws and the surface Newton solve read directly, and because it is the form of the reference codes (§1.6); the inventory form is the choice when the energy ledger must close to round-off.

### 1.4.3 The source terms and the surface gas enthalpy

The terms of (1.1) and (1.3) are (full derivations in §7, §8, §6):

- $-\nabla\!\cdot\mathbf{q}_{\text{cond}}$, the conductive flux divergence, with Fourier's law $\mathbf{q}_{\text{cond}}=-k\,\nabla T$ and a distance-weighted harmonic-mean face conductivity (§1.4.4).
- The **gas enthalpy flux** $\nabla\!\cdot(\sum_g h_g\mathbf J_g)$: the sensible enthalpy carried by the through-flow of volatiles relative to the skeleton. It is assembled at the faces in flux form — every face carries $\sum_g J_{g,f}\,h_g(T_f)$, with the same $J_{g,f}$ the mass rows carry and $T_f$ the distance-weighted face temperature between two cells, the surface temperature on a boundary face — and its temperature-form reading is the term $-\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T$ of (1.3). It is active wherever gas moves, and it is the only place the gas flux enters the energy balance.
- $Q_{\text{rad}}$, the volumetric radiation absorption (model-dependent; §8).
- $Q_{\text{rxn}} = -\sum_r h_r r_r$, the reaction heat (sign convention in §1.5.3), with $h_r$ at the local temperature; in (1.1) it is accompanied by the formation term $\sum_j\omega_j h_j(T)$.

The surface balance itself carries no gas-enthalpy term. Under local thermal equilibrium the gas leaves through the exposed face at $T_s$: the outward mass flux of the boundary face carries its sensible enthalpy at $T_s$ in that face's energy flux, and the surface balance reduces to conducted flux equals applied flux, $q_{\text{BC}}(T_s) = q_{\text{cond}}(T_s)$. The sensible heat the escaping gas carries away is thereby counted once: its enthalpy at formation sits in $h_r(T)$, and its heating from the formation temperature to $T_s$ is the enthalpy the face fluxes carry along its path. A surface term $\sum_g \dot m_g \Delta h_g(T_s, T_\infty)$ in the balance would count that heating a second time — for an outflow of 10 g·m⁻²·s⁻¹ with $c_{p,g} = 1.5$ kJ·kg⁻¹·K⁻¹ heated by 500 K, a spurious sink of about 7.5 kW·m⁻² — and would under-predict surface temperature and MLR. The surface balance is therefore a plain convective + radiative balance — the ThermaKin Eq. 17 convention, which Gpyro's thermal-equilibrium form shares. The only boundary use of the gas flux inside the balance is the blowing model of a convective BC, which reads the outward gas flux to reduce the heat-transfer coefficient (§12.3.3).

### 1.4.4 The discrete energy update

Integrating (1.1) over a finite-volume cell $i$ gives the rows the integrator advances (`_energy_row`, `src/Assembly/residual.jl`). With $s_{if} = +1$ when the normal of face $f$ points out of cell $i$ and $-1$ otherwise, the cell's energy right-hand side is

$$
E_i = \big(Q_{\text{rxn},i} + Q_{\text{rad},i}\big)V_i - \sum_{f\in\partial i} s_{if}\,A_f\,F^{E}_{f},
\qquad
F^{E}_{f} = q_f + \sum_j F_{j,f}\,h_j(T_f),
\tag{1.4}
$$

with $F_{j,f}$ the mass flux of component $j$ through face $f$ along its normal and $F^E_f$ the face's energy flux: conduction plus the sensible enthalpy its mass fluxes carry. The energy row is then

$$
\text{temperature form:}\quad C_i\,\frac{dT_i}{dt} = E_i - \sum_j h_j(T_i)\,\dot M_{j,i},
\qquad
\text{inventory form:}\quad \frac{dH_i}{dt} = E_i + \sum_j \omega_{j,i}V_i\,h_j(T_i),
\tag{1.5}
$$

with $\dot M_{j,i} = -\sum_{f\in\partial i} s_{if}A_f F_{j,f}$ the face-flux part of the mass rate and $C_i$ the cell's heat-capacity sum $\sum_j m_{j,i}c_{p,j}$ at the trial temperature, smoothly floored at $\rho c_{p,\text{floor}}V_i(0)$ ($\rho c_{p,\text{floor}} = 10^3$ J·m⁻³·K⁻¹, `RHO_CP_FLOOR`) so that a cell emptied of all mass keeps a finite rate. The conductive face flux between the owner cell $o$ and the neighbour cell $n$ of a face is

$$
q_f = -k_f\,\frac{T_n-T_o}{d_o+d_n},
\qquad
k_f = \frac{(d_o+d_n)\,k_o k_n}{d_o k_n + d_n k_o + \eta\,(d_o+d_n)},
\tag{1.6}
$$

with $d_o$, $d_n$ the distances of the two cell centers from the face. The face conductivity $k_f$ is the **distance-weighted harmonic mean** (resistance-in-series) of the two cell conductivities, derived from flux continuity across the interface (the total thermal resistance $R = d_o/(k_o A) + d_n/(k_n A)$ inverts to it); it reduces to $2k_ok_n/(k_o+k_n)$ on a uniform mesh, and the guard $\eta = 10^{-20}$ W·m⁻¹·K⁻¹ makes it vanish continuously, with a finite derivative, as either side's conductivity does (`guarded_harmonic_mean`, `src/Properties/effective_properties.jl`). The same mean is applied to the gas-transfer coefficient and, with a smaller guard, to the Darcy mobility (§5, §9).

> **Implementation note.** Each kernel adds its face fluxes, per unit area along the face normal, and its cell sources to the contributions: conduction in `ConductionKernel` (`src/Kernels/conduction.jl`), the gas mass fluxes and the enthalpy they carry in `GasTransportKernel` (`src/Kernels/gas_transport.jl`), the boundary energy fluxes with the enthalpy of the outward mass fluxes at $T_s$ in `ThermalBoundaryKernel` (`src/Kernels/thermal_boundary.jl`), and the reaction sources in `ReactionKernel`. The row assembly `assemble_rows!` (`src/Assembly/residual.jl`) is the only writer of the derivative: it scatters each face flux as $-F_fA_f$ into the owner and $+F_fA_f$ into the neighbour and forms (1.5). The kernel order and the full residual sequence are documented in §7.6 and §15.3.

The species equation is solved in conservative form, the cell masses being the unknowns:

$$
\frac{dm_{j,i}}{dt}
= \omega_{j,i}\,V_i + \dot M_{j,i},
\qquad
\omega_{j,i}=\sum_r \nu_{r,j}\,r_{r,i},
\tag{1.7}
$$

where the face fluxes of the gas transport carry gas components only, and where the equidistribution correction moves the mesh relative to the skeleton, every component also crosses the faces at the relative velocity (§3.2.3, §11.6).

---

## 1.5 Units and coordinate conventions

All quantities are SI (see §2 for the complete table). Three conventions are restated here because they fix signs in the governing equations.

### 1.5.1 Coordinate orientation

The single spatial coordinate is the **axial / through-thickness coordinate $z$** [m]:

- $z=0$ is the **bottom / substrate**, the boundary tag `:bottom`, at face 1 and the pinned node 1;
- $z=L$ is the **top / exposed surface**, the boundary tag `:top`, at face $n+1$.

Heat enters from the **top** (the irradiated face); a charring/receding material shrinks **downward** toward the substrate. The lateral plane carries the cross-sectional area $A$. Boundary conditions are keyed by tag; there are no numeric boundary ids. Cells are numbered $i=1,\dots,n$ from bottom to top, and the state vector interleaves node and cell blocks in the order node 1, cell 1, node 2, …, cell $n$, node $n+1$ (§3.1.2, §15.1); a user who needs a field should read it through the solution's accessors (`get_temperatures`, `get_profile`, User Guide §9) rather than index the state.

### 1.5.2 Flux and divergence sign conventions

- **A face flux is the transport per unit face area along the face normal**, and the normal points from the face's owner cell to its neighbour, out of the domain on a boundary face. On a line every interior face and the `:top` face therefore have their normal along $+z$ (positive flux = transport bottom → top), and the `:bottom` face along $-z$ (positive flux = outflow through the back). Fourier's law $q_f=-k_f(T_n-T_o)/(d_o+d_n)$ gives $q_f>0$ on an interior face when $T_o>T_n$: heat flows upward, from hot to cold, opposite the gradient.
- **The discrete divergence is outflow-positive**: $(\nabla\!\cdot F)_i = (1/V_i)\sum_{f\in\partial i} s_{if}A_f F_f$ is positive when the quantity flows *out* of cell $i$. The rows are extensive (W, kg·s⁻¹): the energy right-hand side (1.4) carries $-\sum_f s_{if}A_fF^E_f$, so a net conductive outflow cools the cell, as physically required.

### 1.5.3 Heat-of-reaction sign convention

The internal (storage) convention is **$h>0$ endothermic (cools)** and **$h<0$ exothermic (heats)**, with $h$ measured per kilogram of the (single) reactant [J·kg⁻¹]. The volumetric reaction heat source carries the minus sign:

$$
Q_{\text{rxn}} = -\sum_r h_r\,r_r,
\tag{1.8}
$$

verified verbatim in `reaction_heat_source` (`src/Physics/kinetics.jl`, `Q -= heats[i]*rates[i]`). Hence an endothermic reaction ($h>0$, $r>0$) gives $Q_{\text{rxn}}<0$ (cooling), and an exothermic reaction ($h<0$) gives $Q_{\text{rxn}}>0$ (heating).

> **Overload caution (binding).** ThermaKin *publishes* the opposite sign — in its manual $h>0$ denotes an exothermic reaction — so a heat of reaction taken from ThermaKin input is entered in Pyrolysis.jl with its sign flipped. Gpyro and FDS use the same convention as Pyrolysis.jl ($h>0$ endothermic; §6.9.2), and their heats are entered as published. This is overload resolution H1 of the nomenclature (§2); the minus sign in (1.8) is what turns the stored heat into the volumetric source.

Two further nomenclature overloads that recur in this guide and are worth flagging here (full list in §2): the symbol $A$ denotes **cross-sectional area** in the physics chapters (and $A_f$ the area of face $f$), while the Arrhenius pre-exponential is always written subscripted, $A_i$ (or $A_{\alpha\beta}$); and $\chi$ always denotes **pyrolysis progress**, never a cross-section ratio (the area ratio is $A/A_0=\mathrm{law}(\bar\chi)$).

### 1.5.4 Fixed physical constants

The constants used throughout are defined in `src/Materials/constants.jl` (the air molar mass of the Chapman–Enskog estimate is `M_AIR` in `src/Properties/effective_properties.jl`) and match the nomenclature table exactly:

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | $R_g$ | 8.31446261815324 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | $\sigma$ | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference temperature | $T_{\text{ref}}$ | 298.15 | K |
| Reference pressure | $P_{\text{ref}}$ | 101325 | Pa |
| Reference molar mass (air, Chapman–Enskog estimate) | $M_{\text{ref}}$ | 0.029 | kg·mol⁻¹ |

The numerical constants and defaults — the property bracket, the floors of the state mapping, the deactivation threshold and band — are tabulated in the numerical-constants table of §2 (Nomenclature).

---

## 1.6 Relation to ThermaKin2Ds, Gpyro, and the FDS solid-phase model

Pyrolysis.jl deliberately adopts the established modeling conventions of three reference solvers so that its results are directly comparable and its inputs largely portable. This subsection maps notation and contrasts conventions; detailed equation-by-equation comparisons appear in the relevant later chapters (§3, §6, §7, §9, §15).

### 1.6.1 ThermaKin2Ds

ThermaKin2Ds is a 1-D/quasi-2-D condensed-phase pyrolysis solver widely used for inverse parameter estimation from cone/TGA data. The correspondence:

- **Energy convention.** ThermaKin's energy equation (its Eq. 17) is the same equation as (1.3): the storage term sums $c_p$ over **all** components including gas, its convective term $c_g(\mathbf{J}_g\!\cdot\!\nabla T)$ is the temperature-form reading of the gas enthalpy flux that Pyrolysis.jl assembles at the faces, and its deformation term is the skeleton-following derivative that the ALE mesh frame absorbs (§11). Because the face fluxes carry the gas enthalpy, the surface balance has no gas-enthalpy term (§1.4.3). ThermaKin's Eq. 15 applies the same deformation term to every component, gas included — the frame choice Pyrolysis.jl makes in §11.1.2.
- **Stoichiometry.** ThermaKin writes stoichiometric coefficients as $\theta_i^j$; Pyrolysis.jl uses $\nu_{i,j}$ (negative for reactants, positive for products). The symbol $\theta$ is reserved here for the volumetric strain rate, so $\theta_i^j$ appears only in ThermaKin-comparison passages and is mapped immediately to $\nu_{i,j}$ (overload A2, §2).
- **Heat of reaction.** ThermaKin publishes $h>0$ as exothermic, the opposite of Pyrolysis.jl (§1.5.3).
- **Kinetics.** ThermaKin supports bimolecular rates ($\xi_{\text{COMP1}}\xi_{\text{COMP2}}$); Pyrolysis.jl is single-reactant only (§6).
- **Boundary fluxes and mesh.** ThermaKin uses piecewise external/flame BCs with reflected incident flux, and a dynamically adjusted mesh; Pyrolysis.jl's ALE framework (§11) with cell deactivation (§14) and its surface-temperature Newton solve (§12) are the analogs.

### 1.6.2 Gpyro

Gpyro is a generalized 1-D (with limited multi-D) condensed-phase pyrolysis model built on Patankar-style control-volume discretization.

- **Conservation form.** Gpyro writes control-volume conservation of $\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, gas mass, and energy. Its bulk density $\bar\rho$ and volume fractions $X_\alpha$ map to Pyrolysis.jl's mixture density $\rho$ and volume fractions $v_j$; its pre-exponential $Z$ maps to $A_i$; its permeability $K$ maps to $\kappa$; its conversion variable $\alpha_{\text{conv}}$ is distinct from the absorption coefficient $\alpha$ (overloads A3, A4, A5, §2). Pyrolysis.jl's cell masses $m_{j,i}$ are Gpyro's $\bar\rho Y_i\,\Delta z$ times the cell's area.
- **Energy storage and gas frame.** Gpyro's thermal-equilibrium form drops the gas storage $\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default (its Eq. 65) and restores it with `FULL_QSG`; Pyrolysis.jl always keeps it (1.2). Gpyro writes its gas balance (Eq. 58) on control volumes that deform with the solid, with only the Darcy flux crossing their faces — the same skeleton-bound frame as the ALE mesh (§11.1.2).
- **Gas momentum.** Gpyro uses a Darcian mass flux $\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ with a pressure-evolution equation; Pyrolysis.jl's Darcy–Fick transport uses $v_g=-(\kappa/\mu)\nabla P$ with an ideal-gas pressure closure (§9), omitting the buoyancy term (negligible in a heated thin sample).
- **Discretization and solve.** Gpyro uses Patankar's Peclet-based scheme selection, harmonic interface interpolation, and TDMA within a fully-implicit backward-Euler linearization. Pyrolysis.jl shares the harmonic interface interpolation (§1.4.4) but writes the semi-discrete system as a mass-matrix DAE and advances it with a stiff Rosenbrock integrator and the exact colored Jacobian instead of a hand-linearized backward Euler (§13, §15).
- **Mixing rules.** Gpyro averages properties by volume fraction with power-law $T$-dependence; Pyrolysis.jl offers parallel, series, weighted, Bruggeman, and Carman–Kozeny rules per property (§5).

### 1.6.3 The FDS solid-phase pyrolysis model

FDS (Fire Dynamics Simulator) embeds a 1-D solid-phase pyrolysis model behind each wall/boundary cell of its CFD gas solver.

- **Energy and storage.** FDS solves a 1-D Fourier conduction equation with the condensed-phase heat capacity and Arrhenius source terms; it carries no in-depth gas, so the gas term of (1.2) is identically zero there. The FDS boundary condition $-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$ (convective + radiative) corresponds to Pyrolysis.jl's surface energy balance (§12).
- **Kinetics.** FDS rates can include density, oxygen, and temperature powers ($\rho^n O_2^n T^n$) and provide a recipe for recovering $A_i,E_i$ from a TGA reference temperature and peak rate; Pyrolysis.jl supports the same TGA-recovery recipe (§6, User Guide §4.9) but without the oxygen/density power factors.
- **Swelling/shrinkage.** FDS uses max/sum density-ratio rules for swelling and shrinkage; Pyrolysis.jl uses the swelling factor $\gamma_j$ and the $\chi$-driven lateral-shrinkage law (§10).
- **Radiation.** FDS uses a two-flux (Schuster–Schwarzschild) in-depth radiation treatment; Pyrolysis.jl offers surface absorption and Beer–Lambert in-depth absorption (§8), with no in-depth emission and no scattering.
- **Energy accounting.** The Pyrolysis.jl conservation ledgers and audit (§16) are the analog of FDS's energy-conservation accounting (its Appendix P).

### 1.6.4 Summary contrast

| Aspect | ThermaKin2Ds | Gpyro | FDS solid phase | **Pyrolysis.jl** |
|---|---|---|---|---|
| Dimensionality | 1-D / quasi-2-D | 1-D (limited multi-D) | 1-D per face | **1-D** |
| Gas sensible storage in $\rho c_p$ | included | dropped by default (`FULL_QSG` restores) | no in-depth gas | **included** |
| Gas enthalpy transport | $c_g\mathbf J_g\!\cdot\!\nabla T$ | conservative div | surface BC | **face fluxes $\sum_g J_g h_g(T_f)$ (ThermaKin's term in flux form)** |
| Heat-of-reaction sign | $h>0$ exo | $h>0$ endo | $h>0$ endo | **$h>0$ endo**, $Q=-\sum h r$ |
| Reactants per reaction | bimolecular | nth-order conv. | $\rho^nO_2^nT^n$ | **single reactant** |
| Gas momentum | Boyle-law Darcy | Darcy + buoyancy | — | **Darcy (no buoyancy)** |
| Time integration | Crank–Nicolson | implicit BE + TDMA | implicit | **mass-matrix DAE, `Rodas5P` with the colored Jacobian** |
| Residual differentiable by AD | no | no | no | **yes, with respect to the state** |

The last row is the distinguishing capability: the Pyrolysis.jl residual is automatically differentiable with respect to the state, which gives the integrator its exact Jacobian by colored forward-mode differentiation (§15). Parameter sensitivities, as in the reference codes, come from repeated forward runs (§17).

---

## 1.7 Document roadmap

The Technical Reference Guide proceeds from continuum theory to discretization to numerics to verification:

- **§2 Nomenclature** — the verbatim symbol table and the three sign conventions; binding overload resolutions.
- **§3 Governing Equations of the Continuum Model** — the full coupled PDE system and its discrete unknowns, the gas-vs-solid species transport, the energy balance and its two energy states, the ideal-gas pressure closure, deformation as mesh motion, and side-by-side comparison tables with the reference codes.
- **§4 Materials & Property Functions** and **§5 Effective Properties & Mixing Rules** — the constitutive layer: components, phases, property-function types and the derivatives automatic differentiation takes of them, mixture density and porosity, and the conductivity/transfer/permeability mixing rules.
- **§6 Reaction Kinetics** — the Arrhenius rate law, smooth temperature gates, the one-sided depletion limiter, stoichiometry, and the heat-of-reaction source.
- **§7 Heat Conduction, Gas Convection & Energy Assembly** — Fourier conduction, the gas enthalpy flux, the energy rows, and the residual sequence.
- **§8 Thermal Radiation** — the three radiation models and the Beer–Lambert face intensities.
- **§9 Gas-Phase Mass Transport & Pressure** and **§10 Volume Change & Lateral Shrinkage** — diffusive and Darcy fluxes, the ideal-gas pressure, the strain rate, the conversion-history pair and the cross-section rows.
- **§11 The ALE Framework**, **§12 Boundary & Initial Conditions**, **§13 Finite-Volume Discretization**, **§14 Cell Deactivation and Burnout** — the discretization, the moving mesh and its correction, and burnout.
- **§15 Temporal Integration, the DAE, and the Colored Jacobian**, **§16 Conservation Ledgers and the Audit**, **§17 Sensitivity Analysis**, **§18 Verification & Validation** — the numerical engine, the audits, and the V&V program.

Part II (the User Guide) mirrors this structure at the workflow level, beginning with installation and a quick-start cone-calorimeter run.

---

## 1.8 Limitations and open issues

The model carries several known limitations, recorded here so that users size their expectations correctly (each is revisited in the relevant chapter):

1. **Dimensionality and lateral transport.** Full solves are 1-D (meshes of more dimensions exist at the residual level only, §1.2.2); lateral shrinkage is a column-global area scaling driven by an *empirical* law $A/A_0=\mathrm{law}(\bar\chi)$, not a resolved 2-D field. The physical justification of the shrinkage law and its interaction with multi-dimensional geometry is an open question (§10).
2. **One temperature.** The matrix and the pore gas share one temperature (local thermal equilibrium); there is no two-temperature gas energy equation, so regimes with very fast, hot gas throughput are outside the model (§3.9).
3. **Single-reactant kinetics only.** Multi-reactant (bimolecular, competitive) reactions are unsupported, as are temperature-dependent activation energy, pressure dependence, and distributed-activation-energy tooling (§6).
4. **The heat-of-reaction datum.** The energy balance takes the heats of reaction at the local temperature (the DSC-style, inverse-calibrated convention), which the formation term $\sum_j\omega_j h_j(T)$ of (1.1) converts to the sensible-enthalpy datum exactly. A heat quoted at the datum $T_{\text{ref}}$ differs from the local-temperature heat by Kirchhoff's law, $h_r(T) = h_r(T_{\text{ref}}) + \sum_j\nu_{r,j}h_j(T)$, and must be converted before it is supplied (§6.6.3). The datum decides which heat is released, not whether energy is conserved: the energy ledger closes exactly in both energy forms for any composition change (§16.5.2).
5. **Transport closures are first-order / heuristic in places.** Gas advection uses first-order upwind, diffusion applies one mixture transfer coefficient to every gas species (no per-species discrimination, with only heuristic Lennard-Jones parameter estimates available for sizing it), and the pore-gas viscosity is one curve per material with no composition-dependent (Wilke) mixing (§5, §9).
6. **Burnout inside an envelope; sensitivities across events.** Cell deactivation is verified under joint refinement of the mesh and the tolerance inside the envelope stated in User Guide §10.9 — the threshold, the band, the cells across a vanishing layer, the tolerance — and outside it the end-of-layer transient need not converge with resolution (§14.2.4). Finite-difference parameter sensitivities remain valid across located deactivation events, whose times are smooth functions of the parameters, but an objective that samples an event at the data's times can be rough at small steps, so the step is chosen by a convergence check (§17.7.3).
7. **Imported calibrations must respect the surface convention.** Pyrolysis.jl carries the gas enthalpy in the face fluxes and applies no surface gas-enthalpy term (§1.4.3); parameters calibrated against a solver that applies both will predict higher $T_s$ and MLR here and should be re-validated (§18).

These limitations bound the model's domain of applicability; within that domain — 1-D, single-reactant, atmospheric-to-moderate-pressure condensed-phase pyrolysis — the formulation is rigorous, differentiable, and conservation-audited.

### 1.8.1 Known gaps and missing physics

The physics below is absent from the model. It is ordered by its expected effect on the validation cases; items already listed above are cross-referenced rather than repeated.

**Likely to bias the validation cases**

1. **Flame heat-flux feedback.** There is no ignition-triggered flame heat-flux boundary condition. After ignition a cone-calorimeter flame returns roughly 20–30 kW·m⁻² (PMMA) or 10–15 kW·m⁻² (wood) to the surface, 40–60 % on top of a 50 kW·m⁻² exposure, so parameters calibrated without it absorb that net-flux bias. The natural form is a flame-flux BC switched on by a critical-MLR or critical-$T_s$ callback. Flame feedback and the blowing reduction of the convective coefficient (the `blowing` model of `ConvectiveBC`, §12.3.3) act in opposite directions and partially cancel: enabling a blowing model alone on a cone example calibrated without either raises $T_s$ and MLR. Both belong in place before the cone examples are recalibrated.
2. **Multi-reactant kinetics** (item 3 above). A reaction with more than one reactant is rejected by problem validation; source evaluation also checks the restriction (`_apply_reaction_source!`, `src/Physics/kinetics.jl`). This blocks char oxidation, recondensation, and autocatalysis (Prout–Tompkins fits cannot be represented). The extension point is a 3-ary reactant tuple `(idx, order, mass_stoich)`; it touches the reaction type, the kinetics, the volumetric strain rate, and the reaction kernel (the colored Jacobian differentiates whatever the kernel computes).

**Relevant to near-term applications**

3. **Oxidative pyrolysis and char oxidation with O₂ ingress** (about 32.8 MJ per kg of carbon). It dominates at sub-critical fluxes, after flameout, and in smoldering, and it controls the late-time MLR tail and the residue mass of the Douglas-fir case. It needs multi-reactant kinetics (item 2) and a per-species ambient concentration, which `ConvectiveMassBC` already accepts (`ξ_inf` per species). The upwind Darcy outflow suppresses ingress at high MLR.
4. **Equilibrium (pressure-dependent) moisture evaporation and recondensation.** Drying is a pressure-independent Arrhenius reaction; no saturation pressure $P_{\text{sat}}(T)$ appears anywhere. Vapor recondensation ahead of the drying front (which raises the local moisture content by 30–50 %) is absent, and at the 1.5–2 atm pore pressure of pressure-treated wood the boiling point shifts to 112–121 °C, which Arrhenius drying cannot represent. It needs item 2 and a smooth $r \propto (P_v - P_{\text{sat}}(T))$ closure; the pressure is already evaluated in the state mapping under Darcy–Fick transport.
5. **Backing layer and contact resistance at the back face.** There is no contact-resistance or substrate model. An adiabatic back face over-predicts late-time back-half temperatures (the missing loss is about 0.7–3 kW·m⁻² for 25 mm of wood over 10–20 min); a lumped Robin condition on the bottom face (radiative + convective, as in example `06_cone_calorimeter.jl`) is the available approximation. A contact conductance plus a semi-infinite-substrate $h(t)$ would close the gap.
6. **Kinetics expressiveness.** Normalized-conversion ($\alpha$) tracking, a $T^n$ pre-factor, and DAEM-style stacking are absent. Any real order $n \ge 0$ per reactant ($\xi^n$) is supported and equals $\xi_0^n(1-\alpha)^n$ up to a rescaling of $A_i$, so the gap is $\alpha$-normalization, not the conversion order.
7. **Transport-closure accuracy** (item 5 above): no composition-dependent (Wilke) viscosity mixing, air-default Lennard-Jones parameters for the Chapman–Enskog estimate, and first-order upwind Darcy advection with no TVD limiter; and under Darcy–Fick transport the diffusive flux keeps the isobaric driving potential $-(\lambda/T)\,\partial(\xi T)/\partial z$, dropping the term proportional to $\nabla(\varphi P)$ that a resolved pressure field and sharp char fronts carry (`gas_flux`, `src/Physics/gas_transport.jl`).

**Niche or future**

- Spectral (two-band) radiation and in-depth emission from glowing char; a model of that kind would add further algebraic intensity rows to the state, as the Beer–Lambert face intensities do (§8).
- Substrate reflection; scattering ($\sigma_s$ is not modeled); view-factor edge fall-off.
- Knudsen/Klinkenberg slip, $\kappa_{\text{eff}} = \kappa(1 + b/P)$: $b \approx 28$ kPa for 1 µm pores but about 280 kPa for 0.1 µm, up to about 4× $\kappa$ at 1 atm in sub-µm char.
- Stefan–Maxwell multicomponent diffusion (first relevant with O₂ ingress against the pyrolyzate efflux) and a tortuosity factor (the porosity $\varphi$ over the tortuosity, a factor of 2–10) on the diffusivity and the gas-transfer coefficient.
- Liquid free-water capillarity (relevant only above the fiber-saturation point, about 30 % moisture content) and bound-water diffusion ($\sqrt{Dt} \approx 0.25$ mm over a cone test, negligible at fire heating rates).
- Mechanics (char cracking, delamination, spalling), anisotropic conductivity (about 2× along the wood grain), melting and dripping, and 2-D/3-D edge burning (about a 10–20 % effect on late-time MLR in cone tests; item 1 above).
