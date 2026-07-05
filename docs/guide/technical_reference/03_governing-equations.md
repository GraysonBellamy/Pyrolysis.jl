# 3. Governing Equations of the Continuum Model

This chapter states the full coupled system of partial differential equations that
Pyrolysis.jl solves *before* it is discretized. We introduce the primary state
variables (cell temperature $T$ and component mass concentrations $\xi_j$) and the
auxiliary states carried only in certain simulation modes (node positions $z_i$ for
the Arbitrary Lagrangian–Eulerian mesh, and the pyrolysis-progress scalar $\chi_i$
for variable cross-section), then derive each conservation law: condensed-phase
mass, per-species transport, the non-conservative ("Formulation A") energy equation,
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
  equilibrium at each point. The justification (gas residence time $\ll$ matrix
  thermal diffusion time) is the same quasi-steady-gas argument that lets us drop gas
  thermal storage in §3.3.
- $\xi_j$ [kg·m⁻³] is the **mass concentration** of component $j$ — mass of $j$ per
  unit *bulk* volume of the column. This is the canonical species state variable, *not*
  the mass fraction $Y_j = \xi_j/\rho$. Choosing $\xi_j$ avoids ever forming $\rho$ at
  intermediate points and degrades gracefully as a component vanishes.

In the discrete model the state is cell-centered: $T_i$ and $\xi_{j,i}$ for cells
$i=1,\dots,n$. Two **auxiliary state blocks** are appended only in the modes that need
them (the bundled `SimulationMode{Mesh,Rad,Trans,Chi}` trait selects which; see §15):

- **Node positions** $z_i$ [m], $i=1,\dots,n_{\mathrm{nodes}}$, present only in ALE
  mode (`Mesh = ALE`). The mesh becomes part of the coupled DAE; its evolution is
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
&\text{Energy (Formulation A):}\quad
\rho c_p^{\mathrm{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\cdot\mathbf q_{\mathrm{cond}} + S_{\mathrm{conv}} + S_{\mathrm{gen}}
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
density — and it omits in-pore gas thermal storage; both choices are justified in
§3.3. Each term and its sign are derived in the sections that follow.

---

## 3.2 Condensed-phase mass and total density

### 3.2.1 Mixture (total) density

The bulk **total density** of the mixture is the sum of all component mass
concentrations,

$$
\boxed{\;\rho = \sum_{j=1}^{N_C} \xi_j. \;}
$$

This is exact by the definition of $\xi_j$ as mass-per-bulk-volume and requires no
property evaluation. It is the quantity returned by `total_density` and is used, for
example, to form mass fractions $Y_j = \xi_j/\rho$ wherever those are needed
(`mass_fractions`, with a $\rho<10^{-20}$ guard returning zero fractions).

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

> **Implementation notes.** `total_density(ξ) = sum(ξ)` and `mixture_density` are in
> `src/Properties/mixture_rules.jl` (lines ~170 and ~674). The $\rho_{\mathrm{mix}}$
> swelling form is `_specific_vol_with_swelling_direct` (same file).

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
invariant is checked by `validate_species_source_conservation` (tolerance
$10^{-12}$) and `check_material_mass_balance` (§6, §16).

### 3.2.3 Phase-resolved mass: solids/liquids are immobile, only gases transport

The per-species conservation law is the foundation of the model. In integral
(finite-volume) form over a cell $i$,

$$
\frac{d}{dt}\!\int_{V_i}\!\xi_j\,dV
= -\oint_{\partial V_i}\!\mathbf J_j\cdot\mathbf n\,dA
+ \int_{V_i}\! \left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}\!dV ,
$$

where $\mathbf J_j$ is the **mass flux** of component $j$. The decisive modeling choice
is *which components transport*. In Pyrolysis.jl the matrix is rigid relative to
itself: solids and liquids do **not** flux between cells — they are consumed/produced
in place by reactions and (in ALE) carried by the mesh, never by molecular transport.
Only gas components diffuse and advect through the pore network. Hence the
**gas-vs-solid dispatch**:

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
> `apply_species_residual_unrolled!` in `src/Materials/materials.jl` (≈ line 890).
> Species sources $S_{j,i}=\sum_r\nu_{r,j}r_r$ come from `species_source_terms!` in
> `src/Physics/kinetics.jl`; the per-reaction application enforces the single-reactant
> mass-balance convention in `_apply_reaction_source!` (same file). The species
> divergence $(J_R-J_L)A/V_i$ is `divergence_species_1d!` in
> `src/Discretization/finite_volume.jl`.

---

## 3.3 The energy equation: Formulation A (non-conservative)

### 3.3.1 Statement

Pyrolysis.jl solves the **non-conservative** energy equation written directly for the
temperature $T$, with the time-storage term carried by the **matrix-only** volumetric
heat capacity:

$$
\boxed{\;
\rho c_p^{\mathrm{eff}}(T,\xi)\,\frac{\partial T}{\partial t}
= -\,\nabla\cdot \mathbf q_{\mathrm{cond}}
+ S_{\mathrm{conv}}
+ S_{\mathrm{gen}}
+ Q_{\mathrm{rad}}
+ Q_{\mathrm{rxn}} .
\;}
$$

The terms, with sign conventions, are:

| Term | Meaning | Continuous form | Sign / units |
|---|---|---|---|
| $\rho c_p^{\mathrm{eff}}\,\partial_t T$ | matrix sensible-heat storage | $\bigl(\sum_{j\in S/L}\xi_j c_{p,j}(T)\bigr)\partial_t T$ | J·m⁻³·s⁻¹ |
| $-\nabla\cdot\mathbf q_{\mathrm{cond}}$ | Fourier conduction | $\mathbf q_{\mathrm{cond}} = -k_{\mathrm{eff}}\nabla T$ | net outflow cools |
| $S_{\mathrm{conv}}$ | gas advection of sensible enthalpy | $-\sum_g c_{p,g}(T)\,\mathbf J_g\!\cdot\!\nabla T$ | W·m⁻³ |
| $S_{\mathrm{gen}}$ | optional gas-generation enthalpy sink | $-\sum_g \Delta h_g(T,T_{\mathrm{ref}})\,\nabla\!\cdot\!\mathbf J_g$ | W·m⁻³ |
| $Q_{\mathrm{rad}}$ | volumetric radiation absorption | model-dependent | W·m⁻³ |
| $Q_{\mathrm{rxn}}$ | reaction heat | $-\sum_r h_r r_r$ | W·m⁻³ |

This is the "Formulation A" referred to throughout the codebase, selected by the
default `energy_form = :standard`. The discrete cell update is

$$
\frac{dT_i}{dt}
= \frac{1}{\rho c_{p,i}^{\mathrm{eff}}}
\Bigl(
-(\nabla\cdot\mathbf q)_i + S_{\mathrm{conv},i} + S_{\mathrm{gen},i}
+ Q_{\mathrm{rad},i} + Q_{\mathrm{rxn},i}
\Bigr),
$$

which is precisely the final assembly loop of `compute_rhs!` (the divergence enters
with a minus sign, the four volumetric sources add, the whole is divided by the
matrix $\rho c_p$). The detailed discretization of each term — the distance-weighted
harmonic-mean face conductivity, the second-order cell-centered temperature gradient in
$S_{\mathrm{conv}}$, the eight-step assembly order — belongs to §7; here we justify the
*continuous* form and especially the two non-textbook modeling choices.

> **Implementation notes.** The final RHS assembly
> $dT_i/dt = (-\,\mathrm{div\_q} + S_{\mathrm{conv}} + S_{\mathrm{gen}} + Q_{\mathrm{rad}} + Q_{\mathrm{rxn}})/\rho c_p$
> is in `compute_rhs!`, `src/Discretization/finite_volume.jl` (Step 8, ≈ line 387).
> The matrix heat capacity is `effective_heat_capacity` and the conductive flux is
> `conductive_flux` / `compute_conductive_fluxes!`.

### 3.3.2 The matrix-only heat capacity and the quasi-steady-gas approximation

The storage coefficient is the **matrix** (solid + liquid) volumetric heat capacity,
with gas components excluded:

$$
\rho c_p^{\mathrm{eff}}(T,\xi)
= \sum_{j\in\{\text{solid},\text{liquid}\}} \xi_j\, c_{p,j}(T).
$$

The exclusion of gas storage is deliberate and is the defining feature of Formulation
A. To see what is dropped, start from a one-temperature energy balance for
a porous medium in which both the matrix and the in-pore gas store sensible heat
(the form below is itself **not** assumption-free — see the caveat after the
equation):

$$
\underbrace{\sum_{j\in S/L}\!\xi_j c_{p,j}\,\frac{\partial T}{\partial t}}_{\text{matrix storage}}
+ \underbrace{\sum_{g}\!\xi_g c_{p,g}\,\frac{\partial T}{\partial t}}_{\text{gas storage}}
= -\nabla\cdot\mathbf q_{\mathrm{cond}}
- \nabla\cdot\Bigl(\textstyle\sum_g c_{p,g}(T-T_{\mathrm{ref}})\mathbf J_g\Bigr)
+ Q_{\mathrm{rad}} + Q_{\mathrm{rxn}} .
$$

Three assumptions are already embedded in this "starting" balance: (i) the
reaction heats inside $Q_{\mathrm{rxn}}$ are referenced to the $T_{\mathrm{ref}}$
datum (with DSC-style heats quoted at the reaction temperature, the gas
products' sensible enthalpy is *inside* $h_r$ and the gas-enthalpy divergence
term must not be double-counted — §6.6.3); (ii) the gas enthalpy flux is
written with the *mixture* temperature (local thermal equilibrium and
quasi-steady gas are partly pre-applied in this one-temperature form); and
(iii) the condensed-product sensible-generation term
$-\sum_{j\in S/L}\Delta h_j(T,T_{\mathrm{ref}})\,\dot\omega_j$ is omitted (its
inclusion would be required for a fully conservative-equivalent assembly —
see the caveat in `compute_gas_generation_sink!`).

Formulation A then makes two coupled approximations:

1. **Drop the gas storage term** $\sum_g \xi_g c_{p,g}\,\partial_t T$. This is the
   *quasi-steady-gas approximation*: the gas residence (advection) time through the
   pores is far shorter than the matrix thermal-diffusion time, so the gas does not
   appreciably store transient heat. At atmospheric pressure the in-pore gas mass
   concentration is tiny relative to the matrix, and the dropped term is roughly
   $0.2\%$ of the matrix storage; the error scales with pressure (it grows in
   high-pressure, low-porosity regimes). The cost is not hidden — it is quantified at
   runtime by the diagnostics "physical-form residual," which compares the matrix
   energy ledger against a total (matrix + gas) ledger (§16).
2. **Expand the gas-enthalpy divergence by the product rule** to move the convected
   sensible enthalpy onto the right-hand side. Writing $\Delta h_g = c_{p,g}(T-T_{\mathrm{ref}})$
   for the moment (constant $c_p$),
   $$
   \nabla\cdot\!\Bigl(\sum_g c_{p,g}(T-T_{\mathrm{ref}})\mathbf J_g\Bigr)
   = \underbrace{\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T}_{-\,S_{\mathrm{conv}}}
   + \underbrace{(T-T_{\mathrm{ref}})\sum_g c_{p,g}\,\nabla\!\cdot\!\mathbf J_g}_{-\,S_{\mathrm{gen}}} .
   $$
   The first piece is the **advection source** $S_{\mathrm{conv}}$; the second is the
   **generation sink** $S_{\mathrm{gen}}$. Formulation A keeps $S_{\mathrm{conv}}$
   always and treats $S_{\mathrm{gen}}$ as optional (§3.3.4).

The net result is a temperature equation whose storage coefficient is matrix-only and
whose gas-sensible-enthalpy transport is carried entirely by the volumetric source
$S_{\mathrm{conv}}$ (plus, optionally, $S_{\mathrm{gen}}$). This is the convention used
by FDS and Gpyro; it differs from ThermaKin, whose Eq. 17 keeps the full
$\sum_j\xi_j c_{p,j}$ (all phases) in the storage term (§3.7).

> **Implementation notes.** `effective_heat_capacity` (matrix only) and
> `total_heat_capacity` (matrix + gas, diagnostics only) are in
> `src/Properties/effective_properties.jl` (≈ lines 322–360). The generated summation
> `_sum_heat_capacity` explicitly guards `if !is_gas(comps[i])` so that gases never
> enter the storage coefficient.

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
$I(z) = I_{\mathrm{surface}}\exp(-\int\alpha\,dz')$; for `P1_QUASI_STEADY` it is the net
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
> `reaction_heat_source` in `src/Physics/kinetics.jl` (≈ line 257), with the negative
> sign documented in-source ("Note the negative sign!"). Reaction rates and heats are
> evaluated together by `evaluate_reactions!` (same file).

### 3.3.4 The advection source $S_{\mathrm{conv}}$ and the optional generation sink $S_{\mathrm{gen}}$

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
This term is **always on** in the `:standard` form; it is the
only carrier of gas sensible enthalpy in the temperature equation under Formulation A.

**$S_{\mathrm{gen}}$ (optional).** The generation sink is the second product-rule piece,
written with the *exact* sensible-enthalpy integral rather than a single-point product:

$$
S_{\mathrm{gen}} = -\sum_g \Delta h_g(T,T_{\mathrm{ref}})\,\nabla\cdot\mathbf J_g,
\qquad
\Delta h_g(T,T_{\mathrm{ref}}) = \int_{T_{\mathrm{ref}}}^{T} c_{p,g}(\tau)\,d\tau
\approx c_{p,g}\!\Bigl(\tfrac{T+T_{\mathrm{ref}}}{2}\Bigr)(T-T_{\mathrm{ref}}),
$$

the integral being evaluated by the midpoint rule in `delta_enthalpy` (§4) — the very
same quadrature and datum $T_{\mathrm{ref}}$ used by the optional surface-transpiration
boundary term. Using the true integral (not $(T-T_{\mathrm{ref}})c_{p,g}(T)$) is what lets
the discrete energy ledger telescope *exactly* with the surface term. $S_{\mathrm{gen}}$
is included only when `energy_form = :with_generation_sink`; otherwise it is zeroed.
Enabling it is a *partial* move toward fully-conservative form: full conservation would
additionally restore gases to $\rho c_p^{\mathrm{eff}}$, which Formulation A does not do.

**Mutual exclusivity with surface transpiration.** With constant $c_p$ the volumetric
integral of the two product-rule pieces equals the boundary outflow term,

$$
\int_0^L (S_{\mathrm{conv}} + S_{\mathrm{gen}})\,dz
= -\sum_g \dot m_g\,\Delta h_g(T_s, T_{\mathrm{ref}}),
$$

i.e. the conservative-form gas-enthalpy divergence is recovered at the boundary. The
default convention carries this transport *volumetrically* via $S_{\mathrm{conv}}$ and a
plain convective+radiative surface BC, with **no** surface-transpiration enthalpy term.
If one also enabled the surface-transpiration BC
$\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$, it would apply the *entire* outflow on top
of $\int S_{\mathrm{conv}}$, double-counting roughly $5$–$10$ kW·m⁻² as a spurious sink.
For this reason `use_transpiration_bc = true` is **rejected at problem construction**;
the supported production path is $S_{\mathrm{conv}}$-only (optionally plus
$S_{\mathrm{gen}}$, which then telescopes with transpiration only if transpiration were
permitted — it is not, in production). The surface energy balance and the
transpiration term are detailed in §7 and §12.

> **Implementation notes.** $S_{\mathrm{conv}}$ is `compute_gas_convective_source!`
> and $S_{\mathrm{gen}}$ is `compute_gas_generation_sink!`, both in
> `src/Physics/convection.jl`. The construction-time guard rejecting
> `use_transpiration_bc` is in `to_problem_def`, `src/Problem/adapter.jl`.
> The midpoint-rule `delta_enthalpy` is in `src/Materials/components.jl`
> (≈ line 492).

### 3.3.5 Why temperature, not enthalpy: trade-offs of the non-conservative form

Solving for $T$ (rather than an enthalpy density $\rho h$) keeps the primary unknown
directly observable, makes property evaluation and the surface Newton solve natural,
and is the FDS/Gpyro/ThermaKin tradition for solid-phase pyrolysis. The price is that
the *discrete* energy balance is only approximately conservative: it does not evolve the
composition-change enthalpy term $\sum_j h_j(T)\,\partial_t\xi_j$ as part of a conserved
enthalpy density (§3.3.2) — its reaction and gas-transport parts are re-expressed as
$Q_{\text{rxn}}$ and $S_{\text{conv}}$/$S_{\text{gen}}$, the condensed-product
sensible-generation part is omitted — and (in ALE) the volume-change contribution to
storage is likewise outside the ledger. These
discrete-closure gaps are bounded, are made visible by the conservation diagnostics
(§16), and in practice are dominated by the quasi-steady-gas storage gap quantified in
§3.3.2.

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
species equations reproduce total gas continuity with no spurious source. (The
model's actual quasi-steady approximation is in the *energy* treatment — the
neglected gas sensible-heat storage, §1.4.1.) The closure is built
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
falls back to its bulk density, recovering the legacy "no pores in the pure phase"
behaviour. (A second, volume-fraction porosity $\varphi=\sum_j v_j^{\text{gas}}$ also
exists historically, but the intrinsic-density form above is the one used by production
property updates and by the pressure closure.) The porosity derivatives the Jacobian
needs ($\partial\varphi/\partial\xi_j = -1/\rho_{i,j}$ for solid/liquid, $0$ for gas, in
the clamp-inactive interior) are obtained by forward-mode AD through `porosity`
itself (§5.3.3, §15).

Pressure matters only when Darcy flow is active (`DarcyFick` mode). In the diffusion-only
mode the pressure field is never assembled; the gas flux is purely Fickian and the
pressure closure is inert. The pressure derivatives $\partial P/\partial T$ and
$\partial P/\partial\xi_k$ (gas direct term plus the porosity-coupling term) that the
Jacobian needs are obtained by forward-mode differentiation of the closure itself —
see §9.6.

> **Implementation notes.** $P = N/\max(\varphi,10^{-12})$ is `ideal_gas_pressure` /
> `compute_pressure!` in `src/Physics/pressure.jl` (porosity floor `_PHI_FLOOR = 1e-12`,
> ≈ line 15; gas-only partial-pressure sum `_gas_partial_pressure_sum`, ≈ line 22). The
> canonical porosity is `porosity` in `src/Properties/mixture_rules.jl`; its
> derivatives are taken by AD through the same function (§5.3.3).

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
$\mu$ from Sutherland's law at the cell temperature) and combines the two sides
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
> Sutherland viscosity is `effective_gas_viscosity` (§5).

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
> (≈ line 103, with `if γ > 0` guards). Absolute rates are `compute_cell_volume_rates!`.
> (Some KB glosses write $\theta=-(\sum_j\dot\xi_j)/\rho_{\mathrm{total}}$; the code
> computes the swelling-weighted per-component form above — that is what is
> implemented.)

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
precision. The mesh kinematics, the ALE advection corrections
($+\,w\cdot\nabla T$ for energy, $+\,w\cdot\nabla\xi$ for gas — the lab-frame gas flux
is already in the Eulerian RHS, so the same chain-rule correction applies — zero for
solids since the solid velocity equals the mesh velocity), the solid **dilation correction**
$\Delta(d\xi_j/dt) = -\xi_j\theta_i$ that keeps solids Lagrangian, and the GCL machinery
are derived in full in §11 (The ALE Framework) and §10 (Volume Change & Lateral
Shrinkage). This chapter records only that the deformation physics is *the same* as
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
> `src/Residual/residual.jl` (≈ lines 355–415, with `w[1]=0`). The χ dynamics
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

| Aspect | Pyrolysis.jl (Formulation A) | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| Primary unknown | temperature $T$ | $T$ | $T$ | $T$ |
| Storage coefficient | matrix-only $\sum_{S/L}\xi_j c_{p,j}$ | **all phases** $\sum_j\xi_j c_{p,j}$ (Eq. 17) | matrix (quasi-steady gas) | matrix (quasi-steady gas) |
| Gas sensible transport | $S_{\mathrm{conv}}=-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$ | convective term $c_g(\mathbf J_g\!\cdot\!\nabla T)$ | gas energy eq. / advection | boundary outflow |
| Conduction | Fourier, harmonic-mean face $k$ | Fourier | Fourier, Patankar harmonic interface | 1D Fourier |
| Reaction heat | $-\sum_r h_r r_r$ ($h>0$ endo) | $\sum h_i r_i$ ($h>0$ exo) | $\sum(\Delta H)\dot\omega'''$ | $\sum h_r\dot\omega'''$ |

The cardinal difference is the **storage coefficient**: ThermaKin Eq. 17 keeps gas
sensible heat in the storage term, whereas Pyrolysis.jl (matching FDS and Gpyro) drops
it under the quasi-steady-gas approximation and carries gas enthalpy by the volumetric
advection source $S_{\mathrm{conv}}$ (§3.3.2). Pyrolysis.jl's $S_{\mathrm{conv}}$ is
exactly ThermaKin's convective term $c_g(\mathbf J_g\cdot\nabla T)$ in form; the optional
$S_{\mathrm{gen}}$ recovers the second product-rule piece that a fully-conservative
formulation would carry. The one-dimensional Cartesian energy equation is identical in
all four codes up to these conventions; ThermaKin additionally supports a cylindrical
(radial) Laplacian $\frac1r\partial_r(rk\partial_r T)$, which Pyrolysis.jl approximates
only through the scalar lateral-shrinkage area $A(t)$ rather than resolving a radial
coordinate.

### 3.7.3 Control-volume conservation forms

Gpyro writes its discretization as control-volume balances of $\bar\rho\,\Delta z$ (bulk
mass per area), $\bar\rho Y_i\,\Delta z$ (per-species mass), gas mass, and energy, with
$\Delta z$ itself a state (variable control-volume thickness). Pyrolysis.jl's ALE mesh
motion (§3.6.3) is the moving-mesh analog: instead of rescaling $\Delta z$ inside fixed
control volumes, the nodes move so that $V_i = A|z_{i+1}-z_i|$ tracks the strain, with the
discrete GCL guaranteeing the volume bookkeeping is exact (§11). FDS adjusts solid-layer
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
| Matrix storage $\rho c_p^{\mathrm{eff}}$ | property + mixing | $\sum_{S/L}\xi_j c_{p,j}(T)$ | §4, §5 |
| Conductive flux $\mathbf q_{\mathrm{cond}}$ | Fourier + mixing | $-k_{\mathrm{eff}}(T,\xi)\nabla T$, harmonic-mean face $k$ | §5, §7 |
| Reaction rate $r_r$ | Arrhenius | $A_r e^{-E_r/R_gT}\prod_j\xi_j^{n_{r,j}}$, tanh gates/limiter | §6 |
| Species source $S_{j}$ | stoichiometry | $\sum_r\nu_{r,j}r_r$ | §6 |
| Reaction heat $Q_{\mathrm{rxn}}$ | enthalpy of reaction | $-\sum_r h_r(T,\xi)\,r_r$ | §6 |
| Radiation source $Q_{\mathrm{rad}}$ | radiation model | surface / Beer–Lambert / P1 | §8 |
| Gas advection $S_{\mathrm{conv}}$ | gas flux + $c_{p,g}$ | $-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$ | §7, §9 |
| Gas generation sink $S_{\mathrm{gen}}$ | enthalpy integral | $-\sum_g\Delta h_g\,\nabla\!\cdot\!\mathbf J_g$ (optional) | §7, §9 |
| Gas diffusive flux $\mathbf J_j^{\mathrm{diff}}$ | Fick (volume-fraction gradient) | $-\lambda_{\mathrm{face}}(\xi_R T_R-\xi_L T_L)/(T_{\mathrm{face}}d_{LR})$ | §9 |
| Gas advective flux $\mathbf J_j^{\mathrm{adv}}$ | upwind intrinsic density $\times$ Darcy | $(\xi_j^{\mathrm{upwind}}/\varphi_{\mathrm{face}})\,v_g$ | §9 |
| Pressure $P$ | ideal gas in pores | $N(\xi,T)/\varphi$ | §9 |
| Porosity $\varphi$ | intrinsic density | $\operatorname{clamp}(1-\sum_{S/L}\xi_j/\rho_{i,j},0,1)$ | §5 |
| Gas velocity $v_g$ | Darcy | $-(\kappa/\mu)\nabla P$, distance-weighted harmonic mobility $\kappa/\mu$, Sutherland $\mu$ | §9 |
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

- **Quasi-steady gas thermal storage.** Formulation A drops $\sum_g\xi_g c_{p,g}\,\partial_t T$
  from the storage term (§3.3.2). The error is $\sim 0.2\%$ of matrix storage at 1 atm
  but **scales with pressure**; no automatic switch to include gas storage in
  high-pressure regimes exists. The cost is exposed by the physical-form residual (§16)
  but is never folded back into the solve.
- **One temperature.** Matrix and in-pore gas are assumed in local thermal equilibrium;
  there is no two-temperature gas energy equation. This is consistent with the
  quasi-steady-gas argument but breaks for very fast, hot gas throughput.
- **Algebraic (not evolved) pressure.** Pressure is closed instantaneously by
  $P=N/\varphi$ (§3.4.2) with a $\varphi$ floor of $10^{-12}$; there is no transient
  gas-mass-balance pressure PDE as in Gpyro. Acoustic and fast compressibility transients
  are not represented.
- **First-order upwind advection** for the Darcy–Fick advective flux (§3.5.1); no TVD or
  MUSCL limiter on the gas advection (the ALE *mesh* advection does use a TVD limiter,
  §11). Over many steps, upwinding can push concentrations slightly non-physical if not
  damped by reaction/heat coupling.
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
- **Asymmetric ambients with a permeable back face.** The $S_{\mathrm{gen}}$ form uses a
  single scalar $T_{\mathrm{ref}}$; if the two boundaries have different ambient
  temperatures and the back face is permeable, the energy ledger fails to close by
  exactly $\sum_g\dot m_{g,\mathrm{bottom}}c_{p,g}(T_{\infty,\mathrm{top}}-T_{\infty,\mathrm{bottom}})$.
  Harmless for the usual impermeable back face (§3.3.4, §7).

These limitations are inherited by every downstream chapter; where a later chapter
relaxes or quantifies one of them (e.g. §16 for the energy/mass ledgers, §18 for
regime-validity checks), the cross-reference is given there.
