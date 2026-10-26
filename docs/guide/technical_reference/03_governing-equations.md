# 3. Governing Equations of the Continuum Model

This chapter states the full coupled system of partial differential equations that
Pyrolysis.jl solves *before* it is discretized, and the discrete unknowns that carry
it. We introduce the continuum state (temperature $T$ and component mass
concentrations $\xi_j$), the integrated unknowns of every cell (the component masses
$m_{j,i}$ and one energy state), and the further unknowns carried only when their
physics is on (node positions and velocities on a moving mesh, face intensities under
Beer–Lambert radiation, the conversion-history pair and the cross-section rows under a
lateral-shrinkage law), then derive each conservation law: condensed-phase mass,
per-species transport, the energy balance and the temperature equation it implies, the
ideal-gas pressure closure, and Darcy momentum. The chapter is deliberately a
*statement and justification* of the continuum model; the discrete operators that
realize it are derived in §7 (Heat Conduction, Gas Convection & Energy Assembly), §9
(Gas-Phase Mass Transport & Pressure), §11 (The ALE Framework), and §13 (Finite-Volume
Discretization). Every equation here is verified against the implementing source and,
where the textbook form and the code differ, the chapter documents what the **code**
actually does.

All symbols, units, and sign conventions are those fixed in §2 (Nomenclature). The
three binding sign conventions are repeated here because the entire chapter rests on
them:

1. **Flux positivity.** A face flux is the transport of its quantity per unit face
   area along the face normal, and the normal points from the face's owner cell to its
   neighbour, out of the domain on a boundary face. On a line every interior face and
   the top face ($z=L$, tag `:top`) have their normal along $+z$, so a positive flux
   moves the quantity from the bottom/substrate toward the exposed surface; the bottom
   face ($z=0$, tag `:bottom`) has its normal along $-z$, so a positive flux there is
   an outflow through the back.
2. **Divergence positivity.** The discrete divergence
   $(\nabla\cdot\boldsymbol F)_i = (1/V_i)\sum_{f\in\partial i} s_{if}A_f F_f$, with
   $s_{if}=+1$ when the normal of face $f$ points out of cell $i$, is positive when the
   quantity flows *out* of cell $i$ on net. In the energy and species balances the
   divergence enters with a leading minus sign, so net outflow lowers the stored
   quantity.
3. **Heat-of-reaction sign.** Internally $h_r > 0$ is **endothermic** (the reaction
   cools the matrix) and $h_r < 0$ is **exothermic**. The volumetric reaction heat
   source therefore carries an explicit minus sign,
   $Q_{\mathrm{rxn}} = -\sum_r h_r r_r$, so an endothermic reaction
   ($h_r>0$, $r_r>0$) gives $Q_{\mathrm{rxn}}<0$. (ThermaKin publishes the opposite
   convention, $h>0=$ exothermic; Gpyro and FDS share ours. The mapping is given in
   §3.7 and in overload note H1 of §2.)

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
phase-specific behaviour is decided when the problem is set up, not in the inner loops
(§4). Solids and liquids form the load-bearing *matrix*; gases are decomposition
products that percolate through the matrix pores and ultimately escape. Reactions
convert a single reactant component into one or more product components (§6).

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
  unit *bulk* volume of the column. This is the canonical constitutive species
  variable, *not* the mass fraction $Y_j = \xi_j/\rho$. Choosing $\xi_j$ avoids ever
  forming $\rho$ at intermediate points and degrades gracefully as a component
  vanishes.

In the discrete model the state is cell-centered, and the integrated unknowns of
cell $i=1,\dots,n$ are

- the **component masses** $m_{j,i}=\int_{V_i}\xi_j\,dV$ [kg]; every kernel reads the
  concentrations $\xi_{j,i}=m_{j,i}/V_i$, and
- one **energy state**: the temperature $T_i$ [K] (`energy_state = :temperature`, the
  default) or the sensible-enthalpy inventory $H_i = \sum_j m_{j,i}h_j(T_i)$ [J]
  (`energy_state = :enthalpy`), $h_j(T)=\int_{T_{\mathrm{ref}}}^{T}c_{p,j}\,dT'$
  (§3.3.5).

Further unknowns appear only in the configurations that need them (`setup` builds the
layout from the problem and the solve keywords; §15.1):

- **Node positions** $z_n$ [m] and **node velocities** $w_n$ [m·s⁻¹],
  $n=1,\dots,n+1$, on a moving mesh (`use_ale = true`). The positions are differential
  unknowns, $dz_n/dt = w_n$; the velocities are algebraic, each held by the volume
  balance of one cell (§3.6.3, §11.3).
- **Face intensities** $I_f$ [W·m⁻²], one per face, under `BEER_LAMBERT`: algebraic
  unknowns that satisfy the attenuation of the beam across each cell (§3.3.3, §8).
- The **conversion-history pair** $(G_i, m_{\mathrm{dry},i})$ [kg] of every cell and the
  global **cross-section area** $A$ [m²] with its **lateral strain rate** $\theta_A$
  [s⁻¹] under a lateral-shrinkage law: the pair differential, $A$ and $\theta_A$
  algebraic (§3.6.4, §10).

The state vector interleaves node and cell blocks, so that every row's unknowns lie
near it and the Jacobian is banded:

$$
\boldsymbol u =
\bigl[\,\underbrace{z_1, w_1, I_1}_{\text{node }1},\;
\underbrace{m_{1,1},\dots,m_{N_C,1},\,T_1,\,G_1,\,m_{\mathrm{dry},1}}_{\text{cell }1},\;
\underbrace{z_2, w_2, I_2}_{\text{node }2},\;\dots,\;
\underbrace{z_{n+1}, w_{n+1}, I_{n+1}}_{\text{node }n+1},\;
A,\;\theta_A \,\bigr],
$$

with $T_i$ the energy state (the inventory $H_i$ under `energy_state = :enthalpy`)
and each entry present only in its configuration. On a line its length is

$$
N_{\mathrm{state}} = (n+1)\,(2\delta_{\mathrm{mov}} + \delta_{\mathrm{BL}})
+ n\,(N_C + 1 + 2\delta_{\chi}) + 2\delta_{\chi},
$$

where $\delta_{\mathrm{mov}}$, $\delta_{\mathrm{BL}}$, $\delta_\chi$ equal one on a
moving mesh, under Beer–Lambert radiation and under a lateral-shrinkage law, and zero
otherwise. For three components on four cells the state has 16 rows on a fixed mesh;
36 on a moving mesh under a lateral-shrinkage law, of which 7 are algebraic (the five
node velocities, $A$ and $\theta_A$); and 41 with Beer–Lambert radiation added, of
which 12 are algebraic. The integrator's mass matrix is diagonal — one on the
differential rows, zero on the algebraic ones — so the semi-discrete system is the
semi-explicit DAE $\mathsf M\,d\boldsymbol u/dt = \boldsymbol f(\boldsymbol u, t)$
(§15.4).

**Why masses, not concentrations.** With $m_{j,i}$ as the unknown, the discrete mass
balance is exact by construction: a cell's mass changes only by its reaction sources
and its face fluxes, and every face flux leaves one cell and enters its neighbour. The
volume of a moving cell enters only through $\xi=m/V$, so a change of volume needs no
dilation term in any row (§11.3.4).

> **Implementation notes.** The state layout and its index functions (`mass_index`,
> `energy_index`, `position_index`, `velocity_index`, `intensity_index`,
> `history_index`, `area_index`, `lateral_rate_index`) are
> `StateLayout{NC,E,Mov,Rad,Chi,D,U,Car}` in `src/States/layout.jl`; the state mapping
> `derive!` (`src/States/derived.jl`) is the only reader of $\boldsymbol u$ and fills
> the derived state every kernel reads; the residual is `residual!` in
> `src/Assembly/residual.jl`. The layout, its row kinds and its index arithmetic are
> detailed in §13 and §15.1.

### 3.1.3 Overview of the coupled system

Collecting the laws derived below, the continuum model is the coupled system

$$
\begin{aligned}
&\text{Energy:}\quad
\partial_t\Big(\sum_j\xi_j h_j\Big) + \nabla\cdot\Big(\mathbf v_s\sum_j\xi_j h_j\Big)
+ \nabla\cdot\Big(\sum_g h_g\mathbf J_g\Big)
= -\nabla\cdot\mathbf q_{\mathrm{cond}} + Q_{\mathrm{rad}} + Q_{\mathrm{rxn}}
+ \sum_j\omega_j h_j, \\[4pt]
&\text{Gas species:}\quad
\partial_t\xi_j + \nabla\cdot(\mathbf v_s\xi_j)
= -\nabla\cdot\mathbf J_j + \omega_j, \\[4pt]
&\text{Condensed (solid/liquid) species:}\quad
\partial_t\xi_j + \nabla\cdot(\mathbf v_s\xi_j) = \omega_j,
\end{aligned}
$$

with $\mathbf v_s$ the skeleton velocity (zero on a fixed mesh, which does not
represent the volume change), $\mathbf J_g$ the gas flux relative to the skeleton and
$\omega_j=\sum_r\nu_{r,j}r_r$ the net production by reaction. The system is closed by
the constitutive laws inventoried in §3.8 (Fourier conduction, Arrhenius kinetics,
mixture properties, the ideal-gas pressure $P = N(\xi,T)/\varphi$, Darcy momentum
$v_g = -(\kappa/\mu)\nabla P$, and the volume-change/ALE relations). The energy
balance is written for the conservative sensible-enthalpy density, with every
component, gas included, in it; its product-rule expansion is the temperature
equation, and the integrator carries either the temperature or the cell's enthalpy
inventory (§3.3). Each term and its sign are derived in the sections that follow.

---

## 3.2 Condensed-phase mass and total density

### 3.2.1 Total density and pore density

The bulk total density is $\rho_{tot}=\sum_j\xi_j$, the total inventory per
bulk volume. Each gas has intrinsic pore density $\rho_j^g=\xi_j/\phi$.
These definitions serve different control volumes. Skeletal occupancy uses
condensed intrinsic densities without swelling weights; matrix transport
weights use condensed $\gamma_j\xi_j/\rho_{b,j}$. The volume-change law
uses its own per-component $\gamma_j/\rho_{b,j}$ (§3.6).

---

### 3.2.2 No separate bulk-continuity PDE

Pyrolysis.jl does **not** solve an independent continuity equation for $\rho$. The
cell masses are the unknowns, and a cell's total mass $\sum_j m_{j,i}$ is a derived
quantity whose rate is the sum of its mass rows (§3.2.3). Because every reaction is
mass-balanced at build time (the single reactant is consumed with stoichiometric mass
$1$ and the product yields sum to $1$, §6), the reaction contribution telescopes,

$$
\sum_{j=1}^{N_C} \omega_j
= \sum_r \Bigl(\underbrace{-1}_{\text{reactant}} + \underbrace{\textstyle\sum_{\text{products}}\nu}_{=\,1}\Bigr) r_r
= 0,
$$

so reactions move mass *between* components without creating or destroying it. A
cell's total mass changes only by its face fluxes; the column's changes only by the
fluxes through its boundary faces (and, when cells are deactivated on a fixed mesh, by
the gas a dead cell vents, which the event records; §14.4). This build-time invariant
is enforced by the `Reaction` constructor (`_validate_stoichiometry`, relative
tolerance $10^{-10}$); `mass_balance_error` and `verify_stoichiometry` report it (§6).

### 3.2.3 Phase-resolved mass: solids/liquids are immobile, only gases transport

The per-species conservation law is the foundation of the model. In integral
(finite-volume) form over a cell $i$ whose faces move with the skeleton,

$$
\frac{d}{dt}\!\int_{V_i}\!\xi_j\,dV
= -\oint_{\partial V_i}\!\mathbf J_j\cdot\mathbf n\,dA
+ \int_{V_i}\! \omega_j\,dV ,
$$

where $\mathbf J_j$ is the **mass flux** of component $j$ *relative to the skeleton*.
The decisive modeling choice is *which components transport*. In Pyrolysis.jl the
matrix is rigid relative to itself: solids and liquids do **not** flux between cells
— they are consumed/produced in place by reactions and (on a moving mesh) carried by
the mesh, never by molecular transport. Only gas components diffuse and advect through
the pore network, and they do so *relative to the skeleton*: the wall-friction, molecular-friction and Darcy closures return fluxes relative to the
pore walls, so the laboratory-frame gas flux is $\xi_g\mathbf v_s + \mathbf J_g$ and the
pore gas is carried with the pores it occupies (§11.1.2). Hence the **gas-vs-solid
dispatch**, in the frame of the skeleton:

$$
\boxed{\;
\partial_t\xi_j + \nabla\cdot(\mathbf v_s\xi_j) =
\begin{cases}
-\,\nabla\cdot \mathbf J_j + \omega_j,
& \text{if } j \text{ is gas},\\[6pt]
\omega_j,
& \text{if } j \text{ is solid or liquid.}
\end{cases}
\;}
$$

The continuous gas flux $\mathbf J_j$ follows the selected lumped or dusty-gas closure in §3.4,
§3.5 and detailed in §9; the reaction source $\omega_j = \sum_r \nu_{r,j}\,r_r$ is the
stoichiometry-weighted sum of reaction rates ($\nu_{r,j}<0$ reactant, $>0$ product;
§6). The discrete realization is one mass row per component and cell,

$$
\frac{dm_{j,i}}{dt} = \omega_{j,i}\,V_i + \dot M_{j,i},
\qquad
\dot M_{j,i} = -\sum_{f\in\partial i} s_{if}\,A_f\,F_{j,f},
\qquad \omega_{j,i} = \sum_r \nu_{r,j}\,r_{r,i},
$$

where the face flux $F_{j,f}$ is the gas-transport flux for a gas component and zero
for a condensed one. The gas-transport kernel flags its gas components once, when the
problem is set up (`map(is_gas, material.components)`), and forms face fluxes for
those alone: the divergence is never formed for the matrix phases. Where an
equidistribution correction moves the mesh relative to the skeleton (§3.6.3), every
component — solid, liquid and gas alike — also crosses the faces with the skeleton at
the relative velocity $-\mathbf c\cdot\mathbf n$, carried by a positivity-bounded
upwind face state (`RelativeFluxKernel`, §11.6); without a correction that kernel is
absent.

> **Implementation notes.** The gas mass fluxes are formed by `GasTransportKernel`
> (`src/Kernels/gas_transport.jl`), the fluxes relative to a corrected mesh by
> `RelativeFluxKernel` (`src/Kernels/relative_flux.jl`). The reaction sources
> $\omega_{j,i}V_i$ are added by `ReactionKernel` (`src/Kernels/reactions.jl`) from
> `species_source_terms!` in `src/Physics/kinetics.jl`; the per-reaction application
> enforces the single-reactant mass-balance convention in `_apply_reaction_source!`
> (same file). The row assembly `assemble_rows!` (`src/Assembly/residual.jl`) scatters
> every face flux as $-F_{j,f}A_f$ into its owner and $+F_{j,f}A_f$ into its neighbour
> and writes the mass rows.

---

## 3.3 The energy equation

### 3.3.1 Statement

Pyrolysis.jl solves the energy balance of a one-temperature porous medium for the
**sensible-enthalpy density** $\sum_j\xi_j h_j(T)$, with every component, gas
included, in it (derivation in §1.4.1 and §3.3.2):

$$
\boxed{\;
\partial_t\Big(\sum_j\xi_j h_j\Big) + \nabla\cdot\Big(\mathbf v_s\sum_j\xi_j h_j\Big)
= -\,\nabla\cdot \mathbf q_{\mathrm{cond}}
-\nabla\cdot\Big(\sum_g h_g\mathbf J_g\Big)
+ Q_{\mathrm{rad}}
+ Q_{\mathrm{rxn}}
+ \sum_j\omega_j h_j .
\;}
$$

The terms, with sign conventions, are:

| Term | Meaning | Continuous form | Sign / units |
|---|---|---|---|
| $\partial_t(\sum_j\xi_jh_j)+\nabla\cdot(\mathbf v_s\sum_j\xi_jh_j)$ | sensible-enthalpy storage, every component, following the skeleton | $h_j(T)=\int_{T_{\mathrm{ref}}}^{T}c_{p,j}\,dT'$ | W·m⁻³ |
| $-\nabla\cdot\mathbf q_{\mathrm{cond}}$ | Fourier conduction | $\mathbf q_{\mathrm{cond}} = -k_{\mathrm{eff}}\nabla T$ | net outflow cools |
| $-\nabla\cdot(\sum_g h_g\mathbf J_g)$ | sensible enthalpy carried by the gas flux | assembled at the faces; temperature-form reading $-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$ | W·m⁻³ |
| $Q_{\mathrm{rad}}$ | volumetric radiation absorption | model-dependent | W·m⁻³ |
| $Q_{\mathrm{rxn}}$ | reaction heat, at the local temperature | $-\sum_r h_r r_r$ | W·m⁻³ |
| $\sum_j\omega_jh_j$ | formation term | converts $h_r(T)$ to the sensible-enthalpy datum | W·m⁻³ |

Its product-rule expansion (§3.3.2) is the **temperature equation**

$$
\rho c_p^{\mathrm{eff}}(T,\xi)\,\frac{D_sT}{Dt}
= -\nabla\cdot\mathbf q_{\mathrm{cond}} - \sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T
+ Q_{\mathrm{rad}} + Q_{\mathrm{rxn}},
\qquad
\rho c_p^{\mathrm{eff}} = \sum_j\xi_jc_{p,j}(T).
$$

The discrete rows integrate the balance over a cell. With the cell's energy
right-hand side $E_i = (Q_{\mathrm{rxn},i}+Q_{\mathrm{rad},i})V_i - \sum_{f\in\partial i}
s_{if}A_fF^E_f$ and the face energy flux $F^E_f = q_f + \sum_jF_{j,f}h_j(T_f)$
(conduction plus the sensible enthalpy the face's mass fluxes carry at the face
temperature $T_f$),

$$
\text{temperature form:}\quad C_i\frac{dT_i}{dt} = E_i - \sum_jh_j(T_i)\,\dot M_{j,i},
\qquad
\text{inventory form:}\quad \frac{dH_i}{dt} = E_i + \sum_j\omega_{j,i}V_i\,h_j(T_i),
$$

with $C_i$ the cell's heat-capacity sum (§3.3.2) and $\dot M_{j,i}$ the face-flux part
of its mass rate (§3.2.3). This is precisely the row assembly of the residual: the face
energy fluxes are scattered into the owner and the neighbour, and the row is formed
from $E_i$ in the form the layout carries. The detailed discretization of each term —
the distance-weighted harmonic-mean face conductivity, the face temperatures of the
enthalpy flux, the order in which the kernels run — belongs to §7; here we justify the
*continuous* form.

> **Implementation notes.** The rows are formed by `assemble_rows!` and `_energy_row`
> in `src/Assembly/residual.jl`; $C_i$ is the derived state's `cp_sum` and the
> sensible enthalpies $h_j(T_i)$ its `h` (`derive!`, `src/States/derived.jl`); the
> conductive flux is formed by `ConductionKernel` (`src/Kernels/conduction.jl`).

### 3.3.2 The mixture heat capacity and where the gas enthalpy goes

The storage coefficient of the temperature equation is the **mixture** (solid +
liquid + gas) volumetric heat capacity,

$$
\rho c_p^{\mathrm{eff}}(T,\xi)
= \sum_{j=1}^{N_C} \xi_j\, c_{p,j}(T),
$$

and its discrete counterpart is the cell's heat-capacity sum
$C_i = \sum_j m_{j,i}\,c_{p,j}(T_{\mathrm{trial},i})$ [J·K⁻¹], over every component
and at the trial temperature (the cell temperature smoothly clamped into the property
bracket), smoothly floored at $\rho c_{p,\mathrm{floor}}V_i(0)$ with
$\rho c_{p,\mathrm{floor}} = 10^3$ J·m⁻³·K⁻¹ so that a cell emptied of all mass keeps a
finite rate. The pore gas is at the matrix temperature (one-temperature model,
§3.1.2) and its concentration $\xi_g$ is an integrated state, so its sensible storage
$\xi_g c_{p,g}\,\partial_t T$ is retained. Gpyro's default thermal-equilibrium exchange model
neglects gas storage and gas conduction (`FULL_QSG` retains both); the standard
FDS wall model has no pore-gas inventory. The size of our gas storage, $\varphi\rho_gc_{p,g}$ against the matrix heat capacity, is estimated
in §1.4.1 at about a quarter of a percent of $\rho c_p$ in a porous char at 1 atm.

Where the rest of the gas enthalpy goes follows from the derivation of §1.4.1,
summarized here. Write the enthalpy of component $j$ as $h_j^0 + h_j(T)$, with $h_j^0$
its formation enthalpy at $T_{\mathrm{ref}}$ and $h_j(T)=\int_{T_{\mathrm{ref}}}^T
c_{p,j}\,dT'$ its sensible enthalpy; let every condensed component move at the
skeleton velocity $\mathbf v_s$, and let $\mathbf J_g = \xi_g(\mathbf v_g - \mathbf v_s)$
be the gas flux relative to the skeleton (§3.2.3). The total enthalpy has no reaction
source; subtracting $h_j^0$ times each species balance from its balance leaves the
formation term $-\sum_jh_j^0\omega_j$, and with the heat of reaction at the local
temperature, $h_r(T) = \sum_j \nu_{r,j}\big(h_j^0 + h_j(T)\big)$, that term equals
$Q_{\mathrm{rxn}} + \sum_j\omega_jh_j(T)$: this is the balance of §3.3.1. Expanding
it with the product rule and the species balances
$\partial_t\xi_j + \nabla\!\cdot(\mathbf v_s\xi_j) = \omega_j - \nabla\!\cdot\mathbf J_j$
(with $\mathbf J_j = 0$ for condensed $j$), and using $\nabla h_g = c_{p,g}\nabla T$,
gives

$$
\Big(\sum_j \xi_j c_{p,j}\Big)\Big(\partial_t T + \mathbf v_s\!\cdot\!\nabla T\Big)
= -\nabla\!\cdot\mathbf q_{\mathrm{cond}}
-\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T
\underbrace{-\sum_r h_r(T)\,r_r}_{Q_{\mathrm{rxn}}}
+ Q_{\mathrm{rad}} ,
$$

the temperature equation of §3.3.1. Three things follow:

1. **The gas enthalpy flux is exact**, with $\mathbf J_g$ the skeleton-relative flux
   the transport closures return. It is the only place the gas flux enters the energy
   balance, and its temperature-form reading is the term
   $-\sum_g c_{p,g}\mathbf J_g\!\cdot\!\nabla T$.
2. **The heats of reaction are taken at the local temperature** — the DSC-style,
   inverse-calibrated convention (§6.6.3) — and the formation term
   $\sum_j\omega_jh_j(T)$ converts them to the sensible-enthalpy datum exactly; the
   composition-change enthalpy $\sum_jh_j\,\partial_t\xi_j$ is carried in full, by the
   formation term and by the gas enthalpy flux. A heat quoted instead at the datum
   $T_{\mathrm{ref}}$ differs from the local-temperature heat by Kirchhoff's law,
   $h_r(T) = h_r(T_{\mathrm{ref}}) + \sum_j\nu_{r,j}h_j(T)$; the code does not apply that
   conversion, so heats of reaction are supplied at the reaction temperature (§6.6.3).
3. **The skeleton-following derivative $\partial_t T + \mathbf v_s\!\cdot\!\nabla T$ is
   the mesh-frame rate on a mesh that moves with the skeleton** (§3.6, §11). The
   temperature equation then carries no advection term of any kind, and the statement
   holds unchanged on a fixed mesh with $\mathbf v_s = 0$. Where an equidistribution
   correction moves the mesh relative to the skeleton, the inventories, and the
   enthalpy they carry, cross the faces at the relative velocity (§11.2.2, §11.6).

The storage and gas temperature-advection terms correspond to ThermaKin Eq. 17
in 1D when the other closures match. Gpyro's `FULL_QSG` also retains gas
conduction, so it is not generally the same equation (§3.7).

> **Implementation notes.** The heat-capacity sum is `heat_capacity_sum`
> (`src/States/caloric.jl`), which has no phase branch and skips only a component
> whose mass is exactly zero; `derive!` (`src/States/derived.jl`) evaluates it at the
> trial temperature on the raw masses and floors it smoothly. The sensible enthalpies
> $h_j(T)$ are the components' enthalpy laws (`component_enthalpies`, same file;
> §4.4.10). The same functions feed the energy rows, the colored Jacobian (through
> dual-number evaluation), the enthalpy that a deactivation transfers to the receiving
> cell (§14.3.3) and the energy ledger (§16).

### 3.3.3 Conduction, radiation, and reaction-heat terms

**Conduction.** $\mathbf q_{\mathrm{cond}} = -k_{\mathrm{eff}}(T,\xi)\nabla T$ is
Fourier's law with the composition- and temperature-dependent effective conductivity
$k_{\mathrm{eff}}$ (mixing rules in §5). The minus sign in $-\nabla\cdot\mathbf q_{\mathrm{cond}}$
makes net conductive *outflow* cool the cell. At a face the flux is discretized with a
distance-weighted harmonic-mean conductivity (resistance-in-series), reducing to the
ordinary harmonic mean on a uniform mesh; the full derivation is in §7.

**Radiation.** $Q_{\mathrm{rad}}$ is the volumetric radiation absorption and is
entirely model-dependent. For `NO_RADIATION` and `SURFACE_ABSORPTION` it is zero —
the absorbed flux $a\,q_{\mathrm{inc}}$ of a `RadiativeFluxBC` enters through the
*surface energy balance* (boundary Newton solve, §12.3.5), not as a volumetric
deposit, and the two models give identical residuals. For `BEER_LAMBERT` it follows
the in-depth attenuation $q''' = \alpha I(z)$ with
$I(z) = I_{\mathrm{inc}}\exp(-\int\alpha\,dz')$: the face intensities are algebraic
unknowns, each satisfying $0 = I_{\mathrm{down}} - I_{\mathrm{up}}\,e^{-\alpha_i\Delta z_i}$
across its cell (and $0 = I_e - a\,q_{\mathrm{inc}}$ at the exposed face), and cell $i$
absorbs $I_{\mathrm{up}}A_{\mathrm{up}} - I_{\mathrm{down}}A_{\mathrm{down}}$, the exact
cell integral of $\alpha I$. The three models, their mutual exclusivity, and their
numerics are the subject of §8.

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
the smoothly-gated Arrhenius form of §6, evaluated at the trial temperature.

> **Implementation notes.** $Q_{\mathrm{rxn}} = -\sum_r h_r r_r$ is
> `reaction_heat_source` in `src/Physics/kinetics.jl`, with the negative
> sign documented in-source ("Note the negative sign!"). Reaction rates and heats are
> evaluated together by `evaluate_reactions!` (same file), called per cell by
> `ReactionKernel` (`src/Kernels/reactions.jl`), which adds $Q_{\mathrm{rxn},i}V_i$ to the
> cell's energy source. The Beer–Lambert rows and cell sources are formed by
> `BeerLambertKernel` (`src/Kernels/radiation.jl`).

### 3.3.4 Gas enthalpy flux and surface energy

Every interior gas flux carries $\sum_gF_{g,f}h_g(T_f)$, using the same
species flux as the mass rows and the common interface temperature. The
shared face value makes this contribution conservative. Subtracting
$\sum_jh_j(T_i)\dot M_{j,i}$ in a temperature row leaves the transported
sensible-enthalpy differences from the cell temperature. Both energy
formulations therefore use identical mass and carried-energy fluxes.

At a boundary, outward gas carries $h_g(T_s)$ and inward gas carries its
specified reservoir enthalpy. An LTE wall film has heat $q_s$ into the wall;
the reservoir control surface has $q_r=q_s+H_r-H_s$, where $H_s$ uses surface
enthalpy for every signed species flux. Thus
$E_{out}=H_r-q_r=H_s-q_s$. Without a thermal film, imposed heat and donor
enthalpy enter directly, including cold gas entering an adiabatic wall.
Chapter 9 states the conditional coupling and supported blowing domains.
No additional incoming-gas heating term is added. The internal-energy and
pressure-work storage convention is as stated in §1.4 and is independent of
this carried-enthalpy balance.

---

### 3.3.5 Temperature or enthalpy: the two energy states

The energy row of every cell carries one of two unknowns, chosen by `energy_state`
(`solve` and `setup` keyword; `TemperatureState`, `EnthalpyState` in
`src/States/layout.jl`):

- **Temperature form** (`:temperature`, the default). The unknown is $T_i$ and the row
  is $C_i\,dT_i/dt = E_i - \sum_jh_j(T_i)\dot M_{j,i}$. The temperature is directly
  observable, the property laws and the surface Newton solve read it with no
  inversion. ThermaKin2Ds and FDS use temperature, while the documented Gpyro method solves condensed enthalpy and then recovers temperature (§7.8).
  Energy is conserved to the integrator's order: the identity
  $dH_i/dt = E_i + \sum_j\omega_{j,i}V_ih_j(T_i)$ holds at every state, but a step's
  change of $H_i$ departs from the integral of its rate by the method's local error.
- **Inventory form** (`:enthalpy`). The unknown is $H_i = \sum_jm_{j,i}h_j(T_i)$ and
  the row is $dH_i/dt = E_i + \sum_j\omega_{j,i}V_ih_j(T_i)$. The inventory less its
  integrated rate is a linear invariant of the semi-discrete system, which every
  Runge–Kutta step conserves to round-off. The temperature is recovered in the state
  mapping by inverting the cell's caloric relation, regularized so that it is defined
  for every real state:
  $$
  H_i = \sum_j m_{j,i}^{+}\,h_j(T_i) + C_{\mathrm{reg},i}(m)\,(T_i - T_{\mathrm{ref}}),
  $$
  with $m^{+}_{j,i}$ the one-sided positive part of the mass (exactly zero at and below
  zero, exactly $m_{j,i}$ above the cell's mass tolerance $\delta_i$) and
  $C_{\mathrm{reg},i}$ a C¹ capacity switch that is zero while the cell's guaranteed
  heat capacity is at least twice its floor and rises to the floor as it vanishes. The
  relation is strictly increasing in $T_i$ with slope at least the floor, so it has
  exactly one root for every real $H_i$; a bracketed Newton iteration started at the
  cell's previous temperature finds it, and a dual-number inventory is differentiated
  through the defining relation (one Newton correction per nesting level), so the
  colored Jacobian carries the implicit-function derivative exactly
  (`CaloricExtension`, `invert_caloric`, `src/States/caloric.jl`). The relation agrees
  with the physical caloric relation on every state in which each component is absent
  or holds at least $\delta_i$.

The two rows are the same equation (the product rule with the mass rows turns either
into the other), and both conserve the sensible-enthalpy inventory: the energy ledger
closes exactly in both, to round-off in the inventory form and at the integrator's
order in the temperature form (§16.5.2). On a gasifying PMMA column (5 mm, 20 cells,
a 50 kW·m⁻² heat flux with re-radiation, 120 s, `abstol = 1e-8`, `reltol = 1e-6`),
audited with `ConservationRecorder` and `conservation_audit`, the relative energy
defect is $5.1\times10^{-5}$ in the temperature form and $1.6\times10^{-6}$ in the
inventory form, and the final surface temperatures differ by $4\times10^{-4}$ K; both
defects include the recorder's quadrature of the rates on the step interpolant
(§16.5.3), which is why the inventory form's is not at round-off.

The absolute tolerances follow the unknown: a temperature row takes `abstol` [K]; an
enthalpy row takes $\varepsilon_m\,m^s_i(0)\,c_{p,\mathrm{ref}}\,$`abstol` [J], so that a
cell at the deactivation threshold resolves its temperature to `abstol` as the
temperature form does (`abstol_vector`, `src/States/tolerances.jl`). The inventory
form pays one caloric inversion per cell per residual evaluation; it is the choice when
the energy ledger must close to round-off. In it the residual depends, at the
inversion's round-off, on the temperatures of the workspace's latest evaluation, from
which the inversion starts.

---

## 3.4 Gas-phase mass conservation and the pressure closure

For each gas, $\partial_t\xi_j+\nabla\cdot(\mathbf v_s\xi_j)=-\nabla\cdot\mathbf J_j+\omega_j$.
The skeleton-relative flux uses the material's `LumpedTransfer` or `DustyGas`
model (Chapter 9). A single face flux is shared by adjacent species rows.

### 3.4.1 The pressure potential

$\phi=1-\sum_{S\cup L}\xi_j/\rho_{i,j}(T)$ and
$p_j=\xi_jRT/(M_j\phi)$ define the intrinsic pore state. The gas EOS is
algebraic; gas storage remains in the evolved species inventories. It is
evaluated even with impermeable boundaries and zero transport coefficients.
For lumped transfer, $J_j^D=-\Lambda M_j\nabla p_j/(RT)$.

### 3.4.2 Physical and trial pore states

Skeletal occupancy is independent of the swelling factor. All gas-bearing
initial and accepted states must have resolved positive pore volume, including
uncertainty from mass, temperature and geometry. Trial evaluation uses a C¹
porosity divisor guard of $10^{-12}$, an identity above twice that value.
`pore_state` centralizes pressure, intrinsic densities and fractions. The
pressure-aware controller and current gas guard scale are described in §9.10.

---

## 3.5 Gas-phase momentum: Darcy's law

Both gas models include $v_D=-\kappa\nabla P/\mu$ and transport each species
with its donor intrinsic density $\rho_j^g=\xi_j/\phi$. The donor's own
porosity and temperature define this density. `κ=0` closes this mechanism;
wall transfer can still conduct. Material interfaces match total species flux
through an interface pressure state, including legal wall-only/Darcy-only
layer pairs. The common interface thermal state also supplies carried
enthalpy (§9.6). Darcy velocity is superficial and differs from total mass or
molar flow when wall/molecular transport is present.

---

## 3.6 Deformation: from ThermaKin's single-reference-point source to ALE mesh motion

### 3.6.1 The volumetric strain rate

As reactions consume condensed mass and (optionally) swell the matrix, the material
deforms. The local **volumetric strain rate** (relative rate of cell-volume change) is

$$
\boxed{\;
\theta_i = \frac{1}{V_i}\frac{dV_i}{dt}
= \sum_{j=1}^{N_C} \gamma_j\,\frac{\omega_{j,i}}{\rho_j(T_i)} ,
\;}
$$

derived from the mixture continuity relation: when component $j$ (with bulk density
$\rho_j$) changes its mass concentration by reaction at rate $\omega_j$, it frees or
claims volume at rate $\omega_j/\rho_j$, weighted by the swelling factor $\gamma_j$ that
decides whether the component participates in the matrix skeleton. Solids default to
$\gamma=1$ (depletion shrinks the matrix); gases default to $\gamma=0$ (they escape
without changing matrix volume); intumescent gases may take $\gamma=1$ (gas creates
matrix volume); partial values $\gamma\in(0,1)$ are allowed. The implementing kernel
`volume_change_rate` sums exactly $\sum_j \gamma_j\,\omega_j/\rho_j(T)$, skipping any
component with $\gamma_j\le 0$, at the trial temperature and with the *reaction* mass
sources alone as the driver: the mesh closure forms it from the kernels' reaction
sources before any external source is added, so a manufactured source of a
verification run does not move the mesh. The absolute rate is $dV_i/dt = V_i\theta_i$
on a mesh that moves with the skeleton.

> **Implementation notes.** $\theta = \sum_j \gamma_j\,\omega_j/\rho_j(T)$ is
> `volume_change_rate` (`_volume_change_rate`) in `src/Physics/volume_change.jl`
> (with `if γ > 0` guards), called per cell by `closure_rows!`
> (`src/Closures/balance.jl`). It is the swelling-weighted per-component form, not
> $-(\sum_j\omega_j)/\rho$.

### 3.6.2 The reference models' deformation formulations

ThermaKin2Ds writes its continuum balances in stationary coordinates. Its
deformation term is the complete divergence
$\partial_x[\xi_j\int_0^x\rho^{-1}\partial_t\rho\,dx']$ (TK Eq. 15),
with a stationary reference plane. Numerically, element thicknesses change and
elements are split or combined (TK §3). An Eulerian equation does not imply
fixed numerical element sizes.

Gpyro's 1D balances instead contain $\bar\rho Y_i\Delta z$ directly and
determine $\Delta z$ from mass and bulk-density consistency (GP Eqs. 3, 56–57).
Condensed mass does not cross these deforming cell boundaries. Neither model
is generally identical to our reaction-only strain closure: their density and
volume constitutive laws must also be matched (§10.11).

---

### 3.6.3 Pyrolysis.jl's recasting via ALE mesh motion

Pyrolysis.jl instead makes the strain drive an **Arbitrary Lagrangian–Eulerian** mesh
(`use_ale = true`). The node velocities are algebraic unknowns, and the mesh velocity
is the solid velocity plus a correction, $w = v_s + c$. The solid velocity $w - c$ is
held by one **balance row** per cell, the cell's volume balance at the solid velocity,
placed on the velocity row of its node farther from the pinned back node:

$$
0 = w_1,
\qquad
0 = A_{i+1}\,(w - c)_{i+1} - A_i\,(w - c)_i - (\theta_i - \theta_A)\,V_i
\quad (i = 1,\dots,n),
$$

with $A_i$ the area of the face at node $i$. On a planar line, where every face has the
area $A$ of the cross section, the rows reduce to the cumulative integral of the axial
strain from the pinned bottom node,

$$
v_{s,1} = 0,
\qquad
v_{s,i+1} = v_{s,i} + (\theta_i - \theta_A)\,\Delta z_i,
\qquad \Delta z_i = \frac{V_i}{A},
$$

where $\theta_A$ is the **column-global radial strain rate** subtracted to avoid
double-counting volume change that is instead expressed through the cross-section (see
below). The node positions advance as $dz_n/dt = w_n$, and the cell volumes are
computed from the positions ($V_i = A(t)\,|z_{i+1}-z_i|$ on a planar line), so
$dV_i/dt = \theta_iV_i + s_i$ holds semi-discretely, with $s_i$ the volume the
correction sweeps: the discrete Geometric Conservation Law holds by construction, not
by a correction (§11.5). Because the cell masses are the unknowns, a change of volume
enters only through $\xi = m/V$; no row carries a dilation term (§11.3.4). The balance
rows are bidiagonal from the pinned node, so the Jacobian stays banded.

While $c = 0$ the mesh moves with the skeleton and every component — the pore gas
included, since its transport flux is relative to the skeleton — rides with it, so the
mesh frame introduces **no advection term** in any equation. The correction is zero
unless an equidistribution correction is installed: by default `setup` installs one
where the mesh moves, cells are deactivated and some cell, converted completely, would
fall below the deactivation threshold, so that the cells of a vanishing layer thin
together and leave in one event (§11.7, §14.2.1). Under it every inventory — each
component's mass, the sensible enthalpy it carries and the conversion-history pair —
crosses the faces at the relative velocity $-c$ (§11.6). The mesh kinematics, the
geometric conservation law and the correction are derived in full in §11 (The ALE
Framework) and the strain rate in §10 (Volume Change & Lateral Shrinkage). This chapter
compares deformation strategies; equality of their constitutive volume laws
must be checked for the material and reaction mechanism (§10.11).

> **Implementation notes.** The balance rows are `BalanceRowClosure` and
> `closure_rows!` in `src/Closures/balance.jl`; the position rows $dz_n/dt = w_n$ and
> the velocity rows are written by `assemble_rows!` (`src/Assembly/residual.jl`); the
> correction is `EquidistributionCorrection` (`src/Closures/correction.jl`), installed
> by `_installed_correction` (`src/Integration/setup.jl`). The initial node velocities
> that satisfy the balance rows are computed by `consistent_state!` (same file).

### 3.6.4 Variable cross-section and the pyrolysis-progress state $\chi$

Lateral (radial) shrinkage or swelling is folded into the scalar cross-sectional area

$$
A(t) = A_0\cdot \mathrm{law}(\bar\chi),
$$

driven by the **column conversion**

$$
\bar\chi = \frac{\sum_i G_i}{\sum_i m_{\mathrm{dry},i}},
\qquad
\frac{dG_i}{dt} = S_{\mathrm{dry},i}\,V_i,
\quad
S_{\mathrm{dry},i} = \sum_{j\in\text{dry gas}} \omega_{j,i},
\quad
\frac{dm_{\mathrm{dry},i}}{dt} = 0 ,
$$

where $(G_i, m_{\mathrm{dry},i})$ is the **conversion-history pair** of §3.1.2: $G_i$ the
cumulative dry pyrolysis gas cell $i$ has produced, and $m_{\mathrm{dry},i}$ the cell's
initial solid mass, which has no source. The cell's pyrolysis progress is
$\chi_i = G_i/m_{\mathrm{dry},i}$, so $\bar\chi$ is the $m_{\mathrm{dry}}$-weighted mean
of the $\chi_i$. Water vapour is excluded from the dry-gas sum (it is
sorption/desorption, not pyrolysis; §4.6.4). Where an equidistribution correction moves
the mesh, the pair crosses the faces with the solid's mass flux (§11.6), so that it stays
attached to the material that produced it. The area and the **column-global radial
strain rate** $\theta_A = (1/A)\,dA/dt$ are algebraic unknowns with two rows,

$$
0 = A - A_0\,\mathrm{law}(\bar\chi_c),
\qquad
0 = \theta_A - \frac{d}{dt}\log\mathrm{law}(\bar\chi_c),
$$

the second the logarithmic rate of the first, taken by forward differentiation of
$\log\mathrm{law}$ along the column sums of the history rates
$(\sum_i\dot G_i, \sum_i\dot m_{\mathrm{dry},i})$ — by the chain rule
$\theta_A = \big(\mathrm{law}'(\bar\chi)/\mathrm{law}(\bar\chi)\big)\,d\bar\chi/dt$ — so the
area row's own rate is exactly $\theta_A A$. Here $\bar\chi_c$ is $\bar\chi$ smoothly
clamped into $[-10^{-3}, 1+10^{-3}]$ (the identity on every physical conversion), and
$\sum_i m_{\mathrm{dry},i}$ is smoothly floored at $10^{-8}$ of its initial value. The
law must be positive and finite on that interval and equal one at zero conversion,
which `LateralShrinkage` checks when the problem is set up; a lateral-shrinkage law
needs a moving mesh (`use_ale = true`) and planar coordinates. The two rows read every
cell's history pair (the second also their rates): they are the system's only dense
rows, and their Jacobian rows are assembled from six column sums after the colored pass
(§15.9). The lateral-shrinkage law, its identity/linear-shrink/swell variants, and the
first-principles-vs-empirical discussion are in §10.

> **Implementation notes.** The cross section is `LateralShrinkage`, its rows
> `lateral_rows` and their consistent values `consistent_lateral`
> (`src/Closures/lateral.jl`); the history sources (the dry gas of each cell's reaction
> sources) are `_history_sources!` (same file), the history and cross-section rows are
> written by `_history_rows!` and `_lateral_rows!` in `src/Assembly/residual.jl`.

---

## 3.7 Comparison to ThermaKin2Ds, Gpyro, and FDS

Source scope is given in §1.6.5. The following mappings are conditional on
matching material definitions, boundary conditions and reaction conventions.

### 3.7.1 Species / mass conservation

| Quantity | Pyrolysis.jl | ThermaKin2Ds | Gpyro 1D | FDS standard solid model |
|---|---|---|---|---|
| Species representation | cell masses; $\xi_j=m_j/V$ | continuum concentrations; numerical component masses | condensed mass fractions and mass per area; separate gas inventories | condensed component densities |
| Condensed motion | skeleton-following mesh; relative flux if corrected | deformation term and changing elements | no condensed flux through deforming cells | thickness update and remapping |
| Gas transport | lumped superficial transfer or dusty-gas friction, with optional Darcy | volume-fraction-gradient law | Darcy plus gas mass-fraction diffusion | instantaneous release to CFD gas |

TK Eq. 13 uses a **volume-fraction**, not a mass-fraction, gradient.
Its reduction to $J=-(\lambda/T)\partial_x(\xi T)$ requires component
reference density proportional to $1/T$ at fixed defining pressure.
GP Eq. 60 uses a different mass-fraction diffusion law (§9.10).

### 3.7.2 Energy equation

ThermaKin2Ds stores all components at one temperature (TK Eq. 17).
Pyrolysis.jl's temperature equation has the corresponding storage and gas
temperature-advection terms in 1D. Full equivalence also requires matching
radiation, deformation, property laws and heat-of-reaction basis.

Gpyro solves condensed sensible enthalpy, recovers temperature, and optionally
solves a separate gas energy equation (GP Eqs. 61, 70 and §6.8).
`FULL_QSG` restores gas storage **and gas conduction** in thermal equilibrium
(Eq. 65), so keeping gas heat capacity alone does not reproduce that option.
FDS has no pore-gas inventory but retains reaction and gas/solid enthalpy
exchange accounting (§7.8.4).

### 3.7.3 Control-volume conservation forms

Gpyro's $\bar\rho Y_i\Delta z$ times a fixed area corresponds to a
condensed species inventory, not to the total solid-plus-gas mixture density.
All four methods account for changing material inventory and geometry in
different ways. Those shared conservation requirements do not imply identical
volume-change constitutive laws (§10.11) or identical remapping errors.

### 3.7.4 Reaction kinetics

ThermaKin2Ds permits two concentration factors. Gpyro's heterogeneous rates
use conversion history and optional local oxygen dependence. FDS uses a
reactant density power with optional oxygen and temperature powers.
Pyrolysis.jl uses one reactant concentration power, smooth temperature gates
and a depletion limiter. Only the matched special cases of §6.9 have the same
unregularized rate law; there is no general automatic parameter conversion.

---

## 3.8 Closure inventory: which constitutive law closes each equation

The PDE system of §3.1.3 is incomplete until every flux and source is expressed in terms
of the state — the masses and the energy state, and the node positions and velocities,
face intensities, history pair and cross-section unknowns where they are carried. The
table below inventories the closures and points to the chapter where each is derived.
This is the contract that makes the residual a *pure function* of the state.

| Governing equation / term | Closure | Constitutive law | Chapter |
|---|---|---|---|
| Total density $\rho$ | exact | $\rho=\sum_j\xi_j$, $\xi_j = m_j/V$ | §3.2, §5 |
| Storage $\rho c_p^{\mathrm{eff}}$ | property + mixing | $\sum_j\xi_j c_{p,j}(T)$, every component; discretely $C_i$, smoothly floored | §4, §5 |
| Sensible enthalpy $h_j$ | caloric law | $\int_{T_{\mathrm{ref}}}^{T}c_{p,j}\,dT'$, closed form or tabulated | §4.4.10 |
| Conductive flux $\mathbf q_{\mathrm{cond}}$ | Fourier + mixing | $-k_{\mathrm{eff}}(T,\xi)\nabla T$, guarded harmonic-mean face $k$ | §5, §7 |
| Reaction rate $r_r$ | Arrhenius | $A_r e^{-E_r/R_gT}\prod_j\xi_j^{n_{r,j}}$, tanh gates, one-sided limiter, at the trial temperature | §6 |
| Species source $\omega_j$ | stoichiometry | $\sum_r\nu_{r,j}r_r$ | §6 |
| Reaction heat $Q_{\mathrm{rxn}}$ | enthalpy of reaction at the local temperature | $-\sum_r h_r(T,\xi)\,r_r$ | §6 |
| Radiation source $Q_{\mathrm{rad}}$ | radiation model | surface balance / Beer–Lambert face intensities | §8 |
| Gas enthalpy flux | conservative face flux | interior face-temperature closure; outward surface enthalpy and specified incoming reservoir enthalpy | §7, §12 |
| Wall and molecular gas flux | lumped transfer or dusty-gas friction | partial-pressure driving force, one consistent face pressure, total-flux interface matching | §9 |
| Gas advective flux | donor pore density × Darcy velocity | $p_{j,up}M_jv/(RT_{up})$ in each half-cell | §9 |
| Pressure $P$ | ideal gas in pores | $N(\xi,T)/\varphi$ | §9 |
| Porosity $\varphi$ | intrinsic density | $1-\sum_{S/L}\xi_j/\rho_{i,j}$, smoothly floored at $10^{-12}$ where divided by | §5 |
| Gas velocity $v_g$ | Darcy | $-(\kappa/\mu)\nabla P$, half-cell mobility $\kappa/\mu$, interface pressure from total species-flux continuity | §9 |
| Effective $k_{\mathrm{eff}},\lambda_{\mathrm{eff}},\kappa_{\mathrm{eff}},\varepsilon_{\mathrm{eff}},a_{\mathrm{eff}},\alpha_{\mathrm{eff}}$ | mixing rules | parallel/series/weighted/Bruggeman/Carman–Kozeny | §5 |
| Volume-change rate $\theta$ | swelling continuity | $\sum_j\gamma_j\omega_j/\rho_j(T)$, reaction sources only | §10 |
| Node velocity $w$ | balance rows | $0=A_{i+1}(w-c)_{i+1}-A_i(w-c)_i-(\theta_i-\theta_A)V_i$, $0 = w_1$ | §11.3.2, (11.9)–(11.10) |
| Node correction $c$ | equidistribution | installed by `setup`'s rule, or given | §11.7 |
| Cross-section $A$, $\theta_A$ | algebraic rows | $A=A_0\,\mathrm{law}(\bar\chi_c)$, $\theta_A = \frac{d}{dt}\log\mathrm{law}(\bar\chi_c)$ | §10, §11.3.3 |
| Boundary fluxes $q_{\mathrm{BC}}, J_{j,s}, v_{g,s}$ | thermal/mass/pressure BCs | surface Newton balance, film resistance, outward Darcy velocity | §12 |
| Initial state | initial conditions | $T_{\mathrm{initial}}$, $\xi_{j,\mathrm{initial}}$ at the cell centres; $m_{j,i}=\xi_jV_i$; energy state from $T$; algebraic unknowns made consistent | §12.8 |

With every entry of this table substituted, the right-hand side of the system is a
fully determined function $\boldsymbol f(\boldsymbol u, t)$, and the algebraic rows are
fully determined constraints. The unified `residual!` evaluates exactly this function in
one sequence for every configuration: the state mapping `derive!` (concentrations,
temperatures, their trial clamp, heat-capacity sums, sensible enthalpies, cell
properties), the node corrections, the kernels in their fixed order, the mesh closure's
rows (strain rates, velocity rows), the external sources, and the row assembly; options
that are off contribute no kernel and no rows (§1.3). The discretization of the
operators (finite-volume face fluxes, face means and face temperatures), the temporal
integration (the Rosenbrock method `Rodas5P` on the mass-matrix DAE, with the stage
matrix $W = J - \mathsf M/(\gamma_{\mathrm{RK}}\Delta t)$), and the colored Jacobian are
the subjects of §13–§15.

---

## 3.9 Limitations and open issues

The continuum model rests on several approximations that are exact only in limiting
regimes. The most important, with the section that quantifies or discusses each:

- **One temperature.** Matrix and in-pore gas are assumed in local thermal equilibrium;
  there is no two-temperature gas energy equation. This breaks for very fast, hot gas
  throughput, where the pore gas does not equilibrate with the walls within a cell.
- **Darcy momentum closure.** Pressure is evaluated from the evolving gas
  masses through $P=N/\varphi$ (§3.4.2), with a porosity floor. Compressible
  gas storage is retained. A separate pressure-evolution equation, as in
  Gpyro, is an alternative choice of dependent variable; neither Darcy model
  resolves acoustic momentum waves (§9.10).
- **First-order upwind advection** for the gas advective flux (§3.5.1); no TVD or
  MUSCL limiter on the gas advection. The scheme is monotone but numerically
  diffusive where the cell Péclet number exceeds one, and its upwind switch is a kink
  in the residual (§17.7.6).
- **Single reactant per reaction.** Multi-reactant (bimolecular) kinetics, supported by
  ThermaKin, is not implemented; the per-reactant tuple would need extending to carry a
  mass stoichiometry (§6). There is also no temperature-dependent $A$ or pressure
  dependence, and no distributed-activation-energy tooling.
- **Scalar lateral shrinkage.** Radial deformation is captured only through the single
  scalar $A(t)=A_0\,\mathrm{law}(\bar\chi)$ (§3.6.4); there is no resolved transverse
  coordinate, so genuinely 2D/axisymmetric effects (ThermaKin's cylindrical Laplacian
  applied to a shrinking sample) are not represented. The functional form of
  $\mathrm{law}(\chi)$ is user-supplied and largely empirical; a first-principles
  closure (e.g. percolation-based) is an open question (§10). Cross-section collapse is
  excluded only by requiring the law to be positive on $[-10^{-3}, 1+10^{-3}]$.

These limitations are inherited by every downstream chapter; where a later chapter
relaxes or quantifies one of them (e.g. §16 for the mass and energy ledgers, §18 for
regime-validity checks), the cross-reference is given there.
