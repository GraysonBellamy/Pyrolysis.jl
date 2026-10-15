# 11. The Arbitrary Lagrangian–Eulerian (ALE) Framework

Pyrolysis is, intrinsically, a *moving-domain* problem: the condensed phase recedes
as it gasifies, intumescent chars swell, and the position of the exposed surface is
itself an unknown. This chapter develops the Arbitrary Lagrangian–Eulerian (ALE)
formulation that Pyrolysis.jl uses to track that motion without solving a moving-
boundary PDE. We derive the ALE time-derivative identity from the kinematics of a
moving sample point, specialize it to the energy and species equations (every
component rides the skeleton, so the mesh frame adds only the dilation), construct the
cumulative mesh velocity from the per-cell volumetric strain rate, prove the
discrete Geometric Conservation Law (GCL) and its time-varying-cross-section
generalization, document the TVD-limited upwind discretization of the
mesh-relative temperature advection, and the conservative remapping operators used
when the mesh topology changes during adaptive refinement (§14). Throughout, the
governing relations are tied line-by-line to the implementation in
`src/Physics/{mesh_motion,gcl,volume_change,remapping}.jl`,
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
> AMR-inside-ALE — is *not representable* without first adding an advection
> operator for the slip velocity $w - v_{\text{solid}}$. The verification
> hook `Workspace.w_override` (which replaces the reaction-driven $w$ with a
> prescribed $w(z,t)$ for MMS runs, §18.2) inherits exactly this constraint:
> under prescribed non-material motion a spatially varying lab-frame-static
> field has a $w\,\partial_z\xi$ term with no discrete operator to represent
> it, so the manufactured-solution driver assembles every species source in
> the material-following frame (see test/verification/exact/mms.jl and
> §11.10 item 7). The temperature field, by contrast, *does* get its slip
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
(`ale_fluxes.jl:8–10`). Here $\partial/\partial t|_{\chi}$ is the **ALE
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
$d(\xi V)/dt = SV - [JA]_{\text{faces}}$ exact (§11.3).

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
surface transpiration term is rejected at construction because it would count
the gas heating a second time (Chapters 7, 12).

---

## 11.3 Cumulative mesh velocity

The ALE freedom is fixed by *choosing* $w(z,t)$. Pyrolysis.jl chooses it so that the
mesh exactly absorbs the local volume change of the material, with the substrate
pinned. This section derives the cumulative-integral formula, its discrete form, and
the cell-centered interpolation used by the advection kernels.

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
(`volume_change.jl:93–105`). Here:

- $\gamma_j$ is the **swelling factor** of component $j$ (component field
  `swelling_factor`): $\gamma=1$ for solids by default (when a solid depletes, the
  matrix loses that volume), $\gamma=0$ for gases by default (escaping gas does not
  change the skeleton volume), and $\gamma\in[0,1]$ for intumescent/partial
  contributors. Only components with $\gamma_j>0$ enter the sum (the implementation
  early-exits on $\gamma_j\le 0$, `volume_change.jl:17,28–36,48–61`).
- $\rho_j(T)$ is the **bulk** (pure-phase) density of component $j$, evaluated at the
  cell temperature; the temperature dependence is honored through the property
  function `comp.density(T)` (`volume_change.jl:19,29,…`).
- $d\xi_j/dt$ is the **species reaction source term** $S_j$ — the *reaction* part of
  the species RHS, supplied as the pre-computed `dξ_rxn` matrix
  (`residual.jl:389–399`). It is negative for a consumed solid and positive for a
  produced char/gas.

**Physical reading and sign.** When a solid of density $\rho_j$ is consumed
($d\xi_j/dt<0$, $\gamma_j=1$), the term $\gamma_j(1/\rho_j)(d\xi_j/dt)<0$, so
$\theta<0$: the cell shrinks. A net char-producing or swelling reaction with the
appropriate $\gamma$ can yield $\theta>0$ (expansion/intumescence). Note that
(11.6) is the *swelling-factor-weighted* form; the loose statement "$\theta =
-(\sum_j d\xi_j/dt)/\rho_{\text{total}}$" sometimes used as shorthand
(`problem_residual.md`) is only correct in the special case of a single solid with
$\gamma=1$ and constant $\rho$ — the code always uses (11.6).

> **Implementation note.** Relative strain rate: `volume_change_rate` /
> `_volume_change_rate` (dispatch-unrolled for $N_C\in\{1,2,3\}$, `@generated` for
> $N_C>3$), `src/Physics/volume_change.jl`. Batch absolute rate
> $dV_i/dt = V_i\theta_i$: `compute_cell_volume_rates!`, same file (lines 126–149).

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
strain is $\theta_L = \theta - \theta_A$ (`residual.jl:407`,
`volume_change.md:77`). Discretizing the integral cell-by-cell with the bottom node
pinned gives the **cumulative mesh-velocity recursion**
$$
\boxed{\;
w_1 = 0,\qquad
w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i
\;}
\tag{11.8}
$$
with cell thickness $\Delta z_i = V_i/A_{\text{section}}$
(`residual.jl:387,407–411`). Equivalently, $w$ at the right node of cell $i$ is the
exclusive prefix sum
$$
w_{i+1} = \sum_{k=1}^{i}(\theta_k - \theta_A)\,\Delta z_k
= \sum_{k\le i}\theta_{L,k}\,\Delta z_k .
\tag{11.9}
$$

The bottom-up sweep (loop `cell_idx = 1:n_cells`, `residual.jl:385`) starts at
`w[1]=0` (substrate pinned, `residual.jl:382`) and accumulates
`cumulative_velocity += θ_L * Δz_cell`, assigning the running total to the cell's
right node (`residual.jl:407–411`). Because the strain rates of all cells below a
node are summed, the top node moves at $w_{n+1} = \sum_k \theta_{L,k}\Delta z_k$,
i.e. the surface velocity equals the total axial volume-change rate of the column,
exactly tracking the receding (or expanding) surface.

**Substrate boundary condition.** Pinning $w_1=0$ enforces zero normal mesh velocity
at the impermeable substrate: the bottom node never moves, all motion is taken up
above it, and the material "stacks up" toward the surface as it is consumed
(`ale.md §6.2`, `volume_change.md` design subtlety 5). This is the natural choice
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
(`residual.jl:519–564`, `_compute_theta_A`). The column-average progress rate is the
volume-weighted dry-gas production normalized by the total initial dry-solid mass,
$$
\frac{d\bar\chi}{dt} = \frac{1}{M_{\text{dry}}}\sum_{i=1}^{n}
\Big(\sum_{j\in\text{dry gas}} S_{j,i}\Big)\,V_i ,
\qquad
M_{\text{dry}} = \sum_i m_{\text{dry},i},
\tag{11.11}
$$
computed at `residual.jl:549–557`, with $\mathrm{law}'(\bar\chi)$ obtained by
automatic differentiation `ForwardDiff.derivative(law, χ̄)` (`residual.jl:562`).
Guards: $\theta_A=0$ for identity-law materials, for missing `dξ_rxn`, for
$M_{\text{dry}}\le\epsilon$, and for $\mathrm{law}(\bar\chi)\le\epsilon$ (preventing
division by a vanishing cross-section, `residual.jl:539–541,560`). For the default
constant-cross-section material the dispatch returns exactly `0.0`
(`residual.jl:293–302`, `NoChi`), and (11.8) reduces to the pure cumulative-strain
integral. The interaction of $\theta_A$ with $\chi$-block dynamics and the lateral
law is developed fully in Chapter 10; here it enters only as the scalar subtracted
in (11.8).

### 11.3.4 Cell-centered velocity and the dilation correction

The advection kernels need a *cell-centered* mesh velocity. It is the arithmetic
mean of the bounding node velocities:
$$
\boxed{\;
w_{\text{cell},i} = \tfrac12\big(w_{\text{left}} + w_{\text{right}}\big)
= \tfrac12\big(w_i + w_{i+1}\big)
\;}
\tag{11.12}
$$
(`ale_fluxes.jl:227,449`). Node velocities are also interpolated to faces (faces
co-locate with nodes in 1D) by $w_{\text{face},f} = w_{\text{node}[f.\text{node\_id}]}$
(`residual.jl:417–423`, `_interpolate_node_to_face!`).

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
`geometry.md §10`, `state_layout`), packed in the $z$-block at offset $n(1+N_C)$.
Their evolution equation is simply the mesh velocity,
$$
\boxed{\;
\frac{dz_i}{dt} = w_i
\;}
\tag{11.15}
$$
written into the residual as `du[zoff+i] = ws.w_work[i]` (`residual.jl:204–207`).
The time integrator (Chapter 15) advances $z_i$ together with $(T,\xi,\chi)$ as one
coupled stiff system; the structured Jacobian carries the $z$-*column* couplings
in the geometry AD block and the $z$-*row* dependence on $(T,\xi)$ as embedded
prefix-sum couplings (plus the WithChi $\theta_A$ rank-1 channels) emitted by
the mesh-velocity provider (Chapter 15, §15.8.4 and §15.11.4).

### 11.4.1 Explicit node advance (mesh-motion utility)

For the diagnostic / mesh-management path (`move_mesh!`), node positions are
advanced by an explicit forward-Euler step,
$$
z_i^{\text{new}} = z_i^{\text{old}} + \Delta t\, w_i ,
\tag{11.16}
$$
realized by constructing a fresh immutable `Node1D` per node (`mesh_motion.jl:44–56`;
`ale.md §1.1`). The reference position $z_{\text{ref}}$ is left untouched so that
displacements $z_i - z_{\text{ref},i}$ remain available (`mesh_motion.jl:267–276`).
Within the implicit ODE solve, by contrast, $z_i$ is evolved through (11.15) by the
integrator's own (higher-order, implicit) scheme — the explicit (11.16) is the
formula used by the standalone `move_mesh!` and by the predictive GCL check.

### 11.4.2 Geometric volume recompute

After the nodes are at their new positions, all cell volumes are recomputed
**analytically from geometry**, not carried forward by a flux formula:
$$
\boxed{\;
V_i^{\text{new}} = A(t)\,\big|z_{\text{right}} - z_{\text{left}}\big|
= A(t)\,|z_{i+1}-z_i|
\;}
\tag{11.17}
$$
with the current (possibly $\chi$-updated) cross-section $A(t)$. Cell centers are the
node midpoints,
$$
z_i^c = \tfrac12\big(z_i + z_{i+1}\big),
\tag{11.18}
$$
and in 1D all face areas equal the global cross-section, $A_f = A(t)$
(`mesh_motion.jl:62–83,100–106,113–134`; `geometry.md §6.1`,
`geometric_quantities.jl:20–46`). During residual evaluation the same geometry is
rebuilt directly from the state-vector node positions, type-stably (carrying
`eltype(u)` so AD propagates through geometry), by
`update_state_geometry!(…, ::ALE)` (`state_cache.jl:336–354`,
`volume_change.md:153–166`) and `compute_geometry_from_state!`
(`geometry.md §7.2`). This keeps the residual a pure function of $u$ that never
mutates the mesh.

**Why geometric recompute rather than a flux update?** Recomputing $V = A(t)\,L$
from current node positions is robust to a simultaneous change in the cross-section
$A(t)$ (variable-area mode): carrying $V_{\text{old}}$ across an $A(t)$ change via a
pure axial-flux formula would freeze the old radial area into the volume. The
geometric path bakes the current $A(t)$ into every volume automatically. As proved
in §11.5, when geometry is computed correctly this path *still* satisfies the GCL to
machine precision (`mesh_motion.jl:68–71` comment; `ale.md §6.1`).

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
\tag{11.19}
$$
(`gcl.jl:1–9`; `geometry.md §8.1`, `compute_gcl_volume_rate`). Integrating over
$[t,t+\Delta t]$ gives the **discrete GCL**:
$$
\boxed{\;
V_i^{\text{new}} - V_i^{\text{old}} = \Delta t\,(w_{\text{right}} - w_{\text{left}})\,A
\;}
\tag{11.20}
$$
(`gcl.jl:24`). Why this matters: a GCL violation $\delta V$ acts as a *spurious mass
and energy source*. If the discrete volume does not change by exactly the swept
volume, the cell appears to gain or lose conserved quantity from nothing. Even tiny
violations $\delta V/V\sim 10^{-6}$ accumulate over thousands of time steps and can
dominate the genuine reaction-driven mass loss (`ale.md §7.2`). Exact GCL is
therefore non-negotiable for long-time accuracy.

### 11.5.2 Exactness of the geometric path

Pyrolysis.jl's `move_mesh!` (with `use_gcl_volumes=true`, the default) does **not**
apply the flux formula (11.20) directly; it recomputes $V_i = A|z_{i+1}-z_i|$ from
the advanced nodes (11.17). We show these coincide to machine precision for constant
$A$. After the explicit advance $z_{\text{left}}^{\text{new}} =
z_{\text{left}}+\Delta t\,w_{\text{left}}$ and $z_{\text{right}}^{\text{new}} =
z_{\text{right}}+\Delta t\,w_{\text{right}}$,
$$
V_i^{\text{new}} = A\,(z_{\text{right}}^{\text{new}} - z_{\text{left}}^{\text{new}})
= A\big[(z_{\text{right}}-z_{\text{left}}) + \Delta t\,(w_{\text{right}}-w_{\text{left}})\big]
= V_i^{\text{old}} + \Delta t\,(w_{\text{right}}-w_{\text{left}})\,A ,
$$
which is exactly (11.20). The only deviation is floating-point round-off in the
subtraction and multiplication, of order $10^{-15}$ relative (`ale.md §2.3, §6.1`).
Hence the geometric path is GCL-exact by algebraic identity, while remaining robust
to $A(t)$ changes. The explicit flux update is also available as the inline helper
$$
V_{\text{new}} = V_{\text{old}} + \Delta t\,(w_{\text{right}} - w_{\text{left}})\,A
$$
(`gcl_volume_update`, `gcl.jl:183–191`) for custom schemes, but it is *not* the path
`move_mesh!` takes.

### 11.5.3 GCL diagnostics

Several verification functions monitor GCL compliance:

- **Predictive check** `gcl_error(mesh, w, dt)` (`gcl.jl:34–77`): assuming the motion
  has *not yet* occurred, computes for each cell the GCL-predicted change $\Delta
  V_{\text{GCL}} = \Delta t\,(w_R-w_L)A$ and the would-be geometric change
  $\Delta V_{\text{actual}} = |z_R^{\text{new}}-z_L^{\text{new}}|A - V_i$, and
  returns the maximum relative error $\max_i |\Delta V_{\text{actual}} - \Delta
  V_{\text{GCL}}|/|V_i|$. Expected value $\approx 10^{-15}$.

- **Post-motion check** `gcl_error_post_motion(mesh, V_old, w, dt)`
  (`gcl.jl:96–131`): after `move_mesh!`, compares the *actual* recorded $V_i^{\text
  {new}}-V_i^{\text{old}}$ against $\Delta t\,(w_R-w_L)A$, returning the maximum
  relative error.

- **Boolean** `check_gcl(mesh, w, dt; tol=1e-12)` (`gcl.jl:152–155`): returns
  `gcl_error(...) < tol`; the default tolerance is $10^{-12}$.

- **Wrapped motion** `move_mesh_gcl!` (`mesh_motion.jl:354–381`): snapshots
  $V_{\text{old}}$, calls `move_mesh!`, then `gcl_error_post_motion`, and `@warn`s
  if the error exceeds `gcl_tol` (default $10^{-10}$).

### 11.5.4 Generalized GCL for time-varying cross-section

When the cross-section changes during the step ($A_{\text{old}}\ne A_{\text{new}}$,
WithChi), the swept volume has two contributions: the axial face motion *and* the
radial area change at constant length. The generalized discrete GCL is
$$
\boxed{\;
V_i^{\text{new}} - V_i^{\text{old}}
\approx
\underbrace{\Delta t\,(w_{\text{right}} - w_{\text{left}})\,A_{\text{new}}}_{\text{axial face flux}}
\;+\;
\underbrace{(A_{\text{new}} - A_{\text{old}})\,L_{\text{old}}}_{\text{radial area change}}
\;}
\tag{11.21}
$$
with $L_{\text{old}} = V_i^{\text{old}}/A_{\text{old}}$ (`gcl.jl:226–278`,
`gcl_residual`; `ale.md §2.7`). The helper `gcl_residual(mesh, V_old, w, dt;
A_old, A_new)` returns the maximum relative residual of (11.21) across active cells;
pass `A_old`, `A_new` explicitly when verifying a step that crossed an $A(t)$ change.
For $A_{\text{old}}=A_{\text{new}}$ the second term vanishes and (11.21) collapses to
(11.20).

### 11.5.5 GCL tracking over a run

The `GCLTracker` (`gcl.jl:296–376`) accumulates the GCL error across the whole
simulation:
$$
\text{max\_error} = \max_n e_n, \qquad
\text{cumulative\_error} = \sum_n |e_n|\,\Delta t_n = \int_0^T|e(t)|\,dt ,
$$
counts steps with $e_n>$ tolerance (default $10^{-10}$), and records a
`(time, error)` history (`update!`, `gcl.jl:328–339`; `report`, `gcl.jl:354–363`).
The time-integrated norm $\int|e|\,dt$ is sensitive to *systematic* (rather than
transient) bias, which is precisely what threatens long-time conservation
(`ale.md §6.7`). In a healthy run per-step errors sit at $10^{-12}$–$10^{-14}$;
sustained errors above $10^{-10}$ signal geometric degeneracy or an ill-conditioned
mesh (Chapter 16, Conservation Diagnostics).

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
\tag{11.22}
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
\tag{11.23}
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
The default for ALE advection is **Van Leer** (`limiter::Limiter = VANLEER`;
`ale.md §4.4`). Guards: if there is no donor-donor cell, or the limiter is
`NONE`, or $|\text{grad}_0|<10^{-20}$ (flat far slope — the ratio is undefined),
the scheme falls back to first-order $\text{grad}_1$.

### 11.6.3 Boundary treatment: zero-gradient ghost and the boundary gradient distance

At a domain boundary the "neighbor" cell does not exist (id $=0$). Two choices make
the stencil well-posed:

1. **Zero-gradient ghost value.** A boundary ghost takes the boundary cell's own
   value, $f_{\text{ghost}} = f_{\text{cell}}$ (`_get_boundary_temperature`,
   `ale_fluxes.jl:257–263`; analogous for $\xi$ via `col_upwind = i` when the upwind
   neighbor is the boundary, `ale_fluxes.jl:487`). This is a Neumann (no-flux)
   condition for scalar advection — physically appropriate, since there is no
   advective transport across the closed domain edge.

2. **Boundary gradient distance.** When the upwind neighbor is a boundary, the
   distance returned by `_upwind_gradient_distance` is the full cell thickness,
   $$
   d = \frac{V_i}{A_{\text{section}}} = \Delta z_i,
   $$
   (`ale_fluxes.jl:288–297`). This is precisely `face_distance` evaluated at the
   upwind face: the **center-to-center** spacing $|z_i^c - z_{\text{neighbor}}^c|$ for
   an interior neighbor, and the center-to-ghost-cell spacing
   $\Delta z_i = 2\times(\text{center-to-face})$ for a boundary. Reusing the genuine
   center-to-center distance (rather than a naive cell width) preserves first-order
   consistency on non-uniform meshes — a mismatched width would introduce a
   zeroth-order $d/\Delta z$ bias that does not vanish as $h\to 0$. At a boundary the
   value of $d$ is in any case **immaterial**: the closed edge uses a zero-gradient
   ghost (the upwind value is set equal to the cell's own value), so the first-order
   gradient is identically zero there and the TVD branch is disabled. Consequently
   $d$ divides a zero numerator in the residual and the $+w/d$ and $-w/d$ face-local
   AD partials cancel exactly in the structured-Jacobian scatter — no residual,
   Jacobian, or dense-reference entry depends on it.

> **Implementation note.** Energy advection (override-only):
> `compute_ale_advection_term!` with the shared limited gradient
> `_limited_gradient`, in `src/Discretization/ale_fluxes.jl`. The species
> dilation is `_apply_ale_dilation_correction!` in `src/Residual/residual.jl`
> (§11.3.4).

---

## 11.7 Conservative remapping (for AMR / mesh changes)

When the mesh *topology* changes — adaptive refinement or coarsening (Chapter 14),
or surface-cell depletion — the cell-centered state must be transferred from the old
mesh to the new one. ALE node motion (§11.4) handles *continuous* deformation; it
does not handle a change in the *number* of cells. The remapping operators in
`src/Physics/remapping.jl` perform that transfer while preserving the conserved
integrals (mass per species, thermal energy).

### 11.7.1 State layout for remapping

The remap operates on a **cell-interleaved** flat vector,
$$
u = [\,T_1, \xi_{1,1}, \dots, \xi_{N_C,1},\; T_2, \xi_{1,2},\dots,\; \dots\,],
\qquad \text{idx}(i) = (i-1)(1+N_C),
$$
with $T$ first then the $N_C$ concentrations *per cell* (`remapping.jl:108,148–152`).
This differs from the *block-major* ODE state vector of Chapter 13 (all $T$, then
all $\xi_1$, …); the remap routines use their own cell-interleaved convention and
the caller is responsible for packing/unpacking accordingly.

### 11.7.2 Overlap geometry

The unifying primitive is the **overlap volume** of a new cell $i_{\text{new}}$ and
an old cell $i_{\text{old}}$. With cell extents $[z_L,z_R]$ on each mesh,
$$
z_L^{\text{ov}} = \max(z_L^{\text{new}},z_L^{\text{old}}),\quad
z_R^{\text{ov}} = \min(z_R^{\text{new}},z_R^{\text{old}}),
$$
and, when $z_R^{\text{ov}}>z_L^{\text{ov}}$ (overlap exists),
$$
V_{\text{ov}} = (z_R^{\text{ov}} - z_L^{\text{ov}})\,A
\tag{11.26}
$$
All conserved-quantity accumulation is weighted by $V_{\text{ov}}$. The remap
is an instantaneous re-gridding at fixed time, so both meshes must describe
the same physical configuration: the routines **validate** that the two
meshes share one cross-section area $A$ (throwing `ArgumentError` otherwise)
— mixing $A_{\text{old}}$ in the overlap volumes with $A_{\text{new}}$ in the
division volumes would silently scale every concentration by
$A_{\text{old}}/A_{\text{new}}$. $A(t)$ evolution belongs to the ODE (the
$\chi$ swelling law, Chapter 10), never to the remap.

### 11.7.3 First-order conservative remap

**Temperature — energy-conserving.** When `material` is supplied, each overlap
contributes its heat capacity $\rho c_p$, so the new-cell temperature conserves
thermal energy:
$$
\boxed{\;
T_{\text{new},i} = \frac{\sum_{\text{ov}}\rho c_p\,T_{\text{old}}\,V_{\text{ov}}}
{\sum_{\text{ov}}\rho c_p\,V_{\text{ov}}}
\;}
\tag{11.27}
$$
with $\rho c_p = $ `effective_heat_capacity(ξ_old, T_old, material)`
(`remapping.jl:155–159,180`). When `material === nothing` it degrades gracefully to
plain volume weighting,
$$
T_{\text{new},i} = \frac{\sum_{\text{ov}} T_{\text{old}}\,V_{\text{ov}}}{\sum_{\text{ov}} V_{\text{ov}}}
\tag{11.28}
$$
(`remapping.jl:161–163`). Energy-conserving is the recommended default for pyrolysis,
where solid and char differ markedly in $\rho c_p$ (`ale.md §6.3`).

**Concentration — mass-conserving.** Each component accumulates mass $\sum_{\text{ov}}
\xi_{\text{old}}V_{\text{ov}}$, then is divided by the **new** cell volume so that the
new-cell mass equals the summed contributions:
$$
\boxed{\;
\xi_{\text{new},j,i} = \frac{\sum_{\text{ov}} \xi_{\text{old},j}\,V_{\text{ov}}}{V_{\text{new},i}}
\;}
\tag{11.29}
$$
(`remapping.jl:167–168,184–186`). Dividing by $V_{\text{new}}$ (not the overlap
total) is exactly what makes $\xi_{\text{new}}V_{\text{new}} = \sum
\xi_{\text{old}}V_{\text{ov}}$, i.e. mass-exact for partitioning meshes.

> **Implementation note.** `remap_solution!(…, ::FirstOrderRemap; material)`,
> `src/Physics/remapping.jl:97–214`.

### 11.7.4 Second-order remap with gradient reconstruction

The second-order method (`SecondOrderRemap`, default limiter `:minmod`,
`remapping.jl:38–41`) reconstructs the old field linearly about each old cell center
before integrating over overlaps. The reconstructed value at an overlap center
$z_{\text{ov}}$ is
$$
T_{\text{recon}} = T_{\text{old},i} + \nabla T_{\text{old},i}\,(z_{\text{ov}} - z_{\text{old},i}^c),
\tag{11.30}
$$
with cell-centered gradients computed by central, one-sided, or zero rules over the
*actual* center-to-center distances on the (non-uniform, possibly partly inactive)
old mesh (`compute_gradients_for_remap`, `remapping.jl:342–408`). The reconstructed
contributions are accumulated and averaged exactly as in first order, then two
**extrema-clamping limiters** are applied to suppress reconstruction overshoot:
$$
T_{\text{new}} \leftarrow \text{clamp}\big(T_{\text{new}},\,T_{\min}^{\text{old}},\,T_{\max}^{\text{old}}\big),
\tag{11.31}
$$
$$
\xi_{\text{new},j} \leftarrow \min\!\big(\max(\xi_{\text{new},j},0),\,\xi_{\max,j}^{\text{old}}\big),
\tag{11.32}
$$
where the old-mesh extrema come from `extrema_values` (`remapping.jl:414–431`). The
clamps enforce a discrete maximum principle (no new extrema, no negative
concentrations) at the cost of *breaking* exact mass conservation — which the next
step repairs.

### 11.7.5 Limiter mass error

Because the extrema clamps (11.31)–(11.32) perturb the integrals, the
second-order remap carries a mass error of order $10^{-3}$–$10^{-5}$. (If a
caller needs exact second-order conservation, a global per-component rescale
$\xi \leftarrow \xi\,M_{\text{old}}/M_{\text{new}}$ against
`check_remap_conservation`'s ledger is three lines.) First-order remap is
mass-exact and is what the production depletion path effectively uses.

### 11.7.6 Conservation check and fallback

`check_remap_conservation` returns, per component, the relative mass error
$$
\text{mass\_error}_j = \frac{|M_{\text{new},j} - M_{\text{old},j}|}{|M_{\text{old},j}|},
\quad M_j = \sum_i\xi_{j,i}V_i
\tag{11.34}
$$
Expected: $\approx 10^{-14}$ for first-order without limiters;
$10^{-3}$–$10^{-5}$ for second-order with limiting. When a new cell overlaps
*no* active old cell (rare edge case), both methods fall back to copying the
state of the nearest old cell by center-to-center distance
(`find_nearest_cell`).

---

## 11.8 Mesh quality and safety

ALE motion can in principle tangle the mesh ($V_i\le 0$) if a node overtakes its
neighbor, e.g. under an aggressive time step near complete burnout. Three utilities
guard against this.

- **`check_mesh_quality(mesh)`** (`mesh_motion.jl`) returns a named tuple:
  `valid`, `min_volume`, `max_aspect_ratio` (largest ratio of spatially
  adjacent active-cell sizes), and `tangled_cells`. Tangling is detected from
  the **signed** cell width $z_R - z_L$ in node order, not from the cached
  `cell_volumes` — those pass through `abs()`, so an inverted cell would
  otherwise look like a perfectly valid thin cell and the rollback below
  could never fire. Aspect-ratio adjacency is taken in spatial order (sorted
  by cell center), which diverges from index order after h-refinement.

- **`mesh_is_valid(mesh)`** (`mesh_motion.jl`) is the cheap boolean: `false`
  if any active cell has non-positive signed width.

- **`move_mesh_safe!(mesh, w, dt; max_aspect_ratio=10.0)`**
  (`mesh_motion.jl:229–256`) snapshots nodes/volumes/time, performs the motion,
  checks quality, and **rolls back** (restoring nodes, volumes, time, and
  recomputing geometry) returning `false` if the result is invalid or the aspect
  ratio exceeds the cap; otherwise returns `true`. This is the transactional motion
  used where a step must be all-or-nothing.

Within the implicit solve the integrator's own step-size control (Chapter 15) is the
first line of defense; these utilities back the explicit mesh-management path and
the AMR/depletion machinery (Chapter 14). The minimum-thickness floor `min_thickness`
carried on the `Workspace` (`problem_residual.md`) likewise prevents surface cells
from collapsing to zero during deep burns.

---

## 11.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

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
pressure gradient gas-transport formulation (Chapters 7, 9;
`discretization.md §3,§4`), so on a constant cross-section with $w=0$ the two
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

## 11.10 Limitations and open issues

The following limitations are documented in the source and KB notes
(`ale.md §11`, `volume_change.md` Open Questions, `discretization.md §15`,
`jacobian.md`):

1. **Dimensionality.** The framework is implemented only for 1D ($D=1$). The mesh
   data structures (`MeshTopology{D,…}`, `MeshGeometry{D}`) are parameterized by
   dimension, and the GCL is written as a face-flux sum, but a 2D/3D extension would
   require vector node positions, per-face normals and variable areas, the full
   divergence-theorem GCL, and surface-integral fluxes (`geometry.md §13`,
   `ale.md §11`). The current code assumes a spatially uniform (though time-varying)
   cross-section.

2. **Separation of axial and radial motion.** Variable cross-section is handled by
   splitting axial ALE motion (with $\theta_A$ subtracted) from a separate radial
   cross-section update via $\mathrm{law}(\bar\chi)$, requiring careful sequencing
   (compute axial velocity, then update area). A unified space-time ALE treatment of
   both channels simultaneously is not implemented (`ale.md §11`, item 1).

3. **Order of the explicit mesh advance.** The standalone `move_mesh!` advances nodes
   with first-order forward Euler (11.16). Higher-order *explicit* mesh motion would
   need substeps or a space-time formulation; within the implicit ODE solve, $z_i$ is
   integrated to the solver's order, so this limitation affects only the diagnostic
   mesh-management path (`ale.md §11`, item 3).

4. **Remap limiter placement and order selection.** The second-order remap applies
   TVD/extrema limiters *after* volume-weighted averaging rather than *during*
   reconstruction, which can permit small violations of the discrete maximum
   principle and a $10^{-3}$–$10^{-5}$-level mass error (§11.7.5). The
   choice between first- and second-order remap is left to the caller; there
   is no automatic, error-estimate-driven order selection
   (`ale.md §11`, items 4–5).

5. **Density evaluation in the strain rate.** The volume-change rate (11.6) uses
   $\rho_j(T)$ at the instantaneous cell temperature; whether an implicit or
   time-averaged density would improve consistency near steep thermal transients is
   an open modeling question, as is the extension of the volume-change accounting to
   multi-reactant stoichiometry (`volume_change.md` open questions 1–2).

6. **Jacobian profile trade-offs.** The structured Jacobian (Chapter 15) represents
   the cumulative mesh velocity as a `PrefixSumCoupling`, the column-global $\theta_A$
   back-reaction as a `Rank1Coupling`, and the $z$-columns as a dedicated geometry
   AD block. Under the `:ale_reduced` profile (used by `ApproxSolve`) the cumulative
   widening, geometry block, and $\chi$-column rank-1 are *dropped*, trading exact
   coupling for speed; the WithChi$\times$Darcy interaction is exercised mainly with
   Fickian transport (`jacobian.md`, open issues). These are solver-accuracy/cost
   choices rather than physical limitations of the ALE formulation itself.

7. **Solid-Lagrangian capability boundary: no general (non-material) mesh
   motion.** As stated in §11.1.2, the formulation identifies mesh velocity
   with solid material velocity and gives no component a discrete advection
   operator. Consequently mesh motion decoupled from the solid — remeshing,
   mesh smoothing/redistribution, prescribed node motion, AMR performed
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
  Finite-Volume Discretization**; the remap operators support **Chapter 14, Adaptive
  Mesh Refinement**.
- The structured Jacobian couplings for ALE mesh velocity, advection, dilation, and
  the geometry block, and the time integration of the coupled $(T,\xi,z,\chi)$ system
  are in **Chapter 15, Temporal Integration & Structured Jacobian**.
- GCL tracking as a conservation diagnostic is in **Chapter 16, Conservation
  Diagnostics**.
