# 11. The Arbitrary Lagrangian–Eulerian (ALE) Framework

Pyrolysis is, intrinsically, a *moving-domain* problem: the condensed phase recedes
as it gasifies, intumescent chars swell, and the position of the exposed surface is
itself an unknown. This chapter develops the Arbitrary Lagrangian–Eulerian (ALE)
formulation that Pyrolysis.jl uses to track that motion without solving a
moving-boundary PDE. The node positions of the mesh are differential unknowns of the
state, $dz_i/dt = w_i$, and the node velocities are algebraic unknowns fixed by one
row per cell. The mesh velocity is the velocity of the solid plus a correction,

$$
\mathbf w = \mathbf v_s + \mathbf c ,
$$

where the solid velocity $\mathbf v_s$ satisfies the volume balance of every cell
(the *balance rows*, §11.3) and the correction $\mathbf c$ is zero unless the mesh is
asked to move relative to the solid: by the equidistribution correction that carries a
vanishing layer to burnout in one event (§11.7), or by a prescribed mesh velocity for
verification (§11.7.8). Wherever $\mathbf c \ne 0$ every inventory crosses the faces
of the mesh at the relative velocity $-\mathbf c$ (the *relative fluxes*, §11.6). The
per-cell component masses are the states, so mass is conserved exactly, and the
mesh-frame concentration change that a shrinking cell produces (the dilation) follows
from $\xi_{j,i} = m_{j,i}/V_i$ without a term of its own (§11.3.4). The cell volumes are
recomputed from the node positions at every residual evaluation, so the discrete
Geometric Conservation Law (GCL) holds identically in the semi-discrete system
(§11.5). Throughout, the relations are tied to the implementation in
`src/Closures/balance.jl` (the balance rows), `src/Closures/correction.jl` (the
correction and its layer time), `src/Closures/lateral.jl` (the variable cross
section), `src/Kernels/relative_flux.jl` (the relative fluxes),
`src/Physics/volume_change.jl` (the strain rate), `src/Meshes/geometry.jl` (the
metrics) and `src/Assembly/residual.jl` (the residual).

The coordinate convention is fixed by the nomenclature (Chapter 2): $z$ is the
through-thickness axial coordinate, $z=0$ is the back/substrate (tag `:bottom`),
$z=L$ is the exposed surface (tag `:top`), heat enters from the top, and the
material shrinks downward. Cell index $i$ increases from bottom to top; on a line
cell $i$ lies between nodes $i$ and $i+1$, and face $f$ sits at node $f$, so the
*left* node and face of a cell are at lower $z$, the *right* at higher $z$. The flux
sign convention is **positive flux = transport in the $+z$ direction**. In this
chapter $i$ denotes the **cell index** (and, on a line, the node index of the node
below cell $i$) and $r$ the reaction index (overload note I1).

---

## 11.1 Motivation: tracking a receding/expanding surface

![Aligned before and after strips for fixed, Lagrangian and corrected ALE meshes, with distinct mesh nodes and material markers.](../../assets/figures/03_mesh_motion.svg)

**Figure 3. Three mesh descriptions (schematic).** All panels share the same initial length and pinned back. A fixed mesh represents conversion as condensed-density loss in unchanged volumes, shown by lighter shading. With `c = 0`, nodes follow the solid. Corrected ALE lets interior nodes move relative to the same material motion: $w=v_s+c$ and the inventory transport velocity relative to the mesh is $v_s-w=-c$. The hollow material markers are separate from the filled mesh nodes; their redistribution is illustrative, not a computed trajectory.

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
   surface collapse one after another, each at its own time, and the time step is
   driven by the smallest cell. For deep burns this is the dominant cost.

3. **Arbitrary Lagrangian–Eulerian (ALE).** Let the mesh move with a velocity
   $w(z,t)$ that is *neither* the material velocity (pure Lagrangian) *nor* zero
   (pure Eulerian) in general, but is chosen to keep the mesh well-conditioned while
   still tracking the receding surface. Pyrolysis.jl adopts this approach.

The defining freedom of ALE is that the mesh velocity $w$ is a *modeling choice*,
constrained only by boundary conditions and the requirement that the mesh follow the
domain boundaries. Pyrolysis.jl splits the choice in two. The **solid velocity**
$v_s$ is physics: it is the velocity at which the skeleton moves when its cells
change volume, fixed by the volume balance of every cell with the back node pinned
(§11.3). The **correction** $c = w - v_s$ is numerics: zero by default, so the mesh
rides the skeleton and is Lagrangian for the condensed phase; nonzero only where the
equidistribution correction is active in a layer that is about to vanish (§11.7) or
where a verification run prescribes the mesh velocity (§11.7.8). The surface node
always moves with the solid, so the mesh tracks the receding (or expanding) surface
exactly.

### 11.1.2 Every component rides the skeleton

The frame of each component follows from the frame of its transport closure:

- **Solid and liquid components** are bound to the deforming matrix skeleton: their
  material velocity is $\mathbf v_s$. They have no transport flux of their own.

- **Gas components** percolate through the porous skeleton. Darcy's law is a drag
  balance between the gas and the pore walls. The lumped and dusty-gas wall
  laws use that same reference frame, so the flux the
  closures return, $\mathbf J_g = \xi_g(\mathbf v_{\text{gas}} - \mathbf v_s)$,
  is the gas flux **relative to the skeleton** (§3.2.3, Chapter 9). The
  laboratory-frame gas flux is $\xi_g\mathbf v_s + \mathbf J_g$: the pore gas is
  carried with the pores it occupies and slips through them only by $\mathbf J_g$.

On a mesh that moves with the skeleton ($\mathbf c = 0$) every component is therefore
**Lagrangian** apart from $\mathbf J_g$: the mass of every component in a cell,
gas included, changes only by reaction and by $\mathbf J_g$ through the cell's
faces, and the faces are crossed by nothing else. Where the mesh moves relative to
the skeleton ($\mathbf c \ne 0$), the skeleton and the pore gas it holds cross the
faces together at the relative velocity $-\mathbf c$, the gas then slipping through
the skeleton by $\mathbf J_g$ as before (§11.2.3, §11.6).

ThermaKin writes Eq. 15 in stationary coordinates, including a deformation
term for every component. Transforming that balance to a skeleton-following
frame removes the corresponding mesh-relative term, provided the velocity and
density closures are matched (§3.6.2). Gpyro Eq. 58 instead uses control
volumes that deform with the solid, with $\dot m''$ crossing their faces. Treating the gas as fixed in the laboratory instead would count the
skeleton's motion twice: once inside the frame of $\mathbf J_g$ and once as a sweep
of the gas by the moving mesh.

> **Formulation-level fact: the mesh velocity and the solid velocity are distinct
> unknowns.** The solid velocity is $\mathbf w - \mathbf c$, and the balance rows
> (§11.3.2) hold *it*, not the mesh velocity, to the cells' volume balance. Any mesh
> motion relative to the solid — the equidistribution correction, a prescribed node
> motion for a manufactured solution — is represented by $\mathbf c$, and its effect
> on the fields is carried by one operator, the relative fluxes of §11.6, which move
> every inventory (each component's mass, its sensible enthalpy and the
> conversion-history pair) across the faces at the relative velocity. With
> $\mathbf c = 0$ that operator vanishes identically and is not assembled: the
> mesh-frame rate of every field is then its material rate. No additional
> skeleton-advection term is needed; gas transport relative to the skeleton remains.

---

## 11.2 The ALE balances on a moving cell

### 11.2.1 Kinematics of a moving sample point and of a moving cell

Let $f(z,t)$ be any field (temperature, concentration). Introduce a moving sample
point whose position $z=Z(\chi,t)$ is labeled by a material-like coordinate $\chi$
that is constant along the mesh trajectory, with mesh velocity
$$
w \;=\; \left.\frac{\partial Z}{\partial t}\right|_{\chi}.
$$
(Here and only here $\chi$ is a mesh label; elsewhere $\chi$ is the pyrolysis
progress, overload note A1.) The value of $f$ seen by an observer riding on that
point is $\hat f(\chi,t) = f(Z(\chi,t),t)$. Differentiating with the chain rule
gives the **ALE time-derivative identity**
$$
\boxed{\;
\left.\frac{\partial f}{\partial t}\right|_{\chi}
=
\left.\frac{\partial f}{\partial t}\right|_{z}
+ \mathbf w\cdot\nabla f
\;}
\tag{11.1}
$$
Here $\partial/\partial t|_{\chi}$ is the **ALE (mesh-point)** rate — the rate at
which $f$ changes for a point carried by the mesh — and $\partial/\partial t|_{z}$ is
the **Eulerian (fixed-point)** rate that the conservation laws of Chapter 3 are written
in. In 1D, $\mathbf w\cdot\nabla f = w\,\partial f/\partial z$.

A finite-volume cell $\Omega_i(t)$ whose boundary moves with the mesh velocity obeys
the integral form of (11.1), the Reynolds transport theorem for a moving control
volume:
$$
\frac{d}{dt}\int_{\Omega_i(t)} f\,dV
= \int_{\Omega_i(t)} \left.\frac{\partial f}{\partial t}\right|_{z} dV
+ \oint_{\partial\Omega_i(t)} f\,\mathbf w\cdot\hat{\mathbf n}\,dA .
\tag{11.2}
$$
Pyrolysis.jl integrates cell inventories (masses, and in the inventory form the
sensible enthalpy), so (11.2) is the form every row of the residual is derived from:
substituting the Eulerian conservation law for $\partial f/\partial t|_z$ turns the
volume integral of a flux divergence into face fluxes, and the last term combines
with the laboratory-frame face flux into a flux *relative to the moving face*.

### 11.2.2 The energy balance — the relative-velocity form

Let $H_i = \sum_j m_{j,i}\,h_j(T_i)$ be the sensible-enthalpy inventory of cell $i$,
with $h_j(T) = \int_{T_{\text{ref}}}^{T} c_{p,j}\,dT'$ on a common datum. The
laboratory-frame enthalpy carried by component $j$ is
$(\xi_j\mathbf v_s + \mathbf J_j)\,h_j$, and (11.2) with $\mathbf w = \mathbf v_s +
\mathbf c$ gives the cell balance
$$
\frac{dH_i}{dt}
= -\oint_{\partial\Omega_i}\Big[\mathbf q + \sum_j \big(\mathbf J_j - \xi_j\,\mathbf c\big)\,h_j\Big]\cdot\hat{\mathbf n}\,dA
+ \int_{\Omega_i}\Big(Q_{\text{rxn}} + Q_{\text{rad}} + \sum_j S_j\,h_j\Big)\,dV ,
\tag{11.3}
$$
with $\mathbf q$ the conductive flux, $Q_{\text{rxn}} = -\sum_r h_r r_r$ the reaction
heat at the local temperature, $Q_{\text{rad}}$ the absorbed radiation and $S_j$ the
reaction mass source of component $j$ (the last term converts the heat of reaction,
quoted at the local temperature, to the sensible-enthalpy datum; Chapter 7). This is
the row of the **inventory form** (`energy_state = :enthalpy`). The **temperature
form** (the default) subtracts $\sum_j h_j(T_i)\,dm_{j,i}/dt$, using the mass balance
(11.5) below, and divides by the cell's heat capacity
$(\rho c_p^{\text{eff}})_i V_i = \sum_j m_{j,i}\,c_{p,j}(T_i)$. Its continuous
counterpart, per unit volume, is
$$
\boxed{\;
\rho c_p^{\text{eff}}\left.\frac{\partial T}{\partial t}\right|_{\chi}
= -\nabla\cdot\mathbf q - \sum_g c_{p,g}\,\mathbf J_g\cdot\nabla T + Q_{\text{rxn}} + Q_{\text{rad}}
\;+\; \rho c_p^{\text{eff}}\,(\mathbf w - \mathbf v_s)\cdot\nabla T
\;}
\tag{11.4}
$$
with $\mathbf w - \mathbf v_s = \mathbf c$. The advecting velocity is the velocity of
the mesh *relative to the material that carries the thermal mass*, and on a mesh that
rides the skeleton it vanishes identically: **with $\mathbf c = 0$ there is no
temperature advection at all**, and the mesh-frame rate of the matrix-carried
temperature *is* its material rate. (Temperature is intensive, so there is also no
dilation contribution: a static temperature field is unchanged by uniform stretching
of the mesh.)

Where $\mathbf c \ne 0$ the last term of (11.4) is not assembled pointwise. It is
the temperature-form reading of the relative enthalpy flux in (11.3): every face
carries $\sum_j F_{j,f}\,h_j(\tilde T_f)$, the sensible enthalpy of the same relative
mass fluxes $F_{j,f}$ that the mass rows carry, at the face temperature $\tilde T_f$
reconstructed from the same cells (§11.6.4). The energy row is then conservative in
enthalpy, and in the temperature form the difference
$h_j(\tilde T_f) - h_j(T_i)$ is the discrete advection of the temperature relative
to the mesh.

> **Why the distinction matters.** Applying the *full* chain-rule term
> $+w\,\partial T/\partial z$ of (11.1) to the energy equation while the mass rows
> treat the solid as mesh-bound would advect the temperature field relative to the
> very material carrying it: a spurious energy sink
> $\int \rho c_p\, w\,\partial_z T\,dz \approx \dot m''\,c_s\,(T_s - \bar T_{\text{form}})$
> — a few percent of the effective heat of gasification in every ALE run with volume
> change and thermal gradients, and a systematic low-MLR bias against
> a consistently discretized material-following energy balance. This is a
> formulation check, not a measured bias against every reference code. The Landau ablation case V-CPL-02 verifies the traveling-wave identity
> that would expose such a term (Chapter 18).

**Sign of the relative energy advection.** With $c>0$ (mesh sliding in $+z$ relative
to the material) and $\partial T/\partial z>0$, the term $c\,\partial T/\partial z>0$:
the node moves into warmer material, so the temperature it carries increases — the
physically correct convective reading of the chain rule.

### 11.2.3 The mass balance of every component

For component $j$ with laboratory-frame flux $\xi_j\mathbf v_s + \mathbf J_j$
($\mathbf J_j \equiv 0$ for condensed components, $\mathbf J_g$ the gas flux
for gases), the Eulerian balance is
$\partial\xi_j/\partial t|_z + \nabla\cdot(\xi_j\mathbf v_s + \mathbf J_j) = S_j$, and
(11.2) gives, for the mass $m_{j,i} = \int_{\Omega_i}\xi_j\,dV$,
$$
\boxed{\;
\frac{dm_{j,i}}{dt}
= S_{j,i}\,V_i - \oint_{\partial\Omega_i}\big(\mathbf J_j + \xi_j(\mathbf v_s - \mathbf w)\big)\cdot\hat{\mathbf n}\,dA
= S_{j,i}\,V_i - \oint_{\partial\Omega_i}\mathbf J_j\cdot\hat{\mathbf n}\,dA
+ \oint_{\partial\Omega_i}\xi_j\,\mathbf c\cdot\hat{\mathbf n}\,dA .
\;}
\tag{11.5}
$$
Discretized over the faces of the cell, with $s_{i,f} = +1$ when cell $i$ owns face
$f$ (the face normal points out of it) and $-1$ otherwise,
$$
\frac{dm_{j,i}}{dt} = S_{j,i}V_i - \sum_{f\in\partial i} s_{i,f}\,A_f\,
\big(J_{j,f} + v_{\text{rel},f}\,\xi_{j,f}\big),
\qquad v_{\text{rel},f} = (\mathbf v_s - \mathbf w)_f\cdot\hat{\mathbf n}_f = -\mathbf c_f\cdot\hat{\mathbf n}_f ,
\tag{11.6}
$$
which is exactly the mass row the assembler writes (Chapter 13): the reaction source
plus the scattered face fluxes of the gas-transport kernel ($J_{j,f}$) and of the
relative-flux kernel ($v_{\text{rel},f}\,\xi_{j,f}$, §11.6). On a mesh that rides the
skeleton the relative term vanishes for **every** component, and the mass of every
component in a cell changes only by reaction and by $\mathbf J_j$ through its faces,
so the column total changes only through the boundary faces, exactly as on a fixed
mesh. In a sealed collapsing column the pore gas is compressed cell by cell
($\xi_g = m_g/V \uparrow$ as $V \downarrow$, V-CPL-03 P4), and under gas
transport the pressure closure turns that compression into a venting response.

**Why not a face-sweep flux for the gas.** Treating the gas as fixed in the
laboratory while the mesh sweeps past it would add $\xi_g\,\mathbf w\cdot\hat{\mathbf n}$
to the gas's face flux. With $\mathbf J_g$ defined relative to the skeleton that term
counts the skeleton's motion a second time: in a uniformly collapsing sealed column it
leaves every interior cell's gas untouched and piles the whole compression into the
top cell, and at an open regressing surface it releases pore gas at $\xi_g\,w$ that
the diffusive and Darcy fluxes already account for.

### 11.2.4 Why the gas is in $\rho c_p^{\text{eff}}$ (cross-reference)

The advection term of (11.4) is the mesh-frame transport of *sensible* enthalpy, and
it vanishes on a mesh that rides the skeleton because the thermal mass rides with it.
That thermal mass is the mixture heat capacity $\rho c_p^{\text{eff}} =
\sum_j\xi_jc_{p,j}$ over every component (Chapter 7): the pore gas is carried with the
pores, so its sensible heat moves with the mesh like the matrix's and needs no
separate advection term either. The complementary gas enthalpy carried by the *flux
relative to the skeleton* is assembled at the faces, in flux form: every interior face
carries $\sum_g J_{g,f}\,h_g(T_f)$ with the same $J_{g,f}$ the mass rows carry, and its
temperature-form reading is the $-\sum_g c_{p,g}\,\mathbf J_g\cdot\nabla T$ term of
(11.4) (Chapter 7). This corresponds to the 1D transport terms of ThermaKin Eq. 17, including
its deformation term, which is the $\mathbf v_s\cdot\nabla T$ piece the mesh frame
absorbs. At a boundary face the outward mass flux carries its sensible enthalpy at the
surface temperature, while incoming species carry their specified reservoir
enthalpy. Wall-film heat and carried enthalpy are converted to the same energy
control surface before assembly (§9.9); no additional outgoing-gas sink is added.

---

## 11.3 The solid velocity: balance rows

The ALE freedom is fixed by *choosing* $w(z,t)$. The solid part of it is not a
choice: it is the velocity at which the skeleton moves when its cells change volume,
with the back node pinned. This section derives the strain rate, the discrete volume
balance that determines the solid velocity, the lateral channel of a variable cross
section, and why the masses as states need no dilation term.

### 11.3.1 The volumetric strain rate $\theta$

For a fixed-mass control volume of the deforming matrix, the relative volume change
rate (volumetric strain rate, dilation) is obtained from the mixture continuity
relation by summing the volume freed or created per unit mass change of each
matrix-contributing component:
$$
\boxed{\;
\theta \;\equiv\; \frac{1}{V}\frac{dV}{dt}
\;=\; \sum_{j=1}^{N_C} \gamma_j\,\frac{1}{\rho_j}\,S_j
\;}
\tag{11.7}
$$
(`volume_change_rate`, `src/Physics/volume_change.jl`). Here:

- $\gamma_j$ is the **swelling factor** of component $j$ (component field
  `swelling_factor`): $\gamma=1$ for solids by default (when a solid depletes, the
  matrix loses that volume), $\gamma=0$ for gases by default (escaping gas does not
  change the skeleton volume), and $\gamma\in[0,1]$ for partial contributors. Only
  components with $\gamma_j>0$ enter the sum.
- $\rho_j(T)$ is the **bulk** (pure-phase) density of component $j$, evaluated at the
  cell's trial temperature through the property function `comp.density(T)`.
- $S_j$ is the **reaction mass source** of component $j$ per unit volume. The
  closure reads it from the reaction kernel's cell sources, $S_{j,i} =
  (\text{reaction mass source})_{j,i}/V_i$, after every kernel has run and before any
  external source is added (`closure_rows!`), so the strain rate is the reactions'
  alone: a manufactured source of a verification run does not move the mesh. It is
  negative for a consumed solid and positive for a produced char or gas.

**Physical reading and sign.** When a solid of density $\rho_j$ is consumed
($S_j<0$, $\gamma_j=1$), the term $\gamma_j S_j/\rho_j<0$, so $\theta<0$: the cell
shrinks. A net char-producing or swelling reaction with the appropriate $\gamma$ can
yield $\theta>0$ (expansion). Note that (11.7) is the *swelling-factor-weighted*
form; the loose statement "$\theta = -(\sum_j S_j)/\rho_{\text{total}}$" is only
correct in the special case of a single solid with $\gamma=1$ and constant $\rho$.
The swelling-factor semantics and the endpoints they imply are developed in
Chapter 10.

### 11.3.2 The discrete volume balance and the balance rows

The solid velocity must reproduce the volume change of every cell. In the continuous
planar picture it is the integral of the *axial* strain rate from the pinned back
face,
$$
v_s(z,t) = \int_0^{z}\big(\theta(z',t) - \theta_A(t)\big)\,dz',
\qquad v_s(0,t)=0,
\tag{11.8}
$$
where $\theta_A$ is the **column-global lateral strain rate** of a variable cross
section, subtracted so the lateral part of the volume change is not counted twice
(§11.3.3); for a constant cross section $\theta_A\equiv 0$. The axial strain rate is
$\theta_L = \theta - \theta_A$.

The discrete statement is the volume balance of each cell at the solid velocity: the
rate at which the faces of cell $i$, moving with the solid, sweep volume equals the
cell's axial volume change,
$$
\sum_{f\in\partial i} s_{i,f}\,A_f\,(\mathbf v_s\cdot\hat{\mathbf n})_f = (\theta_i - \theta_A)\,V_i .
$$
On a line, with face areas $A_f$ that follow the coordinate system (Chapter 13) and
$\mathbf v_s = \mathbf w - \mathbf c$, the algebraic rows of the node velocities are
$$
\boxed{\;
0 = w_1, \qquad
0 = A_{i+1}\,(w - c)_{i+1} - A_i\,(w - c)_i - (\theta_i - \theta_A)\,V_i
\quad (i = 1,\dots,n)
\;}
\tag{11.9}
$$
(`BalanceRowClosure`, `src/Closures/balance.jl`): the pinned back node is at rest,
and cell $i$'s balance sits on the velocity row of its node farther from the pinned
node, node $i+1$. On a planar line, where every face has the area $A$ of the cross
section, the rows reduce to the recursion
$$
v_{s,1} = 0,\qquad v_{s,i+1} = v_{s,i} + (\theta_i - \theta_A)\,\Delta z_i ,
\tag{11.10}
$$
the discrete form of (11.8): the solid velocity at node $i+1$ is the exclusive prefix
sum $\sum_{k\le i}\theta_{L,k}\,\Delta z_k$ of the axial volume change below it, and
the surface node moves at the total axial volume-change rate of the column, exactly
tracking the receding (or expanding) surface. In cylindrical and spherical coordinates
the face areas grow with the radius and (11.9) carries them: a uniform strain rate
gives $v_s = (\theta - \theta_A)\,z$ on a planar line, $\theta z/2$ on a cylinder and
$\theta z/3$ on a sphere (with $z$ the radius), which the tests assert at every node,
with and without a correction.

**Why rows, not a sweep.** The prefix sum (11.10) couples the velocity of every node to
the strain rate of every cell beneath it, a dense lower-triangular block in the
Jacobian. Each row of (11.9) reads one cell and its two nodes (and, under a
correction, what the corrections of those two nodes read, §11.7), so the Jacobian
stays banded; the integrator solves the bidiagonal rows together with the rest of the
system. The rows are algebraic (zero rows of the mass matrix) and determine the
velocities whenever every cell's coefficient on its row node, the area of its face
there, is nonzero (`rank_check`); the system is then a semi-explicit DAE of index one
(Chapter 15). The coefficient of the centre node of a sphere is zero, which is why the
pinned node sits there: the closure requires exactly one pinned node, at an end of the
line, and a line mesh pins its first node.

**Consistent velocities.** The integrator starts from a state whose algebraic rows
hold. `consistent_velocities!` sweeps (11.9) once from the pinned node: the
corrections read positions and masses only, so one pass is exact. `setup` calls it for
the initial state (through `consistent_state!`, which also sets the face intensities
and the cross-section rows), and every deactivation event calls it again after its
transfer, a dead cell's row node taking zero velocity (Chapter 14). A velocity row left
above round-off by the sweep is an error at setup; at an event it ends the solve with
return code `:Failure` and a warning that no consistent node velocities exist.

**Substrate boundary condition.** Pinning $w_1=0$ enforces zero normal mesh velocity
at the impermeable substrate: the back node never moves, all motion is taken up above
it, and the material "stacks up" toward the surface as it is consumed. This is the
natural choice for a slab on a non-receding backing.

### 11.3.3 The lateral strain rate $\theta_A$ (variable cross section)

When the material carries a *lateral shrinkage law* (Chapter 10), the column cross
section evolves as $A(t) = A_0\,\mathrm{law}(\bar\chi)$ with $\bar\chi$ the column's
conversion. Part of the total volume change is then lateral and must be removed from
the axial balance (11.9) so it is not counted twice. The area and its lateral strain
rate are two algebraic unknowns of the state, the *cross-section rows*
(`LateralShrinkage`, `src/Closures/lateral.jl`):
$$
\boxed{\;
0 = A - A_0\,\mathrm{law}(\bar\chi_c),
\qquad
0 = \theta_A - \frac{d}{dt}\log\mathrm{law}(\bar\chi_c)
\;}
\tag{11.11}
$$
with
$$
\bar\chi = \frac{\sum_i G_i}{\sum_i m_{\text{dry},i}},
$$
$G_i$ the cumulative dry pyrolysis gas cell $i$ has produced (sourced by the dry-gas
reaction rates) and $m_{\text{dry},i}$ its initial dry solid mass (no source): the
*conversion-history pair*, carried per cell as two differential unknowns. The time
derivative in the second row is taken along the column sums of the history rows'
rates, $(\sum_i\dot G_i, \sum_i\dot m_{\text{dry},i})$, by forward-mode automatic
differentiation of $\log\mathrm{law}$, so on the solution the area's own rate is
exactly $dA/dt = \theta_A A$. $\bar\chi_c$ is $\bar\chi$ smoothly clamped into
$[-10^{-3}, 1 + 10^{-3}]$ (the identity on $[0, 1]$) and the denominator is smoothly
floored at $10^{-8}$ of the initial dry solid mass, so both rows stay finite with
finite derivatives on any trial state. Both rows read every cell's history pair, the
lateral one also its rates: they are the system's two dense rows. `setup` refuses a
law that is not positive and finite on the clamp band or does not equal one at zero
conversion, and a lateral shrinkage law needs `use_ale = true` and planar
coordinates. In (11.9) $\theta_A$ enters as the scalar subtracted from every cell's
strain rate; the law itself and the history pair are developed in Chapter 10.

### 11.3.4 Masses as states: the dilation is implicit

A formulation that integrates concentrations on a moving mesh needs a dilation term:
in the mesh frame a shrinking cell must concentrate the mass it carries, so
$d\xi_{j,i}/dt$ gains $-\xi_{j,i}\theta_i$. Pyrolysis.jl integrates the **masses**
$m_{j,i}$, whose rows (11.6) have no geometric term at all, and the state mapping
recovers the concentration from the current volume, $\xi_{j,i} = m_{j,i}/V_i$
(Chapter 13). The dilation then holds as an identity: with $\mathbf c = 0$ and no
transport flux,
$$
\frac{d\xi_{j,i}}{dt} = \frac{1}{V_i}\frac{dm_{j,i}}{dt} - \frac{m_{j,i}}{V_i^2}\frac{dV_i}{dt}
= S_{j,i} - \xi_{j,i}\,\theta_i ,
$$
using $dV_i/dt = \theta_i V_i$ (§11.5.2). Every component, the pore gas included, is
concentrated by a shrinking cell and diluted by an expanding one without a term that
could drift from the volume change, and the mass of each component is conserved to
round-off whatever the mesh does.

### 11.3.5 Meshes of more than one dimension (residual level)

The mesh closure is defined on meshes of quadrilaterals and hexahedra as well; these
are building blocks at the residual level, not a supported path for full solves
(`setup` refuses gas transport, in-depth radiation, the variable cross section, the
equidistribution correction and deactivation on them). Two closures exist there:

- **`PotentialClosure`**, one multiplier $\phi_i$ per cell: the node velocities of
  minimum discrete norm $\sum_n W_n|\mathbf w_n|^2$, with $W_n$ the node's lumped
  volume over the live cells, that satisfy every live cell's volume balance
  $\sum_f s_{i,f}\sum_q \mathbf s_{q,f}\cdot\mathbf w_q = (\theta_i - \theta_A)V_i$, with
  $\mathbf s_{q,f}$ the sweep vectors of the face's nodes (Chapter 13). Stationarity
  gives $\mathbf w = W^{-1}D^{\mathsf T}\phi$, a local function of the multipliers of
  the cells around each node, written before the kernels run; each cell's balance is
  the algebraic row of its multiplier. The global system $D\,W^{-1}D^{\mathsf T}\phi =
  (\theta - \theta_A)V$ is solved, by one sparse Cholesky factorization, only for the
  consistent initial state (`consistent_multipliers`). Pinned nodes are at rest, and a
  `slip` wall fixes one velocity component on the faces of a tag.
- **`BalanceRowClosure` in column mode**: a single column of cells along the last
  axis, the lateral components and the pinned nodes at rest, cell $k$'s balance on
  the axial velocity of one node of the face it adds and that face's other nodes tied
  to it, so the column's levels stay level and the rows are the line's.

On such meshes the cells can also carry the deformation gradient
$\mathbf F$, with rows $d\mathbf F_i/dt = (\nabla\mathbf v_s)_i\,\mathbf F_i$ and the
cell-mean solid velocity gradient from the faces' sweep vectors, and the relative
fluxes are first-order upwind (§11.6.5).

---

## 11.4 Node positions and the geometry from the state

On a moving mesh (`use_ale = true`) the node positions $z_i$ are differential unknowns
of the state, interleaved with the cell unknowns (Chapter 13). Their rows are simply
$$
\boxed{\;
\frac{dz_i}{dt} = w_i
\;}
\tag{11.12}
$$
written by the row assembly from the node-velocity unknowns, which are themselves
fixed by (11.9). The time integrator (Chapter 15) advances the positions together with
the masses, the energy states and the algebraic rows as one coupled stiff system; the
colored Jacobian carries every coupling (Chapter 15). Without `use_ale` every node is
pinned: the state carries neither positions nor velocities, and the geometry is
computed once.

### 11.4.1 Geometric volume recompute

At every residual evaluation the state mapping (`derive!`) reads the node positions
from the state and recomputes every metric **from geometry**, not from a flux formula
(`geometry!`, `src/Meshes/geometry.jl`):
$$
\boxed{\;
V_i = \begin{cases}
A\,(z_{i+1} - z_i) & \text{planar},\\
\pi A\,(z_{i+1}^2 - z_i^2) & \text{cylindrical},\\
A\,(z_{i+1}^3 - z_i^3)/3 & \text{spherical},
\end{cases}
\qquad
z_i^c = \tfrac12\big(z_i + z_{i+1}\big),
\;}
\tag{11.13}
$$
with $A$ the cross-section measure (the planar face area, the axial length of a
cylinder, the solid angle of a sphere), under a variable cross section the area
unknown of (11.11). Face areas follow the same coordinate system (Chapter 13). The
volume is signed: an inverted trial mesh gives a non-positive volume, which the state
mapping floors smoothly for evaluation and the admissibility check reads unfloored
(§11.8). The metrics carry the element type of the positions, so automatic
differentiation propagates through the geometry; the residual is a pure function of
the state that never mutates the mesh, whose stored positions are only the reference
positions at which the initial conditions are evaluated.

**Why geometric recompute rather than an evolved volume?** Recomputing the volume
from the current node positions and the current cross section bakes a simultaneous
change of $A(t)$ into every volume automatically; carrying a volume forward by a flux
formula across an $A(t)$ change would keep the area of an earlier time in it. There is one
volume per cell, so there is no second copy that could drift from the positions. As
shown in §11.5, this path satisfies the GCL as an identity of the semi-discrete
system.

---

## 11.5 The Geometric Conservation Law (GCL)

### 11.5.1 Statement and discrete form

The Geometric Conservation Law is the requirement that a *uniform* field be
preserved exactly under mesh motion — equivalently, that the discrete rate of change
of each cell volume equal the net flux of mesh velocity through its faces. In
continuous form,
$$
\frac{dV_i}{dt} = \oint_{\partial\Omega_i}\mathbf w\cdot\hat{\mathbf n}\,dA .
$$
In 1D, cell $i$ has a left face (outward normal $-\hat z$, velocity $w_i$, area
$A_i$) and a right face (outward normal $+\hat z$, velocity $w_{i+1}$, area
$A_{i+1}$), so
$$
\frac{dV_i}{dt} = A_{i+1}\,w_{i+1} - A_i\,w_i
\tag{11.14}
$$
at a constant cross-section measure. Why this matters: a GCL violation acts as a
*spurious source* in every concentration the properties are evaluated at. If the
discrete volume does not change by exactly the swept volume, a uniform field is no
longer uniform on the moving mesh, and the error accumulates over thousands of steps.

### 11.5.2 Exactness of the geometric path

Pyrolysis.jl never evolves a cell volume: the node positions are state (11.12) and
every residual evaluation recomputes $V_i$ from them (11.13). The derivative of each
shell volume in (11.13) with respect to the position of a face is the area of that
face, in every coordinate system,
$$
\frac{\partial V_i}{\partial z_{i+1}} = A_{i+1},\qquad
\frac{\partial V_i}{\partial z_i} = -A_i ,
$$
which the mesh tests assert against automatic differentiation for planar,
cylindrical and spherical lines. Differentiating along the semi-discrete trajectory,
with $dz_i/dt = w_i$ and, under a variable cross section, $dA/dt = \theta_A A$ from
(11.11) and $\partial V_i/\partial A = V_i/A$,
$$
\frac{dV_i}{dt} = A_{i+1}\,w_{i+1} - A_i\,w_i + \theta_A V_i ,
$$
which is (11.14) plus the lateral area change at constant length: the GCL for a
time-varying cross section holds by construction, with no separate correction. With
the balance rows (11.9),
$$
\boxed{\;
\frac{dV_i}{dt} = \theta_i\,V_i + \Gamma_i ,
\qquad \Gamma_i = A_{i+1}\,c_{i+1} - A_i\,c_i ,
\;}
\tag{11.15}
$$
exactly: a cell changes volume by its strain rate and by the sweep $\Gamma_i$ of the
correction, which is zero on a mesh that rides the skeleton. Combined with the mass
rows (11.6), a uniform concentration $\bar\xi$ with no sources and no transport flux is
preserved: the relative flux of a constant field carries $\bar\xi$ through every face
(§11.6.2), so $dm_{j,i}/dt = \bar\xi\,dV_i/dt$ and $m_{j,i} = \bar\xi V_i$ holds for all
time in the semi-discrete system.

**After time discretization.** On a planar line with a constant cross section the
volume is linear in the node positions, so the relation $m_{j,i} = \bar\xi V_i$ is
linear in the state and survives the integrator's step: a uniform state stays uniform
to round-off at any step. The free-stream tests drive a 12-cell column through a
prescribed compression to 70 % of its thickness and assert the uniform concentration
to $10^{-13}$ relative and the temperature to $10^{-10}$ K at fixed steps of 0.2,
0.1 and 0.05 s, in both energy forms (`test/kernels/test_relative_flux.jl`;
V-MMS-02 is the suite's case). In cylindrical and spherical coordinates the volume is
quadratic or cubic in the radii, and with a variable cross section $V = A\,\Delta z$
is not jointly linear in the state; there the free-stream error is the integrator's
truncation error, which the same tests assert to converge at better than order 3.5 in
the step. Carrying $V_i$ as a state would make free stream exact there at the price of
a drift between $V_i$ and the positions, and is not done.

**What verifies it.** The induced-motion manufactured solution V-MMS-04, whose mesh
velocity is reaction-driven (no prescribed motion), closes the kinetics → strain rate
→ balance rows → geometry chain: it asserts second-order convergence of the solid
concentration on the reacting, compressing mesh (an observed order of at least 1.7
and a mean of at least 1.9 over 16, 32 and 64 cells), and its surface-position error
is $1.6\times10^{-5}$ of the exact recession at 32 cells. The prescribed-motion cases
V-MMS-02 and V-MMS-03 move the mesh over a solid at rest, so they check the geometric
path and the relative fluxes (a uniform field stays uniform, smooth fields converge at
the design order) rather than the balance rows (Chapter 18).

---

## 11.6 Relative fluxes: the ALE advection discretization

One kernel implements mesh-motion advection, `RelativeFluxKernel`
(`src/Kernels/relative_flux.jl`). `setup` installs it whenever the correction is not
`NoCorrection` (the equidistribution correction or a prescribed mesh velocity); under
`NoCorrection` it is absent and nothing crosses the faces of the mesh but the
gas-transport fluxes, and a face whose relative velocity is exactly zero is skipped,
so a bound correction that is zero everywhere carries nothing either. It runs after the thermal boundary kernel, so the surface
energy balance and its blowing models see the gas-transport fluxes alone.

### 11.6.1 The relative velocity and the upwind side

Through each face, along its normal, every component moves with the solid at the
relative velocity
$$
v_{\text{rel},f} = (\mathbf v_s - \mathbf w)\cdot\hat{\mathbf n}_f = -c_f\,\hat n_f ,
\tag{11.16}
$$
with $c_f$ the correction of the face's node (on a line face $f$ sits at node $f$) and
$\hat n_f = \pm1$ its normal, $+1$ at every face but the bottom boundary face. Each
component's mass flux is $F_{j,f} = v_{\text{rel},f}\,\xi_{j,f}$, with the face state
$\xi_{j,f}$ reconstructed from the upwind side. The upwind cell $U$ is the face's
owner when $v_{\text{rel},f} \ge 0$ and its neighbour otherwise; on an interior face
with $\hat n_f = +1$ the owner is the lower cell, so material moving in $+z$ relative
to the mesh is taken from below. The side is read from the *value* of
$v_{\text{rel},f}$ alone: a dual number of zero value compares by its partials, which
differ between derivative seeds, so at $v_{\text{rel},f} = 0$ every derivative pass
takes the $v_{\text{rel},f}\ge 0$ side. The flux $v_{\text{rel},f}\,\xi_{j,f}$ is
continuous where the side switches, with a kink. A face whose relative velocity is
exactly zero carries nothing and is skipped.

### 11.6.2 Second-order reconstruction with the van Albada slope

Each face state is the upwind cell's value plus a slope times its distance to the
face. With $U$ the upwind cell, $D$ the downwind cell, $X$ the live cell beyond $U$
on the far side, $d_{U,f}$ the distance from $U$'s centre to the face and $d_{UD}$,
$d_{XU}$ the centre-to-centre distances,
$$
q_f = q_U + d_{U,f}\;\psi_{\text{VA}}(a, b;\epsilon_q),
\qquad
a = \frac{q_U - q_X}{d_{XU}},\quad b = \frac{q_D - q_U}{d_{UD}},
\tag{11.17}
$$
with the smooth van Albada slope
$$
\boxed{\;
\psi_{\text{VA}}(a, b;\epsilon) = \frac{a\,(b^2 + \epsilon) + b\,(a^2 + \epsilon)}{a^2 + b^2 + 2\epsilon}
\;}
\tag{11.18}
$$
(`van_albada`). It equals $a$ when $a = b$, is zero on a uniform field and at
$a = -b$, stays below $0.21\max(|a|,|b|)$ in magnitude when the two slopes differ in
sign, and is $C^\infty$ for $\epsilon > 0$; a limiter with a kink (minmod, van Leer)
would put that kink into the Jacobian and into finite-difference sensitivities. The
slope scales are, per cell, the square of a tolerance over the cell's reference
thickness: $\epsilon_\xi = (\xi_{\text{tol},U}/\Delta z_U)^2$ for the concentrations,
with $\xi_{\text{tol},i}$ the cell's mass tolerance over its initial volume,
$\epsilon_T = (\texttt{abstol}/\Delta z_U)^2$ for the temperature and
$\epsilon_G = (\texttt{reltol}\cdot\varepsilon_m/\Delta z_U)^2$ for the history values.

Three cases complete the stencil:

- **No live cell beyond $U$.** The gradient from beyond is taken as zero, so the face
  state is $q_U + d_{U,f}\,\psi_{\text{VA}}(0, b;\epsilon) = q_U + d_{U,f}\,b\epsilon/(b^2 +
  2\epsilon)$: beside a resolved gradient ($|b| \gg \sqrt\epsilon$) it is the upwind
  cell's own value, so a boundary cell holding almost none of a component is not
  drained.
- **Boundary face.** The face state is the one-sided linear extrapolation of the
  boundary cell, $q_U + d_{U,f}\,(q_U - q_X)/d_{XU}$, whichever way the mesh crosses
  the face, so a boundary flux has no switch.
- **Distances.** $d_{U,f}$ is the centre-to-face distance along the normal, smoothly
  floored at $10^{-8}$ of its reference value; the centre-to-centre distances are sums
  of two such distances. A reconstruction that took the cell widths for the centre
  distances would leave a zeroth-order error on a mesh whose non-uniformity does not
  vanish under refinement; the tests advect a front at rest under a sliding mesh and
  assert convergence above order 1.6 on uniform cells and above order 1.5 on cells
  alternating between $h$ and $2h$.

A uniform field gives $a = b = 0$ and $q_f = q_U$ on every interior face and every
boundary face (and, for a concentration above its tolerance, (11.19) below returns it
unchanged): the relative flux of a constant field is the constant times the
relative velocity, which is what the free-stream property of §11.5.2 needs.

### 11.6.3 Positivity of the inventory face states

A concentration's face state is bounded by the upwind cell's inventory. With
$p = \mathcal P(\xi_U, \xi_{\text{tol},U})$ the upwind concentration smoothly floored at zero
below its tolerance (the one-sided positive part of Chapter 13: exactly zero at and
below zero, exactly $\xi_U$ above the tolerance) and $\delta = \xi_f - \xi_U$ the
reconstruction's increment from (11.17), the face state used is
$$
\xi_{j,f} = p + \frac{\delta\,p}{\sqrt{p^2 + \delta^2 + \xi_{\text{tol},U}^2}} .
\tag{11.19}
$$
It lies in $[0, 2p]$, differs from the reconstruction by $O(\delta^3/p^2)$ where
$|\delta| \ll p$, and is exactly zero for an empty cell, so no face draws a component
out of a cell that holds none of it, whatever the slopes around it: the inventories'
semi-discrete equations keep them nonnegative. The temperature and the history values
are not inventories and keep the reconstruction (11.17).

### 11.6.4 The energy and history fluxes carried with the mass fluxes

The energy flux of a face gains
$$
F^{E}_{\text{rel},f} = \sum_j F_{j,f}\,h_j(\tilde T_f),
\tag{11.20}
$$
the sensible enthalpy of the same mass fluxes at the face temperature $\tilde T_f$
reconstructed by (11.17) on the same stencil. Energy is never reconstructed
separately: an isothermal composition profile advected across the mesh stays
isothermal, which the tests assert to $10^{-9}$ K under the equidistribution
correction in both energy forms, and a uniform composition stays uniform. The
conversion-history pair moves with the solid: with $F_{s,f} = \sum_{j\,\text{solid}}
F_{j,f}$ the solid mass flux, the history fluxes are $F_{s,f}\,(G/m_s)_f$ and
$F_{s,f}\,(m_{\text{dry}}/m_s)_f$, with $(G/m_s)_f$ and $(m_{\text{dry}}/m_s)_f$ the
reconstructions of the per-solid-mass values $G_i/m_{s,i}$ and $m_{\text{dry},i}/m_{s,i}$ (the solid mass smoothly floored at
the cell's mass scale). Under a prescribed mesh motion that keeps the end nodes in
place the column sums of the pair are then conserved to round-off, which a test
asserts.

These fluxes are a continuous Lagrange–remap: the numerical diffusion a discrete
remap would pay at each remeshing is paid here continuously, in smooth form, and only
where $\mathbf c \ne 0$.

### 11.6.5 More dimensions

On a mesh of quadrilaterals or hexahedra (residual level, §11.3.5) the relative fluxes
are first-order upwind: through each face every inventory moves at the solid's
velocity relative to the face, $v_{\text{rel},f}A_f = -\sum_q\mathbf s_{q,f}\cdot\mathbf
c_q$ over the face's nodes, carrying the upwind cell's positive-part concentrations,
their enthalpy at its temperature and its per-solid-mass history values; a boundary
face carries its own cell's whichever way the mesh crosses it.

---

## 11.7 The equidistribution correction

![Uneven cells relax within a reacting layer above an inert substrate; the implementation activation curve rises as remaining thickness falls from 0.3 to 0.1.](../../assets/figures/04_equidistribution.svg)

**Figure 4. Equidistribution within a vanishing layer.** The left schematic isolates redistribution toward equal cell widths. Correction is zero in the inert substrate and at layer boundaries; those boundaries can still move with the solid. The curve evaluates the implemented thickness-activation switch. Its input is remaining **thickness**, distinct from the **condensed-mass** criterion for burnout in Figure 5.

### 11.7.1 Purpose

On a mesh that rides the skeleton, the cells of a gasifying layer thin at rates set
by their own heating: the cell at the surface first, the deeper ones later. Under
deactivation at the mass-fraction threshold $\varepsilon_m$ (Chapter 14) they then
reach the threshold one at a time, each an event that restarts the integrator, and the
last cells of the layer become very thin and stiff before they go. The
equidistribution correction moves the mesh relative to the solid inside a layer that
deactivation can remove, so that the layer's cells thin together, reach the threshold
together and, with the near-threshold band, leave in one event
(`EquidistributionCorrection`, `src/Closures/correction.jl`).

### 11.7.2 The node correction

At an interior node $i$ of a layer of $n_{\text{layer}}$ cells, between cell $i-1$
(inward) and cell $i$ (outward),
$$
\boxed{\;
c_i = a_i\,\frac{n_{\text{layer}}^2}{\tau_{\text{eq}}}\,
\big(\tilde{\mathcal M}_i\,\Delta z_i - \tilde{\mathcal M}_{i-1}\,\Delta z_{i-1}\big)
\;}
\tag{11.21}
$$
(`_correction!`), with $\Delta z_i = z_{i+1} - z_i$ the current cell thickness,
$a_i$ the activation of §11.7.3, $\tilde{\mathcal M}$ the normalized monitor of
§11.7.4 and $\tau_{\text{eq}}$ [s] the layer's equidistribution time (the keyword
`τ`, §11.7.6). A positive $c_i$ moves the node in $+z$ relative to the solid, into the
outer cell: the correction moves each node into the thicker (or heavier) of its two
cells, so $\tilde{\mathcal M}\Delta z$ relaxes toward a common value over the layer.

With $a = 1$ and $\tilde{\mathcal M} = 1$, the correction's share of the thickness rates
is $d\Delta z_i/dt|_c = c_{i+1} - c_i = (n_{\text{layer}}^2/\tau_{\text{eq}})
(\Delta z_{i+1} - 2\Delta z_i + \Delta z_{i-1})$, a discrete diffusion of the
thickness along the layer with zero flux at its bounding nodes ($c = 0$ there,
§11.7.5). Its slowest non-uniform mode decays at
$(n_{\text{layer}}^2/\tau_{\text{eq}})\cdot2\big(1 - \cos(\pi/n_{\text{layer}})\big)
\approx \pi^2/\tau_{\text{eq}}$: the layer equalizes in about $\tau_{\text{eq}}/\pi^2$
whatever its cell count, which is why the gain carries $n_{\text{layer}}^2$. A test
integrates three inert cells under a fully active correction and matches the cell
thicknesses to the matrix exponential of that linear system to $10^{-7}$, the solid
at rest and isothermal and its mass conserved to $10^{-12}$.

### 11.7.3 The activation

The correction acts only where cells have thinned. Each cell has a switch
$$
a_i = \mathcal S\!\left(\frac{\ln r_{\text{hi}} - \ln\big(\Delta z_i/\Delta z_i(0)\big)}{\ln r_{\text{hi}} - \ln r_{\text{lo}}}\right),
\qquad \mathcal S(x) = x^2(3 - 2x)\ \text{clamped to}\ [0,1],
\tag{11.22}
$$
with `activation = (r_hi, r_lo)`, by default `(0.3, 0.1)`: zero while a cell keeps at
least $r_{\text{hi}}$ of its initial thickness, one at and below $r_{\text{lo}}$, and
$C^1$ in $\ln\Delta z$ between. The thickness is floored smoothly at $10^{-8}$ of its
initial value before the logarithm, so the switch is one, with a finite derivative, at
zero and negative thicknesses. A node's activation is the smooth "or" of its two
cells', $a_i \leftarrow 1 - (1 - a_{i-1})(1 - a_i)$, so a node moves as soon as either
neighbour has thinned. Where a layer that keeps part of its thickness by design is
moved (`layers = :all`, §11.7.5), the thresholds must sit below that thinning: a char
that keeps 20 % of its thickness would otherwise activate the correction, which is
why they are keywords and not constants.

### 11.7.4 The monitor

The monitor is $\mathcal M = \alpha_M + \beta_M\,\rho_c/\rho_{\text{ref}}$, with
`monitor = (α_M, β_M)`, $\rho_c$ a cell's condensed mass per volume and
$\rho_{\text{ref}}$ the largest initial one over the mesh. It is divided by its mean
over the layer's cells at the initial state, so $\alpha_M$ and $\beta_M$ set its
shape and not the relaxation time. The default `(1.0, 0.0)` gives
$\tilde{\mathcal M} = 1$: the layer's cells relax toward equal thickness. A mass term
($\beta_M > 0$) keeps the nodes closer to Lagrangian where the composition varies.

### 11.7.5 Layers and the nodes that do not move

A layer is bounded by the pinned node, the mesh's **Lagrangian nodes** (the layer
interfaces) and the end nodes. The correction is zero at each of them, at every node
beside a dead cell (so the outermost live node moves with the solid), and wherever both
cells keep at least $r_{\text{hi}}$ of their initial thickness. The exposed surface
node, an end node, therefore always moves with the solid.

The Lagrangian nodes are those of the mesh plus the solve's `lagrangian` keyword.
Under an equidistribution correction every interior node between cells that start
with different condensed components must be among them, or the correction would mix
the layers' compositions: `setup` refuses a correction that would move such a node,
and by default (`lagrangian = nothing`) takes exactly those nodes.

With `layers = :vanishing` (the default) only the layers that deactivation can remove
move; the others keep their mesh with the solid. A layer can vanish when some of its
cells hold condensed mass and each that does keeps less than the fraction
$\varepsilon_m$ of it once its components are converted completely (`vanishes`). The
fraction each condensed component leaves (`residues`) is the condensed yield of a
reaction that consumes it, each product weighted by its own residue, the largest over
the reactions that consume it (competing reactions: the least gasifying bounds what is
left); it is zero for a component every consuming reaction gasifies and one for an
inert component. `layers = :all` moves every layer, for decks that exercise the
correction on layers that cannot vanish.

### 11.7.6 The layer time

The equidistribution time is one value for every layer when `τ` is a number. By
default it is `ActivationTime()`: each layer's time is set once, when its first cell
activates (falls below $r_{\text{hi}}$ of its initial thickness), from the layer's
thickness $L_{\text{layer}}$ (from its inner bounding node to its outermost live node)
and its thinning rate $\dot L_{\text{layer}}$ (their velocity difference) at that time,
$$
\boxed{\;
\tau_{\text{eq}} = \lambda_\tau\,\frac{L_{\text{layer}}}{|\dot L_{\text{layer}}|}
\;}
\tag{11.23}
$$
with $\lambda_\tau$ the keyword `λ` of `ActivationTime`, by default
`ACTIVATION_TIME_FACTOR` = 1.8. The rate scales with the physical rates and
$\lambda_\tau$ is a pure number, so a deck whose rates are all multiplied by $s$ gets
$\tau_{\text{eq}}/s$, which a test asserts in dimensionless time. With it the 5.9 mm
PMMA layer of V-BO-01 at 50 kW/m² gets about 720 s, and a layer equalizes in about
$\tau_{\text{eq}}/\pi^2 \approx 0.18$ of the life $L_{\text{layer}}/|\dot L_{\text{layer}}|$
the estimate implies, about a fifth of it. Before the estimate the layer's correction is zero,
since no cell of it is active. A layer thinning no faster than `floor` (default 1)
times its initial thickness over the solve's horizon keeps the `fallback` time
(`EQUIDISTRIBUTION_TIME` = 720 s, with which V-BO-01 ends in one deactivation at 30,
60 and 120 cells across the layer) until a cell crosses $r_{\text{lo}}$, where the
estimate is taken once more; a layer that thickens keeps the fallback. A layer is
estimated once: a cell that re-activates later changes nothing.

The estimate is an event of the integrator. A continuous callback watches the smallest
margin $\Delta z_i/(r\,\Delta z_i(0)) - 1$ over the live cells of the layers still
waiting ($r = r_{\text{hi}}$ before the first estimate, $r_{\text{lo}}$ while a floored
rate waits for the second), locates its first root inside the step, and sets the time
of every waiting layer whose margin is within `reltol` of zero
(`estimate_layer_times!`). The state is not modified, and the layer's correction is
zero up to its activation, so the event changes the right-hand side only where the
activation has just begun; the integrator re-evaluates its derivative there under the
policy's flag rule (Chapter 15). Every `integrate` call starts with every layer back
before its first activation, so repeated integrations of one problem are
bit-identical.

### 11.7.7 Installation, keywords and outputs

`setup` installs `EquidistributionCorrection()` when no correction is given and three
conditions hold: the mesh moves (`use_ale = true`), cells are deactivated
(`deactivation = true`), and deactivation at `threshold` can remove some cell (its
components, converted completely, leave less condensed mass than $\varepsilon_m$ of
it). Otherwise the solve gets `NoCorrection()`. When the bound correction moves some
layer, the near-threshold band defaults to 0.1 (Chapter 14). A correction passed as
`correction` takes precedence over the rule; `NoCorrection()` keeps the mesh with the solid; any
correction other than `NoCorrection` needs `use_ale = true`. The numerical
parameters of the burnout path and the envelope over which it is verified are listed
in User Guide Chapter 10.

```julia
using Pyrolysis

pmma = SolidComponent(:PMMA, ρ = 1190.0, φ=.3, c = 1420.0, k = 0.19, ε = 0.9, λ = 2e-5)
board = SolidComponent(:board, ρ = 256.0, φ=.8, c = 1100.0, k = 0.06, ε = 0.9, λ = 2e-5)
mma = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02,)
material = Material(
    name = :PMMA_on_board,
    components = (pmma, board, mma),
    reactions = (Reaction(:gasification, 1 => 3, A = 8.5e12, E = 188e3, h = 8.7e5),),
)
bcs = BoundaryConditionSet(
    thermal = (
        top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)
# 2 mm of PMMA (30 cells) on 2 mm of an inert board (30 cells).
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.004, 60),
    material = material,
    bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> z > 0.002 ? 1190.0 : 0.0, z -> z > 0.002 ? 0.0 : 256.0, z -> 0.0],
    tspan = (0.0, 200.0),
    output_times = collect(0.0:5.0:200.0),
)

# The rule installs the equidistribution correction: the mesh moves, cells are
# deactivated, and the PMMA gasifies completely.
sol = solve(problem; use_ale = true, deactivation = true)

# The same correction, every parameter spelled out.
corr = EquidistributionCorrection(
    τ = ActivationTime(λ = 1.8, fallback = 720.0),
    activation = (0.3, 0.1),
    monitor = (1.0, 0.0),
    layers = :vanishing,
)
sol = solve(problem; use_ale = true, deactivation = true, correction = corr)

sol.statistics.invocations    # deactivation invocations
sol.correction_velocity       # c at every node and saved time [m/s]
sol.extras.min_cell_thickness # the thinnest live cell at every saved time [m]
sol.definition.closure.correction.time    # each moving layer's τ_eq [s]
sol.definition.closure.correction.set_at  # the time it was set [s]
```

Measured with this script: the default call installs the correction with one moving
layer, the PMMA, bounded by the board–PMMA interface (node 31), which `setup` takes
as a Lagrangian node from the initial composition; the layer's equidistribution time
is set to about 172 s at $t \approx 35$ s, and the 30 PMMA cells leave in two
invocations, at $t \approx 91.9$ s and $96.6$ s, the solve then continuing on the
board to the end of its time span (return code `:Success`). The explicit correction
reproduces the default bit for bit. At `threshold = 3e-3` the layer leaves in one
invocation. With `correction = NoCorrection()` a PMMA cell falls below the threshold
inward of a live one at $t \approx 69$ s and the solve ends with outcome `:interior`
(Chapter 14).

### 11.7.8 The prescribed mesh velocity (`w_override`)

For verification by manufactured solutions the `w_override` keyword of `solve` (a
function `w(z, t)`) prescribes the velocity of every node instead
(`PrescribedCorrection`): every velocity row is $0 = w_i - w(z_i, t)$ and the correction
is the node velocity itself, $\mathbf c = \mathbf w$, so the mesh moves over a solid at
rest, $\mathbf v_s = 0$. Every inventory then crosses every face, the boundary faces
included, at $-\mathbf w$ through the relative fluxes. The material must not react (a
reacting solid would move), cells are not deactivated, no other correction may be
given, and the keyword needs `use_ale = true`. V-MMS-02 and V-MMS-03 use it
(§11.5.2).

---

## 11.8 Mesh safety

ALE motion can in principle tangle the mesh ($V_i\le 0$) if a node overtakes its
neighbour, for instance under an aggressive step near complete burnout. Four
mechanisms keep a solve inside the physical domain:

- **Smooth guards on trial states.** Newton iterates leave the physical domain before
  any step can be rejected, so every quantity the residual evaluates is finite with a
  finite derivative at any real state: the volume used for $m/V$ is smoothly floored at
  $10^{-8}$ of the initial volume, every centre-to-face distance at $10^{-8}$ of its
  reference value, and the activation's thickness likewise (Chapter 13).
- **The domain test.** A proposal is inadmissible when a live cell's nodes have
  crossed (a non-positive volume), when the cross-section area is non-positive, when a
  live component mass is more negative than its tolerance, or when a live cell is at or
  below absolute zero, all read from the raw state without a guard; the integrator's
  domain test rejects such a step (`DomainCheck`, Chapter 15).
- **Deactivation.** A cell whose condensed mass falls below $\varepsilon_m$ of its
  initial condensed mass is deactivated before it collapses, its inventories and its
  volume going to the live cell beneath (Chapter 14); the correction of §11.7 keeps the
  cells of a vanishing layer from thinning one at a time.
- **The thickness stop.** The `min_thickness` keyword of `solve` (default
  $10^{-6}$ m) ends the solve with return code `:Terminated` and outcome
  `:min_thickness` through a continuous callback when the distance from the first node
  to the outermost live one falls below it.

---

## 11.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

All three reference models can change their 1D material geometry. Describing
them as fixed-grid codes with no physical surface motion is incorrect.

**ThermaKin2Ds.** TK §3 describes deformation of element thicknesses followed
by size control before each timestep: large elements split, and smaller ones
receive part of a neighboring element with species and energy conservation.
Its Eulerian continuum deformation term does not imply fixed numerical
elements. Deformation is restricted to the primary coordinate.

**Gpyro.** GP Eqs. 56–57 conserve condensed inventories on cells whose
thicknesses change. Only gas crosses the material cell boundaries in this
1D formulation. Nonuniform initial spacing and evolving spacing are distinct
features. Setting `use_ale=false` here does not reproduce Gpyro's deforming
grid. Multidimensional transport must be distinguished from general mesh
deformation; the inspected V2 `UPDATE_MESH` includes multiple deformation
modes and explicitly flags limitations in its 2D/3D distortion treatment.

**FDS.** FDS 6.11.0 updates solid-cell thicknesses and coordinates and remeshes
layers (`REGRID_FACTOR`, `X_S_NEW` and `REMESH_LAYER` in `wall.f90`). This
is motion of the internal solid grid; it is separate from the surrounding
CFD mesh and obstruction representation.

**Pyrolysis.jl.** Node positions are differential unknowns and node velocities
are algebraic unknowns in the same coupled system as mass and energy. The
equidistribution correction transports inventories relative to that mesh;
the geometry obeys the discrete volume identity of §11.5. These are specific
implementation differences, not proof that the reference codes lack geometric
conservation. Volume laws, remapping, radiation and gas transport still have
to be matched before comparing solutions (§§8–10).

Sources: TK §3; GP §4.1; [FDS wall source](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90#L2388-L2679);
[V2 mesh source](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_mesh.f90).

---

## 11.10 Limitations and open issues

The principal limitations and open issues are:

1. **Dimensionality.** Full solves run on a line. The mesh, the state layout, the
   closures (`PotentialClosure`, the column mode of `BalanceRowClosure`), the
   first-order relative fluxes and the deformation-gradient rows exist on meshes of
   quadrilaterals and hexahedra at the residual level (§11.3.5), but gas transport,
   in-depth radiation, the variable cross section, the equidistribution correction and
   deactivation are a line's, and `setup` refuses them elsewhere.

2. **Separation of axial and lateral motion.** A variable cross section is a
   column-global area following an empirical law of the column's conversion, closed by
   two algebraic rows (11.11) that read every cell: the system's two dense rows. A
   resolved lateral deformation is not modelled, and the law needs planar coordinates.

3. **Density evaluation in the strain rate.** The strain rate (11.7) uses $\rho_j(T)$
   at the cell's trial temperature; whether a time-averaged density would improve
   consistency near steep thermal transients is an open modeling question, as is the
   extension of the volume-change accounting to multi-reactant stoichiometry.

4. **Free stream beyond the planar line.** The uniform-field property holds to
   round-off only on a planar line with a constant cross section; in cylindrical and
   spherical coordinates and under a variable cross section it holds to the
   integrator's order (§11.5.2).

5. **The correction is a numerical device with a verified envelope.** The
   equidistribution time is set from one estimate per layer, the activation
   thresholds and the monitor are user parameters, and the relative fluxes it switches
   on add smooth numerical diffusion wherever $\mathbf c \ne 0$. The burnout path is
   verified on the decks listed in User Guide Chapter 10; other materials, thicker or
   charring layers and higher heating lie outside that envelope until checked by the
   same refinement study.

---

### Cross-references

- The volumetric strain rate $\theta$, the swelling factor $\gamma_j$, the lateral
  shrinkage law and the conversion-history pair are developed in **Chapter 10, Volume
  Change & Lateral Shrinkage**.
- The energy assembly, $\rho c_p^{\text{eff}}$ (every component), the gas enthalpy flux
  at the faces and the harmonic-mean face conductivity are in **Chapter 7, Heat
  Conduction, Convection & Energy Assembly**.
- The partial-pressure gas flux and the Darcy velocity are in
  **Chapter 9, Gas-Phase Mass Transport & Pressure**.
- The face-based mesh, the metrics, the state layout and mapping, the kernels and the
  row assembly are in **Chapter 13, Finite-Volume Discretization**.
- Cell deactivation, the near-threshold band and the transfer into the cell beneath
  are in **Chapter 14**.
- The mass-matrix DAE, the integrator policies, the domain test and the colored
  Jacobian of the coupled system are in **Chapter 15**.
- The mass and energy ledgers are in **Chapter 16, Conservation Ledgers and the Audit**.
