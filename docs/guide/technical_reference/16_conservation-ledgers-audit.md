# 16. Conservation Ledgers and the Audit

This chapter derives the conservation ledger that Pyrolysis.jl keeps of a
solve: the mass and sensible-enthalpy **inventories** of the live cells, the
**rates** that change them, the **recorder** that integrates those rates over
the accepted steps, the **audit** that closes the inventories against the
integrated rates and the deactivation events, and the **event records** that
carry what an event takes out of the live region. For every quantity we give
the continuous balance law, its discrete counterpart as the code assembles it,
the sign conventions, and the conditions under which the closure is exact,
holds at the integrator's order, or is limited by a quadrature.

All symbols follow the unified nomenclature of §2 (Nomenclature). In this
chapter `i` is the **cell index**, `j` the **component index** (1…$N_C$), `r`
the **reaction index** and `f` the **face index**. The three binding sign
conventions are in force: a face flux is the transport per unit area along the
face normal, which points from the owner cell to the neighbour ($+z$ on every
interior face and on `:top`, out of the domain on `:bottom`); the discrete
divergence $(1/V_i)\sum_{f\in\partial i}s_{if}A_fF_f$, with $s_{if}=+1$ when the
normal of face $f$ points out of cell $i$ and $-1$ otherwise, is positive for
*outflow*; and $h>0$ is endothermic, with the volumetric heat source
$Q_\mathrm{rxn}=-\sum_r h_r r_r$ (§6, Reaction Kinetics). The ledger itself
counts the mass **outflow** positive when mass leaves and the enthalpy **rate**
positive when energy enters; both are restated where they are used.

The implementing code is the `Outputs` module (`src/Outputs/conservation.jl`:
`ledger_rates`, `live_inventories`, `ConservationRecorder`,
`conservation_audit`), the recorder's binding to the integrator in
`Integration` (`src/Integration/observers.jl`), and the event records of
`Deactivation` (`src/Deactivation/events.jl`). `ConservationRecorder` and
`conservation_audit` are exported by `Pyrolysis`; `ledger_rates` and
`live_inventories` are reached through `Pyrolysis.Internal`.

---

## 16.1 Purpose and design philosophy

A pyrolysis solver couples a stiff reaction network, in-depth radiation, gas
transport, a moving mesh and the deactivation of burnt-out cells, all advanced
by an implicit time integrator. With this many interacting mechanisms a silent
error in flux or source assembly can masquerade as physics. The ledger makes
conservation **observable on every solve**: a solve can carry a recorder that
integrates the ledger's rates over its accepted steps, and every solution
carries the saved series from which an audit can be computed without one.

Two properties of the discretization make the ledger both cheap and sharp.
The **component masses of every cell are states** (§3, §13), so mass enters
and leaves a cell only through its face fluxes and its reaction source, and
the stoichiometry makes the reaction sources of each cell sum to zero over the
components (§6.5.3). The **energy row is assembled from face energy fluxes**
scattered with opposite signs into the two cells of each face, so the interior
faces cancel in any sum over cells. The inventories of the live region
therefore change only by what crosses its boundary faces, by the volumetric
sources of its cells, and by what the deactivation events remove. The ledger
is built from exactly those terms; the audit asks whether the integrator's
trajectory, the integrated boundary and source rates, and the event records
add up.

A non-closing audit points either to the integration (a tolerance too loose
for the quantity audited, a step the recorder could not account for) or to a
**bookkeeping defect**: a flux scattered with one sign, a source missing from
a row, an event that moved an inventory without recording it. What the
ledger cannot prove is that its terms are the right physics. It reads the same
kernel contributions the residual assembles, so **a mechanism counted twice in
the residual is counted twice in the ledger rate as well, and the audit still
closes.** That class of error is caught by checks with an independent
reference: the term-by-term test of the energy row (§16.5.5) and the
exact-solution cases of §18.

> **Implementation note.** The ledger functions are pure functions of a state
> and of the derived state and kernel contributions of one residual evaluation
> at it: `ledger_rates(def, c, ds)` and `live_inventories(def, ext, ds, u)`.
> They read which cells are live from the derived state (the layout's mask)
> and never modify the state. The recorder evaluates them in a workspace of
> its own (§16.8), so the integrator's buffers and trajectory are untouched.

---

## 16.2 The continuous balance laws being audited

Let $\Omega(t)$ be the **live region** of the column: the cells not yet
deactivated, bounded by the back face $z=0$ and by the face of the outermost
live cell (the exposed surface, $z=L(t)$), with cross-sectional area $A(t)$. On
a moving mesh the end nodes and the outermost live node move with the solid
(the mesh correction vanishes there, §11), so $\partial\Omega$ is a material
boundary of the condensed phase; on a fixed mesh it is fixed between events.
(A prescribed mesh velocity, the verification hook `w_override`, moves the
mesh over a solid at rest instead; the faces' relative fluxes then carry
inventory across $\partial\Omega$, and the ledger counts them with the other
boundary fluxes.)

### 16.2.1 Mass

Integrating the component conservation law (§3.2) over $\Omega(t)$, for each
component $j$:

$$
\frac{d}{dt}\int_{\Omega(t)} \xi_j \,dV
= -\oint_{\partial\Omega} J_j\,\hat{\mathbf n}\cdot d\mathbf A
\;+\; \int_{\Omega(t)} \omega_j\,dV ,
$$

with $J_j$ the mass flux of component $j$ relative to the solid (zero for a
condensed component) and $\omega_j=\sum_r\nu_{r,j}r_r$ its net production by
reaction. The single-reactant stoichiometry conserves mass exactly,
$\sum_j \omega_j = 0$ (§6.5.3), so the
**total** mass changes only through the boundary:

$$
\frac{d}{dt}\sum_j\int_{\Omega(t)} \xi_j \,dV
= -\oint_{\partial\Omega} \sum_j J_j\,\hat{\mathbf n}\cdot d\mathbf A .
$$

### 16.2.2 Energy

Every component's sensible enthalpy is taken from one datum,

$$
h_j(T)=\int_{T_\mathrm{ref}}^{T} c_{p,j}(T')\,dT',
\qquad T_\mathrm{ref}=298.15\ \mathrm K ,
$$

and the sensible enthalpy of the live region is
$H=\int_{\Omega}\sum_j \xi_j\,h_j(T)\,dV$. Its balance is

$$
\frac{dH}{dt}
= -\oint_{\partial\Omega}\Big(q_\mathrm{cond}+\sum_j h_j(T)\,J_j\Big)\hat{\mathbf n}\cdot d\mathbf A
\;+\; \int_{\Omega}\Big(Q_\mathrm{rxn}+Q_\mathrm{rad}+\sum_j \omega_j\,h_j(T)\Big)dV ,
$$

with $q_\mathrm{cond}=-k_\mathrm{eff}\,\partial T/\partial z$ the conducted
flux, $\sum_j h_j J_j$ the sensible enthalpy the mass fluxes carry,
$Q_\mathrm{rxn}=-\sum_r h_r r_r$ the reaction heat with each $h_r$ quoted at
the local temperature (§6.6.3), $\sum_j \omega_j\,h_j(T)$ the **formation
term** (the sensible enthalpy a reaction moves from its reactant to its
products at the local temperature), and $Q_\mathrm{rad}$ the in-depth
radiation the cells absorb (nonzero under `BEER_LAMBERT` only; under
`SURFACE_ABSORPTION` the absorbed flux is part of the boundary condition's
flux). At the exposed surface the gas leaves at the surface temperature $T_s$,
so the carried term there is $\sum_j h_j(T_s)J_j$; the surface energy
balance itself has no gas-enthalpy term (§7, §12.2).

This is the law the energy ledger discretizes, and it is a balance of the
conserved quantity: the change of composition at fixed temperature (the term
$\sum_j h_j\,\partial\xi_j/\partial t$ of the product-rule expansion of
$\partial_t\sum_j\xi_j h_j$) is part of it through the formation and
carried terms, and nothing is approximated.

---

## 16.3 Mass ledger

### 16.3.1 The live mass inventory

The mass inventory of the live region is the sum of the mass states of the
live cells,

$$
m_\mathrm{live}(t) = \sum_{i\in\mathrm{live}}\ \sum_{j=1}^{N_C} m_{j,i}(t),
$$

with $m_{j,i}$ [kg] the mass of component $j$ in cell $i$. No quadrature
enters: the inventory is a sum of states. A deactivated cell's masses are
zeroed by the event that took it (they went to its receiver or were recorded
as lost, §16.7), and its rows are masked; the sum runs over the cells the mask
holds live.

### 16.3.2 The outflow rate

The mass leaving the live region per unit time is

$$
\dot m_\mathrm{out}(t) = \sum_{f\in\mathcal F_\mathrm{live}}\ \sum_{j} F_{j,f}\,A_f ,
$$

positive when mass leaves. $\mathcal F_\mathrm{live}$ is the set of **live
boundary faces**: for each boundary face of the mesh, the face that bounds the
live cells on its side (the face itself while its cell is live; after cells
beyond it are deactivated, the face of the outermost live cell; none once no
live cell borders it). $F_{j,f}$ [kg·m⁻²·s⁻¹] is the mass flux of component $j$
through $f$ along the outward normal, as the kernels wrote it: the gas
transport's boundary flux under the face's mass condition (§9, §12.4, §12.5),
and, under a prescribed mesh velocity, the relative flux of the moving face.
$A_f$ is the face's area in the current state, so a variable cross section is
integrated at its current value.

**The discrete identity.** The mass row of cell $i$ is
$dm_{j,i}/dt = \omega_{j,i}\tilde V_i - \sum_{f\in\partial i} s_{if}F_{j,f}A_f$,
with $\omega_{j,i}\tilde V_i$ [kg/s] its reaction source ($\tilde V_i$ the smoothly
floored cell volume, §13.3.2) and $s_{if}F_{j,f}$ the face's flux along the normal
pointing out of cell $i$. Every interior face flux is scattered as $-F A$ into its
owner and $+F A$ into its neighbour, and no flux crosses a face between a live and a
dead cell (the boundary conditions move to the face of the outermost live
cell, and the mesh correction vanishes at its node). Summing over the live
cells and the components, the interior faces cancel and the reaction sources
vanish cell by cell, so

$$
\frac{dm_\mathrm{live}}{dt} = -\dot m_\mathrm{out}
$$

holds at **every state**, not only on the trajectory: $m_\mathrm{live}+\int\dot
m_\mathrm{out}\,dt$ is a linear invariant of the semi-discrete system, in both
energy forms. Carried as an extra state of a fixed-step solve, the integrated
outflow closes the mass rows to the round-off of the inventory (at most
$10^{-13}$ of it in the module test `test/outputs/test_energy_ledger.jl`).

> **Implementation note.** `ledger_rates(def, c, ds)` returns the pair
> `(energy, outflow)`. It visits the boundary faces of the topology
> (`topo.faces[f].cells[2] == 0`), maps each to its live face with
> `Meshes.live_face`, and sums `c.F_mass[j, f] * g.face_area[f]`. The saved
> series `sol.extras.mass_outflow` [kg/s] is this rate at every saved state,
> and `sol.extras.live_mass` [kg] the inventory.

### 16.3.3 The closure and what deactivation takes out

Between two times the inventory changes by the integrated outflow and by the
masses the deactivation events removed from the live region without passing
them to a receiver:

$$
m_\mathrm{live}(t_1)-m_\mathrm{live}(t_0)
= -\int_{t_0}^{t_1}\dot m_\mathrm{out}\,dt \;-\; \sum_{e} m^\mathrm{lost}_e ,
\qquad m^\mathrm{lost}_e=\sum_j m^\mathrm{lost}_{j,e} .
$$

On a fixed mesh an event vents the dead cell's gas (the receiver keeps its
volume, and its pores could hold that gas only at a higher pressure), so
$m^\mathrm{lost}_e$ is the vented gas; on a moving mesh the receiver keeps every
mass and $m^\mathrm{lost}_e=0$; the last live cell, which has no receiver,
loses its whole inventory (§16.7).

The audit writes the mass closure as

$$
\begin{aligned}
\text{change} &= m_\mathrm{live}(t_1)-m_\mathrm{live}(t_0),\qquad
\text{inflow} = -\int_{t_0}^{t_1}\dot m_\mathrm{out}\,dt,\qquad
\text{lost} = \sum_e m^\mathrm{lost}_e,\\
d_m &= \text{change}-\text{inflow}+\text{lost},\qquad
d_m^\mathrm{rel} = \frac{d_m}{\max\big(|\text{change}|,|\text{inflow}|,|\text{lost}|\big)} .
\end{aligned}
$$

The relative defect is scaled by the largest term of the balance, so a
component that starts at zero mass needs no separate treatment; for a solve in
which nothing changes, enters or is lost, the relative defect is undefined and
the defect itself is the quantity to read.

### 16.3.4 Sources of (small) mass-ledger discrepancy

Because $dm_\mathrm{live}/dt=-\dot m_\mathrm{out}$ holds at every state, the
mass defect of a recorded solve is the difference between the integrator's
update of the mass states and the recorder's quadrature of the outflow on the
step's interpolant (§16.5.3). It is set by the interpolant's accuracy and
falls with the tolerance: across deactivation events the module test
(`test/outputs/test_conservation.jl`) asserts a fall of more than tenfold from
`reltol = 1e-6` to `1e-8` and a relative mass defect of at most $10^{-7}$ at the
tighter tolerance, on fixed and moving meshes. Three other sources are of a
different kind:

1. **Saved-series integration.** Without a recorder the audit integrates the
   saved rates by the trapezoid (§16.8.3), whose error is the saving
   interval's; it is orders of magnitude above the recorder's (§16.8.4).
2. **Unrecorded movements.** An event that moved mass without recording it, or
   a boundary flux written on a face the ledger does not read, appears as a
   defect that does **not** fall with the tolerance. That is the bookkeeping
   failure the ledger exists to expose.
3. **External sources.** The volumetric sources of the `sources` keyword (the
   manufactured sources of the verification driver, §18.2) add mass the ledger
   does not count; the audit of such a solve is off by their integral.

---

## 16.4 Energy inventory

The energy inventory of the live region is

$$
H_\mathrm{live}(t)=\sum_{i\in\mathrm{live}} H_i(t),
$$

with $H_i$ [J] the cell's **sensible-enthalpy inventory**. In the inventory
form (`energy_state = :enthalpy`) $H_i$ is the state of the cell's energy row.
In the temperature form it is computed from the cell's temperature and masses
by the same **extended caloric relation** the inventory form inverts:

$$
H_i = \sum_j m^{+}_{j,i}\,h_j(T_i) \;+\; C_{\mathrm{reg},i}(m_{\cdot,i})\,\big(T_i-T_\mathrm{ref}\big).
$$

Here $m^{+}_{j,i}$ is a one-sided C¹ positive part of the mass at the cell's
mass scale (exactly zero at and below zero, exactly $m_{j,i}$ at and above the
scale), and $C_{\mathrm{reg},i}$ [J·K⁻¹] a C¹ regularizing capacity that is zero
while the guaranteed heat capacity of the cell's masses is at least twice the
cell's capacity floor and rises to the floor as that capacity vanishes (§13).
The relation is strictly increasing in $T_i$, so it has one root for every real
inventory. Away from empty cells (each component absent or holding at least
the mass scale, and the guaranteed capacity at least twice the floor) it is the
cell's sensible enthalpy, $H_i=\sum_j m_{j,i}\,h_j(T_i)$.

Both energy forms are therefore audited on the same function of the state,
and their audits are comparable term by term. The inventory is measured from
the common datum $T_\mathrm{ref}$; since every component shares it, a transfer of
mass between cells or between components moves no datum energy. The extended
relation, rather than the bare $\sum_j m_{j,i}h_j$, is what a
deactivation event conserves: the merged temperature is the root of this
relation at the summed sensible enthalpy, and the regularization part, which
only a cell far below its capacity floor carries, is written to the event
(§16.7).

> **Implementation note.** `live_inventories(def, ext, ds, u)`,
> `src/Outputs/conservation.jl`, with the caloric extension of the solve
> (`CaloricExtension`, `src/States/caloric.jl`, built from the cells' mass
> scales and capacity floors). It returns `(energy, mass)`.
> `sol.extras.live_energy` [J] is the energy inventory at every saved state.

---

## 16.5 The energy ledger

### 16.5.1 The enthalpy rate, from the energy row

The energy flux through face $f$ along its normal is conduction plus the
sensible enthalpy its mass fluxes carry,

$$
F^E_f = q_f + \sum_j F_{j,f}\,h_j(T_f)
\qquad [\mathrm{W\,m^{-2}}],
$$

with $T_f$ the face temperature of the flux that carries the mass: the
distance-weighted mean of the two cell temperatures for the gas transport
between cells, the limited upwind face value for the relative flux of a moving
mesh, and the surface temperature $T_s$ on a boundary face. On a boundary face
$q_f$ is minus the heat conducted from the face into its cell (§16.6). With
$S_{E,i}$ [W] the cell's energy sources (the reaction heat
$Q_{\mathrm{rxn},i}\tilde V_i$ and the absorbed in-depth radiation, §13.3.4), the
cell's energy right-hand side is
$E_i=S_{E,i}-\sum_{f\in\partial i}s_{if}F^E_fA_f$. The energy row is (§7, §13)

- inventory form: $\displaystyle \frac{dH_i}{dt} = E_i + \sum_j \omega_{j,i}\tilde V_i\,h_j(T_i)$;
- temperature form: $\displaystyle C_i\,\frac{dT_i}{dt}
  = E_i - \sum_j h_j(T_i)\,\dot M_{j,i}$,

where $\dot M_{j,i}=-\sum_{f\in\partial i}s_{if}F_{j,f}A_f$ is the face-flux part
of the mass rate and $C_i$, the heat-capacity sum $\sum_j m_{j,i}c_{p,j}$, is
smoothly floored and evaluated at the trial temperature. By the product rule,
the temperature form gives
$\frac{d}{dt}\sum_j m_{j,i}h_j(T_i)
=\sum_j m_{j,i}c_{p,j}\,\dot T_i+\sum_j h_j\,(\omega_{j,i}\tilde V_i+\dot M_{j,i})
=E_i+\sum_j \omega_{j,i}\tilde V_i\,h_j(T_i)$: the same rate. Summed over the
live cells, the interior faces cancel and the **enthalpy rate** of the live
region is

$$
\dot H_\mathrm{in} = \sum_{f\in\mathcal F_\mathrm{live}}\big(-F^E_f A_f\big)
\;+\;\sum_{i\in\mathrm{live}}\Big[S_{E,i}+\sum_j \omega_{j,i}\tilde V_i\,h_j(T_i)\Big],
$$

positive when energy enters: the heat conducted in through the live boundary
faces less the enthalpy carried out through them, plus the reaction heat, the
absorbed in-depth radiation and the formation terms of the live cells. It is
the discrete form of the law of §16.2.2 term by term.

> **Implementation note.** `ledger_rates(def, c, ds)` accumulates
> `E -= c.F_energy[f] * g.face_area[f]` over the live boundary faces and
> `E += c.S_energy[i] + Σⱼ c.S_mass[j, i] * ds.h[j, i]` over the live cells
> (`S_mass` holds $\omega_{j,i}\tilde V_i$ [kg/s], `ds.h` the sensible enthalpies
> $h_j(T_i)$ at the cell temperatures). `sol.extras.energy_rate` [W] is this
> rate at every saved state.

### 16.5.2 Exactness in the two energy forms

**Inventory form.** $\dot H_\mathrm{in}$ is the sum of the live energy rows, so
$H_\mathrm{live}-\int\dot H_\mathrm{in}\,dt$ is a linear invariant of the
semi-discrete system, exactly as the mass ledger's. Carried as an extra state
of a fixed-step solve, the ledger integral closes the inventory to round-off
(at most $10^{-13}$ of the integrated rate in `test/outputs/test_energy_ledger.jl`,
on a gasifying column whose integrated rate exceeds $10^5$ J, under Fickian and
Darcy–Fick transport).

**Temperature form.** $H_\mathrm{live}$ is a nonlinear function of the state,
and the integrator advances $T$ and $m$, not $H$. The identity
$dH_\mathrm{live}/dt=\dot H_\mathrm{in}$ holds at every state at which the floored
capacity sum equals the inventory's slope $\partial H_i/\partial T_i$ (away from
empty cells and with the temperature inside the property bracket, where the
trial temperature is the temperature), but a step's change of
$H_\mathrm{live}$ departs from the integral of the rate by the method's local
error. The closure therefore **converges at the integrator's order**: in the
same module test, with fixed steps of a fourth-order method halved three
times, the defect falls at observed orders rising toward four (2.1, 2.5, 3.2
under Fickian transport; 2.7, 3.3, 4.3 under Darcy–Fick), the rise being the
pre-asymptotic approach.

### 16.5.3 Per-step accumulation: the recorder

Over each accepted step $[t_n,t_{n+1}]$ of length $\Delta t$ the recorder
integrates both rates on the step's interpolant $\mathbf u_h(t)$ by
three-point Gauss–Legendre quadrature,

$$
\int_{t_n}^{t_{n+1}}\dot H_\mathrm{in}\,dt \approx
\frac{\Delta t}{2}\sum_{k=1}^{3} w_k\,\dot H_\mathrm{in}\big(\mathbf u_h(t_k),t_k\big),
\qquad
t_k=\frac{t_n+t_{n+1}}{2}+\frac{\Delta t}{2}x_k ,
$$

with $x=(-\sqrt{3/5},\,0,\,\sqrt{3/5})$ and $w=(5/9,\,8/9,\,5/9)$, and the
outflow likewise. Each rate evaluation is one residual evaluation at the
interpolated state in the recorder's own workspace. The rule is exact for
polynomials of degree five in $t$, so the accumulated integral carries the
interpolant's error, which the integrator's tolerance controls.

**Steps that an event ends early.** A located event ends a step at the first
root of its condition (§14, §15). The recorder reads every step's end before
the step's events are applied, through a continuous callback whose condition
never crosses zero, so that it sees each step and never fires. When one of the
event conditions it knows crosses zero downward inside the step (the threshold
margin of a located deactivation, the activation margin of the layer times,
the thickness margin of the `min_thickness` stop), it locates the first
crossing by 60 bisections on the interpolant and integrates up to it as well.
After the step's events it **commits**: the whole step's integrals if the step
ran to its end; the integrals up to the crossing if an event ended the step
there (the step's end within $10^{-10}$ relative of the crossing); `NaN`
otherwise, which is the case of a step that a callback of the caller's ended
early. It then records the inventories at the state after the events. The
initial inventories are those of the set-up state at the start of the time
span.

This ordering keeps an event's state change out of every step's integral: the
inventories after the event include it, and the event record carries what it
removed (§16.7).

> **Implementation note.** `_BoundRecorder`, `_step_integrals`,
> `_read_step!`, `_commit_step!` and `_reader` in
> `src/Integration/observers.jl`. The commit runs in the last callback of the
> solve's callback set, after every other callback of the step.

### 16.5.4 The residual

The audit writes the energy closure as

$$
\begin{aligned}
\text{change} &= H_\mathrm{live}(t_1)-H_\mathrm{live}(t_0),\qquad
\text{inflow} = \int_{t_0}^{t_1}\dot H_\mathrm{in}\,dt,\qquad
\text{lost} = \sum_e\big(H^\mathrm{lost}_e+H^\mathrm{reg}_e\big),\\
d_H &= \text{change}-\text{inflow}+\text{lost},\qquad
d_H^\mathrm{rel}=\frac{d_H}{\max\big(|\text{change}|,|\text{inflow}|,|\text{lost}|\big)},
\end{aligned}
$$

with the event terms of §16.7. The module tests assert, at `abstol = 1e-8`,
`reltol = 1e-6` on a gasifying column whose escaping gas carries its enthalpy
out, $|d_H^\mathrm{rel}|\le10^{-6}$ in both energy forms and a mass defect of at
most $10^{-8}$ of the inventory; across deactivation events, on fixed and moving
meshes, a fall of both defects by more than tenfold from `reltol = 1e-6` to
`1e-8` and $|d_H^\mathrm{rel}|\le10^{-6}$ at the tighter tolerance; and on a venting
column whose cross section shrinks, relative defects of at most $10^{-7}$ (mass)
and $10^{-6}$ (energy) with the outflow integrated at the current face area.

### 16.5.5 Why it closes — and what it does *not* prove

The closure is exact for the semi-discrete system (§16.3.2, §16.5.2): its
defect is the integration error of the trajectory, the recorder's quadrature
of the rates, and nothing else. A defect that does not fall with the tolerance
is a **bookkeeping defect**: a flux computed but not scattered, a source
missing from a row or from the rate, an event that moved an inventory without
recording it, a sign error in a term.

The ledger **cannot detect a double count.** If the residual added a surface
gas-enthalpy term on top of the enthalpy the face flux already carries, the
ledger rate, read from the same contributions, would carry both and the audit
would still close. The energy row is therefore checked term by term against
an independent construction in `test/outputs/test_energy_ledger.jl`, in both
energy forms, under Fickian and Darcy–Fick transport, on uniform and stretched
meshes: every face energy flux equals conduction (Fourier between cells with
the distance-weighted harmonic-mean conductivity; the condition's flux
conducted out of the cell on the boundary) plus the enthalpy its mass fluxes
carry at the face temperature (the distance-weighted mean between cells, the
surface temperature on the boundary) to $10^{-12}$ relative; the carried term is
not small and the surface temperature departs from the cell temperature by
more than 1 K, so an omission, a double count or enthalpy taken at the cell
temperature would show; and every cell's energy row equals the sum of its
conduction, carried, source and formation terms to $10^{-11}$ of their scale.
Separately, `test/boundary_conditions/test_boundary_temperature.jl` asserts
that under a prescribed flux the solved surface temperature and conducted
flux are identical with and without a gas flux across the face: the surface
balance carries no gas-enthalpy term. The exact-solution cases of §18 are the
further, physical anchors.

---

## 16.6 Boundary-flux semantics

On a live boundary face the ledger's term is

$$
-F^E_f A_f = \Big(q_{\mathrm{in},f}-\sum_j F_{j,f}\,h_j(T_s)\Big)A_f ,
$$

where $q_{\mathrm{in},f}=k\,(T_s-T_\mathrm{cell})/d$ is the heat conducted from the
face into its cell at the surface temperature $T_s$ (§12.2). For a condition
that goes through the surface energy balance, $T_s$ is the root of

$$
q_\mathrm{BC,applied}(T_s) = q_\mathrm{cond}(T_s),
$$

so the conducted flux **is** the applied flux of the condition at $T_s$: the
absorbed incident radiation less the re-radiation and convective exchange for
a `RadiativeFluxBC`-based condition, the prescribed $q(t)$ for a `HeatFluxBC`,
the flux the prescribed temperature conducts for a `DirichletThermalBC`, zero
for an `AdiabaticBC`. The surface balance carries no gas-enthalpy term: the gas
leaving through the face takes its enthalpy at $T_s$ through the carried term
of the face's energy flux. A blowing model (§12.3.3) lowers the convective
coefficient at the face's outward mass flux; it changes $q_{\mathrm{in},f}$ and
adds no enthalpy term.

The solution exposes the components of the exposed face's flux at every saved
time (`sol.extras`, §9 of the User Guide; per unit current face area):
`surface_heat_flux_net` (the applied flux of the `:top` condition, positive
into the domain; the conducted flux under a prescribed temperature),
`surface_heat_flux_absorbed`, `surface_heat_flux_reradiation` and
`surface_heat_flux_convection` (its members), `surface_heat_flux_conduction`
($k(T_s-T_\mathrm{cell})/d$), and `back_heat_flux` for the `:bottom` condition.
Under `BEER_LAMBERT` the incident flux enters in depth and is not part of the
net surface flux; the intensity the cells absorb, $\sum_i Q_{\mathrm{rad},i}V_i$ per
unit exposed area, is `volumetric_absorption`, and the intensity leaving
through the last face of the beam is `radiation_transmitted`. The ledger counts
the absorbed part as a volumetric source of the live cells; the transmitted
part leaves the column without entering the inventory, and the incident flux
itself is in neither. Under `SURFACE_ABSORPTION` the absorbed incident flux is
part of the condition's applied flux at the face.

---

## 16.7 Deactivation events in the ledger

A cell whose condensed mass falls below the threshold is deactivated (§14):
its inventories go to the live cell inward of it, the **receiver**, and its
rows are masked. The event changes the state between two steps; the recorder
reads the inventories after it (§16.5.3), so whatever the event takes out of
the live region without passing it to the receiver must be on the event's
record. Each deactivated cell has one `DeactivationEvent`:

| Field | Meaning | Units |
|---|---|---|
| `t` | time of the event | s |
| `invocation` | the deactivation callback's invocation it belongs to, numbered from one over the solve and shared by every cell one invocation deactivates | – |
| `cell` | the dead cell | – |
| `receiver` | the cell that received its inventories; zero when none was left | – |
| `fraction` | the dead cell's condensed mass over its reference when it went | – |
| `rows` | the state rows masked | – |
| `created_volume` | the volume the merge created | m³ |
| `lost_mass` | the masses lost, per component, $m^\mathrm{lost}_{j,e}$ | kg |
| `lost_enthalpy` | the sensible enthalpy lost, $H^\mathrm{lost}_e$ | J |
| `regularization` | the caloric regularization the merged cells held that the receiver does not, $H^\mathrm{reg}_e$ | J |

A merge of several cells into one receiver records one event per dead cell,
the outermost first; the merge's created volume, losses and regularization are
written to the first of them and the others carry zeros, so the sums over the
events are the merges'.

**The energy of a merge.** With $k$ the receiver and $\mathcal D$ the dead cells
that go into it (§14.3.3), the receiver takes the summed sensible enthalpy of
what it receives,
$H_m=\sum_{i\in\{k\}\cup\mathcal D}\sum_j m_{j,i}\,h_j(T_i)-H^\mathrm{lost}_e$,
and its merged temperature $T_m$ is the root of the extended caloric relation
(§16.4) at $H_m$ with the merged masses, so the receiver's inventory after the
merge is $H_m$. The extended inventories of the merged cells before the event
less the receiver's after it and the enthalpy lost is the regularization they
held,

$$
H^\mathrm{reg}_e=\sum_{i\in\{k\}\cup\mathcal D} H_i - H_m - H^\mathrm{lost}_e ,
$$

which is zero for cells away from empty and nonzero only for a cell far below
its capacity floor. The live energy inventory therefore falls at the event by
exactly $H^\mathrm{lost}_e+H^\mathrm{reg}_e$, and its mass inventory by
$m^\mathrm{lost}_e$.

**What is lost, by configuration.**

- *Fixed mesh.* The receiver keeps its volume and takes the dead cells'
  condensed inventories only; their gas is vented: `lost_mass` holds the
  vented gas masses and `lost_enthalpy` their sensible enthalpy at the dead
  cells' temperatures.
- *Moving mesh.* The receiver absorbs the dead cells' volume and every
  inventory: its far face moves out to enclose the merged volume, where every
  node of the dead cells then coincides. The merge never destroys pore space:
  the merged volume is the larger of the summed volumes and the merged
  skeletal volume plus the summed pre-merge pore volumes, the difference being
  `created_volume`. Nothing is lost, and only the regularization can be
  nonzero.
- *The last live cell.* It goes with no receiver (`receiver = 0`): `lost_mass`
  holds all its masses, `lost_enthalpy` their sensible enthalpy, and
  `regularization` the rest of its extended inventory. The solve then ends
  with outcome `:burned_out`.

The audit's event terms are $\sum_e m^\mathrm{lost}_e$ for the mass and
$\sum_e(H^\mathrm{lost}_e+H^\mathrm{reg}_e)$ for the energy (§16.3.3, §16.5.4); the
module test asserts that they are exactly these sums of the records. The other
records of the events are `sol.events` (all of them, in order),
`sol.extras.deactivation_times` (the time of every event),
`sol.extras.live_cells` (the live-cell count at every saved time),
`sol.statistics.invocations` (the invocations that deactivated cells) and
`sol.outcome`.

---

## 16.8 The recorder and the audit in a solve

### 16.8.1 `ConservationRecorder`

`ConservationRecorder()` is a read-only observer passed as `solve`'s
`callback`, alone or in a tuple with other callbacks (several recorders may
observe one solve). Its fields after the solve are

- `mass_outflow` [kg]: $\int\dot m_\mathrm{out}\,dt$ over the accepted steps,
- `energy_inflow` [J]: $\int\dot H_\mathrm{in}\,dt$ over the accepted steps,
- `initial`, `final`: the inventories `(energy, mass)` [J, kg] at the start and
  after the last step,
- `steps`: the steps it committed, equal to `sol.statistics.accepted`.

`solve` separates the recorders from the caller's other callbacks, binds each
to the set-up problem (resetting it, so a recorder can be reused), and
installs its reader among the observers and its commit last. The recorder
changes nothing in the integration: a solve with it and one without give
bit-identical trajectories (asserted in `test/outputs/test_conservation.jl`).
Its cost is four residual evaluations per accepted step (three quadrature
nodes and the inventories), and three more on a step inside which an event
condition crosses zero.

### 16.8.2 `conservation_audit`

`conservation_audit(sol[, recorder])` returns

```julia
(; method,                          # :recorder or :trapezoid
   mass   = (; initial, final, change, inflow, lost, defect, relative),
   energy = (; initial, final, change, inflow, lost, defect, relative))
```

with the terms of §16.3.3 and §16.5.4 in kilograms and joules. With the
recorder that observed the solve, it closes from the solve's start to its end
over the accepted steps, with every event's record.

### 16.8.3 Without a recorder: the saved series

Without a recorder the audit closes from the first to the last saved state: the
inventories are `sol.extras.live_mass` and `sol.extras.live_energy` at those
states, the inflows the trapezoid over the saved states of
`sol.extras.mass_outflow` and `sol.extras.energy_rate`, and the events those with
times in $(t_\mathrm{first},t_\mathrm{last}]$ (a state saved at an event's time is
the state after the event). The trapezoid's error is the saving interval's,
and a rate that jumps at an event makes it first order there. The saved-series
audit is a quick check of a solution in hand, not a closure at the integrator's
tolerance.

### 16.8.4 Example

The ablating solid below burns from its heated top on a fixed mesh; each
deactivation vents the dead cell's gas.

```julia
using Pyrolysis

R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
ρ, c, k, v = 1000.0, 1000.0, 1.0, 2.0e-3
A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, λ = 1e-6)
gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
material = Material(
    name = :ablator,
    components = (solid, gas),
    reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(3e-3, 12),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(2.8e6), bottom = AdiabaticBC()),
        mass = (top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0), bottom = ImpermeableBC()),
    ),
    T_initial = z -> T0 + 1e5 * z,
    ξ_initial = [z -> ρ, z -> 0.0],
    tspan = (0.0, 0.8),
    output_times = collect(0.0:0.05:0.8),
)

rec = ConservationRecorder()
sol = solve(problem; deactivation = true, callback = rec, abstol = 1e-10, reltol = 1e-8)
audit = conservation_audit(sol, rec)        # over the accepted steps
saved = conservation_audit(sol)             # trapezoid over the saved states
audit.mass.relative, audit.energy.relative
```

Run as written (1,156 accepted steps, five events of one cell each), the
recorded audit closes the mass to a relative defect of $-1.9\times10^{-8}$ and
the energy to $-3.2\times10^{-7}$. The five events vent $1.9\times10^{-4}$ kg of gas
and 101.5 J of its enthalpy, against an inventory change of 1.40 kg and
$1.77\times10^{5}$ J: without the event terms the defects would be the lost
terms, $1.4\times10^{-4}$ and $5.7\times10^{-4}$ of the change. The saved-series
audit of the same solution, at the 0.05 s saving interval, leaves relative
defects of $-1.2\times10^{-2}$ (mass) and $-5.9\times10^{-2}$ (energy).

---

## 16.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

The ledger adds an explicit, per-solve conservation audit, which none of the
three reference codes offers in this form.

**Energy accounting.** **FDS** keeps an energy budget of its simulations and
reports cumulative boundary, radiative and chemical contributions against the
stored enthalpy; Pyrolysis.jl's energy ledger is the analog for the condensed
phase, at the level of the live cells: it reconstructs the change of the
stored sensible enthalpy from the integrated boundary and source rates and the
event records and reports the closure defect.

**Mass accounting.** **Gpyro** advances control-volume conservation of
$\bar\rho\,\Delta z$ and $\bar\rho Y_i\,\Delta z$ with a fully implicit
backward-Euler scheme, so conservation is enforced by the discretization and
checked through its solid-density bookkeeping. Pyrolysis.jl is conservative in
the same sense, its cell masses being states and its face fluxes scattered
antisymmetrically, and it *audits* the result, including what deactivation
vents, with the mass ledger.

**Deformation.** ThermaKin represents condensed-phase deformation through a
single-reference-point density-integral source term and Gpyro through a
variable-$\Delta z$ control volume; Pyrolysis.jl moves the mesh with the solid
(§10, §11). The ledger needs no deformation term of its own: it sums states
over the live cells, reads the current face areas, and integrates over the live
boundary faces, which move with the solid; what a deactivation event removes
is on the event's record.

---

## 16.10 Worked summary of the ledger equations

| Ledger | Inventory | Rate | Event terms | Defect | Exactness |
|---|---|---|---|---|---|
| Mass | $m_\mathrm{live}=\sum_{i\in\mathrm{live}}\sum_j m_{j,i}$ | $\dot m_\mathrm{out}=\sum_{f\in\mathcal F_\mathrm{live}}\sum_j F_{j,f}A_f$ (positive leaving) | $\sum_e\sum_j m^\mathrm{lost}_{j,e}$ | $\Delta m_\mathrm{live}+\int\dot m_\mathrm{out}\,dt+\sum_e m^\mathrm{lost}_e$ | linear invariant in both forms |
| Energy | $H_\mathrm{live}=\sum_{i\in\mathrm{live}}H_i$ (extended caloric relation) | $\dot H_\mathrm{in}=\sum_{f\in\mathcal F_\mathrm{live}}(-F^E_fA_f)+\sum_{i\in\mathrm{live}}\big[S_{E,i}+\sum_j \omega_{j,i}\tilde V_i\,h_j(T_i)\big]$ (positive entering) | $\sum_e(H^\mathrm{lost}_e+H^\mathrm{reg}_e)$ | $\Delta H_\mathrm{live}-\int\dot H_\mathrm{in}\,dt+\sum_e(H^\mathrm{lost}_e+H^\mathrm{reg}_e)$ | linear invariant in the inventory form; at the integrator's order in the temperature form |

With the recorder the integrals are three-point Gauss–Legendre on each
accepted step's interpolant, and the defect falls with the tolerance; without
it they are the trapezoid on the saved series. The relative defect divides by
the largest of $|\text{change}|$, $|\text{inflow}|$ and $|\text{lost}|$.

Sign-convention reminders: $F_{j,f}$ and $F^E_f$ along a boundary face's outward
normal (positive leaving), so the boundary contributes $-F^E_fA_f$ to the energy
inflow; $F^E_f=q_f+\sum_j F_{j,f}h_j(T_f)$;
$Q_\mathrm{rxn}=-\sum_r h_r r_r$ with $h>0$ endothermic.

---

## 16.11 Limitations and open issues

The following are the known gaps of the ledger. They are caveats on
*interpretation*, not defects of the kernels.

1. **External sources are outside the ledger.** The volumetric mass and
   energy sources of the `sources` keyword are not in its rates; an audit of a
   solve that uses them is off by their integral.
2. **Steps it cannot account for.** A step that a callback of the caller's
   ends early, at a time no event condition the recorder knows crosses zero,
   leaves the recorder's integrals `NaN`, and the audit with them.
3. **The saved-series audit is coarse.** Without a recorder the integrals are
   trapezoids over the saving interval (§16.8.3).
4. **One aggregate per quantity.** The ledger closes the live region's total
   mass and total sensible enthalpy; it does not decompose the boundary term
   into radiation, convection and conduction, or the mass into components. The
   solution's `extras` give the exposed and back faces' flux components at the
   saved times (§16.6).
5. **Bookkeeping, not physics.** The ledger is built from the terms the
   residual assembles and cannot detect a mechanism counted twice
   (§16.5.5).
6. **Temperature-form closure.** In the temperature form the closure holds at
   the integrator's order, not to round-off (§16.5.2). In a cell close to
   empty, the floored capacity sum of its energy row departs from the slope of
   its extended inventory, which adds to the defect until the cell is
   deactivated.
7. **Relative defect of a quiet solve.** When nothing changes, enters or is
   lost, the relative defect is undefined; read the defect itself.

---

*Cross-references.* The governing equations these ledgers audit are stated in
§3 (Governing Equations); the energy row, its face fluxes and the absence of a
surface gas-enthalpy term in §7 (Heat Conduction, Gas Convection, and Energy
Assembly); the heat-of-reaction sign, its temperature datum and stoichiometric
closure in §6 (Reaction Kinetics); the gas fluxes and the boundary venting in
§9 (Gas-Phase Mass Transport and Pressure); the surface energy balance and the
boundary conditions in §12 (Boundary and Initial Conditions); the moving mesh
and its correction in §11 (The ALE Framework); the state mapping and the
assembly in §13 (Finite-Volume Discretization); the deactivation events in §14;
the integrator and its callbacks in §15. The ledgers are themselves a
verification tool used in §18 (Verification and Validation).
