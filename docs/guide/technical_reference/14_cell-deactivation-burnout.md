# 14. Cell Deactivation and Burnout

This chapter derives how Pyrolysis.jl takes a cell out of the system when its
condensed mass is nearly gone. We treat the reason a cell has to leave, the
correction and the near-threshold band that let a vanishing layer leave in one
event, the mass-fraction threshold, the masks that turn a dead cell's rows into
identity rows, the conservative transfer of its inventories to the live cell
beneath (with the pore rule on a moving mesh and the venting of its gas on a
fixed one), the block-and-band rule that decides which cells go together, the
event that locates the crossing inside a step, the outcomes a solve can end
with, the interior void, and element death on meshes of more than one
dimension. Every rule below is grounded in the implementing function
(`src/Deactivation/`, `src/Integration/setup.jl`, `src/Integration/solve.jl`);
where the code makes a choice that differs from an idealized form, we say so
and trust the code.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$ (axial, through-thickness). **$z=0$ is the
  bottom/substrate** (tag `:bottom`); **$z=L$ is the top/exposed surface** (tag
  `:top`). Heat enters from the top and cells die from the top. On a line, cell
  $i$ lies between nodes $i$ and $i+1$, which are also its faces; the pinned
  node of a moving mesh is node 1, the back face.
- $i$ is a cell index, $j$ a component index, $n$ a node index. The receiver of
  a merge is cell $k$, and $\mathcal{D}$ is the set of cells that go into it,
  ordered outermost first.
- The **condensed** components are the solid and liquid ones (every component
  that is not a gas). $m_{j,i}$ [kg] is the mass of component $j$ in cell $i$,
  a state of the integrator; $m^s_i=\sum_{j\ \rm condensed}m_{j,i}$ is the
  cell's condensed mass and $m^s_i(0)$ its initial value.
- $h_j(T)=\int_{T_{\rm ref}}^{T}c_{p,j}\,dT'$ is the sensible enthalpy of
  component $j$ (§2), $T_{\rm ref}=298.15$ K, and
  $H_i=\sum_j m_{j,i}\,h_j(T_i)$ is the sensible-enthalpy inventory of
  cell $i$.
- $\varepsilon_m$ is the deactivation threshold and $\varepsilon_{\rm band}$ the
  near-threshold band, both fractions of a cell's own initial condensed mass.
  The subscripts distinguish them from the emissivity $\varepsilon$.

---

## 14.1 Why a cell leaves the system

The integrator carries per-cell component masses and an energy state (§3,
§13). Because the masses are the states, mass is conserved exactly and a
gasifying cell's condensed mass falls smoothly toward zero; nothing in the
residual stops it there. What happens to the cell on the way depends on the
mesh:

- **Moving mesh** (`use_ale = true`). The balance rows hold the node velocities
  to every cell's volume balance (§11), so a cell whose solid gasifies into a
  gas that takes no matrix volume (swelling factor $\gamma=0$, §10) shrinks
  with its condensed volume: $\Delta z_i\to 0$.
- **Fixed mesh.** The cell keeps its volume and empties; its concentrations
  fall instead of its thickness.

Either way the cell's heat capacity $C_i=\sum_j m_{j,i}c_{p,j}(T_i)$ falls
toward that of its remaining gas. Its temperature row divides the net energy
rate by $C_i$, and its exchange with the neighbouring cells runs through
conductances of order $kA/\Delta z_i$, so the row's time scale is of order
$C_i\Delta z_i/(kA)$. On a moving mesh $C_i\propto\Delta z_i$ and the time
scale falls as $\Delta z_i^2$; on a fixed mesh it falls with $C_i$. Long before
the cell is empty it carries the fastest time scale of the system while its
contents no longer matter to the rest of the column.

The residual stays finite at every such state through its guards (§13): the
volume a concentration divides by is smoothly floored at $10^{-8}$ of the cell's
initial volume, the heat-capacity sum at $C_{{\rm floor},i}=\rho c_{p,{\rm floor}}V_i(0)$ with
`RHO_CP_FLOOR` $=\rho c_{p,{\rm floor}}=10^{3}$ J m⁻³ K⁻¹, and the masses enter the
property laws through a one-sided positive part below the cell's mass
tolerance. These guards keep a trial state evaluable; they are not physics,
and a cell deactivates long before it reaches them.

**Deactivation** takes a cell out of the system once its condensed mass falls
below the fraction $\varepsilon_m$ of its own initial condensed mass: its
inventories go to the live cell inward of it, its rows become identity rows,
and the boundary conditions move to the face of the outermost live cell. The
state length and the Jacobian's sparsity pattern never change, so the
integrator continues on the same vector and the same matrices. Deactivation is
off by default and is turned on by `solve(problem; deactivation = true)`, with
the keywords `threshold` ($\varepsilon_m$, default $10^{-2}$), `band`
($\varepsilon_{\rm band}$) and `located_events` (default `true`).

---

## 14.2 The correction and the band: one invocation per layer

Under surface heating the cells of a layer reach any mass threshold one after
another, the exposed cell first. Deactivating them one at a time means one
event per cell, each a discrete rewrite of the state followed by a restart of
the integrator. Two mechanisms let a layer that burns out leave in as few
events as possible:

1. **The equidistribution correction** (`EquidistributionCorrection`, §11)
   moves the nodes of a vanishing layer relative to the solid, $w=v_s+c$, so
   that once its first cell activates the layer's cells thin together. The
   inventories cross the moving faces at the relative velocity $-c$ (the
   relative-flux kernel, §11), and the balance rows hold the solid velocity
   $w-c$ to the volume balance.
2. **The near-threshold band** $\varepsilon_{\rm band}$: when a cell crosses
   $\varepsilon_m$, the cells of the contiguous run from the exposed side
   inward that hold less than $\varepsilon_{\rm band}$ of their initial
   condensed mass go with it, in one merge (§14.3.4).

A callback invocation that deactivates cells is counted once
(`sol.statistics.invocations`, and each event's `invocation` field), however
many cells it takes.

### 14.2.1 Where the correction is installed

When `solve` (through `setup`) is given no `correction`, it installs
`EquidistributionCorrection()` if and only if the mesh moves, cells are
deactivated, and deactivation can remove some cell; otherwise it uses
`NoCorrection()` (`_installed_correction`, `src/Integration/setup.jl`). A
correction passed explicitly is used as given; `correction = NoCorrection()`
keeps the mesh with the solid.

Whether a cell can be removed follows from the material's **residues**
(`residues(material)`, `src/Closures/correction.jl`): for each condensed
component, the condensed mass it leaves, as a fraction of its own mass, once
its reactions have converted it completely. A reaction that consumes component
$j$ leaves the condensed yields of its products, each weighted by that
product's own residue, renormalized when the reaction yields part of $j$ back;
competing reactions are bounded by the least gasifying one, so the residue is
the largest over the reactions that consume $j$:

$$
\mathrm{res}_j=\max_{\text{reactions consuming } j}\;
\frac{\sum_{p\neq j}y_p\,\mathrm{res}_p}{1-y_j},
$$

with $y_p$ the reaction's yield of product $p$ ($\mathrm{res}_p=0$ for a gas)
and $y_j$ the share it returns to $j$. The residues are the greatest fixed
point of this relation reached from $\mathrm{res}=1$: exact along a chain of
reactions, geometric around a cycle that loses mass to the gas. A component
that every consuming reaction gasifies has residue zero; an inert component, or
one no reaction consumes, has residue one. A cell **vanishes** when it starts
with condensed mass and

$$
\sum_{j\ \rm condensed}\mathrm{res}_j\,m_{j,i}(0)<\varepsilon_m\,m^s_i(0)
$$

(`vanishes`). A charring material whose char yield exceeds $\varepsilon_m$
therefore never gets the correction, and its cells never reach the threshold.

The correction acts on **layers**, bounded by the pinned node, the Lagrangian
nodes and the end nodes. Under an equidistribution correction the Lagrangian
nodes default to the interior nodes between cells that start with different
condensed components (the layer interfaces), together with the mesh's own;
the `lagrangian` keyword names others, and `setup` refuses an equidistribution
correction that would move the mesh across a layer interface not among them.
With `layers = :vanishing` (the default) only the layers deactivation can
remove move; a substrate keeps its mesh with the solid, its correction
velocity exactly zero (V-BO-08).

### 14.2.2 The layer time

The correction relaxes a layer over its equidistribution time $\tau_{\rm eq}$.
By default each layer sets it once, at its first activation
(`ActivationTime()`): when a cell of the layer first falls below the activation
threshold $r_{\rm hi}$ of its initial thickness, the layer's thickness $L_{\rm layer}$
and thinning rate $\dot L_{\rm layer}$ give

$$
\tau_{\rm eq}=\lambda_\tau\,\frac{L_{\rm layer}}{|\dot L_{\rm layer}|},\qquad \lambda_\tau=1.8
\ (\texttt{ACTIVATION\_TIME\_FACTOR}),
$$

where $\lambda_\tau$ is the `λ` keyword of `ActivationTime`. Since $\dot L_{\rm layer}$
scales with the physical rates and $\lambda_\tau$ is a pure number, a deck
whose rates are all multiplied by $s$ gets $\tau_{\rm eq}/s$. Before the
estimate the layer's correction is zero. A layer thinning no faster than its
initial thickness over the solve's horizon keeps the fallback time
(`EQUIDISTRIBUTION_TIME` $=720$ s) until a cell crosses $r_{\rm lo}$, where the
estimate is taken once more; a thickening layer keeps the fallback. The
activation thresholds default to $(r_{\rm hi},r_{\rm lo})=(0.3,0.1)$.

The estimate is taken by an event of its own: a continuous callback on the
smallest margin $\Delta z_i/(r\,\Delta z_i(0))-1$ over the live cells of the layers still
waiting ($r=r_{\rm hi}$ before the first activation, $r_{\rm lo}$ while a slow
layer waits for the second), located inside the step and widened by the
relative tolerance as the deactivation event is (§14.3.5). It does not modify
the state; the layer's correction is zero up to its activation, so the event
changes the right-hand side only where the activation has just begun. The
callback is present on every moving mesh, with an infinite margin when no
layer waits, so the integrator's type does not depend on the correction
(§15.5.7). The correction's node formula, its smooth activation switch and its
monitor are derived in §11.

### 14.2.3 The band

The band defaults to $\max(\varepsilon_m,\,0.1)$ (`CORRECTION_BAND` $=0.1$) when
the bound correction moves some layer, and to $\varepsilon_m$ otherwise, in
which case only the cells below the threshold go. A given `band` must satisfy
$\varepsilon_m\le\varepsilon_{\rm band}<1$, or `setup` throws an
`ArgumentError`. The band is read on a moving mesh only.

The band must exceed what the slowest cell of a layer still holds when the
first reaches $\varepsilon_m$. A narrower band (or none) stops the run at a cell
above it; the crossing cell below the threshold then lies inward of a live one,
and the solve ends with outcome `:interior` (§14.6). On the five burnout decks
of the supported envelope that fraction is 0.4 to 3.9 % at
$\varepsilon_m=3\times10^{-3}$ and 1.4 to 9.9 % at $\varepsilon_m=10^{-2}$; on
the mixed-composition layer of V-BO-07 the other cells still hold up to about
6 % of their initial mass when the first crosses $\varepsilon_m$.

### 14.2.4 The supported envelope

The burnout path (the equidistribution correction with its layer time set at
activation, the near-threshold band, located deactivation events) is verified
under joint refinement of the mesh and the tolerance on five decks built on
the 5.9 mm PMMA layer on Kaowool under 50 kW/m² of V-BO-01: a uniform mesh, a
mesh graded threefold toward the surface (V-BO-06), a layer with a second,
slower solid mixed in (V-BO-07), a 1 mm layer (V-BO-08), and a layer with a
convective loss in place of re-radiation and an adiabatic back. The numerical
parameters are $\varepsilon_m$, $\varepsilon_{\rm band}$, $(r_{\rm hi},r_{\rm
lo})$, the layer time, the tolerance and the cells across the layer; the
material parameters are those of the decks. Inside the envelope every quantity
of interest (complete burnout, the mass lost, the peak mass-loss rate and its
time, the burnout time, the mass ledger) is within tolerance at every
resolution and the burst at the end of a layer converges under refinement:
$\varepsilon_m=3\times10^{-3}$ with $\varepsilon_{\rm band}=0.1$ on all five
decks from 30 cells across the layer (10 on the 1 mm layer), and the default
$\varepsilon_m=10^{-2}$ with $\varepsilon_{\rm band}=0.1$ on the single-solid
layers from 30 cells (60 on the convective-loss deck) but not on the mixed
layer, whose burst stops converging between 120 and 240 cells. The number of
invocations is reported, not a criterion: one on every deck inside the
envelope except the convective-loss deck at $10^{-2}$, which takes two from 60
cells. The full statement, with its limits, is User Guide chapter 10.

---

## 14.3 Deactivation

![Surface-connected burnout transfers inventory to an inward receiver while retaining inactive state slots; insets show fixed-mesh gas venting and interior-void termination.](../../assets/figures/05_deactivation.svg)

**Figure 5. Burnout and retained state (schematic).** Here $r_i$ is the remaining condensed-mass ratio. A threshold crossing selects a surface-connected band; on a moving mesh its inventories and volume transfer to the inward receiver. Dead nodes coincide at the merged face, while state slots and topology remain allocated. The insets show fixed-mesh gas venting and an interior void blocked by live material. A callback can re-test the receiver and emit multiple per-cell records; gas venting and a final cell without a receiver enter the loss ledger.

The deactivation state of a solve is a `CellDeactivation{P,NC,E}`
(`src/Deactivation/events.jl`), built by `setup` when `deactivation = true`:
the policy (`ConservativeTransfer()`), the threshold $\varepsilon_m$, the band
$\varepsilon_{\rm band}$, the reference masses $m^s_i(0)$, the condensed
flags, the caloric extension that merged temperatures are found with, the
events so far, and the outcome. `integrate` resets it at the start of every
integration (`reset_cells!`: every row live, no events, outcome `:running`),
so repeated integrations of one set-up problem are bit-identical.

### 14.3.1 Live rows and dead rows: the mask

The state layout carries a mask over its rows, `layout.mask::BitVector`
(§15.1): `true` on a live row, `false` on a dead one. A cell is live while its
energy row is; the state mapping copies that into `ds.active` at every residual
evaluation (`derive!`). Deactivating cell $i$ masks every row of its cell block
(masses, energy state, history pair) and every row of its far node, the node
away from the receiver (position, velocity and, under `BEER_LAMBERT`, the face
intensity).

The row assembly writes a dead row as (`_dead_rows!`,
`src/Assembly/residual.jl`)

$$
\frac{du_r}{dt}=0\ \ \text{(differential row)},\qquad
0=u_r\ \ \text{(algebraic row)} .
$$

In the integrator's matrix $W=J-\mathsf{M}/(\gamma_{\rm RK}\Delta t)$ (§15.5.2)
the first is the row $-e_r^{\mathsf T}/(\gamma_{\rm RK}\Delta t)$ and the
second the row $e_r^{\mathsf T}$: each is a scaled identity row. A dead
differential unknown keeps its value and contributes no error; a dead algebraic
unknown is held at zero. The kernels read no dead cell either: an interior face
contributes only when both of its cells are live, every boundary condition
applies on the face that bounds the live cells on its side (`live_face`,
`src/Meshes/topology.jl`), the in-depth beam enters through the upstream face
of the first live cell, the balance row of a dead cell carries no velocity, the
correction is zero at every node beside a dead cell, and the cross-section rows
sum the live cells' history pairs only. The live rows therefore do not depend on
the dead unknowns, and the system splits into the live problem and a block of
identity rows.

What the dead rows hold after an event is set by the transfer (§14.3.3): zero
masses; the merged temperature in the temperature form, zero in the enthalpy
form; on a moving mesh node positions collapsed onto the receiver's far face
with zero velocities; face intensities of zero, written by the consistent sweep
that follows every event. The outputs read the mask of each saved time:
`sol.T`, `sol.ξ` and `sol.cell_centers` are `NaN` in a dead cell, `sol.m` is
zero there, and node positions beyond the outermost live node are `NaN`
(§14.6).

### 14.3.2 The threshold

A live cell is **below the threshold** when its condensed mass is less than the
fraction $\varepsilon_m$ of its own initial condensed mass:

$$
m^s_i<\varepsilon_m\,m^s_i(0)\,(1+s),
$$

with the slack $s=0$ except at a located event (§14.3.5)
(`below_threshold`). Three properties of this criterion are deliberate.

- **Mass, not thickness.** Masses are the states, the criterion reads them
  directly, and it applies unchanged on a fixed mesh, where a cell does not
  thin. With the correction active, the thickness of a vanishing layer's cells
  is controlled relative to the layer, not absolutely.
- **Each cell against its own reference.** The reference $m^s_i(0)$ is
  the cell's own initial condensed mass. On a graded mesh a single column-wide
  reference fails: with cell sizes spanning a ratio above $1/\varepsilon_m$, the
  fine cells would start below it or the coarse ones never reach it (the test
  uses a ratio of 27 against $1/\varepsilon_m=20$). A receiver keeps its own
  reference, so a merged cell deactivates once its whole merged inventory is
  below $\varepsilon_m$ of what the receiver started with.
- **A cell that starts with no condensed mass is never tested.** A gas-only cell
  has $m^s_i(0)=0$ and stays live whatever its contents.

The **threshold margin** is the smallest relative excess over the threshold,

$$
g(u)=\min_{i\ \text{live},\ m^s_i(0)>0}
\frac{m^s_i-\varepsilon_m m^s_i(0)}{\varepsilon_m m^s_i(0)},
$$

`Inf` when no live cell started with condensed mass (`threshold_margin`). It is
zero exactly when some cell holds exactly its threshold and negative exactly
when one is below it; its first root in a step is the located event's time.

The threshold also sets the mass rows' absolute tolerances: every component
mass of cell $i$ carries $\mathrm{atol}=\mathrm{rtol}\,\varepsilon_m\,m^{\rm
c}_{i,0}$ (§15.5.6; for a cell that starts with no condensed mass its total
initial mass, or 1 kg m⁻³ times its volume, stands in for $m^s_i(0)$),
so a cell at its threshold still resolves its own condensed mass to the
relative tolerance.

### 14.3.3 The conservative transfer

`ConservativeTransfer` moves the inventories of the dead cells $\mathcal{D}$ to
the receiver $k$, the live cell inward of them (`_merge!`). It conserves mass
and sensible enthalpy and finds the merged temperature from them; on a moving
mesh the receiver also absorbs the dead cells' volume, and on a fixed mesh the
dead cells' gas is vented.

**Masses.** On a moving mesh every component is transferred:

$$
\hat m_{j,k}=m_{j,k}+\sum_{d\in\mathcal D}m_{j,d}\qquad(\text{every }j).
$$

On a fixed mesh the receiver keeps its volume, and its pores could hold the
dead cells' gas only at a higher pressure, so only the condensed components are
transferred and the gas is vented:

$$
\hat m_{j,k}=m_{j,k}+\sum_{d\in\mathcal D}m_{j,d}\ \ (j\ \text{condensed}),\qquad
\hat m_{j,k}=m_{j,k},\quad m^{\rm lost}_{j,e}=\sum_{d\in\mathcal D}m_{j,d}\ \ (j\ \text{gas}),
$$

the vented masses $m^{\rm lost}_{j,e}$ written to the event as lost.

**Sensible enthalpy and the merged temperature.** The receiver takes the summed
sensible enthalpy of what it receives, each cell's share at its own temperature,

$$
H_m=\sum_j m_{j,k}\,h_j(T_k)+\sum_{d\in\mathcal D}\sum_j (m_{j,d}-m^{\rm vent}_{j,d})\,h_j(T_d),
\qquad
H^{\rm lost}_e=\sum_{d\in\mathcal D}\sum_{j\ \rm gas}m^{\rm vent}_{j,d}\,h_j(T_d),
$$

with $m^{\rm vent}_{j,d}=0$ on a moving mesh, and the merged temperature $T_m$ is the root
of the receiver's extended caloric relation at the merged masses,

$$
\mathcal{E}_k(T_m;\,\hat m)=H_m,\qquad
\mathcal{E}_k(T;m)=\sum_j (m_j)^{+}\,h_j(T)+C_{{\rm reg},k}(m)\,(T-T_{\rm ref}),
$$

solved by the bracketed Newton iteration of `invert_caloric` (§13). Here
$(m_j)^{+}$ is the one-sided positive part of $m_j$ below the cell's mass
tolerance and $C_{{\rm reg},k}$ the capacity switch of the caloric extension,
zero while the cell's guaranteed heat capacity is at least twice its floor
$C_{{\rm floor},k}$. On every state the tolerances resolve (each component
absent or above its mass tolerance, the receiver above twice its floor) the
relation is the physical one, $\sum_j \hat m_{j,k}\,h_j(T_m)=H_m$.
The receiver's energy state becomes $T_m$ in the temperature form and
$\mathcal{E}_k(T_m;\hat m)$ in the enthalpy form.

Two choices matter. First, the enthalpy is the integral $\int c_p\,dT$, not the
product $c_p(T)\,T$, which for a temperature-dependent heat capacity is not an
enthalpy; the integral is in closed form for constant, linear, polynomial and
tabulated heat-capacity laws, and from a cubic Hermite interpolant of its values
on a 1 K grid for any other law (`sensible_enthalpy`). Second, every component
the receiver takes is in the balance, gas included, the same heat capacity the
energy row stores (§5, §7). For constant heat capacities the merged
temperature reduces to the heat-capacity-weighted mean,
$T_m=\sum m c_p T/\sum m c_p$ over what the receiver holds: merging a cell
holding 2000 kJ K⁻¹ at 800 K with one holding 500 kJ K⁻¹ at 400 K gives 720 K,
where the volume-weighted mean would give 600 K.

**The caloric regularization.** The extended energies the cells held before the
event may exceed their sensible enthalpies by the regularization
$C_{\rm reg}(T-T_{\rm ref})$, which only a cell far below its heat-capacity
floor carries. It is not transferred; the merge writes to the event

$$
H^{\rm reg}_e=\sum_{c\in\{k\}\cup\mathcal D}\mathcal{E}_c^{\rm before}-\mathcal{E}_k(T_m;\hat m)-H^{\rm lost}_e,
$$

which is zero to round-off unless a cell sat below its floor. In the
temperature form $\mathcal{E}_c^{\rm before}$ is the extended energy of the
cell's masses at its temperature; in the enthalpy form it is the stored
inventory.

**The history pair.** Under a variable cross section each cell carries the
conversion-history pair (the cumulative dry gas it produced and its initial dry
solid mass, §10). Both are extensive and are summed into the receiver; the dead
cells' pairs are zeroed. The column's conversion, the ratio of the two sums, is
unchanged by the event, and so is the cross section.

**Volume and the pore rule (moving mesh).** The receiver encloses the dead
cells' volume: its far node moves out to the position that bounds the merged
volume $\hat V_k$ from its inner node (`_enclose!`; on a planar mesh
$z_{\rm far}=z_{\rm inner}+\hat V_k/A$, on a radial mesh the corresponding shell),
and every far node of the dead cells is placed on it with zero velocity, the
moving-boundary collapse. A merge never destroys pore space. With
$V_{\rm s}(m,T)=\sum_{j\ \rm condensed}m_j/\rho_{i,j}(T)$ the skeletal volume of
the masses $m$ (intrinsic densities, §4; a component without one falls back to
its bulk density and contributes no pore volume),

$$
\hat V_k=\max\Bigl(V_k+\sum_{d\in\mathcal D}V_d,\;
V_{\rm s}(\hat m_{\cdot,k},T_m)+\sum_{c\in\{k\}\cup\mathcal D}\bigl[V_c-V_{\rm s}(m_{\cdot,c},T_c)\bigr]\Bigr),
$$

the larger of the summed volumes and the merged skeletal volume plus the summed
pre-merge pore volumes. The excess over the summed volumes,
$\hat V_k-V_k-\sum_d V_d\ge 0$, is written to the event as created volume. The
rule exists because the summed volume alone can be inadmissible: two unit cells
of a planar column with unit area, one holding 1 kg of a solid whose intrinsic
density is $300/T$ kg m⁻³ at 300 K, the other 1 kg of a solid of intrinsic
density 1 kg m⁻³ at 500 K, merge at equal heat capacities to 400 K, where their
skeletal volume is $4/3+1=7/3$ m³ against the summed 2 m³: a porosity of
$-1/6$. The rule gives $\hat V_k=7/3$ m³ and a created volume of $1/3$ m³, the
skeletal thermal expansion at the merged temperature. When the cells hold gas,
the retained gas keeps at least the pore volume it had, at a finite pressure
that Darcy transport then relaxes. On a fixed mesh the receiver keeps its
volume and no volume is created.

**Conservation.** Mass is conserved to round-off with the vented mass, and
sensible enthalpy with the vented enthalpy and the regularization; the tests
check $\sum m+\text{lost}$ to $10^{-14}$ and the energy ledger to $10^{-12}$ in
both energy forms, on both meshes, for single merges, runs of cells, a run
through the innermost cell and a run through a layer's innermost cell into the
substrate (`test/deactivation/test_deactivation.jl`).

### 14.3.4 The deactivation pass: the block-and-band rule

`deactivate!(u, d, def, t; slack)` deactivates, in one pass, the cells below the
threshold at the exposed side, in place on the state $u$. Every event of one
call shares one invocation number. The pass repeats:

1. Let `top` be the outermost live cell (`top_live_cell`); if none is left,
   stop.
2. **Moving mesh.** Find the run of live cells from `top` inward whose condensed
   mass is below $\varepsilon_{\rm band}(1+s)$ of their reference
   (`_run_start`); let `lo` be its inward end. **Fixed mesh.** The run is the
   outermost live cell alone, `lo = top`.
3. If no cell of the run is below $\varepsilon_m(1+s)$, stop.
4. If `top` is the innermost cell (cell 1), it goes with no receiver
   (`_lose!`): its masses, its sensible enthalpy and its regularization are
   written to the event as lost, its rows are masked, and the pass returns
   `:burned_out`.
5. Otherwise the receiver is the live cell beneath the run, $k=\max(\texttt{lo}-1,1)$,
   and the cells $\mathcal D=\{k+1,\dots,\texttt{top}\}$ merge into it in one
   transfer (§14.3.3). A run that reaches the innermost cell leaves that cell,
   holding the run's inventories, as the survivor.
6. Return to step 1: the receiver, holding the run's inventories against its own
   reference, is tested in turn, and when it starts a run below the band with a
   cell below the threshold, that run merges on in the same invocation. On a
   fixed mesh the receiver is re-tested with the dead cell's condensed
   inventories added.

After the loop the pass looks for any live cell still below the threshold. If
one remains, it lies inward of a live cell that the pass leaves, which no
policy here handles: the pass sets the outcome `:interior` and returns it.
Otherwise it returns `:transferred` (`:burned_out` from step 4). A crossing
inward of a cell above the band is therefore never merged across: on a moving
mesh the cells between it and the exposed face would have to go too, and they
hold more than the band allows.

The rule takes a contiguous block of dead cells touching the exposed boundary
and a band around the crossing cell; this is what keeps the live cells a
contiguous range $1,\dots,\texttt{top}$ ordered from the back, the receiver a
neighbour of every dead cell it takes, and the exposed face the far face of the
receiver. A positive `slack` widens both tests to $1+s$ times $\varepsilon_m$
and $\varepsilon_{\rm band}$ (§14.3.5).

### 14.3.5 The event in the integrator: location and the safety net

`integrate` installs the deactivation as two callbacks
(`_deactivation_callback`, `src/Integration/solve.jl`):

- **The located event** (`located_events = true`, the default). A continuous
  callback on the threshold margin $g(u)$, tested for downward crossings only
  (`affect_neg!`), with the root found inside the step on the integrator's
  interpolant and the state taken just past it (`RightRootFind`). At the root
  the crossing cell holds exactly its threshold, so the event runs
  `deactivate!` with the slack $s=\mathrm{rtol}$, the mass rows' relative error
  scale there: the crossing cell goes, and with it the cells of its run within
  the tolerance of their thresholds, whose crossing the solve does not resolve
  from the located one. A cell within `reltol` of its threshold inward of a
  cell above the band stays.
- **The safety net.** A discrete callback at the end of every accepted step
  whose condition is that some live cell is below the threshold (`any_below`),
  running `deactivate!` with $s=0$. Without `located_events` it is the only
  event; with it, it catches a crossing the location does not see. Location
  takes a cell at its first crossing, so a cell under a correction that would
  have recovered by the step's end goes there; only the safety net, which sees
  the step's end alone, misses a dip and recovery inside one step.

After a pass that transferred inventories, the event restores the algebraic
rows and checks the result:

1. the face intensities by the attenuation recurrence from the exposed face
   (`_consistent_intensities!`), which writes zero on the faces upstream of the
   first live cell;
2. the cross-section area and lateral rate, then the node velocities by one
   sweep of the balance rows from the pinned node (`_consistent_motion!`,
   §15.4); a velocity row left above round-off returns `:inconsistent`;
3. admissibility (`inadmissible`, §15.5.4) with every live cell, the receiver
   included, held to its own mass tolerances, which are those the domain test
   applies to every later proposal, and under Darcy–Fick transport a porosity in
   $[0,1]$ within the skeletal volume of the mass tolerances; a failure returns
   `:inadmissible`.

Each status has one consequence:

| Status | Integrator | Warning |
|---|---|---|
| `:transferred` | continues, with the policy's restart step | — |
| `:burned_out` | `terminate!` (`:Terminated`) | — |
| `:interior` | `terminate!` (`:Terminated`) | names the cell below its threshold inward of a live cell |
| `:inconsistent` | `terminate!(…, ReturnCode.Failure)` | no consistent node velocities after deactivating cells |
| `:inadmissible` | `terminate!(…, ReturnCode.Failure)` | the state after deactivating cells is inadmissible |

The restart after a transfer is the integrator policy's (§15.5.1): the default
`:heuristic` restart proposes the integrator's own initial-step estimate at the
state the event left (`auto_dt_reset!`, two right-hand-side evaluations);
`:forced_dt` proposes `DT_AFTER_EVENT` $=10^{-6}$ s; `:solver_chosen` keeps the
integrator's proposal. A transfer moves the receiver off the quasi-steady
profiles of its neighbours, so a restart from the step the integrator was
taking before the event would start far from the solution.

The located event writes no extra saved state (`save_positions = (false,
false)`); the integrator applies an event before it saves the state of the step
that triggered it, so a state saved at an event's own time is the one after the
event.

### 14.3.6 The Jacobian through deactivation

The sparsity pattern is built once, from the connectivity, for every row
(§15.11); deactivation masks rows without changing it. A dead differential row
of the Jacobian is zero, a dead algebraic row is the unit row $e_r^{\mathsf T}$,
and the entries coupling live rows to dead unknowns are zero because the
kernels read no dead cell. The coloring, the prepared differentiation and the
linear solver's structure therefore stay valid through every event. After
`deactivate!` the pattern contains the dense reference Jacobian and the colored
Jacobian equals it per entry to $10^{-12}$ (§15.12, §15.13).

---

## 14.4 Fixed-mesh deactivation: venting and the interior void

On a fixed mesh only the outermost live cell goes (§14.3.4): it merges into the
cell beneath it, which keeps its volume, takes the condensed inventories and
none of the dead cell's gas, and becomes the surface. The vent keeps the
receiver's pore gas where its own transport left it: in the test that pins
this, the receiver would hold three times its gas without the vent, and with it
its pressure moves only with the condensed volume and the temperature. The
vented masses and their sensible enthalpy are written to the event and enter
the conservation audit as losses (§16).

A fixed mesh therefore burns a layer **top down**, which requires the exposed
cell to empty first. Under surface absorption it does, since the exposed cell
is the hottest: on the PMMA burnout deck of V-BO-02 (the V-BO-01 deck on a
fixed mesh with the incident flux absorbed at the surface, 60 PMMA cells) every
PMMA cell deactivates as the top live cell, the run loses exactly the layer's
mass and ends on the bare substrate, and the maximum live-cell temperature over
every accepted step is 915 K.

Under in-depth absorption the temperature maximum of a regressing
non-charring polymer lies below the surface, so the exposed cell reacts slowest.
On a fixed mesh it keeps its volume and loses concentration: it persists as a
low-density skin while the cells beneath it empty, and the layer gasifies from
within. The first cell below its threshold then lies inward of a live cell, the
pass of §14.3.4 cannot take it, and the solve ends with return code
`:Terminated`, outcome `:interior` and a warning naming the cell. On the deck of
V-BO-05 (the V-BO-01 deck, in-depth absorption included, on a fixed mesh) this
happens before any deactivation, when the lowest PMMA cell, next to the
Kaowool, crosses $\varepsilon_m$: at 260.35, 259.79, 259.62 and 259.57 s at 30,
60, 120 and 240 cells, both energy forms within 5 ms, with the maximum live-cell
temperature below 1000 K.

The stop is deliberate. A void cell conducts with the gas conductivity and
absorbs nothing, and nothing in the model carries heat across a void by
radiation: the in-depth model absorbs and does not emit, and the only emission
is the surface balance's at the exposed face. The beam transmitted through a
void column lands on the first condensed cell below it, which can return that
heat upward only by conduction through the gas, so continued past its first
void the deck would heat the substrate beneath it far outside the property
bracket. On a moving mesh the configuration rarely arises: with a
non-swelling gas an emptying cell collapses, and the correction brings the
layer to the threshold together. The same outcome ends a moving-mesh solve
whose crossing cell lies inward of a cell above the band (§14.2.3).

---

## 14.5 Element death on meshes of more dimensions

On a conforming mesh of more than one dimension (`box_mesh`, `cell_mesh`) a
cell cannot collapse onto a neighbour without distorting the cells of the
neighbouring column, so it leaves by **element death** instead
(`element_death!(layout, topology, cells)`, `src/Deactivation/element_death.jl`):

- the cells keep their geometry and leave the system: their rows are masked, as
  are the rows of every node that no live cell keeps; a node a dead cell shares
  with a live one stays live and moves with it;
- a face between a live cell and a dead one becomes a boundary face of the live
  cell, its normal turned out of it, carrying the condition of the exposed tag
  (the thermal boundary kernel's `exposed` condition);
- the mesh closure balances the live cells alone; the potential closure's rank
  check takes the live cells (`rank_check(…; active)`);
- no node moves and no inventory is transferred: the dead cells' inventories
  leave the system;
- the state length and the Jacobian's pattern never change.

Element death is a residual-level building block. `setup` refuses
`deactivation = true` on a mesh of more than one dimension (deactivation, like
gas transport, in-depth radiation, the variable cross section and the
equidistribution correction, is a line's), so `solve` never calls it, and full
solves on such meshes are not a supported path. The test applies it to a
3 × 3 mesh, a corner cell and then an edge cell sharing its nodes, in both
energy forms: the live rows equal those of the mesh without the dead cells to
$10^{-12}$, with the faces it exposes under the exposed condition, and the unchanged
pattern contains the Jacobian (`test/multid/test_element_death.jl`).

---

## 14.6 Outcomes and event records

Every deactivated cell writes one `DeactivationEvent`:

| Field | Meaning |
|---|---|
| `t` | the event's time [s] |
| `invocation` | the callback invocation it belongs to, numbered from one over the solve, shared by every cell one invocation deactivates |
| `cell` | the dead cell |
| `receiver` | the receiver; zero when none was left (the cell's inventories are then lost) |
| `fraction` | the dead cell's condensed mass as a fraction of its reference when it went |
| `rows` | the rows masked |
| `created_volume` | the volume the merge created [m³] (§14.3.3) |
| `lost_mass` | the masses lost, per component [kg] (vented gas; a last cell's inventories) |
| `lost_enthalpy` | the sensible enthalpy lost [J] |
| `regularization` | the caloric regularization the cells held that the receiver does not [J] |

A merge of several cells into one receiver records one event per dead cell, the
outermost first; the merge's created volume, losses and regularization are
written to that first event and the others carry zeros, so the sums over the
events are the merges'. The solution carries them as `sol.events` and their
times as `sol.extras.deactivation_times`; `sol.extras.live_cells` counts the live
cells at each saved time and `sol.statistics.invocations` the invocations.

`replay_events!(d, layout, t)` reconstructs the mask of a state saved at $t$
from the events at or before $t$; the solution's reconstruction of every saved
state, `sol(t)` and `interpolate_solution` use it. Since the interpolant between
two saved states that an event separates is not the trajectory, `sol(t)` refuses
a time strictly inside such an interval, naming the event.

A solve ends in one of these ways:

| End | `sol.retcode` | `sol.outcome` | `succeeded(sol)` |
|---|---|---|---|
| the time span reached | `:Success` | `:running` | `true` |
| the last live cell went (§14.3.4, step 4) | `:Terminated` | `:burned_out` | `true` |
| the live thickness fell below `min_thickness` (moving mesh) | `:Terminated` | `:min_thickness` | `true` |
| a cell below the threshold inward of a live cell | `:Terminated` | `:interior` | `false` |
| no consistent node velocities after an event, or an inadmissible state after it | `:Failure` | `:running` | `false` |
| any other integrator failure | its return code | `:running` | `false` |

`completed(sol)` is the strict test, `:Success` only. The thickness stop is a
continuous callback on the distance from the first node to the outermost live
one, $z_{\rm top\,live}-z_1-h_{\min}$ with $h_{\min}=$ `min_thickness` (default
$10^{-6}$ m), so it ends the solve at the located crossing.

**Example.** A 1 mm slab of a solid that gasifies completely ($\rho=1190$
kg m⁻³, $c_p=1420$ J kg⁻¹ K⁻¹, $k=0.21$ W m⁻¹ K⁻¹, one first-order reaction
$A=8.5\times10^{12}$ s⁻¹, $E=188$ kJ mol⁻¹, $h=870$ kJ kg⁻¹) on 20 uniform cells,
under 100 kW/m² of incident radiation absorbed at the surface, with re-radiation
and convection at the top and an adiabatic back:

```julia
using Pyrolysis

solid = SolidComponent(:solid, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, λ = 1e-6)
gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5)
material = Material(
    name = :slab,
    components = (solid, gas),
    reactions = (Reaction(:pyrolysis, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(1e-3, 20),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 100e3) + RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    ),
    tspan = (0.0, 80.0),
    output_times = collect(0.0:10.0:80.0),
)

moving = solve(problem; use_ale = true, deactivation = true,
               radiation_model = SURFACE_ABSORPTION)
fixed = solve(problem; deactivation = true, radiation_model = SURFACE_ABSORPTION)
for sol in (moving, fixed)
    println((sol.retcode, sol.outcome, sol.statistics.invocations, length(sol.events)))
end
```

The script prints

```
(:Terminated, :min_thickness, 1, 19)
(:Terminated, :burned_out, 20, 20)
```

On the moving mesh `setup` installs the equidistribution correction (the solid
leaves no residue) with the band 0.1, and the whole layer goes in one
invocation at $t=25.12$ s: cells 20 to 2 merge into cell 1, the survivor, their
`fraction` rising from 0.0100 for the crossing cell 20 to 0.0483 for cell 2.
The survivor burns on until the live thickness reaches `min_thickness`. On the
fixed mesh the band is the threshold, the cells go one per invocation from the
top, and the last one goes with no receiver.

---

## 14.7 Interaction with the solver

The deactivation is coupled to the time integration only through its callbacks
(§15.5.3) and the read-only residual; the points where the two meet are these.

1. **The state is rewritten in place, then made consistent.** A transfer changes
   the differential unknowns of the receiver and the dead cells; the event then
   recomputes the algebraic unknowns (§15.4) and checks admissibility before the
   integrator continues, with no initialization of its own (`NoInit()`).
2. **The receiver is held to its own tolerances.** After a merge the receiver
   holds the merged cells' masses against its own absolute tolerance, which is
   what the domain test holds every later proposal to (§15.5.4); the event's
   admissibility check uses the same tolerances, so a state the event accepts is
   one the integrator can continue from.
3. **The restart is the policy's** (§14.3.5, §15.5.1).
4. **Events and saved states.** A saved state is reconstructed under the mask of
   its time (§14.6); a solve that ends at an event saves its last state.
5. **Conservation.** The conservation recorder integrates the boundary fluxes
   over each accepted step up to the time an event ends the step early (the
   threshold margin, the layer-time margin and the thickness margin are the
   conditions it watches), and the audit closes mass and energy with the events'
   losses and regularization (§16).
6. **Sensitivities.** Sensitivities are finite differences of whole solves. A
   located event's time is the first root of the threshold margin, a smooth
   function of the parameters, so the derivative exists through an event; an
   objective sampled at fixed times sees the event's burst move across its
   sample times, which the step choice of §17.2.2 accounts for.
7. **Types.** The deactivation callbacks reach the problem through
   fixed-signature function wrappers, so whether cells are deactivated is part
   of the solver configuration the integrator compiles for, while the
   material, the boundary conditions and the deck are not (§15.5.7).

---

## 14.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

The reference codes for condensed-phase pyrolysis take different stances on how
to represent a receding surface. Mapping their notation to ours on first use:

- **ThermaKin2Ds** (Stoliarov–Lyon) represents volume change by a dynamic-mesh
  adjustment driven by the deformation source
  $\int(1/\rho)(\partial\rho/\partial t)\,dx$ about a single moving reference
  point, rather than by a topology change. The Pyrolysis.jl equivalent of that
  motion is the moving mesh of §11, whose balance rows move the nodes with the
  solid. What Pyrolysis.jl adds is the deactivation of a cell once its condensed
  mass falls below $\varepsilon_m$ of its initial value, with its inventories
  transferred conservatively to the cell beneath, and the correction that lets a
  vanishing layer leave in one event.
- **Gpyro** (Lautenberger–Fernandez-Pello) uses a control-volume mesh whose
  thicknesses shrink as the solid is consumed, the front-most control volume
  allowed to vanish. Its "cell disappears" behaviour corresponds to
  deactivation here; the enthalpy-conserving transfer of §14.3.3 is the careful
  form of folding a vanishing control volume's energy into its neighbour, and
  the pore rule keeps the receiver's porosity admissible when the intrinsic
  densities depend on temperature.
- **FDS** (McGrattan et al.) solves the solid phase on a through-thickness grid
  whose surface-cell removal as the surface regresses is the analogue of the
  top-down deactivation of a fixed mesh here (§14.4), with swelling and
  shrinking applied through a density-ratio rule on the cell sizes. Pyrolysis.jl
  folds the removed cell's sensible enthalpy into its interior neighbour and
  writes the vented gas and its enthalpy to the event.

All four codes reconcile a numerical grid with a physically receding boundary,
and all four remove resolution at the surface as material is consumed.
Pyrolysis.jl's choices are a threshold on mass rather than thickness, a
conservative transfer with every loss recorded, a state length and Jacobian
pattern that never change, and an event located inside the step.

> Note: the reference comparison here is drawn from the general published
> descriptions of these codes, so the mapping above is at the level of modeling
> strategy rather than equation-by-equation correspondence. The Pyrolysis.jl side
> of every statement is grounded in the source as cited throughout this chapter.

---

## 14.9 Worked summary of the key relations

**Threshold and margin** (`src/Deactivation/events.jl`):

$$
m^s_i<\varepsilon_m\,m^s_i(0)(1+s),\qquad
g(u)=\min_i\frac{m^s_i-\varepsilon_m m^s_i(0)}{\varepsilon_m m^s_i(0)},
\qquad \varepsilon_m=10^{-2}\ \text{by default},\quad s=\mathrm{rtol}\ \text{at a located event}.
$$

**Band** (`src/Integration/setup.jl`): $\varepsilon_{\rm band}=\max(\varepsilon_m,0.1)$
when the bound correction moves some layer, $\varepsilon_m$ otherwise;
$\varepsilon_m\le\varepsilon_{\rm band}<1$.

**Installation of the correction**: `use_ale`, `deactivation` and some cell with
$\sum_j\mathrm{res}_j m_{j,i}(0)<\varepsilon_m m^s_i(0)$.

**Layer time**: $\tau_{\rm eq}=\lambda_\tau L_{\rm layer}/|\dot L_{\rm layer}|$, $\lambda_\tau=1.8$, at the
layer's first activation ($\Delta z_i<0.3\,\Delta z_i(0)$), fallback 720 s.

**Transfer** (`ConservativeTransfer`):

$$
\hat m_{j,k}=m_{j,k}+\sum_{d\in\mathcal D}(m_{j,d}-m^{\rm vent}_{j,d}),\qquad
\mathcal{E}_k(T_m;\hat m)=\sum_j m_{j,k}h_j(T_k)+\sum_{d\in\mathcal D}\sum_j(m_{j,d}-m^{\rm vent}_{j,d})h_j(T_d),
$$

with $m^{\rm vent}_{j,d}$ the gas masses of the dead cells on a fixed mesh and zero on a
moving one; on a moving mesh

$$
\hat V_k=\max\Bigl(\sum_{c}V_c,\;V_{\rm s}(\hat m,T_m)+\sum_{c}\bigl[V_c-V_{\rm s}(m_c,T_c)\bigr]\Bigr),
\qquad c\in\{k\}\cup\mathcal D .
$$

**Dead rows** (`src/Assembly/residual.jl`): $du_r/dt=0$ on a differential row,
$0=u_r$ on an algebraic row; identity rows of $W$.

**Pass**: moving mesh, the run from the outermost live cell inward below
$\varepsilon_{\rm band}$, if it holds a cell below $\varepsilon_m$, into the cell
beneath it, repeated on the receiver; fixed mesh, the outermost live cell while
below $\varepsilon_m$; the last cell lost; a remaining crossing inward of a live
cell ends the solve with `:interior`.

---

## 14.10 Limitations and open issues

1. **Interior voids.** The model carries no conduction or radiation model for a
   gas gap, so a cell below the threshold inward of a live cell ends the solve
   with outcome `:interior` (§14.4). On a fixed mesh under in-depth absorption a
   non-charring polymer reaches this before its first deactivation.
2. **One invocation per layer is not guaranteed.** The correction and the band
   bring a layer to the threshold together on the decks of the envelope
   (§14.2.4), but graded meshes, heterogeneous layers and non-uniform reaction
   work against it; the number of invocations is reported, and a band narrower
   than what the slowest cell keeps ends the solve with `:interior`.
3. **The burst after a layer's deactivation.** At 30 cells across the layer the
   largest per-step rate just after a layer's deactivation can exceed the main
   peak, up to threefold on the graded mesh; it falls below the main peak by 60
   to 120 cells, and the peak sampled at the output times is within tolerance
   from 30 cells. On the mixed layer of V-BO-07 at the default
   $\varepsilon_m=10^{-2}$ the burst stops converging between 120 and 240 cells;
   $\varepsilon_m=3\times10^{-3}$ converges there.
4. **The envelope's extent.** Materials outside the five decks (other materials,
   thicker or charring layers, higher heating) are not inside the envelope by
   that study; the V-BO-01 ladder at 25 and 100 kW/m², at half and twice its
   pre-exponential factor and with its flux stepped after activation, at 30 to
   120 cells, passes, the hot decks with two invocations.
5. **Element death is residual-level.** Deactivation on a mesh of more than one
   dimension is refused by `setup`; element death is tested at the residual and
   Jacobian level only (§14.5).
6. **Fixed-mesh losses.** On a fixed mesh the dead cells' gas is vented to the
   event, not to the boundary flux; the conservation audit closes with it as a
   loss (§16).
7. **Sensitivities across events.** An objective sampled at fixed times is rough
   on the scale of one sample interval in the event time; the finite-difference
   step must exclude such steps (§17.2.2).
