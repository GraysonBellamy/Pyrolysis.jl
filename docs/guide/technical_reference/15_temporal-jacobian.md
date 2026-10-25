# 15. Temporal Integration, the DAE, and the Colored Jacobian

This chapter describes how the finite-volume rows of §13 become a system the
integrator advances: the state vector and its layout, the immutable definition
and the buffers one residual evaluation writes, the residual as a pure function
of the state, the algebraic rows and the consistent state they start from, the
mass-matrix DAE and the integrator policies that advance it, the domain test
that keeps accepted states physical, the per-row tolerances, and the Jacobian:
colored sparse forward-mode automatic differentiation of the whole residual on
a pattern built from the mesh connectivity, with the linear solve of the
integrator's matrix and the dense reference the pattern and the values are
tested against. Every rule is grounded in the implementing function
(`src/States/`, `src/Assembly/`, `src/Integration/`, `src/ColoredJacobian/`);
where the code makes an approximation, a guard or a sign choice, we state it and
trust the code over any idealized form.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$. **$z=0$ is the bottom/substrate** (tag
  `:bottom`); **$z=L$ is the top/exposed surface** (tag `:top`). On a line, cell
  $i$ lies between nodes $i$ and $i+1$, which are also its faces.
- A face flux $F_f$ is the transport per unit area along the face normal, which
  points from the face's owner cell to its neighbour ($+z$ on every interior
  face and on `:top`, out of the domain, $-z$, on `:bottom`); it is scattered as
  $\mp F_fA_f$ into the owner and the neighbour, and $s_{if}=\pm1$ is the
  orientation of face $f$ relative to cell $i$ ($+1$ when $i$ owns it).
- Heat of reaction: $h>0$ endothermic, $Q_{\rm rxn}=-\sum_r h_r r_r$ (overload
  resolution H1).
- $i$ is a cell index, $j$ a component index, $n$ a node index, $f$ a face
  index; $r$ indexes a row of the state vector where stated. $N_C$ is the number
  of components and $n$ the number of cells where it is a count.
- $\boldsymbol u$ is the state vector, $\boldsymbol f(\boldsymbol u,t)$ the
  residual (the right-hand side the residual function writes into
  $\boldsymbol{du}$), $\mathsf M$ the diagonal mass matrix and
  $J=\partial\boldsymbol f/\partial\boldsymbol u$ the Jacobian. The sans-serif
  $\mathsf M$ distinguishes the mass matrix from a molar mass $M_j$.

---

## 15.1 The state vector: layout and row kinds

![Interleaved node and cell state strip and three actual six-cell Jacobian patterns, distinguishing algebraic rows, local coupling, wider correction reach and global area coupling.](../../assets/figures/15_state_jacobian.svg)

**Figure 15. State layout and Jacobian structure.** The schematic strip shows the interleaved ordering of a six-cell, two-component, variable-area setup. The matrices are actual stored patterns, with Beer–Lambert intensity unknowns. Grey entries are allocated dependencies; blue entries are nonzero initial-state AD derivatives. Correction reserves a wider stencil even while initially inactive. Variable area adds history and the two global algebraic rows outlined in purple. Deactivation preserves this allocation. The CSV records each stored entry and its initial derivative.

### 15.1.1 Interleaved blocks

The state is described by a `StateLayout{NC,E,Mov,Rad,Chi,D,U,Car}`
(`src/States/layout.jl`): a sequence of node blocks and cell blocks, then the
global scalars. On a line the order is

$$
\boldsymbol u=\bigl[\,\text{node }1,\ \text{cell }1,\ \text{node }2,\ \text{cell }2,\ \dots,\ \text{cell }n,\ \text{node }n{+}1,\ \text{globals}\,\bigr],
$$

so every row's unknowns and those of its neighbours lie within a fixed distance
of it and the Jacobian is banded in every geometry. In more dimensions each
cell's block follows the blocks of its nodes not yet placed, cell by cell
(`block_order`).

| Block | Unknowns, in order | Present |
|---|---|---|
| node $n$ | position $z_n$; velocity $w_n$; intensity $I_n$ of the face at node $n$ | position and velocity on a moving mesh (`use_ale`, the velocities when the closure's unknowns are the node velocities); intensity under `BEER_LAMBERT` (a line) |
| cell $i$ | masses $m_{1,i},\dots,m_{N_C,i}$ [kg]; the energy state; the carried fields; the closure's multiplier | always; the energy state is $T_i$ [K] (`energy_state = :temperature`) or the enthalpy inventory $H_i$ [J] (`:enthalpy`); carried fields: the conversion-history pair (under a variable cross section) and, in more dimensions, the deformation gradient; the multiplier under the potential closure of more dimensions |
| globals | cross-section area $A$; its lateral strain rate $\theta_A$ | under a variable cross section (a line) |

With $N_{\rm node}=2\delta_{\rm mov}+\delta_{\rm BL}$ the width of a node block and
$N_{\rm cell}=N_C+1+2\delta_\chi$ that of a cell block on a line ($\delta_{\rm mov}$,
$\delta_{\rm BL}$, $\delta_\chi$ equal to one on a moving mesh, under in-depth
radiation and under a variable cross section, zero otherwise, as in §3.1.2), and
$N_{\rm blk}=N_{\rm node}+N_{\rm cell}$ the stride, the indices are

$$
\begin{aligned}
z_n&: (n-1)N_{\rm blk}+1, & w_n&: (n-1)N_{\rm blk}+2, & I_n&: (n-1)N_{\rm blk}+2\delta_{\rm mov}+1,\\
m_{j,i}&: (i-1)N_{\rm blk}+N_{\rm node}+j, & \text{energy}_i&: (i-1)N_{\rm blk}+N_{\rm node}+N_C+1, & \text{history}_{i,1..2}&: (i-1)N_{\rm blk}+N_{\rm node}+N_C+1+(1,2),\\
A&: nN_{\rm blk}+N_{\rm node}+1, & \theta_A&: nN_{\rm blk}+N_{\rm node}+2, & N_{\rm state}&=nN_{\rm blk}+N_{\rm node}+2\delta_\chi,
\end{aligned}
$$

which `test/multid/test_layout.jl` checks against every index function. Code
never computes an index: kernels, the closure, the deactivation and the outputs
go through `mass_index`, `energy_index`, `history_index`, `position_index`,
`velocity_index`, `intensity_index`, `area_index`, `lateral_rate_index`,
`cell_rows` and `node_rows`.

For two components on 12 cells the state has 36 rows on a fixed mesh, 49 under
in-depth radiation, 62 on a moving mesh and 75 on a moving mesh under in-depth
radiation (measured with the script of §15.13).

### 15.1.2 Row kinds and the mask

Every row is **differential** (masses, energy state, node positions, carried
fields) or **algebraic** (node velocities, multipliers, face intensities, the
area and the lateral rate); `is_differential(layout)` flags them and
`mass_matrix(layout)` is the diagonal matrix with one on the differential rows
and zero on the algebraic ones (§15.5.1).

`layout.mask` marks the live rows. Every row is live at the start of an
integration; deactivation (§14.3.1) and element death (§14.5) clear rows, and
the assembly writes each dead row as an identity row of the system: a zero
derivative on a differential row, $0=u_r$ on an algebraic one. The state length
never changes during a solve.

---

## 15.2 The definition and the workspace

The residual is a pure function of the state and the time, parameterized by an
immutable definition and a set of buffers:

```julia
residual!(du, u, def::CoreDefinition, t, ws::CoreWorkspace) -> nothing
```

- `def::CoreDefinition` (`src/Assembly/definition.jl`) holds the state layout;
  the state setup (the topology, the material, the caloric laws of its
  components, the per-cell reference data, the caloric extension of the
  enthalpy form, and the gas transport the problem needs); the kernels, a tuple
  in assembler order; the mesh closure (`nothing` on a fixed mesh); and the
  external sources (`NoSources()`, or `VolumetricSources` for manufactured
  solutions).
- `ws::CoreWorkspace` holds the buffers of one evaluation: the derived state
  (`DerivedState`, the quantities the state mapping produces) and the kernel
  contributions (`Contributions`) with `Float64` elements, their twins for every
  dual-number element type automatic differentiation evaluates the residual at
  (built on first use, keyed by the element type), and `T_warm`, the latest
  temperatures of a plain evaluation.
- `t` is a plain number. A dual-number state evaluates in its own buffers,
  chosen by the state's element type; a forward-mode derivative in time alone
  (a dual `t` with a `Float64` state) is not supported, and every integrator
  policy takes the time derivative of the residual by finite differences
  (§15.5.1).

In the temperature form the residual is a pure function of $(\boldsymbol u,t)$
given `def`. In the enthalpy form the caloric inversion starts from `T_warm`,
so the result also depends on the evaluation history at the inversion's
round-off; `integrate` copies the set-up temperatures into `T_warm` at the start
of every integration, so repeated integrations are bit-identical.

Every physics option is a type parameter of the definition (the layout's flags,
the kernel tuple, the closure's correction), not a method of the residual:
`residual!` has one method, and an evaluation after warm-up allocates nothing
(`test/jacobian/test_allocations.jl`). The definition also holds the two pieces
of state a solve mutates, the layout's mask and the correction's layer times,
so one set-up problem serves one integration at a time; `integrate` resets both
at its start.

`setup(problem; …)` builds the definition, the buffers and everything else the
integrator needs into a `CoreProblem` (`src/Integration/setup.jl`): the problem,
`def`, `ws`, the initial state `u0` (the masses $m_{j,i}=\xi_j(z^c_i)V_i$ and the
energy state from $T(z^c_i)$ at the cell centres, the algebraic unknowns made
consistent, §15.4), the relative tolerance and the absolute-tolerance vector
(§15.5.6), the colored Jacobian (§15.6), the initial temperatures and masses,
and the deactivation state (§14.3; `nothing` without deactivation). `solve` is
`setup` followed by `integrate`; a custom loop calls the two itself.

---

## 15.3 The residual

### 15.3.1 One evaluation, in order

`residual!` runs, in this order (`src/Assembly/residual.jl`):

1. **The state mapping** (`derive!`, §13), the only reader of $\boldsymbol u$:
   the live cells from the mask; on a moving mesh the geometry from the node
   positions, with the node velocities (or the multipliers) and, under a
   variable cross section, the area and its lateral rate; the carried fields;
   the face intensities; the masses, volumes and concentrations; the
   temperatures (the stored temperature, or the inversion of the enthalpy
   inventory); their trial clamp; the heat-capacity sum; the sensible
   enthalpies; the cell properties; the pore-gas properties the gas transport
   needs.
2. **The node corrections** $c_n$ of the mesh closure (`correction!`): zero
   without a correction, the node velocities under a prescribed mesh velocity,
   the layers' relaxed equidistribution under the equidistribution correction
   (§11).
3. **The kernels**, in assembler order (`_kernels`, `src/Integration/setup.jl`):
   the reactions (if the material has any); the gas transport (if gas can
   move); the in-depth radiation (`BEER_LAMBERT`); conduction; the thermal
   boundary; and, under a correction, the inventory fluxes relative to the mesh.
   Each adds its face fluxes and cell sources to the contributions.
4. **The mesh closure** (`closure_rows!`): the cells' strain rates $\theta_i$
   from the reactions' mass sources, the history sources under a variable cross
   section, then the residuals of the velocity rows (§15.4).
5. **The external sources** (`add_sources!`).
6. **The row assembly** (`assemble_rows!`), the only writer of
   $\boldsymbol{du}$.

### 15.3.2 The rows

The assembly scatters every face flux into its owner and neighbour and writes:

$$
\frac{dm_{j,i}}{dt}=\omega_{j,i}\tilde V_i+\dot M_{j,i},\qquad
\dot M_{j,i}=-\sum_{f\in\partial i}s_{if}\,A_f\,F_{j,f},
$$

with $\omega_{j,i}\tilde V_i$ [kg s⁻¹] the cell's reaction source ($\tilde V_i$
the smoothly floored cell volume, §13.10.3; an external mass source, §13.5.3,
adds to it here and in the formation term below), $\dot M_{j,i}$ the face-flux
part of its mass rate and $s_{if}=\pm1$ the orientation of face $f$ relative to
cell $i$. With $E_i$ [W] the cell's energy sources plus its face energy fluxes
(conduction, the enthalpy the mass fluxes carry, the boundary fluxes; §7, §13.3.4)
and $C_i=\sum_j m_{j,i}c_{p,j}(\tilde T_i)$ the heat-capacity sum at the trial
temperature, smoothly floored at $C_{{\rm floor},i}$, the energy row is

$$
\text{temperature form:}\quad \frac{dT_i}{dt}=\frac{E_i-\sum_j h_j(T_i)\,\dot M_{j,i}}{C_i},
\qquad
\text{enthalpy form:}\quad \frac{dH_i}{dt}=E_i+\sum_j \omega_{j,i}\tilde V_i\,h_j(T_i),
$$

with $h_j(T)$ the sensible enthalpy of component $j$ on the $T_{\rm ref}$ datum.
The two are the same balance: $dH_i/dt=C_i\,dT_i/dt+\sum_j h_j(T_i)\,dm_{j,i}/dt$.
The enthalpy form is a linear invariant every Runge–Kutta step conserves to
round-off; the temperature form conserves energy to the integrator's order. An
external temperature-rate source $S_T$ adds $S_T$ to the temperature row and
$C_iS_T$ to the enthalpy row.

On a moving mesh the position rows are $dz_n/dt=w_n$ and the velocity rows
carry the closure's residuals; the face-intensity rows carry the residuals of
their attenuation equations; the history rows take their sources and face
fluxes; the area and lateral-rate rows take the residuals of the lateral rows
at the column sums of the live cells' history pairs and their rates; the
deformation-gradient rows (more dimensions) are $dF_i/dt=(\nabla v_s)_iF_i$.
The dead rows come last (§15.1.2).

### 15.3.3 Order dependencies

The order is fixed for correctness:

- **The corrections before the kernels.** The relative-flux kernel reads the
  node corrections.
- **The gas transport before the thermal boundary.** The surface energy balance
  of a boundary face reads the face's outward mass fluxes, which feed the
  blowing models of its convective conditions and the sensible enthalpy the
  gas carries out at the face temperature (§12.2).
- **The closure after the kernels and before the sources.** The strain rates
  are the reactions' alone, read from the mass sources the reaction kernel
  wrote; a manufactured mass source does not move the mesh.
- **The assembly last.** It is the only writer of $\boldsymbol{du}$, so every
  kernel only adds to the contributions.

---

## 15.4 The algebraic rows and the consistent state

The algebraic rows carry every non-local dependence of the physics as a local
equation (§15.9):

- **The balance rows** of the node velocities (`BalanceRowClosure`, §11): the
  pinned node at rest, $0=w_1$, and for each cell, on the velocity row of its
  node farther from the pinned node,

  $$
  0=A_{i+1}(w-c)_{i+1}-A_i(w-c)_i-(\theta_i-\theta_A)V_i ,
  $$

  the solid velocity $w-c$ satisfying the cell's volume balance. The rows are
  bidiagonal from the pinned node. Under a prescribed mesh velocity every row
  is $0=w_n-w(z_n,t)$.
- **The face-intensity rows** (`BeerLambertKernel`, §8.3): at the exposed face
  $0=I_e-a\,q_{\rm inc}(t)$, and at every other face, the face of its cell
  farther from the exposed one, $0=I_{\rm down}-I_{\rm up}\exp(-\alpha_i\Delta_i)$,
  with $\alpha_i$ the cell's absorption coefficient and $\Delta_i$ its thickness
  along the beam (§8.3.3).
- **The cross-section rows** (`LateralShrinkage`, §10):
  $0=A-A_0\,\mathrm{law}(\bar\chi_c)$ and
  $0=\theta_A-\frac{d}{dt}\log\mathrm{law}(\bar\chi_c)$, with $\bar\chi_c$ the
  column's conversion smoothly clamped into a neighbourhood of $[0,1]$.
- **The multipliers** of the potential closure on meshes of more dimensions,
  one per cell (residual level, §14.5).

Each block is solvable for its own unknowns given the differential ones: the
balance rows by one sweep from the pinned node, the intensity rows by the
attenuation recurrence from the exposed face, the cross-section rows
explicitly, the multipliers by one symmetric positive-definite solve. The
system is therefore a semi-explicit DAE of index one.

`consistent_state!(u, def, ws, t)` overwrites the algebraic unknowns of
$\boldsymbol u$ with the values that satisfy their rows at $t$, given its
differential unknowns: the face intensities by the recurrence; under a variable
cross section, one residual pass for the strain and history rates, the area and
lateral rate from the column sums, and a second pass at the consistent lateral
rate; then the node velocities by one sweep (one pass is exact, since the
corrections read positions and masses only) or the multipliers by one solve. A
closure row left above $10^{-10}$ of its scale is an error. `setup` calls it at
$t_0$, and the deactivation event calls its parts after every transfer
(§14.3.5); the integrator itself runs no initialization (`initializealg =
NoInit()`). `integrate` refuses an inadmissible initial state (§15.5.4) with an
`ArgumentError`.

---

## 15.5 Time integration

### 15.5.1 The mass-matrix DAE and the integrator policies

The integrator advances

$$
\mathsf M\,\frac{d\boldsymbol u}{dt}=\boldsymbol f(\boldsymbol u,t),\qquad
\mathsf M=\mathrm{diag}(\text{1 on differential rows, 0 on algebraic rows}),
$$

handed to OrdinaryDiffEq as

```julia
f = ODEFunction{true,FullSpecialize}(ErasedRHS(rhs!);
    mass_matrix = M, jac = ErasedJacobian{typeof(J0)}(jac!), jac_prototype = J0)
```

with $\mathsf M$ the identity when every row is differential (a fixed mesh
without in-depth radiation, an ODE) and the diagonal `mass_matrix(layout)`
otherwise; `jac!` the colored Jacobian (§15.6) and `J0 = jac_prototype(cj)` a
zero matrix of its pattern (`integrate`, `src/Integration/solve.jl`). The
integrator is initialized with the per-row absolute tolerances, the relative
tolerance, `maxiters`, `dtmin`, `dtmax`, the first step `dt` (the fixed step
with `adaptive = false`), `saveat`, the callback set of §15.5.3, the domain test
of §15.5.4 as `isoutofdomain`, `initializealg = NoInit()` and `save_end = true`.

An **integrator policy** (`IntegratorPolicy`, `src/Integration/policy.jl`) is an
algorithm together with the two rules that are part of its contract: the
**restart** after an event that modifies the state, and the **observer flag**,
what a per-step observer does with the derivative-discontinuity flag that the
callback machinery raises before every discrete callback's `affect!`.

| Policy | Method | Restart after an event | Observer flag |
|---|---|---|---|
| `rodas5p_policy(ls)` (the default, `default_policy`) | `Rodas5P`, a Rosenbrock method: one linear solve per stage with the Jacobian evaluated at every step | `:heuristic` (keyword `restart`; also `:forced_dt`, `:solver_chosen`) | cleared |
| `kencarp4_policy(ls)` | `KenCarp4`, an ESDIRK method with Newton stage solves, the raw (unsmoothed) error estimate | `:solver_chosen` | kept |
| `fbdf_policy(ls)` | `FBDF`, variable-order backward differentiation | `:forced_dt` | cleared |

`ls` is the linear solver (§15.10). Every policy takes the Jacobian from the
colored one and the time derivative of the right-hand side, where the method
needs it, by finite differences (`autodiff = AutoFiniteDiff()`). The restarts:
`:solver_chosen` keeps the integrator's own proposal; `:forced_dt` proposes
`DT_AFTER_EVENT` $=10^{-6}$ s; `:heuristic` proposes the integrator's
initial-step estimate at the state the event left (two right-hand-side
evaluations). For `Rodas5P` the heuristic restart rejects fewer steps after an
event than keeping the step and keeps the mass ledger closer than either other
kind. The flag rules: `:keep` leaves the flag raised, so the integrator
re-evaluates the derivative at the accepted state; `:clear` clears it, so a
read-only observer is numerically inert. `Rodas5P` rebuilds its Jacobian every
step and is indifferent to the flag; `FBDF` with the flag raised restarts at
order one on every step, and after an event it restarts at order one, hence the
forced small step. `solve` takes a policy through `integrator_policy` and a
linear solver through `linear_solver`; without either it uses
`default_policy(StructuredKLU())`.

### 15.5.2 The W matrix and the stage solves

Each stage solves a linear system with the matrix

$$
W=J-\frac{\mathsf M}{\gamma_{\rm RK}\Delta t},
$$

OrdinaryDiffEq's sign convention (the negative of the textbook
$I/(\gamma_{\rm RK}\Delta t)-J$, with the mass matrix in place of the identity;
the sign changes neither the structure nor the solve). `Rodas5P` solves one
linear system per stage with $W$ at the step's Jacobian and no Newton iteration;
`KenCarp4` and `FBDF` solve each stage by a Newton iteration on $W$, reusing the
Jacobian across steps as the method decides. On an algebraic row
$\mathsf M_{rr}=0$, so the row of $W$ is the Jacobian's row: the constraint is
linearized exactly in every stage. On a dead row $W$ is a scaled identity row
(§14.3.1). Every row of the Jacobian's pattern stores its diagonal, so $W$'s
stored entries are a subset of the pattern's whatever $\Delta t$.

### 15.5.3 Output times, callbacks and observers

States are saved at `saveat` (default `problem.output_times`); the state a solve
ends on is always saved. A solve its controller stopped saves its last
proposal, which it rejected, at the last accepted time; `integrate` reports the
accepted state there instead. Every callback reaches the problem through a
fixed-signature function wrapper (§15.5.7) and is assembled into one
`CallbackSet` (`integrate`):

| Callback | Kind | Role |
|---|---|---|
| step observer | discrete, every accepted step | records the largest constraint residual of the accepted states; ends the solve with `ReturnCode.Failure` and a warning once the domain test gave up (§15.5.4) or a fixed step left the domain; applies the Newton-failure step cap (§15.5.5); applies the flag rule |
| located deactivation | continuous, on the threshold margin (`deactivation`, `located_events`) | §14.3.5 |
| deactivation safety net | discrete, every accepted step (`deactivation`) | §14.3.5 |
| layer times | continuous, on the activation margin (every moving mesh) | sets the layers' equidistribution times (§14.2.2) |
| thickness stop | continuous, on the live thickness minus `min_thickness` (moving mesh) | ends the solve with outcome `:min_thickness` |
| conservation recorders | continuous, never firing (one per `ConservationRecorder` in `callback`) | read each step's end before its events (§16) |
| progress bar | discrete (`progress = true`) | a read-only observer |
| the caller's `callback` | as given | a read-only observer: the flag rule follows each discrete callback's `affect!` |
| accepted-state recorder | discrete, every accepted step, last | records the accepted state and commits the recorders' step integrals; applies the flag rule |

Continuous callbacks are located inside the step on the integrator's
interpolant, and the earliest crossing ends the step there; the discrete ones
run at each accepted step's end in the order listed. The step observer is
always installed, so a further read-only observer leaves the trajectory
unchanged under every policy: a probe callback, alone or in a `CallbackSet`,
leaves the saved times, the states, the accepted steps and the right-hand-side
evaluations bit-identical (`test/integration/test_policies.jl`). A continuous
callback of the caller's raises the flag at its event as any event does.

### 15.5.4 The domain test

`DomainCheck(cp)` (`src/Integration/checks.jl`) is the integrator's
`isoutofdomain` test. A proposed state is rejected, and the step cut by the
controller's smallest factor, when

- it is **inadmissible** (`inadmissible`): a live component mass more negative
  than its absolute tolerance; a live cell at or below absolute zero (its
  energy state at or below the value of zero kelvin at its masses,
  `below_absolute_zero`); on a moving mesh a live cell whose nodes have crossed
  (a non-positive volume); or a non-positive cross-section area. The test reads
  the raw state, never through a guard: the trial clamp keeps such a state
  evaluable, and this test keeps it from being accepted;
- under Darcy–Fick transport a live cell's porosity
  $1-\sum_j m_j/(\rho_{i,j}(T)V)$, read from the raw masses and volume, leaves
  $[0,1]$ by more than the skeletal volume of its mass tolerances;
- an algebraic row's scaled residual exceeds one (`constraint_residual`):

  $$
  \max_{r\ \text{algebraic, live}}\frac{|f_r|}{\mathrm{rtol}\,s_r+\mathrm{atol}_r}>1,
  $$

  with $s_r$ the row's own equation scale and $\mathrm{atol}_r$ the absolute tolerance of
  its unknown in the residual's units:

  | Row | Scale $s_r$ | $\mathrm{atol}_r$ |
  |---|---|---|
  | balance row of cell $k$ | $(\lvert\theta_k\rvert+\lvert\theta_A\rvert)V_k+\sum_f A_f\lvert w_f\rvert+\sum_f A_f\lvert c_f\rvert$ | the node's velocity tolerance times its face area |
  | pinned or prescribed velocity | $\lvert w_n\rvert$ | the node's velocity tolerance |
  | face intensity | the incident-intensity scale | its tolerance |
  | area | $\lvert A\rvert$ | its tolerance |
  | lateral rate | $\lvert\theta_A\rvert$ plus the volume-weighted mean $\lvert\theta\rvert$ | its tolerance |

The algebraic rows are also in the integrator's error norm through their
absolute tolerances; the constraint test is a separate check of the
constraint's accuracy in its own equation's units, which the state tolerance of
a velocity, in length per time, does not provide.

The test reads only the proposals the integrator's error test accepts: it is
bound to the integrator's step controller (`bind_error_test!`) and returns
"inside the domain" for a proposal the error test rejects, which that test and
its controller's own cut then handle. The residual is evaluated once per tested
state. After more than `MAX_CONSTRAINT_REJECTIONS` $=50$ consecutive constraint
rejections the state is accepted and the monitor marked failed, and the step
observer ends the solve with `ReturnCode.Failure` and a warning. A proposal
rejected for a cell at or below absolute zero is recorded, and a solve that
stops there warns naming the cell and the time.

A fixed-step solve (`adaptive = false`) never reaches the domain test; the step
observer holds each accepted state to admissibility and ends the solve with
`ReturnCode.Failure` at the first inadmissible one. The counts are in
`sol.statistics`: `domain_rejections` and `largest_constraint_residual` (one or
below passes; zero for a fixed-step solve), beside the integrator's accepted
and rejected steps, right-hand-side evaluations, Jacobians, linear solves and
nonlinear failures.

### 15.5.5 The Newton-failure step cap

After an accepted step whose nonlinear solves needed retries, the proposed step
is held at or below $1.25\,\Delta t_{\rm acc}$ (`FAILURE_MEMORY_GROWTH`), with
$\Delta t_{\rm acc}$ the step just accepted, for
the next $20$ accepted steps (`FAILURE_MEMORY_STEPS`); a further failure
restarts the memory from that step (`cap_step!`). The error controller alone
regrows the step to its accuracy limit after a failure, and on a strongly
nonlinear front that limit can lie beyond the range in which the stage Newton
converges, which sets up a fail–halve–regrow rhythm. A run without nonlinear
failures never engages the cap and keeps its trajectory exactly. Rosenbrock
methods have no nonlinear solve, so under the default policy the cap never
engages.

### 15.5.6 Per-row absolute tolerances

`abstol_vector` (`src/States/tolerances.jl`) gives every row its own absolute
tolerance, from the cell's initial condensed mass $m^{\rm s}_i(0)$ (for a cell
that starts with none, its total initial mass, or 1 kg m⁻³ times its volume),
the deactivation threshold $\varepsilon_m$ (`threshold`), `reltol` and `abstol`:

| Row | Absolute tolerance |
|---|---|
| mass of any component in cell $i$ | $\mathrm{rtol}\,\varepsilon_m\,m^{\rm s}_i(0)$, never a per-species inventory, which is zero for a product |
| temperature | `abstol` [K] |
| enthalpy inventory | $\varepsilon_m\,m^{\rm s}_i(0)\,c_{p,\rm ref}\,$`abstol`, $c_{p,\rm ref}$ the cell's mass-weighted heat capacity at its initial temperature, so a cell at the threshold resolves its temperature to `abstol` as the temperature form does |
| face intensity | $\mathrm{rtol}\,I_{\rm scale}$, $I_{\rm scale}$ the largest incident intensity over the solve window (at least 1 W m⁻²) |
| node position | $\mathrm{rtol}\,\varepsilon_m\,\Delta z_n^{\min}$, $\Delta z_n^{\min}$ the smallest initial size of the cells around node $n$ |
| node velocity | $\mathrm{rtol}\,\varepsilon_m\,\Delta z_n^{\min}/\Delta t_0$, $\Delta t_0$ the first step `dt` |
| multiplier (more dimensions) | $\mathrm{rtol}\,\varepsilon_m\,\Delta z_i(0)^2/\Delta t_0$, $\Delta z_i(0)=V_i(0)^{1/D}$ the cell's initial size |
| history pair | the cell's mass tolerance |
| deformation gradient | $\mathrm{rtol}$ |
| area | $\mathrm{rtol}\,A_0$ |
| lateral rate | the smallest velocity tolerance over the initial thickness |

The defaults are `reltol` $=10^{-6}$ and `abstol` $=10^{-8}$ K. A position error
within its tolerance moves a cell's concentrations by the relative amount its
mass tolerance allows. The same mass tolerance is the scale below which the
state mapping's positive part and the caloric extension act (§13).

### 15.5.7 Problem-independent integrator types

OrdinaryDiffEq compiles its stepper, nonlinear solver and factorization for the
types of the functions it is handed. A closure over a problem's definition and
buffers carries the material, boundary-condition and kernel types, so every
problem would compile that whole stack again. The right-hand side and the
Jacobian therefore reach the integrator behind fixed-signature
`FunctionWrapper`s (`src/Integration/erased.jl`):

- `ErasedRHS(rhs!)`: calls with `Vector{Float64}` state and `Float64` time go
  through the wrapper; any other argument types go to `rhs!` by dynamic
  dispatch;
- `ErasedJacobian{M}(jac!)`: for the matrix type `M` the integrator allocates
  from the prototype (`SparseMatrixCSC{Float64,Int}`); a call with any other
  matrix type is a `MethodError`, since a converted copy would hand the
  integrator a matrix it never sees.

The domain test and every callback reach the problem through wrappers in the
same way; only the actions on the integrator itself (`terminate!`,
`set_proposed_dt!`, `auto_dt_reset!`, `derivative_discontinuity!`, the reads of
its counters) stay in the callbacks. `solve` converts the tolerances, step
bounds and `saveat` to `Float64`, since they enter the integrator's type. The
integrator's type then depends on the solver configuration only: whether the
mesh moves, whether in-depth radiation adds algebraic rows on a fixed mesh (a
mass matrix), whether cells are deactivated, the policy and its linear solver,
and the optional callbacks installed. Two decks that differ in component count,
property types, reactions, boundary conditions, radiation model and gas
transport share one integrator type and one solution type in each such
configuration (`test/integration/test_policies.jl`). A problem with other
material or boundary-condition types compiles only the Pyrolysis kernels, the
state mapping, the colored Jacobian and the solution's reconstruction. The
precompile workloads (`src/precompile.jl`) cache the integrator for the fixed
and moving meshes under surface absorption and Beer–Lambert absorption, Fickian
and Darcy–Fick transport, with and without deactivation, and the enthalpy form
once.

---

## 15.6 The colored Jacobian

The Jacobian $J=\partial\boldsymbol f/\partial\boldsymbol u$ is the forward-mode
derivative of the whole residual, evaluated on a sparsity pattern built from
the mesh connectivity (§15.11) with its columns colored (§15.7). There is one
Jacobian and no option that selects another: every solve uses it, and the dense
forward-mode Jacobian of §15.13 is the reference it is tested against.

`ColoredJacobianCache(def, ws, u, t)` (`src/ColoredJacobian/colored.jl`), built
once by `setup` at the initial state, holds

- the residual as `f!(du, u, t)` (`ResidualFunction(def, ws)`), with the time a
  constant of the differentiation;
- the backend, DifferentiationInterface's `AutoSparse(AutoForwardDiff(;
  chunksize))` with the connectivity pattern as a known sparsity
  (`KnownJacobianSparsityDetector`) and SparseMatrixColorings' greedy column
  coloring;
- the preparation (`prepare_jacobian`), done once per problem structure and
  reused at every state and time;
- the full pattern, the colored pattern (the dense rows reduced to their
  diagonal, §15.9), the decompression target and its scatter into the full
  pattern, the dense-row assembler (`LateralRows`, or `nothing`), a residual
  buffer and the color count.

`jacobian!(J, cj, u, t)` writes the Jacobian at $(\boldsymbol u,t)$ into `J`,
whose structure must be the pattern's: the decompression writes stored values
in order and checks only the matrix size and stored-entry count, so the cache's
construction asserts that the prepared pattern and the colored pattern share
`colptr` and `rowval`, and the integrator's `jac!` checks
`same_structure(J, pattern)` on every call and errors if the matrix lost the
structure. The integrator allocates its Jacobian from `jac_prototype(cj)`, a
zero matrix with exactly the pattern's structure. After warm-up a Jacobian
evaluation allocates nothing (`test/jacobian/test_allocations.jl`, 100 cells,
every physics configuration, both energy forms).

---

## 15.7 Coloring, batching and decompression

Two columns may share a color when no row holds both; one dual-number residual
evaluation then carries the derivatives of every column of a color at once,
each row reading its entry from the column of that color it stores.

**Coloring.** SparseMatrixColorings colors the columns greedily in each of four
vertex orders (`COLORING_ORDERS`: natural, largest-first, smallest-last and
incidence-degree) and keeps the coloring with the fewest colors $n_c$. The
largest number of stored entries in a row is a lower bound on any column
coloring; on every configuration of the test suite the chosen coloring reaches
it, and the count is constant in the cell count
(`test/jacobian/test_pattern.jl`). On a line it follows

$$
n_c=\begin{cases}
3N_C+9+2\delta_{\rm BL}+7\delta_\chi-n_{\rm mesh}, & \text{no relative fluxes},\\
5N_C+13+2\delta_{\rm BL}+11\delta_\chi-n_{\rm mesh}, & \text{relative fluxes (a correction)},
\end{cases}
\qquad
n_{\rm mesh}=\begin{cases}6 & \text{fixed mesh},\\ 2 & \text{moving mesh},\\ 0 & \text{prescribed mesh velocity},\end{cases}
$$

with $\delta_{\rm BL}$ and $\delta_\chi$ the flags of §15.1.1. For two
components on 12 cells this gives 9 colors on a fixed mesh, 11 under in-depth
radiation, 13 on a moving mesh, 15 on a moving mesh under in-depth radiation
and 21 on a moving mesh with deactivation, which installs the equidistribution
correction and its relative fluxes (measured with the script of §15.13).

**Batching.** One dual-number evaluation carries at most `MAX_CHUNK` $=24$
partials. The colors are split into the fewest batches of at most 24, of equal
width $\lceil n_c/\lceil n_c/24\rceil\rceil$, so up to 24 colors cost one
evaluation of the residual at a dual type of that width; the dense reference
(§15.13) uses the same width, so both evaluate the residual at one dual-number
type at every mesh size.

**Decompression.** DifferentiationInterface decompresses the colored result
into a `SparseMatrixCSC` with the colored pattern. Without dense rows that
matrix is the Jacobian. With them, the colored result is scattered into the
full pattern through a precomputed index map and the two dense rows are then
written by `LateralRows` (§15.9).

---

## 15.8 What the Jacobian differentiates

The colored Jacobian is the derivative of the residual as executed: every
kernel, every guard and every iteration inside one evaluation is differentiated
by the dual numbers that pass through it. No derivative is written by hand, and
a term added to a kernel enters the Jacobian with it. Three kinds of code inside
the residual need care to keep that derivative correct and finite.

- **Smooth guards.** Every trial-state guard is $C^1$ and exactly the identity
  on the states the tolerances can distinguish (`src/States/guards.jl`): the
  trial temperature, the cell temperature smoothly clamped into the property
  bracket $[200,3000]$ K with 1 K transitions, which keeps $\exp(-E/RT)$ finite
  with a finite derivative at zero and negative trial temperatures; the
  one-sided positive part of the masses, exactly zero at and below zero and the
  mass itself above the cell's mass tolerance; the smooth floors of the volume,
  the heat-capacity sum, the conductivity and the porosity; the activation
  switch of the correction. The derivatives are therefore finite at every real
  state, and they are the physical ones on every physical state.
- **The surface energy balance.** Each boundary face solves
  $q_{\rm BC}(T_s)=k(T_s-T_{\rm cell})/d$ by a Newton iteration
  (`solve_boundary_temperature!`, §12.2) whose iterates carry the dual numbers.
  In dual mode the first iterate is never accepted, so at least one Newton
  update propagates the partials, and the converged dual fixed point carries
  the implicit-function derivative of $T_s$ with respect to the cell state and
  the face's outward mass fluxes; the blowing models are evaluated at every
  iterate from those fluxes, so their derivatives enter exactly.
- **The caloric inversion** (enthalpy form). The temperature of a cell is the
  root of its extended caloric relation (§13). For dual inputs the primal root
  is found by the bracketed Newton iteration on the values, then corrected by
  one Newton step per level of dual nesting, which carries the implicit
  derivatives exactly (`invert_caloric`).

The lateral-rate row differentiates the log-area along the history rates by a
nested `ForwardDiff.derivative` inside the residual; its Jacobian row is
assembled separately (§15.9). Nothing else in the residual iterates: the
non-local recurrences of the physics are algebraic unknowns of the state
(§15.4).

The dense reference equals central differences of the residual over every
column, per entry against each row's largest entry to $10^{-6}$, at the
containment states of §15.13 (`test/jacobian/test_containment.jl`), so the
derivative of the code is also the derivative of the function it computes.

---

## 15.9 Non-local dependence: algebraic rows and the two dense rows

Figure 15 in §15.1 compares the actual sparsity patterns with and without correction and variable cross section. Its purple outlines identify the two global algebraic rows; pale entries separate structural dependencies from the derivatives that are nonzero at the initial state.

A moving mesh's node velocity is the swept integral of the strain rates of every
cell between it and the pinned node, and the intensity reaching a cell is the
incident flux attenuated by every cell above it. Eliminated, either would couple
every node or cell to every other one, filling the Jacobian. Both are instead
algebraic unknowns of the state whose rows are local (§15.4): a balance row
reads one cell and its two nodes, besides what the corrections of those nodes
read, and an intensity row reads one cell and its two faces. The pattern stays
banded on a line, and the coloring stays bounded independently of the cell
count.

What remains dense are the two rows of a variable cross section: the area row
and the lateral-rate row read every cell's history pair, the lateral one also
their rates. One dense row would make every pair of columns conflict, so the
coloring sees them reduced to their diagonal (`dense_rows`, `_colored_pattern`)
and `LateralRows` (`src/ColoredJacobian/colored.jl`) writes them after the
colored pass. Every dependence of the two rows on the state passes through six
column sums

$$
\boldsymbol s=\Bigl(\textstyle\sum_i G_i,\ \sum_i m_{{\rm dry},i},\ \sum_i\dot G_i,\ \sum_i\dot m_{{\rm dry},i},\ \theta_A,\ A\Bigr)
$$

over the live cells, with $(G_i,m_{{\rm dry},i})$ cell $i$'s history pair and
$(\dot G_i,\dot m_{{\rm dry},i})$ its rates (`lateral_rows`), so each entry is

$$
\frac{\partial f_\ell}{\partial u_c}=\sum_{q=1}^{6}\frac{\partial f_\ell}{\partial s_q}\,\frac{\partial s_q}{\partial u_c},
\qquad \ell\in\{A,\theta_A\},
$$

with $f_\ell$ the residual of row $\ell$, $\partial f_\ell/\partial s_q$ by forward
differentiation of the six-number function, $\partial(\sum\dot G)/\partial u_c$ and
$\partial(\sum\dot m_{\rm dry})/\partial u_c$ the column sums of the live history rows of
the colored Jacobian, which the colored pass gets exactly, and the remaining
$\partial s_q/\partial u_c$ unit vectors. The cell rows read the area and the
lateral rate as ordinary columns of the pattern.

---

## 15.10 The linear solve of W: StructuredKLU

The default linear solver of every layout is `StructuredKLU()`
(`default_linear_solver()`, `src/Integration/linear_solve.jl`), KLU applied to
$W$ on the entries of $W$ that have been nonzero at some step of the solve.

The Jacobian's prototype holds every entry the connectivity allows, most of
which are zero at any one step, and KLU factors a stored zero like any other
entry, so factoring $W$ on the prototype's structure carries the fill of entries
that never matter. For a mass matrix that is not a multiple of the identity,
OrdinaryDiffEq forms $W$ by a sparse broadcast that stores only the entries
that compute to nonzero, so $W$'s own structure changes whenever an entry passes
through zero, and a KLU on that structure would analyse afresh at every such
step.

The cache keeps the union of the prototype's entries that have been nonzero in
this solve, with the diagonal, and copies each fresh $W$ into a matrix of that
structure, where the inner KLU refactors it numerically. A nonzero outside the
union grows it, and the inner KLU analyses the grown structure afresh; the union
only grows, so the analyses stop once every entry the solve reaches has
appeared. $W$'s stored entries must be a subset of the prototype's (the
integrator hands the solver's initialization an all-ones $W$ of the prototype);
an entry outside the prototype is an error. The tests check that the union is
smaller than the prototype, that the inner KLU analyses once per growth of the
union, and that repeated solves are bit-identical under `StructuredKLU()`,
LinearSolve's `KLUFactorization()` and `UMFPACKFactorization()`, on fixed and
moving meshes and under in-depth radiation (`test/integration/test_policies.jl`).

Another LinearSolve factorization is selected with `linear_solver` (for example
`UMFPACKFactorization()`), which `solve` passes to the default policy, or by
passing a policy built with it (`integrator_policy = fbdf_policy(KLUFactorization())`).
The linear solve is a sparse direct factorization; the package configures no
iterative solver or preconditioner.

---

## 15.11 The sparsity pattern from connectivity

`jacobian_pattern(def)` (`src/ColoredJacobian/pattern.jl`) builds the
structural pattern of $J$ from the connectivity alone: the cells within a reach
of each cell, their nodes, and the nodes' cells. It is not traced from the code;
§15.13 tests that it contains the dense reference. The reach of a cell row is two
faces where a kernel reads a neighbour's neighbour (the inventory fluxes relative
to the mesh, whose face states are reconstructed from two cells upwind; the
cross-diffusion part of tensor conduction in more dimensions, which reads the
neighbours' least-squares gradients), one otherwise.

- **A cell row** (a mass, the energy state, the history pair) reaches the
  masses, energy states and history pairs of the cells within the reach and, on
  a moving mesh, the positions of their nodes: conduction reads the neighbouring
  cell centres, the concentrations the neighbours' volumes, the relative fluxes
  everything the corrections of the cell's nodes read; under a prescribed mesh
  velocity with relative fluxes, also the velocities of its nodes. On a line it
  reads the intensities of its two faces and, under a variable cross section,
  the area.
- **A face-intensity row** reaches the intensities of nodes $f-1$ to $f+1$, the
  unknowns of the cells on either side of its face (the attenuation of the cell
  it leaves, the absorptivity of the exposed cell), on a moving mesh the
  positions of those nodes, and the area.
- **A node-position row** reaches its node's velocity, or, when the closure's
  unknowns are the cells' multipliers, the multipliers of the cells around the
  node and the positions of their nodes.
- **A velocity row** reaches its node's velocity and position; on a line a
  balance row also reaches the velocity and position of its cell's other node,
  its cell's masses and energy state (through the strain rate), the area and the
  lateral rate, and under an equidistribution correction what the corrections
  of its cell's two nodes read: the positions of the nodes of the cells beside
  it and, when the monitor reads the masses, those cells' masses.
- **A multiplier row** reaches the multipliers of the cells that share a node
  with its cell, the positions of their nodes, and its cell's masses and energy
  state; a **deformation-gradient row** also reaches the velocities of its
  cell's nodes, or the multipliers and node positions those velocities read.
- **The cross-section rows**: the area row reaches every cell's history pair and
  the area, the lateral-rate row everything the history rows reach, the history
  pair, the area and the lateral rate.

Every row holds its diagonal. The pattern is fixed for a solve: dead rows keep
their entries and take zero values (§15.12). Radial meshes (cylindrical and
spherical coordinates) share the planar mesh's pattern, colors and bandwidth.
For two components on 12 cells the half-bandwidth is 5 on a fixed mesh, 6 under
in-depth radiation, 9 on a moving mesh, 11 on a moving mesh under in-depth
radiation and 14 with the equidistribution correction (measured with the script
of §15.13); it is constant in the cell count.

---

## 15.12 The Jacobian through deactivation

Deactivation (§14.3) and element death (§14.5) mask rows of the layout and leave
the state length and the pattern unchanged. The Jacobian of a dead differential
row is a zero row and that of a dead algebraic row the unit row $e_r^{\mathsf
T}$; the entries of live rows in dead columns are zero, since the kernels skip
every face beyond the live cells and the boundary conditions act on the live
faces. The preparation, the coloring, the integrator's matrices and the linear
solver's union structure all carry through every event, and $W$ splits into the
live block and scaled identity rows.

The tests apply `deactivate!` to the state and compare the colored Jacobian
with the dense reference afterwards: a single transfer on fixed and moving
meshes, a run merging into the cell beneath, a run through a layer's innermost
cell into the substrate, a merge under a variable cross section, and an exposed
cell deactivated under in-depth radiation, where the exposed face becomes the
receiver's far face. In every case the pattern contains the dense Jacobian and
the two agree per entry to $10^{-12}$ (`test/deactivation/test_deactivation.jl`,
`test/jacobian/test_containment.jl`). Element death on a 3 × 3 mesh leaves the
pattern unchanged and containing the Jacobian (`test/multid/test_element_death.jl`).

---

## 15.13 The dense reference and the containment tests

`dense_jacobian(def, ws, u, t)` (`src/ColoredJacobian/oracle.jl`) differentiates
every column of the residual by forward mode, without the pattern or the
coloring, in batches of the colored Jacobian's width. It is the reference: it
differentiates the same `residual!` the integrator integrates.

- `missing_entries(J_ref, pattern)` lists the nonzero entries of the reference
  that the pattern does not contain, empty when the pattern contains it.
- `compare_entries(J, J_ref; rtol = 1e-12, atol = 0.0)` compares each entry (row $r$,
  column $c$) with a row-magnitude floor,

  $$
  |J_{rc}-J^{\rm ref}_{rc}|\le \mathrm{atol}+\mathrm{rtol}\,\max\bigl(|J^{\rm ref}_{rc}|,\ \max_{c'}|J^{\rm ref}_{rc'}|\bigr),
  $$

  and returns `(passed, max_abs_dev, worst_ratio, worst_index)`. The rows of a
  pyrolysis Jacobian span many orders of magnitude (an Arrhenius mass row above
  a velocity row), so a single global normalization would let a wholly missing
  small-magnitude row pass behind the matrix maximum; with the row floor a
  missing entry cannot hide in its own row, while round-off at a true zero of a
  stiff row stays within $\mathrm{rtol}$ times the row's maximum.

The containment tests (`test/jacobian/test_containment.jl`) require, at every
state they visit, a finite residual and reference, a pattern that contains the
reference, a colored Jacobian with the pattern's structure, and per-entry
agreement to $10^{-12}$. They cover fixed and moving columns of a PMMA-like and a
permeable wood-like material without radiation, with surface absorption and with
in-depth absorption, in both energy forms, at uniform hot, hot and gas-heavy,
partially decomposed and mixed reacting states and at trial states outside the
physical domain; the chain-rule regimes of the smooth reaction gates and of
permeability mixing; the cross section's back-reaction; radial meshes and a
prescribed mesh velocity; the states after a deactivation; the boundary mass
flux on the energy row; and a surface absorptivity and an emissivity that depend
on temperature. The dense reference itself equals central differences
(§15.8). The pattern tests (`test/jacobian/test_pattern.jl`) check the color
counts and the per-entry comparison's behaviour on a stiff row above a small
one.

The same checks run from a script through `Pyrolysis.Internal`:

```julia
using Pyrolysis
using Pyrolysis.Internal: jac_prototype, jacobian!, dense_jacobian, missing_entries,
    compare_entries, state_length, is_differential
using SparseArrays

solid = SolidComponent(:solid, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, α = 2.0, λ = 1e-6)
gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5)
material = Material(
    name = :slab,
    components = (solid, gas),
    reactions = (Reaction(:pyrolysis, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(1e-3, 12),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 100e3) + RadiativeBC(T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 600.0 + 2e5 * z,
    tspan = (0.0, 10.0),
)

for (name, kw) in (
    ("fixed", (;)),
    ("fixed, Beer-Lambert", (radiation_model = BEER_LAMBERT,)),
    ("moving", (use_ale = true,)),
    ("moving, Beer-Lambert", (use_ale = true, radiation_model = BEER_LAMBERT)),
    ("moving, deactivation", (use_ale = true, deactivation = true)),
)
    cp = setup(problem; kw...)
    J = jac_prototype(cp.jacobian)
    jacobian!(J, cp.jacobian, cp.u0, 0.0)
    Jd = dense_jacobian(cp.def, cp.ws, cp.u0, 0.0)
    rows = rowvals(J)
    bandwidth = maximum(abs(rows[k] - c) for c in axes(J, 2) for k in nzrange(J, c))
    println(rpad(name, 22), state_length(cp.def.layout), " rows (",
            count(!, is_differential(cp.def.layout)), " algebraic), half-bandwidth ",
            bandwidth, ", ", cp.jacobian.n_colors, " colors, ",
            length(missing_entries(Jd, cp.jacobian.pattern)), " missing, passed = ",
            compare_entries(Matrix(J), Jd; rtol = 1e-12).passed)
end
```

It prints

```
fixed                 36 rows (0 algebraic), half-bandwidth 5, 9 colors, 0 missing, passed = true
fixed, Beer-Lambert   49 rows (13 algebraic), half-bandwidth 6, 11 colors, 0 missing, passed = true
moving                62 rows (13 algebraic), half-bandwidth 9, 13 colors, 0 missing, passed = true
moving, Beer-Lambert  75 rows (26 algebraic), half-bandwidth 11, 15 colors, 0 missing, passed = true
moving, deactivation  62 rows (13 algebraic), half-bandwidth 14, 21 colors, 0 missing, passed = true
```

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

- **Time scheme.** Gpyro is first-order backward Euler with an outer
  source-linearization sweep. Pyrolysis.jl integrates the mass-matrix DAE with
  an adaptive higher-order method, `Rodas5P` by default, with error control on
  every row and the restart after an event part of the integrator's contract.
- **Jacobian.** Gpyro never forms a global Jacobian; its Patankar linearization
  produces a tridiagonal coefficient matrix per equation, and the coupling
  between equations is handled by the outer iteration. Pyrolysis.jl forms the
  full coupled Jacobian of every row with respect to every unknown, by colored
  automatic differentiation, so cross-equation coupling is inside each stage
  solve. The banded block structure of our pattern on a line is the multi-field
  analogue of Gpyro's per-equation tridiagonal matrix.
- **Mesh.** Gpyro uses fixed (or prescribed variable-$\Delta z$) control volumes;
  Pyrolysis.jl moves the mesh with the solid and carries the node velocities as
  algebraic unknowns, so the volume-change kinematics stay local rows.

### 15.14.2 ThermaKin2Ds

ThermaKin2Ds (Stoliarov et al.) solves the 1D (and 2D) condensed phase with a
**Crank–Nicolson** scheme in the primary (through-thickness, here $x$/our $z$)
direction and explicit treatment in the transverse directions, on a
**dynamically adjusted mesh**. Mapping: ThermaKin's stoichiometric coefficients
$\theta_i^j$ are our $\nu_{i,j}$ (we reserve $\theta$ for strain rate, overload
resolution A2); its component concentrations are our $\xi_j$.

Contrasts:

- **Time scheme.** Crank–Nicolson is second-order and A-stable but not L-stable,
  so it can ring on the stiff Arrhenius transients at ignition; `Rodas5P`,
  `KenCarp4` and `FBDF` are L-stable or stiffly accurate methods with embedded
  error control.
- **Energy convention.** ThermaKin (Eq. 17) carries the gas enthalpy by the
  convective term $c_g(\mathbf J_g\cdot\partial T/\partial x)$ with every
  component in the storage term. Pyrolysis.jl stores every component as well and
  carries the gas enthalpy in flux form, as the sensible enthalpy the face mass
  fluxes carry (§7). The heat-of-reaction sign differs in publication ($h>0$
  exothermic in ThermaKin, endothermic in Gpyro and internally here) and is mapped
  by overload resolution H1.
- **Mesh.** ThermaKin's dynamic mesh adjustment is the closest analogue to the
  moving mesh here; both track a receding surface. In Pyrolysis.jl the mesh
  motion is part of the coupled system: the node positions and velocities are
  unknowns, and the Jacobian carries their coupling to the physics exactly.

### 15.14.3 FDS

The FDS (McDermott et al.) solid-phase pyrolysis model integrates a 1D Fourier
conduction equation with Arrhenius source terms; the surface boundary couples to
the gas-phase CFD. FDS uses a node-based 1D solid grid and an implicit/explicit
split appropriate to its multiphysics coupling. Pyrolysis.jl is a standalone
condensed-phase solver: it does not solve a gas-phase flow field, but it carries
multi-component effective properties, in-depth radiation, Darcy gas transport
and a moving mesh with a variable cross section. FDS carries no in-depth gas, so
its energy equation is the condensed-phase limit of ours; it shares the
Arrhenius kinetic form. On the temporal and Jacobian axis the comparison is
methodological (implicit conduction with explicit or implicit coupling) rather
than a mapping of Jacobian internals.

---

## 15.15 Worked summary of the discrete system

The integrator advances, for each live cell $i$, component $j$ and node $n$,

$$
\frac{dm_{j,i}}{dt}=\omega_{j,i}\tilde V_i+\dot M_{j,i},\qquad
\frac{dT_i}{dt}=\frac{E_i-\sum_j h_j(T_i)\,\dot M_{j,i}}{C_i}
\quad\Bigl(\text{or }\ \frac{dH_i}{dt}=E_i+\sum_j \omega_{j,i}\tilde V_i\,h_j(T_i)\Bigr),
$$

$$
\frac{dz_n}{dt}=w_n,\qquad
0=A_{i+1}(w-c)_{i+1}-A_i(w-c)_i-(\theta_i-\theta_A)V_i,\qquad
0=I_{\rm down}-I_{\rm up}e^{-\alpha_i\Delta_i},
$$

with the cross-section rows under a variable cross section and identity rows on
the dead rows, as $\mathsf M\,d\boldsymbol u/dt=\boldsymbol f(\boldsymbol u,t)$.
Its Jacobian $J=\partial\boldsymbol f/\partial\boldsymbol u$ is the colored
forward-mode derivative of the whole residual on the connectivity pattern, with
the two dense rows of a variable cross section assembled from six column sums;
each stage solves a linear system in
$W=J-\mathsf M/(\gamma_{\rm RK}\Delta t)$, factored by KLU on the union of its
nonzeros; accepted states pass the domain test; and the whole construction is
checked against the dense forward-mode reference.

---

## 15.16 Limitations and open issues

- **Full solves on meshes of more dimensions.** The layout, the pattern, the
  coloring and the residual handle `Mesh{2}` and `Mesh{3}` meshes with the
  potential closure, but gas transport, in-depth radiation, the variable cross
  section, the equidistribution correction and deactivation are a line's, and
  full solves on such meshes are not a supported path.
- **The time derivative.** No policy differentiates the residual in time
  alone by dual numbers; the time derivative is a finite difference.
- **Jacobian cost.** One Jacobian costs one dual-number residual per batch of up
  to 24 colors. The color count grows with the number of components and with
  the relative fluxes of a correction ($5N_C$ against $3N_C$), and on meshes of
  more dimensions the second-neighbour terms make it larger; it is measured, not
  assumed.
- **The union structure only grows.** `StructuredKLU` analyses afresh at each
  growth of its union, so the first steps of a solve, and any step at which a
  fresh entry becomes nonzero, carry an analysis.
- **The enthalpy form's history dependence.** Its residual depends on the
  evaluation history at the caloric inversion's round-off (§15.2).
- **Fixed-step solves** bypass the domain test; only the admissibility of each
  accepted state is checked.
- **No derivative through the solve.** Sensitivities are finite differences of
  whole solves (§17); the Jacobian serves the integrator only.
