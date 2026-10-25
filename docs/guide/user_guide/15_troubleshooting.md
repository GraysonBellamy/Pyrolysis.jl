# 15. Troubleshooting and Performance

This chapter is a field guide for when a Pyrolysis.jl simulation fails, runs
slowly, or produces results you do not trust. It is organised by symptom: each
section names a class of failure, explains where in the solver it arises, and
gives the `solve` keywords, the solution fields or the model changes that
address it. Most failures can be localised before any input is changed, from
four things every solve leaves behind: the return code, the outcome, the solver
statistics and the warnings `solve` prints. For the keywords themselves see
Chapter 8 and the reference tables of Chapter 16; for the numerics behind them,
the Technical Reference chapters named in each section.

Coordinate and sign conventions follow the unified nomenclature: `z = 0` is the
back face (tag `:bottom`), `z = L` is the exposed surface (tag `:top`);
positive flux is transport in `+z`; boundary heat and mass fluxes are positive
into the domain; heat of reaction `h > 0` is endothermic with
`Q_rxn = −Σ h_r r_r`.

---

## 15.1 First steps: the return code, the outcome and the statistics

Once integration has started, `solve` always returns a `PyrolysisSolution`; a
solve that fails returns the solution up to the last state the integrator
accepted. Read how it ended before changing anything.

### 15.1.1 The return code and the outcome

Two fields say how a solve ended: `retcode`, the integrator's return code as a
`Symbol`, and `outcome`, which the deactivation and the thickness stop set when
they end the solve.

```julia
sol = solve(problem; use_ale = true, deactivation = true)

sol.retcode        # :Success, :Terminated, :Failure, :DtLessThanMin, :MaxIters, …
sol.outcome        # :running, :burned_out, :min_thickness or :interior
succeeded(sol)     # :Success, or :Terminated by a complete burnout or the thickness stop
completed(sol)     # :Success only: the solve reached tspan[2]
```

| `retcode` | Meaning | Typical cause | See |
|---|---|---|---|
| `:Success` | reached `tspan[2]` | — | — |
| `:Terminated` | an event ended the solve | the last live cell deactivated (outcome `:burned_out`); a moving mesh thinner than `min_thickness` (`:min_thickness`); a cell below the deactivation threshold inward of a live one (`:interior`); a `terminate!` in a callback of yours (outcome `:running`) | §15.4 |
| `:Failure` | the solve's own checks ended it | more than 50 consecutive proposals rejected on a constraint residual; an inadmissible state, or no consistent node velocities, after a deactivation event; an inadmissible state accepted by a fixed-step solve | §15.3, §15.4.3 |
| `:DtLessThanMin` | the step fell below `dtmin` (default `1e-12` s) | every proposal rejected: a cell driven to absolute zero, a trajectory leaving the physical domain | §15.3.2 |
| `:MaxIters` | the step count reached `maxiters` (default `1_000_000`) | a step size that collapsed long before the end | §15.3 |
| `:Unstable`, `:DtNaN` | a non-finite state or step | a property law, rate or boundary function of yours returning `NaN` or `Inf` | §15.6.3 |

| `outcome` | Set when | `succeeded(sol)` |
|---|---|---|
| `:running` | neither the deactivation nor the thickness stop ended the solve (it reached the end, failed, or a callback of yours ended it) | `true` only with `:Success` |
| `:burned_out` | the last live cell deactivated | `true` |
| `:min_thickness` | on a moving mesh, the distance from the back node to the outermost live node fell below `min_thickness` | `true` |
| `:interior` | a cell below its threshold lies inward of a live cell | `false` |

`succeeded(sol)` is the test to use for a burnout run, which normally ends with
`:Terminated`; `completed(sol)` is the strict test that the solve reached the
end of its time span.

### 15.1.2 The solver statistics

`sol.statistics` is a `SolverStatistics`; `sol.solve_time` is the wall time of
the solve [s].

| Field | Counts | A large value points at |
|---|---|---|
| `accepted`, `rejected` | accepted and rejected steps | a stiff transient, or proposals the domain test rejects |
| `rhs_evaluations`, `jacobians`, `linear_solves` | residual evaluations, colored Jacobian evaluations, linear solves | the cost of the solve |
| `nonlinear_failures` | nonlinear-solver failures; always zero under the default Rosenbrock policy, which has no nonlinear iteration | a front the Newton of `fbdf_policy` or `kencarp4_policy` cannot converge on (§15.5) |
| `domain_rejections` | proposals the domain test rejected | §15.3.1 |
| `largest_constraint_residual` | the largest scaled residual of the algebraic rows over the accepted states; one or below passes; zero without algebraic rows or with a fixed step | §15.3.3 |
| `invocations` | deactivation callback invocations that deactivated cells | §15.4.1 |

`show(sol)` in the REPL prints the return code, the outcome, the saved times
and their range, the live cells, the invocations, the accepted steps and the
solve time.

### 15.1.3 The warnings `solve` prints

| Warning (abridged) | Meaning | See |
|---|---|---|
| `solve(): stopped at t = …: the last proposed state, at t = …, takes cell i to or below absolute zero` | the solve stopped while rejecting proposals that cool a cell to 0 K | §15.3.2 |
| `solve(): a constraint residual stayed above its tolerance for N consecutive steps at t = … (scaled residual r)` | the domain test gave up; the solve ends with `:Failure` | §15.3.3 |
| `solve(): the fixed step accepted at t = … is inadmissible …` | a fixed-step solve left the physical domain; `:Failure` | §15.3.4 |
| `solve(): cell i is below its threshold inward of a live cell at t = …; the solve ends with outcome :interior` | an interior void | §15.4.2 |
| `solve(): the state after deactivating cell(s) at t = … is inadmissible` | an event left the receiver outside its tolerances; `:Failure` | §15.4.3 |
| `solve(): no consistent node velocities after deactivating cell(s) at t = …` | the mesh closure has no solution after an event; `:Failure` | §15.4.3 |

The surface energy balance prints nothing when its Newton iteration does not
converge (§15.6.2).

### 15.1.4 Inspecting the last state

The state a solve ends on is always saved, whether or not `saveat` names its
time, and a solve that fails ends on the last state the integrator accepted.
`sol.t[end]` is therefore where the solve stopped, and the fields there show
why:

```julia
k = n_timesteps(sol)                 # the last saved state
sol.t[k]                             # where the solve stopped [s]
get_profile(sol, :T, k)              # cell temperatures [K], NaN in a dead cell
get_profile(sol, :m, k)              # component masses, NC × n_cells [kg]
get_cell_positions(sol; at = k)      # cell centres [m], NaN in a dead cell
sol.extras.live_cells[k]             # cells still live
```

The conservation audit separates a modelling problem from a numerical one.
With a `ConservationRecorder` as the `callback`, `conservation_audit(sol, rec)`
closes mass and energy over the accepted steps against the live inventories
and the deactivation events; a relative defect far above the tolerance points
at the integration, not at the model (Chapter 9; Technical Reference §16).

```julia
rec = ConservationRecorder()
sol = solve(problem; use_ale = true, deactivation = true, callback = rec)
audit = conservation_audit(sol, rec)
audit.mass.relative, audit.energy.relative
```

---

## 15.2 Errors before integration starts

Errors that stop `solve` before its first step are thrown as exceptions whose
message names the remedy.

### 15.2.1 Keywords

`solve` accepts exactly the keywords of its reference table (§16.7) and refuses
any other with an `ArgumentError`:

```julia
solve(problem; tstops = [10.0])
# ArgumentError: solve(): unknown keyword `tstops`
```

It also refuses an `energy_state` other than `:temperature` or `:enthalpy`, a
non-positive `min_thickness`, and a `band` outside `[threshold, 1)`.

### 15.2.2 Problem validation

Validation checks mesh tags, reaction indices, the initial field count and
time span. Reservoir fractions must name gases in the material and sum to
one; non-vacuum reservoirs need incoming temperature and composition.
Every gas-bearing state needs physical pore volume, including sealed gas.
Setup then checks its uncertainty margin, the matrix support path and any
endpoint domain. See §7.9 and §7.5 for the governing conditions.

---

### 15.2.3 What `setup` refuses

`solve` hands the problem to `setup`, which refuses, with an `ArgumentError`,
a combination it does not implement:

| Message (abridged) | Cause | Remedy |
|---|---|---|
| `material … has a lateral_shrinkage_law, which needs a moving mesh: solve with use_ale = true` | a variable cross section on a fixed mesh | `use_ale = true` |
| `a lateral shrinkage law needs planar coordinates`; `in-depth radiation needs planar coordinates` | a `Cylindrical()` or `Spherical()` mesh with a lateral law or `BEER_LAMBERT` | a planar mesh, or `SURFACE_ABSORPTION` |
| `the lateral shrinkage law must be positive and finite on [-0.001, 1.001]`; `… must equal one at zero conversion` | a law that vanishes or is not one at `χ̄ = 0` | a law such as `χ̄ -> 1 - 0.3χ̄` |
| `a correction moves the mesh relative to the solid, which needs use_ale = true` | a `correction` other than `NoCorrection()` on a fixed mesh | `use_ale = true` |
| `a prescribed mesh velocity needs use_ale = true`, and its siblings | `w_override` on a fixed mesh, on a reacting material, with `deactivation`, or with a `correction` | `w_override` is a verification hook (§16.7) |
| `the cells beside node(s) … start with different condensed components: give the layer interfaces as lagrangian nodes …` | an explicit `lagrangian` list that omits a layer interface under the equidistribution correction | leave `lagrangian` at its default, which lists the interfaces, or add them |
| `the mass conditions name the tag :x, which the mesh does not have` | a mass condition keyed by a tag the mesh lacks | key the conditions by the mesh's tags |

`integrate`, which `solve` calls next, refuses an inadmissible initial state ("a
live component mass below its tolerance, a live cell at or below absolute zero,
or crossed nodes"): a `T_initial` at or below 0 K, or a negative `ξ_initial`.

### 15.2.4 Constructor checks

The constructors reject a deck before a problem exists:

- **Property laws.** `SolidComponent`, `LiquidComponent` and `GaseousComponent`
  sample every law on the property bracket `[200, 3000]` K at 1 K steps: `ρ`,
  the intrinsic density, `c` and `k` must be positive (a gas may have `c = 0`),
  `λ` and `κ` non-negative and `ε` within `[0, 1]`. The message names the first
  failing temperature and suggests holding the fit constant outside its measured
  range, `T -> law(clamp(T, T_lo, T_hi))`, or `validate_properties = false`
  (§15.6.3).
- **Component keywords.** `absorptivity` and `γ` must lie in `[0, 1]`; a gas
  needs `M > 0`; `ρ_intrinsic` and `φ` exclude each other, and `φ` needs a
  constant `ρ`.
- **Material.** `gas_viscosity` must be positive and finite at 200, 300, 500,
  800, 1200 and 1800 K.
- **Reaction.** The order `n ≥ 0`, `gate_width > 0`, yields `≥ 0` summing to one
  (`validate_mass = false` drops the sum check only).
- **Mixing.** `MixingSpec(CARMAN_KOZENY)` needs a positive `length_scale`;
  `TableProperty` needs strictly increasing temperatures.

---

## 15.3 Domain rejections and a collapsing step

### 15.3.1 What the domain test rejects

The integrator hands every proposal its error test accepts to `DomainCheck`,
the domain test. A proposal is rejected, and the step cut by the controller's
smallest factor, when

- it is **inadmissible**: a live component mass below minus its absolute
  tolerance (`reltol · threshold · mˢᵢ(0)`, with `mˢᵢ(0)` the cell's initial
  condensed mass), a live cell at or below absolute zero, and on a moving mesh a
  cell whose nodes have crossed or a non-positive cross-section area;
- under **lumped gas transport** transport, a live cell's porosity
  `1 − Σⱼ mⱼ/(ρⱼ(T) V)` lies outside `[0, 1]` by more than the skeletal volume of
  its mass tolerances;
- an **algebraic row's** scaled residual exceeds one: the node-velocity balance
  rows on a moving mesh, the face intensities under `BEER_LAMBERT`, the cross
  section and its lateral rate under a lateral shrinkage law (§15.3.3).

`sol.statistics.domain_rejections` counts them. A few are normal on a fast
transient: the trial states of an implicit step leave the physical domain, and
the smooth guards of the state mapping keep the residual finite there
(Technical Reference §15). Many rejections together with a shrinking step mean
the trajectory itself is heading out of the domain.

### 15.3.2 `:DtLessThanMin` and the absolute-zero warning

The step falls below `dtmin` when the integrator cannot find a proposal that
both tests accept. When the last rejection was a cell at or below absolute
zero, `solve` names it:

```text
solve(): stopped at t = …: the last proposed state, at t = …, takes cell i to or below absolute zero
```

A cell reaches absolute zero only under a heat sink that does not weaken as it
cools: a `HeatFluxBC` that draws heat out at a fixed rate, or an endothermic
reaction whose rate does not fall with temperature (`E = 0`, or a heat of
reaction given in kJ/kg where J/kg is meant, a thousandfold sink). Inspect the
last saved state around the named cell (§15.1.4) and the boundary conditions
and heats of reaction that act on it. When the trajectory itself leaves the
domain, a smaller `dtmin` or other tolerances move the time at which the solve
stops, not whether it stops.

### 15.3.3 Constraint rejections and `:Failure`

The algebraic rows carry no time derivative: on a moving mesh the node
velocities satisfy the cells' volume balance, under `BEER_LAMBERT` the face
intensities satisfy the attenuation recurrence from the exposed face, and under
a lateral shrinkage law the cross section and its lateral rate satisfy their
closure. The domain test holds each row to `rtol · S + a`, with `S` the scale of
its own equation and `a` the absolute tolerance of its unknown, so tightening
`reltol` tightens this test with the error test.

After more than 50 consecutive constraint rejections
(`Pyrolysis.Integration.MAX_CONSTRAINT_REJECTIONS`) the solve ends with
`:Failure` and the warning of §15.1.3. A custom loop sets the bound with
`integrate(cp; max_constraint_rejections = …)` (Chapter 11). On a solve that
completes, `sol.statistics.largest_constraint_residual` reports the largest
scaled residual over the accepted states; one or below passes.

### 15.3.4 Fixed-step solves

With `adaptive = false` the integrator steps at `dt` and the domain test never
runs, since it reads only the proposals of an adaptive error test. Each
accepted state is held to admissibility instead, and the first inadmissible one
ends the solve with `:Failure` and the warning "the fixed step accepted at t = …
is inadmissible". `largest_constraint_residual` is zero for a fixed-step solve.
Fixed steps serve convergence studies of the time discretization; production
runs use the adaptive default.

---

## 15.4 Burnout runs

With `deactivation = true`, cells whose condensed mass falls below the fraction
`threshold` (`ε_m`, default `1e-2`) of their initial condensed mass deactivate,
and their inventories pass to the live cell beneath. Chapter 10 states the
rule, the band, the equidistribution correction, located events and the
supported envelope; this section covers what goes wrong.

### 15.4.1 Reading a burnout run

```julia
sol.statistics.invocations          # deactivation callback invocations
sol.extras.deactivation_times       # one time per deactivated cell [s]
sol.extras.live_cells               # live cells at each saved time
[(e.t, e.cell, e.receiver, e.invocation) for e in sol.events]
```

Each `DeactivationEvent` in `sol.events` records one dead cell: its time, the
invocation it belongs to, the cell, the receiver (zero when the last live cell
went), the condensed-mass fraction it held, and the masses and enthalpy the
event lost. The number of invocations per vanishing layer is reported, not a
criterion of correctness.

The integrator's interpolant on the step that ends at an event is not the
trajectory, so `sol(t)` and `interpolate_solution(sol, t)` refuse a time
strictly between two saved states that an event separates, with an
`ArgumentError` naming the event. Read such a time from the saved states, or
save it with `saveat`.

### 15.4.2 The interior void (`:interior`)

```text
solve(): cell i is below its threshold inward of a live cell at t = …; the solve ends with outcome :interior
```

Cells deactivate from the exposed side only, so a cell that empties beneath a
live cell ends the solve with `:Terminated`, outcome `:interior`, and
`succeeded(sol) == false`. The solution up to the stop is valid.

- **On a fixed mesh** only the outermost live cell can go. A layer that empties
  from within, as one does when the temperature maximum lies below the surface
  under in-depth absorption, ends the solve at the first interior void. The
  verification case V-BO-05 runs the PMMA-on-Kaowool deck of V-BO-01 with
  in-depth absorption on a fixed mesh and stops at 259.6 s with outcome
  `:interior`: the surface cell persists as a low-density skin while the cells
  beneath it empty. A material that regresses belongs on a moving mesh
  (`use_ale = true`), where its cells thin instead of emptying; V-BO-01 burns
  the same layer out on a moving mesh in one invocation.
- **On a moving mesh** the run of cells from the exposed side whose condensed
  mass is below the band `ε_band` goes together with the cell that crosses the
  threshold. The band must exceed what the slowest cell of a layer keeps when
  the first reaches the threshold; a narrower band, or none, ends the solve with
  `:interior`. The default band is `0.1` where the equidistribution correction
  moves a layer and the threshold otherwise, so `correction = NoCorrection()`
  with deactivation leaves no band beyond the threshold. Widen `band`, or lower
  `threshold`, which lowers what the slowest cell keeps at the first crossing
  (Chapter 10 gives the measured fractions).

### 15.4.3 Inadmissible or inconsistent states after an event

Two warnings end a burnout run with `:Failure`, the outcome left `:running`:

- "the state after deactivating cell(s) at t = … is inadmissible": the receiver,
  holding the merged inventories, fails admissibility against its own
  tolerances, or a lumped gas transport porosity leaves its range;
- "no consistent node velocities after deactivating cell(s) at t = …": the sweep
  of the balance rows leaves a velocity row above round-off.

Neither occurs on the decks of the supported envelope (Chapter 10). Report the
deck and the time the warning names.

### 15.4.4 The thickness stop

On a moving mesh the solve ends with `:Terminated` and outcome `:min_thickness`
when the distance from the back node to the outermost live node falls below
`min_thickness` (default `1e-6` m, which must be positive); `succeeded` counts it
as a finished solve. On a fixed mesh the thickness does not change and the stop
never fires. Without `deactivation`, a cell that gasifies completely thins
toward zero width and the step collapses as it approaches it; turn on
`deactivation`, or raise `min_thickness` to end the run first.

### 15.4.5 The end-of-layer burst

Just after a layer deactivates, the per-step mass-loss rate can spike: at 30
cells across the layer the largest per-step rate can exceed the main peak, up
to threefold on a mesh graded toward the surface, and it falls below the main
peak by 60 to 120 cells. The peak sampled at the output times is within
tolerance from 30 cells. Inside the supported envelope the burst converges
under refinement; a quantity read at every accepted step (with a callback)
needs the finer mesh (Chapter 10).

---

## 15.5 Closed channels and disappearing matrix support

Zero `λ` and zero `κ` are legal explicit channel choices. A positive-weight
zero in an exact series mixture closes that channel; a disappearing barrier
can make its support discontinuous. Setup checks reachable support changes
instead of adding a numerical conductance floor. Choose a continuous physical
mixing closure or a persistent barrier when the material trajectory requires it.
An exactly closed gas-bearing material still computes its pore EOS and enforces
physical pore volume. A matrix-free endpoint is separate data with a declared
skeletal-occupancy domain (§7.5).

---

## 15.6 Temperatures at the edges of the property bracket

### 15.6.1 The bracket and the trial temperature

Every property law is evaluated at its temperature argument clamped to the
property bracket `[200, 3000]` K (`T_PROPERTY_MIN`, `T_PROPERTY_MAX` in
`Pyrolysis.Materials`), and every temperature-sensitive kernel (Arrhenius
factors, gas densities, radiation) evaluates at the trial temperature, the
cell temperature smoothly clamped into the bracket with 1 K transitions. The
trial temperature is the cell temperature on every physical state; it keeps the
residual finite with a finite derivative on the integrator's trial states.

The domain test rejects only a cell at or below absolute zero. An accepted state
with a live cell above 3000 K is not rejected, but the model there evaluates
the cell at the bracket's end: its properties, its rates and its emission are
those at 3000 K. Check the live temperatures of a solve that heats strongly:

```julia
T_live = filter(isfinite, sol.T)       # dead cells are NaN
extrema(T_live)                        # inside (200, 3000) K
```

### 15.6.2 The surface balance at the trial clamp

The gas-free surface balance solves `q_bc(T_s) = k (T_s − T_cell)/d`.
Gas-bearing boundaries solve the mass and thermal equations together when
coupled, with qualified roots and implicit derivatives. A failed coupled
root raises a diagnostic; it is not accepted as a boundary flux. The saved
`get_boundary_exchange(sol, tag)` fields include its residual and iteration
count. Material properties use their declared trial-temperature bracket
(§5.4; Technical Reference §12).

On an accepted state, a surface temperature
within 1 K of either end of the bracket, or beyond it, means the condition was
evaluated at the bracket's end: the re-radiation and convection are those of
3000 K (or 200 K), not of `T_s`. Check the saved surface temperatures:

```julia
Ts = filter(isfinite, sol.surface_T)
all(T -> 201.0 <= T <= 2999.0, Ts)     # the surface balance stayed inside the bracket
```

A surface pinned at the clamp points at the deck: an incident flux on a face
with no loss term (no `RadiativeBC` or `ConvectiveBC`), a flux in the wrong
unit, or an exposed cell emptied of its solid on a fixed mesh without
`deactivation`, which conducts like its pore gas. The verification cases count
the accepted states whose surface balance did not converge or sits on the clamp
(V-BO-01's surface quantity); a custom loop reads the per-face state of the
thermal boundary kernel (`converged`, `iterations`, `T_surface`) through
`Pyrolysis.Internal` (Chapter 11).

### 15.6.3 Property laws outside their measured range

A law is evaluated on the whole bracket, not only where it was measured:

- `TableProperty` holds its end values outside its table; a table that stops at
  600 K returns its 600 K value at 900 K.
- `LinearProperty`, `PolynomialProperty` and `CallableProperty` extrapolate their
  formula to the bracket's ends. A fit that turns non-positive on the bracket is
  rejected by the component constructor; hold it constant outside its range
  instead, as the PMMA example does (`examples/06_cone_calorimeter.jl`):

  ```julia
  k_fit = CallableProperty(T -> T < 378.0 ? 0.15 : 0.34 - 3.9e-4 * min(T, 700.0))
  ```

- `validate_properties = false` skips the check, for laws built from search
  parameters where an unphysical value should reach the objective function; a
  law that is negative on the bracket then reaches the residual.

A law that returns `NaN` or `Inf` (a logarithm or a fractional power of a
negative argument) produces a non-finite state; §15.1.1 lists the return codes.

---

## 15.7 Concentrations near zero

The component masses are the state, and the face fluxes move mass between cells
exactly, so mass is conserved by construction. A component mass may dip below
zero by up to its absolute tolerance, `reltol · threshold · mˢᵢ(0)`; the domain
test rejects anything below that. Small negative values of that order in
`sol.m` or `sol.ξ` are expected, not an error.

- **The depletion limiter.** Each reactant's rate factor is exactly zero at and
  below zero concentration and C¹ through it, for every reaction order:
  `DepletionLimiter(; threshold = 1.0, join = 1.0, base_scale = 1e-8, enabled = true)`,
  the `threshold` [kg·m⁻³] the concentration scale of the roll-off. With
  `enabled = false` the factor is one, an absent reactant reacts at the rate of
  a `base_scale` concentration, and the admissibility check is the only guard
  (§4.10; Technical Reference §6.3).
- **Upwind advection.** The lumped gas transport advective gas flux is first-order upwind
  of the intrinsic gas density `ξ/φ`, stable and diffusive (Technical Reference
  §9).
- **A reactant initialized at zero.** A reaction consumes only what `ξ_initial`
  seeds; check that every reactant a scheme needs starts with mass.

The smooth temperature gates `½(1 + tanh((T − T_min)/w)) · ½(1 + tanh((T_max − T)/w))`,
`w` the reaction's `gate_width` (default 1 K), switch over about ten widths
around `T_min` and `T_max` (Technical Reference §6.2); a reaction that appears
to fire slightly early or late is that switch, not a defect.

---

## 15.8 Moving meshes

With `use_ale = true` the back node stays at rest and every other node moves
with the solid velocity plus the correction, `w = v_s + c` (Technical Reference
§11).

### 15.8.1 Crossed nodes and vanishing cells

A proposal whose nodes have crossed (a cell of non-positive width) or whose
cross-section area is non-positive is inadmissible, so the domain test rejects
it: a mesh that tangles shows as domain rejections and, when the trajectory
itself inverts a cell, as `:DtLessThanMin`. The node positions at the saved
times are `sol.z` (`NaN` beyond the outermost live node):

```julia
any(d -> d <= 0, filter(isfinite, diff(sol.z[:, k])))   # a crossed cell at saved state k
```

A cell that gasifies completely thins toward zero width; burnout belongs to
`deactivation` (§15.4). A swelling factor `γ` or a lateral shrinkage law that
drives a cell's volume through zero is a modelling error; `setup` already
refuses a lateral law that is not positive on the conversion range (§15.2.3).

### 15.8.2 The correction and the layer interfaces

Where the mesh moves, deactivation is on and some layer can vanish, `setup`
installs the equidistribution correction, which moves the nodes of a vanishing
layer relative to its solid so that its cells thin together; the node
corrections `c = w − v_s` at the saved times are `sol.correction_velocity`
(zero where no correction acts). The correction never moves a node between
cells that start with different condensed components: those layer interfaces
are Lagrangian nodes, listed by default. A layered deck should therefore put
nodes on its interfaces, so that each cell starts with one composition;
`create_mesh_from_nodes` builds such a mesh (Chapter 6, Chapter 10).

### 15.8.3 Mass-loss rate under a variable cross section

`sol.mass_loss_rate` is per unit **current** area of the `:top` faces, and
`sol.extras.MLR_gauge` per unit **initial** area (the cone-gauge convention).
They coincide on a constant cross section. An MLR that looks off by a factor
that tracks `A(t)/A_0` is being read from the wrong one.

---

## 15.9 Radiation surprises

The radiation model is the `radiation_model` keyword of `solve` (default
`NO_RADIATION`).

### 15.9.1 Model selection

| Model | An incident `RadiativeFluxBC` | Absorption coefficients `α_j` | Rows added |
|---|---|---|---|
| `NO_RADIATION` | its absorbed share enters the surface energy balance | not read | none |
| `SURFACE_ABSORPTION` | the same | not read | none |
| `BEER_LAMBERT` | absorbed in depth through the face intensities; it contributes nothing to the surface balance | read | one face intensity per node (algebraic); on a fixed mesh these make the system a mass-matrix DAE |

`NO_RADIATION` and `SURFACE_ABSORPTION` treat an incident flux identically: the
absorbed share `a·q` (the condition's `absorptivity`, or the exposed cell's
mixture absorptivity) enters the surface balance. `BEER_LAMBERT` needs planar
coordinates.

### 15.9.2 "Nothing absorbs", or "it never burns out"

`solve` prints no warning for a model that does not match the material:

- Under `BEER_LAMBERT` with every `α_j = 0` the beam crosses the slab and leaves
  through the back face. `sol.extras.volumetric_absorption` and
  `sol.extras.radiation_transmitted` show where the absorbed flux went (both per
  unit exposed area).
- Under `SURFACE_ABSORPTION` the `α_j` are not read at all; switch to
  `BEER_LAMBERT` for depth penetration.
- Under `BEER_LAMBERT`, a semi-transparent layer on a non-absorbing substrate or
  an insulated back absorbs a falling share of the beam as it thins, so its last
  fraction burns slowly; V-BO-03 measures the share that leaves through the back
  face of a transparent substrate rising toward one.

### 15.9.3 Optical-thickness regime

Choose the model by the optical thickness `τ = α L`:

| `τ = α·L` | Regime | Model |
|---|---|---|
| > 10 | optically thick | `SURFACE_ABSORPTION` |
| 0.1–10 | semi-transparent | `BEER_LAMBERT` |
| < 0.1 | transparent | `NO_RADIATION` |

`SURFACE_ABSORPTION` on a semi-transparent sample concentrates the absorbed flux
at the surface and overheats it; `BEER_LAMBERT` on an opaque char adds the
intensity rows and their coupling for no benefit (Technical Reference §8).

### 15.9.4 Time-varying incident flux

A flux given as a function of time (`RadiativeFluxBC(flux = t -> 25e3 * (1 - exp(-t / 10)))`)
is evaluated at every residual evaluation. `solve` has no list of times the
integrator must step to, so a flux feature shorter than the steps the error
test would take across it (a short pulse) can be stepped over: bound the step
with `dtmax` (§15.12.1). Under `BEER_LAMBERT` the face-intensity tolerances
scale with the largest incident flux over the time span, sampled at 101 times.

---

## 15.10 Gas pressure and surface-solve diagnostics

Use `get_pressure(sol)` and `get_porosity(sol)` for current pore states,
`get_boundary_exchange(sol,tag)` for signed species exchange and local residuals,
and `gas_species_audit(sol,recorder)` to separate generation, escape, ingress
and vented species. Do not rebuild boundary fluxes from initial cell widths.

A local `TransportFailure` identifies an unresolved interface or surface root.
Check coefficient support, complete reservoir data, physical pore volume and
the reported failure reason. Very strong nearly isobaric coupling may exceed
Float64 pressure resolution. A failed root is not an accepted boundary flux.
Pressure-aware error control tracks proposed temperature and pore geometry;
validate solution accuracy against a tighter solve, particularly during shrinkage.

---

## 15.11 Parameter sweeps and parallelism

`parameter_sweep(problem, grid; modify, extract = identity, parallel = nworkers() > 1, solve_kwargs...)`
solves one problem per grid point (Chapter 12). The contract that keeps it
correct and fast:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]

# `modify` returns a fresh problem and leaves the base problem alone; it may
# share the base problem's mesh, which a solve never changes.
@everywhere function modify(problem, p)
    m = problem.material
    reaction = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 2e5)
    material = Material(; name = m.name, components = m.components,
                        reactions = (reaction,), mixing = m.mixing)
    return PyrolysisProblem(; mesh = problem.mesh, material, bc_set = problem.bc_set,
                            T_initial = problem.T_initial, ξ_initial = problem.ξ_initial,
                            tspan = problem.tspan, output_times = problem.output_times)
end

# `extract` keeps the result small, and records how each solve ended.
@everywhere extract(sol) = (ok = succeeded(sol), peak_T = maximum(sol.surface_T))

results = parameter_sweep(problem, grid; modify, extract)
```

- **`modify` is required.** Without it `parameter_sweep` throws an
  `ArgumentError`. There is no copy-with-overrides constructor of `Material` or
  `PyrolysisProblem`: rebuild what changes from the base problem's fields.
- **Return small results.** `extract` defaults to `identity`, the whole
  solution; on workers that is costly to send back.
- **Workers need the definitions.** `modify`, `extract` and anything they call
  go to every worker with `@everywhere`.
- **Failures.** A solve that runs and fails returns its solution, so record
  `sol.retcode` or `succeeded(sol)` in `extract`; an exception in `modify` or
  before integration (a validation error) propagates and ends the sweep.

`parallel = false` forces a serial `map` for debugging. `solve` draws no
progress bar unless `progress = true`.

---

## 15.12 Performance

### 15.12.1 Checklist

Work top to bottom; each item is cheap to try and the early ones dominate.

1. **Cells.** Each cell carries its `NC` component masses and one energy
   unknown; a moving mesh adds a position and a velocity per node, `BEER_LAMBERT`
   one face intensity per node, and a lateral shrinkage law two history unknowns
   per cell and two global rows that read every cell. Use the smallest mesh
   that converges the quantities you need (Technical Reference §18.5 and
   `examples/convergence_douglas_fir.jl`), and grade cells toward the surface
   with `create_stretched_mesh(L, n, stretch)` rather than refining uniformly.
2. **Tolerances.** `reltol` (default `1e-6`) is the most direct speed knob: it
   scales every row's tolerance but the temperature's, and the constraint test
   with them. `abstol` (default `1e-8`) is the temperature tolerance in kelvin.
   `threshold` scales the mass tolerances (`reltol · threshold · mˢᵢ(0)`), so a
   lower threshold also tightens them. Loosen for exploration; tighten for
   verification.
3. **The integrator policy.** The default is `Rodas5P`, a Rosenbrock method
   that evaluates the colored Jacobian at every step and has no nonlinear
   iteration. `integrator_policy = fbdf_policy(StructuredKLU())` (variable-order
   BDF) and `integrator_policy = kencarp4_policy(StructuredKLU())` (an SDIRK
   method) are the alternatives. Which is fastest depends on the deck: time a
   representative solve with each (`sol.solve_time`, `sol.statistics`).
   Chapter 11 describes their restart and flag rules.
4. **The linear solver.** `StructuredKLU()`, the default, factors the
   integrator's `W` on the entries that have been nonzero at some step of the
   solve rather than on every entry the connectivity allows.
   `linear_solver = UMFPACKFactorization()` (`using LinearSolve`) is the
   alternative direct solver.
5. **`dtmax`.** The step-size ceiling (default `Inf`) does not speed a solve up.
   Set it when a time-dependent condition has features shorter than the steps
   the error test would take across them (§15.9.4).
6. **The radiation model.** `BEER_LAMBERT` adds one algebraic row per face, each
   intensity coupled to every cell outward of it. Do not pay for in-depth
   absorption the optical thickness does not call for (§15.9.3).
7. **Observers.** A `ConservationRecorder` evaluates the residual four times per
   accepted step (three Gauss–Legendre points and the inventories after the
   step) in a workspace of its own: use it for audits, not for production
   sweeps. Saving more times does not change the steps; each saved state costs
   one residual evaluation when `solve` reconstructs the fields. `progress = true`
   draws a progress bar at the cost of one update per accepted step.
8. **Parallelise sweeps** (§15.11) and **keep one Julia session for many solves**
   (§15.12.2).

### 15.12.2 Compilation and precompile workloads

Integrator code is shared across problems with the same solver configuration.
Kernels and colored Jacobians compile for concrete material, boundary,
transport and energy-state types. Changing numerical values within these
types reuses compiled methods; creating a new callable property type can
require compilation.

`src/precompile.jl` exercises representative fixed and moving meshes,
radiation, deactivation, lumped and dusty gas, and both energy formulations.
`test/integration/test_precompile_workloads.jl` executes these workloads.
Separate first-call compilation from warmed timing, and compare complete
solves at matched temperature, pressure and species accuracy. The pressure
error budget and local nonlinear work are part of that comparison.

---

## 15.13 Where to look next

- The `solve` keywords and what each does: Chapter 8; every keyword and default:
  Chapter 16.
- The solution, its accessors and the conservation audit: Chapter 9; the ledger
  and what it proves: Technical Reference §16.
- Burnout, the band, the correction and the supported envelope: Chapter 10;
  deactivation as implemented: Technical Reference §14.
- The integrator policies, the linear solver, `setup` and `integrate` for custom
  loops, and `Pyrolysis.Internal`: Chapter 11; the time integration, the domain
  test and the colored Jacobian: Technical Reference §15.
- Sweeps and finite-difference sensitivities: Chapters 12 and 13.

The escape hatch for anything not on the curated top-level surface is the
`Pyrolysis.Internal` umbrella (`using Pyrolysis.Internal`), which re-exports
every submodule's public symbols: the residual, the colored Jacobian, the state
layout, the kernels and the closures. Symbols reached that way are not part of
the stable public API and may move between minor releases.
