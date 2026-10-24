# 8. Running the Solver

Once you have a `PyrolysisProblem` (Chapters 5–7), the single function that
integrates it forward in time is `solve`. This chapter is the complete option
reference for `solve(problem; ...)`: how to pick the integrator policy and the
linear solver, how the tolerances reach every row of the state, what the
Jacobian is, how to choose the energy state, how to switch on the optional
physics (a moving mesh, radiation, cell deactivation), how callbacks and the
conservation recorder are wired in, how output is saved, and how to read the
return code and the outcome. Every keyword shown here is a real `solve`
keyword: `solve` refuses an unknown keyword with an `ArgumentError`, so a typo
never passes silently. The physics behind these options is derived in the
Technical Reference; this chapter tells you which knob to turn and what each
one does. The temporal integration, the differential-algebraic system and the
Jacobian are derived in Technical Reference §15.

---

## 8.1 The `solve` entry point

`solve` takes the problem as a positional argument and everything else as
keywords:

```julia
using Pyrolysis

solution = solve(problem)   # all defaults: Rodas5P, StructuredKLU, fixed mesh, no radiation
```

That call validates the problem (`validate_problem`), builds the problem the
integrator runs (`setup`: the initial inventories at the cell centres, the
kernels, the per-row tolerances and the colored Jacobian), integrates it
(`integrate`) and reconstructs a `PyrolysisSolution` at the saved times. The
two halves are exported, so a script can set a problem up once and integrate
it several times (§8.11). `solve` leaves `problem` unchanged: the mesh is never
modified, every node position the solve moves lives in the state vector, and
the solution holds the geometry at each saved time (`solution.z`,
`solution.cell_centers`). The same problem can therefore be solved again, or
shared by the points of a parameter sweep, and always starts from the mesh you
built. The problem is kept as `solution.problem`.

`solve` is a method of the SciML ecosystem's common `solve` function, so
`using Pyrolysis` together with `using OrdinaryDiffEq` names one function and
raises no ambiguity.

The full keyword set, with defaults, is:

| Keyword | Default | What it controls |
|---|---|---|
| `integrator_policy` | `nothing` → `default_policy(linear_solver)` (`Rodas5P`) | Integrator with its restart and flag rules; §8.2 |
| `linear_solver` | `nothing` → `StructuredKLU()` | Linear solve of the default policy; §8.4 |
| `abstol` | `1e-8` | Temperature tolerance [K]; the other rows derive theirs (§8.3) |
| `reltol` | `1e-6` | Relative tolerance of every row; §8.3 |
| `threshold` | `1e-2` | `ε_m`, the deactivation threshold as a fraction of a cell's initial condensed mass; also scales the mass, position and velocity tolerances; §8.3, §8.9 |
| `maxiters` | `1_000_000` | Maximum number of integrator steps; §8.3 |
| `dt` | `problem.dt_initial` | First step [s] (the fixed step with `adaptive = false`); also the time scale of the node-velocity tolerances; §8.3 |
| `adaptive` | `true` | Adaptive step-size control; §8.3 |
| `dtmin` | `1e-12` | Minimum step [s]; §8.3 |
| `dtmax` | `Inf` | Maximum step [s]; §8.3 |
| `energy_state` | `:temperature` | Energy variable per cell, `:temperature` or `:enthalpy`; §8.6 |
| `use_ale` | `false` | Moving mesh; §8.7 |
| `min_thickness` | `1e-6` | End a moving-mesh solve when the live thickness falls below it [m]; §8.7, §8.9 |
| `correction` | `nothing` → by rule | Mesh correction (`NoCorrection()`, `EquidistributionCorrection(...)`); §8.7, Chapter 10 |
| `lagrangian` | `nothing` → the layer interfaces | Node ids that always move with the solid; §8.7 |
| `w_override` | `nothing` | Prescribed mesh velocity `w(x, t)` over a solid at rest (verification); §8.7 |
| `radiation_model` | `NO_RADIATION` | `NO_RADIATION`, `SURFACE_ABSORPTION` or `BEER_LAMBERT`; §8.8 |
| `deactivation` | `false` | Deactivate burned-out cells during the solve; §8.9 |
| `band` | `nothing` → 0.1 under a moving correction, else `threshold` | `ε_band`, the near-threshold band; §8.9 |
| `located_events` | `true` | Locate each deactivation inside the step; §8.9 |
| `callback` | `nothing` | Your callbacks and `ConservationRecorder`s; §8.10 |
| `progress` | `false` | Draw a progress bar; §8.10 |
| `saveat` | `problem.output_times` | Saved times (a vector, or a spacing); §8.12 |
| `validate` | `true` | Run `validate_problem` before solving; §8.11 |
| `sources` | `NoSources()` | External volumetric sources (verification); §8.11 |

The authoritative list is `Pyrolysis.Integration._SOLVE_KEYWORDS`; anything
else is refused before any work is done:

```julia
solve(problem; jacobian = 1)
# ERROR: ArgumentError: solve(): unknown keyword `jacobian`
```

`solve` always returns a `PyrolysisSolution`; see Chapter 9 for its fields and
accessors.

---

## 8.2 Choosing the integrator policy

Pyrolysis problems are **stiff**: Arrhenius rates change by orders of magnitude
across a narrow temperature band, and conduction over fine cells imposes a
small explicit stability limit. With a moving mesh or in-depth radiation they
are also differential-algebraic: the node velocities and the face intensities
are algebraic rows of the state, so the integrator receives a mass-matrix
system. `solve` hands the integrator that system together with the colored
Jacobian (§8.5), the per-row absolute tolerances (§8.3), a state whose
algebraic rows are already consistent (so the integrator runs no
initialization of its own) and the domain test (§8.3).

The integrator is chosen through an **integrator policy**: the
OrdinaryDiffEq algorithm together with two rules that belong to its contract.

- **The restart rule** fixes the step proposed after an event that modifies the
  state (a deactivation, §8.9): `:solver_chosen` keeps the integrator's own
  proposal, `:forced_dt` proposes `DT_AFTER_EVENT = 1e-6` s, and `:heuristic`
  proposes the integrator's initial-step estimate at the state the event left
  (two right-hand-side evaluations).
- **The observer-flag rule** fixes what happens to the integrator's
  derivative-discontinuity flag, which the callback machinery raises before
  every callback: `:keep` leaves it raised, `:clear` clears it. Under either
  rule a read-only observer (the progress bar, the conservation recorder, a
  discrete callback of yours) leaves the steps unchanged (§8.10).

Four policy constructors are exported. Each takes the linear solver as its
argument:

| Constructor | Method | Restart | Observer flag |
|---|---|---|---|
| `default_policy(ls)` | `rodas5p_policy(ls)` | `:heuristic` | `:clear` |
| `rodas5p_policy(ls; restart = :heuristic)` | `Rodas5P`, a Rosenbrock method | any of the three | `:clear` |
| `fbdf_policy(ls)` | `FBDF`, variable-order backward differentiation | `:forced_dt` | `:clear` |
| `kencarp4_policy(ls)` | `KenCarp4`, an ESDIRK method, raw error estimate | `:solver_chosen` | `:keep` |

```julia
solve(problem)                                                     # default_policy(StructuredKLU())
solve(problem; integrator_policy = fbdf_policy(StructuredKLU()))
solve(problem; integrator_policy = kencarp4_policy(StructuredKLU()))
solve(problem; integrator_policy = rodas5p_policy(StructuredKLU(); restart = :forced_dt))
```

The policies are exported by `using Pyrolysis`; no `using OrdinaryDiffEq` is
needed to select one.

How the methods differ in use:

- **`Rodas5P`** (the default) is linearly implicit: it evaluates the colored
  Jacobian at every step and solves linear systems, with no Newton iteration
  and therefore no Newton failures. Its heuristic restart rejects fewer steps
  after an event than keeping the step and keeps the mass ledger closer than
  either other restart. The verification suite runs on it, and the precompile
  workloads compile it, so it is also the policy with the shortest time to
  first solve.
- **`FBDF`** reuses its Jacobian across steps and iterates Newton on each. It
  restarts at order one after an event, hence its forced small restart step;
  with the flag raised it would restart at order one on every step, hence the
  flag is cleared.
- **`KenCarp4`** iterates Newton on each stage, so the Newton-failure step cap
  of §8.3 applies to it, as it does to `FBDF`.

In every policy the time derivative of the right-hand side, where the method
needs one, is a finite difference; the Jacobian with respect to the state is
the colored one. A policy other than the default compiles its integrator on
first use in a session. Chapter 11 describes building an `IntegratorPolicy`
around another OrdinaryDiffEq algorithm, and Technical Reference §15 the
methods and the matrix `W` they factor.

---

## 8.3 Tolerances and step control

The integrator adapts its step to keep the local error estimate of every row
within `abstol_row + reltol · |u_row|`. `reltol` applies to every row; the
absolute tolerances are set **per row** by `setup` from `abstol`, `reltol`,
`threshold` and each cell's initial condensed mass `m_s,i(0)`:

| Row | Absolute tolerance |
|---|---|
| Temperature (temperature form) | `abstol` [K] |
| Enthalpy inventory (enthalpy form) | `threshold · m_s,i(0) · c_p,ref · abstol` [J], so a cell at the threshold resolves its temperature to `abstol` |
| Component mass of cell `i` | `reltol · threshold · m_s,i(0)` [kg] (every component of the cell, a product included) |
| Node position (moving mesh) | `reltol · threshold · h_n`, `h_n` the smallest initial size of the cells around node `n` |
| Node velocity (moving mesh) | the position tolerance over the time scale `dt` |
| Face intensity (`BEER_LAMBERT`) | `reltol` times the largest incident intensity over the time span |
| Cross-section area and its lateral rate (lateral shrinkage) | `reltol · A₀`; the velocity tolerance over the initial thickness |

A cell that starts with no condensed mass takes its total initial mass for
`m_s,i(0)`, or 1 kg/m³ times its volume when it starts empty. Three
consequences follow:

- **`abstol` is a temperature tolerance in kelvin.** It is never compared with
  a concentration or a mass; those rows take their tolerances from `reltol`
  and `threshold`.
- **The mass tolerances follow `reltol` and `threshold`.** A cell keeps its
  masses resolved to the relative accuracy `reltol` of what it holds at the
  deactivation threshold, so tightening `reltol` tightens the mass rows too.
  Changing `threshold` changes the tolerances even when `deactivation` is off.
- **`dt` scales the node-velocity tolerances** on a moving mesh, besides being
  the first step.

The defaults are:

```julia
solve(problem; abstol = 1e-8, reltol = 1e-6)
```

Guidance:

- **Loosen** `reltol` (to `1e-5` or `1e-4`) for exploratory runs and parameter
  sweeps where you need trends.
- **Tighten** (`reltol = 1e-7`, `abstol = 1e-9`) for convergence studies and
  comparisons against reference codes or exact solutions. The Douglas-fir
  example (`examples/07_douglas_fir_cone.jl`) uses `abstol = 1e-9,
  reltol = 1e-7`.
- Drive accuracy through `reltol`. On a temperature row the error allowance is
  `abstol + reltol · T`, and at pyrolysis temperatures `reltol · T` (3e-4 K at
  300 K with the default `reltol`) dominates the default `abstol`.

### Step-size limits

```julia
solve(problem;
      dt       = 1e-4,        # first step [s] (default problem.dt_initial)
      dtmin    = 1e-14,       # allow tiny steps through a rapid reaction front
      dtmax    = Inf,         # no upper cap (default)
      maxiters = 2_000_000)
```

If a solve stops with `retcode == :DtLessThanMin`, read the warning `solve`
prints with it (below) before lowering `dtmin`. If it stops with `:MaxIters`,
raise `maxiters` or loosen the tolerances. `dtmax` caps the step; use it when
a boundary condition changes abruptly (a step or a short pulse in a flux
function), which an adaptive integrator can otherwise step over.

With `adaptive = false` the solve takes fixed steps of `dt`. The domain test
below does not run on fixed steps; instead the solve ends with return code
`:Failure` at the first accepted state that leaves the physical domain.

### The domain test

Besides its error test, the integrator applies a **domain test**
(`DomainCheck`) to every proposal its error test accepts. A proposal is
rejected, and the step cut by the controller's smallest factor, when

- it is inadmissible: a live component mass more negative than its absolute
  tolerance, a live cell at or below absolute zero, or on a moving mesh
  crossed nodes (a non-positive cell volume) or a non-positive cross-section
  area;
- under Darcy–Fick transport a live cell's porosity leaves `[0, 1]` by more
  than the skeletal volume of its mass tolerances;
- an algebraic row (a node velocity, a face intensity, the cross-section area
  or its lateral rate) misses its own equation by more than `reltol` times the
  equation's scale plus the unknown's absolute tolerance.

The test is part of every policy and has no keyword. Its counters are in
`solution.statistics`: `domain_rejections` (the proposals it rejected) and
`largest_constraint_residual` (the largest scaled algebraic-row residual of the
accepted states, where one or below passes; zero for a fixed-step solve). After
more than `MAX_CONSTRAINT_REJECTIONS = 50` consecutive rejections on an
algebraic row the solve ends with return code `:Failure` and a warning naming
the residual. When a solve stops while its last proposal took a cell to or
below absolute zero, `solve` warns

```
solve(): stopped at t = …: the last proposed state, at t = …, takes cell 12 to or below absolute zero
```

which names the cell whose energy balance the step could not keep positive;
Chapter 15 discusses the causes.

### The Newton-failure step cap

Policies that solve their stages by Newton iteration (`FBDF`, `KenCarp4`) carry
a step cap: when an accepted step needed Newton retries, the next
`FAILURE_MEMORY_STEPS = 20` accepted steps are held at or below
`FAILURE_MEMORY_GROWTH = 1.25` times that step, so the controller cannot regrow
the step straight back into the range where Newton failed. A run without
Newton failures is unaffected, and `Rodas5P` has none. Compare
`solution.statistics.nonlinear_failures` with `solution.statistics.accepted`
to see how often Newton fails.

---

## 8.4 The linear solver

Each implicit step solves sparse linear systems in the matrix `W = J − M/(γΔt)`
built from the Jacobian `J` and the mass matrix `M`. The default linear solver
of every policy is `StructuredKLU()`: KLU on the entries of `W` that have been
nonzero at some step of the solve. The Jacobian's prototype holds every entry
the residual's connectivity allows, most of which are zero at any one step;
`StructuredKLU` keeps the union of the entries that have been nonzero so far,
copies each fresh `W` into a matrix of that structure, and refactors it
numerically, analysing the structure afresh only when the union grows. The
union only grows, so the analyses stop once every entry the solve reaches has
appeared.

Two other sparse direct factorizations from LinearSolve.jl are accepted and
tested: `UMFPACKFactorization()` and `KLUFactorization()` (plain KLU on `W`'s
own structure). Pass one either through `linear_solver`, which configures the
default policy, or as the argument of a policy constructor:

```julia
using LinearSolve: UMFPACKFactorization

solve(problem; linear_solver = UMFPACKFactorization())                     # default policy
solve(problem; integrator_policy = fbdf_policy(UMFPACKFactorization()))   # another policy
```

A policy carries its own linear solver, so `linear_solver` is read only when
`integrator_policy` is not given. `StructuredKLU` is exported by
`using Pyrolysis`; the LinearSolve factorizations are imported from
LinearSolve. Under each of the three, repeating a solve gives bit-identical
results.

---

## 8.5 The Jacobian

The Jacobian of the residual is assembled by **colored sparse forward-mode
automatic differentiation** on a sparsity pattern built from the mesh
connectivity (`ColoredJacobian`): every row reads its cell, its faces and the
neighbours they touch, the pattern holds exactly those entries, a greedy
column coloring groups columns that have no row in common, and one
dual-number evaluation of the residual per batch of up to 24 colors fills the
matrix (the two dense rows of a variable cross section are assembled from
column sums afterwards). It covers every configuration
exactly, to round-off: the reactions, conduction, gas transport, the boundary
balances, the node velocities and the correction on a moving mesh, the face
intensities under `BEER_LAMBERT`, the cross-section rows under a lateral
shrinkage law, and the identity rows of deactivated cells. There is no Jacobian
keyword and no strategy to choose; every policy receives this Jacobian, and
`Rodas5P` evaluates it at every step.

The test suite holds the colored Jacobian against a dense forward-mode
reference (`dense_jacobian`). You can run the same check on your own problem at
its initial state through `Pyrolysis.Internal`:

```julia
using Pyrolysis.Internal: jac_prototype, jacobian!, dense_jacobian, compare_entries

cp = setup(problem; use_ale = true, radiation_model = BEER_LAMBERT)  # the solve's set-up keywords
t0 = cp.problem.tspan[1]
J = jac_prototype(cp.jacobian)
jacobian!(J, cp.jacobian, cp.u0, t0)
J_ref = dense_jacobian(cp.def, cp.ws, cp.u0, t0)
compare_entries(J, J_ref)    # (passed = true, max_abs_dev = …, worst_ratio = …, worst_index = …)
```

`compare_entries` tests each entry against a floor of its row's largest
magnitude, so a missing small block cannot hide behind a large one. Chapter 11
describes these building blocks; Technical Reference §15 derives the pattern
and the coloring.

---

## 8.6 The energy state

Each cell carries its component masses and one energy variable. The keyword
`energy_state` selects which:

```julia
solve(problem)                            # energy_state = :temperature (default)
solve(problem; energy_state = :enthalpy)  # the sensible-enthalpy inventory
```

| `energy_state` | Per-cell unknown | Conservation |
|---|---|---|
| `:temperature` | the cell temperature `Tᵢ` | energy conserved to the integrator's order |
| `:enthalpy` | the sensible-enthalpy inventory `Hᵢ = Σⱼ mᵢⱼ hⱼ(Tᵢ)` | a linear invariant every Runge–Kutta step conserves to round-off; `Tᵢ` recovered by a regularized caloric inversion |

Both forms solve the same equations and agree within the tolerances; the
solution reports temperatures either way. The temperature form is the default
and the one most configurations are precompiled for; the enthalpy form suits
energy-ledger studies where the energy closure must not depend on the
tolerance. Any other symbol is refused:

```julia
solve(problem; energy_state = :energy)
# ERROR: ArgumentError: energy_state must be :temperature or :enthalpy, got :energy
```

The energy equation behind both forms is described in Chapter 7 and derived in
Technical Reference §3 and §7.

---

## 8.7 Moving-mesh integration

By default every node is pinned and the mesh never moves. For charring,
ablating, swelling or intumescent materials whose surface recedes or expands,
`use_ale = true` moves the mesh with the material:

```julia
solve(problem; use_ale = true)
```

With `use_ale = true`:

- The node positions join the state, and the node velocities with them as
  algebraic rows (the balance rows): the back node (`z = 0`) stays at rest,
  and every other node moves with the solid velocity plus a correction,
  `w = v_s + c`. The solution reports the node positions in `solution.z`, the
  corrections `c` in `solution.correction_velocity` and the live thickness in
  `solution.extras.thickness_history`.
- The solve ends with return code `:Terminated` and outcome `:min_thickness`
  when the live thickness (from the back node to the outermost live node)
  falls below `min_thickness` (default `1e-6` m). Set it to your physical
  stopping criterion, for example `min_thickness = 0.5e-3`. It must be
  positive, and it is read on a moving mesh only.
- If the material carries a lateral-shrinkage law, the conversion-history pair
  and the cross-section area with its lateral rate join the state, and the
  area `A(t)` varies (§8.12, Technical Reference §10). A lateral-shrinkage law
  requires `use_ale = true` and planar coordinates; `solve` refuses it
  otherwise.

A typical opaque charring cone run:

```julia
solve(problem;
      use_ale         = true,
      radiation_model = SURFACE_ABSORPTION,
      min_thickness   = 0.5e-3)          # stop below 0.5 mm
```

Three keywords shape the motion:

- **`correction`** sets how the mesh moves relative to the solid. With the
  default `nothing`, `setup` installs the equidistribution correction
  (`EquidistributionCorrection()`) where the mesh moves, `deactivation` is on
  and some cell can vanish (its components, converted completely, leave less
  condensed mass than `threshold` of it), and `NoCorrection()` (the mesh moves
  with the solid, `c = 0`) otherwise. Pass either explicitly to override the
  rule. The correction exists to make a burning layer's cells thin together
  and leave together; Chapter 10 describes it and its parameters. A
  correction other than `NoCorrection()` requires `use_ale = true`.
- **`lagrangian`** lists node ids that always move with the solid, beside the
  mesh's own. Under an equidistribution correction every interior node
  between cells that start with different condensed components (a layer
  interface) must be among them, and by default they are exactly those nodes;
  a list that leaves one out is refused with an `ArgumentError` naming it.
- **`w_override`**, a function `w(x, t)` [m/s] of the node position and time,
  prescribes every node's velocity instead, moving the mesh over a solid at
  rest. It is a verification hook: the material must not react, no other
  correction may be given, and cells are not deactivated under it.

For fixed-mesh problems (`use_ale = false`) the thickness is constant until
cells are deactivated, and `solution.z` is a `0 × 0` matrix. Mesh-motion
theory and the geometric conservation law are in Technical Reference §11.

---

## 8.8 Radiation model selection

The `radiation_model` keyword selects where the incident radiant flux of a
`RadiativeFluxBC` (Chapter 5) is absorbed:

| Value | Physics | State and Jacobian |
|---|---|---|
| `NO_RADIATION` (default) | The absorbed incident flux enters the surface energy balance; no in-depth absorption | no extra rows |
| `SURFACE_ABSORPTION` | The same: the absorbed flux `a·q` enters the surface balance (opaque materials); component absorption coefficients are not read | no extra rows |
| `BEER_LAMBERT` | Exponential in-depth absorption through the cells' absorption coefficients (semi-transparent materials) | one algebraic intensity row per face; the Jacobian stays banded |

```julia
solve(problem; radiation_model = SURFACE_ABSORPTION)   # opaque char, wood
solve(problem; radiation_model = BEER_LAMBERT)         # semi-transparent polymer
```

`NO_RADIATION` and `SURFACE_ABSORPTION` treat a `RadiativeFluxBC` the same way;
name `SURFACE_ABSORPTION` when the incident flux is the physics you mean. The
absorptivity `a` is the condition's own `absorptivity` when given, otherwise
the effective absorptivity of the exposed cell (Kirchhoff's law unless a
component sets its own).

Under `BEER_LAMBERT`:

- The intensity entering the exposed face is `a·q(t)`; each cell absorbs the
  exact cell average of the attenuated beam, with the attenuation from the
  cells' absorption coefficients `α` (m²/kg, at the cell's concentrations).
  What leaves through the back face is reported in
  `solution.extras.radiation_transmitted`, what the cells absorbed in
  `solution.extras.volumetric_absorption`. A material with no absorption
  coefficient lets the beam cross the slab unabsorbed.
- The incident flux enters through the `:top` faces only. A
  `RadiativeFluxBC` on any other boundary is refused:
  `ArgumentError: in-depth radiation enters through the :top faces only; the
  condition on :bottom carries an incident radiative flux`. Use
  `SURFACE_ABSORPTION` or a `HeatFluxBC` for radiant heating of the back face.
- The model needs planar coordinates; on a cylindrical or spherical mesh it is
  refused.
- After surface cells are deactivated, the beam enters the first live cell
  through its upstream face.

### Time-varying radiation flux is always tracked

The incident flux `q(t)` of a `RadiativeFluxBC` is evaluated at every
right-hand-side evaluation, so a ramped, delayed or modulated cone flux is
tracked whenever it starts varying. Make the `flux` a function of `t` for a
time-dependent source:

```julia
top_bc = RadiativeFluxBC(flux = t -> 25e3 + 250.0 * t, absorptivity = 0.9) +
         RadiativeBC(ε = 0.9, T_inf = 300) +
         ConvectiveBC(h = 10, T_inf = 300)
```

Under `BEER_LAMBERT` the face-intensity tolerances scale with the largest
incident intensity over the time span (§8.3). Re-radiation losses
(`RadiativeBC`) and convective losses (`ConvectiveBC`) are ordinary thermal
boundary conditions, applied under every `radiation_model`. The radiation
physics and the surface energy balance are derived in Technical Reference §8
and §12.

---

## 8.9 Deactivation and termination

When a layer burns out, its cells empty (and on a moving mesh, thin) toward a
state whose equations the integrator cannot resolve. **Deactivation** removes
such a cell from the live mesh: when its condensed mass falls below the
fraction `threshold` (`ε_m`) of its initial condensed mass, it goes, its
inventories pass to the live cell beneath it, its rows become identity rows,
and the boundary conditions move to the face of the outermost live cell. The
state length and the Jacobian's pattern never change.

```julia
solve(problem;
      use_ale       = true,
      deactivation  = true,      # deactivate burned-out surface cells
      threshold     = 1e-2,      # ε_m: 1 % of the cell's initial condensed mass (default)
      min_thickness = 0.01e-3)   # end when the live thickness < 0.01 mm
```

The keywords:

| Keyword | Default | Effect |
|---|---|---|
| `deactivation` | `false` | Deactivate cells during the solve |
| `threshold` | `1e-2` | `ε_m`; a cell that starts with no condensed mass never deactivates |
| `band` | `max(threshold, 0.1)` when the equidistribution correction moves a layer, else `threshold` | `ε_band` in `[threshold, 1)`: on a moving mesh, the run of cells from the exposed side whose condensed mass is below this fraction goes together with a cell that crosses `ε_m` |
| `located_events` | `true` | Locate each event inside the step at the first time a live cell reaches its threshold, and apply it at the interpolated state there; at the end of every step, cells below the threshold go in one event as a safety net (the only events with `located_events = false`) |
| `correction` | by rule (§8.7) | On a moving mesh with `deactivation`, the equidistribution correction makes a layer's cells thin together and leave in one event |

How the two meshes differ:

- **Moving mesh.** The run of cells goes in one merge into the live cell
  beneath it, which also absorbs their volume; mass and sensible enthalpy are
  transferred, and nothing is lost.
- **Fixed mesh.** The outermost live cell goes while it is below the threshold.
  The receiver keeps its volume and takes the condensed inventories; the dead
  cell's gas is vented, and its masses and sensible enthalpy are recorded on
  the event as lost (`solution.events`, Chapter 9), so the conservation audit
  still closes (§9.8). A material that recedes belongs on the moving mesh.

The ways a deactivating solve ends:

| Outcome (`solution.outcome`) | Return code | Meaning |
|---|---|---|
| `:running` | `:Success` | Reached the end of the time span; cells may have been deactivated on the way |
| `:burned_out` | `:Terminated` | The last live cell went |
| `:min_thickness` | `:Terminated` | A moving mesh thinned below `min_thickness` |
| `:interior` | `:Terminated` | A cell below the threshold lies inward of a live cell; `solve` warns, naming the cell |

The interior case is the material gasifying from within (for example a slab
heated from its back face): deactivation removes cells from the exposed side
only, so `solve` stops with the warning
`solve(): cell 1 is below its threshold inward of a live cell at t = …; the
solve ends with outcome :interior`. A band narrower than what the slowest cell
of a layer keeps when the first reaches the threshold also ends a solve this
way; Chapter 10 states the band the verified decks need.

An event can also end the solve with return code `:Failure`, with a warning:
when the state after the merge is inadmissible (the receiver held to its own
mass tolerances, as the domain test holds every later proposal), or when no
consistent node velocities exist after it.

Deactivated cells read `NaN` in `solution.T`, `solution.ξ`,
`solution.cell_centers`, and beyond the outermost live node in `solution.z`;
`solution.extras.live_cells` gives the live count at each saved time,
`solution.extras.deactivation_times` the time of each deactivated cell, and
`solution.events` the full record of every event.

`examples/06_cone_calorimeter.jl` uses `deactivation = true` with the default
threshold. `examples/07_douglas_fir_cone.jl` leaves deactivation off and stops
on `min_thickness`. Chapter 10 is the full treatment: the threshold, the band,
the correction and its layer time, event location, the outcomes and the
supported envelope of the burnout path. Technical Reference §14 derives the
transfer and the event.

---

## 8.10 Callbacks and observers

`solve` assembles its own callback set: a per-step observer (it records the
domain test's residuals, ends the solve when the domain test gave up, applies
the Newton-failure step cap and the policy's flag rule), the deactivation
events (`deactivation`), the layer-time event of the equidistribution
correction (moving mesh), the thickness stop (moving mesh), the readers of any
conservation recorders, the progress bar (`progress = true`), your callback,
and last a recorder of the accepted state. You supply only the last two kinds
of input, through `callback` and `progress`.

### Your callback

`callback` accepts any OrdinaryDiffEq callback (a `DiscreteCallback`, a
`ContinuousCallback` or a `CallbackSet`), a `ConservationRecorder`, or a tuple
of recorders and callbacks.

- **A discrete callback is a read-only observer.** It may read the integrator
  (`integrator.u`, `integrator.t`, `integrator.dt`), record what it sees, and
  call `terminate!`; it must not modify the state. `solve` applies the policy's
  flag rule after its `affect!`, so adding it leaves the trajectory unchanged
  under every policy: the steps, the saved states and the statistics are
  identical with and without it.
- **A continuous callback is an event.** Its root is located inside the step
  and its `affect!` runs there; the integrator re-evaluates its derivative at
  the event as at any event, so the steps after it differ from a solve without
  it.

Recording the accepted step sizes, a read-only observer:

```julia
using OrdinaryDiffEq: DiscreteCallback, terminate!   # the callback types

dts = Float64[]
record_dt = DiscreteCallback((u, t, integrator) -> true,
                             integrator -> push!(dts, integrator.dt);
                             save_positions = (false, false))
solve(problem; callback = record_dt)
```

Stopping a run after a wall-clock budget, for example in a sweep:

```julia
started = time()
budget = DiscreteCallback((u, t, integrator) -> time() - started > 600, terminate!)
solution = solve(problem; callback = budget)
```

A solve your callback terminates returns `:Terminated` with outcome
`:running`, and `succeeded(solution)` is `false` for it (§8.14), so a
termination of yours is never read as a burnout. A callback that tests a
condition on a cell needs the state's layout; §8.11 shows how to reach it.

### The conservation recorder

`ConservationRecorder()` is a read-only observer that integrates the live
region's boundary mass outflow and enthalpy rate over each accepted step
(three-point Gauss–Legendre quadrature on the step's interpolant, up to the
time where an event ends the step early) and records the live inventories at
the start and after every step. Pass it as the callback, alone or with others,
and close the ledger afterwards with `conservation_audit`:

```julia
recorder = ConservationRecorder()
solution = solve(problem; use_ale = true, deactivation = true,
                 callback = (recorder, record_dt))
audit = conservation_audit(solution, recorder)
audit.mass.relative, audit.energy.relative      # the closures' relative defects
```

A step a callback of yours ends early is one the recorder cannot account for;
its integrals become `NaN`. §9.8 explains the audit and its fields.

### Progress bar

`progress = true` draws a percent-of-time-span progress bar, advanced after
accepted steps; it is a read-only observer and leaves the steps unchanged. It
is off by default. Leave it off for batch runs and parameter sweeps.

---

## 8.11 Setting up once: `setup` and `integrate`

`solve(problem; kwargs...)` is `validate_problem`, then `setup`, then
`integrate`. The two halves are exported for scripts that integrate one set-up
problem several times, or that need the state layout:

```julia
cp = setup(problem; use_ale = true, radiation_model = SURFACE_ABSORPTION)  # a CoreProblem
sol1 = integrate(cp; saveat = collect(0.0:1.0:300.0))
sol2 = integrate(cp; saveat = collect(0.0:1.0:300.0))   # bit-identical to sol1
```

`setup` takes the keywords that define the problem the integrator runs
(`energy_state`, `reltol`, `abstol`, `threshold`, `use_ale`, `radiation_model`,
`correction`, `lagrangian`, `w_override`, `deactivation`, `band`, `sources`,
and `time_scale`, the time scale of the velocity tolerances that `solve` sets
from `dt`); `integrate` takes those of the integration (`maxiters`, `dt`,
`adaptive`, `dtmin`, `dtmax`, `saveat`, `callback`, `min_thickness`,
`linear_solver`, `integrator_policy`, `located_events`, `progress`, and an
initial state `u0`, by default the set-up one). Every call of `integrate`
starts with every cell live and from the set-up temperatures, so repeated
integrations of one problem are bit-identical. An initial state `u0` of your
own must have consistent algebraic rows; an inadmissible one is refused.

The set-up problem exposes the layout of the state, which a callback testing a
cell's value needs. In the temperature form the energy row of a cell is its
temperature; this stops the integration when the back cell reaches 500 K:

```julia
using OrdinaryDiffEq: ContinuousCallback, terminate!
using Pyrolysis.Internal: energy_index

cp = setup(problem; use_ale = true, radiation_model = SURFACE_ABSORPTION)
k = energy_index(cp.def.layout, 1)                    # cell 1, at the back face
stop = ContinuousCallback((u, t, integrator) -> u[k] - 500.0, terminate!)
sol = integrate(cp; callback = stop)                  # ends where T₁ = 500 K
```

Two further `solve` keywords belong here:

- **`validate`** (default `true`) runs `validate_problem(problem)` first, which
  refuses an inconsistent problem with a list of its errors and warns about
  degenerate gas transport (Chapter 7). `validate = false` skips it; `setup`
  and `integrate` never validate.
- **`sources`** is the external-source hook for verification by manufactured
  solutions: `VolumetricSources(mass, energy)` (from `Pyrolysis.Internal`),
  with `mass(z, t)` returning a tuple of volumetric mass sources per component
  [kg/(m³·s)] and `energy(z, t)` a temperature-rate source [K/s], either
  `nothing`, evaluated at the cell centres. The default `NoSources()` leaves
  the physics alone.

Chapter 11 describes `setup`, `integrate` and the building blocks
`Pyrolysis.Internal` exposes (the residual, the Jacobian cache, the layout).

---

## 8.12 Output times, saving, and MLR for variable cross-section

By default `solve` saves the solution at `problem.output_times` (set when you
built the problem; generated over the time span when you gave none). Override
with `saveat`:

```julia
solve(problem; saveat = collect(0.0:1.0:600.0))   # every second over a 600 s run
solve(problem; saveat = 5.0)                      # every 5 s from the start
```

`saveat` accepts a vector of times or a spacing. The state a solve ends on is
always saved, whether or not `saveat` names its time, so a terminated solve's
last saved time is its termination time. A vector `saveat` saves the start
time only if it lists it. `saveat = Float64[]` saves every accepted step.
Each saved state costs one residual evaluation when the solution is
reconstructed, and its fields' memory; for values between saved times, use
`solution(t)` (Chapter 9).

### Mass-loss-rate gauge normalization

The mass-loss rate is reported two ways:

- `solution.mass_loss_rate[i]`: the mass leaving through every boundary face,
  per unit **current** area `A(t)` of the `:top` faces, in kg·m⁻²·s⁻¹ (the back
  face's share is zero for the default impermeable holder).
- `solution.extras.MLR_gauge[i]`: the same mass per unit **initial** area
  `A₀`, that is `mass_loss_rate × A(t)/A₀`. This is the cone-calorimeter
  convention, since a cone gauge measures mass loss per fixed sample
  footprint.

Without a lateral-shrinkage law the two coincide, because `A(t) = A₀`. Use
`MLR_gauge` when comparing to cone-calorimeter data on a shrinking or swelling
sample. The total mass in `solution.extras.total_mass` is likewise per unit
initial area, and `solution.extras.cross_section_area_history` reports `A(t)`.
See Technical Reference §10 for the gauge convention.

---

## 8.13 Performance tips

A short checklist for fast, reliable solves:

1. **Keep the default policy and linear solver.** `Rodas5P` with
   `StructuredKLU` is what the precompile workloads compile, so it starts
   fastest, and it is what the verification suite exercises (§8.2, §8.4).
2. **Right-size the mesh.** Smaller cells force smaller diffusion-limited
   steps, so cost grows faster than the cell count. Use a mesh graded toward
   the surface (`create_stretched_mesh(L, n, 1.1)`, finest cells at the
   exposed top by default) instead of a uniformly fine mesh; see Chapter 6 and
   the `convergence_douglas_fir` example.
3. **Loosen `reltol` for sweeps** and exploratory runs, and confirm the trend at
   the default tolerance before you rely on it.
4. **Save only what you need.** Fewer saved times mean fewer reconstructions
   (§8.12); a `ConservationRecorder` gives an accurate ledger without dense
   saving (§9.8).
5. **Disable the progress bar** (`progress = false`, the default) and keep the
   `extract` closures small in parameter sweeps (Chapter 12).

Chapter 15 covers slow and failing solves in more depth.

---

## 8.14 Return codes and outcomes

`solve` always returns a `PyrolysisSolution`. Two fields say how it ended:
`solution.retcode`, the integrator's return code, and `solution.outcome`, what
the deactivation or the thickness stop made of it. `succeeded` and `completed`
combine them:

```julia
solution = solve(problem; use_ale = true, deactivation = true)
solution.retcode        # e.g. :Terminated
solution.outcome        # e.g. :burned_out
succeeded(solution)     # true: :Success, or a burnout / min_thickness termination
completed(solution)     # false: true only for :Success
```

Common return codes:

| `retcode` | Meaning | Typical response |
|---|---|---|
| `:Success` | Integrated to `tspan[2]` | — |
| `:Terminated` | Stopped by the deactivation (`:burned_out`, `:interior`), the thickness stop (`:min_thickness`) or a callback of yours (`:running`) | Read `solution.outcome`; widen `min_thickness` if the stop is premature |
| `:Failure` | The domain test gave up (§8.3), a fixed step left the physical domain, or a deactivation event left an inadmissible or inconsistent state (§8.9); `solve` warns with the reason | Read the warning; see Chapter 15 |
| `:MaxIters` | Hit `maxiters` before the end | Raise `maxiters` or loosen the tolerances |
| `:DtLessThanMin` | The step fell below `dtmin` | Read the warning that names the cell; check the problem before lowering `dtmin` |

`succeeded(solution)` is `true` for `:Success` and for `:Terminated` with the
outcome `:burned_out` or `:min_thickness`; it is `false` for `:interior`, for
`:Failure`, and for a termination by your own callback. `completed(solution)`
is the strict test, `:Success` only. A failed solve ends on the last state it
accepted, which is saved.

`solution.statistics` holds the solve's counters: `accepted` and `rejected`
steps, `rhs_evaluations`, `jacobians`, `linear_solves`, `nonlinear_failures`,
`domain_rejections`, `largest_constraint_residual` and `invocations` (the
deactivation callback's invocations that deactivated cells).
`solution.solve_time` is the wall time in seconds. For `:Failure`,
`:MaxIters` or `:DtLessThanMin`, see the troubleshooting guidance in
Chapter 15.

---

## 8.15 A complete annotated `solve` call

Putting the pieces together, here is a fully specified solve for an opaque,
receding cone-calorimeter sample whose surface cells deactivate as they burn
out, with the conservation ledger recorded:

```julia
using Pyrolysis

recorder = ConservationRecorder()
solution = solve(
    problem;
    # accuracy
    abstol    = 1e-8,                     # temperature tolerance [K]
    reltol    = 1e-6,
    threshold = 1e-2,                     # ε_m; also scales the mass tolerances
    dtmin     = 1e-14,                    # allow tiny steps through the reaction front
    # integrator and linear solver (the defaults, spelled out)
    integrator_policy = default_policy(StructuredKLU()),
    # physics
    use_ale         = true,               # moving mesh: the surface recedes
    radiation_model = SURFACE_ABSORPTION, # opaque: absorb the incident flux at the surface
    deactivation    = true,               # burned-out surface cells leave the live mesh
    min_thickness   = 0.01e-3,            # stop when the live thickness is < 0.01 mm
    # observers and output
    callback = recorder,                  # conservation ledger over the accepted steps
    saveat   = collect(0.0:1.0:600.0),    # save once per second
    progress = true,
)

@assert succeeded(solution)
audit = conservation_audit(solution, recorder)
```

From here, Chapter 9 covers the `PyrolysisSolution` object: profiles, surface
and back-face temperatures, the mass-loss rate, the conservation audit and
plotting.
