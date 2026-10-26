# 11. Solver and Jacobian Configuration; Extensions

This chapter is the reference for *how* Pyrolysis.jl solves a problem, as
opposed to *what* problem it solves (Chapters 3–7) or *how to read the
result* (Chapter 9). It covers the Jacobian (colored sparse forward-mode
automatic differentiation, with nothing to configure), the integrator
policies and their restart and flag rules, the linear solve of the
integrator's Newton matrix (`StructuredKLU`), the plotting extension,
checking the Jacobian against its dense reference, driving the integration
yourself through `setup` and `integrate`, and what `Pyrolysis.Internal`
exposes (the residual, the Jacobian cache, the state layout). The
mathematics behind the system the integrator solves, its Jacobian and the
Newton matrix is derived in Technical Reference §15; this chapter is about
the knobs.

Three pieces of the solve are separate, and keeping them apart is the key to
this chapter:

1. **The Jacobian** `J = ∂f/∂u` of the residual. `setup` builds it once per
   problem from the mesh connectivity; there is one Jacobian, no keyword and
   no strategy to choose.
2. **The integrator policy** (`integrator_policy = ...`): the time-stepping
   method together with its rule for the step after an event and its rule
   for the derivative-discontinuity flag after a read-only observer.
3. **The linear solver** (`linear_solver = ...`): how the integrator's
   Newton matrix `W = J − M/(γΔt)` is factored, `StructuredKLU()` by default.

The defaults work for every configuration the package supports, and you can
ignore this chapter until a solve is too slow, fails, or needs to be driven
from your own loop.

---

## 11.1 What the defaults do

`solve(problem)` hands the integrator a mass-matrix system

```
M du/dt = f(u, t)
```

with `M` diagonal: one on the differential rows (the component masses, the
energy state, the node positions, the conversion-history pair), zero on the
algebraic rows (the node velocities of a moving mesh, the face intensities
under `BEER_LAMBERT`, the cross-section area and its lateral rate). On a
fixed mesh without in-depth radiation every row is differential and `M = I`.
With it go (see `src/Integration/solve.jl`):

- the **colored Jacobian** and its sparsity prototype (§11.2);
- the **policy** `default_policy(StructuredKLU())`, which is
  `rodas5p_policy(StructuredKLU())`: `Rodas5P`, a Rosenbrock method that
  evaluates the exact Jacobian at every step, with the heuristic restart
  after an event (§11.3);
- per-row **absolute tolerances**, `abstol = 1e-8` [K] on the temperature
  rows and tolerances derived from each cell's initial condensed mass and
  the deactivation `threshold` on the others, with `reltol = 1e-6` (§11.5);
- the **domain test** `DomainCheck`, which rejects a proposal the error test
  accepts when it is inadmissible or its algebraic rows' constraint residuals
  are above their tolerance (§11.5);
- no initialization of the integrator's own: `setup` makes the algebraic
  rows of the initial state consistent.

```julia
using Pyrolysis

sol = solve(problem)                                            # the defaults
sol = solve(problem; integrator_policy = default_policy(StructuredKLU()))  # the same, written out
```

> Defaults at a glance: `integrator_policy = nothing` →
> `rodas5p_policy(linear_solver)`, `linear_solver = nothing` →
> `StructuredKLU()`, `abstol = 1e-8`, `reltol = 1e-6`, `dtmin = 1e-12`,
> `dtmax = Inf`, `maxiters = 1_000_000`. The Jacobian is always the colored
> one.

The integrator never sees the problem's own types: the residual and the
Jacobian reach it behind problem-independent function wrappers (`ErasedRHS`,
`ErasedJacobian`), and so do the callbacks. The integrator therefore compiles
once per solver configuration (fixed or moving mesh, a mass matrix or not,
deactivation or not), and a second problem of the same configuration reuses
it; the package's precompile workloads cover the common configurations
(Chapter 15).

---

## 11.2 The Jacobian: colored sparse forward-mode differentiation

The Jacobian is computed by forward-mode automatic differentiation of the
residual itself, on a sparsity pattern built from the mesh connectivity, with
the columns grouped by a coloring so that one dual-number evaluation of the
residual yields many columns at once. It is exact to round-off (the test
suite holds it to the dense reference at a relative `1e-12`, §11.7), it
covers every coupling the residual has, including the moving-mesh and
in-depth-radiation couplings, and nothing about it is configurable: `solve`
has no Jacobian keyword, and an unknown keyword raises an `ArgumentError`.

### 11.2.1 The pattern

`setup` builds the pattern once from the connectivity alone
(`jacobian_pattern(cp.def)`): a cell's rows reach the masses, energy and
history of the cells within one face of it (two where a kernel reads a
neighbour's neighbour: the relative fluxes of the equidistribution
correction, whose face states are reconstructed from two cells upwind) and,
on a moving mesh, the positions of their nodes; a node's velocity row reaches
its cell's other node, the cell's masses and energy (through the strain
rate) and, under the correction, what the corrections of the cell's nodes
read; a face-intensity row reaches the neighbouring intensities and the
cells on either side of its face. Every row holds its diagonal. The pattern
is fixed for the solve: deactivation masks rows without changing it
(Chapter 10).

### 11.2.2 The coloring and the evaluation

The columns are colored greedily, keeping the best of four vertex orders,
and the colors are evaluated in batches of at most 24 per dual-number
residual (`MAX_CHUNK`), so a Jacobian costs a few residual evaluations
whatever the mesh size. The preparation runs once, in `setup`; each
evaluation then reuses it at any state and time. For scale: the
PMMA-on-Kaowool deck of Chapter 10 at 15 PMMA cells (46 cells, 3
components) takes 12 colors on a fixed mesh (184 unknowns), and 28 on a
moving mesh with in-depth radiation and the correction (325 unknowns).

Because the derivatives of the material properties enter through dual
numbers, every property function is evaluated at dual-number arguments: a
property given as a function (`CallableProperty`, `StateDependentProperty`)
must be written generically, without restricting its argument to `Float64`
(Chapter 3, Technical Reference §4).

Under a variable cross section the area and lateral-rate rows depend on every
cell; they are colored by their diagonal alone, and their entries are
assembled after the colored pass from column sums of the history rows, which
the colored pass gets exactly.

### 11.2.3 How the integrator uses it

`Rodas5P` evaluates the Jacobian at every step; `FBDF` and `KenCarp4` reuse
it across steps as their Newton iterations allow. Where a method needs the
time derivative of the right-hand side, it is a finite difference
(`AutoFiniteDiff()` in every policy); the Jacobian itself is never a finite
difference. The integrator's Jacobian must be allocated from the prototype
(`jac_prototype`), which the integrator does from what `solve` hands it.

---

## 11.3 Integrator policies

A policy, an `IntegratorPolicy`, is a method together with the two rules that
are part of its contract with the package's callbacks:

| Constructor | Method | Restart after an event | Observer flag |
|---|---|---|---|
| `rodas5p_policy(ls; restart = :heuristic)` (the default) | `Rodas5P`, Rosenbrock, the exact Jacobian every step | `:heuristic` (or `:forced_dt`, `:solver_chosen`) | `:clear` |
| `fbdf_policy(ls)` | `FBDF`, variable-order backward differentiation | `:forced_dt` | `:clear` |
| `kencarp4_policy(ls)` | `KenCarp4`, singly diagonally implicit Runge–Kutta, raw error estimate | `:solver_chosen` | `:keep` |

`ls` is the linear solver (§11.4). All three run with the colored Jacobian and
`AutoFiniteDiff()` for the time derivative; the test suite runs each on a
moving mesh in both energy forms and through deactivation events.

```julia
sol = solve(problem; integrator_policy = fbdf_policy(StructuredKLU()))
sol = solve(problem; integrator_policy = kencarp4_policy(StructuredKLU()))
sol = solve(problem; integrator_policy = rodas5p_policy(StructuredKLU(); restart = :forced_dt))
```

### 11.3.1 The restart rule

An event that modifies the state (a deactivation, Chapter 10) leaves the
integrator at a state its last step size knows nothing about. The policy
names the step it proposes there:

- `:heuristic` proposes the integrator's own initial-step estimate at the
  state the event left (two right-hand-side evaluations). It rejects fewer
  steps after an event than keeping the step, and keeps the mass ledger
  closer than either other rule; it is the default policy's.
- `:forced_dt` proposes `DT_AFTER_EVENT` (1e-6 s), the FBDF policy's, since
  FBDF restarts at order one after an event.
- `:solver_chosen` keeps the integrator's own proposal, the KenCarp4
  policy's.

### 11.3.2 The observer-flag rule

The callback machinery raises the integrator's derivative-discontinuity flag
before every callback. A per-step observer (the progress bar, the
conservation recorder, a `DiscreteCallback` of yours passed as `callback`)
changes nothing, so the policy decides what the flag means after it:

- `:clear` clears it, so a read-only observer is numerically inert. FBDF
  with the flag raised would restart at order one on every step; Rodas5P
  rebuilds its Jacobian every step, so the flag changes nothing for it and
  the policy clears it.
- `:keep` leaves it raised: the integrator re-evaluates the derivative at the
  accepted state, which is the KenCarp4 policy's contract.

`solve` applies the rule inside each discrete callback you pass, after its
`affect!`, so a read-only callback leaves the trajectory unchanged under
every policy (a test asserts identical step sequences with and without one).
A callback must therefore not modify `integrator.u`: observers are read-only
(Chapter 8). A continuous callback's event raises the flag as any event does.

### 11.3.3 The Newton-failure step cap

For the Newton-based policies (FBDF, KenCarp4), after an accepted step whose
nonlinear solves needed retries, the proposed step stays at or below
`FAILURE_MEMORY_GROWTH` (1.25) times that step for the next
`FAILURE_MEMORY_STEPS` (20) accepted steps. A Rosenbrock method has no
nonlinear iteration, so the cap never engages under the default policy.

### 11.3.4 Choosing a policy

Leave the default unless you have a reason. `fbdf_policy` and
`kencarp4_policy` are there for comparison and for problems where a
multistep or an SDIRK method is preferred; on the 46-cell fixed-mesh deck of
§11.2.2 under surface absorption over 300 s, at the default tolerances, the
three policies agree on the final temperatures within 1e-3 K (measured for
this chapter). An `IntegratorPolicy` of your
own is built as `Pyrolysis.Integration.IntegratorPolicy(alg, restart,
observer_flag, restart_dt)`, with any OrdinaryDiffEq method that accepts a
mass matrix and a user Jacobian; only the three constructors above are
tested.

---

## 11.4 The linear solver

The integrator solves its Newton (Rosenbrock stage) systems with
`W = J − M/(γΔt)`; the `linear_solver` keyword names the factorization.

### 11.4.1 `StructuredKLU`, the default

`StructuredKLU(inner = KLUFactorization())` is KLU on the entries of `W` that
have been nonzero at some step of the solve. The Jacobian's prototype holds
every entry the connectivity allows, most of which are zero at any one step,
and KLU factors a stored zero like any other entry. With a mass matrix that
is not a multiple of the identity, OrdinaryDiffEq forms `W` by a sparse
broadcast that stores only the entries that compute to nonzero, so `W`'s own
structure changes whenever an entry passes through zero, and a plain KLU on
that structure would analyse afresh at every such step. `StructuredKLU` keeps
the union of the prototype's entries that have been nonzero in the solve,
with the diagonal, copies each fresh `W` into a matrix of that structure and
refactors it numerically; a nonzero outside the union grows it and triggers
one fresh analysis. The union only grows, so the analyses stop once every
entry the solve reaches has appeared.

### 11.4.2 Other factorizations

Any LinearSolve factorization of a `SparseMatrixCSC` is passed directly; no
extension is involved. `UMFPACKFactorization()` and plain
`KLUFactorization()` are the ones the test suite runs (repeated solves under
each are bit-identical, and a deactivating solve under UMFPACK deactivates
the same cells at the same times to `1e-6`):

```julia
using LinearSolve: UMFPACKFactorization, KLUFactorization

sol = solve(problem; linear_solver = UMFPACKFactorization())
```

> **Pairing rule.** `linear_solver` is used only when `integrator_policy` is
> left at its default. With an explicit policy, pass the linear solver to the
> policy's constructor instead:
> `solve(problem; integrator_policy = rodas5p_policy(UMFPACKFactorization()))`.
> A `linear_solver` given beside an `integrator_policy` is ignored.

### 11.4.3 Linear solve and deactivation

Deactivation keeps the state size and the pattern, and a dead row is an
identity row (Chapter 10), so nothing about the linear solve changes through
burnout: the `StructuredKLU` union keeps its entries, and the factorization
its structure.

---

## 11.5 Tolerances, step limits and the domain test

The tolerances and step limits are independent of the policy and the linear
solver:

| Keyword | Default | Meaning |
|---|---|---|
| `abstol` | `1e-8` | The temperature rows' absolute tolerance [K]; the other rows take theirs from it, `reltol` and `threshold`. |
| `reltol` | `1e-6` | Relative tolerance, every row. |
| `threshold` | `1e-2` | The deactivation threshold, which also scales the mass tolerances. |
| `maxiters` | `1_000_000` | Maximum integration steps. |
| `dt` | `problem.dt_initial` | The first step [s]; the fixed step with `adaptive = false`. Also scales the velocity tolerances. |
| `adaptive` | `true` | Adaptive stepping. |
| `dtmin` | `1e-12` | Minimum step [s]. |
| `dtmax` | `Inf` | Maximum step [s]. |

The absolute tolerances are per row: a mass row `reltol · threshold` times
the cell's initial condensed mass; a temperature row `abstol`; an enthalpy
row (`energy_state = :enthalpy`) `threshold` times the cell's initial
condensed mass, a reference heat capacity and `abstol`, so a cell at the
threshold resolves its temperature as the temperature form does; a
face-intensity row `reltol` times the largest incident intensity; a node
position `reltol · threshold` times the smallest initial size of the cells
around it, and its velocity that over `dt`. `setup` builds the vector
(`cp.abstol`); Chapter 8 discusses the choice.

**The domain test.** `DomainCheck` is the integrator's `isoutofdomain`
test. It reads only the proposals the error test accepts, and rejects one,
cutting the step by the controller's smallest factor, when it is
inadmissible (a live component mass more negative than its tolerance, a live
cell at or below absolute zero, crossed nodes, a non-positive cross section),
when under gas transport a cell's porosity leaves `[0, 1]` by more
than its mass tolerances allow, or when an algebraic row's residual exceeds
`reltol` times its own equation's scale plus its unknown's absolute
tolerance. After more than 50 consecutive constraint rejections
(`MAX_CONSTRAINT_REJECTIONS`, the `max_constraint_rejections` keyword of
`integrate`) the solve ends with `:Failure` and a warning. The largest
scaled constraint residual of the accepted states is reported as
`sol.statistics.largest_constraint_residual` (one or below passes). A
fixed-step solve (`adaptive = false`) runs no domain test; its first accepted
state that is inadmissible ends it with `:Failure`.

Loosening tolerances (`reltol = 1e-4`) speeds up exploratory runs; tighten
them for convergence and validation studies (Technical Reference §18). The
mass is conserved exactly whatever the tolerance, since the masses are the
states; the energy closes to the integrator's tolerance in the temperature
form and to round-off in the enthalpy form (Technical Reference §16).

---

## 11.6 Loading and using extensions

Optional functionality is delivered through Julia's package-extension system.
One extension ships:

| Extension | Trigger (`using ...`) | Attaches | Optional dependency |
|---|---|---|---|
| `PyrolysisPlotsExt` | `Plots` (or `RecipesBase`) | `plot(sol; kind = ...)` recipes | `RecipesBase` |

```julia
using Pyrolysis, Plots

sol = solve(problem)
plot(sol)                      # the surface temperature (kind = :surface_T, the default)
plot(sol; kind = :mlr)         # the mass-loss rate, sol.mass_loss_rate
plot(sol; kind = :heatmap)     # T over time and cell index
```

Any other `kind` raises an `ArgumentError` naming the three. Chapter 9 shows
the plots and how to draw depth profiles by hand.

Linear solvers need no extension: `LinearSolve` is a dependency, so any of
its factorizations is passed as `linear_solver` or to a policy's constructor
(§11.4). Sensitivity analysis needs none either: Chapter 13 takes finite
differences of ordinary `solve` calls, and Technical Reference §17 covers the
theory.

---

## 11.7 Checking the Jacobian against the dense reference

The colored Jacobian has a reference: `dense_jacobian(def, ws, u, t)`
differentiates every column of the residual by dense forward-mode automatic
differentiation, without the pattern or the coloring. Two checks follow from
it, both through `Pyrolysis.Internal`: that the pattern contains every
structurally nonzero entry (`missing_entries`), and that the values agree
(`compare_entries`, per entry with a row-magnitude floor, so a missing
coupling in a small-magnitude row cannot hide behind the matrix maximum).
The test suite runs both on every configuration; you need them only for a
residual you have extended or a state you suspect.

```julia
using Pyrolysis, Pyrolysis.Internal

cp = setup(problem; use_ale = true, deactivation = true, radiation_model = BEER_LAMBERT)
sol = integrate(cp; min_thickness = 1e-3)

# A saved state after the burnout: the mask of its time first.
k = lastindex(sol.t)
u, t = sol.raw_solution.u[k], sol.t[k]
replay_events!(cp.deactivation, cp.def.layout, t)

J = jac_prototype(cp.jacobian)            # the pattern's structure, zero values
jacobian!(J, cp.jacobian, u, t)           # the colored Jacobian
J_ref = dense_jacobian(cp.def, cp.ws, u, t)

isempty(missing_entries(J_ref, jacobian_pattern(cp.def)))   # true
r = compare_entries(J, J_ref; rtol = 1e-12)  # (passed, max_abs_dev, worst_ratio, worst_index)
r.passed                                     # true
```

`compare_entries` reports; it does not throw. `replay_events!` sets the
layout's mask to the one a state saved at `t` had (every event at or before
`t` applied); a state from a solve without deactivation needs no mask. The
dense reference costs one dual-number pass per batch of columns over the
whole state, so keep it to small meshes.

---

## 11.8 Custom loops: `setup` and `integrate`

`solve(problem; kwargs...)` is three steps: `validate_problem(problem)`
(unless `validate = false`), `cp = setup(problem; ...)`, and
`integrate(cp; ...)`. The two halves are exported, so you can set a problem
up once and integrate it several times, or drive the integration yourself.

| `setup` (what the problem is) | `integrate` (how it is run) |
|---|---|
| `energy_state`, `reltol`, `abstol`, `threshold`, `pressure_tolerance`, `use_ale`, `radiation_model`, `w_override`, `correction`, `lagrangian`, `deactivation`, `band`, `sources`, `time_scale` (the first step, `solve`'s `dt`) | `u0`, `maxiters`, `dt`, `adaptive`, `dtmin`, `dtmax`, `saveat`, `callback`, `progress`, `min_thickness`, `linear_solver`, `integrator_policy`, `located_events`, `max_constraint_rejections` |

`setup` returns a `CoreProblem`: the problem, the definition the residual
reads (`cp.def`), its buffers (`cp.ws`), the consistent initial state
(`cp.u0`), the relative tolerance and the absolute-tolerance vector, the
colored Jacobian (`cp.jacobian`), the initial temperatures and masses, and
the deactivation state (`cp.deactivation`, `nothing` without deactivation).
Building the colored Jacobian is the expensive part of `setup`, so a study
that integrates one problem several ways sets it up once:

```julia
cp = setup(problem; use_ale = true, deactivation = true, radiation_model = BEER_LAMBERT)
sol_rodas = integrate(cp; min_thickness = 1e-3)
sol_fbdf  = integrate(cp; min_thickness = 1e-3,
                      integrator_policy = fbdf_policy(StructuredKLU()))
```

Every call starts with every cell live, the layer times before their first
activation and the caloric inversion from the set-up temperatures, so
repeated integrations of one `CoreProblem` are bit-identical. `integrate(cp;
u0 = u)` starts from another state; its algebraic rows must be consistent,
which `consistent_state!(u, cp.def, cp.ws, t)` ensures, and an inadmissible
`u0` is refused with an `ArgumentError`.

**Driving OrdinaryDiffEq yourself.** The residual and the Jacobian are
plain functions of the state, reachable through `Internal`:
`residual!(du, u, cp.def, t, cp.ws)` and `jacobian!(J, cp.jacobian, u, t)`.
A hand-written integration wraps them in the `(du, u, p, t)` shape
OrdinaryDiffEq expects, with the layout's mass matrix:

```julia
using Pyrolysis, Pyrolysis.Internal
using OrdinaryDiffEq: ODEFunction, ODEProblem, Rodas5P, AutoFiniteDiff

cp = setup(problem; radiation_model = BEER_LAMBERT)      # a fixed mesh: no events, no motion
layout = cp.def.layout

rhs!(du, u, p, t) = residual!(du, u, cp.def, t, cp.ws)
jac!(J, u, p, t) = jacobian!(J, cp.jacobian, u, t)
f = ODEFunction(rhs!; mass_matrix = mass_matrix(layout),   # the intensity rows are algebraic
                jac = jac!, jac_prototype = jac_prototype(cp.jacobian))
ode = ODEProblem(f, copy(cp.u0), problem.tspan)
raw = solve(ode, Rodas5P(autodiff = AutoFiniteDiff(), linsolve = StructuredKLU());
            abstol = cp.abstol, reltol = cp.reltol)

T_end = [raw.u[end][energy_index(layout, i)] for i in 1:n_cells(layout)]
```

`solve` is the SciML ecosystem's common `solve` function, so `using
Pyrolysis` beside `using OrdinaryDiffEq` names one function. A loop of your
own has none of what `integrate` adds: no domain test, no deactivation
events, no thickness stop, no observers and no `PyrolysisSolution`; it suits
a fixed mesh without events, or experiments on the residual. The energy row
holds the temperature in the temperature form and the enthalpy inventory
under `energy_state = :enthalpy`; for the temperatures and every other
derived field, evaluate the residual and read `cp.ws.derived` (`T`, `ξ`, `m`,
`V`, `x`, …), or let `integrate` build the solution.

The residual is a pure function of `(u, t)` given the definition: it writes
only into the buffers of `cp.ws`. The definition is immutable apart from the
layout's mask (which events set) and the layer times of a bound correction
(which the layer-time event sets); `integrate` resets both at its start.

---

## 11.9 The `Internal` umbrella

### 11.9.1 What it re-exports

The top-level `Pyrolysis` namespace exports a curated surface (Chapter 16).
`Pyrolysis.Internal` re-exports every symbol the submodules export (Materials,
Properties, BoundaryConditions, Physics, Meshes, States, Kernels, Closures,
Assembly, Deactivation, ColoredJacobian, Problem, Outputs, Integration) in one
namespace; a name two submodules export as different objects is left out and
reached through its submodule (`Pyrolysis.Meshes.X`).

```julia
using Pyrolysis.Internal   # everything the submodules export, in one place
```

| Area | Names |
|---|---|
| The residual (Assembly) | `residual!`, `CoreDefinition`, `CoreWorkspace`, `buffers`, `NoSources`, `VolumetricSources` |
| The Jacobian (ColoredJacobian) | `ColoredJacobianCache`, `jac_prototype`, `jacobian!`, `jacobian_pattern`, `color_count`, `dense_jacobian`, `missing_entries`, `compare_entries` |
| The layout (States) | `StateLayout`, `state_length`, `mass_index`, `energy_index`, `history_index`, `position_index`, `velocity_index`, `intensity_index`, `area_index`, `lateral_rate_index`, `cell_rows`, `node_rows`, `is_differential`, `mass_matrix`, `DerivedState`, `derive!` |
| The integration (Integration) | `consistent_state!`, `inadmissible`, `DomainCheck`, `IntegratorPolicy`, `DT_AFTER_EVENT`, `FAILURE_MEMORY_GROWTH`, `FAILURE_MEMORY_STEPS`, `default_linear_solver`, `CORRECTION_BAND` |
| Deactivation | `CellDeactivation`, `DeactivationEvent`, `deactivate!`, `threshold_margin`, `top_live_cell`, `replay_events!` |
| The mesh motion (Closures) | `BalanceRowClosure`, `LayerCorrection`, `PrescribedCorrection`, `residues`, `vanishes`, `EQUIDISTRIBUTION_TIME`, `ACTIVATION_TIME_FACTOR` |
| Meshes | `line_mesh`, `box_mesh`, `cell_mesh`, `MeshTopology`, `MeshGeometry`, `live_face` |

Symbols reached through `Internal` are not part of the stable public API and
may move or change between minor releases; prefer the top-level exports in
production code. The meshes of two and three dimensions (`box_mesh`,
`cell_mesh`) and their closure (`PotentialClosure`) are residual-level
building blocks: `setup` accepts them, but a full solve on them is not a
supported user path, and gas transport, in-depth radiation, the variable
cross section, the correction and deactivation are refused there.

### 11.9.2 Reading a set-up problem

A `CoreProblem`'s definition `cp.def` holds the state layout (`layout`), the
state setup (`setup`: the mesh topology with its pinned and Lagrangian
nodes, the material, its caloric laws, the per-cell reference data, the
caloric extension and the gas transport `setup` chose from the material),
the kernels in assembler order (`kernels`), the mesh closure (`closure`,
`nothing` on a fixed mesh; its `correction` is the bound correction of
Chapter 10) and the external sources (`sources`).

The layout places every unknown. On a line the state is a sequence of node
and cell blocks, node 1, cell 1, node 2, …, cell N, node N + 1, so the
Jacobian is banded: a node block holds the node's position and velocity on a
moving mesh and the intensity of the face at that node under `BEER_LAMBERT`;
a cell block holds the component masses, the energy state and, under a
variable cross section, the conversion-history pair; the cross-section area
and its lateral rate close the vector. Read a raw state only through the
index functions:

```julia
layout = cp.def.layout
u = cp.u0
state_length(layout)                 # total unknowns
u[mass_index(layout, i, j)]          # mass of component j in cell i [kg]
u[energy_index(layout, i)]           # T [K] (or H [J] under energy_state = :enthalpy)
u[position_index(layout, n)]         # node n's position [m] (moving mesh)
u[velocity_index(layout, n)]         # node n's velocity [m/s] (moving mesh)
is_differential(layout)              # true on the rows with a time derivative
layout.mask                          # the live rows
```

Which configuration a problem gets follows from the problem and the
keywords, never from a mode argument: `use_ale` moves the mesh,
`radiation_model` adds the intensity rows, a permeable condensed component
selects lumped gas transport gas transport, a `lateral_shrinkage_law` adds the
cross-section rows, and `energy_state` picks the energy unknown (Chapter 7).

---

## 11.10 Quick configuration recipes

```julia
using Pyrolysis
using LinearSolve: UMFPACKFactorization

# (a) The default: Rodas5P, the colored Jacobian, StructuredKLU.
sol = solve(problem)

# (b) Another integrator policy.
sol = solve(problem; integrator_policy = fbdf_policy(StructuredKLU()))
sol = solve(problem; integrator_policy = kencarp4_policy(StructuredKLU()))

# (c) Rodas5P with a forced restart step after each event.
sol = solve(problem; integrator_policy = rodas5p_policy(StructuredKLU(); restart = :forced_dt))

# (d) Another LinearSolve factorization, with the default policy or inside a policy.
sol = solve(problem; linear_solver = UMFPACKFactorization())
sol = solve(problem; integrator_policy = rodas5p_policy(UMFPACKFactorization()))

# (e) Burnout on a moving mesh (Chapter 10); nothing about the solver changes.
sol = solve(problem; use_ale = true, deactivation = true, radiation_model = BEER_LAMBERT)

# (f) Set up once, integrate twice.
cp = setup(problem; use_ale = true, deactivation = true, radiation_model = BEER_LAMBERT)
sol1 = integrate(cp)
sol2 = integrate(cp; integrator_policy = fbdf_policy(StructuredKLU()))
```

---

## 11.11 Cross-references

- The mass-matrix system, the policies and their rules, the domain test,
  `StructuredKLU`, the colored Jacobian and its pattern, and the dense
  reference: **Technical Reference §15**.
- Deactivation, the correction and their defaults: **Chapter 10** and
  **Technical Reference §14** and **§11**.
- Building a problem and the configurations `setup` derives from it:
  **Chapter 7**.
- The `solve` keyword reference, callbacks and return codes: **Chapter 8**;
  the solution and post-processing: **Chapter 9**.
- Sensitivity workflows (finite differences of whole solves): **Chapter 13**
  and **Technical Reference §17**.
- Failure modes and performance: **Chapter 15**; the complete API and
  default-value tables: **Chapter 16**.
