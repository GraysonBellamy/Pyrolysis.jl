# 8. Running the Solver

Once you have a `PyrolysisProblem` (Chapters 5–7), the single function that
integrates it forward in time is `solve`. This chapter is the complete option
reference for `solve(problem; ...)`: how to pick the time integrator and its
tolerances, how to choose and configure the structured Jacobian and the linear
solver behind it, how to switch on optional physics (ALE mesh motion, radiation,
the gas-generation energy sink), how callbacks are wired in, how warm-starting
and output saving work, and how to read mass-loss-rate output for variable
cross-sections. Every keyword shown here is a real `solve` keyword — `solve`
rejects unknown keywords with an error, so there is no room for typos. The
physics behind these options is derived in the Technical Reference; this chapter
tells you which knob to turn and what each one does. For the underlying temporal
integration and Jacobian theory, see Technical Reference §15.

---

## 8.1 The `solve` entry point

`solve` takes the problem as a positional argument and everything else as
keywords:

```julia
using Pyrolysis

solution = solve(problem)   # all defaults: KenCarp4 + structured Jacobian, no radiation
```

That call already does a great deal: it validates the problem, builds the
immutable problem definition and a mutable workspace, packs the initial state
vector, constructs the default structured Jacobian, assembles the callback set,
runs an implicit ODE integrator, and extracts a `PyrolysisSolution`. It leaves
`problem` unchanged: the run works on a copy with its own mesh, so a problem can
be solved again (or shared by parameter-sweep points) and always starts from the
geometry you built. The copy is kept as `solution.problem`; its mesh holds the
final geometry, live cells and boundary states. The remaining sections explain
how to override each part.

The full keyword signature, with defaults, is:

| Keyword | Default | What it controls |
|---|---|---|
| `integrator` | `nothing` → `KenCarp4` | ODE integrator (constructor or built instance); §8.2 |
| `linear_solver` | `nothing` → `KLUFactorization()`, or `UMFPACKFactorization()` for the exact ALE Jacobian | Sparse factorization paired with a constructor `integrator`; §8.4 |
| `abstol` | `1e-8` | Absolute error tolerance; §8.3 |
| `reltol` | `1e-6` | Relative error tolerance; §8.3 |
| `maxiters` | `1000000` | Maximum integrator steps |
| `dtmin` | `1e-12` | Minimum time step [s] |
| `dtmax` | `Inf` | Maximum time step [s] |
| `jacobian` | `nothing` → `Structured(DirectSolve())` | Pluggable Jacobian backend; §8.5 |
| `verify_jacobian` | `nothing` | Compare active Jacobian to a dense AD reference; §8.6 |
| `use_ale` | `false` | Enable ALE moving-mesh integration; §8.7 |
| `radiation_model` | `NO_RADIATION` | `NO_RADIATION` / `SURFACE_ABSORPTION` / `BEER_LAMBERT`; §8.8 |
| `energy_form` | `:standard` | `:standard` or `:with_generation_sink`; §8.9 |
| `use_transpiration_bc` | `false` | Surface transpiration enthalpy term — must stay `false`; §8.9 |
| `min_thickness` | `1e-6` | Terminate ALE solve when sample thins below this [m]; §8.7, §8.10 |
| `handle_depletion` | `false` | Merge depleted surface cells (ALE only); §8.10 |
| `depletion_threshold` | `0.05` | Cell-thickness fraction that triggers a merge; §8.10 |
| `min_cells` | `2` | Stop merging once this many cells remain; §8.10 |
| `diagnostics` | `false` | Snapshot the conservation ledger during the solve; §8.11 |
| `diagnostics_stride` | `10` | Record a diagnostics snapshot every N steps; §8.11 |
| `callback` | `nothing` | One extra user callback; §8.11 |
| `warm_start` | `nothing` | Reuse a prior solution's final state as the IC; §8.12 |
| `saveat` | `problem.output_times` | Times at which to store output; §8.13 |
| `sensealg` | `nothing` | Forwarded to the integrator (used by the sensitivity extension) |
| `progress` | `true` | Show a progress bar |
| `validate` | `true` | Run `validate_problem` before solving |

`solve` always returns a `PyrolysisSolution`; see Chapter 9 for its fields and
post-processing accessors.

---

## 8.2 Choosing the time integrator

Pyrolysis problems are **stiff**: Arrhenius reaction rates change by orders of
magnitude across a narrow temperature band, and conduction over fine cells
imposes a small explicit stability limit. The solver therefore defaults to an
implicit, stiffly stable method. When `integrator = nothing` (the default),
`solve` builds

```julia
KenCarp4(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization())
```

with `UMFPACKFactorization()` in place of KLU for the exact ALE Jacobian
(§8.4).

`KenCarp4` is a 4th-order ESDIRK method — a good general-purpose choice for the
moderately stiff systems that pyrolysis produces. The `autodiff = AutoFiniteDiff()`
setting tells the integrator to finite-difference its *own* Jacobian when no
pluggable Jacobian is supplied; when you let the default structured Jacobian take
over (§8.5), that pluggable Jacobian is used instead.

You can pass the integrator in three forms:

**1. A constructor** (a type or function). `solve` instantiates it for you and
attaches the chosen linear solver:

```julia
using OrdinaryDiffEq            # solver types; `using Pyrolysis` already loads OrdinaryDiffEq

solve(problem; integrator = KenCarp4)              # rebuilt with default linsolve
solve(problem; integrator = TRBDF2)                # another SDIRK alternative
solve(problem; integrator = FBDF)                  # multistep BDF, good for very stiff TGA-like runs
```

**2. A constructor with a paired `linear_solver`.** Use this when you want a
specific sparse factorization but still want `solve` to wire it together:

```julia
solve(problem; integrator = KenCarp4,
               linear_solver = SparspakFactorization())
```

**3. A fully constructed algorithm instance.** Build it yourself, including the
linear solver, and pass it as `integrator`. In this case do **not** also pass
`linear_solver` — `solve` raises an error if you pair `linear_solver` with an
already-built instance:

```julia
solve(problem;
      integrator = KenCarp4(linsolve = SparspakFactorization()))
```

This third form is what the shipped Douglas-fir example uses
(`examples/07_douglas_fir_cone.jl`).

Recommended starting points:

| Situation | Integrator |
|---|---|
| Default cone / TGA / general use | `KenCarp4` (the default) |
| Very stiff, slow to converge | `FBDF` or `TRBDF2` |
| Tight accuracy near a sharp reaction front | `KenCarp4` with tighter `reltol` |

All integrators come from OrdinaryDiffEq.jl; any stiff solver from that package
works. See Technical Reference §15 for the implicit-RK formulation and the
`W = I/(γΔt) − J` matrix the integrator factors each step.

---

## 8.3 Tolerances and step-size limits

The integrator adapts its time step to keep the local error estimate within
`abstol` (absolute) and `reltol` (relative). The defaults are:

```julia
solve(problem; abstol = 1e-8, reltol = 1e-6)
```

Guidance:

- **Loosen** (`reltol = 1e-4`, `abstol = 1e-6`) for quick exploratory runs or
  parameter sweeps where you only need trends.
- **Tighten** (`reltol = 1e-7`, `abstol = 1e-9`) for convergence studies,
  Jacobian verification, or when comparing against reference codes. The
  Douglas-fir example uses `abstol = 1e-9, reltol = 1e-7`.
- Tolerances mix temperature (≈ 10²–10³ K) and concentration (≈ 0–10³ kg/m³) in
  one state vector, so prefer driving accuracy through `reltol` and keep
  `abstol` comfortably below the smallest concentration you care about.

Step-size limits:

```julia
solve(problem;
      dtmin = 1e-14,   # allow tiny steps through a rapid reaction front
      dtmax = Inf,     # no upper cap (default)
      maxiters = 2_000_000)
```

`dtmin` matters in ALE/charring runs: as cells thin, the diffusive stability
scale (`dt ∝ Δz²`) shrinks, and the integrator may need very small steps to push
through the front. If a solve stops with `retcode == :DtLessThanMin`, lower
`dtmin`. If it stops with `:MaxIters`, raise `maxiters` (or loosen tolerances).

`solve` also caps the step after Newton failures. When an accepted step needed
Newton retries, the next 20 steps are held at or below 1.25 times that step,
so the controller cannot regrow the step straight back into the range where
the Newton failed. A run without Newton failures is unaffected. Compare
`sol.raw_solution.stats.nnonlinconvfail` with the number of accepted steps
(`naccept`) to see how often the Newton fails.

---

## 8.4 The linear solver

Each implicit step solves a sparse linear system in the `W` matrix. With
`linear_solver = nothing`, `solve` picks a sparse direct factorization from the
Jacobian's sparsity pattern:

- **`UMFPACKFactorization()`** for the exact ALE Jacobian
  (`Structured(strategy = DirectSolve())` under ALE). Its geometry block couples
  every cell to every node below it, and UMFPACK's supernodal factorization is
  the fastest on that dense block.
- **`KLUFactorization()`** for every other sparse Jacobian: every
  `ApproxSolve` Jacobian (including the ALE default) and every Eulerian one. KLU refactors the fixed
  pattern in place without allocating, and on these patterns it runs several
  times faster end to end than Sparspak.

The choice depends on the mode and the Jacobian strategy, not on the mesh size.

You select the linear solver either by pairing it with a *constructor*
`integrator`:

```julia
using LinearSolve: SparspakFactorization, UMFPACKFactorization

solve(problem; integrator = KenCarp4,
               linear_solver = UMFPACKFactorization())
solve(problem; integrator = KenCarp4,
               linear_solver = SparspakFactorization())   # pure-Julia sparse LU
```

or by baking it into a fully built integrator:

```julia
solve(problem; integrator = KenCarp4(linsolve = KLUFactorization()))
```

KLU, UMFPACK and Sparspak are interchangeable sparse direct solvers; the
default only picks the one suited to the pattern. Iterative ILU-preconditioned
solvers for large ALE problems are covered in §8.14.

> Note: a dense Jacobian backend such as `DenseAD_ForwardDiff` has no sparse
> Jacobian prototype for a sparse factorization to consume. `solve` then
> defaults to `LUFactorization()` (dense), and replaces an explicitly requested
> sparse factorization with it — with a warning, keeping your integrator type
> and `autodiff` setting.

---

## 8.5 Jacobian backend selection

The single most important performance lever is the Jacobian. By default, `solve`
builds the **structured analytical Jacobian**:

```julia
Structured(strategy = DirectSolve())
```

This is selected automatically when `jacobian = nothing` (the default) for both
Eulerian and ALE problems. It assembles the exact Jacobian from per-operator
analytical and AD contributions (reaction, conduction, gas transport, radiation,
boundary conditions, and the ALE structured couplings) and hands a sparse matrix
to the linear solver. It is far faster than letting the integrator
finite-difference the whole RHS. See Technical Reference §15 for the
operator-based architecture.

You override it through the `jacobian` keyword, passing a `JacobianSpec`. The
spec is built from a *backend* and a *linear-solve strategy*. The constructor
`JacobianSpec(backend, problem; ...)` builds the sparsity patterns and AD caches
once; **do not** construct it inside a hot loop.

### Linear-solve strategies

| Strategy | Behavior | When to use |
|---|---|---|
| `DirectSolve()` | Materialize the full Jacobian (block + densified structured couplings) and factor it. Exact. | Fixed-mesh default; the exact ALE Jacobian when asked for |
| `StructuredSolve()` | Sparse LU on the local block plus Richardson iterations that apply the structured couplings as fast O(N) operators. Exact, converges in 2–5 iterations. | Typical ALE problems |
| `ApproxSolve()` | Drop the structured ALE couplings (Schur-style approximation); solve only the local block. Approximate but much cheaper. | ALE default |
| `MatrixFreeSolve()` | GMRES on the operator action without materializing. | Research / very large problems |

The verifier-only backend `DenseAD_ForwardDiff()` runs ForwardDiff over the
entire residual to produce a dense reference Jacobian. It is correct by
construction but slow — use it only for verification (§8.6), never for
production.

### Examples

Production default (no keyword needed):

```julia
solve(problem)                                  # Structured(DirectSolve())
```

Opt into the exact structured/Richardson solve:

```julia
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem)
solve(problem; jacobian = spec)
```

Opt into the exact ALE Jacobian instead of the reduced default. Because the
spec must match the simulation mode, pass `use_ale = true` to **both** the
`JacobianSpec` constructor and `solve`:

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true)
solve(problem; jacobian = spec, use_ale = true)
```

The exact Jacobian carries the dense mesh-velocity coupling, which the reduced
default drops; both follow the same trajectory, and the exact one costs several
times more per step (Technical Reference §15.10.3).

The mode keywords accepted by `JacobianSpec(backend, problem; ...)` mirror what
`solve` selects: `use_ale`, `radiation_model`, `radiation_intensity`,
`energy_form`, `use_transpiration_bc`, `min_thickness`, `min_density`. They
default to Eulerian / no radiation, so a plain Eulerian problem needs none of
them.

> Depletion handling (`handle_depletion = true`) keeps the state size and
> works with any `jacobian`. See §8.10 and Technical Reference §14.3.

---

## 8.6 Verifying the Jacobian during a solve

To check that your chosen backend agrees with the dense AD reference at runtime,
pass `verify_jacobian` a NamedTuple. At each listed sample time the solver
compares the active Jacobian to a `DenseAD_ForwardDiff()` reference and **warns**
(it does not throw) if the deviation exceeds tolerance:

```julia
solve(problem;
      verify_jacobian = (rtol = 1e-6, atol = 1e-10,
                         sample_times = [0.0, 1.0, 10.0]))
```

Allowed keys are exactly `:rtol`, `:atol`, and `:sample_times`; any other key is
rejected. `rtol` defaults to `1e-6` and `atol` to `1e-10` if omitted. Use this
when developing a custom backend or diagnosing a "Jacobian mismatch" symptom
(slow Newton convergence, unexpected step rejections).

Verification is skipped — with a warning — if no Jacobian spec is active, which
happens under `handle_depletion = true`. For an out-of-solve check at a single
state point, call `verify_jacobian_at` directly (see Technical Reference §15 and
Chapter 11).

---

## 8.7 ALE moving-mesh integration

By default the mesh is **Eulerian** (fixed). For charring, ablating, swelling, or
intumescent materials whose surface recedes or expands, enable the Arbitrary
Lagrangian–Eulerian (ALE) formulation, which adds the node positions `z` to the
state vector and moves the mesh with the material:

```julia
solve(problem; use_ale = true)
```

With `use_ale = true`:

- Node positions evolve as state variables, so the domain thickness
  `L = z[end] − z[1]` changes in time and is reported in
  `solution.extras.thickness_history`.
- A thickness-termination callback stops the integration when the sample thins
  below `min_thickness` (default `1e-6` m), preventing a singular `ρc_p → 0`. Set
  it to your physical stopping criterion, e.g. `min_thickness = 0.5e-3`.
- If the material carries a lateral-shrinkage law, the WithChi cross-section
  state `χ` is integrated too, and the cross-sectional area `A(t)` varies (see
  §8.13 and Technical Reference §10).

A typical opaque charring cone run:

```julia
solve(problem;
      use_ale = true,
      radiation_model = SURFACE_ABSORPTION,
      min_thickness = 0.5e-3,                # stop below 0.5 mm
      integrator = KenCarp4(linsolve = SparspakFactorization()))
```

For Eulerian problems (`use_ale = false`), the mesh never moves, the thickness is
constant, and `solution.z` is left empty. ALE node-motion theory and the
Geometric Conservation Law are in Technical Reference §11.

---

## 8.8 Radiation model selection

The `radiation_model` keyword selects how incident radiant flux (supplied by a
`RadiativeFluxBC` on the top boundary, see Chapter 5) is deposited into the
domain. The three core-path values are:

| Value | Physics | Cost / Jacobian |
|---|---|---|
| `NO_RADIATION` (default) | No volumetric radiative heating | O(1), tridiagonal |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the surface (opaque) | O(1), tridiagonal |
| `BEER_LAMBERT` | Exponential in-depth absorption `q''' = α·I` (semi-transparent) | O(n), non-local upper-triangular Jacobian |

```julia
solve(problem; radiation_model = SURFACE_ABSORPTION)   # opaque char, wood
solve(problem; radiation_model = BEER_LAMBERT)         # semi-transparent polymer
```

Notes and guards (all enforced by `solve`):

- **`SURFACE_ABSORPTION` with absorption coefficients set.** If your material
  components carry non-zero mass absorption `α`, `SURFACE_ABSORPTION` *ignores*
  them and applies the absorbed flux at the surface; `solve` warns. Use
  `BEER_LAMBERT` if you want the in-depth profile.
- **`BEER_LAMBERT` without absorption coefficients.** If no component has an `α`,
  `BEER_LAMBERT` lets all radiation pass through unabsorbed; `solve` warns. Set
  `α` (in m²/kg) on the components.
- **`BEER_LAMBERT` with a bottom `RadiativeFluxBC`.** The in-depth sweep models
  top-incident radiation only, and the boundary dispatch zeroes a
  `RadiativeFluxBC`'s face flux under `BEER_LAMBERT` — back-side heating would
  be silently lost, so `solve` raises an `ArgumentError`. Use
  `SURFACE_ABSORPTION` (face deposit) or a `HeatFluxBC` for back-side radiant
  heating.
- **`P1_QUASI_STEADY` is not available on the core path.** It lives in
  `Pyrolysis.Experimental.RadiationP1`, and `solve` raises an error
  if you pass it. Drive P1 through the Experimental primitives directly.

### Time-varying radiation flux is always tracked

Whenever `radiation_model` is `SURFACE_ABSORPTION` or `BEER_LAMBERT`, `solve`
re-evaluates the incident flux `I(t)` of the top `RadiativeFluxBC` at **every**
RHS call — there is no up-front constancy probe, so a ramped, delayed, or
modulated cone flux is tracked exactly no matter when it starts varying. The
cost is one call of your flux function per RHS evaluation, which is negligible.
Just make the `RadiativeFluxBC` flux a function of `t` for a time-dependent
source:

```julia
top_bc = RadiativeFluxBC(flux = t -> 25e3 + 250.0 * t, absorptivity = 0.9) +
         RadiativeBC(ε = 0.9, T_inf = 300) +
         ConvectiveBC(h = 10, T_inf = 300)
```

Re-radiation losses (`RadiativeBC`) and convective losses (`ConvectiveBC`) are
ordinary thermal boundary conditions and are applied regardless of
`radiation_model`. The radiation physics and the surface energy balance are
derived in Technical Reference §8 and §12.

---

## 8.9 Energy formulation and the transpiration caveat

The bulk energy equation uses **Formulation A** (the ThermaKin/Gpyro/FDS
convention): gas sensible storage is excluded from the effective heat capacity
and the gas enthalpy is carried by the volumetric advection source `S_conv`,
which is always active. The `energy_form` keyword selects whether an additional
gas-generation enthalpy sink is included:

```julia
solve(problem; energy_form = :standard)              # default: S_conv only
solve(problem; energy_form = :with_generation_sink)  # adds −Σ Δh_g (∇·J_g) cell source
```

- `:standard` (default) matches ThermaKin Eq. 17 — `S_conv` alone carries gas
  enthalpy. This is the recommended convention for almost all work.
- `:with_generation_sink` adds the per-cell gas-generation sink term, which
  telescopes with surface transpiration to recover the strict conservative form.
  Use it only if you specifically need the conservative-form energy balance to
  close (see Technical Reference §7 and §16).

Any other symbol is rejected.

### Do not enable `use_transpiration_bc`

The surface transpiration enthalpy term `Σṁ_g·Δh_g` would **double-count** gas
enthalpy because `S_conv` already carries it. `solve` therefore hard-rejects
`use_transpiration_bc = true`:

```julia
solve(problem; use_transpiration_bc = false)   # the only accepted value (default)
```

Leave it at its default `false`. Mass transport across the boundary remains fully
active either way — this flag concerns only the *enthalpy* sink, not the mass
flux. (Parameter sets calibrated against solvers that apply a surface transpiration
term *in addition to* a volumetric advection source embed the double-count; they
predict lower surface temperature and MLR than the `S_conv`-only balance and should
be re-validated when imported.)

---

## 8.10 Depletion handling and termination

When a layer burns out, its cells shrink (ALE) or empty (fixed mesh) toward a
state the implicit integrator cannot resolve — a vanishing cell's equations are
conditioned like `1/Δz²`. Depletion handling **parks** such cells: it merges
each one into a neighbour, conserving mass and condensed enthalpy, and moves
its state slot above the live surface as an inert ghost. The state vector,
mesh and Jacobian keep their size; the boundary conditions act at the live
surface.

```julia
solve(problem;
      use_ale = true,
      handle_depletion = true,    # park depleted cells
      depletion_threshold = 0.05, # park a cell once it is < 5% of its initial thickness
      min_cells = 2,              # keep at least 2 live cells; otherwise terminate
      min_thickness = 0.01e-3)    # terminate when the live thickness < 0.01 mm
```

How it works:

- After every accepted step a callback scans the live cells. Under ALE a cell
  is depleted when its thickness is below `depletion_threshold ×` its **own**
  initial thickness (default `0.05`, range `(0,1)`; the reference is per cell, so
  stretched meshes are judged cell-by-cell). On a fixed mesh (`use_ale = false`)
  the top live cell is depleted when its condensed concentration has fallen
  below the same fraction of its initial value.
- On a fixed mesh only the top live cell parks. When a cell beneath it becomes
  a void first, the material has gasified from within: the solve stops with a
  warning and retcode `:Terminated`, and a regressing material belongs on the
  ALE mesh (§10.2.1).
- The depleted cell is merged into its interior neighbour with an
  **energy-conserving** temperature merge and a volume-weighted concentration
  merge, and its slot is parked; one cell per step, surface first.
- Parking stops once only `min_cells` live cells remain (default `2`); the
  thickness-termination callback ends an ALE solve when the live thickness drops
  below `min_thickness`.
- Ghost slots read as `NaN` in `sol.T`, `sol.ξ` and `sol.z`;
  `sol.extras.live_cells` gives the live count at each output time and
  `sol.extras.park_times` the event times.

`examples/06_cone_calorimeter.jl` uses this configuration
(`handle_depletion = true`, `depletion_threshold = 0.05`).
`examples/07_douglas_fir_cone.jl` instead lets ALE carry the thin cells and
stops purely on `min_thickness` (`handle_depletion = false`).

> Jacobian: depletion handling keeps the pluggable Jacobian. The ALE default,
> `Structured(strategy = ApproxSolve())`, runs through the parks, and
> `DirectSolve` stays available through `jacobian`. Theory in Technical
> Reference §14.3.

---

## 8.11 Callbacks and runtime diagnostics

`solve` assembles a `CallbackSet` internally; you do not manage most of these. In
order, they are: a step-update callback (always present; syncs the integrator
state into the mutable mesh), a thickness-termination callback (ALE), a depletion
callback (`handle_depletion`), the Newton-failure step cap (§8.3; integrators
that solve their stages by Newton), a progress bar (`progress = true`), an
optional diagnostics callback, an optional Jacobian-verification callback, and
finally any user callback you supply.

### Conservation diagnostics

Set `diagnostics = true` to snapshot the conservation ledger during the solve.
Snapshots are recorded every `diagnostics_stride` steps (default `10`) and are
appended to `solution.extras.diagnostics_log`:

```julia
solution = solve(problem; diagnostics = true, diagnostics_stride = 20)
rows = solution.extras.diagnostics_log    # Vector{DiagnosticsRow}
```

Each row captures the matrix sensible energy `E_matrix`, total energy `E_total`,
and the gas-storage gap. This is **snapshot** logging — not a continuous time
history — and is intended for auditing energy/mass balance, not dense
post-processing. `diagnostics_stride` must be ≥ 1. See Chapter 9 for reading
diagnostics and Technical Reference §16 for the ledgers.

### A custom user callback

You may pass one additional callback (any OrdinaryDiffEq callback type) via
`callback`. It runs last in the set:

```julia
using OrdinaryDiffEq: ContinuousCallback, terminate!  # the callback types

stop_at_T = ContinuousCallback(
    (u, t, integ) -> u[1] - 700.0,        # back-face cell (index 1) reaches 700 K
    integ -> terminate!(integ))

solve(problem; callback = stop_at_T)
```

For automated convergence checks against the dense AD reference, use
`verify_jacobian` (§8.6) rather than a hand-rolled callback.

### Progress bar

`progress = true` (default) shows a percent-of-time-span progress bar, advanced
after accepted steps; it leaves the step sequence unchanged. Set
`progress = false` for batch runs and parameter sweeps (the sweep harness does
this automatically).

---

## 8.12 Warm-starting from a prior solution

To start a new solve from the final state of a previous one — useful for
continuation studies and for warm-starting nearby parameter-sweep points — pass
the prior `PyrolysisSolution` as `warm_start`:

```julia
sol1 = solve(problem1)
sol2 = solve(problem2; warm_start = sol1)
```

The prior solution's final raw state vector overwrites the new initial
condition. Compatibility is enforced by a state-vector length check, which
implicitly requires matching component count, mesh size, ALE mode, and χ mode;
mismatches raise a clear error (and the prior solve must have stored its raw
solution). A solve that parked cells cannot seed a warm start: its final state
holds ghost slots, while a new solve starts with every cell live. Note that
`tspan[1]` of the *new* problem still sets the integrator's
start time — warm-starting treats the prior final state as fresh initial data, it
does not continue the clock.

> Caution: warm-starting can introduce subtle artifacts in sensitivity studies
> because the perturbed solve no longer begins from the canonical initial
> condition. Prefer cold starts when computing gradients (Chapter 13).

---

## 8.13 Output times, saving, and MLR for variable cross-section

By default `solve` saves the solution at `problem.output_times` (set when you
built the problem). Override with `saveat`:

```julia
# Save every second over a 600 s run
solve(problem; saveat = collect(0.0:1.0:600.0))
```

`saveat` accepts a vector of times or a step. Fewer save points reduce
interpolation overhead and memory; for the raw dense interpolant use
`interpolate_solution(solution, t)` (Chapter 9).

### Mass-loss-rate gauge normalization

For variable cross-section (WithChi) runs, the mass-loss rate (top-face flux
plus any back-face venting; the bottom share is zero for the default
impermeable holder) is reported two ways in the solution:

- `solution.mass_loss_rate[i]` — the **local** MLR per *current* area `A(t)`,
  in kg·m⁻²·s⁻¹ (both faces; the bottom term is zero for an impermeable
  holder).
- `solution.extras.MLR_gauge[i]` — the **gauge** MLR normalized to the *initial*
  area `A_0`, i.e. `MLR_local × A(t)/A_0`. This is the cone-calorimeter
  convention, since a cone gauge measures mass loss per fixed sample footprint.

For constant cross-section (no lateral-shrinkage law) the two coincide because
`A(t) = A_0`. Use `MLR_gauge` when comparing to cone-calorimeter data on a
shrinking or swelling sample. The total domain mass in
`solution.extras.total_mass` is likewise normalized by the initial area. See
Technical Reference §10 and §16 for the gauge convention.

---

## 8.14 Performance tips and iterative linear solvers

A short checklist for fast, reliable solves:

1. **Use the default structured Jacobian.** It is dramatically faster than the
   integrator's finite-difference Jacobian.
2. **Right-size the mesh.** Resolution costs scale super-linearly (smaller cells
   force smaller diffusion-limited time steps). Use a stretched mesh refined near
   the surface (`create_stretched_mesh(L, n, Val(NC), 1.1; stretch_from = :top)`)
   instead of a uniformly fine mesh; see Chapter 6 and the
   `convergence_douglas_fir` example.
3. **Under ALE, keep the reduced default.** `ApproxSolve` drops the dense ALE
   couplings and the Newton iteration still converges; pass `DirectSolve` only
   to check a configuration against the exact Jacobian (§8.5), or verify with
   `verify_jacobian`.
4. **Loosen tolerances for sweeps.** `reltol = 1e-4` is plenty for trend studies.
5. **Disable the progress bar** (`progress = false`) and use small `extract`
   closures in parameter sweeps (Chapter 12).

### ILU-preconditioned iterative solvers

For very large ALE problems where direct factorization dominates the cost, two
ILU-preconditioned Krylov solvers ship with the package and plug in as the
integrator's linear solver:

```julia
using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization

# ILU-preconditioned GMRES
solve(problem;
      integrator = KenCarp4(linsolve = ILUGMRESFactorization(
          τ = 0.1,            # ILU drop tolerance
          gmres_rtol = 1e-8,  # GMRES relative tolerance
          restart = 30)))     # GMRES restart frequency

# ILU-preconditioned BiCGSTAB (fixed memory, no restart)
solve(problem;
      integrator = KenCarp4(linsolve = ILUBiCGSTABFactorization(
          τ = 0.1, rtol = 1e-8)))
```

`ILUGMRESFactorization` and `ILUBiCGSTABFactorization` live in the `Solver`
submodule — bring them in with
`using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization`
(or `using Pyrolysis.Internal`), or qualify as `Pyrolysis.ILUGMRESFactorization`.
If the ILU factorization becomes singular during a sharp transient, both solvers
fall back to diagonal preconditioning automatically.

> Trade-off: inexact iterative solves can cause more Newton-step rejections on
> stiff Arrhenius transients, and `solve` warns when you select ILU-GMRES.
> Benchmark an iterative configuration end-to-end against the default direct
> factorization, and reserve it for large ALE meshes where the dense structured
> coupling makes direct factorization the bottleneck.

External factorizations are available through extensions: `using KLU` enables
`Pyrolysis.klu_factorization()`, `using Pardiso` enables
`Pyrolysis.pardiso_factorization()`, and `using AlgebraicMultigrid` enables
`Pyrolysis.amg_preconditioner(A)`. These are covered in Chapter 11. The
structured-solve strategies and ILU internals are described in Technical
Reference §15.

---

## 8.15 Reading the return code

`solve` always returns a `PyrolysisSolution`; check `solution.retcode` to confirm
the run finished as intended:

```julia
solution = solve(problem)
println("retcode    = ", solution.retcode)   # e.g. :Success
println("solve time = ", solution.solve_time, " s")
```

Common return codes (from OrdinaryDiffEq):

| `retcode` | Meaning | Typical fix |
|---|---|---|
| `:Success` | Integrated to `tspan[2]` cleanly | — |
| `:Terminated` | A callback stopped it (thickness/depletion/user) | Expected for charring; widen `min_thickness` if premature |
| `:MaxIters` | Hit `maxiters` before the end | Raise `maxiters` or loosen tolerances |
| `:DtLessThanMin` | Step fell below `dtmin` | Lower `dtmin`; check for unphysical states / over-fine mesh |

A `:Terminated` code is normal and expected for ALE charring runs that deplete
the sample. For `:MaxIters` or `:DtLessThanMin`, see the troubleshooting guidance
in Chapter 15.

---

## 8.16 A complete annotated `solve` call

Putting the pieces together, here is a fully specified solve for an opaque
charring cone-calorimeter case with ALE mesh motion (mirroring the shipped
Douglas-fir example):

```julia
using Pyrolysis
using LinearSolve: SparspakFactorization

solution = solve(
    problem;
    # accuracy
    abstol      = 1e-9,
    reltol      = 1e-7,
    dtmin       = 1e-14,                    # allow tiny steps through the reaction front
    # integrator + linear solver (built-instance form)
    integrator  = KenCarp4(linsolve = SparspakFactorization()),
    # physics
    use_ale          = true,               # moving mesh: surface recedes
    radiation_model  = SURFACE_ABSORPTION, # opaque: deposit all flux at the surface
    energy_form      = :standard,          # S_conv-only (ThermaKin) convention
    min_thickness    = 0.5e-3,             # stop when the sample is < 0.5 mm
    handle_depletion = false,              # let ALE carry thin cells (no merging)
    # diagnostics + output
    diagnostics = true,                    # snapshot the energy ledger
    saveat      = collect(0.0:1.0:600.0),  # save once per second
    progress    = true,
)

@assert solution.retcode in (:Success, :Terminated)
```

The default `jacobian = nothing` gives the reduced `ApproxSolve` Jacobian for
this ALE problem; for the exact one, add
`jacobian = Structured(strategy = DirectSolve())`.

From here, Chapter 9 covers the `PyrolysisSolution` object — extracting
temperature and concentration profiles, surface and back-face temperatures,
mass-loss rate, and plotting.
