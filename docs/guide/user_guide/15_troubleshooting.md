# 15. Troubleshooting and Performance

This chapter is a field guide for when a Pyrolysis.jl simulation fails, runs
slowly, or produces results you do not trust. It is organised by symptom: each
section names a class of failure, explains the underlying cause, and gives the
concrete `solve` options, diagnostic calls, or model changes that fix it. The
emphasis throughout is on *root-causing* — most of the tools below let you
localise a problem to a single subsystem (kinetics, conduction, gas transport,
radiation, ALE, or the linear solver) before you start changing inputs. For the
physics and numerics behind the knobs mentioned here, follow the "see Technical
Reference §N" pointers; for the full option list see the *Running the Solver*
chapter (UG-8) and the *API and Parameter Reference* (UG-16).

Coordinate and sign conventions follow the unified nomenclature: `z = 0` is the
substrate (bottom, boundary id 2, tag `:bottom`), `z = L` is the exposed surface
(top, id 1, tag `:top`); positive flux is transport in `+z`; heat of reaction
`h > 0` is endothermic with `Q_rxn = −Σ h_r r_r`.

---

## 15.1 First steps: read the return code and run the diagnostics

Before changing anything, look at two things: the solver's return code and the
diagnostics log.

### 15.1.1 The return code

`solve` stores the integrator's exit status on the solution:

```julia
solution = solve(problem; progress = false)
solution.retcode      # Symbol mirrored from raw_solution.retcode
```

The common values and what they mean:

| `retcode` | Meaning | Typical cause |
|---|---|---|
| `:Success` | Reached `tspan[2]` | — |
| `:Terminated` | A callback stopped the run | Sample fully depleted (ALE thickness termination), an interior void on a fixed mesh (depletion handling), or a user callback called `terminate!` |
| `:MaxIters` | Hit `maxiters` before `t_end` | Stiff transient; step size collapsed; under-resolved front |
| `:DtLessThanMin` | Step fell below `dtmin` | Discontinuity, NaN in the RHS, or genuine stiffness wall |
| `:Unstable` / `:DtNaN` | NaN/Inf in the state | Unphysical state (negative concentration, `ρc_p → 0`), bad property table extrapolation |

A `:Terminated` run is *not* an error in ALE mode — it usually means the sample
burned away. Check `solution.extras.thickness_history[end]` against
`min_thickness` (default `1e-6` m).
On a fixed mesh with `handle_depletion = true`, a `:Terminated` run with the
warning "…the material has gasified from within…" stopped at the first
interior void (§15.4).

### 15.1.2 Turn on the diagnostics log

`diagnostics = true` installs a callback that records the sensible energy every
`diagnostics_stride` accepted steps (default 10) and stores typed rows in
`solution.extras.diagnostics_log`:

```julia
solution = solve(problem;
    diagnostics        = true,
    diagnostics_stride = 5,     # snapshot every 5 accepted steps
    progress           = false,
)

for row in solution.extras.diagnostics_log
    println("t=$(row.t)  E=$(row.E_sensible)")
end
```

Each `DiagnosticsRow` carries `t`, `E_sensible`, `energy_residual`,
`energy_relative`, and `mass_max_relative_error`; in the rows `solve` records,
only `t` and `E_sensible` are filled and the three cumulative-ledger fields are
`NaN` (§9.8). The sensible-energy history shows where the stored energy changes
abruptly — a park event, a stiff transient — and can be set against the
boundary-flux extras (`surface_heat_flux_*`, §9.3). The cumulative energy and
mass ledgers (`EnergyConservationTracker`, `MassConservationTracker`, driven
step by step through `record_substep!`) are not part of `solve`; see
Technical Reference §16.

`diagnostics` is **snapshot-based**, not continuous — it gives you a coarse time
history, not a dense trace. The overhead is one extra energy sweep every
`diagnostics_stride` steps.

---

## 15.2 Stiff or slow solves

Arrhenius kinetics make the system stiff during ignition and the early reaction
transient; conduction adds a parabolic time scale. If a solve is slow, hits
`:MaxIters`, or grinds to `:DtLessThanMin`, work through the following.

### 15.2.1 Confirm the default integrator is appropriate

The default integrator is `KenCarp4` (an SDIRK method) with a sparse direct
linear solve; `solve` always supplies the structured analytical Jacobian
(§15.2.3), so the `autodiff = AutoFiniteDiff()` setting never builds `J`. When
both `integrator` and `linear_solver` are left at `nothing`, `solve` builds:

```julia
# Internal default (UMFPACKFactorization() for the exact ALE Jacobian, §8.4):
# KenCarp4(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization())
```

`KenCarp4` is a good default for moderately stiff pyrolysis. For *very* stiff
problems (high pre-exponential `A_i`, sharp ignition) a BDF method often takes
larger steps. Integrators live in `OrdinaryDiffEq` and its sublibraries, which
are **not** re-exported by `using Pyrolysis`; load them yourself (`FBDF` and
`Rodas5P` come with `using OrdinaryDiffEq`; `KenCarp4` and `TRBDF2` need
`using OrdinaryDiffEqSDIRK`) and pass a constructor or a fully-built algorithm
to the `integrator` kwarg:

```julia
using OrdinaryDiffEq            # for FBDF, Rodas5P, etc.

# (a) pass a constructor — solve() injects autodiff + linsolve for you:
solution = solve(problem; integrator = FBDF)

# (b) pass a constructor and a linear solver together:
solution = solve(problem;
    integrator    = FBDF,
    linear_solver = SparspakFactorization(),   # using LinearSolve: SparspakFactorization
)

# (c) pass a fully-constructed algorithm — then DO NOT also pass linear_solver:
solution = solve(problem; integrator = FBDF(linsolve = SparspakFactorization()))
```

The pairing rule is strict: `linear_solver` may only be combined with
`integrator` when `integrator` is a **constructor** (a `Type`/`Function`). If
you hand `solve` an already-constructed algorithm *and* a `linear_solver`, it
raises an error telling you to fold the `linsolve` into the algorithm yourself.

### 15.2.2 Adjust tolerances

The error tolerances are `abstol = 1e-8` and `reltol = 1e-6` by default. They
trade accuracy for speed directly:

```julia
# Looser — faster, for exploratory runs:
solution = solve(problem; abstol = 1e-6, reltol = 1e-4)

# Tighter — for convergence/validation studies:
solution = solve(problem; abstol = 1e-10, reltol = 1e-8)
```

If you see `:DtLessThanMin`, first try *loosening* tolerances; a too-tight
`reltol` against a stiff transient drives the step size into the floor. The step
floor itself is `dtmin = 1e-12` s and the ceiling `dtmax = Inf`; raising
`dtmin` is rarely the right fix — it usually just converts a slow run into a
failed one. `maxiters` defaults to `1_000_000`; a `:MaxIters` exit at that count
almost always means the step size collapsed, so chase the stiffness, not the
iteration cap.

### 15.2.3 Use the structured Jacobian (it is already the default)

When `jacobian = nothing` (the default), `solve` builds a
`Structured(strategy = DirectSolve())` analytical Jacobian — or
`Structured(strategy = ApproxSolve())` under ALE —
and hands a sparse prototype to the integrator (see Technical Reference §15).
You normally do not need to touch it. The cases where you *do* are covered in
§15.3 (verification), §15.4 (depletion), and §15.7 (large ALE).

### 15.2.4 Reduce mesh size or use warm starts

Cost scales with the number of cells (the state dimension is `n·(1+NC)` plus the
ALE `z` block and the WithChi `χ` block when present). If a 400-cell run is slow,
check whether 100–150 cells already give a converged surface temperature and MLR
(run a quick convergence check — see Technical Reference §18.5 and the
`convergence_douglas_fir` example). For a *sequence* of related solves
(continuation, sweeps), reuse the prior final state with `warm_start`:

```julia
sol1 = solve(problem1; progress = false)
sol2 = solve(problem2; warm_start = sol1, progress = false)
```

The state-vector length, component count `NC`, ALE mode, and χ mode of the two
problems must match (enforced by a length check), and the prior solve must not have
parked cells (its final state holds ghost slots). `tspan[1]` still sets the
integrator's start time, so the prior final state is treated as fresh initial
data.

---

## 15.3 Jacobian mismatch warnings

If the integrator behaves erratically — many rejected steps, repeated Newton
failures, a sudden step-size collapse — suspect a Jacobian that disagrees with
the residual.

### 15.3.1 Compare against the dense reference at sample times

The `verify_jacobian` kwarg installs a callback that, at the listed sample
times, compares the active backend against a `DenseAD_ForwardDiff` reference
(full-residual ForwardDiff — correct by construction, slow). It **warns**, it
does not throw, so a single transient outlier will not abort a good run:

```julia
solution = solve(problem;
    verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 1.0, 10.0]),
    progress        = false,
)
```

The NamedTuple keys are restricted to `:rtol`, `:atol`, and `:sample_times`
(unknown keys raise an error). A warning reports the backend type, the time, the
max absolute and relative deviations, and the worst index. The callback fires
when the integrator crosses any sample time within half the current step.

### 15.3.2 Verify at a single state point

For a focused check (e.g., debugging an operator), call `verify_jacobian_at`
directly. You need the `ProblemDef`/`Workspace` pair and a state, which
`initialize_state` builds from the problem (with its initial conditions
applied):

```julia
def, ws, u = initialize_state(problem; use_ale = false, radiation_model = NO_RADIATION)

spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
report = verify_jacobian_at(spec, def, ws, u, 0.0;
    reference = DenseAD_ForwardDiff(),   # default reference backend
    rtol = 1e-6, atol = 1e-10,
)

report.passed          # Bool
report.max_abs_dev     # max absolute deviation
report.max_rel_dev     # max relative deviation
report.worst_index     # (row, col) of the worst entry
report.mode            # :strict (exact backends) or :masked (approximate)
```

Pass the same `use_ale` / `radiation_model` keywords the solve uses, so the
spec has the solve's mode (§11.2.2).

### 15.3.3 Support-pattern (sparsity) issues

The comparison `mode` matters:

- **`:strict`** — for the exact backend (`DirectSolve`). Both
  in-support disagreements *and* reference entries the backend's sparsity
  *omits* count as failures. A strict failure on an outside-support entry means
  the sparsity pattern is too small — the backend is structurally missing a
  coupling.
- **`:masked`** — for approximate backends (`ApproxSolve`, where
  `is_approximate(backend)` is `true`). Only in-support disagreements are
  asserted; omitted entries are informational. A masked pass with large
  outside-support reference magnitude is *expected* for `ApproxSolve` (it
  deliberately drops the structured ALE couplings); a masked run that struggles
  to converge is the signal to switch back to an exact strategy (§15.7).

If `verify_jacobian_at` reports a strict, outside-support failure on a
production backend, that is a package bug — report it with the `worst_index` and
the problem configuration. The structured operators are themselves checked
against the dense reference in the test suite (Technical Reference §15.13,
§18.3).

---

## 15.4 Depletion runs

Depletion handling (`handle_depletion = true`) parks a cell when it thins below
`depletion_threshold` of its initial thickness (ALE) or when the top live cell's
condensed concentration falls below that fraction of its initial value (fixed
mesh); the cell is merged into its neighbour and its state slot moved above the
live surface (Technical Reference §14.3). The state keeps its size, so the
structured Jacobian and the sparse direct solver run unchanged through burnout.

```julia
solution = solve(problem;
    use_ale             = true,
    handle_depletion    = true,
    depletion_threshold = 0.05,   # park below 5% of initial thickness
    min_cells           = 2,      # stop parking at 2 live cells
    min_thickness       = 1e-6,   # terminate when the live thickness < this
    progress            = false,
)
```

A fixed-mesh depletion run that stops early with the "gasified from within"
warning is not a solver failure: a cell below the surface emptied first, and a
fixed mesh parks only its top live cell. Use `use_ale = true` for a material
that regresses (§10.2.1).

If a depletion run is slow or stalls:

- **Check the property laws first.** Sample `k(T)`, `c(T)` of every component on
  [200, 3000] K. The component constructors reject a law that turns
  non-positive on that bracket, but a law passed with `validate_properties =
  false` or a very small positive conductivity at high temperature can still
  turn a cell into an insulator that heats its neighbours to thousands of
  kelvin.
- **Δt collapsing while every live cell is above the threshold** points at a
  stiff transient in the physics, not at depletion. Look at the surface-Newton
  warnings and at `solution.extras.radiation_transmitted` (a transparent
  substrate under Beer–Lambert radiation loses the incident flux through the
  back face).
- **Newton failures right after a park** (visible in the integrator stats).
  The depletion callback restarts the integrator from `DT_AFTER_PARK`
  (10⁻⁶ s), from which the Newton converges on the post-park transient. A run
  that still fails several times at every park has a kernel coupling a live
  cell to a ghost; the equivalence tests
  `test/residual/test_park_equivalence.jl` and
  `test/jacobian/test_park_jacobian.jl` pin every kernel and operator against a
  compacted problem.
- **Newton failures every few steps on a moving front, away from any park.**
  Look for the `solve` warning about a zero gas-transfer coefficient first: a
  condensed component with `λ = 0` stops gas transport in gas-free cells and
  leaves a front the Newton cannot converge on (§3.2.1; Technical Reference
  §5.6). The remaining failures come from the step controller regrowing the
  step beyond the Newton's range after each failure, which the Newton-failure
  step cap of `solve` holds back (§8.3).

Under ALE the default Jacobian for depletion runs is
`Structured(strategy = ApproxSolve())`. It follows the same trajectory as
`jacobian = Structured(strategy = DirectSolve())` at a fraction of the cost,
including on fast ablation fronts, so switching Jacobians is not a remedy for
Newton failures.

## 15.5 Negative or unphysical states

Symptoms: NaN in the temperature or concentration field, `:Unstable` retcode,
concentrations going slightly negative, or oscillations behind a sharp front.

### 15.5.1 Concentrations going negative

The kinetics use a **depletion limiter** — a smooth `tanh(ξ/threshold)` factor
on each reaction rate (default threshold 1.0 kg/m³) — so reactions switch off
as a reactant is exhausted, for every reaction order n ≥ 0; this is what keeps
the system differentiable and stiffness-bounded
(see Technical Reference §6.3). If you still see meaningful negative
concentrations, the usual causes are:

- **Too-coarse a mesh at the reaction front.** Refine near the surface with a
  stretched mesh (§15.5.4) or increase `n_cells`.
- **First-order upwind overshoot in Darcy–Fick advection.** The advective gas
  flux uses first-order upwinding (stable but diffusive; no TVD limiter on the
  Darcy advection — see Technical Reference §9.3). Small negative excursions of a
  near-zero gas concentration are expected and are damped by reaction/transport
  coupling; large ones indicate the pressure field or permeability is driving
  too strong an advective flux (check §15.9).
- **A reactant initialised at exactly zero that a reaction consumes.** Make sure
  `ξ_initial` seeds every reactant the scheme needs.

The smooth temperature gates `½(1+tanh((T−T_min)/w))·½(1+tanh((T_max−T)/w))`, with
`w` the reaction's `gate_width` (default 1 K), leak by roughly ±10·w around `T_min`/`T_max` (Technical Reference §6.2); if a reaction
appears to fire slightly early/late, that leakage — not a bug — is the reason.

### 15.5.2 `ρc_p → 0` blow-up (ALE)

As a cell decomposes to gas and the gas leaves, its `ρc_p^eff` shrinks toward
the pore-gas thermal mass, and `dT/dt = RHS/ρc_p` can blow up. The intended
remedies are surface depletion handling (§15.4) to park the spent surface cell, and
the thickness-termination callback that stops the run when total thickness drops
below `min_thickness`. If you run ALE *without* depletion and hit NaNs near full
burn-out, enable `handle_depletion = true`.

### 15.5.3 Property-table extrapolation

`TableProperty` interpolates linearly and holds its end values outside its range.
A table that stops at 600 K queried at 900 K returns the 600 K value, which can
badly misstate `c_p` or `k` at high temperature. Make sure every `TableProperty` covers
the full temperature span the solve will reach (the Newton boundary solve clamps
surface temperature to `[200, 3000]` K — see Technical Reference §12.2 — so size
tables accordingly).

### 15.5.4 Refine near the surface

Steep thermal and reaction fronts sit near the exposed surface. A geometric
stretch concentrates cells there at fixed total count:

```julia
# Fine cells near z = L (top/surface), coarsening toward z = 0 (substrate):
mesh = create_stretched_mesh(0.01, 100, Val(2), 1.05; stretch_from = :top)
#                            L      n   NC     stretch  direction
```

`create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from = :top, ...)`
uses geometric grading; `stretch > 1` grades the cell sizes across the domain.
Compare against the uniform `create_uniform_mesh(L, n_cells, Val(NC))` to see
whether near-surface refinement removes the overshoot. `solve` does not refine
the mesh during a run; the error indicators and `relocate_nodes!` of Chapter 10
are standalone building blocks.

---

## 15.6 Energy or mass ledger not closing

The cumulative ledgers report `energy_relative` and `mass_max_relative_error`
when they are driven step by step through `record_substep!` (§9.8; the rows
`solve` logs leave these fields `NaN`). The two failure modes below have
*different* causes.

### 15.6.1 `energy_relative` large

This audits assembly and should close to integrator tolerance. If it does not:

- **A callback is mutating state outside the residual.** The mesh is read-only
  during assembly and is synced only by the step-update callback; a custom
  callback that edits `u` or geometry without informing the integrator
  (`derivative_discontinuity!(integrator, true)`) will break the ledger.
- **Depletion masking.** All domain-total sums (`compute_domain_mass`,
  `compute_sensible_energy`) iterate the **active-cell mask**, not
  `eachindex`. If you wrote custom post-processing that loops over all cells
  (including parked ones), your totals — not the solver's — are wrong.
  Always respect `mesh.active_cell_mask` (Technical Reference §16.1).

### 15.6.2 `energy_relative` drifts during active reaction

A small, composition-rate-dependent residual during active pyrolysis is
expected: the ledger omits the heat-capacity-change and ALE volume-change terms
of its `ρc_p·T` proxy (Technical Reference §16.7). A residual that grows with the
*boundary input* rather than with the reaction rate points at the assembly.

---

## 15.7 ALE issues: tangling, GCL, vanishing cross-section

ALE (`use_ale = true`) lets the mesh follow the receding/expanding surface
(Technical Reference §11). Three things go wrong with moving meshes.

### 15.7.1 Mesh tangling and quality

If a swelling/intumescent material moves nodes past each other, cells invert
(the signed width `z_{i+1} − z_i` turns non-positive, which the volume
`A·|z_{i+1} − z_i|` the residual uses does not show). To locate it, check the
node spacing of the saved times, `diff(solution.z[:, i])`, for non-positive
entries (parked slots read `NaN`). For a `Mesh1D` you hold — such as
`solution.problem.mesh`, which carries the final geometry — the quality tools
of §10.6 apply:

```julia
using Pyrolysis.Internal          # Adaptivity exports are not top-level exports

check_quality(mesh, JacobianMetric())        # every live cell has a positive signed width?
quality_violations(mesh, JacobianMetric())   # indices of inverted cells
q = mesh_quality_summary(mesh)               # NamedTuple
q.min_cell_size
q.max_aspect_ratio
```

In a solve the cumulative mesh velocity is computed for you; tangling there
points to an unphysical swelling factor `γ` or shrinkage law producing a
too-large per-step displacement — refine the time step (lower `dtmax`) or
reduce the swelling magnitude.

### 15.7.2 Geometric conservation

The discrete Geometric Conservation Law holds by construction: every residual
evaluation recomputes each cell volume as `V = A·|z_{i+1} − z_i|` from the node
positions in the state (Technical Reference §11.5), so mesh motion injects no
"phantom" mass. A cell driven to
near-zero width is still a degenerate state; address the underlying cause
(§15.7.3).

### 15.7.3 Vanishing cross-section and thickness

In WithChi mode the cross-section evolves as `A(t) = A_0·law(χ̄)` (Technical
Reference §10.5). A shrinkage law `1 − k χ̄` with too large `k` can drive
`A(t) → 0`, at which point volumes vanish and the residual is singular. The
guards keep `law(χ̄) > ε`, but the right fix is a physical shrinkage coefficient.
Likewise, when the *axial* thickness collapses, the thickness-termination
callback fires and the run ends with `:Terminated` — check
`solution.extras.thickness_history` and `min_thickness`.

Note that MLR is reported two ways for variable cross-section:
`solution.mass_loss_rate` is per *current* area `A(t)`, while
`solution.extras.MLR_gauge` is normalised by the initial area `A_0` (the
cone-gauge convention). They coincide in non-ALE / constant-area runs. If your
MLR looks off by a factor that tracks `A(t)/A_0`, you are reading the wrong one.

### 15.7.4 Choosing a linear-solve strategy for large ALE

For large ALE problems the structured ALE couplings (cumulative mesh velocity,
χ back-reaction) dominate Jacobian assembly/solve cost. The two strategies trade
exactness for speed (Technical Reference §15.10):

| Strategy | Exactness | When to use |
|---|---|---|
| `DirectSolve()` (fixed-mesh default) | Exact; densify + factor | Fixed mesh; exact ALE Jacobian; verification anchor |
| `ApproxSolve()` (ALE default) | Approximate; drops structured couplings (`is_approximate == true`) | ALE: the modified Newton iteration still converges |

Select one by passing the bare backend; it works with or without
`handle_depletion`:

```julia
solution = solve(problem;
    use_ale  = true,
    jacobian = Structured(strategy = DirectSolve()),   # the exact ALE Jacobian
)
```

The two strategies follow the same trajectory, and the exact one costs several
times more per step, so `ApproxSolve` is the default; use
`DirectSolve` to check a configuration against the exact Jacobian. Use
`verify_jacobian` (§15.3) with the masked mode to confirm the approximate
backend still tracks the reference inside its support.

---

## 15.8 Radiation surprises

The four radiation models are mutually exclusive and selected by the
`radiation_model` kwarg (default `NO_RADIATION`). Three of them are reachable
from `solve`; the fourth is quarantined.

### 15.8.1 Model selection and mutual exclusivity

```julia
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)  # opaque, O(1)
solution = solve(problem; radiation_model = BEER_LAMBERT)        # in-depth absorption
# radiation_model = NO_RADIATION  is the default
```

| Model | Behaviour | Cost / Jacobian |
|---|---|---|
| `NO_RADIATION` | No volumetric source | O(1), tridiagonal |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the surface | O(1), tridiagonal |
| `BEER_LAMBERT` | Exponential in-depth attenuation `I(z)` | O(n), upper-triangular (non-local) Jacobian |
| `P1_QUASI_STEADY` | P1 diffusion with re-radiation | **Experimental — rejected by `solve`** |

`NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT` are exported by
`using Pyrolysis`. `P1_QUASI_STEADY` is exported from
`Pyrolysis.Experimental` (its residual lives in `Pyrolysis.Experimental.RadiationP1`)
and is **not** part of the core unified residual; passing it to `solve` errors with a message pointing you to the
Experimental module. Do not expect to use it through the standard `solve` path
(Technical Reference §8.12).

### 15.8.2 "I selected a radiation model but nothing absorbs"

`solve` warns about the two common mismatches:

- **`BEER_LAMBERT` with no absorption coefficients.** Beer–Lambert needs a
  non-zero mass absorption `α_j` on the absorbing components (the volumetric
  `α_eff = Σ_j ξ_j α_j`; Technical Reference §8.3). If every `α_j` is zero, the
  optical thickness `τ_i = α_i Δz_i` is ~0, `1 − exp(−τ) ≈ 0`, and effectively
  no radiation is deposited. Set `α` on the solid (e.g.
  `SolidComponent(:char, …, α = 5.0)`).
- **`SURFACE_ABSORPTION` with non-zero `α`.** The coefficients are *ignored* —
  the absorbed flux is applied at the surface. The warning tells you the in-depth profile
  you specified is not being used; switch to `BEER_LAMBERT` if you want depth
  penetration.

### 15.8.3 Optical-thickness regime

Choose the model by optical thickness `τ = α L` (Technical Reference §8.9):

| `τ = α·L` | Regime | Recommended |
|---|---|---|
| > 10 | Optically thick | `SURFACE_ABSORPTION` (or P1 if re-radiation matters) |
| 1–10 | Intermediate | `BEER_LAMBERT` |
| 0.1–1 | Optically thin | `BEER_LAMBERT` |
| < 0.1 | Transparent | `NO_RADIATION` |

Using `SURFACE_ABSORPTION` on a semi-transparent sample concentrates the
absorbed flux at the surface and overheats it; using `BEER_LAMBERT` on an opaque char
just adds non-local Jacobian cost for no benefit.

### 15.8.4 Time-varying incident flux is always tracked

For a **time-varying** incident flux (e.g. a cone ramp
`RadiativeFluxBC(flux = t -> 25e3 * (1 - exp(-t/10)))`), the flux is re-sampled
at every RHS evaluation whenever `radiation_model` is `SURFACE_ABSORPTION` or
`BEER_LAMBERT` — it is never snapshotted at `t = 0`, and there is no constancy
probe. A flux that is flat early and varies later (a delayed step at `t = 5`,
say) is tracked correctly too. The cost is one call of your flux function per
RHS evaluation, which is negligible; no special arrangement is needed.

### 15.8.5 Marshak ε coefficient (P1 only)

If you do reach into `Pyrolysis.Experimental.RadiationP1`, note the Marshak
boundary relation is `q = c·(G_ext + 4·I_ext − G)` with `c = ε/(2(2−ε))`
(Modest Eq. 16.48): the exchange coefficient is **not** `0.5ε`, and the external
flux carries the half-range weight `4c = 2ε/(2−ε)` — the same weight as the
ambient term — **not** 1. Either simplification silently under-injects radiation
for `ε < 1`. This is internal to the experimental P1 path; it does not affect
the core models reachable from `solve`.

---

## 15.9 Gas-transport and pressure surprises (Darcy mode)

Darcy–Fick transport activates when the material carries permeability (any
solid/liquid component with `κ > 0`); see UG-7 and Technical Reference §9.

- **No pressure-driven flow.** The face Darcy velocity is
  `v_g = −mob_face·∇P` with `mob_face` the distance-weighted harmonic mean of
  the two cells' mobilities `κ/μ`; the face is clamped to zero when either
  side's mobility falls below `1e-30` (an impermeable neighbor blocks the
  face) or a cell's viscosity below `1e-20` Pa·s. If you expect advection and
  see none, check that the *solid/liquid* components carry `κ` (gases never
  do) and that a pressure gradient actually develops.
- **Pressure spikes in nearly-solid cells.** The ideal-gas pressure carries a
  porosity divisor `P = (Σ_g ξ_j R_g T / M_j)/φ` with a floor `φ ≥ 1e-12`. In a
  cell that is almost entirely solid, the same gas mass yields a high pressure;
  that is physical (confined gas), not a bug.
- **Pressure BCs vent gas advectively — and only advectively.**
  `DirichletPressureBC`, `NeumannPressureBC`, and `ConvectivePressureBC` set
  the boundary-face Darcy velocity, which carries the advective venting flux
  `J = v_b·ξ/φ` out of the surface cell (outflow only — inflow advects the
  zero-concentration ambient). They contribute no diffusive film-law flux:
  that channel is governed by the diffusive mass BCs (`ConvectiveMassBC`,
  `DiffusiveMassBC`, …). Combine the two with `+`
  (`DiffusiveMassBC(h_m) + DirichletPressureBC()`) to get both escape paths;
  a pressure BC alone vents at a rate set by the permeability and the
  pressure gradient.
- **All gas species move together.** One effective `λ` serves every gas
  species, and advection carries them all at the same bulk velocity — there is
  no per-species molar-mass scaling (Technical Reference §9.2–9.3).

---

## 15.10 Parameter sweeps and parallelism

`parameter_sweep(problem, grid; modify, extract, parallel, solve_kwargs...)`
runs a grid of solves. The contract that keeps it correct and fast:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]

# modify MUST return a fresh problem and MUST NOT mutate the base problem.
# Material and PyrolysisProblem are built with keyword constructors — there is
# no copy/update constructor, so rebuild from the base problem's fields:
function modify(prob, p)
    base = prob.material
    new_rxn = Reaction(:decomp, :virgin => :gas; A = p.A, E = p.E, h = 200e3)
    new_mat = Material(;
        name                  = base.name,
        components            = base.components,
        reactions             = (new_rxn,),
        mixing                = base.mixing,
        lateral_shrinkage_law = base.lateral_shrinkage_law,
        depletion_limiter     = base.depletion_limiter,
        gas_viscosity         = base.gas_viscosity,
    )
    return PyrolysisProblem(;
        mesh         = prob.mesh,
        material     = new_mat,
        bc_set       = prob.bc_set,
        T_initial    = prob.T_initial,
        ξ_initial    = prob.ξ_initial,
        tspan        = prob.tspan,
        output_times = prob.output_times,
    )
end

# extract pulls small scalars/curves to minimise worker→head transfer:
extract(sol) = (peak_T = maximum(sol.surface_T),
                final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```

Failure modes:

- **`modify` mutating the base problem.** The immutable `ProblemDef`/residual
  are shared across workers; the mutable `Workspace` is per-worker. If `modify`
  edits the shared base in place, you corrupt every other point. Always build a
  new problem.
- **Returning full solutions from a big grid.** `extract` defaults to `identity`
  (the whole solution). For large grids that floods the head process — return a
  small NamedTuple.
- **Warm-start in a sweep.** You can pass `warm_start` through `solve_kwargs`,
  but reusing a neighbour's converged state across perturbed grid points can
  introduce artifacts in *sensitivity* studies (Technical Reference §17.7).
  Prefer cold starts when computing gradients.

`parallel` defaults to `nworkers() > 1` (uses `Distributed.pmap`); set
`parallel = false` to force sequential `map` for debugging. `parameter_sweep`
already runs each solve with `progress = false`.

---

## 15.11 Performance checklist

Work top to bottom; each item is cheap to try and the early ones dominate.

1. **Mesh size.** Cost scales with `n·(1+NC)` (+ the ALE `z` block and the
   WithChi `χ` block). Run a quick convergence check and use the smallest
   `n_cells` that resolves the surface temperature and MLR. Grade toward the
   surface with `create_stretched_mesh` (§15.5.4) rather than refining
   uniformly.
2. **Keep the default structured Jacobian.** `jacobian = nothing` →
   `Structured(strategy = DirectSolve())` (`ApproxSolve()` under ALE, which
   drops the dense geometry block and keeps the pattern local), sparse and
   analytical — far faster than finite-differencing the right-hand side.
3. **Keep the default linear solver.** `linear_solver = nothing` picks UMFPACK
   for the exact ALE Jacobian and KLU for every other sparse Jacobian (§8.4).
   Iterative
   `ILUGMRESFactorization(τ = 0.1, gmres_rtol = 1e-8, restart = 30)` and
   `ILUBiCGSTABFactorization(τ = 0.1, rtol = 1e-8)` exist for very large
   systems, but inexact solves can cause Newton-step rejections on stiff
   transients, and `solve` warns when you select ILU-GMRES. Treat ILU as a last
   resort for problems where direct factorization truly dominates.
4. **Loosen tolerances for exploration.** `abstol`/`reltol` are the most direct
   speed knob; tighten only for validation (§15.2.2).
5. **Choose the integrator.** `KenCarp4` is the default; a BDF (`FBDF`) can take
   larger steps on very stiff problems (§15.2.1).
6. **Right-size the radiation model.** `SURFACE_ABSORPTION` is O(1);
   `BEER_LAMBERT` adds an O(n) sweep and a non-local (upper-triangular) Jacobian
   coupling. Do not pay for in-depth absorption you do not need (§15.8.3).
7. **Warm-start continuation runs** (§15.2.4) and **parallelise sweeps**
   (§15.10).
8. **Disable extras you are not using.** `progress = false` removes the progress
   bar in batch/sweep contexts; leave `diagnostics = false` in production runs —
   it costs an extra energy sweep every `diagnostics_stride` steps and is only
   needed for debugging (§15.1.2).
9. **Keep one Julia session for many solves.** The first solve of a problem
   compiles the residual and Jacobian kernels for its types: the component and
   property types (every `CallableProperty` closure is a type of its own), the
   reactions, the boundary-condition combination, the component count and the
   radiation and transport modes. That takes 15–20 s, or about 2 s when the
   types match a precompiled configuration (constant properties with the cone
   boundary conditions). The time integrator compiles once per solver
   configuration and comes precompiled for the defaults: a fixed or ALE mesh,
   with or without `handle_depletion`, `KenCarp4` with the default
   factorization. Another integrator, a user `callback`, `diagnostics = true` or
   `verify_jacobian` compiles it once more per session, at about the cost of a
   new problem type. Changing values (fluxes, kinetic parameters, mesh size,
   tolerances) never recompiles. Re-evaluating a property closure at the REPL
   does, because the new closure is a new type; a sweep's `modify` function
   builds the same closure type at every point.

The package's allocation-free residual caches, `@generated` component loops, and
per-operator AD workspaces (Technical Reference §15.7, §15.11) are automatic; you
do not configure them. The levers above are the ones you control.

---

## 15.12 Where to look next

- Conservation ledgers and what each residual proves: Technical Reference §16;
  the post-processing/reading side in UG-9.
- Jacobian backends, strategies, linear solvers, and the extension stubs
  (`forward_sensitivity`, `adjoint_sensitivity`): UG-11 and Technical
  Reference §15.
- Depletion handling, and the standalone mesh tools (error indicators,
  monitors, `relocate_nodes!`, quality metrics): UG-10 and Technical
  Reference §14.
- The full keyword and default table for `solve`: UG-16.

The escape hatch for anything not on the curated top-level surface is the
`Pyrolysis.Internal` umbrella (`using Pyrolysis.Internal`), which re-exports
every submodule's public symbols — the mesh-quality and diagnostics helpers in
this chapter live there. Symbols reached that way are not part of the stable
public API and may move between minor releases.
