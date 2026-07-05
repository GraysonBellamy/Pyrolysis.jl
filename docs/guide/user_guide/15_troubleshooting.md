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
conservation ledgers.

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
| `:Terminated` | A callback stopped the run | Sample fully depleted (ALE thickness termination) or a user callback called `terminate!` |
| `:MaxIters` | Hit `maxiters` before `t_end` | Stiff transient; step size collapsed; under-resolved front |
| `:DtLessThanMin` | Step fell below `dtmin` | Discontinuity, NaN in the RHS, or genuine stiffness wall |
| `:Unstable` / `:DtNaN` | NaN/Inf in the state | Unphysical state (negative concentration, `ρc_p → 0`), bad property table extrapolation |

A `:Terminated` run is *not* an error in ALE mode — it usually means the sample
burned away. Check `solution.extras.thickness_history[end]` against
`min_thickness` (default `1e-6` m).

### 15.1.2 Turn on the conservation ledgers

The single most useful debugging switch is `diagnostics = true`. It installs a
callback that snapshots the mass and energy ledgers every `diagnostics_stride`
accepted steps (default 10) and stores typed rows in
`solution.extras.diagnostics_log`:

```julia
solution = solve(problem;
    diagnostics        = true,
    diagnostics_stride = 5,     # snapshot every 5 accepted steps
    progress           = false,
)

for row in solution.extras.diagnostics_log
    println("t=$(row.t)  E_matrix=$(row.E_matrix)  ",
            "r_disc(rel)=$(row.discrete_relative)  ",
            "r_phys(rel)=$(row.physical_relative)  ",
            "ΔE_gas=$(row.gas_storage_gap)  ",
            "m_err=$(row.mass_max_relative_error)")
end
```

Each `DiagnosticsRow` carries `t`, `E_matrix`, `E_total`, `discrete_residual`,
`discrete_relative`, `physical_residual`, `physical_relative`,
`gas_storage_gap`, and `mass_max_relative_error`. How to read them (see
Technical Reference §16):

- **`discrete_relative`** audits the *bookkeeping*. It is assembled from exactly
  the terms the RHS computes, so it closes to integrator tolerance (~`1e-6`
  relative) **by construction**. If it is large, you have a bug in flux/source
  assembly or a callback — not a physics problem. (Caveat: with strongly
  time-varying composition it cannot close exactly, because the ledger omits the
  `∫ T ∂(ρc_p)/∂t` and ALE volume-change terms; expect a small floor.)
- **`physical_relative`** audits the *quasi-steady-gas approximation*. For
  canonical 1-atm pyrolysis it should stay **< 1 %** of initial energy. A larger
  value means the regime stresses the gas-storage assumption.
- **`gas_storage_gap`** (`ΔE_gas = E_total − E_matrix`) is the sensible heat
  carried by gas that the matrix-only `ρc_p^eff` does not store. Canonically
  < 0.2 % of `E_matrix`.
- **`mass_max_relative_error`** is the worst per-component mass-balance error;
  it should be `< 1e-10`.

`diagnostics` is **snapshot-based**, not continuous — it gives you a coarse time
history adequate for spotting where a ledger drifts, not a dense trace. The
overhead is one extra energy/mass sweep every `diagnostics_stride` steps.

---

## 15.2 Stiff or slow solves

Arrhenius kinetics make the system stiff during ignition and the early reaction
transient; conduction adds a parabolic time scale. If a solve is slow, hits
`:MaxIters`, or grinds to `:DtLessThanMin`, work through the following.

### 15.2.1 Confirm the default integrator is appropriate

The default integrator is `KenCarp4` (an SDIRK method) with a sparse direct
linear solve; its `autodiff = AutoFiniteDiff()` setting only takes effect on
the no-Jacobian (depletion) path, since `solve` normally supplies the
structured analytical Jacobian (§15.2.3). When both `integrator` and
`linear_solver` are left at `nothing`, `solve` builds:

```julia
# Internal default:
# KenCarp4(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
```

`KenCarp4` is a good default for moderately stiff pyrolysis. For *very* stiff
problems (high pre-exponential `A_i`, sharp ignition) a BDF method often takes
larger steps. Integrators live in `OrdinaryDiffEq`, which is **not** re-exported
by `using Pyrolysis`; load it (or the relevant sublibrary) yourself and pass a
constructor or a fully-built algorithm to the `integrator` kwarg:

```julia
using OrdinaryDiffEq            # for FBDF, Rodas5P, etc.

# (a) pass a constructor — solve() injects autodiff + linsolve for you:
solution = solve(problem; integrator = FBDF)

# (b) pass a constructor and a linear solver together:
solution = solve(problem;
    integrator    = FBDF,
    linear_solver = SparspakFactorization(),   # re-exported by Pyrolysis.Solver
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

When `jacobian = nothing` (the default) and `handle_depletion = false`, `solve`
builds a `Structured(strategy = DirectSolve())` analytical Jacobian and hands a
sparse prototype to the integrator. This is much faster than the integrator's
own colored finite-difference Jacobian because it is exact and sparsity-aware
(see Technical Reference §15). You normally do not need to touch it. The cases
where you *do* are covered in §15.3 (verification), §15.4 (depletion), and
§15.7 (large ALE).

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
problems must match (enforced by a length check). `tspan[1]` still sets the
integrator's start time, so the prior final state is treated as fresh initial
data.

---

## 15.3 Jacobian mismatch warnings

If the integrator behaves erratically only when the structured Jacobian is
active — many rejected steps, sudden step-size collapse that goes away when you
fall back to the dense reference — suspect a Jacobian that disagrees with the
residual.

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
directly. You need the `ProblemDef`/`Workspace` pair, which you build from the
problem:

```julia
def = to_problem_def(problem; use_ale = false, radiation_model = NO_RADIATION)
ws  = build_workspace(problem, def)
u   = create_state_vector(ws.layout, ws.mesh)

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

`JacobianSpec` also accepts a `PyrolysisProblem` directly
(`JacobianSpec(backend, problem)`) if you do not already have a `def`/`ws` pair.

### 15.3.3 Support-pattern (sparsity) issues

The comparison `mode` matters:

- **`:strict`** — for exact backends (`DirectSolve`, `StructuredSolve`). Both
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

## 15.4 Depletion and the pluggable-Jacobian incompatibility

Surface-cell depletion (`handle_depletion = true`, ALE only) resizes the state
vector mid-solve via `resize!(integrator, new_length)` when a surface cell
thins below `depletion_threshold` of its initial thickness (default 5 %) and
merges it into its interior neighbour (see Technical Reference §14.3).

The pluggable `JacobianSpec` cache is **immutable** and built once at setup, so
it cannot survive a resize. The two are therefore mutually exclusive:

```julia
# ERROR: solve() rejects this combination up front.
solve(problem;
    use_ale          = true,
    handle_depletion = true,
    jacobian         = JacobianSpec(Structured(strategy = DirectSolve()), problem),
)
# → "the `jacobian` kwarg is incompatible with `handle_depletion = true`"
```

The supported way to run depletion is to **let `solve` pick the Jacobian**: when
`handle_depletion = true`, the auto-default returns `nothing` and the integrator
falls back to its own colored finite-difference Jacobian, which *does* survive a
`resize!`:

```julia
solution = solve(problem;
    use_ale             = true,
    handle_depletion    = true,
    depletion_threshold = 0.05,   # merge below 5% of initial thickness
    min_cells           = 2,      # stop merging at 2 cells
    min_thickness       = 1e-6,   # terminate when total thickness < this
    progress            = false,
)
```

The merge is energy-conserving (it does not blow up when a near-empty cell with
`ρc_p → 0` is absorbed; Technical Reference §14.3). Consequence: with depletion
you give up the fast structured Jacobian. If a depletion run is slow, see
§15.4.1 and the linear-solver choice in §15.11.

### 15.4.1 The MKL `ipiv` resize bug

With `handle_depletion = true` and a large state vector (> ~500 elements),
passing a fully built integrator whose `linsolve` is left unset lets
LinearSolve.jl's auto-selection pick an MKL LU factorization whose pivot array
(`ipiv`) is not resized correctly on `resize!`, producing a
`DimensionMismatch: ipiv has length X, but needs to be Y`. `solve` warns when
it detects this risky combination.

The fix is to let `solve` resolve the linear solver: on the depletion path it
defaults to Julia's native `LUFactorization()`, which handles `resize!`
correctly (a sparse factorization cannot run here anyway — there is no sparse
Jacobian prototype — so an explicit `SparspakFactorization`/`KLUFactorization`
is rebuilt to `LUFactorization` with a warning):

```julia
solution = solve(problem;
    use_ale          = true,
    handle_depletion = true,
    progress         = false,
)
```

For the optional `KLU`, `Pardiso`, and AMG extensions and how to load them, see
UG-11. The `klu_factorization`, `pardiso_factorization`, and `amg_preconditioner`
stubs error with a `MethodError` until you `import` the matching package — that
error *is* the prompt to load the dependency.

---

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

The smooth temperature gates `½(1+tanh(T−T_min))·½(1+tanh(T_max−T))` leak by
roughly ±10 K around `T_min`/`T_max` (Technical Reference §6.2); if a reaction
appears to fire slightly early/late, that leakage — not a bug — is the reason.

### 15.5.2 `ρc_p → 0` blow-up (ALE)

As a cell decomposes to gas, its matrix-only `ρc_p^eff` (gas is excluded by
Formulation A) shrinks, and `dT/dt = RHS/ρc_p` can blow up. The intended
remedies are ALE surface depletion (§15.4) to merge the spent surface cell, and
the thickness-termination callback that stops the run when total thickness drops
below `min_thickness`. If you run ALE *without* depletion and hit NaNs near full
burn-out, enable `handle_depletion = true`.

### 15.5.3 Property-table extrapolation

`TableProperty` interpolates linearly and extrapolates beyond its range. A table
that stops at 600 K queried at 900 K will extrapolate — possibly to a negative
`c_p` or `k`, which then produces NaN. Make sure every `TableProperty` covers
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
whether near-surface refinement removes the overshoot. For dynamic refinement
options (error indicators, `relocate_nodes!`), see UG-10.

---

## 15.6 Energy or mass ledger not closing

Use the diagnostics rows (§15.1.2) to localise. The two energy residuals fail
for *different* reasons.

### 15.6.1 `discrete_relative` large

This audits assembly and should close to integrator tolerance. If it does not:

- **A callback is mutating state outside the residual.** The mesh is read-only
  during assembly and is synced only by the step-update callback; a custom
  callback that edits `u` or geometry without informing the integrator
  (`u_modified!`) will break the ledger.
- **AMR / depletion masking.** All domain-total sums (`compute_domain_mass`,
  `compute_sensible_energy`) iterate the **active-cell mask**, not
  `eachindex`. If you wrote custom post-processing that loops over all cells
  (including merged-away ones), your totals — not the solver's — are wrong.
  Always respect `mesh.active_cell_mask` (Technical Reference §16.1).

### 15.6.2 `physical_relative` large (> 1 %)

This audits the quasi-steady-gas approximation. A large value means the regime
genuinely stresses the assumption (gas residence time approaching the thermal
time scale — unusual at 1 atm). Reconsider the setup rather than the numerics;
`gas_storage_gap` will corroborate (it will be a non-trivial fraction of
`E_matrix`).

### 15.6.3 The S_conv vs transpiration double-count

This is the most common "my energy balance is off by a constant rate" trap. The
default energy form (`energy_form = :standard`) carries the gas sensible
enthalpy through the **always-on** volumetric convection source
`S_conv = −Σ_g c_{p,g} J̄_g ∂T/∂z` (ThermaKin Eq. 17 convention; Technical
Reference §7.4). The surface **transpiration** enthalpy term would carry the
*same* enthalpy at the boundary — enabling both double-counts by exactly
`∫ S_conv`. The package forbids the combination:

```julia
# ERROR: solve() rejects this.
solve(problem; use_transpiration_bc = true)
# → "use_transpiration_bc=true double-counts gas convective enthalpy
#    (S_conv is always active ...). Use use_transpiration_bc=false"
```

So `use_transpiration_bc` is `false` by default and there is no supported way to
turn it on through `solve`. Mass transport is unaffected either way — gas still
leaves the domain through the mass BCs; only the *enthalpy bookkeeping* is at
issue. (Parameter sets calibrated against a solver that applies a transpiration term
on top of a volumetric advection source embed the double-count and will predict
higher `T_s` and MLR here; re-validate imported calibrations.)

### 15.6.4 `:standard` vs `:with_generation_sink`

The optional gas-generation sink `S_gen = −Σ_g Δh_g (∇·J_g)` (Technical Reference
§7.5, §16.4) can be added:

```julia
solution = solve(problem; energy_form = :with_generation_sink)
```

`energy_form` accepts only `:standard` (default) or `:with_generation_sink`
(anything else errors). **Datum warning (see §7.7):** the sink is consistent only
with `T_ref`-datum (formation-style) reaction heats. With DSC-style heats quoted at
the reaction temperature — the usual meaning of calibrated heats — the gas sensible
enthalpy is already inside `h`, and the sink subtracts it a *second* time
(≈0.7–1.2 MJ per kg of volatiles at 900 K); `solve` emits a one-time warning when it
is selected. Even with the right datum it is only a *partial* step toward
conservative form (gas heat storage and condensed-product sensible generation remain
outside the assembly). The midpoint-rule `Δh_g` makes `∫(S_conv + S_gen)` telescope
exactly against the surface outflow, which is what the ledger audit checks. For most
cone/TGA work `:standard` is correct and matches FDS/Gpyro/ThermaKin.

---

## 15.7 ALE issues: tangling, GCL, vanishing cross-section

ALE (`use_ale = true`) lets the mesh follow the receding/expanding surface
(Technical Reference §11). Three things go wrong with moving meshes.

### 15.7.1 Mesh tangling and quality

If a swelling/intumescent material moves nodes past each other, cells invert
(`V ≤ 0`) and the geometry recompute produces NaNs. For diagnosis, the
mesh-motion utilities expose quality checks. These operate on a `Mesh1D` and
live in the `Pyrolysis.Internal` umbrella (they are not top-level exports):

```julia
using Pyrolysis.Internal          # mesh-motion helpers are not top-level exports

q = check_mesh_quality(mesh)       # NamedTuple
q.valid                            # all cells V > 0 ?
q.min_volume
q.max_aspect_ratio
q.tangled_cells                    # indices with V ≤ 0

mesh_is_valid(mesh)                # quick boolean
```

`move_mesh_safe!` rejects a motion (and rolls back) if it would tangle or exceed
an aspect-ratio bound (default 10). If you are driving the mesh yourself, prefer
it over the raw `move_mesh!`. For production solves the cumulative mesh velocity
is computed for you; tangling there points to an unphysical swelling factor `γ`
or shrinkage law producing a too-large per-step displacement — refine the time
step (lower `dtmax`) or reduce the swelling magnitude.

### 15.7.2 GCL violations

The discrete Geometric Conservation Law,
`V_new − V_old = Δt (w_R − w_L) A`, must hold to machine precision or the mesh
motion injects "phantom" mass (Technical Reference §11.5). Per-step GCL error
should be `O(1e-12)`–`O(1e-15)`; a value above `1e-10` signals geometric
degeneracy. To check explicitly:

```julia
using Pyrolysis.Internal

err = gcl_error(mesh, w, dt)        # predictive relative error (max over cells)
check_gcl(mesh, w, dt; tol = 1e-12) # Bool

# For a time-varying cross-section (WithChi mode) pass A_old, A_new:
gcl_residual(mesh, V_old, w, dt; A_old = A0, A_new = A1)
```

A `GCLTracker` accumulates the error over a run if you wire it into a custom
driver. In normal use, a GCL violation is not something you tune — it is a
symptom of a degenerate mesh (a cell driven to near-zero width); address the
underlying cause (§15.7.3).

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
χ back-reaction) dominate Jacobian assembly/solve cost. The strategies trade
exactness for speed (Technical Reference §15.10):

| Strategy | Exactness | When to use |
|---|---|---|
| `DirectSolve()` (default) | Exact; densify + factor | Small/medium problems; verification anchor |
| `StructuredSolve()` | Exact; sparse LU on the local block + Richardson (2–5 iters) | Typical ALE problems |
| `ApproxSolve()` | Approximate; drops structured couplings (`is_approximate == true`) | Large ALE where Newton tolerance can absorb a Schur-style approximation |
| `MatrixFreeSolve()` | GMRES on the operator action | Research/experimental |

Select via a `JacobianSpec` (remember: **not** with `handle_depletion = true`):

```julia
solution = solve(problem;
    use_ale  = true,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
)
```

If `ApproxSolve` increases Newton rejections on Arrhenius transients, fall back
to `StructuredSolve` (exact, still O(N) per coupling apply) or the default
`DirectSolve`. Use `verify_jacobian` (§15.3) with the masked mode to confirm an
approximate backend still tracks the reference inside its support.

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
`using Pyrolysis`. `P1_QUASI_STEADY` has moved to
`Pyrolysis.Experimental.RadiationP1` and is **not** part of the core unified
residual; passing it to `solve` errors with a message pointing you to the
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
   `Structured(DirectSolve())`, exact and sparse — far faster than colored FD.
   Only depletion forces the FD fallback (§15.4).
3. **Pick the right linear solver.** `SparspakFactorization` (the default) is
   ~28 % faster than KLU for the block-tridiagonal pattern. For large ALE, try
   `Structured(StructuredSolve())`, then `ApproxSolve()` (§15.7.4). Iterative
   `ILUGMRESFactorization(τ = 0.1, gmres_rtol = 1e-8, restart = 30)` and
   `ILUBiCGSTABFactorization(τ = 0.1, rtol = 1e-8)` exist for very large
   systems, but inexact solves cause Newton-step rejections on stiff transients
   — `solve` warns that Sparspak is usually 1.7–3.3× faster in practice. Treat
   ILU as a last resort for problems where direct factorization truly dominates.
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
   it costs an extra ledger sweep every `diagnostics_stride` steps and is only
   needed for debugging (§15.1.2).

The package's allocation-free residual caches, `@generated` component loops, and
per-operator AD workspaces (Technical Reference §15.7, §15.11) are automatic; you
do not configure them. The levers above are the ones you control.

---

## 15.12 Where to look next

- Conservation ledgers and what each residual proves: Technical Reference §16;
  the post-processing/reading side in UG-9.
- Jacobian backends, strategies, and the extension stubs
  (`klu_factorization`, `pardiso_factorization`, `amg_preconditioner`,
  `forward_sensitivity`, `adjoint_sensitivity`, `verify_mixing_derivatives`):
  UG-11 and Technical Reference §15.
- Adaptivity (error indicators, monitors, `relocate_nodes!`, depletion):
  UG-10 and Technical Reference §14.
- The full keyword and default table for `solve`: UG-16.

The escape hatch for anything not on the curated top-level surface is the
`Pyrolysis.Internal` umbrella (`using Pyrolysis.Internal`), which re-exports
every submodule's public symbols — the mesh-quality, GCL, and remapping helpers
in this chapter live there. Symbols reached that way are not part of the stable
public API and may move between minor releases.
