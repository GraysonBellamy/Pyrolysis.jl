# 10. Adaptive Meshing Options

Pyrolysis.jl gives you several ways to let the mesh follow the physics instead of
fixing it once at problem setup. The most important of these — and the only one
wired directly into the production `solve` path — is **surface-cell depletion
handling**, which merges the receding surface cell into its neighbor as the solid
gasifies and the Arbitrary Lagrangian–Eulerian (ALE) mesh shrinks. Beyond that,
the `Pyrolysis.Adaptivity` submodule exposes standalone tools for **r-refinement**
(moving nodes to concentrate them where the solution is steep), **error
indicators** and **monitor functions** that drive those moves, and **mesh-quality
metrics** to keep the mesh well-conditioned. A fully topology-changing h-AMR
cluster (cell splitting and sibling coarsening) also ships, but it lives in
`Pyrolysis.Experimental` and is not connected to `solve`. This chapter shows how
to use each, what every option means and defaults to, and where the sharp edges
are. For the underlying theory (equidistribution, energy-conserving merges, the
geometric conservation law) see Technical Reference §14 (AMR), §11 (ALE), and
§13 (finite-volume discretization).

> **Coordinate reminder.** `z = 0` is the bottom / substrate (boundary tag
> `:bottom`, id 2); `z = L` is the top / exposed surface (tag `:top`, id 1). Heat
> enters from the top and the material recedes downward, so the "surface cell" is
> the last (highest-`z`) active cell. See Chapter 6 (Meshes and Geometry).

---

## 10.1 When to use adaptivity

Adaptivity costs assembly time and, in some configurations, forces a less
efficient Jacobian. Reach for it only when the physics warrants it:

| Situation | Symptom | Recommended tool |
|---|---|---|
| Surface recession (charring/ablating solid in ALE) | Surface cell thins toward zero; timestep collapses; `ρc_p → 0` trouble | `handle_depletion = true` (§10.2) |
| Sample fully consumed | Total thickness approaches zero | thickness termination (automatic in ALE; §10.3) |
| Steep, *moving* thermal or reaction front on a fixed cell count | Front under-resolved between cells, but you don't want more cells | r-refinement via `relocate_nodes!` (§10.5) |
| You suspect mesh distortion after ALE motion | Large cell-size jumps, near-zero or inverted volumes | quality metrics + `enforce_quality!` (§10.6) |
| You want true local h-refinement | — | Experimental h-AMR (§10.7), with caveats |

A good default workflow is: start with a **uniform or surface-stretched fixed
mesh** (Chapter 6) sized from a convergence study (Technical Reference §18.5),
and turn on depletion handling only once you run an ALE case where the surface
actually recedes. Node relocation and h-AMR are advanced, opt-in tools; most
cone-calorimeter and TGA studies never need them.

---

## 10.2 Surface-cell depletion handling (the production path)

When a charring solid pyrolyzes under ALE mesh motion (`use_ale = true`), the
exposed surface cell loses mass and its node moves downward. Eventually that cell
becomes so thin that its discrete heat capacity collapses and the integrator is
forced into vanishingly small steps. **Depletion handling** detects such thin
cells and *merges* them into a neighbor, removing the singularity while
conserving mass and energy.

### 10.2.1 Turning it on through `solve`

This is the recommended way to use depletion handling — you do not build the
callback yourself. Pass `handle_depletion = true` (which **requires**
`use_ale = true`) to `solve`:

```julia
using Pyrolysis

sol = solve(problem;
    use_ale             = true,    # node positions become state; mesh can shrink
    handle_depletion    = true,    # merge thin surface cells instead of stalling
    depletion_threshold = 0.05,    # merge a cell once it drops below 5% of its
                                   #   initial thickness (default; must be in (0,1))
    min_cells           = 2,       # never merge below this many cells; terminate instead
    min_thickness       = 1e-6,    # [m] terminate when total thickness < this (default)
    radiation_model     = SURFACE_ABSORPTION,
)
```

Defaults and meaning of each depletion-related keyword (verified against
`solve(...)` in `src/Solver/solve.jl`):

| Keyword | Default | Meaning |
|---|---|---|
| `use_ale` | `false` | Enable ALE mesh motion (node positions enter the state vector). Prerequisite for depletion. |
| `handle_depletion` | `false` | Detect thin cells and merge them mid-solve. |
| `depletion_threshold` | `0.05` | Cell is "thin" when current thickness `< threshold ×` its own initial thickness (per-cell reference). Must satisfy `0 < threshold < 1`. |
| `min_cells` | `2` | Once the active cell count reaches this, the run terminates rather than merging further. Must be `≥ 1`. |
| `min_thickness` | `1e-6` | [m] Total-sample thickness below which the integrator terminates (separate from per-cell merging). Must be `> 0`. |

Each callback invocation performs **at most one merge** so the integrator can
re-adapt its step between topology changes; thin cells are processed
surface-first (highest `z` first). The merge is **energy-conserving** — see
§10.2.3.

### 10.2.2 What `solve` wires up under the hood

When `use_ale = true && handle_depletion = true`, `solve` assembles a
`CallbackSet` containing (among others):

1. a **step-update callback** that mirrors the integrator's accepted state back
   into the mutable mesh (so subsequent callbacks see current geometry);
2. a **thickness-termination callback** (§10.3);
3. a **depletion callback** that finds thin cells, merges one, resizes the
   workspace caches, and calls `resize!(integrator, …)` to shrink the ODE state.

You can build these manually if you are driving the integrator yourself, but for
ordinary use the `solve` keywords above are all you need. The public constructors
(all exported by `using Pyrolysis`) are:

```julia
# Build depletion callback state (tracks layout + per-cell initial thicknesses + merge log)
cb_state = DepletionCallbackState(layout, mesh)         # convenience: thicknesses from mesh
# or  DepletionCallbackState(layout, initial_cell_thickness::Vector{Float64})  # one entry per cell

# Build the depletion callback itself
dep_cb = create_depletion_callback(ws, mesh, material, cb_state;
    threshold = 0.05,   # same meaning as depletion_threshold above
    min_cells = 2)

# Termination when the whole sample gets too thin
term_cb = create_thickness_termination_callback(cb_state, 1e-6)  # min_thickness [m]
```

`create_depletion_callback` accepts either a bare `StateLayout` or a
`DepletionCallbackState` as its fourth argument; passing the shared
`DepletionCallbackState` lets the termination callback see the resized layout
after each merge. The `material` argument is what makes the temperature merge
energy-conserving (omit it and the merge falls back to a volume-weighted
average — see §10.2.3).

### 10.2.3 The merge is mass- and energy-conserving

Concentrations are merged volume-weighted (exactly mass-conserving):

```
ξ_merged[j] = (V_remove · ξ_remove[j] + V_keep · ξ_keep[j]) / (V_remove + V_keep)
```

Temperature is merged by **conserving the condensed-phase sensible enthalpy**,
not by naive volume averaging. With `m_j = V · ξ_j` the per-component condensed
mass of each cell, the merged temperature solves

```
Σ_{j ∉ gas} m_j ∫_{T_ref}^{T_merged} c_p,j(T′) dT′  =  H_remove + H_keep
```

(a bracketed Newton solve — the root always lies between the two cell
temperatures). The integral form is exact for any temperature-dependent
`c_p,j(T)`, and gas heat capacity is excluded because Formulation A stores
thermal energy in the condensed matrix only. For constant `c_p` this reduces to
the familiar `ρc_p`-weighted mean.

This matters precisely in the depletion regime: merging a nearly-depleted
(low-`ρc_p`) surface cell with a dense interior cell can introduce a 10–20 %
temperature error if you average naively (e.g. ~120 K error when merging a 1 m³
solid at 800 K with a 1 m³ gas-filled cell at 400 K). The dry-solid mass used as
the pyrolysis-progress denominator (`m_dry_solid_initial`) is summed across the
merged cells, so progress accounting stays consistent. See Technical Reference
§14.3 for the derivation and §10 for the χ progress variable.

The merge prefers the **interior** neighbor (toward the substrate) for thermal
stability; this is `find_merge_neighbor` from `Pyrolysis.Adaptivity`, used
internally by the callback.

### 10.2.4 Inspecting what happened

After a solve that used depletion handling, query the merge history with the
exported helper (note the slightly fragile internals — it returns sentinel
values rather than throwing if it cannot reach the callback state):

```julia
stats = get_depletion_callback_stats(dep_cb)
# NamedTuple:
#   stats.n_merges     :: Int               number of cell merges performed
#   stats.merge_times  :: Vector{Float64}   simulation times [s] of each merge
#   stats.merged_cells :: Vector{Int}       original indices of merged cells
```

In ALE mode the solution arrays carry their own evidence of merges: the cell
count shrinks over time and removed-cell entries are NaN-padded in
`sol.T`, `sol.ξ`, and `sol.z` (see Chapter 9). `sol.extras.thickness_history`
records total sample thickness per output time.

### 10.2.5 Critical caveat: depletion is incompatible with a pluggable Jacobian

Depletion resizes the ODE state vector mid-solve via `resize!(integrator, …)`.
The structured/dense Jacobian backends carry an **immutable cache built once at
problem setup** that cannot survive a resize. Consequently:

- Passing **both** `jacobian = …` and `handle_depletion = true` to `solve`
  raises an error at validation time.
- With `handle_depletion = true` and the default `jacobian = nothing`, `solve`
  supplies **no** analytic Jacobian: the integrator falls back to its own
  **colored finite-difference** Jacobian. Because that dense FD Jacobian is
  incompatible with the sparse default linear solver, `solve` also switches the
  linear solver to a dense `LUFactorization`.

In short: **you cannot combine depletion handling with `Structured(...)` /
`DenseAD_ForwardDiff()` backends.** If you need the structured backend (e.g. for
a large ALE problem), run without depletion and rely on thickness termination
instead, or accept the FD fallback. See §10.8 and Technical Reference §14.7 and
§15.13.

> **MKL note.** The dense `LUFactorization` that `solve` selects under depletion
> is Julia's native LU and survives `resize!`. Do not pass a sparse
> factorization (`SparspakFactorization()`, KLU) together with depletion — there
> is no sparse Jacobian prototype in this mode, and `solve` will warn and swap
> it for dense LU. The MKL `ipiv` resize bug can only bite if you construct your
> own `integrator` without an explicit `linsolve`; the solver warns if it
> detects that combination on a large problem. See Chapter 11 and
> Technical Reference §15.12.

---

## 10.3 Thickness-termination callback

Independently of merging, a shrinking sample can reach (near-)zero thickness,
which makes the volumetric heat capacity singular. In **ALE mode `solve`
automatically installs a thickness-termination callback** that stops the
integration cleanly once the total sample thickness drops below
`min_thickness` (default `1e-6` m):

```julia
sol = solve(problem; use_ale = true, min_thickness = 1e-5)  # terminate at 10 µm
@show sol.retcode    # :Terminated  (vs :Success when the full tspan is reached)
```

The callback measures `thickness = z[end] - z[1]` from the node-position block of
the state vector, so it is meaningful only in ALE mode (`use_ale = true`); in
Eulerian mode the mesh does not move and the callback is not installed. To build
it directly:

```julia
term_cb = create_thickness_termination_callback(layout, Float64(min_thickness))
```

Pass a `DepletionCallbackState` instead of a bare `layout` when you also use
depletion handling, so the callback tracks the resized layout after each merge.

---

## 10.4 Error indicators (driving refinement decisions)

Error indicators estimate where the discretization error is largest. They are the
inputs to refinement/relocation decisions: high indicator value ⇒ resolve more
there. All indicator types live in `Pyrolysis.Adaptivity` and are re-exported, so
`using Pyrolysis` makes them available directly. Each implements
`compute_indicator(ind, mesh)` (returning a per-cell `Vector{Float64}`, with `0`
for inactive cells).

| Indicator | Constructor | Formula (η for cell *i*) | Use |
|---|---|---|---|
| Gradient | `GradientIndicator(:T)` or `GradientIndicator(:T, weight)` | `η = w · hᵢ · \|∇φ\|ᵢ` | Steep fronts (thermal/reaction). Default `weight = 1.0`. |
| Curvature | `CurvatureIndicator(:T)` or `(:T, weight)` | `η = w · hᵢ² · \|∇²φ\|ᵢ` | Sharper features than gradient; zero at boundaries (needs 3 points). |
| Reaction rate | `ReactionRateIndicator()` or `(weight)` | `η = w · hᵢ · maxⱼ \|rⱼ\|ᵢ / max(ξⱼ, ε)` | Active reaction zones. Requires `rates = r` (throws `ArgumentError` without it); standalone-only — it cannot run inside the h-AMR controller cycle. |
| Interface proximity | `InterfaceProximityIndicator(pos, h_target)` or `(pos, h_target, weight)` | `η = w / (1 + dᵢ/h_target)` | Force resolution near a tracked interface position, regardless of solution. Purely geometric. |
| Composite | `CompositeIndicator(inds)` or `(inds, combination)` | `:max`, `:l2`, or `:sum` of components, each first normalized by its own maximum over active cells | Blend several indicators. Default `combination = :max`. |

`φ` is the tracked field: `:T` for temperature, or a component index for a
concentration. `hᵢ` is the cell size and `dᵢ` the distance from the cell center
to the interface. The reaction-rate indicator's `ε = 1e-10` guards the division.

```julia
using Pyrolysis

# A single indicator
η_T = compute_indicator(GradientIndicator(:T), mesh)        # Vector{Float64}, length n_cells

# Blend a temperature gradient with temperature curvature, take the cell-wise max
ind = CompositeIndicator(
    AbstractErrorIndicator[GradientIndicator(:T, 1.0), CurvatureIndicator(:T, 0.5)],
    :max,
)
η = compute_indicator(ind, mesh)

# Reaction-rate indicator: standalone only, with explicit rates (sized like ξ)
η_r = compute_indicator(ReactionRateIndicator(), mesh, T, ξ; rates = r)
```

> **Note.** `compute_indicator(ind, mesh)` reads `T` and `ξ` from the mesh's cell
> states, so make sure the mesh is in sync with the current solution (after a
> solve, the step-update callback keeps `mesh` consistent; if you are
> post-processing a `PyrolysisSolution`, rebuild the field arrays yourself and
> use the four-argument form `compute_indicator(ind, mesh, T, ξ)`). The
> reaction-rate indicator requires a rate matrix that is not available from the
> mesh alone and **throws** in the convenience form — pass `rates = r` (sized
> like `ξ` on the current mesh) via the four-argument form.

These indicators are **diagnostic building blocks** — computing one does not
change the mesh. To act on an indicator you either feed it into the experimental
h-AMR machinery (§10.7) or, more commonly, use a *monitor function* for
r-refinement (§10.5).

---

## 10.5 r-refinement: node relocation by equidistribution

r-refinement keeps the **same number of cells** but slides interior nodes so they
cluster where the solution is hard to resolve. It is driven by a **monitor
function** `ω(z)`; nodes are placed so that each cell carries roughly equal
"monitor integral" (the equidistribution principle, Technical Reference §14.2).
Boundary nodes never move.

### 10.5.1 Monitor functions

| Monitor | Constructor | Formula | Notes |
|---|---|---|---|
| Gradient | `GradientMonitor()` or `GradientMonitor(; α=1.0, variable=:T)` | `ω = √(1 + α·\|∇φ\|²)` | Clusters nodes at steep gradients. Defaults `α = 1.0`, `variable = :T`. |
| Interface | `InterfaceMonitor(z)` or `InterfaceMonitor(z; σ=0.001, ω_base=1.0, ω_peak=10.0)` | `ω = ω_base + ω_peak·exp(−(z−z_interface)²/(2σ²))` | Gaussian peak at interface position `z`. |
| Composite | `CompositeMonitor(monitors, combination)` | `:max`, `:sum`, or `:product` | Combine monitors. |

### 10.5.2 Relocating nodes

```julia
using Pyrolysis

n_iter = relocate_nodes!(mesh;
    monitor    = GradientMonitor(; α = 1.0, variable = :T),
    max_iter   = 10,     # max equidistribution sweeps
    relaxation = 0.5,    # step fraction toward target: z_new = z + β·(z_target − z), β ∈ (0,1]
    tol        = 1e-3,   # converged when max relative node movement < tol
    h_min      = 1e-8,   # [m] minimum cell size enforced during moves
)
```

Defaults (verified in `src/Adaptivity/r_refinement.jl`): `monitor =
GradientMonitor()`, `max_iter = 10`, `relaxation = 0.5`, `tol = 1e-3`,
`h_min = 1e-8`. The function returns the number of iterations performed and
mutates `mesh` in place (updating node positions and re-deriving geometry).

What the algorithm does each sweep: evaluate the monitor, compute target node
positions by equidistributing the monitor integral, then take a relaxed step
toward each target while enforcing `z_{n+1} − z_n ≥ h_min`. Because the
minimum-spacing constraint can prevent nodes from reaching their targets, the
returned iteration count may equal `max_iter` without full convergence.

> **Caveats.** Monitor functions are problem-dependent and have arbitrary units;
> you are responsible for scaling so the equidistribution solve is
> well-conditioned (huge `ω` disparities cause ill-conditioning). Default
> tolerances (`tol = 1e-3`, `relaxation = 0.5`) are untested on realistic
> steep-front pyrolysis problems — treat them as starting points. Boundary nodes
> are pinned. r-refinement is a **standalone mesh operation**: it is not invoked
> automatically by `solve`. Use it between solve segments (relocate, then warm-
> start the next segment) or in a custom driver.

---

## 10.6 Mesh-quality metrics and enforcement

After any mesh motion (ALE or r-refinement) you may want to confirm the mesh is
still well-conditioned, and optionally repair it. The quality types and functions
are in `Pyrolysis.Adaptivity` (re-exported).

### 10.6.1 Metrics

| Metric | Constructor (default field) | Constraint per adjacent pair / cell |
|---|---|---|
| Aspect ratio | `AspectRatioMetric()` → `max_ratio = 2.0` | `hᵢ/hⱼ ≤ max_ratio` |
| Smoothness | `SmoothnessMetric()` → `max_size_jump = 0.5` | `\|hᵢ − hⱼ\| / max(hᵢ, hⱼ) ≤ max_size_jump` |
| Jacobian (positivity) | `JacobianMetric()` → `min_jacobian = 0.0` | `Vᵢ > min_jacobian` (signed width — inverted or collapsed cells fail) |
| Minimum cell size | `MinCellSizeMetric()` → `h_min = 1e-8` | `hᵢ ≥ h_min` |
| Composite | `CompositeQualityMetric(metrics)` | All component metrics must pass |

(Defaults verified in `src/Adaptivity/mesh_quality.jl`.)

### 10.6.2 Checking and summarizing

```julia
using Pyrolysis

metric = AspectRatioMetric(2.0)

check_quality(mesh, metric)        # Bool: does the whole mesh pass?
quality_violations(mesh, metric)   # Vector{Int}: cell indices that violate it

compute_aspect_ratios(mesh)        # Vector{Float64}: per-cell aspect ratio
compute_smoothness(mesh)           # Vector{Float64}: per-cell smoothness

summary = mesh_quality_summary(mesh)
# NamedTuple with:
#   n_active_cells, min_cell_size, max_cell_size,
#   max_aspect_ratio, max_smoothness, all_positive_volume
```

### 10.6.3 Enforcing quality

```julia
remaining = enforce_quality!(mesh, metric; max_iter = 10)
# Returns the number of remaining violations after Laplacian smoothing.
```

`enforce_quality!` applies **Laplacian (neighbor-averaging) smoothing** to the
interior nodes of violating cells, repeating until either no violations remain or
`max_iter` is hit, and returns the count of remaining violations (zero means the
mesh now passes). It is a local, incremental fix and is **not guaranteed to clear
clustered or severe violations** — check the return value. Boundary nodes are
not moved.

A typical post-motion guard:

```julia
relocate_nodes!(mesh; monitor = GradientMonitor(; variable = :T))
metric = CompositeQualityMetric(AbstractQualityMetric[
    AspectRatioMetric(2.0), MinCellSizeMetric(1e-7),
])
if !check_quality(mesh, metric)
    left = enforce_quality!(mesh, metric; max_iter = 20)
    left == 0 || @warn "mesh still has $left quality violation(s) after smoothing"
end
```

---

## 10.7 Experimental h-AMR (cell splitting and coarsening)

True topology-changing h-adaptive refinement — splitting cells into children and
merging sibling cells back — ships with the package but is **quarantined in
`Pyrolysis.Experimental` and is *not* connected to `solve`.** Unlike the
`Adaptivity` exports, `Experimental` is **not re-exported at the top level**, so
you must reach into it explicitly:

```julia
using Pyrolysis

# Single-cell split: returns (child_a_id, child_b_id)
ca, cb = Pyrolysis.Experimental.refine_cell!(mesh, cell_id)

# Threshold-based marking + batch refinement with 2:1 level balance
marked = Pyrolysis.Experimental.mark_cells_for_refinement(
    mesh, η, θ_refine, h_min, max_level)
marked = Pyrolysis.Experimental.enforce_2to1_balance!(mesh, marked, max_level)
Pyrolysis.Experimental.perform_h_refinement!(mesh, marked)

# Sibling coarsening (the inverse, for cells created by a previous refinement)
Pyrolysis.Experimental.coarsen_cells!(mesh, cell_a_id, cell_b_id; material = material)

# Orchestrator
ctrl = Pyrolysis.Experimental.AMRController(;
    error_indicator = GradientIndicator(:T),
    quality_metric  = AspectRatioMetric(),
    θ_refine    = 0.7,    # refine if η > θ_refine · max(η)
    θ_coarsen   = 0.1,    # coarsen if η < θ_coarsen · max(η)
    h_min       = 1e-6,
    h_max       = 1e-2,
    max_level   = 4,
    check_every = 10,
)
result = Pyrolysis.Experimental.adapt_mesh!(mesh, ctrl)   # (n_refined, n_coarsened, quality_ok)
```

The full experimental surface (all under `Pyrolysis.Experimental`) includes
`refine_cell!`, `perform_h_refinement!`, `mark_cells_for_refinement`,
`enforce_2to1_balance!`, `verify_refinement_conservation`, `coarsen_cells!`,
`perform_h_coarsening!`, `mark_cells_for_coarsening`, `can_coarsen`,
`get_sibling`, `verify_coarsening_conservation`, `hr_adapt!` (combined h-then-r),
`AMRController`, `adapt_mesh!`, `should_adapt`, `maybe_adapt!`,
`amr_statistics`, `reset_statistics!`, and `create_default_controller`.

**Caveats — read before using:**

- **Not in the solve loop.** There is no `solve` keyword to enable h-AMR during
  integration. Calling these between solve segments is the only supported pattern,
  and even then the state-vector resize logic that depletion handling uses is not
  shared here.
- **Bookkeeping is leaky.** h-refinement appends children and h-coarsening marks
  cells inactive *without compacting the arrays*, so `length(mesh.cells) >
  n_active_cells` is normal. Any code consuming the mesh afterward must iterate
  the `active_cell_mask`, not raw index ranges. `assert_amr_invariants(mesh)` (a
  core `Adaptivity` export) is called after each operation to catch drift.
- **Conservation is exact, but the temperature split is zeroth-order.** Children
  inherit the parent concentration with equal volume split (mass-exact); the
  temperature split is piecewise-constant — both children inherit the parent
  temperature (the coded gradient correction evaluates to zero because the split
  point coincides with the parent cell center) — exactly energy-conserving but
  adding no sub-cell resolution. Use
  `verify_refinement_conservation` / `verify_coarsening_conservation` to check.
- **Defaults are untuned.** Controller thresholds (`θ_refine = 0.7`,
  `θ_coarsen = 0.1`, `max_level = 4`, `check_every = 10`) have not been validated
  on realistic pyrolysis fronts.

For production surface-recession work, prefer the depletion path (§10.2), which
*is* integrated with the solver and uses the energy-conserving adjacent-cell
merge (`merge_adjacent_cells!`, `merge_surface_cell!`) that lives in core
`Adaptivity`. See Technical Reference §14.4 for the h-AMR algorithms and §14.3
for the production merge.

---

## 10.8 Choosing the Jacobian when you use adaptivity

The single most common configuration mistake with adaptivity is combining
depletion handling with an explicit Jacobian backend. The rule:

| Scenario | Jacobian to use |
|---|---|
| ALE, **no** depletion handling | Default `jacobian = nothing` ⇒ `Structured(strategy = DirectSolve())`; or opt into `ApproxSolve` for large ALE problems (Chapter 11). |
| ALE **with** `handle_depletion = true` | **Leave `jacobian = nothing`.** `solve` supplies no analytic Jacobian and the integrator uses colored FD with a dense `LUFactorization`. Do **not** pass a `JacobianSpec` — it errors. |
| Verifying a backend during an adaptive run | `solve(...; verify_jacobian = (rtol=1e-6, atol=1e-10, sample_times=[…]))` — but only meaningful when a backend is active (i.e. not with depletion). |

Concretely, this is correct:

```julia
# Large ALE problem WITHOUT depletion — opt into the cheaper approximate solve
sol = solve(problem;
    use_ale = true,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
)
```

and this **errors** at validation:

```julia
# WRONG: pluggable Jacobian + depletion are mutually exclusive
solve(problem;
    use_ale = true,
    handle_depletion = true,
    jacobian = JacobianSpec(Structured(strategy = DirectSolve()), problem),  # ← rejected
)
```

If you need both the analytic structured Jacobian *and* surface recession, the
practical recourse is to forgo per-cell merging and rely on the
thickness-termination callback alone (run until the sample is "thin enough,"
then stop), keeping the state-vector length fixed so the structured backend's
cache stays valid. See Chapter 11 for the full Jacobian/linear-solver decision
guide and Technical Reference §15.10–§15.13 for the strategy internals.

---

## 10.9 Summary of options and defaults

| Option | Where | Default | Section |
|---|---|---|---|
| `use_ale` | `solve` kwarg | `false` | 10.2 |
| `handle_depletion` | `solve` kwarg | `false` | 10.2 |
| `depletion_threshold` | `solve` kwarg | `0.05` | 10.2 |
| `min_cells` | `solve` kwarg | `2` | 10.2 |
| `min_thickness` | `solve` kwarg | `1e-6` m | 10.3 |
| `GradientIndicator(:T, w)` | `Adaptivity` | `w = 1.0` | 10.4 |
| `CompositeIndicator(inds, c)` | `Adaptivity` | `c = :max` | 10.4 |
| `GradientMonitor(; α, variable)` | `Adaptivity` | `α=1.0, variable=:T` | 10.5 |
| `InterfaceMonitor(z; σ, ω_base, ω_peak)` | `Adaptivity` | `σ=1e-3, ω_base=1.0, ω_peak=10.0` | 10.5 |
| `relocate_nodes!(; max_iter, relaxation, tol, h_min)` | `Adaptivity` | `10, 0.5, 1e-3, 1e-8` | 10.5 |
| `AspectRatioMetric(max_ratio)` | `Adaptivity` | `2.0` | 10.6 |
| `SmoothnessMetric(max_size_jump)` | `Adaptivity` | `0.5` | 10.6 |
| `MinCellSizeMetric(h_min)` | `Adaptivity` | `1e-8` m | 10.6 |
| `enforce_quality!(; max_iter)` | `Adaptivity` | `10` | 10.6 |
| h-AMR (`refine_cell!`, `AMRController`, …) | `Experimental` (not re-exported) | — experimental | 10.7 |

**Key correctness reminders.** Depletion handling requires `use_ale = true`; it
cannot be combined with a pluggable `jacobian` backend (it falls back to colored
FD with dense LU); the surface merge conserves species mass exactly and condensed-phase sensible
enthalpy via the `∫c_p dT` temperature merge (gas excluded); r-refinement,
quality enforcement, and h-AMR
are standalone tools not driven by `solve`; and the h-AMR cluster is experimental
and reached only through `Pyrolysis.Experimental`. For the physics behind these
operations, see Technical Reference §11 (ALE), §14 (AMR), and §15 (Jacobian).
