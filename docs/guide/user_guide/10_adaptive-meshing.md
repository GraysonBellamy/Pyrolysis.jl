# 10. Adaptive Meshing Options

The only mesh adaptation `solve` performs is **surface-cell depletion
handling**: as the solid gasifies, a depleted cell is *parked* — merged into its
neighbour and its state slot retired above the live surface — on a fixed mesh or
on the shrinking Arbitrary Lagrangian–Eulerian (ALE) mesh. Beyond that, the
`Pyrolysis.Adaptivity` submodule exposes standalone building blocks that `solve`
does not use: **error indicators**, **monitor functions** and **r-refinement**
(moving nodes to concentrate them where a field is steep), and **mesh-quality
metrics**. Cell splitting and sibling coarsening (h-refinement) ship in
`Pyrolysis.Experimental`, and the meshes they produce cannot be passed to
`solve`. This chapter shows how to use each, what every option means and defaults
to, and where the sharp edges are. For the underlying theory (equidistribution,
energy-conserving merges, the geometric conservation law) see Technical Reference
§14, §11 (ALE), and §13 (finite-volume discretization).

> **Coordinate reminder.** `z = 0` is the bottom / substrate (boundary tag
> `:bottom`, id 2); `z = L` is the top / exposed surface (tag `:top`, id 1). Heat
> enters from the top and the material recedes downward, so the "surface cell" is
> the last (highest-`z`) active cell. See Chapter 6 (Meshes and Geometry).

---

## 10.1 When to use adaptivity

Reach for depletion handling only when the physics warrants it, and size the
mesh up front for everything else:

| Situation | Symptom | Recommended tool |
|---|---|---|
| Surface recession or burnout (charring/ablating solid) | Surface cell thins or empties toward zero; timestep collapses; `ρc_p → 0` trouble | `handle_depletion = true` (§10.2) |
| Sample fully consumed | Total thickness approaches zero | thickness termination (automatic in ALE; §10.3) |
| Steep thermal or reaction front | Front under-resolved between cells | a surface-stretched mesh or more cells (Chapter 6) |
| Grading or distortion of a mesh you built | Large cell-size jumps, near-zero volumes | quality metrics (§10.6) to inspect it before solving |

A good default workflow is: start with a **uniform or surface-stretched fixed
mesh** (Chapter 6) sized from a convergence study (Technical Reference §18.5),
and turn on depletion handling once you run a case where the surface actually
recedes or burns out. The indicators, r-refinement and quality tools of
§10.4–§10.6 are building blocks for your own mesh studies; `solve` does not call
them, and most cone-calorimeter and TGA studies never need them.

---

## 10.2 Depletion handling: parked cells (the production path)

When a solid pyrolyzes away, the cells that held it shrink (ALE mesh motion,
`use_ale = true`) or empty (fixed mesh). Either way the integrator ends up
with cells it cannot resolve: a thin cell's equations are conditioned like
`1/Δz²`, so the Newton iteration loses every digit and the step size collapses.
**Depletion handling** detects such cells and *parks* them: the cell's mass and
sensible enthalpy are merged into a neighbour, and its state slot is moved
above the live surface as an inert ghost. Nothing changes size — the state
vector, the mesh arrays and the Jacobian keep their allocation for the whole
solve — and every kernel runs over the live cells `1:mesh.n_active_cells` with
the boundary conditions applied at the live surface.

### 10.2.1 Turning it on through `solve`

Pass `handle_depletion = true` to `solve`, in either mesh mode:

```julia
using Pyrolysis

sol = solve(problem;
    use_ale             = true,    # node positions become state; mesh can shrink
    handle_depletion    = true,    # park depleted cells instead of stalling
    depletion_threshold = 0.05,    # park a cell once it drops below 5% of its
                                   #   initial thickness (default; must be in (0,1))
    min_cells           = 2,       # never park below this many live cells; terminate instead
    min_thickness       = 1e-6,    # [m] terminate when the live thickness < this (default)
    radiation_model     = SURFACE_ABSORPTION,
)
```

| Keyword | Default | Meaning |
|---|---|---|
| `use_ale` | `false` | Enable ALE mesh motion (node positions enter the state vector). |
| `handle_depletion` | `false` | Park depleted cells mid-solve. |
| `depletion_threshold` | `0.05` | ALE: a cell is depleted when its thickness `< threshold ×` its own initial thickness. Fixed mesh: the top live cell is depleted when its condensed concentration `< threshold ×` its initial value. Must satisfy `0 < threshold < 1`. |
| `min_cells` | `2` | Once the live cell count reaches this, the run terminates rather than parking further. Must be `≥ 1`. |
| `min_thickness` | `1e-6` | [m] Live thickness below which an ALE solve terminates. Must be `> 0`. |

Leave `depletion_threshold` at its default unless a case shows otherwise. A
larger value merges the hottest cell into a colder one earlier and disturbs the
surface temperature; on a fast ablation front it does not remove the Newton
failures either (Technical Reference §14.3.2).

Choose the mesh mode by the material, not by cost. A fixed mesh keeps every
cell's volume, so a gasifying layer thins in *density*, and only the top live
cell can park. The layer burns top down when the surface cell is the hottest,
as under surface absorption. Under in-depth absorption the surface cell of a
regressing polymer is the coldest: it stays in place as a porous skin while the
cells beneath it empty, and the model carries no heat across the void they
leave. A fixed-mesh solve therefore stops at the first cell below the surface
that becomes a void (matrix volume fraction below 10⁻⁴), warns that the
material has gasified from within, and returns retcode `:Terminated`
(Technical Reference §14.3.2). A fixed mesh with depletion suits a char that
keeps its volume and a layer that burns top down; use `use_ale = true` for
materials that regress.

Each callback invocation parks **at most one cell** so the integrator can
re-adapt its step between events; candidates are processed surface-first.
After a park the integrator restarts from a step of `DT_AFTER_PARK`
(10⁻⁶ s): the merged cell is off the quasi-steady profiles of its neighbours
and the pre-park step would fail its Newton repeatedly. The step-size
controller grows the step back within a few steps, so a park costs a handful
of extra steps and no failures; expect the step history to show a spike to
10⁻⁶ s at every entry of `sol.extras.park_times`.

### 10.2.2 What `solve` wires up under the hood

With `handle_depletion = true`, `solve` assembles a `CallbackSet` containing
(among others):

1. a **step-update callback** that mirrors the integrator's accepted state back
   into the mutable mesh (so subsequent callbacks see current geometry);
2. a **thickness-termination callback** (ALE, §10.3) on the live thickness;
3. a **depletion callback** that finds the next depleted cell and parks it.

You can build these yourself if you are driving the integrator directly:

```julia
using Pyrolysis.Solver: DepletionCallbackState, create_depletion_callback,
                        create_thickness_termination_callback

# Shared bookkeeping: live count, per-slot depletion reference, park log
cb_state = DepletionCallbackState(layout, mesh, material, u0)

dep_cb  = create_depletion_callback(ws, mesh, material, cb_state;
    threshold = 0.05,   # same meaning as depletion_threshold above
    min_cells = 2)

term_cb = create_thickness_termination_callback(cb_state, 1e-6)  # live thickness [m]
```

### 10.2.3 The park conserves mass and energy

Concentrations are merged volume-weighted (exactly mass-conserving):

```
ξ_keep[j] = (V_park · ξ_park[j] + V_keep · ξ_keep[j]) / V_keep′
```

with `V_keep′ = V_park + V_keep` under ALE (the keeper spans both cells) and
`V_keep′ = V_keep` on a fixed mesh (the keeper's volume is fixed, so its
concentration rises by the parked mass).

Temperature is merged by **conserving the sensible enthalpy**, not by naive
volume averaging. With `m_j = V · ξ_j` the per-component mass of each cell, the
merged temperature solves

```
Σ_j m_j ∫_{T_ref}^{T_merged} c_p,j(T′) dT′  =  H_park + H_keep
```

(a bracketed Newton solve — the root always lies between the two cell
temperatures). The integral form is exact for any temperature-dependent
`c_p,j(T)`, and every component is included, gas and condensed alike — the
same `ρc_p` the temperature equation stores. For constant `c_p` this reduces to
the familiar `ρc_p`-weighted mean.

This matters precisely in the depletion regime: merging a nearly-depleted
(low-`ρc_p`) cell with a dense interior cell can introduce a 10–20 %
temperature error if you average naively (e.g. ~120 K error when merging a 1 m³
solid at 800 K with a 1 m³ gas-filled cell at 400 K). The dry-solid mass used as
the pyrolysis-progress denominator (`m_dry_solid_initial`) is summed across the
merged cells and the merged `χ` is weighted by it, so progress accounting stays
consistent. See Technical Reference §14.3 for the derivation and §10 for the χ
progress variable.

The park prefers the **interior** neighbour (toward the substrate) for thermal
stability; this is `find_merge_neighbor` from `Pyrolysis.Adaptivity`, used
internally by the callback.

### 10.2.4 Inspecting what happened

Ghost slots read as `NaN` in `sol.T`, `sol.ξ` and (ALE) `sol.z`;
`sol.extras.thickness_history` is the live thickness;
`sol.extras.live_cells[k]` is the number of live cells at output time `k`, and
`sol.extras.park_times` lists the event times:

```julia
n_parks = length(sol.extras.park_times)
sol.extras.live_cells[end]          # live cells at the end of the run
```

### 10.2.5 Depletion and the Jacobian

Parking keeps the state size, so the structured Jacobian and its sparse linear
solver run unchanged through burnout; `verify_jacobian` covers
parked states. Under ALE the default backend is
`Structured(strategy = ApproxSolve())`, with or without parking: the full ALE
pattern couples every node to every cell below it and makes the factorized
matrix effectively dense, while the couplings the reduced pattern drops do not
change the trajectory. Pass `jacobian = Structured(strategy = DirectSolve())`
to keep the exact assembly.

A park is a discrete rewrite of the state: parameter sensitivities are valid on
horizons without park events only (Chapter 13).

---

## 10.3 Thickness-termination callback

Independently of parking, a shrinking sample can reach (near-)zero thickness,
which makes the volumetric heat capacity singular. In **ALE mode `solve`
automatically installs a thickness-termination callback** that stops the
integration cleanly once the total sample thickness drops below
`min_thickness` (default `1e-6` m):

```julia
sol = solve(problem; use_ale = true, min_thickness = 1e-5)  # terminate at 10 µm
@show sol.retcode    # :Terminated  (vs :Success when the full tspan is reached)
```

The callback measures the live thickness `z[n_live+1] - z[1]` from the
node-position block of the state vector, so it is meaningful only in ALE mode
(`use_ale = true`); in Eulerian mode the mesh does not move and the callback is
not installed. To build it directly:

```julia
term_cb = create_thickness_termination_callback(layout, Float64(min_thickness))
```

Pass a `DepletionCallbackState` instead of a bare `layout` when you also use
depletion handling, so the callback reads the live cell count.

---

## 10.4 Error indicators (driving refinement decisions)

Error indicators estimate where the discretization error is largest. They are the
inputs to refinement/relocation decisions: high indicator value ⇒ resolve more
there. All indicator types are exported from `Pyrolysis.Adaptivity` (not from the
top level), so bring them in with `using Pyrolysis.Internal`. Each implements
`compute_indicator(ind, mesh)` (returning a per-cell `Vector{Float64}`, with `0`
for inactive cells).

| Indicator | Constructor | Formula (η for cell *i*) | Use |
|---|---|---|---|
| Gradient | `GradientIndicator(:T)` or `GradientIndicator(:T, weight)` | `η = w · hᵢ · \|∇φ\|ᵢ` | Steep fronts (thermal/reaction). Default `weight = 1.0`. |
| Curvature | `CurvatureIndicator(:T)` or `(:T, weight)` | `η = w · hᵢ² · \|∇²φ\|ᵢ` | Sharper features than gradient; zero at boundaries (needs 3 points). |
| Reaction rate | `ReactionRateIndicator()` or `(weight)` | `η = w · hᵢ · maxⱼ \|rⱼ\|ᵢ / max(ξⱼ, ε)` | Active reaction zones. Requires `rates = r` (throws `ArgumentError` without it). |
| Interface proximity | `InterfaceProximityIndicator(pos, h_target)` or `(pos, h_target, weight)` | `η = w / (1 + dᵢ/h_target)` | Force resolution near a tracked interface position, regardless of solution. Purely geometric. |
| Composite | `CompositeIndicator(inds)` or `(inds, combination)` | `:max`, `:l2`, or `:sum` of components, each first normalized by its own maximum over active cells | Blend several indicators. Default `combination = :max`. |

`φ` is the tracked field: `:T` for temperature, or a component index for a
concentration. `hᵢ` is the cell size and `dᵢ` the distance from the cell center
to the interface. The reaction-rate indicator's `ε = 1e-10` guards the division.

```julia
using Pyrolysis, Pyrolysis.Internal

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
> states, so make sure the mesh holds the state you mean (after a solve,
> `solution.problem.mesh` holds the final state; to evaluate a saved time of a
> `PyrolysisSolution`, rebuild the field arrays yourself and use the
> four-argument form `compute_indicator(ind, mesh, T, ξ)`). The
> reaction-rate indicator requires a rate matrix that is not available from the
> mesh alone and **throws** in the convenience form — pass `rates = r` (sized
> like `ξ` on the current mesh) via the four-argument form.

These indicators are **diagnostic building blocks** — computing one does not
change the mesh, and `solve` never computes one.

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
using Pyrolysis, Pyrolysis.Internal

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
> are pinned. r-refinement is a **standalone mesh operation** that `solve` does
> not call: it moves the nodes of a `Mesh1D` you hold, and the cell values stay
> attached to their cells (they are not remapped).

---

## 10.6 Mesh-quality metrics and enforcement

To confirm that a mesh you built or relocated is well-conditioned, and optionally
repair it, use the quality types and functions exported from
`Pyrolysis.Adaptivity` (bring them in with `using Pyrolysis.Internal`). `solve`
does not run them.

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
using Pyrolysis, Pyrolysis.Internal

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

A typical guard after relocating nodes:

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

Topology-changing h-adaptive refinement — splitting cells into children and
merging sibling cells back — ships with the package in `Pyrolysis.Experimental`,
and **its meshes cannot be passed to `solve`**. `Experimental` is not part of
`Pyrolysis.Internal`, so you reach into it explicitly:

```julia
using Pyrolysis, Pyrolysis.Internal   # GradientIndicator, AspectRatioMetric

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

- **Not usable with `solve`.** Refinement leaves inactive parent cells among the
  active ones, while the solver's kernels assume the active cells come first, so
  a refined or coarsened mesh cannot be solved. There is no `solve` keyword for
  h-refinement either.
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

For surface-recession work use the depletion path (§10.2), which is integrated
with the solver and parks depleted cells with the energy-conserving merge
(`park_cell!`) that lives in core `Adaptivity`. See Technical Reference §14.4 for
the h-refinement algorithms and §14.3 for the park.

---

## 10.8 Choosing the Jacobian with depletion handling

Depletion handling keeps the state size, so either Jacobian strategy works with
it. The defaults:

| Scenario | Default Jacobian |
|---|---|
| ALE, with or without `handle_depletion` | `Structured(strategy = ApproxSolve())` — the reduced ALE pattern; it follows the same trajectory as `DirectSolve` at a fraction of the cost per step. Without parking, surface cells that burn away shrink toward zero thickness and stall either Jacobian alike: use `handle_depletion = true` to run through burnout. |
| Fixed mesh | `Structured(strategy = DirectSolve())` (the strategies coincide unless Beer–Lambert couplings are present). |
| Verifying a backend during a depletion run | `solve(...; verify_jacobian = (rtol=1e-6, atol=1e-10, sample_times=[…]))`, with or without depletion. |

Override with a bare backend, built for the solve's own mode:

```julia
sol = solve(problem;
    use_ale = true,
    handle_depletion = true,
    jacobian = Structured(strategy = DirectSolve()),   # exact assembly through burnout
)
```

See Chapter 11 for the full Jacobian/linear-solver decision guide and Technical
Reference §15.6–§15.10 for the strategy internals.

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

**Key correctness reminders.** Depletion handling works in both mesh modes
and with either Jacobian strategy; the park conserves species mass exactly and
sensible enthalpy (every component, gas included) via the `∫c_p dT`
temperature merge; error indicators, r-refinement and quality enforcement are
standalone tools `solve` does not use; and the h-refinement cluster is
experimental, reached only through `Pyrolysis.Experimental`, and its meshes
cannot be passed to `solve`. For the physics behind these operations, see
Technical Reference §11 (ALE), §14 (mesh adaptation), and §15 (Jacobian).
