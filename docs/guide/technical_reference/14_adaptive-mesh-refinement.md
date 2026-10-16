# 14. Adaptive Mesh Refinement and Surface Depletion

This chapter develops the adaptive-mesh-refinement (AMR) machinery of Pyrolysis.jl
and, in particular, the *surface-depletion* mechanism that is the only AMR path
exercised in the production solver. We treat the local error indicators that flag
where the mesh is under-resolved, the equidistribution theory behind
$r$-refinement (node relocation), the energy-conserving cell-merge algebra that
removes the $\rho c_p\!\to\!0$ singularity of a depleting surface cell, the
experimental $h$-refinement/$h$-coarsening (cell split/merge) topology operators,
the mesh-quality metrics that keep the discretization well-conditioned, the
active-cell bookkeeping that all AMR consumers must respect, and the way the
depletion callback parks a burned-out cell without changing the size of the
state. Every discrete rule, threshold, guard, and sign convention below is grounded
in the implementing function so that the mathematics is traceable to the code;
where the code makes an approximation, a guard, or a choice that differs from the
idealized form (including a temperature-interpolation term that vanishes
identically), we say so and trust the code.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$ (axial, through-thickness). **$z=0$ is the
  bottom/substrate** (boundary id 2, tag `:bottom`); **$z=L$ is the top/exposed
  surface** (boundary id 1, tag `:top`). Heat enters from the top; the material
  recedes downward. Cell indices increase from bottom to top, so the "left" cell
  of a face is at lower $z$ and the "right" cell is at higher $z$.
- The cell index is $i$ ($i=1,\dots,n$); $j$ is a component index
  ($j=1,\dots,N_C$); $r$ a reaction index; $f$ a face index. This follows
  overload resolution I1 of §2.
- Following the adaptivity subsystem we write the **cell size** $h_i\equiv\Delta z_i$
  for the through-thickness width of cell $i$, to match the AMR literature; this
  is the same quantity called $\Delta z_i$ elsewhere in the guide. We reserve
  $\eta_i$ for an error indicator value and $\omega(z)$ for a monitor function, as
  in §2.
- $\phi$ in this chapter is a **generic scalar field** (temperature $T$ or a
  concentration $\xi_j$) being tracked by an indicator or a monitor, *not* the
  porosity of §5. The porosity does not appear in this chapter.

The subsystem splits along a maturity line. The **core** `Adaptivity` module
(`src/Adaptivity/`) holds the depleted-cell park (`park_cell!`) with its
energy-conserving temperature merge and neighbor finders — the only adaptivity
`solve` uses — and building blocks that `solve` does not call: the error
indicators, $r$-refinement, the mesh-quality metrics, and the invariant checker.
The **experimental** $h$-AMR cluster (`src/Experimental/`) holds cell splitting
(`refine_cell!`), sibling coarsening (`coarsen_cells!`), the `AMRController`
orchestrator, and the combined `hr_adapt!` driver. It is not wired into `solve`,
and a mesh it refines cannot be passed to `solve`: refinement leaves inactive
parent cells among the active ones, and the solver's kernels assume the active
cells come first (`src/Experimental/Experimental.jl`). This chapter covers both,
flagging the status of each throughout.

---

## 14.1 Error indicators

An *error indicator* $\eta_i\ge 0$ assigns to each active cell a scalar estimate
of the local discretization error, used to decide where to refine (large $\eta$)
or coarsen (small $\eta$). Pyrolysis.jl implements four physically motivated
indicators plus a composite, all subtypes of the abstract type
`AbstractErrorIndicator` (declared in `src/Materials/Materials.jl`, so that the
type is visible to both the `Adaptivity` and `Experimental` modules). All
indicator math lives in `src/Adaptivity/error_indicators.jl`.

The unifying design is **dimensional consistency through the cell size**: an
$O(h^p)$ truncation-error term is mimicked by multiplying a discrete derivative of
the field by the appropriate power of $h_i$. This makes $\eta_i$ behave like a
local error contribution, so that thresholding $\eta_i$ against the global maximum
$\max_j\eta_j$ adapts the mesh toward equidistributed error.

### 14.1.1 Gradient indicator

The gradient indicator targets steep solution fronts (thermal waves, reaction
zones). For a tracked field $\phi$,

$$
\eta_i^{\mathrm{grad}} = w\, h_i \,\bigl|\nabla\phi\bigr|_i ,
$$

where $w$ is a user weight (default $1$), $h_i$ is the cell size, and
$|\nabla\phi|_i$ is the magnitude of the discrete first derivative at cell $i$.
The factor $h_i$ makes $\eta^{\mathrm{grad}}$ scale like the first-order term in a
Taylor expansion: a fixed jump across a cell produces a fixed indicator regardless
of how the cell is sized, so the indicator flags *resolution relative to the
feature*, not absolute slope.

The gradient itself (`compute_gradient_1d`) is the
opposite-distance-weighted mean of the two one-sided slopes
(`opposite_distance_gradient`, §13.6.1) — second order on non-uniform
meshes — and degrades gracefully at boundaries:

$$
(\nabla\phi)_i =
\begin{cases}
\dfrac{d_R\,\frac{\phi_i-\phi_L}{d_L} + d_L\,\frac{\phi_R-\phi_i}{d_R}}{d_L+d_R}, & \text{both neighbors present},\\[2.2ex]
\dfrac{\phi_R-\phi_i}{z^c_R-z^c_i}, & \text{left is a boundary (one-sided right)},\\[2.2ex]
\dfrac{\phi_i-\phi_L}{z^c_i-z^c_L}, & \text{right is a boundary (one-sided left)},\\[2.2ex]
0, & \text{a lone cell (no neighbors)} .
\end{cases}
$$

Here $L,R$ index the cells immediately below/above cell $i$, found through the
cell's left/right faces, $z^c$ are cell-center positions, and
$d_L = z^c_i - z^c_L$, $d_R = z^c_R - z^c_i$ are the center-to-center
distances. A boundary cell takes the one-sided difference toward its only
neighbor, so it never reads a missing value, and a zero center-to-center distance
returns a zero gradient.

### 14.1.2 Curvature indicator

The curvature indicator is *more* sensitive to sharp features than the gradient
indicator, because it carries two powers of the cell size:

$$
\eta_i^{\mathrm{curv}} = w\, h_i^2 \,\bigl|\nabla^2\phi\bigr|_i .
$$

The discrete second derivative (`compute_curvature_1d`) uses the
three-point **non-uniform mesh** stencil:

$$
(\nabla^2\phi)_i =
\frac{2}{h_L+h_R}
\left[
  \frac{\phi_R-\phi_i}{h_R}
  -
  \frac{\phi_i-\phi_L}{h_L}
\right],
\qquad
h_L=z^c_i-z^c_L,\quad h_R=z^c_R-z^c_i .
$$

For $h_L=h_R=h$, this reduces to the standard centered formula
$(\phi_R-2\phi_i+\phi_L)/h^2$. On graded meshes it preserves zero curvature for
linear fields, avoiding false refinement from mesh stretching alone.

Because three collinear points are required, the curvature is returned as
**exactly zero at any cell that lacks an interior neighbor on either side**
(boundary cells); the indicator therefore never drives refinement of the first or
last cell. This is a deliberate restriction: a one-sided second difference at a
boundary is ill-conditioned and would falsely flag the surface cell, which is
already handled by the depletion mechanism (§14.3).

> **Implementation note.** Gradient and curvature indicators are
> `compute_indicator(::GradientIndicator, …)` /
> `compute_indicator(::CurvatureIndicator, …)`, with the derivative kernels
> `compute_gradient_1d` and `compute_curvature_1d`, all in
> `src/Adaptivity/error_indicators.jl`. Both skip cells with
> `active_cell_mask[i] == false`, returning $\eta_i=0$ there. The variable to
> track is the `variable::Union{Symbol,Int}` field — `:T` for temperature, or an
> integer component index `j` to track $\xi_{j,\cdot}$; any other symbol throws
> (component names are not resolvable without a material in scope).

### 14.1.3 Reaction-rate indicator

Conduction- and diffusion-based indicators cannot see a reaction front that has
not yet imprinted a steep temperature or concentration gradient. The
reaction-rate indicator refines wherever chemistry is locally fast relative to the
available reactant:

$$
\eta_i^{\mathrm{react}} = w\, h_i \,\max_j \frac{|r_j|_i}{\max(\xi_{j,i},\,\varepsilon)},
\qquad \varepsilon = 10^{-10}\ \mathrm{kg\,m^{-3}} .
$$

The normalization by $\xi_{j,i}$ converts an absolute rate (kg m$^{-3}$ s$^{-1}$)
into an inverse time scale (s$^{-1}$), i.e. a *fractional* consumption rate; the
$\max_j$ takes the worst-case over the rate matrix rows, so a single fast species
suffices to flag the cell. The floor $\varepsilon$ prevents division by zero in a
depleted cell.

This indicator is special in that it **requires** the per-cell, per-species
reaction rate matrix `rates`, passed as a keyword:
`compute_indicator(ind, mesh, T, ξ; rates = r)` with `r` sized like ξ for the
*current* mesh. Calling it without rates throws an `ArgumentError`, and a `rates` matrix whose column count does not match
the mesh throws `DimensionMismatch`. Because `adapt_mesh!` recomputes
indicators *after* topology changes — when any precomputed rates are stale —
this indicator cannot run inside the controller cycle and is standalone-only;
the default controller uses a gradient + curvature composite instead.
Supplying live rates inside the time-stepping loop would mean recomputing the
Arrhenius source at the indicator-evaluation point (see §6, Reaction
Kinetics) — this coupling gap is an open issue (§14.10).

> **Implementation note.** `compute_indicator(::ReactionRateIndicator, mesh, T, ξ;
> rates=nothing)` in `src/Adaptivity/error_indicators.jl`.

### 14.1.4 Interface-proximity indicator

The previous three indicators are *solution-driven*. The interface-proximity
indicator is purely *geometric*: it forces refinement in a neighborhood of a
tracked interface position $z_{\mathrm{int}}$ regardless of the current solution,
which is useful for a moving pyrolysis front that the solution has not yet
resolved:

$$
\eta_i^{\mathrm{iface}} = \frac{w}{1 + d_i/h_{\mathrm{target}}},
\qquad
d_i = |z^c_i - z_{\mathrm{int}}| ,
$$

with $h_{\mathrm{target}}$ the desired cell size near the interface. The indicator
is a Lorentzian bump: it equals $w$ exactly at the interface ($d_i=0$), falls to
$w/2$ one target-size away, and decays as $w\,h_{\mathrm{target}}/d_i$ far from it.
Because $\eta^{\mathrm{iface}}$ is bounded and smooth, it composes cleanly with the
solution-driven indicators.

> **Implementation note.** `compute_indicator(::InterfaceProximityIndicator, …)`
> in `src/Adaptivity/error_indicators.jl`. The struct carries
> `interface_position`, `h_target`, and `weight`; no derivative is computed.

### 14.1.5 Composite indicator

Physical problems usually need several indicators at once (e.g. a thermal-gradient
indicator together with a curvature indicator). `CompositeIndicator` evaluates
each component indicator $\eta_k(i)$, **normalizes each by its own maximum
over active cells**, $\hat\eta_k(i) = \eta_k(i)/\max_{j\in\text{active}}
\eta_k(j)$, and combines the normalized values per-cell by one of three rules:

$$
\eta_i^{\max} = \max_k \hat\eta_k(i),
\qquad
\eta_i^{\ell_2} = \sqrt{\sum_k \hat\eta_k(i)^2},
\qquad
\eta_i^{\Sigma} = \sum_k \hat\eta_k(i) .
$$

The normalization is what makes the combination meaningful: component
indicators carry incommensurable dimensions (a temperature-gradient indicator
is in K, a reaction-rate indicator in s⁻¹), so combining raw values would let
whichever indicator has the largest unit-system numbers dominate. The $\max$
rule (the default) is the most conservative: a cell is flagged if *any*
component demands it. The $\ell_2$ and $\sum$ rules blend contributions,
useful when several weak indicators should reinforce one another. The
per-component weights $w_k$ are baked into each $\eta_k$ *before*
normalization, so they tilt the within-indicator shape, while the relative
scale across indicators is set by the normalization.

> **Implementation note.** `compute_indicator(::CompositeIndicator, …)` in
> `src/Adaptivity/error_indicators.jl` dispatches on the
> `combination::Symbol` field (`:max`, `:l2`, `:sum`) and errors on any other
> value. A convenience overload `compute_indicator(ind, mesh)` pulls $T$ and $\xi$ from the mesh cell states via `extract_fields` and forwards
> to the field-taking method; this dense extraction happens on every call and is
> not cached.

### 14.1.6 Indicator interface and the active-cell contract

Every indicator computes over `n_cells(mesh)` cells but **writes $\eta_i=0$ for any
inactive cell** (`active_cell_mask[i] == false`). This is the first instance of a
pervasive AMR contract: because $h$-refinement appends children to `mesh.cells`
without ever compacting (§14.6), `length(mesh.cells)` can exceed
`n_active_cells`, and *every* consumer must guard on the mask rather than iterate a
raw index range. The indicators honor this; downstream marking routines
(§14.4.2) compute $\max\eta$ only over masked cells.

---

## 14.2 $r$-refinement: node relocation by equidistribution

$r$-refinement ("relocation refinement") moves the interior nodes of a *fixed
topology* mesh to concentrate resolution where it is needed, without changing the
number of cells. Unlike $h$-refinement it adds no degrees of freedom — it is a
pure remeshing of the existing cells. `solve` does not call it: the ALE mesh moves
with the solid (§11.1.2), and node relocation decoupled from the solid is not
representable there (§11.9 item 5). The governing principle is **equidistribution** of a
user-chosen *monitor function* $\omega(z)>0$. All $r$-refinement code is in
`src/Adaptivity/r_refinement.jl`.

### 14.2.1 The equidistribution principle

Given a monitor $\omega(z)$ that is large where mesh resolution is desired, the
equidistribution principle places node $z_n^\ast$ so that each subinterval carries
an equal share of the total "monitor mass":

$$
\int_{0}^{z_n^\ast}\!\omega(z)\,dz
= \frac{n}{N+1}\int_{0}^{L}\!\omega(z)\,dz,
\qquad n = 1,\dots,N,
$$

where $N$ is the number of interior nodes (so there are $N+1$ cells). Equivalently,
each cell $[z_n^\ast,z_{n+1}^\ast]$ contains the same $\int\omega\,dz$. Where
$\omega$ is large the cells must be small to keep the integral constant — exactly
the desired refinement. The continuous monitor integral is discretized as a sum
over cells,

$$
\Omega_{\mathrm{total}} = \int_0^L\omega\,dz \approx \sum_i \omega_i\,h_i ,
$$

with $\omega_i$ the monitor evaluated at cell center $z^c_i$ and $h_i$ the cell
size.

### 14.2.2 Monitor functions

**Gradient monitor.** The most common monitor is built from the solution
gradient, regularized so it never vanishes:

$$
\omega^{\mathrm{grad}}(z) = \sqrt{1 + \alpha\,|\nabla\phi|^2}\, .
$$

The "$1+$" floor guarantees $\omega\ge 1$ everywhere (a uniform mesh in the
gradient-free limit), while $\alpha$ tunes the sensitivity to the gradient. This
is the classical *arc-length monitor* of moving-mesh PDE methods: in 1D,
$\int\sqrt{1+\alpha\,\phi_z^2}\,dz$ is (for $\alpha=1$) the arc length of the
solution curve, so equidistributing $\omega^{\mathrm{grad}}$ places nodes at equal
arc-length intervals along the solution — naturally clustering them on steep
fronts. The gradient is the same distance-weighted center difference used by the
indicators (`compute_gradient_1d`, `src/Adaptivity/error_indicators.jl`), with
one-sided forms at boundaries.

**Interface monitor.** A Gaussian bump concentrates nodes near a prescribed
interface $z_{\mathrm{int}}$ irrespective of the solution:

$$
\omega^{\mathrm{iface}}(z) = \omega_{\mathrm{base}}
+ \omega_{\mathrm{peak}}\,\exp\!\left(-\frac{(z-z_{\mathrm{int}})^2}{2\sigma^2}\right).
$$

The defaults are $\sigma=10^{-3}\,$m, $\omega_{\mathrm{base}}=1$,
$\omega_{\mathrm{peak}}=10$, so the monitor is ten-fold larger at the interface
than in the far field, with a Gaussian width $\sigma$. Here $\sigma$ is the
interface-monitor width $\sigma_{\mathrm{mon}}$ of §2.

**Composite monitor.** `CompositeMonitor` combines several monitors by `:max`,
`:sum`, or `:product`. The product rule is multiplicative reinforcement (both
monitors must be large); the max rule is the conservative union.

> **Implementation note.** `evaluate_monitor` (`src/Adaptivity/r_refinement.jl`)
> dispatches on `GradientMonitor`, `InterfaceMonitor`, and
> `CompositeMonitor`; monitors are subtypes of the abstract `MonitorFunction`.
> Like the indicators, monitor evaluation skips inactive cells. Note that monitor
> functions are problem-dependent and carry arbitrary units; the user is
> responsible for keeping the disparity in $\omega$ bounded, because a large
> dynamic range in $\omega$ makes the equidistribution solve ill-conditioned.

### 14.2.3 Computing target node positions

`compute_target_positions(mesh, ω)` inverts the discrete
equidistribution relation by piecewise-linear interpolation of the cumulative
monitor integral. The algorithm:

1. **Sort active cells by $z$**, building the list
   $\{(i, z^c_i, \omega_i, h_i)\}$ ordered by center position. This step is what
   makes $r$-refinement robust to the non-monotone cell ordering that AMR can
   leave behind.
2. **Form the cumulative monitor mass.** With $C_0=0$ and
   $C_k = C_{k-1} + \omega_k h_k$ for $k=1,\dots,M$ over the $M$ sorted cells,
   $C_M = \Omega_{\mathrm{total}}$. A companion array $z_k$ holds running cell
   right-boundary estimates, with $z_0 = z_{\min}$ and $z_M = z_{\max}$.
3. **Invert for each interior node.** For interior node $n=1,\dots,N$ (with
   $N$ the number of *active* interior nodes — nodes referenced by an active
   cell, minus the two boundaries; $h$-coarsening orphans are excluded via
   `_active_node_mask`), set the target cumulative value
   $\Omega_n^\ast = \tfrac{n}{N+1}\,\Omega_{\mathrm{total}}$, locate the cumulative
   interval $[C_{k},C_{k+1}]$ containing it by `searchsortedlast`, and
   linearly interpolate within it:

$$
z_n^\ast = z_k + \frac{\Omega_n^\ast - C_k}{C_{k+1}-C_k}\,(z_{k+1}-z_k).
$$

   If a cumulative interval is degenerate ($C_{k+1}=C_k$, a zero-monitor cell) the
   fraction defaults to $\tfrac12$.

The routine returns the empty vector if there are no interior nodes or if
$\Omega_{\mathrm{total}}\le 0$ (which then disables relocation for that step).

The inversion is *node-count* based ($\tfrac{n}{N+1}\Omega_{\mathrm{total}}$), not
*cell-count* based, which is the correct choice for placing the $N$ interior nodes.

### 14.2.4 Relaxed relocation with quality constraints

`relocate_nodes!` drives the interior nodes toward their targets
by damped fixed-point iteration. Per iteration it:

1. Evaluates the monitor and computes target positions.
2. For each *active* interior node $n$ (sorted by position; nodes orphaned by
   $h$-coarsening are skipped via `_active_node_mask`) takes a **relaxed step**

$$
z_n^{\mathrm{new}} = z_n^{\mathrm{old}} + \beta\,(z_n^\ast - z_n^{\mathrm{old}}),
\qquad \beta\in(0,1],
$$

   with default relaxation $\beta=0.5$. Relaxation prevents node overshoot and
   tangling when the monitor changes sharply between iterations.

3. **Enforces a minimum cell size** by clamping the new position against its
   already-updated neighbors:

$$
z_{n-1}^{\mathrm{new}} + h_{\min} \;\le\; z_n^{\mathrm{new}} \;\le\; z_{n+1}^{\mathrm{old}} - h_{\min},
$$

   with the domain boundaries $z_{\min},z_{\max}$ standing in for absent neighbors
   at the ends. The default $h_{\min}=10^{-8}\,$m. This clamp can prevent a node
   from reaching its target if the equidistribution would demand a cell thinner
   than $h_{\min}$.

4. Calls `update_geometry!(mesh)` to refresh centers, volumes, and face positions.

5. **Checks convergence** on the maximum *relative* node movement,

$$
\max_n \frac{|z_n^{\mathrm{new}} - z_n^{\mathrm{old}}|}{L} < \mathrm{tol},
$$

   default $\mathrm{tol}=10^{-3}$. The function returns the iteration count at
   convergence, or `max_iter` (default $10$) if not converged.

**Boundary nodes never move** — only interior nodes are relocated. The boundary
positions are physical domain limits and are fixed by construction.

A critical caveat: `relocate_nodes!` relocates nodes but does **not remap the
solution onto the moved mesh**. It reads the cell states once at entry
(`extract_fields`) to evaluate the monitor and then moves nodes; the cell
*state* (the $T$ and $\xi$ stored per cell) is left attached to its resized
cell, so a cell's mass $\xi_j V_i$ changes with its volume. The package has no
conservative remap; `relocate_nodes!` is a *geometry* operation, conservation of
the transported fields is the caller's responsibility, and `solve` never calls
it.

---

## 14.3 Depletion handling: parked cells (the production path)

The single AMR capability wired into `solve` is **depletion handling**, enabled by
`solve(…; handle_depletion = true)` (default `false`). As
pyrolysis consumes the solid, an ALE cell (§11) thins toward zero thickness; on a
fixed mesh a cell empties instead. A vanishing or emptied cell is numerically
catastrophic for two reasons. First, its matrix heat capacity
$\rho c_p^{\mathrm{eff}} V_i \to 0$ as the solid gasifies, so the energy equation
$\rho c_p^{\mathrm{eff}}\,\dot T_i = (\text{fluxes})$ becomes singular. Second, the
transport rows of a thin cell scale as $1/\Delta z_i^2$, so the condition number
of the implicit stage matrix $W = I/(\gamma\Delta t) - J$ grows as
$\Delta z_i^{-2}$: at $\Delta z \approx 10^{-7}$ m, $\operatorname{cond}(W)$
exceeds $10^{13}$ for any $\Delta t \gtrsim 10^{-3}$ s and Newton's linear solves
carry no usable digits, so the integrator halves $\Delta t$ until the cell's own
time scale is resolved. No Jacobian accuracy and no regularisation of the
properties changes this; the cell must leave the nonlinear system.

Depletion handling **parks** the cell: it merges the cell's mass and sensible
enthalpy into a neighbour, then moves the freed state slot above the live
surface as an inert ghost. The state vector, the mesh arrays, the
`StateLayout`, the sparsity pattern and the Jacobian caches keep their size for
the whole solve; every kernel runs over the live range, and the boundary
conditions act at the live surface.

### 14.3.1 Live range and ghost slots

With $n$ the (fixed) number of state slots and $n_{\rm live}$ =
`mesh.n_active_cells`:

- live cells are `1:n_live`, live nodes and faces `1:n_live+1`; face
  `n_live+1` is the exposed boundary face (`Face1D` with `cell_ids =
  (n_live, 0)`, `boundary_id = 1`; `mesh.boundary_faces[1] == [n_live+1]`) and
  cell `n_live` the boundary cell;
- ghost slots are `n_live+1:n` (cells) and `n_live+2:n+1` (nodes);
  `active_cell_mask[i] == (i ≤ n_live)`;
- a ghost carries $\xi_j = 0$ for every component, the temperature it was
  parked with, $\chi = 0$ and `m_dry_solid_initial = 0`. Under ALE its nodes
  sit on the live surface (zero thickness) and move with it,
  $\dot z_k = w_{n_{\rm live}+1}$ for $k > n_{\rm live}+1$; on a fixed mesh
  they keep their geometry;
- every residual kernel and every Jacobian operator loops over
  `1:mesh.n_active_cells` (cells) or `1:mesh.n_active_cells+1` (faces) and
  detects the top boundary through the face connectivity, so no kernel reads a
  ghost. The pack step writes exactly zero into the ghost $T$, $\xi$ and $\chi$
  rows. Their Jacobian rows are therefore zero and $W$ is block-diagonal with
  $I/(\gamma\Delta t)$ in the ghost block — perfectly conditioned — while the
  error estimate of a ghost component is identically zero;
- the live thickness is $z_{n_{\rm live}+1} - z_1$; every ledger, mass and
  energy sum runs over the live cells.

The correctness statement is an **equivalence**: the residual (and the
Jacobian) of the live rows of a parked state equal, to roundoff, those of the
same physical configuration built as an $n_{\rm live}$-cell problem
(`test/residual/test_park_equivalence.jl`, `test/jacobian/test_park_jacobian.jl`,
for the Eulerian, ALE, Darcy, Beer–Lambert and variable-cross-section modes).

### 14.3.2 Trigger

Under ALE a live cell is depleted when its current thickness falls below a
fixed fraction of its *initial* thickness,

$$
\Delta z_i = z_{i+1}-z_i \;<\; \tau_{\mathrm{dep}}\,\Delta z_{i,\mathrm{init}},
$$

with the default `depletion_threshold` $\tau_{\mathrm{dep}} = 0.05$. Node
positions are read from the ODE state, so detection sees the live integrator
state. The reference thicknesses are stored per state slot in
`DepletionCallbackState.reference` (computed once as $V_{i,\mathrm{init}}/A_0$),
so each cell is judged against its *own* initial thickness — the only detector
that is correct on a stretched mesh, where a single shared reference would, for
any cell-size ratio above $1/\tau_{\mathrm{dep}} = 20$, flag the fine cells at
$t = 0$ or never flag the coarse ones. When an interior cell parks, the
reference entries shift with the cells. `find_thin_cells` returns the thin live
cells surface first (descending $z$) and throws on a non-positive signed
thickness: an inverted cell is a mesh-motion failure, not depletion.

The threshold is not a lever for Newton failures. Raising $\tau_{\rm dep}$
merges the hottest cell into a colder one earlier, which moves the surface
temperature: on the Landau case V-CPL-02 (2 mm/s on 62.5 µm cells), going
from 0.05 to 0.3 shifts the Landau identity residual from $7\times10^{-4}$ to
$2\times10^{-3}$. The Newton failures on a fast ablation front have two
other causes, degenerate gas transport (§5.6) and a step controller that
regrows the step beyond the Newton's range (§15.5.3), and are removed there.
The default stays 0.05.

On a fixed mesh cells cannot thin; they empty. The trigger is the **top live
cell's condensed concentration**,
$\sum_{j\notin\mathrm{gas}}\xi_{j,n_{\rm live}} < \tau_{\mathrm{dep}}\,
\sum_{j\notin\mathrm{gas}}\xi_{j,n_{\rm live}}^{\mathrm{init}}$, with the same
threshold (a cell that starts with no condensed mass never parks). Only the
top live cell parks, so a fixed mesh burns a layer top down: the surface cell
empties and parks, and the cell beneath it becomes the surface. That requires
the surface cell to be the hottest, as under surface absorption. On the PMMA
burnout of V-BO-02 all 60 cells park one by one from 63 to 452 s, and the
maximum live-cell temperature stays at 915 K.

Under in-depth absorption the temperature maximum of a regressing
non-charring polymer lies below the surface, so the surface cell reacts
slowest. Under ALE it shrinks and parks, and the surface advances into the
hotter material. On a fixed mesh the surface cell keeps its volume and loses
concentration: it persists as a low-density skin while the cells beneath it
empty, and the layer gasifies from within. A void cell conducts with the gas
conductivity and absorbs nothing, and nothing in the model carries heat across
a void by radiation. The flux transmitted through the void column is deposited
in the first condensed cell below it, which can return that heat upward only
by conduction through the gas, and the temperatures run away. Continued past its first void, the PMMA
deck heats the top substrate cell past 1000 K within about 1.5 s and past
3000 K, the end of the property bracket, about 60 s later. From the first
interior void on, the model's equations have no solution with valid
properties.

A fixed-mesh solve therefore stops at its first interior void. After every
accepted step the depletion callback tests each live cell below the top live
cell for a matrix volume fraction
$\varphi_m = \sum_j \gamma_j\xi_j/\rho_j < 10^{-4}$ (`MIXING_BLEND_LO`), the
property model's own definition of a void (§5.1); a cell that starts with no
condensed mass is not tested. On a hit the callback warns that the material
has gasified from within, that the model carries no heat across a void and
that a regressing material belongs on the ALE mesh, and terminates the
integration (retcode `:Terminated`). The in-depth PMMA deck stops at
337.0 s, when the lowest PMMA cell becomes a void, with the maximum live-cell
temperature at 912 K (V-BO-05). A fixed
mesh with depletion is the porous-residue model (a char that keeps its
volume, burning top down); a material that regresses belongs on the ALE
mesh.

> **Implementation note.** `is_cell_thin`, `find_thin_cells`,
> `depletion_candidate` and `interior_void_cell` in
> `src/Solver/callbacks.jl`; `matrix_volume_fraction` in
> `src/Properties/mixture_rules.jl`. The callback *condition* fires whenever a
> candidate or (fixed mesh) an interior void exists; the *affect* terminates
> the integration on a void, parks when `n_live > min_cells` (default
> `min_cells = 2`) and **terminates the integration** otherwise.

### 14.3.3 Energy-conserving temperature merge

When cells $a,b$ of volumes $V_a,V_b$ are merged into one, the concentrations
combine by **volume-weighted averaging**, which conserves species mass exactly
because mass is the extensive product $\xi_j V$:

$$
\xi_{\mathrm{merged},j} = \frac{V_a\,\xi_{a,j} + V_b\,\xi_{b,j}}{V_{\mathrm{keep}}},
$$

with $V_{\mathrm{keep}} = V_a + V_b$ under ALE (the keeper spans both cells) and
$V_{\mathrm{keep}} = V_b$ on a fixed mesh (the keeper's volume is fixed; the
parked cell's mass raises its concentration). Either way
$\xi_{\mathrm{merged},j}V_{\mathrm{keep}} = \xi_{a,j}V_a + \xi_{b,j}V_b$.

Temperature is **intensive**, and a naive volume-weighted average of $T$ does
*not* conserve internal energy when the two cells have different heat capacities —
the very situation in depletion, where a hot, dense solid cell is merged with a
cool, mostly-converted (low-capacity) one. Pyrolysis.jl therefore conserves the
**sensible enthalpy** of every component

$$
H(T) \;=\; \sum_{j=1}^{N_C} m_j \int_{T_{\mathrm{ref}}}^{T}
c_{p,j}(T')\,\mathrm{d}T',
\qquad m_j = V\,\xi_j ,
$$

and `T_merged` solves $H_{\mathrm{merged}}(T_{\mathrm{merged}}) = H_a(T_a) +
H_b(T_b)$. Two choices are deliberate:

- **$\int c_p\,\mathrm{d}T$, not $c_p(T)\cdot T$.** For temperature-dependent
  $c_p$ the product form is not an enthalpy: with $c_p = a + bT$ it carries the
  $b$-term with coefficient 1 instead of $\tfrac12$, a few-K systematic bias
  per merge (always the same sign when a hot surface cell merges into a cooler
  interior one). The integral is evaluated per component by 4-point
  Gauss–Legendre quadrature — exact for polynomial $c_p$ up to degree 7, and
  park events are rare so the cost is irrelevant. The datum
  $T_{\mathrm{ref}}$ cancels because each component's total mass is identical
  before and after the merge.
- **Every component is included.** The merged enthalpy sums the gas with the
  condensed components, the same $\rho c_p$ the temperature equation stores
  (§5, §7).

The solve is a bracketed Newton iteration on
$[\min(T_a,T_b),\,\max(T_a,T_b)]$: $H$ is strictly monotone ($c_p > 0$), and
$H_{\mathrm{merged}}(T_{\min}) \le H_{\mathrm{total}} \le
H_{\mathrm{merged}}(T_{\max})$ guarantees the root lies inside the bracket, so
the iteration cannot fail. If neither cell carries any heat capacity (both
empty) the energy balance is degenerate and the routine falls back to the
volume-weighted average.

The accuracy stakes: merging a $1\,$m$^3$ condensed-rich cell at $800\,$K with
$\rho c_p = 2000\,$kJ m$^{-3}$K$^{-1}$ and a $1\,$m$^3$
mostly-converted cell at $400\,$K with
$\rho c_p = 500\,$kJ m$^{-3}$K$^{-1}$ gives, by the wrong volume-weighted rule,
$T=(800+400)/2 = 600\,$K, whereas conserving the enthalpy (constant
$c_p$, where the integral form reduces to the familiar
$\rho c_p$-weighted mean) gives

$$
T_{\mathrm{merged}} = \frac{2000\cdot 800 + 500\cdot 400}{2000 + 500} = 720\ \mathrm{K},
$$

a $120\,$K (20 % of the temperature difference) error avoided.

> **Implementation note.** `merge_temperatures_conservatively(V_a, state_a, V_b,
> state_b, material)` in `src/Adaptivity/cell_merge.jl`. The component heat
> capacity is `material.components[j].heat_capacity`, a callable $c_{p,j}(T)$.
> This is the **single** canonical merge: the park primitive `park_cell!` and
> the experimental `coarsen_cells!` both call it.

### 14.3.4 The park event

`park_cell!(mesh, i_park, i_keep; material, ale)` (`src/Adaptivity/cell_park.jl`)
performs the mesh side of a park; `park!` (`src/Solver/callbacks.jl`) drives it
and rewrites the state vector. The sequence:

1. **Neighbour selection.** `find_merge_neighbor(mesh, i_park)`
   (`src/Adaptivity/surface_depletion.jl`): a boundary cell takes its only
   neighbour; an interior cell takes the **interior** (lower-$z$) one. The
   substrate side has the more slowly varying thermal state, so absorbing the
   thin cell into the cooler, larger interior cell produces the less abrupt
   merged temperature. Neighbours are found through the cell's faces
   (`find_interior_neighbor` / `find_exterior_neighbor`; either is 0 at a
   domain boundary, which is how the live top cell reports no cell above it).
2. **Conservative merge** of the two `CellState`s per §14.3.3, with
   $V_{\mathrm{keep}}$ as defined there; `m_dry_solid_initial` is summed into the
   keeper so the per-cell pyrolysis-progress denominator stays consistent
   (§10). The freed slot becomes the ghost state ($\xi = 0$, $T = T_{\rm merged}$,
   zero volume under ALE, its own volume on a fixed mesh, `m_dry_solid_initial = 0`).
3. **ALE slot permutation.** Cells `i_park+1:n_live` move down one slot and
   the nodes above the shared node `max(i_park, i_keep)` move down one slot
   (per-cell and per-node *values*: `cell_states`, `cell_volumes`,
   `m_dry_solid_initial`, node positions and velocities). The 1-D connectivity
   is regular — cell $i$ always owns nodes and faces $(i, i+1)$ — so no cell,
   face or node reference changes; the freed cell slot `n_live` (the last live
   index) becomes the ghost and the freed node slot lands on the surface. On a
   fixed mesh only the top live cell parks and nothing shifts.
4. **Boundary flags.** Face `n_live + 1` (`n_live` the new count) becomes the exposed
   boundary face, `cell_ids = (n_live, 0)`, `boundary_id = 1`; the previous top
   face becomes an inert ghost face; the matching node flags follow;
   `mesh.boundary_faces[1]` is updated; `n_active_cells` is decremented and
   `active_cell_mask` cleared for the ghost; the mesh geometry is recomputed
   (`update_geometry!`) and the AMR invariants (§14.6) asserted.
5. **State vector.** $T$, $\xi$ and (ALE) $z$ are read back from the mesh for
   every slot — the mesh mirrors the accepted state when the callback runs, so
   live cells shift exactly as the mesh did and the ghost carries the merged
   temperature, zero concentrations and the surface position. The $\chi$ block
   is shifted the same way with
   $\chi_{\rm merged} = (m_a\chi_a + m_b\chi_b)/(m_a+m_b)$ weighted by the
   **initial dry-solid masses** (captured before step 2 sums them), because
   $\chi_i = m_{\mathrm{released},i}/m_{\mathrm{dry,init},i}$ makes
   $\sum_i \chi_i\,m_{\mathrm{dry,init},i}$ the extensive invariant; the ghost
   gets $\chi = 0$. The workspace's cached top boundary face id follows the
   mesh (it is shared with the Jacobian's AD workspace, so one update serves
   both), and `derivative_discontinuity!(integrator, true)` forces a derivative
   recomputation.

Mass, sensible enthalpy, $\chi$-weighted dry mass and (under ALE) live volume are
conserved to roundoff at every event (`test/adaptivity/test_cell_park.jl`).

### 14.3.5 The depletion callback

The live machinery is a `DiscreteCallback` checked at the end of every
accepted integrator step (`create_depletion_callback`,
`src/Solver/callbacks.jl`). Its condition is `depletion_candidate(...) != 0`,
or on a fixed mesh an interior void (`interior_void_cell(...) != 0`); its
`affect!` terminates the integration with a warning on an interior void
(§14.3.2), and otherwise parks **exactly one cell per step**, surface first,
so the integrator adapts its step between events, or terminates the
integration when only `min_cells` live cells remain. It records every event in the shared
`DepletionCallbackState` (`park_times`, `parked_cells`, `n_live`), which the
thickness-termination callback (live thickness) and solution extraction also
read.

After the park the `affect!` flags the state as modified and proposes
`DT_AFTER_PARK` ($10^{-6}$ s) for the next step. The merged keeper cell
conserves mass and enthalpy but sits off the quasi-steady gas-concentration and
temperature profiles of its neighbours, so the residual at the parked state is
orders of magnitude above that of the accepted pre-park state, and the
integrator's extrapolated stage predictor, built from the pre-park derivative,
starts far from the post-park solution. At the pre-park step size the Newton
fails until the controller has cut the step far enough; from the proposed step
it converges at once, and the controller grows the step back within a few
steps. The reset is a cost setting, not a condition for completing the run: on
the 60-park PMMA burnout (V-BO-01) the run takes 858 steps and 8 Newton
failures with it and 1686 steps and 649 failures without it, with the same
mass-loss history. A cascade (several cells reaching the threshold together)
proceeds one park per $10^{-6}$-s step.

The integrator drains `saveat` with the pre-modification interpolant before
any `DiscreteCallback.affect!` runs, so a state saved at exactly the park time
carries the pre-park layout: `n_live_cells(cb, t) = n − #{t_park < t}`.
Extraction NaN-pads the ghost slots of `T`, `ξ` and `z`, reports the live
thickness $z_{n_{\rm live}+1}-z_1$, sums mass and gas generation over the live
cells, and exposes `solution.extras.live_cells` and
`solution.extras.park_times`. No dense output is lost and `saveat`
interpolation is unaffected because the state never changes size.

### 14.3.6 The Jacobian through burnout

Because the state keeps its size, the structured Jacobian (§15.6) stays valid
through every park:

- the sparsity pattern is built once for $n$ slots and is a superset of every
  live pattern;
- every operator assembles the live range only, so ghost rows and columns of
  the sparse block are zero;
- the per-cell payload vectors of the structured couplings (the prefix sums of
  the ALE mesh velocity, the Beer–Lambert suffix sums, the $\chi$ rank-1) are
  zero above $n_{\rm live}$, so the sums stop at the live surface and every
  ghost node row reproduces the surface node row, exactly as the residual's
  $\dot z_{\rm ghost} = w_{n_{\rm live}+1}$ requires;
- the geometry AD pass (§15.7) differentiates the masked residual and is
  consistent by construction.

`verify_jacobian` against the dense ForwardDiff reference passes strictly on
parked states for `Structured(strategy = DirectSolve())`, and the parked
`ApproxSolve` assembly equals that of the compacted problem on the live block.

The ALE default Jacobian, `Structured(strategy = ApproxSolve())` (§15.10.2),
stays the default through parking. The `:full` ALE pattern carries the
prefix-sum couplings of every node to every cell below it, which makes the
factorized $W$ effectively dense and the sparse direct solve 60–90× more
expensive per step than the reduced pattern; parking keeps every live cell
above $\tau_{\rm dep}$ of its initial thickness, and the two strategies follow
the same trajectory. The post-park transient is not a
Jacobian effect: `DirectSolve` restarted from $5\times10^{-2}$ s after every
park fails its Newton more than 500 times on the 60-park PMMA burnout, while
the reduced strategy with the reset of §14.3.5 needs 858 steps and 8
failures at a hundredth of the wall time.

The same holds on fast ablation. On a surface regressing at mm/s through a
reaction front one cell wide (the Landau case V-CPL-02, $v = 2$ mm/s,
$\Delta z = 62.5$ µm) the dropped mesh-velocity couplings are larger, but the
reduced modified Newton still contracts by a factor of 0.03 to 0.18 per
iteration at the controller's steps of 2.5 to 5 ms, and the two strategies
follow the same trajectory: without the step cap of §15.5.3, 1088 steps and
894 Newton failures for the reduced strategy against 1064 and 823 for
`DirectSolve` with UMFPACK, which takes several times the wall time. Those
Newton failures are not a Jacobian effect either; a Newton that re-evaluates
the exact Jacobian at every iterate stalls on the same states. They come from
degenerate gas transport (§5.6) and from the step controller (§15.5.3).

---

## 14.4 $h$-refinement and $h$-coarsening (experimental)

The $h$-adaptivity cluster changes mesh *topology* by splitting a cell into two
children ($h$-refinement) or merging two refinement siblings back into their
parent ($h$-coarsening). It lives in `src/Experimental/` and is **not invoked by
`solve`**; it is exercised by tests and the experimental `AMRController`. A mesh
it refines cannot be passed to `solve`: refinement leaves inactive parent cells
among the active ones, while the solver's kernels assume the active cells are the
contiguous prefix `1:n_active_cells`; it also changes the number of degrees of
freedom, which the fixed-length state vector of `solve` cannot follow.

### 14.4.1 Single-cell refinement

`refine_cell!(mesh, cell_id)` (`src/Experimental/h_refinement.jl`)
splits a leaf cell $p$ into a left child $a$ and a right child $b$:

1. **New center node.** Place a node at the geometric midpoint
   $z_{\mathrm{center}} = \tfrac12(z_{\mathrm{left}} + z_{\mathrm{right}})$, with
   its reference coordinate the midpoint of the parent's node references.
2. **Two children.** Child $a$ spans `(left_node, center_node)` and inherits the
   parent's left face; child $b$ spans `(center_node, right_node)` and inherits
   the right face; a new interior face joins them. Both children carry
   `level = parent.level + 1` and `parent = cell_id`.
3. **Repoint the parent's boundary faces** to the children, and record the
   children on the parent (`children = (a,b)`), which marks the parent a non-leaf.
4. **Interpolate state to the children** (see below).
5. **Append geometry caches** (centers, half volumes, the new face position/area)
   and split the dry-solid denominator evenly,
   $m_{\mathrm{dry},a}=m_{\mathrm{dry},b}=\tfrac12 m_{\mathrm{dry},p}$.
6. **Update the active mask**: mark the parent inactive, append two active
   children, increment `n_active_cells` by one (net $-1+2$), and raise
   `max_level`.
7. **Assert AMR invariants.**

**Conservation of mass** is exact. Concentrations are *not* gradient-interpolated;
the children simply inherit the parent value, $\xi_a = \xi_b = \xi_p$, with equal
volumes $V_a = V_b = V_p/2$. Then

$$
V_a\,\xi_a + V_b\,\xi_b = \tfrac12 V_p\,\xi_p + \tfrac12 V_p\,\xi_p = V_p\,\xi_p ,
$$

so each species mass is conserved exactly.

**Temperature interpolation** is where the *code differs from its docstring*, a
discrepancy that must be recorded for "no mistakes." The docstring advertises a
linearly interpolated split,
$\phi_{a,b} = \phi_p \mp (\nabla\phi\cdot\Delta x)/4$. The code computes a gradient
$T_{\mathrm{grad}}$ and forms

$$
T_a = T_p - T_{\mathrm{grad}}\,\frac{|\Delta z|}{2},
\qquad
T_b = T_p + T_{\mathrm{grad}}\,\frac{|\Delta z|}{2},
\qquad
\Delta z = z_{\mathrm{center}} - z^c_{p} .
$$

But $z_{\mathrm{center}} = \tfrac12(z_{\mathrm{left}}+z_{\mathrm{right}})$ is the
*same point* as the parent cell center $z^c_p = \tfrac12(z_{\mathrm{left}} +
z_{\mathrm{right}})$. Hence **$\Delta z = 0$ identically**, the correction
vanishes, and in the current implementation

$$
T_a = T_b = T_p .
$$

The split is thus a piecewise-constant (zeroth-order) injection of the parent
temperature into both children. This is exactly energy-conserving for the
$T\!\cdot\!V$ proxy ($V_a T_a + V_b T_b = V_p T_p$), so `verify_refinement_conservation`
reports zero error, but it does not realize the advertised second-order
reconstruction. With variable heat capacity the constant-$T$ split is also
energy-conserving in the true $\rho c_p T V$ sense because $\rho c_p$ and $T$ are
identical in both children. The gradient routine `estimate_gradient` is correct
and ready; only the geometric $\Delta z$ argument nullifies it.

> **Implementation note.** `refine_cell!` and `estimate_gradient` in
> `src/Experimental/h_refinement.jl`. The conservation check
> `verify_refinement_conservation` compares parent vs. summed
> children for mass (per species, against `volume_ref`) and a $T\!\cdot\!V$ energy
> proxy.

### 14.4.2 Marking and 2:1 balance

`mark_cells_for_refinement(mesh, η, θ_refine, h_min, max_level)`
flags a cell when **all three** criteria hold:

$$
\eta_i > \theta_{\mathrm{refine}}\,\max_j\eta_j,
\qquad
h_i > 2\,h_{\min},
\qquad
\mathrm{level}_i < \mathrm{max\_level}.
$$

The first is the relative-error threshold (default $\theta_{\mathrm{refine}}=0.7$);
the second ensures the children will still exceed $h_{\min}$; the third caps the
refinement depth. The maximum $\eta$ is taken only over active cells, and the
routine returns empty if $\max\eta\le 0$.

`enforce_2to1_balance!(mesh, marked, max_level)` augments the
marked set so that, after refinement, adjacent cells differ by at most one
refinement level. It is a **fixed-point iteration**: for each currently-marked
cell it inspects both neighbors; if a neighbor is coarser
($\mathrm{level}_{\mathrm{nb}} < \mathrm{level}_i$) and below `max_level`, it marks
the neighbor; the loop repeats until no new cell is marked. The 2:1 balance keeps
level transitions gradual, which bounds the interpolation error and the
conditioning of the discretization across refinement jumps. Two known
limitations: the balance is enforced only on *newly* marked cells (an
already-coarse cell that becomes unbalanced after a complex refinement pattern is
not retroactively un-marked), and the fixed-point approach is conservative and can
over-refine.

`perform_h_refinement!(mesh, marked)` refines the marked cells in
decreasing index order (for determinism) and calls `update_geometry!` once at the
end.

### 14.4.3 Sibling coarsening

`coarsen_cells!(mesh, a, b; material)` (`src/Experimental/h_coarsening.jl`) merges two *siblings* (same parent, both leaves) back into their parent.
Its preconditions distinguish it sharply from the depletion park's adjacent-cell merge:
it requires a parent–child relationship (`can_coarsen` checks parent $\ne 0$, both
leaves), whereas the park (`park_cell!`, §14.3.4) works on any
adjacent live pair. The merge
algebra is identical to §14.3.3 — energy-conserving $T$ (or volume-weighted
without material), volume-weighted $\xi$, summed $m_{\mathrm{dry}}$ — but the
*topology* operation restores the parent: the parent's faces are repointed back to
the parent, the parent is marked active with `children = (0,0)`, the two children
are marked inactive, and `n_active_cells` decrements by one.

Crucially, `coarsen_cells!` **does not remove the orphaned intermediate node and
interior face from the arrays** — it only flips the mask. This is the deliberate
"mark-inactive" strategy (§14.6): removing elements would invalidate indices
throughout the code, so the experimental $h$-AMR path leaves orphans in place,
trading memory for index stability. (Contrast the depletion park, which keeps the
arrays' length and shifts the live values down so the live cells stay a contiguous
prefix, §14.3.4.)

`mark_cells_for_coarsening(mesh, η, θ_coarsen, h_max)` marks a
sibling pair when **both** siblings have low error and the merged cell is not too
large:

$$
\eta_i < \theta_{\mathrm{coarsen}}\,\max_j\eta_j,
\qquad
\eta_{\mathrm{sib}} < \theta_{\mathrm{coarsen}}\,\max_j\eta_j,
\qquad
h_i + h_{\mathrm{sib}} \le h_{\max},
$$

with default $\theta_{\mathrm{coarsen}}=0.1$. `perform_h_coarsening!` iterates the
pairs (re-checking both are still active) and refreshes geometry.

### 14.4.4 The AMR controller and combined $hr$-adaptation

`AMRController{EI,QM}` (`src/Experimental/amr_controller.jl`) bundles an error
indicator, a quality metric, the refine/coarsen thresholds, the size and level
bounds, and a check cadence. Its defaults are $\theta_{\mathrm{refine}}=0.7$,
$\theta_{\mathrm{coarsen}}=0.1$, $h_{\min}=10^{-6}\,$m, $h_{\max}=10^{-2}\,$m,
`max_level=4`, `check_every=10`. `adapt_mesh!(mesh, controller)` runs the full cycle:

1. Compute indicators on the current state.
2. Mark for refinement, enforce 2:1 balance, perform $h$-refinement.
3. *Recompute* indicators on the refined mesh.
4. Mark sibling pairs for coarsening, perform $h$-coarsening.
5. Check the quality metric; if it fails, apply Laplacian smoothing
   (`enforce_quality!`, §14.5) and report whether all violations cleared.

It returns `(n_refined, n_coarsened, quality_ok)`. `should_adapt`/`maybe_adapt!`
gate the cycle to fire every `check_every` steps; `maybe_adapt!` (like
`adapt_mesh!`) accepts a `material` keyword that is forwarded to the coarsening pass for the energy-conserving
temperature merge (without it, sibling merges fall back to volume-weighted-T
averaging). `create_default_controller(; h_target)` builds a controller tuned
for pyrolysis: a `CompositeIndicator` of a temperature gradient ($w=1$) and a
temperature curvature ($w=0.5$) combined by `:max` (a `ReactionRateIndicator`
cannot run inside the cycle — see §14.1.3), a `CompositeQualityMetric` of
aspect-ratio ($\le 2$), smoothness ($\le 0.5$), and minimum-size
($h_{\mathrm{target}}/10$), with $h_{\min}=h_{\mathrm{target}}/10$ and
$h_{\max}=100\,h_{\mathrm{target}}$.

`hr_adapt!(mesh; …)` (`src/Experimental/hr_adapt.jl`) runs $h$-adaptivity *first*
(to capture sharp features by topology change) then $r$-adaptivity *second* (to
smooth and optimize node positions via `relocate_nodes!`). The ordering reflects
the design principle that topology changes should precede node relocation so that
relocation operates on the final cell count.

---

## 14.5 Mesh-quality metrics

Quality metrics guard against discretizations that, while topologically valid, are
numerically dangerous: extreme size jumps degrade interpolation accuracy and
matrix conditioning; non-positive volumes signal tangling. `solve` does not call
them; they serve the experimental `AMRController` (§14.4.4) and diagnostics. All
metrics are subtypes of `AbstractQualityMetric` (`src/Materials/Materials.jl`) and live in
`src/Adaptivity/mesh_quality.jl`. Each metric supports `check_quality` (a boolean)
and `quality_violations` (the offending cell indices).

**Aspect-ratio metric.** Bounds the ratio of adjacent active cell sizes:

$$
\frac{\max(h_i,h_j)}{\min(h_i,h_j)} \le \mathrm{max\_ratio},
\qquad \text{for adjacent } i,j,
$$

default $\mathrm{max\_ratio}=2$. Steep size transitions amplify the leading
truncation-error term across a face and worsen the conditioning of the flux
stencil, so a ratio near unity is desirable. The check sorts active cells by
position and compares each consecutive pair.

**Smoothness metric.** Bounds the *relative* size jump:

$$
\frac{|h_j - h_i|}{\max(h_i,h_j)} \le \mathrm{max\_size\_jump},
$$

default $0.5$. This is closely related to the aspect ratio but measures the jump
as a fraction of the larger cell; it is the more natural bound when cell sizes vary
smoothly. (Note the violation test inside `quality_violations` for smoothness uses
$|h_j-h_i|/h_i$ and $|h_j-h_i|/h_j$ per cell — a per-cell normalization — while the
summary statistic `compute_smoothness` reports $|h_j-h_i|/h$ per neighbor.)

**Jacobian metric.** Ensures positive cell volumes (no inverted/tangled cells):

$$
V_i > \mathrm{min\_jacobian}, \qquad \text{default } \mathrm{min\_jacobian}=0 .
$$

The name follows finite-element usage where a non-positive element Jacobian marks
an inverted element; in 1D it is simply the cell volume.

**Minimum-cell-size metric.** Prevents pathologically small cells:

$$
h_i \ge h_{\min}, \qquad \text{default } h_{\min}=10^{-8}\ \mathrm{m} .
$$

**Composite quality metric.** `CompositeQualityMetric` passes only if *all*
sub-metrics pass; its violation set is the union.

**Enforcement by Laplacian smoothing.** `enforce_quality!(mesh, metric;
max_iter=10)` attempts to repair violations by relaxed Laplacian
smoothing of the interior nodes of violating cells:

$$
z_{\mathrm{new}} = \tfrac12\bigl(z_{\mathrm{old}} + \bar z_{\mathrm{neigh}}\bigr),
$$

where $\bar z_{\mathrm{neigh}}$ averages the positions of the nodes of all active
cells sharing the node. Boundary nodes are never moved. The routine iterates until
no violations remain or `max_iter` is reached, refreshing geometry each pass, and
returns the residual violation count. This is a *local, incremental* fix: it can
fail to converge for clustered or severe violations (e.g. after several large
refinements), in which case the AMR controller reports `quality_ok = false`.

**Summary statistics.** `mesh_quality_summary(mesh)` returns the active-cell count,
min/max cell sizes, max aspect ratio, max smoothness deviation, and whether all
volumes are positive — a one-call diagnostic of mesh health.

> **Implementation note.** All metric types, `check_quality`/`quality_violations`,
> `compute_aspect_ratios`/`compute_smoothness`, `enforce_quality!`, and
> `mesh_quality_summary` in `src/Adaptivity/mesh_quality.jl`; all of them mask out
> inactive cells.

---

## 14.6 AMR bookkeeping and invariants

The AMR data model is split for performance (see §13):
**immutable connectivity** (`Cell1D`, `Face1D` are immutable structs; the topology
graph does not change during a residual evaluation) and **mutable geometry**
(positions, volumes, areas, the active mask). Topology changes — the park inside
`solve`, refinement and coarsening outside it — happen only between residual
evaluations; the residual never mutates the mesh.

The central bookkeeping object is the **active-cell mask** `active_cell_mask::
BitVector`, with `n_active_cells` its population count. The hard rule — call it the
*active-cell contract* — is that **every AMR consumer must iterate via the mask (or
the masked accessors), never over a raw `eachindex(mesh.cells)` range**. The reason
is that the experimental $h$-refinement *appends* children to `mesh.cells` and
*marks-inactive* on coarsening without compacting, so `length(mesh.cells)` can
exceed `n_active_cells`. Iterating the raw range would touch stale parent/child
cells — the failure mode the contract (and the masked accessors) exists to
prevent.

`assert_amr_invariants(mesh)` (`src/Adaptivity/Adaptivity.jl`) is the
cheap consistency check called after every refine/coarsen/merge step. It asserts
that all five parallel per-cell arrays — `active_cell_mask`, `cell_states`,
`cell_volumes`, `cell_centers`, and the $\chi$ denominator `m_dry_solid_initial`
(empty allowed only for hand-built meshes that never sized it) — have length
equal to `length(mesh.cells)`, and
that `n_active_cells == count(active_cell_mask)`. It does not check connectivity
correctness; it catches the most common drift, a counter that disagrees with the
mask or a parallel array of the wrong length.

Two distinct slot strategies coexist, by design:

| Operation | Strategy | Arrays | Index validity |
|---|---|---|---|
| $h$-refinement / $h$-coarsening (experimental) | **mark-inactive** | grow / left in place | preserved (no shift) |
| depletion (production) | **park**: shift live values down, ghost at the end | fixed length | live range renumbered (shift), connectivity regular |

The production depletion path shifts because the live cells must stay a
contiguous prefix `1:n_active_cells` in spatial order (§14.3.1); the arrays keep
their length and the 1-D connectivity never changes. The experimental $h$-AMR
path leaves orphans (so its arrays grow monotonically across refinement cycles)
precisely to avoid invalidating indices held elsewhere. A periodic
full-compaction pass for the
$h$-AMR path is an open design question (§14.10).

No workspace buffer changes size on a park: the caches are allocated once for
$n$ slots, and the Jacobian's sparsity pattern is a superset of every live
pattern (§14.3.6).

---

## 14.7 Interaction with the solver

The relationship between depletion parking and the time integrator is mediated
entirely by callbacks (§15) and the read-only-mesh discipline of the residual (§3,
§13); the $h$- and $r$-adaptivity of §14.2 and §14.4 has no such link. Three
facts complete the picture:

1. **Mesh sync is callback-driven.** The residual and Jacobian assembly treat the
   mutable mesh as read-only; the `StepUpdateCallback` mirrors the integrator's
   accepted state (including ALE node positions and the
   $\chi$-driven cross-section area) back into the mesh at the end of each step,
   so that the depletion check, the thickness-termination check, and
   post-processing all see a consistent mesh. The depletion callback runs *after*
   the step-update callback in the `CallbackSet` (see §15).

2. **Termination guards the singular limit.** Alongside depletion, under ALE, the
   thickness-termination callback (`create_thickness_termination_callback`, a
   `ContinuousCallback` on the live thickness
   $z_{n_{\mathrm{live}}+1}-z_1 - \mathrm{min\_thickness}$, default
   $\mathrm{min\_thickness}=10^{-6}\,$m) stops the integration when the whole
   sample has receded to near-zero thickness, before $\rho c_p^{\mathrm{eff}}V$
   becomes singular globally. Depletion handles *per-cell* thinning; termination
   handles *whole-domain* exhaustion. Depletion also self-terminates when only
   `min_cells` remain.

3. **Depletion keeps the Jacobian.** A park rewrites the state in place and
   leaves its size unchanged (§14.3.6), so the structured Jacobian, its sparse
   linear solver and `verify_jacobian` all run through burnout. A park event is
   nonetheless a discrete rewrite of the state: outputs are smooth in the
   parameters between events and only piecewise smooth across them (§17.7.3).

---

## 14.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

The reference codes for condensed-phase pyrolysis take markedly different stances
on mesh adaptivity, because each makes a different choice about how to represent a
receding surface. Mapping their notation to ours on first use:

- **ThermaKin2Ds** (Stoliarov–Lyon) uses a *single moving reference point* and a
  dynamic-mesh adjustment to represent volume change: the deformation source
  $\int (1/\rho)(\partial\rho/\partial t)\,dx$ drives a redistribution of the mesh
  rather than a topology change. There is no error-indicator-driven AMR. The
  Pyrolysis.jl equivalent of ThermaKin's volume-change-driven motion is the ALE
  mesh velocity of §11; what Pyrolysis.jl adds on top is the *surface-depletion
  merge* (§14.3), which removes a cell once it thins past
  $\tau_{\mathrm{dep}}$ — an operation ThermaKin handles by its reference-point
  rescaling rather than by deleting cells. Both share the goal of keeping the
  discretization well-conditioned as the sample recedes; Pyrolysis.jl's merge is
  the discrete, mass/energy-conserving realization of that goal.

- **Gpyro** (Lautenberger–Fernandez-Pello) uses a variable-$\Delta z$
  control-volume mesh: the cell thicknesses themselves are state-like and shrink as
  the solid is consumed, with the front-most control volume allowed to vanish.
  Gpyro's grid stretching toward the heated face is the static analogue of
  Pyrolysis.jl's geometrically stretched mesh (§13.9) plus the interface-proximity
  indicator (§14.1.4) and the interface monitor (§14.2.2), which dynamically pull
  resolution toward the front. Gpyro's "cell disappears" behavior corresponds to
  Pyrolysis.jl's depleted-cell park; the energy-conserving temperature combine
  (§14.3.3) is the careful version of folding a vanishing control volume's energy
  into its neighbor.

- **FDS** (McGrattan et al.) solves the 1D solid-phase heat equation on a *fixed*
  through-thickness grid with a node-collapsing rule as the surface regresses (the
  cell nearest the surface is removed once its mass drops below a threshold), and
  swelling/shrinking is applied through a density-ratio rule on the cell sizes. The
  FDS surface-node removal is conceptually identical to Pyrolysis.jl's
  surface-depletion merge, with the difference that Pyrolysis.jl folds the removed
  cell's sensible enthalpy into its interior neighbor by the enthalpy-conserving
  merge of §14.3.3, whereas FDS's regression bookkeeping is tied to its fixed-grid
  layout.
  Neither FDS nor Gpyro carries an indicator-driven $h$-refinement of the kind in
  Pyrolysis.jl's experimental cluster (§14.4).

The common thread: all four codes must reconcile a fixed (or slowly varying)
numerical grid with a physically receding boundary, and all four ultimately
*remove* resolution at the surface as material is consumed. Pyrolysis.jl's
contribution is to make that removal a clean, conservation-checked,
energy-conserving cell merge with a documented thinness criterion, decoupled from
the rest of the physics by the callback architecture.

> Note: the reference comparison here is drawn from the general published
> descriptions of these codes, so the mapping above is at the level of modeling
> strategy rather than equation-by-equation correspondence. The Pyrolysis.jl side
> of every statement is grounded in the source as cited throughout this chapter.

---

## 14.9 Worked summary of the key relations

For reference, the binding discrete relations of this chapter, all verified against
the source:

**Error indicators** (`src/Adaptivity/error_indicators.jl`):

$$
\eta_i^{\mathrm{grad}} = w\,h_i|\nabla\phi|_i,
\quad
\eta_i^{\mathrm{curv}} = w\,h_i^2|\nabla^2\phi|_i,
\quad
\eta_i^{\mathrm{react}} = w\,h_i\max_j\frac{|r_j|_i}{\max(\xi_{j,i},10^{-10})},
\quad
\eta_i^{\mathrm{iface}} = \frac{w}{1+d_i/h_{\mathrm{target}}} .
$$

**Composite combination**:
$\eta^{\max}=\max_k\hat\eta_k,\ \eta^{\ell_2}=\sqrt{\textstyle\sum_k\hat\eta_k^2},\ \eta^{\Sigma}=\sum_k\hat\eta_k$, with $\hat\eta_k=\eta_k/\max_{\text{active}}\eta_k$.

**Equidistribution** (`src/Adaptivity/r_refinement.jl`):

$$
\int_0^{z_n^\ast}\!\omega\,dz = \frac{n}{N+1}\int_0^L\!\omega\,dz,
\qquad
\omega^{\mathrm{grad}}=\sqrt{1+\alpha|\nabla\phi|^2},
\qquad
\omega^{\mathrm{iface}}=\omega_{\mathrm{base}}+\omega_{\mathrm{peak}}e^{-(z-z_{\mathrm{int}})^2/2\sigma^2}.
$$

**Relaxed node step**:
$z_n^{\mathrm{new}} = z_n^{\mathrm{old}} + \beta(z_n^\ast - z_n^{\mathrm{old}})$,
clamped to $[z_{n-1}^{\mathrm{new}}+h_{\min},\,z_{n+1}^{\mathrm{old}}-h_{\min}]$,
$\beta=0.5$, $\mathrm{tol}=10^{-3}$.

**Thinness criterion** (`src/Solver/callbacks.jl`):
$\Delta z_i < \tau_{\mathrm{dep}}\,\Delta z_{i,\mathrm{init}}$, $\tau_{\mathrm{dep}}=0.05$.

**Conservative merge** (`src/Adaptivity/cell_merge.jl`):

$$
\xi_{\mathrm{merged},j} = \frac{V_a\xi_{a,j}+V_b\xi_{b,j}}{V_{\mathrm{keep}}},
\qquad
\sum_{j} m_j\!\int_{T_{\mathrm{ref}}}^{T_{\mathrm{merged}}}\! c_{p,j}\,dT'
= H_a(T_a) + H_b(T_b),
\qquad
m_{\mathrm{dry},\mathrm{keep}}\mathrel{+}=m_{\mathrm{dry},\mathrm{remove}},
$$

with $V_{\mathrm{keep}} = V_a+V_b$ under ALE and $V_b$ on a fixed mesh, $m_j = V\,\xi_j$ the mass of component $j$, the sum running over every component ($T_{\mathrm{merged}}$ by
bracketed Newton; reduces to the $\rho c_p$-weighted mean for constant $c_p$;
§14.3.3), and $\chi_{\mathrm{merged}} = (m_a\chi_a+m_b\chi_b)/(m_a+m_b)$ with $m$
the initial dry-solid mass (§14.3.4).

**Refinement marking** (`src/Experimental/h_refinement.jl`):
$\eta_i>\theta_{\mathrm{refine}}\max\eta \ \wedge\ h_i>2h_{\min}\ \wedge\ \mathrm{level}_i<\mathrm{max\_level}$,
$\theta_{\mathrm{refine}}=0.7$.

**Coarsening marking** (`src/Experimental/h_coarsening.jl`):
both siblings $\eta<\theta_{\mathrm{coarsen}}\max\eta \ \wedge\ h_i+h_{\mathrm{sib}}\le h_{\max}$,
$\theta_{\mathrm{coarsen}}=0.1$.

**Quality bounds** (`src/Adaptivity/mesh_quality.jl`):
$\max(h_i,h_j)/\min(h_i,h_j)\le 2$ (aspect),
$|h_j-h_i|/\max(h_i,h_j)\le 0.5$ (smoothness),
$V_i>0$ (Jacobian),
$h_i\ge 10^{-8}$ (min size).

**AMR invariant** (`src/Adaptivity/Adaptivity.jl`):
$\mathrm{len}(\text{mask})=\mathrm{len}(\text{states})=\mathrm{len}(V)=\mathrm{len}(\text{centers})=\mathrm{len}(\text{cells})$
and $\mathrm{n\_active\_cells}=\mathrm{count}(\text{mask})$.

---

## 14.10 Limitations and open issues

The adaptivity subsystem is functional in the depletion-parking form `solve` uses
and experimental in its full $h$-/$r$-adaptive form, which `solve` does not use. The following are the principal
open issues, several of which are discrepancies or design questions surfaced by
the code itself.

1. **$h$-refinement temperature split is zeroth-order, not the advertised
   first-order.** Because the split point coincides with the parent cell center,
   $\Delta z = 0$ and the gradient correction in `refine_cell!` vanishes
   identically: both children inherit $T_a=T_b=T_p$ (§14.4.1). This is
   energy-conserving but loses the second-order reconstruction the docstring
   describes. A correct fix would interpolate to the *child* centers
   ($z^c_p \pm \Delta z_p/4$) rather than to the split point.

2. **Mesh-compaction strategy for $h$-AMR.** The experimental refine/coarsen path
   marks inactive without compacting, so `mesh.cells` grows monotonically across
   adaptation cycles. Whether to add a periodic full-compaction pass (reclaiming
   memory at the cost of a global reindex) or to accept the growth for
   long-running simulations is unresolved (§14.6).

3. **Reaction-rate indicator coupling.** The reaction-rate indicator requires
   the per-cell rate matrix (`rates = r`) and therefore cannot run inside the
   `adapt_mesh!` cycle, where indicators are recomputed after topology changes
   and any precomputed rates are stale (it throws rather than running silently).
   Integrating it with the implicit time-stepping — sampling the rates at the
   adaptation point — is not yet done (§14.1.3).

4. **Sensitivities across park events.** A park is a discrete rewrite of the
   state, so trajectories are only piecewise smooth in the parameters across an
   event (§17.7.3), and a finite difference whose perturbed solves park at
   different times measures the jump.

5. **$r$-refinement defaults are untuned, and relocation does not remap.**
   The convergence tolerance ($10^{-3}$) and relaxation ($0.5$) defaults have not
   been validated on realistic steep-front pyrolysis problems, and
   `relocate_nodes!` moves nodes without conservatively remapping $T,\xi$ onto the
   new layout — that remap is the caller's responsibility (§14.2.4). A non-smooth
   relocation or a missing remap would also make outputs non-smooth in the
   parameters (§17).

6. **2:1 balance is one-sided and Laplacian smoothing is fragile.** Balance is
   enforced only on newly marked cells and may leave residual imbalances for
   complex patterns (§14.4.2); the Laplacian quality enforcer is a local
   incremental fix that can fail to converge for clustered or severe violations
   (§14.5).
