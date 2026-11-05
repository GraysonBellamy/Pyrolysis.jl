# V-HC-05 — steady Kirchhoff profile with k(T) = k₀(1+βT)

**Tags:** `heat_conduction`, `temperature_dependent_k`, `kirchhoff`, `steady`

**References:**

- Kirchhoff transform (exact steady flux/profile), standard texts

## Problem statement

Steady conduction across a slab with k(T) = k₀(1 + βT), fixed T at both
faces. The Kirchhoff transform gives the exact flux and implicit profile:

  q = (k₀/L)·[(T_top − T_bot) + β(T_top² − T_bot²)/2]
  z(T) = k₀[(T − T_bot) + β(T² − T_bot²)/2] / q          (bottom cold face)

The implicit form is compared invert-free: evaluate z(T_num[i]) and compare
to the cell position. k varies 2.25× across 300–800 K with β = 0.01 (the
linear form caps the ratio at T_top/T_bot; 2.25× is char-vs-virgin scale).

What is verified: the T-dependent property path (LinearProperty) and the
distance-weighted harmonic-mean face interpolation of cell-centered k(T)
under strong nonlinearity. Note k(T) is evaluated at cell centers and
harmonic-averaged to faces — NOT the integrated Kirchhoff face conductance —
so unlike V-HC-04 this is a genuine discretization-error measurement, with
a convergence check. Transient T-dependent {k,c} has no exact solution
(that regime is B-02, code-to-code).

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ Kirchhoff-profile error (z(T) map, /L) | 2.464e-06 | — | 2.464e-06 | 0.00012 | yes | implicit profile compared invert-free; n=80 uniform |
| steady flux vs Kirchhoff-exact | 3.25e+04 | 3.25e+04 | 1.096e-06 | 6e-05 | yes | q = (k₀/L)[ΔT + βΔ(T²)/2] = 32500.0 W/m² |
| steadiness (last half-interval drift) | 2.274e-13 | — | 2.274e-13 | 1e-06 | yes | t_end = 10·L²/α_min |

![V-HC-05 T_profile_steady](assets/V-HC-05/V-HC-05_T_profile_steady.png)

![V-HC-05 T_profile_steady_error](assets/V-HC-05/V-HC-05_T_profile_steady_error.png)

![V-HC-05 flux_history](assets/V-HC-05/V-HC-05_flux_history.png)

![V-HC-05 flux_history_error](assets/V-HC-05/V-HC-05_flux_history_error.png)

![V-HC-05 convergence](assets/V-HC-05/V-HC-05_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-05/V-HC-05_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-05 T_profile_steady_vs_codes_error](assets/V-HC-05/V-HC-05_T_profile_steady_vs_codes_error.png)

#### Cost at matched accuracy

Each external code appears at its matched mesh h and time step Δt, at
coarser steps (2Δt, 4Δt) and at refinements (Δt/2, Δt/4, h/2 at Δt/4,
h/2 at Δt/2); this code at its report ladder, the external meshes and any
finer cost-only meshes, each at its report tolerances and at 10, 100 and
1000 times looser ones, the counterpart of the coarser steps. A coarse
configuration a code cannot express or complete is left out; for the
external codes the reason is in the matched provenance. The error is the
largest analytical error over
the curves and references the external run is scored on. Times are FDS's
internal elapsed time, which excludes its MPI start-up; Gpyro's and
ThermaKin's process time; and this code's warm problem construction and
solve. Each is the median of repeated runs on one pinned core
(`scripts/timing/README.md`; data and provenance in
`test/verification/reference/matched/timing/`). On the smallest problems
the external times are mostly input, setup and output. In the figure,
filled markers are each external code's matched mesh and step and this
code's report tolerances, open markers are the other runs, and this
code's line joins its cheapest configuration at each accuracy; lower left
is better. Errors at or below 1e-12 count as exact
when matching accuracy, errors at or below 1e-14 are drawn at 1e-14, and
runs every code solves exactly are not drawn.

![V-HC-05 cost](assets/V-HC-05/V-HC-05_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.03431 | 0.14 | 0.054 | 80, 1000× tol | 0.00823 | 6.56 |
| FDS 6.11.1 | h, 2Δt | 0.03431 | 0.097 | 0.054 | 80, 1000× tol | 0.00823 | 6.56 |
| FDS 6.11.1 | h, 4Δt | 0.03431 | 0.054 | 0.054 | 80, 1000× tol | 0.00823 | 6.56 |
| FDS 6.11.1 | h, Δt/2 | 0.03431 | 0.23 | 0.054 | 80, 1000× tol | 0.00823 | 6.56 |
| FDS 6.11.1 | h, Δt/4 | 0.03431 | 0.375 | 0.054 | 80, 1000× tol | 0.00823 | 6.56 |
| FDS 6.11.1 | h/2, Δt/4 | 0.02887 | 0.555 | 0.365 | 160, 1000× tol | 0.01704 | 21.4 |
| FDS 6.11.1 | h/2, Δt/2 | 0.02887 | 0.365 | 0.365 | 160, 1000× tol | 0.01704 | 21.4 |
| Gpyro V2 | h, Δt | 0.5236 | 0.08788 | 0.02438 | 20, 1000× tol | 0.002405 | ≥ 10.1 |
| Gpyro V2 | h, 2Δt | 0.5236 | 0.04426 | 0.02438 | 20, 1000× tol | 0.002405 | ≥ 10.1 |
| Gpyro V2 | h, 4Δt | 0.5236 | 0.02438 | 0.02438 | 20, 1000× tol | 0.002405 | ≥ 10.1 |
| Gpyro V2 | h, Δt/2 | 0.5236 | 0.1714 | 0.02438 | 20, 1000× tol | 0.002405 | ≥ 10.1 |
| Gpyro V2 | h, Δt/4 | 0.5236 | 0.3417 | 0.02438 | 20, 1000× tol | 0.002405 | ≥ 10.1 |
| Gpyro V2 | h/2, Δt/4 | 0.2636 | 0.5986 | 0.2991 | 40, 1000× tol | 0.004292 | 69.7 |
| Gpyro V2 | h/2, Δt/2 | 0.2636 | 0.2991 | 0.2991 | 40, 1000× tol | 0.004292 | 69.7 |
| ThermaKin2Ds 7 | h, Δt | 5.046 | 1.391 | 0.347 | 20, 1000× tol | 0.002405 | ≥ 144.0 |
| ThermaKin2Ds 7 | h, 2Δt | 5.046 | 0.6946 | 0.347 | 20, 1000× tol | 0.002405 | ≥ 144.0 |
| ThermaKin2Ds 7 | h, 4Δt | 5.046 | 0.347 | 0.347 | 20, 1000× tol | 0.002405 | ≥ 144.0 |
| ThermaKin2Ds 7 | h, Δt/2 | 5.046 | 2.779 | 0.347 | 20, 1000× tol | 0.002405 | ≥ 144.0 |
| ThermaKin2Ds 7 | h, Δt/4 | 5.046 | 5.574 | 0.347 | 20, 1000× tol | 0.002405 | ≥ 144.0 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 2.531 | 11.09 | 5.556 | 20, 1000× tol | 0.002405 | ≥ 2310.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 2.531 | 5.556 | 5.556 | 20, 1000× tol | 0.002405 | ≥ 2310.0 |

Solution overlays including the other codes' points: [T_profile_steady_vs_codes](assets/V-HC-05/V-HC-05_T_profile_steady_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.00642 | 0.0002818 | 0.0005775 |
| 40 | 0.025 | 0.0129 | 7.051e-05 | 0.0001511 |
| 80 | 0.0125 | 0.02577 | 1.763e-05 | 3.87e-05 |
| 160 | 0.00625 | 0.0555 | 4.408e-06 | 9.794e-06 |
| 320 | 0.003125 | 0.1185 | 1.102e-06 | 2.464e-06 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
