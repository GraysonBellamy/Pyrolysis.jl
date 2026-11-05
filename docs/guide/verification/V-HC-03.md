# V-HC-03 — finite slab, constant flux, adiabatic back

**Tags:** `heat_conduction`, `heat_flux_bc`, `energy_balance`

**References:**

- Carslaw & Jaeger (1959) §3.8 (constant-flux slab series)

## Problem statement

Complements the semi-infinite cases at LONG times: the thermal wave reaches
the back face — the regime cone-calorimeter simulations actually live in.
Constant q″ at the top face, insulated back; exact solution
T = T₀ + (q″L/k)·φ(x̃, Fo) with φ from exact/slab_series.jl (C&J §3.8),
x̃ = z/L from the adiabatic back face.

Two parameter-free exact targets beyond the profiles:
- front–back plateau: ΔT → q″L/2k for Fo ≫ 1 (100 K here),
- global energy: mean temperature rises exactly at q″/(ρcL) (0.25 K/s here).

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.1 | 0.0001229 | — | 0.0001229 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| L∞ T-profile error, Fo=0.5 | 9.176e-05 | — | 9.176e-05 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| L∞ T-profile error, Fo=1.5 | 8.138e-05 | — | 8.138e-05 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| front–back ΔT plateau, Fo=1.5 | 99.69 | 99.69 | 3.535e-09 | 1e-05 | yes | steady-shape plateau q″L/2k = 100 K between faces |
| mean-temperature rise (global energy) | 200 | 200 | 1.307e-12 | 1e-05 | yes | exact constant q″/(ρcL) = 0.25 K/s; integrator-tolerance limited |

![V-HC-03 T_profiles](assets/V-HC-03/V-HC-03_T_profiles.png)

![V-HC-03 T_profiles_error](assets/V-HC-03/V-HC-03_T_profiles_error.png)

![V-HC-03 front_back_dT_history](assets/V-HC-03/V-HC-03_front_back_dT_history.png)

![V-HC-03 front_back_dT_history_error](assets/V-HC-03/V-HC-03_front_back_dT_history_error.png)

![V-HC-03 mean_T_history](assets/V-HC-03/V-HC-03_mean_T_history.png)

![V-HC-03 mean_T_history_error](assets/V-HC-03/V-HC-03_mean_T_history_error.png)

![V-HC-03 convergence](assets/V-HC-03/V-HC-03_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-03/V-HC-03_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Gpyro V2 revision `bff0126c68f2` uses 321 boundary nodes (320 intervals, 0.0625 mm spacing) at dt = 0.1 s, with the flux absorbed at the surface node, re-radiation off and h = 0. Its front and back are surface nodes (x̃ = 1 and 0); its volume mean integrates the nodes over half end cells. ThermaKin2Ds 7 uses 320 uniform elements (0.0625 mm) at a fixed dt = 0.1 s; the absorbed flux is incident radiation q″/ε on the exposed element with ε = 1e-8.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-03 T_profiles_vs_codes_error](assets/V-HC-03/V-HC-03_T_profiles_vs_codes_error.png)

![V-HC-03 mean_T_history_vs_codes_error](assets/V-HC-03/V-HC-03_mean_T_history_vs_codes_error.png)

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

![V-HC-03 cost](assets/V-HC-03/V-HC-03_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Gpyro V2 | h, Δt | 0.01157 | 0.4957 | 0.4957 | 40, 1000× tol | 0.002071 | 239.0 |
| Gpyro V2 | h, 2Δt | 0.02269 | 0.2482 | 0.2482 | 40, 1000× tol | 0.002071 | 120.0 |
| Gpyro V2 | h, 4Δt | 0.04494 | 0.1293 | 0.1293 | 20, 1000× tol | 0.001539 | ≥ 84.0 |
| Gpyro V2 | h, Δt/2 | 0.006001 | 0.9843 | 0.9843 | 80, 1000× tol | 0.004 | 246.0 |
| Gpyro V2 | h, Δt/4 | 0.003218 | 1.966 | 1.966 | 80, 1000× tol | 0.004 | 492.0 |
| Gpyro V2 | h/2, Δt/4 | 0.002892 | 3.708 | 3.708 | 80, 1000× tol | 0.004 | 927.0 |
| Gpyro V2 | h/2, Δt/2 | 0.005674 | 1.854 | 1.854 | 80, 1000× tol | 0.004 | 464.0 |
| ThermaKin2Ds 7 | h, Δt | 0.000171 | 1.075 | 0.5399 | 320, 1000× tol | 0.01447 | 37.3 |
| ThermaKin2Ds 7 | h, 2Δt | 0.0001667 | 0.5399 | 0.5399 | 320, 1000× tol | 0.01447 | 37.3 |
| ThermaKin2Ds 7 | h, 4Δt | 0.0004303 | 0.2779 | 0.2779 | 320, 1000× tol | 0.01447 | 19.2 |
| ThermaKin2Ds 7 | h, Δt/2 | 0.000171 | 2.141 | 0.5399 | 320, 1000× tol | 0.01447 | 37.3 |
| ThermaKin2Ds 7 | h, Δt/4 | 0.000171 | 4.279 | 0.5399 | 320, 1000× tol | 0.01447 | 37.3 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 9.19e-05 | 8.523 | 4.261 | 640, 1000× tol | 0.02887 | 148.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 9.19e-05 | 4.261 | 4.261 | 640, 1000× tol | 0.02887 | 148.0 |

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-03/V-HC-03_T_profiles_vs_codes.png), [mean_T_history_vs_codes](assets/V-HC-03/V-HC-03_mean_T_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.002745 | 0.02373 | 0.0314 |
| 40 | 0.025 | 0.004601 | 0.005935 | 0.007862 |
| 80 | 0.0125 | 0.008721 | 0.001484 | 0.001966 |
| 160 | 0.00625 | 0.01686 | 0.000371 | 0.0004916 |
| 320 | 0.003125 | 0.03414 | 9.275e-05 | 0.0001229 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
