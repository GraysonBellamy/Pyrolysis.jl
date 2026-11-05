# V-HC-02 — finite slab, convective BC, Bi = 0.1

**Tags:** `heat_conduction`, `convective_bc`, `biot_sweep`

**References:**

- Carslaw & Jaeger (1959) §3.11 (eigenfunction series)
- FDS Verification Guide `heat_conduction_a–d`

## Problem statement

0.1 m slab initially at Tᵢ = 300 K, convective front (T∞ = 500 K), adiabatic
back, Bi = hL/k ∈ {0.1, 1, 10, 100} — four regimes from one definition:
Bi = 0.1 is lumped-like (BC-dominated), Bi = 100 is surface-Dirichlet-like
(interior-dominated). Exact eigenfunction series with eigencondition
λₙ tan λₙ = Bi (exact/slab_series.jl).

Mesh mapping: adiabatic face is z = 0, so the series coordinate is x̃ = z/L.
Comparison times are Fo = αt/L² ∈ {0.05, 0.1, 0.2, 0.5}, all well above the
60-root series validity floor (Fo ≳ 1.1e-3).

Convergence spot check runs on the Bi = 10 member only (the mixed regime);
the other three assert QoIs at the CI mesh.

## V-HC-02a — finite slab, convective BC, Bi = 0.1

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.05 | 1.442e-05 | — | 1.442e-05 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.05 | 304.9 | 304.9 | 5.69e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.1 | 1.177e-05 | — | 1.177e-05 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.1 | 306.9 | 306.9 | 3.896e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.2 | 1.362e-05 | — | 1.362e-05 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.2 | 309.7 | 309.7 | 2.798e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.5 | 6.881e-06 | — | 6.881e-06 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.5 | 315.6 | 315.6 | 2.603e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |

![V-HC-02a T_profiles](assets/V-HC-02a/V-HC-02a_T_profiles.png)

![V-HC-02a T_profiles_error](assets/V-HC-02a/V-HC-02a_T_profiles_error.png)

![V-HC-02a surface_T_history](assets/V-HC-02a/V-HC-02a_surface_T_history.png)

![V-HC-02a surface_T_history_error](assets/V-HC-02a/V-HC-02a_surface_T_history_error.png)

![V-HC-02a back_T_history](assets/V-HC-02a/V-HC-02a_back_T_history.png)

![V-HC-02a back_T_history_error](assets/V-HC-02a/V-HC-02a_back_T_history_error.png)

![V-HC-02a convergence](assets/V-HC-02a/V-HC-02a_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-02a/V-HC-02a_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-02a T_profiles_vs_codes_error](assets/V-HC-02a/V-HC-02a_T_profiles_vs_codes_error.png)

![V-HC-02a surface_T_history_vs_codes_error](assets/V-HC-02a/V-HC-02a_surface_T_history_vs_codes_error.png)

![V-HC-02a back_T_history_vs_codes_error](assets/V-HC-02a/V-HC-02a_back_T_history_vs_codes_error.png)

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

![V-HC-02a cost](assets/V-HC-02a/V-HC-02a_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.002464 | 0.104 | 0.104 | 80, 1000× tol | 0.002989 | 34.8 |
| FDS 6.11.1 | h, 2Δt | 0.004893 | 0.059 | 0.059 | 40, 1000× tol | 0.001767 | 33.4 |
| FDS 6.11.1 | h, 4Δt | 0.009754 | 0.035 | 0.035 | 40, 1000× tol | 0.001767 | 19.8 |
| FDS 6.11.1 | h, Δt/2 | 0.001251 | 0.195 | 0.195 | 80, 1000× tol | 0.002989 | 65.2 |
| FDS 6.11.1 | h, Δt/4 | 0.0006436 | 0.389 | 0.389 | 160, 1000× tol | 0.005359 | 72.6 |
| FDS 6.11.1 | h/2, Δt/4 | 0.0006156 | 0.42 | 0.42 | 160, 1000× tol | 0.005359 | 78.4 |
| FDS 6.11.1 | h/2, Δt/2 | 0.001223 | 0.212 | 0.212 | 80, 1000× tol | 0.002989 | 70.9 |
| Gpyro V2 | h, Δt | 0.001987 | 0.06888 | 0.06888 | 80, 1000× tol | 0.002989 | 23.0 |
| Gpyro V2 | h, 2Δt | 0.003062 | 0.03615 | 0.03615 | 80, 1000× tol | 0.002989 | 12.1 |
| Gpyro V2 | h, 4Δt | 0.005211 | 0.01996 | 0.01996 | 40, 1000× tol | 0.001767 | 11.3 |
| Gpyro V2 | h, Δt/2 | 0.00145 | 0.1359 | 0.1359 | 80, 1000× tol | 0.002989 | 45.5 |
| Gpyro V2 | h, Δt/4 | 0.001181 | 0.2687 | 0.2263 | 80, 1000× tol | 0.002989 | 75.7 |
| Gpyro V2 | h/2, Δt/4 | 0.0004966 | 0.4531 | 0.4531 | 160, 1000× tol | 0.005359 | 84.5 |
| Gpyro V2 | h/2, Δt/2 | 0.0007653 | 0.2263 | 0.2263 | 160, 1000× tol | 0.005359 | 42.2 |
| ThermaKin2Ds 7 | h, Δt | 0.118 | 0.07167 | 0.01956 | 20, 1000× tol | 0.001133 | ≥ 17.3 |
| ThermaKin2Ds 7 | h, 2Δt | 0.118 | 0.03696 | 0.01956 | 20, 1000× tol | 0.001133 | ≥ 17.3 |
| ThermaKin2Ds 7 | h, 4Δt | 0.118 | 0.01956 | 0.01956 | 20, 1000× tol | 0.001133 | ≥ 17.3 |
| ThermaKin2Ds 7 | h, Δt/2 | 0.118 | 0.1408 | 0.01956 | 20, 1000× tol | 0.001133 | ≥ 17.3 |
| ThermaKin2Ds 7 | h, Δt/4 | 0.118 | 0.2789 | 0.01956 | 20, 1000× tol | 0.001133 | ≥ 17.3 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 0.05922 | 0.5552 | 0.2791 | 20, 1000× tol | 0.001133 | ≥ 246.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 0.05922 | 0.2791 | 0.2791 | 20, 1000× tol | 0.001133 | ≥ 246.0 |

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-02a/V-HC-02a_T_profiles_vs_codes.png), [surface_T_history_vs_codes](assets/V-HC-02a/V-HC-02a_surface_T_history_vs_codes.png), [back_T_history_vs_codes](assets/V-HC-02a/V-HC-02a_back_T_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.001933 | 0.00213 | 0.003479 |
| 40 | 0.025 | 0.003297 | 0.0005325 | 0.000871 |
| 80 | 0.0125 | 0.006166 | 0.0001331 | 0.0002178 |
| 160 | 0.00625 | 0.01133 | 3.328e-05 | 5.446e-05 |
| 320 | 0.003125 | 0.02245 | 8.321e-06 | 1.362e-05 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2 (L2).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-HC-02b — finite slab, convective BC, Bi = 1.0

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.05 | 0.0001504 | — | 0.0001504 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.05 | 341.9 | 341.9 | 0.0002724 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.1 | 0.0001246 | — | 0.0001246 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.1 | 355.3 | 355.3 | 0.0001343 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.2 | 8.822e-05 | — | 8.822e-05 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.2 | 371.3 | 371.3 | 6.552e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.5 | 4.768e-05 | — | 4.768e-05 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.5 | 399.1 | 399.1 | 4.457e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |

![V-HC-02b T_profiles](assets/V-HC-02b/V-HC-02b_T_profiles.png)

![V-HC-02b T_profiles_error](assets/V-HC-02b/V-HC-02b_T_profiles_error.png)

![V-HC-02b surface_T_history](assets/V-HC-02b/V-HC-02b_surface_T_history.png)

![V-HC-02b surface_T_history_error](assets/V-HC-02b/V-HC-02b_surface_T_history_error.png)

![V-HC-02b back_T_history](assets/V-HC-02b/V-HC-02b_back_T_history.png)

![V-HC-02b back_T_history_error](assets/V-HC-02b/V-HC-02b_back_T_history_error.png)

![V-HC-02b convergence](assets/V-HC-02b/V-HC-02b_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-02b/V-HC-02b_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-02b T_profiles_vs_codes_error](assets/V-HC-02b/V-HC-02b_T_profiles_vs_codes_error.png)

![V-HC-02b surface_T_history_vs_codes_error](assets/V-HC-02b/V-HC-02b_surface_T_history_vs_codes_error.png)

![V-HC-02b back_T_history_vs_codes_error](assets/V-HC-02b/V-HC-02b_back_T_history_vs_codes_error.png)

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

![V-HC-02b cost](assets/V-HC-02b/V-HC-02b_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.02011 | 0.105 | 0.105 | 40, 1000× tol | 0.002031 | 51.7 |
| FDS 6.11.1 | h, 2Δt | 0.03806 | 0.058 | 0.058 | 40, 1000× tol | 0.002031 | 28.6 |
| FDS 6.11.1 | h, 4Δt | 0.07402 | 0.035 | 0.035 | 20, 1000× tol | 0.001337 | ≥ 26.2 |
| FDS 6.11.1 | h, Δt/2 | 0.01119 | 0.194 | 0.194 | 80, 1000× tol | 0.003493 | 55.5 |
| FDS 6.11.1 | h, Δt/4 | 0.006762 | 0.382 | 0.382 | 80, 1000× tol | 0.003493 | 109.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.005028 | 0.421 | 0.421 | 80, 1000× tol | 0.003493 | 121.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.009511 | 0.213 | 0.213 | 80, 1000× tol | 0.003493 | 61.0 |
| Gpyro V2 | h, Δt | 0.0066 | 0.06948 | 0.06948 | 80, 1000× tol | 0.003493 | 19.9 |
| Gpyro V2 | h, 2Δt | 0.01274 | 0.03603 | 0.03603 | 80, 1000× tol | 0.003493 | 10.3 |
| Gpyro V2 | h, 4Δt | 0.025 | 0.02015 | 0.02015 | 40, 1000× tol | 0.002031 | 9.92 |
| Gpyro V2 | h, Δt/2 | 0.005579 | 0.1361 | 0.1361 | 80, 1000× tol | 0.003493 | 38.9 |
| Gpyro V2 | h, Δt/4 | 0.005276 | 0.2689 | 0.2291 | 80, 1000× tol | 0.003493 | 65.6 |
| Gpyro V2 | h/2, Δt/4 | 0.001651 | 0.4509 | 0.4509 | 160, 1000× tol | 0.006649 | 67.8 |
| Gpyro V2 | h/2, Δt/2 | 0.003187 | 0.2291 | 0.2291 | 160, 1000× tol | 0.006649 | 34.5 |
| ThermaKin2Ds 7 | h, Δt | 0.7707 | 0.07194 | 0.02023 | 20, 1000× tol | 0.001337 | ≥ 15.1 |
| ThermaKin2Ds 7 | h, 2Δt | 0.7706 | 0.03743 | 0.02023 | 20, 1000× tol | 0.001337 | ≥ 15.1 |
| ThermaKin2Ds 7 | h, 4Δt | 0.7706 | 0.02023 | 0.02023 | 20, 1000× tol | 0.001337 | ≥ 15.1 |
| ThermaKin2Ds 7 | h, Δt/2 | 0.7707 | 0.1409 | 0.02023 | 20, 1000× tol | 0.001337 | ≥ 15.1 |
| ThermaKin2Ds 7 | h, Δt/4 | 0.7707 | 0.2793 | 0.02023 | 20, 1000× tol | 0.001337 | ≥ 15.1 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 0.3855 | 0.5563 | 0.2779 | 20, 1000× tol | 0.001337 | ≥ 208.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 0.3855 | 0.2779 | 0.2779 | 20, 1000× tol | 0.001337 | ≥ 208.0 |

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-02b/V-HC-02b_T_profiles_vs_codes.png), [surface_T_history_vs_codes](assets/V-HC-02b/V-HC-02b_surface_T_history_vs_codes.png), [back_T_history_vs_codes](assets/V-HC-02b/V-HC-02b_back_T_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.002445 | 0.01429 | 0.02246 |
| 40 | 0.025 | 0.004375 | 0.003572 | 0.005634 |
| 80 | 0.0125 | 0.007964 | 0.0008929 | 0.00141 |
| 160 | 0.00625 | 0.01569 | 0.0002232 | 0.0003528 |
| 320 | 0.003125 | 0.03132 | 5.581e-05 | 8.822e-05 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2 (L2).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-HC-02c — finite slab, convective BC, Bi = 10.0

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.05 | 0.001126 | — | 0.001126 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.05 | 453.5 | 453.5 | 0.0003044 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.1 | 0.0005972 | — | 0.0005972 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.1 | 465.9 | 465.9 | 0.0001471 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.2 | 0.0002715 | — | 0.0002715 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.2 | 475.5 | 475.5 | 5.104e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.5 | 0.0002269 | — | 0.0002269 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.5 | 487.1 | 487.1 | 1.973e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |

![V-HC-02c T_profiles](assets/V-HC-02c/V-HC-02c_T_profiles.png)

![V-HC-02c T_profiles_error](assets/V-HC-02c/V-HC-02c_T_profiles_error.png)

![V-HC-02c surface_T_history](assets/V-HC-02c/V-HC-02c_surface_T_history.png)

![V-HC-02c surface_T_history_error](assets/V-HC-02c/V-HC-02c_surface_T_history_error.png)

![V-HC-02c back_T_history](assets/V-HC-02c/V-HC-02c_back_T_history.png)

![V-HC-02c back_T_history_error](assets/V-HC-02c/V-HC-02c_back_T_history_error.png)

![V-HC-02c convergence](assets/V-HC-02c/V-HC-02c_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-02c/V-HC-02c_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-02c T_profiles_vs_codes_error](assets/V-HC-02c/V-HC-02c_T_profiles_vs_codes_error.png)

![V-HC-02c surface_T_history_vs_codes_error](assets/V-HC-02c/V-HC-02c_surface_T_history_vs_codes_error.png)

![V-HC-02c back_T_history_vs_codes_error](assets/V-HC-02c/V-HC-02c_back_T_history_vs_codes_error.png)

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

![V-HC-02c cost](assets/V-HC-02c/V-HC-02c_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.06147 | 0.104 | 0.104 | 80, 1000× tol | 0.004359 | 23.9 |
| FDS 6.11.1 | h, 2Δt | 0.1063 | 0.058 | 0.058 | 40, 1000× tol | 0.002409 | 24.1 |
| FDS 6.11.1 | h, 4Δt | 0.1964 | 0.035 | 0.035 | 40, 1000× tol | 0.002409 | 14.5 |
| FDS 6.11.1 | h, Δt/2 | 0.03938 | 0.199 | 0.199 | 80, 1000× tol | 0.004359 | 45.7 |
| FDS 6.11.1 | h, Δt/4 | 0.02851 | 0.379 | 0.218 | 80, 1000× tol | 0.004359 | 50.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.01537 | 0.424 | 0.424 | 160, 1000× tol | 0.008 | 53.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.02655 | 0.218 | 0.218 | 80, 1000× tol | 0.004359 | 50.0 |
| Gpyro V2 | h, Δt | 0.04356 | 0.06976 | 0.06976 | 80, 1000× tol | 0.004359 | 16.0 |
| Gpyro V2 | h, 2Δt | 0.08326 | 0.03633 | 0.03633 | 40, 1000× tol | 0.002409 | 15.1 |
| Gpyro V2 | h, 4Δt | 0.1646 | 0.01982 | 0.01982 | 40, 1000× tol | 0.002409 | 8.23 |
| Gpyro V2 | h, Δt/2 | 0.02422 | 0.1363 | 0.1363 | 80, 1000× tol | 0.004359 | 31.3 |
| Gpyro V2 | h, Δt/4 | 0.01455 | 0.2697 | 0.2697 | 160, 1000× tol | 0.008 | 33.7 |
| Gpyro V2 | h/2, Δt/4 | 0.01089 | 0.4553 | 0.4553 | 160, 1000× tol | 0.008 | 56.9 |
| Gpyro V2 | h/2, Δt/2 | 0.02081 | 0.2285 | 0.2285 | 80, 1000× tol | 0.004359 | 52.4 |
| ThermaKin2Ds 7 | h, Δt | 2.489 | 0.07214 | 0.07214 | 20, 1000× tol | 0.001528 | ≥ 47.2 |
| ThermaKin2Ds 7 | h, 2Δt | 2.49 | 0.03726 | 0.03726 | 20, 1000× tol | 0.001528 | ≥ 24.4 |
| ThermaKin2Ds 7 | h, 4Δt | 2.49 | 0.02002 | 0.02002 | 20, 1000× tol | 0.001528 | ≥ 13.1 |
| ThermaKin2Ds 7 | h, Δt/2 | 2.489 | 0.1408 | 0.07214 | 20, 1000× tol | 0.001528 | ≥ 47.2 |
| ThermaKin2Ds 7 | h, Δt/4 | 2.489 | 0.281 | 0.07214 | 20, 1000× tol | 0.001528 | ≥ 47.2 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 1.247 | 0.5557 | 0.2803 | 20, 1000× tol | 0.001528 | ≥ 183.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 1.247 | 0.2803 | 0.2803 | 20, 1000× tol | 0.001528 | ≥ 183.0 |

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-02c/V-HC-02c_T_profiles_vs_codes.png), [surface_T_history_vs_codes](assets/V-HC-02c/V-HC-02c_surface_T_history_vs_codes.png), [back_T_history_vs_codes](assets/V-HC-02c/V-HC-02c_back_T_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.003196 | 0.05208 | 0.06959 |
| 40 | 0.025 | 0.005479 | 0.013 | 0.01739 |
| 80 | 0.0125 | 0.01043 | 0.00325 | 0.004345 |
| 160 | 0.00625 | 0.02055 | 0.0008124 | 0.001086 |
| 320 | 0.003125 | 0.04156 | 0.0002031 | 0.0002715 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2.001 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-HC-02d — finite slab, convective BC, Bi = 100.0

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.05 | 0.001275 | — | 0.001275 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.05 | 495 | 495 | 6.105e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.1 | 0.0006256 | — | 0.0006256 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.1 | 496.4 | 496.4 | 2.152e-05 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.2 | 0.0002811 | — | 0.0002811 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.2 | 497.5 | 497.5 | 6.205e-06 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |
| L∞ T-profile error, Fo=0.5 | 0.0002559 | — | 0.0002559 | 0.06 | yes | observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K |
| surface temperature, Fo=0.5 | 498.8 | 498.8 | 2.813e-06 | 0.015 | yes | front face from the surface Newton solve vs series at x̃=1; observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K |

![V-HC-02d T_profiles](assets/V-HC-02d/V-HC-02d_T_profiles.png)

![V-HC-02d T_profiles_error](assets/V-HC-02d/V-HC-02d_T_profiles_error.png)

![V-HC-02d surface_T_history](assets/V-HC-02d/V-HC-02d_surface_T_history.png)

![V-HC-02d surface_T_history_error](assets/V-HC-02d/V-HC-02d_surface_T_history_error.png)

![V-HC-02d back_T_history](assets/V-HC-02d/V-HC-02d_back_T_history.png)

![V-HC-02d back_T_history_error](assets/V-HC-02d/V-HC-02d_back_T_history_error.png)

![V-HC-02d convergence](assets/V-HC-02d/V-HC-02d_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-HC-02d/V-HC-02d_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-HC-02d T_profiles_vs_codes_error](assets/V-HC-02d/V-HC-02d_T_profiles_vs_codes_error.png)

![V-HC-02d surface_T_history_vs_codes_error](assets/V-HC-02d/V-HC-02d_surface_T_history_vs_codes_error.png)

![V-HC-02d back_T_history_vs_codes_error](assets/V-HC-02d/V-HC-02d_back_T_history_vs_codes_error.png)

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

![V-HC-02d cost](assets/V-HC-02d/V-HC-02d_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.06782 | 0.103 | 0.103 | 80, 1000× tol | 0.005075 | 20.3 |
| FDS 6.11.1 | h, 2Δt | 0.116 | 0.058 | 0.058 | 40, 1000× tol | 0.002669 | 21.7 |
| FDS 6.11.1 | h, 4Δt | 0.2129 | 0.035 | 0.035 | 40, 1000× tol | 0.002669 | 13.1 |
| FDS 6.11.1 | h, Δt/2 | 0.0439 | 0.195 | 0.195 | 80, 1000× tol | 0.005075 | 38.4 |
| FDS 6.11.1 | h, Δt/4 | 0.03206 | 0.382 | 0.217 | 80, 1000× tol | 0.005075 | 42.8 |
| FDS 6.11.1 | h/2, Δt/4 | 0.01695 | 0.421 | 0.421 | 160, 1000× tol | 0.009476 | 44.4 |
| FDS 6.11.1 | h/2, Δt/2 | 0.02898 | 0.217 | 0.217 | 80, 1000× tol | 0.005075 | 42.8 |
| Gpyro V2 | h, Δt | 0.06758 | 0.06907 | 0.06907 | 80, 1000× tol | 0.005075 | 13.6 |
| Gpyro V2 | h, 2Δt | 0.1284 | 0.03635 | 0.03635 | 40, 1000× tol | 0.002669 | 13.6 |
| Gpyro V2 | h, 4Δt | 1.354e+07 | 0.04129 | 0.03635 | 20, 1000× tol | 0.001814 | ≥ 20.0 |
| Gpyro V2 | h, Δt/2 | 0.03728 | 0.1358 | 0.1358 | 80, 1000× tol | 0.005075 | 26.8 |
| Gpyro V2 | h, Δt/4 | 0.02232 | 0.2683 | 0.2683 | 80, 1000× tol | 0.005075 | 52.9 |
| Gpyro V2 | h/2, Δt/4 | 0.0169 | 0.4533 | 0.4533 | 160, 1000× tol | 0.009476 | 47.8 |
| Gpyro V2 | h/2, Δt/2 | 0.03214 | 0.2271 | 0.2271 | 80, 1000× tol | 0.005075 | 44.7 |
| ThermaKin2Ds 7 | h, Δt | 3.144 | 0.07284 | 0.0202 | 20, 1000× tol | 0.001814 | ≥ 11.1 |
| ThermaKin2Ds 7 | h, 2Δt | 3.144 | 0.03701 | 0.0202 | 20, 1000× tol | 0.001814 | ≥ 11.1 |
| ThermaKin2Ds 7 | h, 4Δt | 3.144 | 0.0202 | 0.0202 | 20, 1000× tol | 0.001814 | ≥ 11.1 |
| ThermaKin2Ds 7 | h, Δt/2 | 3.144 | 0.1406 | 0.0202 | 20, 1000× tol | 0.001814 | ≥ 11.1 |
| ThermaKin2Ds 7 | h, Δt/4 | 3.144 | 0.279 | 0.0202 | 20, 1000× tol | 0.001814 | ≥ 11.1 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 1.572 | 0.558 | 0.2799 | 20, 1000× tol | 0.001814 | ≥ 154.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 1.572 | 0.2799 | 0.2799 | 20, 1000× tol | 0.001814 | ≥ 154.0 |

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-02d/V-HC-02d_T_profiles_vs_codes.png), [surface_T_history_vs_codes](assets/V-HC-02d/V-HC-02d_surface_T_history_vs_codes.png), [back_T_history_vs_codes](assets/V-HC-02d/V-HC-02d_back_T_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.003507 | 0.0578 | 0.07198 |
| 40 | 0.025 | 0.006092 | 0.01444 | 0.01799 |
| 80 | 0.0125 | 0.01195 | 0.003609 | 0.004497 |
| 160 | 0.00625 | 0.02449 | 0.0009023 | 0.001124 |
| 320 | 0.003125 | 0.04892 | 0.0002256 | 0.0002811 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2 (L2).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
