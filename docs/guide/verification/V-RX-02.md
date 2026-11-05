# V-RX-02 — zero-order decay (exact burnout time)

**Tags:** `kinetics`, `reaction_order`, `isothermal`

**References:**

- Khawam & Flanagan (2006), J. Phys. Chem. B 110:17315

## Problem statement

Khawam & Flanagan (2006) exact isothermal integrations (FDS verifies its
`part_baking_soda` case against the same forms). Sub-cases:

  a) n = 0    — linear decay, exact finite burnout time. Compared only up to
                0.9·t_b with the depletion limiter DISABLED: the power is
                taken on the finite base of the conversion, so at ξ = 0 the
                n = 0 rate stays k ξ₀ and would integrate ξ negative
                (`DepletionLimiter(enabled = false)`).
  b) n = 2    — single-species second-order decay (a C+C kinetic analog,
                not a reproduced ThermaKin run):
                ξ = ξ₀/(1 + kt).
  c) n = 2/3  — contracting-sphere R3 model with FDS's NaHCO₃ parameters
                A = 3.4e11 1/s, E = 103 kJ/mol at T = 450 K; the marker
                is m^(1/3) linear in t.
  d) n = 0 with the limiter ENABLED (default) — non-negativity guard: run to
                2× burnout and require 0 ≤ x ≤ limiter threshold at the end.
                This is the case the limiter exists for.

The rate is k·ξ₀·xⁿ on the conversion x = ξ/ξ₀, k in s⁻¹ for every order, so
published constants for normalized m = ξ/ξ₀ are the solver's as they stand
(see exact/kinetics.jl header).

## V-RX-02a — zero-order decay (exact burnout time)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=0.0) | 7.39e-15 | — | 7.39e-15 | 1e-12 | yes | linear decay; compared only to 0.9·t_b, limiter off |
| total-mass conservation | 7.105e-16 | — | 7.105e-16 | 1e-13 | yes | A→B closed system |

![V-RX-02a xi_histories](assets/V-RX-02a/V-RX-02a_xi_histories.png)

![V-RX-02a xi_histories_error](assets/V-RX-02a/V-RX-02a_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-RX-02a/V-RX-02a_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02a xi_histories_vs_codes_error](assets/V-RX-02a/V-RX-02a_xi_histories_vs_codes_error.png)

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

![V-RX-02a cost](assets/V-RX-02a/V-RX-02a_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 1.421e-14 | 0.183 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h, 2Δt | 1.421e-14 | 0.098 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h, 4Δt | 1.599e-14 | 0.055 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h, Δt/2 | 1.421e-14 | 0.355 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h, Δt/4 | 1.421e-14 | 0.702 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h/2, Δt/4 | 1.421e-14 | 0.705 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| FDS 6.11.1 | h/2, Δt/2 | 1.421e-14 | 0.361 | 0.055 | 1, 10× tol | 0.0005471 | 101.0 |
| Gpyro V2 | h, Δt | 1.421e-14 | 0.06201 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h, 2Δt | 1.421e-14 | 0.03327 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h, 4Δt | 1.421e-14 | 0.01918 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h, Δt/2 | 1.421e-14 | 0.1217 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h, Δt/4 | 1.421e-14 | 0.2395 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h/2, Δt/4 | 1.421e-14 | 0.2558 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |
| Gpyro V2 | h/2, Δt/2 | 1.421e-14 | 0.129 | 0.01918 | 1, 10× tol | 0.0005471 | 35.1 |

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02a/V-RX-02a_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02b — second-order decay (ThermaKin C+C analog)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=2.0) | 2.845e-12 | — | 2.845e-12 | 1e-11 | yes | ξ = ξ₀/(1+kt); k in s⁻¹ |
| total-mass conservation | 2.7e-15 | — | 2.7e-15 | 1e-13 | yes | A→B closed system |

![V-RX-02b xi_histories](assets/V-RX-02b/V-RX-02b_xi_histories.png)

![V-RX-02b xi_histories_error](assets/V-RX-02b/V-RX-02b_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-RX-02b/V-RX-02b_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02b xi_histories_vs_codes_error](assets/V-RX-02b/V-RX-02b_xi_histories_vs_codes_error.png)

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

![V-RX-02b cost](assets/V-RX-02b/V-RX-02b_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.01837 | 0.197 | 0.197 | 1, 1000× tol | 0.0007648 | ≥ 258.0 |
| FDS 6.11.1 | h, 2Δt | 0.03678 | 0.111 | 0.111 | 1, 1000× tol | 0.0007648 | ≥ 145.0 |
| FDS 6.11.1 | h, 4Δt | 0.07365 | 0.069 | 0.069 | 1, 1000× tol | 0.0007648 | ≥ 90.2 |
| FDS 6.11.1 | h, Δt/2 | 0.009184 | 0.367 | 0.367 | 1, 1000× tol | 0.0007648 | ≥ 480.0 |
| FDS 6.11.1 | h, Δt/4 | 0.004591 | 0.712 | 0.712 | 1, 1000× tol | 0.0007648 | ≥ 931.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.004591 | 0.716 | 0.712 | 1, 1000× tol | 0.0007648 | ≥ 931.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.009184 | 0.371 | 0.367 | 1, 1000× tol | 0.0007648 | ≥ 480.0 |
| Gpyro V2 | h, Δt | 0.01837 | 0.0624 | 0.0624 | 1, 1000× tol | 0.0007648 | ≥ 81.6 |
| Gpyro V2 | h, 2Δt | 0.03677 | 0.03288 | 0.03288 | 1, 1000× tol | 0.0007648 | ≥ 43.0 |
| Gpyro V2 | h, 4Δt | 0.07365 | 0.01877 | 0.01877 | 1, 1000× tol | 0.0007648 | ≥ 24.5 |
| Gpyro V2 | h, Δt/2 | 0.009184 | 0.1216 | 0.1216 | 1, 1000× tol | 0.0007648 | ≥ 159.0 |
| Gpyro V2 | h, Δt/4 | 0.004591 | 0.2415 | 0.2415 | 1, 1000× tol | 0.0007648 | ≥ 316.0 |
| Gpyro V2 | h/2, Δt/4 | 0.004591 | 0.2562 | 0.2415 | 1, 1000× tol | 0.0007648 | ≥ 316.0 |
| Gpyro V2 | h/2, Δt/2 | 0.009184 | 0.13 | 0.1216 | 1, 1000× tol | 0.0007648 | ≥ 159.0 |
| ThermaKin2Ds 7 | h, Δt | 4.944e-06 | 0.005601 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h, 2Δt | 4.944e-06 | 0.004461 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h, 4Δt | 4.944e-06 | 0.00397 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h, Δt/2 | 4.944e-06 | 0.008469 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h, Δt/4 | 4.944e-06 | 0.01299 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 4.944e-06 | 0.02773 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 4.944e-06 | 0.01533 | 0.00397 | 1, 1000× tol | 0.0007648 | ≥ 5.19 |

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02b/V-RX-02b_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02c — contracting-sphere R3 decay (FDS NaHCO₃, 450 K)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=0.6666666666666666) | 2.251e-10 | — | 2.251e-10 | 1e-09 | yes | m^(1/3) linear in t; k=0.3765 1/s at 450 K; exercises live Arrhenius (E ≠ 0) at pinned T |
| total-mass conservation | 2.501e-15 | — | 2.501e-15 | 1e-13 | yes | A→B closed system |

![V-RX-02c xi_histories](assets/V-RX-02c/V-RX-02c_xi_histories.png)

![V-RX-02c xi_histories_error](assets/V-RX-02c/V-RX-02c_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-RX-02c/V-RX-02c_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02c xi_histories_vs_codes_error](assets/V-RX-02c/V-RX-02c_xi_histories_vs_codes_error.png)

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

![V-RX-02c cost](assets/V-RX-02c/V-RX-02c_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.04966 | 0.182 | 0.182 | 1, 1000× tol | 0.0007584 | ≥ 240.0 |
| FDS 6.11.1 | h, 2Δt | 0.09934 | 0.099 | 0.099 | 1, 1000× tol | 0.0007584 | ≥ 131.0 |
| FDS 6.11.1 | h, 4Δt | 0.1987 | 0.055 | 0.055 | 1, 1000× tol | 0.0007584 | ≥ 72.5 |
| FDS 6.11.1 | h, Δt/2 | 0.02483 | 0.356 | 0.356 | 1, 1000× tol | 0.0007584 | ≥ 469.0 |
| FDS 6.11.1 | h, Δt/4 | 0.01242 | 0.693 | 0.693 | 1, 1000× tol | 0.0007584 | ≥ 914.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.01242 | 0.703 | 0.693 | 1, 1000× tol | 0.0007584 | ≥ 914.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.02483 | 0.356 | 0.356 | 1, 1000× tol | 0.0007584 | ≥ 469.0 |
| Gpyro V2 | h, Δt | 0.05556 | 0.06236 | 0.06236 | 1, 1000× tol | 0.0007584 | ≥ 82.2 |
| Gpyro V2 | h, 2Δt | 0.1052 | 0.03285 | 0.03285 | 1, 1000× tol | 0.0007584 | ≥ 43.3 |
| Gpyro V2 | h, 4Δt | 0.2046 | 0.01909 | 0.01909 | 1, 1000× tol | 0.0007584 | ≥ 25.2 |
| Gpyro V2 | h, Δt/2 | 0.03141 | 0.122 | 0.122 | 1, 1000× tol | 0.0007584 | ≥ 161.0 |
| Gpyro V2 | h, Δt/4 | 0.01972 | 0.2416 | 0.2416 | 1, 1000× tol | 0.0007584 | ≥ 319.0 |
| Gpyro V2 | h/2, Δt/4 | 0.01972 | 0.2567 | 0.2416 | 1, 1000× tol | 0.0007584 | ≥ 319.0 |
| Gpyro V2 | h/2, Δt/2 | 0.03141 | 0.1295 | 0.122 | 1, 1000× tol | 0.0007584 | ≥ 161.0 |

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02c/V-RX-02c_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02d — zero-order decay, depletion limiter ON (guard)

_limiter roll-off characterization_

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| non-negativity through 2× burnout | 5.369e-05 | 0 | 0 | 1e-09 | yes | the depletion limiter's job for n = 0 kinetics |
| residual mass within limiter threshold | 5.369e-05 | 0 | 5.369e-05 | 0.1 | yes | quadratic onset near 0 ⇒ algebraic tail; observed 5.4e-5 |

![V-RX-02d xi_histories](assets/V-RX-02d/V-RX-02d_xi_histories.png)

![V-RX-02d xi_histories_error](assets/V-RX-02d/V-RX-02d_xi_histories_error.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
