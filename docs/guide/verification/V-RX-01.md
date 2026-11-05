# V-RX-01 — two-step series chain A→B→C, isothermal

**Tags:** `kinetics`, `multi_step`, `isothermal`

**References:**

- FDS Verification Guide `two_step_solid_reaction`
- Lautenberger dissertation (2007) §3.4.4

## Problem statement

Identical problem in FDS (`two_step_solid_reaction`) and Gpyro §3.4.4:
constant rates k₁ = 0.389 s⁻¹, k₂ = 0.262 s⁻¹ (published values, so numbers
are cross-checkable against both codes), isothermal, exact Bateman solution
in exact/kinetics.jl.

Exercises multi-step reaction coupling: intermediate B is produced by
reaction 1 and consumed by reaction 2 SIMULTANEOUSLY (the existing
isothermal test is single-step).

Setup notes (solver conventions that would corrupt the case):
- constant rates realized as E = 0 (rate = A exactly, independent of T, so
  the Dirichlet-pinned temperature field cannot perturb the kinetics);
- depletion limiter DISABLED (tanh roll-off below ~1 kg/m³ would distort
  the late decay: ξ_A < 1 kg/m³ for t > ln(ξ₀)/k₁ ≈ 17.8 s);
- temperature gates left at defaults (fully open);
- h = 0 on both reactions keeps the single cell exactly isothermal.

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error | 4.773e-13 | — | 4.773e-13 | 1e-08 | yes | Bateman solution, k₁=0.389, k₂=0.262 (FDS/Gpyro values) |
| max relative ξ_B error | 1.207e-12 | — | 1.207e-12 | 1e-08 | yes | Bateman solution, k₁=0.389, k₂=0.262 (FDS/Gpyro values) |
| max relative ξ_C error | 7.295e-13 | — | 7.295e-13 | 1e-08 | yes | Bateman solution, k₁=0.389, k₂=0.262 (FDS/Gpyro values) |
| total-mass conservation (closed chain) | 1.251e-15 | — | 1.251e-15 | 1e-13 | yes | Σξ must equal ξ₀; unit yields |
| isothermality drift | 0 | — | 0 | 1e-12 | yes | h = 0 both reactions; Dirichlet-pinned single cell |

![V-RX-01 xi_histories](assets/V-RX-01/V-RX-01_xi_histories.png)

![V-RX-01 xi_histories_error](assets/V-RX-01/V-RX-01_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-RX-01/V-RX-01_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

FDS 6.11.1 and ThermaKin2Ds 7 use one insulated 1 mm cell (element); Gpyro V2 revision `bff0126c68f2` uses one ZEROD isothermal unknown. All three use the same rate constants with E = 0 and unit condensed yields, at dt = 0.001 s.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-01 xi_histories_vs_codes_error](assets/V-RX-01/V-RX-01_xi_histories_vs_codes_error.png)

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

![V-RX-01 cost](assets/V-RX-01/V-RX-01_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.08905 | 0.365 | 0.365 | 1, 1000× tol | 0.0007989 | ≥ 457.0 |
| FDS 6.11.1 | h, 2Δt | 0.1781 | 0.194 | 0.194 | 1, 1000× tol | 0.0007989 | ≥ 243.0 |
| FDS 6.11.1 | h, 4Δt | 0.3565 | 0.101 | 0.101 | 1, 1000× tol | 0.0007989 | ≥ 126.0 |
| FDS 6.11.1 | h, Δt/2 | 0.04452 | 0.723 | 0.723 | 1, 1000× tol | 0.0007989 | ≥ 905.0 |
| FDS 6.11.1 | h, Δt/4 | 0.02226 | 1.428 | 1.428 | 1, 1000× tol | 0.0007989 | ≥ 1790.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.02226 | 1.442 | 1.428 | 1, 1000× tol | 0.0007989 | ≥ 1790.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.04452 | 0.723 | 0.723 | 1, 1000× tol | 0.0007989 | ≥ 905.0 |
| Gpyro V2 | h, Δt | 0.09983 | 0.03667 | 0.03667 | 1, 1000× tol | 0.0007989 | ≥ 45.9 |
| Gpyro V2 | h, 2Δt | 0.2091 | 0.02096 | 0.02096 | 1, 1000× tol | 0.0007989 | ≥ 26.2 |
| Gpyro V2 | h, 4Δt | 0.4457 | 0.01259 | 0.01259 | 1, 1000× tol | 0.0007989 | ≥ 15.8 |
| Gpyro V2 | h, Δt/2 | 0.04833 | 0.06916 | 0.06916 | 1, 1000× tol | 0.0007989 | ≥ 86.6 |
| Gpyro V2 | h, Δt/4 | 0.02363 | 0.1368 | 0.1351 | 1, 1000× tol | 0.0007989 | ≥ 169.0 |
| Gpyro V2 | h/2, Δt/4 | 0.02363 | 0.1351 | 0.1351 | 1, 1000× tol | 0.0007989 | ≥ 169.0 |
| Gpyro V2 | h/2, Δt/2 | 0.04833 | 0.06964 | 0.06916 | 1, 1000× tol | 0.0007989 | ≥ 86.6 |
| ThermaKin2Ds 7 | h, Δt | 5.586e-05 | 0.01056 | 0.01056 | 1, 1000× tol | 0.0007989 | ≥ 13.2 |
| ThermaKin2Ds 7 | h, 2Δt | 7.548e-05 | 0.006872 | 0.006872 | 1, 1000× tol | 0.0007989 | ≥ 8.6 |
| ThermaKin2Ds 7 | h, 4Δt | 0.0001755 | 0.005276 | 0.005276 | 1, 1000× tol | 0.0007989 | ≥ 6.6 |
| ThermaKin2Ds 7 | h, Δt/2 | 4.966e-05 | 0.01718 | 0.01718 | 1, 1000× tol | 0.0007989 | ≥ 21.5 |
| ThermaKin2Ds 7 | h, Δt/4 | 4.966e-05 | 0.03038 | 0.01718 | 1, 1000× tol | 0.0007989 | ≥ 21.5 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 4.966e-05 | 0.07336 | 0.01718 | 1, 1000× tol | 0.0007989 | ≥ 21.5 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 4.966e-05 | 0.03824 | 0.01718 | 1, 1000× tol | 0.0007989 | ≥ 21.5 |

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-01/V-RX-01_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
