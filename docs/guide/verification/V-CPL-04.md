# V-CPL-04 — Frank-Kamenetskii slab: subcritical steady state + criticality bracket

**Tags:** `reaction_conduction`, `bifurcation`, `exothermic`, `coupled_jacobian`

**References:**

- D.A. Frank-Kamenetskii, Diffusion and Heat Exchange in Chemical Kinetics

## Problem statement

Zero-order exothermic Arrhenius source in a slab with Dirichlet walls
(Frank-Kamenetskii). The only case where reaction heat release feeds back
on the temperature field with a sharply known stability boundary
δ_c = 0.878 (slab), δ = (E/RT₀²)·S₀L²/k, L = half-width,
S₀ = |h|·A·e^(−E/RT₀).

Realization: the volumetric heat source is a SELF-REGENERATING reaction
A → A (order n = 0, h < 0) — product index = reactant index nets ξ̇ ≡ 0
exactly (verified: drift 0.0) while depositing |h|·A·ρ·e^(−E/RT) per unit
volume, i.e. a pure zero-order Arrhenius source with no reactant-depletion
clock. The conversion stays at one, where the depletion limiter's onset is
exactly 1.0, and the default gates (T_min = 0) are exactly open, so the
kinetics are pure Arrhenius. A T_max gate at T₀ + 100 K caps the
SUPERCRITICAL branch so the runaway saturates instead of chasing e^θ to
integrator death; the subcritical branch (rise ≤ 15 K) never feels it.

References (exact/frank_kamenetskii.jl, self-tested): the FK closed-form
profile at E/RT₀ → ∞, AND a full-Arrhenius RK4-shooting steady solution at
the case's finite E/RT₀ = 40 — the primary comparison, since the FK-limit
truncation error (~2–5% at E/RT₀ = 40) is physics, not discretization.

Assert-mode protocol:
  main solve  δ = 0.44 ≈ 0.50·δ_c → steady; center/interior T vs shooting
  bracket     δ = 0.70 ≈ 0.80·δ_c → must stay bounded (near-critical, stable)
              δ = 1.32 ≈ 1.50·δ_c → must run away past the detection rise
The δ_c bisection itself (turning-point hunt) is report-mode work.

### Quantities of interest (n = 160)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| steady center temperature vs full-Arrhenius shooting | 503.5 | 503.5 | 1.851e-05 | 0.0025 | yes | δ = 0.44 ≈ 0.5 δ_c, E/RT₀ = 40, relative to the rise |
| steady center rise vs FK closed form (limit check) | 0.2767 | 0.2772 | 0.001877 | 0.006 | yes | quantifies the FK-limit error at finite E/RT₀ |
| steady profile L∞ vs shooting (relative to center rise) | 3.059e-05 | — | 3.059e-05 | 0.0025 | yes | includes the symmetry check |
| steady-state reached (last-output drift) | 0.0003834 | 0 | 0.0003834 | 0.001 | yes | t_end = 5 diffusion times |
| near-critical δ = 0.80 δ_c stays on the stable branch | 0.552 | 0.5575 | 0.009852 | 0.03 | yes | bounded, near the FK lower branch (θ ≈ 0.558) |
| supercritical δ = 1.50 δ_c ignites (thermal runaway) | 104.2 | — | 0 | 0.5 | yes | max rise 104.2 K must exceed 50 K (θ_c ≈ 1.19 ⇒ 14.9 K is the largest stable rise) |

![V-CPL-04 steady_profiles](assets/V-CPL-04/V-CPL-04_steady_profiles.png)

![V-CPL-04 steady_profiles_error](assets/V-CPL-04/V-CPL-04_steady_profiles_error.png)

![V-CPL-04 center_T_histories](assets/V-CPL-04/V-CPL-04_center_T_histories.png)

![V-CPL-04 bifurcation_diagram](assets/V-CPL-04/V-CPL-04_bifurcation_diagram.png)

![V-CPL-04 bifurcation_diagram_error](assets/V-CPL-04/V-CPL-04_bifurcation_diagram_error.png)

![V-CPL-04 convergence](assets/V-CPL-04/V-CPL-04_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-CPL-04/V-CPL-04_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-CPL-04 steady_profiles_vs_codes_error](assets/V-CPL-04/V-CPL-04_steady_profiles_vs_codes_error.png)

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

![V-CPL-04 cost](assets/V-CPL-04/V-CPL-04_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FDS 6.11.1 | h, Δt | 0.0004645 | 0.276 | 0.078 | 80, 1000× tol | 0.002165 | 36.0 |
| FDS 6.11.1 | h, 2Δt | 0.0004638 | 0.144 | 0.078 | 80, 1000× tol | 0.002165 | 36.0 |
| FDS 6.11.1 | h, 4Δt | 0.0004629 | 0.078 | 0.078 | 80, 1000× tol | 0.002165 | 36.0 |
| FDS 6.11.1 | h, Δt/2 | 0.0004645 | 0.537 | 0.078 | 80, 1000× tol | 0.002165 | 36.0 |
| FDS 6.11.1 | h, Δt/4 | 0.0004645 | 1.068 | 0.078 | 80, 1000× tol | 0.002165 | 36.0 |
| FDS 6.11.1 | h/2, Δt/4 | 0.0001102 | 1.451 | 0.729 | 160, 1000× tol | 0.003703 | 197.0 |
| FDS 6.11.1 | h/2, Δt/2 | 0.0001102 | 0.729 | 0.729 | 160, 1000× tol | 0.003703 | 197.0 |
| Gpyro V2 | h, Δt | 0.001745 | 0.2624 | 0.06661 | 80, 1000× tol | 0.002165 | 30.8 |
| Gpyro V2 | h, 2Δt | 0.001743 | 0.1307 | 0.06661 | 80, 1000× tol | 0.002165 | 30.8 |
| Gpyro V2 | h, 4Δt | 0.00174 | 0.06661 | 0.06661 | 80, 1000× tol | 0.002165 | 30.8 |
| Gpyro V2 | h, Δt/2 | 0.001746 | 0.5271 | 0.06661 | 80, 1000× tol | 0.002165 | 30.8 |
| Gpyro V2 | h, Δt/4 | 0.001746 | 1.045 | 0.06661 | 80, 1000× tol | 0.002165 | 30.8 |
| Gpyro V2 | h/2, Δt/4 | 0.0007278 | 1.815 | 0.9063 | 80, 1000× tol | 0.002165 | 419.0 |
| Gpyro V2 | h/2, Δt/2 | 0.0007275 | 0.9063 | 0.9063 | 80, 1000× tol | 0.002165 | 419.0 |
| ThermaKin2Ds 7 | h, Δt | 0.1105 | 0.1545 | 0.04016 | 20, 1000× tol | 0.000944 | ≥ 42.5 |
| ThermaKin2Ds 7 | h, 2Δt | 0.1105 | 0.07779 | 0.04016 | 20, 1000× tol | 0.000944 | ≥ 42.5 |
| ThermaKin2Ds 7 | h, 4Δt | 0.1105 | 0.04016 | 0.04016 | 20, 1000× tol | 0.000944 | ≥ 42.5 |
| ThermaKin2Ds 7 | h, Δt/2 | 0.1105 | 0.3059 | 0.04016 | 20, 1000× tol | 0.000944 | ≥ 42.5 |
| ThermaKin2Ds 7 | h, Δt/4 | 0.1105 | 0.6091 | 0.04016 | 20, 1000× tol | 0.000944 | ≥ 42.5 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 0.0558 | 1.217 | 0.6127 | 20, 1000× tol | 0.000944 | ≥ 649.0 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 0.0558 | 0.6127 | 0.6127 | 20, 1000× tol | 0.000944 | ≥ 649.0 |

Solution overlays including the other codes' points: [steady_profiles_vs_codes](assets/V-CPL-04/V-CPL-04_steady_profiles_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | norm |
| --- | --- | --- | --- |
| 20 | 0.05 | 0.001344 | 0.01094 |
| 40 | 0.025 | 0.002202 | 0.002556 |
| 80 | 0.0125 | 0.003678 | 0.0004599 |
| 160 | 0.00625 | 0.006448 | 6.403e-05 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2.472 (norm), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-10 |
| reltol | 1.0e-8 |
