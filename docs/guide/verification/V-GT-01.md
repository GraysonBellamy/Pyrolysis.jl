# V-GT-01 — Gaussian gas-pulse diffusion (plane source)

**Tags:** `gas_transport`, `lumped_transfer`, `mass_conservation`

**References:**

- J. Crank, The Mathematics of Diffusion (plane-source solution)
- Stoliarov & Lyon, FSS 9:1141 (2008), gas-transport verification

## Problem statement

Pure lumped wall transport of gas released inside the solid (Crank/ThermaKin
plane source). Instead of approximating a delta function, the initial
condition is a finite-width Gaussian — itself an exact solution — so the
comparison is IC-discretization-free:

  ξ(z,t) = M_tot/√(2πσ²(t)) · exp(−(z−z₀)²/(2σ²(t))),  σ²(t) = σ₀² + 2dt

At constant physical porosity φ and temperature, J=−(Λ/φ)∂zξ.
The Gaussian variance therefore grows with diffusivity d=Λ/φ. Adiabatic
faces, no reaction/source, and internally carried enthalpy preserve the
initially uniform temperature; this is checked in the QoIs.
The sealed image-source correction is below the spatial error on this domain.

QoIs: profile L∞ (relative to the evolving peak), total-mass conservation,
and the parameter-free variance growth law Var(t) = σ₀² + 2dt.

### Quantities of interest (n = 640)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| isothermal invariant | 0 | — | 0 | 1e-08 | yes | sealed, source-free uniform-temperature state |
| L∞ ξ-profile error / peak, t=2.0 s | 7.809e-05 | — | 7.809e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| L∞ ξ-profile error / peak, t=4.0 s | 4.826e-05 | — | 4.826e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| L∞ ξ-profile error / peak, t=8.0 s | 2.704e-05 | — | 2.704e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| gas-mass conservation (sealed) | 3.244e-15 | — | 3.244e-15 | 1e-13 | yes | antisymmetric face fluxes telescope exactly |
| variance growth, t=2.0 s | 5e-06 | 5e-06 | 8.809e-15 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |
| variance growth, t=4.0 s | 9e-06 | 9e-06 | 1.012e-10 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |
| variance growth, t=8.0 s | 1.7e-05 | 1.7e-05 | 4.597e-06 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |

![V-GT-01 xi_profiles](assets/V-GT-01/V-GT-01_xi_profiles.png)

![V-GT-01 xi_profiles_error](assets/V-GT-01/V-GT-01_xi_profiles_error.png)

![V-GT-01 variance_history](assets/V-GT-01/V-GT-01_variance_history.png)

![V-GT-01 variance_history_error](assets/V-GT-01/V-GT-01_variance_history_error.png)

![V-GT-01 convergence](assets/V-GT-01/V-GT-01_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-GT-01/V-GT-01_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-GT-01 xi_profiles_vs_codes_error](assets/V-GT-01/V-GT-01_xi_profiles_vs_codes_error.png)

![V-GT-01 variance_history_vs_codes_error](assets/V-GT-01/V-GT-01_variance_history_vs_codes_error.png)

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

![V-GT-01 cost](assets/V-GT-01/V-GT-01_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ThermaKin2Ds 7 | h, Δt | 5.566e-05 | 0.05031 | 0.02707 | 320, 10× tol | 0.01542 | 1.75 |
| ThermaKin2Ds 7 | h, 2Δt | 5.54e-05 | 0.02707 | 0.02707 | 320, 10× tol | 0.01542 | 1.75 |
| ThermaKin2Ds 7 | h, Δt/2 | 5.573e-05 | 0.09509 | 0.02707 | 320, 10× tol | 0.01542 | 1.75 |
| ThermaKin2Ds 7 | h, Δt/4 | 5.575e-05 | 0.1851 | 0.02707 | 320, 10× tol | 0.01542 | 1.75 |
| ThermaKin2Ds 7 | h/2, Δt/4 | 1.396e-05 | 0.3681 | 0.1865 | 640, 1000× tol | 0.02962 | 6.3 |
| ThermaKin2Ds 7 | h/2, Δt/2 | 1.395e-05 | 0.1865 | 0.1865 | 640, 1000× tol | 0.02962 | 6.3 |

Solution overlays including the other codes' points: [xi_profiles_vs_codes](assets/V-GT-01/V-GT-01_xi_profiles_vs_codes.png), [variance_history_vs_codes](assets/V-GT-01/V-GT-01_variance_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 40 | 0.025 | 0.002544 | 0.003984 | 0.01214 |
| 80 | 0.0125 | 0.004499 | 0.0009693 | 0.003074 |
| 160 | 0.00625 | 0.00836 | 0.0002408 | 0.0007709 |
| 320 | 0.003125 | 0.01531 | 6.011e-05 | 0.000193 |
| 640 | 0.001563 | 0.02983 | 1.502e-05 | 4.826e-05 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 2.011 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
