# V-GT-08 — imposed-inflow Darcy slab (bottom to top)

**Tags:** `gas_transport`, `darcy`, `mass_flux_bc`, `analytical`

**References:**

- [Gpyro Verification and Validation Manual v1.0 (2026), §4.3, p. 29, Eqs. (4.5)–(4.6), Appendix A.2](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/documentation/gpyro_verification_manual.pdf#page=31)

## Problem statement

Adapted from Gpyro Verification and Validation Manual v1.0 (2026), §4.3,
p. 29, Eqs. (4.5)–(4.6) and Appendix A.2: a 0.1 m isothermal porous
slab with permeability 5e-10 m², T = 300 K, an outlet at 1 bar and an
imposed incoming pure-gas flux of 10 kg/(m² s). Use constant dynamic
viscosity μ = 2e-5 Pa s instead of Gpyro's Chapman–Enskog correlation.

Darcy's law and the ideal-gas EOS give J = -(κM/μRT) P dP/dz, so
P² = Pout² + 2μRT J x/(κM), where x is distance inward from the outlet.
There is no reaction or diffusive wall transport. The initial pressure is
uniform at 1 bar; the solve must relax to this independently derived profile.
Member a feeds the bottom and vents the top; b mirrors both boundaries.
MassFluxBC is outward-positive, so the inlet flux is -J in both cases.
An adiabatic matrix, gas and incoming gas all at 300 K preserve T exactly.

## V-GT-08a — imposed-inflow Darcy slab (bottom to top)

_The steady mass-balance plot shows the saved samples at 0.5 and 1 s. Markers should lie on the reference lines: inlet −10, outlet +10, and their sum zero kg/(m² s). Flux is outward-positive, so negative means gas enters the slab. The initial sample at t = 0 is omitted from this steady comparison: its outlet flux is zero while pressure builds._

### Quantities of interest (n = 80)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| linear P² profile / pressure-squared rise | 0.0008129 | — | 0.0008129 | 0.0025 | yes | Gpyro Eq. (4.5), μRT/M in place of α; n=80 observed 8.13e-4, ~3× margin |
| prescribed signed inlet flux | -10 | -10 | 0 | 1e-12 | yes |  |
| outlet flux vs imposed inflow | 10 | 10 | 1.799e-14 | 1e-08 | yes | independent steady mass balance after pressure relaxation |
| steady signed boundary balance | 1.794e-13 | — | 1.794e-14 | 1e-08 | yes |  |
| gas inventory steadiness | 0 | — | 0 | 1e-08 | yes | relative inventory change between 0.5 and 1 s |
| isothermal invariant | 2.842e-13 | — | 2.842e-13 | 1e-08 | yes |  |

![V-GT-08a pressure](assets/V-GT-08a/V-GT-08a_pressure.png)

![V-GT-08a pressure_error](assets/V-GT-08a/V-GT-08a_pressure_error.png)

![V-GT-08a boundary_flux](assets/V-GT-08a/V-GT-08a_boundary_flux.png)

![V-GT-08a boundary_flux_error](assets/V-GT-08a/V-GT-08a_boundary_flux_error.png)

![V-GT-08a convergence](assets/V-GT-08a/V-GT-08a_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-GT-08a/V-GT-08a_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Gpyro V2 (`bff0126c68f2`) uses 81 boundary-centered nodes (80 intervals) and Δt = 0.005 s, with temperature held at 300 K. Its molecular collision parameters match the Darcy coefficient νP = μRT/M used in this problem. Native gauge pressure is converted to absolute pressure, depth to z = L − depth, and boundary fluxes to outward-positive kg/(m² s). Both codes are compared with the same steady solution.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-GT-08a pressure_vs_codes_error](assets/V-GT-08a/V-GT-08a_pressure_vs_codes_error.png)

![V-GT-08a boundary_flux_vs_codes_error](assets/V-GT-08a/V-GT-08a_boundary_flux_vs_codes_error.png)

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

![V-GT-08a cost](assets/V-GT-08a/V-GT-08a_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Gpyro V2 | h, Δt | 43.45 | 0.009489 | 0.005602 | 40, 1000× tol | 0.004261 | 1.31 |
| Gpyro V2 | h, 2Δt | 43.45 | 0.006626 | 0.005602 | 40, 1000× tol | 0.004261 | 1.31 |
| Gpyro V2 | h, 4Δt | 43.45 | 0.005602 | 0.005602 | 40, 1000× tol | 0.004261 | 1.31 |
| Gpyro V2 | h, Δt/2 | 43.45 | 0.01428 | 0.005602 | 40, 1000× tol | 0.004261 | 1.31 |
| Gpyro V2 | h, Δt/4 | 43.45 | 0.02496 | 0.005602 | 40, 1000× tol | 0.004261 | 1.31 |
| Gpyro V2 | h/2, Δt/4 | 21.69 | 0.04079 | 0.02311 | 80, 1000× tol | 0.008389 | 2.75 |
| Gpyro V2 | h/2, Δt/2 | 21.69 | 0.02311 | 0.02311 | 80, 1000× tol | 0.008389 | 2.75 |

Solution overlays including the other codes' points: [pressure_vs_codes](assets/V-GT-08a/V-GT-08a_pressure_vs_codes.png), [boundary_flux_vs_codes](assets/V-GT-08a/V-GT-08a_boundary_flux_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | P2 |
| --- | --- | --- | --- |
| 20 | 0.05 | 0.00514 | 0.00319 |
| 40 | 0.025 | 0.008976 | 0.001616 |
| 80 | 0.0125 | 0.01798 | 0.0008129 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 0.9861 (P2), expected 1.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| reltol | 1.0e-9 |
| abstol | 1.0e-11 |
| pressure_tolerance | 0.001 |

## V-GT-08b — imposed-inflow Darcy slab (top to bottom)

_The steady mass-balance plot shows the saved samples at 0.5 and 1 s. Markers should lie on the reference lines: inlet −10, outlet +10, and their sum zero kg/(m² s). Flux is outward-positive, so negative means gas enters the slab. The initial sample at t = 0 is omitted from this steady comparison: its outlet flux is zero while pressure builds._

### Quantities of interest (n = 80)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| linear P² profile / pressure-squared rise | 0.0008129 | — | 0.0008129 | 0.0025 | yes | Gpyro Eq. (4.5), μRT/M in place of α; n=80 observed 8.13e-4, ~3× margin |
| prescribed signed inlet flux | -10 | -10 | 0 | 1e-12 | yes |  |
| outlet flux vs imposed inflow | 10 | 10 | 1.843e-14 | 1e-08 | yes | independent steady mass balance after pressure relaxation |
| steady signed boundary balance | 1.847e-13 | — | 1.847e-14 | 1e-08 | yes |  |
| gas inventory steadiness | 1.721e-16 | — | 1.721e-16 | 1e-08 | yes | relative inventory change between 0.5 and 1 s |
| isothermal invariant | 2.274e-13 | — | 2.274e-13 | 1e-08 | yes |  |

![V-GT-08b pressure](assets/V-GT-08b/V-GT-08b_pressure.png)

![V-GT-08b pressure_error](assets/V-GT-08b/V-GT-08b_pressure_error.png)

![V-GT-08b boundary_flux](assets/V-GT-08b/V-GT-08b_boundary_flux.png)

![V-GT-08b boundary_flux_error](assets/V-GT-08b/V-GT-08b_boundary_flux_error.png)

![V-GT-08b convergence](assets/V-GT-08b/V-GT-08b_convergence.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Costs are compared at matched accuracy only.

Matched native-input runs at seven configurations (coarser and finer time steps at fixed mesh, and a finer mesh) are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at every configuration](assets/V-GT-08b/V-GT-08b_comparison_accuracy.csv) are retained separately, and the same runs are timed for the cost comparison.

Gpyro V2 (`bff0126c68f2`) uses 81 boundary-centered nodes (80 intervals) and Δt = 0.005 s, with temperature held at 300 K. Its molecular collision parameters match the Darcy coefficient νP = μRT/M used in this problem. Native gauge pressure is converted to absolute pressure, depth to z = L − depth, and boundary fluxes to outward-positive kg/(m² s). Both codes are compared with the same steady solution.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-GT-08b pressure_vs_codes_error](assets/V-GT-08b/V-GT-08b_pressure_vs_codes_error.png)

![V-GT-08b boundary_flux_vs_codes_error](assets/V-GT-08b/V-GT-08b_boundary_flux_vs_codes_error.png)

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

![V-GT-08b cost](assets/V-GT-08b/V-GT-08b_cost.png)

For each external mesh and step: the external code's cheapest run at least as accurate (often itself; a coarser or finer run of the same code when that costs less), this code's cheapest configuration at least as accurate on every scored reference (its cells, and its tolerance multiplier when looser than the report's), and the ratio of the two. A ratio marked ≥ uses this code's coarsest configuration, fewest cells at the loosest tolerance, whose accuracy already exceeds the external run's.

| External run | Mesh, step | External error | External time (s) | Cheapest external at this accuracy (s) | Our cells, tolerance | Our time (s) | Cheapest external/our time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Gpyro V2 | h, Δt | 43.21 | 0.009623 | 0.005453 | 40, 1000× tol | 0.004057 | 1.34 |
| Gpyro V2 | h, 2Δt | 43.21 | 0.007086 | 0.005453 | 40, 1000× tol | 0.004057 | 1.34 |
| Gpyro V2 | h, 4Δt | 43.21 | 0.005453 | 0.005453 | 40, 1000× tol | 0.004057 | 1.34 |
| Gpyro V2 | h, Δt/2 | 43.21 | 0.01465 | 0.005453 | 40, 1000× tol | 0.004057 | 1.34 |
| Gpyro V2 | h, Δt/4 | 43.21 | 0.02475 | 0.005453 | 40, 1000× tol | 0.004057 | 1.34 |
| Gpyro V2 | h/2, Δt/4 | 21.63 | 0.03976 | 0.02287 | 80, 1000× tol | 0.008018 | 2.85 |
| Gpyro V2 | h/2, Δt/2 | 21.63 | 0.02287 | 0.02287 | 80, 1000× tol | 0.008018 | 2.85 |

Solution overlays including the other codes' points: [pressure_vs_codes](assets/V-GT-08b/V-GT-08b_pressure_vs_codes.png), [boundary_flux_vs_codes](assets/V-GT-08b/V-GT-08b_boundary_flux_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | P2 |
| --- | --- | --- | --- |
| 20 | 0.05 | 0.004796 | 0.00319 |
| 40 | 0.025 | 0.008986 | 0.001616 |
| 80 | 0.0125 | 0.01833 | 0.0008129 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order 0.9861 (P2), expected 1.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| reltol | 1.0e-9 |
| abstol | 1.0e-11 |
| pressure_tolerance | 0.001 |
