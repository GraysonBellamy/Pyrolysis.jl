# V-BO-02 — PMMA-on-Kaowool burnout on a fixed mesh under surface absorption (top-cell parking)

**Tags:** `eulerian`, `burnout`, `depletion`, `surface_absorption`, `coupled`

**References:**

- MaCFP-2 PMMA gasification case (NIST), 50 kW/m²

## Problem statement

The V-BO-01 deck on a fixed mesh, with the incident flux absorbed at the
surface. A fixed mesh keeps every cell's volume, so a gasifying cell empties
instead of thinning, and only the top live cell parks: once its condensed
concentration falls below 5 % of its initial value it merges into the cell
beneath, which becomes the surface. With the flux absorbed at the surface,
the surface cell is the hottest and empties first, so the layer burns top
down and no cell below the surface empties before it.

The run must park all 60 PMMA cells, lose exactly the mass of the PMMA
layer, end on the bare substrate, close the mass ledger, keep the maximum
live-cell temperature below 1000 K, and take no step below the post-park
reset. The peak mass loss rate is informational: surface absorption is a
different problem from the in-depth reference of V-BO-01.

### Quantities of interest (n = 120)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| complete burnout: live thickness at t_end is the substrate | 0.012 | 0.012 | 3.333e-06 | 1e-05 | yes | every PMMA cell parked; 12.00 ± 0.01 mm |
| mass loss of the gasified PMMA layer | 69.92 | 69.91 | 0.005843 | 0.3 | yes | percent of the initial mass; exact = ρ_PMMA L_PMMA / m₀ |
| peak MLR | 24.25 | 24.61 | 0.01478 | 0.25 | yes | g m⁻² s⁻¹; surface absorption on a fixed mesh, a different problem from the in-depth reference (informational) |
| mass ledger closure | 0.0008933 | 0 | 0.0008933 | 0.001 | yes | total-mass drop vs ∫ ṁ″ dt over the run |
| surface-Newton warnings | 0 | 0 | 0 | 0.5 | yes | boundary temperature solve converged at every evaluation |
| park events | 120 | 120 | 0 | 0.5 | yes | every PMMA cell parks as the top live cell, top down |
| maximum live-cell temperature | 917 | — | 917 | 1000 | yes | K over every accepted step (observed 915); the surface cell is the hottest and empties first, so no interior void forms |
| smallest accepted step | 1e-06 | 1e-06 | 2.525e-09 | 0.001 | yes | s; no step below the post-park reset (observed 1e-6 s) |
| Newton failures | 102 | 0 | 102 | 100 | **no** | budget at the CI resolution (60 PMMA cells); observed 56, about one per park interval |
| accepted steps | 4198 | 0 | 4198 | 2000 | **no** | budget at the CI resolution (60 PMMA cells); observed 1572, about 26 per park interval |

![V-BO-02 mlr_history](assets/V-BO-02/V-BO-02_mlr_history.png)

![V-BO-02 thickness_history](assets/V-BO-02/V-BO-02_thickness_history.png)

![V-BO-02 max_temperature_history](assets/V-BO-02/V-BO-02_max_temperature_history.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| use_ale | false |
| handle_depletion | true |
| depletion_threshold | 0.05 |
| min_cells | 2 |
| radiation_model | SURFACE_ABSORPTION |
| abstol | 1.0e-8 |
| reltol | 1.0e-6 |

