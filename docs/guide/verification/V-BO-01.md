# V-BO-01 — PMMA-on-Kaowool burnout under ALE (parked depletion)

**Tags:** `ale`, `burnout`, `depletion`, `beer_lambert`, `coupled`, `performance`

**References:**

- MaCFP-2 PMMA gasification case (NIST), 50 kW/m²
- Gpyro V2 PMMA deck property fits

## Problem statement

A 5.9 mm PMMA slab on 12 mm Kaowool under 50 kW/m² burns out completely:
every PMMA cell shrinks to the depletion threshold and is parked, the
exposed boundary walks down the stack, and the run ends on the bare
substrate. The deck is the MaCFP PMMA gasification case with the Gpyro V2
property fits (power-law k(T), c(T) above 394 K — positive on the whole
property bracket) and in-depth (Beer–Lambert) absorption.

The run must burn out completely, lose exactly the mass of the PMMA layer,
match the reference peak mass loss rate, its time and the burnout time,
close the mass ledger, and converge the surface temperature at every
evaluation, within a wall-time budget at 0.098 mm cells.

### Quantities of interest (n = 120)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| complete burnout: live thickness at t_end is the substrate | 0.012 | 0.012 | 3.333e-06 | 1e-05 | yes | every PMMA cell parked; 12.00 ± 0.01 mm |
| mass loss of the gasified PMMA layer | 69.92 | 69.91 | 0.005843 | 0.3 | yes | percent of the initial mass; exact = ρ_PMMA L_PMMA / m₀ |
| peak MLR | 32.47 | 32.39 | 0.002449 | 0.03 | yes | g m⁻² s⁻¹, gas generation rate |
| mass ledger closure | 0.002557 | 0 | 0.002557 | 0.001 | **no** | total-mass drop vs ∫ ṁ″ dt over the run |
| surface-Newton warnings | 0 | 0 | 0 | 0.5 | yes | boundary temperature solve converged at every evaluation |
| time of peak MLR | 335.5 | 335 | 0.001493 | 0.03 | yes | s |
| burnout time (live thickness < 12.02 mm) | 363.5 | 363.5 | 0 | 0.01 | yes | s |

![V-BO-01 mlr_history](assets/V-BO-01/V-BO-01_mlr_history.png)

![V-BO-01 thickness_history](assets/V-BO-01/V-BO-01_thickness_history.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| use_ale | true |
| handle_depletion | true |
| depletion_threshold | 0.05 |
| min_cells | 2 |
| min_thickness | 0.001 |
| radiation_model | BEER_LAMBERT |
| abstol | 1.0e-8 |
| reltol | 1.0e-6 |

