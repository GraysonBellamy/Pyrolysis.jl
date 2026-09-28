# V-BO-05 — PMMA-on-Kaowool burnout on a fixed mesh under in-depth absorption stops at the first interior void

**Tags:** `eulerian`, `burnout`, `depletion`, `beer_lambert`, `coupled`

**References:**

- MaCFP-2 PMMA gasification case (NIST), 50 kW/m²

## Problem statement

The V-BO-01 deck, in-depth absorption included, on a fixed mesh. The
temperature maximum of the PMMA lies below the surface, so the surface cell
reacts slowest: it persists as a low-density skin while the cells beneath it
empty, and the layer gasifies from within. A void cell conducts with the gas
conductivity and absorbs nothing, and nothing in the model carries heat
across a void, so the flux transmitted through the void column would heat
the substrate beneath it far outside the property bracket.

The depletion callback stops a fixed-mesh solve when a live cell below the
top live cell becomes a void (matrix volume fraction below 1e-4). The run
must stop between 330 and 345 s with the interior-void warning, before any
park, with the maximum live-cell temperature below 1000 K. The solve ends
with retcode Terminated.

### Quantities of interest (n = 120)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| stop time at the first interior void | 336.7 | 337.5 | 0.7651 | 7.5 | yes | s; window 330–345 s (observed 337.0): the lowest PMMA cell reaches matrix volume fraction < MIXING_BLEND_LO under the skin |
| interior-void warning | 1 | 1 | 0 | 0.5 | yes | the stop names the cause and the ALE mesh |
| park events before the stop | 0 | 0 | 0 | 0.5 | yes | the surface cell persists as a skin while the layer gasifies from within |
| maximum live-cell temperature up to the stop | 890.8 | — | 890.8 | 1000 | yes | K over every accepted step (observed 912); the solve stops before the void column runs away |

![V-BO-05 mlr_history](assets/V-BO-05/V-BO-05_mlr_history.png)

![V-BO-05 max_temperature_history](assets/V-BO-05/V-BO-05_max_temperature_history.png)

![V-BO-05 pmma_profile](assets/V-BO-05/V-BO-05_pmma_profile.png)

![V-BO-05 temperature_profile](assets/V-BO-05/V-BO-05_temperature_profile.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| use_ale | false |
| handle_depletion | true |
| depletion_threshold | 0.05 |
| min_cells | 2 |
| radiation_model | BEER_LAMBERT |
| abstol | 1.0e-8 |
| reltol | 1.0e-6 |

