# V-BO-04 — Uniform zero-order burnout: exact thickness through park cascades

**Tags:** `ale`, `burnout`, `depletion`, `kinetics`, `exact`

**References:**

- closed form L(t) = L₀ exp(−r t/ρ_s)

## Problem statement

An adiabatic slab at uniform T₀ with one zero-order reaction solid → gas,
γ_gas = 0 and no heat of reaction. The temperature stays at T₀ exactly (no
heat source, no gradient), so the volumetric rate r = A·exp(−E/RT₀) is a
constant and every cell loses the same fraction r/ρ_s of its current
volume per unit time. Under ALE the solid concentration stays at ρ_s
(dilation removes volume, not concentration) and the exact thickness is

  L(t) = L₀ · exp(−r t / ρ_s).

Every cell reaches the depletion threshold at the same instant, so the
parking machinery meets its worst case: a cascade of park events, one per
step, surface-first, halving the live cell count at t_k = ln(2^k·20) ρ_s/r.
L(t) and the solid mass must follow the exact law through every cascade,
and ξ_s must stay at ρ_s to roundoff.

### Quantities of interest (n = 40)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| thickness follows L₀ exp(−r t/ρ) through the park cascades | 8.23e-05 | 8.23e-05 | 1.932e-07 | 0.0001 | yes | max relative error over all output times |
| solid concentration stays at ρ_s (dilation removes volume only) | 1000 | 1000 | 2.274e-16 | 1e-08 | yes | max \|ξ_s − ρ_s\|/ρ_s over live cells and times |
| isothermal: T stays at T₀ | 600 | 600 | 1.137e-13 | 1e-06 | yes | K; no heat of reaction, no gradient |
| park cascade: 40 → 5 live cells | 5 | 5 | 0 | 0.5 | yes | 35 park events; every cell hits the threshold at once |

![V-BO-04 thickness_history](assets/V-BO-04/V-BO-04_thickness_history.png)

![V-BO-04 thickness_history_error](assets/V-BO-04/V-BO-04_thickness_history_error.png)

![V-BO-04 live_cells_history](assets/V-BO-04/V-BO-04_live_cells_history.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| use_ale | true |
| handle_depletion | true |
| depletion_threshold | 0.05 |
| min_cells | 2 |
| min_thickness | 1.0e-6 |
| abstol | 1.0e-9 |
| reltol | 1.0e-7 |

