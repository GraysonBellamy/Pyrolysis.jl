# V-HC-03 — finite slab, constant flux, adiabatic back

**Tags:** `heat_conduction`, `heat_flux_bc`, `energy_balance`

**References:**

- Carslaw & Jaeger (1959) §3.8 (constant-flux slab series)

## Problem statement

Complements the semi-infinite cases at LONG times: the thermal wave reaches
the back face — the regime cone-calorimeter simulations actually live in.
Constant q″ at the top face, insulated back; exact solution
T = T₀ + (q″L/k)·φ(x̃, Fo) with φ from exact/slab_series.jl (C&J §3.8),
x̃ = z/L from the adiabatic back face.

Two parameter-free exact targets beyond the profiles:
- front–back plateau: ΔT → q″L/2k for Fo ≫ 1 (100 K here),
- global energy: mean temperature rises exactly at q″/(ρcL) (0.25 K/s here).

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T-profile error, Fo=0.1 | 0.0001225 | — | 0.0001225 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| L∞ T-profile error, Fo=0.5 | 9.172e-05 | — | 9.172e-05 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| L∞ T-profile error, Fo=1.5 | 8.139e-05 | — | 8.139e-05 | 0.006 | yes | observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5 |
| front–back ΔT plateau, Fo=1.5 | 99.69 | 99.69 | 1.155e-08 | 1e-05 | yes | steady-shape plateau q″L/2k = 100 K between faces |
| mean-temperature rise (global energy) | 200 | 200 | 3.411e-13 | 1e-05 | yes | exact constant q″/(ρcL) = 0.25 K/s; integrator-tolerance limited |

![V-HC-03 T_profiles](assets/V-HC-03/V-HC-03_T_profiles.png)

![V-HC-03 T_profiles_error](assets/V-HC-03/V-HC-03_T_profiles_error.png)

![V-HC-03 front_back_dT_history](assets/V-HC-03/V-HC-03_front_back_dT_history.png)

![V-HC-03 front_back_dT_history_error](assets/V-HC-03/V-HC-03_front_back_dT_history_error.png)

![V-HC-03 mean_T_history](assets/V-HC-03/V-HC-03_mean_T_history.png)

![V-HC-03 mean_T_history_error](assets/V-HC-03/V-HC-03_mean_T_history_error.png)

![V-HC-03 convergence](assets/V-HC-03/V-HC-03_convergence.png)

### Comparison with other codes

The same case was solved with Gpyro 0.8200 at 201 nodes (Δz=0.1 mm); decks, outputs, and
run provenance are committed under `test/verification/reference/`.
Each code's error against the same exact solution is drawn below
on a log scale, muted gray behind this solver's series — the
signed linear-scale panel above shows where the error lives, this
one compares magnitudes across codes.

![V-HC-03 T_profiles_vs_codes_error](assets/V-HC-03/V-HC-03_T_profiles_vs_codes_error.png)

![V-HC-03 front_back_dT_history_vs_codes_error](assets/V-HC-03/V-HC-03_front_back_dT_history_vs_codes_error.png)

Wall times at every ladder rung against the reference runs.
Resolutions and simulated spans differ where noted (details in
`reference/timings.csv`), so cross-code timings are indicative
rather than a controlled benchmark; rungs at a matched
resolution are directly comparable.

![V-HC-03 solve_time](assets/V-HC-03/V-HC-03_solve_time.png)

Solution overlays including the other codes' points: [T_profiles_vs_codes](assets/V-HC-03/V-HC-03_T_profiles_vs_codes.png), [front_back_dT_history_vs_codes](assets/V-HC-03/V-HC-03_front_back_dT_history_vs_codes.png).

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 20 | 0.05 | 0.0152 | 0.02373 | 0.0314 |
| 40 | 0.025 | 0.02775 | 0.005935 | 0.007862 |
| 80 | 0.0125 | 0.05262 | 0.001484 | 0.001966 |
| 160 | 0.00625 | 0.1024 | 0.000371 | 0.0004916 |
| 200 | 0.005 | 0.1186 | 0.0002372 | 0.0003141 |
| 320 | 0.003125 | 0.19 | 9.254e-05 | 0.0001225 |

Observed order 2 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

