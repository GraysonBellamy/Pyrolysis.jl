# V-GT-07 — transient binary diffusion through a mass film (h_m = 0.0001)

**Tags:** `gas_transport`, `diffusion`, `mass_film`, `transient`, `analytical`

**References:**

- [Gpyro Verification and Validation Manual v1.0 (2026), §4.4, pp. 30–31, Eqs. (4.7)–(4.11)](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/documentation/gpyro_verification_manual.pdf#page=32)

## Problem statement

Adapted from Gpyro Verification and Validation Manual v1.0 (2026), §4.4,
pp. 30–31, Eqs. (4.7)–(4.11): a 1 m porous slab initially containing
carrier gas, exposed at its top to a reservoir with tracer fraction 0.22;
the bottom is sealed. T = 300 K, P = 1 bar, φ = 0.6, D = 2e-5 m²/s.
Members a/b use h_m = 1e-4 / 1e-2 kg/(m² s), respectively.

Like the upstream input decks, both tracer gases have M = 0.032 kg/mol;
equal molecular masses, heat capacities and wall coefficients preserve
constant total pressure, density and temperature during counterdiffusion.
LumpedTransfer with Λ = φD and κ = 0 realizes the constant-D equation;
Chapman–Enskog diffusivity is replaced by the guide's specified constant.
The zero-total-flow MassTransferFilm law reduces exactly to the Robin BC.

At depth x = L-z, Y = Y∞ Φ(x,t), where Φ is the semi-infinite convective
erfc solution with h/k = b = h_m/(φρ_gD). We stop at 400 s so the back
face is over five diffusion lengths away; the reference tail is checked.
Surface uptake is h_m Y∞ erfcx(b√(Dt)); its time integral independently
predicts the stored tracer mass. All references use original inputs.

## V-GT-07a — transient binary diffusion through a mass film (h_m = 0.0001)

### Quantities of interest (n = 160)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| transient tracer profile / reservoir fraction | 0.0001283 | — | 0.0001283 | 0.0004 | yes | n=160: observed 1.28e-4 / 1.28e-3 for h_m=1e-4 / 1e-2; ~3× margin |
| surface uptake vs analytical flux | 0.0003341 | — | 0.0003341 | 0.0011 | yes | h_m Y∞ erfcx(b√(Dt)); observed 3.34e-4 / 2.45e-3; ~3× margin |
| stored tracer vs integrated analytical influx | 0.0009245 | — | 0.0009245 | 0.003 | yes | integrated exact influx; observed 9.25e-4 / 2.55e-3; ~3× margin |
| total gas inventory conservation | 3.892e-13 | — | 3.892e-13 | 1e-10 | yes | equal-mass binary counterdiffusion has zero net mass influx |
| zero net boundary mass flux | 9.599e-16 | — | 4.363e-11 | 1e-08 | yes | two species counterdiffuse; flux scale h_m Y∞ |
| sealed back face | 0 | — | 0 | 1e-14 | yes |  |
| isobaric invariant | 8.644e-05 | — | 8.644e-10 | 4e-09 | yes | observed relative drift ≤1.05e-9 (0.000105 Pa); pressure tolerance 0.001 Pa |
| isothermal invariant | 0 | — | 0 | 1e-08 | yes |  |
| semi-infinite reference tail at back | 2.444e-16 | — | 2.444e-16 | 1e-12 | yes | bounds the neglected far-boundary penetration |

![V-GT-07a profiles](assets/V-GT-07a/V-GT-07a_profiles.png)

![V-GT-07a profiles_error](assets/V-GT-07a/V-GT-07a_profiles_error.png)

![V-GT-07a uptake](assets/V-GT-07a/V-GT-07a_uptake.png)

![V-GT-07a uptake_error](assets/V-GT-07a/V-GT-07a_uptake_error.png)

![V-GT-07a inventory](assets/V-GT-07a/V-GT-07a_inventory.png)

![V-GT-07a inventory_error](assets/V-GT-07a/V-GT-07a_inventory_error.png)

![V-GT-07a convergence](assets/V-GT-07a/V-GT-07a_convergence.png)

### Comparison with other codes

Gpyro V2 (`bff0126c68f2`) uses 161 boundary-centered nodes (160 intervals) and Δt = 0.05 s. Molecular collision parameters give D = 2e-5 m²/s at 300 K and 1 bar. The far-side reservoir contains only carrier gas so Gpyro's positive floor on the zero film coefficient contributes negligible exchange. Both codes are compared with the analytical Robin solution, with their native gas constants retained. Inputs, raw outputs, refinement checks and extraction are in `test/verification/reference/gpyro/v2_accuracy/`.

Accuracy comparison only; timing measurements are omitted.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-GT-07a profiles_vs_codes_error](assets/V-GT-07a/V-GT-07a_profiles_vs_codes_error.png)

![V-GT-07a uptake_vs_codes_error](assets/V-GT-07a/V-GT-07a_uptake_vs_codes_error.png)

![V-GT-07a inventory_vs_codes_error](assets/V-GT-07a/V-GT-07a_inventory_vs_codes_error.png)

Solution overlays including the other codes' points: [profiles_vs_codes](assets/V-GT-07a/V-GT-07a_profiles_vs_codes.png), [uptake_vs_codes](assets/V-GT-07a/V-GT-07a_uptake_vs_codes.png), [inventory_vs_codes](assets/V-GT-07a/V-GT-07a_inventory_vs_codes.png).

### Convergence

| n_cells | h | Linf |
| --- | --- | --- |
| 40 | 0.025 | 0.001024 |
| 80 | 0.0125 | 0.0002564 |
| 160 | 0.00625 | 6.402e-05 |

Observed order 2 (Linf), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| reltol | 1.0e-9 |
| abstol | 1.0e-11 |
| pressure_tolerance | 0.001 |

## V-GT-07b — transient binary diffusion through a mass film (h_m = 0.01)

### Quantities of interest (n = 160)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| transient tracer profile / reservoir fraction | 0.001282 | — | 0.001282 | 0.004 | yes | n=160: observed 1.28e-4 / 1.28e-3 for h_m=1e-4 / 1e-2; ~3× margin |
| surface uptake vs analytical flux | 0.002445 | — | 0.002445 | 0.0075 | yes | h_m Y∞ erfcx(b√(Dt)); observed 3.34e-4 / 2.45e-3; ~3× margin |
| stored tracer vs integrated analytical influx | 0.002549 | — | 0.002549 | 0.008 | yes | integrated exact influx; observed 9.25e-4 / 2.55e-3; ~3× margin |
| total gas inventory conservation | 2.891e-13 | — | 2.891e-13 | 1e-10 | yes | equal-mass binary counterdiffusion has zero net mass influx |
| zero net boundary mass flux | 8.235e-16 | — | 3.743e-13 | 1e-08 | yes | two species counterdiffuse; flux scale h_m Y∞ |
| sealed back face | 0 | — | 0 | 1e-14 | yes |  |
| isobaric invariant | 0.0001047 | — | 1.047e-09 | 4e-09 | yes | observed relative drift ≤1.05e-9 (0.000105 Pa); pressure tolerance 0.001 Pa |
| isothermal invariant | 0 | — | 0 | 1e-08 | yes |  |
| semi-infinite reference tail at back | 2.427e-15 | — | 2.427e-15 | 1e-12 | yes | bounds the neglected far-boundary penetration |

![V-GT-07b profiles](assets/V-GT-07b/V-GT-07b_profiles.png)

![V-GT-07b profiles_error](assets/V-GT-07b/V-GT-07b_profiles_error.png)

![V-GT-07b uptake](assets/V-GT-07b/V-GT-07b_uptake.png)

![V-GT-07b uptake_error](assets/V-GT-07b/V-GT-07b_uptake_error.png)

![V-GT-07b inventory](assets/V-GT-07b/V-GT-07b_inventory.png)

![V-GT-07b inventory_error](assets/V-GT-07b/V-GT-07b_inventory_error.png)

![V-GT-07b convergence](assets/V-GT-07b/V-GT-07b_convergence.png)

### Comparison with other codes

Gpyro V2 (`bff0126c68f2`) uses 161 boundary-centered nodes (160 intervals) and Δt = 0.05 s. Molecular collision parameters give D = 2e-5 m²/s at 300 K and 1 bar. The far-side reservoir contains only carrier gas so Gpyro's positive floor on the zero film coefficient contributes negligible exchange. Both codes are compared with the analytical Robin solution, with their native gas constants retained. Inputs, raw outputs, refinement checks and extraction are in `test/verification/reference/gpyro/v2_accuracy/`.

Accuracy comparison only; timing measurements are omitted.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-GT-07b profiles_vs_codes_error](assets/V-GT-07b/V-GT-07b_profiles_vs_codes_error.png)

![V-GT-07b uptake_vs_codes_error](assets/V-GT-07b/V-GT-07b_uptake_vs_codes_error.png)

![V-GT-07b inventory_vs_codes_error](assets/V-GT-07b/V-GT-07b_inventory_vs_codes_error.png)

Solution overlays including the other codes' points: [profiles_vs_codes](assets/V-GT-07b/V-GT-07b_profiles_vs_codes.png), [uptake_vs_codes](assets/V-GT-07b/V-GT-07b_uptake_vs_codes.png), [inventory_vs_codes](assets/V-GT-07b/V-GT-07b_inventory_vs_codes.png).

### Convergence

| n_cells | h | Linf |
| --- | --- | --- |
| 40 | 0.025 | 0.002551 |
| 80 | 0.0125 | 0.0006366 |
| 160 | 0.00625 | 0.000159 |

Observed order 2.002 (Linf), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| reltol | 1.0e-9 |
| abstol | 1.0e-11 |
| pressure_tolerance | 0.001 |
