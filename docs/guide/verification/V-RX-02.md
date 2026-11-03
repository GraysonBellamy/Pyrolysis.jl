# V-RX-02 — zero-order decay (exact burnout time)

**Tags:** `kinetics`, `reaction_order`, `isothermal`

**References:**

- Khawam & Flanagan (2006), J. Phys. Chem. B 110:17315

## Problem statement

Khawam & Flanagan (2006) exact isothermal integrations (FDS verifies its
`part_baking_soda` case against the same forms). Sub-cases:

  a) n = 0    — linear decay, exact finite burnout time. Compared only up to
                0.9·t_b with the depletion limiter DISABLED: the power is
                taken on the finite base of the conversion, so at ξ = 0 the
                n = 0 rate stays k ξ₀ and would integrate ξ negative
                (`DepletionLimiter(enabled = false)`).
  b) n = 2    — single-species second-order decay (a C+C kinetic analog,
                not a reproduced ThermaKin run):
                ξ = ξ₀/(1 + kt).
  c) n = 2/3  — contracting-sphere R3 model with FDS's NaHCO₃ parameters
                A = 3.4e11 1/s, E = 103 kJ/mol at T = 450 K; the marker
                is m^(1/3) linear in t.
  d) n = 0 with the limiter ENABLED (default) — non-negativity guard: run to
                2× burnout and require 0 ≤ x ≤ limiter threshold at the end.
                This is the case the limiter exists for.

The rate is k·ξ₀·xⁿ on the conversion x = ξ/ξ₀, k in s⁻¹ for every order, so
published constants for normalized m = ξ/ξ₀ are the solver's as they stand
(see exact/kinetics.jl header).

## V-RX-02a — zero-order decay (exact burnout time)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=0.0) | 7.39e-15 | — | 7.39e-15 | 1e-12 | yes | linear decay; compared only to 0.9·t_b, limiter off |
| total-mass conservation | 7.105e-16 | — | 7.105e-16 | 1e-13 | yes | A→B closed system |

![V-RX-02a xi_histories](assets/V-RX-02a/V-RX-02a_xi_histories.png)

![V-RX-02a xi_histories_error](assets/V-RX-02a/V-RX-02a_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Runtime comparisons require a shared achieved-accuracy target.

Matched native-input runs, independent time refinement at fixed mesh (dt, dt/2, dt/4), and spatial refinement at fixed dt/4 are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at all retained resolutions](assets/V-RX-02a/V-RX-02a_comparison_accuracy.csv) are retained separately. These runs provide accuracy measurements.

Accuracy comparison only; timing measurements are omitted.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02a xi_histories_vs_codes_error](assets/V-RX-02a/V-RX-02a_xi_histories_vs_codes_error.png)

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02a/V-RX-02a_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02b — second-order decay (ThermaKin C+C analog)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=2.0) | 2.845e-12 | — | 2.845e-12 | 1e-11 | yes | ξ = ξ₀/(1+kt); k in s⁻¹ |
| total-mass conservation | 2.7e-15 | — | 2.7e-15 | 1e-13 | yes | A→B closed system |

![V-RX-02b xi_histories](assets/V-RX-02b/V-RX-02b_xi_histories.png)

![V-RX-02b xi_histories_error](assets/V-RX-02b/V-RX-02b_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Runtime comparisons require a shared achieved-accuracy target.

Matched native-input runs, independent time refinement at fixed mesh (dt, dt/2, dt/4), and spatial refinement at fixed dt/4 are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at all retained resolutions](assets/V-RX-02b/V-RX-02b_comparison_accuracy.csv) are retained separately. These runs provide accuracy measurements.

Accuracy comparison only; timing measurements are omitted.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02b xi_histories_vs_codes_error](assets/V-RX-02b/V-RX-02b_xi_histories_vs_codes_error.png)

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02b/V-RX-02b_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02c — contracting-sphere R3 decay (FDS NaHCO₃, 450 K)

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| max relative ξ_A error (n=0.6666666666666666) | 2.251e-10 | — | 2.251e-10 | 1e-09 | yes | m^(1/3) linear in t; k=0.3765 1/s at 450 K; exercises live Arrhenius (E ≠ 0) at pinned T |
| total-mass conservation | 2.501e-15 | — | 2.501e-15 | 1e-13 | yes | A→B closed system |

![V-RX-02c xi_histories](assets/V-RX-02c/V-RX-02c_xi_histories.png)

![V-RX-02c xi_histories_error](assets/V-RX-02c/V-RX-02c_xi_histories_error.png)

### Comparison with other codes

Physical-problem agreement does not establish equal discretization or cost. See the [mesh and time comparison basis](comparison_coverage.md#comparison-basis-audit) and [per-pair inventory](assets/comparison_basis.csv). Runtime comparisons require a shared achieved-accuracy target.

Matched native-input runs, independent time refinement at fixed mesh (dt, dt/2, dt/4), and spatial refinement at fixed dt/4 are retained under `test/verification/reference/matched/`. See the [coverage and matching notes](comparison_coverage.md). The plotted external series use the base resolution; [analytical errors at all retained resolutions](assets/V-RX-02c/V-RX-02c_comparison_accuracy.csv) are retained separately. These runs provide accuracy measurements.

Accuracy comparison only; timing measurements are omitted.

Error panels show absolute error at each code's output points. Zeros are displayed at 1e-12 on the logarithmic axis; the reference CSV precision limits the discernible error.

![V-RX-02c xi_histories_vs_codes_error](assets/V-RX-02c/V-RX-02c_xi_histories_vs_codes_error.png)

Solution overlays including the other codes' points: [xi_histories_vs_codes](assets/V-RX-02c/V-RX-02c_xi_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |

## V-RX-02d — zero-order decay, depletion limiter ON (guard)

_limiter roll-off characterization_

### Quantities of interest (n = 1)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| non-negativity through 2× burnout | 5.369e-05 | 0 | 0 | 1e-09 | yes | the depletion limiter's job for n = 0 kinetics |
| residual mass within limiter threshold | 5.369e-05 | 0 | 5.369e-05 | 0.1 | yes | quadratic onset near 0 ⇒ algebraic tail; observed 5.4e-5 |

![V-RX-02d xi_histories](assets/V-RX-02d/V-RX-02d_xi_histories.png)

![V-RX-02d xi_histories_error](assets/V-RX-02d/V-RX-02d_xi_histories_error.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
