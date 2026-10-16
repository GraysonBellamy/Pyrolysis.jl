# White Pine UMD CAPA N₂ test measurements

Canonical record of pre- and post-test measurements provided 2026-04-30 for
White Pine cone-calorimeter tests in nitrogen. Used as calibration anchors for
the variable-cross-section model and the inverse-fit conductivities.

Sample provenance: MaCFP / UMD CAPA chamber, N₂ atmosphere, dry-conditioned.

## 1. Pre-test data (10 replicates, 3 fluxes)

| File | flux<br/>(kW/m²) | cone T<br/>(K) | D₀<br/>(mm) | L₀<br/>(mm) | m₀<br/>(g) | Δm<br/>(g) | col-8<br/>(K) | col-9<br/>(K) |
|---|---|---|---|---|---|---|---|---|
| UMD_Wood_CAPA_N2_30kW_R1 | 30 | 873.4 | 69.51 | 8.50 | 11.43 | 6.26 | 392.4 | 302.9 |
| UMD_Wood_CAPA_N2_30kW_R2 | 30 | 870.9 | 69.42 | 8.52 | 12.17 | 6.39 | 394.6 | 303.8 |
| UMD_Wood_CAPA_N2_30kW_R3 | 30 | 873.1 | 69.60 | 8.49 | 11.16 | 5.55 | 390.9 | 301.9 |
| UMD_Wood_CAPA_N2_50kW_R1 | 50 | 1002.9 | 69.48 | 8.45 | 12.60 | 9.12 | 447.7 | 307.9 |
| UMD_Wood_CAPA_N2_50kW_R2 | 50 | 1001.9 | 69.48 | 8.52 | 11.97 | 8.60 | 447.7 | 307.6 |
| UMD_Wood_CAPA_N2_50kW_R3 | 50 | 1004.5 | 69.27 | 8.50 | 12.10 | 8.97 | 449.7 | 308.1 |
| UMD_Wood_CAPA_N2_50kW_R4 | 50 | 1004.5 | 69.50 | 8.49 | 11.05 | 8.00 | 449.3 | 307.9 |
| UMD_Wood_CAPA_N2_70kW_R1 | 70 | 1100.2 | 69.62 | 8.50 | 11.48 | 8.69 | 497.0 | 311.7 |
| UMD_Wood_CAPA_N2_70kW_R2 | 70 | 1100.8 | 69.53 | 8.50 | 10.90 | 8.39 | 498.8 | 311.7 |
| UMD_Wood_CAPA_N2_70kW_R3 | 70 | 1100.3 | 69.63 | 8.51 | 11.26 | 8.71 | 497.5 | 312.1 |

Columns 8 and 9 (labels TBD by user — best guesses: peak back temperature
during test (K), and initial sample temperature at start of test (K)).

## 2. Post-test thickness L_f (mm)

Replicate-level measurements:

| flux | replicates | n | mean | spread |
|---|---|---|---|---|
| 30 kW | TBD | 3 | **7.60** | TBD |
| 50 kW | TBD | 4 | **8.80** | TBD |
| 70 kW | TBD | 3 | **9.30** | TBD |

## 3. Post-test diameter D_f (mm)

| flux | D_f | n |
|---|---|---|
| 30 kW | 65.00 | 1 |
| 50 kW | 64.40 | 1 |
| 70 kW | 64.40 | 1 (R1 anchor) |

(More replicates of D_f post-test would be valuable; currently the radial
channel is identified from a single sample per flux.)

## 4. Per-flux pre-test diameter averages

| flux | n | mean D₀ (mm) | A₀ (mm²) |
|---|---|---|---|
| 30 | 3 | 69.51 | 3796.0 |
| 50 | 4 | 69.43 | 3787.3 |
| 70 | 3 | 69.59 | 3804.8 |
| all | 10 | **69.50** | **3794.9** |

## 5. Derived per-flux means

Using replicate means and per-flux pre-test averages:

| flux | m_f/m_0<br/>(meas) | A_f/A_0<br/>(meas) | L_f/L_0<br/>(meas) | V_f/V_0<br/>(meas) | implied ρ_f<br/>(kg/m³) |
|---|---|---|---|---|---|
| 30 | 0.477 ± 0.026 | 0.875 | 0.947 | 0.829 | **207** |
| 50 | 0.273 ± 0.010 | 0.860 | 1.065 | 0.916 | **107** |
| 70 | 0.233 ± 0.008 | 0.856 | 1.133 | 0.971 | **86** |

m_f/m_0 from `(m₀ − Δm)/m₀` averaged over replicates:
- 30 kW: 0.452, 0.475, 0.503 → 0.477
- 50 kW: 0.276, 0.282, 0.259, 0.276 → 0.273
- 70 kW: 0.243, 0.230, 0.227 → 0.233

Implied ρ_f assumes L₀ = 8.50 mm and uses ρ₀ = 360 kg/m³.

## 6. Calibration anchors used in `wp_properties.jl`

The simulation pulls these from the per-flux means above:

```
WP_D0          = 0.06950 m       (overall pre-test mean)
WP_A0          = π·(WP_D0/2)²    ≈ 3.795 × 10⁻³ m²
WP_D0_30       = 0.06951 m       (per-flux mean, n=3)
WP_D0_50       = 0.06943 m       (per-flux mean, n=4)
WP_D0_70       = 0.06959 m       (per-flux mean, n=3)

WP_DF_30_MM    = 65.00           (1 sample)
WP_DF_50_MM    = 64.40           (1 sample)
WP_DF_70_MM    = 64.40           (R1 original anchor)

WP_LF_30_MM    = 7.60            (mean of 3 replicates)
WP_LF_50_MM    = 8.80            (mean of 4 replicates)
WP_LF_70_MM    = 9.30            (mean of 3 replicates)

WP_MASS_RATIO_30 = 0.477         (mean of 3 replicates, ±0.026)
WP_MASS_RATIO_50 = 0.273         (mean of 4 replicates, ±0.010)
WP_MASS_RATIO_70 = 0.233         (mean of 3 replicates, ±0.008)
```

## 7. Calibration constants derived from this data

### Radial shrinkage law
`A/A_0 = 1 − (1 − s_r²) · χ̄^p`, with **p = 1**, **s_r² = 0.807** (LSQ over all
three flux anchors). The data does not require saturation in the χ̄ range
observed; an exponential-saturation form was tested and degenerated to the
linear limit.

### Fixed dilatometry component densities
The inverse-analysis component densities come from dilatometry (the scaled
density column below); `wp_properties.jl` holds the values the model uses:

| component | T_peak (°C) | Y_peak | ρ_raw exp (kg/m³) | ρ_scaled (kg/m³) |
|---|---:|---:|---:|---:|
| virgin | 27 | 1.00 | 362 | 360 |
| int1 | 283 | 0.71 | 330 | 321 |
| int2 | 337 | 0.58 | 280 | 259 |
| int3 | 386 | 0.20 | 168 | 122 |
| int4 | 476 | 0.16 | 142 | 90 |
| char | 812 | 0.15 | 146 | 95 |

Density sensitivity means rerunning with an edited component table.

### Pre-test sample geometry
A_0 = 3.795 × 10⁻³ m² (overall mean), or per-flux from the table above.

## 8. Open questions / outstanding measurements

1. **Columns 8 and 9 of section 1.** Best guesses are peak back temperature
   (K) and initial sample temperature (K). Confirm.
2. **Post-test diameter at 30 / 50 kW** — only 1 replicate currently.
   Replicate spread on D_f would tighten the radial calibration.
3. **Post-test envelope density at 30 and 50 kW** (= m_f / V_f). Currently
   inferred via section 5 from m_f/m_0 and L_f/D_f. Direct measurement
   (weigh + caliper) would decouple from any kinetic-error feedback.
4. **Char swelling / void volume submodel.** The 70 kW residue density
   (86 kg/m³) is below the int4 anchor (91), implying ~5% void / swelling.
   No model for this yet.

## 9. Cross-references

- Calibration code: [`examples/inverse_analysis/wp_properties.jl`](../wp_properties.jl)
- Forward driver: [`examples/inverse_analysis/wp_forward.jl`](../wp_forward.jl)
- Validation driver: [`examples/inverse_analysis/validate_thickness.jl`](../validate_thickness.jl)
