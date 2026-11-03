# MaCFP-2 PMMA community model submissions

Downloaded 2026-07-07 from `MaCFP/matl-db` (branch `master`),
`PMMA/Computational_Results/{TGA_Predictions,Gasification_Predictions}/` —
the idealized "target simulation" predictions submitted by 9 institutions for
the 2021 MaCFP-2 condensed-phase workshop. Row-downsampled 5× (first/last data
rows kept, both header rows preserved) to keep the repo light; fetch the raw
URLs for full resolution.

Scenarios (guidelines v1.3, 2021-02-06): TGA 300→1000 K at 10 and 100 K/min,
5 mg; 1D anaerobic gasification at 10/25/65 kW/m² incident flux, **no
convection**, perfectly insulated back, T₀ = 293 K, thicknesses 6 and 12 mm,
scaled to 10 cm × 10 cm. Gasification columns: time [s], mass [g], back-surface
T [K], top-surface T [K]; TGA: time [s], T [K], mass [mg]. Header text varies
by institution (NIST swaps the name/unit rows; UMD has a trailing comma) —
parse tolerantly.

The authors' [2021 material-property files](https://github.com/MaCFP/matl-db/tree/master/PMMA/Material_Properties/2021)
record calibration methods as well as codes. For example, Aalto records
Gpyro for kinetics and FDS/PROPTI for thermal-property fitting; BUW-FZJ and
NIST record FDS; GIDAZE+ records Gpyro; UCLAN and UMD record ThermaKin;
Sandia records Sierra. A calibration code is not necessarily the generator
of every submitted curve. Do not infer an exact solver/version from an
institution name or these calibration fields alone.
`MaCFP_TGA_10K_*Stats.csv` are the workshop's supplied ensemble statistics.

**Interpretation:** the [authors' computational-results description](https://github.com/MaCFP/matl-db/blob/master/PMMA/Computational_Results/README.md)
explicitly excludes these idealized scenarios from verification and validation.
Each institution uses its own calibrated material model, so the envelope
combines parameter and implementation differences. `MaCFP_PMMA_UMD.json` is
one submitted property set, not a designated universal central model. Using
it here does not isolate solver differences against other submissions with
different inputs. A numerical comparison requires matching the complete
material model, constitutive laws, conditions and discretization study.
No verification case currently reads these files.
