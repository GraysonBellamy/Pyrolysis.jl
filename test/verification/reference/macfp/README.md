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

Codes behind the submissions: Aalto = FDS; BUW-FZJ = FDS (PROPTI, sets A/B);
DBI = FDS + Gpyro; GIDAZE+ = Gpyro; NIST = FDS; Sandia = Sierra
Thermal/Fluids; UCLan = ThermaKin; UMD = ThermaKin2D; UMET = ThermaKin (`_TK`)
and Gpyro (`_GP`). `MaCFP_TGA_10K_*Stats.csv` are the workshop's precomputed
ensemble statistics.

**Caveat for plots:** each institution used its own calibrated property set,
so the spread mixes code differences with calibration differences. For a
comparison that isolates our solver, run our PMMA case with the
community-designated "central" model `MaCFP_PMMA_UMD.json` (committed here;
from `PMMA/Material_Properties/2021/`) under the v1.3 scenario and present
these curves as the community envelope. Intended consumer: the planned B-01
case (`docs/Verification_Suite_Plan_2026-07-06.md`, Tier 4).
