# Gpyro V2 conduction, radiation, kinetics and Stefan comparisons

All Gpyro verification overlays use unmodified V2 revision
`bff0126c68f2a35c084ef04490ae327185f281a3`, built with GNU Fortran 14.4.0
and the upstream `gnu_linux` flags. This directory contains the HC-01,
HC-03, RAD-01, CPL-01 route A and RX-01 comparison datasets.
The executable hash agrees with the `v2_accuracy/` and `matched/` runs.

| Case | Base native unknowns | Intervals / spacing | Requested initial dt |
|:--|--:|:--|--:|
| V-HC-01 | 501 boundary nodes | 500 / 0.1 mm | 0.1 s |
| V-HC-03 | 321 boundary nodes | 320 / 0.0625 mm | 0.1 s |
| V-RAD-01 | 401 boundary nodes | 400 / 0.025 mm | 0.1 s |
| V-CPL-01 route A | 801 boundary nodes | 800 / 0.05 mm | 0.02 s |
| V-RX-01 | one 0D unknown | homogeneous kinetics | 0.001 s |

The four spatial grids match the finest plotted Julia meshes. Boundary
control volumes are half width; unknown counts differ by one from Julia's
cell counts. Full native profiles verify the coordinates and counts. There
is no half-cell shift. Stefan front interpolation uses every native node
in the full profiles.

Each case has `base` (h,dt), `refined` (h/2,dt/2), `time_half` (h,dt/2),
`time_quarter` (h,dt/4), and `space_half` (h/2,dt/4) runs: 25 successful
native solves. The 0D problem keeps one unknown throughout; its last spatial
level repeats the finest temporal configuration and is not spatial evidence.
`DT0` is a requested initial step, not a guarantee of fixed stepping. All 25
logs report successful completion without failed-step messages. The console
step log is sparse and does not replace an accepted-step history.

The generated inputs use half boundary control volumes, fixed geometry,
tight nonlinear solver tolerances and endpoint probes. The gas-species
index is `IO2=1`. RX-01's gas-species solve is disabled because the material
has no porosity and the chain produces no gas.

Native `TSTOP` includes one requested timestep beyond the common observation
horizon solely to retain the final dump. Extracted curves end at the common
horizon. These are accuracy runs without timing measurements; any future
timed comparison must remove or account for that output guard.

`accuracy.csv` records analytical errors at all five levels, evaluated at
native observation coordinates/times. The analytical expressions
in the parent extractor provide the erfc, constant-flux slab, absorption,
Stefan, and Bateman references. `analytical_checks.txt` retains the readable
checks. Native inputs, logs and numerical outputs are under `raw/` and
hashed in `provenance.json`, together with the extracted CSVs and scripts.

Base maximum sampled temperature errors are 0.0552 K (HC-01), 0.0116 K
(HC-03), and 0.0502 K (RAD-01). RX-01's maximum species-fraction error is
9.98e-5, with no printed mass or temperature drift. CPL-01 has a 0.0795 K
maximum error at the selected probe times and a 0.1092 mm maximum front
error at 100, 200 and 299 s. Its finite Gaussian latent-heat peak and front
interpolation differ from the sharp-interface reference; temporal refinement
does not remove that offset. These are sampled errors, not whole-history
error bounds or matched-accuracy performance qualifications.

Reproduce with the same pinned V2 executable:

```sh
python3 test/verification/reference/gpyro/v2_migration/run_migration.py \
  /path/to/gpyro /tmp/gpyro-v2-migration
```

`--extract-only` validates generated inputs, successful termination and grids
before re-extracting an existing workspace.
