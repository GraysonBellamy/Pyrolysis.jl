# Third-party reference data for verification report overlays

Numerical results from other pyrolysis/heat-transfer solvers, used only by
report mode (`scripts/verification/`) to overlay competing codes' curves and
errors on our verification figures. Assert mode (`Pkg.test()`) never reads
this directory.

Each subdirectory records the available version, input, run-date and extraction
information. FDS files are native outputs; Gpyro reference CSVs are extracted
and converted to SI; MaCFP files are row-downsampled. Plotting the retained
curves is reproducible without the external codes. Reproducing every original
run is a stronger claim: the V2 CSV lacks a pinned revision. Repeated FDS
and Gpyro 0.8200 timing measurements have separate provenance records.

| Directory | Code | What it feeds |
|---|---|---|
| `fds/` | FDS 6.11.1 (NIST), run locally on the shipped verification inputs | V-HC-01/02, V-RX-01/02, V-CPL-03 overlays |
| `gpyro/` | Gpyro 0.8200 (Lautenberger), built from source, dissertation §3.4 decks | V-HC-01/03, V-RAD-01, V-CPL-01, V-RX-01 overlays |
| `macfp/` | MaCFP-2 PMMA community submissions (9 institutions: FDS/Gpyro/ThermaKin/Sierra) | B-01 community-envelope comparison (case pending) |
| `timings.csv` | internal elapsed and full process wall times; timing protocols recorded with the data | per-case `solve_time` figure and supporting timing table |

## Timing convention

The primary figures compare **setup and simulation wall time**. Julia times
`case.build(n)` and the solve after a discarded warm-up. Its process startup,
compilation warm-up, and report/plot generation are outside the measurement.
FDS and Gpyro use their **internal total elapsed wall time**, not their
time-stepping-only or per-routine timing statistics. The `wall_s` column in
`timings.csv` holds this internal elapsed time; `process_wall_s` separately
holds full process wall time. Both reference columns are medians of five
fresh-process runs after one warm-up, on one pinned CPU.

The published Julia timing figures and wall-time columns retain measurements
from source revision `a3016a55a68f7e1a2a5baf451f7cfd7dbdf9cb74`. They do not
measure the current gas transport implementation. Updated numerical errors and
QoIs come from separate verification runs. Configurations without a retained
timing have an empty `wall_s` field in CSV and a dash in the displayed table.

| Code | Primary timer starts | Primary timer ends |
|---|---|---|
| Pyrolysis.jl | before model construction in a warmed Julia session | after solve returns |
| FDS 6.11.1 | after MPI initialization; before input reading | final timing report before MPI finalization |
| Gpyro 0.8200 | before input reading and allocation | shutdown timing print after simulation and final timing-file output |

Full process wall times include launch and shutdown and appear in supporting
tables, not in ratios against warmed Julia. No complete-process Julia
measurement is supplied here. FDS/Gpyro output remains inside their internal
timers; Julia report output is outside its timer. Discretizations, adaptive
steps and output workloads also differ, so these are indicative costs, not
pure solver speedups. The matched timing decks use the same simulated spans
as the corresponding Julia cases; per-row details are in the note column.

See [FDS provenance](fds/README.md) and [Gpyro provenance](gpyro/README.md) for
build details, timer boundaries, repeated measurements and input hashes.

ThermaKin2Ds has no retained verification output in this directory. The local
executable in `docs/refs/ThermaKin/` is a Windows PE32+ console binary; this
observation does not characterize every ThermaKin distribution. Stoliarov and
Lyon's [original 2008 paper](https://publications.iafss.org/publications/fss/9/1141/view/fss_9-1141.pdf),
pp. 1146–1149, verifies semi-infinite conduction, absorption, four coupled
reactions and diffusion from a thin layer. These are analytical precedents,
not evidence that a ThermaKin2Ds run has been performed for this report.
