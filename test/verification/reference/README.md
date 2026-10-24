# Third-party reference data for verification report overlays

Numerical results from other pyrolysis/heat-transfer solvers, used only by
report mode (`scripts/verification/`) to overlay competing codes' curves and
errors on our verification figures. Assert mode (`Pkg.test()`) never reads
this directory.

Each subdirectory has its own README recording exact provenance (code
version, inputs, executable identifiers, and how to regenerate). Data files are committed
verbatim from the runs (MaCFP files are row-downsampled; see `macfp/README.md`)
so plots are reproducible without installing the other codes.

| Directory | Code | What it feeds |
|---|---|---|
| `fds/` | FDS 6.11.1 (NIST), run locally on the shipped verification inputs | V-HC-01/02, V-RX-01/02, V-CPL-03 overlays |
| `gpyro/` | Gpyro 0.8200 (Lautenberger), built from source, dissertation §3.4 decks | V-HC-01/03, V-RAD-01, V-CPL-01, V-RX-01 overlays |
| `macfp/` | MaCFP-2 PMMA community submissions (9 institutions: FDS/Gpyro/ThermaKin/Sierra) | B-01 community-envelope comparison (case pending) |
| `timings.csv` | wall times of the FDS/Gpyro runs above; timing protocols recorded with the data | per-case `solve_time` comparison figure |

Timing caveat: the reference runs use each code's own resolution, time step,
and (for the shipped FDS heat_conduction cases) a longer simulated span — the
solve-time figures are indicative, not a controlled benchmark; per-row
details are in `timings.csv`'s note column.

ThermaKin2Ds is not yet represented: the distribution in `docs/refs/ThermaKin/`
is a Windows-only console binary (PE32+), and no wine/mono is assumed on dev
machines. Its verification cases (Stoliarov & Lyon, Fire Safety Science 9,
2008, pp. 1146–1149 — semi-infinite conduction, Beer–Lambert absorption,
4-reaction system, Gaussian plane source ≙ V-GT-01) are fully parameterized in
the paper; run the .exe on any Windows box (or wine) and drop `t`-vs-QoI CSVs
into `thermakin/` following the same provenance pattern.
