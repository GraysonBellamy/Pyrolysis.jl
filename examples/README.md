# Pyrolysis.jl examples

The catalogue of these scripts — what each one models and how to run it — is
the documentation's [Examples page](../docs/src/examples.md).

Run any script from the repository root in the examples environment:

```sh
julia --project=examples examples/06_cone_calorimeter.jl
```

Scripts that write files write them into gitignored locations: the repository
root (`cone_calorimeter_*.png`, `douglas_fir_*`, `ptw_*`), `output/`, and the
`results/` and `*_results/` directories beside the scripts.

## Smoke runs

Every simulation script and the optimizer script take the environment variable
`PYROLYSIS_SMOKE=1`: a short run that checks the script works. It shortens the
time horizon to 10–20 s with a handful of output times, and where a script
loops over many cases it runs a subset. Without the variable a script runs in
full. A smoke run writes the same files as the full run, so it replaces the
outputs of an earlier full run of the same script (the figures of the two
`plot_*` scripts excepted, which take a `_smoke` suffix).

```sh
PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter.jl
```

## Simulation scripts

- **`06_cone_calorimeter.jl`** — PMMA, 5.8 mm, in the cone calorimeter at
  25 kW/m² on a moving (ALE) mesh with cell deactivation; plots MLR,
  temperatures, mass and thickness. Fixtures: none. It overlays ThermaKin
  outputs from `docs/thermakin_output/` when that directory exists (it is not in
  the repository, so the overlay is skipped). Smoke: 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter.jl`
- **`06_cone_calorimeter_surface_absorption.jl`** — the same PMMA case with
  `radiation_model = SURFACE_ABSORPTION` set explicitly (the radiation is
  absorbed at the exposed face, as under 06's default);
  writes figures and a CSV to `output/`. Fixtures: none; the ThermaKin overlay
  reads `docs/thermakin_surface_absorption_output/` when present (not in the
  repository). Smoke: 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter_surface_absorption.jl`
- **`07_douglas_fir_cone.jl`** — Douglas fir, 12 mm, at 75 kW/m²: five-step
  pyrolysis with moisture evaporation, lumped wall and Darcy gas transport with background
  pore air, surface absorption, ALE mesh; writes figures and the
  `douglas_fir_results_<q>kWm2_fit.csv` / `douglas_fir_spatial_<q>kWm2_fit.csv`
  tables to the working directory. Fixtures: none. Smoke: 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/07_douglas_fir_cone.jl`
- **`08_gas_transport.jl`** — a charring wood-like slab at 75 kW/m² with
  pressure-driven (Darcy) and diffusive gas transport; prints MLR, surface
  temperature and the mass budget. Fixtures: none. Smoke: 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/08_gas_transport.jl`
- **`08_pressure_treated_wood_cone.jl`** — pressure-treated wood, 12 mm, at
  50 kW/m² with lumped wall and Darcy gas transport; writes figures and
  `ptw_results_50kWm2.csv` / `ptw_spatial_50kWm2.csv` to the working directory.
  Fixtures: none. Smoke: 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/08_pressure_treated_wood_cone.jl`
- **`convergence_douglas_fir.jl`** — grid convergence of the Douglas fir cone
  case at 480, 240, 120, 60, 30 and 15 cells (finest is the reference); writes
  tables and figures to `examples/convergence_results/`. Runs in parallel under
  `julia -p N`. Fixtures: none. Smoke: the two coarsest grids over 10 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/convergence_douglas_fir.jl`
- **`moisture_hos_comparison.jl`** — Douglas fir at four initial moisture
  contents (3, 6, 9, 12 % dry basis) under an exponential and a constant
  heat-of-sorption model; writes CSVs and figures to
  `examples/moisture_hos_comparison_results/`. Fixtures: none. Smoke: every case
  over 10 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/moisture_hos_comparison.jl`
- **`literature_comparison/literature_comparison.jl`** — the Douglas fir cone
  case with literature heat capacities, heats of reaction and heat of sorption
  substituted one category at a time (baseline plus eight variants; method and
  sources in `literature_comparison_plan.md`); writes CSVs to
  `examples/literature_comparison/results/`. Runs in parallel under
  `julia -p N`. Fixtures: none. Smoke: the baseline and the first variant over
  10 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/literature_comparison/literature_comparison.jl`
- **`macfp/macfp_gasification.jl`** — the MaCFP-4 idealized white-pine
  gasification targets: 10, 30 and 60 kW/m² on 6 and 18 mm samples, with the
  calibrated white-pine model of `inverse_analysis/wp_forward.jl`; writes CSVs
  and figures to `examples/macfp/results/macfp_<timestamp>/`. Runs the cases on
  threads under `julia -t N`. Fixtures: none. Smoke: every case over 10 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/macfp/macfp_gasification.jl`
- **`sensitivity_analysis/sensitivity_douglas_fir.jl`** — propagates the
  thermophysical-property uncertainties (16 parameters, 8192 Sobol samples)
  through the Douglas fir cone case; resumable, writes per-sample curves,
  scalar QoIs, SRRC rankings and prediction envelopes to
  `examples/sensitivity_analysis/results/` (findings in
  `sensitivity_analysis_report.md`). Runs in parallel under `julia -p N`.
  Fixtures: none (an existing `results/` directory is reloaded and completed). Smoke: one forward solve at the nominal parameters over 10 s,
  no files written.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/sensitivity_analysis/sensitivity_douglas_fir.jl`
- **`inverse_analysis/validate_thickness.jl`** — the calibrated white-pine
  model with variable cross section at 30, 50 and 70 kW/m²: final thickness,
  area and mass ratios against the CAPA measurements. Fixtures: the UMD CAPA
  and MaCFP data in `examples/inverse_analysis/data/` (committed). Smoke: every
  scenario over 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/inverse_analysis/validate_thickness.jl`
- **`inverse_analysis/plot_stage2_theta.jl`** — forward-solves the fit and
  validation scenarios at a given Stage-2 conductivity vector and plots them
  against the experiments (`results/fit_plots_stage2_manual.png`). Fixtures:
  `examples/inverse_analysis/data/` (committed). Smoke: every scenario over 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/inverse_analysis/plot_stage2_theta.jl`
- **`inverse_analysis/plot_diagnostics.jl`** — forward-solves each flux at a
  given Stage-2 conductivity vector and plots back and front temperature, mass,
  MLR and thickness against the experiments
  (`results/diagnostics_stage2_manual.png`). Fixtures:
  `examples/inverse_analysis/data/` (committed). Smoke: every scenario over 20 s.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/inverse_analysis/plot_diagnostics.jl`

## Optimizer script

- **`inverse_analysis/inverse_white_pine.jl`** — recovers white-pine component
  conductivities by fitting CAPA back-surface temperatures at several cone heat
  fluxes, with an optional HFM/LFA conductivity prior and Hessian or bootstrap
  covariance (`--stage`, `--fit`, `--val`, `--prior-weight`, `--bootstrap`; see
  the script's docstring); writes `trace.csv`, `best.json` and figures to
  `examples/inverse_analysis/results/run_<timestamp>_s<stage>_pw<weight>/`.
  Runs the scenarios on threads under `julia -t N`. Fixtures: the FSRI, UMD and
  MaCFP CAPA and gasification series, `hfm_dry_summary.csv` and
  `lfa_summary.csv` in `examples/inverse_analysis/data/` (committed). Smoke: one
  objective evaluation at the initial θ, the fit scenarios over 20 s, then exit.
  `PYROLYSIS_SMOKE=1 julia --project=examples examples/inverse_analysis/inverse_white_pine.jl`

The white-pine model lives in three included files: `wp_properties.jl`
(properties, kinetics and constants), `wp_forward.jl` (the forward model) and
`wp_objective.jl` (data loading and the loss). The Python helpers
`derive_cone_view_factor.py`, `derive_shrinkage.py` and
`derive_surface_emissivity.py` derive constants of `wp_properties.jl` from the
data in `examples/inverse_analysis/data/`.

## Analysis scripts

These read result files and have no smoke mode.

- **`df_sensitivity_comparison.jl`** — Douglas fir at 50 kW/m²: Pyrolysis.jl
  against ThermaKin (MLR, temperatures, mass, per-component MLR); writes
  figures to `output/`. Fixtures: `docs/DF_sensitivity_Cone_output.txt`
  (committed) and `douglas_fir_results_50kWm2_fit.csv` in the repository root,
  generated by `07_douglas_fir_cone.jl` run with `q_external = 50e3`.
- **`ptw_comparison.jl`** — pressure-treated wood at 50 kW/m²: Pyrolysis.jl
  against ThermaKin; writes figures to `output/`. Fixtures:
  `docs/PTW_Cone_output.txt` (committed) and `ptw_results_50kWm2.csv` in the
  repository root, generated by `08_pressure_treated_wood_cone.jl`.
- **`sensitivity_analysis/sensitivity_sobol.jl`** — group Sobol indices of the
  Douglas fir sensitivity study through a polynomial chaos expansion, with
  bootstrap confidence intervals and convergence checks. Fixtures:
  `results/parameters.csv` and `results/scalar_qois.csv`, generated by
  `sensitivity_douglas_fir.jl`.
- **`sensitivity_analysis/compare_degrees.jl`** — compares PCE degrees 2, 3 and
  4 for the Sobol analysis (method in `pce_degree_selection.md`). Fixtures: the
  same `results/` files as `sensitivity_sobol.jl`.

## Excluded scripts

- `sensitivity_analysis/sensitivity_sobol.jl` and
  `sensitivity_analysis/compare_degrees.jl` do not load in the examples
  environment: they need `PolyChaos`, which `examples/Project.toml` does not
  list. They also need the 8192-sample result set of a full
  `sensitivity_douglas_fir.jl` run.
- `df_sensitivity_comparison.jl` needs `douglas_fir_results_50kWm2_fit.csv`,
  which is not in the repository; `07_douglas_fir_cone.jl` writes it when run
  at 50 kW/m².
- `ptw_comparison.jl` needs `ptw_results_50kWm2.csv`, which is not in the
  repository; a full run of `08_pressure_treated_wood_cone.jl` writes it.
