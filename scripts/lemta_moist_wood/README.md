# LEMTA moist-wood case

Runs the LEMTA cone-heater drying tests in Pyrolysis.jl and compares them with
the measurements:
- wood blocks at 0–30 % moisture (dry basis), 10 × 10 × 5 cm, grain across the
  heat flux;
- 4 kW/m² for 2.5 h, then cooling;
- thermocouples at 2–20 mm (every 2 mm), 30 and 40 mm, and mass, at 1 Hz.

The comparison scores the evaporation laws against the in-depth temperatures and
the mass.

Sources:
- Diarra et al., Fire Saf. J. 165 (2026) 104951 (the tests);
- Flity et al., Int. J. Therm. Sci. 210 (2025) 109545 (boundary conditions and
  optics on the same rig);
- Flity et al., Int. J. Therm. Sci. 195 (2024) 108661 (conductivity and heat
  capacity).

## Data

The data are Terrei, Diarra & Parent, doi:10.57745/CBLZ8P, v1.0, licence Etalab
2.0 (attribution). They stay outside git:
- `LEMTA_DATA` points at a directory holding the published files; the default is
  `~/Documents/data/lemta-CBLZ8P`.
- `fetch_data.sh` downloads them through the repository's API and checks the MD5s
  the repository states.
- `check_data()` checks each file against the SHA-256 pinned in `data_files.toml`.

```sh
bash scripts/lemta_moist_wood/fetch_data.sh
julia --project=scripts/verification scripts/lemta_moist_wood/check_data.jl
julia --project=scripts/verification scripts/lemta_moist_wood/calibrate_dry.jl
julia --project=scripts/verification scripts/lemta_moist_wood/run_cases.jl dry dry_cal gated_cal chan_cal arrhenius_cal
julia --project=scripts/verification scripts/lemta_moist_wood/compare.jl
```

## Layout

| path | what |
|---|---|
| `decks/<species>.toml` | per species: dry density, absorptivity, ε(T), convection coefficients, permeability, each with its source |
| `LemtaMoistWood.jl` | data loading and the defect list; decks to Pyrolysis.jl problems; the evaporation laws; outputs at the thermocouple depths; metrics; the variants run |
| `data_files.toml`, `fetch_data.sh` | the pinned dataset and its download |
| `check_data.jl` | every run with its defects applied, plotted per species; `results/data_summary.csv` |
| `calibrate_dry.jl`, `calibration.toml` | the factors on the absorbed flux and on k that each species' dry runs need |
| `run_cases.jl` | runs variants: `results/<variant>/<run>.csv` |
| `compare.jl` | metrics (`results/metrics_*.csv`), tables (`results/tables/*.md`), figures (`results/figures/`) |
| `results/` | generated (ignored) |

## Scenario

- **Geometry and mesh.** A 1D slab of 50 mm across the grain, 100 equal cells, on a
  fixed mesh.
- **Initial and ambient state.** The initial temperature is each run's first row,
  interpolated in depth. The ambient is the shallowest thermocouple's first reading.
- **Exposed face:**
  - absorbs α × 4.0 kW/m² while the cone is on (α from the deck);
  - re-radiates at ε(T) to the ambient;
  - convects at the deck's heating coefficient, then its cooling coefficient. The
    switch is a 3 s smooth step at the withdrawal time read from the data.
- **Back face:** the deck's `h_rear` is the whole loss there (Flity et al. 2025 fit no
  radiation on it).
- **Mass exchange:** both faces exchange gas with dry ambient air through the analogy
  films of their convective coefficients.
- **Dry wood:**
  - its density is the filename's;
  - k across the grain is Flity et al. 2024 Eq. 9, k = 0.0386 + 0.1708 (ρ/1000)
    √(T/433);
  - c = 1131 + 4.67 (T − 273);
  - skeletal density 1500 kg/m³.
- **Water.**
  - Bound water is a liquid component at 1000 kg/m³ and 0.6 W/(m K), mixed in
    parallel. That gives the conductivity rise with moisture that Flity et al. 2024
    measured.
  - Initial water is MC/(1 + MC) × m₀, or MC ρ_dry V for the two fir files whose first
    mass includes the holder.
  - Evaporation takes the latent heat at the cell temperature plus the heat of
    sorption of `examples/07`.
- **Gases:** water vapour and pore air. Darcy permeability from the deck, lumped
  transfer λ = 5e-7 m²/s.
- **Evaporation laws (`LAWS`):**
  - `:arrhenius`, examples/07's set;
  - `:gated`, the same switched on at 373.15 K over 2 K;
  - `:chan`, A = 5.13e10 s⁻¹, E = 88 kJ/mol;
  - the `:onset*` laws, a scan of the switch temperature.
- **Dry calibration.** `calibrated = true` applies the species' factors from
  `calibration.toml`. Fir has no dry run.

## Defects

- `oak_test2_0%` is excluded: its channels are scrambled.
- `fir_17,67%` and `fir_9,1%`: only their mass loss is read.
- `spruce_test2_0%` starts with a 20 → 40 °C gradient and cools through its sides. It
  is reported, but left out of the calibration.
- The balance has single-sample spikes: measured mass losses are ±30 s medians.
