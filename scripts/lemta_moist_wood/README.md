# LEMTA moist-wood case

Runs the LEMTA cone-heater drying tests in Pyrolysis.jl and compares them with
the measurements:
- wood blocks at 0–30 % moisture (dry basis), 10 × 10 × 5 cm, grain across the
  heat flux;
- 4 kW/m² for 2.5 h, then cooling;
- thermocouples at 2–20 mm (every 2 mm), 30 and 40 mm, and mass, at 1 Hz.

The comparison scores the moisture model, the phase change of the bound water
and its vapour toward sorption equilibrium, and the one-way evaporation laws
beside it, against the in-depth temperatures, the mass, the mass regained or
lost after the cone is withdrawn, and the sliced moisture profiles.

Sources:
- Diarra et al., Fire Saf. J. 165 (2026) 104951 (the tests; Fig. 11, the
  moisture profiles);
- Flity et al., Int. J. Therm. Sci. 210 (2025) 109545 (boundary conditions and
  optics on the same rig);
- Flity et al., Int. J. Therm. Sci. 195 (2024) 108661 (conductivity and heat
  capacity);
- Dietenberger, Glass & Boardman 2024 (the CAST isotherms).

## Data

The data are Terrei, Diarra & Parent, doi:10.57745/CBLZ8P, v1.0, licence Etalab
2.0 (attribution). They stay outside git:
- `LEMTA_DATA` points at a directory holding the published files; the default is
  `~/Documents/data/lemta-CBLZ8P`.
- `fetch_data.sh` downloads them through the repository's API and checks the MD5s
  the repository states.
- `check_data()` checks each file against the SHA-256 pinned in `data_files.toml`.

`data/diarra_2026_fig11.csv` holds the 22 points of Diarra et al.'s Fig. 11(a)
(spruce at about 10 %, 5 mm slices after 30, 60 and 120 min), read from the
figure to ±0.15 point of moisture content and ±0.5 mm.

```sh
bash scripts/lemta_moist_wood/fetch_data.sh
julia --project=scripts/verification scripts/lemta_moist_wood/check_data.jl
julia --project=scripts/verification scripts/lemta_moist_wood/calibrate_dry.jl
julia --project=scripts/verification scripts/lemta_moist_wood/run_cases.jl dry dry_cal sorption_cal gated_cal chan_cal arrhenius_cal
julia --project=scripts/verification scripts/lemta_moist_wood/compare.jl
julia --project=scripts/verification scripts/lemta_moist_wood/profiles.jl
```

The transfer coefficient's ladder and the other sensitivities are variants like
any other, e.g.

```sh
julia --project=scripts/verification scripts/lemta_moist_wood/run_cases.jl sorption_cal_lambda5e-7 sorption_cal_lambda5e-6
julia --project=scripts/verification scripts/lemta_moist_wood/run_cases.jl sorption_cal_rh0 sorption_cal_rh65 sorption_cal_freewater
```

A run of one variant takes a few seconds to half a minute (the 22 moist runs two to
four minutes); the variants can run side by side as separate processes.

## Layout

| path | what |
|---|---|
| `decks/<species>.toml` | per species: dry density, absorptivity, ε(T), convection coefficients, permeability, isotherm, each with its source |
| `LemtaMoistWood.jl` | data loading and the defect list; decks to Pyrolysis.jl problems; the moisture model and the one-way laws; outputs at the thermocouple depths; metrics; the variants run |
| `data_files.toml`, `fetch_data.sh` | the pinned dataset and its download |
| `data/diarra_2026_fig11.csv` | the measured moisture profiles |
| `check_data.jl` | every run with its defects applied, plotted per species; `results/data_summary.csv` |
| `calibrate_dry.jl`, `calibration.toml` | the factors on the absorbed flux and on k that each species' dry runs need |
| `run_cases.jl` | runs variants: `results/<variant>/<run>.csv` and, for a moist run, `<run>_profiles.csv` |
| `compare.jl` | metrics (`results/metrics_*.csv`), tables (`results/tables/*.md`, the phase change's printed), figures (`results/figures/`) |
| `profiles.jl` | the moisture profiles of `spruce_11,42%` against Fig. 11(a); `results/tables/profiles_diarra.md` |
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
- **Mass exchange:** both faces exchange gas with the ambient air through the analogy
  films of their convective coefficients, in both directions. The air's vapour mole
  fraction is `ambient_rh × p_sat(T_a)/P` (`WATER_SATURATION`).
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
- **Gases:** water vapour and pore air. Darcy permeability from the deck, lumped
  transfer coefficient `λ`, or the dusty-gas closure with `tortuosity`.
- **Dry calibration.** `calibrated = true` applies the species' factors from
  `calibration.toml`. Fir has no dry run.

## The moisture model

The default law, `:sorption`, is a `PhaseChange` of the water `:MOISTURE` and its
vapour `:WATER_VAPOR` with the wood as the sorbing solid. In every cell the water and
the pore vapour exchange toward sorption equilibrium,
p_v = a(X, T) p_sat(T), at the rate K (ξ_v^eq − ξ_v): evaporation where the vapour
is below it, condensation and re-sorption where it is above. The step takes the
latent heat built from the water's and the vapour's heat capacities plus the
isotherm's heat of sorption; nothing about the heat is set in the suite. The pore
vapour starts in equilibrium with the water and the air makes up the ambient
pressure. Its inputs, each a keyword of `problem` and `simulate`:

| input | default | why |
|---|---|---|
| isotherm | the deck's `isotherm`: `CAST_SPRUCE` for spruce, fir and larch, `CAST_BEECH` for oak and poplar | the two CAST sets, softwood and hardwood |
| `ambient_rh` | 0.5 at the ambient temperature | the laboratory humidity is not recorded; the mass regained or lost during cooling is scored at 50 % |
| `λ` | 1.5 × 10⁻⁶ m²/s | the dry wood's lumped transfer coefficient stands for bound-water diffusion across the grain; measured moisture diffusion coefficients put it near 5 × 10⁻⁷ at 10–12 % and two to four times that at 15–20 % (Siau 1984; Keey et al. 2000). The samples span 7–26 %. Here the heating-rate ratio at 100/70 °C matches the measured 0.30; at 5 × 10⁻⁶ the rms is 0.2 K lower and the mass closer, but the ratio is twice the measured |
| `exchange_rate` | 10³ s⁻¹ | the inverse relaxation time of the pore vapour; between 10² and 10⁴ s⁻¹ the scores move by under 0.01 K and 0.001 in mass |
| `tortuosity` | none (lumped transfer) | a value selects the dusty-gas closure, the vapour diffusing through the pore air at 2.6 × 10⁻⁵ m²/s (25 °C, 1 atm) reduced by φ/tortuosity |

`:freewater` is the same step with the activity of free water (`FreeLiquid`, no
isotherm): it cannot hold water below saturation and shows what the isotherm does.

The phase change's runs add to the series the gross condensation, the largest
p_v/p_eq where the wood holds water and the largest moisture content, and a line
with the accepted steps and the temperatures at which the water evaporates. Every
moist run also writes its moisture profiles (`<run>_profiles.csv`: X, p_v, p_eq and
P by cell at 0, 0.5, 1, 1.5, 2, 2.5 h and the end).

## The one-way laws

The entries of `LAWS` make the water evaporate by a `Reaction`, one way, whatever
the vapour pressure: they cannot condense or regain water. They run at λ = 5 × 10⁻⁷
m²/s in dry air (`ONE_WAY`), the settings at which they are scored; a one-way step's
temperatures and mass hardly depend on either. What each is for:
- `:arrhenius`, examples/07's TGA set: a rate law that evaporates at room
  temperature over hours, ruled out by these data;
- `:gated`, the same switched on at 373.15 K over 2 K: the best one-way law, the
  reference the phase change is judged against;
- `:chan`, A = 5.13e10 s⁻¹, E = 88 kJ/mol, the drying set quoted in the literature;
- the `:onset*` laws, a scan of where a one-way sink would have to sit.

Their heat is the latent heat of free water at the cell temperature plus the heat
of sorption of `examples/07` (`heat_of_desorption`).

## Variants

`VARIANTS` in `LemtaMoistWood.jl`; `results/<variant>/`:
- `dry`, `dry_cal` and the dry runs' diagnostics (`dry_cells200`, `dry_tight`,
  `dry_handoff_bc`, `dry_probe_flity`, `dry_k_table1`, `dry_flux3.8`);
- `sorption` (the published deck) and `sorption_cal` (the calibrated deck, the
  primary); on the calibrated deck `sorption_cal_<tag>` with the tags
  `lambda5e-7`, `lambda5e-6` (the transfer coefficient's ladder), `rh0`, `rh65` (the
  ambient humidity), `vent` (κ = 10⁻¹⁵ m², side venting along the grain), `K1e2`,
  `K1e4` (the exchange rate), `cells200`, `tight` (reltol 10⁻⁸, dtmax 10 s),
  `dusty4`, `dusty12` (the dusty-gas closure at tortuosity 4 and 12) and
  `freewater`;
- for each one-way law, `<law>`, `<law>_cal` and on the calibrated deck
  `<law>_cal_<tag>` with `kappa0.1`, `kappa10`, `vent`, `sorption0`, `sorption2` (the
  heat of sorption ×0 and ×2), `kseries`, `cells200`, `tight`; `onset*_cal`.

## Defects

- `oak_test2_0%` is excluded: its channels are scrambled.
- `fir_17,67%` and `fir_9,1%`: only their mass loss is read.
- `spruce_test2_0%` starts with a 20 → 40 °C gradient and cools through its sides. It
  is reported, but left out of the calibration.
- The balance has single-sample spikes: measured mass losses are ±30 s medians.
