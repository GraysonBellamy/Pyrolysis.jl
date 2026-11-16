# LEMTA dry argon case

Decks for the LEMTA tests of seven dry wood species under argon:
- oven-dry blocks, 10 × 10 × 5 cm, grain across the heat flux, in a calcium-silicate
  holder inside a stainless chamber with a sapphire window, argon at 20 L/min;
- a cone heater at 25, 35 and 45 kW/m², 45 mm from the sample;
- the mass of one block per species and flux (0–698 s), and in-depth temperatures of
  another at 35 and 45 kW/m² at 2–20 mm (every 2 mm), 30 and 40 mm.

Every deck value is measured on separate specimens, independently of these cone tests:
TGA kinetics, DSC heats and heat capacities, hot-wire conductivity of wood and char,
FTIR optics, and the chamber's front-face convection and back temperature from
calcium-silicate blocks run in the same chamber. Sources:
- H. Flity, PhD thesis, Université de Lorraine (2023): kinetics, heats, c_p, porosity,
  the chamber;
- Flity et al., Int. J. Therm. Sci. 195 (2024) 108661: wood conductivity and c_p;
- Flity et al., Int. J. Therm. Sci. 210 (2025) 109545: wood optics;
- Flity, Jannot, Acem & Parent, "Thermal conductivity of wood char up to 800 °C"
  (SSRN preprint): char conductivity;
- L. Terrei, PhD thesis, Université de Lorraine (2020): the flux behind the window;
- Jannot & Degiovanni, Int. J. Therm. Sci. 142 (2019) 379: the calcium silicate;
- Thomas & Joseph, Johns Hopkins APL Tech. Dig. 9 (1988) 328: sapphire.

## Data

The data are Flity, Terrei, Acem & Parent, doi:10.57745/QRYAJ9, v1.0, licence Etalab
2.0 (attribution). They stay outside git:
- `LEMTA_ARGON_DATA` points at a directory holding the published files; the default is
  `~/Documents/data/lemta-QRYAJ9`.
- `fetch_data.sh` downloads each spreadsheet as uploaded, checked against the MD5 the
  repository states, and the repository's tab-delimited export of it, checked against
  the SHA-256 pinned in `data_files.toml`. The export holds the spreadsheet's cell
  values unchanged; the loader reads it.
- The comparison also scores the char depths of the same chamber's cone tests: Terrei,
  Erez, Acem, Suzanne, Boulet & Parent, doi:10.57745/I9CL2R, v1.1, licence Etalab 2.0,
  in `LEMTA_CHAR_DATA` (default `~/Documents/data/lemta-I9CL2R`), pinned in
  `char_data_files.toml` and fetched by the same script.
- `check_data()` checks every file of both against the SHA-256s of the two manifests.

```sh
bash scripts/lemta_argon/fetch_data.sh
julia --project=scripts/verification scripts/lemta_argon/check_data.jl
julia --project=scripts/verification scripts/lemta_argon/calibrate_front.jl
julia --project=scripts/verification scripts/lemta_argon/smoke.jl poplar
julia --project=scripts/verification scripts/lemta_argon/run_cases.jl deck
julia --project=scripts/verification scripts/lemta_argon/compare.jl
```

`run_cases.jl` takes one or more variants of `variants.jl` (the decks; 200 cells and
tighter tolerances; one lever each; a diagnostic with the thesis's fitted
conductivity, which is fitted to these tests and is not a deck value), optionally test
or species names and `--part=i/n` to split the runs over processes.

## Layout

| path | what |
|---|---|
| `decks/chamber.toml` | the rig and the laws every species shares: c_p of wood and char, k of wood and char, char optics, the heat-only step, gas transport, the thesis's front convection, the calcium silicate, the cone and window |
| `decks/silcal.csv` | the calcium-silicate runs in the chamber (thesis Fig. 5.5, read from its vector paths): temperatures at 2–40 mm and the rear face, at 25, 35 and 45 kW/m² |
| `decks/<species>.toml` | kinetics, heats, optics and char law, each with its source |
| `LemtaArgon.jl` | data loading and the defect list; decks to Pyrolysis.jl materials and problems; outputs at the thermocouple depths with the conservation ledgers |
| `data_files.toml`, `char_data_files.toml`, `fetch_data.sh` | the pinned datasets and their download |
| `chamber_radiation.jl` | the cone, the sapphire window and the face's radiation exchange |
| `decks/wood_absorptivity.csv` | the woods' spectral absorptivity (Flity et al. 2025 Fig. B.18, read from its vector paths) |
| `calibrate_front.jl`, `front.toml` | the front-face h per flux, fitted to the calcium-silicate runs in 1D, under the cone and window and with a gray face |
| `check_data.jl` | every test with its cuts, drawn per species (`results/figures/data_<species>.png`); `results/data_summary.csv` |
| `smoke.jl` | one species at each flux: solves, ledgers, a sanity check against the data, and sensitivities on its 35 kW/m² tests |
| `variants.jl`, `run_cases.jl` | the comparison's variants, and their runs on every test (`results/<variant>/<test>.csv`; mass tests to 720 s, the char-depth tests' exposure) |
| `compare.jl` | scores against the data: m/m₀, mass loss rates, temperatures by depth, char depth and thickness lost (`results/metrics_*.csv`, `results/tables/*.md`, `results/figures/compare_<species>.png`, `levers.png`) |
| `results/` | generated (ignored) |

## Scenario

- **Geometry and mesh.** A 1D block of 50 mm across the grain, 100 equal cells. The
  chars are denser than their reactants so that a fully converted layer keeps
  1 − `shrinkage` of its thickness (default 0.20, charred blocks across the grain;
  `:cone`, the species' share in this chamber's cone tests; `nothing`, none on a fixed
  mesh). The mesh moves with the solid, and each thermocouple is read at its material
  point.
- **Material.** Three pseudo-components of the dry wood decompose independently at
  first order: S1 → volatiles, S2 → γ2 C2 + volatiles, S3 → γ3 C3 + volatiles, and the
  heat-only C2 → C4 keeps C2's mass. Each pseudo-component has the test's density, so
  the three fill the volume in their shares α; each char holds its reactant's volume.
  - Wood: c = 1131 + 4.67 T(°C); k = 0.0386 + 0.1708 (ρ/1000) √(T/433 K) across the
    grain; ε(T) and α of the species.
  - Char: c = 693 + 3.44 T(°C); k = aT + b + 16σT³/(3β) of the species' char, or the
    law pooled over three chars; ε = α = 0.95.
  - Conductivity mixes in parallel over the condensed volume shares; emissivity and
    absorptivity by those shares.
  - Heats are per g of initial dry wood in the decks, H1/α1, H2/α2 and H4/(α2γ2) per kg
    of reactant in the model, positive endothermic; H3 = 0.
- **Gas.** One lumped volatile species and argon in the pores; every solid at 1e-12 m²
  so the volatiles leave through the exposed face at near-ambient pressure (the block
  vents along its vertical grain, which 1D cannot); the back is impermeable.
- **Exposed face** (`front`, default `:window`). Under the cone and the sapphire window
  (`chamber_radiation.jl`): the coil, gray on the ISO 5660 frustum, at the
  temperature that sends the nominal flux through the cold window to the sample's
  centre; the window's spectral transmittance, reflectance and absorptance from its
  phonon parameters and multiphonon edge (Thomas & Joseph, Johns Hopkins APL Tech.
  Dig. 9 (1988) 328); its temperature from its own balance; the face absorbs the
  transmitted coil radiation at its absorptivity for that spectrum (the woods' from
  their measured spectra, `decks/wood_absorptivity.csv`; char 0.95), and the window's
  emission and its own reflected emission at its emissivity; window and sample are
  solved in turn. Convection to T_a(t) at the h fitted on the calcium-silicate runs
  under the same window (`calibrate_front.jl`, `front.toml`). `front = :silcal`
  instead absorbs the nominal flux at the deck's gray absorptivity and re-radiates
  to T_a(t), with the gray h of the same calibration; `front = :thesis` takes the
  thesis's h.
- **T_a(t).** The test's initial 40 mm temperature plus the rise of the
  calcium-silicate rear face at the same flux (`decks/silcal.csv`), continued past
  600 s at its last slope.
- **Back face** (`back`): `:chamber` holds it at T_a(t); `:tc40` ends the domain at
  40 mm and holds it at the test's own 40 mm record; `:adiabatic`.
- **Initial state.** A temperature test starts from its first row, interpolated in
  depth; a mass test starts uniform at the chamber's initial gas temperature.
- **Char depth.** The char/virgin boundary is where the conversion of the three
  pseudo-components in a material cell, X = 1 − (S1 + S2 + S3)/(their initial mass),
  crosses 0.5 (0.2 and 0.8 bracket it; `CHAR_THRESHOLDS`), measured from the current
  face; the thickness lost is the initial thickness less the current one.

## Defects

- `Mass_Balsa_35kW.m-2_143kg.m-3` and `Temperature_Ipe_35kW.m-2_842kg.m-3` are the
  45 kW/m² tests: the loader reads them at 45.
- 12 thermocouples give no reading; crosstalk windows on balsa 45, beech 35, ipé 842 and
  larch 45 are cut; beech 35's exposure ends at 1078 s.
- Balsa, oak and beech: the glue line opened and the heater reached the 2–6 mm
  junctions, which then read above the wood; those records are cut above 400 °C.
- The balance is disturbed for the first 4–29 s of several mass tests (`t_settle`).
