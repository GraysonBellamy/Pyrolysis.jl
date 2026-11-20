# Bilbao et al. 1996 wet-wood case

*Pinus pinaster* slabs, 110 × 110 × 19 mm, ρ_dry = 430 kg/m³, at 0, 9 and 20 % moisture
(dry basis), heated from one face by a radiant panel 35 mm above them in a furnace purged
with N₂ (O₂ < 1 %); the other faces in fibreglass. The flux rises with the panel from cold
(five settings A–E, measured with a Gardon gauge); temperatures at 1, 9 and 19 mm (the back
face) and the mass every 60 s. Source: R. Bilbao, J.F. Mastral, J. Ceamanos, M.E. Aldea,
*J. Anal. Appl. Pyrolysis* 36 (1996) 81–97.

The case scores the moisture model of `scripts/lemta_moist_wood`, the phase change of the
bound water and its vapour toward sorption equilibrium, and the one-way evaporation laws
beside it (the `examples/07` Arrhenius step, the same step gated at 373 K, and the Chan,
Kelbon & Krieger set) where the LEMTA case cannot: under pyrolysis, with the interior
heating at 10–30 K/min and vapour driven inward against a sealed back, against a dry run at
the same flux.

## Data

The paper prints only figures. `data/` holds them digitised from the 300 dpi page images
(frame fitted as four lines through the rotated scan; markers found by template matching in
a strip at each whole minute and traced by a Viterbi path; checked by overlays and by hand;
±3 K, ±0.003 in X_T):

| file | content |
|---|---|
| `flux.csv` | incident flux of curves A–E [kW/m²] on a 1-min grid (Fig. 2, 1-min medians, 3-point mean) |
| `flux_points.csv` | the digitised Fig. 2 points |
| `temperatures.csv` | `dry_E` (Fig. 3), `y09_E` (Fig. 11), `y20_E` (Fig. 12), `y20_E_repeat` (Fig. 4, a second 0.20 test), `y09_B` (Fig. 8, 19 mm only); `match` is the template score |
| `conversion.csv` | total conversion X_T (mass lost / initial mass, water included) of the first three (Fig. 7) |
| `table1.csv` | the paper's ΔT_max for A–E × 0/0.09/0.20; it does not match the plotted runs (card §2) and is not scored |

The dataset card is `docs/audits/validation/moist_wood/card_bilbao_1996.md` (local).

## Running

```sh
julia --project=scripts/verification scripts/bilbao_wet_wood/smoke.jl y20_E front=tc1 kinetics=larch shrinkage=0.2 k_scale=1.366
julia --project=scripts/verification scripts/bilbao_wet_wood/calibrate_dry.jl tc1
julia --project=scripts/verification scripts/bilbao_wet_wood/run_cases.jl tc1/sorption y09_E y20_E y20_E_repeat
julia --project=scripts/verification scripts/bilbao_wet_wood/run_cases.jl tc1/gated y09_E y20_E y20_E_repeat
julia --project=scripts/verification scripts/bilbao_wet_wood/compare.jl tc1 [family ...]
```

The transfer coefficient's ladder is three variants, `tc1/sorption/lambda1`, `lambda3` and
`lambda10`, each run on every run; the other tags (`variants.jl`) go the same way. A run
takes about 100 s (with the 1 mm layer's solve); variants can run side by side as separate
processes.

`calibrate_dry.jl FAMILY` fits the family's factor to the dry run and writes it to
`calibration.toml`; `run_cases.jl` takes a variant `family/law[/tag]` (`variants.jl`; the
`:tc1` families need runs with a 1 mm record, so not `y09_B`);
`compare.jl` writes `results/metrics.csv`, `results/tables/<family>.md`,
`<family>_phase_change.md` (the columns of the moisture model: rms, the 9 mm reading at the
1 mm record's 440 and 530 K marks, crossing times, pressure, p_v/p_eq, moisture content,
condensation, steps) and `<family>_plateau.md` (the 365→400 K plateau and the crossing times
at 9 and 19 mm), prints the last two, and draws `results/figures/` (`results/` is generated
and ignored).

## Scenario

- **Slab.** 1D, 19 mm across the grain, 100 equal cells; uniform start at 300 K with
  Y₀ρ_dry of bound water.
- **Wood** (`deck.toml`): ρ_dry 430 kg/m³ (measured); k the Flity et al. 2024 law scaled to
  the paper's measured 0.103 W/(m K); c the Flity et al. 2024 law; ε = α = 0.78 (the paper).
- **Pyrolysis.** Most families use the three pseudo-components fitted to larch TGA in N₂ (the
  LEMTA argon decks, Flity 2023), char yield 0.217 against the paper's 0.24; `tc1_bilbao` and
  `bilbao` use the Arrhenius part of the paper's own TGA kinetics (A = 8.31 × 10⁵ s⁻¹,
  E = 101.3 kJ/mol, char 0.24) with its DSC heats (+274 kJ/kg to X_s = 0.6, −353 after). The
  paper's kinetics add a low-temperature clamp and heating-rate terms that a first-order step
  cannot carry.
- **Char**: the LEMTA argon decks' laws (c, the pooled k law, ε = α = 0.95). A fully
  converted layer keeps 80 % of its thickness (`shrinkage = 0.2`, measured across the grain
  on other species, Flity 2023); the mesh moves with the solid and each thermocouple is read
  at its material point.
- **Water**: the LEMTA moist decks' water, vapour and conductivity convention. The law is the
  variant's:
  - `sorption` (the default), the moisture model: a `PhaseChange` of the water and its
    vapour with the three wood components as the sorbing solids (the char holds none),
    exchanging toward sorption equilibrium on the deck's isotherm (`CAST_SPRUCE`) at the
    rate `exchange_rate` (10³ s⁻¹), evaporating where the pore vapour is below
    equilibrium and condensing where it is above; its heat is built from the water's and
    the vapour's heat capacities and the isotherm. The pore vapour starts in equilibrium
    with the water and the N₂ makes up the ambient pressure.
  - `freewater`: the same step with the activity of free water (`FreeLiquid`, no isotherm).
  - the one-way laws of `LAWS` (`gated`, `chan`, `arrhenius`), with the LEMTA decks' heat of
    desorption: they evaporate whatever the vapour pressure and cannot condense. They run at
    the deck's transfer coefficient (5 × 10⁻⁷ m²/s) at every moisture content.
- **Gas**: volatiles, vapour and N₂; permeability 1e-15 m² (wood), 1e-11 m² (char); out
  through the exposed face to dry N₂ at 1 atm (the reservoir names the vapour at zero); the
  back is impermeable. The wood's transfer coefficient stands for bound-water diffusion
  across the grain: 5 × 10⁻⁷ m²/s, times 3 for a run at 15 % moisture or more (`deck.toml`,
  `wet_moisture`, `transfer_wet_scale`; the moisture diffusion coefficients of Siau 1984 and
  Keey et al. 2000). The 9 % run takes 5 × 10⁻⁷, the 20 % runs 1.5 × 10⁻⁶; the tags
  `lambda1`, `lambda3`, `lambda10` set the factor on 5 × 10⁻⁷ for every run. `tortuosity`
  (tags `dusty4`, `dusty12`) replaces the lumped transfer by the dusty-gas closure.
- **Top boundary**, by `front`:
  - `:tc1` (the comparison's main mode): the domain is the 18 mm below the 1 mm
    thermocouple, its top held at that thermocouple's measured record; no flux enters. The
    1 mm layer above is solved alongside at the same record for its mass (`front = :top1`).
  - `:flux`: the whole 19 mm; the exposed face absorbs `flux_factor` × the curve's flux at the
    surface absorptivity and re-radiates to 300 K; no convection (a stably stratified gap
    under a hotter panel; the authors' model has none).
- **Back**: adiabatic (the paper), or `h_back`.
- **Calibration on the dry run** (`calibrate_dry.jl`, 9 and 19 mm over 0–24 min,
  `calibration.toml`), per thermal family (`variants.jl`):
  - under `:tc1`, one factor on the solids' conductivity (wood and char): 1.37 (larch,
    shrinkage), 1.53 (the paper's single step), 1.48 (fixed mesh);
  - under `:flux`, the flux factor: 2.4–2.8. The measured flux cannot heat the slab as fast as
    the thermocouples and the mass show (the 1 mm thermocouple also reads above the flux's
    radiative-equilibrium temperature), and the factor does not carry between runs.
- **Scoring window**: a run is scored until its measured 19 mm record passes 626 K (the dry
  run's at 24 min, the end of the calibrated range), or to its end.

Results: `docs/audits/validation/moist_wood/Results_Bilbao.md` (local). In short, the phase
change brings the back face (19 mm) through 373 K at 10.9 and 13.9 min in the 9 and 20 %
runs (measured 12.7 and 12.9; the gated law 18.0 and 27.0) and the later 9 mm plateau to
the measured level (380 and 382 K at the 530 K mark of the 1 mm record, measured 375 and
378, gated 361–362), with the pore vapour at equilibrium and the moisture rising ahead of the
front to 1.3–1.4 times its start. The 9 mm rms, set after dry-out by the thermal and
pyrolysis model, is not improved, and the early heating of the 20 % run's interior (9 mm at
375 K by 8 min; 339–343 K under every variant of the phase change) is not reproduced.
Among the one-way laws gated <
Chan < Arrhenius in rms in every run, and the Arrhenius step chills the back 10–17 K below its
start.
