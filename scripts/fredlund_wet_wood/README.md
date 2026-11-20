# Fredlund wet-wood cubes

Cubes of spruce and pine, 135 mm, heated from one face across the grain under a cone at
90 kW/m² nominal, in air, with thermocouples at 2.5–130 mm and gauge pore-pressure needles at
10–70 mm, at several moisture contents; and a 150 mm pine cube at 14 % with resistive
moisture sensors. Sources: B. Fredlund, *A model for heat and mass transfer in timber
structures during fire*, Report LUTVDG/(TVBB-1003), Lund University (1988), and *A computer
program for the analysis of timber structures exposed to fire*, Report LUTVDG/(TVBB-3020),
Lund (1985).

The case scores the moisture model of `scripts/lemta_moist_wood` (the phase change of the
bound water and its vapour toward sorption equilibrium) and the gated one-way law beside it
on a wood of low permeability (0.88 × 10⁻¹⁶ m²) at 14.5 and 24 % moisture: the time the
interior spends between 100 and 150 °C, the pore pressures, and the moisture rising ahead of
the drying front. Fredlund's own model (WOOD1, equilibrium drying of free water at the
saturation pressure) is not reproduced; the deck takes its inputs and is scored against the
measurements.

## Data

The reports print figures only (the raw data were kept on floppy discs, 1988 p. 193). The
1988 report is a phone scan. `data/` holds what was read by hand from 200 dpi renders,
zoomed four times, each window's axes mapped through its own grid lines:

| file | content |
|---|---|
| `temperatures.csv` | the measured thermocouples ('+') at 10–70 mm: spruce at 14.5 % (Fig. 9.2), oven-dry spruce (Fig. 9.8, which started at about 50 °C), pine at 24 % (Fig. 9.11); ±5 °C, ±0.3 min |
| `tc10.csv` | the 10 mm thermocouple of the three tests, the boundary of the column below it |
| `pressures.csv` | the peak gauge pressures at each needle depth and their times, spruce at 14.5 % (Fig. 9.3) and pine at 24 % (Fig. 9.12); ±0.3 mWG; the author's cracked records flagged |
| `moisture_1985.csv` | the 1985 cube's resistive moisture readings (Fig. 5.4) against its thermocouples (Fig. 5.3): the peak, its time, the local temperature then, the time the local 100 °C is reached; ±0.3 %, ±0.5 min. No temperature correction of the resistance is described: the peaks fall at 40–85 °C, where an uncorrected meter reads high |
| `flux.csv` | the corrected flux of the 1988 tests (Table 9.1) and the radiation meter's record of the 1985 test (Fig. 5.2) |

The 1988 figures run on the data clock; the exposure starts at about 8 min on it (Fig. 6.26).
The dataset cards are `docs/audits/validation/moist_wood/card_fredlund_1988_thesis.md` and
Card 1 of `card_fredlund_lee_bryden.md` (local).

## Running

```sh
julia --project=scripts/verification scripts/fredlund_wet_wood/smoke.jl spruce14 sorption
julia --project=scripts/verification scripts/fredlund_wet_wood/calibrate_dry.jl 0.6:4.5 0.7:3.5
julia --project=scripts/verification scripts/fredlund_wet_wood/calibrate_dry.jl pick
julia --project=scripts/verification scripts/fredlund_wet_wood/run_cases.jl sorption spruce14 pine24 pine14_1985 spruce_dry
julia --project=scripts/verification scripts/fredlund_wet_wood/run_cases.jl gated
julia --project=scripts/verification scripts/fredlund_wet_wood/compare.jl
```

`calibrate_dry.jl` solves the oven-dry spruce at pairs of factors (`k:kc`, the solids'
conductivity and the char's alone; pairs can run as separate processes) and appends their rms
to `results/calibration_scan.csv`; `pick` writes the best pair to `calibration.toml`, which
every variant then takes. `run_cases.jl` takes a variant `law[/tag]` (`variants.jl`; tags join
with `+`) and the cases (`spruce14`, `pine24`, `pine14_1985`; `spruce_dry` on request) and
writes `results/<variant>/<case>.csv`. A case takes 1–3 min. `compare.jl` writes
`results/metrics.csv` and `results/tables/fredlund.md` and prints the tables: per case the
rms and bias at each depth, the times 100, 150 and 200 °C are reached and the time from 100 to
150 °C; the peak gauge pressures and their times; for the 1985 cube the moisture rise at each
depth, its lead before the local 100 °C and the width of the raised zone (`results/` is
generated and ignored).

## Scenario

- **Column.** 1D across the grain, 0.5 mm cells, from the exposed face to the back of the
  cube; uniform start at 20 °C (the oven-dry cube at 50 °C) with the moisture's water.
- **Wood and char** (`deck.toml`): WOOD1's inputs (1988 pp. 219–233): ρ_dry 386 (spruce),
  495 (pine), 440 (oven-dry spruce); char 120 kg/m³ holding the wood's volume; k and c as
  tables against temperature; one pyrolysis step with Fredlund's in-situ kinetics
  (0.54 s⁻¹, 26.3 kJ/mol, no heat), gated off below 417 K, the lower end of the fit; the
  cell-wall density 1500 kg/m³ sets the pore fraction.
- **Conductivity, calibrated.** WOOD1's conductivities leave the oven-dry interior far too
  cold under its own 10 mm record. Fitted to it (`calibrate_dry.jl`): the wood's ×0.6 and the
  char's ×2.7 (×0.6 on both, ×4.5 on the char), 30 K rms over 20–60 mm; carried unchanged to
  the wet cubes (tag `uncal` drops them).
- **Water**, by the variant's law:
  - `sorption` (the default), the moisture model: a `PhaseChange` of the water and its vapour
    with the wood as the sorbing solid (the char holds none), toward equilibrium on
    `CAST_SPRUCE` at 10³ s⁻¹; the pore vapour starts in equilibrium and air makes up 1 atm.
  - `freewater`: the same with the activity of free water.
  - `gated`, the one-way laws of `LAWS`: they evaporate whatever the vapour pressure.
- **Gas transport.** Permeability 0.88 × 10⁻¹⁶ m² (wood), 5.3 × 10⁻¹² m² (char), WOOD1's
  (1988 eq. 7.40). The wood's transfer coefficient 5 × 10⁻⁷ m²/s, ×3 at 15 % moisture or
  more (the Bilbao deck's rule), with the ladder `lambda1`, `lambda3`, `lambda10`.
- **Top**, by the case's front:
  - `tc10` (the 1988 cubes): the column below the 10 mm thermocouple, its top held at that
    thermocouple's measured record and open to the air. The heating the deck cannot carry
    near the face (the flame in air, the char's oxidation and recession) is then in the
    boundary. The 10 mm pressure is not scored.
  - `flux` (the 1985 cube; tag `fluxfront` for the others): the whole column under the
    measured flux, absorbed at the surface absorptivity (wood 0.78, char 0.95), with
    re-radiation and convection (10 W/(m² K)) to 20 °C. Without flame or oxidation the
    interior heats about half as fast as measured.
- **Air** at 1 atm carrying the vapour that leaves the conditioned wood at rest (the test
  humidity is not recorded).
- **Back**: adiabatic, impermeable (tag `openback`: open to the air).

Results: `docs/audits/validation/moist_wood/Results_Fredlund.md` (local).
