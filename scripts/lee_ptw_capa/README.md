# Lee et al. 2025 treated-pine CAPA case

Disks of southern pine treated with micronized copper azole, 7 cm across and 9.0 mm thick,
dried density 498 kg/m³, conditioned to 6.2 % moisture (dry basis), heated from above in the
Controlled Atmosphere Pyrolysis Apparatus II under nitrogen at 29 and 63 kW/m², three tests
each. The bottom face carries a painted copper foil (ε 0.92) read by an IR camera; the mass is
weighed. Source: S. Lee, J.A. De Beer, P.B. Sunderland, M.J. Gollner, S.I. Stoliarov,
"Pyrolysis models for pressure treated wood and wood–plastic composite", *Fire and Materials*
49 (2025) 85–100.

At both fluxes the bottom face holds a shoulder at 375–385 K, which the authors' model, with
a one-way water step, misses by 30–50 K early. The case scores the moisture model (the phase
change of the bound water and its vapour) on that shoulder, with the paper's step beside it.

## Data

The paper prints figures only; the authors' records are not public.

| file | content |
|---|---|
| `data/fig6_bottom_temperature.csv` | bottom temperature, mean of three tests with its band (markers and bars) and the paper's model line, Fig. 6, both fluxes |
| `data/fig7_mlr.csv` | mass loss rate per area, the same, Fig. 7 |

Read from the 300 dpi figures embedded in the publisher's PDF: the frame and tick marks
located in the image, then in every pixel column the band of black pixels and the red model
line, binned to 5 s (63 kW/m²) and 10 s (29 kW/m²). Reading uncertainty ±2 K where the band
is thin; rows where the red line covers part of the band, or that sit in the axis margin, are
flagged and not scored. The red line at 63 kW/m² equals the authors' ThermaKin output to
about 1 K, which checks the axis mapping.

Two inputs of the deck are not in the paper: the boundary conditions of the authors'
ThermaKin2Ds condition file at 62.7 kW/m² (the nitrogen temperature history, the convection
coefficient, the radiant set point's position dependence, the back face's surroundings), and
the sidewall temperature of the two 2023 tests at 29 kW/m² from the rig's records (S. Lee,
personal communication; the records are not distributed). `deck.toml` gives each value
with its source.

The dataset card is `docs/audits/validation/moist_wood/card_lee_2025.md` (local).

## Running

```sh
julia --project=scripts/verification scripts/lee_ptw_capa/smoke.jl 63 sorption
julia --project=scripts/verification scripts/lee_ptw_capa/run_cases.jl sorption
julia --project=scripts/verification scripts/lee_ptw_capa/run_cases.jl thermakin
julia --project=scripts/verification scripts/lee_ptw_capa/run_cases.jl sorption/kappa10+lambda1e-7 29
julia --project=scripts/verification scripts/lee_ptw_capa/compare.jl
```

`run_cases.jl` takes a variant `law[/tag]` (`variants.jl`; tags join with `+`) and the runs
(`63`, `29`; both by default) and writes `results/<variant>/<run>.csv`. The transfer
coefficient's ladder is `sorption/lambda1e-7`, `lambda2e-7`, the deck's 5e-7, `lambda1.5e-6`
and `lambda2e-5`; the permeability's `kappa0.1`, `kappa10`, `kappa100`. A run takes 1–3 min;
variants can run side by side as separate processes. `compare.jl` writes
`results/metrics.csv`, `results/tables/<run>.md` (the rms over the record and over the
moisture window, the readings at marks, the times 340, 360, 375 and 400 K are reached, the
time between 365 and 395 K, the mass loss rate, the mass lost, the pore pressure, the bottom
cell's moisture content, p_v/p_eq and the steps), and `results/figures/<run>.png`, and prints
the tables (`results/` is generated and ignored). The paper's model is scored the same way
from the figures' red lines.

## Scenario

- **Column.** 1D, the 9.0 mm of wood over the copper foil (50 µm, ρ 8933, c 385, k 401,
  carried in the two bottom cells at half their volume), 180 cells; uniform start at 305.2 K
  with 6.2 % of water.
- **Material** (`deck.toml`): examples/08's treated pine (the ThermaKin components of the
  paper: five solids, four pyrolysis steps with the paper's kinetics and heats, their gases),
  with its parallel weight 0.5 for k, λ and κ. The water holds its own volume at 1000 kg/m³
  and adds to k by that share.
- **Water**, by the variant's law:
  - `sorption` (the default), the moisture model: a `PhaseChange` of the water and its vapour
    with PTW, PTW1, PTW2 and PTW3 as the sorbing solids (the char holds none), toward sorption
    equilibrium on `CAST_SPRUCE` at the exchange rate 10³ s⁻¹. The pore vapour starts in
    equilibrium and nitrogen makes up 1 atm.
  - `freewater`: the same with the activity of free water.
  - `thermakin`: the paper's one-way step (A 1.55e4 s⁻¹, E 43.5 kJ/mol, 2.78 MJ/kg, from
    oriented strand board): it evaporates whatever the vapour pressure and cannot condense.
- **Gas transport.** The stages that are still wood (PTW, PTW1) take the transfer coefficient
  of bound-water diffusion across the grain, 5e-7 m²/s (about 1e-7 is the measured size at
  6 %; examples/08's 2e-5 has no source); the later stages and the char 2e-5. Permeability
  1e-15 m² for the wood stages, 1e-12 for the char.
- **Exposed face** (63 kW/m², the condition file at the centre of its one radial layer):
  61.85 kW/m² incident, absorbed at the surface absorptivity; background 290 K; convection
  6.96 W/(m² K) to nitrogen at 333 K rising to 409 K. At 29 kW/m² the same factor, h and
  nitrogen history (tag `n2bg`: the nitrogen at 290 K). Gas leaves to dry nitrogen.
- **Back**: the foil at ε 0.92, impermeable, exchanging radiation and convection with the
  water-cooled sidewall below it: at 63 kW/m² the condition file's surroundings warming from
  295 to 329.6 K over 460 s with h ramped to 4 W/(m² K); at 29 kW/m² the measured sidewall
  (290.5 → 308.3 K) with h = 4.

Results: `docs/audits/validation/moist_wood/Results_Lee.md` (local). In short: the paper's
step on this deck reproduces the paper's model (each crossing within 6 s), so the boundary
conditions are those the authors used. The phase change brings the bottom face up 30–50 s
earlier than the step. At the deck's permeability, though, the water condensed at the back
boils off under 2–3 bar of pore pressure at 405–410 K, above the measured shoulder. With the
permeability ×10 it fits 29 kW/m² (rms 4.4 K, shoulder 118 s against 129) but heats the back
too late at 63 kW/m². There only a transfer coefficient of 2e-5 m²/s reaches 340 K by 35 s
(measured 34 s), and that value runs 36 s early at 29 kW/m². No constant setting fits both
fluxes.
