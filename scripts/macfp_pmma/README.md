# MaCFP PMMA comparison suite

Runs the MaCFP condensed-phase PMMA property sets in Pyrolysis.jl, FDS,
Gpyro and ThermaKin2Ds, and compares each run with the institution's own
submission, with the other codes' runs of the same deck, and with the NIST
gasification apparatus experiments. `report.jl` writes the validation pages: the
overview `docs/guide/validation/macfp_pmma.md` and its detail pages in
`docs/guide/validation/macfp_pmma/`.

- **MaCFP-2 targets**: TGA at 10 and 100 K/min; 1D anaerobic gasification
  at 10, 25 and 65 kW/m² on 6 and 12 mm slabs, no convection, perfectly
  insulated back, 293 K (guidelines v1.3).
- **MaCFP-3 NIST gasification apparatus**: the eleven PMMA tests at 25 and
  50 kW/m² and the two inert tests, against the experiments and the groups'
  own predictions.

Data come from a checkout of [MaCFP/matl-db](https://github.com/MaCFP/matl-db)
(MIT). Point `MATL_DB` at it; the default is `~/Documents/git/matl-db`.
`check_matl_db()` warns if its PMMA tree differs from the pinned commit
(`MATL_DB_COMMIT` in `MacfpPMMA.jl`).

## Layout

| path | what |
|---|---|
| `decks/<set>.toml` | one canonical, solver-neutral deck per property set |
| `MacfpPMMA.jl` | decks → Pyrolysis.jl problems; scenarios; MaCFP output quantities; comparison and cost metrics; tolerant reader for matl-db CSVs |
| `run_pyrolysis.jl` | Pyrolysis.jl on every deck × gasification case × variant, and the TGA targets |
| `run_apparatus.jl` | Pyrolysis.jl on every deck × apparatus test, and the inert tests |
| `external/` | deck writers, runners and extractors for FDS, Gpyro and ThermaKin2Ds |
| `compare.jl` | figures and metrics of the runs against the submissions and the apparatus experiments (`results/figures/`, `results/metrics.csv`, `results/apparatus_scores.csv`) |
| `report.jl` | the validation overview, its detail pages and their figures |
| `time_pyrolysis.jl`, `run_timing.py`, `timing_report.jl` | cost at matched accuracy on a quiet machine |
| `results/` | generated (ignored): `results/<code>/<set>/<variant>/*.csv`, `results/timing/` |

```sh
julia --project=scripts/macfp_pmma scripts/macfp_pmma/run_pyrolysis.jl [SET ...] [--cells=N]
julia --project=scripts/macfp_pmma scripts/macfp_pmma/run_apparatus.jl [SET ...] [--inert]
python3 scripts/macfp_pmma/run_timing.py scripts/macfp_pmma/results/timing --steps pyrolysis --sets NIST,Aalto_I,UMD,GIDAZE+
julia --project=scripts/macfp_pmma scripts/macfp_pmma/compare.jl
julia --project=scripts/macfp_pmma scripts/macfp_pmma/report.jl
```

The external codes' runs are described beside their runners
(`external/README_fds.md`, `external/README_gpyro_thermakin.md`).

Every code writes one CSV per case with the columns
`time_s,mass_g,mlr_g_m2s,T_top_K,T_back_K,T_1mm_K,thickness_mm` (mass for a
10 cm × 10 cm sample; `T_1mm_K` is 1 mm below the current surface, the
"top" temperature UMD and BUW-FZJ reported).

## Deck format

```toml
id = "NIST"                       # set name
calibration_code = "FDS"          # the code the group calibrated and predicted with
predictions = "NIST"              # file prefix in Computational_Results/Gasification_Predictions
tga_predictions = { 10 = "...", 100 = "..." }
absorption = "surface"            # or "in-depth"
kappa_kernel = "beer-lambert"     # in-depth only: the closure the group fitted κ in ("two-flux" for FDS); recorded, not read
order_basis = "own"               # n ≠ 1: "own" (normalised by the reactant's own mass), "init", "density"
top_temperature = "surface"       # the submission's "top": "surface", "avg_1mm" (UMD), "avg_3pt" (BUW-FZJ)
kinetics_only = false             # true: TGA only

[gas]
cp = 0.0                          # pyrolyzate cp the calibration code used (0: none)

[[components]]                    # condensed components; residues and intermediates start at y0 = 0
id = "PMMA"
rho = 1100.0                      # kg/m³
y0 = 1.0                          # initial mass fraction
emissivity = 0.9
kappa = 2870.0                    # 1/m, in-depth decks
k = { form = "constant", value = 0.20 }
cp = { form = "constant", value = 2200.0 }

[[reactions]]
reactant = "PMMA"
A = 2.85e13                       # 1/s, normalised form for n ≠ 1
E = 1.91e5                        # J/mol
n = 1.0
dH = 8.7e5                        # J per kg of reactant, positive endothermic
residue = ""                      # solid product (component id), the rest is gas
residue_yield = 0.0
```

Property forms: `constant {value}`, `linear {a, b}` (a + bT), `power {c0, Tr, n}`
(c0 (T/Tr)^n), `table {T, v}` (piecewise linear, held constant outside),
`switch {T, below, above}`. Temperatures are in K throughout.

Each deck's header and `notes` record the translation decisions that the
JSON leaves open (the κ unit and kernel, the n ≠ 1 normalisation, the
residue's properties, digits taken from a group's own input file) and the
defects found in the JSON and in matl-db's own decks of the set.

## Scenario choices (MaCFP-2 targets)

- The incident flux is absorbed at the surface mixture emissivity (Kirchhoff),
  re-radiation goes to 293 K, and there is no convection.
- A deck that absorbs at its surface runs under `SURFACE_ABSORPTION`. A deck
  with an absorption coefficient runs under `RADIATIVE_TRANSFER`: every cell
  absorbs and emits (a gray two-stream model), and the exposed face crosses
  radiation at the surface cell's emissivity both ways. A deck's κ is the
  material's absorption coefficient, read as written.
- The heater is `Diffuse()`. The MaCFP-2 targets state no distribution, and
  a diffuse and a collimated heater both satisfy "absorbed at the surface
  emissivity"; the collimated heater runs as a variant.
- The insulated back is a mirror: it returns the radiation reaching it into
  the sample, as FDS's insulated back does (`EMISSIVITY_BACK = 0`). Neither
  Gpyro nor ThermaKin2Ds puts that radiation into the sample: Gpyro drops what
  its march leaves after the last cell, and ThermaKin2Ds deposits none of the
  share that reaches the back. It matters once the sample is thinner than a
  few absorption lengths.
- Variants of the deck as written (`base`):

  | variant | decks | what it changes |
  |---|---|---|
  | `collimated` | in-depth | the heater as one beam along the normal, attenuated at κ rather than at the diffuse field's 1.5κ |
  | `fds_closure` | in-depth | every component's κ × ln(1 + 2κΔz)/(1.5κΔz) at Δz = 50 µm: FDS's two-flux attenuation on its own cells, which matches FDS's sweep cell by cell in transmission and emission |
  | `fds_matched` | UMD, UMET_TK, UMET_GP | `fds_closure` with the one item FDS cannot express aligned with FDS: gas cp 0 (FDS removes the pyrolyzate where it forms); UMET_GP's reaction 2 first order with A (ρ_B/M_ref)^(n−1), FDS's rate on a shrinking cell of pure B |
  | `surface` | in-depth | `SURFACE_ABSORPTION`: no radiation in depth |
  | `h5.33` | BUW-FZJ_A, BUW-FZJ_B | convection h = 5.33 W/(m² K) to 293 K at every flux, which their runs' steady states need |
  | `tablesC` | DBI_1, DBI_calibrated | the kelvin-valued cp and k tables read as °C, as their FDS RAMPs evidently read them |
  | `T303` | Aalto_I | start and surroundings at 303.15 K, where their runs start |

- Pyrolysis.jl gas: one MMA-like gas with the deck's gas cp. The cp is 0 when
  the calibration code carries no gas enthalpy; Gpyro's default cpg = 1000
  J/(kg K) enters only its gas–solid exchange, which changes the result by
  ≤ 0.8 %. Every solid has a pore fraction of 0.01 and a lumped transfer
  coefficient of 2e-5, so the gas leaves promptly.
- Pyrolysis.jl runs on a moving mesh of 20 cells/mm with deactivation. Its MLR
  is −dm/dt of the condensed mass. The interval containing a deactivation event
  is blanked, because the event vents the layer's last 1 % in one step. An
  in-depth deck without a residue ends at burnout (`min_thickness`) unless an
  n ≠ 1 tail outlasts the run (UMET_GP at 10 and 25 kW/m²). At n ≠ 1
  a run can end in a step failure once the reactant is spent: the TGA of the
  n ≠ 1 sets after the mass loss, and DBI's gasification and apparatus runs
  with at most 0.013 % of the mass left.
- The surface emissivity is the volume-weighted value of the surface cell.
  When a residue with a different ε accumulates there (BUW-FZJ_A, ε 0.407),
  the result depends on the cell size in every code and, in Pyrolysis.jl, on
  how long the moving mesh keeps a burnt-out layer in the surface cell:
  Pyrolysis.jl keeps the residue as a film in its Lagrangian surface cell, and
  FDS smears it over its first cell, with the effect converging only slowly as
  the cell shrinks.
- TGA: one 50 µm cell whose faces follow T = 300 K + βt.
- `order_basis` for n ≠ 1: A is in the normalised form, except for
  `"density"` (Sandia_1), where A multiplies ρ^n and is converted against
  ρ₀ = 1180 kg/m³.

## The apparatus (MaCFP-3)

- One shared deck for every set (`apparatus_problem`): the PMMA at the
  thickness that holds the test's measured mass at the deck's density on the
  6.99 cm disc, bonded to 28.52 mm of Kaowool PM (ρ 256, c 1070, the
  README's k(T)); convection h = 10 W/(m² K) to gas at 298 K; re-radiation to
  walls at 289 K; adiabatic beneath the board; the initial temperature of the
  test's first back-thermocouple reading (20.9 °C for the mass tests). Epoxy
  is ignored.
- The heater is `APPARATUS_HEATER = Cone(0, 50°)`: uniform intensity from the
  normal out to 50°. A coaxial disc fitted to the README's radial flux profile
  is seen out to 33–62° at 50 kW/m² within the profile's error bars; at
  25 kW/m² no disc fits within them. The in-depth decks also run the
  `diffuse` and `collimated` heaters, which bound it.
- The Kaowool cell under the PMMA is opaque (κ = 8e6 1/m) and absorbs what
  the PMMA transmits; the board beneath it is transparent (α = 0). Its
  measured k(T) already carries its own radiation, and opaque cells emitting
  at their centre temperature would add about 4σT³Δz of conductivity on the
  board's graded cells (up to about 2.7 mm wide).
- The inert tests (`inert_problem`): five black-painted Kaowool discs with
  thermocouples at 5.72, 11.44 and 17.16 mm, and a black-painted copper disc
  on Kaowool, both at 50 kW/m² on the shared deck's convection and walls.
  They test the shared h and the Kaowool's k.
- Scores: MaCFP's NRMS (matl-db `Scripts/NIST_Gasification_Validation.py`):
  the mean of the runs of the repeats against the mean of the measured
  repeats over the shortest repeat, RMS/(max − min) of the measured mean.
