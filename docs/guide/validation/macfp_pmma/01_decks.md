# Decks and scenario

The property sets, the scenario and the variants behind the [MaCFP PMMA validation](../macfp_pmma.md).

Data come from [MaCFP/matl-db](https://github.com/MaCFP/matl-db) at commit
`2772b0eb`. The decks, the runners and the generator of these pages are
in `scripts/macfp_pmma/` (its `README.md` describes each choice); the translation of each
deck for FDS, Gpyro and ThermaKin2Ds is described beside their runners in
`scripts/macfp_pmma/external/`.

## Scenario

The MaCFP-2 gasification targets expose a 6 or 12 mm slab at 293 K to 10, 25 or
65 kW/m² with no convection; the slab re-radiates to 293 K and its back is perfectly
insulated. The flux is absorbed at the surface emissivity. Pyrolysis.jl runs each case
on a moving mesh of 20 cells/mm with deactivation, at the default tolerances. A deck
that absorbs at its surface runs under `SURFACE_ABSORPTION`. A deck with an absorption
coefficient κ runs under `RADIATIVE_TRANSFER`: every cell absorbs and emits, the
heater is `Diffuse()`, and the insulated back is a mirror. The TGA targets heat one
50 µm cell from 300 K at 10 and 100 K/min.

## Property sets

| set | year | calibrated in | absorption | κ (1/m) | reactions | MaCFP-2 submission |
|---|---|---|---|---|---|---|
| Aalto_I | 2021 | FDS | in-depth | 3421 | 2 | yes |
| BUW-FZJ_A | 2021 | FDS | surface | – | 2 | yes |
| BUW-FZJ_B | 2021 | FDS | surface | – | 2 | yes |
| BUW-FZJ_C | 2023 | FDS | in-depth | 3434, 50000 | 4 | no |
| DBI_1 | 2021 | FDS | surface | – | 1 | yes |
| DBI_calibrated | 2023 | FDS | in-depth | 2250 | 1 | no |
| GIDAZE+ | 2021 | Gpyro | in-depth | 2870 | 1 | yes |
| NIST-StMU | 2023 | FDS | in-depth | 2620 | 2 | no |
| NIST | 2021 | FDS | surface | – | 1 | yes |
| UMD | 2021 | ThermaKin | in-depth | 2870 | 2 | yes |
| UMET_GP | 2021 | Gpyro | in-depth | 2700 | 2 | yes |
| UMET_TK | 2021 | ThermaKin | in-depth | 2700 | 2 | yes |

The sets Sandia_1, Sandia_2, Sandia_3, Sandia_4, Sandia_5, Sandia_6, UCLAN give kinetics only and run the TGA targets alone.

FDS 6.11.1 runs all 12 sets,
Gpyro V2 runs Aalto_I, BUW-FZJ_A, BUW-FZJ_B, GIDAZE+, NIST, UMET_GP, UMET_TK, and ThermaKin2Ds 7 runs Aalto_I, BUW-FZJ_A, BUW-FZJ_B, GIDAZE+, NIST, UMD, UMET_TK in its
`MAX` mode (the absorbed radiation deposited in the element that absorbs most) and
Aalto_I, GIDAZE+, UMD, UMET_TK also in its
`RAND` mode (in an element drawn from the absorption distribution).

## Variants

Each in-depth deck also runs variants that move one choice, and three variants run
what a group evidently ran:

| variant | what it changes |
|---|---|
| `collimated` | the heater as one beam along the normal, attenuated at κ rather than the diffuse field's 1.5κ |
| `fds_closure` | every component's κ × ln(1 + 2κΔz)/(1.5κΔz) at Δz = 50 µm: FDS's two-flux attenuation on its own cells |
| `fds_matched` | `fds_closure` with the one item FDS cannot express aligned with FDS: no pyrolyzate enthalpy in the solid (UMD, UMET_TK); UMET_GP's second reaction first order, as FDS runs it on a shrinking cell |
| `surface` | `SURFACE_ABSORPTION`: no radiation in depth |
| `h5.33` | BUW-FZJ_A/B: convection h = 5.33 W/(m²·K) to 293 K, which their submissions include |
| `tablesC` | DBI_1, DBI_calibrated: the kelvin-valued c and k tables read as °C, as DBI_1's FDS runs read them (DBI_calibrated's are not published) |
| `T303` | Aalto_I: start and surroundings at 303.15 K, where their runs start |
