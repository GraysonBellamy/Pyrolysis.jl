# FDS 6.11 on the MaCFP-2 PMMA gasification targets

`fds.py` turns each canonical deck (`../decks/*.toml`) and each MaCFP-2 1D
gasification case (`GASIFICATION_CASES` in `../MacfpPMMA.jl`, read at run
time) into an FDS `SOLID_PHASE_ONLY` input, runs FDS, and writes the common
CSV `results/fds/<set>/<variant>/gasification_<q>kW_<L>mm.csv`
(`time_s,mass_g,mlr_g_m2s,T_top_K,T_back_K,T_1mm_K,thickness_mm`).
`fds_checks.py` reproduces every number below. Both are Python 3.11+,
standard library only.

FDS is the official Linux bundle `FDS-6.11.1-0-gff928db-release`
(`~/FDS/FDS6/bin/fds`, override with `FDS_EXE`). Source line references are
to the `FDS-6.11.1` tag of github.com/firemodels/fds.

## Running

```sh
cd scripts/macfp_pmma/external
python3 fds.py run                       # the nine MaCFP-2 sets x six cases, 6 FDS processes
python3 fds.py run NIST --cases=25x6     # one set, one case
python3 fds.py run --all                 # also the 2023 sets (BUW-FZJ_C, DBI_calibrated, NIST-StMU)
python3 fds.py run --tag=_c2.5e-5 --cell=2.5e-5   # refinement: results/fds/<set>/base_c2.5e-5/
python3 fds.py generate NIST --out=/tmp/x         # inputs only
python3 fds.py extract                   # re-extract CSVs from kept device files
```

`run` sets `OMP_NUM_THREADS=1`, `I_MPI_FABRICS=shm`, `LD_LIBRARY_PATH` and
`FI_PROVIDER_PATH` to the bundle's Intel MPI (without them FDS crashes in the
OFI provider), and pins each process to CPUs 0-19 (`--cpus`, avoids the
2.5 GHz cores 20-21). Each run directory `run_<q>kW_<L>mm/` keeps
`gasification.fds`, `gasification_devc.csv` and `gasification.out`, the last
trimmed to its header (version, solid grid, devices) and final timing summary
(the per-step diagnostics run to several MB); everything else FDS writes is
deleted. Every run appends a row to
`results/fds/<set>/runs_<variant>.csv`: wall time measured around the process
(includes the ~0.3 s MPI start), FDS's own time-stepping and total wall clock
(from the `.out`), its CPU time (`_cpu.csv`, `Total T_USED`), the FDS revision
and the first 16 hex digits of the deck file's SHA-256 (a changed deck shows
up as a changed hash).

## The scenario in FDS

MaCFP-2 guidelines v1.3: incident radiant flux q = 10, 25, 65 kW/m² absorbed at
the surface emissivity, re-radiation to 293 K, no convection, perfectly
insulated back, initial and ambient temperature 293 K, 6 and 12 mm, outputs for
a 10 cm x 10 cm sample.

- One 3x3x3 gas mesh, `SOLID_PHASE_ONLY=T`; the sample is a 0.1 m x 0.1 m
  `VENT` on the z = 0 boundary. Only its 1D wall is solved.
- **Flux.** `EXTERNAL_FLUX = q + σT∞⁴` with `SKIP_INRAD=T` and radiation on.
  FDS adds `EMISSIVITY x EXTERNAL_FLUX` to the incoming radiation
  (radi.f90:4703-4714; `SKIP_INRAD` drops the gas-phase part) and emits
  `EMISSIVITY σ T_s⁴` (wall.f90:2922, 2964), so the net radiative flux is
  ε q + ε σ (T∞⁴ − T_s⁴) exactly, with T∞ = 293 K and FDS's
  σ = 5.670373e-8 (Pyrolysis.jl: 5.670374419e-8, 2.5e-7 relative).
  For in-depth decks the same ε (q + σT∞⁴) is the incoming stream of FDS's
  two-flux model and the emission comes from the interior (below).
- **No convection.** `HEAT_TRANSFER_COEFFICIENT=0`; FDS uses any value ≥ 0 as
  given (func.f90:3016-3034).
- **Insulated back.** `BACKING='INSULATED'`: no conduction, convection or
  radiation through the back (wall.f90:2152-2157 set `Q_RAD_IN_B = 0` and
  `EMISSIVITY_BACK = 0`). For in-depth decks, radiation that reaches the back
  face is reflected: the two-flux outward stream starts from the full inward
  stream there (`RFLUX_UP = Q_RAD_IN_B + (1 − EMISSIVITY_BACK) RFLUX_UP`,
  wall.f90:2284). Nothing leaves through the back and nothing is deposited at
  the back face itself; the reflected part is absorbed on its way back through
  the slab or leaves through the top. This is FDS's native perfectly insulated
  back, and Pyrolysis.jl's insulated back is the same mirror.
- **Temperatures.** `TMPA = TMP_INNER = 19.85 °C` (293.0 K).
- **End times** are those of `GASIFICATION_CASES` (8000, 12000, 2000, 2400,
  1500, 2000 s for 10/6, 10/12, 25/6, 25/12, 65/6, 65/12).

## Translation, all sets

- One `MATL` per deck component: `DENSITY`, `EMISSIVITY`, `CONDUCTIVITY` or
  `CONDUCTIVITY_RAMP`, `SPECIFIC_HEAT` (kJ/kg/K) or `SPECIFIC_HEAT_RAMP`,
  `ABSORPTION_COEFFICIENT = kappa` for `absorption = "in-depth"` (none for
  surface absorption, FDS's default), and `ADJUST_H=F` on every material so the
  heat of reaction is the deck's constant (FDS ≥ 6.8 otherwise re-references it
  to a virtual-TGA peak with the gas and residue enthalpies).
- Each reaction on its reactant's `MATL`: `A`, `E`, `N_S = n`,
  `HEAT_OF_REACTION = dH/1000` (kJ/kg of reactant, positive endothermic in both
  conventions), residue `MATL_ID(1,j)`/`NU_MATL(1,j) = yield`, gas
  `SPEC_ID(1,j)='MMA'`/`NU_SPEC(1,j) = 1 − yield`. `E` is multiplied by
  R_FDS/R_Pyrolysis.jl = 8.314472/8.31446261815324 so E/R, the quantity in the
  rate, is the same in both codes (a 1.1e-6 change).
- The gas is one surrogate species, `&SPEC ID='MMA', FORMULA='C5H8O2'`. FDS
  needs a tracked species for the gas yield; MMA is not in its table, so FDS
  prints a warning. With `SOLID_PHASE_ONLY` and `ADJUST_H=F` its properties do
  not enter the solid. FDS has no gas flow inside the solid: pyrolyzate leaves
  the cell where it forms, at once.
- The `SURF` layer lists the virgin components with
  `MATL_MASS_FRACTION = y0`; FDS forms the layer density 1/Σ(Y/ρ), the same
  volume-additive density as `MacfpPMMA`. Residues and intermediates are added
  by FDS from the reactions.
- Mixing in FDS: k volume-fraction weighted, ρc mass weighted
  (wall.f90:2818-2832), κ volume weighted (wall.f90:2261-2273), surface ε
  volume weighted over the first cell. All components of a deck share k, c, ε
  and κ except the residues of BUW-FZJ_A/B/C, so the rule matters only there.
- Cells shrink as mass leaves (`ALLOW_SHRINKING` default), keeping
  Σρ_i/ρ̂_i = 1 (wall.f90:2390-2425), and FDS re-nodes the layer when a cell
  changes size by 25 % (`REMESH_RATIO` default). A layer thinner than 1 µm
  (`MINIMUM_LAYER_THICKNESS` default) is removed in one step (at most
  ~0.012 g of a 10 cm x 10 cm sample).

### Temperature-dependent properties

FDS tabulates each material's k and c at integer kelvin
(`ML%K_S(0:5000)`, read.f90:7587-7604) by evaluating the `RAMP`, and looks
them up at `INT(T)` (wall.f90:2820), so inside FDS a property is a staircase of
1 K steps whatever the RAMP. The generator therefore writes each T-dependent law
as a RAMP with a breakpoint at every integer kelvin from 250 to 1500 K (1251
points, T in °C) and `NUMBER_INTERPOLATION_POINTS=1250`, so that FDS's own
resampling (nearest node, func.f90:844) lands exactly on those breakpoints. The
entry at I holds the law at I + 0.5 K, the centre of the bin [I, I+1) that
FDS maps to it. The remaining error is the staircase itself: at most 0.5 K x
|dlaw/dT|, zero on average over each bin. Identical laws share one RAMP.

- Smooth laws (GIDAZE+ k and c power laws, UMET_GP c power law, UMET_TK c
  linear, UMD's linear pieces): ≤ 0.2 % (UMD c below 395 K ≤ 0.4 %).
- Tables (DBI_1, DBI_calibrated, BUW-FZJ_C): exact at the bin centres; the one
  bin holding a table breakpoint carries the kink (DBI_1 c ≤ 0.6 % at 312 K;
  BUW-FZJ_C k 4 % inside the 573-574 K bin). Held constant outside the table,
  as the deck says (DBI k beyond 337 K, c beyond 413 K).
- Switch (UMD k and c at 395 K): bin 394 holds the law below, bin 395 the law
  above, so FDS switches at exactly 395.0 K, reversibly, as the deck reads it.
- Non-positive laws are floored as `MacfpPMMA.positive` floors them: UMD's
  k = 0.34 − 4.2e-4 T above 395 K is floored at 1 % of k(300 K) = 0.0016 W/m/K
  from 806 K (only the char film at burnout gets there). FDS would replace a
  non-positive k by 10000 W/m/K (wall.f90:2831).

### Non-unit order (DBI_1, UMET_GP, DBI_calibrated)

FDS ≥ 6.7.8 evaluates r = A_F ρ_s^n e^(−E/RT) per unit cell volume with ρ_s
the component's partial density in its current cell (wall.f90:3463). The deck's
law (Pyrolysis.jl; Gpyro "own") is r = A M_ref (m/M_ref)^n at the material
point, M_ref the reactant's lineage inventory: Y0 ρ0 for a virgin component,
yield x the parent's inventory for an intermediate. The generator uses
A_F = A / M_ref^(n−1), which makes the two identical in an unshrunk cell:

| Set | Reaction | n | M_ref (kg/m³) | A (1/s) | A_F (FDS) |
|---|---|---|---|---|---|
| DBI_1 | PMMA → gas | 1.03 | 1181.5 | 6.31e8 | 5.1034e8 |
| UMET_GP | A → 0.962 B + gas | 3 | 1206.7 | 1e12 | 6.8675e5 |
| UMET_GP | B → gas | 1.8 | 0.962 x 1206.7 = 1160.8 | 2.51e17 | 8.8686e14 |

They part as the cell shrinks. A pure reactant's partial density stays at its
bulk density while it burns, so FDS's order acts on nothing and the reaction is
first order in conversion; the deck's (1 − α)^n is not reproduced. FDS cannot
express the normalised law on a shrinking grid (`ALLOW_SHRINKING=F` would
keep it, at the cost of a slab that never regresses). `fds_checks.py order`
integrates one shrinking FDS cell and the deck's law through a TGA:

| Set | β (K/min) | T_peak deck / FDS (K) | T at half mass loss deck / FDS (K) | peak rate deck / FDS (%/s) | max \|Δm/m0\| at equal T |
|---|---|---|---|---|---|
| DBI_1 | 10 | 633.7 / 633.8 | 627.3 / 627.1 | 0.259 / 0.263 (+1 %) | 0.008 |
| DBI_1 | 100 | 692.0 / 692.0 | 684.5 / 684.2 | 2.182 / 2.215 (+1 %) | 0.008 |
| UMET_GP | 10 | 634.7 / 634.8 | 633.8 / 629.7 | 0.312 / 0.432 (+38 %) | 0.173 |
| UMET_GP | 100 | 667.8 / 668.0 | 667.1 / 662.5 | 2.823 / 3.912 (+39 %) | 0.173 |

DBI_1's approximation is negligible. UMET_GP's is not: reaction 2 (n = 1.8)
runs first order in FDS, which raises the peak mass-loss rate by 38 % and
moves the half-mass temperature 4-5 K lower; reaction 1 (n = 3) keeps most of
its order because B fills the cell as A converts. In gasification the effect
shows as FDS burning out cleanly where the deck's n = 1.8 tail keeps a few
tenths of a percent of the mass for thousands of seconds (UMET_GP 10 kW/m²
6 mm: FDS below 0.1 % of the initial mass at 4137 s, Pyrolysis.jl not by
8000 s).

### In-depth absorption

FDS's in-depth radiation is a diffuse two-flux model (wall.f90:2260-2291):
each cell attenuates both streams with optical thickness 2κΔx, emits σT⁴ from
its own temperature, and the front face transmits ε of the incoming flux and
reflects 1 − ε of the outgoing stream back in. Every in-depth deck runs with
`ABSORPTION_COEFFICIENT = kappa` as written. That is the model the
FDS-calibrated decks are fitted in (Aalto_I, and the 2023 BUW-FZJ_C,
DBI_calibrated, NIST-StMU: `kappa_kernel = "two-flux"`); the decks whose κ is
a one-way Beer-Lambert coefficient (GIDAZE+ and UMET_GP from Gpyro, UMD and
UMET_TK from ThermaKin: `kappa_kernel = "beer-lambert"`) run the same way.

The implicit sweep crosses a cell with transmission 1/(1 + 2κΔx), so on its
50 µm cells FDS attenuates the diffuse stream at ln(1 + 2κΔx)/Δx, 1.72-1.80κ
for these decks (0.72κ for BUW-FZJ_C's residue). Pyrolysis.jl's diffuse two-stream model attenuates at 1.5κ
and integrates each cell exactly; its `fds_closure` variant takes
κ × ln(1 + 2κΔx)/(1.5κΔx) at Δx = 50 µm, which makes the two codes' cells
agree in transmission and emission, and `fds_matched` also aligns the items
below that FDS cannot express.

### Per-set notes

| Set | FDS materials and reactions | Properties | Absorption | Notes |
|---|---|---|---|---|
| NIST | PMMA → gas | constants | surface (ε 0.9) | |
| GIDAZE+ | PMMA → gas | k, c power laws → RAMPs | in-depth κ 2870 (Beer-Lambert in the deck) | |
| Aalto_I | PMMA_a, PMMA_b → gas (4/96 %) | constants | in-depth κ 3421.27 (two-flux, as calibrated) | |
| BUW-FZJ_A | PMMA_a, PMMA_b → 0.001 residue + gas | constants; residue k 0.263, c 1214.9, ε 0.407 | surface | residue ε dilution, below |
| BUW-FZJ_B | PMMA_a, PMMA_b → 0.001 residue + gas | constants | surface | E1 = 610.8 kJ/mol leaves PMMA_a (2.5 %) unreacted, kept as written |
| DBI_1 | PMMA → gas, n = 1.03 | k, c tables → RAMPs | surface (ε 1.0) | order approximation negligible |
| UMD | PMMA → 0.98 PMMA_int + gas → 0.002 char + gas | switch at 395 K → RAMPs; k floored above 806 K | in-depth κ 2870 (Beer-Lambert in the deck) | deck gas c_p 2000 and λ 2e-5 not expressible |
| UMET_GP | A → 0.962 B + gas (n 3); B → gas (n 1.8) | c power law → RAMP; ρ held at 1206.7 by the deck | in-depth κ 2700 (Beer-Lambert in the deck) | order approximation significant |
| UMET_TK | A → 0.962 B + gas; B → gas | c linear → RAMP; ρ held at 1204.2 by the deck | in-depth κ 2700 (Beer-Lambert in the deck) | deck gas c_p 2000 not expressible |
| BUW-FZJ_C (2023) | 4 virgin → 0.0101 residue + gas | k, c tables → RAMPs; residue κ 50000 | in-depth two-flux | |
| DBI_calibrated (2023) | PMMA → gas, n = 1.03 | k, c tables → RAMPs | in-depth κ 2250 two-flux | |
| NIST-StMU (2023) | PMMA1, PMMA2 → gas; PMMA3 inert (1 %) | constants | in-depth κ 2620 two-flux | ADJUST_H off here (the group ran 6.8 with it on) |

## What FDS cannot express exactly

1. **Non-unit order on a shrinking grid** (DBI_1, UMET_GP, DBI_calibrated):
   sized above. Negligible for DBI; significant for UMET_GP.
2. **Beer-Lambert absorption** (GIDAZE+, UMD, UMET_GP, UMET_TK): FDS only has
   its two-flux model, with interior emission and interface reflection; the
   deck's one-way κ is used as written.
3. **Gas inside the solid** (UMD, UMET_TK give gas c_p 2000 J/kg/K and a
   transfer coefficient 2e-5): FDS removes the pyrolyzate where it forms, so
   the gas carries no enthalpy through the hotter layers above it. With
   in-depth absorption the reaction zone is a fraction of a millimetre below a
   surface some 10-30 K hotter, so the omitted cooling is of order
   ṁ″ c_p,g ΔT ≈ 0.01-0.03 kg/m²/s x 2000 x 10-30 K ≈ 0.2-2 kW/m²
   (an estimate, not a run).
4. **Surface emissivity of an accumulating residue** (BUW-FZJ_A, ε_residue
   0.407): FDS's surface ε is the volume average over its first cell, so the
   few-µm residue film left by a 0.1 % yield is diluted by however much virgin
   material that cell holds. FDS's BUW-FZJ_A answer therefore depends on its
   cell size and is **not converged at 0.05 mm**: halving the cell moves the
   mass by up to 19 g and T_back by up to 138 K during burning, and the burning
   ends later (65 kW/m², 12 mm: 452 s at 0.05 mm, 496 s at 0.025 mm;
   25 kW/m², 12 mm: 1399 / 1454 s). Pyrolysis.jl keeps the residue as a film
   in its surface cell, so its answer depends on its cell size and on how long
   its moving mesh keeps a burnt-out layer in that cell. At 10 kW/m² the two
   codes agree (≤ 0.06 g at 6 mm, ≤ 1.9 g at 12 mm); at 25 and 65 kW/m² the
   mass gaps reach 10, 44, 28 and 55 g at 25/6, 25/12, 65/6 and 65/12
   (`fds_checks.py pyrolysis`). Treat both codes' BUW-FZJ_A curves at 25 and
   65 kW/m² as discretization dependent. BUW-FZJ_B (residue ε 0.989 ≈ virgin 0.990) is
   converged and agrees with Pyrolysis.jl to 0.2 g.
5. **Property staircase** (all RAMP decks): ≤ 0.2 % on smooth laws, above.
6. **Constants**: E/R preserved; σ differs by 2.5e-7. Both negligible.

Nothing else needed approximating: constant properties, n = 1 kinetics,
parallel and series reactions with residues, surface absorption, constant
heats of reaction, the 395 K switch and the held-constant table ends are
native in FDS. Temperature-dependent density (UMET_GP, UMET_TK) is already held
constant by the decks.

## Outputs

Devices on the vent (`TIME_AVERAGED=F`, instantaneous, `DT_DEVC = 1 s`):

| Column | From |
|---|---|
| `time_s` | device time (0, 1, 2, … s) |
| `mass_g` | `SURFACE DENSITY` [kg/m²] x 0.01 m² x 1000 |
| `mlr_g_m2s` | −d(SURFACE DENSITY)/dt x 1000, centred differences (one-sided at the ends), the `MacfpPMMA` definition. FDS's own `MASS FLUX` of MMA agrees to 5e-4 g/m²/s and integrates to the mass lost |
| `T_top_K` | `WALL TEMPERATURE` + 273.15 (front face; for in-depth decks with h = 0 this is the first cell's temperature) |
| `T_back_K` | `BACK WALL TEMPERATURE` + 273.15 |
| `T_1mm_K` | four `INSIDE WALL TEMPERATURE` devices at 1 mm ± 0.5 and ± 1.5 cells below the moving front (FDS measures `DEPTH` from the current surface, dump.f90:8931-8944), each with `INSIDE WALL DEPTH` giving its cell centre; linear interpolation to 1 mm. A wall thinner than 1 mm reports its back temperature (Pyrolysis.jl: the deepest cell centre) |
| `thickness_mm` | `WALL THICKNESS` x 1000 |

Once the wall has burnt away completely (thickness 0) FDS resets its
temperatures to the back gas temperature (wall.f90:2677-2682); those rows carry NaN
temperatures. Decks with a residue never reach zero thickness; their residue
film then sits at radiative equilibrium, ((q + σT∞⁴)/σ)^¼ (818 K at 25 kW/m²,
1036 K at 65 kW/m²), the same in both codes.

## Numerics and convergence

Base: uniform cells of 0.05 mm (`CELL_SIZE(1)=5e-5`, `STRETCH_FACTOR(1)=1`;
120 cells at 6 mm, 240 at 12 mm, coordinates listed in each `.out`), fixed
time step 0.01 s (`LOCK_TIME_STEP=T`, `WALL_INCREMENT=1`, the wall is updated
every step).

Largest change from base over the six cases, half cell (0.025 mm),
`python3 fds_checks.py refine <SET> base_c2.5e-5`. Mass over the whole run; the rest
to the end of the burning phase (the mass within 1 % of the initial mass of its
final value):

| Set | cases | mass (g) | MLR (g/m²/s) | T_top (K) | T_back (K) | T_1mm (K) |
|---|---|---|---|---|---|---|
| NIST | 6 | 0.103 | 0.07 | 1.80 | 2.76 | 2.76 |
| GIDAZE+ | 6 | 1.107 | 13.89 | 11.50 | 14.74 | 14.74 |
| Aalto_I | 6 | 0.434 | 3.16 | 3.32 | 2.86 | 2.86 |
| BUW-FZJ_A | 6 | 19.446 | 30.32 | 52.86 | 138.13 | 122.45 |
| BUW-FZJ_B | 6 | 0.090 | 0.34 | 0.64 | 1.59 | 1.59 |
| DBI_1 | 6 | 0.113 | 2.49 | 2.09 | 3.23 | 3.23 |
| UMD | 6 | 0.893 | 7.26 | 6.66 | 9.87 | 9.87 |
| UMET_GP | 6 | 0.285 | 1.86 | 2.09 | 1.28 | 1.27 |
| UMET_TK | 6 | 0.308 | 1.88 | 2.13 | 1.57 | 1.57 |

Largest change from base, half time step (0.005 s); cases 10/6, 25/12, 65/6;
same limits (`refine <SET> base_dt0.005`):

| Set | cases | mass (g) | MLR (g/m²/s) | T_top (K) | T_back (K) | T_1mm (K) |
|---|---|---|---|---|---|---|
| NIST | 3 | 0.050 | 0.05 | 0.71 | 1.16 | 1.16 |
| GIDAZE+ | 3 | 0.012 | 0.03 | 0.64 | 0.07 | 0.14 |
| Aalto_I | 3 | 0.001 | 0.02 | 0.37 | 0.02 | 0.07 |
| BUW-FZJ_A | 3 | 0.084 | 0.31 | 1.26 | 2.23 | 2.23 |
| BUW-FZJ_B | 3 | 0.018 | 0.15 | 0.40 | 0.51 | 0.51 |
| DBI_1 | 3 | 0.054 | 0.23 | 1.77 | 2.34 | 2.34 |
| UMD | 3 | 0.002 | 0.02 | 0.80 | 0.04 | 0.13 |
| UMET_GP | 3 | 0.001 | 0.01 | 0.21 | 0.01 | 0.08 |
| UMET_TK | 3 | 0.012 | 0.01 | 0.21 | 0.01 | 0.08 |

NIST at three cell sizes (dt 0.01 s), change from base (same limits):

| Case | cells (mm) | mass (g) | T_top (K) | T_back (K) |
|---|---|---|---|---|
| 65kW_6mm | 0.1 | 0.112 | 1.30 / 2.80 | 0.63 / 2.24 |
| 25kW_12mm | 0.1 | 0.051 | 0.27 / 0.60 | 0.09 / 0.72 |
| 25kW_12mm | 0.025 | 0.013 | 0.05 / 0.15 | 0.02 / 0.18 |

## NIST sanity check

NIST calibrated its set in FDS and submitted FDS predictions, so FDS should
reproduce them. `fds_checks.py nist` (temperatures to 20 s before the earlier
burnout, burnout = first time below 0.1 % of the initial mass):

| Case (kW/m², mm) | burnout NIST / FDS (s) | mass max / RMS (g) | T_top max / RMS (K) | T_back max / RMS (K) |
|---|---|---|---|---|
| 10, 6 | 3853 / 3904 | 0.62 / 0.39 | 4.6 / 0.76 | 2.3 / 0.78 |
| 10, 12 | 7268 / 7348 | 1.34 / 0.79 | 4.6 / 0.96 | 4.4 / 1.02 |
| 25, 6 | 767 / 768 | 0.13 / 0.07 | 11.6 / 0.77 | 2.1 / 0.54 |
| 25, 12 | 1564 / 1567 | 0.41 / 0.21 | 11.6 / 0.59 | 3.9 / 0.46 |
| 65, 6 | 244 / 244 | 0.30 / 0.20 | 30.1 / 5.39 | 0.7 / 0.54 |
| 65, 12 | 490 / 489 | 0.25 / 0.20 | 30.1 / 5.06 | 0.4 / 0.23 |

At 25 and 65 kW/m² mass and burnout agree within 0.4 g and 3 s. At
25 kW/m² T_top agrees to 0.8 K RMS; at 65 kW/m² NIST's surface is cooler in
the first 10 s (360 K against 390 K at 1.1 s) and then about 5 K hotter while
burning (RMS 4.7 K after 20 s), with the same mass loss: the signature of a
coarser surface cell in their run (here, going from 0.1 to 0.05 mm cells moves
T_top by 1.6 K at 65 kW/m²; inferred: their grid is not published). The other
large differences are in the last seconds before burnout. At 10 kW/m² NIST's run burns out 1.1–1.3 % earlier, a mass gap
growing linearly to 1.3 g with temperatures 0.6 K lower: about 0.5 % more net
heat. At 10 kW/m² the net flux is only ~2.4 kW/m² after re-radiation, so a
0.5 % difference in how NIST's run realised the radiative boundary (e.g.
re-radiation exchanged through the FDS radiation solver with the enclosure
instead of `SKIP_INRAD`; inferred: their input is not published) suffices.
FDS here and Pyrolysis.jl agree to 0.004 g and 0.03 K (T_back) on the same
case.

## Comparison with the groups' own submissions and with Pyrolysis.jl

FDS base against each group's own MaCFP-2 submission (`fds_checks.py submission`);
worst case over the six cases. Mass over the common time span; temperatures to 20 s
before the earlier burnout (below 0.1 % of the initial mass). The submitted 'top'
column is compared with FDS's T_1mm for BUW-FZJ_A/B and UMD, whose top is a mean over
the top millimetre, and with T_top otherwise.

| Set | max \|Δmass\| (g) | T_top RMS (K) | T_back RMS (K) |
|---|---|---|---|
| Aalto_I | 5.0 | 42.5 | 16.6 |
| BUW-FZJ_A | 67.2 | 171.9 | 189.1 |
| BUW-FZJ_B | 51.6 | 83.1 | 62.5 |
| DBI_1 | 17.6 | 20.6 | 38.5 |
| GIDAZE+ | 7.6 | 13.5 | 199.2 |
| NIST | 1.3 | 5.4 | 1.0 |
| UMD | 16.9 | 129.4 | 94.2 |
| UMET_GP | 19.0 | 20.2 | 20.4 |
| UMET_TK | 29.7 | 8.5 | 28.5 |

Notes on the submissions: Aalto starts at 303.15 K, not 293 K; BUW-FZJ_A/B
and DBI at 293.15 K; BUW-FZJ_B's runs keep 2.6 % unreacted mass and never burn
out. GIDAZE+'s submitted 12 mm back temperatures equal its top temperatures
from mid-test on, with half the mass left (e.g. 65 kW/m² from 285 s), so that
back-temperature statistic says nothing about FDS.

FDS base against Pyrolysis.jl base (`results/pyrolysis/<set>/base`,
`fds_checks.py pyrolysis`); worst case over the six cases. Mass over the
common time span; temperatures to 10 s before the earlier burnout (below
0.1 % of the initial mass) or the end of the shorter run.

| Set | max \|Δmass\| (g) | T_top RMS (K) | T_back RMS (K) |
|---|---|---|---|
| Aalto_I | 1.25 | 2.1 | 3.8 |
| BUW-FZJ_A | 54.93 | 101.4 | 154.2 |
| BUW-FZJ_B | 0.17 | 0.6 | 0.6 |
| DBI_1 | 0.81 | 1.0 | 1.0 |
| GIDAZE+ | 4.00 | 3.6 | 11.1 |
| NIST | 0.18 | 0.6 | 0.3 |
| UMD | 3.34 | 15.7 | 16.1 |
| UMET_GP | 14.35 | 20.1 | 5.4 |
| UMET_TK | 3.03 | 3.0 | 2.8 |

On the surface decks the two codes solve the same problem (BUW-FZJ_A aside,
above). On the in-depth decks the differences are FDS's cell attenuation and
the items FDS cannot express (UMD and UMET_TK's pyrolyzate enthalpy, UMET_GP's
order on a shrinking cell; the UMD and UMET_GP windows include their long
residue tails); Pyrolysis.jl's `fds_closure` and `fds_matched` variants align
them (above).

## Runs and timings

| Variant | runs | successful | sum of wall times (s) | sum of FDS CPU (s) | longest run (s) |
|---|---|---|---|---|---|
| base | 72 | 72 | 2548 | 2521 | 152 |
| base_c1e-4 | 2 | 2 | 22 | 21 | 15 |
| base_c2.5e-5 | 54 | 54 | 2544 | 2524 | 171 |
| base_dt0.005 | 27 | 27 | 1548 | 1538 | 126 |

## Files

- `fds.py`: generator, runner, extractor.
- `fds_checks.py`: the checks above (`order`, `nist`, `submission`, `refine`,
  `pyrolysis`, `runs`, `stale`).
- `results/fds/<set>/base/gasification_<q>kW_<L>mm.csv`: the nine MaCFP-2
  sets x six cases, and the 2023 sets BUW-FZJ_C, DBI_calibrated, NIST-StMU.
- `results/fds/<set>/base_c2.5e-5/`, `base_dt0.005/` (and `NIST/base_c1e-4/`):
  refinement runs.
- `results/fds/<set>/<variant>/run_<q>kW_<L>mm/`: the FDS input, device CSV
  and `.out` of each run; `results/fds/<set>/runs_<variant>.csv`: the run log.
