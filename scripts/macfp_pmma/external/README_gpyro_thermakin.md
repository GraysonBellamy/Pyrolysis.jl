# Gpyro and ThermaKin2Ds on the MaCFP-2 PMMA gasification targets

Native Gpyro V2 and ThermaKin2Ds 7 runs of the MaCFP-2 PMMA 1D gasification
targets for every property set each code can express, generated from the
canonical decks (`../decks/*.toml`) and the scenario in `../MacfpPMMA.jl`
(`GASIFICATION_CASES`, read at run time). Python 3.12+, standard library only.

```sh
# ThermaKin: standard build (MAX runs) + a seeded copy (RAND runs), both in ignored docs/refs/
test/verification/reference/matched/build_thermakin.sh
scripts/macfp_pmma/external/build_thermakin_seeded.sh
python3 scripts/macfp_pmma/external/thermakin.py run --variants base,rand --jobs 3
python3 scripts/macfp_pmma/external/gpyro.py run --jobs 1
# refinement and seed studies (variant directory = <variant><suffix>)
python3 scripts/macfp_pmma/external/thermakin.py run --sets UMD,NIST --cases 25_6 --variants base,rand \
    --dx 2.5e-5 --dt 0.005 --suffix _half
python3 scripts/macfp_pmma/external/thermakin.py run --sets UMD --cases 25_6 --variants rand --seed 2 --suffix _seed2
python3 scripts/macfp_pmma/external/gpyro.py run --cases 25_6 --dz 2.5e-5 --dt 0.01 --suffix _half
# checks
python3 scripts/macfp_pmma/external/sanity.py --table   # against the institutions' own submissions
python3 scripts/macfp_pmma/external/sanity.py --runs    # every run: wall time, end, flags
```

`common.py` holds what both drivers share (deck loading, the case list,
post-processing); `gpyro.py` and `thermakin.py` write the native inputs, run
them (at most `--jobs` processes, one thread each) and extract.

## Outputs

Everything goes to the ignored `../results/` tree:

* `results/{gpyro,thermakin}/<set>/<variant>/gasification_<q>kW_<L>mm.csv`,
  header exactly `time_s,mass_g,mlr_g_m2s,T_top_K,T_back_K,T_1mm_K,thickness_mm`,
  one row per second of simulated time;
* `results/{gpyro,thermakin}/<set>/<variant>/run_<q>kW_<L>mm/`: the native
  inputs (`case.data`; `case.cmp`, `case.cnd`), the native outputs (Gpyro
  point summary `case_summary_01_0001.csv`, gzipped temperature and cell-size
  profiles, `gpyro_case_0001.out`; ThermaKin `case.out.gz`),
  `run.log`, `meta.json` (wall time, last output time, flags) and `aux.csv`
  (native extras: Gpyro's instantaneous MLR; ThermaKin's top-boundary gas
  outflow, retained gas and three candidate readings of UMD's top-temperature
  column).

Variants: `base` = Gpyro as generated, ThermaKin with its default MAX
absorption; `rand` = ThermaKin RAND absorption with seed 1 (in-depth decks
only: for surface decks both modes select the top element); `_half` = half
cell/element size and half time step; `rand_seed2`, `rand_seed3` = other seeds.

Column definitions, identical for both codes:

* `mass_g`: condensed sample mass per m² × 0.01 m² (10 cm × 10 cm), in g.
  ThermaKin: printed total mass minus the backing (below) and the pyrolyzate
  held in the elements. Gpyro: its `TOTAL_MASS` point output (no gas solve).
* `mlr_g_m2s`: minus the centred difference (one-sided at the ends) of the
  areal condensed mass over the 1 s outputs, in g/(m² s) — the definition
  `MacfpPMMA.jl` uses for Pyrolysis.jl. The native instantaneous rates are in
  `aux.csv`: Gpyro's `MLR` agrees with it to 1e-4 g/(m² s); ThermaKin's top
  gas outflow agrees to 0.2 % of the peak rms while burning, with single-output
  spikes of up to 6 % when elements are merged.
* `T_top_K`: Gpyro node 1, which sits on the exposed face; ThermaKin's top
  element, whose centre is half an element (≤ 25 µm) below the face.
* `T_back_K`: Gpyro's last node (on the back face); ThermaKin interpolated at
  the sample/backing interface.
* `T_1mm_K`: interpolated 1 mm below the current surface from the node
  depths (Gpyro: rebuilt from the dumped cell sizes `DLTZN`; ThermaKin: the
  printed element centres); the back-face value once the sample is thinner
  than 1 mm (as `MacfpPMMA.jl`).
* `thickness_mm`: Gpyro `THICKNESS` (sum of cell sizes); ThermaKin printed
  total minus the backing.
* After burnout (sample mass below 1e-5 of the initial) there is no sample:
  the temperatures are NaN. Gpyro ends its own run at burnout (below).

## Scenario mapping (both codes)

| Scenario item | Gpyro V2 | ThermaKin2Ds 7 |
|---|---|---|
| Incident flux q, absorbed at ε | `QE = q`; Gpyro absorbs `EMIS·QE` | external radiation `I = q + σT∞⁴`; ThermaKin deposits `ε(I − σT⁴)` at the absorbing element |
| Re-radiation to 293 K | `RERADIATION=.TRUE., TINF=293` (σ = 5.67e-8 hard-coded) | the `σT∞⁴` added to `I` (native σ = 5.6704e-8) |
| No convection | `HC=0` | `CONVECTION COEFF: 0` |
| Insulated, impermeable back | `QE=0, HC=0, RERADIATION=.FALSE.` | two-element non-emitting backing (below), `MASS TRANSPORT: NO` |
| T0 = 293 K, L = 6 / 12 mm | `TMP_INITIAL`, `ZDIM` | layer `TEMPERATURE`, `THICKNESS` |
| Activation energy | E/R kept: `E[kJ/mol] = E·8.314e-3/8.314462618` (Gpyro's R0) | E/R kept with ThermaKin's R = 8.31447 |
| Heat of reaction, J/kg of reactant, endothermic + | `DHS = DHV = dH` (per kg reactant, eq. 45) | `HEAT: −dH` (exothermic +; reactant stoichiometry 1) |
| Pyrolyzate | no gas solve (instant escape) | component `gas`, deck `[gas]` cp if > 0 else 1 J/(kg K); λ = 2e-5 m²/s on every condensed component (deck `D` if given); top boundary `gas LIN 0.05 0`; gases take no volume |
| Mixing of k, ε, κ | volume fractions (native) | `PARALL CONDUCTIVITY: 1` (volume-weighted, as FDS, Gpyro, Pyrolysis.jl); ε volume-weighted; α_i = κ_i/ρ_i |
| Surface absorption (NIST, BUW-FZJ_A/B) | `KAPPA = 9D9` (≥ 1e6 m⁻¹: the net flux goes to the surface node) | α = 1e6 m²/kg: the top element absorbs all |
| Output | 1 s point and profile dumps, `TSTOP = t_end + 2·dt` (dumps skip TSTOP) | 1 s output, `DURATION = t_end − dt/2`, every element |

**Radiation that reaches the back face.** Neither code puts it into the
sample. Gpyro V2 marches Beer–Lambert from the surface cell inward
(`IN_DEPTH_RADIATION_ABSORPTION`, `gpyro_thermal.f90`) and drops what is left
after the last cell. ThermaKin deposits the radiation of a step in one
element, chosen from the distribution of what each element absorbs (the
reference manual's maximum and random absorption algorithms); the share that
the sample transmits, and that the ε = 0 backing would absorb, is never
deposited, and the deposit and its re-radiation scale with the fraction the
sample absorbs. Both matter only once the sample is thinner than a few
absorption lengths (in-depth decks: the last ~0.5–1 mm), where they slow the
final burnout and keep the remnant cooler than a reflecting back would.
Pyrolysis.jl's insulated back is a mirror: it returns that radiation into
the sample.

## Gpyro V2

Executable `docs/refs/gpyro_build/gpyro`, revision
`bff0126c68f2a35c084ef04490ae327185f281a3` (the matched-comparison build; its
banner hash is checked), one OpenMP thread. Generated inputs follow the
harness's parametric deck (`test/verification/reference/matched/run_comparisons.py`)
with these settings: half cells at both faces (node 1 on the exposed face),
`SOLVER_DEFORMATION_MODE=2` (cells shrink with their mass and the conduction
geometry is rebuilt every step; V2's default), no pressure/gas solve,
`THERMAL_EQUILIBRIUM`, `CONSTANT_DHVOL`, `SHYI_CORRECTION=.TRUE.` (the proper
enthalpy-form conduction term; it is zero unless species of different cp
coexist, i.e. only BUW-FZJ), tolerances `TMPTOL=1D-6, HTOL=1D-10, YITOL=1D-8`,
`EPS=1D-7`.

Gpyro specifics that shape the series:

* **In-depth absorption deposits the net flux.** In V2 (`BC_FOR_HALF_CELLS`)
  the surface net flux `ε·QE − εσ(T₁⁴ − T∞⁴)` is what is spread by
  Beer–Lambert when positive and κ_surface < 1e6 m⁻¹; the technical
  reference's eq. 66–67 (incident flux in depth, emission at the surface)
  describe a scheme V2 does not use.
* **Burnout.** V2 keeps a fully reacted cell with a residual mass `EPS` of its
  initial (no reactions, near-zero size) instead of removing it. It ends a run
  itself when the thickness falls below 5e-4 of the initial
  (`gpyro_main.f90`), or, when the last cell empties within one step, with its
  "Mass in cell has dropped to zero" shutdown; `gpyro.py` accepts the latter
  as burnout when less than 0.5 % of the mass is left or the rest would go in
  3 s at the last MLR (recorded in `meta.json`). Series of fully gasifying
  sets therefore end at burnout; a final dump that V2 writes after jumping its
  clock past `TSTOP` is dropped.
* **Power laws only**, `φ0 (T/Tr)^n` with one global `TREF = 300 K`; every
  power-law deck uses Tr = 300 K.
* **Solid yield = 1 + (ρ_B/ρ_A − 1)χ.** A residue of yield y gets
  χ = (1 − y)/(1 − ρ_B/ρ_A), which keeps both the deck's yield and the
  residue's own density (volume): BUW-FZJ_A χ = 14.15, BUW-FZJ_B χ = 3.299
  (χ > 1 is accepted). Equal reactant and product densities (UMET's A → 0.962 B)
  give no gas at all, so B's density is set 1e-4 below A's and χ = 380.
* **Reaction order.** `CONVENTIONAL_RXN_ORDER=.FALSE.`: Gpyro's own
  normalisation by the species' cumulative formed mass (`order_basis = "own"`,
  UMET_GP's calibration); identical to every other reading for n = 1.

## ThermaKin2Ds 7

Built from the non-redistributable source in `docs/refs/ThermaKin2Ds Code/`
(ignored; nothing from it is tracked). MAX runs use
`docs/refs/thermakin_build/thermakin` (`build_thermakin.sh`, banner-only
portability patch). RAND runs use `docs/refs/thermakin_build_seeded/thermakin`
from `build_thermakin_seeded.sh`, which applies the same banner patch and
replaces the clock seed of the random absorption algorithm by a seed read
from the environment variable `THERMAKIN_SEED` (unset: the clock seed); the
numerics are untouched and the source hash recorded in its `build_info.txt`
equals the standard build's. Element size 5e-5 m, time step 0.01 s (UMD's own
settings, Fiola et al. 2021 §3), Crank–Nicolson, output every second.

ThermaKin specifics that shape the series:

* **The back face always radiates.** The bottom boundary is defined like the
  top one and re-radiates at the emissivity of an element it selects, with
  or without external radiation; the reference manual's cone example puts a
  zero-emissivity layer at the bottom to prevent that loss. To make the back
  adiabatic here, two backing elements (thickness Δx each) sit under the sample: a
  non-reacting thermal copy of the first virgin component (same cp and k laws,
  glass transition included) at a tenth of its density (a tenth of its heat
  capacity per volume, 0.17 %/0.08 % of the 6/12 mm sample's), ε = 0, gas
  transfer 1e-30. Their absorption (in-depth decks: optical depth set so the
  bottom element absorbs 1.05·s/(1+s), s = one sample element's share) wins
  the bottom search against the sample element above, so the bottom
  boundary's emitter is the ε = 0 element. Two elements because ThermaKin
  merges a top element thinner than 0.2·Δx into the one below: a last residue
  layer then mixes into the upper backing element, never into the emitter.
  The backing is excluded from mass and thickness. It costs this: once an
  in-depth sample is about one element thick, MAX selects the backing from
  the top and deposits nothing, so in MAX runs of in-depth decks the last
  ≈ 0.02–0.03 mm (up to 0.4 g) lingers instead of burning out. RAND draws the
  backing only with probability transmittance × 0.12 and leaves ≤ 0.12 g at
  10 kW/m², nothing at 25 and 65 kW/m².
* **MAX vs RAND.** MAX (default) puts the whole net radiation into the one
  element with the largest Beer–Lambert share — for a uniform sample the top
  element — so in-depth decks behave as surface absorbers with the top
  element's temperature emitting. RAND draws that element every step from the
  Beer–Lambert distribution (seed fixed here), which on average reproduces
  in-depth absorption; re-radiation uses the drawn element's temperature.
* **Glass transition (UMD).** The deck's reversible 395 K switch of cp and k
  is reproduced as UMD implemented it (Fiola et al. 2021 §3.2; the ThermaKin
  reference manual's PMMA example): `PMMA_glass` (below-395 K laws) →
  `PMMA` (above-395 K laws), A = 100 s⁻¹, E = 0, h = 0, `TEMP LIMIT: L 395`;
  `PMMA_int` and `char` carry the above-395 K laws. Below 395 K the glassy
  polymer does not decompose (the deck's reaction 1 there: < 0.5 % conversion
  in 3600 s at 380 K, 2 % of it gas).
* **Property forms**: `p0 + p1·T + p2·T^n`: constant, linear (UMD, UMET_TK cp)
  and power laws (GIDAZE+: k = 0.2·300^1.2·T^−1.2, cp = 1200·300^−1.2·T^1.2)
  are exact. Densities are constant (the decks hold ρ constant).
* **UMD's k turns negative.** `0.34 − 4.2e-4·T` crosses zero at 809.5 K.
  ThermaKin takes the law as written; once a thin remnant passes 809.5 K the
  integration diverges ("integration diverged"). `thermakin.py` keeps such a
  run up to its last complete output and flags it (`diverged` in
  `meta.json`). UMD's own submissions end the same way: their last
  temperatures are 809.5–809.6 K (25 kW 6/12 mm, 65 kW 12 mm), i.e. the
  submitted files stop where ThermaKin diverged. (Pyrolysis.jl floors k at
  1 % of k(300 K).)
* **Gas.** The pyrolyzate is a G-state component carrying the deck's gas cp
  (UMD, UMET_TK: 2000 J/(kg K); otherwise 1 J/(kg K), negligible), transported
  with λ = 2e-5 m²/s and removed at the top by `J = 0.05·ξ_top` (the
  reference manual's rapid-removal value). At most 0.015–0.035 g of it is held
  in the sample (< 0.05 % of the initial mass); `mass_g` excludes it.

## Sets per code

All six cases (10/25/65 kW/m² × 6/12 mm, end times from `GASIFICATION_CASES`)
run for every set listed as run.

| Set | Gpyro V2 | ThermaKin2Ds |
|---|---|---|
| NIST | run (constants, surface) | run, `base` only (MAX = RAND for a surface absorber) |
| GIDAZE+ | run (native power laws, κ 2870 m⁻¹) | run, `base` + `rand` (power laws exact as `p2·T^n`) |
| Aalto_I | run (κ 3421.27 m⁻¹ as written) | run, `base` + `rand` |
| BUW-FZJ_A | run (residue via χ = 14.15) | run, `base` |
| BUW-FZJ_B | run (χ = 3.299; component 1 inert, as written) | run, `base` |
| UMET_GP | run (n = 3, 1.8, Gpyro's own normalisation; ρ held at 1206.7 as the deck) | not expressible: ThermaKin rates are first order (n = 3, 1.8) |
| UMET_TK | run; cp 7.5 + 3.7T replaced by 1117.21 (T/300)^0.995861, max error 0.032 % over 293–900 K | run, `base` + `rand` (all native) |
| UMD | skipped: best single power laws are off by 1.1 % (glass cp, 293–395 K), 0.9 % (melt cp, 395–900 K) and 16–34 % (melt k 0.34 − 4.2e-4 T, 395–700/750 K) | run, `base` + `rand` (glass transition as UMD's pseudo-reaction) |
| DBI_1 | skipped: cp table 1800 → 1400 J/(kg K) held above 413 K has no power law within 9.6 % (k table: 1.4 %) | not expressible (n = 1.03; tables) |

Aalto_I runs the deck's literal κ in both codes, as a one-way coefficient,
though Aalto fits it in FDS's two-flux model. Densities are constant in every
`base` run, as the decks hold them; `rhoT` (extra, below) restores UMET's ρ(T)
in UMET's own codes.

## Sanity: the institutions' own codes

`sanity.py --table` compares with the submissions in matl-db
(`PMMA/Computational_Results/Gasification_Predictions/`). Differences are
ours − theirs over the submission's span, burnout tail excluded (submission
mass > 1 % of initial); t50 = time to half the initial mass.

**UMD in ThermaKin: MAX at 10 kW/m², RAND at 25 and 65 kW/m².**

| Case | t50 UMD | t50 MAX | t50 RAND | mass rms MAX / RAND (g) | T_back rms MAX / RAND (K) |
|---|---|---|---|---|---|
| 10 kW 6 mm | 2322 | **2341** | 2124 | **0.24** / 4.2 | **2.1** / 5.2 |
| 10 kW 12 mm | 4868 | **4889** | 4437 | **0.30** / 9.4 | **1.6** / 9.1 |
| 25 kW 6 mm | 535 | 589 | **520** | 5.5 / **1.6** | 24.7 / **5.1** |
| 25 kW 12 mm | 1034 | 1153 | **1004** | 11.1 / **3.1** | 26.1 / **4.6** |
| 65 kW 6 mm | 159 | 169 | **155** | 4.0 / **1.2** | 47.5 / **4.6** |
| 65 kW 12 mm | 290 | 316 | **284** | 8.6 / **1.7** | 54.7 / **4.7** |

At 10 kW/m² MAX reproduces UMD's mass and back temperature to within the
initial-mass offset (UMD's samples are 0.17–0.33 % heavier: 72.84, 72.72 and
145.44 g, i.e. 6.02, 6.01 and 12.02 mm at 1210 kg/m³) and RAND burns 9 %
too fast; at 25 and 65 kW/m² RAND reproduces them (t50 2–3 % early) and MAX
is 6–12 % slow. The documents do not say which mode UMD used for MaCFP-2 (Fiola et al. 2021 is
silent; the group's MaCFP-3 UMD-FSRI README says "Random absorption
algorithm, a Monte-Carlo approach"); the numbers say MAX for the 10 kW/m²
runs and RAND for the others.

UMD's "Top Temperature (1mm)" column is neither the top element (rms 10–12 K
off in MAX at 10 kW, 25–56 K in RAND at 25–65 kW) nor a probe fixed 1 mm
below the initial surface that reads the top element once the surface has
passed it (the reference manual's visualizer convention; 9 K and 24–53 K).
It matches the mean temperature of the top 1 mm of sample — equivalently the
temperature about 0.5 mm below the moving surface (the interpolated depth that reproduces each value is 0.44–0.60 mm at
25 kW/m²): rms 1.3/2.0 K (MAX, 10 kW 6/12 mm), 4.7/2.8 K (RAND, 25 kW),
7.7/11.8 K (RAND, 65 kW). `aux.csv` carries these candidates
(`T_mean_top_1mm_K`, `T_0p5mm_K`, `T_fixed_height_L_minus_1mm_K`). It is
not the temperature 1 mm below the moving surface (`T_1mm_K`). UMD's files
end where their runs evidently stopped: the last temperatures are
809.5–809.6 K, the zero of their k law (above).

**GIDAZE+ in Gpyro: reproduced.** t50 within 1 s in all six cases (1786/1787,
3699/3701, 484/485, 924/924, 148/148, 267/267 s), mass rms 0.14–0.28 g, top
temperature rms 1.0–4.8 K. The differences sit in the last seconds of
burnout: GIDAZE+'s remnant heats to ~790 K while V2 loses the transmitted
radiation (above). Their "Back surface temperature" for the 12 mm cases is
the temperature of the node initially 6 mm deep (matches our profile there to
~1 K; the back face itself is up to 41 K colder at 10 kW, 160 K at 25 kW and
350 K at 65 kW) — a `POINT_Z = 0.006` left over from the 6 mm input.

**UMET_GP in Gpyro: reproduced once ρ(T) is restored.** With the deck's
constant density (`base`): t50 within 0.6–2.3 %, mass rms 0.4–1.6 g,
T_back rms 2–14 K. With UMET's power-law density 1200 (T/300)^−0.236
(`rhoT`, extra variant): t50 within 0.3–1.1 %, temperatures within 0.2–0.8 K rms
(≤ 2.9 K max), mass rms 0.3–0.9 g, mostly the initial-mass offset: UMET's
initial masses (72.25/144.50 g) are those of 1204.2 kg/m³ (UMET_TK's density
at 293 K), 0.2 % below the
power law's 1206.7.

**UMET_TK in ThermaKin: not reproduced.** Neither mode nor ρ(T) matches the
submission within a few per cent: t50 ours/theirs MAX 2255/2374, 4565/4816,
555/539, 1122/1088, 178/171, 335/273 s; RAND 2145, 4344, 525, 1060, 171, 321 s;
MAX with ρ = 1380 − 0.6T (`rhoT`) 2304, 4691, 563, 1136, 178, 334 s. Their
65 kW 12 mm case burns much faster than any variant. The unknowns are UMET's
absorption mode, κ unit (2700 m⁻¹ or m²/kg — the latter is a surface
absorber, which MAX already is), gas cp and swelling settings.

## Convergence

Base resolution: Gpyro Δz = 5e-5 m (121/241 nodes), dt0 = 0.02 s; ThermaKin
Δx = 5e-5 m, dt = 0.01 s. Each `_half` run halves both. Differences are base
− half over the run (sample mass > 2 % of initial; max / rms):

| Code, set, case | mass (g) | T_top (K) | T_back (K) | T_1mm (K) | t50 base/half (s) | peak MLR |
|---|---|---|---|---|---|---|
| Gpyro NIST 25 kW 6 mm | 0.03 / 0.02 | 0.36 / 0.04 | 0.44 / 0.06 | 0.44 / 0.06 | 531/531 | +0.0 % |
| Gpyro GIDAZE+ 25/6 | 0.03 / 0.01 | 0.06 / 0.02 | 0.13 / 0.05 | 0.13 / 0.04 | 484/485 | +0.0 % |
| Gpyro Aalto_I 25/6 | 0.01 / 0.01 | 0.02 / 0.01 | 0.04 / 0.02 | 0.04 / 0.02 | 465/465 | +0.0 % |
| Gpyro BUW-FZJ_A 25/6 | 0.04 / 0.02 | 0.54 / 0.06 | 0.60 / 0.09 | 0.60 / 0.08 | 531/531 | +0.0 % |
| Gpyro BUW-FZJ_B 25/6 | 0.05 / 0.02 | 4.5 / 0.24 | 4.9 / 0.26 | 4.9 / 0.26 | 517/517 | +0.1 % |
| Gpyro UMET_GP 25/6 | 0.01 / 0.00 | 0.03 / 0.01 | 0.04 / 0.02 | 0.04 / 0.01 | 534/534 | +0.0 % |
| Gpyro UMET_TK 25/6 | 0.01 / 0.00 | 0.02 / 0.00 | 0.02 / 0.01 | 0.02 / 0.01 | 529/529 | +0.0 % |
| Gpyro NIST 10/12 | 0.01 / 0.00 | 0.04 / 0.00 | 0.02 / 0.00 | 0.02 / 0.00 | 4336/4336 | +0.0 % |
| Gpyro GIDAZE+ 10/12 | 0.01 / 0.01 | 0.01 / 0.00 | 0.02 / 0.01 | 0.02 / 0.01 | 3699/3700 | +0.0 % |
| ThermaKin NIST MAX 25/6 | 0.15 / 0.08 | 1.6 / 0.74 | 2.1 / 0.29 | 2.1 / 0.30 | 525/526 | +0.2 % |
| ThermaKin UMD MAX 10/6 | 0.44 / 0.28 | 3.0 / 0.73 | 3.2 / 0.81 | 3.2 / 0.80 | 2341/2354 | +0.7 % |
| ThermaKin UMD MAX 25/6 | 1.66 / 0.89 | 8.1 / 3.6 | 15.1 / 3.3 | 15.1 / 3.1 | 589/598 | +2.9 % |
| ThermaKin GIDAZE+ MAX 25/6 | 1.19 / 0.62 | 6.8 / 2.8 | 12.4 / 2.7 | 12.4 / 2.6 | 519/524 | +1.8 % |
| ThermaKin UMD RAND 25/6 | 0.36 / 0.19 | 6.2 / 1.6 | 3.2 / 1.0 | 3.2 / 1.1 | 520/518 | +0.7 % |
| ThermaKin GIDAZE+ RAND 25/6 | 0.50 / 0.24 | 7.2 / 1.5 | 3.3 / 1.3 | 3.3 / 1.3 | 458/456 | −0.4 % |

Gpyro is converged to plotting accuracy everywhere; its largest differences
are in the last seconds of burnout (e.g. NIST 6 K in the final 2 s) and at
the instant BUW-FZJ_B's last reacting component is used up. ThermaKin is
converged for surface absorption and for RAND; MAX on an in-depth deck
converges more slowly (t50 1–1.5 %, peak MLR 2–3 % between the two
resolutions), because MAX deposits the net flux in the top element and
emits at its centre temperature, both of which depend on Δx (the limit is a
surface absorber). The MAX-vs-RAND reading of UMD's submission holds at half
size (10 kW 6 mm MAX t50 2354 vs UMD 2322; 25 kW 6 mm RAND 518 vs 535).

## RAND seed spread

UMD, 25 kW/m², 6 mm, seeds 1, 2, 3 (`rand`, `rand_seed2`, `rand_seed3`); and
65 kW/m², 6 mm, seeds 1, 2. Spread (max − min across seeds) while the sample
mass is above 1 % of its initial, max / rms:

| Case | mass (g) | MLR (g/(m² s)) | T_top (K) | T_back (K) | T_1mm (K) | t50 (s) | peak MLR |
|---|---|---|---|---|---|---|---|
| 25 kW 6 mm, 3 seeds | 0.08 / 0.03 | 0.71 / 0.20 | 7.4 / 2.4 | 6.1 / 0.85 | 6.1 / 1.1 | 520, 520, 520 | 14.74, 14.74, 14.67 |
| 65 kW 6 mm, 2 seeds | 0.27 / 0.08 | 4.2 / 1.3 | 15.9 / 4.8 | 8.2 / 1.3 | 8.2 / 2.3 | 155, 155 | 48.6, 48.3 |

The draw is per 0.01 s step, so the stochastic part averages out over the
1 s outputs; the spread is far below the MAX–RAND difference. After burnout
the seeds differ more (the thin remnant is heated only in steps whose draw
lands in it), and the 65 kW runs diverge at different times (274, 300 s).

## Runs and wall times

Production: Gpyro 42 runs (7 sets × 6 cases), ThermaKin 66 (7 sets × 6 cases
MAX + 4 in-depth sets × 6 cases RAND). Extras: Gpyro 9 `base_half` and 6
`rhoT`; ThermaKin 6 `_half`, 3 seed and 6 `rhoT` runs. No run fails outright.
Wall times (one thread each, measured with three or four runs and unrelated
Julia jobs sharing a 22-core machine; `meta.json` per run): Gpyro base
1234 s for 42 runs (0.5–244 s; longest BUW-FZJ_A 10 kW 12 mm), ThermaKin MAX
2101 s for 42 runs (2–511 s; longest UMD 10 kW 12 mm), RAND 1415 s for 24
runs (2–591 s); extras 1.5 ks in total.

Flags (`sanity.py --runs`):

* ThermaKin, "integration diverged" after burnout (UMD's k < 0 above
  809.5 K), kept to the last output: UMD `rand` 25 kW 12 mm (at 1686 s),
  65 kW 6 mm (300 s); `rand_seed2` 65 kW 6 mm (274 s); `rand_half` 25 kW 6 mm
  (920 s).
* Gpyro, stopped at burnout (series end there): every NIST, GIDAZE+, Aalto_I
  and UMET_TK case, by V2's zero-mass shutdown with < 0.7 % of the mass left;
  UMET_GP at 25 and 65 kW (thickness criterion at 25 kW 12 mm, zero-mass
  shutdown otherwise). BUW-FZJ_A/B (residue) and UMET_GP at 10 kW (slow
  n = 1.8 tail) run to the end time.
* ThermaKin MAX on in-depth decks leaves 0–0.39 g unburnt at the end time
  (Aalto_I 0.30–0.35, GIDAZE+ 0–0.22, UMET_TK 0.32–0.39, UMD 0.21–0.35 g
  including 0.14/0.29 g of char), RAND at 10 kW/m² 0.11–0.12 g; both from
  the backing and the transmitted radiation (above).

## Approximations and their size

Everything else is the deck as written.

* Both codes: densities constant (decks); transmitted radiation leaves through
  the back (Pyrolysis.jl's mirror returns it into the sample) — matters only
  in the last ~0.5–1 mm of in-depth decks. σ: Gpyro 5.67e-8, ThermaKin 5.6704e-8
  W/(m² K⁴) against 5.670374e-8 (−0.007 %, +0.0005 %); E/R preserved exactly.
* Gpyro: UMET_TK cp power-law fit (0.032 %); UMET's intermediate B 1e-4
  lighter than A (B's volume +0.01 %); in-depth decks deposit the net flux in
  depth (V2's scheme, not a choice); residual mass `EPS` = 1e-7 per consumed
  cell.
* ThermaKin: the backing (0.17 %/0.08 % extra heat capacity at the back face
  of a 6/12 mm sample; the MAX burnout tail above); `T_top` half an element
  (≤ 25 µm) below the surface; gas cp 1 J/(kg K) where the deck has none
  (≤ 0.035 g of gas in the sample); UMD's glass transition as a fast
  pseudo-reaction (A = 100 s⁻¹, so ~0.01 s lag at 395 K; the glassy polymer
  does not react below 395 K, < 0.5 % conversion of the deck's reaction 1);
  UMD's negative k above 809.5 K left as written (divergence after burnout);
  MAX on in-depth decks is ThermaKin's own surface-like deposition, RAND its
  sampled Beer–Lambert deposition — both are native modes, neither is Gpyro's
  continuous Beer–Lambert source or Pyrolysis.jl's two-stream transfer.

## Files

`scripts/macfp_pmma/external/`: `common.py`, `gpyro.py`, `thermakin.py`,
`sanity.py`, `build_thermakin_seeded.sh`, this README. Results:
`scripts/macfp_pmma/results/{gpyro,thermakin}/` (ignored), with
`runs.jsonl` (append-only run log) per code. ThermaKin source and executables
stay in the ignored `docs/refs/ThermaKin2Ds Code/`, `docs/refs/thermakin_build/`
and `docs/refs/thermakin_build_seeded/`.
