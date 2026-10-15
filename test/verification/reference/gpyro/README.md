# Gpyro reference-code runs for the Pyrolysis.jl verification suite

Date: 2026-07-07

## Code

- **Gpyro 0.8200** (banner of every `.out` file: `Gpyro 0.8200`,
  `https://github.com/reaxfire/gpyro`). Source checkout: commit
  `6d369dd0fa5a3621334272d0806f48c56e608c01` (2020-09-01), the current master
  of the maintained Gpyro repository (reaxfire is the maintained home of
  lautenberger/gpyro; identical v0.8200 source).
- Serial standalone build (`build/linux/gpyro/gpyro`), compiled 2026-07-07
  with **GNU Fortran 15.2.0 (conda-forge gcc 15.2.0-19)**.
- Smoke-tested against the shipped `samples/01-dtg` case before use.
- Extraction/verification script: `extract_results.py` (same directory);
  it converts Gpyro's output (temperatures in deg C) to SI (K, m, s) and
  cross-checks every run against the analytical solution.

## Conventions common to all decks

- Fixed-format namelist decks derived from the shipped samples
  (`samples/01-dtg`, `samples/02-thermoplastic`).
- Gpyro's 1D grid places node 1 exactly at the exposed surface z = 0 and the
  last node at z = L, so profile dumps need no half-cell shift; `x_m` in the
  result CSVs is depth below the irradiated/cold face.
- Output CSVs written by Gpyro are in deg C; `*_result.csv` files here are
  converted to K. All other quantities already SI.
- Gas-phase solvers disabled (`SOLVE_GAS_YJ/ENERGY/PRESSURE = .FALSE.` for the
  1D cases), `EMIS = 1` with `RERADIATION = .FALSE.` so the absorbed flux is
  exactly the prescribed `QE`, no melting/reactions unless the case needs them.

## Cases

### vhc01 — V-HC-01 (dissertation §3.4.1, Fig. 3.9)

Semi-infinite inert solid, absorbed radiant flux + convective loss
(erfc exact solution).

| Parameter | Value | Same as V-HC-01.jl? |
|---|---|---|
| k | 0.2 W/m-K | yes |
| rho | 1000 kg/m3 | yes |
| c | 1400 J/kg-K | yes |
| q_e (absorbed) | 25 kW/m2 | yes (EMIS=1, KAPPA=9e9 -> surface absorption) |
| h_c | 20 W/m2-K, T_inf = T0 = 300 K | yes |
| L | 0.05 m, adiabatic back | yes (dissertation used 0.04 m) |
| Grid / dt | dz = 0.1 mm (500 nodes), dt = 0.1 s | dissertation numerics |

Deltas from the Julia case: **none** in physics.
Result files: `vhc01_result.csv` (t_s, x_m, T_K profiles at t = 30/60/90/180 s),
`vhc01_surfaceT.csv` (t_s, T_surface_K, T_back_K at 1 s intervals).
Measured agreement with the exact erfc solution: L-inf 0.138 / 0.091 / 0.071 /
0.045 K at t = 30/60/90/180 s; back face undisturbed (300.00 K).

### vrad01 — V-RAD-01 (dissertation §3.4.2 + Appendix A, Fig. 3.10)

Finite slab with in-depth radiation absorption, convective front loss,
insulated back (eigenfunction exact solution; sub-surface temperature peak).

| Parameter | Value | Same as V-RAD-01.jl? |
|---|---|---|
| k, rho, c | 0.2 / 1000 / 1400 | yes |
| delta | 0.01 m | yes |
| kappa | 1000 1/m (`KAPPA(1)=1000` activates Gpyro's in-depth Beer-Lambert sweep) | yes (Julia: mass coefficient alpha = kappa/rho = 1.0 m2/kg) |
| q_e | 25 kW/m2, EMIS = 1, no reradiation | yes |
| h_c | 20 W/m2-K at irradiated face, T_inf = T0 = 300 K | yes |
| Grid / dt | dz = 0.1 mm (101 nodes), dt = 0.1 s | dissertation numerics |

Deltas from the Julia case: **none**. Gpyro transmits the unabsorbed remainder
exp(-kappa*delta)*q_e out of the slab, matching the exact model and the Julia
solver.
Result files: `vrad01_result.csv` (profiles at 30/60/90/180 s),
`vrad01_surfaceT.csv`. Measured agreement with the eigenfunction solution:
L-inf 0.091 / 0.074 / 0.063 / 0.040 K; sub-surface peak present and migrating
(0.2 mm at 30 s -> 0.4 mm at 180 s), exactly the coupled-physics marker
V-RAD-01 checks.

### vcpl01 — V-CPL-01 (dissertation §3.4.3, Eqs. 3.137-3.138, Fig. 3.11)

Neumann two-phase Stefan (solidification) problem via apparent heat capacity:
Gaussian melt peak `TMELT=433 K, DHMELT=50 kJ/kg, SIGMA2MELT=0.1 K^2`
(Gpyro Eq. 3.5c; the deck's latent heat enters through Gpyro's analytic erf
melt-enthalpy term).

| Parameter | Value | Same as V-CPL-01.jl (route A)? |
|---|---|---|
| k_s = k_l | 0.2 W/m-K | yes |
| rho | 900 kg/m3 | yes |
| c_s = c_l | 1700 J/kg-K | yes |
| T_m / DH_m / sigma_m^2 | 433 K / 50 kJ/kg / 0.1 K^2 | yes |
| T0 / T_wall | 438 K / 428 K (TFIXED pins the z=0 surface node) | yes |
| L | 0.04 m, adiabatic back | yes |
| Grid / dt | uniform dz = 0.05 mm (801 nodes), dt0 = 0.02 s | **delta**: Julia case uses a stretched mesh (~14 um at wall); Gpyro uses the dissertation's uniform 0.05 mm. dt: Gpyro's own step controller cut dt to ~0.01 s after early convergence failures at the cp spike (logged; run completes cleanly). |

Both the dissertation and V-CPL-01.jl use this exact parameter set, so a
single deck serves both.
Result files: `vcpl01_result.csv` (t_s, x_m, T_K histories at the four
dissertation depths x = 0.5/1.0/1.5/2.0 mm, ~1 s sampling — the Fig. 3.11
quantities), `vcpl01_front.csv` (t_s, s_front_m: T_m-isotherm crossing).
Eigenvalue check: lambda = 0.2317 (matches dissertation). Front position vs
s(t) = 2 lambda sqrt(alpha_s t): -0.87% / -0.00% / +0.25% at t = 100/200/299 s;
depth-history errors <= 0.17 K. Back face undisturbed (437.9999 K).
Note: `PROFILE_ISKIP=2` (profile dumped every other node, 0.1 mm) because
Gpyro v0.8200 truncates profile CSV rows at 10000 characters (801 nodes would
overflow and corrupt the row); the solve grid itself is unchanged and all four
QoI depths lie exactly on dumped nodes.

### vhc03 — V-HC-03 (not a dissertation case; C&J §3.8 finite-slab Neumann)

Finite slab, constant net surface flux, adiabatic back — the pure Neumann BC
that FDS's 1D solid cannot express without eps-scaling or sigma*T^4
re-emission. Gpyro can: `EMIS=1`, `RERADIATION=.FALSE.`, `HC=0`, surface
absorption via `KAPPA=9D9` (> 1e6 deposits the full `QE` in the surface node;
the in-depth kappa route was not needed).

| Parameter | Value | Same as V-HC-03.jl? |
|---|---|---|
| k / rho / c | 0.5 W/m-K / 1000 kg/m3 / 1000 J/kg-K (alpha = 5e-7 m2/s) | yes |
| q'' (absorbed) | 5000 W/m2 at top face, h = 0, no reradiation | yes |
| L | 0.02 m, adiabatic back, T0 = 300 K | yes |
| Grid / dt | dz = 0.1 mm (201 nodes), dt = 0.1 s | — |
| TSTOP | 1202 s (not 1200) so the t = 1200 s profile dump fires before exit | output timing only |

Deltas from the Julia case: **none** in physics.
Result files: `vhc03_result.csv` (t_s, x_m, T_K profiles at t = 80/400/1200 s,
Fo = 0.1/0.5/1.5), `vhc03_frontback.csv` (t_s, T_front_K, T_back_K at 2 s
intervals). Measured agreement with the exact Fourier-series solution:
L-inf 0.012 / 0.001 / 0.000 K at Fo = 0.1/0.5/1.5; mean-T anchor exact to
< 1 mK at all three times (q/(rho c L) = 0.25 K/s: 320/400/600 K); front-back
dT = 100.000 K at Fo = 1.5 (plateau q''L/2k = 100 K), 69.769 K vs exact
69.788 K at Fo = 0.1.

### vrx01 — V-RX-01 (dissertation §3.4.4, Eqs. 3.139-3.140, Fig. 3.12)

Two-step series chain A -> B -> C, isothermal (Bateman exact solution);
0D (ZEROD) run, the analog of the FDS `two_step_solid_reaction` benchmark.

| Parameter | Value | Same as V-RX-01.jl? |
|---|---|---|
| k_AB | 0.389 1/s (Z=0.389, E=0 -> rate == Z at any T) | yes |
| k_BC | 0.262 1/s (Z=0.262, E=0) | yes |
| Hold temperature | 400 K (TAMB=400, BETA=0 -> particle held at exactly 400 K) | yes (inert with E=0) |
| Densities | R0 = 1000 kg/m3 for A, B, C (the Julia xi0) -> no gas evolved, unit condensed yields, YI == Bateman mass fractions | yes |
| Reaction heats | DHS = DHV = 0 | yes (h = 0) |
| dt | 0.001 s | **delta**: finer than the dissertation-era 0.1 s; Gpyro's species update is first-order in time (measured error 6e-4 in Y_A at dt = 0.01, 7e-5 at dt = 0.001) |

Result file: `vrx01_result.csv` (t_s, Y_A, Y_B, Y_C at 0.1 s intervals).
Measured agreement with Bateman: max |dY| = 7.2e-5 / 6.1e-5 / 1.1e-4 (A/B/C);
M/M0 = 1 to machine output precision; T = 400.0000 K throughout.

## Correspondence to dissertation figures

| Deck | Julia case | Dissertation figure |
|---|---|---|
| vhc01 | V-HC-01 | Fig. 3.9 |
| vrad01 | V-RAD-01 | Fig. 3.10 |
| vcpl01 | V-CPL-01 | Fig. 3.11 |
| vrx01 | V-RX-01 | Fig. 3.12 |
| vhc03 | V-HC-03 | — (not in the dissertation; Carslaw & Jaeger §3.8 exact series) |

## Run logs

`<case>_run.log` (Gpyro stdout) each end with
"End of simulation reached successfully"; the vcpl01 log additionally shows
the 20 early adaptive-dt reductions at the apparent-cp spike (expected).

### vbo01 — V-BO-01 (MaCFP PMMA gasification, 50 kW/m²)

`vbo01_gpyro_v2_pmma.csv` is the **Gpyro V2** (SoftwareX release, GPL-3
GPYRO-V2 repository) run made here on the authors' MaCFP PMMA deck
(`docs/gpyro_v2_comparison/reference/pmma/PMMA.data`: 5.9 mm PMMA on 12 mm
Kaowool, 50 kW/m², power-law k(T) and c(T) fits), copied verbatim from
`docs/gpyro_v2_comparison/reference/pmma/gpyro_v2_pmma_localrun.csv`. Gpyro V2
was built with the docker recipe in `docs/gpyro_v2_comparison/README.md`. The
first header row names the columns; column 2 (`001_MLR`) is the mass loss rate
in **g m⁻² s⁻¹**, column 6 the remaining thickness in m, columns 3–5 the
temperatures (°C) at the surface, 12 mm and 17.9 mm. Feeds the V-BO-01 MLR
overlay in report mode.
