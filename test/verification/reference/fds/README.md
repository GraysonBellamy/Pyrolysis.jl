# FDS 6.11.1 verification outputs

These device outputs and the in-depth profile come from the official FDS
6.11.1 Linux release bundle. The twenty supplied verification decks are the
unmodified files under its `Examples/Heat_Transfer/` and `Examples/Pyrolysis/`
directories. The eight matched decks in `custom/` use the same executable.
Every run completed successfully. Input and output SHA-256 hashes, executable
identity and timing samples are recorded in `provenance.json`.

- Executable: `FDS-6.11.1-0-gff928db-release`, Intel Fortran 2025.3.2.
- [Official release](https://github.com/firemodels/fds/releases/tag/FDS-6.11.1):
  `FDS-6.11.1_SMV-6.11.2_lnx.sh`.
- Bundle SHA-256: `ba8793b974150fdb778b3db0b69ab8db4e4668d5c52987a317e7f5ed58102ea0`.
- Executable SHA-256: `fe2b123dfc759dbebb161cc91912ed60861f3fe2a141602634b53f417f731b1d`.
- One MPI rank and one OpenMP thread, with `I_MPI_FABRICS=shm`.
- Units are FDS-native (first header row): time s, temperature **°C**, mass
  kg or kg/m², thickness m. Convert T with +273.15 when overlaying.

| File(s) | FDS case | Overlays our case |
|---|---|---|
| `heat_conduction_{a..d}_devc.csv` | finite slab, convective BC; FDS a/b/c/d = Bi 100/10/1/0.1 | V-HC-02 d/c/b/a (labels reversed; nondimensional θ–Fo overlay) |
| `convective_cooling_devc.csv` | semi-infinite-like convective heating | V-HC-01 |
| `insulated_steel_plate_prof_1.csv` | constant flux, insulated back | no matched FDS overlay; see the boundary-condition limitation below |
| `two_step_solid_reaction_devc.csv` | A→B→C series kinetics | V-RX-01 |
| `part_baking_soda_{420K,450K,500K}_devc.csv` | isothermal NaHCO₃ decomposition | V-RX-02 |
| `shrink_swell_devc.csv` | 6 shrink/swell variants, M_i/L_i/T_i | V-CPL-03 |
| `matl_e_cons_{1..9}_devc.csv` | enthalpy conservation ledger | reference outputs; no report overlay |

Run each deck in its own output directory with:

```sh
LD_LIBRARY_PATH=<bundle>/bin/intelmpi/lib \
FI_PROVIDER_PATH=<bundle>/bin/intelmpi/prov \
I_MPI_FABRICS=shm OMP_NUM_THREADS=1 <bundle>/bin/fds <case>.fds
```

The FDS rows in `../timings.csv` are medians of five fresh-process runs of
each matched deck after one warm-up, pinned to one CPU. The primary `wall_s`
value is FDS's **Total Elapsed Wall Clock Time**: input reading, setup,
simulation and output after MPI initialization, ending at the final timing
report before MPI finalization. It is not the time-stepping-only timer.
`process_wall_s` retains the full process measurement, including launch and
shutdown, for the supporting table. All samples of both timers and the
time-stepping timer remain in `provenance.json`.

Timer boundaries are defined in FDS's
[main routine](https://github.com/firemodels/fds/blob/FDS-6.11.1/Source/main.f90)
and [timing report](https://github.com/firemodels/fds/blob/FDS-6.11.1/Source/dump.f90).
These are indicative cross-code costs: Julia times warmed model construction
and solve, and the codes retain different output and time-stepping workloads.

## `custom/` — our decks, run with the same binary

Where the shipped FDS case differs from ours, `custom/` holds minimal decks
(committed next to their outputs) that make FDS solve *our* case definition:

- `vhc01_semi_infinite.fds` — V-HC-01 exactly, via the same Robin-BC
  equivalence the case itself uses (absorbed 25 kW/m² + h = 20 W/m²K toward
  300 K ≡ pure convection toward T∞_eff = 1550 K; `EMISSIVITY=0` so radiation
  is inert — note FDS multiplies `EXTERNAL_FLUX` by ε, so the flux form can't
  be used with ε = 0). Front face matches the erfc solution to ≲0.1 K.
- `vcpl03_shrink_swell_paths.fds` — V-CPL-03 P1/P2 thickness *paths*:
  the shipped `shrink_swell.fds` heats with an Arrhenius transient (only the
  endpoints are comparable), so this deck runs the isothermal first-order
  conversion (A = 2 s⁻¹, E = 0, ΔH = 0) with our densities (P2: 1000→2000,
  P1: 1000→500) and L₀ = 2 mm. Endpoints exact; devc values are averages over
  the 0.1 s output interval (plot at t − 0.05 s).
- `heat_conduction_a_ext.fds` — shipped case a (Bi = 100) with `T_END`
  extended 2000 → 5000 s so the Fo = 0.5 profile exists for V-HC-02d.
- `heat_conduction_{b,c,d}_matched.fds`, `two_step_solid_reaction_matched.fds`,
  `part_baking_soda_450K_matched.fds` — shipped cases with `T_END` set to our
  case's simulated span (5000 s ⇒ Fo = 0.5; 20 s; 2.4 s ≈ 0.9·t_burnout).
  The shipped originals simulate up to 10× longer (mostly the featureless
  asymptote), which would inflate their wall times in the solve-time
  comparison and stretch the overlay axes; all V-HC-02/V-RX-01/V-RX-02c
  overlays and `timings.csv` rows therefore use these matched runs. The
  supplied-deck outputs are retained alongside the matched-deck outputs.

Not representable in FDS: V-HC-03's pure Neumann BC (constant net flux, no
convection, no re-radiation). `NET_HEAT_FLUX`/`CONVECTIVE_HEAT_FLUX` don't
work with `MATL_ID`, and `EXTERNAL_FLUX` (even with `SKIP_INRAD=T`) always
carries the ε·σ·Ts⁴ re-emission term. V-HC-03 gets its cross-code overlay
from Gpyro instead.
