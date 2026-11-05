# FDS 6.11.1 verification outputs

Every FDS comparison on the verification pages is a matched native-input run
documented in [`../matched/README.md`](../matched/README.md), with its
refinements and cost. The files here are outputs of FDS's own supplied
verification decks and of two custom decks; no verification page reads them.

These device outputs and the in-depth profile come from the official FDS
6.11.1 Linux release bundle. The twenty supplied verification decks are the
unmodified files under its `Examples/Heat_Transfer/` and `Examples/Pyrolysis/`
directories. The two decks in `custom/` use the same executable.
Every run completed successfully. Input and output SHA-256 hashes and the
executable identity are recorded in `provenance.json`.

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
| `heat_conduction_{a..d}_devc.csv` | finite slab, convective BC; FDS a/b/c/d = Bi 100/10/1/0.1 | Same nondimensional problem as V-HC-02 d/c/b/a; not read by an accuracy overlay |
| `convective_cooling_devc.csv` | semi-infinite-like convective heating | Not read by an accuracy overlay |
| `insulated_steel_plate_devc/prof_1.csv` | heated steel plate | retained; not the V-HC-03 overlay |
| `two_step_solid_reaction_devc.csv` | A→B→C series kinetics | Same chain as V-RX-01; not read by an accuracy overlay |
| `part_baking_soda_{420K,450K,500K}_devc.csv` | isothermal NaHCO₃ decomposition | Not read by an accuracy overlay |
| `shrink_swell_devc.csv` | 6 shrink/swell variants, M_i/L_i/T_i | Endpoint check; comparison paths are in `../matched/` |
| `matl_e_cons_{1..9}_devc.csv` | enthalpy conservation ledger | retained; not read by a current overlay |

Run each deck in its own output directory with:

```sh
LD_LIBRARY_PATH=<bundle>/bin/intelmpi/lib \
FI_PROVIDER_PATH=<bundle>/bin/intelmpi/prov \
I_MPI_FABRICS=shm OMP_NUM_THREADS=1 <bundle>/bin/fds <case>.fds
```

## `custom/` — our decks, run with the same binary

- `vcpl03_shrink_swell_paths.fds` — V-CPL-03 P1/P2 thickness *paths*:
  the shipped `shrink_swell.fds` heats with an Arrhenius transient (only the
  endpoints are comparable), so this deck runs the isothermal first-order
  conversion (A = 2 s⁻¹, E = 0, ΔH = 0) with our densities (P2: 1000→2000,
  P1: 1000→500) and L₀ = 2 mm. Endpoints exact; devc values are averages over
  the 0.1 s output interval (plot at t − 0.05 s).
- `part_baking_soda_450K_matched.fds` — the shipped case with `T_END` set to
  2.4 s ≈ 0.9·t_burnout. V-RX-02c uses the direct native-input run in
  `../matched/` instead.

V-HC-03 has no FDS comparison: its face takes a constant net flux with no
convection or reradiation. In FDS 6.11.0's standard `MATL_ID` wall input path,
[`READ_SURF`](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90#L8290-L8300)
rejects `NET_HEAT_FLUX` and `CONVECTIVE_HEAT_FLUX`; `EXTERNAL_FLUX` is
absorptivity-scaled and the material reradiates. This identifies the limitation
of the tested input path, not every possible FDS configuration or release.
