# Matched comparisons with FDS, Gpyro and ThermaKin2Ds

These accuracy comparisons use native inputs to unmodified FDS 6.11.1 and
Gpyro V2, and to ThermaKin2Ds 7 built from its source with a banner-only
portability change. There are 73 native decks (21 FDS, 28 Gpyro, 24
ThermaKin), each run at seven configurations (504 successful runs, 7 recorded
unavailable), spanning 28 verification case IDs.
V-HC-04c and V-HC-04d share the same physical problem and explicitly reuse
the same native runs: only Julia's mesh differs. P1/P2/P3 are counted
separately within V-CPL-03. V-GT-08 is the only pressure-driven case. The
solve times of these runs and of Pyrolysis.jl on the same curves are in
[`timing/`](timing/), measured under `scripts/timing/README.md`.
See the [coverage matrix](../../../../docs/guide/verification/comparison_coverage.md)
for partial comparisons and exclusions.

## Reproduce and audit

Python 3.12+ and the standard library are sufficient:

```sh
python3 test/verification/reference/matched/run_comparisons.py fds \
  /path/to/fds /tmp/verification-comparisons
python3 test/verification/reference/matched/run_comparisons.py gpyro \
  /path/to/gpyro /tmp/verification-comparisons
test/verification/reference/matched/build_thermakin.sh
python3 test/verification/reference/matched/run_comparisons.py thermakin \
  docs/refs/thermakin_build/thermakin /tmp/verification-comparisons
```

The ThermaKin2Ds source is not redistributable and is not committed. Place
the authors' C++ source in `docs/refs/ThermaKin2Ds Code/` (ignored);
`build_thermakin.sh` copies it to the ignored `docs/refs/thermakin_build/`,
replaces the MSVC-only `localtime_s`/`asctime_s` start-time banner calls
with POSIX `localtime_r`/`asctime_r`, and compiles with `g++ -std=c++11 -O2`.
`thermakin_provenance.json` records the source and executable hashes and
the compiler.

`--extract-only` validates an existing workspace without running the solver.
The default levels are `base` (h,dt), `refined` (h/2,dt/2), `time_half`
(h,dt/2), `time_quarter` (h,dt/4), `space_half` (h/2,dt/4), and the coarser
steps `time_double` (h,2dt) and `time_quadruple` (h,4dt), which look for a
cheaper configuration of each code at its accuracy. A coarse configuration a
code cannot express or complete is recorded under `unavailable` in its
provenance with the reason: ThermaKin's output interval must be a whole number
of steps (V-HC-01, V-RX-04a and V-GT-01 at 4dt); Gpyro's V-RX-04a at 4dt
stops one step past the requested end and its V-CPL-01 at 4dt skips output
times; FDS's accumulated time at V-HC-01's 2dt and 4dt misses the 30 s
profile and writes it one step later. The coupled refinement is not used to
infer isolated spatial or temporal order.
Use `--cases=vhc04a --levels=base,time_half,time_quarter,space_half
--destination=/tmp/comparison-pilot` for a separate pilot; partial studies
cannot overwrite this directory's complete manifests.
The script verifies generated inputs, successful termination (not merely an
exit status), pinned version banners, finite numeric outputs, monotone times
and the final output time. It writes base samples in `fds/` and `gpyro/`, finer
samples in `refined/`, and native inputs/output CSVs in `raw/`. Native logs
containing elapsed times are not published; the timing campaign's summary is. Each run retains
a hashed `grid.json` containing actual initial coordinates, indexing convention,
interval count and spatial unknown count. FDS coordinates come from its native
mesh table, not probe count; Gpyro coordinates come from its native profile,
except V-RX-01, whose `ZEROD` run has one lumped unknown.
The TOML manifests
bind each series to an exact case, curve and named analytical reference.
Report generation fails on a missing binding, missing samples or nonfinite
analytical errors. The extracted coordinates are never fitted or shifted to
improve agreement.

The provenance JSON files hash the executable, manifests, retained decks and
outputs. FDS is the official Linux release
`FDS-6.11.1-0-gff928db-release`, compiled with Intel Fortran 2025.3.2. Gpyro is
[revision bff0126c68f2a35c084ef04490ae327185f281a3](https://github.com/Gedeben/GPYRO-V2/tree/bff0126c68f2a35c084ef04490ae327185f281a3),
built with GNU Fortran 14.4.0 and the upstream `gnu_linux` flags
(`-O2 -unroll -frecord-marker=4 -ffree-line-length-none` with the upstream
OpenMP flags) and no solver-source changes. The source archive SHA-256 is
`061e033d05501e65874738641b441d078245e00f3ef5b30d3f8fe19ac7a3648b`. From the
extracted source's `build/linux/gpyro` directory:

```sh
GPYRO_FCOMPL_GNU=gfortran make -f ../../Makefile_gpyro gnu_linux \
  GIT_HASH=bff0126c68f2a35c084ef04490ae327185f281a3 \
  GIT_STAT= GIT_BRANCH=archive GIT_DATE=unspecified
```

Builds from the same archive with this recipe differ only by their embedded
build date and reproduce every run's outputs byte for byte. Gpyro can exit
with status zero after a failed solve, so the script requires its
successful-termination message rather than the exit status.
Each run requests one OpenMP thread and FDS uses one MPI rank. Gpyro's banner
reports OpenMP disabled in this build.

Regenerate accuracy reports, including per-case `comparison_accuracy.csv`
with maximum and RMS error against the shared analytical functions:

```sh
julia --project=scripts/verification scripts/verification/run_reports.jl \
  --cases=V-HC-01,V-HC-02,V-HC-03,V-HC-04,V-HC-05,V-RX-01,V-RX-02,V-RX-04,V-RX-05,V-CPL-01,V-CPL-03,V-CPL-04,V-GT-01,V-GT-07,V-GT-08 \
  --timings=false --comparisons-only=true --outdir=/tmp/verification-comparison-report
```

Only base external curves are plotted. Every configuration enters the
comparison-accuracy CSV and the cost comparison, which pairs each external
run with each code's cheapest configuration at least as accurate. The command
above produces partial chapters without wall times. Omit
`--comparisons-only=true` to render every curve, including unrelated sweeps.

## Native-input mappings

Coordinates map from native depth to `z = L - depth`. Temperature outputs in
degrees Celsius gain 273.15 K. FDS devices use `TIME_AVERAGED=F`; native output
times are preserved. Gpyro point-output time precision is 0.0001 s, and its
thickness precision is 1e-7 m. Initial boundary values are omitted from the
convective-slab histories, where the eigenfunction expansion is singular at
the initial corner. FDS can select a different internal wall grid, so the
extractor checks its actual count and domain against the request. `grid.json`
stores native mesh-table faces; `PROF CELL_CENTERED=T` explicitly stores
thermal cell-center values. FDS profiles for V-HC-01 and V-HC-02 are written at
every output (`FORMAT_INDEX=1`); the others store only the final profile.

Most Gpyro decks are generated from the case parameters. Physics that
generator does not express (a surface-absorbed flux, in-depth absorption,
apparent-heat-capacity latent heat, a two-step chain, retained reaction gas,
mass-transfer films and Darcy inflow) uses the native V2 inputs in
[`gpyro_decks/`](gpyro_decks/), each at its base level; a level changes only
`DT0`, `NCELLZ(1)` and `TSTOP(1)`, and the script checks the base values
against the case parameters. Those decks follow the [Gpyro Verification and
Validation Manual v1.0](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/documentation/gpyro_verification_manual.pdf)
(§4.4 diffusion, §4.3 and Appendix A.2 Darcy flow, §3.6 reaction enthalpy)
and Lautenberger's [2007 dissertation](https://escholarship.org/uc/item/7wz5m7dg)
§§3.4.1–3.4.4, adapted to the problems below; they are not unchanged upstream
runs. Gpyro's 1D grid places node 1 at the exposed face and the last node at
z = L, with half control volumes at both, so its profiles need no half-cell
shift.

* **V-HC-01:** L = 0.05 m, k = 0.2, ρ = 1000, c = 1400, absorbed flux
  25 kW/m² with h = 20 W/(m² K) to 300 K, insulated back. Base: 500 intervals,
  dt = 0.1 s; profiles at 30, 60, 90 and 180 s and the surface history
  (5 s FDS outputs, 1 s Gpyro outputs). FDS has no pure net-flux face on a
  conducting material, so the face convects to T∞ + q″/h = 1550 K with
  ε = 0, the same Robin condition. Gpyro absorbs `QE` at the surface node
  (`EMIS=1`, `KAPPA=9D9`, no re-radiation).
* **V-HC-02a–d:** L = 0.1 m, ρ = c = 1000, k = 0.5,
  T₀ = 300 K, T∞ = 500 K and h = 0.5/5/50/500 W/(m² K), insulated back.
  Base: 80 intervals (81 Gpyro boundary nodes), dt = 2 s. Profile times are
  1000, 2000, 4000 and 10000 s.
* **V-HC-03 (Gpyro):** L = 0.02 m, k = 0.5, ρ = c = 1000, absorbed flux
  5000 W/m² at the surface node, h = 0, no re-radiation, insulated back.
  Base: 321 boundary nodes, dt = 0.1 s; profiles at Fo = 0.1, 0.5 and 1.5.
  The volume mean integrates the node values over half end cells (80 s
  outputs).
* **V-HC-04a–d:** both solvers use L = 0.02 m, 300/400 K boundaries,
  ρ = c = 1000, bottom k = 1 and top k = 0.1/0.01/0.1 W/(m K),
  interface z = 0.01/0.01/0.0102 m. Native Gpyro obstacle orientation is
  explicit (`GEOMETRY_IS_UPSIDE_DOWN=F`); both blocks cover the domain.
  Base: 80 requested intervals; dt = 10/100/10 s; final time 40000/400000/40000 s.
  Both FDS layers have explicit `CELL_SIZE(layer)` and `STRETCH_FACTOR(layer)`.
  The native grid is retained;
  off-center layer thicknesses can give slightly different widths per layer.
  Each numerical interface is compared against the prescribed physical
  interface, with no coordinate correction. Off-face quantization remains.
  V-HC-04d uses the same physical 10.2 mm interface as V-HC-04c, but Julia
  places a mesh face there using 20 cells in each layer. Its manifest entries
  explicitly reuse `vhc04c.csv`, including the native refinement runs; no
  extra native runs or improved native interface alignment are implied.
  Julia reports both prescribed-geometry and represented-geometry errors,
  retaining V-HC-04c as the off-face geometry-resolution test.
* **V-HC-05:** L = 0.02 m, boundaries 300/800 K, k(T) = 0.2 + 0.002T,
  ρ = c = 1000. FDS uses a linear conductivity ramp. Gpyro uses equal
  amounts of two inert constituents of identical density and heat capacity:
  k₁ = 0.4 and k₂ = 1.2(T/300). Its volume average is exactly the same
  affine law. Base: 80 intervals, dt = 1 s, final time 5000 s.
* **V-RX-02a–c:** A converts entirely to an equal-density, equal-heat-capacity
  solid B, with no reaction heat, transport or volume change. Initial A is
  100/100/1000 kg/m³, order is 0/2/(2/3), and normalized rate coefficient is
  0.05/0.1/[3.4e11 exp(-103000/(8.314462618 × 450))]. FDS multiplies this
  coefficient by ρ^(1−n), because its rate uses concentration powers;
  Gpyro uses normalized concentration. The last coefficient evaluates the
  Arrhenius law at the prescribed 450 K: it matches this isothermal case,
  not a varying-temperature Arrhenius test. No empirical rate rescaling is
  used. Base: one interval, dt = final time/10000, 100 output intervals.
  All runs stop before depletion, where limiter conventions would differ.
* **V-RX-01:** A → B → C with k = 0.389 and 0.262 s⁻¹, E = 0, unit
  condensed yields, ξ₀ = 1000 kg/m³, held at 400 K. Base: dt = 0.001 s,
  0.1 s outputs. FDS uses one insulated 1 mm cell with three materials of
  equal density. Gpyro runs `ZEROD` with one lumped unknown at
  `TAMB=400`, `BETA=0`; its mass fractions are scaled by ξ₀, and the script
  requires M/M₀ = 1 and T = 400 K throughout.
* **V-RX-04a:** L = 1e-4 m, k = 100, ρ = c = 1000, T₀ = 900 K,
  ε = 0.9, radiation environment 300 K, no convection, insulated back.
  Base: four intervals, dt = 0.002 s, final time 8 s. Slab means use
  native control-volume weights: trapezoidal weights on Gpyro boundary nodes
  and equal cell-volume weights on FDS probes placed at its thermal cell
  centers (four base / eight refined cells). The analytical
  reference is the thermally thin cooling limit, as in the Julia case.
* **V-RX-04b:** L = 0.005 m, k = 1, ρ = c = 1000, T₀ = T∞ = 300 K,
  ε = 0.85, h = 15 and absorbed irradiance 20000 W/m², insulated back.
  Base: 20 intervals, dt = 0.2 s, final time 2000 s. Gpyro incident flux
  is q_abs/ε. FDS uses `SKIP_INRAD=T` and external incident flux
  q_abs/ε + σT∞⁴, with radiation enabled and no prescribed `TMP_GAS_FRONT`
  (that option bypasses external flux). The analytical line is the steady
  root, not an exact heating transient. Both native solvers use
  σ = 5.670373e-8 versus Pyrolysis.jl's 5.670374419e-8 W/(m² K⁴).
* **V-RX-05a–b:** FDS uses the insulated A → B problem, rate
  0.2 s⁻¹, initial density 1000 kg/m³, L = 0.001 m, T₀ = 1000 K,
  heat capacities 1200/1900 J/(kg K), and reaction heats +1e6/−1e5 J/kg.
  `ADJUST_H=F` on both materials disables FDS's default reference-temperature
  adjustment of reaction enthalpy. This matches the specified constant local
  reaction heat (see [FDS input processing](https://github.com/firemodels/fds/blob/FDS-6.11.1/Source/read.f90)
  and [solid energy equation](https://github.com/firemodels/fds/blob/FDS-6.11.1/Source/wall.f90)).
  Base: one interval, dt = 0.001 s, final time 20 s. Sensible enthalpy is
  reconstructed from native concentrations and temperature using the common
  298.15 K datum and unit area, rather than comparing different native datums.
* **V-RX-05a–c (Gpyro):** the same kinetics, heats and heat capacities
  (c_G = 1000 J/(kg K)); member c produces 0.4 B + 0.6 G, encoded by
  `R0(B)=400`. Two boundary nodes span the insulated 1 mm slab with no
  deformation, zero permeability and no boundary mass flux; Gpyro's `ZEROD`
  mode prescribes temperature, so the insulated energy balance needs this
  uniform slab. Both species' true solid densities are 1000/0.7 kg/m³,
  matching the evolving pore volume of our fixed mesh, and the pressure
  solve retains all generated gas. A 1e-6 Pa seed pressure avoids a zero
  native diffusivity at vacuum; its gas mass is about 1e-15 kg/m². `CI` and
  `CJ` give the condensed and gas concentrations, and the enthalpy is
  L (c_A ξ_A + c_B ξ_B + c_G ξ_G)(T − 298.15). Member c retains a native
  total-mass offset of about 3.27e-5 of the initial mass that persists when
  dt is halved: the [pressure equation](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_gas.f90)
  uses R = 8.314462 while [`RHOGOFT`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_funcs.f90)
  uses 8.314. The comparison keeps it; no rescaling is applied.
* **V-CPL-01 (Gpyro, route A):** the Neumann two-phase solidification
  problem with apparent heat capacity: k = 0.2, ρ = 900, c = 1700,
  T_m = 433 K, ΔH_m = 50 kJ/kg, σ_m² = 0.1 K², T₀ = 438 K and a 428 K fixed
  face (`TFIXED`), L = 0.04 m, insulated back. Base: 801 boundary nodes,
  dt = 0.02 s, 1 s profiles. Histories are the nodes at 0.5, 1.0, 1.5 and
  2.0 mm; the front is the T_m isotherm interpolated linearly between the
  first bracketing nodes. Route B (different temperatures, latent heat and a
  gated reaction) has no native deck.
* **V-CPL-03 P1/P2/P3:** both solvers use solid A → solid B at 400 K, rate 2 s⁻¹,
  zero reaction heat, initial density 1000 and product density 500/2000/1000
  kg/m³, L = 0.002 m. Gpyro uses `CHI=0`, `SOLVER_DEFORMATION_MODE=3`;
  there is no gas product or pressure solve. Base: eight intervals,
  dt = 0.001 s, final time 8 s. FDS outputs instantaneous thickness samples
  for all three paths.
  The native decks are insulated: uniform initial temperature, equal heat
  capacities and zero reaction heat keep T exactly constant, making the
  Julia case's prescribed 400 K front boundary redundant for these paths.
  The full thickness history is compared with 1 + (Vf/V₀−1)(1−exp(-2t)).
* **V-GT-07a–b (Gpyro):** both gases have M = 32 g/mol; T = 300 K,
  P = 1e5 Pa, porosity 0.6, D = 2e-5 m²/s and external tracer fraction 0.22,
  film coefficients 1e-4 and 1e-2 kg/(m² s). The Chapman–Enskog parameters
  `EPSOK=100 K` and `SIGMA=3.588275763111019 Å` give D exactly; the script
  checks the native D12 profile. Temperature is held by
  `FIX_DOMAIN_TEMPERATURE` with `BETA=0`. Base: 161 boundary nodes,
  dt = 0.05 s; profiles at 50, 150 and 400 s. Gpyro clamps `HM=0` to
  1.19e-7 kg/(m² s) in [`INTERPOLATE_BCS`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_bc.f90),
  so the far reservoir holds only the carrier; its exchange stays below
  2e-22 kg/(m² s). The influx is h_m (0.22 − Y_surface) and the inventory
  integrates the nodes over half end cells at porosity times the native gas
  density, whose `RHOGOFT` uses R = 8314 J/(kmol K) against the benchmark's
  8314.46; this small difference is retained.
* **V-GT-08a–b (Gpyro):** L = 0.1 m, permeability 5e-10 m², porosity 0.6,
  T = 300 K, outlet pressure 1e5 Pa and inflow 10 kg/(m² s)
  (`MDOTPP=10000` g/(m² s)); member a feeds the bottom, b the top. Gpyro uses
  its Chapman–Enskog D12 as kinematic viscosity in the Darcy coefficient:
  M = 29 g/mol, `EPSOK=100 K` and `SIGMA=3.965473809316032 Å` give
  D12 P = μRT/M for μ = 2e-5 Pa s, which the script checks. Base: 81
  boundary nodes, dt = 0.005 s; the steady profile is compared at 1 s and
  the fluxes at 0.5 and 1 s. Native gauge pressure gains 1e5 Pa. The top
  outward flux is minus the native z flux there; the bottom outward flux is
  the native integrated MLR minus the top one, both in kg/(m² s). The t = 0
  dump is omitted: Gpyro has not yet applied the inflow in it.
* **V-CPL-04:** both solvers use self-regenerating A → A, zero order,
  reaction heat −1e5 J/kg, L = 0.02 m, k = 0.2, ρ = c = 1000,
  initial/boundary T = 500 K, E/(RT₀) = 40 and δ = 0.44. The pre-exponential
  is derived from the case's volumetric source. Native activation energies
  preserve E/R using FDS R = 8.314472 J/(mol K) and Gpyro's compiled
  R = 0.008314000152051449 kJ/(mol K). This is a unit/constant conversion,
  not a fit. FDS also uses `ADJUST_H=F`. Base: 80 intervals, dt = 0.25 s,
  final time 2500 s. The compared temperatures stay near 500–504 K, where
  the Julia upper-temperature gate at 600 K is inactive to floating-point
  precision. No runaway or critical-threshold comparison is implied.

## ThermaKin2Ds mappings

ThermaKin2Ds uses uniform, cell-centered elements, so N elements match N
Julia cells of the same width; the first input layer is the top (front)
face and z = L − depth. Each run takes fixed Crank–Nicolson steps with one
linearization per step. `DURATION` is half a step short of the horizon so
that the loop, which runs while k·dt < `DURATION`, ends exactly on it.
Outputs print seven significant digits (about 1e-4 K at 500 K) and every
element at every output. The retained `.out` files are gzipped with the
start-date and total-runtime lines removed, so a rerun reproduces them
byte for byte.

The mappings below are native inputs, not corrections:

* **Boundary location.** Boundary heat and mass exchange act on the
  boundary element's center value. Surface and back-face histories are
  therefore element-center values Δx/2 inside the face, and a fixed or
  film boundary is applied there; the resulting O(Δx) offsets stay in the
  errors.
* **Fixed temperature** (V-HC-04, V-HC-05, V-CPL-04): convective coefficient
  1e9 W/(m² K), with the boundary element initialized at the fixed value.
  Crank–Nicolson does not damp the large-coefficient mode, so an initial
  jump at that element would oscillate indefinitely.
* **Radiation.** Each boundary deposits external radiation in, and emits
  εσT⁴ from, the element selected by its absorption search. The exposed
  element is a separately named copy of the material (`A_surface`) carrying
  the case emissivity; every other element has ε = 0. A mass absorption
  coefficient of 1e6 m²/kg saturates each boundary element, so deposition
  stays in the exposed element and the back boundary's search stops at its
  own non-emitting element. V-RX-04a/b use incident flux q_abs/ε + σT∞⁴ with
  the native σ = 5.6704e-8 W/(m² K⁴).
* **Absorbed flux without emission** (V-HC-01, V-HC-03): ThermaKin has no
  net-flux boundary, so the incident flux is q″/ε on an exposed element with
  ε = 1e-8; its emission loss stays below 1e-3 W/m² below 1200 K.
* **Reactions.** The rate is A·exp(−E/RT)·ξ₁·ξ₂, first order in each of at
  most two reactants. V-RX-02b is A + A with stoichiometry ½, ½ and
  A = k ξ₀^(1−n). V-CPL-04 is a self-regenerating A → A, so A·ξ_A with
  constant ξ_A is the zero-order source; its 600 K upper gate is
  `TEMP LIMIT: U 600`. Activation energies preserve E/R with the native
  R = 8.31447 J/(mol K). ThermaKin's positive reaction heat is exothermic,
  so the sign is reversed.
* **Gases.** ThermaKin's element volume is Σm/ρ over the condensed phases
  plus γ times the gas volumes. V-GT-01 and V-GT-07 use solid swelling 0
  and gas swelling limit 1e-30, so gases occupy no volume and ξ is per
  element volume. V-RX-05c instead uses solid swelling 1 (γ = 1 exactly)
  and gives the retained gas the solid density, so converting solid to gas
  leaves the element volume fixed, as for gas held in the pores. The
  transport law J = −ρ_g λ ∂(ξ/ρ_g)/∂x with constant ρ_g is
  Fickian diffusion of ξ with coefficient λ, set to the case diffusivity on
  every component. ThermaKin assembles gas transport only when a boundary
  lists a mass flux, so the sealed V-GT-01 slab carries a `LIN 0 0` entry.
  V-GT-01's Gaussian initial state is one layer per element at the element
  center; the 30000-character input limit allows at most 320 such layers,
  so its base mesh is 160 elements. V-GT-07 sets the tracer density to
  φρ_g, making ξ/ρ the tracer mass fraction and the linear boundary law
  a·ρ·(ξ/ρ − b) the film flux with a = h_m/(φρ_g) and b = Y∞; the inert
  carrier only fixes the pressure, which ThermaKin does not model. The
  printed boundary mass flow is the uptake flux on the unit area.
* **Conductivity mixing** is series (`PARALL CONDUCTIVITY: 0`), exact for
  layers normal to the heat flow; only the V-HC-04c element that straddles
  the 10.2 mm interface contains both materials.
* **Steps.** Base steps equal the Gpyro steps, except V-HC-05 at 0.05 s:
  with k(T), ThermaKin's linearized step diverges at 0.1 s and above.

ThermaKin cannot express V-CPL-01
(properties are p₀ + p₁T + p₂Tⁿ, which has no Gaussian latent-heat peak),
V-RX-02a/c (zero and two-thirds order), or the pressure-driven cases
V-GT-02 to V-GT-06 and V-GT-08.

The refinement audit is not an assertion of monotone convergence. Steady
profile errors include any residual approach to equilibrium at the fixed
output time; output precision and native constants also impose floors.
No convergence order is inferred automatically. V-RX-04b's whole-history error includes the
physical approach to its steady reference, so it is not a transient
discretization-error estimate.
