# Matched conduction, reaction, deformation and gas-diffusion comparisons

These accuracy comparisons use native inputs to unmodified FDS 6.11.1 and
Gpyro V2, and to ThermaKin2Ds 7 built from its source with a banner-only
portability change. There are 56 solver/subcase pairs, each run at five
configurations (280 successful runs), spanning 25 verification case IDs.
V-HC-04c and V-HC-04d share the same physical problem and explicitly reuse
the same native runs: only Julia's mesh differs. P1/P2/P3 are counted
separately within V-CPL-03. This dataset contains no pressure cases or
timing measurements.
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
(h,dt/2), `time_quarter` (h,dt/4), and `space_half` (h/2,dt/4). The coupled refinement is not used to infer isolated spatial or temporal order.
Use `--cases=vhc04a --levels=base,time_half,time_quarter,space_half
--destination=/tmp/comparison-pilot` for a separate pilot; partial studies
cannot overwrite this directory's complete manifests.
The script verifies generated inputs, successful termination (not merely an
exit status), pinned version banners, finite numeric outputs, monotone times
and the final output time. It writes base samples in `fds/` and `gpyro/`, finer
samples in `refined/`, and native inputs/output CSVs in `raw/`. Native timing
files and logs containing elapsed times are not published. Each run retains
a hashed `grid.json` containing actual initial coordinates, indexing convention,
interval count and spatial unknown count. FDS coordinates come from its native
mesh table, not probe count; Gpyro coordinates come from its native profile.
The TOML manifests
bind each series to an exact case, curve and named analytical reference.
Report generation fails on a missing binding, missing samples or nonfinite
analytical errors. The extracted coordinates are never fitted or shifted to
improve agreement.

The provenance JSON files hash the executable, manifests, retained decks and
outputs. FDS is the official Linux release
`FDS-6.11.1-0-gff928db-release`, compiled with Intel Fortran 2025.3.2. Gpyro is
[revision bff0126c68f2a35c084ef04490ae327185f281a3](https://github.com/Gedeben/GPYRO-V2/tree/bff0126c68f2a35c084ef04490ae327185f281a3),
built with GNU Fortran 14.4.0 and upstream `gnu_linux` flags; build details and
the source archive hash are in [the V2 build documentation](../gpyro/v2_accuracy/README.md).
Each run requests one OpenMP thread and FDS uses one MPI rank. Gpyro's banner
reports OpenMP disabled in this build.

Regenerate accuracy reports, including per-case `comparison_accuracy.csv`
with maximum and RMS error against the shared analytical functions:

```sh
julia --project=scripts/verification scripts/verification/run_reports.jl \
  --cases=V-HC-02,V-HC-04,V-HC-05,V-RX-02,V-RX-04,V-RX-05,V-CPL-03,V-CPL-04 \
  --timings=false --comparisons-only=true --outdir=/tmp/verification-comparison-report
```

Only base external curves are plotted. The other runs are retained to audit
resolution sensitivity without crowding the plots. The command above
produces accuracy-only partial chapters. Runtime comparisons require a
shared achieved-accuracy target. Omit
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
thermal cell-center values.

* **V-HC-02a–d:** Gpyro uses L = 0.1 m, ρ = c = 1000, k = 0.5,
  T₀ = 300 K, T∞ = 500 K and h = 0.5/5/50/500 W/(m² K).
  Base: 81 boundary nodes, dt = 2 s. Profile times are 1000, 2000,
  4000 and 10000 s. FDS curves use the equivalent
  dimensionless slab problem.
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

ThermaKin cannot express V-RAD-01 (absorbed radiation is deposited in one
element per step, not distributed by Beer–Lambert absorption), V-CPL-01
(properties are p₀ + p₁T + p₂Tⁿ, which has no Gaussian latent-heat peak),
V-RX-02a/c (zero and two-thirds order), or the pressure-driven cases
V-GT-02 to V-GT-06 and V-GT-08.

The refinement audit is not an assertion of monotone convergence. Steady
profile errors include any residual approach to equilibrium at the fixed
output time; output precision and native constants also impose floors.
No convergence order is inferred automatically. V-RX-04b's whole-history error includes the
physical approach to its steady reference, so it is not a transient
discretization-error estimate.
