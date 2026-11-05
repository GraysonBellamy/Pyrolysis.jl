# Comparison coverage

External curves are included when native solver inputs represent the same
equations, material laws, boundary conditions and observable. A similar example
name or matching endpoint is insufficient. Comparisons use the shared analytical
reference at each solver's own output coordinates, without fitting simulation
results. Different discretizations remain visible as numerical error.

Comparisons use unmodified FDS 6.11.1, Gpyro V2 revision
`bff0126c68f2a35c084ef04490ae327185f281a3`, and ThermaKin2Ds 7 built from
source with a start-time-banner portability change only. Every external curve
comes from one matched study, run at seven configurations: time-only runs at
2dt, 4dt, dt/2 and dt/4, spatial refinement at fixed dt/4, and coupled
mesh/time refinement. Native inputs, outputs, executable hashes, extraction
code and equation mappings are under `test/verification/reference/matched/`.
Each case's comparison section links its analytical-error CSV; the same runs'
solve times give its cost at matched accuracy. Matching equations does **not** establish
matching discretization, observation operators, or achieved accuracy.

## Comparison basis audit

The [per-pair inventory](assets/comparison_basis.csv) records the actual Julia
mesh used in the published overlay, the external intervals and spacing, and
the time refinement retained for each pair. The
[time-refinement differences](assets/comparison_time_changes.csv) isolate
temporal sensitivity at fixed mesh. These differences have native observable
units; they are not convergence orders.

The inventory covers **77 case/code/subcase pairs across 29 case IDs**:
19 have matching plotted spacing, 25 have different plotted meshes, 27 are
spatially homogeneous ODE comparisons, and 6 compare a material interface
that each code represents differently. There are 511 retained native runs (74 native decks × seven
configurations, 7 recorded unavailable).
The Julia study contains 108 solves across 36 case/mesh/variant workloads;
all three tolerance levels reached their expected return codes and the
accepted-step observer counts agreed with the solver statistics.

The comparison basis and its limits are:

| Comparisons | Mesh and observations |
|:--|:--|
| HC-01 / Gpyro V2 | 501 boundary nodes span 500 intervals at 0.1 mm, matching Julia's plotted mesh. Independent time/space refinements are retained. |
| HC-01, HC-03 / ThermaKin | 500 and 320 uniform elements match Julia's plotted cells exactly. The surface value is the boundary-element center, Δx/2 below the face. |
| HC-01 / FDS | 500 cells at 0.1 mm, matching Julia's plotted mesh, with instantaneous devices and full cell-center profiles. The absorbed flux with convective loss is convection to the equivalent 1550 K gas. |
| HC-02a–d | Julia's plotted 320 cells differ from the 80 external intervals of every code. Julia's resolution ladder includes an 80-cell mesh; the plots use the finest mesh. |
| HC-03 / Gpyro V2 | 321 nodes span 320 intervals at 0.0625 mm, matching Julia's plotted mesh. Independent refinements are retained. |
| HC-04a–d / FDS | Both layers have explicit `CELL_SIZE` and `STRETCH_FACTOR` controls. The base grid has 80 cells; Julia's plotted mesh has 40, and its cost ladder includes 80 cells. |
| HC-04c/d / both codes | Equal total count is insufficient: prescribed and represented material interfaces differ. FDS follows physical layer thicknesses; Gpyro samples obstacle composition on boundary nodes. c is a geometry-quantization experiment; d aligns a Julia mesh face with the interface and reports both prescribed- and represented-geometry errors. |
| HC-05 | Julia's plotted 320 cells versus 80 external intervals; the cost ladder includes the 80-cell rung. |
| RAD-01 / Gpyro V2 | 401 nodes span 400 intervals at 0.025 mm, matching Julia's plotted mesh. Independent time/space refinements are retained. |
| RX-01/02 and RX-05 | Spatially homogeneous ODEs: one Julia cell, two Gpyro boundary nodes, or native 0D are physically comparable, but not equal work. Cost compares time at equal error on the same time observations; no spatial convergence is claimed. RX-05c retains its native gas-constant error floor. |
| RX-04a/b | Base spatial spacing matches (4/20 Julia cells versus 4/20 native intervals). Slab-mean quadrature stays fixed across time refinement. Cooling uses a thermally thin approximation; heating has only a steady exact reference. A whole-history error against the steady line is not transient numerical error. |
| CPL-01 route A | The plotted 800-cell Julia mesh is uniform, matching V2's 801 boundary nodes. Every native node contributes to front interpolation. Independent refinements retain the finite Gaussian-peak and sharp-interface-reference distinction. Route B remains excluded. |
| CPL-03 P1/P2/P3 | Base spacing matches 8 Julia cells and 8 external intervals. Timing covers the same paths and output cadence, without the report's unrelated QoI sub-solves. |
| CPL-04 stable branch | Julia's plotted 160 cells versus 80 native intervals; the cost ladder includes the 80-cell rung. Timing covers the δ = 0.44 solve only, not the report's bifurcation sweep. Profiles are compared at the common final time. |
| GT-07a/b and GT-08a/b | Base spacing matches (160/80 intervals). Independent time and space refinements isolate their respective sensitivities. Cell-average and point-value references stay separate, and native constant differences are retained. |
| GT-01 / ThermaKin | 160 elements versus Julia's plotted 640 cells; Julia's ladder has a 160-cell rung. ThermaKin's 30000-character input limit caps the one-layer-per-element Gaussian initial state at 320 elements. |
| ThermaKin boundary conditions | Fixed temperatures, films and fluxes act on boundary-element centers, so fixed-temperature and film cases carry an O(Δx) boundary offset that refinement reduces but the base comparison retains. |

### Spatial comparison rule

Use the same domain, boundary locations, material interfaces, initial state,
physical horizon and output quantities. For a uniform domain of length L,
Julia's N cells have N+1 faces and h = L/N. Gpyro with boundary-centered
unknowns needs N+1 nodes to have the same interval spacing; its two boundary
control volumes are half width. Equal spacing does not mean equal degrees of
freedom or identical boundary discretization. Record both counts.

There is a pinned-V2 two-node exception: the matched RX-02 native profile
prints coordinates `[0, 3L/4]`. In `COMPUTE_Z_CELL_POSITION`, the first/last
position formulas overlap when `NCELLZ=2`; `DZN` still contains two widths
L/2, so the total control volume remains L. The extractor preserves these
coordinates and flags the exception instead of shifting them. These retained
runs are homogeneous ODE checks only. Do not use a two-node V2 deck as evidence
of a matched spatial transport mesh. The V2 RX-05 two-node decks print the
same coordinates and likewise do not qualify a spatial transport comparison.

FDS's `PROF` output defaults to temperatures interpolated to **faces**, even
though the thermal unknowns are cell centered. Its printed "Solid Phase Node"
coordinates are the cell faces; N+1 such coordinates describe N cells.
This is established by the pinned [FDS profile writer](https://raw.githubusercontent.com/firemodels/fds/FDS-6.11.1/Source/dump.f90)
and [mesh setup](https://raw.githubusercontent.com/firemodels/fds/FDS-6.11.1/Source/read.f90).
Never count temperature probes as cells or shift all native coordinates by
half a cell. Retain actual grid coordinates and specify whether each value
is a center value, face reconstruction, point sample or cell average.

The matched FDS runner requests `CELL_CENTERED=T`. HC-04a's 80-cell mesh
has a maximum native-center profile error of 4.55e-6 K, near CSV precision.
Face-interpolated values at a conductivity jump use a different observation
operator and must be evaluated separately. Radiative cooling uses a
volume-weighted mean from native FDS cell-center probes, cross-checked
against the final native center profile.

For stretched/layered meshes compare the full face list, widths near the
heated surface and reaction front, widths in each layer, and the represented
interface. For moving meshes match the initial geometry and separately
compare subsequent geometry; equal initial N is not equal resolution for all
time. Evaluate conserved inventories using each control volume's measure.
Analytical cell-average references and point references are different
observations and should not be pooled into one error norm.

### Temporal comparison rule

Keep each code's native integrator. An adaptive average step
`(t_end-t_start)/accepted_steps` is useful work metadata, but is **not** a
matching accuracy criterion: short ignition, diffusion-startup or front
events may require far smaller steps. `diff(saved_output_times)` is not the
integrator's step history. Likewise, native `DT`/`DT0` is a request: FDS can
substep the solid calculation and Gpyro can reduce failed steps.

For example, the 80-cell Julia HC-05 baseline took 588 accepted steps over
5000 s: mean dt = 8.50 s, median dt = 0.0815 s, minimum dt = 0.0001 s and
maximum dt = 2113 s. The large late steady-state steps dominate the average;
a fixed 8.50 s reference run would not reproduce its startup resolution.

The native runner freezes each external code's mesh, observation times,
output cadence, quadrature, iteration tolerances and equations, and runs
`(h,dt)`, the coarser `(h,2dt)` and `(h,4dt)`, the finer `(h,dt/2)` and
`(h,dt/4)`, and the finer mesh at `(h/2,dt/4)` and `(h/2,dt/2)`. Native
step histories are not retained: output intervals are not integration steps.
Every code is sampled instantaneously and scored against the analytical
reference at its own output points. Profiles come only from an output at the
case's observation time, never from another output and never by
extrapolation, so a coarse configuration that misses an observation time is
recorded as unavailable. Julia's time study runs each pair's native interval
count at three tolerance levels, with step caps from its own baseline run,
and records its step and work counters (see
[Reproduce the audit](#reproduce-the-audit)).

The [time-refinement differences](assets/comparison_time_changes.csv)
record each code's change from dt to dt/2 and from dt/2 to dt/4 at fixed
mesh, as observed sensitivity, not an order or an upper error bound. No
temporal allowance is subtracted: each run's measured error includes its
temporal error, and the cost comparison below selects runs by that measured
error. Approach to a steady reference is not transient discretization error.

### Runtime at matched accuracy

Each case's cost section plots solve time against analytical error for every
configuration of each external code and of this code. This code's
configurations are its resolution ladder, each pair's native interval count
and any finer cost-only meshes, each at the report tolerances and at 10, 100
and 1000 times looser ones. The looser tolerances are the counterpart of the
external codes' coarser steps; a level that fails is left out. A run's error
is its largest analytical error over the curves and references it is scored
on, at its own output coordinates. For each external configuration, the
table lists that code's cheapest configuration at least as accurate, this
code's cheapest configuration at least as accurate on every scored
reference, and their ratio. Errors at or below 1e-12 count as exact. A ratio
marked ≥ uses this code's coarsest configuration (fewest cells at the
loosest tolerance), which is already more accurate; "not reached" means no
configuration of this code is as accurate.

Times are medians of repeated runs, serial, on one pinned core of an idle
machine (`scripts/timing/README.md`):
FDS's internal elapsed time, which excludes its MPI start-up; Gpyro's and
ThermaKin's process time; and this code's warm problem construction and
solve after a discarded compiling solve. Every timed external run reproduced
its committed native outputs byte for byte, so each time belongs to the run
whose error the page reports. Times include native output, and Gpyro's
include its one step past the last observation. Plotting, analytical
references and the reports' unrelated sweeps are not timed. The machine,
commits, executable hashes, timer boundaries and disturbance record are in
`test/verification/reference/matched/timing/provenance.json`.

A ratio is a cost at equal accuracy for that problem, not a pure solver
speedup. On the smallest problems the external times are mostly input,
setup and output, and the codes' discretizations, step control and output
workloads differ. Approximate references (the thermally thin cooling limit,
the sharp-interface Stefan solution, steady profiles) put the same model
error into every code's error.

### Reproduce the audit

```sh
julia --project=scripts/verification scripts/verification/export_comparison_meshes.jl \
  test/verification/reference/comparison_meshes.toml
python3 scripts/verification/audit_comparison_basis.py
julia --project=scripts/verification scripts/verification/run_comparison_time_study.jl
python3 -m unittest discover -s scripts/verification -p test_comparison_basis.py
```

`matched/run_comparisons.py` runs all seven levels by default for FDS, Gpyro
V2 and ThermaKin2Ds; `matched/build_thermakin.sh` builds the ThermaKin
executable from a local, uncommitted source copy. Its `--cases`,
`--levels` and `--destination` options support isolated pilots; a partial
study cannot replace the complete reference manifest. Each matched run stores a hashed
`grid.json`; the audit checks every retained native-input/output hash.

The Julia time-study runner uses each pair's external interval count (separate
P1/P2/P3 problems for CPL-03), three tolerance levels, and maximum accepted-step
caps derived from its own baseline run. Its read-only callback records actual
accepted dt and verifies its count against solver statistics. Results are
[step/work counters](assets/comparison_julia_steps.csv) and
[changes at fixed mesh](assets/comparison_julia_time_changes.csv). These use
the case's saved times; they do not by themselves establish dense common-time
observation accuracy or repair the interface exceptions. No wall times are
measured by this runner.

## Gpyro V2 accuracy

**All Gpyro verification comparisons use V2 revision
`bff0126c68f2a35c084ef04490ae327185f281a3`.**

The HC-01, HC-03, RAD-01 and CPL-01 spatial base grids match the plotted Julia interval spacing:
501 Gpyro boundary nodes versus 500 Julia cells for HC-01; 321 versus 320
for HC-03; 401 versus 400 for RAD-01; and 801 versus 800 for CPL-01.
RX-01 is a `ZEROD` chain with one lumped unknown, and RX-05 a two-node
homogeneous slab. Stefan front interpolation uses every node in the full
native profiles. Each run stops one requested step past its last
observation so that the final output is written; the timed runs include that
step.

The finite Gaussian latent-heat peak and spatial front interpolation keep
CPL-01's front offset from the sharp-interface reference at every
configuration. RX-05c retains a native total-mass offset of about 3.27e-5
of the initial mass from inconsistent gas constants; halving dt leaves it
unchanged. Each case's comparison-accuracy CSV records the analytical errors
at all seven configurations.

## ThermaKin2Ds accuracy

ThermaKin2Ds 7 contributes 24 native cases (165 runs, 3 recorded
unavailable) to the matched dataset. Its elements are cell centered, so HC-01, HC-03 and RX-04a/b use
exactly Julia's plotted cells. Halving its fixed Crank–Nicolson step changes
the outputs by at most a few units in the seventh printed digit, so its
errors are spatial or observational.

The homogeneous cases sit at its output precision: 5.6e-5 kg/m³ for
RX-01, 4.9e-6 kg/m³ for RX-02b, 4.9e-4 K and 5e-5 kg/m³ for RX-05a–c and
2.8e-7 relative thickness for CPL-03. The flux-driven HC-03 profiles
(1.7e-4 K) and the sealed GT-01 profiles (5.6e-5 kg/m³) show second-order
spatial error of the same size as Julia's at the same mesh. ThermaKin applies fixed temperatures, convection and mass-transfer
films at the boundary element's center, Δx/2 inside the face. That location
gives first-order errors which halve under spatial refinement: 1.1–1.2 K for
HC-04a–d, 5.1 K for HC-05, 0.11 K for CPL-04, up to 3.1 K for the HC-02
profiles, 5.1 K for the HC-01 surface history and 1.2e-3 to 1.2e-2 in tracer
fraction for GT-07a/b. The [time-refinement
differences](assets/comparison_time_changes.csv) and each case's
comparison-accuracy CSV record all seven configurations.

## Solver coverage

| Verification problem | FDS | Gpyro V2 | ThermaKin2Ds | Scope |
|:--|:--|:--|:--|:--|
| V-HC-01: semi-infinite solid, flux and convection | Available | Available | Available | Profiles and surface history |
| V-HC-02a–d: convective slab / Biot sweep | Profiles and histories | Profiles and histories | Profiles and histories | Four profiles and front/back histories |
| V-HC-03: finite slab, constant flux | Excluded | Available | Available | Profiles and volume-mean history |
| V-HC-04a–d: layered conduction | Available | Available | Available | Steady profile against the prescribed interface; c/d share native runs, with only Julia's mesh aligned in d |
| V-HC-05: temperature-dependent conductivity | Available | Available | Available | Steady Kirchhoff profile; same affine conductivity law |
| V-RX-01: series chain | Available | Available | Available | Three concentration histories |
| V-RX-02a–c: reaction order | Direct inputs | Direct inputs | b only | Zero-, second-, and two-thirds-order concentration histories |
| V-RX-04a–b: radiative boundaries | Available | Available | Available | Slab-mean cooling and surface approach to radiative/convective equilibrium |
| V-RX-05a–c: unequal heat capacities / reaction heat | a, b | a–c | a–c | Temperature, concentrations and sensible enthalpy; FDS requires `ADJUST_H=F` |
| V-CPL-03: volume change | P1/P2/P3 | P1/P2/P3 | P1/P2/P3 | Pure solid swelling, shrinking and constant-volume paths only |
| V-CPL-04: thermal feedback | Partial | Partial | Partial | Subcritical steady profile at δ = 0.44 only; no claim about the gated runaway branch |
| V-GT-01: Gaussian gas diffusion | Excluded | Not run | Available | Profiles and variance growth |
| V-GT-07a–b: diffusion through a mass film | Excluded | Available | Available | Profiles, uptake and inventory |

Other retained comparisons include V-RAD-01 (Gpyro V2), the
apparent-heat-capacity route of V-CPL-01 (Gpyro V2) and Gpyro V-GT-08.

V-HC-04c remains an off-face geometry-resolution test. V-HC-04d places a
Julia mesh face at the same physical 10.2 mm interface; both report prescribed
and represented geometry errors. Reusing the native runs is explicit in the
manifests because the materials, geometry and boundary conditions are identical.

V-GT-07 reports pointwise errors and analytical cell-average errors for
Julia's finite-volume inventories. Gpyro boundary-node outputs
retain point references. Resolution figures show Julia's maximum profile,
surface-flux and inventory errors over 50, 150 and 400 s against mesh
spacing. The external codes' errors at every configuration are in the case's
comparison-accuracy CSV and cost figure; their `time_half`, `time_quarter`
and `space_half` runs separate temporal and spatial sensitivity.

## Exclusions and limits

| Cases / solver | Reason for exclusion |
|:--|:--|
| Gas transport / FDS | No matching porous-solid gas storage and transport equations; gas-production capability does not supply them |
| V-HC-03 / FDS | A matching pure net-flux boundary on a conducting material has not been established: `NET_HEAT_FLUX` is incompatible with `MATL_ID`, and external irradiation also emits thermal radiation |
| V-RAD-01 / FDS | Native internal radiation includes emission; this benchmark specifies absorption-only deposition |
| V-RX-02d | Package-specific smooth depletion limiter; ordinary zero-order burnout is a different equation |
| V-RX-03 and V-RX-05c / FDS | Retained pore gas contributes storage and energy; FDS releases generated gas |
| V-CPL-01 / additional FDS run | FDS samples heat capacity on a 1 K property table, which does not resolve this 0.316 K-wide Gaussian latent-heat peak; a tabulated trial is not published as the same closure |
| V-CPL-01b | Different temperatures, latent heat and gated reaction closure from the retained Gpyro apparent-heat-capacity case |
| V-CPL-02/02b, remaining V-CPL-03 permutations, V-CPL-05 | Full matching of gas retention/transport, reaction gates, deformation and lateral geometry laws has not been established; matching only an endpoint or an isolated TGA rate would not verify the full case |
| V-CPL-04 runaway / bifurcation | Native Arrhenius self-heating alone lacks the benchmark's upper-temperature gate |
| V-MMS and burnout regressions | Custom forcing, mesh/deactivation semantics or coupled gas laws are not represented by a matched native deck |
| V-RAD-01 / ThermaKin2Ds | External radiation is deposited in one element per step (maximum or random absorption), not distributed by Beer–Lambert absorption |
| V-CPL-01 / ThermaKin2Ds | Properties are p₀ + p₁T + p₂Tⁿ; there is no Gaussian latent-heat peak |
| V-RX-02a/c / ThermaKin2Ds | Rates are first order in each of at most two reactants |
| V-GT-02 to V-GT-06, V-GT-08 / ThermaKin2Ds | Gas transport is volume-fraction diffusion; there is no pressure, Darcy or dusty-gas law |

The radiative cooling reference is the thermally thin limit, and the heating
reference is a steady equilibrium, not an exact transient. The external error
panels use these same interpretations as the Pyrolysis.jl panels. Native
Stefan–Boltzmann constants and CSV rounding impose small accuracy floors.
Neither a refinement run nor a small error proves a higher convergence order:
the retained refinement results document the observed change.
