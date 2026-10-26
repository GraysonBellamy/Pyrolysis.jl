# Comparison coverage

External curves are included when native solver inputs represent the same
equations, material laws, boundary conditions and observable. A similar example
name or matching endpoint is insufficient. Comparisons use the shared analytical
reference at each solver's own output coordinates, without fitting simulation
results. Different discretizations remain visible as numerical error.

Comparisons use unmodified FDS 6.11.1 and Gpyro V2 revision
`bff0126c68f2a35c084ef04490ae327185f281a3`. The native refinement studies
include time-only runs at dt/2 and dt/4, spatial refinement at fixed dt/4,
and coupled mesh/time refinement. Native inputs, outputs, executable hashes,
extraction code and equation mappings are under `test/verification/reference/`.
Each case's comparison section links its analytical-error CSV. These studies
provide accuracy measurements. Matching equations does **not** establish
matching discretization, observation operators, or achieved accuracy.

## Comparison basis audit

The [per-pair inventory](assets/comparison_basis.csv) records the actual Julia
mesh used in the published overlay, external intervals and spacing, time
qualification, and requirements for a matched comparison. The
[time-refinement differences](assets/comparison_time_changes.csv) isolate
temporal sensitivity at fixed mesh. These differences have native observable
units; they are not convergence orders or runtime qualifications.

The inventory covers **52 case/code/subcase pairs across 28 case IDs**:
12 have matching plotted spacing, 17 have different plotted meshes, 19 are
spatially homogeneous ODE comparisons, and 4 require a separate material
interface audit. There are 220 retained native accuracy runs across the
three runners (44 native solver/subcase combinations × five configurations).
The Julia study contains 105 solves across 35 case/mesh/variant workloads;
all three tolerance levels reached their expected return codes and the
accepted-step observer counts agreed with the solver statistics.

The comparison basis and its limits are:

| Comparisons | Mesh, observations and qualification |
|:--|:--|
| HC-01 / Gpyro V2 | 501 boundary nodes span 500 intervals at 0.1 mm, matching Julia's plotted mesh. Independent time/space refinements are retained. |
| HC-01 / FDS | Published Julia uses 500 cells versus the reported native 133; native probes are time averaged. A matched comparison requires instantaneous probes (or identical averaging), a retained native grid, and the same mesh in Julia. |
| HC-02a–d | Julia's plotted 320 cells differ from Gpyro's 80 intervals and FDS's 101 cells. Julia's resolution ladder includes 80- and 101-cell meshes; the plots use the finest mesh. FDS uses t_Julia = 2 t_FDS, so its 0.5 s outer step represents 1 s on our axis. A matched comparison requires instantaneous FDS outputs and separate overlays at each native mesh spacing. |
| HC-03 / Gpyro V2 | 321 nodes span 320 intervals at 0.0625 mm, matching Julia's plotted mesh. Independent refinements are retained. The 1200 s observation horizon uses a one-requested-step native output guard; exclude that guard from timed workloads. |
| HC-04a–d / FDS | Both layers have explicit `CELL_SIZE` and `STRETCH_FACTOR` controls. The base grid has 80 cells; Julia's published mesh has 40. A matched comparison requires an 80-cell Julia solve and matching layer boundaries. |
| HC-04c/d / both codes | Equal total count is insufficient: prescribed and represented material interfaces differ. FDS follows physical layer thicknesses; Gpyro samples obstacle composition on boundary nodes. Keep c as a geometry-quantization experiment. For d, compare explicit interface-aligned meshes and report each layer's widths and assigned-material boundary. |
| HC-05 | Julia's plotted 320 cells versus 80 external intervals. Use the existing 80-cell Julia rung. |
| RAD-01 / Gpyro V2 | 401 nodes span 400 intervals at 0.025 mm, matching Julia's plotted mesh. Independent time/space refinements are retained. |
| RX-01/02 and RX-05 | Spatially homogeneous ODEs: one Julia cell, two Gpyro boundary nodes, or native 0D are physically comparable, but not equal work. Match time observations and error targets; do not manufacture a spatial convergence claim. RX-05c retains its native gas-constant error floor. |
| RX-04a/b | Base spatial spacing matches (4/20 Julia cells versus 4/20 native intervals). Slab-mean quadrature must remain fixed during time refinement. Cooling uses a thermally thin approximation; heating has only a steady exact reference. A whole-history error against the steady line is not transient numerical error. |
| CPL-01 route A | The plotted 800-cell Julia mesh is uniform, matching V2's 801 boundary nodes. Every native node contributes to front interpolation. Independent refinements retain the finite Gaussian-peak and sharp-interface-reference distinction. Route B remains excluded. |
| CPL-03 P1/P2/P3 | Base spacing matches 8 Julia cells and 8 external intervals. Time only the same paths and output cadence, excluding unrelated Julia QoI sub-solves. |
| CPL-04 stable branch | Julia's plotted 160 cells versus 80 native intervals; use the existing 80-cell rung. Time the δ = 0.44 solve only, excluding the report's unrelated bifurcation sweep. Establish settling at the common final time. |
| GT-07a/b and GT-08a/b | Base spacing matches (160/80 intervals). Independent time and space refinements isolate their respective sensitivities. Preserve cell-average versus point-value distinctions and native constant differences. |

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
of a matched spatial transport mesh. The V2 RX-05 two-node decks use
the same native implementation but retain only point histories, not a grid
profile; they likewise do not qualify a spatial transport comparison.

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

1. Freeze the mesh and physical observation times. For the reference codes
   run dt, dt/2 and dt/4 with unchanged output cadence, quadrature, solver
   iteration tolerances and equations. Retain accepted/rejected-step or
   substep information when available; explicitly mark unavailable counters.
2. For Julia run the baseline tolerances and two tighter settings (e.g.
   rtol/10 and rtol/100, with matching atol scaling). Also reduce `dtmax`
   when there is a fast event or a tolerance plateau. Keep the mesh fixed.
   Record accepted/rejected steps and RHS/Jacobian/linear-solve counts;
   record min/median/max accepted dt if a step observer is enabled.
3. Compare the same physical quantities at common times. Use instantaneous
   samples, or apply the same averaging window to both codes. Never silently
   select the nearest output or extrapolate a missing endpoint. Bound any
   interpolation error by increasing output density separately.
4. Allocate a temporal-sensitivity allowance of **10% of the chosen total
   error target**. Require both refinement changes to fall within it; inspect
   nonmonotone changes, roundoff and reference-model floors rather than
   inferring an order automatically. This is an observed sensitivity check,
   not a mathematical upper error bound.
5. Then refine space with the qualified time controls held fixed and recheck
   temporal sensitivity on the finer mesh. A dt/4 spatial pair alone does
   not establish time independence at the fine mesh. For steady comparisons
   also check late-time drift (or a longer run); do not interpret approach to
   equilibrium as transient discretization error.

The native runners retain the five configurations `(h,dt)`, `(h,dt/2)`,
`(h,dt/4)`, `(h/2,dt/4)` and `(h/2,dt/2)`. No unobserved step
history is reconstructed from output spacing.

### Runtime at matched accuracy

Publish a work-versus-error curve, with matched-spacing errors as a separate
result. Choose the observation set, normalization and error targets **before**
tuning either solver. A useful initial ladder is normalized maximum error
10⁻², 10⁻³ and 10⁻⁴, supplemented by physical absolute tolerances. Normalize
temperature by the imposed temperature change, concentrations by initial
bulk density, inventory/enthalpy by their common initial scales, and flux
by the prescribed or analytical scale. For a zero invariant use an explicit
nonzero absolute scale. These are proposed study targets, not retroactive
claims that the existing runs meet them. CI regression tolerances need not
be suitable cross-code accuracy budgets.

For each target select the least costly **qualified** mesh/time combination
for each solver; different cell counts and timesteps are allowed in this
performance comparison. Require all selected QoIs to meet the target and
the temporal/output-sampling allowances. A code whose native constant,
closure or output precision imposes a larger floor is "target not attained",
not assigned a misleading speedup. Approximate thermally thin, sharp-front
and steady references need separately verified limits; their model error
cannot be charged as time-integrator error.

Benchmark the selected configurations on the same idle machine and pinned
CPU, serially, with one warm-up and at least five repetitions. Report median
and spread, executable/source hashes, thread count, hardware, output cadence,
actual grid, tolerances and actual work counters. Compare warmed Julia
construction-plus-solve with native internal setup-plus-simulation timers;
list native output costs explicitly. Full process times are a separate
metric and require a full-process Julia measurement too. Do not time plotting,
analytical references, unrelated parameter sweeps or qualification searches.

**No runtime pair is qualified under this policy.** The reports provide
accuracy comparisons. Matched-accuracy runtime ratios require
interval-matched meshes, common observations, per-QoI accuracy qualification,
and measurements of the qualified workloads.

### Reproduce the audit

```sh
julia --project=scripts/verification scripts/verification/export_comparison_meshes.jl \
  test/verification/reference/comparison_meshes.toml
python3 scripts/verification/audit_comparison_basis.py
julia --project=scripts/verification scripts/verification/run_comparison_time_study.jl
python3 -m unittest discover -s scripts/verification -p test_comparison_basis.py
```

`matched/run_comparisons.py` runs all five levels by default. Its `--cases`,
`--levels` and `--destination` options support isolated pilots; a partial
study cannot replace the complete reference manifest. The Gpyro V2 accuracy
runner and `gpyro/v2_migration/run_migration.py` also retain the five levels.
The latter covers HC-01, HC-03, RAD-01, RX-01 and CPL-01 route A with the
pinned V2 executable. Each matched run stores a hashed
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
RX-01 is a spatially homogeneous, one-unknown 0D chain. Stefan front
interpolation uses every node in the full native profiles.

The dataset contains 25 successful native runs (five cases × five
configurations). Inputs, logs, full profiles, extracted curves and hashes are in
`test/verification/reference/gpyro/v2_migration/`.

The [analytical-error table](assets/comparison_v2_migration_accuracy.csv)
records sampled errors at all five configurations. Base maximum profile
errors are 0.0552 K for HC-01, 0.0116 K for HC-03 and 0.0502 K for RAD-01.
RX-01's maximum species-fraction error is 9.98e-5, with no printed mass or
temperature drift. CPL-01's maximum front error at the selected 100, 200 and
299 s samples is 0.1092 mm; the finite Gaussian latent-heat peak and spatial
front interpolation retain an offset from the sharp-interface reference.
These sampled checks do not establish whole-history error bounds.

## Solver coverage

| Verification problem | FDS | Gpyro V2 | Scope |
|:--|:--|:--|:--|
| V-HC-02a–d: convective slab / Biot sweep | Equivalent dimensionless problem | Profiles and histories | Four profiles and front/back histories |
| V-HC-04a–d: layered conduction | Available | Available | Steady profile against the prescribed interface; c/d share native runs, with only Julia's mesh aligned in d |
| V-HC-05: temperature-dependent conductivity | Available | Available | Steady Kirchhoff profile; same affine conductivity law |
| V-RX-02a–c: reaction order | Direct inputs | Direct inputs | Zero-, second-, and two-thirds-order concentration histories |
| V-RX-04a–b: radiative boundaries | Available | Available | Slab-mean cooling and surface approach to radiative/convective equilibrium |
| V-RX-05a–b: unequal heat capacities / reaction heat | Available | Available | Temperature, solid concentrations and sensible enthalpy; FDS requires `ADJUST_H=F` |
| V-CPL-03: volume change | P1/P2/P3 | P1/P2/P3 | Pure solid swelling, shrinking and constant-volume paths only |
| V-CPL-04: thermal feedback | Partial | Partial | Subcritical steady profile at δ = 0.44 only; no claim about the gated runaway branch |

Other retained comparisons include V-HC-01, V-HC-03 (Gpyro V2), V-RAD-01
(Gpyro V2), V-RX-01, the apparent-heat-capacity route of V-CPL-01 (Gpyro V2),
and Gpyro V-GT-07/08 and V-RX-05c runs.

V-HC-04c remains an off-face geometry-resolution test. V-HC-04d places a
Julia mesh face at the same physical 10.2 mm interface; both report prescribed
and represented geometry errors. Reusing the native runs is explicit in the
manifests because the materials, geometry and boundary conditions are identical.

V-GT-07 reports pointwise errors and analytical cell-average errors for
Julia's finite-volume inventories. Gpyro boundary-node outputs
retain point references. Resolution figures compare maximum profile, surface
flux and inventory errors over 50, 150 and 400 s against actual mesh spacing.
The Gpyro `base`/`refined` pair also halves its timestep, so those
two plotted samples do not establish an isolated spatial convergence order.
The retained `time_half`, `time_quarter` and `space_half` runs separate
temporal and spatial sensitivity.

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
| ThermaKin2Ds | No runnable local executable and retained verification run available; local material files alone do not establish comparison data |

The radiative cooling reference is the thermally thin limit, and the heating
reference is a steady equilibrium, not an exact transient. The external error
panels use these same interpretations as the Pyrolysis.jl panels. Native
Stefan–Boltzmann constants and CSV rounding impose small accuracy floors.
Neither a refinement run nor a small error proves a higher convergence order:
the retained refinement results document the observed change.
