# Gpyro V2 accuracy comparisons

These runs compare all seven members of V-GT-07, V-GT-08 and V-RX-05 with
Gpyro V2 and the independent analytical solutions used by Pyrolysis.jl.
They were generated on 2026-10-26. **No timing measurements are retained or
used for this comparison.** The machine was loaded during the runs.

The problem sources are the [Gpyro Verification and Validation Manual v1.0
(15 January 2026)](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/documentation/gpyro_verification_manual.pdf):
§4.4, pp. 30–31, Eqs. (4.7)–(4.11) for diffusion; §4.3, p. 29,
Eqs. (4.5)–(4.6) and Appendix A.2 for Darcy flow; and §3.6, p. 23,
Eq. (3.11) for reaction enthalpy. The adaptations below match our verification
cases; these are new runs, rather than curves digitized from the manual.

## Upstream case files

The actual upstream input decks were inspected before preparing these adapted
runs. Their adjacent `run.sh` and `plot_results.py` files were also reviewed.
All links below identify the same pinned revision as the executable:

- Film diffusion: [`hm_0001/gas_diffusion.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Gas_Transport/Gas_Diffusion/hm_0001/gas_diffusion.data)
  and [`hm_01/gas_diffusion.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Gas_Transport/Gas_Diffusion/hm_01/gas_diffusion.data).
- Prescribed-inflow Darcy flow: [`darcy_mfixed/darcy.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Gas_Transport/darcy_mfixed/darcy.data).
- Reaction enthalpy: [`endothermic_reaction.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Species_Conservation/Reaction_Enthalpy/endothermic_reaction/endothermic_reaction.data),
  [`exothermic_reaction.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Species_Conservation/Reaction_Enthalpy/exothermic_reaction/exothermic_reaction.data)
  and [`sensible_enthalpy.data`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/Verification/Species_Conservation/Reaction_Enthalpy/sensible_enthalpy/sensible_enthalpy.data).

The diffusion decks assign both gases a molecular mass of 32 g/mol; our
equal-mass tracer follows those inputs. The Darcy deck supplies `MDOTPP=10000`
in g/(m² s), matching our 10 kg/(m² s) imposed inflow. The upstream reaction
cases use externally heated slabs and temperature-dependent kinetics, and
their post-processing checks the chemical heat source reconstructed from
reaction rates. Our insulated, constant-rate reaction cases instead allow
independent closed-form temperature and mass checks. The retained decks are
adaptations to the problems described below, not unchanged upstream runs.

## Executable and reproduction

- Source: [Gedeben/GPYRO-V2, revision
  `bff0126c68f2a35c084ef04490ae327185f281a3`](https://github.com/Gedeben/GPYRO-V2/tree/bff0126c68f2a35c084ef04490ae327185f281a3).
- Source archive SHA-256:
  `061e033d05501e65874738641b441d078245e00f3ef5b30d3f8fe19ac7a3648b`.
- GNU Fortran 14.4.0, local `gcc:14` image `9188ac751ca2`, upstream
  `gnu_linux` flags: `-O2 -unroll -frecord-marker=4 -ffree-line-length-none`,
  with the upstream OpenMP flags. Runs set `OMP_NUM_THREADS=1`.
- No solver-source changes. Build metadata identifies the pinned source archive.
  The executable SHA-256 and every retained deck/output hash are in
  [provenance.json](provenance.json).

From the extracted source's `build/linux/gpyro` directory, build with:

```sh
GPYRO_FCOMPL_GNU=gfortran make -f ../../Makefile_gpyro gnu_linux \
  GIT_HASH=bff0126c68f2a35c084ef04490ae327185f281a3 \
  GIT_STAT= GIT_BRANCH=archive GIT_DATE=unspecified
```

Then run the retained decks, both base and refined resolutions, and extract:

```sh
python3 test/verification/reference/gpyro/v2_accuracy/run_accuracy.py \
  /path/to/gpyro /tmp/gpyro-accuracy-comparison \
  --compiler 'GNU Fortran 14.4.0; upstream gnu_linux flags'
```

`--extract-only` rechecks an existing workspace without rerunning the solver.
The script requires successful termination, the pinned revision banner,
finite outputs, the requested final time and the expected transport
coefficients. Gpyro can exit with status zero after a failed solve, so exit
status alone is not accepted. [decks/](decks/) contains the base inputs;
[raw/](raw/) retains native numerical CSVs and the actual inputs for both
resolutions. Native timing files are excluded. Extracted CSVs beside this
README feed report mode; [refined/](refined/) and [accuracy.csv](accuracy.csv)
retain the refinement audit.

Regenerate the Julia report figures without measuring or publishing wall time:

```sh
julia --project=scripts/verification scripts/verification/run_reports.jl \
  --cases=V-GT-07,V-GT-08,V-RX-05 --timings=false
```

## Matching the problems

### V-GT-07: film diffusion

Both gases have molecular mass 32 g/mol, T = 300 K, P = 100000 Pa,
porosity 0.6, D = 2e-5 m²/s and external tracer fraction 0.22.
The two exposed-film coefficients are 1e-4 and 1e-2 kg/(m² s).
Gpyro's Chapman–Enskog expression is evaluated with `EPSOK=100 K`
and `SIGMA=3.588275763111019 Å` for both gases, giving the specified D.
These parameters are calculated from the transport formula, without fitting
simulation output. Temperature is fixed at 300 K (`FIX_DOMAIN_TEMPERATURE`
with `BETA=0`). The slab is 1 m deep; the base run has 161 boundary-centered
nodes (160 intervals) and dt = 0.05 s. The refined run has 321 nodes and
dt = 0.025 s. Profiles are compared at 50, 150 and 400 s.

Gpyro clamps `HM=0` to the single-precision machine epsilon,
1.1920928955078125e-7 kg/(m² s), in
[`INTERPOLATE_BCS`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_bc.f90).
The far-side reservoir therefore contains the initial carrier only, rather
than the exposed tracer mixture. At this distant boundary the largest
computed tracer fraction is 1.66e-15 and the residual exchange is below
2e-22 kg/(m² s). Its effect is negligible over the reported interval.
Leaving the exposed composition at the back would create a second tracer
source and would not represent the sealed-back benchmark.

The influx is reconstructed from the native surface fraction and Robin law,
`h_m (0.22 − Y_surface)`. Inventory uses the native boundary-node control
volumes (trapezoidal weights), porosity, and native gas density. Gpyro's
`RHOGOFT` uses R = 8314 J/(kmol K), while the shared analytical benchmark
and Pyrolysis.jl use 8314.46261815324 J/(kmol K). This small density
and film-Biot difference is retained and is not fitted away.

### V-GT-08: prescribed-inflow Darcy flow

The slab is 0.1 m deep, permeability 5e-10 m², porosity 0.6,
T = 300 K, outlet pressure 100000 Pa, and incoming flux 10 kg/(m² s).
Member a feeds the bottom; member b feeds the top. The initial pressure is
uniform at the outlet value. Gpyro uses 81 nodes (80 intervals), dt = 0.005 s;
the refined run uses 161 nodes and dt = 0.0025 s. The shared steady analytical
profile is checked at 1 s, and signed fluxes at 0.5 and 1 s.

Gpyro uses its Chapman–Enskog `D12` as kinematic viscosity in the Darcy
coefficient. Set M = 29 g/mol, `EPSOK=100 K`, and
`SIGMA=3.965473809316032 Å` so that `D12 × P = μRT/M` for the benchmark's
μ = 2e-5 Pa s and gas constant. This directly matches the steady pressure
equation. Native gas-density rounding affects transient storage, which is
not compared with a steady reference. Temperature is held at 300 K with
`FIX_DOMAIN_TEMPERATURE` and `BETA=0`.

Native gauge pressure is shifted by 100000 Pa. Native depth becomes
`z = 0.1 − depth` in Julia coordinates. Gpyro takes incoming `MDOTPP=10000`
in g/(m² s). The top outward flux is minus the native z flux at the top;
the bottom outward flux is native integrated MLR minus the top outward flux.
All are divided by 1000 for kg/(m² s). This uses actual boundary quantities,
not the last interior face. The t = 0 flux dump is omitted because Gpyro has
not yet initialized the imposed inflow in that dump.

### V-RX-05: insulated reaction

The cases use A = 0.2 s⁻¹, E = 0, first-order kinetics, T_initial = 1000 K,
initial condensed concentration 1000 kg/m³, and c_A/c_B/c_G =
1200/1900/1000 J/(kg K). The reaction heats are +1e6, −1e5 and 0 J/kg.
Members a/b retain all product as solid B; member c produces 0.4 B + 0.6 G.

Two boundary nodes span a uniform 1 mm insulated slab with no mesh deformation,
zero permeability and zero prescribed boundary mass flux. Base dt is 0.001 s;
refined dt is 0.0005 s. Gpyro's `ZEROD` mode prescribes temperature, so the
insulated energy calculation uses this spatially uniform slab instead.
`R0(B)=400` for member c encodes the 0.4 condensed yield in Gpyro. Both
species' true solid densities are 1000/0.7 kg/m³, matching the evolving pore
volume in our fixed mesh. All generated gas is retained by the pressure solve.

A seed pressure of 1e-6 Pa avoids division by zero in native diffusivity at
vacuum; its initial gas mass is about 1.15e-15 kg/m². It has negligible
thermal capacity in a/b. The zero-heat case has no temperature change,
so Gpyro's default thermal-equilibrium treatment gives the same temperature
invariant without requiring a gas-storage heat term.

Native temperatures are converted from °C to K. `CI` and `CJ` supply the
condensed and gaseous bulk concentrations. Sensible enthalpy is reconstructed
as `L (c_A C_A + c_B C_B + c_G C_G) (T − 298.15)`, using the common datum;
it is not Gpyro's native condensed-only specific-enthalpy output.

Member c reveals a persistent native mass discrepancy: approximately
3.271e-5 of initial total mass. The
[pressure equation](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_gas.f90)
uses the `GPYRO_VARS` gas constant (declared `8.314462`, represented as
8.314461708068848), but
[`RHOGOFT`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_funcs.f90)
uses 8.314 in molar units. The reported retained-gas concentration agrees with
`(1000 − C_A − C_B) × 8.314461708068848 / 8.314`
within 3.3e-6 kg/m³ at the saved times. Halving dt leaves this total-mass offset
unchanged. The comparison preserves that discrepancy; no rescaling of the
reported gas mass is applied.

## Accuracy and resolution

Errors use the same analytical functions as the Julia cases, evaluated at
each code's own output locations and times; no interpolation between codes
is needed. This compares the stated numerical configurations, not equal
error budgets. A Gpyro boundary node is not a Julia cell center.
Native point outputs have seven decimal places (temperature resolution
1e-7 K; boundary-flux resolution 1e-10 kg/(m² s)); zero errors mean agreement
at that output precision. Logarithmic report panels floor zero errors at
1e-12 for display.
| Case | Error metric | Pyrolysis.jl | Gpyro base | Gpyro refined |
| --- | --- | --- | --- | --- |
| V-GT-07a | max tracer-fraction error / 0.22 | 0.0001283 | 0.00033449 | 9.5296e-05 |
| V-GT-07b | max tracer-fraction error / 0.22 | 0.001282 | 0.00051712 | 0.00013634 |
| V-GT-08a | max P² error / P² rise | 0.00081294 | 0.0016411 | 0.00081933 |
| V-GT-08b | max P² error / P² rise | 0.00081294 | 0.0016314 | 0.00081692 |
| V-RX-05a | max temperature error (K) | 3.1624e-08 | 0.013802 | 0.0069006 |
| V-RX-05b | max temperature error (K) | 5.4374e-09 | 0.0013802 | 0.00069002 |
| V-RX-05c | max species concentration error / initial concentration | 2.9656e-12 | 4.5978e-05 | 3.6955e-05 |

Julia uses its finest report resolution: 160 cells for diffusion, 80 cells
for Darcy, and one cell for the uniform reaction. Gpyro refinement improves
the diffusion and Darcy profiles and halves the a/b reaction-temperature
errors. The retained-gas mass offset persists with refinement; its zero
temperature error is only an output-precision statement.
