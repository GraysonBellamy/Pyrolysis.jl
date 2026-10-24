# 18. Verification and Validation

This chapter sets out the verification and validation (V&V) methodology for
Pyrolysis.jl. Following Roache's standard distinction, *verification* asks
whether the discrete equations are solved correctly (a mathematics question,
internal to the code) and *validation* asks whether the right equations are
solved (a physics question, answered against experiment and trusted reference
codes). We catalogue the exact-solution and manufactured-solution checks of the
verification suite, the Jacobian-against-reference tests, the conservation
ledgers, the spatial, temporal and tolerance convergence studies, the
experimental and code-to-code validation cases, the regime-validity checks that
tell the user when a modelling approximation is being stressed, and the known
discrepancies that still require attention. Every tolerance, sign convention and
discrete formula stated here has been traced to the implementing source or test;
where a formula appears it is the one the code evaluates.

All symbols conform to the nomenclature of §2; the binding conventions are reused
throughout: (i) a face flux is the transport per unit area along the face normal,
which points from the owner cell to the neighbour ($+z$ on interior faces and on
`:top`, out of the domain on `:bottom`); (ii) the discrete divergence
$(1/V_i)\sum_fs_{i,f}A_fF_f$ is positive for net outflow; (iii) heat of reaction
$h>0$ is endothermic, with the volumetric source $Q_\mathrm{rxn}=-\sum_rh_rr_r$.

---

## 18.1 The V&V philosophy

### 18.1.1 Verification versus validation

A pyrolysis solver is a chain of approximations: a continuum PDE system (§3) is
closed by constitutive laws (§4–§10), projected onto a finite-volume mesh (§13),
advanced by a stiff integrator with a colored sparse Jacobian (§15), and read out
through post-processing (§16). Errors can enter at every link, and they have
qualitatively different remedies. We therefore separate the two questions
explicitly.

- **Verification** — *Are we solving the equations right?* This is purely
  internal. It compares the code's discrete output against (a) closed-form
  solutions of the *same* continuous equations, or manufactured ones, (b) the
  dense forward-mode reference of its own Jacobian, and (c) the conservation laws
  the discretization is designed to satisfy. A verification failure is a *bug* (a
  wrong sign, a missing factor, a flux scattered with one sign, a Jacobian entry
  outside the pattern), not a modelling deficiency. Verification can, in
  principle, be made exact: the component masses are states, so a conservative
  scheme conserves mass to the integrator's tolerance; the colored Jacobian equals
  the dense reference per entry to round-off; and cell volumes recomputed from the
  node positions satisfy the geometric conservation law as an identity.

- **Validation** — *Are we solving the right equations?* This is external. It
  compares the code's predictions against measured cone-calorimeter mass-loss and
  surface-temperature histories, against thermogravimetric (TGA) curves, and
  against the output of established codes (ThermaKin2Ds, Gpyro, the FDS solid
  phase). A validation discrepancy can be a modelling choice (one-temperature gas,
  gray-body radiation, single-reactant kinetics), a parameter-calibration issue,
  or genuinely missing physics. Validation is never exact: it is bounded by
  measurement uncertainty and by the irreducible approximations of a 1D
  condensed-phase model.

The two are kept methodologically distinct because conflating them is the most
common way to draw a wrong conclusion. A model that agrees beautifully with a TGA
curve after the activation energy has been tuned to it has been *calibrated*, not
validated; conversely, a perfectly verified code can still disagree with
experiment because its physics is incomplete. Throughout this chapter we keep the
ledger checks (verification) separate from the experimental comparisons
(validation), and we report the regime checks (§18.7) that quantify how hard a
given run leans on an approximation.

### 18.1.2 The verification hierarchy in Pyrolysis.jl

Pyrolysis.jl layers its verification from cheap, fundamental checks to experimental
validation:

1. **Kernel-level and exact-solution checks** (§18.2): the kernels tested in
   isolation (`test/kernels/`, `test/physics/`), and the cases of the verification
   suite compared against closed-form or manufactured solutions — conduction,
   kinetics, in-depth absorption, gas transport, coupled fronts, the geometric
   conservation law, burnout through a deactivation event.
2. **Jacobian verification** (§18.3): the colored Jacobian compared per entry with
   the dense forward-mode reference of the whole residual, in the test suite. There
   is no in-solve Jacobian check; a user's problem is checked with the same
   functions from a script.
3. **Conservation ledgers** (§18.4): the mass and energy ledgers of §16, which any
   solve can carry (`callback = ConservationRecorder()`, then
   `conservation_audit(sol, rec)`), driven by the module tests and by the burnout
   cases' own per-step ledgers.
4. **Convergence studies** (§18.5): spatial, temporal and tolerance refinement,
   confirming that the discretization error shrinks at the designed rate.
5. **Code- and experiment-level validation** (§18.6): cone calorimeter, Douglas fir,
   pressure-treated wood, moisture transport, against measurements and other codes.

The first three run in the test suite, and the verification suite asserts
convergence orders as well (§18.1.3); the convergence studies of §18.5 and the
validation cases are scripted examples in `examples/`.

### 18.1.3 The verification suite: test mode and report mode

The exact-solution and manufactured-solution cases form a verification suite that
runs in two modes from a single case definition. The suite has to stay cheap enough
for `Pkg.test()`, while publishable evidence (fine meshes, resolution ladders,
figures) is far too expensive for it and needs plotting dependencies the test
environment should not carry.

- **Test (assert) mode** runs inside `Pkg.test()`. `test/runtests.jl` runs the
  module suites directory by directory, each file in its own module, and then
  `test/verification/runcases.jl`, the suite's entry point. It first runs the
  self-tests of the exact-solution library (`test/verification/exact/selftests.jl`),
  then every case: a discarded warm-up solve, so the wall-time budget measures the
  case and not compilation; the timed solve at the case's CI resolution, asserting
  `sol.retcode == case.expected_retcode` and each quantity of interest (QoI) below
  its pinned tolerance; where the case defines one, a spatial-convergence spot check
  asserting every observed order between consecutive resolutions above a minimum
  and their mean above a second one; and the wall time below five times the case's
  budget. It plots nothing and writes no files. Three environment variables steer
  it: `PYROLYSIS_VERIFICATION_PRINT_QOIS=1` prints every observed QoI and
  convergence rate (how tolerances are re-pinned, at about three times the
  observed error); `PYROLYSIS_TEST_FAST=1` skips the convergence spot checks;
  `PYROLYSIS_ENERGY_STATE=enthalpy` runs every solve of the suite in the enthalpy
  form (a run of `runtests.jl` in that form runs only the verification cases, since
  the module suites cover both energy forms themselves). `PYROLYSIS_TEST_GROUP`
  selects a part of the tests (`jacobian`, `burnout`, `other`, `domains` for the
  verification cases, or `all`).
- **Continuous integration** (`.github/workflows/tests.yml`) runs a fast matrix on
  every pull request (the three module groups and the verification cases side by
  side, in the temperature form, without the spot checks) and a full matrix weekly
  (`Pkg.test()` as a user runs it, with the spot checks, and the verification cases
  in the enthalpy form with their spot checks).
- **Report mode** runs in the isolated environment `scripts/verification/` (its
  own `Project.toml` with CairoMakie, so the test manifest stays free of plotting
  dependencies):

  ```
  julia --project=scripts/verification scripts/verification/run_reports.jl [--cases=… --outdir=… --chapters=…]
  ```

  `run_reports.jl` runs each case over its report resolution ladder with its report
  solver settings and writes the artifact set under `<outdir>/assets/<case id>/`:
  overlay figures (solver values as markers over densely sampled exact curves),
  pointwise-error figures that show *where* the discretization struggles, a log–log
  convergence figure with the fitted observed order and an expected-order reference
  slope, and CSVs of every curve, the resolution ladder, the QoIs at the finest
  rung, the ladder QoI errors and the solver configuration (`report_harness.jl`,
  with all plotting in `render.jl`). `chapters.jl` then writes one Markdown chapter
  per case family (problem statement, figures, error and order tables) and an index.
  Report mode does not assert the QoIs (a failed solve still aborts the run): coarse
  ladder rungs may exceed the pinned tolerances, and the tables record
  within-tolerance yes/no instead. The default output directory is the committed
  `docs/guide/verification/`; the documentation build publishes it, and the manual
  `verification-report` workflow regenerates it and opens a pull request with the
  diff.

**The shared case definition.** `test/verification/casedefs.jl` declares the types
both modes consume and depends on nothing from `Test`, so the report runner can
include it from its isolated environment. Each file in `test/verification/cases/`
defines a `VerificationCase` with

- identity and provenance: `id`, `title`, physics `tags`, `references`;
- `build(n_cells) -> PyrolysisProblem`, optional `solve_kwargs`, a `solve_fn` when
  the case solves through its own path (the manufactured-solution cases use
  `solve_mms`, the burnout cases a solve with their step recorder), and
  `expected_retcode` (`:Success` by default, `:Terminated` when the stop is the
  answer);
- `qois(sol) -> Vector{QoIResult}`: each QoI's value, exact value, error and
  tolerance, with a `comment` recording the tolerance's provenance;
- `ci::CIConfig`: the CI resolution, the wall-time budget and an optional
  `ConvergenceSpec` (resolutions, error functional, expected order, minimum
  observed order and minimum mean order);
- `report::NamedTuple`: the report-mode keys — the resolution ladder
  `n_cells_ladder`, `curves(sol) -> Vector{ReportCurve}` (plain data: `CurveSeries`
  samples, each with an optional `exact` callable on the same abscissa from which
  the pointwise error is derived), `error_norms`, `expected_order`, report
  `solve_kwargs`, and optional temporal or parameter sweeps. A case with
  `report = (;)` is skipped in report mode.

Every solve a case makes goes through `suite_solve`, which applies the energy state
of `PYROLYSIS_ENERGY_STATE`, so all of a case's quantities come from one energy form.
The case file's leading comment banner is the problem statement the report chapter
prints, so the published chapter cannot drift from what the test runs. Reference
solutions are computed in the run at high precision (series summed to machine
tolerance, transcendental roots by bisection), never hard-coded; they live in the
shared library `test/verification/exact/` (series, roots, the manufactured-solution
driver), which carries its own self-tests.

**Case catalogue.** Each case family has one file in `test/verification/cases/`;
lettered members (V-HC-02a–d, V-HC-04a–c, V-RX-02a–d, V-RX-04a/b, V-CPL-01b,
V-CPL-02b) share the family's file and its report chapter. Most families compare
against an exact or manufactured solution; the burnout cases V-BO-01, 02, 03 and
05 to 08 check mass budgets, ledgers and a pinned reference run instead.

| Family | Cases | Reference solution | What it verifies |
|---|---|---|---|
| Heat conduction | V-HC-01 – V-HC-05 | Carslaw–Jaeger erfc and eigenfunction series, resistance network, Kirchhoff transform | the convective surface Newton (V-HC-01; Biot sweep V-HC-02a–d), the constant-flux slab to its back-face plateau (V-HC-03), the harmonic-mean face conductivity at a layer jump on and off a face (V-HC-04a–c), $k(T)$ (V-HC-05) |
| Kinetics and surface radiation | V-RX-01 – V-RX-04 | Bateman chain, Khawam–Flanagan isothermal integrals, parallel-reaction splits, implicit radiative-cooling law | multi-step coupling, reaction order and the depletion limiter (V-RX-02d keeps it on as a guard), two-product yields and heat superposition, `RadiativeBC` through the surface Newton |
| In-depth radiation | V-RAD-01 | Lautenberger eigenfunction series | Beer–Lambert absorption with conduction and a convective face, including the sub-surface temperature peak |
| Gas transport | V-GT-01, V-GT-02 | Gaussian plane source; steady compressible Darcy $P^2$ profile | the Fickian flux and mass conservation; the pressure closure, Darcy mobility, upwinding and the pressure outflow condition |
| Coupled physics | V-CPL-01 – V-CPL-05 | Neumann two-phase Stefan solution, Landau ablation wave, exact shrink/swell endpoints, Frank-Kamenetskii slab, 0D TGA curve | a latent-heat front (apparent heat capacity, V-CPL-01; gated reaction, V-CPL-01b), surface regression on a moving mesh (V-CPL-02, 02b), volume change and the lateral-shrinkage law, the reaction–conduction bifurcation, the boundary and volume coupling of the kinetics |
| Manufactured solutions | V-MMS-01 – V-MMS-04 | the manufactured-solution driver of §18.2 | design order of the coupled residual on a fixed mesh in space and time, free-stream preservation under prescribed motion, order under prescribed and reaction-driven motion |
| Burnout | V-BO-01 – V-BO-08 | exact exponential thickness (V-BO-04); otherwise mass budgets, ledgers and a pinned reference run | deactivation through complete burnout on a moving mesh under in-depth absorption (V-BO-01) and on a fixed mesh under surface absorption (V-BO-02); the transmitted-flux series (V-BO-03); the interior-void stop (V-BO-05, `expected_retcode = :Terminated`); a graded mesh (V-BO-06), a layer of mixed composition (V-BO-07) and a thin layer on a deep substrate (V-BO-08) |

The burnout cases share V-BO-01's deck (`V-BO-01.jl` defines it, its step
recorder and its QoIs): a 5.9 mm PMMA slab on 12 mm of Kaowool under 50 kW m⁻², the
MaCFP PMMA gasification case with the Gpyro V2 property fits, solved at 60 PMMA
cells (0.098 mm). V-BO-01 requires complete burnout, the exact mass of the PMMA
layer lost, the peak gas-generation rate (32.39 g m⁻² s⁻¹ within 3 %), its time
(335 s within 3 %) and the burnout time (363.5 s within 1 %) of a reference run of
the deck, a mass ledger closed within $3\times10^{-5}$ (§18.4.1), and no accepted state
whose surface balance did not converge or sits on the trial clamp; the number of
deactivation invocations is reported, not asserted. V-BO-02 runs the deck on a
fixed mesh under surface absorption, where only the top cell deactivates, one cell
at a time (60 events), with the live-cell temperature below 1000 K. V-BO-03 makes
the substrate transparent and requires the transmitted flux
(`extras.radiation_transmitted`) to exceed 1 % of the entering flux. V-BO-05 runs the
deck on a fixed mesh under in-depth absorption, where the sub-surface temperature
maximum empties the cells beneath a surface skin; it must stop with outcome
`:interior` and return code `:Terminated` at 259.6 s, before any deactivation
(§14.4). V-BO-06 grades the PMMA cells toward the surface (each 1/1.02 of the one
beneath at 60 cells) under V-BO-01's references; V-BO-07 mixes in a second solid of
half the density and a tenth of the pre-exponential factor, graded from none at
the substrate to half at the surface, against its own reference run at 240 PMMA
cells; V-BO-08 burns a 1 mm layer off the 12 mm substrate and requires the correction
velocity to be exactly zero at every Kaowool node at every saved time.

Kernel-level checks sit beside the suite in the module tests: the Beer–Lambert beam's
energy telescoping and the incident intensity under Kirchhoff's law and an explicit
absorptivity (`test/kernels/test_radiation.jl`); the depletion limiter, the
constant-heating-rate peak temperature and the Kissinger slope
(`test/kernels/test_kinetics.jl`); the gas-transport fluxes, the vent and the
upwind side at zero velocity (`test/kernels/test_gas_transport.jl`); the relative
fluxes' free stream, isothermal advection and convergence order
(`test/kernels/test_relative_flux.jl`); and the closures of `test/physics/`. Two
conventions matter when a case transcribes literature parameters: $h>0$ is
endothermic here, the opposite of ThermaKin's published sign (Gpyro and FDS use the
same sign as Pyrolysis.jl, §6.9.2); and most exact-kinetics cases disable the
depletion limiter (`DepletionLimiter(enabled = false)`) and leave the temperature
gates open (§18.2.2), while V-RX-02d keeps the limiter on as a guard.

**Cross-code overlays.** Where a case has an analogue in another code, that code's
numerical output is committed under `test/verification/reference/<code>/` with a
provenance README: FDS 6.11.0 and Gpyro 0.8200 runs of the same problems, a Gpyro V2
run of the V-BO-01 deck, and the MaCFP-2 PMMA community submissions (which no case
reads). They are read only inside a case's report `curves`; assert mode never
touches them. Overlays are drawn for V-HC-01 (FDS and Gpyro), V-HC-02 (FDS), V-HC-03
(Gpyro; FDS cannot run a pure-Neumann slab), V-RAD-01 (Gpyro), V-RX-01 (FDS and
Gpyro), V-RX-02 (FDS), V-CPL-01 (Gpyro), V-CPL-03 (FDS) and V-BO-01 (the Gpyro V2
gas-generation rate). A series marked `crosscode = true` is drawn as a muted
context layer in a separate `_vs_codes` variant of each affected figure, so the base
figure keeps the solver's own error on its own scale; for exact-solution cases the
other codes' *errors* against the same exact solution are compared.
`reference/timings.csv` feeds the wall-time comparison: warmed Julia model
construction and solve versus FDS/Gpyro internal total elapsed time, including
setup and simulation. Full process times appear separately in supporting
tables. Output workloads, resolutions and time steps differ, so these are
indicative costs rather than pure solver speedups. The reference README
defines the timer boundaries and provenance. ThermaKin2Ds is not represented:
its distribution is a Windows console binary.

---

## 18.2 Method of manufactured solutions and analytical checks

The most rigorous verification tool is the *method of manufactured solutions*
(MMS): one chooses analytic fields, substitutes them into the governing equations
to obtain the exact forcing, supplies that forcing to the discrete solver, and
confirms that the numerical solution converges to the chosen fields at the designed
order. The suite's driver (`test/verification/exact/mms.jl`) provides `MMSFields`,
exact-source generators that evaluate every coefficient through the package's own
property and kinetics functions (mixing rules, floors, gates and limiter included)
and take the spatial derivatives of the fields by nested ForwardDiff, manufactured
Neumann boundary conditions (`HeatFluxBC` and `MassFluxBC` with time-dependent
fluxes), and `solve_mms`. `solve_mms` passes the sources to `solve` through its
`sources` keyword as `VolumetricSources`: the species sources on the mass rows as
$S_\xi V$, the temperature-rate source on the energy row (times the heat-capacity
sum in the enthalpy form), both functions of the cell-centre position and the time,
which the colored Jacobian differentiates with the rest of the residual, so the
Jacobian stays exact with the sources active (§13.5.3). The cases are V-MMS-01
(the coupled system on a fixed mesh), V-MMS-02 (free-stream preservation under
prescribed motion), V-MMS-03 (order under prescribed motion) and V-MMS-04 (order
under reaction-driven motion).

Mesh motion enters in two ways. Under **prescribed motion**, `solve_mms` passes the
mesh velocity to `solve` through its `w_override` keyword: the mesh moves over a
solid at rest and every inventory crosses the faces at the mesh-relative velocity
(the relative fluxes of §11.6 with the prescribed correction, §11.7.8), so the
manufactured sources are the fixed-frame ones. Under **reaction-driven motion**
(`MMSMotion`, V-MMS-04) the mesh is the solid and every inventory rides it, so each
source is assembled in the material-following frame,

$$
S_j=\partial_t\xi^*_j+w^*\partial_z\xi^*_j+[\partial_zJ^*_j]_\text{gas}-\dot\omega^*_j+\xi^*_j\theta^*,
$$

with $w^*$ the integral of the manufactured strain rate $\theta^*$; the temperature
source is assembled in the same frame. Because the boundary values a condition fixes
enter through a one-sided difference over the half-cell (§13.12), the driver imposes
fluxes rather than face temperatures, so the boundary does not cap the observed order.

> **What MMS cannot see: wrong-equation defects.** MMS verifies that the
> discretization solves the code's *own* target equations at the designed order —
> the manufactured source is built from those same equations, so a term the
> formulation wrongly includes (or omits) is compensated in the source and is
> invisible at every order and every resolution. A spurious advection term on a
> mesh that rides the solid, for example, would pass every MMS case untouched. Only
> checks with an *independent physical identity* that the code does not define can
> catch a wrong-equation defect: an exact solution (the Landau travelling wave of
> V-CPL-02, the exact thickness of V-BO-04), a conservation endpoint (the
> shrink/swell endpoints and the sealed column of V-CPL-03), or the term-by-term
> construction of the energy row (§18.4.3). This is why the suite pairs MMS
> (discretization order) with exact-solution and conservation anchors (formulation
> correctness): the two families have disjoint blind spots.

Pyrolysis.jl also checks each physics kernel in a parameter regime where an exact
solution of its equation is already known and no manufactured forcing is needed.
These are the building blocks; each anchors a specific term of the residual.

### 18.2.1 Pure conduction

With gas transport, reaction and radiation absent, the energy equation of §7
collapses to the linear heat equation

$$
\rho c_p^\mathrm{eff}\,\frac{\partial T}{\partial t}=-\nabla\!\cdot\mathbf q_\mathrm{cond},
\qquad\mathbf q_\mathrm{cond}=-k\,\nabla T .
$$

For constant $\rho c_p^\mathrm{eff}$ and $k$ on a fixed mesh this is
$\partial_tT=\alpha_\mathrm{th}\,\partial_{zz}T$ with thermal diffusivity
$\alpha_\mathrm{th}=k/\rho c_p^\mathrm{eff}$, which admits a wealth of closed-form
solutions: the error-function response of a semi-infinite solid, with its
convective-surface (Biot-coupled) variant (V-HC-01), the eigenfunction series of a
finite slab with a convective front face (V-HC-02a–d over Bi = 0.1 to 100) or a
constant flux (V-HC-03, with its back-face plateau $\Delta T\to q''L/2k$ and the exact
mean-temperature rise $q''/\rho cL$), and the decaying-cosine eigenmode

$$
T(z,t)=T_\infty+(T_0-T_\infty)\cos\!\Big(\frac{\pi z}{2L}\Big)
\exp\!\Big(-\alpha_\mathrm{th}\frac{\pi^2}{4L^2}t\Big)
$$

on $[0,L]$ with an insulated bottom (`AdiabaticBC`) and a fixed top
(`DirichletThermalBC`). The single decay constant $\alpha_\mathrm{th}\pi^2/4L^2$
makes this an exacting check on the whole conduction pathway: the face
conductivity, the scatter, the heat-capacity sum and the time integration must all
be correct for the numerical decay rate to match.

The discrete face flux verified here is (§13.7.1)

$$
q_f=-k_f\,\frac{T_n-T_o}{d_o+d_n},\qquad
k_f=\frac{(d_o+d_n)\,k_ok_n}{d_ok_n+d_nk_o+\eta(d_o+d_n)},
$$

the guarded distance-weighted harmonic mean. On a *uniform* mesh
$d_o=d_n=\Delta z/2$, and for constant $k$ the scheme is the textbook second-order central
difference; V-HC-01 observes spatial orders 1.991 and 1.998 between 40, 80 and 160
cells. The harmonic mean (rather than the arithmetic one) is verified by the
two-layer slab V-HC-04: with the conductivity jump on a face, the steady cell
temperatures and the boundary fluxes match the series-resistance network to
integrator precision for ratios 10:1 and 100:1; with the interface off a face
(V-HC-04c) the numerics match the network of the face-quantized interface to
$10^{-11}$, and the whole error against the true interface is that quantization
(§13.12). V-HC-05 measures the genuine discretization error of a $k(T)$ slab against
the Kirchhoff-transform solution, with a convergence check.

> **Implementation note.** `ConductionKernel` (`src/Kernels/conduction.jl`) writes
> the flux on the interior faces between two live cells and skips every other face;
> the boundary faces belong to the thermal boundary kernel. `assemble_rows!`
> (`src/Assembly/residual.jl`) scatters $\mp q_fA_f$ into the owner and neighbour
> energy rows, and the temperature row divides by the floored heat-capacity sum
> (§13.3.4).

### 18.2.2 Single-reaction TGA (0D-idealized cell)

A thermogravimetric experiment heats a small, thermally thin sample at a constant
rate $\beta=dT/dt$ and records the residual mass fraction. In the limit of a single,
well-mixed cell with a uniform temperature ramp, all transport vanishes and the
species balance for a one-step first-order decomposition becomes the 0D ODE

$$
\frac{d\xi_1}{dt}=-A_1\exp\!\Big(-\frac{E_1}{R_gT}\Big)\,\xi_1,\qquad T(t)=T_0+\beta t .
$$

Changing the independent variable to temperature and integrating gives the classic
non-isothermal conversion law

$$
\frac{\xi_1(T)}{\xi_1(0)}=\exp\!\left[-\frac{A_1}{\beta}\int_{T_0}^{T}\exp\!\Big(-\frac{E_1}{R_g\,T'}\Big)\,dT'\right],
$$

whose integral is the Arrhenius temperature integral. Its inflection (maximum
reaction rate) sits at the peak temperature $T_p$ fixed by

$$
\frac{E_1\,\beta}{R_g\,T_p^2}=A_1\exp\!\Big(-\frac{E_1}{R_g\,T_p}\Big),
$$

the FDS-style recipe for recovering $(A_1,E_1)$ from a measured $T_p$ and peak rate
(§6.9, §18.6.2). The module tests drive a 0.1 mm cell between two faces whose
temperature follows the ramp: the mass-loss rate peaks within 2 K of the exact $T_p$
at 10 K/min, and the Kissinger slope over 5, 10, 20 and 40 K/min recovers $E_1$ within
5 % (`test/kernels/test_kinetics.jl`). V-CPL-05 makes the chain exact: the 0D curve
against the temperature integral computed to high precision, then 1D slabs under the
same ramp on both faces, whose volume-mean mass must converge to the 0D curve at
second order in the thickness (the interior lags the faces by the quasi-steady
parabola $\Delta T\approx\beta(L/2)^2/2\alpha$).

Two implementation details must be accounted for when the check is made
*quantitatively* exact:

- **Mass-concentration basis.** The rate law is evaluated on the mass concentration
  $\xi_j$ (kg m⁻³), not a mole fraction or a normalized conversion; for a first-order
  single reactant it is linear in $\xi_1$, and the solution above holds with $\xi_1(0)$
  the initial concentration. Published constants for a normalized mass convert as
  $A_\text{ours}=k_\text{pub}\,\xi_0^{1-n}$ (the header of
  `test/verification/exact/kinetics.jl`).
- **Smooth gates and the depletion limiter perturb the tails.** The temperature gates
  are smooth tanh ramps of width `gate_width` (default 1 K) that leak outside the
  nominal window (§6.2), and the one-sided depletion limiter multiplies each
  reactant's rate by $g(\xi/\xi_\mathrm{th})$, a cubic onset that is zero with zero
  slope at $\xi=0$ joined to $\tanh$ (default $\xi_\mathrm{th}=1$ kg m⁻³), with the
  reaction-order power taken on a finite base (§6.3). Far above the threshold the
  factor is one to within $2e^{-2\xi/\xi_\mathrm{th}}$ (below $5\times10^{-9}$ at
  $10\,\xi_\mathrm{th}$). To compare against the bare Arrhenius solution one either
  leaves the gates fully open and keeps $\xi_1$ well above the threshold, or disables
  the limiter with `DepletionLimiter(enabled = false)`; the exact-kinetics cases do
  the latter (V-RX-01, V-RX-02a–c, V-RX-03, V-CPL-01b, V-CPL-03, V-CPL-05, V-GT-02),
  and V-RX-02d keeps the limiter
  on to require $0\le\xi\le\xi_\mathrm{th}$ at twice the zero-order burnout time.

### 18.2.3 Beer–Lambert transmitted fraction

The in-depth radiation model (§8.3) carries the intensity of every face as an
algebraic unknown (§13.5.2): at the face the beam enters,
$0=I_e-a\,q_\mathrm{inc}(t)$; at every other face,
$0=I_\text{down}-I_\text{up}e^{-\tau_i}$ with $\tau_i=\alpha_i\Delta z_i$ the cell's
optical thickness; and each cell absorbs $I_\text{up}A_\text{up}-I_\text{down}A_\text{down}$,
the exact cell average of the attenuated beam,

$$
q'''_i\,\Delta z_i=I_\text{up}\bigl(1-e^{-\tau_i}\bigr)\quad\text{(planar)} .
$$

Two exact energy identities verify the kernel:

- **Total absorbed flux.** The absorbed energy telescopes,

  $$
  \sum_iq'''_i\,\Delta z_i=a\,q_\mathrm{inc}\big(1-e^{-\tau_\mathrm{total}}\big),
  \qquad\tau_\mathrm{total}=\sum_i\alpha_i\Delta z_i,
  $$

  exact at the discrete level and independent of the cell count, because each cell
  deposits $I_\text{up}(1-e^{-\tau_i})$ and passes $I_\text{up}e^{-\tau_i}$ on. Here $a$
  is the effective absorptivity of the exposed cell (Kirchhoff's law unless a
  component sets its own) or the condition's explicit absorptivity (§12.3.5). The
  solution reports it as `extras.volumetric_absorption`.
- **Transmitted fraction.** The intensity leaving the last face is
  $a\,q_\mathrm{inc}e^{-\tau_\mathrm{total}}$, reported as `extras.radiation_transmitted`.
  Absorbed plus transmitted equals $a\,q_\mathrm{inc}$ to round-off; the remaining
  $(1-a)\,q_\mathrm{inc}$ is reflected at the surface and never enters the medium.

`test/kernels/test_radiation.jl` ("the beam's energy telescopes, Kirchhoff and
explicit absorptivity") asserts at a consistent state that the entering intensity is
$a\,q_\mathrm{inc}$ to $10^{-14}$, that every cell absorbs $I_\text{up}(1-e^{-\tau_i})$ and
the last face carries the transmitted intensity to $10^{-12}$, and that absorbed plus
transmitted equals the entering intensity to $10^{-12}$, under both absorptivity
paths. A transparent cell ($\alpha_i=0$) absorbs nothing and passes the whole
intensity; in the opposite limit $\tau_1\gg1$ the exposed cell absorbs the entering
intensity $a\,q_\mathrm{inc}$ entirely, the volumetric counterpart of
`SURFACE_ABSORPTION`, which applies the same $a\,q_\mathrm{inc}$ at the face (§8.2).
V-RAD-01 couples the kernel with conduction and a convective face against
Lautenberger's eigenfunction solution, and V-BO-03 exercises the transmitted series
through a burnout.

The incident flux is sampled at every residual evaluation:
`test/kernels/test_boundary_conditions.jl` ("a time-varying incident flux is sampled
at every evaluation") solves under $q_0(1+\sin2\pi t)$ and under $q_0$ and requires
the absorbed intensity at $t=0.25$ s to be twice the constant run's, to $10^{-4}$
(a flux frozen at its initial value would make the two runs identical, since
$\sin 0=\sin2\pi=0$).

### 18.2.4 Geometric Conservation Law to machine precision

The single most decisive moving-mesh verification is the Geometric Conservation Law
(GCL, §11.5). A moving-mesh scheme that does not reproduce the swept cell volume
exactly injects *phantom* mass and energy sources proportional to the volume error;
over hundreds of stiff steps these can swamp the genuine reaction-driven mass loss.
Pyrolysis.jl satisfies it *geometrically*: the node positions are states
($dz_n/dt=w_n$) and every residual evaluation recomputes $V_i$ from them and the
cross section, so the volume the scheme uses is always the swept one and no
separately evolved volume can drift (§11.5.2). With the balance rows this gives
$dV_i/dt=\theta_iV_i+\Gamma_i$ identically, with $\Gamma_i$ the sweep of the
correction, and since the masses are the states a uniform field with no sources and
no transport fluxes stays uniform.

V-MMS-02 asserts exactly that: a spatially uniform state with no reactions and
zero-flux boundaries, on a mesh compressed to 70 % of its thickness through solve's
`w_override` keyword with the uniform-stretch velocity $w^*(z,t)=(\dot s/s)\,z$, must
stay uniform. The temperature drift is pinned at $10^{-10}$ K (observed exactly zero
at 20 cells) and the relative concentration drift at $10^{-12}$ (observed zero for the
solid and $1.5\times10^{-16}$ for the gas): the relative flux of the constant field must
balance the change of the cell volumes discretely. The node trajectories
$z_n(t)=z_n(0)\,s(t)$ are known exactly, so the top-cell-centre position doubles as a
check on the prescribed-motion path ($4\times10^{-7}$ relative at `reltol = 1e-8`).
The module tests repeat the free-stream check in both energy forms at fixed steps
(`test/kernels/test_relative_flux.jl`): exact to round-off on a planar line, and
convergent at the integrator's order on cylindrical and spherical lines, whose
volumes are not linear in the radii (§11.5.2). V-MMS-03 then measures the observed
order of a smooth manufactured solution under prescribed mesh motion, and V-MMS-04
under reaction-driven motion, which closes the kinetics → strain rate → balance rows
→ geometry chain. A drift in V-MMS-02 is a moving-mesh transport defect that the
smooth-field order studies could hide inside their truncation error.

### 18.2.5 Exact burnout through a deactivation event

Deactivation (§14.3) is a discrete event: cells whose condensed mass falls below the
threshold $\varepsilon_m$ of their reference mass leave the system, their inventories
transferred to the live cell beneath. V-BO-04 checks it against an exact solution that
passes through the event. An adiabatic 10 mm slab at a uniform 600 K carries one
zero-order reaction solid → gas with no heat of reaction and $\gamma_\text{gas}=0$, on
40 cells of a moving mesh. The temperature stays at $T_0$ exactly (no source, no
gradient), so the volumetric rate $r$ is constant and every cell loses the same
fraction $r/\rho_s$ of its volume per unit time. On a mesh that rides the solid the
solid concentration stays at $\rho_s$ (the shrinkage removes volume, not
concentration) and the thickness is

$$
L(t)=L_0\exp\!\Big(-\frac{r\,t}{\rho_s}\Big).
$$

Every cell reaches the threshold at the same instant,
$t=\ln(1/\varepsilon_m)\,\rho_s/r$ (460.5 s at $\varepsilon_m=10^{-2}$), the worst case for
the event: every cell is in the near-threshold band when the first crosses, so the
whole column merges into the innermost cell in one invocation of 39 events, and that
survivor holds the column's mass (about $40\,\varepsilon_m$ of its own initial mass) at
$\rho_s$ and burns on at the exact rate. The QoIs require the thickness to follow the
exact law through the event within $10^{-6}$ relative (observed $1.5\times10^{-9}$ in
the temperature form, $1.2\times10^{-9}$ in the enthalpy form), the solid
concentration to stay at $\rho_s$ within $10^{-8}$, the temperature at $T_0$ within
$10^{-6}$ K, one invocation leaving one live cell, and the located event within
$10^{-6}$ of the exact crossing time (observed $3.1\times10^{-10}$ and
$1.8\times10^{-10}$). The case checks the event location (§14.3.5), the conservative
transfer (§14.3.3) and the restart after the event together.

---

## 18.3 Jacobian verification

### 18.3.1 Why and against what

Every integrator policy linearizes the system with the Jacobian
$J=\partial\boldsymbol f/\partial\boldsymbol u$ through the matrix
$W=J-\mathsf M/(\gamma_\mathrm{RK}\Delta t)$ (§15.5.2). The default, `Rodas5P`, is a
Rosenbrock method: each stage solves one linear system with $W$, with no Newton
iteration, and the method's order conditions are built on the Jacobian of the
right-hand side, so an inexact $J$ lowers the accuracy of the step and corrupts its
embedded error estimate. The Newton-based policies (`fbdf_policy`,
`kencarp4_policy`) iterate on $W$: a wrong $J$ leaves the fixed point of the Newton
iteration unchanged, so a converged step still solves the discrete equations, but it
degrades or destroys convergence, inflating step rejections and, at worst, stalling
the solve on an Arrhenius transient. Under either kind a Jacobian entry missing from
the sparsity pattern is silently dropped.

The Jacobian is colored sparse forward-mode automatic differentiation of the whole
residual on a pattern built from the mesh connectivity (§15.6–§15.11). The reference
against which "true" is defined is `dense_jacobian(def, ws, u, t)`
(`src/ColoredJacobian/oracle.jl`): every column of the same `residual!` the
integrator evaluates, differentiated by forward mode without the pattern or the
coloring. It is correct by construction and slow. Verification compares the two,
and the reference itself is compared with central differences, an oracle
independent of automatic differentiation.

### 18.3.2 The comparison functions

Three functions of `Pyrolysis.Internal` (from `ColoredJacobian`) carry the check:

- `dense_jacobian(def, ws, u, t) -> Matrix{Float64}`, the reference;
- `missing_entries(J_ref, pattern) -> Vector{Tuple{Int,Int}}`, the nonzero entries
  of the reference that the pattern does not contain, empty when the pattern
  contains it;
- `compare_entries(J, J_ref; rtol = 1e-12, atol = 0.0)
  -> (passed, max_abs_dev, worst_ratio, worst_index)`, the per-entry comparison.

### 18.3.3 The per-entry criterion and the containment tests

`compare_entries` requires, for every entry,

$$
|\Delta J_{ij}|\;\le\;\texttt{atol}+\texttt{rtol}\cdot
\max\big(|J^{\rm ref}_{ij}|,\ \max\nolimits_{j'}|J^{\rm ref}_{ij'}|\big).
$$

The row-magnitude floor keeps honest round-off at a structural zero of a stiff row
from failing (a $10^{-9}$ residue in a $10^6$-scale Arrhenius row passes), while a
wholly missing block in a small-magnitude row (a velocity row lies many orders of
magnitude below a mass row) is measured against its own row's scale and fails
regardless of the problem's size. A global normalization by $\max|J^{\rm ref}|$
would make pass and fail scale-dependent: a small-magnitude row could lose an entire
block behind the Arrhenius rows. `test/jacobian/test_pattern.jl` checks the
criterion's behaviour on a stiff row above a small one.

There is one Jacobian and one regime: the colored Jacobian must contain every nonzero
of the reference and equal it per entry at `rtol = 1e-12`. The containment tests
(`test/jacobian/test_containment.jl`; §15.13) require, at every state they visit, a
finite residual and reference, an empty `missing_entries`, a colored Jacobian with
the pattern's structure, and that agreement. They cover fixed and moving columns of a
PMMA-like and a permeable wood-like material, without radiation, with surface
absorption and with in-depth absorption, under Fickian and Darcy–Fick transport, in
both energy forms and at 6 and 8 cells; at the set-up state, at uniform hot, hot and
gas-heavy, partially decomposed and mixed reacting states, and at trial states
outside the physical domain (perturbed, zero and negated states); the chain-rule
regimes of the smooth temperature gates within a kelvin of both edges, of a gas heat
capacity that depends on temperature, of the gas viscosity of a hot Darcy column, of
permeability mixing, of a single cell and of 200 cells; the cross section's
back-reaction under a lateral-shrinkage law; cylindrical and spherical meshes; a
prescribed mesh velocity; the states after `deactivate!` with in-depth radiation; the
boundary mass flux on the energy row; and a surface absorptivity and an emissivity
that depend on temperature. The dense reference equals central differences to
$10^{-6}$ per entry with the same row floor (the testset "the dense Jacobian equals
central differences"). Each kernel's tests repeat the containment check for its own
physics (`test/kernels/test_gas_transport.jl`, `test_radiation.jl`,
`test_relative_flux.jl`), and `test/assembly/test_residual_jacobian.jl` checks the
assembled residual and Jacobian on a stretched mesh under closed-form, convective and
combined boundary conditions, on and off the physical states. A residual and a
colored Jacobian allocate nothing after warm-up (`test/jacobian/test_allocations.jl`),
and JET's `@test_call` analyses both entry points at the concrete types `setup` builds
(`test/type_analysis/`).

### 18.3.4 Checking a user's problem

No solve option compares Jacobians; the same functions check any problem from a
script, at any state, through `setup` and `integrate`:

```julia
using Pyrolysis
using Pyrolysis.Internal: jac_prototype, jacobian!, dense_jacobian, missing_entries,
    compare_entries

# A charring slab on a moving mesh, integrated to a reacting state.
virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1500.0, k = 0.2, λ = 1e-6)
char = SolidComponent(:char, ρ = 200.0, c = 1200.0, k = 0.1, λ = 1e-6)
gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5)
rxn = Reaction(:decomp, 1 => (2, 3); A = 1e13, E = 1.88e5, h = 8.7e5, yields = (0.2, 0.8))
material = Material(name = :slab, components = (virgin, char, gas), reactions = (rxn,))
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(5e-3, 20),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(Inf), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, 120.0),
    output_times = [0.0, 60.0, 120.0],
)

cp = setup(problem; use_ale = true)
sol = integrate(cp)
u, t = sol.raw_solution.u[end], sol.t[end]           # the state the solve ended on

J = jac_prototype(cp.jacobian)                       # the colored Jacobian's sparse matrix
jacobian!(J, cp.jacobian, u, t)
J_ref = dense_jacobian(cp.def, cp.ws, u, t)          # dense forward-mode reference
println(sol.retcode, ", ", size(J_ref), ", ", count(!iszero, J_ref), " nonzero reference entries")
println(isempty(missing_entries(J_ref, cp.jacobian.pattern)))
println(compare_entries(Matrix(J), J_ref; rtol = 1e-12))
```

prints

```
Success, (122, 122), 880 nonzero reference entries
true
(passed = true, max_abs_dev = 0.0, worst_ratio = 0.0, worst_index = (1, 1))
```

The colored Jacobian evaluates the residual on dual numbers in batches of its color
count and the reference in batches of the same width, so the two agree to the last
bit here. The dense reference costs one residual evaluation per batch of columns, so
it is for checks, not for production.

### 18.3.5 Deactivated states

Deactivation (§14.3) keeps the state's size and turns the rows of the dead cells and
nodes into identity rows, so the colored Jacobian is the one the integrator uses
through burnout and the containment check covers it.
`test/deactivation/test_deactivation.jl` ("masked rows are identity rows and the live
block is the compacted problem's") builds a six-cell problem and a problem of only its
live cells, on fixed and moving meshes, under in-depth radiation, Darcy–Fick transport
with a surface balance, and a lateral-shrinkage law, in both energy forms, with one and
two dead cells, at the set-up state and a perturbed one. It requires the live rows'
residual to equal the compacted problem's to $10^{-11}$ of the magnitude of their
terms; the live block of the dense Jacobian to equal the compacted problem's per entry
to $10^{-11}$, with no live row reading a dead unknown; every dead differential row to
be zero with a zero derivative and every dead algebraic row to equal its unknown, an
identity row; and the colored Jacobian to equal the dense reference to $10^{-12}$.
`test/jacobian/test_containment.jl` repeats the containment check after `deactivate!`
under in-depth radiation.

---

## 18.4 Conservation-based verification

The ledgers of §16 serve as verification: the module tests and the burnout cases
drive them per accepted step, and any solve can carry them
(`callback = ConservationRecorder()`, then `conservation_audit(sol, rec)`; without a
recorder the audit integrates the saved series by the trapezoid, §16.8). The checks
divide into those that must close by construction (a failure is a bug) and those that
quantify an approximation (a large value is a regime warning, not a bug).

### 18.4.1 Mass ledger

The component masses are the states and every face flux is scattered antisymmetrically
(§13.3.3), so the live mass changes only by the outflow through the live boundary faces
and by the masses deactivation removes (§16.3). The recorder integrates the outflow over
each accepted step by three-point Gauss–Legendre quadrature on the step's interpolant,
up to an event that ends a step early (§16.5.3), and the audit writes the defect

$$
d_m=\big[m_\text{live}(t_1)-m_\text{live}(t_0)\big]+\int_{t_0}^{t_1}\dot m_\text{out}\,dt+\sum_em^\text{lost}_e ,
$$

which closes to the integrator's tolerance; a defect that does not fall with the
tolerance is a flux-assembly or event bookkeeping bug. The module tests
(`test/outputs/test_conservation.jl`) assert at `abstol = 1e-8`, `reltol = 1e-6` a mass
defect of at most $10^{-8}$ of the inventory on a gasifying column; across deactivation
events on fixed and moving meshes a fall by more than tenfold from `reltol = 1e-6` to
`1e-8`; and on a venting column whose cross section shrinks a relative defect of at most
$10^{-7}$ with the outflow integrated at the current face area (§16.5.4). The burnout
cases carry their own per-step ledger of the same construction (the live mass drop
against the integrated outflow and the vented mass): V-BO-01 and its variants require it
within $3\times10^{-5}$ (it closes to at most $5.2\times10^{-7}$ in the temperature form
and $3.2\times10^{-6}$ in the enthalpy form on the moving-mesh cases at their CI
resolution), and V-BO-02 within $3\times10^{-4}$, since each deactivation on a fixed mesh
leaves a gas transient (observed $6.1\times10^{-5}$ at 60 cells, $1.9\times10^{-8}$ at
`reltol = 1e-8`). V-GT-01 asserts total-mass conservation of a diffusing gas pulse in a
sealed column, and V-CPL-03 total-mass conservation in every shrink/swell permutation,
the sealed char-and-gas column among them.

### 18.4.2 Energy ledger

The energy ledger balances the live sensible-enthalpy inventory
$H_\text{live}=\sum_i\sum_jm_{j,i}h_j(T_i)$ (the extended caloric relation, §16.4)
against the enthalpy rate read from the energy rows (the boundary energy fluxes of the
live boundary faces, the reaction heat, the absorbed radiation and the formation terms)
and against what the deactivation events remove (§16.5, §16.7). It closes **exactly in
both energy forms**: in the enthalpy form $H_\text{live}-\int\dot H_\text{in}\,dt$ is a
linear invariant of the semi-discrete system, which a fixed-step solve closes to
round-off (at most $10^{-13}$ of the integrated rate in
`test/outputs/test_energy_ledger.jl`); in the temperature form the closure converges at
the integrator's order (observed orders 2.1, 2.5 and 3.2 under Fickian transport, 2.7,
3.3 and 4.3 under Darcy–Fick, with fixed steps of a fourth-order method halved three
times; §16.5.2). With the recorder, the module tests assert a relative energy defect of
at most $10^{-6}$ in both energy forms on a gasifying column whose escaping gas carries
its enthalpy out, and at the tighter tolerance across deactivation events (§16.5.4). A
defect that does not fall with the tolerance is a bookkeeping bug.

**The boundary flux.** The ledger's boundary term is the energy flux of the live
boundary faces: the condition's flux conducted into the cell (§13.7.5) plus the
enthalpy the outward mass fluxes carry at the surface temperature. The ledger *audits
bookkeeping, not double counting*: it reads the same kernel contributions the residual
assembles, so if the residual counted a mechanism twice the ledger would still close.
Detecting a double count requires an *independent* check (§18.4.3).

### 18.4.3 Independent checks of the gas-enthalpy coupling

The gas enthalpy is carried by the face energy flux: on an interior face the
gas-transport kernel adds $\sum_jJ_{j,f}h_j(T_f)$ at the distance-weighted face
temperature, on a boundary face the thermal boundary kernel adds the enthalpy of the
outward fluxes at the surface temperature, and the surface balance carries no
gas-enthalpy term (§7.4.3, §13.7.4–§13.7.5). A surface term
$\sum_g\dot m_g\,\Delta h_g(T_s,T_\infty)$ added on top would count the escaping
enthalpy twice, and the ledger of §18.4.2 would record both and still close. The
coupling is therefore pinned by checks built independently of the residual:

- **The energy row term by term** (`test/outputs/test_energy_ledger.jl`, "energy row
  term by term against the ledger"): in both energy forms, under Fickian and Darcy–Fick
  transport, on uniform and stretched meshes, every face energy flux equals the
  conduction (Fourier between cells with the distance-weighted harmonic-mean
  conductivity; the condition's conducted flux on the boundary) plus the enthalpy its
  mass fluxes carry at the face temperature, to $10^{-12}$ relative, and every cell's
  energy row equals the sum of its conduction, carried, source and formation terms to
  $10^{-11}$ of their scale. The carried term is not small and the surface temperature
  departs from the cell temperature by more than 1 K, so an omission, a double count or
  enthalpy taken at the cell temperature would show.
- **No gas enthalpy in the surface balance**
  (`test/boundary_conditions/test_boundary_temperature.jl`, "Outward gas flux enters
  only through blowing"): under a prescribed flux, the solved surface temperature and
  the conducted flux are identical with and without a gas flux across the face.
- **Isothermal transport** (`test/kernels/test_gas_transport.jl`, "isothermal transport
  leaves the temperature at rest"): gas diffusing through an isothermal column leaves
  every temperature at rest in both energy forms.
- **An exact solution that would see a spurious term.** V-CPL-02 sets the gas heat
  capacity equal to the solid's, which makes the Landau travelling-wave identity
  $v=v_{ss}(T_s)$ exact for any kinetics; a spurious advection of the temperature at
  the mesh velocity on a mesh that rides the collapsing solid would break it by about
  $c\,(T_s-\bar T_\text{form})\approx R T_s^2c/E$ per kilogram, a few percent, and the
  case pins the identity at $2.1\times10^{-3}$ (observed $7.1\times10^{-4}$ at 160 cells).

---

## 18.5 Convergence studies

A scheme can be conservative and have a correct Jacobian yet still be *under-resolved*.
Convergence studies confirm that the discretization error shrinks at the designed rate
as the mesh, the time step or the integrator tolerance is refined, and they let the user
pick a resolution that is converged for the quantity of interest (QoI) at acceptable
cost.

### 18.5.1 Spatial convergence (uniform and stretched)

The cell-centered finite-volume scheme of §13 is second-order accurate in space for
smooth solutions: the two-point face gradient and the harmonic-mean face coefficient
are $O(\Delta z^2)$ on a uniform mesh, and the distance-weighted forms keep second order
on smoothly stretched meshes (§13.6). The verification suite asserts the orders against
exact solutions where it can: expected order 2 (each observed order above 1.5 and their
mean above 1.8) for V-HC-01, 02c, 03 and 05, V-RAD-01, V-GT-01, V-CPL-04 and V-MMS-01
and 03, and above 1.7 with a mean above 1.9 for V-MMS-04; V-HC-01 observes 1.991 and
1.998. V-GT-02 expects first order (observed 0.89 and 0.94), the signature of the
first-order upwind Darcy advection (§13.7.3, §18.7.3). The boundary faces are first-order
(§13.12); the manufactured-solution cases impose fluxes there so that they do not cap
the observed order.

With an exact solution unavailable, the finest mesh is the reference
(Richardson-style self-convergence): `examples/convergence_douglas_fir.jl` solves the
Douglas-fir cone case of §18.6.2 on 480, 240, 120, 60, 30 and 15 uniform cells (finest
first, the reference) at `abstol = 1e-9`, `reltol = 1e-7`, and reports per-QoI relative
errors

$$
\Delta f(n)=\frac{|f(n)-f_\mathrm{ref}|}{|f_\mathrm{ref}|}\times100\%
$$

for the scalar QoIs (peak mass-loss rate and its time, the time to ignition, the total
mass loss, the peak surface and back-face temperatures) and curve-level $L_2$ relative
errors

$$
\frac{\lVert f(n)-f_\mathrm{ref}\rVert_2}{\lVert f_\mathrm{ref}\rVert_2}\times100\%
$$

for the mass-loss rate, the surface temperature and the back-face temperature
histories. Plotting $\Delta f$ against the cell count on log–log axes recovers the slope
$-p$; a slope near $-2$ confirms second-order convergence, a shallower slope flags a
first-order contribution (the Darcy advection, a surface-absorbed flux resolved by the
boundary half-cell).

Stretched and graded meshes (`create_stretched_mesh`, `create_mesh_from_nodes`, §13.9)
put the resolution where the steep thermal and reaction fronts live. V-CPL-01 runs the
Stefan problem on a mesh stretched toward the cold wall (about 14 μm at the wall), the
energy-row test runs on stretched meshes (§18.4.3), the relative fluxes converge above
order 1.5 on cells alternating between $h$ and $2h$ (`test/kernels/test_relative_flux.jl`),
and V-BO-06 burns the V-BO-01 layer on a mesh graded threefold toward the surface under
V-BO-01's references. An incorrect non-uniform treatment shows up as a loss of order on
the stretched run that is absent on the uniform one.

### 18.5.2 Temporal convergence

The default integrator `Rodas5P` is a fifth-order Rosenbrock method with an embedded
error estimate; `fbdf_policy` (variable-order BDF) and `kencarp4_policy` (fourth-order
ESDIRK) are the alternatives (§15.5.1). A temporal study fixes the mesh and sweeps the
step: V-MMS-01's report mode runs a fixed-step sweep (`adaptive = false`, steps 0.5,
0.25, 0.125 and 0.0625 s at 80 cells) against the manufactured solution, and the
energy-ledger test measures the temperature form's closure at fixed steps halved three
times (§18.4.2). In practice the Arrhenius transients, the sharp onset of decomposition
as the temperature rises, dominate the temporal error and force the adaptive controller
to small steps; the smooth temperature gates and the one-sided depletion limiter
(§6.2–§6.3) keep the residual and its Jacobian differentiable through them, so the error
estimator behaves. A deactivation event is located inside its step at the first root of
the threshold margin and the integration restarts from the state the event leaves
(§14.3.5); V-BO-04 locates it within $3.1\times10^{-10}$ of the exact crossing
(§18.2.5), and the module tests require a located event's time to converge with the
tolerance on fixed and moving meshes.

### 18.5.3 Integrator-tolerance study

The step size is chosen adaptively to meet `abstol`/`reltol` (defaults $10^{-8}$ K on the
temperature rows, with the mass and enthalpy rows' tolerances derived from the cell's
reference condensed mass and `threshold`, and `reltol = 1e-6`; §15.5.6), so the tolerance
pair is the practical temporal-error control. A tolerance study sweeps `reltol` (and
`abstol` in proportion) and tracks a QoI to its tolerance-converged plateau; the
convergence example uses tightened tolerances (`abstol = 1e-9`, `reltol = 1e-7`) so that
the *spatial* error dominates the comparison and is not masked by temporal error. The
ledger audit of §18.4 provides an independent confirmation that the tolerance is tight
enough: its defect falls with the tolerance, so a defect that does not shrink as `reltol`
is tightened indicates a bookkeeping issue rather than under-resolution. The three knobs
($\Delta z$, the step, the tolerance) are studied separately so that the dominant error
source is unambiguous; only when all three are converged is a quoted QoI a property of
the *model* rather than of the *discretization*.

**Burnout: joint refinement of the mesh and the tolerance.** Through a burnout the
deactivation threshold is a numerical parameter too, and the end-of-layer burst (the
largest per-step mass-loss rate just after a layer's deactivation) is the sensitive
quantity. The burnout path (the equidistribution correction with its layer time set at
activation, the near-threshold band, located events) is verified under joint refinement
of the mesh and the tolerance on five decks built on V-BO-01: the uniform mesh, the mesh
graded threefold toward the surface (V-BO-06), the mixed layer (V-BO-07), the 1 mm layer
(V-BO-08), and the 5.9 mm layer with a convective loss in place of re-radiation and an
adiabatic back. Inside the supported envelope every QoI (complete burnout, the mass lost,
the peak and its time, the burnout time, the mass ledger) is within tolerance at every
resolution and the burst converges (the largest per-step excess of $N$ cells over $2N$ and
the successive $L_1$ differences falling at least 1.5× per doubling): $\varepsilon_m=3\times10^{-3}$
with $\varepsilon_\text{band}=0.1$ on all five decks from 30 cells across the layer (10 on
the 1 mm layer); the default $\varepsilon_m=10^{-2}$ with $\varepsilon_\text{band}=0.1$ on the
single-solid layers from 30 cells (60 on the convective-loss deck) but not on the mixed
layer (§18.8.2); at `reltol` $10^{-6}$ and $10^{-7}$, in either energy form, and with
either activation bound pair $(0.3,0.1)$ and $(0.7,0.3)$. At 30 cells the largest per-step
rate just after a layer's deactivation can exceed the main peak, up to threefold on the
graded mesh; it falls below the main peak by 60 to 120 cells, and the peak sampled at the
output times is within tolerance from 30 cells. V-BO-01 at 25 and 100 kW m⁻², at half and
twice its pre-exponential factor and with its flux stepped after activation, at 30 to 120
cells, passes the same QoIs, the hot decks with two invocations. Materials outside these
decks (other materials, thicker or charring layers, higher heating) are not inside the
envelope by this study. The full statement with its numerical parameters is §14.2.4 and
User Guide chapter 10.

---

## 18.6 Validation against experiments and reference codes

Validation compares converged Pyrolysis.jl predictions against measured data and
against the output of established condensed-phase pyrolysis codes. The notation of
the reference codes is mapped to ours on first use; the comparisons draw on the
governing-form, kinetics, transport and radiation contrasts documented throughout
§3–§9.

### 18.6.1 Cone calorimeter: HRR/MLR/$T_s$

The cone calorimeter is the canonical bench-scale validation experiment: a sample is
exposed to a constant radiant flux (typically 25–75 kW m⁻²) and the mass-loss rate
(MLR), the surface temperature and the oxygen-consumption heat-release rate (HRR) are
recorded. Pyrolysis.jl, being a condensed-phase solver, predicts the MLR and $T_s$
directly and an HRR *proxy* through an effective heat of combustion times the
gas-generation rate; the gas-phase flame is not modelled. `examples/06_cone_calorimeter.jl`
(PMMA under a 25 kW m⁻² cone, an incident flux with absorptivity 0.96 plus
re-radiation and convection on the top face, a convective mass condition, back-face
re-radiation and convection with an impermeable back, single-step Arrhenius
kinetics, temperature-dependent properties, a moving mesh with deactivation) is the
reference workflow; its variant `06_cone_calorimeter_surface_absorption.jl` names
`radiation_model = SURFACE_ABSORPTION` explicitly, which gives the same residual as
the default `NO_RADIATION` for an incident flux absorbed at the face (§8.2).

The validation QoIs are the peak MLR and its timing, the residual (char) mass
fraction, a time-to-ignition surrogate (a surface-temperature or gas-generation
threshold), and the $T_s(t)$ history. Two implementation points materially affect
quantitative agreement and are themselves verified:

- **Surface temperature in a cell-centered scheme.** Temperatures live at cell
  centres, but a $T^4$ radiation condition must be evaluated at the *surface*. The
  thermal boundary kernel solves the surface energy balance
  $q_\mathrm{bc}(T_s)=k\,(T_s-T_\mathrm{cell})/d$ for $T_s$ by Newton's method at
  every residual evaluation (§12.2, §13.7.5) and stores it; the solution's
  `surface_T` reads it at every saved time. Validation against a measured $T_s$
  requires this solved value, not the boundary cell's temperature; V-HC-01, V-HC-02
  and V-RX-04b verify it against exact solutions.
- **MLR normalization.** `sol.extras.surface_mass_flux` is the outward mass flux of
  the `:top` faces per unit current face area, and `sol.extras.MLR_gauge` the mass
  leaving the live boundary faces per unit *initial* area of the `:top` faces, the
  cone gauge convention. Under a lateral-shrinkage law the two differ by $A(t)/A_0$
  (§10), and a prediction compared with a measurement reported per initial sample
  area uses the gauge series.

### 18.6.2 Douglas fir and the convergence case

`examples/07_douglas_fir_cone.jl` is the most fully featured validation case:
multi-step wood decomposition kinetics (five sequential steps through four
intermediates to char, and moisture evaporation), surface absorption of the incident flux (wood is opaque),
Darcy–Fick gas transport with the pore gas initialized at ambient pressure
(`background_gas_concentration`; the species molar masses enter the pressure
closure, and all gases share one diffusive transfer coefficient, §5.11), and a
moving mesh tracking the shrinking sample, solved at `abstol = 1e-9`,
`reltol = 1e-7` with the default integrator and the colored Jacobian. The companion
`convergence_douglas_fir.jl` (§18.5.1) establishes the mesh resolution at which the
validation QoIs are converged, so that any residual disagreement with experiment is
attributable to physics or parameters rather than to under-resolution. Validation
here checks the MLR$(t)$ shape (the characteristic two-peak signature of
hemicellulose/cellulose and lignin decomposition), the char yield, and the back-face
temperature rise (a sensitive probe of the in-depth conduction and of the effective
conductivity's mixing rule, §5.5).

The FDS-style kinetic recovery of §18.2.2 is the practical link between TGA data and
the cone model: $(A_1,E_1)$ are recovered from a measured TGA peak temperature and
peak rate through the inflection relation, then carried into the cone run. This is a
*calibration* step (the kinetics are fitted to TGA), so the cone comparison that
follows is a genuine *validation* only insofar as the cone data were not used in the
fit, the distinction §18.1.1 insists on.

### 18.6.3 Pressure-treated wood and moisture transport

`examples/08_pressure_treated_wood_cone.jl` (and its ThermaKin comparison
`ptw_comparison.jl`) validates the Darcy-driven transport path: a permeable charring
solid develops an internal pressure $P=\sum_g\xi_jR_gT/(M_j\varphi)$ (§9.4), and the
Darcy velocity $v_g=-(\kappa/\mu)\,\partial P/\partial z$ (§9.5) advects the pyrolysis
gases outward. The validation QoIs are the MLR augmentation and the internal-pressure
history relative to a Fickian-only run; `08_darcy_fick_demo.jl` isolates the
pressure-driven transport for a controlled check. Moisture transport
(`moisture_hos_comparison.jl`) validates the state-dependent heat-of-sorption path:
the evaporation reaction carries a `StateDependentProperty` heat $h(T,\xi)$ that
depends on the local moisture concentration (§6.6), and the moisture front, the
evaporation plateau in $T(t)$ (the constant-temperature stall while water boils
off) and the dried-mass profile are the comparison targets.

### 18.6.4 Comparison to ThermaKin2Ds, Gpyro, and FDS

Code-to-code comparison (`examples/literature_comparison/`,
`examples/macfp/macfp_gasification.jl`, and the cross-code overlays of the
verification suite, §18.1.3) is run alongside experimental validation because it
isolates *modelling-convention* differences from *numerical* ones. The three
reference codes differ from Pyrolysis.jl in documented, specific ways; mapping their
notation to ours:

- **ThermaKin2Ds.** Maps its stoichiometric coefficients $\theta_i^{\,j}$ to our
  $\nu_{i,j}$ and its component concentrations to our $\xi_j$, and publishes its
  heats of reaction with $h>0$ exothermic (the opposite of ours). ThermaKin's energy
  equation (its Eq. 17) stores sensible heat in every component, gas included, and
  carries the gas enthalpy by the convective term $c_g\,(J_g\cdot\partial T/\partial x)$;
  Pyrolysis.jl stores every component in $C_i$ and carries the gas enthalpy as the face
  flux $\sum_jJ_{j,f}h_j(T_f)$, the conservative form of the same term (§13.3.4), so the
  convective coupling is a clean cross-check. ThermaKin supports *bimolecular*
  reactions; Pyrolysis.jl enforces single-reactant kinetics (§6.5), so multi-reactant
  ThermaKin schemes must be re-expressed as reaction sequences before comparison.
- **Gpyro.** Maps its pre-exponential $Z$ to our $A_i$, its permeability $K$ to our
  $\kappa$, its volume fractions $X_\alpha$ to our $v_j$, and its volumetric rate
  $\dot\omega'''$ to our $r$. Gpyro is fully conservative (control-volume balances of
  $\bar\rho\Delta z$, $\bar\rho Y_i\Delta z$, gas mass and energy) and uses
  volume-fraction property averaging with power-law temperature dependence;
  Pyrolysis.jl's WEIGHTED/PARALLEL/SERIES/BRUGGEMAN mixing rules (§5.5) span and
  generalize Gpyro's averaging. Gpyro's Darcian flux
  $\dot m''=-(K/\nu)(\partial P/\partial z-g\rho_g)$ matches our Darcy form (§9.5) up
  to the buoyancy term, which Pyrolysis.jl omits (negligible for thin condensed-phase
  slabs). Gpyro selects swelling and shrinkage mechanisms by solid-fraction state;
  Pyrolysis.jl's lateral-shrinkage law $A(t)=A_0\,\mathrm{law}(\bar\chi)$ (§10.5) is the
  1D analogue. Gpyro's solver is fully implicit backward Euler with Patankar
  linearization and TDMA; Pyrolysis.jl integrates with `Rodas5P` and the colored
  Jacobian, so code-to-code agreement on a steady or slowly varying case validates the
  *physics* independent of the time integrator. Gpyro 0.8200 runs of V-HC-01, V-HC-03,
  V-RAD-01, V-RX-01 and V-CPL-01, and a Gpyro V2 run of the V-BO-01 deck, are overlaid
  in report mode.
- **FDS solid phase.** Maps its absorption coefficient $\kappa_s$ to our $\alpha$ and its
  pre-exponential $Z$ to our $A_i$. FDS solves 1D Fourier conduction with Arrhenius
  kinetics (with $\rho^n$, $O_2^n$ and $T^n$ rate forms) and a two-flux in-depth
  radiation model; Pyrolysis.jl's in-depth model is the one-directional Beer–Lambert
  beam of the incident flux (§8.3), so the two agree where in-depth emission is
  negligible. FDS carries no in-depth gas, so its energy equation is the
  condensed-phase limit of ours; its TGA-based kinetic recovery (peak temperature and
  rate → $A$, $E$) is the recipe of §18.2.2 and §18.6.2. FDS's surface energy balance
  $-k\,\partial T/\partial x=q_c''+q_r''$ is the multi-mechanism analogue of our surface
  Newton balance (§12.2). FDS 6.11.0 runs of V-HC-01, V-HC-02, V-RX-01, V-RX-02 and
  V-CPL-03 are overlaid in report mode.

In all three comparisons, agreement *after* matching conventions (gas treatment,
radiation model, kinetic form) is the verification that the shared physics is
implemented consistently; *residual* disagreement after convention matching is
attributable to the documented modelling differences, not to error.

---

## 18.7 Regime-validity checks

A model is only valid in the regime where its approximations hold. The checks below
tell the user when a run is leaving the validated regime; they are validation aids,
not verification checks, because they assess the *applicability of the equations*,
not the *correctness of the solve*.

### 18.7.1 Optical thickness and the radiation model

The optical thickness of the column,

$$
\tau_\mathrm{total}=\int\alpha\,dz=\sum_i\alpha_i\,\Delta z_i,\qquad
\alpha_i=\sum_j\alpha_j\,\xi^+_{j,i},
$$

with $\alpha_j$ the components' mass absorption coefficients (m² kg⁻¹), selects the
radiation model (§8): for $\tau_\mathrm{total}\gg1$ the incident flux is absorbed within a
few cells of the surface and `SURFACE_ABSORPTION` is adequate; for an intermediate or
thin column `BEER_LAMBERT` distributes it in depth; for $\tau_\mathrm{total}\ll1$ the
medium is effectively transparent. The package computes no regime verdict; $\tau$
follows by hand from the components' absorption coefficients and the cell
thicknesses, and a solve under `BEER_LAMBERT` reports the absorbed and the transmitted
intensity at every saved time (`sol.extras.volumetric_absorption`,
`sol.extras.radiation_transmitted`). A transmitted share that is not small means that
flux leaves through the back face without heating the column, a regime that a back-face
condition or a substrate layer must represent; V-BO-03 exercises it on a transparent
substrate (its transmitted share exceeds 1 % of the entering flux as the PMMA thins).
The absorbed-plus-transmitted identity of §18.2.3 is the underlying conservation check;
the regime judgement layers on top of it.

### 18.7.2 Pore-gas thermal equilibrium

The energy equation stores the pore gas at the matrix temperature (§7.3). That
one-temperature assumption holds when the gas equilibrates with the pore walls within
a cell, which is the case for pore-scale heat-transfer lengths far below the cell size —
every validated regime here (1 atm, pore sizes well below the mesh spacing, gas
residence times of milliseconds). A run with coarse open pores or very fast, hot
throughput, where the gas does not equilibrate with the walls, is outside the model:
there is no two-temperature gas option (§3.9), and the prediction should be treated
with caution.

### 18.7.3 Péclet number and upwind diffusion

The Darcy–Fick species transport uses first-order upwinding for the advective part
(§9.3, §13.7.3). First-order upwind introduces a numerical diffusion of order
$|v_g|\Delta z/2$, so its importance is measured by the cell Péclet number
$\mathrm{Pe}=|v_g|\Delta z/\lambda$; when $\mathrm{Pe}\gg1$ the scheme stays robust but the
artificial diffusion can smear a gas-composition front. The convergence study (§18.5.1)
exposes it as a first-order contribution to the spatial error of pressure-driven cases;
V-GT-02 measures it directly (observed orders 0.89 and 0.94 against the exact $P^2$
profile). The practical regime check is to confirm that the cell Péclet number is
moderate over the resolved front; if it is large, refine the mesh near the front (a
stretched mesh) until the upwind diffusion is below the physical diffusion. The
conduction analogue is benign (the conductive flux is a two-point central difference,
second-order), so this regime concern is specific to the advective gas transport.

### 18.7.4 Threshold sensitivity

Two thresholds govern near-zero-concentration behaviour, and a validated run should be
insensitive to both within their physically negligible range.

- The **depletion limiter** threshold $\xi_\mathrm{th}$ (`DepletionLimiter(threshold)`,
  default 1 kg m⁻³, §6.3) sets where the rate roll-off begins; the factor is one to
  within $2e^{-2\xi/\xi_\mathrm{th}}$ above it.
- The **deactivation** threshold $\varepsilon_m$ (`threshold`, a fraction of a cell's
  reference condensed mass, default $10^{-2}$, §14.3.2) sets when a cell leaves the
  system, together with the near-threshold band (`band`, default $\max(\varepsilon_m,0.1)$
  where the equidistribution correction moves a layer, otherwise $\varepsilon_m$) that
  takes the cells close to the threshold in the same event. The events are located
  inside the step (`located_events = true`).

A sensitivity check varies each threshold by a factor of a few to ten and confirms that
the QoIs (peak MLR, char yield, burnout time) move by less than their tolerance. A
strong dependence means the limiter is acting in a physically active region rather than
only damping the stiffness near true depletion, or that the deactivation numerics sit
outside the supported envelope (§18.5.3); lower the threshold in either case (the
envelope holds at $\varepsilon_m=3\times10^{-3}$ on every deck it covers). Both
thresholds act smoothly or through located events: the limiter is $C^1$ (§17.6), and a
located event's time is a smooth function of the parameters, so finite-difference
sensitivities remain valid across a deactivation, with the caveats of §18.8.3.

---

## 18.8 Known discrepancies and re-validation needs

The following are documented incompletenesses that affect, or require re-checking of,
validated results.

### 18.8.1 Surface gas-enthalpy convention and calibrated parameters

The gas leaving through a face carries its sensible enthalpy at the face temperature in
the face energy flux, and the surface balance has no gas-enthalpy term (§7.4.3,
§18.4.3), since a surface term would count that enthalpy a second time. A balance that
adds a surface term on top of a carried enthalpy artificially depresses the surface
temperature, so parameter sets calibrated in a model that does so read a *higher* $T_s$
and mass-loss rate here (§1.8, item 7). Imported calibrations must therefore be
re-validated against data under this convention. The fitted kinetic and property
parameters of the bundled calibrated examples (the inverse-analysis cases and the cone
runs) are in that position, together with the flame heat-flux feedback (absent,
§1.8.1) and the blowing reduction of the convective coefficient (`ConvectiveBC`'s
blowing model, §12.3.3) that a recalibration should include. This is a validation-state
issue, not a verification failure.

### 18.8.2 The end-of-layer burst outside the supported envelope

The burnout path is verified inside the envelope of §18.5.3 (§14.2.4). Outside it the
QoIs can still pass while the end-of-layer burst (the largest per-step mass-loss rate
just after a layer's deactivation) does not converge with resolution. V-BO-07, the layer
of mixed composition, is the documented case: at the default $\varepsilon_m=10^{-2}$ its
QoIs pass at the CI resolution, but under refinement the largest per-step excess of $N$
cells over $2N$ stops falling between 120 and 240 cells, while at
$\varepsilon_m=3\times10^{-3}$ it falls by a factor of two per doubling. When the first
cell of that layer crosses $\varepsilon_m$, the others still hold up to about 6 % of their
initial mass (the slow solid concentrated near the surface), so one invocation merges
cells far from their threshold. The convective-loss deck at $\varepsilon_m=10^{-2}$ takes
two invocations from 60 cells. A run whose layer is not of the decks' kind, or whose
invocation count is above one, should be checked by refinement at a lower threshold
before its burst is trusted; the integrated QoIs (mass lost, burnout time) are the robust
ones.

### 18.8.3 Sensitivities across deactivation events

Deactivation is a sequence of discrete state changes (§14.3). The Jacobian runs through
them (§18.3.5), and a located event's time is a smooth function of the parameters (it is
the first root of the threshold margin), so the trajectory after the event moves
smoothly with the parameters as long as the event takes the same cells. Two things break
that smoothness (§17.7.3): an objective that samples the solution at fixed times sees
the event's burst move across its sample times, which adds a roughness of the period of
one sample interval in the event time (on V-BO-01 the deactivation moves 0.46 s per
$10^{-3}$ of flux against the reference's 1 s samples, and the $10^{-3}$ difference
departs by 6 %; §17.2.2); and a perturbation that changes which cells an invocation
takes, or splits one invocation into two, changes the merged state by a finite amount.
The finite-difference convergence check of §17.2.2 selects the usable steps;
`sol.statistics.invocations` and `sol.events` compared between the perturbed solves
detect a change of decision.

### 18.8.4 Heuristic and unvalidated closures

Several constitutive closures are heuristic and carry their own validation debt: the
Lennard-Jones parameter estimates for the Chapman–Enskog diffusivity
($\sigma\approx3.5+0.5\log_{10}(M/0.016)$, $\varepsilon/k\approx150+100(M/0.044)$, §5.10)
are explicitly flagged as unvalidated; the lateral-shrinkage laws (§10.5) are empirical
rather than first-principles; and the single gray radiation model is a simplification.
Validation cases that depend sensitively on these closures should report that
dependence.

### 18.8.5 Stale-flux regression

One silent failure mode carries a standing regression guard. The incident radiative flux
is evaluated at *every* residual evaluation with no constancy detection (§8.3.4, §13.5.2),
so a time-varying flux can never be frozen at its initial value. A regression would be
silent, so the modulated-flux comparison under $q_0(1+\sin2\pi t)$ against a constant
$q_0$ (§18.2.3) is kept as a permanent check in `test/kernels/test_boundary_conditions.jl`.

---

## 18.9 Reproducibility

Verifiable results must be *reproducible*, which requires a controlled software
environment and an automated quality gate.

### 18.9.1 Quality tooling

Static and dynamic quality tools run in isolated environments, so the package's own
manifest is not pinned to their dependencies (`.github/workflows/quality.yml`):

- **JuliaFormatter** — `julia --project=scripts scripts/format.jl` enforces one
  canonical source format (`--check` in CI, gating), so diffs are semantic and review is
  reliable.
- **Aqua.jl** — `julia --project=test/quality test/quality/aqua.jl` audits package
  hygiene; the compat bounds, undefined exports, method ambiguities and unbound type
  parameters gate the build. A clean Aqua run is a precondition for trusting that the
  public API surface (User Guide chapter 16) is consistent.
- **JET.jl** — `julia --project=test/quality test/quality/jet.jl` runs the whole-package
  analysis restricted to the package's own modules, which the script gates to zero
  reports (the CI job is advisory). The test suite's type analysis
  (`test/type_analysis/`, run in a fresh process by `Pkg.test()`) applies JET's
  `@test_call` to `residual!` and the colored `jacobian!` at the concrete types `setup`
  builds (a fixed mesh; a moving mesh with in-depth radiation and Darcy–Fick transport; a
  moving mesh under the equidistribution correction; both energy forms). A type
  instability or a dynamic dispatch in the hot paths would slow the solve and threaten the
  dual-number evaluation the colored Jacobian relies on.
- **The code scan** — `test/code_scan/` asserts that no file of the source, test,
  example, script, extension, benchmark and workflow trees names an identifier from its
  list of names the package does not define (`removed_names.txt`), so a stale name in a
  comment, a docstring or an unexercised script cannot go unnoticed.

### 18.9.2 Test layout

The test suite mirrors the source module tree (`test/meshes/`, `states/`, `kernels/`,
`closures/`, `assembly/`, `deactivation/`, `jacobian/`, `outputs/`, `integration/`,
`problem/`, `multid/`, and the shared-layer `boundary_conditions/`, `materials/`,
`physics/`), and `test/runtests.jl` walks the directories, so a new `test_*.jl` file is
picked up without registration; each file runs in its own module. The full suite runs
through `Pkg.test()`. The verification machinery of this chapter has direct coverage:

- Jacobian verification: `test/jacobian/test_containment.jl` (the colored Jacobian
  against the dense reference over the configurations, states and regimes of §18.3.3,
  and the dense reference against central differences), `test_pattern.jl` (the per-entry
  criterion and the color counts), `test_allocations.jl` (zero bytes per residual and per
  Jacobian), the kernels' containment testsets, and `test/deactivation/test_deactivation.jl`
  for the masked rows (§18.3.5).
- Conservation: `test/outputs/test_conservation.jl` (the recorder and the audit on
  gasifying, deactivating and shrinking columns), `test/outputs/test_energy_ledger.jl`
  (the energy row term by term and the ledger's exactness in both energy forms),
  `test/assembly/test_residual_jacobian.jl` (the antisymmetric scatter), and
  `test/boundary_conditions/test_boundary_temperature.jl` (no gas enthalpy in the surface
  balance).
- The verification suite: `test/verification/runcases.jl` runs every case of
  `test/verification/cases/` in assert mode (§18.1.3).
- Radiation: `test/kernels/test_radiation.jl` (the beam's energy, the absorptivity
  paths, the moving mesh, planar coordinates only) and
  `test/kernels/test_boundary_conditions.jl` (the time-varying incident flux).
- Determinism: `test/integration/test_policies.jl` and `test/closures/test_correction.jl`
  require repeated solves to be bit-identical under each linear solver, with in-depth
  radiation and with the correction and its events.

### 18.9.3 Determinism and environment

Reproducible numbers require a pinned `Manifest.toml` (the integrator, linear-solver and
automatic-differentiation package versions all affect the last digits) and a fixed choice
of the integrator policy, the tolerances and the linear solver. The defaults — the
`Rodas5P` policy (`default_policy`), the `StructuredKLU()` linear solver, the colored
Jacobian (the only one), `abstol = 1e-8`, `reltol = 1e-6`, `energy_state = :temperature`,
`threshold = 1e-2` and located events — are the reference configuration. `integrate`
resets the mutable pieces of a set-up problem at its start (the layout's mask, the
correction's layer times, the enthalpy form's warm-start temperatures), so repeated
integrations of one set-up problem are bit-identical, and the determinism tests hold
repeated solves to that under `StructuredKLU()` and `UMFPACKFactorization()`. The
alternatives (`fbdf_policy`, `kencarp4_policy`, `UMFPACKFactorization()`, the enthalpy
form) trade cost against robustness on particular problems and should be re-checked
against the reference configuration on a new problem when first adopted: the same QoIs
at both, a ledger audit (§18.4), and, if a modification of the residual is involved, the
Jacobian check of §18.3.4.

---

## 18.10 Limitations and open issues

The V&V coverage above has gaps that a user should weigh:

1. **MMS coverage.** The manufactured-solution driver (§18.2) covers the coupled residual
   on a fixed mesh (conduction with temperature-dependent properties, two-step kinetics
   with reaction heat, Fickian gas transport) and the moving mesh under prescribed and
   reaction-driven motion. In-depth radiation, Darcy–Fick transport, the variable cross
   section and deactivation are verified by exact solutions and ledgers (V-RAD-01,
   V-GT-02, V-CPL-03, V-BO-04) rather than by a manufactured order study, and the boundary
   faces enter the manufactured cases as imposed fluxes only.

2. **Convergence reference is self-consistent, not exact, for the full cases.** The
   convergence studies of the validation cases use the finest mesh as the reference
   (Richardson self-convergence); they confirm *grid convergence* and its *rate*, but not
   the absolute error against a true solution. For QoIs with a known exact value (the
   verification suite) the absolute error is bounded; for the full nonlinear cone case it
   is not.

3. **The burnout envelope is bounded by its decks** (§18.5.3, §18.8.2): the deactivation
   numerics are verified under joint refinement on five PMMA-on-Kaowool decks and the
   V-BO-01 flux and kinetics ladder; other materials, thicker or charring layers and
   higher heating are outside it by that study, and outside the envelope the end-of-layer
   burst need not converge.

4. **Sensitivities across events** (§18.8.3): finite differences remain valid across a
   located deactivation, but objectives sampled at fixed times are rough at small steps,
   and a perturbation that changes an event's decision measures a jump.

5. **Validation is parameter-entangled.** Several cone and inverse-analysis cases were
   calibrated against the data they are compared to, and those calibrations need
   re-fitting under the surface convention of §18.8.1. Re-fitting against held-out data
   (calibration and validation data disjoint) is the outstanding validation task, and is
   the discipline §18.1.1 requires for a comparison to count as validation rather than
   calibration.

6. **Heuristic closures carry validation debt** (§18.8.4): the Lennard-Jones diffusivity
   estimates and the empirical shrinkage laws are unvalidated or simplified, and cases that
   lean on them inherit that uncertainty.

7. **Cross-code coverage.** ThermaKin2Ds is not among the overlays (its distribution is a
   Windows console binary), and the MaCFP-2 PMMA community submissions, committed under
   `test/verification/reference/macfp/`, are read by no case.
