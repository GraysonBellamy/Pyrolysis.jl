# 18. Verification and Validation

This chapter sets out the verification and validation (V&V) methodology for
Pyrolysis.jl. Following Roache's now-standard distinction, *verification* asks
whether the discrete equations are solved correctly (a mathematics question,
internal to the code) and *validation* asks whether the right equations are
solved (a physics question, answered against experiment and trusted reference
codes). We catalogue the analytical and manufactured-solution checks, the
Jacobian-against-reference machinery, the runtime conservation ledgers used as
continuous self-checks, the spatial/temporal/tolerance convergence studies, the
experimental and code-to-code validation cases, the regime-validity diagnostics
that tell the user when a modeling approximation is being stressed, and the known
discrepancies that still require attention. Every tolerance, sign convention, and
discrete formula stated here has been traced to the implementing source; where a
formula appears it is the one the code evaluates, not an idealized textbook form.

All symbols conform to the unified nomenclature of §2; the three binding sign
conventions are reused throughout: (i) positive flux is transport in the $+z$
direction (bottom→top); (ii) the discrete divergence $(\mathrm F_R-\mathrm
F_L)\,A/V_i$ is positive for net outflow; (iii) heat of reaction $h>0$ is
endothermic, with the volumetric source $Q_\mathrm{rxn}=-\sum_r h_r r_r$.

---

## 18.1 The V&V philosophy

### 18.1.1 Verification versus validation

A pyrolysis solver is a chain of approximations: a continuum PDE system (§3) is
closed by constitutive laws (§4–§10), projected onto a finite-volume mesh (§13),
advanced by a stiff implicit integrator with a structured Jacobian (§15), and
read out through post-processing (§16). Errors can enter at every link, and they
have qualitatively different remedies. We therefore separate the two questions
explicitly.

- **Verification** — *Are we solving the equations right?* This is purely
  internal. It compares the code's discrete output against (a) closed-form
  solutions of the *same* discrete or continuous equations, (b) the
  finite-difference / automatic-differentiation reference of its own Jacobian,
  and (c) the conservation laws the discretization is designed to satisfy. A
  verification failure is a *bug* (a wrong sign, a missing factor, a stale cache,
  an incorrect limiter), not a modeling deficiency. Verification can, in
  principle, be made exact: a correctly implemented conservative scheme conserves
  mass to round-off, a correctly assembled Jacobian matches its AD reference to
  AD precision, and cell volumes recomputed from the node positions satisfy the
  geometric conservation law as an identity.

- **Validation** — *Are we solving the right equations?* This is external. It
  compares the code's predictions against measured cone-calorimeter mass-loss
  and surface-temperature histories, against thermogravimetric (TGA) curves, and
  against the published output of established codes (ThermaKin2Ds, Gpyro, the FDS
  solid phase). A validation discrepancy can be a modeling choice
  (one-temperature gas, gray-body radiation, single-reactant kinetics), a
  parameter-calibration issue, or a genuine missing physics. Validation is never
  exact: it is bounded by measurement uncertainty and by the irreducible
  approximations of a 1D condensed-phase model.

The two are kept methodologically distinct because conflating them is the most
common way to draw a wrong conclusion. A model that agrees beautifully with a
TGA curve after the activation energy has been tuned to it has been *calibrated*,
not validated; conversely, a perfectly verified code can still disagree with
experiment because its physics is incomplete. Throughout this chapter we keep the
ledger checks (verification) separate from the experimental comparisons
(validation), and we report the regime diagnostics (§18.7) that quantify how hard
a given run leans on an approximation.

### 18.1.2 The verification hierarchy in Pyrolysis.jl

Pyrolysis.jl layers its verification from cheap, fundamental checks to experimental
validation; of these, the Jacobian comparison can also be armed during a
production solve (§18.3.4):

1. **Unit-level analytical checks** (§18.2): isolated kernels — conduction,
   single-reaction kinetics, Beer–Lambert attenuation, the GCL, the P1 sign —
   compared against hand-derivable solutions.
2. **Jacobian verification** (§18.3): the production structured Jacobian compared
   against a whole-residual ForwardDiff reference, either as a one-shot call or
   as an in-solve callback.
3. **Conservation ledgers** (§18.4): the mass and energy trackers of §16, driven
   by the tests and by the verification cases that carry their own per-step
   ledgers; in a production solve `diagnostics = true` logs sensible-energy
   snapshots (§16.8.3).
4. **Convergence studies** (§18.5): spatial, temporal, and integrator-tolerance
   refinement, confirming that the discretization error shrinks at the designed
   rate.
5. **Code- and experiment-level validation** (§18.6): cone calorimeter, Douglas
   fir, pressure-treated wood, moisture transport, against measurements and
   reference-code output.

The first three run in the test suite (the verification suite of §18.1.3 also
checks convergence orders), while the convergence studies of §18.5 and the validation cases are scripted
examples shipped in `examples/`.

### 18.1.3 The verification suite: test mode and report mode

The exact-solution and manufactured-solution cases (catalogued below; the MMS
driver is described in §18.2) form a
verification suite that runs in two modes from a single case definition. The
whole suite has to stay cheap enough for `Pkg.test()`, while publishable
evidence — fine meshes, resolution ladders, figures — is far too expensive for it
and needs plotting dependencies the test environment should not carry.

- **Test (assert) mode** runs inside `Pkg.test()`. `test/runtests.jl` picks up the
  `.jl` files directly under `test/verification/` (non-recursively), and
  `runcases.jl` among them is the suite's entry point: it runs every case at its CI resolution, asserts each quantity of interest (QoI)
  against a pinned tolerance, runs a short spatial-convergence check where the
  case defines one, and asserts the wall time (after a discarded warm-up solve)
  against five times the case's budget. It plots nothing and writes no files.
  Setting `PYROLYSIS_VERIFICATION_PRINT_QOIS=1` prints every observed QoI and
  convergence rate, which is how tolerances are re-pinned.
- **Report mode** runs in the isolated environment `scripts/verification/` (its
  own `Project.toml` with CairoMakie, so the test manifest stays free of plotting
  dependencies; the same pattern as `scripts/profile/` and `test/quality/`):

  ```
  julia --project=scripts/verification scripts/verification/run_reports.jl [--cases=… --outdir=… --chapters=…]
  ```

  `run_reports.jl` runs each case over its report resolution ladder with its
  report solver settings and writes the artifact set under
  `<outdir>/assets/<case id>/`: overlay figures (solver values as markers over
  densely sampled exact curves), pointwise-error figures that show *where* the
  discretization struggles, a log–log convergence figure with the fitted observed
  order and an expected-order reference slope, and CSVs of every curve, the
  resolution ladder, the QoIs at the finest rung, the ladder QoI errors, and the
  solver configuration (`report_harness.jl`, with all plotting in `render.jl`).
  `chapters.jl` then writes one Markdown chapter per case family — problem
  statement, figures, error and order tables — plus an index dashboard. Report
  mode does not assert the QoIs (a failed solve still aborts the run): coarse ladder rungs may exceed the pinned tolerances,
  and the tables record within-tolerance yes/no instead. The default output
  directory is the committed `docs/guide/verification/`; the documentation build
  copies it into the published site, and the manual `verification-report`
  workflow regenerates it and opens a pull request with the diff.

**The shared case definition.** `test/verification/casedefs.jl` declares the
types both modes consume and depends on nothing from `Test`, so the report runner
can include it from its isolated environment. Each file in
`test/verification/cases/` defines a `VerificationCase` with

- identity and provenance: `id`, `title`, physics `tags`, `references`;
- `build(n_cells) -> PyrolysisProblem`, optional `solve_kwargs`, a `solve_fn`
  when the case solves through its own path (the MMS cases use the wrapped-RHS
  `solve_mms`), and `expected_retcode` (`:Success` by default, `:Terminated`
  when a callback stopping the solve is the answer);
- `qois(sol) -> Vector{QoIResult}`: each QoI's value, exact value, error, and
  tolerance, with a `comment` recording the tolerance provenance (tolerances are
  pinned at about three times the error observed when pinned);
- `ci::CIConfig`: the CI resolution, the wall-time budget, and an optional
  `ConvergenceSpec` (resolutions, error functional, expected and minimum observed
  order);
- `report::NamedTuple`: the report-mode keys — the resolution ladder
  `n_cells_ladder`, `curves(sol) -> Vector{ReportCurve}` (plain data:
  `CurveSeries` samples, each with an optional `exact` callable on the same
  abscissa from which the pointwise error is derived), `error_norms`,
  `expected_order`, report `solve_kwargs`, and optional temporal or parameter
  sweeps. A case with `report = (;)` is skipped in report mode.

The case file's leading comment banner is the problem statement the report
chapter prints, so the published chapter cannot drift from what the test runs.
Reference solutions are computed in the run at high precision (series summed to
machine tolerance, transcendental roots by bisection), never hard-coded; they
live in the shared library `test/verification/exact/` (series, roots, the MMS
driver), which carries its own self-tests.

**Case catalogue.** Each case family has one file in `test/verification/cases/`;
lettered members such as V-HC-02a–d share the family file and its report chapter.
Most families compare against an exact or manufactured solution; the burnout
family and V-CPL-06 check budgets and a pinned reference run instead.

| Family | Cases | Reference solution | What it verifies |
|---|---|---|---|
| Heat conduction | V-HC-01 – V-HC-05 | Carslaw–Jaeger erfc and eigenfunction series, resistance network, Kirchhoff transform | the convective-BC surface Newton (V-HC-01; Biot sweep V-HC-02), the constant-flux slab to its back-face plateau (V-HC-03), the harmonic-mean face conductivity at a layer jump (V-HC-04), $k(T)$ (V-HC-05) |
| Kinetics and surface radiation | V-RX-01 – V-RX-04 | Bateman chain, Khawam–Flanagan isothermal integrals, parallel-reaction splits, implicit radiative-cooling law | multi-step coupling, reaction order and the depletion limiter, two-product yields and heat superposition, `RadiativeBC` through the surface Newton |
| In-depth radiation | V-RAD-01 | Lautenberger eigenfunction series | Beer–Lambert absorption with conduction and a convective face, including the sub-surface temperature peak |
| Gas transport | V-GT-01, V-GT-02 | Gaussian plane source; steady compressible Darcy $P^2$ profile | the Fickian operator and mass conservation; the pressure closure, Darcy mobility and outflow BC |
| Coupled physics | V-CPL-01 – V-CPL-06 | Neumann two-phase Stefan solution, Landau ablation wave, exact shrink/swell endpoints, Frank-Kamenetskii slab, 0D TGA curve, discrete gas-enthalpy ledger | a latent-heat front, ALE surface regression, volume change and the lateral-shrinkage law, the reaction–conduction bifurcation, the BC/volume coupling of kinetics, $S_\mathrm{conv}$ (§18.4.3) |
| Manufactured solutions | V-MMS-01 – V-MMS-04 | the MMS driver of §18.2 | design order of the coupled Eulerian residual in space and time, free-stream preservation under prescribed motion, ALE order under prescribed and reaction-driven motion |
| Burnout | V-BO-01 – V-BO-05 | exact exponential thickness (V-BO-04); otherwise mass budgets and a pinned reference run | parked-cell depletion through complete burnout on ALE and fixed meshes (§14.3), the transmitted-flux warning, the interior-void stop (`expected_retcode = :Terminated`, V-BO-05) |

Kernel-level analytic checks are in `test/verification/test_analytical.jl`, which
`Pkg.test()` runs beside the suite: semi-infinite conduction, the Beer–Lambert
profile and energy split, isothermal and constant-heating-rate kinetics, and the
mass and energy ledgers. Two conventions matter when a case transcribes
literature parameters: $h>0$ is endothermic here, the opposite of ThermaKin,
Gpyro and FDS; and most exact-kinetics cases disable the depletion limiter
(`DepletionLimiter(enabled = false)`) and open the temperature gates (§18.2.2),
while V-RX-02d keeps the limiter on as a guard.

**Cross-code overlays.** Where a case has a published analog in another code,
that code's numerical output is committed under
`test/verification/reference/<code>/` with a provenance README (FDS and Gpyro runs
on the same problems; MaCFP-2 PMMA community submissions, which no case reads) and
read only inside a case's
report `curves`. A series marked `crosscode = true` is drawn as a muted context
layer in a separate `_vs_codes` variant of each affected figure, so the base
figure keeps the solver's own error on its own scale; for exact-solution cases the
other codes' *errors* against the same exact solution are compared.
`reference/timings.csv` feeds a wall-time comparison at every ladder rung, and a
case's ladder includes a rung matched to the reference run's resolution where
needed.

---

## 18.2 Method of manufactured solutions and analytical checks

The most rigorous verification tool is the *method of manufactured solutions*
(MMS): one chooses an analytic field, substitutes it into the governing PDE to
obtain the exact forcing, supplies that forcing to the discrete solver, and
confirms that the numerical solution converges to the chosen field at the
designed order. The verification suite ships a coupled MMS driver
(`test/verification/exact/mms.jl`: `MMSFields`, a ForwardDiff exact-source
generator, manufactured-BC builders, and the wrapped-RHS `solve_mms` path)
exercised by the V-MMS cases in `test/verification/cases/` — coupled Eulerian
MMS, ALE free-stream/GCL preservation, ALE MMS under prescribed mesh motion via
the `Workspace.w_override` hook, and ALE MMS under reaction-driven motion
(V-MMS-04, `MMSMotion`). Two constraints of that driver
follow from the formulation: the structured geometry Jacobian assumes
reaction-driven motion, so `w_override` runs use finite-difference Jacobians;
and because the scheme is a **solid-Lagrangian ALE** (mesh velocity ≡ solid
material velocity; no component carries a discrete advection operator —
§11.1.2, §11.9 item 5), every manufactured species source on a moving mesh
is assembled in the material-following frame, with the prescribed or the
reaction-driven mesh velocity as its $w$.

> **What MMS cannot see: wrong-PDE defects.** MMS verifies that the
> discretization solves the code's *own* target PDE at the designed order —
> the manufactured source is built from that same PDE, so any term the
> formulation wrongly includes (or omits) is self-consistently compensated
> in the source and invisible at every order and every resolution. A
> spurious ALE advection term, for example, would pass every MMS case
> untouched. Only checks with an *independent physical identity* the code
> does not get to define — an exact solution (the Landau traveling wave of
> V-CPL-02), a conservation ledger (V-CPL-03, V-CPL-06) — can catch a
> wrong-PDE defect. This is why the suite pairs MMS (discretization order)
> with exact-solution and conservation anchors (formulation correctness):
> the two families have disjoint blind spots.

In addition, Pyrolysis.jl supports the *analytical-solution* special case of
MMS for each isolated physics kernel, in which a parameter regime is chosen
so that an exact solution of the governing equation is already known and no
manufactured forcing is needed. These are the building blocks; each is
independently verifiable and each anchors a specific term of `compute_rhs!`.

### 18.2.1 Pure conduction

With all gas transport, reaction, and radiation terms disabled, the energy
equation of §7 collapses to the linear heat equation

$$
\rho c_p^\mathrm{eff}\,\frac{\partial T}{\partial t}
   = -\nabla\!\cdot\mathbf q_\mathrm{cond},
\qquad
\mathbf q_\mathrm{cond} = -k\,\nabla T .
$$

For constant $\rho c_p^\mathrm{eff}$ and $k$ on a fixed (Eulerian) mesh this is
$\partial_t T = \alpha_\mathrm{th}\,\partial_{zz}T$ with thermal diffusivity
$\alpha_\mathrm{th}=k/\rho c_p^\mathrm{eff}$, which admits a wealth of
closed-form solutions: the error-function semi-infinite response to a step
surface temperature, the Fourier-series response to a constant surface flux, and
the decaying-cosine eigenmode

$$
T(z,t) = T_\infty + (T_0-T_\infty)\,
   \cos\!\Big(\tfrac{\pi z}{2L}\Big)\,
   \exp\!\Big(-\alpha_\mathrm{th}\,\tfrac{\pi^2}{4L^2}\,t\Big)
$$

on $[0,L]$ with an insulated bottom ($\partial_z T|_{z=0}=0$, an `AdiabaticBC`)
and a fixed top ($T(L,t)=T_\infty$, a `DirichletThermalBC`). The single decay
constant $\lambda=\alpha_\mathrm{th}\pi^2/4L^2$ makes this an exacting check on
the whole conduction pathway: the harmonic-mean face conductivity, the divergence
operator, the diffusivity assembly, and the implicit time integration must all be
correct for the numerical decay rate to match.

The discrete face flux verified here is (§7.1, `conductive_flux`,
`src/Physics/heat_conduction.jl`)

$$
q_f = -k_f\,\frac{T_R-T_L}{d_L+d_R},
\qquad
k_f = \frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
$$

(distance-weighted harmonic mean). On a *uniform* mesh, $d_L=d_R=\tfrac12\Delta
z$ and $k_f$ reduces to the ordinary harmonic mean $2k_Lk_R/(k_L+k_R)$, and for
constant $k$ the scheme is the textbook second-order central difference; the
numerical decay rate then matches $\lambda$ to $O(\Delta z^2)+O(\Delta t^p)$ where
$p$ is the integrator order. The harmonic mean (rather than arithmetic) is
verified by a *two-layer* slab with a conductivity jump: the steady-state
solution has a continuous flux and a kinked temperature profile, and only the
harmonic (resistance-in-series) face value reproduces the analytic interface
temperature.

**Implementation note.** Interior conductive fluxes are assembled by
`compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`), which zeroes
boundary faces (the BC solver overwrites them); the divergence
$(\nabla\!\cdot\mathbf q)_i=(q_{f_R}-q_{f_L})A/V_i$ is `divergence_1d!`
(`src/Discretization/finite_volume.jl`). The RHS assembles
$dT_i/dt = -(\nabla\!\cdot\mathbf q)_i/\rho c_p^\mathrm{eff}_{,i}$ in step 8 of
`compute_rhs!`.

### 18.2.2 Single-reaction TGA (0D-idealized cell)

A thermogravimetric experiment heats a small, thermally-thin sample at a constant
rate $\beta = dT/dt$ and records the residual mass fraction. In the limit of a
single, well-mixed cell with a uniform temperature ramp, all transport vanishes
and the species balance for a one-step first-order decomposition $A_1\!\to\!
\text{products}$ becomes the 0D ODE

$$
\frac{d\xi_1}{dt}
   = -A_1\exp\!\Big(-\frac{E_1}{R_g T}\Big)\,\xi_1 ,
\qquad T(t)=T_0+\beta t .
$$

Changing the independent variable to temperature and integrating gives the
classic non-isothermal conversion law

$$
\frac{\xi_1(T)}{\xi_1(0)}
   = \exp\!\left[
       -\frac{A_1}{\beta}\int_{T_0}^{T}
         \exp\!\Big(-\frac{E_1}{R_g\,T'}\Big)\,dT'
     \right],
$$

whose integral is the Arrhenius temperature integral (approximated by the
Coats–Redfern or Senum–Yang series, or evaluated numerically). The conversion
$\chi_\mathrm{conv}=1-\xi_1/\xi_1(0)$ traces the familiar sigmoid; its inflection
(maximum reaction rate) sits at the temperature where
$d^2\chi_\mathrm{conv}/dT^2=0$, giving the reference relation between the peak
temperature $T_p$ and the kinetic pair $(A_1,E_1)$,

$$
\frac{E_1\,\beta}{R_g\,T_p^2}
  = A_1\exp\!\Big(-\frac{E_1}{R_g\,T_p}\Big),
$$

which is exactly the FDS-style recipe for recovering $(A_1,E_1)$ from a measured
$T_p$ and peak rate (see §6.9 and §18.6.2). Reproducing $T_p$, the peak rate, and
the residual fraction in a single-cell Pyrolysis run verifies the Arrhenius
kernel, the species stoichiometry, and the depletion limiter together.

Two implementation subtleties must be accounted for when this check is made
*quantitatively* exact rather than merely qualitative:

- **Mass-concentration basis.** The rate law is evaluated on mass concentration
  $\xi_j$ (kg m$^{-3}$), not mole fraction or normalized conversion:
  $r_i=A_i\exp(-E_i/R_gT)\prod_j\xi_j^{\,n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ (Eq. 1 of §6,
  `reaction_rate`, `src/Physics/kinetics.jl`). For a first-order single
  reactant this is linear in $\xi_1$ and the analytic solution above holds with
  $\xi_1(0)$ the *initial concentration*.
- **Smooth gates and the depletion limiter perturb the tails.** The hard
  Heaviside gates $H(T-T_\mathrm{min})H(T_\mathrm{max}-T)$ are replaced by tanh
  ramps $\tfrac12(1+\tanh(T-T_\mathrm{min}))\cdot\tfrac12(1+\tanh(T_\mathrm{max}
  -T))$ (§6.2, `kinetics.jl`), which leak $\sim\!10$ K outside the nominal
  window; and the depletion limiter multiplies the rate by
  $\tanh(\xi_1/\text{threshold})$ (default threshold $1.0$ kg m$^{-3}$, §6.3,
  `kinetics.jl`). To compare against the bare-Arrhenius analytic solution one
  must either set $T_\mathrm{min},T_\mathrm{max}$ far outside the active range
  (so the gates read $\approx 1$) and keep $\xi_1$ well above the limiter
  threshold throughout, or disable the limiter via
  `DepletionLimiter(enabled=false)`. With those settings the numerical curve
  matches the analytic conversion to integrator tolerance; the residual gap *is*
  the smoothing error and is itself a quantitative verification of the gate
  width.

### 18.2.3 Beer–Lambert transmitted fraction

The in-depth radiation model (§8.3) integrates the Beer–Lambert law
$\partial I/\partial z=-\alpha_\mathrm{eff}I$, giving
$I(z)=I_\mathrm{surface}\exp\!\big(-\!\int_z^L\alpha_\mathrm{eff}\,dz'\big)$ and
the exact cell-averaged volumetric source

$$
q'''_i = I_\mathrm{in}\,\frac{1-\exp(-\tau_i)}{\Delta z_i},
\qquad
\tau_i = \alpha_i\,\Delta z_i,
$$

assembled by `compute_beer_lambert_source!` in a top-to-bottom sweep
($i=n\!:\!-1\!:\!1$) so that absorbed flux upstream correctly attenuates the
intensity reaching downstream cells (`src/Physics/radiation.jl`). Two exact
energy-balance identities verify this kernel:

- **Total absorbed flux.** Summing the volumetric source times cell thickness,

  $$
  \sum_i q'''_i\,\Delta z_i
     = \alpha_\mathrm{surf}\,I_\mathrm{surface}\big(1-e^{-\tau_\mathrm{total}}\big),
  \qquad
  \tau_\mathrm{total}=\sum_i\alpha_i\Delta z_i,
  $$

  where $\alpha_\mathrm{surf}$ is the surface transmissivity of §8.3.2
  ($\varepsilon_{\text{eff}}$ on the Kirchhoff path; 1 when the BC's explicit
  absorptivity is already folded into the intensity), is a telescoping sum that
  is *exact at the discrete level* — independent of cell count — because each
  cell deposits $I_\mathrm{in}(1-e^{-\tau_i})$ and passes
  $I_\mathrm{in}e^{-\tau_i}$ to the next. The solution reports it as
  `extras.volumetric_absorption`, and the sum of the deposited source over the
  cells must equal this analytic value.

- **Transmitted fraction.** The flux escaping the back face is
  $\alpha_\mathrm{surf}\,I_\mathrm{surface}\,e^{-\tau_\mathrm{total}}$, reported as
  `extras.radiation_transmitted` (§8.3.2). Conservation requires
  $\text{absorbed}+\text{transmitted}=\alpha_\mathrm{surf}\,I_\mathrm{surface}$
  to round-off — the remaining $(1-\alpha_\mathrm{surf})\,I_\mathrm{surface}$ is
  reflected at the surface and never enters the medium — which is the strongest
  single check on the Beer–Lambert path.

The transparent-cell guard ($\alpha>10^{-20}$, else $q'''_i=0$,
`radiation.jl`) is verified by the $\alpha\to0$ limit, where $1-e^{-\tau}\to
0$ and all flux is transmitted. The opposite limit $\tau_1\gg1$ in the surface
cell deposits the entire entering flux
$\alpha_\mathrm{surf}\,I_\mathrm{surface}$ in the top cell — consistent with
`SURFACE_ABSORPTION`, which applies the same
$\varepsilon_{\text{eff}}\,I_\mathrm{surface}$ at the face (§8.2), providing a
cross-check between the two non-P1 radiation models.

### 18.2.4 Geometric Conservation Law to machine precision

The single most decisive ALE verification is the Geometric Conservation Law
(GCL, §11.5). A moving-mesh scheme that does not reproduce the swept cell volume
exactly injects *phantom* mass and energy sources proportional to the volume
error; over hundreds of stiff steps these can swamp the genuine reaction-driven
mass loss. The continuous 1D GCL is $dV_i/dt = A\,(w_R-w_L)$ (plus the
area-change term $|z_R-z_L|\,dA/dt$ under a lateral shrinkage law). Pyrolysis.jl
satisfies it *geometrically*: the node positions are state ($dz_i/dt=w_i$) and
every residual evaluation recomputes $V_i=A(t)\,|z_R-z_L|$ from them, so the
volume the scheme uses is always the swept one and no separately evolved volume
can drift (§11.5.2). With the cumulative mesh velocity this gives
$dV_i/dt=V_i\theta_i$ identically, and with the dilation correction a uniform
field with no sources and no fluxes stays uniform.

The verification case V-MMS-02 asserts exactly that: a spatially uniform state
with no reactions and zero-flux boundaries, on a mesh compressed to 70 % of its
thickness through the `Workspace.w_override` hook with the uniform-stretch
velocity $w^*(z,t)=(\dot s/s)\,z$, must stay uniform to machine precision. The
node trajectories $z_i(t)=z_i(0)\,s(t)$ are known exactly, so the top-cell-center
position doubles as a check on the prescribed-motion path. V-MMS-03 then measures
the observed order of a smooth manufactured solution under prescribed smooth
mesh motion. A drift in V-MMS-02 is a moving-mesh transport defect that the
smooth-field order study could hide inside its truncation error.

### 18.2.5 P1 sign test

The P1 quasi-steady radiation model (§8.4, Experimental) solves
$\nabla\!\cdot(D_\mathrm{rad}\nabla G)=\alpha G-4\sigma\alpha T^4$ for the
incident field $G$ and feeds the net source

$$
S_\mathrm{rad} = \alpha\,(G-4\sigma T^4)
$$

into the energy equation. The sign is the part most easily gotten wrong, and it
is verified directly: $S_\mathrm{rad}>0$ when $G>4\sigma T^4$ (the local field
exceeds the blackbody level → net absorption → heating), $S_\mathrm{rad}<0$ when
$T^4$ dominates (net emission → cooling), and $S_\mathrm{rad}=0$ at local
equilibrium. Two corollary checks pin the boundary and the volumetric coupling:

- **Marshak equilibrium.** With no external source ($I_\mathrm{ext}=0$) and the
  surface in equilibrium with its surroundings ($T_s=T_\infty$, so
  $G_\mathrm{surf}=G_\mathrm{ext}=4\sigma T_\infty^4$), the Marshak boundary flux
  $q_\mathrm{BC}=\frac{\varepsilon}{2(2-\varepsilon)}
  (G_\mathrm{ext}+4I_\mathrm{ext}-G_\mathrm{surf})$ must vanish identically. This
  verifies the exact $\frac{\varepsilon}{2(2-\varepsilon)}$ coefficient (Modest
  Eq. 16.48) rather than the incorrect simplification $0.5\varepsilon$. A second
  regression pins the **external-flux weight**: the $I_\mathrm{ext}$ response is
  $2\varepsilon/(2-\varepsilon)\cdot I_\mathrm{ext}$ — the same half-range weight
  the ambient term carries — and a cold optically-thick absorber at
  $\varepsilon=0.9$ absorbs the continuum-P1 fraction $4c/(1+\sqrt3\,c)\approx
  0.96$ of the incident flux (a weight-1 injection would absorb only
  $\approx 0.59$).
- **No volumetric emissivity factor.** The volumetric emission is the blackbody
  value $4\sigma\alpha T^4$ with *no* $\varepsilon$ factor (gray-body local
  equilibrium, $\varepsilon=1$ in the interior); $\varepsilon<1$ enters *only* at
  the boundary through the Marshak term. A test that compares interior emission
  against $4\sigma\alpha T^4$ verifies this convention.

The Newton solve for $G$ uses an allocation-free Thomas (TDMA) factorization with
relative convergence $\|R\|_\infty<\text{tol}\cdot\|G\|_\infty$ (G fields are
$10^4$–$10^6$ W m$^{-2}$, so absolute tolerance would be unphysical) and a
non-negativity clip $G\ge0$; these are verified against analytic optically-thick
and optically-thin limits in `test/radiation/test_p1_validation.jl`. Because P1
lives in `Pyrolysis.Experimental` and uses a *separate* residual
(`RadiationP1`), the production `solve` does not route
`radiation_model=P1_QUASI_STEADY` through the unified residual path; it rejects
the setting outright at keyword validation with an error directing the user to
the `Pyrolysis.Experimental.RadiationP1` primitives, which are invoked directly
(`_validate_solve_kwargs`, src/Solver/solve.jl). P1's verification is therefore confined to
the Experimental test set (§8.12).

---

## 18.3 Jacobian verification

### 18.3.1 Why and against what

The implicit integrators used for these stiff systems (default `KenCarp4`, an
SDIRK method) solve a Newton system at every stage with the iteration matrix
$W=I/(\gamma_\mathrm{RK}\Delta t)-J$ (§15.5). A wrong Jacobian $J$ does not change
the *fixed point* of the Newton iteration — so a converged step is still a
correct solution of the discrete ODE — but it degrades or destroys convergence,
inflating step rejections and, at worst, stalling the solve on an Arrhenius
transient. The production Jacobian is the *structured* operator assembly of §15
(a sparse block $A$ plus O($N$) prefix/suffix/rank-1 couplings); its analytical
and AD-per-operator pieces are optimizations of a quantity that must equal the
true Jacobian of the residual. The reference against which "true" is defined is
`DenseAD_ForwardDiff`: a single ForwardDiff pass over the *entire* residual
vector, which is correct by construction (it differentiates exactly the function
the integrator evaluates) and slow (dense, full-width). Verification compares the
two.

### 18.3.2 `verify_jacobian_at`

The standalone comparator is

```julia
verify_jacobian_at(spec_active, def, ws, u, t;
    reference = DenseAD_ForwardDiff(), rtol = 1e-6, atol = 1e-10) -> NamedTuple
```

(`src/Jacobian/verify.jl`). It assembles the active backend's Jacobian and the
reference Jacobian at the same state $(u,t)$ and reports the worst absolute and
relative deviation and the $(\text{row},\text{col})$ where it occurs.

### 18.3.3 Strict versus masked comparison

The comparison mode is selected by the `is_approximate` trait of the active
backend (§15.6, §15.10):

- **Strict / support-complete** (exact backends, `is_approximate == false`):
  every entry of the dense reference is compared against the active backend's
  dense view, with the backend's view zero outside its declared support, under
  the **per-entry** criterion

  $$
  |\Delta J_{ij}| \;\le\; \texttt{atol} + \texttt{rtol}\cdot
  \max\big(|J^{\rm ref}_{ij}|,\ \max\nolimits_j |J^{\rm ref}_{i,:}|\big).
  $$

  The row-magnitude floor keeps honest roundoff at structural zeros of a
  stiff row from failing (a $10^{-9}$ residue in a $10^6$-scale Arrhenius row
  passes), while a wholly missing block in a *small*-magnitude row — the ALE
  z-rows are ~10 orders below the ξ-rows — contributes to its own row's
  reference scale and fails regardless of problem size. (A global
  max-$|J^{\rm ref}|$ normalization would make pass/fail scale-dependent —
  small-magnitude rows such as the ALE z-rows could lose an entire block
  behind the Arrhenius rows' magnitude — which is why the criterion is
  per-entry with a row floor.) A reference non-zero that falls outside the backend's
  support surfaces as a deviation — a *missing-entry* failure. The verifier
  separately reports `max_abs_dev_in_support` (genuine disagreement on a
  claimed entry) and `max_abs_dev_outside_support` (magnitude of reference
  entries the support omits). The exact backend, `Structured(DirectSolve())`,
  must pass this regime.

- **Masked** (approximate backends, `is_approximate == true`): entries outside
  the active backend's support are zeroed in the difference, so only in-support
  disagreement is asserted, under the coarse global-scale criterion
  $\texttt{max\_abs\_dev}\le\texttt{atol} \lor
  \texttt{max\_abs\_dev}/\max_{ij}|J^{\rm ref}_{ij}|\le\texttt{rtol}$. This is
  the regime for `Structured(ApproxSolve())`, which *deliberately* drops the
  structured ALE couplings (the cumulative mesh-velocity prefix sum, the z-row
  prefix blocks, the Beer–Lambert suffix, the $\theta_A$ rank-1 back-reaction,
  the geometry $z$-columns) under the `:ale_reduced` sparsity profile. Its
  in-support values legitimately deviate by the dropped coupling magnitudes,
  so masked comparison is a sanity check
  rather than an exactness gate; the dropped couplings are by design (they
  cost accuracy in the Newton step, not in the converged solution).

### 18.3.4 In-solve verification callback

The same comparator can be installed as a `DiscreteCallback` that fires when the
integrator crosses requested sample times:

```julia
solve(problem; verify_jacobian = (rtol=1e-6, atol=1e-10, sample_times=[0.0, 1.0, 10.0]))
```

(`src/Solver/solve.jl`). The callback compares the active backend against the
`DenseAD_ForwardDiff` reference at each sample time and `@warn`s on a failure but
*does not throw*: a single transient outlier (e.g. an extreme Arrhenius spike
where the analytical and AD paths differ at the last AD digit) should not abort
an otherwise-good production run. This lets a long calibration sweep be run with
verification armed at negligible extra cost — the reference Jacobian is built
only at the sample times.

### 18.3.5 Parked states

Depletion handling (§14.3) parks cells without changing the state size, so the
structured Jacobian is the one the integrator uses through burnout and the
verifier covers it: `test/jacobian/test_park_jacobian.jl` builds the specs on
the fresh problem, parks surface, interior and bottom cells in every mode
(Eulerian, ALE, Darcy, Beer–Lambert, variable cross-section), and checks
`Structured(strategy = DirectSolve())` strictly against the dense reference, the
zero ghost rows, the ghost node rows equal to the surface node row, and the
live block of the parked Jacobian against the compacted problem's — for
`ApproxSolve` too, since that strategy drops couplings by design and is not
compared against the dense reference on a reacting state.

---

## 18.4 Conservation-based verification

The conservation ledgers of §16 serve as verification where the tests and the
verification cases drive them per step (`record_substep!`, or a case's own
per-step ledger callback); a production `solve(problem; diagnostics=true)` logs
only sensible-energy snapshots, with the ledger columns NaN (§16.8.3). The checks
divide cleanly into those that must close by construction (so a failure is a bug)
and those that quantify an approximation (so a large value is a regime warning,
not a bug).

### 18.4.1 Mass ledger

The per-component mass tracker (`MassConservationTracker`,
`src/Diagnostics/conservation_checks.jl`) maintains
$m_j(t)=\sum_{i\in\mathrm{active}}\xi_j^{(i)}V_i$ (masked: only active cells,
§16.3.1) and accumulates boundary efflux
$\Delta m_{\mathrm{out},j}=-(J_{\mathrm{top},j}+J_{\mathrm{bot},j})\,A\,\Delta t$
with the convention that positive $J$ is *inflow*, so positive flux *reduces*
efflux. The balance check is

$$
\max_j
  \frac{\big|m_{\mathrm{current},j}
          -(m_{\mathrm{initial},j}-m_{\mathrm{efflux},j})\big|}
       {|m_{\mathrm{initial},j}|}
  < 10^{-10}
$$

(`check_mass_conservation`, default tol $10^{-10}$). When a component starts at
zero or near-zero mass (e.g. char before any reaction), the check switches from
relative to *absolute* error to avoid division by a vanishing baseline. For ALE
runs with a time-varying cross-section the area in the efflux term is the
trapezoidal substep average $A=\tfrac12(A^\mathrm{pre}+A^\mathrm{post})$, which
keeps the contribution symmetric and consistent with the geometric volume update.
A mass-balance failure indicates a flux-assembly or callback bug.

### 18.4.2 Energy ledger

Energy verification uses the sensible-energy ledger of §16.4–16.5. The stored
quantity is $E(t)=\sum_i\rho c_p^\mathrm{eff}\,T_i\,V_i$ with the mixture
$\rho c_p^\mathrm{eff}$ (every component, gas included, §16.4), balanced against
the *exact terms the RHS assembles*:

$$
E(t) = E(0)
  + \int_0^t\!\Big[
      q_\mathrm{cond,boundary}\,A
      + Q_\mathrm{rxn}
      + S_\mathrm{conv}
    \Big]\,d\tau .
$$

The residual $r=E(t)-E(0)-\sum(\text{sources})$ should close to integrator
tolerance ($\sim10^{-6}$ relative, `check_energy_conservation` default tol
$10^{-6}$) because it is built from exactly the terms `compute_rhs!` evaluates. A
failure here is a bookkeeping bug. The one caveat (§16.7, §18.8) is that this
ledger *cannot* close exactly when composition varies in time, because it omits
the $\int T\,\partial(\rho c_p)/\partial t\,dV\,dt$ heat-capacity-change term of its
$\rho c_p\,T$ proxy density (§16.7) and the ALE volume-change term; the residual
then carries a small, composition-rate-dependent offset that is a known
incompleteness of the *audit*, not of the *solver*.

**The boundary flux.** The ledger uses $q_\mathrm{cond,boundary}$ (the conducted
flux from the surface Newton solve). The surface balance has no gas-enthalpy
term, so it equals the applied boundary flux $q_\mathrm{BC,applied}$
(the externally driven flux, which may include radiation/convection; §16.6).
Crucially, the discrete ledger
*audits bookkeeping, not double-counting*: if the RHS itself double-counts a
mechanism, the discrete ledger will still close, because it is built from the same
(buggy) terms. Detecting a double-count requires an *independent* check — an
exact solution or a ledger built from outside the RHS (§18.4.3).

### 18.4.3 Independent closure of the gas-enthalpy coupling

The gas enthalpy is carried entirely by the always-on volumetric advection source

$$
S_{\mathrm{conv},i}
  = -\sum_g c_{p,g}(T_i)\,\bar J_{g,i}\,
      \Big(\frac{\partial T}{\partial z}\Big)_i
$$

(§7.4, `compute_gas_convective_source!`), and the surface balance carries no
gas-enthalpy term (§7.4.3). A surface term $\sum_g\dot m_g\,\Delta h_g(T_s,T_\infty)$
added on top of $S_\mathrm{conv}$ would count the advected enthalpy twice, yet
the energy ledger of §18.4.2, built from the same RHS terms, would record both
and still close. The coupling is therefore pinned by a ledger built independently
of the RHS: the verification case V-CPL-06 (a char-layer column whose gas forms in
a cool band and vents through a hot layer, $h=0$ so only the sensible ledger acts)
checks

- **V-CPL-06a (Eulerian):** the discrete summation-by-parts identity
  $\int S_\mathrm{conv}\,dV\,dt = -c_g\big(\int J_\mathrm{top}T_{f,\mathrm{top}}A\,dt
  - \int\sum_i(\mathrm{gen}_i-\dot\xi_{g,i})\,T_i V_i\,dt\big)$ — the energy
  $S_\mathrm{conv}$ removes from the field equals the energy the species ledger
  says the gas picks up between formation and exit — and the physical reheat
  reference $-c_g\sum\mathrm{gen}\,(T_\mathrm{exit}-T_\mathrm{form})\,V$ to
  quasi-steady accuracy;
- **V-CPL-06b (ALE):** the same column with the fuel band collapsing; the ALE
  energy budget closes with $E=\sum\rho c_p T V$ over every component, and the
  mesh-relative $T$-advection accumulator is exactly zero (§11.2.2).

Its ledger runs in a per-accepted-step callback around the pure `residual!` with
trapezoidal time quadrature (a rectangle rule leaves an $O(\Delta t\,\dot X)$
residual that is pure quadrature). V-CPL-06a shows that $S_\mathrm{conv}$ alone
removes the whole gas reheat, so a surface gas-enthalpy term would remove it a
second time; a mesh-velocity $T$-advection term would appear in V-CPL-06b as a
sink of about 66 kJ per kg of volatiles. The absence of the surface term is
asserted directly by `test/boundary_conditions/test_boundary_temperature.jl`
(testset "Outward gas flux enters only through blowing"): under a prescribed-flux
BC, the solved surface temperature and conducted flux are identical with and
without a gas flux across the face.

---

## 18.5 Convergence studies

A scheme can be conservative and have a correct Jacobian yet still be *under-
resolved*. Convergence studies confirm that the discretization error shrinks at
the designed rate as the mesh, time step, or integrator tolerance is refined, and
they let the user pick a resolution that is converged for the quantity of
interest (QoI) at acceptable cost.

### 18.5.1 Spatial convergence (uniform and stretched)

The cell-centered FVM of §13 is second-order accurate in space for smooth
solutions: the harmonic-mean face conductivity and the central-difference
gradient are each $O(\Delta z^2)$ on a uniform mesh, and the distance-weighted
generalizations preserve formal second order on smooth, slowly-stretched meshes.
The error of a QoI $f$ should therefore scale as $f(\Delta z)-f_\mathrm{exact}
=C\,\Delta z^p+\cdots$ with $p\approx2$. With an exact solution unavailable, the
finest mesh is used as the reference (Richardson-style self-convergence): the
shipped `examples/convergence_douglas_fir.jl` solves a Douglas-fir cone case at a
sequence of cell counts (finest first = reference), then reports per-QoI relative
errors

$$
\Delta f(n) = \frac{|f(n)-f_\mathrm{ref}|}{|f_\mathrm{ref}|}\times100\%
$$

for scalar QoIs — peak mass-loss rate, time-to-ignition, final mass, peak surface
temperature — and *curve-level* $L_2$ relative errors

$$
\frac{\lVert f(n)-f_\mathrm{ref}\rVert_2}{\lVert f_\mathrm{ref}\rVert_2}\times100\%
$$

for the full MLR$(t)$, $T_s(t)$, and back-face $T(t)$ histories. Plotting
$\Delta f$ against cell count on log–log axes recovers the slope $-p$; a slope
near $-2$ confirms second-order spatial convergence, while a shallower slope
flags a first-order contamination (e.g. the first-order upwind in the Darcy–Fick
advection of §9.3, or the surface-cell radiation deposition).

Stretched meshes (`create_stretched_mesh`, geometric ratio, §13.9) are tested
separately because the boundary-refinement they provide near the heated surface
is where the steep thermal and reaction fronts live; the distance-weighted
harmonic mean is what preserves accuracy there (a 1–5% flux-accuracy gain near
boundaries versus the naive arithmetic mean, §7.1). A convergence study on a
stretched mesh verifies that the non-uniform-mesh corrections in the gradient and
face-value reconstructions (§13.6) are correct: an *incorrect* non-uniform
treatment shows up as a loss of order on the stretched run that is absent on the
uniform run.

### 18.5.2 Temporal convergence

The default integrator `KenCarp4` is a 4th-order SDIRK method; an integrator-order
study fixes the mesh, sweeps the maximum step (via `dtmax` or by tightening
tolerances), and confirms the temporal error scales at the method's design order
for smooth intervals. In practice the *Arrhenius transients* — the sharp onset of
decomposition as a reaction gate opens — dominate the temporal error and force
the adaptive controller to small steps; the tanh-smoothed gates and the
tanh depletion limiter (§6.2–6.3) are essential here precisely because they keep
the RHS and its Jacobian differentiable, so the error estimator and the Newton
solver behave. A temporal study that holds the mesh fixed and tightens `dtmax`
isolates the time-integration error from the spatial error.

### 18.5.3 Integrator-tolerance study

Because the step size is chosen adaptively to meet `abstol`/`reltol` (defaults
$10^{-8}$ and $10^{-6}$, `src/Solver/solve.jl`), the practical temporal-error
control knob is the tolerance pair. A tolerance study sweeps `reltol` (and
`abstol` in proportion) and tracks a QoI to its tolerance-converged plateau; the
convergence example uses tightened tolerances (`abstol=1e-9`, `reltol=1e-7`) for
its reference runs so that the *spatial* error dominates the comparison and is not
masked by temporal error. The energy-ledger residual (§18.4.2) provides an
independent confirmation that the tolerance is tight enough: it closes to roughly
`reltol`, so a residual that fails to shrink as `reltol` is tightened
indicates a bookkeeping issue rather than under-resolution. The three knobs —
$\Delta z$, integrator order/`dtmax`, and tolerance — are studied separately so
that the dominant error source is unambiguous; only when all three are converged
is a quoted QoI a property of the *model* rather than of the *discretization*.

---

## 18.6 Validation against experiments and reference codes

Validation compares converged Pyrolysis.jl predictions against measured data and
against the published output of established condensed-phase pyrolysis codes. The
notation of the reference codes is mapped to ours on first use; the comparisons
draw on the governing-form, kinetics, transport, and radiation contrasts
documented throughout §3–§9.

### 18.6.1 Cone calorimeter: HRR/MLR/$T_s$

The cone calorimeter is the canonical bench-scale validation experiment: a sample
is exposed to a constant radiant flux (typically 25–75 kW m$^{-2}$) and the
mass-loss rate (MLR), surface temperature, and oxygen-consumption heat-release
rate (HRR) are recorded. Pyrolysis.jl, being a condensed-phase solver, predicts
MLR and $T_s$ directly and an HRR *proxy* via the effective heat of combustion
times the gas-generation rate; the gas-phase flame is not modeled. The shipped
`examples/06_cone_calorimeter.jl` (PMMA, 25 kW m$^{-2}$, top convection +
re-radiation, back-face re-radiation and convection with an impermeable back,
single-step Arrhenius, temperature-dependent properties) is the reference workflow; its surface-
absorption variant `06_cone_calorimeter_surface_absorption.jl` exercises the
opaque-material radiation path.

The validation QoIs are the peak MLR and its timing, the residual (char) mass
fraction, the time-to-ignition surrogate (a surface-temperature or
gas-generation threshold), and the $T_s(t)$ history. Two implementation points
materially affect quantitative agreement and are themselves verified:

- **Surface temperature in cell-centered FV.** Temperature lives at cell
  *centers*, but a $T^4$ radiation BC must be evaluated at the *surface*. Using
  the boundary-cell-center temperature directly produces an $\sim$20% energy
  imbalance for strong radiative BCs; `solve_boundary_temperature` instead Newton-
  iterates for the true $T_s$ satisfying $q_\mathrm{applied}(T_s)=
  q_\mathrm{conducted}(T_s,T_\mathrm{cell})$ and stores it for read-out
  (`get_surface_temperature`, §12.2). Validation against measured $T_s$ requires
  this solved value, not the cell-center temperature.
- **MLR gauge normalization.** For a fixed cross-section MLR is per unit area; for
  a shrinking variable cross-section (WithChi) the *gauge* MLR is normalized to
  the initial area, $\mathrm{MLR}_\mathrm{gauge}=\mathrm{MLR}_\mathrm{local}\cdot
  A(t)/A_0$ (§10, `src/Solver/solve.jl`), matching the cone gauge convention so
  the prediction is comparable to a measurement reported per *initial* sample
  area.

### 18.6.2 Douglas fir and the convergence case

`examples/07_douglas_fir_cone.jl` is the most fully-featured validation case:
multi-step wood decomposition kinetics, surface-absorption radiation (wood is
opaque), Darcy–Fick gas transport (species molar masses enter the pressure
closure; all gases share one diffusive transfer coefficient, §5.11), ALE
tracking of the shrinking sample, and the `ApproxSolve` reduced Jacobian for a
$\sim$8× speedup over the full structured Jacobian. The companion
`convergence_douglas_fir.jl` (§18.5.1) establishes the mesh resolution at which
the validation QoIs are converged, so that any residual disagreement with
experiment is attributable to physics or parameters rather than to
under-resolution. Validation here checks the MLR$(t)$ shape (the characteristic
two-peak signature of hemicellulose/cellulose and lignin decomposition), the
char yield, and the back-face temperature rise (a sensitive probe of the in-depth
conduction and the effective conductivity mixing rule, §5.5).

The FDS-style kinetic recovery of §18.2.2 is the practical link between TGA data
and the cone model: $(A_1,E_1)$ are recovered from a measured TGA peak temperature
and peak rate via the inflection relation, then carried into the cone run. This is
a *calibration* step (the kinetics are fit to TGA), so the cone comparison that
follows is a genuine *validation* only insofar as the cone data were not used in
the fit — a distinction §18.1.1 insists on.

### 18.6.3 Pressure-treated wood and moisture transport

`examples/08_pressure_treated_wood_cone.jl` (and `ptw_comparison.jl`) validate the
Darcy-driven transport path: a permeable charring solid develops an internal
pressure $P=\sum_g\xi_jR_gT/M_j/\varphi$ (§9.4), and the Darcy velocity
$v_g=-(\kappa_\mathrm{face}/\mu_\mathrm{face})\nabla P$ (§9.5) advects the
pyrolysis gases outward. The validation QoI is the MLR augmentation and the
internal-pressure history relative to a Fickian-only run; the Darcy–Fick demo
`08_darcy_fick_demo.jl` isolates the pressure-driven transport for a controlled
check. Moisture transport (`moisture_hos_comparison.jl`) validates the
state-dependent heat-of-sorption path: the evaporation "reaction" uses a
`StateDependentProperty` heat $h(T,\xi)$ that depends on the local moisture
concentration (§6.6), and the moisture front, the evaporation-plateau in $T(t)$
(the constant-temperature stall while water boils off), and the dried-mass profile
are the comparison targets.

### 18.6.4 Comparison to ThermaKin2Ds, Gpyro, and FDS

Code-to-code comparison (`examples/literature_comparison/`) is run alongside
experimental validation because it isolates *modeling-convention* differences
from *numerical* ones. The three reference codes differ from Pyrolysis.jl in
documented, specific ways; mapping their notation to ours:

- **ThermaKin2Ds.** Maps its stoichiometric coefficients $\theta_i^{\,j}$ to our
  $\nu_{i,j}$ and its component concentrations to our $\xi_j$. ThermaKin's energy
  equation (its Eq. 17) is Pyrolysis.jl's energy equation term for term: the
  storage term $\sum_j\xi_j c_{p,j}$ summed over *all* components including gas,
  the convective term, and a deformation term that the ALE mesh frame absorbs
  (§3.3.2). ThermaKin's convective energy term
  $c_g(J_g\cdot\partial T/\partial x)$ is identical in form to our $S_\mathrm{conv}$
  (the central-difference gradient matches ThermaKin's discretization, §7.4), so
  the convective coupling is a clean cross-check.
  ThermaKin supports *bimolecular* reactions ($\xi_\mathrm{COMP1}\,
  \xi_\mathrm{COMP2}$); Pyrolysis.jl enforces single-reactant kinetics
  ($N_\mathrm{Rin}=1$, §6.5), so multi-reactant ThermaKin schemes must be
  re-expressed as reaction sequences before comparison.

- **Gpyro.** Maps its pre-exponential $Z$ to our $A_i$, its permeability $K$ to
  our $\kappa$, its volume fractions $X_\alpha$ to our $v_j$, and its volumetric
  rate $\dot\omega'''$ to our $r$. Gpyro is fully conservative (control-volume
  conservation of $\bar\rho\Delta z$, $\bar\rho Y_i\Delta z$, gas mass, and
  energy) and uses volume-fraction property averaging with power-law
  temperature-dependence; Pyrolysis.jl's WEIGHTED/PARALLEL/SERIES/BRUGGEMAN
  mixing rules (§5.5) span and generalize Gpyro's averaging. Gpyro's Darcian mass
  flux $\dot m''=-(K/\nu)(\partial P/\partial z-g\rho_g)$ matches our Darcy form
  (§9.5) up to the buoyancy term $g\rho_g$, which Pyrolysis.jl omits (negligible
  for thin condensed-phase slabs). Gpyro selects swelling/shrinkage mechanisms by
  solid-fraction state; Pyrolysis.jl's lateral-shrinkage law $A(t)=A_0\cdot
  \mathrm{law}(\bar\chi)$ (§10.5) is the 1D analog. Gpyro's solver is a
  fully-implicit backward-Euler / Patankar linearization with TDMA; Pyrolysis.jl
  uses higher-order SDIRK with a structured Jacobian — so code-to-code agreement
  on a steady or slowly-varying case validates the *physics* independent of the
  time integrator.

- **FDS solid phase.** Maps its absorption $\kappa_s$ to our $\alpha$ and its
  pre-exponential $Z$ to our $A_i$. FDS solves 1D Fourier conduction with
  Arrhenius kinetics ($\rho^n O_2^n T^n$ rate forms) and a two-flux
  (Schuster–Schwarzschild) in-depth radiation model; Pyrolysis.jl's Beer–Lambert
  and P1 models bracket the FDS two-flux treatment. FDS carries no in-depth gas,
  so its energy equation is the condensed-phase limit of ours; the FDS
  TGA-reference kinetic recovery (peak temperature and
  rate → $A,E$) is the same recipe used in §18.2.2 and §18.6.2. FDS's surface
  energy balance $-k\,\partial T/\partial x=q_c''+q_r''$ is the multi-mechanism
  analog of our `solve_boundary_temperature` Newton balance (§12.2).

In all three comparisons, agreement *after* matching conventions (gas treatment,
radiation model, kinetic form) is the verification that the shared physics is
implemented consistently; *residual* disagreement after convention-matching is
attributable to the documented modeling differences, not to error.

---

## 18.7 Regime-validity checks

A model is only valid in the regime where its approximations hold. Pyrolysis.jl
exposes runtime diagnostics that tell the user when a run is leaving the
validated regime — these are validation aids, not verification checks, because
they assess the *applicability of the equations*, not the *correctness of the
solve*.

### 18.7.1 Optical-thickness and P1 validity

The optical thickness $\tau_\mathrm{total}=\int\alpha\,dz=\sum_i\alpha_i\Delta z_i$
selects the appropriate radiation model and is computed by `optical_thickness`
(Experimental). The validated regimes (§8): $\tau>10$ optically thick (P1
excellent, surface absorption acceptable); $1<\tau<10$ intermediate (P1 good,
Beer–Lambert recommended); $0.1<\tau<1$ optically thin (Beer–Lambert); $\tau<0.1$
transparent (use `NO_RADIATION`). `check_p1_validity` flags a run whose optical
thickness places it outside the P1-accurate band, and `radiation_energy_balance`
confirms the radiation source integrates to the absorbed-minus-emitted flux. The
absorbed-plus-transmitted identity of §18.2.3 is the underlying conservation
check; the validity query layers a *regime* judgement on top of it.

### 18.7.2 Pore-gas thermal equilibrium

The energy equation stores the pore gas at the matrix temperature (§7.3). That
one-temperature assumption holds when the gas equilibrates with the pore walls
within a cell, which is the case for pore-scale heat-transfer lengths far below
the cell size — every validated regime here (1 atm, pore sizes well below the
mesh spacing, gas residence times of milliseconds). A run with coarse open pores
or very fast, hot throughput, where the gas does not equilibrate with the
walls, is outside the model: there is no two-temperature gas option (§3.9), and
the prediction should be treated with caution.

### 18.7.3 Péclet number and upwind diffusion

The Darcy–Fick species transport uses first-order upwinding for the advective
part (§9.3). First-order upwind introduces numerical diffusion proportional to the
cell Péclet number $\mathrm{Pe}=|v_g|\Delta z/\lambda$; when $\mathrm{Pe}\gg1$
the scheme is robust but the artificial diffusion can smear the gas-composition
front. The convergence study (§18.5.1) exposes this as a first-order contribution
to the spatial error of pressure-driven cases. The practical regime check is to
confirm that the cell Péclet number is moderate over the resolved front — if it is
large, the user should refine the mesh near the front (a stretched mesh) until the
upwind diffusion is below the physical diffusion. The conduction analog is benign
(the conductive flux is central-differenced and unconditionally second-order), so
this regime concern is specific to the advective gas transport.

### 18.7.4 Depletion-threshold sensitivity

Two thresholds govern near-zero-concentration behavior, and a validated run should
be insensitive to both within their physically-negligible range. The *depletion
limiter* threshold (default $1.0$ kg m$^{-3}$, min $0.01$, §6.3) sets where the
tanh rate-suppression begins; the *surface-depletion* threshold (default 5% of
each cell's own initial thickness, `depletion_threshold=0.05`, §14.3,
`src/Solver/solve.jl`) sets when a thinned (ALE) or emptied (fixed mesh) cell is
parked. A sensitivity check varies
each threshold by an order of magnitude and confirms the QoIs (peak MLR, char
yield) move by less than the tolerance; a strong dependence means the threshold is
acting in a physically-active region rather than only suppressing
numerical stiffness near true depletion, and the value should be lowered. These
are AD-relevant: the limiter is tanh-smooth (so sensitivities are well-defined,
§17.6), but a park is a discrete rewrite of the state (whose size stays fixed)
and is *not* differentiable — a caveat for sensitivity studies (§17.7, §18.8.3).

---

## 18.8 Known discrepancies and re-validation needs

Honest V&V records the open items. The following are documented incompletenesses
that affect, or require re-checking of, validated results.

### 18.8.1 Surface gas-enthalpy convention and calibrated examples

The energy balance follows the $S_\mathrm{conv}$-only convention (ThermaKin
Eq. 17 / FDS): $S_\mathrm{conv}$ carries the escaping-gas enthalpy and the
surface balance has no gas-enthalpy term (§7.4.3), since a surface term would
count that enthalpy a second time. A balance that adds a surface term alongside
$S_\mathrm{conv}$ artificially depresses the surface temperature, so parameter
sets calibrated in a solver that does so read a *higher* $T_s$ and mass-loss rate
here (§1.8, item 7). The fitted kinetic and property parameters of the bundled
calibrated examples — the inverse-analysis cases and the cone runs — *need
re-validation* against data, together with the flame heat-flux feedback and blowing that a
recalibration should include (§1.8.1). This is a validation-state issue, not a
verification failure.

### 18.8.2 Composition-varying discrete energy residual does not close exactly

The energy ledger (§18.4.2) is missing the
$\int T\,\partial(\rho c_p)/\partial t\,dV\,dt$ heat-capacity-change term of its
$\rho c_p\,T$ proxy density (heat capacity changing at constant temperature as
composition reacts — §16.7) and an ALE volume-change term. As a result it cannot close to integrator tolerance when
composition varies in time; the residual carries a small composition-rate-
dependent offset. A ledger that closes exactly would require per-substep
accumulation of $\rho c_p(t_i)\,\Delta T_i\,V(t_i)$ per cell, which is not
implemented. This is a limitation of the *audit*, not of the solver: the solver's
energy update is correct; the bookkeeping formula used to *check* it is
incomplete. Until it is extended, a non-zero residual during active
reaction should be interpreted against this known offset, not as a bug.

### 18.8.3 Sensitivities across park events

Depletion handling is a sequence of discrete state rewrites (§14.3.4, §17.7).
The Jacobian and the verifier run through them (§18.3.5), but smooth
sensitivity analysis does not: forward sensitivities are valid only between
park events and adjoint sensitivity through an event is unsupported. Runs that
need smooth sensitivities must stay on horizons without park events.

### 18.8.4 Heuristic and unvalidated closures

Several constitutive closures are heuristic and carry their own validation debt:
the Lennard-Jones parameter estimates for Chapman–Enskog diffusivity
($\sigma\approx3.5+0.5\log_{10}(M/0.016)$,
$\varepsilon/k\approx150+100(M/0.044)$, §5.10) are explicitly flagged unvalidated;
the lateral-shrinkage laws (§10.5) are empirical rather than first-principles; and
the single global radiation model is a simplification. Validation cases that depend
sensitively on these closures should report that dependence.

### 18.8.5 Stale-flux and Marshak-coefficient regressions

Two silent-failure modes carry standing regression guards. The incident
radiation flux is re-sampled at *every* RHS evaluation with no constancy
detection (§8.3.4), so a time-varying flux can never be frozen at its initial
value; and the Marshak BC uses the exact
$\frac{\varepsilon}{2(2-\varepsilon)}$ coefficient rather than
$0.5\varepsilon$, with the external flux at its half-range weight
$2\varepsilon/(2-\varepsilon)$ rather than 1 (§18.2.5). Both would be silent
errors if regressed, so the equilibrium P1 sign test, the cold-absorber
absorbed-fraction test, and a modulated-flux comparison
($q_0(1+\sin 2\pi t)$) are kept as permanent checks.

---

## 18.9 Reproducibility

Verifiable results must be *reproducible*, which requires a controlled software
environment and an automated quality gate.

### 18.9.1 Quality tooling

Three static/dynamic quality tools run in isolated environments (so the
package's own Manifest is not pinned to their dependencies):

- **JuliaFormatter** — enforces a single canonical source format, so diffs are
  semantic and review is reliable.
- **Aqua.jl** — audits package hygiene: undefined exports, method ambiguities,
  stale dependencies, unbound type parameters, piracy. A clean Aqua run is a
  precondition for trusting that the public API surface (UG-16) is consistent.
- **JET.jl** — static type-stability and dispatch analysis. Because the hot paths
  rely on `@generated` unrolling over heterogeneous component tuples and on
  end-to-end `Tu`-typed (Float64 or `ForwardDiff.Dual`) caches for AD-cleanliness
  (§4.6.5, §7.7), a JET regression that introduces a type instability or a dynamic
  dispatch would both slow the solve and threaten the AD/sensitivity paths.

These run via dedicated project environments (`--project=scripts`,
`--project=test/quality`); the main test manifest is intentionally separate
because the quality tools' pins are fragile.

### 18.9.2 Test layout

The test suite mirrors the source module tree (`test/` mirrors `src/`), and
`runtests.jl` auto-walks the directories so a new test file is picked up without
manual registration. The full suite is run via `Pkg.test()`. The verification
machinery of this chapter has direct test coverage:

- Jacobian verification: `test/perf/test_verify_jacobian.jl` drives
  `verify_jacobian_at` over the canonical-configuration × backend matrix;
  per-operator tests (`test/jacobian/test_operator_*.jl`) check each operator's
  Jacobian row against AD.
- Conservation: the tracker tests exercise the mass and energy ledgers, and the
  verification case V-CPL-06 pins the gas-enthalpy coupling, and
  `test/boundary_conditions/test_boundary_temperature.jl` asserts that the surface
  balance carries no gas-enthalpy term (§18.4.3).
- The verification suite: `test/verification/runcases.jl` runs every case of
  `test/verification/cases/` in assert mode (§18.1.3).
- Radiation: `test/radiation/test_p1_sign.jl` (sign), `test_p1_validation.jl`
  (analytic limits), `test/jacobian/test_operator_radiation.jl` (Beer–Lambert
  upper-triangular structure).
- Allocation discipline: `test/perf/test_zero_allocation_rhs.jl` asserts the RHS
  assembly allocates zero bytes — a reproducibility-relevant performance
  invariant, since allocations in the hot path would couple results to GC timing.

### 18.9.3 Determinism and environment

Reproducible numbers require a pinned `Manifest.toml` (the integrator,
linear-solver, and AD package versions all affect the last digits) and a fixed
choice of integrator, tolerances, Jacobian backend, and linear solver. The
defaults — `KenCarp4` with `AutoFiniteDiff`, the pattern-chosen factorization
(`KLUFactorization` or `UMFPACKFactorization`), `abstol=1e-8`/`reltol=1e-6`,
`Structured(DirectSolve())` (`ApproxSolve` under ALE) — are the reference
configuration; deviations (e.g. `ApproxSolve` on a fixed mesh, `ILUGMRESFactorization`
for large ALE runs)
trade a controlled amount of accuracy for performance and should be re-checked
against the reference configuration with `verify_jacobian` armed when first
adopted on a new problem.

---

## 18.10 Limitations and open issues

The V&V coverage above has gaps that a user should weigh:

1. **MMS coverage is bounded by the solid-Lagrangian formulation.** The
   coupled MMS driver (§18.2, `test/verification/exact/mms.jl`) closes the
   design-order convergence gap for the fully coupled residual, but mesh
   motion in the MMS cases is either reaction-driven or prescribed via the
   verification-only `w_override` hook; under prescribed (non-material)
   motion the species sources are assembled in the material-following frame
   and the Jacobian falls back to finite differences (§11.9 item 5). General
   ALE (a species advection operator for mesh motion decoupled from the
   solid) remains out of scope.

2. **Convergence reference is self-consistent, not exact.** The convergence
   studies use the finest mesh as reference (Richardson self-convergence); they
   confirm *grid convergence* and the *rate*, but not the absolute error against
   a true solution. For QoIs with a known analytic value (the conduction
   eigenmode, the Beer–Lambert transmitted fraction) the absolute error can be
   bounded; for the full nonlinear cone case it cannot.

3. **Audit ledgers are incomplete for varying composition and ALE volume change**
   (§18.8.2). The discrete energy residual does not close exactly during active
   reaction; the missing heat-capacity-change and volume-change terms must be
   accounted for when interpreting it.

4. **Depletion is non-differentiable** (§18.8.3): the structured Jacobian and its
   verifier run through park events (§18.3.5), but sensitivities are valid only
   between them, so the most strongly-shrinking cases (deep char, ablation) have
   the narrowest sensitivity horizons.

5. **P1 radiation is Experimental and separately verified** (§18.2.5, §8.12): it
   has its own residual outside the production `solve` path, so its V&V does not
   ride the unified test matrix and must be exercised through the Experimental
   tests.

6. **Validation is parameter-entangled.** Several cone/inverse cases were
   calibrated against the data they are compared to, and those calibrations need
   re-fitting (§18.8.1). Re-fitting against held-out data
   — keeping calibration and validation data disjoint — is the outstanding
   validation task, and is the discipline §18.1.1 requires for a comparison to
   count as validation rather than calibration.

7. **Heuristic closures carry validation debt** (§18.8.4): the LJ diffusivity
   estimates and the empirical shrinkage laws are unvalidated or simplified, and cases that lean on them inherit that
   uncertainty.
