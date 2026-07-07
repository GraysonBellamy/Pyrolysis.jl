# Verification Suite Plan — 2026-07-06

Design document for a comprehensive verification (exact-solution) and benchmark
(code-to-code) suite for Pyrolysis.jl, grounded in what the established
condensed-phase pyrolysis codes were verified against:

- **FDS Verification Guide** (NIST SP 1018-2), Heat Conduction + Pyrolysis chapters;
  case sources in `firemodels/fds` `Verification/Pyrolysis/`.
- **Gpyro** — Lautenberger dissertation §3.4 (*A Generalized Pyrolysis Model for
  Combustible Solids*, UC Berkeley 2007, [eScholarship](https://escholarship.org/uc/item/7wz5m7dg))
  and the Gpyro3D paper (Fire Safety Science 11:193–207, 2014).
- **ThermaKin** — Stoliarov & Lyon, Fire Safety Science 9:1141–1152 (2008)
  (peer-reviewed companion to FAA report DOT/FAA/AR-TN08/17), "Model Verification" section.
- **MaCFP Condensed Phase WG** — `MaCFP/matl-db` target simulations (MaCFP-2 PMMA,
  MaCFP-4 pine wood, 2026) — these are code-to-code, not exact-solution.

Terminology follows the standard split: **verification** = solving the equations
right (exact/manufactured solutions, convergence order, conservation);
**validation** = solving the right equations (experiments). This plan is about
verification plus code-to-code benchmarks; experimental validation (cone
examples, MaCFP experimental data) is the guide's Part 2 and comes later.

## Status of existing coverage

`test/verification/test_analytical.jl` already covers: semi-infinite constant-flux
conduction (erf) with 2nd-order spatial convergence, Beer–Lambert absorption
profile + energy split, isothermal first-order decay, constant-heating-rate TGA
vs high-precision temperature integral + DTG peak + Kissinger E-recovery, a
decaying-sine exact solution, inert/reacting mass & energy ledgers, ALE mass
conservation + GCL, and steady Darcy velocity. Targeted tests exist for
Kirchhoff properties, P1 limits (Experimental), analytical Jacobians, blowing
correction, and non-uniform gradients.

The cases below are **additions**, chosen to (a) match what the reference codes
verified, (b) exercise couplings none of our current tests reach (moving
boundary with exact solution, latent-heat front, pressure-driven outflow,
reaction–conduction feedback), and (c) close the documented MMS gap
(`docs/src/guide/technical_reference/18_verification-validation.md` §18.2).

Case IDs: `V-HC` heat conduction, `V-RX` kinetics, `V-RAD` radiation, `V-GT` gas
transport, `V-CPL` coupled/moving-boundary, `V-MMS` manufactured solutions,
`B-` code-to-code benchmarks. Every case must report an error norm **and** an
observed convergence order where a spatial/temporal discretization is exercised
(matching FDS's practice for `convective_cooling` and `ht3d_sphere`).

---

## Tier 1 — Single-physics, closed-form exact solutions

### V-HC-01 — Semi-infinite solid, convective surface BC (Biot-coupled erfc)

The one verification case shared by **all three** reference codes
(FDS `ht3d_slab`, Gpyro §3.4.1, ThermaKin Eq. 25). Carslaw & Jaeger:

```
(T(x,t) − T₀)/(T∞ − T₀) = erfc(η) − exp(hx/k + h²αt/k²)·erfc(η + h√(αt)/k),   η = x/(2√(αt))
```

Unlike our existing constant-flux erf case, this exercises `ConvectiveBC` and
the **surface Newton solve** (the surface temperature is an unknown coupled to
the interior), which is where BC bugs live. Gpyro's variant adds a constant
absorbed radiant flux (25 kW/m²) plus convective loss (h = 20 W/m²K) — the exact
solution generalizes (superposition, constant properties); use Gpyro's
parameters (k = 0.2, ρ = 1000, c = 1400) so results are directly comparable to
dissertation Fig. 3.x profiles at t = 30/60/90/180 s.

### V-HC-02 — Finite slab, convective BC, Biot sweep (FDS `heat_conduction_a–d`)

0.1 m slab, convective front, adiabatic back, Bi = hL/k ∈ {0.1, 1, 10, 100}.
Exact eigenfunction series, eigencondition `λₙ tan λₙ = Bi`:

```
θ/θᵢ = Σₙ Cₙ exp(−λₙ² Fo) cos(λₙ x/L),   Cₙ = 4 sin λₙ / (2λₙ + sin 2λₙ)
```

The Biot sweep is the point: Bi = 0.1 is lumped-like (BC-dominated), Bi = 100 is
surface-Dirichlet-like (interior-dominated). One case definition, four regimes,
front- and back-face temperature histories as targets.

### V-HC-03 — Finite slab, constant flux, adiabatic back (exact series)

Complements the existing semi-infinite erf test at **long times** (thermal wave
reaches the back face — the regime cone-calorimeter simulations actually live in):

```
T(x̃,t) = T₀ + (q″L/k)·[Fo + (3x̃² − 1)/6 − (2/π²) Σₙ ((−1)ⁿ/n²) e^(−n²π²Fo) cos(nπx̃)],  x̃ from back face
```

Target both the front/back ΔT plateau `q″L/2k` and the late-time linear-in-time
mean rise (a global energy check with an exact constant).

### V-HC-04 — Two-layer slab: steady piecewise-linear profile + interface flux

FDS `insulated_steel_plate` / `ht_pile` analog. Since laminates enter via
composition step ICs (`step_profile`), this doubles as a verification of the
**distance-weighted harmonic-mean face conductivity** at a property jump:
steady profile is exactly piecewise linear, interface temperature from the
resistance network, and the numerical steady flux must be exact (harmonic mean
is nodally exact for piecewise-constant k). Run with k-ratio 10:1 and 100:1,
including a mesh where the interface does **not** coincide with a face
(quantifies the smearing error users will actually see).

### V-HC-05 — Temperature-dependent k(T): steady Kirchhoff-exact profile

Steady conduction across a slab with k(T) = k₀(1 + βT), fixed T at both faces.
Kirchhoff transform gives the exact implicit profile
`k₀(T − T₁) + k₀β(T² − T₁²)/2 = C·x` and exact flux
`q″ = (k₀/L)[(T₁ − T₂) + β(T₁² − T₂²)/2]`. Exercises the T-dependent property
path and face interpolation under strong nonlinearity (β chosen so k varies 3×
across the slab — comparable to char vs virgin wood). Transient T-dependent
{k,c} has no exact solution (FDS resorts to HEATING/ABAQUS code-to-code — see B-02).

### V-RX-01 — Two-step series chain A→B→C (the FDS/Gpyro shared benchmark)

Both FDS (`two_step_solid_reaction`) and Gpyro (§3.4.4) use the identical
problem: constant rates k₁ = 0.389 s⁻¹, k₂ = 0.262 s⁻¹, exact

```
ξ_A = e^(−k₁t),  ξ_B = k₁/(k₂−k₁)·(e^(−k₁t) − e^(−k₂t)),  ξ_C = 1 − ξ_A − ξ_B
```

Exercises multi-step `ReactionSet` coupling (intermediate B is produced and
consumed simultaneously — our existing isothermal test is single-step). Use the
published rate constants so the numbers are cross-checkable against both codes.
Run isothermal (Dirichlet T, reactions gated fully open). The tanh depletion
limiter must be disabled (`DepletionLimiter(enabled=false)`) or its threshold
subtracted from the tolerance budget — document which.

### V-RX-02 — nth-order and contracting-geometry isothermal decay

FDS `part_baking_soda` verifies against Khawam & Flanagan (2006) exact isothermal
integrations. For n ≠ 1: `ξ(t) = [ξ₀^(1−n) − (1−n)k't]^(1/(1−n))` (finite
burnout time for n < 1). Sub-cases: n = 0 (linear, exact burnout time —
stresses the depletion limiter deliberately), n = 2 (ThermaKin's second-order
C+C reaction, representable here as single-reactant n = 2), n = 2/3
(contracting-sphere R3 model, FDS's NaHCO₃ case: m^(1/3) linear in t; use their
A = 3.4×10¹¹ 1/s, E = 103 kJ/mol at T = 420/450/500 K).

### V-RX-03 — Competing parallel reactions A→νB·B, A→νC·C

Exact: `ξ_A = ξ₀e^(−(k₁+k₂)t)`; yields split as k₁:k₂ scaled by ν. Verifies
stoichiometric bookkeeping and heat-release superposition when both reactions
carry different h (exact energy ledger target, analogous to FDS `matl_e_cons`).

### V-RX-04 — Thermally-thin radiative cooling (exact implicit solution)

A thin slab (Bi ≪ 1) cooling purely by radiation: `dT/dt = −K(T⁴ − T∞⁴)` with
K = εσA/(mc) has the closed form t = [G(Tᵢ) − G(T)]/K where the antiderivative
(partial fractions, valid on the cooling branch T > T∞) is

```
G(T) = (1/(4T∞³))·ln((T − T∞)/(T + T∞)) − (1/(2T∞³))·arctan(T/T∞)
```

(Re-derive and unit-test G against numerical quadrature before use — standard
but sign-error-prone.)

This is the only closed-form transient involving the T⁴ nonlinearity — it
verifies `RadiativeBC` + surface Newton against something other than steady
algebra. Companion steady case: radiative–convective surface equilibrium
`ε q″ = εσ(Tₛ⁴ − T∞⁴) + h(Tₛ − T∞)` with adiabatic back → exact algebraic Tₛ
(polynomial root), directly checking the surface Newton residual.

### V-RAD-01 — In-depth absorption with conduction: eigenfunction exact solution

Lautenberger dissertation §3.4.2 + Appendix A (Green's function derivation):
finite slab δ = 1 cm, in-depth source `q″κe^(−κz)` with κ = 1000 m⁻¹,
convective loss at the front, insulated back — exact eigenfunction series with
eigencondition `cot(λδ) = kλ/h`. The physics marker is the **sub-surface
temperature peak**, which our existing Beer–Lambert test (profile + energy
split, no conduction coupling) never checks. This is the definitive
`BEER_LAMBERT` + conduction + `ConvectiveBC` verification. Use Lautenberger's
parameters for cross-code comparability.

### V-GT-01 — Instantaneous plane source: Gaussian diffusion (ThermaKin case)

Gas species released from a plane at x₀ inside the solid, pure Fickian
transport: `[C](x,t) = M_C/(2√(πdt)) · exp(−(x−x₀)²/(4dt))`. Verifies the
`Fickian` gas-transport operator + `ImpermeableBC` (mass conserved to
tolerance; boundary influence quantifiable via image sources). ThermaKin's
verification set includes exactly this (Crank plane-source solution).

### V-GT-02 — Steady compressible Darcy outflow with uniform generation (P² profile)

Isothermal ideal gas, uniform volumetric generation ġ (zero-order gas-producing
reaction at clamped T), sealed back at z = 0, `DirichletPressureBC(P∞)` at
z = L. Steady mass balance `d/dz[(κM/μRT)·P dP/dz] = −ġ` integrates exactly:

```
P²(z) = P∞² + (μRT ġ/(κM))·(L² − z²),    exit mass flux = ġL  (exact)
```

This is the canonical porous-flow exact solution (gas-reservoir form) and is a
real test — the equation is nonlinear in P, exercises the algebraic pressure
closure, harmonic-mean mobility, Darcy upwinding, and the outflow mass BC
simultaneously. Our current Darcy test only checks the velocity–gradient
relation pointwise. Parameter choice should give a back-face overpressure of
order 10–100 kPa (wood-pyrolysis-relevant magnitudes per Ch. 9 of the technical
reference).

---

## Tier 2 — Coupled physics with exact solutions (the non-trivial core)

These are the cases that make the suite "not too trivial": each couples two or
more operators and has a sharp, quantitative exact target.

### V-CPL-01 — Neumann two-phase Stefan problem (latent-heat front)

Gpyro's flagship verification (§3.4.3) — the only reference-code case with an
exact **moving front**. Two-phase solidification into a semi-infinite melt
(solid 0 < x < s(t) at wall temperature T_w < T_m, liquid x > s initially at
T_l0 > T_m); front position `s(t) = 2λ√(α_s t)` with λ from the classical
Neumann transcendental equation (Carslaw & Jaeger §11.2):

```
e^(−λ²)/erf(λ) − [k_l √α_s (T_l0 − T_m)]/[k_s √α_l (T_m − T_w)] · e^(−λ²α_s/α_l)/erfc(λ√(α_s/α_l))
   = λ √π · ΔH_m / (c_s (T_m − T_w))
```

Cross-check the transcribed form against dissertation Eq. 3.138d before use;
with Gpyro's setup (polypropylene-like: T_m = 433 K, ΔH_m = 50 kJ/kg) the root
must reproduce λ = 0.2317, which serves as the unit test of the reference
implementation itself.

Implementation in our framework: latent heat via a narrow `StateDependentProperty`
c_p peak (apparent-heat-capacity method — exactly Gpyro's Gaussian-peak
treatment) or via a reaction (solid_A → solid_B, γ chosen for zero volume
change, h = ΔH_m, tight tanh temperature gate). Both routes should be run; the
peak-width/gate-width convergence study (exact solution recovered as width → 0
with mesh refinement) is itself a valuable result and belongs in the guide.
Targets: front trajectory (√t slope), interface temperature plateau, and both
similarity profiles.

### V-CPL-02 — Landau ablation problem (exact steadily-regressing surface) ⭐

**The flagship ALE case, and a differentiator** — none of FDS/Gpyro/ThermaKin
verify against it (confirmed: ablation appears only in Gpyro's literature
review). Landau (1950): semi-infinite solid, constant net surface flux q″,
material removed at fixed ablation temperature T_p with heat of ablation ΔH_p.
The exact steady traveling wave:

```
v_ss = q″ / (ρ[c(T_p − T₀) + ΔH_p])           (regression rate)
T(ζ) = T₀ + (T_p − T₀)·exp(−v_ss ζ/α)          (ζ = depth below moving surface)
```

with thermal-wave thickness α/v_ss and exact steady MLR = ρ·v_ss. This
exercises, in one problem: ALE mesh motion driven by volume change, surface
regression, kinetics-controlled front (a high-activation-energy solid→gas
reaction with γ_gas = 0 approximates the fixed-T_p ablation limit), the surface
energy balance under mass loss, and `S_conv` advection of the escaping gas.

Verification protocol: (1) fix T_p by choosing (A, E) pairs along a constant-
"ablation temperature" locus and confirm v_ss is insensitive to the pair
(kinetics-limit convergence to the ablation limit as E ↑); (2) check MLR
plateau against v_ss to <1 %; (3) check the exponential in-depth profile in
surface-attached coordinates; (4) check total-thickness history slope. Also run
with the blowing correction OFF (the exact solution assumes prescribed net
flux) — a second variant with `ConvectiveBC` + blowing has no exact solution
and is a regime check only.

### V-CPL-03 — Shrink/swell exact endpoints (FDS `shrink_swell` analog)

Reaction A→B with prescribed intrinsic densities; FDS verifies exact final
thickness/density for six ρ_product ≷ ρ_reactant permutations. For our ALE +
swelling-factor γ model the exact final thickness is **not** the FDS formula —
it must be derived by integrating *our* volumetric-strain model
(θ = Σⱼ γⱼ ξ̇ⱼ/ρⱼ, technical reference Ch. 10/11) over complete conversion for
each permutation; the derivation belongs in the case docstring and its algebra
should be spot-checked against a brute-force fine-mesh run before the value is
pinned. Mass must be exactly conserved in all permutations. Extends the existing ALE
mass-conservation test to **exact geometric endpoints**, including the swelling
branch and the χ-driven lateral-area law (WithChi): final A/A₀ exact from the
prescribed law at χ = 1. (Note: there is no "γ > 1 intumescent branch" — γ is
validated to [0, 1] because the same weight feeds the density/conductivity
mixing rules; swelling is exercised through low-bulk-density products,
ρ_B < ρ_A, as V-CPL-03 P1 does. See tech-ref §4.1.5.)

### V-CPL-04 — Frank-Kamenetskii critical slab (reaction–conduction bifurcation)

Zero-order exothermic reaction in a slab with Dirichlet walls. In the FK limit
(E/RT₀ → ∞), criticality at δ_c = 0.878 with
δ = (E/RT₀²)·(Q A ρ L² / k)·e^(−E/RT₀) (L = half-width), and the exact
subcritical steady profile (x̃ = x/L ∈ [−1,1], θ = (E/RT₀²)(T − T₀)):

```
θ(x̃) = θ_m − 2 ln cosh( x̃ · √(δ e^(θ_m)/2) ),   θ_m fixed by θ(±1) = 0
```

(two θ_m roots below δ_c — take the lower, stable branch; at δ_c the roots
merge at θ_m ≈ 1.19). Protocol: run at E/RT₀ ≈ 30–50 where FK is accurate to
<2 %, bracket δ_c by bisection on q (or L), and compare the subcritical steady
center temperature to the exact FK profile; additionally pin the reference with
an independent high-precision BVP solve of the full-Arrhenius steady equation
(BoundaryValueDiffEq in the test, tolerances 1e-10). This is the only case in
the suite where reaction heat release feeds back on the temperature field with
a sharp, exactly-known stability boundary — a genuinely hard test of the
coupled Jacobian and of stiff-integrator robustness near a turning point. Note:
the smooth tanh temperature gates must be fully open across the relevant T
range so the kinetics are pure Arrhenius.

### V-CPL-05 — TGA-to-gasification consistency ladder (0D→1D exactness chain)

Not a new exact solution but an exactness *chain* used by MaCFP: the same
kinetics run (a) 0D isothermal (exact, V-RX cases), (b) 0D constant heating
rate (existing high-precision temperature-integral reference), (c) 1D thin-limit
gasification (as L → 0 and q″ balance → prescribed β, the 1D solve must
converge to the 0D TGA curve — an internal consistency limit with a quantified
convergence rate). Catches BC/volume-coupling errors that pure 0D tests can't.

---

## Tier 3 — Method of Manufactured Solutions (infrastructure)

### V-MMS-01 — General source-injection MMS driver

The documented gap (tech ref §18.2): we currently verify only against exact
solutions that happen to have zero residual source. Build a test-only wrapper
around the pure-residual layer (`ProblemDef`/`residual!`) that injects
manufactured volumetric sources S_T(z,t), S_ξⱼ(z,t) computed symbolically (or
via ForwardDiff on the manufactured fields) so that a chosen smooth
T*(z,t), ξ*ⱼ(z,t) is the exact solution of the **full coupled system** —
T-dependent properties, multi-step Arrhenius, Darcy–Fick transport all active.
This is the only way to get formal order-of-accuracy verification of the
coupled operators (FDS uses MMS for exactly this reason in the fluid phase).

Deliverables: (a) driver utility in `test/verification/mms/`; (b) MMS case with
all couplings on, log–log convergence with observed order ≥ 1.9 spatial /
matching the integrator's temporal order; (c) an **ALE-MMS** variant with a
prescribed smooth mesh motion z*(t) verifying GCL + moving-mesh transport terms
at design order (the hardest and most valuable single test in the suite).

---

## Tier 4 — Code-to-code benchmarks (community-anchored)

### B-01 — MaCFP-2 PMMA target simulations

The exact prescribed scenarios every major code ran (MaCFP 2020 guidelines):
0D TGA at 10 and 100 K/min (300–800 K, spatially uniform T); 1D gasification of
a 6 mm slab at 25 and 65 kW/m², anaerobic, perfectly insulated back; outputs
MLR and back-surface temperature. Adopt one published parameter set (e.g. the
NIST set from `matl-db/Verification/Model_predictions/`) verbatim and overlay
our curves on the submitted FDS/Gpyro/ThermaKin results. Watch the
heat-of-reaction **sign convention** (ours: h > 0 endothermic — opposite
ThermaKin/Gpyro) when transcribing parameters.

### B-02 — FDS `heat_conduction_kc` analog (transient T-dependent k, c)

No exact solution exists (FDS compares to HEATING 7.3 / ABAQUS). Reproduce the
FDS case definition and compare to the published FDS/HEATING curves; cheap
insurance on the transient nonlinear-property path that V-HC-05 only checks at
steady state.

### B-03 — MaCFP-4 pine wood (2026, timely)

0D TGA at 10/100 K/min; 1D gasification at 6 and 18 mm × 10/30/60 kW/m²,
anaerobic. Directly in our wheelhouse (multi-step wood kinetics + Darcy–Fick +
ALE, cf. `examples/07_douglas_fir_cone.jl`). **Consider submitting Pyrolysis.jl
as a participant** — the workshop is active in 2026 and participation would put
the code on the community map alongside UC Boulder / FSRI-UMD / CORIA / NIST
submissions. This is a strategic item, not just a test.

---

## Explicitly out of scope / known solver limits the suite cannot cover

- **Oxidative kinetics** (no O₂ dependence in `Reaction`): MaCFP-5 oxidation
  cases and FDS oxidative TGA cases are unrepresentable until the kinetics
  model grows an oxidizer term. Note in the guide as a limitation.
- **Catalytic reactions** (rate depends on a species not consumed — ThermaKin
  verification case 4): not representable with single-reactant stoichiometry.
- **Bimolecular cross-reactions** (A+B): single-reactant enforced.
- **P1 radiation**: verifiable only through `Pyrolysis.Experimental`
  (optically thick/thin limit tests already exist there); the E₃ plane-layer
  emission exact solution (FDS `wall_internal_radiation`, Zel'dovich & Raizer)
  becomes relevant if/when P1 or two-flux graduates to core.
- **First-class laminates with contact resistance**: V-HC-04 covers property
  jumps via composition; a true contact-resistance model has no home yet.

## Cross-cutting requirements (every case)

1. **Convergence order reported**, not just an error at one resolution
   (spatial via mesh sweep ×4; temporal via tol/dt sweep where relevant).
2. **Conservation ledgers on** (mass, energy, species) with case-specific
   budgets — reuse the existing ledger machinery.
3. **Both Jacobian strategies** where feasible: `DirectSolve` (exact) and
   `ApproxSolve` for ALE cases (documents the approximation's effect against a
   known-exact answer — currently only self-referenced).
4. **Reference solutions computed in-test at high precision** (series summed to
   machine tolerance, transcendental roots via bisection to 1e-14), never
   hard-coded floats without provenance; store derivations in the case file's
   docstring with the citation.
5. **Fixed random-free setup**: cases are deterministic; tolerances pinned with
   headroom (~3× observed error at the CI mesh) and a comment stating the
   observed value at pin time.
6. **Diagnostics are first-class, not just asserts**: every case must define
   its quantities of interest (QoIs), the exact-solution callables, and the
   comparison curves it wants plotted — so the same case definition can drive
   both a fast pass/fail test and a publishable report (see dual-mode
   architecture below). A case that can only say "pass" is incomplete.

## Dual-mode architecture: CI asserts vs. publishable reports

The existing test suite already takes 10+ minutes, and publication-quality
diagnostics (fine meshes, multi-resolution sweeps, figure rendering) are far
too expensive — and pull in plotting dependencies — to live in `Pkg.test()`.
So every case runs in one of two modes from a **single case definition**:

- **Assert mode** (in `Pkg.test()`, budgeted): coarse-but-sufficient mesh,
  endpoint/QoI checks against pinned tolerances, a 2–3 point convergence-order
  spot check where cheap. No plotting, no artifact I/O. Target: the whole
  verification group adds ≤ 2–3 min to the suite; per-case budgets recorded in
  the case metadata and enforced (a case that blows its budget moves its
  expensive parts to report mode, it doesn't slow CI).
- **Report mode** (manual / scheduled, isolated env — same pattern as
  `scripts/profile/` and `test/quality/`): full mesh-refinement ladder (×4–5
  resolutions), tight integrator tolerances, and generation of the standard
  artifact set below. This is what populates the guide; CI never runs it.

Standard artifact set per case (FDS-Verification-Guide style — each case in
their guide is exactly "problem statement → figure(s) overlaying expected vs
computed → error statement"):

1. **Overlay plot(s)**: solver solution vs exact/reference solution for each
   QoI — e.g. T(x) profiles at several times with analytical curves as lines
   and solver values as markers; T(t) or MLR(t) histories; front position s(t)
   vs 2λ√(αt); P(z) vs the exact P² profile. Reference curves drawn from the
   in-test high-precision evaluation, never digitized/hard-coded.
2. **Error diagnostics plot**: pointwise error profile e(x) or error-vs-time
   e(t) at the finest resolution — shows *where* the discretization struggles
   (boundary cells, fronts, interfaces), which a norm alone hides.
3. **Convergence plot**: log–log error norm (L₂ and L∞) vs Δz (and Δt/tol
   where relevant) with fitted observed order and a reference-slope triangle.
4. **Summary table** (machine-written): resolution ladder, error norms,
   observed order, expected order, tolerance, pass/fail, wall time — plus
   solver configuration (integrator, tolerances, Jacobian strategy) for
   reproducibility.
5. All numeric outputs also written as CSV next to the figures so the guide
   tables regenerate from data, and so future runs can be diffed against a
   committed baseline (regression detection beyond pass/fail).
6. **Cross-code overlays** (2026-07-07): where a case has a published analog
   in another solver, that solver's *numerical* result is committed under
   `test/verification/reference/<code>/` (provenance READMEs per directory)
   and drawn as an extra `CurveSeries` on the case's overlay/error figures.
   Sources verified 2026-07-07: FDS 6.11.0 run locally on its shipped
   verification inputs (FDS never commits its own outputs — firebot
   regenerates them nightly), Gpyro 0.8200 built from source with decks
   reconstructed from Lautenberger's dissertation §3.4, and the MaCFP-2 PMMA
   community submissions (downloadable CSVs; feeds the future B-01).
   ThermaKin ships as a Windows-only binary (`docs/refs/ThermaKin/`) — its
   FSS9 verification cases (incl. the V-GT-01 Gaussian) are fully
   parameterized in the paper, so either run the `.exe` on Windows/wine or
   request a rerun from the Stoliarov group (stolia@umd.edu); slot the CSVs
   into `reference/thermakin/`. Since every correct code lies on the exact
   curve at plot scale, exact-solution overlays compare *error* (their series
   gets the same `exact` callable, so the renderer derives their pointwise
   error panel next to ours); solution-curve overlays are reserved for
   genuinely-disagreeing comparisons (MaCFP, B-02 HEATING/ABAQUS).
   Presentation (`CurveSeries.crosscode = true` → render.jl): every affected
   figure is emitted twice — the base variant stays ours-vs-exact only (our
   error legible on its own scale), and a `_vs_codes` variant adds the other
   codes as a muted context layer (gray shades + per-code marker shapes,
   thinned, always drawn *under* our colored series). Cases with reference
   runs also get a `solve_time` figure: log-x lollipop of our wall time at
   every ladder rung vs the other codes' runs (`reference/timings.csv`,
   same-machine 2026-07-07, resolution recorded per row). A discarded warm-up
   solve precedes the timed ladder so Julia JIT compile is not billed to the
   coarsest rung. Where a reference run's resolution falls outside our
   ladder, the case's `n_cells_ladder` carries a matched rung (V-HC-01:
   133/500, V-HC-02: 101, V-HC-03: 200, V-CPL-01: 800) so error and solve
   time compare like-for-like; the committed decks also allow the converse
   (rerun their code at our resolution) if ever needed.

## Proposed layout & guide scaffolding (Phase 2, after cases land)

```
test/verification/
  cases/                 # one file per case: V-HC-01.jl, ... (definition + assert mode)
  exact/                 # reusable exact-solution library (series, roots, MMS driver)
  runcases.jl            # assert-mode harness used by Pkg.test()
scripts/verification/    # report mode — ISOLATED env (own Project.toml with
                         # CairoMakie/CSV; keeps plotting deps out of the test
                         # manifest, same rationale as scripts/profile/)
  run_reports.jl         # runs cases/ in report mode → figures, CSVs, tables
docs/src/guide/verification/
  index.md               # suite overview, methodology, tolerance policy
  V-HC-01.md ...         # one generated chapter per case: statement, exact
                         # solution + derivation/citation, figures, error/order
                         # table, tolerance provenance, solver config
  assets/                # generated figures + CSVs (regenerated, never hand-edited)
```

Each case file carries a `VerificationCase` metadata block (id, title, physics
tags, references, QoIs, exact-solution callables, CI resolution + budget,
report resolution ladder, tolerances, expected orders) that all three
consumers read — the assert harness, the report runner, and the Documenter
chapter generator — the same single-source pattern as `scripts/split_guides.jl`.
Chapter prose (problem statement, solution derivation) lives in the case file's
docstring/metadata so the guide never drifts from what the test actually runs.

Workflow integration: report mode runs as a manual `workflow_dispatch` (and
optionally on release tags) that regenerates `docs/src/guide/verification/`
and can be reviewed as a PR diff of figures/CSVs; the docs deploy
(`docs.yml` → Cloudflare Pages) then publishes it. `Pkg.test()` stays
plot-free and within budget.

## Implementation handoff notes (read before writing code)

Orientation for whoever picks this up fresh:

**Read first.** `docs/src/guide/technical_reference/` chapters 03 (governing
equations), 06 (kinetics), 10–11 (volume change/ALE), 12 (BCs), 18 (existing
V&V methodology and its stated MMS gap). Then skim
`test/verification/test_analytical.jl` end-to-end — it is the house style for
analytical tests (fixture helpers, tolerance pinning, convergence-rate
assertions) and the thing this suite extends. `test/runtests.jl` auto-walks
test directories, so new files under `test/verification/` are picked up
without registration.

**Case setup skeleton.** A case is: `SolidComponent`/`GaseousComponent` →
`Reaction(...)` → `Material` → `create_uniform_mesh`/`create_stretched_mesh` →
`BoundaryConditionSet(thermal=(top=…, bottom=…), mass=(…))` →
`PyrolysisProblem` → `solve(problem; …)`. Extract results via `sol.T`, `sol.ξ`,
`get_cell_positions(sol)`, `get_surface_temperature`, and
`sol.extras.{total_mass, gas_generation_rate, thickness_history}` (MLR =
`gas_generation_rate`). Working examples: the module docstring in
`src/Pyrolysis.jl` and every testset in `test_analytical.jl`.

**Solver conventions that will silently corrupt a case if missed:**

- **Heat of reaction sign: h > 0 is ENDOTHERMIC** — opposite to
  ThermaKin/Gpyro/FDS conventions. Flip signs when transcribing any literature
  parameter set (bites hardest in B-01/B-03 and V-CPL-01/02).
- **Depletion limiter**: `Reaction` applies a tanh depletion limiter by
  default (threshold ~1 kg/m³). Exact-kinetics cases (all V-RX, V-CPL-02/04)
  must pass `DepletionLimiter(enabled=false)` or account for it in the
  tolerance budget — verify the current default in `src/Materials/reactions.jl`
  before assuming.
- **Temperature gates are smooth tanh**, not Heaviside: for pure-Arrhenius
  cases set `T_min=0, T_max=Inf` (gates fully open); for gate-dependent cases
  (Landau ablation-limit, Stefan-via-reaction) the gate width is a convergence
  parameter and must appear in the report-mode sweep.
- **`use_transpiration_bc=true` is rejected at construction** (would
  double-count S_conv, which is always on). Blowing correction is opt-in and
  should stay OFF for every exact-solution case (the solutions assume
  prescribed net flux / plain h).
- **ALE**: `solve(...; use_ale=true, min_thickness=…)`; γ (swelling factor)
  defaults: solid = 1, gas = 0. Mesh is pinned at the bottom; the surface
  moves — build regression cases (Landau) with the heated face on `top`.
- **Reference numerics**: default integrator KenCarp4 with
  `abstol=1e-8, reltol=1e-6`; report-mode runs should tighten these until
  temporal error is subdominant to spatial (verify, don't assume).

**Formula provenance rule.** The closed forms in this document were
transcribed from the cited sources and standard texts but have NOT been
independently re-derived — treat each as a claim. Before pinning any
tolerance: re-derive (or verify against the citation), and unit-test the
reference implementation itself against an independent check (numerical
quadrature for integrals, high-precision ODE/BVP solve for profiles,
`λ = 0.2317` for Stefan, published FDS/Gpyro curve values where given). The
exact-solution library in `test/verification/exact/` should carry these
self-tests permanently.

**Suggested first PR** (small, proves the architecture end-to-end): the
`exact/` library skeleton + `VerificationCase` metadata type + assert-mode
harness + **V-HC-01 only**, wired into `Pkg.test()` within budget. Second PR:
report-mode runner in `scripts/verification/` (isolated Project.toml, CairoMakie)
rendering V-HC-01's artifact set + generated guide chapter. Everything after
that is adding cases to a working pipeline, roughly in the priority order
below. Keep quality tooling in mind: JuliaFormatter/Aqua/JET run in isolated
envs (`--project=scripts`, `--project=test/quality`).

## Suggested implementation order

| Priority | Cases | Rationale |
|---|---|---|
| P0 | V-HC-01/02, V-RX-01, V-RAD-01, V-GT-02, V-CPL-02 (Landau) | Highest new-coverage per effort; Landau is the flagship ALE proof |
| P1 | V-CPL-01 (Stefan), V-MMS-01 (+ALE-MMS), V-CPL-03, V-RX-02/03, V-RX-04, V-HC-03/04/05, V-GT-01 | Completes the exact-solution core + closes the MMS gap |
| P2 | V-CPL-04 (FK), V-CPL-05, B-01, B-02 | Hard coupled physics + community anchoring |
| P3 | B-03 (MaCFP-4 participation decision), guide scaffolding | Strategic / documentation phase |

## Primary references

- H.S. Carslaw & J.C. Jaeger, *Conduction of Heat in Solids*, 2nd ed. (1959) — V-HC-01/02/03, sphere/series solutions.
- C. Lautenberger, dissertation §3.4 + Appendix A (2007) — V-HC-01 params, V-RAD-01, V-CPL-01, V-RX-01.
- S.I. Stoliarov & R.E. Lyon, FSS 9:1141 (2008) — V-GT-01, V-RX-02 (n=2), verification methodology.
- FDS Verification Guide (NIST SP 1018-2) + `firemodels/fds/Verification/Pyrolysis/` — V-HC-02, V-RX-01/02, V-CPL-03, B-02, case-per-chapter guide format.
- A. Khawam & D.R. Flanagan, J. Phys. Chem. B 110:17315 (2006) — solid-state reaction models (V-RX-02).
- H.G. Landau, Q. Appl. Math. 8:81 (1950) — V-CPL-02.
- D.A. Frank-Kamenetskii, *Diffusion and Heat Exchange in Chemical Kinetics* — V-CPL-04.
- J. Crank, *The Mathematics of Diffusion* — V-GT-01.
- MaCFP condensed-phase guidelines (2020) + `MaCFP/matl-db` + MaCFP-4 modeling spec (2026) — B-01/B-03.
- P.M. Muskat-form compressible Darcy (gas reservoir P² solution, standard porous-media texts) — V-GT-02.
