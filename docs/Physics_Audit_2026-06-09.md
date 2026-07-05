# Pyrolysis.jl — Roadmap & Known Gaps

The physics audit (2026-06-09 through 2026-06-11) found and fixed 14 critical/major/minor
groups covering 48 numbered items plus 21 documentation findings. All fix-checklist items
except **item 9** have been resolved. This document is the forward-looking record:
remaining open work, known approximations accepted by design, and missing physics ordered
by impact.

The evidence archive for all resolved items is in `docs/verification_2026-06-10/`.
Verified-sound physics is listed at the bottom — do not "fix" those.

---

## Open checklist item

### Gap-1 (flame feedback) — required before recalibration

Identical to P1.1 in the Missing Physics section below. Gap-2 (blowing) was
**implemented 2026-06-11** (see "Resolved since the audit" below); flame
feedback is now the only gate. **Complete it before any recalibration pass of
the cone examples** — blowing alone raises T_s/MLR and the two partially
cancel.

---

## Resolved since the audit

**P1.2 Blowing correction — IMPLEMENTED 2026-06-11 (opt-in).**
`use_blowing_correction = true` (solve/to_problem_def) rescales every
`ConvectiveBC` in the surface Newton by h·B/(e^B−1), B = ṁ″c_p,g/h, with ṁ″
the true outward gas flux at the face (film law + Darcy venting). Stable
`blowing_factor` (Taylor at |B|<1e-5, overflow-safe e^{−B} form), conductance
rebuilt at every Newton iterate with cp_g(T_s), frozen-B Newton derivative
(first-order consistent, house precedent). Jacobian: `BCThermalOperator` and
`BCConvectionOperator` rebuild ṁ″ from the stencil via the shared
`_bc_boundary_face_fluxes` kernel — support-complete strict verification
passes with blowing on (Fickian/Darcy–Fick/ALE). Diagnostics:
`BoundaryState.blowing_mcp` + blowing-aware `compute_boundary_heat_fluxes`.
Default OFF; examples not recalibrated (gated on Gap-1). Guide §12.3.3, §7.6;
tests in `test/boundary_conditions/test_blowing_correction.jl`.
Known approximation: under ALE the operators omit the z-channel of B (through
the face distance in the flux law), consistent with the existing
frozen-geometry convention — Newton-robustness only, residual exact.

---

## Known approximations (accepted or documented)

These are not bugs; they are physics limitations or accepted design trade-offs, documented
in the guide and/or code. Listed here as a record of what the model does not do.

**Minor 6 — Midpoint Δh is single-interval (accepted by design)**
`src/Materials/components.jl:515-518`: c_p((T₁+T₂)/2)·ΔT; truncation error −c_p″ΔT³/24
≈ −0.8% over 300→800 K for typical curved char c_p. The AD value/derivative-consistency
rationale is sound. Composite midpoint would require N× c_p evaluations in Jacobian-hot
paths for < 1% gain. Documented in the `delta_enthalpy` docstring.

**Minor 9 — Diffusive driving potential ignores ∇P and ∇φ**
`gas_transport.jl:36-47`: J = −(λ/T)∂(ξT)/∂z equals Fick's −ρ_gD∇Y iff ∇(φP) = 0;
the dropped term −Y·∇(φP)/(φP)·λ matters exactly where DarcyFick resolves ∇P and at
sharp char fronts. Inherited ThermaKin convention; documented in `gas_flux` docstring.
Future fix: add the ∇(φP) diffusive correction in DarcyFick mode.

**Minor 28 — Beer–Lambert Jacobian omits z-channel derivatives under ALE**
`operators/radiation.jl:371-406, 428-469` fill only T/ξ columns; the residual uses live
geometry under ALE. The `is_approximate` flag is currently false — this is an unflagged
inexact Jacobian. Impact: Newton convergence slightly degraded under ALE + Beer–Lambert.
Future fix: add the z-channel column derivatives; or set `is_approximate = true`.

**Minor 44 — ALE gas-species mass closure is truncation-order, not telescoping**
The +w·∇ξ ALE correction is pointwise (non-conservative form), so gas-species mass
closure attains truncation-order accuracy, not the 1e-10-level that solid closure
achieves. Documented in `check_mass_conservation` docstring and guide §16.3.3/§16.3.4.

**Minor 47 — All-phase gas volume fraction uses ρ_g(T, P_REF), not local P**
`components.jl:261`: gas volume is over-weighted by ~P/P_REF in pressurized cells,
affecting λ/κ mixing and swelling-gas mixture density (porosity and P = N/φ unaffected).
Documented in `GaseousComponent` constructor.

**Minor 48 — Bruggeman warm start tested on value only**
In AD passes the Dual-typed cache carries stale partials in the warm start, but ≥1
Dual-propagating Newton update always executes before the tolerance break, annihilating
stale partials to O(tol). The pre-update |df| < 1e-20 break on degenerate inputs is the
only true pass-through. Documented in guide §5.15.2.

---

## Missing physics (prioritized)

### P1 — likely biasing current validation cases

**1. Flame heat-flux feedback**
No ignition-triggered flame flux BC exists. Post-ignition cone flames return
~20–30 kW/m² (PMMA; Beaulieu & Dembsey, Stoliarov's ThermaKin cone studies) / ~10–15
kW/m² (wood) — 40–60% on top of a 50 kW/m² exposure. If the shipped examples match data
without it, calibrated kinetics/properties have absorbed a 30–50% net-flux bias.

*Implementation:* low — ignition-triggered `FlameFluxBC` (critical-MLR or critical-T_s
callback) in `src/BoundaryConditions/` + `src/Solver/callbacks.jl`.
**Complete before any recalibration pass of the cone examples.**

**2. Blowing correction h·B/(e^B−1)** — **DONE 2026-06-11**, see "Resolved
since the audit" above. Opt-in `use_blowing_correction`; enable together with
(1) before recalibration. (Original finding, for the record: B = 4 → h/h₀ =
0.075; PMMA peak B ≈ 4, wood B ≈ 0.8–1.3 → 35–52% h reduction; shipped
examples' h = 7.2–10 W/m²K increases B.)

**3. Multi-reactant kinetics** (N_in = 1 enforced, `kinetics.jl:183-185`)
Hard blocker for char oxidation, recondensation, and autocatalysis (Prout–Tompkins fits
structurally unrepresentable). The 3-ary reactant tuple is already anticipated
(`reactions.jl:43-49`).

*Implementation:* medium — touches reactions, kinetics, volume_change, and the
reaction_source Jacobian operator.

### P2 — matters for stated near-term applications

**4. Oxidative pyrolysis / char oxidation with O₂ ingress** (32.8 MJ/kg C)
Dominant at sub-critical fluxes, after flameout, and for smolder; controls the late-time
MLR tail and residue mass in the Douglas-fir case. Requires (3) plus per-species ambient
concentration (`ConvectiveMassBC.ξ_inf` is already per-species after Minor 20).
Upwind Darcy outflow naturally suppresses ingress at high MLR (emergent blowing).

**5. Equilibrium (pressure-dependent) moisture evaporation + recondensation**
Drying is P-independent Arrhenius today; no P_sat anywhere. Vapor recondensation ahead
of the drying front (Perré & Turner; raises local MC 30–50%); at PTW's 1.5–2 atm pore
pressure the boiling point shifts to 112–121 °C (T_sat verified from steam tables) —
unrepresentable with Arrhenius drying. Needs (3) + a smooth r ∝ (P_v − P_sat(T))
closure; P already lives in the cache.

**6. Backing layer / contact resistance at the back face**
Default adiabatic back (package default + wood examples; example 06 already has a
back-loss Robin BC) over-predicts late-time back-half temperatures (~0.7–3 kW/m² for
25 mm wood over 10–20 min). Low cost as a lumped Robin BC with contact conductance +
semi-infinite substrate h(t).

**7. Kinetics expressiveness**
Normalized conversion (α) tracking and the T^n prefactor are absent. Note: arbitrary
real order ξ^n per reactant is already supported and equals ξ₀^n(1−α)^n up to rescaling
A — the gap is α-normalization and DAEM-style stacking, not conversion order per se.

**8. Transport-closure accuracy debt** (documented)
- Wilke viscosity mixing absent (single Sutherland air law)
- Air-default Lennard-Jones parameters for Chapman–Enskog
- First-order upwind Eulerian gas advection (no TVD)
- `permeability_derivatives` (test-only) hard-codes ∂κ/∂T = 0 — the live operator
  correctly differentiates user κ(T) via Duals; this is dead-code hygiene, not a
  production bug

### P3 — niche / future (brief)

- Transient gas energy storage / auto-switch at high P (~0.2%/atm)
- Spectral two-band radiation; in-depth emission for glowing char (the regime where the
  guide points to P1)
- Substrate reflection; scattering (σ_s plumbed in Experimental P1, never exercised);
  view-factor edge falloff
- Knudsen/Klinkenberg slip: κ_eff = κ(1 + b/P); b ≈ 28 kPa at 1 µm pores but
  ≈ 280 kPa at 0.1 µm → up to ~4× κ at 1 atm for sub-µm char
- Stefan–Maxwell multicomponent diffusion (first real customer: O₂ ingress vs
  pyrolyzate efflux); tortuosity factor (φ/τ, factor 2–10) on D and λ_eff
- Liquid free-water capillarity (only above FSP ~30% MC); bound-water diffusion
  (√(Dt) ≈ 0.25 mm per cone test → negligible at fire heating rates)
- Mechanics (char cracking, delamination, intumescence, spalling); anisotropic k
  (wood grain ~2×); melting/dripping; 2D/3D edge burning (documented out of scope;
  ~10–20% late-MLR effect in cone tests)

---

## Verified sound — do not "fix"

These were re-derived from first principles and adversarially checked; they are
correct. A fixing agent should treat deviations from these as regressions.

- **Formulation A core**: with h_rxn at local T, `:standard` drops *only* gas sensible
  storage (~0.2%); matrix-only ρc_p exactly consistent; ξ-based bulk concentrations
  make porosity implicit — no missing φ in the energy equation.
- **Conduction**: face flux −k_f(T_R−T_L)/(d_L+d_R) with distance-weighted harmonic k_f
  is the exact two-resistance form; divergence telescopes; boundary insertion signs
  correct at both faces.
- **S_conv** (sign, units, the two ½ factors, compile-time gas skip, Step ordering —
  matches ThermaKin) and **S_gen** (sign/units; the c_p(T) split ∇Δh_g = c_p,g(T)∇T is
  exact). **Q_rxn = −Σh_r·r_r** (h > 0 endothermic cools; ThermaKin sign-flip caution
  correctly documented).
- **Kinetics machinery**: Arrhenius factorization, C^∞ temperature gates, Σν = 1
  enforcement, ΣS_j ≡ 0 by construction, per-reaction mass/energy consistency.
  **Property derivatives** verified to machine precision.
- **Mixing**: mass/volume fractions + derivatives; mixture density with per-component
  swelling; porosity (intrinsic-density form, clamped); matrix-only ρc_p; conductivity
  parallel/series/weighted on the renormalized solid basis; ε_eff and
  α_eff = Σξ_jα_j (m⁻¹, correct depletion thinning).
- **Chapman–Enskog**: formula, unit conversions (no Pa/atm error), combining rules, all
  8 Neufeld coefficients — matches literature.
- **Beer–Lambert discretization**: exact cell-integrated source
  q_i = I_in(1−e^(−τ_i))/Δz_i; strict conservation by telescoping; backward sweep
  causal (surface = cell N at z = L); live-geometry Δz (ALE-correct). Non-local
  Jacobian suffix-sum factorization exact (T/ξ channels).
- **P1 interior**: ∇·(1/(3(α+σ_s))∇G) = α(G−4σT⁴) re-derived from RTE moments — no
  4π/factor-4 errors; blackbody volumetric emission is the correct LTE closure; Marshak
  G-exchange coefficient ε/(2(2−ε)) correct; Thomas assembly correct; G = 4σT⁴ cold start.
- **Gas diffusion kernel**: exact discretization of −(λ/T)∂(ξT)/∂z; reduces to Fick at
  uniform T; positivity-safe; λ_face distance-weighted harmonic on live geometry.
- **Ideal-gas closure**: P = N/max(φ,1e-12), consistent in both production sites;
  production Jacobian's ∂P/∂ξ_k includes the porosity-coupling term with the right
  sign (via Duals). **Darcy sign conventions** interior + boundary consistent.
  **Surface diffusive flux** series-resistance law re-derived; MLR uses the identical
  formula.
- **Species conservation**: summed species equations reproduce total gas continuity
  with no spurious source from the algebraic pressure closure.
- **ALE core**: θ = Σγ_j(1/ρ_j)dξ_j/dt from the additive-volume constraint
  (reaction-only, no transport contamination); dilation −ξθ_cell counts axial and
  lateral channels exactly once (dm/dt = SV identically, correct for γ = 0 rigid
  char); θ_A/dA/dt/dχ/dt mutually consistent (water excluded by name filter);
  bottom-anchored cumulative w; solids get exactly zero ALE advection;
  advection/dilation mutually exclusive via `is_gas`; production lab-frame flux +
  +w·∇ξ split is an exact identity (no double count).
- **GCL**: volumes recomputed as V = A(t)|Δz| from state z every RHS — discrete GCL
  holds as an algebraic identity; uniform fields preserved exactly.
- **Residual assembly**: every term of the energy and gas-species equations appears
  exactly once with the right sign in BOTH Eulerian and ALE paths; plain ODEProblem,
  no mass matrix, ρc_p re-evaluated every call; state geometry has no redundant
  drifting copy; AD-clean (no Float64 stripping of state-dependent quantities).
- **Surface energy balance**: signs (positive-into-domain), Kirchhoff gray closure,
  F on both emission and absorption, Newton derivative complete for the implemented
  terms (dk/dT correctly absent — k_eff frozen at the cell value consistently),
  q_BC = q_cond at convergence, sound primal/Dual AD design. Boundary sign flips
  and Δz_half = V/(2A) correct; BoundaryConditionSet machinery correct.
  Transpiration: hard-zeroed; `to_problem_def` throws on `use_transpiration_bc=true`.
- **Conservation ledger algebra** (§16.3–16.8) implemented as documented — note: no
  production caller audits closure at runtime; all ledger residual columns are NaN in
  `snapshot_diagnostics_row`.
- **Correctly negligible** (estimates verified): gravity in Darcy (~0.14 Pa vs kPa),
  Soret/Dufour (leading drift already in the ξT-gradient flux), local thermal
  non-equilibrium (pore equilibration ~0.02–0.13 ms).
- `_upwind_gradient_distance` boundary full-Δz is correct and inert — no 0.5 factor
  needed.

---

## Previously removed claims — do not reinstate

Two originally-reported sub-claims failed adversarial verification and are **not**
defects:

1. The claim that `radiation_types.jl`'s "+O(n²) dense fill" memory comment was wrong.
   The comment is **correct**: the default `:full` sparsity profile emits NC·n(n+1)/2
   slots for the Beer–Lambert coupling; the O(n) `:ale_reduced` profile drops the
   coupling.
2. The claim that production Darcy faces use an *unweighted* harmonic mean. Production
   (`finite_volume.jl`) is distance-weighted; the unweighted form exists only in the
   deleted dead `darcy_velocity_face` family.
