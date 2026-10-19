# Mesh framework: the Phase 1 review

**Scope.** An adversarial review of all Phase 1 work on branch `mesh-framework-phase1` at `8484153` (the diff `17bc477..8484153`): the nine new modules, their tests, the gate scripts and tables, `docs/Mesh_Framework_Phase1_Gate.md`, plan §12, the suite routing and the old-core edits, against the design (`docs/Mesh_Framework_Design_2026-10-17.md`), the plan and the repository rules. The review ran as nine area reviews (states and meshes, the conduction/reaction/boundary kernels, transport/radiation/closures, assembly/deactivation/outputs, the colored Jacobian and integration, the tests, the documents' traceability, the repository rules, and measurements of the gate's soft spots), each returning findings with `file:line` evidence, verified against the code and, where a claim was numerical, re-measured in persistent Julia sessions on the same machine (load 17 to 21 of 22 logical CPUs throughout; the owner's other runs). Its fixes are on branch `phase1-review`; plan §12 records them as D38 to D44 and C32 to C36.

**Verdict.** No finding changes a verification result or a gate verdict. Two defects of the solve's bookkeeping were real (a terminated solve lost its last state; a state saved at an event's own time was reconstructed under the mask before the event), both fixed. Three claims of the gate document were wrong or unsupported and are corrected: the cause given for V-MMS-03's temperature order (the limiter; measured to be the upwind-biased face temperature's pre-asymptotic term), the statement that C31's remedy is upstream (a package-side remedy exists and is measured), and a profile share and a wall time that no committed table holds. The old core's bit-identity under the §2.3 edits, asserted but never committed as evidence, was recomputed on both trees and is identical (C35). The rest is robustness, test coverage, wording and traceability.

## 1. Findings

Severity: `blocker` none; `major` a real defect, an unsupported gate claim or a missing proof; `minor` robustness, coverage, pins, wording; `note` recorded, no action. Status: `fixed` (on `phase1-review`), `followup` (yours to do next, listed in the handoff prompt), `owner` (a decision), `upstream`, `none`.

| ID | Sev. | Where | Finding | Status |
|---|---|---|---|---|
| D-1 | major | `src/Deactivation/events.jl` `replay_events!` | A state saved at an event's own time is the post-event state (the integrator applies a discrete callback before `savevalues!`), but the replay applied events strictly before `t`: the dead cell showed a temperature instead of `NaN`, `live_cells` was one too many, the surface series read a massless cell. | fixed (D38) |
| D-2 | major | `src/Integration/solve.jl` `integrate` | A solve ended by deactivation or `min_thickness` saved no state at or after the termination unless `saveat` named that time (`saveat = [0, 20]` on a burn-out at 0.59 s returned the initial state alone); the old core behaves the same. | fixed (D38) |
| X-1 | major | `docs/Mesh_Framework_Phase1_Gate.md` item 3, plan §12.6, `V-MMS-03.jl` header | The gate's cause for V-MMS-03's order (van Albada clipping) is refuted: the species converge at 2.00, the temperature's 1.63 to 1.92 (to 512 cells) is a pre-asymptotic O(h³) term of the upwind-biased face temperature; a centred face temperature gives 1.98 to 2.00 at 4 to 6× the error. | fixed (documents, C33); owner (face temperature) |
| M-1 | major | gate item 2, plan C31 | "The remedy is upstream" is wrong: a LinearSolve algorithm of the package's own keeps KLU's analysis on the prototype's structure (prototype measured: 90 → 1 analyses, temperatures within 6e-11 K, ratios 2.2/2.6/1.6/2.1 → 1.7/2.0/1.1/1.4). | fixed (documents, C32); owner (adopt) |
| R-1, K-1 | major | plan §12.1 M3, commit `73200f8` | The bit-identity of the old core under the §2.3 edits rested on a session hash; nothing committed re-checks it. | fixed (recomputed, C35); followup (commit the script) |
| G-1 | major | gate item 5 | "16 s against 3.6 s" on V-BO-02 contradicts `case_accuracy.csv` (28.4 against 3.73 s). | fixed |
| G-2 | major | gate item 1, plan M9 row | "188 quantities" per form was in no table (a static count gave 183); counted in a run: 33 cases, 188 quantities, 11 convergence checks, so the gate is right and now says so. | fixed (confirmed) |
| G-3 | major | gate item 2 | The 3 to 22× accuracy gap at equal tolerance is attributed to the old core's smoothed estimator; `audit_breakdown.md` shows the new core's smoothed `KenCarp4` misses the same target, so the cause is open. | fixed (softened) |
| G-4 | major | gate item 2, C31 | "A profile puts 72 % of that solve in fresh factorizations" is uncommitted; the committed substitution saves 12 to 38 % of the solve. | fixed |
| TS-1 | major | `src/Integration/solve.jl` callback set | A read-only user callback raised the derivative-discontinuity flag on every step, which the policy's flag rule cannot clear from another callback, so under `fbdf_policy` it restarted the method at order one every step (34 → 83,735 steps). | fixed (D40) + test |
| TS-2 | major | tests | No test covered a terminated solve on the normal output path (see D-2). | fixed (test) |
| TS-3 | major | tests | Plan §6.5 item 8 (repeated solves bit-identical under KLU and UMFPACK) had no new-core test; the property holds. | fixed (test) |
| TS-4 | major | tests | The three checks of design §3.8 are tested at unit level only; no solve is driven into a constraint rejection, the give-up path or the post-event admissibility refusal. | followup |
| J-1 | major | tests | The `:forced_dt` restart and KenCarp4's `:keep` flag are never exercised through a real deactivation event in `Pkg.test()`. | followup |
| T-1 | major | `radiation.jl`, `setup.jl` | Radial Beer–Lambert would create energy (intensity rows on `I`, sources on `I A`); no case uses it. | fixed (refused, D45); owner (C37) |
| D-3 | major | `events.jl` fixed-mesh transfer | The dead cell's pore gas is added to the receiver's unchanged pore volume (a pressure step per event under Darcy–Fick); nothing ledgers the pore volume or checks the pressure. | owner (C34, C4) + followup test |
| X-3 | major | gate item 8 | V-BO-02's event-driven 26 K persists through whole park intervals (median 0.85 K), and before any event its temperature error is 0.065 K, nine times the budget unit; the gate said neither. | fixed (documents); owner (open) |
| D-4, S-5 | minor | `setup.jl`, `tolerances.jl` | A cell starting without condensed mass had the tolerance scale (total mass or unit density) as its deactivation reference, so it was below the threshold at t = 0. | fixed (D39) |
| D-6 | minor | `Outputs/solution.jl` | A radial mesh from r = 0 made every `:bottom` output `0/0`. | fixed (D42) |
| S-2 | minor | `caloric.jl` `_invert_primal` | The bracket rests on a slope bound that holds only up to a law's variation between its 1 K samples; a collapsed bracket was returned without a residual check. | fixed (D43) |
| S-3 | minor | plan D28 | The regularization boundary is `4 RHO_CP_FLOOR`, not 2 (the guaranteed capacity uses half the law's minimum). | fixed (plan) |
| TS-5 | minor | `test/runtests.jl` | The three-way `Pkg.test()` ran `test/newcore` three times (no file there reads the energy form); `test_analytical.jl` (old-core solves) runs under the new-core label. | fixed (D44, the first); followup (the second) |
| TS-6 | minor | `thermal_boundary.jl` | The surface Newton received the raw cell temperature, so wild trial states ran it to its cap with warnings. | fixed (D41) |
| J-3 | minor | `checks.jl`, `solve.jl` docstrings | A domain rejection shrinks the step by the controller's `qmin` (one fifth), not half; the Newton-failure cap is inert under Rodas5P. | fixed |
| J-7 | minor | `checks.jl` | A non-finite constraint residual passed the domain test (`NaN > 1` is false). | fixed |
| K-2 | minor | `Materials/reactions.jl` | The `OneSidedLimiter` docstring misstated the disabled behaviour at order zero. | fixed |
| K-3 | minor | `test_kinetics.jl:73` | The calibration guard is exact by construction; the test pinned it at 1e-14. | fixed (`==`) |
| K-5 | minor | `Physics/kinetics.jl` | Stale comments after the limiter dispatch. | fixed |
| S-14 | minor | `Meshes/Meshes.jl` | The docstring put the metrics on `Face`. | fixed |
| T-3, T-7 | minor | `relative_flux.jl`, `gas_transport.jl` | A wrong claim about opposite-sign van Albada slopes; a wrong orientation word. | fixed |
| T-4 | minor | tests | No radial Beer–Lambert, `:bottom` mass-flux sign or two-cell face-flux value test. | followup |
| R-2, R-3, R-4, D-7, TS-13 | minor | docstrings, comments | Conflict labels and historical wording in code comments; "agent summary" in the design's §1.3 title. | fixed |
| R-6 | minor | plan C20 | The guides' stale references also cover `boundary_fluxes.jl` and the `Material` signature. | fixed |
| D-8 | minor | plan §3.4 | The `S_energy` comment carried the formation term the row assembly adds. | fixed |
| G-5, G-6, G-7, G-9, G-12, G-13 | minor | gate document, plan §12.6 | Item 1's scope attribution; the 7e-6 K caveat; the 256-cell session number; minima against medians; superseded §12.6 bullets; which item 2 verdicts lie within noise. | fixed |
| G-8, G-10 | minor | gate document, scripts | Session-only numbers (`Pkg.test()` counts, test counts, report-mode counts, C30's kelvin differences) are now labelled; the gate scripts record no commit. | fixed (labels); followup (commit stamp) |
| G-11 | minor | `V-BO-02.jl` | The enthalpy-form step budget of 4,000 extends a ruling that covered 3,000. | owner |
| X-2 | minor | gate item 8, `case_accuracy.jl` | V-CPL-05's references stop at depletion, but over the common saved times the error is 0.2 of the budget unit; the script compares nothing when the reference stops early. | fixed (documents); followup (script) |
| K-4 | minor | `gas_transport.jl` `outward_darcy_velocity` | The convective pressure vent is gated on `mob > 0`, a jump in a state-derived quantity. | owner (C36) |
| K-6 | minor | `surface_balance.jl` | The surface Newton's hard clamp and warnings inside the residual are plan §3.8's choice but not a recorded departure from design §3.10. | fixed (recorded, C36); owner |
| S-1 | minor | plan §3.2 | The "two functions" contract is not what is built: Deactivation, the lateral-row Jacobian, the domain test and the thickness terminator read the raw state by index. | followup (plan text) |
| S-4 | minor | plan §3.2 | `T_warm` carries the last evaluated temperature, not the last accepted; harmless, the plan sentence is wrong. | followup (plan text) |
| S-6 | minor | tests | `abstol_vector` and the moving-mesh branch of `derive!` have no direct test. | followup |
| J-2 | minor | tests | No containment run for radial × Beer–Lambert, radial × Darcy, radial × relative flux, or after `deactivate!` with radiation; the refusal of deactivation under `w_override` is untested. | followup |
| J-5 | minor | `setup.jl` `_consistent_motion!` | A failing velocity sweep after an event throws from inside the callback instead of ending the solve with `:Failure`. | followup |
| TS-7, TS-8, TS-9, TS-10 | minor | tests | A column-scaled tolerance in `test_radiation.jl`; the ledger order pinned at 2 with 2.17 observed; the five loosest and five most fragile pins. | followup |
| TS-11 | minor | tests | Untested behaviours: radial deactivation, Beer–Lambert through a deactivating solve, the switch's effect on the rates, `:interior`/`:burned_out` through a solve, `UMFPACKFactorization`, several refusals. | followup (two added) |
| TS-12, R-7 | minor | gate scripts | `Dates` and `Printf` are not in `test/Project.toml`; the scripts run outside `Pkg.test()`. | followup |
| R-5 | minor | design document | Five companion documents it cites are gitignored. | owner |
| J-4 | minor | `solve.jl` | `progress` and `jacobian` are accepted and ignored for the routing's sake: a keyword shim in a new module. | owner |
| D-5 | minor | `residual.jl` `buffers` | A dual `t` with a `Float64` state (forward-mode `tgrad`) would throw; unreachable with the shipped policies. | followup (document) |
| D-9, D-10, J-6, J-8, K-7, K-8, K-9, S-7 to S-13, T-5, T-6, TS-13, X-4, G-14, G-15, R-8 | note | various | See §2. | none |

## 2. Findings by area

### 2.1 States and meshes (S)

Verified correct: the interleaved layout and every index function tile the state exactly; `differential` and the mass matrix match design §3.2; `positive_part`, `smooth_floor`, `smooth_clamp` and the capacity switch are C¹ as stated; the extended caloric relation is strictly increasing with slope at least the floor; the inversion's convergence criterion and implicit derivatives are as §12.5 says; D32 and D33 hold; one datum `T_REF` for every component; the temperature-form row is design §3.2's expansion term by term; every tolerance scale is positive; the geometry writer is unique and C16 holds algebraically.

S-2 (fixed, D43): the guaranteed capacity is half the law's sampled minimum, a bound only up to the variation between samples; the collapse exit of `_invert_primal` now returns only when the residual is small and otherwise doubles the bracket away from the datum. S-3 (fixed): D28 now states the `4 RHO_CP_FLOOR` bound. S-5 (fixed, D39): the deactivation reference is the initial condensed mass itself. Notes: S-7 (the Carman–Kozeny cap and the volume-fraction cut inside Properties are hard pieces the new core reaches under Darcy–Fick), S-8 (two homogeneous tuples are indexed in loops; harmless), S-9 (`sensible_enthalpy(::TableProperty)` recomputes the datum integral per call), S-11 (the enthalpy tolerance is 100× tighter than the temperature form's near `T_REF`), S-12 (`@boundscheck` trait checks are elided under `@inbounds`), S-13 (`live_face` and the volume merge support deactivation from the exposed side only, as the plan says).

### 2.2 Conduction, reactions, boundary kernels and the old-core edits (K)

Verified correct: the one-sided limiter is exactly zero at and below zero for every order, C¹ at the join, bit-for-bit equal to the old limiter above it; the finite base is formed in log space as D4 says; the formation term and the reaction heat enter the two energy forms exactly as design §3.2; the reactivity switch is C¹, per cell, inert without deactivation; the face conductivity, gradients and boundary flux ownership; the surface Newton's code is byte-identical to the old core's (moved); every closed form and every sign along the outward normal; every old-core edit is additive and the default limiter is the old one; `scripts/mesh_phase0/` untouched after M0.

K-1 (with R-1; C35): the fingerprint was recomputed on `17bc477` and on the review tree in the same Julia: the hash of 3,072 reaction rates and heats (orders 0, 0.5, 1, 2; thresholds 1e-3, 0.1 and 10 kg/m³; limiter on and off; 16 concentrations from −1 to 1e3 kg/m³ at four temperatures) and the `(T, ξ, t)` hashes of eight old-core solves (V-RX-01, 02a, 02c, 02d, 03, 04a, V-CPL-01, V-CPL-04 at their CI resolutions) are identical on `17bc477` and on the review tree (rates `0xc8fea394cd668f5a`; every solve hash equal), so the §2.3 edits left the old core bit-identical. K-4 (owner, C36): `outward_darcy_velocity(::ConvectivePressureBC)` is `mob > 0 ? h_P (P − P_ext) : 0`, a jump. K-6 (recorded, C36): the surface Newton's clamp and warnings. K-8 (note): the "full tensor" is a 1×1 signature in 1D. K-9 (note): `ε` and `a` are the cell's at its trial temperature, not at `T_s`, as on the old core. K-7 (note): `Contributions.T_face` is not reset; only live faces are read.

### 2.3 Transport, radiation, closures (T)

Verified correct: the flux form conserves (one scatter per face, mass and energy sharing `F`); the `:bottom` normal and every boundary flux along the outward normal; `outward_darcy_velocity` reproduces the old `boundary_darcy_velocity` on every condition; the outflow-only vent and the first-order upwind are the old core's (C21); the Fick coefficient, mobility mean, face porosity and the centred face temperature of the gas enthalpy; the boundary enthalpy once, at `T_s`, after the gas kernel and before the relative kernel; Beer–Lambert on a planar mesh (order, rows, telescoping source, consistent sweep, the exposed face after deactivation); the closure's rows, sweep and rank check with radial areas; the lateral rows (a C¹ clamp, not a hard one). The transport reviewer's own hypothesis for V-MMS-03 (the van Albada slope near the travelling extremum, with the thermal Péclet number explaining why the gas is immune) was the one the measurements of §2.9 tested: replacing the slope by the unlimited central one gave 1.54, 1.78, 1.90, so it is not the dominant term; the centred face temperature is.

T-1 (fixed by refusal, D45; physics the owner's, C37): on a cylinder or sphere the intensity rows attenuate the intensity while the sources read the face areas, so a transparent medium would absorb `I (A_up − A_down)` per cell; no case uses the combination; `setup` now refuses it. T-2 is K-4 (C36). T-3 (fixed): the `van_albada` docstring claimed the slope nearly vanishes for opposite signs; it reaches 0.21 of the larger gradient. T-7 (fixed): a wrong orientation word in the gas kernel's docstring. T-4 (followup): no radial Beer–Lambert test, no residual-level `:bottom` mass-flux sign assertion, no two-cell value check of the new face flux against `gas_flux`. T-5 (note): `mass_bc_flux(::DiffusiveMassBC)`'s `ρ <= 0` branch is a cut reachable only in an empty trial cell. T-6 (note): `consistent_velocities!` divides by the row node's face area, guarded by `rank_check`, outside the residual.

### 2.4 Assembly, deactivation, outputs (D)

Verified correct: the assembler order; both energy rows; face scatter with each face's own area; dead rows; the cascade loop; the merged temperature; the volume rule equals design §3.5 term by term; admissibility after an event; every output name's semantics against the old core's (`mass_loss_rate`, `MLR_gauge`, `total_mass`, `gas_generation_rate`, `thickness_history`, the surface decomposition, `volumetric_absorption`, `radiation_transmitted`, `z`, `park_times`, `live_cells`); the reconstruction is a pure function of the saved state.

D-1 and D-2 are fixed (D38) with a regression test (`test_deactivation.jl`: a four-cell burn-out with `saveat = [0, 20]` ends on its last state with every cell `NaN`, and a state saved exactly at the first event's time carries the post-event mask). D-3 is the owner's (C34). D-9 (note): the dead cell's caloric regularization is dropped without an event entry. D-10 (note): the transfer's inversion warm-starts at `T_REF`.

### 2.5 Colored Jacobian and integration (J)

Verified correct: the pattern matches every kernel's reads (the 4,780 entries of `ale` at 100 cells reproduced by hand), the exposed face after deactivation included; D17's dense rows are exact; the coloring and scatter; the setup's refusals and the gas-model rule; `NoInit` with the sweep; the policies against design §3.8; the domain test runs on every adaptive proposal; the callback mechanics in the installed DiffEqBase.

On C31: W's structure is the numerically nonzero part of J (SparseArrays' `_broadcast_zeropres!` stores only nonzero results); `Wfact_t` is not a route because `calc_W!` returns `new_W = false` after it and the Rosenbrock step then never marks the linear cache fresh (a stale factorization, not a slow one); `KLUFactorization(check_pattern = false)` still re-analyses on an nnz change and is unsafe otherwise; LinearSolve's nonstructural-zero reduction is inert for OrdinaryDiffEq-built caches. The package-side remedy (M-1, C32) is a `LinearSolve.AbstractSparseFactorization` subtype whose `init_cacheval` captures the all-ones full-structure W that `build_J_W` hands `init`, and whose `solve!` scatters each fresh W into it and delegates to an inner KLU cache; the prototype is `review_structured_klu.jl` in the handoff's scratch directory (not package code). Measured on the gate's decks in one session, both cores back to back, warm minima of three, load 18 to 19:

| deck | new, default [s] | new, structured [s] | old [s] | ratio, default | ratio, structured | ΔT [K] |
|---|---|---|---|---|---|---|
| eulerian 100 | 0.0313 | 0.0346 | 0.0388 | 0.81 | 0.89 | 0 |
| ale 25 | 0.0292 | 0.0229 | 0.0133 | 2.20 | 1.72 | 5.5e-11 |
| ale 100 | 0.117 | 0.0919 | 0.0451 | 2.60 | 2.04 | 5.9e-11 |
| ale_darcy 25 | 0.0536 | 0.0366 | 0.0330 | 1.62 | 1.11 | 7.1e-6 |
| ale_darcy 100 | 0.279 | 0.190 | 0.132 | 2.11 | 1.44 | 5.2e-11 |

On `ale` 100: 90 steps, one KLU analysis against 90, identical step counts and statistics. With it, three of the six matched-accuracy decks are under 1.5×; the two plain `ale` decks still fail, as `w_formation.md` predicted. J-1, J-2, J-5 are followups; J-3 and J-7 fixed; J-6 (note): `rhs_evaluations` omits the domain test's residual per proposal and the dense rows' residual per Jacobian; J-8 (note): the velocity tolerance scales as `1/dt`.

### 2.6 Tests (TS)

Verified correct: every refusal of §12.2 but one is tested; the §3.2, §3.5, §3.6 and §6.5 lists are covered except item 8 (now added) and the tensor test that has no kernel to run on; the dense oracle is pattern-free; the old-test modifications are import moves only; the test-count arithmetic (4,113 + 7,526 = 11,639; 8,212 − 7,526 = 686 per form).

TS-1 (fixed, D40, test added): the policy's flag rule now runs inside each discrete callback the user passes, after its `affect!`. The review's first fix, the user callback placed before the step observer, was inert, as the queued test showed (`fbdf_policy`: 57 steps became 999,996, the step limit): the installed DiffEqBase raises the flag before every discrete callback's `affect!`, captures each callback's own verdict right after it, and folds the verdicts whatever their order. TS-3 (test added): two solves bit-identical under KLU and UMFPACK. TS-5 (fixed, D44): the enthalpy-form run skips `test/newcore`. TS-6 (fixed, D41). TS-4, TS-7 to TS-12 are followups; the loosest pins are 400 to 29,000× above their observed values (`test_mesh_motion.jl:347`, `test_policies.jl:225`, `test_relative_flux.jl:88`, `test_solve.jl:117`), the most fragile the ledger order `>= 2` against 2.17 observed and the exact color counts fitted to the code. TS-13 (note): `missing_entries` is blind to couplings that vanish at a state (the limiter at the zero and negative states).

### 2.7 The documents and the gate tables (G)

Every number of the gate document was traced to a committed table or marked a session number; the matched-reltol rule, the equal-tolerance and matched ratios, the degenerate-deck ratios, items 3, 6, 7, 8, 9 and the provisional-choices numbers reproduce from the tables; every §12.5 constant matches the code; every C-number cited in §12 is defined; the gate scripts write only under `measurements/`. G-1 to G-9, G-12, G-13 are corrected in the gate document; G-10 (no commit stamp in the tables; two scripts edited after their tables) and G-14 (V-RAD-01's constraint residual prints exactly zero; a run of the case gives 63 accepted steps, 0 rejected, 0 domain rejections and 0.0 again, so a probe of the intensity-row residual at accepted states is queued as `jobs1/09_vrad_residual.jl`) are followups; G-11 is the owner's; G-15 is a wording of D37 (the margin bounds the total L2, not the time error).

### 2.8 Repository rules (R)

Clean: no AI attribution in the 42 commits, the PR body, the workflow or the diff; no dates, milestone names or historical words in any changed `src/` file; no dual path or shim in `src/` apart from the two ignored keywords (J-4); the move is byte-identical apart from a docstring path and a `using` line; every old-core edit is in the plan's list; all Manifests consistent; the formatter check clean; every commit signed with a body. Fixed: the conflict labels and historical wording in comments (R-3, R-4), the design's §1.3 title (R-2), C20 (R-6). R-5 (the design cites gitignored companions) and the file names `leftover1_jacobian_cost.jl`/`leftover2_integrators.jl` are the owner's. R-8 (note): three commits mix themes.

### 2.9 Measurements (X and the coordinator)

X-1 (V-MMS-03), X-2 (V-CPL-05) and X-3 (V-BO-02 event timing) are in §1 and in the gate document; the server state after the override experiments is tainted and was not reused. X-4 confirms the C31 mechanism in the installed sources. The coordinator's measurements: the C31 structure numbers reproduce exactly; the structured-KLU prototype (§2.5); a terminated solve's saved states on both cores (D-2).

## 3. Remediation plan

1. **Fixed on `phase1-review`** (one theme per commit, listed in §4): the solve's bookkeeping (D-1, D-2), the deactivation reference (D-4), the user callback's flag (TS-1), the surface Newton's trial temperature (TS-6), the inversion bracket (S-2), zero-area outputs (D-6), the domain test's non-finite residual (J-7), the test-runner skip (TS-5), the regression tests (D-1/D-2, TS-1, TS-3), the exact calibration guard (K-3), docstrings and comments (K-2, K-5, S-14, J-3, R-2 to R-4, D-7), and the gate document and plan §12 corrections (G-1 to G-9, G-12, G-13, X-1 to X-3, M-1, S-3, R-6, D-8, C32 to C36, D38 to D44).
2. **For the owner** (§7): adopt the structure-keeping linear solve (C32) and re-measure item 2; V-MMS-03's face temperature (C33); the fixed-mesh transfer of the dead cell's gas (C34, C4); V-BO-02's enthalpy-form step budget (G-11); the non-smooth pieces of C36; whether V-BO-02's pre-event error (X-3) is investigated now or in Phase 2; the design's gitignored companions (R-5); the ignored keywords (J-4).
3. **Followups** (the handoff prompt): commit the fingerprint script (C35); the tests of TS-4, J-1, J-2, S-6, TS-11; `case_accuracy.jl` comparing over common times (X-2); a commit stamp in the gate scripts (G-10); `Dates`/`Printf` in `test/Project.toml` (TS-12); `test_analytical.jl` out of the new-core domain (TS-5); J-5; the plan sentences of S-1 and S-4; TS-7/TS-8's criteria.
4. **Upstream**: the structure-keeping `jacobian2W!` is still worth proposing (C31), but it is no longer the only route; the text for the issue is in the handoff prompt.

## 4. What was implemented, with gates

Source: `src/Integration/solve.jl` (save_end, the user callback's flag rule, docstrings), `src/Deactivation/events.jl` (replay rule, reference guard), `src/Integration/setup.jl` (condensed reference), `src/Integration/checks.jl` (non-finite residual, docstring), `src/Kernels/thermal_boundary.jl` (trial temperature), `src/Outputs/solution.jl` (zero-area faces), `src/States/caloric.jl` (bracket widening), `src/States/derived.jl`, `src/States/guards.jl`, `src/Meshes/Meshes.jl`, `src/Materials/reactions.jl`, `src/Physics/kinetics.jl` (docstrings and comments). Tests: `test/newcore/test_deactivation.jl` (replay rule, the terminated-solve testset), `test/newcore/test_policies.jl` (observer inertness, determinism), `test/newcore/test_kinetics.jl` (`==`), comment edits in four files, `test/runtests.jl` (D44). Documents: the gate document, plan §12 and §3.4, the design's §1.3 title. The old core's source is untouched (the two edited old-core files, `src/Materials/reactions.jl` and `src/Physics/kinetics.jl`, changed in docstrings and comments only; the fingerprint of C35 covers the tree).

Gates run on the review tree (persistent sessions, load 17 to 21):

- Formatter: `julia --project=scripts scripts/format.jl` applied (it rewrapped `scripts/verification/run_reports.jl` and nothing else of note), then `--check`: "Format check passed: all files already formatted."
- Old core: `src/Physics/kinetics.jl` and `src/Materials/reactions.jl` changed in comments and a docstring only; the fingerprint of C35 (rates and eight solves) is identical on `17bc477` and on this tree.
- New-core measurements that back the fixes: the structured-KLU prototype (§2.5); a four-cell burn-out with `saveat = [0, 20]` now ends on its last state (checked by the new testset, written and queued); the C31 structure numbers reproduce; the QoI count (188) and V-RAD-01's statistics from a run of every routed case.
- Targeted test files: `test_deactivation.jl`, `test_policies.jl`, `test_states.jl`, `test_kinetics.jl` (server 1) and `test_solve.jl`, `test_boundary_conditions.jl`, `test_mesh_motion.jl`, `test_residual_jacobian.jl`, `test_energy_ledger.jl`, `test_radiation.jl`, `test_gas_transport.jl`, `test_relative_flux.jl` (server 4) were queued, with `NEW_CORE_CASES` in both forms after them; the session reached its limit before they finished. Their logs are `/tmp/claude-1000/-home-gbellamy-Documents-git-Pyrolysis-jl/65209d8c-094c-4cf2-8307-4394f4e714e7/scratchpad/jobs1/08_tests_a.log` and `/tmp/claude-1000/-home-gbellamy-Documents-git-Pyrolysis-jl/65209d8c-094c-4cf2-8307-4394f4e714e7/scratchpad/jobs4/01_tests_b.log`, `02_cases_T.log`, `03_cases_H.log` (each written when its job ends; a `.done` file marks it). The first task of the handoff prompt is to read them and fix anything they show before committing.
- Not run: Aqua and JET (`julia --project=test/quality test/quality/aqua.jl`, `jet.jl`); the containment suite; the three-way `Pkg.test()`.

## 5. What was not implemented, and why

- The structure-keeping linear solve (C32): a new numerical component in the solve path; the owner decides the default (§7, question 1). The prototype and its measurement are complete.
- V-MMS-03's face temperature (C33): a discretization choice with a trade-off (order against error); the owner's.
- The fixed-mesh transfer (C34): C4 is the owner's.
- The followup tests (TS-4, J-1, J-2, S-6, TS-11) and script changes (X-2, G-10): time; they are specified in the handoff prompt.
- The fingerprint script as a committed test: it needs the pre-Phase-1 tree to prove identity, so it is a script with a recorded digest, not a `Pkg.test()` item; the handoff prompt has it.

## 6. The three-way `Pkg.test()`

Not run: the session that performed the review reached its usage limit after the fixes were written and the targeted tests queued. The handoff prompt runs it on the final commit (old core about 85 min; `PYROLYSIS_CORE=new` temperature form about 60 min; the enthalpy form now the verification suite alone), the three in parallel in detached worktrees with the development servers stopped, and records the counts here and in plan §12.1.

## 7. Owner questions

See the review's final message and `docs/HANDOFF_Phase1_Review_Followup.md` (gitignored), which carry the same list with the options and recommendations.
