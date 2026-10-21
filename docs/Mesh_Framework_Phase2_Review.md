# Mesh framework: the Phase 2 review

**Scope.** An adversarial review of all Phase 2 work on branch `mesh-framework-phase2` at `8166e77` (the diff `bb2860b..8166e77`: 70 files, 9,504 insertions, 849 deletions; 96 commits over `mesh-framework-phase1`): the correction and its binding, the block-and-band rule, event location, the structured factorization, the restart policy, the old-core edits (the surface Newton's trial clamp, the conductivity skip), the tests, the eight V-BO cases, the seven gate scripts and their tables, `docs/Mesh_Framework_Phase2_Gate.md` and plan §12 to §14, against the design (`docs/Mesh_Framework_Design_2026-10-17.md`), the plan (`docs/Mesh_Framework_Phase2_Plan.md`) and the repository rules. The review ran as one coordinating read of every changed source, test and case file, plus two independent reads with no shared conclusions: the gate scripts and tables against the gate document and §7.3, and the deactivation, event and correction code against its contracts. It ran the format check, Aqua and JET on `8166e77`; it ran no solve, because the machine carried fourteen OpenFOAM processes at a load of 14.5 throughout, which is also why the gate's item 8 could not be re-measured here.

**Verdict.** No finding changes a verification result, a test count or a gate verdict. One defect of the kernels is real: after a merge, the relative-flux stencil on the face just inward of the outermost live node reads the dead cell beyond it (its zeroed masses, its frozen temperature and a distance built on its frozen node) whenever that node's correction is nonzero, which every class B deck meets after its first event; the effect is bounded by the limiter and face-conservative, so no ledger term is wrong, but live rows depend on masked state, the two energy forms part at that face beyond round-off, and no test reaches the path (K-1). The gate document misreads one table (the largest residual is 0.993, not 0.998; G-2), overrides by hand a ladder column whose criterion passes a run with no invocation (S-1), rests item 9's recommended resolution on numbers no committed table holds (G-3), and, with every table and plan §12, stamps a commit that does not exist in the repository (G-1). The phase is not closed at its gate: item 8 is recorded and undecided, and on the recorded numbers the phase's own stop rule of §5.7 and §7.4 is reached on five rows (X-1); items 1 and 9 carry four and two owner exceptions still to be granted (C91 to C94, V-BO-06 at 3e-2); MF waits for the upstream release; and the pull request still describes the planning milestone (P-1). What the phase set out to prove, one invocation per layer with the layer's QoIs intact, holds at the default threshold on every deck but V-BO-07; what it set out to keep, the old core's order of wall time, it has not shown.

**Contents.** 1 Findings · 2 Findings by area · 3 What remains · 4 What this review ran · 5 What the followup implemented, with gates · 6 What was not implemented, and why · 7 The three-way `Pkg.test()` · 8 The owner's answers

---

## 1. Findings

Severity: `blocker` none; `major` a real defect, an unsupported gate claim or a missing proof; `minor` robustness, coverage, pins, wording; `note` recorded, no action. Status: `followup` (to do next, §3), `owner` (a decision), `none`; after the followup (§5), `fixed` with its record and commit.

| ID | Sev. | Where | Finding | Status |
|---|---|---|---|---|
| K-1 | major | `src/Kernels/relative_flux.jl:135-158, :179-181`; `src/Closures/correction.jl:265-273`; `src/Deactivation/events.jl:393-397, :427-437` | After a merge into receiver `r`, `_correction!` zeroes `c[r+1]` (cell `r+1` is dead) but not `c[r]`, whose node lies between two live cells of a layer with gain. When cell `r` is thicker than cell `r−1`, the usual state after it absorbed the run's volume, `c[r] > 0`, so the face at node `r` has `u = −c[r] < 0` and `_stencil` takes side 2: `U = r`, `D = r−1`, `X = far[2, f] = r+1`, the dead cell. `_face_value` then reads `ξ(X) = 0`, `T(X)` frozen at the merge (`T_m` in the temperature form, the inversion of a zero inventory in the enthalpy form) and `d_XU` built on the dead node `x[r+2]`, which stays where `_enclose!` left it while the live node `r+1` regresses. Live rows read masked state, the two forms differ at that face's enthalpy flux beyond round-off, and `derived.jl:145-146` ("the kernels skip every face and patch beyond" the live cells) does not hold for this kernel. Unreachable where the receiver's layer has zero gain (V-BO-01/03/06/08, whose receiver is Kaowool; V-BO-04, whose survivor is cell 1); reachable on every class B deck after its first event (V-CPL-02, V-CPL-02b, `ale` and `ale_darcy` with deactivation). Bounded: the van Albada slope tends to the downwind gradient when the upwind one is large, and the flux is face-conservative, so no ledger term is wrong and no QoI verdict changes. The pattern (reach 2) holds the dead columns, so the colored and dense Jacobians agree. The masked-rows test keeps every cell above `r_hi`, so `c ≡ 0` there and the path is untested. | fixed (D90, `46059c2`, with its test); the path is not reached on V-CPL-02/02b, whose digests stay (C97) |
| G-1 | major | `test/newcore/gates/measurements/*.md:3` (all seven tables), gate §Scope, plan §12.1 M6c and M7 | Stamped "commit 5f08134", which does not exist in the repository: `stamp.jl` read HEAD before the commit queue was re-created by `git am -S`, which changed every hash. The signed commit with the same tree (`543bbb8`) is `78ce699`, written nowhere in a tracked file; the mapping lives only in the gitignored `docs/refs/phase2_commit_queue/TREES.txt`. The tables do describe HEAD's code (`git diff 78ce699..HEAD` touches documents and tables only), so no number is in doubt, but D52's purpose, a reader reproducing a table at its commit, fails. | fixed (D98, `25dac70`; the mapping in the gate and plan §12.1, C95) |
| G-2 | major | gate §Summary item 6 and §6, plan §12.1 M6c | "The largest accepted scaled residual 0.998 under C52's scale, none above one under Phase 1's": the ladder's maximum is 0.993 (`burnout_ladder.csv`, column `c52`: V-BO-01, 60 cells, temperature form, 3e-3) and 0.994 under Phase 1's scale; 0.998 is the thinnest-cell ratio column. The verdict stands. | fixed (gate, plan §12.1) |
| S-1 | major | `test/newcore/gates/burnout_ladder.jl:416-417`, `burnout_ladder.md` | The item 1 column passes `length(inv) <= layers && all(largest ≤ 5 ε_m)`, so a run with no invocation passes: V-BO-08 at 3e-2 (`:interior` at 67.04 s) and V-BO-08's reference row without the correction read "pass" in the table; the gate document corrects them by hand (§1, C94). §7.3 asks for one invocation per layer in every run. | fixed (D96, `c36760a`; the column recomputed, `d3fc995`) |
| X-1 | major | gate item 8, `burnout_timing.md`, plan §5.7, §7.3 item 8, §7.4 | Item 8 is recorded, not decided, and on the recorded numbers the phase's stop rule is reached: V-BO-01 and V-BO-03 at 120 cells 2.11 and 2.09× at matched accuracy, the installed correction 3.73× on `ale` 100 and 2.38 and 4.62× on `ale_darcy` 25 and 100, all at a load of 17 to 21. The installed-correction rows are at equal tolerance (reltol 1e-6), not matched on the window RMS, so they are outside item 8's protocol either way. The V-BO-01 and V-BO-03 ratios at 60 and 120 cells have stood above 1.5× in every session since M4 (1.62 and 2.28×, 1.97 and 2.74× at M4; 1.85 and 2.11×, 1.88 and 2.09× at M6c), so a quiet machine is unlikely to bring them under the pass line alone; the lever named is `FBDF` (0.80× modeled and 0.83× wall on V-BO-01 on #4662's merge commit, Appendix C.7), which is MF's. R41's deck is the correction's `N²/τ` gain on a 100-cell single layer (1,620 steps against 1,107), which no release of `FBDF` removes. | fixed: the protocol (D101, `9729b62`); the quiet-machine run pending (§6); owner: the §7.4 decision (open) |
| G-3 | major | gate item 9 and §9, plan C93 | The numbers the recommended resolution rests on (spreads 0.0039 and 0.0079, σ 2.8e-5 and 2.2e-5 s at reltol 1e-7) were taken by "the same script, scratch": no committed table holds them, against D52, and the gate reads item 9 as failing on two of four configurations at the case's tolerance. | fixed (D100, `c174629`; the table, `b675f82`) |
| P-1 | major | pull request #14, plan §9.1 | Title "Plan Phase 2 (burnout) of the mesh framework" and body "Planning only: no package code changes", unchanged since M0, where §9.1 says the body is §11 and §13; 47 commits (from `0b14741`, M2) are unpushed, so the draft shows M1. | owner: the push and both pull requests (texts drafted, followup plan §4) |
| I-1 | minor | `src/Integration/setup.jl:395-399`, `src/Closures/correction.jl:203-209` | The band defaults to `CORRECTION_BAND` whenever the correction given is an `EquidistributionCorrection`, before `bind_correction` decides whether any layer moves. An explicit `EquidistributionCorrection()` on a deck whose layers cannot vanish, with deactivation on, runs the band at 0.1 with zero gain everywhere: C58's configuration, which changed V-CPL-02b (35 cells against 34, its energy-balance QoI failing). The installation rule itself cannot reach it; the public keyword can. | fixed (D91, `aa8fa06`) |
| I-2 | minor | `src/Integration/solve.jl:72-74, :403-405` | The `solve` docstring says every cell within `reltol` of its threshold goes with the crossing one; `deactivate!` takes the contiguous run from the top only, as its own docstring says. The callback comment says a cell under a correction can dip below its threshold and recover and that location fires at the first crossing, which removes it: recovery is invisible to the end-of-step net only. | fixed (D95, `6bb9f18`) |
| I-3 | minor | `src/Integration/solve.jl:422`, `src/Integration/checks.jl:18-23`, `src/Deactivation/events.jl:372` | After a merge the receiver holds the sum of `k` cells' component masses but is tested against its own per-cell `abstol`; `k` admissible negatives of `−abstol` (a near-depleted gas) sum to `−k·abstol` and end the solve `:Failure` from states the domain test accepted. Not observed on any gate run. | fixed: documented (D92, `4efeef8`) + test |
| X-2 | minor | `test/verification/cases/V-BO-01.jl:322-325`, D86 | The mass-ledger QoI is a trapezoid over the saved times, so V-BO-08's verdict depended on `save_every` (1.08e-3 at 0.5 s against 1e-3; 1.5e-4 at 0.1 s; the old core 3.2e-4 and 4.6e-5) and the case saves every 0.1 s to pass. The burst recorder of `burnout_burst.jl` already integrates the rate per accepted step. | owner (followup plan §7.1, open) |
| X-3 | minor | `src/States/derived.jl:365`, `src/Assembly/residual.jl:11`, C79 | In the enthalpy form `T_warm` carries the latest plain evaluation's temperatures, a step recorder's included and the previous `integrate(cp)` of the same `CoreProblem` included (nothing resets it), so the trajectory depends on the evaluation history at round-off; `residual!`'s docstring claims a pure function of `(u, t)` given `def`. The determinism test builds a fresh problem per solve and cannot see it. | fixed (D94, `705db0b`) + test; C79 stands |
| S-2 | minor | `test/newcore/gates/burnout_ladder.jl:206-214`, D73 | The probe is a continuous callback evaluated at the proposed step end, before the deactivation root is located, so on a located step the thinnest cell (item 10) and the scaled residual (item 6) are read on a state the located event then discards: conservative for item 10 (the cell is thinner there than anywhere on the trajectory), off-trajectory for item 6. The gate's 0.865 `ε_m h₀` is such a reading. | fixed (the gate, the ladder's preamble, `5111807`) |
| S-3 | minor | `test/newcore/gates/burnout_timing.jl:185-190` against `case_accuracy.jl:114-127` | `event_shift` is a drifted copy of `case_accuracy.jl`'s (per-cell event times through `unique`, no single-invocation rule, D77): on the new core one invocation is `N` events at one time, so the column is blank on every new-core row of `burnout_timing.md` while the script header and the table's preamble claim parity. | fixed (D97, `7587326`, `08de603`) |
| S-4 | minor | plan §7.3 items 3, 5, 6, 8 | Sub-criteria no script measures: item 5(A)'s smallest `h/h₀` per deck; item 5(B)'s largest activation, active nodes and V-BO-04's exactly-zero correction (only the equal step counts with and without it, indirectly); item 3's "no larger than the old core's closure" is read by eye, the script prints no verdict; item 8's 1.5 and 2 thresholds are not encoded; item 6's "before C52" exists only as `phase1_above` (D71). | fixed (D99: `dca69ae`, `706aa52`, `263262f`; plan §7.3 item 6) |
| G-4 | minor | gate document | "The same steps the review measured (1,620 on `ale_darcy` 100)": the review measured 1,623, 1,620 is this pass's; "1 of 39 cells" is one invocation of 39 cells; V-BO-03 "(enthalpy alike)" where its (c) column differs (0.591 and 0.237 against 0.727 and 0.488, both passing); "3,022.5 to 3,023.0 s" cited to the ladder, taken from the burst table; item 7 "every deactivating case" against the script's "with events" (16 of 20 rows); item 4 has no section of its own. | fixed (gate) |
| C-1 | minor | `src/Closures/correction.jl:133-151` | `vanishing` reads a self-yield reaction `A → αA + (1 − α) G` and a reversible chain as non-vanishing (conservative: no correction, where the layer does gasify), and counts a reaction that can never run (`A = 0`) as a consumer (the correction installed on a residue it never moves: harmless, the activation stays zero). | fixed (D93, `9902314`) + test |
| TS-1 | minor | `test/newcore/test_deactivation.jl` ("masked rows are identity rows"), `test_correction.jl` | No test evaluates the residual or the Jacobian on a post-merge state with the correction active at the node inward of the outermost live node; every masked-row state keeps its cells above `r_hi`. K-1 is the gap. | fixed (with K-1, `46059c2`) |
| R-1 | minor | `scripts/mesh_phase2_prototype/` | `baselines.jl` and `cases.jl` remain under §9.2's "until their measurements have gate scripts of their own": `vbo_old`, `vbo_new`, `vbo_ladder_old`, `vbo01_qois_old` and `df_cores` now have (`burnout_ladder.jl`, `burnout_identity.jl`, `burnout_timing.jl`), `vbo02_restart` is D80's measurement and a test, `integrators` and `tolweights` have none. | owner (open) |
| K-2 | owner | C73, C74 | Carried from M3: the unlimited slope beside a boundary cell (not positivity-preserving; a test deck drained a cell to `DtLessThanMin`) and the upwind side at an exact tie (taken from the dual's partials, so the colored and dense passes can differ on a tie). Neither showed on the ladder. | owner (open) |
| R-3 | owner | design §3.3, §3.5, §3.7, §6 | The errata of plan §7.2 are not applied: the design still multiplies the rates by the switch (line 139) and lists the "smooth reactivity switch-off" in Phase 2 (line 259). | owner (open) |
| D-1 | note | `src/Deactivation/events.jl:258-270` | `deactivate!`'s loop re-tests the receiver after a moving-mesh merge, so a receiver that is itself below the band and the threshold after absorbing the run merges on in the same invocation; §3.4 and the docstring describe one merge. Reachable only when the run's whole inventory is below `ε_m` of the receiver's reference. | fixed (D95, `6bb9f18`) |
| D-2 | note | `src/Deactivation/events.jl:359-375`, `src/Integration/solve.jl:423-424` | On a fixed mesh under Darcy–Fick a receiver with little pore space takes the dead cells' condensed mass into its unchanged volume, and `_porosity_out` then ends the solve `:Failure` rather than refusing the merge. D47's design; a hard failure where a rule could stand. | none |
| C-2 | note | `src/Closures/correction.jl:230-236` | `activation`'s docstring says the smoothstep's argument is clamped to `[0, 1]`; the code clamps above only. Below `r_hi` the argument is negative only through the smooth floor's excess (of order 1e-10 at most), where `x²(3 − 2x)` is 1e-20. | none |
| P-2 | note | plan §2.2, §3.9 | Still list `switch_width` and "Switch width 1.5 `ε_m`" as Phase 2 content; §10 and §12 record the deletion. The plan's convention. | none |
| R-2 | note | `docs/refs/phase2_commit_queue/` (local, gitignored) | The eight patches are applied; `TREES.txt` is the only record of G-1's mapping. | done (the local directory deleted after G-1's mapping) |
| S-5 | note | gate scripts | `burnout_timing.jl:66-76` includes `scripts/mesh_phase0/baselines/jacobian_common.jl` by relative path; `burnout_identity.jl:369-407` scrapes `examples/07_douglas_fir_cone.jl` between two text anchors and evaluates it; `burnout_ladder.jl:620` hard-codes "activation 0.3/0.1"; `_g` is defined in all seven scripts with three NaN behaviours; under location the crossing cell sits exactly at `ε_m`, so `e.fraction < ε_m` (`V-BO-01.jl:456`) counts it in the band ("0 below" on V-BO-03 and V-BO-04's rows). | none |

## 2. Findings by area

### 2.1 The correction (C)

Verified correct: the node correction is §3.2's with the `N²/τ` gain per layer, zero at the pinned, end and Lagrangian nodes and beside a dead cell (`_correction!`); the "or" activation and its `C¹` joins; the monitor normalized by the initial mean (C72), so `β = 0` gives `M̃ = 1` exactly; the exact zeros with zero partials above `r_hi` (D66), which `iszero` of a dual reads, so the kernel's skip is seeding-independent and the colored Jacobian equals the dense one; `bind_correction`'s line-mesh guard and layer bounds; the `vanishing` fixpoint on chains, char formers, competing paths and cycles as the tests state; the balance rows' sweep and the one-pass consistent velocities (`w − c` the solid velocity on planar, cylindrical and spherical meshes); the velocity rows' widened reach in the pattern (positions `k − 1 … k + 2`, the masses only with `β ≠ 0`), two colors below the design's count (D67), with the test's formula handling the prescribed and equidistribution cases apart. Findings: C-1, C-2; I-1 is the binding's.

### 2.2 Relative fluxes (K)

Verified correct: the kernel reads `ds.c` alone and the solid velocity never appears; the face-conservative mass and enthalpy fluxes; the `NaN` propagation of C47 (`T * NaN` in `_tabulated` and `enthalpy_slope`, tested through a tabulated law). Finding K-1 is the one defect of this review: `_stencil` has no notion of a dead cell. Every other kernel reaches the live cells through `live_face` or `ds.active`; this one was written for a prescribed velocity over a mesh that never deactivates (Phase 1's refusal), and the correction brought it into the deactivating decks without that guard. K-2 carries C73 and C74.

### 2.3 Deactivation and events (D)

Verified correct, beyond the module's own tests: the node indices of `_nodes`, `_enclose!` and `_mask!` (the shared node stays live and holds the merged volume's outer face, the dead cells' far nodes collapse onto it and are masked); the merge's bookkeeping over several dead cells (masses, the vented gas on a fixed mesh, sensible enthalpy, the history pair, created volume, the regularization `E − E_m − H_vented`, all on the first event so the sums are the merge's, D53); the run reaching cell 1 and the survivor below `ε_m` lost in the same invocation; a receiver with zero reference never divided by and never deactivated; `threshold_margin`'s sign exactly the threshold test's (D76); two crossings in one step located at their own roots; the initial-time root left to the safety net; the slack's effect on the final `:interior` test benign (a cell inward of the run is above `ε_band (1 + slack)`); the mutable state reset in `integrate`. Findings: D-1, D-2, I-3.

### 2.4 Integration, setup, policies, the linear solve (I)

Verified correct: the installation rule (`use_ale`, deactivation, some cell vanishes) with the composition interfaces Lagrangian by default (C85) and `_require_layers` on an explicit list; the refusals (a correction on a fixed mesh, beside `w_override`, a prescribed velocity on a reacting or deactivating deck); the band's range check; the constraint scale with the sweep (C52) and the round-off check of `_consistent_motion!`; the heuristic restart (`auto_dt_reset!`, then the proposal) and the policy's three restart kinds; `StructuredKLU`'s union (grows on a nonzero outside it, errors outside the prototype, the diagonal always in), with the tests' one analysis per growth; the located callback's constructor (`affect!` positional `nothing`, `affect_neg!` the event, `RightRootFind`), the safety net after it, and the terminated outcomes saved. The trial clamp of C36 is handed to the shared Newton as a keyword whose default keeps the old core's clamp and warnings, and the fingerprint, `test_surface_bc_guards.jl` and the old core's `Pkg.test()` were run at M1; nothing in the shared files changed after. Findings: I-1, I-2, I-3, X-3.

### 2.5 Tests (TS)

The seventeen tests of §6.1 are present: items 1 to 5, 11, 14 and 16 in `test_correction.jl`, 6 to 10 in `test_deactivation.jl`, 12 in `test_states.jl` and `test_relative_flux.jl`, 13 in `test_boundary_conditions.jl`, 15 in `test_allocations.jl`, item 11's color count in `test_pattern.jl`; item 17 (the fallback) was not built. The located-event tests assert what M5 records (the shift falling with the tolerance, the same cells as the end-of-step event, the threshold-exact crossing, the safety net's crossing, two roots in one step, the admissible located state). `runtests.jl` runs `test/newcore` in the temperature form only (D44); the correction's enthalpy-form coverage is `test_correction.jl`'s explicit loops and the suite. Finding TS-1.

### 2.6 The cases (X)

V-BO-01's QoIs read the main peak before the last event (C86), the one-invocation QoI with the merged-fraction clause, and the clamped-balance count on accepted states (D72); V-BO-03's transmitted-flux predicate on the series; V-BO-04's plain law at 1e-6 and the located crossing at 1e-6 (§14.6 rows 1 and 2); V-BO-05's `:interior` stop at the 240-cell time with a ladder-derived tolerance (C45); V-BO-06's grading held across resolutions (D89); V-BO-07's references from the old core at 240 cells (D87) and its invocation asserted without the bound (C91); V-BO-08's zero correction on the substrate and its 0.1 s saves (D86). Every case names its re-pin's old value and reason (R37). Findings: X-2, X-3; X-1 is the gate's.

### 2.7 The gate scripts and tables (S, G)

Every number of the gate document's items 1, 2, 3, 5, 7, 8 and 9 was found in its table (the item 1 counts and fractions, the QoI triples, the burst measures and windows, the ledger closures, the digests and step counts, the sixteen shift and gain pairs, the matched and equal-tolerance ratios, the profiler shares, the gradient spreads); the scripts' criteria match §7.3 where §7.3 states one (the burst's 1.5× per doubling and 1.5× against the old core, item 7's 0.05 and threefold gain, item 10's `0.5 ε_m h₀`, item 9's 0.1 and the constant count). Findings: G-1 to G-4, S-1 to S-5.

### 2.8 The documents and the rules (P, R)

The plan's §12 is the status of record and records every departure and conflict with its evidence; §14 records the decision review's implementation; the gate document reads every item on the final defaults and names its open decisions. No AI attribution anywhere; every commit signed; one theme per commit; no dates or history in code or docstrings (the one grep hit, `EQUIDISTRIBUTION_TIME`'s docstring, states a justification, not a history); the guides untouched; the formatter, Aqua and JET green. Findings: P-1, P-2, R-1 to R-3.

## 3. What remains

In the order to do it.

1. **The review's fixes** (one theme per commit): K-1 with its test (TS-1), the class B digests re-recorded and `burnout_identity.jl` and the two V-CPL rows of `case_accuracy.jl` re-measured; I-1; I-2; I-3; C-1; X-3's reset and docstring; S-1's criterion and the regenerated column; S-3; G-1's mapping in the gate and plan §12.1 and `stamp.jl`'s tree hash; G-2 and G-4; S-2's sentence; D-1's sentence; then R-2.
2. **The measurements still owed to the gate:** `burnout_timing.jl` on a quiet machine with the installed-correction rows matched on the window RMS (X-1), `burnout_gradients.jl` with a reltol 1e-7 configuration (G-3), the sub-criteria of S-4 either measured or restated. None changes code.
3. **The owner's decisions, which this review does not take:** item 8 under §7.4, on the quiet-machine table (X-1); C91 (V-BO-07's merged fraction), C92 (the old core's ledger on V-BO-06), C93 (item 9's noise), C94 and V-BO-06 at 3e-2 (the owner exceptions of item 1); C73 and C74 (K-2); the design errata (R-3); the prototype directory (R-1); X-2's ledger.
4. **Design work the gate names, for the cutover, not Phase 2:** the layer time `τ`. It is one constant, 720 s, tuned on V-BO-01's 5.9 mm layer; the gain `N²/τ` makes a 100-cell single layer stiff (R41: 1,620 against 1,107 steps, 3.7 to 4.6× the old core on `ale` 100 and `ale_darcy`) and a 1 mm layer under-synchronized (C94: `τ/π²` = 73 s against an activated life of 14 s). The plan's own heuristic (`T_sync` at most a tenth of the layer's remaining life after full activation, §3.2) is not a rule in the code. R27's lever, a `τ` per layer from the deck, is the one item that addresses both; whether it is set from the layer's thickness and a regression-rate scale at setup or from the observed thinning at activation is the design question.
5. **MF** on the upstream release, as §14.3 lists it, re-running items 5, 8 and 9 and the suite, with item 12.
6. **Process:** push the branch and update pull request #14 (P-1) before the owner's review of the gate.

## 4. What this review ran

On `8166e77`, whose source and test tree equals the measured tree (`git diff 78ce699..HEAD` outside `measurements/` and the two documents is empty; `78ce699`'s tree is `543bbb8`, the gate's): `julia --project=scripts scripts/format.jl --check` (all files formatted), `test/quality/aqua.jl` (9 of 9), `test/quality/jet.jl` (no Pyrolysis-owned reports). No solve and no `Pkg.test()`: the machine carried fourteen `foamRun` processes at a load of 14.5 for the whole session, and the three-way `Pkg.test()` of M7 ran on this code. The two independent reads were given the files and the contracts, not the coordinator's conclusions; their findings were re-derived from the code and the tables before entering §1 (K-1 traced through `_stencil`'s `far` table and `_enclose!`; G-2 against `burnout_ladder.csv`'s `c52` column; S-1 against `burnout_ladder.jl:416-417`; S-3 against `burnout_timing.md`'s new-core rows; S-2 against the probe's `ContinuousCallback`).

## 5. What the followup implemented, with gates

On branch `phase2-review` (`docs/Mesh_Framework_Phase2_Review_Followup_Plan.md`, the owner's answers in its §7), 26 signed commits, one theme each, in the followup plan's §5 order but for three departures: S-2's sentence in the ladder's preamble is a commit of its own (`5111807`), the event-shift guard that the first event-timing run showed missing is another (`08de603`), and the order follows the checks' completion: item 8's thresholds (`263262f`) before items 5 and 3 (`dca69ae`, `706aa52`), item 9's configurations (`c174629`) before the installed decks' protocol (`9729b62`), and the event-timing table (`f073c3b`) after the ledger and gradient tables. Plan §12.2 and §12.3 record D90 to D101 and C95 to C97.

### 5.1 Code

| Finding | Commit | Change | Gate |
|---|---|---|---|
| K-1, TS-1 | `46059c2` | `_stencil` treats a dead cell beyond the upwind cell as absent (D90) | the new testset fails on `8166e77` (three of four perturbations and the dead-column test, each form) and passes; `test_correction.jl` 719, `test_relative_flux.jl` 410, `test_deactivation.jl` 536, containment 2,046, pattern 319, allocations 35 (zero bytes with the correction active); the twelve V-BO digests equal to the parent's; V-CPL-02/02b bit-identical in both forms (C97); `ale_darcy` 25 and 100 change (585 to 583, 1,620 to 1,623 steps) |
| I-1 | `aa8fa06` | the band's default from the bound correction (D91) | `test_correction.jl` 721, `test_deactivation.jl` 536 |
| I-3 | `4efeef8` | the per-cell post-event test kept and documented (D92) | `test_deactivation.jl` 548, `test_policies.jl` 202 |
| C-1 | `9902314` | `vanishing` ignores a self-yield (D93) | `test_correction.jl` 724; the identity classes unchanged |
| X-3 | `705db0b` | `integrate` resets `T_warm` (D94) | the new test fails on `8166e77` and passes; `test_policies.jl` 209, `test_solve.jl` 75; the twelve V-BO digests (six in the enthalpy form) equal to the parent's |
| I-2, D-1 | `6bb9f18` | docstrings and a comment (D95) | on this last source commit: the format check clean, Aqua 9 of 9, JET with no Pyrolysis-owned reports, the old-core fingerprint equal to 17bc477's digests |

### 5.2 Gate scripts

| Finding | Commit | Change |
|---|---|---|
| G-1 | `25dac70` | `stamp.jl` prints the tree beside the commit (D98) |
| S-3 | `7587326`, `08de603` | `event_shift.jl`, shared by `case_accuracy.jl` and `burnout_timing.jl`; the second commit restores the guard for a solution without events, which the first dropped (the MMS rows) (D97) |
| S-1 | `c36760a`, `d3fc995` | `item1_verdict` (exactly one invocation per layer) and `--item1`; the column recomputed: four rows fail (D96) |
| S-2 | `5111807` | the ladder's preamble says what the probe reads on a located step |
| S-4 | `dca69ae`, `706aa52`, `263262f` | item 5's thinnest cell, activation and V-BO-04; item 3's verdict; item 8's thresholds (D99) |
| X-1 | `9729b62` | the installed decks under item 8's protocol, the cost split at the first activation, the equal-gain row (D101) |
| G-3 | `c174629` | item 9's configurations at reltol 1e-7 on the finite-volume core (D100) |

### 5.3 Tables

Measured in a detached worktree of `6bb9f18` (tree `44e66df`), `case_accuracy.jl` on `08de603` (the same worktree moved to it; the two trees differ in one gate script), one script at a time after one precompile, at loads of 17 to 19 (fourteen OpenFOAM processes ran throughout), and stamped with the commit and its tree.

| Table | Commit | Result |
|---|---|---|
| `burnout_identity.md`, `burnout_identity_reference.csv` | `946c1fb` | compare mode on M6c's reference: every deck equal, V-CPL-02/02b included (C97), V-BO-04 without a reference; `--record` adds V-BO-04, and compare mode reads all nine equal |
| `case_accuracy.{csv,md}` | `f073c3b` | on `08de603`: all 100 rows equal M6c's but the wall time (the first run, on `6bb9f18`, lost the eight MMS rows to the event-shift guard, D97); item 7 passes on every case with an event |
| `energy_ledger.md` | `b4a0b16` | every number M6c's; item 3 passes in both forms |
| `burnout_gradients.{csv,md}` | `b675f82` | the case's tolerance rows M6c's to every digit; at 1e-7 every configuration passes |
| `burnout_ladder.{csv,md}` | `d3fc995` | recomputed from the csv on `c36760a`, no solve: four rows fail item 1 |
| `burnout_timing.{csv,md}` | none yet | waits for a quiet machine (§6) |
| `burnout_burst.{csv,md}` | none | no fix reaches the V-BO decks (the twelve digests); M6c's, stamped `5f08134` (C95) |

### 5.4 Documents

The gate document (G-1's mapping, G-2, G-4, S-1, S-2, item 4's section, the re-measured numbers, C96 and C97), plan §7.3 item 6, §9.3, §12.1 to §12.4 and §14.6 row 6, and this review. R-2: the local `docs/refs/phase2_commit_queue/` deleted once G-1's mapping was recorded.

## 6. What was not implemented, and why

- **The owner's decisions** (followup plan §7.1): item 8 under §7.4, C91, C92, C93, C94, V-BO-06 at 3e-2, C73 and C74 (K-2), the design errata (R-3), the prototype directory (R-1), the ledger QoI (X-2), the per-layer `τ` and MF. Nothing in the followup presumes an answer.
- **S-2's probe at the located state:** not taken (followup plan §2.4): it needs an observer hook in `integrate` and a ladder re-run, to tighten a reading that is already conservative.
- **I-1's refusal and I-3's summed tolerance:** the alternatives of followup plan §1.2 and §1.3, not taken for the reasons given there (D91, D92).
- **The class B identity reference for D90:** not re-recorded, since D90 leaves V-CPL-02/02b bit-identical (C97); the re-record adds V-BO-04 alone.
- **Item 8's quiet-machine run** (followup plan §3.2): the machine carried fourteen OpenFOAM processes at a load of 15 to 19 for the whole followup, so `burnout_timing.jl`'s new protocol (D101) has not been run beyond the scratch check of `9729b62`. By question 10 the followup goes to review with item 8 pending; its table and the gate's item 8 text follow as one table commit and one gate commit on `phase2-review` before the merge.
- **P-1:** the push, the retitle and the bodies of #14 and the followup's pull request are the owner's (followup plan §4); their texts are drafted.

## 7. The three-way `Pkg.test()`

In three detached worktrees of `6bb9f18`, the last commit that changes source or tests (the later commits change gate scripts, which `runtests.jl` skips, tables and documents), started one at a time after the previous one's precompile and run side by side with the followup's measurements at loads of 16 to 19:

| run | tests | wall |
|---|---|---|
| old core | 13,668/13,668 | 109 min (106 min of tests) |
| `PYROLYSIS_CORE=new`, temperature form | 10,074/10,074 | 82 min (79 min of tests) |
| `PYROLYSIS_CORE=new`, enthalpy form (the verification suite alone, D44) | 543/543 | 23 min (19 min of tests) |

Against M7 (13,598, 10,004 and 543 at `5f08134`), the old-core and temperature-form runs grow by the followup's 70 tests, which both run (`test_correction.jl` 673 to 724, `test_deactivation.jl` 536 to 548, `test_policies.jl` 202 to 209); the enthalpy-form run is the verification suite alone and stays at 543.

## 8. The owner's answers

| Question (followup plan §7.2) | Answer | Where |
|---|---|---|
| 1. I-1: the band from the bound correction | as recommended | D91 |
| 2. I-3: the per-cell test documented | as recommended | D92 |
| 3. X-3: the reset to `T0` | as recommended | D94 |
| 4. S-1: recomputed from the csv | as recommended | D96 |
| 5. S-2: a sentence | as recommended | the gate's §6 and 10, the ladder's preamble |
| 6. S-3: one shared function | as recommended | D97 |
| 7. S-4: as its table | as the table | D99 |
| 8. G-3: all four configurations at 1e-7 on the finite-volume core | as recommended | D100 |
| 9. X-1: the split at the first activation and the equal-gain row | as recommended | D101 |
| 10. The pull request and item 8 | as recommended: the machine was not quiet, so the pull request carries item 8 pending | §6 |

The phase's decisions of followup plan §7.1 stay open: item 8 under §7.4, C91, C92, C93, C94, V-BO-06 at 3e-2, C73 and C74, the design errata (R-3), the prototype directory (R-1), the ledger QoI (X-2), the per-layer `τ` and MF; the followup adds C96 (V-BO-04's correction at round-off).
