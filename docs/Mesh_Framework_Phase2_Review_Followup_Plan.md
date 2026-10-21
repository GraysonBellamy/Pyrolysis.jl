# Mesh framework: the Phase 2 review's followups

**Scope.** The plan for the followups of `docs/Mesh_Framework_Phase2_Review.md` (the review of `mesh-framework-phase2` at `8166e77`): its findings with status `followup` (§1 and §3 items 1, 2 and 6 of the review), grouped as the review's §3 groups them, each with the change, the files, the test that proves it, the gate it must keep, the record it writes and an estimate; then the commit list, the re-measurements and their machine conditions, the owner's inputs and the followup's closing gate. Nothing here re-reviews the phase, and nothing starts MF or the design of a per-layer `τ`.

**Branch.** `phase2-review`, off `mesh-framework-phase2` at `8166e77`; its pull request targets `mesh-framework-phase2`, as Phase 1's `phase1-review` targeted `mesh-framework-phase1` (#13). Nothing goes to `main`, and #14 stays against `mesh-framework-phase1`. Records continue at D90 and C95 in plan §12.2 and §12.3. No old-core file is touched by any item (`src/BoundaryConditions`, `src/Properties`, `src/Materials`, `src/Solver`, `src/Discretization`, `src/Jacobian`, `src/Physics`); the code fixes are in `src/Kernels`, `src/Closures`, `src/Integration`, `src/Assembly` and `src/Deactivation`.

**What was verified for this plan.** `5f08134` is not an object of the repository; `78ce699` ("Time the correction where the installation rule installs it") has tree `543bbb8`, and `git diff 78ce699..HEAD` touches the two documents and the tables under `measurements/` only, so every table describes `8166e77`'s code. The `file:line` references of the review's §1 hold on `8166e77` except S-3's: `event_shift` is `burnout_timing.jl:205-210` (the review cites 185-190) against `case_accuracy.jl:114-129`. I-3's `events.jl:372` is the line that sums the run's masses into the receiver (written to the state at `:389`). Two facts the plan rests on: the receiver of every V-BO deck lies in a zero-gain layer (the Kaowool under V-BO-01, 03, 06, 07 and 08; V-BO-04's survivor is cell 1, beside the pinned node), and in the enthalpy form `setup` leaves `T_warm` equal to `T0` exactly (`_invert_primal` returns its starting point when the residual there is zero, which it is for the energy `setup` computed from `T0`).

**Contents.** 1 Code fixes · 2 Gate scripts and tables · 3 Re-measurements · 4 Process · 5 Commits · 6 Records · 7 The owner's inputs · 8 Closing gate · 9 Estimate

---

## 1. Code fixes

### 1.1 K-1 and TS-1: the relative-flux stencil reads a dead cell

- **Change.** `_stencil` (`src/Kernels/relative_flux.jl:135-145`) takes the live-cell mask and treats `X`, the cell beyond the upwind cell, as absent when `!active[X]`: `X = kernel.far[side, f]`, then `X = X != 0 && !active[X] ? 0 : X`. `contribute!` passes `ds.active`. On the face at node `r` inward of the outermost live node after a merge, the face state is then the upwind cell's value with the downwind gradient alone, exactly the stencil of the outermost cell of a mesh that never had the dead cells. For the concentrations the van Albada slope already approached it (the upwind gradient against a zeroed cell is large); for the temperature it did not, and the two forms differed: the frozen `T_m` gives an upwind gradient near zero in the temperature form, the inversion of a zero inventory (`T_REF`) a large one in the enthalpy form. `D` needs no guard: the only face with a dead downwind cell is the one at the outermost live node, and `_correction!` zeroes `c` there, so the kernel skips it. The kernel's docstring says "where no live cell lies beyond"; `derived.jl:145-146` ("the kernels skip every face and patch beyond them") becomes true and stays as written.
- **Files.** `src/Kernels/relative_flux.jl`; `test/newcore/test_correction.jl`.
- **Test.** A testset in `test_correction.jl`, both energy forms, "the relative fluxes after a merge read no dead cell": `_ablator(6)` with `_ablation` (moving mesh, deactivation, `EquidistributionCorrection(τ = 2.0)`), its cells at `h/h₀ = (1, 1, 1, 0.2, 0.25, 0.02)` by `_with_thicknesses`, cell 6's condensed mass set to 0.5 % of its reference; `deactivate!` by hand merges cell 6 into cell 5 (volume 0.27 `h₀`, below `r_hi`). Asserted: the precondition, `c[5] > 0` (node 5 active between two live cells, so the face at node 5 takes side 2 with `X = 6`); the residual's live rows are `==` the unperturbed ones after each of cell 6's rows (its masses, its energy) and node 7's position row is perturbed in turn; the dense Jacobian's live rows are zero in every dead column; `_contained` at that state (the pattern contains the dense Jacobian and the colored one equals it, the same helper as `jacobian/test_containment.jl`'s). The perturbation assertions fail on `8166e77` (the face reads `ξ(X)`, `T(X)` and `d_XU` through node 7) and pass after the fix. `reset_cells!` restores the mask.
- **Gate.** `test_correction.jl`, `test_relative_flux.jl`, `test_deactivation.jl`, `jacobian/test_containment.jl`, `jacobian/test_pattern.jl`, `jacobian/test_allocations.jl` (zero bytes per residual with the correction active: the mask read allocates nothing). The containment proof: digests of `(T, m, z, t)` of V-BO-01, 03, 04, 06, 07 and 08 at their CI resolutions as the suite solves them, both forms, on the parent and on the fix, twelve pairs equal (a scratch script, as M3's `identity.jl`; the digests go into D90). The class B decks change: `burnout_identity.jl` in compare mode reads class A equal and V-CPL-02 and V-CPL-02b different, then `--record` re-records them (§14.6 row 6, again); `NEW_CORE_CASES` passes V-CPL-02/02b in both forms with the QoI values recorded as §14.6 row 3 records them; `case_accuracy.jl` is re-measured for item 7 (§3.1). Not re-run, because nothing they solve reaches the path: `burnout_ladder.jl`, `burnout_burst.jl` and `burnout_gradients.jl`'s 1e-6 rows (V-BO decks only, receivers in zero-gain layers, proved by the digests), `energy_ledger.jl` (V-BO-04, whose survivor is cell 1, and V-CPL-06's char-layer column, solved without deactivation), and the V-BO rows of `burnout_timing.jl`. The timing table's installed-correction rows (`ale_darcy` and V-CPL-02/02b after their first event) do reach it, and item 8's run (§3.2) re-measures them after the fix.
- **Record.** D90 (old: `X = far[side, f]` whatever its state; new: absent when dead; reason: K-1; the twelve V-BO digests; the class B digests' old and new values; V-CPL-02/02b's QoI values). The face at the receiver now takes C73's slope (the downwind gradient without a limiter), which the owner's answer to C73 then covers too (§7).
- **Estimate.** 1.5 days (the fix an hour; the test half a day; the digests, the re-record and `NEW_CORE_CASES`' two cases half a day; `case_accuracy.jl` is §3.1's).

### 1.2 I-1: the band default ignores whether the correction moves anything

- **Change.** Take the band's default from the bound correction: `CORRECTION_BAND` (or the threshold if larger) when the closure's correction is a `LayerCorrection` with some nonzero gain, the threshold otherwise. The `something(band, …)` of `setup.jl:395-399` moves below the closure's construction (`CellDeactivation` is built after it, so this is a local reorder), with a helper `_moves(closure)`; the docstrings of `setup` (`:341-347`) and `CORRECTION_BAND` (`:141-149`) say "when the bound correction moves some layer". The refusal (an equidistribution correction that moves nothing while cells deactivate) is not taken: it would break item 5(A)'s explicit-correction check on any class A deck with deactivation on, and a correction that moves nothing synchronizes nothing, which is all the band is for (C58's own reading of §3.4).
- **What C58's configuration does afterwards.** An explicit `EquidistributionCorrection()` on a deck whose layers cannot vanish (or a one-cell vanishing layer, which has no interior node and so no gain), with deactivation on, runs the band at the threshold: its deactivation is `NoCorrection`'s (the contiguous block below `ε_m` from the exposed side), the relative-flux kernel stays installed with `c ≡ 0` and skips every face, and the trajectory agrees with `NoCorrection`'s to round-off through the wider pattern (C56's 1e-10), not bit for bit. The V-CPL-02b change of C58 cannot recur.
- **Files.** `src/Integration/setup.jl`; `test/newcore/test_correction.jl`.
- **Test.** In the setup testset ("setup binds the correction to the layers…"): `_inert(6)` with `use_ale`, `deactivation` and `EquidistributionCorrection()` has band `1e-2`; with `layers = :all`, `CORRECTION_BAND`; the existing assertions (the two-layer deck at `CORRECTION_BAND`, `NoCorrection` at the threshold, `threshold = 0.2` at 0.2) unchanged.
- **Gate.** `test_correction.jl`, `test_deactivation.jl`. No table moves: no verification case passes an explicit correction, the ladder passes one only on V-BO decks (their PMMA layers move), and the identity script's explicit-correction decks (V-CPL-03, the Douglas fir deck) run without deactivation.
- **Record.** D91 (old: D54's default keyed on the correction's type; new: keyed on the bound gain; reason: I-1).
- **Estimate.** Half a day.

### 1.3 I-3: a merge's receiver against its own tolerance

- **Change.** Keep the per-cell test and document it. The post-event test (`solve.jl:422` → `inadmissible`, `checks.jl:13-30`) holds the receiver's merged masses to the receiver's own `abstol`, which is also what the domain test holds every later proposal to. Summing the merged cells' tolerances in the post-event test alone would pass an event whose next proposal the domain test must reject, so a receiver at `−k·abstol` would end in `DtLessThanMin` rather than in the named `:Failure`; summing them in the integrator's tolerance vector as well would change every trajectory after an event. The docstrings say so: `inadmissible`'s (each live cell, a merge's receiver included, against its own tolerances), the callback's comment in `solve.jl`, and `solve`'s sentence on `:Failure`.
- **Files.** `src/Integration/checks.jl`, `src/Integration/solve.jl` (docstrings and a comment); `test/newcore/test_deactivation.jl`.
- **Test.** A moving deck whose two outermost cells (below `ε_m`) and receiver each hold a gas mass of `−0.6` times its tolerance (equal cells, equal tolerances): after `deactivate!` the receiver's merged `−1.8` makes `inadmissible(u, layout, cp.abstol)` true, and at `−0.3` each the merged `−0.9` leaves it false; `DomainCheck(cp).abstol === cp.abstol` pins that the domain test applies the same tolerances to the proposals after the event.
- **Gate.** `test_deactivation.jl`, `test_policies.jl` (the domain test's testset). No code path changes, so nothing else.
- **Record.** D92 (the post-event test is the domain test's, by design; reason: I-3).
- **Estimate.** Half a day.

### 1.4 C-1: `vanishing` and a self-yield

- **Change.** In `vanishing` (`src/Closures/correction.jl:133-151`) a product equal to the reactant is not counted among the reaction's products, and a reaction whose yield of its own reactant is one consumes nothing (without that rule `A → A` would read as vanishing). The reversible chain and the reaction that can never run (`A = 0`) are recorded as they are, conservative and harmless. The constructors accept a self-product (no validation forbids it; checked on the first test).
- **Files.** `src/Closures/correction.jl`; `test/newcore/test_correction.jl`.
- **Test.** In the vanishing testset ("a solve gets the correction where…"): `1 => (1, 3)` with yields `(0.4, 0.6)` gives `(true, false, false)`; `1 => (1, 2)` with `(0.4, 0.6)` gives `(false, false, false)` (the product `b` is a residue); `1 => 1` alone gives `(false, false, false)`.
- **Gate.** `test_correction.jl`; `burnout_identity.jl`'s classes unchanged on every deck (no current deck has a self-yield); the suite at the closing gate.
- **Record.** D93 (old: a self-yield read as non-vanishing; new: as the rest of its products decide; reason: C-1).
- **Estimate.** Half a day.

### 1.5 X-3: the warm start across integrations

- **Change.** `integrate` resets the shared warm start, `copyto!(cp.ws.T_warm, cp.T0)`, beside `reset_cells!`, so every `integrate(cp)` starts its inversions from the set-up temperatures. The warm start inside a solve stays (C79, the owner's). `residual!`'s docstring (`src/Assembly/residual.jl:11`) replaces "a pure function of `(u, t)` given `def`" by the caveat: in the inventory form the caloric inversion starts from the workspace's latest plain evaluation, so the result depends on the evaluation history at the inversion's round-off; in the temperature form it is a pure function of `(u, t)`. No record labels in the docstring.
- **Files.** `src/Integration/solve.jl`, `src/Assembly/residual.jl`; `test/newcore/test_policies.jl`.
- **Test.** Beside "repeated solves are bit-identical…": one set-up problem in the enthalpy form on a deck with events, `integrate(cp)`, then `residual!` at a perturbed state (which moves `T_warm`), then `integrate(cp)` again: the two solutions equal, and equal to a fresh setup's. It fails on `8166e77`.
- **Gate.** `test_policies.jl`, `test_solve.jl`. No table moves: the temperature form never reads `T_warm`; the first `integrate` after `setup` starts where it did, since `setup` leaves `T_warm == T0` (the six enthalpy-form digests of §1.1 confirm it); no gate script integrates one set-up problem twice in the enthalpy form (the ladder sets up per solve; the timing and gradient scripts run the temperature form).
- **Record.** D94 (old: `T_warm` carried from the previous `integrate`; new: reset to `T0`; reason: X-3; C79 stands).
- **Estimate.** Half a day.

### 1.6 I-2 and D-1: docstrings that describe more than the code does

- **Change.** `solve`'s docstring (`solve.jl:72-74`): at a located root the threshold and the band are widened by `reltol` (`deactivate!`'s `slack`), and the cells that go are the contiguous run from the exposed side that the widened tests take, which holds the crossing cell; a cell within `reltol` of its threshold inward of a cell above the band stays. The callback's comment (`solve.jl:403-405`): the margin is tested for downward crossings, and location takes a cell at its first crossing, so a cell that would have recovered by the step's end goes there; only the end-of-step net, which sees the step's end alone, misses a dip and recovery inside one step. `deactivate!`'s docstring (`events.jl:227-248`): after a merge the receiver, holding the run's inventories against its own reference, is tested in turn and, when below the band with a cell below `ε_m`, merges on in the same invocation.
- **Files.** `src/Integration/solve.jl`, `src/Deactivation/events.jl`.
- **Test.** None (text only); the formatter and JET.
- **Record.** D95 (the receiver's re-test is the rule, against §3.4's "one merge"; the located event's run is the contiguous one).
- **Estimate.** Two hours.

## 2. Gate scripts and tables

The script changes land before any re-measurement (§5), so each table is measured once, on the followup's final code.

### 2.1 G-1: a stamp that survives re-signing

- **Change.** `test/newcore/gates/stamp.jl` reads `git rev-parse --short HEAD^{tree}` beside the commit: `"commit 78ce699 (tree 543bbb8)"`, with the uncommitted-changes suffix as now. The mapping goes into the gate's scope ("the tables are stamped `5f08134`, the pre-signing hash of `78ce699`, whose tree `543bbb8` they measured") and into plan §12.1's M6c and M7 rows. Plan §9.3 gains the rule: a gate table is measured in a detached worktree of a signed commit, never of a commit still to be signed or rebased.
- **Files.** `test/newcore/gates/stamp.jl`; the gate document; the plan.
- **Test.** The next table written carries both hashes (§3.1's first run).
- **Record.** C95 (old: seven tables and plan §12.1 stamped with a hash that does not exist; new: the mapping recorded, the tree printed; resolution D98).
- **Estimate.** Two hours.

### 2.2 S-1: item 1's column passes a run without an invocation

- **Change.** One function, `item1_verdict(invocations, largest, ε_m; layers = 1)`, returns pass when `invocations == layers` and the largest merged fraction is at most `5 ε_m` (§7.3: one invocation per layer in every run); `measure` uses it (`burnout_ladder.jl:416-417`). The column is recomputed without solving: the criterion reads only `invocations`, `largest_fraction` and `eps_m`, which `burnout_ladder.csv` holds per row, so an `--item1` mode reads the csv, rewrites its `item1` column and the item 1 cell of the matching rows of `burnout_ladder.md` (matched on case, cells, form, `ε_m`, `τ`, correction, `β`), and adds to the stamp line "item 1 recomputed from this table's csv at commit … (tree …)". A full re-run would cost one to two hours at load and reproduce the same numbers, since no fix reaches the V-BO decks (§1.1).
- **What changes.** Four rows go from pass to fail: V-BO-08 at `3e-2` with the correction (`:interior` at 67.04 s) and V-BO-08's reference row without it (`:interior` at 68.86 s), in both forms. The first is C94's row, which the gate already reads as failing. The second is not a matrix run (§7.3 names V-BO-08 at its three thresholds with the correction); the gate says so beside it instead of the criterion exempting it.
- **Files.** `test/newcore/gates/burnout_ladder.jl`; `measurements/burnout_ladder.{csv,md}`.
- **Test.** The recompute's diff touches the `item1` column of exactly those four csv rows and their md cells, plus the stamp line.
- **Record.** D96.
- **Estimate.** Half a day.

### 2.3 S-3: one event-shift measure

- **Change.** A shared file `test/newcore/gates/event_shift.jl`, included as `stamp.jl` is, with `invocation_times(sol)` (the new core's invocations, each at its first event's time; the old core's park times) and `event_shift(sol, ref) -> (shift, interval)` under D77's rule (the median interval between the reference's invocations, or between its saved times for a single invocation). `case_accuracy.jl` and `burnout_timing.jl` drop their own copies, and the timing table's column divides the two. On the new core the function is `case_accuracy.jl`'s, so item 7's numbers do not move; the timing table's column fills on its new-core rows at item 8's run.
- **Files.** `test/newcore/gates/event_shift.jl` (new), `case_accuracy.jl`, `burnout_timing.jl`.
- **Test.** `case_accuracy.jl`'s re-run (§3.1) reproduces every shift of the rows K-1 does not touch.
- **Record.** D97.
- **Estimate.** Half a day.

### 2.4 S-2: the probe reads the un-located step end

- **Change.** A sentence in the gate's items 6 and 10, and in the ladder's preamble: the probe reads each step's end before the step's event is located, so on a located step it reads a state after the event's time that the trajectory does not pass through; the thinnest cell there is thinner than at the event (item 10's reading is conservative), and the residual is that of a proposal the domain test accepted rather than of a state on the trajectory (item 6's 0.993 is such a reading where its step was located).
- **Not taken: a probe at the located state.** The deactivation callback applies its event at the root before any user callback runs, so the located pre-event state is visible only through a hook in `integrate` (an observer keyword with its tests, about a day, an API addition to a module that is not yet user-facing), followed by a ladder re-run (one to two hours at load), all to tighten a conservative reading.
- **Record.** The gate document only.
- **Estimate.** In §2.7's half day.

### 2.5 S-4: the sub-criteria no script measures

| §7.3 sub-criterion | Today | Choice | What changes |
|---|---|---|---|
| 3: the energy ledger "no larger than the old core's closure" | printed beside it, judged by eye | measure | `energy_ledger.jl` prints a verdict per form on V-BO-04 before its event (the new core's closure at most the old core's). Across the event the old core has no located event, so that row stays reported |
| 5(A): each deck's smallest `h/h₀` (C56) | not recorded | measure | `burnout_identity.jl` records, per moving deck, the thinnest live `h/h₀` over the saved states of both solves (the rule's and the explicit correction's) |
| 5(B): the largest activation and the active nodes | only the largest `\|c\|` | measure | `burnout_identity.jl` evaluates `Closures.activation` per node from the saved positions and the bound correction, and records its largest value and the most active nodes (a node with nonzero gain and activation) at any saved time |
| 5(B): V-BO-04's correction exactly zero | indirect (equal steps with and without it in the ladder) | measure | V-BO-04 joins `burnout_identity.jl` as a class B deck at its CI resolution: its digest in the reference, its largest `\|c\|` (exactly 0.0) and its steps without the correction |
| 6: the domain rejections "before and after C52" | `phase1_above` (D71) | restate | §7.3 item 6 becomes D71's measure: on one trajectory, the largest accepted scaled residual under C52's scale and under Phase 1's, and the number of accepted states Phase 1's scale would have rejected (none is a pass). A second scale in the package would be a legacy mode |
| 8: the 1.5 and 2 thresholds | ratios only | measure | `burnout_timing.jl` gives every matched row a verdict column: pass at or below 1.5, the owner's between 1.5 and 2, stop above 2 |

- **Files.** `burnout_identity.jl`, `energy_ledger.jl`, `burnout_timing.jl`; plan §7.3 item 6.
- **Record.** D99 (the measured sub-criteria and item 6's restatement).
- **Estimate.** One day (the identity additions half of it).

### 2.6 G-3: item 9 at reltol 1e-7

- **Change.** `CONFIGS` in `burnout_gradients.jl` gains a relative tolerance per configuration: the four of the case's tolerance stay, and the same four at reltol 1e-7 are added for the finite-volume core only. The old core is not run at 1e-7: item 9 judges the finite-volume core alone, and the old core's rows stay at the case's tolerance as the reported comparison. The md's section titles name the tolerance.
- **Why all four.** C93's recommended answer (a) rests on the absorption coefficient's two configurations; answer (c), reading item 9 at 1e-7, needs all four. The other two cost about ten minutes more, so the table holds whatever the owner answers.
- **Runtime.** 36 more finite-volume solves of V-BO-01 (27 at 60 cells, 9 at 120), each about one and a half to three times its 1e-6 time: under half an hour at the review's load, and the whole script under an hour. No machine condition (spreads and σ, not times).
- **Record.** D100 (old: C93's 1e-7 numbers from a scratch run; new: a committed table).
- **Estimate.** Half a day.

### 2.7 G-2 and G-4: the gate document's numbers and wording

- G-2: item 6's largest accepted scaled residual is 0.993 under C52's scale (V-BO-01, 60 cells, temperature form, `ε_m = 3e-3`; `burnout_ladder.csv` column `c52`) and 0.994 under Phase 1's (column `phase1`), not 0.998; in the summary, §6 and plan §12.1's M6c row. The verdict stands.
- §8: "1,620 steps on `ale_darcy` 100, against the decision review's 1,623 (R41)" for "the same steps the review measured", unless item 8's run changes it.
- §1: V-BO-04 "one invocation of 39 cells (1.0 %)" for "1 of 39 cells".
- §2's burst table: V-BO-03 in its two forms' rows (column (c) 0.727 and 0.488 in the temperature form, 0.591 and 0.237 in the enthalpy form, both passing), not "(enthalpy alike)".
- §1: V-BO-03's invocation times (3,022.46 to 3,023.0 s, at 30 to 240 cells) cited to `burnout_burst.md`, where they are, not to the ladder.
- Item 7 in the wording of its script: it passes on every deactivating case that has an event (16 of 20 rows); the four rows of V-BO-03 and V-BO-05 have no event to time and lie outside the criterion, rather than "the rows the table counts as not passing".
- Item 4 gets a section of its own (the bullet under §1 moves there).
- With these, the items of §2.1 (scope), §2.2 (item 1's rows), §2.4 (the probe sentence) and every number §3 re-measures.
- **Estimate.** Half a day, with §2.4.

## 3. Re-measurements

### 3.1 The tables after the fixes

Every run below is made in a detached worktree of the last code commit of §5 (commit 17), signed before it is measured (G-1's rule), one script at a time, after one precompile, so no two processes precompile Pyrolysis at once.

| Table | Script | Why | Expected | Machine | Runtime |
|---|---|---|---|---|---|
| `burnout_identity.md`, `burnout_identity_reference.csv` | `burnout_identity.jl`, then `--record` | K-1 changes class B; S-4 adds V-BO-04 and the activation columns | compare mode: class A equal, V-CPL-02/02b different; after `--record`, every deck equal and V-BO-04 a new reference row | any | about 15 min |
| `case_accuracy.{csv,md}` | `case_accuracy.jl` (all of `NEW_CORE_CASES`, both forms; its case filter would write a table of those cases alone) | K-1 on V-CPL-02/02b (item 7, and its largest constraint residual, 0.999 on V-CPL-02) | V-CPL-02/02b rows new and passing item 7; every other row's steps, error and shift reproduced | any (walls reported only) | 1 h at a load of 15 |
| `energy_ledger.md` | `energy_ledger.jl` | S-4's verdict for item 3 | the numbers of M6c, with the verdict | any | 15 to 30 min |
| `burnout_gradients.{csv,md}` | `burnout_gradients.jl` | G-3 | the 1e-6 rows reproduced to their printed digits; the 1e-7 rows new | any | under 1 h |
| `burnout_ladder.{csv,md}` | `burnout_ladder.jl --item1` | S-1 | four item 1 cells flip (§2.2) | none (no solve) | seconds |
| `burnout_timing.{csv,md}` | `burnout_timing.jl` | X-1, S-3, S-4, and K-1 on the installed rows | §3.2 | quiet | §3.2 |
| `burnout_burst.{csv,md}` | none | no fix reaches the V-BO decks | stamped `5f08134`, read through the mapping | | |

### 3.2 Item 8 on a quiet machine (X-1)

- **The script change.** `deck_rows` takes a deck descriptor (name, cells, the problem, the old core's keywords, `setup`'s keywords, the default relative tolerance, the window's field, `min_thickness`) instead of reading `CASES` alone, so the installed-correction decks (`ale` at 25 and 100 cells to 400 s and `ale_darcy` at 25 and 100 to 300 s, deactivation on, default reltol 1e-6; V-CPL-02 and V-CPL-02b at their CI resolutions and tolerances) run the V-BO rows' protocol: each core's ladder (reltol 1e-5 to 1e-7) against its own reference at 1e-9, the window before the first event of any solve of the deck (the whole run where no solve has an event), the matched tolerance the loosest rung whose window RMS is at or below the old default's, the old default timed over nine solves. The installed section then shows the matched ratio and verdict (§2.5), the equal-tolerance ratio, and the same deck with `NoCorrection()` (steps, outcome, end).
- **R41's cost apart from the integrator's.** Before the first activation `c` is exactly zero (the relative-flux faces are skipped, and V-BO-08's substrate showed the trajectory bit-identical to the uncorrected one up to that step), so the new core there solves the problem without a correction on the wider pattern, and its ratio to the old core over that interval is the integrator's and the per-step cost's. After it the gain `N²/τ` acts. Per matched row, on the installed decks and on V-BO-01 and V-BO-03: `t_a`, the first saved time with a nonzero correction velocity (the old core has none); each core's accepted steps before and after `t_a`; the split of each core's time at `t_a`, taken from a read-only step recorder in a separate solve and applied as a fraction to the timed median, so the recorder's own cost stays out of the ratio; the ratio before and after `t_a`; and the gain's share of the new core's time, `(T_after − r_before · T_old,after) / T_new`. One diagnostic row, reported and not gated: `ale_darcy` 100 at `τ = 16 × 720` s, whose per-node gain equals `ale_darcy` 25's at the default, with its steps, outcome and invocations, which says whether R41's steps are the `N²` itself. A per-layer `τ` is not designed here.
- **Machine condition.** A load average below 2 at the start and end of the run with no `foamRun` or other job on the machine (`uptime`, `pgrep -c foamRun`), one thread, nothing else running in parallel, the load at start and end recorded in the table as now; a run whose end load exceeds 3 is repeated. The run waits for the owner to free the machine.
- **Runtime.** About one hour on a quiet machine (the V-BO part's timed solves were 22 minutes at a load of 17 to 21, where times run about three times the quiet ones; the references at 1e-9 and 1e-10 and the profiler are most of the rest; the installed decks add about half). Two to three hours at the review's load, which would not decide the item.
- **What the owner gets.** For each deck the matched ratio with its verdict and, beside it, how much of it is the integrator (the before-activation ratio, R40's, which MF's `FBDF` addresses) and how much the gain (R41's, which a layer-scaled `τ` addresses): §7.4's decision on that table.
- **Record.** D101 (old: the installed rows at equal tolerance, outside item 8's protocol; new: matched on the window RMS, the cost split at the first activation; reason: X-1).
- **Estimate.** One and a half days for the script, plus the run.

## 4. Process

### 4.1 P-1: pull request #14

The owner pushes and edits; nothing here runs `git push` or `gh pr edit`.

- **Push** `mesh-framework-phase2` (47 signed commits from `0b14741`) now: independent of the followup, and the draft then shows M2 to M7.
- **Retitle** #14 "Phase 2 of the mesh framework: burnout on the finite-volume core", still a draft against `mesh-framework-phase1`.
- **Replace the body** once the followup is merged into `mesh-framework-phase2`, so that it carries the corrected numbers: plan §11's first paragraph ("What Phase 2 builds"); the gate document's summary table as the followup leaves it, with the three-way `Pkg.test()` counts; the open decisions of §7 below, answered or pending; links to the plan, the gate, the review and this plan; the line that nothing goes to `main` from it. The text is assembled at commit 25 of §5 and handed to the owner as a file.

### 4.2 R-2

After commits 23 and 24 record G-1's mapping, delete the local `docs/refs/phase2_commit_queue/` (gitignored; no commit).

### 4.3 The followup's pull request

`phase2-review` → `mesh-framework-phase2`, titled "Follow up the Phase 2 review". The body: the review's verdict in two sentences; the fixes by finding with their D and C numbers; the gate numbers that changed (item 1's four rows, item 5's class B digests, item 6's 0.993, item 7's V-CPL rows, item 9 at 1e-7, item 8's table); the owner's decisions still open. The owner pushes the branch and opens it.

## 5. Commits

One theme each, subject imperative, body stating what, why and the evidence, signed with `git commit -S`; no AI attribution. Before each: `julia --project=scripts scripts/format.jl` and the checks named. Short commits, since signing has timed out on long sessions.

| # | Subject | Finding | Files | Checks before committing |
|---|---|---|---|---|
| 1 | Add the Phase 2 review | | `docs/Mesh_Framework_Phase2_Review.md`, unchanged | |
| 2 | Plan the Phase 2 review's followups | | this document, with the owner's answers | |
| 3 | Stamp the gate tables with their tree | G-1 | `gates/stamp.jl` | the stamp printed in a scratch session |
| 4 | Share the event-shift measure between the gate scripts | S-3 | `gates/event_shift.jl`, `case_accuracy.jl`, `burnout_timing.jl` | both scripts loaded; the function on one V-CPL-02 and one old-core V-BO-01 solve |
| 5 | Count exactly one invocation per layer in item 1 | S-1 | `burnout_ladder.jl` | the script loaded |
| 6 | Recompute the burnout ladder's item 1 column | S-1 | `burnout_ladder.{csv,md}` | the diff is the four rows and the stamp line |
| 7 | Record item 5's activation, thinnest cell and V-BO-04 | S-4 | `burnout_identity.jl` | compare mode on the current reference: the eight decks equal, the new columns filled |
| 8 | Give item 3's energy ledger a verdict | S-4 | `energy_ledger.jl` | the script loaded |
| 9 | Encode item 8's thresholds | S-4 | `burnout_timing.jl` | the script loaded |
| 10 | Match the installed correction on the window RMS and split its cost at the first activation | X-1 | `burnout_timing.jl` | one installed deck (`ale` 25) and V-BO-01 at 30 cells run through the new path in a scratch session |
| 11 | Measure item 9 at reltol 1e-7 as well | G-3 | `burnout_gradients.jl` | one configuration at 1e-7, two perturbations, in a scratch session |
| 12 | Read no dead cell in the relative-flux stencil | K-1, TS-1 | `src/Kernels/relative_flux.jl`, `test_correction.jl` | §1.1's test files; the twelve V-BO digests; `burnout_identity.jl` compare mode (class B differs, nothing else) |
| 13 | Take the band default from the bound correction | I-1 | `src/Integration/setup.jl`, `test_correction.jl` | `test_correction.jl`, `test_deactivation.jl` |
| 14 | Hold a merge's receiver to its own mass tolerance, as the domain test does | I-3 | `src/Integration/checks.jl`, `src/Integration/solve.jl`, `test_deactivation.jl` | `test_deactivation.jl`, `test_policies.jl` |
| 15 | Read a self-yield as gasifying in `vanishing` | C-1 | `src/Closures/correction.jl`, `test_correction.jl` | `test_correction.jl` |
| 16 | Restart the caloric warm start in every integration | X-3 | `src/Integration/solve.jl`, `src/Assembly/residual.jl`, `test_policies.jl` | `test_policies.jl`, `test_solve.jl`; the enthalpy-form digests of commit 12 unchanged |
| 17 | Correct the solve and deactivation docstrings to the code | I-2, D-1 | `src/Integration/solve.jl`, `src/Deactivation/events.jl` | Aqua and JET on this commit (commits 12 to 16 change internal signatures); the old-core fingerprint once |
| 18 | Re-record the identity reference on the followup's code | K-1, S-4 | `burnout_identity_reference.csv`, `burnout_identity.md` | measured on 17 |
| 19 | Re-measure event timing | K-1 | `case_accuracy.{csv,md}` | measured on 17 |
| 20 | Re-measure the energy ledgers with item 3's verdict | S-4 | `energy_ledger.md` | measured on 17 |
| 21 | Measure item 9's gradients at reltol 1e-7 | G-3 | `burnout_gradients.{csv,md}` | measured on 17 |
| 22 | Measure item 8 on a quiet machine | X-1 | `burnout_timing.{csv,md}` | measured on 17, §3.2's condition |
| 23 | Correct the Phase 2 gate document | G-1, G-2, G-4, S-1, S-2 | `docs/Mesh_Framework_Phase2_Gate.md` | every number traced to its table |
| 24 | Record the review's followups in the Phase 2 plan | | `docs/Mesh_Framework_Phase2_Plan.md` (§7.3 item 6, §9.3, §12.1 to §12.4, §14.6 row 6) | |
| 25 | Record what the followup implemented in the review | | `docs/Mesh_Framework_Phase2_Review.md` (the status column; §5 what was implemented, with gates; §6 what was not, and why; §7 the three-way `Pkg.test()`; §8 the owner's answers) | the closing gate (§8) |

If the machine is not quiet when commits 18 to 21 are done, commit 22 and item 8's part of 23 follow when it is (question 10 of §7.2).

## 6. Records

| Record | Old | New | Reason |
|---|---|---|---|
| D90 | `_stencil` takes `X = far[side, f]` whatever its state | `X` absent when dead; the receiver's face takes the outermost cell's stencil | K-1; the twelve V-BO digests equal, the class B digests and V-CPL-02/02b's QoI values listed |
| D91 | the band defaults to `CORRECTION_BAND` for any `EquidistributionCorrection` (D54) | for a bound correction with some nonzero gain; the threshold otherwise | I-1; C58's configuration as §1.2 states |
| D92 | the post-event test's tolerance for a receiver unstated | the receiver's own, as the domain test's, by design | I-3 |
| D93 | a self-yield reads as non-vanishing | the reactant's own product ignored; a full self-yield consumes nothing | C-1 |
| D94 | `T_warm` carried between `integrate` calls | reset to `T0` at each | X-3; C79 unchanged |
| D95 | `solve`'s and `deactivate!`'s docstrings: every cell within `reltol` goes; one merge | the contiguous run; the receiver re-tested | I-2, D-1 |
| D96 | item 1's column `invocations ≤ layers` | `== layers`, recomputed from the csv; four rows fail | S-1 |
| D97 | two `event_shift` functions, the timing table's blank on the new core | one, shared | S-3 |
| D98 | the stamp names the commit | commit and tree; measure on signed commits (§9.3) | G-1 |
| D99 | five sub-criteria unmeasured | four measured, item 6 restated to D71 | S-4 |
| D100 | item 9 at 1e-7 from a scratch run | a committed configuration | G-3 |
| D101 | the installed rows at equal tolerance | matched on the window RMS; the cost split at the first activation | X-1 |
| C95 | seven tables, plan §12.1 and the gate stamped `5f08134` | the mapping to `78ce699` (tree `543bbb8`) recorded | G-1; resolution D98 |

Plan §14.6 row 6 (the class B digests, re-recorded at W1, W5 and MF) gains D90. A conflict found while implementing takes the next C number, from C96.

## 7. The owner's inputs

### 7.1 Decisions the plan waits on, with the review's recommendation

*Owner:* open. None is answered with the followup, and nothing in it presumes an answer; each stays open in the gate document and in plan §12.4.

| Input | Recommendation (the review or the gate) | What waits on it |
|---|---|---|
| Item 8 under §7.4 (X-1): continue with R40 and R41 carried, or stop | decide on the quiet-machine table of §3.2; the V-BO-01 and V-BO-03 ratios at 60 and 120 cells have stood above 1.5 in every session since M4, so the lever is `FBDF` at MF (R40) for the integrator's part and a layer-scaled `τ` (R41, R27) for the gain's | the gate's item 8 verdict; whether MF and the `τ` design come before the cutover |
| C91: V-BO-07's merged fraction (items 1 and 2) | (a) an owner exception, the mechanism named | the gate's items 1 and 2 |
| C92: the old core's ledger on V-BO-06, bounded at 5e-3 | confirm | nothing |
| C93: item 9's noise | (a) pass with the noise named, reltol 1e-7 for a user who differences the invocation time; decided on §2.6's committed table | the gate's item 9 |
| C94: V-BO-08 at `3e-2` ends `:interior` | (a) an owner exception for Phase 2, R27 at the cutover | the gate's item 1 |
| V-BO-06 at `3e-2`: two invocations | extend V-BO-01's `3e-2` exception to it | the gate's item 1 |
| C73 and C74 (K-2), carried from M3 | no recommendation; neither showed on the ladder. After D90, C73's slope also serves the receiver's face after a merge | nothing in this followup |
| The design errata of plan §7.2 (R-3) | apply them at the design's next revision (§13, question 13) | nothing in this followup |
| `scripts/mesh_phase2_prototype/` (R-1) | delete the covered experiments and the README's claims, or the directory (`integrators` and `tolweights`, Appendix A.7 and A.9, have no gate script) | nothing in this followup |
| The mass-ledger QoI (X-2) | a per-step ledger, or the output-grid dependence stated in every V-BO case's comment | V-BO re-pins if the ledger changes (R37) |
| The layer time `τ` per layer (review §3 item 4) | design work for the cutover | not this followup |
| MF | `FBDF` on the upstream release, re-running items 5, 8, 9 and 12 and the suite (§14.3) | not this followup |

### 7.2 Choices this plan makes, which the owner may overturn

1. **I-1.** The band from the bound correction (§1.2). *Recommended.* The alternative, refusing an equidistribution correction that moves nothing while cells deactivate, breaks item 5(A)'s explicit-correction check wherever deactivation is on. *Changes:* with the refusal, a `setup` error and its test instead of D91's default. *Owner:* as recommended.
2. **I-3.** The per-cell test documented (§1.3). *Recommended.* *Changes:* with the summed tolerance, the post-event test sums the merged cells' tolerances, and a receiver between `−abstol` and `−k·abstol` stalls in the domain test instead of failing at the event. *Owner:* as recommended.
3. **X-3.** The reset to `T0` (§1.5). *Recommended:* bit-identical for every first integration, and every later one starts at the same point. *Owner:* as recommended.
4. **S-1.** The column recomputed from the csv, and the uncorrected V-BO-08 reference row reading "fail" with the gate saying it is outside the matrix (§2.2). *Recommended.* *Changes:* with a re-run, one to two hours at load and the ladder's stamp renewed. *Owner:* as recommended.
5. **S-2.** The sentence (§2.4). *Recommended.* *Changes:* with the probe, a day for an `integrate` hook and a ladder re-run. *Owner:* as recommended.
6. **S-3.** One shared function (§2.3). *Recommended.* *Changes:* with the column dropped, the timing table loses its shift and the header its parity claim. *Owner:* as recommended.
7. **S-4.** The table of §2.5: four measured, item 6 restated to D71. *Owner:* as the table.
8. **G-3.** All four configurations at 1e-7, on the finite-volume core only (§2.6). *Recommended.* *Changes:* with the two absorption-coefficient configurations alone, ten minutes saved and C93's answer (c) left without numbers. *Owner:* as recommended.
9. **X-1.** The cost split at the first activation and the equal-gain diagnostic row (§3.2). *Recommended.* *Changes:* without them, item 8's table cannot say which lever its ratios need. *Owner:* as recommended.
10. **The followup's pull request and item 8.** *Recommended:* if the machine is quiet within a week of commit 21, the pull request carries item 8's table; otherwise it goes to review with item 8 pending, and the quiet run follows as one table commit and one gate commit on `phase2-review` before the merge. *Changes:* with "wait", the merge waits for the machine. *Owner:* as recommended.

## 8. Closing gate

1. **Scope of the change.** `git diff --stat mesh-framework-phase2..phase2-review -- src` lists only `src/Kernels/relative_flux.jl`, `src/Integration/setup.jl`, `solve.jl`, `checks.jl`, `src/Closures/correction.jl`, `src/Assembly/residual.jl` and `src/Deactivation/events.jl`; no old-core file.
2. **Quality.** The format check clean, Aqua 9 of 9 and JET with no Pyrolysis-owned reports, on the last source commit.
3. **The old core.** `oldcore_fingerprint.jl` equal to 17bc477's digests, once (no old-core file changes, so no old-core `Pkg.test()` beyond the three-way run).
4. **The targeted tests** per commit as §5 lists them, their counts in the review's §5 (test_correction.jl from 673, test_deactivation.jl from 536, test_policies.jl from 202, each grown by its new testsets).
5. **The tables.** `burnout_identity.jl` in compare mode equal to the re-recorded reference on every deck, class A unchanged since M6c; `case_accuracy.md` passing item 7 on every row with an event, V-CPL-02/02b's rows new and every other row's steps and errors as at M6c; `energy_ledger.md` with item 3's verdict; `burnout_gradients.md`'s 1e-6 rows as at M6c and its 1e-7 rows; `burnout_ladder.md`'s recomputed item 1 column; `burnout_timing.md` on a quiet machine (or pending, per question 10). Each stamped with commit 17 and its tree.
6. **The three-way `Pkg.test()`** in three detached worktrees of the last commit that changes source or tests, precompiled once before the three start: the old core, `PYROLYSIS_CORE=new` in the temperature form (the verification suite and every `test/newcore` file), and in the enthalpy form (the verification suite); counts and wall times in the review's §7. M7's counts were 13,598, 10,004 and 543; the first two grow by the new tests, the third stays.
7. **Records.** Plan §12 (D90 to D101, C95, a §12.1 row for the followup with its commits and gate evidence, §12.4's next step), §7.3 item 6, §9.3, §14.6 row 6; the gate document; the review's status column and §5 to §8.

## 9. Estimate

| Group | Days |
|---|---|
| Code fixes (§1; commits 12 to 17) | 3.5 |
| Gate scripts (§2; commits 3 to 11) | 4 |
| Re-measurements (§3; commits 18 to 22), work only | 1.5 |
| Documents and records (commits 23 to 25, §4's texts) | 2 |
| Closing gate | 0.5 |
| **Total** | **11.5, about two and a half weeks** |

Machine time beside it: about three hours at a load of 15 for commits 18 to 21 and the three-way `Pkg.test()` (about two hours side by side, the old core's run the longest), and about an hour on a quiet machine for item 8, whose wait is not in the estimate. Plan §7 gave the review and its followups one to two weeks; the followup is larger because item 8's protocol is extended to the installed decks (§3.2) and four sub-criteria gain measurements (§2.5).
