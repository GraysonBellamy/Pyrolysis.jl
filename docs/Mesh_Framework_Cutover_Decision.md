# Mesh framework: the cutover decision

**Scope.** The evidence pack of the Phase 3 plan (`docs/Mesh_Framework_Phase3_Plan.md` §2): what the owner reads to decide whether the finite-volume core replaces the old core. It is assembled at M0 on branch `mesh-framework-phase3`, stacked on `phase2-review` at `3bf19f2`, from the Phase 1 and Phase 2 gates (`docs/Mesh_Framework_Phase1_Gate.md`, `docs/Mesh_Framework_Phase2_Gate.md`), their measurement tables under `test/newcore/gates/measurements/`, the decision review of 2026-10-20 and the independent review of 2026-10-21 (both local, under `docs/`), extended with the multi-D probe's result after M1, and kept as the record of the decision. Every number names its table; nothing is measured here that a gate script did not write, except where a section says so.

**Status.** Sections 1 to 5 are complete. Section 6 (item 8 on a quiet machine) waits for the machine, by the owner's answer. Section 7 holds the owner's M0 answers; question 1, the cutover itself, waits for item 8. Section 8 records the minimal numerical workflow. Section 9 is written after M1. Section 10 is the decision.

---

## 1. What the migration buys

The design argues the migration on its own terms (design §6, Phase 2: "the exact Jacobian at reduced cost, one-callback deactivation and the multi-D readiness"), not on a speed-up.

| claim | evidence | source |
|---|---|---|
| An exact Jacobian, cheaper than the old exact one on moving meshes | one colored Jacobian against one old `DirectSolve` assembly: 0.32 to 0.36× at 25 cells, 0.11 to 0.12× at 100, 0.02 to 0.03× at 400 (`ale`, `ale_darcy`, both forms); on the Eulerian deck 1.11 to 1.44×. Against the old reduced default `ApproxSolve`, which is not exact: 1.29 to 2.01× | `leftover1_jacobian_cost.md` |
| One invocation per vanishing layer | every V-BO-01, 03, 04, 06 and 08 run at `ε_m` of 3e-3 and 1e-2, largest merged fraction 0.57 to 4.53 %; the exceptions of §4 | Phase 2 item 1, `burnout_ladder.md` |
| Located event times | the located invocation-time shift 2.7e-8 to 6.0e-3 of the interval between invocations, 10 to 4.9e7 times smaller than without location, on every deactivating case with an event | Phase 2 item 7, `case_accuracy.md` |
| Smooth event times under parameter changes | the invocation time's central-difference spread at most 0.145 at the case's tolerance and 0.0101 at reltol 1e-7; the old core's first park 1.4 to 109 | Phase 2 item 9, `burnout_gradients.md` |
| A smaller end-of-layer burst | V-BO-01's per-step excess of N over 2N 5 to 10× smaller than the old core's at `ε_m` = 1e-2, 1.45 to 1.7× at 3e-3 | Phase 2 item 2, `burnout_burst.md` |
| Decks that finish | `ale` 100, `ale_darcy` 25 and 100 complete with the installed correction and end `:interior` without it | Phase 2 item 8, `burnout_timing.md` |
| Exact mass conservation | masses are the states: V-CPL-03's total mass to 1.6e-15 against the old core's 6.6e-7; every burnout run's mass ledger over the accepted steps at most 3.1e-6 against a tolerance of 3e-5 | Phase 1 plan §12 (M6), Phase 2 item 3 |
| A tighter energy balance | the temperature form's closure 15 to 38× smaller than the old core's on V-CPL-06's column; the inventory form's exact to round-off with the integrator's quadrature | Phase 1 item 7, `energy_ledger.md` |
| A trial-state contract | no hard clamp in the residual; inadmissible proposals rejected by the domain test; the independent review's four defects fixed with regression tests (D105 to D108) | Phase 2 gate, decisions |
| Dimension-general interfaces | `Face{D,K}`, `Cell{D,F,K}`, `MeshTopology{D}`, `Mesh{D}`, `MeshGeometry{D,T}` are general; the layout, the kernels, the closure and the pattern are 1D by construction (plan §1.2, row 5). Not evidence until M1's probe | §9 |

## 2. What it costs

- **Wall time on the burnout path.** At matched accuracy on a loaded machine (load 17.5 to 21): V-BO-01 1.00, 1.85 and 2.11× and V-BO-03 1.57, 1.88 and 2.09× the old core at 30, 60 and 120 cells; V-BO-04 0.97 to 1.09× (`burnout_timing.md`, M6c's table). On V-BO-01 at 60 cells the colored Jacobian takes 39 % of the time (1,040 µs per evaluation against the old structured Jacobian's 368 µs), the right-hand side 32 %, the factorization and solves 11 %; `Rodas5P` builds an exact Jacobian every step (R40). The quiet-machine re-measurement is §6.
- **The correction's cost on long single layers.** At equal tolerance, 3.73× the old core on `ale` 100 and 2.38 and 4.62× on `ale_darcy` 25 and 100: the `N²/τ` gain (R41). These decks do not finish without the correction (§1).
- **The audit decks.** Phase 1 item 2 at matched accuracy: 0.99 to 1.88×; `ale` 100 fails 1.5× at 1.88, an owner exception (C40). At equal tolerance the cores take about the same time (0.75 to 1.19×), but `Rodas5P`'s error is 3 to 22× the old default's on five of six decks.
- **The inventory form.** On the V-HC cases 2 to 12× the temperature form's steps; on the audit decks 2.9 to 6.3× the old core at matched accuracy (Phase 1 gate). It is not the default.
- **`FBDF` is not counted.** The decision review measured the burnout path under `FBDF`, with the upstream domain-rejection bug patched in process, at about half `Rodas5P`'s cost and below the old core on every deck run. No OrdinaryDiffEqBDF release contains the fix (SciML/OrdinaryDiffEq.jl#4662 merged, unreleased; #4695 open); by the independent review's condition, an anticipated speed-up is not a passed gate.
- **Open accuracy limits.** V-BO-07's burst criterion (a) falls 1.19 and 1.15× per doubling against 1.5 (C91); V-BO-08 at `ε_m` = 3e-2 ends `:interior` (C94); one constant `τ` = 720 s keeps a layer in one invocation up to about twice its calibrated layer time (C98, R27).
- **For users.** The public API changes without a shim (plan Appendix A.4): `solve`'s keywords, the mesh constructors, the solution's accessors, `DepletionLimiter`'s behaviour. About 6,400 User Guide and 4,900 Technical Reference lines are rewritten and 3,000 edited (plan Appendix A.3).
- **Phase 3 itself.** 20 to 24 weeks of one developer with every workstream; 15 to 18 with the layer time and the envelope's comparator deferred (plan §6).

## 3. The independent review's six conditions

| # | condition | status at M0 | where Phase 3 answers it |
|---|---|---|---|
| 1 | the four defects have regression tests; the cases pass in both forms | done: `test_relative_flux.jl` (face state in [0, 2p]; an empty cell not drained), `test_policies.jl` (at or below absolute zero inadmissible), `test_gas_transport.jl` (`_pressure_tie`; `NoGasTransport` structural), `test_correction.jl` (`tied` containment); 543 of 543 in each form on `f54a6fb` | plan §7 item 12(1) |
| 2 | a stated burnout accuracy envelope under joint refinement | open | W9, M5; plan §7 item 12(2) |
| 3 | item 8 on suitable hardware with cost breakdowns | open: the machine has not been quiet | §6 here; W5 at M6 |
| 4 | moving 2D and 3D demonstrations before deletion | open | W4, M1; plan §7 item 6; §9 here |
| 5 | a qualified parameter-to-objective derivative, or finite differences stated as the route | open | W10, M5; plan §7 item 12(5) |
| 6 | an automated numerical test gate in CI | started: §8 | W7; plan §7 item 10 |

## 4. The owner's open exceptions

- **Phase 1 item 2** (C40): `ale` 100 at 1.88× at matched accuracy; Phase 1 closed with it recorded.
- **Phase 2 item 1**, each with its mechanism named: V-BO-01 and V-BO-06 at `ε_m` = 3e-2, two invocations (C78); V-BO-07, 5.82 and 6.06 % merged (C91); V-BO-08 at 3e-2, `:interior` (C94); V-BO-01 at `τ` = 3600 s, two invocations (C98).
- **Phase 2 item 2**: V-BO-07's burst criterion (a) on its rate (C91).
- **Phase 2 item 9**: passes with the noise named; reltol 1e-7 for a user who differences the invocation time at a relative step below 2e-3 (C93).
- **Phase 2 item 8**: undecided (§6).

## 5. The stop conditions

From plan §2. The cutover is not approved, and Phase 3 does not proceed past M0, if the quiet-machine item 8 reads above 2× on V-BO-01 or V-BO-03 at the CI resolution and the owner grants no exception (P2 §7.4), or if the owner reads the burnout path's cost as "the bad sign for the architecture" of P2 §5.7. Phase 3 stops after M1 if the probe shows that the layout or the assembly must change in a way that invalidates the 1D core's gates (plan §7 item 6). The alternatives are then P2 §7.4's: R40's levers (ii) to (v), waiting for MF, or reconsidering the cutover.

## 6. Item 8 on a quiet machine

**Pending: the owner waits for a quiet window.** D101's protocol: `burnout_timing.jl` in a detached worktree of `phase2-review`'s head, at a load below 2 with nothing else running, both cores in one process, the installed decks matched on the window RMS and each matched row's cost split at the first activation into the integrator's part (R40) and the gain's (R41). At M0's start the load was 14.6 (fourteen OpenFOAM processes); the owner expects the machine free from late morning on 2026-10-23. The run takes about an hour.

```
git worktree add --detach ../Pyrolysis.jl-item8 phase2-review
cd ../Pyrolysis.jl-item8
julia --project=test --startup-file=no test/newcore/gates/burnout_timing.jl
```

The table is Phase 2's, so it lands on `phase2-review` (#15): `test/newcore/gates/measurements/burnout_timing.md` in one commit, the Phase 2 gate's item 8 in a second; `mesh-framework-phase3` is then rebased onto it and the reading entered here. M0 changes no source, so the measurement is the same on either branch.

## 7. The owner's answers at M0

The questions of plan §12 that W2's ports depend on, with the plan's recommendation and the owner's answer.

| # | question | recommendation | answer |
|---|---|---|---|
| 1 | the cutover | go on this pack, confirmed on M1's probe | waits for item 8 (§6); M1 starts before it, by the owner's choice |
| 2 | branching (C106) | merge the stack into `main` once the minimal workflow runs, both cores present, the old one the default | **the alternative**: `main` stays at `17bc477`; Phase 3 stacks on `phase2-review` as `mesh-framework-phase3`, and everything lands through #12 at the end |
| 3 | `Experimental` (C107) | delete P1 and h-AMR; P1 described as a limitation | as recommended |
| 4 | `Diagnostics` (C108) | delete; a user-facing `ConservationRecorder` and `conservation_audit` on `Outputs`; TRG chapter 16 rewritten | as recommended |
| 8 | progress bar and `parameter_sweep` (C111) | port both onto the new core | as recommended |
| 9 | `DepletionLimiter` (C110) | the one-sided limiter takes the name | as recommended |
| 10 | the mesh API (C102) | constructors without `Val(NC)`, coordinates on the mesh, wired into `setup` at M2 | as recommended |
| 18 | the removed keywords | `warm_start`, `verify_jacobian`, `diagnostics`, `min_cells`, `jacobian`, `integrator`, `extra_kwargs` removed without replacement | as recommended |
| 20 | the dependency removals | `Sparspak`, `IncompleteLU`, `Krylov`, `CommonSolve`, `RecursiveFactorization`, `FiniteDiff` removed by `Pkg.rm` | as recommended |
| 21 | the burnout path's hard gate (W9) | accuracy under joint refinement; the invocation count reported | as recommended |

## 8. The minimal numerical workflow

`.github/workflows/tests.yml` (W7's first step) runs on pull requests to `main` and to the stack's branches (`mesh-framework-*`, `phase*-review`), and by dispatch once it is on `main`: every `test/newcore` file and every case of `NEW_CORE_CASES` at its CI resolution in the temperature form, without the convergence spot checks (`PYROLYSIS_CORE=new`, `PYROLYSIS_TEST_FAST=1`), as two jobs side by side (plan D110), each printing the commit, its tree, the pull request's head and the Julia and manifest identifiers before it runs. Since the stack stays off `main` (question 2), the workflow gates the pull requests into `mesh-framework-phase3` and that branch's own pull request, and #12 when it carries everything.

Local measurement of the partition: each job's command (`Pkg.test()` with its group) on `3bf19f2` with the M0 changes, the two side by side at a load of 16 to 17. The verification cases pass 521 of 521 in 14.7 min of tests (15.0 min wall); the `test/newcore` files pass 9,696 of 9,696 in 70.7 min (71.0 min wall). Together they are the temperature-form new-core run's 10,239 tests less the 22 assertions of the 11 convergence checks. The `newcore` job is the long pole: its largest files, timed one by one in a separate process at the same load, are `test_policies.jl` 12.8 min, `jacobian/test_containment.jl` 10.9, `test_mesh_motion.jl` 7.1, `test_deactivation.jl` 6.6 and `test_correction.jl` 5.9.

On the runner, the first run (#16, run 39789492297, on `a1bb27d`, the dependency cache cold) passed both jobs. The cases job took 15 min wall (8.5 min instantiating and precompiling, 4.7 min of tests). The `newcore` job took 36 min wall (9.3 min instantiating, 24.8 min of tests). The runner ran the tests about three times faster than the loaded local machine. The `newcore` job is 6 min over M3's 30-minute budget (W7's second step), before the cache warms.

## 9. The multi-D probe

Written at M1's gate: the nine criteria of plan §7 item 6, each with its result, and whether the decision is confirmed.

## 10. The decision

Recorded by the owner.
