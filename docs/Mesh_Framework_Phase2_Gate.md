# Mesh framework: the Phase 2 gate

**Scope.** The Phase 2 gate of the mesh-framework design (`docs/Mesh_Framework_Design_2026-10-17.md` §6, Phase 2; `docs/Mesh_Framework_Phase2_Plan.md` §7.3, items 1 to 12 as the decision review restated them, §14.5), measured on branch `mesh-framework-phase2` at `5f08134` (the M6c measurement pass, plan §14.4). That hash is the commit as prepared before signing, and the tables are stamped with it; the signed commit "Time the correction where the installation rule installs it" carries the same tree, `543bbb8`. Each item gives what was measured, the numbers, and whether it passes. The tables behind the numbers are under `test/newcore/gates/measurements/`, each stamped with the commit it was measured on, and the plan's §12 holds the implementation record (D53 to D89, C57 to C94).

**The integrator.** No registered OrdinaryDiffEqBDF release contains the upstream fix of the domain-rejection bug (SciML/OrdinaryDiffEq.jl#4687: #4662 merged and unreleased, #4695 open; the latest release is 2.4.13). By the owner's instruction Phase 2 closes on `Rodas5P` (plan §14.8, question 17's option (a)). Every table is taken on `Rodas5P` with the heuristic restart and event location on, item 12 does not apply, and MF, the `FBDF` switch, is the first milestone after the release, before the cutover, re-running items 5, 8 and 9 and the suite (§14.3).

**Machine.** One shared machine (Intel Core Ultra 7 155H, 22 logical CPUs), Julia 1.13.0, one thread per process. Fourteen OpenFOAM processes held the load at about 15 throughout, and the pass ran its gate scripts and the three-way `Pkg.test()` side by side, so each table records its own load (17 to 33). Every comparison of the two cores ran in one process. No quiet machine was available, so item 8's wall times are recorded and not decided.

**Decisions taken for the gate** (the owner; plan §12.4, §13, §14.8):

- The burnout path: the block-and-band rule with the equidistribution correction; the slot-shifting interior merge stays the unbuilt fallback (item 11).
- The correction is installed by `setup` where the mesh moves, deactivation is on and a layer can vanish (C85, C87), with zero gain on a layer that cannot vanish (D83); the material residual-fraction field is dropped (D69).
- The reactivity switch is deleted (§14, W1).
- Event location is on by default (D78); the heuristic restart after an event is kept (D61, D80).
- The integrator stays `Rodas5P` until the upstream release; the gas rows keep the solid's tolerance scale until MF (M2 decision 1, C89).
- The sparse factorization learns the union of `W`'s nonzeros (D62).
- Item 8 matches on the window RMS (C70); the burnout path was accepted at M2 at about 1.5× (R40).
- The peak QoI reads the main peak (C86), and the burst is judged by refinement on both cores (§14, W4).
- V-BO-01 at `ε_m = 3e-2` (two invocations) is an owner exception of item 1 (M4 decision 1).
- The caloric inversion keeps its warm start (C79).

## Summary

| item | what | result |
|---|---|---|
| 1 | one invocation per layer, the largest merged fraction at most 5 `ε_m` | passes on every V-BO-01, 03, 04, 06 and 08 run at `ε_m` of 3e-3 and 1e-2 and at `τ = 360` s; V-BO-01 at 3e-2 is the owner's exception; **fails** on V-BO-06 at 3e-2 (two invocations), V-BO-01 at `τ = 3600` s (two, C78), V-BO-07 (5.82 and 6.07 %, C91) and V-BO-08 at 3e-2 (`:interior`, C94): owner's decisions |
| 2 | the cases' QoIs on every run; the burst by refinement | every run passes its case except the two-invocation rows (their invocation QoI) and V-BO-08's `:interior` rows; the burst criteria pass on every gated family in both forms; reported: V-BO-07's per-step excess stops falling at 240 cells (C91) |
| 3 | ledgers | passes: every V-BO run's mass ledger at its case's tolerance; V-BO-04's energy ledger 3.65e-6 and 3.98e-6 before its event (old core 7.1e-5), 1.85e-5 and 1.78e-5 across it |
| 4 | admissibility at every event | passes: no inadmissible result on 50 runs |
| 5 | charring decks identical, split by installation | passes: the eight decks equal their references; class A unchanged under an explicit correction |
| 6 | constraint residuals and the domain test | passes: the largest accepted scaled residual 0.998, none above one under Phase 1's scale |
| 7 | event timing | passes on every deactivating case with events (16 of 20 rows; V-BO-03 and V-BO-05 have none) |
| 8 | wall time against the old core | measured on a loaded machine, recorded and not decided: at matched accuracy V-BO-01 1.00, 1.85 and 2.11× and V-BO-03 1.57, 1.88 and 2.09× at 30, 60 and 120 cells, V-BO-04 0.97 to 1.09×; the installed correction 0.80 to 1.32× on V-CPL-02/02b and `ale` 25, 2.4 to 4.6× on `ale` 100 and `ale_darcy` (R41); owner's decision after a quiet-machine run |
| 9 | smooth finite-difference gradients | the flux and the pre-exponential factor pass; **fails** on the absorption coefficient at 60 and 120 cells (the invocation time's spread 0.143 and 0.128, noise that falls to 0.004 and 0.008 at reltol 1e-7, C93): owner's decision |
| 10 | minimum cell thickness | passes: the thinnest live cell at least 0.865 `ε_m h₀`, no `DtLessThanMin` |
| 11 | the fallback | not triggered: items 1 and 2 pass on V-BO-01 at the defaults |
| 12 | Phase 1's items 2 and 9 under `FBDF` | does not apply (MF has not landed) |

The full three-way `Pkg.test()` (session counts; no log is committed), each in a detached worktree of `5f08134`: the old core 13,598 of 13,598 (117 min), `PYROLYSIS_CORE=new` 10,004 of 10,004 in the temperature form (83 min: the verification suite and every `test/newcore` file) and 543 of 543 in the enthalpy form (19 min: the verification suite; `runtests.jl` runs the `test/newcore` files in the temperature form only), side by side at loads of 17 to 33.

## 1. One invocation per layer

Table: `measurements/burnout_ladder.md` (script `burnout_ladder.jl`), 25 configurations in each energy form, the band at 0.1, the monitor without its mass term unless stated.

| deck | cells | `ε_m` | invocations (largest merged fraction), both forms | item 1 |
|---|---|---|---|---|
| V-BO-01 | 30, 60, 120 | 3e-3 | 1 (0.61 to 0.64 %) | pass |
| V-BO-01 | 30, 60, 120 | 1e-2 | 1 (3.84, 4.08, 4.19 %) | pass |
| V-BO-01 | 30, 60, 120 | 3e-2 | 2 (9.48 to 9.94 %, then 4.35 to 4.48 %) | owner's exception |
| V-BO-01, `τ = 360` s | 60 | 1e-2 | 1 (1.22 %) | pass |
| V-BO-01, `τ = 3600` s | 60 | 1e-2 | 2 (9.05 %, then 1.29 %) | fail (C78) |
| V-BO-03, to 3,100 s | 30, 60, 120 | 1e-2 | 1 (1.0 %, at 3,022.5 to 3,023.0 s) | pass |
| V-BO-04, with and without the correction | 40 | 1e-2 | 1 of 39 cells (1.0 %) | pass |
| V-BO-06 | 60 | 3e-3, 1e-2 | 1 (0.57, 1.6 %) | pass |
| V-BO-06 | 60 | 3e-2 | 2 (9.8 %, then 6.8 %) | fail |
| V-BO-07, `β = 0` and `β = 1` | 60 | 1e-2 | 1 (5.82, 6.07 %) | fail (C91) |
| V-BO-08 | 10 | 3e-3, 1e-2 | 1 (1.48, 4.53 %) | pass |
| V-BO-08 | 10 | 3e-2 | 0: `:interior` at 67.04 s | fail (C94) |

- **The 3e-2 rows.** At that threshold the band (0.1) is 3.3 `ε_m`, so the first cell crosses before the main peak while most of the layer holds more than the band; the first invocation takes the cells within the band and the second the rest (C78). V-BO-06 behaves as V-BO-01 there. The owner's exception names V-BO-01's row; V-BO-06's is the same mechanism on the graded mesh.
- **V-BO-07 (C91).** The burn runs down from the surface. The upper cells lose their PMMA first and keep the slow second solid, while the nearly pure PMMA beneath burns last and fast, so when the first cell crosses `ε_m` the others still hold up to 5.8 % (`β = 0`) or 6.1 % (`β = 1`) of their mass. Neither `τ` (90 to 720 s) nor the monitor's mass term brings it within the clause. With the second solid at PMMA's kinetics, so that only the density varies, the layer goes at 2.92 % and 2.37 %, within it.
- **V-BO-08 at 3e-2 (C94).** The 1 mm layer's activated life (its first cell below `r_hi` at about 55 s, the layer gone at 69 s) is short against the correction's synchronization time `τ/π²` = 73 s, so at the largest threshold a cell beneath the surface reaches it first. The table's item 1 column counts at most one invocation per layer and reads this row as passing; it fails item 1 and its case's QoIs. At the default threshold V-BO-08 passes, at 4.53 % against 5 % (R46's margin).
- **Item 4.** No run's deactivation was refused as inadmissible, and every invocation's events name its cells.

## 2. The cases' QoIs and the burst

Tables: `measurements/burnout_ladder.md` (QoIs per run) and `measurements/burnout_burst.md` (script `burnout_burst.jl`).

**The QoIs.** Every run of the matrix passes all its case's QoIs in assert mode, except the two-invocation rows above, which fail only their invocation QoI (V-BO-07 asserts one invocation without the bound and passes, C91), and V-BO-08's two `:interior` rows (at 3e-2, and the reference run without the correction), which end before the layer burns out. V-BO-01's main peak, its time and the burnout time against the references 32.39 (3 %), 335.0 s (3 %) and 363.5 s (1 %):

| cells | peak, `t_peak`, `t_bo` at 1e-2 (old core) |
|---|---|
| 30 | 32.248, 333.0, 363.0 (32.254, 333.0, 364.0) |
| 60 | 32.391, 335.0, 363.0 (32.392, 335.0, 363.5) |
| 120 | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) |

V-BO-06 reads 32.288 at 334.0 s and 364.0 s (the old core 32.292, 334.0, 364.0), V-BO-07 30.542 at 287.5 s and 324.0 s against its references 30.64, 288.0 and 325.0 (D87), and V-BO-08 29.712 at 47.2 s.

**The burst by refinement.** The gas-generation rate at every accepted step and post-event state, at 30, 60, 120 and 240 PMMA cells on both cores, from each deck's old-core 240-cell peak time to its end.

| family | (a) per-step excess of N over 2N [g m⁻² s⁻¹] | (b) successive L1 [g/m²] | (c) L1 against the old core's 240 | verdict |
|---|---|---|---|---|
| V-BO-01, old core | 31.9, 7.73, 3.18 | 25.0, 9.34, 4.40 | | |
| V-BO-01, temperature, 1e-2 | 3.16, 1.36, 0.642 | 22.0, 10.2, 5.0 | 39.7, 19.4, 10.0, 5.93 | pass (gated) |
| V-BO-01, temperature, 3e-3 | 19.1, 5.14, 2.19 | 22.7, 11.7, 5.86 | 36.7, 18.4, 8.21, 3.21 | pass (gated) |
| V-BO-01, enthalpy, 1e-2 | 3.17, 1.35, 0.636 | 22.0, 10.2, 5.01 | 39.8, 19.4, 10.0, 5.92 | pass (gated) |
| V-BO-01, enthalpy, 3e-3 | 18.7, 4.98, 2.10 | 22.7, 11.6, 5.81 | 36.8, 18.4, 8.33, 3.23 | pass (gated) |
| V-BO-01, temperature (enthalpy alike), 3e-2 | 4.59, 1.80, 0.909 | 21.8, 9.95, 5.21 | 38.5, 18.4, 8.9, 4.71 | pass (reported) |
| V-BO-03 extended, temperature (enthalpy alike) | 0.0069, 0.0021, 0.00055 | 8.68, 2.29, 0.739 | 11.1, 2.6, 0.727, 0.488 | pass (reported) |
| V-BO-06, temperature (enthalpy alike) | 25.6, 5.23, 2.08 | 45.0, 20.9, 9.91 | 64.1, 29.9, 14.5, 7.85 | pass (reported) |
| V-BO-07, temperature | 1.61, 0.799, 1.99 | 23.5, 11.8, 6.11 | 54.4, 31.6, 20.9, 16.3 | (a) fails (reported, C91) |
| V-BO-07, enthalpy | 1.60, 0.797, 1.18 | 23.5, 11.8, 6.10 | 54.4, 31.6, 20.9, 16.3 | (a) fails (reported, C91) |

On V-BO-01 the new core's excess is smaller than the old core's at every resolution (5 to 10× at 1e-2, 1.45 to 1.7× at 3e-3), and each family's L1 against its own 240-cell run is within 1.33× of the old core's. V-BO-07's excess stops falling at 120 → 240 cells and exceeds the old core's there (1.99 against 0.90): the merge lumps the slow solid's residue (C91). By item 2's text that run goes to the owner.

## 3. Ledgers

Tables: `measurements/burnout_ladder.md` (every V-BO run's mass ledger, a QoI of its case) and `measurements/energy_ledger.md` (`energy_ledger.jl`).

Every V-BO run's mass ledger is within its case's tolerance (1e-3, V-BO-02's 3e-3 on the new core, V-BO-06's 5e-3 on the old core, C92). V-BO-04's energy ledger:

| | temperature form | enthalpy form | old core |
|---|---|---|---|
| before the event, at the case's tolerance | 3.65e-6 | 3.98e-6 | 7.1e-5 |
| across the event, location on | 1.85e-5 | 1.78e-5 | |
| across the event, location off | 1.17e-2 | 1.11e-2 | |
| fixed-step defect over the scale (dt 0.1 to 0.025 s) | 4.2e-15 to 2.7e-13 | 8.2e-15 | |

The column stays at its initial temperature, so its fixed-step defect is round-off and shows no order (D75); V-CPL-06's Eulerian and moving columns converge at fourth order in the temperature form and are exact in the enthalpy form. Item 3 passes.

## 5. The charring decks, split by installation

Table: `measurements/burnout_identity.md` (`burnout_identity.jl`).

Every deck equals its reference digest. Class A (the rule installs nothing): V-CPL-03, V-CPL-06a, V-HC-05, V-RAD-01, V-RX-02d and the Douglas fir deck over 60 s. Solved again with an explicit `EquidistributionCorrection()`, V-CPL-03 and the Douglas fir deck take the same steps, within 0.0 and 2.0e-14 relative, with the correction exactly zero. Class B (the rule installs the correction): V-CPL-02 and V-CPL-02b, 1,134 and 638 steps against 2,092 and 1,222 without it, the largest correction velocity 5.7e-4 and 3.4e-4 m/s. Item 5 passes.

## 6 and 10. Residuals, the domain test and the thinnest cell

Table: `measurements/burnout_ladder.md`.

The largest scaled constraint residual of the accepted states is 0.998 under C52's scale on any run, and no accepted state lies above one under Phase 1's scale, so C52 removed no rejection (D71). The domain test rejects 0 to 174 proposals per run. The thinnest live cell at the end of any step, before its event, is at least 0.865 `ε_m h₀` on every run (V-BO-07: 1.7 and 1.95), and no run ends in `DtLessThanMin`. Items 6 and 10 pass.

## 7. Event timing

Table: `measurements/case_accuracy.md` (`case_accuracy.jl`).

Every case of `NEW_CORE_CASES` in both forms, the deactivating ones with location off and on. The located invocation-time shift, over the reference's median interval between invocations (or its saved interval for a single invocation):

| case | temperature form | enthalpy form | gain (off over on) |
|---|---|---|---|
| V-BO-01 | 6.0e-3 | 1.0e-5 | 76, 51,800 |
| V-BO-02 | 3.2e-4 | 1.1e-5 | 378, 7,930 |
| V-BO-04 | 2.7e-8 | 2.9e-7 | 4.9e7, 6.5e5 |
| V-BO-06 | 5.3e-3 | 5.3e-5 | 99, 10,700 |
| V-BO-07 | 3.3e-3 | 1.5e-5 | 8.0, 22,800 |
| V-BO-08 | 3.1e-4 | 5.7e-4 | 3,600, 2,280 |
| V-CPL-02 | 1.8e-4 | 1.8e-5 | 173, 1,700 |
| V-CPL-02b | 5.1e-5 | 2.7e-5 | 670, 1,170 |

Every one is within 0.05 and at least three times smaller than without location. V-BO-03 (no invocation by its 500 s end, C75) and V-BO-05 (its `:interior` stop comes before any event) have no event to time and are the four rows the table counts as not passing. The largest constraint residual over every case is 0.999 (V-CPL-02). Item 7 passes.

## 8. Wall time

Table: `measurements/burnout_timing.md` (`burnout_timing.jl`), at a load of 17.5 to 21.

Matched on the window RMS before the first event (D58, C70), each core against its own reference at reltol 1e-9, the old default timed over nine solves; the ratio is the new core's warm time at its matched tolerance over the old default's.

| deck | 30 cells | 60 cells | 120 cells | at equal tolerance |
|---|---|---|---|---|
| V-BO-01 | 1.00 (reltol 1e-5) | 1.85 (1e-5) | 2.11 (3e-6) | 1.44, 1.60, 2.45 |
| V-BO-03 | 1.57 (1e-5) | 1.88 (1e-5) | 2.09 (3e-6) | 3.44, 2.72, 2.32 |
| V-BO-04 | 1.06 (1e-5) | 1.09 (1e-5) | 0.97 (1e-5) | 2.46, 2.40, 2.18 |

Where V-BO-01's time goes at 60 cells (the profiler's shares at the matched tolerance): the colored Jacobian 39 % (one evaluation 1,040 µs against the old core's structured Jacobian at 368 µs, 25,000 stored entries against 8,887), the right-hand side 32 % (39 against 32 µs per call at this load), the factorization and the solves 11 %. `Rodas5P` builds an exact Jacobian every step (R40).

The installed correction, at the case's tolerance (reltol 1e-6 on the audit decks with deactivation on), against the old core and against the same deck with `NoCorrection()`:

| deck | old core [s] (steps) | the rule [s] (steps; invocations) | rule / old | without the correction |
|---|---|---|---|---|
| `ale` 25 | 0.071 (331) | 0.070 (334; 0) | 0.99 | completes, 440 steps |
| `ale` 100 | 0.46 (699) | 1.71 (1,257; 0) | 3.73 | `:interior` at 161.8 s |
| `ale_darcy` 25 | 0.12 (569) | 0.28 (585; 1) | 2.38 | `:interior` at 84.8 s |
| `ale_darcy` 100 | 0.54 (1,107) | 2.49 (1,620; 1) | 4.62 | `:interior` at 83.1 s |
| V-CPL-02 | 1.82 (1,037) | 1.45 (1,134; 19) | 0.80 | completes, 2,092 steps |
| V-CPL-02b | 0.38 (615) | 0.50 (638; 11) | 1.32 | completes, 1,222 steps |

Without the correction three of the audit decks stop `:interior`, so the correction is what lets them finish; its cost there is the `N²` gain on a long single-material layer (R41), the same steps the review measured (1,620 on `ale_darcy` 100).

The load (17.5 to 21, fourteen OpenFOAM processes and the gate's own runs) matters here: M4's tables at a load of 15 to 16 read V-BO-01 at 1.16 to 2.28×, and two runs at that load differed by up to 0.9 in a ratio. The item is therefore recorded, not decided. On these numbers V-BO-01 and V-BO-03 lie above 1.5× at 60 and 120 cells and slightly above 2 at 120, and three installed decks above 2. Under §7.3 item 8, ratios between 1.5 and 2 are the owner's call and above 2 stop the phase for the owner. The levers are R40's (`FBDF` at MF first, then forming `W` in the union, the absent-work skips and element-local Jacobians) and R41's (the gain on long layers).

## 9. Smooth gradients

Table: `measurements/burnout_gradients.md` (`burnout_gradients.jl`).

V-BO-01 through the case's solve, nine perturbations per configuration and central differences at four steps, both cores in one process. On the finite-volume core every run ends in one invocation, none falls back, and the peak, its time and the burnout time are smooth on every configuration (spreads at most 0.041).

| configuration | invocation time: spread (σ, h10) | item 9 |
|---|---|---|
| absorption coefficient, 60 cells | 0.143 (5.5e-4 s, 1.9e-3) | fail |
| incident flux, 60 cells | 0.00097 | pass |
| pre-exponential factor, 60 cells | 0.0061 | pass |
| absorption coefficient, 120 cells | 0.128 (3.7e-4 s, 1.3e-3) | fail |

The invocation time is the least sensitive quantity under the absorption coefficient (about 3 s per unit `ln α`), so the step noise is 10 to 20 % of a central difference at h = 1e-3 and 2e-3. At reltol 1e-7 the same configurations give spreads of 0.0039 and 0.0079 (σ 2.8e-5 and 2.2e-5 s) and pass (C93). The old core's first park is erratic on every configuration (spreads 1.4 to 109). Item 9 fails on two configurations at the case's tolerance; the owner decides between reading it with the noise named, restating it on h10, or reading it at reltol 1e-7.

## 11 and 12

Item 11, the fallback, is conditioned on items 1 or 2 failing on V-BO-01 at the defaults, or item 8 failing without an exception. Items 1 and 2 pass on V-BO-01 at the defaults, and item 8 is not decided (above), so the fallback is not built. Item 12 applies only once MF lands.

## The new cases (M6b)

- **V-BO-05** is routed with its own stop (C45): the first interior `ε_m` crossing at 259.57 s ± 0.7 s, the outcome `:interior`, no event before it, the largest live temperature 693 K.
- **V-BO-06** puts V-BO-01 on a mesh graded toward the surface (D89). Its peak, timing and burnout fall within V-BO-01's references; it takes one invocation at 1.6 %.
- **V-BO-07** mixes a second solid into the layer, with references taken from the old core at 240 cells (D87). Its one invocation leaves 5.8 % (C91).
- **V-BO-08** is a 1 mm layer on the 12 mm Kaowool, run to 150 s with saves every 0.1 s (D86). It takes one invocation at 4.53 %, its correction is exactly zero at every Kaowool node, and its substrate is bit-identical to the uncorrected solve until the first activation.

`NEW_CORE_CASES` passes 543 of 543 in each form, and the V-BO cases pass 63 of 63 on the old core.

## Open questions for the owner

- **C91:** V-BO-07 fails item 1's merged-fraction clause and item 2's reported burst measure on its graded kinetics. Recommended: an owner exception with the mechanism named.
- **C92:** the old core's ledger on V-BO-06 is bounded at 5e-3. Confirm.
- **C93:** item 9's invocation time under the absorption coefficient is noise-limited at the case's tolerance. Recommended: pass, with the noise named and reltol 1e-7 for users who difference the invocation time.
- **C94:** V-BO-08 at 3e-2 ends `:interior`. Recommended: an owner exception for Phase 2, with a layer-scaled `τ` (R27) as cutover design work.
- **V-BO-06 at 3e-2:** two invocations, the same mechanism as V-BO-01's excepted 3e-2 row. Recommended: extend that exception to it.
- **Item 8:** re-run `burnout_timing.jl` on a quiet machine, then decide: V-BO-01 and V-BO-03 read 1.85 to 2.11× at 60 and 120 cells here, and the installed correction 2.4 to 4.6× on three audit decks (R41). Recommended: take the decision on MF's `FBDF` numbers, the first of R40's levers, which the review measured at 0.56× (modeled) on V-BO-01.
- **Carried from M3:** C73 (the boundary slope of the relative fluxes) and C74 (the upwind side at an exact tie).

## Carried beyond Phase 2

- **MF:** the `FBDF` switch on the release (plan §14.3), re-running items 5, 8 and 9 and the suite, with item 12.
- **R40's levers** stay at the phase gate unless MF lands (§12.4).
- **R27:** a layer time scaled to the layer (C94), if the owner wants V-BO-08's 3e-2 row or its margins at the default improved.
- **W8:** element-local Jacobians, at the start of the multi-D phase.
