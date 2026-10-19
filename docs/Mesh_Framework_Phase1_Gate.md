# Mesh framework: the Phase 1 gate

**Scope.** The Phase 1 gate of the mesh-framework design (`docs/Mesh_Framework_Design_2026-10-17.md` §6, Phase 1; `docs/Mesh_Framework_Phase1_Plan.md` §7.3, items 1 to 10), measured on branch `mesh-framework-phase1` at the M9 commits, and items 2, 6 and 8 re-measured on branch `phase1-review` at `603e863`, after the review's fixes and the owner's answers to it (the structure-keeping linear solve, the fixed-mesh vent; `docs/Mesh_Framework_Phase1_Review.md`). Each item gives what was measured, the numbers, and whether it passes; the tables behind the numbers are under `test/newcore/gates/measurements/`, written by the gate scripts beside them, each stamped with the commit it was measured on. Plan §12 is the implementation record (departures D1 to D52, conflicts C13 to C39).

**Machine.** One shared machine (Intel Core Ultra 7 155H, 22 logical CPUs), Julia 1.13.0, one thread per process, the owner's other runs holding the load at 17 to 24 throughout; every comparison of the two cores was run in one session, and each table records its load. Timings are warm medians of three (`w_formation.jl`: warm minima of three), and the same configuration's time differs by up to 1.9× between the gate scripts' runs.

**Decisions taken for the gate** (the owner, plan §12.4):

- Integrator (C28): `Rodas5P` with the solver-chosen restart is the default; the finite-slab criterion is read at the spatial floor (D37).
- Reactivity switch (C27): kept at 1.5 `ε_m` for Phase 1; revisited on the V-BO ladder in Phase 2.
- V-MMS-02/03 (C24): run with their pins kept as exactness checks. The owner chose tighter tolerances; under `Rodas5P` every pin holds at the cases' own tolerances, so none were applied (D35).
- Darcy upwind (C21): the first-order upwind stays. Non-finite inventory (C22): applied (D33).
- Burnout cases (C29): V-BO-02 joins the new core with its step budget re-pinned; V-BO-01, 03, 04 and 05 move to Phase 2.
- Linear solve (C32): `StructuredKLU`, keeping KLU's analysis across a solve, is the default on a layout with algebraic rows (D46).
- Fixed-mesh transfer (C34): the dead cell's gas is vented to the event's ledger (D47); V-BO-02's ledger quantity on the new core is held to 3e-3 (C39).
- V-MMS-03's face temperature (C33): the upwind-biased reconstruction stays.

## Summary

| item | what | result |
|---|---|---|
| 1 | verification cases within their tolerances | passes for every case of §6.2; the six cases of C7 and C29 lie outside its scope |
| 2 | audit decks within 1.5× of the old core at matched accuracy | with the structure-keeping linear solve, temperature form 0.99 to 1.88× over two runs (Eulerian 1.16 to 1.45, ale_darcy 0.99 to 1.38, ale 25 1.48 and 1.51, ale 100 1.39 and 1.88), against 1.9 to 2.8× on the moving-mesh and Darcy decks before; **fails** on ale 100 and sits at the threshold on ale 25 |
| 3 | MMS orders and the V-MMS-01 temporal sweep | passes; V-MMS-03's temperature converges below second order on the new core (1.6 to 1.7), its species at second order (C33) |
| 4 | radial bandwidth and colors equal to planar | passes |
| 5 | the suite with `H` as the state | passes |
| 6 | constraint residuals below tolerance on every case | passes (largest 0.42; V-RAD-01's exact zero is exact) |
| 7 | temperature-form energy balance no worse than the old core's, at the integrator's order; inventory form exact | passes |
| 8 | achieved error against requested tolerance on every case | 60 of 66 solves within budget; the deactivating cases' event times follow the step, not the tolerance; V-CPL-05 within budget over the times its reference reaches |
| 9 | V-HC-03 and V-HC-02c at 320 cells with the chosen integrator | passes (at the spatial floor) |
| 10 | leftovers 1 and 2 recorded | recorded; the affordability decision is open (item 2) |

The full three-way `Pkg.test()` (session counts; no log is committed), each in a detached worktree of the commit: the old core 11,639/11,639 (83 min, at `8788d48`), `PYROLYSIS_CORE=new` 8,212/8,212 in the temperature form (58 min) and 8,212/8,212 in the enthalpy form (59 min, both at `6c0346a`). At `8788d48` the two new-core runs stopped at `test/newcore/test_policies.jl`, which imported SciMLBase, a package the test environment does not carry (the development sessions do); `6c0346a` imports it through OrdinaryDiffEq and changes nothing else, and test_policies passes 113/113 in the test environment. The three-way run on `phase1-review`'s final code is recorded in the review document's §6.

## 1. The verification cases

`NEW_CORE_CASES` (22 case families: V-RX-01..04, V-BO-02, V-CPL-01/01b, V-CPL-02/02b, V-CPL-03, V-CPL-04, V-CPL-05, V-GT-01/02, V-MMS-01..04, V-RAD-01, V-HC-01..05) through the suite's own routing (`test/verification/casedefs.jl`), in assert mode with the convergence spot checks and the wall-time budgets: every quantity of interest and every convergence check passes in both energy forms, 188 quantities and 11 convergence checks per form (counted in a run over the 33 routed cases). The old core, in the same session, passes every case of the suite through the same files (the step recorder and the routing changed).

What the suite's routing (plan §6.1, C7 and C29) keeps off the new core; none of these cases is in §6.2, item 1's list:

- V-CPL-06a/b read the old core's internals (C7). Their obligation, the energy row against the ledger with gas transport on, is `test/newcore/test_energy_ledger.jl` (102/102), and item 7 compares the two cores on V-CPL-06's own column.
- V-BO-01, 03, 04 and 05 end at an interior crossing (C29): a cell below the threshold inward of a live one, which design §3.5 ends as `:interior` and the old core parks. V-BO-01 and 03 (in-depth absorption on a moving mesh) end at 242 s before any event; V-BO-04 (uniform zero-order burnout, every cell at once) at 468 s after one event; V-BO-05 ends as intended at an interior void, at 263.5 s against the pinned 330 to 345 s, because `ε_m` is a mass fraction (C9). They need an interior-merge or cascade policy, which is Phase 2's burnout work.

Case-level changes the new core needed: V-BO-02's step budget (3,000 and 4,000 accepted steps in the temperature and enthalpy forms against the observed 2,354 and 3,163, and 1,706 and 2,456 with the fixed-mesh vent, D47; the old core keeps 2,000 against 1,572) and, with the vent, its mass-ledger tolerance on the new core (3e-3, observed 1.5e-3 and 8.5e-4: the time error of the gas transient after each park, which its 0.5 s outputs do not resolve, 4e-5 at reltol 1e-8; C39), V-BO's step recorder through the layout (C8), and the caloric extension admitting a gas that carries no heat (D32). Rerun on `phase1-review` with every change of the review and its followup: every quantity and convergence check passes in both forms, 188 and 11 per form.

Report mode (`scripts/verification/run_reports.jl` with `PYROLYSIS_CORE=new`) runs the 21 routed families unchanged, 31 case reports with their figures, CSVs and chapters, written to a scratch directory (the committed guide location is not touched in Phase 1).

## 2. The audit decks at matched accuracy

Table: `measurements/audit_decks.md` (script `audit_decks.jl`); breakdown `measurements/audit_breakdown.md` (`audit_breakdown.jl`).

The audit decks (`eulerian`, `ale`, `ale_darcy` at 25 and 100 cells, 80 s) are timed on each core's default integrator on a tolerance ladder; a run's error is its largest temperature difference against the same core at reltol 1e-11, so it is the time error alone. The target is the old core's error at its default tolerances; the new core's matched time is its cheapest ladder solve at or below it.

The two PMMA decks have no converged reference as defined (C30): the solid's gas-transfer coefficient is zero, gas transport vanishes in the solid, and neither core's solution converges with tolerance (the new core at reltol 1e-10 and 1e-11 differs by 0.02 to 0.05 K, the old core at its default and at 1e-10 by up to 0.21 K). The matched comparison uses them with a solid coefficient of 1e-6 m²/s, which converges both cores to 1e-8 K; the originals are timed at equal tolerance.

| deck | old [s] | target [K] | temperature form: matched reltol, ratio | enthalpy form: matched reltol, ratio | equal tolerance (temperature) |
|---|---|---|---|---|---|
| eulerian 25 | 0.0042 | 9.9e-6 | 1e-7, **1.16** | 1e-7, 4.23 | 0.75 |
| eulerian 100 | 0.013 | 1.4e-5 | 1e-7, **1.32** | 1e-7, 7.78 | 0.83 |
| ale 25 | 0.0046 | 2.3e-5 | 1e-7, **1.48** | 3e-7, 3.56 | 0.95 |
| ale 100 | 0.015 | 1.0e-5 | 1e-7, 1.88 | 3e-8, 13.0 | 1.19 |
| ale_darcy 25 | 0.011 | 1.3e-4 | 1e-6, **1.06** | 1e-6, 2.22 | 1.06 |
| ale_darcy 100 | 0.042 | 6.0e-5 | 3e-7, **1.33** | 3e-7, 3.36 | 0.99 |

The table is the second of two runs at `603e863` in one worktree, the load at 12 to 14; the first gave 1.45, 1.35, 1.51, 1.39, 0.99 and 1.38 in the temperature form. The new core's matched times agree between the runs (ale 100: 0.0291 and 0.0281 s); the old core's moved (ale 100: 0.021 and 0.015 s, and 0.016 s in `audit_breakdown.md`), so the second run's 1.88 on ale 100 is the representative figure. Before the structure-keeping linear solve the temperature form stood at 1.18, 1.44, 2.09, 2.78, 1.91 and 2.32 (the M9 measurement, load 20 to 24).

On the original PMMA decks at the default tolerances the temperature form takes 0.58, 0.80, 1.22 and 0.69× the old core's time (eulerian 25, 100; ale 25, 100), and the enthalpy form 1.67, 2.36, 2.44 and 1.64×.

The item passes on the Eulerian and Darcy decks in both runs, sits at the threshold on ale 25 (1.48 and 1.51) and fails on ale 100 (1.88). What closed most of the gap is the linear solve (C31, C32; D46): the moving mesh's algebraic rows make the mass matrix a `Diagonal` with zeros, OrdinaryDiffEq then forms `W = J − M/(γΔt)` by a sparse broadcast that keeps only the entries computing to nonzero, so `W`'s structure changed from step to step (2,414 to 3,623 of `J`'s 4,780 stored entries on ale 100) and KLU analysed it afresh at every step (90 analyses in 90 steps); `StructuredKLU` scatters each `W` into the prototype's fixed structure, so one analysis serves the solve, with the same steps and temperatures within 6e-11 K (7e-6 K on ale_darcy 25, where the step sequence shifts). In `measurements/audit_breakdown.md` the share of the new core's solve outside the Jacobians, right-hand sides, factorizations and solves fell from 0.55 to 0.63 to 0.37 to 0.44 on the moving-mesh decks. The old core's moving mesh has no algebraic rows, so its iteration matrix kept its structure throughout.

What remains is common to every deck. At equal tolerance the two cores now take about the same time (0.75 to 1.19), but `Rodas5P`'s error is 3 to 22× the old default's on five of the six decks (the old core's `KenCarp4` with the smoothed estimate; on ale_darcy 25 it is half), so matching it needs reltol 1e-7 to 3e-7 and 1.3 to 1.8× the old core's steps; whether that gap is the estimator's or the residual's is open (`KenCarp4` with the smoothed estimate on the new core's residual does not reach the old default's error on two decks either, `measurements/audit_breakdown.md`). And each step forms and factors a new Jacobian, where the old core reuses its (approximate, on the moving mesh) Jacobian over many steps: the colored Jacobians are a fifth to a quarter of the new core's solve.

## 3. The manufactured-solution orders

Table: `measurements/mms_orders.md` (`mms_orders.jl`).

| case | old core | new core (both forms) |
|---|---|---|
| V-MMS-01 (20, 40, 80 cells) | 1.98, 1.99 | 1.98, 2.00 |
| V-MMS-03 (16, 32, 64) | 2.21, 2.22 | 3.27, 1.70 |
| V-MMS-04 (16, 32, 64) | 2.00, 2.01 | 2.00, 2.00 |

V-MMS-02 is an exactness case: its free-stream drifts hold in both forms (gas 3.6e-14 and 5.0e-14 against the pin 1e-12). Every quantity of the four cases passes on both cores. V-MMS-01's temporal sweep (80 cells, fixed steps 0.5 to 0.0625 s, against an eighth of the smallest step): the old core's KenCarp4 3.99, 4.01, 3.91 in L2; the new core's `Rodas5P` 3.92, 4.28, 4.43 (temperature) and 4.45, 4.03, 3.78 (enthalpy).

V-MMS-03 on the new core passes its check (each order above 1.5, the mean above 1.8) and its error is 5 to 10× smaller than the old core's on every mesh, but its temperature does not reach second order: refined to 128 cells the orders are 3.27, 1.70, 1.63 (reproduced in the review; a 256-cell value of 1.82 was a session number). The species do: the gas converges at 1.99, 2.00, 2.00 in both norms and the uniform solid is exact, tighter time tolerances change nothing, and the temperature error at 64 cells is a smooth global profile, so the van Albada reconstruction's clipping at the manufactured extrema, which would show in the species first, is not the cause. The cause is the upwind-biased face temperature of the relative energy flux, a pre-asymptotic term: to 512 cells the temperature orders rise 1.63, 1.83, 1.92, and a centred face temperature gives 1.98 to 2.00 at 4 to 6× the error (C33). Nothing is wrong in the formula; the choice between the two face temperatures is the owner's, and it matters once the correction `c` makes the relative fluxes production paths (Phase 2).

## 4. Radial bandwidth and colors

`test/newcore/jacobian/test_pattern.jl` (254/254): the radial patterns, bandwidths and color counts equal the planar ones at four cell counts.

## 5. The suite with `H` as the state

`NEW_CORE_CASES` passes in the inventory form, 188 quantities and 11 convergence checks, the same as the temperature form. Its cost: on the V-HC cases 2 to 12× the temperature form's steps (V-HC-02a 665 against 56) and 9 to 60× its wall time; on V-BO-02 3,163 against 2,354 steps and 28 s against 3.7 s at M9, and with the fixed-mesh vent 2,456 against 1,706 steps and 7.8 s against 1.6 s (`case_accuracy.md`); on the audit decks 2.9 to 6.3× the old core at matched accuracy against 1.2 to 2.8× for the temperature form (item 2).

## 6. Constraint residuals

Table: `measurements/case_accuracy.md` (`case_accuracy.jl`). Every case of `NEW_CORE_CASES` at its CI resolution and solve keywords, in both forms, records the largest scaled constraint residual of its accepted states (`statistics.largest_constraint_residual`; one or below passes): every solve passes. The largest are on the deactivating cases (V-CPL-02 0.24 and 0.42 in the temperature and enthalpy forms, V-CPL-02b 0.05 and 0.08), the moving-mesh MMS cases next (V-MMS-02 1.8e-3 and 2.0e-3), and a fixed mesh without in-depth radiation has no algebraic rows (zero); re-measured at `603e863`, unchanged. V-RAD-01's zero is exact, not a missed row: its deck is inert (one solid, no reactions or gas, a constant flux and absorption coefficient), so its intensity rows do not depend on the temperature, every stage's intensity increment is zero and the accepted states keep the consistent intensities bit for bit; on a reacting translucent slab the statistic reads 5.3e-4. The domain test rejects a step whose residual exceeds one (C5), so these are the accepted states' margins.

## 7. The energy balance against the ledger

Table: `measurements/energy_ledger.md` (`energy_ledger.jl`), on V-CPL-06's char-layer column (80 cells, the window from 120 s), at the case's tolerance (abstol 1e-8, reltol 1e-6), |ΔE − ∫ ledger| over the integrated boundary heat with the same trapezoidal quadrature over the accepted steps on both cores:

| deck | old core | new, temperature form | new, enthalpy form |
|---|---|---|---|
| Eulerian | 2.27e-5 | 5.9e-7 | 1.03e-6 |
| moving mesh | 4.69e-5 | 3.08e-6 | 4.32e-6 |

The temperature form's closure is 15 to 38× smaller than the old core's. With the ledger integral carried as a state, so the integrator's quadrature accounts for it (20 cells, fixed steps 2, 1, 0.5 s over 240 s, against a ledger of 1.24e6 J), the temperature form's defect is 2.97e-4, 1.70e-5, 8.7e-7 J on the Eulerian column (orders 4.12, 4.29) and 2.97e-4, 1.71e-5, 9.1e-7 J on the moving one (4.12, 4.22), the integrator's order; the inventory form's is 4.1e-9 and 3.3e-9 J, round-off. On the inventory form the trapezoidal closure above is the quadrature's error, not the form's.

## 8. Achieved error against requested tolerance

Table: `measurements/case_accuracy.md`. Each case at its tolerance against itself at both tolerances divided by 100: the largest temperature difference over the saved times and the cells live in both, over reltol times the largest temperature, is within budget at 10. Re-measured at `603e863`: 60 of the 66 solves are, most by orders of magnitude (the V-HC and V-RX cases at 1e-8 to 0.05, V-CPL-01 0.96, V-MMS-03 1.69 the largest). Outside it, the six solves of the deactivating cases:

- V-CPL-02, V-CPL-02b and V-BO-02, both forms. A cell is deactivated at the end of the step on which it crossed the threshold (a discrete callback, plan §12.2), so an event's time is resolved to the step, not to the tolerance: between the two tolerances the event times move by 5 to 12 % of the interval between events, and the temperatures follow them through whole park intervals (V-BO-02 26 K just after the last park, a median of 0.85 K between events over cells live in both; V-CPL-02b 0.13 and 0.61 K, V-CPL-02 0.04 K, which alone is within the temperature budget). Before its first event V-BO-02's temperature error is 0.065 K, which the event timing does not explain; measured on `phase1-review` (temperature form, CI resolution, t < 69 s against the same case at tolerances divided by 100): it lies in the heated PMMA layer (cells 155 to 163, about 684 K, at t = 36.5 s) and follows `reltol` alone, 0.065, 0.0072 and 0.0004 K at reltol 1e-6, 1e-7 and 1e-8, while dividing the temperature rows', the mass rows' or the energy rows' absolute tolerance by 100 leaves it at 0.06 to 0.10 K. It is the temperature rows' time error at the requested relative tolerance, not a tolerance-vector scale. The old core on the same deck errs by 0.41 K over t ≤ 60 s at its default (0.20 K at reltol 1e-7), six times the new core's. Nothing is changed for it; both cores exceed the budget unit on this heated layer. Their quantities of interest pass at their tolerances.

V-CPL-05, in both forms, is within budget since `case_accuracy.jl` compares over the saved times both solves reach (review X-2): its depletion limiter is disabled by design (C3), so no reference completes at a tenth or a hundredth of its tolerances (`DtLessThanMin` at the depletion point); over the 72 times to 720 s the error is 1.45e-6 and 2.86e-6 K, 0.22 and 0.43 of the budget. Locating the crossing (a continuous callback on the condensed mass) would tie the event times to the tolerance; it belongs with Phase 2's event work.

## 9. The finite slab

Table: `measurements/finite_slab.md` (`finite_slab.jl`), 320 cells, largest L2 over the output times [K]:

| case | reltol | old default | spatial floor | `Rodas5P` (default) | passes |
|---|---|---|---|---|---|
| V-HC-03 | 1e-6 | 5.89e-3 | 9.27e-5 | 9.27e-5 | yes |
| V-HC-03 | 1e-8 | 1.07e-4 | 9.27e-5 | 9.27e-5 | yes |
| V-HC-02c | 1e-6 | 2.20e-4 | 5.932e-4 | 5.933e-4 | yes, at the floor |
| V-HC-02c | 1e-8 | 5.928e-4 | 5.932e-4 | 5.932e-4 | yes, at the floor |

The floor is the default at reltol 1e-11, abstol 1e-13. On V-HC-02c the old default's L2 at reltol 1e-6 lies below the floor, its time error cancelling part of the spatial error (C28); `Rodas5P` is within 0.01 % of the floor at both tolerances. KenCarp4 with the raw estimate, the former provisional default, is at 5.93e-2 there on both cores.

## 10. Leftovers 1 and 2

Leftover 1 (`measurements/leftover1_jacobian_cost.md`, M7) and leftover 2 (`measurements/leftover2_integrators.md`, M8) are recorded with their tables; plan §12.6 has the numbers. Leftover 1's affordability decision (AD alone, or AD with analytical fast paths for conduction and gas transport, §7.1) was left for this item's solve times, which item 2 now gives.

## The provisional choices of plan §10

These are the owner's; the measurements and a recommendation:

- **Energy state.** Temperature, as the default. The inventory form passes everything the temperature form does, but at 2 to 12× its steps on the V-HC cases and 2.9 to 6.3× the old core's time at matched accuracy on the audit decks (the temperature form 1.2 to 2.8×). Its advantage, the exact ledger, shows only with the integrator's own quadrature, where the temperature form's defect converges at fourth order (item 7); with the trapezoidal ledger the two forms close alike (5.9e-7 and 1.0e-6 Eulerian, 3.1e-6 and 4.3e-6 moving). The inventory form stays the alternative, as the relative-flux work of Phase 2 plans.
- **Integrator.** Decided: `Rodas5P`, solver-chosen restart.
- **Lifting against bordering.** Lifting. The pattern is not what costs: 9 to 15 colors, below the lifted model's 15 to 17 (D23), and the colored Jacobian 13 to 25 % of a matched-accuracy solve. On a moving mesh the cost that failed item 2 was the iteration matrix's structure (C31). The eliminated system of the bordered form has no algebraic rows, so it would have escaped C31 as the old core does, at the price of every coupling of the eliminated rows enumerated by hand; the structure-keeping linear solve (C32, now the default) removes the same cost with the lifted form kept.
- **The correction `c`.** Phase 2, as planned. Phase 1 runs at `c = 0`; the relative-flux machinery it will drive is exact under a prescribed mesh velocity (isothermal advection, the free stream of V-MMS-02) and converges below second order on V-MMS-03 (item 3), which Phase 2 should look at first.
- **Analytical fast paths.** Not yet. The Jacobian is 13 to 24 % of a matched-accuracy solve, so fast paths for conduction and gas transport could take at most that from it. The order of the levers by their measured size: the structure-keeping `W` (C31, C32), now the default and measured in item 2 (the moving-mesh and Darcy decks from 1.9 to 2.8× down to 1.06 to 1.88×), then the integrator's accuracy per tolerance (the new core needs reltol 1e-7 where the old core's smoothed KenCarp4 meets the error at 1e-6), then reuse of the Jacobian across steps (a Rosenbrock-W method, an integrator choice), and only then the Jacobian's own cost. Item 2 still fails on ale 100 and sits at the threshold on ale 25; the next lever is the integrator's accuracy per tolerance, not the Jacobian.
- **`ε_m`.** 1e-2 stays for Phase 1; Phase 2 tunes it on the V-BO ladder, which C29 also needs.
- **Tolerance-vector scales.** As listed (§12.5). Item 8 finds no case whose per-row tolerances let the temperature error leave its budget; the two exceptions are the event timing and a case without a limiter, which no scale changes.

## What the gate leaves open

- Item 2 on ale 100 (1.88×) and ale 25 (at the threshold), and the affordability decision with it; the next lever is the integrator's accuracy per tolerance.
- C29: interior crossings, the V-BO ladder.
- V-MMS-03's temperature order on the relative-flux path (C33; the upwind-biased face temperature kept for Phase 1).
- Radial Beer–Lambert's physics (C37; refused until a radial case needs it).
- The surface Newton's clamp, to be replaced by the trial clamp in Phase 2 (C36).
- The design's §3.7 color count at `c = 0` (C25) and the masked-row wording (C26), to correct at the next revision of the design.
