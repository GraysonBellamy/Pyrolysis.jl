# Mesh framework: Phase 2 implementation plan

**Scope.** The implementation plan for Phase 2 of the mesh-framework design (`docs/Mesh_Framework_Design_2026-10-17.md` §6, "Burnout"), written against branch `mesh-framework-phase1` at `bb2860b` (the merge of Phase 1, its review and the review's followups) and verified against that code. It takes the shape of the Phase 1 plan (`docs/Mesh_Framework_Phase1_Plan.md`): what was verified, the module layout, the contracts, the kernel and solver plans, the tests and cases, the milestones with their gates, the risks, the branch rules, the gated choices and the owner's questions. Every number a decision rests on was measured on this tree during planning (§1.5 and Appendix A, with the prototype's standing against each gate item stated in §7.3: the burnout cases on both cores, the old core's burnout wall times, a scratch prototype of the correction on V-BO-01, event location by a continuous callback, the restart policy after an event, the integrator candidates' accuracy per tolerance, and the Douglas fir deck on both cores). The prototype is planning evidence, not package code; it is committed beside this plan so its tables can be reproduced (§9.2).

**Conventions.** `file:line` refers to `bb2860b`. "§3.x" without a document name is a section of the design document; "plan §x" is the Phase 1 plan; "gate §x" its gate document (`docs/Mesh_Framework_Phase1_Gate.md`). D1 to D52 and C1 to C40 are the Phase 1 plan's departures and conflicts; this plan's conflicts continue the series from C41, and its implementation record will continue the departures from D53. Case ids (V-BO-01, …) are the verification cases under `test/verification/cases/`. "The prototype" is `scripts/mesh_phase2_prototype/` (§9.2).

**Contents.** 1 What was verified · 2 Module layout · 3 Data structures and contracts · 4 Kernel-by-kernel plan · 5 Solver plan · 6 Test and verification plan · 7 Milestones, gates, estimates · 8 Risk register · 9 Branch and process · 10 Decisions carried as gated choices · 11 Summary · 12 Implementation status · 13 Owner's questions · Appendix A Measurements · Appendix B What the reviews changed

---

## 1. What was verified

### 1.1 What Phase 1 built, and what it did not

The Phase 2 planning brief (the owner's, kept outside the repository) lists deactivation as built and the correction `c` as not built. Both hold, with these precisions.

**Deactivation exists as the terminal event, from the exposed side only.** `src/Deactivation/events.jl:174-192` (`deactivate!`) deactivates the outermost live cell while it is below `ε_m` of its own initial condensed mass, re-testing the receiver after each transfer, and returns `:interior` as soon as any other live cell is below the threshold (`:187`); the solve then ends `:Terminated` (D27, `src/Integration/solve.jl:378-420`). The transfer (`:242-326`) is the conservative volume merge of §3.5 with the pore-preserving rule, on a fixed mesh the condensed inventories only with the gas vented (D47), the regularization ledgered (D48), and the dead rows masked (`_mask!`, `:357-361`). Every crossing at a step is handled in one `DiscreteCallback` at the end of the step; the restart step follows the policy (`rodas5p_policy` is solver-chosen, `src/Integration/policy.jl:62-74`). There is no cascade rule beyond "the outermost cell while below": a block of cells that cross together deactivates one cell, whose receiver then rises above the threshold, and the rest end the solve as `:interior`. That is what ends V-BO-04 after one event (§1.5).

**The relative-flux kernel exists, under a prescribed mesh velocity only.** `src/Kernels/relative_flux.jl` reconstructs every face state with the van Albada slope and carries the enthalpy of the same mass fluxes at the reconstructed face temperature, but its relative velocity is `u = −w·n` (`:180`): it assumes a solid at rest (`:37`), and `setup` installs it only under `w_override` (`src/Integration/setup.jl:350-361`), refuses an override on a reacting material (`:76-80`) and refuses deactivation under an override (`:81-84`). The kernel is therefore a verification path (V-MMS-02/03/04, isothermal advection, free stream, history transport), not a production path.

**The correction `c` does not exist.** No monitor, floor, activation switch or response time `τ` is in `src/` (a search of `Closures`, `Kernels`, `States` and `Integration` for a monitor, an activation or an equidistribution finds nothing); the balance rows carry no sweep term `s_k` (`src/Closures/balance.jl:139-148`), the consistent sweep none (`:178-194`), and the Jacobian pattern's velocity rows reach their own cell only (`src/ColoredJacobian/pattern.jl:96-115`). Plan §10's "soft activation and relative fluxes implemented, inactive" overstates what landed: the relative fluxes exist, the activation does not. The design's §3.3 and §6 Phase 2 items 2 and 3 are the whole of what this plan builds.

**The pieces the correction needs are in place.** `MeshTopology{1}` carries a `lagrangian` node set (`src/Meshes/topology.jl:109-119`; `line_mesh(…; lagrangian = …)`, `src/Meshes/generators.jl:7-29`) that nothing sets from a problem; `CellReference` holds the initial node positions (`src/States/derived.jl:21-29`), so the initial thickness of every cell is available for the activation; `DerivedState` holds `x`, `w` and `active` (`:145-172`); the pattern is built from connectivity with a `reconstruction` flag (`pattern.jl:52-58`) and verified by containment; the domain test scales a balance row's residual by `|θ_k| V_k + |θ_A| V_k + Σ_f A_f |w_f|` (`src/Integration/checks.jl:203-210`), which the sweep term joins.

### 1.2 The burnout cases on both cores (re-run)

Every number the planning brief quotes about the V-BO cases was re-run on `bb2860b` (Appendix A.1, A.2; old core and new core in one session each, warm medians of three, load 11 of 22 CPUs):

| case | old core | new core (temperature form, Phase 1 as built) |
|---|---|---|
| V-BO-01 (moving, in-depth) | `Success`, 0.80 s, 853 steps, 60 parks, every QoI | `:interior` at 241.6 s before any event (cells 179 and 180 below the threshold under live cells 181 and 182); 0.45 s |
| V-BO-02 (fixed, surface) | `Success`, 1.29 s, 1,572 steps, 60 parks | `Success`, 1.32 s, 1,706 steps (953 rejected), 60 events, every QoI (ledger 1.5e-3 against its 3e-3, C39) |
| V-BO-03 (moving, transparent substrate) | `Success`, 0.58 s, 681 steps, 48 parks | `:interior` at 241.6 s, 0 events |
| V-BO-04 (uniform zero-order) | `Success`, 0.034 s, 353 steps, 35 parks | `:interior` at 468.3 s after one event (the receiver rose above the threshold, 38 cells stayed below) |
| V-BO-05 (fixed, in-depth) | `Terminated` at 337.0 s as pinned, 1.42 s | `:interior` at 263.9 s, 0.34 s (the mass-fraction threshold against the pinned matrix-volume-fraction criterion, C9) |

C29 holds exactly as recorded. The old core's wall times at the report ladder (Appendix A.2; 30, 60 and 120 PMMA cells): V-BO-01 0.25, 0.58, 2.11 s; V-BO-02 0.45, 1.11, 3.70 s; V-BO-03 0.21, 0.54, 1.79 s. The same configuration's time moves by up to 40 % between runs on this machine (V-BO-01 at 60 cells: 0.58 and 0.80 s in two sessions), so every comparison below is of the two cores back to back in one process, and the gate's ratios are judged against that spread.

### 1.3 Pointer audit

Every `file:line` and library claim the design's Phase 2 text, plan §12.7 and the planning brief rely on, checked against `bb2860b`.

| Claim | Location checked | Verdict |
|---|---|---|
| `deactivate!` handles the outermost live cell only and ends `:interior` otherwise (C29, D27) | `events.jl:174-192` | holds |
| Conservative transfer by volume merge, pore-preserving; fixed mesh transfers condensed inventories and vents the gas (D47) | `events.jl:242-326` (`_transfer!`, `_merge_volume!`) | holds |
| Dead rows are identity rows, the pattern fixed (C26) | `src/Assembly/residual.jl:227-236` (`_dead_rows!`); `pattern.jl` docstring | holds |
| The reactivity switch is one above 1.5 `ε_m` (C27) | `src/Kernels/reactions.jl:12` (`SWITCH_WIDTH`), `:40-52` | holds; `C¹` smoothstep |
| Relative fluxes active only under `w_override`; the solid at rest | `relative_flux.jl:37, :180`; `setup.jl:350-361, :76-84` | holds |
| Van Albada reconstruction, upwind-biased face temperature (C33) | `relative_flux.jl:9-17, :150-160, :183-185` | holds |
| Pattern with reconstruction: cells `i−2..i+2`, positions `i−2..i+3`, the cell's two velocities; velocity rows reach their own cell | `pattern.jl:52-58, :70-77, :96-115` | holds; the velocity rows must widen under the correction (C43) |
| Constraint scale of a balance row omits any sweep term | `checks.jl:203-210` | holds (C52) |
| The deactivation callback: one `DiscreteCallback`, the consistent sweep after the transfer, admissibility, the policy's restart | `solve.jl:378-420` | holds |
| The integrator raises the derivative-discontinuity flag before each discrete callback's `affect!` and folds every callback's verdict, so a policy's flag rule runs inside the wrapped `affect!` (D40) | `solve.jl:301-324` (`_observing`) | holds |
| `min_thickness` ends a moving mesh's solve by a `ContinuousCallback` (D18) | `solve.jl:350-364` | holds |
| Event times follow the step, not the tolerance (gate item 8) | `solve.jl:378-420`; measured again in A.5 (0.88 s shift between reltol 1e-6 and 1e-8, 0.133 of the median event interval on V-BO-02) | holds |
| `MeshTopology.lagrangian` exists and nothing sets it from a problem | `topology.jl:109-119`; `generators.jl:29`; `setup.jl` (no `lagrangian` keyword) | holds |
| The surface Newton clamps its iterate and warns (C36) | `src/BoundaryConditions/surface_balance.jl:174, :227-232, :241-256` | holds |
| The Darcy advection's first-order upwind switch (C21) | `src/Kernels/gas_transport.jl:105, :126` | holds |
| `NEW_CORE_CASES` routes V-BO-02 and not V-BO-01, 03, 04, 05 | `test/verification/casedefs.jl:161-184` | holds |
| V-BO-05's pinned stop window 330 to 345 s | `test/verification/cases/V-BO-05.jl:80-82` | holds; C45 |
| The V-BO step recorder reads the new core through the layout (C8) | `V-BO-01.jl` `_vbo_core_solve`, `_vbo_recorder` | holds |
| `StructuredKLU` is the default linear solve on a layout with algebraic rows (D46) | `src/Integration/linear_solve.jl:97` (`default_linear_solver`) | holds |
| Tabulated sensible enthalpy at a non-finite temperature | `src/States/caloric.jl:112-116` (`floor(Int, …)`) | throws `InexactError` instead of returning `NaN` (C47; hit by the prototype) |
| The Douglas fir deck (14 components with liquid moisture, a state-dependent heat of sorption, Darcy–Fick, ALE, surface absorption) runs on the new core | Appendix A.8 | holds: 49 colors, surface temperature within 0.052 K of the old core over 60 s |

Library facts, in the pinned versions of the repository `Manifest.toml`: `ContinuousCallback` with `affect_neg!` and `rootfind = RightRootFind` locates a downward crossing and applies the affect at the state just past the root (SciMLBase 3.56.1; exercised in A.5); the Rosenbrock methods `Rodas5Pe`, `Rodas5Pr`, `Rodas4P2`, `Rodas4P`, `Rodas3P`, `ROS34PW2` and `Rodas23W` are reachable through the umbrella's `OrdinaryDiffEqRosenbrock` (exercised in A.7; the umbrella re-exports `Rodas5P` and `Rosenbrock23` only).

### 1.4 Conflicts between the design, the planning brief and the code, with proposed resolutions

Each entry records both locations; none is resolved silently, and the owner's answer changes the scope where marked.

**C41. The relative velocity is `v_s − w = −c`, not `−w`, once the correction exists.** `relative_flux.jl:180` reads the node velocity and assumes the solid at rest, which `setup` enforces by refusing reactions under an override. With `w = v_s + c` the inventories cross a face at `(v_s − w)·n = −c·n`, so the kernel must read the correction, not the mesh velocity, and the balance rows (which produce `v_s` implicitly through `w`) must carry the correction's sweep. *Resolution:* the derived state gains the per-node correction `c` (§3.2), filled by the closure's correction object before the kernels run; the kernel reads `ds.c`; under `w_override` the correction object is the prescribed velocity itself with `v_s = 0`, so the verification path is the same code (§3.3). The prototype implements exactly this (Appendix A.3) and its pattern is contained in the dense Jacobian.

**C42. The cascade rule.** Design §3.5 says every cell below the threshold at a step deactivates in the same callback; `deactivate!` deactivates the outermost cell while it is below and stops once the receiver, holding the dead cell's inventory, is above. A block of cells crossing together therefore ends the solve (`:interior`, V-BO-04). *Resolution:* the block rule of §3.4: the contiguous run of cells from the exposed side that lie within a near-threshold band `κ ε_m` and contain at least one cell below `ε_m` merges, cell by cell, into the first live cell beneath the run (or into the lowest cell of the column, which survives), in one invocation. Measured on V-BO-01 (A.3): with `κ = 1` the band leaves cells at 1.3 to 3.6 % of their mass under a vanished top and the solve ends `:interior`; with `κ = 10` the whole PMMA layer ends in one invocation and every QoI holds. `κ` is a gated parameter (§10).

**C43. The Jacobian pattern under the correction.** `pattern.jl:96-115` gives a balance row its node, its cell's other node and its cell's masses and energy. The sweep term `s_k = Σ_f A_f (c·n)_f` reads the correction at the cell's two nodes, and the correction at a node reads the positions of the three nodes around it and the masses of the two cells sharing it (the monitor), so a balance row reaches the positions of nodes `k−1..k+2` and the masses of cells `k−1..k+1` (the energy only if the optional reaction term is used). The cell rows already reach that far through the reconstruction. *Resolution:* the pattern widens for the velocity rows when a correction is present (§3.6); measured at 30 colors on V-BO-01 (in-depth radiation, three components) against 18 at `c = 0`, contained in the dense Jacobian (A.3).

**C44. The design's activation is per cell; the correction is per node.** §3.3 defines the switch in `log(h_i/h_i(0))` per cell and says nothing about the node between a thin cell and a healthy one. A product of the two cells' switches leaves such a node still, so a single thin cell between two healthy ones is never fed and crosses the threshold alone, which is V-BO-01's failure (a sub-surface cell under the absorption peak). *Resolution:* the node switch is the smooth "or" of its cells', `a_n = 1 − (1 − a_i)(1 − a_{i+1})`; the node then moves whenever either cell is thin, feeds the thin cell from the healthy one, and the healthy cell's thinning activates the next node, so the activated region spreads through the vanishing layer (A.3; the "and" form is measured in A.4 for the record).

**C45. V-BO-05's pinned stop.** The case pins the interior-void stop at 330 to 345 s under the old core's criterion (a matrix volume fraction below 1e-4 beneath a surviving skin); the new core ends `:interior` at 263.9 s when a cell's condensed mass falls below `ε_m` of its initial value, with no warning and no "park events" counterpart. There is no exact mapping (C9). *Resolution:* V-BO-05 is re-pinned for the new core at its own criterion (the stop time at the first `ε_m` crossing beneath a live cell on a fixed mesh, with a tolerance set from a ladder), the warning QoI replaced by the outcome `:interior`; on a fixed mesh no correction acts (every node pinned), so the case's purpose, a documented stop, is unchanged. Owner question 8.

**C46. V-BO-04's QoIs are the old core's.** "Park cascade: 40 → 5 live cells" and the exact thickness law through pairwise cascades describe the old core's one-cell-per-invocation parking at a 5 % thickness threshold. Under the block rule all 40 cells cross `ε_m` together and merge into the lowest cell in one invocation; the survivor holds `40 ε_m` of a cell's reference mass, so the exact law `L(t) = L₀ exp(−r t/ρ_s)` continues to hold (its thickness is `40 ε_m h₀ = ε_m L₀` at the crossing, which is the law's value) until it crosses alone and the solve ends `:burned_out`. *Resolution:* re-pin the live-cell QoI to one invocation and one survivor, keep the exact-law, concentration and isothermal QoIs, and pin the final outcome; the one-callback result of Phase 0 item 7 (three cascades) is the old core's and is not reproduced. Owner question 8.

**C47. A non-finite face temperature throws.** `_tabulated` (`caloric.jl:112-116`) takes `floor(Int, T)` of the trial temperature; at `NaN` this throws `InexactError` instead of returning `NaN`, so a wild trial state in the relative-flux kernel (which evaluates the enthalpies at a reconstructed face temperature) aborts the solve where the design's trial-state contract (§3.2, D33) wants a `NaN` residual and a rejected step. Hit once by the prototype at 120 cells (A.4). *Resolution:* `_tabulated` and `enthalpy_slope` (`caloric.jl:140`, the same `floor(Int, …)`, reached through the caloric inversion) return `NaN` for a non-finite argument (two lines, trial-state contract); a test drives the relative-flux kernel with a `NaN` face temperature.

**C48. "Surface node moves with the physical surface velocity" against equidistribution.** §3.3 keeps `c = 0` at the surface node. The top cell is then fed only through its lower node, so the top cell of an activated region is thinner than the equidistributed value and crosses the threshold first; the band of C42 absorbs the difference (its measured width is 1.3 to 9.2 % of a cell's mass at the two invocations of the `τ = 1` run, A.3). Kept as designed: the surface is the physical surface, and the band is the cheaper remedy.

**C49. The wall-time clause of the Phase 2 gate inherits Phase 1's item 2.** The gate requires wall time no worse than the old core's burnout runs. The old core's V-BO-01 at 60 PMMA cells takes 0.58 to 0.80 s; the prototype takes 1.6 to 2.3 s (`τ = 1`, two invocations) or 2.8 to 3.5 s (`τ = 0.2`, one invocation) at the same tolerance (A.3), 3.5 to 6× and measured in another process than the old core's runs (A.2), and Phase 1's moving-mesh decks already stand at 1.5 to 1.9× at matched accuracy (gate item 2, C40). Of the prototype's cost about half is the colored Jacobian `Rodas5P` builds every step (30 colors, 1,091 steps against the old core's 164 Jacobians over 853 steps) and a fifth the rejected steps (246 of 1,337 proposals, 52 of them by the domain test). *Resolution:* the clause is kept as written and measured at matched accuracy on the V-BO ladder; the levers are listed with their measured size (§5.5) and the exception of C40 is extended to it only by the owner (question 6). No integrator candidate closes Phase 1's accuracy-per-tolerance gap on the moving-mesh deck (A.7).

**C50. The energy state for the relative-flux work.** The design chose the inventory form for the work where `c ≠ 0`; the gate chose the temperature form as the default (C40) because the inventory form costs 2 to 13×. *Resolution:* the V-BO ladder runs in both forms (the prototype's enthalpy run is in A.4); the temperature form stays the default unless its ledger closure on the burnout cases exceeds the old core's, which is a gate item (§7.3 item 3).

**C51. The reactivity switch under a synchronized layer.** With the correction active the whole layer approaches the threshold together, so the switch (C27: `s = x²(3 − 2x)`, `x = m/(1.5 ε_m m_ref)`, 0.74 at the threshold) throttles every cell of the layer at once; on the old core it throttled one cell at a time. The band then deactivates cells inside the switch's range, not above it: the synchronized runs merge cells at 1.0 to 1.75 % of their mass (30 cells), 1.0 to 1.2 % (`τ_p = 0.1`) and 0.3 to 0.6 % at `ε_m = 3e-3`. The receiver's switch and deactivation reference are its own initial condensed mass, so a substrate receiver has `s = 1` and the deposited residue gasifies unthrottled. *Resolution:* the switch's fate under a correction is decided by R25's decomposition measurement (§7.3 item 2); the structural alternatives are no switch under an `EquidistributionCorrection` (the band pre-empts the stiff approach the switch was for, since cells leave at or above `ε_m`) or a band trigger at the switch's onset; question 12.

**C52. The constraint scale omits the sweep.** `checks.jl:203-210` scales a balance residual by `|θ_k| V_k + |θ_A| V_k + Σ_f A_f |w_f|`; with the correction the row's own scale gains `Σ_f A_f |c_f|`, and without it the domain test rejects steps spuriously where the correction is large (the prototype's 52 domain rejections at `τ = 0.2` with the largest accepted residual at 0.69 of the bound, A.3). *Resolution:* add the term (§3.7).

**C53. Events and invocations.** `DeactivationEvent` is one transfer; the design's "one event per layer" counts callback invocations (one restart, one consistent sweep). *Resolution:* every event records the invocation it belongs to, the outputs and the gate count invocations, `park_times` keeps one entry per cell as the old core's.

**C54. The interior policy on a moving mesh at `c = 0`.** The design's fallback "the core ships with `c = 0` and deactivation" leaves V-BO-01 and V-BO-03 at `:interior` (A.1), which the old core handles by parking the interior cell and shifting the live cells down one slot (`src/Adaptivity/cell_park.jl:90-143`). A slot-shifting interior merge on the new core is possible in 1D (the state is permuted inside the callback, the pattern never changes) but every per-slot quantity would have to shift with it: the six per-cell vectors of `CellReference` (`derived.jl:21-29`), the per-row tolerance vector the integrator holds (`tolerances.jl`), `CellDeactivation.reference` and the reactivity switch's reference, the relative-flux kernel's per-cell slope scales, the Lagrangian node indices, the layout mask, `replay_events!` and `park_times`, and the integrator's interpolant after a state permutation (which also touches located events, §3.5). The planning brief asked for a scratch measurement of this policy on V-BO-01 (events per layer with interior cells merged); it was not run, so whether V-BO-01 and V-BO-03 pass their QoIs on the new core at `c = 0` with it is not known, and its estimate below is unverified. *Resolution:* not built in Phase 2 as planned; the block-and-band rule plus the correction are the primary path, and the slot-shifting merge is the fallback the owner schedules only if the correction fails its gate (question 1; about one and a half weeks, unverified, with the contracts above and the tests of §6.1 item 17).

**C55. Thresholds as material parameters (design §7 question 9).** The activation thresholds of §3.3 (30 % and 10 % of the initial thickness) must sit below a layer's design shrinkage; a char former that keeps 20 % of its thickness must not activate. Nothing in `Material` declares a residual thickness. *Resolution:* `setup` takes the activation thresholds and the band width as keywords with package defaults, and `Material` gains an optional residual-thickness fraction (additive, old core untouched, fingerprint-checked) from which the defaults are derived; question 4.

**C56. The gate's "charring decks identical to Phase 1 (same code path with `c = 0`)".** When the correction object is present but never activates, the relative-flux kernel adds exact zeros and the residual is bit-identical; the Jacobian pattern is wider (C43), so the colored passes and the factorization change round-off and the trajectory agrees to round-off only. *Resolution:* two checks: a deck solved without a correction object is bit-identical to Phase 1; the same deck with the correction object present and inactive agrees to 1e-10 relative in every output (§7.3 item 5).

### 1.5 What the prototype established (summary; Appendix A.3, A.4)

The correction was prototyped as method overrides on the new core in a dedicated session (no package code changed): the node correction `c_n = a_n (M_{n} h_{n} − M_{n−1} h_{n−1})/τ` with `M = α + β·(condensed mass per volume)/ρ_ref`, the smooth "or" activation of C44 in `log(h/h₀)` between 30 % and 10 %, `c = 0` at the pinned, surface and layer-interface nodes; the balance rows and the consistent sweep with the sweep term; the relative-flux kernel at `u = −c·n`; the pattern of C43 (contained in the dense Jacobian, colored equal to dense); the block-and-band rule of C42. On V-BO-01 at its CI resolution (60 PMMA cells, in-depth radiation, temperature form, the case's tolerances):

| configuration | invocations (cells) | outcome | peak MLR [g m⁻² s⁻¹] (ref 32.39) | t_peak [s] (335) | t_bo [s] (363.5) | mass loss [%] (exact 69.915; old core 69.921) | wall [s] (old 0.58 to 0.80, another process) |
|---|---|---|---|---|---|---|---|
| `c = 0`, block rule | 0 | `:interior` at 241.6 s | 23.58 (−27 %) | 241.6 | — | 35.5 | 0.40 |
| `τ = 1`, `κ = 1` | 3 (12) | `:interior` at 355.9 s | 32.371 | 334.5 | — | 67.7 | 1.26 |
| `τ = 0.2`, `κ = 1` | 1 (2) | `:interior` at 359.5 s | 32.391 | 335.0 | — | 68.7 | 2.03 |
| `τ = 0.2`, `κ = 3` | 1 (41) | `:interior` at 362.4 s | 32.391 | 335.0 | — | 69.4 | 2.35 to 2.75 |
| `τ = 1`, `κ = 10` | 2 (34 + 26) | `Success` | 32.392 | 335.0 | 363.5 | 69.921 | 1.6 to 2.3 |
| `τ = 0.2`, `κ = 10` | **1 (60)** | `Success` | 32.391 | 335.0 | 363.5 | 69.921 | 2.8 to 3.5 |

With `τ = 0.2` s and `κ = 10` the PMMA layer ends in one invocation at 358.4 s (the merged cells at 1 to 3.6 % of their mass), every V-BO-01 quantity is within its tolerance (peak within 0.01 %, burnout time exact, ledger 7e-5), the smallest live cell over the run is 1.0 µm (the design's bound relative to the layer: the 60 cells shrink with it), and the cost is 3.5 to 6× the old core's at equal tolerance. The one-event-per-layer hypothesis therefore holds on the flagship case at 60 cells for two thresholds with a band of ten thresholds; at 120 cells, at `ε_m = 3e-2` and at `τ_p = 1` the layer takes two invocations, which the `N²` scaling of §3.2 explains and A.10 removes at a constant layer time; the gate's graded grid, heterogeneous layer and the rest of the matrix remain to be run (§7.3). At 30 PMMA cells the layer still ends in one invocation but the peak reads 36.7 against the old core's 32.25 at the same resolution, an end-of-layer spike that the ladder must resolve (A.4, R27).

---

## 2. Module layout

### 2.1 Rules

Unchanged from plan §2.1: the two cores stay separate modules; the old core is not edited except additively (C55's material field, if the owner takes it, and the trial-clamp keyword of C36 in the shared surface Newton, `src/BoundaryConditions/surface_balance.jl`), and any such edit re-runs `test/newcore/gates/oldcore_fingerprint.jl`, whose digests must still equal 17bc477's (C35), and the old core's `Pkg.test()`. Nothing is deprecated; a replaced rule is deleted. The guides are not touched (the new core is not user-facing until cutover). The new modules are reached as `Pyrolysis.<Module>.name` until cutover (C10).

### 2.2 What changes in each new module

| Module | Phase 2 change | Why |
|---|---|---|
| `Closures` | The correction interface (`Correction`, `NoCorrection`, `EquidistributionCorrection`, `PrescribedCorrection`), `correction!` filling the per-node `c`, the monitor and the activation, the sweep term in the balance rows and in `consistent_velocities!`, the Lagrangian-node rule (`c = 0`) | design §3.3; C41, C44 |
| `Kernels` | `RelativeFluxKernel` reads `ds.c` and is installed whenever a correction object is present; its face temperature stays upwind-biased (C33) unless the ladder says otherwise (§10); `ReactivitySwitch` gains a width field set by `setup`; `ThermalBoundaryKernel` passes the trial-clamp keyword of C36 to the shared surface Newton | C41, C51, C36 |
| `States` | `DerivedState.c` per node (zero on a fixed mesh); `CellReference.thickness` (`h_i(0)`); the tabulated enthalpy and its slope return `NaN` at a non-finite temperature (C47); no new state rows | C41, C44, C47 |
| `Deactivation` | The block-and-band rule in `deactivate!`; the invocation index on `DeactivationEvent`; the layer-end rule (receiver beneath the block, or the lowest cell survives); the fixed-mesh rule unchanged | C42, C46, C53 |
| `ColoredJacobian` | The pattern's velocity rows widen under a correction; the color-count test gains the correction configurations | C43 |
| `Integration` | `setup` keywords `correction` (default `NoCorrection`), `lagrangian`, `activation`, `band`, `switch_width`, each added to `_SOLVE_KEYWORDS` and forwarded by `_solve`; the constraint scale with the sweep (C52); event location by a continuous callback with the discrete predicate as its safety net, behind an `integrate`/`solve` keyword `located_events` (the callback is built in `integrate`, §5.2); the restart policy measured (§5.3) | §3.5, §3.7 |
| `Outputs` | `correction_velocity` per node and `min_cell_thickness` as series; events with their invocation; `invocations` as a statistic | C53, the gate |
| `Assembly` | `residual!` calls `correction!` after `derive!` and before the kernels (the kernels read `c`, the closure's rows read it again) | C41 |
| `Materials` (old core shared, additive) | Optional residual-thickness fraction (owner question 4) | C55 |
| `BoundaryConditions` (old core shared, additive) | A keyword of `solve_boundary_temperature!` selecting the smooth trial clamp in place of the hard clamp and the warnings; its default is today's clamp and warnings, which the old core (Solver, Jacobian operators, Discretization) keeps | C36 |
| `test/newcore`, `test/verification` | New contracts and tests (§6); V-BO-01, 03, 04, 05 routed with re-pinned QoIs; three new burnout cases (§6.2) | §7.3 |

### 2.3 Not changed

`Meshes`, the energy forms, the tolerance vectors, the kinetics and property layers, the Beer–Lambert and gas-transport kernels, the MMS driver (its `w_override` path becomes `PrescribedCorrection` behind the same keyword, bit-identical), the harness, the Phase 1 gate tables (`case_accuracy.jl` and `energy_ledger.jl` are extended, §6.3).

---

## 3. Data structures and contracts

Each entry: a Julia sketch (names provisional), the invariant it keeps, where the prototype exercised it, and the contract test it must keep passing.

### 3.1 The correction interface

```julia
abstract type Correction end
struct NoCorrection <: Correction end                       # c ≡ 0; the kernel is not installed; Phase 1 bit for bit
struct PrescribedCorrection{F} <: Correction; w::F; end      # the verification override: the velocity rows are w_n − w(x_n, t), c is the state's w, v_s = 0
struct EquidistributionCorrection{M} <: Correction
    monitor::M                 # M(ds, i) = α + β·(condensed mass per volume)/ρ_ref (+ γ·reaction term, optional)
    τ::Float64                 # response time [s]
    r_hi::Float64; r_lo::Float64   # activation in log(h/h₀): 0 above r_hi (0.3), 1 below r_lo (0.1), C¹ between
    lagrangian::Vector{Int}    # nodes held at c = 0 (layer interfaces), besides the pinned and surface nodes
end
correction!(c::AbstractVector, corr::Correction, ds::DerivedState, topo, reference, t)   # per node, pure in (x, m, t)
sweep(corr, ds, topo, k) = Σ_f σ_f A_f n_f c_{node(f)}                                    # the volume rate cell k gains from c
```

Invariants: `c = 0` at the pinned node, at the surface (outermost live) node, at every Lagrangian node and at every node both of whose cells are above `r_hi` of their initial thickness, so a deck that never activates is unaffected; `c` is a `C¹` function of the state (the activation is a clamped cubic in `log(h/h₀)`, `C¹` and not `C²` like the reactivity switch and the smooth floors the design already accepts under §3.10, with the quintic smoothstep as the `C²` upgrade if the ladder shows step restarts at the joins; the monitor smooth; the thickness floored only to keep the logarithm finite, since the activation is already one for every `h ≤ r_lo h₀`, negative included, and at a crossed node `c = (h_a − h_b)/τ` pushes the node back); with the balance row `Σ_f A_f (w·n)_f = (θ_k − θ_A) V_k + s_k` the cell volume obeys `dV_k/dt = θ_k V_k + s_k` semi-discretely (the `balance.jl` docstring's `dV_k/dt = θ_k V_k` no longer holds under a correction and changes), a uniform concentration crosses the faces at `−ξ s_k` because the van Albada slope of a uniform field is zero, so `d(ξV)/dt = ωV − ξ s_k` gives `V dξ/dt = (ω − ξθ)V`: the sweep cancels exactly and the free stream holds for any `c`; the row stays bidiagonal in `w` because `c` reads positions and masses only, so the one-pass consistent sweep stays exact; with `NoCorrection` nothing is installed and the residual is Phase 1's bit for bit; with `PrescribedCorrection` the velocity rows are replaced by `0 = w_n − w(x_n, t)` and `correction!` copies the state's velocity, `ds.c .= ds.w`, which is today's override path bit for bit: the algebraic rows hold only to the domain test's tolerance between sweeps, so evaluating the prescribed function into `c` instead would change every face flux by that residual and move the cell rows' dependence from the velocity columns to the position columns. Prototype: `corr_c` and the overridden `_velocity_rows!`, `_consistent_velocities!` (A.3). Tests: §6.1 items 1 to 4.

### 3.2 The node correction

For interior node `n` between cells `n−1` (inward) and `n` (outward), in 1D:

`c_n = a_n · N² · (M̃_n h_n − M̃_{n−1} h_{n−1}) / τ`, `a_n = 1 − (1 − a_{n−1})(1 − a_n)`, `a_i = S((ln r_hi − ln(h_i/h_i(0))) / (ln r_hi − ln r_lo))`, `S(x) = x²(3 − 2x)` clamped to `[0, 1]`, `M̃ = M / ⟨M⟩` the monitor normalized by its mean over the layer, and `N` the number of cells between the layer's bounding nodes (pinned, Lagrangian or surface).

The sign moves the node into the thicker (or heavier) cell, which is the relaxed equidistribution of `M h`; at exact equidistribution `M̃_i h_i` is constant over the activated region and the bounds of §3.3 (design) hold relative to the region's total. The factor `N²` is MMPDE5's `1/Δξ²` with `Δξ = 1/N`: without it (the prototype's form, `c_n = a_n (h_n − h_{n−1})/τ_p`) a fully activated region obeys `ḣ_i = (h_{i+1} − 2h_i + h_{i−1})/τ_p`, a diffusion in the cell index whose layer-scale mode relaxes in `T_sync ≈ τ_p N²/π²` and whose stiffest mode decays at `4/τ_p`, so the same `τ_p` means 18, 73 and 290 s of synchronization time at 30, 60 and 120 cells against a layer life of about 60 to 70 s once the layer is below 30 % of its thickness. That reproduces every row of A.3 and A.4 without tuning: 120 cells cannot synchronize at `τ_p = 0.2` (two invocations, the band saturated), 60 is marginal (the band of ten thresholds takes cells at 1 to 3.6 %), 30 over-synchronizes (cells at 1 to 1.75 %, R25), `τ_p = 1` gives two invocations and `τ_p = 5` none. With `N²` in the formula `τ` is the layer's equidistribution time (`T_sync = τ/π²`), the same `τ` holds across the ladder, and the prototype's runs correspond to `τ = 0.2 N²` s, that is 180, 720 and 2,880 s at 30, 60 and 120 cells; the default is set by the ladder at constant `τ` (A.10 tests the scaling), and the heuristic is `T_sync` at most about a tenth of the layer's remaining life after full activation. The price of a small `τ` is the stiffest mode (the prototype's 1,091 to 1,286 steps for `τ_p` from 0.2 to 0.1). The monitor is normalized because an unnormalized `M = α + β ρ/ρ_ref` only rescales the relaxation (`M = 2` on PMMA doubles the correction, which is the `τ/2` duplicate of A.4), so `α` and `β` are hidden `τ` scalings and only the monitor's shape over the activated region acts; on a layer of uniform density `β` changes nothing at all, and on a layer of pure cells of different densities each cell's mass fraction equals its thickness ratio, so `β` acts only where cells hold mixed compositions (V-BO-07, §6.2). The activation reads each cell's initial thickness from `CellReference`; a receiver that absorbed a dead cell's volume keeps its own reference, so its activation drops after a merge, as it should, and when the receiver is a layer cell the "or" at its lower node then redistributes the merged residue down the layer at about `h₀/τ_p`, twenty times the burning rate, which lifts the cells beneath out of the band and is why a first band merge inside the layer is followed by a second invocation (A.3's `τ_p = 1` row). The derived field `ds.c` is filled once per residual evaluation and read by the relative-flux kernel and the closure; it is differentiated by the colored AD like every other derived quantity.

### 3.3 Relative fluxes under the correction

`RelativeFluxKernel` is unchanged except that its relative velocity is `u = −c_f · n_f` (C41), so every inventory, the pore gas included, crosses faces at the correction alone; the solid velocity never appears because the balance rows carry it. The enthalpy flux stays `Σ_j F_j h_j(T_f)` with the same `F_j` as the mass rows (isothermal advection exact); the face temperature stays the upwind-biased reconstruction (C33), the smaller error in the practical range, revisited only if V-MMS-03's order matters for a burnout QoI (§10). A gas species' flux across one face then carries enthalpy at two face temperatures, the upwind-biased one for the `−c` part and the centred one for the Darcy–Fick part, each part isothermally exact. The upwind side of the reconstruction switches at `c = 0` (`relative_flux.jl:137`), a kink the flux of a uniform composition does not see but a mixed activated region does, and which gate item 9 may detect; the temperature slope scale is the solve's `abstol` over the cell thickness, effectively unregularized where the two gradients oppose (Phase 1's choice). The prototype ran without the history pair, so its transport under a correction with reactions on rests on Phase 1's prescribed-velocity tests until §6.1 item 5 runs. Under `PrescribedCorrection` the kernel reads the state's velocity through `ds.c` (§3.1), so its inputs are exactly today's and V-MMS-02/03/04 and the isothermal, free-stream and history tests run unchanged; the pattern keeps the cell rows' two velocity entries in that case. Prototype: the overridden `contribute!` (A.3). Tests: `test_relative_flux.jl` as it stands, plus §6.1 item 5 (isothermal advection with reactions off under an `EquidistributionCorrection`).

### 3.4 Deactivation: the block-and-band rule

```julia
# on CellDeactivation: band::Float64 — ε_band, the near-threshold band as an absolute mass fraction (default 0.1, ≥ ε_m)
deactivate!(u, d, def, t) -> :transferred | :burned_out | :interior
```

Rule, on a moving mesh only (a fixed mesh keeps the Phase 1 rule, D47, so V-BO-02 and V-BO-05 are unchanged by construction; the prototype applied the band regardless of the mesh and never ran it on V-BO-02): let `top` be the outermost live cell and `lo` the inward end of the maximal run `lo:top` of live cells whose condensed mass is below `ε_band` of their reference. If the run holds at least one cell below `ε_m`, every cell of the run transfers, from `top` inward, into the cell beneath it, so the run's whole inventory lands in cell `lo − 1` (the first live cell beneath the run) in one invocation; if `lo = 1` the run is the whole column and cell 1 survives holding it; if the run is the single last cell it is lost to the ledger (`_lose!`). After the pass, any live cell still below `ε_m` ends the solve as `:interior`, as now. With `ε_band = ε_m` the rule is today's extended to a contiguous block; the band exists because the surface cell crosses first (C48) and the synchronization at the threshold is finite (A.3), and it deactivates cells at up to ten thresholds, inside the reactivity switch's range, which is where R25's mechanism lives. Invariants: every transfer is the conservative volume merge (nothing is lost but the last cell's inventory); the run's inventories are summed once and the merged temperature inverted once, rather than chained through up to sixty sequential inversions and pore-volume rules as the prototype did (the chain is conservative, but each link re-inverts the caloric relation and applies the pore rule afresh, which with temperature-dependent densities creates volume differently from one merge); admissibility after the pass; the state size and pattern unchanged; the invocation recorded on every event, with the number of cells below `ε_m` and the number the band took (C53). Prototype: the overridden `deactivate!` with `BAND`, a multiple of `ε_m` there (A.3). Tests: §6.1 items 6 to 8.

The fixed-mesh rule is unchanged (D47); a fixed mesh has no correction (every node pinned), so V-BO-02 and V-BO-05 are unaffected by this plan except C45.

### 3.5 Event location

```julia
located_events::Bool          # integrate/solve keyword (the callback is built in integrate); default decided by the gate (§10)
```

With location on, the deactivation callback is a `ContinuousCallback` on the smallest margin `min_i (m_i^cond / (ε_m m_i^ref) − 1)` over the live cells, `affect_neg!` only, `rootfind = RightRootFind`, followed by the same `DiscreteCallback` as a safety net; the affect at a root nudges the crossing cell below the strict threshold by four ulps (the root is the threshold) before `deactivate!`, then runs the consistent sweep, the admissibility check and the restart as now. Measured on V-BO-02 (A.5): the event-time shift between reltol 1e-6 and 1e-8 falls from 0.88 s (0.133 of the median interval) to 0.15 s (0.022), at a cost of 2 to 4 % in wall time; without the safety net a root that leaves the cell at the threshold is missed and 21 of 60 events are lost. The located event is applied at an interpolated state of the interpolant's order (`Rodas5P`'s stiff-aware fourth-order interpolant, below the method's fifth), whose algebraic rows the consistent sweep rebuilds; the discrete safety net guarantees that no crossing is missed, and two cells crossing in one step are not missed because the margin has one root per step and the condition is re-evaluated after the state changes (asserted by §6.1 item 9). Under a correction a cell's mass is not monotone (the relative flux feeds it), so a cell can dip below `ε_m` inside a step and recover: location fires at the first crossing, the discrete predicate at a step's end, and the plan's semantics is the first crossing. The four-ulp nudge of the prototype is replaced by passing the located cell to `deactivate!` (non-strict at a root), the strict predicate kept for the safety net. Tests: §6.1 items 9 and 10.

### 3.6 The Jacobian pattern under a correction

A balance row of cell `k` reaches the positions of nodes `max(1, k−1)..min(N+1, k+2)` and the masses and energy of cells `max(1, k−1)..min(N, k+1)` besides what it reaches now (`pattern.jl:96-115`); the cell rows are unchanged (the reconstruction already reaches `i−2..i+2`). The design's §3.7 count with the correction active, `5NC + 13 + 2·bl + 11·χ`, is measured again with these rows (30 on V-BO-01 at `bl = 1`, `NC = 3`, against the formula's 30); the containment test runs with every node partially activated (`r_hi > 1 > r_lo`, so the activation is strictly inside `(0, 1)` at `t = 0` and its derivatives are exercised; a saturated switch has zero derivative) and on perturbed, zero and negative states, as in Phase 1. The energy columns are a structural superset the correction never reads (`V_eval` and the masses do not depend on `T`); they stay only if the optional reaction term is wanted. Prototype: the overridden `jacobian_pattern`, `missing_entries` empty, colored equal to dense (A.3), with the caveat that its first containment run saturated the activation and used `β = 0`, so the mass columns it added were trivially contained; the check was repeated with a partial activation and `β = 1` (A.10).

### 3.7 The constraint scale and the domain test

A balance row's scale becomes `(|θ_k| + |θ_A|) V_k + Σ_f A_f |w_f| + Σ_f A_f |c_f|` (C52); the velocity tolerance is unchanged. The admissibility predicate is unchanged (crossed nodes, negative masses, porosity, area). In an activated region `|c_f|` (about 5e-4 m/s) exceeds `|w_f|` (about 3e-5 m/s of burning) by an order, so after C52 the scaled test is loosest exactly where the correction acts; the gate reports the unscaled residual beside the scaled one. The round-off scale of the consistent sweep is a global sum over cells and faces and gains the sweep term the same way. The prototype's accepted residuals reach 0.69 of the bound and 52 proposals of 1,337 were rejected by the domain test at `τ_p = 0.2`; the gate records both after C52 (§7.3 item 6).

### 3.8 Trial-state contract additions

| Quantity | Guard | Where |
|---|---|---|
| Thickness ratio in the activation | `h` floored smoothly at `DISTANCE_FLOOR_FRACTION · h₀` before the logarithm; `a = 1` at or below zero thickness | `correction!` |
| Monitor density | the floored volume `V_eval` | `correction!` |
| Tabulated enthalpy at a non-finite temperature | returns `NaN` (C47) | `caloric.jl` |
| Surface Newton iterate | the trial clamp of §3.8 (plan) replaces the hard clamp and the warnings (C36) | `surface_balance.jl` |

The quality grep of plan §3.8 (no `clamp`, `max`, `min` on state-derived quantities in the new modules) covers `correction!`; the smoothstep's clamping of its argument is the one-sided `C¹` join the design allows for the kinetics limiter, and is reviewed as such.

### 3.9 Parameters and their defaults (provisional, gate-measured)

| Parameter | Default | From |
|---|---|---|
| `τ` (the layer's equidistribution time, `N²` in the formula) | set by the ladder at constant `τ`; the prototype's runs are `τ = 0.2 N²` s (720 s at 60 cells, `T_sync` 73 s) | A.3, A.4 and A.10 (the `N²` scaling); the heuristic `T_sync` at most a tenth of the layer's remaining life |
| `r_hi`, `r_lo` | 0.3, 0.1 | design §3.3; A.4 measures 0.5/0.2 |
| Monitor `α`, `β` | 1, 0, normalized by the mean over the layer (`β` acts only on mixed compositions; unnormalized it is a hidden `τ`, A.4) | design §3.3 |
| Band `ε_band` | 0.1 (absolute; the prototype's 10 `ε_m` at `ε_m = 1e-2`) | A.3 (3 `ε_m` fails, 10 `ε_m` passes at 60 cells); unmeasured at the other thresholds; fixed before the ladder, not tuned on it (§7.3 item 1) |
| `ε_m` | 1e-2 | Phase 1; 3e-3 and 3e-2 measured (A.4) |
| Switch width | 1.5 `ε_m` | C27; measured on and off (§7.3 items 1 and 2); R25 argues for off under a correction |
| Lagrangian nodes | an explicit list given to `setup` (the burnout decks know their interface node); a composition rule is a cutover convenience, not Phase 2's | §3.1 |
| `located_events` (an `integrate`/`solve` keyword) | off until the gate | §3.5 |
| Restart after an event | solver-chosen (`Rodas5P`) | A.6: forced 1e-6 s costs 30 % on V-BO-02 and cuts the ledger quantity 20× |


---

## 4. Kernel-by-kernel plan

In build order. Each entry: the equations, the current implementation to read first, what changes, and the case the item lands with. Every item lands with its dense-AD containment (the configurations of §6.1 item 11 at 6 and 8 cells) before the next starts, as in Phase 1.

### 4.1 The correction (Closures)

*Equations.* §3.2 of this plan; the balance row of cell `k` becomes `Σ_f A_f (w·n)_f − (θ_k − θ_A) V_k − s_k = 0` with `s_k = Σ_f A_f (c·n)_f`, so `w − c` satisfies the volume balance and is the solid velocity; the pinned node's row `w_p = 0`; a Lagrangian node has `c = 0` and no row of its own (its velocity follows the balance rows as now). *Read first.* `src/Closures/balance.jl:139-200` (`_velocity_rows!`, `_consistent_velocities!`, `_coefficient`, `_sweep_rate`), `src/States/derived.jl` (`CellReference`, `DerivedState`, `derive!`), `src/Meshes/topology.jl` (`lagrangian`), `src/Integration/setup.jl:300-420`. *Changes.* The `Correction` types and `correction!` (§3.1); `ds.c`; `CellReference.thickness`; the sweep term in both velocity functions; `setup` builds the correction from its keywords, takes the layer-interface nodes as an explicit list (§3.9) and passes them to `line_mesh`; `PrescribedCorrection` replaces the `override` field of `BalanceRowClosure` (one path for both; the field is read in `balance.jl`, `pattern.jl` and `checks.jl`, and `test/newcore/test_mesh_motion.jl` constructs a closure with the keyword and changes with it). The round-off check of `_consistent_motion!` (`setup.jl`) gains the sweep in its scale as C52's does. *Lands with.* The unit tests of §6.1 items 1 to 4; V-MMS-02/03/04 and the Phase 1 relative-flux tests unchanged under `PrescribedCorrection` (bit-identical to Phase 1, §7.3 item 5).

### 4.2 Relative fluxes under the correction (Kernels)

*Equations.* §3.3. *Read first.* `src/Kernels/relative_flux.jl` (`_stencil`, `_face_value`, `contribute!`), `test/newcore/test_relative_flux.jl`. *Changes.* `u = −c_f n_f` from `ds.c`; the kernel installed whenever the correction is not `NoCorrection`; the `setup` refusals of an override on a reacting material and of deactivation under an override stay for `PrescribedCorrection` and do not apply to the equidistribution correction. *Lands with.* Isothermal advection under an `EquidistributionCorrection` with reactions off (§6.1 item 5); V-BO-01 at `c ≠ 0` (§6.2).

### 4.3 The block-and-band rule (Deactivation)

*Equations.* §3.4. *Read first.* `src/Deactivation/events.jl` (`deactivate!`, `_transfer!`, `_lose!`, `_mask!`), `test/newcore/test_deactivation.jl` ("every cell below the threshold goes in one pass", "a void under a live cell ends the solve"). *Changes.* The run from the top within the band; the chain of transfers into the cell beneath the run; the survivor rule; the invocation index on `DeactivationEvent`; `replay_events!` unchanged (it reads times). The band parameter is an absolute mass fraction `ε_band ≥ ε_m` (default 0.1), not a multiple of `ε_m`: at `ε_m = 3e-2` a band of ten thresholds took cells at 29 % of their mass (A.4), which is not a near-threshold band. *Lands with.* V-BO-04 re-pinned (C46: one invocation, one survivor, the exact law through the event); the tests of §6.1 items 6 to 8; V-BO-02 and V-CPL-02/02b unchanged (a fixed mesh and a one-cell-at-a-time ablation front never form a block).

### 4.4 The Jacobian pattern (ColoredJacobian)

*Equations.* §3.6. *Read first.* `src/ColoredJacobian/pattern.jl:96-115`, `test/newcore/jacobian/test_pattern.jl`, `test_containment.jl`. *Changes.* The velocity rows' reach under a correction; the color-count formula gains the correction configurations; containment with every node activated (§3.6). *Lands with.* Containment on V-BO-01's configuration (in-depth radiation, three components, moving mesh) and on the V-BO-04 and heterogeneous-layer configurations, at 6 and 8 cells, five states, both energy forms.

### 4.5 Setup, the domain test and the outputs (Integration, Outputs)

*Read first.* `src/Integration/setup.jl` (`_require_supported`, `setup`), `checks.jl` (`constraint_residual`), `solve.jl` (`_deactivation_callback`, `_observe!`), `src/Outputs/solution.jl` (`CoreExtras`, `SolverStatistics`). *Changes.* Keywords `correction`, `lagrangian`, `activation`, `band`, `switch_width`, `located_events`; the sweep in the constraint scale (C52); `correction_velocity` and `min_cell_thickness` series, `invocations` in the statistics; the trial clamp in the surface Newton (C36: the iterate is clamped smoothly into the property bracket as every other temperature-sensitive kernel's input is, the warnings removed, the old core unchanged because the Newton body is shared and the clamp is a keyword the new core sets). *Lands with.* The keywords and outputs at M3; V-BO-01 and V-BO-03 are routed at M4 (§6.2, §7); C36 lands at M1 with its test (a wild trial state runs the Newton without a warning and with finite derivatives on the new core; the old core keeps its clamp and warnings).

### 4.6 Event location (Integration)

*Equations.* §3.5. *Read first.* `solve.jl:378-420`, the prototype's continuous variant (A.5), SciMLBase `ContinuousCallback` (`affect_neg!`, `RightRootFind`). *Changes.* Behind `located_events`; the discrete safety net kept; the nudge at the root; the measurement of §7.3 item 7 decides the default. *Lands with.* §6.1 items 9 and 10; V-BO-02 at both tolerances (the event-time shift as a QoI of the gate script).

### 4.7 What is not built

The slot-shifting interior merge (C54), the Péclet weighting of the Darcy upwind (C21: the default integrator has no Newton to fail at the kink, and the moving-mesh decks of Phase 2 do not reverse their flow; left to the owner's call), a centred face temperature for the relative energy flux (C33: kept upwind-biased), the monitor's reaction-rate term `γ` (design §3.3 marks it optional; the normalized mass term already acts only on mixed compositions, §3.2), and any change to the energy forms.

---

## 5. Solver plan

### 5.1 The integrator and its policies

`Rodas5P` with the solver-chosen restart stays the default (C28, C40). On the event-dense fixed-mesh deck the alternatives measured in this session (A.6): the forced restart (`DT_AFTER_EVENT = 1e-6` s) takes 2,790 accepted steps against 1,706, 30 % more wall time, and cuts the ledger quantity of C39 from 1.5e-3 to 7.6e-5 (the post-event gas transient is then resolved); `KenCarp4` with the raw estimate takes 10,021 steps, 1,784 nonlinear failures and 3.75× the time. The restart policy is therefore a gated choice for deactivating solves (§10): the forced restart buys the old core's ledger tolerance back at a measured price.

### 5.2 Event location

§3.5 and §4.6. The continuous callback's cost is 2 to 4 % of the solve; its benefit is the event time tied to the tolerance (shift 0.022 of the median interval between reltol 1e-6 and 1e-8, against 0.133 for the discrete callback). The gate item 8 of Phase 1 (achieved error against requested tolerance) is re-run on the deactivating cases with location on and off (§7.3 item 7); the default follows.

### 5.3 The constraint check and the domain test

C52's scale; `MAX_CONSTRAINT_REJECTIONS` unchanged; the largest accepted residual and the domain rejections per solve are gate quantities (the prototype's 0.69 and 52 of 1,091 proposals on V-BO-01 are the numbers to beat, A.3).

### 5.4 Minimum cell size

The design asks Phase 2 to measure, not assume, the smallest cell the correction produces. Measured on the prototype (A.3, A.4): 1.0 µm at 60 PMMA cells, 2.0 µm at 30, 0.49 µm at 120, 0.31 µm at `ε_m = 3e-3`, all about `ε_m`-to-`ε_band` of the initial cell, which is the design's "bounded relative to the layer": the layer's cells shrink with the layer until the band takes them. No step fell below 2e-5 s and no solve stalled; the conduction coefficient `k/h²` of a micrometre cell is 2e11 W m⁻³ K⁻¹ and `StructuredKLU` factored every iteration matrix. The `min_thickness` termination (D18) measures the whole slab and does not bind. The gate records `min_cell_thickness` on every run (§7.3 item 9).

### 5.5 Wall time: where the prototype's cost is, and the levers

At 60 PMMA cells the prototype takes 2.8 to 3.5 s against the old core's 0.58 to 0.80 s (A.3: 3.5 to 6×), at 120 cells 9.3 s against 2.1 s and at 30 cells 1.46 against 0.25 s (A.4; the old core's times are A.2's medians from another process, so the ratios carry the machine's spread); V-BO-02, which has no correction, takes 1.32 s against 1.29 s in one process (A.1). The correction's cost has three measured parts:

1. **A Jacobian every step at 30 colors.** `Rodas5P` builds 1,091 Jacobians over 1,091 accepted steps, where the old core built 164 over 853; the correction's pattern colors in 30 against 18 at `c = 0` (C43), so one Jacobian costs about 1.7× more. Roughly half of the solve is Jacobians. A Rosenbrock-W method would allow reuse, but none of the installed W-methods (`ROS34PW2`, `Rodas23W`) reaches the old core's accuracy within 3× of its time even on the smooth audit decks (A.7), and OrdinaryDiffEq rebuilds their Jacobian every step as well; `KenCarp4`, which reuses its Jacobian, is 3.75× slower than `Rodas5P` on the event-dense deck (A.6). The integrator's accuracy per tolerance, Phase 1's open lever, does not close on the moving-mesh deck with any candidate (A.7: the best, `Rodas5P` at reltol 3e-7, is 1.14× the old core's error at 1.52× its time; `Rodas5Pe` at 1e-7 is 0.28× at 1.77×). The tolerance weights of the algebraic rows are probed in A.9.
2. **Rejected steps.** 246 of 1,337 proposals, 52 of them by the domain test on the constraint residual (C52 is expected to remove most of those).
3. **The steps themselves.** 1,091 against 853: the correction adds stiffness at `τ = 0.2` s (at `τ = 1`, 691 steps; at `τ = 0.1`, 1,286).

The gate's wall-time clause is measured at matched accuracy (the Phase 1 protocol: a run's error is its largest temperature difference against the same core at reltol 1e-11) on V-BO-01, 03 and 04 at 30, 60 and 120 cells, both cores in one session, three repeats; the result is reported as a ratio with the machine's spread, and the owner decides whether C40's exception extends to it (question 6). The plan does not promise parity: the prototype is 3.5 to 6× at equal tolerance before C52 and before any tuning of `τ`, matched accuracy raises the ratio, and the levers above do not reach 2×; the honest expectation is 4× or more, and whether any ratio is accepted is question 6.

### 5.6 Energy form

The temperature form stays the default; the enthalpy form reproduces every V-BO-01 quantity with the correction at 3.1 to 3.9× the temperature form's time (A.4: 10.9 s against 2.8 to 3.5 s), consistent with Phase 1's 2 to 13×. The gate runs the V-BO ladder in both forms (§7.3 item 3).

### 5.7 Performance: why the burnout path is slower, and the milestone that must close it

**The owner's rule.** The new core must run the burnout path within the same order of wall time as the old core at matched accuracy, or faster; a path that cannot is a bad sign for the architecture, not a detail to carry. Milestone M2 exists to settle that before the correction is built, with a stop at its gate: within 1.5× passes, 1.5 to 2× is the owner's call, above 2× the phase stops for the owner's decision (continue with the levers below, ship `c = 0` with the fallback, or reconsider the cutover itself).

**Where the time goes** (V-BO-01 at 60 PMMA cells, equal tolerance; A.2, A.3):

| | old core | new core with the correction |
|---|---|---|
| accepted steps (rejected) | 853 (3) | 1,091 (246) |
| Jacobians | 164, the reduced approximate one, reused across steps | 1,091, exact, 30 colors (about 15 to 25 residual-equivalents each) |
| right-hand sides | 19,076 | 11,787, plus 1,337 residuals for the domain test |
| residual cost | 1 | 1.15 to 1.54 (Phase 1 leftover 1) |
| wall | 0.58 to 0.80 s | 2.8 to 3.5 s |

Four causes, none intrinsic to the design: (1) `Rodas5P` needs an exact Jacobian every step by construction, about half the solve, where `KenCarp4` on the old core reuses an approximate one over about five steps; it was chosen because the Newton-based candidates behaved badly on the new residual (`KenCarp4` on V-BO-02: 10,021 steps and 1,784 nonlinear failures against the old core's 1,572 and 56; on the Eulerian audit deck its error does not fall with the tolerance, A.7), and that behaviour is undiagnosed; (2) the correction's reconstruction widens the pattern from 18 to 30 colors, 1.7× per Jacobian; (3) more steps and rejections: the correction's stiffness (691, 1,091 and 1,286 steps at `τ_p` of 1, 0.2 and 0.1), an 18 % rejection rate and one residual per proposal for the domain test, and at matched accuracy `Rodas5P`'s reltol 1e-7 to 3e-7 against the old default's 1e-6; (4) the residual itself at 1.15 to 1.54× with no kernel ever profiled or tuned. The old core's path carries ten months of performance work (KLU, the structural diagonal, the ring trim, the reduced Jacobian reused across steps, the step reset after a park); the new core has had none.

**The levers, in the order the measurements rank them.** (a) Diagnose the Newton-based integrators on the new residual: the suspects are the error weights of the mass rows (`reltol · ε_m · m(0)`, which A.9 shows carry the accuracy), the smoothed estimator on the new variables, the one-sided limiter and smooth floors near depletion, and the restart after events; a reusing integrator that qualifies cuts the Jacobian count five-fold and is the one lever that can reach parity alone. (b) Failing that, Jacobian reuse for a Rosenbrock-W method (`RosenbrockW6S4OS` is fourth order and tolerates an inexact Jacobian); OrdinaryDiffEq rebuilds the Jacobian of its W-methods every step, so this is a small stepper change, upstream or local. (c) The per-Jacobian cost: analytical fast paths for the conduction and gas-transport kernels inside the dual passes (declined in Phase 1 when the Jacobian was a quarter of a solve, worth revisiting at half) and a compact smooth stencil for the relative flux. (d) Steps and rejections: C52, the controller's gains, the largest `τ` the band permits. (e) A kernel performance pass on the residual (the boundary-condition kernels led the Phase 1 profiles). Stacked without (a) these plausibly reach 1.5 to 3×; parity or better needs (a) or (b).

**What is measured.** `burnout_timing.jl` with a cost breakdown per solve (Jacobians, right-hand sides, factorizations, solves, the domain test, the remainder) on V-BO-02, on the prototype's correction run of V-BO-01 (the prototype is enough for M2: the lever is the integrator and the Jacobian, not the correction's code) and on the audit decks, both cores in one process, before and after each lever, so the gate's ratio has its anatomy beside it.


---

## 6. Test and verification plan

### 6.1 New tests (`test/newcore/`)

1. **The correction's zeros and smoothness** (`test_correction.jl`). `c = 0` at the pinned, surface and Lagrangian nodes and wherever both cells are above `r_hi`; the activation is `C¹` in `log(h/h₀)` (finite differences across both joins); the sign moves a node into the thicker cell; at or below zero thickness the activation is one and finite.
2. **Relaxed equidistribution.** On a three-cell inert mesh with every node activated (`r_hi` above one), `M_i h_i` relaxes to a common value with a rate set by `τ`, against the closed form of the linearized recurrence.
3. **The balance rows with the sweep.** After the consistent sweep `w − c` satisfies every cell's volume balance to round-off; the scaled constraint residual (C52) is at round-off; a uniform strain rate with `c ≠ 0` gives `w = v_s + c` exactly on planar, cylindrical and spherical meshes.
4. **`PrescribedCorrection` is the Phase 1 override.** V-MMS-03's trajectory hash under `w_override` equals Phase 1's (bit for bit); `test_relative_flux.jl` unchanged.
5. **Isothermal advection and free stream under an `EquidistributionCorrection`** with reactions off: temperature constant to 1e-9 K, total mass to round-off, a uniform concentration uniform.
6. **The block rule** (`test_deactivation.jl`). Every cell of a column crossing at once ends in one invocation with one survivor whose inventories are the column's; a run of `k` cells under a live cell merges into that cell in one invocation; every inventory conserved through the chain; admissibility after the pass; invocation ids on the events; `park_times` one per cell.
7. **The band.** A cell at 1.5 `ε_m` beneath a top cell at 0.9 `ε_m` merges with it; a cell above the band stays; the band is an absolute fraction (`ε_band = 0.1` takes the same cells at `ε_m = 1e-2` and `3e-3`).
8. **An interior void beneath a healthy skin still ends the solve** (`:interior`).
9. **Event location.** On a two-cell deck the located event time's shift between tolerances falls with the tolerance; the safety net deactivates a cell a root left at the threshold (the nudge); location and the discrete predicate agree on the cell.
10. **A located event's state** is admissible and leaves the velocity and intensity rows at round-off after the sweep.
11. **Containment** with the correction active (every node activated and the perturbed, zero and negative states) on the V-BO-01 (in-depth, three components), V-BO-04 and heterogeneous-layer configurations, at 6 and 8 cells, both forms; colored equal to dense; the color count against the formula.
12. **Trial states.** `sensible_enthalpy` of a non-finite temperature is `NaN` (C47); the relative-flux kernel on a `NaN` state returns a `NaN` residual without throwing; the activation at negative thickness.
13. **The surface Newton's trial clamp** (C36): a wild trial state runs without a warning, with finite values and derivatives; physical states unchanged bit for bit.
14. **Outputs.** `correction_velocity`, `min_cell_thickness`, the events' invocations and the `invocations` statistic; NaN beyond the live cells after a block merge.
15. **Allocations.** Zero bytes per residual and per colored Jacobian with the correction present (the Phase 1 gates extended).
16. **Determinism.** Repeated solves bit-identical with the correction and events, under `StructuredKLU` and UMFPACK.
17. **The fallback, only if built** (C54, §7.4). An interior cell's merge and the slot shift conserve every inventory, shift every per-slot quantity listed in C54, keep the mask contiguous from the top, leave the state admissible and the algebraic rows consistent, and V-BO-01 and V-BO-03 pass their QoIs at `c = 0`.

### 6.2 Verification cases

| Case | Routing on the new core | Acceptance |
|---|---|---|
| V-BO-01 | routed with the correction (`correction = :equidistribution`, the defaults, `lagrangian` the interface node) | every current QoI at its tolerance (the case's fixed references) at 30, 60 and 120 PMMA cells, the 30-cell row an owner exception if R25 stands; invocations exactly one for the PMMA layer as a new QoI, with the largest merged fraction at most 5 `ε_m`; the surface-Newton warning count replaced by the count of clamped trial iterates once C36 removes the warnings (question 8) |
| V-BO-02 | as now (fixed mesh, no correction) | unchanged; the ledger QoI at C39's 3e-3 unless the restart policy changes (question 7) |
| V-BO-03 | routed with the correction | its four QoIs as pinned (the transmitted-flux warning becomes the `radiation_transmitted` series' peak, the surface-Newton count the clamped-iterate count; owner question 8) |
| V-BO-04 | routed, no correction needed (all cells cross together; the block rule) | re-pinned per C46: one invocation, one survivor, the exact law through the event, isothermal, `ξ_s = ρ_s` on live cells, final outcome |
| V-BO-05 | routed (fixed mesh) | re-pinned per C45: `:interior` at the first `ε_m` crossing beneath a live cell, the stop time from a ladder |
| V-BO-06 (new) | V-BO-01 on a mesh graded toward the surface: `_vbo_problem` gains a grading keyword (a geometric ratio of 1.02 per cell across the 60 PMMA cells, so the surface cell is a third of the interface cell; the Kaowool at the interface pitch), the same materials, conditions and 500 s | one invocation; the V-BO-01 case's fixed references at their tolerances (they are mesh-independent); the old core at the same graded mesh reported beside |
| V-BO-07 (new) | a heterogeneous layer of mixed composition: every PMMA cell also holds a second solid of the same properties but half the density (605 kg m⁻³) and a pre-exponential factor a tenth of PMMA's, with a mass fraction graded linearly from 0 at the interface to 0.5 at the surface, so the density varies inside one activated region and the normalized monitor's mass term acts (two pure sub-layers separated by a Lagrangian node would equidistribute separately and `β` would change nothing at equilibrium); the Kaowool substrate, conditions and 500 s of V-BO-01 | one invocation for the layer; mass loss equal to the layer's initial mass (exact); ledger at 1e-3; `β = 1` against `β = 0` on invocations, merged fractions and MLR; the old core at the same deck reported (peak and timing have no fixed reference until the CI run pins them) |
| V-BO-08 (new) | localized burnout: a 1 mm PMMA layer (ten cells at the V-BO-01 pitch) on the 12 mm substrate, the V-BO-01 materials and conditions | one invocation; `correction_velocity` exactly zero at every Kaowool node on every saved step; the Kaowool temperatures within 1e-10 of the `NoCorrection` run until the first activation; mass loss exact; ledger at 1e-3; the MLR curve against the old core at the same deck |
| V-CPL-02/02b, V-CPL-03, V-MMS-02/03/04 and every Phase 1 case | unchanged routing | bit-identical to Phase 1 with `NoCorrection` (§7.3 item 5) |
| The Douglas fir deck (`examples/07_douglas_fir_cone.jl`, 60 s and 600 s) | a gate script, not a case | bit-identical with `NoCorrection`; to 1e-10 with an inactive correction; wall time recorded on both cores (A.8); the matched-accuracy comparison is Phase 3's (question 9) |

### 6.3 Gate scripts (`test/newcore/gates/`)

`burnout_ladder.jl` (item 1 of §7.3: the matrix of cases, grids, thresholds, monitor settings and `τ`; invocations, outcome, QoIs, wall time, `min_cell_thickness`, constraint residual, domain rejections), `burnout_identity.jl` (item 5: hashes of the charring decks and the Douglas fir deck with `NoCorrection` and with an inactive correction), `burnout_timing.jl` (item 8: both cores at matched accuracy, the Phase 1 protocol, with the cost breakdown of §5.7), `burnout_gradients.jl` (item 9: the Phase 0 `fd_gradients.jl` protocol on the new core), `case_accuracy.jl` extended to the deactivating cases with location on and off (item 7), `energy_ledger.jl` extended to a burnout deck (item 3). Every table stamped with its commit (D52); the measurements under `measurements/`.

---

## 7. Milestones, gates, estimates

Single developer; the old core untouched (except C55 and C36, fingerprint-checked); targeted gates per milestone (plan §9.3, D21); the full three-way `Pkg.test()` at M7. The design estimated 4 to 6 weeks for Phase 2; the estimate below is 11.5 to 14.5, because the burnout path's wall time must be brought within the old core's order first (M2, §5.7), and because the cascade rule, the event location, the three new cases and the gradient re-measurement are concrete items the design listed in one line each, and because the wall-time clause needs the matched-accuracy protocol on a ladder.

| Milestone | Builds | Gate to pass before the next | Estimate |
|---|---|---|---|
| M0 Branch and plan | This plan, the prototype and its tables (§9.2) | Format check, Aqua and JET green locally (the Quality workflow runs only for `main`, §9.1); the owner's answers to questions 1 to 4 gate M1's start, not M0's length | 2 days |
| M1 Deactivation | The block-and-band rule, invocation ids, the survivor and layer-end rules, C47, C36 | §6.1 items 6 to 8, 12, 13; V-BO-04 re-pinned and routed; V-BO-02, V-CPL-02/02b bit-identical; for C36 (an edit of shared old-core code) the fingerprint, `test/solver/test_surface_bc_guards.jl`, `test/boundary_conditions/test_boundary_temperature.jl` and the full old-core `Pkg.test()`, the old core keeping its clamp and warnings (keyword off) | 1 week |
| M2 Performance of the burnout path | The diagnosis and the levers of §5.7, in the order the measurements rank them: why the Newton-based integrators fail on the new residual (the erratic `KenCarp4` error and its nonlinear failures), Jacobian reuse (a qualifying reusing integrator, or a reuse policy for a Rosenbrock-W method), the per-Jacobian cost (fast paths for the conduction and gas-transport kernels in the dual passes, a compact relative-flux stencil), the step and rejection count (C52, the controller), and a kernel performance pass on the residual; `burnout_timing.jl` with its cost breakdown | V-BO-02 (the burnout path without a correction) and the prototype's correction run on V-BO-01 at matched accuracy against the old core, both cores in one process: within 1.5× passes (the same order); 1.5 to 2× is the owner's call; above 2× is the bad sign of §5.7 and the phase stops at this point for the owner's decision | 2 to 3 weeks |
| M3 The correction | `Correction` types, `correction!`, `ds.c`, `CellReference.thickness`, the sweep in the balance rows and the consistent sweep, the Lagrangian nodes, `RelativeFluxKernel` on `ds.c`, the pattern (C43), C52, `setup` keywords, outputs | §6.1 items 1 to 5, 11, 14 to 16; V-MMS-02/03/04 and `test_relative_flux.jl` unchanged under `PrescribedCorrection`; V-BO-01 at 60 cells ends in one invocation with every QoI within its tolerance (the prototype's A.3 run, not within round-off of it: C52 changes the rejected steps) | 2 weeks |
| M4 The ladder | V-BO-01 and V-BO-03 routed; the matrix of §7.3 item 1 at the first thresholds; the end-of-layer tail diagnostic (R25); both energy forms; `burnout_ladder.jl`, `burnout_timing.jl`, the `energy_ledger.jl` extension; `min_cell_thickness`; a tuning allowance of half a week for `τ`, the band and the activation | §7.3 items 1 to 3, 6, 8, 10 on V-BO-01/03/04; the decision point of §7.4 | 2 to 2.5 weeks |
| M5 Events | Event location behind `located_events` (prototyped, A.5), the restart-policy measurement on V-BO-01 | §6.1 items 9, 10; §7.3 item 7; the defaults of questions 5 and 7 chosen | 0.5 week |
| M6 The gate's checks | V-BO-06/07/08 (a graded mesh with its old-core reference, a new two-solid material, a thin layer); `burnout_gradients.jl` (nine runs); the Douglas fir identity and timing; V-BO-05 re-pinned from a ladder; the rest of item 1's matrix | §7.3 items 1, 4, 5, 9 complete | 2 to 2.5 weeks |
| M7 Gate document | `docs/Mesh_Framework_Phase2_Gate.md`, §12 of this plan, the three-way `Pkg.test()` in detached worktrees | The Phase 2 gate (§7.3) | 1 week |
| Review | The adversarial review of the phase and its followups, as Phase 1 experienced it (D38 to D52, C32 to C39, a re-measured gate) | The review's findings closed or owned | 1 to 2 weeks |

Total: 11.5 to 14.5 weeks (the design said 4 to 6); 13 to 16 if the fallback of §7.4 is built. Phase 1's record is an aggregate (about twelve and a half of thirteen to fourteen weeks, plus an unbudgeted review); its per-milestone actuals were never recorded, so these rows are estimates against the deliverable count (about 25 code-level items in §2.2 and §4, 17 tests, four new gate scripts, three new cases, four re-pins), not against measured velocity.

### 7.1 What the owner decides before M1

1. Question 1 (the interior policy and whether the slot-shifting merge is insurance now or the fallback later).
2. Question 2 (the band as an absolute fraction, its default).
3. Question 4 (thresholds as material parameters).
4. Question 8's V-BO-04 re-pin (M3 needs it).

### 7.2 Design errata to apply at the design's next revision

§3.3: the activation is per node, the smooth "or" of its two cells' switches (C44); §3.5: the cascade is the block-and-band rule, the band an absolute mass fraction (C42); §3.7: the color count with the correction active includes the velocity rows' reach (C43) and the `c = 0` count is `3NC + 7 + 2·bl + 7·χ` (C25); §3.5 and §6 Phase 2: "one event per layer" counts callback invocations (C53). No design content changes beyond these; the owner applies them.

### 7.3 The Phase 2 gate, as written in the design, made measurable

Every item is measured and reported in `docs/Mesh_Framework_Phase2_Gate.md` with its script and table.

1. **One event per layer or fewer across the ladder.** `burnout_ladder.jl`, table `measurements/burnout_ladder.md`: V-BO-01 at 30, 60 and 120 PMMA cells × `ε_m ∈ {3e-3, 1e-2, 3e-2}` at the default `τ` and `ε_band`, plus `τ ∈ {0.1, 0.2, 1}` at 60 cells and `ε_m = 1e-2`; V-BO-07 with `β ∈ {0, 1}` (the monitor's mass term acts only where the density varies, so the "two monitor settings" are carried by the heterogeneous layer; on V-BO-01 the prototype's `β = 1` was a `τ/2` duplicate, A.4); V-BO-06 and V-BO-08 at their CI resolutions and the three thresholds; V-BO-04 with and without the correction. Per run the table records the invocations, and per invocation the cells below `ε_m` and the cells the band took with the largest merged mass fraction. Pass: one invocation per layer in every run, with the largest merged fraction at most 5 `ε_m` (so the band, fixed before the ladder, cannot be widened to force the count; the prototype's 3.6 % at 60 cells passes, its 10 % at 120 cells fails). The prototype already fails this strict criterion on three rows (A.3, A.4): 120 cells (two invocations, the band saturated at 10 %), `ε_m = 3e-2` (two) and `τ = 1` (two); the 120-cell and `τ_p = 1` rows are the `N²` scaling of §3.2 (a constant layer time `τ` across the ladder is the lever, A.10), the `3e-2` row the absolute band's 3.3 `ε_m` there, which 3 `ε_m` failed at 60 cells and any run left at two invocations is an owner exception recorded per run, not a pass. Every run of the matrix is also judged on all of its case's QoIs (item 2).
2. **Peak MLR and its timing.** `burnout_ladder.jl` runs every V-BO case of the matrix in the suite's assert mode against the cases' own fixed references (V-BO-01: peak 32.39 within 3 %, its time 335.0 s within 3 %, the burnout time 363.5 s within 1 %, the mass loss within 0.3 points, the final thickness, the ledger, the surface-Newton count), and reports the old core at the same resolution in the same session beside each run. Pass: every run of the matrix passes its case, the 30-cell row included unless the owner excludes it (question 8); the prototype fails it at 30 cells (peak +13.7 %, `t_peak` +8.4 %) and at `τ = 0.1` (`t_peak` 361.5 s, +7.9 %), which the tail diagnostic of M4 (R25) must resolve. Because the merge burst can sit below the main peak, the item adds a curve criterion the case already carries (`:mlr_history`): the largest difference from the old core's MLR after `t_peak` at most 3 % of the peak.
3. **Ledgers closed.** The mass ledger of every V-BO run at its case's tolerance in the temperature form (V-BO-01: 1e-3, the prototype 7e-5 to 2.5e-4; V-BO-02: C39's 3e-3 unless question 7 changes the restart); the energy ledger (`energy_ledger.jl` extended to V-BO-04 at 40 cells, whose event is one and exactly timed, over the event-free window before it and, with located events, across it, the event's vented gas, lost inventory, regularization and created volume entering both closures as the event records them, the sum over live cells only) no larger than the old core's closure on the same deck and converging at the integrator's order under fixed-step refinement on the event-free window; exact in the enthalpy form.
4. **Admissibility at every event.** The deactivation callback asserts it (Phase 1); `burnout_ladder.jl` counts events, invocations and admissibility failures (zero) over the matrix.
5. **Charring decks identical to Phase 1.** `burnout_identity.jl`, run first at M0 on `bb2860b` to record the Phase 1 digests (as `oldcore_fingerprint.jl` holds 17bc477's): the hash of `(T, m, x, t)` of V-CPL-02/02b, V-CPL-03, V-CPL-06's char-layer column, V-RX-02d (a residue former), V-HC-05, V-RAD-01 and the Douglas fir deck (60 s) with `NoCorrection` equals that digest; with the correction present every output agrees to 1e-10 relative on the same accepted-step sequence (a changed step sequence is reported, with the temperature difference, as Phase 1's 7e-6 K on ale_darcy 25 was), the largest activation over the run is exactly zero and each deck's smallest `h/h₀` is recorded (C56).
6. **Constraint residuals and the domain test.** `burnout_ladder.jl` records per run the largest accepted scaled residual under both the Phase 1 scale and C52's, and the domain rejections before and after C52 on the same commit (the domain test rejects above one, so "at or below one" cannot fail and is not the criterion; the prototype's 0.69 and 52 of 1,337 proposals are the Phase 1-scale numbers).
7. **Event timing.** `case_accuracy.jl` on V-BO-01, 02, 04 and V-CPL-02/02b with location off and on: the achieved temperature error over the cells live in both solves, and the event-time shift between the tolerance and the tolerance divided by 100 as a fraction of the median event interval. Phase 1's budget (the shift at most ten times `reltol`) is unreachable by three orders of magnitude (A.5: 0.133 off, 0.022 on), so the criterion is restated and the script's budget changed: pass when the located shift is at most 0.05 of the median interval and at least three times smaller than the discrete one on every deactivating case (A.5 passes both on V-BO-02); the default follows (question 5).
8. **Wall time no worse than the old core's burnout runs.** `burnout_timing.jl`, table `measurements/burnout_timing.md`: V-BO-01, 03 and 04 at 30, 60 and 120 cells, both cores in one process with the same statistic (median of three), at matched accuracy with an error measure defined for decks with events: the temperature error over the saved times and cells live in both solves plus the event-time shift (the `case_accuracy.jl` measure), each core against its own reference at reltol 1e-9 (a reltol 1e-11 reference of a 60-event solve is not shown to complete; the script checks that the reference it uses does); and, for the cost of a correction present but inactive (the design's "checked on the PMMA cone cases"), the audit decks (`ale`, `ale_darcy` at 25 and 100 cells) and V-CPL-02/03 with the correction object installed, against Phase 1's numbers. Pass: the same order as the old core, that is a ratio at most 1.5 on every deck at matched accuracy; 1.5 to 2 is the owner's call; above 2 is the bad sign of §5.7. Standing before M2, from the tables: the prototype is 3.5 to 6× at equal tolerance (A.3, A.4: 2.76 to 3.5 s against 0.58 to 0.80 s at 60 cells, 9.3 against 2.1 s at 120, 1.46 against 0.25 s at 30) and matched accuracy raises the ratio (Phase 1 item 2: `Rodas5P` needs reltol 1e-7 to 3e-7 to match the old core's error); M2's levers (§5.7) must bring it within the pass before M3 is built, and the item is re-measured on the finished correction at M4.
9. **Smooth finite-difference gradients.** `burnout_gradients.jl`, table `measurements/burnout_gradients.md`, the Phase 0 item 1 protocol on the new core with the correction (nine runs, the PMMA absorption coefficient perturbed by up to ±0.8 %, central differences at four steps). Phase 0 already found the integrated quantity smooth through 60 parks (spread 0.025) and the event-timed ones not (`t_peak` interpolated 0.11, `t_bo` interpolated 0.20); the clause's target is the event-timed quantities. Pass: `peak_mlr` (interpolated), `t_peak` (interpolated), `t_bo` (interpolated) and the invocation time each with a spread below 0.1 across the four steps, and the invocation count the same for every perturbation; a count that changes fails the item, with the design's bracketed differencing as the documented consequence.
10. **Minimum cell thickness.** `burnout_ladder.jl` records the smallest live cell per run; pass when it is at least `0.5 ε_m h₀` on every run (the prototype: 1.03 µm against `ε_m h₀ = 0.98` µm at 60 cells, 2.0/1.97 at 30, 0.49/0.49 at 120, 0.31/0.29 at `ε_m = 3e-3`) and no run ends in `DtLessThanMin`.
11. **The fallback.** If items 1 or 2 fail on V-BO-01 at the defaults after M4's tuning, or item 8 fails and the owner grants no exception, the core ships with `c = 0`, the block-and-band rule and the slot-shifting interior merge (C54): that departs from the design's "nothing else changes", because at `c = 0` the new core cannot finish V-BO-01 or V-BO-03 without an interior policy (A.1). The fallback's own acceptance: V-BO-01 to V-BO-05 pass their cases on the new core at `c = 0` (one invocation per cell on V-BO-01 and 03, as the old core's parks), the ledgers and admissibility of items 3 and 4, item 8's timing with the old core's one-cell events as the comparison. The cutover case is then argued on the exact Jacobian, one-invocation cascades and multi-D readiness alone (design §6).

### 7.4 The decision points at M2 and M4

M2 ends on §7.3 item 8 measured on V-BO-02 and the prototype's V-BO-01 run with the levers of §5.7 applied: within 1.5× the phase continues; above 2× it stops for the owner's decision (continue, ship `c = 0` with the fallback, or reconsider the cutover), and between the two the owner decides.


M4 ends with the matrix of item 1 at 60 and 120 cells. If V-BO-01 ends in one invocation at every threshold at both resolutions with the QoIs of item 2 within tolerance (the prototype's standing: one invocation at 60 cells for `ε_m` of 3e-3 and 1e-2, two at 120 cells and at 3e-2, A.3 and A.4), M5 to M7 proceed as planned. If not, the owner chooses between tuning (`τ`, the band, the activation) within M4's allowance and the fallback of item 11, which adds the slot-shifting interior merge (about one and a half weeks, unverified, C54) before M5. Under the fallback M5 survives unchanged (event location applies to any deactivation), M6's three cases assert the documented per-layer invocation count of item 1's fallback clause instead of one, V-BO-07's `β` measurement is dropped, and the gradient re-measurement runs through the many events as Phase 0 did.

---

## 8. Risk register

Carried from design §5 and plan §8 (R20 to R22 are Phase 2's) and extended with what this session measured.

| # | Risk | Where it shows | Mitigation / owner |
|---|---|---|---|
| R20 | Event ordering under a parameter perturbation makes finite-difference gradients non-smooth | M6, item 9 | With one invocation per layer there is one event time to differentiate; measured by the Phase 0 protocol; the bracketed differencing of design §5 if it persists |
| R21 | `ε_m` small enough for one event per layer is too small for the Newton | M4 | Measured: 3e-3 gives one invocation and 0.31 µm cells with no stall (A.4); `Rodas5P` has no Newton; the domain test and the step counts are the gate quantities |
| R22 | The correction fails to carry a layer, or needs problem-dependent tuning (Furzeland) | M4, §7.4 | The prototype carries V-BO-01 at `τ = 0.2` with a band of 10 `ε_m` in one invocation at 60 cells for `ε_m` of 3e-3 and 1e-2, in two at 120 cells and at 3e-2 (A.3, A.4); the absolute band is unmeasured at the other thresholds; V-BO-06/07/08 test the tuning; the fallback is item 11 |
| R25 | **An end-of-layer MLR spike.** After the band merge the last per cent of the layer gasifies in a burst: at 30 PMMA cells the peak reads 36.7 against the old core's 32.25 at that resolution (+14 %), and at 60 cells with `τ_p = 0.1` the peak moves to the merge time (A.4); at 60 and 120 cells with `τ_p = 0.2` there is no spike. The mechanism from the rule as written: with first-order kinetics the merge changes the total rate `Σ m_i k(T_i) s_i` only through the switch and the temperature, and it changes both: the merged residue was throttled (`s` 0.74 to 1 inside the band) and the substrate receiver's switch is one, releasing up to 1.35×; the thin film transmits most of the incident flux into a receiver that absorbs it in a fraction of a micrometre, so the merged temperature is above the film's, and 5 K is +27 % in the rate at 650 K; the 30-cell receiver has twice the heat capacity and pulls the merged temperature further toward the substrate's, which is why it spikes more than `τ_p = 0.1` at 60 cells despite the tighter throttling there. The peak QoI is a point sample of a 2 to 3 s burst on a 0.5 s output grid | M4, item 2 | The tail diagnostic of M4: at the last pre-merge step record per cell `m_i`, `T_i`, `s_i`, at the post-event state the residue and `T_m`, decompose `Σ m_i k(T_i) s_i → m_res k(T̄) s̄ → m_res k(T_m)`, compare the burst's integrated mass with the merged residue, and repeat with the switch off (a surviving spike is the temperature part). Merging into the lowest layer cell as the survivor removes only the switch part and trades a second invocation against item 1; the structural options are C51's. The owner decides the switch (question 12) and the 30-cell re-pin (question 8) on that table |
| R26 | Wall time 3.5 to 6× the old core's at equal tolerance, more at matched accuracy (C49, item 8); the owner's rule is the same order or better | M2, item 8 | The levers of §5.7 in their measured order, the integrator diagnosis first; a stop at M2's gate above 2× (§7.4) |
| R27 | `τ` is a time and the burnout time scale varies by material and flux (Furzeland) | M4, M6 | The `N²` form makes `τ` the layer's equidistribution time, mesh-independent (A.10); the ladder's `τ` row; V-BO-07 and V-BO-08 measure `T_sync` against the layer's life; a `setup` keyword |
| R28 | The band semantics: a multiple of `ε_m` takes cells at 29 % of their mass when `ε_m = 3e-2` (A.4) | M3 | The band is an absolute mass fraction (§4.3) |
| R29 | The correction's stiffness raises the step count and the rejections (1,091 steps at `τ = 0.2`, 1,286 at 0.1, against 853 on the old core) | M4 | `τ` default from the ladder; C52; the restart policy (A.6) |
| R30 | The widened pattern costs 1.7× per Jacobian (30 against 18 colors) | M1 | Measured and accepted; a narrower stencil would be first-order and kinked (C21) and is not pursued |
| R31 | Trial-state gaps on the relative-flux path (C47: a `NaN` face temperature throws) | M1, M3 | C47 fixed; §6.1 item 12 drives the kernel with `NaN`, zero and negative states |
| R32 | A deck's layer interface is not declared, so the correction mixes two materials across it | M3 | The explicit `lagrangian` list is required whenever a correction is selected on a deck with more than one condensed composition; `setup` refuses the combination otherwise |
| R33 | The fallback needs the slot-shifting merge, which is a cross-cutting contract (C54) | §7.4 | Scheduled only if the gate fails; its cost estimated at one and a half weeks |
| R34 | The three new cases have no old-core reference for peak and timing (the old core parks interior cells and may differ at the grading) | M6 | Compared against the old core at the same mesh in the same session; a difference is reported with both values, not pinned silently |
| R35 | The enthalpy form at 3 to 4× the temperature form's time on the burnout decks (A.4) | M4 | The temperature form stays the default (C50); the enthalpy form runs the gate only |
| R36 | The Douglas fir deck at 49 colors is 1.9× the old core at equal tolerance (A.8) | Phase 3 | Recorded; the matched-accuracy comparison belongs to the cutover table (question 9) |
| R37 | A re-pinned QoI hides a regression | M3, M4, M6 | Every re-pin lists the old value, the new value and the reason (C45, C46, R25), and is the owner's |
| R38 | The located event is applied at an interpolated state of the method's order | M5 | The safety net; the admissibility check; item 7 measures the error against the tolerance with location on |
| R39 | Two prototype runs in one session once showed no events where a clean session shows two invocations (A.3's note) | M1 | The runs were discarded; the package tests of §6.1 item 16 (repeated solves bit-identical with events) guard the package implementation |

---

## 9. Branch and process

### 9.1 Branches and the `main` rule

`mesh-framework-phase2` is branched off `mesh-framework-phase1` at `bb2860b` and returns to it by pull request; this plan opens the draft pull request. **Nothing goes to `main` until the owner commits to this architecture at the cutover decision:** no push to `main`, no change to pull request #12 (`mesh-framework-phase1` → `main`), which the owner alone handles. This supersedes plan §9.1 and C40's "the draft pull request becomes a real one": #12 stays as the owner left it, and Phase 2's work is reviewed against `mesh-framework-phase1`. Every commit is GPG-signed, one theme per commit, the subject imperative, the body stating what, why and the evidence; no AI attribution of any kind anywhere (commits, pull requests, issues, files); the draft pull request's body is this plan's summary (§11) and its owner questions (§13). The Quality workflow (`.github/workflows/quality.yml`) triggers only on `main`, so it never runs on the draft pull request: the format check (`julia --project=scripts scripts/format.jl --check`), Aqua and JET are run locally before every push, as plan §9.3 does, and the workflow can be run by hand on the branch (`workflow_dispatch`).

### 9.2 The prototype as evidence

`scripts/mesh_phase2_prototype/` holds five files: `common.jl` (the method overrides, in a dedicated Julia session only, never in a session that runs the package's tests), `run.jl` (the prototype's experiments by name: `check`, `sweep`, `band`, `repeat`, `ladder`, `order`, `continuous`, `scaling`), `baselines.jl` (the measurements that override nothing: `vbo_old`, `vbo_new`, `vbo_ladder_old`, `vbo01_qois_old`, `vbo02_restart`, `df_cores`, `integrators`, `tolweights`), `cases.jl` (the case runner the two `vbo_*` baselines use) and `README.md` (what is overridden, how to run, that it is planning evidence and not package code). Every script prints its table to standard output; no measurement file is written, and the tables are Appendix A, so the "scripts and tables" commit holds the scripts alone. The directory is formatted with the repository's formatter and committed in its own commit before the plan. `common.jl` and `run.jl` are deleted when M3 lands, since M3 reproduces the one-invocation run inside the package and M2 measures on the prototype's run (§5.7, §7); `baselines.jl` and `cases.jl` stay until their measurements have gate scripts of their own (§6.3).

### 9.3 Commit and review discipline

As plan §9.3: one theme per commit; replaced rules deleted; no dates or historical notes in code or docstrings; the guides untouched; `julia --project=scripts scripts/format.jl` before every commit (the formatter covers `scripts/`); Aqua and JET after any change to exports, dependencies or signatures; targeted gates per milestone and the full three-way `Pkg.test()` at M7; the old-core fingerprint and the old core's `Pkg.test()` after any edit of Materials or BoundaryConditions; gate tables stamped with their commit. The development loop of Phase 1 (a scratch environment mirroring the root manifest, Revise job servers, never two processes precompiling Pyrolysis at once, measurements inside functions, both cores in one process, timings repeated; recorded in the owner's local Phase 1 handoff) applies unchanged, with one addition from this planning: a server in which a job has overridden package methods is retired before any package test runs in it, and a job that fails while defining an override leaves the server unusable (R39).

---

## 10. Decisions carried as gated choices

| Decision | Phase 2 default | Alternative | Decided by | When |
|---|---|---|---|---|
| The correction | `EquidistributionCorrection` with the layer time `τ` of §3.2, `α = 1`, `β = 0`, activation 0.3/0.1, `ε_band = 0.1`, selected per solve by the `correction` keyword (the burnout cases pass it; `setup`'s default is `NoCorrection`, since a correction object present but inactive still widens the pattern to 30 colors, C43, C56) | `c = 0` plus the block-and-band rule and the slot-shifting interior merge (item 11); or the correction as `setup`'s default, decided at the cutover | §7.3 items 1, 2 and 8 on the matrix | M4 (§7.4), M6 |
| `τ` (layer time, `N²` form) | set by the ladder at constant `τ` (the prototype's 720 s at 60 cells) | the prototype's per-cell form | the `τ` row of item 1, A.10 and R25's tail diagnostic | M4 |
| The band | 0.1 (absolute) | 0.03 to 0.3 | item 1's merged-band fractions | M4 |
| `ε_m` | 1e-2 | 3e-3, 3e-2 | item 1 across thresholds; the Phase 1 V-CPL-02 numbers | M4 |
| The monitor's mass term `β` (normalized) | 0 | 1 | V-BO-07 (mixed composition, the only deck where it acts) | M6 |
| The reactivity switch | on, 1.5 `ε_m` (C27) | off | item 2 and the tail diagnostic | M4 |
| Event location | off | on | item 7 | M5 |
| Restart after an event | solver-chosen | forced `DT_AFTER_EVENT` | item 7's ledger quantities against the 30 % cost (A.6) | M5 |
| The face temperature of the relative energy flux | upwind-biased (C33) | centred | only if a burnout QoI's convergence needs it | M4 |
| Energy state for `c ≠ 0` | temperature | enthalpy | item 3 | M4 |
| The wall-time clause | the same order as the old core: ratio ≤ 1.5 at matched accuracy | 1.5 to 2 the owner's call; above 2 the phase stops | item 8 at M2 and M4 | M2, M4 |
| The burnout path's integrator | `Rodas5P` with the exact Jacobian every step | a Newton-based integrator with Jacobian reuse once its behaviour on the new residual is understood, or a Rosenbrock-W method with a reuse policy | M2's diagnosis and the cost breakdown of §5.7 | M2 |

---

## 11. Summary

**What Phase 2 builds.** The correction `c` of design §3.3 on the new core: a per-node relaxed equidistribution of a monitor with a smooth per-node activation in the cells' thickness ratio, carried by the balance rows as a sweep term so that `w − c` is the solid velocity, with the existing relative-flux kernel reading `−c` as the relative velocity and every inventory, the enthalpy and the history pair crossing faces with the same fluxes; the Jacobian pattern widened for it and verified by containment; the cascade rule that the design's one-callback-per-step needs, a block of near-threshold cells from the exposed side merging into the cell beneath in one invocation; event location by a continuous callback with the discrete predicate as its safety net; the constraint scale, the outputs, the surface Newton's trial clamp and a trial-state gap closed. The old core is untouched except for two additive, fingerprint-checked edits (the clamp keyword, an optional material field). V-BO-01, 03, 04 and 05 join the new core with re-pins the owner approves; three burnout cases are added (graded, heterogeneous, localized).

**What the evidence says.** A scratch prototype of exactly this design on the Phase 1 code carries the V-BO-01 PMMA layer through burnout in one deactivation at 60 PMMA cells for `ε_m` of 1e-2 and 3e-3 with every quantity of interest within its tolerance (peak within 0.01 %, burnout time exact, ledger 7e-5), in two invocations at 120 cells and at `ε_m = 3e-2`, with cells shrinking to a micrometre and no stall; the pieces that made it work and were not in the design's text are the per-node "or" activation and the near-threshold band, both recorded as conflicts. The cost is 3.5 to 6× the old core's at equal tolerance, mostly the per-step Jacobian at 30 colors, and more at matched accuracy; no integrator candidate closes Phase 1's accuracy-per-tolerance gap; event location ties the event times to the tolerance for 2 to 4 % of the solve; a forced restart buys the old ledger tolerance back for 30 %. At 30 cells the end of the layer shows an MLR spike that the ladder must settle.

**What it costs.** Eleven and a half to fourteen and a half weeks of one developer, the review and a performance milestone included (the design said 4 to 6); thirteen to sixteen if the correction fails its gate and the slot-shifting merge is built; seventeen new tests, three new cases, four new gate scripts and two extended; no change to the old core beyond the clamp keyword (default unchanged) and an optional material field.

**What its gate proves.** One invocation per layer across grids, thresholds and monitor settings, with a graded grid, a heterogeneous layer and a localized layer among the checks; the peak and its timing against the old core; closed ledgers; admissibility at every event; the charring decks and the Douglas fir deck bit-identical without a correction; constraint residuals and domain rejections; event times tied to the tolerance; wall time at matched accuracy against the old core, reported with its spread; smooth gradients through the one event; the minimum cell size.

**What remains open after it.** Whether M2 brings the burnout path within the old core's order (the owner's rule; above 2× the phase stops); the interior merge on a moving mesh at `c = 0` (built only as the fallback); the design's items for later phases (radial Beer–Lambert, P1, the variable cross section's fate, multi-D deactivation by element death); the cutover's end-to-end comparison on the real wood decks, for which the Douglas fir deck now runs on the new core.

---

## 12. Implementation status

Branch `mesh-framework-phase2`, draft pull request #14 against `mesh-framework-phase1`. M0: this plan, the prototype (§9.2) and the owner's answers (§13, every question as recommended). Implementation starts with M1 (§7) and must pass M2's performance gate (§5.7) before the correction is built; its departures continue the series at D53, its conflicts at C57, and this section records each milestone's commits and gate evidence as plan §12.1 did.


---

## 13. Owner's questions

Each with a recommendation and what the answer changes. Questions 1, 2, 4 and the V-BO-04 part of 8 are needed before M1 starts (§7.1); the rest by the milestone named. **Owner: every question decided as recommended** (the interior policy by the block-and-band rule with the slot-shifting merge as the fallback only; the band an absolute fraction of 0.1; `τ` the layer's equidistribution time set by the ladder; thresholds as `setup` keywords with the optional material field; event location on if item 7 passes; the wall-time ratio decided on item 8's table; the solver-chosen restart kept; the V-BO-04, V-BO-05, V-BO-03 and surface-Newton re-pins approved and the 30-cell V-BO-01 row decided after the tail diagnostic; the Douglas fir comparison at the cutover; the accuracy-per-tolerance investigation closed; the temperature form; the switch decided on R25's decomposition; the design errata at the design's next revision; radial Beer–Lambert deferred to Phase 4a; the Péclet weighting not in Phase 2).

1. **The interior policy on a moving mesh.** The block-and-band rule plus the correction is the primary path; the slot-shifting interior merge (C54, the old core's rule) is either insurance built now or the fallback built only if the correction fails its gate. *Recommendation:* the fallback only; the prototype carries V-BO-01, 03's deck and 04 without it, and it is a cross-cutting contract (every per-slot constant shifts) worth one and a half weeks. *Changes:* built now, M3 grows by one and a half weeks and the fallback of item 11 is already in hand; later, the schedule risk of §7.4 stays.
2. **The band.** An absolute near-threshold mass fraction `ε_band` (default 0.1) or a multiple of `ε_m`. *Recommendation:* absolute, 0.1, set by the ladder (a multiple took cells at 29 % of their mass at `ε_m = 3e-2`, A.4; three thresholds leave cells behind at 60 cells, A.3). *Changes:* the parameter's name and default in `setup` and in the design errata (§7.2).
3. **`τ`.** The layer's equidistribution time with the `N²` factor in the formula (§3.2), or the prototype's per-cell form. *Recommendation:* the layer time, a `setup` keyword whose default the ladder sets at constant `τ` across 30, 60 and 120 cells (the prototype's one-invocation run is `τ = 720` s at 60 cells, `T_sync` 73 s), with the heuristic that `T_sync = τ/π²` is at most a tenth of the layer's remaining life after full activation; the per-cell form makes `τ` mesh-dependent and is what produced the ladder's rows. *Changes:* R27 reduces to one physical ratio, `T_sync` against the layer's life, measured on V-BO-07 and V-BO-08.
4. **Thresholds as material parameters** (design question 9). `setup` keywords with package defaults plus an optional `Material` residual-thickness fraction (additive to the old core, fingerprint-checked), or package constants. *Recommendation:* the keywords and the optional field, since a char former that keeps 20 % of its thickness must not activate and only the material knows that. *Changes:* with the field, M1 edits Materials and runs the fingerprint and the old core's `Pkg.test()`; without it, charring decks rely on the activation never reaching 30 % (true of the current decks) and a new material could activate wrongly.
5. **Event location by default.** On (2 to 4 % of the solve; the event time tied to the tolerance, A.5) or off. *Recommendation:* on, if §7.3 item 7 passes with it; the discrete predicate stays as the safety net either way. *Changes:* the default of `located_events` and whether item 7's budget is met by the deactivating cases.
6. **The wall-time clause.** As written (no worse than the old core) or an exception for Phase 2. *Owner's decision:* the burnout path must run within the same order as the old core or faster; a path that cannot is a bad sign for the architecture. Written into the plan as milestone M2 (§5.7, §7): the integrator diagnosis and the Jacobian levers first, measured at matched accuracy with a cost breakdown, a pass at 1.5×, the owner's call between 1.5 and 2×, and a stop above 2× for the owner's decision on the phase and the cutover. *Changes:* M2 precedes the correction's build; the estimate grows by two to three weeks; the rest of the phase depends on M2's gate.
7. **The restart after an event** on deactivating solves. Solver-chosen (the default; V-BO-02's ledger quantity 1.5e-3 against C39's 3e-3) or forced `DT_AFTER_EVENT` (30 % more time, ledger 7.6e-5, A.6). *Recommendation:* keep solver-chosen and C39's tolerance; measure both on V-BO-01 with the correction at M5, and switch only if a QoI needs it. *Changes:* V-BO-02's ledger tolerance (back to 1e-3 under the forced restart) and the step budgets.
8. **Re-pins.** V-BO-04 per C46 (one invocation, one survivor); V-BO-05 per C45 (`:interior` at the first `ε_m` crossing, the time from a ladder); V-BO-03's warning QoI becomes the transmitted series' peak; the surface-Newton warning count of V-BO-01/03 becomes a count of clamped trial iterates once C36 removes the warnings, or is dropped; V-BO-01 at 30 PMMA cells if the end-of-layer spike (R25) persists after M4's diagnostic. *Recommendation:* approve the first three now (M3 needs V-BO-04's); decide V-BO-01's after the diagnostic, with the mechanism named. *Changes:* which QoIs the gate asserts and which it reports.
9. **The Douglas fir deck.** The both-cores comparison in Phase 2 or in Phase 3. *Recommendation:* Phase 2 runs it as the "charring deck identical" check (bit-identical with `NoCorrection`) and records both cores' times (A.8: 1.9× at equal tolerance, 49 colors); the matched-accuracy comparison on the real wood decks is the cutover's end-to-end table (design §6 Phase 3 item 6). *Changes:* whether M6 carries a matched-accuracy run of a 600 s deck (about a day of measurement) or Phase 3 does.
10. **The accuracy-per-tolerance investigation** carried from gate item 2. Close it with A.7 and A.9 (no candidate meets the old core's error within 1.5× of its time on the moving-mesh deck; the tolerance weights are not a lever; `Rodas5Pe` is the closest on the Darcy and Eulerian decks), or spend more. *Recommendation:* close it and carry the exception to Phase 3, where the cutover decision weighs it. *Changes:* a week of M5 if the owner wants the Rosenbrock error estimators studied further.
11. **The energy form for `c ≠ 0`.** Temperature (the gate's default) or the inventory form the design named for this work. *Recommendation:* temperature; the prototype's enthalpy run reproduces every QoI at 3 to 4× the time (A.4), and the ledgers close in both. *Changes:* which form the V-BO cases run in by default; the ladder runs both regardless.
12. **The reactivity switch** (C27). Keep at 1.5 `ε_m`, drop it under an `EquidistributionCorrection`, or trigger the band at its onset. *Recommendation:* decide on R25's decomposition measurement at M4; the hypothesis is that the switch is the release part of the end-of-layer spike and is no longer needed once the band deactivates cells at or above `ε_m`. *Changes:* a `setup` keyword either way; with the switch off C27's V-CPL-02 numbers are re-run.
13. **The design errata of §7.2.** Apply them to the design document in M0 (a design revision commit) or at the design's next revision. *Recommendation:* at the next revision, listed here, since the owner owns the design text.
14. **Radial Beer–Lambert** (C37, carried open). *Recommendation:* deferred to Phase 4a, where the first radial case needs it; no burnout case is radial and `setup` refuses the combination (D45). *Changes:* nothing in Phase 2.
15. **The Péclet weighting of the Darcy upwind** (C21). *Recommendation:* not in Phase 2: the default integrator has no Newton to fail at the kink, and no Phase 2 deck reverses its flow. *Changes:* a week of M6 if the owner wants it measured on the audit decks now.

---

## Appendix A. Measurements

All on `bb2860b`, Julia 1.13.0, one thread per process, the machine shared (load 11 of 22 logical CPUs throughout, lower than Phase 1's 17 to 24), warm times unless stated; the two cores always in one process back to back. The scripts are in `scripts/mesh_phase2_prototype/` (§9.2): `baselines.jl <name>` for A.1, A.2 and A.6 to A.9, `run.jl <experiment>` for A.3 to A.5. The prototype's overrides persist in a session, so each `run.jl` experiment is its own process; the A.3 and A.4 tables below come from the clean sequences `sweep`, `band`, `repeat`, `ladder` and `order` (an earlier sequence in a session whose first job had failed while defining an override gave two zero-event runs that no clean session reproduces; they are discarded, R39).

### A.1 The burnout cases on both cores (`vbo_old`, `vbo_new`; CI resolutions, the cases' tolerances, median of three)

| case | core | retcode, outcome | t_end [s] | wall [s] | steps (rejected) | Jacobians | events or parks | QoIs |
|---|---|---|---|---|---|---|---|---|
| V-BO-01 | old | Success | 500 | 0.80 | 853 (3) | 164 | 60 parks | 7/7 |
| V-BO-01 | new | Terminated, `:interior` | 241.6 | 0.45 | 333 (144) | 333 | 0 | 2/7 (ledger, surface warnings) |
| V-BO-02 | old | Success | 500 | 1.29 | 1,572 (0), 56 Newton failures | 341 | 60 | 10/10 |
| V-BO-02 | new | Success | 500 | 1.32 | 1,706 (953) | 1,706 | 60 | 10/10 (ledger 1.5e-3 of 3e-3; peak 23.67, informational) |
| V-BO-03 | old | Success | 500 | 0.58 | 681 (3) | 126 | 48 | 4/4 |
| V-BO-03 | new | Terminated, `:interior` | 241.6 | 0.44 | 333 (144) | 333 | 0 | 1/4 |
| V-BO-04 | old | Success | 480 | 0.034 | 353 (0) | 77 | 35 | 4/4 |
| V-BO-04 | new | Terminated, `:interior` | 468.3 | 0.039 | 304 (2) | 304 | 1 | 1/4 |
| V-BO-05 | old | Terminated (as pinned) | 337.0 | 1.42 | 169 (9) | 86 | 0 | 4/4 |
| V-BO-05 | new | Terminated, `:interior` | 263.9 | 0.34 | 348 (143) | 348 | 0 | 2/4 (stop time, warning) |

The new core's V-BO-01 ends with cells 179 and 180 below the threshold under live cells 181 and 182 (the absorption peak depletes a sub-surface cell first); its peak MLR to that point is 23.58 at 241.6 s.

### A.2 The old core's burnout ladder (`vbo_ladder_old`, `vbo01_qois_old`; median of three)

| case | PMMA cells (total) | wall [s] | steps | Newton failures | parks | V-BO-01 peak [g m⁻² s⁻¹], t_peak, t_bo [s] |
|---|---|---|---|---|---|---|
| V-BO-01 | 30 (91) | 0.25 | 648 | 6 | 30 | 32.254, 333.0, 364.0 (mass loss 69.9206 %) |
| V-BO-01 | 60 (182) | 0.58 | 853 | 8 | 60 | 32.392, 335.0, 363.5 (69.9206 %) |
| V-BO-01 | 120 (364) | 2.11 | 1,461 | 62 | 120 | 32.469, 335.5, 363.5 (69.9206 %) |
| V-BO-02 | 30 / 60 / 120 | 0.45 / 1.11 / 3.70 | 1,013 / 1,572 / 4,198 | 45 / 56 / 102 | 30 / 60 / 120 | — |
| V-BO-03 | 30 / 60 / 120 | 0.21 / 0.54 / 1.79 | 515 / 681 / 1,259 | 2 / 2 / 62 | 24 / 48 / 96 | — |

V-BO-01 at 60 cells took 0.80 s in A.1 and 0.58 s here, the machine's spread.

### A.3 The prototype on V-BO-01 at 60 PMMA cells (`sweep`, `band`, `repeat`, `order`; in-depth radiation, temperature form, reltol 1e-6, abstol 1e-8, `ε_m = 1e-2`, `min_thickness = 1e-3` m)

The prototype: `c_n = a_n (h_n − h_{n−1})/τ_p` (`α = 1`, `β = 0`; the per-cell form without the `N²` factor of §3.2, so its `τ_p` is the layer time divided by `N²`), the "or" activation between 30 % and 10 % of the initial thickness, `c = 0` at the pinned, surface and PMMA–Kaowool interface nodes; the balance rows and the consistent sweep with the sweep term; `RelativeFluxKernel` at `u = −c·n`; the widened pattern (contained in the dense Jacobian at 12 cells, colored equal to dense to round-off; the first check saturated the activation and used `β = 0`, so the added mass columns were trivially contained, see A.10); the block-and-band rule. Wall is the smaller of two warm solves (the `repeat` runs: three solves, identical events and outcomes each time).

| configuration | invocations (cells; largest merged mass fraction) | outcome, t_end [s] | peak, t_peak, t_bo | mass loss [%] | ledger | wall [s] | steps (rejected; by the domain test) | largest constraint residual | smallest cell [m] | colors |
|---|---|---|---|---|---|---|---|---|---|---|
| `c = 0`, block rule | 0 | `:interior`, 241.6 | 23.58, 241.6, — | 35.5 | 4.5e-7 | 0.40 | 333 (144; 0) | 2.8e-3 | 9.6e-7 | 18 |
| `τ = 1`, band `ε_m` | 3 (2, 4, 6) | `:interior`, 355.9 (cells 161 to 164 below) | 32.371, 334.5, — | 67.7 | 2.8e-6 | 1.26 | 590 (180; 5) | 0.18 | 9.7e-7 | 30 |
| `τ = 0.2`, band `ε_m` | 1 (2) | `:interior`, 359.5 (172, 173) | 32.391, 335.0, — | 68.7 | 4.1e-7 | 2.03 | 1,007 (213; 52) | 0.69 | 9.8e-7 | 30 |
| `τ = 5`, band `ε_m` | 1 (2) | `:interior`, 276.4 (175, 176) | 26.74, 276.4, — | 44.0 | 5.1e-4 | 0.79 | 353 (143; 0) | 4.8e-2 | 9.7e-7 | 30 |
| `τ = 0.2`, band 3 `ε_m` | 1 (41; 3.0 %) | `:interior`, 362.4 (cell 123) | 32.391, 335.0, — | 69.4 | 3.6e-6 | 2.35 to 2.75 | 1,047 (223; 54) | 0.69 | 9.8e-7 | 30 |
| `τ = 1`, band 3 `ε_m` | 2 (50) | `:interior`, 363.1 | — | — | — | 1.41 | 621 (203; 6) | 0.22 | 9.8e-7 | 30 |
| `τ = 1`, band 10 `ε_m` | 2 (34, 9.2 %; 26, 1.3 %) | Success, 500 | 32.392, 335.0, 363.5 | 69.921 | 2.0e-4 | 1.57 to 2.33 | 691 (200; 10) | 0.19 | 9.8e-7 | 30 |
| `τ = 0.2`, band 10 `ε_m` | **1 (60; 3.6 %)** | Success, 500 | **32.391, 335.0, 363.5** | 69.921 | 7.2e-5 | 2.76 to 3.5 | 1,091 (246; 52) | 0.69 | 1.03e-6 | 30 |

Reference values: the case's fixed references 32.39, 335.0, 363.5 and the exact mass loss 69.915 % (`VBO_MASS_LOSS_EXACT`); the old core at 60 cells 32.392, 335.0, 363.5 and 69.9206 % (A.2, `vbo01_qois_old`); the case's tolerances 3 %, 3 %, 1 %, 0.3 percentage points, ledger 1e-3. The walls here are the smaller of two solves in the prototype's process; the old core's are A.2's medians from another process. The `τ = 1` run's first invocation at 302.7 s takes 34 cells of which the largest holds 9.2 % of its mass: the band of `κ = 10` is what makes it one invocation, and the surface cell crossing first (C48) is what the band absorbs.

### A.4 The prototype across grids, thresholds, monitor, activation and energy form (`ladder`; `τ = 0.2`, band 10 `ε_m` unless stated)

| run | invocations (cells; largest fraction) | outcome | peak, t_peak, t_bo | wall [s] (old core, same resolution) | steps (rejected; domain) | largest residual | smallest cell [m] |
|---|---|---|---|---|---|---|---|
| 120 PMMA cells | 2 (73, 10.0 %; 47, 1.7 %) | Success | 32.469, 335.5, 363.5 (old 32.469, 335.5, 363.5) | 9.3 (2.11) | 1,303 (460; 46) | 0.90 | 4.9e-7 |
| 30 PMMA cells | 1 (30; 1.75 %) | Success | **36.676, 361.0**, 363.5 (old 32.254, 333.0, 364.0) | 1.46 (0.25) | 867 (254; 41) | 0.60 | 2.0e-6 |
| 60 cells, `τ = 0.1` | 1 (60; 1.2 %) | Success | **32.987, 361.5**, 363.5 | 4.2 | 1,286 (253; 73) | 0.92 | 1.0e-6 |
| 60 cells, `β = 1` (mass term) | 1 (60; 1.2 %) | Success | 32.987, 361.5, 363.5 (identical to the row above: the term is inert at constant density, and the run differs from the `β = 0` run only through its `τ`-equivalent) | 4.4 | 1,301 (259; 77) | 0.88 | 1.0e-6 |
| 60 cells, "and" activation (C44) | 1 (56; 9.8 %) | `:interior` at 362.4 (cell 123) | 32.391, 335.0, — | 4.6 | 1,416 (271; 45) | 0.65 | 9.8e-7 |
| 60 cells, `ε_m = 3e-3` | 1 (60; 0.6 %) | Success | 32.391, 335.0, 363.5 | 3.6 | 1,169 (202; 70) | 0.97 | 3.1e-7 |
| 60 cells, `ε_m = 3e-2` | 2 (47, 28.7 %; 13, 4.7 %) | Success | 32.542, 333.0, 363.0 | 3.1 | 984 (229; 46) | 0.12 | 3.0e-6 |
| 60 cells, activation 0.5/0.2 | 1 (60; 1.85 %) | Success | 32.391, 335.0, 363.5 | 3.3 | 1,102 (221; 60) | 0.46 | 1.0e-6 |
| 60 cells, enthalpy form | 1 (60; 3.6 %) | Success | 32.391, 335.0, 363.5 | 10.9 | 1,585 (636; 55) | 0.45 | 1.0e-6 |

Every Success run loses 69.921 % of the mass and ends at 11.997 mm with a ledger between 7e-5 and 2.5e-4; three of them (30 cells, `τ = 0.1`, `β = 1`) fail the case's `t_peak` tolerance (361.0 and 361.5 s against 335.0, +8 %) and the 30-cell run its peak tolerance. The two bold rows are R25: the peak moves to the merge time with a spike of +14 % (30 cells) and +1.8 % (`τ = 0.1`). The `β = 1` row's identity with the `τ = 0.1` row is the unnormalized monitor (`M = 1 + ρ/ρ_ref = 2` on the PMMA layer doubles the correction, which is `τ_p = 0.1`): `β` is a hidden `τ` until the monitor is normalized, and does nothing on a uniform-density layer, as §3.2 says. The `ε_m = 3e-2` row shows why the band is an absolute fraction (R28). The 120-cell run at `τ = 1`, band 3 `ε_m`, threw `InexactError` in the tabulated enthalpy at a non-finite face temperature (C47).

### A.5 Event location on V-BO-02 (`continuous`; fixed mesh, 60 events)

| callback | reltol | wall [s] | steps (rejected) | events | first, last event [s] | ledger quantity |
|---|---|---|---|---|---|---|
| discrete (Phase 1) | 1e-6 | 1.19 to 1.33 | 1,706 (953) | 60 | 69.201, 465.254 | 1.50e-3 |
| discrete | 1e-8 | 3.54 to 3.58 | 5,626 (2,276) | 60 | 69.125, 464.649 | 4.1e-5 |
| continuous + discrete safety net | 1e-6 | 1.24 to 1.36 | 1,733 (937) | 60 | 69.035, 464.354 | 1.37e-3 |
| continuous + safety net | 1e-8 | 3.63 | 5,644 (2,278) | 60 | 69.033, 464.362 | 7.0e-6 |

Event-time shift between the two tolerances, largest over the 60 events: discrete 0.88 s (0.133 of the median interval of 6.6 s), continuous 0.15 s (0.022). The experiment ran on the fixed mesh without the correction, so event location was never exercised together with it; M5 does that. Between the discrete and the located events at reltol 1e-8 the largest difference is 0.38 s. Without the safety net and the nudge at the root, the located run at 1e-8 missed 21 of the 60 events (a root left the cell at the threshold, the strict predicate never fired again, and the margin stayed negative). The continuous callback's cost is 2 to 4 % of the solve.

### A.6 The restart policy after an event on V-BO-02 (`vbo02_restart`; new core)

| policy | wall [s] | steps (rejected) | Jacobians | RHS | peak | ledger quantity |
|---|---|---|---|---|---|---|
| `Rodas5P`, solver-chosen (default) | 1.37 | 1,706 (953) | 1,706 | 22,978 | 23.67 | 1.50e-3 |
| `Rodas5P`, forced `DT_AFTER_EVENT = 1e-6` s | 1.77 | 2,790 (286) | 2,790 | 27,398 | 23.71 | 7.6e-5 |
| `KenCarp4` raw, flag kept | 5.14 | 10,021 (43), 1,784 nonlinear failures | 1,025 | 158,084 | — | — |

### A.7 Integrator candidates on the audit decks (`integrators`; new core, 100 cells, solid gas-transfer coefficient 1e-6 m²/s on the PMMA decks (C30), 80 s; error = largest temperature difference against the new core's `Rodas5P` at reltol 1e-10; target = the old core's error at its default against its own reltol 1e-10 reference; time = smaller of two warm solves)

| deck | old core default: error [K], wall [s] | best candidates (reltol): error / target, wall / old | `Rodas5P` (1e-6, 3e-7, 1e-7): error / target, wall / old |
|---|---|---|---|
| ale 100 | 1.03e-5, 0.0163 | `Rodas5P` 3e-7: 1.14, 1.52; `Rodas5Pe` 1e-7: 0.28, 1.77; `Rodas4P2` 3e-7: 1.49, 2.13 | 9.3, 1.19; 1.14, 1.52; 0.63, 1.92 |
| eulerian 100 | 1.37e-5, 0.0129 | `Rodas5Pe` 1e-7: 0.95, 1.26; `Rodas5Pr` 3e-7: 0.80, 1.52; `KenCarp4` raw 1e-7: 0.87, 1.45 | 21.6, 0.85; 2.06, 1.07; 0.71, 1.38 |
| ale_darcy 100 | 5.96e-5, 0.0419 | `Rodas5Pe` 1e-6: 1.30, 1.01; `Rodas4P2` 1e-6: 0.73, 1.30; `Rodas5P` 3e-7: 0.24, 1.38 | 3.9, 1.03; 0.24, 1.38; 0.09, 1.81 |

Every candidate's three rows are in the script's output; the W-methods (`ROS34PW2`, `Rodas23W`) and the lower-order Rosenbrock methods reach the target only at 3 to 5× the old core's time, `KenCarp4` and `Kvaerno5` on the new core's residual miss it by 7 to 36× at reltol 1e-6 on the moving-mesh deck (the old core's smoothed `KenCarp4` meets it there), and `KenCarp47` misses it everywhere. No candidate meets the target within 1.5× of the old core's time on ale 100; `Rodas5Pe` is the closest on the other two decks.

### A.8 The Douglas fir deck on both cores (`df_cores`; `examples/07_douglas_fir_cone.jl`'s material, mesh and boundary conditions, 240 cells, 14 components, 6 reactions, ALE, surface absorption, Darcy–Fick, reltol 1e-7, abstol 1e-9, the first 60 s of the example's 600)

| core | retcode | wall [s] | steps (rejected) | Jacobians | colors | surface temperature at 60 s [K] | thickness at 60 s [m] |
|---|---|---|---|---|---|---|---|
| old (`ApproxSolve`, KLU) | Success | 2.72 | 168 | 47 | — | 764.90 | 0.011520 |
| new (`Rodas5P`, `StructuredKLU`) | Success | 5.10 | 313 (32) | 313 | 49 | 764.92 | 0.011520 |

Largest surface-temperature difference over the 61 outputs 0.052 K; thickness difference 7e-8 m; the new core's largest constraint residual 1.8e-3. The deck runs on the new core without change (liquid moisture, the state-dependent heat of sorption, the Darcy–Fick transport with seven gases); at equal tolerance it is 1.9× the old core, and the color count is 49 against the design's estimate of 42 to 48.

### A.9 The tolerance weights of the algebraic and inventory rows (`tolweights`; new core `Rodas5P`, 100 cells, error against the new core at reltol 1e-10)

Scaling the absolute tolerance of the velocity rows (×100, ×10⁴), the position rows (×100) or the energy rows (×100, ×0.01) changes the error by less than 2 % on the ale and eulerian decks at reltol 1e-6 and 3e-7 (the relative tolerance dominates those rows); loosening the mass rows' absolute tolerance by 100 or 10⁴ raises the error 1.3 to 2.9× (26.8 to 31.9× the target on ale 100 against 9.3× at the baseline). The mass rows' tight absolute tolerance (`reltol · ε_m · m_i(0)`) is therefore doing useful work, and no tolerance weight is a lever for Phase 1's accuracy-per-tolerance gap; the gap is the formulation's (masses and positions as states), not the error weights'.

### A.10 The `N²` scaling of the relaxation and the non-saturated containment check (`scaling`, `check`)

The numerics review derived that the prototype's per-cell form relaxes a layer in `T_sync ≈ τ_p N²/π²`, so one layer time corresponds to `τ_p ∝ 1/N²`; tested by scaling the 60-cell run's `τ_p = 0.2` s to the other grids (band 10 `ε_m`, `ε_m = 1e-2`, in-depth radiation, temperature form):

| run | invocations (cells; largest merged fraction) | outcome | peak, t_peak, t_bo (old core at that grid) | wall [s] | steps (rejected; domain) | smallest cell [m] |
|---|---|---|---|---|---|---|
| 120 cells, `τ_p = 0.05` (= 0.2 · (60/120)²) | **1 (120; 3.4 %)** | Success | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 11.0 | 2,111 (525; 153) | 5.0e-7 |
| 120 cells, `τ_p = 0.1` | 2 (90, 9.8 %; 30, 2.6 %) | Success | 32.455, 335.5, 363.5 | 9.0 | 1,735 (477; 104) | 4.9e-7 |
| 120 cells, `τ_p = 0.2` (A.4) | 2 (73, 10.0 %; 47, 1.7 %) | Success | 32.469, 335.5, 363.5 | 9.3 | 1,303 (460; 46) | 4.9e-7 |
| 30 cells, `τ_p = 0.8` (= 0.2 · (60/30)²) | 1 (30; 3.35 %) | Success | **35.600, 359.0**, 363.0 (32.254, 333.0, 364.0) | 0.86 | 691 (240; 25) | 2.0e-6 |
| 30 cells, `τ_p = 0.2` (A.4) | 1 (30; 1.75 %) | Success | **36.676, 361.0**, 363.5 | 1.46 | 867 (254; 41) | 2.0e-6 |

At one layer time (`τ = 720` s, `T_sync` 73 s) the three grids merge their layer in one invocation with the same largest merged fraction (3.35, 3.6 and 3.4 % at 30, 60 and 120 cells), which is the scaling the formula of §3.2 predicts; the 120-cell runs at `τ_p` of 0.1 and 0.2 are its under-relaxed side (two invocations, the first band saturated at 10 %). The 30-cell spike shrinks from +13.7 % to +9.9 % at the layer time and remains, so over-synchronization is not its cause and R25's mechanism (the merge's release of the switch and the receiver's temperature) stands as the hypothesis to measure. The cost of the layer time at 120 cells is the stiffest mode: 2,111 steps against 1,303, with 153 domain rejections (C52).

The containment check was repeated with a partial activation (`r_hi = 1.5 > 1 > r_lo = 0.5`, so the activation is strictly inside `(0, 1)` at `t = 0` and its derivatives are exercised) and the monitor's mass term on (`β = 1`): no entry of the dense Jacobian lies outside the widened pattern and the colored Jacobian equals the dense one to round-off, at 12 cells on the perturbed state of A.3's check.

## Appendix B. What the reviews changed

The adversarial reviews of this draft (five areas, as the planning brief asks): the findings, what each changed in the plan, and what was verified against the code or a measurement.

Five reviews were asked for (feasibility against the code, the numerics of the correction and the relative fluxes, the gate's measurability, scope and estimates, the repository rules). All five returned and are folded in below.

**Repository rules and deliverable hygiene.** Eight findings, all verified and applied: `baselines.jl` imported `ADTypes`, which `test/Project.toml` does not carry, so the `integrators` baseline could not run under the documented command (fixed: it reads `AutoFiniteDiff` from `Pyrolysis.Integration`); the Quality workflow triggers only on `main` and never on a draft pull request against `mesh-framework-phase1`, so M0's gate and §9.1 now name the local format, Aqua and JET checks; C36 edits shared old-core code (`src/BoundaryConditions/surface_balance.jl`) and is now counted beside C55 in §2.1, §2.2, §7 and §11 with the fingerprint and the old core's tests on M1's gate; §9.2 now lists the five files, the seven experiments and the eight baselines, says the tables are Appendix A and no measurement file is written, and keeps `baselines.jl` and `cases.jl` past M1; the citations of gitignored handoff files were replaced by "the planning brief"; "this planning session" became "during planning". Verified clean by the review: no AI marker, no dates or history in the scripts, formatted and parsing, nothing written into the tree, no scratch paths, the branch and `main` rules as the brief states them, overrides only by `@eval`, the guides untouched.

**Scope and estimates.** Eleven findings, applied as follows. The fallback (C54) had no measurement and an unsupported estimate: the plan now enumerates what the slot-shifting merge touches, states that the brief's scratch measurement of an interior-merge policy was not run, marks the estimate unverified, adds §6.1 item 17 and says what survives of M5 and M6 under the fallback (§7.4). M1's gate depended on M3's band rule: the milestones are reordered (M1 deactivation, M3 the correction) and the reproduction criterion is one invocation with every QoI within tolerance, not round-off agreement (C52 changes the rejected steps). The default correction's Jacobian cost on non-burnout decks was untimed: `setup`'s default is `NoCorrection` with the burnout cases selecting the correction, and §7.3 item 8 times the audit decks and V-CPL-02/03 with an inactive correction. The estimates were not credible as a range: M4 and M6 grew to 2 to 2.5 weeks each with a tuning allowance, M5 shrank to half a week, the review is its own line, the total is 9.5 to 11.5 weeks, and the plan says that Phase 1's per-milestone actuals were never recorded. Question numbers in C54 and C55 were fixed, the pull request body rule added to §9.1, radial Beer–Lambert and the Péclet weighting given owner lines (questions 14 and 15), M0's deliverable no longer includes the owner's answers, and the composition rule for Lagrangian nodes was dropped from Phase 2 as scope creep (an explicit list only; R32 reworded).

**Feasibility against the code.** Seventeen findings; the load-bearing constructions were verified feasible (`ds.c` filled between `derive!` and the kernels with everything the closure needs; the block-and-band rule as a chain of today's `_transfer!` with the static adjacency and the masks; the widened balance-row reach exactly the correction's dependencies). Two could not be built as written and are changed: `PrescribedCorrection` is bit-identical to Phase 1 only if it copies the state's velocity into `c` (the algebraic rows hold only to tolerance between sweeps), so §3.1 and §3.3 now say so and the pattern keeps the cell rows' velocity entries under it; C36 is an edit of the shared `BoundaryConditions` module, now filed there with its default unchanged and `Kernels` passing the keyword. Also applied: `located_events` is an `integrate`/`solve` keyword and every new keyword joins `_SOLVE_KEYWORDS`; the band is an absolute fraction everywhere (§3.4, §3.9); the band rule is gated on a moving mesh, since the prototype never ran it on V-BO-02 and the rule as prototyped was not; `_consistent_motion!`'s round-off scale gains the sweep; C47 covers `enthalpy_slope` too; `test_mesh_motion.jl` is listed as a changed test; the gate scripts were struck from §2.3; the containment configuration is partially activated so the activation's derivatives are exercised; the `_observing` and `rodas5p_policy` pointers were corrected.

**The numerics of the correction and the relative fluxes.** Twelve findings, three major, all verified and applied. The prototype's node formula omits MMPDE5's `1/Δξ²`, so its `τ_p` is a per-cell time whose layer-scale relaxation is `τ_p N²/π²`: the 30/60/120-cell rows of A.4 were that scaling, not tuning; §3.2 now carries the `N²` factor, `τ` is the layer's equidistribution time, and A.10 confirms the scaling (one invocation at 120 cells once `τ_p` is scaled by `(60/120)²`, the same merged fraction on all three grids). The reactivity switch's premise in C51 was false: the band merges cells inside the switch's range and the substrate receiver's switch is one, which with first-order kinetics is one of the two ways the merge can move the MLR (the other the receiver's temperature); C51, R25 and question 12 now carry that mechanism with the decomposition measurement, and A.10's 30-cell run at the layer time (spike +9.9 %, down from +13.7 %) shows over-synchronization is not the cause. The monitor is normalized by its mean over the layer, since unnormalized `α` and `β` are hidden `τ` scalings (the `β = 1` row of A.4), and V-BO-07 is redefined as a mixed-composition layer, the only configuration where the mass term acts. Also applied: the containment check's saturated activation and `β = 0` made the added mass columns trivially contained (repeated non-saturated in A.10; the package test requires it); the GCL derivation under the sweep and the bidiagonal sweep recorded in §3.1; the activation's floor described as value-inert with the restoring behaviour at crossed nodes and its `C¹` class; the two face temperatures of a gas flux, the upwind-side kink at `c = 0` and the untested history pair in §3.3; one inversion per invocation instead of a chain, and the redistribution after a merge inside the layer, in §3.4; the interpolant's order, the first-crossing semantics under a non-monotone mass and the replacement of the nudge in §3.5; the unscaled residual reported beside C52's in §3.7; the energy columns dropped from C43's reach.

**The gate's measurability.** Nineteen findings, most major, applied: item 1's band criterion was tautological and gameable (the band is now fixed before the ladder, the largest merged fraction bounded at 5 `ε_m`, cells below `ε_m` and cells taken by the band recorded per invocation, and the rows the prototype already fails at the strict criterion listed: 120 cells, `ε_m = 3e-2`, `τ = 1`); the monitor's two settings are carried by V-BO-07, since `β = 1` on V-BO-01 is a `τ/2` duplicate that fails the timing QoI; item 2 judges every run in the case's assert mode against its fixed references (the 30-cell row included unless the owner excludes it) and adds an MLR-curve criterion for the merge burst; item 3 holds V-BO-02 at C39's 3e-3 and names V-BO-04, the window and the event terms for the energy ledger; item 5 records the Phase 1 digests at M0 on `bb2860b`, defines the 1e-10 check on the same step sequence, records the largest activation and smallest thickness ratio, and names the charring decks; item 6 reports both scales and the rejections before and after C52; item 7's budget, unreachable by three orders of magnitude, is restated (located shift at most 0.05 of the median interval and at least three times smaller than the discrete one); item 8 quotes 3.5 to 6× at equal tolerance, states that matched accuracy raises the ratio and that no lever reaches 2×, defines the error measure for decks with events, requires both cores in one process, and times the audit decks with an inactive correction; item 9 passes on the event-timed quantities and the invocation count; item 10 bounds the smallest cell at `0.5 ε_m h₀`; item 11 states the fallback's acceptance and its departure from the design's "nothing else changes"; the three new cases received builder changes, data and references; C53 records the band's reach per invocation; the surface-Newton warning QoI is re-pinned (question 8); §7.4 judges at 60 and 120 cells; C49's "246 of 1,091", §1.5's single walls, the mass-loss reference, the 3.75× of A.6 and the 0.052 K of A.8 were corrected.

