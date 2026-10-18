# Mesh framework: Phase 1 implementation plan

**Scope.** The implementation plan for Phase 1 of the mesh-framework design (`docs/Mesh_Framework_Design_2026-10-17.md` §6, "New 1D core beside `Mesh1D`"), reviewed against `main` at 17bc477. Implementation is under way on branch `mesh-framework-phase1`; §12 records what has landed, the departures from this plan and the conflicts found while building. The design document is the standing design; the harness is its executable contract set; Phase 0's numbers are measurements on proxies wherever the design says so. The four provisional choices (energy state, integrator with its policies, lifting against bordering, the correction) are not resolved here; §10 carries each with the measurement that decides it.

**Conventions.** `file:line` refers to 17bc477. "§3.x" without a document name is a section of the design document. "Harness" is `scripts/mesh_phase0/`. Case ids (V-HC-03, …) are the verification cases under `test/verification/cases/`.

**Contents.** 1 What was verified · 2 Module layout · 3 Data structures and contracts · 4 Kernel-by-kernel plan · 5 Solver plan · 6 Test and verification plan · 7 Milestones, gates, estimates · 8 Risk register · 9 Branch and process · 10 Decisions carried as gated choices · 11 Summary · 12 Implementation status

---

## 1. What was verified

### 1.1 The two checks the plan depends on

**`initialize_state(problem)` is the canonical setup path and returns `(def, ws, u0)`.** Confirmed. `src/Problem/adapter.jl:287-303` defines it as `apply_initial_conditions!` → `to_problem_def` → `build_workspace` → `create_state_vector`, with keywords `use_ale`, `radiation_model`, `radiation_intensity`, `radiation_kirchhoff`; `solve` calls it at `src/Solver/solve.jl:627`; it is exported at `src/Pyrolysis.jl:173`; `test/solver/test_initialize_state.jl` pins the contract; the verification helpers (`test/helpers.jl`, `test/verification/exact/mms.jl:412`) and every Jacobian fixture use it. The new core's adapter (§2.4) reuses `apply_initial_conditions!` and `PyrolysisProblem` and replaces the other three steps.

**The derivative-discontinuity flag is set before every `DiscreteCallback` `affect!`.** Confirmed for the pinned DiffEqBase 7.21.2 (`~/.julia/packages/DiffEqBase/wJ4a3/src/callbacks.jl`): `apply_discrete_callback!` sets `integrator.derivative_discontinuity = true` immediately before `callback.affect!(integrator)` (line 772), captures the value afterwards (779) and, when it is still true, calls `reeval_internals_due_to_modification!` (which on a DAE re-runs `initialize_dae!` with the callback's `initializealg`, `OrdinaryDiffEqCore/8VHiN/src/integrators/integrator_interface.jl:66-91`). Continuous callbacks do the same (line 699). Consequences for the new core, all confirmed by the Phase 0 measurement (`results/integrator_flag_summary.md`): a per-step callback that calls nothing leaves the flag set (the `flag_off` variant is bit-identical to `current`); an observer must call `derivative_discontinuity!(integrator, false)` explicitly; `u_modified!` is a deprecated alias of the same function (SciMLBase 3.56.1 `integrator_interface.jl:341-345`). The package's own observers already do this (`src/Solver/solve.jl:283, 785`, `src/Solver/callbacks.jl:604`, `src/Diagnostics/callback.jl:147`); the step-update callback deliberately keeps the flag (`callbacks.jl:540`).

### 1.2 Which finite-slab case the callbacks docstring refers to

The docstring (`src/Solver/callbacks.jl:515-521`: "on the finite-slab heat-conduction case at 320 cells the L2 error grows 15× and the step count 2.6×") entered with commit a2dce2f, whose message generalizes it to "6-15x accuracy on stiff fine meshes". Both finite-slab families carried a 320-cell rung in their report ladders at that commit (`V-HC-02.jl:119`, `V-HC-03.jl:128`), and no test or script at a2dce2f records the measurement. Phase 0 reproduced the effect on both: V-HC-03 at 320 cells and the case's tolerance (reltol 1e-8) gives 51.6× / 2.04× (smoothed estimate) and 44× / 2.59× (raw estimate); V-HC-02c gives 88× / 2.01× at reltol 1e-6 and 1.00× at reltol 1e-8. The 2.6× step ratio matches V-HC-03 with the raw estimate, the 15× error ratio matches neither exactly, and the solver defaults changed after a2dce2f (KLU default, smoothed estimate). **Conclusion: the case cannot be identified uniquely; the docstring describes a range.** The Phase 1 gate therefore names both cases (V-HC-03 and V-HC-02c at 320 cells, both tolerances), as the design already does, and the new core's own documentation names the case id and tolerance of every number it quotes.

### 1.3 Pointer audit

Every `file:line` claim of the design that the plan relies on, checked against 17bc477.

| Design claim | Location checked | Verdict |
|---|---|---|
| Block-major state layout | `src/Problem/state_layout.jl:17-22` | holds |
| Residual reads the state from `u` | `residual.jl:54-58` (`populate!`) | holds |
| Prefix-sum mesh velocity "`residual.jl:401-460`" | `_compute_ale_mesh_velocity!` is at `residual.jl:329-373` (file is 537 lines) | content holds, range wrong |
| θ_A subtracted in the velocity | `residual.jl:359-361` | holds |
| θ from `volume_change_rate` | `src/Physics/volume_change.jl:106` | holds |
| Park/merge "`callbacks.jl:490-650`" | depletion callback `420-487`, step update `524-545` (file is 613 lines) | content holds, range wrong |
| Merge conserves enthalpy "`cell_park.jl:90`" | `cell_park.jl:88-92` | holds; see conflict C4 for the fixed-mesh branch |
| Structured walk "`structured_backend.jl:498-728`" | `_compute_jacobian!` is `300-590` (file is 655 lines) | content holds, range wrong |
| Dense geometry block | `operators/geometry.jl:226-262, 268-322` | holds |
| Dense n² sparsity lookup | `fast_sparsity_map.jl:27-55` | holds |
| Divergence reads one face area for both faces | `finite_volume.jl:44-88` (`A = state.A_face[f_left]`) | holds |
| Kinetics clamps ξ at zero; limiter multiplies the rate | `kinetics.jl:63`; structure `57-66` | holds |
| Porosity clamped to [0, 1] | `mixture_rules.jl:400` | holds |
| Surface Newton termination test; Dual-safe Newton | `boundary_temperature.jl:215`; body `130-330` (function at 62) | holds |
| `MeshTopology{D}` / `MeshGeometry{D}` scaffolds | `topology.jl:24` (struct), `76` (`mesh_topology`, not 83); `mesh_geometry.jl:32` (struct), `45` (constructor) | holds, two lines off |
| BC consumers hard-code `top` and `bottom` | 58 `Val(:top)`/`Val(:bottom)` sites in 7 files (Discretization, Jacobian operators, solve.jl, BoundaryConditions) | holds |
| Global scalars have no trial-state guard | `state_cache.jl` `update_cross_section_area!` throws on non-positive `A` | holds |
| Size: `src/` 31.8k lines, Jacobian 8.6k, Discretization 3.4k, Experimental 2.9k, Solver 2.7k, Problem 2.6k, Adaptivity 1.8k, Geometry 1.2k, Physics 1.1k, Residual 0.6k | `wc`: 31,823; 8,579; 3,387; 2,935; 2,749; 2,603; 1,772; 1,217; 1,131; 587 | holds |
| "17 value-dependent branches" in Discretization | not reproduced; three upwind switches (`ale_fluxes.jl:105`, `gas_transport.jl:107, 231`) plus value conditionals in the surface Newton, blowing and Darcy code, order of magnitude consistent | not verified exactly |
| Park threshold 5% of thickness, `DT_AFTER_PARK = 1e-6`, one cell per callback | `solve.jl:530`, `callbacks.jl:391`, docstring `394-418` | holds |
| KLU default, UMFPACK with the geometry block | `solve.jl:311-315` | holds |

Library facts the design cites, checked in the pinned versions of the repository `Manifest.toml`:

| Claim | Where | Verdict |
|---|---|---|
| `OrdinaryDiffEqCore` defaults initialization to `CheckInit` | `OrdinaryDiffEqCore` 4.18.0 `initialize_dae.jl:27-50` (`DefaultInit` → `CheckInit`) | holds |
| `KenCarp4` builds the `algebraic_vars` mask only for `ImplicitEuler`; the ESDIRK estimate sums every row | `OrdinaryDiffEqSDIRK` 2.9.6 `sdirk_caches.jl:188-191` | holds |
| The nonlinear solver's convergence test uses `internalnorm` | `OrdinaryDiffEqNonlinearSolve` 2.9.8 `newton.jl:551-555` (out-of-place step; the in-place step is `639-650`) | holds |
| SparseMatrixColorings and DifferentiationInterface are "installed" | present in `Manifest.toml` (0.4.28, 0.7.21) as indirect dependencies only; not in `Project.toml` `[deps]` | holds as stated; Phase 1 must add them as direct dependencies with `[compat]` entries (Aqua gates this) |
| `FBDF`, `Rodas5P` reachable | through the `OrdinaryDiffEq` umbrella (7.8.1); `KenCarp4` through the direct dependency `OrdinaryDiffEqSDIRK` | holds |

### 1.4 Conflicts between the design and the code, with proposed resolutions

Each entry records both locations. None is resolved silently; the owner's answer changes the scope of Phase 1 where marked.

**C1. The Phase 1 gate requires deactivation that Phase 1 does not build (scope decision).** Design §6 Phase 1 step 6 limits deactivation to "a deactivation callback stub (masking only)" and leaves "deactivation as the terminal event" to Phase 2 step 1, but the Phase 1 gate reads "all verification cases within their existing tolerances". V-CPL-02 and V-CPL-02b (Tier 2, Landau ablation) run with `handle_depletion = true, depletion_threshold = 0.05` (`V-CPL-02.jl:258-266, 408-413`) and park cells throughout; V-BO-01 to V-BO-05 are burnout cases by definition. *Proposed resolution:* bring the 1D terminal event into Phase 1 as step 6b, the design's conservative transfer by volume merge plus moving-boundary collapse, masked rows, the mass-fraction threshold and every crossing in one callback, which is "exactly the current merge" (§3.5) and has its Phase 0 prototype (`merge_state`) and contract (`test_transfer.jl`); leave the correction `c`, the threshold tuning and the gradient re-measurement in Phase 2. Cost: about two weeks (§7). The alternative keeps the design's split and removes the seven cases from the Phase 1 gate (they join Phase 2's); the plan is written for the first option and marks what moves if the owner chooses the second.

**C2. "Properties unchanged" against "no hard clamp remains in the residual" (§3.2).** The property layer the new core is told to import unchanged applies hard clamps that the trial-state contract forbids in the residual: `max(·, RHO_CP_FLOOR)` on the heat-capacity sum (`src/Properties/effective_properties.jl:314`), `clamp(1 − solid_vol, 0, 1)` on porosity (`mixture_rules.jl:400`), `max(k_eff, 0)` (`effective_properties.jl:487`), `clamp(ε, 0, 1)` and `clamp(a, 0, 1)` (`:614, :639`), `max(T, T_GAS_PROPERTY_FLOOR)` (`:241`), `max(μ, MU_FLOOR)` (`:66`), the face-mean cut band (`harmonic_mean_weighted`, zero below the cut). *Proposed resolution:* Properties gains smooth counterparts beside the hard ones as new functions (a smooth floor on the capacity sum, a smooth porosity floor, a smooth face-mean guard); the new core calls only those; the hard ones stay for the old core and are deleted at cutover. The property-law argument bracket `[T_PROPERTY_MIN, T_PROPERTY_MAX]` (`src/Materials/constants.jl:22-23`) is not a conflict: it is the accepted bracket the design's smooth trial-temperature clamp maps into, so the law's own clamp never binds.

**C3. "Physics unchanged" against the kinetics limiter.** The one-sided `C¹` limiter and finite base of §3.2 replace `max(ξ[idx], 0)` and the tanh `_depletion_factor` inside `Physics.reaction_rate` (`kinetics.jl:63-66`), which the old core shares; `Material.depletion_limiter` is the concrete `DepletionLimiter` (`materials.jl:202`, `reactions.jl:53-56`). *Proposed resolution:* an abstract limiter supertype in Materials, the current `DepletionLimiter` and a new one-sided limiter (`threshold`, join point `x₀`, finite-base scale `δ`) as its subtypes, `reaction_rate` dispatching on the limiter type, and the material field widened to the supertype. The old core's behaviour is unchanged; the new core validates at setup that the material carries the one-sided limiter. The `enabled = false` convention of the exact-solution cases (V-CPL-01/03/05, V-GT-02, V-RX-01/02/03) is kept: the finite base still applies, the vector field may then point out of the domain, and `isoutofdomain` is the guard. The alternative, a second Arrhenius kernel in the new core, is a dual code path and is rejected.

**C4. Fixed-mesh deactivation is not one of the three policies.** On a fixed mesh the current park keeps the receiver's volume and leaves the ghost's (`V_keep = ale ? V_p + V_k : V_k`, `cell_park.jl:90-98`): mass is transferred, volume is not, concentrations rise. §3.5's conservative transfer moves the receiver's far face; element death keeps the geometry but ledgers the residual instead of transferring it. V-BO-02 and V-BO-05 are fixed-mesh burnout cases with pinned QoIs. *Proposed resolution:* in the new core a fixed mesh is a mesh with every node pinned, and conservative transfer on a pinned mesh transfers inventories without volume (today's rule). Owner question for §3.5; if the owner prefers element death with a ledgered residual, V-BO-02's QoIs are re-pinned in Phase 2.

**C5. Constraint-residual rejection from a callback is not available (§3.8).** "The step-update callback evaluates the algebraic rows at the accepted state and rejects the step" cannot be done: a `DiscreteCallback` runs after acceptance (DiffEqBase `apply_discrete_callback!`), and OrdinaryDiffEq rejects a step only through its error estimate or the `isoutofdomain` predicate. *Proposed resolution:* the scaled constraint residual is part of the `isoutofdomain` predicate (a true rejection, the step is halved), the observer records the residual at every accepted step for the gate, and a residual that stays above tolerance across a bounded number of consecutive rejections terminates the solve with a return code. Cost: one residual evaluation per proposed step, which on `KenCarp4` is the evaluation the per-step flag already pays for.

**C6. The MMS driver injects sources by index arithmetic on the old layout.** `solve_mms` adds `du[i] += S_T` (K/s) and `du[n+(j-1)n+i] += S_ξ` (kg/(m³·s)) (`test/verification/exact/mms.jl:437-448`) and requires `jacobian = :fd` under `w_override` (`:416-424`). On the new core the sources enter through the layout (mass rows `S_ξ V_i`; the energy row `S_T` in the temperature form, `(Σ_j m_ij c_pj) S_T` in the inventory form, both pure functions of `(z, t)` so the Jacobian stays exact) and the `:fd` requirement disappears because the velocity override is a row the AD Jacobian differentiates. *Resolution:* the new assembler has a `sources` hook (per-cell volumetric sources per inventory and energy, and a nodal velocity override), and `solve_mms` gains a `core` keyword that routes through it.

**C7. V-CPL-06 reads old-core internals.** Its ledger reads `ws.rhs_cache.S_conv`, `ws.ale_advection`, `ws.θ_work`, `ws.prop_cache.ρcp` and `ws.top_face_ids` (`V-CPL-06.jl:88-215`). It cannot be run through an adapter. *Resolution:* it is re-authored as the "energy row checked term by term against the ledger with gas transport on" test of Phase 1 step 8, on the kernel-contribution interface (§3.4); the old case stays with the old core until cutover.

**C8. The burnout cases' step recorder reads temperatures by old-layout index.** `_vbo_solve_fn` reads `maximum(view(integ.u, 1:n))` (`V-BO-01.jl:150-156`). It is layout-specific and moves with the V-BO cases to Phase 2, where the recorder reads through the layout accessors.

**C9. Threshold semantics.** The new threshold is a mass fraction `ε_m` (§3.5, starting at 1e-2); V-CPL-02's tolerances were pinned under a 5% thickness threshold. There is no exact mapping. A QoI that leaves its tolerance under the new threshold is a gate finding, re-pinned only with the owner's agreement (the suite's 3× observed rule, `PYROLYSIS_VERIFICATION_PRINT_QOIS=1`).

**C10. One `Internal` namespace, two cores.** `src/Internal.jl` re-exports every submodule's exports; new modules exporting `solve`, `residual!` or `StateLayout` would collide. *Resolution:* the new modules are not added to `Internal` or to the top-level exports until cutover and are reached as `Pyrolysis.<Module>.name`; nothing in the old core imports them.

**C11. The test runner's directory list is fixed.** `test/runtests.jl:37-49` walks a literal list and isolates `jacobian/` files in modules. The new core's tests need an entry (`newcore`, isolated like `jacobian/` because they define types).

**C12. Default error estimate.** The old core's default is `KenCarp4()` with the smoothed estimate; the design's provisional default is the raw estimate (`smooth_est = false`). Not a conflict, but the Phase 1 gate compares the new core against the old core's actual default, so both numbers are reported.

**Design errata (line pointers).** The three ranges in the table above that exceed their files and the two off-by-few lines are corrected in the design document in the first commit (§9); no design content changes.

### 1.5 Harness status on the current tree

`julia --project=. scripts/mesh_phase0/runtests.jl --case=kinetics` passes on 17bc477 (374 assertions, exit 0), so the harness runs against the committed project environment. JuliaFormatter in check mode reports `scripts/mesh_phase0/` already formatted, so bringing it under version control adds no formatting churn and the Quality workflow's format gate (which covers `scripts/`) stays green.

---

## 2. Module layout

### 2.1 Rules

The two cores are separate modules until Phase 3. No module contains both an old and a new code path, no shim adapts one to the other, and the old core is not edited except for the four additive changes listed in §2.3 (a limiter type in Materials, smooth property variants in Properties, direct dependencies in `Project.toml`, and a pure move of the face-local boundary code from Discretization into BoundaryConditions). The new modules load after Physics and before `Problem`; they do not import `Geometry` (beyond the two scaffold types), `Discretization`, `Residual`, `Jacobian`, `Adaptivity`, `Solver` or `Experimental`, and nothing in those imports them. They are not added to `Internal` or to the top-level exports until cutover.

### 2.2 New modules

| Module (directory) | Owns | Imports unchanged from the existing layers | Phase 0 prototype it grows from |
|---|---|---|---|
| `Meshes` (`src/Meshes/`) | Face, cell and node types with tags; `MeshTopology{D}` and `MeshGeometry{D}` taken over from `src/Geometry/topology.jl`, `mesh_geometry.jl` as the connectivity and the per-step metric buffers; generators for 1D planar and 1D radial (cylinder, sphere) from node positions; face metrics (area, unit normal, tensor contraction) and cell volumes computed from node positions; the pinned-node and Lagrangian-node sets; a builder from a `Mesh1D` (node positions, cross-section area, boundary ids → tags) so every existing `PyrolysisProblem` is accepted | `BoundaryConditions` tag names; nothing else | `structured_mesh`, `cell_volume`, `sweep_operator` (`prototypes.jl`) |
| `States` (`src/States/`) | The state layout (block accessors, masks, tolerance vectors), the state mapping (the only reader of the raw state: derived cache of `T`, `ξ`, `V`, mass fractions, porosity, pressure, face metrics), the energy-state switch, the solid-carried fields (`χ` pair, reference orientation, deformation gradient), the extended caloric map for the inventory form, the Dual-typed derived caches for AD | `Materials` (components, `Material`, constants), `Properties` (mixture rules, property cache evaluation at a state) | `positive_mass`, `trial_temperature`, `capacity_floor`, `temperature_from_inventory`, `NonlinearCaloric.invert_temperature`, `material_axes` |
| `Kernels` (`src/Kernels/`) | The kernel interface (contributions: face fluxes, cell sources, strain rates, algebraic rows) and the 1D kernel set: face-scatter divergence, full-tensor conduction, Darcy–Fick gas transport in flux form with the enthalpy face flux, Beer–Lambert with lifted face intensities, reactions with the one-sided limiter, the strain rate `θ`, the per-cell volume balance with `θ_A`, the global cross-section rows, the surface balance per boundary face with blowing and absorptivity, tag-driven thermal, mass and pressure boundary conditions, the relative-flux reconstruction (van Albada) for `c ≠ 0` | `Physics` (`volume_change_rate`, `ideal_gas_pressure`, the pointwise `gas_flux`, `gas_flux_darcy_fick` face formulas, `reaction_rate` through the new limiter type, `evaluate_reactions!`, `species_source_terms!`, `reaction_heat_source`), `Properties` (face averaging, effective properties, Darcy permeability, radiative properties, gas viscosity, the smooth variants of C2), `BoundaryConditions` (every BC type, `BoundaryConditionSet`, `get_thermal_bc`, `get_mass_bc`, and, after the move of §2.3, `BlowingFace`, `blowing_state`, `blowing_ratio`, `mass_bc_flux`, `evaluate_bc_flux_and_derivative` and the per-face surface Newton) | the `rhs!` closures of `test_integration.jl` (shared inventory fluxes), `reaction_rate` (`prototypes.jl`) |
| `Closures` (`src/Closures/`) | The closure interface (unknowns, operators, boundary treatment, rank conditions) and its 1D implementation: the balance-row closure on nodal velocities, pinned and Lagrangian nodes, the consistent sweep, the `w_override` row; the potential closure's interface (cell multipliers) with a 1D implementation that must reproduce the balance rows, so the interface is exercised before Phase 4c | `Meshes` | `sweep_operator`, `closure_velocity`, the DAE closure of `test_integration.jl` |
| `Deactivation` (`src/Deactivation/`) | The three policies as types, the event (cells and nodes it changes, the inventory transfer, the ledger entry), the admissibility assertions, the mask; Phase 1 implements conservative transfer and moving-boundary collapse in 1D and the masking (C1); element death is an interface with a 1D pinned-mesh implementation (C4) | `Materials` (intrinsic densities, `c_p` integrals), `Properties` | `merge_state` |
| `Assembly` (`src/Assembly/`) | The residual: state mapping → kernels → row assembly (differential rows from the contributions, algebraic rows from the closure, intensity and global kernels, masked rows), the `sources` hook (C6), the ledger-term access the energy-row test needs; a pure function of `(du, u, def, t, ws)` like today's `residual!` | `Materials`, `BoundaryConditions` | — |
| `ColoredJacobian` (`src/ColoredJacobian/`) | The hand-built sparsity pattern from connectivity (memoized per problem structure), column coloring (SparseMatrixColorings), the colored ForwardDiff closure (DifferentiationInterface), the `jac!` and `jac_prototype`, the hand-written dense rows for the global scalars, the dense-AD oracle and the containment test helper | `ForwardDiff`, `SparseMatrixColorings`, `DifferentiationInterface`, `SparseArrays` | `structural_pattern`, `column_colors`, `ColorCache`, `colored_jac!` (`restart_probe.jl`); `baselines/lifted_containment.jl` |
| `Integration` (`src/Integration/`) | The integrator interface (candidates with their restart and per-step-flag policies), the `ODEFunction` with the diagonal mass matrix, the erased RHS and Jacobian types, `NoInit` with the project's consistent sweep, the full-norm error control with the tolerance vectors, the `isoutofdomain` predicate with the constraint check, the observer callbacks, the deactivation callback, the Newton-failure cap and progress callbacks (carried over), `solve` for the new core | `OrdinaryDiffEq`, `OrdinaryDiffEqSDIRK`, `LinearSolve`, `SciMLBase`, `FunctionWrappers`, `ProgressMeter` | `run_ladder`, `run_restarts`, `callback_control` (`solver_probe.jl`, `restart_probe.jl`) |
| `Outputs` (`src/Outputs/`) | The solution type with reconstructed fields (concentrations, temperatures, node positions, per-face surface temperatures and fluxes, area-weighted totals, mass-loss rate, thickness history, ledgers read from inventories, event log), the harness adapter surface (§6.1) | — | — |

Phase 1 deliverable per module is listed in §7. The `Experimental` namespace (P1, h-AMR) stays on the old core; its fate is the owner's question 5 and belongs to Phase 3.

### 2.3 What changes in the existing layers (additive only)

- `Materials`: abstract limiter supertype, the one-sided limiter type, `Material.depletion_limiter` field widened (C3). `Physics.reaction_rate` gains a method on the new type.
- `Properties`: smooth floors beside the hard ones (C2); an unclamped accessor for the per-component heat-capacity sum and the solid specific volume so the state mapping can apply its own smooth floors.
- `BoundaryConditions`: the mesh-free, face-local boundary code moves here from Discretization as a pure move with no behaviour change: `src/Discretization/blowing.jl` (`BlowingFace`, `blowing_state`, `blowing_ratio`, `blowing_factor`; it imports only Materials and Properties), `boundary_fluxes.jl` (`mass_bc_flux` and the thermal flux evaluators; no mesh reference), and from `boundary_temperature.jl` the surface Newton body (`solve_boundary_temperature`, lines 62-142, a pure function of its arguments) with the `evaluate_bc_flux_and_derivative` overloads (lines 311-660); the mesh-bound wrappers (`solve_boundary_temperature!`, `update_boundary_state!`, the mesh accessors) stay in Discretization and import the moved functions. Both cores then call one implementation. If the Newton's diagnostics sink turns out to be typed on `BoundaryState`, that plain struct moves with it.
- `Project.toml`: `SparseMatrixColorings` and `DifferentiationInterface` as direct dependencies with `[compat]`; `test/Project.toml` the same where the tests call them directly.
- `test/runtests.jl`: the `newcore` directory (C11).

### 2.4 The adapter

`PyrolysisProblem` stays the user-facing problem. `Integration.setup(problem; options)` replaces `initialize_state` for the new core: `apply_initial_conditions!(problem)` (unchanged), a `Meshes` mesh built from `problem.mesh` (node positions, cross-section area, boundary ids → `:top`, `:bottom`), the material and `bc_set` taken as they are, initial inventories `m_ij = ξ_j(z_i) V_i`, the energy state from `T(z_i)`, nodal velocities, face intensities and global scalars from the consistent sweep. The options are the current `solve` keywords the new core honours (`use_ale` → free surface node or all nodes pinned; `radiation_model`; `handle_depletion`, `depletion_threshold` mapped to the policy and `ε_m`; `abstol`, `reltol`, `maxiters`, `dtmin`, `dtmax`, `saveat`, `callback`, `progress`, `min_thickness`, `integrator`, `linear_solver`) plus the new ones (`energy_state`, `integrator_policy`, `closure`, `deactivation`). `jacobian` is accepted and ignored with a warning in Phase 1 (one Jacobian path).

### 2.5 Deleted at cutover (Phase 3, one commit)

`src/Jacobian/` entirely; `src/Discretization/` entirely; `src/Residual/`; `src/Adaptivity/`; `src/Solver/`; the `Mesh1D` files of `src/Geometry/` (`nodes.jl`, `cells.jl`, `faces.jl`, `mesh.jl`, `mesh_generation.jl`, `geometric_quantities.jl`) with `MeshTopology`/`MeshGeometry` already moved; `src/Problem/{state_layout,state_cache,workspace,adapter}.jl`; `src/Physics/convection.jl` and the suffix sweep in `src/Physics/radiation.jl`; the hard-clamped property variants and the tanh limiter; `src/Experimental/` ported or deleted per the owner's answer; the `jacobian/`, `discretization/`, `ale/`, `adaptivity/`, `residual/` test directories after their obligations are ported (§6.3).

---

## 3. Data structures and contracts

Each entry gives a type signature (Julia sketches, names provisional), the invariant it keeps, the Phase 0 prototype it generalizes and the contract test it must keep passing. Harness tests are referenced by file; they are ported into `test/newcore/` against the production types (§6.4) and the harness copies stay as the record of the Phase 0 gate.

### 3.1 Geometry

```julia
struct Face{D}
    id::Int
    nodes::NTuple{2^(D-1),Int}      # 1 node in 1D
    cells::Tuple{Int,Int}           # (owner, neighbour); 0 = boundary
    tag::Symbol                     # :interior, :top, :bottom, later :axis, :side, …
end
struct Cell{D}
    id::Int
    faces::NTuple{2D,Int}
    nodes::NTuple{2^D,Int}
end
struct MeshTopology{D}              # immutable: faces, cells, nodes, boundary_tags::Dict{Symbol,Vector{Int}},
                                    # pinned::Vector{Int}, lagrangian::Vector{Int}
struct MeshGeometry{D,T}            # per step, element type T so Duals pass: face_area, face_normal (SVector{D,T}),
                                    # face_center, cell_volume, cell_center, node_position (SVector{D,T})
geometry!(g::MeshGeometry, x::AbstractVector, topo::MeshTopology, cs::CrossSection)   # the only writer; reads x only
```

Invariants: `Σ_f A_f n_f = 0` over every cell's faces; `V_i = (1/D) Σ_f A_f (x_f · n_f)` equals the generator's volume formula (planar `A h`, cylinder `π(r₊² − r₋²)`, sphere `(4π/3)(r₊³ − r₋³)`) to round-off; `∂V_i/∂x_n` by AD of `geometry!` equals the face-sweep operator `D` of the closure (§3.5) to round-off on every generator; face normals point from `cells[1]` to `cells[2]`. The 1D generator reads node positions, so `Mesh1D`'s uniform and stretched meshes map one to one. Prototypes: `structured_mesh`, `cell_volume`, `sweep_operator`. Tests kept passing: `test_closure.jl` ("independent volume derivative matches face sweeps", the skewed and single-vertex meshes); new: the radial generators' areas and volumes against the formulas.

### 3.2 State layout and the state mapping

```julia
struct StateLayout{NC,E<:EnergyState,R,C}         # E ∈ {TemperatureState, EnthalpyState}; R radiation, C cross-section traits
    n_cells::Int; n_nodes::Int; n_faces::Int
    # per-node block:  x_n (diff), w_n (alg);  per-face block: I_f (alg, Beer–Lambert only)
    # per-cell block:  m_{i,1..NC} (diff), E_i (diff), G_i, M_i (diff, WithChi only)
    # trailing:        A, θ_A (alg, WithChi only)
    mask::BitVector                                # live rows; dead rows are identity rows
end
index(l, ::Val{:mass}, i, j); index(l, ::Val{:energy}, i); index(l, ::Val{:x}, n); index(l, ::Val{:w}, n); index(l, ::Val{:I}, f); …
mass_matrix(l) -> Diagonal            # 1 on differential rows, 0 on algebraic rows
abstol_vector(l, rtol, abstol, scales) -> Vector{Float64}
```

Ordering is interleaved by position (node block, cell block, node block, …, globals last) so that the Jacobian is banded in every geometry and "radial bandwidth equals planar" is a measurable statement; column coloring does not depend on the ordering. Kernels never compute an index: every access goes through `index`.

```julia
mutable struct DerivedState{NC,T}     # the derived cache; T = Float64 or a Dual
    T::Vector{T}; ξ::Matrix{T}; V::Vector{T}; Y::Matrix{T}; φ::Vector{T}; P::Vector{T}
    T_trial::Vector{T}                # smooth clamp of T into [T_min, T_max] (temperature form) — what every T-sensitive kernel reads
    cp_sum::Vector{T}                 # smoothly floored Σ_j m_ij c_pj
    χ::Vector{T}; A::Base.RefValue{T}; θ_A::Base.RefValue{T}
    geometry::MeshGeometry{1,T}
    props::PropertyCache{T}           # existing type, filled from (T, ξ)
    T_warm::Vector{Float64}           # last accepted temperatures: warm start of the inventory inversion
end
derive!(ds::DerivedState, u, layout, topo, material)     # the only reader of u
assemble_rows!(du, contributions, ds, layout)            # the only writer of du
```

Contract: `derive!` and `assemble_rows!` are the only two functions that know whether `E_i` is `T_i` or `H_i` and whether the inventories are masses; the Phase 1 gate's "suite passes once with `H` as the state" is the executable form of this contract, and a grep for `index(…, :energy` outside these two functions is part of the quality check. Inventory form: `T_i` is the root of `H_i = Σ_j m_j⁺ h_j(T) + κ(m)(T − T_ref)` with the smoothstep positive part and the capacity switch of §3.2, solved by the bracketed Newton with warm start and differentiated by the implicit-function corrections (the `invert_temperature` construction). Prototypes: `temperature_from_inventory`, `positive_mass`, `capacity_floor`, `NonlinearCaloric`. Tests kept passing: `test_caloric.jl`, `test_nonlinear_caloric.jl` (256-bit bisection oracle, implicit first and second derivatives, floor activation boundary), `test_caloric_comparison.jl` (both forms fourth order, inventory ledger below 1e-12), `test_energy.jl` (temperature form conserves energy at the integrator's order, inventory form to round-off).

### 3.3 Tolerance vectors

`abstol_k` per row: mass rows `rtol · ε_m · m_i^s(0)`, energy row in the inventory form `rtol · ε_m · m_i^s(0) · c_p,ref · ΔT_ref`, temperature rows the user's `abstol` in kelvin, node positions `abstol_x`, velocities `abstol_x / τ_x` with `τ_x` the first-step scale, intensities `rtol · I_inc`, global rows `rtol · A₀` and `rtol · |θ|_ref`. These are implementation parameters measured at the gate ("achieved error against requested tolerance"), not design choices; `m_i^s(0)` is the cell's total initial solid mass, never a per-species or reference-state inventory (§3.2).

### 3.4 Kernel interface

```julia
abstract type Kernel end
struct Contributions{NC,T}
    F_mass::Matrix{T}        # NC × n_faces: mass flux per unit face area along the face normal
    F_energy::Vector{T}      # n_faces: conduction + carried enthalpy + radiative (face) energy flux
    F_carried::Matrix{T}     # carried inventories (G, M, Q): flux per unit area
    S_mass::Matrix{T}        # NC × n_cells: ω_ij V_i
    S_energy::Vector{T}      # Q_rxn V_i + Q_rad V_i + Σ_j ω_ij h_j V_i (inventory form adds the formation term)
    S_carried::Matrix{T}
    θ::Vector{T}             # per-cell strain rate
    r_alg::Vector{T}         # algebraic rows written by the closure, intensity and global kernels
end
contribute!(c::Contributions, ::ConductionKernel, ds::DerivedState, topo, material, bcs, t)
```

Contract: a kernel reads the derived state, the topology and metrics, the material and the boundary conditions; it never reads `u`, `du` or another kernel's output except through `c` (the enthalpy flux reads `F_mass` of the same face, so gas transport runs before the energy kernel, and the order is fixed by the assembler, not by the kernels). Face fluxes are scattered as `±F_f A_f` by the assembler into both cells; boundary faces receive their flux from the boundary-condition kernel on the same arrays. The design's trial-state contract (§3.8 of this plan) applies to every kernel: finite value and finite derivative at any real state. Prototypes: the `rhs!` closures of `test_integration.jl` ("isothermal composition transport with common inventory fluxes", "weighted moving-mesh closure through AD and a DAE restart"). Tests kept passing: those two, plus `test_fields.jl` ("reaction on, transport off", "reaction off, transport on").

### 3.5 Closure interface

```julia
abstract type Closure end
struct BalanceRowClosure <: Closure end          # 1D and column mode: unknowns w_n
struct PotentialClosure  <: Closure end          # unknowns φ_i; 1D implementation must reproduce BalanceRowClosure
unknowns(::Closure, topo) -> block description
residual_rows!(r_alg, ::Closure, u_alg, ds, topo, s)       # per cell: Σ_f A_f (w·n)_f − (θ_i − θ_A) V_i − s_i
velocity!(w, ::Closure, u_alg, topo)                        # identity for balance rows; w = W⁻¹ Dᵀ φ for the potential form
consistent_init!(u_alg, ::Closure, ds, topo, s)             # one bidiagonal (1D) or SPD solve
rank_check(::Closure, topo) -> Bool                         # D has full row rank for this boundary set
override!(r_alg, w_fun, x, t)                               # w_override: rows become w_n − w_fun(x_n, t)
```

Invariants: at every accepted step every live cell's balance residual, scaled by `|θ_i| V_i + Σ_f A_f |w_f|`, is below `rtol` (the constraint check of §5.4); a uniform `θ` on a cylinder gives `w = θ r/2` and on a sphere `θ r/3`; the pinned back node has `w = 0`; the surface node moves at the physical surface velocity; with the lumped-volume weights two unit squares under uniform expansion move uniformly `(1, 1, 1)` (the unweighted form gives `(2/3, 4/3, 2/3)`). The cell-multiplier form is in the interface so that Phase 4c changes the implementation, not the callers. Prototypes: `sweep_operator`, `closure_velocity`. Tests kept passing: `test_closure.jl` (two-square negative control, connected 2D/3D meshes, rank under each boundary set, fully pinned body throws), `test_integration.jl` (closure DAE through AD and a restart, constraint residual below 1e-8 on every accepted state).

### 3.6 Deactivation policies

```julia
abstract type DeactivationPolicy end
struct ConservativeTransfer <: DeactivationPolicy end    # receiver absorbs volume (free mesh) or inventories only (pinned mesh, C4)
struct MovingBoundaryCollapse <: DeactivationPolicy end  # dead nodes coincide with the receiver's new far face
struct ElementDeath <: DeactivationPolicy end            # geometry kept, residual ledgered (interface; 1D pinned-mesh instance only)
struct DeactivationEvent; cells::Vector{Int}; nodes::Vector{Int}; ledger::LedgerEntry; end
detect(policy, ds, layout, ε_m) -> Vector{Int}           # every cell below ε_m of its initial solid mass, same step
apply!(u, mask, topo, policy, cells, material) -> DeactivationEvent
admissible(ds) -> Bool                                   # porosity in [0, 1], positive volumes, finite pressure
```

Invariants (§3.5): a merge never destroys pore space (merged volume = max of the summed geometric volumes and the merged skeletal volume plus the summed pre-merge pore volumes); the created volume is ledgered; total inventories are conserved to round-off; the merged temperature follows the merged enthalpy; admissibility holds after every event; the state length and the sparsity pattern never change; every crossing at a step is handled in one callback, re-testing the predicate before each further transfer (Phase 0 item 7). Prototype: `merge_state`. Tests kept passing: `test_transfer.jl` (pore-free counterexample, gas-containing merges with finite pressure, a pore-free donor with a porous receiver); the one-callback parking result of Phase 0 (`results/one_callback_parking.md`) is reproduced on the new core in Phase 2 when V-BO-04 runs.

### 3.7 Solid-carried fields

```julia
struct CarriedField{S}
    name::Symbol                     # :G, :M, :F, :orientation, …
    frame::Symbol                    # :reference (constant per cell) or :material (convected)
    packed::Bool                     # intensive a carried as Q = m_s a
    source::S                        # per-cell source of the packed inventory: a (dm_s/dt)_rxn, and for F: m_s (∇v_s) F
end
material_axes(F, reference) -> orthonormal frame        # priority: grain, then radial projected, then cross product
physical_tensor(k_principal, axes) -> SMatrix{D,D}
```

Invariants: a constant intensive field stays constant while its carrier reacts away (the packed row carries `a (dm_s/dt)_rxn`); with `c = 0` nothing moves relative to the cells; `χ_i = G_i/M_i` and the column mean `χ̄ = Σ G_i / Σ M_i` (the current mass-weighted mean, `state_cache.jl:_chi_bar_from_state`); a rigid 90° rotation of `diag(4, 1)` gives `diag(1, 4)`; a 45° fibre under `diag(2, 1)` ends at `atan(1/2)` with unchanged principal values; the declared axis priority decides the transverse direction. In Phase 1's 1D core `F` is the axial stretch and the orientation is the axis; the full-frame functions are implemented and tested now so no kernel changes in Phase 4. Prototype: `material_axes`. Tests kept passing: `test_fields.jl` (all three testsets).

### 3.8 Trial-state contract

Every residual evaluation at any real state returns finite values and finite derivatives; the physical vector field never points out of the non-negative domain; the admissibility predicate reads unguarded quantities. The guarded quantities, listed:

| Quantity | Guard (smooth, `C¹`) | Applied where | Predicate reads |
|---|---|---|---|
| Reactant concentration in a rate | one-sided limiter `g(ξ/ξ_th)`: 0 at and below 0, joins `tanh` at `x₀` (default 1); finite base `δ·softplus(ξ/δ)` for the power | `reaction_rate` (new limiter type) | unguarded `m_ij < 0` |
| Temperature in Arrhenius, gas density, radiation, property laws | smooth clamp into `[T_PROPERTY_MIN, T_PROPERTY_MAX]` with a 1 K transition | `ds.T_trial`, read by every temperature-sensitive kernel; properties then never hit their own hard clamp | unguarded `T` |
| Heat-capacity sum (temperature form) | smooth floor at `RHO_CP_FLOOR`-equivalent per cell, inactive above twice the floor | `ds.cp_sum` | — |
| Inventory inversion (enthalpy form) | positive part `m⁺` (zero below 0, exact above `δ`, smoothstep between) and capacity switch `κ` | `derive!` | unguarded masses |
| Porosity | smooth floor replacing `clamp(·, 0, 1)` (C2) | `ds.φ`, pressure, Darcy | unguarded `1 − Σ m_j/(ρ_j V)` outside `[0, 1]` |
| Cell volume | smooth floor for `m/V`, face metrics and the activation switch | `geometry!` | unguarded `V ≤ 0` |
| Column mean `χ̄` | smooth clamp into `[0, 1]` before the lateral law | global rows | — |
| Face conductivity / transfer means | smooth guard band replacing the cut (`harmonic_mean_weighted`) | Properties smooth variant | — |
| Surface Newton iterate | starts from the clamped cell temperature, iterate clamped (today's rule) | per-face surface balance | — |

Prototypes: `depletion_factor`, `reaction_rate`, `trial_temperature`, `positive_mass`, `capacity_floor`, `softplus`. Tests kept passing: `test_kinetics.jl` (finite rates and derivatives at absent reactants for orders 0 to 2, `C¹` onset, rates above the join equal the production kernel to 1e-14, trial temperature safe in representative kernels, absent reactant through a reacting step, per-order approach to depletion), `test_caloric.jl`. No hard clamp remains in the new core's residual; the quality check greps the new modules for `clamp(`, `max(`, `min(` on state-derived quantities and the review signs off on each remaining use (array bounds and index arithmetic only).

---

## 4. Kernel-by-kernel plan

In the design's order. Each entry: the equations (§3 of the design), the current implementation to read first, what changes, and the exact-solution case the kernel lands with. Every kernel lands with its dense-AD check (the containment and value oracle of §4.9 run at 6 to 8 cells for the kernel's physics options) before the next starts.

### 4.1 Geometry

*Equations.* Faces carry `A_f`, `n_f`; cells carry `V_i(x)`; planar `A h`, cylinder and sphere from the radius; `dV_i/dt = Σ_f A_f (w·n)_f` semi-discretely (§3.3, discrete GCL). *Read first.* `src/Geometry/topology.jl`, `mesh_geometry.jl` (the scaffolds), `mesh.jl` (`Mesh1D`, `BoundaryState`), `mesh_generation.jl`, `geometric_quantities.jl` (`update_geometry!`, the opposite-distance gradient), `src/Problem/state_cache.jl:update_state_geometry!` (the Dual-typed geometry recompute the residual uses today). *Changes.* Face-based types with tags; metrics computed from `x` into a `MeshGeometry{1,T}` whose element type follows `u` so AD passes through; the radial generators; the `Mesh1D` builder; the shared-area divergence disappears with the per-face area (`finite_volume.jl:44-88`). *Lands with.* No verification case of its own; the unit tests of §3.1, the radial area and volume formulas, and the free-stream convergence test of §6.5 (which needs the velocity rows of §4.8, so it closes there).

### 4.2 State

*Equations.* §3.2: `m_ij`, `E_i`, `x_n`, `w_n`, `I_f`, `G_i`, `M_i`, `A`, `θ_A`; derived `V_i`, `ξ_ij = m_ij/V_i`, `T_i`, mass fractions, porosity, ideal-gas pressure; the temperature-form row `(Σ_j m_ij c_pj) dT_i/dt = RHS_i − Σ_j h_j(T_i) dm_ij/dt`; the inventory row `dH_i/dt = RHS_i`; sensible enthalpies `h_j(T) = ∫_{T_ref}^T c_pj` from one reference. *Read first.* `src/Problem/state_layout.jl`, `state_cache.jl` (`populate!`, `update_cross_section_area!`, `_chi_bar_from_state`), `workspace.jl`, `adapter.jl` (`to_problem_def`, `build_workspace`, `build_ad_workspace`), `src/Adaptivity/cell_merge.jl` (the `∫c_p dT` quadrature the merge already uses: reuse it for `h_j`), `src/Properties/property_cache.jl`, `effective_properties.jl:300-320` (`effective_heat_capacity`). *Changes.* The layout of §3.2 with masks and tolerance vectors; `derive!` as the only reader; both energy forms behind it; the Dual-typed derived caches per chunk width (today's `build_ad_workspace` pattern); initial inventories from the problem's `ξ` callables times the initial volumes. *Lands with.* V-RX-01, V-RX-02a-d, V-RX-03 (one-cell kinetics, no transport: state plus reactions), run in both energy forms.

### 4.3 Conduction

*Equations.* `q_f = −n_f · K_f ∇T|_f`, two-point normal part plus the cross-diffusion part evaluated always (§3.6 Transport); in 1D `q_f = −k_f (T_R − T_L)/d_LR` with the distance-weighted harmonic face conductivity; the energy row takes `−Σ_f q_f A_f`. *Read first.* `src/Physics/heat_conduction.jl` (38 lines), `src/Discretization/finite_volume.jl:225-243` (step 1), `src/Properties/effective_properties.jl` (face averaging, `harmonic_mean_weighted`), `src/Geometry/geometric_quantities.jl` (gradient stencils), Technical Reference §7.1-7.3, §13.6-13.7. *Changes.* A face loop over the topology with per-face area; the tensor signature from the first line (an `SMatrix{1,1}` in 1D); the smooth face-mean guard (C2); boundary faces take the surface-balance flux (§4.7). *Lands with.* V-HC-01 (semi-infinite, convective), V-HC-03 (finite slab, flux), V-HC-02a-d (Biot sweep), V-HC-04a/b (two-layer), V-HC-05 (`k(T)`, Kirchhoff). The minimal tag-driven boundary set these need (heat flux, adiabatic, convective with the per-face Newton, Dirichlet) lands here; the radiative, blowing and absorptivity members land in §4.7.

### 4.4 Reactions

*Equations.* `dm_ij/dt ∋ ω_ij V_i`, `ω` from the Arrhenius rate with the one-sided limiter and finite base (§3.2), the smooth temperature gates as now; `Q_rxn(T_i) V_i` with heats quoted at the local temperature; in the inventory form the additional `Σ_j ω_ij h_j(T_i) V_i` that makes that convention exact; reaction rates multiplied by the smooth remaining-mass-fraction switch-off near `ε_m` (§3.5; implemented now, inactive without deactivation). *Read first.* `src/Physics/kinetics.jl` (`reaction_rate`, the generated `evaluate_reactions!`, `species_source_terms!`, `reaction_heat_source`), `src/Materials/reactions.jl` (`Reaction`, `DepletionLimiter`, stoichiometry), `finite_volume.jl:244-264` (step 2), Technical Reference §6. *Changes.* The limiter type (C3); evaluation at `T_trial`; sources into `S_mass`, `S_energy`, `S_carried` (the `G` row sums the dry-gas sources, `Materials.dry_pyrolysis_gas_indices`). *Lands with.* V-RX-01..04 (re-run with the state), V-CPL-04 (Frank-Kamenetskii), V-CPL-05 (TGA ladder); the depletion contract tests at orders 0, 0.5, 1; V-RX-02d (limiter on through burnout) is the calibration guard: rates above `x₀ ξ_th` equal the production kernel to 1e-14.

### 4.5 Gas transport and the enthalpy face flux

*Equations.* `F_{ij,f}` for every gas `j`: the volume-fraction-gradient diffusive flux `−(λ_f/T_f) (ξ_R T_R − ξ_L T_L)/d` and, under Darcy–Fick, the advective `(ξ/φ)_upwind v_g` with `v_g = −(κ/μ) ∇P` and `P = Σ_j ξ_j R T/(M_j φ)`; mass rows `dm_ij/dt ∋ −Σ_f F_{ij,f} A_f`; the energy row carries `Σ_f Σ_j F_{ij,f} h_j(T_f) A_f` with `T_f` the distance-weighted face temperature (centred), which replaces `S_conv`; the surface keeps no gas-enthalpy term (the gas leaves at `h_g(T_s)` inside the flux sum); mass and pressure boundary fluxes at boundary faces including Darcy venting. *Read first.* `src/Physics/gas_transport.jl` (`gas_flux`, `gas_flux_darcy_fick`, `compute_gas_fluxes!`, boundary venting at `:231`), `src/Physics/pressure.jl`, `src/Discretization/darcy_velocity.jl`, `finite_volume.jl:266-290` (step 3), `:322-337` (step 6, what is replaced), `:629-750` (`apply_mass_bc_flux!`), `:777-870` (`update_darcy_properties!`), `src/Physics/convection.jl` (the `S_conv` being replaced and its boundary blend), `src/Discretization/boundary_fluxes.jl`, Technical Reference §7.4, §9, the V-CPL-06 ledger construction. *Changes.* Flux form into `F_mass`; the enthalpy flux into `F_energy` from the same `F_mass`; porosity and pressure from the derived state with smooth floors; per-face BC fluxes by tag. The relative solid flux for `c ≠ 0` (van Albada reconstruction, upwind face state, enthalpy with the same flux) is implemented here and exercised in Phase 1 only under `w_override`. *Lands with.* V-GT-01 (Gaussian gas pulse, Fickian), V-GT-02 (steady compressible Darcy outflow), V-MMS-01 (coupled, gas transport on), the energy-row-against-ledger test (§6.5, replacing V-CPL-06a/b), the boundary-venting and Darcy obligations of `test/discretization/test_darcy*.jl` and `test/jacobian/test_darcy_venting_jacobian.jl` as combined-physics containment states.

### 4.6 Radiation

*Equations.* Lifted face intensities: per cell `0 = I_{f⁻} − I_{f⁺} exp(−α_i Δ_i)`, the exposed-face row `0 = I_top − a_eff I_inc(t)` (Kirchhoff `a_eff = ε_eff` of the surface cell by default, the explicit absorptivity otherwise), the cell source `Q_rad,i V_i = (I_{f⁺} − I_{f⁻}) A` (the exact cell-averaged absorption, identical to today's `I_in (1 − e^{−τ})/Δz`), and the transmitted flux is `I` at the back face, a state the outputs read directly; surface absorption stays a surface-balance term. *Read first.* `src/Physics/radiation.jl` (`compute_beer_lambert_source!`), `src/Properties/radiation_cache.jl`, `effective_properties.jl` (`effective_emissivity`, `effective_absorptivity`), `src/Jacobian/operators/radiation.jl` and `sparsity.jl:_emit_radiation_pattern!` (the suffix coupling the lifting removes), `src/Solver/solve.jl:971` (`_get_radiation_flux_function`), Technical Reference §8.3, §8.10. *Changes.* Intensities as algebraic unknowns with local rows; the consistent sweep fills them by the recurrence at `t₀` and after events; the energy row reads two faces. *Lands with.* V-RAD-01 (in-depth absorption with conduction, sub-surface peak); the Kirchhoff and `ε(T)` chains of `test_kirchhoff_eps_eff_jacobian.jl` and `test_surface_radiation_eps_T_jacobian.jl` as containment states; V-BO-03's transmitted-flux guard in Phase 2.

### 4.7 Boundary conditions

*Equations.* Per boundary face with tag: the surface balance `q_bc(T_s) = k (T_s − T_cell)/Δz_half` solved by Newton per face (Dual-safe, started from the clamped cell temperature), `q_bc` from the thermal BC members (flux, convective with the blowing ratio of the face's outward mass flux, re-radiation `ε_eff σ T_s⁴`, incident `a_eff q_inc` under surface absorption), blowing conductance and ratio as per-face diagnostics; mass BCs (impermeable, prescribed, convective, diffusive with the λρ formulation) and pressure BCs per face. *Read first.* `src/Discretization/boundary_temperature.jl:62-330` (the Newton), `boundary_fluxes.jl`, `blowing.jl` (`BlowingFace`, `blowing_state`, `blowing_ratio`), `finite_volume.jl:292-307` and `:508-600` (`_apply_thermal_bc_static!`), `src/BoundaryConditions/boundary_conditions.jl` (types, `BoundaryConditionSet`, `get_thermal_bc`, `has_radiative_flux_bc`), `src/Jacobian/operators/bc_thermal.jl`, `bc_mass.jl`, `bc_convection.jl` (which couplings the `∂T_s` chain creates), Technical Reference §12. *Changes.* Consumers iterate `boundary_tags`; one `FaceBalanceState` per boundary face (replacing `mesh.boundary_states[id]`); no `Val(:top)`; the surface Newton's math is kept verbatim (it is already a pure function of its arguments, and after the move of §2.3 it lives in BoundaryConditions) and re-hosted on the face. *Lands with.* V-HC-02 (convective Newton), V-RX-04a/b (radiative cooling, radiative-convective equilibrium), the ports of `test/boundary_conditions/` (blowing models, blowing correction, boundary temperature) and `test/solver/test_surface_bc_guards.jl`, the boundary-flux energy-row obligation of `test_jacobian_coupling_regressions.jl` as a combined-physics containment state (diffusive mass BC at both faces, Fickian and Darcy, nonlinear surface balance).

### 4.8 Mesh motion

*Equations.* `dx_n/dt = w_n`; the balance rows `A_{k+1} w_{k+1} − A_k w_k = (θ_k − θ_A) V_k + s_k` from the pinned back node (`s = 0` in Phase 1), the surface node moving with the physical surface; `θ_i = Σ_j (γ_j/ρ_j(T)) ω_ij` (`volume_change_rate` on the reaction mass rates); `G`, `M` rows; the global rows `0 = A − A₀ law(χ̄)`, `0 = θ_A − law'(χ̄) χ̄̇ / law(χ̄)` with `χ̄ = Σ G/Σ M`, `χ̄̇ = Σ Ġ/Σ M`; no dilation term (masses are the states); the `w_override` row `0 = w_n − w(x_n, t)`; relative fluxes active only when `w ≠ v_s`. *Read first.* `src/Residual/residual.jl:329-373` (velocity), `:374-537` (`θ_A`, χ rows, dilation, cross-section), `src/Physics/volume_change.jl`, `src/Discretization/ale_fluxes.jl` (the mesh-relative temperature advection the inventory fluxes replace), `state_cache.jl` (`update_cross_section_area!`), `src/Jacobian/operators/ale_mesh_velocity.jl`, `chi.jl`, `sparsity.jl:_emit_geometry_pattern!` (the prefix-sum and rank-1 couplings that become local), Technical Reference §10, §11, `test/verification/exact/mms.jl:240-300` (`mms_mesh_velocity`, `MMSMotion`). *Changes.* Velocities as algebraic unknowns through the closure interface; the strain rate from the mass sources; the cross-section scalars as two dense rows; the GCL holds semi-discretely by construction (volumes are read from `x`, sweep rates from the same `w`). *Lands with.* V-CPL-03 (shrink/swell endpoints, six permutations and the WithChi law), V-MMS-02 (free-stream preservation under prescribed motion), V-MMS-03 (prescribed smooth motion), V-MMS-04 (induced motion, material-frame sources), the isothermal advection test and the radial free-stream convergence test (§6.5), the history-transport and rotation tests, and, with step 6b of C1, V-CPL-02/02b (Landau ablation). V-CPL-01/01b (Stefan, Eulerian) close under §4.3-4.4.

### 4.9 Jacobian

*Equations.* Pattern from connectivity: a cell row reaches the unknowns of cells `i−1..i+1` and nodes `i−1..i+2` with `c = 0` (the containment finding: conduction reads neighbour centres and volumes), cells `i−2..i+2` and nodes `i−2..i+3` with the five-cell reconstruction active, its two face intensities, the `A` column; a velocity row reaches `w_k, w_{k+1}, x_k, x_{k+1}`, its cell's unknowns (through `θ_k`) and `θ_A`; an intensity row reaches `I_f, I_{f+1}` and its cell; the boundary cell's rows reach the neighbour cell through the face flux the surface Newton reads; the two global rows are dense and written by hand after the colored pass; every row stores its diagonal; masked rows are diagonal only. Expected colors `3NC + 9 + 2·bl + 7·χ` (15 to 24 for two components, 20 for V-BO-01) and `5NC + 13 + 2·bl + 11·χ` with the reconstruction active (Phase 0 item 4). *Read first.* `src/Jacobian/sparsity.jl`, `fast_sparsity_map.jl`, `dense_ad_reference.jl`, `verify.jl` (`_strict_per_entry_compare`: the per-entry criterion with the row-magnitude floor is kept), `operators/geometry.jl` (the Dual workspace pattern), `scripts/mesh_phase0/baselines/{jacobian_colors,lifted_containment,ad_overhead}.jl`, `restart_probe.jl` (`ColorCache`, `colored_jac!`). *Changes.* `pattern(layout, topo, options)` memoized; `SparseMatrixColorings.coloring` (column partition, greedy, largest-first fallback); `DifferentiationInterface.prepare_jacobian` with `AutoSparse(AutoForwardDiff(); sparsity_detector = KnownJacobianSparsityDetector(pattern), coloring_algorithm = …)` and `jacobian!` writing into the integrator's `SparseMatrixCSC` in place (pattern identity with `jac_prototype` asserted by `colptr`/`rowval` equality at build time); the dense rows added after; the containment helper: dense ForwardDiff of the production residual at 6 to 8 cells for every kernel family and physics option, `nnz(dense) ⊆ pattern` and values equal to round-off. No tracing. *Lands with.* The containment suite over the configurations of `test/perf/perf_helpers.jl` (baseline Eulerian, surface absorption, Darcy, ALE no radiation, ALE Beer–Lambert Darcy, Beer–Lambert moderate) plus WithChi and the boundary-flux states; V-MMS-01 at its tight tolerance as the first end-to-end case on the colored Jacobian.

### 4.10 Solver

§5.

### 4.11 Outputs

*Read first.* `src/Problem/pyrolysis_problem.jl:226-282` (`PyrolysisSolution`, the `extras` schema), `src/Solver/solve.jl:989-1300` (`extract_solution!`), `:1312-1600` (`extract_solution_ale!`, NaN-padding above the surface, `park_times`, `live_cells`), `src/Diagnostics/conservation_checks.jl` (ledgers), `ext/PyrolysisPlotsExt.jl` (what the recipes read), `test/verification/casedefs.jl`, the solution fields the cases read (§6.1). *Changes.* A new solution type with reconstructed fields; ledgers read inventories; per-face surface temperatures and fluxes with area-weighted totals; the mass-loss rate from the exposed face's mass flux; the event log. *Lands with.* The harness adapter and every case above.

---

## 5. Solver plan

### 5.1 One interface, three candidates

```julia
struct IntegratorPolicy{Alg}
    alg::Alg                       # KenCarp4(smooth_est = false), Rodas5P(), FBDF()
    restart::Symbol                # :solver_chosen or :forced_dt (value: DT_AFTER_EVENT)
    observer_flag::Symbol          # :keep (KenCarp4) or :clear (FBDF, Rodas5P indifferent)
    needs_exact_jacobian::Bool     # Rodas5P
end
```

Defaults and their Phase 0 evidence (§3.8 of the design; `results/integrator_flag_summary.md`, `restart_comparison.csv`): `KenCarp4` with the raw estimate, the per-step flag kept and the solver-chosen restart is the default (qualified at 100 events and production tolerances: error 8.1e-7, 0.207 s on the lifted ladder); `Rodas5P` is the favourite (qualified under both restart policies, 3.4e-8, 0.183 s; bit-identical under every flag variant; needs the exact Jacobian, which the colored AD provides; its cost scales with the production color count, 2.2 to 4.8× the ladder's per-step Jacobian); `FBDF` needs the observer flag cleared (with the flag it stalls at 3 of 1,200 s) and a forced small restart step at events (2.2e-5 against the 1e-5 budget with the solver-chosen step), and is the fastest at tight tolerance on the finite slab (0.0284 s against 0.06 s for `KenCarp4`, same L2 at reltol 1e-8). Every observer callback asks the policy which flag value to set; the deactivation callback always sets the flag, runs the consistent sweep after the transfer and, under `:forced_dt`, proposes the restart step. `Rodas5P` needs `∂f/∂t`: the erased RHS exposes it through finite differences (`autodiff = AutoFiniteDiff()` as today) because the incident flux is a function of `t`.

### 5.2 Problem construction

`ODEFunction{true,FullSpecialize}(ErasedRHS(rhs!); mass_matrix = mass_matrix(layout), jac = ErasedJacobian{SparseMatrixCSC{Float64,Int}}(jac!), jac_prototype)`; `KLUFactorization()` default, `UMFPACKFactorization()` selectable (the pivot-growth failure on the Eulerian boundary-layer cases recurs: `solve.jl:300-310`); `initializealg = NoInit()`; `abstol` the vector of §3.3; `reltol` as now; `isoutofdomain` of §5.4; `maxiters`, `dtmin`, `dtmax`, `saveat` as now; the erased types, wrapped callbacks, progress callback and the Newton-failure cap carried over from `src/Solver/integrator_functions.jl`, `callbacks.jl:584-613`, `solve.jl:775-789`. The precompile workloads for the new core are added at cutover, not in Phase 1.

### 5.3 Consistent initialization

`consistent_sweep!(u)`: face intensities by the recurrence from the exposed face; nodal velocities by the closure's `consistent_init!` (one bidiagonal solve in 1D, after one kernel pass for `θ`); the global scalars from `χ̄`; then the algebraic residual is evaluated and must be below the constraint tolerance before the integrator starts. The sweep runs at `t₀` and inside every callback that modifies `u` (deactivation), where the stepper re-enters initialization with `NoInit` and would otherwise continue from an inconsistent state (`reeval_internals_due_to_modification!` → `initialize_dae!`, §1.1).

### 5.4 Three separate checks

1. **Step error.** The integrator's own estimate on every row, algebraic rows included, in the default `internalnorm`; the raw estimate for `KenCarp4`.
2. **Nonlinear convergence.** The integrator's Newton, with the same norm (`newton.jl:551-555`, `639-650`); never masked.
3. **Constraint accuracy.** `isoutofdomain(u, p, t)` returns true when any live balance residual exceeds `rtol × (|θ_i| V_i + Σ_f A_f |w_f|)`, any intensity residual exceeds `rtol × I_inc`, or a cross-section residual exceeds `rtol × A₀` or `rtol × |θ|_ref` (C5); the same predicate returns true on negative mass, physical porosity outside `[0, 1]` or non-positive volume, read unguarded. A step rejected this way is halved by the controller; more than a configured number of consecutive constraint rejections ends the solve with `ReturnCode.Failure` and a diagnostic. The step-update observer records the scaled residuals at every accepted step; their maximum is a gate quantity.

### 5.5 Running the finite-slab gate

The protocol of `scripts/mesh_phase0/baselines/integrator_flag.jl` is re-hosted as `test/newcore/gates/finite_slab.jl` (a script, not part of `Pkg.test()`): V-HC-03 and V-HC-02c at 320 cells, the two tolerance pairs (reltol 1e-6/abstol 1e-8 and 1e-8/1e-10), the three candidates each under its own policy, L2 and L∞ against the series at every output time, steps, RHS calls, Jacobians, warm time (median of three). Pass criterion: for the chosen candidate and policy, L2 max ≤ the old core's `KenCarp4` default at the same tolerance (V-HC-03: 1.07e-4 K at reltol 1e-8, 5.89e-3 K at 1e-6; V-HC-02c: 5.93e-4 K and 2.20e-4 K), with the old numbers re-measured on the same machine in the same session. The old core is driven by the existing script unchanged.

---

## 6. Test and verification plan

### 6.1 The adapter

The verification suite reads a fixed set of solution fields: `sol.t`, `sol.T`, `sol.z` (node positions, NaN above the surface after events), `sol.surface_T`, `sol.mass_loss_rate`, `sol.retcode`, `sol.problem`, and from `extras`: `thickness_history`, `total_mass`, `gas_generation_rate`, `surface_mass_flux`, `cross_section_area_history`, `park_times`, `live_cells`, `radiation_transmitted`; the report harness also reads `back_temperature` and the surface-flux decomposition. The new `Outputs` solution type provides these names with the same semantics (the type is new; no shim converts one solution into the other).

Routing: `test/verification/casedefs.jl` gets a one-line switch in `_case_solve` on the environment variable `PYROLYSIS_CORE` (`old` default, `new`); under `new`, a case without a `solve_fn` is solved by `Integration.solve(problem; progress = false, kwargs...)` with the keyword mapping of §2.4, the MMS cases go through `solve_mms(…; core = :new)` (C6), the V-CPL-06 pair is replaced by the energy-row test of §6.5 under `new`, and the V-BO cases are skipped under `new` until Phase 2 (C8) unless C1's step 6b lands, in which case V-BO-01/03/04 run and V-BO-02/05 wait for the pinned-mesh policy decision (C4). `Pkg.test()` runs the old core; a second invocation with `PYROLYSIS_CORE=new` runs the new one, and the Phase 1 gate requires both green. The convergence spot checks, the wall-time budgets (`wall < 5 × budget_s`) and the pinned tolerances apply unchanged under `new`.

The report mode (`scripts/verification/run_reports.jl`, Gpyro and FDS overlays from `test/verification/reference/`) consumes the same case definitions and the same solution fields, so it runs under `new` without change; the overlays are data and do not move.

### 6.2 What runs through the adapter unchanged

Every Tier 1 case (V-HC-01..05, V-RX-01..04, V-RAD-01, V-GT-01..02), every Tier 2 case without events (V-CPL-01/01b, V-CPL-03, V-CPL-04, V-CPL-05), the MMS driver (V-MMS-01..04) through its `core` keyword, and the exact-solution self-tests (`test/verification/exact/selftests.jl`, which are solver-free). With step 6b: V-CPL-02/02b.

### 6.3 Which `test/jacobian/` obligations are ported now

The design places this port in Phase 3; the Phase 1 gate needs the new core's Jacobian tests, so the obligations are written now under `test/newcore/jacobian/` and the old tests stay until cutover. AD removes hand differentiation, not wrong sparsity, state mapping, boundary coupling or residual equations, so each obligation becomes a containment or combined-physics test:

| Old test | Obligation | New form |
|---|---|---|
| `test_jacobian_coupling_regressions.jl` (strict per-entry criterion) | a missing small-magnitude block must fail regardless of problem size | the comparison helper keeps the per-entry rule with the row-magnitude floor; the synthetic unit test is ported as is |
| same (z-row × (T, ξ) block; WithChi θ_A back-reaction) | velocity rows depend on the cell's state through `θ` and on `θ_A` | containment at 6 cells with hot Arrhenius; WithChi probe with the `m_dry_solid_initial = 0` trap guarded |
| same (boundary mass flux on the energy row) | the boundary cell's energy row carries the face mass flux and the `∂T_s` chain | combined-physics containment: diffusive mass BC at both faces, Fickian and Darcy, convective plus re-radiative surface balance |
| `test_park_jacobian.jl` | exactness on masked states; ghost rows identity; masked-state Jacobian equals the compacted problem's on the live block | with step 6b, the same assertions on the new mask; otherwise Phase 2 |
| `test_regressions.jl` (22 chain-rule regimes) | each regime drives a specific coupling | the state builders are ported as containment states; the operator-level asserts are dropped |
| `test_sparsity.jl` | diagonal present on every row; profile shapes | pattern tests: diagonal on every row, masked rows diagonal only, bandwidth equality planar vs radial, color count against the Phase 0 formula |
| `test_geometry_ad_sparse_scatter.jl` | the pattern misses no geometry entry | subsumed by containment (node columns) |
| `test_kirchhoff_eps_eff_jacobian.jl`, `test_surface_radiation_eps_T_jacobian.jl` | the `a_eff`, `ε(T)` chains reach the energy rows | containment states with `ε(T)` components and per-component absorptivity; the surface-temperature column against central differences is kept |
| `test_darcy_venting_jacobian.jl` | the venting chain (φ_face, solid-basis κ, zero mobility, overpressure) | one containment state with all four features |
| `test_face_mean_guard.jl` | residual and Jacobian share one guard band | moot under AD; replaced by a trial-state test that the smooth face-mean guard is finite with finite derivative on both sides of the band |
| `test_chi_and_reduced.jl`, `test_structured_couplings.jl`, `test_structured_backend.jl`, `test_interface.jl`, `test_operator_*.jl` | operator-layer behaviour | not ported (the layer is deleted at cutover); the hand-written dense rows for `A`, `θ_A` get their own dense-AD comparison |

The performance invariants of `test/perf/` are ported in kind: zero-allocation residual after warm-up for the four configurations, a measured allocation budget for the colored Jacobian (DifferentiationInterface's prepared `jacobian!` is expected near zero; the budget is set from the measurement, with the ALE budgets of today, 6 KiB, as the ceiling), method count, JET type stability on the new modules, and the integrator-type independence test (`test/solver/test_integrator_type.jl`) on the new erased types.

### 6.4 Harness contracts moved into `Pkg.test()`

The harness stays as the Phase 0 record. Its contracts are re-expressed against the production types in `test/newcore/contracts/`: kinetics (`test_kinetics.jl`: finite rates and derivatives, `C¹` onset, rates equal the production kernel above the join, trial temperature in every kernel, absent reactant through a reacting step, approach to depletion at orders 0, 0.5, 1), caloric (`test_caloric.jl`, `test_nonlinear_caloric.jl`, `test_caloric_comparison.jl`), closure (`test_closure.jl`), transfer (`test_transfer.jl`), fields (`test_fields.jl`), energy (`test_energy.jl`), the DAE restart and shared-flux tests of `test_integration.jl`, and the callback controls of `test_restart_callbacks.jl` (an observer is numerically inert: bit-identical trajectory and counters with and without it). The solver and restart ladders (`test_solvers.jl`, `test_restarts.jl`) are measurements, not unit tests; they stay in the harness and are re-run on the new residual as gate scripts (§7, leftovers 2 and 3).

### 6.5 New tests Phase 1 adds

1. **Isothermal advection.** An isothermal composition profile with two species of different `c_p` under `w_override` (so `w ≠ v_s` and the relative fluxes are active) stays isothermal to round-off at every accepted step; the negative control (energy flux omitted) changes the temperature by more than 100 K. Generalizes `test_integration.jl`'s two-cell test and Phase 0 item 9 to the production kernels.
2. **Radial free-stream convergence.** A uniform concentration on a shrinking sphere and cylinder with a prescribed radial velocity stays uniform with an error that converges at the integrator's order under step refinement (Phase 0 item 8: 3.4e-5, 2.2e-6, 1.4e-7 at steps 0.2, 0.1, 0.05 on the sphere); planar stays uniform to round-off.
3. **Rigid rotation and history transport.** `diag(4, 1)` rotated 90° gives `diag(1, 4)`; a 45° fibre under `diag(2, 1)` ends at `atan(1/2)`; the perpendicular pair follows the declared priority; history transport under `w_override` with reactions off conserves mass-weighted means; a constant packed field under reaction with transport off stays constant.
4. **Energy row against the ledger.** With gas transport on (Fickian and Darcy, Eulerian and moving), the energy row expanded term by term (conduction, carried enthalpy, radiation, reaction heat, formation term) equals the per-step ledger built from inventories: exact to round-off in the inventory form, converging at the integrator's order in the temperature form; this is V-CPL-06's obligation on the kernel-contribution interface (C7).
5. **Both energy states.** The suite under `PYROLYSIS_CORE=new` runs twice, `energy_state = :temperature` and `:enthalpy`, with the §3.3 tolerance vectors.
6. **Constraint residuals.** The maximum scaled balance, intensity and cross-section residual over every accepted step of every case is recorded and must be below `rtol`.
7. **Containment.** §4.9, over every kernel family and physics option, at 6 and 8 cells.
8. **Determinism.** Repeated solves bit-identical for KLU and UMFPACK (`test/solver/test_repeated_solve_determinism.jl` in kind).
9. **Trial states.** Zero and negative trial temperatures, negative masses and volumes give finite values and derivatives in every kernel (the design's list, §3.2).

---

## 7. Milestones, gates, estimates

Single developer; the 1D suite stays green on the old core throughout (nothing in the old core changes except §2.3). Weeks are elapsed estimates against the design's 8 to 12; the total below is 12 to 14 because step 6b (C1) and the adapter are added and because the three Phase 0 leftovers are measured inside the phase.

| Milestone | Builds | Gate to pass before the next | Estimate |
|---|---|---|---|
| M0 Branch | Branch, first commit (§9), this plan, the direct dependencies, `test/newcore/` scaffold in the runner | Quality workflow green on a draft PR; harness quick cases pass from the committed tree | 2 days |
| M1 Geometry and state | `Meshes`, `States` (both energy forms, tolerance vectors, carried-field types), `Assembly` skeleton with the `sources` hook, the dense-AD oracle, a minimal `Integration.solve` on `KenCarp4` | §3.1 and §3.2 contracts; V-RX-01..03 in both energy forms through the adapter | 1.5 weeks |
| M2 Conduction and basic BCs | Conduction kernel, flux/adiabatic/convective/Dirichlet faces with the per-face Newton | V-HC-01..05 at their tolerances and convergence orders; containment for conduction and BCs | 1.5 weeks |
| M3 Reactions | Limiter type (C3), reactions kernel, formation term, switch-off | V-RX-01..04, V-CPL-04, V-CPL-05; depletion contracts; calibration guard (rates equal above the join) | 1 week |
| M4 Gas transport | Fickian and Darcy–Fick flux form, enthalpy face flux, mass and pressure BCs, venting, smooth property variants (C2) | V-GT-01/02, V-MMS-01; energy row against the ledger (test 4) in both forms; Darcy containment states | 1.5 weeks |
| M5 Radiation and full BCs | Lifted intensities, surface absorption, Kirchhoff and explicit absorptivity, blowing, transpiration options | V-RAD-01, V-RX-04a/b; `test/boundary_conditions/` ports; radiation containment states | 1 week |
| M6 Mesh motion | Closure interface and balance rows, `θ`, `G`/`M`, global rows, `w_override`, relative-flux machinery at `c = 0`, radial generators, carried `F` | V-CPL-03, V-MMS-02/03/04; tests 1, 2, 3; free-stream exact on planar; velocity containment incl. WithChi | 2 weeks |
| M7 Colored Jacobian | Pattern, coloring, DI closure, dense rows, containment suite, allocation gates | Containment green for every configuration; color counts match the Phase 0 formula; **leftover 1**: the lifted assembly timed against the proxy (§7.1) | 1.5 weeks |
| M8 Solver and deactivation | Three candidates behind the policy interface, `NoInit` sweep, `isoutofdomain` with the constraint check, observers, Newton cap, step 6b (1D transfer, collapse, mask, one-callback cascade) | V-CPL-02/02b; harness ladders re-run on the new residual; **leftover 2**: the integrator decision (§7.2); **leftover 3**: the finite-slab gate (§5.5) | 2 weeks (0.5 of which is step 6b beyond the stub) |
| M9 Outputs and the gate | `Outputs`, the adapter complete, report mode under `new`, both energy states, audit decks timed | The Phase 1 gate (§7.3) | 1 week |

Total: 12 to 14 weeks.

### 7.1 Leftover 1: the lifted assembly's measured cost

Phase 0 projected the colored Jacobian from the current residual's dual cost (0.32 to 0.73 residual-equivalents per color) and the lifted color count (15 to 36): 0.9 to 1.24× the structured assembly on Eulerian decks, 1.1 to 1.3× the ALE reduced assembly, 0.02 to 0.09× the ALE exact assembly with `c = 0`. The new residual has more rows per cell and two global rows, so the projection is replaced by a measurement at M7: on the audit decks (`eulerian`, `ale`, `ale_darcy` at 25, 100 and 400 cells, the `problem_for` decks of `scripts/mesh_phase0/baselines/integrator_observer_audit.jl`, which copies `scripts/performance_audit/audit.jl`), time one colored Jacobian of the new core against one `Structured` assembly of the old core (`DirectSolve` and `ApproxSolve`) on the same machine in the same session, with the residual-equivalent per color and the color count reported separately so the projection's two factors are each checked. The affordability decision (AD alone, or AD plus analytical fast paths for conduction and gas transport verified against AD) is made on this number and on the matched-accuracy solve times of M9.

### 7.2 Leftover 2: the integrator on the new residual

At M8 the harness ladders are not enough: the three candidates run the audit decks and the event-dense case (V-BO-04 if step 6b lands, otherwise V-CPL-02) on the new core with the colored Jacobian, each under its own policy, and are ranked by achieved error against a tight reference, constraint residual, steps, Jacobians and warm time. The candidate that qualifies on every deck with the best time is the Phase 1 default; `KenCarp4` raw with the flag and the solver-chosen restart is the fallback if none dominates. This is a gated decision (§10), recorded with the table in the gate report.

### 7.3 The Phase 1 gate, as written in the design

Every item is measured and reported in a gate document (`docs/guide/verification/` is not touched; the report is a plan-level document beside this one):

1. All verification cases of §6.2 within their existing tolerances on the new core, in assert mode, with the convergence spot checks (and V-CPL-02/02b with step 6b).
2. The audit decks within 1.5× of the old core's wall time at matched accuracy (target 1× by cutover): warm medians, same machine, same session, same tolerances, error against a tight reference.
3. MMS orders reproduced: V-MMS-01..04 spatial orders and the V-MMS-01 temporal sweep.
4. The radial generator's Jacobian has the same bandwidth and color count as the planar one at every cell count.
5. The suite passes once with `H` as the state through the state mapping.
6. Constraint residuals below tolerance at every accepted step on every case (test 6).
7. The energy-balance error in the temperature form no larger than the old core's ledger closure on the same decks at the same tolerance, and converging at the integrator's order under refinement; exact in the inventory form (test 4).
8. Achieved error against requested tolerance on every case, within its budget (a success return code is not evidence).
9. V-HC-03 and V-HC-02c at 320 cells within the current default's L2 error with the chosen integrator and its flag policy (§5.5).
10. Leftovers 1 and 2 recorded with their tables.

---

## 8. Risk register

Carried from design §5 (condensed to what bears on Phase 1) and extended with what the code review found.

| # | Risk | Where it shows | Mitigation / owner |
|---|---|---|---|
| R1 | Colored AD assembly slower than the hand-derived operators on the real residual (projection 1.1 to 1.8× the reduced assembly, measured on a proxy) | M7, leftover 1 | Measure both factors separately; analytical fast paths for conduction and gas transport verified against AD if the 1.5× wall-time gate fails; AD stays the definition of correctness |
| R2 | The hand-built pattern misses an entry (an upwind switch, the surface Newton's neighbour reach, the WithChi dense rows) | M7 | Containment against dense AD per kernel family and physics option at 6 to 8 cells, including trial states off the physical manifold (`restart_probe.jl`'s `:perturbed` and `:zero` states) |
| R3 | Algebraic rows upset the SDIRK estimate, or the Newton stops with the constraint unconverged | M8 | Rows stay in the norm; raw estimate; constraint check in `isoutofdomain` (C5); consistent sweep at `t₀` and every event; ladders re-run on the new residual |
| R4 | The integrator chosen on the ladder loses on the real decks, or a restart policy hides a miss behind a success code | M8, leftover 2 | Gate on achieved error and constraint residual, never on the return code; each candidate under its own policy; `KenCarp4` raw as the fallback |
| R5 | Inventory tolerances or the enthalpy inversion cause step-size trouble (82 Newton iterations and 2 to 8× the temperature form's time in the probe) | M1, M9 (gate item 5) | Warm start from the last accepted temperature; bracketed Newton; tolerance vectors of §3.3; the inventory form stays the alternative, not the default |
| R6 | A trial state outside the domain produces a non-finite residual or derivative | every kernel | The guarded-quantity list of §3.8; test 9; the quality grep for hard clamps |
| R7 | The vector field consumes an absent reactant | M3 | One-sided limiter; absent-reactant and per-order tests; V-RX-02d |
| R8 | Phase 1 scope grows (C1 step 6b, the adapter, the ledger re-authoring) | M8, M9 | Owner's call on C1 before M1; the adapter is the only consumer-facing code written in Phase 1; everything else is kernels and contracts |
| R9 | Properties and Materials edits (C2, C3) change old-core results | M3, M4 | Additive only; the old core's `Pkg.test()` runs after each; the calibration guard compares rates above the join to 1e-14 |
| R10 | `ε_m` as a mass fraction changes V-CPL-02's tolerances pinned under a thickness threshold | M8 (C9) | Report, do not re-pin silently; the owner decides |
| R11 | The pinned-mesh deactivation rule is not in the design (C4) | Phase 2 for V-BO-02/05; M8 for the interface | Owner question; the transfer-without-volume rule is the default proposal |
| R12 | The MMS source hook and the velocity override change the manufactured-solution driver's contract (C6) | M6 | `solve_mms(core = :new)` keeps the old path intact; sources are pure functions of `(z, t)` so the Jacobian stays exact; the `:fd` requirement is dropped only for the new core |
| R13 | V-CPL-06's ledger cannot be ported (C7) | M4 | Re-authored on the contribution interface; the old case stays with the old core |
| R14 | `Rodas5P` needs `∂f/∂t` and an exact Jacobian every step; cost scales with the production color count (2.2 to 4.8× the ladder's) | M8 | Finite-difference `tgrad` through the erased RHS; counted in leftover 2 |
| R15 | `FBDF` stalls when any observer leaves the flag set; a user callback that forgets to clear it reproduces the stall | M8 | The policy interface sets the flag for every package observer; user callbacks are wrapped with the policy's flag value |
| R16 | `NoInit` with an inconsistent state after an event silently continues | M8 | The sweep asserts the constraint residual after every event and throws above tolerance |
| R17 | DifferentiationInterface writes into a matrix with a different pattern than `jac_prototype`, or allocates per call | M7 | Pattern identity asserted at build; allocation gate; fall back to a hand-rolled colored ForwardDiff pass (the `restart_probe.jl` construction) if DI allocates |
| R18 | The two-core period lets the suites drift (a fix lands in one core only) | whole phase | Both suites run on every milestone; the adapter's keyword mapping is tested against `solve`'s signature |
| R19 | The `Internal` namespace and the top-level exports collide with the new modules (C10) | M0 | New modules are reached by module path until cutover; Aqua's undefined-exports and ambiguity checks run on every milestone |
| R20 | Event ordering under parameter perturbation makes finite-difference gradients non-smooth | Phase 2 | Carried; Phase 0 item 1 bounds what inverse analysis gets today |
| R21 | `ε_m` small enough for one event per layer is too small for the Newton | Phase 2 | Carried; start at 1e-2 |
| R22 | The correction fails to carry a layer through burnout or needs tuning | Phase 2 | Carried; `c = 0` plus deactivation is the fallback and is what Phase 1 ships |
| R23 | Verification tolerances do not transfer because the variables changed | M2 onward | QoIs compare quantities, not states; a miss is reported against the pinned tolerance with the observed value |
| R24 | 2D scope creep | whole phase | Only the interfaces are dimension-general; no 2D kernel, generator or test in Phase 1 |

---

## 9. Branch and process

### 9.1 Branch

`mesh-framework-phase1`, off `main` at 17bc477. Nothing from it goes to `main` until the owner is satisfied; a draft pull request exists from M0 so the Quality workflow (format, Aqua, JET) runs on every push, and it is converted to a real pull request at the Phase 1 gate. `Pkg.test()` and the `PYROLYSIS_CORE=new` run are local, as the test suite is today.

### 9.2 The first commit

It brings the design document and the Phase 0 harness under version control and nothing else. The `.gitignore` excludes both (`docs/*_YYYY-MM-DD.md` at line 49 covers the design document; `**/results/` at line 31 covers `scripts/mesh_phase0/results/`). Two ways were weighed: renaming (an undated design file and a `results/` directory renamed, which edits the harness's `RESULTS_DIR` and every `joinpath(@__DIR__, "results", …)`), or negation patterns. The negations are chosen because they leave the harness untouched and keep the design's own path, which the harness, the results summaries and this plan cite. Verified in a scratch repository with the same patterns: the three lines

```
!docs/Mesh_Framework_Design_2026-10-17.md
!scripts/mesh_phase0/results/
scripts/mesh_phase0/results/*.log
```

re-include the design document and the results directory, keep every other dated note ignored, and keep the run logs out (the logs are reproducible, 140 KB in total and up to 45 KB each; the `.md`, `.csv` and `.txt` evidence is committed, about 360 KB). The harness needs no formatting (JuliaFormatter check mode reports it clean) and runs on the committed project environment (§1.5). The literature survey the design cites (`reports/Pyrolysis mesh motion alternatives.md` and `research_notes/`, untracked and not ignored) is included in the same commit unless the owner keeps it local; the design document's line-pointer errata (§1.4, last paragraph) are applied in this commit, and nothing else in it changes.

Commit message, by the repository's convention (subject in the imperative, a body that says what and why, no attribution of any kind):

```
Bring the mesh-framework design and the Phase 0 harness under version control

The design document is the standing design for the new 1D core and the
harness is its executable contract set; both were local notes excluded by
the dated-note and results patterns. Negations re-include exactly these two
paths; run logs stay ignored.
```

### 9.3 Commit and review discipline

- One theme per commit; each milestone ends on a commit the owner can review, with the gate numbers of §7 in the milestone's message or in the gate report.
- No AI attribution anywhere: no `Co-Authored-By` trailer, no "generated with" line, in commits, the pull request, issues or files. This overrides any tool default.
- Replaced features are deleted, never deprecated: inside the new modules nothing is kept "for compatibility"; the old core is not edited except §2.3, and its deletion is Phase 3's one commit.
- No dates or historical notes in code, docstrings or the guides. The new modules' docstrings describe what the code does now. The guides are not edited in Phase 1 (the new core is not user-facing until cutover); when they are, the monoliths are edited and `scripts/split_guides.jl` regenerates `docs/guide/`.
- Formatting: `julia --project=scripts scripts/format.jl` before every commit (the formatter targets `src`, `test`, `ext`, `examples`, `benchmark`, `scripts`); the workflow runs `--check`.
- Quality: `julia --project=test/quality test/quality/aqua.jl` and `jet.jl` after any change to exports, dependencies or signatures; `test/perf/runtests.jl` in a fresh process for the allocation gates.
- Dependencies: `SparseMatrixColorings` and `DifferentiationInterface` added with `[compat]` pins to the installed versions (0.4.28, 0.7.21); `OrdinaryDiffEqBDF` and `OrdinaryDiffEqRosenbrock` are reached through the umbrella as today; the test and quality manifests are regenerated and committed when they change.
- The harness is not modified by this effort. Re-running its ladders on the new residual is done by new scripts under `test/newcore/gates/` that `include` the harness probes, so the Phase 0 evidence is never overwritten.
- Gates per milestone are targeted, not the full suites. (1) No old-core `Pkg.test()` when the milestone changes no old-core source file (`git diff --stat <last milestone> -- src` lists only the new modules and the include line); the old-core cases whose test files changed are run through the case runner on the old core instead, in the same session as their new-core numbers. A milestone that edits the old core (§2.3) runs its before/after fingerprint and the old-core tests of the touched layer, and the full old-core `Pkg.test()`. (2) The new core's tests and `NEW_CORE_CASES` run in persistent sessions that revise the tree instead of recompiling it, split across a few processes, both energy forms. (3) The containment matrices of `test/newcore/` are trimmed to what a milestone's kernels need (both sizes and every state for the configurations the milestone adds; one size and the physical and trial states for the rest). (4) The full three-way run (old core, `PYROLYSIS_CORE=new` in each form) is the Phase 1 gate's (M9), and runs before the draft pull request becomes a real one. Quality (format, Aqua, JET) runs at every milestone.

### 9.4 What the owner decides before M1

1. C1: step 6b in Phase 1 (plan as written) or the design's split (seven cases move to the Phase 2 gate).
2. C4: the pinned-mesh deactivation rule.
3. The literature survey in the first commit or not.
4. Module names of §2.2 (provisional).

---

## 10. Decisions carried as gated choices

None of these is resolved by this plan. Each is carried with its Phase 1 default, the alternative, the measurement that decides it and the milestone where that measurement is taken.

| Decision | Phase 1 default | Alternative | Decided by | When |
|---|---|---|---|---|
| Energy state | Temperature (`TemperatureState`) | Enthalpy inventory (`EnthalpyState`) | Gate items 5 and 7: both forms run the suite; the temperature form's ledger error must be no larger than the old core's closure and converge at the integrator's order; the inventory form's cost (Newton iterations, wall time) on the audit decks; the default is picked by the owner from that table, and the inventory form is what the relative-flux work of Phase 2 uses | M9 |
| Integrator with restart and per-step-flag policy | `KenCarp4`, raw estimate, flag kept, solver-chosen restart | `Rodas5P` (either policy, exact Jacobian), `FBDF` (observer flag cleared, forced restart step) | Leftover 2: qualification on every audit deck and the event-dense case by achieved error and constraint residual, then time; the finite-slab gate (§5.5) with the chosen policy | M8 |
| Lifting against bordering | Lifting (algebraic unknowns, local rows, colored AD) | Bordered solve of the eliminated system (exact, every coupling enumerated by hand) | Leftover 1 and gate item 2: colored assembly cost and matched-accuracy solve time within 1.5× of the old core; if lifting fails on cost and the pattern is the cause, the bordered form is prototyped on the same decks before any analytical fast path is written | M7, M9 |
| The correction `c` | `c = 0` (soft activation and relative fluxes implemented, inactive) | Activated correction with the MMPDE monitor | Phase 2's gate (one event per layer or fewer on V-BO-01 across the ladder); Phase 1 only proves the machinery is inert at `c = 0` (charring decks reproduce) and exact under `w_override` (isothermal advection) | Phase 2 |
| Analytical fast paths | None | Conduction and gas transport operators verified against AD | Leftover 1 | M7 |
| Deactivation threshold `ε_m` | 1e-2 | Tightened under the Phase 2 gate | V-CPL-02 and the V-BO ladder | M8 (first value), Phase 2 |
| Tolerance-vector scales (§3.3) | As listed | Per-case overrides | Gate item 8 | M9 |

---

## 11. Summary

**What Phase 1 builds.** A second 1D core beside the existing one, in new modules (`Meshes`, `States`, `Kernels`, `Closures`, `Deactivation`, `Assembly`, `ColoredJacobian`, `Integration`, `Outputs`) that import the material, property, boundary-condition and physics layers unchanged except for four additive edits (a one-sided depletion limiter type, smooth property floors beside the hard ones, two direct dependencies, and a pure move of the face-local boundary code into BoundaryConditions). The integrator carries per-cell component masses, an energy state (temperature by default, the enthalpy inventory behind the same two-function state mapping), node positions, algebraic nodal velocities, lifted Beer–Lambert face intensities and the two cross-section scalars; kernels are face-based from the first line and return flux and source contributions; the Jacobian is colored sparse automatic differentiation on a hand-built, containment-verified pattern; the integrator candidates sit behind one interface with their restart and per-step-flag policies, `NoInit` with the project's consistent sweep, the full error norm and a constraint check that rejects steps through `isoutofdomain`. With the proposed step 6b it also carries the 1D terminal deactivation event (transfer, collapse, masked rows, one-callback cascades) so that the Landau ablation cases run; the correction stays off. Every existing problem definition is accepted through an adapter; the verification suite runs on either core by one environment variable.

**What it costs.** Twelve to fourteen weeks of one developer (the design said 8 to 12; the difference is the deactivation event, the adapter and taking the three Phase 0 leftovers inside the phase), the ported Jacobian obligations of §6.3, the harness contracts of §6.4 and the nine new tests of §6.5 under `test/newcore/`, and no change to the old core's behaviour. The first commit brings the design document and the Phase 0 harness under version control with two gitignore negations.

**What its gate proves.** Every exact-solution and coupled case of the suite within its existing tolerance on the new core, in both energy states; the manufactured-solution orders; constraint residuals and achieved errors measured, not inferred from return codes; the energy balance no worse than today's ledger closure and converging at the integrator's order; the radial Jacobian as banded as the planar one; the 320-cell finite-slab accuracy reproduced with the chosen integrator and policy; the audit decks within 1.5× of the old core's wall time; and the three things Phase 0 could not settle, measured on the real residual: the colored assembly's cost against the proxy projection, the integrator with its two policies, and the finite-slab accuracy.

**What remains unproven after it.** Burnout: whether the correction `c` brings a vanishing layer to one terminal event (or none for residue formers), its tuning, and the finite-difference gradient smoothness across former park times (Phase 2). Cost at 1× rather than 1.5× (cutover target). The inventory form's cost on production decks beyond the audit set. Rosenbrock cost at production color counts on event-dense decks. The pinned-mesh deactivation rule and the fixed-mesh burnout cases. Everything beyond 1D: the potential and mechanics closures, element death on conforming meshes, the cross-diffusion terms, cone geometry, and the 3D preconditioner (Phase 4). The old core's deletion and the guides (Phase 3).

---

## 12. Implementation status

Branch `mesh-framework-phase1`, draft pull request #12. This section is the record of what has landed against §7; the sections above are the plan as written, and where an item below departs from them, this section governs.

### 12.1 Milestones

| Milestone | Commits | Gate evidence |
|---|---|---|
| M0 | `5a2e5b3` (design, harness, survey, errata), `92e37ab` (this plan) | Quality workflow green on #12; harness quick cases (kinetics, caloric, closure, transfer, fields, energy, integration, nonlinear_caloric) pass from a clean checkout of the committed tree, 1,126 assertions; old-core `Pkg.test()` 4113/4113 in 36 min 28 s on the reference machine |
| M1 | `cfeb6f1` (guarded face mean), `edbb40a` (Meshes, States), `af0d861` (Kernels, Assembly, ColoredJacobian, dependencies), `832a26d` (Integration, Outputs, suite routing) | V-RX-01, V-RX-02a-d, V-RX-03 pass every QoI in both energy forms; containment of the dense-AD pattern and colored = dense to 1e-12 per entry with the row floor, at 6 and 8 cells, on the physical, perturbed, zero and negative states; 12 colors (3 × (NC + 1)) on the fixed mesh; zero-allocation residual in both forms; enthalpy form conserves summed enthalpy to round-off and the temperature form converges at fourth order (two-cell closed test); old-core `Pkg.test()` 5263/5263 (41 min 27 s); `PYROLYSIS_CORE=new` 1591/1591 in each energy form; Aqua 9/9, JET clean; Quality workflow green |
| M2 | `a75a9dc` (move of the face-local boundary code, byte-for-byte apart from docstring references), `45ffd4e` (unclamped mixture ε, a), `f65d848` (surface balance per face, outputs, V-HC routing) | V-HC-01, V-HC-02a-d, V-HC-03, V-HC-04a-c, V-HC-05 pass every QoI and spatial convergence check in both energy forms (orders 1.99 to 2.00); containment and colored = dense with the surface Newton on convective, radiative and combined faces; zero-allocation residual with those faces; old-core `Pkg.test()` 5508/5508 (41 min 18 s); `PYROLYSIS_CORE=new` 1926/1926 (temperature), 1926/1926 (enthalpy); Aqua 9/9, JET clean |
| M3 | `73200f8` (one-sided limiter in Materials and Physics), `0e85ee7` (the reaction kernel on it, the kinetics contracts), `b2c9228` (`core_solve`, case routing) | Old-core rates bit-identical: the hash of 1,170 rate values and gradients (orders 0 to 2, three thresholds, limiter on and off, concentrations from −1 to 1e3) and of eight reacting solves (V-RX-01, 02a, 02c, 02d, 03, 04a, V-CPL-01, 04) unchanged by the limiter commit; kinetics contracts 1,209/1,209 (calibration guard to 1e-14 above the join, the kernel's sources equal to the material limiter's, absent reactant exactly absent through a reacting solve at orders 0, 0.5, 1 in both forms); V-RX-01 to V-RX-04b and V-CPL-04 with its routed criticality sub-solves pass every QoI in both energy forms, and on the old core in the same session; the old-core `Pkg.test()` of the M4 row covers this tree |
| M4 | `82ecf3f` (face-mean guard argument, skeletal volume fraction), `36c1917` (gas transport, enthalpy face flux, outputs, tests), `17ca282` (`solve_mms(core = :new)`, case routing) | V-GT-01 (orders 2.01, 2.00), V-GT-02 (P² error 7.08e-4 against the pinned 7.1e-4; orders 0.89, 0.94), V-CPL-05 (thin-limit deviation 1.42e-3, L-rates 2.004, 2.002) and V-MMS-01 (T error 0.0187 and 0.0135 K against a 0.05 K tolerance; orders 1.98, 2.00) pass every QoI in both forms; containment and colored = dense on Fick (diffusive, film and flux mass conditions) and Darcy (pressure venting with zero-mobility faces; combined film, convective-pressure and Neumann-pressure conditions) configurations with blowing on the surface, at 6 and 8 cells, both forms, five states; zero-allocation residual; exact isothermal transport; energy row against the ledger term by term to 1e-11 and over solves (inventory form 1e-13 of the ledger; temperature form orders 2.1, 2.5 Fick and 2.7, 3.3 Darcy under fixed-step refinement); new-core tests 3,162/3,162; old-core `Pkg.test()` 7275/7275 (49 min 08 s, on a shared machine); `PYROLYSIS_CORE=new` 3746/3746 (temperature), 3746/3746 (enthalpy); Aqua 9/9, JET clean |
| M5 | `ce9ea32` (lifted Beer–Lambert intensities, surface absorption, outputs, tests, V-RAD-01 routing) | V-RAD-01 passes every QoI in both forms (L∞ at 90 s 0.0649 K, the old core's 0.0649; orders 1.996, 1.999); V-RX-04a/b unchanged; containment and colored = dense with lifted intensities under Kirchhoff (ε(T), a component absorptivity) and explicit time-dependent absorptivity, with and without gas transport, at 6 and 8 cells, both forms, five states, at most 3(NC+1) + 2 colors; consistent initial intensities; the beam's cell sources equal `I_up (1 − e^{−τ})` and telescope to the transmitted intensity to 1e-12; the exposed cell's temperature and mass columns equal central differences; ported surface-condition guards (no state carried between evaluations, ambient ingress, two-face venting closing against the inventories to 1e-3, Dirichlet face temperature, bottom incident flux refused under Beer–Lambert, surface absorption equal to Beer–Lambert's total when nothing is transmitted, a time-varying flux sampled every evaluation, blowing lowering the convective loss); new-core tests 3,653/3,653; old-core `Pkg.test()` 7766/7766 (51 min 47 s, shared machine); `PYROLYSIS_CORE=new` 4252/4252 (temperature), 4252/4252 (enthalpy); Aqua 9/9, JET clean |
| M6 | `53e2031` (closure, moving mesh, cross section, relative fluxes, outputs, tests), `01a2ee3` (`solve_mms(core = :new)` on moving meshes, case routing) | V-CPL-03 (thickness ratios within 2.6e-7 of the exact endpoints, total mass to 1.6e-15) and V-MMS-04 (T 2.05e-5 K, solid 6.61e-7, gas 2.299e-3, surface position 1.633e-5, orders 2.03, 2.00; the old core 2.03e-5, 6.57e-7, 2.299e-3, 1.633e-5 in the same session) pass every QoI in both forms; V-MMS-02 and V-MMS-03 pass all but one QoI each (C24, open); velocity containment and colored = dense on planar, cylindrical, spherical, variable-cross-section (the dense rows against dense AD), Beer–Lambert, Fick and Darcy moving meshes and under a prescribed velocity with the relative fluxes, at 6 and 8 cells, both forms, five states; zero-allocation residual; test 1 (isothermal advection to 1e-9 K in both forms, more than 100 K without the enthalpy flux), test 2 (planar free stream within 1e-13 at fixed steps in both forms; radial at orders 3.84, 3.99 on the cylinder and 3.86, 3.99 on the sphere), test 3 (rotation, fibre and priority from M1; history transport conserving the history sums to 1e-12; the packed field per D15); new-core tests 5,291/5,291; old-core `Pkg.test()` 9404/9404 (80 min 42 s, shared machine); `PYROLYSIS_CORE=new` 5918/5918 (temperature, 47 min 13 s), 5918/5918 (enthalpy, 44 min 28 s); Aqua 9/9, JET clean |

Effort to date against §7: M0 to M6 plus about half of M7 (D3, D17), about eight and a half of the thirteen to fourteen weeks.

### 12.2 What is built

- **Module order.** `Materials → Properties → BoundaryConditions → Geometry → Physics → Meshes → States → Kernels → Closures → Assembly → ColoredJacobian → Adaptivity → Discretization → Problem → Outputs → Integration → Residual → …`. The kernel-level modules load before Discretization and Problem, so they cannot import them; Outputs and Integration load after Problem (they need `PyrolysisProblem`) and before Residual. `Deactivation` does not exist yet (M8).
- **Meshes.** `Face{D,K}`, `Cell{D,F,K}`, `MeshTopology{D,C,FT,CT}` (coordinate system, faces with owner/neighbour/tag, cells, node count, boundary tags, pinned and Lagrangian nodes, cross-section measure), `Mesh{D}` (topology plus reference positions), `line_mesh` for planar, cylindrical and spherical 1D, `MeshGeometry{D,T}` written only by `geometry!`. Normals point from owner to neighbour, outward on boundary faces. In 1D cell `i` lies between faces (and nodes) `i` and `i + 1`; face 1 is `:bottom` (owner cell 1, normal −1), face `N + 1` is `:top`.
- **States.** `StateLayout{NC,E,Mov,Rad,Chi}` with the interleaved node/cell/global ordering; a fixed mesh without intensities carries no node block. Index functions `mass_index`, `energy_index`, `history_index`, `position_index`, `velocity_index`, `intensity_index` (face `f` sits at node `f`), `area_index`, `lateral_rate_index`, `cell_rows`, `node_rows`; `mask` and `differential` BitVectors; `mass_matrix`. `StateSetup` carries the `GasTransport` model (`NoGasTransport`, `FickTransport`, `DarcyFickTransport`). `DerivedState{NC,T}` filled only by `derive!`: masses, concentrations and positive-part concentrations, temperatures and their trial clamp, the floored heat-capacity sum, sensible enthalpies, conductivity (smoothly floored), ε, a, and as the physics needs them face intensities `I` and absorption coefficients `α`, gas-transfer coefficients `λ`, the floored porosity `φ`, pressure `P` and mobility `mob`; on a moving mesh the node positions and velocities `x`, `w`, and under a variable cross section the area `A`, its lateral rate `θ_A` and the history pair `G`, `M`. `CaloricLaws`, `CaloricExtension`, `invert_caloric`, `CellReference`, `abstol_vector` (mass, energy, intensity, position, velocity, history and cross-section rows), `CarriedField`, `CONVERSION_HISTORY`, `material_axes`, `physical_tensor`.
- **Kernels.** `Contributions{NC,T}` (`F_mass`, `F_energy`, `S_mass`, `S_energy`, `θ`, `M_flux`, `T_face`, `r_intensity`, `r_velocity`, `F_carried`, `S_carried`, reaction scratch). In assembler order: `ReactionKernel` (a `OneSidedLimiter` from the material's threshold and flag, trial temperature); `GasTransportKernel{M}` (Fick and Darcy–Fick fluxes in flux form, mass-condition patches per tag, the carried enthalpy at the distance-weighted face temperature, `outward_darcy_velocity`); `BeerLambertKernel` (exposed-face and per-cell intensity rows, cell sources `I_up A_up − I_down A_down`, `sweep_intensities!`, `incident_intensity`, `transmitted_face`); `ConductionKernel`; `ThermalBoundaryKernel` (closed forms, the surface Newton under the radiation model, the enthalpy the boundary mass flux carries at the face temperature); `RelativeFluxKernel` (last, under a prescribed mesh velocity only: every component at `(vₛ − w)·n` times its face concentration, the history pair at the solid mass flux times its per-solid-mass face values, the enthalpy `Σⱼ Fⱼ hⱼ(T_f)` of the same fluxes at the reconstructed face temperature; face states the upwind value plus the smooth `van_albada` slope, the downwind gradient where no cell lies beyond, one-sided extrapolation on boundary faces).
- **Closures.** `Closure`; `BalanceRowClosure{O,X}`: the pinned node's row `w_p = 0`, cell `k`'s balance `Σ_f A_f (w·n)_f − (θ_k − θ_A) V_k` on the velocity row of its node farther from the pinned node, `θ_k` by `Physics.volume_change_rate` on the reaction mass sources at the trial temperature; the `override` form `w_n − w(x_n, t)`; `closure_rows!` (θ, the history sources, the velocity residuals; after the kernels and before the external sources, so θ is the reactions' alone), `consistent_velocities!` (one bidiagonal sweep), `rank_check`. `LateralShrinkage` (law, `A₀`, dry-gas flags, `Σ M` floor), `chi_bar`, `lateral_rows` (both rows from the six column sums `Σ G, Σ M, Σ Ġ, Σ Ṁ, θ_A, A`), `consistent_lateral`.
- **Assembly.** `CoreDefinition` (with the closure, `nothing` on a fixed mesh), `CoreWorkspace`, `residual!` (derive → kernels → closure → external sources → rows), `assemble_rows!` (the energy row in either form, the intensity rows from `r_intensity`, the position rows `dxₙ/dt = wₙ`, the velocity rows from `r_velocity`, the history rows, the area and lateral-rate rows, dead rows), `NoSources`, `VolumetricSources`.
- **ColoredJacobian.** `jacobian_pattern(def)` (cell rows: cells i−1..i+1, the intensities of their two faces, the positions of nodes i−1..i+2 on a moving mesh, with relative fluxes cells i−2..i+2, positions i−2..i+3 and their two velocities, the area; intensity rows: intensities and positions of nodes f−1..f+1 and cells f−1..f; position rows: their velocity; velocity rows: their node, and a balance row its cell's other node, masses and energy, the area and the lateral rate; the area and lateral-rate rows dense; diagonal everywhere), `dense_rows`, `LateralRows`, `ColoredJacobianCache` (the dense rows colored by their diagonal, the colored result scattered into the full pattern, the dense rows written after), `jac_prototype`, `jacobian!`, `dense_jacobian`, `missing_entries`, `compare_entries`.
- **Integration.** `setup` (the §2.4 adapter: gas model from the material and mass conditions, kernels, the moving mesh under `use_ale` with the back node pinned, the cross section under a lateral law, a prescribed `w_override`, radial `coordinates`, the tolerance vector with `time_scale`), `consistent_state!` (intensities by the recurrence, then the lateral rows from the column sums and the velocity sweep after a kernel pass, the velocity rows checked at round-off), `integrate(cp)` and `solve = setup + integrate` (adds `dt`, `adaptive`, `min_thickness`: a moving mesh's solve ends `:Terminated` below it; a mass-matrix DAE once algebraic rows are present), `IntegratorPolicy`, `default_policy`. `inadmissible` adds crossed nodes and a non-positive area. Refused with an `ArgumentError`: `P1_QUASI_STEADY`, an incident flux on any tag but `:top` under `BEER_LAMBERT`, a lateral law without `use_ale` or off planar coordinates, a law not equal to one at zero conversion or not positive on `[−w, 1 + w]`, a lateral law without dry solid, a prescribed velocity without `use_ale` or on a reacting material.
- **Outputs.** `CoreSolution` adds `boundary_mass_flux` and `boundary_blowing_ratio` per tag and `extras` (`surface_mass_flux`, `surface_blowing_ratio`, `radiation_transmitted`, `volumetric_absorption`, `thickness_history`, `cross_section_area_history`, `total_mass` per initial `:top` area) under the old schema's names; `z` holds the node positions on a moving mesh.
- **Old core (the §2.3 changes so far).** `guarded_harmonic_mean` (now with a guard argument), `mixture_emissivity`, `mixture_absorptivity`, `skeletal_volume_fraction` in Properties; `AbstractDepletionLimiter`, `OneSidedLimiter` and the limiter type parameter of `Material` in Materials, `reaction_rate`'s per-reactant dispatch and `evaluate_reactions!`'s limiter argument in Physics; the two direct dependencies; the move of §2.3.
- **Suite.** `PYROLYSIS_CORE=new` routes the cases of `NEW_CORE_CASES` (V-RX-01..04, V-CPL-03, V-CPL-04, V-CPL-05, V-GT-01/02, V-MMS-01, V-MMS-04, V-RAD-01, V-HC-01..05) through `Integration.solve` in the energy state of `PYROLYSIS_ENERGY_STATE`; `core_solve` routes every case-level sub-solve the same way; `solve_mms(…; core, energy_state)` solves on either core, on the new one with `use_ale`, `w_override` (Eulerian sources) and `motion` (material-frame sources). V-MMS-02 and V-MMS-03 pass the suite's core but are not in `NEW_CORE_CASES` (C24).

### 12.3 Departures from this plan

- **D1. Dependencies landed in M1, not M0.** Aqua's `test_stale_deps` fails on a declared dependency nothing loads, so SparseMatrixColorings and DifferentiationInterface came with the first code that imports them.
- **D2. Conduction and closed-form boundary faces landed in M1.** V-RX-03's QoI is the quasi-steady drift `T_c − T_bc = Q̇ Δz²/(4k)` of a cell between Dirichlet faces, so M1 needed conduction to Dirichlet faces. M2 added the surface Newton.
- **D3. The pattern, coloring and colored closure landed in M1.** M7 keeps the moving-mesh, global rows of the pattern, the dense hand-written rows, the allocation budget of the colored Jacobian, the color-count check against the Phase 0 formula and leftover 1. The intensity rows landed in M5.
- **D4. C3 refined (applied in M3).** `Material` gains a limiter type parameter instead of an abstract field, and the new core converts the material's limiter to a `OneSidedLimiter` with the same threshold and enabled flag. The finite base is formed in log space below zero: the handoff draft's `δ log1p(exp(ξ/δ))` underflows to zero before its slope near `ξ/δ ≈ −730`, which made the power's derivative `0 · ∞` for fractional orders; with the limiter enabled the rate short-circuits to zero at and below zero, where it and every derivative vanish.
- **D5. The move of §2.3 is narrower and has one load-order change.** See C13 to C15.
- **D6. Case-level sub-solves (applied in M3).** `core_solve` in `casedefs.jl`; V-CPL-03, V-CPL-04 and V-CPL-05 call it in their QoIs and report curves.
- **D7. Milestone membership (applied).** V-RX-04a/b and V-CPL-04 joined in M3, V-CPL-05 in M4.
- **D8. Porosity admissibility moves to M8.** §12.5 put the porosity check in M4. In the inventory form the porosity needs the temperature (intrinsic densities are functions of it), which `inadmissible` can only get from the inversion; the M8 predicate evaluates the residual at the proposed state for the constraint check (C5) and gets it there. M4 floors the porosity smoothly in the residual and adds no admissibility test.
- **D9. The reaction switch-off near `ε_m` lands with deactivation (M8).** It does nothing without deactivation.
- **D10. No transpiration option to port.** `use_transpiration_bc` was deleted from the old core before Phase 1; M5's "transpiration options" are empty.
- **D11. `SURFACE_ABSORPTION` is the surface balance of `NO_RADIATION`.** Both put the incident flux into the surface Newton at the exposed cell's absorptivity (or the condition's), so the new core accepts both and their residuals are identical.
- **D12. A conductivity floor.** The cell conductivity is floored smoothly at 1e-6 W/(m·K): a trial state with every mass at or below zero has none, and the closed-form face temperature of a heat-flux face divides by it (an infinite face temperature then made the carried enthalpy `0 · ∞`). Not in the §3.8 list.
- **D13. The potential closure is not built.** §2.2 and §3.5 give the closure interface a 1D implementation of the potential form that must reproduce the balance rows. Its unknowns are cell multipliers, a layout the 1D core does not have, and in 1D the balance rows alone fix `w`; M6 builds the interface (`Closure`, `closure_rows!`, `consistent_velocities!`, `rank_check`) and `BalanceRowClosure`, and the potential form lands with Phase 4c, where its layout is needed.
- **D14. Under a prescribed mesh velocity the solid is at rest.** The override rows replace the balance rows that would give the solid velocity, and without them `vₛ` is the swept strain integral, a dense coupling. `setup` therefore refuses an override on a reacting material, `vₛ = 0`, and every inventory, gas included, crosses the faces at `−w·n` through `RelativeFluxKernel`. The old core under its override kept the species on the mesh and advected the temperature alone; on the new core every V-MMS-02/03 source is the Eulerian one (`solve_mms` header).
- **D15. The deformation gradient is the axial stretch, with no state.** At `c = 0` a cell is a material volume, so `F = h/h(0)`; under the override the solid is at rest and `F = 1`. Nothing in 1D reads `F`, and no packed (intensive) carried field exists in Phase 1, so the "constant packed field under reaction" part of test 3 has nothing to run on; it lands with the first packed field. The rest of test 3 runs (rotation, fibre, axis priority from M1; history transport under the override in M6).
- **D16. The lateral rows.** The conversion is clamped into `[−w, 1 + w]` (`w = 1e-3`), the identity on `[0, 1]`, as every other guard is the identity on physical states; the handoff's "clamped into [0, 1]" would distort the law near zero conversion, where every column starts. The lateral rate is the logarithmic rate of the area row along `(Σ Ġ, Σ Ṁ)` (the quotient rule with `Ṁ` and the clamp's slope), so the area row's own rate is `θ_A A` exactly; at `c = 0` (`Ṁ = 0`) and on physical states it is §4.8's `law'(χ̄) χ̄̇ / law(χ̄)`. `setup` requires `law(0) = 1` (the old core accepted any law and started at `A₀ law(0)`), a law positive on the clamp band, and dry solid in the column (`Σ M > 0`; the `m_dry_solid_initial = 0` trap of §6.3 is refused, and the containment states set `χ̄ > 0`).
- **D17. The dense rows by the chain rule through six sums.** Both rows depend on the state only through `Σ G, Σ M, Σ Ġ, Σ Ṁ, θ_A, A`; `LateralRows` differentiates the six-number function, takes `∂Σ Ġ/∂u` and `∂Σ Ṁ/∂u` as column sums of the history rows of the colored pass (exact there), and writes both rows into the full pattern, against dense AD in the containment tests. The `A` and `θ_A` columns stay in the colored pass (a dense column costs one color; only a dense row breaks the coloring), where the handoff had them by hand too.
- **D18. `min_thickness` ends the solve.** As on the old core, a continuous callback terminates a moving mesh's solve (`:Terminated`) when the thickness falls below it (default 1e-6 m), rather than an inert keyword waiting for M8.
- **D19. Radial meshes through the adapter.** `setup(…; coordinates = Cylindrical() | Spherical())` reads the problem's node positions as radii with the default cross section (unit axial length, 4π), so the radial generators run through the same adapter as the planar cases (tests 2 and the bandwidth check).
- **D20. V-MMS-04's ApproxSolve companion runs on the old core only.** The new core has one Jacobian; its two companion QoIs would repeat the main solve.
- **D21. Targeted gates per milestone (owner, after M6).** The full old-core and two-form new-core suites took 81, 47 and 44 minutes at M6 on the shared machine, besides recompiles of tens of minutes, while the old core was not edited. From M7 the per-milestone gate is §9.3's targeted procedure; the full three-way run belongs to the Phase 1 gate. The R9 rule (the old core's `Pkg.test()` after an edit of Materials or Properties) still holds for old-core edits.

### 12.4 Conflicts found during implementation

- **C13. Load order.** `BoundaryConditions` loaded before `Properties` (`src/Pyrolysis.jl` at 17bc477, lines 63-74), and `blowing.jl` imports `Properties.chapman_enskog_diffusivity` (`blowing.jl:12`). `Properties` imported `BoundaryConditions` and `Geometry` without using either (`Properties.jl:14-15`). Resolved: those two imports deleted and `Properties` loaded before `BoundaryConditions`, the one non-additive edit to the old core so far.
- **C14. `boundary_fluxes.jl` is not mesh-free in whole.** `_bc_boundary_face_fluxes` calls `Discretization.boundary_darcy_velocity` with `Val(tag)`. It stays in Discretization; the new core has its own `Kernels.outward_darcy_velocity` (M4), tested equal to the tagged kernel along the outward normal on every pressure condition and a combination.
- **C15. The surface Newton's diagnostics argument is typed on `BoundaryState`.** `BoundaryState` moved from `Geometry/mesh.jl` to `BoundaryConditions/boundary_state.jl`.
- **C16. §3.1's invariant `Σ_f A_f n_f = 0` holds only for planar 1D.** The volume identity `V = (1/d) Σ_f A_f (x_f · n_f)` holds in all three geometries and is what the tests assert.
- **C17. SparseMatrixColorings checks only size and stored-entry count before decompressing.** `ColoredJacobianCache` asserts pattern identity with `jac_prototype`; the integrator's `jac!` checks it on every call.
- **C18. C10 in practice.** `using Pyrolysis` exports old-core names the new modules also define; test files import the new names explicitly.
- **C19. `get_cell_positions` is typed on `PyrolysisSolution`.** The cases call the suite's `cell_positions` (both solution types).
- **C20. Guides cite moved paths.** Technical Reference chapters 7 and 12 name `src/Discretization/boundary_temperature.jl` and `Discretization/blowing.jl`; fix when the guides are next regenerated.
- **C21. The upwind switch of the Darcy advection against "no non-smooth limiters in the residual".** Design §3.10 excludes "hard clamps and non-smooth limiters anywhere in the residual"; plan §4.5 and the M4 handoff specify the first-order upwind `(ξ/φ)_upwind v_g`, which the new core implements as the old one does (`src/Physics/gas_transport.jl:107`, `ξ_upwind = v_g >= 0 ? ξ_L : ξ_R`; the vent's outflow-only rule at `:231-232`; new core `src/Kernels/gas_transport.jl`, `_advective_flux` and `_vent_rate`). The flux is continuous with a kink in its derivative at zero velocity. Measured on a Darcy char-layer column whose flow reverses near t = 0.4 s (cold back face): fixed-step `KenCarp4` fails Newton at t = 0.40 to 0.425 s for every step in {0.05, 0.025, 0.0125} and tolerance in {1e-8, 1e-10}, on the old core and the new one alike; the same column with an adiabatic back, whose flow never reverses, converges at every step from 0.05 down; adaptive solves recover by rejecting steps (the reversing column: 302 accepted, 13 rejected, 1 nonlinear failure). Not resolved: the plan's scheme is kept, the ledger test uses the non-reversing column, and the alternative, a smooth upwind weight set by the face Péclet number `v d / λ_f` (central where diffusion dominates, upwind where advection does, with no new scale), changes the discretization and is for the owner, best decided with the M8 integrator measurements, which count nonlinear failures.
- **C22. A NaN trial inventory aborts an inventory-form solve.** `invert_caloric` throws on a non-finite inventory (`src/States/caloric.jl`, `_invert_primal`; pinned by `test/newcore/test_states.jl:244`), so a NaN stage value ends the solve with an exception, where the temperature form returns a NaN residual and the integrator rejects the step. The ported guard test (`test_surface_bc_guards.jl:56` → `test/newcore/test_boundary_conditions.jl`) drives the inventory form with a wild finite inventory only. Proposed for M8 (step rejection and the Newton-failure cap): return a NaN temperature on a non-finite inventory so both forms reject alike.
- **C24. Two moving-mesh QoIs are pinned at the old core's state-level exactness.** V-MMS-02 (`test/verification/cases/V-MMS-02.jl`, "free-stream ξ_gas relative drift", tolerance 1e-12, observed 1.5e-16 on the old core) and V-MMS-03 (`V-MMS-03.jl`, "relative L∞ ξ_solid (mesh-bound) error", tolerance 1e-10, observed 2.1e-11) were pinned on a core whose states are the concentrations: a uniform concentration has a zero rate there and stays fixed under every Newton iterate. On the new core the states are the masses and the positions (design §3.2), so a uniform concentration is the relation `m = ξ V(x)` between two kinds of state. Measured on the new core in this session: V-MMS-02's gas drift is 1.8e-12 (temperature form) and 1.2e-10 (enthalpy form), and with converged stage solves (`NLNewton(κ = 1e-6)`) 2.5e-15, so the discrete GCL is exact and the drift is the Newton leftover the default `κ` allows on the gas rows, whose tolerance is set by the cell's solid mass (§12.5; 4e-8 of the gas itself here). It depends on the Newton path (κ = 1e-3 gave 2.2e-9). V-MMS-03's uniform solid misses its pin at 5.6e-10 (temperature) and 2.3e-10 (enthalpy) against the solve's reltol of 1e-8: with converged stage solves it is 1.03e-10 and at reltol 1e-10 (225 steps against 64) 1.3e-11, so here the time error of the moving positions, which the concentration `m/V(x)` carries, is part of it. Every other QoI of both cases passes on the new core, in both forms. At reltol 1e-10 and abstol 1e-12 (167 steps against 47 on V-MMS-02, 225 to 239 against 64 to 71 on V-MMS-03) every pin holds in both forms: gas drift 1.4e-14 (temperature) and 4.3e-15 (enthalpy), uniform solid 1.3e-11 and 7.0e-12; with converged stage solves alone V-MMS-02 holds (2.5e-15, 2.1e-15) and V-MMS-03's temperature form stays at 1.03e-10. Not resolved: either these two cases run at the tighter tolerances on the new core, keeping their pins as exactness checks, or the two QoIs are re-pinned for the new core at the integrator's tolerance (a pin at three times the observed gas drift would be fragile, since the drift follows the Newton path); both are the owner's call. Until then the two cases are routed but kept out of `NEW_CORE_CASES`, and tests 1 and 2 hold the property (planar free stream within 1e-13 in both forms at fixed steps; isothermal advection within 1e-9 K).
- **C23. Property-layer cuts the new core still reaches.** `volume_fractions` and `mixture_density` return zero below a total of 1e-20 (`src/Properties/mixture_rules.jl`), `_carman_kozeny_permeability` uses `max` and `min` (`src/Properties/effective_properties.jl`), and `ConvectivePressureBC` vents only from a permeable cell (kept from the old core in `outward_darcy_velocity`). All are finite with finite derivatives, so the trial-state contract holds; they are kinks the quality grep of §3.8 (new modules only) does not see. Smooth variants are added beside them if a gate traces a failure to one.

### 12.5 Implementation parameters chosen (provisional, gate-measured)

- Mass tolerance `reltol · ε_m · mˢᵢ(0)` with `ε_m = 1e-2`; enthalpy tolerance `ε_m · mˢᵢ(0) · c_p,ref · abstol`; face-intensity tolerance `reltol · I_scale`, `I_scale` the largest incident intensity on 101 samples of the window, at least 1 W/m².
- Positive-part scale `δᵢ` equals the cell's mass tolerance. Properties are evaluated at `ξ⁺` and the trial temperature; transport and kinetics read the raw concentrations.
- Smooth floors: volume `1e-8 Vᵢ(0)`, center-to-face distances `1e-8` of their reference values, heat-capacity sum `RHO_CP_FLOOR · Vᵢ(0)`, conductivity 1e-6 W/(m·K), porosity 1e-12 (the old core's `_PHI_FLOOR`), viscosity `MU_FLOOR`; face-mean guards 1e-20 (conductivity, gas transfer) and 1e-30 m²/(Pa·s) (mobility). Trial temperature: 1 K transitions inside `[T_PROPERTY_MIN, T_PROPERTY_MAX]`.
- `OneSidedLimiter` defaults: `join = 1` (in units of the threshold), `base_scale = 1e-8` kg/m³.
- Gas transport model: none without a gas component; Darcy–Fick when `has_darcy_flow`; Fick when some component's gas-transfer coefficient is nonzero at any of four probe temperatures or a mass condition is not impermeable; none otherwise.
- Caloric extension and inversion as before (guaranteed capacity from half the sampled minimum; `κ = Cᵢ (1 − C/(2Cᵢ))²`; convergence at `4 eps · max(|H|, C |T − T_REF|, C)`).
- Admissibility (`inadmissible`): a mass below minus its absolute tolerance; on a moving mesh crossed nodes (a non-positive volume) and a non-positive area. The porosity check with the M8 predicate (D8).
- Moving-mesh tolerances: node positions `reltol · ε_m · hₙ` (the thinnest initial cell beside the node), velocities that over the first step (`time_scale`, the solve's `dt`), the history pair the cell's mass tolerance, the area `reltol · A₀`, the lateral rate the smallest velocity tolerance over the initial thickness.
- Lateral rows: conversion clamp band `w = 1e-3` outside `[0, 1]`; `Σ M` floored at `1e-8 Σ M(0)`.
- Relative fluxes: van Albada `ε` per cell the square of a tolerance over the reference cell thickness (concentrations `δᵢ/Vᵢ(0)`, temperature the solve's `abstol` [K], per-solid-mass history values `reltol · ε_m`); the per-solid-mass values divide by the solid mass smoothly floored at `δᵢ`.

### 12.6 Measurements so far

- Inventory-form cost (gate item 5 input): on V-HC the enthalpy form takes about 2× the accepted steps and 10 to 20× the wall time of the temperature form (V-HC-01: 135 against 69 steps, 0.22 s against 0.013 s).
- New core against old on the M3 cases (same session, accepted steps; the new core's raw estimate against the old core's smoothed one, C12): V-RX-01 373 against 185, V-RX-02b 150 against 104, V-RX-03 315 against 143, V-RX-04a 58 against 32, V-RX-04b 141 against 122, V-CPL-04 24 against 29; V-RAD-01 73 against 61.
- Temperature-form energy ledger under fixed-step refinement on the char-layer column (10 cells, 10 s, defect against a ledger of 3.1e5 to 3.5e5 J): Fick 0.63, 0.14, 0.026, 0.0029 J at steps 0.1, 0.05, 0.025, 0.0125 s (orders 2.1, 2.5, 3.2); Darcy 0.019, 0.0029, 0.00029 J at 0.05, 0.025, 0.0125 s (orders 2.7, 3.3); the inventory form at 1e-9 J, round-off.
- The fixed-step Newton failure at a Darcy flow reversal (C21).
- Moving-mesh Jacobian (warm, function scope, the reacting slab of `test_mesh_motion.jl`; 25 and 100 cells): colors 12 fixed, 17 moving, 24 with the cross section; one colored Jacobian costs 11 to 13 residuals on the fixed mesh, 21 to 24 moving and 28 to 33 with the cross section (the dense rows add a residual and the six-sum assembly); the residual allocates nothing, the colored Jacobian 14 to 54 KiB fixed, 26 to 95 KiB moving, 29 to 108 KiB with the cross section, growing with the cell count (the M7 allocation budget, D3).
- New core against old on the M6 cases (same session): V-CPL-03's P1 solve 141 accepted steps, 4 rejected (old 85); total-mass error 1.6e-15 against the old core's 6.6e-7 (masses are the states); V-MMS-03's T error 0.0088 K against 0.071 K, gas 1.68e-3 against 1.78e-3.
- The gas free-stream drift of V-MMS-02 on the new core against the Newton tolerance: 1.7e-12 (default `κ`), 2.2e-9 (`κ = 1e-3`), 2.5e-15 (`κ = 1e-6`); without gas diffusion 1.9e-14; at reltol 1e-11 1.1e-14; at a fixed step of 0.05 s 9.5e-15 (C24).

### 12.7 Next

The owner's call on C24 closes M6 (V-MMS-02 and V-MMS-03 into `NEW_CORE_CASES`). Then M7 as in §7 with its remainder (D3, D17): the containment suite over the `test/perf/perf_helpers.jl` configurations, the colored Jacobian's allocation budget, the color counts against the Phase 0 formula, and leftover 1; then M8 (with D8, D9, C21, C22) and M9.
