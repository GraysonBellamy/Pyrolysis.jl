# Mesh framework: test port ledger

The old core's tests are deleted with the old core. This document records, before the deletion, every obligation (a
property a test asserts that matters) of every test file the deletion removes or edits, and its disposition. Its scope:
the old core's test directories `test/jacobian`, `test/discretization`, `test/solver`, `test/perf`, `test/adaptivity`,
`test/ale`, `test/residual` and `test/radiation`; `test/verification/test_analytical.jl`; `test/helpers.jl`; and the
kept shared-layer directories `test/boundary_conditions`, `test/materials`, `test/physics` and `test/geometry`.

Each obligation has one disposition:

- **covered**: a test of the new core (`test/newcore/`), a verification case that runs on the new core
  (`test/verification/cases/`), or a kept shared-layer test asserts it.
- **ported**: the property is now asserted by a new test, written for it and named in the Where column.
- **dropped**: the obligation tests a deleted feature, or a definition that has no caller once the old core is gone.
  The Where column names it and the decision that deletes it.
- **kept** (the shared-layer directories only): the assertion stays in its file, on the new core or a surviving
  shared-layer function.

For covered, ported and kept rows the Where column names the file and the testset (for an interpolated name, its
literal prefix and the configuration; `>` marks nesting) and what it asserts. One row stays **port**: its test is
written in the deletion commit (Section 4).

The deletion commit's message gives each deleted directory's counts from Section 1 and cites this document (C125).
Rows are numbered per file. Tests are cited by testset name, not by line.

## 1. Counts

| Directory | Files | Obligations | Covered | Ported | Port | Dropped | Kept |
|---|---|---|---|---|---|---|---|
| `test/jacobian` | 26 | 136 | 69 | 14 | 0 | 53 | 0 |
| `test/discretization` | 6 | 81 | 23 | 38 | 0 | 20 | 0 |
| `test/solver` | 14 | 75 | 32 | 23 | 1 | 19 | 0 |
| `test/perf` | 8 | 12 | 9 | 1 | 0 | 2 | 0 |
| `test/adaptivity` | 4 | 53 | 21 | 9 | 0 | 23 | 0 |
| `test/ale` | 4 | 15 | 7 | 2 | 0 | 6 | 0 |
| `test/residual` | 2 | 10 | 2 | 5 | 0 | 3 | 0 |
| `test/radiation` | 4 | 30 | 1 | 3 | 0 | 26 | 0 |
| `test/verification/test_analytical.jl` | 1 | 25 | 20 | 5 | 0 | 0 | 0 |
| `test/helpers.jl` | 1 | 2 | 0 | 0 | 0 | 2 | 0 |
| `test/boundary_conditions` | 3 | 41 | 5 | 1 | 0 | 4 | 31 |
| `test/materials` | 5 | 25 | 0 | 0 | 0 | 1 | 24 |
| `test/physics` | 3 | 18 | 5 | 0 | 0 | 2 | 11 |
| `test/geometry` | 2 | 14 | 13 | 0 | 0 | 0 | 1 |
| **Total** | **83** | **537** | **207** | **101** | **1** | **161** | **67** |

## 2. The ledger

### test/jacobian

#### test_chi_and_reduced.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | With a lateral law on a moving mesh the Jacobian's χ block is present, finite and non-zero. | covered | `test/newcore/jacobian/test_containment.jl` "the cross section's back-reaction: …": a lateral law with non-degenerate histories; colored equals dense per entry, and the area, history and θ_A entries are non-zero. |
| 2 | Without a lateral law a moving-mesh layout has no χ states. | ported | `test/newcore/test_mesh_motion.jl` "moving mesh: …" (all nine configurations, `use_ale = true`): `has_chi(layout) == (name in (:chi, :darcy))`, so a moving mesh without a lateral law has no history, area or lateral-rate states. |
| 3 | `ApproxSolve`'s reduced sparsity has entries. | dropped | The structured Jacobian (`ApproxSolve`). |
| 4 | The reduced Jacobian evaluates to finite values. | dropped | The structured Jacobian (`ApproxSolve`). |
| 5 | An ALE solve with `jacobian = Structured(strategy = ApproxSolve())` succeeds. | dropped | The structured Jacobian; the `jacobian` keyword is among the removed keywords (owner question 18). |

#### test_darcy_venting_jacobian.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The species rows equal dense AD per entry on a venting state with all four features: φ_face from intrinsic densities across a char front, solid-basis κ mixing, zero mobility beside an impermeable solid, overpressure against a top pressure condition. | covered | `test/newcore/test_gas_transport.jl` "gas transport: …" (`darcy_venting`, `darcy_combined`): containment and colored equal to dense on every row, at 6 and 8 cells, both energy forms, five states. |
| 2 | The venting term is load-bearing: the top gas row's diagonal is non-zero. | covered | `test/newcore/test_gas_transport.jl` "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out": overpressure in the top cell, outward gas flux at the top face, zero mobility in the impermeable cells. |

#### test_face_mean_guard.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Residual and Jacobian share one guard band: below the cut both return a zero face value and zero derivatives. | dropped | The structured Jacobian's hand derivative `_harmonic_mean_with_derivs`, and the cut face mean: deleted with the old core, `harmonic_mean_weighted` has no caller outside it. Under AD one function is differentiated. |
| 2 | The face mean is finite and non-negative, with a finite derivative on both sides of the guard, and continuous across it. | covered | `test/newcore/jacobian/test_containment.jl` "the smooth face mean is finite with a finite derivative through its guard": `guarded_harmonic_mean` and its ForwardDiff gradient are finite on a grid through the guard, the value is non-negative, and there is no jump at the guard. |
| 3 | In the regular regime the face value is the distance-weighted harmonic mean. | covered | Same testset: far above the guard it equals the distance-weighted harmonic mean to 1e-14 (its reference is the closed form of `test/discretization/test_gas_flux.jl` row 8). |
| 4 | The operator's derivatives equal the analytic harmonic-mean sensitivities. | dropped | The structured Jacobian (the hand derivative `_harmonic_mean_with_derivs`). |
| 5 | One side at zero gives a zero face value, so a non-conducting or impermeable side blocks the face. | ported | `test/newcore/jacobian/test_containment.jl` "the smooth face mean is finite with a finite derivative through its guard": `guarded(0.0, b, d_a, d_b) == 0` and `guarded(b, 0.0, d_a, d_b) == 0` for `b in (η, 0.45)`. |
| 6 | The guard branch returns a zero of the promoted type (Dual, Float32). | dropped | The cut face mean's branch: deleted with the old core, `harmonic_mean_weighted` has no caller outside it. The smooth form has no branch and takes Dual arguments (row 2). |

#### test_geometry_ad.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The Jacobian's node-position columns equal dense AD (ALE, Fickian). | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_no_radiation` and the three ALE Beer–Lambert decks): no dense non-zero outside the pattern and colored equal to dense per entry over every column, position columns included; 8 states, 6 and 8 cells. |
| 2 | `apply_geometry_jacobian!` with a `GeometryADCache` fills the geometry block as a separate pass. | dropped | The structured Jacobian (its geometry AD pass). |
| 3 | On the ALE Beer–Lambert Darcy deck the geometry block holds at least n_nodes non-zeros, so row 1 is not vacuous. | ported | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (the four moving configurations, 6 and 8 cells, both forms, state `mixed`): every node's position column has a non-zero in the cells' energy or mass rows, so containment covers every position column. At the uniform `hot` state the position columns are nearly empty, which is why the check runs at `mixed`. |
| 4 | A fixed mesh has no node block: the state length is n(1 + NC). | covered | `test/newcore/test_states.jl` "state layout: interleaved blocks and row kinds": `state_length == n(NC+1+2χ) + (n+1)(2·moving + intensities) + 2χ`, so a fixed mesh has no node rows. |

#### test_geometry_ad_sparse_scatter.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The dense and sparse scatters of the geometry block agree bit for bit on the structural entries. | dropped | The structured Jacobian (the scatter in `src/Jacobian`). |
| 2 | The pattern misses no node-position entry, including on the optically moderate Beer–Lambert decks whose energy rows read every node above, and the block equals dense. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_no_radiation`, `ale_beer_lambert_darcy`, `ale_beer_lambert_darcy_moderate`, `ale_beer_lambert_fickian_moderate`): `missing_entries` empty and colored equal to dense per entry over all columns. |

#### test_interface.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The backend interface's public surface exists, and `JacobianSpec(backend, def, ws)` is inferred. | dropped | The structured Jacobian (`JacobianBackend`, `JacobianSpec`). |
| 2 | A backend without `_compute_jacobian!` or `_build_jacobian_spec` fails with a clean message. | dropped | The structured Jacobian. |
| 3 | `FastSparsityMap` and `fast_sparse_add!` add in-pattern writes and drop out-of-pattern ones. | dropped | The structured Jacobian (`src/Jacobian`). |

#### test_jacobian_coupling_regressions.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The comparison is per entry with a row-magnitude floor, so a missing small block fails whatever the problem size. | covered | `test/newcore/jacobian/test_pattern.jl` "the comparison is per entry, with a row-magnitude floor": the same synthetic matrix and cases on `compare_entries`, the helper every containment test uses. |
| 2 | Node velocities depend on the state of every upstream cell (the z-row × (T, ξ) block), at 6 cells with hot Arrhenius. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (moving configurations at 6 and 8 cells, state `hot`): containment, and the top velocity row's entries in the last cell's energy and mass columns are non-zero. |
| 3 | The θ_A back-reaction under a lateral law, with the zero-χ̄ trap guarded. | covered | `test/newcore/jacobian/test_containment.jl` "the cross section's back-reaction: …": non-degenerate histories guarded, ∂v/∂θ_A and the area and lateral-rate entries non-zero, 6 and 8 cells, four states. |
| 4 | The boundary cell's energy row carries the boundary mass flux and the ∂T_s chain: Fick and Darcy, a diffusive top or both faces, a convective plus reradiative top. | covered | `test/newcore/jacobian/test_containment.jl` "boundary mass flux on the energy row: …": the same five cases at 8 cells, same ramp and gas state; containment, and the top energy row's gas and neighbour entries and the permeable bottom's entry are non-zero. |

#### test_kirchhoff_eps_eff_jacobian.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Under Beer–Lambert with Kirchhoff absorptivity, ∂a_eff/∂T_n (a Kirchhoff ε(T) component, or a T-dependent density) reaches the energy rows, on fixed and moving meshes. | ported | `test/newcore/test_radiation.jl`, on the material `absorptivity_of_temperature` (virgin solid: Kirchhoff `ε = (0.95, -1.5e-4)`, no absorptivity, a density falling with T; char: `ε = 0.85, absorptivity = 0.6`): containment in "Beer–Lambert: …" under both tops and in "Beer–Lambert on a moving mesh: …" (five states, `J[e, c] != 0` at `u0` under the Kirchhoff top); in "Beer–Lambert: exposed-cell columns against central differences" the entering-intensity row's entry in the exposed cell's energy column is non-zero under the Kirchhoff top and the column matches central differences. |
| 2 | Under an explicit absorptivity the entering intensity does not depend on T_n. | ported | `test/newcore/test_radiation.jl` "Beer–Lambert: exposed-cell columns against central differences" and "Beer–Lambert on a moving mesh: …": on `absorptivity_of_temperature` under the explicit-absorptivity top, `abs(J[e, c]) <= 1e-12 * maximum(abs, J[e, :])`. The same entry is non-zero under the Kirchhoff top, so the check fails if the explicit branch reads the cell's absorptivity. |
| 3 | The exposed cell's ξ columns (attenuation, ∂a_eff/∂ξ) match dense per entry for Kirchhoff, explicit and component absorptivity, on fixed and moving meshes. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: …" (Kirchhoff and explicit tops, component absorptivity, 6 and 8 cells, five states) and "Beer–Lambert: exposed-cell columns against central differences" (the mass column is non-zero under Kirchhoff and matches central differences); moving: `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`beer_lambert`). |
| 4 | The T_n column against central differences in a thin medium, and no T_n chain at all when ε and ρ are constant. | ported | `test/newcore/test_radiation.jl` "Beer–Lambert: exposed-cell columns against central differences": the exposed cell's energy and virgin-mass columns against central differences on both materials under both tops, the T-dependent one included; on `component_absorptivity`, whose absorptivity does not depend on T, `abs(J[e, c]) <= 1e-12 * maximum(abs, J[e, :])` under the Kirchhoff top. The chain runs through the exposed face's intensity state, so the check holds at any optical thickness, not only in a thin medium. |
| 5 | The Jacobian reads the current incident flux, not the one at construction, on fixed and moving meshes. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: …": the flux and the absorptivity are functions of t, evaluated at every call, and colored and dense agree per entry at t = 1; `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`beer_lambert`) on a moving mesh. |

#### test_operator_ale_mesh_velocity.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The operator's `apply!` reproduces the residual's `w_node`, `θ`, `θ_A` bit for bit. | dropped | The structured Jacobian (its operator layer). |
| 2 | At a hot reacting ALE state the strain rate θ is non-zero, so the mesh velocity reads the cells' state. | covered | `test/newcore/test_mesh_motion.jl` "balance rows: a uniform strain rate gives θ x, θ r/2, θ r/3" (θ > 0 and w = θ x/d on all three coordinates); `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_no_radiation`, `hot`: the top velocity row's energy and mass entries are non-zero). |
| 3 | Without a lateral shrinkage law the lateral strain rate is zero. | covered | `test/newcore/test_mesh_motion.jl` "balance rows: a uniform strain rate gives θ x, θ r/2, θ r/3" (a law-free slab has axial w = θ x/d with no θ_A share) and "solve: the swelling column's thickness follows its conversion" (the cross-section history stays 1). |

#### test_operator_bc_convection.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `apply!` equals the residual's boundary-cell `S_conv/ρcp` slice, Fickian and Darcy. | dropped | `S_conv` and the structured Jacobian's operator layer. |
| 2 | Gas leaving through a boundary face carries its sensible enthalpy out of the boundary cell's energy row. | covered | At uniform temperature: `test/newcore/test_gas_transport.jl` "gas transport: isothermal transport leaves the temperature at rest (…)" (a Darcy vent at the top, both forms: the temperature rows are at rest to 1e-12, which holds only if the boundary face carries F_mass·h). |
| 3 | The same enthalpy is carried at the solved surface temperature T_s, not the cell's, when the two differ. | ported | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …" (Fick and Darcy, both forms, uniform and stretched meshes): every boundary face's energy flux less Σ_j F_mass[j, f]·h_j(T_face[f]) equals the conducted flux to 1e-12. At the top face T_face is more than 1 K from the cell temperature and the carried enthalpy is above 1e-3 of the conducted flux, so enthalpy taken at the cell temperature would miss by 1e-3 to 3e-2 of the conducted flux. |
| 4 | `contribute_to_jacobian!` equals dense AD on the boundary cell and its neighbour. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 5 | The boundary energy row's T_s chain reaches the neighbour cell and the boundary cell's gas (the reference's non-vacuity guard). | dropped | The guard of row 4's hand-written operator check (the structured Jacobian), and the `S_conv` boundary gradient that produced the neighbour reach. In the flux form T_s reads only the boundary cell, and AD makes the chain exact. |
| 6 | Support completeness: the boundary rows reach nothing beyond the boundary cell and its neighbour. | covered | `test/newcore/jacobian/test_containment.jl` "boundary mass flux on the energy row: …" (no dense non-zero outside the pattern, same fixture); `test/newcore/test_gas_transport.jl` "gas transport: …" (`fick_diffusive`, `fick_film`, `darcy_venting`, `darcy_combined`: 3(NC+1) colors under a boundary mass flux). |

#### test_operator_bc_thermal_mass.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The thermal and mass operators' `apply!` equals the residual's boundary-condition slice. | dropped | The structured Jacobian (its operator layer). |
| 2 | The thermal condition (prescribed flux plus convection, adiabatic back) enters the boundary cells' energy rows as the surface-balance face flux. | covered | `test/verification/cases/V-HC-01.jl` QoIs "L∞ T-profile error, t=… s" and "surface temperature, t=… s" (`HeatFluxBC + ConvectiveBC` top, adiabatic bottom); `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (T_s = T_cell + q d/k, conducted flux = q); V-HC-02 (convective, Biot sweep). |
| 3 | The `ConvectiveMassBC` film flux enters only the boundary cell's gas rows, with the right sign. | covered | Sign and solid rows: `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition" (the top gas row is positive, the solid rows and solid fluxes are zero). |
| 4 | The film flux's value, h_m(ξ_inf − ξ), inward and outward. | ported | `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition": at `u0` the top face's `F_mass[2, top] ≈ -0.05·0.5` and the top gas row ≈ `A·0.05·0.5`; with a uniform gas mass of `1.0·V_i` the flux is `0.05·(1.0 - 0.5)`, the top row `-A·0.05·0.5` and every other gas row exactly 0; rtol 1e-12. |
| 5 | `contribute_to_jacobian!` equals dense AD on the cell-local support. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 6 | Support completeness: the boundary rows reach only the boundary cell's (T, ξ). | covered | `test/newcore/test_gas_transport.jl` "gas transport: …" (`fick_film`: `HeatFluxBC + ConvectiveMassBC` top, convective plus `MassFluxBC` bottom, five states): no missing entry, colored equal to dense; `test/newcore/jacobian/test_containment.jl` "boundary mass flux on the energy row: …". |
| 7 | ∂T_s/∂ξ through k_eff(ξ) reaches the exposed cell's energy row exactly. | covered | `test/newcore/jacobian/test_containment.jl` "surface absorption with ε(T) and an absorptivity: …": the exposed energy row against central differences in its energy and both mass columns, with composition- and T-dependent k under a convective plus reradiative balance; `test/newcore/test_radiation.jl` "Beer–Lambert: exposed-cell columns against central differences". |

#### test_operator_conduction.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `apply!` equals −div(q_cond)/ρcp computed from the residual's own kernels. | dropped | The structured Jacobian (its operator layer). |
| 2 | Conduction on a temperature ramp with composition-dependent k: Fourier face fluxes with the distance-weighted harmonic-mean conductivity, summed into the cell rows. | covered | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …": each interior face's conducted part is the Fourier flux with the guarded harmonic mean, and each energy row is its signed face sum plus sources, both forms; `test/verification/cases/V-HC-04.jl` QoI "steady …-face flux vs network"; V-HC-05 (k(T)). |
| 3 | A uniform temperature conducts nothing. | covered | `test/newcore/test_gas_transport.jl` "gas transport: isothermal transport leaves the temperature at rest (…)": uniform 650 K, adiabatic, temperature rows within 1e-12·T. |
| 4 | `contribute_to_jacobian!` equals dense AD on the T row × three-cell stencil. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 5 | Support completeness: conduction rows reach only T and ξ of cells i−1 to i+1. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`: no missing entry); `test/newcore/test_residual_jacobian.jl` "residual and colored Jacobian: …" (fixed-mesh rows reach three cells: 3(NC+1) colors). |

#### test_operator_conduction_analytical.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The hand-derived conduction rows equal the operator's AD path to 1e-12 on 4 configurations × 5 states. | dropped | The structured Jacobian (hand-written operator Jacobians). |
| 2 | The conduction couplings (T and composition of the three-cell stencil, fixed and moving meshes) are exact on the canonical configurations at hot, ramp and uniform states. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …": the same four configurations and three more, eight states, colored equal to dense per entry to 1e-12. |

#### test_operator_darcy_transport.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `apply!` equals the residual's interior Darcy–Fick slice. | dropped | The structured Jacobian (its operator layer). |
| 2 | The interior gas energy term is `S_conv`. | dropped | `S_conv` (`Physics/convection.jl`). |
| 3 | Each interior face carries the sensible enthalpy of its mass flux at the face temperature, under Darcy–Fick. | covered | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …" (Darcy): each interior face's energy flux minus Σ F_mass·h(T_f) equals the Fourier flux to 1e-12, and the carried term is not negligible; `test/newcore/test_gas_transport.jl` "gas transport: isothermal transport leaves the temperature at rest (…)". |
| 4 | Pure Darcy interior transport: face mobility, pressure gradient, upwinding, venting outflow. | covered | `test/verification/cases/V-GT-02.jl` (λ = 0) QoIs "P² profile error / overpressure² scale", "back-face pressure", "exit mass flux vs ġL", and its first-order ConvergenceSpec. |
| 5 | Under Darcy–Fick an interior face carries the Fick flux plus the upwind Darcy flux, both active. | ported | `test/newcore/test_gas_transport.jl` "gas transport: an interior Darcy–Fick face carries the Fick flux plus the upwind advection of ξ/φ, and no advection across the char front" (`_gas_slab(8; darcy = true, …)` at `u0`): every interior face and gas has `F_mass[j, f] ≈ fick + v·ξ_up/φ_f` to 1e-12, with v, the Fick and the advective parts non-zero on the faces between char-bearing cells, and the solid fluxes zero. The reference face means are the closed-form distance-weighted harmonic means (zero when either side is zero), not the kernel's guarded mean. |
| 6 | `contribute_to_jacobian!` equals dense AD on the three-cell support. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 7 | Support completeness of the Darcy–Fick rows. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian_darcy`); `test/newcore/test_gas_transport.jl` "gas transport: …" (`darcy_venting`, `darcy_combined`). |

#### test_operator_darcy_transport_analytical.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The analytical Darcy rows equal the AD path to 1e-12 on 3 configurations × 5 states. | dropped | The structured Jacobian (hand-written operator Jacobians). |
| 2 | The hand-coded mobility chain dμ/dT follows the material's own viscosity curve, not a built-in Sutherland derivative. | dropped | The structured Jacobian (the hand-coded chain). AD differentiates whatever the residual evaluates, and the new core reads `material.gas_viscosity`. |
| 3 | The Darcy couplings (pressure, mobility, porosity, upwind) are exact on the Eulerian and the ALE Beer–Lambert Darcy configurations. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian_darcy`, `ale_beer_lambert_darcy`, `ale_beer_lambert_darcy_moderate`) and "chain-rule regime …" (`darcy_viscosity_hot_gas`); `test/newcore/test_gas_transport.jl` "gas transport: …" (`darcy_venting`, `darcy_combined`). |

#### test_operator_gas_transport.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `apply!` equals the residual's interior Fickian slice. | dropped | The structured Jacobian (its operator layer). |
| 2 | The interior Fick flux value, and its divergence into the gas rows. | covered | `test/newcore/test_gas_transport.jl` "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, and gas leaving through :bottom leaves the bottom cell": the interior `F_mass` equals `Physics.gas_flux` to 1e-14, and the bottom gas row changes by −A·F (rtol 1e-10); `test/verification/cases/V-GT-01.jl` QoIs "L∞ ξ-profile error / peak, t=… s", "gas-mass conservation (sealed)", "variance growth, t=… s". |
| 3 | The interior gas energy term is `S_conv`. | dropped | `S_conv` (`Physics/convection.jl`). |
| 4 | Each interior face carries the enthalpy of its Fickian mass flux at the face temperature. | covered | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …" (Fick). |
| 5 | `contribute_to_jacobian!` equals dense AD on the three-cell support. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 6 | Support completeness, with the solid rows receiving no transport. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`); `test/newcore/test_gas_transport.jl` "gas transport: …" and "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out" (the solid fluxes are all zero). |

#### test_operator_gas_transport_analytical.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The analytical Fickian rows equal the AD path to 1e-12 on 3 configurations × 5 states. | dropped | The structured Jacobian (hand-written operator Jacobians). |
| 2 | The Fickian gas-transport couplings are exact on fixed and moving PMMA and under surface absorption. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`, `eulerian_surface_absorption`, `ale_no_radiation`, eight states each); `test/newcore/test_gas_transport.jl` "gas transport: …" (`fick_diffusive`, `fick_film`). |

#### test_operator_radiation.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | NO_RADIATION and SURFACE_ABSORPTION deposit nothing in depth. | covered | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate" (under SURFACE_ABSORPTION no volumetric absorption, the surface absorbs q·ε); `test/newcore/test_boundary_conditions.jl` "surface absorption is the surface balance of the incident flux" (the NO_RADIATION residual equals SURFACE_ABSORPTION's); `test/newcore/test_radiation.jl` "Beer–Lambert is refused off planar coordinates" (no intensity unknowns under SURFACE_ABSORPTION). |
| 2 | P1_QUASI_STEADY is a no-op. | dropped | P1 radiation, deleted with Experimental (owner question 3). |
| 3 | The Beer–Lambert `apply!` and the workspace's `q_rad` equal `compute_beer_lambert_source!` for both conventions. | dropped | The structured Jacobian's operator wrapper, and deleted with the old core: `compute_beer_lambert_source!` has no caller outside it. |
| 4 | The beam entering the exposed cell is the surface absorptivity times the incident flux: the cell's ε_eff under Kirchhoff, the condition's absorptivity when explicit. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity" (I_in ≈ a·q, and the exposed cell's intensity equals it); `test/materials/test_emissivity_absorptivity.jl` "a_eff: Kirchhoff invariant and per-component override" (a_eff is ε_eff bit for bit without a component absorptivity). |
| 5 | Each cell absorbs I_up(1 − e^(−τ_i)) with τ = αΔz, so the Kirchhoff profile is the bare profile scaled by the surface absorptivity. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": the per-cell energy source, the transmitted face, and absorbed plus transmitted equal to the entering beam. |
| 6 | Non-Beer–Lambert modes zero the workspace's `q_rad` and `I_in`. | dropped | The structured Jacobian (the operator workspace `build_beer_lambert_workspace`). |

#### test_operator_reaction_source.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `apply!` equals the per-cell (Q_rxn/ρcp, species sources) slice under `TanhLimiter`. | dropped | The structured Jacobian's operator layer, and deleted with the old core: `TanhLimiter` has no caller outside it. |
| 2 | At a hot, gas-bearing state the reaction kernel's energy and species sources are the material's rates, heats and yields, per cell. | covered | `test/newcore/test_kinetics.jl` "reaction kernel: sources equal the material's rates" (S_mass = ω·V and S_energy = −Σ h·r·V, two reactions with yields, under the material's `DepletionLimiter`); the division by the heat capacity in the temperature form: `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …". |
| 3 | `contribute_to_jacobian!` equals dense AD on the cell-local (1+NC)² block. | dropped | The structured Jacobian (a hand-written operator Jacobian). |
| 4 | Support completeness: the reaction couplings are cell-local and reach T and every ξ. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`, states `hot` and `gas_heavy`: no missing entry) and "chain-rule regime …" (the window regimes: the mass row's temperature entry is non-zero). |
| 5 | The cold initial state (300 K, no gas) gives a finite source with no blow-up. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`, state `u0`, the same PMMA at 300 K without gas: residual and Jacobian finite) and "chain-rule regime …" (`cold_below_window`). |

#### test_park_jacobian.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | On a masked state the exact Jacobian equals dense AD (fixed and ALE Fickian, ALE Beer–Lambert). | covered | `test/newcore/test_deactivation.jl` "masked rows are identity rows and the live block is the compacted problem's: …" (`fixed`, `moving`, `beer_lambert`; both forms; set-up and perturbed states): the masked Jacobian equals dense per entry. After a real deactivation: "a transfer conserves the inventories and the receiver absorbs the volume: …" and `test/newcore/jacobian/test_containment.jl` "after deactivate! with in-depth radiation: …". |
| 2 | Dead-cell rows are identically zero. | covered | Same testset: dead differential rows (masses, energy, position) are zero, dead algebraic rows (velocity, intensity) are identity rows. |
| 3 | Ghost node rows copy the surface node row. | dropped | Parking (its ghost nodes), replaced by deactivation, where a dead node's position row is a zero differential row. |
| 4 | The live block equals the compacted problem's Jacobian. | covered | Same testset: the live block equals a 5-cell problem of the live thickness to rtol 1e-11, and the live-to-dead block is zero. |
| 5 | Rows 1, 2 and 4 under Darcy–Fick transport (fixed mesh with two dead cells; ALE with Beer–Lambert and a convective top). | ported | `test/newcore/test_deactivation.jl` "masked rows are identity rows and the live block is the compacted problem's: …" (`darcy`: fixed mesh, a permeable solid with `ρ_intrinsic`, a `DirichletPressureBC` top; `darcy_surface_balance`: moving mesh, BEER_LAMBERT, a radiative, convective and reradiative top; gas in every live cell; one and two dead cells; both forms; set-up and perturbed states): the set-up uses `DarcyFickTransport`, the live block equals the compacted problem's Jacobian, the live-to-dead block is zero, dead differential rows are zero and dead algebraic rows identity rows, and colored equals dense. |
| 6 | Rows 1, 2 and 4 under a lateral shrinkage law (ALE, two dead cells, a mass condition on top). | ported | Same testset, `lateral` (moving mesh, `χ -> 1 - 0.2χ`, a `DiffusiveMassBC(0.01)` top, histories (0.3 + 0.05i)·M and zero in the dead cells): the compact rows mapped through `cell_rows`, `node_rows`, `area_index` and `lateral_rate_index`, with the same assertions. |
| 7 | Parking an interior cell (cell 5 alone; cells 8, 4, 1) keeps the Jacobian exact. | dropped | Parking (interior parking). Deactivation takes only the exposed run; a sub-threshold cell under a live one ends the solve with outcome `:interior`. |
| 8 | `ApproxSolve`'s masked live block equals the compacted problem's, and its ghost rows are zero. | dropped | The structured Jacobian (`ApproxSolve`). |

#### test_regressions.jl

Every covered row relies on the containment check of `test/newcore/jacobian/test_containment.jl`: a finite residual and
dense Jacobian, no dense non-zero outside the pattern, and colored equal to dense per entry at rtol 1e-12, in both energy
forms, at 6 and 8 cells unless the regime says otherwise.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Baseline Eulerian PMMA at its default state. | covered | "perf configuration …" (`baseline_eulerian`, state `u0`). |
| 2 | Reaction window, lower edge (1 K above T_min): the gate's sech² slope. | covered | "chain-rule regime …" (`window_lower_edge`), which also asserts a non-zero temperature entry in the mass row. |
| 3 | Reaction window, upper edge (1 K below T_max). | covered | "chain-rule regime …" (`window_upper_edge`), with the same assertion. |
| 4 | A narrow gate (`gate_width = 0.1`): no Jacobian path hard-codes the 1 K ramp. | dropped | The structured Jacobian's hand paths. Under AD the width reaches the Jacobian only through `reaction_rate`, whose narrow-gate value and slope `test/physics/test_kinetics.jl` "temperature gate width" keeps. |
| 5 | Sutherland μ(T) in Darcy's ∂v_gas/∂T (ALE Beer–Lambert Darcy, 1000 K, ξ = (400, 5)). | covered | "perf configuration …" (`ale_beer_lambert_darcy`, state `gas_heavy`, the old state exactly) and "chain-rule regime …" (`darcy_viscosity_hot_gas`); μ(T) is the material's default `gas_viscosity`. |
| 6 | A T-dependent gas c_p, with nothing outside the support. | covered | "chain-rule regime …" (`gas_heat_capacity_of_temperature`): the same `LinearProperty(1100, 0.8)` and gas-rich state; no missing entry is the outside-support check. |
| 7 | Cell-varying reaction rates (the per-cell dS/dξ slice) on an Eulerian temperature ramp. | covered | "perf configuration …" (`baseline_eulerian`, state `mixed`: T from 450 to 650 K with a composition ramp). |
| 8 | `apply_limiter` accepts Duals (AD through the ALE residual at a T ramp). | dropped | The old ALE temperature-advection kernel and its limiters. |
| 9 | Eulerian PMMA hot (1000 K, (1100, 0)). | covered | "perf configuration …" (`baseline_eulerian`, state `hot`). |
| 10 | Eulerian PMMA hot and gas-heavy (1200 K, (400, 8)). | covered | "perf configuration …" (`baseline_eulerian`, state `gas_heavy`). |
| 11 | Eulerian PMMA partially decomposed (800 K, (700, 2)). | covered | "perf configuration …" (`baseline_eulerian`, state `partial`). |
| 12 | Eulerian surface absorption, hot with gas (1000 K, (700, 5)). | covered | "perf configuration …" (`eulerian_surface_absorption`, states `hot`, `gas_heavy`, `partial`): the same coupling, the surface balance with gas present, though not the exact state. |
| 13 | ALE without radiation, hot. | covered | "perf configuration …" (`ale_no_radiation`, state `hot`, with the non-zero velocity-row entries). |
| 14 | ALE without radiation, hot and gas-heavy. | covered | "perf configuration …" (`ale_no_radiation`, state `gas_heavy`). |
| 15 | ALE Beer–Lambert Darcy, hot (1000 K, (400, 0)). | covered | "perf configuration …" (`ale_beer_lambert_darcy`, state `hot`). |
| 16 | ALE Beer–Lambert Darcy, hot with gas (1000 K, (400, 5)). | covered | "perf configuration …" (`ale_beer_lambert_darcy`, state `gas_heavy`). |
| 17 | ALE Beer–Lambert Darcy, partially decomposed (800 K, (200, 2)). | covered | "perf configuration …" (`ale_beer_lambert_darcy`, state `partial`). |
| 18 | A single-cell mesh, hot. | covered | "chain-rule regime …" (`single_cell`, and `single_moving_cell`). |
| 19 | 200 cells, hot with gas (1000 K, (1100, 0.5)). | covered | "chain-rule regime …" (`two_hundred_cells`, the same state). |
| 20 | A zero-rate cold state (T_min = 800 K, T = 400 K). | covered | "chain-rule regime …" (`cold_below_window`: the same reaction and state, at 8 cells). |
| 21 | Three components, the κ-mixing chain rule (two permeable solids; ALE with Beer–Lambert in the old deck). | covered | "chain-rule regime …" (`permeability_mixing`: κ 1e-14 and 1e-12, fixed mesh); three components on moving meshes: `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`darcy`, `beer_lambert`). |

#### test_sparsity.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `build_sparsity` is the source of truth: both profiles' sizes and Bool type, the specs' patterns equal it, a geometry block if and only if ALE. | dropped | The structured Jacobian (`build_sparsity`, its profiles, `JacobianSpec`). The new pattern is checked by containment. |
| 2 | On the optically moderate ALE Beer–Lambert decks the pattern holds every far-node coupling. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_beer_lambert_darcy_moderate`, `ale_beer_lambert_fickian_moderate`): no dense non-zero outside the pattern in eight states. |
| 3 | The old reach: energy rows to the top node, species rows to i + 3, energy rows to node 4 without Beer–Lambert. | dropped | The structured Jacobian's operator pattern. The beam is carried by intensity states, and the reach is bounded by the color count and bandwidth (`test/newcore/jacobian/test_pattern.jl` "colors: …", "radial meshes: the planar pattern, bandwidth and colors"). |
| 4 | Threaded pattern builds share one cache safely. | dropped | The structured Jacobian (`_SPARSITY_CACHE`). The new core builds the pattern per `setup` and keeps no global cache. |
| 5 | The `sparsity_profile` trait. | dropped | The structured Jacobian. |
| 6 | Every row stores its diagonal, so W has J's pattern. | covered | `test/newcore/jacobian/test_pattern.jl` "colors: …" (fixed, moving, Beer–Lambert, cross section, reconstruction and equidistribution cases; NC = 2, 3); also `test/newcore/test_residual_jacobian.jl` "residual and colored Jacobian: …", `test/newcore/test_radiation.jl` "Beer–Lambert: …", `test/newcore/test_mesh_motion.jl` "moving mesh: …". |
| 7 | The reduced profile's structural-zero diagonals stay zero after assembly. | dropped | The structured Jacobian (`ApproxSolve`, `:ale_reduced`). |
| 8 | `materializes_geometry_block` (the default-factorization rule) agrees with the pattern. | dropped | The structured Jacobian (the rule chose UMFPACK for the ALE DirectSolve geometry block). The default is `StructuredKLU`: `test/newcore/test_policies.jl` "one KLU analysis per growth of the union, …". |
| 9 | The reduced ALE pattern stays within a one-cell stencil. | dropped | The structured Jacobian (`:ale_reduced`). The new bandwidth is constant in n and equal planar and radial (`test/newcore/jacobian/test_pattern.jl` "radial meshes: the planar pattern, bandwidth and colors"). |
| 10 | `:ale_reduced` gating of specific entries (χ columns, χ × z, the Beer–Lambert suffix, the Kirchhoff T_n slots of `:full`). | dropped | The structured Jacobian (the profile). The Kirchhoff T_n coupling it located is `test_kirchhoff_eps_eff_jacobian.jl` row 1. |

#### test_structured_backend.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `is_approximate`: DirectSolve exact (the default), ApproxSolve approximate. | dropped | The structured Jacobian. |
| 2 | The default backend agrees strictly with dense AD on six canonical configuration and off-default state pairs. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`baseline_eulerian`, `eulerian_surface_absorption`, `ale_no_radiation`, `ale_beer_lambert_darcy`; states `hot`, `gas_heavy`), per entry at rtol 1e-12. The old surface-absorption state (1000 K, (700, 5)) is not among the states; the same coupling with gas present is. |
| 3 | Strict agreement at a ramped ALE state (where the old limiter was active). | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_no_radiation`, state `mixed`) and `test/newcore/test_mesh_motion.jl` "moving mesh: …" (state `u0`, a temperature ramp). The limiter it exercised is deleted. |
| 4 | ApproxSolve routes the verifier to the masked comparison. | dropped | The structured Jacobian (`ApproxSolve`). |
| 5 | The ApproxSolve masked round trip on an ALE lateral-law Darcy surface-absorption deck. | dropped | The structured Jacobian (`ApproxSolve`). The exact Jacobian of a lateral-law Darcy deck is "the cross section's back-reaction: …". |
| 6 | The ApproxSolve masked round trip on ALE Beer–Lambert Darcy. | dropped | The structured Jacobian (`ApproxSolve`). |

#### test_structured_couplings.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `densify!` of the node prefix-sum, scaled suffix-sum and rank-1 couplings, and `compose_into!`, match dense references. | dropped | The structured Jacobian (its couplings). The couplings they encoded are contained per entry (`test_jacobian_coupling_regressions.jl` rows 2 and 3; the Beer–Lambert decks of "perf configuration …"). |

#### test_structured_regressions.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Ramp-ALE strict agreement: nothing outside the support, in-support deviation within 1e-7 relative. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_no_radiation`, state `mixed`): no missing entry and per-entry equality at rtol 1e-12; the velocity-row couplings at state `hot`. The dilation rows are a deleted term. |
| 2 | A local temperature minimum gives a finite Jacobian that agrees with dense AD through the van Leer limiter at r = 0. | dropped | The old ALE temperature-advection kernel and its limiters. The surviving van Albada slope has a finite gradient at zero slopes (`test/newcore/test_relative_flux.jl` "van Albada slope") and is contained at sign-changing velocities ("relative fluxes: containment, …"). |

#### test_surface_radiation_eps_T_jacobian.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Under surface absorption, an ε(T) component's emission chain and a component absorptivity enter the exposed energy row exactly (fixed mesh). | covered | `test/newcore/jacobian/test_containment.jl` "surface absorption with ε(T) and an absorptivity: …" (`component_absorptivity`, fixed mesh: virgin ε = (0.8, 5e-5) with absorptivity 0.7, Kirchhoff char, radiative plus convective top, both forms): containment and per-entry equality. |
| 2 | The absorption slope q_inc·da/dT of a Kirchhoff ε(T) component, also on a moving mesh. | ported | `test/newcore/jacobian/test_containment.jl` "surface absorption with ε(T) and an absorptivity: …" (`absorptivity_of_temperature`: virgin Kirchhoff `ε = (0.95, -1.5e-4)` without an absorptivity, char `ε = 0.85, absorptivity = 0.6`; fixed and moving meshes, both forms): containment, and the exposed energy row in its energy and two mass columns against central differences, rtol 1e-5. |
| 3 | The surface-temperature column of the exposed energy row against central differences. | covered | Same testset: the exposed row in its energy and two mass columns against central differences of `residual!`, rtol 1e-5, both forms. |
| 4 | The same check on the conduction neighbour's row n − 1, and on the old material, fixed and ALE. | ported | Same testset: the conduction neighbour's row `energy_index(layout, n - 1)` in the exposed cell's energy column against central differences, rtol 1e-5, on both materials, fixed and moving, both forms. |
| 5 | Load-bearing: ε frozen at its value reproduces the residual but moves J[n, n] off the central difference by more than 20× the tolerance. | ported | `test/newcore/jacobian/test_containment.jl` "surface absorption: the ε(T) slope is load-bearing in the exposed energy row: …" (fixed and moving meshes, both forms): with the virgin ε frozen at its exposed-cell value the residual equals the live material's (rtol 1e-12), the live diagonal is within 1e-5·abs(CD) of the central difference and the frozen one misses it by more than 20 times that. The state (400 + 40i K) is cooler than the main testset's, where the reaction heat swamps the slope (C129). |

### test/discretization

#### test_conservation.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `MassConservationTracker` and `compute_domain_mass` start from the inventories ρ·V·A per component. | dropped | Diagnostics (owner question 4): the trackers are replaced by `ConservationRecorder` and `conservation_audit`. |
| 2 | The tracker accumulates the boundary efflux J·A·dt over substeps. | dropped | Diagnostics (owner question 4). The recorder's integrated outflow is checked against the trapezoid of `ledger_rates` in `test/newcore/test_conservation.jl` "the recorder on the gasifying column: …". |
| 3 | A mass-conserving evolution passes the mass check. | covered | Total mass on a solve: `test/newcore/test_conservation.jl` "the recorder on the gasifying column: …" (a venting deck's mass defect within 1e-8 of the initial mass) and "the recorder across deactivation events: …" (relative defect within 1e-7). The per-component closure belongs to the tracker (Diagnostics). |
| 4 | The mass books close at round-off between the residual's mass rows and the boundary outflow (the old check's 1e-10). | ported | `test/newcore/test_energy_ledger.jl` "energy ledger over a solve: …" (Fick and Darcy): the solve carries the integrated boundary mass outflow of `ledger_rates` as a second augmented unknown; the total mass plus that outflow stays at the initial inventory M0 within 1e-13·M0 in the inventory form and at every temperature-form step, with the outflow above 1e-4·M0. The bound is relative to the inventory: 1e-13 of the outflow lies below the round-off of the mass sum (C129). |
| 5 | The mass and energy reports print their headings and PASSED or FAILED. | dropped | Diagnostics (owner question 4). |
| 6 | `EnergyConservationTracker` and `compute_sensible_energy` start from Σ ξ_j c_j T·V·A. | dropped | Diagnostics (owner question 4). |
| 7 | Pore gas stores sensible heat: the capacity sums every component, gas included. | covered | `test/newcore/test_states.jl` "the heat-capacity sum skips the absent components": `heat_capacity_sum` is Σ m_j c_j over all three components, the gas among them. |
| 8 | The energy tracker accumulates boundary conduction, reaction heat and the `S_conv` sink. | dropped | Diagnostics (owner question 4); `S_conv` is deleted. |
| 9 | A duplicate save time does not advance the diagnostics stride. | dropped | Diagnostics (owner question 4); `diagnostics_stride` is among the removed keywords (owner question 18). |
| 10 | An energy-conserving evolution passes the energy check (relative defect below 1e-6). | covered | `test/newcore/test_conservation.jl` "the recorder on the gasifying column: …" (energy relative defect within 1e-6 on a gasifying, venting deck); `test/newcore/test_energy_ledger.jl` "energy ledger over a solve: …" (round-off in the inventory form, the integrator's order in the temperature form). |
| 11 | Freshly built trackers report conservation at t = 0. | dropped | Diagnostics (owner question 4). |
| 12 | Under a time-varying cross section the boundary efflux, integrated at the current area, still closes the books. | ported | `test/newcore/test_conservation.jl` "the recorder on a venting column whose cross section shrinks: …" (both forms; `use_ale = true`, a lateral law `1 - 0.3χ`, a `DirichletPressureBC` top): `:Success`, `rec.steps == sol.statistics.accepted`, the area falls by more than 1 %, `audit.mass.inflow < 0`, the mass relative defect within 1e-7 and the energy one within 1e-6, and the inflow times the area's relative fall above 1e3 times the mass defect, so an efflux integrated at the initial area would fail. |

#### test_darcy.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `ideal_gas_pressure`: one gas at φ = 1 gives ξRT/M; P is linear in T; a solid-only material and a zero gas concentration give exactly 0. | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "ideal-gas pressure: P = Σ ξ R T/(M φ) over the gases": one gas at φ = 1 gives ξRT/M; a solid-only material and a zero gas give exactly 0; P(600 K)/P(300 K) = 2. |
| 2 | The pressure carries the porosity divisor P = N/φ. | covered | `test/physics/test_background_air.jl` "EOS round-trip at ambient": `ideal_gas_pressure` at φ = 2/3 returns 101325 Pa for the inverted air concentration, to 1e-12. |
| 3 | `porosity` sums the skeletal volumes of two solids (φ = 0.25) and the pressure divides by it. | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "the pore fraction sums the skeletal volumes and the pressure divides by it": φ = 1 − `skeletal_volume_fraction` = 0.25 over two solids, and P = N/φ. The clamped `porosity` is asserted in "porosity is the pore fraction on physical states" (C127, Section 4). |
| 4 | ForwardDiff through the porosity and `ideal_gas_pressure` matches central differences in ξ_solid, ξ_gas and T. | ported | Same parent > "pressure: ForwardDiff through the pore fraction equals central differences": ∂/∂ξ_s, ∂/∂ξ_g and ∂/∂T to 1e-6, at φ = 0.5 so the solid column is load-bearing. |
| 5 | The pressure and its derivatives stay finite when the skeleton over-fills the cell with gas present (φ at its floor, P = N/φ_floor). | ported | `test/newcore/jacobian/test_containment.jl` "perf configuration …" and "radial …" (the wood configurations, state `overfilled = uniform(800.0, (550.0, 1.0))`, a skeleton 10 % over its volume with gas): contained; under `DarcyFickTransport` (`baseline_eulerian_darcy`, `ale_beer_lambert_darcy`, `ale_beer_lambert_darcy_moderate`) every `ds.φ` equals `POROSITY_FLOOR`, every `ds.P` is finite and the smallest exceeds 1e12. |
| 6 | `effective_permeability`: one solid returns its κ; PARALLEL on two solids lies between them; gases add nothing. | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "permeability: one solid, two solids in parallel, gases add nothing". |
| 7 | Solid-basis mixing: a uniform κ₀ with 30 % gas gives κ₀ under PARALLEL and SERIES; an impermeable phase makes SERIES exactly 0, PARALLEL (0.2/0.7)κ and WEIGHTED their mean; AD equals FD across WEIGHTED. | ported | Same parent > "permeability mixes on the solid basis; an impermeable phase zeroes the series bound": the values to 1e-12, and AD against finite differences across WEIGHTED to 1e-5. |
| 8 | `has_darcy_flow` is true for κ > 0 and false for κ = 0. | covered | `test/newcore/test_gas_transport.jl` "gas transport: …": the set-up transport is `DarcyFickTransport` exactly for the Darcy decks (chosen by `has_darcy_flow`), and the Fick decks have κ = 0 on every solid. |
| 9 | Sutherland: 1.716e-5 at 273.15 K and 1.846e-5 at 300 K, rising with T, within 1 % of air at 400 and 500 K. | covered | `test/materials/test_gas_viscosity_property.jl` "SutherlandProperty is Sutherland's law": the closed form at six temperatures to 1e-14, both anchors. |
| 10 | A custom Sutherland triple differs from the default. | covered | `test/materials/test_gas_viscosity_property.jl` "end-to-end: doubling μ doubles the Darcy overpressure²" (`SutherlandProperty(2μ_ref, T_ref, S)`, on the new solve). |
| 11 | `mixture_gas_viscosity` equals `effective_gas_viscosity` equals the default Sutherland curve. | covered | `test/materials/test_gas_viscosity_property.jl` "Material defaults to air": the material's default `gas_viscosity` is the Sutherland curve (the two accessors are deleted with the old core). |
| 12 | Tripwire: μ is independent of ξ, which the analytical Darcy Jacobian assumes. | dropped | The structured Jacobian (its hand-written Darcy operator). The colored AD Jacobian differentiates any ξ dependence. |
| 13 | Sutherland at a negative trial temperature is finite and equals its value at the bracket. | covered | `test/materials/test_property_validation.jl` "argument clamp is the identity inside the bracket": `SutherlandProperty()(-5.0) == SutherlandProperty()(T_PROPERTY_MIN)`. |
| 14 | `chapman_enskog_diffusivity` at −43 K is finite and equals its value at 0 K. | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "Chapman–Enskog diffusivity" > "a nonphysical trial temperature evaluates at the floor": at −43 K the diffusivity is finite and equals its value at 0 K and at `T_GAS_PROPERTY_FLOOR`. |
| 15 | Chapman–Enskog: Neufeld Ω_D at T* = 1, 2, 5 and falling; D rising between T^1.5 and T^1.75; D ∝ 1/P; lighter gases faster; literature magnitudes; custom Lennard-Jones parameters; the air constants. | ported | Same parent's "Neufeld collision integral at T* = 1, 2, 5, falling", "rises with T between T^1.5 and T^1.75", "inverse in pressure, faster for lighter gases", "literature magnitudes in air at 1 atm" and "custom Lennard-Jones parameters and the air constants". |
| 16 | `gas_flux_darcy_fick` pointwise: v = 0 is `gas_flux`; λ = 0 is pure upwind; advection carries ξ/φ_face; the combined sum; the Péclet limits. | dropped | Deleted with the old core: `gas_flux_darcy_fick` has no caller outside it. |
| 17 | The face loop at zero Darcy velocity gives the diffusion-only flux. | ported | `test/newcore/test_gas_transport.jl` "gas transport: an interior Darcy–Fick face carries the Fick flux plus the upwind advection of ξ/φ, and no advection across the char front" (its last block, `_pressure_tie(8; outlet = false)`): v == 0 at every interior face, the Fick flux is non-zero, and `F_mass` equals it to rtol 1e-14 (the reference's closed-form mean differs from the kernel's guarded one by round-off). |
| 18 | A positive Darcy velocity adds upwind advection at every interior face. | ported | Same testset (`_gas_slab(8; darcy = true, …)`): every interior face has `F_mass[j, f] ≈ fick + v·ξ_up/φ_f` to 1e-12, with v = −m_f ΔP/d and ξ_up taken by the sign of v; the advective part is non-zero on every char–char face. |
| 19 | Solid components carry no Darcy–Fick flux. | covered | `test/newcore/test_gas_transport.jl` "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out": the solid rows of `F_mass` are all zero under Darcy venting. |
| 20 | `S_conv` values: −c_g J ∂T/∂z on a linear profile with the boundary blend, proportional to the gradient. | dropped | `S_conv` (`Physics/convection.jl`). The carried enthalpy is a face flux, checked term by term in `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …". |
| 21 | Gas flowing at a uniform temperature does not change the temperature. | covered | `test/newcore/test_gas_transport.jl` "gas transport: isothermal transport leaves the temperature at rest (…)": the temperature rate within 1e-12·T, the inventory rate Σ dm·h. |
| 22 | Pressure-condition constructors: Dirichlet default 101325 Pa and custom; Neumann default 0 and custom; the convective fields. | ported | `test/boundary_conditions/test_mass_bcs.jl` "Mass boundary conditions" > "pressure conditions: atmospheric and impermeable defaults, given values kept": `DirichletPressureBC().P == 101325` and a given P, `NeumannPressureBC().∂P∂n == 0` and a given one, and the fields of `ConvectivePressureBC(1e-6, 101325)`. |
| 23 | Overpressure vents through a Dirichlet top at outward velocity, and an impermeable bottom stays closed. | covered | `test/newcore/test_gas_transport.jl` "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out": overpressure in the top cell, outward gas flux at the top face, zero flux at the bottom face. |
| 24 | `NeumannPressureBC(0)` on both faces gives zero boundary velocity. | ported | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form": the Neumann velocity is −m·∂P∂n in both signs, and `NeumannPressureBC(0.0)` gives exactly 0, at three pressures. |
| 25 | `ConvectivePressureBC` gives v = h_P (P_int − P_ext). | covered | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form": `outward_darcy_velocity(film, 2e5, 2e-7, d) == 1e-6·(2e5 - 1e5)`. |
| 26 | Bruggeman: a single component, two equal components and a zero volume fraction return the property. | ported | `test/physics/test_gas_closures.jl` "Bruggeman effective medium" > "a single phase, two equal phases and an empty phase return the property". |
| 27 | Bruggeman lies strictly between series and parallel, equals the root 2.2612078, is symmetric under relabelling and satisfies its defining equation. | ported | Same parent > "strictly between the series and parallel bounds" (the root 2.2612078, three phases, a 10:1 contrast), "symmetric under relabelling, not under swapping the fractions" and "the root satisfies Σ v(k − k_eff)/(k + 2k_eff) = 0" (below 1e-8). |
| 28 | Bruggeman with an insulating phase: k_c(3p − 1)/2 above p = 1/3, 0 below, and the three-phase residual. | ported | Same parent > "an insulating phase: k_c(3p − 1)/2 above p = 1/3, zero below", with the three-phase residual. |
| 29 | `effective_conductivity` under BRUGGEMAN lies between the series and parallel bounds. | ported | Same parent > "the conductivity rule lies between the bounds; the gas-transfer rule refuses it": `effective_conductivity` under BRUGGEMAN lies between the series and parallel bounds. |
| 30 | BRUGGEMAN for the gas-transfer coefficient throws `ArgumentError`. | ported | Same testset: `effective_gas_transfer` under BRUGGEMAN throws `ArgumentError`. |
| 31 | `MixingSpec(CARMAN_KOZENY)` requires a positive finite `length_scale`; its fields; other rules unaffected. | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "Carman–Kozeny: κ = min(ℓ²φ³/(180(1 − φ)²), ℓ²)" > "the closure needs a positive finite length scale". |
| 32 | Carman–Kozeny κ = min(ℓ²φ³/(180(1−φ)²), ℓ²) across four porosities, monotone in φ, with a finite negative AD derivative at 800 kg/m³. | ported | Same parent > "closed form across porosities, increasing with porosity" (the 4-argument `effective_permeability(ξ, T, mat, φ)` at the state's pore fraction against the closed form at a hand-written φ = 1 − ρ_s/1500, to 1e-12) and "below the cap more solid lowers κ, with a finite AD slope". |
| 33 | The Carman–Kozeny free-flow cap: κ == ℓ² at φ ≈ 0.933 and φ = 1, with an AD derivative exactly 0 on the cap. | ported | Same parent > "the free-flow cap binds at high porosity with zero slope". |
| 34 | `_permeability_from_v` rejects CARMAN_KOZENY. | ported | Same parent > "the volume-fraction entry point refuses the closure": `_permeability_from_v` throws for CARMAN_KOZENY. |

#### test_darcy_boundary_venting.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The surface cell's venting rate grows with κ (more than 10× per two decades, 5× per half decade), unlike a film law alone. | ported | `test/newcore/test_gas_transport.jl` "gas transport: a pressure condition vents its cell at v ξ/φ beside the film flux, faster through a more permeable solid": the top outflow over κ in (1e-14, 1e-12, 1e-11) grows by more than 10× and 5× (about 100× and 10×). |
| 2 | Uniform state: the surface gas tendency is exactly −(v_b ξ/φ + λhξ/(λ + hd))/Δz with v_b = (κ/μ)(P − P_bc)/d > 0, and the advective part is 2v_bξ at φ = 0.5. | ported | Same testset, on a uniform char column (φ = 0.5, 600 K) under `DiffusiveMassBC(0.05) + DirichletPressureBC(101325)`: `ds.φ ≈ 0.5`, P and the mobility against independent closed forms, v = m(P − 101325)/d > 0, `F_mass[2, top] ≈ 2vξ + λhξ/(λ + hd)` to 1e-12 (the film term well above the tolerance), and the top gas row −A·F (rtol 1e-10). |
| 3 | Without a pressure condition only the film law acts. | ported | Same testset: with `DiffusiveMassBC(0.05)` alone the top flux is the film term λhξ/(λ + hd) (rtol 1e-12). |
| 4 | A char (κ > 0) and virgin (κ = 0) face carries exactly zero Darcy velocity while a char–char face flows. | ported | `test/newcore/test_gas_transport.jl` "gas transport: an interior Darcy–Fick face carries the Fick flux plus the upwind advection of ξ/φ, and no advection across the char front": at the char-front face (one side with `ds.mob == 0`) v == 0 and the flux is the Fick flux (rtol 1e-12), where an arithmetic face mobility would advect more than 10× the Fick flux; a char–char face has v ≠ 0. |
| 5 | Neumann ∂P/∂n is along the outward normal: outflow at both ends, magnitude (κ/μ)·abs(∂P/∂n). | ported | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form" (the Neumann velocity −m·∂P∂n along the outward normal) and "gas transport: a pressure condition on :bottom vents the bottom cell along the outward normal": `NeumannPressureBC(-1e5)` and `DirichletPressureBC(101325)` at `:bottom` give an outward v > 0, `F_mass[2, bottom] ≈ v·ξ/φ` (1e-12) and a bottom gas row that differs from the sealed deck's by −A·F (1e-10). |
| 6 | The convective pressure condition gives h_P(P − P_ext), zero on an impermeable cell, outward at the bottom. | covered | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form": zero at zero mobility, exact at a physical mobility, the smooth gate's derivative; the outward form makes the bottom sign intrinsic. |
| 7 | A Dirichlet overpressure vents outward at both ends, and an impermeable cell does not vent. | ported | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form": the Dirichlet velocity is m(P − P_bc)/d and 0 at m = 0; the vent is outward at the top in "gas transport: a pressure condition vents its cell …" and "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out", and at the bottom in "gas transport: a pressure condition on :bottom …". |
| 8 | A combined mass condition sums its members' velocities; non-pressure members drive none. | ported | `test/newcore/test_gas_transport.jl` "outward Darcy velocity of each mass condition in closed form": `DiffusiveMassBC` gives exactly 0, and `DiffusiveMassBC + DirichletPressureBC(9e4)` equals its Dirichlet member exactly. |

#### test_discretization.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | minmod, van Leer and superbee values, type stability, and the limiter enum dispatch. | dropped | The old ALE temperature-advection kernel and its limiters. The new reconstruction's slope is `test/newcore/test_relative_flux.jl` "van Albada slope". |
| 2 | `divergence_1d!` values: a constant flux gives 0, F = z gives 1. | dropped | The old `residual!` (`divergence_1d!` in `src/Discretization`). The new core scatters face fluxes into mass-unit rows. |
| 3 | The divergence telescopes: Σ div·V equals the boundary flux difference times A. | covered | `test/newcore/test_residual_jacobian.jl` "row assembly conserves mass and scatters fluxes antisymmetrically" (a closed slab; conduction cancels pairwise to 1e-12); `test/newcore/test_gas_transport.jl` "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, and gas leaving through :bottom leaves the bottom cell" (a boundary face's flux leaves its cell exactly). |
| 4 | `RHSCache` sizes. | dropped | The old `residual!` (`RHSCache`). |
| 5 | A linear steady profile has a constant interior conductive flux −kΔT/L and zero divergence. | covered | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …" (every interior face's conducted part is −k_f ΔT/d to 1e-12); `test/verification/cases/V-HC-04.jl` QoIs "L∞ steady profile vs (quantized) network" and "steady …-face flux vs network". |
| 6 | The blowing model reads the outward gas flux at both the top and the bottom face: outflow χ < 1, inflow χ > 1, a hotter surface under outflow. | ported | `test/newcore/test_boundary_conditions.jl` "blowing reads the outward gas flux of its face, at :top and at :bottom": at each tag a `HeatFluxBC + ConvectiveBC(…, blowing = CouetteFilm())` face with an outward gas flux of +0.01, 0 or −0.01 (`MassFluxBC(∓0.01)`, `ImpermeableBC`) has a `blowing_ratio` below 1, equal to 1 and above 1, and face temperatures ordered outflow > none > inflow, the same at both tags. |
| 7 | Manufactured sinusoidal conduction: interior fluxes within 5 %. | covered | The second-order conduction checks on the new core: `test/verification/cases/V-HC-03.jl` ConvergenceSpec (order 2) and `test/verification/cases/V-MMS-01.jl` QoI "L∞ T error, t=… s" with its ConvergenceSpec. |
| 8 | `update_darcy_properties!` fills P > 0, φ in [0, 1], κ ≥ 0, μ > 0 and the sized gradient and velocity arrays. | dropped | The old `residual!` and its Darcy sub-cache (`update_darcy_properties!`; `PropertyCache` is deleted with the old core, no caller outside it). The new state mapping's P, φ and mobility are exercised by V-GT-02 and "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out". |
| 9 | `compute_rhs!` with Darcy transport on a uniform state: no temperature change and zero solid tendency. | covered | `test/newcore/test_gas_transport.jl` "gas transport: isothermal transport leaves the temperature at rest (…)" (uniform T with a non-uniform gas: the temperature exactly at rest) and "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out" (zero solid fluxes). |
| 10 | At κ = 0 the Darcy path reduces to diffusion. | covered | `test/newcore/test_gas_transport.jl` "gas transport: …": decks with κ = 0 set up `FickTransport`, so the Darcy path never runs at zero permeability. |
| 11 | Darcy flow runs from high to low pressure (a hot top gives v < 0). | covered | `test/verification/cases/V-GT-02.jl` QoIs "exit mass flux vs ġL" (to 5e-5) and "P² profile error / overpressure² scale" (to 2.5e-3): a reversed v = −m∇P cannot reach the outflow parabola. |

#### test_gas_flux.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `gas_flux` at equal ξ and unequal T is non-zero, flows hot to cold, and equals −λ(ξ_R T_R − ξ_L T_L)/(T̄ d). | ported | `test/physics/test_gas_closures.jl` "Gas closures" > "gas flux: J = −λ(ξ_R T_R − ξ_L T_L)/(T̄ d)" > "equal concentrations at unequal temperatures flow hot to cold". |
| 2 | At uniform T `gas_flux` is Fick's −λΔξ/d. | ported | Same parent > "at a uniform temperature it is Fick's law". |
| 3 | The product form equals the expanded form −λ[∂ξ + (ξ̄/T̄)∂T] to 1e-14. | ported | Same parent > "the product form equals −λ[∂ξ + (ξ̄/T̄)∂T]" (1e-14). |
| 4 | `gas_flux` equals the partial-pressure form −(λM/(R T̄))∂P. | ported | Same parent > "it is the partial-pressure form −(λM/(R T̄))∂P". |
| 5 | The face loop gives non-zero, hot-to-cold interior gas fluxes under a temperature ramp, and zero solid fluxes. | covered | `test/newcore/test_gas_transport.jl` "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, and gas leaving through :bottom leaves the bottom cell" (the interior face equals `Physics.gas_flux` at the trial temperatures of a ramp to 1e-14); `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition" (the solid fluxes are zero). |
| 6 | `gas_flux_darcy_fick` is v·ξ_up plus diffusion, upwinded on the sign of v. | dropped | Deleted with the old core: `gas_flux_darcy_fick` has no caller outside it. |
| 7 | The harmonic face mean at equal distances is 2ab/(a + b), dominated by the small value. | covered | The mean survives as `guarded_harmonic_mean`: `test/newcore/jacobian/test_containment.jl` "the smooth face mean is finite with a finite derivative through its guard" (far above the guard it equals the distance-weighted harmonic mean to 1e-14). |
| 8 | The covering face-mean references do not use the deleted cut mean, and the equal-distance case is pinned. | ported | `test/newcore/jacobian/test_containment.jl` "the smooth face mean is finite with a finite derivative through its guard" (the reference is the closed form `(d_a + d_b)/(d_a/0.15 + d_b/0.45)`; `guarded(1e-6, 1e-4, d, d) ≈ 2e-10/(1e-6 + 1e-4)`, rtol 1e-10, and lies between 1.9e-6 and 2e-6) and `test/newcore/test_gas_transport.jl` "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, …" (the reference is the closed form `(d1 + d2)/(d1/λ1 + d2/λ2)`, rtol 1e-14, which the arithmetic mean misses by 2 %). |
| 9 | The face flux uses the harmonic mean of λ, not of λρ. | covered | `test/newcore/test_gas_transport.jl` "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, and gas leaving through :bottom leaves the bottom cell": the interior flux equals `gas_flux` with λ_f the distance-weighted harmonic mean of the cells' `ds.λ`, to 1e-14 (its reference is row 8). |
| 10 | `mass_bc_flux(DiffusiveMassBC)` uses λ directly: −λξ/d for ρ down to 1e-30, and −λhξ/(λ + hd) with a film. | ported | `test/boundary_conditions/test_mass_bcs.jl` "Mass boundary conditions" > "the diffusive film flux reads λ, not λρ": `mass_bc_flux(DiffusiveMassBC(), …) ≈ -λξ/d` for ρ in (1, 1e-20, 1e-30), and −λhξ/(λ + hd) with a film. |

#### test_nonuniform_gradients.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `opposite_distance_gradient` is exact on quadratics at any spacing, its coefficients reproduce it and sum to zero, and the wide and half-mean stencils are not exact. | dropped | The old core's `src/Geometry`; its only users are `S_conv`, the error indicators and the structured Jacobian, all deleted. |
| 2 | The `S_conv` gradient is exact on quadratics on a stretched mesh, boundary blend included. | dropped | `S_conv` (`Physics/convection.jl`). |
| 3 | Its successor on a stretched mesh: the enthalpy a face carries at the distance-weighted face temperature, with distance-weighted face means. | ported | `test/newcore/test_energy_ledger.jl` "energy row term by term against the ledger: …, stretched mesh" (Fick and Darcy, both forms, `create_stretched_mesh(L, 10, 1.2)`): the uniform mesh's 1e-12 face and 1e-11 cell bounds, with every interior face's two cell distances differing by more than 10 %. |
| 4 | The error-indicator gradient is second order on a stretched mesh. | dropped | Error indicators and h-AMR, deleted with Experimental (owner question 3). |
| 5 | `_limited_gradient` blends with β = d₀/(d₀ + d₁), exact on quadratics, and the VANLEER result stays within the TVD bound. | dropped | The old ALE temperature-advection kernel and its limiters. The new reconstruction is second order on a non-uniform moving mesh (`test/newcore/test_relative_flux.jl` "the relative fluxes converge at second order on a front at rest"). |
| 6 | The structured Jacobian equals dense AD on a stretched mesh (Fick, Darcy, ALE Fick). | dropped | The structured Jacobian. The coupling is contained on stretched meshes: `test/newcore/test_gas_transport.jl` "gas transport: …" (`create_stretched_mesh`, Fick and Darcy) and `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`fick`, `darcy`). |

### test/solver

#### test_boundary_flux_reporting.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The surface-gas reporter `_surface_gas_flux` equals the residual's boundary face flux, on fixed and ALE meshes. | dropped | The old `PyrolysisSolution` and its reporter (`src/Solver`). `CoreSolution` reads the boundary mass flux from the residual's own face fluxes, so only one path is left. |
| 2 | The ambient inflow of a species absent from the cell enters the boundary flux at the film value, at zero and at a slightly negative concentration. | ported | `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition": the top face's flux is `-0.05·0.5` at zero gas and `-0.05·(0.5 + 1e-12)` at a gas mass of `-1e-12·V` (the film law at the unclamped concentration), rtol 1e-12, and `Integration.solve(problem; use_ale).boundary_mass_flux[:top][1] ≈ -0.025` with `:Success` on fixed and moving meshes. |
| 3 | Under SURFACE_ABSORPTION the reported net flux of a combined condition equals the applied flux, with a_eff absorption. | covered | `test/boundary_conditions/test_boundary_temperature.jl` "surface Newton: a_eff absorption vs. two-piece split": the decomposition absorbs a_eff·q_inc and its net equals the solved applied flux, also through the face method. |
| 4 | `compute_boundary_heat_fluxes(...).net` equals the applied flux of `evaluate_bc_flux_and_derivative` for a combined and a standalone condition under NO_RADIATION, SURFACE_ABSORPTION and BEER_LAMBERT, with the pinned values: standalone under SURFACE_ABSORPTION net = absorbed = ε·q; combined under BEER_LAMBERT absorbed 0 and net = reradiation + convection, the parts summing to net; combined under SURFACE_ABSORPTION the absorbed flux included. | ported | `test/boundary_conditions/test_boundary_temperature.jl` "net is the applied face flux, per BC and per model": the same loop and pinned values through the face methods with a zero-outflow `BlowingFace`. |
| 5 | At every saved time the solved extras close the surface ledger (net ≈ conduction); under Beer–Lambert the deposit is positive and the face absorbs nothing; under surface absorption there is no deposit. | covered | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate": the decomposition sums, net ≈ conduction (rtol 1e-8), surface absorption absorbs ε·q with no deposit, Beer–Lambert absorbs nothing at the face with a positive deposit. |
| 6 | Under Beer–Lambert the deposit is bounded by ε·q_ext. | ported | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate" (its BEER_LAMBERT branch): `all(x.volumetric_absorption .<= x.surface_emissivity .* 4e4 .+ 1e-9)`. |
| 7 | With a ≠ ε the beam enters at a_eff·q_ext, and the cells' absorption plus the transmitted intensity equals it (residual level). | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": absorbed plus transmitted equals the entering beam, on a fixed mesh. The value a ≠ ε itself is pinned only by row 8. |
| 8 | The solution's extras report it: the deposit plus the transmitted beam ≈ a·q_ext, `surface_emissivity` ≈ ε and a positive deposit, on fixed and moving meshes. | ported | `test/newcore/test_solve.jl` "outputs: the in-depth deposit and the transmitted beam sum to the absorbed share: use_ale = …" (`false`, `true`): a solid with `ε = 0.95, absorptivity = 0.6` and a gas that takes no volume, under BEER_LAMBERT: `surface_emissivity ≈ 0.95`, the face absorbs nothing, the deposit is positive and at most a·q, the transmitted beam exceeds 0.1·a·q, deposit plus transmitted ≈ `0.6·q` against the literal a (rtol 1e-12), net ≈ conduction, and on the moving mesh the column thins by more than 0.5 %. |

#### test_default_linear_solver.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The default Jacobian is DirectSolve on a fixed mesh and ApproxSolve under ALE. | dropped | The structured Jacobian (`DirectSolve`, `ApproxSolve`); the colored Jacobian is the only one. |
| 2 | The default factorization follows the Jacobian pattern: UMFPACK for the ALE DirectSolve geometry block, KLU otherwise. | dropped | The structured Jacobian (the pattern rule of `_resolve_integrator`); the extra solvers go with the dependency removals (owner question 20). |
| 3 | A solve without `linear_solver` factors with the default, on fixed and moving meshes. | covered | `test/newcore/test_policies.jl` "one KLU analysis per growth of the union, use_ale = …": `default_linear_solver() isa StructuredKLU`, and the integrator cache of a `linear_solver = nothing` solve is a `StructuredKLUCache`, on both meshes. |
| 4 | An explicit `linear_solver` is the one the integrator uses. | ported | `test/newcore/test_policies.jl` "one KLU analysis per growth of the union, use_ale = …" (both meshes): with `KLUFactorization()` and `UMFPACKFactorization()` the integrator's `alg.linsolve` is the given type, the KLU cache is not a `StructuredKLUCache` and the UMFPACK one is an `UmfpackLU`. |

#### test_ilu_gmres.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | ILU-GMRES solves an ALE deck, matches the KLU baseline to 1 K, and succeeds at τ = 0.01 and 0.5. | dropped | ILU/GMRES, among the dependency removals (owner question 20); the `integrator` keyword is among the removed keywords (owner question 18). |
| 2 | Constructor defaults, `show`, and the LinearSolve interface of both Krylov algorithms. | dropped | ILU/GMRES (owner question 20). |
| 3 | An unconverged Krylov solve reports `ConvergenceFailure` to the Newton. | dropped | ILU/GMRES (owner question 20). |

#### test_initialize_state.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The initial state is the problem's callable initial fields at the cell centres, not a mesh placeholder. | ported | `test/newcore/test_problem.jl` "the set-up state holds the initial fields at the cell centres": on a stretched mesh with non-uniform callable fields, `cp.T0`, the energy rows and the mass rows ξ·V (rtol 1e-12) equal the file's recorded cell-centre values (not `initial_fields`, which `setup` itself calls), and so do the solve's first saved T and ξ. |
| 2 | The keywords reach the definition: mode, radiation model, radiation intensity, the ALE node block from the live mesh. | dropped | `initialize_state` and the old `ProblemDef`/`Workspace`. The incident intensity is a `RadiativeFluxBC`, and the new solve refuses unknown keywords. |
| 3 | Under a lateral shrinkage law the dry-solid denominator is never zero (the `m_dry_solid_initial = 0` trap). | covered | `test/newcore/test_mesh_motion.jl` "lateral rows: the lateral rate is the area row's logarithmic rate": `LateralShrinkage(mat, 2.0, 0.0)` throws `ArgumentError`, so a zero denominator is refused. |
| 4 | That denominator is V_i Σ_solid ξ_j(0) from the problem's initial fields. | ported | Same testset: under a lateral law with `use_ale = true`, `has_chi(layout)`, and in every cell G = 0 and M = ξ_solid·V (rtol 1e-12), the gas the cells hold excluded. |

#### test_integrator_type.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Decks that differ in components, properties, reactions, conditions and transport run on one integrator type and one solution type per solver configuration. | covered | `test/newcore/test_policies.jl` "the integrator's type does not depend on the problem: …": plain and rich decks on a fixed mesh, a fixed mesh with Beer–Lambert, a moving mesh and a moving deactivating mesh, asserting `===` on the integrator and solution types. |
| 2 | The same across radiation models on a fixed mesh (plain NO_RADIATION and rich BEER_LAMBERT on one type). | dropped | By design, pending the owner's decision: on a fixed mesh Beer–Lambert adds algebraic intensity rows and a mass matrix, so it compiles a second integrator (Section 4). |

#### test_jacobian_kwarg.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | A bare `Structured` backend is built for the solve's mode. | dropped | The `jacobian` keyword, among the removed keywords (owner question 18), and the structured Jacobian. |
| 2 | Prebuilt specs and the dense AD reference are rejected with `TypeError`. | dropped | The same keyword and the structured Jacobian. |

#### test_newton_failure_cap.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The cap's constants are pinned: growth 1.25, 20 steps. | ported | `test/newcore/test_policies.jl` "Newton-failure step cap": `FAILURE_MEMORY_GROWTH == 1.25` and `FAILURE_MEMORY_STEPS == 20`. |
| 2 | After a Newton failure the proposed step stays at or below growth × the failing step for the memory's accepted steps, the cap binds, and then it releases. | covered | `test/newcore/test_policies.jl` "Newton-failure step cap": the same three properties with the cap's step record. |
| 3 | Without a failure the capped trajectory is the free one. | covered | Same testset: no nonlinear failures, and `t` and `u` identical with and without the cap. |

#### test_precompile_workloads.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Every precompile workload runs to its end, the names are unique, and the depletion workloads deactivate cells. | port | Written in the deletion commit with the workloads it tests (C121; Section 4): `test/newcore/test_precompile_workloads.jl`, in which the names of `Pyrolysis._precompile_workloads()` are unique, each runs without a warning to `:Success` or `:Terminated`, and each deactivation workload has events. |

#### test_problem.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The keyword and positional constructors keep their arguments (mesh, material, tspan, component count). | ported | `test/newcore/test_problem.jl` "problem construction: the arguments kept, the conditions and output times filled": on `create_uniform_mesh(0.01, 20)` the keyword and positional constructors keep `mesh` and `material` (`===`), `tspan` and the component count. |
| 2 | `BoundaryConditionSet` fills missing conditions (thermal bottom `AdiabaticBC`, mass `ImpermeableBC`), and the problem validates. | ported | Same testset: a set given only the top's thermal condition gets an `AdiabaticBC` bottom and `ImpermeableBC` faces, and the problem passes `validate_problem` and `is_valid_problem`. |
| 3 | Custom callable initial fields are stored and evaluated at the cell centres. | covered | `test/newcore/test_problem.jl` "initial fields: the profiles at the cell centres": `initial_fields` equals the recorded cell values. |
| 4 | Empty `output_times` are generated from `tspan[1]` to `tspan[2]`. | ported | Same testset: for `tspan = (0, 100)` the generated output times are sorted, start at 0 and end at 100. |
| 5 | `show` of the problem, both forms, names the type, the material and its components. | ported | Same testset: both `show` forms contain "PyrolysisProblem" and the material name, the compact one "20 cells", the full one "Material: inert" and "Components: 2". |
| 6 | A valid problem passes `validate_problem` and `is_valid_problem`. | covered | `test/newcore/test_problem.jl` "validation on the native mesh: the same refusals". |
| 7 | An inverted `tspan` is refused. | covered | Same testset. |
| 8 | A permeable matrix with φ₀ = 0 is refused, and the message names `ρ_intrinsic`. | covered | Same testset. |
| 9 | The rest of the Darcy rule: `is_valid_problem` false; the message names `background_gas_concentration` and "porosity"; `ρ_intrinsic = 1500` restores validity; a zero-porosity sub-region is refused; a material whose initial κ_eff is 0 passes. | ported | `test/newcore/test_problem.jl` "validation on the native mesh: the same refusals": `is_valid_problem` is false; the message names `background_gas_concentration` and "porosity"; `ρ_intrinsic = 1500` restores validity; a region with no pore space is refused; a matrix permeable only in its char passes. |
| 10 | An impermeable non-porous (Fickian) material passes. | covered | `test/newcore/test_problem.jl` "validation on the native mesh: the same refusals": its valid problem is Fickian, without `ρ_intrinsic`. |
| 11 | Degenerate Fickian gas transport (a gas-producing reaction, λ_solid = 0, no κ) warns "gas-transfer coefficient"; λ_solid > 0, κ > 0 or no reaction keeps validation silent. | ported | `test/newcore/test_problem.jl` "validation warns when Fickian gas transport lacks a condensed transfer coefficient": the warning for a gas-producing reaction with λ_solid = 0 and no κ, and silent validation with λ_solid > 0, κ > 0 or no reaction (`@test_logs`). |
| 12 | `apply_initial_conditions!(problem)` writes the cell states and `m_dry_solid_initial` into the mesh. | dropped | Deleted with the old core: `apply_initial_conditions!` has no caller outside it (it writes the old `Mesh1D`'s cell states). The values are row 3; the denominator is `test_initialize_state.jl` rows 3 and 4. |
| 13 | `apply_initial_conditions!(mesh, material; T, ξ)` writes uniform states. | dropped | Deleted with the old core: this `apply_initial_conditions!` method takes the old `Mesh1D` and has no caller outside it. |
| 14 | The step, linear, gaussian and tanh profile helpers give their values. | covered | `test/newcore/test_problem.jl` "initial fields: the profiles at the cell centres": the fields are built from all four helpers and equal the recorded values bit for bit. |
| 15 | An empty `PyrolysisSolution(problem)` and its `show`. | dropped | The old `PyrolysisSolution`. `CoreSolution`'s `show` is in `test/newcore/test_outputs.jl` "accessors on a moving deactivating deck". |
| 16 | A basic solve succeeds, saves the requested times, and the surface heats under a flux. | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (success and a field per saved time); the heating: `test/verification/cases/V-HC-01.jl` QoI "surface temperature, t=… s" (against Carslaw–Jaeger within 1.5 K). |
| 17 | The reported `solve_time` is positive. | ported | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances": `sol.solve_time > 0`. |
| 18 | `abstol` and `reltol` are accepted. | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" passes `abstol = 1e-10, reltol = 1e-8` and succeeds. |
| 19 | `get_profile(sol, :T or :ξ1, k)` has one value per cell. | covered | `test/newcore/test_outputs.jl` "accessors on a moving deactivating deck": `get_profile` equals the saved arrays for `:T`, `:ξ` and `:ξ2`. |
| 20 | `get_cell_positions` has one ascending entry per cell. | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (the exact centres) and `test/newcore/test_outputs.jl` "accessors on a moving deactivating deck". |
| 21 | The full workflow: the solid is consumed, gas is produced, the surface heats. | covered | `test/newcore/test_solve.jl` "solve: kinetics through the adapter, both energy forms" (exact first-order decay) and "outputs: the surface-flux decomposition and the mass-loss rate" (the gas leaves: a positive final mass-loss rate); the heating: `test/verification/cases/V-HC-01.jl` QoI "surface temperature, t=… s". |

#### test_progress_bar.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `progress = true` takes exactly the steps of `progress = false`, on fixed and moving meshes. | covered | `test/newcore/test_outputs.jl` "the progress bar is inert: …" (a fixed deck and a moving deactivating deck, under the default, Rodas5P, FBDF and KenCarp4 policies): the hash of T, m, z, t, the event times and the statistics, and the retcode, are equal with and without the bar. |

#### test_repeated_solve_determinism.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Repeated solves of one problem in one process are bit-identical under each sparse factorization, on a moving mesh and with deactivation. | covered | Without radiation: `test/newcore/test_policies.jl` "repeated solves are bit-identical under each linear solver" (StructuredKLU, KLU, UMFPACK; moving; `t`, `u` and accepted counts equal) and `test/newcore/test_correction.jl` "repeated solves with the correction and events are bit-identical under each linear solver" (with events). |
| 2 | The same under Beer–Lambert, whose intensity rows change W's pattern, on fixed and moving meshes. | ported | `test/newcore/test_policies.jl` "repeated solves with in-depth radiation are bit-identical under each linear solver: use_ale = …" (`false`, `true`): `_slab(12; radiating = true)` under BEER_LAMBERT, solved twice under each of StructuredKLU, KLU and UMFPACK, with `raw_solution.t`, `raw_solution.u` and the accepted counts equal. |

#### test_solve_leaves_problem.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | After a moving solve that deactivates cells, the caller's mesh positions are unchanged. | covered | `test/newcore/test_outputs.jl` "accessors on a moving deactivating deck": after the moving, deactivating solve, the midpoints of `problem.mesh.positions` equal the saved initial cell centres exactly. |
| 2 | The topology arrays the solve shares with the caller are unchanged too. | ported | `test/newcore/test_solve.jl` "solve leaves the caller's mesh unchanged on a moving, deactivating deck": after a 12-cell ablating solve with events and a thinned column, the mesh positions and every `MeshTopology` field are `isequal` to a `deepcopy` taken before it. |
| 3 | A second solve of the same problem repeats the first bit for bit, on a moving mesh and with events. | covered | `test/newcore/test_correction.jl` "repeated solves with the correction and events are bit-identical under each linear solver" (one problem, events, `t`, `u` and event times equal); `test/newcore/test_policies.jl` "repeated solves are bit-identical under each linear solver". |
| 4 | `sol.problem !== problem`: the solve integrates a copy. | dropped | The old solve's problem copy (`src/Solver`). The new solve writes nothing into the problem and `CoreSolution` keeps the caller's; rows 1 and 2 are the contract. |
| 5 | A warm start from a parked solve is rejected. | dropped | `warm_start`, among the removed keywords (owner question 18), and parking, replaced by deactivation. |
| 6 | A `parameter_sweep` whose `modify` reuses the base mesh leaves that mesh unchanged; each point equals a fresh solve, and different points differ. | ported | `test/newcore/test_parameter_sweep.jl` "parameter sweep: points that share the base mesh leave it unchanged": a serial sweep whose `modify` builds every point on `problem.mesh` (an 8-cell ablator, `use_ale = true, deactivation = true`, the grid `[A, 2A, A]`, `extract` returning t, T, m, z and the event times): the base mesh's positions and topology are unchanged, each point `isequal`s a fresh solve on a fresh mesh and has events, and two different points differ. |

#### test_surface_bc_guards.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | When its iteration budget runs out the surface Newton warns "did not converge". | dropped | Deleted with the old core: the allocating `solve_boundary_temperature` and the warnings of the `trial === nothing` path of `solve_boundary_temperature!` have no caller outside it. The surviving path records an exhausted budget on the state without a warning (row 2). |
| 2 | It records `converged == false` and the iteration count, and converges on the default budget. | ported | `test/boundary_conditions/test_boundary_temperature.jl` "surface Newton: an exhausted iteration budget is recorded": `solve_boundary_temperature!` with `trial = trial_temperature` and `max_iter = 1` on a radiative plus convective face records `converged == false` and `iterations == 1` with an imbalance above 1 W/m² and no warning; the default budget converges to 1e-8. |
| 3 | The residual carries no surface-Newton state between evaluations (a wild or NaN trial, then a good state, equals a fresh workspace), fixed mesh. | covered | `test/newcore/test_boundary_conditions.jl` "the residual carries no surface-balance state between evaluations", both energy forms. |
| 4 | The same on a moving mesh. | ported | `test/newcore/test_boundary_conditions.jl` "the residual carries no surface-balance state between evaluations": the loop runs for `use_ale in (false, true)` in both forms, comparing the residual's vectors at rtol 1e-12 (Section 4: the enthalpy form's warm start). |
| 5 | A final iterate pinned at the [200, 3000] K clamp warns and returns T_s = 3000. | dropped | Deleted with the old core: the iterate clamp and its warning belong to the `trial === nothing` path of `solve_boundary_temperature!`, which has no caller outside it. The surviving path neither clamps nor warns (`test/newcore/test_boundary_conditions.jl` "the surface balance reads the trial clamp: no clamp, no warning, finite derivatives"). |
| 6 | The mass conditions write only gas rows of the boundary face flux; the solid rows carry no film flux. | covered | `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition" (zero solid flux under a film law) and `test/newcore/test_gas_transport.jl` "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out". |
| 7 | `ambient_concentration` takes a scalar for every species, a tuple per species, or a function of time; a negative ambient is refused at construction. | ported | `test/boundary_conditions/test_mass_bcs.jl` "Mass boundary conditions" > "ambient concentrations: one for every species, one per species, or a function of time": `ambient_concentration` of a scalar, a per-species tuple and a function of time, and `ConvectiveMassBC(h_m = 0.05, ξ_inf = (-0.1, 0.0))` throwing `ArgumentError`. |
| 8 | A tuple of the wrong length is refused at validation. | covered | `test/newcore/test_problem.jl` "validation on the native mesh: the same refusals". |
| 9 | Ambient gas enters a gas-free column, and the solid sees no ambient solid. | covered | `test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition": the gas row is positive and every solid row is zero. |
| 10 | Venting through both faces closes against the mass lost. | covered | `test/newcore/test_boundary_conditions.jl` "mass leaving through both faces closes against the inventories": bottom ≈ top, and the integrated outflow equals the mass lost (rtol 1e-3). |
| 11 | `mass_loss_rate` counts the bottom face. | ported | `test/newcore/test_boundary_conditions.jl` "mass leaving through both faces closes against the inventories": `sol.mass_loss_rate` equals the top plus the bottom outflow (rtol 1e-12), with the bottom flux equal to the top one (rtol 1e-6), so a top-only rate would be off by half. |
| 12 | A standalone Dirichlet face is the reported surface temperature. | covered | `test/newcore/test_boundary_conditions.jl` "a prescribed face temperature is the reported surface temperature": `surface_T` is 600 K at every saved time. |
| 13 | Its reported net flux is finite and, from a 300 K start, flows in. | ported | `test/newcore/test_boundary_conditions.jl` "a prescribed face temperature is the reported surface temperature": `all(isfinite, sol.extras.surface_heat_flux_net)` and `sol.extras.surface_heat_flux_net[2] > 0`. |
| 14 | Under BEER_LAMBERT a bottom `RadiativeFluxBC` is refused; SURFACE_ABSORPTION accepts it. | covered | `test/newcore/test_boundary_conditions.jl` "in-depth radiation enters through the top faces only". |
| 15 | A time-varying incident flux equal at t = 0 and t = 1 is resampled at every evaluation. | covered | `test/newcore/test_boundary_conditions.jl` "a time-varying incident flux is sampled at every evaluation": the varying flux's surface is more than 1 K hotter at t = 0.5. |

#### test_time_integration.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The old block layout: T and ξ blocks, the z block under ALE, no χ block. | dropped | The old `StateLayout`. The new interleaved layout is `test/newcore/test_states.jl` "state layout: interleaved blocks and row kinds". |
| 2 | The residual is finite at the initial state, on fixed and moving meshes. | covered | `test/newcore/test_residual_jacobian.jl` "residual and colored Jacobian: …" (fixed, at `u0`) and `test/newcore/jacobian/test_containment.jl` "perf configuration …" (the moving configurations at `u0`). |
| 3 | `compute_jacobian!` of `Structured(DirectSolve)` is finite and square, Eulerian and ALE. | dropped | The structured Jacobian (`compute_jacobian!`). |

### test/perf

#### perf_helpers.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | No assertions: it provides the canonical configurations and the mixed reacting state the invariants run on. | covered | Rebuilt as problems in `test/newcore/jacobian/test_allocations.jl` (`_configurations`, run by "allocations: …") and `test/newcore/jacobian/test_containment.jl` (`_perf_configurations` and the `mixed` state, run by "perf configuration …"). |

#### runtests.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The allocation probes run hermetically and still measure zero. | covered | `test/newcore/jacobian/test_allocations.jl` "allocations: …": measured inside a function in the shared test process; process state can only add allocations, never hide them. |
| 2 | The quality invariants run in a fresh Julia process and are skipped when JET is unavailable. | covered | `test/newcore/quality/runtests.jl`, run by `test/runtests.jl` in a fresh process and excluded from the in-process walk, with the JET-availability skip; `test/newcore/quality/type_analysis.jl` "type analysis: …". The method count still runs in-process (Section 4). |

#### test_ad_jacobian_correctness.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The AD reference Jacobian equals a centred finite-difference Jacobian of the residual (max error over max abs(J) below 1e-6) on the Eulerian and ALE baselines. | ported | `test/newcore/jacobian/test_containment.jl` "the dense Jacobian equals central differences: …" (`baseline_eulerian`, `ale_no_radiation`, 8 cells): the dense AD Jacobian equals central differences of `residual!` per entry with the row-magnitude floor at rtol 1e-6, which implies the obligation's global criterion, in the temperature form at `u0` and `mixed` and in the enthalpy form at `mixed`. At `u0` the gas masses sit on the joint of their positive part, which the difference step straddles (C129). |

#### test_jacobian_allocations.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Jacobian assembly allocates within a budget after warm-up on the five canonical configurations. | covered | `test/newcore/jacobian/test_allocations.jl` "allocations: …": `jacobian!` allocates at most `JACOBIAN_BUDGET = 0` bytes on those five and three more, both forms. |

#### test_method_count.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The residual's method table stays bounded (no trait explosion, at most 16). | covered | `test/newcore/jacobian/test_allocations.jl` "one residual method": `length(methods(Pyrolysis.Assembly.residual!)) == 1`, stricter; it runs in-process (Section 4). |

#### test_type_stability.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | JET's `@test_call` on the concrete core residual call reports nothing in the core's modules. | covered | `test/newcore/quality/type_analysis.jl` "type analysis: …": `@test_call` on `residual!` and the colored `jacobian!` at the concrete types `setup` builds (a fixed PMMA column, a moving Beer–Lambert Darcy–Fick column, a moving equidistribution column), both forms, in a fresh process. |

#### test_verify_jacobian.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The dense AD reference passes against itself. | dropped | The structured Jacobian (`verify_jacobian_at`, `DenseAD_ForwardDiff`). |
| 2 | The production Jacobian equals the dense AD reference on the structural support, on every canonical configuration. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …": all seven configurations at 6 and 8 cells, both forms, eight states; containment and colored equal to dense per entry at rtol 1e-12. |
| 3 | ApproxSolve's masked comparison on the ALE configurations. | dropped | The structured Jacobian (`ApproxSolve`, `is_approximate`). |
| 4 | On an optically moderate Beer–Lambert column at a reacting state, the energy rows' dependence on every node above is in the pattern. | covered | `test/newcore/jacobian/test_containment.jl` "perf configuration …" (`ale_beer_lambert_darcy_moderate`, `ale_beer_lambert_fickian_moderate`, state `mixed`): no missing entry, per-entry equality. |

#### test_zero_allocation_rhs.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The residual allocates zero bytes after warm-up on the baseline, ALE, surface-absorption and ALE Beer–Lambert Darcy configurations. | covered | `test/newcore/jacobian/test_allocations.jl` "allocations: …": the residual allocates 0 bytes on all of them and four more, both forms, measured inside a function. |

### test/adaptivity

#### test_adaptivity.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `GradientIndicator` is non-negative and nearly constant on a linear field. | dropped | Error indicators, deleted with Experimental (owner question 3). |
| 2 | Indicator variable resolution: an `Int` tracks that component; unknown symbols and out-of-range indices throw. | dropped | Error indicators (owner question 3). |
| 3 | `CurvatureIndicator` is small on a linear field and exact (0 and 2) on a non-uniform 3-cell stencil. | dropped | Error indicators (owner question 3). |
| 4 | `ReactionRateIndicator` needs rates, peaks at the reacting cell, and rejects stale rates. | dropped | Error indicators (owner question 3). |
| 5 | `CompositeIndicator`'s max-normalized `:max` and `:sum` combinations. | dropped | Error indicators (owner question 3). |
| 6 | `InterfaceProximityIndicator` peaks at the interface. | dropped | Error indicators (owner question 3). |
| 7 | h-refinement of one cell: topology, hierarchy, level, volume conservation. | dropped | h-AMR, deleted with Experimental (owner question 3). |
| 8 | h-refinement conserves each component's mass. | dropped | h-AMR (owner question 3). |
| 9 | Batch refinement marks high-indicator cells, and 2:1 balance keeps neighbour level differences at most one. | dropped | h-AMR (owner question 3). |
| 10 | h-coarsening: topology, mass conservation, `can_coarsen`, pair marking. | dropped | h-AMR (owner question 3). |
| 11 | r-refinement: the monitors are at least one and peak at the interface; relocation moves interior nodes and keeps the end nodes and positive volumes. | dropped | r-refinement, deleted with Experimental (owner question 3). |
| 12 | Mesh quality metrics (aspect ratio, smoothness, Jacobian, composite, summary) on a uniform mesh. | dropped | Mesh quality metrics, deleted with Experimental (owner question 3). |
| 13 | `AMRController`: defaults, `adapt_mesh!` result and quality, `should_adapt` cadence, `maybe_adapt!`, statistics, the default controller. | dropped | h-AMR (owner question 3). |
| 14 | The refine-coarsen cycle and a full AMR cycle conserve mass and keep mesh quality. | dropped | h-AMR (owner question 3). |

#### test_cell_park.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Parking the surface cell on a moving mesh conserves solid and gas mass, sensible enthalpy, live volume and thickness. | covered | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (moving): masses plus the lost mass (zero) at 1e-14, enthalpy plus the lost enthalpy at 1e-12, no created volume, the top node position kept. |
| 2 | After a park the live range is 1:n_live, the exposed boundary moves to the new top face, and dead cells hold nothing. | covered | Same testset: the top live cell, zero dead masses, dead rows masked and the receiver live, and the `:top` condition acting on the receiver's far face. |
| 3 | The dead node sits on the surface and the dead cell has zero thickness. | covered | Same testset: the dead node at the old top position, coincident with its neighbour, and its rows masked. |
| 4 | Nodes inward of the merge are untouched, and the live nodes stay strictly increasing. | ported | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (moving) and "a run below the threshold merges into the cell beneath it in one invocation: …": the positions of the receiver's inward node and of every node beneath it are bitwise unchanged by `deactivate!`. |
| 5 | A park records its time and cell, and a saved state's live count follows the park times. | covered | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (the event's time, invocation, cell and receiver) and "a terminated solve saves its last state, and a state saved at an event's time is the one after it" (the state saved at an event's time is the post-event one; Section 4). |
| 6 | Parking an interior cell on a moving mesh shifts the cells above it down one slot and compacts the reference. | dropped | Parking (interior parking). Cells die from the exposed side only; a cell below the threshold under a live one gives `:interior`. |
| 7 | Parking the bottom cell merges it into the cell above. | dropped | Parking. The replacement, a run reaching the innermost cell leaving it as the survivor, is `test/newcore/test_deactivation.jl` "a column below the band leaves its innermost cell as the survivor: …". |
| 8 | A one-cell mesh cannot be parked (the termination floor). | dropped | Parking, and `min_cells` among the removed keywords (owner question 18). The replacement: the last cell goes with no receiver ("a column below the band leaves its innermost cell as the survivor: …") and the solve ends `:burned_out` ("a terminated solve saves its last state, and a state saved at an event's time is the one after it"). |
| 9 | With a lateral law, parking conserves Σχ_i m_dry,i and Σ m_dry_init, and the dead cell's χ is zero. | ported | `test/newcore/test_deactivation.jl` "a merge carries the history pair to the receiver and keeps the cross section: …" (both forms; histories Gᵢ = (0.1 + 0.05i)Mᵢ, so dropping a cell's pair would move the area): ΣG and ΣM unchanged (1e-14), the receiver holds both cells' pairs exactly and the dead cell's pair is zero, the area unchanged (1e-12) with the area and lateral-rate residuals at round-off, and containment with colored equal to dense. |
| 10 | On a fixed mesh an interior cell cannot be parked. | covered | `test/newcore/test_deactivation.jl` "a void under a live cell ends the solve" (fixed and moving: `:Terminated`, outcome `:interior`, no events) and "a run below the threshold merges into the cell beneath it in one invocation: …" (`deactivate!` returns `:interior`). |
| 11 | On a fixed mesh a surface park conserves masses and enthalpy and moves no volume. | covered | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (fixed: masses plus the lost mass and enthalpy plus the lost enthalpy conserved, the dead cell's gas vented with the event; Section 4) and "on a fixed mesh the receiver keeps its volume and none of the dead cell's gas". |
| 12 | The depletion reference is per cell (the initial thickness on a moving mesh, the initial condensed concentration on a fixed one). | ported | `test/newcore/test_deactivation.jl` "each cell is held to its own initial condensed mass: …" (fixed and moving; 12 cells stretched at 1.35, a cell-size ratio of 27, ε_m = 0.05): the reference equals each cell's initial condensed mass and nothing is below at `u0` for both `stretch_from`; with `:top`, cell 1 at 0.04 of its own reference is the only cell below, and at 0.06 none is. |
| 13 | Candidates come surface-most first, and dead cells are never scanned. | covered | `test/newcore/test_deactivation.jl` "a run below the threshold merges into the cell beneath it in one invocation: …" (events outermost first; masked cells skipped on the next invocation) and "a transfer conserves the inventories and the receiver absorbs the volume: …". |
| 14 | `find_thin_cells` rejects a reference vector of the wrong length. | dropped | Parking (`find_thin_cells`). |
| 15 | A tangled (inverted) cell is refused, not taken as depleted. | covered | `test/newcore/test_policies.jl` "domain check: admissibility, porosity and constraint residuals" (crossed nodes are inadmissible); the new trigger is a mass fraction and never reads a cell's thickness. |
| 16 | On a fixed mesh only the top live cell is a candidate, and only below the threshold. | covered | `test/newcore/test_deactivation.jl` "the near-threshold band" (on a fixed mesh only the outermost cell is taken and the pass stops) and "a void under a live cell ends the solve". |
| 17 | The interior-void test is the matrix volume fraction below `MIXING_BLEND_LO` (a thin porous cell is not a void; the lowest void is reported). | dropped | Parking (`interior_void_cell` and its criterion). The new criterion is the mass threshold under a live cell ("a void under a live cell ends the solve"; `test/verification/cases/V-BO-05.jl` QoI "outcome :interior"). |
| 18 | An empty top live cell, also after a park, takes the park path, not the void path, and dead cells are not scanned. | covered | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (`:transferred` for a top cell below the threshold) and "a run below the threshold merges into the cell beneath it in one invocation: …" (the pass after an event). |
| 19 | A cell that started with no condensed mass is never tested. | ported | `test/newcore/test_deactivation.jl` "a cell that started with no condensed mass is never tested: …" (fixed and moving): cell 3 starts with gas only and `reference[3] == 0`; with a trace below zero in its solid it is not below the threshold and the margin is (1 − ε_m)/ε_m; with cell 6 depleted the margin is −0.5 and `deactivate!` returns `:transferred` (6 into 5); on the moving mesh a solve succeeds with every event outward of cell 3. |
| 20 | The mesh mirrors the state (cell temperatures and concentrations) and keeps its live topology (each cell's node pair). | dropped | The old `Mesh1D` and `CellState`; the new topology is immutable. |
| 21 | After a park exactly one face and one node carry the top boundary tag. | covered | `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …": the live top face is the receiver's far face (`live_face`), on both meshes. |

#### test_mesh_merge_guards.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Constant c_p: the merged temperature is the heat-capacity-weighted mean. | covered | `test/newcore/test_deactivation.jl` "a pore-free donor and a porous receiver" (inside "a merge never destroys pore space": the merged temperature is 300 + 0.951·200/1.951) and "on a fixed mesh the receiver keeps its volume and none of the dead cell's gas" (the mass-weighted mean), with equal heat capacities. |
| 2 | The pore gas stores sensible heat in the merge. | covered | `test/newcore/test_deactivation.jl` "a pore-free donor and a porous receiver": the receiver's 1e-3 kg of gas enters the weight (0.951, not 0.95), atol 1e-10. |
| 3 | Components with different heat capacities are weighted by their own c_p (gas at 5× the solid's). | ported | `test/newcore/test_deactivation.jl` "a merge conserves the sensible enthalpy of every component" > "unequal heat capacities: each component weighted by its own": with the gas at 5× the solid's heat capacity, a moving-mesh pair merges at Σ m_j c_j T/Σ m_j c_j (atol 1e-10), and Σ_j m_j h_j is conserved (1e-12). |
| 4 | With a linear c_p(T) the merge conserves Σ m ∫c_p dT (the exact root), not the product form (off by more than 1 K). | ported | Same parent > "a linear heat capacity: the enthalpy's root, not the product form": with `c = LinearProperty(500, 2)` and equal masses at 800 and 400 K the merged temperature is the root of the summed enthalpy (1e-6 K), Σ m h is conserved (1e-12), and the product form is more than 1 K off. |
| 5 | When both cells are empty the merge falls back to the volume-weighted mean. | dropped | Mesh merging, replaced by deactivation (`merge_temperatures_conservatively`). The new merge inverts the caloric extension, whose degenerate regime is "a transfer writes the regularization the receiver does not hold: …" and `test/newcore/test_states.jl` "extended caloric map: inversion, implicit derivatives, exact region". |
| 6 | `find_thin_cells` refuses an inverted cell. | covered | `test/newcore/test_policies.jl` "domain check: admissibility, porosity and constraint residuals" (crossed nodes are inadmissible); the deactivation trigger reads masses, not thickness. |
| 7 | The depletion condition keeps firing at the `min_cells` floor. | dropped | `min_cells`, among the removed keywords (owner question 18). The last cell goes with no receiver ("a column below the band leaves its innermost cell as the survivor: …"). |

#### test_surface_depletion.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The depletion trigger flags a cell below the threshold of its reference, and not a healthy one. | covered | Now a condensed-mass fraction: `test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …" (below at 0.5 ε_m), "a run below the threshold merges into the cell beneath it in one invocation: …" (the receiver is not below), "the near-threshold band" (nothing below at 5 ε_m) and "the end-of-step safety net takes a crossing location cannot see: …" (a zero margin at the threshold). |
| 2 | The reference is per cell: on a stretched mesh with a cell ratio of 27 nothing is flagged at t = 0, and a cell at 4 % of its own reference is flagged, at 6 % not. | ported | `test/newcore/test_deactivation.jl` "each cell is held to its own initial condensed mass: …" (the test of `test_cell_park.jl` row 12). |
| 3 | `find_thin_cells` rejects a reference vector of the wrong length. | dropped | Parking (`find_thin_cells`). |
| 4 | The thickness termination ends a shrinking solve at `min_thickness`. | covered | `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion": `:Terminated`, with the thickness and the time against the closed form at 1e-6. |
| 5 | After a park the termination reads the live thickness (the outermost live node). | ported | `test/newcore/test_outputs.jl` "succeeded and completed by return code and outcome": the deactivating deck stopped at `min_thickness = 2e-3` has events, fewer live cells at the end than at the start, and `thickness_history[end] == z[live+1, end] - z[1, end]`, the live span. |
| 6 | A saved state's live cell count follows the park times. | covered | `test/newcore/test_deactivation.jl` "a solve through deactivations: …" (`live_cells` follows the event count) and "a terminated solve saves its last state, and a state saved at an event's time is the one after it" (the post-event state; Section 4). |
| 7 | The step after a park starts from a reset step size: under a forced restart. | covered | `test/newcore/test_deactivation.jl` "the restart after an event: …" (the forced-restart policy: the step after each event equals `DT_AFTER_EVENT`). |
| 8 | The same on the default solve path. | ported | `test/newcore/test_deactivation.jl` "the heuristic restart replaces the step the solver would keep after an event": `_ablator(30)` on a moving mesh under `rodas5p_policy(KLUFactorization())`, whose restart is the default heuristic, and under the same policy with `restart = :solver_chosen` takes identical steps up to the same first event and a different step after it. The KLU policy is used because the default StructuredKLU moves the event by 1e-15 relative between the two runs. |
| 9 | Park times are recorded on a fixed mesh, ordered, inside the span, with the live count. | covered | `test/newcore/test_deactivation.jl` "a solve through deactivations: …" (fixed: the event times, cells counting down from 30 outermost first, `live_cells`); `test/verification/cases/V-BO-02.jl` QoI "park events" (the count). |
| 10 | An interior void stops a fixed-mesh solve before any park, with the bottom cell empty and the top still solid. | covered | `test/newcore/test_deactivation.jl` "a void under a live cell ends the solve" (`:Terminated`, `:interior`, no events, t below the end, both meshes); `test/verification/cases/V-BO-05.jl` QoIs "stop time at the first interior crossing of ε_m", "outcome :interior", "deactivation events before the stop". |
| 11 | The interior-void stop warns ("gasified from within"). | covered | The stop is reported instead: `test/newcore/test_outputs.jl` "succeeded and completed by return code and outcome" (`:Terminated`, outcome `:interior`, neither succeeded nor completed) and `test/newcore/test_deactivation.jl` "a void under a live cell ends the solve". The new core does not warn (Section 4). |

### test/ale

#### test_advection.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | An inert column on the moving-mesh layout (no volume change): the mesh stays put and the temperatures and surface temperature equal the fixed-mesh solve's. | ported | `test/newcore/test_mesh_motion.jl` "solve: an inert column on a moving mesh stays put and matches the fixed mesh": an inert solid on `create_uniform_mesh(0.02, 40)` under a `DirichletThermalBC(500.0)` top, solved for 50 s with and without `use_ale` in both forms (abstol 1e-10, reltol 1e-8): both succeed on the same saved times, the bottom heats by more than 1 K, the temperatures agree within 1e-4 K and the surface temperatures within 1e-3 K, the nodes stay put to 1e-12·L and `thickness_history` stays 0.02 (rtol 1e-12). |
| 2 | Conduction matches the analytic error-function solution within 2 K on both layouts. | covered | Fixed: `test/verification/cases/V-HC-01.jl` QoI "L∞ T-profile error, t=… s" against erfc, with its ConvergenceSpec (minimum order 1.5); moving: `test/verification/cases/V-MMS-03.jl` QoI "L∞ T error on moving mesh, t=… s" (within 0.25 K). |
| 3 | A reacting, shrinking run on the moving mesh carries no spurious or wrong-sign temperature advection. | covered | `test/verification/cases/V-CPL-02.jl` QoI "steady energy balance: v vs Landau v_ss(measured T_s)" (2.1e-3; a spurious mesh-velocity advection breaks it by a few percent); `test/verification/cases/V-MMS-04.jl` QoIs on T, ξ and the surface position under reaction-driven motion. |
| 4 | Conduction on the moving-mesh layout converges at second order (rates above 1.5, mean above 1.8). | covered | `test/verification/cases/V-MMS-03.jl` ConvergenceSpec (`l2_T_final` on 16, 32 and 64 cells, minimum order 1.5, mean 1.8) on a compressing mesh. |

#### test_donor_upwind.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The first-order donor side follows c = −w (forward for w > 0, backward for w < 0, zero at the ghost side). | dropped | The old ALE temperature-advection kernel and its upwind/donor-cell limiters. |
| 2 | A step translated by the old kernel obeys the discrete maximum principle and does not increase the total variation (NONE, VANLEER). | dropped | The old ALE temperature-advection kernel and its limiters. The new relative flux has positivity tests only (Section 4). |
| 3 | Mesh-relative transport moves the field the right way (a front at rest in the material while the mesh slides). | covered | `test/newcore/test_relative_flux.jl` "the relative fluxes converge at second order on a front at rest": the cell masses against the exact integral over each cell's moving nodes; a reversed transport is an O(1) error. |
| 4 | The TVD-limited gradient uses the near/far slope ratio (bounded, no sign flip, smooth and degenerate limits) for minmod, van Leer and superbee. | dropped | The old ALE kernel's limiters. The new van Albada slope has its own test (`test/newcore/test_relative_flux.jl` "van Albada slope"). |

#### test_unified.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The moving-mesh state carries the node positions, and the thickness read from it is L. | covered | `test/newcore/test_states.jl` "state layout: interleaved blocks and row kinds" (position and velocity rows per node, positions differential); `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion" (the initial thickness and its law). |
| 2 | The residual is finite on the position rows. | covered | `test/newcore/test_mesh_motion.jl` "moving mesh: …": residual finite, velocity rows at round-off at the start, eight configurations. |
| 3 | A moving-mesh solve runs without an assembled Jacobian (`jacobian = nothing`). | dropped | The `jacobian` keyword, among the removed keywords (owner question 18). |
| 4 | A moving-mesh solve runs with `Structured(strategy = DirectSolve())`. | dropped | The structured Jacobian. |

#### test_upwind.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The old upwind gradient uses the centre-to-centre distance (error below 10 % away from z = 0 on an alternating h/2h mesh). | dropped | The old ALE temperature-advection kernel. |
| 2 | Mesh-relative transport has no zeroth-order bias on a mesh whose non-uniformity does not vanish under refinement (rate at least one). | ported | `test/newcore/test_relative_flux.jl` "the relative fluxes converge on a front at rest over alternating cell sizes": the front-at-rest error on meshes of alternating cells h, 2h, … (`create_mesh_from_nodes`) at 32, 64 and 128 cells falls at observed orders above 1.5 (1.67 and 1.78). |
| 3 | Uniform-mesh accuracy of the transport. | covered | `test/newcore/test_relative_flux.jl` "the relative fluxes converge at second order on a front at rest" (rate above 1.6 on 32, 64 and 128 uniform cells). |

### test/residual

#### test_ale_dilation.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Every component's concentration row, gas included, carries the dilation term −ξ_j θ_i. | dropped | The ALE dilation term (masses are the states in the new core). |
| 2 | On a sealed, uniformly reacting, collapsing column each cell's inventory of every component changes by its reaction source only. | covered | `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion" (a sealed, uniformly reacting column on a moving mesh: total mass with the gas held at rtol 1e-12); `test/newcore/test_kinetics.jl` "reaction kernel: sources equal the material's rates" (S_mass = ω·V per cell and component); `test/verification/cases/V-BO-04.jl` QoI "solid concentration stays at ρ_s (dilation removes volume only)". |
| 3 | Each cell's volume changes at θ_i V_i. | covered | `test/newcore/test_mesh_motion.jl` "balance rows: a uniform strain rate gives θ x, θ r/2, θ r/3"; `test/newcore/test_correction.jl` "w − c is the solid velocity: θ x, θ r/2, θ r/3 with c ≠ 0" (per cell, the volume rate of w equals θ_k V_k plus that of c at 1e-11). |

#### test_park_equivalence.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The live rows' residual on a masked state equals the compacted problem's residual (fixed Fickian, moving Fickian, moving Beer–Lambert). | ported | `test/newcore/test_deactivation.jl` "masked rows are identity rows and the live block is the compacted problem's: …" (`fixed`, `moving`, `beer_lambert` and the configurations of rows 2, 3 and 7; both forms; one and two dead cells against 5- and 4-cell compacted problems; set-up and perturbed states): with gas in every live cell, written identically into the compact state, each live row's residual equals the compacted problem's within 1e-11 of the magnitude of its terms, abs.(J_c)·abs.(u_c) + abs.(du_c). |
| 2 | The same equivalence with Darcy transport (fixed mesh, two dead cells). | ported | Same testset, `darcy` (fixed mesh, the solid with `κ = 1e-12` and `ρ_intrinsic = 1250`, a `DirichletPressureBC(101325.0)` top) and `darcy_surface_balance` (moving, BEER_LAMBERT), with one and two dead cells: the set-up uses `DarcyFickTransport`, and the residual and Jacobian equalities hold on the live rows. |
| 3 | The same equivalence with the lateral shrinkage law (moving, two dead cells, history rows). | ported | Same testset, `lateral` (moving, `χ -> 1 - 0.2χ`, histories (0.3 + 0.05i)·M and zero in the dead cells): the compact rows mapped through `cell_rows`, `node_rows`, `area_index` and `lateral_rate_index`, with the residual and Jacobian equalities. |
| 4 | Equivalence after interior and bottom parks. | dropped | Parking (interior parking); cells die from the exposed side only. |
| 5 | The dead rows' residual is exactly zero (T, ξ, χ). | ported | Same testset: every dead differential row's residual is 0 and every dead algebraic row's equals its state value, the dead algebraic rows set to 0.5 so that the second check is not 0 == 0. |
| 6 | ALE ghost nodes carry the surface velocity. | dropped | Parking (its ghost nodes). A dead node is frozen at the receiver's far face with zero velocity (`test/newcore/test_deactivation.jl` "a transfer conserves the inventories and the receiver absorbs the volume: …"; Section 4). |
| 7 | The equivalence under a surface-balance top (the surface Newton and its ∂T_s chain on the new exposed face), with Darcy and Beer–Lambert on a moving mesh. | ported | Same testset, `darcy_surface_balance`: the top `RadiativeFluxBC(flux = 2.8e6) + ConvectiveBC(h = 10.0, T_inf = 300.0) + RadiativeBC(T_inf = 300.0)` on a moving mesh under BEER_LAMBERT, with the porous, permeable, gas-seeded solid of row 2, and the residual and Jacobian equalities on the live rows. |

### test/radiation

Every file tests the P1 quasi-steady model, deleted with Experimental (owner question 3). The rows that survive are
the shared parameter helpers, the volumetric absorption coefficient ξ·α_mass and the emissivity mixture.

#### test_p1.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The `RadiationModel` enum's integer codes 0–3, including `P1_QUASI_STEADY = 3`, and the members' type. | dropped | The `P1_QUASI_STEADY` member, deleted with P1 (Experimental, owner question 3). Nothing reads the other members' codes, and every new-core radiation test constructs them. |
| 2 | The radiation diffusion coefficient D_rad = 1/(3(α + σ_s)), its transparent clamp, type stability. | dropped | P1, deleted with Experimental (owner question 3). |
| 3 | Marshak, reflective and Dirichlet radiation conditions: construction, defaults, validation, `default_radiation_bcs`. | dropped | P1's radiation conditions (owner question 3). |
| 4 | `evaluate_param` (a Real or a function of t) and `evaluate_param_T` (a Real or a function of T) return the right value with an inferred type. | ported | `test/materials/test_property_validation.jl` "parameter evaluation: numbers, functions and property laws": `evaluate_param` on a Real and a function of t; `evaluate_param_T` on a Real, a function of T, a `LinearProperty` and a `ConstantProperty`; a Dual argument giving a Dual with zero partials for both Real methods and a Dual through a function of t; `@inferred` on every method. |
| 5 | The Marshak flux-derivative coefficient c = ε/(2(2−ε)) against central differences, with a T-dependent emissivity. | dropped | P1's Marshak closure (owner question 3). |
| 6 | The half-cell effective Marshak coefficient c_eff = 1/(1/c + δ/D) and its limits. | dropped | P1's Marshak closure (owner question 3). |
| 7 | `RadiationCache`: construction, resize, `initialize_radiation_intensity!`, energy density, net source. | dropped | P1 storage (owner question 3); deleted with the old core, `RadiationCache` has no caller outside it. |
| 8 | The Marshak condition's flux q = c(G_ext + 4I_ext − G), the half-range weight 4c, zero at equilibrium. | dropped | P1's Marshak condition (owner question 3). |
| 9 | The reflective condition's flux is zero. | dropped | P1's reflective condition (owner question 3). |
| 10 | The volumetric absorption coefficient is ξ·α_mass. | ported | `test/physics/test_radiation.jl` "Beer–Lambert surface absorptivity" > "the absorption coefficient is Σ ξ_k α_k of each cell's composition": `effective_absorption` for one component and for two solids plus a gas, and after `setup(…; radiation_model = BEER_LAMBERT)` and one `residual!` the `ds.α` of a graded state equals the literal sums in all 8 cells, all distinct. |
| 11 | The cell emissivity is the component's. | covered | The volume-weighted sum: `test/materials/test_emissivity_absorptivity.jl` "ε_eff = Σ v_j ε_j(T) at the cell temperature" (rtol 1e-14). The per-cell value: "the derived state fills ε and a per cell at the cell's composition and temperature" in the same file (`ds.ε[i] ≈ mixture_emissivity(ξ_i, T_i, mat)` per cell, rtol 1e-14; the rewrite of that file's row 6). |
| 12 | The P1 sources Q_emission = 4σαT⁴ and Q_absorption = αG. | dropped | P1's volumetric emission and absorption (owner question 3). |
| 13 | The radiation property update: D_rad per cell and at faces. | dropped | P1 (owner question 3). The absorption coefficient it also reads is row 10. |
| 14 | The sign of the P1 flux, and a zero flux at a reflective face. | dropped | P1 (owner question 3). |
| 15 | Optical thickness τ = αL, per-cell τ, and the `check_p1_validity` regimes. | dropped | P1's diagnostics (owner question 3). The absorption coefficient they also read is row 10. |
| 16 | The keys and signs of `radiation_energy_balance`. | dropped | P1's diagnostic (owner question 3). |

#### test_p1_experimental.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | G is not a state variable: the state has no G block, and G is re-solved at every evaluation. | dropped | `P1Residual`, deleted with Experimental (owner question 3). |
| 2 | The constructor requires the P1 radiation model (a guard against counting the external flux twice). | dropped | `P1Residual` (owner question 3). |
| 3 | `RadiativeFluxBC + RadiativeBC` contribute nothing at the boundary under P1. | dropped | P1's boundary handling (owner question 3). The new core refuses P1 at setup. |
| 4 | The P1 source enters the energy rows: a hot Marshak field heats the surface cell faster. | dropped | `P1Residual` (owner question 3). |

#### test_p1_sign.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The P1 source S_rad = κ(G − 4σT⁴) is positive in a hot field, negative in a cold one and zero at equilibrium. | dropped | P1's flux divergence, deleted with Experimental (owner question 3). |

#### test_p1_validation.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The Rosseland-limit equilibrium: G starts at 4σT⁴ and the hot-to-cold G ratio is 16. | dropped | P1 (owner question 3). The absorption coefficient it also reads is `test_p1.jl` row 10. |
| 2 | The P1 source cools hot regions, and an external flux heats the surface cell. | dropped | P1 (owner question 3). |
| 3 | The quasi-steady integral balance closes, and in a closed reflective system G stays finite and positive. | dropped | P1's quasi-steady solve (owner question 3). |
| 4 | The cold absorber's absorbed fraction matches continuum P1: q_abs/I_ext = 4c/(1 + √3 c). | dropped | P1's Marshak closure (owner question 3). |
| 5 | The cell absorption coefficient is ξ·α_mass = 10000. | ported | The test of `test_p1.jl` row 10. |
| 6 | Optical-thickness diagnostics: τ scaling, per-cell τ, validity. | dropped | P1's diagnostics (owner question 3). The absorption coefficient they also read is `test_p1.jl` row 10. |
| 7 | Reflective, zero-flux Marshak and Dirichlet-G radiation conditions. | dropped | P1's conditions (owner question 3). |
| 8 | The quasi-steady solver's G is finite and positive for four temperature profiles. | dropped | P1's Newton solve (owner question 3). |
| 9 | The P1 residual differs from the no-radiation residual (S_rad is non-zero). | dropped | P1 (owner question 3). |

### test/verification

Only `test_analytical.jl` is deleted; the verification cases stay.

#### test_analytical.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | A Dirichlet-wall transient matches erf, max error below 2 K (100 cells). | covered | `test/verification/cases/V-CPL-01.jl` QoI "T-history L∞ error at x = … mm" (a Dirichlet top wall, tolerances 0.05 to 0.4 K); `test/verification/cases/V-HC-01.jl` QoI "L∞ T-profile error, t=… s" (within 0.75 K at 100 cells); `test/newcore/test_boundary_conditions.jl` "a prescribed face temperature is the reported surface temperature". |
| 2 | Transient Dirichlet conduction converges at order 2 (every rate above 1.5, mean above 1.8). | covered | Transient order: the ConvergenceSpec of `test/verification/cases/V-HC-01.jl` (40, 80, 160 cells, L2 at 60 s); Dirichlet faces: the ConvergenceSpecs of `V-HC-05.jl` and `V-CPL-04.jl` (both faces Dirichlet), all at 1.5 and 1.8; full runs only (Section 4). |
| 3 | Beer–Lambert: each cell absorbs I_up(1 − e^(−αΔz)), within 1 %. | covered | The attenuation: `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity" (the per-cell source I_up(1 − e^(−τ)) to rtol 1e-12, on both absorptivity paths). |
| 4 | The cell absorption coefficient ρ·α_mass that the expected source used. | ported | `test/physics/test_radiation.jl` "Beer–Lambert surface absorptivity" > "the absorption coefficient is Σ ξ_k α_k of each cell's composition" (the test of `test/radiation/test_p1.jl` row 10). |
| 5 | Beer–Lambert energy conservation: absorbed I₀(1 − e^(−αL)) within 5 %, transmitted I₀e^(−αL) within 1 %. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": the transmitted intensity and absorbed plus transmitted against the entering beam, rtol 1e-12; `test/verification/cases/V-RAD-01.jl` couples absorption with conduction. |
| 6 | First-order Arrhenius decay ξ₀e^(−kt) in a Dirichlet-pinned cell, within 1 %. | covered | `test/verification/cases/V-RX-01.jl` QoI "max relative ξ_… error" (within 1e-8); `V-RX-02.jl` variant c (live Arrhenius at a pinned temperature, within 1e-9); `test/newcore/test_solve.jl` "solve: kinetics through the adapter, both energy forms". |
| 7 | Mass conservation in the reaction: ξ_A + ξ_B = ξ₀ to 1e-8. | covered | `test/verification/cases/V-RX-01.jl` QoI "total-mass conservation (closed chain)" and `V-RX-02.jl` QoI "total-mass conservation", within 1e-13. |
| 8 | TGA at a constant heating rate: conversion matches the quadrature solution (0.02 between 5 % and 95 %, 0.05 overall), reaching significant conversion. | covered | `test/verification/cases/V-CPL-05.jl` QoIs "thin-limit deviation from 0D TGA (L = … mm)" (within 3e-3 of the exact curve) and "complete conversion at t_end (coarsest rung)" (within 1e-6). |
| 9 | The mass-loss-rate peak lies within 2 K of the exact first-order peak temperature (βE/(RT_p²) = A e^(−E/RT_p)). | ported | `test/newcore/test_kinetics.jl` "constant heating rate: the mass-loss rate peaks at the first-order peak temperature": a 0.1 mm one-cell column between faces at `300 + βt` (β = 10 K/min, A = 1e13, E = 180 kJ/mol, no reaction heat) solves to `:Success`, and the temperature at the middle of the saved interval of the largest mass-loss rate lies within 2 K of the T_p solving βE/(RT_p²) = A e^(−E/RT_p) (0.2 K). |
| 10 | Over β = 5, 10, 20, 40 K/min the Kissinger slope recovers E within 5 %. | ported | `test/newcore/test_kinetics.jl` "constant heating rates: the Kissinger slope recovers the activation energy": on row 9's deck at β = 5, 10, 20 and 40 K/min the least-squares slope of ln(β/T_p²) against 1/T_p gives E within 5 % (1.2 %). |
| 11 | The decaying sine mode with Dirichlet faces converges spatially (rates above 1.0, mean above 1.5). | covered | The ConvergenceSpec of `test/verification/cases/V-HC-02.jl` (Bi = 10 variant: a transient eigen-series, 40, 80, 160 cells, 1.5 and 1.8); `V-HC-05.jl` and `V-CPL-04.jl` for the Dirichlet faces; full runs only (Section 4). |
| 12 | An inert material keeps `extras.total_mass` constant under a heat flux. | covered | `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion" (`total_mass` constant, rtol 1e-12); on a fixed mesh the solid mass rows are exactly zero (`test/newcore/test_boundary_conditions.jl` "an ambient gas enters a gas-free column through a film condition"). |
| 13 | The reported `total_mass` is the absolute per-area inventory ξ₀·L. | ported | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (inert, fixed mesh, a `HeatFluxBC` top, adiabatic bottom): `sol.extras.total_mass ≈ fill(800.0 * L, 3)`, rtol 1e-12. |
| 14 | A closed system of 50 Dirichlet-pinned reacting cells conserves A + B to 1e-8. | covered | `test/newcore/test_residual_jacobian.jl` "row assembly conserves mass and scatters fluxes antisymmetrically" (the mass rates of a closed reacting slab sum to within 1e-12 of their scale); `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion" (total mass over a solve, rtol 1e-12); V-RX-01, V-RX-02, V-RX-03 total-mass QoIs within 1e-13. |
| 15 | An adiabatic closed slab conserves ∫ρcT dV and its mean temperature to 1e-6. | covered | `test/verification/cases/V-HC-03.jl` QoI "mean-temperature rise (global energy)" (within 1e-5 K of q″Δt/(ρcL)); `test/newcore/test_solve.jl` "energy forms: invariant to round-off, temperature form at fourth order" (a closed adiabatic pair conserves enthalpy to 1e-12 in the enthalpy form). |
| 16 | In an adiabatic exothermic slab ΔE_thermal = −h·Δm_A within 5 %. | covered | `test/verification/cases/V-RX-03.jl` QoI "reaction-heat superposition (quasi-steady cell drift)" (within 1e-4 K, pinning the reaction heat's magnitude); `test/newcore/test_residual_jacobian.jl` "row assembly conserves mass and scatters fluxes antisymmetrically" (Σ dH/dt equals the reaction heat plus the formation terms, rtol 1e-12); `test/newcore/test_energy_ledger.jl` "energy ledger over a solve: …". |
| 17 | Under ALE the initial mass equals the final mass plus ∫ gas generation dt, within 10 %. | covered | `test/verification/cases/V-CPL-02.jl` QoI "mass ledger closure over plateau window" (within 4.9e-3); `test/newcore/test_conservation.jl` "the recorder across deactivation events: …" (moving: relative defect within 1e-7). |
| 18 | Under ALE with depletion `thickness_history` starts at L, decreases monotonically, ends below L, and has one value per saved time. | covered | `test/newcore/test_deactivation.jl` "a solve through deactivations: …" (moving: non-increasing within 1e-12·L, ending at the live span); `test/newcore/test_mesh_motion.jl` "solve: the swelling column's thickness follows its conversion" (starting at L₀, the exact L(t)); every extras series has one value per saved time by construction. |
| 19 | Under ALE an inert column keeps its thickness L and its mass (1e-6). | covered | `test/newcore/test_correction.jl` "isothermal advection and free stream under the equidistribution correction: …": an inert two-solid ALE column keeps its top node and Σm to rtol 1e-12. It is isothermal, but the velocity rows read only the reaction mass sources, so heat cannot move an inert column. |
| 20 | Darcy–Fick transport leaves the solid rows untouched, and sealed faces carry no gas. | covered | `test/newcore/test_gas_transport.jl` "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out" (zero solid fluxes, zero flux through the sealed back face); `test/verification/cases/V-GT-01.jl` QoI "gas-mass conservation (sealed)" (Fickian, within 1e-13). |
| 21 | The gas mass rates of a sealed Darcy–Fick column (κ > 0, a pressure gradient) sum to zero. | ported | `test/newcore/test_gas_transport.jl` "gas transport: a sealed Darcy–Fick column keeps its gas and its solid (… form)" (both forms): `DarcyFickTransport`, every mobility positive, the largest pressure more than 10× the smallest, non-zero face fluxes, exactly zero flux through both sealed faces, zero solid rows, and the gas rows summing to zero within 1e-12 of their absolute sum; a solve succeeds with the gas inventory constant to 1e-12. |
| 22 | Under a prescribed flux T_s = T₀ + (2q/k)√(αt/π), within 2 %. | covered | `test/verification/cases/V-HC-01.jl` QoI "surface temperature, t=… s" (within 1.5 K through the surface balance, flux plus convection); `test/verification/cases/V-HC-03.jl` QoI "L∞ T-profile error, Fo=…" (a pure `HeatFluxBC`, within 0.006 K); `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (under a pure flux T_s = T_cell + q d/k, rtol 1e-12). |
| 23 | The top boundary balance (applied against conducted flux) closes to 1e-6. | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances": the conducted flux equals q (rtol 1e-12) at every saved time. |
| 24 | T_s at 1, 5, 10 and 20 s within 3 %, increasing. | covered | `test/verification/cases/V-HC-01.jl` QoI "surface temperature, t=… s" at 30, 60, 90 and 180 s, within 1.5 K of the exact value; the exact increments are far larger than 1.5 K, so the result is monotone. |
| 25 | With `HeatFluxBC + RadiativeBC` on top the top imbalance stays below 1e-4 and the adiabatic bottom's below 1e-10. | covered | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate" (net ≈ conduction, rtol 1e-8, with the decomposition summing to 1e-12, at every saved time) and "outputs: face temperatures from the boundary balances" (adiabatic bottom: `boundary_T[:bottom] == T[1, :]`, zero back flux). |

### test

Only `helpers.jl` at the top level. Every file that includes it is deleted, and
`test/geometry/test_variable_cross_section.jl` builds its problem inline, so the file goes whole (Section 3.1).

#### helpers.jl

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `fixture_material(; lateral, darcy)` gives a solid-to-gas material with an optional lateral law and solid permeability, and `fixture_problem(…)` puts it on a uniform 1 cm mesh with a 10 kW/m² flux top, an adiabatic bottom and a Fickian top mass condition. | dropped | No surviving consumer: every file that includes `helpers.jl` or calls a `fixture_*` helper is deleted, and `test/geometry/test_variable_cross_section.jl` builds its problem inline. The file goes whole (Section 3.1). |
| 2 | `fixture_state(problem; use_ale, radiation_model, radiation_intensity)` returns `(def, ws, u, du)` from `initialize_state`. | dropped | As row 1; `fixture_state` also reads `initialize_state`, deleted with the old core. |

### test/boundary_conditions

#### test_blowing_correction.jl

Old-only testsets, deleted whole (Section 3.1): "Blowing correction" > "evaluate_bc_flux_and_derivative: face
dispatch" > "evaluate_bc_flux_and_derivative(::ConvectiveBC) without a face" and "Blowing correction" >
"Structured Jacobian (verify_jacobian_at): strict agreement with blowing". Outside them the file uses no old-core
name: the surface Newton runs in place with `trial = trial_temperature`, the convective reference is a closed form,
and the solves are `Pyrolysis.Integration.solve` on the native mesh.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `blowing_factor` values and limits, suction enhancement, monotone and positive. | kept | "Blowing correction" > "blowing_factor: values and limits". |
| 2 | `blowing_factor` is C¹ across the Taylor and exact seams, with f′(0) = −1/2. | kept | "Blowing correction" > "blowing_factor: branch seams (value + derivative)". |
| 3 | `blowing_factor` and its ForwardDiff partials stay finite at large ±B. | kept | "Blowing correction" > "blowing_factor: large-B overflow safety". |
| 4 | `evaluate_bc_flux_and_derivative` face dispatch: only convective members rescale, zero outflow reproduces the face-less flux bit for bit, user rules apply. | kept | "Blowing correction" > "evaluate_bc_flux_and_derivative: face dispatch": the face-less reference is the closed form h(T_∞ − T_s), reproduced bit for bit at zero outflow. |
| 5 | `CouetteFilm` λ semantics, validation, bit-identity of the default, and the decomposition applying the same factor. | kept | "Blowing correction" > "CouetteFilm: film constant". |
| 6 | The surface Newton with blowing: a converged balance, the `blowing_mcp` and `blowing_ratio` diagnostics, the on, off, zero and suction ordering. | kept | "Blowing correction" > "surface Newton with blowing: converged balance" (`solve_boundary_temperature!` with `trial = trial_temperature`). |
| 7 | ∂T_s/∂T_cell and ∂T_s/∂ṁ″ by ForwardDiff match central differences. | kept | "Blowing correction" > "AD: T_s derivatives match finite differences" (with `trial`). |
| 8 | The structured backend's Jacobian equals dense AD under the strict criterion with a blowing model on the top face. | dropped | The structured Jacobian (`verify_jacobian_at`). Its assertion sits in the old-only "Blowing correction" > "Structured Jacobian (verify_jacobian_at): strict agreement with blowing" (Section 3.1). |
| 9 | The coupling row 8 protected, under `CouetteFilm()`: Fickian and Darcy–Fick on a fixed mesh, Fickian on a moving one. | covered | `test/newcore/test_gas_transport.jl` "gas transport: …" (`fick_diffusive`, `darcy_venting`, `darcy_combined` carry `ConvectiveBC(…, blowing = CouetteFilm())` with outflow): containment and per-entry equality at rtol 1e-12, five states, both forms; `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`fick`, a `CouetteFilm` top). |
| 10 | The same under `CouetteFilm(0.31)`, `PropertyRatioFilm(NITROGEN)` and `PoolingFilm(NITROGEN; …)` (fixed Fick, fixed Darcy–Fick, moving Fick). | ported | `test/newcore/test_gas_transport.jl` "gas transport: …" (`fick_couette_film`, `fick_property_ratio`, `fick_pooling`, `darcy_couette_film`, `darcy_property_ratio`, `darcy_pooling`: `CouetteFilm(0.31)`, `PropertyRatioFilm(NITROGEN)`, `PoolingFilm(NITROGEN; diffusion_volume = (gas = 60.0,))`): containment, a finite residual and Jacobian, colored equal to dense at rtol 1e-12 at five states, both forms, 6 and 8 cells, and a residual that allocates nothing. The top is `ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = m) + RadiativeBC(T_inf = 300.0)` over a `DiffusiveMassBC(0.02)` (Fick) or a gentle `ConvectivePressureBC(1e-8, 101325.0)` (Darcy–Fick): with T_inf = 1000 K the face stays below the free stream, so the pooling root solve never runs, and a `DirichletPressureBC` vent drives the blowing factor to 0 (C129). "gas transport: … blows through a face hotter than its free stream" asserts at `u0` that gas leaves, T_face > T_inf, 0 < χ < 1 and, for pooling, ψ < 1. On a moving mesh: `test/newcore/test_mesh_motion.jl` "moving mesh: pooling, …" (`HeatFluxBC(3e4) + ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = PoolingFilm(…)) + RadiativeBC(T_inf = 300.0)`), with ψ < 1 at `u0`, `hot` and `perturbed`. |
| 11 | A solve with `CouetteFilm` ends hotter and loses mass faster than without; the ratio is identically 1 off and below 1 on. | covered | `test/newcore/test_boundary_conditions.jl` "blowing lowers the convective loss of an outgassing surface": the same deck shape, a hotter surface and a larger surface mass flux with blowing, the ratio 1 off and below 1 on. Also asserted on the new solve in "Blowing correction" > "integration: solve() with blowing". |
| 12 | The saved decomposition stays consistent under blowing: net ≈ conduction on and off, and the convective loss is positive and below h(T_s − T_∞). | kept | "Blowing correction" > "integration: solve() with blowing" (`Pyrolysis.Integration.solve`): net ≈ conduction on and off (1e-6), and the convective loss in (0, h(T_s − T_∞)). |
| 13 | `PoolingFilm` runs through the integrator: `:Success`, ψ in [0.59, 1], the face fraction in [0, 1], a ratio below 1. | kept | Same testset: `PoolingFilm` through the integrator, `:Success`, ψ in [0.59, 1], the face fraction in [0, 1], a ratio below 1. |

#### test_blowing_models.jl

Old-core names: none; every `solve_boundary_temperature!` call passes `trial`.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Gas data: NASA-7, conductivities, Eucken, Fuller (values, P and T scaling), `GasSpecies` validation and the range clamp. | kept | "Blowing models" > "gas data". |
| 2 | `PropertyRatioFilm` limits, the Couette reduction at n_Γ = 0, the calculator's reference ratios, monotone in the outflow. | kept | "Blowing models" > "PropertyRatioFilm". |
| 3 | `PoolingFilm` reference ratios, ψ, Y, Λ, χ = ψ f(φ/ψ), the limits and the pooling floor. | kept | "Blowing models" > "PoolingFilm: reference values". |
| 4 | AD through the pooling face solve and through the surface Newton matches central differences on the switch. | kept | "Blowing models" > "PoolingFilm: AD through the face solve" (the surface Newton with `trial`). |
| 5 | The per-species mixture path: one component equals the lumped description, the purge is excluded from the face, Wilke bounds, inflow excluded. | kept | "Blowing models" > "per-species mixture path". |

#### test_boundary_temperature.jl

Old-only testsets, deleted whole (Section 3.1): "surface Newton: a_eff absorption vs. two-piece split" >
"effective_absorptivity on the line" and > "compute_boundary_heat_fluxes(::CombinedThermalBC) without a face";
"surface Newton ignores the stored state" > "solve_boundary_temperature! without a trial clamp: a NaN cell warns";
"P1_QUASI_STEADY: RadiativeFluxBC contributes nothing at the face"; "evaluate_bc_flux_and_derivative without a
face: combined and convective"; "Discretization accessors: get_surface_temperature,
check_boundary_energy_balance". Outside them the surface Newton runs in place with `trial = trial_temperature`, the
combined and convective conditions go through the face methods with a zero-outflow `BlowingFace`, the mixtures are
`mixture_emissivity` and `mixture_absorptivity`, and the solves are `Pyrolysis.Integration.solve` on the native
mesh.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | A prescribed flux (a function): T_s = T_cell + q d/k within 5 iterations, with the stored inputs. | kept | "solve_boundary_temperature!" > "Pure Heat Flux BC". |
| 2 | A constant prescribed flux: applied = conducted = q. | kept | "solve_boundary_temperature!" > "Constant Heat Flux BC". |
| 3 | A convective face: the closed-form T_s, between T_∞ and T_cell. | kept | "solve_boundary_temperature!" > "Convective BC". |
| 4 | A radiative face: the balance and εσ(T_∞⁴ − T_s⁴) at the solved T_s. | kept | "solve_boundary_temperature!" > "Radiative BC". |
| 5 | The combined cone condition converges within 20 iterations to a balance. | kept | "solve_boundary_temperature!" > "Combined BC (flux + radiative + convective)". |
| 6 | An adiabatic face: one iteration, T_s = T_cell. | kept | "solve_boundary_temperature!" > "Adiabatic BC". |
| 7 | A Dirichlet face returns the prescribed temperature and its conducted flux. | kept | "solve_boundary_temperature!" > "Dirichlet BC". |
| 8 | A time-dependent flux is sampled at t. | kept | "solve_boundary_temperature!" > "Time-dependent BC". |
| 9 | An outward gas flux changes nothing without a convective member: the balance carries no gas-enthalpy term. | kept | "solve_boundary_temperature!" > "Outward gas flux enters only through blowing"; its comment points to the gas-transport kernel's enthalpy flux. |
| 10 | Convergence robustness: a hot radiative face, very thin and very thick half-cells. | kept | "Convergence Robustness" (three nested testsets). |
| 11 | `RadiativeFluxBC` under NO_RADIATION and SURFACE_ABSORPTION absorbs ε_eff·q with zero slope; under BEER_LAMBERT it contributes nothing at the face. | kept | "RadiativeFluxBC with Radiation Models" > "NO_RADIATION: flux applied at boundary", "SURFACE_ABSORPTION: flux applied at boundary" and "BEER_LAMBERT: flux handled volumetrically". |
| 12 | `RadiativeFluxBC` under P1_QUASI_STEADY contributes nothing at the face. | dropped | P1 radiation, deleted with Experimental (owner question 3). Its assertion sits in the old-only "P1_QUASI_STEADY: RadiativeFluxBC contributes nothing at the face" (Section 3.1). |
| 13 | An explicit absorptivity overrides ε_eff under the surface models and is zero under BEER_LAMBERT. | kept | "RadiativeFluxBC with Radiation Models" > "RadiativeFluxBC with explicit absorptivity", without the P1 lines. |
| 14 | A combined `RadiativeFluxBC + ConvectiveBC`: sums and slope per model. | kept | "RadiativeFluxBC with Radiation Models" > "Combined BC with RadiativeFluxBC", through the face method. |
| 15 | The `a_eff` keyword: absorbs at a_eff and emits at ε_eff, an explicit α ignores it, combined and face-aware forms forward it, volumetric models give zero. | kept | "RadiativeFluxBC with Radiation Models" > "absorptivity = nothing absorbs at a_eff, emits at ε_eff", the combined condition through the face method. |
| 16 | `get_surface_temperature` and `get_back_face_temperature` read the mesh's stored boundary states. | dropped | The old `Mesh1D`'s boundary states and the mesh-bound wrappers of the surface Newton in `src/Discretization`. Their successors `sol.surface_T` and `sol.boundary_T` are asserted in `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances". The assertions sit in the old-only "Discretization accessors: get_surface_temperature, check_boundary_energy_balance" (Section 3.1). |
| 17 | `check_boundary_energy_balance` computes the relative imbalance of the stored states. | dropped | The mesh-bound wrappers in `src/Discretization`. The balance itself is row 19. The assertion sits in the same old-only testset as row 16. |
| 18 | The flux and slope of `HeatFluxBC`, `ConvectiveBC`, `RadiativeBC` and `AdiabaticBC`. | kept | "BC Flux and Derivative Evaluation", the `ConvectiveBC` case through the face method. |
| 19 | During a solve the surface temperature satisfies the face balance at both faces and exceeds the cell under heating. | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances" (T_s = T_cell + q d/k, net = q, conduction ≈ q to 1e-12, an adiabatic back at the cell temperature with zero flux) and "outputs: the surface-flux decomposition and the mass-loss rate" (a nonlinear radiative plus convective top: net ≈ conduction to 1e-8). Also asserted on the new solve in "Integration with Solver" > "the face balances hold at both faces during a solve, the heated one converged". |
| 20 | After the solve the stored states carry `converged = true` at both faces. | covered | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate" (an independent decomposition at the solved T_s equals the conducted flux, rtol 1e-8) and "outputs: face temperatures from the boundary balances" (the adiabatic back). The flag: "Integration with Solver" > "the face balances hold at both faces during a solve, the heated one converged" reads `ThermalBoundaryKernel.states[top].converged`; an adiabatic face keeps `converged = false` (Section 4). |
| 21 | `solution.surface_T` is the solved face temperature, not the top cell's, and differs from it under 75 kW/m². | covered | `test/newcore/test_solve.jl` "outputs: face temperatures from the boundary balances": `surface_T == boundary_T[:top]` and ≈ T_cell + q d/k. Also asserted on the new solve in "Integration with Solver" > "Solution output uses boundary temperature". |
| 22 | Per-component absorptivity on a line a_j = a + bε_j: the single condition equals the two-piece split; the decomposition is absorbed(a) + emitted(ε) + convective; the Kirchhoff default differs; the diagnostic decomposition applies a_eff. | kept | "surface Newton: a_eff absorption vs. two-piece split": the surface Newton in place with `trial`, `mixture_emissivity` and `mixture_absorptivity`, and the face decomposition. |
| 23 | The surface Newton starts from the cell temperature, never the stored state, and stores only finite iterates. | kept | "surface Newton ignores the stored state", on the trial path: a NaN cell gives a NaN face, not stored, without a warning. |

### test/materials

#### test_emissivity_absorptivity.jl

Old-only testsets, deleted whole (Section 3.1): "effective_absorptivity: Kirchhoff invariant and per-component
override" and "PropertyCache and update_properties! fill ε and a per cell". Outside them the file runs on `setup`,
`residual!` and `buffers` and on the `mixture_*` sums; each `effective_emissivity` assertion (Section 4) has a
`mixture_emissivity` assertion of the same property beside it.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | A component's ε takes every property-function form, with defaults, and Duals propagate. | kept | "component ε accepts every property-function form". |
| 2 | Surface absorptivity: construction, the static `Nothing`/`Float64` parameter, the Kirchhoff fallback, validation. | kept | "surface absorptivity: construction, validation, accessors". |
| 3 | ε_eff = Σ v_j ε_j(T) at the cell temperature, the constant case, the AD slope. | kept | "ε_eff = Σ v_j ε_j(T) at the cell temperature": each `effective_emissivity` assertion (two temperatures, the T dependence, the constant case, the AD slope) has a `mixture_emissivity` assertion of the same property beside it. |
| 4 | The clamp still applies: an unphysical ε = 1.3 law gives ε_eff == 1. | dropped | The [0, 1] clamp on emissivity, one of the hard-clamped property variants; the new core reads the unclamped sums. Its assertion sits in "effective_emissivity clamps an unphysical law at 1", deleted with the clamp (Section 4). |
| 5 | a_eff ≡ ε_eff bit for bit when no component sets one; an override replaces ε_j(T) in the same sum; slopes; a depleted cell gives 0. | kept | "a_eff: Kirchhoff invariant and per-component override", on `mixture_absorptivity` and `mixture_emissivity`, the bitwise identity included. |
| 6 | The per-cell ε and a the solver uses are the property functions at each cell's (ξ, T): a == ε for a Kirchhoff material, distinct and each correct with one override, varying with T. | kept | "the derived state fills ε and a per cell at the cell's composition and temperature": `setup`, masses and temperatures set per cell, `residual!` and `buffers`; `ds.a == ds.ε` for a Kirchhoff material, `ds.ε[i]` and `ds.a[i]` against the `mixture_*` sums (1e-14), distinct values, and the override case. |

#### test_gamma_semantics.jl

Old-core names: none. `effective_emissivity` (Section 4) is asserted beside `mixture_emissivity`; no expected value
changes.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | One γ semantics: the density, the mixing weights and θ agree for γ = 0, 1 and partial. | kept | "one γ semantics across θ, ρ_mix, and mixing weights". |
| 2 | γ is validated to [0, 1] at construction for each phase. | kept | Same testset. |
| 3 | A fully converted cell falls back to the gas conductivity, and a pure-gas cell's ε is 0. | kept | "depleted-cell floor and liquid ε default", with `mixture_emissivity` beside `effective_emissivity` (0.90, and 0 for a pure-gas cell). |
| 4 | Liquids default to ε = 0.95 and raise ε_eff of a wet surface. | kept | Same testset (`mixture_emissivity` 0.9·0.9 + 0.1·0.95 for the wet surface). |

#### test_gas_viscosity_property.jl

Old-only testset, deleted whole (Section 3.1): "effective_gas_viscosity and mixture_gas_viscosity" (the two
accessors and their hard floor). Outside it the file reads `material.gas_viscosity`, the new core's mobility
`ds.mob` and `Pyrolysis.Outputs.get_cell_positions`, and takes the pore fraction as 1 − `skeletal_volume_fraction`.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `SutherlandProperty()` is Sutherland's law with the documented anchors. | kept | "Material gas viscosity" > "SutherlandProperty is Sutherland's law". |
| 2 | The property-function contract: the 2-argument form, AD-clean. | kept | "Material gas viscosity" > "property-function contract". |
| 3 | `Material` defaults to air bit for bit. | kept | "Material gas viscosity" > "Material defaults to air", on `gas_viscosity`. |
| 4 | The material's curve is consulted (a constant, a number, a function). | kept | "Material gas viscosity" > "material curve is actually consulted". |
| 5 | The surrogate's quartic fit is accepted and lies below air. | kept | "Material gas viscosity" > "polynomial curve (the surrogate's native form)". |
| 6 | Guards: a low-temperature evaluation is finite, a curve going negative is floored at `MU_FLOOR`, and the floor never binds physically. | kept | "Material gas viscosity" > "guards: μ only ever divides": the −50 K equality through `T_PROPERTY_MIN`; the new core's mobility `ds.mob ≈ κ/μ` for the default curve at 650 K and for a falling user fit at 1000 K and at 3000 K, where the curve is negative and `MU_FLOOR` binds through the smooth floor. |
| 7 | End to end: doubling μ doubles the steady Darcy overpressure², each solve matching V-GT-02's exact profile. | kept | "Material gas viscosity" > "end-to-end: doubling μ doubles the Darcy overpressure²": `Pyrolysis.Integration.solve` on the native mesh, `Pyrolysis.Outputs.get_cell_positions`, φ = 1 − `skeletal_volume_fraction`; the profile errors (1.3e-3 and 2.6e-3) within 5e-3 and the ratio within 1e-2 (4.0e-3). It is the test that `gas_viscosity` reaches the new core's mobility. |
| 8 | The constructor rejects a non-positive or NaN curve. | kept | "Material gas viscosity" > "constructor rejects an unusable curve". |

#### test_kirchhoff.jl

Old-core names: none.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The Kirchhoff-lifted heats equal the DSC heats at the DSC peaks, and untouched stages are unchanged. | kept | "White pine Kirchhoff heat lift". |
| 2 | The lifted heats lie below DSC off-peak, finite with a negative ForwardDiff slope. | kept | Same testset. |

#### test_property_validation.jl

Old-core names: none.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | Laws evaluate on the property bracket: the identity inside, clamped outside with zero slope, for every law type. | kept | "property-law bracket" > "argument clamp is the identity inside the bracket". |
| 2 | A zero-crossing conductivity is rejected, naming the temperature and the opt-out. | kept | "property-law bracket" > "zero-crossing conductivity is rejected with the crossing temperature". |
| 3 | Non-positive ρ, c, k and out-of-range ε are rejected for every phase; good laws build. | kept | "property-law bracket" > "other laws and component types". |
| 4 | State-dependent laws are not sampled. | kept | "property-law bracket" > "state-dependent laws are not sampled". |
| 5 | `identically_zero` decides by structure, never by samples. | kept | "a law vanishes by its structure, never by samples". |

### test/physics

#### test_background_air.jl

Old-core names: none. The fields are `initial_fields(problem, problem.mesh)` at
`Pyrolysis.Problem.cell_centers(problem.mesh)` on the native mesh, and the pore fraction is
1 − `skeletal_volume_fraction`. `background_gas_concentration` reads the clamped `porosity` internally, with the
same values either way (Section 4).

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `background_gas_concentration` inverts the pressure closure at ambient: a hand number, exactly 1 atm back, the index form equal to the name form, linear in P. | kept | "Background pore air" > "EOS round-trip at ambient": φ = 1 − `skeletal_volume_fraction`, the same hand numbers. |
| 2 | A uniform composition with air at 1 atm gives a flat 1 atm pressure field, with no driving gradient. | kept | "Background pore air" > "equilibrium hold: uniform air ⇒ no Darcy driving force": the closure on `initial_fields` at the cell centres of the native mesh, and the core's own `ds.P` after `setup` and one `residual!` (`DarcyFickTransport`, every cell at 101325 Pa to 1e-12). |
| 3 | A non-gas component or an unknown name is refused. | kept | "Background pore air" > "argument validation". |
| 4 | Without `ρ_intrinsic` there is no pore volume, and the helper returns exactly 0. | kept | "Background pore air" > "no-pore-space degeneracy is visible, not silent". |

#### test_kinetics.jl

Old-only testset, deleted whole (Section 3.1): "TanhLimiter". The other testsets run on `DepletionLimiter`. Four
groups of expected values differ from the tanh limiter's, because below the join the onset is the one-sided cubic
g(x) = (s − 2y)x³ + (3y − s)x² (y = tanh 1, s = sech² 1) instead of tanh x, and the power reads the finite base
b(ξ) = ξ + δ log(1 + e^(−ξ/δ)) instead of max(ξ, 0).

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The limited rate is continuous with a continuous derivative across the threshold, and has the stated closed form. | kept | "depletion limiter" > "continuous and smooth across the threshold": below the join the closed form k ξ g(ξ/ξ_th) to 1e-12, and the C¹ check at the join. |
| 2 | Far above the threshold the bare Arrhenius rate returns; the saturated path is bit-identical to the unlimited rate with slope k. | kept | "depletion limiter" > "bare Arrhenius rate recovered far above threshold" (values unchanged). |
| 3 | Order-0 reactions stop at depletion: zero at 0 and below, a finite slope at 0⁺, the unlimited rate k at ξ = 0. | kept | "depletion limiter" > "zero-order reactions stop at depletion": zero at 0 and below, the unlimited rate k, r₀(1e-3) = k g(1e-3), and a slope at 1e-12 of k g′(1e-12)/ξ_th, vanishing at depletion. |
| 4 | For n < 0.5 the slope stays finite and falls toward depletion. | kept | "depletion limiter" > "bounded ∂r/∂ξ at depletion for n < 0.5": the slope at 1e-12 by the chain rule k[n bⁿ⁻¹ b′ g + bⁿ g′] (1e-9), and the slopes at 1e-12, 1e-8 and 1e-4 increasing. |
| 5 | The temperature gate width: the default bare tanh ramp, narrow gates that sharpen without moving, T_max mirrored, the AD slope, invalid widths refused, the width threaded through `Material`. | kept | "depletion limiter" > "temperature gate width" (values unchanged under `DepletionLimiter(enabled = false)` at ξ = 1). |
| 6 | Negative orders are refused in every constructor form; n = 0 and fractional orders are accepted. | kept | "depletion limiter" > "construction rejects negative orders". |

#### test_radiation.jl

Old-only testset, deleted whole (Section 3.1): "compute_beer_lambert_source!, compute_radiation_source!,
extract_radiative_flux_from_bcs", which holds every old-core name the file uses (the `Mesh1D` fixture,
`initialize_state`, `populate!`, `update_properties!`, the two sources, `extract_radiative_flux_from_bcs`,
`P1_QUASI_STEADY`). The rest of the file runs on `setup(…; radiation_model = BEER_LAMBERT)`, `residual!`, `buffers`
and the lifted kernel.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | The flux entering the medium is the exposed cell's absorptivity times the incident intensity. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": the incident intensity equals `ds.a`·q at the exposed cell, and the exposed cell's intensity equals it. |
| 2 | Under Kirchhoff's law (no component absorptivity) that absorptivity is the exposed cell's emissivity, the one the surface balance emits with, strictly between the components'. | kept | "Beer–Lambert surface absorptivity" > "Kirchhoff path: the exposed cell absorbs at the emissivity it emits with": after `setup(…; radiation_model = BEER_LAMBERT)` and one `residual!`, `ds.a[8] === ds.ε[8]`, `0.85 < ds.ε[8] < 0.95`, `ds.ε[8]` against `mixture_emissivity` (1e-13), and `incident_intensity` and the exposed intensity ≈ ε·50e3. |
| 3 | The absorbed power is a·I·(1 − e^(−τ_tot)) with exact per-cell averages, not the 1/ε-jumped I(1 − e^(−τ)). | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": the per-cell source I_up(1 − e^(−τ)), absorbed plus transmitted equal to the entering beam, which is a·q with a < 1. |
| 4 | In the optically thick limit Beer–Lambert absorbs what surface absorption does. | covered | `test/newcore/test_boundary_conditions.jl` "surface absorption is the surface balance of the incident flux": at α = 1e3 the energy rate under BEER_LAMBERT equals SURFACE_ABSORPTION's (1e-10). |
| 5 | `extract_radiative_flux_from_bcs` folds an explicit α into the intensity and reports it; the Kirchhoff default reports none and the bare flux. | dropped | Deleted with the old core: `extract_radiative_flux_from_bcs` has no caller outside it (the Beer–Lambert kernel reads the `RadiativeFluxBC` itself). Its assertions sit in the old-only "compute_beer_lambert_source!, compute_radiation_source!, extract_radiative_flux_from_bcs" (Section 3.1). |
| 6 | The kernel does not apply an explicit absorptivity twice. | covered | `test/newcore/test_radiation.jl` "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity": with an explicit `absorptivity = t -> 0.9 - 0.01t` the entering beam is (0.9 − 0.01t)·4e4. |
| 7 | No volumetric source under NO_RADIATION and SURFACE_ABSORPTION. | covered | `test/newcore/test_solve.jl` "outputs: the surface-flux decomposition and the mass-loss rate" (no volumetric absorption under SURFACE_ABSORPTION); `test/newcore/test_radiation.jl` "Beer–Lambert is refused off planar coordinates" (no intensity rows under SURFACE_ABSORPTION); `test/newcore/test_boundary_conditions.jl` "surface absorption is the surface balance of the incident flux" (the NO_RADIATION residual equals SURFACE_ABSORPTION's). |
| 8 | No volumetric source under P1_QUASI_STEADY. | dropped | P1 radiation, deleted with Experimental (owner question 3). Its assertion sits in the same old-only testset as row 5. |

### test/geometry

#### test_mesh_generation.jl

Old-core names: `create_stretched_mesh(L, n, Val(2), s; stretch_from)` and the `Node1D` fields `z` and
`boundary_id`. The file goes whole (Section 3.1): every row is covered.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | `stretch_from = :top`: n positive cells summing to L and ending at L, finest at the top, with adjacent ratio s. | covered | `test/newcore/test_meshes.jl` "stretched meshes: the ratio, the refined side and the ends": the ends exactly 0 and L, the adjacent ratio, the finest cell L(r−1)/(rⁿ−1); positivity through the strictly-increasing check of "line meshes: connectivity and tags". |
| 2 | `stretch_from = :bottom` refines the back. | covered | `test/newcore/test_meshes.jl` "stretched meshes: the ratio, the refined side and the ends". |
| 3 | The two orientations are mirror images. | covered | `test/newcore/test_meshes.jl` "stretched meshes: the ratio, the refined side and the ends". |
| 4 | After the reflection the positions are sorted, `:bottom` at z = 0 and `:top` at z = L. | covered | `test/newcore/test_meshes.jl` "line meshes: connectivity and tags" (face 1 `:bottom`, the last face `:top`, non-increasing positions rejected), with "stretched meshes: …"' ends, and "uniform meshes: the nodes, the tags and the volumes" (the tags on a constructor's mesh). |
| 5 | A unit ratio falls back to the uniform mesh. | covered | `test/newcore/test_meshes.jl` "stretched meshes: the ratio, the refined side and the ends": a unit ratio equals the uniform range, also at 1 + 1e-12. |
| 6 | A bad `stretch_from` and a ratio below one are refused. | covered | `test/newcore/test_meshes.jl` "stretched meshes: the ratio, the refined side and the ends": `:side` and 0.9 throw `ArgumentError`. |

#### test_variable_cross_section.jl

Old-core names: none. The problem is built inline on `create_uniform_mesh(0.01, 5)` with a lateral-law material
(`:cross_section_test`), and runs on `setup`, `Pyrolysis.Assembly.residual!` and
`Pyrolysis.Integration.solve`. Besides row 6, the file asserts rows 1 to 5, 7 and 8 on the new core as well.

| # | Obligation | Disposition | Where |
|---|---|---|---|
| 1 | A lateral law adds the cross-section rows to the moving layout, with the stated state length. | covered | `test/newcore/test_states.jl` "state layout: interleaved blocks and row kinds" (the length formula; the history pair and the area and lateral-rate rows partition the state) and "absolute tolerances of every row kind: …" (a lateral-law moving layout has χ rows); `test/newcore/jacobian/test_containment.jl` "the cross section's back-reaction: …" (`has_chi(layout)` after `setup` of a lateral-law deck with `use_ale = true`). Also asserted in "Variable Cross Section" > "a lateral law adds the history pair, the area and the lateral rate" (a state length of 39). |
| 2 | The area follows the conversion through the law (χ = 0.5 gives A = 0.9 A₀). | covered | `test/newcore/test_mesh_motion.jl` "lateral rows: the lateral rate is the area row's logarithmic rate": A = A₀·law(ΣG/ΣM), dA/dt = A·θ_A, the rows at round-off and finite on trial states. Also asserted in the testset of row 1 (at χ̄ = 0.5 the area row is non-zero at A₀ and within 1e-14·A₀ at 0.9 A₀). |
| 3 | The residual is finite, the cross-section rows included. | covered | `test/newcore/test_mesh_motion.jl` "moving mesh: …" (`chi`, `darcy`: residual finite at five states, the area and lateral-rate rows within 1e-14 at the start, containment). Also asserted in the testset of row 1. |
| 4 | Without a lateral law the layout carries no cross-section rows. | covered | `test/newcore/test_mesh_motion.jl` "moving mesh: …": `has_chi(layout) == (name in (:chi, :darcy))` on all nine configurations under `use_ale = true` (`test/jacobian/test_chi_and_reduced.jl` row 2); `test/newcore/test_states.jl` "state layout: interleaved blocks and row kinds" (the law-free lengths). Also asserted in "Variable Cross Section" > "without a lateral law the moving layout carries no cross-section rows". |
| 5 | A lateral law on a fixed mesh is refused. | covered | `test/newcore/test_mesh_motion.jl` "closure: rank, consistent sweep and refusals": `setup` of a lateral-law deck on a fixed mesh throws `ArgumentError`. |
| 6 | The refusal names `lateral_shrinkage_law` and `use_ale = true`, through `solve` with and without `validate`. | kept | "Variable Cross Section" > "a lateral law needs use_ale = true": for `setup(problem; use_ale = false)`, `Integration.solve(problem)` and `Integration.solve(problem; validate = false)` the `ArgumentError` message contains `lateral_shrinkage_law` and `use_ale = true`, and `has_chi` holds under `use_ale = true`. |
| 7 | With `use_ale = true` the cross-section rows are built. | covered | `test/newcore/test_mesh_motion.jl` "balance rows: a uniform strain rate gives θ x, θ r/2, θ r/3" (`setup(…; law, use_ale = true)` reads `lateral_rate_index`, and θ_A matches the closed form); `test/newcore/jacobian/test_containment.jl` "the cross section's back-reaction: …". Also asserted in the testsets of rows 1 and 6. |
| 8 | A variable-cross-section solve succeeds and records a finite area history. | covered | `test/verification/cases/V-CPL-03.jl` QoIs "P7: final lateral area ratio A_f/A₀ = law(χ̄ = 1/2)" (a non-finite ratio fails its tolerance) and "…: retcode" (fails unless `:Success`). Also asserted in "Variable Cross Section" > "a variable-cross-section solve records a finite area history". |

## 3. Testsets and shared definitions deleted with the old core

### 3.1 Testsets and files deleted whole

The kept directories isolate every assertion of a definition that goes with the old core in a testset named after it,
so the deletion removes each of these testsets whole. What they asserted that survives is asserted beside them on the
new core (the kept rows of Section 2). Names used only inside them are qualified there, so no top-level import becomes
unused.

| # | File | Testset |
|---|---|---|
| 1 | `test/physics/test_kinetics.jl` | "TanhLimiter" |
| 2 | `test/physics/test_radiation.jl` | "compute_beer_lambert_source!, compute_radiation_source!, extract_radiative_flux_from_bcs" |
| 3 | `test/materials/test_emissivity_absorptivity.jl` | "effective_absorptivity: Kirchhoff invariant and per-component override" |
| 4 | `test/materials/test_emissivity_absorptivity.jl` | "PropertyCache and update_properties! fill ε and a per cell" |
| 5 | `test/materials/test_gas_viscosity_property.jl` | "effective_gas_viscosity and mixture_gas_viscosity" |
| 6 | `test/boundary_conditions/test_blowing_correction.jl` | "Blowing correction" > "evaluate_bc_flux_and_derivative: face dispatch" > "evaluate_bc_flux_and_derivative(::ConvectiveBC) without a face" |
| 7 | `test/boundary_conditions/test_blowing_correction.jl` | "Blowing correction" > "Structured Jacobian (verify_jacobian_at): strict agreement with blowing" |
| 8 | `test/boundary_conditions/test_boundary_temperature.jl` | "surface Newton: a_eff absorption vs. two-piece split" > "effective_absorptivity on the line" |
| 9 | `test/boundary_conditions/test_boundary_temperature.jl` | "surface Newton: a_eff absorption vs. two-piece split" > "compute_boundary_heat_fluxes(::CombinedThermalBC) without a face" |
| 10 | `test/boundary_conditions/test_boundary_temperature.jl` | "surface Newton ignores the stored state" > "solve_boundary_temperature! without a trial clamp: a NaN cell warns" |
| 11 | `test/boundary_conditions/test_boundary_temperature.jl` | "P1_QUASI_STEADY: RadiativeFluxBC contributes nothing at the face" |
| 12 | `test/boundary_conditions/test_boundary_temperature.jl` | "evaluate_bc_flux_and_derivative without a face: combined and convective" |
| 13 | `test/boundary_conditions/test_boundary_temperature.jl` | "Discretization accessors: get_surface_temperature, check_boundary_energy_balance" |

Two more files go whole: `test/helpers.jl` (no surviving file includes it; its rows are dropped) and
`test/geometry/test_mesh_generation.jl` (every row is covered by `test/newcore/test_meshes.jl`).

Three more testsets go only on the owner's decision (Section 4): `test/physics/test_gas_closures.jl` "Gas closures" >
"porosity is the pore fraction on physical states" and "Gas closures" > "Carman–Kozeny: κ = min(ℓ²φ³/(180(1 − φ)²),
ℓ²)" > "porosity: the three-argument form reads the clamped pore fraction" (C127), and
`test/materials/test_emissivity_absorptivity.jl` "effective_emissivity clamps an unphysical law at 1".

The kept directories run in the old-core group of `test/runtests.jl` until the deletion moves them (D134).

### 3.2 Shared definitions deleted with the old core

These definitions sit in the kept shared layers but have no caller outside the old core, so they are deleted with it
(C127) and their tests are dropped. Each lists the surviving test files that use it, which the deletion then edits;
the files deleted whole (the deleted directories, `test/verification/test_analytical.jl`, `test/helpers.jl` and
`test/geometry/test_mesh_generation.jl`) are not listed. Where a use sits in a testset of Section 3.1, its number
there is given.

**Materials**

- `AbstractNode`, `AbstractCell`, `AbstractFace`, `AbstractMesh`, `AbstractErrorIndicator`, `AbstractQualityMetric`,
  `AbstractPyrolysisSolution`: no surviving test.
- `has_absorption_coefficients`: no surviving test.
- `TanhLimiter` and its two conversion constructors: `test/newcore/test_kinetics.jl` (its import, and "calibration
  guard: rates above the join equal the tanh limiter"), `test/newcore/gates/oldcore_fingerprint.jl`,
  `test/physics/test_kinetics.jl` (testset 1).
- `compute_gas_fluxes_unrolled!`, `compute_gas_fluxes_darcy_fick_unrolled!`, `set_boundary_advective_flux_unrolled!`,
  `apply_species_residual_unrolled!`: no surviving test.
- The enum member `P1_QUASI_STEADY` and the P1 branches that test for it in shared code:
  `test/newcore/test_boundary_conditions.jl` ("in-depth radiation enters through the top faces only", the refusal),
  `test/boundary_conditions/test_boundary_temperature.jl` (testset 11), `test/physics/test_radiation.jl` (testset 2).
- `requires_radiation_cache`, `AbstractRadiationBC`: no surviving test.

**Properties**

- `effective_gas_viscosity`: `test/materials/test_gas_viscosity_property.jl` (testset 5),
  `test/newcore/test_kinetics.jl` ("trial temperature is safe in representative kernels"),
  `test/verification/casedefs.jl` (an import), `test/verification/cases/V-GT-02.jl`.
- `mixture_gas_viscosity`: `test/materials/test_gas_viscosity_property.jl` (testset 5).
- `effective_absorptivity`, `_absorptivity_from_v`: `test/boundary_conditions/test_boundary_temperature.jl` (testset 8),
  `test/materials/test_emissivity_absorptivity.jl` (testsets 3 and 4).
- `harmonic_mean_weighted` (the cut face mean): no surviving test.
- `DarcyCache`, `resize_darcy_cache!`, `PropertyCache` (all four constructors), `_build_property_cache`,
  `resize_cache!`, `update_properties!`: `test/materials/test_emissivity_absorptivity.jl` (testset 4),
  `test/physics/test_radiation.jl` (testset 2).
- The whole of `radiation_cache.jl`: `RadiationCache` (all three constructors), `resize_radiation_cache!`,
  `initialize_radiation_intensity!`, `radiation_energy_density`, `net_radiation_source`, `total_radiation_flux`: no
  surviving test.
- The 4-argument forwarder of `_conductivity_sums`, and the `k_init` keyword of `_conductivity_from_v` and
  `bruggeman_effective_property`: no surviving test.

**BoundaryConditions**

- `extract_radiative_flux_from_bcs`: `test/physics/test_radiation.jl` (testset 2).
- `has_radiative_flux_bc`: no surviving test.
- The allocating `solve_boundary_temperature`: no surviving test.
- `_evaluated(::Nothing, T)` and the `trial === nothing` path of `solve_boundary_temperature!` (the iterate clamp and
  the two warnings): `test/newcore/test_boundary_conditions.jl` ("the surface balance reads the trial clamp: no clamp,
  no warning, finite derivatives", its `trial = nothing` cases), `test/boundary_conditions/test_boundary_temperature.jl`
  (testset 10).
- The face-less `compute_boundary_heat_fluxes(::CombinedThermalBC, …)`:
  `test/boundary_conditions/test_boundary_temperature.jl` (testset 9).
- The face-less `compute_boundary_heat_fluxes(::ConvectiveBC, …)`: no surviving test.
- The face-less `evaluate_bc_flux_and_derivative(::CombinedThermalBC, …, ::RadiationModel)`:
  `test/boundary_conditions/test_boundary_temperature.jl` (testset 12).
- `evaluate_bc_flux_and_derivative(::ConvectiveBC, T_s, t, ε)`: `test/boundary_conditions/test_boundary_temperature.jl`
  (testset 12), `test/boundary_conditions/test_blowing_correction.jl` (testset 6).
- `evaluate_bc_flux_and_derivative(::RadiativeFluxBC, T_s, t, ε)`: no surviving test.

**Physics**

- `_depletion_factor` and `_reactant_factor(conc, ::TanhLimiter, ξ, order)`: reached only through `TanhLimiter` (see
  Materials).
- `conductive_flux`, `compute_conductive_fluxes!` (the whole of `heat_conduction.jl`): no surviving test.
- `compute_beer_lambert_source!`, `compute_radiation_source!` (the whole of `radiation.jl`):
  `test/physics/test_radiation.jl` (testset 2).
- `gas_flux_darcy_fick`, `compute_gas_fluxes!`, `compute_gas_fluxes_darcy_fick!` (`gas_transport.jl` keeps only
  `gas_flux`): no surviving test (a comment in `test/verification/cases/V-GT-01.jl` names `compute_gas_fluxes!`).
- `compute_gas_convective_source!` (`S_conv`) and `_gas_convective_source_unrolled` (the whole of `convection.jl`): no
  surviving test.
- `volume_change_rate(dξ::AbstractVector, T, material)`: no surviving test (the tests use the tuple method, which
  stays).

**Problem**

- `PyrolysisSolution` (the struct, `PyrolysisSolution(problem)`, the call `(sol)(t)` and both `show` methods),
  `_extras`, `_empty_extras`, `_allocate_extras!`, and the six accessor methods on it (`get_times`,
  `get_temperatures`, `get_concentrations`, `get_gas_generation_rate`, `n_timesteps`, `succeeded`; the generic
  functions stay with the `CoreSolution` methods): `test/newcore/gates/burnout_burst.jl` (comments in
  `test/verification/cases/V-GT-02.jl` and `test/verification/exact/mms.jl` name it).
- `get_mesh`, `get_material`, `n_components(::PyrolysisProblem)`, `tspan(::PyrolysisProblem)`: no caller anywhere.
- `apply_initial_conditions!` (both methods) and `populate_m_dry_solid_initial!`:
  `test/verification/cases/V-CPL-06.jl` (an old-core case).
- The `Mesh1D` methods of `_n_cells`, `cell_centers`, `cell_volumes`, `validate_mesh` and `_validate_boundaries!`: no
  surviving test.

## 4. Open points

The points that need the owner's decision or a reviewer's confirmation before the deletion, and the findings that
bear on the rows above.

- `test/solver/test_precompile_workloads.jl` row 1, the one row still **port**: `src/precompile.jl` holds the old
  core's workloads until the deletion, which writes the new workloads and their test together (C121).
- `test/adaptivity/test_surface_depletion.jl` row 11: the new core reports an interior void as `outcome === :interior`
  without a warning. Dropping the user-facing message needs the owner's decision; if it stays, the port is a
  `@test_logs (:warn, r"interior")` in `test/newcore/test_deactivation.jl` "a void under a live cell ends the solve".
- `test/solver/test_integrator_type.jl` row 2: on a fixed mesh Beer–Lambert adds algebraic intensity rows and a mass
  matrix, so a fixed-mesh Beer–Lambert deck compiles a second integrator. The owner should confirm that this is
  accepted.
- `test/solver/test_jacobian_kwarg.jl` rows 1 and 2: the new solve still accepts the `jacobian` keyword and ignores it
  with a warning, although owner question 18 removes it. Once it is removed, the refusal of unknown keywords covers it.
- `test/perf/test_method_count.jl` row 1 and `test/perf/runtests.jl` row 2: the method count runs in the shared test
  process (D132), where it can fail spuriously but never pass wrongly, while plan §5.1 puts it in the fresh-process
  quality test. The owner should confirm the location; if it moves, the one assertion goes into
  `test/newcore/quality/runtests.jl`.
- `porosity` (C127): the clamped pore fraction is still reached, through `validate_darcy_initial_porosity`, the
  Carman–Kozeny branch of the 3-argument `effective_permeability` and `background_gas_concentration`, and the guide
  documents it. The kept and ported state-level assertions read 1 − `skeletal_volume_fraction`
  (`test/physics/test_gas_closures.jl` pressure and Carman–Kozeny, `test/physics/test_background_air.jl` rows 1 and 2,
  `test/materials/test_gas_viscosity_property.jl` row 7), so the decision touches only the two testsets of
  `test/physics/test_gas_closures.jl` named after it (Section 3.1), `test/newcore/test_gas_transport.jl` "Carman–Kozeny
  permeability at the smooth porosity" (the 3-argument permeability of an overfilled state is zero), and
  `test/verification/casedefs.jl` and `test/verification/cases/V-GT-02.jl`, which import it.
  `background_gas_concentration`'s values are the same either way.
- `effective_emissivity` and `effective_heat_capacity` (C127; C2 against the guide): neither has a caller after the
  deletion, but the User Guide documents both for inspecting effective properties. Each kept assertion of
  `effective_emissivity` (`test/materials/test_emissivity_absorptivity.jl` "ε_eff = Σ v_j ε_j(T) at the cell
  temperature", `test/materials/test_gamma_semantics.jl` "depleted-cell floor and liquid ε default") has a
  `mixture_emissivity` assertion of the same property beside it, so deleting the function deletes those lines only;
  `test/verification/exact/mms.jl` reads `effective_heat_capacity`. The clamp goes either way, and with it
  `test/materials/test_emissivity_absorptivity.jl` "effective_emissivity clamps an unphysical law at 1" (row 4).
- `test/ale/test_donor_upwind.jl` rows 2 and 4: the new mesh-relative transport has positivity tests only
  (`test/newcore/test_relative_flux.jl` "a boundary cell beside a steep front is not drained: …" keeps each fraction
  of a near-step binary inside (0, 1)); nothing bounds total variation or overshoot on a translated step. Whether to
  add such a bound is the owner's call.
- Semantic changes counted as covered: on a fixed mesh the dead cell's gas is vented with the event
  (`test/adaptivity/test_cell_park.jl` row 11); a state saved at an event's time is the post-event one
  (`test_cell_park.jl` row 5, `test/adaptivity/test_surface_depletion.jl` row 6); a dead node is frozen with zero
  velocity (`test/residual/test_park_equivalence.jl` row 6). A reviewer should confirm all three.
- `test/verification/test_analytical.jl` rows 2 and 11: the convergence checks run only when `PYROLYSIS_TEST_FAST` is
  not 1, that is in the full matrix, not on pull requests. After the deletion the order behind these rows is asserted
  only there.
- The surface balance's `converged` flag (`test/boundary_conditions/test_boundary_temperature.jl` row 20):
  "Integration with Solver" > "the face balances hold at both faces during a solve, the heated one converged" reads it
  at the heated face. An `AdiabaticBC` face keeps `BoundaryState()`'s default after a solve (`converged = false`,
  `iterations = 0`, T_s = T_cell): the kernel sets its temperature without the Newton, while row 20's obligation has
  `converged = true` at both faces. A prescribed-flux face is closed the same way, so `test/newcore/test_solve.jl`
  "outputs: face temperatures from the boundary balances" never writes its faces' states; the flag means something
  only where the Newton runs, as in "outputs: the surface-flux decomposition and the mass-loss rate". Is the default
  deliberate (no balance to solve), or should the kernel mark such faces converged?
- The enthalpy form's residual depends on the evaluation history through the caloric inversion's warm start (C79,
  deliberate): on a fresh workspace a second evaluation at the same state differs from the first by up to 3e-10
  relative on interior energy rows and 1 ulp on mass rows, on fixed and moving meshes alike; the temperature form is
  bit-identical. `test/newcore/test_boundary_conditions.jl` "the residual carries no surface-balance state between
  evaluations" (`test/solver/test_surface_bc_guards.jl` rows 3 and 4) therefore compares the residual vectors at
  rtol 1e-12 by norm; a per-entry comparison at 1e-12 in the enthalpy form would fail.
- `PoolingFilm`: the pooling configurations' zero-allocation checks hold (`test/newcore/test_gas_transport.jl` "gas
  transport: …", `fick_pooling` and `darcy_pooling`; `test/newcore/test_mesh_motion.jl` "moving mesh: pooling, …"),
  at the initial state, where the pooling root solve runs (D133).
- Outside the configuration of `test/solver/test_boundary_flux_reporting.jl` row 8: with faster kinetics
  (E = 6.5e4 J/mol) its deck fully converts the top cell on a fixed mesh, whose mixture absorptivity and emissivity
  then fall to about 2e-5 (a cell holding only gas that takes no volume); on a moving mesh without deactivation the
  same deck stops with `DtLessThanMin` near 2e5 K. Neither is asserted; whether either needs a guard is the owner's
  call.
- Tests the ledger credits that still use names the deletion removes. Each is rewritten before or with the deletion,
  or the rows citing it lose their coverage; the verification suite's routing (`core_solve`, `core_mesh`,
  `solve_mms`'s `core`) is the deletion's (C126):
  - `test/newcore/test_kinetics.jl`: the `TanhLimiter` import and "calibration guard: rates above the join equal the
    tanh limiter", and `effective_gas_viscosity` in "trial temperature is safe in representative kernels".
  - `test/newcore/test_boundary_conditions.jl`: the `P1_QUASI_STEADY` refusal in "in-depth radiation enters through the
    top faces only", and the `trial = nothing` cases of "the surface balance reads the trial clamp: no clamp, no
    warning, finite derivatives".
  - `test/verification/casedefs.jl` and `test/verification/cases/V-GT-02.jl`: `effective_gas_viscosity` and the clamped
    `porosity`, for the exact profile.
  - `test/verification/exact/mms.jl`: `Pyrolysis.Solver: KenCarp4, AutoFiniteDiff`, `initialize_state` (the old-core
    branch of `solve_mms`) and `effective_heat_capacity`.
  - `test/verification/cases/V-BO-02.jl`: `Pyrolysis.Solver.DT_AFTER_PARK`, its old-core budget. `V-CPL-06.jl` is an
    old-core case (`Pyrolysis.Solver`, `apply_initial_conditions!`) and goes with the old core.
  - The gate scripts `test/newcore/gates/burnout_burst.jl` (`initialize_state`, `PyrolysisSolution`,
    `Pyrolysis.Solver`), `audit_breakdown.jl` (`Pyrolysis.Solver: KenCarp4`) and `oldcore_fingerprint.jl`
    (`TanhLimiter`); no row cites them.
- Comments in the new core's tests that name old-core files or describe the old core, to reword at the deletion: in
  `test/newcore/jacobian/test_containment.jl` the comments above `_pmma` and above `_incident` and
  `_perf_configurations` (citing `test/perf/perf_helpers.jl` and the old core's `radiation_intensity` keyword) and
  above `_states` (citing `test/jacobian/test_regressions.jl`), and the `hot`-state comment in "perf configuration …"
  ("the block the old core's prefix sum carried"); the comment above the configurations of
  `test/newcore/jacobian/test_allocations.jl` (citing `test/perf/perf_helpers.jl`); and "fail on either core" in the
  comment above `test/newcore/test_energy_ledger.jl`'s `_column`.
