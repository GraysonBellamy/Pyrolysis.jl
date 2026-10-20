# The energy balance against the ledger

Produced by `test/newcore/gates/energy_ledger.jl` at commit 7ceaf74. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [18.2, 17.1, 17.8] at the start and [17.3, 17.4, 17.6] at the end. Phase 1's gate item 7 on V-CPL-06's char-layer column, and Phase 2's item 3 on V-BO-04's burning column before and across its event.

## Closure at the case's tolerance

80 cells, abstol 1.0e-8, reltol 1.0e-6, the window t ≥ 120.0 s; |ΔE − ∫ ledger| over the integrated boundary heat, trapezoidal over the accepted steps on both cores.

| deck | run | retcode | closure |
|---|---|---|---|
| Eulerian | old core | Success | 2.27e-5 |
| Eulerian | new core, temperature form | Success | 5.91e-7 |
| Eulerian | new core, enthalpy form | Success | 7.62e-7 |
| moving | old core | Success | 4.69e-5 |
| moving | new core, temperature form | Success | 3.08e-6 |
| moving | new core, enthalpy form | Success | 4.26e-6 |

## Refinement with the ledger integrated by the integrator

20 cells, 0 to 240 s at fixed steps; the largest |Σ H − H₀ − E| over the solve [J] against the largest |E|.

| deck | form | dt [s] | retcode | defect [J] | scale [J] | order |
|---|---|---|---|---|---|---|
| Eulerian | temperature | 2.0 | Success | 0.000297 | 1.24e6 |  |
| Eulerian | temperature | 1.0 | Success | 1.7e-5 | 1.24e6 | 4.12 |
| Eulerian | temperature | 0.5 | Success | 8.71e-7 | 1.24e6 | 4.29 |
| Eulerian | enthalpy | 2.0 | Success | 4.07e-9 | 1.24e6 |  |
| moving | temperature | 2.0 | Success | 0.000297 | 1.24e6 |  |
| moving | temperature | 1.0 | Success | 1.71e-5 | 1.24e6 | 4.12 |
| moving | temperature | 0.5 | Success | 9.15e-7 | 1.24e6 | 4.22 |
| moving | enthalpy | 2.0 | Success | 3.26e-9 | 1.24e6 |  |

## V-BO-04 before its first event

40 cells on a moving mesh, abstol 1.0e-9, reltol 1.0e-7 (the case's), 0 to 290.0 s without deactivation; |ΔE − ∫ ledger| over |ΔE|, each core's own, trapezoidal over the accepted steps (the old core's budget is V-CPL-06b's).

| run | retcode | closure |
|---|---|---|
| old core | Success | 7.1e-5 |
| new core, temperature form | Success | 3.65e-6 |
| new core, enthalpy form | Success | 3.98e-6 |

20 cells, 0 to 290.0 s at fixed steps; the largest |Σ H − H₀ − E| over the solve [J] against the largest |E|, and their ratio. The column stays at its initial temperature, so its enthalpy is linear in the masses in both forms and the ledger an invariant the integrator keeps to round-off; an order is given only between defects above 1.0e-12 of the scale.

| form | dt [s] | retcode | defect [J] | scale [J] | ratio | order |
|---|---|---|---|---|---|---|
| temperature | 0.1 | Success | 1.19e-8 | 2.85e6 | 4.16e-15 |  |
| temperature | 0.05 | Success | 2.85e-7 | 2.85e6 | 9.97e-14 |  |
| temperature | 0.025 | Success | 7.77e-7 | 2.85e6 | 2.72e-13 |  |
| enthalpy | 0.1 | Success | 2.33e-8 | 2.85e6 | 8.16e-15 |  |

## V-BO-04 across its event

The case's run, 40 cells on a moving mesh with deactivation, abstol 1.0e-9, reltol 1.0e-7, event location off and on. Inventories over the live cells (the extended energy, the inventory itself in the enthalpy form, and the mass), the boundary terms at the live boundary faces, the sources over the live cells, trapezoidal over the accepted steps with a node where the smallest threshold margin is 1.0e-8; each invocation's lost enthalpy and regularization and its lost masses leave the inventories at the event. Closure: |Δ − ∫ ledger + lost| over |Δ|, energy and mass. Event: the first invocation's time [s]; lost energy, regularization: the events' sums over |ΔE|; created volume: the events' sum [m³]. The old core's budget does not park cells; its closure before the event is the table above.

| form | location | retcode | outcome | invocations | events | event [s] | energy closure | mass closure | lost energy | regularization | created volume | steps |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| temperature | off | Success | running | 1 | 39 | 473.0847 | 0.00661 | 0.00661 | 1.22e-18 | 1.22e-18 | 0.0 | 325 |
| temperature | on | Success | running | 1 | 39 | 465.4896 | 1.58e-5 | 1.58e-5 | 1.22e-18 | 1.22e-18 | 0.0 | 325 |
| enthalpy | off | Success | running | 1 | 39 | 471.6878 | 0.00555 | 0.00555 | 2.43e-18 | 2.43e-18 | 2.71e-20 | 328 |
| enthalpy | on | Success | running | 1 | 39 | 465.4896 | 1.46e-5 | 1.46e-5 | 3.65e-18 | 3.65e-18 | 5.42e-20 | 328 |
