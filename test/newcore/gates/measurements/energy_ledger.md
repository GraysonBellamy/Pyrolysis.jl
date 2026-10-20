# The energy balance against the ledger

Produced by `test/newcore/gates/energy_ledger.jl` at commit 675b4f0. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [15.2, 15.8, 16.2] at the start and [16.8, 17.0, 16.7] at the end. Phase 1's gate item 7 on V-CPL-06's char-layer column, and Phase 2's item 3 on V-BO-04's burning column before its first event.

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
