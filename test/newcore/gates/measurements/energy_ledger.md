# Gate item 7: the energy balance against the ledger

Produced by `test/newcore/gates/energy_ledger.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [20.7, 19.9, 20.3] at the start and [23.8, 22.1, 21.2] at the end. V-CPL-06's char-layer column.

## Closure at the case's tolerance

80 cells, abstol 1.0e-8, reltol 1.0e-6, the window t ≥ 120.0 s; |ΔE − ∫ ledger| over the integrated boundary heat, trapezoidal over the accepted steps on both cores.

| deck | run | retcode | closure |
|---|---|---|---|
| Eulerian | old core | Success | 2.27e-5 |
| Eulerian | new core, temperature form | Success | 5.9e-7 |
| Eulerian | new core, enthalpy form | Success | 1.03e-6 |
| moving | old core | Success | 4.69e-5 |
| moving | new core, temperature form | Success | 3.08e-6 |
| moving | new core, enthalpy form | Success | 4.32e-6 |

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
