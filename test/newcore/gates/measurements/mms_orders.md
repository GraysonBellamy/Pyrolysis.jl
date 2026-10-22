# Gate item 3: the manufactured-solution orders

Produced by `test/newcore/gates/mms_orders.jl` at commit f54a6fb (tree b544f02). Julia 1.13.0, 1 thread, 22 logical CPUs; load average [19.7, 19.0, 18.9] at the start and [19.3, 19.3, 19.2] at the end. Each core on its default integrator.

## Spatial orders (each case's convergence check)

| case | run | cells | errors | observed orders |
|---|---|---|---|---|
| V-MMS-01 | old core | 20, 40, 80 | 0.0276, 0.00699, 0.00176 | 1.98, 1.99 |
| V-MMS-03 | old core | 16, 32, 64 | 0.234, 0.0508, 0.0109 | 2.21, 2.22 |
| V-MMS-04 | old core | 16, 32, 64 | 0.000753, 0.000188, 4.67e-5 | 2.0, 2.01 |
| V-MMS-01 | new core, temperature form | 20, 40, 80 | 0.0253, 0.00641, 0.00161 | 1.98, 2.0 |
| V-MMS-03 | new core, temperature form | 16, 32, 64 | 0.0504, 0.00443, 0.00145 | 3.51, 1.61 |
| V-MMS-04 | new core, temperature form | 16, 32, 64 | 0.000753, 0.000188, 4.7e-5 | 2.0, 2.0 |
| V-MMS-01 | new core, enthalpy form | 20, 40, 80 | 0.0253, 0.00641, 0.00161 | 1.98, 2.0 |
| V-MMS-03 | new core, enthalpy form | 16, 32, 64 | 0.0504, 0.00443, 0.00145 | 3.51, 1.61 |
| V-MMS-04 | new core, enthalpy form | 16, 32, 64 | 0.000753, 0.000188, 4.7e-5 | 2.0, 2.0 |

## V-MMS-01 temporal sweep

| run | cells | dt [s] | L2 [K] | L∞ [K] | L2 orders | L∞ orders |
|---|---|---|---|---|---|---|
| old core | 80 | 0.5, 0.25, 0.125, 0.0625 | 3.17e-6, 1.98e-7, 1.23e-8, 8.21e-10 | 4.6e-6, 2.88e-7, 1.82e-8, 1.21e-9 | 4.0, 4.01, 3.91 | 4.0, 3.98, 3.91 |
| new core, temperature form | 80 | 0.5, 0.25, 0.125, 0.0625 | 2.53e-7, 1.67e-8, 8.62e-10, 3.94e-11 | 4.84e-7, 3.33e-8, 2.14e-9, 1.71e-10 | 3.92, 4.28, 4.45 | 3.86, 3.96, 3.65 |
| new core, enthalpy form | 80 | 0.5, 0.25, 0.125, 0.0625 | 8.67e-8, 3.97e-9, 2.43e-10, 1.77e-11 | 2.06e-7, 1.68e-8, 1.35e-9, 1.34e-10 | 4.45, 4.03, 3.78 | 3.61, 3.64, 3.34 |

## Quantities of interest at the CI resolution

| case | run | retcode | quantity | error | tolerance | passes |
|---|---|---|---|---|---|---|
| V-MMS-01 | old core | Success | L∞ T error, t=10.0 s | 0.0165 | 0.05 | true |
| V-MMS-01 | old core | Success | L∞ T error, t=20.0 s | 0.0122 | 0.05 | true |
| V-MMS-01 | old core | Success | relative L∞ ξ_A error, t=20.0 s | 6.81e-5 | 0.0002 | true |
| V-MMS-01 | old core | Success | relative L∞ ξ_B error, t=20.0 s | 2.19e-5 | 0.0001 | true |
| V-MMS-01 | old core | Success | relative L∞ ξ_gas error, t=20.0 s | 0.00991 | 0.03 | true |
| V-MMS-02 | old core | Success | free-stream T drift (all cells, all snapshots) | 5.68e-14 | 1.0e-10 | true |
| V-MMS-02 | old core | Success | free-stream ξ_solid relative drift | 1.89e-16 | 1.0e-12 | true |
| V-MMS-02 | old core | Success | free-stream ξ_gas relative drift | 1.48e-16 | 1.0e-12 | true |
| V-MMS-02 | old core | Success | top cell-center trajectory vs z₀·s(t_end) | 4.04e-7 | 1.5e-6 | true |
| V-MMS-03 | old core | Success | L∞ T error on moving mesh, t=10.0 s | 0.0711 | 0.25 | true |
| V-MMS-03 | old core | Success | relative L∞ ξ_solid (mesh-bound) error, t=10.0 s | 2.11e-11 | 1.0e-10 | true |
| V-MMS-03 | old core | Success | relative L∞ ξ_gas error, t=10.0 s | 0.00178 | 0.006 | true |
| V-MMS-04 | old core | Success | L∞ T error (DirectSolve), t=20.0 s | 2.03e-5 | 6.0e-5 | true |
| V-MMS-04 | old core | Success | relative L∞ ξ_solid (z-varying, material frame) error, t=20.0 s | 6.57e-7 | 2.0e-6 | true |
| V-MMS-04 | old core | Success | relative L∞ ξ_gas error, t=20.0 s | 0.0023 | 0.007 | true |
| V-MMS-04 | old core | Success | relative surface-position error, t=20.0 s | 1.63e-5 | 5.0e-5 | true |
| V-MMS-04 | old core | Success | L∞ T error (ApproxSolve), t=20.0 s | 2.04e-5 | 6.0e-5 | true |
| V-MMS-04 | old core | Success | relative L∞ ξ_solid error (ApproxSolve), t=20.0 s | 6.57e-7 | 2.0e-6 | true |
| V-MMS-01 | new core, temperature form | Success | L∞ T error, t=10.0 s | 0.0187 | 0.05 | true |
| V-MMS-01 | new core, temperature form | Success | L∞ T error, t=20.0 s | 0.0135 | 0.05 | true |
| V-MMS-01 | new core, temperature form | Success | relative L∞ ξ_A error, t=20.0 s | 7.51e-5 | 0.0002 | true |
| V-MMS-01 | new core, temperature form | Success | relative L∞ ξ_B error, t=20.0 s | 2.41e-5 | 0.0001 | true |
| V-MMS-01 | new core, temperature form | Success | relative L∞ ξ_gas error, t=20.0 s | 0.0111 | 0.03 | true |
| V-MMS-02 | new core, temperature form | Success | free-stream T drift (all cells, all snapshots) | 3.41e-13 | 1.0e-10 | true |
| V-MMS-02 | new core, temperature form | Success | free-stream ξ_solid relative drift | 1.84e-14 | 1.0e-12 | true |
| V-MMS-02 | new core, temperature form | Success | free-stream ξ_gas relative drift | 2.92e-14 | 1.0e-12 | true |
| V-MMS-02 | new core, temperature form | Success | top cell-center trajectory vs z₀·s(t_end) | 7.32e-9 | 1.5e-6 | true |
| V-MMS-03 | new core, temperature form | Success | L∞ T error on moving mesh, t=10.0 s | 0.00909 | 0.25 | true |
| V-MMS-03 | new core, temperature form | Success | relative L∞ ξ_solid (mesh-bound) error, t=10.0 s | 5.32e-12 | 1.0e-10 | true |
| V-MMS-03 | new core, temperature form | Success | relative L∞ ξ_gas error, t=10.0 s | 0.00168 | 0.006 | true |
| V-MMS-04 | new core, temperature form | Success | L∞ T error (DirectSolve), t=20.0 s | 2.07e-5 | 6.0e-5 | true |
| V-MMS-04 | new core, temperature form | Success | relative L∞ ξ_solid (z-varying, material frame) error, t=20.0 s | 6.59e-7 | 2.0e-6 | true |
| V-MMS-04 | new core, temperature form | Success | relative L∞ ξ_gas error, t=20.0 s | 0.0023 | 0.007 | true |
| V-MMS-04 | new core, temperature form | Success | relative surface-position error, t=20.0 s | 1.63e-5 | 5.0e-5 | true |
| V-MMS-01 | new core, enthalpy form | Success | L∞ T error, t=10.0 s | 0.0187 | 0.05 | true |
| V-MMS-01 | new core, enthalpy form | Success | L∞ T error, t=20.0 s | 0.0135 | 0.05 | true |
| V-MMS-01 | new core, enthalpy form | Success | relative L∞ ξ_A error, t=20.0 s | 7.51e-5 | 0.0002 | true |
| V-MMS-01 | new core, enthalpy form | Success | relative L∞ ξ_B error, t=20.0 s | 2.41e-5 | 0.0001 | true |
| V-MMS-01 | new core, enthalpy form | Success | relative L∞ ξ_gas error, t=20.0 s | 0.0111 | 0.03 | true |
| V-MMS-02 | new core, enthalpy form | Success | free-stream T drift (all cells, all snapshots) | 7.39e-13 | 1.0e-10 | true |
| V-MMS-02 | new core, enthalpy form | Success | free-stream ξ_solid relative drift | 2.5e-14 | 1.0e-12 | true |
| V-MMS-02 | new core, enthalpy form | Success | free-stream ξ_gas relative drift | 5.0e-14 | 1.0e-12 | true |
| V-MMS-02 | new core, enthalpy form | Success | top cell-center trajectory vs z₀·s(t_end) | 6.33e-9 | 1.5e-6 | true |
| V-MMS-03 | new core, enthalpy form | Success | L∞ T error on moving mesh, t=10.0 s | 0.00908 | 0.25 | true |
| V-MMS-03 | new core, enthalpy form | Success | relative L∞ ξ_solid (mesh-bound) error, t=10.0 s | 3.8e-12 | 1.0e-10 | true |
| V-MMS-03 | new core, enthalpy form | Success | relative L∞ ξ_gas error, t=10.0 s | 0.00168 | 0.006 | true |
| V-MMS-04 | new core, enthalpy form | Success | L∞ T error (DirectSolve), t=20.0 s | 2.07e-5 | 6.0e-5 | true |
| V-MMS-04 | new core, enthalpy form | Success | relative L∞ ξ_solid (z-varying, material frame) error, t=20.0 s | 6.59e-7 | 2.0e-6 | true |
| V-MMS-04 | new core, enthalpy form | Success | relative L∞ ξ_gas error, t=20.0 s | 0.0023 | 0.007 | true |
| V-MMS-04 | new core, enthalpy form | Success | relative surface-position error, t=20.0 s | 1.63e-5 | 5.0e-5 | true |
