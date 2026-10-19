# Gate item 2: the audit decks at matched accuracy

Produced by `test/newcore/gates/audit_decks.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [24.7, 22.9, 21.0] at the start and [20.7, 21.9, 21.4] at the end. Each core's default integrator; error: the largest temperature difference [K] over the saved times and cells against the same core at reltol 1.0e-11 (abstol a hundredth of reltol); spread: that reference against one at a tenth of its tolerances. Target: the old core's error at its default tolerances (reltol 1.0e-6). Matched: the new core's cheapest ladder solve at or below the target. Warm medians of 3 solves [s].

| deck | new core | target [K] | spread [K] | old [s] | equal-tolerance ratio | new error at default [K] | matched reltol | matched error [K] | matched [s] | matched ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | temperature form | 9.9e-6 | 1.38e-9 | 0.00926 | 0.719 | 5.09e-5 | 1.0e-7 | 2.95e-6 | 0.011 | 1.18 |
| eulerian 25, solid λ 1.0e-6 | enthalpy form | 9.9e-6 | 1.38e-9 | 0.00926 | 3.1 | 9.17e-5 | 1.0e-7 | 3.69e-6 | 0.0429 | 4.63 |
| eulerian 100, solid λ 1.0e-6 | temperature form | 1.37e-5 | 5.42e-9 | 0.0406 | 0.884 | 0.000297 | 1.0e-7 | 9.77e-6 | 0.0584 | 1.44 |
| eulerian 100, solid λ 1.0e-6 | enthalpy form | 1.37e-5 | 5.42e-9 | 0.0406 | 3.28 | 6.12e-5 | 1.0e-7 | 7.36e-6 | 0.257 | 6.33 |
| ale 25, solid λ 1.0e-6 | temperature form | 2.31e-5 | 1.22e-9 | 0.0139 | 1.42 | 6.82e-5 | 1.0e-7 | 3.09e-6 | 0.0292 | 2.09 |
| ale 25, solid λ 1.0e-6 | enthalpy form | 2.31e-5 | 1.22e-9 | 0.0139 | 2.82 | 0.000124 | 3.0e-7 | 6.91e-6 | 0.047 | 3.37 |
| ale 100, solid λ 1.0e-6 | temperature form | 1.03e-5 | 3.14e-9 | 0.0448 | 1.69 | 9.66e-5 | 1.0e-7 | 6.54e-6 | 0.124 | 2.78 |
| ale 100, solid λ 1.0e-6 | enthalpy form | 1.03e-5 | 3.14e-9 | 0.0448 | 3.15 | 0.000162 | 1.0e-7 | 4.76e-6 | 0.264 | 5.9 |
| ale_darcy 25 | temperature form | 0.000126 | 3.03e-8 | 0.0192 | 1.91 | 6.77e-5 | 1.0e-6 | 6.77e-5 | 0.0367 | 1.91 |
| ale_darcy 25 | enthalpy form | 0.000126 | 3.03e-8 | 0.0192 | 2.87 | 7.55e-5 | 1.0e-6 | 7.55e-5 | 0.0552 | 2.87 |
| ale_darcy 100 | temperature form | 5.96e-5 | 3.0e-8 | 0.0775 | 1.67 | 0.000232 | 3.0e-7 | 1.44e-5 | 0.18 | 2.32 |
| ale_darcy 100 | enthalpy form | 5.96e-5 | 3.0e-8 | 0.0775 | 2.82 | 7.2e-5 | 3.0e-7 | 2.78e-5 | 0.291 | 3.76 |

## The original PMMA decks (solid λ = 0) at the default tolerances

No tolerance converges these decks (see the script's header); warm times only, and the ratio to the old core's.

| deck | run | retcode | steps | warm [s] | ratio |
|---|---|---|---|---|---|
| eulerian 25 | old core | Success | 292 | 0.0397 | 1.0 |
| eulerian 25 | new core, temperature form | Success | 239 | 0.0227 | 0.572 |
| eulerian 25 | new core, enthalpy form | Success | 272 | 0.0589 | 1.48 |
| eulerian 100 | old core | Success | 941 | 0.405 | 1.0 |
| eulerian 100 | new core, temperature form | Success | 989 | 0.336 | 0.83 |
| eulerian 100 | new core, enthalpy form | Success | 1037 | 0.85 | 2.1 |
| ale 25 | old core | Success | 196 | 0.0311 | 1.0 |
| ale 25 | new core, temperature form | Success | 254 | 0.0543 | 1.74 |
| ale 25 | new core, enthalpy form | Success | 261 | 0.0896 | 2.88 |
| ale 100 | old core | Success | 1660 | 0.775 | 1.0 |
| ale 100 | new core, temperature form | Success | 1118 | 0.872 | 1.13 |
| ale 100 | new core, enthalpy form | Success | 837 | 1.13 | 1.46 |

## Ladders

| deck | run | reltol | retcode | error [K] | steps | warm [s] |
|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000329 | 23 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000159 | 28 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 3.42e-5 | 34 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 6.43e-5 | 42 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 9.9e-6 | 52 | 0.00926 |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 1.79e-5 | 64 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 3.77e-6 | 80 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 9.99e-7 | 103 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 9.73e-7 | 128 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.00108 | 25 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00121 | 29 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000205 | 33 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000194 | 40 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 5.09e-5 | 48 | 0.00666 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 1.72e-5 | 64 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 2.95e-6 | 84 | 0.011 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 7.2e-7 | 113 | 0.014 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 7.43e-7 | 153 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.00183 | 33 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.000373 | 38 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.000763 | 44 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.00023 | 54 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 9.17e-5 | 66 | 0.0287 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 1.77e-5 | 84 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 3.69e-6 | 114 | 0.0429 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 1.37e-6 | 164 | 0.0629 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 2.19e-7 | 229 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000204 | 27 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000107 | 32 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 3.52e-5 | 38 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 1.4e-5 | 47 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 1.37e-5 | 56 | 0.0406 |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 8.36e-6 | 69 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 2.0e-6 | 85 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 1.66e-6 | 106 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 5.8e-7 | 133 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0115 | 31 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00142 | 36 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000719 | 42 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000305 | 49 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 0.000297 | 58 | 0.0359 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 2.83e-5 | 75 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 9.77e-6 | 98 | 0.0584 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 3.21e-6 | 129 | 0.0746 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 1.08e-6 | 165 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.00665 | 42 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00184 | 49 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.000808 | 57 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.000344 | 70 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 6.12e-5 | 85 | 0.133 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 1.39e-5 | 113 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 7.36e-6 | 168 | 0.257 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 1.61e-6 | 291 | 0.452 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 2.48e-6 | 492 |  |
| ale 25, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000553 | 23 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000236 | 27 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 9.82e-5 | 33 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 5.35e-5 | 41 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 2.31e-5 | 51 | 0.0139 |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 4.91e-6 | 63 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 2.46e-6 | 78 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 8.29e-7 | 99 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 3.26e-7 | 124 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0027 | 25 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.000464 | 29 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000703 | 33 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000124 | 39 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 6.82e-5 | 46 | 0.0198 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 2.76e-5 | 59 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 3.09e-6 | 79 | 0.0292 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 5.26e-7 | 108 | 0.0402 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 2.01e-7 | 141 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.000542 | 32 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00154 | 36 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.000419 | 42 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 3.51e-5 | 52 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 0.000124 | 63 | 0.0393 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 6.91e-6 | 80 | 0.047 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 4.06e-6 | 107 | 0.0617 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 3.72e-7 | 151 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 3.73e-7 | 213 |  |
| ale 100, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000608 | 28 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000199 | 32 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 0.000106 | 37 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 3.09e-5 | 46 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 1.03e-5 | 54 | 0.0448 |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 7.31e-6 | 67 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 4.43e-6 | 81 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 9.1e-7 | 104 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 5.61e-7 | 129 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0159 | 30 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00144 | 34 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.00109 | 39 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000164 | 47 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 9.66e-5 | 54 | 0.0759 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 1.18e-5 | 70 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 6.54e-6 | 90 | 0.124 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 5.45e-7 | 117 | 0.16 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 3.45e-7 | 153 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.0446 | 40 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00256 | 47 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.00146 | 54 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.000445 | 66 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 0.000162 | 80 | 0.141 |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 4.59e-5 | 105 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 4.76e-6 | 151 | 0.264 |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 2.9e-6 | 252 | 0.445 |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 3.79e-7 | 427 |  |
| ale_darcy 25 | old core | 0.0001 | Success | 0.000278 | 79 |  |
| ale_darcy 25 | old core | 3.0e-5 | Success | 0.000131 | 59 |  |
| ale_darcy 25 | old core | 1.0e-5 | Success | 0.000132 | 63 |  |
| ale_darcy 25 | old core | 3.0e-6 | Success | 0.000154 | 69 |  |
| ale_darcy 25 | old core | 1.0e-6 | Success | 0.000126 | 75 | 0.0192 |
| ale_darcy 25 | old core | 3.0e-7 | Success | 3.74e-5 | 90 |  |
| ale_darcy 25 | old core | 1.0e-7 | Success | 3.1e-5 | 107 |  |
| ale_darcy 25 | old core | 3.0e-8 | Success | 5.31e-5 | 122 |  |
| ale_darcy 25 | old core | 1.0e-8 | Success | 6.22e-6 | 166 |  |
| ale_darcy 25 | new core, temperature form | 0.0001 | Success | 0.00773 | 45 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-5 | Success | 0.00181 | 54 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-5 | Success | 0.000685 | 62 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-6 | Success | 0.00053 | 76 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-6 | Success | 6.77e-5 | 100 | 0.0367 |
| ale_darcy 25 | new core, temperature form | 3.0e-7 | Success | 5.53e-5 | 135 | 0.0439 |
| ale_darcy 25 | new core, temperature form | 1.0e-7 | Success | 9.75e-6 | 186 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-8 | Success | 1.73e-6 | 259 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-8 | Success | 5.77e-7 | 353 |  |
| ale_darcy 25 | new core, enthalpy form | 0.0001 | Success | 0.0137 | 50 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-5 | Success | 0.00177 | 61 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-5 | Success | 0.0011 | 75 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-6 | Success | 0.000157 | 102 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-6 | Success | 7.55e-5 | 128 | 0.0552 |
| ale_darcy 25 | new core, enthalpy form | 3.0e-7 | Success | 9.35e-6 | 180 | 0.0779 |
| ale_darcy 25 | new core, enthalpy form | 1.0e-7 | Success | 3.85e-6 | 240 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-8 | Success | 1.19e-6 | 337 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-8 | Success | 3.25e-7 | 479 |  |
| ale_darcy 100 | old core | 0.0001 | Success | 0.00022 | 589 |  |
| ale_darcy 100 | old core | 3.0e-5 | Success | 0.000287 | 180 |  |
| ale_darcy 100 | old core | 1.0e-5 | Success | 0.000124 | 90 |  |
| ale_darcy 100 | old core | 3.0e-6 | Success | 0.000162 | 74 |  |
| ale_darcy 100 | old core | 1.0e-6 | Success | 5.96e-5 | 80 | 0.0775 |
| ale_darcy 100 | old core | 3.0e-7 | Success | 0.000255 | 92 |  |
| ale_darcy 100 | old core | 1.0e-7 | Success | 8.63e-5 | 98 |  |
| ale_darcy 100 | old core | 3.0e-8 | Success | 4.38e-5 | 126 |  |
| ale_darcy 100 | old core | 1.0e-8 | Success | 9.52e-6 | 169 |  |
| ale_darcy 100 | new core, temperature form | 0.0001 | Success | 0.0116 | 48 |  |
| ale_darcy 100 | new core, temperature form | 3.0e-5 | Success | 0.00192 | 56 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-5 | Success | 0.00175 | 66 |  |
| ale_darcy 100 | new core, temperature form | 3.0e-6 | Success | 0.000225 | 80 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-6 | Success | 0.000232 | 102 | 0.13 |
| ale_darcy 100 | new core, temperature form | 3.0e-7 | Success | 1.44e-5 | 140 | 0.18 |
| ale_darcy 100 | new core, temperature form | 1.0e-7 | Success | 5.48e-6 | 188 | 0.241 |
| ale_darcy 100 | new core, temperature form | 3.0e-8 | Success | 4.67e-6 | 260 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-8 | Success | 4.61e-7 | 350 |  |
| ale_darcy 100 | new core, enthalpy form | 0.0001 | Success | 0.0101 | 55 |  |
| ale_darcy 100 | new core, enthalpy form | 3.0e-5 | Success | 0.00305 | 65 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-5 | Success | 0.00222 | 78 |  |
| ale_darcy 100 | new core, enthalpy form | 3.0e-6 | Success | 0.000911 | 96 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-6 | Success | 7.2e-5 | 126 | 0.218 |
| ale_darcy 100 | new core, enthalpy form | 3.0e-7 | Success | 2.78e-5 | 170 | 0.291 |
| ale_darcy 100 | new core, enthalpy form | 1.0e-7 | Success | 5.8e-6 | 251 | 0.429 |
| ale_darcy 100 | new core, enthalpy form | 3.0e-8 | Success | 1.8e-6 | 405 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-8 | Success | 8.48e-7 | 643 |  |
