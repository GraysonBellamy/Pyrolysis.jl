# Gate item 2: the audit decks at matched accuracy

Produced by `test/newcore/gates/audit_decks.jl` at commit 603e863. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [13.5, 14.3, 13.5] at the start and [12.2, 13.5, 13.4] at the end. Each core's default integrator; error: the largest temperature difference [K] over the saved times and cells against the same core at reltol 1.0e-11 (abstol a hundredth of reltol); spread: that reference against one at a tenth of its tolerances. Target: the old core's error at its default tolerances (reltol 1.0e-6). Matched: the new core's cheapest ladder solve at or below the target. Warm medians of 3 solves [s].

| deck | new core | target [K] | spread [K] | old [s] | equal-tolerance ratio | new error at default [K] | matched reltol | matched error [K] | matched [s] | matched ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | temperature form | 9.9e-6 | 1.38e-9 | 0.00418 | 0.747 | 5.09e-5 | 1.0e-7 | 2.95e-6 | 0.00484 | 1.16 |
| eulerian 25, solid λ 1.0e-6 | enthalpy form | 9.9e-6 | 1.38e-9 | 0.00418 | 2.52 | 9.17e-5 | 1.0e-7 | 3.69e-6 | 0.0177 | 4.23 |
| eulerian 100, solid λ 1.0e-6 | temperature form | 1.37e-5 | 5.42e-9 | 0.0132 | 0.833 | 0.000297 | 1.0e-7 | 9.77e-6 | 0.0174 | 1.32 |
| eulerian 100, solid λ 1.0e-6 | enthalpy form | 1.37e-5 | 5.42e-9 | 0.0132 | 3.92 | 6.12e-5 | 1.0e-7 | 7.36e-6 | 0.103 | 7.78 |
| ale 25, solid λ 1.0e-6 | temperature form | 2.31e-5 | 1.22e-9 | 0.00455 | 0.954 | 6.82e-5 | 1.0e-7 | 3.09e-6 | 0.00675 | 1.48 |
| ale 25, solid λ 1.0e-6 | enthalpy form | 2.31e-5 | 1.22e-9 | 0.00455 | 2.9 | 0.000108 | 3.0e-7 | 1.44e-5 | 0.0162 | 3.56 |
| ale 100, solid λ 1.0e-6 | temperature form | 1.03e-5 | 3.14e-9 | 0.0149 | 1.19 | 9.66e-5 | 1.0e-7 | 6.54e-6 | 0.0281 | 1.88 |
| ale 100, solid λ 1.0e-6 | enthalpy form | 1.03e-5 | 3.14e-9 | 0.0149 | 4.05 | 0.000188 | 3.0e-8 | 1.84e-6 | 0.194 | 13.0 |
| ale_darcy 25 | temperature form | 0.000126 | 3.03e-8 | 0.011 | 1.06 | 7.28e-5 | 1.0e-6 | 7.28e-5 | 0.0116 | 1.06 |
| ale_darcy 25 | enthalpy form | 0.000126 | 3.03e-8 | 0.011 | 2.22 | 0.000121 | 1.0e-6 | 0.000121 | 0.0243 | 2.22 |
| ale_darcy 100 | temperature form | 5.96e-5 | 3.0e-8 | 0.0423 | 0.994 | 0.000232 | 3.0e-7 | 1.44e-5 | 0.0564 | 1.33 |
| ale_darcy 100 | enthalpy form | 5.96e-5 | 3.0e-8 | 0.0423 | 2.5 | 7.77e-5 | 3.0e-7 | 2.74e-5 | 0.142 | 3.36 |

## The original PMMA decks (solid λ = 0) at the default tolerances

No tolerance converges these decks (see the script's header); warm times only, and the ratio to the old core's.

| deck | run | retcode | steps | warm [s] | ratio |
|---|---|---|---|---|---|
| eulerian 25 | old core | Success | 292 | 0.0288 | 1.0 |
| eulerian 25 | new core, temperature form | Success | 239 | 0.0168 | 0.583 |
| eulerian 25 | new core, enthalpy form | Success | 272 | 0.0479 | 1.67 |
| eulerian 100 | old core | Success | 941 | 0.288 | 1.0 |
| eulerian 100 | new core, temperature form | Success | 989 | 0.229 | 0.797 |
| eulerian 100 | new core, enthalpy form | Success | 1037 | 0.679 | 2.36 |
| ale 25 | old core | Success | 196 | 0.0229 | 1.0 |
| ale 25 | new core, temperature form | Success | 254 | 0.0279 | 1.22 |
| ale 25 | new core, enthalpy form | Success | 260 | 0.0557 | 2.44 |
| ale 100 | old core | Success | 1660 | 0.551 | 1.0 |
| ale 100 | new core, temperature form | Success | 1004 | 0.381 | 0.692 |
| ale 100 | new core, enthalpy form | Success | 1099 | 0.904 | 1.64 |

## Ladders

| deck | run | reltol | retcode | error [K] | steps | warm [s] |
|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000329 | 23 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000159 | 28 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 3.42e-5 | 34 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 6.43e-5 | 42 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 9.9e-6 | 52 | 0.00418 |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 1.79e-5 | 64 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 3.77e-6 | 80 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 9.99e-7 | 103 |  |
| eulerian 25, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 9.73e-7 | 128 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.00108 | 25 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00121 | 29 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000205 | 33 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000194 | 40 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 5.09e-5 | 48 | 0.00312 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 1.72e-5 | 64 |  |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 2.95e-6 | 84 | 0.00484 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 7.2e-7 | 113 | 0.00622 |
| eulerian 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 7.43e-7 | 153 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.00183 | 33 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.000373 | 38 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.000763 | 44 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.00023 | 54 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 9.17e-5 | 66 | 0.0105 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 1.77e-5 | 84 |  |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 3.69e-6 | 114 | 0.0177 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 1.37e-6 | 164 | 0.0251 |
| eulerian 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 2.19e-7 | 229 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000204 | 27 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000107 | 32 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 3.52e-5 | 38 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 1.4e-5 | 47 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 1.37e-5 | 56 | 0.0132 |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 8.36e-6 | 69 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 2.0e-6 | 85 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 1.66e-6 | 106 |  |
| eulerian 100, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 5.8e-7 | 133 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0115 | 31 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00142 | 36 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000719 | 42 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000305 | 49 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 0.000297 | 58 | 0.011 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 2.83e-5 | 75 |  |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 9.77e-6 | 98 | 0.0174 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 3.21e-6 | 129 | 0.0224 |
| eulerian 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 1.08e-6 | 165 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.00665 | 42 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00184 | 49 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.000808 | 57 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.000344 | 70 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 6.12e-5 | 85 | 0.0519 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 1.39e-5 | 113 |  |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 7.36e-6 | 168 | 0.103 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 1.61e-6 | 291 | 0.181 |
| eulerian 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 2.48e-6 | 492 |  |
| ale 25, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000553 | 23 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000236 | 27 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 9.82e-5 | 33 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 5.35e-5 | 41 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 2.31e-5 | 51 | 0.00455 |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 4.91e-6 | 63 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 2.46e-6 | 78 |  |
| ale 25, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 8.29e-7 | 99 |  |
| ale 25, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 3.26e-7 | 124 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0027 | 25 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.000464 | 29 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.000703 | 33 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000124 | 39 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 6.82e-5 | 46 | 0.00434 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 2.76e-5 | 59 |  |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 3.09e-6 | 79 | 0.00675 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 5.26e-7 | 108 | 0.00916 |
| ale 25, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 2.01e-7 | 141 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.000539 | 32 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00149 | 36 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 8.48e-5 | 43 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 5.1e-5 | 52 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 0.000108 | 63 | 0.0132 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 1.44e-5 | 80 | 0.0162 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 1.54e-6 | 107 | 0.0213 |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 7.68e-7 | 153 |  |
| ale 25, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 1.88e-7 | 212 |  |
| ale 100, solid λ 1.0e-6 | old core | 0.0001 | Success | 0.000608 | 28 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-5 | Success | 0.000199 | 32 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-5 | Success | 0.000106 | 37 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-6 | Success | 3.09e-5 | 46 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-6 | Success | 1.03e-5 | 54 | 0.0149 |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-7 | Success | 7.31e-6 | 67 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-7 | Success | 4.43e-6 | 81 |  |
| ale 100, solid λ 1.0e-6 | old core | 3.0e-8 | Success | 9.1e-7 | 104 |  |
| ale 100, solid λ 1.0e-6 | old core | 1.0e-8 | Success | 5.61e-7 | 129 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 0.0001 | Success | 0.0159 | 30 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-5 | Success | 0.00144 | 34 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-5 | Success | 0.00109 | 39 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-6 | Success | 0.000164 | 47 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-6 | Success | 9.66e-5 | 54 | 0.0178 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-7 | Success | 1.18e-5 | 70 |  |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-7 | Success | 6.54e-6 | 90 | 0.0281 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 3.0e-8 | Success | 5.45e-7 | 117 | 0.0356 |
| ale 100, solid λ 1.0e-6 | new core, temperature form | 1.0e-8 | Success | 3.45e-7 | 153 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 0.0001 | Success | 0.0449 | 40 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-5 | Success | 0.00257 | 47 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-5 | Success | 0.00155 | 54 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-6 | Success | 0.000439 | 66 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-6 | Success | 0.000188 | 80 | 0.0605 |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-7 | Success | 5.83e-5 | 105 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-7 | Success | 1.28e-5 | 150 |  |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 3.0e-8 | Success | 1.84e-6 | 251 | 0.194 |
| ale 100, solid λ 1.0e-6 | new core, enthalpy form | 1.0e-8 | Success | 1.94e-7 | 426 | 0.341 |
| ale_darcy 25 | old core | 0.0001 | Success | 0.000278 | 79 |  |
| ale_darcy 25 | old core | 3.0e-5 | Success | 0.000131 | 59 |  |
| ale_darcy 25 | old core | 1.0e-5 | Success | 0.000132 | 63 |  |
| ale_darcy 25 | old core | 3.0e-6 | Success | 0.000154 | 69 |  |
| ale_darcy 25 | old core | 1.0e-6 | Success | 0.000126 | 75 | 0.011 |
| ale_darcy 25 | old core | 3.0e-7 | Success | 3.74e-5 | 90 |  |
| ale_darcy 25 | old core | 1.0e-7 | Success | 3.1e-5 | 107 |  |
| ale_darcy 25 | old core | 3.0e-8 | Success | 5.31e-5 | 122 |  |
| ale_darcy 25 | old core | 1.0e-8 | Success | 6.22e-6 | 166 |  |
| ale_darcy 25 | new core, temperature form | 0.0001 | Success | 0.00788 | 45 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-5 | Success | 0.00184 | 54 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-5 | Success | 0.000594 | 62 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-6 | Success | 0.000521 | 76 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-6 | Success | 7.28e-5 | 100 | 0.0116 |
| ale_darcy 25 | new core, temperature form | 3.0e-7 | Success | 5.27e-5 | 135 | 0.0152 |
| ale_darcy 25 | new core, temperature form | 1.0e-7 | Success | 6.72e-6 | 186 |  |
| ale_darcy 25 | new core, temperature form | 3.0e-8 | Success | 1.14e-6 | 259 |  |
| ale_darcy 25 | new core, temperature form | 1.0e-8 | Success | 4.47e-7 | 353 |  |
| ale_darcy 25 | new core, enthalpy form | 0.0001 | Success | 0.0132 | 50 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-5 | Success | 0.00176 | 61 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-5 | Success | 0.000939 | 75 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-6 | Success | 0.000154 | 102 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-6 | Success | 0.000121 | 128 | 0.0243 |
| ale_darcy 25 | new core, enthalpy form | 3.0e-7 | Success | 2.22e-5 | 179 | 0.033 |
| ale_darcy 25 | new core, enthalpy form | 1.0e-7 | Success | 9.67e-6 | 234 |  |
| ale_darcy 25 | new core, enthalpy form | 3.0e-8 | Success | 7.69e-7 | 337 |  |
| ale_darcy 25 | new core, enthalpy form | 1.0e-8 | Success | 7.09e-7 | 480 |  |
| ale_darcy 100 | old core | 0.0001 | Success | 0.00022 | 589 |  |
| ale_darcy 100 | old core | 3.0e-5 | Success | 0.000287 | 180 |  |
| ale_darcy 100 | old core | 1.0e-5 | Success | 0.000124 | 90 |  |
| ale_darcy 100 | old core | 3.0e-6 | Success | 0.000162 | 74 |  |
| ale_darcy 100 | old core | 1.0e-6 | Success | 5.96e-5 | 80 | 0.0423 |
| ale_darcy 100 | old core | 3.0e-7 | Success | 0.000255 | 92 |  |
| ale_darcy 100 | old core | 1.0e-7 | Success | 8.63e-5 | 98 |  |
| ale_darcy 100 | old core | 3.0e-8 | Success | 4.38e-5 | 126 |  |
| ale_darcy 100 | old core | 1.0e-8 | Success | 9.52e-6 | 169 |  |
| ale_darcy 100 | new core, temperature form | 0.0001 | Success | 0.0116 | 48 |  |
| ale_darcy 100 | new core, temperature form | 3.0e-5 | Success | 0.00192 | 56 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-5 | Success | 0.00175 | 66 |  |
| ale_darcy 100 | new core, temperature form | 3.0e-6 | Success | 0.000225 | 80 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-6 | Success | 0.000232 | 102 | 0.0421 |
| ale_darcy 100 | new core, temperature form | 3.0e-7 | Success | 1.44e-5 | 140 | 0.0564 |
| ale_darcy 100 | new core, temperature form | 1.0e-7 | Success | 5.48e-6 | 188 | 0.0757 |
| ale_darcy 100 | new core, temperature form | 3.0e-8 | Success | 4.67e-6 | 260 |  |
| ale_darcy 100 | new core, temperature form | 1.0e-8 | Success | 4.61e-7 | 350 |  |
| ale_darcy 100 | new core, enthalpy form | 0.0001 | Success | 0.01 | 55 |  |
| ale_darcy 100 | new core, enthalpy form | 3.0e-5 | Success | 0.0054 | 65 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-5 | Success | 0.00225 | 78 |  |
| ale_darcy 100 | new core, enthalpy form | 3.0e-6 | Success | 0.000718 | 96 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-6 | Success | 7.77e-5 | 126 | 0.106 |
| ale_darcy 100 | new core, enthalpy form | 3.0e-7 | Success | 2.74e-5 | 170 | 0.142 |
| ale_darcy 100 | new core, enthalpy form | 1.0e-7 | Success | 7.92e-6 | 250 | 0.21 |
| ale_darcy 100 | new core, enthalpy form | 3.0e-8 | Success | 2.14e-6 | 411 |  |
| ale_darcy 100 | new core, enthalpy form | 1.0e-8 | Success | 3.98e-7 | 641 |  |
