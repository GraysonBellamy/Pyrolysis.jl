# The end-to-end comparison of the two cores

Produced by `cutover_comparison.jl` at commit b8950cc (tree aac0390). Julia 1.13.0, one thread per process. Accuracy: 4 processes at a one-minute load of 2.6, 2.6, 2.6, 2.6 at their starts. Timing: one process on CPUs `1`, one-minute load 1.5 at the start and 1.7 at the end, 623 s. Memory: one process per core and row. The script's header defines every row, measure and budget.

## The rows

Window: the saved times before the first event of any solve of the row [s] (`Inf`: none). Target: the old default's window RMS against its own reltol 1e-9 solve [K; relative thickness on V-BO-04]; matched: the finite-volume core's loosest rung at or below it whose own closure is within its budget (on the window alone: the loosest at or below it). Accuracy: every verdict of the common quantities (below). Timing: warm medians [s]; ratio: matched over the old default.

| row | cells | reltol | window [s] | target | spreads old, new | matched reltol (on the window alone) | accuracy | old [s] | new matched [s] | ratio | equal-tolerance ratio | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | 25 | 1.0e-6 | Inf | 8.89e-5 | 1.2e-6, 6.65e-9 | 1.0e-5 (1.0e-5) | pass | 0.00372 | 0.00264 | 0.71 | 0.911 | pass |
| eulerian 100, solid λ 1.0e-6 | 100 | 1.0e-6 | Inf | 8.42e-5 | 1.35e-6, 1.59e-8 | 3.0e-6 (3.0e-6) | pass | 0.0133 | 0.0111 | 0.83 | 1.0 | pass |
| ale 25, solid λ 1.0e-6 | 25 | 1.0e-6 | Inf | 8.49e-5 | 1.36e-6, 6.39e-9 | 1.0e-5 (1.0e-5) | pass | 0.00444 | 0.00437 | 0.984 | 1.21 | pass |
| ale 100, solid λ 1.0e-6 | 100 | 1.0e-6 | Inf | 7.41e-5 | 1.74e-6, 5.63e-9 | 3.0e-6 (3.0e-6) | pass | 0.0142 | 0.0186 | 1.31 | 1.48 | pass |
| ale_darcy 25 | 25 | 1.0e-6 | Inf | 7.1e-5 | 1.41e-6, 1.81e-8 | 1.0e-6 (1.0e-6) | pass | 0.0101 | 0.0113 | 1.12 | 1.12 | pass |
| ale_darcy 100 | 100 | 1.0e-6 | Inf | 6.69e-5 | 1.57e-6, 1.76e-8 | 1.0e-6 (1.0e-6) | pass | 0.0385 | 0.0443 | 1.15 | 1.15 | pass |
| V-BO-01 | 182 | 1.0e-6 | 167.8 | 0.00866 | 0.000323, 1.39e-6 | 1.0e-5 (1.0e-5) | pass | 0.535 | 0.66 | 1.23 | 1.76 | pass |
| V-BO-02 | 182 | 1.0e-6 | 62.42 | 0.015 | 0.000638, 1.01e-6 | 1.0e-5 (1.0e-5) | **fail** | 0.934 | 0.575 | 0.616 | 0.905 | not read (accuracy) |
| V-BO-03 | 182 | 1.0e-6 | 167.9 | 0.00866 | 0.000341, 7.6e-7 | 1.0e-5 (1.0e-5) | pass | 0.408 | 0.627 | 1.54 | 2.22 | owner |
| V-BO-04 | 40 | 1.0e-7 | 299.8 | 6.44e-9 | 2.9e-11, 3.71e-14 | 1.0e-5 (1.0e-5) | pass | 0.0265 | 0.0261 | 0.983 | 2.09 | pass |
| V-BO-05 | 182 | 1.0e-6 | 259.8 | 0.00651 | 0.000269, 4.59e-7 | 3.0e-6 (3.0e-6) | pass | 1.01 | 0.2 | 0.197 | 0.214 | pass |
| V-BO-06 | 132 | 1.0e-6 | 167.6 | 0.0123 | 0.000317, 1.05e-6 | 1.0e-5 (1.0e-5) | pass | 0.388 | 0.546 | 1.41 | 1.99 | pass |
| V-BO-07 | 182 | 1.0e-6 | 307.8 | 0.00905 | 0.000395, 9.5e-7 | 1.0e-5 (1.0e-5) | pass | 0.765 | 0.707 | 0.924 | 1.48 | pass |
| V-BO-08 | 132 | 1.0e-6 | 65.26 | 0.00766 | 0.000132, 7.33e-7 | 3.0e-6 (1.0e-5) | pass | 0.156 | 0.132 | 0.848 | 0.965 | pass |
| Douglas fir, 600 s | 240 | 1.0e-7 | Inf | 1.23e-5 | 1.71e-6, 1.43e-7 | 3.0e-7 (3.0e-7) | pass | 5.4 | 7.59 | 1.41 | 1.87 | pass |
| white pine, 30kW | 160 | 1.0e-7 | Inf | 0.000187 | 5.68e-5, 2.31e-6 | 1.0e-5 (1.0e-5) | pass | 3.0 | 4.55 | 1.52 | 2.36 | owner |
| white pine, 50kW | 160 | 1.0e-7 | Inf | 0.000135 | 6.52e-5, 1.85e-6 | 1.0e-5 (1.0e-5) | pass | 3.11 | 3.53 | 1.13 | 2.57 | pass |
| PMMA cone | 72 | 1.0e-6 | 348.5 | 0.0156 | 0.000314, 0.000799 | 1.0e-5 (1.0e-5) | pass | 0.434 | 0.991 | 2.28 | 3.14 | **over 2** |

## The common quantities

Between the old default and the finite-volume core's matched solve. Temperatures [K] at the old core's live cell centres, over the window (the whole run on a row without events) and over the whole run; the mass loss [%] and its difference in points; the burnout time [s], the main peak [g m⁻² s⁻¹] and its time [s], the tail [g/m²], each with its difference relative to the old core's value; each core's closure against its case's ledger tolerance.

| row | quantity | old | new | difference | budget | verdict |
|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | T_rms |  |  | 0.0105 | 1.0 | pass |
| eulerian 25, solid λ 1.0e-6 | T_max |  |  | 0.0615 | 5.0 | pass |
| eulerian 25, solid λ 1.0e-6 | mass_loss | 2.8159 | 2.8116 | 0.00429 | 0.3 | pass |
| eulerian 25, solid λ 1.0e-6 | closure_old |  |  | 1.32e-5 |  |  |
| eulerian 25, solid λ 1.0e-6 | closure_new |  |  | 1.1e-6 |  |  |
| eulerian 100, solid λ 1.0e-6 | T_rms |  |  | 0.000434 | 1.0 | pass |
| eulerian 100, solid λ 1.0e-6 | T_max |  |  | 0.00321 | 5.0 | pass |
| eulerian 100, solid λ 1.0e-6 | mass_loss | 2.9262 | 2.9261 | 0.000155 | 0.3 | pass |
| eulerian 100, solid λ 1.0e-6 | closure_old |  |  | 7.31e-6 |  |  |
| eulerian 100, solid λ 1.0e-6 | closure_new |  |  | 1.5e-6 |  |  |
| ale 25, solid λ 1.0e-6 | T_rms |  |  | 0.00779 | 1.0 | pass |
| ale 25, solid λ 1.0e-6 | T_max |  |  | 0.0353 | 5.0 | pass |
| ale 25, solid λ 1.0e-6 | mass_loss | 2.8494 | 2.847 | 0.0024 | 0.3 | pass |
| ale 25, solid λ 1.0e-6 | closure_old |  |  | 2.48e-5 |  |  |
| ale 25, solid λ 1.0e-6 | closure_new |  |  | 4.19e-6 |  |  |
| ale 100, solid λ 1.0e-6 | T_rms |  |  | 0.00021 | 1.0 | pass |
| ale 100, solid λ 1.0e-6 | T_max |  |  | 0.00127 | 5.0 | pass |
| ale 100, solid λ 1.0e-6 | mass_loss | 3.0036 | 3.0035 | 4.4e-5 | 0.3 | pass |
| ale 100, solid λ 1.0e-6 | closure_old |  |  | 2.82e-5 |  |  |
| ale 100, solid λ 1.0e-6 | closure_new |  |  | 4.84e-7 |  |  |
| ale_darcy 25 | T_rms |  |  | 0.105 | 1.0 | pass |
| ale_darcy 25 | T_max |  |  | 0.317 | 5.0 | pass |
| ale_darcy 25 | mass_loss | 32.775 | 32.795 | 0.0196 | 0.3 | pass |
| ale_darcy 25 | closure_old |  |  | 1.45e-5 |  |  |
| ale_darcy 25 | closure_new |  |  | 1.08e-6 |  |  |
| ale_darcy 100 | T_rms |  |  | 0.00269 | 1.0 | pass |
| ale_darcy 100 | T_max |  |  | 0.0158 | 5.0 | pass |
| ale_darcy 100 | mass_loss | 32.42 | 32.419 | 0.000921 | 0.3 | pass |
| ale_darcy 100 | closure_old |  |  | 1.69e-5 |  |  |
| ale_darcy 100 | closure_new |  |  | 1.19e-6 |  |  |
| V-BO-01 | T_rms_window |  |  | 0.0119 | 1.0 | pass |
| V-BO-01 | T_max_window |  |  | 0.592 | 5.0 | pass |
| V-BO-01 | T_rms |  |  | 0.545 | 1.0 | no common reference |
| V-BO-01 | T_max |  |  | 33.4 | 5.0 | no common reference |
| V-BO-01 | mass_loss | 69.921 | 69.921 | 3.57e-6 | 0.3 | pass |
| V-BO-01 | burnout | 364.0 | 363.5 | 0.00137 | 0.01 | pass |
| V-BO-01 | peak | 32.392 | 32.391 | 4.05e-5 | 0.03 | pass |
| V-BO-01 | t_peak | 335.0 | 335.0 | 0.0 | 0.03 | pass |
| V-BO-01 | tail | 883.01 | 885.16 | 0.00244 | 0.1 | no common reference |
| V-BO-01 | closure_old |  |  | 0.000478 | 0.001 | pass |
| V-BO-01 | closure_new |  |  | 1.59e-7 | 3.0e-5 | pass |
| V-BO-02 | T_rms_window |  |  | 0.0173 | 1.0 | pass |
| V-BO-02 | T_max_window |  |  | 0.379 | 5.0 | pass |
| V-BO-02 | T_rms |  |  | 9.59 | 1.0 | no common reference |
| V-BO-02 | T_max |  |  | 185.0 | 5.0 | no common reference |
| V-BO-02 | mass_loss | 69.921 | 69.921 | 1.25e-5 | 0.3 | pass |
| V-BO-02 | burnout | 452.0 | 464.0 | 0.0265 | 0.01 | no common reference |
| V-BO-02 | peak | 24.623 | 23.755 | 0.0352 | 0.03 | **fail** |
| V-BO-02 | t_peak | 429.5 | 438.5 | 0.021 | 0.03 | no common reference |
| V-BO-02 | tail | 507.62 | 558.94 | 0.101 | 0.1 | no common reference |
| V-BO-02 | closure_old |  |  | 2.3e-5 | 0.001 | pass |
| V-BO-02 | closure_new |  |  | 9.58e-5 | 0.0003 | pass |
| V-BO-03 | T_rms_window |  |  | 0.0119 | 1.0 | pass |
| V-BO-03 | T_max_window |  |  | 0.592 | 5.0 | pass |
| V-BO-03 | T_rms |  |  | 0.0117 | 1.0 | no common reference |
| V-BO-03 | T_max |  |  | 0.592 | 5.0 | no common reference |
| V-BO-03 | mass_loss | 68.879 | 68.879 | 0.000103 | 0.3 | no common reference |
| V-BO-03 | burnout | 421.0 | 421.0 | 0.0 | 0.01 | no common reference |
| V-BO-03 | peak | 30.387 | 30.385 | 4.43e-5 | 0.03 | no common reference |
| V-BO-03 | t_peak | 317.0 | 317.0 | 0.0 | 0.03 | no common reference |
| V-BO-03 | tail | 1367.6 | 1368.0 | 0.000282 | 0.1 | no common reference |
| V-BO-03 | closure_old |  |  | 0.000233 |  |  |
| V-BO-03 | closure_new |  |  | 1.91e-9 |  |  |
| V-BO-04 | T_rms_window |  |  | 2.85e-13 | 1.0 | pass |
| V-BO-04 | T_max_window |  |  | 1.48e-12 | 5.0 | pass |
| V-BO-04 | T_rms |  |  | 3.27e-13 | 1.0 | no common reference |
| V-BO-04 | T_max |  |  | 1.48e-12 | 5.0 | no common reference |
| V-BO-04 | mass_loss | 99.177 | 99.177 | 1.17e-5 | 0.3 | no common reference |
| V-BO-04 | closure_old |  |  | 0.00155 |  |  |
| V-BO-04 | closure_new |  |  | 0.00156 |  |  |
| V-BO-05 | T_rms_window |  |  | 0.00703 | 1.0 | pass |
| V-BO-05 | T_max_window |  |  | 0.278 | 5.0 | pass |
| V-BO-05 | T_rms |  |  | 0.00703 | 1.0 | no common reference |
| V-BO-05 | T_max |  |  | 0.278 | 5.0 | no common reference |
| V-BO-05 | mass_loss | 43.069 | 43.07 | 0.00118 | 0.3 | no common reference |
| V-BO-05 | closure_old |  |  | 1.62e-5 |  |  |
| V-BO-05 | closure_new |  |  | 3.76e-7 |  |  |
| V-BO-06 | T_rms_window |  |  | 0.0138 | 1.0 | pass |
| V-BO-06 | T_max_window |  |  | 0.537 | 5.0 | pass |
| V-BO-06 | T_rms |  |  | 0.404 | 1.0 | no common reference |
| V-BO-06 | T_max |  |  | 35.8 | 5.0 | no common reference |
| V-BO-06 | mass_loss | 69.948 | 69.948 | 4.73e-6 | 0.3 | pass |
| V-BO-06 | burnout | 364.0 | 364.0 | 0.0 | 0.01 | pass |
| V-BO-06 | peak | 32.292 | 32.288 | 0.000132 | 0.03 | pass |
| V-BO-06 | t_peak | 334.0 | 334.0 | 0.0 | 0.03 | pass |
| V-BO-06 | tail | 918.54 | 919.83 | 0.0014 | 0.1 | no common reference |
| V-BO-06 | closure_old |  |  | 0.00383 | 0.005 | pass |
| V-BO-06 | closure_new |  |  | 5.24e-7 | 3.0e-5 | pass |
| V-BO-07 | T_rms_window |  |  | 0.016 | 1.0 | pass |
| V-BO-07 | T_max_window |  |  | 0.303 | 5.0 | pass |
| V-BO-07 | T_rms |  |  | 0.853 | 1.0 | no common reference |
| V-BO-07 | T_max |  |  | 37.0 | 5.0 | no common reference |
| V-BO-07 | mass_loss | 65.338 | 65.338 | 1.43e-5 | 0.3 | pass |
| V-BO-07 | burnout | 325.0 | 324.0 | 0.00308 | 0.01 | pass |
| V-BO-07 | peak | 30.542 | 30.542 | 7.91e-7 | 0.03 | pass |
| V-BO-07 | t_peak | 287.5 | 287.5 | 0.0 | 0.03 | pass |
| V-BO-07 | tail | 938.75 | 939.05 | 0.00032 | 0.1 | no common reference |
| V-BO-07 | closure_old |  |  | 3.06e-5 | 0.001 | pass |
| V-BO-07 | closure_new |  |  | 1.05e-7 | 3.0e-5 | pass |
| V-BO-08 | T_rms_window |  |  | 0.00863 | 1.0 | pass |
| V-BO-08 | T_max_window |  |  | 0.346 | 5.0 | pass |
| V-BO-08 | T_rms |  |  | 0.562 | 1.0 | no common reference |
| V-BO-08 | T_max |  |  | 18.3 | 5.0 | no common reference |
| V-BO-08 | mass_loss | 27.924 | 27.924 | 1.25e-5 | 0.3 | pass |
| V-BO-08 | burnout | 69.0 | 69.5 | 0.00725 | 0.01 | no common reference |
| V-BO-08 | peak | 29.712 | 29.712 | 4.71e-6 | 0.03 | no common reference |
| V-BO-08 | t_peak | 47.0 | 47.0 | 0.0 | 0.03 | no common reference |
| V-BO-08 | tail | 634.38 | 634.86 | 0.000758 | 0.1 | no common reference |
| V-BO-08 | closure_old |  |  | 0.000317 | 0.001 | pass |
| V-BO-08 | closure_new |  |  | 1.19e-5 | 3.0e-5 | pass |
| Douglas fir, 600 s | T_rms |  |  | 0.12 | 1.0 | pass |
| Douglas fir, 600 s | T_max |  |  | 0.753 | 5.0 | pass |
| Douglas fir, 600 s | mass_loss | 84.131 | 84.128 | 0.00285 | 0.3 | pass |
| Douglas fir, 600 s | peak | 33.079 | 33.062 | 0.000535 | 0.03 | pass |
| Douglas fir, 600 s | t_peak | 284.0 | 284.0 | 0.0 | 0.03 | pass |
| Douglas fir, 600 s | tail | 920.22 | 921.78 | 0.00169 | 0.1 | pass |
| Douglas fir, 600 s | closure_old |  |  | 0.00254 |  |  |
| Douglas fir, 600 s | closure_new |  |  | 4.43e-8 |  |  |
| white pine, 30kW | T_rms |  |  | 0.0293 | 1.0 | pass |
| white pine, 30kW | T_max |  |  | 0.0897 | 5.0 | pass |
| white pine, 30kW | mass_loss | 72.756 | 72.73 | 0.026 | 0.3 | pass |
| white pine, 30kW | peak | 4.3241 | 4.3226 | 0.000356 | 0.03 | pass |
| white pine, 30kW | t_peak | 69.407 | 69.407 | 0.0 | 0.03 | pass |
| white pine, 30kW | tail | 2060.2 | 2059.4 | 0.000373 | 0.1 | pass |
| white pine, 30kW | closure_old |  |  | 1.91e-8 |  |  |
| white pine, 30kW | closure_new |  |  | 5.22e-8 |  |  |
| white pine, 50kW | T_rms |  |  | 0.0341 | 1.0 | pass |
| white pine, 50kW | T_max |  |  | 0.128 | 5.0 | pass |
| white pine, 50kW | mass_loss | 78.95 | 78.931 | 0.0186 | 0.3 | pass |
| white pine, 50kW | peak | 9.5342 | 9.5327 | 0.000155 | 0.03 | pass |
| white pine, 50kW | t_peak | 18.113 | 18.113 | 0.0 | 0.03 | pass |
| white pine, 50kW | tail | 2335.1 | 2334.6 | 0.00024 | 0.1 | pass |
| white pine, 50kW | closure_old |  |  | 3.19e-7 |  |  |
| white pine, 50kW | closure_new |  |  | 1.02e-6 |  |  |
| PMMA cone | T_rms_window |  |  | 0.0177 | 1.0 | pass |
| PMMA cone | T_max_window |  |  | 0.345 | 5.0 | pass |
| PMMA cone | T_rms |  |  | 0.0488 | 1.0 | no common reference |
| PMMA cone | T_max |  |  | 0.345 | 5.0 | no common reference |
| PMMA cone | mass_loss | 92.671 | 92.695 | 0.0234 | 0.3 | no common reference |
| PMMA cone | burnout | 1592.0 | 1592.0 | 0.0 | 0.01 | no common reference |
| PMMA cone | peak | 5.4843 | 5.4866 | 0.00042 | 0.03 | no common reference |
| PMMA cone | t_peak | 650.67 | 653.33 | 0.0041 | 0.03 | no common reference |
| PMMA cone | tail | 4310.5 | 4296.9 | 0.00317 | 0.1 | no common reference |
| PMMA cone | closure_old |  |  | 3.27e-6 |  |  |
| PMMA cone | closure_new |  |  | 8.27e-8 |  |  |

## Closures and case QoIs

Each core's own closure (relative): the finite-volume core's mass and energy from its conservation recorder (mass from the case's ledger on the cases), the old core's mass from the saved-time trapezoid. QoIs: passing of counted, the failing named.

| row | old: retcode, mass | old QoIs | new: outcome, mass, energy | new QoIs |
|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | Success, 1.32e-5 |  | running, 1.1e-6, 9.35e-7 |  |
| eulerian 100, solid λ 1.0e-6 | Success, 7.31e-6 |  | running, 1.5e-6, -1.98e-7 |  |
| ale 25, solid λ 1.0e-6 | Success, 2.48e-5 |  | running, -4.19e-6, 2.02e-7 |  |
| ale 100, solid λ 1.0e-6 | Success, 2.82e-5 |  | running, 4.84e-7, 9.09e-7 |  |
| ale_darcy 25 | Success, 1.45e-5 |  | running, -1.08e-6, -1.88e-6 |  |
| ale_darcy 100 | Success, 1.69e-5 |  | running, -1.19e-6, -2.23e-6 |  |
| V-BO-01 | Success, 0.000478 | 7/7 | running, 1.59e-7,  | 8/8 |
| V-BO-02 | Success, 2.3e-5 | 10/10 | running, 9.58e-5,  | 10/10 |
| V-BO-03 | Success, 0.000233 | 3/4 (transmitted-flux warning fired) | running, 1.91e-9,  | 4/4 |
| V-BO-04 | Success, 0.00155 | 4/4 | running, 0.00156,  | 5/5 |
| V-BO-05 | Terminated, 1.62e-5 | 3/4 (interior-void warning) | interior, 3.76e-7,  | 4/4 |
| V-BO-06 | Success, 0.00383 | 7/7 | running, 5.24e-7,  | 8/8 |
| V-BO-07 | Success, 3.06e-5 | 7/7 | running, 1.05e-7,  | 8/8 |
| V-BO-08 | Success, 0.000317 | 4/4 | running, 1.19e-5,  | 6/6 |
| Douglas fir, 600 s | Success, 0.00254 |  | running, 4.43e-8, -5.65e-8 |  |
| white pine, 30kW | Success, 1.91e-8 |  | running, 5.22e-8, 1.7e-6 |  |
| white pine, 50kW | Success, 3.19e-7 |  | running, 1.02e-6, -6.89e-6 |  |
| PMMA cone | Success, 3.27e-6 |  | running, -8.27e-8, -0.000105 |  |

## The counters

Of the timed solves: accepted (rejected) steps, right-hand sides, Jacobians, linear solves, nonlinear failures, domain rejections (the old core counts none), events and invocations; the old default and the finite-volume core at its matched and default rungs.

| row | solve | reltol | warm [s] (range, n) | steps (rejected) | RHS | Jacobians | linear solves | nonlinear failures | domain rejections | events, invocations |
|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | old default | 1.0e-6 | 0.00372 (0.00357–0.00396, 9) | 52 (3) | 1311 | 10 | 1314 | 0 |  | 0, 0 |
| eulerian 25, solid λ 1.0e-6 | matched | 1.0e-5 | 0.00264 (0.00259–0.00272, 3) | 33 (0) | 297 | 33 | 264 | 0 | 0 | 0, 0 |
| eulerian 25, solid λ 1.0e-6 | new default | 1.0e-6 | 0.00339 (0.00325–0.00375, 3) | 48 (1) | 440 | 48 | 392 | 0 | 0 | 0, 0 |
| eulerian 100, solid λ 1.0e-6 | old default | 1.0e-6 | 0.0133 (0.0128–0.0137, 9) | 56 (3) | 1362 | 12 | 1365 | 1 |  | 0, 0 |
| eulerian 100, solid λ 1.0e-6 | matched | 3.0e-6 | 0.0111 (0.011–0.0115, 3) | 49 (2) | 457 | 49 | 408 | 0 | 0 | 0, 0 |
| eulerian 100, solid λ 1.0e-6 | new default | 1.0e-6 | 0.0133 (0.0132–0.0137, 3) | 58 (2) | 538 | 58 | 480 | 0 | 0 | 0, 0 |
| ale 25, solid λ 1.0e-6 | old default | 1.0e-6 | 0.00444 (0.00419–0.00607, 9) | 51 (2) | 1337 | 12 | 1339 | 0 |  | 0, 0 |
| ale 25, solid λ 1.0e-6 | matched | 1.0e-5 | 0.00437 (0.00432–0.00455, 3) | 33 (2) | 313 | 33 | 280 | 0 | 0 | 0, 0 |
| ale 25, solid λ 1.0e-6 | new default | 1.0e-6 | 0.00536 (0.00511–0.00556, 3) | 46 (1) | 422 | 46 | 376 | 0 | 0 | 0, 0 |
| ale 100, solid λ 1.0e-6 | old default | 1.0e-6 | 0.0142 (0.0135–0.015, 9) | 54 (3) | 1283 | 15 | 1286 | 0 |  | 0, 0 |
| ale 100, solid λ 1.0e-6 | matched | 3.0e-6 | 0.0186 (0.0179–0.019, 3) | 47 (1) | 431 | 47 | 384 | 0 | 0 | 0, 0 |
| ale 100, solid λ 1.0e-6 | new default | 1.0e-6 | 0.021 (0.0209–0.021, 3) | 54 (0) | 486 | 54 | 432 | 0 | 0 | 0, 0 |
| ale_darcy 25 | old default | 1.0e-6 | 0.0101 (0.00962–0.0104, 9) | 75 (2) | 2428 | 49 | 2430 | 2 |  | 0, 0 |
| ale_darcy 25 | matched | 1.0e-6 | 0.0113 (0.011–0.0119, 3) | 100 (3) | 924 | 100 | 824 | 0 | 0 | 0, 0 |
| ale_darcy 25 | new default | 1.0e-6 | 0.0113 (0.011–0.0119, 3) | 100 (3) | 924 | 100 | 824 | 0 | 0 | 0, 0 |
| ale_darcy 100 | old default | 1.0e-6 | 0.0385 (0.0382–0.0399, 9) | 80 (2) | 2726 | 51 | 2728 | 2 |  | 0, 0 |
| ale_darcy 100 | matched | 1.0e-6 | 0.0443 (0.0425–0.0465, 3) | 102 (2) | 934 | 102 | 832 | 0 | 0 | 0, 0 |
| ale_darcy 100 | new default | 1.0e-6 | 0.0443 (0.0425–0.0465, 3) | 102 (2) | 934 | 102 | 832 | 0 | 0 | 0, 0 |
| V-BO-01 | old default | 1.0e-6 | 0.535 (0.527–0.549, 9) | 853 (3) | 19076 | 164 | 19079 | 8 |  | 60, 60 |
| V-BO-01 | matched | 1.0e-5 | 0.66 (0.657–0.666, 3) | 647 (200) | 7427 | 649 | 6776 | 0 | 2 | 60, 1 |
| V-BO-01 | new default | 1.0e-6 | 0.943 (0.937–0.945, 3) | 985 (242) | 10803 | 985 | 9816 | 0 | 0 | 60, 1 |
| V-BO-02 | old default | 1.0e-6 | 0.934 (0.924–0.953, 9) | 1572 (0) | 38574 | 341 | 38574 | 56 |  | 60, 60 |
| V-BO-02 | matched | 1.0e-5 | 0.575 (0.572–0.593, 3) | 1355 (99) | 13107 | 1355 | 11632 | 0 | 0 | 60, 60 |
| V-BO-02 | new default | 1.0e-6 | 0.845 (0.843–0.883, 3) | 1978 (288) | 20226 | 1978 | 18128 | 0 | 0 | 60, 60 |
| V-BO-03 | old default | 1.0e-6 | 0.408 (0.405–0.418, 9) | 681 (3) | 14538 | 126 | 14541 | 2 |  | 48, 48 |
| V-BO-03 | matched | 1.0e-5 | 0.627 (0.626–0.629, 3) | 601 (195) | 6972 | 604 | 6368 | 0 | 3 | 0, 0 |
| V-BO-03 | new default | 1.0e-6 | 0.904 (0.9–0.908, 3) | 928 (231) | 10200 | 928 | 9272 | 0 | 0 | 0, 0 |
| V-BO-04 | old default | 1.0e-7 | 0.0265 (0.0257–0.0272, 9) | 353 (0) | 7309 | 77 | 7309 | 13 |  | 35, 35 |
| V-BO-04 | matched | 1.0e-5 | 0.0261 (0.0258–0.0274, 3) | 124 (0) | 1118 | 124 | 992 | 0 | 0 | 39, 1 |
| V-BO-04 | new default | 1.0e-7 | 0.0555 (0.0552–0.0557, 3) | 322 (0) | 2900 | 322 | 2576 | 0 | 0 | 39, 1 |
| V-BO-05 | old default | 1.0e-6 | 1.01 (1.0–1.03, 9) | 169 (9) | 5824 | 86 | 5833 | 8 |  | 0, 0 |
| V-BO-05 | matched | 3.0e-6 | 0.2 (0.199–0.203, 3) | 294 (185) | 4126 | 294 | 3832 | 0 | 0 | 0, 0 |
| V-BO-05 | new default | 1.0e-6 | 0.217 (0.217–0.225, 3) | 347 (143) | 4267 | 347 | 3920 | 0 | 0 | 0, 0 |
| V-BO-06 | old default | 1.0e-6 | 0.388 (0.385–0.392, 9) | 896 (0) | 19556 | 169 | 19556 | 9 |  | 60, 60 |
| V-BO-06 | matched | 1.0e-5 | 0.546 (0.541–0.552, 3) | 684 (234) | 8032 | 686 | 7344 | 0 | 2 | 60, 1 |
| V-BO-06 | new default | 1.0e-6 | 0.772 (0.772–0.78, 3) | 1040 (236) | 11250 | 1040 | 10208 | 0 | 0 | 60, 1 |
| V-BO-07 | old default | 1.0e-6 | 0.765 (0.762–0.786, 9) | 794 (0) | 16962 | 181 | 16962 | 10 |  | 60, 60 |
| V-BO-07 | matched | 1.0e-5 | 0.707 (0.707–0.71, 3) | 466 (132) | 5253 | 467 | 4784 | 0 | 1 | 60, 1 |
| V-BO-07 | new default | 1.0e-6 | 1.13 (1.13–1.15, 3) | 785 (227) | 8883 | 785 | 8096 | 0 | 0 | 60, 1 |
| V-BO-08 | old default | 1.0e-6 | 0.156 (0.154–0.158, 9) | 291 (1) | 7355 | 65 | 7356 | 7 |  | 10, 10 |
| V-BO-08 | matched | 3.0e-6 | 0.132 (0.131–0.132, 3) | 187 (54) | 2117 | 187 | 1928 | 0 | 0 | 10, 1 |
| V-BO-08 | new default | 1.0e-6 | 0.151 (0.149–0.151, 3) | 221 (51) | 2399 | 221 | 2176 | 0 | 0 | 10, 1 |
| Douglas fir, 600 s | old default | 1.0e-7 | 5.4 (5.31–5.42, 3) | 621 (11) | 17502 | 141 | 17513 | 0 |  | 0, 0 |
| Douglas fir, 600 s | matched | 3.0e-7 | 7.59 (7.59–7.61, 3) | 704 (49) | 6728 | 704 | 6024 | 0 | 0 | 0, 0 |
| Douglas fir, 600 s | new default | 1.0e-7 | 10.1 (10.1–10.1, 3) | 936 (78) | 9048 | 936 | 8112 | 0 | 0 | 0, 0 |
| white pine, 30kW | old default | 1.0e-7 | 3.0 (2.94–3.01, 3) | 802 (2) | 16184 | 80 | 16186 | 0 |  | 0, 0 |
| white pine, 30kW | matched | 1.0e-5 | 4.55 (4.53–4.55, 3) | 744 (0) | 6696 | 744 | 5952 | 0 | 0 | 0, 0 |
| white pine, 30kW | new default | 1.0e-7 | 7.09 (6.86–7.54, 3) | 1048 (69) | 9984 | 1048 | 8936 | 0 | 0 | 0, 0 |
| white pine, 50kW | old default | 1.0e-7 | 3.11 (3.07–3.16, 9) | 679 (0) | 14775 | 136 | 14775 | 1 |  | 0, 0 |
| white pine, 50kW | matched | 1.0e-5 | 3.53 (3.45–3.54, 3) | 528 (2) | 4768 | 528 | 4240 | 0 | 0 | 0, 0 |
| white pine, 50kW | new default | 1.0e-7 | 8.01 (7.98–8.14, 3) | 1167 (108) | 11367 | 1167 | 10200 | 0 | 0 | 0, 0 |
| PMMA cone | old default | 1.0e-6 | 0.434 (0.425–0.436, 9) | 2026 (12) | 46103 | 408 | 46115 | 59 |  | 58, 58 |
| PMMA cone | matched | 1.0e-5 | 0.991 (0.981–1.01, 3) | 2351 (621) | 26138 | 2362 | 23776 | 0 | 11 | 0, 0 |
| PMMA cone | new default | 1.0e-6 | 1.36 (1.36–1.43, 3) | 3312 (898) | 36996 | 3316 | 33680 | 0 | 4 | 0, 0 |

## Memory

Each core on each row in a process of its own after one warm-up solve, at the timing's tolerance: the measured solve's allocated bytes [MiB] and the process's peak resident set [MiB], and their ratios.

| row | old allocated | new allocated | ratio | old peak RSS | new peak RSS | ratio |
|---|---|---|---|---|---|---|
| eulerian 25, solid λ 1.0e-6 | 0.907 | 0.866 | 0.955 | 974.0 | 1240.0 | 1.27 |
| eulerian 100, solid λ 1.0e-6 | 3.77 | 3.45 | 0.914 | 974.0 | 1220.0 | 1.25 |
| ale 25, solid λ 1.0e-6 | 1.02 | 1.55 | 1.52 | 974.0 | 1220.0 | 1.25 |
| ale 100, solid λ 1.0e-6 | 4.57 | 6.84 | 1.5 | 975.0 | 1260.0 | 1.3 |
| ale_darcy 25 | 1.04 | 1.96 | 1.89 | 974.0 | 1230.0 | 1.27 |
| ale_darcy 100 | 4.62 | 8.46 | 1.83 | 976.0 | 1230.0 | 1.26 |
| V-BO-01 | 76.2 | 79.0 | 1.04 | 1150.0 | 1240.0 | 1.08 |
| V-BO-02 | 71.4 | 32.1 | 0.45 | 1140.0 | 1260.0 | 1.1 |
| V-BO-03 | 77.3 | 76.2 | 0.986 | 1160.0 | 1330.0 | 1.15 |
| V-BO-04 | 3.29 | 5.42 | 1.65 | 1090.0 | 1260.0 | 1.15 |
| V-BO-05 | 59.0 | 30.8 | 0.522 | 1160.0 | 1250.0 | 1.08 |
| V-BO-06 | 51.7 | 60.2 | 1.16 | 1140.0 | 1250.0 | 1.1 |
| V-BO-07 | 151.0 | 86.9 | 0.575 | 1140.0 | 1250.0 | 1.09 |
| V-BO-08 | 20.7 | 22.0 | 1.06 | 1160.0 | 1330.0 | 1.15 |
| Douglas fir, 600 s | 2540.0 | 1320.0 | 0.52 | 1350.0 | 1530.0 | 1.13 |
| white pine, 30kW | 1990.0 | 316.0 | 0.158 | 1450.0 | 1430.0 | 0.992 |
| white pine, 50kW | 2140.0 | 368.0 | 0.172 | 1310.0 | 1470.0 | 1.13 |
| PMMA cone | 23.7 | 47.2 | 1.99 | 1120.0 | 1310.0 | 1.17 |

## The white-pine objective

The data term of `loss` at the fitted stage-2 conductivities over the 30 and 50 kW/m² scenarios: as fitted (the old core's `KenCarp4` with Sparspak; `loss` itself gives 0.667351), the old core at its default, the finite-volume core at its default and matched rungs; the back-face temperature of the finite-volume default against the fitted solve over the data times [K].

| as fitted | old default | new default | new matched | drift (new default − as fitted) |
|---|---|---|---|---|
| 0.667351 | 0.667351 | 0.666991 | 0.667005 | -0.00036 |

Back-face temperature, 30kW: RMS 0.0283 K, largest 0.0454 K.

Back-face temperature, 50kW: RMS 0.0513 K, largest 0.069 K.

## The recheck configurations

The finite-volume core at the deck's tolerance, warm median, range and count, and allocated bytes [MiB].

| configuration | warm [s] | range [s] | n | allocated [MiB] | steps |
|---|---|---|---|---|---|
| V-BO-01 | 0.967 | 0.958–0.989 | 5 | 91.0 | 985 |
| ale 100, solid λ 1.0e-6 | 0.0221 | 0.0213–0.0229 | 5 | 7.25 | 54 |
| Douglas fir, 60 s | 3.98 | 3.95–4.05 | 5 | 789.0 | 313 |
