# Gate items 1, 2, 4, 6 and 10 at M4: the burnout ladder

Produced by `test/newcore/gates/burnout_ladder.jl` at commit 675b4f0. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [15.5, 15.8, 16.2] at the start and [16.0, 15.8, 16.3] at the end. The finite-volume core with the equidistribution correction (band 0.1, activation 0.3/0.1, `β = 0`), the PMMA–Kaowool interface Lagrangian; τ [s] is the layer's equidistribution time, with the prototype's per-cell `τ_p = τ/N²` [s] in parentheses; switch: the reactivity switch's width in ε_m (0: none). The old core runs each case at the same resolution with its own keywords (its threshold is a thickness fraction).

## Outcomes and QoIs (items 1, 2 and 4)

Invocations: per invocation its cells (below ε_m, taken by the band; the largest merged mass fraction). Item 1 passes with at most one invocation per layer and the largest fraction at most 5 ε_m. QoIs: the case's own, in assert mode against its fixed references; the failing ones are named. Curve: the largest difference of the mass loss rate from the old core's at the saved times from the old core's peak on, over that peak (at most 3 %), and the time of that difference [s]. Inadmissible: results of a deactivation the solve refused.

| case | cells | form | ε_m | τ (τ_p) | switch | corrected | retcode, outcome, end [s] | invocations | item 1 | QoIs | failing | peak, t_peak, t_bo (old core) | curve | events | inadmissible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 30 | temperature | 0.003 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 30 (6 below, 24 band; 0.525 %) | pass | 8/8 |  | 32.248, 333.0, 364.0 (32.254, 333.0, 364.0) | 38.0 % (364.0) | 30 | 0 |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 30 (4 below, 26 band; 3.06 %) | pass | 6/8 | peak MLR; time of peak MLR | 35.419, 359.5, 363.0 (32.254, 333.0, 364.0) | 36.1 % (365.0) | 30 | 0 |
| V-BO-01 | 30 | temperature | 0.03 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 17 (1 below, 16 band; 9.03 %); 13 (3 below, 10 band; 4.16 %) | fail | 5/8 | peak MLR; time of peak MLR; one invocation for the PMMA layer | 38.66, 360.5, 363.5 (32.254, 333.0, 364.0) | 35.8 % (365.0) | 30 | 0 |
| V-BO-01 | 60 | temperature | 0.003 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 60 (5 below, 55 band; 0.586 %) | pass | 8/8 |  | 32.391, 335.0, 363.5 (32.392, 335.0, 363.5) | 24.4 % (363.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 60 (6 below, 54 band; 3.38 %) | pass | 8/8 |  | 32.391, 335.0, 363.5 (32.392, 335.0, 363.5) | 17.0 % (365.0) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.03 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 35 (1 below, 34 band; 9.95 %); 25 (6 below, 19 band; 4.17 %) | fail | 7/8 | one invocation for the PMMA layer | 32.393, 335.0, 363.5 (32.392, 335.0, 363.5) | 14.3 % (365.0) | 60 | 0 |
| V-BO-01 | 120 | temperature | 0.003 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 120 (8 below, 112 band; 0.597 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 7.19 % (363.5) | 120 | 0 |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 120 (14 below, 106 band; 3.36 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 8.77 % (364.5) | 120 | 0 |
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 69 (1 below, 68 band; 9.83 %); 51 (14 below, 37 band; 4.1 %) | fail | 7/8 | one invocation for the PMMA layer | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 12.0 % (360.5) | 120 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | 1.5 | yes | Success, running, 500.0 | 60 (1 below, 59 band; 1.18 %) | pass | 7/8 | time of peak MLR | 33.007, 361.5, 363.5 (32.392, 335.0, 363.5) | 23.6 % (361.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | 1.5 | yes | Success, running, 500.0 | 34 (1 below, 33 band; 9.32 %); 26 (3 below, 23 band; 1.26 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 18.8 % (363.5) | 60 | 0 |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | 0.0 | yes | Success, running, 500.0 | 30 (4 below, 26 band; 3.46 %) | pass | 6/8 | peak MLR; time of peak MLR | 35.109, 359.0, 363.0 (32.254, 333.0, 364.0) | 38.1 % (364.5) | 30 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | 0.0 | yes | Success, running, 500.0 | 60 (5 below, 55 band; 3.83 %) | pass | 8/8 |  | 32.391, 335.0, 363.0 (32.392, 335.0, 363.5) | 17.3 % (365.0) | 60 | 0 |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | 0.0 | yes | Success, running, 500.0 | 120 (12 below, 108 band; 3.85 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 10.2 % (364.5) | 120 | 0 |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 0.119 % (359.0) | 0 | 0 |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 0.0315 % (356.5) | 0 | 0 |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 0.0256 % (377.0) | 0 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | 1.5 | no | Success, running, 480.0 | 39 (39 below, 0 band; 0.98 %) | pass | 5/5 |  |  |  % () | 39 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | 1.5 | yes | Success, running, 480.0 | 39 (39 below, 0 band; 0.98 %) | pass | 5/5 |  |  |  % () | 39 | 0 |
| V-BO-01 | 30 | enthalpy | 0.003 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 30 (4 below, 26 band; 0.557 %) | pass | 8/8 |  | 32.248, 333.0, 364.0 (32.254, 333.0, 364.0) | 33.0 % (364.0) | 30 | 0 |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 30 (3 below, 27 band; 3.2 %) | pass | 6/8 | peak MLR; time of peak MLR | 34.514, 359.5, 363.0 (32.254, 333.0, 364.0) | 36.1 % (365.0) | 30 | 0 |
| V-BO-01 | 30 | enthalpy | 0.03 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 17 (1 below, 16 band; 9.04 %); 13 (2 below, 11 band; 4.27 %) | fail | 5/8 | peak MLR; time of peak MLR; one invocation for the PMMA layer | 37.195, 360.5, 363.5 (32.254, 333.0, 364.0) | 35.9 % (365.0) | 30 | 0 |
| V-BO-01 | 60 | enthalpy | 0.003 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 60 (10 below, 50 band; 0.546 %) | pass | 8/8 |  | 32.391, 335.0, 364.0 (32.392, 335.0, 363.5) | 28.2 % (363.5) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 60 (6 below, 54 band; 3.39 %) | pass | 8/8 |  | 32.391, 335.0, 363.5 (32.392, 335.0, 363.5) | 17.1 % (365.0) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.03 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 35 (1 below, 34 band; 9.95 %); 25 (3 below, 22 band; 4.45 %) | fail | 7/8 | one invocation for the PMMA layer | 32.393, 335.0, 363.5 (32.392, 335.0, 363.5) | 15.3 % (365.0) | 60 | 0 |
| V-BO-01 | 120 | enthalpy | 0.003 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 120 (19 below, 101 band; 0.563 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 8.46 % (363.5) | 120 | 0 |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 120 (7 below, 113 band; 3.65 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 9.52 % (364.5) | 120 | 0 |
| V-BO-01 | 120 | enthalpy | 0.03 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 69 (2 below, 67 band; 9.83 %); 51 (4 below, 47 band; 4.42 %) | fail | 7/8 | one invocation for the PMMA layer | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 7.61 % (360.0) | 120 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 360 (0.1) | 1.5 | yes | Success, running, 500.0 | 60 (20 below, 40 band; 1.1 %) | pass | 8/8 |  | 32.39, 335.0, 363.5 (32.392, 335.0, 363.5) | 19.2 % (362.0) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | 1.5 | yes | Success, running, 500.0 | 34 (2 below, 32 band; 9.19 %); 26 (1 below, 25 band; 1.3 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 17.1 % (363.5) | 60 | 0 |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | 0.0 | yes | Success, running, 500.0 | 30 (3 below, 27 band; 3.61 %) | pass | 6/8 | peak MLR; time of peak MLR | 34.481, 359.0, 363.0 (32.254, 333.0, 364.0) | 39.5 % (364.5) | 30 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | 0.0 | yes | Success, running, 500.0 | 60 (6 below, 54 band; 3.8 %) | pass | 8/8 |  | 32.391, 335.0, 363.0 (32.392, 335.0, 363.5) | 17.3 % (365.0) | 60 | 0 |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | 0.0 | yes | Success, running, 500.0 | 120 (8 below, 112 band; 4.05 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 10.7 % (364.5) | 120 | 0 |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 0.113 % (359.0) | 0 | 0 |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 0.0248 % (356.5) | 0 | 0 |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | 1.5 | yes | Success, running, 500.0 | 0 | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 0.0263 % (377.0) | 0 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | 1.5 | no | Success, running, 480.0 | 39 (39 below, 0 band; 0.995 %) | pass | 5/5 |  |  |  % () | 39 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | 1.5 | yes | Success, running, 480.0 | 39 (39 below, 0 band; 0.995 %) | pass | 5/5 |  |  |  % () | 39 | 0 |

## Cells, residuals and cost (items 6 and 10)

Smallest cell: the thinnest live cell at the end of any step, before its deactivation event [m], and the smallest ratio of a live cell to its initial thickness over ε_m (item 10: at least 0.5, and no `DtLessThanMin`). C52 and Phase 1: the largest scaled constraint residual of the accepted states under C52's scale and under Phase 1's, which omits the sweep of the node corrections; Phase 1 > 1: the accepted states Phase 1's scale would have rejected. Domain: the domain test's rejections (C52's scale, admissibility included). Wall: the warm median of 3 solves without the recorders [s], the old core's at the same resolution in parentheses; same: whether those solves took the instrumented solve's steps and rejections (in the enthalpy form the recorders' residual evaluations move the warm start of the caloric inversion, so the two solves part at round-off).

| case | cells | form | ε_m | τ (τ_p) | switch | corrected | smallest cell [m] | ratio / ε_m | item 10 | C52 | Phase 1 | Phase 1 > 1 | domain | steps (rejected) | wall (old core) | same |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 30 | temperature | 0.003 | 720 (0.8) | 1.5 | yes | 5.56e-7 | 0.942 | pass | 0.369 | 0.372 | 0 | 24 | 746 (225) | 1.51 (0.724) | true |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | 1.5 | yes | 1.89e-6 | 0.96 | pass | 0.265 | 0.265 | 0 | 21 | 687 (221) | 1.03 (0.724) | true |
| V-BO-01 | 30 | temperature | 0.03 | 720 (0.8) | 1.5 | yes | 5.56e-6 | 0.942 | pass | 0.0652 | 0.0652 | 0 | 19 | 727 (262) | 1.46 (0.724) | true |
| V-BO-01 | 60 | temperature | 0.003 | 720 (0.2) | 1.5 | yes | 2.92e-7 | 0.99 | pass | 0.993 | 0.994 | 0 | 68 | 1188 (201) | 4.44 (1.76) | true |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | 1.5 | yes | 9.59e-7 | 0.975 | pass | 0.609 | 0.609 | 0 | 61 | 1126 (225) | 4.35 (1.76) | true |
| V-BO-01 | 60 | temperature | 0.03 | 720 (0.2) | 1.5 | yes | 2.69e-6 | 0.914 | pass | 0.191 | 0.192 | 0 | 44 | 1094 (227) | 4.13 (1.76) | true |
| V-BO-01 | 120 | temperature | 0.003 | 720 (0.05) | 1.5 | yes | 1.47e-7 | 0.993 | pass | 0.942 | 0.971 | 0 | 174 | 2323 (571) | 18.6 (7.27) | true |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | 1.5 | yes | 4.74e-7 | 0.965 | pass | 0.967 | 0.976 | 0 | 152 | 2120 (489) | 16.3 (7.27) | true |
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | 1.5 | yes | 1.34e-6 | 0.906 | pass | 0.716 | 0.716 | 0 | 129 | 1983 (452) | 15.7 (7.27) | true |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | 1.5 | yes | 9.82e-7 | 0.999 | pass | 0.954 | 0.955 | 0 | 71 | 1296 (239) | 4.93 (1.76) | true |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | 1.5 | yes | 9.4e-7 | 0.956 | pass | 0.194 | 0.194 | 0 | 2 | 704 (167) | 2.76 (1.76) | true |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | 0.0 | yes | 1.88e-6 | 0.954 | pass | 0.265 | 0.265 | 0 | 21 | 686 (221) | 1.37 (0.724) | true |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | 0.0 | yes | 9.57e-7 | 0.973 | pass | 0.609 | 0.609 | 0 | 61 | 1126 (224) | 4.17 (1.76) | true |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | 0.0 | yes | 4.74e-7 | 0.964 | pass | 0.967 | 0.976 | 0 | 151 | 2117 (493) | 15.7 (7.27) | true |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | 1.5 | yes | 2.24e-6 | 1.14 | pass | 0.123 | 0.124 | 0 | 21 | 621 (211) | 1.19 (0.708) | true |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | 1.5 | yes | 1.1e-6 | 1.12 | pass | 0.802 | 0.802 | 0 | 60 | 1064 (213) | 3.86 (1.48) | true |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | 1.5 | yes | 5.48e-7 | 1.11 | pass | 0.833 | 0.834 | 0 | 166 | 2110 (501) | 10.6 (5.06) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | 1.5 | no | 2.45e-6 | 0.98 | pass | 5.03e-8 | 5.03e-8 | 0 | 0 | 324 (2) | 0.104 (0.0771) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | 1.5 | yes | 2.45e-6 | 0.98 | pass | 5.04e-8 | 5.04e-8 | 0 | 0 | 324 (2) | 0.161 (0.0771) | true |
| V-BO-01 | 30 | enthalpy | 0.003 | 720 (0.8) | 1.5 | yes | 5.76e-7 | 0.977 | pass | 0.179 | 0.181 | 0 | 22 | 963 (316) | 5.01 (0.724) | false |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | 1.5 | yes | 1.92e-6 | 0.974 | pass | 0.119 | 0.119 | 0 | 21 | 919 (323) | 4.71 (0.724) | false |
| V-BO-01 | 30 | enthalpy | 0.03 | 720 (0.8) | 1.5 | yes | 5.71e-6 | 0.968 | pass | 0.0364 | 0.0364 | 0 | 15 | 925 (321) | 4.82 (0.724) | false |
| V-BO-01 | 60 | enthalpy | 0.003 | 720 (0.2) | 1.5 | yes | 2.8e-7 | 0.949 | pass | 0.909 | 0.929 | 0 | 56 | 1635 (570) | 16.8 (1.76) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | 1.5 | yes | 9.6e-7 | 0.977 | pass | 0.469 | 0.47 | 0 | 49 | 1579 (612) | 16.0 (1.76) | false |
| V-BO-01 | 60 | enthalpy | 0.03 | 720 (0.2) | 1.5 | yes | 2.89e-6 | 0.979 | pass | 0.11 | 0.11 | 0 | 43 | 1569 (596) | 15.8 (1.76) | false |
| V-BO-01 | 120 | enthalpy | 0.003 | 720 (0.05) | 1.5 | yes | 1.41e-7 | 0.959 | pass | 0.984 | 0.986 | 0 | 160 | 2947 (928) | 54.4 (7.27) | false |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | 1.5 | yes | 4.88e-7 | 0.992 | pass | 0.945 | 0.945 | 0 | 143 | 2788 (860) | 47.0 (7.27) | false |
| V-BO-01 | 120 | enthalpy | 0.03 | 720 (0.05) | 1.5 | yes | 1.45e-6 | 0.985 | pass | 0.53 | 0.53 | 0 | 112 | 2648 (854) | 54.1 (7.27) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 360 (0.1) | 1.5 | yes | 9.17e-7 | 0.932 | pass | 0.895 | 0.895 | 0 | 69 | 1770 (583) | 13.1 (1.76) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | 1.5 | yes | 9.68e-7 | 0.984 | pass | 0.142 | 0.142 | 0 | 0 | 1168 (567) | 13.1 (1.76) | false |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | 0.0 | yes | 1.91e-6 | 0.972 | pass | 0.119 | 0.119 | 0 | 22 | 919 (317) | 4.46 (0.724) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | 0.0 | yes | 9.54e-7 | 0.97 | pass | 0.469 | 0.47 | 0 | 50 | 1581 (607) | 16.1 (1.76) | false |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | 0.0 | yes | 4.84e-7 | 0.985 | pass | 0.945 | 0.945 | 0 | 143 | 2789 (860) | 51.0 (7.27) | false |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | 1.5 | yes | 2.24e-6 | 1.14 | pass | 0.0977 | 0.0977 | 0 | 22 | 872 (329) | 4.69 (0.708) | false |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | 1.5 | yes | 1.1e-6 | 1.12 | pass | 0.306 | 0.306 | 0 | 53 | 1542 (576) | 15.8 (1.48) | false |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | 1.5 | yes | 5.48e-7 | 1.11 | pass | 0.837 | 0.839 | 0 | 143 | 2666 (758) | 53.6 (5.06) | false |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | 1.5 | no | 2.49e-6 | 0.995 | pass | 2.31e-8 | 2.31e-8 | 0 | 0 | 327 (4) | 0.178 (0.0771) | true |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | 1.5 | yes | 2.49e-6 | 0.995 | pass | 2.31e-8 | 2.31e-8 | 0 | 0 | 327 (4) | 0.232 (0.0771) | true |

The script's C52 scale against `Integration.constraint_residual` at every probed state: largest relative difference 0.0.

## R25: the end of the layer (V-BO-01)

At the state that ends the step of the layer's invocation, before the merge, and after it: the cells merged, the residue `m_res` [g/m²] (after the merge, in the receiver), the mass-weighted temperature `T̄` and switch `s̄` of the merged PMMA [K], the receiver's temperature before the merge and the merged temperature `T_m` [K] with its switch, and the receiver's weight in it, `(T_m − T̄)/(T_receiver − T̄)`; the PMMA gasification rate [g m⁻² s⁻¹] `R0 = Σ mᵢ k(Tᵢ) sᵢ` before the merge, `Ra = m_res k(T̄) s̄` lumped, `Rb = m_res k(T̄)` with the switch released and `Rc = m_res k(T_m) s_m` after the merge; the mass gasified from the event to the end of the run [g/m²] (from `Rc` at the event to the first saved time after it, then the saved times); the largest mass loss rate at those saved times and its time, the old core's largest over the ten saved times from the event's.

| cells | form | ε_m | τ (τ_p) | switch | event [s] | merged | m_res (after) | T̄ | s̄ | T receiver | T_m (s_m) | weight | R0 | Ra | Rb | Rc | Rb/Ra | Rc/Rb | Rc/R0 | gasified after | largest after (t) | old core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | temperature | 0.003 | 720 (0.8) | 1.5 | 363.77 | 30 | 29.3 (29.3) | 759.27 | 0.93 | 824.38 | 784.47 (1.0) | 0.387 | 19.4 | 19.2 | 20.7 | 47.6 | 1.07 | 2.3 | 2.46 | 30.8 | 29.94 (364.0) | 17.67 |
| 30 | temperature | 0.01 | 720 (0.8) | 1.5 | 359.17 | 30 | 140.0 (140.0) | 722.55 | 0.962 | 782.19 | 729.75 (1.0) | 0.121 | 26.2 | 25.3 | 26.3 | 34.4 | 1.04 | 1.31 | 1.31 | 140.0 | 35.42 (359.5) | 29.54 |
| 30 | temperature | 0.03 | 720 (0.8) | 1.5 | 312.16 | 17 | 217.0 (217.0) | 678.95 | 0.946 | 681.79 | 678.95 (1.0) | 0.00101 | 6.8 | 6.7 | 7.08 | 7.08 | 1.06 | 1.0 | 1.04 | 1600.0 | 38.66 (360.5) | 31.51 |
| 60 | temperature | 0.003 | 720 (0.2) | 1.5 | 363.35 | 60 | 32.1 (32.1) | 758.16 | 0.95 | 793.59 | 766.13 (1.0) | 0.225 | 20.8 | 20.7 | 21.8 | 28.6 | 1.05 | 1.31 | 1.37 | 32.5 | 27.98 (363.5) | 20.07 |
| 60 | temperature | 0.01 | 720 (0.2) | 1.5 | 358.71 | 60 | 149.0 (149.0) | 721.69 | 0.967 | 758.2 | 723.91 (1.0) | 0.0609 | 27.3 | 26.3 | 27.2 | 29.6 | 1.03 | 1.09 | 1.08 | 149.0 | 32.18 (359.0) | 28.26 |
| 60 | temperature | 0.03 | 720 (0.2) | 1.5 | 310.21 | 35 | 218.0 (218.0) | 678.88 | 0.947 | 681.73 | 678.88 (1.0) | 0.001 | 6.82 | 6.72 | 7.1 | 7.1 | 1.06 | 1.0 | 1.04 | 1670.0 | 32.39 (335.0) | 31.33 |
| 120 | temperature | 0.003 | 720 (0.05) | 1.5 | 363.13 | 120 | 32.5 (32.5) | 758.48 | 0.951 | 778.42 | 760.98 (1.0) | 0.126 | 21.4 | 21.2 | 22.3 | 24.3 | 1.05 | 1.09 | 1.14 | 32.4 | 23.0 (363.5) | 20.67 |
| 120 | temperature | 0.01 | 720 (0.05) | 1.5 | 358.64 | 120 | 148.0 (148.0) | 722.32 | 0.965 | 747.06 | 723.11 (1.0) | 0.0321 | 27.7 | 26.6 | 27.6 | 28.4 | 1.04 | 1.03 | 1.03 | 148.0 | 30.74 (359.0) | 28.12 |
| 120 | temperature | 0.03 | 720 (0.05) | 1.5 | 309.28 | 69 | 206.0 (206.0) | 678.68 | 0.944 | 681.64 | 678.68 (1.0) | 0.000946 | 6.36 | 6.27 | 6.64 | 6.64 | 1.06 | 1.0 | 1.04 | 1700.0 | 32.47 (335.5) | 31.24 |
| 60 | temperature | 0.01 | 360 (0.1) | 1.5 | 361.43 | 60 | 78.4 (78.4) | 739.49 | 0.823 | 776.14 | 743.45 (1.0) | 0.108 | 23.0 | 22.7 | 27.6 | 31.8 | 1.21 | 1.15 | 1.38 | 78.4 | 33.01 (361.5) | 25.37 |
| 60 | temperature | 0.01 | 3600 (1.0) | 1.5 | 301.95 | 34 | 137.0 (137.0) | 676.78 | 0.982 | 679.84 | 676.78 (1.0) | 0.000721 | 4.04 | 4.01 | 4.08 | 4.08 | 1.02 | 1.0 | 1.01 | 1920.0 | 32.39 (335.0) | 30.46 |
| 30 | temperature | 0.01 | 720 (0.8) | 0.0 | 358.66 | 30 | 153.0 (153.0) | 719.94 | 1.0 | 778.7 | 726.51 (1.0) | 0.112 | 26.8 | 26.1 | 26.1 | 33.4 | 1.0 | 1.28 | 1.25 | 153.0 | 35.11 (359.0) | 29.54 |
| 60 | temperature | 0.01 | 720 (0.2) | 0.0 | 358.2 | 60 | 163.0 (163.0) | 719.15 | 1.0 | 755.48 | 721.19 (1.0) | 0.0561 | 27.8 | 27.0 | 27.0 | 29.2 | 1.0 | 1.08 | 1.05 | 163.0 | 32.05 (358.5) | 28.26 |
| 120 | temperature | 0.01 | 720 (0.05) | 0.0 | 358.09 | 120 | 163.0 (163.0) | 719.53 | 1.0 | 744.54 | 720.27 (1.0) | 0.0294 | 28.3 | 27.4 | 27.4 | 28.2 | 1.0 | 1.03 | 0.997 | 163.0 | 30.78 (358.5) | 28.19 |
| 30 | enthalpy | 0.003 | 720 (0.8) | 1.5 | 363.7 | 30 | 30.9 (30.9) | 757.3 | 0.943 | 822.72 | 781.83 (1.0) | 0.375 | 19.4 | 19.2 | 20.4 | 46.1 | 1.06 | 2.26 | 2.38 | 32.4 | 28.32 (364.0) | 17.67 |
| 30 | enthalpy | 0.01 | 720 (0.8) | 1.5 | 358.99 | 30 | 144.0 (144.0) | 721.72 | 0.965 | 781.01 | 728.68 (1.0) | 0.117 | 26.4 | 25.4 | 26.3 | 34.2 | 1.04 | 1.3 | 1.3 | 144.0 | 34.51 (359.5) | 29.54 |
| 30 | enthalpy | 0.03 | 720 (0.8) | 1.5 | 312.14 | 17 | 217.0 (217.0) | 678.95 | 0.947 | 681.79 | 678.95 (1.0) | 0.00101 | 6.81 | 6.71 | 7.08 | 7.08 | 1.06 | 1.0 | 1.04 | 1600.0 | 37.2 (360.5) | 31.51 |
| 60 | enthalpy | 0.003 | 720 (0.2) | 1.5 | 363.46 | 60 | 30.2 (30.2) | 759.25 | 0.935 | 794.63 | 767.6 (1.0) | 0.236 | 20.0 | 19.9 | 21.3 | 28.2 | 1.07 | 1.33 | 1.41 | 30.7 | 29.21 (363.5) | 20.07 |
| 60 | enthalpy | 0.01 | 720 (0.2) | 1.5 | 358.69 | 60 | 150.0 (150.0) | 721.6 | 0.967 | 758.1 | 723.82 (1.0) | 0.0607 | 27.3 | 26.3 | 27.2 | 29.6 | 1.03 | 1.09 | 1.08 | 150.0 | 32.18 (359.0) | 28.26 |
| 60 | enthalpy | 0.03 | 720 (0.2) | 1.5 | 310.26 | 35 | 218.0 (218.0) | 678.88 | 0.947 | 681.73 | 678.89 (1.0) | 0.001 | 6.82 | 6.72 | 7.09 | 7.09 | 1.06 | 1.0 | 1.04 | 1670.0 | 32.39 (335.0) | 31.33 |
| 120 | enthalpy | 0.003 | 720 (0.05) | 1.5 | 363.22 | 120 | 30.9 (30.9) | 759.34 | 0.94 | 779.17 | 761.95 (1.0) | 0.131 | 20.6 | 20.5 | 21.8 | 23.8 | 1.06 | 1.09 | 1.15 | 30.8 | 23.41 (363.5) | 20.67 |
| 120 | enthalpy | 0.01 | 720 (0.05) | 1.5 | 358.3 | 120 | 158.0 (158.0) | 720.73 | 0.971 | 745.59 | 721.48 (1.0) | 0.0303 | 27.9 | 26.9 | 27.7 | 28.5 | 1.03 | 1.03 | 1.02 | 157.0 | 30.69 (358.5) | 28.19 |
| 120 | enthalpy | 0.03 | 720 (0.05) | 1.5 | 309.3 | 69 | 206.0 (206.0) | 678.68 | 0.944 | 681.64 | 678.68 (1.0) | 0.000946 | 6.36 | 6.27 | 6.64 | 6.64 | 1.06 | 1.0 | 1.04 | 1700.0 | 32.47 (335.5) | 31.24 |
| 60 | enthalpy | 0.01 | 360 (0.1) | 1.5 | 361.64 | 60 | 73.6 (73.6) | 742.32 | 0.768 | 778.78 | 746.48 (1.0) | 0.114 | 22.2 | 22.0 | 28.7 | 33.2 | 1.3 | 1.16 | 1.5 | 73.5 | 30.43 (362.0) | 24.21 |
| 60 | enthalpy | 0.01 | 3600 (1.0) | 1.5 | 302.56 | 34 | 136.0 (136.0) | 676.76 | 0.982 | 679.82 | 676.76 (1.0) | 0.000715 | 3.98 | 3.95 | 4.02 | 4.03 | 1.02 | 1.0 | 1.01 | 1900.0 | 32.39 (335.0) | 30.57 |
| 30 | enthalpy | 0.01 | 720 (0.8) | 0.0 | 358.49 | 30 | 158.0 (158.0) | 719.23 | 1.0 | 777.59 | 725.58 (1.0) | 0.109 | 26.9 | 26.2 | 26.2 | 33.3 | 1.0 | 1.27 | 1.24 | 157.0 | 34.48 (359.0) | 29.54 |
| 60 | enthalpy | 0.01 | 720 (0.2) | 0.0 | 358.23 | 60 | 162.0 (162.0) | 719.29 | 1.0 | 755.64 | 721.34 (1.0) | 0.0564 | 27.8 | 27.0 | 27.0 | 29.2 | 1.0 | 1.08 | 1.05 | 162.0 | 32.0 (358.5) | 28.26 |
| 120 | enthalpy | 0.01 | 720 (0.05) | 0.0 | 357.87 | 120 | 169.0 (169.0) | 718.63 | 1.0 | 743.68 | 719.34 (1.0) | 0.0284 | 28.4 | 27.5 | 27.5 | 28.2 | 1.0 | 1.03 | 0.995 | 169.0 | 30.71 (358.5) | 28.3 |

## The old core

Curve against the next resolution: the curve criterion of the old core at these cells against the old core at twice as many, and its time [s].

| case | cells | retcode | parks | wall [s] | curve against the next resolution |
|---|---|---|---|---|---|
| V-BO-01 | 30 | Success | 30 | 0.724 | 18.4 % (365.0) |
| V-BO-01 | 60 | Success | 60 | 1.76 | 6.79 % (365.0) |
| V-BO-01 | 120 | Success | 120 | 7.27 |  |
| V-BO-03 | 30 | Success | 24 | 0.708 | 0.725 % (357.0) |
| V-BO-03 | 60 | Success | 48 | 1.48 | 0.184 % (356.5) |
| V-BO-03 | 120 | Success | 96 | 5.06 |  |
| V-BO-04 | 40 | Success | 35 | 0.0771 |  |
