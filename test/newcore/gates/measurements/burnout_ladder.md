# Gate items 1, 2, 4, 6 and 10: the burnout ladder

Produced by `test/newcore/gates/burnout_ladder.jl` at commit 5f08134; item 1 recomputed from this table's csv at commit c36760a (tree 3307204). Julia 1.13.0, 1 thread, 22 logical CPUs; load average [24.0, 20.7, 18.4] at the start and [18.7, 18.9, 19.8] at the end. The finite-volume core with the equidistribution correction (band 0.1, activation 0.3/0.1, the monitor's mass term `β` per row), the layer interfaces Lagrangian; τ [s] is the layer's equidistribution time, with the prototype's per-cell `τ_p = τ/N²` [s] in parentheses. The old core runs each case at the same resolution with its own keywords (its threshold is a thickness fraction).

## Outcomes and QoIs (items 1, 2 and 4)

Invocations: per invocation its cells (below ε_m, taken by the band; the largest merged mass fraction). Item 1 passes with exactly one invocation per layer and the largest fraction at most 5 ε_m. QoIs: the case's own, in assert mode against its fixed references; the failing ones are named. Inadmissible: results of a deactivation the solve refused.

| case | cells | form | ε_m | τ (τ_p) | correction | retcode, outcome, end [s] | invocations | item 1 | QoIs | failing | peak, t_peak, t_bo (old core) | events | inadmissible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 30 | temperature | 0.003 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 30 (1 below, 29 band; 0.614 %) | pass | 8/8 |  | 32.248, 333.0, 364.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 30 (1 below, 29 band; 3.84 %) | pass | 8/8 |  | 32.248, 333.0, 363.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 30 | temperature | 0.03 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 17 (1 below, 16 band; 9.71 %); 13 (1 below, 12 band; 4.48 %) | fail | 7/8 | one invocation for the PMMA layer | 32.254, 333.0, 363.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.632 %) | pass | 8/8 |  | 32.391, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 4.08 %) | pass | 8/8 |  | 32.391, 335.0, 363.0 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 33 (1 below, 32 band; 9.48 %); 27 (1 below, 26 band; 4.35 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 120 | temperature | 0.003 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 120 (1 below, 119 band; 0.64 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 120 (1 below, 119 band; 4.19 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 67 (0 below, 67 band; 9.94 %); 53 (1 below, 52 band; 4.43 %) | fail | 7/8 | one invocation for the PMMA layer | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.22 %) | pass | 8/8 |  | 32.39, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | β = 0.0 | Success, running, 500.0 | 33 (1 below, 32 band; 9.05 %); 27 (1 below, 26 band; 1.29 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 3100.0 | 30 (0 below, 30 band; 1.0 %) | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 30 | 0 |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 3100.0 | 60 (0 below, 60 band; 1.0 %) | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 60 | 0 |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 3100.0 | 120 (0 below, 120 band; 1.0 %) | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 120 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | none | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | β = 0.0 | Success, running, 480.0 | 39 (0 below, 39 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-06 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.569 %) | pass | 8/8 |  | 32.288, 334.0, 364.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.6 %) | pass | 8/8 |  | 32.288, 334.0, 364.0 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 43 (1 below, 42 band; 9.8 %); 17 (1 below, 16 band; 6.8 %) | fail | 7/8 | one invocation for the PMMA layer | 32.36, 334.0, 363.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 5.82 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 1.0 | Success, running, 500.0 | 60 (1 below, 59 band; 6.07 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-08 | 10 | temperature | 0.003 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 1.48 %) | pass | 6/6 |  | 29.712, 47.2, 69.3 (29.712, 47.2, 69.0) | 10 | 0 |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 4.53 %) | pass | 6/6 |  | 29.712, 47.2, 69.2 (29.712, 47.2, 69.0) | 10 | 0 |
| V-BO-08 | 10 | temperature | 0.03 | 720 (7.2) | β = 0.0 | Terminated, interior, 67.04 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.2,  (29.712, 47.2, 69.0) | 0 | 0 |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | none | Terminated, interior, 68.86 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.2,  (29.712, 47.2, 69.0) | 0 | 0 |
| V-BO-01 | 30 | enthalpy | 0.003 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 30 (1 below, 29 band; 0.613 %) | pass | 8/8 |  | 32.248, 333.0, 364.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 30 (1 below, 29 band; 3.84 %) | pass | 8/8 |  | 32.248, 333.0, 363.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 30 | enthalpy | 0.03 | 720 (0.8) | β = 0.0 | Success, running, 500.0 | 17 (1 below, 16 band; 9.71 %); 13 (1 below, 12 band; 4.48 %) | fail | 7/8 | one invocation for the PMMA layer | 32.254, 333.0, 363.0 (32.254, 333.0, 364.0) | 30 | 0 |
| V-BO-01 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.631 %) | pass | 8/8 |  | 32.391, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 4.08 %) | pass | 8/8 |  | 32.391, 335.0, 363.0 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 33 (1 below, 32 band; 9.48 %); 27 (1 below, 26 band; 4.35 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 120 | enthalpy | 0.003 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 120 (1 below, 119 band; 0.639 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 120 (1 below, 119 band; 4.19 %) | pass | 8/8 |  | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 120 | enthalpy | 0.03 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 67 (1 below, 66 band; 9.94 %); 53 (1 below, 52 band; 4.43 %) | fail | 7/8 | one invocation for the PMMA layer | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 360 (0.1) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.22 %) | pass | 8/8 |  | 32.39, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | β = 0.0 | Success, running, 500.0 | 33 (1 below, 32 band; 9.05 %); 27 (1 below, 26 band; 1.29 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 3100.0 | 30 (0 below, 30 band; 1.0 %) | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 30 | 0 |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 3100.0 | 60 (0 below, 60 band; 1.0 %) | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 60 | 0 |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 3100.0 | 120 (0 below, 120 band; 1.0 %) | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 120 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | none | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | β = 0.0 | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-06 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.568 %) | pass | 8/8 |  | 32.288, 334.0, 364.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.6 %) | pass | 8/8 |  | 32.288, 334.0, 364.0 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 43 (1 below, 42 band; 9.8 %); 17 (1 below, 16 band; 6.8 %) | fail | 7/8 | one invocation for the PMMA layer | 32.362, 334.0, 363.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 5.82 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 1.0 | Success, running, 500.0 | 60 (1 below, 59 band; 6.07 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-08 | 10 | enthalpy | 0.003 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 1.48 %) | pass | 6/6 |  | 29.712, 47.2, 69.3 (29.712, 47.2, 69.0) | 10 | 0 |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 4.53 %) | pass | 6/6 |  | 29.712, 47.2, 69.2 (29.712, 47.2, 69.0) | 10 | 0 |
| V-BO-08 | 10 | enthalpy | 0.03 | 720 (7.2) | β = 0.0 | Terminated, interior, 67.04 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.2,  (29.712, 47.2, 69.0) | 0 | 0 |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | none | Terminated, interior, 68.86 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.2,  (29.712, 47.2, 69.0) | 0 | 0 |

## Cells, residuals and cost (items 6 and 10)

Smallest cell: the thinnest live cell at the end of any step, before its deactivation event [m], and the smallest ratio of a live cell to its initial thickness over ε_m (item 10: at least 0.5, and no `DtLessThanMin`). C52 and Phase 1: the largest scaled constraint residual of the accepted states under C52's scale and under Phase 1's, which omits the sweep of the node corrections; Phase 1 > 1: the accepted states Phase 1's scale would have rejected. Domain: the domain test's rejections (C52's scale, admissibility included). Wall: the warm median of 3 solves without the recorders [s], the old core's at the same resolution in parentheses; same: whether those solves took the instrumented solve's steps and rejections (in the enthalpy form the recorders' residual evaluations move the warm start of the caloric inversion, so the two solves part at round-off).

| case | cells | form | ε_m | τ (τ_p) | correction | smallest cell [m] | ratio / ε_m | item 10 | C52 | Phase 1 | Phase 1 > 1 | domain | steps (rejected) | wall (old core) | same |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 30 | temperature | 0.003 | 720 (0.8) | β = 0.0 | 5.29e-7 | 0.897 | pass | 0.369 | 0.372 | 0 | 24 | 745 (223) | 0.68 (0.674) | true |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | 1.88e-6 | 0.954 | pass | 0.265 | 0.265 | 0 | 21 | 686 (220) | 0.659 (0.674) | true |
| V-BO-01 | 30 | temperature | 0.03 | 720 (0.8) | β = 0.0 | 5.18e-6 | 0.878 | pass | 0.0593 | 0.0593 | 0 | 18 | 725 (262) | 0.704 (0.674) | true |
| V-BO-01 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | 2.62e-7 | 0.887 | pass | 0.993 | 0.994 | 0 | 68 | 1188 (201) | 2.06 (1.58) | true |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 9.57e-7 | 0.973 | pass | 0.609 | 0.609 | 0 | 61 | 1127 (225) | 1.97 (1.58) | true |
| V-BO-01 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | 2.89e-6 | 0.978 | pass | 0.27 | 0.27 | 0 | 45 | 1098 (224) | 2.18 (1.58) | true |
| V-BO-01 | 120 | temperature | 0.003 | 720 (0.05) | β = 0.0 | 1.33e-7 | 0.902 | pass | 0.942 | 0.971 | 0 | 174 | 2323 (570) | 10.3 (3.56) | true |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | 4.74e-7 | 0.964 | pass | 0.967 | 0.976 | 0 | 151 | 2118 (493) | 9.11 (3.56) | true |
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | β = 0.0 | 1.47e-6 | 0.994 | pass | 0.919 | 0.921 | 0 | 128 | 1981 (452) | 7.48 (3.56) | true |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | β = 0.0 | 9.42e-7 | 0.958 | pass | 0.954 | 0.955 | 0 | 71 | 1294 (239) | 2.33 (1.58) | true |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | β = 0.0 | 8.51e-7 | 0.865 | pass | 0.222 | 0.222 | 0 | 2 | 704 (166) | 1.32 (1.58) | true |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | 1.96e-6 | 0.997 | pass | 0.123 | 0.124 | 0 | 21 | 678 (211) | 0.834 (0.456) | true |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 9.8e-7 | 0.997 | pass | 0.802 | 0.802 | 0 | 60 | 1121 (213) | 2.05 (1.29) | true |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | 4.9e-7 | 0.997 | pass | 0.833 | 0.834 | 0 | 166 | 2168 (501) | 9.45 (7.31) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | none | 2.34e-6 | 0.936 | pass | 1.13e-11 | 1.13e-11 | 0 | 0 | 322 (0) | 0.0922 (0.076) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | β = 0.0 | 2.34e-6 | 0.936 | pass | 6.15e-10 | 6.15e-10 | 0 | 0 | 322 (0) | 0.153 (0.076) | true |
| V-BO-06 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | 2.51e-7 | 0.987 | pass | 0.948 | 0.96 | 0 | 64 | 1251 (233) | 3.15 (1.2) | true |
| V-BO-06 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 5.37e-7 | 0.874 | pass | 0.664 | 0.664 | 0 | 50 | 1137 (214) | 2.84 (1.2) | true |
| V-BO-06 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | 1.64e-6 | 0.869 | pass | 0.23 | 0.23 | 0 | 43 | 1116 (205) | 2.05 (1.2) | true |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 1.92e-6 | 1.95 | pass | 0.684 | 0.687 | 0 | 9 | 732 (215) | 2.52 (2.64) | true |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 1.0 | 1.67e-6 | 1.7 | pass | 0.496 | 0.499 | 0 | 8 | 731 (214) | 2.06 (2.64) | true |
| V-BO-08 | 10 | temperature | 0.003 | 720 (7.2) | β = 0.0 | 2.76e-7 | 0.934 | pass | 0.145 | 0.145 | 0 | 0 | 235 (56) | 0.61 (0.594) | true |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | β = 0.0 | 8.92e-7 | 0.907 | pass | 0.0498 | 0.0498 | 0 | 0 | 213 (54) | 0.409 (0.594) | true |
| V-BO-08 | 10 | temperature | 0.03 | 720 (7.2) | β = 0.0 | 2.88e-6 | 0.976 | pass | 0.0411 | 0.0411 | 0 | 0 | 109 (41) | 0.23 (0.594) | true |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | none | 8.67e-7 | 0.882 | pass | 0.0431 | 0.0431 | 0 | 0 | 120 (45) | 0.196 (0.594) | true |
| V-BO-01 | 30 | enthalpy | 0.003 | 720 (0.8) | β = 0.0 | 5.87e-7 | 0.995 | pass | 0.179 | 0.181 | 0 | 23 | 965 (313) | 4.2 (0.674) | false |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | 1.91e-6 | 0.972 | pass | 0.119 | 0.119 | 0 | 22 | 919 (319) | 2.43 (0.674) | false |
| V-BO-01 | 30 | enthalpy | 0.03 | 720 (0.8) | β = 0.0 | 5.7e-6 | 0.965 | pass | 0.0364 | 0.0364 | 0 | 15 | 927 (323) | 2.43 (0.674) | false |
| V-BO-01 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | 2.83e-7 | 0.958 | pass | 0.909 | 0.929 | 0 | 56 | 1634 (569) | 8.73 (1.58) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 9.54e-7 | 0.97 | pass | 0.469 | 0.47 | 0 | 49 | 1580 (613) | 8.1 (1.58) | false |
| V-BO-01 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | 2.71e-6 | 0.919 | pass | 0.137 | 0.137 | 0 | 43 | 1565 (586) | 11.1 (1.58) | false |
| V-BO-01 | 120 | enthalpy | 0.003 | 720 (0.05) | β = 0.0 | 1.43e-7 | 0.97 | pass | 0.984 | 0.986 | 0 | 160 | 2946 (929) | 61.2 (3.56) | false |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | 4.84e-7 | 0.985 | pass | 0.945 | 0.945 | 0 | 143 | 2789 (860) | 68.4 (3.56) | false |
| V-BO-01 | 120 | enthalpy | 0.03 | 720 (0.05) | β = 0.0 | 1.47e-6 | 0.998 | pass | 0.612 | 0.612 | 0 | 111 | 2644 (846) | 56.8 (3.56) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 360 (0.1) | β = 0.0 | 9.13e-7 | 0.929 | pass | 0.895 | 0.895 | 0 | 69 | 1765 (580) | 11.6 (1.58) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | β = 0.0 | 9.14e-7 | 0.93 | pass | 0.16 | 0.16 | 0 | 0 | 1165 (565) | 14.7 (1.58) | false |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | 1.96e-6 | 0.997 | pass | 0.0977 | 0.0977 | 0 | 22 | 940 (329) | 6.67 (0.456) | false |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 9.81e-7 | 0.997 | pass | 0.306 | 0.306 | 0 | 53 | 1610 (576) | 16.6 (1.29) | false |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | 4.9e-7 | 0.997 | pass | 0.837 | 0.839 | 0 | 143 | 2734 (758) | 67.7 (7.31) | false |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | none | 2.5e-6 | 0.998 | pass | 1.65e-8 | 1.65e-8 | 0 | 0 | 322 (0) | 0.186 (0.076) | true |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | β = 0.0 | 2.5e-6 | 0.998 | pass | 1.65e-8 | 1.65e-8 | 0 | 0 | 322 (0) | 0.28 (0.076) | true |
| V-BO-06 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | 2.41e-7 | 0.938 | pass | 0.988 | 0.989 | 0 | 58 | 1709 (648) | 14.4 (1.2) | false |
| V-BO-06 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 5.55e-7 | 0.926 | pass | 0.419 | 0.419 | 0 | 48 | 1639 (671) | 15.0 (1.2) | false |
| V-BO-06 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | 1.64e-6 | 0.876 | pass | 0.104 | 0.104 | 0 | 44 | 1643 (644) | 15.0 (1.2) | false |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 1.92e-6 | 1.95 | pass | 0.572 | 0.574 | 0 | 7 | 1150 (524) | 18.8 (2.64) | false |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 1.0 | 1.67e-6 | 1.7 | pass | 0.537 | 0.54 | 0 | 9 | 1149 (526) | 22.1 (2.64) | false |
| V-BO-08 | 10 | enthalpy | 0.003 | 720 (7.2) | β = 0.0 | 2.95e-7 | 0.998 | pass | 0.128 | 0.128 | 0 | 1 | 306 (72) | 2.8 (0.594) | true |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | β = 0.0 | 8.66e-7 | 0.881 | pass | 0.105 | 0.105 | 0 | 0 | 279 (87) | 2.63 (0.594) | false |
| V-BO-08 | 10 | enthalpy | 0.03 | 720 (7.2) | β = 0.0 | 2.83e-6 | 0.96 | pass | 0.048 | 0.048 | 0 | 0 | 181 (88) | 1.25 (0.594) | false |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | none | 9.79e-7 | 0.996 | pass | 0.105 | 0.105 | 0 | 0 | 188 (81) | 1.05 (0.594) | false |

The script's C52 scale against `Integration.constraint_residual` at every probed state: largest relative difference 0.0.

## R25: the end of the layer (V-BO-01)

At the state that ends the step of the layer's invocation, before the merge, and after it: the cells merged, the residue `m_res` [g/m²] (after the merge, in the receiver), the mass-weighted temperature `T̄` of the merged PMMA [K], the receiver's temperature before the merge and the merged temperature `T_m` [K], and the receiver's weight in it, `(T_m − T̄)/(T_receiver − T̄)`; the PMMA gasification rate [g m⁻² s⁻¹] `R0 = Σ mᵢ k(Tᵢ)` before the merge, `Ra = m_res k(T̄)` lumped and `Rc = m_res k(T_m)` after the merge; the mass gasified from the event to the end of the run [g/m²] (from `Rc` at the event to the first saved time after it, then the saved times); the largest mass loss rate at those saved times and its time, the old core's largest over the ten saved times from the event's.

| cells | form | ε_m | τ (τ_p) | event [s] | merged | m_res (after) | T̄ | T receiver | T_m | weight | R0 | Ra | Rc | Rc/Ra | Rc/R0 | gasified after | largest after (t) | old core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | temperature | 0.003 | 720 (0.8) | 363.55 | 30 | 29.8 (33.6) | 757.9 | 823.37 | 778.46 | 0.314 | 20.1 | 20.0 | 45.0 | 2.25 | 2.24 | 35.6 | 26.2 (364.0) | 17.67 |
| 30 | temperature | 0.01 | 720 (0.8) | 358.2 | 30 | 153.0 (165.0) | 719.94 | 778.7 | 724.15 | 0.0717 | 26.8 | 26.1 | 33.1 | 1.27 | 1.23 | 165.0 | 35.09 (358.5) | 29.54 |
| 30 | temperature | 0.03 | 720 (0.8) | 307.55 | 17 | 229.0 (229.0) | 678.53 | 681.35 | 678.53 | 0.00107 | 7.39 | 7.35 | 7.35 | 1.0 | 0.994 | 1740.0 | 34.82 (360.0) | 31.11 |
| 60 | temperature | 0.003 | 720 (0.2) | 363.25 | 60 | 29.9 (34.2) | 759.09 | 794.53 | 763.63 | 0.128 | 21.0 | 20.9 | 28.0 | 1.34 | 1.33 | 34.4 | 26.5 (363.5) | 20.07 |
| 60 | temperature | 0.01 | 720 (0.2) | 357.92 | 60 | 163.0 (171.0) | 719.15 | 755.48 | 719.98 | 0.0227 | 27.8 | 27.0 | 29.2 | 1.08 | 1.05 | 171.0 | 31.94 (358.5) | 28.26 |
| 60 | temperature | 0.03 | 720 (0.2) | 305.59 | 33 | 204.0 (204.0) | 678.12 | 681.16 | 678.12 | 0.00117 | 6.45 | 6.41 | 6.42 | 1.0 | 0.995 | 1810.0 | 32.39 (335.0) | 30.89 |
| 120 | temperature | 0.003 | 720 (0.05) | 363.04 | 120 | 30.7 (34.5) | 759.13 | 779.0 | 758.84 | -0.0147 | 21.6 | 21.5 | 24.0 | 1.11 | 1.11 | 34.1 | 22.63 (363.5) | 20.67 |
| 120 | temperature | 0.01 | 720 (0.05) | 357.71 | 120 | 163.0 (174.0) | 719.53 | 744.54 | 718.71 | -0.033 | 28.3 | 27.4 | 28.3 | 1.03 | 1.0 | 174.0 | 30.83 (358.0) | 28.3 |
| 120 | temperature | 0.03 | 720 (0.05) | 304.63 | 67 | 204.0 (204.0) | 678.09 | 681.15 | 678.1 | 0.000958 | 6.46 | 6.43 | 6.43 | 1.0 | 0.994 | 1840.0 | 32.47 (335.5) | 30.78 |
| 60 | temperature | 0.01 | 360 (0.1) | 361.33 | 60 | 76.9 (79.7) | 737.06 | 774.1 | 740.17 | 0.0841 | 25.0 | 24.8 | 28.7 | 1.16 | 1.15 | 79.6 | 31.06 (361.5) | 25.37 |
| 60 | temperature | 0.01 | 3600 (1.0) | 298.79 | 33 | 132.0 (133.0) | 676.48 | 679.52 | 676.49 | 0.00493 | 3.88 | 3.86 | 3.89 | 1.01 | 1.0 | 2010.0 | 32.39 (335.0) | 30.12 |
| 30 | enthalpy | 0.003 | 720 (0.8) | 363.55 | 30 | 33.4 (33.6) | 754.85 | 820.55 | 778.12 | 0.354 | 20.3 | 20.2 | 44.5 | 2.2 | 2.19 | 35.4 | 26.31 (364.0) | 17.67 |
| 30 | enthalpy | 0.01 | 720 (0.8) | 358.2 | 30 | 158.0 (165.0) | 719.23 | 777.59 | 724.15 | 0.0844 | 26.9 | 26.2 | 33.1 | 1.26 | 1.23 | 165.0 | 35.08 (358.5) | 29.54 |
| 30 | enthalpy | 0.03 | 720 (0.8) | 307.55 | 17 | 228.0 (229.0) | 678.53 | 681.36 | 678.53 | 0.00121 | 7.35 | 7.31 | 7.35 | 1.01 | 0.999 | 1740.0 | 34.86 (360.0) | 31.11 |
| 60 | enthalpy | 0.003 | 720 (0.2) | 363.25 | 60 | 32.6 (34.2) | 756.63 | 792.35 | 763.27 | 0.186 | 21.0 | 21.0 | 27.6 | 1.32 | 1.31 | 34.4 | 26.52 (363.5) | 20.07 |
| 60 | enthalpy | 0.01 | 720 (0.2) | 357.91 | 60 | 162.0 (171.0) | 719.29 | 755.64 | 719.98 | 0.0189 | 27.8 | 27.0 | 29.2 | 1.08 | 1.05 | 171.0 | 31.94 (358.5) | 28.26 |
| 60 | enthalpy | 0.03 | 720 (0.2) | 305.59 | 33 | 204.0 (204.0) | 678.12 | 681.16 | 678.12 | 0.00145 | 6.44 | 6.4 | 6.42 | 1.0 | 0.997 | 1810.0 | 32.39 (335.0) | 30.89 |
| 120 | enthalpy | 0.003 | 720 (0.05) | 363.04 | 120 | 33.3 (34.5) | 756.79 | 776.95 | 758.47 | 0.0835 | 21.6 | 21.6 | 23.7 | 1.1 | 1.09 | 34.1 | 22.67 (363.5) | 20.67 |
| 120 | enthalpy | 0.01 | 720 (0.05) | 357.71 | 120 | 169.0 (174.0) | 718.63 | 743.68 | 718.71 | 0.00294 | 28.4 | 27.5 | 28.3 | 1.03 | 0.996 | 174.0 | 30.83 (358.0) | 28.3 |
| 120 | enthalpy | 0.03 | 720 (0.05) | 304.63 | 67 | 204.0 (204.0) | 678.09 | 681.15 | 678.1 | 0.00132 | 6.45 | 6.41 | 6.43 | 1.0 | 0.996 | 1840.0 | 32.47 (335.5) | 30.78 |
| 60 | enthalpy | 0.01 | 360 (0.1) | 361.33 | 60 | 75.0 (79.7) | 737.66 | 774.68 | 740.17 | 0.068 | 24.9 | 24.7 | 28.7 | 1.16 | 1.15 | 79.6 | 31.05 (361.5) | 25.37 |
| 60 | enthalpy | 0.01 | 3600 (1.0) | 298.78 | 33 | 131.0 (133.0) | 676.47 | 679.5 | 676.49 | 0.00912 | 3.85 | 3.83 | 3.89 | 1.02 | 1.01 | 2010.0 | 32.39 (335.0) | 30.12 |

## V-BO-08: the substrate before the first activation

The temperatures of the cells that start without PMMA at the saved times up to the end of the last accepted step before the first whose state has an active node (the step that crosses the activation moves the saved times its interpolant spans): the corrected run at ε_m = 1e-2 against the run without the correction, the largest difference and that over the largest temperature (§6.2: at most 1e-10 on the same accepted steps); the accepted steps of the two solves up to the last saved time compared, and the time of the first step at which they part (a changed sequence is reported with the difference, as item 5 does). The run without the correction is this comparison's reference only: its layer ends `:interior`, and it is not a run of items 1 and 2.

| form | first active step [s] | compared to [s] | saved times | largest difference [K] | relative | steps (without) | parted at [s] | verdict |
|---|---|---|---|---|---|---|---|---|
| temperature | 56.05 | 54.884 | 549 | 0.0 | 0.0 | 95 (95) |  | pass |
| enthalpy | 56.215 | 55.157 | 552 | 0.0 | 0.0 | 163 (163) |  | pass |

## The old core

| case | cells | retcode | parks | wall [s] |
|---|---|---|---|---|
| V-BO-01 | 30 | Success | 30 | 0.674 |
| V-BO-01 | 60 | Success | 60 | 1.58 |
| V-BO-01 | 120 | Success | 120 | 3.56 |
| V-BO-03 | 30 | Success | 26 | 0.456 |
| V-BO-03 | 60 | Success | 52 | 1.29 |
| V-BO-03 | 120 | Success | 103 | 7.31 |
| V-BO-04 | 40 | Success | 35 | 0.076 |
| V-BO-06 | 60 | Success | 60 | 1.2 |
| V-BO-07 | 60 | Success | 60 | 2.64 |
| V-BO-08 | 10 | Success | 10 | 0.594 |
