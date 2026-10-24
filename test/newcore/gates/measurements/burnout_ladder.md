# Gate items 1, 2, 4, 6 and 10: the burnout ladder

Produced by `test/newcore/gates/burnout_ladder.jl` at commit f54a6fb (tree b544f02). Julia 1.13.0, 1 thread, 22 logical CPUs; load average [20.0, 18.3, 18.3] at the start and [19.8, 21.6, 20.5] at the end. The finite-volume core with the equidistribution correction (band 0.1, activation 0.3/0.1, the monitor's mass term `β` per row), the layer interfaces Lagrangian; τ [s] is the layer's equidistribution time, with the prototype's per-cell `τ_p = τ/N²` [s] in parentheses. The old core runs each case at the same resolution with its own keywords (its threshold is a thickness fraction).

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
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | β = 0.0 | Success, running, 500.0 | 67 (1 below, 66 band; 9.94 %); 53 (1 below, 52 band; 4.43 %) | fail | 7/8 | one invocation for the PMMA layer | 32.469, 335.5, 363.5 (32.469, 335.5, 363.5) | 120 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.22 %) | pass | 8/8 |  | 32.39, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | β = 0.0 | Success, running, 500.0 | 33 (1 below, 32 band; 9.05 %); 27 (1 below, 26 band; 1.29 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 3100.0 | 30 (0 below, 30 band; 1.0 %) | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 30 | 0 |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 3100.0 | 60 (0 below, 60 band; 1.0 %) | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 60 | 0 |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 3100.0 | 120 (0 below, 120 band; 1.0 %) | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 120 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | none | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | β = 0.0 | Success, running, 480.0 | 39 (0 below, 39 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-06 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.569 %) | pass | 8/8 |  | 32.288, 334.0, 364.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.6 %) | pass | 8/8 |  | 32.288, 334.0, 364.0 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 43 (1 below, 42 band; 9.8 %); 17 (1 below, 16 band; 6.8 %) | fail | 7/8 | one invocation for the PMMA layer | 32.359, 334.0, 363.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 5.82 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 1.0 | Success, running, 500.0 | 60 (1 below, 59 band; 6.06 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-08 | 10 | temperature | 0.003 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 1.48 %) | pass | 6/6 |  | 29.712, 47.0, 69.5 (29.712, 47.0, 69.0) | 10 | 0 |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 4.53 %) | pass | 6/6 |  | 29.712, 47.0, 69.5 (29.712, 47.0, 69.0) | 10 | 0 |
| V-BO-08 | 10 | temperature | 0.03 | 720 (7.2) | β = 0.0 | Terminated, interior, 67.04 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.0,  (29.712, 47.0, 69.0) | 0 | 0 |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | none | Terminated, interior, 68.86 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.0,  (29.712, 47.0, 69.0) | 0 | 0 |
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
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | β = 0.0 | Success, running, 500.0 | 33 (0 below, 33 band; 9.05 %); 27 (1 below, 26 band; 1.29 %) | fail | 7/8 | one invocation for the PMMA layer | 32.392, 335.0, 363.5 (32.392, 335.0, 363.5) | 60 | 0 |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | Success, running, 3100.0 | 30 (0 below, 30 band; 1.0 %) | pass | 4/4 |  | 30.394, 316.5,  (30.4, 316.5, ) | 30 | 0 |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 3100.0 | 60 (0 below, 60 band; 1.0 %) | pass | 4/4 |  | 30.385, 317.0,  (30.387, 317.0, ) | 60 | 0 |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | Success, running, 3100.0 | 120 (0 below, 120 band; 1.0 %) | pass | 4/4 |  | 30.383, 317.0,  (30.383, 317.0, ) | 120 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | none | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | β = 0.0 | Success, running, 480.0 | 39 (1 below, 38 band; 1.0 %) | pass | 5/5 |  |  | 39 | 0 |
| V-BO-06 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 0.568 %) | pass | 8/8 |  | 32.288, 334.0, 364.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 1.6 %) | pass | 8/8 |  | 32.288, 334.0, 364.0 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-06 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 43 (1 below, 42 band; 9.8 %); 17 (1 below, 16 band; 6.8 %) | fail | 7/8 | one invocation for the PMMA layer | 32.359, 334.0, 363.5 (32.292, 334.0, 364.0) | 60 | 0 |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | Success, running, 500.0 | 60 (1 below, 59 band; 5.82 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 1.0 | Success, running, 500.0 | 60 (1 below, 59 band; 6.06 %) | fail | 8/8 |  | 30.542, 287.5, 324.0 (30.542, 287.5, 325.0) | 60 | 0 |
| V-BO-08 | 10 | enthalpy | 0.003 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 1.48 %) | pass | 6/6 |  | 29.712, 47.0, 69.5 (29.712, 47.0, 69.0) | 10 | 0 |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | β = 0.0 | Success, running, 150.0 | 10 (1 below, 9 band; 4.53 %) | pass | 6/6 |  | 29.712, 47.0, 69.5 (29.712, 47.0, 69.0) | 10 | 0 |
| V-BO-08 | 10 | enthalpy | 0.03 | 720 (7.2) | β = 0.0 | Terminated, interior, 67.04 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.0,  (29.712, 47.0, 69.0) | 0 | 0 |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | none | Terminated, interior, 68.86 | 0 | fail | 3/6 | complete burnout: live thickness at t_end is the substrate; mass loss of the gasified PMMA layer; one invocation for the PMMA layer | 29.712, 47.0,  (29.712, 47.0, 69.0) | 0 | 0 |

## Cells, residuals and cost (items 6 and 10)

Smallest cell: the thinnest live cell at the end of any step, before its deactivation event [m], and the smallest ratio of a live cell to its initial thickness over ε_m (item 10: at least 0.5, and no `DtLessThanMin`). On a step whose event is located, the step's end is read before the event is located, a state after the event's time that the trajectory does not pass through: the thinnest cell there is thinner than at the event (item 10's reading is conservative), and the residuals are those of a proposal the domain test accepted rather than of a state on the trajectory. C52 and Phase 1: the largest scaled constraint residual of the accepted states under C52's scale and under Phase 1's, which omits the sweep of the node corrections; Phase 1 > 1: the accepted states Phase 1's scale would have rejected. Domain: the domain test's rejections (C52's scale, admissibility included). Wall: the warm median of 3 solves without the recorders [s], the old core's at the same resolution in parentheses; same: whether those solves took the instrumented solve's steps and rejections (in the enthalpy form the recorders' residual evaluations move the warm start of the caloric inversion, so the two solves part at round-off).

| case | cells | form | ε_m | τ (τ_p) | correction | smallest cell [m] | ratio / ε_m | item 10 | C52 | Phase 1 | Phase 1 > 1 | domain | steps (rejected) | wall (old core) | same |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 30 | temperature | 0.003 | 720 (0.8) | β = 0.0 | 5.71e-7 | 0.968 | pass | 0.381 | 0.381 | 0 | 24 | 744 (227) | 1.41 (0.492) | true |
| V-BO-01 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | 1.96e-6 | 0.994 | pass | 0.184 | 0.184 | 0 | 19 | 685 (230) | 1.28 (0.492) | true |
| V-BO-01 | 30 | temperature | 0.03 | 720 (0.8) | β = 0.0 | 5.64e-6 | 0.956 | pass | 0.0624 | 0.0624 | 0 | 21 | 733 (264) | 1.37 (0.492) | true |
| V-BO-01 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | 2.76e-7 | 0.934 | pass | 0.987 | 0.988 | 0 | 65 | 1177 (201) | 4.24 (1.18) | true |
| V-BO-01 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 9.37e-7 | 0.953 | pass | 0.559 | 0.563 | 0 | 59 | 1125 (225) | 4.23 (1.18) | true |
| V-BO-01 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | 2.9e-6 | 0.982 | pass | 0.27 | 0.271 | 0 | 49 | 1110 (229) | 4.0 (1.18) | true |
| V-BO-01 | 120 | temperature | 0.003 | 720 (0.05) | β = 0.0 | 1.47e-7 | 1.0 | pass | 0.973 | 1.03 | 1 | 178 | 2339 (577) | 20.1 (9.62) | true |
| V-BO-01 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | 4.76e-7 | 0.968 | pass | 0.974 | 0.982 | 0 | 142 | 2083 (483) | 9.26 (9.62) | true |
| V-BO-01 | 120 | temperature | 0.03 | 720 (0.05) | β = 0.0 | 1.32e-6 | 0.893 | pass | 0.97 | 0.97 | 0 | 127 | 1976 (448) | 18.8 (9.62) | true |
| V-BO-01 | 60 | temperature | 0.01 | 360 (0.1) | β = 0.0 | 8.55e-7 | 0.87 | pass | 0.878 | 0.878 | 0 | 74 | 1307 (241) | 6.83 (1.18) | true |
| V-BO-01 | 60 | temperature | 0.01 | 3600 (1.0) | β = 0.0 | 9.58e-7 | 0.974 | pass | 0.191 | 0.191 | 0 | 2 | 713 (173) | 4.13 (1.18) | true |
| V-BO-03 | 30 | temperature | 0.01 | 720 (0.8) | β = 0.0 | 1.96e-6 | 0.997 | pass | 0.21 | 0.21 | 0 | 21 | 674 (201) | 2.34 (1.26) | true |
| V-BO-03 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 9.8e-7 | 0.997 | pass | 0.539 | 0.539 | 0 | 64 | 1131 (210) | 4.97 (2.7) | true |
| V-BO-03 | 120 | temperature | 0.01 | 720 (0.05) | β = 0.0 | 4.9e-7 | 0.997 | pass | 0.998 | 0.999 | 0 | 158 | 2141 (492) | 19.6 (5.83) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | none | 2.34e-6 | 0.936 | pass | 1.13e-11 | 1.13e-11 | 0 | 0 | 322 (0) | 0.0796 (0.0522) | true |
| V-BO-04 | 40 | temperature | 0.01 | 720 (0.45) | β = 0.0 | 2.34e-6 | 0.936 | pass | 5.38e-10 | 5.38e-10 | 0 | 0 | 322 (0) | 0.137 (0.0522) | true |
| V-BO-06 | 60 | temperature | 0.003 | 720 (0.2) | β = 0.0 | 2.37e-7 | 0.918 | pass | 0.956 | 0.97 | 0 | 62 | 1241 (229) | 2.19 (1.62) | true |
| V-BO-06 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 5.58e-7 | 0.935 | pass | 0.873 | 0.874 | 0 | 49 | 1132 (214) | 3.31 (1.62) | true |
| V-BO-06 | 60 | temperature | 0.03 | 720 (0.2) | β = 0.0 | 1.64e-6 | 0.954 | pass | 0.196 | 0.196 | 0 | 46 | 1129 (206) | 3.53 (1.62) | true |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 0.0 | 1.89e-6 | 1.92 | pass | 0.668 | 0.671 | 0 | 5 | 724 (218) | 3.68 (1.66) | true |
| V-BO-07 | 60 | temperature | 0.01 | 720 (0.2) | β = 1.0 | 1.66e-6 | 1.68 | pass | 0.608 | 0.611 | 0 | 9 | 734 (213) | 3.43 (1.66) | true |
| V-BO-08 | 10 | temperature | 0.003 | 720 (7.2) | β = 0.0 | 2.76e-7 | 0.934 | pass | 0.145 | 0.145 | 0 | 0 | 235 (56) | 0.365 (0.548) | true |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | β = 0.0 | 8.92e-7 | 0.907 | pass | 0.0498 | 0.0498 | 0 | 0 | 213 (54) | 0.512 (0.548) | true |
| V-BO-08 | 10 | temperature | 0.03 | 720 (7.2) | β = 0.0 | 2.88e-6 | 0.976 | pass | 0.0411 | 0.0411 | 0 | 0 | 109 (41) | 0.293 (0.548) | true |
| V-BO-08 | 10 | temperature | 0.01 | 720 (7.2) | none | 8.67e-7 | 0.882 | pass | 0.0431 | 0.0431 | 0 | 0 | 120 (45) | 0.287 (0.548) | true |
| V-BO-01 | 30 | enthalpy | 0.003 | 720 (0.8) | β = 0.0 | 5.83e-7 | 0.988 | pass | 0.35 | 0.35 | 0 | 21 | 959 (328) | 4.93 (0.492) | false |
| V-BO-01 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | 1.94e-6 | 0.985 | pass | 0.0947 | 0.0947 | 0 | 20 | 918 (345) | 4.83 (0.492) | false |
| V-BO-01 | 30 | enthalpy | 0.03 | 720 (0.8) | β = 0.0 | 5.26e-6 | 0.892 | pass | 0.0367 | 0.0367 | 0 | 12 | 925 (352) | 4.64 (0.492) | false |
| V-BO-01 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | 2.86e-7 | 0.971 | pass | 0.951 | 0.952 | 0 | 58 | 1658 (617) | 17.2 (1.18) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 9.45e-7 | 0.961 | pass | 0.434 | 0.437 | 0 | 51 | 1583 (581) | 14.8 (1.18) | false |
| V-BO-01 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | 2.56e-6 | 0.867 | pass | 0.122 | 0.122 | 0 | 41 | 1575 (627) | 14.7 (1.18) | false |
| V-BO-01 | 120 | enthalpy | 0.003 | 720 (0.05) | β = 0.0 | 1.41e-7 | 0.957 | pass | 0.964 | 0.997 | 0 | 161 | 2912 (877) | 53.7 (9.62) | false |
| V-BO-01 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | 4.79e-7 | 0.974 | pass | 0.886 | 0.887 | 0 | 135 | 2710 (807) | 45.0 (9.62) | false |
| V-BO-01 | 120 | enthalpy | 0.03 | 720 (0.05) | β = 0.0 | 1.45e-6 | 0.981 | pass | 0.572 | 0.573 | 0 | 109 | 2634 (830) | 32.7 (9.62) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 360 (0.1) | β = 0.0 | 9.12e-7 | 0.928 | pass | 0.988 | 0.989 | 0 | 67 | 1769 (576) | 11.5 (1.18) | false |
| V-BO-01 | 60 | enthalpy | 0.01 | 3600 (1.0) | β = 0.0 | 9.1e-7 | 0.925 | pass | 0.136 | 0.136 | 0 | 0 | 1172 (547) | 16.9 (1.18) | false |
| V-BO-03 | 30 | enthalpy | 0.01 | 720 (0.8) | β = 0.0 | 1.96e-6 | 0.997 | pass | 0.0813 | 0.0814 | 0 | 22 | 927 (303) | 10.2 (1.26) | false |
| V-BO-03 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 9.82e-7 | 0.999 | pass | 0.331 | 0.331 | 0 | 54 | 1611 (580) | 13.0 (2.7) | false |
| V-BO-03 | 120 | enthalpy | 0.01 | 720 (0.05) | β = 0.0 | 4.9e-7 | 0.997 | pass | 0.993 | 0.996 | 0 | 138 | 2810 (910) | 48.0 (5.83) | false |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | none | 2.5e-6 | 0.998 | pass | 1.65e-8 | 1.65e-8 | 0 | 0 | 322 (0) | 0.228 (0.0522) | true |
| V-BO-04 | 40 | enthalpy | 0.01 | 720 (0.45) | β = 0.0 | 2.5e-6 | 0.998 | pass | 1.65e-8 | 1.65e-8 | 0 | 0 | 322 (0) | 0.316 (0.0522) | true |
| V-BO-06 | 60 | enthalpy | 0.003 | 720 (0.2) | β = 0.0 | 2.35e-7 | 0.906 | pass | 0.962 | 0.963 | 0 | 58 | 1715 (654) | 8.95 (1.62) | false |
| V-BO-06 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 5.67e-7 | 0.962 | pass | 0.374 | 0.375 | 0 | 51 | 1620 (606) | 8.12 (1.62) | false |
| V-BO-06 | 60 | enthalpy | 0.03 | 720 (0.2) | β = 0.0 | 1.64e-6 | 0.858 | pass | 0.151 | 0.151 | 0 | 42 | 1624 (634) | 11.5 (1.62) | false |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 0.0 | 1.9e-6 | 1.93 | pass | 0.384 | 0.386 | 0 | 7 | 1179 (596) | 14.9 (1.66) | false |
| V-BO-07 | 60 | enthalpy | 0.01 | 720 (0.2) | β = 1.0 | 1.66e-6 | 1.69 | pass | 0.443 | 0.445 | 0 | 6 | 1174 (592) | 15.0 (1.66) | false |
| V-BO-08 | 10 | enthalpy | 0.003 | 720 (7.2) | β = 0.0 | 2.89e-7 | 0.98 | pass | 0.128 | 0.128 | 0 | 1 | 308 (88) | 1.77 (0.548) | false |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | β = 0.0 | 9.49e-7 | 0.965 | pass | 0.0687 | 0.0687 | 0 | 0 | 273 (72) | 1.63 (0.548) | false |
| V-BO-08 | 10 | enthalpy | 0.03 | 720 (7.2) | β = 0.0 | 2.77e-6 | 0.939 | pass | 0.0379 | 0.0379 | 0 | 0 | 170 (65) | 1.53 (0.548) | false |
| V-BO-08 | 10 | enthalpy | 0.01 | 720 (7.2) | none | 9.07e-7 | 0.922 | pass | 0.0913 | 0.0913 | 0 | 0 | 183 (67) | 1.21 (0.548) | false |

The script's C52 scale against `Integration.constraint_residual` at every probed state: largest relative difference 0.0.

## R25: the end of the layer (V-BO-01)

At the state that ends the step of the layer's invocation, before the merge, and after it: the cells merged, the residue `m_res` [g/m²] (after the merge, in the receiver), the mass-weighted temperature `T̄` of the merged PMMA [K], the receiver's temperature before the merge and the merged temperature `T_m` [K], and the receiver's weight in it, `(T_m − T̄)/(T_receiver − T̄)`; the PMMA gasification rate [g m⁻² s⁻¹] `R0 = Σ mᵢ k(Tᵢ)` before the merge, `Ra = m_res k(T̄)` lumped and `Rc = m_res k(T_m)` after the merge; the mass gasified from the event to the end of the run [g/m²] (from `Rc` at the event to the first saved time after it, then the saved times); the largest mass loss rate at those saved times and its time, the old core's largest over the ten saved times from the event's.

| cells | form | ε_m | τ (τ_p) | event [s] | merged | m_res (after) | T̄ | T receiver | T_m | weight | R0 | Ra | Rc | Rc/Ra | Rc/R0 | gasified after | largest after (t) | old core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | temperature | 0.003 | 720 (0.8) | 363.55 | 30 | 32.4 (33.6) | 755.96 | 821.49 | 778.46 | 0.343 | 20.4 | 20.4 | 45.0 | 2.21 | 2.2 | 35.6 | 26.2 (364.0) | 17.67 |
| 30 | temperature | 0.01 | 720 (0.8) | 358.2 | 30 | 164.0 (165.0) | 718.35 | 776.21 | 724.15 | 0.1 | 27.1 | 26.3 | 33.1 | 1.26 | 1.22 | 165.0 | 35.09 (358.5) | 29.54 |
| 30 | temperature | 0.03 | 720 (0.8) | 307.55 | 17 | 229.0 (229.0) | 678.53 | 681.35 | 678.53 | 0.00117 | 7.36 | 7.32 | 7.35 | 1.0 | 0.998 | 1740.0 | 34.82 (360.0) | 31.11 |
| 60 | temperature | 0.003 | 720 (0.2) | 363.25 | 60 | 31.7 (34.2) | 757.75 | 793.31 | 763.63 | 0.166 | 21.3 | 21.2 | 28.0 | 1.32 | 1.32 | 34.4 | 26.5 (363.5) | 20.07 |
| 60 | temperature | 0.01 | 720 (0.2) | 357.92 | 60 | 157.0 (171.0) | 720.0 | 756.46 | 719.98 | -0.00071 | 27.7 | 26.9 | 29.2 | 1.09 | 1.05 | 171.0 | 31.94 (358.5) | 28.26 |
| 60 | temperature | 0.03 | 720 (0.2) | 305.59 | 33 | 204.0 (204.0) | 678.12 | 681.16 | 678.12 | 0.00134 | 6.45 | 6.41 | 6.42 | 1.0 | 0.996 | 1810.0 | 32.39 (335.0) | 30.89 |
| 120 | temperature | 0.003 | 720 (0.05) | 363.04 | 120 | 34.5 (34.5) | 756.42 | 776.63 | 758.84 | 0.119 | 22.1 | 22.1 | 24.0 | 1.09 | 1.08 | 34.1 | 22.63 (363.5) | 20.67 |
| 120 | temperature | 0.01 | 720 (0.05) | 357.71 | 120 | 164.0 (174.0) | 719.36 | 744.37 | 718.71 | -0.026 | 28.3 | 27.4 | 28.3 | 1.03 | 1.0 | 174.0 | 30.83 (358.0) | 28.3 |
| 120 | temperature | 0.03 | 720 (0.05) | 304.63 | 67 | 204.0 (204.0) | 678.09 | 681.15 | 678.1 | 0.0012 | 6.45 | 6.42 | 6.43 | 1.0 | 0.996 | 1840.0 | 32.47 (335.5) | 30.78 |
| 60 | temperature | 0.01 | 360 (0.1) | 361.33 | 60 | 71.0 (79.7) | 738.94 | 775.9 | 740.17 | 0.0334 | 24.7 | 24.5 | 28.7 | 1.17 | 1.16 | 79.6 | 31.06 (361.5) | 25.37 |
| 60 | temperature | 0.01 | 3600 (1.0) | 298.79 | 33 | 131.0 (133.0) | 676.47 | 679.51 | 676.49 | 0.00652 | 3.87 | 3.85 | 3.89 | 1.01 | 1.01 | 2010.0 | 32.39 (335.0) | 30.12 |
| 30 | enthalpy | 0.003 | 720 (0.8) | 363.55 | 30 | 33.2 (33.6) | 754.99 | 820.7 | 778.12 | 0.352 | 20.2 | 20.2 | 44.5 | 2.2 | 2.2 | 35.4 | 26.31 (364.0) | 17.67 |
| 30 | enthalpy | 0.01 | 720 (0.8) | 358.2 | 30 | 161.0 (165.0) | 718.72 | 776.8 | 724.15 | 0.0936 | 27.0 | 26.2 | 33.1 | 1.26 | 1.22 | 165.0 | 35.08 (358.5) | 29.54 |
| 30 | enthalpy | 0.03 | 720 (0.8) | 307.55 | 17 | 229.0 (229.0) | 678.53 | 681.35 | 678.53 | 0.00112 | 7.38 | 7.33 | 7.35 | 1.0 | 0.996 | 1740.0 | 34.85 (360.0) | 31.11 |
| 60 | enthalpy | 0.003 | 720 (0.2) | 363.25 | 60 | 33.1 (34.2) | 756.31 | 792.05 | 763.27 | 0.195 | 21.1 | 21.1 | 27.6 | 1.31 | 1.31 | 34.4 | 26.52 (363.5) | 20.07 |
| 60 | enthalpy | 0.01 | 720 (0.2) | 357.91 | 60 | 160.0 (171.0) | 719.66 | 756.07 | 719.98 | 0.00867 | 27.8 | 27.0 | 29.2 | 1.08 | 1.05 | 171.0 | 31.94 (358.5) | 28.26 |
| 60 | enthalpy | 0.03 | 720 (0.2) | 305.59 | 33 | 204.0 (204.0) | 678.12 | 681.16 | 678.12 | 0.000999 | 6.46 | 6.42 | 6.42 | 1.0 | 0.994 | 1810.0 | 32.39 (335.0) | 30.89 |
| 120 | enthalpy | 0.003 | 720 (0.05) | 363.04 | 120 | 32.8 (34.5) | 757.14 | 777.25 | 758.47 | 0.0663 | 21.5 | 21.5 | 23.7 | 1.1 | 1.1 | 34.1 | 22.67 (363.5) | 20.67 |
| 120 | enthalpy | 0.01 | 720 (0.05) | 357.71 | 120 | 166.0 (174.0) | 719.12 | 744.15 | 718.71 | -0.0166 | 28.3 | 27.4 | 28.3 | 1.03 | 0.999 | 174.0 | 30.83 (358.0) | 28.3 |
| 120 | enthalpy | 0.03 | 720 (0.05) | 304.63 | 67 | 204.0 (204.0) | 678.09 | 681.15 | 678.1 | 0.00148 | 6.44 | 6.41 | 6.43 | 1.0 | 0.997 | 1840.0 | 32.47 (335.5) | 30.78 |
| 60 | enthalpy | 0.01 | 360 (0.1) | 361.33 | 60 | 75.0 (79.7) | 737.67 | 774.69 | 740.17 | 0.0676 | 24.9 | 24.7 | 28.7 | 1.16 | 1.15 | 79.6 | 31.05 (361.5) | 25.37 |
| 60 | enthalpy | 0.01 | 3600 (1.0) | 298.78 | 33 | 131.0 (133.0) | 676.47 | 679.5 | 676.49 | 0.00912 | 3.85 | 3.83 | 3.89 | 1.02 | 1.01 | 2010.0 | 32.39 (335.0) | 30.12 |

## V-BO-08: the substrate before the first activation

The temperatures of the cells that start without PMMA at the saved times up to the end of the last accepted step before the first whose state has an active node (the step that crosses the activation moves the saved times its interpolant spans): the corrected run at ε_m = 1e-2 against the run without the correction, the largest difference and that over the largest temperature (§6.2: at most 1e-10 on the same accepted steps); the accepted steps of the two solves up to the last saved time compared, and the time of the first step at which they part (a changed sequence is reported with the difference, as item 5 does). The run without the correction is this comparison's reference only: its layer ends `:interior`, and it is not a run of items 1 and 2.

| form | first active step [s] | compared to [s] | saved times | largest difference [K] | relative | steps (without) | parted at [s] | verdict |
|---|---|---|---|---|---|---|---|---|
| temperature | 56.05 | 54.884 | 110 | 0.0 | 0.0 | 95 (95) |  | pass |
| enthalpy | 55.64 | 54.525 | 110 | 0.0 | 0.0 | 157 (157) |  | pass |

## The old core

| case | cells | retcode | parks | wall [s] |
|---|---|---|---|---|
| V-BO-01 | 30 | Success | 30 | 0.492 |
| V-BO-01 | 60 | Success | 60 | 1.18 |
| V-BO-01 | 120 | Success | 120 | 9.62 |
| V-BO-03 | 30 | Success | 26 | 1.26 |
| V-BO-03 | 60 | Success | 52 | 2.7 |
| V-BO-03 | 120 | Success | 103 | 5.83 |
| V-BO-04 | 40 | Success | 35 | 0.0522 |
| V-BO-06 | 60 | Success | 60 | 1.62 |
| V-BO-07 | 60 | Success | 60 | 1.66 |
| V-BO-08 | 10 | Success | 10 | 0.548 |
