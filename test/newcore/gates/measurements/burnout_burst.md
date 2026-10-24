# Gate item 2: the end-of-layer burst by refinement

Produced by `test/newcore/gates/burnout_burst.jl` at commit f54a6fb (tree b544f02). Julia 1.13.0, 1 thread; load average [20.0, 18.3, 18.3] at the start and [19.3, 19.3, 19.2] at the end. The finite-volume core through the case's solve and the old core with the case's keywords, the gas-generation rate at every accepted step and post-event state [g m⁻² s⁻¹]; each deck's window from the old core's 240-cell peak time to its end (V-BO-01 [336.0, 500.0] s, V-BO-03 [317.0, 3100.0] s, V-BO-06 [336.0, 500.0] s, V-BO-07 [288.0, 500.0] s; V-BO-03 is run to 3100.0 s). (a) The largest per-step excess of N cells over 2N after the window's start, falling at least 1.5× per doubling and no larger than the old core's; (b) the L1 of successive grids [g/m²], falling at least 1.5× per doubling, and the L1 against the family's 240-cell run within 1.5× of the old core's; (c) the L1 against the old core's 240-cell run falling with refinement. Gated: V-BO-01 at ε_m of 0.01 and 0.003 in both forms; the rest is reported.

## The criteria

Each measure lists 30→60, 60→120 and 120→240 (the L1 against the old core's 240-cell run: 30, 60, 120 and 240), then its ratios per doubling.

| family | (a) excess of N over 2N; ratios | (a) | (b) successive L1; ratios | (b) self-L1 over the old core's | (b) | (c) L1 against the old core's 240 | (c) | gated |
|---|---|---|---|---|---|---|---|---|
| V-BO-01, old core | 31.9, 7.73, 3.18; 4.13, 2.43 | | 25.0, 9.34, 4.4; 2.68, 2.12 | | | | | |
| V-BO-03, old core | 0.0115, 0.00395, 0.000778; 2.92, 5.08 | | 9.75, 2.53, 0.661; 3.86, 3.82 | | | | | |
| V-BO-06, old core | 114.0, 25.5, 7.29; 4.47, 3.5 | | 51.7, 19.9, 8.44; 2.6, 2.36 | | | | | |
| V-BO-07, old core | 12.0, 3.67, 0.896; 3.28, 4.1 | | 21.0, 9.05, 5.16; 2.32, 1.76 | | | | | |
| V-BO-01, temperature, ε_m = 0.01 | 3.16, 1.36, 0.642; 2.33, 2.12 | pass | 22.0, 10.2, 5.0; 2.15, 2.05 | 1.03, 1.15, 1.14 | pass | 39.7, 19.4, 10.0, 5.92 | pass | yes |
| V-BO-01, temperature, ε_m = 0.003 | 19.1, 5.14, 2.19; 3.72, 2.35 | pass | 22.7, 11.7, 5.85; 1.94, 1.99 | 1.02, 1.27, 1.33 | pass | 36.7, 18.4, 8.21, 3.21 | pass | yes |
| V-BO-01, temperature, ε_m = 0.03 | 4.59, 1.8, 0.908; 2.55, 1.98 | pass | 21.8, 9.95, 5.21; 2.19, 1.91 | 1.01, 1.14, 1.18 | pass | 38.5, 18.4, 8.9, 4.71 | pass | reported |
| V-BO-01, enthalpy, ε_m = 0.01 | 3.17, 1.35, 0.636; 2.34, 2.13 | pass | 22.0, 10.2, 5.01; 2.15, 2.04 | 1.03, 1.15, 1.14 | pass | 39.8, 19.4, 10.0, 5.92 | pass | yes |
| V-BO-01, enthalpy, ε_m = 0.003 | 18.7, 4.98, 2.1; 3.75, 2.37 | pass | 22.7, 11.6, 5.81; 1.95, 2.0 | 1.02, 1.26, 1.32 | pass | 36.8, 18.4, 8.33, 3.23 | pass | yes |
| V-BO-01, enthalpy, ε_m = 0.03 | 4.6, 1.8, 0.901; 2.55, 2.0 | pass | 21.8, 9.93, 5.18; 2.19, 1.92 | 1.01, 1.14, 1.18 | pass | 38.5, 18.5, 8.91, 4.69 | pass | reported |
| V-BO-03, temperature, ε_m = 0.01 | 0.00804, 0.00189, 0.000504; 4.26, 3.74 | pass | 8.74, 2.21, 0.768; 3.96, 2.87 | 0.887, 0.893, 1.16 | pass | 11.1, 2.6, 0.727, 0.49 | pass | reported |
| V-BO-03, enthalpy, ε_m = 0.01 | 0.00602, 0.00184, 0.000545; 3.27, 3.37 | pass | 8.67, 2.27, 0.648; 3.83, 3.5 | 0.894, 0.907, 0.981 | pass | 11.4, 2.75, 0.591, 0.233 | pass | reported |
| V-BO-06, temperature, ε_m = 0.01 | 25.6, 5.23, 2.08; 4.89, 2.52 | pass | 45.1, 20.9, 9.91; 2.15, 2.11 | 0.974, 1.12, 1.17 | pass | 64.1, 29.9, 14.5, 7.85 | pass | reported |
| V-BO-06, enthalpy, ε_m = 0.01 | 25.6, 5.21, 2.08; 4.91, 2.51 | pass | 45.1, 20.9, 9.84; 2.16, 2.12 | 0.973, 1.11, 1.17 | pass | 64.1, 29.8, 14.6, 7.82 | pass | reported |
| V-BO-07, temperature, ε_m = 0.01 | 1.64, 0.809, 0.68; 2.03, 1.19 | **fail** | 23.6, 11.8, 6.15; 2.0, 1.92 | 1.33, 1.47, 1.19 | pass | 54.3, 31.5, 20.9, 16.3 | pass | reported |
| V-BO-07, enthalpy, ε_m = 0.01 | 1.62, 1.41, 0.666; 1.15, 2.12 | **fail** | 23.6, 11.8, 6.14; 2.0, 1.93 | 1.33, 1.47, 1.19 | pass | 54.4, 31.6, 20.9, 16.3 | pass | reported |

## The runs

| case | core | form | ε_m | cells | retcode | steps (rejected) | invocations or parks | largest per-step rate before the last event (t) | peak QoI (t) | last event [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | old | old | NaN | 30 | Success | 648 (1) | 30 | 32.647 (359.17) | 32.254 (333.0) | 364.584 |
| V-BO-01 | old | old | NaN | 60 | Success | 853 (3) | 60 | 32.392 (334.67) | 32.392 (335.0) | 364.564 |
| V-BO-01 | old | old | NaN | 120 | Success | 1461 (0) | 120 | 32.469 (335.6) | 32.469 (335.5) | 364.608 |
| V-BO-01 | old | old | NaN | 240 | Success | 3277 (0) | 240 | 32.51 (336.1) | 32.51 (336.0) | 364.649 |
| V-BO-01 | new | temperature | 0.01 | 30 | Success | 685 (230) | 1 | 32.247 (333.51) | 32.248 (333.0) | 358.205 |
| V-BO-01 | new | temperature | 0.01 | 60 | Success | 1125 (225) | 1 | 32.391 (334.9) | 32.391 (335.0) | 357.915 |
| V-BO-01 | new | temperature | 0.01 | 120 | Success | 2083 (483) | 1 | 32.469 (335.63) | 32.469 (335.5) | 357.711 |
| V-BO-01 | new | temperature | 0.01 | 240 | Success | 4239 (1053) | 1 | 32.51 (336.02) | 32.51 (336.0) | 357.593 |
| V-BO-01 | new | temperature | 0.003 | 30 | Success | 744 (227) | 1 | 32.247 (332.83) | 32.248 (333.0) | 363.548 |
| V-BO-01 | new | temperature | 0.003 | 60 | Success | 1177 (201) | 1 | 32.391 (335.01) | 32.391 (335.0) | 363.25 |
| V-BO-01 | new | temperature | 0.003 | 120 | Success | 2339 (577) | 1 | 32.469 (335.62) | 32.469 (335.5) | 363.036 |
| V-BO-01 | new | temperature | 0.003 | 240 | Success | 4688 (1082) | 1 | 32.51 (335.98) | 32.51 (336.0) | 362.911 |
| V-BO-01 | new | temperature | 0.03 | 30 | Success | 733 (264) | 2 | 32.254 (333.18) | 32.254 (333.0) | 359.999 |
| V-BO-01 | new | temperature | 0.03 | 60 | Success | 1110 (229) | 2 | 32.392 (334.72) | 32.392 (335.0) | 359.763 |
| V-BO-01 | new | temperature | 0.03 | 120 | Success | 1976 (448) | 2 | 32.469 (335.56) | 32.469 (335.5) | 359.664 |
| V-BO-01 | new | temperature | 0.03 | 240 | Success | 3897 (1039) | 2 | 32.51 (335.95) | 32.51 (336.0) | 359.553 |
| V-BO-01 | new | enthalpy | 0.01 | 30 | Success | 918 (345) | 1 | 32.248 (332.99) | 32.248 (333.0) | 358.202 |
| V-BO-01 | new | enthalpy | 0.01 | 60 | Success | 1583 (581) | 1 | 32.391 (334.94) | 32.391 (335.0) | 357.912 |
| V-BO-01 | new | enthalpy | 0.01 | 120 | Success | 2710 (807) | 1 | 32.469 (335.7) | 32.469 (335.5) | 357.71 |
| V-BO-01 | new | enthalpy | 0.01 | 240 | Success | 4763 (1035) | 1 | 32.51 (336.02) | 32.51 (336.0) | 357.594 |
| V-BO-01 | new | enthalpy | 0.003 | 30 | Success | 959 (328) | 1 | 32.248 (333.32) | 32.248 (333.0) | 363.555 |
| V-BO-01 | new | enthalpy | 0.003 | 60 | Success | 1658 (617) | 1 | 32.391 (334.78) | 32.391 (335.0) | 363.255 |
| V-BO-01 | new | enthalpy | 0.003 | 120 | Success | 2912 (877) | 1 | 32.469 (335.62) | 32.469 (335.5) | 363.042 |
| V-BO-01 | new | enthalpy | 0.003 | 240 | Success | 5154 (1049) | 1 | 32.51 (336.0) | 32.51 (336.0) | 362.919 |
| V-BO-01 | new | enthalpy | 0.03 | 30 | Success | 925 (352) | 2 | 32.254 (333.33) | 32.254 (333.0) | 359.998 |
| V-BO-01 | new | enthalpy | 0.03 | 60 | Success | 1575 (627) | 2 | 32.392 (334.75) | 32.392 (335.0) | 359.759 |
| V-BO-01 | new | enthalpy | 0.03 | 120 | Success | 2634 (830) | 2 | 32.469 (335.56) | 32.469 (335.5) | 359.663 |
| V-BO-01 | new | enthalpy | 0.03 | 240 | Success | 4410 (896) | 2 | 32.51 (335.95) | 32.51 (336.0) | 359.554 |
| V-BO-03 | old | old | NaN | 30 | Success | 579 (1) | 26 | 30.398 (315.98) | 30.4 (316.5) | 649.317 |
| V-BO-03 | old | old | NaN | 60 | Success | 769 (3) | 52 | 30.386 (317.2) | 30.387 (317.0) | 1037.46 |
| V-BO-03 | old | old | NaN | 120 | Success | 1721 (0) | 103 | 30.383 (316.95) | 30.383 (317.0) | 2039.06 |
| V-BO-03 | old | old | NaN | 240 | Success | 3997 (0) | 207 | 30.383 (316.97) | 30.383 (317.0) | 2728.03 |
| V-BO-03 | new | temperature | 0.01 | 30 | Success | 674 (201) | 1 | 30.393 (316.21) | 30.394 (316.5) | 3022.46 |
| V-BO-03 | new | temperature | 0.01 | 60 | Success | 1131 (210) | 1 | 30.385 (316.83) | 30.385 (317.0) | 3022.87 |
| V-BO-03 | new | temperature | 0.01 | 120 | Success | 2141 (492) | 1 | 30.383 (317.0) | 30.383 (317.0) | 3022.97 |
| V-BO-03 | new | temperature | 0.01 | 240 | Success | 4301 (1051) | 1 | 30.383 (316.95) | 30.383 (317.0) | 3023.0 |
| V-BO-03 | new | enthalpy | 0.01 | 30 | Success | 927 (303) | 1 | 30.394 (316.84) | 30.394 (316.5) | 3022.46 |
| V-BO-03 | new | enthalpy | 0.01 | 60 | Success | 1611 (580) | 1 | 30.385 (316.92) | 30.385 (317.0) | 3022.87 |
| V-BO-03 | new | enthalpy | 0.01 | 120 | Success | 2810 (910) | 1 | 30.383 (317.02) | 30.383 (317.0) | 3022.97 |
| V-BO-03 | new | enthalpy | 0.01 | 240 | Success | 4909 (1079) | 1 | 30.383 (316.99) | 30.383 (317.0) | 3023.0 |
| V-BO-06 | old | old | NaN | 30 | Success | 609 (5) | 30 | 49.998 (360.28) | 36.649 (360.5) | 365.307 |
| V-BO-06 | old | old | NaN | 60 | Success | 896 (0) | 60 | 32.29 (334.66) | 32.292 (334.0) | 364.894 |
| V-BO-06 | old | old | NaN | 120 | Success | 1543 (0) | 120 | 32.415 (334.99) | 32.415 (335.0) | 364.765 |
| V-BO-06 | old | old | NaN | 240 | Success | 3678 (0) | 240 | 32.482 (335.73) | 32.481 (336.0) | 364.724 |
| V-BO-06 | new | temperature | 0.01 | 30 | Success | 754 (294) | 1 | 32.072 (331.77) | 32.072 (331.5) | 361.847 |
| V-BO-06 | new | temperature | 0.01 | 60 | Success | 1132 (214) | 1 | 32.288 (334.04) | 32.288 (334.0) | 361.164 |
| V-BO-06 | new | temperature | 0.01 | 120 | Success | 2178 (515) | 1 | 32.414 (335.18) | 32.414 (335.0) | 360.802 |
| V-BO-06 | new | temperature | 0.01 | 240 | Success | 4394 (1166) | 1 | 32.481 (335.79) | 32.481 (336.0) | 360.616 |
| V-BO-06 | new | enthalpy | 0.01 | 30 | Success | 934 (350) | 1 | 32.072 (331.79) | 32.072 (331.5) | 361.847 |
| V-BO-06 | new | enthalpy | 0.01 | 60 | Success | 1620 (606) | 1 | 32.288 (334.13) | 32.288 (334.0) | 361.161 |
| V-BO-06 | new | enthalpy | 0.01 | 120 | Success | 2941 (1052) | 1 | 32.414 (335.14) | 32.414 (335.0) | 360.801 |
| V-BO-06 | new | enthalpy | 0.01 | 240 | Success | 4852 (935) | 1 | 32.481 (335.74) | 32.481 (336.0) | 360.617 |
| V-BO-07 | old | old | NaN | 30 | Success | 604 (3) | 30 | 30.41 (286.8) | 30.412 (286.0) | 325.659 |
| V-BO-07 | old | old | NaN | 60 | Success | 794 (0) | 60 | 30.542 (287.45) | 30.542 (287.5) | 325.758 |
| V-BO-07 | old | old | NaN | 120 | Success | 1093 (2) | 120 | 30.609 (288.18) | 30.609 (288.0) | 325.818 |
| V-BO-07 | old | old | NaN | 240 | Success | 1848 (0) | 240 | 30.644 (288.34) | 30.644 (288.0) | 326.391 |
| V-BO-07 | new | temperature | 0.01 | 30 | Success | 506 (151) | 1 | 30.412 (286.0) | 30.412 (286.0) | 314.468 |
| V-BO-07 | new | temperature | 0.01 | 60 | Success | 724 (218) | 1 | 30.542 (287.36) | 30.542 (287.5) | 314.115 |
| V-BO-07 | new | temperature | 0.01 | 120 | Success | 1262 (432) | 1 | 30.609 (288.03) | 30.609 (288.0) | 313.945 |
| V-BO-07 | new | temperature | 0.01 | 240 | Success | 2256 (769) | 1 | 30.644 (288.24) | 30.643 (288.0) | 313.84 |
| V-BO-07 | new | enthalpy | 0.01 | 30 | Success | 717 (305) | 1 | 30.412 (285.96) | 30.412 (286.0) | 314.465 |
| V-BO-07 | new | enthalpy | 0.01 | 60 | Success | 1179 (596) | 1 | 30.542 (287.37) | 30.542 (287.5) | 314.113 |
| V-BO-07 | new | enthalpy | 0.01 | 120 | Success | 1859 (755) | 1 | 30.609 (287.95) | 30.609 (288.0) | 313.945 |
| V-BO-07 | new | enthalpy | 0.01 | 240 | Success | 2810 (774) | 1 | 30.644 (288.23) | 30.643 (288.0) | 313.838 |
