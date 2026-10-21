# Gate item 2: the end-of-layer burst by refinement

Produced by `test/newcore/gates/burnout_burst.jl` at commit 5f08134. Julia 1.13.0, 1 thread; load average [20.4, 18.8, 17.5] at the start and [20.4, 20.9, 22.6] at the end. The finite-volume core through the case's solve and the old core with the case's keywords, the gas-generation rate at every accepted step and post-event state [g m⁻² s⁻¹]; each deck's window from the old core's 240-cell peak time to its end (V-BO-01 [336.0, 500.0] s, V-BO-03 [317.0, 3100.0] s, V-BO-06 [336.0, 500.0] s, V-BO-07 [288.0, 500.0] s; V-BO-03 is run to 3100.0 s). (a) The largest per-step excess of N cells over 2N after the window's start, falling at least 1.5× per doubling and no larger than the old core's; (b) the L1 of successive grids [g/m²], falling at least 1.5× per doubling, and the L1 against the family's 240-cell run within 1.5× of the old core's; (c) the L1 against the old core's 240-cell run falling with refinement. Gated: V-BO-01 at ε_m of 0.01 and 0.003 in both forms; the rest is reported.

## The criteria

Each measure lists 30→60, 60→120 and 120→240 (the L1 against the old core's 240-cell run: 30, 60, 120 and 240), then its ratios per doubling.

| family | (a) excess of N over 2N; ratios | (a) | (b) successive L1; ratios | (b) self-L1 over the old core's | (b) | (c) L1 against the old core's 240 | (c) | gated |
|---|---|---|---|---|---|---|---|---|
| V-BO-01, old core | 31.9, 7.73, 3.18; 4.13, 2.43 | | 25.0, 9.34, 4.4; 2.68, 2.12 | | | | | |
| V-BO-03, old core | 0.0115, 0.00395, 0.000778; 2.92, 5.08 | | 9.75, 2.53, 0.661; 3.86, 3.82 | | | | | |
| V-BO-06, old core | 114.0, 25.5, 7.29; 4.47, 3.5 | | 51.7, 19.9, 8.44; 2.6, 2.36 | | | | | |
| V-BO-07, old core | 12.0, 3.67, 0.896; 3.28, 4.1 | | 21.0, 9.05, 5.16; 2.32, 1.76 | | | | | |
| V-BO-01, temperature, ε_m = 0.01 | 3.16, 1.36, 0.642; 2.33, 2.12 | pass | 22.0, 10.2, 5.0; 2.15, 2.05 | 1.03, 1.15, 1.14 | pass | 39.7, 19.4, 10.0, 5.93 | pass | yes |
| V-BO-01, temperature, ε_m = 0.003 | 19.1, 5.14, 2.19; 3.72, 2.35 | pass | 22.7, 11.7, 5.86; 1.94, 1.99 | 1.02, 1.27, 1.33 | pass | 36.7, 18.4, 8.21, 3.21 | pass | yes |
| V-BO-01, temperature, ε_m = 0.03 | 4.59, 1.8, 0.909; 2.55, 1.98 | pass | 21.8, 9.95, 5.21; 2.19, 1.91 | 1.01, 1.14, 1.18 | pass | 38.5, 18.4, 8.9, 4.71 | pass | reported |
| V-BO-01, enthalpy, ε_m = 0.01 | 3.17, 1.35, 0.636; 2.34, 2.13 | pass | 22.0, 10.2, 5.01; 2.15, 2.04 | 1.03, 1.15, 1.14 | pass | 39.8, 19.4, 10.0, 5.92 | pass | yes |
| V-BO-01, enthalpy, ε_m = 0.003 | 18.7, 4.98, 2.1; 3.75, 2.38 | pass | 22.7, 11.6, 5.81; 1.95, 2.0 | 1.02, 1.26, 1.32 | pass | 36.8, 18.4, 8.33, 3.23 | pass | yes |
| V-BO-01, enthalpy, ε_m = 0.03 | 4.6, 1.8, 0.901; 2.55, 2.0 | pass | 21.8, 9.93, 5.18; 2.19, 1.92 | 1.01, 1.14, 1.18 | pass | 38.5, 18.5, 8.91, 4.7 | pass | reported |
| V-BO-03, temperature, ε_m = 0.01 | 0.00693, 0.00213, 0.000545; 3.26, 3.9 | pass | 8.68, 2.29, 0.739; 3.8, 3.1 | 0.888, 0.885, 1.12 | pass | 11.1, 2.6, 0.727, 0.488 | pass | reported |
| V-BO-03, enthalpy, ε_m = 0.01 | 0.0068, 0.00214, 0.000546; 3.18, 3.91 | pass | 8.68, 2.26, 0.574; 3.85, 3.93 | 0.895, 0.898, 0.869 | pass | 11.4, 2.75, 0.591, 0.237 | pass | reported |
| V-BO-06, temperature, ε_m = 0.01 | 25.6, 5.23, 2.08; 4.89, 2.52 | pass | 45.0, 20.9, 9.91; 2.15, 2.11 | 0.974, 1.12, 1.17 | pass | 64.1, 29.9, 14.5, 7.85 | pass | reported |
| V-BO-06, enthalpy, ε_m = 0.01 | 25.6, 5.21, 2.08; 4.91, 2.51 | pass | 45.1, 20.9, 9.84; 2.16, 2.12 | 0.973, 1.11, 1.17 | pass | 64.1, 29.8, 14.6, 7.83 | pass | reported |
| V-BO-07, temperature, ε_m = 0.01 | 1.61, 0.799, 1.99; 2.02, 0.402 | **fail** | 23.5, 11.8, 6.11; 1.99, 1.93 | 1.33, 1.47, 1.18 | pass | 54.4, 31.6, 20.9, 16.3 | pass | reported |
| V-BO-07, enthalpy, ε_m = 0.01 | 1.6, 0.797, 1.18; 2.01, 0.677 | **fail** | 23.5, 11.8, 6.1; 1.99, 1.93 | 1.33, 1.46, 1.18 | pass | 54.4, 31.6, 20.9, 16.3 | pass | reported |

## The runs

| case | core | form | ε_m | cells | retcode | steps (rejected) | invocations or parks | largest per-step rate before the last event (t) | peak QoI (t) | last event [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | old | old | NaN | 30 | Success | 648 (1) | 30 | 32.647 (359.17) | 32.254 (333.0) | 364.584 |
| V-BO-01 | old | old | NaN | 60 | Success | 853 (3) | 60 | 32.392 (334.67) | 32.392 (335.0) | 364.564 |
| V-BO-01 | old | old | NaN | 120 | Success | 1461 (0) | 120 | 32.469 (335.6) | 32.469 (335.5) | 364.608 |
| V-BO-01 | old | old | NaN | 240 | Success | 3277 (0) | 240 | 32.51 (336.1) | 32.51 (336.0) | 364.649 |
| V-BO-01 | new | temperature | 0.01 | 30 | Success | 686 (220) | 1 | 32.247 (332.77) | 32.248 (333.0) | 358.205 |
| V-BO-01 | new | temperature | 0.01 | 60 | Success | 1127 (225) | 1 | 32.391 (334.88) | 32.391 (335.0) | 357.915 |
| V-BO-01 | new | temperature | 0.01 | 120 | Success | 2118 (493) | 1 | 32.469 (335.52) | 32.469 (335.5) | 357.711 |
| V-BO-01 | new | temperature | 0.01 | 240 | Success | 4203 (1042) | 1 | 32.51 (335.99) | 32.51 (336.0) | 357.593 |
| V-BO-01 | new | temperature | 0.003 | 30 | Success | 745 (223) | 1 | 32.248 (333.24) | 32.248 (333.0) | 363.548 |
| V-BO-01 | new | temperature | 0.003 | 60 | Success | 1188 (201) | 1 | 32.391 (335.04) | 32.391 (335.0) | 363.25 |
| V-BO-01 | new | temperature | 0.003 | 120 | Success | 2323 (570) | 1 | 32.469 (335.53) | 32.469 (335.5) | 363.036 |
| V-BO-01 | new | temperature | 0.003 | 240 | Success | 4692 (1089) | 1 | 32.51 (336.05) | 32.51 (336.0) | 362.911 |
| V-BO-01 | new | temperature | 0.03 | 30 | Success | 725 (262) | 2 | 32.254 (332.99) | 32.254 (333.0) | 359.999 |
| V-BO-01 | new | temperature | 0.03 | 60 | Success | 1098 (224) | 2 | 32.392 (334.92) | 32.392 (335.0) | 359.763 |
| V-BO-01 | new | temperature | 0.03 | 120 | Success | 1981 (452) | 2 | 32.469 (335.69) | 32.469 (335.5) | 359.664 |
| V-BO-01 | new | temperature | 0.03 | 240 | Success | 3904 (1041) | 2 | 32.51 (335.98) | 32.51 (336.0) | 359.553 |
| V-BO-01 | new | enthalpy | 0.01 | 30 | Success | 919 (319) | 1 | 32.247 (332.84) | 32.248 (333.0) | 358.202 |
| V-BO-01 | new | enthalpy | 0.01 | 60 | Success | 1580 (613) | 1 | 32.391 (334.92) | 32.391 (335.0) | 357.912 |
| V-BO-01 | new | enthalpy | 0.01 | 120 | Success | 2789 (860) | 1 | 32.469 (335.66) | 32.469 (335.5) | 357.71 |
| V-BO-01 | new | enthalpy | 0.01 | 240 | Success | 4743 (1006) | 1 | 32.51 (335.97) | 32.51 (336.0) | 357.594 |
| V-BO-01 | new | enthalpy | 0.003 | 30 | Success | 965 (313) | 1 | 32.248 (333.29) | 32.248 (333.0) | 363.555 |
| V-BO-01 | new | enthalpy | 0.003 | 60 | Success | 1634 (569) | 1 | 32.391 (334.88) | 32.391 (335.0) | 363.255 |
| V-BO-01 | new | enthalpy | 0.003 | 120 | Success | 2946 (929) | 1 | 32.469 (335.7) | 32.469 (335.5) | 363.042 |
| V-BO-01 | new | enthalpy | 0.003 | 240 | Success | 5066 (984) | 1 | 32.51 (335.99) | 32.51 (336.0) | 362.919 |
| V-BO-01 | new | enthalpy | 0.03 | 30 | Success | 927 (323) | 2 | 32.254 (333.32) | 32.254 (333.0) | 359.997 |
| V-BO-01 | new | enthalpy | 0.03 | 60 | Success | 1565 (586) | 2 | 32.392 (334.65) | 32.392 (335.0) | 359.759 |
| V-BO-01 | new | enthalpy | 0.03 | 120 | Success | 2644 (846) | 2 | 32.469 (335.55) | 32.469 (335.5) | 359.663 |
| V-BO-01 | new | enthalpy | 0.03 | 240 | Success | 4437 (955) | 2 | 32.51 (336.04) | 32.51 (336.0) | 359.553 |
| V-BO-03 | old | old | NaN | 30 | Success | 579 (1) | 26 | 30.398 (315.98) | 30.4 (316.5) | 649.317 |
| V-BO-03 | old | old | NaN | 60 | Success | 769 (3) | 52 | 30.386 (317.2) | 30.387 (317.0) | 1037.46 |
| V-BO-03 | old | old | NaN | 120 | Success | 1721 (0) | 103 | 30.383 (316.95) | 30.383 (317.0) | 2039.06 |
| V-BO-03 | old | old | NaN | 240 | Success | 3997 (0) | 207 | 30.383 (316.97) | 30.383 (317.0) | 2728.03 |
| V-BO-03 | new | temperature | 0.01 | 30 | Success | 678 (211) | 1 | 30.394 (316.48) | 30.394 (316.5) | 3022.46 |
| V-BO-03 | new | temperature | 0.01 | 60 | Success | 1121 (213) | 1 | 30.385 (316.62) | 30.385 (317.0) | 3022.87 |
| V-BO-03 | new | temperature | 0.01 | 120 | Success | 2168 (501) | 1 | 30.383 (317.04) | 30.383 (317.0) | 3022.97 |
| V-BO-03 | new | temperature | 0.01 | 240 | Success | 4294 (1047) | 1 | 30.383 (317.01) | 30.383 (317.0) | 3023.0 |
| V-BO-03 | new | enthalpy | 0.01 | 30 | Success | 940 (329) | 1 | 30.394 (316.59) | 30.394 (316.5) | 3022.46 |
| V-BO-03 | new | enthalpy | 0.01 | 60 | Success | 1610 (576) | 1 | 30.385 (317.03) | 30.385 (317.0) | 3022.87 |
| V-BO-03 | new | enthalpy | 0.01 | 120 | Success | 2734 (758) | 1 | 30.383 (317.01) | 30.383 (317.0) | 3022.97 |
| V-BO-03 | new | enthalpy | 0.01 | 240 | Success | 4888 (1062) | 1 | 30.383 (316.99) | 30.383 (317.0) | 3023.0 |
| V-BO-06 | old | old | NaN | 30 | Success | 609 (5) | 30 | 49.998 (360.28) | 36.649 (360.5) | 365.307 |
| V-BO-06 | old | old | NaN | 60 | Success | 896 (0) | 60 | 32.29 (334.66) | 32.292 (334.0) | 364.894 |
| V-BO-06 | old | old | NaN | 120 | Success | 1543 (0) | 120 | 32.415 (334.99) | 32.415 (335.0) | 364.765 |
| V-BO-06 | old | old | NaN | 240 | Success | 3678 (0) | 240 | 32.482 (335.73) | 32.481 (336.0) | 364.724 |
| V-BO-06 | new | temperature | 0.01 | 30 | Success | 753 (294) | 1 | 32.072 (331.7) | 32.072 (331.5) | 361.848 |
| V-BO-06 | new | temperature | 0.01 | 60 | Success | 1137 (214) | 1 | 32.288 (334.06) | 32.288 (334.0) | 361.164 |
| V-BO-06 | new | temperature | 0.01 | 120 | Success | 2153 (511) | 1 | 32.414 (335.15) | 32.414 (335.0) | 360.802 |
| V-BO-06 | new | temperature | 0.01 | 240 | Success | 4388 (1167) | 1 | 32.481 (335.73) | 32.481 (336.0) | 360.616 |
| V-BO-06 | new | enthalpy | 0.01 | 30 | Success | 936 (351) | 1 | 32.072 (331.7) | 32.072 (331.5) | 361.847 |
| V-BO-06 | new | enthalpy | 0.01 | 60 | Success | 1639 (671) | 1 | 32.288 (333.91) | 32.288 (334.0) | 361.161 |
| V-BO-06 | new | enthalpy | 0.01 | 120 | Success | 2915 (1019) | 1 | 32.414 (335.16) | 32.414 (335.0) | 360.801 |
| V-BO-06 | new | enthalpy | 0.01 | 240 | Success | 4656 (802) | 1 | 32.481 (335.79) | 32.481 (336.0) | 360.617 |
| V-BO-07 | old | old | NaN | 30 | Success | 604 (3) | 30 | 30.41 (286.8) | 30.412 (286.0) | 325.659 |
| V-BO-07 | old | old | NaN | 60 | Success | 794 (0) | 60 | 30.542 (287.45) | 30.542 (287.5) | 325.758 |
| V-BO-07 | old | old | NaN | 120 | Success | 1093 (2) | 120 | 30.609 (288.18) | 30.609 (288.0) | 325.818 |
| V-BO-07 | old | old | NaN | 240 | Success | 1848 (0) | 240 | 30.644 (288.34) | 30.644 (288.0) | 326.391 |
| V-BO-07 | new | temperature | 0.01 | 30 | Success | 509 (154) | 1 | 30.412 (286.0) | 30.412 (286.0) | 314.415 |
| V-BO-07 | new | temperature | 0.01 | 60 | Success | 732 (215) | 1 | 30.542 (287.46) | 30.542 (287.5) | 314.091 |
| V-BO-07 | new | temperature | 0.01 | 120 | Success | 1277 (436) | 1 | 30.609 (288.04) | 30.609 (288.0) | 313.938 |
| V-BO-07 | new | temperature | 0.01 | 240 | Success | 2248 (767) | 1 | 30.644 (288.24) | 30.643 (288.0) | 313.842 |
| V-BO-07 | new | enthalpy | 0.01 | 30 | Success | 719 (322) | 1 | 30.412 (286.26) | 30.412 (286.0) | 314.413 |
| V-BO-07 | new | enthalpy | 0.01 | 60 | Success | 1150 (524) | 1 | 30.542 (287.33) | 30.542 (287.5) | 314.09 |
| V-BO-07 | new | enthalpy | 0.01 | 120 | Success | 1862 (785) | 1 | 30.609 (287.97) | 30.609 (288.0) | 313.937 |
| V-BO-07 | new | enthalpy | 0.01 | 240 | Success | 2806 (757) | 1 | 30.644 (288.23) | 30.643 (288.0) | 313.84 |
