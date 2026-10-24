# The verification cases on both cores

Produced by `case_comparison.jl` at commit e4ac985 (tree facf693). Julia 1.13.0, one thread per process. Accuracy: the old core at a one-minute load of 0.91, the finite-volume core at 0.91. Timing: one process on CPUs `1`, one-minute load 1.27 at the start and 1.11 at the end, 345 s. The script's header defines every measure.

## Per case

Quantities passing on each core; on the exact-solution cases, the geometric mean of the finite-volume core's error over the old core's on the quantities both report with both errors above round-off (1.0e-12, relative to the exact value where the case gives one; below one: the finite-volume core closer to the exact solution), and their count (the burnout cases' accuracy is `cutover_comparison.md`'s); the convergence check's observed orders; each core's accepted (rejected) steps and Jacobians (blank where a case's solve returns no integrator statistics); the warm medians [s] and their ratio.

| case | cells | passing old, new | error ratio (n) | orders old | orders new | steps old | steps new | Jacobians old, new | old [s] | new [s] | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | 60 | 7/7, 8/8 |  (0) |  |  | 853 (3) | 985 (242) | 164, 985 |  |  |  |
| V-BO-02 | 60 | 10/10, 10/10 |  (0) |  |  | 1572 (0) | 1978 (288) | 341, 1978 |  |  |  |
| V-BO-03 | 60 | 3/4, 4/4 |  (0) |  |  | 681 (3) | 928 (231) | 126, 928 |  |  |  |
| V-BO-04 | 40 | 4/4, 5/5 |  (0) |  |  | 353 (0) | 322 (0) | 77, 322 |  |  |  |
| V-BO-05 | 60 | 3/4, 4/4 |  (0) |  |  | 169 (9) | 347 (143) | 86, 347 |  |  |  |
| V-BO-06 | 60 | 7/7, 8/8 |  (0) |  |  | 896 (0) | 1040 (236) | 169, 1040 |  |  |  |
| V-BO-07 | 60 | 7/7, 8/8 |  (0) |  |  | 794 (0) | 785 (227) | 181, 785 |  |  |  |
| V-BO-08 | 10 | 4/4, 6/6 |  (0) |  |  | 291 (1) | 221 (51) | 65, 221 |  |  |  |
| V-CPL-01 | 150 | 9/9, 9/9 | 1.0 (9) |  |  | 319 (0) | 246 (24) | 157, 246 | 0.0513 | 0.0368 | 0.717 |
| V-CPL-01b | 150 | 9/9, 9/9 | 1.0 (9) |  |  | 640 (69) | 1141 (241) | 167, 1142 | 0.228 | 0.251 | 1.1 |
| V-CPL-02 | 160 | 6/6, 6/6 | 1.18 (6) |  |  | 1037 (0) | 1053 (41) | 497, 1054 | 0.494 | 0.587 | 1.19 |
| V-CPL-02b | 100 | 2/2, 2/2 | 1.88 (2) |  |  | 615 (0) | 639 (2) | 288, 639 | 0.194 | 0.227 | 1.17 |
| V-CPL-03 | 8 | 18/18, 18/18 | 0.909 (9) |  |  | 85 (0) | 88 (0) | 15, 88 | 0.00238 | 0.00267 | 1.12 |
| V-CPL-04 | 40 | 6/6, 6/6 | 1.0 (5) | 2.1, 2.47 | 2.1, 2.47 | 29 (0) | 28 (0) | 5, 28 | 0.00123 | 0.00163 | 1.33 |
| V-CPL-05 | 16 | 4/4, 4/4 | 0.999 (3) |  |  | 343 (26) | 237 (31) | 17, 237 | 0.00976 | 0.00841 | 0.862 |
| V-GT-01 | 160 | 7/7, 7/7 | 1.0 (5) | 2.01, 2.0 | 2.01, 2.0 | 76 (0) | 46 (0) | 23, 46 | 0.0253 | 0.012 | 0.476 |
| V-GT-02 | 80 | 4/4, 4/4 | 1.0 (3) | 0.892, 0.938 | 0.892, 0.938 | 44 (0) | 52 (4) | 2, 52 | 0.0075 | 0.01 | 1.34 |
| V-HC-01 | 100 | 9/9, 9/9 | 1.0 (9) | 1.99, 2.0 | 1.99, 2.0 | 56 (0) | 60 (0) | 11, 60 | 0.00642 | 0.00588 | 0.916 |
| V-HC-02a | 80 | 8/8, 8/8 | 1.0 (8) |  |  | 81 (1) | 56 (0) | 5, 56 | 0.00522 | 0.00421 | 0.806 |
| V-HC-02b | 80 | 8/8, 8/8 | 1.0 (8) |  |  | 115 (0) | 75 (0) | 5, 75 | 0.0091 | 0.00529 | 0.582 |
| V-HC-02c | 80 | 8/8, 8/8 | 1.0 (8) | 2.0, 2.0 | 2.0, 2.0 | 149 (0) | 98 (0) | 11, 98 | 0.0122 | 0.00676 | 0.554 |
| V-HC-02d | 80 | 8/8, 8/8 | 1.0 (8) |  |  | 164 (0) | 113 (0) | 20, 113 | 0.0141 | 0.00764 | 0.543 |
| V-HC-03 | 80 | 5/5, 5/5 | 0.646 (4) | 2.0, 2.0 | 2.0, 2.0 | 86 (0) | 86 (0) | 7, 86 | 0.00709 | 0.00643 | 0.908 |
| V-HC-04a | 40 | 4/4, 4/4 |  (0) |  |  | 112 (0) | 99 (0) | 14, 99 | 0.00675 | 0.00626 | 0.928 |
| V-HC-04b | 40 | 4/4, 4/4 |  (0) |  |  | 155 (0) | 97 (0) | 8, 97 | 0.00833 | 0.006 | 0.72 |
| V-HC-04c | 40 | 5/5, 5/5 | 1.0 (1) |  |  | 112 (0) | 99 (0) | 14, 99 | 0.00677 | 0.00618 | 0.912 |
| V-HC-05 | 80 | 3/3, 3/3 | 1.0 (2) | 1.97, 1.98 | 1.97, 1.98 | 138 (0) | 243 (0) | 26, 243 | 0.0129 | 0.0145 | 1.12 |
| V-MMS-01 | 40 | 5/5, 5/5 | 1.11 (5) | 1.98, 1.99 | 1.98, 2.0 |  |  | ,  | 0.118 | 0.177 | 1.51 |
| V-MMS-02 | 20 | 4/4, 4/4 | 0.0181 (1) |  |  |  |  | ,  | 0.00746 | 0.00576 | 0.773 |
| V-MMS-03 | 32 | 3/3, 3/3 | 0.312 (3) | 2.21, 2.22 | 3.51, 1.61 |  |  | ,  | 0.0545 | 0.0366 | 0.671 |
| V-MMS-04 | 32 | 6/6, 4/4 | 1.0 (4) | 2.0, 2.01 | 2.0, 2.0 |  |  | ,  | 0.148 | 0.138 | 0.935 |
| V-RAD-01 | 100 | 11/11, 11/11 | 0.999 (10) | 2.0, 2.0 | 2.0, 2.0 | 61 (1) | 63 (0) | 10, 63 | 0.0183 | 0.0069 | 0.377 |
| V-RX-01 | 1 | 5/5, 5/5 | 0.000434 (1) |  |  | 185 (0) | 175 (5) | 30, 175 | 0.00211 | 0.00123 | 0.582 |
| V-RX-02a | 1 | 2/2, 2/2 |  (0) |  |  | 6 (0) | 13 (0) | 5, 13 | 0.000361 | 0.000656 | 1.82 |
| V-RX-02b | 1 | 2/2, 2/2 | 0.00017 (1) |  |  | 104 (0) | 242 (0) | 33, 242 | 0.00105 | 0.00128 | 1.22 |
| V-RX-02c | 1 | 2/2, 2/2 | 0.816 (1) |  |  | 86 (11) | 183 (11) | 13, 183 | 0.0009 | 0.00108 | 1.2 |
| V-RX-02d | 1 | 2/2, 2/2 |  (0) |  |  | 144 (15) | 306 (29) | 18, 306 | 0.00118 | 0.00152 | 1.29 |
| V-RX-03 | 1 | 6/6, 6/6 | 0.0279 (2) |  |  | 143 (0) | 160 (0) | 17, 160 | 0.00189 | 0.00174 | 0.921 |
| V-RX-04a | 4 | 2/2, 2/2 | 1.04 (2) |  |  | 32 (0) | 78 (7) | 14, 78 | 0.000613 | 0.000872 | 1.42 |
| V-RX-04b | 20 | 3/3, 3/3 | 0.168 (2) |  |  | 122 (0) | 159 (0) | 10, 159 | 0.00332 | 0.0031 | 0.934 |

Over the exact-solution cases: the geometric mean of the error ratios 0.768 over 131 quantities (median 1.0), and of the time ratios 0.91 over 32 cases (from 0.377 to 1.82).

## The convergence checks

Each case's own error norm at the check's resolutions on each core, and the ratio of the finite-volume core's to the old core's.

| case | cells | old errors | new errors | ratios | expected order |
|---|---|---|---|---|---|
| V-CPL-04 | 20, 40, 80 | 0.0109, 0.00256, 0.00046 | 0.0109, 0.00256, 0.00046 | 1.0, 1.0, 1.0 | 2.0 |
| V-GT-01 | 80, 160, 320 | 3.23e-5, 8.02e-6, 2.0e-6 | 3.23e-5, 8.02e-6, 2.0e-6 | 1.0, 1.0, 1.0 | 2.0 |
| V-GT-02 | 40, 80, 160 | 0.00131, 0.000708, 0.00037 | 0.00131, 0.000708, 0.00037 | 1.0, 1.0, 1.0 | 1.0 |
| V-HC-01 | 40, 80, 160 | 0.312, 0.0785, 0.0197 | 0.312, 0.0785, 0.0197 | 1.0, 1.0, 1.0 | 2.0 |
| V-HC-02c | 40, 80, 160 | 0.013, 0.00325, 0.000812 | 0.013, 0.00325, 0.000812 | 1.0, 1.0, 1.0 | 2.0 |
| V-HC-03 | 40, 80, 160 | 0.00594, 0.00148, 0.000371 | 0.00594, 0.00148, 0.000371 | 1.0, 1.0, 1.0 | 2.0 |
| V-HC-05 | 40, 80, 160 | 0.000151, 3.87e-5, 9.79e-6 | 0.000151, 3.87e-5, 9.79e-6 | 1.0, 1.0, 1.0 | 2.0 |
| V-MMS-01 | 20, 40, 80 | 0.0276, 0.00699, 0.00176 | 0.0253, 0.00641, 0.00161 | 0.918, 0.917, 0.915 | 2.0 |
| V-MMS-03 | 16, 32, 64 | 0.234, 0.0508, 0.0109 | 0.0504, 0.00443, 0.00145 | 0.215, 0.0873, 0.133 | 2.0 |
| V-MMS-04 | 16, 32, 64 | 0.000753, 0.000188, 4.67e-5 | 0.000753, 0.000188, 4.7e-5 | 1.0, 1.0, 1.01 | 2.0 |
| V-RAD-01 | 50, 100, 200 | 0.0646, 0.0162, 0.00405 | 0.0646, 0.0162, 0.00405 | 1.0, 1.0, 1.0 | 2.0 |

## Every quantity

As each case reports it: the value on each core, the exact or reference value (blank where the case has none), each core's error and the case's tolerance; failing values in bold.

| case | quantity | old | new | exact | old error | new error | tol |
|---|---|---|---|---|---|---|---|
| V-BO-01 | complete burnout: live thickness at t_end is the substrate | 0.0119967 | 0.0119967 | 0.012 | 3.33e-6 | 3.33e-6 | 1.0e-5 |
| V-BO-01 | mass loss of the gasified PMMA layer | 69.9206 | 69.9206 | 69.9148 | 0.00584 | 0.00584 | 0.3 |
| V-BO-01 | peak MLR | 32.3921 | 32.3908 | 32.39 | 6.37e-5 | 2.35e-5 | 0.03 |
| V-BO-01 | mass ledger closure | 0.000477815 | 6.1682e-8 | 0 | 0.000478 | 6.17e-8 | 0.001 |
| V-BO-01 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-01 | time of peak MLR | 335.0 | 335.0 | 335.0 | 0 | 0 | 0.03 |
| V-BO-01 | burnout time (live thickness < 12.02 mm) | 363.5 | 363.0 | 363.5 | 0 | 0.00138 | 0.01 |
| V-BO-01 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-01 | deactivation invocations for the PMMA layer (reported) |  | 1.0 | 1.0 |  | 0 | Inf |
| V-BO-02 | complete burnout: live thickness at t_end is the substrate | 0.0119967 | 0.0119967 | 0.012 | 3.33e-6 | 3.33e-6 | 1.0e-5 |
| V-BO-02 | mass loss of the gasified PMMA layer | 69.9206 | 69.9206 | 69.9148 | 0.00584 | 0.00583 | 0.3 |
| V-BO-02 | peak MLR | 24.6226 | 23.755 | 24.62 | 0.000104 | 0.0351 | 0.25 |
| V-BO-02 | mass ledger closure | 2.29638e-5 | 6.12771e-5 | 0 | 2.3e-5 | 6.13e-5 | 0.001 |
| V-BO-02 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-02 | park events | 60.0 | 60.0 | 60.0 | 0 | 0 | 0.5 |
| V-BO-02 | maximum live-cell temperature | 915.418 | 912.628 |  | 915.0 | 913.0 | 1000.0 |
| V-BO-02 | smallest accepted step | 1.0e-6 | 4.98607e-5 | 1.0e-6 | 2.52e-9 | 0 | 0.001 |
| V-BO-02 | Newton failures | 56.0 | 0 | 0 | 56.0 | 0 | 100.0 |
| V-BO-02 | accepted steps | 1572.0 | 1978.0 | 0 | 1570.0 | 1980.0 | 2000.0 |
| V-BO-02 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-03 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-03 | transmitted-flux warning fired | 0 |  | 1.0 | **1.0** |  | 0.5 |
| V-BO-03 | peak transmitted fraction of the incident flux | 0.777013 | 0.776995 | 1.0 | 0.223 | 0.223 | 0.3 |
| V-BO-03 | mass loss (informational) | 68.8788 | 68.8787 | 69.9148 | 1.04 | 1.04 | 10.0 |
| V-BO-03 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-03 | transmitted flux reported (series peak above 1 % of the entering flux) |  | 0.776995 | 1.0 |  | 0 | 0.5 |
| V-BO-04 | thickness follows L₀ exp(−r t/ρ) through the park cascades | 8.22975e-5 |  | 8.22975e-5 | 1.84e-7 |  | 0.0001 |
| V-BO-04 | solid concentration stays at ρ_s (dilation removes volume only) | 1000.0 | 1000.0 | 1000.0 | 2.27e-16 | 8.35e-12 | 1.0e-8 |
| V-BO-04 | isothermal: T stays at T₀ | 600.0 | 600.0 | 600.0 | 1.14e-13 | 1.25e-12 | 1.0e-6 |
| V-BO-04 | park cascade: 40 → 5 live cells | 5.0 |  | 5.0 | 0 |  | 0.5 |
| V-BO-04 | thickness follows L₀ exp(−r t/ρ) through the event |  | 8.22975e-5 | 8.22975e-5 |  | 1.27e-9 | 1.0e-6 |
| V-BO-04 | one invocation deactivates the column into one survivor |  | 1.0 | 1.0 |  | 0 | 0.5 |
| V-BO-04 | the event is located at the exact crossing of ε_m |  | 460.517 | 460.517 |  | 1.57e-10 | 1.0e-6 |
| V-BO-05 | stop time at the first interior void | 336.988 |  | 337.5 | 0.512 |  | 7.5 |
| V-BO-05 | interior-void warning | 0 |  | 1.0 | **1.0** |  | 0.5 |
| V-BO-05 | park events before the stop | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-05 | maximum live-cell temperature up to the stop | 911.988 | 693.379 |  | 912.0 | 693.0 | 1000.0 |
| V-BO-05 | stop time at the first interior crossing of ε_m |  | 259.791 | 259.57 |  | 0.221 | 0.7 |
| V-BO-05 | outcome :interior |  | 1.0 | 1.0 |  | 0 | 0.5 |
| V-BO-05 | deactivation events before the stop |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-06 | complete burnout: live thickness at t_end is the substrate | 0.011981 | 0.011981 | 0.011981 | 2.59e-12 | 4.03e-10 | 1.0e-5 |
| V-BO-06 | mass loss of the gasified PMMA layer | 69.9481 | 69.9481 | 69.9481 | 4.75e-8 | 4.78e-6 | 0.3 |
| V-BO-06 | peak MLR | 32.2925 | 32.2882 | 32.39 | 0.00301 | 0.00314 | 0.03 |
| V-BO-06 | mass ledger closure | 0.00382521 | 9.43722e-8 | 0 | 0.00383 | 9.44e-8 | 0.005 |
| V-BO-06 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-06 | time of peak MLR | 334.0 | 334.0 | 335.0 | 0.00299 | 0.00299 | 0.03 |
| V-BO-06 | burnout time (live thickness < 12.0 mm) | 364.0 | 364.0 | 363.5 | 0.00138 | 0.00138 | 0.01 |
| V-BO-06 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-06 | deactivation invocations for the PMMA layer (reported) |  | 1.0 | 1.0 |  | 0 | Inf |
| V-BO-07 | complete burnout: live thickness at t_end is the substrate | 0.0119967 | 0.0119967 | 0.0119967 | 5.09e-11 | 1.9e-9 | 1.0e-5 |
| V-BO-07 | mass loss of the gasified PMMA layer | 65.3383 | 65.3383 | 65.3383 | 4.21e-7 | 1.48e-5 | 0.3 |
| V-BO-07 | peak MLR | 30.5418 | 30.5418 | 30.64 | 0.0032 | 0.0032 | 0.03 |
| V-BO-07 | mass ledger closure | 3.06029e-5 | 6.09075e-8 | 0 | 3.06e-5 | 6.09e-8 | 0.001 |
| V-BO-07 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-07 | time of peak MLR | 287.5 | 287.5 | 288.0 | 0.00174 | 0.00174 | 0.03 |
| V-BO-07 | burnout time (live thickness < 12.02 mm) | 325.0 | 324.0 | 325.0 | 0 | 0.00308 | 0.01 |
| V-BO-07 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-07 | deactivation invocations for the PMMA layer (reported) |  | 1.0 | 1.0 |  | 0 | Inf |
| V-BO-08 | complete burnout: live thickness at t_end is the substrate | 0.0119967 | 0.0119967 | 0.0119967 | 9.32e-12 | 4.48e-10 | 1.0e-5 |
| V-BO-08 | mass loss of the gasified PMMA layer | 27.9239 | 27.9239 | 27.9239 | 2.61e-7 | 1.27e-5 | 0.3 |
| V-BO-08 | mass ledger closure | 0.000317492 | 5.73086e-6 | 0 | 0.000317 | 5.73e-6 | 0.001 |
| V-BO-08 | surface-Newton warnings | 0 |  | 0 | 0 |  | 0.5 |
| V-BO-08 | clamped surface balances on accepted states |  | 0 | 0 |  | 0 | 0.5 |
| V-BO-08 | deactivation invocations for the PMMA layer (reported) |  | 1.0 | 1.0 |  | 0 | Inf |
| V-BO-08 | the correction is zero at every Kaowool node |  | 0 | 0 |  | 0 | 2.22e-16 |
| V-CPL-01 | T-history L∞ error at x = 0.5 mm | 0.0137082 | 0.0135896 |  | 0.0137 | 0.0136 | 0.05 |
| V-CPL-01 | T-history L∞ error at x = 1.0 mm | 0.0272804 | 0.0270525 |  | 0.0273 | 0.0271 | 0.1 |
| V-CPL-01 | T-history L∞ error at x = 1.5 mm | 0.0376655 | 0.0373394 |  | 0.0377 | 0.0373 | 0.12 |
| V-CPL-01 | T-history L∞ error at x = 2.0 mm | 0.133802 | 0.133687 |  | 0.134 | 0.134 | 0.4 |
| V-CPL-01 | front position at t = 100.0 s | 0.00169157 | 0.00169168 | 0.00167516 | 0.0098 | 0.00986 | 0.03 |
| V-CPL-01 | front position at t = 200.0 s | 0.00239072 | 0.00239084 | 0.00236903 | 0.00916 | 0.00921 | 0.03 |
| V-CPL-01 | front position at t = 300.0 s | 0.00292913 | 0.00292927 | 0.00290145 | 0.00954 | 0.00959 | 0.03 |
| V-CPL-01 | √t front-growth coefficient | 0.000169109 | 0.000169117 | 0.000167516 | 0.00951 | 0.00956 | 0.03 |
| V-CPL-01 | back-face disturbance (semi-infinite validity) | 9.23744e-5 | 9.23867e-5 | 0 | 9.24e-5 | 9.24e-5 | 0.01 |
| V-CPL-01b | T-history L∞ error at x = 0.5 mm | 0.169686 | 0.169689 |  | 0.17 | 0.17 | 0.55 |
| V-CPL-01b | T-history L∞ error at x = 1.0 mm | 0.351478 | 0.351479 |  | 0.351 | 0.351 | 1.1 |
| V-CPL-01b | T-history L∞ error at x = 1.5 mm | 0.608528 | 0.608539 |  | 0.609 | 0.609 | 1.9 |
| V-CPL-01b | T-history L∞ error at x = 2.0 mm | 1.23142 | 1.23142 |  | 1.23 | 1.23 | 3.8 |
| V-CPL-01b | front position at t = 50.0 s | 0.00156779 | 0.00156779 | 0.00158957 | 0.0137 | 0.0137 | 0.045 |
| V-CPL-01b | front position at t = 100.0 s | 0.00222408 | 0.00222408 | 0.00224799 | 0.0106 | 0.0106 | 0.045 |
| V-CPL-01b | front position at t = 150.0 s | 0.00272619 | 0.00272619 | 0.00275322 | 0.00981 | 0.00982 | 0.045 |
| V-CPL-01b | √t front-growth coefficient | 0.000222188 | 0.000222188 | 0.000224799 | 0.0116 | 0.0116 | 0.035 |
| V-CPL-01b | back-face disturbance (semi-infinite validity) | 5.64121e-8 | 5.64137e-8 | 0 | 5.64e-8 | 5.64e-8 | 0.01 |
| V-CPL-02 | steady energy balance: v vs Landau v_ss(measured T_s) | 0.00194225 | 0.00194535 | 0.00194445 | 0.00113 | 0.00187 | 0.0021 |
| V-CPL-02 | regression consistency: thickness slope vs gas generation | 0.00194225 | 0.00194535 | 0.00194191 | 0.000179 | 0.000291 | 0.0021 |
| V-CPL-02 | mass ledger closure over plateau window | 0.777824 | 0.779089 | 0.776618 | 0.00155 | 0.00165 | 0.0049 |
| V-CPL-02 | in-depth profile decay rate (surface-attached frame) | 0.00196461 | 0.0019677 | 0.00194225 | 0.0115 | 0.0115 | 0.027 |
| V-CPL-02 | design-locus regression rate (loose) | 0.00194225 | 0.00194535 | 0.002 | 0.0289 | 0.0273 | 0.1 |
| V-CPL-02 | back-face rise (semi-infinite validity) | 0.000353589 | 0.000353623 | 0 | 0.000354 | 0.000354 | 0.01 |
| V-CPL-02b | steady energy balance: v vs Landau v_ss(measured T_s) | 0.00189854 | 0.00189991 | 0.00190147 | 0.00154 | 0.0018 | 0.002 |
| V-CPL-02b | regression consistency: thickness slope vs gas generation | 0.00189854 | 0.00189991 | 0.00189834 | 0.000108 | 0.000325 | 0.0017 |
| V-CPL-03 | P1: final thickness ratio | 2.0 | 2.0 | 2.0 | 7.12e-7 | 5.63e-8 | 3.0e-6 |
| V-CPL-03 | P1: total mass conservation | 2.0 | 2.0 | 2.0 | 6.55e-7 | 4.44e-16 | 2.0e-6 |
| V-CPL-03 | P2: final thickness ratio | 0.5 | 0.5 | 0.5 | 2.87e-8 | 1.13e-7 | 3.0e-6 |
| V-CPL-03 | P2: total mass conservation | 2.0 | 2.0 | 2.0 | 1.41e-7 | 8.88e-16 | 2.0e-6 |
| V-CPL-03 | P3: final thickness ratio | 1.0 | 1.0 | 1.0 | 0 | 0 | 3.0e-6 |
| V-CPL-03 | P3: total mass conservation | 2.0 | 2.0 | 2.0 | 0 | 3.33e-16 | 2.0e-6 |
| V-CPL-03 | P4: final thickness ratio | 0.3 | 0.3 | 0.3 | 1.01e-7 | 2.63e-7 | 3.0e-6 |
| V-CPL-03 | P4: total mass conservation | 2.0 | 2.0 | 2.0 | 1.61e-7 | 9.99e-16 | 2.0e-6 |
| V-CPL-03 | P4: final char mass (solid ledger) | 0.6 | 0.6 | 0.6 | 2.74e-7 | 1.13e-7 | 2.0e-6 |
| V-CPL-03 | P5: final thickness ratio | 1.0 | 1.0 | 1.0 | 0 | 0 | 3.0e-6 |
| V-CPL-03 | P5: total mass conservation | 2.0 | 2.0 | 2.0 | 6.66e-16 | 8.88e-16 | 2.0e-6 |
| V-CPL-03 | P5: final char mass (solid ledger) | 1.0 | 1.0 | 1.0 | 1.13e-7 | 1.13e-7 | 2.0e-6 |
| V-CPL-03 | P6: final thickness ratio | 0.5 | 0.5 | 0.5 | 2.87e-8 | 1.13e-7 | 3.0e-6 |
| V-CPL-03 | P6: total mass conservation | 2.0 | 2.0 | 2.0 | 1.41e-7 | 8.88e-16 | 2.0e-6 |
| V-CPL-03 | P7: final thickness ratio | 1.33333 | 1.33333 | 1.33333 | 7.12e-8 | 3.75e-8 | 3.0e-6 |
| V-CPL-03 | P7: total mass conservation | 2.0 | 2.0 | 2.0 | 2.05e-8 | 2.22e-16 | 2.0e-6 |
| V-CPL-03 | P7: final char mass (solid ledger) | 1.0 | 1.0 | 1.0 | 1.33e-7 | 1.13e-7 | 2.0e-6 |
| V-CPL-03 | P7: final lateral area ratio A_f/A₀ = law(χ̄ = 1/2) | 0.75 | 0.75 | 0.75 | 5.07e-8 | 3.75e-8 | 3.0e-6 |
| V-CPL-04 | steady center temperature vs full-Arrhenius shooting | 503.46 | 503.46 | 503.457 | 0.000738 | 0.000739 | 0.0025 |
| V-CPL-04 | steady center rise vs FK closed form (limit check) | 0.276761 | 0.276762 | 0.277242 | 0.00173 | 0.00173 | 0.006 |
| V-CPL-04 | steady profile L∞ vs shooting (relative to center rise) | 0.000737019 | 0.000738805 |  | 0.000737 | 0.000739 | 0.0025 |
| V-CPL-04 | steady-state reached (last-output drift) | 0.000380642 | 0.000385454 | 0 | 0.000381 | 0.000385 | 0.001 |
| V-CPL-04 | near-critical δ = 0.80 δ_c stays on the stable branch | 0.552044 | 0.552044 | 0.557537 | 0.00985 | 0.00985 | 0.03 |
| V-CPL-04 | supercritical δ = 1.50 δ_c ignites (thermal runaway) | 104.179 | 104.179 |  | 0 | 0 | 0.5 |
| V-CPL-05 | thin-limit deviation from 0D TGA (L = 0.5 mm) | 0.00141962 | 0.00141962 |  | 0.00142 | 0.00142 | 0.003 |
| V-CPL-05 | L-convergence rate, rung 1 → 2 | 2.0045 | 2.0045 | 2.0 | 0.0045 | 0.0045 | 0.1 |
| V-CPL-05 | L-convergence rate, rung 2 → 3 | 2.00153 | 2.00153 | 2.0 | 0.00153 | 0.00153 | 0.1 |
| V-CPL-05 | complete conversion at t_end (coarsest rung) | -1.15839e-15 | -3.9452e-11 | 0 | 1.16e-15 | 3.95e-11 | 1.0e-6 |
| V-GT-01 | L∞ ξ-profile error / peak, t=2.0 s | 0.00124671 | 0.00124672 |  | 0.00125 | 0.00125 | 0.004 |
| V-GT-01 | L∞ ξ-profile error / peak, t=4.0 s | 0.000770828 | 0.000770832 |  | 0.000771 | 0.000771 | 0.004 |
| V-GT-01 | L∞ ξ-profile error / peak, t=8.0 s | 0.000432363 | 0.000432358 |  | 0.000432 | 0.000432 | 0.004 |
| V-GT-01 | gas-mass conservation (sealed) | 1.9464e-15 | 1.08134e-15 |  | 1.95e-15 | 1.08e-15 | 1.0e-13 |
| V-GT-01 | variance growth, t=2.0 s | 5.0e-6 | 5.0e-6 | 5.0e-6 | 1.52e-15 | 2.54e-15 | 2.0e-5 |
| V-GT-01 | variance growth, t=4.0 s | 9.0e-6 | 9.0e-6 | 9.0e-6 | 1.66e-10 | 1.66e-10 | 2.0e-5 |
| V-GT-01 | variance growth, t=8.0 s | 1.69999e-5 | 1.69999e-5 | 1.7e-5 | 4.99e-6 | 4.99e-6 | 2.0e-5 |
| V-GT-02 | P² profile error / overpressure² scale | 0.000708073 | 0.000708073 |  | 0.000708 | 0.000708 | 0.0025 |
| V-GT-02 | back-face pressure | 110669.0 | 110669.0 | 110676.0 | 0.000716 | 0.000716 | 0.0025 |
| V-GT-02 | exit mass flux vs ġL | 0.0499996 | 0.0499996 | 0.05 | 7.25e-6 | 7.25e-6 | 5.0e-5 |
| V-GT-02 | steadiness of exit flux | 2.55824e-8 | 1.19071e-13 |  | 2.56e-8 | 1.19e-13 | 1.0e-6 |
| V-HC-01 | L∞ T-profile error, t=30.0 s | 0.240386 | 0.240404 |  | 0.24 | 0.24 | 0.75 |
| V-HC-01 | L∞ T-profile error, t=60.0 s | 0.205884 | 0.205903 |  | 0.206 | 0.206 | 0.75 |
| V-HC-01 | L∞ T-profile error, t=90.0 s | 0.182133 | 0.182126 |  | 0.182 | 0.182 | 0.75 |
| V-HC-01 | L∞ T-profile error, t=180.0 s | 0.138575 | 0.13857 |  | 0.139 | 0.139 | 0.75 |
| V-HC-01 | surface temperature, t=30.0 s | 546.251 | 546.251 | 545.747 | 0.503 | 0.503 | 1.5 |
| V-HC-01 | surface temperature, t=60.0 s | 625.757 | 625.757 | 625.502 | 0.255 | 0.255 | 1.5 |
| V-HC-01 | surface temperature, t=90.0 s | 680.072 | 680.072 | 679.912 | 0.159 | 0.159 | 1.5 |
| V-HC-01 | surface temperature, t=180.0 s | 784.967 | 784.967 | 784.91 | 0.0571 | 0.0571 | 1.5 |
| V-HC-01 | back-face rise (semi-infinite validity) | 1.15108e-9 | 1.15904e-9 | 0 | 1.15e-9 | 1.16e-9 | 1.0e-6 |
| V-HC-02a | L∞ T-profile error, Fo=0.05 | 0.000228527 | 0.000230697 |  | 0.000229 | 0.000231 | 0.06 |
| V-HC-02a | surface temperature, Fo=0.05 | 304.949 | 304.949 | 304.948 | 0.000916 | 0.00091 | 0.015 |
| V-HC-02a | L∞ T-profile error, Fo=0.1 | 0.000187165 | 0.000188381 |  | 0.000187 | 0.000188 | 0.06 |
| V-HC-02a | surface temperature, Fo=0.1 | 306.942 | 306.942 | 306.941 | 0.000625 | 0.000623 | 0.015 |
| V-HC-02a | L∞ T-profile error, Fo=0.2 | 0.000217297 | 0.000217815 |  | 0.000217 | 0.000218 | 0.06 |
| V-HC-02a | surface temperature, Fo=0.2 | 309.716 | 309.716 | 309.716 | 0.000448 | 0.000448 | 0.015 |
| V-HC-02a | L∞ T-profile error, Fo=0.5 | 0.000110151 | 0.000110097 |  | 0.00011 | 0.00011 | 0.06 |
| V-HC-02a | surface temperature, Fo=0.5 | 315.645 | 315.645 | 315.644 | 0.000417 | 0.000417 | 0.015 |
| V-HC-02b | L∞ T-profile error, Fo=0.05 | 0.00239577 | 0.00239599 |  | 0.0024 | 0.0024 | 0.06 |
| V-HC-02b | surface temperature, Fo=0.05 | 341.929 | 341.929 | 341.925 | 0.00436 | 0.00436 | 0.015 |
| V-HC-02b | L∞ T-profile error, Fo=0.1 | 0.00198947 | 0.00198926 |  | 0.00199 | 0.00199 | 0.06 |
| V-HC-02b | surface temperature, Fo=0.1 | 355.287 | 355.287 | 355.285 | 0.00215 | 0.00215 | 0.015 |
| V-HC-02b | L∞ T-profile error, Fo=0.2 | 0.00141015 | 0.00141028 |  | 0.00141 | 0.00141 | 0.06 |
| V-HC-02b | surface temperature, Fo=0.2 | 371.323 | 371.323 | 371.322 | 0.00105 | 0.00105 | 0.015 |
| V-HC-02b | L∞ T-profile error, Fo=0.5 | 0.00076202 | 0.000761883 |  | 0.000762 | 0.000762 | 0.06 |
| V-HC-02b | surface temperature, Fo=0.5 | 399.096 | 399.096 | 399.096 | 0.000713 | 0.000713 | 0.015 |
| V-HC-02c | L∞ T-profile error, Fo=0.05 | 0.0180188 | 0.0180182 |  | 0.018 | 0.018 | 0.06 |
| V-HC-02c | surface temperature, Fo=0.05 | 453.53 | 453.53 | 453.535 | 0.00487 | 0.00487 | 0.015 |
| V-HC-02c | L∞ T-profile error, Fo=0.1 | 0.00955771 | 0.00955729 |  | 0.00956 | 0.00956 | 0.06 |
| V-HC-02c | surface temperature, Fo=0.1 | 465.883 | 465.883 | 465.885 | 0.00235 | 0.00235 | 0.015 |
| V-HC-02c | L∞ T-profile error, Fo=0.2 | 0.00434463 | 0.00434493 |  | 0.00434 | 0.00434 | 0.06 |
| V-HC-02c | surface temperature, Fo=0.2 | 475.503 | 475.503 | 475.504 | 0.000817 | 0.000817 | 0.015 |
| V-HC-02c | L∞ T-profile error, Fo=0.5 | 0.00362978 | 0.00362964 |  | 0.00363 | 0.00363 | 0.06 |
| V-HC-02c | surface temperature, Fo=0.5 | 487.134 | 487.134 | 487.134 | 0.000316 | 0.000316 | 0.015 |
| V-HC-02d | L∞ T-profile error, Fo=0.05 | 0.0204031 | 0.0204033 |  | 0.0204 | 0.0204 | 0.06 |
| V-HC-02d | surface temperature, Fo=0.05 | 494.958 | 494.958 | 494.959 | 0.000978 | 0.000978 | 0.015 |
| V-HC-02d | L∞ T-profile error, Fo=0.1 | 0.0100111 | 0.0100111 |  | 0.01 | 0.01 | 0.06 |
| V-HC-02d | surface temperature, Fo=0.1 | 496.433 | 496.433 | 496.434 | 0.000344 | 0.000344 | 0.015 |
| V-HC-02d | L∞ T-profile error, Fo=0.2 | 0.00449762 | 0.00449742 |  | 0.0045 | 0.0045 | 0.06 |
| V-HC-02d | surface temperature, Fo=0.2 | 497.508 | 497.508 | 497.508 | 9.93e-5 | 9.93e-5 | 0.015 |
| V-HC-02d | L∞ T-profile error, Fo=0.5 | 0.00409437 | 0.00409411 |  | 0.00409 | 0.00409 | 0.06 |
| V-HC-02d | surface temperature, Fo=0.5 | 498.818 | 498.818 | 498.819 | 4.5e-5 | 4.5e-5 | 0.015 |
| V-HC-03 | L∞ T-profile error, Fo=0.1 | 0.00196561 | 0.00196616 |  | 0.00197 | 0.00197 | 0.006 |
| V-HC-03 | L∞ T-profile error, Fo=0.5 | 0.00146635 | 0.0014682 |  | 0.00147 | 0.00147 | 0.006 |
| V-HC-03 | L∞ T-profile error, Fo=1.5 | 0.00130224 | 0.00130211 |  | 0.0013 | 0.0013 | 0.006 |
| V-HC-03 | front–back ΔT plateau, Fo=1.5 | 98.75 | 98.75 | 98.75 | 3.22e-7 | 5.62e-8 | 1.0e-5 |
| V-HC-03 | mean-temperature rise (global energy) | 200.0 | 200.0 | 200.0 | 1.14e-13 | 1.53e-12 | 1.0e-5 |
| V-HC-04a | L∞ steady profile vs (quantized) network | 1.54046e-11 | 5.68434e-14 |  | 1.54e-11 | 5.68e-14 | 1.0e-8 |
| V-HC-04a | steady bottom-face flux vs network | 909.091 | 909.091 | 909.091 | 3.39e-11 | 9.08e-14 | 1.0e-8 |
| V-HC-04a | steady top-face flux vs network | 909.091 | 909.091 | 909.091 | 3.34e-13 | 9.25e-15 | 1.0e-8 |
| V-HC-04a | steadiness (last half-interval drift) | 8.46001e-10 | 1.13687e-13 |  | 8.46e-10 | 1.14e-13 | 1.0e-6 |
| V-HC-04b | L∞ steady profile vs (quantized) network | 3.02407e-11 | 5.68434e-14 |  | 3.02e-11 | 5.68e-14 | 1.0e-8 |
| V-HC-04b | steady bottom-face flux vs network | 99.0099 | 99.0099 | 99.0099 | 4.08e-10 | 5.91e-13 | 1.0e-8 |
| V-HC-04b | steady top-face flux vs network | 99.0099 | 99.0099 | 99.0099 | 9.77e-14 | 5.88e-15 | 1.0e-8 |
| V-HC-04b | steadiness (last half-interval drift) | 3.11189e-9 | 1.13687e-13 |  | 3.11e-9 | 1.14e-13 | 1.0e-6 |
| V-HC-04c | L∞ steady profile vs (quantized) network | 1.54046e-11 | 5.68434e-14 |  | 1.54e-11 | 5.68e-14 | 1.0e-8 |
| V-HC-04c | steady bottom-face flux vs network | 909.091 | 909.091 | 909.091 | 3.39e-11 | 9.08e-14 | 1.0e-8 |
| V-HC-04c | steady top-face flux vs network | 909.091 | 909.091 | 909.091 | 3.34e-13 | 9.25e-15 | 1.0e-8 |
| V-HC-04c | interface-quantization smearing (documented) | 0.0163636 | 0.0163636 | 0 | 0.0164 | 0.0164 | 0.0208 |
| V-HC-04c | steadiness (last half-interval drift) | 8.46001e-10 | 1.13687e-13 |  | 8.46e-10 | 1.14e-13 | 1.0e-6 |
| V-HC-05 | L∞ Kirchhoff-profile error (z(T) map, /L) | 3.86964e-5 | 3.86964e-5 |  | 3.87e-5 | 3.87e-5 | 0.00012 |
| V-HC-05 | steady flux vs Kirchhoff-exact | 32500.6 | 32500.6 | 32500.0 | 1.75e-5 | 1.75e-5 | 6.0e-5 |
| V-HC-05 | steadiness (last half-interval drift) | 3.5888e-8 | 2.27374e-13 |  | 3.59e-8 | 2.27e-13 | 1.0e-6 |
| V-MMS-01 | L∞ T error, t=10.0 s | 0.0165159 | 0.0186741 |  | 0.0165 | 0.0187 | 0.05 |
| V-MMS-01 | L∞ T error, t=20.0 s | 0.0122428 | 0.013515 |  | 0.0122 | 0.0135 | 0.05 |
| V-MMS-01 | relative L∞ ξ_A error, t=20.0 s | 6.80988e-5 | 7.50646e-5 |  | 6.81e-5 | 7.51e-5 | 0.0002 |
| V-MMS-01 | relative L∞ ξ_B error, t=20.0 s | 2.18523e-5 | 2.40685e-5 |  | 2.19e-5 | 2.41e-5 | 0.0001 |
| V-MMS-01 | relative L∞ ξ_gas error, t=20.0 s | 0.00990746 | 0.0110621 |  | 0.00991 | 0.0111 | 0.03 |
| V-MMS-02 | free-stream T drift (all cells, all snapshots) | 5.68434e-14 | 3.41061e-13 | 0 | 5.68e-14 | 3.41e-13 | 1.0e-10 |
| V-MMS-02 | free-stream ξ_solid relative drift | 1.89478e-16 | 1.83794e-14 | 0 | 1.89e-16 | 1.84e-14 | 1.0e-12 |
| V-MMS-02 | free-stream ξ_gas relative drift | 1.4803e-16 | 2.91619e-14 | 0 | 1.48e-16 | 2.92e-14 | 1.0e-12 |
| V-MMS-02 | top cell-center trajectory vs z₀·s(t_end) | 0.006825 | 0.006825 | 0.006825 | 4.04e-7 | 7.32e-9 | 1.5e-6 |
| V-MMS-03 | L∞ T error on moving mesh, t=10.0 s | 0.0710876 | 0.0090851 |  | 0.0711 | 0.00909 | 0.25 |
| V-MMS-03 | relative L∞ ξ_solid (mesh-bound) error, t=10.0 s | 2.10811e-11 | 5.31656e-12 |  | 2.11e-11 | 5.32e-12 | 1.0e-10 |
| V-MMS-03 | relative L∞ ξ_gas error, t=10.0 s | 0.00178289 | 0.00168 |  | 0.00178 | 0.00168 | 0.006 |
| V-MMS-04 | L∞ T error (DirectSolve), t=20.0 s | 2.03406e-5 | 2.06606e-5 |  | 2.03e-5 | 2.07e-5 | 6.0e-5 |
| V-MMS-04 | relative L∞ ξ_solid (z-varying, material frame) error, t=20.0 s | 6.57437e-7 | 6.59427e-7 |  | 6.57e-7 | 6.59e-7 | 2.0e-6 |
| V-MMS-04 | relative L∞ ξ_gas error, t=20.0 s | 0.00229895 | 0.00229964 |  | 0.0023 | 0.0023 | 0.007 |
| V-MMS-04 | relative surface-position error, t=20.0 s | 1.63289e-5 | 1.63342e-5 | 0.00896269 | 1.63e-5 | 1.63e-5 | 5.0e-5 |
| V-MMS-04 | L∞ T error (ApproxSolve), t=20.0 s | 2.04286e-5 |  |  | 2.04e-5 |  | 6.0e-5 |
| V-MMS-04 | relative L∞ ξ_solid error (ApproxSolve), t=20.0 s | 6.57287e-7 |  |  | 6.57e-7 |  | 2.0e-6 |
| V-RAD-01 | L∞ T-profile error, t=30.0 s | 0.0438306 | 0.0434579 |  | 0.0438 | 0.0435 | 0.25 |
| V-RAD-01 | L∞ T-profile error, t=60.0 s | 0.056841 | 0.0566988 |  | 0.0568 | 0.0567 | 0.25 |
| V-RAD-01 | L∞ T-profile error, t=90.0 s | 0.0649277 | 0.0648546 |  | 0.0649 | 0.0649 | 0.25 |
| V-RAD-01 | L∞ T-profile error, t=180.0 s | 0.0796417 | 0.0796245 |  | 0.0796 | 0.0796 | 0.25 |
| V-RAD-01 | front-face temperature, t=30.0 s | 468.542 | 468.541 | 468.622 | 0.0803 | 0.0806 | 0.25 |
| V-RAD-01 | front-face temperature, t=60.0 s | 548.123 | 548.123 | 548.201 | 0.0775 | 0.0776 | 0.25 |
| V-RAD-01 | front-face temperature, t=90.0 s | 604.294 | 604.294 | 604.368 | 0.074 | 0.0741 | 0.25 |
| V-RAD-01 | front-face temperature, t=180.0 s | 717.608 | 717.608 | 717.673 | 0.0648 | 0.0648 | 0.25 |
| V-RAD-01 | sub-surface peak value, t=90 s | 608.944 | 608.944 | 608.932 | 0.0124 | 0.0123 | 0.05 |
| V-RAD-01 | sub-surface peak location, t=90 s | 0.00035 | 0.00035 | 0.000319 | 3.1e-5 | 3.1e-5 | 0.00015 |
| V-RAD-01 | peak strictly interior, t=90 s | 97.0 | 97.0 |  | 0 | 0 | 0.5 |
| V-RX-01 | max relative ξ_A error | 1.25588e-9 | 4.77257e-13 |  | 1.26e-9 | 4.77e-13 | 1.0e-8 |
| V-RX-01 | max relative ξ_B error | 2.78142e-9 | 1.20696e-12 |  | 2.78e-9 | 1.21e-12 | 1.0e-8 |
| V-RX-01 | max relative ξ_C error | 1.52554e-9 | 7.29472e-13 |  | 1.53e-9 | 7.29e-13 | 1.0e-8 |
| V-RX-01 | total-mass conservation (closed chain) | 2.27374e-15 | 1.25056e-15 |  | 2.27e-15 | 1.25e-15 | 1.0e-13 |
| V-RX-01 | isothermality drift | 5.68434e-14 | 0 |  | 5.68e-14 | 0 | 1.0e-12 |
| V-RX-02a | max relative ξ_A error (n=0.0) | 1.42109e-16 | 7.4607e-15 |  | 1.42e-16 | 7.46e-15 | 1.0e-12 |
| V-RX-02a | total-mass conservation | 2.84217e-16 | 7.10543e-16 |  | 2.84e-16 | 7.11e-16 | 1.0e-13 |
| V-RX-02b | max relative ξ_A error (n=2.0) | 1.72038e-8 | 2.92971e-12 |  | 1.72e-8 | 2.93e-12 | 5.0e-8 |
| V-RX-02b | total-mass conservation | 1.13687e-15 | 5.68434e-16 |  | 1.14e-15 | 5.68e-16 | 1.0e-13 |
| V-RX-02c | max relative ξ_A error (n=0.6666666666666666) | 2.77351e-10 | 2.26259e-10 |  | 2.77e-10 | 2.26e-10 | 1.0e-9 |
| V-RX-02c | total-mass conservation | 4.54747e-16 | 6.82121e-16 |  | 4.55e-16 | 6.82e-16 | 1.0e-13 |
| V-RX-02d | non-negativity through 2× burnout | -8.46971e-13 | 0.00539209 | 0 | 8.47e-13 | 0 | 1.0e-9 |
| V-RX-02d | residual mass within limiter threshold | -8.46971e-13 | 0.00539209 | 0 | 8.47e-13 | 0.00539 | 1.0 |
| V-RX-03 | max relative ξ_A error | 2.26927e-9 | 1.7705e-12 |  | 2.27e-9 | 1.77e-12 | 1.0e-8 |
| V-RX-03 | max relative ξ_B error | 1.05899e-9 | 8.26674e-13 |  | 1.06e-9 | 8.27e-13 | 1.0e-8 |
| V-RX-03 | max relative ξ_C error | 3.02569e-10 | 2.36142e-13 |  | 3.03e-10 | 2.36e-13 | 1.0e-8 |
| V-RX-03 | max relative ξ_G error | 9.07707e-10 | 7.08553e-13 |  | 9.08e-10 | 7.09e-13 | 1.0e-8 |
| V-RX-03 | total-mass conservation (incl. retained gas) | 1.36424e-15 | 1.02318e-15 |  | 1.36e-15 | 1.02e-15 | 1.0e-13 |
| V-RX-03 | reaction-heat superposition (quasi-steady cell drift) | 3.23527e-5 | 3.23527e-5 |  | 3.24e-5 | 3.24e-5 | 0.0001 |
| V-RX-04a | max mean-T error vs implicit exact | 0.00389327 | 0.0041815 |  | 0.00389 | 0.00418 | 0.012 |
| V-RX-04a | internal ΔT (thin-limit validity) | 0.00685372 | 0.00685374 |  | 0.00685 | 0.00685 | 0.05 |
| V-RX-04b | steady surface temperature vs algebraic root | 732.797 | 732.797 | 732.797 | 3.75e-9 | 2.54e-10 | 1.0e-6 |
| V-RX-04b | steadiness (last half-interval drift) | 0.000107416 | 0.000107343 |  | 0.000107 | 0.000107 | 0.001 |
| V-RX-04b | steady uniformity (adiabatic back) | 8.7698e-10 | 2.467e-11 |  | 8.77e-10 | 2.47e-11 | 1.0e-7 |
