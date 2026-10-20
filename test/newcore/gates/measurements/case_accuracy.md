# Gate items 6 and 8 (Phase 1) and 7 (Phase 2): constraint residuals, achieved error and event timing

Produced by `test/newcore/gates/case_accuracy.jl` at commit 7ceaf74. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [18.2, 17.1, 17.8] at the start and [16.7, 17.2, 17.5] at the end. The cases V-BO-01, V-BO-02, V-BO-04, V-CPL-02, V-CPL-02b, each at its CI resolution and solve keywords, on the default integrator. Error: the largest temperature difference [K] over the saved times and the cells live in both solves, against the same case at both tolerances divided by 100, or by 10 where that does not complete (tightened) (under deactivation the saved times with the same live cells in both; a case saved at every step its final state; when one solve stops before the other, the saved times both reached; compared: the saved times compared, zero for the final state; to: the last time compared [s]); ratio: that error over reltol times the reference's largest temperature, within budget at 10 or less. Under deactivation each case is solved with event location off and on (location; its reference likewise) and judged on item 7 below: invocations, the deactivation invocations; event shift, the largest difference of the invocation times from the reference's over the median interval between the reference's invocations (between its saved times for a single invocation). Constraint: the largest scaled constraint residual of the accepted states (one or below passes; zero for a fixed-step solve or a configuration without algebraic rows). Wall: one warm solve [s].

| case | form | location | retcode | reltol | tightened | error [K] | compared | to [s] | ratio | invocations | event shift | within budget | constraint | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | off | Success | 1.0e-6 | 100 | 2.49 | 1001 | 500.0 | 2580.0 | 1 | 0.218 | true | 0.609 | 1126 | 4.64 |
| V-BO-01 | temperature | on | Success | 1.0e-6 | 100 | 0.234 | 1001 | 500.0 | 243.0 | 1 | 0.00601 | true | 0.609 | 1127 | 4.58 |
| V-BO-02 | temperature | off | Success | 1.0e-6 | 100 | 25.1 | 938 | 500.0 | 27500.0 | 60 | 0.138 | true | 0.0 | 2014 | 3.84 |
| V-BO-02 | temperature | on | Success | 1.0e-6 | 100 | 0.137 | 1000 | 500.0 | 150.0 | 60 | 0.000345 | true | 0.0 | 2035 | 3.98 |
| V-BO-04 | temperature | off | Success | 1.0e-7 | 100 | 1.25e-12 | 241 | 480.0 | 2.08e-8 | 1 | 0.133 | true | 5.03e-8 | 324 | 0.0754 |
| V-BO-04 | temperature | on | Success | 1.0e-7 | 100 | 1.25e-12 | 241 | 480.0 | 2.08e-8 | 1 | 1.21e-5 | true | 5.03e-8 | 324 | 0.0821 |
| V-CPL-02 | temperature | off | Success | 1.0e-5 | 100 | 0.0427 | 41 | 2.0 | 5.78 | 58 | 0.0754 | true | 0.256 | 2049 | 1.76 |
| V-CPL-02 | temperature | on | Success | 1.0e-5 | 100 | 0.00919 | 41 | 2.0 | 1.24 | 58 | 0.000189 | true | 0.372 | 2094 | 1.87 |
| V-CPL-02b | temperature | off | Success | 1.0e-5 | 100 | 0.128 | 41 | 2.0 | 16.5 | 34 | 0.0772 | true | 0.0599 | 1192 | 0.619 |
| V-CPL-02b | temperature | on | Success | 1.0e-5 | 100 | 0.00766 | 41 | 2.0 | 0.992 | 34 | 8.63e-5 | true | 0.0577 | 1221 | 0.667 |
| V-BO-01 | enthalpy | off | Success | 1.0e-6 | 100 | 3.85 | 1000 | 500.0 | 3990.0 | 1 | 0.407 | true | 0.382 | 1571 | 15.0 |
| V-BO-01 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0308 | 1001 | 500.0 | 31.9 | 1 | 5.84e-6 | true | 0.382 | 1571 | 19.3 |
| V-BO-02 | enthalpy | off | Success | 1.0e-6 | 100 | 22.6 | 967 | 500.0 | 24800.0 | 60 | 0.0862 | true | 0.0 | 2794 | 20.6 |
| V-BO-02 | enthalpy | on | Success | 1.0e-6 | 100 | 0.00907 | 1001 | 500.0 | 9.94 | 60 | 2.76e-5 | true | 0.0 | 2846 | 14.7 |
| V-BO-04 | enthalpy | off | Success | 1.0e-7 | 100 | 3.99e-10 | 240 | 480.0 | 6.65e-6 | 1 | 0.159 | true | 2.31e-8 | 327 | 0.199 |
| V-BO-04 | enthalpy | on | Success | 1.0e-7 | 100 | 3.99e-10 | 241 | 480.0 | 6.65e-6 | 1 | 5.24e-6 | true | 2.31e-8 | 327 | 0.2 |
| V-CPL-02 | enthalpy | off | Success | 1.0e-5 | 100 | 0.126 | 40 | 2.0 | 17.0 | 58 | 0.0689 | true | 0.375 | 2052 | 6.51 |
| V-CPL-02 | enthalpy | on | Success | 1.0e-5 | 100 | 0.03 | 41 | 2.0 | 4.06 | 58 | 7.27e-5 | true | 0.51 | 2101 | 4.98 |
| V-CPL-02b | enthalpy | off | Success | 1.0e-5 | 100 | 0.0851 | 41 | 2.0 | 11.0 | 34 | 0.0626 | true | 0.0791 | 1261 | 2.03 |
| V-CPL-02b | enthalpy | on | Success | 1.0e-5 | 100 | 0.018 | 41 | 2.0 | 2.33 | 34 | 2.29e-5 | true | 0.0791 | 1291 | 2.46 |

20 of 20 solves within budget; 20 of 20 with every constraint residual at or below one.

## Item 7: event timing with location off and on

Per deactivating case and form, the largest shift of the invocation times between the solve and its reference [s] with location off and on, the median interval between the reference's invocations, or between its saved times for a single invocation [s], and the shifts as fractions of it, and gain, the shift off over the shift on. Passes when the located shift is at most 0.05 of the interval and the gain at least 3.0, both solves completing as the case expects.

| case | form | invocations | shift off [s] | shift on [s] | interval [s] | off / interval | on / interval | gain | passes |
|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | 1 | 0.109 | 0.00301 | 0.5 | 0.218 | 0.00601 | 36.2 | true |
| V-BO-02 | temperature | 60 | 0.925 | 0.0023 | 6.68 | 0.138 | 0.000345 | 402.0 | true |
| V-BO-04 | temperature | 1 | 0.266 | 2.43e-5 | 2.0 | 0.133 | 1.21e-5 | 10900.0 | true |
| V-CPL-02 | temperature | 58 | 0.00243 | 6.1e-6 | 0.0323 | 0.0754 | 0.000189 | 399.0 | true |
| V-CPL-02b | temperature | 34 | 0.00409 | 4.56e-6 | 0.0528 | 0.0775 | 8.63e-5 | 898.0 | true |
| V-BO-01 | enthalpy | 1 | 0.204 | 2.92e-6 | 0.5 | 0.407 | 5.84e-6 | 69800.0 | true |
| V-BO-02 | enthalpy | 60 | 0.574 | 0.000185 | 6.68 | 0.0859 | 2.76e-5 | 3110.0 | true |
| V-BO-04 | enthalpy | 1 | 0.318 | 1.05e-5 | 2.0 | 0.159 | 5.24e-6 | 30300.0 | true |
| V-CPL-02 | enthalpy | 58 | 0.00225 | 2.35e-6 | 0.0323 | 0.0696 | 7.27e-5 | 958.0 | true |
| V-CPL-02b | enthalpy | 34 | 0.00333 | 1.21e-6 | 0.0528 | 0.063 | 2.29e-5 | 2750.0 | true |

10 of 10 pass.
