# Gate items 6 and 8 (Phase 1) and 7 (Phase 2): constraint residuals, achieved error and event timing

Produced by `test/newcore/gates/case_accuracy.jl` at commit f54a6fb (tree b544f02). Julia 1.13.0, 1 thread, 22 logical CPUs; load average [20.2, 18.2, 18.3] at the start and [19.7, 19.9, 19.5] at the end. Every case of `NEW_CORE_CASES`, each at its CI resolution and solve keywords, on the default integrator. Error: the largest temperature difference [K] over the saved times and the cells live in both solves, against the same case at both tolerances divided by 100, or by 10 where that does not complete (tightened) (under deactivation the saved times with the same live cells in both; a case saved at every step its final state; when one solve stops before the other, the saved times both reached; compared: the saved times compared, zero for the final state; to: the last time compared [s]); ratio: that error over reltol times the reference's largest temperature, within budget at 10 or less. Under deactivation each case is solved with event location off and on (location; its reference likewise) and judged on item 7 below: invocations, the deactivation invocations; event shift, the largest difference of the invocation times from the reference's over the median interval between the reference's invocations (between its saved times for a single invocation). Constraint: the largest scaled constraint residual of the accepted states (one or below passes; zero for a fixed-step solve or a configuration without algebraic rows). Wall: one warm solve [s].

| case | form | location | retcode | reltol | tightened | error [K] | compared | to [s] | ratio | invocations | event shift | within budget | constraint | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | off | Success | 1.0e-6 | 100 | 12.1 | 1000 | 500.0 | 12500.0 | 1 | 0.822 | true | 0.559 | 1124 | 3.49 |
| V-BO-01 | temperature | on | Success | 1.0e-6 | 100 | 0.248 | 1001 | 500.0 | 257.0 | 1 | 0.00601 | true | 0.559 | 1125 | 5.59 |
| V-BO-02 | temperature | off | Success | 1.0e-6 | 100 | 26.1 | 940 | 500.0 | 28600.0 | 60 | 0.12 | true | 0.0 | 1944 | 2.21 |
| V-BO-02 | temperature | on | Success | 1.0e-6 | 100 | 0.179 | 1001 | 500.0 | 197.0 | 60 | 0.000316 | true | 0.0 | 1977 | 2.68 |
| V-BO-03 | temperature | off | Success | 1.0e-6 | 100 | 0.0912 | 1001 | 500.0 | 133.0 | 0 |  | false | 0.539 | 1073 | 3.43 |
| V-BO-03 | temperature | on | Success | 1.0e-6 | 100 | 0.0912 | 1001 | 500.0 | 133.0 | 0 |  | false | 0.539 | 1073 | 2.51 |
| V-BO-04 | temperature | off | Success | 1.0e-7 | 100 | 1.25e-12 | 240 | 480.0 | 2.08e-8 | 1 | 1.31 | true | 5.38e-10 | 321 | 0.143 |
| V-BO-04 | temperature | on | Success | 1.0e-7 | 100 | 1.25e-12 | 241 | 480.0 | 2.08e-8 | 1 | 2.7e-8 | true | 5.38e-10 | 322 | 0.228 |
| V-BO-05 | temperature | off | Terminated | 1.0e-6 | 100 | 0.138 | 520 | 260.0 | 199.0 | 0 |  | false | 0.0184 | 347 | 0.542 |
| V-BO-05 | temperature | on | Terminated | 1.0e-6 | 100 | 0.138 | 520 | 260.0 | 199.0 | 0 |  | false | 0.0184 | 347 | 0.592 |
| V-BO-06 | temperature | off | Success | 1.0e-6 | 100 | 5.51 | 1001 | 500.0 | 5580.0 | 1 | 0.213 | true | 0.873 | 1133 | 2.38 |
| V-BO-06 | temperature | on | Success | 1.0e-6 | 100 | 0.268 | 1001 | 500.0 | 271.0 | 1 | 0.00533 | true | 0.873 | 1132 | 2.39 |
| V-BO-07 | temperature | off | Success | 1.0e-6 | 100 | 5.35 | 1001 | 500.0 | 5540.0 | 1 | 0.231 | true | 0.668 | 724 | 2.72 |
| V-BO-07 | temperature | on | Success | 1.0e-6 | 100 | 0.177 | 1001 | 500.0 | 184.0 | 1 | 0.00333 | true | 0.668 | 724 | 4.73 |
| V-BO-08 | temperature | off | Success | 1.0e-6 | 100 | 1.97 | 301 | 150.0 | 2060.0 | 1 | 0.222 | true | 0.0498 | 213 | 0.559 |
| V-BO-08 | temperature | on | Success | 1.0e-6 | 100 | 0.00708 | 301 | 150.0 | 7.41 | 1 | 6.18e-5 | true | 0.0498 | 213 | 1.01 |
| V-CPL-01 | temperature |  | Success | 1.0e-6 | 100 | 0.000421 | 31 | 300.0 | 0.96 | 0 |  | true | 0.0 | 246 | 0.171 |
| V-CPL-01b | temperature |  | Success | 1.0e-6 | 100 | 2.83e-5 | 31 | 150.0 | 0.0586 | 0 |  | true | 0.0 | 1141 | 1.31 |
| V-CPL-02 | temperature | off | Success | 1.0e-5 | 100 | 0.629 | 40 | 2.0 | 85.5 | 19 | 0.0307 | true | 0.856 | 1123 | 2.59 |
| V-CPL-02 | temperature | on | Success | 1.0e-5 | 100 | 0.0284 | 41 | 2.0 | 3.86 | 19 | 0.000184 | true | 0.999 | 1134 | 2.59 |
| V-CPL-02b | temperature | off | Success | 1.0e-5 | 100 | 1.29 | 40 | 2.0 | 167.0 | 11 | 0.033 | true | 0.189 | 636 | 0.983 |
| V-CPL-02b | temperature | on | Success | 1.0e-5 | 100 | 0.0161 | 41 | 2.0 | 2.09 | 11 | 5.13e-5 | true | 0.197 | 638 | 1.02 |
| V-CPL-03 | temperature |  | Success | 1.0e-8 | 100 | 0.0 | 3 | 8.0 | 0.0 | 0 |  | true | 1.05e-9 | 88 | 0.00884 |
| V-CPL-04 | temperature |  | Success | 1.0e-6 | 100 | 1.13e-6 | 11 | 2500.0 | 0.00224 | 0 |  | true | 0.0 | 28 | 0.0063 |
| V-CPL-05 | temperature |  | Success | 1.0e-8 | 10 | 2.38e-6 | 72 | 720.0 | 0.361 | 0 |  | true | 0.0 | 237 | 0.0294 |
| V-GT-01 | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 4 | 8.0 | 0.0 | 0 |  | true | 0.0 | 46 | 0.0381 |
| V-GT-02 | temperature |  | Success | 1.0e-8 | 100 | 0.0 | 3 | 1.0 | 0.0 | 0 |  | true | 0.0 | 52 | 0.0324 |
| V-HC-01 | temperature |  | Success | 1.0e-8 | 100 | 3.03e-7 | 5 | 180.0 | 0.0395 | 0 |  | true | 0.0 | 60 | 0.0205 |
| V-HC-02a | temperature |  | Success | 1.0e-8 | 100 | 4.51e-8 | 5 | 10000.0 | 0.0143 | 0 |  | true | 0.0 | 56 | 0.0144 |
| V-HC-02b | temperature |  | Success | 1.0e-8 | 100 | 9.55e-8 | 5 | 10000.0 | 0.024 | 0 |  | true | 0.0 | 75 | 0.019 |
| V-HC-02c | temperature |  | Success | 1.0e-8 | 100 | 1.24e-7 | 5 | 10000.0 | 0.0256 | 0 |  | true | 0.0 | 98 | 0.0228 |
| V-HC-02d | temperature |  | Success | 1.0e-8 | 100 | 1.28e-7 | 5 | 10000.0 | 0.0258 | 0 |  | true | 0.0 | 113 | 0.0265 |
| V-HC-03 | temperature |  | Success | 1.0e-8 | 100 | 2.43e-8 | 4 | 1200.0 | 0.00366 | 0 |  | true | 0.0 | 86 | 0.0201 |
| V-HC-04a | temperature |  | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 | 0 |  | true | 0.0 | 99 | 0.0158 |
| V-HC-04b | temperature |  | Success | 1.0e-8 | 100 | 3.98e-13 | 3 | 400000.0 | 1.0e-7 | 0 |  | true | 0.0 | 97 | 0.0153 |
| V-HC-04c | temperature |  | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 | 0 |  | true | 0.0 | 99 | 0.0153 |
| V-HC-05 | temperature |  | Success | 1.0e-8 | 100 | 1.14e-13 | 3 | 5000.0 | 1.43e-8 | 0 |  | true | 0.0 | 243 | 0.0446 |
| V-MMS-01 | temperature |  | Success | 1.0e-8 | 100 | 2.2e-8 | 3 | 20.0 | 0.00326 | 0 |  | true | 0.0 | 127 | 0.536 |
| V-MMS-02 | temperature |  | Success | 1.0e-8 | 100 | 3.98e-13 | 41 | 10.0 | 7.96e-8 | 0 |  | true | 0.00178 | 26 | 0.0152 |
| V-MMS-03 | temperature |  | Success | 1.0e-8 | 100 | 1.56e-5 | 41 | 10.0 | 2.52 | 0 |  | true | 2.05e-6 | 77 | 0.617 |
| V-MMS-04 | temperature |  | Success | 1.0e-8 | 100 | 8.98e-9 | 41 | 20.0 | 0.00146 | 0 |  | true | 9.02e-9 | 63 | 0.527 |
| V-RAD-01 | temperature |  | Success | 1.0e-8 | 100 | 2.08e-7 | 5 | 180.0 | 0.0286 | 0 |  | true | 0.0 | 63 | 0.0197 |
| V-RX-01 | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 7 | 20.0 | 0.0 | 0 |  | true | 0.0 | 175 | 0.00278 |
| V-RX-02a | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 | 0 |  | true | 0.0 | 13 | 0.000989 |
| V-RX-02b | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 | 0 |  | true | 0.0 | 242 | 0.00322 |
| V-RX-02c | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 7 | 7.17 | 0.0 | 0 |  | true | 0.0 | 183 | 0.00292 |
| V-RX-02d | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 21 | 40.0 | 1.42e-6 | 0 |  | true | 0.0 | 306 | 0.00418 |
| V-RX-03 | temperature |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 15.0 | 2.84e-6 | 0 |  | true | 0.0 | 160 | 0.00268 |
| V-RX-04a | temperature |  | Success | 1.0e-8 | 100 | 3.02e-7 | 6 | 8.0 | 0.0336 | 0 |  | true | 0.0 | 78 | 0.00238 |
| V-RX-04b | temperature |  | Success | 1.0e-8 | 100 | 5.97e-9 | 3 | 2000.0 | 0.000814 | 0 |  | true | 0.0 | 159 | 0.00767 |
| V-BO-01 | enthalpy | off | Success | 1.0e-6 | 100 | 6.47 | 1001 | 500.0 | 6720.0 | 1 | 0.55 | true | 0.475 | 1585 | 20.1 |
| V-BO-01 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0366 | 1001 | 500.0 | 38.0 | 1 | 3.36e-6 | true | 0.475 | 1584 | 21.0 |
| V-BO-02 | enthalpy | off | Success | 1.0e-6 | 100 | 18.6 | 947 | 500.0 | 20400.0 | 60 | 0.0942 | true | 0.0 | 2771 | 14.5 |
| V-BO-02 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0147 | 1001 | 500.0 | 16.2 | 60 | 1.03e-5 | true | 0.0 | 2785 | 14.2 |
| V-BO-03 | enthalpy | off | Success | 1.0e-6 | 100 | 0.0323 | 1001 | 500.0 | 47.2 | 0 |  | false | 0.441 | 1542 | 27.8 |
| V-BO-03 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0323 | 1001 | 500.0 | 47.2 | 0 |  | false | 0.441 | 1542 | 17.3 |
| V-BO-04 | enthalpy | off | Success | 1.0e-7 | 100 | 3.19e-7 | 241 | 480.0 | 0.00532 | 1 | 0.176 | true | 1.65e-8 | 322 | 0.441 |
| V-BO-04 | enthalpy | on | Success | 1.0e-7 | 100 | 3.19e-7 | 241 | 480.0 | 0.00532 | 1 | 2.11e-7 | true | 1.65e-8 | 322 | 0.268 |
| V-BO-05 | enthalpy | off | Terminated | 1.0e-6 | 100 | 0.0146 | 520 | 260.0 | 21.1 | 0 |  | false | 0.00359 | 837 | 12.7 |
| V-BO-05 | enthalpy | on | Terminated | 1.0e-6 | 100 | 0.0146 | 520 | 260.0 | 21.1 | 0 |  | false | 0.00359 | 837 | 11.7 |
| V-BO-06 | enthalpy | off | Success | 1.0e-6 | 100 | 1.36 | 1001 | 500.0 | 1380.0 | 1 | 0.0655 | true | 0.348 | 1626 | 13.4 |
| V-BO-06 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0109 | 1001 | 500.0 | 11.1 | 1 | 5.54e-5 | true | 0.348 | 1626 | 9.93 |
| V-BO-07 | enthalpy | off | Success | 1.0e-6 | 100 | 1.37 | 1001 | 500.0 | 1420.0 | 1 | 0.0818 | true | 0.354 | 1166 | 31.9 |
| V-BO-07 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0242 | 1001 | 500.0 | 25.1 | 1 | 3.54e-5 | true | 0.354 | 1166 | 19.7 |
| V-BO-08 | enthalpy | off | Success | 1.0e-6 | 100 | 0.00982 | 301 | 150.0 | 10.3 | 1 | 0.00119 | true | 0.0681 | 278 | 2.16 |
| V-BO-08 | enthalpy | on | Success | 1.0e-6 | 100 | 0.00467 | 301 | 150.0 | 4.89 | 1 | 0.000117 | true | 0.0681 | 276 | 2.08 |
| V-CPL-01 | enthalpy |  | Success | 1.0e-6 | 100 | 0.000165 | 31 | 300.0 | 0.376 | 0 |  | true | 0.0 | 391 | 0.901 |
| V-CPL-01b | enthalpy |  | Success | 1.0e-6 | 100 | 1.43e-5 | 31 | 150.0 | 0.0296 | 0 |  | true | 0.0 | 1144 | 1.21 |
| V-CPL-02 | enthalpy | off | Success | 1.0e-5 | 100 | 1.21 | 41 | 2.0 | 164.0 | 19 | 0.0325 | true | 0.976 | 1223 | 5.46 |
| V-CPL-02 | enthalpy | on | Success | 1.0e-5 | 100 | 0.00604 | 41 | 2.0 | 0.821 | 19 | 1.8e-5 | true | 0.858 | 1280 | 3.18 |
| V-CPL-02b | enthalpy | off | Success | 1.0e-5 | 100 | 0.433 | 41 | 2.0 | 56.3 | 11 | 0.0321 | true | 0.341 | 663 | 0.964 |
| V-CPL-02b | enthalpy | on | Success | 1.0e-5 | 100 | 0.0224 | 41 | 2.0 | 2.91 | 11 | 2.74e-5 | true | 0.341 | 673 | 1.05 |
| V-CPL-03 | enthalpy |  | Success | 1.0e-8 | 100 | 5.68e-13 | 3 | 8.0 | 1.42e-7 | 0 |  | true | 1.28e-9 | 89 | 0.00606 |
| V-CPL-04 | enthalpy |  | Success | 1.0e-6 | 100 | 6.27e-7 | 11 | 2500.0 | 0.00124 | 0 |  | true | 0.0 | 31 | 0.00567 |
| V-CPL-05 | enthalpy |  | Success | 1.0e-8 | 10 | 1.89e-6 | 72 | 720.0 | 0.287 | 0 |  | true | 0.0 | 453 | 0.0781 |
| V-GT-01 | enthalpy |  | Success | 1.0e-10 | 100 | 2.07e-11 | 4 | 8.0 | 0.000519 | 0 |  | true | 0.0 | 47 | 0.0395 |
| V-GT-02 | enthalpy |  | Success | 1.0e-8 | 100 | 6.82e-13 | 3 | 1.0 | 1.36e-7 | 0 |  | true | 0.0 | 52 | 0.0399 |
| V-HC-01 | enthalpy |  | Success | 1.0e-8 | 100 | 4.98e-9 | 5 | 180.0 | 0.000651 | 0 |  | true | 0.0 | 294 | 0.206 |
| V-HC-02a | enthalpy |  | Success | 1.0e-8 | 100 | 1.69e-8 | 5 | 10000.0 | 0.00536 | 0 |  | true | 0.0 | 665 | 0.723 |
| V-HC-02b | enthalpy |  | Success | 1.0e-8 | 100 | 2.22e-8 | 5 | 10000.0 | 0.00557 | 0 |  | true | 0.0 | 547 | 0.451 |
| V-HC-02c | enthalpy |  | Success | 1.0e-8 | 100 | 1.47e-8 | 5 | 10000.0 | 0.00302 | 0 |  | true | 0.0 | 487 | 0.485 |
| V-HC-02d | enthalpy |  | Success | 1.0e-8 | 100 | 3.17e-8 | 5 | 10000.0 | 0.00637 | 0 |  | true | 0.0 | 468 | 0.449 |
| V-HC-03 | enthalpy |  | Success | 1.0e-8 | 100 | 5.83e-9 | 4 | 1200.0 | 0.000876 | 0 |  | true | 0.0 | 511 | 0.527 |
| V-HC-04a | enthalpy |  | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 | 0 |  | true | 0.0 | 504 | 0.222 |
| V-HC-04b | enthalpy |  | Success | 1.0e-8 | 100 | 1.51e-8 | 3 | 400000.0 | 0.0038 | 0 |  | true | 0.0 | 655 | 0.453 |
| V-HC-04c | enthalpy |  | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 | 0 |  | true | 0.0 | 504 | 0.229 |
| V-HC-05 | enthalpy |  | Success | 1.0e-8 | 100 | 3.8e-7 | 3 | 5000.0 | 0.0477 | 0 |  | true | 0.0 | 518 | 0.241 |
| V-MMS-01 | enthalpy |  | Success | 1.0e-8 | 100 | 1.38e-8 | 3 | 20.0 | 0.00205 | 0 |  | true | 0.0 | 126 | 0.517 |
| V-MMS-02 | enthalpy |  | Success | 1.0e-8 | 100 | 7.96e-13 | 41 | 10.0 | 1.59e-7 | 0 |  | true | 0.00196 | 27 | 0.13 |
| V-MMS-03 | enthalpy |  | Success | 1.0e-8 | 100 | 7.81e-6 | 41 | 10.0 | 1.26 | 0 |  | true | 4.35e-7 | 88 | 0.189 |
| V-MMS-04 | enthalpy |  | Success | 1.0e-8 | 100 | 1.38e-8 | 41 | 20.0 | 0.00224 | 0 |  | true | 1.94e-7 | 67 | 1.04 |
| V-RAD-01 | enthalpy |  | Success | 1.0e-8 | 100 | 4.16e-8 | 5 | 180.0 | 0.00573 | 0 |  | true | 0.0 | 366 | 0.351 |
| V-RX-01 | enthalpy |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 20.0 | 2.84e-6 | 0 |  | true | 0.0 | 174 | 0.00284 |
| V-RX-02a | enthalpy |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 | 0 |  | true | 0.0 | 14 | 0.000776 |
| V-RX-02b | enthalpy |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 | 0 |  | true | 0.0 | 243 | 0.00318 |
| V-RX-02c | enthalpy |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 7.17 | 2.53e-6 | 0 |  | true | 0.0 | 184 | 0.00292 |
| V-RX-02d | enthalpy |  | Success | 1.0e-10 | 100 | 2.27e-13 | 21 | 40.0 | 5.68e-6 | 0 |  | true | 0.0 | 313 | 0.00555 |
| V-RX-03 | enthalpy |  | Success | 1.0e-10 | 100 | 1.71e-13 | 7 | 15.0 | 4.26e-6 | 0 |  | true | 0.0 | 167 | 0.00379 |
| V-RX-04a | enthalpy |  | Success | 1.0e-8 | 100 | 2.16e-7 | 6 | 8.0 | 0.024 | 0 |  | true | 0.0 | 87 | 0.00278 |
| V-RX-04b | enthalpy |  | Success | 1.0e-8 | 100 | 9.46e-9 | 3 | 2000.0 | 0.00129 | 0 |  | true | 0.0 | 297 | 0.0445 |

92 of 100 solves within budget; 100 of 100 with every constraint residual at or below one.

## Item 7: event timing with location off and on

Per deactivating case and form, the largest shift of the invocation times between the solve and its reference [s] with location off and on, the median interval between the reference's invocations, or between its saved times for a single invocation [s], and the shifts as fractions of it, and gain, the shift off over the shift on. Passes when the located shift is at most 0.05 of the interval and the gain at least 3.0, both solves completing as the case expects.

| case | form | invocations | shift off [s] | shift on [s] | interval [s] | off / interval | on / interval | gain | passes |
|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | 1 | 0.411 | 0.003 | 0.5 | 0.822 | 0.00601 | 137.0 | true |
| V-BO-02 | temperature | 60 | 0.799 | 0.00211 | 6.68 | 0.12 | 0.000316 | 378.0 | true |
| V-BO-03 | temperature | 0 |  |  |  |  |  |  | false |
| V-BO-04 | temperature | 1 | 2.62 | 5.4e-8 | 2.0 | 1.31 | 2.7e-8 | 4.86e7 | true |
| V-BO-05 | temperature | 0 |  |  |  |  |  |  | false |
| V-BO-06 | temperature | 1 | 0.106 | 0.00266 | 0.5 | 0.213 | 0.00533 | 39.9 | true |
| V-BO-07 | temperature | 1 | 0.115 | 0.00166 | 0.5 | 0.231 | 0.00333 | 69.3 | true |
| V-BO-08 | temperature | 1 | 0.111 | 3.09e-5 | 0.5 | 0.222 | 6.18e-5 | 3600.0 | true |
| V-CPL-02 | temperature | 19 | 0.00297 | 1.78e-5 | 0.0966 | 0.0307 | 0.000184 | 167.0 | true |
| V-CPL-02b | temperature | 11 | 0.00522 | 8.11e-6 | 0.158 | 0.033 | 5.13e-5 | 643.0 | true |
| V-BO-01 | enthalpy | 1 | 0.275 | 1.68e-6 | 0.5 | 0.55 | 3.36e-6 | 164000.0 | true |
| V-BO-02 | enthalpy | 60 | 0.629 | 6.89e-5 | 6.68 | 0.0941 | 1.03e-5 | 9120.0 | true |
| V-BO-03 | enthalpy | 0 |  |  |  |  |  |  | false |
| V-BO-04 | enthalpy | 1 | 0.352 | 4.22e-7 | 2.0 | 0.176 | 2.11e-7 | 833000.0 | true |
| V-BO-05 | enthalpy | 0 |  |  |  |  |  |  | false |
| V-BO-06 | enthalpy | 1 | 0.0328 | 2.77e-5 | 0.5 | 0.0655 | 5.54e-5 | 1180.0 | true |
| V-BO-07 | enthalpy | 1 | 0.0409 | 1.77e-5 | 0.5 | 0.0818 | 3.54e-5 | 2310.0 | true |
| V-BO-08 | enthalpy | 1 | 0.000595 | 5.83e-5 | 0.5 | 0.00119 | 0.000117 | 10.2 | true |
| V-CPL-02 | enthalpy | 19 | 0.00314 | 1.74e-6 | 0.0966 | 0.0325 | 1.8e-5 | 1800.0 | true |
| V-CPL-02b | enthalpy | 11 | 0.00509 | 4.34e-6 | 0.158 | 0.0322 | 2.74e-5 | 1170.0 | true |

16 of 20 pass.
