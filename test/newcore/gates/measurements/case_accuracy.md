# Gate items 6 and 8 (Phase 1) and 7 (Phase 2): constraint residuals, achieved error and event timing

Produced by `test/newcore/gates/case_accuracy.jl` at commit 08de603 (tree 8edc104). Julia 1.13.0, 1 thread, 22 logical CPUs; load average [18.1, 17.7, 17.6] at the start and [17.4, 16.7, 16.6] at the end. Every case of `NEW_CORE_CASES`, each at its CI resolution and solve keywords, on the default integrator. Error: the largest temperature difference [K] over the saved times and the cells live in both solves, against the same case at both tolerances divided by 100, or by 10 where that does not complete (tightened) (under deactivation the saved times with the same live cells in both; a case saved at every step its final state; when one solve stops before the other, the saved times both reached; compared: the saved times compared, zero for the final state; to: the last time compared [s]); ratio: that error over reltol times the reference's largest temperature, within budget at 10 or less. Under deactivation each case is solved with event location off and on (location; its reference likewise) and judged on item 7 below: invocations, the deactivation invocations; event shift, the largest difference of the invocation times from the reference's over the median interval between the reference's invocations (between its saved times for a single invocation). Constraint: the largest scaled constraint residual of the accepted states (one or below passes; zero for a fixed-step solve or a configuration without algebraic rows). Wall: one warm solve [s].

| case | form | location | retcode | reltol | tightened | error [K] | compared | to [s] | ratio | invocations | event shift | within budget | constraint | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | off | Success | 1.0e-6 | 100 | 4.33 | 1000 | 500.0 | 4490.0 | 1 | 0.454 | true | 0.609 | 1126 | 4.31 |
| V-BO-01 | temperature | on | Success | 1.0e-6 | 100 | 0.248 | 1001 | 500.0 | 258.0 | 1 | 0.00601 | true | 0.609 | 1127 | 2.69 |
| V-BO-02 | temperature | off | Success | 1.0e-6 | 100 | 26.1 | 940 | 500.0 | 28600.0 | 60 | 0.12 | true | 0.0 | 1944 | 1.9 |
| V-BO-02 | temperature | on | Success | 1.0e-6 | 100 | 0.179 | 1001 | 500.0 | 197.0 | 60 | 0.000316 | true | 0.0 | 1977 | 2.12 |
| V-BO-03 | temperature | off | Success | 1.0e-6 | 100 | 0.0912 | 1001 | 500.0 | 133.0 | 0 |  | false | 0.802 | 1064 | 2.74 |
| V-BO-03 | temperature | on | Success | 1.0e-6 | 100 | 0.0912 | 1001 | 500.0 | 133.0 | 0 |  | false | 0.802 | 1064 | 2.56 |
| V-BO-04 | temperature | off | Success | 1.0e-7 | 100 | 2.16e-12 | 240 | 480.0 | 3.6e-8 | 1 | 1.31 | true | 6.15e-10 | 321 | 0.106 |
| V-BO-04 | temperature | on | Success | 1.0e-7 | 100 | 2.16e-12 | 241 | 480.0 | 3.6e-8 | 1 | 2.7e-8 | true | 6.15e-10 | 322 | 0.11 |
| V-BO-05 | temperature | off | Terminated | 1.0e-6 | 100 | 0.138 | 520 | 260.0 | 199.0 | 0 |  | false | 0.0184 | 347 | 0.815 |
| V-BO-05 | temperature | on | Terminated | 1.0e-6 | 100 | 0.138 | 520 | 260.0 | 199.0 | 0 |  | false | 0.0184 | 347 | 0.549 |
| V-BO-06 | temperature | off | Success | 1.0e-6 | 100 | 5.2 | 1000 | 500.0 | 5270.0 | 1 | 0.528 | true | 0.664 | 1138 | 2.31 |
| V-BO-06 | temperature | on | Success | 1.0e-6 | 100 | 0.268 | 1001 | 500.0 | 272.0 | 1 | 0.00534 | true | 0.664 | 1137 | 3.44 |
| V-BO-07 | temperature | off | Success | 1.0e-6 | 100 | 0.457 | 1001 | 500.0 | 474.0 | 1 | 0.0263 | true | 0.684 | 733 | 6.04 |
| V-BO-07 | temperature | on | Success | 1.0e-6 | 100 | 0.177 | 1001 | 500.0 | 184.0 | 1 | 0.00329 | true | 0.684 | 732 | 4.14 |
| V-BO-08 | temperature | off | Success | 1.0e-6 | 100 | 14.1 | 1500 | 150.0 | 14800.0 | 1 | 1.11 | true | 0.0498 | 213 | 0.619 |
| V-BO-08 | temperature | on | Success | 1.0e-6 | 100 | 0.0262 | 1501 | 150.0 | 27.4 | 1 | 0.000309 | true | 0.0498 | 213 | 0.582 |
| V-CPL-01 | temperature |  | Success | 1.0e-6 | 100 | 0.000421 | 31 | 300.0 | 0.96 | 0 |  | true | 0.0 | 246 | 0.109 |
| V-CPL-01b | temperature |  | Success | 1.0e-6 | 100 | 2.83e-5 | 31 | 150.0 | 0.0586 | 0 |  | true | 0.0 | 1141 | 0.785 |
| V-CPL-02 | temperature | off | Success | 1.0e-5 | 100 | 0.32 | 41 | 2.0 | 43.5 | 19 | 0.0316 | true | 0.855 | 1127 | 2.06 |
| V-CPL-02 | temperature | on | Success | 1.0e-5 | 100 | 0.0284 | 41 | 2.0 | 3.86 | 19 | 0.000184 | true | 0.999 | 1134 | 2.53 |
| V-CPL-02b | temperature | off | Success | 1.0e-5 | 100 | 1.52 | 40 | 2.0 | 197.0 | 11 | 0.0343 | true | 0.189 | 636 | 0.776 |
| V-CPL-02b | temperature | on | Success | 1.0e-5 | 100 | 0.0161 | 41 | 2.0 | 2.09 | 11 | 5.13e-5 | true | 0.197 | 638 | 0.85 |
| V-CPL-03 | temperature |  | Success | 1.0e-8 | 100 | 0.0 | 3 | 8.0 | 0.0 | 0 |  | true | 1.05e-9 | 88 | 0.00688 |
| V-CPL-04 | temperature |  | Success | 1.0e-6 | 100 | 1.13e-6 | 11 | 2500.0 | 0.00224 | 0 |  | true | 0.0 | 28 | 0.0046 |
| V-CPL-05 | temperature |  | Success | 1.0e-8 | 10 | 2.38e-6 | 72 | 720.0 | 0.361 | 0 |  | true | 0.0 | 237 | 0.0246 |
| V-GT-01 | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 4 | 8.0 | 0.0 | 0 |  | true | 0.0 | 46 | 0.0254 |
| V-GT-02 | temperature |  | Success | 1.0e-8 | 100 | 0.0 | 3 | 1.0 | 0.0 | 0 |  | true | 0.0 | 52 | 0.0193 |
| V-HC-01 | temperature |  | Success | 1.0e-8 | 100 | 3.03e-7 | 5 | 180.0 | 0.0395 | 0 |  | true | 0.0 | 60 | 0.0182 |
| V-HC-02a | temperature |  | Success | 1.0e-8 | 100 | 4.51e-8 | 5 | 10000.0 | 0.0143 | 0 |  | true | 0.0 | 56 | 0.0171 |
| V-HC-02b | temperature |  | Success | 1.0e-8 | 100 | 9.55e-8 | 5 | 10000.0 | 0.024 | 0 |  | true | 0.0 | 75 | 0.0218 |
| V-HC-02c | temperature |  | Success | 1.0e-8 | 100 | 1.24e-7 | 5 | 10000.0 | 0.0256 | 0 |  | true | 0.0 | 98 | 0.0297 |
| V-HC-02d | temperature |  | Success | 1.0e-8 | 100 | 1.28e-7 | 5 | 10000.0 | 0.0258 | 0 |  | true | 0.0 | 113 | 0.0335 |
| V-HC-03 | temperature |  | Success | 1.0e-8 | 100 | 2.43e-8 | 4 | 1200.0 | 0.00366 | 0 |  | true | 0.0 | 86 | 0.00982 |
| V-HC-04a | temperature |  | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 | 0 |  | true | 0.0 | 99 | 0.00906 |
| V-HC-04b | temperature |  | Success | 1.0e-8 | 100 | 3.98e-13 | 3 | 400000.0 | 1.0e-7 | 0 |  | true | 0.0 | 97 | 0.0085 |
| V-HC-04c | temperature |  | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 | 0 |  | true | 0.0 | 99 | 0.0089 |
| V-HC-05 | temperature |  | Success | 1.0e-8 | 100 | 1.14e-13 | 3 | 5000.0 | 1.43e-8 | 0 |  | true | 0.0 | 243 | 0.0272 |
| V-MMS-01 | temperature |  | Success | 1.0e-8 | 100 | 2.2e-8 | 3 | 20.0 | 0.00326 | 0 |  | true | 0.0 | 127 | 0.433 |
| V-MMS-02 | temperature |  | Success | 1.0e-8 | 100 | 3.41e-13 | 41 | 10.0 | 6.82e-8 | 0 |  | true | 0.00178 | 26 | 0.0204 |
| V-MMS-03 | temperature |  | Success | 1.0e-8 | 100 | 1.09e-5 | 41 | 10.0 | 1.75 | 0 |  | true | 1.43e-6 | 80 | 0.847 |
| V-MMS-04 | temperature |  | Success | 1.0e-8 | 100 | 8.98e-9 | 41 | 20.0 | 0.00146 | 0 |  | true | 9.02e-9 | 63 | 0.564 |
| V-RAD-01 | temperature |  | Success | 1.0e-8 | 100 | 2.08e-7 | 5 | 180.0 | 0.0286 | 0 |  | true | 0.0 | 63 | 0.0148 |
| V-RX-01 | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 7 | 20.0 | 0.0 | 0 |  | true | 0.0 | 175 | 0.0033 |
| V-RX-02a | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 | 0 |  | true | 0.0 | 13 | 0.000931 |
| V-RX-02b | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 | 0 |  | true | 0.0 | 242 | 0.00275 |
| V-RX-02c | temperature |  | Success | 1.0e-10 | 100 | 0.0 | 7 | 7.17 | 0.0 | 0 |  | true | 0.0 | 183 | 0.00226 |
| V-RX-02d | temperature |  | Success | 1.0e-10 | 100 | 5.68e-14 | 21 | 40.0 | 1.42e-6 | 0 |  | true | 0.0 | 306 | 0.00343 |
| V-RX-03 | temperature |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 15.0 | 2.84e-6 | 0 |  | true | 0.0 | 160 | 0.00251 |
| V-RX-04a | temperature |  | Success | 1.0e-8 | 100 | 3.02e-7 | 6 | 8.0 | 0.0336 | 0 |  | true | 0.0 | 78 | 0.00224 |
| V-RX-04b | temperature |  | Success | 1.0e-8 | 100 | 5.97e-9 | 3 | 2000.0 | 0.000814 | 0 |  | true | 0.0 | 159 | 0.0087 |
| V-BO-01 | enthalpy | off | Success | 1.0e-6 | 100 | 5.06 | 1000 | 500.0 | 5250.0 | 1 | 0.533 | true | 0.382 | 1574 | 14.7 |
| V-BO-01 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0308 | 1001 | 500.0 | 31.9 | 1 | 1.03e-5 | true | 0.382 | 1572 | 15.6 |
| V-BO-02 | enthalpy | off | Success | 1.0e-6 | 100 | 16.2 | 958 | 500.0 | 17800.0 | 60 | 0.0856 | true | 0.0 | 2778 | 16.7 |
| V-BO-02 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0147 | 1001 | 500.0 | 16.1 | 60 | 1.08e-5 | true | 0.0 | 2797 | 17.0 |
| V-BO-03 | enthalpy | off | Success | 1.0e-6 | 100 | 0.0093 | 1001 | 500.0 | 13.6 | 0 |  | false | 0.423 | 1576 | 17.1 |
| V-BO-03 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0093 | 1001 | 500.0 | 13.6 | 0 |  | false | 0.423 | 1576 | 21.5 |
| V-BO-04 | enthalpy | off | Success | 1.0e-7 | 100 | 5.89e-7 | 241 | 480.0 | 0.00981 | 1 | 0.189 | true | 1.65e-8 | 322 | 0.216 |
| V-BO-04 | enthalpy | on | Success | 1.0e-7 | 100 | 5.89e-7 | 241 | 480.0 | 0.00981 | 1 | 2.93e-7 | true | 1.65e-8 | 322 | 0.221 |
| V-BO-05 | enthalpy | off | Terminated | 1.0e-6 | 100 | 0.0332 | 520 | 260.0 | 47.9 | 0 |  | false | 0.00351 | 825 | 9.67 |
| V-BO-05 | enthalpy | on | Terminated | 1.0e-6 | 100 | 0.0332 | 520 | 260.0 | 47.9 | 0 |  | false | 0.00351 | 825 | 9.77 |
| V-BO-06 | enthalpy | off | Success | 1.0e-6 | 100 | 20.6 | 1001 | 500.0 | 20800.0 | 1 | 0.56 | true | 0.436 | 1627 | 12.4 |
| V-BO-06 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0137 | 1001 | 500.0 | 13.9 | 1 | 5.25e-5 | true | 0.436 | 1627 | 13.6 |
| V-BO-07 | enthalpy | off | Success | 1.0e-6 | 100 | 6.2 | 1001 | 500.0 | 6420.0 | 1 | 0.339 | true | 0.384 | 1150 | 19.0 |
| V-BO-07 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0133 | 1001 | 500.0 | 13.8 | 1 | 1.48e-5 | true | 0.384 | 1150 | 18.9 |
| V-BO-08 | enthalpy | off | Success | 1.0e-6 | 100 | 5.46 | 1499 | 150.0 | 5720.0 | 1 | 1.3 | true | 0.0982 | 284 | 2.41 |
| V-BO-08 | enthalpy | on | Success | 1.0e-6 | 100 | 0.0133 | 1501 | 150.0 | 13.9 | 1 | 0.000571 | true | 0.0982 | 283 | 2.41 |
| V-CPL-01 | enthalpy |  | Success | 1.0e-6 | 100 | 0.000165 | 31 | 300.0 | 0.376 | 0 |  | true | 0.0 | 391 | 0.663 |
| V-CPL-01b | enthalpy |  | Success | 1.0e-6 | 100 | 1.43e-5 | 31 | 150.0 | 0.0296 | 0 |  | true | 0.0 | 1144 | 1.5 |
| V-CPL-02 | enthalpy | off | Success | 1.0e-5 | 100 | 1.11 | 41 | 2.0 | 151.0 | 19 | 0.0305 | true | 0.976 | 1223 | 4.44 |
| V-CPL-02 | enthalpy | on | Success | 1.0e-5 | 100 | 0.00604 | 41 | 2.0 | 0.821 | 19 | 1.8e-5 | true | 0.858 | 1280 | 4.66 |
| V-CPL-02b | enthalpy | off | Success | 1.0e-5 | 100 | 0.31 | 41 | 2.0 | 40.4 | 11 | 0.0321 | true | 0.341 | 663 | 1.45 |
| V-CPL-02b | enthalpy | on | Success | 1.0e-5 | 100 | 0.0224 | 41 | 2.0 | 2.91 | 11 | 2.74e-5 | true | 0.341 | 673 | 1.5 |
| V-CPL-03 | enthalpy |  | Success | 1.0e-8 | 100 | 5.68e-13 | 3 | 8.0 | 1.42e-7 | 0 |  | true | 1.28e-9 | 89 | 0.0088 |
| V-CPL-04 | enthalpy |  | Success | 1.0e-6 | 100 | 6.27e-7 | 11 | 2500.0 | 0.00124 | 0 |  | true | 0.0 | 31 | 0.00904 |
| V-CPL-05 | enthalpy |  | Success | 1.0e-8 | 10 | 1.89e-6 | 72 | 720.0 | 0.287 | 0 |  | true | 0.0 | 453 | 0.1 |
| V-GT-01 | enthalpy |  | Success | 1.0e-10 | 100 | 2.07e-11 | 4 | 8.0 | 0.000519 | 0 |  | true | 0.0 | 47 | 0.0598 |
| V-GT-02 | enthalpy |  | Success | 1.0e-8 | 100 | 6.82e-13 | 3 | 1.0 | 1.36e-7 | 0 |  | true | 0.0 | 52 | 0.0313 |
| V-HC-01 | enthalpy |  | Success | 1.0e-8 | 100 | 4.98e-9 | 5 | 180.0 | 0.000651 | 0 |  | true | 0.0 | 294 | 0.363 |
| V-HC-02a | enthalpy |  | Success | 1.0e-8 | 100 | 1.69e-8 | 5 | 10000.0 | 0.00536 | 0 |  | true | 0.0 | 665 | 0.726 |
| V-HC-02b | enthalpy |  | Success | 1.0e-8 | 100 | 2.22e-8 | 5 | 10000.0 | 0.00557 | 0 |  | true | 0.0 | 547 | 0.542 |
| V-HC-02c | enthalpy |  | Success | 1.0e-8 | 100 | 1.47e-8 | 5 | 10000.0 | 0.00302 | 0 |  | true | 0.0 | 487 | 0.444 |
| V-HC-02d | enthalpy |  | Success | 1.0e-8 | 100 | 3.17e-8 | 5 | 10000.0 | 0.00637 | 0 |  | true | 0.0 | 468 | 0.409 |
| V-HC-03 | enthalpy |  | Success | 1.0e-8 | 100 | 5.83e-9 | 4 | 1200.0 | 0.000876 | 0 |  | true | 0.0 | 511 | 0.488 |
| V-HC-04a | enthalpy |  | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 | 0 |  | true | 0.0 | 504 | 0.314 |
| V-HC-04b | enthalpy |  | Success | 1.0e-8 | 100 | 1.51e-8 | 3 | 400000.0 | 0.0038 | 0 |  | true | 0.0 | 655 | 0.393 |
| V-HC-04c | enthalpy |  | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 | 0 |  | true | 0.0 | 504 | 0.288 |
| V-HC-05 | enthalpy |  | Success | 1.0e-8 | 100 | 3.8e-7 | 3 | 5000.0 | 0.0477 | 0 |  | true | 0.0 | 518 | 0.406 |
| V-MMS-01 | enthalpy |  | Success | 1.0e-8 | 100 | 1.38e-8 | 3 | 20.0 | 0.00205 | 0 |  | true | 0.0 | 126 | 0.684 |
| V-MMS-02 | enthalpy |  | Success | 1.0e-8 | 100 | 1.19e-12 | 41 | 10.0 | 2.39e-7 | 0 |  | true | 0.00196 | 27 | 0.0165 |
| V-MMS-03 | enthalpy |  | Success | 1.0e-8 | 100 | 3.46e-6 | 41 | 10.0 | 0.557 | 0 |  | true | 5.84e-7 | 91 | 0.396 |
| V-MMS-04 | enthalpy |  | Success | 1.0e-8 | 100 | 1.38e-8 | 41 | 20.0 | 0.00224 | 0 |  | true | 1.94e-7 | 67 | 0.434 |
| V-RAD-01 | enthalpy |  | Success | 1.0e-8 | 100 | 4.16e-8 | 5 | 180.0 | 0.00573 | 0 |  | true | 0.0 | 366 | 0.741 |
| V-RX-01 | enthalpy |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 20.0 | 2.84e-6 | 0 |  | true | 0.0 | 174 | 0.00313 |
| V-RX-02a | enthalpy |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 | 0 |  | true | 0.0 | 14 | 0.001 |
| V-RX-02b | enthalpy |  | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 | 0 |  | true | 0.0 | 243 | 0.00324 |
| V-RX-02c | enthalpy |  | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 7.17 | 2.53e-6 | 0 |  | true | 0.0 | 184 | 0.00278 |
| V-RX-02d | enthalpy |  | Success | 1.0e-10 | 100 | 2.27e-13 | 21 | 40.0 | 5.68e-6 | 0 |  | true | 0.0 | 313 | 0.00414 |
| V-RX-03 | enthalpy |  | Success | 1.0e-10 | 100 | 1.71e-13 | 7 | 15.0 | 4.26e-6 | 0 |  | true | 0.0 | 167 | 0.00407 |
| V-RX-04a | enthalpy |  | Success | 1.0e-8 | 100 | 2.16e-7 | 6 | 8.0 | 0.024 | 0 |  | true | 0.0 | 87 | 0.00333 |
| V-RX-04b | enthalpy |  | Success | 1.0e-8 | 100 | 9.46e-9 | 3 | 2000.0 | 0.00129 | 0 |  | true | 0.0 | 297 | 0.0544 |

92 of 100 solves within budget; 100 of 100 with every constraint residual at or below one.

## Item 7: event timing with location off and on

Per deactivating case and form, the largest shift of the invocation times between the solve and its reference [s] with location off and on, the median interval between the reference's invocations, or between its saved times for a single invocation [s], and the shifts as fractions of it, and gain, the shift off over the shift on. Passes when the located shift is at most 0.05 of the interval and the gain at least 3.0, both solves completing as the case expects.

| case | form | invocations | shift off [s] | shift on [s] | interval [s] | off / interval | on / interval | gain | passes |
|---|---|---|---|---|---|---|---|---|---|
| V-BO-01 | temperature | 1 | 0.227 | 0.00301 | 0.5 | 0.454 | 0.00601 | 75.6 | true |
| V-BO-02 | temperature | 60 | 0.799 | 0.00211 | 6.68 | 0.12 | 0.000316 | 378.0 | true |
| V-BO-03 | temperature | 0 |  |  |  |  |  |  | false |
| V-BO-04 | temperature | 1 | 2.62 | 5.4e-8 | 2.0 | 1.31 | 2.7e-8 | 4.85e7 | true |
| V-BO-05 | temperature | 0 |  |  |  |  |  |  | false |
| V-BO-06 | temperature | 1 | 0.264 | 0.00267 | 0.5 | 0.528 | 0.00534 | 98.9 | true |
| V-BO-07 | temperature | 1 | 0.0131 | 0.00165 | 0.5 | 0.0263 | 0.00329 | 7.99 | true |
| V-BO-08 | temperature | 1 | 0.111 | 3.09e-5 | 0.1 | 1.11 | 0.000309 | 3600.0 | true |
| V-CPL-02 | temperature | 19 | 0.00306 | 1.78e-5 | 0.0966 | 0.0317 | 0.000184 | 173.0 | true |
| V-CPL-02b | temperature | 11 | 0.00543 | 8.11e-6 | 0.158 | 0.0343 | 5.13e-5 | 670.0 | true |
| V-BO-01 | enthalpy | 1 | 0.267 | 5.15e-6 | 0.5 | 0.533 | 1.03e-5 | 51800.0 | true |
| V-BO-02 | enthalpy | 60 | 0.571 | 7.21e-5 | 6.68 | 0.0855 | 1.08e-5 | 7930.0 | true |
| V-BO-03 | enthalpy | 0 |  |  |  |  |  |  | false |
| V-BO-04 | enthalpy | 1 | 0.378 | 5.85e-7 | 2.0 | 0.189 | 2.93e-7 | 645000.0 | true |
| V-BO-05 | enthalpy | 0 |  |  |  |  |  |  | false |
| V-BO-06 | enthalpy | 1 | 0.28 | 2.62e-5 | 0.5 | 0.56 | 5.25e-5 | 10700.0 | true |
| V-BO-07 | enthalpy | 1 | 0.169 | 7.42e-6 | 0.5 | 0.339 | 1.48e-5 | 22800.0 | true |
| V-BO-08 | enthalpy | 1 | 0.13 | 5.71e-5 | 0.1 | 1.3 | 0.000571 | 2280.0 | true |
| V-CPL-02 | enthalpy | 19 | 0.00296 | 1.74e-6 | 0.0966 | 0.0306 | 1.8e-5 | 1700.0 | true |
| V-CPL-02b | enthalpy | 11 | 0.00509 | 4.34e-6 | 0.158 | 0.0322 | 2.74e-5 | 1170.0 | true |

16 of 20 pass.
