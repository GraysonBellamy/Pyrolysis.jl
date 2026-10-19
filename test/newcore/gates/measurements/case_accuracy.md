# Gate items 6 and 8: constraint residuals and achieved error

Produced by `test/newcore/gates/case_accuracy.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [20.8, 21.2, 22.2] at the start and [21.3, 21.6, 21.9] at the end. Every case of `NEW_CORE_CASES` at its CI resolution and solve keywords, on the default integrator. Error: the largest temperature difference [K] over the saved times and the cells live in both solves, against the same case at both tolerances divided by 100, or by 10 where that does not complete (tightened) (under deactivation the saved times with the same live cells in both; a case saved at every step its final state; compared: the saved times compared, zero for the final state); ratio: that error over reltol times the reference's largest temperature, within budget at 10 or less; event shift: under deactivation the largest event-time difference over the reference's median event interval, over reltol, on which those cases' budget is judged. Constraint: the largest scaled constraint residual of the accepted states (one or below passes; zero for a fixed-step solve or a configuration without algebraic rows). Wall: one warm solve [s].

| case | form | retcode | reltol | tightened | error [K] | compared | ratio | event shift | within budget | constraint | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-02 | temperature | Success | 1.0e-6 | 100 | 26.1 | 939 | 28600.0 | 117000.0 | false | 0.0 | 2354 | 3.73 |
| V-CPL-01 | temperature | Success | 1.0e-6 | 100 | 0.000421 | 31 | 0.96 |  | true | 0.0 | 246 | 0.0634 |
| V-CPL-01b | temperature | Success | 1.0e-6 | 100 | 2.85e-5 | 31 | 0.0589 |  | true | 0.0 | 1141 | 0.529 |
| V-CPL-02 | temperature | Success | 1.0e-5 | 100 | 0.0379 | 41 | 5.13 | 6480.0 | false | 0.244 | 1904 | 3.52 |
| V-CPL-02b | temperature | Success | 1.0e-5 | 100 | 0.128 | 41 | 16.5 | 6740.0 | false | 0.0495 | 1057 | 0.803 |
| V-CPL-03 | temperature | Success | 1.0e-8 | 100 | 0.0 | 3 | 0.0 |  | true | 1.05e-9 | 88 | 0.00525 |
| V-CPL-04 | temperature | Success | 1.0e-6 | 100 | 1.13e-6 | 11 | 0.00224 |  | true | 0.0 | 28 | 0.00288 |
| V-CPL-05 | temperature | Success | 1.0e-8 | 10 |  | 0 |  |  | false | 0.0 | 237 | 0.0155 |
| V-GT-01 | temperature | Success | 1.0e-10 | 100 | 0.0 | 4 | 0.0 |  | true | 0.0 | 46 | 0.0213 |
| V-GT-02 | temperature | Success | 1.0e-8 | 100 | 0.0 | 3 | 0.0 |  | true | 0.0 | 52 | 0.0185 |
| V-HC-01 | temperature | Success | 1.0e-8 | 100 | 3.03e-7 | 5 | 0.0395 |  | true | 0.0 | 60 | 0.00958 |
| V-HC-02a | temperature | Success | 1.0e-8 | 100 | 4.51e-8 | 5 | 0.0143 |  | true | 0.0 | 56 | 0.00711 |
| V-HC-02b | temperature | Success | 1.0e-8 | 100 | 9.55e-8 | 5 | 0.024 |  | true | 0.0 | 75 | 0.0089 |
| V-HC-02c | temperature | Success | 1.0e-8 | 100 | 1.24e-7 | 5 | 0.0256 |  | true | 0.0 | 98 | 0.0113 |
| V-HC-02d | temperature | Success | 1.0e-8 | 100 | 1.28e-7 | 5 | 0.0258 |  | true | 0.0 | 113 | 0.0127 |
| V-HC-03 | temperature | Success | 1.0e-8 | 100 | 2.43e-8 | 4 | 0.00366 |  | true | 0.0 | 86 | 0.00999 |
| V-HC-04a | temperature | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 4.29e-8 |  | true | 0.0 | 99 | 0.00964 |
| V-HC-04b | temperature | Success | 1.0e-8 | 100 | 3.98e-13 | 3 | 1.0e-7 |  | true | 0.0 | 97 | 0.00937 |
| V-HC-04c | temperature | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 4.29e-8 |  | true | 0.0 | 99 | 0.00949 |
| V-HC-05 | temperature | Success | 1.0e-8 | 100 | 1.14e-13 | 3 | 1.43e-8 |  | true | 0.0 | 243 | 0.0267 |
| V-MMS-01 | temperature | Success | 1.0e-8 | 100 | 2.31e-8 | 3 | 0.00342 |  | true | 0.0 | 127 | 0.355 |
| V-MMS-02 | temperature | Success | 1.0e-8 | 100 | 3.41e-13 | 41 | 6.82e-8 |  | true | 0.00178 | 26 | 0.0114 |
| V-MMS-03 | temperature | Success | 1.0e-8 | 100 | 1.09e-5 | 41 | 1.75 |  | true | 1.43e-6 | 80 | 0.42 |
| V-MMS-04 | temperature | Success | 1.0e-8 | 100 | 8.6e-9 | 41 | 0.0014 |  | true | 9.13e-9 | 63 | 0.29 |
| V-RAD-01 | temperature | Success | 1.0e-8 | 100 | 2.08e-7 | 5 | 0.0286 |  | true | 0.0 | 63 | 0.0127 |
| V-RX-01 | temperature | Success | 1.0e-10 | 100 | 0.0 | 7 | 0.0 |  | true | 0.0 | 175 | 0.00236 |
| V-RX-02a | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 1.42e-6 |  | true | 0.0 | 13 | 0.000728 |
| V-RX-02b | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 1.42e-6 |  | true | 0.0 | 242 | 0.00225 |
| V-RX-02c | temperature | Success | 1.0e-10 | 100 | 0.0 | 7 | 0.0 |  | true | 0.0 | 183 | 0.00183 |
| V-RX-02d | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 21 | 1.42e-6 |  | true | 0.0 | 306 | 0.00288 |
| V-RX-03 | temperature | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 2.84e-6 |  | true | 0.0 | 160 | 0.00203 |
| V-RX-04a | temperature | Success | 1.0e-8 | 100 | 3.02e-7 | 6 | 0.0336 |  | true | 0.0 | 78 | 0.00179 |
| V-RX-04b | temperature | Success | 1.0e-8 | 100 | 5.97e-9 | 3 | 0.000814 |  | true | 0.0 | 159 | 0.00616 |
| V-BO-02 | enthalpy | Success | 1.0e-6 | 100 | 26.6 | 955 | 29200.0 | 93100.0 | false | 0.0 | 3163 | 28.4 |
| V-CPL-01 | enthalpy | Success | 1.0e-6 | 100 | 0.000165 | 31 | 0.376 |  | true | 0.0 | 391 | 0.478 |
| V-CPL-01b | enthalpy | Success | 1.0e-6 | 100 | 1.43e-5 | 31 | 0.0296 |  | true | 0.0 | 1144 | 0.947 |
| V-CPL-02 | enthalpy | Success | 1.0e-5 | 100 | 0.0362 | 41 | 4.9 | 7250.0 | false | 0.422 | 1901 | 4.82 |
| V-CPL-02b | enthalpy | Success | 1.0e-5 | 100 | 0.614 | 41 | 79.5 | 5390.0 | false | 0.0791 | 1114 | 2.61 |
| V-CPL-03 | enthalpy | Success | 1.0e-8 | 100 | 4.55e-13 | 3 | 1.14e-7 |  | true | 1.28e-9 | 89 | 0.00822 |
| V-CPL-04 | enthalpy | Success | 1.0e-6 | 100 | 6.27e-7 | 11 | 0.00124 |  | true | 0.0 | 31 | 0.00578 |
| V-CPL-05 | enthalpy | Success | 1.0e-8 | 10 |  | 0 |  |  | false | 0.0 | 455 | 0.0619 |
| V-GT-01 | enthalpy | Success | 1.0e-10 | 100 | 2.07e-11 | 4 | 0.000517 |  | true | 0.0 | 47 | 0.0386 |
| V-GT-02 | enthalpy | Success | 1.0e-8 | 100 | 7.39e-13 | 3 | 1.48e-7 |  | true | 0.0 | 52 | 0.0216 |
| V-HC-01 | enthalpy | Success | 1.0e-8 | 100 | 4.98e-9 | 5 | 0.000651 |  | true | 0.0 | 294 | 0.214 |
| V-HC-02a | enthalpy | Success | 1.0e-8 | 100 | 1.69e-8 | 5 | 0.00536 |  | true | 0.0 | 665 | 0.42 |
| V-HC-02b | enthalpy | Success | 1.0e-8 | 100 | 2.22e-8 | 5 | 0.00557 |  | true | 0.0 | 547 | 0.316 |
| V-HC-02c | enthalpy | Success | 1.0e-8 | 100 | 1.47e-8 | 5 | 0.00302 |  | true | 0.0 | 487 | 0.258 |
| V-HC-02d | enthalpy | Success | 1.0e-8 | 100 | 3.17e-8 | 5 | 0.00637 |  | true | 0.0 | 468 | 0.243 |
| V-HC-03 | enthalpy | Success | 1.0e-8 | 100 | 5.83e-9 | 4 | 0.000876 |  | true | 0.0 | 511 | 0.288 |
| V-HC-04a | enthalpy | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 0.0267 |  | true | 0.0 | 504 | 0.168 |
| V-HC-04b | enthalpy | Success | 1.0e-8 | 100 | 1.51e-8 | 3 | 0.0038 |  | true | 0.0 | 655 | 0.23 |
| V-HC-04c | enthalpy | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 0.0267 |  | true | 0.0 | 504 | 0.169 |
| V-HC-05 | enthalpy | Success | 1.0e-8 | 100 | 3.8e-7 | 3 | 0.0477 |  | true | 0.0 | 518 | 0.234 |
| V-MMS-01 | enthalpy | Success | 1.0e-8 | 100 | 1.32e-8 | 3 | 0.00195 |  | true | 0.0 | 126 | 0.492 |
| V-MMS-02 | enthalpy | Success | 1.0e-8 | 100 | 1.19e-12 | 41 | 2.39e-7 |  | true | 0.00196 | 27 | 0.0129 |
| V-MMS-03 | enthalpy | Success | 1.0e-8 | 100 | 3.46e-6 | 41 | 0.557 |  | true | 5.84e-7 | 91 | 0.249 |
| V-MMS-04 | enthalpy | Success | 1.0e-8 | 100 | 1.35e-8 | 41 | 0.0022 |  | true | 1.95e-7 | 67 | 0.293 |
| V-RAD-01 | enthalpy | Success | 1.0e-8 | 100 | 4.16e-8 | 5 | 0.00573 |  | true | 0.0 | 366 | 0.267 |
| V-RX-01 | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 2.84e-6 |  | true | 0.0 | 174 | 0.00226 |
| V-RX-02a | enthalpy | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 1.42e-6 |  | true | 0.0 | 14 | 0.000871 |
| V-RX-02b | enthalpy | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 1.42e-6 |  | true | 0.0 | 243 | 0.00263 |
| V-RX-02c | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 2.53e-6 |  | true | 0.0 | 184 | 0.00234 |
| V-RX-02d | enthalpy | Success | 1.0e-10 | 100 | 2.27e-13 | 21 | 5.68e-6 |  | true | 0.0 | 313 | 0.00332 |
| V-RX-03 | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 2.84e-6 |  | true | 0.0 | 167 | 0.00294 |
| V-RX-04a | enthalpy | Success | 1.0e-8 | 100 | 2.16e-7 | 6 | 0.024 |  | true | 0.0 | 87 | 0.00229 |
| V-RX-04b | enthalpy | Success | 1.0e-8 | 100 | 9.46e-9 | 3 | 0.00129 |  | true | 0.0 | 297 | 0.0335 |

58 of 66 solves within budget; 66 of 66 with every constraint residual at or below one.
