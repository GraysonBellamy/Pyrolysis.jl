# Gate items 6 and 8: constraint residuals and achieved error

Produced by `test/newcore/gates/case_accuracy.jl` at commit 603e863. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [13.4, 12.4, 12.2] at the start and [14.3, 14.5, 13.5] at the end. Every case of `NEW_CORE_CASES` at its CI resolution and solve keywords, on the default integrator. Error: the largest temperature difference [K] over the saved times and the cells live in both solves, against the same case at both tolerances divided by 100, or by 10 where that does not complete (tightened) (under deactivation the saved times with the same live cells in both; a case saved at every step its final state; when one solve stops before the other, the saved times both reached; compared: the saved times compared, zero for the final state; to: the last time compared [s]); ratio: that error over reltol times the reference's largest temperature, within budget at 10 or less; event shift: under deactivation the largest event-time difference over the reference's median event interval, over reltol, on which those cases' budget is judged. Constraint: the largest scaled constraint residual of the accepted states (one or below passes; zero for a fixed-step solve or a configuration without algebraic rows). Wall: one warm solve [s].

| case | form | retcode | reltol | tightened | error [K] | compared | to [s] | ratio | event shift | within budget | constraint | steps | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V-BO-02 | temperature | Success | 1.0e-6 | 100 | 26.4 | 943 | 500.0 | 28900.0 | 133000.0 | false | 0.0 | 1706 | 1.56 |
| V-CPL-01 | temperature | Success | 1.0e-6 | 100 | 0.000421 | 31 | 300.0 | 0.96 |  | true | 0.0 | 246 | 0.0399 |
| V-CPL-01b | temperature | Success | 1.0e-6 | 100 | 2.85e-5 | 31 | 150.0 | 0.0589 |  | true | 0.0 | 1141 | 0.293 |
| V-CPL-02 | temperature | Success | 1.0e-5 | 100 | 0.0379 | 41 | 2.0 | 5.13 | 6480.0 | false | 0.244 | 1904 | 1.1 |
| V-CPL-02b | temperature | Success | 1.0e-5 | 100 | 0.128 | 41 | 2.0 | 16.6 | 6750.0 | false | 0.0495 | 1057 | 0.367 |
| V-CPL-03 | temperature | Success | 1.0e-8 | 100 | 0.0 | 3 | 8.0 | 0.0 |  | true | 9.86e-10 | 88 | 0.00397 |
| V-CPL-04 | temperature | Success | 1.0e-6 | 100 | 1.13e-6 | 11 | 2500.0 | 0.00224 |  | true | 0.0 | 28 | 0.00185 |
| V-CPL-05 | temperature | Success | 1.0e-8 | 10 | 1.45e-6 | 72 | 720.0 | 0.221 |  | true | 0.0 | 237 | 0.0111 |
| V-GT-01 | temperature | Success | 1.0e-10 | 100 | 0.0 | 4 | 8.0 | 0.0 |  | true | 0.0 | 46 | 0.0143 |
| V-GT-02 | temperature | Success | 1.0e-8 | 100 | 0.0 | 3 | 1.0 | 0.0 |  | true | 0.0 | 52 | 0.0153 |
| V-HC-01 | temperature | Success | 1.0e-8 | 100 | 3.03e-7 | 5 | 180.0 | 0.0395 |  | true | 0.0 | 60 | 0.00843 |
| V-HC-02a | temperature | Success | 1.0e-8 | 100 | 4.51e-8 | 5 | 10000.0 | 0.0143 |  | true | 0.0 | 56 | 0.00808 |
| V-HC-02b | temperature | Success | 1.0e-8 | 100 | 9.55e-8 | 5 | 10000.0 | 0.024 |  | true | 0.0 | 75 | 0.00781 |
| V-HC-02c | temperature | Success | 1.0e-8 | 100 | 1.24e-7 | 5 | 10000.0 | 0.0256 |  | true | 0.0 | 98 | 0.0106 |
| V-HC-02d | temperature | Success | 1.0e-8 | 100 | 1.28e-7 | 5 | 10000.0 | 0.0258 |  | true | 0.0 | 113 | 0.0123 |
| V-HC-03 | temperature | Success | 1.0e-8 | 100 | 2.43e-8 | 4 | 1200.0 | 0.00366 |  | true | 0.0 | 86 | 0.00678 |
| V-HC-04a | temperature | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 |  | true | 0.0 | 99 | 0.00812 |
| V-HC-04b | temperature | Success | 1.0e-8 | 100 | 3.98e-13 | 3 | 400000.0 | 1.0e-7 |  | true | 0.0 | 97 | 0.00831 |
| V-HC-04c | temperature | Success | 1.0e-8 | 100 | 1.71e-13 | 3 | 40000.0 | 4.29e-8 |  | true | 0.0 | 99 | 0.00857 |
| V-HC-05 | temperature | Success | 1.0e-8 | 100 | 1.14e-13 | 3 | 5000.0 | 1.43e-8 |  | true | 0.0 | 243 | 0.0243 |
| V-MMS-01 | temperature | Success | 1.0e-8 | 100 | 2.31e-8 | 3 | 20.0 | 0.00342 |  | true | 0.0 | 127 | 0.204 |
| V-MMS-02 | temperature | Success | 1.0e-8 | 100 | 3.98e-13 | 41 | 10.0 | 7.96e-8 |  | true | 0.00178 | 26 | 0.00601 |
| V-MMS-03 | temperature | Success | 1.0e-8 | 100 | 1.05e-5 | 41 | 10.0 | 1.69 |  | true | 1.58e-6 | 80 | 0.0427 |
| V-MMS-04 | temperature | Success | 1.0e-8 | 100 | 8.67e-9 | 41 | 20.0 | 0.00141 |  | true | 9.1e-9 | 63 | 0.133 |
| V-RAD-01 | temperature | Success | 1.0e-8 | 100 | 2.08e-7 | 5 | 180.0 | 0.0286 |  | true | 0.0 | 63 | 0.00879 |
| V-RX-01 | temperature | Success | 1.0e-10 | 100 | 0.0 | 7 | 20.0 | 0.0 |  | true | 0.0 | 175 | 0.00132 |
| V-RX-02a | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 |  | true | 0.0 | 13 | 0.000583 |
| V-RX-02b | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 |  | true | 0.0 | 242 | 0.000993 |
| V-RX-02c | temperature | Success | 1.0e-10 | 100 | 0.0 | 7 | 7.17 | 0.0 |  | true | 0.0 | 183 | 0.00106 |
| V-RX-02d | temperature | Success | 1.0e-10 | 100 | 5.68e-14 | 21 | 40.0 | 1.42e-6 |  | true | 0.0 | 306 | 0.00129 |
| V-RX-03 | temperature | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 15.0 | 2.84e-6 |  | true | 0.0 | 160 | 0.00152 |
| V-RX-04a | temperature | Success | 1.0e-8 | 100 | 3.02e-7 | 6 | 8.0 | 0.0336 |  | true | 0.0 | 78 | 0.00117 |
| V-RX-04b | temperature | Success | 1.0e-8 | 100 | 5.97e-9 | 3 | 2000.0 | 0.000814 |  | true | 0.0 | 159 | 0.00346 |
| V-BO-02 | enthalpy | Success | 1.0e-6 | 100 | 23.1 | 964 | 500.0 | 25300.0 | 82900.0 | false | 0.0 | 2456 | 7.79 |
| V-CPL-01 | enthalpy | Success | 1.0e-6 | 100 | 0.000165 | 31 | 300.0 | 0.376 |  | true | 0.0 | 391 | 0.381 |
| V-CPL-01b | enthalpy | Success | 1.0e-6 | 100 | 1.43e-5 | 31 | 150.0 | 0.0296 |  | true | 0.0 | 1144 | 0.637 |
| V-CPL-02 | enthalpy | Success | 1.0e-5 | 100 | 0.0359 | 41 | 2.0 | 4.86 | 6580.0 | false | 0.422 | 1901 | 3.08 |
| V-CPL-02b | enthalpy | Success | 1.0e-5 | 100 | 0.614 | 41 | 2.0 | 79.5 | 5390.0 | false | 0.0791 | 1114 | 1.04 |
| V-CPL-03 | enthalpy | Success | 1.0e-8 | 100 | 5.12e-13 | 3 | 8.0 | 1.28e-7 |  | true | 1.61e-9 | 89 | 0.00418 |
| V-CPL-04 | enthalpy | Success | 1.0e-6 | 100 | 6.27e-7 | 11 | 2500.0 | 0.00124 |  | true | 0.0 | 31 | 0.00411 |
| V-CPL-05 | enthalpy | Success | 1.0e-8 | 10 | 2.86e-6 | 72 | 720.0 | 0.434 |  | true | 0.0 | 455 | 0.0458 |
| V-GT-01 | enthalpy | Success | 1.0e-10 | 100 | 2.07e-11 | 4 | 8.0 | 0.000517 |  | true | 0.0 | 47 | 0.0263 |
| V-GT-02 | enthalpy | Success | 1.0e-8 | 100 | 7.39e-13 | 3 | 1.0 | 1.48e-7 |  | true | 0.0 | 52 | 0.0121 |
| V-HC-01 | enthalpy | Success | 1.0e-8 | 100 | 4.98e-9 | 5 | 180.0 | 0.000651 |  | true | 0.0 | 294 | 0.182 |
| V-HC-02a | enthalpy | Success | 1.0e-8 | 100 | 1.69e-8 | 5 | 10000.0 | 0.00536 |  | true | 0.0 | 665 | 0.364 |
| V-HC-02b | enthalpy | Success | 1.0e-8 | 100 | 2.22e-8 | 5 | 10000.0 | 0.00557 |  | true | 0.0 | 547 | 0.27 |
| V-HC-02c | enthalpy | Success | 1.0e-8 | 100 | 1.47e-8 | 5 | 10000.0 | 0.00302 |  | true | 0.0 | 487 | 0.222 |
| V-HC-02d | enthalpy | Success | 1.0e-8 | 100 | 3.17e-8 | 5 | 10000.0 | 0.00637 |  | true | 0.0 | 468 | 0.205 |
| V-HC-03 | enthalpy | Success | 1.0e-8 | 100 | 5.83e-9 | 4 | 1200.0 | 0.000876 |  | true | 0.0 | 511 | 0.246 |
| V-HC-04a | enthalpy | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 |  | true | 0.0 | 504 | 0.138 |
| V-HC-04b | enthalpy | Success | 1.0e-8 | 100 | 1.51e-8 | 3 | 400000.0 | 0.0038 |  | true | 0.0 | 655 | 0.19 |
| V-HC-04c | enthalpy | Success | 1.0e-8 | 100 | 1.06e-7 | 3 | 40000.0 | 0.0267 |  | true | 0.0 | 504 | 0.137 |
| V-HC-05 | enthalpy | Success | 1.0e-8 | 100 | 3.8e-7 | 3 | 5000.0 | 0.0477 |  | true | 0.0 | 518 | 0.191 |
| V-MMS-01 | enthalpy | Success | 1.0e-8 | 100 | 1.32e-8 | 3 | 20.0 | 0.00195 |  | true | 0.0 | 126 | 0.324 |
| V-MMS-02 | enthalpy | Success | 1.0e-8 | 100 | 7.39e-13 | 41 | 10.0 | 1.48e-7 |  | true | 0.00196 | 27 | 0.00664 |
| V-MMS-03 | enthalpy | Success | 1.0e-8 | 100 | 3.46e-6 | 41 | 10.0 | 0.558 |  | true | 5.81e-7 | 91 | 0.166 |
| V-MMS-04 | enthalpy | Success | 1.0e-8 | 100 | 1.37e-8 | 41 | 20.0 | 0.00223 |  | true | 1.94e-7 | 67 | 0.232 |
| V-RAD-01 | enthalpy | Success | 1.0e-8 | 100 | 4.56e-8 | 5 | 180.0 | 0.00627 |  | true | 0.0 | 367 | 0.207 |
| V-RX-01 | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 20.0 | 2.84e-6 |  | true | 0.0 | 174 | 0.00139 |
| V-RX-02a | enthalpy | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 18.0 | 1.42e-6 |  | true | 0.0 | 14 | 0.00056 |
| V-RX-02b | enthalpy | Success | 1.0e-10 | 100 | 5.68e-14 | 6 | 100.0 | 1.42e-6 |  | true | 0.0 | 243 | 0.00122 |
| V-RX-02c | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 7.17 | 2.53e-6 |  | true | 0.0 | 184 | 0.00106 |
| V-RX-02d | enthalpy | Success | 1.0e-10 | 100 | 2.27e-13 | 21 | 40.0 | 5.68e-6 |  | true | 0.0 | 313 | 0.00178 |
| V-RX-03 | enthalpy | Success | 1.0e-10 | 100 | 1.14e-13 | 7 | 15.0 | 2.84e-6 |  | true | 0.0 | 167 | 0.00203 |
| V-RX-04a | enthalpy | Success | 1.0e-8 | 100 | 2.16e-7 | 6 | 8.0 | 0.024 |  | true | 0.0 | 87 | 0.00175 |
| V-RX-04b | enthalpy | Success | 1.0e-8 | 100 | 9.46e-9 | 3 | 2000.0 | 0.00129 |  | true | 0.0 | 297 | 0.0251 |

60 of 66 solves within budget; 66 of 66 with every constraint residual at or below one.
