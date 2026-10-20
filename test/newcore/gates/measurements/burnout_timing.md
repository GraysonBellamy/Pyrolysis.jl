# Gate item 8 at M2: the burnout path against the old core

Produced by `test/newcore/gates/burnout_timing.jl` at commit 616bcd4 with uncommitted changes. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [10.3, 11.3, 11.5] at the start and [10.3, 10.4, 10.8] at the end. 60 PMMA cells; V-BO-01 on the finite-volume core with the prototype's correction (τ 0.2 s per cell, band 10 ε_m). Window error: the largest temperature difference [K] over the cells at the saved times before the first event of any solve of the deck, against the same core at reltol 1.0e-9; run error: over the whole run at the saved times with the same live cells; shift: the largest event-time difference over the reference's median interval; MLR: the L1 error of the mass loss rate history over the reference's integral. Warm medians of 3 solves (a solve slower than 15 times the old default is timed once and its configuration not taken to tighter tolerances); ratio: over the old core's default.

## V-BO-02

Window: the saved times before 62.4 s. Reference spread (reltol 1.0e-10 against 1.0e-9): old 0.0309 K, new 3.34e-5 K. Target: the old core's window error at reltol 1.0e-6, 0.375 K.

| run | reltol | retcode | window error [K] | run error [K] | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-5 | Success | 0.438 | 50.3 | 0.266 | 0.0141 | 24.539 | 452.0 | 0.00126 | 0.815 | 0.811 | 1276 | 0 | 32497 | 293 | 47 |
| old core default (KenCarp4) | 3.0e-6 | Success | 0.418 | 46.9 | 0.239 | 0.013 | 24.55 | 452.0 | 0.000217 | 0.942 | 0.937 | 1434 | 0 | 36411 | 329 | 57 |
| old core default (KenCarp4) | 1.0e-6 | Success | 0.375 | 30.4 | 0.185 | 0.0101 | 24.623 | 452.0 | 2.3e-5 | 1.0 | 1.0 | 1572 | 0 | 38574 | 341 | 56 |
| old core default (KenCarp4) | 3.0e-7 | Success | 0.375 | 44.1 | 0.21 | 0.0118 | 24.519 | 452.0 | 2.13e-5 | 1.08 | 1.08 | 1761 | 11 | 41739 | 345 | 11 |
| old core default (KenCarp4) | 1.0e-7 | Success | 0.199 | 39.7 | 0.223 | 0.0111 | 24.512 | 452.0 | 0.000105 | 1.23 | 1.23 | 2029 | 44 | 48530 | 357 | 4 |
| Rodas5P (default) | 1.0e-5 | Success | 0.531 | 17.2 | 0.149 | 0.00377 | 23.722 | 465.0 | 0.00535 | 0.521 | 0.518 | 818 | 276 | 9570 | 818 | 0 |
| Rodas5P (default) | 3.0e-6 | Success | 0.206 | 24.4 | 0.134 | 0.00897 | 23.718 | 465.5 | 0.00444 | 0.783 | 0.78 | 1183 | 542 | 14983 | 1183 | 0 |
| Rodas5P (default) | 1.0e-6 | Success | 0.065 | 29.1 | 0.147 | 0.00981 | 23.666 | 465.5 | 0.0015 | 1.18 | 1.18 | 1706 | 953 | 22978 | 1706 | 0 |
| Rodas5P (default) | 3.0e-7 | Success | 0.0146 | 20.8 | 0.114 | 0.00758 | 23.688 | 465.5 | 0.000955 | 1.77 | 1.76 | 2494 | 1619 | 35398 | 2494 | 0 |
| Rodas5P (default) | 1.0e-7 | Success | 0.00445 | 16.9 | 0.0677 | 0.00422 | 23.73 | 465.0 | 9.13e-5 | 2.84 | 2.82 | 3788 | 3010 | 58172 | 3788 | 0 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-5 | Success | 0.365 | 12.1 | 0.0734 | 0.00277 | 23.728 | 465.0 | 0.0238 | 7.36 | 7.33 | 18551 | 0 | 215217 | 1584 | 4255 |
| KenCarp4 raw (kencarp4_policy) | 3.0e-6 | Success | 0.308 | 1.09 | 0.032 | 0.000717 | 23.733 | 465.0 | 0.00177 | 7.22 | 7.19 | 17161 | 3 | 220267 | 1583 | 3157 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-6 | Success | 0.255 | 6.01 | 0.0878 | 0.00168 | 23.724 | 465.0 | 0.0252 | 5.11 | 5.09 | 10021 | 43 | 158084 | 1025 | 1784 |
| KenCarp4 raw (kencarp4_policy) | 3.0e-7 | Success | 0.279 | 10.4 | 0.0698 | 0.00206 | 23.722 | 465.0 | 0.0247 | 3.02 | 3.01 | 4695 | 262 | 93485 | 530 | 980 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-7 | Success | 0.0559 | 19.9 | 0.0765 | 0.00507 | 23.73 | 465.0 | 0.0124 | 2.83 | 2.82 | 3923 | 510 | 87316 | 417 | 843 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-5 | Success | 0.448 | 42.3 | 0.27 | 0.016 | 23.724 | 466.0 | 5.76e-5 | 0.724 | 0.72 | 1052 | 0 | 18186 | 377 | 51 |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-6 | Success | 0.432 | 47.8 | 0.318 | 0.0211 | 23.651 | 466.5 | 1.72e-5 | 0.769 | 0.766 | 1122 | 0 | 19394 | 417 | 39 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-6 | Success | 0.374 | 48.1 | 0.267 | 0.0173 | 23.666 | 466.0 | 0.000606 | 0.858 | 0.854 | 1239 | 0 | 21904 | 405 | 40 |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-7 | Success | 0.357 | 30.9 | 0.195 | 0.0137 | 23.669 | 466.0 | 8.48e-5 | 1.0 | 0.997 | 1385 | 7 | 26778 | 400 | 47 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-7 | Success | 0.24 | 24.5 | 0.125 | 0.00826 | 23.697 | 465.5 | 4.58e-5 | 1.2 | 1.19 | 1545 | 54 | 32631 | 432 | 49 |

Matched accuracy (the loosest ladder tolerance at or below the target):

| run | reltol | window error [K] | warm [s] | ratio |
|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.375 | 1.0 | 1.0 |
| Rodas5P (default) | 3.0e-6 | 0.206 | 0.783 | 0.78 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-5 | 0.365 | 7.36 | 7.33 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-6 | 0.374 | 0.858 | 0.854 |

Where the time goes at the matched tolerance (shares of the profiler's samples) and unit costs at 200.0 s [µs]:

| run | reltol | jacobian | domain | events | factorization | linsolve | rhs | wform | observers | integrator | RHS | Jacobian | refactor + solve | solve | states | nnz(W) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.047 | 0.0 | 0.00049 | 0.06 | 0.31 | 0.46 | 0.022 | 0.0017 | 0.095 | 13.4 | 132.0 | 20.5 | 6.48 | 728 | 3655 |
| Rodas5P (default) | 3.0e-6 | 0.23 | 0.042 | 0.0027 | 0.071 | 0.15 | 0.41 | 0.026 | 0.0 | 0.074 | 19.3 | 145.0 | 16.7 | 5.42 | 728 | 2794 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-5 | 0.044 | 0.049 | 0.00089 | 0.04 | 0.23 | 0.56 | 0.017 | 3.4e-5 | 0.065 | 20.0 | 148.0 | 17.7 | 5.97 | 728 | 2806 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-6 | 0.089 | 0.029 | 0.0029 | 0.056 | 0.23 | 0.51 | 0.022 | 0.0 | 0.068 | 19.1 | 146.0 | 16.7 | 5.73 | 728 | 2810 |

## V-BO-01

Window: the saved times before 167.8 s. Reference spread (reltol 1.0e-10 against 1.0e-9): old 0.0158 K, new 7.0e-5 K. Target: the old core's window error at reltol 1.0e-6, 0.263 K.

| run | reltol | retcode | window error [K] | run error [K] | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-5 | Success | 0.247 | 0.547 | 0.539 | 0.000271 | 32.392 | 363.5 | 0.00353 | 0.994 | 1.77 | 1304 | 0 | 31398 | 232 | 120 |
| old core default (KenCarp4) | 3.0e-6 | Success | 0.271 | 1.65 | 0.589 | 0.000343 | 32.392 | 363.5 | 0.00106 | 0.675 | 1.2 | 890 | 0 | 20802 | 185 | 37 |
| old core default (KenCarp4) | 1.0e-6 | Success | 0.263 | 0.807 | 1.08 | 0.000296 | 32.392 | 363.5 | 0.000478 | 0.561 | 1.0 | 853 | 3 | 19076 | 164 | 8 |
| old core default (KenCarp4) | 3.0e-7 | Success | 0.206 | 1.79 | 1.06 | 6.47e-5 | 32.392 | 363.5 | 2.96e-5 | 0.632 | 1.13 | 979 | 10 | 21951 | 151 | 5 |
| old core default (KenCarp4) | 1.0e-7 | Success | 0.247 | 1.95 | 0.644 | 3.18e-5 | 32.392 | 363.5 | 0.000124 | 0.754 | 1.35 | 1141 | 38 | 26199 | 158 | 5 |
| Rodas5P (default) | 1.0e-5 | Success | 0.453 | 6.17 |  | 0.000389 | 32.391 | 363.5 | 3.04e-5 | 1.24 | 2.21 | 733 | 187 | 8093 | 733 | 0 |
| Rodas5P (default) | 3.0e-6 | Success | 0.305 | 5.76 |  | 0.000849 | 32.391 | 363.5 | 0.000161 | 1.56 | 2.77 | 928 | 248 | 10336 | 928 | 0 |
| Rodas5P (default) | 1.0e-6 | Success | 0.0911 | 10.5 |  | 0.000121 | 32.391 | 363.5 | 6.87e-5 | 1.74 | 3.1 | 1090 | 245 | 11770 | 1090 | 0 |
| Rodas5P (default) | 3.0e-7 | Success | 0.0225 | 2.05 |  | 0.000263 | 32.391 | 363.5 | 7.08e-5 | 2.68 | 4.78 | 1581 | 571 | 18797 | 1581 | 0 |
| Rodas5P (default) | 1.0e-7 | Success | 0.00956 | 1.82 |  | 1.42e-5 | 32.391 | 363.5 | 3.96e-5 | 3.44 | 6.14 | 1999 | 831 | 24639 | 1999 | 0 |
| KenCarp4 raw (kencarp4_policy) | 1.0e-5 | Success | 0.461 | 12.3 |  | 0.00015 | 32.391 | 363.5 | 0.000193 | 127.0 | 226.0 | 130225 | 6 | 1678590 | 13994 | 36009 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-5 | Success | 0.467 | 3.11 |  | 0.000316 | 32.391 | 363.5 | 3.94e-5 | 0.834 | 1.49 | 464 | 15 | 10899 | 155 | 22 |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-6 | Success | 0.297 | 4.29 |  | 0.000333 | 32.391 | 363.5 | 1.8e-6 | 1.04 | 1.85 | 557 | 56 | 14206 | 164 | 19 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-6 | Success | 0.296 | 9.4 |  | 0.000166 | 32.391 | 363.5 | 7.33e-5 | 1.26 | 2.25 | 691 | 83 | 18397 | 152 | 16 |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-7 | Success | 0.183 | 1.57 |  | 0.000259 | 32.391 | 363.5 | 8.55e-5 | 1.62 | 2.89 | 886 | 166 | 24468 | 153 | 5 |
| KenCarp4 smoothed, forced restart, flag cleared | 1.0e-7 | Success | 0.197 | 4.23 |  | 0.000323 | 32.391 | 363.5 | 8.55e-6 | 2.14 | 3.81 | 1267 | 227 | 34006 | 141 | 3 |

Matched accuracy (the loosest ladder tolerance at or below the target):

| run | reltol | window error [K] | warm [s] | ratio |
|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.263 | 0.561 | 1.0 |
| Rodas5P (default) | 1.0e-6 | 0.0911 | 1.74 | 3.1 |
| KenCarp4 raw (kencarp4_policy) | none on the ladder | | | |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-7 | 0.183 | 1.62 | 2.89 |

Where the time goes at the matched tolerance (shares of the profiler's samples) and unit costs at 200.0 s [µs]:

| run | reltol | jacobian | domain | events | factorization | linsolve | rhs | wform | observers | integrator | RHS | Jacobian | refactor + solve | solve | states | nnz(W) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.063 | 0.0 | 0.0012 | 0.046 | 0.27 | 0.45 | 0.031 | 0.0058 | 0.13 | 18.1 | 176.0 | 20.9 | 7.34 | 911 | 3891 |
| Rodas5P (default) | 1.0e-6 | 0.36 | 0.023 | 0.00027 | 0.12 | 0.13 | 0.22 | 0.057 | 0.00014 | 0.097 | 27.7 | 518.0 | 40.7 | 10.9 | 1277 | 6731 |
| KenCarp4 smoothed, forced restart, flag cleared | 3.0e-7 | 0.061 | 0.019 | 0.00043 | 0.058 | 0.31 | 0.44 | 0.026 | 0.0 | 0.091 | 27.1 | 496.0 | 38.0 | 10.3 | 1277 | 6706 |

