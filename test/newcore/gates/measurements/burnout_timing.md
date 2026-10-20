# Gate item 8 at M2: the burnout path against the old core

Produced by `test/newcore/gates/burnout_timing.jl` at commit 3efdc60 with uncommitted changes. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [2.2, 1.6, 1.7] at the start and [2.4, 2.0, 1.8] at the end. 60 PMMA cells; V-BO-01 on the finite-volume core with the prototype's correction (τ 0.2 s per cell, band 10 ε_m). Window RMS and largest: the RMS and the largest of the temperature differences [K] over the cells at the saved times before the first event of any solve of the deck, against the same core at reltol 1.0e-9; the new core is matched on the RMS. Run error: the largest difference over the whole run at the saved times with the same live cells; shift: the largest event-time difference over the reference's median interval; MLR: the L1 error of the mass loss rate history over the reference's integral. Warm medians of 3 solves, of 9 for the old default (a solve slower than 15 times the old default is timed once and its configuration not taken to tighter tolerances); ratio: over the old core's default.

## V-BO-02

Window: the saved times before 62.4 s. Reference spread (reltol 1.0e-10 against 1.0e-9, window RMS): old 0.000638 K, new 1.01e-6 K. Target: the old core's window RMS at reltol 1.0e-6, 0.015 K.

| run | reltol | retcode | window RMS [K] | window largest [K] | run error [K] | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-5 | Success | 0.0221 | 0.438 | 50.3 | 0.266 | 0.0141 | 24.539 | 452.0 | 0.00126 | 0.827 | 0.839 | 1276 | 0 | 32497 | 293 | 47 |
| old core default (KenCarp4) | 3.0e-6 | Success | 0.0161 | 0.418 | 46.9 | 0.239 | 0.013 | 24.55 | 452.0 | 0.000217 | 0.924 | 0.937 | 1434 | 0 | 36411 | 329 | 57 |
| old core default (KenCarp4) | 1.0e-6 | Success | 0.015 | 0.375 | 30.4 | 0.185 | 0.0101 | 24.623 | 452.0 | 2.3e-5 | 0.986 | 1.0 | 1572 | 0 | 38574 | 341 | 56 |
| old core default (KenCarp4) | 3.0e-7 | Success | 0.0109 | 0.375 | 44.1 | 0.21 | 0.0118 | 24.519 | 452.0 | 2.13e-5 | 1.12 | 1.14 | 1761 | 11 | 41739 | 345 | 11 |
| old core default (KenCarp4) | 1.0e-7 | Success | 0.00561 | 0.199 | 39.7 | 0.223 | 0.0111 | 24.512 | 452.0 | 0.000105 | 1.28 | 1.3 | 2029 | 44 | 48530 | 357 | 4 |
| Rodas5P (default) | 1.0e-5 | Success | 0.0101 | 0.531 | 44.7 | 0.244 | 0.0157 | 23.621 | 466.0 | 2.39e-5 | 0.615 | 0.624 | 1360 | 129 | 13392 | 1360 | 0 |
| Rodas5P (default) | 3.0e-6 | Success | 0.00622 | 0.206 | 26.8 | 0.154 | 0.0101 | 23.709 | 465.5 | 5.76e-5 | 0.74 | 0.751 | 1636 | 177 | 16260 | 1636 | 0 |
| Rodas5P (default) | 1.0e-6 | Success | 0.00245 | 0.065 | 28.0 | 0.141 | 0.00992 | 23.71 | 465.5 | 4.05e-5 | 0.929 | 0.942 | 2014 | 350 | 21046 | 2014 | 0 |
| Rodas5P (default) | 3.0e-7 | Success | 0.000264 | 0.0145 | 18.2 | 0.102 | 0.00727 | 23.721 | 465.5 | 0.00011 | 1.3 | 1.32 | 2671 | 815 | 30679 | 2671 | 0 |
| Rodas5P (default) | 1.0e-7 | Success | 7.09e-5 | 0.00441 | 14.2 | 0.0501 | 0.00303 | 23.708 | 465.0 | 0.000115 | 2.09 | 2.12 | 3797 | 2261 | 52381 | 3797 | 0 |

Matched accuracy (the loosest ladder tolerance at or below the target):

| run | reltol | window RMS [K] | window largest [K] | warm [s] | ratio |
|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.015 | 0.375 | 0.986 | 1.0 |
| Rodas5P (default) | 1.0e-5 | 0.0101 | 0.531 | 0.615 | 0.624 |

Where the time goes at the matched tolerance (shares of the profiler's samples) and unit costs at 200.0 s [µs]:

| run | reltol | jacobian | domain | events | factorization | linsolve | rhs | wform | observers | integrator | RHS | Jacobian | states | nnz(J) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.053 | 0.0 | 0.0016 | 0.079 | 0.3 | 0.46 | 0.024 | 0.0042 | 0.081 | 13.3 | 118.0 | 728 | 8704 |
| Rodas5P (default) | 1.0e-5 | 0.32 | 0.029 | 0.0031 | 0.091 | 0.11 | 0.34 | 0.045 | 0.0 | 0.066 | 13.5 | 147.0 | 728 | 8704 |

## V-BO-01

Window: the saved times before 167.8 s. Reference spread (reltol 1.0e-10 against 1.0e-9, window RMS): old 0.000323 K, new 1.35e-6 K. Target: the old core's window RMS at reltol 1.0e-6, 0.00866 K.

| run | reltol | retcode | window RMS [K] | window largest [K] | run error [K] | shift | MLR | peak | t_bo [s] | ledger | warm [s] | ratio | steps | rejected | RHS | Jacobians | nonlinear failures |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-5 | Success | 0.0104 | 0.247 | 0.547 | 0.539 | 0.000271 | 32.392 | 363.5 | 0.00353 | 0.898 | 1.61 | 1304 | 0 | 31398 | 232 | 120 |
| old core default (KenCarp4) | 3.0e-6 | Success | 0.0112 | 0.271 | 1.65 | 0.589 | 0.000343 | 32.392 | 363.5 | 0.00106 | 0.604 | 1.08 | 890 | 0 | 20802 | 185 | 37 |
| old core default (KenCarp4) | 1.0e-6 | Success | 0.00866 | 0.263 | 0.807 | 1.08 | 0.000296 | 32.392 | 363.5 | 0.000478 | 0.559 | 1.0 | 853 | 3 | 19076 | 164 | 8 |
| old core default (KenCarp4) | 3.0e-7 | Success | 0.00648 | 0.206 | 1.79 | 1.06 | 6.47e-5 | 32.392 | 363.5 | 2.96e-5 | 0.633 | 1.13 | 979 | 10 | 21951 | 151 | 5 |
| old core default (KenCarp4) | 1.0e-7 | Success | 0.00602 | 0.247 | 1.95 | 0.644 | 3.18e-5 | 32.392 | 363.5 | 0.000124 | 0.75 | 1.34 | 1141 | 38 | 26199 | 158 | 5 |
| Rodas5P (default) | 1.0e-5 | Success | 0.00591 | 0.453 | 7.57 |  | 0.000477 | 32.391 | 363.5 | 4.53e-5 | 0.915 | 1.64 | 766 | 179 | 8329 | 767 | 0 |
| Rodas5P (default) | 3.0e-6 | Success | 0.00285 | 0.177 | 13.2 |  | 0.000557 | 32.391 | 363.5 | 7.83e-5 | 1.09 | 1.95 | 930 | 219 | 10125 | 931 | 0 |
| Rodas5P (default) | 1.0e-6 | Success | 0.0033 | 0.0912 | 11.1 |  | 0.000527 | 32.391 | 363.5 | 7.41e-5 | 1.28 | 2.29 | 1115 | 226 | 11845 | 1115 | 0 |
| Rodas5P (default) | 3.0e-7 | Success | 0.000499 | 0.0225 | 4.64 |  | 0.000403 | 32.391 | 363.5 | 1.18e-5 | 1.87 | 3.35 | 1581 | 564 | 18743 | 1581 | 0 |
| Rodas5P (default) | 1.0e-7 | Success | 0.000117 | 0.00886 | 5.55 |  | 6.07e-5 | 32.391 | 363.5 | 4.28e-5 | 2.4 | 4.29 | 2017 | 799 | 24547 | 2017 | 0 |

Matched accuracy (the loosest ladder tolerance at or below the target):

| run | reltol | window RMS [K] | window largest [K] | warm [s] | ratio |
|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.00866 | 0.263 | 0.559 | 1.0 |
| Rodas5P (default) | 1.0e-5 | 0.00591 | 0.453 | 0.915 | 1.64 |

Where the time goes at the matched tolerance (shares of the profiler's samples) and unit costs at 200.0 s [µs]:

| run | reltol | jacobian | domain | events | factorization | linsolve | rhs | wform | observers | integrator | RHS | Jacobian | states | nnz(J) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| old core default (KenCarp4) | 1.0e-6 | 0.055 | 0.0 | 0.0016 | 0.071 | 0.29 | 0.45 | 0.025 | 0.0059 | 0.1 | 15.1 | 157.0 | 911 | 8887 |
| Rodas5P (default) | 1.0e-5 | 0.39 | 0.024 | 0.0 | 0.092 | 0.068 | 0.22 | 0.075 | 0.00022 | 0.12 | 22.7 | 470.0 | 1277 | 27904 |

