# Gate item 2: where the time goes at matched accuracy

Produced by `test/newcore/gates/audit_breakdown.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [23.8, 21.4, 20.8] at the start and [23.6, 22.3, 21.2] at the end. Each run at the cheapest tolerance of its ladder whose error is at or below the old default's (error against the same core at reltol 1.0e-11). Unit costs of the new core at the initial state [µs]: one residual, one colored Jacobian, one KLU refactorization of the iteration matrix with a solve, one solve; the old core's row repeats them. Parts: Jacobians × Jacobian, right-hand sides × residual, factorizations × factorization, solves × solve, each as a share of the warm time; the remainder is the rest.

| deck | run | reltol | error [K] | warm [s] | steps | rejected | RHS | Jacobians | factorizations | solves | residual [µs] | Jacobian [µs] | factorization [µs] | solve [µs] | Jacobians share | RHS share | factorizations share | solves share | rest |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 100, solid λ 1.0e-6 | old core default | 1.0e-6 | 1.37e-5 | 0.0245 | 56 | 3 | 1362 | 12 | 31 | 1365 | 10.4 | 79.0 | 14.8 | 2.99 |  |  |  |  |  |
| eulerian 100, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 9.77e-6 | 0.0329 | 98 | 4 | 914 | 98 | 102 | 816 | 10.4 | 79.0 | 14.8 | 2.99 | 0.24 | 0.29 | 0.046 | 0.074 | 0.36 |
| eulerian 100, solid λ 1.0e-6 | new KenCarp4 raw | 1.0e-7 | 1.19e-5 | 0.061 | 74 | 0 | 1651 | 13 | 37 | 1577 | 10.4 | 79.0 | 14.8 | 2.99 | 0.017 | 0.28 | 0.009 | 0.077 | 0.62 |
| eulerian 100, solid λ 1.0e-6 | new KenCarp4 smoothed | NaN |  |  | 0 | 0 | 0 | 0 | 0 | 0 | 10.4 | 79.0 | 14.8 | 2.99 |  |  |  |  |  |
| ale 25, solid λ 1.0e-6 | old core default | 1.0e-6 | 2.31e-5 | 0.0138 | 51 | 2 | 1337 | 12 | 25 | 1339 | 6.77 | 52.5 | 12.1 | 2.36 |  |  |  |  |  |
| ale 25, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 3.09e-6 | 0.0314 | 79 | 4 | 743 | 79 | 83 | 664 | 6.77 | 52.5 | 12.1 | 2.36 | 0.13 | 0.16 | 0.032 | 0.05 | 0.63 |
| ale 25, solid λ 1.0e-6 | new KenCarp4 raw | 3.0e-8 | 7.04e-6 | 0.0264 | 76 | 1 | 1783 | 11 | 29 | 1707 | 6.77 | 52.5 | 12.1 | 2.36 | 0.022 | 0.46 | 0.013 | 0.15 | 0.36 |
| ale 25, solid λ 1.0e-6 | new KenCarp4 smoothed | 3.0e-8 | 2.26e-5 | 0.0185 | 51 | 0 | 1222 | 15 | 26 | 1222 | 6.77 | 52.5 | 12.1 | 2.36 | 0.043 | 0.45 | 0.017 | 0.16 | 0.34 |
| ale 100, solid λ 1.0e-6 | old core default | 1.0e-6 | 1.03e-5 | 0.0446 | 54 | 3 | 1283 | 15 | 34 | 1286 | 25.2 | 249.0 | 58.5 | 9.6 |  |  |  |  |  |
| ale 100, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 6.54e-6 | 0.124 | 90 | 3 | 834 | 90 | 93 | 744 | 25.2 | 249.0 | 58.5 | 9.6 | 0.18 | 0.17 | 0.044 | 0.057 | 0.55 |
| ale 100, solid λ 1.0e-6 | new KenCarp4 raw | 3.0e-8 | 5.69e-6 | 0.106 | 89 | 1 | 2047 | 14 | 39 | 1958 | 25.2 | 249.0 | 58.5 | 9.6 | 0.033 | 0.49 | 0.022 | 0.18 | 0.28 |
| ale 100, solid λ 1.0e-6 | new KenCarp4 smoothed | NaN |  |  | 0 | 0 | 0 | 0 | 0 | 0 | 25.2 | 249.0 | 58.5 | 9.6 |  |  |  |  |  |
| ale_darcy 25 | old core default | 1.0e-6 | 0.000126 | 0.0359 | 75 | 2 | 2428 | 49 | 62 | 2430 | 9.8 | 108.0 | 16.1 | 3.5 |  |  |  |  |  |
| ale_darcy 25 | new Rodas5P (default) | 1.0e-6 | 6.77e-5 | 0.0586 | 100 | 3 | 924 | 100 | 103 | 824 | 9.8 | 108.0 | 16.1 | 3.5 | 0.18 | 0.15 | 0.028 | 0.049 | 0.58 |
| ale_darcy 25 | new KenCarp4 raw | 1.0e-5 | 7.43e-5 | 0.156 | 341 | 0 | 7171 | 71 | 140 | 6830 | 9.8 | 108.0 | 16.1 | 3.5 | 0.049 | 0.45 | 0.014 | 0.15 | 0.33 |
| ale_darcy 25 | new KenCarp4 smoothed | 3.0e-6 | 5.74e-5 | 0.114 | 191 | 0 | 4904 | 51 | 96 | 4904 | 9.8 | 108.0 | 16.1 | 3.5 | 0.048 | 0.42 | 0.014 | 0.15 | 0.36 |
| ale_darcy 100 | old core default | 1.0e-6 | 5.96e-5 | 0.14 | 80 | 2 | 2726 | 51 | 69 | 2728 | 40.1 | 470.0 | 79.7 | 15.6 |  |  |  |  |  |
| ale_darcy 100 | new Rodas5P (default) | 3.0e-7 | 1.44e-5 | 0.294 | 140 | 3 | 1284 | 140 | 143 | 1144 | 40.1 | 470.0 | 79.7 | 15.6 | 0.22 | 0.18 | 0.039 | 0.061 | 0.5 |
| ale_darcy 100 | new KenCarp4 raw | 3.0e-7 | 1.66e-5 | 0.606 | 270 | 1 | 8076 | 63 | 126 | 7806 | 40.1 | 470.0 | 79.7 | 15.6 | 0.049 | 0.53 | 0.017 | 0.2 | 0.2 |
| ale_darcy 100 | new KenCarp4 smoothed | 3.0e-7 | 3.22e-5 | 0.596 | 268 | 3 | 7969 | 64 | 124 | 7972 | 40.1 | 470.0 | 79.7 | 15.6 | 0.05 | 0.54 | 0.017 | 0.21 | 0.19 |
