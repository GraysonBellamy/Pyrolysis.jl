# Gate item 2: where the time goes at matched accuracy

Produced by `test/newcore/gates/audit_breakdown.jl` at commit 603e863. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [11.3, 11.9, 12.0] at the start and [13.4, 12.4, 12.2] at the end. Each run at the cheapest tolerance of its ladder whose error is at or below the old default's (error against the same core at reltol 1.0e-11). Unit costs of the new core at the initial state [µs]: one residual, one colored Jacobian, one KLU refactorization of the iteration matrix with a solve, one solve; the old core's row repeats them. Parts: Jacobians × Jacobian, right-hand sides × residual, factorizations × factorization, solves × solve, each as a share of the warm time; the remainder is the rest.

| deck | run | reltol | error [K] | warm [s] | steps | rejected | RHS | Jacobians | factorizations | solves | residual [µs] | Jacobian [µs] | factorization [µs] | solve [µs] | Jacobians share | RHS share | factorizations share | solves share | rest |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 100, solid λ 1.0e-6 | old core default | 1.0e-6 | 1.37e-5 | 0.0139 | 56 | 3 | 1362 | 12 | 31 | 1365 | 5.48 | 43.8 | 12.8 | 2.44 |  |  |  |  |  |
| eulerian 100, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 9.77e-6 | 0.019 | 98 | 4 | 914 | 98 | 102 | 816 | 5.48 | 43.8 | 12.8 | 2.44 | 0.23 | 0.26 | 0.069 | 0.1 | 0.34 |
| eulerian 100, solid λ 1.0e-6 | new KenCarp4 raw | 1.0e-7 | 1.19e-5 | 0.0203 | 74 | 0 | 1651 | 13 | 37 | 1577 | 5.48 | 43.8 | 12.8 | 2.44 | 0.028 | 0.45 | 0.023 | 0.19 | 0.31 |
| eulerian 100, solid λ 1.0e-6 | new KenCarp4 smoothed | NaN |  |  | 0 | 0 | 0 | 0 | 0 | 0 | 5.48 | 43.8 | 12.8 | 2.44 |  |  |  |  |  |
| ale 25, solid λ 1.0e-6 | old core default | 1.0e-6 | 2.31e-5 | 0.00505 | 51 | 2 | 1337 | 12 | 25 | 1339 | 1.94 | 20.8 | 4.62 | 1.11 |  |  |  |  |  |
| ale 25, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 3.09e-6 | 0.00751 | 79 | 4 | 743 | 79 | 83 | 664 | 1.94 | 20.8 | 4.62 | 1.11 | 0.22 | 0.19 | 0.051 | 0.098 | 0.44 |
| ale 25, solid λ 1.0e-6 | new KenCarp4 raw | 3.0e-8 | 7.04e-6 | 0.00876 | 76 | 1 | 1783 | 11 | 29 | 1707 | 1.94 | 20.8 | 4.62 | 1.11 | 0.026 | 0.4 | 0.015 | 0.22 | 0.35 |
| ale 25, solid λ 1.0e-6 | new KenCarp4 smoothed | 3.0e-8 | 2.26e-5 | 0.00693 | 51 | 0 | 1222 | 15 | 26 | 1222 | 1.94 | 20.8 | 4.62 | 1.11 | 0.045 | 0.34 | 0.017 | 0.2 | 0.4 |
| ale 100, solid λ 1.0e-6 | old core default | 1.0e-6 | 1.03e-5 | 0.0155 | 54 | 3 | 1283 | 15 | 34 | 1286 | 6.93 | 78.1 | 18.0 | 4.12 |  |  |  |  |  |
| ale 100, solid λ 1.0e-6 | new Rodas5P (default) | 1.0e-7 | 6.54e-6 | 0.0286 | 90 | 3 | 834 | 90 | 93 | 744 | 6.93 | 78.1 | 18.0 | 4.12 | 0.25 | 0.2 | 0.058 | 0.11 | 0.39 |
| ale 100, solid λ 1.0e-6 | new KenCarp4 raw | 3.0e-8 | 5.69e-6 | 0.0323 | 89 | 1 | 2047 | 14 | 39 | 1958 | 6.93 | 78.1 | 18.0 | 4.12 | 0.034 | 0.44 | 0.022 | 0.25 | 0.26 |
| ale 100, solid λ 1.0e-6 | new KenCarp4 smoothed | NaN |  |  | 0 | 0 | 0 | 0 | 0 | 0 | 6.93 | 78.1 | 18.0 | 4.12 |  |  |  |  |  |
| ale_darcy 25 | old core default | 1.0e-6 | 0.000126 | 0.011 | 75 | 2 | 2428 | 49 | 62 | 2430 | 2.88 | 25.4 | 5.78 | 1.27 |  |  |  |  |  |
| ale_darcy 25 | new Rodas5P (default) | 1.0e-6 | 7.28e-5 | 0.0121 | 100 | 3 | 924 | 100 | 103 | 824 | 2.88 | 25.4 | 5.78 | 1.27 | 0.21 | 0.22 | 0.049 | 0.087 | 0.43 |
| ale_darcy 25 | new KenCarp4 raw | 1.0e-5 | 7.43e-5 | 0.0476 | 341 | 0 | 7172 | 71 | 140 | 6831 | 2.88 | 25.4 | 5.78 | 1.27 | 0.038 | 0.43 | 0.017 | 0.18 | 0.33 |
| ale_darcy 25 | new KenCarp4 smoothed | 3.0e-6 | 5.74e-5 | 0.0335 | 191 | 0 | 4904 | 51 | 96 | 4904 | 2.88 | 25.4 | 5.78 | 1.27 | 0.039 | 0.42 | 0.017 | 0.19 | 0.34 |
| ale_darcy 100 | old core default | 1.0e-6 | 5.96e-5 | 0.0424 | 80 | 2 | 2726 | 51 | 69 | 2728 | 10.3 | 106.0 | 22.4 | 5.24 |  |  |  |  |  |
| ale_darcy 100 | new Rodas5P (default) | 3.0e-7 | 1.44e-5 | 0.0593 | 140 | 3 | 1284 | 140 | 143 | 1144 | 10.3 | 106.0 | 22.4 | 5.24 | 0.25 | 0.22 | 0.054 | 0.1 | 0.37 |
| ale_darcy 100 | new KenCarp4 raw | 3.0e-7 | 1.66e-5 | 0.183 | 270 | 1 | 8076 | 63 | 126 | 7806 | 10.3 | 106.0 | 22.4 | 5.24 | 0.036 | 0.46 | 0.015 | 0.22 | 0.27 |
| ale_darcy 100 | new KenCarp4 smoothed | 3.0e-7 | 3.22e-5 | 0.186 | 268 | 3 | 7969 | 64 | 124 | 7972 | 10.3 | 106.0 | 22.4 | 5.24 | 0.036 | 0.44 | 0.015 | 0.22 | 0.28 |
