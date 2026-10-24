# Leftover 2: the integrator candidates on the finite-volume core

Produced by `test/newcore/gates/leftover2_integrators.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [18.0, 18.2, 18.3] at the start and [18.5, 18.4, 18.4] at the end. Error: the largest final-time temperature difference against KenCarp4 at abstol 1e-10, reltol 1e-8 [K] (V-CPL-02: the final surface temperature); spread: Rodas5P's at the same tight tolerances against that reference. Constraint: the largest scaled constraint residual of the accepted states (at most one passes). Warm medians of 3 solves. `qualifies` is a provisional rule for this table, not the decision (plan §10).

| deck | candidate | retcode | error | spread | constraint | steps | rejected | Jacobians | nonlinear failures | warm [s] | QoI failures | qualifies |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eulerian 25 | KenCarp4 | Success | 0.0339 | 2.86e-5 | 0.0 | 268 | 1 | 57 | 40 | 0.0352 |  | false |
| eulerian 25 | Rodas5P | Success | 0.0527 | 2.86e-5 | 0.0 | 239 | 45 | 241 | 0 | 0.0256 |  | false |
| eulerian 25 | Rodas5P forced restart | Success | 0.0527 | 2.86e-5 | 0.0 | 239 | 45 | 241 | 0 | 0.0254 |  | false |
| eulerian 25 | FBDF | DtLessThanMin |  | 2.86e-5 | 0.0 | 191 | 23 | 13 | 14 | 0.00756 |  | false |
| eulerian 25 | old core default | Success |  | 2.86e-5 |  | 292 | 24 | 73 | 14 | 0.0438 |  | false |
| eulerian 100 | KenCarp4 | Success | 0.0629 | 0.0227 | 0.0 | 866 | 3 | 147 | 141 | 0.323 |  | false |
| eulerian 100 | Rodas5P | Success | 0.0617 | 0.0227 | 0.0 | 989 | 235 | 992 | 0 | 0.372 |  | false |
| eulerian 100 | Rodas5P forced restart | Success | 0.0617 | 0.0227 | 0.0 | 989 | 235 | 992 | 0 | 0.373 |  | false |
| eulerian 100 | FBDF | DtLessThanMin |  | 0.0227 | 0.0 | 365 | 17 | 47 | 33 | 0.0519 |  | false |
| eulerian 100 | old core default | Success |  | 0.0227 |  | 941 | 18 | 200 | 62 | 0.447 |  | false |
| ale 25 | KenCarp4 | Success | 0.0148 | 0.0611 | 0.0188 | 281 | 3 | 60 | 40 | 0.054 |  | false |
| ale 25 | Rodas5P | Success | 0.0146 | 0.0611 | 3.69e-5 | 254 | 46 | 254 | 0 | 0.0608 |  | false |
| ale 25 | Rodas5P forced restart | Success | 0.0146 | 0.0611 | 3.69e-5 | 254 | 46 | 254 | 0 | 0.0613 |  | false |
| ale 25 | FBDF | DtLessThanMin |  | 0.0611 | 0.00124 | 258 | 21 | 18 | 20 | 0.0165 |  | false |
| ale 25 | old core default | Success |  | 0.0611 |  | 196 | 23 | 54 | 8 | 0.0355 |  | false |
| ale 100 | KenCarp4 | Success | 0.0635 | 0.0281 | 0.00695 | 1162 | 3 | 230 | 245 | 0.718 |  | false |
| ale 100 | Rodas5P | Success | 0.0507 | 0.0281 | 1.85e-6 | 1118 | 260 | 1124 | 0 | 0.993 |  | false |
| ale 100 | Rodas5P forced restart | Success | 0.0507 | 0.0281 | 1.85e-6 | 1118 | 260 | 1124 | 0 | 0.995 |  | false |
| ale 100 | FBDF | DtLessThanMin |  | 0.0281 | 5.58e-6 | 299 | 13 | 45 | 29 | 0.0913 |  | false |
| ale 100 | old core default | Success |  | 0.0281 |  | 1660 | 10 | 330 | 143 | 0.857 |  | false |
| ale_darcy 25 | KenCarp4 | Success | 6.06e-5 | 7.85e-7 | 0.000368 | 73 | 0 | 49 | 2 | 0.0382 |  | true |
| ale_darcy 25 | Rodas5P | Success | 6.09e-5 | 7.85e-7 | 0.0135 | 100 | 3 | 100 | 0 | 0.0379 |  | true |
| ale_darcy 25 | Rodas5P forced restart | Success | 6.09e-5 | 7.85e-7 | 0.0135 | 100 | 3 | 100 | 0 | 0.0379 |  | true |
| ale_darcy 25 | FBDF | Success | 8.43e-5 | 7.85e-7 | 0.0217 | 165 | 4 | 25 | 0 | 0.0207 |  | true |
| ale_darcy 25 | old core default | Success |  | 7.85e-7 |  | 75 | 2 | 49 | 2 | 0.0214 |  | false |
| ale_darcy 100 | KenCarp4 | Success | 2.05e-5 | 1.3e-6 | 0.00281 | 503 | 0 | 79 | 62 | 0.561 |  | true |
| ale_darcy 100 | Rodas5P | Success | 9.37e-5 | 1.3e-6 | 0.0048 | 102 | 2 | 102 | 0 | 0.156 |  | true |
| ale_darcy 100 | Rodas5P forced restart | Success | 9.37e-5 | 1.3e-6 | 0.0048 | 102 | 2 | 102 | 0 | 0.156 |  | true |
| ale_darcy 100 | FBDF | Success | 0.000177 | 1.3e-6 | 0.00441 | 184 | 6 | 26 | 0 | 0.0751 |  | true |
| ale_darcy 100 | old core default | Success |  | 1.3e-6 |  | 80 | 2 | 51 | 2 | 0.088 |  | false |
| V-CPL-02 (58 events) | KenCarp4 | Success | 0.00221 | 0.00028 | 0.4 | 2473 | 43 | 648 | 581 | 2.69 | 0 | true |
| V-CPL-02 (58 events) | Rodas5P | Success | 0.000599 | 0.00028 | 0.244 | 1904 | 741 | 1904 | 0 | 2.87 | 0 | true |
| V-CPL-02 (58 events) | Rodas5P forced restart | Success | 0.00132 | 0.00028 | 0.256 | 2049 | 2 | 2049 | 0 | 2.73 | 0 | true |
| V-CPL-02 (0 events) | FBDF | DtLessThanMin |  | 0.00028 | 0.322 | 30 | 5 | 2 | 9 | 0.024 |  | false |
| V-CPL-02 | old core default | Success |  | 0.00028 |  | 1037 | 0 | 497 | 78 | 1.48 |  | false |
