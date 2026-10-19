# Leftover 3: the finite-slab gate

Produced by `test/newcore/gates/finite_slab.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [17.7, 18.2, 18.3] at the start and [19.6, 18.7, 18.5] at the end. 320 cells; warm medians of 3 solves; L2 and L∞ the temperature error against the exact series [K], the largest over the output times. `passes`: L2 no larger than the old core's default at the same tolerance.

| case | reltol | run | retcode | steps | rejected | RHS | Jacobians | warm [s] | L2 max | L∞ max | passes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V-HC-03 | 1.0e-6 | old KenCarp4 (default) | Success | 52 | 0 | 325 | 5 | 0.0133 | 0.00589141 | 0.0100849 | true |
| V-HC-03 | 1.0e-6 | old KenCarp4 raw | Success | 57 | 0 | 355 | 5 | 0.0121 | 0.00490972 | 0.0084185 | true |
| V-HC-03 | 1.0e-6 | new KenCarp4 raw | Success | 57 | 0 | 355 | 5 | 0.0154 | 0.00490972 | 0.0084185 | true |
| V-HC-03 | 1.0e-6 | new Rodas5P | Success | 54 | 0 | 486 | 54 | 0.0252 | 9.26857e-5 | 0.000122921 | true |
| V-HC-03 | 1.0e-6 | new Rodas5P forced restart | Success | 54 | 0 | 486 | 54 | 0.0254 | 9.26857e-5 | 0.000122921 | true |
| V-HC-03 | 1.0e-6 | new FBDF | Success | 148 | 9 | 219 | 1 | 0.0155 | 0.000106597 | 0.000181666 | true |
| V-HC-03 | 1.0e-8 | old KenCarp4 (default) | Success | 92 | 1 | 1499 | 7 | 0.0377 | 0.000107262 | 0.000193635 | true |
| V-HC-03 | 1.0e-8 | old KenCarp4 raw | Success | 110 | 1 | 1805 | 6 | 0.044 | 9.11372e-5 | 0.000121087 | true |
| V-HC-03 | 1.0e-8 | new KenCarp4 raw | Success | 110 | 1 | 1805 | 6 | 0.0514 | 9.11372e-5 | 0.000121087 | true |
| V-HC-03 | 1.0e-8 | new Rodas5P | Success | 97 | 0 | 873 | 97 | 0.043 | 9.27401e-5 | 0.000122897 | true |
| V-HC-03 | 1.0e-8 | new Rodas5P forced restart | Success | 97 | 0 | 873 | 97 | 0.0435 | 9.27401e-5 | 0.000122897 | true |
| V-HC-03 | 1.0e-8 | new FBDF | Success | 301 | 7 | 393 | 1 | 0.0262 | 9.27449e-5 | 0.000123049 | true |
| V-HC-02c | 1.0e-6 | old KenCarp4 (default) | Success | 75 | 1 | 900 | 5 | 0.0241 | 0.000220394 | 0.000339036 | true |
| V-HC-02c | 1.0e-6 | old KenCarp4 raw | Success | 70 | 0 | 431 | 5 | 0.0139 | 0.0593386 | 0.0794116 | false |
| V-HC-02c | 1.0e-6 | new KenCarp4 raw | Success | 70 | 0 | 431 | 5 | 0.0177 | 0.0593375 | 0.0794102 | false |
| V-HC-02c | 1.0e-6 | new Rodas5P | Success | 61 | 0 | 549 | 61 | 0.0297 | 0.000593298 | 0.00112661 | false |
| V-HC-02c | 1.0e-6 | new Rodas5P forced restart | Success | 61 | 0 | 549 | 61 | 0.0302 | 0.000593298 | 0.00112661 | false |
| V-HC-02c | 1.0e-6 | new FBDF | Success | 169 | 3 | 243 | 1 | 0.0184 | 0.000671528 | 0.00134139 | false |
| V-HC-02c | 1.0e-8 | old KenCarp4 (default) | Success | 156 | 0 | 3423 | 12 | 0.0813 | 0.000592804 | 0.00112516 | true |
| V-HC-02c | 1.0e-8 | old KenCarp4 raw | Success | 127 | 0 | 2884 | 16 | 0.0698 | 0.000589916 | 0.00111724 | true |
| V-HC-02c | 1.0e-8 | new KenCarp4 raw | Success | 127 | 0 | 2884 | 16 | 0.0809 | 0.000589916 | 0.00111724 | true |
| V-HC-02c | 1.0e-8 | new Rodas5P | Success | 112 | 0 | 1008 | 112 | 0.0516 | 0.000593245 | 0.00112637 | false |
| V-HC-02c | 1.0e-8 | new Rodas5P forced restart | Success | 112 | 0 | 1008 | 112 | 0.0508 | 0.000593245 | 0.00112637 | false |
| V-HC-02c | 1.0e-8 | new FBDF | Success | 356 | 9 | 510 | 1 | 0.0312 | 0.000593121 | 0.00112509 | false |
