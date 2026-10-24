# Leftover 3: the finite-slab gate

Produced by `test/newcore/gates/finite_slab.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [20.7, 20.3, 19.7] at the start and [22.8, 21.2, 20.2] at the end. 320 cells; warm medians of 3 solves; L2 and L∞ the temperature error against the exact series [K], the largest over the output times. Spatial floor (the default at reltol 1.0e-11, abstol 1.0e-13): V-HC-03 9.27497e-5 K, V-HC-02c 0.000593244 K. `passes`: L2 no larger than the old core's default at the same tolerance or than the floor plus 1 %, whichever is larger.

| case | reltol | run | retcode | steps | rejected | RHS | Jacobians | warm [s] | L2 max | L∞ max | passes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V-HC-03 | 1.0e-6 | old KenCarp4 (default) | Success | 52 | 0 | 325 | 5 | 0.0119 | 0.00589141 | 0.0100849 | true |
| V-HC-03 | 1.0e-6 | old KenCarp4 raw | Success | 57 | 0 | 355 | 5 | 0.0124 | 0.00490972 | 0.0084185 | true |
| V-HC-03 | 1.0e-6 | new KenCarp4 raw | Success | 57 | 0 | 355 | 5 | 0.0152 | 0.00490972 | 0.0084185 | true |
| V-HC-03 | 1.0e-6 | new Rodas5P (default) | Success | 54 | 0 | 486 | 54 | 0.0236 | 9.26857e-5 | 0.000122921 | true |
| V-HC-03 | 1.0e-6 | new Rodas5P forced restart | Success | 54 | 0 | 486 | 54 | 0.0238 | 9.26857e-5 | 0.000122921 | true |
| V-HC-03 | 1.0e-6 | new FBDF | Success | 148 | 9 | 219 | 1 | 0.0169 | 0.000106597 | 0.000181666 | true |
| V-HC-03 | 1.0e-8 | old KenCarp4 (default) | Success | 92 | 1 | 1499 | 7 | 0.0682 | 0.000107262 | 0.000193635 | true |
| V-HC-03 | 1.0e-8 | old KenCarp4 raw | Success | 110 | 1 | 1805 | 6 | 0.0787 | 9.11372e-5 | 0.000121087 | true |
| V-HC-03 | 1.0e-8 | new KenCarp4 raw | Success | 110 | 1 | 1805 | 6 | 0.0864 | 9.11372e-5 | 0.000121087 | true |
| V-HC-03 | 1.0e-8 | new Rodas5P (default) | Success | 97 | 0 | 873 | 97 | 0.0736 | 9.27401e-5 | 0.000122897 | true |
| V-HC-03 | 1.0e-8 | new Rodas5P forced restart | Success | 97 | 0 | 873 | 97 | 0.0737 | 9.27401e-5 | 0.000122897 | true |
| V-HC-03 | 1.0e-8 | new FBDF | Success | 301 | 7 | 393 | 1 | 0.0428 | 9.27449e-5 | 0.000123049 | true |
| V-HC-02c | 1.0e-6 | old KenCarp4 (default) | Success | 75 | 1 | 900 | 5 | 0.0408 | 0.000220394 | 0.000339036 | true |
| V-HC-02c | 1.0e-6 | old KenCarp4 raw | Success | 70 | 0 | 431 | 5 | 0.023 | 0.0593386 | 0.0794116 | false |
| V-HC-02c | 1.0e-6 | new KenCarp4 raw | Success | 70 | 0 | 431 | 5 | 0.0285 | 0.0593375 | 0.0794102 | false |
| V-HC-02c | 1.0e-6 | new Rodas5P (default) | Success | 61 | 0 | 549 | 61 | 0.0475 | 0.000593298 | 0.00112661 | true |
| V-HC-02c | 1.0e-6 | new Rodas5P forced restart | Success | 61 | 0 | 549 | 61 | 0.0475 | 0.000593298 | 0.00112661 | true |
| V-HC-02c | 1.0e-6 | new FBDF | Success | 169 | 3 | 243 | 1 | 0.0281 | 0.000671528 | 0.00134139 | false |
| V-HC-02c | 1.0e-8 | old KenCarp4 (default) | Success | 156 | 0 | 3423 | 12 | 0.139 | 0.000592804 | 0.00112516 | true |
| V-HC-02c | 1.0e-8 | old KenCarp4 raw | Success | 127 | 0 | 2884 | 16 | 0.126 | 0.000589916 | 0.00111724 | true |
| V-HC-02c | 1.0e-8 | new KenCarp4 raw | Success | 127 | 0 | 2884 | 16 | 0.133 | 0.000589916 | 0.00111724 | true |
| V-HC-02c | 1.0e-8 | new Rodas5P (default) | Success | 112 | 0 | 1008 | 112 | 0.0832 | 0.000593245 | 0.00112637 | true |
| V-HC-02c | 1.0e-8 | new Rodas5P forced restart | Success | 112 | 0 | 1008 | 112 | 0.0873 | 0.000593245 | 0.00112637 | true |
| V-HC-02c | 1.0e-8 | new FBDF | Success | 356 | 9 | 510 | 1 | 0.0517 | 0.000593121 | 0.00112509 | true |
