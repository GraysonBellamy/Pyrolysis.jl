# Mesh phase 0 baselines: step-update flag × integrator

Scripts: `scripts/mesh_phase0/baselines/integrator_flag.jl` (shared helpers in `common.jl`); the one-callback parking measurement is reported in `one_callback_parking.md`. Data: `integrator_flag.csv` in this directory; log alongside. Run as `julia --project=test --startup-file=no --threads=1 <script>` from the repository root. Environment: Julia 1.13.0, OrdinaryDiffEq 7.8.1 (Core 4.18.0, SDIRK 2.9.6, BDF 2.4.11, Rosenbrock 2.7.4), DiffEqBase 7.21.2, SciMLBase 3.56.1, LinearSolve 5.18.1 (KLU everywhere).

## 1. Step-update callback variants × integrator × tolerance

Variants of the per-step callback (`create_step_update_callback`, replaced in-process): `current` = mirror + `derivative_discontinuity!(integrator, true)` (package); `flag_off` = mirror, no flag call; `observer` = mirror + `derivative_discontinuity!(integrator, false)` (`u_modified!` is a deprecated alias of the same function). Every integrator uses `autodiff = AutoFiniteDiff()`, `linsolve = KLUFactorization()`; `KenCarp4_raw` is `KenCarp4(smooth_est = false)`. Tolerance pairs (reltol, abstol): (1e-6, 1e-8) = `solve` default, (1e-8, 1e-10) = the cases' own `solve_kwargs`. Errors are against the exact series at every output time after t = 0: L2 = RMS over cells, Linf = max over cells; "max" = maximum over output times, "final" = last output time (Fo = 1.5 for V-HC-03, Fo = 0.5 for V-HC-02c). Warm time = median of 3 solves after one warm-up; `maxiters = 20 000`; wall-clock guard 300 s per solve (never triggered). Cells in the tables read `reltol 1e-6 / reltol 1e-8`.

### V-HC-03 (finite slab, constant flux, adiabatic back), 320 cells

| integrator | variant | retcode (t reached) | steps | rej | nf | njacs | nlfail | warm s | L2 max (K) | Linf max (K) | L2 final (K) | QoI fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KenCarp4 | current | ok / ok | 52 / 92 | 0 / 1 | 325 / 1499 | 5 / 7 | 0 / 0 | 0.0175 / 0.06 | 0.00589 / 1.07e-04 | 0.0101 / 1.94e-04 | 8.14e-05 / 8.14e-05 | 1 / 0 |
| KenCarp4 | flag_off | ok / ok | 52 / 92 | 0 / 1 | 325 / 1499 | 5 / 7 | 0 / 0 | 0.0173 / 0.0602 | 0.00589 / 1.07e-04 | 0.0101 / 1.94e-04 | 8.14e-05 / 8.14e-05 | 1 / 0 |
| KenCarp4 | observer | ok / ok | 58 / 188 | 0 / 0 | 304 / 956 | 5 / 5 | 0 / 0 | 0.0173 / 0.0458 | 0.0106 / 0.00553 | 0.018 / 0.00938 | 8.14e-05 / 8.14e-05 | 1 / 1 |
| KenCarp4_raw | current | ok / ok | 57 / 110 | 0 / 1 | 355 / 1805 | 5 / 6 | 0 / 0 | 0.0184 / 0.0717 | 0.00491 / 9.11e-05 | 0.00842 / 1.21e-04 | 8.14e-05 / 8.14e-05 | 1 / 0 |
| KenCarp4_raw | flag_off | ok / ok | 57 / 110 | 0 / 1 | 355 / 1805 | 5 / 6 | 0 / 0 | 0.0185 / 0.0714 | 0.00491 / 9.11e-05 | 0.00842 / 1.21e-04 | 8.14e-05 / 8.14e-05 | 1 / 0 |
| KenCarp4_raw | observer | ok / ok | 56 / 285 | 0 / 0 | 294 / 1440 | 5 / 5 | 0 / 0 | 0.0179 / 0.0642 | 0.00515 / 0.00401 | 0.00876 / 0.00681 | 8.14e-05 / 8.14e-05 | 1 / 1 |
| FBDF | current | MaxIters (t=3.1876831450039638) / MaxIters (t=0.007679514116566112) | 19999 / 19997 | 1 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | flag_off | MaxIters (t=3.1876831450039638) / MaxIters (t=0.007679514116566112) | 19999 / 19997 | 1 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | observer | ok / ok | 148 / 301 | 9 / 7 | 219 / 393 | 1 / 1 | 0 / 0 | 0.0172 / 0.0284 | 1.07e-04 / 9.27e-05 | 1.82e-04 / 1.23e-04 | 9.13e-05 / 8.14e-05 | 1 / 0 |
| Rodas5P | current | ok / ok | 54 / 97 | 0 / 0 | 486 / 873 | 54 / 97 | 0 / 0 | 0.0309 / 0.0538 | 9.27e-05 / 9.27e-05 | 1.23e-04 / 1.23e-04 | 8.14e-05 / 8.14e-05 | 0 / 0 |
| Rodas5P | flag_off | ok / ok | 54 / 97 | 0 / 0 | 486 / 873 | 54 / 97 | 0 / 0 | 0.0307 / 0.0536 | 9.27e-05 / 9.27e-05 | 1.23e-04 / 1.23e-04 | 8.14e-05 / 8.14e-05 | 0 / 0 |
| Rodas5P | observer | ok / ok | 54 / 97 | 0 / 0 | 486 / 873 | 54 / 97 | 0 / 0 | 0.0307 / 0.0542 | 9.27e-05 / 9.27e-05 | 1.23e-04 / 1.23e-04 | 8.14e-05 / 8.14e-05 | 0 / 0 |

Ratios against `current` for the same integrator and tolerance (V-HC-03, 320 cells):

| integrator | variant | reltol 1e-6: L2max× / Linfmax× / steps× | reltol 1e-8: L2max× / Linfmax× / steps× |
|---|---|---|---|
| KenCarp4 | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4 | observer | 1.80 / 1.79 / 1.12 | 51.60 / 48.46 / 2.04 |
| KenCarp4_raw | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4_raw | observer | 1.05 / 1.04 / 0.98 | 44.00 / 56.24 / 2.59 |
| FBDF | flag_off | n/a / n/a / 1.00 | n/a / n/a / 1.00 |
| FBDF | observer | n/a / n/a / 0.01 | n/a / n/a / 0.02 |
| Rodas5P | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| Rodas5P | observer | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |

### V-HC-03 (finite slab, constant flux, adiabatic back), 80 cells

| integrator | variant | retcode (t reached) | steps | rej | nf | njacs | nlfail | warm s | L2 max (K) | Linf max (K) | L2 final (K) | QoI fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KenCarp4 | current | ok / ok | 54 / 86 | 0 / 0 | 334 / 1817 | 5 / 7 | 0 / 0 | 0.00473 / 0.0181 | 0.0126 / 0.00148 | 0.0226 / 0.00197 | 0.0013 / 0.0013 | 1 / 0 |
| KenCarp4 | flag_off | ok / ok | 54 / 86 | 0 / 0 | 334 / 1817 | 5 / 7 | 0 / 0 | 0.00476 / 0.0182 | 0.0126 / 0.00148 | 0.0226 / 0.00197 | 0.0013 / 0.0013 | 1 / 0 |
| KenCarp4 | observer | ok / ok | 56 / 86 | 0 / 0 | 291 / 1733 | 5 / 7 | 0 / 0 | 0.00458 / 0.0178 | 0.00903 / 0.00148 | 0.0163 / 0.00197 | 0.0013 / 0.0013 | 2 / 0 |
| KenCarp4_raw | current | ok / ok | 49 / 94 | 0 / 0 | 305 / 1960 | 5 / 8 | 0 / 0 | 0.00432 / 0.0197 | 0.00773 / 0.00148 | 0.0143 / 0.00197 | 0.0013 / 0.0013 | 3 / 0 |
| KenCarp4_raw | flag_off | ok / ok | 49 / 94 | 0 / 0 | 305 / 1960 | 5 / 8 | 0 / 0 | 0.0044 / 0.0195 | 0.00773 / 0.00148 | 0.0143 / 0.00197 | 0.0013 / 0.0013 | 3 / 0 |
| KenCarp4_raw | observer | ok / ok | 51 / 94 | 0 / 0 | 266 / 1868 | 5 / 8 | 0 / 0 | 0.00423 / 0.0194 | 0.00474 / 0.00148 | 0.00886 / 0.00197 | 0.0013 / 0.0013 | 2 / 0 |
| FBDF | current | MaxIters (t=3.2420762047817844) / MaxIters (t=0.013347140571771929) | 20000 / 19997 | 0 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | flag_off | MaxIters (t=3.2420762047817844) / MaxIters (t=0.013347140571771929) | 20000 / 19997 | 0 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | observer | ok / ok | 131 / 276 | 4 / 7 | 194 / 392 | 1 / 1 | 0 / 0 | 0.00418 / 0.00781 | 0.00151 / 0.00148 | 0.00195 / 0.00197 | 0.0013 / 0.0013 | 1 / 0 |
| Rodas5P | current | ok / ok | 48 / 86 | 0 / 0 | 432 / 774 | 48 / 86 | 0 / 0 | 0.00709 / 0.0121 | 0.00148 / 0.00148 | 0.00197 / 0.00197 | 0.0013 / 0.0013 | 0 / 0 |
| Rodas5P | flag_off | ok / ok | 48 / 86 | 0 / 0 | 432 / 774 | 48 / 86 | 0 / 0 | 0.00705 / 0.0121 | 0.00148 / 0.00148 | 0.00197 / 0.00197 | 0.0013 / 0.0013 | 0 / 0 |
| Rodas5P | observer | ok / ok | 48 / 86 | 0 / 0 | 432 / 774 | 48 / 86 | 0 / 0 | 0.00724 / 0.0122 | 0.00148 / 0.00148 | 0.00197 / 0.00197 | 0.0013 / 0.0013 | 0 / 0 |

Ratios against `current` for the same integrator and tolerance (V-HC-03, 80 cells):

| integrator | variant | reltol 1e-6: L2max× / Linfmax× / steps× | reltol 1e-8: L2max× / Linfmax× / steps× |
|---|---|---|---|
| KenCarp4 | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4 | observer | 0.71 / 0.72 / 1.04 | 1.00 / 1.00 / 1.00 |
| KenCarp4_raw | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4_raw | observer | 0.61 / 0.62 / 1.04 | 1.00 / 1.00 / 1.00 |
| FBDF | flag_off | n/a / n/a / 1.00 | n/a / n/a / 1.00 |
| FBDF | observer | n/a / n/a / 0.01 | n/a / n/a / 0.01 |
| Rodas5P | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| Rodas5P | observer | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |

### V-HC-02c (finite slab, convective BC, Bi = 10), 320 cells

| integrator | variant | retcode (t reached) | steps | rej | nf | njacs | nlfail | warm s | L2 max (K) | Linf max (K) | L2 final (K) | QoI fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KenCarp4 | current | ok / ok | 75 / 156 | 1 / 0 | 900 / 3423 | 5 / 12 | 0 / 0 | 0.0386 / 0.133 | 2.20e-04 / 5.93e-04 | 3.39e-04 / 0.00113 | 1.66e-04 / 1.74e-04 | 0 / 0 |
| KenCarp4 | flag_off | ok / ok | 75 / 156 | 1 / 0 | 900 / 3423 | 5 / 12 | 0 / 0 | 0.0387 / 0.133 | 2.20e-04 / 5.93e-04 | 3.39e-04 / 0.00113 | 1.66e-04 / 1.74e-04 | 0 / 0 |
| KenCarp4 | observer | ok / ok | 151 / 156 | 0 / 0 | 766 / 3284 | 5 / 11 | 0 / 0 | 0.0374 / 0.13 | 0.0194 / 5.93e-04 | 0.0369 / 0.00113 | 0.019 / 1.74e-04 | 1 / 0 |
| KenCarp4_raw | current | ok / ok | 70 / 127 | 0 / 0 | 431 / 2884 | 5 / 16 | 0 / 0 | 0.0219 / 0.117 | 0.0593 / 5.90e-04 | 0.0794 / 0.00112 | 0.0593 / 1.74e-04 | 2 / 0 |
| KenCarp4_raw | flag_off | ok / ok | 70 / 127 | 0 / 0 | 431 / 2884 | 5 / 16 | 0 / 0 | 0.0216 / 0.117 | 0.0593 / 5.90e-04 | 0.0794 / 0.00112 | 0.0593 / 1.74e-04 | 2 / 0 |
| KenCarp4_raw | observer | ok / ok | 74 / 127 | 0 / 0 | 382 / 2756 | 5 / 16 | 0 / 0 | 0.0216 / 0.114 | 0.0607 / 5.90e-04 | 0.0812 / 0.00112 | 0.0607 / 1.74e-04 | 3 / 0 |
| FBDF | current | MaxIters (t=4.344456554001141) / MaxIters (t=0.01690266246315271) | 20000 / 19997 | 0 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | flag_off | MaxIters (t=4.344456554001141) / MaxIters (t=0.01690266246315271) | 20000 / 19997 | 0 / 3 | 40005 / 40004 | 1 / 1 | 0 / 0 | NaN / NaN | NaN / NaN | NaN / NaN | NaN / NaN | -1 / -1 |
| FBDF | observer | ok / ok | 169 / 356 | 3 / 9 | 243 / 510 | 1 / 1 | 0 / 0 | 0.0181 / 0.0347 | 6.72e-04 / 5.93e-04 | 0.00134 / 0.00113 | 1.81e-04 / 1.74e-04 | 0 / 0 |
| Rodas5P | current | ok / ok | 61 / 112 | 0 / 0 | 549 / 1008 | 61 / 112 | 0 / 0 | 0.0348 / 0.0625 | 5.93e-04 / 5.93e-04 | 0.00113 / 0.00113 | 1.74e-04 / 1.74e-04 | 0 / 0 |
| Rodas5P | flag_off | ok / ok | 61 / 112 | 0 / 0 | 549 / 1008 | 61 / 112 | 0 / 0 | 0.035 / 0.0624 | 5.93e-04 / 5.93e-04 | 0.00113 / 0.00113 | 1.74e-04 / 1.74e-04 | 0 / 0 |
| Rodas5P | observer | ok / ok | 61 / 112 | 0 / 0 | 549 / 1008 | 61 / 112 | 0 / 0 | 0.0349 / 0.0627 | 5.93e-04 / 5.93e-04 | 0.00113 / 0.00113 | 1.74e-04 / 1.74e-04 | 0 / 0 |

Ratios against `current` for the same integrator and tolerance (V-HC-02c, 320 cells):

| integrator | variant | reltol 1e-6: L2max× / Linfmax× / steps× | reltol 1e-8: L2max× / Linfmax× / steps× |
|---|---|---|---|
| KenCarp4 | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4 | observer | 88.00 / 108.92 / 2.01 | 1.00 / 1.00 / 1.00 |
| KenCarp4_raw | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4_raw | observer | 1.02 / 1.02 / 1.06 | 1.00 / 1.00 / 1.00 |
| FBDF | flag_off | n/a / n/a / 1.00 | n/a / n/a / 1.00 |
| FBDF | observer | n/a / n/a / 0.01 | n/a / n/a / 0.02 |
| Rodas5P | flag_off | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| Rodas5P | observer | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |

### Accuracy and cost relative to the package default (KenCarp4, current, same tolerance), V-HC-03 at 320 cells

| integrator | variant | reltol 1e-6: L2max× / steps× / time× | reltol 1e-8: L2max× / steps× / time× |
|---|---|---|---|
| KenCarp4 | current | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| KenCarp4 | flag_off | 1.00 / 1.00 / 0.99 | 1.00 / 1.00 / 1.00 |
| KenCarp4 | observer | 1.80 / 1.12 / 0.99 | 51.60 / 2.04 / 0.76 |
| KenCarp4_raw | current | 0.83 / 1.10 / 1.06 | 0.85 / 1.20 / 1.19 |
| KenCarp4_raw | flag_off | 0.83 / 1.10 / 1.06 | 0.85 / 1.20 / 1.19 |
| KenCarp4_raw | observer | 0.87 / 1.08 / 1.03 | 37.38 / 3.10 / 1.07 |
| FBDF | current | n/a / 384.60 / n/a | n/a / 217.36 / n/a |
| FBDF | flag_off | n/a / 384.60 / n/a | n/a / 217.36 / n/a |
| FBDF | observer | 0.02 / 2.85 / 0.99 | 0.86 / 3.27 / 0.47 |
| Rodas5P | current | 0.02 / 1.04 / 1.77 | 0.86 / 1.05 / 0.90 |
| Rodas5P | flag_off | 0.02 / 1.04 / 1.76 | 0.86 / 1.05 / 0.89 |
| Rodas5P | observer | 0.02 / 1.04 / 1.76 | 0.86 / 1.05 / 0.90 |

### Does the docstring's 15× / 2.6× effect reproduce?

- `flag_off` is identical to `current` in every counter and every error on all 24 (case, cells, integrator, tolerance) pairs: DiffEqBase 7.21.2 sets `derivative_discontinuity = true` before every `affect!`, so not calling the function leaves the flag on. Only the explicit `false` (`observer`) changes anything; the docstring's "removing the flag" corresponds to `observer`.
- V-HC-03, 320 cells, the case's tolerance (reltol 1e-8): KenCarp4 current → observer: L2 max 1.07e-04 → 0.00553 K (51.6×), steps 92 → 188 (2.04×); KenCarp4 smooth_est=false: 9.11e-05 → 0.00401 K (44.0×), steps 110 → 285 (2.59×). The direction and magnitude of the docstring's claim reproduce (error ratio larger than 15×, step ratio 2.0–2.6×).
- At the `solve` default tolerance (reltol 1e-6) the same comparison gives 1.80× / 1.12× on V-HC-03 (both variants sit at L2 max ≈ 0.00589–0.0106 K, failing the case's 0.006 K L∞ tolerance), while on V-HC-02c the effect appears at reltol 1e-6 (88× L2 max, 2.01× steps) and vanishes at reltol 1e-8 (1.00×). At 80 cells the observer variant is as accurate or more accurate (0.6–1.0×). The benefit of the flag for KenCarp4 is therefore real but depends on tolerance and mesh; the two finite-slab cases cannot be told apart as the docstring's source.
- The final-time error (L2 final 8.14e-5 K at 320 cells, 1.30e-3 K at 80 cells on V-HC-03) is the same for every finished configuration: it is the spatial discretisation error of the steady-shape profile; the variants differ only in the transient outputs.

### FBDF

- With `current` or `flag_off`, FBDF never finishes: every accepted step is followed by a flagged discontinuity, FBDF's `perform_step!` resets to order 1 and clears its difference history (`iters_from_event` is 1 at every step), and its step controller never grows the step. `dt` stays at 8.2e-5 s (the initial 1e-4 s after one early rejection) at reltol 1e-6 and at 3.8e-7 s at reltol 1e-8, so 20 000 accepted steps reach t = 3.19 s and 0.0077 s of 1200 s (`MaxIters`; one Jacobian, one Newton iteration per step). A probe with `maxiters = 100 000` stalled the same way. The same happens at 80 cells and on V-HC-02c.
- With `observer`, FBDF finishes with one Jacobian: 148 steps / 219 f-evals and L2 max 1.07e-04 K at reltol 1e-6, which is 55× below KenCarp4 current at the same tolerance (0.00589 K) and equal to what KenCarp4 current reaches at reltol 1e-8 (1.07e-04 K); at reltol 1e-8, 301 steps and L2 max 9.27e-05 K (the 320-cell spatial floor, equal to Rodas5P). Warm time 0.0172 s / 0.0284 s against 0.0175 s / 0.06 s for KenCarp4 current. FBDF therefore requires the observer variant.

### Rodas5P

- Bit-identical under all three variants (a Rosenbrock method rebuilds W every step: njacs = steps). L2 max 9.27e-05 K at both tolerances (the spatial floor), 54 / 97 steps, warm 0.0309 / 0.0538 s against 0.0175 / 0.06 s for KenCarp4 current; on V-HC-02c Rodas5P, FBDF-observer and KenCarp4 at reltol 1e-8 agree at L2 max 5.9e-4 K (the KenCarp4-current value of 2.2e-4 K at reltol 1e-6 lies below that floor through cancellation between temporal and spatial error, not through a more accurate time integration).

## 2. One-callback parking

Moved to `one_callback_parking.md` (generated by `one_callback_parking.jl`; data in `one_callback_parking.csv`, `one_callback_parking_qois.csv`, `one_callback_parking_ledger.csv`). That file holds the corrected loop (it re-tests the callback's own event predicate before every further `deplete!`, so a cascade that leaves exactly `min_cells` live cells no longer terminates the solve), the `min_cells = 20` V-BO-04 termination case, the QoI tables, and the output-spacing check of the sampled mass ledger.

## Caveats

- Warm times are medians of 3 solves on a shared 22-core machine under a load average of about 17 from other processes; repeated runs of the conduction matrix varied by up to ~40% in absolute warm time (0.013–0.018 s for KenCarp4 current at 320 cells, reltol 1e-6) while every ratio between configurations within a run stayed within ±0.02. Counters and errors are identical between runs. Driver wall time: integrator_flag 146–170 s.
- The 12 FBDF `current`/`flag_off` rows stopped at `maxiters = 20 000` (reduced from 100 000 after a probe stalled identically); they have no warm timing and no errors (`repeats_consistent = false` only there). No configuration timed out or threw.
- The step-update override adds a `time()` deadline check per accepted step to every variant, including `current`; the counters and errors of `current` match the package code run without it (V-HC-03 at 320/80 cells: 92/86 steps, identical QoIs).
- "QoI fails" in the conduction tables counts the cases' own QoIs at the given tolerance; at reltol 1e-6 several finished configurations exceed the cases' 0.006 K L∞ tolerance because the cases themselves run at reltol 1e-8.
