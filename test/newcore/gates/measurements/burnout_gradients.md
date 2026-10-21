# Gate item 9: smooth finite-difference gradients through burnout

Produced by `test/newcore/gates/burnout_gradients.jl` at commit 5f08134. Julia 1.13.0, 1 thread; load average [20.4, 18.8, 17.5] at the start and [26.3, 23.9, 20.3] at the end. V-BO-01 (temperature form), the finite-volume core through the case's solve at the package defaults, the old core with the case's keywords. Each parameter perturbed by ε ∈ {0, ±1, ±2, ±4, ±8}·1e-3; D(h) = (Q(h) − Q(−h))/(2h) at h ∈ {1, 2, 4, 8}·1e-3. Spread: (max − min)/|mean| of the four differences (smooth below 0.1, moderate below 1); σ: the RMS residual of a quadratic in ε over the nine runs; h10 = 10σ/|slope|. Item 9 passes when the interpolated peak, its time, the interpolated burnout time and the invocation time are smooth on the finite-volume core and the invocation count is the same in every run. Misfits: mean squares against the Gpyro V2 MLR [(g m⁻² s⁻¹)²] and back temperature [K²], reported only.

## PMMA absorption coefficient, 60 PMMA cells

**finite-volume core.** Every run Success; invocations [1], safety-net fallbacks 0 of 9; 1127 steps at ε = 0 (1104 to 1160).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3909 | -2.48 | -2.49 | -2.48 | -2.49 | 0.000761 | smooth | 3.44e-6 | 1.39e-5 |
| t_peak_interp | 334.846 | 33.6 | 33.2 | 33.8 | 34.1 | 0.0243 | smooth | 0.00105 | 0.00031 |
| t_bo_interp | 362.976 | 18.6 | 18.4 | 19.1 | 19.2 | 0.0409 | smooth | 0.001 | 0.000525 |
| invocation time | 357.915 | 2.98 | 2.57 | 2.97 | 2.96 | 0.143 | moderate | 0.000554 | 0.00188 |
| misfit_mlr | 19.4302 | -24.8 | -24.5 | -25.7 | -29.0 | 0.174 | moderate | 0.00865 | 0.00307 |
| misfit_Tback | 7263.28 | 1400.0 | 1380.0 | 1400.0 | 1400.0 | 0.0132 | smooth | 0.0257 | 0.000184 |

**old core.** Every run Success; parks [60]; 853 steps at ε = 0 (829 to 880).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3922 | -2.49 | -2.48 | -2.49 | -2.48 | 0.00391 | smooth | 7.58e-6 | 3.05e-5 |
| t_peak_interp | 334.828 | 31.7 | 35.3 | 33.7 | 34.1 | 0.106 | moderate | 0.00286 | 0.000841 |
| t_bo_interp | 363.353 | 18.1 | 22.1 | 19.8 | 19.6 | 0.199 | moderate | 0.00387 | 0.00196 |
| first park | 169.788 | 867.0 | -127.0 | 17.3 | -43.8 | 5.57 | erratic | 0.812 | 0.319 |
| misfit_mlr | 18.7587 | -102.0 | -78.0 | -71.9 | -73.1 | 0.372 | moderate | 0.0283 | 0.00385 |
| misfit_Tback | 7254.77 | 1260.0 | 1460.0 | 1370.0 | 1370.0 | 0.144 | moderate | 0.159 | 0.00116 |

Item 9 on this configuration: **fail**.

## incident flux, 60 PMMA cells

**finite-volume core.** Every run Success; invocations [1], safety-net fallbacks 0 of 9; 1127 steps at ε = 0 (1107 to 1154).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3909 | 42.6 | 42.6 | 42.6 | 42.6 | 7.77e-5 | smooth | 3.81e-6 | 8.94e-7 |
| t_peak_interp | 334.846 | -399.0 | -398.0 | -398.0 | -398.0 | 0.00283 | smooth | 0.000961 | 2.41e-5 |
| t_bo_interp | 362.976 | -457.0 | -457.0 | -456.0 | -457.0 | 0.000951 | smooth | 0.00129 | 2.82e-5 |
| invocation time | 357.915 | -469.0 | -469.0 | -469.0 | -469.0 | 0.000968 | smooth | 0.000571 | 1.22e-5 |
| misfit_mlr | 19.4302 | 902.0 | 956.0 | 948.0 | 927.0 | 0.0588 | smooth | 0.0974 | 0.00104 |
| misfit_Tback | 7263.28 | -16600.0 | -16600.0 | -16600.0 | -16600.0 | 0.00101 | smooth | 0.0183 | 1.1e-5 |

**old core.** Every run Success; parks [60]; 853 steps at ε = 0 (819 to 888).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3922 | 42.6 | 42.6 | 42.6 | 42.6 | 0.000171 | smooth | 8.99e-6 | 2.11e-6 |
| t_peak_interp | 334.828 | -399.0 | -402.0 | -398.0 | -398.0 | 0.00927 | smooth | 0.00414 | 0.000104 |
| t_bo_interp | 363.353 | -456.0 | -459.0 | -456.0 | -455.0 | 0.0089 | smooth | 0.00507 | 0.000111 |
| first park | 169.788 | -1100.0 | -655.0 | -127.0 | -157.0 | 1.91 | erratic | 0.814 | 0.0438 |
| misfit_mlr | 18.7587 | 949.0 | 873.0 | 862.0 | 882.0 | 0.0977 | smooth | 0.0956 | 0.00109 |
| misfit_Tback | 7254.77 | -16600.0 | -16800.0 | -16700.0 | -16700.0 | 0.0115 | smooth | 0.153 | 9.19e-5 |

Item 9 on this configuration: **pass**.

## pre-exponential factor, 60 PMMA cells

**finite-volume core.** Every run Success; invocations [1], safety-net fallbacks 0 of 9; 1127 steps at ε = 0 (1092 to 1127).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3909 | 1.79 | 1.79 | 1.79 | 1.79 | 0.00105 | smooth | 2.77e-6 | 1.55e-5 |
| t_peak_interp | 334.846 | -24.2 | -24.1 | -24.3 | -24.4 | 0.0139 | smooth | 0.000584 | 0.000239 |
| t_bo_interp | 362.976 | -26.6 | -27.0 | -27.4 | -27.4 | 0.03 | smooth | 0.00114 | 0.000416 |
| invocation time | 357.915 | -28.7 | -28.5 | -28.7 | -28.6 | 0.00614 | smooth | 0.000247 | 8.64e-5 |
| misfit_mlr | 19.4302 | 32.1 | 32.1 | 33.1 | 42.2 | 0.29 | moderate | 0.0226 | 0.00565 |
| misfit_Tback | 7263.28 | -651.0 | -642.0 | -651.0 | -647.0 | 0.0138 | smooth | 0.0111 | 0.000172 |

**old core.** Every run Success; parks [60]; 853 steps at ε = 0 (830 to 873).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.3922 | 1.79 | 1.79 | 1.79 | 1.79 | 0.0011 | smooth | 2.86e-6 | 1.6e-5 |
| t_peak_interp | 334.828 | -24.0 | -24.1 | -24.8 | -25.1 | 0.0446 | smooth | 0.00151 | 0.000605 |
| t_bo_interp | 363.353 | -27.1 | -27.0 | -27.0 | -28.0 | 0.0352 | smooth | 0.00224 | 0.000809 |
| first park | 169.788 | -250.0 | 516.0 | -241.0 | -53.0 | 109.0 | erratic | 0.812 | 0.127 |
| misfit_mlr | 18.7587 | 104.0 | 85.5 | 101.0 | 96.5 | 0.195 | moderate | 0.0248 | 0.00256 |
| misfit_Tback | 7254.77 | -644.0 | -668.0 | -660.0 | -670.0 | 0.0403 | smooth | 0.0405 | 0.000607 |

Item 9 on this configuration: **pass**.

## PMMA absorption coefficient, 120 PMMA cells

**finite-volume core.** Every run Success; invocations [1], safety-net fallbacks 0 of 9; 2118 steps at ε = 0 (2085 to 2128).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.4691 | -2.54 | -2.53 | -2.53 | -2.53 | 0.000612 | smooth | 3.32e-6 | 1.31e-5 |
| t_peak_interp | 335.625 | 34.8 | 34.6 | 34.3 | 34.1 | 0.0208 | smooth | 0.000899 | 0.000263 |
| t_bo_interp | 363.008 | 20.6 | 20.2 | 19.9 | 19.9 | 0.0371 | smooth | 0.000953 | 0.000478 |
| invocation time | 357.711 | 3.31 | 3.08 | 2.92 | 2.96 | 0.128 | moderate | 0.00037 | 0.00125 |
| misfit_mlr | 19.4212 | -38.8 | -38.8 | -38.3 | -38.6 | 0.0146 | smooth | 0.00167 | 0.000432 |
| misfit_Tback | 7351.37 | 1380.0 | 1370.0 | 1360.0 | 1360.0 | 0.0153 | smooth | 0.0171 | 0.000126 |

**old core.** Every run Success; parks [120]; 1461 steps at ε = 0 (1431 to 1534).

| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |
|---|---|---|---|---|---|---|---|---|---|
| peak_interp | 32.4694 | -2.53 | -2.53 | -2.53 | -2.53 | 0.00138 | smooth | 7.54e-6 | 2.97e-5 |
| t_peak_interp | 335.621 | 35.0 | 35.3 | 34.9 | 33.9 | 0.041 | smooth | 0.00261 | 0.000763 |
| t_bo_interp | 363.227 | -29.4 | 21.4 | 20.9 | 20.4 | 6.09 | erratic | 0.0374 | 0.0188 |
| first park | 167.588 | -69.3 | -41.4 | -175.0 | -99.7 | 1.39 | erratic | 0.305 | 0.0275 |
| misfit_mlr | 19.0708 | 113.0 | -67.1 | -64.2 | -56.5 | 9.65 | erratic | 0.131 | 0.0231 |
| misfit_Tback | 7345.1 | 2680.0 | 1370.0 | 1350.0 | 1330.0 | 0.801 | moderate | 0.973 | 0.0072 |

Item 9 on this configuration: **fail**.

