# Inverse Analysis Plan — White Pine Component Thermal Conductivities

## Goal

Recover the temperature-independent thermal conductivities of the virgin and char
solid components of white pine (and, if needed, the four pyrolysis intermediates)
by fitting CAPA back-surface temperature time-series. All other thermophysical
and kinetic properties are held fixed at their independently measured values;
the cone heat flux conditions used for fitting are **30 kW/m² and 50 kW/m²**.
The **70 kW/m²** record is held out for validation.

---

## Status (high-level)

| Area                                | Status     |
|-------------------------------------|------------|
| Project scaffolding & files                   | ✅ done     |
| White pine properties from JSON               | ✅ wired in |
| Forward model (Stages 1/2/3)                  | ✅ done     |
| BC scenario container                         | ✅ done     |
| Replicate-averaging data loader               | ✅ done     |
| BC time-series (Inner/Outer wall T)           | ✅ wired in |
| Top-gas analytic model                        | ✅ wired in |
| Loss + trace log (data/prior split)           | ✅ done     |
| HFM/LFA k(T) prior penalty                    | ✅ done     |
| Multi-start NelderMead via Fminbox            | ✅ done     |
| Sobol multi-start + grid-minima seeding       | ✅ done     |
| Per-loss parallel scenarios (`Threads.@spawn`) | ✅ done     |
| Hessian covariance (rederived) + bootstrap    | ✅ done     |
| Grid scan (Stage 1)                           | ✅ done     |
| End-to-end smoke test                         | ✅ pass     |
| **First fit run**                             | ⏳ pending  |
| **Validation against 70 kW/m²**               | ⏳ pending  |

Per-forward-solve wall time observed in trace: **~25 s** per scenario at
160 cells, 601 output samples (`KenCarp4` + Sparspak). With `julia -t auto`
the two fit scenarios run in parallel inside each loss call → effective
**~25 s/eval**. NelderMead at 2 params (Stage 1) typically converges in
30–60 evals → **12–25 min/start**. Multi-start (3 starts) Stage 1: roughly
**30–60 minutes total**. Stage 3 (4 params): 80–150 evals/start, expect
30–60 min/start.

---

## Architecture

```
examples/inverse_analysis/
├── PLAN.md                         (this document)
├── wp_properties.jl                ✅ all white pine constants from MaCFP JSON
├── wp_forward.jl                   ✅ build_and_solve(θ, scenario; stage)
├── wp_objective.jl                 ✅ replicate loader + loss + trace log
├── inverse_white_pine.jl           ✅ multi-start NelderMead + covariance + plots
├── data/
│   ├── macfp_white_pine_property_file.json   ✅ source of all WP properties
│   ├── UMD_Wood_CAPA_N2_30kW_R{1,2,3}.csv     ✅ 3 replicates (back-temp target)
│   ├── UMD_Wood_CAPA_N2_50kW_R{1,2,3,4}.csv   ✅ 4 replicates (R4 ~95% NaN)
│   ├── UMD_Wood_CAPA_N2_70kW_R{1,2,3}.csv     ✅ 3 replicates (validation only)
│   ├── MaCFP_White_Pine_CAPA_30kW_260404_R{4,5,6}_detailed.csv  ✅ BC source (Inner/Outer wall T)
│   ├── MaCFP_White_Pine_CAPA_50kW_260404_R{3,7,8}_detailed.csv  ✅ BC source
│   ├── MaCFP_White_Pine_CAPA_70kW_260404_R{1,2,9}_detailed.csv  ✅ BC source
│   └── bcs_<flux>kW.csv                       ⏳ optional override (top-gas T-series will go here)
└── results/
    ├── trace.csv                   per-evaluation log (auto-rolled per run)
    ├── grid_scan.csv               5×5 residual surface (Stage 1 only)
    ├── best.json                   best θ + std error + correlation matrix
    └── fit_plots.png               experimental vs simulated traces
```

---

## Component schema (NC = 14)

The white pine JSON adds an inert **ash** component (0.4% by mass, same ρ/cp/ε/k
as char). The full ordering used everywhere in the code:

```
1 = Virgin (WP)        7 = Ash (inert)        12 = WP_gas4
2 = WP_int1            8 = Moisture           13 = WP_gas5
3 = WP_int2            9 = WP_gas1            14 = WATER_VAPOR
4 = WP_int3           10 = WP_gas2
5 = WP_int4           11 = WP_gas3
6 = Char
```

Ash sits in the char matrix throughout the simulation and shares the char's
conductivity parameter (so optimizing `k_char` automatically updates ash's `k`).

Initial mass concentrations: `0.996 × 360 = 358.56 kg/m³` virgin + `0.004 ×
360 = 1.44 kg/m³` ash. Moisture concentration is zero per user spec.

---

## Forward model parameterization

Three escalating parameter schemas, switchable via `--stage=1|2|3`. Always
start at Stage 1 — only escalate if residuals demand it.

| Stage | θ                                              | Dim | Intermediate `k_i` rule |
|-------|------------------------------------------------|-----|--------------------------|
| 1     | `[k_virgin, k_char]`                           | 2   | `k_virgin + (k_char−k_virgin)·progress_i` (cumulative-yield) |
| 2     | `[k_virgin, k_int1..4, k_char]`                | 6   | independent per intermediate |
| 3     | `[k_v_lo, Δk_v, k_c_lo, Δk_c]`                 | 4   | Each phase's `k(T)` linear, parameterized over its physically active window: virgin on `[298, 500] K` (HFM/LFA prior range), char on `[500, 1100] K` (pyrolysis-onset to peak surface T). `Δk_x = k_x(T_HI) − k_x(T_LO)` per phase. Intermediates: `k_i(T) = (1−p_i)·k_v(T) + p_i·k_c(T)` — convex combination of the two phase lines by yield-progress. `Δk` lower bounds slightly negative so init=0 reproduces Stage 1 from the interior of the Fminbox box |

`progress_i` derived from white pine cumulative solid yields:
```
ν = (0.914206, 0.785307, 0.386929, 0.778861, 0.697933)
final char yield = 0.151004  (matches JSON exactly)
```

If Stage 1 fits poorly, escalate to **Stage 3 first** (more physically
interpretable, fewer parameters), then Stage 2 only if Stage 3 also fails.

---

## Boundary conditions (final values)

```julia
# top — exposed to cone heater
top_thermal = RadiativeFluxBC(flux = q_ext) +                # uses material ε
              RadiativeBC(T_inf = T_top_surr) +              # uses material ε
              ConvectiveBC(h = 8.6, T_inf = T_top_gas)       # h_top = 8.6 W/(m²·K)
top_mass    = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

# bottom — free to surroundings, impermeable
bottom_thermal = RadiativeBC(ε = 0.92, T_inf = T_bot_surr) + # ε_bot pinned
                 ConvectiveBC(h = 4.0, T_inf = T_bot_gas)    # h_bot = 4.0 W/(m²·K)
bottom_mass    = ImpermeableBC()
```

- **Top emissivity is computed from components** (decomposes virgin → char
  through the run), so `RadiativeFluxBC` and front `RadiativeBC` both pass
  `ε = nothing` and pick up the effective material emissivity.
- **Bottom emissivity is fixed at 0.92** (high-emissivity coating / measurement
  target on back face), independent of decomposition state.
- **`T_inf` callables** (`T_top_gas`, `T_top_surr`, `T_bot_gas`, `T_bot_surr`)
  are resolved per-flux in this order, with later sources overriding earlier:
  1. **Analytic top-gas model** for `T_top_gas` only (user-supplied
     two-exponential fit, see "Top-gas analytic model" below). Coefficients
     interpolate linearly in nominal heat flux.
  2. **MaCFP detailed CSVs** (replicate-averaged across the matching
     `MaCFP_White_Pine_CAPA_<flux>kW_*_detailed.csv` files):
     - `Inner Wall Temperature (K)` → `T_bot_gas` AND `T_bot_surr`
     - `Outer Wall Temperature (K)` → `T_top_surr`
  3. **`bcs_<flux>kW.csv`** if present: any subset of the four columns
     (`T_top_gas`, `T_top_surr`, `T_bot_gas`, `T_bot_surr`) overrides the
     corresponding default — drop-in slot for measured replacements.
  4. **Constant `T_AMBIENT_DEFAULT = 300 K`** for any series still unset
     (warning emitted; with the current setup all four are resolved).

### Top-gas analytic model

`T_top_gas(t; q̇") = T1·exp(T2·t) + T3·exp(T4·t) + T_HFG`, with coefficients
tabulated against nominal flux in `wp_objective.jl::T_TOP_GAS_TABLE` for
{30, 40, 50, 60, 70, 80} kW/m². Linearly interpolated in `q̇"` between rows;
flat extrapolation outside [30, 80]. The 50 and 70 kW/m² rows are themselves
linear interpolations of (40, 60) and (60, 80) respectively, matching the
user's interpolated profile spec exactly. **Note**: the model takes the
**nominal** flux (30, 50, 70) — `HEAT_FLUX_CAL` is not applied here, since the
fit was made against nominal. The cone radiative input `q_ext` still uses the
calibrated value.

---

## Objective function

```
loss(θ) = Σ_q  w_q · MSE_normalized_q  +  w_prior · MSE_prior
MSE_normalized_q = (1/N_q) · Σ_i  ((T̂_back(t_i; θ, q) − T_back^exp(t_i; q)) / σ_q,i)²
σ_q,i  = max( std_replicates_q,i , 0.5 K )
MSE_prior = (1/N_prior) · Σ_j ((k_model(T_j; θ) − k_meas,j) / σ_eff,j)²
```

Implementation:
- **Replicate averaging**: for each flux, all `UMD_Wood_CAPA_N2_<flux>kW_R*.csv`
  are auto-discovered, NaN rows dropped per file, and at each timestep the mean
  is taken across whichever replicates are still valid. Replicate times are
  rounded to integer seconds (1 Hz grid) for robust dict keying.
- **Resampling**: linear interpolation of simulated `T_back` (on
  `solution.times`) onto experimental times. Never the other direction.
- **Per-point σ**: between-replicate population std with a 0.5 K floor (TC
  noise lower bound). At single-replicate timestamps σ falls back to the
  floor.
- **Per-flux weighting** equal (`w_q = 0.5` each) by default.
- **Parallel scenarios**: each loss call dispatches the two fit scenarios via
  `Threads.@spawn` when `Threads.nthreads() > 1` (run with `julia -t auto`).
  Single-thread mode falls back to a sequential loop to avoid spawn overhead.
- **k(T) priors (HFM/LFA, parallel)**: Gaussian penalty on `k_v(T)` using per-
  point uncertainties from the CSVs. Effective prior σ combines the reported σ
  with dataset-level floors in `wp_objective.jl`:
  `σ_eff = sqrt((σ_meas·scale)^2 + max(σ_abs, σ_rel·|k|)^2)`. Weight is
  controlled by the `--prior-weight=W` CLI flag, `W ∈ [0, 1]` (default 0.5):
  0 disables the prior, 0.5 weights it equal to one flux scenario, 1.0
  weights it equal to the combined data term.
- **Time window**: skip first 1 s (`t_skip=1.0`); tail follows experimental
  CSV duration.
- **Failure handling**: non-`Success` retcode or NaN/Inf in back T → that
  scenario contributes `FAILURE_PENALTY = 1e6 × w_q` to the loss. Box-bound
  enforcement is delegated entirely to `Fminbox`; no in-loss bounds penalty.
- **Logging**: every evaluation appends one row to `results/trace.csv` with
  columns `iter,timestamp,stage,theta,loss,data_loss,prior_loss,wallclock_s,
  success`. Each fresh run rolls the existing trace to a timestamped file.

---

## Optimizer

**Optim.jl `Fminbox(NelderMead())`** with multi-start, Hessian-based covariance
(default) or residual bootstrap (`--bootstrap=N`).

- 2–6 parameters; derivative-free inner solver — robust to the loss-surface
  noise from solver tolerances + linear interpolation onto experimental times.
  No FD gradient evaluations (LBFGS+FD made every line search probe-noisy and
  burned ~5 evals/step before).
- Bounds enforced via `Fminbox`. The loss only returns a `FAILURE_PENALTY` for
  failed forward solves, never for bound violations.
- **Multi-start initial guesses** (`n_starts=3` default):
    1. nominal `θ_init` from `wp_properties.jl`,
    2. up to `k=2` lowest-loss points from the Stage 1 grid scan
       (`results/grid_scan.csv`) when available,
    3. remaining slots filled by a Sobol sequence over the box.

After convergence, parameter covariance for the σ-normalized weighted MSE
loss `L(θ) = (1/N)·Σ((r_i/σ_i)²)`:

```
Σ_θ̂ ≈ (2 · L(θ̂) / (N − p)) · H_L⁻¹
```

Derived from Gauss-Newton (H_L ≈ (2/N)·JᵀWJ) and the standard weighted-NLS
result `Σ = (JᵀWJ)⁻¹·ŝ²` with `ŝ² = N·L(θ̂)/(N−p)`. `H_L` is computed by
`FiniteDiff.finite_difference_hessian` central differences. Note that this
assumes the σ_i correctly capture noise — `χ²_red ≈ L(θ̂)` is a sanity check
(should be ~1; >1 means σ underestimated and σ_θ correspondingly underestimated).

The `--bootstrap=N` flag opts into a residual bootstrap that resamples per-
time residuals with replacement, refits N times, and reports the empirical
covariance + median. Slower (N × refit cost) but doesn't assume Gauss-Newton
or a noise-free Hessian. Best.json records `covariance_method` so downstream
consumers know which was used.

For Stage 1, a 5×5 grid scan runs first → free identifiability check + warm
start (top-2 minima seed the multi-start list).

---

## Outstanding — what I (Claude) need to do next

These are pending the next "go" signal from the user; nothing is blocked
on the user *today*:

1. **First fit run + diagnosis** of Stage 1 result quality. Will bring back:
   plots, recovered (k_virgin, k_char) ± std, validation residual at 70 kW/m².
2. **Stage 3 escalation** if Stage 1 validation is poor.
3. **Front-`RadiativeFluxBC` absorptivity sensitivity check** — currently
   tracks effective material emissivity (decomposes during run). Worth a
   one-eval comparison against fixing α at the virgin value, since this
   choice affects the early-time front-side energy balance.
4. **Optional** richer reporting in `best.json`: solve diagnostics (final
   thickness, peak surface T, peak MLR per scenario) so the writeup has
   secondary QoIs without rerunning.

---

## Outstanding — what you (user) need to provide / decide

| Item | Status | Impact |
|---|---|---|
| **Inner / Outer wall T** time-series | ✅ provided | Loaded from `MaCFP_White_Pine_CAPA_<flux>kW_*_detailed.csv` and replicate-averaged. Inner Wall → `T_bot_gas`/`T_bot_surr`; Outer Wall → `T_top_surr`. |
| **Top-surface gas T** | ✅ provided (analytic) | User-supplied two-exponential fit `T(t; q̇") = T1·exp(T2·t) + T3·exp(T4·t) + T_HFG` with table coefficients at {30, 40, 50, 60, 70, 80} kW/m², linear interpolation between. Coefficients may be updated; if a measured time-series becomes available later it can be dropped in via `bcs_<flux>kW.csv` to override. |
| **Optim.jl** added to env | ⏳ user step | `]add Optim` in the Pyrolysis.jl env (or use a separate examples project). |
| **TC mounting** clarification | ⏳ pending | CSV header is `TC back 1 (K)` (thermocouple, not IR as originally stated). If the TC is bonded to the back surface, current setup is correct. If embedded at depth d, we'd need to compare to the simulated T at depth d rather than the back surface. |
| **Whether to fit `h_bot`** as a 3rd parameter | open | If Stage 1 shows systematic long-time bias, we may want to admit `h_bot` (currently fixed at 4.0) into θ. Defer until first fit is in. |

---

## Properties to verify (currently using placeholders / defaults)

These have no measured WP value and use either DF defaults or simple
interpolation. They are unlikely to dominate the back-temperature fit but worth
a sanity pass before publication:

| Property | Current value | Source | Sensitivity |
|---|---|---|---|
| Virgin emissivity | 0.81 | DF default (softwood) | Low–moderate (front-side balance) |
| Char emissivity | 0.70 | DF default | Low–moderate (late-time front T) |
| Ash emissivity | 0.70 | = char (assumed) | Negligible (0.4% mass) |
| Pyrolysis gas M, c, k, λ (gas1–gas5) | DF speciation | DF placeholder | Low (only affects in-depth gas energy balance) |
| Solid-phase λ (gas transport through pores) | DF defaults (5e-7 → 1e-5) | DF placeholder | Low–moderate (late-time mass loss) |
| Permeability κ (Darcy flow) | DF defaults (1e-15 → 1e-12) | DF placeholder | Very low for back T |
| Convective h_top | 8.6 W/(m²·K) | user-provided | Moderate (front-side balance) |
| Convective h_bot | 4.0 W/(m²·K) | user-provided | **High** (long-time back T) |
| Bottom ε | 0.92 | user-provided, fixed | **High** (back-side radiation balance) |
| Initial T (uniform) | T_bot_surr(0) | scenario-driven | Low |

The two **High-sensitivity** items are user-provided. The user-provided
front-side numbers (`h_top`, top emissivity tracking) are also locked in.
Everything else is either physically minor or improvable later.

---

## Pitfalls (mitigations in place)

1. **Identifiability of (k_virgin, k_char)** — likely OK because 30 kW/m²
   exercises virgin-dominated regime longer while 50 kW/m² drives deeper
   conversion. The 5×5 grid scan visually confirms this before committing
   to optimization.
2. **Confounding with kinetics, ρ, cp** — held fixed at calibrated WP values.
3. **Back BC characterization** — h_bot, ε_bot are user-supplied; remaining
   risk is the BC time-series (currently constant 300 K placeholder).
4. **Solver-noise vs gradient methods** — handled via Fminbox(LBFGS); if noise
   dominates, revert to a derivative-free method or tighten solver tolerances.
5. **Pre-ignition transients** — first 1 s skipped in residual (`t_skip=1.0`).
6. **T-dependent k** — addressed in Stage 3 if needed.
7. **Anisotropy** — recovered `k` is transverse-to-grain effective value
   (cone calorimeter geometry). Documented in output.
8. **Multi-start** — 3 random starts within bounds + nominal at start #1.
9. **Validation discipline** — 70 kW/m² evaluated only after optimizer
   terminates.
10. **ALE thickness threshold** — `min_thickness=0.5e-3`; 8 mm samples have
    plenty of margin.
11. **Memory** — only `(θ, loss, retcode, t)` logged per evaluation; full
    solutions saved only for the final best parameters.
12. **Replicate quality** — 50 kW/m² R4 is ~95% NaN (failed test); loader
    automatically falls back to remaining 3 replicates after t≈43 s.
    Between-replicate std reported alongside the average.

---

## Confirmed inputs (from user, locked in)

- 5-step kinetic mechanism, calibrated WP values from
  `data/macfp_white_pine_property_file.json` (A·log10 → s⁻¹, E kJ→J, h kJ→J)
- Fixed dilatometry densities: ρ_virgin = 360, ρ_int1..4 = (321, 259, 122, 90),
  ρ_char = 95, ρ_ash = 95 kg/m³
- Virgin cp 3-segment piecewise linear (HFM/bridge/STA from JSON), char cp
  linear from JSON
- Composition: 99.6% virgin wood, 0.4% inert ash
- Same 5 pyrolysis-gas + water-vapor speciation as DF (placeholder for now)
- Initial moisture = 0
- Sample thickness 8 mm
- Back BC: free to surroundings (radiation + convection), impermeable
- h_top = 8.6 W/(m²·K), h_bot = 4.0 W/(m²·K)
- Top emissivity = effective material (decomposing); bottom emissivity = 0.92 fixed
- 30 + 50 kW/m² for fitting; 70 kW/m² for validation
- **Heat flux calibration**: actual delivered cone flux is **97.3% of nominal**
  (set in `inverse_white_pine.jl` as `HEAT_FLUX_CAL = 0.973`); applied when
  constructing each scenario, so the model sees 29.19 / 48.65 / 68.11 kW/m²
- No optimization-budget constraint
- Hessian-based covariance estimates desired

---

## Execution order (updated)

1. ✅ Stand up `wp_forward.jl` with white pine baseline; verify smoke test.
2. ✅ Wire in MaCFP detailed CSVs as Inner/Outer wall BC source
   (replicate-averaged).
3. ✅ Wire in user's analytic two-exponential `T_top_gas(t; q̇")` model.
4. ⏳ **Run Stage 1 first fit** with full BC stack
   (all four BC series now time-varying).
   Report (k_virgin, k_char), validation overlay, gut-check on residual.
5. ⏳ When user updates top-gas coefficients (or supplies a measured
   `bcs_<flux>kW.csv`), rerun and compare.
6. ⏳ If validation at 70 kW/m² is poor → **Stage 3** (T-dependent linear).
7. ⏳ If Stage 3 still poor → **Stage 2** (per-intermediate independent).
8. ⏳ Final reporting: best θ + std errors + correlation matrix, residual
   plots, validation overlay, copy of trace.csv.

---

## Dependency notes

- `Optim.jl` — needs `]add Optim` (not in main package deps).
- `FiniteDiff.jl` — already in `Project.toml`.
- `Plots.jl` — already in `Project.toml`.
- No JSON parser dependency: WP properties are translated into `wp_properties.jl`
  by hand (ten-line block, easy to re-sync if the JSON changes).
