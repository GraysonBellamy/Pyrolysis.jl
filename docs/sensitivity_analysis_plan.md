# Sensitivity Analysis Implementation Plan for ThermaKin3

## Context

We need a sensitivity analysis framework for the ThermaKin3 pyrolysis simulation package. The goal is to propagate uncertainties in thermophysical input parameters (heat capacities, heat of sorption, heats of reaction) through the douglas fir cone calorimeter simulation and quantify their impact on key outputs (MLR curve, surface/back temperature, scalar summaries).

Each simulation takes ~6-8 minutes. We will use **direct Latin Hypercube Sampling** with **500-1000 samples** run in parallel via `Distributed.jl`. No surrogate model is needed at this simulation cost. Results will be analyzed with `GlobalSensitivity.jl` for Sobol indices.

---

## Uncertain Parameters (16 total)

### Heat Capacity Regression Coefficients — DF Virgin + Char (8 params)

The intermediate solid cp's (DF_int1 through DF_int4) are **not independent** — they are
linearly interpolated between virgin DF and char cp based on cumulative solid yield through
the reaction chain. Therefore only the DF and char cp coefficients are uncertain; the
intermediates are derived.

**Interpolation formula:**
```
final_char_yield = 0.95 × 0.80 × 0.37 × 0.79 × 0.74 = 0.1644
progress_i = (1 - cumulative_solid_yield_i) / (1 - final_char_yield)
cp_int_i(T) = (1 - progress_i) × cp_DF(T) + progress_i × cp_char(T)
```

**Cumulative solid yields and interpolation weights:**
| Intermediate | Cumulative Yield | Progress (α) |
|-------------|-----------------|---------------|
| DF_int1 | 0.95 | 0.0598 |
| DF_int2 | 0.76 | 0.2872 |
| DF_int3 | 0.2812 | 0.8602 |
| DF_int4 | 0.2222 | 0.9309 |

**Uncertain cp parameters:**

| # | Parameter | Nominal | Component | Form |
|---|-----------|---------|-----------|------|
| 1 | `df_cp_a1` | -53.05 | DF (T≤323) | intercept |
| 2 | `df_cp_b1` | 4.245 | DF (T≤323) | slope |
| 3 | `df_cp_a2` | -31.03 | DF (323<T≤425) | intercept |
| 4 | `df_cp_b2` | 4.177 | DF (323<T≤425) | slope |
| 5 | `df_cp_a3` | -2131.0 | DF (T>425) | intercept |
| 6 | `df_cp_b3` | 9.119 | DF (T>425) | slope |
| 7 | `char_cp_a` | 347.2 | Char | intercept |
| 8 | `char_cp_b` | 1.784 | Char | slope |

### Heat of Sorption Equation Parameters (3 params)

| # | Parameter | Nominal | Description |
|---|-----------|---------|-------------|
| 9 | `hos_baseline` | 2275.0 | Baseline HoS [kJ/kg] |
| 10 | `hos_amplitude` | 860.0 | Exponential amplitude [kJ/kg] |
| 11 | `hos_decay` | 0.120 | Exponential decay rate [1/%MC] |

### Heats of Reaction (5 params)

| # | Parameter | Nominal | Reaction |
|---|-----------|---------|----------|
| 12 | `h_rxn1` | -20067.0 | pyrolysis_step1 [J/kg] |
| 13 | `h_rxn2` | -8264.0 | pyrolysis_step2 [J/kg] |
| 14 | `h_rxn3` | 273999.0 | pyrolysis_step3 [J/kg] |
| 15 | `h_rxn4` | -224430.0 | pyrolysis_step4 [J/kg] |
| 16 | `h_rxn5` | 555.6 | pyrolysis_step5 [J/kg] |

**Total: 16 uncertain parameters** (8 intermediate cp coefficients are derived, not independent)

---

## Outputs to Capture Per Simulation

### Scalar QoIs
- Peak MLR [kg/(m²·s)]
- Time to peak MLR [s]
- Time to ignition (first time MLR > threshold, e.g. 3e-3 kg/(m²·s)) [s]
- Total mass loss fraction [%]
- Peak surface temperature [K]
- Peak back-face temperature [K]
- Final thickness [m] (from ALE)

### Time-Series Curves (stored per run)
- `surface_mass_flux` vs `times` (MLR curve)
- `surface_temperature` vs `times`
- `back_temperature` vs `times`

---

## File Structure

Create a single new example/analysis script:

```
examples/sensitivity_analysis/
├── sensitivity_douglas_fir.jl      # Main driver script
└── results/                        # Output directory (created at runtime)
    ├── parameters.csv              # All sampled parameter sets
    ├── scalar_qois.csv             # Scalar QoIs for each sample
    ├── mlr_curves/                 # One CSV per run: time, mlr
    │   ├── run_0001.csv
    │   └── ...
    ├── surface_temp_curves/        # One CSV per run: time, T_surface
    │   ├── run_0001.csv
    │   └── ...
    ├── back_temp_curves/           # One CSV per run: time, T_back
    │   ├── run_0001.csv
    │   └── ...
    ├── sobol_indices.csv           # First-order and total Sobol indices
    └── summary_plots.png           # Envelope plots, histograms, tornado chart
```

All in a standalone script — **no modifications to the core Pyrolysis.jl package**.

---

## Implementation Steps

### Step 1: Define the Parameter Space

```julia
using Distributions

# Parameter definition struct
struct UncertainParameter
    name::Symbol
    nominal::Float64
    distribution::Distribution  # From Distributions.jl
end

# User specifies ±X as 95% CI → Normal(nominal, X/1.96)
# Or can use any Distributions.jl distribution (Uniform, LogNormal, etc.)
function param_from_ci(name, nominal, half_width_95)
    σ = half_width_95 / 1.96
    return UncertainParameter(name, nominal, Normal(nominal, σ))
end

# Define all 16 parameters with their uncertainties
# Example (user fills in actual ± values):
params = [
    # DF heat capacity (piecewise linear, 3 segments)
    param_from_ci(:df_cp_a1, -53.05, 5.0),    # ±5.0 at 95% CI
    param_from_ci(:df_cp_b1, 4.245, 0.4),
    param_from_ci(:df_cp_a2, -31.03, 3.0),
    param_from_ci(:df_cp_b2, 4.177, 0.4),
    param_from_ci(:df_cp_a3, -2131.0, 200.0),
    param_from_ci(:df_cp_b3, 9.119, 0.9),
    # Char heat capacity (linear)
    param_from_ci(:char_cp_a, 347.2, 35.0),
    param_from_ci(:char_cp_b, 1.784, 0.18),
    # Heat of sorption equation
    param_from_ci(:hos_baseline, 2275.0, 227.0),
    param_from_ci(:hos_amplitude, 860.0, 86.0),
    param_from_ci(:hos_decay, 0.120, 0.012),
    # Heats of reaction
    param_from_ci(:h_rxn1, -20067.0, 2007.0),
    param_from_ci(:h_rxn2, -8264.0, 826.0),
    param_from_ci(:h_rxn3, 273999.0, 27400.0),
    param_from_ci(:h_rxn4, -224430.0, 22443.0),
    param_from_ci(:h_rxn5, 555.6, 55.6),
]
```

### Step 2: Latin Hypercube Sampling

```julia
using QuasiMonteCarlo

n_samples = 500  # User-adjustable
n_params = length(params)

# Generate LHS samples in [0,1]^d
lhs_unit = QuasiMonteCarlo.sample(n_samples, zeros(n_params), ones(n_params),
                                   LatinHypercubeSample())

# Transform to physical space using inverse CDF (quantile function)
samples = Matrix{Float64}(undef, n_params, n_samples)
for i in 1:n_params
    for j in 1:n_samples
        samples[i, j] = quantile(params[i].distribution, lhs_unit[i, j])
    end
end
```

### Step 3: Build the Simulation Function

This is the core piece — a function that takes a parameter vector and returns QoIs.

```julia
function build_and_solve(θ::Vector{Float64})
    # Unpack parameter vector (order matches params definition — 16 total)
    df_cp_a1, df_cp_b1, df_cp_a2, df_cp_b2, df_cp_a3, df_cp_b3,
    char_cp_a, char_cp_b,
    hos_baseline, hos_amplitude, hos_decay,
    h_rxn1, h_rxn2, h_rxn3, h_rxn4, h_rxn5 = θ

    # --- Build parametric property functions ---
    df_cp = CallableProperty(T -> begin
        if T <= 323
            df_cp_a1 + df_cp_b1 * T
        elseif T <= 425
            df_cp_a2 + df_cp_b2 * T
        else
            df_cp_a3 + df_cp_b3 * T
        end
    end)

    char_cp = CallableProperty(T -> char_cp_a + char_cp_b * T)

    # --- Derive intermediate cp's via linear interpolation ---
    # Intermediates are interpolated between virgin DF and char based on
    # cumulative solid yield through the reaction chain:
    #   progress_i = (1 - cumulative_yield_i) / (1 - final_char_yield)
    #   cp_int_i(T) = (1 - progress_i) * cp_DF(T) + progress_i * cp_char(T)
    #
    # Cumulative solid yields:  int1=0.95, int2=0.76, int3=0.2812, int4=0.2222
    # Final char yield: 0.95 * 0.80 * 0.37 * 0.79 * 0.74 = 0.16439
    const FINAL_CHAR_YIELD = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
    progress = [(1 - y) / (1 - FINAL_CHAR_YIELD)
                for y in [0.95, 0.95*0.80, 0.95*0.80*0.37, 0.95*0.80*0.37*0.79]]

    function make_interp_cp(α)
        CallableProperty(T -> (1 - α) * df_cp(T) + α * char_cp(T))
    end

    int1_cp = make_interp_cp(progress[1])  # α ≈ 0.0598
    int2_cp = make_interp_cp(progress[2])  # α ≈ 0.2872
    int3_cp = make_interp_cp(progress[3])  # α ≈ 0.8602
    int4_cp = make_interp_cp(progress[4])  # α ≈ 0.9309

    hos = StateDependentProperty((T, ξ) -> begin
        MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS
        return (hos_baseline + hos_amplitude * exp(-hos_decay * MC)) * 1e3
    end)

    # --- Build components (same as 07_douglas_fir_cone.jl but with parametric cp) ---
    # All other properties (ρ, k, ε, λ, κ, etc.) are FIXED at nominal values
    components = (
        SolidComponent(:DF, ρ=530.0, c=df_cp, k=0.145, ...),        # uses parametric df_cp
        SolidComponent(:DF_int1, ρ=507.2, c=int1_cp, k=..., ...),   # uses parametric int1_cp
        SolidComponent(:DF_int2, ..., c=int2_cp, ...),
        SolidComponent(:DF_int3, ..., c=int3_cp, ...),
        SolidComponent(:DF_int4, ..., c=int4_cp, ...),
        SolidComponent(:Char, ..., c=char_cp, ...),
        LiquidComponent(:MOISTURE, ...),                              # cp FIXED (excluded)
        GaseousComponent(:DF_gas1, ...), ...,                         # all FIXED
        GaseousComponent(:WATER_VAPOR, ...),
    )

    # --- Build reactions (same structure, parametric h values and HoS) ---
    reactions = (
        Reaction(:pyrolysis_step1, 1 => (2, 8), A=74816.95, E=71.59e3,
                 h=h_rxn1, n=1.0, yields=(0.95, 0.05)),
        Reaction(:pyrolysis_step2, 2 => (3, 9), A=23550492839.0, E=140.9e3,
                 h=h_rxn2, n=1.0, yields=(0.80, 0.20)),
        Reaction(:pyrolysis_step3, 3 => (4, 10), A=7.744618e12, E=181.0e3,
                 h=h_rxn3, n=1.0, yields=(0.37, 0.63)),
        Reaction(:pyrolysis_step4, 4 => (5, 11), A=103276.14, E=98.94e3,
                 h=h_rxn4, n=1.0, yields=(0.79, 0.21)),
        Reaction(:pyrolysis_step5, 5 => (6, 12), A=6.397348, E=56.01e3,
                 h=h_rxn5, n=1.0, yields=(0.74, 0.26)),
        Reaction(:moisture_evaporation, 7 => 13, A=2.997e5, E=5.0e4,
                 h=hos, n=1.0),
    )

    # --- Assemble material ---
    material = Material(
        name = :DouglasFir,
        components = components,
        reactions = reactions,
        conductivity_mixing = WEIGHTED, β = 0.5,
        gas_transfer_mixing = WEIGHTED, β_λ = 0.5,
        permeability_mixing = WEIGHTED, β_κ = 0.5,
    )

    # --- Build problem (mesh, BCs, ICs all FIXED at nominal) ---
    # Copy the mesh, BC, IC setup from 07_douglas_fir_cone.jl verbatim
    problem = PyrolysisProblem(
        mesh = ...,
        material = material,
        bc_set = ...,
        T_initial = ...,
        ξ_initial = ...,
        tspan = (0.0, 600.0),
        output_times = collect(range(0.0, 600.0, length=3601)),
        ...
    )

    # --- Solve ---
    solution = solve(problem;
        abstol = 1e-9, reltol = 1e-7,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        use_analytical_jacobian = true,
        solver = KenCarp4(linsolve = KLUFactorization()),
        progress = false,   # disable progress bars for batch runs
    )

    return solution
end
```

### Step 4: Extract QoIs from Solution

```julia
struct SensitivityResult
    # Scalars
    peak_mlr::Float64
    time_to_peak_mlr::Float64
    time_to_ignition::Float64    # first time MLR > threshold
    total_mass_loss_pct::Float64
    peak_surface_temp::Float64
    peak_back_temp::Float64
    final_thickness::Float64
    retcode::Symbol
    solve_time::Float64

    # Curves (stored as vectors aligned with common time grid)
    times::Vector{Float64}
    mlr_curve::Vector{Float64}
    surface_temp_curve::Vector{Float64}
    back_temp_curve::Vector{Float64}
end

const MLR_IGNITION_THRESHOLD = 3e-3  # kg/(m²·s), adjustable

function extract_qois(sol::PyrolysisSolution)
    t = sol.times
    mlr = sol.surface_mass_flux
    T_surf = sol.surface_temperature
    T_back = sol.back_temperature

    peak_mlr = maximum(mlr)
    peak_idx = argmax(mlr)
    time_to_peak = t[peak_idx]

    tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
    time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf

    total_mass_loss = (1.0 - sol.total_mass[end] / sol.total_mass[1]) * 100.0

    return SensitivityResult(
        peak_mlr, time_to_peak, time_to_ignition,
        total_mass_loss,
        maximum(T_surf), maximum(T_back),
        sol.thickness_history[end],
        sol.retcode, sol.solve_time,
        t, mlr, T_surf, T_back
    )
end
```

### Step 5: Parallel Execution with Distributed.jl

```julia
using Distributed

# Add worker processes (user sets based on available cores)
n_workers = 8
addprocs(n_workers)

# Load package on all workers
@everywhere begin
    using Pyrolysis  # or however the package is loaded
    # Include the build_and_solve and extract_qois functions
    # (either via include() or defining them here)
end

# Run all samples in parallel with checkpointing
results_dir = "results"
mkpath(results_dir)
mkpath(joinpath(results_dir, "mlr_curves"))
mkpath(joinpath(results_dir, "surface_temp_curves"))
mkpath(joinpath(results_dir, "back_temp_curves"))

# Save parameter samples first
# (CSV with columns: sample_id, param1, param2, ..., param24)
save_parameters_csv(joinpath(results_dir, "parameters.csv"), params, samples)

# Parallel map with progress tracking and checkpointing
results = pmap(1:n_samples; retry_delays=zeros(0)) do i
    θ = samples[:, i]

    # Check if this run already completed (for resumability)
    scalar_file = joinpath(results_dir, "scalar_qois.csv")
    # ... checkpointing logic ...

    try
        sol = build_and_solve(θ)
        result = extract_qois(sol)

        # Save curves immediately (incremental output)
        save_curve_csv(joinpath(results_dir, "mlr_curves",
                       @sprintf("run_%04d.csv", i)), result.times, result.mlr_curve)
        save_curve_csv(joinpath(results_dir, "surface_temp_curves",
                       @sprintf("run_%04d.csv", i)), result.times, result.surface_temp_curve)
        save_curve_csv(joinpath(results_dir, "back_temp_curves",
                       @sprintf("run_%04d.csv", i)), result.times, result.back_temp_curve)

        return (i, result)
    catch e
        @warn "Sample $i failed: $e"
        return (i, nothing)  # Record failure, don't crash the batch
    end
end

# Collect scalar QoIs into summary CSV
save_scalar_qois_csv(joinpath(results_dir, "scalar_qois.csv"), results)
```

### Step 6: Sensitivity Analysis (Sobol Indices)

```julia
using GlobalSensitivity

# For regression-based Sobol from LHS data:
# Extract scalar QoI matrix: [n_samples × n_qois]
qoi_matrix = hcat([
    [r.peak_mlr, r.time_to_peak_mlr, r.time_to_ignition,
     r.total_mass_loss_pct, r.peak_surface_temp, r.peak_back_temp]
    for (_, r) in results if r !== nothing
]...)'

param_matrix = samples[:, [i for (i, r) in results if r !== nothing]]'

# Regression-based sensitivity (works directly on LHS samples)
# Uses rank correlation coefficients as sensitivity measure
for (qoi_idx, qoi_name) in enumerate(["peak_mlr", "t_peak_mlr", "t_ignition",
                                       "mass_loss", "peak_T_surf", "peak_T_back"])
    y = qoi_matrix[:, qoi_idx]

    # Standardized Rank Regression Coefficients (SRRC)
    # These approximate first-order Sobol indices for monotonic relationships
    srrc = compute_srrc(param_matrix, y)

    println("$qoi_name sensitivity ranking:")
    sorted = sortperm(abs.(srrc), rev=true)
    for j in sorted
        println("  $(params[j].name): SRRC = $(round(srrc[j], digits=4))")
    end
end
```

**Alternative: For exact Sobol indices**, use the Saltelli sampling scheme:

```julia
# Saltelli scheme for exact Sobol indices
# Requires N*(2d+2) = N*(2*16+2) = 34*N evaluations
# At N=16: 544 simulations; at N=32: 1088 simulations
# More expensive but gives exact first-order and total-order indices

using GlobalSensitivity

lb = [quantile(p.distribution, 0.005) for p in params]  # Lower bounds
ub = [quantile(p.distribution, 0.995) for p in params]  # Upper bounds

# Define model function for GlobalSensitivity.jl interface
function model_for_gsa(param_matrix)
    # param_matrix is [n_samples × n_params]
    n = size(param_matrix, 1)
    qois = Vector{Float64}(undef, n)
    # Run in parallel...
    for i in 1:n
        sol = build_and_solve(param_matrix[i, :])
        result = extract_qois(sol)
        qois[i] = result.peak_mlr  # or whichever scalar QoI
    end
    return qois
end

sobol_result = gsa(model_for_gsa, Sobol(),
                   [[lb[i], ub[i]] for i in 1:n_params],
                   samples = 32)  # N=32 → 1600 model evaluations
```

### Step 7: Visualization

```julia
using Plots

# --- 1. MLR Envelope Plot ---
# Load all MLR curves, compute percentiles at each time point
all_mlr = hcat([r.mlr_curve for (_, r) in results if r !== nothing]...)
t = results[1][2].times

p5  = [quantile(all_mlr[i, :], 0.05) for i in 1:length(t)]
p25 = [quantile(all_mlr[i, :], 0.25) for i in 1:length(t)]
p50 = [quantile(all_mlr[i, :], 0.50) for i in 1:length(t)]
p75 = [quantile(all_mlr[i, :], 0.75) for i in 1:length(t)]
p95 = [quantile(all_mlr[i, :], 0.95) for i in 1:length(t)]

plot(t, p50, label="Median", lw=2, color=:blue)
plot!(t, p5, fillrange=p95, alpha=0.2, label="90% CI", color=:blue)
plot!(t, p25, fillrange=p75, alpha=0.3, label="50% CI", color=:blue)
xlabel!("Time [s]")
ylabel!("MLR [kg/(m²·s)]")
title!("Mass Loss Rate — Uncertainty Envelope")

# --- 2. Scalar QoI Histograms ---
# Histograms of peak_mlr, time_to_ignition, etc.

# --- 3. Tornado Chart (parameter importance) ---
# Bar chart of |SRRC| values sorted by magnitude for each QoI

# --- 4. Surface/Back Temperature Envelopes ---
# Same percentile envelope approach as MLR
```

---

## Key Design Decisions

### Why LHS over Saltelli sampling
- LHS with 500 samples + regression-based sensitivity (SRRC) gives both **output distributions** and **parameter importance rankings** in one batch
- Saltelli would require ~1600 samples for exact Sobol but only gives sensitivity, not full output distributions
- SRRC is a well-accepted approximation for monotonic input-output relationships, which thermophysical parameters typically have
- If exact Sobol indices are needed later, the script supports switching to Saltelli scheme

### Checkpointing strategy
- Each completed simulation saves its curves to disk immediately
- A master CSV accumulates scalar QoIs row by row
- If the batch is interrupted, completed runs are detected and skipped on restart
- This is critical for multi-day runs — no work is lost

### Error handling
- Failed simulations (solver divergence, NaN, etc.) are logged but don't crash the batch
- The `retcode` field from each solution is saved to identify problematic parameter combinations
- Summary statistics exclude failed runs

### Distributed.jl approach
- Each worker gets its own Julia process with full Pyrolysis.jl loaded
- Workers are stateless — each `pmap` call builds the full material/problem from scratch
- No shared mutable state, no thread-safety concerns
- `pmap` handles load balancing automatically

---

## Dependencies to Add

```julia
# In the script (not the core package Project.toml):
using Pkg
Pkg.add(["Distributions", "QuasiMonteCarlo", "GlobalSensitivity", "CSV", "DataFrames", "Printf"])
```

| Package | Purpose |
|---------|---------|
| `Distributions.jl` | Define parameter uncertainty distributions |
| `QuasiMonteCarlo.jl` | Latin Hypercube Sampling |
| `GlobalSensitivity.jl` | Sobol indices, sensitivity analysis |
| `CSV.jl` + `DataFrames.jl` | Save/load results |
| `Printf` | Formatted file names |
| `Distributed` | Parallel execution (stdlib) |

---

## Critical Files to Reference During Implementation

| File | Purpose |
|------|---------|
| `examples/07_douglas_fir_cone.jl` | Base simulation to parameterize — copy the full component/reaction/BC/IC setup |
| `src/Core/property_functions.jl` | `CallableProperty`, `StateDependentProperty` constructors |
| `src/Core/components.jl` | `SolidComponent`, `GaseousComponent`, `LiquidComponent` constructors |
| `src/Core/reactions.jl` | `Reaction` constructor with yields syntax |
| `src/Core/materials.jl` | `Material` constructor, mixing rules |
| `src/Problem/pyrolysis_problem.jl` | `PyrolysisProblem` constructor, `PyrolysisSolution` fields |
| `src/Solvers/solve.jl` | `solve()` function signature and keyword arguments |

---

## Verification

1. **Nominal validation**: Run with all parameters at nominal values and compare against the existing `07_douglas_fir_cone.jl` output — results should be identical
2. **Single-parameter sweep**: Perturb one parameter (e.g., `h_rxn3`) at ±10% with 5 values to verify the script produces sensible, monotonic trends
3. **Small batch test**: Run 10-20 LHS samples to verify parallelism, checkpointing, and CSV output before committing to a full 500-sample run
4. **Sobol sanity check**: Parameters known to have large effects (like the large endothermic `h_rxn3 = 273999 J/kg`) should show high sensitivity indices
