# 14. Worked Examples and Tutorials

This chapter walks through the example scripts shipped in the `examples/`
directory of Pyrolysis.jl. Each example is a complete, runnable workflow that
exercises a different combination of physics (radiation models, lumped gas transport
transport, moisture, ALE shrinkage), solver options, and post-processing. For
each one we state what is being modeled, which simulation modes and boundary
conditions are configured, how to run it, what outputs to look at, and which
chapters of the Technical Reference derive the underlying physics. The code
fragments here are drawn from the shipped scripts; where a snippet is
abbreviated, the original file is cited so you can read it in full.

All example code uses the package's coordinate convention: `z = 0` is the back
face (tag `:bottom`) and `z = L` is the exposed surface (tag `:top`); heat enters
from the top and the material recedes downward. SI units are used throughout
(see Technical Reference §2).

The examples live in `examples/` with their own `Project.toml` environment that
adds `Plots` (and a few helper packages) on top of Pyrolysis. Run them from the
repository root:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=examples examples/06_cone_calorimeter.jl
```

or from a Julia session started in the repository root:

```julia
using Pkg
Pkg.activate("examples")
Pkg.instantiate()         # first time only

include("examples/06_cone_calorimeter.jl")
```

Scripts that write files write them into the working directory (the repository
root: `cone_calorimeter_*.png`, `douglas_fir_*`, `ptw_*`), into `output/`, or
into `results/` directories beside the scripts. Every simulation script takes
the environment variable `PYROLYSIS_SMOKE=1`: a short run (10–20 s of simulated
time, and a subset of the cases where a script loops over many) that checks the
script works and writes the same files as a full run:

```bash
PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter.jl
```

`examples/README.md` lists every script with its fixtures, its outputs and its
smoke command.

A roadmap of what each example demonstrates:

| Example file | Physics highlighted | Modes | §14 |
|---|---|---|---|
| `06_cone_calorimeter.jl` | PMMA cone, single-step kinetics, incident flux at the surface (`NO_RADIATION`) | ALE, deactivation | 14.1 |
| `06_cone_calorimeter_surface_absorption.jl` | PMMA cone, surface-only radiation | ALE, SURFACE_ABSORPTION, deactivation | 14.1 |
| `07_douglas_fir_cone.jl` | Wood, multi-step kinetics, lumped gas transport, moisture | ALE, SURFACE_ABSORPTION, LumpedTransfer | 14.2 |
| `convergence_douglas_fir.jl` | Grid-convergence study | as 07 at 50 kW/m², swept over `n_cells` | 14.3 |
| `08_pressure_treated_wood_cone.jl` | Treated wood + moisture, lumped gas transport | ALE, SURFACE_ABSORPTION, LumpedTransfer | 14.4 |
| `08_gas_transport.jl` | Lumped and dusty-gas transport, finite film | Moving mesh, both transport models | 14.5 |
| `moisture_hos_comparison.jl` | State-dependent heat of sorption | ALE, SURFACE_ABSORPTION, LumpedTransfer | 14.6 |
| `macfp/macfp_gasification.jl` | MaCFP-4 white-pine gasification targets (10/30/60 kW/m², 6/18 mm), the calibrated model of §14.7 at a constant cross section | ALE, SURFACE_ABSORPTION, threads | — |
| `inverse_analysis/inverse_white_pine.jl` (with `validate_thickness.jl`, `plot_diagnostics.jl`, `plot_stage2_theta.jl`) | Parameter estimation from data | ALE, optimization | 14.7 |
| `sensitivity_analysis/sensitivity_douglas_fir.jl`, `sensitivity_analysis/sensitivity_sobol.jl`, `sensitivity_analysis/compare_degrees.jl` | Uncertainty / Sobol | parameter sweep, post-process | 14.8 |
| `literature_comparison/literature_comparison.jl` | Property-variant comparison | as 07 at 50 kW/m², swept over property sets | 14.9 |

---

## 14.1 PMMA cone calorimeter (`06_cone_calorimeter.jl`)

The executable example defines virgin polymer, residual char, MMA and inert
air, with temperature-dependent thermal properties and decomposition into
MMA. The condensed phases declare nominal pore volume and superficial
lumped transfer. `PYROLYSIS_PMMA_POROSITY` selects the nominal porosity for
sensitivity studies. The atmospheric reservoir contains air, and the initial
air inventory follows the pore EOS at 300 K.

The top receives radiation, reradiates, and exchanges with a 330 K heat and
mass film. Incoming air carries that reservoir enthalpy. The bottom is
insulated and impermeable. The solve uses ALE, in-depth absorption and
surface deactivation. `06_cone_calorimeter_surface_absorption.jl` uses the
surface-absorption radiation choice with the same explicit gas closure.

```bash
PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter.jl
PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter_surface_absorption.jl
```

Inspect `get_pressure`, `get_porosity` and `get_boundary_exchange` alongside
temperatures, mass loss and thickness. Mesh, nominal porosity and pressure
error budget must be varied jointly when assessing dense-polymer results.

---

## 14.2 Douglas fir cone calorimeter (`07_douglas_fir_cone.jl`)

The wood example includes multi-step pyrolysis, absorbed moisture, water
vapor and inert pore air. Condensed phases provide skeletal densities,
lumped matrix transfer and permeability. The complete atmospheric reservoir
has incoming air at the same temperature as the external thermal film.
The initial air inventory follows the actual pore fraction; zero swelling
weight on absorbed water does not remove its skeletal occupancy.

```bash
PYROLYSIS_SMOKE=1 julia --project=examples examples/07_douglas_fir_cone.jl
```

The moving mesh, radiation and deactivation choices are explicit solve
keywords. Use the species boundary diagnostics to separate combustible-gas
escape from vapor transport and air ingress. Reconstructing flux from the
initial mesh or only the pressure gradient would omit parts of this balance.

---

## 14.3 Grid-convergence study (`convergence_douglas_fir.jl`)

### What it models

The Douglas fir cone case at 50 kW/m² (13 components: the material of §14.2
without the background air), solved at six spatial resolutions
(480, 240, 120, 60, 30, 15 cells → 0.025–0.8 mm cells) to quantify how coarse
the mesh can be before the quantities of interest (peak MLR, time to ignition,
total mass loss, peak surface temperature, solve time) drift from the finest-grid
reference. This is the spatial half of the verification workflow in Technical
Reference §18.

### How it is structured

The whole build-and-solve is wrapped in a function taking `n_cells`, so the same
material and BCs are re-instantiated at each resolution. The material and problem
are rebuilt inside the function (rather than mutating a shared problem) so workers
do not share mutable state:

```julia
@everywhere function build_and_solve(n_cells::Int)
    # ... build material (13 components, 6 reactions) and BCs ...
    mesh = create_uniform_mesh(12.0e-3, n_cells)
    ξ_initial_funcs = [z -> 530.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
                       z -> 31.8,                                 # 6 % moisture, dry basis
                       z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0]
    problem = PyrolysisProblem(mesh = mesh, material = material, bc_set = bc_set,
        T_initial = z -> 300.0, ξ_initial = ξ_initial_funcs,
        tspan = (0.0, 600.0), dt_initial = 1e-3,
        output_times = collect(range(0.0, 600.0, length = 3601)))
    return solve(problem;
        abstol = 1e-9, reltol = 1e-7, dtmin = 1e-14, progress = false,
        use_ale = true, radiation_model = SURFACE_ABSORPTION, min_thickness = 0.5e-3)
end
```

### How to run it

The grid sizes are independent, so the study parallelizes trivially with
`Distributed.pmap`. The script loads packages on all workers and calls
`pmap(build_and_solve, GRID_SIZES)`:

```bash
julia --project=examples examples/convergence_douglas_fir.jl            # sequential (1 worker)
julia -p auto --project=examples examples/convergence_douglas_fir.jl    # parallel (all cores)
julia -p 4 --project=examples examples/convergence_douglas_fir.jl       # 4 workers
```

```julia
using Distributed
@everywhere begin
    using Pyrolysis
    using LinearAlgebra
end
solutions = pmap(build_and_solve, [480, 240, 120, 60, 30, 15])
```

The script writes its tables and figures to `examples/convergence_results/`; its
smoke run solves the two coarsest grids over 10 s.

### What to look at

The script tabulates scalar QoIs and their relative error versus the 480-cell
reference, and computes curve-level `L2` errors for the MLR, surface, and back
temperature traces. The takeaway you want from any convergence study is the
coarsest mesh whose QoI error is below your tolerance (e.g. 1%). Because output
is saved to a fixed `output_times` grid, the curve `L2` comparison is
element-wise across resolutions:

```julia
mlr_l2 = norm(r.mlr_curve .- ref.mlr_curve) / norm(ref.mlr_curve) * 100   # [%]
```

Cross-reference: the same kind of distributed sweep is also available through the
package's `parameter_sweep` helper, with a `modify` that replaces the mesh
(§12.6 and §14.8.3); this script calls `pmap` over its own
`build_and_solve(n_cells)`, which builds the whole problem at each resolution.

---

## 14.4 Pressure-treated wood (`08_pressure_treated_wood_cone.jl`)

### What it models

Cone pyrolysis of fire-retardant pressure-treated wood (PTW) with moisture, at
50 kW/m², 12 mm on 240 cells. Structurally similar to Douglas fir but with a
different component ordering (water and water vapor come first), three solid
intermediates between the virgin wood and the char instead of four, and five
reactions. It exercises the same lumped gas transport + surface-absorption + ALE stack.

### Material highlights

Twelve components: `WATER` (liquid), `WATER_VAPOR` (gas), `PTW`/`PTW1`/`PTW2`/
`PTW3`/`Char` (solids, with skeletal densities `ρ_intrinsic`), `PTW_gas1`…`PTW_gas4`
(gases), and the background `AIR` (gas, initialized to 1 atm in the pores by
`background_gas_concentration`). The liquid water is bound water absorbed into
the cell walls, not a separate phase:

```julia
LiquidComponent(:WATER;
    ρ = 1e5,           # pseudo-density [kg/m³]
    c = water_cp,
    k = 0.145,
    ε = 0.81,
    α = 1000.0,        # not used with SURFACE_ABSORPTION
    κ = 0.0,           # liquid blocks Darcy flow
)
```

A liquid's default `γ=0` removes its direct strain and matrix transport
weight. Its mass, heat capacity, evaporation energy and skeletal occupancy
remain. The imported pseudo-density is a constitutive input to the latter;
it does not prove an unobstructed void. Set `γ=1` for a matrix-forming melt
and provide its transport coefficients explicitly.

### Reaction heats and rates

Pyrolysis.jl uses `h > 0` for endothermic reactions and `h < 0` for
exothermic reactions. Water evaporation and the first decomposition step
absorb heat; the third decomposition step releases it:

```julia
Reaction(:water_evaporation, 1 => 2; A = 15532.2434, E = 43528.868,
         h = 2.78e6, n = 1.0),                     # evaporation: endothermic
Reaction(:ptw_step1, 3 => (4, 8); A = 19175814.8751, E = 106116.6114,
         h = 6.82e3, n = 1.0, yields = (0.73075, 0.26925)),    # endothermic
Reaction(:ptw_step3, 5 => (6, 10); A = 179027.4062, E = 100963.197,
         h = -1.82e5, n = 1.0, yields = (0.81744, 0.18256)),   # exothermic
```

These reactions use Arrhenius rates without explicit temperature cutoffs.
`Reaction` also supports smooth `tanh` gates through `T_min`, `T_max` and
`gate_width` (Technical Reference §6).

### Solve and outputs

The solve call is the same as Douglas fir's (ALE + `SURFACE_ABSORPTION`,
`abstol = 1e-9`, `reltol = 1e-7`, `min_thickness = 0.5e-3`, the default
integrator). The script writes two CSVs — a per-time summary
(`ptw_results_50kWm2.csv`) and a spatial dump (`ptw_spatial_50kWm2.csv`) sampling
`solution.T` and `solution.ξ` at selected times and cell positions — which is the
pattern to copy when exporting profiles for external plotting.

---

## 14.5 Gas transport example (`08_gas_transport.jl`)

The example solves a charring slab with `LumpedTransfer` and `DustyGas`, using
independent wall-transfer data for the two constitutive closures. Both use the
same matrix permeability, skeletal pore volume, atmospheric air initialization,
complete reservoir and optional species film. It prints pressure, product
escape and species ledger defects using saved diagnostics.

```sh
PYROLYSIS_SMOKE=1 julia --project=examples examples/08_gas_transport.jl
```

The full example uses a moving mesh. The reported Darcy velocity is distinct
from total species flow, which includes wall/molecular transport.

---

## 14.6 Moisture heat-of-sorption study (`moisture_hos_comparison.jl`)

### What it models

A controlled study of how the *form* of the moisture heat of sorption affects a
Douglas fir cone simulation (50 kW/m², the solve of §14.3 on 240 cells). Two HoS
models are compared at four initial moisture contents (3, 6, 9, 12% dry
basis):

1. **Exponential (state-dependent):** `HoS(MC) = (a + b·exp(−c·MC))·1e3` J/kg —
   the energy cost rises as moisture depletes (tightly-bound water).
2. **Constant average:** a single value equal to the analytical integral average
   of the exponential model over `[0, MC₀]` — same total energy, uniform in time.

This isolates a single modeling choice (8 runs = 4 MC × 2 models, solved one
after another) and is the template for any "does this property form matter?"
investigation. The script writes CSVs and figures to
`examples/moisture_hos_comparison_results/`.

### The one thing that changes

Everything is identical between the two models except the `h` passed to the
moisture-evaporation reaction. The state-dependent case uses a
`StateDependentProperty`; the constant case passes a plain `Float64`, which is
auto-converted to a `ConstantProperty`:

```julia
if hos_model == :exponential
    h_evap = StateDependentProperty((T, ξ) -> begin
        MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS
        return (HOS_A + HOS_B * exp(-HOS_C * MC)) * 1e3   # J/kg
    end)
elseif hos_model == :constant
    avg_kJ = compute_average_hos_kJ(HOS_A, HOS_B, HOS_C, mc_pct)   # mc_pct = 100 mc_db
    h_evap = avg_kJ * 1e3        # Float64 -> ConstantProperty automatically
end

Reaction(:moisture_evaporation, 7 => 13; A = 2.997e5, E = 5.0e4, h = h_evap, n = 1.0)
```

This is the general rule for reaction heats and property fields: pass a number
for a constant, or any of the property-function types (`CallableProperty`,
`LinearProperty`, `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`,
`StateDependentProperty`) for richer dependence (Technical Reference §4).

### What to look at

The script extracts per-model scalar metrics (peak MLR, time to ignition, time to
moisture depletion, back temperature at 100 s) and percent differences, plus
difference curves `exp − const` for MLR and temperatures. The physically
interesting result is that matching the *integrated* energy is not enough — the
temporal distribution of the sorption sink shifts the drying-front passage and
hence the early MLR. The back temperature at a fixed time (`T_back_100s`) is a
sensitive indicator of how fast the drying front advances. Technical Reference §6
(state-dependent heats) and §16 (energy bookkeeping) cover the underlying physics.

---

## 14.7 Inverse analysis: parameter estimation (`inverse_analysis/inverse_white_pine.jl`)

### What it models

Recovery of unknown material parameters (component thermal conductivities of a
white-pine model) by fitting simulated back-surface temperature traces to
experimental CAPA cone data at multiple heat-flux levels. This is the practical
calibration workflow: wrap the forward solve in a loss function, minimize it over
the parameters, and quantify the uncertainty of the fit. The white-pine model
lives in three included files: `wp_properties.jl` (properties, kinetics and
constants), `wp_forward.jl` (the forward model) and `wp_objective.jl` (data
loading and the loss).

### The pattern

The driver is built from three pieces (the heavy lifting is in
`wp_objective.jl`, included at the top):

1. **`build_and_solve(θ, scenario; stage)`** (in `wp_forward.jl`) — builds a
   problem from a parameter vector `θ` and a scenario (flux level, BCs,
   thickness, cross section), then `solve`s it. This is the only place the
   package API appears; the rest is optimization plumbing. The mesh is
   `create_uniform_mesh(L, n_cells; cross_section = A_0)` with the sample's
   pre-test area, and the solve keywords are

   ```julia
   (abstol = 1e-7, reltol = 1e-7,
    dtmin = 1e-8, dtmax = 2.0, maxiters = 200_000,   # fail fast instead of grinding
    progress = verbose,
    use_ale = true, radiation_model = SURFACE_ABSORPTION, min_thickness = 0.5e-3)
   ```

   It returns a `PyrolysisSolution` (or `nothing` when the solve throws), from
   which the fit reads `sol.extras.back_temperature` (and `sol.surface_T` for
   the FSRI front-surface panels). A residual is formed only from a solve with
   `retcode == :Success`; a failed scenario adds a fixed penalty to the loss.

2. **A loss function** that, for each scenario `k`, interpolates the simulated
   back temperature onto the experimental time grid and accumulates a
   `σ`-weighted mean-square residual with the scenario weight `w[k]` (`1/N` for
   `N` scenarios by default), plus an optional conductivity prior:

   ```julia
   T_sim = interpolate_at(t_exp[mask], sol.t, sol.extras.back_temperature)
   r     = T_sim .- T_back_exp[mask]
   data_contrib += w[k] * mean((r ./ σ) .^ 2)
   ```

   The scenarios of one loss evaluation are solved concurrently, one
   `Threads.@spawn` task each, when Julia runs with more than one thread.

3. **A box-constrained optimizer.** The example uses `Optim.jl`'s
   `Fminbox(NelderMead())` (derivative-free, robust to the small loss-surface
   kinks that come from solver tolerances and interpolation), with multi-start
   initialization seeded by a Sobol sequence and an optional Stage-1 grid scan:

   ```julia
   using Optim
   res = optimize(loss_fn, lower, upper, θ0, Fminbox(NelderMead()),
                  Optim.Options(iterations = max_iter))
   θ_hat = Optim.minimizer(res)
   ```

### Uncertainty quantification

After the fit, parameter covariance is estimated either from the finite-
difference Hessian of the loss (default) or by residual bootstrap
(`--bootstrap=N`). The Hessian path uses `FiniteDiff.finite_difference_hessian`
and the weighted-NLS covariance `Σ_θ = (2·L(θ̂)/(N−p))·H⁻¹`, with `N` the number
of data points and `p` the number of parameters. Standard errors are
`sqrt.(diag(Σ))`.

### How to run it

```bash
julia -t auto --project=examples examples/inverse_analysis/inverse_white_pine.jl
julia -t auto --project=examples examples/inverse_analysis/inverse_white_pine.jl --stage=3
julia -t auto --project=examples examples/inverse_analysis/inverse_white_pine.jl --prior-weight=0
julia -t auto --project=examples examples/inverse_analysis/inverse_white_pine.jl --bootstrap=200
```

Each invocation writes a timestamped results directory
(`examples/inverse_analysis/results/run_<timestamp>_s<stage>_pw<weight>/`) with
the per-evaluation trace, the best parameters (`best.json` with standard errors
and the correlation matrix), and fit-vs-experiment plots. Its smoke run
(`PYROLYSIS_SMOKE=1`) makes one objective evaluation at the initial `θ` over
20 s and exits. Beside it, `validate_thickness.jl` runs the calibrated model
with its lateral shrinkage law (a variable cross section) at 30, 50 and
70 kW/m² and compares the final thickness, area and mass ratios with the CAPA
measurements, reading the area from `sol.extras.cross_section_area_history`;
`plot_diagnostics.jl` and `plot_stage2_theta.jl` forward-solve the scenarios at a
given conductivity vector and plot them against the experiments.

### Gradient-based alternative

This example uses a derivative-free optimizer on top of bare `solve` calls. A
gradient-based optimizer can use a central-difference gradient of the same loss,
at `2p` solves per gradient for `p` parameters (Chapter 13, §13.7; §13.4 for the
step):

```julia
using FiniteDiff, Optim

g!(G, θ) = FiniteDiff.finite_difference_gradient!(G, loss_fn, θ, Val(:central);
                                                  relstep = 1e-4)
res = optimize(loss_fn, g!, lower, upper, θ0, Fminbox(LBFGS()))
```

> Caveat: non-smooth lateral-shrinkage laws, deactivation events sampled at the
> data's times, and the solve's own error at small steps make the loss rough
> on the scale of a finite-difference step, which can confuse gradient-based
> methods; the derivative-free NelderMead route is robust to them. See §13.4,
> §13.9 and Technical Reference §17.

---

## 14.8 Sensitivity studies

Three scripts under `examples/sensitivity_analysis/` build a full
uncertainty-quantification pipeline around the Douglas fir cone case.

### 14.8.1 Forward uncertainty propagation (`sensitivity_douglas_fir.jl`)

Propagates uncertainty in 16 thermophysical inputs (instrument biases, heat-
capacity regression coefficients, heat-of-sorption parameters, and five heats of
reaction) through the simulation using Sobol quasi-random sampling. The inputs are
grouped into seven correlated blocks, each sampled jointly from a multivariate
Normal or scaled Student-t distribution with the appropriate covariance.

The structure mirrors §14.3 and §14.7: a single `build_and_solve(θ::Vector)`
(here `θ` has length 16) that maps the sampled parameters into property closures
and reaction heats, then solves with the standard ALE + surface-absorption stack
(`abstol = 1e-9`, `reltol = 1e-7`, `dtmin = 1e-14`, `min_thickness = 0.5e-3`,
`progress = false`, the default integrator). The 8192 samples are evaluated in
parallel, over the sample ids that have no result on disk yet, each solve in a
`try`/`catch` so that a failed sample does not stop the study:

```julia
using Distributed
@everywhere begin
    using Pyrolysis
    using Pyrolysis: get_boundary_exchange
end

samples = generate_samples(CORR_GROUPS, N_SAMPLES)   # 16 × 8192 parameter matrix
remaining_ids = [i for i = 1:N_SAMPLES if i ∉ completed_ids]
results = pmap(remaining_ids) do i
    try
        result = extract_qois(build_and_solve(samples[:, i]))
        # ... write the run's curves to results/ at once ...
        (i, result)
    catch err
        (i, nothing)
    end
end
```

```bash
julia -p 8 --project=examples examples/sensitivity_analysis/sensitivity_douglas_fir.jl   # 8 workers
```

An existing `results/` directory is reloaded and completed, so an interrupted
study resumes. The smoke run makes one forward solve at the nominal parameters
over 10 s and writes no files.

Each run is reduced to scalar QoIs (peak MLR, time to ignition, surface/back
temperature rise times, combustible-only burning metrics, etc.) and a few
curves; the script then computes Standardized Rank Regression Coefficients
(SRRC) for the sensitivity ranking and prediction-envelope percentiles. Note the
public `get_boundary_exchange` accessor: the script sums the assembled
species fluxes over combustible gases, with the current face area and
geometry already accounted for.

### 14.8.2 Global Sobol indices (`sensitivity_sobol.jl`, `compare_degrees.jl`)

`sensitivity_sobol.jl` consumes the `parameters.csv` and `scalar_qois.csv` written
by §14.8.1, fits an adaptive-degree (2–4) Polynomial Chaos Expansion surrogate to
each QoI, and reads off first-order and total group Sobol indices analytically
from the PCE coefficients, with BCa-bootstrap confidence intervals. It does not
re-run the solver — it is pure post-processing of the sampled outputs.
`compare_degrees.jl` is a companion that settles the PCE degree choice (2 vs 3 vs
4) by a convergence rule on the inter-degree change in the Sobol indices.

These two are statistical post-processing of the §14.8.1 outputs; the only
Pyrolysis-specific dependency is the QoI/parameter CSVs produced by the sampling
run. They need `PolyChaos`, which `examples/Project.toml` does not list: add it
to your environment before running them. The variance-decomposition theory is
Technical Reference §17.

### 14.8.3 The `parameter_sweep` helper

For sweeps that do not need a fully bespoke driver, the package provides
`parameter_sweep(problem, grid; modify, extract)`. The `modify` closure builds a
fresh problem from each grid point (it must not mutate the base), and `extract`
reduces each solution to a small summary so workers ship little data back:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]

@everywhere function modify(prob, p)
    m = prob.material
    new_rxn = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 200e3)
    new_mat = Material(; name = m.name, components = m.components,
                       reactions = (new_rxn,), mixing = m.mixing)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
                            T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
                            tspan = prob.tspan, output_times = prob.output_times)
end

@everywhere extract(sol) = (peak_T = maximum(sol.surface_T),
                            final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```

Here `problem` is a base problem whose material has a component 1 that
decomposes into a component 2. The sweep runs in parallel via
`Distributed.pmap` when `nworkers() > 1`, otherwise sequentially; the workers
need `modify` and `extract` (`@everywhere`). Every keyword besides `modify`,
`extract` and `parallel` goes to `solve`, the same at every point. See
Chapter 12 for the full contract.

---

## 14.9 Literature property comparison (`literature_comparison/literature_comparison.jl`)

### What it models

Nine runs of the Douglas fir cone case at 50 kW/m² (one baseline plus eight
variants), each substituting one category of property — virgin `c_p`, char
`c_p`, heats of reaction, or heat of sorption — with a value sourced from the
literature, holding everything else fixed (method and sources in
`literature_comparison_plan.md` beside the script). It quantifies how much each
property choice moves the predicted MLR and temperatures, and is a model for any
"which published property set should I trust?" study.

### The pattern: serializable run specs

Because the variants run on parallel workers, the script avoids serializing
closures across processes. It encodes each run as a plain struct of *symbols*
that index into worker-local dispatch dictionaries of property functions:

```julia
struct LitCompRunSpec
    label::String
    virgin_cp_key::Symbol      # :baseline / :harada / :dupont
    char_cp_key::Symbol        # :baseline / :hostikka / :edrich
    hrxn_scale::Float64        # uniform multiplier on all nominal heats of reaction
    hrxn_scale_secondary::Float64  # two-zone secondary scale; NaN -> use hrxn_scale
    hos_key::Symbol            # :baseline / :free_water
end

const VIRGIN_CP_DISPATCH = Dict{Symbol,Function}(:baseline => baseline_virgin_cp,
                                                :harada   => harada_virgin_cp,
                                                :dupont   => dupont_virgin_cp)
# build_and_solve(spec) looks up the function by key, wraps it in CallableProperty,
# builds the material on create_uniform_mesh(12.0e-3, 240), and solves with the
# standard ALE + SURFACE_ABSORPTION stack.
```

This symbol-dispatch idiom is the recommended way to vary functional properties
across a distributed sweep, since `Function` values do not serialize cleanly to
workers.

### How to run it

```bash
julia -p 9 --project=examples examples/literature_comparison/literature_comparison.jl
julia --project=examples examples/literature_comparison/literature_comparison.jl   # serial
```

### What to look at

Per-run CSVs (`{label}_results.csv`: time, MLR, surface/back temperature, total
mass, thickness) and a summary QoI table (`literature_comparison_qois.csv`), in
`examples/literature_comparison/results/`. The comparison reveals which property categories dominate the spread in predicted
burning behavior — typically the heats of reaction and the high-temperature
char/virgin `c_p` — which feeds back into prioritizing measurement effort.
Technical Reference §18 frames this as validation against alternative inputs.

---

## 14.10 Common patterns across the examples

Reading the examples together, a few idioms recur and are worth internalizing:

- **Build, do not mutate.** For sweeps and parallel work, wrap material/problem
  construction in a function that takes the varied parameters and returns a fresh
  `PyrolysisProblem`. Each `solve` sets up its own definition, workspace and
  Jacobian cache from the problem, so nothing is shared between concurrent
  solves, whether on `Distributed` workers or on threads.
- **Temperature-dependent properties are closures.** `CallableProperty(T -> ...)`
  for `c_p(T)`, `k(T)`; `StateDependentProperty((T, ξ) -> ...)` for properties or
  reaction heats that depend on local composition (e.g. heat of sorption). A bare
  number becomes a `ConstantProperty`.
- **Combine BCs with `+`.** `RadiativeFluxBC(...) + RadiativeBC(...) + ConvectiveBC(...)`
  forms a `CombinedThermalBC`. Assign one complete mass closure to each tag.
- **Declare the model and geometry choices.** `Material` selects
  `LumpedTransfer` or `DustyGas`, each with independent Darcy permeability. The moving mesh, the radiation model and deactivation
  are explicit `solve` keywords (`use_ale`, `radiation_model`, `deactivation`).
  A variable cross section turns on when the `Material` is given a
  `lateral_shrinkage_law` (a function `χ̄ -> A/A_0` of the column's conversion,
  e.g. `chi -> 1 - 0.2chi`), which needs `use_ale = true`. See Chapter 7 and
  Technical Reference §10.
- **Pore air and complete reservoirs.** In `07_douglas_fir_cone.jl`,
  `08_pressure_treated_wood_cone.jl` and `08_gas_transport.jl` the solids
  declare skeletal densities (`ρ_intrinsic`) and an inert gas component starts
  at a 1 atm partial pressure (`background_gas_concentration`), so the pore
  space is not a vacuum against the atmospheric pressure boundary (Technical
  Reference §9).
- **Read the right MLR.** `solution.extras.surface_mass_flux` is the physical
  surface efflux; `solution.extras.gas_generation_rate` is the gas the
  reactions produce; `solution.mass_loss_rate` is the outflow per current area
  and `solution.extras.MLR_gauge` the outflow per initial area, as a cone
  calorimeter reports it (they differ under a variable cross section). The
  efflux and the production coincide only at quasi-steady state.
- **Use the boundary-flux extras for energy accounting.** `surface_heat_flux_*` and
  `surface_emissivity`/`surface_conductivity` in `extras` give the surface energy
  balance without re-deriving it; under `BEER_LAMBERT` add
  `extras.volumetric_absorption`, the in-depth deposit the face fluxes exclude.

### Caveats to keep in mind

- **Radiation models.** The radiation models are `NO_RADIATION`,
  `SURFACE_ABSORPTION` and `BEER_LAMBERT`. The shipped examples use
  `NO_RADIATION` (flux via the BC, applied in the surface balance exactly as
  under `SURFACE_ABSORPTION`) or `SURFACE_ABSORPTION`; `BEER_LAMBERT` is one
  keyword away in `06_cone_calorimeter.jl`.
- **Boundary energy.** Incoming gas carries reservoir enthalpy and outgoing
  gas carries surface enthalpy. The coupled LTE film converts heat to the
  reservoir control surface; adiabatic and prescribed-heat boundaries use the
  incoming enthalpy directly. `get_boundary_exchange` reports the assembled
  species flux, carried enthalpy, heat and total energy (§5.4).
- **Experimental validation.** The cone examples need validation of their
  material parameters against experimental data. Numerical convergence and
  conservation alone do not establish that calibration (§14.1).

For deeper troubleshooting (rejected steps, ledger imbalance, crossed mesh
nodes, radiation surprises), see the Troubleshooting and Performance chapter
(Chapter 15); for the full keyword reference, see the Running the Solver and
API Reference chapters (Chapters 8 and 16).
