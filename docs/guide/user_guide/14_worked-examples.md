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
| `06_cone_calorimeter.jl` | PMMA cone, single-step kinetics, incident flux at the surface (`SURFACE_ABSORPTION`) | ALE, deactivation | 14.1 |
| `07_douglas_fir_cone.jl` | Wood, multi-step kinetics, lumped gas transport, the phase change of the moisture | ALE, SURFACE_ABSORPTION, LumpedTransfer | 14.2 |
| `convergence_douglas_fir.jl` | Grid-convergence study | as 07 at 50 kW/m², swept over `n_cells` | 14.3 |
| `08_pressure_treated_wood_cone.jl` | Treated wood transcribed from ThermaKin, its one-way drying step, lumped gas transport | ALE, SURFACE_ABSORPTION, LumpedTransfer | 14.4 |
| `08_gas_transport.jl` | Lumped and dusty-gas transport, finite film | Moving mesh, both transport models | 14.5 |
| `macfp/macfp_gasification.jl` | MaCFP-4 white-pine gasification targets (10/30/60 kW/m², 6/18 mm), the calibrated model of §14.6 at a constant cross section | ALE, SURFACE_ABSORPTION, threads | — |
| `inverse_analysis/inverse_white_pine.jl` (with `validate_thickness.jl`, `plot_diagnostics.jl`, `plot_stage2_theta.jl`) | Parameter estimation from data | ALE, optimization | 14.6 |

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
insulated and impermeable. The solve uses ALE, surface absorption and
surface deactivation.

```bash
PYROLYSIS_SMOKE=1 julia --project=examples examples/06_cone_calorimeter.jl
```

Inspect `sol.gas.pressure`, `sol.gas.porosity` and `get_boundary_exchange` alongside
temperatures, mass loss and thickness. Mesh, nominal porosity and pressure
error budget must be varied jointly when assessing dense-polymer results.

---

## 14.2 Douglas fir cone calorimeter (`07_douglas_fir_cone.jl`)

The wood example includes multi-step pyrolysis, bound moisture and its vapour,
and inert pore air. The moisture and its vapour exchange through a two-way phase
change on the CAST spruce isotherm, held by the wood and its four intermediates and
not by the char (§4.14):

```julia
phase_changes = (
    PhaseChange(:moisture, :MOISTURE => :WATER_VAPOR;
                activity = CAST_SPRUCE,
                sorbents = (:DF, :DF_int1, :DF_int2, :DF_int3, :DF_int4)),
),
```

The water evaporates where the wood heats and condenses where the vapour driven
inward meets cooler wood, so the moisture rises ahead of the drying front and the
back face holds near the temperature at which the bound water's equilibrium pressure
meets the pore pressure. The pore vapour starts in equilibrium with the water and the
air at the rest of 101325 Pa (`equilibrium_vapour_pressure`,
`equilibrium_vapour_concentration`, `background_gas_concentration`), and the ambient
air carries the humidity in equilibrium with the initial moisture at 300 K (30.7 %),
so the conditioned sample is at rest before the heater acts. Condensed phases provide
skeletal densities, lumped matrix transfer and permeability. The complete atmospheric
reservoir has incoming air at the same temperature as the external thermal film.

```bash
PYROLYSIS_SMOKE=1 julia --project=examples examples/07_douglas_fir_cone.jl
```

The moving mesh, radiation and deactivation choices are explicit solve
keywords. Use the species boundary diagnostics to separate combustible-gas
escape from vapor transport and air ingress. Reconstructing flux from the
initial mesh or only the pressure gradient would omit parts of this balance.

The script reports the water evaporated and condensed over the run from
`extras.evaporation_rate` and `extras.condensation_rate` (§9.3), and its results CSV
carries the net water-vapour production (`mlr_gen_water`, evaporation less
condensation) and the net exchange at the drying front (`phase_change_rate_front`).
At the default exchange rate about 40 % of the water evaporated condenses again, and
the vapour runs at most 11 % above its equilibrium pressure; the results are flat in
the exchange rate from `1e2` to `1e4` 1/s (Technical Reference §6.7.5).

---

## 14.3 Grid-convergence study (`convergence_douglas_fir.jl`)

### What it models

The Douglas fir cone case of §14.2 at 50 kW/m² (14 components, five reactions and
the phase change of the moisture), solved at six spatial resolutions
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
    # ... build material (14 components, 5 reactions, 1 phase change) and BCs,
    # the reservoir's humidity in equilibrium with the initial moisture ...
    mesh = uniform_mesh(12.0e-3, n_cells)
    ξ_DF = 530.0 / (1 + 0.06 * 530.0 / 1000.0)   # wood and 6 % bound water fill the cell (§3.2.2)
    ξc = (ξ_DF, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06 * ξ_DF, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    p_v = equilibrium_vapour_pressure(material, ξc, 300.0)
    ξ_v = equilibrium_vapour_concentration(material, ξc, 300.0)
    ξ_air = background_gas_concentration(material, :AIR, ξc, 300.0; P = 101325.0 - p_v)
    ξ_initial_funcs = [z -> ξ_DF, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
                       z -> 0.06 * ξ_DF,                          # 6 % moisture, dry basis
                       z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0, z -> 0.0,
                       z -> ξ_v, z -> ξ_air]                      # vapour, air
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
(§12.6 and §14.7); this script calls `pmap` over its own
`build_and_solve(n_cells)`, which builds the whole problem at each resolution.

---

## 14.4 Pressure-treated wood (`08_pressure_treated_wood_cone.jl`)

### What it models

Cone pyrolysis of fire-retardant pressure-treated wood (PTW) with moisture, at
50 kW/m², 12 mm on 240 cells, transcribed from a ThermaKin deck (`PTW.cmp`) for
comparison with it. Structurally similar to Douglas fir but with a different
component ordering (water and water vapor come first), three solid intermediates
between the virgin wood and the char instead of four, and five reactions. Its water
leaves by ThermaKin's one-way drying reaction, not by a phase change, so the vapour it
drives inward does not condense (§4.14.5). It exercises the same lumped gas transport +
surface-absorption + ALE stack.

### Material highlights

Twelve components: `WATER` (liquid), `WATER_VAPOR` (gas), `PTW`/`PTW1`/`PTW2`/
`PTW3`/`Char` (solids, with skeletal densities `ρ_intrinsic`), `PTW_gas1`…`PTW_gas4`
(gases), and the background `AIR` (gas, initialized to 1 atm in the pores by
`background_gas_concentration`). The liquid water is bound water absorbed into
the cell walls, not a separate phase:

```julia
LiquidComponent(:WATER;
    ρ = 1000.0,        # liquid water [kg/m³]
    c = water_cp,
    k = 0.6,           # liquid water [W/(m·K)]
    ε = 0.81,
    α = 1000.0,        # not used with SURFACE_ABSORPTION
    κ = 0.0,           # liquid blocks Darcy flow
)
```

The water holds its own volume, so drying shrinks the wood, and it conducts
and emits by its volume share. It declares neither `λ` nor `κ`, so it is no
part of the gas-transport matrix; its mass, heat capacity, evaporation energy
and skeletal occupancy remain. A liquid that should join the matrix, such as
a melt, declares its transport coefficients explicitly.

### Reaction heats and rates

Pyrolysis.jl uses `h > 0` for endothermic reactions and `h < 0` for
exothermic reactions. Water evaporation (ThermaKin's step, with its heat sign
flipped) and the first decomposition step absorb heat; the third decomposition step
releases it:

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

## 14.6 Inverse analysis: parameter estimation (`inverse_analysis/inverse_white_pine.jl`)

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
   `uniform_mesh(L, n_cells; cross_section = A_0)` with the sample's
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

## 14.7 The `parameter_sweep` helper

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
                       reactions = (new_rxn,), phase_changes = m.phase_changes,
                       mixing = m.mixing,
                       gas_viscosity = m.gas_viscosity, gas_transport = m.gas_transport)
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

## 14.8 Common patterns across the examples

Reading the examples together, a few idioms recur and are worth internalizing:

- **Build, do not mutate.** For sweeps and parallel work, wrap material/problem
  construction in a function that takes the varied parameters and returns a fresh
  `PyrolysisProblem`. Each `solve` sets up its own definition, workspace and
  Jacobian cache from the problem, so nothing is shared between concurrent
  solves, whether on `Distributed` workers or on threads.
- **Temperature-dependent properties are closures.** `CallableProperty(T -> ...)`
  for `c_p(T)`, `k(T)`; `StateDependentProperty((T, ξ) -> ...)` for properties or
  reaction heats that depend on local composition (e.g. a conversion-dependent heat). A bare
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
  balance without re-deriving it; under `RADIATIVE_TRANSFER` add
  `extras.radiation_absorbed − radiation_emitted − radiation_back`, the radiation the
  cells take in, which the face fluxes exclude.

### Caveats to keep in mind

- **Radiation models.** The radiation models are `SURFACE_ABSORPTION` and
  `RADIATIVE_TRANSFER`. The shipped examples use `SURFACE_ABSORPTION` (flux via the BC,
  applied in the surface balance). `06_cone_calorimeter.jl` runs under
  `RADIATIVE_TRANSFER` once its `RadiativeFluxBC` drops `absorptivity` and its
  `RadiativeBC` members drop `ε` (§2.6.2).
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
