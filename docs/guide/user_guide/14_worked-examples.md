# 14. Worked Examples and Tutorials

This chapter walks through the example scripts shipped in the `examples/`
directory of Pyrolysis.jl. Each example is a complete, runnable workflow that
exercises a different combination of physics (radiation models, Darcy–Fick
transport, moisture, ALE shrinkage), solver options, and post-processing. For
each one we state what is being modeled, which simulation modes and boundary
conditions are configured, how to run it, what outputs to look at, and which
chapters of the Technical Reference derive the underlying physics. The code
fragments here are drawn from the shipped scripts; where a snippet is
abbreviated, the original file is cited so you can read it in full.

All example code uses the package's coordinate convention: `z = 0` is the
bottom/substrate (boundary id 2, tag `:bottom`) and `z = L` is the top/exposed
surface (boundary id 1, tag `:top`); heat enters from the top and the material
recedes downward. SI units are used throughout (see Technical Reference §2).

The examples live in `examples/` with their own `Project.toml` environment that
adds `Plots` (and a few helper packages) on top of Pyrolysis. To run them:

```julia
cd("examples")
using Pkg
Pkg.activate(".")
Pkg.instantiate()         # first time only

include("06_cone_calorimeter.jl")
```

Or from a shell:

```bash
cd examples
julia --project=. -e 'include("06_cone_calorimeter.jl")'
```

A roadmap of what each example demonstrates:

| Example file | Physics highlighted | Modes | §14 |
|---|---|---|---|
| `06_cone_calorimeter.jl` | PMMA cone, single-step kinetics, incident flux at the surface (`NO_RADIATION`) | ALE, depletion | 14.1 |
| `06_cone_calorimeter_surface_absorption.jl` | PMMA cone, surface-only radiation | ALE, SURFACE_ABSORPTION | 14.1 |
| `07_douglas_fir_cone.jl` | Wood, multi-step kinetics, Darcy–Fick, moisture | ALE, SURFACE_ABSORPTION, DarcyFick | 14.2 |
| `convergence_douglas_fir.jl` | Grid-convergence study | as 07, swept over `n_cells` | 14.3 |
| `08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl` | Treated wood + moisture, Darcy–Fick | ALE, SURFACE_ABSORPTION, DarcyFick | 14.4 |
| `08_darcy_fick_demo.jl` | Pressure-driven gas transport, permeability | ALE, DarcyFick | 14.5 |
| `moisture_hos_comparison.jl` | State-dependent heat of sorption | ALE, SURFACE_ABSORPTION, DarcyFick | 14.6 |
| `macfp/macfp_gasification.jl` | MaCFP-4 white-pine gasification targets (10/30/60 kW/m², 6/18 mm) | ALE, SURFACE_ABSORPTION | — |
| `inverse_analysis/inverse_white_pine.jl` | Parameter estimation from data | ALE, optimization | 14.7 |
| `sensitivity_analysis/sensitivity_douglas_fir.jl`, `sensitivity_sobol.jl`, `compare_degrees.jl`, `df_sensitivity_comparison.jl` | Uncertainty / Sobol / inter-code | parameter sweep, post-process | 14.8 |
| `literature_comparison/literature_comparison.jl` | Property-variant comparison | as 07, swept over property sets | 14.9 |

---

## 14.1 Cone calorimeter: PMMA (`06_cone_calorimeter.jl`)

### What it models

Standard cone-calorimeter pyrolysis of poly(methyl methacrylate) (PMMA), a
non-charring thermoplastic that depolymerizes almost entirely to gaseous MMA
monomer. The sample is a 5.8 mm slab heated from the top by a radiant cone
heater; it recedes as the material decomposes and gas escapes. This is the
canonical "your first cone run" workflow: one solid that becomes a small char
residue plus a gas, temperature-dependent `c_p` and `k`, combined thermal
boundary conditions, and ALE mesh motion so the shrinking surface is tracked
exactly. The physics is derived in Technical Reference §3 (governing
equations), §6 (kinetics), §7 (energy assembly), §8 (radiation), and §10–§11
(volume change and ALE).

### Material

The material has three components — virgin PMMA, char, and MMA gas — and one
endothermic decomposition reaction. Temperature-dependent properties are
supplied as `CallableProperty` closures (piecewise functions of `T`):

```julia
using Pyrolysis

# Piecewise temperature-dependent heat capacity and conductivity (break at 378 K)
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2 * T) : (1040.0 + 3.08 * T))
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * T))

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin;
            ρ = 1210.0,      # bulk density [kg/m³]
            c = pmma_cp,     # specific heat [J/(kg·K)]
            k = pmma_k,      # thermal conductivity [W/(m·K)]
            ε = 0.96,        # surface emissivity
            α = 1.47,        # mass absorption coefficient [m²/kg] (used by BEER_LAMBERT)
        ),
        SolidComponent(:char;
            ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96, α = 1.47,
        ),
        GaseousComponent(:MMA;
            M = 0.10012,     # molar mass [kg/mol]
            c = 1500.0,      # gas heat capacity [J/(kg·K)]
            k = 0.02,        # gas conductivity [W/(m·K)]
            λ = 1e-5,        # gas transfer (diffusion) coefficient [m²/s]
        ),
    ),
    reactions = (
        # virgin (1) -> char (2) + MMA gas (3); products carry mass yields
        Reaction(:decomposition, 1 => (2, 3);
            A = 1.5e14,            # pre-exponential factor [s⁻¹]
            E = 2.03e5,            # activation energy [J/mol]
            h = 820e3,             # heat of reaction [J/kg] (h > 0 = endothermic)
            n = 1.0,               # reaction order
            yields = (0.002, 0.998),  # 0.2% char, 99.8% gas; must sum to 1.0
        ),
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)
```

Annotations on the non-obvious arguments:

- **`α` (mass absorption coefficient, m²·kg⁻¹).** Only consumed when
  `radiation_model = BEER_LAMBERT`; with `SURFACE_ABSORPTION` it is ignored
  (all incident flux lands in the surface cell). See §8.
- **`h` sign.** Pyrolysis.jl uses the storage convention `h > 0 endothermic`
  (the volumetric source is `Q_rxn = −Σ h_r r_r`). This is the opposite sign
  to ThermaKin/Gpyro publications — see Technical Reference §6 and the
  nomenclature overload note H1. PMMA depolymerization is endothermic, hence
  `h = +820e3`.
- **`yields`.** For a single-reactant reaction with two products, the implicit
  reactant coefficient is `−1` and the product yields must sum to `1.0`
  (mass conservation, enforced at construction).
- **`mixing`.** `MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL))`
  selects the conductivity mixing rule (`WEIGHTED` with `β = 0.5`, a 50/50
  parallel/series blend) and the gas-transfer mixing rule (`PARALLEL`). The
  third slot, permeability `κ`, defaults to `MixingSpec(WEIGHTED, 0.5)`.
  Mixing rules are derived in Technical Reference §5.

### Mesh

A uniform mesh is created with `create_uniform_mesh(L, n_cells, Val(NC); ...)`.
The third positional argument is a `Val` carrying the component count `NC` (a
compile-time constant for type stability):

```julia
sample_thickness = 5.8e-3   # 5.8 mm
n_cells = 72                # ~0.08 mm cells

mesh = create_uniform_mesh(
    sample_thickness, n_cells, Val(3);   # NC = 3 (virgin, char, gas)
    T_initial = 300.0,                   # initial temperature [K]
    ξ_initial = (1210.0, 0.0, 0.0),      # initial concentrations [kg/m³]: full virgin
)
```

`ξ_initial` is a per-component tuple of initial mass concentrations `ξ_j`
(kg·m⁻³), not mass fractions. The cross-sectional area defaults to `A_0 = 1.0`
m² (`cross_section_area` keyword); in 1D it is a gauge normalization for
per-area outputs (see the Meshes and Geometry chapter).

### Boundary conditions

Thermal BCs are composed with the `+` operator, which builds a
`CombinedThermalBC`. The top surface receives a radiative cone flux, re-radiates
to the surroundings, and loses heat by convection; the bottom re-radiates and
convects to ambient and is impermeable:

```julia
q_external = 0.9813 * 25e3   # 25 kW/m² cone flux (with a radial-uniformity factor)
h_conv = 7.2                 # convective coefficient [W/(m²·K)]

top_thermal_bc =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +  # incident cone flux
    RadiativeBC(ε = 0.96, T_inf = 330) +                       # re-radiation loss
    ConvectiveBC(h = h_conv, T_inf = 330)                      # convective loss

top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)        # external mass transfer

bottom_thermal_bc = RadiativeBC(ε = 0.95, T_inf = 307) + ConvectiveBC(h = 4, T_inf = 307)
bottom_mass_bc = ImpermeableBC()

bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal_bc, bottom = bottom_thermal_bc),
    mass    = (top = top_mass_bc,    bottom = bottom_mass_bc),
)
```

Key points:

- **`RadiativeFluxBC`** supplies the incident cone flux. `flux` may be a number
  (constant) or a function of time `t -> q(t)`; `absorptivity` fixes the
  absorbed fraction (default `nothing` → use the composition-weighted effective
  emissivity via Kirchhoff). With `SURFACE_ABSORPTION` the absorbed flux is
  applied at the surface through the boundary energy balance; with
  `BEER_LAMBERT` it is attenuated in depth.
- **`RadiativeBC(; ε = nothing, T_inf, F = 1.0)`** is a *loss* term
  `q = Fεσ(T_∞⁴ − T_s⁴)`, not a heat input. With `ε = nothing` it uses the
  material's effective emissivity (which tracks the virgin→char transition).
- **`ConvectiveMassBC(; h_m, ξ_inf = 0.0)`** lets gas escape with an external
  film coefficient `h_m`; `DiffusiveMassBC(h_m)` is the diffusion-resistance
  variant, and `DiffusiveMassBC()` (no argument) is pure internal diffusion
  (`h_m = ∞`). Mass BCs are derived in Technical Reference §9 and §12.

### Problem and solve

```julia
problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bc_set,
    T_initial = z -> 300.0,                  # initial-condition closures (function of z)
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, 1600.0),
    dt_initial = 1e-3,
    output_times = collect(range(0.0, 1600.0, length = 601)),  # save every ~2.7 s
)

solution = solve(problem;
    abstol = 1e-8, reltol = 1e-6,
    dtmin = 1e-14,                # allow tiny steps through the rapid reaction front
    progress = true,
    use_ale = true,              # moving mesh tracks the receding surface
    min_thickness = 0.01e-3,     # terminate when sample < 0.01 mm
    handle_depletion = true,     # park thin surface cells
    depletion_threshold = 0.05,  # park a cell at 5% of its initial thickness
    min_cells = 1,               # keep at least 1 cell
)
```

Notes on the solve options (full reference: Running the Solver chapter; defaults
from `src/Solver/solve.jl`):

- **`use_ale = true`** adds node positions `z` to the state vector and evolves
  them so the mesh follows the material (Technical Reference §11). Without it the
  mesh is Eulerian (fixed) and the surface does not recede.
- **`radiation_model`** is omitted here, so it defaults to `NO_RADIATION`. The
  incident flux still enters through `RadiativeFluxBC`, but to distribute it in
  depth you would pass `radiation_model = BEER_LAMBERT` (PMMA is
  semi-transparent; the `α = 1.47` values then matter). The surface-absorption
  variant (§14.1.1) passes `SURFACE_ABSORPTION`.
- **`handle_depletion = true`** turns on cell parking (Chapter 10), which
  carries the run through burnout (Technical Reference §14.3).

### What to look at

Outputs live on the returned `PyrolysisSolution`:

```julia
times    = solution.t                                    # saved times [s]
T_surf   = solution.surface_T                            # surface temperature [K]
T_back   = solution.extras.back_temperature              # back-face temperature [K]
total_m  = solution.extras.total_mass                    # mass per initial area [kg/m²]
gen_rate = solution.extras.gas_generation_rate           # volumetric gas production [kg/(m²·s)]
thick    = solution.extras.thickness_history             # current thickness [m]

println("Return code: ", solution.retcode)
println("Peak surface T: ", round(maximum(T_surf), digits = 1), " K")
println("Mass loss: ", round((1 - total_m[end]/total_m[1]) * 100, digits = 1), " %")
```

`surface_T` and `back_temperature` come from the boundary-temperature Newton
solve (not raw cell-center temperatures), which closes the surface energy
balance — important for `T⁴` radiative BCs (Technical Reference §12). The
script prints a side-by-side comparison with ThermaKin values at t ≈ 600 s
(hard-coded in the script) and overlays ThermaKin traces only when
`docs/thermakin_output/` is present; the repository does not ship those traces.

If you have `Plots` loaded, the package's plot recipe gives quick views:

```julia
using Plots
plot(solution; kind = :surface_T)   # surface temperature vs time
plot(solution; kind = :mlr)         # mass-loss rate
plot(solution; kind = :heatmap)     # T(z, t) field
```

> Caveat: this example's parameter set is pending re-validation against measured
> data (Technical Reference §18.8.1), so the ThermaKin comparison is indicative,
> not a regression target; treat the printed peak-surface-T and mass-loss numbers as a
> reference for magnitude and time scale rather than exact values to match.

### 14.1.1 Surface-absorption variant (`06_cone_calorimeter_surface_absorption.jl`)

This companion script is identical in spirit but selects the surface-only
radiation model. The absorption coefficients are set to `α = 0.0` to make
explicit that they are unused, and the solve passes
`radiation_model = SURFACE_ABSORPTION`:

```julia
solution = solve(problem;
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,   # all incident flux in the surface cell
    min_thickness = 0.01e-3,
    handle_depletion = true,
)
```

`SURFACE_ABSORPTION` is `O(1)` per step with a tridiagonal Jacobian, versus the
upper-triangular non-local Jacobian of `BEER_LAMBERT`; it is the right choice
for optically thick (opaque) materials where `τ = αΔz ≫ 1`. Compare the two
scripts to see the surface-temperature and MLR difference between in-depth and
surface deposition (Technical Reference §8).

---

## 14.2 Cone calorimeter: Douglas fir (`07_douglas_fir_cone.jl`)

### What it models

Pyrolysis of Douglas fir wood — an opaque charring solid with moisture — under a
75 kW/m² cone. This is the flagship "everything on" example: multi-step
decomposition kinetics (five solid→solid+gas steps that march virgin wood
through four intermediates to char, plus a moisture-evaporation step), a
state-dependent heat of sorption, Darcy–Fick gas transport with
kinetic-theory-sized gas diffusion coefficients, surface-absorption radiation,
and ALE mesh motion with the
analytical (structured) Jacobian for speed. The physics spans Technical
Reference §5–§12.

### Material highlights

The material declares 14 components — six solids (`DF`, `DF_int1`…`DF_int4`,
`Char`), one liquid (`MOISTURE`), and seven gases (`DF_gas1`…`DF_gas5`,
`WATER_VAPOR`, and the inert background `AIR` initialized at 1 atm) — and six
reactions. A few patterns worth lifting:

**Kinetic-theory sizing of the gas `λ` values.** Each gas component's transfer
coefficient `λ` is *assigned* a value scaled by `√(M_ref/M_j)` from kinetic
theory (lighter molecules get a larger `λ`), where `M_ref = 0.029` kg/mol
(air). These per-component values feed the single mixture coefficient that all
species share at a face — the flux kernels apply no per-species scaling of
their own (Technical Reference §9.2):

```julia
λ_gas = 3.0e-5  # reference free-molecular diffusion at ~500 K [m²/s]

GaseousComponent(:DF_gas3; M = 0.110, c = 1300.0, k = 0.018,
                 λ = λ_gas * sqrt(0.029 / 0.110))   # heavy tar: slowest
GaseousComponent(:WATER_VAPOR; M = 0.018, c = water_vapor_cp, k = 0.025,
                 λ = λ_gas * sqrt(0.029 / 0.018))   # H₂O: fastest
```

**Composition-dependent permeability.** Virgin wood is nearly impermeable
(`κ ≈ 1e-15` m²) and char is open-pored (`κ ≈ 1e-12` m²); the per-component `κ`
values make the effective permeability rise naturally as the matrix chars,
driving Darcy flow (Technical Reference §9). Any non-zero `κ` switches the
problem into Darcy–Fick transport mode automatically.

**State-dependent heat of sorption.** The moisture-evaporation reaction uses a
`StateDependentProperty` whose value depends on the local moisture concentration
(more energy is needed to release the remaining tightly-bound water):

```julia
const MOISTURE_IDX = 7
const ρ_DRY_BASIS = 530.0   # dry wood density for the % moisture-content basis

heat_of_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS    # moisture content [% dry basis]
    return (2450.42 + 875.64 * exp(-0.1122 * MC)) * 1e3   # J/kg
end)

# ...used as the heat of one reaction:
Reaction(:moisture_evaporation, 7 => 13;
    A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0)
```

`StateDependentProperty` closures receive `(T, ξ)` where `ξ` is the full
component tuple. They are AD-clean and may be used for any property field or for
a reaction heat. See Technical Reference §4 (property functions) and §6.

**Multi-step kinetics.** The five wood-decomposition reactions are sequential —
each consumes the previous intermediate. Pyrolysis.jl enforces a single reactant
per reaction (Technical Reference §6), so multi-step schemes are built as a chain
of single-reactant reactions:

```julia
Reaction(:pyrolysis_step1, 1 => (2, 8);  A = 74816.95,      E = 71.59e3,
         h = -19200.66, n = 1.0, yields = (0.95, 0.05)),   # exothermic (h < 0)
Reaction(:pyrolysis_step3, 3 => (4, 10); A = 7.744618e12,   E = 181.0e3,
         h = 264424.79, n = 1.0, yields = (0.37, 0.63)),   # endothermic (h > 0)
# ... steps 2, 4, 5 similarly
```

The full 14-component, six-reaction definition is the `douglas_fir = Material(...)`
block of `examples/07_douglas_fir_cone.jl`.

### Boundary conditions

The top combines incident flux (using the material's effective emissivity by
omitting `ε`), re-radiation, and convection; the back is adiabatic. The top
mass BC combines diffusive escape with a Dirichlet pressure (atmospheric exhaust)
so gas can both diffuse and be pushed out by Darcy flow:

```julia
top_thermal_bc =
    RadiativeFluxBC(flux = 75e3) +            # uses material ε (ε = nothing default)
    RadiativeBC(T_inf = 300) +
    ConvectiveBC(h = 10.0, T_inf = 300)

top_mass_bc = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal_bc, bottom = AdiabaticBC()),
    mass    = (top = top_mass_bc,    bottom = ImpermeableBC()),
)
```

`DirichletPressureBC(P)` sets the surface pressure for the Darcy velocity
computation (Technical Reference §9, §12). Pressure BCs are only meaningful when
the material has permeability (Darcy–Fick mode).

### Solve

```julia
using OrdinaryDiffEq, LinearSolve
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve: SparspakFactorization

solution = solve(problem;
    abstol = 1e-9, reltol = 1e-7,
    dtmin = 1e-14,
    progress = true,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,   # opaque wood: deposit at surface
    min_thickness = 0.5e-3,                 # terminate at 0.5 mm
    handle_depletion = false,               # ALE alone handles recession here
    integrator = KenCarp4(linsolve = SparspakFactorization()),
)
```

Under ALE the default Jacobian is `Structured(strategy = ApproxSolve())`, which
drops the dense ALE mesh-velocity coupling (a Schur-style approximation that
still converges every term through Newton). For the exact analytical Jacobian,
including that coupling, pass `DirectSolve`:

```julia
solution = solve(problem;
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    jacobian = Structured(strategy = DirectSolve()),
    integrator = KenCarp4(linsolve = SparspakFactorization()),
)
```

Both follow the same trajectory; the exact Jacobian costs several times more
per step (Technical Reference §15.10.2; see the Solver/Jacobian configuration
chapter for choosing strategies).

### What to look at

The solver records a full boundary-flux decomposition in `extras`:

```julia
mlr           = solution.extras.surface_mass_flux           # actual MLR at surface [kg/(m²·s)]
q_absorbed    = solution.extras.surface_heat_flux_absorbed  # absorbed radiation [W/m²]
q_reradiation = solution.extras.surface_heat_flux_reradiation
q_convection  = solution.extras.surface_heat_flux_convection
q_conduction  = solution.extras.surface_heat_flux_conduction
q_net         = solution.extras.surface_heat_flux_net
ε_eff_front   = solution.extras.surface_emissivity          # tracks virgin -> char
k_eff_front   = solution.extras.surface_conductivity
```

Two distinct mass-loss measures are available and worth distinguishing:

- **`solution.extras.surface_mass_flux`** — gas actually leaving through the
  surface, summed over species. This is the physical MLR.
- **`solution.extras.gas_generation_rate`** — total volumetric gas production
  integrated over the domain. In quasi-steady burning the two agree; during
  transients they differ because generated gas may be stored in pores before
  escaping.

Component-specific generation rates require post-processing. The script
re-evaluates the reactions at each saved state and weights by the product mass
yield. Note the access pattern for internal helpers and for the reaction's
`.products` field:

```julia
NC, NR = 14, 6
rates, heats = zeros(NR), zeros(NR)
for j in 1:n_cells
    ξ = ntuple(k -> solution.ξ[k, j, i], NC)
    Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, solution.T[j, i], douglas_fir)
    # reaction 1: .products is ((comp_idx1, yield1), (comp_idx2, yield2));
    # products[2][2] is the second product's mass yield (the gas)
    gas1_gen += douglas_fir.reactions[1].products[2][2] * rates[1] * cell_vol
end
```

`evaluate_reactions!` lives in the `Pyrolysis.Physics` submodule (not a
top-level export); reach it via the qualified path as shown, or via
`using Pyrolysis.Internal`. A `Reaction` stores its stoichiometry in two
fields: `.reactants`, an `NTuple` of `(component_idx, reaction_order)` pairs,
and `.products`, an `NTuple` of `(component_idx, mass_yield)` pairs. For a
reaction with two products, `products[2][2]` is the second product's mass yield
and `products[1][2]` the first; the per-product yields sum to `1.0`. This is an
internal API and may change; prefer the solver's own `surface_mass_flux` /
`gas_generation_rate` when they suffice.

The surface-temperature solve is the headline diagnostic here: in cell-centered
finite volume, using `T_cell` directly for the `T⁴` boundary causes a ~20%
energy imbalance, so the solver Newton-iterates for the true `T_s`. For
post-processing, read it directly off the returned solution as
`solution.surface_T` (and the back face as
`solution.extras.back_temperature`) — these are populated from that
boundary-temperature Newton solve. At the lower level the solved value is cached
per boundary in the final mesh's `boundary_states` (`solution.problem.mesh`) and
can be queried with `get_surface_temperature(mesh)` /
`get_back_face_temperature(mesh)` on that mesh (exported
from the `Pyrolysis.Discretization` submodule), but for normal use the
`solution` fields are the intended API.

---

## 14.3 Grid-convergence study (`convergence_douglas_fir.jl`)

### What it models

The same Douglas fir cone case, solved at six spatial resolutions
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
function build_and_solve(n_cells::Int)
    # ... build material (13 components, 6 reactions) and BCs ...
    mesh = create_uniform_mesh(12.0e-3, n_cells, Val(13);
        T_initial = 300.0,
        ξ_initial = (530.0, 0,0,0,0,0, 31.8, 0,0,0,0,0,0))   # DF + 6% moisture
    problem = PyrolysisProblem(mesh = mesh, material = material, bc_set = bc_set,
        T_initial = z -> 300.0, ξ_initial = ξ_initial_funcs,
        tspan = (0.0, 600.0), dt_initial = 1e-3,
        output_times = collect(range(0.0, 600.0, length = 3601)))
    return solve(problem;
        abstol = 1e-9, reltol = 1e-7, dtmin = 1e-14, progress = false,
        use_ale = true, radiation_model = SURFACE_ABSORPTION, min_thickness = 0.5e-3,
        handle_depletion = false,
        integrator = KenCarp4(linsolve = SparspakFactorization()))
end
```

### How to run it

The grid sizes are independent, so the study parallelizes trivially with
`Distributed.pmap`. The script loads packages on all workers and calls
`pmap(build_and_solve, GRID_SIZES)`:

```bash
julia convergence_douglas_fir.jl            # sequential (1 worker)
julia -p auto convergence_douglas_fir.jl    # parallel (all cores)
julia -p 4 convergence_douglas_fir.jl       # 4 workers
```

```julia
using Distributed
@everywhere begin
    using Pyrolysis, LinearSolve, Sparspak, OrdinaryDiffEq
    using OrdinaryDiffEqSDIRK: KenCarp4
    using LinearSolve: SparspakFactorization
end
solutions = pmap(build_and_solve, [480, 240, 120, 60, 30, 15])
```

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
package's `parameter_sweep` helper (see §14.8.4 and the Parameter Sweeps
chapter); this script uses raw `pmap` because each point needs a fully rebuilt
problem.

---

## 14.4 Pressure-treated wood (`08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl`)

### What it models

Cone pyrolysis of fire-retardant pressure-treated wood (PTW) with moisture, at
50 kW/m². Structurally similar to Douglas fir but with a different component
ordering (water and water vapor come first), four solid intermediates instead of
five, and ThermaKin-style temperature-limited reactions. It exercises the same
Darcy–Fick + surface-absorption + ALE stack and is a good template for adapting
a ThermaKin `.cmp` material file into Pyrolysis.jl.

### Material highlights

Eleven components: `WATER` (liquid), `WATER_VAPOR` (gas), `PTW`/`PTW1`/`PTW2`/
`PTW3`/`Char` (solids), and `PTW_gas1`…`PTW_gas4` (gases). The liquid water uses
a large pseudo-density so it occupies negligible volume (bound water absorbed
into cell walls, not a separate phase):

```julia
LiquidComponent(:WATER;
    ρ = 1e5,           # pseudo-density [kg/m³] (matches ThermaKin SWELLING: 0 behavior)
    c = water_cp,
    k = 0.145,
    ε = 0.81,
    κ = 0.0,           # liquid blocks Darcy flow
)
```

This is the practical way to model "non-swelling" species: the swelling factor
defaults to `γ = 1` for solids/liquids, but a very large `ρ_j` keeps the volume
fraction near zero so the bound water does not inflate the bulk volume. (For a
truly volumeless gas you would instead set `γ = 0`; Technical Reference §10.)

### ThermaKin-to-Pyrolysis sign mapping

PTW reactions illustrate the heat-of-reaction sign flip between the codes.
ThermaKin publishes `h > 0 = exothermic`; Pyrolysis.jl stores `h > 0 =
endothermic`. So a ThermaKin `HEAT: -6.82E3 J/kg` (exothermic) is entered as
`h = +6.82e3`, while ThermaKin `HEAT: 1.82E5` (its endothermic) becomes
`h = -1.82e5`:

```julia
Reaction(:water_evaporation, 1 => 2; A = 15532.2434, E = 43528.868,
         h = 2.78e6, n = 1.0),                     # evaporation: endothermic
Reaction(:ptw_step1, 3 => (4, 8); A = 19175814.875, E = 106116.6,
         h = 6.82e3, n = 1.0, yields = (0.73075, 0.26925)),
Reaction(:ptw_step3, 5 => (6, 10); A = 179027.4062, E = 100963.197,
         h = -1.82e5, n = 1.0, yields = (0.81744, 0.18256)),
```

When porting from ThermaKin, always re-read Technical Reference §6 / nomenclature
note H1 and flip the sign deliberately. ThermaKin `TEMP LIMIT` thresholds are not
needed: Arrhenius rates fall off naturally below the onset temperature, and the
smooth `tanh` gates (`T_min`/`T_max` on `Reaction`, defaults `0.0`/`Inf`) provide
the same effect when a hard cutoff is required.

### Solve and outputs

The solve call is the same shape as Douglas fir (ALE + `SURFACE_ABSORPTION` +
Sparspak). The script writes two CSVs — a per-time summary
(`ptw_results_50kWm2.csv`) and a spatial dump (`ptw_spatial_50kWm2.csv`) sampling
`solution.T` and `solution.ξ` at selected times and cell positions — which is the
pattern to copy when exporting profiles for external plotting or for
inter-code comparison. The companion `ptw_comparison.jl` parses ThermaKin text
output and overlays it.

---

## 14.5 Darcy–Fick gas transport demo (`08_darcy_fick_demo.jl`)

### What it models

A minimal demonstration of pressure-driven gas transport in a charring material
under 75 kW/m². The point is the transport physics, not a specific real material:
a single decomposition `virgin → 0.20 char + 0.80 gas`, with the solid matrix
carrying permeability so generated gas builds internal pressure and is pushed
toward the surface by Darcy flow (Technical Reference §9):

- Darcy–Fick flux: `J_j = (ξ_j/φ) v_g + J_diff` (intrinsic density advected by the superficial velocity)
- Darcy velocity: `v_g = −(κ/μ) ∇P`
- Ideal-gas pressure in pores: `P = Σ_g ξ_j R_g T / M_j / φ`

### Enabling Darcy flow

There is no flag for Darcy mode — it is selected automatically when any
component carries `κ > 0`:

```julia
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, ε = 0.9, α = 5.0,
                        κ = 1e-15,            # tight virgin structure
                        ρ_intrinsic = 1500.0) # cell-wall (skeletal) density → φ ≈ 0.67
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, ε = 0.95, α = 50.0,
                        κ = 1e-12,            # porous char
                        ρ_intrinsic = 1950.0) # carbon skeletal density → φ ≈ 0.92
gas    = GaseousComponent(:volatiles; M = 0.044, c = 1200.0, k = 0.03, λ = 2e-5)
air    = GaseousComponent(:air;       M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)

material = Material(
    name = :CharringWood,
    components = (virgin, char, gas, air),
    reactions = (Reaction(:pyrolysis, 1 => (2, 3);
                          yields = (0.20, 0.80), A = 1.0e13, E = 150e3, h = 300e3,
                          T_min = 400.0, T_max = 1500.0),),
    mixing = MaterialMixing(k = MixingSpec(PARALLEL, 0.5), κ = MixingSpec(WEIGHTED, 0.5)),
)

# The pore space must not start as a vacuum: the solids declare skeletal
# (intrinsic) densities so the material has a real pore volume, and an inert
# background air component is initialized to 1 atm partial pressure with the
# single-species inverse of the pressure closure:
ξ_air = background_gas_concentration(material, :air, (500.0, 0.0, 0.0, 0.0), 300.0)
# ξ_initial = (virgin = 500.0, char = 0.0, volatiles = 0.0, air = ξ_air)
```

A Darcy–Fick run started with an evacuated pore space (all gas ξ = 0) against a
101325 Pa boundary spuriously advects early product gas *inward* down the
absolute-pressure slope; the background-air initialization (Technical Reference
§9.5.4) makes the initial state mechanically neutral.

The reaction here uses temperature gates `T_min = 400` and `T_max = 1500`,
implemented as smooth `tanh` ramps (Technical Reference §6.2) — handy when you
want a reaction to be inactive outside a band.

### Solve

```julia
solution_darcy = solve(problem;
    use_ale = true,
    min_thickness = 0.5e-3,
    progress = true,
)
```

The script suggests a comparison experiment: rebuild the same material with all
`κ = 0.0` (diffusion-only Fickian transport) and overlay the MLR curves — Darcy
flow typically shows faster initial gas release and lower peak internal pressure
because bulk flow drains the pores that diffusion alone would leave pressurized.

### What to look at

```julia
peak_mlr   = maximum(solution_darcy.mass_loss_rate)            # MLR per current area
initial_m  = solution_darcy.extras.total_mass[1]
final_m    = solution_darcy.extras.total_mass[end]
println("Mass released: ", round((1 - final_m/initial_m) * 100, digits = 1), " %")
```

Use `solution.mass_loss_rate` (a top-level field, MLR per the current area
`A(t)`) or `solution.extras.surface_mass_flux` for the surface flux; index both
by `solution.t`. Mass conservation across the run is the headline check — the
released mass should equal the integrated surface efflux to integrator tolerance
(Technical Reference §16).

---

## 14.6 Moisture heat-of-sorption study (`moisture_hos_comparison.jl`)

### What it models

A controlled study of how the *form* of the moisture heat of sorption affects a
Douglas fir cone simulation. Two HoS models are compared at four initial moisture
contents (3, 6, 9, 12% dry basis):

1. **Exponential (state-dependent):** `HoS(MC) = (a + b·exp(−c·MC))·1e3` J/kg —
   the energy cost rises as moisture depletes (tightly-bound water).
2. **Constant average:** a single value equal to the analytical integral average
   of the exponential model over `[0, MC₀]` — same total energy, uniform in time.

This isolates a single modeling choice (8 runs = 4 MC × 2 models) and is the
template for any "does this property form matter?" investigation.

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
    avg_kJ = compute_average_hos_kJ(HOS_A, HOS_B, HOS_C, mc_db * 100.0)
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

Recovery of unknown material parameters (component thermal conductivities) by
fitting simulated back-surface temperature traces to experimental
cone-calorimeter data at multiple heat-flux levels. This is the practical
calibration workflow: wrap the forward solve in a loss function, minimize it over
the parameters, and quantify the uncertainty of the fit.

### The pattern

The driver is built from three pieces (the heavy lifting is in
`wp_objective.jl`, included at the top):

1. **`build_and_solve(θ, scenario; stage)`** — builds a problem from a parameter
   vector `θ` and a scenario (flux level, BCs, thickness), then `solve`s it. This
   is the only place the package API appears; the rest is optimization
   plumbing. It returns a `PyrolysisSolution`, from which the fit reads
   `sol.extras.back_temperature` (and `sol.surface_T` for the FSRI front-surface
   panels).

2. **A loss function** that, for each scenario, interpolates the simulated back
   temperature onto the experimental time grid and accumulates a `σ`-weighted
   mean-square residual:

   ```julia
   T_sim = interpolate_at(t_exp[mask], sol.t, sol.extras.back_temperature)
   r     = T_sim .- T_back_exp[mask]
   loss += 0.5 * mean((r ./ σ) .^ 2)
   ```

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
difference Hessian of the loss (Gauss–Newton, default) or by residual bootstrap
(`--bootstrap=N`). The Hessian path uses `FiniteDiff.finite_difference_hessian`
and the weighted-NLS covariance `Σ_θ = (2·L(θ̂)/(N−p))·H⁻¹`. Standard errors are
`sqrt.(diag(Σ))`.

### How to run it

```bash
julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl
julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --stage=3
julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --prior-weight=0
julia -t auto --project=. examples/inverse_analysis/inverse_white_pine.jl --bootstrap=200
```

Each invocation writes a timestamped results directory with the per-evaluation
trace, the best parameters (`best.json` with standard errors and the correlation
matrix), and fit-vs-experiment plots.

### Gradient-based alternative

This example uses a derivative-free optimizer on top of bare `solve` calls. If
you prefer gradient-based optimization, Pyrolysis.jl exposes
`forward_sensitivity` and `adjoint_sensitivity` (extension-attached; load
`SciMLSensitivity` and `Zygote`). They follow the same wrap-the-solve pattern via
a `p_inject` closure that builds a perturbed problem from `θ`:

```julia
using Pyrolysis, SciMLSensitivity, Zygote

p_inject = (base, θ) -> rebuild_problem_with(base, θ)   # your problem builder

# Forward mode (small parameter counts): ForwardDiff Jacobian of an output_fn
output, ∂out_∂θ = forward_sensitivity(base_problem, θ, p_inject;
                                      output_fn = sol -> sol.surface_T[end])

# Adjoint mode (large parameter counts, scalar loss): Zygote gradient
loss, grad = adjoint_sensitivity(base_problem, θ, p_inject,
                                 sol -> sum(abs2, sol.extras.back_temperature .- data))
```

> Caveat: warm starts, non-smooth lateral-shrinkage laws, and depletion
> parking introduce discontinuities that can confuse
> gradient-based methods; the derivative-free NelderMead route is robust to them.
> See the Sensitivity Analysis chapter and Technical Reference §17.

---

## 14.8 Sensitivity studies

Four scripts under `examples/sensitivity_analysis/` (plus the inter-code
comparison `df_sensitivity_comparison.jl`) build a full uncertainty-quantification
pipeline around the Douglas fir cone case.

### 14.8.1 Forward uncertainty propagation (`sensitivity_douglas_fir.jl`)

Propagates uncertainty in 16 thermophysical inputs (instrument biases, heat-
capacity regression coefficients, heat-of-sorption parameters, and five heats of
reaction) through the simulation using Sobol quasi-random sampling. The inputs are
grouped into seven correlated blocks, each sampled jointly from a multivariate
Normal or scaled Student-t distribution with the appropriate covariance.

The structure mirrors §14.3 and §14.7: a single `build_and_solve(θ::Vector)`
(here `θ` has length 16) that maps the sampled parameters into property closures
and reaction heats, then solves with the standard ALE + surface-absorption stack.
The 8192 samples are evaluated in parallel:

```julia
using Distributed
@everywhere begin
    using Pyrolysis
    using Pyrolysis: compute_surface_mass_flux, effective_gas_transfer
    using OrdinaryDiffEq, LinearSolve, Sparspak
    using OrdinaryDiffEqSDIRK: KenCarp4
    using LinearSolve: SparspakFactorization
end

samples = generate_samples(CORR_GROUPS, 8192)        # 16 × 8192 parameter matrix
results = pmap(j -> extract_qois(build_and_solve(samples[:, j])), 1:8192)
```

```bash
julia -p 8 sensitivity_douglas_fir.jl    # 8 parallel workers
```

Each run is reduced to scalar QoIs (peak MLR, time to ignition, surface/back
temperature rise times, combustible-only burning metrics, etc.) and a few
curves; the script then computes Standardized Rank Regression Coefficients
(SRRC) for the sensitivity ranking and prediction-envelope percentiles. Note the
qualified imports of internal helpers (`compute_surface_mass_flux`,
`effective_gas_transfer`) used to reconstruct a combustible-only MLR by summing
only over the fuel-gas species — an example of reaching into the
`Pyrolysis`/`Pyrolysis.Internal` surface for a post-processing computation the
solver does not provide directly.

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
run. The variance-decomposition theory is Technical Reference §17.

### 14.8.3 Inter-code comparison (`df_sensitivity_comparison.jl`)

A pure post-processing/plotting script that parses a ThermaKin text output file
and a Pyrolysis.jl results CSV and overlays MLR, surface/back temperature, mass,
and thickness. It contains no `solve` call — it is the template for validating
Pyrolysis.jl against a reference code once you have both result files (Technical
Reference §18).

### 14.8.4 The `parameter_sweep` helper

For sweeps that do not need a fully bespoke driver, the package provides
`parameter_sweep(problem, grid; modify, extract)`. The `modify` closure builds a
fresh problem from each grid point (it must not mutate the base), and `extract`
reduces each solution to a small summary so workers ship little data back:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]

function modify(prob, p)
    new_rxn = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 200e3)
    new_mat = Material(name = :M, components = prob.material.components,
                       reactions = (new_rxn,))
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
                            T_initial = z -> 300.0, tspan = prob.tspan,
                            output_times = prob.output_times)
end

extract(sol) = (peak_T = maximum(sol.surface_T), final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```

It runs in parallel via `Distributed.pmap` when `nworkers() > 1`, otherwise
sequentially. See the Parameter Sweeps chapter for the full contract.

---

## 14.9 Literature property comparison (`literature_comparison/literature_comparison.jl`)

### What it models

Ten runs of the Douglas fir cone case (one baseline plus nine variants), each
substituting one category of property — virgin `c_p`, char `c_p`, heats of
reaction, or heat of sorption — with a value sourced from the literature, holding
everything else fixed. It quantifies how much each property choice moves the
predicted MLR and temperatures, and is a model for any "which published property
set should I trust?" study.

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

const VIRGIN_CP_DISPATCH = Dict(:baseline => baseline_virgin_cp,
                                :harada   => harada_virgin_cp,
                                :dupont   => dupont_virgin_cp)
# build_and_solve(spec) looks up the function by key, wraps it in CallableProperty,
# builds the material, and solves with the standard ALE + SURFACE_ABSORPTION stack.
```

This symbol-dispatch idiom is the recommended way to vary functional properties
across a distributed sweep, since `Function` values do not serialize cleanly to
workers.

### How to run it

```bash
julia -p 10 --project=. examples/literature_comparison/literature_comparison.jl
julia --project=. examples/literature_comparison/literature_comparison.jl   # serial
```

### What to look at

Per-run CSVs (`{label}_results.csv`: time, MLR, surface/back temperature, total
mass, thickness) and a summary QoI table (`literature_comparison_qois.csv`). The
comparison reveals which property categories dominate the spread in predicted
burning behavior — typically the heats of reaction and the high-temperature
char/virgin `c_p` — which feeds back into prioritizing measurement effort.
Technical Reference §18 frames this as validation against alternative inputs.

---

## 14.10 Common patterns across the examples

Reading the examples together, a few idioms recur and are worth internalizing:

- **Build, do not mutate.** For sweeps and parallel work, wrap material/problem
  construction in a function that takes the varied parameters and returns a fresh
  `PyrolysisProblem`. The immutable `ProblemDef` is safe to share across workers;
  the mutable `Workspace` is allocated per worker by `solve`.
- **Temperature-dependent properties are closures.** `CallableProperty(T -> ...)`
  for `c_p(T)`, `k(T)`; `StateDependentProperty((T, ξ) -> ...)` for properties or
  reaction heats that depend on local composition (e.g. heat of sorption). A bare
  number becomes a `ConstantProperty`.
- **Combine BCs with `+`.** `RadiativeFluxBC(...) + RadiativeBC(...) + ConvectiveBC(...)`
  forms a `CombinedThermalBC`; the same works for mass BCs
  (`DiffusiveMassBC(...) + DirichletPressureBC(...)`).
- **Modes are inferred, then opted into.** Darcy–Fick turns on whenever a
  component has `κ > 0`. ALE and the radiation model are explicit `solve`
  keywords (`use_ale`, `radiation_model`). Variable cross-section (the `χ`/WithChi
  mode) turns on when the `Material` is given a `lateral_shrinkage_law`
  (a function `χ̄ -> A/A_0`, e.g. `chi -> 1 - 0.2chi`). See the Building a Problem
  chapter and Technical Reference §10.
- **Read the right MLR.** `solution.extras.surface_mass_flux` is the physical
  surface efflux; `solution.extras.gas_generation_rate` is volumetric production;
  `solution.mass_loss_rate` is per current area. They coincide only at
  quasi-steady state.
- **Use the boundary-flux extras for energy accounting.** `surface_heat_flux_*` and
  `surface_emissivity`/`surface_conductivity` in `extras` give the surface energy
  balance without re-deriving it; under `BEER_LAMBERT` add
  `extras.volumetric_absorption`, the in-depth deposit the face fluxes exclude.

### Caveats to keep in mind

- **`P1_QUASI_STEADY` radiation is experimental** and quarantined in
  `Pyrolysis.Experimental`; the core `solve` rejects it as a `radiation_model`.
  The shipped examples use `NO_RADIATION` (flux via BC only) or `SURFACE_ABSORPTION`;
  `BEER_LAMBERT` is one keyword away in `06_cone_calorimeter.jl`.
- **The boundary energy balance has no surface gas-enthalpy term**, because
  `S_conv` already carries the gas enthalpy (a surface term would double count
  it; §5.4, Technical Reference §7).
- **The cone examples' parameter sets are pending re-validation against data**;
  treat their experimental and ThermaKin overlays as indicative comparisons rather
  than tight regression targets (Technical Reference §18.8.1). Parameters
  calibrated in a solver that adds a surface gas-enthalpy term on top of `S_conv`
  predict a higher surface temperature and MLR here (Technical Reference §1.8).

For deeper troubleshooting (stiff solves, ledger imbalance, ALE mesh tangling,
radiation surprises), see the Troubleshooting and Performance chapter; for the
full keyword reference, see the Running the Solver and API Reference chapters.
