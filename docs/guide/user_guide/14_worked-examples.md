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
| `07_douglas_fir_cone.jl` | Wood, multi-step kinetics, Darcy–Fick, moisture | ALE, SURFACE_ABSORPTION, DarcyFick | 14.2 |
| `convergence_douglas_fir.jl` | Grid-convergence study | as 07 at 50 kW/m², swept over `n_cells` | 14.3 |
| `08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl` | Treated wood + moisture, Darcy–Fick | ALE, SURFACE_ABSORPTION, DarcyFick | 14.4 |
| `08_darcy_fick_demo.jl` | Pressure-driven gas transport, permeability | ALE, DarcyFick | 14.5 |
| `moisture_hos_comparison.jl` | State-dependent heat of sorption | ALE, SURFACE_ABSORPTION, DarcyFick | 14.6 |
| `macfp/macfp_gasification.jl` | MaCFP-4 white-pine gasification targets (10/30/60 kW/m², 6/18 mm), the calibrated model of §14.7 at a constant cross section | ALE, SURFACE_ABSORPTION, threads | — |
| `inverse_analysis/inverse_white_pine.jl` (with `validate_thickness.jl`, `plot_diagnostics.jl`, `plot_stage2_theta.jl`) | Parameter estimation from data | ALE, optimization | 14.7 |
| `sensitivity_analysis/sensitivity_douglas_fir.jl`, `sensitivity_analysis/sensitivity_sobol.jl`, `sensitivity_analysis/compare_degrees.jl`, `df_sensitivity_comparison.jl` | Uncertainty / Sobol / inter-code | parameter sweep, post-process | 14.8 |
| `literature_comparison/literature_comparison.jl` | Property-variant comparison | as 07 at 50 kW/m², swept over property sets | 14.9 |

---

## 14.1 Cone calorimeter: PMMA (`06_cone_calorimeter.jl`)

### What it models

Standard cone-calorimeter pyrolysis of poly(methyl methacrylate) (PMMA), a
non-charring thermoplastic that depolymerizes almost entirely to gaseous MMA
monomer. The sample is a 5.8 mm slab heated from the top by a radiant cone
heater; it recedes as the material decomposes and gas escapes. This is the
canonical "your first cone run" workflow: one solid that becomes a small char
residue plus a gas, temperature-dependent `c_p` and `k`, combined thermal
boundary conditions, a moving (ALE) mesh that follows the receding surface, and
cell deactivation to remove the surface cells that burn out. The physics is
derived in Technical Reference §3 (governing equations), §6 (kinetics), §7
(energy assembly), §8 (radiation), §10–§11 (volume change and the moving mesh)
and §14 (deactivation).

### Material

The material has three components — virgin PMMA, char, and MMA gas — and one
endothermic decomposition reaction. Temperature-dependent properties are
supplied as `CallableProperty` closures (piecewise functions of `T`):

```julia
using Pyrolysis

# Piecewise temperature-dependent heat capacity and conductivity (break at 378 K).
# The fits are measured on ~300–700 K and held constant outside that range: the
# low-T cp line is negative below 204 K and the high-T k line crosses zero at 872 K.
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T))
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))

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

- **`α` (mass-basis absorption coefficient, m²·kg⁻¹).** Only consumed when
  `radiation_model = BEER_LAMBERT`; with `SURFACE_ABSORPTION` or `NO_RADIATION`
  it is ignored (all incident flux enters the surface energy balance). See
  Technical Reference §8.
- **`h` sign.** Pyrolysis.jl uses the storage convention `h > 0 endothermic`
  (the volumetric source is `Q_rxn = −Σ h_r r_r`). ThermaKin publishes the
  opposite sign (`h > 0` exothermic); Gpyro and FDS use the same sign as
  Pyrolysis.jl (Technical Reference §6.9 and the nomenclature overload note
  H1). PMMA depolymerization is endothermic, hence `h = +820e3`.
- **`yields`.** For a single-reactant reaction with two products, the implicit
  reactant coefficient is `−1` and the product yields must sum to `1.0`
  (mass conservation, enforced at construction).
- **`mixing`.** `MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL))`
  selects the conductivity mixing rule (`WEIGHTED` with `β = 0.5`, a 50/50
  parallel/series blend) and the gas-transfer mixing rule (`PARALLEL`). The
  third slot, permeability `κ`, defaults to `MixingSpec(WEIGHTED, 0.5)`.
  Mixing rules are derived in Technical Reference §5.

### Mesh

A uniform mesh is created with `create_uniform_mesh(L, n_cells)`. The mesh
holds the geometry only — no component count and no initial state:

```julia
sample_thickness = 5.8e-3   # 5.8 mm
n_cells = 72                # ~0.08 mm cells

mesh = create_uniform_mesh(sample_thickness, n_cells)
```

The initial temperature and concentrations are the problem's `T_initial` and
`ξ_initial` (below), functions of the height `z` that `setup` evaluates at the
cell centres. `ξ_initial` gives the initial mass concentrations `ξ_j`
(kg·m⁻³) of the components, not mass fractions. The cross section defaults to
1 m² (`cross_section` keyword); in 1D it is a gauge normalization for
per-area outputs (Chapter 6).

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
  absorbed fraction (default `nothing` → the composition-weighted surface
  absorptivity of the cell under the face, which is the emissivity for a
  component without its own `absorptivity`, by Kirchhoff). With
  `SURFACE_ABSORPTION` or `NO_RADIATION` the absorbed flux is applied at the
  surface through the boundary energy balance; with `BEER_LAMBERT` it is
  attenuated in depth.
- **`RadiativeBC(; ε = nothing, T_inf, F = 1.0)`** is a *loss* term
  `q = Fεσ(T_∞⁴ − T_s⁴)`, not a heat input. With `ε = nothing` it uses the
  mixture emissivity of the cell under the face (which tracks the virgin→char
  transition).
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
    progress = true,              # draw a progress bar
    use_ale = true,               # moving mesh follows the receding surface
    min_thickness = 0.01e-3,      # stop when the sample is thinner than 0.01 mm
    deactivation = true,          # remove burned-out surface cells
)
```

Notes on the solve options (full reference: Chapter 8 and the `solve`
docstring):

- **`use_ale = true`** adds the node positions and velocities to the state
  (the velocities as algebraic balance rows) so the mesh follows the solid
  (Technical Reference §11). Without it the mesh is fixed and the surface does
  not recede.
- **`radiation_model`** is omitted here, so it defaults to `NO_RADIATION`.
  The incident flux of `RadiativeFluxBC` then enters the surface energy
  balance exactly as under `SURFACE_ABSORPTION`: this script and its
  surface-absorption variant (§14.1.1) solve the same equations (their
  surface temperatures agree to the last digit on the first 300 s, measured).
  To distribute the flux in depth, pass `radiation_model = BEER_LAMBERT`
  (PMMA is semi-transparent; the `α = 1.47` values then matter).
- **`deactivation = true`** removes a surface cell once its condensed mass
  falls below 1 % of its initial condensed mass (the default `threshold`),
  in an event located inside the step (Chapter 10, Technical Reference §14).
  PMMA leaves 0.2 % char, so every cell can burn out, and `setup` installs the
  equidistribution correction, which moves the mesh relative to the solid so
  that the layer thins evenly and its cells go together. The solve ends with
  outcome `:burned_out` when the last cell goes, or `:min_thickness` below
  `min_thickness`. At 25 kW/m² neither happens within the 1600 s: the run
  ends with `:Success` at 1600 s, with 0.42 mm and 7 % of the mass left and no
  deactivation event (measured).
- **Gas transport.** The solids keep the default gas-transfer coefficient
  `λ = 0` under Fickian transport, so `solve` warns that the gas transport is
  degenerate and the integrator rejects many steps. A small `λ` on the solids
  (`λ = 1e-6`, as in the condensed version of this deck in §12.9) removes
  the warning (§15.5).

### What to look at

Outputs live on the returned `PyrolysisSolution`:

```julia
times    = solution.t                                    # saved times [s]
T_surf   = solution.surface_T                            # surface temperature [K]
T_back   = solution.extras.back_temperature              # back-face temperature [K]
total_m  = solution.extras.total_mass                    # mass per initial area [kg/m²]
gen_rate = solution.extras.gas_generation_rate           # gas the reactions produce, per initial area [kg/(m²·s)]
thick    = solution.extras.thickness_history             # current thickness [m]

println("Return code: ", solution.retcode, ", outcome: ", solution.outcome)
println("Peak surface T: ", round(maximum(T_surf), digits = 1), " K")
println("Mass loss: ", round((1 - total_m[end]/total_m[1]) * 100, digits = 1), " %")
println("Deactivation events: ", solution.statistics.invocations)
```

`surface_T` and `back_temperature` are the solved face temperatures, not
cell-centre temperatures: at every evaluation the face temperature is the
root of the face energy balance, found by a Newton iteration, which matters
for `T⁴` radiative BCs (Technical Reference §12). After a deactivation they
are the temperatures of the face bounding the live cells, and the dead cells
read `NaN` in `solution.T` and `solution.ξ` (Chapter 9). The
script prints a side-by-side comparison with ThermaKin values at t ≈ 600 s
(hard-coded in the script) and overlays ThermaKin traces only when
`docs/thermakin_output/` is present; the repository does not ship those traces.

If you have `Plots` loaded, the package's plot recipe gives quick views:

```julia
using Plots
plot(solution; kind = :surface_T)   # surface temperature vs time
plot(solution; kind = :mlr)         # mass-loss rate, solution.mass_loss_rate
plot(solution; kind = :heatmap)     # T over time and cell index
```

> Caveat: this example's parameter set is pending re-validation against measured
> data (Technical Reference §18.8.1), so the ThermaKin comparison is indicative,
> not a regression target; treat the printed peak-surface-T and mass-loss numbers as a
> reference for magnitude and time scale rather than exact values to match.

### 14.1.1 Surface-absorption variant (`06_cone_calorimeter_surface_absorption.jl`)

This companion script selects the surface-only radiation model explicitly.
The absorption coefficients are set to `α = 0.0` to make explicit that they
are unused, and the solve passes `radiation_model = SURFACE_ABSORPTION`; it
writes its figures and a CSV to `output/`:

```julia
solution = solve(problem;
    abstol = 1e-8, reltol = 1e-6,
    dtmin = 1e-14,
    progress = true,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,   # all incident flux in the surface balance
    min_thickness = 0.01e-3,
    deactivation = true,
)
```

`SURFACE_ABSORPTION` adds no unknowns: the absorbed flux is a term of the
surface energy balance. `BEER_LAMBERT` adds the face intensities as algebraic
rows of the state, each reading one cell and its two faces, so its Jacobian
stays banded (Technical Reference §8). Surface absorption is the right choice
for optically thick (opaque) materials, where the optical depth of a cell,
`τ = α_eff Δz` with `α_eff = Σ_j α_j ξ_j`, is much larger than one.
Because `06_cone_calorimeter.jl` leaves the radiation model at `NO_RADIATION`,
the two scripts give the same result; to see the difference between in-depth
and surface deposition, run `06_cone_calorimeter.jl` with
`radiation_model = BEER_LAMBERT`. On this deck at 300 s the in-depth case has
a surface 38 K cooler (630 K against 668 K) and a peak mass-loss rate of
6.2 g/(m²·s) over the first 300 s against 3.9 (measured): the absorbed flux
heats a layer beneath the surface instead of the surface itself, so the
surface re-radiates less.

---

## 14.2 Cone calorimeter: Douglas fir (`07_douglas_fir_cone.jl`)

### What it models

Pyrolysis of Douglas fir wood — an opaque charring solid with moisture — under a
75 kW/m² cone. This is the flagship "everything on" example: multi-step
decomposition kinetics (five solid→solid+gas steps that march virgin wood
through four intermediates to char, plus a moisture-evaporation step), a
state-dependent heat of sorption, Darcy–Fick gas transport with
kinetic-theory-sized gas diffusion coefficients, background air in the pores,
surface-absorption radiation, and ALE mesh motion. The physics spans Technical
Reference §5–§12.

### Material highlights

The material declares 14 components — six solids (`DF`, `DF_int1`…`DF_int4`,
`Char`), one liquid (`MOISTURE`), and seven gases (`DF_gas1`…`DF_gas5`,
`WATER_VAPOR`, and the inert background `AIR`, whose initial concentration
`background_gas_concentration(douglas_fir, :AIR, ξ_condensed_initial, T_ambient)`
sets to a 1 atm partial pressure in the pores) — and six reactions. A few
patterns worth lifting:

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
problem into Darcy–Fick transport mode automatically. The solids also declare
skeletal densities (`ρ_intrinsic`, from 1500 kg/m³ for the virgin wood to
1950 kg/m³ for the char),
which give the matrix a real pore volume to hold the background air.

**State-dependent heat of sorption.** The moisture-evaporation reaction uses a
`StateDependentProperty` whose value depends on the local moisture concentration
(more energy is needed to release the remaining tightly-bound water):

```julia
const MOISTURE_IDX = 7
const ρ_DRY_BASIS = 530.0   # dry wood density for the % moisture-content basis

heat_of_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS    # moisture content [% dry basis]
    return (2450.420909987441 + 875.6441221735336 * exp(-0.1121684273316275 * MC)) * 1e3  # J/kg
end)

# ...used as the heat of one reaction:
Reaction(:moisture_evaporation, 7 => 13;
    A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0)
```

`StateDependentProperty` closures receive `(T, ξ)` where `ξ` is the full
component tuple. They are evaluated on the dual numbers of the Jacobian like
every other property, so they must be generic Julia (no `Float64`
conversions), and may be used for any property field or for a reaction heat.
See Technical Reference §4 (property functions) and §6.

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

The top combines incident flux (absorbed at the material's surface
absorptivity, by omitting `absorptivity`), re-radiation (at the material's
emissivity, by omitting `ε`), and convection; the back is adiabatic. The top
mass BC combines diffusive escape with a Dirichlet pressure (atmospheric exhaust)
so gas can both diffuse and be pushed out by Darcy flow:

```julia
top_thermal_bc =
    RadiativeFluxBC(flux = 75e3) +            # absorptivity = nothing: the material's a_eff
    RadiativeBC(T_inf = 300) +                # ε = nothing: the material's ε_eff
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
solution = solve(problem;
    abstol = 1e-9, reltol = 1e-7,
    dtmin = 1e-14,
    progress = true,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,   # opaque wood: deposit at surface
    min_thickness = 0.5e-3,                 # stop at 0.5 mm
)
```

Everything else is the default: the `Rodas5P` integrator with the
`StructuredKLU` linear solver, and the Jacobian by colored sparse forward-mode
differentiation on the pattern of the mesh connectivity (Chapter 11). There is
no Jacobian to choose. Deactivation is off: the char layer stays in place, and
the run ends at the end of `tspan` (600 s) or when the sample is thinner than
0.5 mm. To try another integrator, pass
`integrator_policy = fbdf_policy(StructuredKLU())` or
`kencarp4_policy(StructuredKLU())` (§11.3).

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
- **`solution.extras.gas_generation_rate`** — the gas the reactions produce,
  integrated over the live cells, per initial area. In quasi-steady burning the
  two agree; during transients they differ because generated gas may be stored
  in pores before escaping.

Component-specific generation rates require post-processing. The script
re-evaluates the reactions at each saved state and weights by the product mass
yield. Note the access pattern for internal helpers and for the reaction's
`.products` field (here at saved time `i`, with the cell volumes from the saved
node positions of the moving mesh):

```julia
NC, NR = 14, 6
rates, heats = zeros(NR), zeros(NR)
z = solution.z[:, i]                                 # node positions [m]
gas1_gen = 0.0
for j in axes(solution.T, 1)
    isfinite(solution.T[j, i]) || continue           # a deactivated cell reads NaN
    ξ = ntuple(k -> solution.ξ[k, j, i], NC)
    cell_vol = (z[j+1] - z[j]) * 1.0                 # cross section 1 m²
    Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, solution.T[j, i], douglas_fir)
    # reaction 1: .products is ((comp_idx1, yield1), (comp_idx2, yield2));
    # products[2][2] is the second product's mass yield (the gas)
    gas1_gen += douglas_fir.reactions[1].products[2][2] * rates[1] * cell_vol  # kg/s
end
```

`rates` are the reactant consumption rates [kg/(m³·s)]. The script itself takes
uniform cells of the current thickness, an approximation on this mesh.

`evaluate_reactions!` lives in the `Pyrolysis.Physics` submodule (not a
top-level export); reach it via the qualified path as shown, or via
`using Pyrolysis.Internal`. A `Reaction` stores its stoichiometry in two
fields: `.reactants`, an `NTuple` of `(component_idx, reaction_order)` pairs,
and `.products`, an `NTuple` of `(component_idx, mass_yield)` pairs. For a
reaction with two products, `products[2][2]` is the second product's mass yield
and `products[1][2]` the first; the per-product yields sum to `1.0`. This is an
internal API and may change; prefer the solver's own `surface_mass_flux` /
`gas_generation_rate` when they suffice.

The surface-temperature solve is the headline diagnostic here: in a
cell-centred finite-volume method the cell temperature under the face is not
the face temperature, and using it in the `T⁴` boundary terms would misstate
the radiative loss, so the face temperature is solved from the face energy
balance by a Newton iteration. The solved face temperatures are on the
solution at every saved time: `solution.surface_T` for the exposed face and
`solution.boundary_T[:bottom]` (equal to `solution.extras.back_temperature`)
for the back face, `solution.boundary_T` holding one series per boundary tag.
The script closes by checking the surface balance at the last saved time: the
flux the top condition applies, `solution.extras.surface_heat_flux_net[end]`,
against the flux conducted into the surface cell,
`solution.extras.surface_heat_flux_conduction[end]`, whose relative difference
it expects below `1e-6`, and it prints the jump between each solved face
temperature and the centre temperature of the outermost and innermost live
cells.

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
(§12.6 and §14.8.4); this script calls `pmap` over its own
`build_and_solve(n_cells)`, which builds the whole problem at each resolution.

---

## 14.4 Pressure-treated wood (`08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl`)

### What it models

Cone pyrolysis of fire-retardant pressure-treated wood (PTW) with moisture, at
50 kW/m², 12 mm on 240 cells. Structurally similar to Douglas fir but with a
different component ordering (water and water vapor come first), three solid
intermediates between the virgin wood and the char instead of four, and five
reactions taken from a ThermaKin material file. It exercises the same
Darcy–Fick + surface-absorption + ALE stack and is a good template for adapting
a ThermaKin `.cmp` material file into Pyrolysis.jl.

### Material highlights

Twelve components: `WATER` (liquid), `WATER_VAPOR` (gas), `PTW`/`PTW1`/`PTW2`/
`PTW3`/`Char` (solids, with skeletal densities `ρ_intrinsic`), `PTW_gas1`…`PTW_gas4`
(gases), and the background `AIR` (gas, initialized to 1 atm in the pores by
`background_gas_concentration`). The liquid water is bound water absorbed into
the cell walls, not a separate phase:

```julia
LiquidComponent(:WATER;
    ρ = 1e5,           # pseudo-density [kg/m³] carried over from the ThermaKin file
    c = water_cp,
    k = 0.145,
    ε = 0.81,
    α = 1000.0,        # not used with SURFACE_ABSORPTION
    κ = 0.0,           # liquid blocks Darcy flow
)
```

A liquid's swelling factor defaults to `γ = 0`, which makes it volumeless: it
contributes no strain, no volume to the bulk and no weight in the conductivity
and emissivity mixing — ThermaKin's `L SWELLING: 0` behavior for bound water.
The large pseudo-density is therefore not needed for that; it is kept to match
the ThermaKin input. Set `γ = 1` for a liquid that occupies bulk volume (a melt
layer), as solids do by default (Technical Reference §10).

### ThermaKin-to-Pyrolysis sign mapping

PTW reactions illustrate the heat-of-reaction sign flip between the codes.
ThermaKin publishes `h > 0 = exothermic`; Pyrolysis.jl stores `h > 0 =
endothermic`. So ThermaKin's `HEAT: -2.78E6 J/kg` for the evaporation and
`HEAT: -6.82E3 J/kg` for the first decomposition step (both endothermic) are
entered as `h = +2.78e6` and `h = +6.82e3`, while ThermaKin's `HEAT: 1.82E5`
(exothermic) becomes `h = -1.82e5`:

```julia
Reaction(:water_evaporation, 1 => 2; A = 15532.2434, E = 43528.868,
         h = 2.78e6, n = 1.0),                     # evaporation: endothermic
Reaction(:ptw_step1, 3 => (4, 8); A = 19175814.8751, E = 106116.6114,
         h = 6.82e3, n = 1.0, yields = (0.73075, 0.26925)),    # endothermic
Reaction(:ptw_step3, 5 => (6, 10); A = 179027.4062, E = 100963.197,
         h = -1.82e5, n = 1.0, yields = (0.81744, 0.18256)),   # exothermic
```

When porting from ThermaKin, always re-read Technical Reference §6 / nomenclature
note H1 and flip the sign deliberately (the comments in the script describe
ThermaKin's sign the other way round; the values follow the rule above).
ThermaKin's `TEMP LIMIT` lower bounds (290 K for the evaporation, 378 K for the
decomposition steps) appear in the script as comments only: Arrhenius rates
fall off naturally below the onset temperature, and the smooth `tanh` gates
(`T_min`/`T_max` on `Reaction`, defaults `0.0`/`Inf`, width `gate_width`)
provide the same effect when a cutoff is required.

### Solve and outputs

The solve call is the same as Douglas fir's (ALE + `SURFACE_ABSORPTION`,
`abstol = 1e-9`, `reltol = 1e-7`, `min_thickness = 0.5e-3`, the default
integrator). The script writes two CSVs — a per-time summary
(`ptw_results_50kWm2.csv`) and a spatial dump (`ptw_spatial_50kWm2.csv`) sampling
`solution.T` and `solution.ξ` at selected times and cell positions — which is the
pattern to copy when exporting profiles for external plotting or for
inter-code comparison. The companion `ptw_comparison.jl` (in `examples/`) parses
the ThermaKin text output `docs/PTW_Cone_output.txt` and overlays it on
`ptw_results_50kWm2.csv`, which a full run of the cone script writes.

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

The mesh is `create_uniform_mesh(0.010, 50)` (10 mm, 50 cells); the top
condition is `HeatFluxBC(75e3) + RadiativeBC(ε = 0.9, T_inf = 300.0)` with the
mass condition `DiffusiveMassBC(0.01) + DirichletPressureBC(101325.0)`, whose
pressure anchor gives the face its Darcy venting velocity; the back is
adiabatic and impermeable. The run covers 120 s with the default tolerances and
integrator.

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
released mass should equal the integrated surface efflux to integrator
tolerance. Pass a `ConservationRecorder` and audit the run (§9.8, Technical
Reference §16):

```julia
rec = ConservationRecorder()
sol = solve(problem; use_ale = true, min_thickness = 0.5e-3, callback = rec)
audit = conservation_audit(sol, rec)
audit.mass.relative, audit.energy.relative     # relative defects of the two ledgers
```

On this deck over the 120 s the run releases 38 % of the mass, and the
recorder closes the mass ledger to 2e-7 and the energy ledger to 5e-7 relative
(measured).

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

Three scripts under `examples/sensitivity_analysis/` (plus the inter-code
comparison `examples/df_sensitivity_comparison.jl`) build a full
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
    using Pyrolysis: compute_surface_mass_flux, effective_gas_transfer, get_mass_bc
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
qualified imports of internal helpers (`compute_surface_mass_flux`,
`effective_gas_transfer`, `get_mass_bc`) used to reconstruct a combustible-only
MLR by summing only over the fuel-gas species — an example of reaching into the
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
run. They need `PolyChaos`, which `examples/Project.toml` does not list: add it
to your environment before running them. The variance-decomposition theory is
Technical Reference §17.

### 14.8.3 Inter-code comparison (`df_sensitivity_comparison.jl`)

A pure post-processing/plotting script in `examples/` that parses a ThermaKin
text output file (`docs/DF_sensitivity_Cone_output.txt`, a 50 kW/m² run) and a
Pyrolysis.jl results CSV (`douglas_fir_results_50kWm2_fit.csv` in the repository
root, which `07_douglas_fir_cone.jl` writes when run with `q_external = 50e3`)
and overlays MLR, surface/back temperature, mass, and per-component MLR,
writing its figures to `output/`. It contains no `solve` call — it is the
template for validating Pyrolysis.jl against a reference code once you have
both result files (Technical Reference §18).

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
  forms a `CombinedThermalBC`; the same works for mass BCs
  (`DiffusiveMassBC(...) + DirichletPressureBC(...)`).
- **Modes are inferred, then opted into.** Darcy–Fick turns on whenever a
  component has `κ > 0`. The moving mesh, the radiation model and deactivation
  are explicit `solve` keywords (`use_ale`, `radiation_model`, `deactivation`).
  A variable cross section turns on when the `Material` is given a
  `lateral_shrinkage_law` (a function `χ̄ -> A/A_0` of the column's conversion,
  e.g. `chi -> 1 - 0.2chi`), which needs `use_ale = true`. See Chapter 7 and
  Technical Reference §10.
- **Pore air for Darcy–Fick.** In `07_douglas_fir_cone.jl`,
  `08_pressure_treated_wood_cone.jl` and `08_darcy_fick_demo.jl` the solids
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
- **The boundary energy balance has no surface gas-enthalpy term**, because the
  gas enthalpy is carried in flux form by the mass fluxes across the faces,
  the boundary face included (a surface term would double count it; §5.4,
  Technical Reference §7).
- **The cone examples' parameter sets are pending re-validation against data**;
  treat their experimental and ThermaKin overlays as indicative comparisons rather
  than tight regression targets (Technical Reference §18.8.1). Parameters
  calibrated in a solver that adds a surface gas-enthalpy term on top of the
  enthalpy the gas carries predict a higher surface temperature and MLR here
  (Technical Reference §1.8).

For deeper troubleshooting (rejected steps, ledger imbalance, crossed mesh
nodes, radiation surprises), see the Troubleshooting and Performance chapter
(Chapter 15); for the full keyword reference, see the Running the Solver and
API Reference chapters (Chapters 8 and 16).
