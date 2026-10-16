# 1. Introduction and Installation

Pyrolysis.jl is a Julia package for simulating the thermal decomposition
(pyrolysis) of condensed-phase materials exposed to high heat flux. This chapter
explains what the package is, when to reach for it, how it is organized relative
to the companion Technical Reference Guide, how to install it and its optional
extensions, and how to verify the install with a complete ten-line "hello
pyrolysis" run. By the end you will have a working environment and a mental model
for the rest of the User Guide.

This is a *User Guide*: it teaches you how to *use* the package — defining inputs,
running simulations, and post-processing results. Wherever the physics or
numerical method behind a feature matters, this guide points you to the
**Technical Reference** (cited as "see Technical Reference §N"), which derives and
justifies the mathematics.

---

## 1.1 What Pyrolysis.jl is

Pyrolysis.jl solves the coupled heat-transfer, reaction-kinetics, and
gas-transport equations governing a one-dimensional (through-thickness) slab of
decomposing material. The primary state variables are the cell-centered
temperature `T` (K) and the per-component mass concentrations `ξ_j` (kg/m³), and
the package integrates these forward in time using stiff ODE solvers from the
SciML ecosystem (OrdinaryDiffEq.jl).

The model is built on a cell-centered **finite-volume** spatial discretization
(see Technical Reference §13) with an optional **Arbitrary Lagrangian–Eulerian
(ALE)** moving mesh (see Technical Reference §11) that tracks surface recession,
swelling, and intumescence. Reactions are **Arrhenius** rate laws with smooth
temperature gates (see Technical Reference §6); heat transfer is Fourier
conduction plus gas-advection energy transport (see Technical Reference §7);
radiation, in-depth gas transport, and variable cross-section are all optional,
opt-in physics.

The design philosophy is **"a generalized model with optional physics."** A
minimal run needs only a mesh, a material (one solid plus one gas), one reaction,
and boundary conditions. From there you turn on exactly the extra physics your
problem needs — radiation absorption, pressure-driven Darcy flow, mesh motion,
lateral shrinkage, adaptivity — by setting the corresponding component fields or
`solve` keywords. Nothing you do not ask for is computed.

Pyrolysis.jl is intended for the same class of problems as **ThermaKin2Ds**,
**Gpyro**, and the **FDS** solid-phase pyrolysis model: cone-calorimeter
simulations, TGA/DSC/MCC idealizations, charring and non-charring solids,
intumescent coatings, pressure-treated wood, and moisture transport. The energy
equation is the non-conservative temperature equation of a one-temperature porous
medium: every component, pore gas included, stores sensible heat, and the
enthalpy carried by the gas flux relative to the skeleton enters through the
advection source `S_conv`. This is the ThermaKin convention (see Technical
Reference §3 and §7).

### What is not modeled

To set expectations up front, the current model does **not** support:

- Multi-dimensional (2D/3D) transport — the model is strictly 1D
  through-thickness.
- A full two-temperature gas phase (the pore gas shares the matrix temperature
  and stores sensible heat at it — see Technical Reference §3).
- Multi-reactant reactions — each reaction has a single reactant (`NRin == 1`);
  multi-step schemes are built as *sequences* of single-reactant reactions (see
  Technical Reference §6).
- Temperature- or pressure-dependent pre-exponential factors, and distributed
  activation-energy tooling.

---

## 1.2 Capabilities at a glance

| Capability | How it is enabled | Reference |
|---|---|---|
| 1D condensed-phase heat conduction | Always on | TR §7 |
| Multi-component materials (solid / liquid / gas) | `SolidComponent`, `LiquidComponent`, `GaseousComponent` | TR §4 |
| Temperature-dependent properties | `LinearProperty`, `PolynomialProperty`, `TableProperty`, … | TR §4 |
| Multi-step Arrhenius kinetics | A tuple of `Reaction`s on the `Material` | TR §6 |
| Effective (mixture) properties, mixing rules | `MixingSpec`, `MaterialMixing` | TR §5 |
| Surface / in-depth radiation | `solve(...; radiation_model = …)` + `RadiativeFluxBC` | TR §8 |
| Pressure-driven gas transport (Darcy–Fick) | Set `κ > 0` on a solid/liquid component | TR §9 |
| Volume change, swelling, intumescence | `swelling_factor` (γ) on components | TR §10 |
| Lateral shrinkage / variable cross-section | `lateral_shrinkage_law` on the `Material` | TR §10 |
| Moving mesh (surface recession) | `solve(...; use_ale = true)` | TR §11 |
| Surface-cell depletion (cell parking) | `solve(...; handle_depletion = true)` | TR §14 |
| Error indicators, r-refinement, mesh-quality metrics (standalone; not used by `solve`) | `Pyrolysis.Internal` / submodule access | TR §14 |
| Sensible-energy diagnostics log | `solve(...; diagnostics = true)` | TR §16 |
| Sensitivity analysis (forward / adjoint) | `forward_sensitivity`, `adjoint_sensitivity` (extension) | TR §17 |
| Parameter sweeps (parallel) | `parameter_sweep` | UG §12 |

A few features are explicitly **experimental** or carry caveats; this guide flags
them plainly wherever they appear:

- **P1 quasi-steady radiation** lives in `Pyrolysis.Experimental`, not in the core
  solve path. The exported radiation models you can pass to `solve` are
  `NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT`; requesting
  `P1_QUASI_STEADY` is rejected by `solve` (see Technical Reference §8).
- **Depletion handling** (`handle_depletion = true`) parks depleted cells
  without changing the state size, so it works with every Jacobian backend;
  it is a sequence of discrete state rewrites, so parameter sensitivities are
  valid only between park events (Chapter 10, Technical Reference §14, §17).
- **h-refinement / h-coarsening** (cell splitting / coarsening) lives in
  `Pyrolysis.Experimental`, and its meshes cannot be passed to `solve`. The only
  mesh adaptation `solve` performs is depleted-cell parking; error indicators,
  r-refinement, and mesh-quality metrics are standalone building blocks.

---

## 1.3 Coordinate and unit conventions

Pyrolysis.jl uses SI units throughout, and the geometric conventions below are
binding across both guides (see Technical Reference §1 and §2 for the full
nomenclature).

- The spatial coordinate is **`z`** (axial, through-thickness), in metres.
- **`z = 0` is the bottom / substrate** (boundary id 2, tag `:bottom`).
- **`z = L` is the top / exposed surface** (boundary id 1, tag `:top`).
- Heat enters from the **top**; charring materials shrink **downward**.
- Cell 1 is at `z ≈ 0` (bottom); cell `n` is at `z ≈ L` (top).
- Positive flux means transport in the **+z** direction (bottom → top).

Key state and input units:

| Quantity | Symbol | Units |
|---|---|---|
| Temperature | `T` | K |
| Mass concentration (primary species state) | `ξ_j` | kg·m⁻³ |
| Domain thickness | `L` | m |
| Heat flux | `q` | W·m⁻² |
| Density (bulk) | `ρ_j` | kg·m⁻³ |
| Heat capacity | `c_p` | J·kg⁻¹·K⁻¹ |
| Thermal conductivity | `k` | W·m⁻¹·K⁻¹ |
| Pre-exponential factor | `A_i` | s⁻¹ (first order) |
| Activation energy | `E_i` | J·mol⁻¹ |
| Heat of reaction | `h` | J·kg⁻¹ of reactant |
| Permeability | `κ` | m² |
| Cross-sectional area | `A` | m² |

> **Sign convention for the heat of reaction.** Internally Pyrolysis.jl uses
> **`h > 0` = endothermic** (the reaction absorbs heat and cools the material) and
> **`h < 0` = exothermic**. Note this is the *opposite* of the ThermaKin/Gpyro
> publishing convention (`h > 0` = exothermic). The volumetric heat source is
> `Q_rxn = −Σ h_r r_r`, so an endothermic reaction (`h > 0`) produces
> `Q_rxn < 0`. See Technical Reference §6 (overload note H1 in §2).

The state variable is **mass concentration `ξ_j`** (kg/m³), not mass fraction:
`ξ_j = Y_j ρ`. A component initially present at its bulk density `ρ_j` has
`ξ_j(0) = ρ_j`. This is why the default initial condition places the first solid
component at its reference density (see §1.8).

---

## 1.4 How this guide relates to the Technical Reference

Pyrolysis.jl ships two companion documents that share one unified nomenclature
(Technical Reference §2):

- The **Technical Reference Guide** derives the continuum equations, the
  constitutive laws, the discretization, the Jacobian, and the conservation
  diagnostics. Read it when you need to know *why* a result is what it is, or what
  approximation a feature rests on.
- The **User Guide** (this document) teaches the *workflow*: how to define
  materials, reactions, and boundary conditions; how to build and solve a problem;
  how to post-process; and how to troubleshoot.

The two are cross-referenced. Throughout the User Guide you will see "see
Technical Reference §N"; the Technical Reference chapters are numbered:
1 Introduction & Scope; 2 Nomenclature; 3 Governing Equations; 4 Materials &
Property Functions; 5 Effective Properties & Mixing Rules; 6 Reaction Kinetics;
7 Heat Conduction, Convection & Energy Assembly; 8 Thermal Radiation; 9 Gas-Phase
Mass Transport & Pressure; 10 Volume Change & Lateral Shrinkage; 11 ALE Framework;
12 Boundary & Initial Conditions; 13 Finite-Volume Discretization; 14 Adaptive
Mesh Refinement; 15 Temporal Integration & Structured Jacobian; 16 Conservation
Diagnostics; 17 Sensitivity Analysis; 18 Verification & Validation.

### Roadmap of the rest of the User Guide

| Chapter | Topic |
|---|---|
| UG 2 | Quick-start tutorial: your first cone-calorimeter run |
| UG 3 | Defining materials, components, and properties |
| UG 4 | Defining reactions |
| UG 5 | Boundary and initial conditions |
| UG 6 | Meshes and geometry |
| UG 7 | Building a problem: simulation modes |
| UG 8 | Running the solver |
| UG 9 | The solution object and plotting |
| UG 10 | Adaptive meshing options |
| UG 11 | Solver and Jacobian configuration; extensions |
| UG 12 | Parameter sweeps |
| UG 13 | Sensitivity analysis |
| UG 14 | Worked examples and tutorials |
| UG 15 | Troubleshooting and performance |
| UG 16 | API and parameter reference |

If you are new to the package, work through UG 1 (this chapter) and UG 2 first,
then dip into UG 3–8 as you build your own model.

---

## 1.5 Installation

### 1.5.1 Julia version requirements

Pyrolysis.jl requires **Julia 1.10 or newer** (this is the `julia` compat bound
declared in `Project.toml`). The package uses Julia's weak-dependency *package
extension* mechanism (Julia 1.9+) for its optional features, so a current Julia
release is required. Download Julia from <https://julialang.org/downloads/> or
install it via [`juliaup`](https://github.com/JuliaLang/juliaup), the official
version manager:

```bash
# Install juliaup, then the latest stable Julia
curl -fsSL https://install.julialang.org | sh
juliaup add release
```

Verify your Julia version:

```julia
julia> VERSION
v"1.10.0"     # or newer
```

### 1.5.2 Installing the package

Pyrolysis.jl is a standard Julia package. The recommended practice is to install
it into a **project environment** so its dependency versions are pinned in a
`Project.toml`/`Manifest.toml` pair that travels with your work.

Create or activate a project, then add the package:

```julia
using Pkg
Pkg.activate("my_pyrolysis_study")   # create/activate a local environment
Pkg.add("Pyrolysis")                 # if registered in a registry you have
```

If you are working from a local clone or a Git URL (the typical case while the
package is pre-registration), add it by path or URL:

```julia
using Pkg
Pkg.activate("my_pyrolysis_study")

# From a local checkout:
Pkg.add(path = "/path/to/Pyrolysis.jl")

# Or directly from a Git remote:
Pkg.add(url = "https://github.com/<owner>/Pyrolysis.jl")
```

You can do the same from the Pkg REPL mode (press `]` at the `julia>` prompt):

```
(my_pyrolysis_study) pkg> add /path/to/Pyrolysis.jl
```

### 1.5.3 Developing the package

If you intend to modify the source — fix a bug, add physics, or run the test
suite — use `Pkg.develop` (`dev` in the Pkg REPL) instead of `add`. This places
an editable, tracked checkout under your Julia depot (or at a path you choose):

```julia
using Pkg
Pkg.develop(path = "/path/to/Pyrolysis.jl")
```

```
(@v1.10) pkg> dev /path/to/Pyrolysis.jl
```

A `dev`-ed package picks up your local edits on the next `using Pyrolysis` (with
[`Revise.jl`](https://github.com/timholy/Revise.jl) loaded first, edits are
reflected without restarting Julia).

### 1.5.4 Project.toml and Manifest.toml

A project environment is described by two files:

- **`Project.toml`** lists the *direct* dependencies and version compatibility
  bounds (`[compat]`).
- **`Manifest.toml`** records the *exact* resolved versions of every direct and
  transitive dependency — it is the reproducibility lockfile.

To reproduce someone else's environment, place their `Project.toml` and
`Manifest.toml` in a directory and instantiate it:

```julia
using Pkg
Pkg.activate("path/to/their/environment")
Pkg.instantiate()    # downloads and resolves the exact pinned versions
```

Commit both files to version control for any study you want to be reproducible.

---

## 1.6 Optional extensions (load-on-demand)

Pyrolysis.jl keeps its core lightweight by declaring two optional capabilities —
plotting recipes and SciML sensitivity analysis — as **package extensions**. Each
extension is wired to *weak dependencies* and activates automatically the moment
you `using` the corresponding packages alongside `Pyrolysis` — you never call the
extension explicitly, you just load its trigger packages. The relevant entries
from the package's `Project.toml` are:

```toml
[weakdeps]
RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
SciMLSensitivity = "1ed8b502-d754-442c-8d5d-10ac956f44a1"
Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[extensions]
PyrolysisPlotsExt = "RecipesBase"
PyrolysisSciMLSensitivityExt = ["SciMLSensitivity", "Zygote"]
```

The table below summarizes each extension, the package(s) that trigger it, and
what it adds. Neither is required for a basic solve.

| Extension | Trigger package(s) | What it provides |
|---|---|---|
| `PyrolysisPlotsExt` | `Plots` (or `RecipesBase`) | Plot recipes for `PyrolysisSolution` (`plot(sol; kind = …)`) |
| `PyrolysisSciMLSensitivityExt` | `SciMLSensitivity`, `Zygote` | `forward_sensitivity`, `adjoint_sensitivity` |

Install the trigger packages you need into the same environment, then load them:

```julia
using Pkg
Pkg.add(["Plots", "SciMLSensitivity", "Zygote"])

using Pyrolysis
using Plots                      # activates PyrolysisPlotsExt → plot(sol; kind = …)
using SciMLSensitivity, Zygote   # activates PyrolysisSciMLSensitivityExt
```

Alternative sparse linear solvers need no extension: any LinearSolve factorization
(for example `KLUFactorization()`, or `MKLPardisoFactorize()` with Pardiso/MKL
installed) is passed directly as `linear_solver` after `using LinearSolve` (UG §11.4).

### Stub behavior when an extension is not loaded

The functions `forward_sensitivity` and `adjoint_sensitivity` are exported by
`Pyrolysis` as *stubs* with no methods. Calling one before
`using SciMLSensitivity, Zygote` raises a `MethodError` — that is the signal that
you need to load the trigger packages; add them and the concrete methods appear.
The plotting recipes have no stub: without `Plots` loaded there is no `plot`
function to call.

The two extensions are independent — neither depends on the other — so you load
only the ones a given study requires. Plotting recipes and sensitivity analysis
are covered in detail in UG §9 and UG §13.

---

## 1.7 Verifying the install

The fastest correctness check is to load the package and confirm the public
surface is present:

```julia
using Pyrolysis

@assert isdefined(Pyrolysis, :Material)
@assert isdefined(Pyrolysis, :solve)
@assert isdefined(Pyrolysis, :create_uniform_mesh)
println("Pyrolysis.jl loaded successfully.")
```

If you `dev`-ed the package and want to run its test suite (this requires the
test environment to resolve):

```julia
using Pkg
Pkg.test("Pyrolysis")
```

For a real end-to-end check, run the minimal example in §1.8 — completing it with
`retcode == :Success` is the strongest confirmation that conduction, kinetics, gas
transport, the boundary solve, and the integrator are all working together.

---

## 1.8 A minimal "hello pyrolysis" example

The script below is a complete, runnable simulation: a 1 cm slab of PMMA-like
material, modeled as one solid that decomposes into one gas, heated from the top
by a constant 50 kW/m² flux, insulated at the bottom, integrated for 60 seconds.
Every function, type, and keyword used here is part of the public API, and the
script runs as written, returning `retcode = :Success`.

```julia
using Pyrolysis

# 1. Define the material: one solid component + one gaseous product.
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86, λ = 1e-6),
        GaseousComponent(:MMA,  M = 0.100, c = 1500, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        # :virgin → :MMA, endothermic (h > 0): absorbs 870 kJ/kg of reactant.
        Reaction(:decomposition, :virgin => :MMA, A = 8.5e12, E = 188e3, h = 870e3),
    ),
)

# 2. Build a uniform mesh: 1 cm thick, 40 cells, NC = 2 components.
mesh = create_uniform_mesh(0.01, 40, Val(2))

# 3. Boundary conditions: constant heat flux on top, insulated bottom,
#    gas free to diffuse out the top, impermeable bottom.
bcs = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3),     bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(),    bottom = ImpermeableBC()),
)

# 4. Assemble the problem. T_initial defaults to z -> 300.0; ξ_initial defaults
#    to the first solid at its reference density (1190 kg/m³) and the gas at 0.
problem = PyrolysisProblem(
    mesh      = mesh,
    material  = pmma,
    bc_set    = bcs,
    T_initial = z -> 300.0,
    tspan     = (0.0, 60.0),
)

# 5. Solve and inspect.
solution = solve(problem; progress = false)
println("retcode      = ", solution.retcode)          # :Success
println("surface T(t) = ", solution.surface_T[end], " K")
println("T field size = ", size(solution.T))           # (n_cells, n_times)
```

A few notes on what each step does and what is defaulted:

- **`Material(; name, components, reactions = (), mixing = MaterialMixing(),
  lateral_shrinkage_law = nothing)`** — the components tuple is heterogeneous (a
  `SolidComponent` and a `GaseousComponent`); the reaction's symbolic
  `:virgin => :MMA` pair is resolved to component indices against the components
  tuple inside the constructor. Materials and reactions are covered in UG §3 and
  UG §4.
- **`SolidComponent(:name; ρ, c, k, …)`** requires bulk density `ρ`, heat
  capacity `c`, and conductivity `k`; `ε` (emissivity, a number or a property
  function of temperature) defaults to `0.9`, the surface `absorptivity` for
  incident radiation to `nothing` (Kirchhoff, `a = ε(T)`), the
  swelling factor `γ` to `1.0`, permeability `κ` to `0.0` (impermeable, so no
  Darcy flow), the gas-transfer coefficient `λ` to `0.0`, and the mass-absorption
  coefficient `α` to `0.0`. Any of `ρ`, `c`, `k` may be a number or a property
  function (UG §3).
- **`GaseousComponent(:name; M, c, k, λ)`** requires molar mass `M` (kg/mol),
  heat capacity `c`, conductivity `k`, and gas-transfer coefficient `λ`. Its
  density defaults to the ideal-gas law and its swelling factor `γ` to `0.0` (gas
  escapes without changing the matrix volume).
- **`Reaction(:name, reactant => product; A, E, h)`** with symbol payloads creates
  a single-reactant, single-product reaction. `n` (reaction order) defaults to
  `1.0`, `T_min` to `0.0`, `T_max` to `Inf`, `gate_width` (tanh gate ramp) to
  `1.0` K. Remember `h > 0` is endothermic.
- **`create_uniform_mesh(L, n_cells, Val(NC))`** builds the mesh. The third
  argument is the number of components *as a type parameter*, written `Val(2)` for
  two components — it must match the material's component count. Optional keywords
  are `T_initial = 300.0`, `ξ_initial = nothing`, and `cross_section_area = 1.0`.
  Geometry is covered in UG §6.
- **`BoundaryConditionSet(; thermal, mass)`** takes per-boundary NamedTuples keyed
  by `:top` and `:bottom`. Unspecified slots default to `AdiabaticBC` (thermal)
  and `ImpermeableBC` (mass). `DiffusiveMassBC()` defaults to infinite external
  conductance (`h_m = Inf`, i.e. pure internal-diffusion-limited escape). Boundary
  conditions are covered in UG §5.
- **`PyrolysisProblem(; mesh, material, bc_set, tspan, …)`**: `T_initial` defaults
  to `z -> 300.0`, `ξ_initial` defaults to the first solid component at its
  reference density (others zero), `dt_initial` to `1e-4`, and `output_times` is
  auto-generated (11–101 evenly spaced samples over `tspan`) when left empty.
- **`solve(problem; …)`** returns a `PyrolysisSolution`. Key defaults are
  `integrator = nothing` (→ `KenCarp4`, a stiff SDIRK method),
  `abstol = 1e-8`, `reltol = 1e-6`, `use_ale = false` (fixed Eulerian mesh),
  `radiation_model = NO_RADIATION`, `diagnostics = false`, `validate = true`,
  `progress = true`, and `handle_depletion = false`. Setting
  `progress = false` (as above) suppresses the progress bar, which is convenient in
  scripts and tests. The full option set is documented in UG §8.

### Reading the result

`solve` returns a `PyrolysisSolution`. The fields you will use most are:

- `solution.t` — the vector of saved times (s).
- `solution.T` — temperature field, an `n_cells × n_times` matrix (K).
- `solution.ξ` — concentrations, an `NC × n_cells × n_times` array (kg/m³).
- `solution.surface_T` — the exposed-surface (top) temperature per saved time (K).
- `solution.mass_loss_rate` — mass-loss rate per current area (both faces;
  bottom term zero for an impermeable holder) (kg·m⁻²·s⁻¹).
- `solution.retcode` — the integrator return code (`:Success`, `:Terminated`,
  `:MaxIters`, …).

Convenience accessors for the solution — `get_profile`, `get_cell_positions`, and
`interpolate_solution` — are exported by the `Solver` submodule but **not**
re-exported at the top level (to keep generic names like `get_profile` out of the
default namespace). Reach them via the `Pyrolysis.Internal` umbrella or by their
submodule path:

```julia
using Pyrolysis
using Pyrolysis.Internal           # brings get_profile, get_cell_positions, … into scope

T_profile = get_profile(solution, :T, 1)        # temperature in every cell at time index 1
ξ1_final  = get_profile(solution, :ξ1, length(solution.t))  # ξ of component 1 at the last time
positions = get_cell_positions(solution)        # cell-center z positions (m)

# Or without `using Pyrolysis.Internal`, fully qualified:
T_profile = Pyrolysis.Solver.get_profile(solution, :T, 1)
```

If you have `Plots` installed and loaded, the plotting extension provides recipes
for the common views (covered fully in UG §9):

```julia
using Plots
plot(solution)                       # default: surface temperature trace
plot(solution; kind = :mlr)          # mass-loss rate vs time
plot(solution; kind = :heatmap)      # T(z, t) heatmap
```

---

## 1.9 Public API surface and the `Internal` escape hatch

`using Pyrolysis` brings a deliberately curated set of names into scope — the
stable, user-facing API. Broadly, the top-level exports cover:

- **Materials:** `Material`, `SolidComponent`, `GaseousComponent`,
  `LiquidComponent`, `Reaction`, `DepletionLimiter`, the mixing rules
  (`PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN`, `CARMAN_KOZENY`), `MixingSpec`,
  `MaterialMixing`, and the phase predicates (`is_solid`, `is_gas`, `is_liquid`,
  `phase`, `molar_mass`).
- **Property functions:** `ConstantProperty`, `LinearProperty`,
  `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`,
  `SutherlandProperty`, `StateDependentProperty`.
- **Boundary conditions:** the thermal BCs (`HeatFluxBC`, `ConvectiveBC`,
  `RadiativeBC`, `RadiativeFluxBC`, `CombinedThermalBC`, `DirichletThermalBC`,
  `AdiabaticBC`), the mass BCs (`ImpermeableBC`, `MassFluxBC`, `ConvectiveMassBC`,
  `DiffusiveMassBC`), the pressure BCs (`DirichletPressureBC`, `NeumannPressureBC`,
  `ConvectivePressureBC`), and the container `BoundaryConditionSet`.
- **Radiation models:** `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`,
  `BEER_LAMBERT` (the `P1_QUASI_STEADY` value is exported from
  `Pyrolysis.Experimental`).
- **Geometry:** `create_uniform_mesh`, `create_stretched_mesh`,
  `create_mesh_from_nodes`.
- **Problem and solve:** `PyrolysisProblem`, `PyrolysisSolution`, `solve`,
  `parameter_sweep`.
- **Jacobian configuration:** `Structured` and its two linear-solve strategies
  (`DirectSolve`, `ApproxSolve`), `JacobianSpec`, `DenseAD_ForwardDiff` (a
  verification reference), `verify_jacobian_at`, and more (UG §11).
- **Extension stubs:** `forward_sensitivity`, `adjoint_sensitivity`.

Anything not exported at the top level remains reachable in one of two ways:

```julia
using Pyrolysis.Internal     # re-exports every submodule's public surface in one namespace
# e.g. get_profile, effective_conductivity, divergence_1d!, … all become available

# Or reach a specific submodule directly:
Pyrolysis.Solver.get_profile(solution, :T, 1)
Pyrolysis.Properties.effective_conductivity
```

The submodules are: `Materials`, `BoundaryConditions`, `Geometry`, `Properties`,
`Physics`, `Adaptivity`, `Discretization`, `Problem`, `Residual`, `Jacobian`,
`Diagnostics`, and `Solver`. The `Experimental` namespace (P1 radiation,
h-refinement/coarsening) is *not* part of `Internal` and is accessed only as
`Pyrolysis.Experimental.X`.

A couple of practical notes you will encounter later in the guide:

- The simulation-mode **trait subtypes** (`Eulerian`, `ALE`, `Fickian`,
  `DarcyFick`, `NoChi`, `WithChi`, …) are deliberately *not* exported at the top
  level — they would shadow generic names. When you need one (rarely, for
  low-level work), import it explicitly, e.g.
  `import Pyrolysis.Problem: Eulerian`. In normal use you select these modes
  indirectly through `solve` keywords and the material's fields (UG §7).
- Symbols reached through `Internal` or a submodule are **not** part of the stable
  public API and may change between minor releases. Prefer the top-level exports
  in production code; treat `Internal` as an escape hatch.

> **Caveat — API stability.** Top-level `using Pyrolysis` exports are stable.
> Submodule / `Internal` symbols can move between minor releases. Everything in
> `Pyrolysis.Experimental` is explicitly unstable.

---

## 1.10 Where the examples live

The repository ships a set of runnable example scripts under the `examples/`
directory, each with its own project environment (it adds `Plots` for
visualization). These are the best starting points for real configurations, and
UG §14 walks through them in detail. The most useful to start with:

| File | What it demonstrates |
|---|---|
| `examples/06_cone_calorimeter.jl` | Standard cone-calorimeter PMMA run: single-step kinetics, temperature-dependent properties, combined thermal BCs, gas generation |
| `examples/06_cone_calorimeter_surface_absorption.jl` | Same case with the `SURFACE_ABSORPTION` radiation model |
| `examples/07_douglas_fir_cone.jl` | Douglas fir: multi-step kinetics, surface-absorption radiation, ALE mesh tracking shrinkage, Darcy–Fick transport, `ApproxSolve` for speed |
| `examples/08_darcy_fick_demo.jl` | Pressure-driven (Darcy–Fick) gas transport with permeability |
| `examples/08_pressure_treated_wood_cone.jl` | Pressure-treated wood with moisture |
| `examples/convergence_douglas_fir.jl` | Spatial / temporal convergence study |
| `examples/moisture_hos_comparison.jl` | Moisture transport with state-dependent heat of sorption |
| `examples/ptw_comparison.jl` | Pressure-treated wood comparison |
| `examples/df_sensitivity_comparison.jl` | Sensitivity-analysis workflow |

To run an example, activate its environment and include it:

```julia
using Pkg
Pkg.activate("examples")     # the examples directory has its own Project.toml
Pkg.instantiate()            # first time only: resolve the example dependencies
include("examples/06_cone_calorimeter.jl")
```

The **Examples** page of the documentation site (source `docs/src/examples.md`)
is the authoritative catalogue of what each script covers. With the install
verified and the minimal example running, you are ready for UG §2, which builds a
full cone-calorimeter simulation step by step.
