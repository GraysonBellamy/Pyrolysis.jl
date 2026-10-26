# 1. Introduction and Installation

Pyrolysis.jl is a Julia package for simulating the thermal decomposition
(pyrolysis) of condensed-phase materials exposed to high heat flux. This chapter
explains what the package is, when to reach for it, how it is organized relative
to the companion Technical Reference Guide, how to install it and its optional
extension, and how to verify the install with a complete "hello pyrolysis" run.
By the end you will have a working environment and a mental model for the rest
of the User Guide.

This is a *User Guide*: it teaches you how to *use* the package — defining inputs,
running simulations, and post-processing results. Wherever the physics or
numerical method behind a feature matters, this guide points you to the
**Technical Reference** (cited as "see Technical Reference §N"), which derives and
justifies the mathematics.

---

## 1.1 What Pyrolysis.jl is

Pyrolysis.jl solves the coupled heat-transfer, reaction-kinetics, and
gas-transport equations governing a one-dimensional (through-thickness) slab of
decomposing material. You describe the composition by the per-component mass
concentrations `ξ_j` (kg/m³) and the temperature `T` (K); the integrator carries,
per cell, the component masses and an energy state (the temperature by default,
or the cell's enthalpy inventory), and integrates them forward in time with stiff
solvers from the SciML ecosystem (OrdinaryDiffEq.jl). Because the masses
themselves are the states, mass is conserved exactly.

The model is built on a cell-centered **finite-volume** discretization (see
Technical Reference §13) on a mesh that can move with the solid, the
**Arbitrary Lagrangian–Eulerian (ALE)** framework (see Technical Reference §11),
which tracks surface recession, swelling, and intumescence. Reactions are
**Arrhenius** rate laws with smooth temperature gates (see Technical Reference
§6); heat transfer is Fourier conduction plus the enthalpy the gas carries across
the cell faces (see Technical Reference §7); radiation, in-depth gas transport,
the variable cross-section, and burnout by cell deactivation are all optional,
opt-in physics.

The design philosophy is **"a generalized model with optional physics."** A
minimal run needs only a mesh, a material (one solid plus one gas), one reaction,
and boundary conditions. From there you turn on exactly the extra physics your
problem needs — radiation absorption, pressure-driven Darcy flow, mesh motion,
lateral shrinkage, burnout — by setting the corresponding component fields or
`solve` keywords. Nothing you do not ask for is computed.

Pyrolysis.jl is intended for the same class of problems as **ThermaKin2Ds**,
**Gpyro**, and the **FDS** solid-phase pyrolysis model: cone-calorimeter
simulations, TGA/DSC/MCC idealizations, charring and non-charring solids,
intumescent coatings, pressure-treated wood, and moisture transport. The energy
balance is that of a one-temperature porous medium: every component, the pore
gas included, stores sensible heat at the cell temperature, and the gas fluxes
carry their sensible enthalpy across the faces, leaving the exposed surface at
the surface temperature (see Technical Reference §3 and §7).

### What is not modeled

To set expectations up front, the model does **not** support:

- Multi-dimensional (2D/3D) solves — a problem is one-dimensional through the
  thickness. (The mesh and residual layers also build quadrilateral and
  hexahedral meshes, but these are building blocks, not a supported solve path.)
- A full two-temperature gas phase (the pore gas shares the matrix temperature
  and stores sensible heat at it — see Technical Reference §3).
- Multi-reactant reactions — each reaction has a single reactant (`NRin == 1`);
  multi-step schemes are built as *sequences* of single-reactant reactions (see
  Technical Reference §6).
- Temperature- or pressure-dependent pre-exponential factors, and distributed
  activation-energy tooling.
- In-depth emission and scattering: in-depth radiation is the Beer–Lambert
  absorption of the incident flux (see Technical Reference §8).

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
| Pressure-driven gas transport | Set `κ > 0` on a solid/liquid component | TR §9 |
| Volume change, swelling, intumescence | The swelling factor `γ` on components | TR §10 |
| Lateral shrinkage / variable cross-section | `lateral_shrinkage_law` on the `Material` (with `use_ale = true`) | TR §10 |
| Mesh moving with the solid (surface recession) | `solve(...; use_ale = true)` | TR §11 |
| Planar, cylindrical, or spherical geometry | `coordinates` keyword of the mesh constructors | UG 6 |
| Burnout by cell deactivation | `solve(...; deactivation = true)` | UG 10, TR §14 |
| Temperature or enthalpy as the energy state | `solve(...; energy_state = :enthalpy)` | UG 8, TR §3 |
| Integrator and linear solver | `solve(...; integrator_policy = …, linear_solver = …)` | UG 11, TR §15 |
| Mass and energy conservation audit | `callback = ConservationRecorder()`, `conservation_audit` | UG 9, TR §16 |
| Sensitivity analysis (local finite differences, global Sobol) | finite differences of `solve`; `parameter_sweep` | UG 13 |
| Parameter sweeps (serial or `Distributed`) | `parameter_sweep` | UG 12 |

A few features carry caveats; this guide flags them plainly wherever they
appear:

- The radiation models you can pass to `solve` are `NO_RADIATION`,
  `SURFACE_ABSORPTION`, and `BEER_LAMBERT`. `BEER_LAMBERT` needs planar
  coordinates (see Technical Reference §8).
- **Burnout** (`deactivation = true`) removes cells whose condensed mass falls
  below a fraction of their initial condensed mass, in events located inside
  the step. The state keeps its size; an event is a discrete change of the
  state, so a derivative of an output with respect to a parameter is taken by
  finite differences of whole solves (UG 10, UG 13, Technical Reference §14 and
  §17). Its verified envelope is stated in UG 10.
- Sensitivities are finite differences of whole solves; the package takes no
  derivative through the solve (UG 13).

---

## 1.3 Coordinate and unit conventions

Pyrolysis.jl uses SI units throughout, and the geometric conventions below are
binding across both guides (see Technical Reference §1 and §2 for the full
nomenclature).

- The spatial coordinate is **`z`** (axial, through-thickness), in metres.
- **`z = 0` is the bottom / substrate** (boundary tag `:bottom`).
- **`z = L` is the top / exposed surface** (boundary tag `:top`).
- Heat enters from the **top**; charring materials shrink **downward**.
- Cell 1 is at `z ≈ 0` (bottom); cell `n` is at `z ≈ L` (top).
- Positive flux means transport in the **+z** direction (bottom → top).

Key state and input units:

| Quantity | Symbol | Units |
|---|---|---|
| Temperature | `T` | K |
| Mass concentration (input and output composition) | `ξ_j` | kg·m⁻³ |
| Component mass per cell (the integrator's state) | `m_ij` | kg |
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
> **`h < 0` = exothermic**. Note this is the *opposite* of the ThermaKin
> publishing convention (`h > 0` = exothermic); Gpyro and FDS use this one. The volumetric heat source is
> `Q_rxn = −Σ h_r r_r`, so an endothermic reaction (`h > 0`) produces
> `Q_rxn < 0`. See Technical Reference §6 (overload note H1 in §2).

The composition you specify and read is the **mass concentration `ξ_j`**
(kg/m³), not the mass fraction: `ξ_j = Y_j ρ`. A component initially present at
its bulk density `ρ_j` has `ξ_j(0) = ρ_j`. This is why the default initial
condition places the first solid component at its reference density (see §1.8).
At the start of a solve the concentrations are evaluated at the cell centres
`z_i` and turned into the component masses `m_ij = ξ_j(z_i) V_i` the integrator
carries; the solution reports both.

---

## 1.4 How this guide relates to the Technical Reference

Pyrolysis.jl ships two companion documents that share one unified nomenclature
(Technical Reference §2):

- The **Technical Reference Guide** derives the continuum equations, the
  constitutive laws, the discretization, the Jacobian, and the conservation
  ledgers. Read it when you need to know *why* a result is what it is, or what
  approximation a feature rests on.
- The **User Guide** (this document) teaches the *workflow*: how to define
  materials, reactions, and boundary conditions; how to build and solve a problem;
  how to post-process; and how to troubleshoot.

The two are cross-referenced. Throughout the User Guide you will see "see
Technical Reference §N"; the Technical Reference chapters are numbered:
1 Introduction & Scope; 2 Nomenclature; 3 Governing Equations; 4 Materials &
Property Functions; 5 Effective Properties & Mixing Rules; 6 Reaction Kinetics;
7 Heat Conduction, Convection & Energy Assembly; 8 Thermal Radiation; 9 Gas-Phase
Mass Transport & Pressure; 10 Volume Change & Lateral Shrinkage; 11 The ALE
Framework; 12 Boundary & Initial Conditions; 13 Finite-Volume Discretization;
14 Cell Deactivation and Burnout; 15 Temporal Integration, the DAE, and the
Jacobian; 16 Conservation Ledgers and the Audit; 17 Sensitivity Analysis;
18 Verification & Validation.

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
| UG 10 | Burnout and the correction |
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
extension* mechanism for its plotting recipes, so a current Julia release is
required. Download Julia from <https://julialang.org/downloads/> or install it
via [`juliaup`](https://github.com/JuliaLang/juliaup), the official version
manager:

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

Create or activate a project, then add the package from its Git repository:

```julia
using Pkg
Pkg.activate("my_pyrolysis_study")   # create/activate a local environment
Pkg.add(url = "https://github.com/GraysonBellamy/Pyrolysis.jl")
```

or from a local clone:

```julia
using Pkg
Pkg.activate("my_pyrolysis_study")
Pkg.add(path = "/path/to/Pyrolysis.jl")
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

Pyrolysis.jl keeps its core lightweight by declaring its plotting recipes as a
**package extension**. The extension is wired to a *weak dependency* and
activates automatically the moment you `using` its trigger package alongside
`Pyrolysis` — you never call the extension explicitly. The relevant entries from
the package's `Project.toml` are:

```toml
[weakdeps]
RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[extensions]
PyrolysisPlotsExt = "RecipesBase"
```

| Extension | Trigger package(s) | What it provides |
|---|---|---|
| `PyrolysisPlotsExt` | `Plots` (or `RecipesBase`) | Plot recipes for `PyrolysisSolution` (`plot(sol; kind = …)`) |

It is not required for a basic solve. Install the trigger package into the same
environment, then load it:

```julia
using Pkg
Pkg.add("Plots")

using Pyrolysis
using Plots                      # activates PyrolysisPlotsExt → plot(sol; kind = …)
```

Without `Plots` loaded there is no `plot` function to call; the recipes are
covered in UG 9.

The alternative linear solver needs no extension: the default `StructuredKLU()`
is exported, and LinearSolve's `UMFPACKFactorization()` is passed as
`linear_solver` after `Pkg.add("LinearSolve")` and `using LinearSolve` (UG 8,
UG 11). Sensitivity analysis needs none either: it takes finite differences of
ordinary `solve` calls (UG 13).

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

If you `dev`-ed the package and want to run its test suite:

```julia
using Pkg
Pkg.test("Pyrolysis")
```

The whole suite takes a while. The environment variable `PYROLYSIS_TEST_GROUP`
runs one part of it (`jacobian`, `burnout`, `other`, `domains`, or `all`, the
default), and `PYROLYSIS_TEST_FAST=1` skips the verification cases' convergence
spot checks (the Development page of the documentation).

For a real end-to-end check, run the minimal example in §1.8 — completing it
with `retcode == :Success` is the strongest confirmation that conduction,
kinetics, gas transport, the boundary solve, and the integrator are all working
together.

---

## 1.8 A minimal gas-producing slab

This illustrative porous matrix decomposes into a gas and exchanges with
an explicit vacuum at the exposed face. Its physical porosity is 0.3 and
its lumped superficial coefficient is declared on the condensed component.
The finite endpoint covers complete loss of matrix support.

```julia
using Pyrolysis
material = Material(name=:porous_polymer,
    components=(
        SolidComponent(:virgin; ρ=1190., φ=.3, c=1420., k=.21, λ=1e-6),
        GaseousComponent(:MMA; M=.100, c=1500., k=.02)),
    reactions=(Reaction(:decomposition, :virgin=>:MMA;
                        A=8.5e12, E=188e3, h=870e3),),
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(
        λ=1e-6, κ=0., max_skeletal_occupancy=.001)))
bcs = BoundaryConditionSet(
    thermal=(top=HeatFluxBC(50e3), bottom=AdiabaticBC()),
    mass=(top=AmbientGasBC(GasReservoir(P=0., T=300.)), bottom=ImpermeableBC()))
problem = PyrolysisProblem(mesh=create_uniform_mesh(.01,40), material=material,
    bc_set=bcs, T_initial=z->300., tspan=(0.,60.))
solution = solve(problem; pressure_tolerance=1.)
@assert completed(solution)
println(solution.surface_T[end])
println(maximum(get_pressure(solution)))
```

The default initial concentration fills the first solid at its bulk density;
other species start empty. For a positive-pressure reservoir, declare its
composition and initialize pore gas consistently when equilibrium is intended.
`pressure_tolerance` is an absolute local pressure-error budget in Pa. The
solver also controls the mass, thermal and geometry rows and checks physical
pore-volume admissibility throughout the evolution.

Gas properties are specified by `GaseousComponent(;M,c,k,...)`; transport data
belong to the matrix and its declared `LumpedTransfer` or `DustyGas` model.
The default mass boundary is impermeable and the default thermal boundary is
adiabatic. A finite external mass film is an optional `MassTransferFilm(g)`
with conductance in kg/(m² s). See Chapters 3, 5 and 7 for the complete inputs.

---

## 1.9 Public API surface and the `Internal` escape hatch

`using Pyrolysis` brings a deliberately curated set of names into scope — the
stable, user-facing API. The top-level exports cover:

- **Materials:** `Material`, `SolidComponent`, `GaseousComponent`,
  `LiquidComponent`, `Reaction`, `DepletionLimiter`, the mixing rules
  (`MixingRule`, `PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN`, `CARMAN_KOZENY`),
  `MixingSpec`, `MaterialMixing`, and the phase predicates (`is_solid`, `is_gas`,
  `is_liquid`, `phase`, `molar_mass`).
- **Property functions:** `ConstantProperty`, `LinearProperty`,
  `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`,
  `SutherlandProperty`, `StateDependentProperty`.
- **Boundary conditions:** the thermal conditions (`HeatFluxBC`, `ConvectiveBC`,
  `RadiativeBC`, `RadiativeFluxBC`, `CombinedThermalBC`, `DirichletThermalBC`,
  `AdiabaticBC`), the mass conditions (`ImpermeableBC`, `MassFluxBC`,
  `AmbientGasBC`, `GasReservoir`, `MassTransferFilm`), the pressure
  conditions (`NeumannPressureBC`, `ConvectivePressureBC`),
  the container `BoundaryConditionSet`; the blowing models of `ConvectiveBC`
  (`AbstractBlowingModel`, `NoBlowing`, `CouetteFilm`, `PropertyRatioFilm`,
  `PoolingFilm`, `FunctionBlowing`) with the gas data and helpers they read
  (`GasSpecies`, `NITROGEN`, `nasa7_heat_capacity`, `eucken_viscosity`,
  `fuller_diffusivity`, `BlowingFace`, `blowing_state`, `blowing_ratio`,
  `blowing_factor`); and `background_gas_concentration` for pore-air
  initialization in lumped gas transport problems.
- **Radiation models:** `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`,
  `BEER_LAMBERT`.
- **The mesh:** `create_uniform_mesh`, `create_stretched_mesh`,
  `create_mesh_from_nodes`, the coordinate systems `Planar`, `Cylindrical`,
  `Spherical`, and the mesh's parts `MeshTopology`, `MeshGeometry`.
- **Problem and solve:** `PyrolysisProblem`, `PyrolysisSolution`,
  `SolverStatistics`, `solve`, `setup`, `integrate`, `CoreProblem` (the set-up
  problem `setup` returns and `integrate` runs, for custom loops, UG 11), and
  `parameter_sweep`.
- **The mesh correction and the integrator:** `NoCorrection`,
  `EquidistributionCorrection`, `ActivationTime` (UG 10); the integrator
  policies `default_policy`, `rodas5p_policy`, `fbdf_policy`,
  `kencarp4_policy`, and the linear solver `StructuredKLU` (UG 11).
- **The solution's accessors and the conservation audit:** `get_times`,
  `get_temperatures`, `get_concentrations`, `get_gas_generation_rate`,
  `n_timesteps`, `succeeded`, `completed`, `get_profile`, `get_cell_positions`,
  `interpolate_solution`, `ConservationRecorder`, `conservation_audit` (UG 9).

Anything not exported at the top level remains reachable in one of two ways:

```julia
using Pyrolysis.Internal     # re-exports every submodule's public surface in one namespace
# e.g. residual!, dense_jacobian, initial_fields, effective_conductivity, … all become available

# Or reach a specific submodule directly:
Pyrolysis.Properties.effective_conductivity
Pyrolysis.Problem.initial_fields
```

The submodules, in load order, are: `Materials`, `Properties`,
`BoundaryConditions`, `Physics` (the shared physics layers), `Meshes`, `States`,
`Kernels`, `Closures`, `Assembly`, `Deactivation`, `ColoredJacobian` (the
discretization), `Problem`, `Outputs`, and `Integration`. `Internal` re-exports
the union of their exports: among them the residual (`residual!`), the colored
Jacobian (`jacobian!`, and `dense_jacobian`, the dense forward-mode reference the
tests compare it with), the state layout and its indices, the kernels and
closures of a set-up problem, the deactivation events, and the
multi-dimensional mesh builders (`box_mesh`, `cell_mesh`) at the residual level.
A name that two submodules export as different objects is not re-exported;
reach it through its submodule.

A couple of practical notes you will encounter later in the guide:

- The simulation modes — a fixed or moving mesh, the radiation rows, lumped or
  lumped gas transport gas transport, the energy state, the variable cross-section — are
  not types you pick. `setup` decides them from the material (a permeable
  component selects gas transport, a `lateral_shrinkage_law` the
  cross-section rows) and from the `solve` keywords (`use_ale`,
  `radiation_model`, `energy_state`, `deactivation`) (UG 7).
- Symbols reached through `Internal` or a submodule are **not** part of the stable
  public API and may change between minor releases. Prefer the top-level exports
  in production code; treat `Internal` as an escape hatch.

> **Caveat — API stability.** Top-level `using Pyrolysis` exports are stable.
> Submodule / `Internal` symbols can move between minor releases.

---

## 1.10 Where the examples live

The repository ships a set of runnable example scripts under the `examples/`
directory, with their own project environment (it adds `Plots` and other
post-processing packages). These are the best starting points for real
configurations, and UG 14 walks through them in detail. The most useful to start
with:

| File | What it demonstrates |
|---|---|
| `examples/06_cone_calorimeter.jl` | PMMA, 5.8 mm, at 25 kW/m²: single-step kinetics, temperature-dependent properties, combined thermal conditions, a mesh moving with the solid and cell deactivation; overlays ThermaKin outputs when present |
| `examples/06_cone_calorimeter_surface_absorption.jl` | The same case with `radiation_model = SURFACE_ABSORPTION`; writes figures and a CSV to `output/` |
| `examples/07_douglas_fir_cone.jl` | Douglas fir, 12 mm, at 75 kW/m²: five-step kinetics with moisture evaporation, gas transport with background pore air, surface absorption, a moving mesh |
| `examples/08_gas_transport.jl` | Lumped and dusty-gas transport with a finite reservoir film |
| `examples/08_pressure_treated_wood_cone.jl` | Pressure-treated wood with moisture and multi-step Arrhenius kinetics |
| `examples/convergence_douglas_fir.jl` | Grid-convergence study (uses `Distributed`) |
| `examples/moisture_hos_comparison.jl` | Moisture transport with state-dependent heat of sorption |
| `examples/export_results_csv.jl` | Writes a solution's series and profiles to CSV files |

The subdirectories hold larger studies: `inverse_analysis/` (fitting white-pine
properties to cone data), `sensitivity_analysis/`, `literature_comparison/`, and
`macfp/` (the MaCFP gasification targets).

To run an example, start Julia in the examples environment from the repository
root:

```bash
julia --project=examples examples/06_cone_calorimeter.jl
```

or, from a running session:

```julia
using Pkg
Pkg.activate("examples")     # the examples directory has its own Project.toml
Pkg.instantiate()            # first time only: resolve the example dependencies
include("examples/06_cone_calorimeter.jl")
```

Every simulation script runs a short check of itself with the environment
variable `PYROLYSIS_SMOKE=1` (a 10–20 s horizon). The **Examples** page of the
documentation site (source `docs/src/examples.md`) is the authoritative
catalogue of what each script covers. With the install verified and the minimal
example running, you are ready for UG 2, which builds a full cone-calorimeter
simulation step by step.
