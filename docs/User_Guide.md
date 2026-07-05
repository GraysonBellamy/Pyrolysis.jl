# Pyrolysis.jl — User Guide

A practical guide to using Pyrolysis.jl: installation, materials, reactions, boundary conditions, meshes, building and running simulations, post-processing, adaptivity, solver/Jacobian configuration, parameter sweeps, sensitivity analysis, worked examples, troubleshooting, and a complete API reference. Companion to the **Technical Reference Guide**.

_This document is the canonical source; the per-chapter views in `docs/guide/user_guide` are regenerated from it with `julia scripts/split_guides.jl`. Each chapter is grounded in the Pyrolysis.jl source and was adversarially verified against the implementation._

---

## Contents

- [1. Introduction and Installation](#1-introduction-and-installation)
- [2. Quick-Start Tutorial: Your First Cone-Calorimeter Run](#2-quick-start-tutorial-your-first-cone-calorimeter-run)
- [3. Defining Materials, Components, and Properties](#3-defining-materials-components-and-properties)
- [4. Defining Reactions](#4-defining-reactions)
- [5. Boundary and Initial Conditions](#5-boundary-and-initial-conditions)
- [6. Meshes and Geometry](#6-meshes-and-geometry)
- [7. Building a Problem: Simulation Modes](#7-building-a-problem-simulation-modes)
- [8. Running the Solver](#8-running-the-solver)
- [9. The Solution Object and Plotting](#9-the-solution-object-and-plotting)
- [10. Adaptive Meshing Options](#10-adaptive-meshing-options)
- [11. Solver and Jacobian Configuration; Extensions](#11-solver-and-jacobian-configuration-extensions)
- [12. Parameter Sweeps](#12-parameter-sweeps)
- [13. Sensitivity Analysis](#13-sensitivity-analysis)
- [14. Worked Examples and Tutorials](#14-worked-examples-and-tutorials)
- [15. Troubleshooting and Performance](#15-troubleshooting-and-performance)
- [16. API and Parameter Reference](#16-api-and-parameter-reference)

---

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
equation uses **Formulation A** (the non-conservative convention in which gas
sensible storage is removed from the matrix heat capacity and carried by the
advection source `S_conv`), which matches the FDS/Gpyro/ThermaKin convention (see
Technical Reference §3 and §7).

### What is not modeled

To set expectations up front, the current model does **not** support:

- Multi-dimensional (2D/3D) transport — the model is strictly 1D
  through-thickness.
- A full two-temperature gas phase (the gas is treated quasi-steady; the cost of
  this approximation is typically well under 1% at atmospheric pressure — see
  Technical Reference §3 and §16).
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
| Surface-cell depletion / merging | `solve(...; use_ale = true, handle_depletion = true)` | TR §14 |
| Adaptive mesh refinement (r-refinement, indicators) | `Pyrolysis.Internal` / submodule access | TR §14 |
| Conservation diagnostics | `solve(...; diagnostics = true)` | TR §16 |
| Sensitivity analysis (forward / adjoint) | `forward_sensitivity`, `adjoint_sensitivity` (extension) | TR §17 |
| Parameter sweeps (parallel) | `parameter_sweep` | UG §12 |

A few features are explicitly **experimental** or carry caveats; this guide flags
them plainly wherever they appear:

- **P1 quasi-steady radiation** lives in `Pyrolysis.Experimental`, not in the core
  solve path. The exported radiation models you can pass to `solve` are
  `NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT`; requesting
  `P1_QUASI_STEADY` is rejected by `solve` (see Technical Reference §8).
- The **surface transpiration boundary condition** is incompatible with the
  default `S_conv` energy convention. The volumetric advection source already
  carries the gas enthalpy, so enabling transpiration would double-count it.
  `use_transpiration_bc = true` is rejected at construction; leave it at its
  default `false` (see Technical Reference §7 and §12).
- **Surface-cell depletion handling** (`handle_depletion = true`) resizes the
  state vector mid-solve and is therefore incompatible with a pluggable Jacobian
  backend. When you enable it, the solver falls back to the integrator's own
  colored finite-difference Jacobian (see §1.7 and Technical Reference §14, §15).
- **h-AMR** (cell splitting / coarsening) is quarantined in
  `Pyrolysis.Experimental`; the supported adaptivity in the core package is
  r-refinement (node relocation) and surface-depletion merging.

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

Pyrolysis.jl keeps its core lightweight by declaring several heavy or
vendor-specific capabilities as **package extensions**. Each extension is wired to
a *weak dependency* and activates automatically the moment you `using` the
corresponding package alongside `Pyrolysis` — you never call the extension
explicitly, you just load its trigger package. The relevant entries from the
package's `Project.toml` are:

```toml
[extensions]
PyrolysisAMGExt = "AlgebraicMultigrid"
PyrolysisKLUExt = "KLU"
PyrolysisPardisoExt = "Pardiso"
PyrolysisPlotsExt = "RecipesBase"
PyrolysisSciMLSensitivityExt = ["SciMLSensitivity", "Zygote"]
PyrolysisSymbolicsExt = "Symbolics"
```

The table below summarizes each extension, the package(s) that trigger it, and
what it adds. None of these is required for a basic solve.

| Extension | Trigger package(s) | What it provides |
|---|---|---|
| `PyrolysisKLUExt` | `KLU` | `Pyrolysis.klu_factorization()` — KLU sparse linear solver |
| `PyrolysisPardisoExt` | `Pardiso` | `Pyrolysis.pardiso_factorization()` — Intel MKL Pardiso solver |
| `PyrolysisAMGExt` | `AlgebraicMultigrid` | `Pyrolysis.amg_preconditioner(A)` — AMG preconditioner for Krylov solvers |
| `PyrolysisPlotsExt` | `Plots` (or `RecipesBase`) | Plot recipes for `PyrolysisSolution` (`plot(sol; kind = …)`) |
| `PyrolysisSciMLSensitivityExt` | `SciMLSensitivity`, `Zygote` | `forward_sensitivity`, `adjoint_sensitivity` |
| `PyrolysisSymbolicsExt` | `Symbolics` | `verify_mixing_derivatives` (test/verification only) |

Install the trigger packages you need into the same environment, then load them:

```julia
using Pkg
Pkg.add(["Plots", "KLU", "SciMLSensitivity", "Zygote"])

using Pyrolysis
using Plots           # activates PyrolysisPlotsExt → plot(sol; kind = …)
using KLU             # activates PyrolysisKLUExt   → klu_factorization()
using SciMLSensitivity, Zygote   # activates PyrolysisSciMLSensitivityExt
```

### Stub behavior when an extension is not loaded

The functions `klu_factorization`, `pardiso_factorization`, `amg_preconditioner`,
`forward_sensitivity`, `adjoint_sensitivity`, and `verify_mixing_derivatives` are
exported by `Pyrolysis` as *stubs* with no methods. Calling one before its
extension is loaded raises a `MethodError` — that is the signal that you need to
`using` the corresponding trigger package. For example, calling
`forward_sensitivity(...)` without `SciMLSensitivity` and `Zygote` loaded will
error; add them and the concrete method appears.

The extensions are independent peers — no extension depends on another — so you
load only the ones a given study requires. Plotting recipes, sensitivity, and
alternate linear solvers are each covered in detail in their own chapters
(UG §9, UG §11, UG §13).

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

> **Depletion and the Jacobian.** If you later enable surface-cell depletion
> handling (`solve(...; use_ale = true, handle_depletion = true)`), be aware that
> it resizes the state vector while integrating, which a pluggable `JacobianSpec`
> cannot survive. The solver therefore refuses to combine `handle_depletion = true`
> with an explicit `jacobian` argument and, when depletion is on, falls back to the
> integrator's own colored finite-difference Jacobian. See UG §10 and Technical
> Reference §14–§15.

---

## 1.8 A minimal "hello pyrolysis" example

The script below is a complete, runnable simulation: a 1 cm slab of PMMA-like
material, modeled as one solid that decomposes into one gas, heated from the top
by a constant 50 kW/m² flux, insulated at the bottom, integrated for 60 seconds.
Every function, type, and keyword used here is part of the public API and runs
against the real package (this exact pattern was verified end-to-end while writing
this chapter; it returns `retcode = :Success`).

```julia
using Pyrolysis

# 1. Define the material: one solid component + one gaseous product.
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86),
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
  capacity `c`, and conductivity `k`; `ε` (emissivity) defaults to `0.9`, the
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
  `1.0`, `T_min` to `0.0`, `T_max` to `Inf`. Remember `h > 0` is endothermic.
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
  `radiation_model = NO_RADIATION`, `energy_form = :standard`,
  `diagnostics = false`, `validate = true`, `progress = true`,
  `handle_depletion = false`, and `use_transpiration_bc = false`. Setting
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
plot(solution; kind = :residual)     # energy-conservation residual
```

---

## 1.9 Public API surface and the `Internal` escape hatch

`using Pyrolysis` brings a deliberately curated set of names into scope — the
stable, user-facing API. Broadly, the top-level exports cover:

- **Materials:** `Material`, `SolidComponent`, `GaseousComponent`,
  `LiquidComponent`, `Reaction`, `ReactionSet`, the mixing rules
  (`PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN`, `CARMAN_KOZENY`), `MixingSpec`,
  `MaterialMixing`, and the phase predicates (`is_solid`, `is_gas`, `is_liquid`,
  `phase`, `molar_mass`).
- **Property functions:** `ConstantProperty`, `LinearProperty`,
  `PolynomialProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`,
  `StateDependentProperty`.
- **Boundary conditions:** the thermal BCs (`HeatFluxBC`, `ConvectiveBC`,
  `RadiativeBC`, `RadiativeFluxBC`, `CombinedThermalBC`, `DirichletThermalBC`,
  `AdiabaticBC`), the mass BCs (`ImpermeableBC`, `MassFluxBC`, `ConvectiveMassBC`,
  `DiffusiveMassBC`), the pressure BCs (`DirichletPressureBC`, `NeumannPressureBC`,
  `ConvectivePressureBC`), and the container `BoundaryConditionSet`.
- **Radiation models:** `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`,
  `BEER_LAMBERT` (the `P1_QUASI_STEADY` value is exported from
  `Pyrolysis.Experimental`).
- **Geometry:** `create_uniform_mesh`, `create_stretched_mesh`.
- **Problem and solve:** `PyrolysisProblem`, `PyrolysisSolution`, `solve`,
  `parameter_sweep`.
- **Jacobian configuration:** `JacobianSpec`, `Structured`, the linear-solve
  strategies (`DirectSolve`, `StructuredSolve`, `ApproxSolve`, `MatrixFreeSolve`),
  `DenseAD_ForwardDiff` (a verification reference), `verify_jacobian_at`, and more
  (UG §11).
- **Extension stubs:** `klu_factorization`, `pardiso_factorization`,
  `amg_preconditioner`, `forward_sensitivity`, `adjoint_sensitivity`,
  `verify_mixing_derivatives`.

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
`Diagnostics`, and `Solver`. The `Experimental` namespace (P1 radiation, h-AMR) is
*not* part of `Internal` and is accessed only as `Pyrolysis.Experimental.X`.

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
> `Pyrolysis.Experimental` is explicitly unstable and will be restructured or
> migrated to dedicated extensions in later phases.

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

The `examples/README.md` is the authoritative catalogue of what each script
covers. With the install verified and the minimal example running, you are ready
for UG §2, which builds a full cone-calorimeter simulation step by step.


---

# 2. Quick-Start Tutorial: Your First Cone-Calorimeter Run

This chapter walks you through a complete, runnable cone-calorimeter simulation
from an empty Julia session to plotted results. We model a slab of PMMA
(poly(methyl methacrylate)) heated from above by a radiant cone, decomposing
into a monomer gas plus a small char residue. By the end you will have built a
material, defined boundary conditions, created a mesh, assembled a problem,
solved it, and extracted surface temperature and mass-loss rate. Every code
snippet uses the real public API and can be pasted into the REPL as-is. The
deeper "why" of each physics choice is deferred to the Technical Reference; the
goal here is to get you to a working run and a mental model you can reuse.

The mental model is a five-step pipeline that this chapter follows in order:

> **material → boundary conditions → mesh → problem → solve → solution**

You define *what* decomposes (the `Material`), *how* energy and mass cross the
edges (the `BoundaryConditionSet`), *where* the domain lives (the `Mesh1D`),
bundle them into a `PyrolysisProblem`, hand that to `solve`, and read results
off the returned `PyrolysisSolution`.

---

## 2.1 Before you begin: coordinates, units, and conventions

Two conventions are worth fixing in your head before writing any code, because
they govern how every boundary condition and result is interpreted.

**Coordinate axis.** The single spatial coordinate is `z`, running through the
thickness of the sample:

- `z = 0` is the **bottom / substrate** (boundary tag `:bottom`, internal id 2).
- `z = L` is the **top / exposed surface** (boundary tag `:top`, internal id 1).

Heat enters from the top; the material recedes toward the bottom as it
decomposes. Positive flux is in the `+z` direction (bottom → top). See
Technical Reference §13 for the finite-volume layout and §12 for the boundary
tag/id bridge.

**Units.** Pyrolysis.jl is strict SI throughout — metres, kilograms, seconds,
kelvin, joules, watts, pascals. There are no implicit "kW" or "mm". A 25 kW/m²
cone flux is written `25e3`; a 5.8 mm sample is `5.8e-3`. The primary species
state variable is the **mass concentration** `ξⱼ` in kg/m³ (not a mass
fraction). For the full symbol table see Technical Reference §2.

**One import.** Everything in this tutorial comes from a single package, plus
`Plots` for the figures:

```julia
using Pyrolysis
using Plots          # optional; enables plot(sol) recipes
```

`using Pyrolysis` brings in the curated public surface: the material and
component constructors, the boundary-condition types, the mesh creators,
`PyrolysisProblem`, `PyrolysisSolution`, and `solve`. Loading `Plots`
automatically activates the plotting recipes described in §2.9.

---

## 2.2 Step 1 — Define the material

A `Material` is a tuple of **components** (the chemical species, each a phase
with its own thermophysical properties) plus a tuple of **reactions** (the
Arrhenius kinetics that convert one species into others), plus the
**mixing rules** used to blend per-component properties into effective
mixture properties. PMMA needs three components:

1. **Virgin PMMA** — the solid that decomposes.
2. **Char** — a small solid residue left behind.
3. **MMA gas** — the monomer that escapes through the matrix.

### 2.2.1 Temperature-dependent properties

PMMA's heat capacity and conductivity vary with temperature. The most general
way to encode an arbitrary `f(T)` is a `CallableProperty`, which wraps any
single-argument function:

```julia
# Heat capacity [J/(kg·K)], piecewise about 378 K
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2 * T) : (1040.0 + 3.08 * T))

# Thermal conductivity [W/(m·K)], piecewise about 378 K
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * T))
```

`CallableProperty` is one of several property-function types. You can also pass
a bare number (auto-converted to a `ConstantProperty`), a `LinearProperty(a, b)`
for `a + b·T`, a `PolynomialProperty`, an `ArrheniusProperty`, or a
`TableProperty` for tabulated data. Property types and their derivatives are
covered in full in User Guide §UG-3 and Technical Reference §4.

### 2.2.2 Component constructors

Each component is built with a keyword constructor. Solids and liquids carry a
bulk density, heat capacity, conductivity, emissivity, an in-depth radiation
absorption coefficient, and optionally permeability and a swelling factor.
Gases instead carry a molar mass (for the ideal-gas density) and a gas-transfer
(diffusion) coefficient `λ`.

```julia
virgin = SolidComponent(
    :virgin,
    ρ = 1210.0,     # bulk density [kg/m³]
    c = pmma_cp,    # heat capacity [J/(kg·K)]
    k = pmma_k,     # thermal conductivity [W/(m·K)]
    ε = 0.96,       # surface emissivity [-]   (default 0.9)
    α = 1.47,       # mass absorption coefficient [m²/kg]  (default 0.0)
)

char = SolidComponent(
    :char,
    ρ = 1210.0,
    c = pmma_cp,
    k = pmma_k,
    ε = 0.96,
    α = 1.47,
)

mma = GaseousComponent(
    :MMA,
    M = 0.10012,    # molar mass [kg/mol]  (required)
    c = 1500.0,     # gas heat capacity [J/(kg·K)]
    k = 0.02,       # gas conductivity [W/(m·K)]
    λ = 1e-5,       # gas-transfer (diffusion) coefficient [m²/s]
)
```

The first positional argument is the component **name** (a `Symbol`). All
thermophysical values are keywords. Defaults worth knowing:

| Keyword | Meaning | Default (`SolidComponent`) | Default (`GaseousComponent`) |
|---|---|---|---|
| `ρ` | bulk density [kg/m³] | required | ideal-gas law from `M` |
| `c` | heat capacity [J/(kg·K)] | required | required |
| `k` | thermal conductivity [W/(m·K)] | required | required |
| `ε` | surface emissivity [-] | `0.9` | `0.0` |
| `α` | mass absorption coeff. [m²/kg] | `0.0` | `0.0` |
| `γ` | swelling factor [-] | `1.0` (solid shrinks) | `0.0` (gas escapes) |
| `λ` | gas-transfer coeff. [m²/s] | `0.0` | required |
| `κ` | permeability [m²] | `0.0` (impermeable) | n/a (gases flow through) |
| `M` | molar mass [kg/mol] | n/a | required |

A gas's density defaults to the ideal-gas law `ρ_g = P_ref·M/(R_g·T)` with
`P_ref = 101325` Pa and `R_g = 8.314...` J/(mol·K). The bulk-vs-intrinsic
density distinction (`ρ_intrinsic`, used only in the porosity formula) and the
swelling factor `γ` are explained in Technical Reference §4 and §10; the
defaults above are correct for this opaque, non-swelling case.

### 2.2.3 The decomposition reaction

PMMA depolymerises in a single endothermic step: virgin PMMA becomes a tiny
char fraction plus MMA gas. Reactions reference components by **index** into the
components tuple (1 = virgin, 2 = char, 3 = MMA). The `1 => (2, 3)` syntax means
"reactant 1 produces products 2 and 3", and `yields` gives the per-product mass
fractions, which **must sum to 1.0**:

```julia
decomp = Reaction(
    :decomposition,
    1 => (2, 3),         # virgin → char + MMA
    A = 1.5e14,          # pre-exponential factor [s⁻¹]
    E = 2.03e5,          # activation energy [J/mol]
    h = 820e3,           # heat of reaction [J/kg]; h > 0 = endothermic
    n = 1.0,             # reaction order (default 1.0)
    yields = (0.002, 0.998),  # 0.2 % char, 99.8 % gas
)
```

Two conventions to internalise:

- **Heat-of-reaction sign:** `h > 0` is **endothermic** (the reaction absorbs
  heat and cools the material); `h < 0` is **exothermic**. PMMA pyrolysis is
  endothermic, so `h = 820e3` is positive. (Note: ThermaKin and Gpyro publish
  the opposite sign convention — see Technical Reference §6 and the §2
  nomenclature note H1.)
- **Single-reactant only:** each reaction consumes exactly one reactant species
  (with mass stoichiometry 1.0). Multi-reactant reactions are not supported.
  Multi-step schemes are modelled as a *tuple of single-reactant reactions* (see
  User Guide §UG-4).

For a one-product reaction the simpler form `Reaction(:name, 1 => 2; A, E, h)`
omits `yields` (the single product then has stoichiometry 1.0). Temperature
gates `T_min`/`T_max` default to `0.0` and `Inf` (always active). See Technical
Reference §6 for the full rate law, the smooth `tanh` temperature gates, and the
depletion limiter.

### 2.2.4 Assemble the `Material`

Finally bundle components, reactions, and mixing rules. `MaterialMixing` selects
how per-component thermal conductivity (`k`), gas-transfer (`λ`), and
permeability (`κ`) are averaged into effective mixture properties:

```julia
pmma = Material(
    name = :PMMA,
    components = (virgin, char, mma),
    reactions  = (decomp,),
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),  # k = 0.5·k_parallel + 0.5·k_series
        λ = MixingSpec(PARALLEL),       # volume-weighted gas transfer
    ),
)
```

`MaterialMixing` defaults are `k = MixingSpec(WEIGHTED, 0.5)`,
`λ = MixingSpec(PARALLEL, 0.5)`, and `κ = MixingSpec(WEIGHTED, 0.5)`, so you
only need to pass the channels you want to change. The weight `β` (the second
argument of `MixingSpec`) is **only used by the `WEIGHTED` rule**; for `PARALLEL`
(and the other non-`WEIGHTED` rules) it is carried but ignored, so
`MixingSpec(PARALLEL)` and `MixingSpec(PARALLEL, 0.5)` are equivalent — that is
why the example above writes `λ = MixingSpec(PARALLEL)` with no weight. The
available rules are `PARALLEL`, `SERIES`, `WEIGHTED` (with weight `β`),
`BRUGGEMAN` (conductivity only), and `CARMAN_KOZENY` (permeability only). The
mixing theory lives in Technical Reference §5; how to choose rules is in User
Guide §UG-3.

> **Tip — components reference by index, not name (here).** Because we passed
> `virgin`, `char`, `mma` in that order, the reaction's `1 => (2, 3)` lines up
> with virgin → char + MMA. You can instead write reactions with *symbolic*
> names (`Reaction(:decomp, :virgin => (:char, :MMA); ...)`), which `Material`
> resolves to indices for you. Both styles are valid; the numeric style is used
> in the shipped `examples/06_cone_calorimeter.jl`.

---

## 2.3 Step 2 — Define the boundary conditions

The cone calorimeter heats the top face with a radiant flux while that same
face loses heat by re-radiation and convection and lets gas escape. The bottom
face sits on a sample holder and loses a little heat to the surroundings while
remaining gas-tight. Boundary conditions come in two kinds we need here —
**thermal** and **mass** — each specified independently for `:top` and
`:bottom`, then collected into a `BoundaryConditionSet`.

### 2.3.1 Thermal boundary conditions

Thermal BC types **add** with the `+` operator, building a `CombinedThermalBC`
whose flux is the sum of the parts. The top face combines three effects:

```julia
q_external = 25e3   # incident cone flux [W/m²] = 25 kW/m²
h_conv     = 7.2    # convective coefficient [W/(m²·K)]

top_thermal =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +  # incident cone radiation in
    RadiativeBC(ε = 0.96, T_inf = 330.0) +                     # re-radiation loss out
    ConvectiveBC(h = h_conv, T_inf = 330.0)                    # convective loss out
```

What each piece does:

- **`RadiativeFluxBC(flux, absorptivity)`** is the *incident* external radiation
  from the cone heater. The absorbed flux is `flux × absorptivity`. If you pass
  `absorptivity = nothing` (the default), the effective surface emissivity is
  used instead via Kirchhoff's law (`α = ε_eff`), automatically tracking
  composition as the surface chars. The way this flux is deposited depends on
  the `radiation_model` you choose at solve time (see §2.6).
- **`RadiativeBC(ε, T_inf, F)`** is the surface's own **re-radiation loss**,
  `q = F·ε·σ·(T_∞⁴ − T_s⁴)`. With `ε = nothing` it uses the effective
  emissivity; the view factor `F` defaults to `1.0`. Because `T_s` exceeds
  `T_∞` once heated, this term is a loss.

> **Don't confuse the two radiative BCs.** `RadiativeBC` applies the
> *re-radiation loss* (the `T_s⁴` term) using the surface's **own emissivity**
> `ε` — it has **no `absorptivity` field**. `RadiativeFluxBC` applies *incident*
> external radiation scaled by **`absorptivity`** — it has no `T_s⁴` term. The
> example above is correct: `RadiativeFluxBC(... absorptivity = 0.96)` for the
> cone flux in, and `RadiativeBC(ε = 0.96, ...)` for re-radiation out.
- **`ConvectiveBC(h, T_inf)`** is Newton cooling, `q = h·(T_∞ − T_s)` — also a
  loss once the surface is hot.

The bottom face only loses heat:

```julia
bottom_thermal =
    RadiativeBC(ε = 0.95, T_inf = 307.0) +
    ConvectiveBC(h = 4.0, T_inf = 307.0)
```

Other thermal BCs available (full reference in User Guide §UG-5 / Technical
Reference §12):

| Type | Flux | Notes |
|---|---|---|
| `AdiabaticBC()` | `0` | perfect insulation |
| `HeatFluxBC(q)` | prescribed `q` or `q(t)` | `+` into domain; not subject to absorptivity |
| `ConvectiveBC(h, T_inf)` | `h(T_∞ − T_s)` | `h`, `T_inf` may be functions of `t` |
| `RadiativeBC(ε, T_inf, F)` | `F·ε·σ·(T_∞⁴ − T_s⁴)` | re-radiation; `ε = nothing` ⇒ `ε_eff` |
| `RadiativeFluxBC(flux, absorptivity)` | `flux·α` (model-dependent) | incident external radiation |
| `DirichletThermalBC(T)` | enforces `T_s = T` | fixed-temperature face |
| `bc1 + bc2 + …` | sum | `CombinedThermalBC` |

> **Note — surface temperature is solved, not assumed.** In a cell-centred
> finite-volume method the temperature lives at cell *centres*, but the
> radiative loss `σ T_s⁴` is highly nonlinear in the *surface* temperature.
> Pyrolysis.jl therefore runs a small Newton solve at each boundary to find the
> true `T_s` that balances applied and conducted flux. You never call this
> yourself; the solved value lands in `solution.surface_T`. See Technical
> Reference §12.

### 2.3.2 Mass boundary conditions

The top face lets MMA gas escape; the bottom is impermeable:

```julia
top_mass    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)  # gas escapes; h_m in m/s
bottom_mass = ImpermeableBC()                            # no gas crosses
```

- **`ConvectiveMassBC(h_m, ξ_inf)`** drives a per-species flux
  `J = h_m·(ξ_∞ − ξ_s)`. With ambient concentration `ξ_inf = 0.0`, gas leaves
  the surface at a rate set by the mass-transfer coefficient `h_m` [m/s].
- **`ImpermeableBC()`** sets `J = 0` for all species — the default for a face if
  you omit it.

> **Mind the units: `h_m` [m/s] is not `h` [W/(m²·K)].** The mass-transfer
> coefficient `h_m ≈ 0.05 m/s` used here is typical for a cone heater with
> moderate air flow. It is related to the convective heat-transfer coefficient
> `h ≈ 7.2 W/(m²·K)` by the Chilton–Colburn analogy, roughly
> `h_m ≈ h / (ρ·c_p)` (air `ρ·c_p ≈ 1.2 × 1005 ≈ 1.2e3 J/(m³·K)`). Do not swap
> `h_m` and `h`: they have different units and different scaling.

Alternatives include `MassFluxBC(J)` (prescribed flux), `DiffusiveMassBC(h_m)`
(internal-diffusion-limited with optional external resistance; `DiffusiveMassBC()`
means no external resistance), and the pressure BCs used in Darcy mode. See
User Guide §UG-5.

> **Caveat — do not enable surface transpiration here.** Pyrolysis.jl's default
> energy form (`energy_form = :standard`) already carries gas enthalpy through a
> volumetric convection source `S_conv`. The opt-in surface-transpiration term
> (`use_transpiration_bc = true`) would double-count that enthalpy, and `solve`
> rejects it with an error. Leave it off (the default). See §2.6 and Technical
> Reference §7.

### 2.3.3 Collect into a `BoundaryConditionSet`

The set is tag-keyed: a `NamedTuple` for `thermal` and one for `mass`, each with
`top` and `bottom` entries. Any face you omit defaults to `AdiabaticBC()`
(thermal) or `ImpermeableBC()` (mass).

```julia
bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal, bottom = bottom_thermal),
    mass    = (top = top_mass,    bottom = bottom_mass),
)
```

---

## 2.4 Step 3 — Create the mesh

A 1D finite-volume mesh divides the sample thickness into `n_cells` cells.
`create_uniform_mesh` is the workhorse; its signature is
`create_uniform_mesh(L, n_cells, Val(NC); kwargs...)`, where `L` is the domain
thickness in metres, `n_cells` is the cell count, and `Val(NC)` encodes the
number of components (here 3) as a type parameter for compile-time
specialisation.

```julia
sample_thickness = 5.8e-3   # 5.8 mm
n_cells          = 72       # → 0.08 mm per cell
T_ambient        = 300.0    # initial / ambient temperature [K]

mesh = create_uniform_mesh(
    sample_thickness,
    n_cells,
    Val(3);                          # NC = 3 components (virgin, char, MMA)
    T_initial = T_ambient,           # initial cell temperature [K] (default 300.0)
    ξ_initial = (1210.0, 0.0, 0.0),  # initial concentrations [kg/m³]: all virgin
)
```

The `Val(3)` argument **must match** the number of components in your material.
`ξ_initial` is a tuple of one concentration per component: the sample starts as
solid virgin PMMA at 1210 kg/m³, with no char and no gas. `cross_section_area`
defaults to `1.0` m² — in a 1D model the cross-section is a gauge normalisation
(MLR and mass are reported per unit area), so leaving it at 1.0 is correct and
conventional. See User Guide §UG-6 for when to change it.

To refine the mesh near the heated surface (where the reaction front is
steepest), use `create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from = :top)`
instead, with a `stretch` factor above 1.0. That is optional; a uniform mesh is
fine for a first run.

> **Tip — resolution vs. cost.** 72 cells over 5.8 mm resolves the thermal
> wave and reaction zone comfortably for PMMA. Too few cells smears the surface
> gradient; too many slows the solve (the diffusion time-step scales like
> `Δz²`). A mesh-convergence study is described in Technical Reference §18.

---

## 2.5 Step 4 — Build the problem

`PyrolysisProblem` ties the mesh, material, boundary conditions, initial
conditions, and time span together. Initial conditions are given as
*functions of `z`* (allowing spatially varying profiles), here constant:

```julia
t_final = 1600.0   # simulate 1600 s

problem = PyrolysisProblem(
    mesh      = mesh,
    material  = pmma,
    bc_set    = bc_set,
    T_initial = z -> T_ambient,                         # uniform 300 K
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],      # one closure per component
    tspan     = (0.0, t_final),
    dt_initial = 1e-3,                                  # initial step guess [s]
    output_times = collect(range(0.0, t_final, length = 601)),  # save ~601 points
)
```

Argument notes:

- **`T_initial`** — a closure `z -> T`. Defaults to `z -> 300.0` if omitted.
- **`ξ_initial`** — a `Vector` of closures, one per component (`z -> ξⱼ`). If
  omitted, the first solid component defaults to its reference density and the
  rest to zero. A single closure is broadcast to all components.
- **`tspan`** — the `(t_start, t_end)` time window in seconds (required).
- **`dt_initial`** — the integrator's first-step guess (default `1e-4`); the
  adaptive solver adjusts from there.
- **`output_times`** — the times at which the solution is sampled and stored.
  Fewer points means less interpolation overhead and memory; if left empty, a
  reasonable default grid is generated.

The mesh's `ξ_initial` (Step 3) and the problem's `ξ_initial` (here) describe
the same starting state; the problem-level closures are authoritative for the
solve. Keeping them consistent avoids confusion.

`PyrolysisProblem` does not yet choose the *physics modes* (ALE vs. fixed mesh,
radiation model, etc.) — those are selected as options to `solve`, described
next. The simulation-mode trait system is covered in User Guide §UG-7.

---

## 2.6 Step 5 — Solve

`solve(problem; kwargs...)` integrates the problem forward and returns a
`PyrolysisSolution`. With no options it runs a fixed (Eulerian) mesh, no
volumetric radiation, the `:standard` energy form, and a stiff implicit
integrator. For a cone run we make two deliberate choices: track the receding
surface with ALE mesh motion, and deposit the cone radiation at the surface.

```julia
solution = solve(
    problem;
    abstol = 1e-8,                          # absolute ODE tolerance (default 1e-8)
    reltol = 1e-6,                          # relative ODE tolerance (default 1e-6)
    use_ale = true,                         # move the mesh as the sample shrinks
    radiation_model = SURFACE_ABSORPTION,   # absorb cone flux at the surface
    min_thickness = 0.01e-3,                # terminate when sample < 0.01 mm thick
    handle_depletion = true,                # merge surface cells as they thin out
    depletion_threshold = 0.05,             # merge a cell below 5 % of its initial Δz
    min_cells = 1,                          # never go below 1 cell
    progress = true,                        # show a progress bar
)
```

The first run incurs Julia compilation latency (tens of seconds); subsequent
runs in the same session are fast.

### 2.6.1 The options you just used

| Keyword | Default | What it does here |
|---|---|---|
| `abstol` / `reltol` | `1e-8` / `1e-6` | ODE error tolerances |
| `use_ale` | `false` | enables moving-mesh tracking of surface recession |
| `radiation_model` | `NO_RADIATION` | `SURFACE_ABSORPTION` absorbs the cone flux at the surface |
| `min_thickness` | `1e-6` | ALE termination thickness [m] |
| `handle_depletion` | `false` | merge thin surface cells (requires `use_ale = true`) |
| `depletion_threshold` | `0.05` | cell-thinness fraction that triggers a merge |
| `min_cells` | `2` | floor on cell count before termination |
| `progress` | `true` | progress bar |

This run sets `min_cells = 1` to let the mesh thin all the way down as the
material decomposes fully. The default `min_cells = 2` is more conservative,
keeping at least two cells for stability; use `min_cells = 1` only when you
expect the sample to burn through completely (as for this PMMA slab). The solver
validates `min_cells ≥ 1`, so `1` is the smallest legal value.

### 2.6.2 Choosing the radiation model

The cone heater's incident flux can be deposited four ways; for an opaque solid
like PMMA, surface absorption is the right, fast choice:

| `radiation_model` | Behaviour | Use when |
|---|---|---|
| `NO_RADIATION` (default) | No volumetric absorption; a `RadiativeFluxBC` is still absorbed *at the surface* as a boundary flux | opaque material, radiation handled purely as a surface flux |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the surface (boundary energy balance); O(1), tridiagonal Jacobian | opaque material (optical thickness ≫ 1) |
| `BEER_LAMBERT` | In-depth exponential absorption using component `α`; `RadiativeFluxBC` becomes volumetric (zero boundary flux) | semi-transparent material (polymers, foams) |
| `P1_QUASI_STEADY` | P1 diffusion approximation with re-radiation | optically thick + strong re-radiation (experimental) |

Two important caveats:

- Under `NO_RADIATION` and `SURFACE_ABSORPTION` the `RadiativeFluxBC` is
  absorbed at the surface, so the two give essentially the same surface
  behaviour for an opaque sample; the difference is bookkeeping. Under
  `BEER_LAMBERT` the *same* `RadiativeFluxBC` is instead distributed
  volumetrically (and returns zero as a boundary flux) — switching models
  changes where the energy goes, by design. Set component `α` values for
  `BEER_LAMBERT` or you will get a warning that no absorption occurs. In
  `BEER_LAMBERT` mode the boundary-temperature Newton solver still runs (it
  finds the `T_s` that balances conducted, convective, and re-radiation flux),
  but `RadiativeFluxBC` contributes **zero** flux at the boundary; the absorbed
  fraction of the incident cone radiation — `ε_eff·Φ` by Kirchhoff's law
  (`α = ε_eff` of the live surface cell), or `α·Φ` with an explicit
  `absorptivity` — is deposited volumetrically via
  `compute_beer_lambert_source!` and reported in
  `solution.extras.volumetric_absorption`.
- **`P1_QUASI_STEADY` is experimental** and is *rejected* by the core `solve`
  (it lives in `Pyrolysis.Experimental` with a separate residual path). Do not
  pass it to `solve`.

### 2.6.3 Other solve options worth knowing

These are not needed for this run but round out the picture (full reference in
User Guide §UG-8):

- **`integrator`** (default `KenCarp4`) and **`linear_solver`** (default
  `SparspakFactorization`) — the time integrator and sparse linear solver.
- **`energy_form`** (default `:standard`) — `:standard` carries gas enthalpy via
  the volumetric `S_conv` term; `:with_generation_sink` adds the gas-generation
  enthalpy sink, which is consistent **only with T_ref-datum reaction heats**
  (it double-counts gas sensible enthalpy under the usual DSC-style heats —
  §7.7). Leave at `:standard` unless you have a specific reason.
- **`use_transpiration_bc`** (default `false`) — **leave off.** It double-counts
  gas enthalpy against `S_conv` and `solve` raises an error if set `true`.
- **`use_blowing_correction`** (default `false`) — opt-in Couette-film blowing
  correction: every `ConvectiveBC` is rescaled `h → h·B/(e^B − 1)` with
  `B = ṁ″·c_p,g/h` built from the outward gas flux at the face. Physically real
  and safe with the default energy convention (it modifies the convective
  *coefficient*, not the enthalpy accounting), but enable it on calibrated
  examples only together with flame heat-flux feedback — alone it raises `T_s`
  and MLR relative to the calibration (§5.2.3 of the Technical Reference,
  §12.3.3).
- **`diagnostics`** (default `false`) — snapshot the conservation ledgers every
  `diagnostics_stride` steps (default 10) for energy/mass-balance auditing.
- **`jacobian`** (default `nothing` ⇒ a structured backend with a direct solve)
  — selects the Jacobian backend. **WARNING:** `handle_depletion = true` is a
  **hard incompatibility** with a pluggable `jacobian`, not a soft warning: the
  depletion callback resizes the state vector mid-solve, which the pre-allocated
  Jacobian caches of the structured backend cannot survive. If you specify
  **both** `handle_depletion = true` **and** a `jacobian`, `solve` raises an
  error and refuses to run. To use depletion, **omit the `jacobian` kwarg
  entirely** — the integrator then falls back to its own colored
  finite-difference Jacobian. See Technical Reference §15 and User Guide §UG-11.
- **`warm_start`** — reuse a prior solution's final state as the initial
  condition for a continuation run.

> **Unknown keywords are rejected.** `solve` validates its keyword arguments and
> errors on anything it does not recognise, so a typo is caught immediately
> rather than silently ignored.

---

## 2.7 The solution object

`solve` returns a `PyrolysisSolution`. Its core fields are concretely typed; a
typed `extras` NamedTuple holds opt-in diagnostics. The fields you will use most:

| Field | Type / shape | Meaning |
|---|---|---|
| `solution.t` | `Vector{Float64}` | saved times [s] |
| `solution.T` | `Matrix{Float64}` `(n_cells, n_times)` | temperature profiles [K] |
| `solution.ξ` | `Array{Float64,3}` `(NC, n_cells, n_times)` | concentrations [kg/m³] |
| `solution.z` | `Matrix{Float64}` `(n_nodes, n_times)` | node positions [m] (ALE; `0×0` in Eulerian mode) |
| `solution.surface_T` | `Vector{Float64}` | solved surface temperature per time [K] |
| `solution.mass_loss_rate` | `Vector{Float64}` | MLR per current area (top + permeable-bottom) [kg/(m²·s)] |
| `solution.retcode` | `Symbol` | integrator status (`:Success`, `:Terminated`, …) |
| `solution.solve_time` | `Float64` | wall-clock solve time [s] |

The `extras` NamedTuple adds (among others) `gas_generation_rate`,
`back_temperature`, `total_mass`, `thickness_history`, `MLR_gauge`, and a
per-component surface heat-flux breakdown. A few are central to cone analysis:

| `solution.extras.…` | Meaning |
|---|---|
| `gas_generation_rate` | volumetric gas production rate, per initial area [kg/(m²·s)] |
| `back_temperature` | back-face (`z = 0`) temperature [K] |
| `total_mass` | remaining sample mass per initial area [kg/m²] |
| `thickness_history` | current sample thickness [m] |
| `MLR_gauge` | MLR re-gauged to the *initial* area, `MLR·A(t)/A₀` |

Always check the return code first:

```julia
println("retcode  = ", solution.retcode)     # :Terminated means the sample burned through
println("solve t  = ", round(solution.solve_time, digits = 2), " s")
println("saved    = ", length(solution.t), " time points")
```

For an ALE cone run that consumes the sample, `:Terminated` is the *expected*
outcome (the thickness-termination callback fired); `:Success` means it ran the
full `tspan` without depleting.

### 2.7.1 MLR vs. gas-generation rate (a common gotcha)

There are two related "mass-loss" quantities and they are *not* identical:

- **`solution.mass_loss_rate`** — the actual gas mass flux leaving the sample,
  per current cross-sectional area [kg/(m²·s)]. It counts both faces; with the
  usual impermeable holder the bottom term is zero, so this is the surface
  MLR. (The strictly top-face flux is `solution.extras.surface_mass_flux`.)
- **`solution.extras.gas_generation_rate`** — the *volumetric* gas production
  integrated over the domain, per *initial* area. It is a convenient proxy for
  MLR and is what the shipped cone example plots, but it leads the surface flux
  slightly because gas takes time to escape.

For a typical cone-calorimeter run with a **fixed** cross-sectional area (no
lateral shrinkage), `A(t) = A₀` always, so `mass_loss_rate` and `MLR_gauge` are
identical and the distinction is moot. The distinction becomes important only if
you enable a **lateral-shrinkage law** (the `WithChi` cross-section mode), where
`A(t)` shrinks over time: `MLR_gauge` then automatically re-normalises the
surface flux by `A(t)/A₀`, making it directly comparable to standard
cone-calorimeter data reported per *original* specimen area. For this fixed-area
run (`cross_section_area = 1.0`, no lateral shrinkage law), `mass_loss_rate` and
`MLR_gauge` coincide.

---

## 2.8 Step 6 — Extract surface temperature and mass-loss rate

Everything you need is already on the solution object. The simplest extraction
is direct field access:

```julia
times     = solution.t
T_surface = solution.surface_T                       # [K]
T_back    = solution.extras.back_temperature         # [K]
mlr       = solution.mass_loss_rate .* 1000          # [g/(m²·s)]
gen_rate  = solution.extras.gas_generation_rate .* 1000  # [g/(m²·s)] proxy
total_m   = solution.extras.total_mass               # [kg/m²]

println("peak surface T = ", round(maximum(T_surface), digits = 1), " K")
println("peak MLR       = ", round(maximum(mlr), digits = 2), " g/(m²·s)")
println("mass loss      = ",
        round((1 - total_m[end] / total_m[1]) * 100, digits = 1), " %")
```

To pull a quantity at a particular time, find the nearest saved index and read
the field:

```julia
i600 = argmin(abs.(times .- 600.0))                  # index nearest t = 600 s
println("at t≈600 s:  T_s = ", round(T_surface[i600], digits = 1), " K,",
        "  MLR = ", round(mlr[i600], digits = 2), " g/(m²·s)")
```

### 2.8.1 Temperature profiles through the depth

`solution.T[:, i]` is the temperature of every cell at saved time index `i`
(bottom cell first, surface cell last). For a depth profile at, say, 120 s:

```julia
i = argmin(abs.(times .- 120.0))
T_profile = solution.T[:, i]                         # [K], one value per cell
# Depth axis spans the *current* thickness at this time (ALE mesh shrinks):
L_now  = solution.extras.thickness_history[i]
depths = range(0.0, L_now, length = length(T_profile))  # [m], z = 0 at the back
```

In ALE mode the per-cell node positions are also available in `solution.z[:, i]`
when you need the exact (non-uniform) cell spacing.

### 2.8.2 Convenience accessors (optional)

The Solver submodule provides a few helper accessors —
`get_profile(solution, :T, i)`, `get_cell_positions(solution)`, and
`interpolate_solution(solution, t)` (which interpolates the *raw* state vector
at an arbitrary `t`). These are **not** part of the curated top-level export
list, so `using Pyrolysis` alone does **not** bring them into bare scope. They
live in the `Solver` submodule (re-exported through `Internal`); access them
either by bringing them in with `using Pyrolysis.Internal` and calling them
directly, or with fully-qualified names like `Pyrolysis.Solver.get_profile(...)`:

```julia
using Pyrolysis.Internal     # exposes get_profile, get_cell_positions, interpolate_solution, …

T_at_step = get_profile(solution, :T, i600)   # equivalent to solution.T[:, i600]
state_t   = interpolate_solution(solution, 333.0)  # raw state vector at t = 333 s
```

For everyday use, direct field access (`solution.T`, `solution.surface_T`, …) is
the simplest and most robust path.

---

## 2.9 Step 7 — Plot the results

With `Plots` loaded, the package's recipe gives you quick standard views off the
solution directly:

```julia
using Plots

plot(solution)                      # surface-temperature trace (the default)
plot(solution; kind = :mlr)         # mass-loss-rate curve [kg/(m²·s)]
plot(solution; kind = :heatmap)     # T(z, t) field as a heatmap
plot(solution; kind = :residual)    # energy-conservation residual (zeros unless diagnostics on)
```

The recognised `kind` values are `:surface_T` (default), `:mlr`, `:heatmap`, and
`:residual`; any other value raises an informative error.

For a publication-style multi-panel figure you build the series yourself from the
fields. A typical cone summary — an MLR/gas-generation proxy, the surface and
back temperatures, and temperature profiles at several times:

```julia
# Panel 1: gas-generation-rate proxy for HRR/MLR
p1 = plot(times, gen_rate;
    xlabel = "Time [s]", ylabel = "Gas generation [g/(m²·s)]",
    title = "Gas generation (MLR proxy)", lw = 2, label = "Pyrolysis.jl")

# Panel 2: surface and back temperatures
p2 = plot(times, T_surface;
    xlabel = "Time [s]", ylabel = "Temperature [K]",
    title = "Surface & back temperature", lw = 2, label = "Surface (z = L)")
plot!(p2, times, T_back; lw = 2, label = "Back (z = 0)")

# Panel 3: through-thickness temperature profiles at selected times
p3 = plot(xlabel = "Depth [mm]", ylabel = "Temperature [K]",
          title = "Temperature profiles")
for t_target in (30, 120, 300, 600, 1200)
    i = argmin(abs.(times .- t_target))
    Tp = solution.T[:, i]
    length(Tp) <= 1 && continue                      # skip if the mesh merged to 1 cell
    L_now = solution.extras.thickness_history[i]
    z_mm  = range(0.0, L_now * 1000, length = length(Tp))
    plot!(p3, z_mm, Tp; lw = 1.5, label = "t = $(round(times[i])) s")
end

fig = plot(p1, p2, p3; layout = (1, 3), size = (1200, 380))
savefig(fig, "cone_quickstart.png")
```

> **Note — true HRR needs a combustion-efficiency factor.** The cone HRR is
> usually `HRR = MLR × Δh_c × χ`, where `Δh_c` is the effective heat of
> combustion and `χ` the combustion efficiency. Pyrolysis.jl is a *condensed-
> phase* solver and outputs the mass-loss rate (and a gas-generation proxy); it
> does not model the gas-phase flame, so multiply by your `Δh_c·χ` to get an HRR
> estimate. See Technical Reference §1 for the scope boundary.

---

## 2.10 Complete script

The full tutorial in one runnable block. (Aside from the explanatory comments,
this mirrors `examples/06_cone_calorimeter_surface_absorption.jl`.)

```julia
using Pyrolysis
using Plots

# ---- Step 1: material -------------------------------------------------------
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2 * T) : (1040.0 + 3.08 * T))
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * T))

virgin = SolidComponent(:virgin, ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96, α = 1.47)
char   = SolidComponent(:char,   ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96, α = 1.47)
mma    = GaseousComponent(:MMA,  M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5)

decomp = Reaction(:decomposition, 1 => (2, 3);
                  A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0, yields = (0.002, 0.998))

pmma = Material(
    name = :PMMA,
    components = (virgin, char, mma),
    reactions  = (decomp,),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# ---- Step 2: boundary conditions -------------------------------------------
q_external = 25e3
h_conv     = 7.2
T_ambient  = 300.0

top_thermal =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +
    RadiativeBC(ε = 0.96, T_inf = 330.0) +
    ConvectiveBC(h = h_conv, T_inf = 330.0)
bottom_thermal =
    RadiativeBC(ε = 0.95, T_inf = 307.0) +
    ConvectiveBC(h = 4.0, T_inf = 307.0)

top_mass    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)
bottom_mass = ImpermeableBC()

bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal, bottom = bottom_thermal),
    mass    = (top = top_mass,    bottom = bottom_mass),
)

# ---- Step 3: mesh -----------------------------------------------------------
sample_thickness = 5.8e-3
n_cells          = 72
mesh = create_uniform_mesh(sample_thickness, n_cells, Val(3);
                           T_initial = T_ambient, ξ_initial = (1210.0, 0.0, 0.0))

# ---- Step 4: problem --------------------------------------------------------
t_final = 1600.0
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    T_initial = z -> T_ambient,
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, t_final),
    dt_initial = 1e-3,
    output_times = collect(range(0.0, t_final, length = 601)),
)

# ---- Step 5: solve ----------------------------------------------------------
solution = solve(problem;
    abstol = 1e-8, reltol = 1e-6,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    min_thickness = 0.01e-3,
    handle_depletion = true, depletion_threshold = 0.05, min_cells = 1,
    progress = true,
)

# ---- Step 6: extract --------------------------------------------------------
times     = solution.t
T_surface = solution.surface_T
mlr       = solution.mass_loss_rate .* 1000
gen_rate  = solution.extras.gas_generation_rate .* 1000
println("retcode = ", solution.retcode)
println("peak T_s = ", round(maximum(T_surface), digits = 1), " K")
println("peak MLR = ", round(maximum(mlr), digits = 2), " g/(m²·s)")

# ---- Step 7: plot -----------------------------------------------------------
plot(solution; kind = :mlr)
```

---

## 2.11 What to do next

- **Switch to in-depth radiation.** Set `radiation_model = BEER_LAMBERT` and
  give the solid components a non-zero absorption coefficient `α` to model a
  semi-transparent polymer. Compare the surface-temperature trace to the
  surface-absorption run above. (Technical Reference §8.)
- **Add pressure-driven transport.** Give a component a permeability `κ > 0` to
  switch the solver into Darcy–Fick mode for pressure-treated wood or charring
  systems. (User Guide §UG-7; Technical Reference §9.)
- **Vary the cone flux.** Replace the scalar `q_external` with a function of
  time, `RadiativeFluxBC(flux = t -> ramp(t), absorptivity = 0.96)`; the solver
  re-samples the flux at every RHS evaluation.
- **Sweep parameters.** Use `parameter_sweep` to run a grid of heat fluxes or
  kinetic parameters (User Guide §UG-12), or `forward_sensitivity` /
  `adjoint_sensitivity` to compute gradients of an output with respect to
  `A`, `E`, `h`, or material properties (User Guide §UG-13).
- **Audit conservation.** Re-run with `diagnostics = true` and inspect
  `solution.extras.diagnostics_log` to verify the energy and mass ledgers close.
  (Technical Reference §16.)

You now have the full pipeline — material, boundary conditions, mesh, problem,
solve, solution — and can recombine these pieces for any of the worked cases in
User Guide §UG-14.


---

# 3. Defining Materials, Components, and Properties

A `Material` is the physical heart of every Pyrolysis.jl simulation: it bundles the
condensed-phase and gaseous **components** that make up your sample, the
**reactions** that convert them into one another, the **property functions** that
give each component its temperature- (and optionally composition-) dependent
thermophysical behaviour, and the **mixing rules** that homogenize those component
properties into the effective bulk properties the solver actually integrates. This
chapter shows how to build each of these pieces with the real constructors, what
every field and default means, and how to assemble and validate a complete
material. The physics and derivations behind these objects live in the Technical
Reference (see Technical Reference §4 *Materials & Property Functions* and §5
*Effective Properties & Mixing Rules*); here the focus is on *how to specify* them
correctly.

Throughout, the spatial coordinate is `z` (through-thickness): `z = 0` is the
bottom/substrate, `z = L` is the top/exposed surface (see Technical Reference §1).
The primary species state variable is the **mass concentration** `ξ_j` in kg·m⁻³
(not mass fraction). All inputs are SI: densities in kg·m⁻³, heat capacities in
J·kg⁻¹·K⁻¹, conductivities in W·m⁻¹·K⁻¹, activation energies in J·mol⁻¹.

---

## 3.1 Overview: the component → material hierarchy

The build order is always the same:

1. **Components** (`SolidComponent`, `LiquidComponent`, `GaseousComponent`) — each
   carries its phase as a *type*, so the solver dispatches on phase at compile
   time. Solids and liquids form the matrix; gases flow through it.
2. **Property functions** — every per-component thermophysical property (`ρ`, `c`,
   `k`, `λ`, `κ`, intrinsic density) is a callable of temperature. You may pass a
   bare number, a `(a, b)` tuple, a function, or an explicit property-function
   object; the constructors convert for you.
3. **Reactions** (`Reaction`) — single-reactant Arrhenius decomposition steps with
   product yields, activation energy, and heat of reaction. Components may be
   referenced by **symbol** (resolved automatically) or by integer index.
4. **Mixing rules** (`MixingSpec`, `MaterialMixing`) — pick how component thermal
   conductivity `k`, gas-transfer coefficient `λ`, and permeability `κ` are
   homogenized.
5. **Material** (`Material`) — the assembled object, optionally checked with the
   internal helper `Pyrolysis.Materials.validate_material` (or simply validated by
   `solve` when you build the problem; see §3.7).

A minimal end-to-end example (one solid, one gas, one reaction) appears in the
package's top-level docstring and in Chapter 2; this chapter expands every option.

> **A note on what is exported.** `using Pyrolysis` brings the *constructors* into
> scope: `SolidComponent`, `LiquidComponent`, `GaseousComponent`, `Material`,
> `Reaction`, `ReactionSet`, all property-function types, the `MixingRule` enum
> values, `MixingSpec`, and `MaterialMixing`. The various **helper and inspection
> functions** used in this chapter are *not* part of the curated top-level surface,
> so you reach them by qualifying through the submodule that defines them. This works
> the same way whether or not the function is exported from that submodule — Julia's
> `Module.name` syntax resolves any name defined in the module:
>
> - From the `Materials` submodule: `Pyrolysis.Materials.n_components`,
>   `Pyrolysis.Materials.n_reactions`, `Pyrolysis.Materials.dry_pyrolysis_gas_indices`,
>   `Pyrolysis.Materials.gas_component_indices`, and `Pyrolysis.Materials.has_darcy_flow`
>   are *exported* from `Materials`; the validation helpers
>   `Pyrolysis.Materials.validate_material`,
>   `Pyrolysis.Materials.verify_stoichiometry`, and
>   `Pyrolysis.Materials.check_material_mass_balance` are *internal* (not exported, not
>   part of the stable public API) but are still callable through the same qualified
>   path.
> - From the `Properties` submodule: the effective-property functions
>   (`Pyrolysis.Properties.mixture_density`, etc.).
>
> For convenience you may `using Pyrolysis.Internal`, which re-exports the union of
> every submodule's *exports* and so brings the exported helpers
> (`n_components`, `dry_pyrolysis_gas_indices`, `has_darcy_flow`, …) into scope
> unqualified — but **not** the internal validation helpers above, which still
> require the `Pyrolysis.Materials.` prefix. `Internal` is not the stable public
> surface; prefer qualified access in production code. If you want validation as part
> of the public API, build a `PyrolysisProblem` and let `solve` validate (it runs the
> same checks by default, `validate = true`); see §3.7.

---

## 3.2 Defining components

### 3.2.1 `SolidComponent`

A solid is the canonical matrix-forming component. The keyword constructor is:

```julia
SolidComponent(name::Symbol;
               ρ, c, k,
               ε = 0.9, α = 0.0, γ = 1.0,
               λ = 0.0, κ = 0.0,
               ρ_intrinsic = nothing, φ = nothing)
```

| Keyword | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `name` | — | Component identifier (positional, a `Symbol`) | — | required |
| `ρ` | `ρ_j` | **Bulk** (pure-phase) density | kg·m⁻³ | required |
| `c` | `c_{p,j}` | Specific heat capacity | J·kg⁻¹·K⁻¹ | required |
| `k` | `k_j` | Thermal conductivity | W·m⁻¹·K⁻¹ | required |
| `ε` | `ε_j` | Surface emissivity | – | `0.9` |
| `α` | `α_j` | Mass-basis absorption coefficient | m²·kg⁻¹ | `0.0` |
| `γ` | `γ_j` | Swelling factor (volume-change participation) | – | `1.0` |
| `λ` | `λ_j` | Gas-transfer (diffusion) coefficient contribution | m²·s⁻¹ | `0.0` |
| `κ` | `κ_j` | Darcy permeability | m² | `0.0` |
| `ρ_intrinsic` | `ρ_{i,j}` | Skeletal/intrinsic density (porosity only) | kg·m⁻³ | `= ρ` |
| `φ` | — | Pure-phase porosity (alternative to `ρ_intrinsic`) | – | `nothing` |

```julia
using Pyrolysis

# Virgin wood: temperature-rising heat capacity, opaque, tight matrix
virgin = SolidComponent(:virgin,
    ρ = 530.0,            # bulk density [kg/m³]
    c = (1500.0, 1.0),    # LinearProperty: cp(T) = 1500 + 1.0·T  [J/(kg·K)]
    k = 0.12,             # ConstantProperty 0.12 W/(m·K)
    ε = 0.88,             # emissivity
    λ = 2e-5,             # contributes to mixture gas-transfer coefficient
    κ = 1e-15)            # very low permeability (tight cellular structure)

# Char: lighter, more conductive, more emissive, absorbing, porous
char = SolidComponent(:char,
    ρ = 150.0, c = 1000.0, k = 0.10,
    ε = 0.95,
    α = 5.0,              # mass absorption coeff (used only by in-depth radiation)
    λ = 2e-5,
    κ = 1e-12)            # open, porous char conducts gas readily
```

Notes on the less-obvious fields:

- **`ε` and `α` are plain `Float64`s, not property functions.** Emissivity and the
  absorption coefficient do not carry temperature dependence at the component level.
- **`α` (absorption) is ignored unless an in-depth radiation model is active.** With
  `radiation_model = SURFACE_ABSORPTION` (or `NO_RADIATION`) the absorbed flux `α·Φ`
  is applied at the surface through the boundary energy balance and the component `α`
  values play no role — in
  fact `solve` will emit a warning if you set non-zero `α` together with
  `SURFACE_ABSORPTION`. Use `BEER_LAMBERT` (or P1, experimental) to make `α` matter.
  See §3.6 and Technical Reference §8.
- **`λ` on a solid** does not transport the solid (solids are stationary); it is the
  solid's *contribution* to the effective mixture gas-transfer coefficient through
  the mixing rule (see §3.5). ThermaKin-style setups typically use `λ ≈ 2e-5`.
- **`κ = 0.0` (default) means impermeable.** Setting any solid/liquid component's
  `κ > 0` switches the whole problem into Darcy–Fick (pressure-driven) gas transport
  (§3.6.2 and Technical Reference §9).
- **Swelling factor `γ`** is the fraction of the component's volume that counts as
  bulk (matrix) volume, with **one meaning everywhere**: the dilation rate θ, the
  mixture density, and the conductivity/emissivity volume weights all use the same
  `γ_j·ξ_j/ρ_j` weighting. `γ = 1.0` (default for solids) means the solid shrinks as
  it is consumed; `γ = 0.0` makes the component consistently *volumeless* — no volume
  change, no volume in the bulk density, and **no weight in conductivity/emissivity
  mixing** (use with care for a "rigid char"); intermediate values are for
  intumescent/partial-expansion cases. Values outside [0, 1] are rejected at
  construction. See Technical Reference §5.1–5.2 and §10.

### 3.2.2 `LiquidComponent`

Liquids are matrix-forming like solids (same field set) and are the standard way to
represent **bound/free moisture** in hygroscopic materials. The defaults differ from
solids: a liquid by default is strongly emissive (ε = 0.95, near water's ≈ 0.96),
non-swelling, and blocks gas flow. (For the default γ = 0 the emissivity carries no
weight in the mixture anyway — a γ = 0 liquid is volumeless in the emissivity mixing;
the 0.95 matters as soon as you give the liquid γ > 0.)

```julia
LiquidComponent(name::Symbol;
                ρ, c, k,
                λ = 2e-5, κ = 0.0,
                ε = 0.95, α = 0.0, γ = 0.0,
                ρ_intrinsic = nothing, φ = nothing)
```

| Keyword | Meaning | Units | Default |
|---|---|---|---|
| `ρ`, `c`, `k` | bulk density, heat capacity, conductivity | (as solid) | required |
| `λ` | gas-transfer contribution | m²·s⁻¹ | `2e-5` |
| `κ` | permeability | m² | `0.0` (liquid blocks pores) |
| `ε` | emissivity | – | `0.95` (water ≈ 0.96) |
| `α` | mass absorption | m²·kg⁻¹ | `0.0` |
| `γ` | swelling factor | – | `0.0` |
| `ρ_intrinsic` / `φ` | skeletal density / pure-phase porosity | kg·m⁻³ / – | `= ρ` / `nothing` |

```julia
# Liquid water (moisture) in wood
moisture = LiquidComponent(:moisture,
    ρ = 1000.0, c = 4180.0, k = 0.60)   # λ defaults to 2e-5, κ=0 (blocks flow)
```

The default `γ = 0.0` for liquids means moisture evaporation does **not** by itself
drive bulk shrinkage; set `γ = 1.0` if you want liquid loss to shrink the matrix.

### 3.2.3 `GaseousComponent`

Gases flow *through* the matrix; they are the products of decomposition reactions.
A gaseous component has **no permeability field** (`get_permeability` returns 0 for
gases) and carries a **molar mass** instead, which feeds the ideal-gas density and
the pore-pressure closure `P = Σ_g ξ_g R_g T/(M_g φ)`.

```julia
GaseousComponent(name::Symbol;
                 M, c, k, λ,
                 ε = 0.0, α = 0.0, γ = 0.0)
```

| Keyword | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `M` | `M_j` | Molar mass | kg·mol⁻¹ | required |
| `c` | `c_{p,j}` | Heat capacity | J·kg⁻¹·K⁻¹ | required |
| `k` | `k_j` | Thermal conductivity | W·m⁻¹·K⁻¹ | required |
| `λ` | `λ_j` | Gas-transfer (diffusion) coefficient | m²·s⁻¹ | required |
| `ε` | `ε_j` | Emissivity | – | `0.0` |
| `α` | `α_j` | Mass absorption | m²·kg⁻¹ | `0.0` |
| `γ` | `γ_j` | Swelling factor | – | `0.0` |

```julia
# Pyrolysis gas (modelled as a CO₂-like vapour) and water vapour
gas   = GaseousComponent(:gas, M = 0.044, c = 1800.0, k = 0.03, λ = 1e-4)
vapor = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025, λ = 2e-5)
```

**Density is automatic.** Unlike solids/liquids, you do *not* pass `ρ` to a gas. The
constructor installs the ideal-gas law

```
ρ_g(T) = P_ref · M / (R_g · T)
```

with `P_ref = 101325 Pa` and `R_g = 8.31446261815324 J·mol⁻¹·K⁻¹`. For example a
`M = 0.044` gas gives `ρ_g(300) ≈ 1.787 kg·m⁻³`. The ideal-gas default is correct for
almost all pyrolysis applications. See Technical Reference §4.3.

**Naming and water vapour.** Components named so that the name contains `H2O` or
`water` (case-insensitive, e.g. `:H2O`, `:water_vapor`) are treated as *moisture*
gases and are **excluded** from the dry-pyrolysis-gas set that drives lateral
shrinkage `χ` (Technical Reference §10.3). Name your dry-pyrolysis gas anything else
(`:gas`, `:MMA`, `:volatiles`). You can inspect the classification:

```julia
Pyrolysis.Materials.gas_component_indices(material)        # all gas indices
Pyrolysis.Materials.dry_pyrolysis_gas_indices(material)    # excludes :H2O / :water*
```

### 3.2.4 Bulk vs. intrinsic density and porosity

Two densities coexist, and getting them right matters for porous-media physics:

- **Bulk density `ρ_j`** is the pure-phase density used by *mixing rules*, *mixture
  density*, and *volume-change* kernels. This is what you pass as `ρ`.
- **Intrinsic (skeletal) density `ρ_{i,j}`** is the density of the cell-wall
  material itself, ignoring pores within that component. It is used **only** by the
  porosity formula

  ```
  φ = clamp(1 − Σ_{j∈solid,liquid} ξ_j / ρ_{i,j}, 0, 1)
  ```

  and, through `φ`, by the gas-pressure closure `P = N/φ` and the Carman–Kozeny
  permeability model (Technical Reference §5.3, §9.4).

If you never supply an intrinsic density, it **defaults to the bulk density**, which
makes the pure phase pore-free and gives `φ = 0` for a fully solid cell — the
historical, backward-compatible default. To model a genuinely porous solid you must
supply *either* `ρ_intrinsic` *or* `φ` (not both — passing both raises an
`ArgumentError`):

```julia
# Option A: give the skeletal density directly
porous_char = SolidComponent(:char, ρ = 150.0, c = 1000.0, k = 0.10,
                             ρ_intrinsic = 1450.0)   # cell-wall density of carbon

# Option B: give the pure-phase porosity (only valid with a constant bulk ρ)
porous = SolidComponent(:porous, ρ = 400.0, c = 1000.0, k = 0.10, φ = 0.5)
# → intrinsic_density = ρ / (1 − φ) = 400 / 0.5 = 800 kg/m³
```

The `φ` shortcut requires a **constant** bulk `ρ`; for a temperature-dependent bulk
density you must pass `ρ_intrinsic` explicitly. `φ` must satisfy `0 ≤ φ < 1`.

> **Practical implication.** If you intend to use Carman–Kozeny permeability or the
> gas-pressure formula and you leave intrinsic density at its default, porosity will
> be ~0 and those closures will behave degenerately. Always set `ρ_intrinsic`/`φ`
> on the solids when porosity matters.

### 3.2.5 Inspecting a component

Phase predicates and accessors dispatch statically on the component type:

```julia
is_solid(virgin)        # true
is_gas(gas)             # true
is_liquid(moisture)     # true
phase(virgin)           # SOLID  (a Phase enum value)
molar_mass(gas)         # 0.044  (0.0 for solids/liquids)

# Density / property accessors (qualify through Materials):
Pyrolysis.Materials.get_density(virgin, 300.0)            # bulk ρ at T
Pyrolysis.Materials.get_intrinsic_density(porous, 300.0)  # skeletal ρ
Pyrolysis.Materials.get_permeability(char, 500.0)         # κ at T (0 for gases)
```

---

## 3.3 Property functions

Every per-component thermophysical property is an `AbstractPropertyFunction`: a
callable `p(T)` (and, for some, `p(T, ξ)`). You rarely construct these directly —
the component constructors call `to_property` on whatever you pass — but knowing the
types lets you express any temperature dependence you need.

### 3.3.1 Shorthand conversions (`to_property`)

`to_property` (internal: `Pyrolysis.Materials.to_property`) maps convenient inputs to
property objects, and the component/reaction constructors apply it automatically to
`ρ`, `c`, `k`, `λ`, `κ`, `ρ_intrinsic`, and `h`:

| You pass | You get | Meaning |
|---|---|---|
| a `Real` number, e.g. `0.21` | `ConstantProperty(0.21)` | constant in `T` |
| a 2-tuple `(a, b)` | `LinearProperty(a, b)` | `a + b·T` |
| a 1-arg function `T -> …` | `CallableProperty` | arbitrary `f(T)` |
| an `AbstractPropertyFunction` | passes through unchanged | explicit object |

So `k = 0.21`, `k = (0.20, 1e-4)`, `k = T -> 0.20 + 1e-4*T`, and
`k = LinearProperty(0.20, 1e-4)` are four equivalent ways to specify a conductivity.

### 3.3.2 The property-function types

All are exported by `using Pyrolysis`.

**`ConstantProperty(value)`** — `p(T) = value`.

```julia
k = ConstantProperty(0.21)
```

**`LinearProperty(a, b)`** — `p(T) = a + b·T`. Equivalent to the `(a, b)` shorthand.

```julia
cp = LinearProperty(1420.0, 2.5)   # cp(T) = 1420 + 2.5·T
```

**`PolynomialProperty(coeffs::NTuple)`** — `p(T) = c₀ + c₁T + c₂T² + …`, evaluated
with `evalpoly` (numerically stable, Horner form). Pass the coefficients
low-order-first as a tuple:

```julia
cp = PolynomialProperty((1000.0, 2.0, 0.001))   # 1000 + 2T + 0.001T²
cp(300.0)   # = 1690.0
```

**`ArrheniusProperty(A, E[, T_ref])`** — `p(T) = A · exp(−E / (R_g·T))`, `E` in
J·mol⁻¹. `T_ref` defaults to 298.15 K. Mostly used for rate-constant-shaped
properties; reaction kinetics themselves are specified on the `Reaction` (§3.4), not
as an `ArrheniusProperty`.

```julia
p = ArrheniusProperty(8.5e12, 188e3)   # T_ref defaults to T_REF = 298.15
```

**`TableProperty(temperatures, values)`** — tabulated `(Tᵢ, vᵢ)` with linear
interpolation. Both arguments are tuples; **temperatures must be strictly
increasing** (the constructor errors otherwise). Queries outside the table
**extrapolate flat** (clamp to the nearest endpoint value). A binary search is used
for tables with more than 8 entries.

```julia
k = TableProperty((300.0, 400.0, 500.0), (0.20, 0.22, 0.25))
k(350.0)   # = 0.21  (linear interpolation between 0.20 and 0.22)
k(250.0)   # = 0.20  (flat extrapolation below the table)
```

**`CallableProperty(func)`** — wraps any `f(T) -> value`. Use for closed-form
dependences not covered by the other types:

```julia
k = CallableProperty(T -> 0.10 + 3.0e-4 * (T - 300.0))
```

**`StateDependentProperty(func)`** — a property that depends on both temperature and
the full cell concentration vector, `f(T, ξ) -> value`. This is the mechanism for
**composition-dependent heats of reaction**, most importantly the moisture heat of
sorption, where the sorption enthalpy rises as moisture content falls.

```julia
# Heat of sorption as a function of moisture content (%); ξ is the full
# mass-concentration tuple, so you index the moisture component yourself.
MOISTURE_IDX = 3
heat_sorption = StateDependentProperty((T, ξ) -> begin
    ξ_solid = ξ[1] + ξ[2]                              # dry-solid concentration
    MC = 100.0 * ξ[MOISTURE_IDX] / max(ξ_solid, 1.0)   # moisture content [%]
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3          # [J/kg], endothermic (h > 0)
end)
```

A `StateDependentProperty` **must** be called with state: `p(T, ξ)`. Calling it with
only `T` raises an error. When used as a reaction heat, the kinetics evaluation
passes `ξ` automatically. See §3.4.5 and Technical Reference §6.6.

### 3.3.3 Derivatives (for the implicit Jacobian)

You generally do not call these, but it is useful to know they exist and how they
behave. Temperature derivatives `dp/dT` are exact and closed-form for
`ConstantProperty` (0), `LinearProperty` (`b`), `PolynomialProperty` (unrolled), and
`ArrheniusProperty` (`p·E/(R_g T²)`); they are finite-difference (`δT = 0.1 K`) for
`TableProperty`, `CallableProperty`, and `StateDependentProperty`. Concentration
derivatives `∂p/∂ξ_j` are zero for temperature-only properties and adaptive
finite-difference (`δξ = max(|ξ|·1e-4, 1e-6)`) for `StateDependentProperty`. Because
table/callable/state-dependent properties use finite differences, prefer the
analytic types (constant/linear/polynomial/Arrhenius) where you can — they keep the
Jacobian exact and the solve crisp. See Technical Reference §4.5.

---

## 3.4 Defining reactions

Reactions convert one reactant component into one or two product components by
Arrhenius kinetics. They are usually attached to a material when you build it. The
full physics — smooth temperature gates, the depletion limiter, species sources — is
in **Chapter 4** and Technical Reference §6; here we cover the constructor surface
enough to build a material.

### 3.4.1 Constructor forms

There are four `Reaction` constructors. They differ only in how reactant/product are
named (integer index vs. symbol) and whether there are one or two products.

```julia
# (a) integer index, single product
Reaction(name, reactant::Int => product::Int;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, validate_mass = true)

# (b) symbol, single product (resolved by Material)
Reaction(name, :reactant => :product;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, validate_mass = true)

# (c) integer index, two products
Reaction(name, reactant::Int => (p1::Int, p2::Int);
         A, E, h, yields, n = 1.0, T_min = 0.0, T_max = Inf, validate_mass = true)

# (d) symbol, two products (resolved by Material)
Reaction(name, :reactant => (:p1, :p2);
         A, E, h, yields = (y1, y2), n = 1.0, T_min = 0.0, T_max = Inf,
         validate_mass = true)
```

| Argument | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `name` | — | reaction identifier (`Symbol`, positional) | — | required |
| `A` | `A_i` | Arrhenius pre-exponential | s⁻¹ | required |
| `E` | `E_i` | activation energy | J·mol⁻¹ | required |
| `h` | `h` | heat of reaction (number or property fn) | J·kg⁻¹ | required |
| `yields` | `ν_{i,j}` | per-product mass yields (two-product forms) | – | required for (c)/(d) |
| `n` | `n_{i,j}` | reaction order in the reactant (must be ≥ 0) | – | `1.0` |
| `T_min` | `T_min,i` | lower temperature gate | K | `0.0` |
| `T_max` | `T_max,i` | upper temperature gate | K | `Inf` |
| `validate_mass` | — | enforce `Σ yields = 1` at construction | – | `true` |

Yields must each be ≥ 0 and the order must satisfy `n ≥ 0` — both are enforced
unconditionally at construction; `validate_mass = false` disables only the
`Σ yields = 1` sum check.

### 3.4.2 Symbolic vs. index references

The **symbol** forms ((b), (d)) are the recommended style: you name components, and
the `Material` constructor resolves the names to indices using the `components`
tuple, so you never track integer positions by hand. Internally a symbolic reaction
is stored as a `PendingReaction` until `Material(...)` resolves it.

```julia
pyrolysis = Reaction(:pyrolysis, :virgin => (:char, :gas),
                     A = 1e10, E = 140e3, h = 500e3,   # h > 0: endothermic
                     yields = (0.25, 0.75),            # 25% char + 75% gas
                     T_min = 450.0)                    # gate below 450 K
```

The **index** forms ((a), (c)) reference components by their 1-based position in the
material's `components` tuple. They are handy in programmatically generated schemes:

```julia
rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3)
# component 1 → component 2, single product, default order/gates
```

### 3.4.3 Heat-of-reaction sign convention

Pyrolysis.jl uses the **storage convention `h > 0 = endothermic`** (cools the
material). The volumetric source is `Q_rxn = −Σ h_r r_r`, so an endothermic reaction
produces `Q_rxn < 0`. This is the *opposite* of the ThermaKin/Gpyro publication
convention (`h > 0 = exothermic`); when you transcribe literature parameters, flip
the sign. Examples:

```julia
endo = Reaction(:pyrolysis, :virgin => :gas, A = 1e10, E = 140e3, h =  500e3)  # cools
exo  = Reaction(:char_ox,   :char   => :ash, A = 1e6,  E = 120e3, h = -300e3)  # heats
```

See the overload note H1 in the nomenclature and Technical Reference §6.6.

### 3.4.4 Yields and mass balance

For single-product reactions the product yield is fixed at `1.0` (one unit of
reactant → one unit of product). For two-product reactions you supply
`yields = (y1, y2)`, and **mass conservation requires `y1 + y2 = 1.0`**. With the
default `validate_mass = true`, the constructor enforces this and errors on a
violation; you can bypass the check with `validate_mass = false` (not recommended).

```julia
# Charring: 1 kg virgin → 0.20 kg char + 0.80 kg gas
Reaction(:charring, :virgin => (:char, :gas),
         A = 1e15, E = 200e3, h = 500e3, yields = (0.20, 0.80))
```

Only **single-reactant** reactions are supported (the reactant is consumed with mass
stoichiometry 1.0). Multi-step schemes are modelled as *sequences* of single-reactant
reactions sharing intermediate components — see Chapter 4.

### 3.4.5 State-dependent heat (moisture sorption)

Pass a `StateDependentProperty` (§3.3.2) as `h` to make the heat depend on local
composition. This is the standard way to model the moisture heat of sorption:

```julia
heat_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[3] / max(ξ[1] + ξ[2], 1.0)
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3
end)

drying = Reaction(:drying, :moisture => :H2O,
                  A = 5.6e8, E = 43e3, h = heat_sorption)
```

### 3.4.6 Verifying reactions

Two internal helpers (not exported; qualify through `Materials`) check mass balance
after the fact:

```julia
Pyrolysis.Materials.verify_stoichiometry(rxn)                  # Bool, warns on failure
Pyrolysis.Materials.check_material_mass_balance(mat.reactions) # all reactions
```

Both compare the reactant mass (1.0) against the sum of product yields within a
tolerance of `1e-10` and emit a warning (returning `false`) on imbalance. They are
also invoked automatically by `validate_material` (§3.7) and by `solve`'s
problem-level validation, so for routine use you can rely on those rather than
calling these private helpers directly.

---

## 3.5 Mixing rules: from component to effective properties

The solver integrates **effective** (homogenized) bulk properties, not individual
component properties. Three properties are mixed per material — thermal conductivity
`k`, gas-transfer coefficient `λ`, and permeability `κ` — and **each may use a
different rule**. Heat capacity, mixture density, emissivity, and absorption are not
configurable: they always use fixed weighting rules (concentration- or
volume-weighted; see Technical Reference §5).

### 3.5.1 `MixingRule` and `MixingSpec`

The available rules are the `MixingRule` enum values (all exported):

| Rule | Formula | Valid channels |
|---|---|---|
| `PARALLEL` | `k = Σ_j v_j k_j` (volume-weighted, upper bound) | k, λ, κ |
| `SERIES` | `1/k = Σ_j v_j/k_j` (harmonic, lower bound) | k, λ, κ |
| `WEIGHTED` | `k = β·k_∥ + (1−β)·k_series` | k, λ, κ |
| `BRUGGEMAN` | EMT: `Σ_j v_j (k_j−k)/(k_j+2k) = 0` | **k only** |
| `CARMAN_KOZENY` | `κ = min(ℓ²φ³ / (180(1−φ)²), ℓ²)` (porosity-based) | **κ only** |

A `MixingSpec` bundles a rule with its parameters:

```julia
MixingSpec(rule::MixingRule, β::Real = 0.5; length_scale::Real = NaN)
```

- `β` is the blend weight, used **only** by `WEIGHTED`, and must lie in `[0, 1]`
  (default `0.5`). It is harmlessly carried for other rules.
- `length_scale` (`ℓ`, the characteristic particle/pore diameter in metres) is used
  **only** by `CARMAN_KOZENY` and is **required** for it — the constructor errors if
  you request `CARMAN_KOZENY` without a positive, finite `length_scale`.

```julia
MixingSpec(WEIGHTED, 0.6)                       # 60% parallel / 40% series for k
MixingSpec(PARALLEL)                            # β ignored
MixingSpec(CARMAN_KOZENY; length_scale = 5e-5)  # ℓ = 50 µm
MixingSpec(BRUGGEMAN)                            # EMT (k channel only)
```

### 3.5.2 `MaterialMixing`

`MaterialMixing` builds the `(k, λ, κ)` NamedTuple that `Material` expects. Each
channel accepts either a `MixingSpec` or a bare `MixingRule` (β defaults to 0.5 for
the bare form). The **defaults** are:

```julia
MaterialMixing(; k = MixingSpec(WEIGHTED, 0.5),   # conductivity
                 λ = MixingSpec(PARALLEL, 0.5),   # gas transfer
                 κ = MixingSpec(WEIGHTED, 0.5))   # permeability
```

```julia
mixing = MaterialMixing(
    k = MixingSpec(WEIGHTED, 0.5),                       # blended bounds
    λ = MixingSpec(PARALLEL),                            # parallel
    κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),  # porosity-based κ
)

# bare-rule shorthand (β = 0.5):
mixing2 = MaterialMixing(k = WEIGHTED, λ = PARALLEL, κ = SERIES)
```

If you omit `mixing` from `Material(...)`, the defaults above are used.

### 3.5.3 Choosing rules — guidance

- **Thermal conductivity `k`.** The default `WEIGHTED` with `β = 0.5` (midway between
  the parallel upper bound and series lower bound) is a sensible neutral choice. Use
  `BRUGGEMAN` for char-forming materials where the microstructure is random/
  self-similar — it captures percolation-like drops as a continuous phase breaks up
  (Technical Reference §5.5). `BRUGGEMAN` is solved by a warm-started Newton
  iteration and is valid **only for `k`**; requesting it for `λ` or `κ` raises an
  error at evaluation time.
- **Gas transfer `λ`.** `PARALLEL` (the default) is appropriate; EMT (`BRUGGEMAN`)
  is rejected for diffusivities because there is no physical EMT closure for them.
- **Permeability `κ`.** Two mutually exclusive approaches:
  1. **Component-based** (`PARALLEL`/`SERIES`/`WEIGHTED`): mix the per-component `κ`
     values you set on the solids/liquids.
  2. **Carman–Kozeny** (`CARMAN_KOZENY`): ignore component `κ` and compute
     permeability from **porosity** via `κ = min(ℓ²φ³/(180(1−φ)²), ℓ²)` — the
     packed-bed law capped at the free-flow scale `ℓ²` for very open structures
     (`φ ≳ 0.93`). This captures the steep rise in permeability as virgin → char
     raises porosity.

> **Important caveat for Carman–Kozeny.** Whether the solver runs in Darcy–Fick
> (pressure-driven) mode is decided by `has_darcy_flow(material)`, which checks
> whether any *component* carries `κ > 0` — it does **not** look at the κ mixing
> rule. So choosing `CARMAN_KOZENY` alone does **not** turn on Darcy transport: with
> all component `κ = 0`, `has_darcy_flow` returns `false` and gas transport stays
> Fickian. If you want pressure-driven flow *and* a Carman–Kozeny effective
> permeability, set a non-zero component `κ` (to trip `has_darcy_flow`) while using
> the `CARMAN_KOZENY` rule for the effective value, or use component-based κ mixing.
> Also remember that Carman–Kozeny needs a real porosity, so set `ρ_intrinsic`/`φ`
> on the solids (§3.2.4) — otherwise `φ ≈ 0` and `κ ≈ 0`.

### 3.5.4 Inspecting effective properties

The effective-property functions live in the `Properties` submodule (qualify them).
They take `(ξ, T, material)` where `ξ` is the mass-concentration tuple:

```julia
using Pyrolysis
P = Pyrolysis

ξ = (300.0, 50.0, 0.5)   # (virgin, char, gas) concentrations [kg/m³]
T = 600.0

P.Properties.mixture_density(ξ, T, material)         # ρ with swelling
P.Properties.effective_heat_capacity(ξ, T, material) # ρc_p (matrix only)
P.Properties.effective_conductivity(ξ, T, material)  # k_eff
P.Properties.effective_gas_transfer(ξ, T, material)  # λ_eff
P.Properties.effective_permeability(ξ, T, material)  # κ_eff
P.Properties.effective_emissivity(ξ, T, material)    # ε_eff
P.Properties.effective_absorption(ξ, material)       # α_eff (note: no T arg)
P.Properties.porosity(ξ, T, material)                # φ
```

Two physics notes worth remembering (both in Technical Reference §5):

- **Heat capacity excludes gas storage.** `effective_heat_capacity` sums only
  solid/liquid components: `ρc_p^eff = Σ_{j∈S,L} ξ_j c_{p,j}(T)`. Gas sensible heat
  is carried by the advective source term instead (the quasi-steady-gas
  approximation, ~0.2% of the matrix term at 1 atm). `total_heat_capacity` (which
  includes gas) exists for diagnostics only.
- **Conductivity and emissivity use the swelling-weighted (matrix) volume
  fractions** — γ = 0 components (canonical gases) carry no weight in these
  surface/matrix properties, matching the ThermaKin convention; a fully depleted
  cell falls back to the all-phase occupancy so k_eff degrades to the gas-mixture
  conductivity instead of 0. **Gas transfer uses the all-phase volume fractions**
  (gas included) — diffusion happens through the pore space the matrix basis
  excludes. Component-based permeability mixes on the matrix basis. Absorption uses
  concentration weighting; mixture density uses mass fractions with the swelling
  factor.

---

## 3.6 Assembling a `Material`

### 3.6.1 The constructor

```julia
Material(; name::Symbol,
           components::Tuple,
           reactions::Tuple = (),
           mixing = MaterialMixing(),
           lateral_shrinkage_law = nothing,
           depletion_limiter = DepletionLimiter())
```

| Keyword | Meaning | Default |
|---|---|---|
| `name` | material identifier (`Symbol`) | required |
| `components` | tuple of components | required |
| `reactions` | tuple of reactions (may be empty `()`) | `()` |
| `mixing` | `(k, λ, κ)` NamedTuple from `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `A(t)/A_0 = law(χ̄)` or `nothing` (identity) | `nothing` |
| `depletion_limiter` | `DepletionLimiter` applied to every reaction-rate evaluation (smooth `tanh(ξ/threshold)` roll-off near reactant depletion; §4.10) | `DepletionLimiter()` (`threshold = 1.0` kg·m⁻³, enabled) |

- **Component order defines the index space.** The first component is index 1, the
  second is 2, etc. The state vector's species blocks follow this order, and integer
  reaction references (§3.4.2) use it. Component names must be unique.
- **Symbolic reactions are resolved here.** Any `PendingReaction` built with symbol
  references is matched against the `components` tuple's name→index map. An unknown
  symbol raises a clear error naming the offending reaction and component.
- **`lateral_shrinkage_law`** is the variable-cross-section (intumescence/shrink)
  hook; it must be a `Function` of the column-average pyrolysis progress `χ̄`, or
  `nothing`. Most simulations leave it `nothing`. See Technical Reference §10.5 and
  Chapter 7.

### 3.6.2 Worked example: charring wood with moisture

A wood with virgin solid, char solid, liquid moisture, a dry pyrolysis gas, and water
vapour. The solids carry permeabilities, so the material will run in Darcy–Fick mode;
it uses `WEIGHTED` conductivity and `CARMAN_KOZENY` permeability, with porosity set
on the solids via `φ`.

```julia
using Pyrolysis

# --- Components -----------------------------------------------------------
virgin = SolidComponent(:virgin,
    ρ = 530.0, c = (1500.0, 1.0), k = 0.12,
    ε = 0.88, λ = 2e-5,
    κ = 1e-15,           # tight virgin matrix (trips Darcy mode)
    φ = 0.55)            # pure-phase porosity → intrinsic density set internally

char = SolidComponent(:char,
    ρ = 150.0, c = 1000.0, k = 0.10,
    ε = 0.95, α = 5.0, λ = 2e-5,
    κ = 1e-12,           # open, porous char
    φ = 0.85)

moisture = LiquidComponent(:moisture, ρ = 1000.0, c = 4180.0, k = 0.60)

gas   = GaseousComponent(:gas, M = 0.044, c = 1800.0, k = 0.03, λ = 1e-4)
vapor = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025, λ = 2e-5)

# --- Reactions ------------------------------------------------------------
pyrolysis = Reaction(:pyrolysis, :virgin => (:char, :gas),
                     A = 1e10, E = 140e3, h = 500e3,     # endothermic
                     yields = (0.25, 0.75), T_min = 450.0)

# Moisture evaporation with state-dependent heat of sorption
heat_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[3] / max(ξ[1] + ξ[2], 1.0)
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3
end)
drying = Reaction(:drying, :moisture => :H2O,
                  A = 5.6e8, E = 43e3, h = heat_sorption)

# --- Material -------------------------------------------------------------
wood = Material(
    name = :CharringWood,
    components = (virgin, char, moisture, gas, vapor),   # indices 1..5
    reactions  = (pyrolysis, drying),
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),
        λ = MixingSpec(PARALLEL),
        κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),
    ),
)

# --- Inspect / validate ---------------------------------------------------
# All reached via the Materials submodule (none are top-level Pyrolysis exports).
# n_components / n_reactions / has_darcy_flow / dry_pyrolysis_gas_indices are
# exported from Materials; validate_material is an internal helper — both kinds work
# through the qualified `Pyrolysis.Materials.` path.
Pyrolysis.Materials.n_components(wood)              # 5
Pyrolysis.Materials.n_reactions(wood)               # 2
Pyrolysis.Materials.has_darcy_flow(wood)            # true (virgin/char κ > 0)
Pyrolysis.Materials.dry_pyrolysis_gas_indices(wood) # (4,)  — :H2O (5) excluded
Pyrolysis.Materials.validate_material(wood)         # true (internal helper)
```

### 3.6.3 Worked example: a simple polymer (PMMA)

A non-charring polymer needs only a virgin solid and a single gas, one decomposition
reaction, and no permeability (so transport stays Fickian). Note the positive `h`:
PMMA depolymerization is endothermic (≈ 870 kJ·kg⁻¹ absorbed), so under the
`h > 0 = endothermic` convention it takes a positive value.

```julia
using Pyrolysis

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86),
        GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:depolymerization, :virgin => :MMA,
                 A = 8.5e12, E = 188e3, h = 870e3),   # h > 0 = endothermic
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)

Pyrolysis.Materials.validate_material(pmma)   # true (internal helper, qualified)
```

With `κ = 0` on every component, `has_darcy_flow(pmma)` is `false` and gas leaves by
Fickian diffusion only — the default transport mode. To make MMA absorb in depth you
would set `α` on the solid and choose `radiation_model = BEER_LAMBERT` at solve time.

---

## 3.7 Validating a material

`validate_material` is an **internal** helper (not exported; reach it as
`Pyrolysis.Materials.validate_material`). It is the one call to make if you want an
explicit material check before building a problem. It checks that:

- the material has at least one component,
- every reaction's reactant/product indices are in range `1:NC`, and
- every reaction is mass-balanced (via `verify_stoichiometry`).

```julia
Pyrolysis.Materials.validate_material(wood)   # returns true; warns on problems
```

Mass-balance failures are reported as **warnings** (the function still returns
`true`); index errors **throw**. For a standalone balance check use
`Pyrolysis.Materials.check_material_mass_balance(material.reactions)`. Note that the
two-product `Reaction` constructors already enforce `Σ yields = 1.0` at construction
time (unless `validate_mass = false`), so a material built normally will pass.

Because these helpers are private, the **public** way to validate is to build a
`PyrolysisProblem` and let `solve` validate by default (`validate = true` on
`solve`); the problem-level validation runs the same stoichiometry and index checks.
Calling the internal `validate_material` as you build the material is still useful for
the clearest, earliest error messages during development.

---

## 3.8 Quick reference

**Component defaults** (the values you get if you omit the keyword):

| Field | Solid | Liquid | Gas |
|---|---|---|---|
| `ε` (emissivity) | `0.9` | `0.95` | `0.0` |
| `α` (absorption) | `0.0` | `0.0` | `0.0` |
| `γ` (swelling) | `1.0` | `0.0` | `0.0` |
| `λ` (gas transfer) | `0.0` | `2e-5` | required |
| `κ` (permeability) | `0.0` | `0.0` | n/a (no field) |
| `ρ_intrinsic` | `= ρ` | `= ρ` | n/a |
| density | required `ρ` | required `ρ` | ideal gas from `M` |

**Reaction defaults:** `n = 1.0`, `T_min = 0.0`, `T_max = Inf`,
`validate_mass = true`; single-product yield fixed at `1.0`, two-product `yields`
must sum to `1.0`. Sign: `h > 0` endothermic.

**Mixing defaults:** `k = WEIGHTED (β=0.5)`, `λ = PARALLEL (β=0.5)`,
`κ = WEIGHTED (β=0.5)`. `BRUGGEMAN` is k-only; `CARMAN_KOZENY` is κ-only and requires
`length_scale`.

**Physical constants used by the constructors:** `R_g = 8.31446261815324
J·mol⁻¹·K⁻¹`, `P_ref = 101325 Pa`, `T_ref = 298.15 K`, `σ = 5.670374419e-8
W·m⁻²·K⁻⁴` (all exported as `R_GAS`, `P_REF`, `T_REF`, `STEFAN_BOLTZMANN`).

**Where to go next:** reactions in depth and TGA parameter recovery → Chapter 4;
boundary and initial conditions → Chapter 5; meshes → Chapter 6; building and
running the problem → Chapters 7–8. The physics behind everything here is in
Technical Reference §4 (materials & property functions), §5 (effective properties &
mixing), §6 (kinetics), and §10 (volume change & swelling).


---

# 4. Defining Reactions

Reactions tell Pyrolysis.jl how the solid (and any liquid) phases of a material
decompose into char, intermediates, and gases as the temperature rises. Each
reaction is an Arrhenius rate law with a stoichiometry, a heat of reaction, and
an optional temperature window. This chapter shows how to build reactions with
the `Reaction` constructor, how to specify products and yields, how the heat-of-
reaction sign convention works, how to model multi-step schemes (the single-
reactant restriction), how the depletion limiter affects rates near zero
concentration, and how to validate your stoichiometry. The physics and
derivations behind every formula here are in Technical Reference §6 (Reaction
Kinetics and Stoichiometry).

In this chapter the index `i` denotes a reaction and `j` denotes a component, per
the kinetics convention in the nomenclature.

---

## 4.1 The reaction rate law

Every reaction in Pyrolysis.jl evaluates the same volumetric rate (kg·m⁻³·s⁻¹):

```
r_i = A_i · exp(−E_i / (R_g · T)) · Π_j ξ_j^{n_{i,j}} · tanh(ξ_j/ξ_th) · gate(T)
```

where

- `A_i` is the Arrhenius pre-exponential factor (units depend on order; see §4.4),
- `E_i` is the activation energy (J·mol⁻¹),
- `R_g = 8.314462618` J·mol⁻¹·K⁻¹ is the universal gas constant,
- `T` is the absolute temperature (K),
- `ξ_j` is the mass concentration of reactant `j` (kg·m⁻³) — note
  **mass concentration**, not mole fraction or mass fraction,
- `n_{i,j}` is the reaction order with respect to reactant `j` (must be ≥ 0;
  negative orders are rejected at construction),
- `tanh(ξ_j/ξ_th)` is the smooth depletion roll-off of §4.10, and
- `gate(T)` is a smooth temperature window built from `T_min` / `T_max` (§4.6).

The rate is on a **mass-concentration** basis, which is what makes the species
sources directly additive into the finite-volume residual. See Technical
Reference §6.1 for the derivation.

---

## 4.2 The `Reaction` constructor

`Reaction` is exported by `using Pyrolysis`. There are four call forms,
distinguished by whether you name components by **symbol** or by **integer index**,
and whether the reaction has **one product** or **two products**.

### Symbolic vs. numeric component references

You can refer to components either by the `Symbol` you gave them in their
constructor, or by their 1-based position in the material's `components` tuple.

```julia
using Pyrolysis

# Symbolic — components are referenced by name; resolved when the Material is built.
Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)

# Numeric — :virgin is component 1, :MMA is component 2 in the components tuple.
Reaction(:decomposition, 1 => 2; A = 8.5e12, E = 188e3, h = 870e3)
```

The symbolic form returns a lightweight `PendingReaction`; `Material(...)`
resolves the names to indices automatically using the `components` tuple, so you
never construct a `PendingReaction` yourself. If a reaction references a name that
is not in the components tuple, `Material(...)` raises a clear error. Prefer the
**symbolic** form: it is self-documenting and immune to component reordering.

### Single-product form

```
Reaction(name, reactant => product; A, E, h, n=1.0, T_min=0.0, T_max=Inf, validate_mass=true)
```

The product receives the full mass yield of 1.0 automatically; you do not pass
`yields`.

```julia
# virgin PMMA fully volatilizes to MMA gas, endothermic.
Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3, n = 1.0)
```

### Two-product form

```
Reaction(name, reactant => (product1, product2); A, E, h, yields, n=1.0, T_min=0.0, T_max=Inf, validate_mass=true)
```

Here `yields = (y1, y2)` gives the mass split between the two products and **must
sum to 1.0** (see §4.3).

```julia
# Charring: 1 kg virgin → 0.2 kg char + 0.8 kg gas, endothermic.
Reaction(:charring, :virgin => (:char, :gas); A = 1e15, E = 200e3, h = 500e3,
         yields = (0.2, 0.8))
```

### Keyword arguments and defaults

| Keyword | Meaning | Units | Default |
|---|---|---|---|
| `A` | Arrhenius pre-exponential factor (required) | s⁻¹ (first order) | — |
| `E` | Activation energy (required) | J·mol⁻¹ | — |
| `h` | Heat of reaction (required; number, `(a,b)` tuple, or property function) | J·kg⁻¹ | — |
| `yields` | Mass split for the two-product form (required there) | – | — |
| `n` | Reaction order w.r.t. the reactant | – | `1.0` |
| `T_min` | Lower temperature gate | K | `0.0` |
| `T_max` | Upper temperature gate | K | `Inf` |
| `validate_mass` | Enforce mass balance at construction | – | `true` |

`A`, `E`, and `h` are **required** keyword arguments — there are no defaults, and
omitting any of them is an error. `yields` is required only for the two-product
form. All arguments after the component pair are keyword arguments (note the `;`
in the signatures above), though the examples also work with `,` because Julia
accepts trailing keyword arguments either way.

---

## 4.3 Products, yields, and mass balance

Pyrolysis.jl uses a **single-reactant** convention (see §4.7). Per unit reaction
extent, the one reactant is consumed with mass stoichiometry `−1.0` (this is
implicit — you never type it), and the products are generated with the mass
stoichiometries you supply. Mass is conserved when the product yields sum to one:

```
Σ_products ν_{i,j} = 1.0
```

For the single-product form, the lone product automatically gets yield `1.0`, so
mass balance is trivially satisfied. For the two-product form, your `yields`
tuple must sum to `1.0`:

```julia
yields = (0.2, 0.8)   # 0.2 + 0.8 = 1.0  ✓  (mass conserved)
yields = (0.2, 0.7)   # 0.2 + 0.7 = 0.9  ✗  (validate_mass throws an error)
```

When `validate_mass = true` (the default), the constructor runs an internal
stoichiometry check and **throws an error at construction time** if the masses do
not balance. This catches typos before you ever solve. You can pass
`validate_mass = false` to skip the check (not recommended — it lets mass
conservation silently break).

The species source assembled per cell is the stoichiometric sum
`S_j = Σ_i ν_{i,j} r_i` (kg·m⁻³·s⁻¹), with the reactant carrying `−r_i` and each
product carrying `+yield · r_i`. See Technical Reference §6.5.

---

## 4.4 The Arrhenius parameters `A` and `E`

- `A` (pre-exponential factor) and `E` (activation energy, J·mol⁻¹) are stored as
  plain scalars (`Float64`). They are **temperature-independent**: a
  temperature-dependent `A(T)` is not supported directly (see §4.14).
- The units of `A` depend on the reaction order. For the standard first-order
  case (`n = 1.0`) the rate constant `A·exp(−E/R_gT)` has units s⁻¹ and the rate
  `r` comes out in kg·m⁻³·s⁻¹ because it is multiplied by a concentration. For a
  bimolecular convention `A` would be m³·kg⁻¹·s⁻¹. **The constructor does not
  check these units** — getting them right is your responsibility.

```julia
Reaction(:pyrolysis, :virgin => (:char, :gas);
         A = 1.0e13,    # s⁻¹
         E = 150e3,     # J/mol
         h = 300e3,     # J/kg
         yields = (0.20, 0.80))
```

§4.9 gives a practical recipe for recovering `A` and `E` from a TGA reference
temperature and peak rate, FDS-style.

---

## 4.5 Reaction order `n`

The keyword `n` sets the order with respect to the reactant concentration: the
rate carries the factor `ξ_reactant^n`. The default is first order (`n = 1.0`),
which is the usual choice for solid-state thermal decomposition. Non-integer
orders are allowed.

```julia
# Half-order decomposition
Reaction(:decomp, :virgin => :gas; A = 1e12, E = 180e3, h = 500e3, n = 0.5)
```

The rate evaluation has a fast path for orders 0 and 1 (avoiding a general
power), so first-order reactions are both the most common and the cheapest.

---

## 4.6 Temperature gates `T_min` and `T_max`

`T_min` and `T_max` define the temperature window in which the reaction is active.
They are **not hard cutoffs**: instead of step functions, Pyrolysis.jl multiplies
the rate by smooth `tanh` ramps,

```
r ← r · ½(1 + tanh(T − T_min)) · ½(1 + tanh(T_max − T))
```

so the reaction turns on and off smoothly. The smoothness is deliberate — hard
steps create discontinuities in the Jacobian that break automatic differentiation
and Newton convergence (Technical Reference §6.2).

```julia
Reaction(:pyrolysis, :virgin => (:char, :gas);
         A = 1.0e13, E = 150e3, h = 300e3, yields = (0.20, 0.80),
         T_min = 400.0,    # ramps on near 400 K
         T_max = 1500.0)   # ramps off near 1500 K
```

Defaults are `T_min = 0.0` and `T_max = Inf`. The lower gate is applied only when
`T_min` is finite **and** strictly positive, and the upper gate only when `T_max`
is finite; with the defaults, neither gate fires and the bare Arrhenius rate is
used at all temperatures.

> **Caveat — gate leakage.** Because the ramp uses `tanh(T − T_min)` (argument in
> kelvin), the transition is sharp on the kelvin scale but not instantaneous: the
> reaction is "weakly active" within roughly ±10 K of each bound. If you need a
> hard kinetic threshold, the gates are an approximation; for the AD- and
> Newton-stability reasons above, they are intentional. See Technical Reference
> §6.2.

---

## 4.7 The single-reactant restriction and multi-step schemes

The current implementation supports **single-reactant reactions only**
(`NRin == 1`). Each `Reaction` consumes exactly one component. Attempting to apply
a multi-reactant reaction raises an error at evaluation. This means you cannot
write a bimolecular reaction such as `char + O₂ → ...` as a single `Reaction` with
two reactants.

The way to model realistic decomposition is a **sequence of single-reactant
reactions**, each consuming one component and producing intermediates plus gases.
This is exactly how the shipped wood and pressure-treated-wood models work: a
chain `virgin → int1 → int2 → ... → char`, with a gas split released at every
step. Below is a three-step charring chain (symbolic names; numeric indices would
read identically):

```julia
using Pyrolysis

virgin = SolidComponent(:virgin, ρ = 530.0, c = 1500.0, k = 0.12)
int1   = SolidComponent(:int1,   ρ = 400.0, c = 1400.0, k = 0.11)
int2   = SolidComponent(:int2,   ρ = 250.0, c = 1300.0, k = 0.12)
char   = SolidComponent(:char,   ρ = 130.0, c = 1100.0, k = 0.10, ε = 0.95)
gas1   = GaseousComponent(:gas1, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
gas2   = GaseousComponent(:gas2, M = 0.040, c = 1500.0, k = 0.03, λ = 1e-5)
gas3   = GaseousComponent(:gas3, M = 0.050, c = 1500.0, k = 0.03, λ = 1e-5)

material = Material(
    name = :CharringWood,
    components = (virgin, int1, int2, char, gas1, gas2, gas3),
    reactions = (
        Reaction(:step1, :virgin => (:int1, :gas1);
                 A = 7.48e4,  E = 71.6e3,  h = -1.92e4, yields = (0.95, 0.05)),
        Reaction(:step2, :int1   => (:int2, :gas2);
                 A = 2.36e10, E = 140.9e3, h =  2.64e5, yields = (0.80, 0.20)),
        Reaction(:step3, :int2   => (:char, :gas3);
                 A = 6.40e0,  E = 56.0e3,  h = -3.3e2,  yields = (0.74, 0.26)),
    ),
)
```

Each step is mass-balanced on its own (every `yields` sums to 1.0), and the chain
as a whole conserves solid + gas mass.

> **Note.** O₂-dependent char oxidation (an exothermic surface reaction common in
> flaming-combustion models) cannot currently be expressed as a true bimolecular
> reaction. You can approximate it with a single-reactant exothermic char-→-ash
> step gated by temperature (§4.8), but the explicit oxygen dependence is a
> documented limitation (Technical Reference §6.10).

---

## 4.8 Heat of reaction `h` and its sign convention

The heat of reaction is supplied through the `h` keyword in J·kg⁻¹ (per kg of the
reactant). Pyrolysis.jl uses the **storage (enthalpy) sign convention**:

- **`h > 0` is endothermic** — the reaction absorbs heat and **cools** the
  material.
- **`h < 0` is exothermic** — the reaction releases heat and **heats** the
  material.

The volumetric heat source assembled into the energy equation is
`Q_rxn = −Σ_i h_i r_i` (W·m⁻³), so the minus sign converts the stored enthalpy
into a source: an endothermic reaction (`h > 0`) gives `Q_rxn < 0` (cooling), and
an exothermic reaction (`h < 0`) gives `Q_rxn > 0` (heating). See Technical
Reference §6.6.

> **Important — opposite to ThermaKin / Gpyro.** ThermaKin2Ds and Gpyro publish
> heats of reaction with the *opposite* sign (in those codes `h > 0` is
> exothermic). When you port parameters from a ThermaKin or Gpyro input file,
> **flip the sign of every heat of reaction.** The shipped examples annotate this
> conversion explicitly; e.g. a ThermaKin "+2.78E6 J/kg" (exothermic there)
> becomes `h = -2.78e6` here (exothermic in this convention) — and a ThermaKin
> "−2.78E6 J/kg" (endothermic there) becomes `h = 2.78e6` here (endothermic in
> this convention). In both cases the rule is the same: flip the sign.

### Constant heat (the common case)

Pass a plain number; it is wrapped as a `ConstantProperty`:

```julia
# Endothermic pyrolysis (absorbs heat).
Reaction(:decomposition, :virgin => :gas; A = 1.5e14, E = 2.03e5, h = 820e3)

# Exothermic char-forming step (releases heat) — note the negative h.
Reaction(:exo_step, :int4 => (:char, :gas);
         A = 6.4e0, E = 56e3, h = -3.26e2, yields = (0.74, 0.26))
```

### Temperature-dependent heat

Pass a two-tuple `(a, b)` for an affine `h(T) = a + b·T` (a `LinearProperty`), or
any of the property-function types from Chapter 3:

```julia
# h(T) = 870e3 + 100*T  J/kg  (h > 0 = endothermic)
Reaction(:decomp, :virgin => :gas; A = 8.5e12, E = 188e3, h = (870e3, 100.0))
```

### State-dependent heat (moisture sorption)

For moisture evaporation, the heat of sorption depends on the current moisture
content, not only on temperature. Use a `StateDependentProperty`, whose closure
receives both the temperature `T` and the full concentration tuple `ξ`. The
kinetics evaluator automatically passes `ξ` when the heat is state-dependent.

```julia
const MOISTURE_IDX = 7        # index of moisture in the components tuple
const ρ_DRY_BASIS  = 530.0    # dry-wood density for the moisture-content ratio [kg/m³]

# HoS(MC) = 2450.4 + 875.6 * exp(-0.1122 * MC) [kJ/kg], MC = moisture content (% dry basis)
heat_of_sorption = StateDependentProperty(
    (T, ξ) -> begin
        MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS          # moisture content, % dry basis
        return (2450.420909987441 + 875.6441221735336 *
                exp(-0.1121684273316275 * MC)) * 1e3          # J/kg (endothermic, h > 0)
    end,
)

Reaction(:moisture_evaporation, :moisture => :water_vapor;
         A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0)
```

The closure indexes `ξ` by component position, so `MOISTURE_IDX` must match the
moisture component's slot in the `components` tuple. Evaporation is endothermic,
so the returned value is positive. See Technical Reference §6.6 and Chapter 3 for
the property-function types.

---

## 4.9 Recovering `A` and `E` from TGA data (FDS-style recipe)

A common practical task is fitting Arrhenius parameters to a single-step
thermogravimetric analysis (TGA) curve. The FDS approach uses two readouts from a
constant-heating-rate TGA run: the **peak-rate temperature** `T_p` (K, where the
mass-loss rate is maximal) and the **peak reaction rate** `r_p` (s⁻¹, the maximum
fractional mass-loss rate divided by the heating rate `β_h` in K·s⁻¹). For a
first-order reaction the standard estimates are

```
E ≈ e · r_p · R_g · T_p²
A ≈ e · r_p · exp(E / (R_g · T_p))
```

where `e ≈ 2.71828` is Euler's number and `R_g = 8.314462618` J·mol⁻¹·K⁻¹. As a
runnable example, suppose the TGA peak is at `T_p = 640 K` with a heating rate of
10 K·min⁻¹ and a peak fractional mass-loss rate of 0.012 K⁻¹:

```julia
const R_g = 8.314462618    # J/(mol·K)

β_h = 10.0 / 60.0          # heating rate, K/s  (10 K/min)
T_p = 640.0                # peak-rate temperature, K
r_p = 0.012 * β_h          # peak reaction rate, s⁻¹  (0.012 K⁻¹ × β_h)

E = exp(1) * r_p * R_g * T_p^2          # J/mol
A = exp(1) * r_p * exp(E / (R_g * T_p)) # s⁻¹

rxn = Reaction(:decomposition, :virgin => :gas; A = A, E = E, h = 500e3)
```

These are first-pass estimates; for production fits, refine `A`, `E`, and `h`
against the measured curve using the sensitivity / inverse-analysis tooling (User
Guide Chapter 13). See Technical Reference §6.9 for the relationship to the FDS,
ThermaKin, and Gpyro kinetic formulations.

---

## 4.10 The depletion limiter

As a reactant concentration approaches zero, a bare Arrhenius rate can still try
to react at a finite speed, which makes implicit solvers stiff — and for a
zero-order fit it would keep consuming reactant straight through zero.
Pyrolysis.jl rolls the rate off smoothly by multiplying it with a `tanh` factor
per reactant:

```
r *= tanh(ξ / threshold)
```

with negative concentrations clamped to zero first. The factor is smooth
everywhere (no threshold branch): it is ≈ 1 well above `threshold` (within 4%
at `2·threshold`, indistinguishable from 1 beyond `~20·threshold`) and drives
the rate to zero with a bounded derivative as the reactant runs out — for every
reaction order `n ≥ 0`, including zero-order. The limiter is configured by the
`DepletionLimiter` type, exported by the `Physics` submodule
(`Pyrolysis.Physics.DepletionLimiter`):

```julia
Base.@kwdef struct DepletionLimiter
    threshold::Float64 = 1.0    # kg/m³ — rate roll-off scale
    enabled::Bool      = true
end
```

| Field | Meaning | Units | Default |
|---|---|---|---|
| `threshold` | Concentration scale of the rate roll-off | kg·m⁻³ | `1.0` |
| `enabled` | Turn limiting on/off | – | `true` |

The default limiter (`threshold = 1.0`, `enabled = true`) is applied
automatically during a solve. The limiter is a property of the **material**, not
of individual reactions or of `solve`: pass a custom one with the `Material`
constructor's `depletion_limiter` keyword,

```julia
material = Material(name = …, components = …, reactions = …,
                    depletion_limiter = DepletionLimiter(threshold = 0.1))
```

and the solver's kinetics path reads `material.depletion_limiter` for every rate
evaluation. **When to tune:** the default `threshold = 1.0` kg·m⁻³ is appropriate for
most condensed-phase materials whose concentrations are hundreds to thousands of
kg·m⁻³. If you work with very dilute species (e.g. trace moisture at a few
kg·m⁻³), the default threshold may damp them more than intended; conversely, if a
solve is stiff or produces small negative concentrations near full conversion, a
slightly larger threshold can help. See Technical Reference §6.3.

You can evaluate a single rate with a custom or disabled limiter for testing:

```julia
using Pyrolysis
using Pyrolysis.Physics: reaction_rate, DepletionLimiter   # exported by Physics, not re-exported at the top level

# Disable limiting (validation / verification only — may be stiff in a full
# solve, and removes the only guard that stops zero-order reactions at ξ = 0).
r_raw = reaction_rate(rxn, ξ_cell, T_cell; limiter = DepletionLimiter(enabled = false))
```

where `ξ_cell` is the concentration `NTuple` and `T_cell` the temperature.

---

## 4.11 Assembling reactions into a material

Reactions are passed to `Material(...)` as a tuple via the `reactions` keyword.
The components tuple comes first (so symbolic names can be resolved), and the
reaction tuple may mix heterogeneous reaction types (different orders, product
counts, heat-property types) without penalty:

```julia
using Pyrolysis

material = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86),
        SolidComponent(:char,   ρ = 5.0,  c = 1000, k = 0.10, ε = 0.95),
        GaseousComponent(:MMA,  M = 0.10012, c = 1500, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => (:char, :MMA);
                 A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0,
                 yields = (0.002, 0.998)),   # 0.2% char, 99.8% gas
    ),
)
```

A material can also have **no reactions** — pass `reactions = ()` (the default) for
a pure heat-conduction problem. `Material(...)` resolves symbolic reactions,
validates that every reaction references a real component, and checks each
reaction's mass balance at construction. Chapter 3 covers the rest of the
`Material` constructor (components, mixing rules, lateral-shrinkage law).

### Building a reaction set separately

If you want to construct the reaction collection on its own (for reuse or
inspection), `ReactionSet` wraps a tuple of reactions and is exported:

```julia
rxns = ReactionSet(
    Reaction(:step1, :virgin => (:int1, :gas1);
             A = 7.48e4,  E = 71.6e3,  h = -1.92e4, yields = (0.95, 0.05)),
    Reaction(:step2, :int1   => (:char, :gas2);
             A = 2.36e10, E = 140.9e3, h =  2.64e5, yields = (0.80, 0.20)),
)
```

In normal use you pass the plain tuple straight to
`Material(..., reactions = (...))`; `ReactionSet` is mainly an accessor
convenience.

---

## 4.12 Validating reactions

Two helpers let you check mass balance explicitly, beyond the automatic
check at construction. They are **not** part of the top-level `Pyrolysis`
exports — reach them through their submodule as
`Pyrolysis.Materials.verify_stoichiometry` and
`Pyrolysis.Materials.check_material_mass_balance`, or bring them into scope
with `using Pyrolysis.Materials: verify_stoichiometry, check_material_mass_balance`.

**`verify_stoichiometry(rxn; tol=1e-10) -> Bool`** confirms that a single
reaction's product yields sum to 1.0 (within `tol`). It returns `true` when
balanced and issues a warning and returns `false` otherwise:

```julia
using Pyrolysis.Materials: verify_stoichiometry, check_material_mass_balance

rxn = Reaction(:charring, :virgin => (:char, :gas);
               A = 1e15, E = 200e3, h = 500e3, yields = (0.25, 0.75))
verify_stoichiometry(rxn)   # true: 0.25 + 0.75 = 1.0
```

**`check_material_mass_balance(reactions; tol=1e-10) -> Bool`** runs
`verify_stoichiometry` over every reaction in a material and returns `true` only
if all pass:

```julia
if !check_material_mass_balance(material.reactions)
    @warn "Material has reactions with mass-balance issues"
end
```

Both are diagnostics: mass conservation is enforced by the single-reactant
construction convention (§4.3), not by a projection during the solve, so these
checks are for catching setup mistakes. See Technical Reference §6.8.

---

## 4.13 Worked example: a single-step endothermic decomposition

Putting it together — a one-reaction PMMA-like material, evaluated at a single
cell state to show the rate, heat source, and species sources:

```julia
using Pyrolysis
# `reaction_rate` lives in the Physics submodule and is exported from it, but is
# not re-exported at the top level — bring it into scope explicitly.
using Pyrolysis.Physics: reaction_rate
# `verify_stoichiometry` / `check_material_mass_balance` are not exported at all;
# reach them through the Materials submodule.
using Pyrolysis.Materials: verify_stoichiometry, check_material_mass_balance

# Three components: virgin solid, char residue, gas product.
virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
char   = SolidComponent(:char,   ρ = 5.0,  c = 1000, k = 0.10, ε = 0.95)
gas    = GaseousComponent(:MMA,  M = 0.10012, c = 1500, k = 0.02, λ = 1e-5)

rxn = Reaction(:decomposition, :virgin => (:char, :MMA);
               A = 1.5e14,     # s⁻¹
               E = 2.03e5,     # J/mol
               h = 820e3,      # J/kg  (positive ⇒ endothermic ⇒ cools)
               n = 1.0,
               yields = (0.002, 0.998))

material = Material(name = :PMMA, components = (virgin, char, gas), reactions = (rxn,))

# Verify mass balance.
@assert verify_stoichiometry(rxn)                  # 0.002 + 0.998 = 1.0
@assert check_material_mass_balance(material.reactions)

# Evaluate the rate at one cell state (concentrations in kg/m³, T in K).
T_cell = 700.0
ξ_cell = (1000.0, 0.0, 0.0)                         # (virgin, char, MMA)
r = reaction_rate(rxn, ξ_cell, T_cell)             # kg·m⁻³·s⁻¹

Q_rxn = -rxn.heat(T_cell, ξ_cell) * r              # W·m⁻³  (= -h·r); endothermic ⇒ < 0
S = (-r, 0.002 * r, 0.998 * r)                     # species sources (virgin, char, gas)
@assert isapprox(S[2] + S[3], r; rtol = 1e-12)     # product mass = consumed mass
```

`r` is the volumetric reaction rate, `Q_rxn = −h·r` is the volumetric heat source
(negative here because the reaction is endothermic), and the species sources `S`
conserve mass (the products' generation balances the reactant's consumption).
This is exactly the per-cell computation the solver performs internally; see
Technical Reference §6.7.

---

## 4.14 Summary and limitations

- A `Reaction` is `A·exp(−E/R_gT)·ξ^n` with smooth `T_min`/`T_max` gates and a
  heat of reaction `h`.
- Reference components by **symbol** (recommended) or integer index; product
  `yields` must sum to **1.0**, enforced at construction unless
  `validate_mass = false`.
- **Sign convention: `h > 0` is endothermic, `h < 0` is exothermic** — opposite
  to ThermaKin/Gpyro; flip the sign when porting parameters.
- Only **single-reactant** reactions are supported; build multi-step schemes as
  chains of single-reactant reactions. Multi-reactant (e.g. O₂-dependent)
  kinetics, temperature-dependent `A`, and pressure-dependent rates are not
  supported (Technical Reference §6.10).
- The **depletion limiter** (default `threshold = 1.0` kg·m⁻³) is applied
  automatically and is rarely tuned; lower it only for genuinely dilute species
  and raise it only to relieve stiffness near full conversion.
- Validate with `verify_stoichiometry` and `check_material_mass_balance`.

For the mathematics behind the rate law, the smooth gates, the depletion limiter,
the stoichiometric sources, and the heat-source assembly, see Technical Reference
§6 (Reaction Kinetics and Stoichiometry). For defining the components a reaction
acts on, see Chapter 3 (Defining Materials, Components, and Properties).


---

# 5. Boundary and Initial Conditions

Boundary conditions (BCs) tell the solver how heat and gas enter and leave the
1-D column; initial conditions set the temperature and composition at `t = 0`.
This chapter shows how to build a `BoundaryConditionSet`, the complete catalogue
of thermal, mass, pressure, and mesh-motion BCs with every field and default,
how the surface energy balance is solved, and how to specify initial
temperature and concentration profiles. The physics and derivations behind each
BC live in the Technical Reference §12 (and the radiation models in §8); this
chapter is about *how to specify them correctly*.

> **Coordinate and sign reminders** (see Technical Reference §2). The axial
> coordinate is `z`. `z = 0` is the **bottom / substrate** (tag `:bottom`,
> boundary id 2); `z = L` is the **top / exposed surface** (tag `:top`, boundary
> id 1). Heat normally enters from the top. For *all* BCs in this package the
> sign convention is: **positive flux = into the domain** (heating, or mass
> entering). A gas escaping out of the top therefore has a *negative* boundary
> mass flux in the BC's own convention.

---

## 5.1 The `BoundaryConditionSet`

All boundary conditions are bundled in a single tag-keyed container,
`BoundaryConditionSet`. It holds exactly two kinds of BC — **thermal** and
**mass** — each as a `NamedTuple` keyed by boundary tag (`:top`, `:bottom` in
1-D):

```julia
using Pyrolysis

bc_set = BoundaryConditionSet(
    thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top    = DiffusiveMassBC(0.01),
               bottom = ImpermeableBC()),
)
```

Both keyword arguments default to fully-insulated, impermeable boundaries:

```julia
# Equivalent to the explicit defaults below
bc_set = BoundaryConditionSet()
#   thermal = (top = AdiabaticBC(),   bottom = AdiabaticBC())
#   mass    = (top = ImpermeableBC(), bottom = ImpermeableBC())
```

Key facts about the constructor:

- **Partial specification is allowed.** Whatever you supply is merged over the
  defaults, so `BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3),))`
  leaves the bottom adiabatic and both mass faces impermeable. (Note the
  trailing comma — a one-element `NamedTuple` needs it.)
- **Type checking is enforced at construction.** Passing a mass BC into a
  thermal slot (or vice-versa) throws an `ArgumentError` immediately.
- **It is type-stable.** Each tag compiles to a direct field read via
  `Val{tag}` dispatch, so the residual specializes statically with zero runtime
  overhead.

The set is then handed to `PyrolysisProblem` via the `bc_set` keyword:

```julia
problem = PyrolysisProblem(
    mesh = mesh,
    material = material,
    bc_set = bc_set,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 100.0),
)
```

> **Note on mesh motion and radiation BCs.** The `BoundaryConditionSet` has
> *only* `thermal` and `mass` fields. ALE mesh motion needs no BC at all — it is
> computed automatically from the physics (§5.7) — and the P1 radiation
> BCs (§5.6) are **not** stored in this container: they live in a separate
> `RadiationBCSet` inside the `Experimental` module. External incident radiation
> (cone heater, panel) is specified as a **thermal** BC (`RadiativeFluxBC`,
> §5.2.6), not as a radiation BC.

---

## 5.2 Thermal boundary conditions

Thermal BCs are subtypes of `AbstractThermalBC`. They are the heat-side input to
the **surface energy balance** that the solver closes by a Newton iteration for
the surface temperature `T_s` (§5.4). Every flux is reported with the **positive
= into domain** convention.

| BC type | Physical flux `q_BC` | Key parameters (defaults) |
|---|---|---|
| `AdiabaticBC` | `0` | — |
| `HeatFluxBC` | `q(t)` (prescribed) | `flux` |
| `ConvectiveBC` | `h·(T_∞ − T_s)` | `h`, `T_inf` |
| `RadiativeBC` | `F·ε_eff·σ·(T_∞⁴ − T_s⁴)` | `ε` (`nothing`), `T_inf`, `F` (`1.0`) |
| `RadiativeFluxBC` | `α·Φ_incident` (model-dependent) | `flux`, `absorptivity` (`nothing`) |
| `DirichletThermalBC` | enforces `T_s = T` | `T` |
| `CombinedThermalBC` | `Σ_i q_BC,i` | built via `+` |

All exported by `using Pyrolysis`.

### 5.2.1 `AdiabaticBC` — perfect insulation

```julia
bc = AdiabaticBC()
```

Zero heat flux. This is the default for both faces and the usual choice for the
substrate (back face) in a cone-calorimeter idealization where the sample is
mounted on insulation.

### 5.2.2 `HeatFluxBC` — prescribed flux

```julia
bc = HeatFluxBC(50e3)                            # constant 50 kW/m² into the surface
bc = HeatFluxBC(t -> 50e3 * (1 - exp(-t / 10)))  # ramp-up over ~10 s
```

The single field, `flux`, is the heat flux `q_BC` in W/m², **positive into the
domain**. It may be a constant (`Real`) or a function of time `t`. Use this for
direct-contact heating or for an opaque surface where you want to inject a flux
directly (rather than via the radiation models, see `RadiativeFluxBC`).
`dq_BC/dT_s = 0`, so it does not couple to the surface temperature.

### 5.2.3 `ConvectiveBC` — Newton cooling

```julia
bc = ConvectiveBC(h = 10.0, T_inf = 300.0)
bc = ConvectiveBC(h = 10.0, T_inf = t -> 300 + 0.1 * t)   # warming ambient
```

Newton's cooling law `q_BC = h·(T_∞ − T_s)`. Fields:

- `h` — convective heat-transfer coefficient `h_conv` in W·m⁻²·K⁻¹. Constant or
  `h(t)`. Validated to be non-negative.
- `T_inf` — ambient temperature `T_∞` in K. Constant or `T_∞(t)`. Validated to
  be positive.

The flux is positive (heating) when `T_∞ > T_s` and negative (cooling) when the
surface is hotter than the surroundings — the usual re-radiation/convection loss
direction during pyrolysis. Its temperature derivative `dq_BC/dT_s = −h` makes
the Newton solve converge quickly. The `T_inf` here is also the value used to
resolve the ambient temperature for the surface balance (§5.8).

**Blowing correction (opt-in).** Strong outgassing thickens the boundary layer
and suppresses convective exchange. With `solve(...; use_blowing_correction =
true)` the coefficient is rescaled by the classical film factor
`h → h·B/(e^B − 1)`, `B = ṁ″·c_p,g/h`, using the outward gas mass flux at the
face. At a wood-like peak `ṁ″ ≈ 10 g·m⁻²·s⁻¹` against `h ≈ 10` this is a
35–50 % reduction of `h`; at PMMA-like peaks (`B ≈ 4`) convective cooling is
essentially extinguished. The back face is unaffected automatically (no
outflow → `B = 0`). Off by default because the shipped examples were
calibrated without it (and without flame feedback — the two partially cancel;
see the physics roadmap before recalibrating).

### 5.2.4 `RadiativeBC` — surface re-radiation

```julia
bc = RadiativeBC(T_inf = 300.0)                  # uses material's effective emissivity
bc = RadiativeBC(ε = 0.9, T_inf = 300.0)         # fixed emissivity
bc = RadiativeBC(T_inf = 300.0, F = 0.7)         # partially enclosed surface
bc = RadiativeBC(ε = T -> 0.85 + 0.1*(T-300)/700, T_inf = 300.0)  # T-dependent ε
```

Gray-body radiative exchange with the surroundings,
`q_BC = F·ε_eff·σ·(T_∞⁴ − T_s⁴)`, where `σ = 5.670374419×10⁻⁸` W·m⁻²·K⁻⁴.
Fields:

- `ε` — surface emissivity. Default `nothing`, which means **use the effective
  emissivity** `ε_eff = Σ_j v_j ε_j` (volume-weighted) computed from the live
  composition. This is the recommended setting: it automatically tracks the
  change from virgin polymer to char as the material decomposes. You can also
  pass a constant in `[0, 1]` or a function `ε(T)`.
- `T_inf` — surroundings temperature `T_∞` in K (constant or `T_∞(t)`).
- `F` — view factor in `[0, 1]`, default `1.0` (surface fully sees the
  surroundings). Use `F < 1` for partially enclosed geometries.

Its temperature derivative is `dq_BC/dT_s = −4·F·ε_eff·σ·T_s³`. When `ε` is a
function of temperature, the solver **neglects** the `dε/dT` term in this
derivative; this slightly slows Newton convergence but does not change the
converged root (see Technical Reference §12.2).

### 5.2.5 `DirichletThermalBC` — prescribed temperature

```julia
bc = DirichletThermalBC(500.0)            # hold this face at 500 K
bc = DirichletThermalBC(t -> 300 + 5*t)   # ramped wall temperature
```

Forces the surface to a prescribed temperature `T = T_s`. This BC bypasses the
ordinary surface-balance Newton iteration and is instead enforced with a large
penalty conductance (`h_large = 1e8` W·m⁻²·K⁻¹) plus a special case in the RHS
assembly. Use it for verification problems or when an experiment fixes the
surface temperature directly. The single field `T` may be constant or `T(t)`.

### 5.2.6 `RadiativeFluxBC` — external incident radiation

This is the BC for a **cone heater, furnace, or radiant panel** — an external
radiation source heating the surface. It is distinct from `RadiativeBC` (which
is the surface's *own* re-radiation loss) and from `HeatFluxBC` (which injects a
bare flux with no absorptivity or in-depth treatment).

```julia
bc = RadiativeFluxBC(flux = 25e3)                          # incident 25 kW/m², α = ε_eff
bc = RadiativeFluxBC(flux = 50e3, absorptivity = 0.9)      # fixed absorptivity
bc = RadiativeFluxBC(flux = t -> 25e3*(1 - exp(-t/10)))    # time-varying flux
```

Fields:

- `flux` — incident radiative flux `Φ_incident` in W/m² (constant or `Φ(t)`).
- `absorptivity` — surface absorptivity `α` in `[0, 1]`. Default `nothing`,
  which uses the effective emissivity `ε_eff` as the absorptivity (Kirchhoff's
  law `α = ε_eff`), again tracking composition automatically. You can override
  with a constant or `α(t)`.

**Its behavior depends on the `radiation_model` keyword passed to `solve()`**
(see §5.5 and Technical Reference §8):

| `radiation_model` (in `solve`) | What `RadiativeFluxBC` does |
|---|---|
| `NO_RADIATION` | Absorbed flux `α·Φ` deposited as a surface boundary heat flux |
| `SURFACE_ABSORPTION` | Same — absorbed at the surface (opaque material) |
| `BEER_LAMBERT` | Boundary flux is **zero**; the absorbed `α·Φ` is distributed in-depth volumetrically (`extras.volumetric_absorption`) |
| `P1_QUASI_STEADY` | Boundary flux is **zero**; radiation handled by the P1 solver (Experimental) |

The boundary flux from this BC has `dq_BC/dT_s = 0` (it is the *incoming* source
and does not depend on the surface temperature).

> **Do not double-count radiation.** Do not combine `RadiativeFluxBC` and
> `HeatFluxBC` for the *same* radiation source — that adds the energy twice.
> `RadiativeFluxBC` is the correct, absorptivity-aware choice for a radiant
> source.

### 5.2.7 `CombinedThermalBC` — adding BCs with `+`

A real exposed surface usually has several simultaneous heat exchanges:
incoming radiation, re-radiation loss, and convective loss. Add them with `+`:

```julia
top_thermal =
    RadiativeFluxBC(flux = 50e3, absorptivity = 0.96) +  # cone heater (incoming)
    RadiativeBC(ε = 0.96, T_inf = 330.0) +               # re-radiation loss
    ConvectiveBC(h = 10.0, T_inf = 330.0)                # convective loss
```

The combined flux is the sum of the component fluxes,
`q_BC = Σ_i q_BC,i(T_s, t)`, and the derivative is summed likewise — assembled
by a generated function for type stability. You can chain any number of thermal
BCs. The ambient temperature for the surface balance is taken from the first
component that carries one (preferring `ConvectiveBC`/`RadiativeBC`, see §5.8).

A canonical cone-calorimeter top face:

```julia
bc_set = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.96) +
              RadiativeBC(ε = 0.96, T_inf = 330.0) +
              ConvectiveBC(h = 10.0, T_inf = 330.0),
        bottom = AdiabaticBC(),
    ),
    mass = (
        top    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
        bottom = ImpermeableBC(),
    ),
)
# Pair with:  solve(problem; radiation_model = SURFACE_ABSORPTION)
```

---

## 5.3 Mass boundary conditions

Mass BCs (subtypes of `AbstractMassBC`) control how **gas** escapes (or enters)
through a boundary. Solid and liquid components never transport across the
boundary; only gaseous components do. The sign convention is **positive = mass
into the domain**, so gas escaping at the top has a negative flux in the BC's
own convention (the solver applies the correct face orientation internally).

| BC type | Physical flux `J_j` | Key parameters (defaults) |
|---|---|---|
| `ImpermeableBC` | `0` | — |
| `MassFluxBC` | `J(t)` (prescribed) | `flux` |
| `ConvectiveMassBC` | `h_m·(ξ_∞ − ξ_s)` | `h_m`, `ξ_inf` (`0.0`) |
| `DiffusiveMassBC` | combined internal + external resistance | `h_m` (`Inf`) |
| `CombinedMassBC` | `Σ_i J_j,i` | built via `+` |

Pressure BCs (§5.5) are also `AbstractMassBC` subtypes but only act in Darcy
mode. The first four BC types above (and the pressure BCs) are exported by
`using Pyrolysis`; `CombinedMassBC` is produced by the `+` operator rather than
called directly.

### 5.3.1 `ImpermeableBC` — no gas escape

```julia
bc = ImpermeableBC()
```

Zero mass flux for every component. This is the default for both faces and the
usual choice for the substrate.

### 5.3.2 `MassFluxBC` — prescribed mass flux

```julia
bc = MassFluxBC(0.0)             # constant prescribed flux [kg/(m²·s)]
bc = MassFluxBC(t -> 1e-3)       # time-varying
```

Prescribes the boundary mass flux `J` directly (positive = into domain), as a
constant or `J(t)`. The same value applies to all gas components.

### 5.3.3 `ConvectiveMassBC` — external convective transfer

```julia
bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)   # gas swept away to clean atmosphere
```

External convective mass transfer, `J_j = h_m·(ξ_{j,∞} − ξ_{j,s})`. Fields:

- `h_m` — external (film) mass-transfer coefficient `h_m` in **m/s**
  (velocity-based). Constant or `h_m(t)`; validated non-negative. Typical
  values: natural convection 0.01–0.05 m/s, forced convection 0.05–0.5 m/s. As
  a rule of thumb from the heat–mass analogy,
  `h_m ≈ h_conv / (ρ_air · c_{p,air})`; for `h_conv ≈ 10` W·m⁻²·K⁻¹ this gives
  `h_m ≈ 0.008` m/s.
- `ξ_inf` — ambient mass concentration `ξ_∞` in kg/m³, default `0.0` (gas
  escapes into clean surroundings). Accepts a scalar (broadcast to every gas
  species), a function of time `t -> ξ_∞(t)`, or an `NTuple` with one entry
  per component (`ξ_inf = (0.0, 0.0, 0.27)` — e.g. ambient air ingress for a
  background-air species without injecting pyrolyzate). Entries are validated
  non-negative; a tuple's length must equal the material's component count.

### 5.3.4 `DiffusiveMassBC` — internal diffusion with optional external resistance

```julia
bc = DiffusiveMassBC()            # pure internal diffusion (h_m = Inf, default)
bc = DiffusiveMassBC(0.01)        # add a finite external film resistance, h_m = 0.01 m/s
```

Couples the **internal** Fickian diffusion gradient in the surface cell to an
optional **external** film resistance. The single field is:

- `h_m` — external mass-transfer coefficient in m/s. Default `Inf` (the
  no-argument constructor `DiffusiveMassBC()` sets `h_m = Inf`), meaning *no*
  external resistance — the flux is the pure diffusive limit
  `J_j = −λ_j ξ_{j,s}/d`, where `d` is the distance from the surface cell center
  to the boundary face. A finite `h_m` adds a series resistance:
  `J_j = −λ_j h_m ξ_{j,s}/(λ_j + h_m d)`.

This is the most physically consistent escape BC because it ties the surface
loss to the same internal transport law used everywhere else (it also keeps the
species balance consistent with ALE surface recession). Use it (with `Inf` or a
finite `h_m`) when you want gas escape governed by the material's own gas
transfer coefficient rather than an imposed convective coefficient.

### 5.3.5 `CombinedMassBC` — adding mass BCs with `+`

```julia
bc = DiffusiveMassBC(Inf) + ConvectiveMassBC(h_m = 0.01, ξ_inf = 0.0)
```

Sums component mass fluxes, `J_j = Σ_i J_{j,i}`, via the `+` operator (generated
function, type-stable). For most pyrolysis runs a single mass BC suffices;
combined mass BCs are for modeling competing transport mechanisms or
research/validation studies.

---

## 5.4 The surface energy balance (how thermal BCs are applied)

Thermal BCs are not applied as a raw flux in the cell; instead the solver finds
the **surface (face) temperature** `T_s` that makes the boundary energy balance
close, then conducts the resulting flux into the adjacent cell. With the default
convention (transpiration off), the balance is

```
q_BC(T_s) = k_eff · (T_s − T_cell) / Δz_{1/2}
```

— the applied BC flux equals the heat conducted from the surface face into the
first cell center, with `k_eff` the effective conductivity, `T_cell` the
adjacent cell temperature, and `Δz_{1/2}` the face-to-center distance. This is
solved by a damped Newton iteration on `T_s` (see Technical Reference §12.2). You
do not call this yourself — it runs inside the RHS assembly — but it is worth
knowing because it explains a few user-visible behaviors:

- **Convergence is robust** for physical BCs: convective flux is linear in `T_s`,
  radiative is `T_s⁴`; typical solves converge in 2–5 iterations (default
  tolerance `1e-10`, up to 50 iterations).
- **Surface temperature is clamped to `[200, 3000] K`** during the iteration to
  prevent unphysical excursions.
- The solver is **AD-aware**: it carries `ForwardDiff.Dual` types through the
  Newton iteration so sensitivities (see the sensitivity-analysis chapter) flow
  correctly through `T_s`.
- The converged surface temperature is available afterwards as
  `solution.surface_T`, and the flux decomposition (absorbed, re-radiation,
  convection, conducted) is available in `solution.extras` when diagnostics are
  enabled.

The **surface transpiration enthalpy** term `Σ_g ṁ_g Δh_g(T_s, T_∞)` is an
optional addition to this balance, controlled by `use_transpiration_bc`. It is
**off by default and rejected by `solve()` when set to `true`** — see the
warning in §5.12.

---

## 5.5 Pressure boundary conditions (Darcy mode only)

When the material has a permeability (so the solver runs in Darcy–Fick transport
mode, see Technical Reference §9), gas transport is pressure-driven and you can
prescribe **pressure** at a boundary. Pressure BCs are `AbstractMassBC` subtypes,
so they go in the `mass` slot of the `BoundaryConditionSet`. They do **not**
contribute a diffusive mass flux; instead they set the Darcy gas velocity at the
boundary face (via `apply_pressure_bc_to_darcy_velocity!`). In Fickian (non-Darcy)
mode they have no effect.

| BC type | Effect | Parameters (defaults) |
|---|---|---|
| `DirichletPressureBC` | fixed boundary pressure `P` | `P` (`101325.0`) |
| `NeumannPressureBC` | prescribed pressure gradient `∂P/∂n` | `∂P∂n` (`0.0` = impermeable) |
| `ConvectivePressureBC` | Robin outflow `v_g·n = h_P·(P_int − P_ext)` | `h_P`, `P_ext` |

All three are exported by `using Pyrolysis`.

### 5.5.1 `DirichletPressureBC` — fixed pressure

```julia
bc = DirichletPressureBC()           # default: atmospheric, 101325 Pa
bc = DirichletPressureBC(101325.0)   # explicit atmospheric
```

Holds the boundary at pressure `P` (Pa, validated positive). The Darcy velocity
is computed from the gradient between the boundary and the interior cell,
`v_g = −(κ/μ)·(P_bc − P_interior)/d`. This is the natural exposed-surface
condition: the outer face sits at atmospheric pressure and gas flows out down the
pressure gradient built up by pyrolysis.

### 5.5.2 `NeumannPressureBC` — prescribed gradient

```julia
bc = NeumannPressureBC()        # default: zero gradient (impermeable to gas flow)
bc = NeumannPressureBC(0.0)     # same as above
bc = NeumannPressureBC(1000.0)  # prescribed gradient [Pa/m] forcing an outflow
```

Sets the normal pressure gradient `∂P/∂n` (Pa/m), default `0.0`. A zero gradient
is a no-flow (impermeable) back face for the pressure field — the usual choice for
a sealed substrate in Darcy mode.

### 5.5.3 `ConvectivePressureBC` — Robin (resistive) outflow

```julia
bc = ConvectivePressureBC(1e-6, 101325.0)   # finite surface resistance, atmospheric exterior
```

A Robin-type condition, `v_g·n = h_P·(P_internal − P_external)`, for a surface
with finite resistance to gas flow (a char skin, a coating). Fields:

- `h_P` — pressure-transfer coefficient `h_P` in m·Pa⁻¹·s⁻¹.
- `P_ext` — external pressure `P_external` in Pa.

Limits: `h_P → ∞` approaches `DirichletPressureBC`; `h_P → 0` approaches an
impermeable (zero-flow) boundary.

---

## 5.6 Radiation boundary conditions (P1 model, Experimental)

The four radiation models are selected by the `radiation_model` keyword of
`solve()` (default `NO_RADIATION`), not by a BC in the set. For
`NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT`, the incoming radiation
is specified by a `RadiativeFluxBC` in the thermal slot (§5.2.6) — there is
nothing else to set.

The **P1 quasi-steady** model (`P1_QUASI_STEADY`) is different: it solves a
volumetric radiation field and needs its own boundary treatment (Marshak,
reflective, or Dirichlet on the intensity `G`). These P1 BCs are **experimental**
and live in `Pyrolysis.Experimental`, not in the `BoundaryConditionSet`:

```julia
import Pyrolysis.Experimental: MarshakRadiationBC, ReflectiveRadiationBC,
                               DirichletRadiationBC, RadiationBCSet, P1_QUASI_STEADY
```

The standard P1 configuration is a **Marshak** condition on the exposed top
(carrying the external flux and ambient field) and a **reflective** condition on
the insulated bottom. The Marshak flux is
`q_BC = c·(G_ext + 4·I_ext − G_surf)` with `c = ε/(2(2−ε))` and
`G_ext = 4σT_∞⁴` — the external flux carries the half-range weight
`4c = 2ε/(2−ε)`, the same weight as the ambient term (see Technical Reference
§8.6 for the derivation).

> **Status caveat.** P1 radiation is in `Pyrolysis.Experimental` and runs through
> a separate residual path. For most cone-calorimeter and TGA work,
> `SURFACE_ABSORPTION` (opaque) or `BEER_LAMBERT` (semi-transparent) with a
> plain `RadiativeFluxBC` is the supported, production route. Reach for P1 only
> when you specifically need in-depth re-radiation in an optically thick sample.

---

## 5.7 Mesh motion in ALE simulations (no boundary conditions needed)

ALE mesh motion is **not configured through boundary conditions** — there are no
mesh-motion BC types, and `BoundaryConditionSet` carries only `thermal` and
`mass` slots. The mesh-velocity boundary treatment is fixed by the formulation:
the bottom node is pinned (`w_1 = 0`, hard-coded in the residual's mesh-velocity
integration) and every other node's velocity follows cumulatively from the
volumetric strain as the material shrinks or swells (see Technical Reference
§11.3), so the top surface automatically follows the material (recedes or
expands). You enable the moving mesh simply with `solve(problem; use_ale =
true)`.

---

## 5.8 Choosing ambient temperatures

Several thermal BCs carry an ambient temperature `T_∞`: `ConvectiveBC` and
`RadiativeBC` each have their own `T_inf`. This `T_∞` is used both in the BC's own
flux and as the **reference temperature for the surface energy balance**. The
resolution rules (see Technical Reference §12.3.8) are:

- `ConvectiveBC` and `RadiativeBC` use their own `T_inf`.
- BCs without an ambient (`HeatFluxBC`, `AdiabaticBC`, `RadiativeFluxBC`,
  `DirichletThermalBC`) carry none, and the resolver falls back to the package
  default `T_REF = 298.15 K`.
- For a `CombinedThermalBC`, the first child that carries an ambient wins — so a
  combined top BC with a `RadiativeBC(T_inf = 330)` and a `ConvectiveBC(T_inf = 330)`
  resolves the surface ambient to 330 K. (Resolution keys on whether a BC *has* an
  ambient, not on its value: setting `T_inf = 298.15` exactly is honored, never
  mistaken for "defaulted".)

Practical guidance: set `T_inf` consistently on the loss terms (`RadiativeBC`,
`ConvectiveBC`) to the experimental surroundings temperature. In a cone
calorimeter the cone shroud and chamber are warmer than room temperature; values
in the 300–360 K range are common and used in the shipped examples.

---

## 5.9 Time-dependent parameters

Any BC parameter shown above as "constant or function of `t`" accepts a
single-argument function. The solver evaluates it at the current time. Examples:

```julia
# Heat-flux ramp
HeatFluxBC(t -> 50e3 * min(t / 5.0, 1.0))           # linear ramp to 50 kW/m² over 5 s

# Time-varying incident radiation
RadiativeFluxBC(flux = t -> 25e3 * (1 - exp(-t / 10)))

# Warming ambient in a convective loss term
ConvectiveBC(h = 10.0, T_inf = t -> 300 + 0.5 * t)

# Pulsed mass extraction
MassFluxBC(t -> t < 50 ? 0.0 : -1e-3)
```

> **Note on time-varying radiation flux.** When `radiation_model` is
> `SURFACE_ABSORPTION` or `BEER_LAMBERT`, `solve()` re-samples the incident
> flux of the top `RadiativeFluxBC` at **every** RHS evaluation — constant and
> time-varying fluxes alike, with no up-front constancy detection. A ramped or
> modulated flux is therefore always tracked exactly. You do not need to do
> anything special — just supply `flux = t -> …`.

---

## 5.10 Initial conditions

Initial conditions are supplied to `PyrolysisProblem` as functions of the spatial
coordinate `z` (in meters, with `z = 0` at the bottom, `z = L` at the top). Two
inputs:

- `T_initial` — a function `z -> T` returning the initial temperature in K.
  Default: `z -> 300.0`.
- `ξ_initial` — the initial mass concentrations `ξ_j(z)` in kg/m³, one entry per
  component. Accepts a `Vector` of functions (one per component), a single
  function (broadcast to all components), or you may omit it for the auto-default.

```julia
problem = PyrolysisProblem(
    mesh = mesh,
    material = material,            # NC = 2 here: virgin solid + gas
    bc_set = bc_set,
    T_initial = z -> 300.0,                  # uniform 300 K
    ξ_initial = [z -> 1190.0, z -> 0.0],     # full virgin density, no gas
    tspan = (0.0, 100.0),
)
```

### 5.10.1 Auto-defaults

If you omit `ξ_initial`, the constructor builds a sensible default: the **first
component, if it is a solid, is placed at its reference density** (its density
function evaluated at 300 K) and all other components are set to zero. This is
correct for the common "single virgin solid that decomposes into gas(es) and
char" setup, where component 1 is the virgin solid. Always double-check it
matches your component ordering; if not, pass `ξ_initial` explicitly.

`T_initial` defaults to a uniform 300 K.

### 5.10.2 Uniform vs. depth-varying profiles

A uniform value is just a constant-returning closure:

```julia
T_initial = z -> 300.0
ξ_initial = [z -> 1190.0, z -> 0.0]   # uniform composition
```

For a non-uniform start (e.g. a pre-charred surface layer, or a moisture
gradient) make the closures depend on `z`:

```julia
L = 0.02  # domain thickness [m]
# Warmer at the exposed top, cooler at the substrate
T_initial = z -> 300.0 + 50.0 * (z / L)

# A thin char layer near the surface (component 2 = char), virgin below
ξ_initial = [
    z -> z > 0.9 * L ? 200.0 : 1190.0,   # virgin: depleted near top
    z -> z > 0.9 * L ? 150.0 : 0.0,      # char: present near top
    z -> 0.0,                            # gas: none initially
]
```

The initial conditions are evaluated at each cell center and stored in the mesh.
They also populate the per-cell **dry-solid mass** denominator used by the
lateral-shrinkage (χ) dynamics, computed as
`m_{dry,i} = V_i · Σ_{solid} ξ_{j,i}(0)` (see Technical Reference §12.8 and
§10.3).

### 5.10.3 Profile helper functions

Convenience profile builders are available in the `Problem` submodule (not
top-level exported). Import them explicitly when you want them:

```julia
import Pyrolysis.Problem: step_profile, linear_profile, gaussian_profile, tanh_profile

z_int = 0.9 * L
# Sharp step: char above the interface, virgin below
ξ_char_initial = z -> step_profile(z; z_interface = z_int,
                                       value_left = 0.0, value_right = 150.0)

# Linear temperature ramp through the thickness
T_initial = z -> linear_profile(z; z_min = 0.0, z_max = L,
                                    value_min = 300.0, value_max = 350.0)

# Smooth (tanh) interface instead of a hard step
ξ_char_initial = z -> tanh_profile(z; z_center = z_int, width = 0.001,
                                      value_left = 0.0, value_right = 150.0)

# Gaussian bump (e.g. a localized moisture pocket)
ξ_moist_initial = z -> gaussian_profile(z; z_center = 0.5 * L, width = 0.002,
                                           amplitude = 100.0, baseline = 0.0)
```

Their signatures:

| Helper | Keyword arguments |
|---|---|
| `step_profile(z; z_interface, value_left, value_right)` | step at `z_interface` |
| `linear_profile(z; z_min, z_max, value_min, value_max)` | clamped linear ramp |
| `gaussian_profile(z; z_center, width, amplitude, baseline = 0.0)` | Gaussian bump |
| `tanh_profile(z; z_center, width, value_left, value_right)` | smooth step |

A smooth (`tanh`) interface is generally preferable to a hard `step_profile` for
the solver: sharp discontinuities at `t = 0` can provoke small startup
transients.

---

## 5.11 Common configurations (worked recipes)

### 5.11.1 Cone calorimeter (opaque charring solid)

External cone flux absorbed at the surface, with re-radiation and convective
losses; gas escapes by convective mass transfer; insulated, impermeable back
face:

```julia
bc_set = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.96) +
              RadiativeBC(ε = 0.96, T_inf = 330.0) +
              ConvectiveBC(h = 10.0, T_inf = 330.0),
        bottom = AdiabaticBC(),
    ),
    mass = (
        top    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
        bottom = ImpermeableBC(),
    ),
)
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)
```

### 5.11.2 Semi-transparent polymer (in-depth absorption)

Same incoming radiation, but absorbed in depth via Beer–Lambert (requires
absorption coefficients on the components):

```julia
bc_set = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 50e3) +          # α defaults to ε_eff
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
)
solution = solve(problem; radiation_model = BEER_LAMBERT)
```

### 5.11.3 Furnace temperature program on the exposed face

A prescribed heating program is most directly modeled by holding the surface
temperature with a time-dependent `DirichletThermalBC`:

```julia
β = 10.0 / 60.0   # 10 K/min in K/s
bc_set = BoundaryConditionSet(
    thermal = (top = DirichletThermalBC(t -> 300.0 + β * t),
               bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
)
```

### 5.11.4 Pressure-driven char (Darcy mode)

With a permeable material, set atmospheric pressure on the exposed face and a
sealed (zero-gradient) back face. The transport mode switches to Darcy–Fick
automatically when the material has permeability:

```julia
bc_set = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = DirichletPressureBC(101325.0),    # atmospheric exposed face
               bottom = NeumannPressureBC(0.0)),       # sealed substrate
)
```

---

## 5.12 Caveats and pitfalls

- **Surface transpiration is forbidden in the production path.** `solve()`
  rejects `use_transpiration_bc = true` with an `ArgumentError`. The reason
  (Technical Reference §7.6, §12.3): the volumetric gas-advection source
  `S_conv = −Σ_g c_{p,g} J_g ∂T/∂z` is *always* active and already carries gas
  enthalpy; adding the surface transpiration term `Σ_g ṁ_g Δh_g(T_s, T_∞)`
  double-counts it, producing a spurious ~5–10 kW/m² sink that under-predicts
  surface temperature and mass-loss rate. The supported convention (matching
  ThermaKin, Gpyro, and FDS) is `S_conv` only, i.e. `use_transpiration_bc =
  false` (the default). Gas mass transport still works either way — the
  transpiration term affects only the *energy* balance.

- **`radiation_model` defaults to `NO_RADIATION`.** A `RadiativeFluxBC` only does
  something if you pass an appropriate `radiation_model` to `solve()`. With the
  default `NO_RADIATION`, an incident `RadiativeFluxBC` is still absorbed at the
  surface (the absorbed flux is applied as a boundary heat flux), but choose
  `SURFACE_ABSORPTION` (opaque) or `BEER_LAMBERT` (semi-transparent) explicitly
  to make the intent clear and to enable the in-depth distribution.

- **Only thermal and mass BCs go in the set.** `BoundaryConditionSet` accepts
  only `AbstractThermalBC` and `AbstractMassBC` payloads — a wrong-kind BC
  throws at construction. ALE mesh motion needs no BC at all (it is automatic
  with `use_ale = true`), and P1 radiation BCs (`MarshakRadiationBC` etc.) go
  through the Experimental P1 residual path, not the set.

- **Pressure BCs need Darcy mode.** `DirichletPressureBC`/`NeumannPressureBC`/
  `ConvectivePressureBC` only have an effect when the material has permeability
  (Darcy–Fick transport). In Fickian mode they contribute zero flux and are
  silently inert — use `ImpermeableBC`/`DiffusiveMassBC`/`ConvectiveMassBC` for
  diffusion-only gas escape.

- **Don't double-specify a radiation source.** Use `RadiativeFluxBC` (not also
  `HeatFluxBC`) for the cone/panel; use `RadiativeBC` for the surface's own
  re-radiation loss. They are different physical terms.

- **Match `ξ_initial` to your component order.** The auto-default only fills
  component 1 (if solid). For multi-solid or moisture-bearing materials, supply
  `ξ_initial` explicitly so virgin solid, char, moisture, and gases start with
  the right concentrations.

For deeper background on the boundary physics, the surface Newton solve, and the
radiation models, see Technical Reference §12 (boundary and initial conditions),
§8 (thermal radiation), §7 (energy assembly and the `S_conv` convention), and §9
(gas transport and pressure).


---

# 6. Meshes and Geometry

The mesh is the spatial backbone of every Pyrolysis.jl simulation: it discretizes the
through-thickness coordinate `z` into the finite-volume cells over which temperature and
species concentrations are solved. Pyrolysis.jl is a strictly **one-dimensional** (1D)
through-thickness solver, so "geometry" here means a stack of cells along `z`, plus a single
scalar cross-sectional area `A` that turns those 1D cells into physical volumes. This chapter
shows how to create meshes (uniform and stretched), how to choose resolution, what the
cross-sectional area `A_0` means in a 1D model, how to inspect a mesh, and how the spatial
coordinate is laid out in the raw state vector — including the special case where the mesh
itself becomes part of the solved state (ALE). For the underlying finite-volume theory see
Technical Reference §13; for the moving-mesh (ALE) formulation see Technical Reference §11.

> **Coordinate convention (read this first).** Throughout Pyrolysis.jl, `z = 0` is the
> **bottom / substrate** boundary (tag `:bottom`, internal `boundary_id = 2`) and `z = L` is
> the **top / exposed surface** (tag `:top`, internal `boundary_id = 1`). Heat enters from the
> top; the material recedes toward the bottom. Cell indices increase from bottom (cell 1, near
> `z = 0`) to top (cell `n`, near `z = L`). Positive flux is transport in the `+z` direction.
> (You may see a stray comment in `examples/06_cone_calorimeter.jl` saying "z=0 is the exposed
> TOP surface" — that comment is wrong; the mesh generator and the rest of the package use the
> convention stated here.)

---

## 6.1 The mesh in one dimension

A 1D mesh in Pyrolysis.jl consists of:

- **Nodes** — the vertices that bound cells. For `n` cells there are `n + 1` nodes. Node 1
  sits at `z = 0` (bottom), node `n+1` sits at `z = L` (top).
- **Cells** — the finite-volume control volumes. Cell `i` spans nodes `i` and `i+1`, i.e. the
  interval `[z_i, z_{i+1}]`. Temperature `T_i` and the per-component mass concentrations
  `ξ_{j,i}` are stored at the cell **center** `z_i^c = (z_i + z_{i+1}) / 2`.
- **Faces** — in 1D each face is a single point that co-locates with a node, so there are also
  `n + 1` faces. Fluxes (conduction, gas transport) are evaluated at faces. Faces 1 and `n+1`
  are the boundary faces.

The mesh container type is `Mesh1D{NC}`, parameterized by the number of chemical components
`NC` (the number of `SolidComponent`/`LiquidComponent`/`GaseousComponent` entries in your
`Material`). You never construct `Mesh1D` directly; you call one of the two mesh-generation
functions described in §6.2. The relationship `n_faces = n_nodes = n_cells + 1` always holds.

Each cell carries a `CellState{NC}` (temperature, an `NTuple{NC}` of concentrations, and a
cached volume); these are populated from your initial conditions when you call `solve`. The
volume of cell `i` is computed geometrically from the node positions and the cross-sectional
area:

```
V_i = |z_{i+1} − z_i| · A
```

See Technical Reference §13 for the full description of the cell/node/face layout and the
divergence operator built on it.

---

## 6.2 Creating a mesh

Two mesh creators are exported at the top level (`using Pyrolysis` makes both available):

| Function | Purpose |
|---|---|
| `create_uniform_mesh(L, n_cells, Val(NC); ...)` | Equal cell spacing `Δz = L / n_cells`. |
| `create_stretched_mesh(L, n_cells, Val(NC), stretch; ...)` | Geometrically graded cells (fine near one surface, coarse near the other). |

Both return a `Mesh1D{NC}`.

### 6.2.1 Uniform meshes

```julia
create_uniform_mesh(L::Real, n_cells::Int, ::Val{NC};
                    T_initial::Real = 300.0,
                    ξ_initial = nothing,
                    cross_section_area::Real = 1.0) -> Mesh1D{NC}
```

Positional arguments:

- `L` — domain thickness in metres (`m`). This is the total material length `z = L`.
- `n_cells` — number of finite-volume cells.
- `Val(NC)` — the number of chemical components, passed as a **type parameter** via
  `Val(...)`. `NC` must equal the number of components in the `Material` you intend to use.
  Passing a plain integer (e.g. `3` instead of `Val(3)`) will not compile.

Keyword arguments and their defaults:

| Keyword | Default | Meaning |
|---|---|---|
| `T_initial` | `300.0` | Uniform initial cell temperature `[K]`. Note: `solve` re-applies the problem's `T_initial` closure, so this matters mainly when you inspect the bare mesh. |
| `ξ_initial` | `nothing` | Initial mass concentrations `[kg·m⁻³]`. `nothing` → all zero; a scalar → that value for every component; an `NTuple`/vector of length `NC` → per-component values. |
| `cross_section_area` | `1.0` | Cross-sectional area `A_0` `[m²]` (see §6.4). |

A 1 cm slab with 100 cells and 2 components:

```julia
using Pyrolysis

mesh = create_uniform_mesh(0.01, 100, Val(2))   # L = 1 cm, 100 cells, NC = 2
```

A typical cone-calorimeter sample (5.8 mm, 72 cells, 3 components — virgin solid, char,
pyrolysis gas — initialized with the virgin solid at its bulk density):

```julia
sample_thickness = 5.8e-3      # 5.8 mm  [m]
n_cells          = 72          # 0.08 mm resolution

mesh = create_uniform_mesh(
    sample_thickness,
    n_cells,
    Val(3);                              # NC = 3: virgin, char, gas
    T_initial = 300.0,
    ξ_initial = (1210.0, 0.0, 0.0),      # full virgin density, no char, no gas
)
```

The node positions are `z_i = (i − 1)·Δz` with `Δz = L / n_cells`, so cell `i` occupies
`[(i−1)Δz, iΔz]` and every cell has the same volume `V_i = Δz · A`.

### 6.2.2 Stretched (graded) meshes

In a cone-calorimeter or fire-exposure problem the steepest temperature and reaction-rate
gradients are at the **top** (exposed) surface, where the pyrolysis front lives. A graded mesh
puts small cells where the action is and large cells deep in the cold substrate, giving the
same accuracy for fewer total cells.

```julia
create_stretched_mesh(L::Real, n_cells::Int, ::Val{NC}, stretch::Real;
                      T_initial::Real = 300.0,
                      ξ_initial = nothing,
                      cross_section_area::Real = 1.0,
                      stretch_from::Symbol = :top) -> Mesh1D{NC}
```

The extra arguments beyond the uniform case:

| Argument | Default | Meaning |
|---|---|---|
| `stretch` (positional) | — | Geometric stretching factor. `stretch = 1` ⇒ uniform (any value with `\|stretch − 1\| < 1e-10` is treated exactly as uniform). `stretch > 1` makes successive cells larger by that factor along the stretching direction. Must be `≥ 1` — a sub-unity value throws an `ArgumentError`; choose which end is refined with `stretch_from`, never with `stretch < 1`. |
| `stretch_from` (keyword) | `:top` | Which surface the **fine** end is anchored to: `:top` puts the smallest cells at `z = L` (the exposed surface); `:bottom` puts the smallest cells at `z = 0`. |

The cell thicknesses follow a geometric series. With `stretch_from = :top` the smallest cell is
adjacent to the exposed surface and the spacing grows geometrically away from it; the first cell
thickness `Δz_0` is chosen so the cells sum exactly to `L`:

```
Δz_0 = L · (stretch − 1) / (stretch^n − 1)
```

Fine cells at the exposed surface, coarse toward the cold back:

```julia
# 50 cells over 1 cm, smallest cells at the top (exposed) surface
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from = :top)
```

A modest factor (`1.05`–`1.15`) is usually enough. Very aggressive stretching produces large
cell-to-cell size jumps, which can hurt accuracy and conditioning; the conduction and transport
schemes use distance-weighted harmonic-mean face properties precisely so that graded meshes
stay second-order accurate (Technical Reference §7 and §13), but extreme ratios still degrade
the gradient reconstruction. If you stretch hard, sanity-check the cell-size jump with
`Pyrolysis.Geometry.aspect_ratio` (§6.5).

> **Tip.** When in doubt, start uniform, run a convergence study (§6.3), and only switch to a
> stretched mesh if a uniform mesh with the same near-surface resolution is too expensive.

---

## 6.3 Choosing `n_cells`: resolution versus cost

The number of cells controls both accuracy and run time. There is no universal "correct"
`n_cells`; it depends on how thin the reacting/thermal layer is relative to the sample. Some
practical guidance:

- **Resolve the thermal/reaction front.** For high-flux cone tests the front can be a fraction
  of a millimetre thick. As a rule of thumb you want several cells across that layer. The cone
  example above uses `72` cells over `5.8 mm` (≈ `0.08 mm`/cell). The Douglas-fir convergence
  example (`examples/convergence_douglas_fir.jl`) sweeps `n_cells` and reports cell size in mm,
  peak mass-loss rate, and time-to-ignition so you can see when the answer stops moving.
- **Cost scales with the state size.** The solved state has length
  `n · (1 + NC)` for a fixed (Eulerian) mesh, plus `n_nodes` more for ALE and `n` more for a
  variable cross-section (see §6.6). The stiff implicit solver factorizes a Jacobian whose
  size grows with that state, so doubling `n_cells` more than doubles the per-step cost.
- **Run a convergence study.** Solve the same problem at increasing `n_cells` and watch a
  scalar output of interest (peak MLR, ignition time, back-face temperature) converge. Stop
  refining when successive results agree to your tolerance. This is the recommended way to pick
  a production resolution. See Technical Reference §18 for the verification/convergence
  methodology and the `convergence_douglas_fir` example for a worked spatial-convergence sweep.

A minimal convergence-loop sketch (build the same problem at each resolution, extract one
scalar, compare):

```julia
using Pyrolysis

results = Dict{Int,Float64}()
for n in (24, 48, 96, 192)
    mesh = create_uniform_mesh(5.8e-3, n, Val(3); ξ_initial = (1210.0, 0.0, 0.0))
    problem = PyrolysisProblem(
        mesh      = mesh,
        material  = material,                       # your Material{3,...}
        bc_set    = bcs,                            # your BoundaryConditionSet
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
        tspan     = (0.0, 1600.0),
    )
    sol = solve(problem; radiation_model = SURFACE_ABSORPTION, progress = false)
    results[n] = maximum(sol.mass_loss_rate)        # peak MLR  [kg·m⁻²·s⁻¹]
end
```

When the entries of `results` agree to your target tolerance, you have enough cells.

---

## 6.4 Cross-sectional area `A_0` in a 1D model

A 1D through-thickness model has no real lateral extent, so the cross-sectional area is a
**gauge / normalization constant** rather than a physical width. It enters the model only
through cell volumes and face areas:

```
V_i = |z_{i+1} − z_i| · A          A_face = A   (same A on every face in 1D)
```

Because the governing equations are written per unit volume and per unit area, the area cancels
out of intensive results (temperature profiles, MLR per unit area, surface temperature) **as
long as `A` is constant**. The default `cross_section_area = 1.0 m²` is therefore the natural
choice for almost every simulation, and it makes per-area diagnostics read directly. Set it to
something else only if you have a specific reason to track absolute (extensive) quantities such
as total mass or total energy in physical units.

The mesh stores two area fields:

- `cross_section_area` — the **current** area `A(t)` `[m²]`. For a fixed-cross-section problem
  this never changes.
- `cross_section_area_initial` — the fixed **calibration reference** `A_0` `[m²]`, set once at
  mesh creation.

These differ only when you enable a **lateral-shrinkage law** on the material, which makes the
area evolve as `A(t) = A_0 · law(χ̄)` with `χ̄` the column-averaged pyrolysis progress. That is
the "variable cross-section" (`WithChi`) mode; it is triggered by attaching a
`lateral_shrinkage_law` to the `Material` (covered in Chapter 7 and Technical Reference §10).
When the area varies, mass-loss-rate diagnostics are normalized to the original area `A_0`
(the `MLR_gauge` extra), so they remain comparable to fixed-area runs.

You pass the area at mesh-creation time:

```julia
# 100 cm² sample face, otherwise identical to the default-area mesh
mesh = create_uniform_mesh(0.01, 100, Val(2); cross_section_area = 1e-2)
```

---

## 6.5 Inspecting a mesh

Once you have a `Mesh1D`, you can query its size and geometry. **Note on imports:** the basic
size accessors and geometry queries are exported by the `Geometry` submodule but are *not* part
of the curated top-level `Pyrolysis` surface. Reach them either fully qualified
(`Pyrolysis.Geometry.cell_center`, etc.) or by bringing in the umbrella namespace with
`using Pyrolysis.Internal`, which re-exports every submodule's public symbols. The two mesh
*creators* (`create_uniform_mesh`, `create_stretched_mesh`) are top-level exported; the mesh
*type* (`Mesh1D`) is exported by the `Geometry` submodule — reach it as
`Pyrolysis.Geometry.Mesh1D` or through the `Pyrolysis.Internal` umbrella — not from the
curated top-level surface.

```julia
using Pyrolysis
using Pyrolysis.Internal   # makes the accessors below available unqualified
```

### 6.5.1 Size accessors

| Call | Returns |
|---|---|
| `n_cells(mesh)` | Number of cells `n`. |
| `n_nodes(mesh)` | Number of nodes (`= n + 1`). |
| `n_faces(mesh)` | Number of faces (`= n + 1`). |
| `n_active_cells(mesh)` | Number of active (leaf) cells — equals `n_cells` unless adaptive refinement has been applied (Chapter 10). |

```julia
mesh = create_uniform_mesh(0.01, 100, Val(2))
n_cells(mesh)        # 100
n_nodes(mesh)        # 101
n_faces(mesh)        # 101
```

### 6.5.2 Geometry queries

| Call | Returns |
|---|---|
| `cell_center(mesh, i)` | Center position `z_i^c` of cell `i` `[m]`. |
| `cell_size(mesh, i)` | Thickness `Δz_i = \|z_{i+1} − z_i\|` of cell `i` `[m]`. |
| `cell_volume(mesh, i)` | Cached volume `V_i` of cell `i` `[m³]`. |
| `face_area(mesh, f)` | Area `A` of face `f` `[m²]` (constant in 1D). |
| `total_volume(mesh)` | Sum of active-cell volumes `[m³]`. |
| `domain_bounds(mesh)` | Tuple `(z_min, z_max)` from node positions `[m]`. |
| `boundary_face_ids(mesh, boundary_id)` | Vector of face indices for a boundary id (`1` = top, `2` = bottom). |
| `cross_section_area(mesh)` | Current `A(t)` `[m²]` (top-level exported). |

```julia
domain_bounds(mesh)              # (0.0, 0.01)
cell_center(mesh, 1)             # ≈ 5.0e-5   (first cell, near z = 0, the bottom)
cell_center(mesh, n_cells(mesh)) # ≈ 9.95e-3 (last cell, near z = L, the top/exposed surface)
cell_size(mesh, 1)              # 1.0e-4    (uniform Δz)
total_volume(mesh)              # 0.01 * 1.0  (length × default area = 0.01 m³)

# Face indices for the two boundaries (id 1 = top/exposed, id 2 = bottom/substrate)
boundary_face_ids(mesh, 1)       # [101]  — top face at z = L
boundary_face_ids(mesh, 2)       # [1]    — bottom face at z = 0
```

For the user-facing boundary tags (`:top`, `:bottom`) that you use when specifying boundary
conditions, see Chapter 5; internally `:top` maps to id `1` and `:bottom` to id `2`.

### 6.5.3 Cell-by-cell positions (for plotting)

To get all cell-center positions as a vector — for plotting a temperature profile against
depth — read `mesh.cell_centers` directly, or, after a solve, use the solution accessor
`get_cell_positions` (top-level exported via the Solver module):

```julia
z = mesh.cell_centers                 # Vector{Float64} of length n_cells  [m]

# after solving:
z = get_cell_positions(solution)      # same vector, taken from solution.problem.mesh
```

### 6.5.4 Mesh-quality metrics (graded/adaptive meshes)

For stretched or adaptively refined meshes you can check grading quality with the metrics in
the `Geometry` submodule. These four are defined but **not exported** (not even through
`Pyrolysis.Internal`), so call them **fully qualified** with the `Pyrolysis.Geometry.` prefix:

| Call | Meaning |
|---|---|
| `Pyrolysis.Geometry.aspect_ratio(mesh, i)` | Largest size ratio between cell `i` and its neighbours (`1.0` = uniform locally). |
| `Pyrolysis.Geometry.max_aspect_ratio(mesh)` | Worst cell-to-neighbour size ratio over the active mesh. |
| `Pyrolysis.Geometry.min_cell_size(mesh)` | Smallest active cell thickness `[m]`. |
| `Pyrolysis.Geometry.max_cell_size(mesh)` | Largest active cell thickness `[m]`. |

```julia
using Pyrolysis
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from = :top)
Pyrolysis.Geometry.max_aspect_ratio(mesh)   # ≈ 1.1 for stretch = 1.1
Pyrolysis.Geometry.min_cell_size(mesh)      # the smallest cell, at the exposed surface
```

A `max_aspect_ratio` close to your `stretch` factor is expected and fine. Values far above ~2
indicate aggressive grading; consider more cells or a gentler stretch.

---

## 6.6 The state-vector layout (for users who index raw state)

Most users never touch the raw ODE state vector — `solve` returns a `PyrolysisSolution` with
named fields (`T`, `ξ`, `z`, `surface_T`, `mass_loss_rate`, …; see Chapter 9). But if you write
a callback, a custom diagnostic, or a sensitivity output function, you need to know how the flat
state vector `u` is packed. The packing is **block-major** and is described by a `StateLayout`,
which is itself derived from the problem's `SimulationMode` (Chapter 7). With `n = n_cells` and
`NC` components, the blocks appear in this order:

| Block | Present when | Indices | Contents |
|---|---|---|---|
| `T` | always | `1 .. n` | Cell temperatures `[K]` |
| `ξ_j` | always | `n + (j−1)·n + (1..n)` for `j = 1..NC` | Mass concentrations `ξ_{j,i}` `[kg·m⁻³]` |
| `z` | ALE only | `n·(1+NC) + (1..n_nodes)` | Node positions `[m]` |
| `χ` | variable cross-section only | trailing `n` entries | Pyrolysis progress `χ_i` (dimensionless) |

The total length is `state_length(layout) = n·(1+NC) [+ n_nodes if ALE] [+ n if WithChi]`.

### 6.6.1 Reading state safely with accessors

Rather than computing offsets by hand, use the layout-aware accessors (top-level exported):

```julia
# `layout` is a StateLayout; `u` is a raw state vector (e.g. from a callback)
T_i  = get_temperature(u, layout, i)        # temperature of cell i  [K]
ξ_ji = get_concentration(u, layout, i, j)   # concentration of comp j in cell i  [kg·m⁻³]

# block views (no copies): Val(:T), Val(:ξ) (with component index), Val(:z), Val(:χ)
T_block  = block_view(u, layout, Val(:T))       # view of all temperatures
ξ1_block = block_view(u, layout, Val(:ξ), 1)    # view of component 1's concentrations
```

ALE-only accessors return zero / empty for non-ALE layouts, so they are safe to call
unconditionally:

```julia
z_i        = get_node_position(u, layout, i)         # node position  [m]  (ALE only)
z_nodes    = get_node_positions(u, layout)           # view of the whole z block (ALE only)
thickness  = get_thickness_from_state(u, layout)     # z_top − z_bottom  [m]; 0 for Eulerian
```

You can build a `StateLayout` directly from the component count, mesh size, and mode:

```julia
layout = StateLayout(NC, n_cells(mesh), n_nodes(mesh), mode)   # mode :: SimulationMode
```

In practice the layout is created for you inside the solver workspace; when you write a callback
you receive it (or can read it from the workspace). See Technical Reference §15 for the full
state-layout descriptor and the residual that consumes it, and Chapter 7 for how
`SimulationMode` is selected.

---

## 6.7 Eulerian versus ALE geometry: when `z` becomes state

There are two ways Pyrolysis.jl handles the mesh during a solve, controlled by the `use_ale`
keyword of `solve` (default `false`):

### 6.7.1 Eulerian (fixed mesh) — the default

With `use_ale = false`, the node positions are **fixed** for the whole simulation. Geometry
(cell centers, volumes, face positions) is constant and is mirrored once when the solver builds
its workspace. The state vector contains only the `T` and `ξ` blocks (plus `χ` if you have a
variable cross-section). This is the right choice for almost all cases: TGA-like idealizations,
cone-calorimeter runs without explicit surface recession, and any problem where you do not need
the mesh to physically move with the receding/swelling material. Surface recession in the
Eulerian picture is represented through the density/concentration fields and optional adaptive
surface-cell depletion (Chapter 10), not by moving nodes.

```julia
sol = solve(problem; progress = false)            # Eulerian (use_ale = false, the default)
```

### 6.7.2 ALE (moving mesh) — `use_ale = true`

With `use_ale = true`, the node positions `z_i` become **part of the solved state**. The mesh
moves with the material: solid components are tracked Lagrangian-style (locked to the nodes)
while gas is treated Eulerian, so the mesh contracts as the solid pyrolyzes and shrinks (or
expands for an intumescent/swelling material). The state vector then includes the `z` block
(see §6.6), and the residual computes a mesh velocity `w` cumulatively from the pinned bottom
node upward, advancing nodes as `dz_i/dt = w_i`. The cell volumes are recomputed geometrically
from the moving nodes each evaluation, which keeps the discrete Geometric Conservation Law
satisfied to machine precision (Technical Reference §11).

```julia
sol = solve(problem; use_ale = true, progress = false)   # moving (Lagrangian) mesh
```

In ALE mode the solution carries the node-position history in `solution.z` (shape
`(n_nodes, n_steps)`; empty for Eulerian runs), so you can plot mesh motion and surface
recession over time. Use ALE when surface regression or volume change is a first-class output
you want resolved geometrically; otherwise prefer the cheaper Eulerian default. ALE is also a
prerequisite for the variable-cross-section (lateral shrinkage) physics, since the column
strain and the area law are coupled through the mesh-velocity computation.

> **Caveat — vanishing cross-section / mesh tangling.** ALE meshes can in principle tangle
> (nodes crossing) or shrink toward zero thickness for very aggressive recession. The solver
> guards against this with the `min_thickness` keyword of `solve` (default `1e-6 m`) and
> mesh-quality checks; if a run terminates early or reports a thickness floor, see the
> troubleshooting notes in Chapter 15 and the ALE quality discussion in Technical Reference §11.

### 6.7.3 Geometry derived from state (advanced)

When you do need geometry from a raw ALE state vector (for a custom diagnostic), the
state-layout helpers compute it without touching the mutable mesh:

```julia
thickness = get_thickness_from_state(u, layout)   # current z_top − z_bottom  [m]
```

The lower-level routines `compute_geometry_from_state!` / `update_geometry_cache!`
(in `Pyrolysis.Geometry`) fill pre-allocated arrays of volumes, centers, and face distances
directly from the `z` block; these are what the residual uses internally and are documented in
Technical Reference §13. For variable-cross-section runs, `cross_section_area_from_state` and
`column_chi_bar_from_state` (in `Pyrolysis.Residual`) reconstruct `A(t)` and `χ̄` from the
state — advanced hooks you will rarely need directly.

---

## 6.8 Summary and cross-references

- Build meshes with `create_uniform_mesh(L, n_cells, Val(NC); ...)` or
  `create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from = :top, ...)`. Both return
  `Mesh1D{NC}`; `NC` must match your `Material`.
- `z = 0` is the bottom/substrate, `z = L` is the top/exposed surface; cells number from bottom
  to top.
- Choose `n_cells` by resolving the near-surface front; confirm with a convergence study
  (Technical Reference §18; `convergence_douglas_fir` example). Stretch toward `:top` to
  concentrate cells at the exposed surface.
- `cross_section_area` (`A_0`, default `1.0 m²`) is a gauge constant; it cancels from intensive
  results and only matters as `A(t) = A_0·law(χ̄)` when a lateral-shrinkage law is active
  (Technical Reference §10).
- Inspect meshes with `n_cells`/`n_nodes`/`n_faces`, `cell_center`, `cell_size`, `cell_volume`,
  `total_volume`, `domain_bounds`, and `boundary_face_ids` — qualify with `Pyrolysis.Geometry.`
  or bring them in via `using Pyrolysis.Internal`.
- The state vector is block-major `[T | ξ_1..ξ_NC | z (ALE) | χ (WithChi)]`; use
  `get_temperature`, `get_concentration`, `block_view`, and the ALE accessors instead of raw
  offsets.
- `use_ale = false` (default) keeps the mesh fixed (Eulerian); `use_ale = true` makes node
  positions part of the state so the mesh moves with the material (Technical Reference §11).

For the finite-volume discretization built on this mesh, see Technical Reference §13; for
moving-mesh mechanics and the Geometric Conservation Law, Technical Reference §11; for how the
mesh plugs into a problem and how modes are selected, Chapter 7; for boundary conditions on the
top/bottom faces, Chapter 5; and for adaptive refinement and surface depletion, Chapter 10.


---

# 7. Building a Problem: Simulation Modes

By this point you can build the three ingredients of a simulation — a mesh
(Chapter 6), a material with components, reactions, and properties (Chapters 3–4),
and a `BoundaryConditionSet` (Chapter 5). This chapter shows how to assemble them
into a `PyrolysisProblem` and, critically, how to select the *physics* that the
solver actually integrates. Pyrolysis.jl is a "generalized model with optional
physics": a single problem object can be solved with a fixed or moving mesh, with
or without pressure-driven gas flow, with one of several radiation closures, and
with or without a variable cross-section. These choices form the **simulation
mode**, and most of them are selected by keyword arguments to `solve` (or
automatically inferred from the material). This chapter explains every axis, its
default, and how the axes map onto the four residual paths the solver dispatches
to. The physics and derivations live in the Technical Reference; here we focus on
how to *configure* a run correctly.

---

## 7.1 The `PyrolysisProblem` object

A `PyrolysisProblem` couples a mesh, a material, boundary conditions, initial
conditions, and time-span metadata. It is the canonical user-facing entry point
and the only thing you pass to `solve`.

### 7.1.1 Constructor

The keyword constructor is:

```julia
PyrolysisProblem(;
    mesh,                       # Mesh1D from create_uniform_mesh / create_stretched_mesh
    material,                   # Material with components + reactions
    bc_set::BoundaryConditionSet,
    T_initial = z -> 300.0,     # Initial temperature closure  z -> T [K]
    ξ_initial = nothing,        # Initial concentration closures; auto-default if nothing
    tspan,                      # (t_start, t_end) [s]  — required
    dt_initial = 1e-4,          # Initial timestep hint [s]
    output_times = Float64[],   # Snapshot times [s]; auto-generated if empty
)
```

A complete, minimal example (single charring solid + one gas product, constant
heat flux at the top, adiabatic and impermeable at the bottom):

```julia
using Pyrolysis

# --- Material: one solid that decomposes into a gas ---
virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

decomp = Reaction(:pyrolysis, :virgin => :MMA;
                  A = 8.5e12, E = 188e3, h = 870e3)   # h > 0 = endothermic

pmma = Material(name = :PMMA,
                components = (virgin, gas),
                reactions  = (decomp,))

# --- Mesh: 1 cm thick, 100 cells, NC = 2 components ---
mesh = create_uniform_mesh(0.01, 100, Val(2))

# --- Boundary conditions ---
bc_set = BoundaryConditionSet(
    thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top    = DiffusiveMassBC(0.01),   # gas escapes to atmosphere
               bottom = ImpermeableBC()),
)

# --- Problem ---
problem = PyrolysisProblem(
    mesh      = mesh,
    material  = pmma,
    bc_set    = bc_set,
    T_initial = z -> 300.0,
    tspan     = (0.0, 100.0),
)
```

> The component count passed to the mesh (`Val(2)`) **must** equal the number of
> components in the material. The validator (§7.9) checks this.

### 7.1.2 Initial conditions

`T_initial` and `ξ_initial` are *closures* of the spatial coordinate `z`. They are
evaluated at each cell center when the problem is solved.

- `T_initial :: z -> T` returns the initial temperature [K]. Default
  `z -> 300.0`.
- `ξ_initial` is a `Vector` of closures, one per component, each returning the
  initial mass concentration `ξ_j` [kg/m³]. You may pass:
  - `nothing` (default): the constructor auto-fills the first component, **if it
    is a solid**, at its reference density `ρ(300.0)`, and sets all others to
    zero;
  - a single closure (broadcast to all components);
  - an explicit `Vector` of closures.

```julia
# Explicit, per-component initial concentrations:
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0,   # virgin at full density
                 z -> 0.0],     # gas starts empty
    tspan = (0.0, 100.0),
)
```

Several spatial-profile helpers are exported for non-uniform initial states:
`step_profile`, `linear_profile`, `gaussian_profile`, and `tanh_profile`. For
example, a warm front pre-heated near the top surface (`z = L`):

```julia
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    T_initial = z -> linear_profile(z; z_min = 0.0, z_max = 0.01,
                                       value_min = 300.0, value_max = 400.0),
    tspan = (0.0, 100.0),
)
```

The problem-level initial-condition step also populates the per-cell initial
dry-solid mass `m_dry,i` on the mesh — the denominator used by the
pyrolysis-progress state χ (§7.6). See Technical Reference §12 for the initial-
and boundary-condition treatment.

### 7.1.3 Time span and output times

- `tspan = (t_start, t_end)` in seconds. The validator requires
  `t_start < t_end` and `t_start ≥ 0`.
- `dt_initial` is an initial timestep hint passed to the integrator (default
  `1e-4` s). It must be positive and no larger than the time span.
- `output_times` is the vector of times at which the solution is sampled. If you
  leave it empty (the default), the constructor auto-generates between 11 and 101
  evenly spaced samples over `tspan`. To control output explicitly:

```julia
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    tspan = (0.0, 120.0),
    output_times = collect(0.0:10.0:120.0),   # every 10 s
)
```

By default `solve` saves at exactly these `output_times` (it passes them as its
`saveat`); you can override the save points with the `saveat` keyword of `solve`
(see Chapter 8).

---

## 7.2 Simulation modes: the four trait axes

Internally the residual is specialized on a single bundled trait,
`SimulationMode{Mesh, Rad, Trans, Chi}`, that carries four orthogonal axes. You
never construct this object directly in normal use — `solve` derives it from your
keyword arguments and from the material. Understanding the axes tells you exactly
which physics is active.

| Axis | Trait type | Possible values | Selected by | Default |
|---|---|---|---|---|
| Mesh motion | `MeshMode` | `Eulerian`, `ALE` | `use_ale` kwarg | `Eulerian` |
| Radiation | `RadiationMode` | `NoRadiation`, `SurfaceAbsorption`, `BeerLambert` | `radiation_model` kwarg | `NoRadiation` |
| Gas transport | `TransportMode` | `Fickian`, `DarcyFick` | material permeability (auto) | `Fickian` |
| Cross-section | `ChiMode` | `NoChi`, `WithChi` | material `lateral_shrinkage_law` (auto) | `NoChi` |

Two axes are chosen by you at solve time (mesh motion, radiation); two are
inferred automatically from the material (transport, cross-section). The trait
subtype names (`Eulerian`, `ALE`, `Fickian`, etc.) are *not* exported at the top
level; if you ever need to refer to them, import from `Pyrolysis.Problem`, e.g.

```julia
using Pyrolysis.Problem: ALE, BeerLambert, DarcyFick, WithChi
```

The corresponding compile-time accessors are `mesh_mode`, `radiation_mode`,
`transport_mode`, and `chi_mode`, applied to a `SimulationMode` instance. See
Technical Reference §15.1 for how the bundled trait drives the four residual
paths.

---

## 7.3 Mesh-motion mode: Eulerian vs. ALE (`use_ale`)

### 7.3.1 What it controls

- **Eulerian** (`use_ale = false`, default): the mesh is fixed. Node positions
  `z` are *not* part of the state vector; geometry is constant and is mirrored
  into the workspace once at construction. The state vector holds only
  temperatures and concentrations: `[T | ξ_1..ξ_NC]`.
- **ALE** (`use_ale = true`): the mesh moves. Node positions `z` become part of
  the state vector and evolve with a mesh velocity `w` derived from the
  reaction-driven volumetric strain rate `θ`. This lets the domain shrink as the
  solid is consumed (charring/ablation) or expand (intumescence) without a
  moving-boundary PDE. The state vector grows to `[T | ξ_1..ξ_NC | z]`.

```julia
# Eulerian (default): fixed mesh
sol = solve(problem)

# ALE: moving mesh that tracks recession / expansion
sol = solve(problem; use_ale = true)
```

### 7.3.2 When to use ALE

Use ALE when the through-thickness extent of the material changes appreciably —
charring wood, ablating polymers, intumescent coatings, or any case where surface
recession matters for the surface temperature and mass-loss rate. For a TGA-like
or thin-film case where geometry change is negligible, Eulerian is cheaper and
simpler.

ALE adds work arrays, mesh-velocity integration, advection corrections, and a
geometric-conservation-law-exact volume update; the geometry is recomputed from
the node positions at every residual evaluation. Bottom node 1 is pinned
(`w_1 = 0`): the substrate is immobile and the material moves as it shrinks or
swells. See Technical Reference §11 for the full ALE framework and §13.9 for mesh
generation that pairs well with it (stretching toward the surface).

> ALE is also what surface-depletion handling and thickness-termination callbacks
> build on; those are covered in Chapter 10.

---

## 7.4 Radiation mode (`radiation_model`)

The volumetric radiation closure is selected by the `radiation_model` keyword,
which takes a `RadiationModel` enum value. Three values are wired into the core
solver:

| Value | Meaning | Cost | Jacobian pattern |
|---|---|---|---|
| `NO_RADIATION` (default) | No volumetric absorption. External radiative *boundary* conditions still apply. | O(1) | tridiagonal |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the exposed surface via the boundary energy balance (opaque). Component `α` values are ignored. | O(1) | tridiagonal |
| `BEER_LAMBERT` | In-depth exponential absorption `I(z) = I_surface·exp(−∫α dz)`. Requires component absorption coefficients `α`. | O(n) | upper-triangular |

```julia
# Opaque material — absorb the cone flux at the surface:
sol = solve(problem; radiation_model = SURFACE_ABSORPTION)

# Semi-transparent material — absorb in depth (needs α on components):
sol = solve(problem; radiation_model = BEER_LAMBERT)
```

### 7.4.1 How radiation actually enters the energy equation

The *incident* radiative flux comes from a `RadiativeFluxBC` in your thermal
`bc_set` (e.g., a cone heater). The `radiation_model` keyword decides *where that
flux is deposited*:

- With `SURFACE_ABSORPTION`, the absorbed flux `α·Φ` (`α = ε_eff` by Kirchhoff,
  or the explicit `absorptivity`) enters at the exposed surface through the
  boundary energy balance, raising `T_s`; nothing is deposited volumetrically.
- With `BEER_LAMBERT`, the intensity is attenuated into depth according to the
  cell-local optical thickness `τ_i = α_i·Δz_i`, using each component's mass-basis
  absorption `α` (the effective volumetric absorption is `α_eff = Σ_j ξ_j α_j`).
- With `NO_RADIATION`, no volumetric source is added; a `RadiativeFluxBC` still
  acts — its absorbed flux `α·Φ` is applied at the surface exactly as under
  `SURFACE_ABSORPTION` — and `RadiativeBC` (re-radiation loss,
  `Fεσ(T_∞⁴ − T_s⁴)`) and `ConvectiveBC` act at the boundary as always.

`solve` re-samples the incident flux at every RHS evaluation whenever the model
consumes it (`SURFACE_ABSORPTION` / `BEER_LAMBERT`) — constant and time-varying
fluxes alike — so a ramped cone flux is always tracked exactly.

### 7.4.2 Diagnostics and warnings

`solve` validates your choice against the material:

- `SURFACE_ABSORPTION` with components that have non-zero `α` warns that those
  coefficients are ignored.
- `BEER_LAMBERT` with no `α` on any component warns that radiation passes through
  unabsorbed.

Choose the model by optical thickness `τ = α·L`: `SURFACE_ABSORPTION` for opaque
chars (`τ ≫ 1`), `BEER_LAMBERT` for semi-transparent polymers and foams, and
`NO_RADIATION` for transparent samples or when heating is purely conductive/
convective.

### 7.4.3 P1 is experimental

The fourth enum value, `P1_QUASI_STEADY`, is **not** wired into the core unified
solve path. It lives in `Pyrolysis.Experimental.RadiationP1` as a thin residual
wrapper: the ODE state is the ordinary unified layout, and the radiation field
`G` is a quasi-steady algebraic field re-solved into the radiation cache at
every residual evaluation — it is not part of the state vector.
Passing `radiation_model = P1_QUASI_STEADY` to `solve` raises an error — the core
solve path rejects it at keyword validation rather than silently degrading. Use the
Experimental primitives directly if you need it. The P1-only radiation boundary types (`MarshakRadiationBC`,
`ReflectiveRadiationBC`, `DirichletRadiationBC`) are likewise quarantined in
`Pyrolysis.Experimental`. See Technical Reference §8 for all four models and
§8.12 for the experimental status.

---

## 7.5 Gas-transport mode: Fickian vs. Darcy–Fick (automatic)

The transport axis is **not** a `solve` keyword — it is inferred from the
material. If any solid or liquid component has a non-zero permeability `κ`, the
material reports `has_darcy_flow == true` and the solver selects `DarcyFick`
(pressure-driven Darcy advection plus Fickian diffusion). Otherwise it selects
`Fickian` (diffusion only).

| Mode | Active when | Gas flux |
|---|---|---|
| `Fickian` (default) | all `κ = 0` | diffusive only: `J = −λ_face (ξ_R T_R − ξ_L T_L)/(T_face·d)` |
| `DarcyFick` | any solid/liquid `κ > 0` | advective + diffusive: `J = (ξ_upwind/φ_face)·v_g + J_diff`, with `v_g = −(κ/μ)∇P` (ξ is per bulk volume, v_g superficial — the advected scalar is the intrinsic density ξ/φ) |

### 7.5.1 Enabling Darcy–Fick

Put a permeability on the relevant solid/liquid components:

```julia
# Diffusion-only (Fickian): no κ → all gas transport is diffusive
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15)

# Darcy–Fick: assign permeability and the solver switches automatically
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, κ = 1e-15,
                        ρ_intrinsic = 1500.0)  # skeletal density → real pore volume
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, κ = 1e-12,
                        ρ_intrinsic = 1950.0)
gas    = GaseousComponent(:volatiles; M = 0.030, c = 1100.0, k = 0.02, λ = 1e-5)
air    = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)

wood = Material(name = :CharringWood,
                components = (virgin, char, gas, air),
                reactions  = (decomp,),
                mixing = MaterialMixing(k = MixingSpec(PARALLEL, 0.5),
                                        κ = MixingSpec(WEIGHTED, 0.5)))
```

> **Initialize the pore gas.** A Darcy–Fick problem needs (i) skeletal
> densities on the condensed components (`ρ_intrinsic`, or the `φ` kwarg) so a
> real pore volume exists, and (ii) an inert background air component
> initialized at ambient pressure with
> `background_gas_concentration(material, :air, ξ_condensed₀, T₀)` — otherwise
> the pore space starts as a vacuum against the 1 atm boundary and early
> product gas is spuriously advected inward. See Technical Reference §9.5.4 and
> the shipped `08_darcy_fick_demo.jl`.

Here the effective permeability is computed from the component values via the `κ`
mixing rule. As virgin converts to char, `κ` rises and bulk gas flow becomes more
important. The solver allocates the Darcy sub-cache only when `DarcyFick` is
selected, so a Fickian material pays nothing for the machinery.

Gas-only components carry no permeability (gases flow *through* the matrix); their
permeability is always zero. Only solids and liquids contribute.

### 7.5.2 Pressure boundary conditions

Darcy flow is driven by the ideal-gas pressure field and its gradient. To control
the pressure at a boundary, attach a pressure BC to the **mass** slot of the
boundary using the `+` operator. Pressure BCs are not a separate field of
`BoundaryConditionSet`; they ride along with the mass BC:

```julia
# Top surface: gas diffuses out AND is held at atmospheric pressure
top_mass = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

bc_set = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(75e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = top_mass, bottom = ImpermeableBC()),
)
```

The available pressure BCs are `DirichletPressureBC(P)` (fixed pressure),
`NeumannPressureBC(∂P/∂n)` (prescribed gradient), and `ConvectivePressureBC`
(Robin outflow). They set the boundary-face Darcy velocity, which vents gas
advectively (the venting flux scales with the interior concentration,
permeability, and pressure gradient); they contribute no diffusive film-law
flux of their own — that channel comes from the diffusive mass BC they are
combined with, and the two add. See Technical Reference §9.5.3 for the venting
chain and §12.5 for the pressure BCs.

> A Darcy–Fick material is almost always run with `use_ale = true`, because
> pressure-driven outflow accompanies recession; the shipped pressure-treated-wood
> and Douglas-fir examples do exactly this.

---

## 7.6 Variable cross-section mode: NoChi vs. WithChi (automatic)

The cross-section axis is also inferred from the material. If the material's
`lateral_shrinkage_law` is `nothing` (the default), the column cross-sectional
area is fixed at `A_0` and the mode is `NoChi`. If you supply a law — a function
of the mass-weighted column-average pyrolysis progress `χ̄` — the mode becomes
`WithChi` and the area evolves as `A(t) = A_0 · law(χ̄)`.

### 7.6.1 Defining a shrinkage / swelling law

```julia
# Linear lateral shrinkage: area shrinks to 80% as pyrolysis completes
shrink_law = χ̄ -> 1.0 - 0.2 * χ̄

# Linear lateral swelling (intumescence): area grows by 50%
swell_law  = χ̄ -> 1.0 + 0.5 * χ̄

wood = Material(name = :CharringWood,
                components = (virgin, char, gas),
                reactions  = (decomp,),
                lateral_shrinkage_law = shrink_law)
```

The progress variable `χ_i` in each cell tracks the cumulative dry-pyrolysis-gas
mass released per unit initial dry-solid mass; only "dry" gases (not water vapor)
contribute. `χ̄` is the mass-weighted column average. The law must be a smooth
function: its derivative is taken by automatic differentiation at run time, so
piecewise/non-smooth laws give undefined results. The default identity behavior
(no law) keeps the χ block out of the state vector entirely.

### 7.6.2 Interaction with ALE

`WithChi` is most meaningful together with `use_ale = true`: the radial
(cross-section) strain rate `θ_A = law'(χ̄)/law(χ̄) · dχ̄/dt` is subtracted from
each cell's volumetric strain to get the axial component that drives mesh motion,
so the two channels (radial area change and axial node motion) do not
double-count. With `WithChi`, the state vector carries the χ block as well:
`[T | ξ_1..ξ_NC | z | χ]`. The χ block is *reconstructed* from the dry-gas source
at every residual evaluation rather than being a free integrator state, which
keeps it consistent with the accumulated gas release.

```julia
problem = PyrolysisProblem(mesh = mesh, material = wood, bc_set = bc_set,
                           tspan = (0.0, 120.0))

# WithChi is auto-selected because `wood` carries a lateral_shrinkage_law:
sol = solve(problem; use_ale = true)

# The cross-section history is recorded in the solution extras:
sol.extras.cross_section_area_history
```

See Technical Reference §10 for volume change, swelling, and the lateral-shrinkage
law, and §11 for how `θ_A` couples into ALE mesh motion.

---

## 7.7 Energy form (`energy_form`) and the transpiration caveat

Separate from the four trait axes, a runtime option selects the bulk energy
formulation. It is passed as the `energy_form` keyword to `solve`:

| Value | Meaning |
|---|---|
| `:standard` (default) | Formulation A. The volumetric gas-advection source `S_conv = −Σ_g c_{p,g} J_g ∂T/∂z` is active; the gas-generation enthalpy-divergence term is dropped. This matches the FDS / Gpyro / ThermaKin convention. |
| `:with_generation_sink` | Adds the gas-generation enthalpy sink `S_gen = −Σ_g Δh_g(T_i,T_ref)(∇·J_g)_i` as an extra cell source — a *partial* move toward conservative form. **Requires T_ref-datum reaction heats** (see warning below). |

```julia
sol = solve(problem)                                   # :standard (default)
sol = solve(problem; energy_form = :with_generation_sink)
```

`energy_form` is a runtime field, not a trait — it adds one additive term and is
cheaper to branch on at run time than to compile as a separate method. Any value
other than these two raises an `ArgumentError` at problem construction.

> **Warning — `:with_generation_sink` is datum-dependent.** The sink is
> consistent only when every reaction heat `h` is referenced to the
> `T_ref ≈ 298.15 K` ambient datum (formation-enthalpy style). If your heats
> are DSC-style values quoted at the reaction temperature — the usual meaning
> of calibrated heats, including anything ported from ThermaKin or Gpyro —
> the gas products' sensible enthalpy is already inside `h`, and enabling the
> sink subtracts it a **second** time (≈0.7–1.2 MJ per kg of volatiles at
> 900 K, comparable to the heat of pyrolysis itself). In that case
> `:standard` is the exact form; leave the default alone. A one-time warning
> is emitted when the sink is selected (Technical Reference §6.6.3, §7.5.4).

### 7.7.1 Do not enable the surface transpiration BC with the default convention

There is a fourth `solve` keyword, `use_transpiration_bc` (default `false`),
which would fold a surface transpiration enthalpy term `Σ_g ṁ_g Δh_g(T_s, T_∞)`
into the boundary energy balance. **It is rejected at problem construction**:
setting `use_transpiration_bc = true` throws an `ArgumentError`. The reason is a
double-count. The volumetric advection source `S_conv` is *always* active (in
both energy forms); the surface transpiration term overlaps it by exactly
`∫S_conv`, producing a spurious sink of order 5–10 kW/m² that under-predicts
surface temperature and mass-loss rate. The supported convention is `S_conv` only
(`use_transpiration_bc = false`) with a plain convective + radiative boundary
condition. The transpiration path is retained for closed-ledger experiments only
and must never be combined with the volumetric source. See Technical Reference
§7.6 and §16.4 for the energy-balance accounting.

---

## 7.8 How modes map to the residual paths

The four trait axes plus the energy form fully determine the discrete RHS the
solver integrates. The two mesh-motion paths are the primary fork:

- **Eulerian path** — geometry is constant. The residual populates the typed
  state cache, updates properties, computes the physics RHS (dispatched on
  `Fickian` vs `DarcyFick`), and packs `dT` and `dξ`. No mesh velocity, no
  advection corrections.
- **ALE path** — geometry is recomputed from the node positions. The residual
  additionally updates the cross-section area (for `WithChi`), computes the
  per-cell strain `θ` and the column radial strain `θ_A`, integrates the
  cumulative mesh velocity `w` from the bottom up, applies ALE advection and the
  solid-dilation correction, packs the `z` block, and (for `WithChi`) packs the χ
  derivatives.

Radiation, transport, and cross-section are handled by compile-time dispatch
*within* those two paths. The net effect is a small, type-stable set of
specialized residuals selected entirely by the mode bundle — you do not choose a
"path" directly; you choose the axes and the path follows. See Technical Reference
§15.1–§15.3 for the dispatch and the full ALE assembly sequence.

---

## 7.9 Validating a problem (`validate_problem`)

Before solving, run the validator to catch structural mistakes early:

```julia
validate_problem(problem)                 # throws on any error, returns true if OK
validate_problem(problem; strict = true)  # adds expensive checks (see below)
```

`solve` calls `validate_problem(problem; strict = false)` automatically unless you
pass `validate = false`. The default checks are:

- the mesh has cells, nodes, faces, cell states, and initialized geometry;
- the material has at least one component, and every reaction references valid
  component indices;
- the mesh has the expected top (id 1) and bottom (id 2) boundaries;
- the time span is valid (`t_start < t_end`, `t_start ≥ 0`, positive `dt_initial`
  no larger than the span);
- the number of initial concentration closures equals the component count.

`strict = true` additionally verifies reaction stoichiometric mass balance
(product yields summing to 1.0), that the initial-condition closures are callable
at `z = 0`, and that all cell volumes are positive.

Two convenience helpers exist: `is_valid_problem(problem)` returns a `Bool`
without throwing, and `check_consistency(problem)` returns a NamedTuple
`(valid, messages)` describing component-count mismatches between mesh, material,
and initial conditions.

---

## 7.10 Choosing modes for common scenarios

The table below is a quick recipe lookup. "auto" means the axis is inferred from
the material, not set on `solve`.

| Scenario | `use_ale` | `radiation_model` | transport (auto) | cross-section (auto) | `energy_form` |
|---|---|---|---|---|---|
| TGA / DSC-like (thin, no recession) | `false` | `NO_RADIATION` | Fickian (`κ = 0`) | NoChi (no law) | `:standard` |
| Cone calorimeter, opaque char | `true` | `SURFACE_ABSORPTION` | Fickian or Darcy | NoChi or WithChi | `:standard` |
| Cone calorimeter, semi-transparent polymer | `true` | `BEER_LAMBERT` (set `α`) | Fickian | NoChi | `:standard` |
| Intumescent coating (swelling) | `true` | `SURFACE_ABSORPTION` | Darcy (`κ > 0`) | WithChi (swell law) | `:standard` |
| Pressure-driven char / pressure-treated wood | `true` | `SURFACE_ABSORPTION` | Darcy (`κ > 0`) | NoChi or WithChi | `:standard` |

### 7.10.1 Worked example: TGA-like, minimal mode

A thin sample with no geometry change, no radiation, diffusion-only gas escape:

```julia
virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
decomp = Reaction(:pyrolysis, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)
pmma   = Material(name = :PMMA, components = (virgin, gas), reactions = (decomp,))

mesh = create_uniform_mesh(1e-3, 40, Val(2))          # 1 mm, 40 cells

bc_set = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(t -> 1e3 * t),         # linear heating ramp
               bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = pmma, bc_set = bc_set,
                           tspan = (0.0, 600.0))

sol = solve(problem)        # Eulerian, NoRadiation, Fickian, NoChi, :standard
```

### 7.10.2 Worked example: cone calorimeter, opaque, recession + Darcy

A charring wood under a cone heater, with surface absorption, ALE recession, and
pressure-driven gas escape (Darcy auto-selected by the component permeabilities):

```julia
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, κ = 1e-15,
                        ρ_intrinsic = 1500.0)
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, κ = 1e-12,
                        ρ_intrinsic = 1950.0)
gas    = GaseousComponent(:volatiles; M = 0.030, c = 1100.0, k = 0.02, λ = 1e-5)
decomp = Reaction(:pyrolysis, :virgin => (:char, :gas);
                  A = 1e10, E = 150e3, h = 300e3, yields = (0.25, 0.75))

wood = Material(name = :CharringWood,
                components = (virgin, char, gas),
                reactions  = (decomp,),
                mixing = MaterialMixing(k = MixingSpec(PARALLEL, 0.5),
                                        κ = MixingSpec(WEIGHTED, 0.5)))

mesh = create_uniform_mesh(0.02, 80, Val(3))           # 2 cm, 80 cells, NC = 3

# Top: incident cone flux (surface-absorbed) + re-radiation + convection.
# Mass: diffusive escape + atmospheric pressure (drives Darcy outflow).
bc_set = BoundaryConditionSet(
    thermal = (top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.9) +
                     RadiativeBC(ε = 0.9, T_inf = 300.0) +
                     ConvectiveBC(h = 10.0, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
               bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = wood, bc_set = bc_set,
                           tspan = (0.0, 120.0))

sol = solve(problem;
            use_ale         = true,                 # recession
            radiation_model = SURFACE_ABSORPTION)   # opaque: deposit flux at surface
# transport = DarcyFick (auto, from κ); cross-section = NoChi (no law); energy_form = :standard
```

For a production Darcy run, also initialize the pore gas with a background air
component — see the pore-gas initialization note in §7.5.1.

To add lateral intumescent swelling, give `wood` a `lateral_shrinkage_law`
(e.g. `χ̄ -> 1.0 + 0.5*χ̄`) and the cross-section axis becomes `WithChi`
automatically — no change to the `solve` call is required.

---

## 7.11 Summary

- A `PyrolysisProblem` bundles mesh, material, `bc_set`, initial conditions, and
  `tspan`/`output_times`. Component counts must agree across mesh, material, and
  initial conditions.
- The active physics is the **simulation mode**: four trait axes plus the
  `energy_form` option.
  - *Mesh motion* (`use_ale`, default `false`) and *radiation*
    (`radiation_model`, default `NO_RADIATION`) are chosen at solve time.
  - *Gas transport* (Fickian vs. Darcy–Fick) and *cross-section* (NoChi vs.
    WithChi) are inferred automatically from the material — from component
    permeability `κ` and from `lateral_shrinkage_law`, respectively.
- `energy_form = :standard` (default) is the canonical convention; the surface
  transpiration BC is forbidden because it double-counts `S_conv`.
- `validate_problem` (run automatically by `solve` unless `validate = false`)
  catches structural errors before integration.

The next chapter, "Running the Solver," covers the rest of the `solve` keyword
surface: integrator and tolerance selection, Jacobian backends and linear
solvers, callbacks, warm starts, and output control.


---

# 8. Running the Solver

Once you have a `PyrolysisProblem` (Chapters 5–7), the single function that
integrates it forward in time is `solve`. This chapter is the complete option
reference for `solve(problem; ...)`: how to pick the time integrator and its
tolerances, how to choose and configure the structured Jacobian and the linear
solver behind it, how to switch on optional physics (ALE mesh motion, radiation,
the gas-generation energy sink), how callbacks are wired in, how warm-starting
and output saving work, and how to read mass-loss-rate output for variable
cross-sections. Every keyword shown here is a real `solve` keyword — `solve`
rejects unknown keywords with an error, so there is no room for typos. The
physics behind these options is derived in the Technical Reference; this chapter
tells you which knob to turn and what each one does. For the underlying temporal
integration and Jacobian theory, see Technical Reference §15.

---

## 8.1 The `solve` entry point

`solve` takes the problem as a positional argument and everything else as
keywords:

```julia
using Pyrolysis

solution = solve(problem)   # all defaults: KenCarp4 + structured Jacobian, no radiation
```

That call already does a great deal: it validates the problem, builds the
immutable problem definition and a mutable workspace, packs the initial state
vector, constructs the default structured Jacobian, assembles the callback set,
runs an implicit ODE integrator, and extracts a `PyrolysisSolution`. The
remaining sections explain how to override each part.

The full keyword signature, with defaults, is:

| Keyword | Default | What it controls |
|---|---|---|
| `integrator` | `nothing` → `KenCarp4` | ODE integrator (constructor or built instance); §8.2 |
| `linear_solver` | `nothing` → `SparspakFactorization()` | Sparse factorization paired with a constructor `integrator`; §8.4 |
| `abstol` | `1e-8` | Absolute error tolerance; §8.3 |
| `reltol` | `1e-6` | Relative error tolerance; §8.3 |
| `maxiters` | `1000000` | Maximum integrator steps |
| `dtmin` | `1e-12` | Minimum time step [s] |
| `dtmax` | `Inf` | Maximum time step [s] |
| `jacobian` | `nothing` → `Structured(DirectSolve())` | Pluggable Jacobian backend; §8.5 |
| `verify_jacobian` | `nothing` | Compare active Jacobian to a dense AD reference; §8.6 |
| `use_ale` | `false` | Enable ALE moving-mesh integration; §8.7 |
| `radiation_model` | `NO_RADIATION` | `NO_RADIATION` / `SURFACE_ABSORPTION` / `BEER_LAMBERT`; §8.8 |
| `energy_form` | `:standard` | `:standard` or `:with_generation_sink`; §8.9 |
| `use_transpiration_bc` | `false` | Surface transpiration enthalpy term — must stay `false`; §8.9 |
| `min_thickness` | `1e-6` | Terminate ALE solve when sample thins below this [m]; §8.7, §8.10 |
| `handle_depletion` | `false` | Merge depleted surface cells (ALE only); §8.10 |
| `depletion_threshold` | `0.05` | Cell-thickness fraction that triggers a merge; §8.10 |
| `min_cells` | `2` | Stop merging once this many cells remain; §8.10 |
| `diagnostics` | `false` | Snapshot the conservation ledger during the solve; §8.11 |
| `diagnostics_stride` | `10` | Record a diagnostics snapshot every N steps; §8.11 |
| `callback` | `nothing` | One extra user callback; §8.11 |
| `warm_start` | `nothing` | Reuse a prior solution's final state as the IC; §8.12 |
| `saveat` | `problem.output_times` | Times at which to store output; §8.13 |
| `sensealg` | `nothing` | Forwarded to the integrator (used by the sensitivity extension) |
| `progress` | `true` | Show a progress bar |
| `validate` | `true` | Run `validate_problem` before solving |

`solve` always returns a `PyrolysisSolution`; see Chapter 9 for its fields and
post-processing accessors.

---

## 8.2 Choosing the time integrator

Pyrolysis problems are **stiff**: Arrhenius reaction rates change by orders of
magnitude across a narrow temperature band, and conduction over fine cells
imposes a small explicit stability limit. The solver therefore defaults to an
implicit, stiffly stable method. When `integrator = nothing` (the default),
`solve` builds:

```julia
KenCarp4(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
```

`KenCarp4` is a 4th-order ESDIRK method — a good general-purpose choice for the
moderately stiff systems that pyrolysis produces. The `autodiff = AutoFiniteDiff()`
setting tells the integrator to finite-difference its *own* Jacobian when no
pluggable Jacobian is supplied; when you let the default structured Jacobian take
over (§8.5), that pluggable Jacobian is used instead.

You can pass the integrator in three forms:

**1. A constructor** (a type or function). `solve` instantiates it for you and
attaches the chosen linear solver:

```julia
using OrdinaryDiffEq            # solver types; `using Pyrolysis` already loads OrdinaryDiffEq

solve(problem; integrator = KenCarp4)              # rebuilt with default linsolve
solve(problem; integrator = TRBDF2)                # another SDIRK alternative
solve(problem; integrator = FBDF)                  # multistep BDF, good for very stiff TGA-like runs
```

**2. A constructor with a paired `linear_solver`.** Use this when you want a
specific sparse factorization but still want `solve` to wire it together:

```julia
solve(problem; integrator = KenCarp4,
               linear_solver = SparspakFactorization())
```

**3. A fully constructed algorithm instance.** Build it yourself, including the
linear solver, and pass it as `integrator`. In this case do **not** also pass
`linear_solver` — `solve` raises an error if you pair `linear_solver` with an
already-built instance:

```julia
solve(problem;
      integrator = KenCarp4(linsolve = SparspakFactorization()))
```

This third form is what the shipped Douglas-fir example uses
(`examples/07_douglas_fir_cone.jl`).

Recommended starting points:

| Situation | Integrator |
|---|---|
| Default cone / TGA / general use | `KenCarp4` (the default) |
| Very stiff, slow to converge | `FBDF` or `TRBDF2` |
| Tight accuracy near a sharp reaction front | `KenCarp4` with tighter `reltol` |

All integrators come from OrdinaryDiffEq.jl; any stiff solver from that package
works. See Technical Reference §15 for the implicit-RK formulation and the
`W = I/(γΔt) − J` matrix the integrator factors each step.

---

## 8.3 Tolerances and step-size limits

The integrator adapts its time step to keep the local error estimate within
`abstol` (absolute) and `reltol` (relative). The defaults are:

```julia
solve(problem; abstol = 1e-8, reltol = 1e-6)
```

Guidance:

- **Loosen** (`reltol = 1e-4`, `abstol = 1e-6`) for quick exploratory runs or
  parameter sweeps where you only need trends.
- **Tighten** (`reltol = 1e-7`, `abstol = 1e-9`) for convergence studies,
  Jacobian verification, or when comparing against reference codes. The
  Douglas-fir example uses `abstol = 1e-9, reltol = 1e-7`.
- Tolerances mix temperature (≈ 10²–10³ K) and concentration (≈ 0–10³ kg/m³) in
  one state vector, so prefer driving accuracy through `reltol` and keep
  `abstol` comfortably below the smallest concentration you care about.

Step-size limits:

```julia
solve(problem;
      dtmin = 1e-14,   # allow tiny steps through a rapid reaction front
      dtmax = Inf,     # no upper cap (default)
      maxiters = 2_000_000)
```

`dtmin` matters in ALE/charring runs: as cells thin, the diffusive stability
scale (`dt ∝ Δz²`) shrinks, and the integrator may need very small steps to push
through the front. If a solve stops with `retcode == :DtLessThanMin`, lower
`dtmin`. If it stops with `:MaxIters`, raise `maxiters` (or loosen tolerances).

---

## 8.4 The linear solver

Each implicit step solves a sparse linear system in the `W` matrix. The default
linear solver is `SparspakFactorization()`, a sparse direct factorization chosen
because it is fast on the block-tridiagonal sparsity that 1D pyrolysis produces.

You select the linear solver either by pairing it with a *constructor*
`integrator`:

```julia
using LinearSolve: SparspakFactorization, KLUFactorization, LUFactorization

solve(problem; integrator = KenCarp4,
               linear_solver = SparspakFactorization())   # default
solve(problem; integrator = KenCarp4,
               linear_solver = KLUFactorization())        # general-purpose sparse LU
```

or by baking it into a fully built integrator:

```julia
solve(problem; integrator = KenCarp4(linsolve = KLUFactorization()))
```

Sparspak and KLU are both interchangeable sparse direct solvers; Sparspak is the
default for the pyrolysis sparsity pattern. Iterative ILU-preconditioned solvers
for large ALE problems are covered in §8.14.

> Note: `solve` substitutes `LUFactorization()` (dense) — with a warning, and
> keeping your integrator type and `autodiff` setting — in two cases where a
> sparse solver cannot be used: (1) when no pluggable Jacobian is active
> (`handle_depletion = true`, see §8.10), and (2) when a dense Jacobian
> backend such as `DenseAD_ForwardDiff` is selected. In both cases there is
> no sparse Jacobian prototype for a sparse factorization to consume.

---

## 8.5 Jacobian backend selection

The single most important performance lever is the Jacobian. By default, `solve`
builds the **structured analytical Jacobian**:

```julia
Structured(strategy = DirectSolve())
```

This is selected automatically when `jacobian = nothing` (the default) for both
Eulerian and ALE problems. It assembles the exact Jacobian from per-operator
analytical and AD contributions (reaction, conduction, gas transport, radiation,
boundary conditions, and the ALE structured couplings) and hands a sparse matrix
to the linear solver. It is far faster than letting the integrator
finite-difference the whole RHS. See Technical Reference §15 for the
operator-based architecture.

You override it through the `jacobian` keyword, passing a `JacobianSpec`. The
spec is built from a *backend* and a *linear-solve strategy*. The constructor
`JacobianSpec(backend, problem; ...)` builds the sparsity patterns and AD caches
once; **do not** construct it inside a hot loop.

### Linear-solve strategies

| Strategy | Behavior | When to use |
|---|---|---|
| `DirectSolve()` | Materialize the full Jacobian (block + densified structured couplings) and factor it. Exact. | Default; small to typical problems |
| `StructuredSolve()` | Sparse LU on the local block plus Richardson iterations that apply the structured couplings as fast O(N) operators. Exact, converges in 2–5 iterations. | Typical ALE problems |
| `ApproxSolve()` | Drop the structured ALE couplings (Schur-style approximation); solve only the local block. Approximate but much cheaper. | Large ALE problems where Newton still converges |
| `MatrixFreeSolve()` | GMRES on the operator action without materializing. | Research / very large problems |

The verifier-only backend `DenseAD_ForwardDiff()` runs ForwardDiff over the
entire residual to produce a dense reference Jacobian. It is correct by
construction but slow — use it only for verification (§8.6), never for
production.

### Examples

Production default (no keyword needed):

```julia
solve(problem)                                  # Structured(DirectSolve())
```

Opt into the exact structured/Richardson solve:

```julia
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem)
solve(problem; jacobian = spec)
```

Opt into the approximate ALE solve for a large moving-mesh problem. Because the
spec must match the simulation mode, pass `use_ale = true` to **both** the
`JacobianSpec` constructor and `solve`:

```julia
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
solve(problem; jacobian = spec, use_ale = true)
```

On the shipped Douglas-fir ALE case, `ApproxSolve` drops the dense mesh-velocity
coupling and gives a large assembly speedup, with a solution that matches the
full Jacobian to within solver tolerance.

The mode keywords accepted by `JacobianSpec(backend, problem; ...)` mirror what
`solve` selects: `use_ale`, `radiation_model`, `radiation_intensity`,
`energy_form`, `use_transpiration_bc`, `min_thickness`, `min_density`. They
default to Eulerian / no radiation, so a plain Eulerian problem needs none of
them.

> Caveat: a pluggable `jacobian` is **incompatible** with
> `handle_depletion = true`. Depletion resizes the state vector mid-solve, which
> the immutable Jacobian cache cannot survive, so `solve` errors if you request
> both. With depletion enabled, `solve` falls back to the integrator's own
> colored finite-difference Jacobian. See §8.10 and Technical Reference §14.

---

## 8.6 Verifying the Jacobian during a solve

To check that your chosen backend agrees with the dense AD reference at runtime,
pass `verify_jacobian` a NamedTuple. At each listed sample time the solver
compares the active Jacobian to a `DenseAD_ForwardDiff()` reference and **warns**
(it does not throw) if the deviation exceeds tolerance:

```julia
solve(problem;
      verify_jacobian = (rtol = 1e-6, atol = 1e-10,
                         sample_times = [0.0, 1.0, 10.0]))
```

Allowed keys are exactly `:rtol`, `:atol`, and `:sample_times`; any other key is
rejected. `rtol` defaults to `1e-6` and `atol` to `1e-10` if omitted. Use this
when developing a custom backend or diagnosing a "Jacobian mismatch" symptom
(slow Newton convergence, unexpected step rejections).

Verification is skipped — with a warning — if no Jacobian spec is active, which
happens under `handle_depletion = true`. For an out-of-solve check at a single
state point, call `verify_jacobian_at` directly (see Technical Reference §15 and
Chapter 11).

---

## 8.7 ALE moving-mesh integration

By default the mesh is **Eulerian** (fixed). For charring, ablating, swelling, or
intumescent materials whose surface recedes or expands, enable the Arbitrary
Lagrangian–Eulerian (ALE) formulation, which adds the node positions `z` to the
state vector and moves the mesh with the material:

```julia
solve(problem; use_ale = true)
```

With `use_ale = true`:

- Node positions evolve as state variables, so the domain thickness
  `L = z[end] − z[1]` changes in time and is reported in
  `solution.extras.thickness_history`.
- A thickness-termination callback stops the integration when the sample thins
  below `min_thickness` (default `1e-6` m), preventing a singular `ρc_p → 0`. Set
  it to your physical stopping criterion, e.g. `min_thickness = 0.5e-3`.
- If the material carries a lateral-shrinkage law, the WithChi cross-section
  state `χ` is integrated too, and the cross-sectional area `A(t)` varies (see
  §8.13 and Technical Reference §10).

A typical opaque charring cone run:

```julia
solve(problem;
      use_ale = true,
      radiation_model = SURFACE_ABSORPTION,
      min_thickness = 0.5e-3,                # stop below 0.5 mm
      integrator = KenCarp4(linsolve = SparspakFactorization()))
```

For Eulerian problems (`use_ale = false`), the mesh never moves, the thickness is
constant, and `solution.z` is left empty. ALE node-motion theory and the
Geometric Conservation Law are in Technical Reference §11.

---

## 8.8 Radiation model selection

The `radiation_model` keyword selects how incident radiant flux (supplied by a
`RadiativeFluxBC` on the top boundary, see Chapter 5) is deposited into the
domain. The three core-path values are:

| Value | Physics | Cost / Jacobian |
|---|---|---|
| `NO_RADIATION` (default) | No volumetric radiative heating | O(1), tridiagonal |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the surface (opaque) | O(1), tridiagonal |
| `BEER_LAMBERT` | Exponential in-depth absorption `q''' = α·I` (semi-transparent) | O(n), non-local upper-triangular Jacobian |

```julia
solve(problem; radiation_model = SURFACE_ABSORPTION)   # opaque char, wood
solve(problem; radiation_model = BEER_LAMBERT)         # semi-transparent polymer
```

Notes and guards (all enforced by `solve`):

- **`SURFACE_ABSORPTION` with absorption coefficients set.** If your material
  components carry non-zero mass absorption `α`, `SURFACE_ABSORPTION` *ignores*
  them and applies the absorbed flux at the surface; `solve` warns. Use
  `BEER_LAMBERT` if you want the in-depth profile.
- **`BEER_LAMBERT` without absorption coefficients.** If no component has an `α`,
  `BEER_LAMBERT` lets all radiation pass through unabsorbed; `solve` warns. Set
  `α` (in m²/kg) on the components.
- **`BEER_LAMBERT` with a bottom `RadiativeFluxBC`.** The in-depth sweep models
  top-incident radiation only, and the boundary dispatch zeroes a
  `RadiativeFluxBC`'s face flux under `BEER_LAMBERT` — back-side heating would
  be silently lost, so `solve` raises an `ArgumentError`. Use
  `SURFACE_ABSORPTION` (face deposit) or a `HeatFluxBC` for back-side radiant
  heating.
- **`P1_QUASI_STEADY` is not available on the core path.** It lives in
  `Pyrolysis.Experimental.RadiationP1`, and `solve` raises an error
  if you pass it. Drive P1 through the Experimental primitives directly.

### Time-varying radiation flux is always tracked

Whenever `radiation_model` is `SURFACE_ABSORPTION` or `BEER_LAMBERT`, `solve`
re-evaluates the incident flux `I(t)` of the top `RadiativeFluxBC` at **every**
RHS call — there is no up-front constancy probe, so a ramped, delayed, or
modulated cone flux is tracked exactly no matter when it starts varying. The
cost is one call of your flux function per RHS evaluation, which is negligible.
Just make the `RadiativeFluxBC` flux a function of `t` for a time-dependent
source:

```julia
top_bc = RadiativeFluxBC(flux = t -> 25e3 + 250.0 * t, absorptivity = 0.9) +
         RadiativeBC(ε = 0.9, T_inf = 300) +
         ConvectiveBC(h = 10, T_inf = 300)
```

Re-radiation losses (`RadiativeBC`) and convective losses (`ConvectiveBC`) are
ordinary thermal boundary conditions and are applied regardless of
`radiation_model`. The radiation physics and the surface energy balance are
derived in Technical Reference §8 and §12.

---

## 8.9 Energy formulation and the transpiration caveat

The bulk energy equation uses **Formulation A** (the ThermaKin/Gpyro/FDS
convention): gas sensible storage is excluded from the effective heat capacity
and the gas enthalpy is carried by the volumetric advection source `S_conv`,
which is always active. The `energy_form` keyword selects whether an additional
gas-generation enthalpy sink is included:

```julia
solve(problem; energy_form = :standard)              # default: S_conv only
solve(problem; energy_form = :with_generation_sink)  # adds −Σ Δh_g (∇·J_g) cell source
```

- `:standard` (default) matches ThermaKin Eq. 17 — `S_conv` alone carries gas
  enthalpy. This is the recommended convention for almost all work.
- `:with_generation_sink` adds the per-cell gas-generation sink term, which
  telescopes with surface transpiration to recover the strict conservative form.
  Use it only if you specifically need the conservative-form energy balance to
  close (see Technical Reference §7 and §16).

Any other symbol is rejected.

### Do not enable `use_transpiration_bc`

The surface transpiration enthalpy term `Σṁ_g·Δh_g` would **double-count** gas
enthalpy because `S_conv` already carries it. `solve` therefore hard-rejects
`use_transpiration_bc = true`:

```julia
solve(problem; use_transpiration_bc = false)   # the only accepted value (default)
```

Leave it at its default `false`. Mass transport across the boundary remains fully
active either way — this flag concerns only the *enthalpy* sink, not the mass
flux. (Parameter sets calibrated against solvers that apply a surface transpiration
term *in addition to* a volumetric advection source embed the double-count; they
predict lower surface temperature and MLR than the `S_conv`-only balance and should
be re-validated when imported.)

---

## 8.10 Surface depletion and termination

For ALE charring runs where the surface cell shrinks to near-zero thickness, you
can merge depleted surface cells into their neighbors instead of dragging tiny
cells through the integration. This is the production depletion path, enabled
only together with `use_ale = true`:

```julia
solve(problem;
      use_ale = true,
      handle_depletion = true,    # merge thin surface cells
      depletion_threshold = 0.05, # merge a cell once it is < 5% of its initial thickness
      min_cells = 1,              # keep at least 1 cell; otherwise terminate
      min_thickness = 0.01e-3)    # terminate when total thickness < 0.01 mm
```

How it works:

- A discrete callback scans for **thin cells** — those whose current thickness is
  below `depletion_threshold ×` their **own** initial thickness (default `0.05`,
  range `(0,1)`). The reference is per-cell, so stretched meshes are judged
  cell-by-cell rather than against a single surface-cell thickness.
- A thin surface cell is merged into an interior neighbor with an
  **energy-conserving** temperature merge and a volume-weighted concentration
  merge, then the state vector is resized.
- Merging stops once only `min_cells` remain (default `2`); from there the
  thickness-termination callback ends the solve when total thickness drops below
  `min_thickness`.

The example `examples/06_cone_calorimeter.jl` uses exactly this configuration
(`handle_depletion = true`, `depletion_threshold = 0.05`, `min_cells = 1`). An
alternative is to *not* merge and instead let ALE carry the thin cells, stopping
purely on `min_thickness` — that is what `examples/07_douglas_fir_cone.jl` does
(`handle_depletion = false`).

> Important interaction with the Jacobian: because depletion resizes the state
> vector via `resize!`, it is incompatible with a pluggable `jacobian`. `solve`
> errors if you pass both. With `handle_depletion = true` the solver
> automatically uses the integrator's colored finite-difference Jacobian and a
> dense (resize-safe) LU linear solver; sparse factorizations passed via
> `linear_solver` are swapped to dense LU with a warning in this mode.
> The MKL pivot-array resize bug can only bite if you pass a fully built
> integrator with its `linsolve` left unset, letting LinearSolve auto-select
> MKL LU; `solve` warns about that combination on large depleting states.
> Depletion merge theory is in Technical Reference §14.

---

## 8.11 Callbacks and runtime diagnostics

`solve` assembles a `CallbackSet` internally; you do not manage most of these. In
order, they are: a step-update callback (always present; syncs the integrator
state into the mutable mesh), a thickness-termination callback (ALE), a depletion
callback (ALE + `handle_depletion`), a progress bar (`progress = true`), an
optional diagnostics callback, an optional Jacobian-verification callback, and
finally any user callback you supply.

### Conservation diagnostics

Set `diagnostics = true` to snapshot the conservation ledger during the solve.
Snapshots are recorded every `diagnostics_stride` steps (default `10`) and are
appended to `solution.extras.diagnostics_log`:

```julia
solution = solve(problem; diagnostics = true, diagnostics_stride = 20)
rows = solution.extras.diagnostics_log    # Vector{DiagnosticsRow}
```

Each row captures the matrix sensible energy `E_matrix`, total energy `E_total`,
and the gas-storage gap. This is **snapshot** logging — not a continuous time
history — and is intended for auditing energy/mass balance, not dense
post-processing. `diagnostics_stride` must be ≥ 1. See Chapter 9 for reading
diagnostics and Technical Reference §16 for the ledgers.

### A custom user callback

You may pass one additional callback (any OrdinaryDiffEq callback type) via
`callback`. It runs last in the set:

```julia
using DiffEqCallbacks  # for the callback types

stop_at_T = ContinuousCallback(
    (u, t, integ) -> u[1] - 700.0,        # back-face cell (index 1) reaches 700 K
    integ -> terminate!(integ))

solve(problem; callback = stop_at_T)
```

For automated convergence checks against the dense AD reference, use
`verify_jacobian` (§8.6) rather than a hand-rolled callback.

### Progress bar

`progress = true` (default) shows a percent-of-time-span progress bar. Set
`progress = false` for batch runs and parameter sweeps (the sweep harness does
this automatically).

---

## 8.12 Warm-starting from a prior solution

To start a new solve from the final state of a previous one — useful for
continuation studies and for warm-starting nearby parameter-sweep points — pass
the prior `PyrolysisSolution` as `warm_start`:

```julia
sol1 = solve(problem1)
sol2 = solve(problem2; warm_start = sol1)
```

The prior solution's final raw state vector overwrites the new initial
condition. Compatibility is enforced by a state-vector length check, which
implicitly requires matching component count, mesh size, ALE mode, and χ mode;
mismatches raise a clear error (and the prior solve must have stored its raw
solution). Note that `tspan[1]` of the *new* problem still sets the integrator's
start time — warm-starting treats the prior final state as fresh initial data, it
does not continue the clock.

> Caution: warm-starting can introduce subtle artifacts in sensitivity studies
> because the perturbed solve no longer begins from the canonical initial
> condition. Prefer cold starts when computing gradients (Chapter 13).

---

## 8.13 Output times, saving, and MLR for variable cross-section

By default `solve` saves the solution at `problem.output_times` (set when you
built the problem). Override with `saveat`:

```julia
# Save every second over a 600 s run
solve(problem; saveat = collect(0.0:1.0:600.0))
```

`saveat` accepts a vector of times or a step. Fewer save points reduce
interpolation overhead and memory; for the raw dense interpolant use
`interpolate_solution(solution, t)` (Chapter 9).

> ALE + depletion exception: when `use_ale = true` and `handle_depletion = true`,
> the integrator cannot interpolate across a state-vector resize, so `solve`
> internally switches to "save every step" and then filters down to your
> requested `saveat` times during extraction. You still get output near the times
> you asked for, plus the final time.

### Mass-loss-rate gauge normalization

For variable cross-section (WithChi) runs, the mass-loss rate (top-face flux
plus any back-face venting; the bottom share is zero for the default
impermeable holder) is reported two ways in the solution:

- `solution.mass_loss_rate[i]` — the **local** MLR per *current* area `A(t)`,
  in kg·m⁻²·s⁻¹ (both faces; the bottom term is zero for an impermeable
  holder).
- `solution.extras.MLR_gauge[i]` — the **gauge** MLR normalized to the *initial*
  area `A_0`, i.e. `MLR_local × A(t)/A_0`. This is the cone-calorimeter
  convention, since a cone gauge measures mass loss per fixed sample footprint.

For constant cross-section (no lateral-shrinkage law) the two coincide because
`A(t) = A_0`. Use `MLR_gauge` when comparing to cone-calorimeter data on a
shrinking or swelling sample. The total domain mass in
`solution.extras.total_mass` is likewise normalized by the initial area. See
Technical Reference §10 and §16 for the gauge convention.

---

## 8.14 Performance tips and iterative linear solvers

A short checklist for fast, reliable solves:

1. **Use the default structured Jacobian.** It is dramatically faster than the
   integrator's finite-difference Jacobian. Only fall back to FD when forced
   (depletion).
2. **Right-size the mesh.** Resolution costs scale super-linearly (smaller cells
   force smaller diffusion-limited time steps). Use a stretched mesh refined near
   the surface (`create_stretched_mesh(L, n, Val(NC), 1.1; stretch_from = :top)`)
   instead of a uniformly fine mesh; see Chapter 6 and the
   `convergence_douglas_fir` example.
3. **For large ALE problems, try `ApproxSolve`.** Dropping the dense ALE
   couplings cuts assembly cost substantially when Newton still converges (§8.5).
   Verify with `verify_jacobian` first.
4. **Loosen tolerances for sweeps.** `reltol = 1e-4` is plenty for trend studies.
5. **Disable the progress bar** (`progress = false`) and use small `extract`
   closures in parameter sweeps (Chapter 12).

### ILU-preconditioned iterative solvers

For very large ALE problems where direct factorization dominates the cost, two
ILU-preconditioned Krylov solvers ship with the package and plug in as the
integrator's linear solver:

```julia
using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization

# ILU-preconditioned GMRES
solve(problem;
      integrator = KenCarp4(linsolve = ILUGMRESFactorization(
          τ = 0.1,            # ILU drop tolerance
          gmres_rtol = 1e-8,  # GMRES relative tolerance
          restart = 30)))     # GMRES restart frequency

# ILU-preconditioned BiCGSTAB (fixed memory, no restart)
solve(problem;
      integrator = KenCarp4(linsolve = ILUBiCGSTABFactorization(
          τ = 0.1, rtol = 1e-8)))
```

`ILUGMRESFactorization` and `ILUBiCGSTABFactorization` live in the `Solver`
submodule — bring them in with
`using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization`
(or `using Pyrolysis.Internal`), or qualify as `Pyrolysis.ILUGMRESFactorization`.
If the ILU factorization becomes singular during a sharp transient, both solvers
fall back to diagonal preconditioning automatically.

> Trade-off: inexact iterative solves can cause more Newton-step rejections on
> stiff Arrhenius transients. For typical 1D pyrolysis problems the direct
> Sparspak solver is usually faster end-to-end despite being "more expensive" per
> factorization. `solve` warns when you select ILU-GMRES, recommending Sparspak
> unless the problem is genuinely factorization-bound. Reserve the iterative
> solvers for large ALE meshes where the dense structured coupling makes direct
> factorization the bottleneck.

External factorizations are available through extensions: `using KLU` enables
`Pyrolysis.klu_factorization()`, `using Pardiso` enables
`Pyrolysis.pardiso_factorization()`, and `using AlgebraicMultigrid` enables
`Pyrolysis.amg_preconditioner(A)`. These are covered in Chapter 11. The
structured-solve strategies and ILU internals are described in Technical
Reference §15.

---

## 8.15 Reading the return code

`solve` always returns a `PyrolysisSolution`; check `solution.retcode` to confirm
the run finished as intended:

```julia
solution = solve(problem)
println("retcode    = ", solution.retcode)   # e.g. :Success
println("solve time = ", solution.solve_time, " s")
```

Common return codes (from OrdinaryDiffEq):

| `retcode` | Meaning | Typical fix |
|---|---|---|
| `:Success` | Integrated to `tspan[2]` cleanly | — |
| `:Terminated` | A callback stopped it (thickness/depletion/user) | Expected for charring; widen `min_thickness` if premature |
| `:MaxIters` | Hit `maxiters` before the end | Raise `maxiters` or loosen tolerances |
| `:DtLessThanMin` | Step fell below `dtmin` | Lower `dtmin`; check for unphysical states / over-fine mesh |

A `:Terminated` code is normal and expected for ALE charring runs that deplete
the sample. For `:MaxIters` or `:DtLessThanMin`, see the troubleshooting guidance
in Chapter 15.

---

## 8.16 A complete annotated `solve` call

Putting the pieces together, here is a fully specified solve for an opaque
charring cone-calorimeter case with ALE mesh motion (mirroring the shipped
Douglas-fir example):

```julia
using Pyrolysis
using LinearSolve: SparspakFactorization

solution = solve(
    problem;
    # accuracy
    abstol      = 1e-9,
    reltol      = 1e-7,
    dtmin       = 1e-14,                    # allow tiny steps through the reaction front
    # integrator + linear solver (built-instance form)
    integrator  = KenCarp4(linsolve = SparspakFactorization()),
    # physics
    use_ale          = true,               # moving mesh: surface recedes
    radiation_model  = SURFACE_ABSORPTION, # opaque: deposit all flux at the surface
    energy_form      = :standard,          # S_conv-only (ThermaKin) convention
    min_thickness    = 0.5e-3,             # stop when the sample is < 0.5 mm
    handle_depletion = false,              # let ALE carry thin cells (no merging)
    # diagnostics + output
    diagnostics = true,                    # snapshot the energy ledger
    saveat      = collect(0.0:1.0:600.0),  # save once per second
    progress    = true,
)

@assert solution.retcode in (:Success, :Terminated)
```

The default `jacobian = nothing` gives the structured `DirectSolve` Jacobian
appropriate for this ALE problem; to opt into the approximate solve, add
`jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)`.

From here, Chapter 9 covers the `PyrolysisSolution` object — extracting
temperature and concentration profiles, surface and back-face temperatures,
mass-loss rate, and plotting.


---

# 9. The Solution Object and Plotting

Every call to `solve` returns a single `PyrolysisSolution`. This chapter is a complete
reference for that object: the concrete fields it always carries, the typed `extras`
NamedTuple of opt-in diagnostics, the accessor functions for pulling out profiles and
interpolating at arbitrary times, and the `Plots.jl` recipes shipped with the package.
It closes with practical post-processing recipes — computing surface temperature and
mass-loss rate, reading the conservation ledger, handling moving (ALE) meshes, and
preparing curves for comparison against cone-calorimeter or TGA data. The numerics that
produce these quantities are derived in the Technical Reference (energy assembly in
Technical Reference §7, the boundary-temperature solve in §12, conservation diagnostics
in §16).

A note on coordinates that governs every array below: `z = 0` is the **bottom /
substrate** (boundary id 2, tag `:bottom`) and `z = L` is the **top / exposed surface**
(boundary id 1, tag `:top`). Heat enters from the top and the material recedes downward.
"Surface" always means the top, exposed face; "back face" means the bottom.

---

## 9.1 What `solve` returns

`solve(problem; ...) -> PyrolysisSolution`. The return type is a single mutable struct.
You never construct it yourself — `solve` allocates it, populates the arrays during
solution extraction, and hands it back.

```julia
using Pyrolysis

solution = solve(problem; abstol = 1e-8, reltol = 1e-6)

solution                       # PyrolysisSolution(:Success, 601 timesteps, 12.4s)
typeof(solution)               # PyrolysisSolution{...}
```

The struct is defined in `src/Problem/pyrolysis_problem.jl`. Its layout is a **fixed core
schema** of concretely typed fields plus one typed `extras` NamedTuple. The fixed core is
what the plotting recipes and accessors rely on; `extras` holds secondary diagnostics that
are always present (the keys never change) but whose contents depend on the run.

`PyrolysisProblem`, `PyrolysisSolution`, and `solve` are all top-level exports of
`Pyrolysis`. Many of the convenience accessors discussed in §9.4 are *not* — they live in
submodules. The distinction matters for which `using` line you write, and is called out
explicitly where it applies.

---

## 9.2 Core fields

These fields are always present after a successful `solve`. They are populated for every
saved output time (`solution.t`), so an index `i` into the time vector selects a consistent
snapshot across all of them.

| Field | Type | Shape | Meaning / units |
|---|---|---|---|
| `problem` | `PyrolysisProblem` | scalar | The original problem (mesh, material, BCs) |
| `t` | `Vector{Float64}` | `(n_steps,)` | Saved solution times [s] |
| `T` | `Matrix{Float64}` | `(n_cells, n_steps)` | Cell-center temperatures [K] |
| `ξ` | `Array{Float64,3}` | `(NC, n_cells, n_steps)` | Mass concentrations `ξ_j` [kg·m⁻³] |
| `z` | `Matrix{Float64}` | `(n_nodes, n_steps)` | Node positions [m]; `0×0` in Eulerian mode |
| `surface_T` | `Vector{Float64}` | `(n_steps,)` | Top surface temperature `T_s` [K] |
| `mass_loss_rate` | `Vector{Float64}` | `(n_steps,)` | MLR per current area (top + permeable-bottom) [kg·m⁻²·s⁻¹] |
| `energy_residual` | `Vector{Float64}` | `(n_steps,)` | Energy-conservation residual; zeros if tracker off |
| `retcode` | `Symbol` | scalar | ODE-solver return code |
| `solve_time` | `Float64` | scalar | Wall-clock solve time [s] |
| `raw_solution` | `Any` | — | Raw `DifferentialEquations.jl` solution (interpolation source) |
| `extras` | `NamedTuple` | — | Opt-in diagnostics (see §9.3) |

A few points worth internalizing:

- **`T` and `ξ` are cell-centered.** Row `j` is cell `j`, counting from the bottom
  (`z = 0`) toward the top (`z = L`). So `solution.T[end, i]` is the cell nearest the
  exposed surface and `solution.T[1, i]` is the cell at the substrate. The actual *surface*
  temperature (on the exposed face, not the adjacent cell center) is the separately solved
  `solution.surface_T[i]` — see §9.5.

- **`ξ` is component-major in its first index.** `solution.ξ[k, j, i]` is the
  concentration of component `k` in cell `j` at saved time `i`. The component order matches
  the order you listed components in the `Material` constructor. These are *mass
  concentrations* `ξ_j = Y_j ρ` in kg·m⁻³, the primary species state variable — not mass
  fractions (Technical Reference §3, §4).

- **`z` is empty in Eulerian mode.** When you solve without `use_ale=true`, the mesh is
  fixed, so node positions are constant and `solution.z` is a `0×0` matrix. Use
  `get_cell_positions(solution)` (§9.4) to get the fixed cell centers instead. In ALE mode
  `z` is `(n_nodes, n_steps)`; see §9.7.

- **`retcode`** is `Symbol(raw_solution.retcode)`. Common values:
  `:Success` (completed to `t_end`), `:Terminated` (a callback stopped it — e.g. the ALE
  thickness-termination callback fired because the sample fully depleted),
  `:MaxIters`, `:DtLessThanMin`.

- **`energy_residual`** is all zeros unless a conservation tracker wrote into it. It is the
  series the `:residual` plot recipe consumes (§9.6). For the richer per-snapshot ledger,
  enable `diagnostics=true` and read `extras.diagnostics_log` (§9.8).

### Direct field access vs. accessors

The shipped examples read fields directly, which is the most robust path because the field
names are part of the fixed schema:

```julia
times      = solution.t                 # Vector{Float64}
T_field    = solution.T                 # (n_cells, n_steps)
T_at_t5    = solution.T[:, 5]           # profile at saved time index 5
T_surface  = solution.surface_T         # Vector{Float64}
mlr        = solution.mass_loss_rate    # Vector{Float64}
xi_virgin  = solution.ξ[1, :, end]      # component 1 profile at last time
```

There is also a small set of accessor functions (§9.4) — but be careful which package
namespace exports them.

### Display and quick checks

`PyrolysisSolution` has a pretty `show`, and a couple of boolean helpers:

```julia
solution                     # one-line summary
display(solution)            # multi-line summary (retcode, time range, peaks)

solution.retcode == :Success # success check that needs no import
```

Two helper functions, `succeeded(sol)` (true iff `retcode == :Success`) and
`n_timesteps(sol)` (`== length(sol.t)`), are exported by the `Problem` submodule, **not** by
top-level `Pyrolysis`. With plain `using Pyrolysis` you must qualify them
(`Pyrolysis.succeeded(solution)`) or pull them in via `using Pyrolysis.Internal`. Reading
`solution.retcode` directly is the simplest success check and needs no import.

---

## 9.3 The `extras` NamedTuple

`solution.extras` is a typed NamedTuple whose **keys are fixed at construction** — they are
always all present, even if a particular run leaves some of them at their default (zeros, or
an empty `diagnostics_log`). The numeric vectors are sized `(n_steps,)` and indexed by the
same time index as the core fields. This is where most cone-calorimeter post-processing
quantities live.

| `extras` key | Units | Meaning |
|---|---|---|
| `gas_generation_rate` | kg·m⁻²·s⁻¹ | Volumetric gas production integrated over the column, gauge-normalized |
| `back_temperature` | K | Back-face (bottom) surface temperature at `z=0` |
| `total_mass` | kg·m⁻² | Total domain mass per initial cross-section area |
| `thickness_history` | m | Sample thickness `z[end]-z[1]`; constant in Eulerian mode |
| `surface_mass_flux` | kg·m⁻²·s⁻¹ | Per-species summed **top-face** mass flux (excludes any bottom venting; cf. `mass_loss_rate`) |
| `surface_heat_flux_net` | W·m⁻² | Heat flux the boundary balance applies at the top face (`net = conduction + transpiration`); under `BEER_LAMBERT` it excludes the incident radiation, which is reported in `volumetric_absorption` |
| `surface_heat_flux_absorbed` | W·m⁻² | Absorbed incident (radiative) flux |
| `surface_heat_flux_reradiation` | W·m⁻² | Re-radiation loss `εσ(T_∞⁴−T_s⁴)` |
| `surface_heat_flux_convection` | W·m⁻² | Convective exchange `h_conv(T_∞−T_s)` |
| `surface_heat_flux_conduction` | W·m⁻² | Heat conducted into the slab at the surface |
| `surface_heat_flux_transpiration` | W·m⁻² | Enthalpy carried by escaping gas (zero with default S_conv) |
| `volumetric_absorption` | W·m⁻² | Beer–Lambert in-depth absorbed flux `α_surf·q_ext·(1−e^{−τ_tot})`; zero unless `radiation_model = BEER_LAMBERT` |
| `surface_emissivity` | – | Composition-weighted `ε_eff` at the surface |
| `surface_conductivity` | W·m⁻¹·K⁻¹ | Effective surface conductivity `k_eff` |
| `back_heat_flux` | W·m⁻² | Heat flux at the back face |
| `MLR_gauge` | kg·m⁻²·s⁻¹ | MLR normalized to initial area `A_0` (cone gauge convention) |
| `cross_section_area_history` | m² | Cross-section `A(t)` history (varies only in WithChi mode) |
| `diagnostics_log` | — | `Vector{DiagnosticsRow}`; empty unless `diagnostics=true` (§9.8) |

```julia
gen_rate  = solution.extras.gas_generation_rate    # kg/(m²·s)
T_back    = solution.extras.back_temperature       # K
mass_m2   = solution.extras.total_mass             # kg/m²
thickness = solution.extras.thickness_history      # m
q_net     = solution.extras.surface_heat_flux_net  # W/m²
```

There is one exported accessor for `extras`:

```julia
# Exported by the Problem submodule; qualify or `using Pyrolysis.Internal`.
Pyrolysis.get_gas_generation_rate(solution)   # == solution.extras.gas_generation_rate
```

Note its docstring caveat: `gas_generation_rate` is the *volumetric* production summed over
the domain, **not** the surface mass-loss rate. The canonical MLR is
`solution.mass_loss_rate` (per current area) and `extras.MLR_gauge` (per initial area). See
§9.5 for which one to use.

---

## 9.4 Accessor functions and interpolation

The package ships a handful of convenience functions. **None of them are top-level exports
of `Pyrolysis`.** `get_profile`, `get_cell_positions`, and `interpolate_solution` are
exported by the `Solver` submodule; `get_surface_temperature` and
`get_back_face_temperature` by the `Discretization` submodule. To call them unqualified,
bring them in via the umbrella namespace:

```julia
using Pyrolysis
using Pyrolysis.Internal     # re-exports every submodule's public symbol
```

or qualify each call (`Pyrolysis.get_profile(...)`). The umbrella namespace is documented as
the "escape hatch" for advanced users; symbols reached through it are not part of the stable
public API. For production scripts that only ever need the fixed fields, plain field access
(§9.2) avoids the issue entirely.

### `get_profile(solution, field, time_idx)`

Pulls a single field at one saved time index. `field` is `:T` for temperature, `:ξ` for the
full `(NC, n_cells)` concentration slab, or `:ξ1`, `:ξ2`, … (equivalently `:xi1`, `:xi2`,
…) for an individual component.

```julia
T_prof   = get_profile(solution, :T, 10)     # == solution.T[:, 10]
xi_all   = get_profile(solution, :ξ, 10)     # == solution.ξ[:, :, 10]  (NC × n_cells)
xi_char  = get_profile(solution, :ξ2, 10)    # == solution.ξ[2, :, 10]
```

`time_idx` is an index into `solution.t`, not a time in seconds. To find the index nearest a
target time, use `argmin(abs.(solution.t .- t_target))`.

### `get_cell_positions(solution)`

Returns the cell-center positions from the problem's mesh
(`solution.problem.mesh.cell_centers`). This is the right depth axis for plotting Eulerian
profiles, where the mesh is fixed and `solution.z` is empty.

```julia
z_centers = get_cell_positions(solution)     # Vector, length n_cells [m]
```

In ALE mode the mesh moves, so for moving-mesh plots you should reconstruct positions per
time step from `solution.z` instead (§9.7).

### `interpolate_solution(solution, t)` and the call syntax

Both interpolate the **raw state vector** at an arbitrary continuous time using the
underlying `DifferentialEquations.jl` dense interpolant. The state vector is packed
block-major `[T | ξ | z (ALE) | χ (WithChi)]` (Technical Reference §13, §15) — it is *not*
the convenient `(n_cells, n_steps)` view.

```julia
u_t = interpolate_solution(solution, 42.7)   # raw state vector at t = 42.7 s
u_t = solution(42.7)                          # identical: PyrolysisSolution is callable
```

Both error if `raw_solution === nothing` (e.g. a stub solution that was never solved). To
turn a raw state vector into physical quantities, use the `StateLayout` accessors. These
*are* top-level exports of `Pyrolysis` (and are also re-exported by the `Problem` submodule),
so plain `using Pyrolysis` is enough to call them unqualified:
`get_temperature`, `get_concentration`, and `get_node_positions`.

Mind their signatures — they are **per-cell / per-component scalar accessors**, not bulk
extractors:

- `get_temperature(u, layout, i)` returns the temperature of a *single* cell `i` (one scalar).
- `get_concentration(u, layout, i, j)` returns the concentration of component `j` in a
  *single* cell `i` (one scalar). Note the order: cell index first, component index second.
- `get_node_positions(u, layout)` is the exception — it returns a view of *all* node
  positions (ALE layouts only).

To rebuild a whole profile, loop over the cells. The `layout` comes from a freshly built
workspace:

```julia
def     = to_problem_def(solution.problem)
ws      = build_workspace(solution.problem, def)
layout  = ws.layout
n_cells = layout.n_cells

T_cells = [get_temperature(u_t, layout, i) for i in 1:n_cells]      # cell temperatures [K]
xi_1    = [get_concentration(u_t, layout, i, 1) for i in 1:n_cells] # component 1 [kg/m³]

# A single cell, if that is all you need:
T_cell1 = get_temperature(u_t, layout, 1)          # temperature of cell 1 [K]
```

For most workflows you do not need any of this — saving at the times you care about
(`output_times` / `saveat`) and indexing the core arrays is simpler and avoids rebuilding a
layout. Reach for interpolation, and these scalar accessors, only when you need a value
strictly between saved times.

### `get_surface_temperature(mesh)` / `get_back_face_temperature(mesh)`

These take a **mesh**, not a solution, and read the last-solved boundary state
(`mesh.boundary_states`). They reflect the mesh's state at the *end* of the solve (the
step-update callback syncs the mesh to the final integrator state), so they are a snapshot,
not a time history. For the time history of the surface temperature, use
`solution.surface_T` (top) and `solution.extras.back_temperature` (bottom), which are solved
at every saved time during extraction.

```julia
using Pyrolysis.Internal: get_surface_temperature, get_back_face_temperature
mesh = solution.problem.mesh
T_s_final    = get_surface_temperature(mesh)       # top face, final state [K]
T_back_final = get_back_face_temperature(mesh)     # bottom face, final state [K]
```

---

## 9.5 Mass-loss rate (MLR) and gas generation

Three related but distinct quantities are available. Choosing correctly matters when the
cross-section changes (variable-`A` / WithChi mode) or for cone-gauge comparison.

1. **`solution.mass_loss_rate`** — MLR per *current* area `A(t)`, [kg·m⁻²·s⁻¹]. It is the
   per-species-summed gas flux leaving the slab — the exposed face plus, when a permeable
   mass BC is configured there, the back face (zero for the default impermeable holder).
   The strictly top-face flux is `solution.extras.surface_mass_flux`. In constant-area runs
   (Eulerian, or ALE without a lateral-shrinkage law), this is the MLR.

2. **`solution.extras.MLR_gauge`** — surface MLR normalized to the *initial* area `A_0`, the
   cone-calorimeter "gauge" convention: `MLR_gauge = mass_loss_rate · A(t)/A_0`. With a
   variable cross-section this differs from the local MLR; with constant area the two
   coincide. **For cone-calorimeter comparison, use `MLR_gauge`** (the lab measures mass per
   fixed sample-holder area), unless your area is constant — in which case either works.

3. **`solution.extras.gas_generation_rate`** — *volumetric* gas production integrated over
   the whole column and divided by `A_0`, [kg·m⁻²·s⁻¹]. This is the total gas being
   *created* by reactions, not the gas crossing the surface. It can differ transiently from
   the surface flux because gas takes time to diffuse/advect out and can accumulate in the
   pore space. Use it to diagnose where in the reaction the production peaks; use the surface
   flux quantities for what actually leaves the sample.

```julia
times     = solution.t
mlr_local = solution.mass_loss_rate          # per A(t)
mlr_gauge = solution.extras.MLR_gauge        # per A_0 (cone convention)
gen_rate  = solution.extras.gas_generation_rate

# Peak MLR and the time it occurs:
i_peak    = argmax(mlr_gauge)
t_peak    = times[i_peak]
peak_mlr  = mlr_gauge[i_peak]
```

### An HRR proxy

The package does not compute a heat-release rate — that requires a gas-phase combustion
model the solid solver deliberately omits (Technical Reference §1). The standard *proxy* is
the mass-loss rate scaled by an assumed effective heat of combustion `ΔH_c`:

```julia
ΔH_c      = 25e6                                 # J/kg, effective heat of combustion (assumed)
hrr_proxy = solution.extras.MLR_gauge .* ΔH_c    # W/m²
```

State this assumption plainly in any comparison; the proxy is only as good as `ΔH_c`.

---

## 9.6 Plotting with the `Plots.jl` recipes

The `PyrolysisPlotsExt` extension attaches a `Plots.jl` recipe for `PyrolysisSolution`. It
loads automatically when both `Pyrolysis` and `Plots` (or `RecipesBase`) are present in your
environment — no explicit `using PyrolysisPlotsExt`. The recipe interprets a `kind` keyword
to route to one of four plot types; the default is the surface-temperature trace.

```julia
using Pyrolysis, Plots

solution = solve(problem)

plot(solution)                      # default: surface temperature vs time
plot(solution; kind = :surface_T)   # surface temperature vs time (explicit)
plot(solution; kind = :mlr)         # mass-loss rate vs time
plot(solution; kind = :residual)    # energy-conservation residual vs time
plot(solution; kind = :heatmap)     # temperature field T(z, t) as a heatmap
```

| `kind` | x | y | source field |
|---|---|---|---|
| `:surface_T` (default) | `t` [s] | surface T [K] | `sol.t`, `sol.surface_T` |
| `:mlr` | `t` [s] | MLR [kg·m⁻²·s⁻¹] | `sol.t`, `sol.mass_loss_rate` |
| `:residual` | `t` [s] | energy residual [W·m⁻²] | `sol.t`, `sol.energy_residual` |
| `:heatmap` | `t` [s] | cell index (depth) | `sol.t`, `axes(sol.T,1)`, `sol.T` |

Any unrecognized `kind` raises an `ArgumentError` listing the valid options.

Notes that follow from the recipe implementation:

- The recipes consume only the fixed core schema (`t`, `T`, `surface_T`, `mass_loss_rate`,
  `energy_residual`), so they work regardless of which `extras` you enabled.
- `:mlr` plots `mass_loss_rate` (per current area). For the cone-gauge curve, plot
  `extras.MLR_gauge` yourself (see §9.5).
- `:residual` plots `energy_residual`, which is zeros unless a conservation tracker filled
  it; a flat line at zero simply means the tracker was off.
- `:heatmap` uses **cell index** on the y axis, not physical depth. For a depth-resolved
  heatmap with metric units, build it yourself from `get_cell_positions` (Eulerian) or
  `solution.z` (ALE).

Because these are standard recipes, all the usual `Plots.jl` attributes pass through:

```julia
plot(solution; kind = :mlr, label = "25 kW/m²", lw = 2, color = :firebrick)
plot!(solution2; kind = :mlr, label = "50 kW/m²")   # overlay a second run
savefig("mlr_comparison.png")
```

### Plotting depth profiles by hand

There is no built-in recipe for "temperature vs depth at a fixed time," but it is a couple of
lines. In **Eulerian** mode use the fixed cell centers:

```julia
using Pyrolysis, Plots
using Pyrolysis.Internal: get_cell_positions

z = get_cell_positions(solution)            # fixed cell centers [m]
plt = plot(xguide = "depth z [mm]", yguide = "T [K]", title = "Temperature profiles")
for t_target in (10.0, 30.0, 60.0)
    i = argmin(abs.(solution.t .- t_target))
    plot!(plt, z .* 1e3, solution.T[:, i]; label = "t = $(solution.t[i]) s")
end
plt
```

To plot a species field instead, swap `solution.T[:, i]` for `solution.ξ[k, :, i]` (component
`k`).

---

## 9.7 Moving meshes: post-processing ALE solutions

When you solve with `use_ale=true`, node positions become part of the state and the mesh
recedes/expands in time (Technical Reference §11). This changes how you read profiles:

- **`solution.z`** is now `(n_nodes, n_steps)` — node positions in metres at each saved
  time. `n_nodes = n_cells + 1`.
- **`solution.extras.thickness_history`** is the sample thickness `z[end] − z[1]` at each
  time (it shrinks as the surface recedes). In Eulerian mode this is a constant equal to the
  domain length.
- If surface-cell **depletion** is enabled (`handle_depletion=true`, which requires
  `use_ale=true`), the cell count drops over the run. The output arrays are pre-allocated at
  the *maximum* cell count and rows for cells that no longer exist are filled with `NaN`.
  Always filter `NaN` before plotting or reducing a moving-mesh profile.
- **`solution.extras.cross_section_area_history`** carries `A(t)` if you used a
  lateral-shrinkage law (WithChi mode); otherwise it stays at `A_0`.

A depth profile at a fixed ALE time uses the cell midpoints reconstructed from the node
positions, and must skip the `NaN`-padded entries:

```julia
using Pyrolysis, Plots

i = argmin(abs.(solution.t .- 30.0))         # nearest saved time to 30 s
z_nodes = solution.z[:, i]                    # node positions at this time [m]
T_cells = solution.T[:, i]

valid = .!isnan.(T_cells)                     # drop depleted/padded cells
z_mid = 0.5 .* (z_nodes[1:end-1] .+ z_nodes[2:end])   # cell-center positions
plot(z_mid[valid] .* 1e3, T_cells[valid];
     xguide = "depth z [mm]", yguide = "T [K]",
     title = "ALE temperature profile at t = $(round(solution.t[i], digits=1)) s")
```

Surface recession over the whole run:

```julia
plot(solution.t, solution.extras.thickness_history .* 1e3;
     xguide = "time [s]", yguide = "thickness [mm]", legend = false,
     title = "Sample thinning")
```

To visualize the moving mesh itself, plot each node's trajectory (every column of
`solution.z` is the node set at one time):

```julia
plt = plot(xguide = "time [s]", yguide = "node position z [m]", legend = false)
for j in axes(solution.z, 1)
    plot!(plt, solution.t, solution.z[j, :])   # NaN gaps render as breaks after depletion
end
plt
```

If the sample fully depletes, the run terminates early (the ALE thickness-termination
callback fires when thickness drops below `min_thickness`), so `solution.retcode` will be
`:Terminated` and `solution.t[end]` will be the burn-through time, not `tspan[2]`.

---

## 9.8 Reading the conservation diagnostics

Pass `diagnostics=true` to `solve` to log snapshots of the conservation ledger. The solver
records a `DiagnosticsRow` every `diagnostics_stride` accepted steps (default 10) into
`solution.extras.diagnostics_log` (Technical Reference §16).

```julia
solution = solve(problem; diagnostics = true, diagnostics_stride = 10)

log = solution.extras.diagnostics_log    # Vector{DiagnosticsRow}; empty if diagnostics=false
```

Each row is a typed `NamedTuple` with these fields:

| Field | Units | Meaning |
|---|---|---|
| `t` | s | Simulation time of the snapshot |
| `E_matrix` | J | Matrix-only sensible energy `Σ ρc_p^eff T V` |
| `E_total` | J | Total sensible energy `Σ ρc_p^total T V` (matrix + gas) |
| `discrete_residual` | J | Discrete energy-balance residual |
| `discrete_relative` | – | `discrete_residual` relative to `|E_matrix(0)|` |
| `physical_residual` | J | Physical energy-balance residual |
| `physical_relative` | – | `physical_residual` relative to `|E_total(0)|` |
| `gas_storage_gap` | J | `ΔE_gas = E_total − E_matrix` (sensible heat in gas, not evolved) |
| `mass_max_relative_error` | – | Max relative mass-conservation error over components |

**Important caveat about snapshot mode.** With `diagnostics=true`, `solve` uses the
*snapshot* recorder, which evaluates only the instantaneous energies. In a snapshot row, the
cumulative-ledger fields — `discrete_residual`, `discrete_relative`, `physical_residual`,
`physical_relative`, `mass_max_relative_error` — are filled with `NaN`, because the
cumulative source integrals are not accumulated in this path. The fields you can trust from
`solve(...; diagnostics=true)` are `t`, `E_matrix`, `E_total`, and `gas_storage_gap`. The
full cumulative residuals require wiring `record_substep!` into the residual layer, which is
beyond the standard `solve` interface.

Reading the energies and the gas-storage gap (the quantity that justifies the matrix-only
heat-capacity convention, Technical Reference §16):

```julia
log      = solution.extras.diagnostics_log
ts       = [row.t for row in log]
E_matrix = [row.E_matrix for row in log]
E_total  = [row.E_total for row in log]
gas_gap  = [row.gas_storage_gap for row in log]

# The gas-storage gap should be a small fraction of the matrix energy
# (< ~0.2% at 1 atm; the physical-form residual should stay < 1%).
rel_gap = gas_gap ./ E_matrix
@show maximum(abs, rel_gap)
```

If you want a true running energy balance with non-`NaN` residuals, the lower-level
`Diagnostics` API (`EnergyConservationTracker`, `update_energy_tracker!`,
`check_energy_conservation`, `energy_conservation_report`) is available via
`using Pyrolysis.Internal`; that path is for advanced auditing and is described in the
Technical Reference §16.

A reminder tied to the energy convention (Technical Reference §7, §16): the default
`energy_form = :standard` carries gas sensible enthalpy through the advection source
`S_conv` and **excludes** a separate gas-generation sink. The surface transpiration BC is
mutually exclusive with that convention — `solve` rejects `use_transpiration_bc=true` to
avoid double-counting — so under default settings
`extras.surface_heat_flux_transpiration` is zero. If you switch to
`energy_form = :with_generation_sink`, the conservative-form sink is added back; the
diagnostics ledger reflects whichever form you chose.

---

## 9.9 Post-processing for cone / TGA comparison

Putting the pieces together, here is a typical reduction of a cone-calorimeter run to the
curves you would overlay on experimental data. The setup (material, BCs, mesh) is shown in
full in the worked-examples chapter; here we focus on the *output* side.

```julia
using Pyrolysis, Plots

solution = solve(problem;
                 abstol = 1e-8, reltol = 1e-6,
                 use_ale = true, handle_depletion = true,
                 min_thickness = 0.01e-3)

@assert solution.retcode in (:Success, :Terminated)

times      = solution.t
T_surface  = solution.surface_T                 # K, exposed-face temperature
T_back     = solution.extras.back_temperature   # K, substrate-face temperature
mlr_gauge  = solution.extras.MLR_gauge          # kg/(m²·s), cone-gauge normalized
mass_m2    = solution.extras.total_mass         # kg/m², remaining mass per initial area

# Mass-loss fraction relative to the initial mass:
mass_frac  = mass_m2 ./ mass_m2[1]

# Surface and back temperature, then MLR, stacked:
p1 = plot(times, T_surface; label = "T_surface",
          xguide = "time [s]", yguide = "T [K]")
plot!(p1, times, T_back; label = "T_back")
p2 = plot(times, mlr_gauge; legend = false,
          xguide = "time [s]", yguide = "MLR [kg/(m²·s)]")
plot(p1, p2; layout = (2, 1), size = (700, 600))
```

For a **TGA-style** idealization (a thin, near-isothermal sample), the natural output is
remaining mass fraction versus temperature. Because the package solves a spatially resolved
slab, take the column total mass and pair it with a representative temperature (the surface,
or a column average):

```julia
mass_frac = solution.extras.total_mass ./ solution.extras.total_mass[1]
T_repr    = solution.surface_T                  # or vec(sum(solution.T, dims=1)) ./ size(solution.T, 1)
plot(T_repr, mass_frac; xguide = "T [K]", yguide = "m/m₀", legend = false,
     title = "Mass-loss curve")
```

For overlaying simulated and measured curves on a common time base, interpolate onto the
experiment's sample times rather than the solver's saved times. The cleanest approach is to
make the solver save exactly where you want to compare — set `output_times` on the problem
(or `saveat` on `solve`) to the experiment's times — which avoids any reliance on the raw
interpolant:

```julia
problem = PyrolysisProblem(; mesh, material, bc_set,
                           tspan = (0.0, t_final),
                           output_times = exp_times)
solution = solve(problem)
# Now solution.extras.MLR_gauge is aligned 1:1 with exp_times.
```

---

## 9.10 Quick reference

```julia
# Core fields (always present; index i selects a saved time):
solution.t                      # times [s]
solution.T[:, i]                # cell temperatures at time i [K]   (bottom→top)
solution.ξ[k, :, i]             # component k concentration profile [kg/m³]
solution.z[:, i]                # node positions [m] (ALE; 0×0 in Eulerian)
solution.surface_T[i]           # exposed-face surface temperature [K]
solution.mass_loss_rate[i]      # MLR per current area (top + permeable-bottom) [kg/(m²·s)]
solution.retcode                # :Success / :Terminated / ...
solution.solve_time             # wall-clock [s]

# extras (always-present keys; some default to zeros):
solution.extras.MLR_gauge                 # MLR per initial area (cone gauge)
solution.extras.back_temperature          # back-face T [K]
solution.extras.total_mass                # domain mass per A_0 [kg/m²]
solution.extras.thickness_history         # sample thickness [m]
solution.extras.gas_generation_rate       # volumetric gas production [kg/(m²·s)]
solution.extras.surface_heat_flux_net     # net surface flux [W/m²]
solution.extras.diagnostics_log           # Vector{DiagnosticsRow} (needs diagnostics=true)

# Accessors (Solver/Discretization submodule exports — qualify or `using Pyrolysis.Internal`):
get_profile(solution, :T, i)              # field at one time index
get_cell_positions(solution)              # fixed cell centers (Eulerian) [m]
interpolate_solution(solution, t)         # raw state vector at continuous t  (== solution(t))
get_surface_temperature(mesh)             # final-state top face T [K]
get_back_face_temperature(mesh)           # final-state bottom face T [K]

# Plotting recipes (need `using Plots`):
plot(solution)                  # :surface_T (default)
plot(solution; kind = :mlr)     # mass-loss rate
plot(solution; kind = :residual)# energy residual
plot(solution; kind = :heatmap) # T(z,t) by cell index
```

See the Technical Reference for the physics behind these outputs: energy assembly and the
`S_conv` convention (§7), the boundary-temperature Newton solve that produces `surface_T`
(§12), ALE node motion behind `z` and `thickness_history` (§10, §11), and the dual energy
ledgers behind the diagnostics rows (§16).


---

# 10. Adaptive Meshing Options

Pyrolysis.jl gives you several ways to let the mesh follow the physics instead of
fixing it once at problem setup. The most important of these — and the only one
wired directly into the production `solve` path — is **surface-cell depletion
handling**, which merges the receding surface cell into its neighbor as the solid
gasifies and the Arbitrary Lagrangian–Eulerian (ALE) mesh shrinks. Beyond that,
the `Pyrolysis.Adaptivity` submodule exposes standalone tools for **r-refinement**
(moving nodes to concentrate them where the solution is steep), **error
indicators** and **monitor functions** that drive those moves, and **mesh-quality
metrics** to keep the mesh well-conditioned. A fully topology-changing h-AMR
cluster (cell splitting and sibling coarsening) also ships, but it lives in
`Pyrolysis.Experimental` and is not connected to `solve`. This chapter shows how
to use each, what every option means and defaults to, and where the sharp edges
are. For the underlying theory (equidistribution, energy-conserving merges, the
geometric conservation law) see Technical Reference §14 (AMR), §11 (ALE), and
§13 (finite-volume discretization).

> **Coordinate reminder.** `z = 0` is the bottom / substrate (boundary tag
> `:bottom`, id 2); `z = L` is the top / exposed surface (tag `:top`, id 1). Heat
> enters from the top and the material recedes downward, so the "surface cell" is
> the last (highest-`z`) active cell. See Chapter 6 (Meshes and Geometry).

---

## 10.1 When to use adaptivity

Adaptivity costs assembly time and, in some configurations, forces a less
efficient Jacobian. Reach for it only when the physics warrants it:

| Situation | Symptom | Recommended tool |
|---|---|---|
| Surface recession (charring/ablating solid in ALE) | Surface cell thins toward zero; timestep collapses; `ρc_p → 0` trouble | `handle_depletion = true` (§10.2) |
| Sample fully consumed | Total thickness approaches zero | thickness termination (automatic in ALE; §10.3) |
| Steep, *moving* thermal or reaction front on a fixed cell count | Front under-resolved between cells, but you don't want more cells | r-refinement via `relocate_nodes!` (§10.5) |
| You suspect mesh distortion after ALE motion | Large cell-size jumps, near-zero or inverted volumes | quality metrics + `enforce_quality!` (§10.6) |
| You want true local h-refinement | — | Experimental h-AMR (§10.7), with caveats |

A good default workflow is: start with a **uniform or surface-stretched fixed
mesh** (Chapter 6) sized from a convergence study (Technical Reference §18.5),
and turn on depletion handling only once you run an ALE case where the surface
actually recedes. Node relocation and h-AMR are advanced, opt-in tools; most
cone-calorimeter and TGA studies never need them.

---

## 10.2 Surface-cell depletion handling (the production path)

When a charring solid pyrolyzes under ALE mesh motion (`use_ale = true`), the
exposed surface cell loses mass and its node moves downward. Eventually that cell
becomes so thin that its discrete heat capacity collapses and the integrator is
forced into vanishingly small steps. **Depletion handling** detects such thin
cells and *merges* them into a neighbor, removing the singularity while
conserving mass and energy.

### 10.2.1 Turning it on through `solve`

This is the recommended way to use depletion handling — you do not build the
callback yourself. Pass `handle_depletion = true` (which **requires**
`use_ale = true`) to `solve`:

```julia
using Pyrolysis

sol = solve(problem;
    use_ale             = true,    # node positions become state; mesh can shrink
    handle_depletion    = true,    # merge thin surface cells instead of stalling
    depletion_threshold = 0.05,    # merge a cell once it drops below 5% of its
                                   #   initial thickness (default; must be in (0,1))
    min_cells           = 2,       # never merge below this many cells; terminate instead
    min_thickness       = 1e-6,    # [m] terminate when total thickness < this (default)
    radiation_model     = SURFACE_ABSORPTION,
)
```

Defaults and meaning of each depletion-related keyword (verified against
`solve(...)` in `src/Solver/solve.jl`):

| Keyword | Default | Meaning |
|---|---|---|
| `use_ale` | `false` | Enable ALE mesh motion (node positions enter the state vector). Prerequisite for depletion. |
| `handle_depletion` | `false` | Detect thin cells and merge them mid-solve. |
| `depletion_threshold` | `0.05` | Cell is "thin" when current thickness `< threshold ×` its own initial thickness (per-cell reference). Must satisfy `0 < threshold < 1`. |
| `min_cells` | `2` | Once the active cell count reaches this, the run terminates rather than merging further. Must be `≥ 1`. |
| `min_thickness` | `1e-6` | [m] Total-sample thickness below which the integrator terminates (separate from per-cell merging). Must be `> 0`. |

Each callback invocation performs **at most one merge** so the integrator can
re-adapt its step between topology changes; thin cells are processed
surface-first (highest `z` first). The merge is **energy-conserving** — see
§10.2.3.

### 10.2.2 What `solve` wires up under the hood

When `use_ale = true && handle_depletion = true`, `solve` assembles a
`CallbackSet` containing (among others):

1. a **step-update callback** that mirrors the integrator's accepted state back
   into the mutable mesh (so subsequent callbacks see current geometry);
2. a **thickness-termination callback** (§10.3);
3. a **depletion callback** that finds thin cells, merges one, resizes the
   workspace caches, and calls `resize!(integrator, …)` to shrink the ODE state.

You can build these manually if you are driving the integrator yourself, but for
ordinary use the `solve` keywords above are all you need. The public constructors
(all exported by `using Pyrolysis`) are:

```julia
# Build depletion callback state (tracks layout + per-cell initial thicknesses + merge log)
cb_state = DepletionCallbackState(layout, mesh)         # convenience: thicknesses from mesh
# or  DepletionCallbackState(layout, initial_cell_thickness::Vector{Float64})  # one entry per cell

# Build the depletion callback itself
dep_cb = create_depletion_callback(ws, mesh, material, cb_state;
    threshold = 0.05,   # same meaning as depletion_threshold above
    min_cells = 2)

# Termination when the whole sample gets too thin
term_cb = create_thickness_termination_callback(cb_state, 1e-6)  # min_thickness [m]
```

`create_depletion_callback` accepts either a bare `StateLayout` or a
`DepletionCallbackState` as its fourth argument; passing the shared
`DepletionCallbackState` lets the termination callback see the resized layout
after each merge. The `material` argument is what makes the temperature merge
energy-conserving (omit it and the merge falls back to a volume-weighted
average — see §10.2.3).

### 10.2.3 The merge is mass- and energy-conserving

Concentrations are merged volume-weighted (exactly mass-conserving):

```
ξ_merged[j] = (V_remove · ξ_remove[j] + V_keep · ξ_keep[j]) / (V_remove + V_keep)
```

Temperature is merged by **conserving the condensed-phase sensible enthalpy**,
not by naive volume averaging. With `m_j = V · ξ_j` the per-component condensed
mass of each cell, the merged temperature solves

```
Σ_{j ∉ gas} m_j ∫_{T_ref}^{T_merged} c_p,j(T′) dT′  =  H_remove + H_keep
```

(a bracketed Newton solve — the root always lies between the two cell
temperatures). The integral form is exact for any temperature-dependent
`c_p,j(T)`, and gas heat capacity is excluded because Formulation A stores
thermal energy in the condensed matrix only. For constant `c_p` this reduces to
the familiar `ρc_p`-weighted mean.

This matters precisely in the depletion regime: merging a nearly-depleted
(low-`ρc_p`) surface cell with a dense interior cell can introduce a 10–20 %
temperature error if you average naively (e.g. ~120 K error when merging a 1 m³
solid at 800 K with a 1 m³ gas-filled cell at 400 K). The dry-solid mass used as
the pyrolysis-progress denominator (`m_dry_solid_initial`) is summed across the
merged cells, so progress accounting stays consistent. See Technical Reference
§14.3 for the derivation and §10 for the χ progress variable.

The merge prefers the **interior** neighbor (toward the substrate) for thermal
stability; this is `find_merge_neighbor` from `Pyrolysis.Adaptivity`, used
internally by the callback.

### 10.2.4 Inspecting what happened

After a solve that used depletion handling, query the merge history with the
exported helper (note the slightly fragile internals — it returns sentinel
values rather than throwing if it cannot reach the callback state):

```julia
stats = get_depletion_callback_stats(dep_cb)
# NamedTuple:
#   stats.n_merges     :: Int               number of cell merges performed
#   stats.merge_times  :: Vector{Float64}   simulation times [s] of each merge
#   stats.merged_cells :: Vector{Int}       original indices of merged cells
```

In ALE mode the solution arrays carry their own evidence of merges: the cell
count shrinks over time and removed-cell entries are NaN-padded in
`sol.T`, `sol.ξ`, and `sol.z` (see Chapter 9). `sol.extras.thickness_history`
records total sample thickness per output time.

### 10.2.5 Critical caveat: depletion is incompatible with a pluggable Jacobian

Depletion resizes the ODE state vector mid-solve via `resize!(integrator, …)`.
The structured/dense Jacobian backends carry an **immutable cache built once at
problem setup** that cannot survive a resize. Consequently:

- Passing **both** `jacobian = …` and `handle_depletion = true` to `solve`
  raises an error at validation time.
- With `handle_depletion = true` and the default `jacobian = nothing`, `solve`
  supplies **no** analytic Jacobian: the integrator falls back to its own
  **colored finite-difference** Jacobian. Because that dense FD Jacobian is
  incompatible with the sparse default linear solver, `solve` also switches the
  linear solver to a dense `LUFactorization`.

In short: **you cannot combine depletion handling with `Structured(...)` /
`DenseAD_ForwardDiff()` backends.** If you need the structured backend (e.g. for
a large ALE problem), run without depletion and rely on thickness termination
instead, or accept the FD fallback. See §10.8 and Technical Reference §14.7 and
§15.13.

> **MKL note.** The dense `LUFactorization` that `solve` selects under depletion
> is Julia's native LU and survives `resize!`. Do not pass a sparse
> factorization (`SparspakFactorization()`, KLU) together with depletion — there
> is no sparse Jacobian prototype in this mode, and `solve` will warn and swap
> it for dense LU. The MKL `ipiv` resize bug can only bite if you construct your
> own `integrator` without an explicit `linsolve`; the solver warns if it
> detects that combination on a large problem. See Chapter 11 and
> Technical Reference §15.12.

---

## 10.3 Thickness-termination callback

Independently of merging, a shrinking sample can reach (near-)zero thickness,
which makes the volumetric heat capacity singular. In **ALE mode `solve`
automatically installs a thickness-termination callback** that stops the
integration cleanly once the total sample thickness drops below
`min_thickness` (default `1e-6` m):

```julia
sol = solve(problem; use_ale = true, min_thickness = 1e-5)  # terminate at 10 µm
@show sol.retcode    # :Terminated  (vs :Success when the full tspan is reached)
```

The callback measures `thickness = z[end] - z[1]` from the node-position block of
the state vector, so it is meaningful only in ALE mode (`use_ale = true`); in
Eulerian mode the mesh does not move and the callback is not installed. To build
it directly:

```julia
term_cb = create_thickness_termination_callback(layout, Float64(min_thickness))
```

Pass a `DepletionCallbackState` instead of a bare `layout` when you also use
depletion handling, so the callback tracks the resized layout after each merge.

---

## 10.4 Error indicators (driving refinement decisions)

Error indicators estimate where the discretization error is largest. They are the
inputs to refinement/relocation decisions: high indicator value ⇒ resolve more
there. All indicator types live in `Pyrolysis.Adaptivity` and are re-exported, so
`using Pyrolysis` makes them available directly. Each implements
`compute_indicator(ind, mesh)` (returning a per-cell `Vector{Float64}`, with `0`
for inactive cells).

| Indicator | Constructor | Formula (η for cell *i*) | Use |
|---|---|---|---|
| Gradient | `GradientIndicator(:T)` or `GradientIndicator(:T, weight)` | `η = w · hᵢ · \|∇φ\|ᵢ` | Steep fronts (thermal/reaction). Default `weight = 1.0`. |
| Curvature | `CurvatureIndicator(:T)` or `(:T, weight)` | `η = w · hᵢ² · \|∇²φ\|ᵢ` | Sharper features than gradient; zero at boundaries (needs 3 points). |
| Reaction rate | `ReactionRateIndicator()` or `(weight)` | `η = w · hᵢ · maxⱼ \|rⱼ\|ᵢ / max(ξⱼ, ε)` | Active reaction zones. Requires `rates = r` (throws `ArgumentError` without it); standalone-only — it cannot run inside the h-AMR controller cycle. |
| Interface proximity | `InterfaceProximityIndicator(pos, h_target)` or `(pos, h_target, weight)` | `η = w / (1 + dᵢ/h_target)` | Force resolution near a tracked interface position, regardless of solution. Purely geometric. |
| Composite | `CompositeIndicator(inds)` or `(inds, combination)` | `:max`, `:l2`, or `:sum` of components, each first normalized by its own maximum over active cells | Blend several indicators. Default `combination = :max`. |

`φ` is the tracked field: `:T` for temperature, or a component index for a
concentration. `hᵢ` is the cell size and `dᵢ` the distance from the cell center
to the interface. The reaction-rate indicator's `ε = 1e-10` guards the division.

```julia
using Pyrolysis

# A single indicator
η_T = compute_indicator(GradientIndicator(:T), mesh)        # Vector{Float64}, length n_cells

# Blend a temperature gradient with temperature curvature, take the cell-wise max
ind = CompositeIndicator(
    AbstractErrorIndicator[GradientIndicator(:T, 1.0), CurvatureIndicator(:T, 0.5)],
    :max,
)
η = compute_indicator(ind, mesh)

# Reaction-rate indicator: standalone only, with explicit rates (sized like ξ)
η_r = compute_indicator(ReactionRateIndicator(), mesh, T, ξ; rates = r)
```

> **Note.** `compute_indicator(ind, mesh)` reads `T` and `ξ` from the mesh's cell
> states, so make sure the mesh is in sync with the current solution (after a
> solve, the step-update callback keeps `mesh` consistent; if you are
> post-processing a `PyrolysisSolution`, rebuild the field arrays yourself and
> use the four-argument form `compute_indicator(ind, mesh, T, ξ)`). The
> reaction-rate indicator requires a rate matrix that is not available from the
> mesh alone and **throws** in the convenience form — pass `rates = r` (sized
> like `ξ` on the current mesh) via the four-argument form.

These indicators are **diagnostic building blocks** — computing one does not
change the mesh. To act on an indicator you either feed it into the experimental
h-AMR machinery (§10.7) or, more commonly, use a *monitor function* for
r-refinement (§10.5).

---

## 10.5 r-refinement: node relocation by equidistribution

r-refinement keeps the **same number of cells** but slides interior nodes so they
cluster where the solution is hard to resolve. It is driven by a **monitor
function** `ω(z)`; nodes are placed so that each cell carries roughly equal
"monitor integral" (the equidistribution principle, Technical Reference §14.2).
Boundary nodes never move.

### 10.5.1 Monitor functions

| Monitor | Constructor | Formula | Notes |
|---|---|---|---|
| Gradient | `GradientMonitor()` or `GradientMonitor(; α=1.0, variable=:T)` | `ω = √(1 + α·\|∇φ\|²)` | Clusters nodes at steep gradients. Defaults `α = 1.0`, `variable = :T`. |
| Interface | `InterfaceMonitor(z)` or `InterfaceMonitor(z; σ=0.001, ω_base=1.0, ω_peak=10.0)` | `ω = ω_base + ω_peak·exp(−(z−z_interface)²/(2σ²))` | Gaussian peak at interface position `z`. |
| Composite | `CompositeMonitor(monitors, combination)` | `:max`, `:sum`, or `:product` | Combine monitors. |

### 10.5.2 Relocating nodes

```julia
using Pyrolysis

n_iter = relocate_nodes!(mesh;
    monitor    = GradientMonitor(; α = 1.0, variable = :T),
    max_iter   = 10,     # max equidistribution sweeps
    relaxation = 0.5,    # step fraction toward target: z_new = z + β·(z_target − z), β ∈ (0,1]
    tol        = 1e-3,   # converged when max relative node movement < tol
    h_min      = 1e-8,   # [m] minimum cell size enforced during moves
)
```

Defaults (verified in `src/Adaptivity/r_refinement.jl`): `monitor =
GradientMonitor()`, `max_iter = 10`, `relaxation = 0.5`, `tol = 1e-3`,
`h_min = 1e-8`. The function returns the number of iterations performed and
mutates `mesh` in place (updating node positions and re-deriving geometry).

What the algorithm does each sweep: evaluate the monitor, compute target node
positions by equidistributing the monitor integral, then take a relaxed step
toward each target while enforcing `z_{n+1} − z_n ≥ h_min`. Because the
minimum-spacing constraint can prevent nodes from reaching their targets, the
returned iteration count may equal `max_iter` without full convergence.

> **Caveats.** Monitor functions are problem-dependent and have arbitrary units;
> you are responsible for scaling so the equidistribution solve is
> well-conditioned (huge `ω` disparities cause ill-conditioning). Default
> tolerances (`tol = 1e-3`, `relaxation = 0.5`) are untested on realistic
> steep-front pyrolysis problems — treat them as starting points. Boundary nodes
> are pinned. r-refinement is a **standalone mesh operation**: it is not invoked
> automatically by `solve`. Use it between solve segments (relocate, then warm-
> start the next segment) or in a custom driver.

---

## 10.6 Mesh-quality metrics and enforcement

After any mesh motion (ALE or r-refinement) you may want to confirm the mesh is
still well-conditioned, and optionally repair it. The quality types and functions
are in `Pyrolysis.Adaptivity` (re-exported).

### 10.6.1 Metrics

| Metric | Constructor (default field) | Constraint per adjacent pair / cell |
|---|---|---|
| Aspect ratio | `AspectRatioMetric()` → `max_ratio = 2.0` | `hᵢ/hⱼ ≤ max_ratio` |
| Smoothness | `SmoothnessMetric()` → `max_size_jump = 0.5` | `\|hᵢ − hⱼ\| / max(hᵢ, hⱼ) ≤ max_size_jump` |
| Jacobian (positivity) | `JacobianMetric()` → `min_jacobian = 0.0` | `Vᵢ > min_jacobian` (signed width — inverted or collapsed cells fail) |
| Minimum cell size | `MinCellSizeMetric()` → `h_min = 1e-8` | `hᵢ ≥ h_min` |
| Composite | `CompositeQualityMetric(metrics)` | All component metrics must pass |

(Defaults verified in `src/Adaptivity/mesh_quality.jl`.)

### 10.6.2 Checking and summarizing

```julia
using Pyrolysis

metric = AspectRatioMetric(2.0)

check_quality(mesh, metric)        # Bool: does the whole mesh pass?
quality_violations(mesh, metric)   # Vector{Int}: cell indices that violate it

compute_aspect_ratios(mesh)        # Vector{Float64}: per-cell aspect ratio
compute_smoothness(mesh)           # Vector{Float64}: per-cell smoothness

summary = mesh_quality_summary(mesh)
# NamedTuple with:
#   n_active_cells, min_cell_size, max_cell_size,
#   max_aspect_ratio, max_smoothness, all_positive_volume
```

### 10.6.3 Enforcing quality

```julia
remaining = enforce_quality!(mesh, metric; max_iter = 10)
# Returns the number of remaining violations after Laplacian smoothing.
```

`enforce_quality!` applies **Laplacian (neighbor-averaging) smoothing** to the
interior nodes of violating cells, repeating until either no violations remain or
`max_iter` is hit, and returns the count of remaining violations (zero means the
mesh now passes). It is a local, incremental fix and is **not guaranteed to clear
clustered or severe violations** — check the return value. Boundary nodes are
not moved.

A typical post-motion guard:

```julia
relocate_nodes!(mesh; monitor = GradientMonitor(; variable = :T))
metric = CompositeQualityMetric(AbstractQualityMetric[
    AspectRatioMetric(2.0), MinCellSizeMetric(1e-7),
])
if !check_quality(mesh, metric)
    left = enforce_quality!(mesh, metric; max_iter = 20)
    left == 0 || @warn "mesh still has $left quality violation(s) after smoothing"
end
```

---

## 10.7 Experimental h-AMR (cell splitting and coarsening)

True topology-changing h-adaptive refinement — splitting cells into children and
merging sibling cells back — ships with the package but is **quarantined in
`Pyrolysis.Experimental` and is *not* connected to `solve`.** Unlike the
`Adaptivity` exports, `Experimental` is **not re-exported at the top level**, so
you must reach into it explicitly:

```julia
using Pyrolysis

# Single-cell split: returns (child_a_id, child_b_id)
ca, cb = Pyrolysis.Experimental.refine_cell!(mesh, cell_id)

# Threshold-based marking + batch refinement with 2:1 level balance
marked = Pyrolysis.Experimental.mark_cells_for_refinement(
    mesh, η, θ_refine, h_min, max_level)
marked = Pyrolysis.Experimental.enforce_2to1_balance!(mesh, marked, max_level)
Pyrolysis.Experimental.perform_h_refinement!(mesh, marked)

# Sibling coarsening (the inverse, for cells created by a previous refinement)
Pyrolysis.Experimental.coarsen_cells!(mesh, cell_a_id, cell_b_id; material = material)

# Orchestrator
ctrl = Pyrolysis.Experimental.AMRController(;
    error_indicator = GradientIndicator(:T),
    quality_metric  = AspectRatioMetric(),
    θ_refine    = 0.7,    # refine if η > θ_refine · max(η)
    θ_coarsen   = 0.1,    # coarsen if η < θ_coarsen · max(η)
    h_min       = 1e-6,
    h_max       = 1e-2,
    max_level   = 4,
    check_every = 10,
)
result = Pyrolysis.Experimental.adapt_mesh!(mesh, ctrl)   # (n_refined, n_coarsened, quality_ok)
```

The full experimental surface (all under `Pyrolysis.Experimental`) includes
`refine_cell!`, `perform_h_refinement!`, `mark_cells_for_refinement`,
`enforce_2to1_balance!`, `verify_refinement_conservation`, `coarsen_cells!`,
`perform_h_coarsening!`, `mark_cells_for_coarsening`, `can_coarsen`,
`get_sibling`, `verify_coarsening_conservation`, `hr_adapt!` (combined h-then-r),
`AMRController`, `adapt_mesh!`, `should_adapt`, `maybe_adapt!`,
`amr_statistics`, `reset_statistics!`, and `create_default_controller`.

**Caveats — read before using:**

- **Not in the solve loop.** There is no `solve` keyword to enable h-AMR during
  integration. Calling these between solve segments is the only supported pattern,
  and even then the state-vector resize logic that depletion handling uses is not
  shared here.
- **Bookkeeping is leaky.** h-refinement appends children and h-coarsening marks
  cells inactive *without compacting the arrays*, so `length(mesh.cells) >
  n_active_cells` is normal. Any code consuming the mesh afterward must iterate
  the `active_cell_mask`, not raw index ranges. `assert_amr_invariants(mesh)` (a
  core `Adaptivity` export) is called after each operation to catch drift.
- **Conservation is exact, but the temperature split is zeroth-order.** Children
  inherit the parent concentration with equal volume split (mass-exact); the
  temperature split is piecewise-constant — both children inherit the parent
  temperature (the coded gradient correction evaluates to zero because the split
  point coincides with the parent cell center) — exactly energy-conserving but
  adding no sub-cell resolution. Use
  `verify_refinement_conservation` / `verify_coarsening_conservation` to check.
- **Defaults are untuned.** Controller thresholds (`θ_refine = 0.7`,
  `θ_coarsen = 0.1`, `max_level = 4`, `check_every = 10`) have not been validated
  on realistic pyrolysis fronts.

For production surface-recession work, prefer the depletion path (§10.2), which
*is* integrated with the solver and uses the energy-conserving adjacent-cell
merge (`merge_adjacent_cells!`, `merge_surface_cell!`) that lives in core
`Adaptivity`. See Technical Reference §14.4 for the h-AMR algorithms and §14.3
for the production merge.

---

## 10.8 Choosing the Jacobian when you use adaptivity

The single most common configuration mistake with adaptivity is combining
depletion handling with an explicit Jacobian backend. The rule:

| Scenario | Jacobian to use |
|---|---|
| ALE, **no** depletion handling | Default `jacobian = nothing` ⇒ `Structured(strategy = DirectSolve())`; or opt into `ApproxSolve` for large ALE problems (Chapter 11). |
| ALE **with** `handle_depletion = true` | **Leave `jacobian = nothing`.** `solve` supplies no analytic Jacobian and the integrator uses colored FD with a dense `LUFactorization`. Do **not** pass a `JacobianSpec` — it errors. |
| Verifying a backend during an adaptive run | `solve(...; verify_jacobian = (rtol=1e-6, atol=1e-10, sample_times=[…]))` — but only meaningful when a backend is active (i.e. not with depletion). |

Concretely, this is correct:

```julia
# Large ALE problem WITHOUT depletion — opt into the cheaper approximate solve
sol = solve(problem;
    use_ale = true,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
)
```

and this **errors** at validation:

```julia
# WRONG: pluggable Jacobian + depletion are mutually exclusive
solve(problem;
    use_ale = true,
    handle_depletion = true,
    jacobian = JacobianSpec(Structured(strategy = DirectSolve()), problem),  # ← rejected
)
```

If you need both the analytic structured Jacobian *and* surface recession, the
practical recourse is to forgo per-cell merging and rely on the
thickness-termination callback alone (run until the sample is "thin enough,"
then stop), keeping the state-vector length fixed so the structured backend's
cache stays valid. See Chapter 11 for the full Jacobian/linear-solver decision
guide and Technical Reference §15.10–§15.13 for the strategy internals.

---

## 10.9 Summary of options and defaults

| Option | Where | Default | Section |
|---|---|---|---|
| `use_ale` | `solve` kwarg | `false` | 10.2 |
| `handle_depletion` | `solve` kwarg | `false` | 10.2 |
| `depletion_threshold` | `solve` kwarg | `0.05` | 10.2 |
| `min_cells` | `solve` kwarg | `2` | 10.2 |
| `min_thickness` | `solve` kwarg | `1e-6` m | 10.3 |
| `GradientIndicator(:T, w)` | `Adaptivity` | `w = 1.0` | 10.4 |
| `CompositeIndicator(inds, c)` | `Adaptivity` | `c = :max` | 10.4 |
| `GradientMonitor(; α, variable)` | `Adaptivity` | `α=1.0, variable=:T` | 10.5 |
| `InterfaceMonitor(z; σ, ω_base, ω_peak)` | `Adaptivity` | `σ=1e-3, ω_base=1.0, ω_peak=10.0` | 10.5 |
| `relocate_nodes!(; max_iter, relaxation, tol, h_min)` | `Adaptivity` | `10, 0.5, 1e-3, 1e-8` | 10.5 |
| `AspectRatioMetric(max_ratio)` | `Adaptivity` | `2.0` | 10.6 |
| `SmoothnessMetric(max_size_jump)` | `Adaptivity` | `0.5` | 10.6 |
| `MinCellSizeMetric(h_min)` | `Adaptivity` | `1e-8` m | 10.6 |
| `enforce_quality!(; max_iter)` | `Adaptivity` | `10` | 10.6 |
| h-AMR (`refine_cell!`, `AMRController`, …) | `Experimental` (not re-exported) | — experimental | 10.7 |

**Key correctness reminders.** Depletion handling requires `use_ale = true`; it
cannot be combined with a pluggable `jacobian` backend (it falls back to colored
FD with dense LU); the surface merge conserves species mass exactly and condensed-phase sensible
enthalpy via the `∫c_p dT` temperature merge (gas excluded); r-refinement,
quality enforcement, and h-AMR
are standalone tools not driven by `solve`; and the h-AMR cluster is experimental
and reached only through `Pyrolysis.Experimental`. For the physics behind these
operations, see Technical Reference §11 (ALE), §14 (AMR), and §15 (Jacobian).


---

# 11. Solver and Jacobian Configuration; Extensions

This chapter is the reference for tuning *how* Pyrolysis.jl solves a problem,
as opposed to *what* problem it solves (Chapters 3–7) or *how to read the
result* (Chapter 9). It covers choosing a Jacobian backend and linear-solve
strategy by problem size and simulation mode, configuring the integrator's
linear solver (sparse-direct or ILU-preconditioned Krylov), loading the optional
vendor and sensitivity extensions, debugging Jacobian correctness with
`verify_jacobian_at`, and dropping down to the pure-residual layer
(`ProblemDef` / `Workspace` / `residual!`) when you need to build a custom
solver loop. The mathematics behind the structured Jacobian and the implicit-RK
`W` matrix is derived in Technical Reference §15; this chapter is about the
knobs.

Two distinct linear systems appear in a stiff implicit solve, and Pyrolysis.jl
lets you configure each one separately. Keeping them straight is the key to this
chapter:

1. **The Jacobian backend** (`jacobian = ...` keyword) decides *how the matrix*
   `J = ∂(du/dt)/∂u` *is assembled* — by structured operator assembly (the
   default: `solve` auto-builds `Structured(strategy = DirectSolve())` when
   `jacobian = nothing`), by dense automatic differentiation (verification
   only), or by letting the integrator build it via colored finite differences
   (the fallback used when depletion handling disables the pluggable backend).
2. **The integrator's linear solver** (`linear_solver = ...` keyword) decides
   *how the implicit-stage system* `W x = b`, with `W = I/(γ_RK Δt) − J`, *is
   solved* once `J` exists — by sparse-direct factorization (Sparspak, KLU,
   Pardiso) or by a preconditioned Krylov iteration (ILU-GMRES, ILU-BiCGSTAB).

You can mix and match: a pluggable structured Jacobian with a sparse-direct
linear solver, an integrator-built finite-difference Jacobian with KLU, and so
on. The defaults (auto-built structured Jacobian; `SparspakFactorization`) work well for
the typical 1D cone-calorimeter problem and you can ignore this chapter until a
solve is too slow or you need sensitivities.

---

## 11.1 The two layers, and what the defaults do

When you call `solve(problem)` with no Jacobian or linear-solver keywords, the
following happens (see `src/Solver/solve.jl`):

- The integrator defaults to `KenCarp4(autodiff = AutoFiniteDiff(),
  linsolve = SparspakFactorization())`. `KenCarp4` is an SDIRK (singly
  diagonally implicit Runge–Kutta) method well suited to the moderate stiffness
  of Arrhenius pyrolysis. The `AutoFiniteDiff()` setting only matters on the
  rare paths where no pluggable Jacobian is active (depletion handling): there
  it makes the integrator finite-difference the residual itself.
- The Jacobian defaults to the **structured analytical backend**: with
  `jacobian = nothing` (the default), `solve` constructs
  `Structured(strategy = DirectSolve())` internally — for both Eulerian and
  ALE problems — and hands the integrator an exact sparse Jacobian plus its
  sparsity prototype. The single exception is `handle_depletion = true`, which
  disables the pluggable backend (the spec cannot survive a mid-solve state
  resize) and falls back to the integrator's colored finite differences.
- The linear solver defaults to `SparspakFactorization()`, a sparse-direct LU
  that is typically faster than KLU for the block-tridiagonal sparsity pattern
  of a 1D problem.

In short: by default, **Layer 1 is the structured analytical Jacobian**
(`Structured(strategy = DirectSolve())`, built for you) and **Layer 2 is
Sparspak**. The `jacobian` keyword is for *overriding* the strategy
(`ApproxSolve`, `StructuredSolve`, `MatrixFreeSolve`) or supplying a custom
spec — not for switching the structured machinery on.

```julia
using Pyrolysis

# The default solve: auto-built structured Jacobian + Sparspak linear solve.
sol = solve(problem)
```

> Defaults at a glance: `integrator = nothing` → `KenCarp4`,
> `linear_solver = nothing` → `SparspakFactorization()`,
> `jacobian = nothing` → auto-built `Structured(strategy = DirectSolve())`
> (colored-FD fallback only under `handle_depletion = true`),
> `abstol = 1e-8`, `reltol = 1e-6`.

---

## 11.2 Choosing a Jacobian backend

A **Jacobian backend** is supplied as a `JacobianSpec` to the `jacobian`
keyword of `solve`. Construct it from a backend tag plus your problem:

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)
```

There are two production-relevant backend tags, both exported at top level:

| Backend tag | Purpose | Use in production? |
|---|---|---|
| `Structured(strategy = ...)` | Operator-based structured assembly (fast, exact or approximate per strategy) | Yes |
| `DenseAD_ForwardDiff()` | Whole-residual ForwardDiff; dense Jacobian | No — verification reference only |

`Structured` is parameterized by a **linear-solve strategy** (§11.3) that
determines how the assembled structured Jacobian is inverted. The default,
selected when you pass `jacobian = nothing`, is `Structured(strategy =
DirectSolve())`.

`DenseAD_ForwardDiff()` runs ForwardDiff over the entire residual vector and
returns a dense matrix. It is correct by construction but slow; its only role is
as the reference against which the structured backend is verified (see §11.7).
Do not use it for production solves.

### 11.2.1 The `JacobianSpec` constructor

`JacobianSpec` has two forms (`src/Jacobian/interface.jl`):

```julia
# From a PyrolysisProblem (the common form). Pass use_ale = true for ALE.
JacobianSpec(backend, problem; use_ale = false, radiation_model = nothing,
             radiation_intensity = 0.0, energy_form = :standard,
             use_transpiration_bc = false, min_thickness = 1e-6)

# From an already-built (ProblemDef, Workspace) pair (low-level form).
JacobianSpec(backend, def, ws)
```

Important: a `JacobianSpec` builds its sparsity pattern and AD workspaces once,
at construction time, and it bakes in the simulation mode. **If your solve uses
ALE, you must tell the spec by passing `use_ale = true`** — otherwise the
sparsity pattern omits the geometry (z-column) block and the Jacobian will be
wrong:

```julia
# WRONG for an ALE solve — the spec is built for the Eulerian mode.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem)
solve(problem; jacobian = spec, use_ale = true)   # mode mismatch

# RIGHT — the spec's mode matches the solve.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
solve(problem; jacobian = spec, use_ale = true)
```

If a problem also uses radiation, pass the matching `radiation_model` (e.g.
`radiation_model = BEER_LAMBERT`) so the spec includes the Beer–Lambert
non-local sparsity. Similarly, pass `energy_form = :with_generation_sink` if you
solve with that energy form. For most users the simplest correct approach is to
let `solve` build the spec for you — see §11.2.3.

### 11.2.2 When *not* to supply a Jacobian

A pluggable `JacobianSpec` is **incompatible with surface-depletion handling**.
Depletion (`handle_depletion = true`) resizes the state vector mid-solve via
`resize!(integrator, ...)`, but the spec's cache is built once at setup time and
cannot survive a resize. `solve` rejects the combination at validation:

```julia
# ERROR: handle_depletion + a pluggable jacobian is not allowed.
solve(problem; jacobian = spec, handle_depletion = true, use_ale = true)
```

When you need depletion, leave `jacobian = nothing`; the integrator falls back
to its own colored finite-difference Jacobian, which tolerates the resize. This
is the only supported configuration for depletion-tracking solves (see Chapter
10 and Technical Reference §14.7).

### 11.2.3 Letting `solve` pick the backend

For routine work, the cleanest pattern is to pass `jacobian = nothing` (the
default) and let `solve` construct the default structured Jacobian for you, *or* to
opt into a non-default backend explicitly only when you need it. The auto-default when
`jacobian` is left at its default is `Structured(strategy = DirectSolve())`
constructed internally with the correct mode flags — you do not have to thread
`use_ale` / `radiation_model` through `JacobianSpec` yourself. Supply your own
spec only when you want a non-default strategy (`ApproxSolve`, `StructuredSolve`,
`MatrixFreeSolve`) or a custom reference for verification.

---

## 11.3 Linear-solve strategies for the structured Jacobian

The `strategy` argument to `Structured(...)` is a `LinearSolveStrategy`. It
governs how the structured Jacobian — a sparse local block `A` plus a list of
fast-apply non-local couplings (prefix-sum, suffix-sum, rank-1) — is inverted
inside the Newton step. All four strategies are top-level exports.

| Strategy | What it does | Exact? | Best for |
|---|---|---|---|
| `DirectSolve()` | Densify the structured couplings into `A`, factor the full sparse matrix (KLU/UMFPACK) | Exact | Small/typical problems; the verification anchor |
| `StructuredSolve()` | Sparse-LU on `A` + preconditioned Richardson over the couplings (2–5 iters) | Exact | Large problems where dense assembly is costly but exactness is required |
| `ApproxSolve()` | Drop the structured couplings; solve `A⁻¹ b` only (Schur-style approximation) | Approximate | Large ALE problems where speed dominates (~7.8× faster assembly) |
| `MatrixFreeSolve()` | GMRES over the operator action without densifying | Exact (iterative) | Research / very large `N`; benchmark before relying on it |

Construct any of these as a parameterless instance:

```julia
JacobianSpec(Structured(strategy = DirectSolve()),     problem)               # default
JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
JacobianSpec(Structured(strategy = ApproxSolve()),     problem; use_ale = true)
JacobianSpec(Structured(strategy = MatrixFreeSolve()), problem; use_ale = true)
```

### 11.3.1 `DirectSolve` — the default and the anchor

`DirectSolve` materializes `A + Σⱼ densify(couplingⱼ)` and factors it directly.
It is exact and predictable. For Eulerian problems and modest ALE meshes it is
the right choice and the one you get when you do nothing. It is also the backend
you should verify against first when proving a configuration is correct.

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)
```

### 11.3.2 `StructuredSolve` — exact, cheaper assembly at scale

`StructuredSolve` keeps the solve exact but avoids densifying the non-local
ALE/radiation couplings. It factors only the sparse local block `A`, then runs a
preconditioned Richardson iteration that applies each coupling through its O(N)
`apply!` method:

```
F = lu(A); y = F \ b
repeat:  r = b − A·y − Σⱼ Cⱼ·y;  Δy = F \ r;  y ← y + Δy   (2–5 iterations)
```

Convergence is fast (2–5 iterations to machine precision) when the coupling
magnitudes are small relative to `A`'s diagonal — which is the usual ALE regime.
Choose `StructuredSolve` when an ALE mesh is large enough that the dense
assembly of `DirectSolve` dominates, but you still need an exact Newton step.

### 11.3.3 `ApproxSolve` — fastest, for large ALE

`ApproxSolve` drops the structured couplings entirely and solves `A⁻¹ b`. The
Newton step is then a Schur-style approximation: it omits the dense ALE
mesh-velocity coupling block, takes a few extra Newton iterations, but each
linear solve and the assembly itself are far cheaper (the assembly profile drops
to `:ale_reduced`, skipping the cumulative ALE widening, the geometry AD pass,
the Beer–Lambert suffix, and the χ-column rank-1 term). This is the strategy the
Douglas-fir ALE example uses for roughly an 8× speedup over the full Jacobian.

```julia
# Large ALE cone problem — fast approximate Jacobian.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)
```

`ApproxSolve` is the **only** strategy that deliberately omits known-nonzero
entries; it carries the `is_approximate` trait, which switches
`verify_jacobian_at` (§11.7) into masked-comparison mode so the dropped entries
are reported informationally rather than counted as failures.

### 11.3.4 `MatrixFreeSolve` — research / very large N

`MatrixFreeSolve` runs GMRES over the structured-Jacobian operator action
without ever forming `A + C`. It is exact in the limit of GMRES convergence and
useful when even materializing the coupled operator is expensive. It is marked
experimental in the source; benchmark a specific problem before using it in
production, and prefer `DirectSolve`/`StructuredSolve` unless the matrix-free
path measurably wins.

### 11.3.5 A sizing rule of thumb

- Eulerian (fixed mesh), any size, no depletion → `jacobian = nothing`
  (default) is fine; it already gives you
  `Structured(strategy = DirectSolve())`.
- ALE, small/medium mesh → `DirectSolve`.
- ALE, large mesh, exactness required → `StructuredSolve`.
- ALE, large mesh, speed over exactness → `ApproxSolve`.
- Depletion tracking → no pluggable Jacobian (`jacobian = nothing`).

---

## 11.4 Configuring the integrator's linear solver

Independently of the Jacobian backend, you choose how the per-stage implicit
system `W x = b` (`W = I/(γ_RK Δt) − J`, see Technical Reference §15.5) is
solved. This is the `linear_solver` keyword. There are two families:
sparse-direct factorizations and ILU-preconditioned Krylov iterations.

### 11.4.1 Sparse-direct factorizations

The default is `SparspakFactorization()` (about 28% faster than KLU for the
block-tridiagonal 1D pattern). `KLUFactorization()` is a slightly more
general-purpose alternative with the same interface. You can pass either by
pairing it with a *constructor* integrator:

```julia
using LinearSolve: KLUFactorization, SparspakFactorization

# Default (no need to write this explicitly):
sol = solve(problem; linear_solver = SparspakFactorization())

# Switch to KLU:
sol = solve(problem; linear_solver = KLUFactorization())
```

> Pairing rule: `linear_solver` is combined with `integrator` only when
> `integrator` is a *constructor* such as `KenCarp4` (or left at its default).
> If you want a fully constructed algorithm with a non-default linear solver,
> build it yourself and pass it as the `integrator`:
> `solve(problem; integrator = KenCarp4(linsolve = KLUFactorization()))`.
> Passing both a fully constructed `integrator` and a separate `linear_solver`
> raises an error.

A note on the no-Jacobian path: with `handle_depletion = true` (or a dense
backend such as `DenseAD_ForwardDiff`), there is no sparse `jac_prototype`,
so the sparse `Sparspak`/`KLU` factorizations cannot run. `solve` then
defaults the linear solver to a dense `LUFactorization()` — and if you
explicitly request a sparse factorization it warns and rebuilds your
algorithm with `LUFactorization` instead. On the ordinary default path
(auto-built structured Jacobian) the sparse prototype is present and
Sparspak applies.

### 11.4.2 ILU-preconditioned GMRES

For large ALE problems where direct factorization of `W` dominates the cost,
Pyrolysis.jl ships an ILU-preconditioned GMRES solver,
`ILUGMRESFactorization` (exported by the `Solver` submodule; qualify as
`Pyrolysis.ILUGMRESFactorization` or `using Pyrolysis.Solver`). It builds an incomplete-LU
preconditioner from the concrete sparse `W` and iterates GMRES against it,
falling back to diagonal scaling if the ILU factorization fails on a stiff
transient.

```julia
# ILUGMRESFactorization(; τ = 0.1, gmres_rtol = 1e-8, restart = 30)
alg = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1))

sol = solve(problem;
    integrator = alg,
    jacobian   = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale    = true)
```

Tuning parameters and their defaults (the field name `τ` is the ILU drop
tolerance, written `τ_ILU` in the nomenclature):

| Parameter | Default | Meaning |
|---|---|---|
| `τ` | `0.1` | ILU drop tolerance. **Smaller** → denser ILU → better preconditioner, slower setup. |
| `gmres_rtol` | `1e-8` | GMRES relative convergence tolerance. |
| `restart` | `30` | GMRES restart frequency (memory parameter). |

Caveat: inexact linear solves can increase Newton rejections on Arrhenius
transients, which sometimes *raises* total wall time despite cheaper individual
solves. Always benchmark an ILU-GMRES configuration end-to-end against the
default `SparspakFactorization()` before adopting it for a production study.

### 11.4.3 ILU-preconditioned BiCGSTAB

`ILUBiCGSTABFactorization` (exported by the `Solver` submodule; qualify or
`using Pyrolysis.Solver`) is the BiCGSTAB analogue. It uses
the same ILU preconditioner but a fixed-memory BiCGSTAB iteration (no restart
parameter), which can be faster for some matrix structures:

```julia
# ILUBiCGSTABFactorization(; τ = 0.1, rtol = 1e-8)
alg = KenCarp4(linsolve = ILUBiCGSTABFactorization(τ = 0.1, rtol = 1e-8))
sol = solve(problem; integrator = alg,
    jacobian = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale  = true)
```

| Parameter | Default | Meaning |
|---|---|---|
| `τ` | `0.1` | ILU drop tolerance. |
| `rtol` | `1e-8` | BiCGSTAB relative convergence tolerance. |

### 11.4.4 The MKL `ipiv` caveat and depletion

When surface depletion resizes the state vector mid-solve, the integrator must
also resize its linear-solver workspace. The MKL-backed pivot array (`ipiv`)
used by some factorizations does not resize cleanly across a state-vector
length change. Depletion forces `jacobian = nothing` (§11.2.2), so there is no
sparse Jacobian prototype and the sparse factorizations cannot run on this
path: `solve` defaults the linear solver to Julia's native, resize-safe
`LUFactorization()` and rebuilds an explicitly requested
`SparspakFactorization`/`KLUFactorization` to `LUFactorization`, with a
warning. A robust depletion configuration is therefore simply:

```julia
# Depletion-safe: no pluggable Jacobian; solve picks the resize-safe dense LU.
sol = solve(problem;
    use_ale             = true,
    handle_depletion    = true,
    depletion_threshold = 0.05)
```

The MKL pitfall survives only when you pass a fully constructed integrator
whose `linsolve` is unset (LinearSolve's auto-selection may then pick MKL LU);
`solve` warns when it detects that combination on a large depleting state.

---

## 11.5 Tolerances and integrator selection

The ODE error tolerances are independent of both Jacobian and linear-solver
choices:

| Keyword | Default | Meaning |
|---|---|---|
| `abstol` | `1e-8` | Absolute error tolerance. |
| `reltol` | `1e-6` | Relative error tolerance. |
| `maxiters` | `1000000` | Maximum integration steps. |
| `dtmin` | `1e-12` | Minimum time step [s]. |
| `dtmax` | `Inf` | Maximum time step [s]. |

To use a different stiff integrator, pass its constructor (the linear solver and
`AutoFiniteDiff` autodiff are wired in for you) or a fully built instance:

```julia
using OrdinaryDiffEq: Rodas5P, FBDF

# Constructor form — linear_solver is paired in automatically:
sol = solve(problem; integrator = Rodas5P, linear_solver = SparspakFactorization())

# Fully constructed form (do not also pass linear_solver):
sol = solve(problem; integrator = FBDF(linsolve = KLUFactorization()))
```

Loosening tolerances (`reltol = 1e-4`) speeds up exploratory runs; tighten them
for convergence and validation studies (Technical Reference §18.5). Tolerance
choice does not affect mass/energy *bookkeeping* closure — that closes to the
integrator tolerance by construction (Technical Reference §16) — but it does
affect physical accuracy.

---

## 11.6 Loading and using extensions

Optional functionality is delivered through Julia's package-extension system
(Julia 1.9+). The hooks (`klu_factorization`, `pardiso_factorization`,
`amg_preconditioner`, `forward_sensitivity`, `adjoint_sensitivity`,
`verify_mixing_derivatives`) are **exported stubs** in core Pyrolysis.jl. Each
stub errors with a load hint until you `using` the corresponding optional
package, at which point the extension attaches a concrete method. No extension
depends on another; all are independent peers.

| Extension | Trigger (`using ...`) | Attaches | Optional dep(s) |
|---|---|---|---|
| `PyrolysisKLUExt` | `KLU` | `klu_factorization()` | `KLU` |
| `PyrolysisPardisoExt` | `Pardiso` | `pardiso_factorization()` | `Pardiso` |
| `PyrolysisAMGExt` | `AlgebraicMultigrid` | `amg_preconditioner(A)` | `AlgebraicMultigrid` |
| `PyrolysisSciMLSensitivityExt` | `SciMLSensitivity, Zygote` | `forward_sensitivity`, `adjoint_sensitivity` | `SciMLSensitivity`, `Zygote` |
| `PyrolysisSymbolicsExt` | `Symbolics` | `verify_mixing_derivatives` | `Symbolics` |
| `PyrolysisPlotsExt` | `RecipesBase`/`Plots` | `plot(sol; kind = ...)` recipes | `RecipesBase` |

### 11.6.1 Linear-solver extensions: KLU, Pardiso, AMG

These three extensions each expose a vendor solver under a stable
Pyrolysis-level name so the benchmark harness and out-of-tree code can request
it without core Pyrolysis taking a hard dependency.

```julia
using Pyrolysis
using KLU                              # triggers PyrolysisKLUExt

ls  = Pyrolysis.klu_factorization()    # -> LinearSolve.KLUFactorization()
sol = solve(problem; linear_solver = ls)
```

```julia
using Pyrolysis
using Pardiso                          # triggers PyrolysisPardisoExt

# Returns LinearSolve.MKLPardisoFactorize(); for users with an Intel MKL license.
ls  = Pyrolysis.pardiso_factorization()
sol = solve(problem; integrator = KenCarp4(linsolve = ls))
```

```julia
using Pyrolysis
using AlgebraicMultigrid               # triggers PyrolysisAMGExt
using SparseArrays, LinearAlgebra

# Smoothed-aggregation AMG preconditioner for a sparse matrix A, for use as a
# preconditioner inside a Krylov linsolve.
A    = sprand(100, 100, 0.05) + 100I   # illustrative sparse matrix
prec = Pyrolysis.amg_preconditioner(A)
```

Note that `KLUFactorization` is already reachable as
`LinearSolve.KLUFactorization` (LinearSolve is a hard dependency); the KLU
extension simply confirms KLU itself is loaded and exposes it under the stable
name `Pyrolysis.klu_factorization`. The AMG preconditioner is meant for the
preconditioner slot of a Krylov solver, not as a `linear_solver` on its own.

**Stub behavior without the trigger.** If you call any of these before loading
the optional package, you get a `MethodError` pointing at the unloaded
extension:

```julia
using Pyrolysis
Pyrolysis.klu_factorization()   # MethodError — `using KLU` not yet loaded
```

The fix is always the same: add the corresponding `using` line.

### 11.6.2 Sensitivity extension

`forward_sensitivity` and `adjoint_sensitivity` are stubs until you load
`SciMLSensitivity` and `Zygote`, which activate `PyrolysisSciMLSensitivityExt`.
Both take a `p_inject` closure that builds a fresh problem from a parameter
vector `θ` (it must not mutate the base problem). The closure rebuilds the
problem from the base problem's mesh and BCs with a material that carries the
injected parameters; the exact way you parameterize the material is up to you
(Chapter 13 develops realistic `p_inject` closures in full):

```julia
using Pyrolysis
using SciMLSensitivity, Zygote          # triggers PyrolysisSciMLSensitivityExt

# p_inject(base, θ) -> fresh PyrolysisProblem with θ injected.
# `make_material(θ)` is your own function that builds a Material from the
# parameter vector (e.g. setting reaction A, E). The fresh problem reuses the
# base problem's mesh, BCs, time span, and initial conditions.
function p_inject(base, θ)
    mat = make_material(θ)
    return PyrolysisProblem(base.mesh, mat, base.bc_set, base.tspan;
                            T_initial = base.T_initial,
                            output_times = base.output_times)
end

# Forward mode (ForwardDiff) — small parameter counts:
θ = [1e13, 2.0e5]
output, ∂out_∂θ = forward_sensitivity(base_problem, θ, p_inject;
                                      output_fn = sol -> maximum(sol.surface_T))

# Adjoint mode (Zygote + SciMLSensitivity) — large parameter counts, scalar loss:
loss, grad = adjoint_sensitivity(base_problem, θ, p_inject, sol -> sum(sol.surface_T);
                                 sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()))
```

This chapter only covers the *loading* mechanics; the full sensitivity workflow
(forward vs. adjoint choice, inverse analysis, Sobol global sensitivity) is the
subject of Chapter 13, and the theory is in Technical Reference §17.

### 11.6.3 Symbolics verification extension (test-only)

`verify_mixing_derivatives` is a test-only stub activated by `using Symbolics`.
It compares the package's analytical mixing-rule derivatives against
Symbolics-derived ground truth. Production code never calls it; it exists to
catch drift bugs in the mixing-rule Jacobian terms.

```julia
using Pyrolysis
using Symbolics                         # triggers PyrolysisSymbolicsExt

# Returns a per-rule max relative error (NaN for BRUGGEMAN, which has no closed form).
report = verify_mixing_derivatives(material; rtol = 1e-8, T_test = 600.0)
```

---

## 11.7 Debugging Jacobian correctness with `verify_jacobian_at`

When a custom or approximate Jacobian configuration produces suspicious Newton
behavior (excessive rejections, divergence, or wrong steady state), compare it
against the dense ForwardDiff reference at a state point. There are two ways: a
one-shot call, and an in-solve callback.

### 11.7.1 One-shot verification

`verify_jacobian_at` takes an active spec, a `(def, ws)` pair, a state `u`, and
a time `t`, and reports the deviation from the `DenseAD_ForwardDiff` reference.
Its signature is `verify_jacobian_at(spec, def, ws, u, t; reference =
DenseAD_ForwardDiff(), rtol = 1e-6, atol = 1e-10)`:

```julia
using Pyrolysis

# Build the (def, ws) pair for the mode you intend to solve in.
def = to_problem_def(problem; use_ale = true)
ws  = build_workspace(problem, def)
u   = create_state_vector(ws.layout, ws.mesh)

spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)

report = verify_jacobian_at(spec, def, ws, u, 0.0;
                            reference = DenseAD_ForwardDiff(),
                            rtol = 1e-6, atol = 1e-10)

report.passed                      # Bool — strict: per-entry |ΔJ_ij| ≤ atol + rtol·max(|J_ref[i,j]|, row max);
                                   #        masked: max_abs_dev ≤ atol || max_rel_dev ≤ rtol
report.max_abs_dev                 # largest element-wise deviation under the mode's rule
report.max_rel_dev
report.worst_index                 # (row, col) of the worst deviation
report.mode                        # :strict (exact backends) or :masked (approximate)
report.max_abs_dev_in_support      # in-support disagreement
report.max_abs_dev_outside_support # reference magnitude the backend's support omits
```

Comparison mode is chosen automatically by the backend's `is_approximate`
trait (§11.3.3):

- **Strict / support-complete** (exact backends, e.g. `DirectSolve`,
  `StructuredSolve`): every entry of the dense reference is compared against the
  active backend; a non-zero reference entry that the backend's support omits
  surfaces as a deviation, and its magnitude appears in
  `max_abs_dev_outside_support`. The strict criterion is per-entry with a
  row-magnitude floor, so a missing coupling in a small-magnitude row cannot
  hide behind the global matrix maximum. This is the regime exact backends
  must satisfy.
- **Masked** (`ApproxSolve`): entries outside the backend's declared support are
  zeroed in the diff; the omitted-entry magnitudes are reported in
  `max_abs_dev_outside_support` but not asserted. Use `is_approximate(spec.backend)`
  to confirm which mode you are in.

`verify_jacobian_at` reports deviations; it does not throw. You inspect
`report.passed` and the per-entry fields to localize a discrepancy.

### 11.7.2 In-solve verification callback

To sample Jacobian agreement at several times during an actual solve, pass the
`verify_jacobian` keyword as a `NamedTuple` with keys drawn from `{:rtol, :atol,
:sample_times}`. This installs a callback that compares the active backend
against the dense reference each time the integrator crosses one of the listed
times and **warns** (does not abort) on a mismatch:

```julia
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem;
    jacobian        = spec,
    use_ale         = true,
    verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 1.0, 10.0]))
```

The callback only runs when a pluggable `jacobian` is supplied (there is nothing
to verify against the integrator's own finite-difference Jacobian). It warns
rather than throws so a single transient outlier does not abort an otherwise
healthy production run.

---

## 11.8 The pure-residual layer: building a custom solver loop

Beneath `solve`, Pyrolysis.jl exposes a decoupled residual API so you can drive
the time integration yourself, build a non-standard solver, or feed the residual
to external tooling. The three pieces are all top-level exports:

- `ProblemDef` — an **immutable** frozen problem definition (material, mesh
  topology, BCs, initial conditions, simulation-mode traits, options). Safe to
  share across threads and sweep workers.
- `Workspace` — a **mutable** bundle of per-step caches, work arrays, and the
  geometry-from-state machinery. Allocate **one per worker**.
- `residual!(du, u, def, t, ws)` — the pure ODE right-hand side. Note the
  argument order: `du`, `u`, `def`, `t`, `ws`.

The builders `to_problem_def`, `build_workspace`, and `create_state_vector`
turn a `PyrolysisProblem` into this lower-level form:

```julia
using Pyrolysis
import OrdinaryDiffEq                                  # from your environment

# 1. Freeze the problem definition (pick the mode flags you want).
def = to_problem_def(problem;
    use_ale         = false,
    radiation_model = NO_RADIATION,
    energy_form     = :standard)

# 2. Build a mutable workspace (one per worker/thread).
ws = build_workspace(problem, def)

# 3. Initial state vector, packed [T, ξ, z (ALE), χ (WithChi)].
u0 = create_state_vector(ws.layout, ws.mesh)

# 4. Wrap residual! in the (du, u, p, t) shape OrdinaryDiffEq expects.
#    Note the residual! argument order is (du, u, def, t, ws).
rhs!(du, u, p, t) = residual!(du, u, def, t, ws)

# 5. Hand-roll your own integration.
prob = OrdinaryDiffEq.ODEProblem(rhs!, u0, problem.tspan)
raw  = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.KenCarp4(autodiff = OrdinaryDiffEq.AutoFiniteDiff()))
```

`StateLayout` (also exported) describes the packing and provides accessors so
you can read fields out of a raw state vector without re-deriving the offsets:

```julia
layout = ws.layout

# Index helpers and block accessors:
state_length(layout)                       # total state-vector length
get_temperature(u, layout, i)              # T in cell i
get_concentration(u, layout, i, j)         # ξ of component j in cell i (cell index first)
get_node_positions(u, layout)              # z nodes (ALE only)
get_thickness_from_state(u, layout)        # current sample thickness
block_view(u, layout, ...)                 # a view of one block

# For variable-cross-section runs:
column_chi_bar_from_state(u, layout, material)        # column-average χ̄
cross_section_area_from_state(mesh, material, u, layout)
```

You can also attach a structured Jacobian to your hand-rolled loop. Build a spec
from the same `(def, ws)` pair and call `compute_jacobian!`:

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
J    = zeros(state_length(layout), state_length(layout))
compute_jacobian!(J, spec, u0, def, 0.0, ws)
```

When driving the solver yourself, remember the contract that `solve` normally
enforces for you: **the residual and Jacobian assembly are mesh-read-only**. The
mutable mesh state is only meant to be synchronized to the integrator's `u`
between accepted steps (this is what the built-in `StepUpdateCallback` does
inside `solve`). If your custom loop relies on the mesh geometry being current
for post-processing or AMR, you must call
`update_mesh_from_state!(mesh, u, layout)` between steps yourself.

---

## 11.9 The `Internal` umbrella and importing trait subtypes

The top-level `Pyrolysis` namespace deliberately exports only a curated public
surface. Two escape hatches give you the rest.

### 11.9.1 `Pyrolysis.Internal`

`Pyrolysis.Internal` re-exports **every** symbol exported by the core
submodules (Materials, BoundaryConditions, Geometry, Properties, Physics,
Adaptivity, Discretization, Problem, Residual, Jacobian, Diagnostics, Solver) in
one namespace:

```julia
using Pyrolysis.Internal   # everything the core submodules export, in one place
```

Use this when you need an internal building block — a mixing function, a flux
kernel, a cache type — that is not on the curated top-level list. Caveat:
symbols reached through `Internal` are **not** part of the stable public API and
may move or change between minor releases; prefer top-level exports in
production code. (The `Experimental` namespace — P1 radiation, h-AMR — is *not*
included in `Internal`; reach it explicitly as `Pyrolysis.Experimental.X`.)

### 11.9.2 Importing the simulation-mode trait subtypes

The `SimulationMode` *abstract* trait is exported, but the concrete subtypes —
`Eulerian`, `ALE`, `Fickian`, `DarcyFick`, `NoChi`, `WithChi`, `NoRadiation`,
`SurfaceAbsorption`, `BeerLambert` — are **not** re-exported at top level, to
keep the public surface clean and avoid shadowing generic names like `ALE`. If
you need to dispatch on or construct one (e.g. for a custom residual
specialization), reach it through the `Problem` submodule:

```julia
# Fully qualified:
mode_tag = Pyrolysis.Problem.Eulerian()

# Or import the ones you need explicitly:
using Pyrolysis.Problem: Eulerian, ALE, Fickian, DarcyFick, WithChi
```

In normal usage you never construct these directly — `solve` and
`to_problem_def` select the right mode bundle from your keyword flags
(`use_ale`, `radiation_model`, the presence of permeability for `DarcyFick`, a
lateral-shrinkage law for `WithChi`). You only reach for the trait subtypes when
building custom dispatch at the residual layer (see Chapter 7 for how the trait
axes map to the four residual paths, and Technical Reference §15.1).

> Note: P1 radiation (`P1_QUASI_STEADY`) is **not** part of the unified
> residual; it lives in `Pyrolysis.Experimental` with its own residual wrapper
> and is rejected by `solve` (`radiation_model = P1_QUASI_STEADY` raises an
> error). The supported radiation models for `solve` are `NO_RADIATION`,
> `SURFACE_ABSORPTION`, and `BEER_LAMBERT` (Technical Reference §8.12).

---

## 11.10 Quick configuration recipes

```julia
# (a) Default, simplest correct solve — Eulerian, no depletion.
sol = solve(problem)

# (b) Explicit exact structured Jacobian (Eulerian).
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)

# (c) Large ALE, fast approximate Jacobian (~8× faster assembly).
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)

# (d) Large ALE, exact but cheap assembly via Richardson.
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)

# (e) ILU-GMRES linear solve for a large W matrix.
sol = solve(problem;
    integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
    jacobian   = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale    = true)

# (f) Depletion tracking — no pluggable Jacobian; solve picks the resize-safe dense LU.
sol = solve(problem;
    use_ale          = true,
    handle_depletion = true)

# (g) Vendor solver via extension.
using KLU
sol = solve(problem; linear_solver = Pyrolysis.klu_factorization())

# (h) Verify a custom backend against the dense reference during a solve.
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true,
             verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 5.0]))
```

---

## 11.11 Cross-references

- Structured-Jacobian architecture, the implicit-RK `W = I/(γ_RK Δt) − J`
  matrix, per-operator Jacobian formulas, and the prefix/suffix/rank-1
  couplings: **Technical Reference §15**.
- Surface-depletion handling and the depletion ↔ pluggable-Jacobian
  incompatibility: **Chapter 10** and **Technical Reference §14.7**.
- Building a problem and selecting simulation modes: **Chapter 7**.
- The `solve` keyword reference and post-processing: **Chapter 8** and
  **Chapter 9**.
- Sensitivity workflows that build on the sensitivity extension: **Chapter 13**
  and **Technical Reference §17**.
- The complete API and default-value tables: **Chapter 16**.


---

# 12. Parameter Sweeps

A *parameter sweep* runs the same base problem many times, each time with one or
more inputs changed, and reduces every run to a small summary. Pyrolysis.jl ships
a single harness for this — `parameter_sweep(problem, grid; modify, extract)` —
that materializes a fresh problem at each grid point, solves it, and collects the
result. This chapter explains how to build the grid, write the `modify` and
`extract` closures correctly (Pyrolysis types are immutable, so you *rebuild*
rather than mutate), run the sweep in parallel with `Distributed.jl`, reuse a
warm start between grid points, and tabulate the output. Sweeps give you discrete
maps of how outputs respond to inputs; for *derivatives* of outputs with respect
to inputs (gradients), use the sensitivity machinery in Chapter 13 instead. For
the theory behind global/variance-based (Sobol) sensitivity built on top of this
harness, see Technical Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the bottom/substrate
> (boundary id 2), `z = L` is the top/exposed surface (boundary id 1); heat
> enters from the top. Arrhenius parameters are the pre-exponential `A_i` [s⁻¹],
> the activation energy `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹]
> with the **`h > 0` endothermic** convention (see Technical Reference §2 and §6).

---

## 12.1 The `parameter_sweep` function

The public entry point is exported from `Pyrolysis` and has this signature:

```julia
parameter_sweep(
    problem,
    grid;
    modify::Function = _default_sweep_modify,   # (base_problem, point) -> new_problem
    extract::Function = identity,               # solution -> summary
    parallel::Bool = nworkers() > 1,            # use Distributed.pmap when true
    solve_kwargs...,                            # forwarded to every solve() call
) -> Vector
```

The body is intentionally tiny — the harness only orchestrates; you supply the
two closures that define *what changes* (`modify`) and *what you keep*
(`extract`):

```julia
work = function (point)
    new_problem = modify(problem, point)
    sol = solve(new_problem; progress = false, solve_kwargs...)
    return extract(sol)
end

pts = collect(grid)
return parallel ? pmap(work, pts) : map(work, pts)
```

Key facts that follow directly from this implementation:

| Aspect | Behavior |
|---|---|
| Return value | A `Vector`, one entry per grid point, in `collect(grid)` order. Each entry is whatever `extract` returns (the full `PyrolysisSolution` by default). |
| `progress` | Always forced to `false` inside the sweep (a progress bar per grid point would be noise); you cannot re-enable it via `solve_kwargs`. |
| Base problem | **Never mutated.** Each grid point gets a fresh problem from `modify`. |
| Grid materialization | `collect(grid)` runs first, so lazy iterators (e.g. `Iterators.product`) are realized into a concrete collection before dispatch. |
| Failures | Not retried. A `solve` that errors propagates out of `pmap`/`map` unless you guard inside `extract` (see §12.7). |
| `solve_kwargs` | Forwarded verbatim to each `solve()`. This is where `abstol`, `reltol`, `use_ale`, `radiation_model`, `dtmin`, etc. go (see Chapter 8 for the full `solve` option list). |

The default `modify` is a guard that throws a helpful `ArgumentError` — you must
supply your own:

```julia
function _default_sweep_modify(_, _)
    throw(ArgumentError(
        "parameter_sweep requires a `modify::(problem, point) -> problem` " *
        "closure that materialises a fresh PyrolysisProblem from each grid " *
        "point. See the docstring for an example."))
end
```

The default `extract` is `identity`, i.e. each entry of the returned vector is the
entire `PyrolysisSolution`. For anything beyond a handful of points you should
pass an `extract` that returns a small `NamedTuple` (see §12.4) — full solutions
are expensive to ship across worker boundaries and hold a lot of memory.

---

## 12.2 Building the grid

The `grid` argument is **any iterable of parameter points**. A "point" is whatever
your `modify` closure knows how to consume — a `NamedTuple`, a `Vector`, a `Dict`,
a scalar — the harness does not interpret it. Common shapes:

**1-D sweep (a vector of scalars).** Sweep a single quantity such as the incident
heat flux:

```julia
fluxes = [15e3, 25e3, 35e3, 50e3, 75e3]   # W/m², a Vector{Float64}
```

**N-D Cartesian grid (a vector of `NamedTuple`s).** A comprehension over two axes
produces a `Matrix{NamedTuple}`; `collect`/`pmap` flatten it column-major:

```julia
grid = [(A = a, E = e)
        for a in (1e12, 1e13, 1e14),   # rows
            e in (1.5e5, 2.0e5)]        # columns
# grid isa Matrix{@NamedTuple{A::Float64, E::Float64}}, size (3, 2)
```

**Lazy product (materialized by the harness).** For higher dimensions a
comprehension allocates the whole array up front; `Iterators.product` defers that
and is `collect`-ed internally:

```julia
grid = Iterators.product((1e12, 1e13, 1e14), (1.5e5, 2.0e5))
# each point is a Tuple (a, e); your modify must index it as point[1], point[2]
```

**Named axes via `Dict`.** Useful when `modify` looks up parameters by name and
not every point varies the same set:

```julia
grid = [Dict(:A => 1e13, :E => 1.8e5),
        Dict(:A => 5e13, :E => 1.9e5)]
```

Pick the shape that makes `modify` simplest to write. Because the result vector is
in `collect(grid)` order, keep the grid object around so you can pair each result
back with the input that produced it (see §12.8).

---

## 12.3 The `modify` closure — rebuild, do not mutate

`modify(base_problem, point)` must **return a fresh `PyrolysisProblem`** for the
given point. Pyrolysis's `Material`, `Reaction`, component, and `PyrolysisProblem`
types are immutable `struct`s with no in-place setters and **no copy
constructors** — you rebuild the pieces you want to change and reuse the rest by
reading the struct's fields. The relevant fields:

- `PyrolysisProblem`: `mesh`, `material`, `bc_set`, `T_initial`, `ξ_initial`,
  `tspan`, `dt_initial`, `output_times`.
- `Material`: `name`, `components` (a tuple), `reactions` (a tuple), `mixing`,
  `lateral_shrinkage_law`.
- `Reaction`: `name`, `reactants`, `products`, `A`, `E`, `heat`, `T_min`, `T_max`.

The constructors you rebuild with are the standard keyword constructors (Chapters
3–5):

- `PyrolysisProblem(; mesh, material, bc_set, T_initial, ξ_initial, tspan,
  dt_initial, output_times)`.
- `Material(; name, components, reactions, mixing, lateral_shrinkage_law)`.
- `Reaction(name, reactant => product; A, E, h, n, T_min, T_max)` and
  `Reaction(name, reactant => (p1, p2); A, E, h, yields, n, T_min, T_max)`.

### 12.3.1 Sweeping Arrhenius parameters (A, E)

The cleanest pattern: rebuild the single reaction you are sweeping, splice it back
into a fresh reactions tuple, rebuild the material, and rebuild the problem. The
example below assumes the swept reaction is the first one and reuses everything
else from the base material. Note the `Reaction` fields: `reactants` is a tuple of
`(component_idx, order)` pairs and `products` is a tuple of
`(component_idx, mass_stoich)` pairs (see Chapter 4):

```julia
function modify_AE(prob, point)
    mat = prob.material
    old_rxn = mat.reactions[1]               # the reaction we are sweeping

    # Rebuild reaction 1 with the swept A, E; reuse its other parameters.
    reactant_idx = old_rxn.reactants[1][1]   # (component_idx, order)
    order        = old_rxn.reactants[1][2]
    product_idx  = old_rxn.products[1][1]     # (component_idx, mass_stoich)
    new_rxn = Reaction(
        old_rxn.name,
        reactant_idx => product_idx;
        A     = point.A,                     # swept pre-exponential [s⁻¹]
        E     = point.E,                     # swept activation energy [J/mol]
        h     = old_rxn.heat,                # reuse heat-of-reaction property
        n     = order,
        T_min = old_rxn.T_min,
        T_max = old_rxn.T_max,
    )

    # Splice the rebuilt reaction back into the tuple. Tuples are immutable, so
    # construct a new one (here: replace index 1, keep the rest):
    new_reactions = (new_rxn, mat.reactions[2:end]...)

    new_mat = Material(
        name                  = mat.name,
        components            = mat.components,        # reused unchanged
        reactions             = new_reactions,
        mixing                = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
    )

    return PyrolysisProblem(
        mesh         = prob.mesh,
        material     = new_mat,
        bc_set       = prob.bc_set,
        T_initial    = prob.T_initial,
        ξ_initial    = prob.ξ_initial,
        tspan        = prob.tspan,
        dt_initial   = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

A few notes:

- `old_rxn.heat` is already a property function (an `AbstractPropertyFunction`),
  so passing it straight back as `h` is fine — the `Reaction` constructor accepts
  either a number or a property function and runs `to_property` on it internally.
- If the reaction has **two products**, use the two-product `Reaction`
  constructor and pass the yields you read back from the old reaction:
  `reactant_idx => (old_rxn.products[1][1], old_rxn.products[2][1])` with
  `yields = (old_rxn.products[1][2], old_rxn.products[2][2])`.
- Mass balance is re-validated when the reaction is reconstructed
  (`validate_mass = true` by default), so a malformed point fails fast.

### 12.3.2 Sweeping the heat of reaction (h)

Heat of reaction is a property of the reaction, swept the same way. Remember the
sign convention: **`h > 0` is endothermic** (see Technical Reference §6):

```julia
function modify_h(prob, point)
    mat = prob.material
    old = mat.reactions[1]
    new_rxn = Reaction(
        old.name,
        old.reactants[1][1] => old.products[1][1];
        A = old.A, E = old.E,
        h = point.h,                         # [J/kg], h > 0 endothermic
        n = old.reactants[1][2],
        T_min = old.T_min, T_max = old.T_max,
    )
    new_mat = Material(
        name = mat.name, components = mat.components,
        reactions = (new_rxn, mat.reactions[2:end]...),
        mixing = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
    )
    return PyrolysisProblem(
        mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

### 12.3.3 Sweeping a boundary condition (incident heat flux)

To sweep an *applied* quantity rather than a material property, rebuild only the
`bc_set`. A `BoundaryConditionSet` exposes two `NamedTuple` fields, `thermal` and
`mass`, each keyed by `:top` and `:bottom`. The example below sweeps the incident
radiant flux delivered by a `RadiativeFluxBC` at the top surface (the cone-heater
flux), keeping the bottom thermal BC and both mass BCs unchanged:

```julia
function modify_flux(prob, q_incident)            # point is a scalar flux [W/m²]
    bcs = prob.bc_set
    new_top_thermal =
        RadiativeFluxBC(flux = q_incident, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = 330.0) +
        ConvectiveBC(h = 7.2, T_inf = 330.0)

    new_bc_set = BoundaryConditionSet(
        thermal = (top = new_top_thermal, bottom = bcs.thermal.bottom),
        mass    = (top = bcs.mass.top,    bottom = bcs.mass.bottom),
    )

    return PyrolysisProblem(
        mesh = prob.mesh, material = prob.material, bc_set = new_bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

The `+` operator on thermal BCs builds a `CombinedThermalBC` (Chapter 5). Each
grid point therefore gets its own combined boundary condition with the swept
incident flux while re-radiation and convection stay fixed. Because the swept
quantity is the incident flux of a `RadiativeFluxBC`, run the sweep with a
matching `radiation_model` in `solve_kwargs` (e.g. `SURFACE_ABSORPTION`); see
Chapter 8 and Technical Reference §8 for how the radiation model interprets that
flux.

> **Caveat — keep `modify` self-contained for parallel runs.** When you run with
> workers (§12.5) the closure body executes on a remote process. Anything it
> references — helper functions, constants — must be defined on the workers too
> (use `@everywhere`). Reading fields off the `base_problem` is safe because the
> base problem is shipped to each worker by `pmap`.

---

## 12.4 The `extract` closure — keep results small

`extract(solution)` reduces a `PyrolysisSolution` to whatever you want to keep.
The default is `identity` (the full solution). For real sweeps, return a small
`NamedTuple` of scalars or short vectors. The fields available on a
`PyrolysisSolution` (Chapter 9) include:

| Field | Meaning |
|---|---|
| `sol.t` | Saved times [s] |
| `sol.T` | Temperature field, `(n_cells, n_steps)` [K] |
| `sol.ξ` | Concentrations, `(NC, n_cells, n_steps)` [kg/m³] |
| `sol.surface_T` | Surface temperature at each saved time [K] |
| `sol.mass_loss_rate` | MLR per *current* area `A(t)` (top + permeable-bottom) [kg/(m²·s)] |
| `sol.retcode` | Integrator return code (`:Success`, `:Terminated`, …) |
| `sol.solve_time` | Wall-clock solve time [s] |
| `sol.extras.total_mass` | Domain mass per initial area [kg/m²] |
| `sol.extras.thickness_history` | Sample thickness over time [m] |
| `sol.extras.MLR_gauge` | Gauge MLR (`MLR_local · A(t)/A_0`) [kg/(m²·s)] |
| `sol.extras.gas_generation_rate` | Volumetric gas production rate [kg/(m²·s)] |

A typical scalar-summary `extract`:

```julia
extract_summary(sol) = (
    retcode    = sol.retcode,
    peak_T     = isempty(sol.surface_T)      ? NaN : maximum(sol.surface_T),
    peak_MLR   = isempty(sol.mass_loss_rate) ? NaN : maximum(sol.mass_loss_rate),
    final_mass = isempty(sol.extras.total_mass) ? NaN : sol.extras.total_mass[end],
    solve_time = sol.solve_time,
)
```

If you need a *curve* (say the surface-temperature trace) but not the full field,
return just that vector plus its time base:

```julia
extract_curve(sol) = (t = copy(sol.t), surface_T = copy(sol.surface_T))
```

> Use `copy` when returning arrays from `extract` so the summary does not alias
> the (about-to-be-garbage-collected) solution's internal buffers.

Always carry `sol.retcode` into your summary so you can detect grid points that
did not finish cleanly — `:Success` means the integrator reached `t_end`,
`:Terminated` means a callback stopped it (e.g. ALE thickness termination), and
codes like `:MaxIters` or `:DtLessThanMin` flag a failed solve.

---

## 12.5 Parallel execution with `Distributed.jl`

`parallel` defaults to `true` whenever more than one worker is available
(`nworkers() > 1`); otherwise the sweep runs sequentially with `map`. When
parallel, the harness uses `Distributed.pmap`, which load-balances grid points
across workers. The thread-safety model is built into the design:

- The pure residual `residual!(du, u, def, t, ws)` is **immutable** and reused;
  the mutable `Workspace` is allocated lazily inside each worker's `solve()`, so
  each grid point's solve has its own scratch state — no contention.
- The immutable `ProblemDef`/`Material` data is shipped to workers by `pmap`.

A complete parallel sweep:

```julia
using Distributed
addprocs(4)                       # start 4 worker processes
@everywhere using Pyrolysis       # load the package on every worker

# Make the closures (and any helpers/constants they use) available on workers:
@everywhere function modify_AE(prob, point)
    mat = prob.material
    old = mat.reactions[1]
    new_rxn = Reaction(old.name, old.reactants[1][1] => old.products[1][1];
                       A = point.A, E = point.E, h = old.heat,
                       n = old.reactants[1][2], T_min = old.T_min, T_max = old.T_max)
    new_mat = Material(name = mat.name, components = mat.components,
                       reactions = (new_rxn, mat.reactions[2:end]...),
                       mixing = mat.mixing,
                       lateral_shrinkage_law = mat.lateral_shrinkage_law)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
                            T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
                            tspan = prob.tspan, dt_initial = prob.dt_initial,
                            output_times = prob.output_times)
end

@everywhere extract_summary(sol) = (
    retcode = sol.retcode,
    peak_T  = isempty(sol.surface_T) ? NaN : maximum(sol.surface_T),
)

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

results = parameter_sweep(
    problem, grid;
    modify  = modify_AE,
    extract = extract_summary,
    abstol  = 1e-8, reltol = 1e-6,    # solve_kwargs forwarded to every solve()
)
```

Notes on the parallel path:

- Define `modify` and `extract` with `@everywhere` (or in a module loaded
  `@everywhere`). A bare `function` defined only on the head process will not be
  found on the workers and the sweep will fail with an `UndefVarError`.
- `parallel` is computed once at call time from `nworkers()`. If you call
  `parameter_sweep` before `addprocs`, it runs sequentially even if you add
  workers afterward. Force the path explicitly with `parallel = true` /
  `parallel = false` if you want to override the default.
- There is no benefit to more workers than grid points; `pmap` simply leaves the
  extra workers idle.
- The order of the returned vector matches `collect(grid)`, *not* the order in
  which workers finished — `pmap` preserves input order.

For a sequential sweep (debugging, single core, or small grids where worker
start-up cost dominates) either omit `addprocs` or pass `parallel = false`:

```julia
results = parameter_sweep(problem, fluxes; modify = modify_flux,
                          extract = extract_summary, parallel = false)
```

---

## 12.6 Warm starting across grid points

`solve` accepts `warm_start::PyrolysisSolution`, which seeds the integrator's
initial state with a previous solution's *final* state. On a sweep whose points
form a smooth continuation (e.g. a monotone heat-flux ramp), warm starting can cut
per-point solve time by giving the implicit integrator a better starting guess.

The compatibility rules are strict and enforced by a state-vector length check:
the component count `NC`, the mesh size, the ALE mode (`use_ale`), and the χ mode
must all be identical between the prior solve and the new one — otherwise `solve`
errors. (`tspan[1]` still sets the integrator's start time; warm starting just
reuses the prior final state as fresh initial data.) The prior solution must also
still carry its `raw_solution`, so warm-starting from a heavily reduced `extract`
output will not work — keep the full `PyrolysisSolution` for the run you chain
from.

Because `pmap` runs grid points concurrently and out of any meaningful order,
warm starting is **not compatible with `parallel = true`** — there is no
"previous" point on a parallel worker. Do warm-started sweeps **sequentially**,
chaining each solve to the last:

```julia
function warm_chained_sweep(problem, points; modify, extract, solve_kwargs...)
    results = Vector{Any}(undef, length(points))
    prev = nothing
    for (i, pt) in enumerate(points)
        prob_i = modify(problem, pt)
        sol = prev === nothing ?
              solve(prob_i; progress = false, solve_kwargs...) :
              solve(prob_i; progress = false, warm_start = prev, solve_kwargs...)
        results[i] = extract(sol)
        prev = sol      # chain the FULL solution into the next point
    end
    return results
end

# Sweep a heat-flux ramp, warm-starting each point from the previous one:
fluxes = collect(15e3:5e3:75e3)
results = warm_chained_sweep(problem, fluxes;
                             modify = modify_flux, extract = extract_summary,
                             abstol = 1e-8, reltol = 1e-6)
```

> **Accuracy caveat.** Warm starting changes only the *initial guess*, not the
> equations, so the converged trajectory of each point is unaffected to integrator
> tolerance. But it does couple grid points: a point solved from a warm start may
> take a slightly different step sequence than the same point solved cold. For
> sweeps feeding a sensitivity or calibration study where reproducibility of each
> point in isolation matters, prefer cold solves. See the warm-start discussion in
> Technical Reference §17 (sensitivity pitfalls).

---

## 12.7 Fault tolerance

A `solve` that throws aborts the whole sweep (the exception propagates out of
`pmap`/`map`). For long grids you usually want a failed point to record a marker
and let the rest continue. The harness always calls `solve` itself, so to own the
error handling completely, take over the per-point body: make `modify` a
pass-through and do the rebuild, solve, and reduction inside one `extract`-side
function guarded by `try`/`catch`. The cleanest way is to fold the work into
`extract` while letting `modify` build the problem normally and catching only the
solve inside a small wrapper of your own:

```julia
@everywhere function point_result(prob, point)
    try
        new_prob = modify_AE(prob, point)
        sol = solve(new_prob; progress = false, abstol = 1e-8, reltol = 1e-6)
        return (ok = sol.retcode === :Success, retcode = sol.retcode,
                peak_T = isempty(sol.surface_T) ? NaN : maximum(sol.surface_T))
    catch err
        return (ok = false, retcode = :Error, peak_T = NaN, error = string(err))
    end
end

# Drive it through your own map/pmap rather than the harness, so the try/catch
# wraps the solve:
using Distributed
pts = collect(grid)
results = (nworkers() > 1 ? pmap : map)(p -> point_result(problem, p), pts)
```

If you stick with `parameter_sweep` itself, remember its `extract` runs *after*
`solve` returns — so an `extract`-side `try` catches problems in your reduction
(e.g. `maximum` on an empty vector), but **not** an exception thrown inside
`solve`. To catch solve-time exceptions you must move the `solve` into your own
loop, as above.

---

## 12.8 Collecting and tabulating results

The result is a flat `Vector` in `collect(grid)` order. To turn it into something
you can plot or write out, zip it back against the grid.

**1-D sweep → paired vectors:**

```julia
fluxes  = [15e3, 25e3, 35e3, 50e3, 75e3]
results = parameter_sweep(problem, fluxes; modify = modify_flux,
                          extract = extract_summary, parallel = false)

peak_T = [r.peak_T for r in results]
for (q, pT) in zip(fluxes, peak_T)
    println("q = $(q/1000) kW/m²  ->  peak T_s = $(round(pT, digits=1)) K")
end
```

**2-D grid → reshape into a matrix.** Because a comprehension grid is a `Matrix`
and `pmap`/`map` flatten column-major, reshape the result back to the grid's
shape:

```julia
A_vals = (1e12, 1e13, 1e14)
E_vals = (1.5e5, 2.0e5)
grid   = [(A = a, E = e) for a in A_vals, e in E_vals]   # 3×2 Matrix

results = parameter_sweep(problem, grid; modify = modify_AE,
                          extract = extract_summary)

peak_T = reshape([r.peak_T for r in results], size(grid))  # 3×2 again
# peak_T[i, j] corresponds to (A = A_vals[i], E = E_vals[j])
```

**Write to a delimited file** for external post-processing:

```julia
using DelimitedFiles
open("sweep_results.csv", "w") do io
    println(io, "A,E,peak_T,final_mass,retcode")
    for (p, r) in zip(vec(grid), results)
        println(io, "$(p.A),$(p.E),$(r.peak_T),$(r.final_mass),$(r.retcode)")
    end
end
```

(If you prefer a `DataFrame`, build one from the same `zip` of grid points and
result `NamedTuple`s; that requires the `DataFrames.jl` package, which is not a
dependency of Pyrolysis.jl.)

---

## 12.9 Worked example: heat-flux series on a cone-calorimeter problem

This puts the pieces together on a self-contained PMMA cone problem (the setup is
the Chapter 2 / example 06 problem in condensed form). The sweep varies the
incident cone flux and records peak surface temperature, peak MLR, and final mass.

```julia
using Pyrolysis

# --- Base material: PMMA (virgin -> char + MMA gas) ---
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2T) : (1040.0 + 3.08T))
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4T))

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96),
        SolidComponent(:char,   ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96),
        GaseousComponent(:MMA,  M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, 1 => (2, 3);
                 A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0,
                 yields = (0.002, 0.998)),     # h > 0 endothermic
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# --- Mesh, BCs, and base problem ---
L, n_cells = 5.8e-3, 72
mesh = create_uniform_mesh(L, n_cells, Val(3);
                           T_initial = 300.0, ξ_initial = (1210.0, 0.0, 0.0))

base_bc = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 25e3, absorptivity = 0.96) +
              RadiativeBC(ε = 0.96, T_inf = 330.0) +
              ConvectiveBC(h = 7.2, T_inf = 330.0),
        bottom = RadiativeBC(ε = 0.95, T_inf = 307.0) +
                 ConvectiveBC(h = 4.0, T_inf = 307.0),
    ),
    mass = (
        top    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
        bottom = ImpermeableBC(),
    ),
)

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = base_bc,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, 600.0),
    output_times = collect(range(0.0, 600.0, length = 121)),
)

# --- The sweep: incident cone flux, 15..75 kW/m² ---
function modify_flux(prob, q_incident)
    bcs = prob.bc_set
    new_top =
        RadiativeFluxBC(flux = q_incident, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = 330.0) +
        ConvectiveBC(h = 7.2, T_inf = 330.0)
    new_bc_set = BoundaryConditionSet(
        thermal = (top = new_top, bottom = bcs.thermal.bottom),
        mass    = (top = bcs.mass.top, bottom = bcs.mass.bottom),
    )
    return PyrolysisProblem(
        mesh = prob.mesh, material = prob.material, bc_set = new_bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end

extract_summary(sol) = (
    retcode    = sol.retcode,
    peak_T     = maximum(sol.surface_T),
    peak_MLR   = maximum(sol.mass_loss_rate),
    final_mass = sol.extras.total_mass[end],
)

fluxes  = [15e3, 25e3, 35e3, 50e3, 75e3]
results = parameter_sweep(problem, fluxes;
                          modify = modify_flux, extract = extract_summary,
                          parallel = false,
                          radiation_model = SURFACE_ABSORPTION,  # interprets the cone flux
                          abstol = 1e-8, reltol = 1e-6)

for (q, r) in zip(fluxes, results)
    println("q = $(q/1000) kW/m²: retcode=$(r.retcode), ",
            "peak T_s=$(round(r.peak_T, digits=1)) K, ",
            "peak MLR=$(round(r.peak_MLR*1000, digits=3)) g/(m²·s)")
end
```

This sweep runs in fixed-mesh (Eulerian) mode because no `use_ale = true` is
passed in `solve_kwargs`. To sweep an ALE problem instead — for example to study
how the swept flux changes the final sample thickness — add the ALE options to
the sweep call and read the ALE-specific extras:

```julia
extract_ale(sol) = (
    retcode         = sol.retcode,
    peak_T          = maximum(sol.surface_T),
    final_thickness = sol.extras.thickness_history[end],
    peak_MLR_gauge  = maximum(sol.extras.MLR_gauge),  # gauge-normalized MLR
)

results_ale = parameter_sweep(problem, fluxes;
    modify = modify_flux, extract = extract_ale, parallel = false,
    radiation_model = SURFACE_ABSORPTION,
    use_ale = true, handle_depletion = true,
    min_thickness = 0.01e-3, depletion_threshold = 0.05, min_cells = 1,
    dtmin = 1e-14)
```

In variable-cross-section/ALE runs use `sol.extras.MLR_gauge` rather than
`sol.mass_loss_rate` when comparing across grid points: `mass_loss_rate` is per
*current* area `A(t)`, whereas `MLR_gauge = MLR_local · A(t)/A_0` is normalized to
the initial gauge area `A_0`, which is the apples-to-apples quantity for cone
comparisons (see Chapter 8 and Technical Reference §10). Note also that when
`handle_depletion = true`, the solver falls back to a finite-difference Jacobian
(a pluggable `jacobian` backend is incompatible with mid-solve state-vector
resizing); this is automatic and needs no action in the sweep, but it is why
ALE-with-depletion points are slower per solve (see Chapter 11).

---

## 12.10 Worked example: an activation-energy / pre-exponential grid

A 2-D kinetics grid maps how peak surface temperature and the timing of peak gas
release respond to the Arrhenius pair `(A_1, E_1)`. Reusing `modify_AE` from
§12.3.1:

```julia
A_vals = (5e12, 1e13, 5e13, 1e14)
E_vals = (1.7e5, 1.9e5, 2.1e5)
grid   = [(A = a, E = e) for a in A_vals, e in E_vals]   # 4×3 Matrix

extract_kinetics(sol) = (
    retcode = sol.retcode,
    peak_T  = maximum(sol.surface_T),
    # time of peak gas-generation rate as a crude "ignition/burnout" proxy:
    t_peak  = sol.t[argmax(sol.extras.gas_generation_rate)],
)

results = parameter_sweep(problem, grid;
                          modify = modify_AE, extract = extract_kinetics,
                          radiation_model = SURFACE_ABSORPTION)

peak_T = reshape([r.peak_T for r in results], size(grid))
t_peak = reshape([r.t_peak for r in results], size(grid))
# peak_T[i, j], t_peak[i, j] correspond to (A = A_vals[i], E = E_vals[j])
```

Run this one in parallel for the 12 points by `addprocs`-ing first and defining
`modify_AE`/`extract_kinetics` with `@everywhere` (§12.5). Because the grid is a
`Matrix`, the reshape recovers the `(A, E)` layout for heatmaps or contour plots.

---

## 12.11 From sweeps to sensitivities

A parameter sweep gives you a **discrete map** of output-vs-input: you sample the
input space and read off outputs. When you instead need the **derivative** of an
output with respect to inputs — for gradient-based calibration, uncertainty
propagation, or optimization — use the sensitivity entry points covered in
Chapter 13:

- `forward_sensitivity(base_problem, θ, p_inject; output_fn)` — ForwardDiff
  Jacobian of an output with respect to a small parameter vector `θ`.
- `adjoint_sensitivity(base_problem, θ, p_inject, loss_fn; sensealg)` — Zygote
  gradient of a scalar loss for large parameter counts / inverse problems.

These require loading the `PyrolysisSciMLSensitivityExt` extension
(`using SciMLSensitivity, Zygote`); without it the exported names are stubs that
error with a hint to load the extension. The `p_inject` closure plays exactly the
role `modify` plays here — it rebuilds a fresh problem from a parameter vector —
so a `modify` you write for a sweep is usually one rename away from a `p_inject`
you can hand to `forward_sensitivity` (the only difference is the parameter shape:
`p_inject` takes a numeric vector `θ`, `modify` takes an arbitrary grid point).

The sweep harness is also the substrate for **global / variance-based (Sobol)**
sensitivity: you sample the parameter space (the grid), run the sweep, and feed
the input/output table into a Sobol estimator. See the
`examples/sensitivity_analysis/sensitivity_sobol.jl` study and Technical Reference
§17.5 for that workflow.

---

## 12.12 Quick reference

| Item | Value / note |
|---|---|
| Function | `parameter_sweep(problem, grid; modify, extract, parallel, solve_kwargs...)` |
| `modify` default | guard that throws `ArgumentError` — you must supply one |
| `extract` default | `identity` (returns the full `PyrolysisSolution`) |
| `parallel` default | `nworkers() > 1` |
| `progress` | forced `false` inside the sweep |
| Return | `Vector` in `collect(grid)` order |
| Grid shapes | `Vector{<:NamedTuple}`, `Iterators.product`, `Vector`/`Dict`, comprehension `Matrix` |
| Rebuild, don't mutate | no copy constructors; reuse struct fields and call keyword constructors |
| Parallel safety | immutable residual reused; mutable `Workspace` allocated per worker `solve` |
| Define closures on workers | `@everywhere` for `modify`/`extract` and their helpers |
| Warm start | sequential only; `NC` / mesh / ALE / χ mode must match; prior `raw_solution` required (see §12.6) |
| Cross-grid MLR (ALE) | compare `sol.extras.MLR_gauge`, not `sol.mass_loss_rate` |
| Gradients instead of maps | Chapter 13 (`forward_sensitivity` / `adjoint_sensitivity`) |


---

# 13. Sensitivity Analysis

*Sensitivity analysis* asks how the outputs of a pyrolysis simulation change when
its inputs change. Pyrolysis.jl supports three complementary flavors: **local
forward** sensitivities (`∂output/∂θ` by automatic differentiation through the
whole solve), **local adjoint** sensitivities (`∂loss/∂θ` for a scalar loss, the
efficient mode for inverse problems with many parameters), and **global /
variance-based** (Sobol) sensitivities built on the parameter-sweep harness from
Chapter 12. This chapter shows how to use each, what can legitimately be a
parameter, how to wire your parameters into a fresh problem, and the numerical
pitfalls that make a sensitivity ill-defined. For the mathematics behind the
methods — Dual-number propagation, the augmented adjoint ODE, and the variance
decomposition — see Technical Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the bottom/substrate
> (boundary id 2, tag `:bottom`); `z = L` is the top/exposed surface (boundary id
> 1, tag `:top`); heat enters from the top. The Arrhenius parameters are the
> pre-exponential `A_i` [s⁻¹] (or [m³·kg⁻¹·s⁻¹] for higher-order reactions), the
> activation energy `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹] with
> the **`h > 0` endothermic** storage convention (Technical Reference §2 and §6).
> Mass concentration `ξ_j` [kg·m⁻³] is the primary species state.

---

## 13.1 Overview: when to use which method

| Method | Entry point | Returns | Cost scaling | Best when |
|---|---|---|---|---|
| **Local forward** | `forward_sensitivity` | output value + Jacobian `∂out/∂θ` | linear in `length(θ)` | few parameters (`≲ 20`), one or many outputs |
| **Local adjoint** | `adjoint_sensitivity` | scalar loss + gradient `∂loss/∂θ` | ~independent of `length(θ)` | many parameters, **one scalar** loss (calibration, optimization) |
| **Global (Sobol)** | `parameter_sweep` + a Sobol design | variance indices per parameter/group | `O(N_samples)` solves | uncertainty quantification, importance ranking over wide ranges |

Forward and adjoint are **local**: they linearize about a single parameter point
`θ`. Sobol is **global**: it explores a distribution of inputs and attributes
output variance to inputs, but it gives no derivative.

`forward_sensitivity` and `adjoint_sensitivity` are exported by `Pyrolysis` as
**stubs** — they error with a clear message until you load the
`PyrolysisSciMLSensitivityExt` extension (§13.6). The Sobol path needs no
extension; it is pure `parameter_sweep` plus your choice of sampling package.

---

## 13.2 The `p_inject` pattern: turning a vector into a problem

Both local methods share one idea. You hand them a **base problem**, a
**parameter vector `θ`**, and a **closure** that rebuilds a fresh
`PyrolysisProblem` from `(base_problem, θ)`:

```julia
p_inject(base_problem, θ) -> PyrolysisProblem
```

The extension differentiates the chain `θ -> output_fn(solve(p_inject(base, θ)))`
end-to-end. The whole solver (state cache, geometry, the boundary-temperature
Newton solve, the property derivatives) is type-generic in the scalar element
type, so when ForwardDiff passes a `Dual` vector the sensitivities propagate
through every stage. You never rewrite the residual.

> **Pyrolysis types are immutable: rebuild, do not mutate.** There is no
> `Material(existing; ...)` or `PyrolysisProblem(existing; ...)` copy
> constructor. `p_inject` must construct the components, reactions, material,
> mesh, BCs and problem **from scratch**, substituting the entries of `θ` where
> they belong — exactly like the `modify` closure for a parameter sweep
> (Chapter 12), and exactly like the `build_and_solve(θ)` function in
> `examples/sensitivity_analysis/sensitivity_douglas_fir.jl`. The base problem
> is only there to supply the parts of the configuration you are *not* varying.

A minimal, correct `p_inject` for a single charring reaction whose
pre-exponential `A`, activation energy `E`, and heat of reaction `h` we want to
differentiate:

```julia
using Pyrolysis

# θ = [A, E, h]
function p_inject(base, θ)
    A, E, h = θ[1], θ[2], θ[3]

    virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20)
    char   = SolidComponent(:char;   ρ = 200.0,  c = 1200.0, k = 0.10)
    gas    = GaseousComponent(:gas;   M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

    rxn = Reaction(:decomp, 1 => (2, 3);
                   A = A, E = E, h = h,        # ← the differentiated parameters
                   yields = (0.2, 0.8), n = 1.0)

    material = Material(
        name       = :charring_polymer,
        components = (virgin, char, gas),
        reactions  = (rxn,),
    )

    mesh = create_uniform_mesh(10.0e-3, 40, Val(3); cross_section_area = 1.0)

    bc_set = BoundaryConditionSet(
        thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0)
                            + ConvectiveBC(h = 10.0, T_inf = 300.0),
                   bottom = AdiabaticBC()),
        mass    = (top = DiffusiveMassBC(Inf), bottom = ImpermeableBC()),
    )

    return PyrolysisProblem(;
        mesh     = mesh,
        material = material,
        bc_set   = bc_set,
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1190.0, z -> 0.0, z -> 0.0],
        tspan    = (0.0, 300.0),
        output_times = collect(range(0.0, 300.0, length = 151)),
    )
end
```

A few rules that keep `p_inject` differentiable and consistent:

- **Place `θ` entries directly into the constructor calls.** The element type of
  `θ` flows into the `Reaction` and `SolidComponent` fields. Do not round, clamp,
  or `Float64(...)` a `θ` entry — that would break AD by stripping the `Dual`.
- **Everything not in `θ` is fixed.** Read fixed values from `base` if it is
  convenient, or just hard-code them as above.
- **Constructors that *must* receive `Float64`** (e.g. `tspan`, `output_times`,
  mesh thickness, `n_cells`) should be plain `Float64`/`Int` literals, never
  functions of `θ`. You differentiate physical parameters, not the time grid or
  the mesh count.
- **`yields` must sum to 1.0** for mass conservation (Technical Reference §6.5);
  `Reaction(...)` validates this at construction. If you make a yield a parameter
  you must rebuild the complementary yields so they still sum to 1.

### What can be a parameter `θ`

Anything that enters a constructor as a number can be a parameter. Common, useful
handles:

| Parameter | Where it lives | Notes |
|---|---|---|
| `A_i`, `E_i` | `Reaction(...; A=, E=)` | the dominant kinetic sensitivities |
| `h` (heat of reaction) | `Reaction(...; h=)` | scalar, tuple `(a,b)` → linear, or a closure |
| `ρ_j`, `c_p`, `k_j` | `SolidComponent(...; ρ=, c=, k=)` | pass a number, or wrap a `θ`-dependent `CallableProperty` |
| `λ_j`, `κ_j` | component `λ=`, `κ=` | gas-transfer / permeability for Darcy mode |
| `ε_j`, `α_j` | component `ε=`, `α=` | emissivity / mass-absorption |
| BC parameters | `HeatFluxBC`, `ConvectiveBC(h=,T_inf=)`, `RadiativeFluxBC(flux=)`, … | external flux, `h_conv`, `T_∞` |
| mixing `β` | `MixingSpec(WEIGHTED, β)` | conductivity / transport blend weight |

For a temperature-dependent property whose coefficients are parameters, build the
property as a closure over `θ` so the `Dual` propagates:

```julia
# θ = [a, b]; c_p(T) = a + b*T as a differentiable callable
cp = CallableProperty(T -> θ[1] + θ[2] * T)
virgin = SolidComponent(:virgin; ρ = 1190.0, c = cp, k = 0.20)
```

`Reaction` accepts a scalar `h` (a constant), a tuple `h = (a, b)` (linear in `T`,
`a + b*T`), or a property object such as `StateDependentProperty((T, ξ) -> ...)`
for moisture-dependent heats of sorption (Technical Reference §6.6). The
`build_and_solve(θ)` in `sensitivity_douglas_fir.jl` demonstrates the
state-dependent form for heat of sorption.

---

## 13.3 Forward-mode local sensitivity

`forward_sensitivity` differentiates an output (scalar or vector) with respect to
`θ` using `ForwardDiff.jacobian`. Its signature (attached by the extension to the
exported stub) is:

```julia
forward_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function;
    output_fn::Function = identity,   # solution -> output (scalar or vector)
    solve_kwargs...,                  # forwarded to every internal solve()
) -> (output, ∂output_∂θ)
```

- `output` is `output_fn(solve(p_inject(base, θ)))` evaluated at the base `θ`.
- `∂output_∂θ` is the Jacobian: if `output` is a length-`m` vector and `θ` has
  length `n`, the result is `m × n`. If `output_fn` returns a scalar, wrap it in a
  1-vector (ForwardDiff differentiates vector-valued functions).

`output_fn` reduces a `PyrolysisSolution` to the quantity you care about. Use the
solution's documented fields and accessors (Chapter 9): `sol.surface_T`,
`sol.mass_loss_rate`, `sol.T`, `sol.ξ`, `sol.extras.total_mass`,
`sol.extras.surface_mass_flux`, etc.

> **Cost.** ForwardDiff carries one extra "tangent" per parameter, so the solve
> cost scales roughly linearly in `length(θ)`. Forward mode is the right choice
> for a handful of parameters and is robust because it differentiates the exact
> code path the solver runs. For dozens of parameters and a single scalar loss,
> prefer adjoint mode (§13.5).

### Worked example: sensitivity of peak surface temperature to `(A, E, h)`

```julia
using Pyrolysis
using SciMLSensitivity, Zygote   # loads PyrolysisSciMLSensitivityExt

θ0 = [1.0e13, 1.88e5, 8.70e5]     # [A (1/s), E (J/mol), h (J/kg)]

# Scalar QoI wrapped in a 1-vector so ForwardDiff differentiates it.
peak_Ts(sol) = [maximum(sol.surface_T)]

# `base` is unused inside p_inject here, but it must still be passed
# (and could supply fixed pieces if you wire it in).
base = p_inject(nothing, θ0)

output, J = forward_sensitivity(base, θ0, p_inject;
                                output_fn = peak_Ts,
                                radiation_model = SURFACE_ABSORPTION,
                                abstol = 1e-9, reltol = 1e-7)

peak     = output[1]              # peak surface temperature [K] at θ0
dpeak_dA = J[1, 1]               # ∂(peak T_s)/∂A   [K·s]
dpeak_dE = J[1, 2]               # ∂(peak T_s)/∂E   [K·mol/J]
dpeak_dh = J[1, 3]              # ∂(peak T_s)/∂h   [K·kg/J]
```

`solve_kwargs` (here `radiation_model`, `abstol`, `reltol`) are forwarded
verbatim to each internal `solve` call; `progress = false` is forced internally.

### Differentiating a whole curve

`output_fn` may return a vector — e.g. the surface-temperature history sampled at
the problem's `output_times`. The Jacobian then has one row per saved time:

```julia
Ts_history(sol) = copy(sol.surface_T)             # length = number of saved times

output, J = forward_sensitivity(base, θ0, p_inject;
                                output_fn = Ts_history,
                                radiation_model = SURFACE_ABSORPTION)
# J[i, k] = ∂T_s(t_i)/∂θ_k
```

> **Keep the output dimension fixed.** `output_fn` must return a vector of the
> same length for every `θ` ForwardDiff probes. Returning `sol.surface_T`
> directly is fine because `saveat = output_times` is fixed and independent of
> `θ`. Returning something whose length depends on the trajectory (e.g. "number
> of time steps the integrator took", or — in ALE mode — a node array whose
> length changes after a depletion merge) will break the Jacobian assembly.

---

## 13.4 Manual forward sensitivity (without the extension)

`forward_sensitivity` is a thin convenience wrapper. If you would rather not pull
in `SciMLSensitivity`/`Zygote` just for forward mode, you can call `ForwardDiff`
yourself on the same `θ -> output_fn(solve(p_inject(base, θ)))` chain:

```julia
using Pyrolysis
using ForwardDiff

g(θ) = [maximum(solve(p_inject(base, θ);
                      progress = false,
                      radiation_model = SURFACE_ABSORPTION).surface_T)]

output = g(θ0)
J = ForwardDiff.jacobian(g, θ0)
```

This is exactly what the extension does internally; the extension just packages it
and provides the matching adjoint method. Forward mode needs only `ForwardDiff`
(a hard dependency of Pyrolysis), so this manual route works in any environment.

---

## 13.5 Adjoint-mode local sensitivity

For a **scalar loss** and **many parameters**, reverse mode (adjoint) computes the
full gradient at a cost roughly independent of `length(θ)`. The entry point is:

```julia
adjoint_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function,
    loss_fn::Function;                # solution -> scalar ℝ
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
    solve_kwargs...,
) -> (loss, ∂loss_∂θ)
```

- `loss_fn` reduces a solution to a single real number (a misfit, an objective).
- `sensealg` is a SciMLSensitivity reverse-mode algorithm. The default is
  `InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`. It is forwarded to
  `solve` as the first-class `sensealg` keyword (the only reason `solve` accepts a
  `sensealg` kwarg at all), so the ODE integrator augments the solve with adjoint
  dynamics. The gradient is taken with `Zygote.gradient`.

### Worked example: gradient of a TGA/cone misfit

A typical calibration loss is the sum of squared errors between the simulated and
measured surface-temperature (or MLR) histories on the saved time grid:

```julia
using Pyrolysis
using SciMLSensitivity, Zygote

# Measured surface temperature at the same output_times as the problem.
T_data = load_my_surface_temperature()   # Vector{Float64}, your data

loss_fn(sol) = sum(abs2, sol.surface_T .- T_data)

θ0 = [1.0e13, 1.88e5, 8.70e5]
base = p_inject(nothing, θ0)

loss, grad = adjoint_sensitivity(base, θ0, p_inject, loss_fn;
                                 radiation_model = SURFACE_ABSORPTION,
                                 abstol = 1e-9, reltol = 1e-7)
# grad[k] = ∂loss/∂θ_k  — feed straight into a gradient-based optimizer
```

`grad` is what an optimizer (e.g. `Optim.LBFGS`, `Optimisers.Adam`) consumes to
fit `(A, E, h, …)` to data. This is the inverse-analysis recipe in §13.7.

> **Try forward mode first for small problems.** Adjoint mode carries
> significant machinery (continuous adjoint ODE, VJP construction) and can be
> sensitive to `sensealg` choice on stiff Arrhenius transients. For `≲ 20`
> parameters, forward mode is usually simpler and just as fast in wall-clock
> terms. Reserve adjoint mode for genuinely high-dimensional fits.

### Choosing a `sensealg`

`sensealg` is forwarded straight to the underlying ODE solve. SciMLSensitivity
offers several reverse-mode algorithms; the default `InterpolatingAdjoint` stores
a continuous interpolation of the forward solution and is a good general choice.
If you hit instability or excessive memory use, consult the SciMLSensitivity
documentation for alternatives (e.g. `QuadratureAdjoint`, `GaussAdjoint`, or a
different `autojacvec`); pass them through the `sensealg` keyword unchanged.

---

## 13.6 Loading the extension

`forward_sensitivity` and `adjoint_sensitivity` are exported stubs. Without the
extension they raise a clear error telling you which package to load. The
extension `PyrolysisSciMLSensitivityExt` is triggered automatically when **both**
`SciMLSensitivity` and `Zygote` are present in your environment:

```julia
using Pyrolysis
using SciMLSensitivity     # both are required to trigger the extension
using Zygote
```

The extension is declared in `Project.toml` as
`PyrolysisSciMLSensitivityExt = ["SciMLSensitivity", "Zygote"]`, so both packages
must be installed and loaded. (Forward mode via `ForwardDiff`, §13.4, does not
need the extension.) Once loaded, the concrete methods attach to the exported
stubs and the calls above work.

---

## 13.7 Inverse analysis: fitting parameters to data

The forward/adjoint machinery exists to support **inverse analysis** — recovering
material parameters from TGA, DSC, MCC, or cone-calorimeter measurements. The
pattern is:

1. **Write `p_inject(base, θ)`** that rebuilds the problem with the candidate
   parameters at the right places (§13.2). Make sure `output_times` matches the
   times at which you have data so simulated and measured curves align.
2. **Write a scalar `loss_fn(sol)`** comparing the relevant simulated field
   (`sol.surface_T`, `sol.mass_loss_rate`, `sol.extras.total_mass`, …) to data.
3. **Get the gradient** with `adjoint_sensitivity` (many parameters) or assemble
   it from `forward_sensitivity` on a 1-vector output (few parameters).
4. **Drive an optimizer** with `(loss, grad)`.

A compact LBFGS loop using forward-mode gradients (no extension required):

```julia
using Pyrolysis, ForwardDiff, Optim

T_data = load_my_surface_temperature()
loss(θ) = sum(abs2, solve(p_inject(base, θ);
                          progress = false,
                          radiation_model = SURFACE_ABSORPTION).surface_T .- T_data)

θ0 = [1.0e13, 1.88e5, 8.70e5]
res = optimize(loss, θ0, LBFGS(); autodiff = :forward)   # ForwardDiff under the hood
θ_fit = Optim.minimizer(res)
```

For larger parameter vectors, supply an explicit adjoint gradient instead:

```julia
function g!(G, θ)
    _, grad = adjoint_sensitivity(base, θ, p_inject,
                                  sol -> sum(abs2, sol.surface_T .- T_data);
                                  radiation_model = SURFACE_ABSORPTION)
    G .= grad
end
res = optimize(loss, g!, θ0, LBFGS())
```

> **Scale your parameters.** `A ~ 10¹³`, `E ~ 10⁵`, and `h ~ 10⁶` differ by many
> orders of magnitude. Optimize in log-space for `A` (`θ = log10(A)`) or
> non-dimensionalize so the optimizer sees comparable scales; otherwise the line
> search struggles and gradients are numerically lopsided. Apply the inverse
> transform inside `p_inject` (e.g. `A = 10^θ[1]`).

The shipped `examples/df_sensitivity_comparison.jl` is **not** a programmatic
fitter — it parses a ThermaKin2Ds reference output
(`docs/DF_sensitivity_Cone_output.txt`) and an already-fitted Julia result CSV and
overlays them (MLR, surface/back temperature, mass, thickness, component MLR). Use
it as a template for *validating* a fit against a reference code, not for
generating one. The actual fitting workflows live under
`examples/inverse_analysis/`.

---

## 13.8 Global (Sobol) sensitivity via the sweep harness

Variance-based global sensitivity quantifies how much of an output's variance each
input (or group of inputs) is responsible for, over a *distribution* of inputs
rather than at a single point. Pyrolysis.jl has no dedicated Sobol function;
instead you build a sampling design, run it through `parameter_sweep`
(Chapter 12), and post-process. The shipped example
`examples/sensitivity_analysis/sensitivity_douglas_fir.jl` is a complete,
production-grade realization of this for a 16-parameter Douglas-fir cone model,
and `examples/sensitivity_analysis/sensitivity_sobol.jl` computes the variance
decomposition from its output.

The workflow has three stages:

**1. Sample the parameter space.** The example uses Owen-scrambled Sobol
quasi-random samples from `QuasiMonteCarlo.jl`, transformed to correlated
multivariate-Normal / Student-t groups whose covariances come from the
measurement campaign:

```julia
using QuasiMonteCarlo, Distributions, LinearAlgebra

# n_sobol unit-hypercube points, then map to each group's distribution.
sobol_unit = QuasiMonteCarlo.sample(
    n_samples, zeros(n_sobol), ones(n_sobol),
    SobolSample(R = OwenScramble(base = 2, pad = 32)),
)
# ... Φ⁻¹ to standard normal, Cholesky of the group covariance, t-scaling ...
# samples[:, j] is the j-th 16-parameter realization θ_j
```

**2. Solve every sample.** Each column of `samples` is a full parameter vector
`θ`; a `build_and_solve(θ)` function (the example's `p_inject`-equivalent)
rebuilds the Douglas-fir material and problem and solves it. Because every solve
is independent, this parallelizes cleanly with `Distributed.jl`:

```julia
using Distributed
@everywhere using Pyrolysis

# build_and_solve(θ) rebuilds the material/problem from θ and returns a solution.
# extract_qois(sol) reduces it to scalar QoIs + curves.
results = pmap(j -> extract_qois(build_and_solve(samples[:, j])), 1:n_samples)
```

This is the same separation as `parameter_sweep`: an immutable base, a `modify`/
`build_and_solve` that materializes a fresh problem per point, and an `extract`
that pulls out quantities of interest (peak MLR, time to ignition, total mass
loss, peak surface temperature, MLR/temperature curves, …). The example's solve
uses `use_ale = true`, `radiation_model = SURFACE_ABSORPTION`, and the default
`Structured(strategy = DirectSolve())` Jacobian.

**3. Decompose the variance.** `sensitivity_sobol.jl` fits a Polynomial Chaos
Expansion (PCE) surrogate to each scalar QoI and reads first-order and total
*group* Sobol indices analytically from the PCE coefficients, with bootstrap
confidence intervals. The Douglas-fir study groups the 16 parameters into 7
measurement sources (HFM/DSC calibration biases, three heat-capacity regressions,
heat-of-sorption parameters, and five heats of reaction), and reports how much of
each QoI's variance each *source* explains. A lighter-weight alternative used in
the same example is **Standardized Rank Regression Coefficients (SRRC)**, which
need only the sampled `(θ, QoI)` pairs and no surrogate.

> **You do not need PCE to get useful global sensitivity.** SRRC (rank regression
> of standardized inputs against the QoI) is a few lines on top of the sweep
> output and is robust for monotone responses. Use PCE when you want variance
> *decomposition* with interaction terms. Either way, the expensive part is the
> `N_samples` solves, which the harness parallelizes.

For the theory of variance-based sensitivity and its relation to the sweep
harness, see Technical Reference §17.5.

---

## 13.9 Pitfalls

Sensitivities are only meaningful where the map `θ -> output` is well-behaved.
The package is deliberately built to be AD-clean (smooth tanh temperature gates
and depletion limiter, midpoint-rule enthalpy, smooth flux limiters — Technical
Reference §17.6), but several configurations can still spoil a gradient:

- **Warm start.** Reusing a prior solution as the initial state
  (`warm_start = previous_solution`) makes the result depend on history, which is
  *not* differentiated by `p_inject`. For sensitivity studies, **start each solve
  from the problem's own initial conditions** (do not pass `warm_start`), so the
  differentiated chain is self-contained. This is why the sweep/sensitivity
  examples build a fresh problem per point and do not warm-start.

- **Non-smooth lateral-shrinkage laws (WithChi mode).** The cross-section law
  `A(t) = A_0 · law(χ̄)` is differentiated with ForwardDiff inside the residual
  (Technical Reference §10.6). A `law` with a kink (e.g. a hard `min`/`max` or a
  piecewise definition with a slope discontinuity) injects a non-smooth point into
  the sensitivity. Use smooth shrink/swell laws when you intend to differentiate
  through them.

- **Surface-depletion callback discontinuities.** Cell merging
  (`handle_depletion = true`, requires `use_ale = true`) resizes the state vector
  mid-solve at a discrete event — a genuine discontinuity in the trajectory and
  its length. Differentiating across a merge is not meaningful, and the output
  dimension can change between `θ` probes. **Disable depletion for sensitivity
  studies** (`handle_depletion = false`, the default); if you need long ALE runs,
  stop before the first merge or use a thickness-termination criterion instead.

- **Depletion + pluggable Jacobian are mutually exclusive.** `solve` rejects
  `handle_depletion = true` together with an explicit `jacobian = ...`
  (state-vector resizing is incompatible with a fixed pluggable backend); when
  depletion is enabled it falls back to the integrator's colored finite-difference
  Jacobian. This is a further reason to keep depletion off in sensitivity work,
  where you want the default `Structured` backend and a smooth path.

- **Transpiration BC is forbidden with the default convention.** Leave
  `use_transpiration_bc = false` (the default). The volumetric advection source
  `S_conv` already carries gas enthalpy; enabling the surface transpiration term
  double-counts it, and problem construction / `solve` rejects it outright. There
  is no sensitivity question that requires turning it on.

- **Tanh-gate leakage.** The smooth temperature gates leak a few percent of the
  reaction rate within roughly ±10 K of `T_min`/`T_max` (Technical Reference
  §6.2). Near those edges, `∂r/∂T_min` and `∂r/∂T_max` are smooth but small and
  not physically meaningful at the very edge of the window; interpret
  gate-boundary sensitivities with care.

- **Integrator tolerances set a noise floor.** Finite `abstol`/`reltol` introduce
  small, parameter-dependent integration error. For tight sensitivities (and for
  reproducible adjoints), tighten tolerances (e.g. `abstol = 1e-9`,
  `reltol = 1e-7`, as the Douglas-fir example does) so the sensitivity signal sits
  well above the tolerance noise.

---

## 13.10 Quick reference

| Task | Call |
|---|---|
| Forward Jacobian of an output | `forward_sensitivity(base, θ, p_inject; output_fn)` |
| Forward Jacobian, no extension | `ForwardDiff.jacobian(θ -> out(solve(p_inject(base, θ))), θ)` |
| Adjoint gradient of a scalar loss | `adjoint_sensitivity(base, θ, p_inject, loss_fn; sensealg)` |
| Load the extension | `using SciMLSensitivity, Zygote` |
| Global/Sobol | sample → `pmap`/`parameter_sweep` → PCE or SRRC (see examples) |

**Defaults to remember.** `forward_sensitivity`: `output_fn = identity`.
`adjoint_sensitivity`: `sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`.
For sensitivity-safe solves: keep `warm_start = nothing`,
`handle_depletion = false`, `use_transpiration_bc = false`, the default
`jacobian = nothing` (→ `Structured(strategy = DirectSolve())`), and tighten
`abstol`/`reltol`.

For parameter sweeps (discrete response maps) see Chapter 12; for the solution
fields and accessors used in `output_fn`/`loss_fn` see Chapter 9; for the solver
options forwarded through `solve_kwargs` see Chapter 8; for the underlying theory
see Technical Reference §17.


---

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
| `06_cone_calorimeter.jl` | PMMA cone, single-step kinetics, in-depth radiation | ALE, depletion | 14.1 |
| `06_cone_calorimeter_surface_absorption.jl` | PMMA cone, surface-only radiation | ALE, SURFACE_ABSORPTION | 14.1 |
| `07_douglas_fir_cone.jl` | Wood, multi-step kinetics, Darcy–Fick, moisture | ALE, SURFACE_ABSORPTION, DarcyFick | 14.2 |
| `convergence_douglas_fir.jl` | Grid-convergence study | as 07, swept over `n_cells` | 14.3 |
| `08_pressure_treated_wood_cone.jl`, `ptw_comparison.jl` | Treated wood + moisture, Darcy–Fick | ALE, SURFACE_ABSORPTION, DarcyFick | 14.4 |
| `08_darcy_fick_demo.jl` | Pressure-driven gas transport, permeability | ALE, DarcyFick | 14.5 |
| `moisture_hos_comparison.jl` | State-dependent heat of sorption | ALE, SURFACE_ABSORPTION, DarcyFick | 14.6 |
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
    handle_depletion = true,     # merge thin surface cells (production depletion path)
    depletion_threshold = 0.05,  # merge a cell at 5% of its initial thickness
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
- **`handle_depletion = true`** turns on the surface-cell merge callback. It
  requires `use_ale = true`. Because it resizes the state vector mid-solve, it is
  **incompatible with a pluggable `jacobian` backend** (the solver errors if you
  pass both) and falls back to the integrator's colored finite-difference
  Jacobian (Technical Reference §14, §15).
- **`use_transpiration_bc`** is left at its default `false`. Enabling it is
  rejected: the default `:standard` energy form already carries gas enthalpy via
  the `S_conv` advection term, so a surface-transpiration term would
  double-count it (Technical Reference §7; see also the Troubleshooting chapter).

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
script overlays these against bundled ThermaKin reference traces and prints a
side-by-side comparison at t ≈ 600 s.

If you have `Plots` loaded, the package's plot recipe gives quick views:

```julia
using Plots
plot(solution; kind = :surface_T)   # surface temperature vs time
plot(solution; kind = :mlr)         # mass-loss rate
plot(solution; kind = :heatmap)     # T(z, t) field
```

> Caveat: this example's parameter set was calibrated against a solver that
> applies a surface transpiration term in addition to the volumetric advection
> source. Under the `S_conv`-only convention used here, surface temperature and
> MLR shift slightly, so the bundled ThermaKin comparison is indicative, not a
> regression target; treat the printed peak-surface-T and mass-loss numbers as a
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

The material declares 13 components — six solids (`DF`, `DF_int1`…`DF_int4`,
`Char`), one liquid (`MOISTURE`), and six gases (`DF_gas1`…`DF_gas5`,
`WATER_VAPOR`) — and six reactions. A few patterns worth lifting:

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

The full 13-component, six-reaction definition is in
`examples/07_douglas_fir_cone.jl` (lines 216–445).

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

The default Jacobian is `Structured(strategy = DirectSolve())`, which assembles
the exact analytical Jacobian including the dense ALE mesh-velocity coupling. For
large ALE problems you can drop that coupling (a Schur-style approximation that
still converges every term through Newton) for a substantial speedup:

```julia
solution = solve(problem;
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
    integrator = KenCarp4(linsolve = SparspakFactorization()),
)
```

Legacy benchmarks measured ~7.8× speedup at 240 cells with solutions identical
to within solver tolerance (Technical Reference §15; see the Solver/Jacobian
configuration chapter for choosing strategies).

### What to look at

Because radiation is active, the solver records a full boundary-flux
decomposition in `extras`:

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
NC, NR = 13, 6
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
per boundary in `mesh.boundary_states` and can be queried with
`get_surface_temperature(mesh)` / `get_back_face_temperature(mesh)` (exported
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
    dt_initial = 0.01,
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

> Caveat: warm starts, non-smooth lateral-shrinkage laws, and the
> depletion-merge callback introduce discontinuities that can confuse
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
- **Use the boundary-flux extras for energy accounting.** With a radiation model
  active, `surface_heat_flux_*` and `surface_emissivity`/`surface_conductivity`
  in `extras` give a complete surface energy balance without re-deriving it.

### Caveats to keep in mind

- **`P1_QUASI_STEADY` radiation is experimental** and quarantined in
  `Pyrolysis.Experimental`; the core `solve` rejects it as a `radiation_model`.
  The shipped examples use `NO_RADIATION` (flux via BC only), `SURFACE_ABSORPTION`,
  or `BEER_LAMBERT`.
- **Surface transpiration BC is off by default and incompatible with the default
  `:standard` energy form** — enabling `use_transpiration_bc` is rejected because
  `S_conv` already carries gas enthalpy (double counting). Technical Reference §7.
- **Depletion handling and a pluggable Jacobian are mutually exclusive.**
  `handle_depletion = true` resizes the state vector mid-solve, so the solver
  errors if you also pass a `jacobian = JacobianSpec(...)`, and it falls back to
  the integrator's colored finite-difference Jacobian. Technical Reference §14–§15.
- **Several cone examples carry parameter sets calibrated under a
  transpiration-plus-`S_conv` convention** and may need re-calibration under the
  `S_conv`-only energy balance; treat their bundled experimental/ThermaKin
  overlays as indicative comparisons rather than tight regression targets
  (Technical Reference §18).

For deeper troubleshooting (stiff solves, ledger imbalance, ALE mesh tangling,
radiation surprises), see the Troubleshooting and Performance chapter; for the
full keyword reference, see the Running the Solver and API Reference chapters.


---

# 15. Troubleshooting and Performance

This chapter is a field guide for when a Pyrolysis.jl simulation fails, runs
slowly, or produces results you do not trust. It is organised by symptom: each
section names a class of failure, explains the underlying cause, and gives the
concrete `solve` options, diagnostic calls, or model changes that fix it. The
emphasis throughout is on *root-causing* — most of the tools below let you
localise a problem to a single subsystem (kinetics, conduction, gas transport,
radiation, ALE, or the linear solver) before you start changing inputs. For the
physics and numerics behind the knobs mentioned here, follow the "see Technical
Reference §N" pointers; for the full option list see the *Running the Solver*
chapter (UG-8) and the *API and Parameter Reference* (UG-16).

Coordinate and sign conventions follow the unified nomenclature: `z = 0` is the
substrate (bottom, boundary id 2, tag `:bottom`), `z = L` is the exposed surface
(top, id 1, tag `:top`); positive flux is transport in `+z`; heat of reaction
`h > 0` is endothermic with `Q_rxn = −Σ h_r r_r`.

---

## 15.1 First steps: read the return code and run the diagnostics

Before changing anything, look at two things: the solver's return code and the
conservation ledgers.

### 15.1.1 The return code

`solve` stores the integrator's exit status on the solution:

```julia
solution = solve(problem; progress = false)
solution.retcode      # Symbol mirrored from raw_solution.retcode
```

The common values and what they mean:

| `retcode` | Meaning | Typical cause |
|---|---|---|
| `:Success` | Reached `tspan[2]` | — |
| `:Terminated` | A callback stopped the run | Sample fully depleted (ALE thickness termination) or a user callback called `terminate!` |
| `:MaxIters` | Hit `maxiters` before `t_end` | Stiff transient; step size collapsed; under-resolved front |
| `:DtLessThanMin` | Step fell below `dtmin` | Discontinuity, NaN in the RHS, or genuine stiffness wall |
| `:Unstable` / `:DtNaN` | NaN/Inf in the state | Unphysical state (negative concentration, `ρc_p → 0`), bad property table extrapolation |

A `:Terminated` run is *not* an error in ALE mode — it usually means the sample
burned away. Check `solution.extras.thickness_history[end]` against
`min_thickness` (default `1e-6` m).

### 15.1.2 Turn on the conservation ledgers

The single most useful debugging switch is `diagnostics = true`. It installs a
callback that snapshots the mass and energy ledgers every `diagnostics_stride`
accepted steps (default 10) and stores typed rows in
`solution.extras.diagnostics_log`:

```julia
solution = solve(problem;
    diagnostics        = true,
    diagnostics_stride = 5,     # snapshot every 5 accepted steps
    progress           = false,
)

for row in solution.extras.diagnostics_log
    println("t=$(row.t)  E_matrix=$(row.E_matrix)  ",
            "r_disc(rel)=$(row.discrete_relative)  ",
            "r_phys(rel)=$(row.physical_relative)  ",
            "ΔE_gas=$(row.gas_storage_gap)  ",
            "m_err=$(row.mass_max_relative_error)")
end
```

Each `DiagnosticsRow` carries `t`, `E_matrix`, `E_total`, `discrete_residual`,
`discrete_relative`, `physical_residual`, `physical_relative`,
`gas_storage_gap`, and `mass_max_relative_error`. How to read them (see
Technical Reference §16):

- **`discrete_relative`** audits the *bookkeeping*. It is assembled from exactly
  the terms the RHS computes, so it closes to integrator tolerance (~`1e-6`
  relative) **by construction**. If it is large, you have a bug in flux/source
  assembly or a callback — not a physics problem. (Caveat: with strongly
  time-varying composition it cannot close exactly, because the ledger omits the
  `∫ T ∂(ρc_p)/∂t` and ALE volume-change terms; expect a small floor.)
- **`physical_relative`** audits the *quasi-steady-gas approximation*. For
  canonical 1-atm pyrolysis it should stay **< 1 %** of initial energy. A larger
  value means the regime stresses the gas-storage assumption.
- **`gas_storage_gap`** (`ΔE_gas = E_total − E_matrix`) is the sensible heat
  carried by gas that the matrix-only `ρc_p^eff` does not store. Canonically
  < 0.2 % of `E_matrix`.
- **`mass_max_relative_error`** is the worst per-component mass-balance error;
  it should be `< 1e-10`.

`diagnostics` is **snapshot-based**, not continuous — it gives you a coarse time
history adequate for spotting where a ledger drifts, not a dense trace. The
overhead is one extra energy/mass sweep every `diagnostics_stride` steps.

---

## 15.2 Stiff or slow solves

Arrhenius kinetics make the system stiff during ignition and the early reaction
transient; conduction adds a parabolic time scale. If a solve is slow, hits
`:MaxIters`, or grinds to `:DtLessThanMin`, work through the following.

### 15.2.1 Confirm the default integrator is appropriate

The default integrator is `KenCarp4` (an SDIRK method) with a sparse direct
linear solve; its `autodiff = AutoFiniteDiff()` setting only takes effect on
the no-Jacobian (depletion) path, since `solve` normally supplies the
structured analytical Jacobian (§15.2.3). When both `integrator` and
`linear_solver` are left at `nothing`, `solve` builds:

```julia
# Internal default:
# KenCarp4(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
```

`KenCarp4` is a good default for moderately stiff pyrolysis. For *very* stiff
problems (high pre-exponential `A_i`, sharp ignition) a BDF method often takes
larger steps. Integrators live in `OrdinaryDiffEq`, which is **not** re-exported
by `using Pyrolysis`; load it (or the relevant sublibrary) yourself and pass a
constructor or a fully-built algorithm to the `integrator` kwarg:

```julia
using OrdinaryDiffEq            # for FBDF, Rodas5P, etc.

# (a) pass a constructor — solve() injects autodiff + linsolve for you:
solution = solve(problem; integrator = FBDF)

# (b) pass a constructor and a linear solver together:
solution = solve(problem;
    integrator    = FBDF,
    linear_solver = SparspakFactorization(),   # re-exported by Pyrolysis.Solver
)

# (c) pass a fully-constructed algorithm — then DO NOT also pass linear_solver:
solution = solve(problem; integrator = FBDF(linsolve = SparspakFactorization()))
```

The pairing rule is strict: `linear_solver` may only be combined with
`integrator` when `integrator` is a **constructor** (a `Type`/`Function`). If
you hand `solve` an already-constructed algorithm *and* a `linear_solver`, it
raises an error telling you to fold the `linsolve` into the algorithm yourself.

### 15.2.2 Adjust tolerances

The error tolerances are `abstol = 1e-8` and `reltol = 1e-6` by default. They
trade accuracy for speed directly:

```julia
# Looser — faster, for exploratory runs:
solution = solve(problem; abstol = 1e-6, reltol = 1e-4)

# Tighter — for convergence/validation studies:
solution = solve(problem; abstol = 1e-10, reltol = 1e-8)
```

If you see `:DtLessThanMin`, first try *loosening* tolerances; a too-tight
`reltol` against a stiff transient drives the step size into the floor. The step
floor itself is `dtmin = 1e-12` s and the ceiling `dtmax = Inf`; raising
`dtmin` is rarely the right fix — it usually just converts a slow run into a
failed one. `maxiters` defaults to `1_000_000`; a `:MaxIters` exit at that count
almost always means the step size collapsed, so chase the stiffness, not the
iteration cap.

### 15.2.3 Use the structured Jacobian (it is already the default)

When `jacobian = nothing` (the default) and `handle_depletion = false`, `solve`
builds a `Structured(strategy = DirectSolve())` analytical Jacobian and hands a
sparse prototype to the integrator. This is much faster than the integrator's
own colored finite-difference Jacobian because it is exact and sparsity-aware
(see Technical Reference §15). You normally do not need to touch it. The cases
where you *do* are covered in §15.3 (verification), §15.4 (depletion), and
§15.7 (large ALE).

### 15.2.4 Reduce mesh size or use warm starts

Cost scales with the number of cells (the state dimension is `n·(1+NC)` plus the
ALE `z` block and the WithChi `χ` block when present). If a 400-cell run is slow,
check whether 100–150 cells already give a converged surface temperature and MLR
(run a quick convergence check — see Technical Reference §18.5 and the
`convergence_douglas_fir` example). For a *sequence* of related solves
(continuation, sweeps), reuse the prior final state with `warm_start`:

```julia
sol1 = solve(problem1; progress = false)
sol2 = solve(problem2; warm_start = sol1, progress = false)
```

The state-vector length, component count `NC`, ALE mode, and χ mode of the two
problems must match (enforced by a length check). `tspan[1]` still sets the
integrator's start time, so the prior final state is treated as fresh initial
data.

---

## 15.3 Jacobian mismatch warnings

If the integrator behaves erratically only when the structured Jacobian is
active — many rejected steps, sudden step-size collapse that goes away when you
fall back to the dense reference — suspect a Jacobian that disagrees with the
residual.

### 15.3.1 Compare against the dense reference at sample times

The `verify_jacobian` kwarg installs a callback that, at the listed sample
times, compares the active backend against a `DenseAD_ForwardDiff` reference
(full-residual ForwardDiff — correct by construction, slow). It **warns**, it
does not throw, so a single transient outlier will not abort a good run:

```julia
solution = solve(problem;
    verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 1.0, 10.0]),
    progress        = false,
)
```

The NamedTuple keys are restricted to `:rtol`, `:atol`, and `:sample_times`
(unknown keys raise an error). A warning reports the backend type, the time, the
max absolute and relative deviations, and the worst index. The callback fires
when the integrator crosses any sample time within half the current step.

### 15.3.2 Verify at a single state point

For a focused check (e.g., debugging an operator), call `verify_jacobian_at`
directly. You need the `ProblemDef`/`Workspace` pair, which you build from the
problem:

```julia
def = to_problem_def(problem; use_ale = false, radiation_model = NO_RADIATION)
ws  = build_workspace(problem, def)
u   = create_state_vector(ws.layout, ws.mesh)

spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
report = verify_jacobian_at(spec, def, ws, u, 0.0;
    reference = DenseAD_ForwardDiff(),   # default reference backend
    rtol = 1e-6, atol = 1e-10,
)

report.passed          # Bool
report.max_abs_dev     # max absolute deviation
report.max_rel_dev     # max relative deviation
report.worst_index     # (row, col) of the worst entry
report.mode            # :strict (exact backends) or :masked (approximate)
```

`JacobianSpec` also accepts a `PyrolysisProblem` directly
(`JacobianSpec(backend, problem)`) if you do not already have a `def`/`ws` pair.

### 15.3.3 Support-pattern (sparsity) issues

The comparison `mode` matters:

- **`:strict`** — for exact backends (`DirectSolve`, `StructuredSolve`). Both
  in-support disagreements *and* reference entries the backend's sparsity
  *omits* count as failures. A strict failure on an outside-support entry means
  the sparsity pattern is too small — the backend is structurally missing a
  coupling.
- **`:masked`** — for approximate backends (`ApproxSolve`, where
  `is_approximate(backend)` is `true`). Only in-support disagreements are
  asserted; omitted entries are informational. A masked pass with large
  outside-support reference magnitude is *expected* for `ApproxSolve` (it
  deliberately drops the structured ALE couplings); a masked run that struggles
  to converge is the signal to switch back to an exact strategy (§15.7).

If `verify_jacobian_at` reports a strict, outside-support failure on a
production backend, that is a package bug — report it with the `worst_index` and
the problem configuration. The structured operators are themselves checked
against the dense reference in the test suite (Technical Reference §15.13,
§18.3).

---

## 15.4 Depletion and the pluggable-Jacobian incompatibility

Surface-cell depletion (`handle_depletion = true`, ALE only) resizes the state
vector mid-solve via `resize!(integrator, new_length)` when a surface cell
thins below `depletion_threshold` of its initial thickness (default 5 %) and
merges it into its interior neighbour (see Technical Reference §14.3).

The pluggable `JacobianSpec` cache is **immutable** and built once at setup, so
it cannot survive a resize. The two are therefore mutually exclusive:

```julia
# ERROR: solve() rejects this combination up front.
solve(problem;
    use_ale          = true,
    handle_depletion = true,
    jacobian         = JacobianSpec(Structured(strategy = DirectSolve()), problem),
)
# → "the `jacobian` kwarg is incompatible with `handle_depletion = true`"
```

The supported way to run depletion is to **let `solve` pick the Jacobian**: when
`handle_depletion = true`, the auto-default returns `nothing` and the integrator
falls back to its own colored finite-difference Jacobian, which *does* survive a
`resize!`:

```julia
solution = solve(problem;
    use_ale             = true,
    handle_depletion    = true,
    depletion_threshold = 0.05,   # merge below 5% of initial thickness
    min_cells           = 2,      # stop merging at 2 cells
    min_thickness       = 1e-6,   # terminate when total thickness < this
    progress            = false,
)
```

The merge is energy-conserving (it does not blow up when a near-empty cell with
`ρc_p → 0` is absorbed; Technical Reference §14.3). Consequence: with depletion
you give up the fast structured Jacobian. If a depletion run is slow, see
§15.4.1 and the linear-solver choice in §15.11.

### 15.4.1 The MKL `ipiv` resize bug

With `handle_depletion = true` and a large state vector (> ~500 elements),
passing a fully built integrator whose `linsolve` is left unset lets
LinearSolve.jl's auto-selection pick an MKL LU factorization whose pivot array
(`ipiv`) is not resized correctly on `resize!`, producing a
`DimensionMismatch: ipiv has length X, but needs to be Y`. `solve` warns when
it detects this risky combination.

The fix is to let `solve` resolve the linear solver: on the depletion path it
defaults to Julia's native `LUFactorization()`, which handles `resize!`
correctly (a sparse factorization cannot run here anyway — there is no sparse
Jacobian prototype — so an explicit `SparspakFactorization`/`KLUFactorization`
is rebuilt to `LUFactorization` with a warning):

```julia
solution = solve(problem;
    use_ale          = true,
    handle_depletion = true,
    progress         = false,
)
```

For the optional `KLU`, `Pardiso`, and AMG extensions and how to load them, see
UG-11. The `klu_factorization`, `pardiso_factorization`, and `amg_preconditioner`
stubs error with a `MethodError` until you `import` the matching package — that
error *is* the prompt to load the dependency.

---

## 15.5 Negative or unphysical states

Symptoms: NaN in the temperature or concentration field, `:Unstable` retcode,
concentrations going slightly negative, or oscillations behind a sharp front.

### 15.5.1 Concentrations going negative

The kinetics use a **depletion limiter** — a smooth `tanh(ξ/threshold)` factor
on each reaction rate (default threshold 1.0 kg/m³) — so reactions switch off
as a reactant is exhausted, for every reaction order n ≥ 0; this is what keeps
the system differentiable and stiffness-bounded
(see Technical Reference §6.3). If you still see meaningful negative
concentrations, the usual causes are:

- **Too-coarse a mesh at the reaction front.** Refine near the surface with a
  stretched mesh (§15.5.4) or increase `n_cells`.
- **First-order upwind overshoot in Darcy–Fick advection.** The advective gas
  flux uses first-order upwinding (stable but diffusive; no TVD limiter on the
  Darcy advection — see Technical Reference §9.3). Small negative excursions of a
  near-zero gas concentration are expected and are damped by reaction/transport
  coupling; large ones indicate the pressure field or permeability is driving
  too strong an advective flux (check §15.9).
- **A reactant initialised at exactly zero that a reaction consumes.** Make sure
  `ξ_initial` seeds every reactant the scheme needs.

The smooth temperature gates `½(1+tanh(T−T_min))·½(1+tanh(T_max−T))` leak by
roughly ±10 K around `T_min`/`T_max` (Technical Reference §6.2); if a reaction
appears to fire slightly early/late, that leakage — not a bug — is the reason.

### 15.5.2 `ρc_p → 0` blow-up (ALE)

As a cell decomposes to gas, its matrix-only `ρc_p^eff` (gas is excluded by
Formulation A) shrinks, and `dT/dt = RHS/ρc_p` can blow up. The intended
remedies are ALE surface depletion (§15.4) to merge the spent surface cell, and
the thickness-termination callback that stops the run when total thickness drops
below `min_thickness`. If you run ALE *without* depletion and hit NaNs near full
burn-out, enable `handle_depletion = true`.

### 15.5.3 Property-table extrapolation

`TableProperty` interpolates linearly and extrapolates beyond its range. A table
that stops at 600 K queried at 900 K will extrapolate — possibly to a negative
`c_p` or `k`, which then produces NaN. Make sure every `TableProperty` covers
the full temperature span the solve will reach (the Newton boundary solve clamps
surface temperature to `[200, 3000]` K — see Technical Reference §12.2 — so size
tables accordingly).

### 15.5.4 Refine near the surface

Steep thermal and reaction fronts sit near the exposed surface. A geometric
stretch concentrates cells there at fixed total count:

```julia
# Fine cells near z = L (top/surface), coarsening toward z = 0 (substrate):
mesh = create_stretched_mesh(0.01, 100, Val(2), 1.05; stretch_from = :top)
#                            L      n   NC     stretch  direction
```

`create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from = :top, ...)`
uses geometric grading; `stretch > 1` grades the cell sizes across the domain.
Compare against the uniform `create_uniform_mesh(L, n_cells, Val(NC))` to see
whether near-surface refinement removes the overshoot. For dynamic refinement
options (error indicators, `relocate_nodes!`), see UG-10.

---

## 15.6 Energy or mass ledger not closing

Use the diagnostics rows (§15.1.2) to localise. The two energy residuals fail
for *different* reasons.

### 15.6.1 `discrete_relative` large

This audits assembly and should close to integrator tolerance. If it does not:

- **A callback is mutating state outside the residual.** The mesh is read-only
  during assembly and is synced only by the step-update callback; a custom
  callback that edits `u` or geometry without informing the integrator
  (`u_modified!`) will break the ledger.
- **AMR / depletion masking.** All domain-total sums (`compute_domain_mass`,
  `compute_sensible_energy`) iterate the **active-cell mask**, not
  `eachindex`. If you wrote custom post-processing that loops over all cells
  (including merged-away ones), your totals — not the solver's — are wrong.
  Always respect `mesh.active_cell_mask` (Technical Reference §16.1).

### 15.6.2 `physical_relative` large (> 1 %)

This audits the quasi-steady-gas approximation. A large value means the regime
genuinely stresses the assumption (gas residence time approaching the thermal
time scale — unusual at 1 atm). Reconsider the setup rather than the numerics;
`gas_storage_gap` will corroborate (it will be a non-trivial fraction of
`E_matrix`).

### 15.6.3 The S_conv vs transpiration double-count

This is the most common "my energy balance is off by a constant rate" trap. The
default energy form (`energy_form = :standard`) carries the gas sensible
enthalpy through the **always-on** volumetric convection source
`S_conv = −Σ_g c_{p,g} J̄_g ∂T/∂z` (ThermaKin Eq. 17 convention; Technical
Reference §7.4). The surface **transpiration** enthalpy term would carry the
*same* enthalpy at the boundary — enabling both double-counts by exactly
`∫ S_conv`. The package forbids the combination:

```julia
# ERROR: solve() rejects this.
solve(problem; use_transpiration_bc = true)
# → "use_transpiration_bc=true double-counts gas convective enthalpy
#    (S_conv is always active ...). Use use_transpiration_bc=false"
```

So `use_transpiration_bc` is `false` by default and there is no supported way to
turn it on through `solve`. Mass transport is unaffected either way — gas still
leaves the domain through the mass BCs; only the *enthalpy bookkeeping* is at
issue. (Parameter sets calibrated against a solver that applies a transpiration term
on top of a volumetric advection source embed the double-count and will predict
higher `T_s` and MLR here; re-validate imported calibrations.)

### 15.6.4 `:standard` vs `:with_generation_sink`

The optional gas-generation sink `S_gen = −Σ_g Δh_g (∇·J_g)` (Technical Reference
§7.5, §16.4) can be added:

```julia
solution = solve(problem; energy_form = :with_generation_sink)
```

`energy_form` accepts only `:standard` (default) or `:with_generation_sink`
(anything else errors). **Datum warning (see §7.7):** the sink is consistent only
with `T_ref`-datum (formation-style) reaction heats. With DSC-style heats quoted at
the reaction temperature — the usual meaning of calibrated heats — the gas sensible
enthalpy is already inside `h`, and the sink subtracts it a *second* time
(≈0.7–1.2 MJ per kg of volatiles at 900 K); `solve` emits a one-time warning when it
is selected. Even with the right datum it is only a *partial* step toward
conservative form (gas heat storage and condensed-product sensible generation remain
outside the assembly). The midpoint-rule `Δh_g` makes `∫(S_conv + S_gen)` telescope
exactly against the surface outflow, which is what the ledger audit checks. For most
cone/TGA work `:standard` is correct and matches FDS/Gpyro/ThermaKin.

---

## 15.7 ALE issues: tangling, GCL, vanishing cross-section

ALE (`use_ale = true`) lets the mesh follow the receding/expanding surface
(Technical Reference §11). Three things go wrong with moving meshes.

### 15.7.1 Mesh tangling and quality

If a swelling/intumescent material moves nodes past each other, cells invert
(`V ≤ 0`) and the geometry recompute produces NaNs. For diagnosis, the
mesh-motion utilities expose quality checks. These operate on a `Mesh1D` and
live in the `Pyrolysis.Internal` umbrella (they are not top-level exports):

```julia
using Pyrolysis.Internal          # mesh-motion helpers are not top-level exports

q = check_mesh_quality(mesh)       # NamedTuple
q.valid                            # all cells V > 0 ?
q.min_volume
q.max_aspect_ratio
q.tangled_cells                    # indices with V ≤ 0

mesh_is_valid(mesh)                # quick boolean
```

`move_mesh_safe!` rejects a motion (and rolls back) if it would tangle or exceed
an aspect-ratio bound (default 10). If you are driving the mesh yourself, prefer
it over the raw `move_mesh!`. For production solves the cumulative mesh velocity
is computed for you; tangling there points to an unphysical swelling factor `γ`
or shrinkage law producing a too-large per-step displacement — refine the time
step (lower `dtmax`) or reduce the swelling magnitude.

### 15.7.2 GCL violations

The discrete Geometric Conservation Law,
`V_new − V_old = Δt (w_R − w_L) A`, must hold to machine precision or the mesh
motion injects "phantom" mass (Technical Reference §11.5). Per-step GCL error
should be `O(1e-12)`–`O(1e-15)`; a value above `1e-10` signals geometric
degeneracy. To check explicitly:

```julia
using Pyrolysis.Internal

err = gcl_error(mesh, w, dt)        # predictive relative error (max over cells)
check_gcl(mesh, w, dt; tol = 1e-12) # Bool

# For a time-varying cross-section (WithChi mode) pass A_old, A_new:
gcl_residual(mesh, V_old, w, dt; A_old = A0, A_new = A1)
```

A `GCLTracker` accumulates the error over a run if you wire it into a custom
driver. In normal use, a GCL violation is not something you tune — it is a
symptom of a degenerate mesh (a cell driven to near-zero width); address the
underlying cause (§15.7.3).

### 15.7.3 Vanishing cross-section and thickness

In WithChi mode the cross-section evolves as `A(t) = A_0·law(χ̄)` (Technical
Reference §10.5). A shrinkage law `1 − k χ̄` with too large `k` can drive
`A(t) → 0`, at which point volumes vanish and the residual is singular. The
guards keep `law(χ̄) > ε`, but the right fix is a physical shrinkage coefficient.
Likewise, when the *axial* thickness collapses, the thickness-termination
callback fires and the run ends with `:Terminated` — check
`solution.extras.thickness_history` and `min_thickness`.

Note that MLR is reported two ways for variable cross-section:
`solution.mass_loss_rate` is per *current* area `A(t)`, while
`solution.extras.MLR_gauge` is normalised by the initial area `A_0` (the
cone-gauge convention). They coincide in non-ALE / constant-area runs. If your
MLR looks off by a factor that tracks `A(t)/A_0`, you are reading the wrong one.

### 15.7.4 Choosing a linear-solve strategy for large ALE

For large ALE problems the structured ALE couplings (cumulative mesh velocity,
χ back-reaction) dominate Jacobian assembly/solve cost. The strategies trade
exactness for speed (Technical Reference §15.10):

| Strategy | Exactness | When to use |
|---|---|---|
| `DirectSolve()` (default) | Exact; densify + factor | Small/medium problems; verification anchor |
| `StructuredSolve()` | Exact; sparse LU on the local block + Richardson (2–5 iters) | Typical ALE problems |
| `ApproxSolve()` | Approximate; drops structured couplings (`is_approximate == true`) | Large ALE where Newton tolerance can absorb a Schur-style approximation |
| `MatrixFreeSolve()` | GMRES on the operator action | Research/experimental |

Select via a `JacobianSpec` (remember: **not** with `handle_depletion = true`):

```julia
solution = solve(problem;
    use_ale  = true,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
)
```

If `ApproxSolve` increases Newton rejections on Arrhenius transients, fall back
to `StructuredSolve` (exact, still O(N) per coupling apply) or the default
`DirectSolve`. Use `verify_jacobian` (§15.3) with the masked mode to confirm an
approximate backend still tracks the reference inside its support.

---

## 15.8 Radiation surprises

The four radiation models are mutually exclusive and selected by the
`radiation_model` kwarg (default `NO_RADIATION`). Three of them are reachable
from `solve`; the fourth is quarantined.

### 15.8.1 Model selection and mutual exclusivity

```julia
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)  # opaque, O(1)
solution = solve(problem; radiation_model = BEER_LAMBERT)        # in-depth absorption
# radiation_model = NO_RADIATION  is the default
```

| Model | Behaviour | Cost / Jacobian |
|---|---|---|
| `NO_RADIATION` | No volumetric source | O(1), tridiagonal |
| `SURFACE_ABSORPTION` | Absorbed flux `α·Φ` applied at the surface | O(1), tridiagonal |
| `BEER_LAMBERT` | Exponential in-depth attenuation `I(z)` | O(n), upper-triangular (non-local) Jacobian |
| `P1_QUASI_STEADY` | P1 diffusion with re-radiation | **Experimental — rejected by `solve`** |

`NO_RADIATION`, `SURFACE_ABSORPTION`, and `BEER_LAMBERT` are exported by
`using Pyrolysis`. `P1_QUASI_STEADY` has moved to
`Pyrolysis.Experimental.RadiationP1` and is **not** part of the core unified
residual; passing it to `solve` errors with a message pointing you to the
Experimental module. Do not expect to use it through the standard `solve` path
(Technical Reference §8.12).

### 15.8.2 "I selected a radiation model but nothing absorbs"

`solve` warns about the two common mismatches:

- **`BEER_LAMBERT` with no absorption coefficients.** Beer–Lambert needs a
  non-zero mass absorption `α_j` on the absorbing components (the volumetric
  `α_eff = Σ_j ξ_j α_j`; Technical Reference §8.3). If every `α_j` is zero, the
  optical thickness `τ_i = α_i Δz_i` is ~0, `1 − exp(−τ) ≈ 0`, and effectively
  no radiation is deposited. Set `α` on the solid (e.g.
  `SolidComponent(:char, …, α = 5.0)`).
- **`SURFACE_ABSORPTION` with non-zero `α`.** The coefficients are *ignored* —
  the absorbed flux is applied at the surface. The warning tells you the in-depth profile
  you specified is not being used; switch to `BEER_LAMBERT` if you want depth
  penetration.

### 15.8.3 Optical-thickness regime

Choose the model by optical thickness `τ = α L` (Technical Reference §8.9):

| `τ = α·L` | Regime | Recommended |
|---|---|---|
| > 10 | Optically thick | `SURFACE_ABSORPTION` (or P1 if re-radiation matters) |
| 1–10 | Intermediate | `BEER_LAMBERT` |
| 0.1–1 | Optically thin | `BEER_LAMBERT` |
| < 0.1 | Transparent | `NO_RADIATION` |

Using `SURFACE_ABSORPTION` on a semi-transparent sample concentrates the
absorbed flux at the surface and overheats it; using `BEER_LAMBERT` on an opaque char
just adds non-local Jacobian cost for no benefit.

### 15.8.4 Time-varying incident flux is always tracked

For a **time-varying** incident flux (e.g. a cone ramp
`RadiativeFluxBC(flux = t -> 25e3 * (1 - exp(-t/10)))`), the flux is re-sampled
at every RHS evaluation whenever `radiation_model` is `SURFACE_ABSORPTION` or
`BEER_LAMBERT` — it is never snapshotted at `t = 0`, and there is no constancy
probe. A flux that is flat early and varies later (a delayed step at `t = 5`,
say) is tracked correctly too. The cost is one call of your flux function per
RHS evaluation, which is negligible; no special arrangement is needed.

### 15.8.5 Marshak ε coefficient (P1 only)

If you do reach into `Pyrolysis.Experimental.RadiationP1`, note the Marshak
boundary relation is `q = c·(G_ext + 4·I_ext − G)` with `c = ε/(2(2−ε))`
(Modest Eq. 16.48): the exchange coefficient is **not** `0.5ε`, and the external
flux carries the half-range weight `4c = 2ε/(2−ε)` — the same weight as the
ambient term — **not** 1. Either simplification silently under-injects radiation
for `ε < 1`. This is internal to the experimental P1 path; it does not affect
the core models reachable from `solve`.

---

## 15.9 Gas-transport and pressure surprises (Darcy mode)

Darcy–Fick transport activates when the material carries permeability (any
solid/liquid component with `κ > 0`); see UG-7 and Technical Reference §9.

- **No pressure-driven flow.** The face Darcy velocity is
  `v_g = −mob_face·∇P` with `mob_face` the distance-weighted harmonic mean of
  the two cells' mobilities `κ/μ`; the face is clamped to zero when either
  side's mobility falls below `1e-30` (an impermeable neighbor blocks the
  face) or a cell's viscosity below `1e-20` Pa·s. If you expect advection and
  see none, check that the *solid/liquid* components carry `κ` (gases never
  do) and that a pressure gradient actually develops.
- **Pressure spikes in nearly-solid cells.** The ideal-gas pressure carries a
  porosity divisor `P = (Σ_g ξ_j R_g T / M_j)/φ` with a floor `φ ≥ 1e-12`. In a
  cell that is almost entirely solid, the same gas mass yields a high pressure;
  that is physical (confined gas), not a bug.
- **Pressure BCs vent gas advectively — and only advectively.**
  `DirichletPressureBC`, `NeumannPressureBC`, and `ConvectivePressureBC` set
  the boundary-face Darcy velocity, which carries the advective venting flux
  `J = v_b·ξ/φ` out of the surface cell (outflow only — inflow advects the
  zero-concentration ambient). They contribute no diffusive film-law flux:
  that channel is governed by the diffusive mass BCs (`ConvectiveMassBC`,
  `DiffusiveMassBC`, …). Combine the two with `+`
  (`DiffusiveMassBC(h_m) + DirichletPressureBC()`) to get both escape paths;
  a pressure BC alone vents at a rate set by the permeability and the
  pressure gradient.
- **All gas species move together.** One effective `λ` serves every gas
  species, and advection carries them all at the same bulk velocity — there is
  no per-species molar-mass scaling (Technical Reference §9.2–9.3).

---

## 15.10 Parameter sweeps and parallelism

`parameter_sweep(problem, grid; modify, extract, parallel, solve_kwargs...)`
runs a grid of solves. The contract that keeps it correct and fast:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]

# modify MUST return a fresh problem and MUST NOT mutate the base problem.
# Material and PyrolysisProblem are built with keyword constructors — there is
# no copy/update constructor, so rebuild from the base problem's fields:
function modify(prob, p)
    base = prob.material
    new_rxn = Reaction(:decomp, :virgin => :gas; A = p.A, E = p.E, h = 200e3)
    new_mat = Material(;
        name                  = base.name,
        components            = base.components,
        reactions             = (new_rxn,),
        mixing                = base.mixing,
        lateral_shrinkage_law = base.lateral_shrinkage_law,
    )
    return PyrolysisProblem(;
        mesh         = prob.mesh,
        material     = new_mat,
        bc_set       = prob.bc_set,
        T_initial    = prob.T_initial,
        ξ_initial    = prob.ξ_initial,
        tspan        = prob.tspan,
        output_times = prob.output_times,
    )
end

# extract pulls small scalars/curves to minimise worker→head transfer:
extract(sol) = (peak_T = maximum(sol.surface_T),
                final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```

Failure modes:

- **`modify` mutating the base problem.** The immutable `ProblemDef`/residual
  are shared across workers; the mutable `Workspace` is per-worker. If `modify`
  edits the shared base in place, you corrupt every other point. Always build a
  new problem.
- **Returning full solutions from a big grid.** `extract` defaults to `identity`
  (the whole solution). For large grids that floods the head process — return a
  small NamedTuple.
- **Warm-start in a sweep.** You can pass `warm_start` through `solve_kwargs`,
  but reusing a neighbour's converged state across perturbed grid points can
  introduce artifacts in *sensitivity* studies (Technical Reference §17.7).
  Prefer cold starts when computing gradients.

`parallel` defaults to `nworkers() > 1` (uses `Distributed.pmap`); set
`parallel = false` to force sequential `map` for debugging. `parameter_sweep`
already runs each solve with `progress = false`.

---

## 15.11 Performance checklist

Work top to bottom; each item is cheap to try and the early ones dominate.

1. **Mesh size.** Cost scales with `n·(1+NC)` (+ the ALE `z` block and the
   WithChi `χ` block). Run a quick convergence check and use the smallest
   `n_cells` that resolves the surface temperature and MLR. Grade toward the
   surface with `create_stretched_mesh` (§15.5.4) rather than refining
   uniformly.
2. **Keep the default structured Jacobian.** `jacobian = nothing` →
   `Structured(DirectSolve())`, exact and sparse — far faster than colored FD.
   Only depletion forces the FD fallback (§15.4).
3. **Pick the right linear solver.** `SparspakFactorization` (the default) is
   ~28 % faster than KLU for the block-tridiagonal pattern. For large ALE, try
   `Structured(StructuredSolve())`, then `ApproxSolve()` (§15.7.4). Iterative
   `ILUGMRESFactorization(τ = 0.1, gmres_rtol = 1e-8, restart = 30)` and
   `ILUBiCGSTABFactorization(τ = 0.1, rtol = 1e-8)` exist for very large
   systems, but inexact solves cause Newton-step rejections on stiff transients
   — `solve` warns that Sparspak is usually 1.7–3.3× faster in practice. Treat
   ILU as a last resort for problems where direct factorization truly dominates.
4. **Loosen tolerances for exploration.** `abstol`/`reltol` are the most direct
   speed knob; tighten only for validation (§15.2.2).
5. **Choose the integrator.** `KenCarp4` is the default; a BDF (`FBDF`) can take
   larger steps on very stiff problems (§15.2.1).
6. **Right-size the radiation model.** `SURFACE_ABSORPTION` is O(1);
   `BEER_LAMBERT` adds an O(n) sweep and a non-local (upper-triangular) Jacobian
   coupling. Do not pay for in-depth absorption you do not need (§15.8.3).
7. **Warm-start continuation runs** (§15.2.4) and **parallelise sweeps**
   (§15.10).
8. **Disable extras you are not using.** `progress = false` removes the progress
   bar in batch/sweep contexts; leave `diagnostics = false` in production runs —
   it costs an extra ledger sweep every `diagnostics_stride` steps and is only
   needed for debugging (§15.1.2).

The package's allocation-free residual caches, `@generated` component loops, and
per-operator AD workspaces (Technical Reference §15.7, §15.11) are automatic; you
do not configure them. The levers above are the ones you control.

---

## 15.12 Where to look next

- Conservation ledgers and what each residual proves: Technical Reference §16;
  the post-processing/reading side in UG-9.
- Jacobian backends, strategies, and the extension stubs
  (`klu_factorization`, `pardiso_factorization`, `amg_preconditioner`,
  `forward_sensitivity`, `adjoint_sensitivity`, `verify_mixing_derivatives`):
  UG-11 and Technical Reference §15.
- Adaptivity (error indicators, monitors, `relocate_nodes!`, depletion):
  UG-10 and Technical Reference §14.
- The full keyword and default table for `solve`: UG-16.

The escape hatch for anything not on the curated top-level surface is the
`Pyrolysis.Internal` umbrella (`using Pyrolysis.Internal`), which re-exports
every submodule's public symbols — the mesh-quality, GCL, and remapping helpers
in this chapter live there. Symbols reached that way are not part of the stable
public API and may move between minor releases.


---

# 16. API and Parameter Reference

This chapter is the consolidated reference for the public API surface of
Pyrolysis.jl and every user-facing parameter and its default. Earlier chapters
teach *workflow*; this one is the lookup table you return to when you need the
exact constructor keyword, the precise default value, or the right place to
import a name from. Every type, function, keyword, and default shown here was
checked against the package source — where a name is reachable only through a
submodule (and not from `using Pyrolysis` alone), that is stated explicitly so
you never get an `UndefVarError` from a copied snippet. For the physics behind any
quantity, follow the "see Technical Reference §N" pointers.

Two coordinate/convention reminders that the whole API assumes (see Technical
Reference §2 and §3):

- `z = 0` is the **bottom / substrate** (boundary id 2, tag `:bottom`); `z = L`
  is the **top / exposed surface** (boundary id 1, tag `:top`). Heat enters from
  the top. Boundary heat and mass fluxes are **positive into the domain**.
- The **primary species state is the mass concentration** `ξ_j` in kg·m⁻³, not a
  mass fraction. The heat-of-reaction sign convention is **`h > 0` endothermic
  (cools)**, `h < 0` exothermic (heats) — see Technical Reference §6.

---

## 16.1 What `using Pyrolysis` gives you, and what it does not

The top-level `Pyrolysis` module exports a **curated** surface. Many useful
names live in submodules and are *not* re-exported at the top level. There are
three ways to reach a name:

1. **Top-level export** — available after `using Pyrolysis`. This is the stable
   public API; prefer it in production code.
2. **Submodule** — e.g. `Pyrolysis.Discretization.get_surface_temperature`, or
   import it: `using Pyrolysis.Discretization: get_surface_temperature`.
3. **`Pyrolysis.Internal`** — an umbrella namespace that re-exports *every*
   submodule's public symbol at once:

   ```julia
   using Pyrolysis            # curated top-level surface
   using Pyrolysis.Internal   # everything every submodule exports
   ```

   Symbols reached via `Internal` are **not** part of the stable API and may
   move between minor releases. `Pyrolysis.Experimental` is **not** included in
   `Internal`; access it explicitly (`Pyrolysis.Experimental.X`).

> Caveat that bites copy-pasters: the solution post-processors `get_profile`,
> `get_cell_positions`, `interpolate_solution`; the callbacks
> `create_depletion_callback`, `create_thickness_termination_callback`,
> `create_step_update_callback`; the iterative linear solvers
> `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`; the surface-temperature
> accessors `get_surface_temperature`, `get_back_face_temperature`; and the
> geometry accessors (`n_cells`, `cell_center`, `total_volume`, …) are exported
> by their submodules but **not** at the top level. Use `using Pyrolysis.Internal`
> or a submodule-qualified import to call them.

### 16.1.1 Submodule map

The package is a layered DAG of submodules. The table lists each submodule and
the kind of names it owns. Everything a submodule exports is reachable as
`Pyrolysis.<Submodule>.X` and is bundled into `Pyrolysis.Internal`.

| Submodule | Owns |
|---|---|
| `Materials` | `Material`, components, `Reaction`/`ReactionSet`, property-function types, `MixingRule`/`MixingSpec`/`MaterialMixing`, `RadiationModel` enum, physical constants |
| `BoundaryConditions` | all thermal/mass/pressure/mesh BC types, `BoundaryConditionSet`, `ambient_temperature` |
| `Geometry` | `Mesh1D`, `Node1D`/`Cell1D`/`Face1D`/`CellState`, `MeshTopology`/`MeshGeometry`, `create_uniform_mesh`/`create_stretched_mesh`, geometry accessors (`n_cells`, `cell_center`, `total_volume`, `domain_bounds`, `cell_size`, `boundary_face_ids`, …) |
| `Properties` | mixture/effective-property functions, `PropertyCache`, `RadiationCache` |
| `Physics` | kinetics, conduction, gas transport, radiation (surface/Beer–Lambert), volume change, pressure, ALE mesh motion, GCL, remapping |
| `Adaptivity` | error indicators, monitors/`relocate_nodes!`, quality metrics, surface-depletion merge machinery, `find_*` helpers |
| `Discretization` | finite-volume operators, flux limiters, boundary fluxes, `solve_boundary_temperature`, `get_surface_temperature`/`get_back_face_temperature`, ALE/Darcy fluxes |
| `Problem` | `PyrolysisProblem`, `PyrolysisSolution`, `ProblemDef`, `Workspace`, `SimulationMode` traits, `StateLayout`, validation, initial-condition helpers |
| `Residual` | `residual!`, column-average helpers |
| `Jacobian` | `JacobianBackend`/`JacobianSpec`, `Structured`/`DenseAD_ForwardDiff`, `LinearSolveStrategy` family, operators, structured couplings, `verify_jacobian_at` |
| `Diagnostics` | conservation trackers, `check_*`/`*_report`, `DiagnosticsCallback`, `DiagnosticsState`/`DiagnosticsRow` |
| `Solver` | `solve`, `parameter_sweep`, callbacks, ILU solvers, `get_profile`/`get_cell_positions`/`interpolate_solution` |
| `Experimental` | P1 radiation (`P1_QUASI_STEADY`, Marshak/Reflective/Dirichlet radiation BCs, the P1 residual), h-AMR (refinement/coarsening, `AMRController`) — **unstable**, not in `Internal` |

### 16.1.2 Extension stubs

Three families of functions are *stubs* in core Pyrolysis that gain concrete
methods only when an optional dependency is loaded. Calling a stub without its
extension raises an error telling you which package to add.

| Stub(s) | Triggered by | Provides |
|---|---|---|
| `klu_factorization` | `using KLU` | `LinearSolve.KLUFactorization()` under a stable name |
| `pardiso_factorization` | `using Pardiso` | `LinearSolve.MKLPardisoFactorize()` |
| `amg_preconditioner` | `using AlgebraicMultigrid` | smoothed-aggregation AMG preconditioner |
| `forward_sensitivity`, `adjoint_sensitivity` | `using SciMLSensitivity, Zygote` | local forward/adjoint sensitivities (§16.10) |
| `verify_mixing_derivatives` | `using Symbolics` | symbolic mixing-derivative check (test-only) |

Plotting recipes (`plot(sol; kind = …)`) attach when `using Plots` (or
`RecipesBase`) is loaded; see Chapter 9.

---

## 16.2 Materials, components, and property functions

See Technical Reference §4 (materials and property functions) and §5 (effective
properties). All names in this section are top-level exported.

### 16.2.1 Components

Three concrete component types encode phase at the type level. Each ships a
keyword convenience constructor.

```julia
SolidComponent(name::Symbol; ρ, c, k,
               ε = 0.9, α = 0.0, γ = 1.0, λ = 0.0, κ = 0.0,
               ρ_intrinsic = nothing, φ = nothing)

LiquidComponent(name::Symbol; ρ, c, k,
                λ = 2e-5, κ = 0.0, ε = 0.95, α = 0.0, γ = 0.0,
                ρ_intrinsic = nothing, φ = nothing)

GaseousComponent(name::Symbol; M, c, k, λ, ε = 0.0, α = 0.0, γ = 0.0)
```

| Keyword | Meaning (symbol) | Units | Default |
|---|---|---|---|
| `ρ` | bulk (pure-phase) density `ρ_j` | kg·m⁻³ | required (solid/liquid) |
| `M` | molar mass `M_j` (gas only; sets ideal-gas density `ρ_g = P_ref M/(R_g T)`) | kg·mol⁻¹ | required (gas) |
| `c` | specific heat `c_{p,j}` | J·kg⁻¹·K⁻¹ | required |
| `k` | thermal conductivity `k_j` | W·m⁻¹·K⁻¹ | required |
| `λ` | gas-transfer (diffusion) coefficient `λ_j` | m²·s⁻¹ | `0.0` (solid), `2e-5` (liquid), required (gas) |
| `κ` | permeability `κ_j` (solid/liquid only) | m² | `0.0` |
| `ε` | emissivity `ε_j` | – | `0.9` (solid), `0.95` (liquid), `0.0` (gas) |
| `α` | mass absorption coefficient `α_j` | m²·kg⁻¹ | `0.0` |
| `γ` | swelling factor `γ_j` | – | `1.0` (solid), `0.0` (liquid/gas) |
| `ρ_intrinsic` | skeletal/intrinsic density `ρ_{i,j}` (porosity formula only) | kg·m⁻³ | `= ρ` |
| `φ` | pure-phase porosity; sets `ρ_intrinsic = ρ/(1−φ)` (constant `ρ` only) | – | `nothing` |

`ρ`, `c`, `k`, `λ`, `κ` accept either a plain number or a property-function
object (§16.2.2); they are converted with `to_property` automatically. Pass
**either** `ρ_intrinsic` or `φ`, never both — supplying both throws an
`ArgumentError`. The bulk density `ρ` is what mixing rules and mixture density
use; `ρ_intrinsic` is used **only** in the porosity formula
`φ = clamp(1 − Σ ξ_j/ρ_{i,j}, 0, 1)` (see §5).

```julia
virgin = SolidComponent(:virgin, ρ = 530, c = 2000, k = 0.12, ε = 0.9, κ = 1e-15)
char   = SolidComponent(:char,   ρ = 150, c = 1500, k = 0.20, ε = 0.95, α = 5.0, κ = 1e-12)
water  = LiquidComponent(:H2O_l, ρ = 1000, c = 4180, k = 0.6)
gas    = GaseousComponent(:gas,  M = 0.044, c = 1800, k = 0.03, λ = 1e-4)
vapor  = GaseousComponent(:H2O,  M = 0.018, c = 2000, k = 0.025, λ = 2e-5)
```

Phase predicates (top-level exported): `is_solid(c)`, `is_liquid(c)`,
`is_gas(c)`, `phase(c)`, and `molar_mass(c)` (returns `0.0` for non-gas).

### 16.2.2 Property-function types

Every property field can be a constant number or one of these callable property
objects (all top-level exported). They implement `p(T)` (and `p(T, ξ)` for the
state-dependent one) and carry analytical or finite-difference derivatives used
in Jacobian assembly (see §4.5).

| Type | Constructor | Value `p(T)` |
|---|---|---|
| `ConstantProperty` | `ConstantProperty(value)` | `value` |
| `LinearProperty` | `LinearProperty(a, b)` | `a + b·T` |
| `PolynomialProperty` | `PolynomialProperty((c0, c1, …))` | `evalpoly(T, coeffs)` (low-to-high order) |
| `ArrheniusProperty` | `ArrheniusProperty(A, E, T_ref)` | `A·exp(−E/(R_g·T))` |
| `TableProperty` | `TableProperty((T1,…), (v1,…))` | linear interpolation; nearest-value extrapolation |
| `CallableProperty` | `CallableProperty(f)` | `f(T)` |
| `StateDependentProperty` | `StateDependentProperty(f)` | `f(T, ξ)` (errors if called as `p(T)`) |

`TableProperty` requires strictly increasing temperatures (checked at
construction); it uses a linear search for ≤8 points and a binary search above
that. `StateDependentProperty` is the type to use for moisture-dependent heat of
sorption and similar composition-dependent quantities.

The convenience converter `to_property` (used implicitly by the constructors)
maps `Real → ConstantProperty`, `Tuple → LinearProperty`,
`Function → CallableProperty`, and passes property objects through.

```julia
# Temperature-dependent conductivity that rises with T:
char = SolidComponent(:char, ρ = 150, c = 1500, k = LinearProperty(0.10, 1.5e-4))

# Tabulated heat capacity:
cp_tab = TableProperty((300.0, 500.0, 800.0), (1500.0, 1800.0, 2100.0))
virgin = SolidComponent(:virgin, ρ = 530, c = cp_tab, k = 0.12)
```

### 16.2.3 Reactions

`Reaction` (top-level exported) is a single-reactant Arrhenius reaction. Four
constructor forms cover numeric vs. symbolic component references and one vs.
two products. Symbolic forms return a `PendingReaction` that the `Material`
constructor resolves against the component names.

```julia
# one product, numeric indices
Reaction(name, reactant::Int => product::Int;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, validate_mass = true)

# two products, numeric indices
Reaction(name, reactant::Int => (p1::Int, p2::Int);
         A, E, h, yields::Tuple, n = 1.0, T_min = 0.0, T_max = Inf, validate_mass = true)

# one product, symbolic names
Reaction(name, reactant::Symbol => product::Symbol; A, E, h, n = 1.0, …)

# two products, symbolic names
Reaction(name, reactant::Symbol => (p1::Symbol, p2::Symbol); A, E, h, yields, n = 1.0, …)
```

| Keyword | Meaning (symbol) | Units | Default |
|---|---|---|---|
| `A` | Arrhenius pre-exponential `A_i` | s⁻¹ | required |
| `E` | activation energy `E_i` | J·mol⁻¹ | required |
| `h` | heat of reaction (number, tuple, or property) — **`h > 0` endothermic** | J·kg⁻¹ of reactant | required |
| `yields` | mass yields `(y1, y2)` of the two products (must sum to 1.0) | – | required for two-product form |
| `n` | reaction order `n_{i,reactant}` | – | `1.0` |
| `T_min` | lower temperature gate `T_min,i` (smooth tanh ramp) | K | `0.0` |
| `T_max` | upper temperature gate `T_max,i` (smooth tanh ramp) | K | `Inf` |
| `validate_mass` | enforce mass balance at construction | – | `true` |

Mass-balance convention: the single reactant is consumed with implicit mass 1.0;
product yields must sum to 1.0 (the one-product form uses yield 1.0
automatically). Each yield must additionally be ≥ 0 and the order must satisfy
`n ≥ 0` — both enforced unconditionally at construction; `validate_mass = false`
disables only the sum check. `h` accepts a number (`ConstantProperty`), a `(a, b)` tuple
(`LinearProperty`, i.e. `a + b·T`), or any property object — use
`StateDependentProperty` for moisture sorption. See §6 for the rate law, tanh
gates (≈10 K leakage caveat), and the depletion limiter.

```julia
# endothermic charring: virgin -> 0.2 char + 0.8 gas
rxn = Reaction(:charring, :virgin => (:char, :gas);
               A = 1e15, E = 200e3, h = 500e3, yields = (0.2, 0.8),
               T_min = 473.0, T_max = 1273.0)

# exothermic char oxidation (negative h releases heat)
rxn2 = Reaction(:char_oxidation, 1 => 2; A = 1e10, E = 160e3, h = -25e6)
```

`ReactionSet` is a tuple-based collection type (top-level exported); the common
path is to pass a plain tuple of reactions to `Material`.

### 16.2.4 Mixing rules: `MixingRule`, `MixingSpec`, `MaterialMixing`

Effective conductivity `k`, gas transfer `λ`, and permeability `κ` are each
mixed with their own rule (see §5). The rule enum (top-level exported) is:

```julia
@enum MixingRule PARALLEL SERIES WEIGHTED BRUGGEMAN CARMAN_KOZENY
```

| Value | Formula / use |
|---|---|
| `PARALLEL` | `k = Σ v_j k_j` (upper bound) |
| `SERIES` | `1/k = Σ v_j/k_j` (lower bound) |
| `WEIGHTED` | `k = β k_∥ + (1−β) k_series` |
| `BRUGGEMAN` | effective-medium theory; **`k` only** (errors for `λ` and `κ`) |
| `CARMAN_KOZENY` | porosity-based permeability `κ = min(ℓ²φ³/(180(1−φ)²), ℓ²)`; **`κ` only** |

`MixingSpec` bundles a rule with the weighting `β` and (for Carman–Kozeny) the
characteristic length:

```julia
MixingSpec(rule::MixingRule, β::Real = 0.5; length_scale::Real = NaN)
```

`β` must be in `[0, 1]` (default `0.5`). For `CARMAN_KOZENY` you **must** pass a
positive `length_scale` (the characteristic pore/particle diameter `ℓ` in
metres) via the keyword — note it is `length_scale`, not a positional `ℓ`.

`MaterialMixing` builds the `(k, λ, κ)` NamedTuple that `Material` expects.
Each channel accepts a `MixingSpec` or a bare `MixingRule` (in which case
`β = 0.5`):

```julia
MaterialMixing(; k = MixingSpec(WEIGHTED, 0.5),
                 λ = MixingSpec(PARALLEL, 0.5),
                 κ = MixingSpec(WEIGHTED, 0.5))
```

```julia
mixing = MaterialMixing(
    k = MixingSpec(WEIGHTED, 0.6),
    λ = PARALLEL,                              # bare rule -> β = 0.5
    κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),
)
```

### 16.2.5 Material

`Material` (top-level exported) assembles components, reactions, mixing, and an
optional lateral-shrinkage law.

```julia
Material(; name::Symbol,
           components::Tuple,
           reactions::Tuple = (),
           mixing = MaterialMixing(),
           lateral_shrinkage_law::Union{Function,Nothing} = nothing,
           depletion_limiter::DepletionLimiter = DepletionLimiter())
```

| Keyword | Meaning | Default |
|---|---|---|
| `name` | material identifier | required |
| `components` | tuple of component objects (heterogeneous; order defines index `j`) | required |
| `reactions` | tuple of `Reaction`/`PendingReaction` (symbolic names resolved here) | `()` |
| `mixing` | `(k, λ, κ)` NamedTuple from `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `χ̄ ↦ A(t)/A_0`; `nothing` ⇒ constant cross-section | `nothing` |
| `depletion_limiter` | per-material `DepletionLimiter` for the reaction-rate depletion roll-off (see §16.2.3 note and Technical Reference §6.3) | `DepletionLimiter()` |

Providing a non-`nothing` `lateral_shrinkage_law` switches the simulation into
the variable-cross-section (WithChi) mode under ALE (see §10 and §16.6). The law
is any function of the column-average pyrolysis progress `χ̄`, e.g.
`χ̄ -> 1 - 0.3*χ̄` for 30 % lateral shrink at full conversion.

```julia
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86),
        GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3),
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)
```

`DepletionLimiter(; threshold = 1.0, enabled = true)` — `threshold` [kg·m⁻³] is
the concentration scale of the smooth `tanh(ξ/threshold)` rate roll-off applied
per reactant; `enabled = false` removes the only guard that stops zero-order
reactions at ξ = 0. Exported by the `Physics` and `Materials` submodules (not at
the top level): `using Pyrolysis.Physics: DepletionLimiter`.

> Pyrolysis.jl does **not** provide a "copy-with-overrides" `Material(prob.material; …)`
> or `PyrolysisProblem(prob; …)` constructor. To vary a parameter (e.g. in a
> sweep), rebuild the full object from scratch — see §16.9.

### 16.2.6 `RadiationModel` enum

The volumetric radiation absorption model (top-level exported, except the P1
value) is an enum:

```julia
@enum RadiationModel NO_RADIATION SURFACE_ABSORPTION BEER_LAMBERT P1_QUASI_STEADY
```

| Value | Behaviour | Exported from |
|---|---|---|
| `NO_RADIATION` | no volumetric absorption (boundary radiative BCs still apply) | top level |
| `SURFACE_ABSORPTION` | absorbed flux `α·Φ` applied at the exposed surface via the boundary energy balance (opaque) | top level |
| `BEER_LAMBERT` | exponential in-depth absorption; non-local (upper-triangular) Jacobian | top level |
| `P1_QUASI_STEADY` | P1 diffusion approximation — **experimental** | `Pyrolysis.Experimental` |

`P1_QUASI_STEADY` is **not** wired into the core `solve` path; passing it to
`solve` raises an error (see §16.7 and §8). It lives in `Pyrolysis.Experimental`
with its own residual.

---

## 16.3 Meshes and geometry

See Technical Reference §13 (discretization) and the geometry notes in §13. Only
the two mesh creators and the topology/geometry view types are top-level
exported; the geometry accessors require `Pyrolysis.Internal` or the `Geometry`
submodule.

### 16.3.1 Mesh creators (top-level exported)

```julia
create_uniform_mesh(L::Real, n_cells::Int, ::Val{NC};
                    T_initial = 300.0, ξ_initial = nothing,
                    cross_section_area = 1.0)

create_stretched_mesh(L::Real, n_cells::Int, ::Val{NC}, stretch::Real;
                      T_initial = 300.0, ξ_initial = nothing,
                      cross_section_area = 1.0, stretch_from = :top)
```

| Argument / keyword | Meaning | Default |
|---|---|---|
| `L` | domain thickness `L` | required (m) |
| `n_cells` | number of finite-volume cells `n` | required |
| `Val(NC)` | number of components, as a `Val` type parameter | required |
| `stretch` | geometric stretch ratio `Δz_{i+1}/Δz_i` (positional, stretched mesh only); must be `≥ 1` (sub-unity throws — pick the refined end with `stretch_from`) | required |
| `T_initial` | uniform initial temperature for the seeded mesh | `300.0` K |
| `ξ_initial` | per-cell initial concentrations (number, vector, or `nothing`) | `nothing` |
| `cross_section_area` | initial cross-section `A_0` (gauge normalization in 1D) | `1.0` m² |
| `stretch_from` | `:top` (fine near surface) or `:bottom` | `:top` |

The third argument is a **`Val`**: `Val(NC)` where `NC` is the component count —
e.g. `create_uniform_mesh(0.01, 100, Val(2))` for a 1 cm, 100-cell, 2-component
mesh. The initial conditions seeded here are usually overwritten by the
problem's `T_initial`/`ξ_initial` closures, so you typically leave
`T_initial`/`ξ_initial` at their defaults and set them on `PyrolysisProblem`.

```julia
mesh   = create_uniform_mesh(0.01, 100, Val(2))                 # 1 cm, 100 cells
meshst = create_stretched_mesh(0.02, 80, Val(3), 1.03; stretch_from = :top)
```

### 16.3.2 Geometry accessors (not top-level; use `Internal` or `Geometry`)

After a solve you usually read geometry through the `PyrolysisSolution`
(§16.5). When you need to inspect the mesh directly, import from the submodule:

```julia
using Pyrolysis.Geometry: n_cells, n_nodes, cell_center, cell_size,
                          total_volume, domain_bounds, boundary_face_ids

n = n_cells(mesh)               # number of cells
zc = [cell_center(mesh, i) for i in 1:n]
(zmin, zmax) = domain_bounds(mesh)
```

`MeshTopology` and `MeshGeometry` (the immutable-connectivity and mutable-geometry
views used by the low-level residual layer) *are* top-level exported. The
component count of a material is `n_components(material)`, exported from
`Pyrolysis.Materials` (also in `Internal`).

### 16.3.3 State-vector layout

The flat ODE state is packed block-major as `[T | ξ_1..ξ_NC | z | χ]`, where the
`z` block exists only in ALE mode and the `χ` block only in WithChi mode (see §15
and §16.6). `StateLayout` (top-level exported) is the descriptor; `state_length`,
`pack_state!`, `create_state_vector`, `update_mesh_from_state!`, `block_view`,
and the accessors `get_temperature`, `get_concentration`, `get_node_position`,
`get_node_positions`, `get_thickness_from_state`, `column_chi_bar_from_state`,
`cross_section_area_from_state` are also top-level exported for users who index
raw state. Most users never touch these — `solve` and the solution object handle
packing for you. See §16.11 for the low-level residual surface.

---

## 16.4 Boundary and initial conditions

See Technical Reference §12. All concrete BC types listed here and
`BoundaryConditionSet` are top-level exported. The P1 radiation BCs (Marshak,
Reflective, Dirichlet) are **not**; they live in `Pyrolysis.Experimental`.

### 16.4.1 Thermal BCs

| Type | Constructor | Flux `q_BC` (positive = into domain) |
|---|---|---|
| `AdiabaticBC` | `AdiabaticBC()` | `0` |
| `HeatFluxBC` | `HeatFluxBC(flux)` (number or `t -> q`) | `flux` or `flux(t)` |
| `ConvectiveBC` | `ConvectiveBC(; h, T_inf)` | `h(t)·(T_∞ − T_s)` |
| `RadiativeBC` | `RadiativeBC(; ε = nothing, T_inf, F = 1.0)` | `F·ε_eff·σ·(T_∞⁴ − T_s⁴)` |
| `DirichletThermalBC` | `DirichletThermalBC(T)` | enforces `T_s = T` (penalty in Newton) |
| `RadiativeFluxBC` | `RadiativeFluxBC(flux; absorptivity = nothing)` or `RadiativeFluxBC(; flux, absorptivity = nothing)` | incident external radiation (model-dependent, below) |
| `CombinedThermalBC` | `bc1 + bc2 + …` | sum of components |

Notes and defaults:

- `ConvectiveBC` requires `h` and `T_inf` keywords; both accept a number or a
  function of `t`.
- `RadiativeBC`: `ε = nothing` (default) uses the composition-weighted effective
  emissivity `ε_eff` (Kirchhoff), tracking the surface as it chars; pass a
  number or `T -> ε` to override. `F` is the view factor, default `1.0`.
- `RadiativeFluxBC` is **incident external** radiation (e.g. a cone heater), not
  the surface's own re-radiation. `absorptivity = nothing` (default) uses
  `ε_eff` (Kirchhoff). Its behaviour depends on the `radiation_model` passed to
  `solve`: for `NO_RADIATION`/`SURFACE_ABSORPTION` it is applied as a surface
  flux; for `BEER_LAMBERT`/`P1_QUASI_STEADY` it returns **zero** boundary flux
  because the incident energy is deposited volumetrically (avoids
  double-counting). The flux can also be a time-varying function; `solve`
  re-samples it at every RHS evaluation.
- Combine BCs with `+`: `HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)`.

The surface temperature `T_s` is found by a Newton solve of the surface energy
balance `q_BC(T_s) = k_eff·(T_s − T_cell)/Δz_{1/2}` (clamped to `[200, 3000]` K,
tol `1e-10`, ≤50 iterations) — see §12.

### 16.4.2 Mass BCs

| Type | Constructor | Flux `J_j` (positive = into domain) |
|---|---|---|
| `ImpermeableBC` | `ImpermeableBC()` | `0` (default) |
| `MassFluxBC` | `MassFluxBC(flux)` (number or `t -> J`) | `flux` or `flux(t)` |
| `ConvectiveMassBC` | `ConvectiveMassBC(; h_m, ξ_inf = 0.0)` | `h_m(t)·(ξ_∞ − ξ_s)` |
| `DiffusiveMassBC` | `DiffusiveMassBC()` or `DiffusiveMassBC(h_m)` | internal-diffusion + external-resistance combined |
| `CombinedMassBC` | `bc1 + bc2 + …` | sum of components |

`DiffusiveMassBC()` defaults to `h_m = Inf` (no external resistance → pure
diffusive limit `J = λξ/d`); a finite `h_m` gives the combined-resistance form
`J = λ h_m ξ/(λ + h_m d)`. `ConvectiveMassBC` defaults `ξ_inf = 0.0` (escape to
clean atmosphere); `ξ_inf` also accepts `t -> ξ_∞(t)` or a per-component
`NTuple` (length = NC, entries ≥ 0).

### 16.4.3 Pressure BCs (Darcy mode)

| Type | Constructor | Meaning |
|---|---|---|
| `DirichletPressureBC` | `DirichletPressureBC(P)` or `DirichletPressureBC()` | fixed boundary pressure (default `101325` Pa) |
| `NeumannPressureBC` | `NeumannPressureBC(∂P∂n)` or `NeumannPressureBC()` | prescribed `∂P/∂n` (default `0.0` ⇒ impermeable) |
| `ConvectivePressureBC` | `ConvectivePressureBC(h_P, P_ext)` | Robin: `v_g·n = h_P(P_int − P_ext)` |

`ConvectivePressureBC` takes `h_P` (pressure-transfer coefficient, m·Pa⁻¹·s⁻¹)
and `P_ext` positionally. Pressure BCs only act when the material enables Darcy
flow (any solid/liquid `κ > 0`); see §9.

### 16.4.4 Mesh motion (ALE)

There are no mesh-motion BC types. ALE mesh motion is governed entirely by the
residual: the bottom node is pinned (`w_1 = 0`) and node velocities integrate
the volumetric strain upward, so the surface follows the material automatically
(§7.3; Technical Reference §11.3). Enable it with `use_ale = true` — nothing to
construct.

### 16.4.5 `BoundaryConditionSet`

The container that pairs `:top` and `:bottom` BCs (top-level exported):

```julia
BoundaryConditionSet(; thermal::NamedTuple = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                       mass::NamedTuple    = (top = ImpermeableBC(), bottom = ImpermeableBC()))
```

User-supplied partial NamedTuples are merged with the defaults, then validated.
Tags map to ids `:top → 1`, `:bottom → 2`.

```julia
bcs = BoundaryConditionSet(
    thermal = (
        top    = RadiativeFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0) +
                 ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (
        top    = DiffusiveMassBC(),     # gas escapes at the surface
        bottom = ImpermeableBC(),
    ),
)
```

### 16.4.6 Initial conditions

Initial conditions are set on `PyrolysisProblem` (§16.5.1) via the `T_initial`
and `ξ_initial` closures (functions of `z`). When `ξ_initial` is omitted, the
default seeds the **first solid component at its reference density** and all
others to zero. `T_initial` defaults to `z -> 300.0`. Profile helpers
(`step_profile`, `linear_profile`, `gaussian_profile`, `tanh_profile`) live in
`Pyrolysis.Problem`.

> **Transpiration caveat.** The surface gas-transpiration enthalpy term is
> incompatible with the default `S_conv` gas-convection convention (they
> double-count). `solve` rejects `use_transpiration_bc = true`. Leave it off.

---

## 16.5 The problem and solution objects

See Technical Reference §3 and §15. `PyrolysisProblem` and `PyrolysisSolution`
are top-level exported.

### 16.5.1 `PyrolysisProblem`

```julia
PyrolysisProblem(; mesh, material, bc_set::BoundaryConditionSet,
                   T_initial = z -> 300.0,
                   ξ_initial = nothing,
                   tspan::Tuple{Real,Real},
                   dt_initial::Real = 1e-4,
                   output_times::Vector{<:Real} = Float64[])

# convenience positional form:
PyrolysisProblem(mesh, material, bc_set, tspan; kwargs...)
```

| Keyword | Meaning | Default |
|---|---|---|
| `mesh` | a `Mesh1D` from a creator | required |
| `material` | a `Material` | required |
| `bc_set` | a `BoundaryConditionSet` | required |
| `T_initial` | `z -> T0` initial temperature closure | `z -> 300.0` |
| `ξ_initial` | `z -> ξ_j` closure, a vector of such closures, or `nothing` | `nothing` (first solid at ref density) |
| `tspan` | `(t_start, t_end)` | required (s) |
| `dt_initial` | initial timestep hint; also sets the auto output-time count | `1e-4` |
| `output_times` | snapshot times; empty ⇒ auto 11–101 evenly spaced over `tspan` | `Float64[]` |

```julia
problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],   # virgin density, gas zero
    tspan = (0.0, 100.0),
)
```

`validate_problem(problem; strict = false)` (in `Pyrolysis.Problem`) checks mesh,
material stoichiometry, BC coverage, time span, and initial-condition counts;
`solve` runs it automatically unless `validate = false`. `validate_material(m)`
(in `Pyrolysis.Materials`) validates a material in isolation.

### 16.5.2 `PyrolysisSolution`

`solve` returns a `PyrolysisSolution`. Core fields (all read directly):

| Field | Type | Meaning |
|---|---|---|
| `problem` | `PyrolysisProblem` | the originating problem |
| `t` | `Vector{Float64}` | saved times [s] |
| `T` | `Matrix{Float64}` | temperatures, `(n_cells, n_steps)` [K] |
| `ξ` | `Array{Float64,3}` | concentrations, `(NC, n_cells, n_steps)` [kg·m⁻³] |
| `z` | `Matrix{Float64}` | node positions, `(n_nodes, n_steps)` [m]; empty `0×0` in Eulerian mode |
| `surface_T` | `Vector{Float64}` | solved surface temperature per step [K] |
| `mass_loss_rate` | `Vector{Float64}` | MLR per **current** area `A(t)` (top + permeable-bottom) [kg·m⁻²·s⁻¹] |
| `energy_residual` | `Vector{Float64}` | energy-conservation residual (zeros if tracker disabled) |
| `retcode` | `Symbol` | integrator return code (e.g. `:Success`, `:Terminated`) |
| `solve_time` | `Float64` | wall-clock solve time [s] |
| `raw_solution` | (untyped) | the underlying ODE solution (interpolation source) |
| `extras` | `NamedTuple` | opt-in diagnostics (below) |

`extras` is a NamedTuple of vectors, written in place during extraction:
`gas_generation_rate`, `back_temperature`, `total_mass`, `thickness_history`,
`surface_mass_flux`, `surface_heat_flux_net`, `surface_heat_flux_absorbed`,
`surface_heat_flux_reradiation`, `surface_heat_flux_convection`,
`surface_heat_flux_conduction`, `surface_heat_flux_transpiration`,
`volumetric_absorption`,
`surface_emissivity`, `surface_conductivity`, `back_heat_flux`, `MLR_gauge`,
`cross_section_area_history`, `diagnostics_log`.

```julia
sol = solve(problem)
times    = sol.t
Tsurf    = sol.surface_T                       # K, per saved time
mlr_top  = sol.extras.surface_mass_flux        # kg/(m^2 s), top face only
mlrgauge = sol.extras.MLR_gauge                # normalized to A_0 (variable cross-section)
Tprofile = sol.T[:, end]                       # final temperature profile
ξvirgin  = sol.ξ[1, :, end]                    # final virgin concentration profile
```

`MLR_gauge` equals `mass_loss_rate` scaled by `A(t)/A_0`; for constant
cross-section the two coincide.

### 16.5.3 Solution post-processors (not top-level; use `Internal`/`Solver`)

```julia
using Pyrolysis.Internal   # or: using Pyrolysis.Solver: get_profile, …

get_profile(sol, field::Symbol, time_idx::Int)   # :T, :ξ, or :ξ1/:ξ2/…
get_cell_positions(sol)                          # cell-center positions from the mesh
interpolate_solution(sol, t::Real)               # full raw state vector at arbitrary t
```

`get_profile(sol, :T, i)` returns `sol.T[:, i]`; `get_profile(sol, :ξ, i)`
returns the full `NC × n_cells` slice; `get_profile(sol, :ξ1, i)` returns
component 1. The surface/back-face temperature accessors read the *mesh's* solved
boundary state, not the solution:

```julia
using Pyrolysis.Discretization: get_surface_temperature, get_back_face_temperature,
                                check_boundary_energy_balance
Ts_final    = get_surface_temperature(problem.mesh)   # mesh holds last accepted step
Tback_final = get_back_face_temperature(problem.mesh)
```

For per-step surface temperature, prefer `sol.surface_T`.

---

## 16.6 Simulation modes (traits)

See Technical Reference §15. A simulation's behaviour is a bundled trait
`SimulationMode{Mesh,Rad,Trans,Chi}` derived automatically by `solve`/`to_problem_def`
from your options and material — you rarely construct it by hand. `SimulationMode`
is top-level exported; the **trait subtypes are deliberately not exported** (to
keep the namespace clean). Reach them via `Pyrolysis.Problem`:

| Axis | Subtypes | Selected by |
|---|---|---|
| `MeshMode` | `Eulerian`, `ALE` | `use_ale` kwarg of `solve` |
| `RadiationMode` | `NoRadiation`, `SurfaceAbsorption`, `BeerLambert` (P1 has no trait subtype — it is quarantined in `Experimental`) | `radiation_model` kwarg |
| `TransportMode` | `Fickian`, `DarcyFick` | material has any `κ > 0` (Darcy) vs. not (Fickian) |
| `ChiMode` | `NoChi`, `WithChi` | material has a `lateral_shrinkage_law` (under ALE) |

```julia
using Pyrolysis.Problem: Eulerian, ALE, Fickian, DarcyFick, NoChi, WithChi,
                         mesh_mode, radiation_mode, transport_mode, chi_mode
```

Rules of thumb for choosing modes via `solve` options:

- TGA-like / fixed slab: defaults (Eulerian, Fickian, `NO_RADIATION`).
- Cone calorimeter with surface recession: `use_ale = true`,
  `radiation_model = SURFACE_ABSORPTION` (opaque) or `BEER_LAMBERT`
  (semi-transparent).
- Pressure-driven char / pressure-treated wood: give components `κ > 0` ⇒
  Darcy–Fick automatically; add pressure BCs.
- Intumescent / lateral shrink: give the material a `lateral_shrinkage_law` and
  run with `use_ale = true` ⇒ WithChi.

---

## 16.7 `solve`: full keyword reference

See Technical Reference §15. `solve` is top-level exported. It rejects unknown
keywords (typo protection), so only the keywords below are valid.

```julia
solve(problem::PyrolysisProblem; kwargs...) -> PyrolysisSolution
```

| Keyword | Type | Default | Purpose |
|---|---|---|---|
| `integrator` | type / instance / `nothing` | `nothing` → `KenCarp4(…)` | ODE algorithm (constructor or fully-built instance) |
| `linear_solver` | factorization / `nothing` | `nothing` → `SparspakFactorization()` | sparse/dense linear solver (only paired with a *constructor* `integrator`) |
| `jacobian` | `JacobianSpec` / `nothing` | `nothing` → `Structured(strategy = DirectSolve())` | pluggable Jacobian backend |
| `verify_jacobian` | `NamedTuple` / `nothing` | `nothing` | `(rtol=…, atol=…, sample_times=[…])` to compare backend vs. `DenseAD_ForwardDiff` |
| `abstol` | `Real` | `1e-8` | absolute ODE tolerance |
| `reltol` | `Real` | `1e-6` | relative ODE tolerance |
| `sensealg` | any | `nothing` | forwarded to the ODE solve (used by sensitivity ext) |
| `maxiters` | `Int` | `1000000` | max integration steps |
| `dtmin` | `Real` | `1e-12` | minimum timestep [s] |
| `dtmax` | `Real` | `Inf` | maximum timestep [s] |
| `saveat` | vector | `problem.output_times` | snapshot times |
| `callback` | callback / `nothing` | `nothing` | extra user callback(s) |
| `progress` | `Bool` | `true` | progress bar |
| `validate` | `Bool` | `true` | run `validate_problem` first |
| `use_ale` | `Bool` | `false` | enable ALE moving mesh |
| `min_thickness` | `Real` | `1e-6` | ALE termination thickness [m] |
| `handle_depletion` | `Bool` | `false` | merge depleted surface cells (requires `use_ale = true`) |
| `depletion_threshold` | `Real` | `0.05` | cell-thickness fraction below which a cell merges |
| `min_cells` | `Int` | `2` | minimum cell count before terminating instead of merging |
| `radiation_model` | `RadiationModel` | `NO_RADIATION` | `NO_RADIATION` / `SURFACE_ABSORPTION` / `BEER_LAMBERT` |
| `diagnostics` | `Bool` | `false` | enable conservation-diagnostics snapshots |
| `diagnostics_stride` | `Int` | `10` | snapshot every N steps |
| `energy_form` | `Symbol` | `:standard` | `:standard` or `:with_generation_sink` |
| `use_transpiration_bc` | `Bool` | `false` | **rejected if `true`** (double-counts `S_conv`) |
| `warm_start` | `PyrolysisSolution` / `nothing` | `nothing` | reuse a prior solution's final state as initial data |

Key behaviours and constraints:

- **Default integrator / linear solver**: `KenCarp4(autodiff = AutoFiniteDiff(),
  linsolve = SparspakFactorization())`. KenCarp4 is an SDIRK method suited to the
  moderately stiff Arrhenius system. Pass a constructor (e.g. `integrator = KenCarp4`)
  to combine with `linear_solver`, or a fully-built instance
  (e.g. `integrator = KenCarp4(linsolve = ILUGMRESFactorization())`) on its own.
- **Default Jacobian**: `Structured(strategy = DirectSolve())`. Pass a
  `JacobianSpec` to opt into other strategies (§16.8).
- **`radiation_model = P1_QUASI_STEADY`** is rejected (P1 is experimental).
- **`energy_form`** must be `:standard` (Formulation A, the FDS/Gpyro/ThermaKin
  convention) or `:with_generation_sink` (adds the `−Σ Δh_g ∇·J_g` cell source).
- **`handle_depletion = true` is incompatible with a pluggable `jacobian`**: the
  state vector resizes mid-solve, which a pre-built `JacobianSpec` cannot survive.
  `solve` errors if both are set; with depletion alone it falls back to the
  integrator's colored finite-difference Jacobian and switches the linear solver
  to a dense LU.
- **MLR gauge note**: in variable-cross-section runs, `sol.mass_loss_rate` is per
  current area and `sol.extras.MLR_gauge` is normalized to `A_0`.

```julia
# Cone calorimeter, opaque, ALE surface recession, default Jacobian:
sol = solve(problem;
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            min_thickness = 0.5e-3,
            abstol = 1e-9, reltol = 1e-7,
            integrator = KenCarp4(linsolve = SparspakFactorization()))
```

---

## 16.8 Jacobian and linear-solver configuration

See Technical Reference §15. These names are top-level exported.

### 16.8.1 Backends

| Name | Role |
|---|---|
| `Structured` | production operator-based backend; `Structured(; strategy = DirectSolve())` |
| `DenseAD_ForwardDiff` | dense whole-residual AD reference — verification only, not for production |
| `JacobianBackend` | abstract supertype |
| `JacobianSpec` | the single user-facing handle that bundles a backend + sparsity + caches |

Build a spec from a problem (the common form) or from a `(def, ws)` pair:

```julia
JacobianSpec(backend::JacobianBackend, problem::PyrolysisProblem; kwargs...)
JacobianSpec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...)

spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem)
sol  = solve(problem; jacobian = spec, use_ale = true)
```

### 16.8.2 Linear-solve strategies

`Structured` takes a `LinearSolveStrategy` (all top-level exported):

| Strategy | Behaviour | When |
|---|---|---|
| `DirectSolve` | materialize `A + Σ densify(couplings)`, factor (exact) | default; small/typical problems |
| `StructuredSolve` | sparse LU on `A` + preconditioned Richardson (exact, 2–5 iters) | typical, faster assembly |
| `ApproxSolve` | drop the structured ALE couplings (Schur-style approximation) | large ALE problems (~8× speedup; less-accurate Newton step) |
| `MatrixFreeSolve` | GMRES over the operator action (no materialization) | research/experimental |

`ApproxSolve` is marked approximate via the `is_approximate` trait; the
verifier compares it in masked (in-support) mode. For the Douglas-fir cone case,
`Structured(strategy = ApproxSolve())` gives roughly an 8× assembly speedup over
the full Jacobian at solver-tolerance-equivalent results.

### 16.8.3 Iterative linear solvers (not top-level; use `Internal`/`Solver`)

```julia
using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization

ILUGMRESFactorization(; τ = 0.1, gmres_rtol = 1e-8, restart = 30)
ILUBiCGSTABFactorization(; τ = 0.1, rtol = 1e-8)

sol = solve(problem;
            integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)))
```

| Keyword | Meaning | Default |
|---|---|---|
| `τ` | ILU drop tolerance `τ_ILU` | `0.1` |
| `gmres_rtol` / `rtol` | Krylov relative tolerance | `1e-8` |
| `restart` | GMRES restart frequency | `30` (GMRES only) |

Both fall back to diagonal preconditioning if the ILU factorization is singular
during a transient. Use them for large ALE problems where direct factorization
dominates cost (trade-off: inexact solves can cause Newton rejections on sharp
Arrhenius transients).

### 16.8.4 Verification

```julia
verify_jacobian_at(spec, def, ws, u, t; reference = DenseAD_ForwardDiff(), rtol, atol)
is_approximate(backend)            # trait: true for ApproxSolve
benchmark_jacobian_backends(...)   # timing harness
```

Or run verification *during* a solve at chosen times:

```julia
sol = solve(problem; verify_jacobian = (rtol = 1e-6, atol = 1e-10,
                                        sample_times = [0.0, 1.0, 10.0]))
```

It warns (does not abort) when the active backend deviates from the dense AD
reference.

### 16.8.5 Solver extensions (vendor factorizations)

Load the dependency, then pass the factorization to a constructor `integrator`:

```julia
using Pyrolysis, KLU
sol = solve(problem; integrator = KenCarp4(linsolve = Pyrolysis.klu_factorization()))
```

`Pyrolysis.pardiso_factorization()` (Intel MKL Pardiso, `using Pardiso`) and
`Pyrolysis.amg_preconditioner(A)` (`using AlgebraicMultigrid`) work the same
way. Calling a stub without its extension loaded raises an error naming the
package to add.

> **MKL ipiv caveat.** When `handle_depletion = true` resizes the state vector,
> avoid factorizations with fixed pivot buffers (MKL Pardiso). With depletion
> the solver already uses a dense, resize-safe `LUFactorization()` and a
> finite-difference Jacobian (explicit sparse choices are rebuilt to dense LU
> with a warning), so this is only a concern if you pass a fully built
> integrator with a non-default `linsolve`.

---

## 16.9 Parameter sweeps

See Chapter 12. `parameter_sweep` is top-level exported.

```julia
parameter_sweep(problem, grid;
                modify::Function = _default_sweep_modify,
                extract::Function = identity,
                parallel::Bool = nworkers() > 1,
                solve_kwargs...) -> Vector
```

| Argument | Meaning |
|---|---|
| `problem` | base problem (never mutated) |
| `grid` | any iterable of parameter points (`Vector{NamedTuple}`, `Iterators.product`, …) |
| `modify` | `(base_problem, point) -> new_problem` — must return a **fresh** problem |
| `extract` | `solution -> summary` (default `identity`; prefer a small NamedTuple for big grids) |
| `parallel` | use `Distributed.pmap` when more than one worker is available |
| `solve_kwargs...` | forwarded to each `solve` call |

Because there is no copy-with-overrides constructor, `modify` must rebuild the
objects it changes from scratch:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

function modify(prob, p)
    comps = prob.material.components
    new_rxn = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 200e3)
    new_mat = Material(name = prob.material.name,
                       components = comps,
                       reactions = (new_rxn,),
                       mixing = prob.material.mixing)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat,
                            bc_set = prob.bc_set, tspan = prob.tspan,
                            output_times = prob.output_times)
end

extract(sol) = (peak_T = maximum(sol.surface_T),
                final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify = modify, extract = extract)
```

---

## 16.10 Sensitivity analysis (extension)

See Chapter 13 and Technical Reference §17. `forward_sensitivity` and
`adjoint_sensitivity` are top-level exported **stubs** — they gain methods only
under `using SciMLSensitivity, Zygote` (which loads `PyrolysisSciMLSensitivityExt`).

```julia
forward_sensitivity(base_problem, θ::AbstractVector, p_inject::Function;
                    output_fn::Function = identity, solve_kwargs...)
    -> (output, ∂output_∂θ)

adjoint_sensitivity(base_problem, θ::AbstractVector, p_inject::Function, loss_fn::Function;
                    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
                    solve_kwargs...)
    -> (loss, grad)
```

| Argument | Meaning |
|---|---|
| `base_problem` | the unperturbed `PyrolysisProblem` |
| `θ` | parameter vector to differentiate with respect to |
| `p_inject` | `(base_problem, θ) -> problem'` — injects `θ` into a fresh problem |
| `output_fn` | `solution -> y` reducing a solution to the output(s) (forward; default `identity`) |
| `loss_fn` | `solution -> scalar` loss (adjoint) |
| `sensealg` | reverse-mode algorithm (adjoint; default `InterpolatingAdjoint(…)`) |

Forward mode (ForwardDiff) suits small `θ`; adjoint mode (Zygote +
SciMLSensitivity) suits large `θ` / scalar loss (inverse problems). `p_inject` is
the same pattern as `modify` in §16.9 — rebuild the perturbed problem from
`θ`. Sensitivities are reliable thanks to the AD-clean design (tanh gates,
midpoint enthalpy, smooth limiters); avoid warm starts and the depletion
callback when differentiating (non-smooth).

```julia
using Pyrolysis, SciMLSensitivity, Zygote

θ0 = [8.5e12, 188e3]   # A, E
p_inject = (prob, θ) -> begin
    rxn = Reaction(:decomp, 1 => 2; A = θ[1], E = θ[2], h = 870e3)
    mat = Material(name = :PMMA, components = prob.material.components,
                   reactions = (rxn,), mixing = prob.material.mixing)
    PyrolysisProblem(mesh = prob.mesh, material = mat, bc_set = prob.bc_set,
                     tspan = prob.tspan, output_times = prob.output_times)
end

out, J = forward_sensitivity(problem, θ0, p_inject;
                             output_fn = sol -> [maximum(sol.surface_T)])
```

---

## 16.11 Adaptivity and callbacks

See Technical Reference §14. The depletion path is production-ready;
h-refinement/coarsening and the `AMRController` are in `Pyrolysis.Experimental`.

### 16.11.1 Surface-depletion via `solve`

The simplest path is the `solve` keywords (§16.7): `handle_depletion = true`
(needs `use_ale = true`), `depletion_threshold` (default `0.05`), `min_cells`
(default `2`). Remember this falls back to a finite-difference Jacobian.

### 16.11.2 Callbacks (not top-level; use `Internal`/`Solver`)

```julia
using Pyrolysis.Solver: create_depletion_callback, create_thickness_termination_callback,
                        create_step_update_callback, DepletionCallbackState,
                        get_depletion_callback_stats

create_thickness_termination_callback(sv_ref, min_thickness::Float64)
create_depletion_callback(ws, mesh, material, sv_or_state; threshold = 0.05, min_cells = 2)
create_step_update_callback(ws, material)
get_depletion_callback_stats(cb)                  # merge history
```

You normally do not build these by hand — `solve` assembles them from its
keywords. They are exposed for custom solver loops.

### 16.11.3 Error indicators, monitors, quality metrics (Adaptivity submodule)

For manual r-refinement / mesh-quality work, import from `Pyrolysis.Adaptivity`:

```julia
using Pyrolysis.Adaptivity: GradientIndicator, CurvatureIndicator,
                            ReactionRateIndicator, InterfaceProximityIndicator,
                            CompositeIndicator, compute_indicator,
                            GradientMonitor, InterfaceMonitor, CompositeMonitor,
                            relocate_nodes!, AspectRatioMetric, SmoothnessMetric,
                            JacobianMetric, MinCellSizeMetric, CompositeQualityMetric,
                            check_quality, enforce_quality!, mesh_quality_summary
```

`relocate_nodes!(mesh; monitor = GradientMonitor(), max_iter = 10, relaxation = 0.5,
tol = 1e-3, h_min = 1e-8)` moves interior nodes by equidistribution (boundary
nodes are pinned). Quality metrics carry defaults (e.g. `AspectRatioMetric(2.0)`,
`SmoothnessMetric(0.5)`, `MinCellSizeMetric(1e-8)`).

### 16.11.4 Experimental h-AMR

`refine_cell!`, `perform_h_refinement!`, `coarsen_cells!`, `AMRController`,
`adapt_mesh!`, `create_default_controller`, etc. are in `Pyrolysis.Experimental`
and are unstable. The active-cell bookkeeping requires iterating the mask, not
`eachindex`; treat this path as research-grade.

---

## 16.12 Conservation diagnostics

See Technical Reference §16. Enable snapshots inline with `solve(problem;
diagnostics = true, diagnostics_stride = N)` — rows land in
`sol.extras.diagnostics_log`. The tracker types and check functions are in
`Pyrolysis.Diagnostics` (and in `Internal`).

```julia
using Pyrolysis.Internal   # or: using Pyrolysis.Diagnostics

# trackers
MassConservationTracker, EnergyConservationTracker, ReactionMassFlowTracker
DiagnosticsState, DiagnosticsRow, DiagnosticsCallback

# instantaneous quantities
compute_domain_mass(mesh, material)             # NTuple of per-component mass [kg]
compute_sensible_energy(mesh, material)         # matrix-only E_matrix [J]
compute_total_sensible_energy(mesh, material)   # E_total (matrix + gas) [J]

# checks (return NamedTuples)
check_mass_conservation(tracker, mesh, material; tol = 1e-10)
check_energy_conservation(tracker, mesh, material; tol = 1e-6)
check_reaction_mass_balance(tracker; tol = 1e-10)

# human-readable reports (return String)
mass_conservation_report(tracker, mesh, material)
energy_conservation_report(tracker, mesh, material)
reaction_mass_flow_report(tracker, material)
```

The energy diagnostics keep two ledgers: a **discrete** form that closes to
integrator tolerance by construction (audits bookkeeping) and a **physical**
form that quantifies the gas-storage gap `ΔE_gas = E_total − E_matrix` — a
runtime measure of the quasi-steady-gas approximation cost (canonically < 1 %).
A `DiagnosticsRow` carries `t`, `E_matrix`, `E_total`, `discrete_residual`,
`discrete_relative`, `physical_residual`, `physical_relative`, `gas_storage_gap`,
`mass_max_relative_error`.

The low-level residual surface — `ProblemDef`, `Workspace`, `residual!`,
`to_problem_def`, `build_workspace`, `build_ad_workspace`, `StateLayout` — is
top-level exported for users building custom solver or sensitivity loops. The
canonical call is `residual!(du, u, def, t, ws)`; build `(def, ws)` with
`def = to_problem_def(problem; use_ale, radiation_model, energy_form)` (plus,
for direct Beer–Lambert use, `radiation_intensity` and `radiation_kirchhoff` —
set `radiation_kirchhoff = false` if you have already folded an explicit
absorptivity into the intensity, else the live surface `ε_eff` is applied at
injection) and
`ws = build_workspace(problem, def)`. See Technical Reference §15.

---

## 16.13 Default-value reference

A single place to look up the defaults baked into the package. Physics is in the
noted TR chapter.

### `solve` keyword defaults

| Keyword | Default |
|---|---|
| `abstol` | `1e-8` |
| `reltol` | `1e-6` |
| `maxiters` | `1000000` |
| `dtmin` / `dtmax` | `1e-12` / `Inf` |
| `progress` | `true` |
| `validate` | `true` |
| `use_ale` | `false` |
| `min_thickness` | `1e-6` m |
| `handle_depletion` | `false` |
| `depletion_threshold` | `0.05` |
| `min_cells` | `2` |
| `radiation_model` | `NO_RADIATION` |
| `diagnostics` | `false` |
| `diagnostics_stride` | `10` |
| `energy_form` | `:standard` |
| `use_transpiration_bc` | `false` (rejected if `true`) |
| default integrator | `KenCarp4(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())` |
| default Jacobian | `Structured(strategy = DirectSolve())` |

### Component / material / reaction defaults

| Quantity | Default |
|---|---|
| `SolidComponent` `ε`, `α`, `γ` | `0.9`, `0.0`, `1.0` |
| `SolidComponent` `λ`, `κ` | `0.0`, `0.0` |
| `LiquidComponent` `λ` | `2e-5` m²·s⁻¹ |
| `LiquidComponent`/`GaseousComponent` `ε`, `γ` | `0.0`, `0.0` |
| `ρ_intrinsic` (when omitted) | `= ρ` (no pure-phase pores) |
| `Reaction` `n`, `T_min`, `T_max`, `validate_mass` | `1.0`, `0.0`, `Inf`, `true` |
| `MixingSpec` `β` | `0.5` |
| `MaterialMixing` `k`, `λ`, `κ` | `WEIGHTED`, `PARALLEL`, `WEIGHTED` (all `β = 0.5`) |
| `Material` `reactions`, `lateral_shrinkage_law` | `()`, `nothing` |

### BC defaults

| Quantity | Default |
|---|---|
| `RadiativeBC` `ε`, `F` | `nothing` (⇒ `ε_eff`), `1.0` |
| `RadiativeFluxBC` `absorptivity` | `nothing` (⇒ `ε_eff`, Kirchhoff) |
| `ConvectiveMassBC` `ξ_inf` | `0.0` |
| `DiffusiveMassBC` `h_m` | `Inf` (no external resistance) |
| `DirichletPressureBC` `P` | `101325` Pa |
| `NeumannPressureBC` `∂P∂n` | `0.0` |
| `BoundaryConditionSet` thermal / mass | `(top = AdiabaticBC(), bottom = AdiabaticBC())` / `(top = ImpermeableBC(), bottom = ImpermeableBC())` |

### Problem / mesh defaults

| Quantity | Default |
|---|---|
| `PyrolysisProblem` `T_initial` | `z -> 300.0` |
| `PyrolysisProblem` `ξ_initial` | first solid at ref density, others 0 |
| `PyrolysisProblem` `dt_initial` | `1e-4` s |
| `PyrolysisProblem` `output_times` | auto 11–101 evenly spaced over `tspan` |
| mesh creators `cross_section_area` (`A_0`) | `1.0` m² |
| mesh creators `T_initial` | `300.0` K |
| `create_stretched_mesh` `stretch_from` | `:top` |

### Linear-solver / Jacobian defaults

| Quantity | Default |
|---|---|
| `ILUGMRESFactorization` `τ`, `gmres_rtol`, `restart` | `0.1`, `1e-8`, `30` |
| `ILUBiCGSTABFactorization` `τ`, `rtol` | `0.1`, `1e-8` |
| `Structured` `strategy` | `DirectSolve()` |

### Kernel constants and limits

| Quantity (symbol) | Default | TR |
|---|---|---|
| depletion limiter threshold | `1.0` kg·m⁻³ | §6 |
| surface-temperature Newton clamp `[T_min, T_max]` | `[200, 3000]` K | §12 |
| surface-temperature Newton tol / max-iter | `1e-10` / `50` | §12 |
| porosity floor (ideal-gas pressure) | `1e-12` | §9 |
| matrix ρc_p floor (`RHO_CP_FLOOR`) | `1e3` J·m⁻³·K⁻¹ | §5 |
| Sutherland constant `S` / reference viscosity `μ_ref` | `110.4` K / `1.716e-5` Pa·s (at 273.15 K) | §5 |
| Kozeny constant | `180` | §5 |
| reference temperature `T_ref` | `298.15` K | §2 |
| reference pressure `P_ref` | `101325` Pa | §2 |
| gas constant `R_g` | `8.31446261815324` J·mol⁻¹·K⁻¹ | §2 |
| Stefan–Boltzmann `σ` | `5.670374419e-8` W·m⁻²·K⁻⁴ | §2 |
| GCL tolerance | `1e-12` | §11 |
| mass / reaction conservation tol | `1e-10` | §16 |
| energy conservation tol | `1e-6` | §16 |

---

## 16.14 Quick import cheat-sheet

| You want… | Reach it via |
|---|---|
| `Material`, components, `Reaction`, property types, `MixingSpec`, `RadiationModel` | top level (`using Pyrolysis`) |
| BC types, `BoundaryConditionSet` | top level |
| `create_uniform_mesh`, `create_stretched_mesh` | top level |
| `PyrolysisProblem`, `PyrolysisSolution`, `solve`, `parameter_sweep` | top level |
| `JacobianSpec`, `Structured`, `DirectSolve`/`ApproxSolve`/…, `verify_jacobian_at` | top level |
| `forward_sensitivity`, `adjoint_sensitivity` | top level (stubs; need `SciMLSensitivity, Zygote`) |
| `klu_factorization`, `pardiso_factorization`, `amg_preconditioner` | top level (stubs; need the dep) |
| `ProblemDef`, `Workspace`, `residual!`, `StateLayout`, `to_problem_def`, `build_workspace` | top level |
| `get_profile`, `get_cell_positions`, `interpolate_solution` | `Pyrolysis.Solver` / `Pyrolysis.Internal` |
| `get_surface_temperature`, `get_back_face_temperature` | `Pyrolysis.Discretization` / `Internal` |
| `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`, depletion callbacks | `Pyrolysis.Solver` / `Internal` |
| geometry accessors (`n_cells`, `cell_center`, `total_volume`, …), `n_components` | `Pyrolysis.Geometry` / `Pyrolysis.Materials` / `Internal` |
| conservation trackers, `check_*`, `*_report`, `DiagnosticsCallback` | `Pyrolysis.Diagnostics` / `Internal` |
| trait subtypes (`Eulerian`, `ALE`, `Fickian`, `DarcyFick`, `NoChi`, `WithChi`) | `Pyrolysis.Problem` |
| P1 radiation, h-AMR, Marshak/Reflective radiation BCs | `Pyrolysis.Experimental` (unstable) |
| everything every submodule exports, at once | `using Pyrolysis.Internal` |


---

