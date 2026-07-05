# Pyrolysis.jl Examples

This directory contains example scripts demonstrating the Pyrolysis.jl package.

## Setup

This directory has its own `Project.toml` environment that includes visualization dependencies (Plots) alongside Pyrolysis. This keeps the main package lightweight — users installing Pyrolysis via `Pkg.add()` won't require Plots unless they want to run examples.

To run examples, activate the examples environment:

```julia
cd("examples")
using Pkg; Pkg.activate(".")
Pkg.instantiate()  # first time only

include("06_cone_calorimeter.jl")
```

Or from the shell:
```bash
cd examples
julia --project=. -e 'include("06_cone_calorimeter.jl")'
```

## Available Examples

### 06_cone_calorimeter.jl - Cone Calorimeter Simulation of PMMA

Simulates the pyrolysis of PMMA under standard cone calorimeter conditions:
- **Heat flux**: 50 kW/m² external radiant heat flux
- **Top surface**: Convective and re-radiation losses, mass release
- **Bottom surface**: Adiabatic, impermeable (insulated sample holder)
- **Sample**: 6 mm thick PMMA

This is a fundamental 1D pyrolysis simulation demonstrating:
- Material definition with temperature-dependent properties
- Single-step Arrhenius decomposition kinetics
- Combined thermal boundary conditions (radiation + convection)
- Gas generation and release

**Run with:**
```julia
using Pyrolysis
include("examples/06_cone_calorimeter.jl")
```

### 07_douglas_fir_cone.jl - Douglas Fir Cone Calorimeter with ALE

Comprehensive simulation of wood pyrolysis with advanced features:
- **Material**: Douglas Fir with multi-step decomposition kinetics
- **Radiation**: Surface absorption model for opaque materials
- **Transport**: Darcy-Fick gas transport with species-specific scaling
- **Mesh**: ALE mesh motion tracking material shrinkage
- **Performance**: Reduced Jacobian for ~8× speedup over full ALE Jacobian (see below)

Demonstrates boundary temperature diagnostics showing the difference between
the solved surface temperature and cell-center temperature.

**Run with:**
```julia
using Pyrolysis
include("examples/07_douglas_fir_cone.jl")
```

## Reduced Jacobian for ALE Simulations

ALE (moving mesh) simulations include a mesh velocity equation where each node's
velocity depends on the reaction rates of all cells below it (cumulative sum).
This creates a dense lower-triangular block in the Jacobian that dominates the
sparse factorization cost (69% of non-zeros at 240 cells, O(N^2.3) scaling).

The default `Structured(strategy = DirectSolve())` includes the structured
coupling exactly. Pass `Structured(strategy = ApproxSolve())` to opt into the
Schur-style approximation that drops the dense coupling — Newton converges
all terms to tolerance using the approximate search direction; legacy
benchmarks measured **~7.8× speedup** at 240 cells with identical solutions.

```julia
solution = solve(problem;
    use_ale = true,
    jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem),
    integrator = KenCarp4(linsolve = SparspakFactorization()),
)
```

## Boundary Temperature Behavior

In cell-centered finite volume methods, temperature is computed at cell centers.
For boundary conditions with strong temperature dependence (especially radiation
with T⁴ scaling), Pyrolysis.jl solves for the true surface temperature satisfying
energy balance:

```
q_applied(T_s) = q_conducted(T_s, T_cell)
```

This ensures:
- Energy balance at boundaries to machine precision
- Correct temperature for radiative flux calculation
- Accurate surface temperature output for comparison with experiments

The solved boundary state is stored in `mesh.boundary_states` and can be accessed via:
```julia
using Pyrolysis.Discretization: get_surface_temperature, check_boundary_energy_balance

T_surface = get_surface_temperature(mesh)
balance = check_boundary_energy_balance(mesh)  # Dict{Int, Float64}
```

See `docs/boundary_temperature_enhancement.md` for detailed physics explanation.

## Planned Examples

1. **01_heat_conduction.jl** - Simple verification case (heat conduction only)
2. **02_single_reaction.jl** - Basic pyrolysis with single decomposition reaction
3. **03_multi_step.jl** - Complex multi-step kinetics
4. **04_charring_ablation.jl** - Charring material with surface regression
5. **05_intumescent.jl** - Swelling/intumescent material
7. **07_adaptive_showcase.jl** - AMR demonstration

## Running Examples

```julia
using Pyrolysis
include("examples/01_heat_conduction.jl")
```
