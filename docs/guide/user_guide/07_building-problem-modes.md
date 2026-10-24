# 7. Building a Problem: Simulation Modes

By this point you can build the three ingredients of a simulation: a mesh (Chapter 6), a
material with components, reactions and properties (Chapters 3–4), and a
`BoundaryConditionSet` (Chapter 5). This chapter shows how to assemble them into a
`PyrolysisProblem` and how the *physics* the solver integrates is selected. Pyrolysis.jl is
a generalized model with optional physics: one problem can be solved on a fixed or a moving
mesh, with one of three radiation models, with no gas transport, Fickian transport or
pressure-driven Darcy–Fick transport, with a fixed or a variable cross section, and with the
temperature or the enthalpy inventory as the energy state. Some of these choices are `solve`
keywords; the others follow from the material and the boundary conditions. `setup`, which
`solve` calls first, reads them all and builds the state and the residual they need. This
chapter explains every choice, its default, and what `setup` decides from it. The physics and
the derivations live in the Technical Reference; here the focus is on configuring a run
correctly.

---

## 7.1 The `PyrolysisProblem` object

A `PyrolysisProblem` couples a mesh, a material, boundary conditions, initial conditions and
the time span. It is the only thing you pass to `solve`.

### 7.1.1 Constructor

The keyword constructor is:

```julia
PyrolysisProblem(;
    mesh,                       # a Mesh{1} from create_uniform_mesh and its siblings
    material,                   # Material with components and reactions
    bc_set::BoundaryConditionSet,
    T_initial = z -> 300.0,     # initial temperature, z -> T [K]
    ξ_initial = nothing,        # initial concentrations; defaulted when nothing
    tspan,                      # (t_start, t_end) [s], required
    dt_initial = 1e-4,          # the first step [s]
    output_times = Float64[],   # saved times [s]; generated when empty
)
```

The positional form `PyrolysisProblem(mesh, material, bc_set, tspan; kwargs...)` takes the
same keywords. A complete, minimal example (one solid that decomposes into a gas, a constant
heat flux on the exposed surface, an adiabatic and impermeable back):

```julia
using Pyrolysis

# --- Material: one solid that decomposes into a gas ---
virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20, λ = 1e-6)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

decomp = Reaction(:pyrolysis, :virgin => :MMA;
                  A = 8.5e12, E = 188e3, h = 870e3)   # h > 0 = endothermic

pmma = Material(name = :PMMA,
                components = (virgin, gas),
                reactions  = (decomp,))

# --- Mesh: 1 cm thick, 100 cells ---
mesh = create_uniform_mesh(0.01, 100)

# --- Boundary conditions ---
bc_set = BoundaryConditionSet(
    thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top    = DiffusiveMassBC(0.01),   # gas escapes to the atmosphere
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

sol = solve(problem)
```

The mesh knows nothing about the material: the component count comes from the material
alone, and the validator (§7.9) checks that `ξ_initial` has one closure per component. The
small gas-transfer coefficient on the solid (`λ = 1e-6`) keeps Fickian gas transport from
degenerating where the solid holds no gas (§7.5, §7.9).

### 7.1.2 Initial conditions

`T_initial` and `ξ_initial` are functions of the position. `setup` evaluates them at the
centre `x_i` of every cell of the mesh, at its reference positions: `T0[i] = T_initial(x_i)`
and `ξ0[j, i] = ξ_initial[j](x_i)` (`Pyrolysis.Problem.initial_fields(problem, mesh)` returns
both). The initial state then holds the masses `m_{i,j} = ξ0[j, i]·V_i` and the energy state
from `T0[i]`.

- `T_initial :: z -> T` returns the initial temperature [K]. Default `z -> 300.0`.
- `ξ_initial` holds one function per component, each returning the initial mass
  concentration `ξ_j` [kg/m³]. You may pass:
  - `nothing` (default): the constructor fills the first component, **if it is a solid**, at
    its reference density `ρ(300 K)`, and sets every other component to zero;
  - a single function, used for every component;
  - an explicit `Vector` of functions, one per component.

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

The `Problem` submodule exports profile helpers for non-uniform initial states (§5.10.3):
`step_profile`, `linear_profile`, `gaussian_profile` and `tanh_profile`. A sample preheated
toward the exposed surface (`z = L`):

```julia
using Pyrolysis.Problem: linear_profile

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    T_initial = z -> linear_profile(z; z_min = 0.0, z_max = 0.01,
                                       value_min = 300.0, value_max = 400.0),
    tspan = (0.0, 100.0),
)
```

Under a lateral shrinkage law (§7.6) the initial state also holds each cell's conversion
history, `G_i = 0` and `M_i` the cell's initial solid mass, the denominator of the column
conversion `χ̄`. A Darcy–Fick material's pore gas is initialized through `ξ_initial` as well
(§7.5.1). See Technical Reference §12 for the initial and boundary conditions.

### 7.1.3 Time span and output times

- `tspan = (t_start, t_end)` in seconds. The validator requires `t_start < t_end` and
  `t_start ≥ 0`.
- `dt_initial` [s] is the first step the integrator tries (the default of `solve`'s `dt`
  keyword) and, on a moving mesh, the time scale of the node-velocity tolerances (default
  `1e-4` s). It must be positive and no larger than the time span.
- `output_times` is the vector of times at which the solution is saved. Left empty (the
  default), the constructor generates between 11 and 101 evenly spaced times over `tspan`. To
  control the output explicitly:

```julia
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    tspan = (0.0, 120.0),
    output_times = collect(0.0:10.0:120.0),   # every 10 s
)
```

`solve` saves at these times by default (they are its `saveat`), and always saves the state
it ends on; the `saveat` keyword of `solve` overrides them (Chapter 8).

---

## 7.2 What `setup` decides

`solve(problem; kwargs...)` validates the problem, calls `setup(problem; ...)` to build the
problem the integrator runs (a `CoreProblem`), and integrates it (Chapter 8). `setup` reads
every modelling choice once and builds the state layout and the residual's kernels to match;
nothing is dispatched on during the integration. The choices:

| Choice | Selected by | Values | Default |
|---|---|---|---|
| Mesh motion | `use_ale` keyword | fixed (`false`), moving (`true`) | fixed |
| Radiation | `radiation_model` keyword | `NO_RADIATION`, `SURFACE_ABSORPTION`, `BEER_LAMBERT` | `NO_RADIATION` |
| Gas transport | the material and the mass conditions | none, Fickian, Darcy–Fick | from the problem |
| Cross section | the material's `lateral_shrinkage_law` | fixed, variable | fixed |
| Energy state | `energy_state` keyword | `:temperature`, `:enthalpy` | `:temperature` |
| Deactivation | `deactivation` keyword | off, on | off |
| Mesh correction | `correction` keyword, or installed | none, equidistribution | from the problem |

Three choices are made at solve time (mesh motion, radiation, energy state), two follow from
the problem (gas transport, cross section), and deactivation and the correction combine
both: you turn deactivation on, and `setup` installs the correction where the mesh moves and a
layer can burn out (§7.3.3). The coordinate system and the cross-section measure are the
mesh's (§6.4). A combination the implementation does not support is refused with an
`ArgumentError` before anything is integrated (§7.9.1). §7.8 shows how the choices shape the
state and the residual, and how to inspect what `setup` built.

---

## 7.3 Mesh motion: fixed or moving (`use_ale`)

### 7.3.1 What it controls

- **Fixed** (`use_ale = false`, default): every node is pinned. The node positions are not
  part of the state, the geometry is constant, and a cell carries its component masses and
  its energy state.
- **Moving** (`use_ale = true`): the node positions and velocities are part of the state. The
  back node is pinned; every other node moves with `w = v_s + c`, the solid velocity that the
  cells' reaction strain rates `θ` give through one volume-balance row per cell, plus the
  correction `c` (zero unless one is installed, §7.3.3). This lets the domain shrink as the
  solid is consumed (charring, ablation) or expand (swelling, intumescence) without a
  moving-boundary equation. Every component, pore gas included, rides the mesh; the gas
  moves relative to it through its transport flux.

```julia
# Fixed mesh (default)
sol = solve(problem)

# Moving mesh that tracks recession and expansion
sol = solve(problem; use_ale = true)
```

### 7.3.2 When to use a moving mesh

Use a moving mesh whenever the through-thickness extent of the material changes
appreciably: charring wood, ablating or gasifying polymers, intumescent coatings, any case in
which surface recession matters for the surface temperature and the mass-loss rate, and every
solve that burns cells out at the surface. For a TGA-like or thin-film case whose geometry
change does not matter, the fixed mesh is cheaper: it has no node unknowns (§7.8).

On a moving mesh the cell volumes follow from the node positions at every evaluation, so the
geometric conservation law holds by construction; the solve ends with outcome
`:min_thickness` when the thickness falls below `min_thickness` (default `1e-6` m). See
Chapter 6 for the node positions in the solution (§6.7) and for meshes stretched toward the
surface (§6.2.2), which pair well with it, and Technical Reference §11 for the moving-mesh
formulation.

### 7.3.3 Deactivation and the correction

`deactivation = true` deactivates, during the solve, every cell whose condensed mass falls
below `threshold` (default `1e-2`) of its initial condensed mass, from the exposed side, in
an event located inside the step; its inventories go to the live cell beneath it. Chapter 10
describes the burnout path, its parameters and the envelope in which it is verified.

On a moving mesh with deactivation on, `setup` installs the **equidistribution correction**
(`EquidistributionCorrection()`, its layer time set at each layer's first activation,
`ActivationTime()`) when some layer can vanish: its components, converted completely, leave
less condensed mass than `threshold` of what it started with (`Pyrolysis.Closures.residues`
gives the fraction each component leaves). The correction moves the mesh relative to the
solid, thinning a vanishing layer's cells together so that the layer leaves in one event;
the near-threshold band then defaults to `0.1`. Otherwise (a fixed mesh, deactivation off,
or a charring layer whose residue stays above the threshold) the solve has no correction and
a moving mesh moves with the solid.

| Keyword | Effect |
|---|---|
| `correction = NoCorrection()` | Keeps the mesh with the solid whatever `setup` would install. |
| `correction = EquidistributionCorrection(; ...)` | Configures the correction (Chapter 10); needs `use_ale = true`. |
| `lagrangian = [nodes...]` | Nodes that always move with the solid. Under the equidistribution correction every interior node between cells that start with different condensed components must be among them, and by default they are exactly those nodes. |
| `band` | The near-threshold band of a moving mesh, in `[threshold, 1)`. |

---

## 7.4 Radiation model (`radiation_model`)

Where the incident radiation is absorbed is selected by the `radiation_model` keyword, which
takes a `RadiationModel` value:

| Value | Meaning | Unknowns added |
|---|---|---|
| `NO_RADIATION` (default) | No volumetric absorption. A `RadiativeFluxBC` is absorbed at the exposed surface. | none |
| `SURFACE_ABSORPTION` | The incident flux is absorbed at the exposed surface, through the surface energy balance (opaque material). Component `α` values are not read. | none |
| `BEER_LAMBERT` | In-depth exponential absorption, `I(z) = I_e exp(−∫ α_eff dz)`, with `α_eff = Σ_j ξ_j α_j` from the components' mass absorption coefficients `α` [m²/kg]. | one intensity per face, algebraic |

```julia
# Opaque material: absorb the cone flux at the surface
sol = solve(problem; radiation_model = SURFACE_ABSORPTION)

# Semi-transparent material: absorb in depth (needs α on the components)
sol = solve(problem; radiation_model = BEER_LAMBERT)
```

### 7.4.1 How radiation enters the energy equation

The *incident* radiative flux comes from a `RadiativeFluxBC` in the thermal entry of the
`:top` tag (a cone heater, for instance). The `radiation_model` keyword decides where that
flux is deposited:

- With `SURFACE_ABSORPTION`, the absorbed flux `a·q_inc` enters the surface energy balance at
  the exposed face, raising `T_s`; `a` is the condition's own `absorptivity` when it gives
  one, otherwise the exposed cell's effective absorptivity (Kirchhoff's law, the emissivity,
  unless a component sets its own `absorptivity`). Nothing is deposited in depth.
- With `NO_RADIATION`, the same: a `RadiativeFluxBC` is absorbed at the surface exactly as
  under `SURFACE_ABSORPTION`, and no volumetric source exists. The two values describe the
  same model; name the one that says what you mean.
- With `BEER_LAMBERT`, the exposed face's intensity is `I_e = a·q_inc(t)` and every other
  face's intensity is the one above it attenuated across the cell, `exp(−α_eff Δz)`; each
  cell absorbs the difference between the intensity entering and leaving it, the exact cell
  average of the attenuated beam, and the intensity leaving the back face is transmitted out
  of the domain (`sol.extras.radiation_transmitted`). The face intensities are algebraic rows
  of the state, so the Jacobian stays banded. The `RadiativeFluxBC` then contributes nothing
  to the surface balance.

Under every model `RadiativeBC` (re-radiation, `Fεσ(T_∞⁴ − T_s⁴)`) and `ConvectiveBC` act at
the boundary, and a time-varying incident flux is evaluated at the time of every residual
evaluation, so a ramped or stepped cone flux is tracked exactly.

### 7.4.2 Restrictions and choosing a model

- `BEER_LAMBERT` needs planar coordinates, and the incident flux enters through the `:top`
  faces only: a `RadiativeFluxBC` on another tag is refused with an `ArgumentError`.
- `BEER_LAMBERT` with `α = 0` on every component absorbs nothing in depth: the whole absorbed
  share of the incident flux passes through and leaves through the back face. No warning is
  given; check `sol.extras.volumetric_absorption` and `sol.extras.radiation_transmitted`.
- `SURFACE_ABSORPTION` ignores the components' `α`, without a warning.

Choose the model by the optical thickness `τ = α_eff·L`: `SURFACE_ABSORPTION` for opaque
chars (`τ ≫ 1`), `BEER_LAMBERT` for semi-transparent polymers and foams, and `NO_RADIATION`
for transparent samples or when the heating is purely conductive or convective. See
Technical Reference §8 for the models.

---

## 7.5 Gas transport: none, Fickian or Darcy–Fick (automatic)

The gas-transport model is not a `solve` keyword: `setup` infers it from the material and
the mass boundary conditions.

| Model | Selected when | Gas flux |
|---|---|---|
| none | the material has no gas component, or no component has a gas-transfer coefficient `λ` and every mass condition is `ImpermeableBC` | no gas moves |
| Fickian | gas can move and no condensed component is permeable (every `κ = 0`) | diffusive: `J = −(λ_f/T_f)(ξ_nb T_nb − ξ_o T_o)/d` |
| Darcy–Fick | some solid or liquid component has `κ > 0` (`has_darcy_flow(material)`) | diffusive plus advective: `J = J_diff + v ξ_up/φ_f`, with the superficial velocity `v = −(κ/μ)∇P` |

In the Fickian flux `λ_f` is the face mean of the cells' gas-transfer coefficients and `d`
the distance between the cell centres; in the Darcy–Fick flux `ξ_up` is the concentration
upwind of `v` and `φ_f` the face porosity, so that the advected quantity is the intrinsic
gas density `ξ/φ` (`ξ` is per bulk volume, `v` superficial). The pressure `P` is the
ideal-gas pressure of the pore gas, and it and the Darcy velocity are evaluated only under
Darcy–Fick transport, so a Fickian material carries no pressure computation. See Technical
Reference §9 for the transport model.

Under Fickian transport the cell coefficient is the volume-fraction mix of the component
values, so a condensed component with `λ = 0` stops all transport in a cell that holds no gas.
The validator warns about it (§7.9); give the condensed components a small `λ` (for example
`1e-6` m²/s).

### 7.5.1 Enabling Darcy–Fick

Put a permeability on the relevant solid and liquid components:

```julia
# Diffusion only (Fickian): no κ, all gas transport is diffusive
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, λ = 1e-6)

# Darcy–Fick: a permeability selects it
virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, κ = 1e-15,
                        ρ_intrinsic = 1500.0)  # skeletal density: a real pore volume
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, κ = 1e-12,
                        ρ_intrinsic = 1950.0)
gas    = GaseousComponent(:volatiles; M = 0.030, c = 1100.0, k = 0.02, λ = 1e-5)
air    = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)

decomp = Reaction(:pyrolysis, :virgin => (:char, :volatiles);
                  A = 1e10, E = 150e3, h = 300e3, yields = (0.25, 0.75))

wood = Material(name = :CharringWood,
                components = (virgin, char, gas, air),
                reactions  = (decomp,),
                mixing = MaterialMixing(k = MixingSpec(PARALLEL, 0.5),
                                        κ = MixingSpec(WEIGHTED, 0.5)))
```

> **Initialize the pore gas.** A Darcy–Fick problem needs (i) skeletal densities on the
> condensed components (`ρ_intrinsic`, or the `φ` keyword, §3.2.4) so that a real pore volume
> exists, and (ii) an inert background air component initialized at ambient pressure with
> `background_gas_concentration(material, :air, ξ_condensed₀, T₀)`; otherwise the pore space
> starts as a vacuum against the 1 atm boundary and early product gas is spuriously advected
> inward. A permeable cell with no pore space is refused by the validator (§7.9). See
> Technical Reference §9.5.4 and `examples/08_darcy_fick_demo.jl`.

```julia
# Pore air at 1 atm in the virgin solid at 300 K
ξ_air = background_gas_concentration(wood, :air, (500.0, 0.0, 0.0, 0.0), 300.0)
ξ_initial = [z -> 500.0, z -> 0.0, z -> 0.0, z -> ξ_air]
```

The effective permeability is computed from the component values with the `κ` mixing rule
(§3.5). As virgin converts to char, `κ` rises and the bulk gas flow gains weight. Gas
components carry no permeability (the gas flows *through* the matrix); only solids and
liquids contribute.

### 7.5.2 Pressure boundary conditions

Darcy flow is driven by the gradient of the ideal-gas pore pressure. To control the pressure
at a boundary, attach a pressure condition to the **mass** entry of that tag with the `+`
operator. Pressure conditions are not a separate field of `BoundaryConditionSet`; they ride
along with the mass condition:

```julia
# Top surface: gas diffuses out AND is held at atmospheric pressure
top_mass = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

bc_set = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(75e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = top_mass, bottom = ImpermeableBC()),
)
```

The pressure conditions are `DirichletPressureBC(P)` (fixed pressure),
`NeumannPressureBC(∂P/∂n)` (prescribed gradient) and `ConvectivePressureBC` (a resistive
outflow law). Each sets the outward Darcy velocity of its boundary face, which vents gas
advectively against an ambient that carries none of the tracked gases; they contribute no
diffusive film flux of their own: that comes from the mass condition they are combined with,
and the two add. Under Fickian transport a pressure condition drives nothing. See §5.5,
Technical Reference §9.5.3 for the venting flux and §12.5 for the pressure conditions.

> The shipped Darcy–Fick examples (`08_darcy_fick_demo.jl`, `07_douglas_fir_cone.jl`,
> `08_pressure_treated_wood_cone.jl`) run on a moving mesh, since pressure-driven outflow
> accompanies recession.

---

## 7.6 Variable cross section (automatic)

The cross section also follows from the material. If the material's
`lateral_shrinkage_law` is `nothing` (the default), the column's cross section is fixed at
the mesh's `A₀`. If you supply a law, a function of the column conversion `χ̄`, the area
varies as `A(t) = A₀ · law(χ̄)` (§6.4.3).

### 7.6.1 Defining a shrinkage or swelling law

```julia
# Linear lateral shrinkage: the area shrinks to 80% as pyrolysis completes
shrink_law = χ̄ -> 1.0 - 0.2 * χ̄

# Linear lateral swelling (intumescence): the area grows by 50%
swell_law  = χ̄ -> 1.0 + 0.5 * χ̄

pmma_shrinking = Material(name = :PMMA,
                          components = (virgin, gas),
                          reactions  = (decomp,),
                          lateral_shrinkage_law = shrink_law)
```

The column conversion is `χ̄ = Σ_i G_i / Σ_i M_i`, with `G_i` the cumulative dry pyrolysis
gas cell `i` has produced and `M_i` its initial solid mass; only "dry" gases contribute, every
gas component except water vapour (a component whose name contains `H2O` or `water`, in any
case). The law must be smooth (its derivative is taken by automatic differentiation, so a
piecewise law with a kink gives an ill-defined lateral rate), positive and finite on
`[−0.001, 1.001]`, and equal to one at `χ̄ = 0`; `setup` checks the last two and refuses a law
that fails them. Without a law, no cross-section unknowns enter the state.

### 7.6.2 Interaction with the moving mesh

A variable cross section requires `use_ale = true` and planar coordinates: `setup` refuses a
lateral shrinkage law on a fixed mesh with an `ArgumentError` naming the fix. With the law,
`setup` adds to the state each cell's conversion-history pair `(G_i, M_i)`, carried with the
solid, and two algebraic rows for the column: the area `A = A₀·law(χ̄)` and its lateral strain
rate `θ_A = d ln law(χ̄)/dt`. Each cell's volume balance uses the axial strain `θ_i − θ_A`, so
the lateral and the axial channels share the volume change without counting it twice.

```julia
problem = PyrolysisProblem(mesh = mesh, material = pmma_shrinking, bc_set = bc_set,
                           tspan = (0.0, 100.0))

sol = solve(problem; use_ale = true)

sol.extras.cross_section_area_history     # A(t) at the saved times [m²]
sol.extras.MLR_gauge                      # mass-loss rate per unit initial area
```

See Technical Reference §10 for volume change, swelling and the lateral shrinkage law, and
§11 for how `θ_A` enters the mesh motion.

---

## 7.7 The energy equation and the energy state

The energy equation is the same in every mode; the choices above change how its terms are
evaluated, not which terms appear:

```
ρc_p ∂T/∂t = −∇·q − Σ_g J_g·∇h_g + Q_rad + Q_rxn
```

Here `ρc_p = Σ_j ξ_j c_p,j` sums every component, pore gas included; `q = −k∇T` is the
conducted flux; `J_g` are the gas mass fluxes relative to the mesh and `h_g` the gas sensible
enthalpies, so the second term is the sensible heat the moving gas exchanges with the matrix
it crosses (`−Σ_g c_p,g J_g·∇T`); `Q_rad` is the in-depth radiation source (`BEER_LAMBERT`
only); and `Q_rxn = −Σ_r h_r r_r` is the reaction heat (`h_r > 0` endothermic).

The solver assembles it in flux form. Every face carries the conducted heat and the sensible
enthalpy `Σ_j J_j h_j(T_f)` of the gas crossing it at the face temperature; every cell's row
subtracts the enthalpy `Σ_j h_j(T_i)·Ṁ_{i,j}` that its net gas inflow brings at the cell's own
temperature, which leaves exactly the exchange term above. The gas that leaves through a
boundary face carries its enthalpy at the surface temperature. On a moving mesh under a
correction the mass and enthalpy inventories also cross the faces at the relative velocity
`−c`. This is the ThermaKin convention (Technical Reference §3 and §7).

**Reaction heats are DSC-style.** A heat `h` measured by DSC, or calibrated in ThermaKin or
Gpyro, is quoted at the reaction temperature, so it already contains the gas products'
sensible enthalpy at that temperature; the reaction kernel evaluates it at the cell's
temperature. The energy equation therefore has no separate gas-generation enthalpy term. Such
a term, `−Σ_g Δh_g(T, T_ref)(∇·J_g)`, belongs only with heats referenced to the
`T_ref ≈ 298.15 K` ambient datum (formation-enthalpy style); with DSC-style heats it would
subtract the gas products' sensible enthalpy a second time (≈0.7–1.2 MJ per kg of volatiles
at 900 K, comparable to the heat of pyrolysis itself). If your heats are referenced to the
ambient datum, convert them to the reaction temperature before use (Technical Reference §6).

### 7.7.1 No surface gas-enthalpy term

The surface energy balance contains no gas-enthalpy term: the flux the boundary condition
applies equals the conducted flux (§5.4), and the enthalpy of the escaping gas leaves through
the face's energy flux at the surface temperature, as above. A surface term
`Σ_g ṁ_g Δh_g(T_s, T_∞)` would remove that enthalpy a second time, a spurious sink that
under-predicts the surface temperature and the mass-loss rate. Use a plain convective and
radiative boundary condition; the escaping gas influences the boundary balance only through
the blowing model of a `ConvectiveBC` (§5.2.3). See Technical Reference §7 and §16 for the
energy accounting.

### 7.7.2 The energy state (`energy_state`)

The integrator carries one energy unknown per cell, chosen with the `energy_state` keyword:

| `energy_state` | Unknown | Row | Conservation |
|---|---|---|---|
| `:temperature` (default) | the cell temperature `T_i` [K] | `(Σ_j m_{i,j} c_p,j) dT_i/dt = E_i − Σ_j h_j(T_i) Ṁ_{i,j}` | energy to the integrator's order |
| `:enthalpy` | the sensible-enthalpy inventory `H_i = Σ_j m_{i,j} h_j(T_i)` [J] | `dH_i/dt = E_i + Σ_j ω_{i,j} V_i h_j(T_i)` | the summed inventory, a linear invariant, to round-off |

`E_i` is the cell's energy sources plus the energy its faces bring in, `Ṁ_{i,j}` the face
part of its mass rate and `ω_{i,j}` its reaction mass rates. In the inventory form the
temperature is recovered from `H_i` and the masses by a regularized inversion of the caloric
laws at every evaluation. Both forms run on every path (fixed and moving meshes,
deactivation), and the burnout envelope of Chapter 10 holds for either. On the problem of
§7.1.1 at the default tolerances the two forms' temperatures differ by at most 0.016 K over
the 100 s (measured by solving it in both forms). Any other value is refused with an
`ArgumentError`.

```julia
sol = solve(problem; energy_state = :enthalpy)
```

---

## 7.8 How the choices shape the state and the residual

The choices of §7.2 determine which unknowns the state holds and which kernels the residual
runs. The rows of the state:

| Rows | Present when | Kind |
|---|---|---|
| component masses `m_{i,1} … m_{i,NC}`, per cell | always | differential |
| energy state (`T_i` or `H_i`), per cell | always | differential |
| node position `x_k`, per node | moving mesh | differential |
| node velocity `w_k`, per node | moving mesh | algebraic (the volume-balance rows) |
| face intensity `I_f`, per face | `BEER_LAMBERT` | algebraic (the attenuation rows) |
| conversion history `G_i`, `M_i`, per cell | lateral shrinkage law | differential |
| cross section `A` and lateral rate `θ_A`, one each | lateral shrinkage law | algebraic |

With `n` cells and `NC` components the state has `n(NC + 1)` rows, plus `2(n + 1)` on a
moving mesh, `n + 1` under `BEER_LAMBERT` and `2n + 2` under a lateral shrinkage law. The
system is a mass-matrix differential-algebraic system `M du/dt = f(u, t)`, the algebraic rows
having no time derivative, and `setup` makes the algebraic unknowns consistent with their rows
before the solve starts (Chapter 8, Technical Reference §15). Because the masses are the
states, mass is conserved exactly.

The residual runs its kernels in a fixed order, each present only when the problem needs it:
the reactions (when the material has reactions), the gas transport (Fickian or Darcy–Fick),
the in-depth radiation (`BEER_LAMBERT`), the conduction, the thermal boundary conditions,
and, under a correction, the relative fluxes of the inventories across the moving faces;
then, on a moving mesh, the volume-balance rows. You do not choose a path through the code;
you make the choices and `setup` builds the residual that implements them. To see what it
built:

```julia
cp = setup(problem; use_ale = true)
cp.def.setup.gas                            # FickTransport(): the gas transport chosen
map(k -> nameof(typeof(k)), cp.def.kernels) # the kernels, in order
length(cp.u0)                               # 502 for the problem of §7.1.1
```

`setup` takes the keywords of `solve` that shape the problem (`energy_state`, `reltol`,
`abstol`, `threshold`, `use_ale`, `radiation_model`, `deactivation`, `band`, `correction`,
`lagrangian`, `w_override`, `sources`); `integrate(cp; ...)` then runs it with the
integrator's keywords. Chapter 11 covers this custom-loop path, and Technical Reference §13
and §15 the assembly and the integration.

---

## 7.9 Validating a problem (`validate_problem`)

Before solving, run the validator to catch structural mistakes early:

```julia
using Pyrolysis.Problem: validate_problem, is_valid_problem

validate_problem(problem)                 # throws on any error, returns true if OK
validate_problem(problem; strict = true)  # adds the checks below
is_valid_problem(problem)                 # a Bool, never throws
```

`solve` calls `validate_problem(problem; strict = false)` unless you pass `validate = false`.
The default checks are:

- the mesh has one finite position per node;
- the material has at least one component, and every reaction references valid component
  indices;
- the mesh carries the boundary tags `:top` (the exposed surface) and `:bottom`, every
  boundary tag of the mesh has a thermal condition, and a per-species `ξ_inf` tuple on a
  `ConvectiveMassBC` has one entry per component;
- the time span is valid (`t_start < t_end`, `t_start ≥ 0`, a positive `dt_initial` no larger
  than the span);
- the number of initial concentration functions equals the component count;
- a Darcy–Fick material does not start with (near-)zero porosity, the pore fraction
  `1 − skeletal_volume_fraction` evaluated at the cell centres, in a cell whose effective
  permeability is nonzero (§7.5.1).

It also warns, without failing, when a gas-producing Fickian material has a condensed
component with gas-transfer coefficient `λ = 0`: gas transport then vanishes where the gas
fraction vanishes, the gas field forms a sharp front, and the integrator rejects step
after step on it. Give such components a small `λ`, for example `1e-6` m²/s. Errors are
collected and thrown together in one message.

`strict = true` additionally verifies the reactions' stoichiometric mass balance (product
yields summing to one), that the initial-condition functions evaluate at the first cell
centre, and that every cell volume is positive.

### 7.9.1 What `setup` refuses

Some combinations pass validation but are not implemented; `setup` refuses them with an
`ArgumentError` before anything is integrated, whether or not `solve` validated:

- a lateral shrinkage law without `use_ale = true`, on a non-planar mesh, or a law that is
  not positive and finite on `[−0.001, 1.001]` or not one at zero conversion;
- `BEER_LAMBERT` on a cylindrical or spherical mesh, or with a `RadiativeFluxBC` on a tag
  other than `:top`;
- a `correction` other than `NoCorrection()` without `use_ale = true`;
- a `w_override` (a prescribed mesh velocity, for verification) without `use_ale = true`,
  with a reacting material, with `deactivation`, or together with a correction;
- under the equidistribution correction, a `lagrangian` list that leaves out an interior node
  between cells that start with different condensed components;
- a `band` outside `[threshold, 1)`, and an `energy_state` other than `:temperature` or
  `:enthalpy`.

`solve` itself refuses an unknown keyword and a non-positive `min_thickness`, also with an
`ArgumentError`.

---

## 7.10 Choosing modes for common scenarios

The table below is a quick recipe lookup. "auto" marks a choice that follows from the
material, not from a `solve` keyword.

| Scenario | `use_ale` | `radiation_model` | gas transport (auto) | cross section (auto) | `deactivation` |
|---|---|---|---|---|---|
| TGA-like thin sample, no recession | `false` | `NO_RADIATION` | Fickian | fixed | `false` |
| Cone calorimeter, opaque charring solid | `true` | `SURFACE_ABSORPTION` | Fickian or Darcy–Fick | fixed or variable | `false` |
| Cone calorimeter, semi-transparent polymer to burnout | `true` | `BEER_LAMBERT` (set `α`) | Fickian | fixed | `true` |
| Intumescent coating (swelling) | `true` | `SURFACE_ABSORPTION` | Darcy–Fick (`κ > 0`) | variable (swelling law) | `false` |
| Pressure-driven char, pressure-treated wood | `true` | `SURFACE_ABSORPTION` | Darcy–Fick (`κ > 0`) | fixed or variable | `false` |

### 7.10.1 Worked example: TGA-like, minimal mode

A thin sample whose two faces follow a furnace program of 10 K/min, with no radiation and
diffusion-only gas escape through the top, on a fixed mesh:

```julia
using Pyrolysis

virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20, λ = 1e-6)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
decomp = Reaction(:pyrolysis, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)
pmma   = Material(name = :PMMA, components = (virgin, gas), reactions = (decomp,))

mesh = create_uniform_mesh(1e-4, 10)                  # 0.1 mm, 10 cells

heating = t -> 300.0 + t / 6                          # 10 K/min
bc_set = BoundaryConditionSet(
    thermal = (top = DirichletThermalBC(heating), bottom = DirichletThermalBC(heating)),
    mass    = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = pmma, bc_set = bc_set,
                           tspan = (0.0, 2400.0))

sol = solve(problem)        # fixed mesh, NO_RADIATION, Fickian, fixed cross section

# Mass-loss fraction against the furnace temperature
loss = 1 .- sol.extras.total_mass ./ sol.extras.total_mass[1]
T_furnace = heating.(sol.t)
```

### 7.10.2 Worked example: cone calorimeter, opaque, recession and Darcy flow

A charring wood under a cone heater, with surface absorption, a moving mesh and
pressure-driven gas escape (Darcy–Fick selected by the component permeabilities), its pore
space filled with air at 1 atm:

```julia
using Pyrolysis

virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, κ = 1e-15,
                        ρ_intrinsic = 1500.0)
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, κ = 1e-12,
                        ρ_intrinsic = 1950.0)
gas    = GaseousComponent(:volatiles; M = 0.030, c = 1100.0, k = 0.02, λ = 1e-5)
air    = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
decomp = Reaction(:pyrolysis, :virgin => (:char, :volatiles);
                  A = 1e10, E = 150e3, h = 300e3, yields = (0.25, 0.75))

wood = Material(name = :CharringWood,
                components = (virgin, char, gas, air),
                reactions  = (decomp,),
                mixing = MaterialMixing(k = MixingSpec(PARALLEL, 0.5),
                                        κ = MixingSpec(WEIGHTED, 0.5)))

# Pore air at 1 atm in the virgin solid
ξ_air = background_gas_concentration(wood, :air, (500.0, 0.0, 0.0, 0.0), 300.0)

mesh = create_uniform_mesh(0.02, 80)                  # 2 cm, 80 cells

# Top: incident cone flux (absorbed at the surface), re-radiation, convection.
# Mass: diffusive escape and atmospheric pressure (drives the Darcy outflow).
bc_set = BoundaryConditionSet(
    thermal = (top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.9) +
                     RadiativeBC(ε = 0.9, T_inf = 300.0) +
                     ConvectiveBC(h = 10.0, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
               bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = wood, bc_set = bc_set,
                           ξ_initial = [z -> 500.0, z -> 0.0, z -> 0.0, z -> ξ_air],
                           tspan = (0.0, 120.0))

sol = solve(problem;
            use_ale         = true,                 # the mesh follows the solid
            radiation_model = SURFACE_ABSORPTION)   # opaque: absorb at the surface
# gas transport = Darcy–Fick (from κ); cross section fixed (no law)
```

To add lateral intumescent swelling, give `wood` a `lateral_shrinkage_law`
(`χ̄ -> 1.0 + 0.5 * χ̄`, say): the cross section becomes variable and the `solve` call stays
the same, since it already moves the mesh. This char keeps a quarter of the virgin mass, more
than the deactivation threshold, so a burnout solve of it would get no correction (§7.3.3).

---

## 7.11 Summary

- A `PyrolysisProblem` bundles the mesh, the material, `bc_set`, the initial conditions and
  `tspan`/`output_times`. The mesh carries no component count; the initial conditions,
  evaluated at the cell centres by `setup`, must give one function per component.
- The active physics is a set of choices that `setup` reads once:
  - *mesh motion* (`use_ale`, default `false`), *radiation* (`radiation_model`, default
    `NO_RADIATION`) and the *energy state* (`energy_state`, default `:temperature`) are
    chosen at solve time;
  - *gas transport* (none, Fickian or Darcy–Fick) and the *cross section* (fixed or
    variable) follow from the material: its gas components and transfer coefficients, the
    mass conditions, the component permeabilities `κ`, and its `lateral_shrinkage_law`;
  - *deactivation* is a keyword, and on a moving mesh `setup` installs the equidistribution
    correction when a layer can burn out.
- The choices fix the rows of the state (masses and energy per cell; positions and
  velocities on a moving mesh; face intensities under `BEER_LAMBERT`; the conversion history
  and the cross section under a lateral law) and the kernels of the residual.
- There is one energy equation, assembled in flux form with the gas enthalpy carried across
  the faces and DSC-style reaction heats; the surface balance has no gas-enthalpy term, which
  would count the escaping gas's enthalpy twice.
- `validate_problem` (run by `solve` unless `validate = false`) catches structural errors,
  and `setup` refuses unsupported combinations, before integration.

The next chapter, "Running the Solver", covers the rest of the `solve` keyword surface: the
integrator policies and linear solvers, tolerances and step limits, deactivation, callbacks
and the progress bar, output control and the return codes and outcomes.
