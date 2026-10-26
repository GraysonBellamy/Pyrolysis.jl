# 7. Building a Problem: Simulation Modes

By this point you can build the three ingredients of a simulation: a mesh (Chapter 6), a
material with components, reactions and properties (Chapters 3–4), and a
`BoundaryConditionSet` (Chapter 5). This chapter shows how to assemble them into a
`PyrolysisProblem` and how the *physics* the solver integrates is selected. Pyrolysis.jl is
a generalized model with optional physics: one problem can be solved on a fixed or a moving
mesh, with one of three radiation models, with no gas transport, lumped transport or
pressure-driven gas transport, with a fixed or a variable cross section, and with the
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
virgin = SolidComponent(:virgin; ρ = 1190.0, φ=.3, c = 1500.0, k = 0.20, λ = 1e-6)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02,)

decomp = Reaction(:pyrolysis, :virgin => :MMA;
                  A = 8.5e12, E = 188e3, h = 870e3)   # h > 0 = endothermic

pmma = Material(
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6,κ=0.,max_skeletal_occupancy=.001)),
    name = :PMMA,
                components = (virgin, gas),
                reactions  = (decomp,))

# --- Mesh: 1 cm thick, 100 cells ---
mesh = create_uniform_mesh(0.01, 100)

# --- Boundary conditions ---
bc_set = BoundaryConditionSet(
    thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top    = AmbientGasBC(GasReservoir(P=0., T=300.)),   # gas escapes to the specified vacuum
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
declared coefficient `λ = 1e-6` specifies superficial wall transport through
the illustrative porous solid (§7.5).

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
conversion `χ̄`. A lumped gas transport material's pore gas is initialized through `ξ_initial` as well
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
| Gas transport | material model and complete mass conditions | `LumpedTransfer` or `DustyGas`; EOS for every gas-bearing material, including sealed gas | from the problem |
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

![Matching slabs with surface deposition and Beer–Lambert face intensities, absorbed cell power, transmitted radiation and intensity decreasing with inward depth.](../../assets/figures/07_radiation.svg)

**Figure 7. Where applied radiation is absorbed (schematic).** `NO_RADIATION` and `SURFACE_ABSORPTION` both apply `RadiativeFluxBC` at the surface; the former does **not** disable that boundary condition. Beer–Lambert introduces face-intensity unknowns and deposits the entering-minus-leaving radiant power in each cell. The exponential inset illustrates uniform optical thickness $\alpha L=2$; its depth coordinate $L-z$ increases inward, opposite to the common $+z$ direction. Radiation remaining at the back face is transmitted rather than absorbed.

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

## 7.5 Gas transport and pore-state control

`Material(...;gas_transport=LumpedTransfer())` uses the default superficial
partial-pressure law plus Darcy transport. Participating condensed components
require an explicit `λ` [m²/s]; `κ` [m²] defaults to zero. Gas components have
no matrix coefficients.

```julia
matrix = SolidComponent(:matrix;ρ=500.,ρ_intrinsic=1500.,c=1500.,k=.2,
                        λ=1e-5,κ=1e-14)
air = GaseousComponent(:air;M=.02897,c=1006.,k=.026)
material = Material(name=:PorousSlab,components=(matrix,air))
ξ_air = background_gas_concentration(material,:air,(500.,0.),300.)
reservoir = GasReservoir(P=101325.,T=300.,mole_fractions=(air=1.,))
problem = PyrolysisProblem(mesh=create_uniform_mesh(.01,10),material=material,
    bc_set=BoundaryConditionSet(thermal=(top=HeatFluxBC(1000.),bottom=AdiabaticBC()),
        mass=(top=AmbientGasBC(reservoir),bottom=ImpermeableBC())),
    T_initial=z->300.,ξ_initial=[z->500.,z->ξ_air],tspan=(0.,.1))
sol = solve(problem;pressure_tolerance=.1,reltol=1e-7)
pressure = get_pressure(sol)
exchange = get_boundary_exchange(sol,:top)
```

### 7.5.1 Dusty-gas transport

```julia
dusty = DustyGas(tortuosity=2.,D_ref=2e-5,T_ref=300.,P_ref=101325.,
                 pore_diameter=1e-6)
material = Material(name=:DustySlab,components=(matrix,air),gas_transport=dusty)
```

The common molecular-diffusion reference data and tortuosity are explicit.
An explicit pore diameter may be a temperature law. If omitted, the diameter
comes from effective permeability, `sqrt(32τκ_eff/φ)`. A path crossing zero
permeability needs an explicit diameter or a positive, zero-preserving
`capillary_regularization` whose sensitivity is checked for that problem.
Per-component pore-diameter averaging is not used.

### 7.5.2 Domain, endpoints and tolerances

Pore pressure, intrinsic density and fractions are computed for all gas-bearing
materials, including sealed ones. Physical initial and accepted pore volumes
must remain above the trial floor with a resolved error margin. The default
`pressure_tolerance=1.0` Pa is an absolute local error budget shared among gas
species. Proposed-state control includes changing temperature, condensed masses
and geometry; it is not a bound on accumulated global error. Tighten it along
with mesh/time tolerances when qualifying pressure-sensitive predictions.

Matrix transport weights use condensed `γξ/ρ_bulk`; skeletal volume is a
separate sum without `γ`. Retained cells that lose their matrix need a declared
endpoint, for example `LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-5,κ=0.,
max_skeletal_occupancy=1e-3))`. The finite coefficients and supported occupancy
are physical data. Without an endpoint, the model ends at `:matrix_exhausted`.
A nominal porosity for a dense melt must be tested jointly with mesh and
pressure tolerance; it cannot be inferred from the numerical floor.

The single conservative face flux matches total species flow through the two
half paths. Wall-only and Darcy-only layers may therefore conduct in series.
Gas transport remains limited to one dimension when gases can move.

---

## 7.6 Variable cross section (automatic)

![Initial, axially shrunken, laterally shrunken and swollen one-dimensional samples, showing volume as area times thickness and two MLR normalizations.](../../assets/figures/11_shrinkage.svg)

**Figure 11. Axial versus lateral volume change (schematic).** A planar cell has $V=A\Delta z$ and total strain rate $\theta=\theta_L+\theta_A$. The axial balance uses $\theta-\theta_A$, so lateral volume change is counted once. The drawn height stands for cross-sectional area, not a resolved lateral mesh. For the same total outward mass rate, `mass_loss_rate` divides by the current exposed area and `extras.MLR_gauge` by its initial area. Shrinking area makes the former larger than the latter.

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

pmma_shrinking = Material(
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6,κ=0.,max_skeletal_occupancy=.001)),
    name = :PMMA,
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
boundary outflow carries surface enthalpy, while inflow carries its specified
incoming enthalpy with the control-surface heat conversion of §5.4. On a moving mesh under a
correction the mass and enthalpy inventories also cross the faces at the relative velocity
`−c`. The gas temperature-advection term corresponds to ThermaKin Eq. 17;
matching a complete model also requires its other closures (Technical Reference §7.8).

**Reaction heats use the local reaction-temperature convention.** The kernel
evaluates `h` at the cell temperature. Before importing a DSC or calibrated
ThermaKin/Gpyro value, check its sign, units, reaction basis and temperature datum
(Technical Reference §6.9). A heat consistent with this convention already
accounts for the products' enthalpies at that temperature. The energy equation therefore has no separate gas-generation enthalpy term. Such
a term, `−Σ_g Δh_g(T, T_ref)(∇·J_g)`, belongs only with heats referenced to the
`T_ref ≈ 298.15 K` ambient datum (formation-enthalpy style); with DSC-style heats it would
subtract the gas products' sensible enthalpy a second time (≈0.7–1.2 MJ per kg of volatiles
at 900 K, comparable to the heat of pyrolysis itself). If your heats are referenced to the
ambient datum, convert them to the reaction temperature before use (Technical Reference §6).

### 7.7.1 Boundary energy control surface

Boundary species carry signed donor enthalpy. Outflow uses the surface
caloric enthalpy and inflow the specified reservoir enthalpy. A wall thermal
film is converted to the reservoir control surface as described in §5.4,
so incoming energy is counted exactly once. Thermal-film temperature and its
reservoir stream temperature must agree; radiative surroundings are separate.
The energy-state choice below retains the package's sensible-enthalpy storage
convention. A pressure-work/internal-energy-storage redesign is outside this model.

---

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
the reactions (when the material has reactions), the gas transport (lumped or dusty-gas transport),
the in-depth radiation (`BEER_LAMBERT`), the conduction, the thermal boundary conditions,
and, under a correction, the relative fluxes of the inventories across the moving faces;
then, on a moving mesh, the volume-balance rows. You do not choose a path through the code;
you make the choices and `setup` builds the residual that implements them. To see what it
built:

```julia
cp = setup(problem; use_ale = true)
cp.def.setup.gas                            # LumpedTransfer(): the declared gas model
map(k -> nameof(typeof(k)), cp.def.kernels) # the kernels, in order
length(cp.u0)                               # 502 for the problem of §7.1.1
```

`setup` takes the keywords of `solve` that shape the problem (`energy_state`, `reltol`,
`abstol`, `threshold`, `use_ale`, `radiation_model`, `deactivation`, `band`, `correction`,
`lagrangian`, `w_override`, `sources`); `integrate(cp; ...)` then runs it with the
integrator's keywords. Chapter 11 covers this custom-loop path, and Technical Reference §13
and §15 the assembly and the integration.

---

## 7.9 Validating a problem

```julia
using Pyrolysis.Problem: validate_problem, is_valid_problem
validate_problem(problem)
is_valid_problem(problem)
```

Validation checks geometry, reaction indices, boundary tags, initial fields
and the time span. Gas-bearing initial states must have positive physical
pore volume. Setup strengthens this to a resolved margin using the actual
mass, thermal and geometry tolerances, including sealed gas. Named reservoir
fractions must refer to gases in the material and sum to one; non-vacuum
reservoirs require complete temperature and composition data. Incoming
prescribed species fluxes require their temperatures.

Participating condensed components require explicit lumped `λ`, and
transport support changes must have a valid constitutive path. A matrix
endpoint declares its remaining skeletal-occupancy domain. Step acceptance
and post-event checks enforce pore and endpoint admissibility as the state
evolves; `validate=false` does not disable these physical checks.

---

## 7.10 Choosing modes for common scenarios

The table below is a quick recipe lookup. "auto" marks a choice that follows from the
material, not from a `solve` keyword.

| Scenario | `use_ale` | `radiation_model` | gas transport (auto) | cross section (auto) | `deactivation` |
|---|---|---|---|---|---|
| TGA-like thin sample, no recession | `false` | `NO_RADIATION` | lumped | fixed | `false` |
| Cone calorimeter, opaque charring solid | `true` | `SURFACE_ABSORPTION` | lumped or dusty-gas transport | fixed or variable | `false` |
| Cone calorimeter, semi-transparent polymer to burnout | `true` | `BEER_LAMBERT` (set `α`) | lumped | fixed | `true` |
| Intumescent coating (swelling) | `true` | `SURFACE_ABSORPTION` | lumped gas transport (`κ > 0`) | variable (swelling law) | `false` |
| Pressure-driven char, pressure-treated wood | `true` | `SURFACE_ABSORPTION` | lumped gas transport (`κ > 0`) | fixed or variable | `false` |

### 7.10.1 Worked example: TGA-like, minimal mode

A thin sample whose two faces follow a furnace program of 10 K/min, with no radiation and
diffusion-only gas escape through the top, on a fixed mesh:

```julia
using Pyrolysis

virgin = SolidComponent(:virgin; ρ = 1190.0, φ=.3, c = 1500.0, k = 0.20, λ = 1e-6)
gas    = GaseousComponent(:MMA;  M = 0.100, c = 1500.0, k = 0.02,)
decomp = Reaction(:pyrolysis, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)
pmma   = Material(
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6,κ=0.,max_skeletal_occupancy=.001)),
    name = :PMMA, components = (virgin, gas), reactions = (decomp,))

mesh = create_uniform_mesh(1e-4, 10)                  # 0.1 mm, 10 cells

heating = t -> 300.0 + t / 6                          # 10 K/min
bc_set = BoundaryConditionSet(
    thermal = (top = DirichletThermalBC(heating), bottom = DirichletThermalBC(heating)),
    mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = pmma, bc_set = bc_set,
                           tspan = (0.0, 2400.0))

sol = solve(problem)        # fixed mesh, NO_RADIATION, lumped, fixed cross section

# Mass-loss fraction against the furnace temperature
loss = 1 .- sol.extras.total_mass ./ sol.extras.total_mass[1]
T_furnace = heating.(sol.t)
```

### 7.10.2 Worked example: cone calorimeter, opaque, recession and Darcy flow

A charring wood under a cone heater, with surface absorption, a moving mesh and
pressure-driven gas escape (the default lumped model with nonzero permeability), its pore
space filled with air at 1 atm:

```julia
using Pyrolysis

virgin = SolidComponent(:virgin; ρ = 500.0, c = 2000.0, k = 0.15, λ=2e-5, κ = 1e-15,
                        ρ_intrinsic = 1500.0)
char   = SolidComponent(:char;   ρ = 150.0, c = 1500.0, k = 0.10, λ=2e-5, κ = 1e-12,
                        ρ_intrinsic = 1950.0)
gas    = GaseousComponent(:volatiles; M = 0.030, c = 1100.0, k = 0.02,)
air    = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026,)
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
    mass    = (top = AmbientGasBC(GasReservoir(P=101325., T=300., mole_fractions=(air=1.,))),
               bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(mesh = mesh, material = wood, bc_set = bc_set,
                           ξ_initial = [z -> 500.0, z -> 0.0, z -> 0.0, z -> ξ_air],
                           tspan = (0.0, 120.0))

sol = solve(problem;
            use_ale         = true,                 # the mesh follows the solid
            radiation_model = SURFACE_ABSORPTION)   # opaque: absorb at the surface
# Lumped gas transport uses the condensed-phase λ and κ; cross section is fixed.
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
  - *gas transport* (none, lumped or dusty-gas transport) and the *cross section* (fixed or
    variable) follow from the material: its gas components and transfer coefficients, the
    mass conditions, the component permeabilities `κ`, and its `lateral_shrinkage_law`;
  - *deactivation* is a keyword, and on a moving mesh `setup` installs the equidistribution
    correction when a layer can burn out.
- The choices fix the rows of the state (masses and energy per cell; positions and
  velocities on a moving mesh; face intensities under `BEER_LAMBERT`; the conversion history
  and the cross section under a lateral law) and the kernels of the residual.
- There is one energy equation, assembled in flux form with the gas enthalpy carried across
  the faces and DSC-style reaction heats. Outflow carries surface enthalpy; inflow
  carries its specified incoming enthalpy. Thermal-film heat is converted to the
  same reservoir control surface before assembly (§5.4), accounting for incoming
  energy once.
- `validate_problem` (run by `solve` unless `validate = false`) catches structural errors,
  and `setup` refuses unsupported combinations, before integration.

The next chapter, "Running the Solver", covers the rest of the `solve` keyword surface: the
integrator policies and linear solvers, tolerances and step limits, deactivation, callbacks
and the progress bar, output control and the return codes and outcomes.
