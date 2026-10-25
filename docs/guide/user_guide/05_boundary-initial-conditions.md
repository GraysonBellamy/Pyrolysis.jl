# 5. Boundary and Initial Conditions

Boundary conditions (BCs) tell the solver how heat and gas enter and leave the
1-D column; initial conditions set the temperature and composition at `t = 0`.
This chapter shows how to build a `BoundaryConditionSet`, the complete catalogue
of thermal, mass, and pressure BCs with every field and default,
how the surface energy balance is solved, and how to specify initial
temperature and concentration profiles. The physics and derivations behind each
BC live in the Technical Reference §12 (and the radiation models in §8); this
chapter is about *how to specify them correctly*.

> **Coordinate and sign reminders** (see Technical Reference §2). The axial
> coordinate is `z`. `z = 0` is the **bottom / substrate** (tag `:bottom`);
> `z = L` is the **top / exposed surface** (tag `:top`). Heat normally enters from the top. For *all* BCs in this package the
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
    mass    = (top    = AmbientGasBC(GasReservoir(P=0., T=300.)),
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
- **It is type-stable.** Each tag's condition is a concrete type, so the
  residual specializes on it statically.
- **The tags are the mesh's.** A line mesh has `:top` and `:bottom`; validation
  requires a thermal condition for every tag of the mesh (the merged defaults
  supply them), and `setup` refuses a set that names a tag the mesh does not have
  with an `ArgumentError`.

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

> **Note on mesh motion and radiation.** The `BoundaryConditionSet` has *only*
> `thermal` and `mass` fields. Mesh motion needs no BC at all — it follows from
> the physics (§5.7) — and there is no radiation BC type: external incident
> radiation (cone heater, panel) is specified as a **thermal** BC
> (`RadiativeFluxBC`, §5.2.6), and the radiation model chooses where it is
> absorbed (§5.6).

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
| `RadiativeFluxBC` | `a·q_inc` (model-dependent) | `flux`, `absorptivity` (`nothing` ⇒ `a_eff`) |
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
`dq_BC/dT_s = 0`, so it does not couple to the surface temperature: alone on a
face it is applied in closed form, without the surface Newton of §5.4 (the
recorded face temperature is `T_cell + q·d/k`).

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
the Newton solve converge quickly.

**Blowing.** Strong outgassing thickens the boundary layer and suppresses
convective exchange, and a gas denser than the free stream can pool on the
face. How a `ConvectiveBC` responds is its `blowing` model, `NoBlowing()` by
default (the coefficient stays `h` whatever leaves the face). The surface
balance hands the model the outward mass flux of every gas component at the
face at each iterate. A sealed face has zero species fluxes and is unaffected. Four models ship:

```julia
ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = CouetteFilm())        # h·B/(e^B − 1), B = ṁ″·c_p,g/h
ConvectiveBC(h = 11.2, T_inf = 360.0, blowing = CouetteFilm(0.31))    # film law with a fitted constant
ConvectiveBC(h = 11.2, T_inf = 360.0, blowing = PropertyRatioFilm(NITROGEN))
ConvectiveBC(h = 11.2, T_inf = 360.0, blowing = PoolingFilm(NITROGEN; volatiles = pyrolysate))
ConvectiveBC(h = 11.2, T_inf = 360.0,
             blowing = FunctionBlowing((face, h, T_inf) -> exp(-1500.0 * sum(face.ṁ) / h)))
```

- `CouetteFilm(λ = 0.5)`: the classical film law `h → h·R_b` with the
  blowing ratio `R_b = B_λ/(e^{B_λ} − 1)` and the film parameter
  `B_λ = 2λ·Σ_j J_j c_p,j(T_s)/h`, using signed outward species fluxes.
  This handles suction and counterflow; zero net mass flow need not give a zero
  heat-capacity flux.
  `λ = 0.5` is the Couette film (`B_λ = B`). At a wood-like peak
  `ṁ″ ≈ 10 g·m⁻²·s⁻¹` against `h ≈ 10` this is a 35–50 % reduction of `h`;
  at PMMA-like peaks (`B ≈ 4`) convective cooling is essentially
  extinguished. A smaller `λ` fitted to a resolved flow holds for that gas
  only.
- `PropertyRatioFilm(free_stream; n_Γ = 0.6)`: the film law for any light
  gas, `R_b = B_Γ/(e^{B_Γ} − 1)` with the film parameter
  `B_Γ = (ṁ″·c_p,∞/h)·(Γ_∞/Γ_g)^n_Γ`, `Γ = k/c_p` of the free stream
  and of the leaving mixture (heat capacity by mass fraction, conductivity by
  Wilke's rule). No constant is refitted for a new gas.
- `PoolingFilm(free_stream; n_Γ, n_D, ψ0, Λc, Λw)`: the same `B_Γ` with a
  pooling factor `ψ`, `R_b = B_Γ/(e^{B_Γ/ψ} − 1)`, that falls from 1 to the
  floor `ψ_0` (keyword `ψ0`) when the gas at the face is denser than the free
  stream, switching over a logistic in the face buoyancy `Λ` centred at `Λ_c`
  with width `Λ_w` (keywords `Λc`, `Λw`). The face
  composition uses the solved mass-film state when available; otherwise it
  follows the empirical species film balance with binary diffusivities and
  requires zero free-stream species flux. The three switch constants belong to
  the apparatus flow (the defaults are the CAPA chamber's); the gas enters
  through its properties. Describe the volatiles either through the
  material's gas components, with Fuller diffusion volumes
  `diffusion_volume = (name = V, ...)`, or as one lumped `GasSpecies`
  (`volatiles = ...`).
- `FunctionBlowing(f)`: your own blowing ratio `R_b = f(face, h, T_inf)` from the
  `BlowingFace` (outward species fluxes `face.ṁ`, surface temperature
  `face.T_s`, `face.material`). `f` must return 1 at zero species fluxes, be smooth and
  accept dual numbers (branch only on `ForwardDiff.value`).

The empirical property-ratio and pooling closures support outward products
or pure free-stream suction and reject mixed counterflow. Use the signed
Couette closure for counterflow.

The free stream is a `GasSpecies` (molar mass, `c_p(T)`, `k(T)`, Fuller
diffusion volume, optional viscosity); `NITROGEN` ships, and
`nasa7_heat_capacity` builds a heat capacity from NASA-7 coefficients. Each
convective member of a combined BC applies its own model. The solution
records the applied ratio `R_b` and, for the pooling model, `ψ`, the face fraction
and the face buoyancy in `extras.surface_blowing_ratio`,
`surface_pooling_factor`, `surface_face_fraction` and
`surface_face_buoyancy`. The shipped examples were calibrated with the bare
coefficient (and without flame feedback — the two partially cancel; see Technical
Reference §1.8.1 before recalibrating).

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
  emissivity** `ε_eff = Σ_j v_j ε_j(T)` (volume-weighted, each component
  emissivity evaluated at the boundary-cell temperature) computed from the live
  composition. This is the recommended setting: it automatically tracks the
  change from virgin polymer to char as the material decomposes, and any
  temperature dependence you gave the component `ε`. You can also pass a
  constant in `[0, 1]` or a function `ε(T)`.
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

Forces the surface to a prescribed temperature `T_s = T`. Alone on a face it is
applied in closed form, without the surface Newton of §5.4: the face temperature
is `T` and the flux conducted into the cell is `k_eff (T − T_cell)/d`, which the
solution reports as `extras.surface_heat_flux_net`. Inside a
`CombinedThermalBC` it enters the surface balance as a large conductance
(`1e8` W·m⁻²·K⁻¹) to `T`, so the face then holds `T` only approximately. Use it for verification problems or when an experiment fixes the
surface temperature directly. The single field `T` may be constant or `T(t)`.

### 5.2.6 `RadiativeFluxBC` — external incident radiation

This is the BC for a **cone heater, furnace, or radiant panel** — an external
radiation source heating the surface. It is distinct from `RadiativeBC` (which
is the surface's *own* re-radiation loss) and from `HeatFluxBC` (which injects a
bare flux with no absorptivity or in-depth treatment).

```julia
bc = RadiativeFluxBC(flux = 25e3)                          # incident 25 kW/m², a = a_eff
bc = RadiativeFluxBC(flux = 50e3, absorptivity = 0.9)      # fixed absorptivity
bc = RadiativeFluxBC(flux = t -> 25e3*(1 - exp(-t/10)))    # time-varying flux
```

Fields:

- `flux` — incident radiative flux `q_inc` in W/m² (constant or `q_inc(t)`).
- `absorptivity` — surface absorptivity `a` in `[0, 1]`. Default `nothing`,
  which uses the material's effective surface absorptivity
  `a_eff = Σ_j v_j a_j`: each component's own `absorptivity` when set, else its
  emissivity `ε_j(T)` (Kirchhoff's law). With no component absorptivity set this
  is exactly `ε_eff`; setting one lets the absorbed heater flux differ from the
  re-radiation without splitting the BC. Either way it tracks composition and
  temperature automatically. You can override with a constant or `a(t)`.

**Its behavior depends on the `radiation_model` keyword passed to `solve()`**
(see §7.4 and Technical Reference §8):

| `radiation_model` (in `solve`) | What `RadiativeFluxBC` does |
|---|---|
| `NO_RADIATION` | Absorbed flux `a·q_inc` deposited as a surface boundary heat flux |
| `SURFACE_ABSORPTION` | Same — absorbed at the surface (opaque material) |
| `BEER_LAMBERT` | Boundary flux is **zero**; the absorbed `a·q_inc` enters the first live cell as the intensity of its exposed face and is absorbed in depth (`extras.volumetric_absorption`); what leaves the last face is transmitted (`extras.radiation_transmitted`) |

Under `BEER_LAMBERT` the incident radiation enters through the `:top` face only:
a `RadiativeFluxBC` on any other tag is refused with an `ArgumentError`. With
`absorptivity = nothing` the absorbed fraction is the effective absorptivity of
the first live cell.

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
BCs. Each term keeps its own `T_inf` (§5.8).

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
        top    = AmbientGasBC(GasReservoir(P=0., T=300.)),
        bottom = ImpermeableBC(),
    ),
)
# Pair with:  solve(problem; radiation_model = SURFACE_ABSORPTION)
```

---

## 5.3 Mass boundary conditions

Species fluxes are positive outward, in kg/(m² s). Assign one complete mass
condition to each tag; mass conditions do not combine with `+`.

### 5.3.1 Impermeable and prescribed species flux

```julia
sealed = ImpermeableBC()
prescribed = MassFluxBC(fluxes=(product=1e-4,air=-2e-5),
                       inflow_temperature=(air=350.,))
```

An inward species flux needs its incoming temperature. Fluxes and incoming
temperatures may be functions of time. A prescribed flux requires no internal
transport coefficient; its carried energy is still assembled.

### 5.3.2 Complete reservoirs

```julia
reservoir = GasReservoir(P=101325., T=300., mole_fractions=(air=1.,))
mass = (top=AmbientGasBC(reservoir; film=MassTransferFilm(0.02)),
        bottom=ImpermeableBC())
```

Reservoir mole fractions are named, nonnegative and sum to one. Unknown nonzero
species are rejected. Pressure, temperature and composition may depend on time.
An explicit vacuum `GasReservoir(P=0.,T=300.)` may omit composition. Atmospheric
exchange requires the atmosphere's tracked species and appropriate initial pore
gas. Both escape and ingress follow the same constitutive path.

### 5.3.3 Optional mass film

`MassTransferFilm(g)` is an isobaric common-diffusivity film with conductance
`g=ρD/δ` in kg/(m² s), supplied as a number or temperature law. A finite film
solves surface fractions jointly with the thermal balance. Infinite conductance
fixes composition to the reservoir. Zero conductance permits only donor bulk
flow. The model and signed blowing/reversal limits are in Technical Reference §9.8.

---

## 5.4 Surface temperature and carried energy

The surface temperature is determined from the thermal condition and conduction
through the current half-cell distance. With a mass film or hydraulic resistance,
surface temperature, fractions and pressure are solved together whenever their
equations couple. A bounded, scaled local Newton solve must converge and have a
resolved nonsingular Jacobian; otherwise the solver diagnoses the local failure.

Outward gas carries surface enthalpy; inward gas carries its specified incoming
enthalpy. The reported boundary heat uses the reservoir control surface. For an
LTE wall heat film, `q_reservoir=q_surface+H_reservoir-H_surface`, so total outward
energy is `H_reservoir-q_reservoir=H_surface-q_surface`. An adiabatic boundary
has no hidden incoming-gas heating. A convective film and its reservoir must use
the same incoming stream temperature; radiative surroundings can differ.

`CouetteFilm` uses the signed sum of species heat-capacity fluxes. Empirical
pooling/property-ratio models support outward products or pure freestream
suction and diagnose unsupported mixed counterflow. See Technical Reference
§9.9 for the control surface, approximation and supported domains.

---

## 5.5 Pressure-gradient and hydraulic-resistance boundaries

```julia
gradient = NeumannPressureBC(gradient=-1e5,T=300.,mole_fractions=(air=1.,))
reservoir = GasReservoir(P=101325.,T=300.,mole_fractions=(air=1.,))
resistive = ConvectivePressureBC(reservoir;resistance=1e6,flow=:darcy_velocity)
```

The pressure gradient is outward [Pa/m]. It determines surface pressure through
the current half-cell distance; accepted surface pressure must be nonnegative.
The resistance law is `P_surface-P_reservoir=resistance*flow`. Choose
`:darcy_velocity` [m/s], `:mass` [kg/(m² s)] or `:molar` [mol/(m² s)]; resistance
units are pressure divided by that quantity. A zero resistance fixes pressure.
`ConvectivePressureBC` also accepts `film=MassTransferFilm(...)`. Both variants
supply complete composition and incoming energy data and support either gas model.

---

## 5.6 Radiation models and the boundary conditions

The three radiation models are selected by the `radiation_model` keyword of
`solve()` (default `NO_RADIATION`), not by a BC in the set, and there is no
radiation BC type: the incoming radiation is always a `RadiativeFluxBC` in the
thermal slot (§5.2.6), and the surface's own emission a `RadiativeBC` (§5.2.4).

| `radiation_model` | Where the absorbed incident flux goes | Requirements |
|---|---|---|
| `NO_RADIATION` (default) | the surface balance of the face carrying the `RadiativeFluxBC` | — |
| `SURFACE_ABSORPTION` | the same (opaque material) | — |
| `BEER_LAMBERT` | in depth, through algebraic face intensities attenuated cell by cell with the components' `α` | the `RadiativeFluxBC` on `:top` only; planar coordinates; `α > 0` on the absorbing components |

Under `BEER_LAMBERT` the face intensities are unknowns of the state, each tied to
its neighbour by the exact attenuation across one cell, so the Jacobian stays
banded (Technical Reference §8.3). An in-depth radiative transfer with
re-emission (a P1 model) is not implemented.

---

## 5.7 Mesh motion in ALE simulations (no boundary conditions needed)

Mesh motion is **not configured through boundary conditions** — there are no
mesh-motion BC types, and `BoundaryConditionSet` carries only `thermal` and
`mass` slots. With `solve(problem; use_ale = true)` the mesh's pinned node (on a
line the back node, `z = 0`) stays at rest and every other node moves with the
solid velocity plus a correction, `w = v_s + c`. The solid velocity follows from
the cells' volumetric strain rates through one algebraic balance row per cell,
so the top surface follows the material as it shrinks or swells. The correction
`c` is zero unless an equidistribution correction is installed, which `setup`
does where the mesh moves, cells are deactivated and a layer can vanish; it
moves the nodes inside a vanishing layer so that its cells thin together, and
keeps layer interfaces with the solid (Chapter 10; Technical Reference §11).

---

## 5.8 Choosing ambient temperatures

Two thermal BCs carry an ambient temperature `T_∞`: `ConvectiveBC` and
`RadiativeBC` each have their own `T_inf`, which enters only that BC's own flux.
BCs without an ambient (`HeatFluxBC`, `AdiabaticBC`, `RadiativeFluxBC`,
`DirichletThermalBC`) carry none. In a `CombinedThermalBC` each term keeps its
own `T_inf`: a top BC of `RadiativeBC(T_inf = 330)` plus
`ConvectiveBC(T_inf = 300)` re-radiates to 330 K and convects to 300 K. The
surface energy balance (§5.4) has no ambient reference temperature of its own.

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

> **Note on time-varying radiation flux.** Every condition parameter is
> evaluated at the time of each residual evaluation — constant and time-varying
> alike, with no up-front constancy detection; under `BEER_LAMBERT` the
> exposed-face intensity row is `I_e = a·q_inc(t)` at that time. A ramped or modulated
> flux is therefore always tracked. You do not need to do anything special —
> just supply `flux = t -> …`.

---

## 5.10 Initial conditions

Initial conditions are supplied to `PyrolysisProblem` as functions of the spatial
coordinate `z` (in meters, with `z = 0` at the bottom, `z = L` at the top).
`setup` evaluates them at the cell centres of the mesh's reference positions:
the initial temperature of cell `i` is `T_initial(z_i)` and its initial mass of
component `j` is `m_ij = ξ_initial[j](z_i)·V_i`, since the component masses are
the states (a profile that varies inside a cell is sampled at its centre, not
averaged). Nothing is stored on the mesh. Two inputs:

- `T_initial` — a function `z -> T` returning the initial temperature in K.
  Default: `z -> 300.0`.
- `ξ_initial` — the initial mass concentrations `ξ_j(z)` in kg/m³, one entry per
  component. Accepts a `Vector` of functions (one per component), a single
  function (used for every component), or you may omit it for the auto-default.

To see the values a mesh gets, call `initial_fields` (exported by the `Problem`
submodule), which returns the cell-centre temperatures and the `NC × n`
concentration matrix:

```julia
T0, ξ0 = Pyrolysis.Problem.initial_fields(problem, problem.mesh)
```

The integrator refuses an initial state in which a component mass is more
negative than its tolerance or a cell is at or below absolute zero (an
`ArgumentError` from `solve`), so the profiles must be non-negative.

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
component, if it is a `SolidComponent`, is placed at its reference density** (its density
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

The initial masses also set the per-cell **dry-solid mass**, the denominator of
the pyrolysis progress `χ` whose column average `χ̄` a lateral-shrinkage law
reads, `m_{dry,i} = V_i · Σ_{solid} ξ_{j,i}(0)` over the
`SolidComponent`s (see Technical Reference §12.8 and §10.3).

**Layered samples.** Put a material interface on a mesh node (for example with
`create_mesh_from_nodes`) and give the components a step there, so every cell
starts with one layer's composition. `setup` finds the layer interfaces from
which condensed components each cell starts with, and on a moving mesh with
deactivation it keeps those nodes with the solid (Chapter 10); the `lagrangian`
keyword of `solve` names them explicitly.

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

A smooth (`tanh`) interface spreads a composition change over several cells and
can soften the startup transient of a single-material gradient. For a layered
sample, prefer a step on a node: a smooth profile leaves every component present
in the transition cells, so the interface that `setup` detects moves to wherever
the profile's tail reaches exactly zero.

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
        top    = AmbientGasBC(GasReservoir(P=0., T=300.)),
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
        top = RadiativeFluxBC(flux = 50e3) +          # a defaults to a_eff (= ε_eff by Kirchhoff)
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
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
    mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)
```

### 5.11.4 A pressure reservoir for a gas-bearing char

For a material containing `:air`, an atmospheric boundary is:

```julia
reservoir = GasReservoir(P=101325., T=300., mole_fractions=(air=1.,))
bcs = BoundaryConditionSet(
    thermal=(top=ConvectiveBC(h=10., T_inf=300.), bottom=AdiabaticBC()),
    mass=(top=AmbientGasBC(reservoir), bottom=ImpermeableBC()))
```

This permits both ingress and escape. Add `film=MassTransferFilm(g)` for
finite external composition transfer. Use `ConvectivePressureBC` with a
named flow measure when the exterior instead has a hydraulic resistance.

---

## 5.12 Caveats and pitfalls

- **Incoming energy needs complete gas data.** Outflow carries enthalpy at the
  surface temperature; inflow carries its specified reservoir enthalpy. An LTE
  thermal film accounts for the adjustment between those control surfaces in
  its coupled surface balance. Adiabatic and prescribed-heat boundaries use
  reservoir carried enthalpy directly. Use the assembled boundary diagnostics
  to avoid adding incoming energy twice (§5.4).

- **`radiation_model` defaults to `NO_RADIATION`.** A `RadiativeFluxBC` only does
  something if you pass an appropriate `radiation_model` to `solve()`. With the
  default `NO_RADIATION`, an incident `RadiativeFluxBC` is still absorbed at the
  surface (the absorbed flux is applied as a boundary heat flux), but choose
  `SURFACE_ABSORPTION` (opaque) or `BEER_LAMBERT` (semi-transparent) explicitly
  to make the intent clear and to enable the in-depth distribution. Under
  `BEER_LAMBERT` the `RadiativeFluxBC` must sit on `:top`.

- **Only thermal and mass BCs go in the set.** `BoundaryConditionSet` accepts
  only `AbstractThermalBC` and `AbstractMassBC` payloads — a wrong-kind BC
  throws at construction. Mesh motion needs no BC at all (it is automatic
  with `use_ale = true`), and there is no radiation BC type.

- **Complete gas data.** Pressure boundaries use the chosen gas model and
  require species and incoming-energy data. Both ingress and escape are
  represented. Sealed gas still has a pore EOS.

- **Don't double-specify a radiation source.** Use `RadiativeFluxBC` (not also
  `HeatFluxBC`) for the cone/panel; use `RadiativeBC` for the surface's own
  re-radiation loss. They are different physical terms.

- **Match `ξ_initial` to your component order.** The auto-default only fills
  component 1 (if solid). For lumped gas transport materials, initialize the pore gas at
  ambient pressure with `background_gas_concentration(material, :air, ξ_condensed, T)`
  (exported at the top level): a pore space started as a vacuum against a
  101325 Pa boundary draws early product gas inward. For multi-solid or moisture-bearing materials, supply
  `ξ_initial` explicitly so virgin solid, char, moisture, and gases start with
  the right concentrations.

For deeper background on the boundary physics, the surface Newton solve, and the
radiation models, see Technical Reference §12 (boundary and initial conditions),
§8 (thermal radiation), §7 (energy assembly and the gas enthalpy flux), and §9
(gas transport and pressure).
