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
face at each iterate, so the back face is unaffected automatically (no
outflow). Four models ship:

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
  `B_λ = 2λ·ṁ″·c_p,g/h`, `c_p,g` the heat capacity of the gas leaving the face.
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
  composition comes from a species-by-species film balance with the binary
  diffusivities into the free stream. The three switch constants belong to
  the apparatus flow (the defaults are the CAPA chamber's); the gas enters
  through its properties. Describe the volatiles either through the
  material's gas components, with Fuller diffusion volumes
  `diffusion_volume = (name = V, ...)`, or as one lumped `GasSpecies`
  (`volatiles = ...`).
- `FunctionBlowing(f)`: your own blowing ratio `R_b = f(face, h, T_inf)` from the
  `BlowingFace` (outward species fluxes `face.ṁ`, surface temperature
  `face.T_s`, `face.material`). `f` must return 1 at zero outflow, be smooth and
  accept dual numbers (branch only on `ForwardDiff.value`).

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
mode. All of these types (and the pressure BCs) are exported by
`using Pyrolysis`; `CombinedMassBC` is produced by the `+` operator rather than
called directly. Mass conditions act only when the problem has gas transport:
a material with a gas component whose components carry a gas-transfer
coefficient, or a mass condition that passes mass (Chapter 7).

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
  `J_j = −λ_eff ξ_{j,s}/d`, where `λ_eff` is the surface cell's effective
  gas-transfer coefficient and `d` the distance from its center to the boundary
  face. A finite `h_m` adds a series resistance:
  `J_j = −λ_eff h_m ξ_{j,s}/(λ_eff + h_m d)`.

This is the most physically consistent escape BC because it ties the surface
loss to the same internal transport law used everywhere else. Use it (with `Inf` or a
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

![Surface cell with incoming radiation, inward conduction and outward radiative and convective losses; gas enthalpy appears separately in the outward face flux.](../../assets/figures/06_surface_energy.svg)

**Figure 6. Surface energy transfers (schematic, surface absorption).** Incoming radiation balances conduction into the solid and convective and radiative losses. The surface-temperature calculation uses $T_s$, the adjacent cell temperature and distance $d$. Gas outflow modifies convection through blowing; its sensible enthalpy enters the face energy flux separately, as shown below. The arrows depict losses from a hot surface; signed exchanges can reverse. Beer–Lambert absorption instead deposits radiation within the cells (Figure 7).

Thermal BCs are not applied as a raw flux in the cell; instead the solver finds
the **surface (face) temperature** `T_s` that makes the boundary energy balance
close, then conducts the resulting flux into the adjacent cell. The balance is

```
q_BC(T_s) = k_eff · (T_s − T_cell) / d
```

— the applied BC flux equals the heat conducted from the surface face into the
first cell center, with `k_eff` the cell's effective conductivity, `T_cell` its
trial temperature, and `d` the face-to-center distance. It is solved by a Newton
iteration on `T_s` (see Technical Reference §12.2). You do not call this
yourself — it runs inside the residual — but it explains a few user-visible
behaviors:

- **Closed-form faces skip it.** An `AdiabaticBC`, a lone `HeatFluxBC` and a
  lone `DirichletThermalBC` give the face flux directly; every other condition
  (and every `CombinedThermalBC`) goes through the Newton iteration.
- **It always has a root.** The Newton starts from the cell temperature clamped
  to the property bracket `[200, 3000]` K and evaluates the condition at the
  *trial clamp* of the iterate (the smooth clamp into the bracket with 1 K
  transitions), so the balance has a root for any cell temperature and the
  iterate itself is never clamped. A root outside the bracket belongs to a trial
  state that the integrator's error control rejects.
- **Convergence** is on the flux residual, relative tolerance `1e-10` of the
  larger of the applied and conducted fluxes (with a round-off floor), in at
  most 50 iterations; on the cone recipe of §5.11.1 (20 cells, 120 s, with a
  blowing model) it converged in 4 to 5 iterations at every accepted step. A
  balance that does not converge is not an error: the flux is taken at the last
  iterate, no warning is printed, and the per-face diagnostics record
  `converged` and `iterations`.
- **It is differentiated exactly.** The iteration carries dual numbers, so the
  Jacobian's derivatives with respect to the state flow through `T_s`.
- The converged surface temperature is available afterwards as
  `solution.surface_T` (and per tag in `solution.boundary_T`), and the flux
  decomposition (absorbed, re-radiation, convection, conducted, net) in
  `solution.extras` (§9.3).

The balance has **no surface gas-enthalpy term**. The gas leaving through the
face carries its sensible enthalpy at the face temperature, `Σ_j J_j h_j(T_s)`,
in the face's energy flux, added once the balance has fixed `T_s`; a surface
term `Σ_g ṁ_g Δh_g(T_s, T_∞)` in the balance would count that enthalpy a second
time. The outward gas flux enters the balance only through the blowing model of
a `ConvectiveBC` (§5.2.3).

---

## 5.5 Pressure boundary conditions (Darcy mode only)

When the material has a permeability (so the solver runs in Darcy–Fick transport
mode, see Technical Reference §9), gas transport is pressure-driven and you can
prescribe **pressure** at a boundary. Pressure BCs are `AbstractMassBC` subtypes,
so they go in the `mass` slot of the `BoundaryConditionSet`. They contribute no
film (diffusive) mass flux; instead they set the outward Darcy velocity `v` of
the boundary face from the boundary cell's pressure `P` and mobility `κ/μ`. An
outflow (`v > 0`) carries the cell's gas out at its intrinsic density `ξ_j/φ`;
an inflow carries none of the tracked gases (the ambient is taken to hold none
of them). A mass condition that is not a pressure condition drives no Darcy
velocity, so an `ImpermeableBC` face is sealed under Darcy–Fick as well. In
Fickian (non-Darcy) mode pressure BCs have no effect.

| BC type | Outward Darcy velocity | Parameters (defaults) |
|---|---|---|
| `DirichletPressureBC` | `v = (κ/μ)(P − P_bc)/d` | `P` (`101325.0`) |
| `NeumannPressureBC` | `v = −(κ/μ)·∂P/∂n` | `∂P∂n` (`0.0` = no flow) |
| `ConvectivePressureBC` | `v = h_P·(P − P_ext)` | `h_P`, `P_ext` |

All three are exported by `using Pyrolysis`. A pressure condition combines with
a film condition through `+`, e.g. `DiffusiveMassBC(0.01) +
DirichletPressureBC(101325.0)` (diffusive escape plus pressure-driven venting),
as `examples/08_darcy_fick_demo.jl` does.

### 5.5.1 `DirichletPressureBC` — fixed pressure

```julia
bc = DirichletPressureBC()           # default: atmospheric, 101325 Pa
bc = DirichletPressureBC(101325.0)   # explicit atmospheric
```

Holds the boundary at pressure `P` (Pa). The outward Darcy velocity is the
one-sided gradient between the interior cell and the boundary,
`v = (κ/μ)·(P_interior − P_bc)/d`. This is the natural exposed-surface
condition: the outer face sits at atmospheric pressure and gas flows out down the
pressure gradient built up by pyrolysis.

### 5.5.2 `NeumannPressureBC` — prescribed gradient

```julia
bc = NeumannPressureBC()         # default: zero gradient (no Darcy flow)
bc = NeumannPressureBC(0.0)      # same as above
bc = NeumannPressureBC(-1000.0)  # pressure falling outward at 1 kPa/m: an outflow
```

Sets the pressure gradient `∂P/∂n` (Pa/m) along the **outward** normal, default
`0.0`. A zero gradient is a no-flow (impermeable) back face for the pressure
field — the usual choice for a sealed substrate in Darcy mode. A positive
gradient (pressure rising outward) would drive an inflow, which carries no gas.

### 5.5.3 `ConvectivePressureBC` — Robin (resistive) outflow

```julia
bc = ConvectivePressureBC(1e-6, 101325.0)   # finite surface resistance, atmospheric exterior
```

A Robin-type condition, `v·n = h_P·(P_internal − P_external)`, for a surface
with finite resistance to gas flow (a char skin, a coating). Fields:

- `h_P` — pressure-transfer coefficient `h_P` in m·Pa⁻¹·s⁻¹.
- `P_ext` — external pressure `P_external` in Pa.

The velocity is gated to zero on an impermeable boundary cell. Limits:
`h_P → ∞` approaches `DirichletPressureBC`; `h_P → 0` approaches an impermeable
(zero-flow) boundary.

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
        top = RadiativeFluxBC(flux = 50e3) +          # a defaults to a_eff (= ε_eff by Kirchhoff)
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

- **The escaping gas carries no surface enthalpy term.** The gas enthalpy flux
  `Σ_j J_j h_j` through every face, the boundary face included at the surface
  temperature, already carries the enthalpy of the escaping gas (Technical
  Reference §7, §12.2); a surface term `Σ_g ṁ_g Δh_g(T_s, T_∞)` on top of it
  would double-count that enthalpy, a spurious sink that under-predicts surface
  temperature and mass-loss rate. Parameter sets calibrated with a code that
  adds a surface gas-enthalpy term may need re-calibration.

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

- **Pressure BCs need Darcy mode.** `DirichletPressureBC`/`NeumannPressureBC`/
  `ConvectivePressureBC` only have an effect when the material has permeability
  (Darcy–Fick transport). In Fickian mode they contribute zero flux and are
  silently inert — use `ImpermeableBC`/`DiffusiveMassBC`/`ConvectiveMassBC` for
  diffusion-only gas escape. The pressure conditions vent outflow only; an
  inward Darcy velocity carries no gas.

- **Don't double-specify a radiation source.** Use `RadiativeFluxBC` (not also
  `HeatFluxBC`) for the cone/panel; use `RadiativeBC` for the surface's own
  re-radiation loss. They are different physical terms.

- **Match `ξ_initial` to your component order.** The auto-default only fills
  component 1 (if solid). For Darcy–Fick materials, initialize the pore gas at
  ambient pressure with `background_gas_concentration(material, :air, ξ_condensed, T)`
  (exported at the top level): a pore space started as a vacuum against a
  101325 Pa boundary draws early product gas inward. For multi-solid or moisture-bearing materials, supply
  `ξ_initial` explicitly so virgin solid, char, moisture, and gases start with
  the right concentrations.

For deeper background on the boundary physics, the surface Newton solve, and the
radiation models, see Technical Reference §12 (boundary and initial conditions),
§8 (thermal radiation), §7 (energy assembly and the gas enthalpy flux), and §9
(gas transport and pressure).
