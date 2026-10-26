# 2. Quick-Start Tutorial: Your First Cone-Calorimeter Run

This chapter walks you through a complete, runnable cone-calorimeter simulation
from an empty Julia session to plotted results. We model a slab of PMMA
(poly(methyl methacrylate)) heated from above by a radiant cone, decomposing
into a monomer gas plus a small char residue. By the end you will have built a
material, defined boundary conditions, created a mesh, assembled a problem,
solved it, and extracted surface temperature and mass-loss rate. Every code
snippet uses the public API and can be pasted into the REPL as-is. The deeper
"why" of each physics choice is deferred to the Technical Reference; the goal
here is to get you to a working run and a mental model you can reuse.

The mental model is a pipeline that this chapter follows in order:

> **material → boundary conditions → mesh → problem → solve → solution**

You define *what* decomposes (the `Material`), *how* energy and mass cross the
edges (the `BoundaryConditionSet`), *where* the domain lives (the mesh), bundle
them into a `PyrolysisProblem`, hand that to `solve`, and read results off the
returned `PyrolysisSolution`.

The deck uses the material, atmospheric pore gas and boundary conditions of
`examples/06_cone_calorimeter.jl`.

---

## 2.1 Before you begin: coordinates, units, and conventions

Two conventions are worth fixing in your head before writing any code, because
they govern how every boundary condition and result is interpreted.

**Coordinate axis.** The single spatial coordinate is `z`, running through the
thickness of the sample:

- `z = 0` is the **bottom / substrate** (boundary tag `:bottom`).
- `z = L` is the **top / exposed surface** (boundary tag `:top`).

Heat enters from the top; the material recedes toward the bottom as it
decomposes. Positive flux is in the `+z` direction (bottom → top). See
Technical Reference §13 for the finite-volume layout and §12 for the boundary
conditions on the tagged faces.

**Units.** Pyrolysis.jl is strict SI throughout — metres, kilograms, seconds,
kelvin, joules, watts, pascals. There are no implicit "kW" or "mm". A 25 kW/m²
cone flux is written `25e3`; a 5.8 mm sample is `5.8e-3`. The composition is
given and reported as the **mass concentration** `ξⱼ` in kg/m³ (not a mass
fraction). For the full symbol table see Technical Reference §2.

**One import.** Everything in this tutorial comes from a single package, plus
`Plots` for the figures:

```julia
using Pyrolysis
using Plots          # optional; enables plot(sol) recipes
```

`using Pyrolysis` brings in the curated public surface: the material and
component constructors, the boundary-condition types, the mesh constructors,
`PyrolysisProblem`, `PyrolysisSolution`, `solve`, and the solution's accessors.
Loading `Plots` activates the plotting recipes described in §2.9.

---

## 2.2 Step 1 — Define the material

A `Material` is a tuple of **components** (the chemical species, each a phase
with its own thermophysical properties) plus a tuple of **reactions** (the
Arrhenius kinetics that convert one species into others), plus the
**mixing rules** used to blend per-component properties into effective
mixture properties. This atmospheric PMMA model tracks four components:

1. **Virgin PMMA** — the solid that decomposes.
2. **Char** — a small solid residue left behind.
3. **MMA gas** — the monomer that escapes through the matrix.
4. **Air** — the initial pore gas and external reservoir species.

### 2.2.1 Temperature-dependent properties

PMMA's heat capacity and conductivity vary with temperature. The most general
way to encode an arbitrary `f(T)` is a `CallableProperty`, which wraps any
single-argument function:

```julia
# Heat capacity [J/(kg·K)], piecewise about 378 K, held constant below 300 K
pmma_cp = CallableProperty(
    T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T),
)

# Thermal conductivity [W/(m·K)], piecewise about 378 K, held constant above 700 K
pmma_k = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))
```

The fits are measured on about 300–700 K and are held constant outside that
range. This is required, not cosmetic: property laws are evaluated on the bracket
200–3000 K (their argument is clamped to it), and the component constructors
sample every law over that bracket and refuse one that leaves its admissible
range. The unclamped low-temperature heat-capacity line is negative below 204 K
and the unclamped conductivity line crosses zero at 872 K, so either would be
refused with an `ArgumentError` naming the temperature.

`CallableProperty` is one of several property-function types. You can also pass
a bare number (auto-converted to a `ConstantProperty`), a `LinearProperty(a, b)`
for `a + b·T`, a `PolynomialProperty`, an `ArrheniusProperty`, a
`TableProperty` for tabulated data, or a `SutherlandProperty(μ_ref, T_ref, S)`
for the gas-viscosity correlation μ_ref·(T/T_ref)^{3/2}·(T_ref+S)/(T+S).
Property types and their derivatives are covered in full in UG 3 and Technical
Reference §4.

### 2.2.2 Component constructors

```julia
virgin = SolidComponent(:virgin; ρ=1210., φ=.01, c=pmma_cp, k=pmma_k,
                        ε=.96, λ=1e-5)
char = SolidComponent(:char; ρ=1210., φ=.01, c=pmma_cp, k=pmma_k,
                      ε=.96, λ=1e-5)
mma = GaseousComponent(:MMA; M=.10012, c=1500., k=.02)
air = GaseousComponent(:AIR; M=.02897, c=1006., k=.026)
```

The nominal pore fraction and superficial transfer coefficient are explicit
material inputs. Gas has no matrix transfer coefficient. Here `κ` defaults
to zero, so gas moves through the lumped wall channel. An open atmospheric
problem includes a background gas and complete reservoir composition. The
initial air loading is computed from the pore volume and reservoir pressure.
The cone deck uses the nominal pore fraction φ=0.01 qualified over
φ=0.01–0.03 with joint mesh and pressure-tolerance refinement (TR §9.10).
For the complete component options see §16.2.1.

---

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
  endothermic, so `h = 820e3` is positive. (Note: ThermaKin publishes
  the opposite sign convention — see Technical Reference §6 and the §2
  nomenclature note H1.)
- **Single-reactant only:** each reaction consumes exactly one reactant species
  (with mass stoichiometry 1.0). Multi-reactant reactions are not supported.
  Multi-step schemes are modelled as a *tuple of single-reactant reactions* (see
  UG 4).

For a one-product reaction the simpler form `Reaction(:name, 1 => 2; A, E, h)`
omits `yields` (the single product then has stoichiometry 1.0). Temperature
gates `T_min`/`T_max` default to `0.0` and `Inf` (always active); their smooth
tanh ramp width is `gate_width` (default `1.0` K, UG §4.6). See Technical
Reference §6 for the full rate law, the smooth `tanh` temperature gates, and the
depletion limiter.

### 2.2.4 Assemble the `Material`

Finally bundle components, reactions, and mixing rules. `MaterialMixing` selects
how per-component thermal conductivity (`k`), gas-transfer (`λ`), and
permeability (`κ`) are averaged into effective mixture properties:

```julia
pmma = Material(
    name = :PMMA,
    components = (virgin, char, mma, air),
    reactions  = (decomp,),
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),  # k = 0.5·k_parallel + 0.5·k_series
        λ = MixingSpec(PARALLEL),       # condensed bulk-occupancy-weighted transfer
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
mixing theory lives in Technical Reference §5; how to choose rules is in UG 3.

> **Tip — components reference by index, not name (here).** Because we passed
> `virgin`, `char`, `mma`, `air` in that order, the reaction's `1 => (2, 3)` lines up
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
q_external = 0.9813 * 25e3   # incident cone flux [W/m²]: 25 kW/m², radially averaged
h_conv     = 7.2             # convective coefficient [W/(m²·K)]

top_thermal =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +  # incident cone radiation in
    RadiativeBC(ε = 0.96, T_inf = 330.0) +                     # re-radiation loss out
    ConvectiveBC(h = h_conv, T_inf = 330.0)                    # convective loss out
```

The nominal 25 kW/m² is scaled by 0.9813 for the radial non-uniformity of the
heater's flux over the sample, as the example does.

What each piece does:

- **`RadiativeFluxBC(flux, absorptivity)`** is the *incident* external radiation
  from the cone heater. The absorbed flux is `flux × absorptivity`. If you pass
  `absorptivity = nothing` (the default), the material's effective surface
  absorptivity `a_eff = Σ_j v_j a_j` of the exposed cell is used instead, where
  each component's `a_j` is its own `absorptivity` keyword when set and its
  emissivity `ε_j(T)` otherwise (Kirchhoff's law, so `a_eff = ε_eff` unless a
  component sets one). Either way it tracks composition as the surface chars.
  Where this flux is deposited depends on the `radiation_model` you choose at
  solve time (see §2.6).
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

Other thermal BCs available (full reference in UG 5 / Technical Reference §12):

| Type | Flux | Notes |
|---|---|---|
| `AdiabaticBC()` | `0` | perfect insulation |
| `HeatFluxBC(q)` | prescribed `q` or `q(t)` | `+` into domain; not subject to absorptivity |
| `ConvectiveBC(h, T_inf)` | `h(T_∞ − T_s)` | `h`, `T_inf` may be functions of `t` |
| `RadiativeBC(ε, T_inf, F)` | `F·ε·σ·(T_∞⁴ − T_s⁴)` | re-radiation; `ε = nothing` ⇒ `ε_eff` |
| `RadiativeFluxBC(flux, absorptivity)` | `a·flux`, `a` the absorptivity (model-dependent) | incident external radiation |
| `DirichletThermalBC(T)` | enforces `T_s = T` | fixed-temperature face |
| `bc1 + bc2 + …` | sum | `CombinedThermalBC` |

> **Note — surface temperature is solved, not assumed.** In a cell-centred
> finite-volume method the temperature lives at cell *centres*, but the
> radiative loss `σ T_s⁴` is highly nonlinear in the *surface* temperature.
> Pyrolysis.jl therefore solves, at each boundary face, the surface energy
> balance `q_bc(T_s) = k (T_s − T_cell)/d` for the surface temperature `T_s`.
> Gas-bearing faces solve the mass and thermal state jointly when coupled.
> The solved value lands in
> `solution.surface_T`. See UG §5.4 and Technical Reference §12.

### 2.3.2 Mass boundary conditions

Choose one complete species boundary per tag. The top exchanges gas with an
air reservoir; the bottom is sealed. The reservoir supplies pressure,
temperature and mole fractions:

```julia
reservoir = GasReservoir(P=101325., T=330., mole_fractions=(AIR=1.,))
g_film = 0.05 * 101325. * 0.02897 / (8.31446261815324 * 330.)
top_mass = AmbientGasBC(reservoir; film=MassTransferFilm(g_film))
bottom_mass = ImpermeableBC()
```

The finite external mass film uses a superficial conductance in kg/(m² s),
here a transfer velocity of 0.05 m/s times the reservoir density. Flux can
reverse; incoming gas brings its
specified enthalpy. `MassFluxBC` prescribes named outward species fluxes with
incoming temperatures for negative entries. Pressure-gradient and hydraulic
resistance closures also require complete species and thermodynamic data (§5.5).

---

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
`create_uniform_mesh(L, n_cells; coordinates = Planar(), cross_section = 1.0)`,
where `L` is the domain thickness in metres and `n_cells` is the cell count.

```julia
sample_thickness = 5.8e-3   # 5.8 mm
n_cells          = 72       # → 0.08 mm per cell
T_ambient        = 300.0    # initial / ambient temperature [K]

mesh = create_uniform_mesh(sample_thickness, n_cells)
```

The mesh holds geometry only: node positions, faces, cells, and the boundary
tags (`:bottom` at `z = 0`, `:top` at `z = L`). It carries no component count
and no initial state; the initial temperature and composition are the
problem's (§2.5). `cross_section` defaults to `1.0` m² — in a 1D planar model
the cross-section is a gauge normalisation (MLR and mass are reported per unit
area), so leaving it at 1.0 is correct and conventional. See UG 6 for when to
change it, and for cylindrical and spherical coordinates.

To refine the mesh near the heated surface (where the reaction front is
steepest), use `create_stretched_mesh(L, n_cells, stretch; stretch_from = :top)`
instead, with a `stretch` factor above 1.0, or place the nodes yourself with
`create_mesh_from_nodes(z)`. That is optional; a uniform mesh is fine for a
first run.

> **Tip — resolution vs. cost.** 72 cells over 5.8 mm resolves the thermal
> wave and reaction zone comfortably for PMMA. Too few cells smears the surface
> gradient; too many slows the solve, since every step's work grows with the
> cell count. A grid-convergence study is UG §14.3; Technical Reference §18
> describes the convergence studies of the verification suite.

---

## 2.5 Step 4 — Build the problem

`PyrolysisProblem` ties the mesh, material, boundary conditions, initial
conditions, and time span together. Initial conditions are given as
*functions of `z`* (allowing spatially varying profiles), here constant:

```julia
t_final = 1600.0   # simulate 1600 s
ξ_air = background_gas_concentration(pmma, 4, (1210., 0., 0., 0.), T_ambient; P=101325.)

problem = PyrolysisProblem(
    mesh      = mesh,
    material  = pmma,
    bc_set    = bc_set,
    T_initial = z -> T_ambient,                         # uniform 300 K
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0, z -> ξ_air],      # one closure per component
    tspan     = (0.0, t_final),
    dt_initial = 1e-3,                                  # first step [s]
    output_times = collect(range(0.0, t_final, length = 601)),  # save 601 points
)
```

Argument notes:

- **`T_initial`** — a closure `z -> T`. Defaults to `z -> 300.0` if omitted.
- **`ξ_initial`** — a `Vector` of closures, one per component (`z -> ξⱼ`). If
  omitted, the first solid component defaults to its reference density and the
  rest to zero. A single closure is used for every component.
- **`tspan`** — the `(t_start, t_end)` time window in seconds (required).
- **`dt_initial`** — the integrator's first step (default `1e-4`); the
  adaptive integrator adjusts from there. It is also `solve`'s default `dt`.
- **`output_times`** — the times at which the solution is saved. Fewer points
  means less memory and less reconstruction work after the solve; if left empty,
  11–101 evenly spaced times are generated. The state a solve ends on is always
  saved.

When the problem is set up, `T_initial` and `ξ_initial` are evaluated at the
cell centres, and the concentrations become the component masses of each cell
(`mᵢⱼ = ξⱼ(zᵢ) Vᵢ`), the states the integrator carries.

`PyrolysisProblem` does not choose the *physics modes* (a moving or fixed mesh,
the radiation model, deactivation, the energy state) — those are selected as
options to `solve`, described next, together with what the material implies
(UG 7).

---

## 2.6 Step 5 — Solve

`solve(problem; kwargs...)` integrates the problem forward and returns a
`PyrolysisSolution`. With no options it runs on a fixed mesh, with the incident
flux absorbed at the surface, the `Rodas5P` integrator and the `StructuredKLU()`
linear solver. For a cone run we make three deliberate choices: move the mesh
with the receding solid, absorb the cone radiation at the surface, and let cells
deactivate as the sample burns out.

```julia
solution = solve(
    problem;
    abstol = 1e-11,                         # temperature tolerance [K]
    reltol = 1e-9,                          # resolves thin-cell pore volume
    pressure_tolerance = 1e-4,              # absolute local pressure budget [Pa]
    dtmin = 1e-14,                          # smallest step [s] (default 1e-12)
    use_ale = true,                         # move the mesh with the solid as it shrinks
    radiation_model = SURFACE_ABSORPTION,   # absorb the cone flux at the surface
    min_thickness = 0.01e-3,                # stop when the sample is < 0.01 mm thick
    deactivation = true,                    # deactivate burned-out cells
    progress = true,                        # show a progress bar
)
```

The first solve in a Julia session compiles the kernels for this problem's
types — its property, reaction, and boundary-condition types — unless a
precompile workload of the package shares them; the time integrator itself comes
compiled with the package. Later solves in the same session reuse the compiled
code, also with other parameter values (UG 15).

### 2.6.1 The options you just used

| Keyword | Default | What it does here |
|---|---|---|
| `abstol` / `reltol` | `1e-8` / `1e-6` | error tolerances: `abstol` is the temperature tolerance [K]; condensed mass rows use initial condensed mass and `threshold`; gas rows and proposed-state pressure errors use `pressure_tolerance` |
| `dtmin` | `1e-12` | smallest step the integrator may take [s] |
| `use_ale` | `false` | moves the mesh with the solid (surface recession) |
| `radiation_model` | `NO_RADIATION` | `SURFACE_ABSORPTION` absorbs the cone flux at the surface |
| `min_thickness` | `1e-6` | on a moving mesh, ends the solve when the thickness falls below it [m] |
| `deactivation` | `false` | deactivates a cell once its condensed mass falls below `threshold` (default `1e-2`) of its initial condensed mass |
| `progress` | `false` | progress bar |

With `use_ale = true` and `deactivation = true` on a layer that can burn out
completely (here the char is 0.2 % of the mass, below the threshold), `solve`
also moves the mesh relative to the solid with the equidistribution correction,
which keeps the cells of the vanishing layer thinning together so that the layer
goes in few deactivation events. The threshold, the band of near-threshold cells
that go with a crossing cell, the correction, and the supported envelope are
UG 10's subject.

### 2.6.2 Choosing the radiation model

The cone heater's incident flux can be deposited at the surface or in depth;
for an opaque solid like PMMA, surface absorption is the right, fast choice:

| `radiation_model` | Behaviour | Use when |
|---|---|---|
| `NO_RADIATION` (default) | A `RadiativeFluxBC` is absorbed *at the surface*, in the surface energy balance | radiation handled purely as a surface flux |
| `SURFACE_ABSORPTION` | The same: the absorbed flux enters the surface energy balance; component `α` is not read | opaque material (optical thickness ≫ 1) |
| `BEER_LAMBERT` | In-depth exponential absorption using the component `α`; the `RadiativeFluxBC` contributes no boundary flux | semi-transparent material (polymers, foams) |

Two important caveats:

- `NO_RADIATION` and `SURFACE_ABSORPTION` give the same solve: both apply the
  absorbed `RadiativeFluxBC` flux at the face, through the surface energy
  balance. (The example script relies on the default; this tutorial names
  `SURFACE_ABSORPTION` to make the choice visible.)
- Under `BEER_LAMBERT` the *same* `RadiativeFluxBC` is instead absorbed in depth
  (and returns zero as a boundary flux) — switching models changes where the
  energy goes, by design. The face intensities become algebraic unknowns of the
  state: the intensity entering the exposed face is `a·Φ`, with the condition's
  explicit `absorptivity` or, without one, the exposed cell's effective surface
  absorptivity (`= ε_eff` by Kirchhoff's law unless a component sets
  `absorptivity`), and each cell absorbs what the beam loses across it at its
  absorption coefficient. The surface temperature Newton solve still runs (it
  balances conducted, convective, and re-radiated flux). The absorbed power is
  reported in `solution.extras.volumetric_absorption` and the power leaving
  through the back face in `solution.extras.radiation_transmitted`; with every
  `α` zero the beam crosses the slab unabsorbed. `BEER_LAMBERT` needs planar
  coordinates and at most one `RadiativeFluxBC` on the exposed boundary;
  setup rejects multiple incident flux conditions. Sources sharing one
  absorptivity can be combined in a single flux function. See Technical
  Reference §8.

### 2.6.3 Other solve options worth knowing

These are not needed for this run but round out the picture (full reference in
UG 8):

- **`integrator_policy`** (default `default_policy`, the `Rodas5P` Rosenbrock
  method; `fbdf_policy` and `kencarp4_policy` are the alternatives) and
  **`linear_solver`** (default `StructuredKLU()`, or LinearSolve's
  `UMFPACKFactorization()`) — the time integrator and sparse linear solver
  (UG 11).
- **`energy_state`** (default `:temperature`) — `:enthalpy` integrates each
  cell's enthalpy inventory instead of its temperature.
- **`threshold`** (default `1e-2`) — the deactivation threshold `ε_m`, a fraction
  of each cell's initial condensed mass, which also scales the mass tolerances
  (UG 10).
- **`saveat`** (default `problem.output_times`), **`dtmax`**, **`maxiters`** —
  saving and step control.
- **`callback`** — a read-only observer of the accepted steps; a
  `ConservationRecorder()` there integrates the boundary mass and energy fluxes
  over each step, and `conservation_audit(solution, recorder)` closes mass and
  energy against the inventories (UG 9, Technical Reference §16).
- The effect of outgassing on convective exchange is a property of each
  `ConvectiveBC` (its `blowing` model, UG §5.2.3), not a solve option.

> **Unknown keywords are rejected.** `solve` validates its keyword arguments and
> throws an `ArgumentError` on anything it does not recognise, so a typo is
> caught immediately rather than silently ignored.

---

## 2.7 The solution object

`solve` returns a `PyrolysisSolution`. Its fields are concretely typed; a typed
`extras` NamedTuple holds the secondary series. The fields you will use most:

| Field | Type / shape | Meaning |
|---|---|---|
| `solution.t` | `Vector{Float64}` | saved times [s] |
| `solution.T` | `Matrix{Float64}` `(n_cells, n_times)` | cell temperatures [K], `NaN` in a deactivated cell |
| `solution.ξ` | `Array{Float64,3}` `(NC, n_cells, n_times)` | concentrations [kg/m³], `NaN` in a deactivated cell |
| `solution.m` | `Array{Float64,3}` `(NC, n_cells, n_times)` | component masses per cell [kg], zero in a deactivated cell |
| `solution.z` | `Matrix{Float64}` `(n_nodes, n_times)` | node positions [m] on a moving mesh (`0×0` on a fixed one), `NaN` beyond the outermost live node |
| `solution.cell_centers` | `Matrix{Float64}` `(n_cells, n_times)` | cell-centre positions [m], `NaN` in a deactivated cell |
| `solution.surface_T` | `Vector{Float64}` | solved surface temperature per time [K] |
| `solution.mass_loss_rate` | `Vector{Float64}` | mass leaving through every boundary face, per unit current area of the top face [kg/(m²·s)] |
| `solution.retcode` | `Symbol` | integrator status (`:Success`, `:Terminated`, `:Failure`, …) |
| `solution.outcome` | `Symbol` | `:running` unless burnout or the thickness stop ended the solve |
| `solution.events` | `Vector` | the deactivation events, in order |
| `solution.statistics` | `SolverStatistics` | accepted and rejected steps, Jacobians, deactivation invocations, … |
| `solution.solve_time` | `Float64` | wall-clock solve time [s] |

The `extras` NamedTuple adds (among others) the surface heat-flux breakdown,
the live inventories, and these series central to cone analysis:

| `solution.extras.…` | Meaning |
|---|---|
| `gas_generation_rate` | gas produced by the reactions of the live cells, per initial area [kg/(m²·s)] |
| `back_temperature` | back-face (`z = 0`) temperature [K] |
| `total_mass` | remaining sample mass per initial area [kg/m²] |
| `thickness_history` | current sample thickness [m] |
| `MLR_gauge` | mass-loss rate per *initial* area [kg/(m²·s)] |
| `live_cells` | number of live cells |
| `deactivation_times` | the time of every deactivation, one entry per cell [s] |

Always check how the solve ended first:

```julia
println("retcode  = ", solution.retcode)     # Success: the solve reached t_final
println("outcome  = ", solution.outcome)     # running: neither burnout nor min_thickness ended it
println("solve t  = ", round(solution.solve_time, digits = 2), " s")
println("saved    = ", length(solution.t), " time points")
```

This run reaches `t_final = 1600 s` with `retcode = :Success` and outcome
`:running`: the sample has lost 92.7 % of its mass and 0.42 mm of it remain. The
combinations you can meet are:

| `retcode` | `outcome` | Meaning | `succeeded` | `completed` |
|---|---|---|---|---|
| `:Success` | `:running` | reached the end of `tspan` | true | true |
| `:Terminated` | `:burned_out` | the last cell deactivated | true | false |
| `:Terminated` | `:min_thickness` | the thickness fell below `min_thickness` | true | false |
| `:Terminated` | `:interior` | a cell below the threshold lies inward of a live one; the solve warns and ends | false | false |
| `:Failure` | — | the integrator gave up, or an event left an inadmissible state (with a warning) | false | false |

`succeeded(solution)` and `completed(solution)` read these for you (UG 9). Run
this deck to 2400 s and the layer burns out: in our run the cells deactivated in
two events, at about 1040 s and 1935 s. The char, 0.2 % of the initial mass,
collects in the last live cell, which keeps more than 1 % of its own initial
condensed mass and so stays live, 0.012 mm thick: that solve also ends with
`:Success`.

### 2.7.1 MLR vs. gas-generation rate (a common gotcha)

There are two related "mass-loss" quantities and they are *not* identical:

- **`solution.mass_loss_rate`** — the actual gas mass flux leaving the sample,
  per current cross-sectional area [kg/(m²·s)]. It counts both faces; with the
  usual impermeable holder the bottom term is zero, so this is the surface
  MLR. (The strictly top-face flux is `solution.extras.surface_mass_flux`.)
- **`solution.extras.gas_generation_rate`** — the gas the reactions produce,
  integrated over the live cells, per *initial* area. It is a convenient proxy
  for MLR and is what the shipped cone example plots, but it leads the surface
  flux slightly because gas takes time to escape.

For a typical cone-calorimeter run with a **fixed** cross-sectional area (no
lateral shrinkage), `A(t) = A₀` always, so `mass_loss_rate` and `MLR_gauge` are
identical and the distinction is moot. The distinction becomes important only if
you give the material a **lateral-shrinkage law** (`lateral_shrinkage_law`,
which needs `use_ale = true`), where `A(t)` shrinks over time: `MLR_gauge` then
re-normalises the surface flux by `A(t)/A₀`, making it directly comparable to
standard cone-calorimeter data reported per *original* specimen area. For this
fixed-area run, `mass_loss_rate` and `MLR_gauge` coincide.

---

## 2.8 Step 6 — Extract surface temperature and mass-loss rate

Everything you need is already on the solution object. The simplest extraction
is direct field access:

```julia
times     = solution.t
T_surface = solution.surface_T                           # [K]
T_back    = solution.extras.back_temperature             # [K]
mlr       = solution.mass_loss_rate .* 1000              # [g/(m²·s)]
gen_rate  = solution.extras.gas_generation_rate .* 1000  # [g/(m²·s)] proxy
total_m   = solution.extras.total_mass                   # [kg/m²]

println("peak surface T = ", round(maximum(T_surface), digits = 1), " K")
println("peak MLR       = ", round(maximum(mlr), digits = 2), " g/(m²·s)")
println("mass lost      = ",
        round((1 - total_m[end] / total_m[1]) * 100, digits = 1), " %")
```

For this deck the peak surface temperature is 672.3 K, the peak MLR
5.49 g/(m²·s), and the mass lost 92.7 % (values from a run of this script).
The surface series read `NaN` once no live cell borders the surface, so take
`maximum(filter(isfinite, T_surface))` on a run that burns out completely.

To pull a quantity at a particular time, find the nearest saved index and read
the field:

```julia
i600 = argmin(abs.(times .- 600.0))                  # index nearest t = 600 s
println("at t ≈ 600 s:  T_s = ", round(T_surface[i600], digits = 1), " K,",
        "  MLR = ", round(mlr[i600], digits = 2), " g/(m²·s)")
```

(672.2 K and 5.47 g/(m²·s), with 4.22 mm of the sample left.)

### 2.8.1 Temperature profiles through the depth

`solution.T[:, i]` is the temperature of every cell at saved time index `i`
(bottom cell first, surface cell last). On a moving mesh the cells move with the
solid, so read their positions at the same time:

```julia
i = argmin(abs.(times .- 120.0))
T_profile = solution.T[:, i]                       # [K], one value per cell
z_cells   = get_cell_positions(solution; at = i)   # [m], the cell centres at that time
```

`get_cell_positions(solution; at = i)` is the column `solution.cell_centers[:, i]`;
a deactivated cell reads `NaN` in both arrays. The node positions are
`solution.z[:, i]` (`get_profile(solution, :z, i)`), and the current thickness
is `solution.extras.thickness_history[i]`.

### 2.8.2 Convenience accessors

The solution's accessors are exported with the package:

| Accessor | Returns |
|---|---|
| `get_times(sol)`, `n_timesteps(sol)` | `sol.t`, its length |
| `get_temperatures(sol)`, `get_concentrations(sol)` | `sol.T`, `sol.ξ` |
| `get_gas_generation_rate(sol)` | `sol.extras.gas_generation_rate` |
| `get_profile(sol, field, i)` | the field `:T`, `:ξ`, `:ξj`, `:m`, `:mj` or `:z` at saved index `i` |
| `get_cell_positions(sol; at = lastindex(sol.t))` | the cell centres at saved index `at` |
| `sol(t)`, `interpolate_solution(sol, t)` | the fields `(; t, T, ξ, m, z, cell_centers)` at any `t` in the span |

```julia
T_at_step = get_profile(solution, :T, i600)   # the same as solution.T[:, i600]
fields    = solution(333.0)                   # the fields at t = 333 s, between saved times
fields.T                                      # cell temperatures [K]
```

`solution(t)` reconstructs the fields from the integrator's interpolant between
the saved states; it refuses a `t` inside a saved interval that a deactivation
event separates, since the interpolant there is not the trajectory (UG 9). For
everyday use, direct field access (`solution.T`, `solution.surface_T`, …) is the
simplest path.

---

## 2.9 Step 7 — Plot the results

![Actual tutorial solution: surface and back temperatures, boundary mass-loss rate compared with gas generation, and profiles plotted at current cell centres.](../../assets/figures/08_quickstart.svg)

**Figure 8. Results of the complete tutorial in §2.10.** The 72-cell moving-mesh solve reaches 1,600 s successfully at an incident flux of 24.5325 kW/m². The solid green trace is the actual `solution.mass_loss_rate`; the dashed gas-generation trace is a reaction-rate proxy and can differ because gas is stored before it escapes. Both are displayed in g/(m²·s); this deck has constant cross section. Profiles use `solution.cell_centers[:, k]` at each saved time, and dots mark the current exposed-face temperature and position. The peak actual MLR is approximately 5.49 g/(m²·s). The editable Julia generator and the underlying CSV data are described in the Development page.

With `Plots` loaded, the package's recipe gives you quick standard views off the
solution directly:

```julia
using Plots

plot(solution)                      # surface-temperature trace (the default)
plot(solution; kind = :mlr)         # mass-loss-rate curve [kg/(m²·s)]
plot(solution; kind = :heatmap)     # T(cell, t) field as a heatmap
```

The recognised `kind` values are `:surface_T` (default), `:mlr`, and `:heatmap`
(time against cell index); any other value raises an informative error.

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
p3 = plot(xlabel = "z [mm] (0 = back face)", ylabel = "Temperature [K]",
          title = "Temperature profiles")
for t_target in (30, 120, 300, 600, 1200)
    k = argmin(abs.(times .- t_target))
    z_k  = get_cell_positions(solution; at = k)
    live = .!isnan.(z_k)                             # deactivated cells read NaN
    count(live) > 1 || continue                      # skip if one live cell remains
    plot!(p3, z_k[live] .* 1000, solution.T[live, k];
          lw = 1.5, label = "t = $(round(Int, times[k])) s")
end

fig = plot(p1, p2, p3; layout = (1, 3), size = (1200, 380))
savefig(fig, "cone_quickstart.png")
```

> **Note — true HRR needs a combustion-efficiency factor.** The cone HRR is
> usually `HRR = MLR × Δh_c,eff`, where the effective heat of combustion
> `Δh_c,eff` is the complete heat of combustion times the combustion
> efficiency. Pyrolysis.jl is a *condensed-phase* solver and outputs the
> mass-loss rate (and a gas-generation proxy); it does not model the gas-phase
> flame, so multiply by your `Δh_c,eff` to get an HRR estimate. See Technical
> Reference §1 for the scope boundary.

---

## 2.10 Complete script

The full tutorial in one runnable block uses the same material and atmospheric
exchange as `examples/06_cone_calorimeter.jl`. The example adds a ThermaKin
overlay and a longer analysis of the results.

```julia
using Pyrolysis
using Plots

# ---- Step 1: material -------------------------------------------------------
pmma_cp = CallableProperty(
    T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T),
)
pmma_k = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))

virgin = SolidComponent(:virgin, ρ = 1210.0, φ=.01, c = pmma_cp, k = pmma_k, ε = 0.96, α = 1.47,
                        λ = 1e-5)
char   = SolidComponent(:char,   ρ = 1210.0, φ=.01, c = pmma_cp, k = pmma_k, ε = 0.96, α = 1.47,
                        λ = 1e-5)
mma    = GaseousComponent(:MMA, M = 0.10012, c = 1500.0, k = 0.02)
air    = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026)

decomp = Reaction(:decomposition, 1 => (2, 3);
                  A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0, yields = (0.002, 0.998))

pmma = Material(
    name = :PMMA,
    components = (virgin, char, mma, air),
    reactions  = (decomp,),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# ---- Step 2: boundary conditions -------------------------------------------
q_external = 0.9813 * 25e3
h_conv     = 7.2
T_ambient  = 300.0

top_thermal =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +
    RadiativeBC(ε = 0.96, T_inf = 330.0) +
    ConvectiveBC(h = h_conv, T_inf = 330.0)
bottom_thermal =
    RadiativeBC(ε = 0.95, T_inf = 307.0) +
    ConvectiveBC(h = 4.0, T_inf = 307.0)

reservoir   = GasReservoir(P=101325., T=330., mole_fractions=(AIR=1.,))
g_film      = 0.05 * 101325. * 0.02897 / (8.31446261815324 * 330.)
top_mass    = AmbientGasBC(reservoir; film=MassTransferFilm(g_film))
bottom_mass = ImpermeableBC()

bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal, bottom = bottom_thermal),
    mass    = (top = top_mass,    bottom = bottom_mass),
)

# ---- Step 3: mesh -----------------------------------------------------------
sample_thickness = 5.8e-3
n_cells          = 72
mesh = create_uniform_mesh(sample_thickness, n_cells)

# ---- Step 4: problem --------------------------------------------------------
t_final = 1600.0
ξ_air = background_gas_concentration(pmma, 4, (1210., 0., 0., 0.), T_ambient; P=101325.)
problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = bc_set,
    T_initial = z -> T_ambient,
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0, z -> ξ_air],
    tspan = (0.0, t_final),
    dt_initial = 1e-3,
    output_times = collect(range(0.0, t_final, length = 601)),
)

# ---- Step 5: solve ----------------------------------------------------------
solution = solve(problem;
    abstol = 1e-11, reltol = 1e-9, pressure_tolerance = 1e-4, dtmin = 1e-14,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    min_thickness = 0.01e-3,
    deactivation = true,
    progress = true,
)

# ---- Step 6: extract --------------------------------------------------------
times     = solution.t
T_surface = solution.surface_T
mlr       = solution.mass_loss_rate .* 1000
gen_rate  = solution.extras.gas_generation_rate .* 1000
println("retcode = ", solution.retcode, ", outcome = ", solution.outcome)
println("peak T_s = ", round(maximum(T_surface), digits = 1), " K")
println("peak MLR = ", round(maximum(mlr), digits = 2), " g/(m²·s)")

# ---- Step 7: plot -----------------------------------------------------------
plot(solution; kind = :mlr)
```

---

## 2.11 What to do next

- **Switch to in-depth radiation.** Set `radiation_model = BEER_LAMBERT` (the
  solids already carry `α = 1.47` m²/kg) to model a semi-transparent polymer.
  Compare the surface-temperature trace to the surface-absorption run above.
  (Technical Reference §8.)
- **Add pressure-driven transport.** Give a component a permeability `κ > 0` to
  switch the solver into lumped gas transport mode for pressure-treated wood or charring
  systems. (UG 7; Technical Reference §9.)
- **Vary the cone flux.** Replace the scalar `q_external` with a function of
  time, `RadiativeFluxBC(flux = t -> ramp(t), absorptivity = 0.96)`; the solver
  re-samples the flux at every right-hand-side evaluation.
- **Run to burnout.** Extend `t_final` and read `solution.events`,
  `solution.extras.live_cells` and `solution.outcome`; UG 10 explains the
  threshold, the band, the mesh correction, and the envelope in which burnout is
  verified.
- **Sweep parameters.** Use `parameter_sweep` to run a grid of heat fluxes or
  kinetic parameters (UG 12), or finite differences of the same solves to
  estimate derivatives of an output with respect to `A`, `E`, `h`, or material
  properties (UG 13).
- **Audit the conservation.** Re-run with `rec = ConservationRecorder()` passed
  as `callback = rec`, then call `conservation_audit(solution, rec)` for the
  mass and energy closures over the solve (UG 9; Technical Reference §16).

You now have the full pipeline — material, boundary conditions, mesh, problem,
solve, solution — and can recombine these pieces for any of the worked cases in
UG 14.
