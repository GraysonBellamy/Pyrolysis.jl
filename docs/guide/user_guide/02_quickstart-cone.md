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
