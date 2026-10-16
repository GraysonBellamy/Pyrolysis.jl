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
> components in the material. The validator (§7.9) checks that `ξ_initial` has
> one closure per material component.

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

Several spatial-profile helpers are exported by the `Problem` submodule (import them
with `import Pyrolysis.Problem: linear_profile`, §5.10.3) for non-uniform initial states:
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

ALE adds work arrays, mesh-velocity integration, the dilation correction, and a
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

`WithChi` requires `use_ale = true`: only the ALE residual evolves `χ` and the
cross-section, so `solve` rejects a lateral shrinkage law on the fixed mesh with
an `ArgumentError`. The radial (cross-section) strain rate
`θ_A = law'(χ̄)/law(χ̄) · dχ̄/dt` is subtracted from each cell's volumetric
strain to get the axial component that drives mesh motion, so the two channels
(radial area change and axial node motion) do not double-count. With `WithChi`, the state vector carries the χ block as well:
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

## 7.7 The energy equation

The energy equation is the same in every mode; the four trait axes change how its
terms are computed, not which terms appear:

```
ρc_p ∂T/∂t = −∇·q + S_conv + Q_rad + Q_rxn
```

Here `ρc_p` sums every component, pore gas included; `−∇·q` is conduction;
`S_conv = −Σ_g c_{p,g} J_g ∂T/∂z` carries the sensible enthalpy of the gas flux
relative to the skeleton; `Q_rad` is the volumetric radiation source
(`BEER_LAMBERT` only); and `Q_rxn = −Σ_r h_r r_r` is the reaction heat. This is
the ThermaKin convention (Technical Reference §3 and §7).

**Reaction heats are DSC-style.** A heat `h` measured by DSC, or calibrated in
ThermaKin or Gpyro, is quoted at the reaction temperature, so it already contains
the gas products' sensible enthalpy at that temperature. The energy equation
therefore has no separate gas-generation enthalpy term. Such a term,
`−Σ_g Δh_g(T, T_ref)(∇·J_g)`, belongs only with heats referenced to the
`T_ref ≈ 298.15 K` ambient datum (formation-enthalpy style); with DSC-style heats
it would subtract the gas products' sensible enthalpy a second time (≈0.7–1.2 MJ
per kg of volatiles at 900 K, comparable to the heat of pyrolysis itself). If
your heats are referenced to the ambient datum, convert them to the reaction
temperature before use (Technical Reference §6).

### 7.7.1 No surface gas-enthalpy term

The boundary energy balance contains no gas-enthalpy term: the applied BC flux
equals the conducted flux (§5.4). A surface term `Σ_g ṁ_g Δh_g(T_s, T_∞)` would
overlap the always-active volumetric source by exactly `∫S_conv`, a spurious sink
of order 5–10 kW/m² that under-predicts surface temperature and mass-loss rate.
Use a plain convective + radiative boundary condition; the escaping gas
influences the boundary only through the blowing model of a `ConvectiveBC`
(§5.2.3). See Technical Reference §7 and §16 for the energy-balance accounting.

---

## 7.8 How modes map to the residual paths

The four trait axes fully determine the discrete RHS the solver integrates. The
two mesh-motion paths are the primary fork:

- **Eulerian path** — geometry is constant. The residual populates the typed
  state cache, updates properties, computes the physics RHS (dispatched on
  `Fickian` vs `DarcyFick`), and packs `dT` and `dξ`. No mesh velocity, no
  advection corrections.
- **ALE path** — geometry is recomputed from the node positions. The residual
  additionally updates the cross-section area (for `WithChi`), computes the
  per-cell strain `θ` and the column radial strain `θ_A`, integrates the
  cumulative mesh velocity `w` from the bottom up, applies the dilation
  correction to every component (every component rides the mesh, so there is
  no advection term), packs the `z` block, and (for `WithChi`) packs the χ
  derivatives.

Radiation, transport, and cross-section are handled by compile-time dispatch
*within* those two paths. The net effect is a small, type-stable set of
specialized residuals selected entirely by the mode bundle — you do not choose a
"path" directly; you choose the axes and the path follows. See Technical Reference
§15.1–§15.4 for the dispatch and the full ALE assembly sequence.

---

## 7.9 Validating a problem (`validate_problem`)

Before solving, run the validator to catch structural mistakes early:

```julia
Pyrolysis.Problem.validate_problem(problem)                 # throws on any error, returns true if OK
Pyrolysis.Problem.validate_problem(problem; strict = true)  # adds expensive checks (see below)
```

`solve` calls `validate_problem(problem; strict = false)` automatically unless you
pass `validate = false`. The default checks are:

- the mesh has cells, nodes, faces, cell states, and initialized geometry;
- the material has at least one component, and every reaction references valid
  component indices;
- the mesh has the expected top (id 1) and bottom (id 2) boundaries, and a
  per-species `ξ_inf` tuple on a `ConvectiveMassBC` has one entry per component;
- the time span is valid (`t_start < t_end`, `t_start ≥ 0`, positive `dt_initial`
  no larger than the span);
- the number of initial concentration closures equals the component count;
- a Darcy–Fick material does not start with (near-)zero porosity in a cell whose
  effective permeability is nonzero (§7.5.1).

It also warns, without failing, when a gas-producing Fickian material has a
condensed component with gas-transfer coefficient `λ = 0` (gas transport then
vanishes wherever the gas fraction vanishes; give such components a small `λ`,
e.g. `1e-6`). Errors are collected and thrown together in one message.

`strict = true` additionally verifies reaction stoichiometric mass balance
(product yields summing to 1.0), that the initial-condition closures are callable
at `z = 0`, and that all cell volumes are positive.

`Pyrolysis.Problem.is_valid_problem(problem)` returns a `Bool` without throwing.

---

## 7.10 Choosing modes for common scenarios

The table below is a quick recipe lookup. "auto" means the axis is inferred from
the material, not set on `solve`.

| Scenario | `use_ale` | `radiation_model` | transport (auto) | cross-section (auto) |
|---|---|---|---|---|
| TGA / DSC-like (thin, no recession) | `false` | `NO_RADIATION` | Fickian (`κ = 0`) | NoChi (no law) |
| Cone calorimeter, opaque char | `true` | `SURFACE_ABSORPTION` | Fickian or Darcy | NoChi or WithChi |
| Cone calorimeter, semi-transparent polymer | `true` | `BEER_LAMBERT` (set `α`) | Fickian | NoChi |
| Intumescent coating (swelling) | `true` | `SURFACE_ABSORPTION` | Darcy (`κ > 0`) | WithChi (swell law) |
| Pressure-driven char / pressure-treated wood | `true` | `SURFACE_ABSORPTION` | Darcy (`κ > 0`) | NoChi or WithChi |

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

sol = solve(problem)        # Eulerian, NoRadiation, Fickian, NoChi
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
# transport = DarcyFick (auto, from κ); cross-section = NoChi (no law)
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
- The active physics is the **simulation mode**: four trait axes.
  - *Mesh motion* (`use_ale`, default `false`) and *radiation*
    (`radiation_model`, default `NO_RADIATION`) are chosen at solve time.
  - *Gas transport* (Fickian vs. Darcy–Fick) and *cross-section* (NoChi vs.
    WithChi) are inferred automatically from the material — from component
    permeability `κ` and from `lateral_shrinkage_law`, respectively.
- There is one energy equation, with gas enthalpy carried by `S_conv` and
  DSC-style reaction heats; the boundary balance has no surface gas-enthalpy
  term, which would double-count `S_conv`.
- `validate_problem` (run automatically by `solve` unless `validate = false`)
  catches structural errors before integration.

The next chapter, "Running the Solver," covers the rest of the `solve` keyword
surface: integrator and tolerance selection, Jacobian backends and linear
solvers, callbacks, warm starts, and output control.
