# 16. API and Parameter Reference

This chapter is the consolidated reference for the public API of Pyrolysis.jl:
every name `using Pyrolysis` exports, with its signature, keywords and defaults;
the keyword tables of `solve`, `setup` and `integrate`; the fields of the
solution; and the defaults and constants built into the package. Earlier
chapters teach the workflow; this one is the lookup table you return to for the
exact keyword, the default value, or the module a name lives in. Every
signature and default here is read from the source. A name that `using
Pyrolysis` does not bring into scope is marked with the submodule that holds
it, so a copied snippet never raises an `UndefVarError`. For the physics behind
a quantity, follow the Technical Reference pointers.

The whole API assumes four conventions (Technical Reference §2 and §3):

- `z = 0` is the **back face** (boundary tag `:bottom`), `z = L` the **exposed
  surface** (tag `:top`). Heat enters from the top.
- The fluxes a boundary condition applies are **positive into the domain**. The
  solution reports the mass leaving the sample as a positive outflow
  (`mass_loss_rate`, `boundary_mass_flux`).
- Concentrations are **mass concentrations** `ξ_j` in kg·m⁻³, not mass
  fractions. The solver's states are the component masses per cell
  `m_ij = ξ_j V_i` [kg]; initial conditions and boundary data are given as
  concentrations.
- The heat of reaction is stored **`h > 0` endothermic** (absorbs heat),
  `h < 0` exothermic (Technical Reference §6).

---

## 16.1 What `using Pyrolysis` gives you, and what it does not

The top-level `Pyrolysis` module exports a curated surface of 99 names. The
package is organized in submodules, and most of their names are not
re-exported at the top level. There are three ways to reach a name:

1. **Top-level export**, available after `using Pyrolysis`. This is the stable
   public API; prefer it in production code.
2. **Submodule**, for example `Pyrolysis.Meshes.n_cells`, or imported by name:
   `using Pyrolysis.Meshes: n_cells`.
3. **`Pyrolysis.Internal`**, an umbrella namespace that re-exports every name
   the submodules export, in one place:

   ```julia
   using Pyrolysis            # the curated top-level surface
   using Pyrolysis.Internal   # every name the submodules export
   ```

   A name that two submodules export as different objects is left out of
   `Internal`; reach it through its submodule. Names reached through
   `Internal` are not part of the stable API and may move or change between
   minor releases (Chapter 11, §11.9).

The 99 top-level exports, by group:

| Group | Exports | Section |
|---|---|---|
| Components and materials | `SolidComponent`, `LiquidComponent`, `GaseousComponent`, `Material`, `phase`, `molar_mass`, `is_gas`, `is_solid`, `is_liquid` | 16.2.1, 16.2.5 |
| Reactions | `Reaction`, `DepletionLimiter` | 16.2.3, 16.2.5 |
| Mixing | `MixingRule`, `PARALLEL`, `SERIES`, `WEIGHTED`, `BRUGGEMAN`, `CARMAN_KOZENY`, `MixingSpec`, `MaterialMixing` | 16.2.4 |
| Property functions | `ConstantProperty`, `LinearProperty`, `PolynomialProperty`, `SutherlandProperty`, `ArrheniusProperty`, `TableProperty`, `CallableProperty`, `StateDependentProperty` | 16.2.2 |
| Radiation models | `RadiationModel`, `NO_RADIATION`, `SURFACE_ABSORPTION`, `BEER_LAMBERT` | 16.2.6 |
| Thermal conditions | `AdiabaticBC`, `HeatFluxBC`, `ConvectiveBC`, `RadiativeBC`, `RadiativeFluxBC`, `DirichletThermalBC`, `CombinedThermalBC` | 16.4.1 |
| Mass conditions | `ImpermeableBC`, `MassFluxBC`, `AmbientGasBC`, `GasReservoir`, `MassTransferFilm` | 16.4.2 |
| Pressure conditions | `NeumannPressureBC`, `ConvectivePressureBC` | 16.4.3 |
| Condition set and pore air | `BoundaryConditionSet`, `background_gas_concentration` | 16.4.5, 16.4.6 |
| Blowing models and gas data | `AbstractBlowingModel`, `NoBlowing`, `CouetteFilm`, `PropertyRatioFilm`, `PoolingFilm`, `FunctionBlowing`, `GasSpecies`, `NITROGEN`, `nasa7_heat_capacity`, `eucken_viscosity`, `fuller_diffusivity`, `BlowingFace`, `blowing_state`, `blowing_ratio`, `blowing_factor` | 16.4.7 |
| Meshes | `create_uniform_mesh`, `create_stretched_mesh`, `create_mesh_from_nodes`, `Planar`, `Cylindrical`, `Spherical`, `MeshTopology`, `MeshGeometry` | 16.3 |
| Problem and solve | `PyrolysisProblem`, `solve`, `setup`, `integrate`, `CoreProblem`, `parameter_sweep` | 16.5.1, 16.7, 16.9 |
| Correction | `NoCorrection`, `EquidistributionCorrection`, `ActivationTime` | 16.11 |
| Integrator and linear solver | `default_policy`, `rodas5p_policy`, `fbdf_policy`, `kencarp4_policy`, `StructuredKLU` | 16.8 |
| Solution and accessors | `PyrolysisSolution`, `SolverStatistics`, `get_times`, `get_temperatures`, `get_concentrations`, `get_gas_generation_rate`, `n_timesteps`, `succeeded`, `completed`, `get_profile`, `get_cell_positions`, `interpolate_solution` | 16.5 |
| Conservation | `ConservationRecorder`, `conservation_audit` | 16.12 |

**Names that are not top-level exports.** The names below are often needed
and live in submodules (all of them are also in `Pyrolysis.Internal`). Import
them by name, for example `using Pyrolysis.Problem: validate_problem`.

| Names | Submodule |
|---|---|
| `validate_problem`, `is_valid_problem`, `initial_fields`, `step_profile`, `linear_profile`, `gaussian_profile`, `tanh_profile` | `Pyrolysis.Problem` |
| `n_cells`, `n_nodes`, `n_faces`, `boundary_tags`, `boundary_faces`, `is_fixed`, `geometry!`, `line_mesh`, `box_mesh`, `cell_mesh` | `Pyrolysis.Meshes` |
| `n_components`, `n_reactions`, `R_GAS`, `STEFAN_BOLTZMANN`, `T_REF`, `P_REF`, `T_PROPERTY_MIN`, `T_PROPERTY_MAX`, `RHO_CP_FLOOR`, `Phase` (`SOLID`, `LIQUID`, `GAS`) | `Pyrolysis.Materials` |
| `SolutionExtras`, `ledger_rates`, `live_inventories` | `Pyrolysis.Outputs` |
| `DeactivationEvent`, `CellDeactivation`, `replay_events!`, `top_live_cell`, `threshold_margin` | `Pyrolysis.Deactivation` |
| `residual!`, `CoreDefinition`, `CoreWorkspace`, `NoSources`, `VolumetricSources` | `Pyrolysis.Assembly` |
| `jac_prototype`, `jacobian!`, `jacobian_pattern`, `color_count`, `dense_jacobian`, `missing_entries`, `compare_entries` | `Pyrolysis.ColoredJacobian` |
| `StateLayout`, `state_length`, `mass_index`, `energy_index`, `position_index`, `velocity_index`, `intensity_index`, `is_differential`, `mass_matrix` | `Pyrolysis.States` |
| `IntegratorPolicy`, `DomainCheck`, `consistent_state!`, `default_linear_solver`, `DT_AFTER_EVENT`, `CORRECTION_BAND`, `MAX_CONSTRAINT_REJECTIONS`, `FAILURE_MEMORY_GROWTH`, `FAILURE_MEMORY_STEPS` | `Pyrolysis.Integration` |
| `PrescribedCorrection`, `residues`, `vanishes`, `EQUIDISTRIBUTION_TIME`, `ACTIVATION_TIME_FACTOR` | `Pyrolysis.Closures` |

`to_property` is not exported by any submodule; call it as
`Pyrolysis.Materials.to_property`.

### 16.1.1 Submodule map

The package is a layered set of submodules, loaded in the order of the table,
each using only those before it. Everything a submodule exports is reachable as
`Pyrolysis.<Submodule>.X` and is bundled into `Pyrolysis.Internal`.

| Submodule | Owns |
|---|---|
| `Materials` | components, `Material`, `Reaction`, `DepletionLimiter`, the property-function types, the mixing rules, the blowing-model types and gas data, the `RadiationModel` enum, the physical constants |
| `Properties` | the mixture rules (volume fractions, density, skeletal volume) and the effective properties (conductivity, emissivity, absorptivity, absorption, gas transfer, permeability, gas diffusivity) |
| `BoundaryConditions` | the thermal, mass and pressure condition types, `BoundaryConditionSet`, the face-local surface energy balance, the blowing-model evaluation (`BlowingFace`, `blowing_state`), the mass-condition fluxes |
| `Physics` | reaction kinetics, the diffusive gas flux at a face, volume change, the ideal-gas pressure, `background_gas_concentration` |
| `Meshes` | the face-based mesh (`Mesh`, `MeshTopology`, `MeshGeometry`, `geometry!`), the coordinate systems, the line constructors, and the meshes of two and three dimensions (`box_mesh`, `cell_mesh`) |
| `States` | the state layout and its index functions, the energy states, the caloric laws, the per-cell reference data and tolerances, the state mapping (`derive!`) |
| `Kernels` | the physics kernels: reactions, conduction, gas transport, in-depth radiation, the thermal boundary, the relative fluxes of a correction |
| `Closures` | the mesh motion: the balance rows of a line, the multiplier closure of other meshes, the corrections, the cross-section rows |
| `Assembly` | the residual `residual!(du, u, def, t, ws)`, its definition and buffers, the external sources |
| `Deactivation` | cell deactivation: the threshold test, the merge, the events and their replay |
| `ColoredJacobian` | the sparsity pattern, the coloring, the colored Jacobian, the dense reference |
| `Problem` | `PyrolysisProblem`, the initial fields and profile helpers, `validate_problem` |
| `Outputs` | `PyrolysisSolution`, `SolutionExtras`, `SolverStatistics`, the accessors, the conservation ledger, recorder and audit |
| `Integration` | `setup`, `integrate`, `solve`, `CoreProblem`, the integrator policies, `StructuredKLU`, the domain test, `parameter_sweep` |

### 16.1.2 Extensions

One package extension ships. Loading `Plots` (or `RecipesBase`) attaches a plot
recipe for `PyrolysisSolution`:

```julia
using Pyrolysis, Plots
plot(sol)                    # surface temperature, the default kind = :surface_T
plot(sol; kind = :mlr)       # sol.mass_loss_rate against time
plot(sol; kind = :heatmap)   # sol.T over time and cell index
```

Any other `kind` raises an `ArgumentError` naming the three (Chapter 9, §9.6).
Linear solvers need no extension: `LinearSolve` is a dependency, and any of its
sparse factorizations is passed as the `linear_solver` keyword or to a policy
constructor (§16.8). Sensitivities need none either: they are finite
differences of whole solves (§16.10).

---

## 16.2 Materials, components, and property functions

See Technical Reference §4 (materials and property functions) and §5
(effective properties). Every name in this section is a top-level export unless
stated otherwise.

### 16.2.1 Components

```julia
SolidComponent(name::Symbol; ρ, c, k, ε=.9, absorptivity=nothing,
    α=0., γ=1., λ=nothing, κ=0., ρ_intrinsic=nothing, φ=nothing,
    validate_properties=true)
LiquidComponent(name::Symbol; ρ, c, k, ε=.95, absorptivity=nothing,
    α=0., γ=0., λ=nothing, κ=0., ρ_intrinsic=nothing, φ=nothing,
    validate_properties=true)
GaseousComponent(name::Symbol; M, c, k, ε=0., α=0., γ=0.,
    validate_properties=true)
```

Density is in kg/m³, heat capacity in J/(kg K), conductivity in W/(m K),
molar mass in kg/mol, condensed `λ` in m²/s and permeability in m².
Participating condensed components require explicit `λ` for `LumpedTransfer`;
zero is a deliberate closed wall channel. `γ=0` condensed components may
have neither nonzero `λ` nor nonzero `κ`.

Specify either intrinsic density or a pure-phase porosity, with
`ρ_intrinsic=ρ/(1-φ)` for constant bulk density. Every condensed component
occupies skeletal volume regardless of `γ`. Gas-bearing initialized and
accepted states need resolved positive pore volume. A gas component has no
matrix transport coefficient. Its `γ` may affect volume change but never
matrix transport weights.

Numeric properties become `ConstantProperty`; tuples, callables and explicit
property objects follow the conversion rules in §16.2.2. The phase predicates
`is_solid`, `is_liquid`, `is_gas`, `phase` and `molar_mass` are exported.

---

### 16.2.2 Property-function types

Every property field takes a number or one of these callable objects. They
implement `p(T)`, and `p(T, ξ)` for the state-dependent one; the Jacobian
differentiates through them by forward-mode automatic differentiation
(Technical Reference §4.5).

| Type | Constructor | Value `p(T)` |
|---|---|---|
| `ConstantProperty` | `ConstantProperty(value)` | `value` |
| `LinearProperty` | `LinearProperty(a, b)` | `a + b·T` |
| `PolynomialProperty` | `PolynomialProperty((c0, c1, …))` | `evalpoly(T, coeffs)`, lowest order first |
| `ArrheniusProperty` | `ArrheniusProperty(A, E, T_ref = 298.15)` | `A·exp(−E/(R_g T))` (`T_ref` is stored, not used in the value) |
| `TableProperty` | `TableProperty((T1, …), (v1, …))` | linear interpolation; the end values are held outside the table |
| `CallableProperty` | `CallableProperty(f)` | `f(T)` |
| `SutherlandProperty` | `SutherlandProperty(μ_ref, T_ref, S)`; `SutherlandProperty()` is air (`1.716e-5` Pa·s at `273.15` K, `S = 110.4` K) | `μ_ref (T/T_ref)^1.5 (T_ref + S)/(T + S)` |
| `StateDependentProperty` | `StateDependentProperty(f)` | `f(T, ξ)`, `ξ` the concentration tuple (calling it as `p(T)` is an error) |

Every type except `ConstantProperty` evaluates its law at `T` clamped to the
property bracket `[200, 3000]` K (`T_PROPERTY_MIN`, `T_PROPERTY_MAX` in
`Pyrolysis.Materials`): inside the bracket the law sees `T` unchanged, outside
it the nearest end, so a trial temperature far from the physical one never
reaches a fitted law's extrapolation. `TableProperty` requires strictly
increasing temperatures (checked at construction). `StateDependentProperty` is
the type for a composition-dependent quantity such as a moisture-dependent heat
of sorption.

```julia
# Conductivity rising with temperature:
char = SolidComponent(:char, ρ = 150, c = 1500, k = LinearProperty(0.10, 1.5e-4))

# Tabulated heat capacity:
cp_tab = TableProperty((300.0, 500.0, 800.0), (1500.0, 1800.0, 2100.0))
virgin = SolidComponent(:virgin, ρ = 530, c = cp_tab, k = 0.12)
```

### 16.2.3 Reactions

`Reaction` is a single-reactant Arrhenius reaction with one or two products.
The component references are indices or names; a name is resolved against the
component names when the `Material` is built.

```julia
# one product
Reaction(name, reactant::Int => product::Int;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0, validate_mass = true)
Reaction(name, reactant::Symbol => product::Symbol; A, E, h, n = 1.0, …)

# two products
Reaction(name, reactant::Int => (p1::Int, p2::Int);
         A, E, h, yields, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0,
         validate_mass = true)
Reaction(name, reactant::Symbol => (p1::Symbol, p2::Symbol); A, E, h, yields, n = 1.0, …)
```

| Keyword | Meaning (symbol) | Units | Default |
|---|---|---|---|
| `A` | pre-exponential factor `A_i` | s⁻¹ (first order) | required |
| `E` | activation energy `E_i` | J·mol⁻¹ | required |
| `h` | heat of reaction, **`h > 0` endothermic**; number, `(a, b)` tuple (`a + b·T`) or property object | J per kg of reactant | required |
| `yields` | mass yields `(y1, y2)` of the two products | – | required (two products) |
| `n` | reaction order, `n ≥ 0` | – | `1.0` |
| `T_min`, `T_max` | smooth lower and upper temperature gates | K | `0.0`, `Inf` |
| `gate_width` | width of the tanh gates, positive | K | `1.0` |
| `validate_mass` | require the yields to sum to one | – | `true` |

The reactant is consumed with mass stoichiometry one; the product yields must
sum to one (the one-product form has yield one). Every yield must be
non-negative and the order non-negative, enforced at construction whatever
`validate_mass` says; `validate_mass = false` skips only the sum check. Use a
`StateDependentProperty` for `h` when it depends on composition. Technical
Reference §6 gives the rate law, the gates and the depletion limiter.

```julia
# endothermic charring: virgin -> 0.2 char + 0.8 gas
rxn = Reaction(:charring, :virgin => (:char, :gas);
               A = 1e15, E = 200e3, h = 500e3, yields = (0.2, 0.8),
               T_min = 473.0, T_max = 1273.0)

# exothermic step (negative h releases heat)
rxn2 = Reaction(:char_oxidation, 1 => 2; A = 1e10, E = 160e3, h = -25e6)
```

Reactions are passed to `Material` as a tuple.

### 16.2.4 Mixing rules: `MixingRule`, `MixingSpec`, `MaterialMixing`

The effective conductivity `k`, gas-transfer coefficient `λ` and permeability
`κ` are each mixed by their own rule (Technical Reference §5):

```julia
@enum MixingRule PARALLEL SERIES WEIGHTED BRUGGEMAN CARMAN_KOZENY
```

| Value | Formula and use |
|---|---|
| `PARALLEL` | `k = Σ v_j k_j` (upper bound) |
| `SERIES` | `1/k = Σ v_j/k_j` (lower bound) |
| `WEIGHTED` | `k = β k_parallel + (1 − β) k_series` |
| `BRUGGEMAN` | effective-medium theory; the `k` channel only |
| `CARMAN_KOZENY` | porosity-based permeability `κ = min(ℓ² φ³/(180 (1 − φ)²), ℓ²)`; the `κ` channel only |

`MixingSpec` bundles a rule with its weight and, for Carman–Kozeny, the length:

```julia
MixingSpec(rule::MixingRule, β::Real = 0.5; length_scale::Real = NaN)
```

`β` lies in `[0, 1]`. `CARMAN_KOZENY` requires a positive `length_scale`, the
characteristic pore or particle size `ℓ` in metres, given as the keyword.
`MaterialMixing` builds the `(k, λ, κ)` NamedTuple `Material` expects; each
channel takes a `MixingSpec` or a bare `MixingRule` (then `β = 0.5`):

```julia
MaterialMixing(; k = MixingSpec(WEIGHTED, 0.5),
                 λ = MixingSpec(PARALLEL, 0.5),
                 κ = MixingSpec(WEIGHTED, 0.5))
```

```julia
mixing = MaterialMixing(
    k = MixingSpec(WEIGHTED, 0.6),
    λ = PARALLEL,                                  # bare rule: β = 0.5
    κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),
)
```

### 16.2.5 `Material` and `DepletionLimiter`

```julia
Material(; name::Symbol,
           components::Tuple,
           reactions::Tuple = (),
           mixing = MaterialMixing(),
           lateral_shrinkage_law = nothing,
           depletion_limiter = DepletionLimiter(),
           gas_viscosity = SutherlandProperty(),
           gas_transport = LumpedTransfer())
```

| Keyword | Meaning | Default |
|---|---|---|
| `name` | material identifier | required |
| `components` | tuple of components; the order defines the index `j` | required |
| `reactions` | tuple of `Reaction`s (names resolved here) | `()` |
| `mixing` | the `(k, λ, κ)` NamedTuple of `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `χ̄ ↦ A/A₀`, a function of the column-average conversion; `nothing` keeps the cross section fixed | `nothing` |
| `depletion_limiter` | the limiter applied to every reaction rate of the material | `DepletionLimiter()` |
| `gas_transport` | declared pore-gas model and endpoint domain | `LumpedTransfer()` for gas-bearing materials |
| `gas_viscosity` | pore-gas viscosity `μ(T)` [Pa·s] in the Darcy mobility `κ/μ` (both gas models); number or property, positive and finite at sample temperatures from 200 to 1800 K | `SutherlandProperty()` (air) |

A `lateral_shrinkage_law` adds the variable cross section: the
conversion-history pair and the area with its lateral rate join the state.
It needs `use_ale = true` and planar coordinates, and `setup` refuses it
otherwise (§16.6; Technical Reference §10). The law must be positive and
finite on `[−0.001, 1.001]` and equal one at zero conversion, for example
`χ̄ -> 1 - 0.3χ̄` for 30 % lateral shrinkage at full conversion.

```julia
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, φ=.3, c = 1420, k = 0.21, ε = 0.86, λ = 1e-6),
        GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02,),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3),
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)
```

```julia
DepletionLimiter(; threshold = 1.0, join = 1.0, base_scale = 1e-8, enabled = true)
```

The depletion limiter multiplies each reactant's rate by a factor that is
exactly zero at and below zero concentration and rises to one above
`threshold` [kg·m⁻³]: a cubic onset from zero up to `join · threshold`, joined
with value and slope to `tanh(ξ/threshold)` above. The rate is `C¹` through zero
for every reaction order, so an absent or negative trial concentration is never
consumed. The order's power is taken of `base_scale · softplus(ξ/base_scale)`,
which equals `ξ` to within `base_scale` [kg·m⁻³]. `threshold`, `join` and
`base_scale` must be positive. With `enabled = false` the factor is one and
only the finite power base remains, so a zero-order reaction then runs at its
full rate at any concentration (Chapter 4, §4.10; Technical Reference §6).

> Pyrolysis.jl has no copy-with-overrides constructor: to vary a parameter (in a
> sweep, for example), rebuild the `Material` and the `PyrolysisProblem` from
> their parts (§16.9).

### 16.2.6 The `RadiationModel` enum

```julia
@enum RadiationModel NO_RADIATION SURFACE_ABSORPTION BEER_LAMBERT
```

| Value | Behaviour |
|---|---|
| `NO_RADIATION` | no in-depth absorption; the absorbed incident flux `a·q` of a `RadiativeFluxBC` enters the surface energy balance |
| `SURFACE_ABSORPTION` | the same treatment as `NO_RADIATION` (an opaque material); the component absorption coefficients are not read |
| `BEER_LAMBERT` | exponential in-depth absorption through the components' `α`; the face intensities join the state as algebraic rows; planar coordinates only |

The model is chosen by `solve`'s `radiation_model` keyword (§16.7; Chapter 8,
§8.8). `NO_RADIATION` and `SURFACE_ABSORPTION` give identical residuals; name
`SURFACE_ABSORPTION` when the incident flux is the physics you mean. The
re-radiation (`RadiativeBC`) and convection (`ConvectiveBC`) losses apply under
every model.

---

## 16.3 Meshes and geometry

See Chapter 6 and Technical Reference §13. The mesh constructors and the
coordinate systems are top-level exports; the mesh's own functions are in
`Pyrolysis.Meshes`.

### 16.3.1 Mesh constructors

```julia
create_uniform_mesh(L::Real, n_cells::Integer;
                    coordinates = Planar(), cross_section = <default>)

create_stretched_mesh(L::Real, n_cells::Integer, stretch::Real;
                      stretch_from = :top, coordinates = Planar(),
                      cross_section = <default>)

create_mesh_from_nodes(z::AbstractVector{<:Real};
                       coordinates = Planar(), cross_section = <default>)
```

| Argument or keyword | Meaning | Default |
|---|---|---|
| `L` | domain thickness, or outer radius in cylindrical and spherical coordinates [m] | required |
| `n_cells` | number of cells, at least one | required |
| `stretch` | ratio of successive cell sizes away from the `stretch_from` boundary, `≥ 1` (1 is uniform; a smaller value is refused) | required |
| `stretch_from` | `:top` (finest cell at the exposed surface) or `:bottom` (finest at the back) | `:top` |
| `z` | node positions, starting at zero and strictly increasing [m] | required |
| `coordinates` | `Planar()`, `Cylindrical()` or `Spherical()`; positions are radii in the last two | `Planar()` |
| `cross_section` | the cross-section measure `A₀`: the face area of a planar column [m²], the axial length of a cylinder [m], the solid angle of a sphere | `1.0`, `1.0`, `4π` |

Each constructor returns a `Mesh{1}` (type `Pyrolysis.Meshes.Mesh`) with two
fields: `topology`, a `MeshTopology`, and `positions`, the node positions (one
`SVector{1}` per node). Cell 1 lies at the back; the face at `z = 0` carries the
tag `:bottom` and the face at `z = L` the tag `:top`; the back node is pinned.
A mesh holds no component count and no initial state: `setup` evaluates the
problem's `T_initial` and `ξ_initial` at the cell centres (§16.4.6). `solve`
never modifies the mesh, so one mesh can serve many problems.

```julia
mesh    = create_uniform_mesh(0.01, 100)                            # 1 cm, 100 cells
meshst  = create_stretched_mesh(0.02, 80, 1.03; stretch_from = :top)
rod     = create_uniform_mesh(0.005, 40; coordinates = Cylindrical())  # radius 5 mm
layered = create_mesh_from_nodes([0.0, 0.002, 0.004, 0.006, 0.007, 0.008])
sample  = create_uniform_mesh(0.01, 50; cross_section = 0.01)        # a 100 cm² face
```

### 16.3.2 Mesh types and accessors

| Name | What it is |
|---|---|
| `Planar`, `Cylindrical`, `Spherical` | the coordinate systems, constructed without arguments |
| `MeshTopology` | the immutable connectivity: fields `coordinates`, `faces`, `cells`, `n_nodes`, `boundary_tags` (tag → face ids), `pinned` (nodes at rest), `lagrangian` (nodes that always move with the solid), `cross_section` |
| `MeshGeometry` | the metrics at one set of node positions: `node_position`, `face_area`, `face_normal`, `face_center`, `face_sweep`, `cell_volume`, `cell_center`, `cell_node_volume`; built empty by `MeshGeometry{T}(topology)` and filled by `Pyrolysis.Meshes.geometry!(g, x, cross_section, topology)`, its only writer |

The accessors take the topology (`mesh.topology`) and live in
`Pyrolysis.Meshes`: `n_cells`, `n_nodes`, `n_faces`, `boundary_tags` (sorted),
`boundary_faces(topology, tag)` (an `ArgumentError` for an unknown tag) and
`is_fixed` (every node pinned). The cell centres and volumes at the reference
positions are `cell_centers(mesh)` and `cell_volumes(mesh)` of
`Pyrolysis.Problem` (imported by name; they are not exported).

```julia
using Pyrolysis.Meshes: n_cells, n_nodes, boundary_faces, MeshGeometry, geometry!

topo = mesh.topology
n_cells(topo), n_nodes(topo)          # (100, 101)
boundary_faces(topo, :top)            # [101], the face at z = L
g = MeshGeometry{Float64}(topo)
geometry!(g, first.(mesh.positions), topo.cross_section, topo)
g.cell_volume                         # V_i [m³]
```

`Pyrolysis.Meshes` also builds meshes of quadrilaterals and hexahedra
(`box_mesh`, `cell_mesh`). `setup` accepts them at the residual level, but a
full solve on them is not a supported user path, and gas transport, in-depth
radiation, the variable cross section, the correction and deactivation are
refused there (Chapter 6, §6.8). The component count of a material is
`n_components(material)` (`Pyrolysis.Materials`).

### 16.3.3 The state vector

Most users never index the raw state: `solve` returns named fields (§16.5). For
a custom loop or a callback that reads the integrator's state, the state is a
sequence of node and cell blocks, on a line node 1, cell 1, node 2, …, cell `n`,
node `n + 1`, so the Jacobian is banded:

| Block | Holds | Present when |
|---|---|---|
| node `k` | position `x_k` [m] and velocity `w_k` [m/s] | moving mesh (`use_ale = true`) |
| node `k` | intensity `I_k` of the face at the node [W/m²] | `radiation_model = BEER_LAMBERT` |
| cell `i` | masses `m_i1 … m_i,NC` [kg] and the energy state (`T_i` [K], or `H_i` [J] under `energy_state = :enthalpy`) | always |
| cell `i` | conversion-history pair `G_i`, `M_i` [kg] | variable cross section |
| global | cross-section area `A` [m²] and lateral strain rate `θ_A` [1/s] | variable cross section |

The node velocities, the face intensities, `A` and `θ_A` are algebraic rows;
the others are differential. Read a raw state only through the index functions
of `Pyrolysis.States` on the layout of a set-up problem, `cp.def.layout`
(`state_length`, `mass_index(layout, i, j)`, `energy_index(layout, i)`,
`position_index(layout, k)`, `velocity_index(layout, k)`,
`intensity_index(layout, f)`, `history_index(layout, i, a)`,
`area_index(layout)`, `lateral_rate_index(layout)`, `is_differential(layout)`,
`mass_matrix(layout)`). Chapter 6 (§6.6), Chapter 11 (§11.9.2) and Technical
Reference §15.1 describe the layout.

---

## 16.4 Boundary and initial conditions

See Chapter 5 and Technical Reference §12. Every condition type in this section
and `BoundaryConditionSet` are top-level exports.

### 16.4.1 Thermal conditions

| Type | Constructor | Flux `q_BC` (positive into the domain) |
|---|---|---|
| `AdiabaticBC` | `AdiabaticBC()` | `0` |
| `HeatFluxBC` | `HeatFluxBC(flux)`, a number or `t -> q` | `flux` |
| `ConvectiveBC` | `ConvectiveBC(; h, T_inf, blowing = NoBlowing(), validate = true)`, or `ConvectiveBC(h, T_inf)` | `R_b h (T_∞ − T_s)`, `R_b` the blowing model's ratio (one without outflow) |
| `RadiativeBC` | `RadiativeBC(; ε = nothing, T_inf, F = 1.0, validate = true)` | `F ε σ (T_∞⁴ − T_s⁴)`; `ε = nothing` is the exposed cell's effective emissivity |
| `RadiativeFluxBC` | `RadiativeFluxBC(flux; absorptivity = nothing, validate = true)`, or `RadiativeFluxBC(; flux, absorptivity = nothing, validate = true)` | `a q(t)`, incident external radiation; `absorptivity = nothing` is the exposed cell's effective absorptivity |
| `DirichletThermalBC` | `DirichletThermalBC(T)`, a number or `t -> T` | holds `T_s = T` |
| `CombinedThermalBC` | `bc1 + bc2 + …` | the sum of its members |

Notes:

- `h`, `T_inf`, `F`, `flux` and `absorptivity` take a number or a function of
  `t`; `RadiativeBC`'s `ε` takes a number or a function of the surface
  temperature. With `validate = true` a constant value is checked: `h ≥ 0`,
  `T_inf > 0`, and `ε`, `F` and `absorptivity` in `[0, 1]`.
- `blowing` is a blowing model (§16.4.7) that reduces `h` with the outward gas
  flux at the face.
- The effective absorptivity is the composition-weighted mean of the
  components' `absorptivity` where set and their emissivity otherwise
  (Kirchhoff), so with no component absorptivity it equals the effective
  emissivity.
- How a `RadiativeFluxBC` acts depends on `radiation_model`: under
  `NO_RADIATION` and `SURFACE_ABSORPTION` the absorbed flux `a q` enters the
  surface balance; under `BEER_LAMBERT` the condition contributes zero to the
  surface balance and the absorbed flux enters the first live cell as the
  intensity of its exposed face, through the `:top` faces only (a
  `RadiativeFluxBC` on another tag is refused). In a `CombinedThermalBC` under
  `BEER_LAMBERT` only the first `RadiativeFluxBC` member feeds the beam; a
  second radiative member is dropped silently, a limitation to keep in mind
  when combining two heaters.
- An adiabatic face, a `HeatFluxBC` alone, or a `DirichletThermalBC` alone, is
  applied in closed form. Any other condition sets the surface temperature `T_s` by a Newton
  solve of the surface energy balance `q_BC(T_s) = k_eff (T_s − T_cell)/d`
  (tolerance `1e-10`, at most 50 iterations, the condition evaluated at the
  iterate smoothly clamped into `[200, 3000]` K); inside a combination a
  `DirichletThermalBC` enters as a conductance of `1e8` W·m⁻²·K⁻¹ to `T`
  (Chapter 5, §5.4; Technical Reference §12.2). The flux of a time-dependent
  condition is evaluated at every right-hand-side evaluation.
- Combine conditions with `+`:
  `RadiativeFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)`.

### 16.4.2 Mass and reservoir conditions

| Type | Construction | Meaning |
|---|---|---|
| `ImpermeableBC` | `ImpermeableBC()` | zero species flow |
| `GasReservoir` | `GasReservoir(P=...,T=...,mole_fractions=(air=1.,))` | complete external ideal-gas state |
| `AmbientGasBC` | `AmbientGasBC(reservoir;film=nothing)` | pore/reservoir exchange in either direction |
| `MassTransferFilm` | `MassTransferFilm(g)` | isobaric common-diffusivity conductance [kg/(m² s)] |
| `MassFluxBC` | `MassFluxBC(fluxes=(air=-1e-4,),inflow_temperature=(air=300.,))` | prescribed outward species flux and incoming energy |

Mass conditions are complete closures and do not add with `+`. Reservoir
fractions are normalized named gas fractions; explicit vacuum may omit them.

---

### 16.4.3 Pressure conditions

| Type | Construction | Meaning |
|---|---|---|
| `ConvectivePressureBC` | `ConvectivePressureBC(reservoir;resistance=...,flow=:mass,film=nothing)` | pressure drop proportional to selected flow |

The hydraulic flow selector is `:darcy_velocity`, `:mass` or `:molar`.
Resistance units match that selector and values must be finite/nonnegative.
These conditions work with either material gas model.

---

### 16.4.4 Mesh motion

There are no mesh-motion conditions. With `use_ale = true` the pinned back node
stays at rest and every other node moves with the solid velocity plus a
correction, `w = v_s + c`, the velocities being algebraic rows of the state
(§16.7; Chapter 8, §8.7; Technical Reference §11).

### 16.4.5 `BoundaryConditionSet`

```julia
BoundaryConditionSet(; thermal::NamedTuple = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                       mass::NamedTuple    = (top = ImpermeableBC(), bottom = ImpermeableBC()))
```

The conditions are keyed by the mesh's boundary tags. A partial NamedTuple is
merged with the defaults, so `thermal = (top = HeatFluxBC(50e3),)` leaves the
bottom adiabatic. `validate_problem` requires a thermal condition for every
tag of the mesh, and `setup` refuses a condition keyed by a tag the mesh does
not have.

```julia
bcs = BoundaryConditionSet(
    thermal = (
        top    = RadiativeFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0) +
                 ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (
        top    = AmbientGasBC(GasReservoir(P=0., T=300.)),     # gas escapes at the surface
        bottom = ImpermeableBC(),
    ),
)
```

### 16.4.6 Initial conditions

The initial conditions belong to `PyrolysisProblem` (§16.5.1): `T_initial`, a
function `z -> T` (default `z -> 300.0`), and `ξ_initial`, a vector of
functions `z -> ξ_j`, one per component, or a single function used for every
component. Without `ξ_initial` the first component, if it is a solid, starts at
its density at 300 K and every other component at zero. `setup` evaluates the
functions at the cell centres of the reference mesh, `Tᵢ = T_initial(zᵢ)` and
`m_ij = ξ_j(zᵢ) Vᵢ` (`Pyrolysis.Problem.initial_fields(problem, mesh)` returns
both fields). The profile helpers of `Pyrolysis.Problem` build such functions:

| Helper | Signature |
|---|---|
| `step_profile` | `step_profile(z; z_interface, value_left, value_right)` |
| `linear_profile` | `linear_profile(z; z_min, z_max, value_min, value_max)` |
| `tanh_profile` | `tanh_profile(z; z_center, width, value_left, value_right)` |
| `gaussian_profile` | `gaussian_profile(z; z_center, width, amplitude, baseline = 0.0)` |

`background_gas_concentration(material, species, ξ, T; P = 101325.0)` returns
the concentration [kg·m⁻³] of the gas component `species` (an index or a name)
whose partial pressure is `P` in the pore space of the condensed composition
`ξ` (an `NTuple` with one entry per component) at temperature `T`. Use the
intended initial pore gas for either transport model (§7.5.1).

Incoming gas carries its specified reservoir enthalpy. Outgoing gas carries
surface enthalpy. The heat closure determines whether a coupled LTE-film
adjustment is required (§5.4; Technical Reference §12).

### 16.4.7 Blowing models and gas data

A `ConvectiveBC`'s `blowing` keyword takes one of these models
(`AbstractBlowingModel` is their supertype); each evaluates its declared flux domain and returns
the blowing ratio `R_b = h_eff/h` (Chapter 5, §5.2.3):

| Model | Constructor | Ratio |
|---|---|---|
| `NoBlowing` | `NoBlowing()` (the default) | `R_b = 1` |
| `CouetteFilm` | `CouetteFilm(λ = 0.5)` or `CouetteFilm(; λ)` | `B_λ/(e^{B_λ} − 1)`, `B_λ = 2λ Σ_j N_j c_p,j(T_s)/h`, including signed counterflow |
| `PropertyRatioFilm` | `PropertyRatioFilm(free_stream; n_Γ = 0.6, volatiles = nothing)` | `B_Γ/(e^{B_Γ} − 1)` with `B_Γ = (ṁ″ c_p,∞/h)(Γ_∞/Γ_g)^n_Γ`, `Γ = k/c_p` |
| `PoolingFilm` | `PoolingFilm(free_stream; n_Γ = 0.6, n_D = 0.5, ψ0 = 0.59, Λc = 0.35, Λw = 0.10, volatiles = nothing, diffusion_volume = (;))` | `B_Γ/(e^{B_Γ/ψ} − 1)` with a pooling factor `ψ` falling from 1 to `ψ_0` (keyword `ψ0`) when the gas at the face is denser than the free stream |
| `FunctionBlowing` | `FunctionBlowing(f)` | `R_b = f(face, h, T_inf)`; `f` returns 1 at zero species fluxes, is smooth and accepts dual numbers |

`CouetteFilm` supports signed fluxes. `PropertyRatioFilm` and `PoolingFilm`
require their declared outward-flow domain and reject unsupported counterflow.
The gas data the models read:

| Name | Signature and meaning |
|---|---|
| `GasSpecies` | `GasSpecies(name; M, c, k, μ = nothing, V_D = NaN, D = nothing, T_range = (200.0, 3000.0))`: molar mass [kg/mol], heat capacity and conductivity (numbers or property functions), viscosity (`nothing` is the Eucken estimate), Fuller diffusion volume, binary diffusivity into the free stream (overrides Fuller), and the temperature range the laws are evaluated on |
| `NITROGEN` | nitrogen as a free stream: NASA-7 heat capacity, polynomial conductivity and viscosity, `V_D = 18.5`, range 200–1250 K |
| `nasa7_heat_capacity` | `nasa7_heat_capacity(a, M)` or `nasa7_heat_capacity(a_low, a_high, M; T_common = 1000.0)`: `c_p` [J/(kg·K)] from five NASA-7 coefficients per range |
| `eucken_viscosity` | `eucken_viscosity(k, cp, M)`: `μ = k/(c_p + 5R/(4M))` [Pa·s] |
| `fuller_diffusivity` | `fuller_diffusivity(T, M_a, V_a, M_b, V_b; P = P_REF)`: the Fuller binary diffusivity [m²/s], molar masses in kg/mol, `P` in Pa |
| `BlowingFace` | `BlowingFace(ṁ, T_s, material)`: what a model sees, the outward flux of every component `ṁ` (an `NTuple`, positive leaving), the surface temperature and the material |
| `blowing_state` | `blowing_state(model, face, h, T_inf) -> (; χ, φ, ψ, Y, Λ, states)`: the field `χ` is the ratio `R_b`, `φ` the film parameter `B_λ` or `B_Γ`, `ψ` the pooling factor, `Y` the face mass fraction of the volatiles, `Λ` the face buoyancy, `states` the number of self-consistent states (1 or 2) |
| `blowing_ratio` | `blowing_ratio(model, face, h, T_inf) -> R_b` |
| `blowing_factor` | `blowing_factor(x) = x/(e^x − 1)`, the film factor `F_b`, smooth for every real `x`, one at `x = 0` |

The solution records the applied ratio per boundary tag
(`boundary_blowing_ratio`) and, for the `:top` faces, the ratio, the pooling
factor, the face fraction and the face buoyancy in `extras` (§16.5.4).

---

## 16.5 The problem and solution objects

See Chapters 7 and 9 and Technical Reference §3 and §15. `PyrolysisProblem`,
`PyrolysisSolution`, `SolverStatistics` and the accessors are top-level
exports.

### 16.5.1 `PyrolysisProblem`

```julia
PyrolysisProblem(; mesh, material, bc_set::BoundaryConditionSet,
                   T_initial = z -> 300.0,
                   ξ_initial = nothing,
                   tspan::Tuple{Real,Real},
                   dt_initial::Real = 1e-4,
                   output_times::Vector{<:Real} = Float64[])

PyrolysisProblem(mesh, material, bc_set, tspan; kwargs...)   # positional form
```

| Keyword | Meaning | Default |
|---|---|---|
| `mesh` | a `Mesh{1}` from a constructor of §16.3.1 | required |
| `material` | a `Material` | required |
| `bc_set` | a `BoundaryConditionSet` | required |
| `T_initial` | initial temperature `z -> T` [K] | `z -> 300.0` |
| `ξ_initial` | initial concentrations: a vector of functions `z -> ξ_j` [kg·m⁻³], one per component, or one function for every component | `nothing`: the first component, if a solid, at its density at 300 K, the others zero |
| `tspan` | `(t_start, t_end)` [s] | required |
| `dt_initial` | the first step `solve` takes by default (its `dt`), and the step behind the automatic output count | `1e-4` s |
| `output_times` | the saved times `solve` uses by default | `Float64[]`: `min(101, max(11, round((t_end − t_start)/dt_initial/100)))` evenly spaced times over `tspan`, both ends included |

The fields are `mesh`, `material`, `bc_set`, `T_initial`, `ξ_initial` (always a
vector of functions), `tspan` (a `Float64` tuple), `dt_initial` and
`output_times`.

```julia
problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bcs,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],   # virgin at its density, no gas
    tspan = (0.0, 100.0),
)
```

`validate_problem(problem; strict = false)` (`Pyrolysis.Problem`) checks the
mesh positions, the material's components and reaction indices, the boundary
tags and conditions, the time span, the number of initial-concentration
functions and, under gas transport, the initial porosity; it collects
every error and throws them in one message (an `ErrorException`), and it warns
about invalid material or boundary data.
`strict = true` adds the stoichiometry, an evaluation of the initial-condition
functions and the cell volumes. `is_valid_problem(problem)` returns a `Bool`.
`solve` runs `validate_problem` unless `validate = false` (Chapter 7, §7.9).

### 16.5.2 `PyrolysisSolution`

`solve` and `integrate` return a `PyrolysisSolution`. `n_t` is the number of
saved times, `NC` the number of components; every time-indexed array is
filled at every saved time.

| Field | Type | Meaning |
|---|---|---|
| `problem` | `PyrolysisProblem` | the problem solved; its mesh is unchanged |
| `t` | `Vector{Float64}` | saved times [s] |
| `T` | `Matrix{Float64}` | cell temperatures, `n_cells × n_t` [K]; `NaN` in a deactivated cell |
| `ξ` | `Array{Float64,3}` | concentrations, `NC × n_cells × n_t` [kg·m⁻³]; `NaN` in a deactivated cell |
| `m` | `Array{Float64,3}` | component masses per cell, `NC × n_cells × n_t` [kg]; zero in a deactivated cell |
| `z` | `Matrix{Float64}` | node positions, `n_nodes × n_t` [m] on a moving mesh, `0 × 0` on a fixed one; `NaN` beyond the outermost live node |
| `correction_velocity` | `Matrix{Float64}` | node corrections `c = w − v_s` [m/s], `n_nodes × n_t` on a moving mesh (zero without a correction), `0 × 0` on a fixed one; `NaN` beyond the outermost live node |
| `cell_centers` | `Matrix{Float64}` | cell-centre positions, `n_cells × n_t` [m]; `NaN` in a deactivated cell |
| `boundary_T` | `Dict{Symbol,Vector{Float64}}` | per boundary tag, the area-weighted face temperature [K]; after deactivation, of the face bounding the live cells on that side |
| `surface_T` | `Vector{Float64}` | temperature of the `:top` faces [K] |
| `mass_loss_rate` | `Vector{Float64}` | mass leaving through every boundary face, per unit current area of the `:top` faces [kg·m⁻²·s⁻¹] |
| `boundary_mass_flux` | `Dict{Symbol,Vector{Float64}}` | per tag, the outward mass flux of every component summed, area-weighted [kg·m⁻²·s⁻¹] |
| `boundary_blowing_ratio` | `Dict{Symbol,Vector{Float64}}` | per tag, the area-weighted ratio `h/h₀` its convective condition applied (1 without outflow, without a blowing model, or on a face fixed in closed form) |
| `extras` | `SolutionExtras` | the secondary series (§16.5.4) |
| `events` | `Vector{DeactivationEvent}` | the deactivation events, in order (§16.5.4) |
| `outcome` | `Symbol` | `:running` unless the deactivation or the thickness stop ended the solve: `:burned_out`, `:interior` or `:min_thickness` |
| `retcode` | `Symbol` | the integrator's return code (§16.7.1) |
| `solve_time` | `Float64` | wall time of the solve [s] |
| `statistics` | `SolverStatistics` | the solver's counters (§16.5.4) |
| `raw_solution` | ODE solution | the integrator's solution: the raw states and the interpolant |
| `definition` | `CoreDefinition` | the definition of the solve, through which `sol(t)` reconstructs the fields |

`T`, `ξ` and `m` are cell-centred, cell 1 at the back (`z = 0`); `ξ[k, i, n]`
is component `k` in cell `i` at saved time `n`. An event that leaves the state
inadmissible or the node velocities inconsistent ends the solve with a warning,
`retcode` `:Failure` and outcome `:running`.

```julia
sol = solve(problem)
times    = sol.t
Tsurf    = sol.surface_T                       # K, per saved time
mlr      = sol.mass_loss_rate                  # kg/(m² s), per current top area
mlrgauge = sol.extras.MLR_gauge                # kg/(m² s), per initial top area
Tprofile = sol.T[:, end]                       # final temperature profile
ξvirgin  = sol.ξ[1, :, end]                    # final virgin concentration profile
```

### 16.5.3 Accessors and reconstruction

All top-level exports:

| Function | Returns |
|---|---|
| `get_times(sol)` | `sol.t` |
| `get_temperatures(sol)` | `sol.T` |
| `get_concentrations(sol)` | `sol.ξ` |
| `get_gas_generation_rate(sol)` | `sol.extras.gas_generation_rate`, the volumetric production of the live cells (not the mass-loss rate) |
| `n_timesteps(sol)` | `length(sol.t)` |
| `get_profile(sol, field, i)` | one field at saved time `sol.t[i]`: `:T`; `:ξ` (`NC × n_cells`) or `:ξj` / `:xij` for component `j`; `:m` (`NC × n_cells`) or `:mj`; `:z`, the node positions (on a fixed mesh the mesh's). Deactivated cells, and nodes beyond the outermost live one, are `NaN`; another symbol is an `ArgumentError` |
| `get_cell_positions(sol; at = lastindex(sol.t))` | the cell centres at saved time `sol.t[at]` [m], `NaN` for a deactivated cell; the initial centres are `at = 1` when the first saved time is the start |
| `sol(t)`, `interpolate_solution(sol, t)` | the NamedTuple `(; t, T, ξ, m, z, cell_centers)` at any `t` in the span, reconstructed from the integrator's interpolant through the state mapping, with the mask of the last event at or before `t` (`z` empty on a fixed mesh) |
| `succeeded(sol)` | `true` for `retcode` `:Success`, and for `:Terminated` with outcome `:burned_out` or `:min_thickness`; `false` for `:interior`, `:Failure`, and a termination by a callback of yours |
| `completed(sol)` | `true` for `:Success` only, the end of the time span reached |
| `conservation_audit(sol[, recorder])` | the mass and energy closures (§16.12) |

`sol(t)` refuses a `t` strictly between two saved states that a deactivation
event separates (the interpolant there is not the trajectory), with an
`ArgumentError` naming the event, and a `t` outside the span. At a saved time it
returns the saved fields exactly. `show(sol)` prints the return code, the
outcome, the number of saved times and the solve time; the `text/plain` form
adds the time range, the live cells, the invocations and the accepted steps.

### 16.5.4 The `extras` series, the events and the statistics

`sol.extras` is a `SolutionExtras` (`Pyrolysis.Outputs`), a typed NamedTuple
whose keys are always present, one value per saved time unless noted.
Quantities of the `:top` or `:bottom` faces are `NaN` without that tag or once
no live cell borders it; per-area quantities divide by the current area of the
faces they name unless noted.

| Key | Units | Meaning |
|---|---|---|
| `gas_generation_rate` | kg·m⁻²·s⁻¹ | gas the reactions of the live cells produce, per unit initial `:top` area |
| `back_temperature` | K | temperature of the `:bottom` faces |
| `total_mass` | kg·m⁻² | summed component masses, per unit initial `:top` area |
| `thickness_history` | m | distance from the first node to the outermost live node |
| `min_cell_thickness` | m | thickness of the thinnest live cell |
| `surface_mass_flux` | kg·m⁻²·s⁻¹ | outward mass flux of the `:top` faces |
| `surface_heat_flux_net` | W·m⁻² | heat flux the `:top` condition applies, positive into the domain (under a prescribed temperature, the conducted flux); excludes the incident flux absorbed in depth under `BEER_LAMBERT` |
| `surface_heat_flux_absorbed` | W·m⁻² | incident radiation absorbed at the face |
| `surface_heat_flux_reradiation`, `surface_heat_flux_convection` | W·m⁻² | re-radiative and convective losses, positive leaving |
| `surface_heat_flux_conduction` | W·m⁻² | flux conducted from the `:top` faces into their cells, `k (T_s − T_cell)/d` |
| `surface_blowing_ratio`, `surface_pooling_factor`, `surface_face_fraction`, `surface_face_buoyancy` | – | blowing diagnostics of the `:top` balance (1, 1, 0 and 1 without outflow or without a blowing model) |
| `volumetric_absorption`, `radiation_transmitted` | W·m⁻² | intensity absorbed in depth and leaving through the last face of the beam, per unit exposed area (zero without `BEER_LAMBERT`) |
| `surface_emissivity`, `surface_conductivity` | –, W·m⁻¹·K⁻¹ | emissivity and conductivity of the cell under the `:top` faces |
| `back_heat_flux` | W·m⁻² | heat flux the `:bottom` condition applies, positive into the domain |
| `MLR_gauge` | kg·m⁻²·s⁻¹ | mass-loss rate per unit initial `:top` area, `mass_loss_rate · A(t)/A₀` |
| `cross_section_area_history` | m² | the cross-section measure `A(t)` |
| `live_cells` | – | number of live cells |
| `deactivation_times` | s | time of every deactivation, one entry per deactivated cell, in order (not indexed by saved time) |
| `live_mass`, `live_energy` | kg, J | mass and energy of the live cells |
| `mass_outflow`, `energy_rate` | kg·s⁻¹, W | mass outflow and enthalpy rate of the live region |

`sol.events` lists one `DeactivationEvent` (`Pyrolysis.Deactivation`) per
deactivated cell, with the fields `t` [s], `invocation` (numbered from one,
shared by the cells one invocation deactivates), `cell`, `receiver` (zero when
no live cell was left), `fraction` (the cell's condensed mass as a fraction of
its initial one when it went), `rows` (the state rows masked),
`created_volume` [m³], `lost_mass` (per component, kg), `lost_enthalpy` [J] and
`regularization` [J]. A merge of several cells records one event per dead cell,
the outermost first, and writes the merge's created volume, losses and
regularization to that first event (Chapter 9, §9.3; Technical Reference §14).

`sol.statistics` is a `SolverStatistics` with the fields `accepted`,
`rejected`, `rhs_evaluations`, `jacobians`, `linear_solves`,
`nonlinear_failures`, `domain_rejections` (proposals the domain test
rejected), `largest_constraint_residual` (the largest scaled algebraic-row
residual of the accepted states, one or below passing; zero for a fixed-step
solve) and `invocations` (deactivation invocations that deactivated cells).

---

## 16.6 What `setup` decides

See Chapter 7 and Technical Reference §15. A problem has no mode argument:
`setup` (which `solve` calls) reads every modelling choice once, from the
keywords and from the problem, and builds the state layout and the kernels to
match.

| Choice | Selected by | Values | Default |
|---|---|---|---|
| Mesh motion | `use_ale` | fixed, moving (`w = v_s + c`) | fixed |
| Radiation | `radiation_model` | `NO_RADIATION`, `SURFACE_ABSORPTION`, `BEER_LAMBERT` | `NO_RADIATION` |
| Gas transport | material model and complete mass conditions | `LumpedTransfer` or `DustyGas`; EOS for every gas-bearing material, including sealed gas | from the problem |
| Cross section | the material's `lateral_shrinkage_law` | fixed, variable (needs `use_ale = true` and planar coordinates) | fixed |
| Energy state | `energy_state` | `:temperature`, `:enthalpy` | `:temperature` |
| Deactivation | `deactivation` | off, on | off |
| Correction | `correction`, or by rule | `NoCorrection()`; `EquidistributionCorrection()` where the mesh moves, deactivation is on and some cell can vanish (its components, converted completely, leave less condensed mass than `threshold` of its initial condensed mass) | by rule |
| Layer interfaces | `lagrangian`, or by rule | node ids moving with the solid; under the equidistribution correction the interior nodes between cells that start with different condensed components | by rule |

The coordinate system and the cross-section measure are the mesh's.

Choosing for common scenarios:

- TGA-like or fixed slab: the defaults (a fixed mesh, `NO_RADIATION`).
- Cone calorimeter with recession: `use_ale = true` with
  `radiation_model = SURFACE_ABSORPTION` (opaque) or `BEER_LAMBERT`
  (semi-transparent).
- Burnout of a layer: add `deactivation = true`; on a moving mesh `setup`
  installs the equidistribution correction (Chapter 10).
- Pressure-driven gas flow: give a condensed component `κ > 0` and add pressure
  conditions (§16.4.3).
- Lateral shrinkage or intumescence: give the material a
  `lateral_shrinkage_law` and solve with `use_ale = true`.

A combination that is not implemented is refused by `setup` with an
`ArgumentError` before anything is integrated: a lateral shrinkage law without
`use_ale` or on a non-planar mesh; `BEER_LAMBERT` on a non-planar mesh or with
a `RadiativeFluxBC` on a tag other than `:top`; a correction without
`use_ale`; a `w_override` without `use_ale`, with a reacting material, with
deactivation or with a correction; an explicit `lagrangian` list that omits a
layer interface under the equidistribution correction; with deactivation, a
`band` outside `[threshold, 1)`; an `energy_state` other than `:temperature` or `:enthalpy`
(Chapter 7, §7.9.1; Chapter 15, §15.2).

---

## 16.7 The solve keywords

| Keyword | Type | Default | What it does |
|---|---|---|---|
| `energy_state` | `Symbol` | `:temperature` | The energy unknown of each cell: `:temperature` (the cell temperature) or `:enthalpy` (the sensible-enthalpy inventory) |
| `abstol` | `Real` | `1e-8` | Absolute tolerance of the temperature rows [K]; the enthalpy rows take `threshold · m_s,i(0) · c_p,ref · abstol`, the other rows derive theirs from `reltol` and `threshold` |
| `pressure_tolerance` | `1.0` | Absolute local pressure-error budget [Pa], with proposed-state pore/thermal/geometry sensitivities |
| `reltol` | `Real` | `1e-6` | Relative tolerance of every row; also scales the mass, position, velocity, intensity and cross-section tolerances, the domain test's constraint allowance, and the widening of a located event's threshold |
| `threshold` | `Real` | `1e-2` | `ε_m`, the deactivation threshold as a fraction of a cell's initial condensed mass; also scales the mass and position tolerances, with or without `deactivation` |
| `maxiters` | `Int` | `1_000_000` | Maximum number of integrator steps |
| `dt` | `Real` | `problem.dt_initial` | First step [s]; the fixed step when `adaptive = false`; the time scale of the node-velocity tolerances |
| `adaptive` | `Bool` | `true` | Adaptive step control; `false` takes fixed steps of `dt` without the domain test and ends with `:Failure` at the first inadmissible accepted state |
| `dtmin` | `Real` | `1e-12` | Minimum step [s] |
| `dtmax` | `Real` | `Inf` | Maximum step [s] |
| `saveat` | number or vector | `problem.output_times` | Saved times, a vector or a spacing [s]; the final state is always saved; `Float64[]` saves every accepted step |
| `callback` | callback, recorder, tuple or `nothing` | `nothing` | An OrdinaryDiffEq callback or `CallbackSet` (each discrete callback a read-only observer), a `ConservationRecorder`, or a tuple of recorders and callbacks |
| `progress` | `Bool` | `false` | Draw a progress bar, a read-only observer |
| `validate` | `Bool` | `true` | Run `validate_problem(problem; strict = false)` before setting up |
| `use_ale` | `Bool` | `false` | Move the mesh: the node positions and velocities join the state, `w = v_s + c` |
| `w_override` | function or `nothing` | `nothing` | Prescribed node velocity `w(z, t)` [m/s] over a solid at rest (verification); needs `use_ale = true`, a non-reacting material, no deactivation and no correction |
| `min_thickness` | `Real` | `1e-6` | On a moving mesh, end with `:Terminated` and outcome `:min_thickness` when the live thickness falls below it [m]; must be positive |
| `radiation_model` | `RadiationModel` | `NO_RADIATION` | `NO_RADIATION`, `SURFACE_ABSORPTION` or `BEER_LAMBERT` (§16.2.6) |
| `linear_solver` | factorization or `nothing` | `nothing`, which is `StructuredKLU()` | Linear solver of the default policy; ignored when `integrator_policy` is given |
| `integrator_policy` | `IntegratorPolicy` or `nothing` | `nothing`, which is `default_policy(linear_solver)` | The integrator with its restart and observer-flag rules (§16.8) |
| `sources` | sources object | `NoSources()` | External volumetric sources for verification: `Pyrolysis.Assembly.VolumetricSources(mass, energy)` |
| `deactivation` | `Bool` | `false` | Deactivate cells whose condensed mass falls below `threshold` of their initial one (§16.11) |
| `band` | `Real` or `nothing` | `nothing`, which is `max(threshold, 0.1)` when the equidistribution correction moves a layer and `threshold` otherwise | `ε_band` in `[threshold, 1)`: on a moving mesh, the run of cells from the exposed side below this fraction goes with a cell that crosses `ε_m` |
| `correction` | `Correction` or `nothing` | `nothing`, which is the rule of §16.6 | `NoCorrection()` or `EquidistributionCorrection(...)` (§16.11) |
| `lagrangian` | vector of node ids or `nothing` | `nothing`, which is the layer interfaces | Nodes that always move with the solid, beside the mesh's own |
| `located_events` | `Bool` | `true` | Locate each deactivation inside the step at the first threshold crossing; the end-of-step test stays as a safety net |

`Pyrolysis.Integration._SOLVE_KEYWORDS` is the authoritative list, and `solve`
refuses any other keyword with `ArgumentError: solve(): unknown keyword` before
doing any work. Chapter 8 explains each keyword in context.

### 16.7.1 Behaviour, tolerances and outcomes

```julia
solve(problem::PyrolysisProblem; kwargs...) -> PyrolysisSolution
```

`solve` runs `validate_problem` (unless `validate = false`), then
`setup(problem; ...)` with the keywords of §16.7.2 (its `time_scale` set to
`dt`), then `integrate(cp; ...)` with those of §16.7.3, and reconstructs the
solution at the saved times. It never modifies `problem` or its mesh, so a
problem can be solved again or shared by the points of a sweep.

The absolute tolerances are set per row by `setup`, from `abstol`, `reltol`,
`threshold` and each cell's initial condensed mass `m_s,i(0)` (its total
initial mass when it holds no condensed mass, or 1 kg/m³ times its volume when
it starts empty):

| Row | Absolute tolerance |
|---|---|
| temperature | `abstol` [K] |
| enthalpy inventory (`energy_state = :enthalpy`) | `threshold · m_s,i(0) · c_p,ref · abstol` [J] |
| condensed component mass of cell `i` | `reltol · threshold · m_s,i(0)` [kg] |
| gas component mass | initial `pressure_tolerance · M_j · φ_i · V_i/(N_g R T_i)` [kg], supplemented by the proposed-state pressure norm |
| node position (moving mesh) | `reltol · threshold · h_n`, `h_n` the smallest initial size of the cells around node `n` |
| node velocity (moving mesh) | the position tolerance over `dt` |
| face intensity (`BEER_LAMBERT`) | `reltol` times the largest incident intensity over the time span (at least 1 W/m²) |
| conversion-history pair (lateral shrinkage) | the cell's mass tolerance |
| cross-section area and lateral rate (lateral shrinkage) | `reltol · A₀`; the smallest node-velocity tolerance over the initial thickness |

The pressure norm evaluates the proposed temperature, masses, pore volume
and geometry, including caloric inversion in the enthalpy formulation. Gas
trial-guard scales follow those same evolving quantities. Tighten
`pressure_tolerance` and compare complete trajectories when qualifying a
solution; it is a local-error budget, not a global pressure bound.

The domain test rejects invalid masses, temperatures or geometry and checks
algebraic constraints. Every gas-bearing state, including sealed gas, must
have physical pore volume above twice the divisor floor plus its mass,
thermal and geometric uncertainty margin. Accepted and post-event states
must also remain in the declared transport domain. Repeated constraint
rejections or an invalid fixed step terminate with a diagnosed failure.

How a solve ends:

| `retcode` | `outcome` | Meaning |
|---|---|---|
| `:Success` | `:running` | reached `tspan[2]` (cells may have been deactivated on the way) |
| `:Terminated` | `:burned_out` | the last live cell was deactivated |
| `:Terminated` | `:min_thickness` | a moving mesh thinned below `min_thickness` |
| `:Terminated` | `:interior` | a cell below the threshold lies inward of a live cell; `solve` warns, naming the cell and the time |
| `:Terminated` | `:running` | a callback of yours called `terminate!` |
| `:Failure` | `:running` | the domain test gave up, a fixed step left the physical domain, or a deactivation event left the state inadmissible or the node velocities inconsistent; `solve` warns with the reason |
| `:MaxIters` | `:running` | `maxiters` steps before the end |
| `:DtLessThanMin` | `:running` | the step fell below `dtmin` |

A solve that stops early ends on the last state it accepted, which is saved.
`succeeded` and `completed` combine the two fields (§16.5.3). Chapter 15
covers the failures.

```julia
# Opaque cone sample, receding surface, surface cells deactivated as they burn out:
recorder = ConservationRecorder()
sol = solve(problem;
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            deactivation = true,
            min_thickness = 0.5e-3,
            abstol = 1e-9, reltol = 1e-7,
            callback = recorder)
audit = conservation_audit(sol, recorder)
```

### 16.7.2 `setup` keywords

```julia
setup(problem; kwargs...) -> CoreProblem
```

| `setup` keyword | Default | Meaning |
|---|---|---|
| `energy_state`, `reltol`, `abstol`, `threshold`, `pressure_tolerance`, `use_ale`, `radiation_model`, `w_override`, `correction`, `lagrangian`, `deactivation`, `band`, `sources` | as in §16.7 | as in §16.7 |
| `time_scale` | `problem.dt_initial` | time scale of the node-velocity tolerances [s]; `solve` passes its `dt` |
| `closure` | `nothing` | the mesh closure of a mesh of two or three dimensions, `:potential` (its default) or `:balance`; a line moves by its balance rows and takes no other |
| `slip` | `()` | `tag => k` pairs fixing velocity component `k` to zero on a tag (two or three dimensions) |
| `tensor` | `nothing` | each cell's dimensionless conductivity tensor (two or three dimensions; the identity by default) |
| `history` | `false` | carry the conversion-history pair (two or three dimensions; on a line a lateral shrinkage law adds it) |
| `deformation` | `false` | carry the deformation gradient (two or three dimensions) |

The last five apply to the meshes of two and three dimensions only, a residual
level that is not a supported user path for full solves (§16.3.2); on a line
`slip`, `tensor`, `history` and `deformation` are refused. `setup` returns a `CoreProblem` with the fields `problem`, `def`
(the definition the residual reads: `layout`, `setup`, `kernels`, `closure`,
`sources`), `ws` (its buffers), `u0` (the consistent initial state), `reltol`,
`abstol` (the per-row vector), `jacobian` (the colored Jacobian cache), `T0`
and `m0` (the initial temperatures and masses) and `deactivation` (the
deactivation state, `nothing` without deactivation). Building the colored
Jacobian is the expensive part of `setup`.

### 16.7.3 `integrate` keywords

```julia
integrate(cp::CoreProblem; kwargs...) -> PyrolysisSolution
```

| `integrate` keyword | Default | Meaning |
|---|---|---|
| `u0` | `cp.u0` | initial state; its algebraic rows must be consistent (`Pyrolysis.Integration.consistent_state!(u, cp.def, cp.ws, t)`), and an inadmissible one is refused |
| `maxiters`, `dt`, `adaptive`, `dtmin`, `dtmax`, `saveat`, `callback`, `min_thickness`, `linear_solver`, `integrator_policy`, `located_events`, `progress` | as in §16.7, with `dt = cp.problem.dt_initial` and `saveat = cp.problem.output_times` | as in §16.7 |
| `max_constraint_rejections` | `50` | consecutive constraint rejections of the domain test before the solve ends with `:Failure` |
| `started` | `time()` | the wall-clock time the reported `solve_time` counts from |

Every call of `integrate` starts with every cell live, the layer times unset and
the caloric inversion from the set-up temperatures, so repeated integrations
of one `CoreProblem` are bit-identical (Chapter 11, §11.8).

```julia
cp   = setup(problem; use_ale = true, radiation_model = SURFACE_ABSORPTION)
sol1 = integrate(cp)
sol2 = integrate(cp; integrator_policy = fbdf_policy(StructuredKLU()))
```

---

## 16.8 Integrator, linear solver and Jacobian

See Chapters 8 and 11 and Technical Reference §15. The policy constructors and
`StructuredKLU` are top-level exports.

### 16.8.1 Integrator policies

An integrator policy is an OrdinaryDiffEq method together with two rules of
its contract: the **restart** rule, the step proposed after an event that
modifies the state (`:solver_chosen` keeps the integrator's proposal,
`:forced_dt` proposes `DT_AFTER_EVENT = 1e-6` s, `:heuristic` proposes the
integrator's initial-step estimate at the state the event left), and the
**observer-flag** rule, what a read-only observer does with the integrator's
derivative-discontinuity flag (`:keep` leaves it raised, `:clear` clears it).
Each constructor takes the linear solver as its argument:

| Constructor | Method | Restart | Observer flag |
|---|---|---|---|
| `default_policy(ls)` | `rodas5p_policy(ls)` | `:heuristic` | `:clear` |
| `rodas5p_policy(ls; restart = :heuristic)` | `Rodas5P`, a Rosenbrock method: the colored Jacobian at every step, a linear solve per stage, no Newton iteration | any of the three | `:clear` |
| `fbdf_policy(ls)` | `FBDF`, variable-order backward differentiation with Newton iteration | `:forced_dt` | `:clear` |
| `kencarp4_policy(ls)` | `KenCarp4`, an ESDIRK method with Newton iteration, the raw error estimate | `:solver_chosen` | `:keep` |

Every policy differentiates the right-hand side in time by finite differences
(`AutoFiniteDiff()`) where its method needs that derivative; the Jacobian with
respect to the state is the colored one (§16.8.3). For the Newton-based
policies, after an accepted step whose nonlinear solves needed retries, the
next `FAILURE_MEMORY_STEPS = 20` accepted steps are held at or below
`FAILURE_MEMORY_GROWTH = 1.25` times that step; `Rodas5P` has no Newton
failures. A policy of your own is
`Pyrolysis.Integration.IntegratorPolicy(alg, restart, observer_flag, restart_dt)`
with an OrdinaryDiffEq method that accepts a mass matrix and a user Jacobian;
`restart` must be one of the three symbols and `observer_flag` `:keep` or
`:clear`. Only the three constructors are tested.

```julia
solve(problem)                                                     # default_policy(StructuredKLU())
solve(problem; integrator_policy = fbdf_policy(StructuredKLU()))
solve(problem; integrator_policy = kencarp4_policy(StructuredKLU()))
solve(problem; integrator_policy = rodas5p_policy(StructuredKLU(); restart = :forced_dt))
```

### 16.8.2 The linear solver

`StructuredKLU(inner = KLUFactorization())` is the default linear solver
(`Pyrolysis.Integration.default_linear_solver()`): KLU on the entries of the
integrator's matrix `W = J − M/(γΔt)` that have been nonzero at some step of the
solve, the structure analysed afresh only when that union of entries grows.
Any LinearSolve factorization of a `SparseMatrixCSC` is accepted;
`UMFPACKFactorization()` and plain `KLUFactorization()` are the ones the test
suite runs, and repeated solves under each of the three are bit-identical.

```julia
using LinearSolve: UMFPACKFactorization

solve(problem; linear_solver = UMFPACKFactorization())                      # default policy
solve(problem; integrator_policy = rodas5p_policy(UMFPACKFactorization()))  # inside a policy
```

`linear_solver` configures the default policy only; with an explicit
`integrator_policy` it is ignored, and the policy's own argument is the linear
solver (Chapter 11, §11.4).

### 16.8.3 The Jacobian

The Jacobian is colored sparse forward-mode automatic differentiation on a
sparsity pattern built from the mesh connectivity: a greedy column coloring
groups columns with no row in common, and one dual-number evaluation of the
residual per batch of up to 24 colors fills the matrix. It is exact to
round-off in every configuration, and there is no keyword or strategy to
choose. A dense forward-mode Jacobian is the reference the test suite holds it
against. The building blocks are in `Pyrolysis.ColoredJacobian` (and
`Internal`):

| Function | Meaning |
|---|---|
| `jac_prototype(cj)` | the pattern's sparse matrix with zero values, `cj = cp.jacobian` |
| `jacobian!(J, cj, u, t)` | the colored Jacobian at `(u, t)`, written into `J` |
| `jacobian_pattern(def)` | the Boolean sparsity pattern of a definition, `def = cp.def` |
| `color_count(def)` | the number of colors |
| `dense_jacobian(def, ws, u, t)` | the dense forward-mode reference |
| `missing_entries(J_ref, pattern)` | the structurally nonzero entries of `J_ref` the pattern lacks |
| `compare_entries(J, J_ref; rtol = 1e-12, atol = 0.0)` | `(; passed, max_abs_dev, worst_ratio, worst_index)`, each entry tested against a floor of its row's largest magnitude |

```julia
using Pyrolysis.Internal: jac_prototype, jacobian!, dense_jacobian, compare_entries

cp = setup(problem)
t0 = cp.problem.tspan[1]
J = jac_prototype(cp.jacobian)
jacobian!(J, cp.jacobian, cp.u0, t0)
compare_entries(J, dense_jacobian(cp.def, cp.ws, cp.u0, t0)).passed   # true
```

Chapter 11 (§11.2, §11.7) describes the pattern and the check on a state after
deactivation.

---

## 16.9 Parameter sweeps

See Chapter 12. `parameter_sweep` is a top-level export.

```julia
parameter_sweep(problem, grid; modify, extract = identity,
                parallel = nworkers() > 1, solve_kwargs...) -> Vector
```

| Argument | Meaning |
|---|---|
| `problem` | the base problem, never modified |
| `grid` | any iterable of parameter points (a vector of NamedTuples, `Iterators.product(...)`, …), collected before the sweep |
| `modify` | `(problem, point) -> problem`, returning a fresh problem (it may share the base problem's mesh, which a solve never changes); required, an `ArgumentError` without it |
| `extract` | `solution -> result`; with workers, prefer a small NamedTuple to a whole solution |
| `parallel` | `Distributed.pmap` over the workers, by default whenever there is more than one; otherwise a serial `map` |
| `solve_kwargs...` | passed to every `solve` |

The results come back in the grid's order. The workers need the package and
the definitions `modify` and `extract` use (`@everywhere`). Since there is no
copy-with-overrides constructor, `modify` rebuilds what it changes:

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

@everywhere function modify(problem, p)
    m = problem.material
    reaction = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 2e5)
    material = Material(; name = m.name, components = m.components,
                        reactions = (reaction,), mixing = m.mixing,
                        lateral_shrinkage_law = m.lateral_shrinkage_law,
                        depletion_limiter = m.depletion_limiter,
                        gas_viscosity = m.gas_viscosity, gas_transport = m.gas_transport)
    return PyrolysisProblem(; mesh = problem.mesh, material, bc_set = problem.bc_set,
                            T_initial = problem.T_initial, ξ_initial = problem.ξ_initial,
                            tspan = problem.tspan, dt_initial = problem.dt_initial,
                            output_times = problem.output_times)
end
@everywhere extract(sol) = (peak_T = maximum(sol.surface_T),
                            final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```

---

## 16.10 Sensitivity analysis

See Chapter 13 and Technical Reference §17. There is no sensitivity function:
sensitivities are finite differences of whole solves, on problems that a
closure `p_inject(base, θ)` rebuilds from a parameter vector, the same pattern
as `modify` in §16.9. No derivative is taken through `solve`.

```julia
θ0 = [log(8.5e12), 188e3]                  # ln A, E
function p_inject(prob, θ)
    m = prob.material
    rxn = Reaction(:decomp, 1 => 2; A = exp(θ[1]), E = θ[2], h = 870e3)
    mat = Material(name = m.name, components = m.components,
                   reactions = (rxn,), mixing = m.mixing)
    return PyrolysisProblem(mesh = prob.mesh, material = mat, bc_set = prob.bc_set,
                            ξ_initial = prob.ξ_initial, tspan = prob.tspan,
                            output_times = prob.output_times)
end
y(θ) = solve(p_inject(problem, θ); reltol = 1e-7, abstol = 1e-9).extras.total_mass[end]

h = 2e-3                                   # relative step
dy = [(y(θ0 .* (1 .+ h .* (1:2 .== k))) - y(θ0 .* (1 .- h .* (1:2 .== k)))) /
      (2h * θ0[k]) for k in 1:2]           # central differences ∂y/∂θ_k
```

| Choice | Guidance |
|---|---|
| step and tolerance | choose them per case by a convergence check: central differences at relative steps of 1, 2, 4 and 8e-3 at two tolerances (reltol 1e-7 and 1e-8, `abstol` scaled with it), and take a value where the Richardson extrapolants of successive steps agree within 5 % and the two tolerances agree (Chapter 13, §13.4) |
| parallel | `parameter_sweep` over the perturbed parameter vectors (Chapter 13, §13.5) |
| events | a deactivation's time is a smooth function of the parameters, but an objective that samples the solution at fixed times across it is rough on the scale of one sample interval; the convergence check excludes such steps (Technical Reference §17) |
| qualify | transport support changes, endpoint blends and local interface roots before differentiating outputs |

---

## 16.11 Deactivation and the correction

See Chapter 10 and Technical Reference §11 and §14. The correction types are
top-level exports; deactivation is switched by `solve` keywords.

### 16.11.1 Keywords

| Keyword | Default | Effect |
|---|---|---|
| `deactivation` | `false` | deactivate cells whose condensed mass falls below `threshold` of their initial condensed mass |
| `threshold` | `1e-2` | `ε_m`; a cell that starts without condensed mass never deactivates |
| `band` | `max(threshold, 0.1)` when the equidistribution correction moves a layer, else `threshold` | `ε_band` in `[threshold, 1)`, read on a moving mesh |
| `located_events` | `true` | locate each event inside the step at the first crossing |
| `correction` | the rule of §16.6 | `NoCorrection()` or `EquidistributionCorrection(...)` |
| `lagrangian` | the layer interfaces | node ids that always move with the solid |
| `use_ale` | `false` | a moving mesh, which the correction requires |
| `min_thickness` | `1e-6` m | the thickness stop of a moving mesh |

How a deactivation proceeds: cells go from the exposed (`:top`) side only. On a
moving mesh the run of live cells from the outermost one inward whose condensed
mass is below `ε_band` goes, when it holds a cell below `ε_m`, in one merge into
the live cell beneath it, which takes their inventories and volume; nothing is
lost. On a fixed mesh the outermost live cell goes while it is below `ε_m`: the
receiver keeps its volume and takes the condensed inventories, and the dead
cell's gas is vented and recorded on the event. A dead cell's rows become
identity rows, the boundary conditions move to the face of the outermost live
cell, and the state length and the Jacobian's pattern never change. With
`located_events` the event is located at the first time a live cell's condensed
mass reaches its threshold and applied there, the threshold and the band
widened by `reltol`; at the end of every step the cells below the threshold go
in one event, the only events without location. The outcomes are those of
§16.7.1. The supported envelope of these parameters is Chapter 10, §10.9; for a
layer of mixed composition use `threshold = 3e-3`.

### 16.11.2 The correction types

```julia
NoCorrection()

EquidistributionCorrection(; τ = ActivationTime(), activation = (0.3, 0.1),
                             monitor = (1.0, 0.0), layers = :vanishing)

ActivationTime(; λ = 1.8, fallback = 720.0, floor = 1.0)
```

| Name and keyword | Meaning | Default |
|---|---|---|
| `NoCorrection()` | the mesh moves with the solid, `c = 0` | — |
| `EquidistributionCorrection` `τ` | the layer's equidistribution time `τ_eq` [s]: an `ActivationTime` rule, or one positive number for every layer | `ActivationTime()` |
| `activation` | `(r_hi, r_lo)`, `0 < r_lo < r_hi`: a cell's switch is zero while it keeps at least `r_hi` of its initial thickness and one at and below `r_lo` | `(0.3, 0.1)` |
| `monitor` | `(α, β)`, `α > 0`, `β ≥ 0`: the monitor `α + β ρ/ρ_ref`, normalized over the layer | `(1.0, 0.0)` |
| `layers` | `:vanishing` moves only the layers deactivation can remove; `:all` every layer | `:vanishing` |
| `ActivationTime` `λ` | `λ_τ`: the layer's `τ_eq` is `λ_τ` times its thickness `L_layer` over its thinning rate, both taken at its first activation | `1.8` |
| `fallback` | the time [s] a layer keeps while it thins too slowly to estimate (or thickens) | `720.0` |
| `floor` | the thinning below which the estimate waits: `floor` times the initial thickness over the solve's horizon | `1.0` |

The correction requires `use_ale = true`. Under an equidistribution correction
every layer interface (an interior node between cells that start with different
condensed components) must be a Lagrangian node; the default `lagrangian` lists
exactly those nodes. A layer deactivation cannot remove keeps its mesh with the
solid unless `layers = :all`.

```julia
# The installed correction, written out:
sol = solve(problem; use_ale = true, deactivation = true,
            correction = EquidistributionCorrection(τ = ActivationTime(λ = 1.8)))

# The mesh with the solid (no correction) and a band of the threshold's width:
sol = solve(problem; use_ale = true, deactivation = true,
            correction = NoCorrection(), band = 1e-2)
```

### 16.11.3 Constants and building blocks

| Name | Value | Submodule |
|---|---|---|
| `CORRECTION_BAND` | `0.1`, the default band under a moving correction | `Pyrolysis.Integration` |
| `ACTIVATION_TIME_FACTOR` | `1.8`, the default `λ_τ` | `Pyrolysis.Closures` |
| `EQUIDISTRIBUTION_TIME` | `720.0` s, the default `fallback` | `Pyrolysis.Closures` |
| `DT_AFTER_EVENT` | `1e-6` s, the `:forced_dt` restart step | `Pyrolysis.Integration` |
| `residues(material)` | per component, the condensed mass fraction left after complete conversion | `Pyrolysis.Closures` |
| `top_live_cell(layout)` | the outermost live cell (zero when none) | `Pyrolysis.Deactivation` |
| `replay_events!(cells, layout, t)` | set the layout's mask to the one a state saved at `t` had, `cells = cp.deactivation` | `Pyrolysis.Deactivation` |
| `DeactivationEvent` | one deactivated cell (§16.5.4) | `Pyrolysis.Deactivation` |

---

## 16.12 Conservation

See Chapter 9 (§9.8) and Technical Reference §16. `ConservationRecorder` and
`conservation_audit` are top-level exports.

```julia
ConservationRecorder()                       # pass as solve's `callback`
conservation_audit(sol[, recorder]) -> (; method, mass, energy)
```

`ConservationRecorder()` is a read-only observer: over each accepted step it
integrates the live region's mass outflow and enthalpy rate by three-point
Gauss–Legendre quadrature on the step's interpolant (up to the time where an
event ends the step early), and it records the live inventories at the start
and after every step. Its fields are `mass_outflow` [kg], `energy_inflow` [J],
`initial` and `final` (each `(energy, mass)`) and `steps`. A step a callback of
yours ends early cannot be accounted for; its integrals become `NaN`.

`conservation_audit` returns `method` (`:recorder`, or `:trapezoid` without a
recorder, the trapezoid rule over the saved series from the first to the last
saved state) and, for `mass` and `energy` each, the NamedTuple:

| Field | Meaning (mass in kg, energy in J) |
|---|---|
| `initial`, `final` | the live inventories at the start and the end |
| `change` | `final − initial` |
| `inflow` | the time integral of the rate entering the live region: the negated mass outflow; the enthalpy rate |
| `lost` | what the deactivation events took out: their lost masses; their lost enthalpy and caloric regularization |
| `defect` | `change − inflow + lost`, zero for an exact ledger |
| `relative` | `defect` over the largest magnitude among `change`, `inflow` and `lost` |

```julia
recorder = ConservationRecorder()
sol = solve(problem; use_ale = true, deactivation = true, callback = recorder)
audit = conservation_audit(sol, recorder)
audit.mass.relative, audit.energy.relative
```

The series behind the trapezoid are `sol.extras.live_mass`, `live_energy`,
`mass_outflow` and `energy_rate` (§16.5.4). The masses are the states, so
mass is conserved exactly whatever the tolerance; the energy ledger closes in
both energy forms (Technical Reference §16).

The residual itself is reachable for custom loops:
`Pyrolysis.Assembly.residual!(du, u, cp.def, t, cp.ws)` writes the row
derivatives of the state `u` at time `t` into `du`, and fills the derived state
`cp.ws.derived` (`T`, `ξ`, `m`, `V`, `x`, …); `Pyrolysis.Outputs.ledger_rates`
and `live_inventories` evaluate the ledger's rates and inventories from it
(Chapter 11, §11.8).

---

## 16.13 Default-value reference

The defaults built into the package, in one place. The physics is in the noted
Technical Reference chapter.

### 16.13.1 `solve`, `setup` and `integrate`

| Keyword | Default |
|---|---|
| `energy_state` | `:temperature` |
| `abstol` / `reltol` | `1e-8` (K) / `1e-6` |
| `pressure_tolerance` | `1.0` Pa |
| `threshold` | `1e-2` |
| `maxiters` | `1_000_000` |
| `dt` | `problem.dt_initial` |
| `adaptive` | `true` |
| `dtmin` / `dtmax` | `1e-12` / `Inf` |
| `saveat` | `problem.output_times` |
| `callback` | `nothing` |
| `progress` | `false` |
| `validate` | `true` |
| `use_ale` | `false` |
| `w_override` | `nothing` |
| `min_thickness` | `1e-6` m |
| `radiation_model` | `NO_RADIATION` |
| `linear_solver` | `nothing`, which is `StructuredKLU()` |
| `integrator_policy` | `nothing`, which is `default_policy(linear_solver)`, `Rodas5P` with the heuristic restart |
| `sources` | `NoSources()` |
| `deactivation` | `false` |
| `band` | `max(threshold, 0.1)` under a correction that moves a layer, else `threshold` |
| `correction` | `EquidistributionCorrection()` where the mesh moves, deactivation is on and a cell can vanish; else `NoCorrection()` |
| `lagrangian` | the layer interfaces under the equidistribution correction; else none |
| `located_events` | `true` |
| `setup`'s `time_scale` | `problem.dt_initial` |
| `integrate`'s `max_constraint_rejections` | `50` |

### 16.13.2 Components, materials and reactions

| Quantity | Default |
|---|---|
| `SolidComponent` `ε`, `α`, `γ`, `λ`, `κ` | `0.9`, `0.0`, `1.0`, `nothing`, `0.0` |
| `LiquidComponent` `ε`, `α`, `γ`, `λ`, `κ` | `0.95`, `0.0`, `0.0`, `nothing`, `0.0` |
| `GaseousComponent` `ε`, `α`, `γ` | `0.0`, `0.0`, `0.0` |
| `absorptivity` (solid, liquid) | `nothing` (Kirchhoff, `a_j = ε_j(T)`) |
| `ρ_intrinsic` | `= ρ` (no pure-phase pores) |
| `validate_properties` | `true` |
| `Reaction` `n`, `T_min`, `T_max`, `gate_width`, `validate_mass` | `1.0`, `0.0` K, `Inf`, `1.0` K, `true` |
| `MixingSpec` `β`, `length_scale` | `0.5`, `NaN` |
| `MaterialMixing` `k`, `λ`, `κ` | `WEIGHTED`, `PARALLEL`, `WEIGHTED`, each `β = 0.5` |
| `Material` `reactions`, `lateral_shrinkage_law` | `()`, `nothing` |
| `Material` `depletion_limiter` | `DepletionLimiter()`: `threshold = 1.0` kg·m⁻³, `join = 1.0`, `base_scale = 1e-8` kg·m⁻³, `enabled = true` |
| `Material` `gas_viscosity` | `SutherlandProperty()`: `μ_ref = 1.716e-5` Pa·s, `T_ref = 273.15` K, `S = 110.4` K |
| `ArrheniusProperty` `T_ref` | `298.15` K |

### 16.13.3 Boundary conditions

| Quantity | Default |
|---|---|
| `ConvectiveBC` `blowing`, `validate` | `NoBlowing()`, `true` |
| `RadiativeBC` `ε`, `F` | `nothing` (the effective emissivity), `1.0` |
| `RadiativeFluxBC` `absorptivity` | `nothing` (the effective absorptivity, the effective emissivity unless a component sets `absorptivity`) |
| `BoundaryConditionSet` `thermal` / `mass` | `(top = AdiabaticBC(), bottom = AdiabaticBC())` / `(top = ImpermeableBC(), bottom = ImpermeableBC())` |
| `background_gas_concentration` `P` | `101325.0` Pa |
| `AmbientGasBC` `film` | `nothing` |
| `LumpedTransfer` `endpoint`, `blend` | `nothing`, `(1e-4,1e-3)` |
| `CouetteFilm` `λ` | `0.5` |
| `PropertyRatioFilm` `n_Γ` | `0.6` |
| `PoolingFilm` `n_Γ`, `n_D`, `ψ0`, `Λc`, `Λw` | `0.6`, `0.5`, `0.59`, `0.35`, `0.10` |
| `GasSpecies` `μ`, `V_D`, `D`, `T_range` | `nothing` (Eucken), `NaN`, `nothing` (Fuller), `(200.0, 3000.0)` K |

### 16.13.4 Problem and mesh

| Quantity | Default |
|---|---|
| `PyrolysisProblem` `T_initial` | `z -> 300.0` |
| `PyrolysisProblem` `ξ_initial` | the first component, if a solid, at its density at 300 K; the others zero |
| `PyrolysisProblem` `dt_initial` | `1e-4` s |
| `PyrolysisProblem` `output_times` | `min(101, max(11, round((t_end − t_start)/dt_initial/100)))` evenly spaced times over `tspan` |
| mesh `coordinates` | `Planar()` |
| mesh `cross_section` | `1.0` (planar, m²), `1.0` (cylindrical, m), `4π` (spherical) |
| `create_stretched_mesh` `stretch_from` | `:top` |

### 16.13.5 Correction, integrator and linear solver

| Quantity | Default |
|---|---|
| `EquidistributionCorrection` `τ`, `activation`, `monitor`, `layers` | `ActivationTime()`, `(0.3, 0.1)`, `(1.0, 0.0)`, `:vanishing` |
| `ActivationTime` `λ`, `fallback`, `floor` | `1.8`, `720.0` s, `1.0` |
| `CORRECTION_BAND` | `0.1` |
| `rodas5p_policy` `restart` | `:heuristic` |
| `DT_AFTER_EVENT` (`:forced_dt` restart) | `1e-6` s |
| `FAILURE_MEMORY_GROWTH`, `FAILURE_MEMORY_STEPS` | `1.25`, `20` |
| `MAX_CONSTRAINT_REJECTIONS` | `50` |
| `StructuredKLU` `inner` | `KLUFactorization()` |
| colors per dual-number batch of the Jacobian | up to `24` |

### 16.13.6 Kernel constants and limits

| Quantity | Value | Name (submodule) | TR |
|---|---|---|---|
| property bracket: law argument clamp, constructor validation range, the trial temperature of the kernels | `[200, 3000]` K | `T_PROPERTY_MIN`, `T_PROPERTY_MAX` (`Materials`) | §4 |
| width of the smooth trial-temperature clamp at the bracket ends | `1.0` K | `TRIAL_TEMPERATURE_WIDTH` (`Materials`) | §12 |
| surface-temperature Newton tolerance / iterations | `1e-10` / `50` | keywords of `solve_boundary_temperature!` (`BoundaryConditions`) | §12 |
| `DirichletThermalBC` conductance inside a combination | `1e8` W·m⁻²·K⁻¹ | — | §12 |
| pore-fraction divisor floor (both gas models) | `1e-12` | `POROSITY_FLOOR` (`Physics`) | §9 |
| admissible pore margin | raw porosity greater than twice the divisor floor plus mass/thermal/geometry uncertainty | initial, proposed and post-event checks | §9 |
| smooth floor on a cell's heat capacity `Σ_j m_ij c_p,j` | `RHO_CP_FLOOR · Vᵢ(0)`, `RHO_CP_FLOOR = 1e3` J·m⁻³·K⁻¹ | `RHO_CP_FLOOR` (`Materials`) | §5 |
| cell-conductivity floor | `1e-6` W·m⁻¹·K⁻¹ | `CONDUCTIVITY_FLOOR` (`States`) | §5 |
| cell-volume floor | `1e-8 · Vᵢ(0)` | `VOLUME_FLOOR_FRACTION` (`States`) | §13 |
| gas-viscosity floor | `1e-12` Pa·s | `MU_FLOOR` (`Properties`) | §5 |
| temperature floor of the Chapman–Enskog diffusivity (the blowing models) | `100` K | `T_GAS_PROPERTY_FLOOR` (`Properties`) | §12 |
| thermal-property mixing-weight blend window (matrix volume fraction) | `1e-4 … 1e-3` | `MIXING_BLEND_LO`, `MIXING_BLEND_HI` (`Properties`) | §5 |
| gas-transport endpoint blend window (condensed matrix occupancy) | `(1e-4,1e-3)` by default | `LumpedTransfer` / `DustyGas` keyword `blend` (`Materials`) | §9 |
| Kozeny constant | `180` | — | §5 |
| clamp band of the column conversion before the lateral law | `1e-3` outside `[0, 1]` | `LATERAL_CLAMP_WIDTH` (`Closures`) | §10 |
| reference temperature `T_ref` | `298.15` K | `T_REF` (`Materials`) | §2 |
| reference pressure `P_ref` | `101325` Pa | `P_REF` (`Materials`) | §2 |
| gas constant `R_g` | `8.31446261815324` J·mol⁻¹·K⁻¹ | `R_GAS` (`Materials`) | §2 |
| Stefan–Boltzmann constant `σ` | `5.670374419e-8` W·m⁻²·K⁻⁴ | `STEFAN_BOLTZMANN` (`Materials`) | §2 |

---

## 16.14 Quick import cheat-sheet

| You want | Reach it through |
|---|---|
| components, `Material`, `Reaction`, `DepletionLimiter`, property types, mixing rules, `RadiationModel` | top level (`using Pyrolysis`) |
| condition types, `BoundaryConditionSet`, blowing models, `GasSpecies`, `background_gas_concentration` | top level |
| `create_uniform_mesh`, `create_stretched_mesh`, `create_mesh_from_nodes`, `Planar`, `Cylindrical`, `Spherical` | top level |
| `PyrolysisProblem`, `solve`, `setup`, `integrate`, `CoreProblem`, `parameter_sweep` | top level |
| `PyrolysisSolution`, `SolverStatistics`, the accessors, `ConservationRecorder`, `conservation_audit` | top level |
| `NoCorrection`, `EquidistributionCorrection`, `ActivationTime`, the policies, `StructuredKLU` | top level |
| `UMFPACKFactorization`, `KLUFactorization` | `LinearSolve` (`using LinearSolve: …`) |
| `validate_problem`, `is_valid_problem`, `initial_fields`, the profile helpers | `Pyrolysis.Problem` |
| `n_cells`, `n_nodes`, `n_faces`, `boundary_tags`, `boundary_faces`, `geometry!`, `box_mesh`, `cell_mesh` | `Pyrolysis.Meshes` |
| `n_components`, the physical constants, `T_PROPERTY_MIN`, `T_PROPERTY_MAX` | `Pyrolysis.Materials` |
| `SolutionExtras`, `ledger_rates`, `live_inventories` | `Pyrolysis.Outputs` |
| `DeactivationEvent`, `replay_events!`, `top_live_cell` | `Pyrolysis.Deactivation` |
| `residual!`, `VolumetricSources`, `NoSources` | `Pyrolysis.Assembly` |
| `jac_prototype`, `jacobian!`, `dense_jacobian`, `compare_entries`, `missing_entries` | `Pyrolysis.ColoredJacobian` |
| the state layout and its index functions | `Pyrolysis.States` |
| `IntegratorPolicy`, `consistent_state!`, `DomainCheck`, `CORRECTION_BAND` | `Pyrolysis.Integration` |
| everything every submodule exports, at once | `using Pyrolysis.Internal` |
