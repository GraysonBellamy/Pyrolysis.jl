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
> geometry accessors (`n_cells`, `cell_size`, `domain_bounds`, …) are exported
> by their submodules but **not** at the top level. Use `using Pyrolysis.Internal`
> or a submodule-qualified import to call them.

### 16.1.1 Submodule map

The package is a layered DAG of submodules. The table lists each submodule and
the kind of names it owns. Everything a submodule exports is reachable as
`Pyrolysis.<Submodule>.X` and is bundled into `Pyrolysis.Internal`.

| Submodule | Owns |
|---|---|
| `Materials` | `Material`, components, `Reaction`/`DepletionLimiter`, property-function types, `MixingRule`/`MixingSpec`/`MaterialMixing`, `RadiationModel` enum, physical constants |
| `BoundaryConditions` | all thermal/mass/pressure BC types, `BoundaryConditionSet`, `ambient_concentration` |
| `Geometry` | `Mesh1D`, `Node1D`/`Cell1D`/`Face1D`/`CellState`, `create_uniform_mesh`/`create_stretched_mesh`/`create_mesh_from_nodes`, geometry accessors (`n_cells`, `n_nodes`, `n_faces`, `cell_size`, `domain_bounds`, …); `MeshTopology`/`MeshGeometry` are placeholders for a multi-dimensional extension (the 1D residual reads `Mesh1D` and the `StateCache`) |
| `Properties` | mixture/effective-property functions, `PropertyCache`, `RadiationCache` |
| `Physics` | kinetics, conduction, gas transport, radiation (surface/Beer–Lambert), volume change, ideal-gas pressure |
| `Adaptivity` | depleted-cell parking (`park_cell!`, the conservative temperature merge), `find_*` helpers; standalone error indicators, monitors/`relocate_nodes!`, quality metrics |
| `Discretization` | finite-volume operators, flux limiters, boundary fluxes, `solve_boundary_temperature`, `get_surface_temperature`/`get_back_face_temperature`, ALE/Darcy fluxes |
| `Problem` | `PyrolysisProblem`, `PyrolysisSolution`, `ProblemDef`, `Workspace`, `SimulationMode` traits, `StateLayout`, validation, initial-condition helpers |
| `Residual` | `residual!`, column-average helpers |
| `Jacobian` | `JacobianBackend`/`JacobianSpec`, `Structured`/`DenseAD_ForwardDiff`, `LinearSolveStrategy` family, operators, structured couplings, `verify_jacobian_at` |
| `Diagnostics` | conservation trackers, `check_*`/`*_report`, `DiagnosticsCallback`, `DiagnosticsState`/`DiagnosticsRow` |
| `Solver` | `solve`, `parameter_sweep`, callbacks, ILU solvers, `get_profile`/`get_cell_positions`/`interpolate_solution` |
| `Experimental` | P1 radiation (`P1_QUASI_STEADY`, Marshak/Reflective/Dirichlet radiation BCs, the P1 residual), h-refinement/coarsening (`AMRController`; its meshes cannot be passed to `solve`) — **unstable**, not in `Internal` |

### 16.1.2 Extensions

Linear solvers need no extension: pass any LinearSolve factorization (for
example `KLUFactorization()` or `MKLPardisoFactorize()`) directly (§16.8).
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
               ε = 0.9, absorptivity = nothing, α = 0.0, γ = 1.0, λ = 0.0, κ = 0.0,
               ρ_intrinsic = nothing, φ = nothing, validate_properties = true)

LiquidComponent(name::Symbol; ρ, c, k,
                λ = 2e-5, κ = 0.0, ε = 0.95, absorptivity = nothing, α = 0.0, γ = 0.0,
                ρ_intrinsic = nothing, φ = nothing, validate_properties = true)

GaseousComponent(name::Symbol; M, c, k, λ, ε = 0.0, α = 0.0, γ = 0.0,
                 validate_properties = true)
```

| Keyword | Meaning (symbol) | Units | Default |
|---|---|---|---|
| `ρ` | bulk (pure-phase) density `ρ_j` | kg·m⁻³ | required (solid/liquid) |
| `M` | molar mass `M_j` (gas only; sets ideal-gas density `ρ_g = P_ref M/(R_g T)`) | kg·mol⁻¹ | required (gas) |
| `c` | specific heat `c_{p,j}` | J·kg⁻¹·K⁻¹ | required |
| `k` | thermal conductivity `k_j` | W·m⁻¹·K⁻¹ | required |
| `λ` | gas-transfer (diffusion) coefficient `λ_j` | m²·s⁻¹ | `0.0` (solid), `2e-5` (liquid), required (gas) |
| `κ` | permeability `κ_j` (solid/liquid only) | m² | `0.0` |
| `ε` | emissivity `ε_j(T)` (number or property function) | – | `0.9` (solid), `0.95` (liquid), `0.0` (gas) |
| `absorptivity` | surface absorptivity `a_j` for incident radiation (solid/liquid only; `nothing` = Kirchhoff `ε_j(T)`) | – | `nothing` |
| `α` | mass absorption coefficient `α_j` | m²·kg⁻¹ | `0.0` |
| `γ` | swelling factor `γ_j` | – | `1.0` (solid), `0.0` (liquid/gas) |
| `validate_properties` | sample the laws on `[200, 3000]` K at 1 K steps and reject `ρ`, `c`, `k` ≤ 0, `λ`, `κ` < 0 or `ε` outside `[0, 1]` | – | `true` |
| `ρ_intrinsic` | skeletal/intrinsic density `ρ_{i,j}` (porosity formula only) | kg·m⁻³ | `= ρ` |
| `φ` | pure-phase porosity; sets `ρ_intrinsic = ρ/(1−φ)` (constant `ρ` only) | – | `nothing` |

`ρ`, `c`, `k` (and `λ` of a liquid or gas) accept either a plain number or a
property-function object (§16.2.2); they are converted with `to_property`
automatically. `κ`, and `λ` of a `SolidComponent`, take a number. Pass
**either** `ρ_intrinsic` or `φ`, never both — supplying both throws an
`ArgumentError`. The bulk density `ρ` is what mixing rules and mixture density
use; `ρ_intrinsic` is used **only** in the porosity formula
`φ = clamp(1 − Σ ξ_j/ρ_{i,j}, 0, 1)` (see Technical Reference §5.3).

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
state-dependent one); the Jacobian assembly differentiates them with ForwardDiff
(UG §3.3.3; Technical Reference §4.5).

| Type | Constructor | Value `p(T)` |
|---|---|---|
| `ConstantProperty` | `ConstantProperty(value)` | `value` |
| `LinearProperty` | `LinearProperty(a, b)` | `a + b·T` |
| `PolynomialProperty` | `PolynomialProperty((c0, c1, …))` | `evalpoly(T, coeffs)` (low-to-high order) |
| `ArrheniusProperty` | `ArrheniusProperty(A, E, T_ref)` | `A·exp(−E/(R_g·T))` |
| `TableProperty` | `TableProperty((T1,…), (v1,…))` | linear interpolation; nearest-value extrapolation |
| `CallableProperty` | `CallableProperty(f)` | `f(T)` |
| `SutherlandProperty` | `SutherlandProperty(μ_ref, T_ref, S)`; `SutherlandProperty()` = air | `μ_ref·(T/T_ref)^1.5·(T_ref + S)/(T + S)` |
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
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0, validate_mass = true)

# two products, numeric indices
Reaction(name, reactant::Int => (p1::Int, p2::Int);
         A, E, h, yields::Tuple, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0,
         validate_mass = true)

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
| `gate_width` | width `w` of the smooth tanh temperature gates | K | `1.0` |
| `validate_mass` | enforce mass balance at construction | – | `true` |

Mass-balance convention: the single reactant is consumed with implicit mass 1.0;
product yields must sum to 1.0 (the one-product form uses yield 1.0
automatically). Each yield must additionally be ≥ 0 and the order must satisfy
`n ≥ 0` — both enforced unconditionally at construction; `validate_mass = false`
disables only the sum check. `h` accepts a number (`ConstantProperty`), a `(a, b)` tuple
(`LinearProperty`, i.e. `a + b·T`), or any property object — use
`StateDependentProperty` for moisture sorption. See Technical Reference §6 for the rate law, tanh
gates (≈10 K leakage caveat), and the depletion limiter.

```julia
# endothermic charring: virgin -> 0.2 char + 0.8 gas
rxn = Reaction(:charring, :virgin => (:char, :gas);
               A = 1e15, E = 200e3, h = 500e3, yields = (0.2, 0.8),
               T_min = 473.0, T_max = 1273.0)

# exothermic char oxidation (negative h releases heat)
rxn2 = Reaction(:char_oxidation, 1 => 2; A = 1e10, E = 160e3, h = -25e6)
```

Reactions are passed to `Material` as a plain tuple.

### 16.2.4 Mixing rules: `MixingRule`, `MixingSpec`, `MaterialMixing`

Effective conductivity `k`, gas transfer `λ`, and permeability `κ` are each
mixed with their own rule (see Technical Reference §5). The rule enum (top-level exported) is:

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
           depletion_limiter::DepletionLimiter = DepletionLimiter(),
           gas_viscosity = SutherlandProperty())
```

| Keyword | Meaning | Default |
|---|---|---|
| `name` | material identifier | required |
| `components` | tuple of component objects (heterogeneous; order defines index `j`) | required |
| `reactions` | tuple of `Reaction`/`PendingReaction` (symbolic names resolved here) | `()` |
| `mixing` | `(k, λ, κ)` NamedTuple from `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `χ̄ ↦ A(t)/A_0`; `nothing` ⇒ constant cross-section | `nothing` |
| `depletion_limiter` | per-material `DepletionLimiter` for the reaction-rate depletion roll-off (see §16.2.3 note and Technical Reference §6.3) | `DepletionLimiter()` |
| `gas_viscosity` | pore-gas dynamic viscosity μ(T) used in the Darcy mobility `κ/μ` (Darcy–Fick mode only) | `SutherlandProperty()` (air) |

Providing a non-`nothing` `lateral_shrinkage_law` switches the simulation into
the variable-cross-section (WithChi) mode under ALE (see §16.6 and Technical
Reference §10). The law
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
reactions at ξ = 0. Top-level exported (and also exported by the `Physics` and
`Materials` submodules).

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
the mesh creators and the topology/geometry placeholder types are top-level
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

create_mesh_from_nodes(node_positions::AbstractVector{<:Real}, ::Val{NC};
                       T_initial = 300.0, ξ_initial = nothing,
                       cross_section_area = 1.0)
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
using Pyrolysis.Geometry: n_cells, n_nodes, n_faces, cell_size, domain_bounds

n = n_cells(mesh)               # number of cells
dz = [cell_size(mesh, i) for i in 1:n]
zc = mesh.cell_centers          # cached cell centers (a Mesh1D field)
(zmin, zmax) = domain_bounds(mesh)
```

The cached geometry is held in `Mesh1D` fields: `cell_centers`, `cell_volumes`,
`face_positions`, `face_areas`, `cross_section_area` (the current area `A`), and
`boundary_faces` (boundary id → face
indices). `MeshTopology` and `MeshGeometry` *are* top-level exported, but they are
placeholders for a multi-dimensional extension: the 1D residual reads `Mesh1D`
and the `StateCache`, not them. The
component count of a material is `n_components(material)`, exported from
`Pyrolysis.Materials` (also in `Internal`).

### 16.3.3 State-vector layout

The flat ODE state is packed block-major as `[T | ξ_1..ξ_NC | z | χ]`, where the
`z` block exists only in ALE mode and the `χ` block only in WithChi mode (see §16.6
and Technical Reference §15.1). `StateLayout` (top-level exported) is the descriptor; `state_length`,
`pack_state!`, `create_state_vector`, `update_mesh_from_state!`, `block_view`,
and the accessors `get_temperature`, `get_concentration`, `get_node_position`,
`get_node_positions`, `get_thickness_from_state`, `column_chi_bar_from_state`,
`cross_section_area_from_state` are also top-level exported for users who index
raw state. Most users never touch these — `solve` and the solution object handle
packing for you. See §16.12 for the low-level residual surface.

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
| `ConvectiveBC` | `ConvectiveBC(; h, T_inf, blowing = NoBlowing())` | `h(t)·(T_∞ − T_s)`, `h` reduced by the `blowing` model |
| `RadiativeBC` | `RadiativeBC(; ε = nothing, T_inf, F = 1.0)` | `F·ε_eff·σ·(T_∞⁴ − T_s⁴)` |
| `DirichletThermalBC` | `DirichletThermalBC(T)` | enforces `T_s = T` (penalty in Newton) |
| `RadiativeFluxBC` | `RadiativeFluxBC(flux; absorptivity = nothing)` or `RadiativeFluxBC(; flux, absorptivity = nothing)` | incident external radiation (model-dependent, below) |
| `CombinedThermalBC` | `bc1 + bc2 + …` | sum of components |

Notes and defaults:

- `ConvectiveBC` requires `h` and `T_inf` keywords; both accept a number or a
  function of `t`.
- `blowing` takes a blowing model — `NoBlowing()` (default), `CouetteFilm`,
  `PropertyRatioFilm`, `PoolingFilm`, `FunctionBlowing`, all top-level exported —
  that reduces `h` with the outward gas flux at the face (§5.2.3).
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
tol `1e-10`, ≤50 iterations) — see Technical Reference §12.2.

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
flow (any solid/liquid `κ > 0`); see Technical Reference §9.

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

> The boundary energy balance has no surface gas-enthalpy term: `S_conv`
> carries the gas enthalpy volumetrically, and a surface term would
> double-count it (§5.4).

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

`validate_problem(problem; strict = false)` (in `Pyrolysis.Problem`) checks the
mesh, the material's component count and reaction indices, the boundaries, the
time span, the initial-condition count, and (Darcy–Fick) the initial porosity;
`strict = true` adds the stoichiometry, initial-condition-call and positive-volume
checks (§7.9). `solve` runs it
automatically unless `validate = false`.

### 16.5.2 `PyrolysisSolution`

`solve` returns a `PyrolysisSolution`. Core fields (all read directly):

| Field | Type | Meaning |
|---|---|---|
| `problem` | `PyrolysisProblem` | the copy of the problem the run integrated (its mesh holds the final state) |
| `t` | `Vector{Float64}` | saved times [s] |
| `T` | `Matrix{Float64}` | temperatures, `(n_cells, n_steps)` [K] |
| `ξ` | `Array{Float64,3}` | concentrations, `(NC, n_cells, n_steps)` [kg·m⁻³] |
| `z` | `Matrix{Float64}` | node positions, `(n_nodes, n_steps)` [m]; empty `0×0` in Eulerian mode |
| `surface_T` | `Vector{Float64}` | solved surface temperature per step [K] |
| `mass_loss_rate` | `Vector{Float64}` | MLR per **current** area `A(t)` (top + permeable-bottom) [kg·m⁻²·s⁻¹] |
| `retcode` | `Symbol` | integrator return code (e.g. `:Success`, `:Terminated`) |
| `solve_time` | `Float64` | wall-clock solve time [s] |
| `raw_solution` | (untyped) | the underlying ODE solution (interpolation source) |
| `extras` | `NamedTuple` | opt-in diagnostics (below) |

`extras` is a NamedTuple of vectors, written in place during extraction:
`gas_generation_rate`, `back_temperature`, `total_mass`, `thickness_history`,
`surface_mass_flux`, `surface_heat_flux_net`, `surface_heat_flux_absorbed`,
`surface_heat_flux_reradiation`, `surface_heat_flux_convection`,
`surface_heat_flux_conduction`, `surface_blowing_ratio`,
`surface_pooling_factor`, `surface_face_fraction`, `surface_face_buoyancy`,
`volumetric_absorption`, `radiation_transmitted`,
`surface_emissivity`, `surface_conductivity`, `back_heat_flux`, `MLR_gauge`,
`cross_section_area_history`, `live_cells`, `park_times`, `diagnostics_log`
(§9.3).

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
boundary state, not the solution. `solve` leaves `problem` unchanged and works on
a copy, kept as `sol.problem`, whose mesh holds the last accepted step:

```julia
using Pyrolysis.Discretization: get_surface_temperature, get_back_face_temperature,
                                check_boundary_energy_balance
Ts_final    = get_surface_temperature(sol.problem.mesh)
Tback_final = get_back_face_temperature(sol.problem.mesh)
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
| `linear_solver` | factorization / `nothing` | `nothing` → `KLUFactorization()`, or `UMFPACKFactorization()` for the exact ALE Jacobian | any LinearSolve factorization (only paired with a *constructor* `integrator`) |
| `jacobian` | `Structured` / `nothing` | `nothing` → `Structured(strategy = DirectSolve())` (fixed mesh) or `Structured(strategy = ApproxSolve())` (ALE) | Jacobian strategy; `solve` builds the spec |
| `verify_jacobian` | `NamedTuple` / `nothing` | `nothing` | `(rtol=…, atol=…, sample_times=[…])` to compare backend vs. `DenseAD_ForwardDiff` |
| `abstol` | `Real` | `1e-8` | absolute ODE tolerance |
| `reltol` | `Real` | `1e-6` | relative ODE tolerance |
| `maxiters` | `Int` | `1000000` | max integration steps |
| `dtmin` | `Real` | `1e-12` | minimum timestep [s] |
| `dtmax` | `Real` | `Inf` | maximum timestep [s] |
| `saveat` | vector | `problem.output_times` | snapshot times |
| `callback` | callback / `nothing` | `nothing` | extra user callback(s) |
| `progress` | `Bool` | `true` | progress bar |
| `validate` | `Bool` | `true` | run `validate_problem` first |
| `use_ale` | `Bool` | `false` | enable ALE moving mesh |
| `min_thickness` | `Real` | `1e-6` | ALE termination on the live thickness [m] |
| `handle_depletion` | `Bool` | `false` | park depleted cells (ALE: thin cells; fixed mesh: the emptied top cell, and stop at the first interior void) |
| `depletion_threshold` | `Real` | `0.05` | fraction of the initial thickness (ALE) or initial condensed concentration (fixed mesh) below which a cell parks |
| `min_cells` | `Int` | `2` | minimum live cell count before terminating instead of parking |
| `radiation_model` | `RadiationModel` | `NO_RADIATION` | `NO_RADIATION` / `SURFACE_ABSORPTION` / `BEER_LAMBERT` |
| `diagnostics` | `Bool` | `false` | log sensible-energy snapshots in `extras.diagnostics_log` |
| `diagnostics_stride` | `Int` | `10` | snapshot every N accepted steps |
| `warm_start` | `PyrolysisSolution` / `nothing` | `nothing` | reuse a prior solution's final state as initial data |

Key behaviours and constraints:

- **Default integrator / linear solver**: `KenCarp4(autodiff = AutoFiniteDiff(),
  linsolve = KLUFactorization())`, with `UMFPACKFactorization()` for the
  exact ALE Jacobian (`DirectSolve` under ALE). KenCarp4 is an SDIRK method
  suited to the
  moderately stiff Arrhenius system. Pass a constructor (e.g. `integrator = KenCarp4`)
  to combine with `linear_solver`, or a fully-built instance
  (e.g. `integrator = KenCarp4(linsolve = ILUGMRESFactorization())`) on its own.
- **Default Jacobian**: `Structured(strategy = DirectSolve())` on a fixed mesh,
  `Structured(strategy = ApproxSolve())` under ALE.
  Pass a bare `Structured` backend to opt into the other strategy (§16.8); a
  `DenseAD_ForwardDiff()` backend or a prebuilt `JacobianSpec` raises a
  `TypeError`.
- **`radiation_model = P1_QUASI_STEADY`** is rejected (P1 is experimental).
- **`handle_depletion = true`** keeps the state size (parked cells become
  inert ghost slots), so it combines with either strategy; ghost slots read as
  `NaN` in the solution arrays and `solution.extras.live_cells` /
  `park_times` record the events.
- **MLR gauge note**: in variable-cross-section runs, `sol.mass_loss_rate` is per
  current area and `sol.extras.MLR_gauge` is normalized to `A_0`.

```julia
# Cone calorimeter, opaque, ALE surface recession, default Jacobian:
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve: SparspakFactorization
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
| `Structured` | production operator-based backend; `Structured(; strategy = DirectSolve())` — the only backend the `jacobian` keyword accepts |
| `DenseAD_ForwardDiff` | dense whole-residual AD reference — verification only, not for production |
| `JacobianBackend` | abstract supertype |
| `JacobianSpec` | bundles a backend with its sparsity patterns and caches (fields `backend`, `sparsity_phys`, `sparsity_geom`, `cache`) |

Pass the bare backend to `solve`, which builds the spec for its own mesh copy and
mode. Build a spec yourself only to drive `compute_jacobian!` or
`verify_jacobian_at` outside `solve`, from a `(def, ws)` pair:

```julia
sol = solve(problem; use_ale = true, jacobian = Structured(strategy = DirectSolve()))

def, ws, u0 = initialize_state(problem; use_ale = true)
spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
```

### 16.8.2 Linear-solve strategies

`Structured` takes a `LinearSolveStrategy` (both top-level exported):

| Strategy | Behaviour | When |
|---|---|---|
| `DirectSolve` | materialize `A + Σ densify(couplings)`, factor (exact) | fixed-mesh default; exact ALE Jacobian |
| `ApproxSolve` | drop the structured couplings (Schur-style approximation) | ALE default |

`ApproxSolve` is marked approximate via the `is_approximate` trait; the
verifier compares it in masked (in-support) mode. It is the ALE default
(§11.3.2).

### 16.8.3 Iterative linear solvers (not top-level; use `Internal`/`Solver`)

```julia
using Pyrolysis.Solver: ILUGMRESFactorization, ILUBiCGSTABFactorization
using OrdinaryDiffEqSDIRK: KenCarp4

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
```

Or run verification *during* a solve at chosen times:

```julia
sol = solve(problem; verify_jacobian = (rtol = 1e-6, atol = 1e-10,
                                        sample_times = [0.0, 1.0, 10.0]))
```

It warns (does not abort) when the active backend deviates from the dense AD
reference.

### 16.8.5 Other LinearSolve factorizations

Any LinearSolve factorization is passed directly, as `linear_solver` or inside a
built integrator; an explicit choice always replaces the pattern-based default:

```julia
using Pyrolysis
using LinearSolve: KLUFactorization, SparspakFactorization
sol = solve(problem; linear_solver = SparspakFactorization())
using OrdinaryDiffEqSDIRK: KenCarp4
sol = solve(problem; integrator = KenCarp4(linsolve = KLUFactorization()))
```

Intel MKL Pardiso works the same way once Pardiso.jl is installed and loaded:
`using Pardiso; solve(problem; linear_solver = MKLPardisoFactorize())`, with
`MKLPardisoFactorize` from LinearSolve.

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
                       mixing = prob.material.mixing,
                       lateral_shrinkage_law = prob.material.lateral_shrinkage_law,
                       depletion_limiter = prob.material.depletion_limiter,
                       gas_viscosity = prob.material.gas_viscosity)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat,
                            bc_set = prob.bc_set, T_initial = prob.T_initial,
                            ξ_initial = prob.ξ_initial, tspan = prob.tspan,
                            output_times = prob.output_times)
end

extract(sol) = (peak_T = maximum(sol.surface_T),
                final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify = modify, extract = extract)
```

---

## 16.10 Sensitivity analysis

See Chapter 13 and Technical Reference §17. There is no sensitivity API:
parameter derivatives are central differences of `solve` calls on problems that
a closure `p_inject(base, θ)` rebuilds from a parameter vector — the same
pattern as `modify` in §16.9. The solve path is `Float64` throughout, so
`ForwardDiff.Dual` parameters fail at construction and there is no automatic
differentiation through `solve`.

```julia
using Pyrolysis, FiniteDiff

θ0 = [log10(8.5e12), 188e3]   # log10(A), E
p_inject = (prob, θ) -> begin
    rxn = Reaction(:decomp, 1 => 2; A = 10.0^θ[1], E = θ[2], h = 870e3)
    mat = Material(name = :PMMA, components = prob.material.components,
                   reactions = (rxn,), mixing = prob.material.mixing)
    PyrolysisProblem(mesh = prob.mesh, material = mat, bc_set = prob.bc_set,
                     tspan = prob.tspan, output_times = prob.output_times)
end
y(θ) = [solve(p_inject(problem, θ); progress = false).extras.total_mass[end]]

J = FiniteDiff.finite_difference_jacobian(y, θ0, Val(:central); relstep = 1e-2)
```

| Choice | Guidance |
|---|---|
| `relstep` | `1e-3` to `1e-2`; check by halving (§13.4) |
| outputs | integrated (mass, time averages) over pointwise values (§13.4) |
| parallel | `parameter_sweep` over the perturbed vectors (§13.5) |
| avoid | `warm_start`, `handle_depletion = true` (both non-smooth in `θ`) |

---

## 16.11 Adaptivity and callbacks

See Technical Reference §14. Depleted-cell parking is the only mesh adaptation
`solve` performs. The indicators, monitors and quality metrics are standalone
building blocks, and h-refinement/coarsening and the `AMRController` are in
`Pyrolysis.Experimental`; their meshes cannot be passed to `solve`.

### 16.11.1 Surface-depletion via `solve`

The simplest path is the `solve` keywords (§16.7): `handle_depletion = true`,
`depletion_threshold` (default `0.05`), `min_cells` (default `2`), in either
mesh mode.

### 16.11.2 Callbacks (not top-level; use `Internal`/`Solver`)

```julia
using Pyrolysis.Solver: create_depletion_callback, create_thickness_termination_callback,
                        create_step_update_callback, DepletionCallbackState, n_live_cells

DepletionCallbackState(layout, mesh, material, u0)      # live count, per-slot reference, park log
create_thickness_termination_callback(sv_ref, min_thickness::Float64)
create_depletion_callback(ws, mesh, material, cb_state; threshold = 0.05, min_cells = 2)
create_step_update_callback(ws, material)
n_live_cells(cb_state, t)                                 # live cells in the state saved at t
```

You normally do not build these by hand — `solve` assembles them from its
keywords. They are exposed for custom solver loops.

### 16.11.3 Error indicators, monitors, quality metrics (Adaptivity submodule)

For manual r-refinement / mesh-quality work on a mesh you hold (`solve` does not
call these), import from `Pyrolysis.Adaptivity`:

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
and are unstable. A refined or coarsened mesh cannot be passed to `solve`
(refinement leaves inactive parent cells among the active ones, and the kernels
assume the active cells come first); the active-cell bookkeeping requires
iterating the mask, not `eachindex`. Treat this path as research-grade (§10.7).

---

## 16.12 Conservation diagnostics

See Technical Reference §16. Enable snapshots inline with `solve(problem;
diagnostics = true, diagnostics_stride = N)` — rows land in
`sol.extras.diagnostics_log`. The tracker types and check functions are in
`Pyrolysis.Diagnostics` (and in `Internal`).

```julia
using Pyrolysis.Internal   # or: using Pyrolysis.Diagnostics

# trackers and the step-by-step ledger
MassConservationTracker, EnergyConservationTracker
DiagnosticsState, DiagnosticsRow, DiagnosticsCallback
record_substep!, snapshot_diagnostics_row

# instantaneous quantities
compute_domain_mass(mesh, material)             # NTuple of per-component mass [kg]
compute_sensible_energy(mesh, material)         # E_sensible, every component [J]

# checks (return NamedTuples)
check_mass_conservation(tracker, mesh, material; tol = 1e-10)
check_energy_conservation(tracker, mesh, material; tol = 1e-6)

# human-readable reports (return String)
mass_conservation_report(tracker, mesh, material)
energy_conservation_report(tracker, mesh, material)
```

The energy diagnostics keep one ledger — the sensible energy `Σ ρc_p T V` over
every component — that closes to integrator tolerance by construction for constant composition and volume (it audits
bookkeeping) when it is driven step by step through `record_substep!`. A
`DiagnosticsRow` carries `t`, `E_sensible`, `energy_residual`,
`energy_relative`, `mass_max_relative_error`; the rows `solve` logs
(`snapshot_diagnostics_row`) fill only `t` and `E_sensible` and leave the
cumulative fields `NaN` (§9.8).

The low-level residual surface — `ProblemDef`, `Workspace`, `residual!`,
`to_problem_def`, `build_workspace`, `build_ad_workspace`, `StateLayout` — is
top-level exported for users building custom solver or sensitivity loops. The
canonical call is `residual!(du, u, def, t, ws)`; build `(def, ws, u0)` with
`initialize_state(problem; use_ale, radiation_model)` (plus, for direct
Beer–Lambert use, `radiation_intensity` and `radiation_kirchhoff` — set
`radiation_kirchhoff = false` if you have already folded an explicit
absorptivity into the intensity, else the live surface `ε_eff` is applied at
injection), which applies the initial conditions and runs `to_problem_def` and
`build_workspace` with the same keywords. See Technical Reference §15.

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
| default integrator | `KenCarp4(autodiff = AutoFiniteDiff())` |
| default linear solver | `KLUFactorization()`; `UMFPACKFactorization()` for the exact ALE Jacobian (`DirectSolve` under ALE) |
| default Jacobian | `Structured(strategy = DirectSolve())`; `ApproxSolve` under ALE |

### Component / material / reaction defaults

| Quantity | Default |
|---|---|
| `SolidComponent` `ε`, `α`, `γ` | `0.9`, `0.0`, `1.0` |
| `SolidComponent` `λ`, `κ` | `0.0`, `0.0` |
| `SolidComponent`/`LiquidComponent` `absorptivity` | `nothing` (Kirchhoff, `a_j = ε_j(T)`) |
| `LiquidComponent` `λ` | `2e-5` m²·s⁻¹ |
| `LiquidComponent` `ε`, `γ` | `0.95`, `0.0` |
| `GaseousComponent` `ε`, `γ` | `0.0`, `0.0` |
| `ρ_intrinsic` (when omitted) | `= ρ` (no pure-phase pores) |
| `Reaction` `n`, `T_min`, `T_max`, `gate_width`, `validate_mass` | `1.0`, `0.0`, `Inf`, `1.0` K, `true` |
| `MixingSpec` `β` | `0.5` |
| `MaterialMixing` `k`, `λ`, `κ` | `WEIGHTED`, `PARALLEL`, `WEIGHTED` (all `β = 0.5`) |
| `Material` `reactions`, `lateral_shrinkage_law` | `()`, `nothing` |

### BC defaults

| Quantity | Default |
|---|---|
| `RadiativeBC` `ε`, `F` | `nothing` (⇒ `ε_eff`), `1.0` |
| `RadiativeFluxBC` `absorptivity` | `nothing` (⇒ `a_eff`; `= ε_eff` by Kirchhoff unless a component sets `absorptivity`) |
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
| `Structured()` `strategy` (constructor default) | `DirectSolve()` |

### Kernel constants and limits

| Quantity (symbol) | Default | TR |
|---|---|---|
| depletion limiter threshold | `1.0` kg·m⁻³ | §6 |
| property bracket `[T_PROPERTY_MIN, T_PROPERTY_MAX]` (law argument clamp, constructor validation, surface-temperature Newton clamp) | `[200, 3000]` K | §4, §12 |
| mixing-weight blend window (matrix volume fraction) | `1e-4 … 1e-3` | §5 |
| surface-temperature Newton tol / max-iter | `1e-10` / `50` | §12 |
| porosity floor (ideal-gas pressure) | `1e-12` | §9 |
| matrix ρc_p floor (`RHO_CP_FLOOR`) | `1e3` J·m⁻³·K⁻¹ | §5 |
| pore-gas viscosity `Material.gas_viscosity` | `SutherlandProperty()` — air: `μ_ref` `1.716e-5` Pa·s at `T_ref` `273.15` K, `S` `110.4` K | §5 |
| gas-viscosity floor (`MU_FLOOR`) | `1e-12` Pa·s | §5 |
| Kozeny constant | `180` | §5 |
| reference temperature `T_ref` | `298.15` K | §2 |
| reference pressure `P_ref` | `101325` Pa | §2 |
| gas constant `R_g` | `8.31446261815324` J·mol⁻¹·K⁻¹ | §2 |
| Stefan–Boltzmann `σ` | `5.670374419e-8` W·m⁻²·K⁻⁴ | §2 |
| mass conservation tol | `1e-10` | §16 |
| energy conservation tol | `1e-6` | §16 |

---

## 16.14 Quick import cheat-sheet

| You want… | Reach it via |
|---|---|
| `Material`, components, `Reaction`, property types, `MixingSpec`, `RadiationModel` | top level (`using Pyrolysis`) |
| BC types, `BoundaryConditionSet` | top level |
| `create_uniform_mesh`, `create_stretched_mesh`, `create_mesh_from_nodes` | top level |
| `PyrolysisProblem`, `PyrolysisSolution`, `solve`, `parameter_sweep` | top level |
| `JacobianSpec`, `Structured`, `DirectSolve`/`ApproxSolve`, `verify_jacobian_at` | top level |
| `KLUFactorization`, `UMFPACKFactorization`, `SparspakFactorization`, `MKLPardisoFactorize` | `LinearSolve` (`using LinearSolve: …`; Pardiso/MKL for the last) |
| `ProblemDef`, `Workspace`, `residual!`, `StateLayout`, `initialize_state`, `to_problem_def`, `build_workspace` | top level |
| `get_profile`, `get_cell_positions`, `interpolate_solution` | `Pyrolysis.Solver` / `Pyrolysis.Internal` |
| `get_surface_temperature`, `get_back_face_temperature` | `Pyrolysis.Discretization` / `Internal` |
| `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`, depletion callbacks | `Pyrolysis.Solver` / `Internal` |
| geometry accessors (`n_cells`, `cell_size`, `domain_bounds`, …), `n_components` | `Pyrolysis.Geometry` / `Pyrolysis.Materials` / `Internal` |
| conservation trackers, `check_*`, `*_report`, `DiagnosticsCallback` | `Pyrolysis.Diagnostics` / `Internal` |
| trait subtypes (`Eulerian`, `ALE`, `Fickian`, `DarcyFick`, `NoChi`, `WithChi`) | `Pyrolysis.Problem` |
| P1 radiation, h-refinement/coarsening, Marshak/Reflective radiation BCs | `Pyrolysis.Experimental` (unstable) |
| everything every submodule exports, at once | `using Pyrolysis.Internal` |
