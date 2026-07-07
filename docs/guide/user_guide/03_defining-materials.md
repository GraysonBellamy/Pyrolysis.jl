# 3. Defining Materials, Components, and Properties

A `Material` is the physical heart of every Pyrolysis.jl simulation: it bundles the
condensed-phase and gaseous **components** that make up your sample, the
**reactions** that convert them into one another, the **property functions** that
give each component its temperature- (and optionally composition-) dependent
thermophysical behaviour, and the **mixing rules** that homogenize those component
properties into the effective bulk properties the solver actually integrates. This
chapter shows how to build each of these pieces with the real constructors, what
every field and default means, and how to assemble and validate a complete
material. The physics and derivations behind these objects live in the Technical
Reference (see Technical Reference §4 *Materials & Property Functions* and §5
*Effective Properties & Mixing Rules*); here the focus is on *how to specify* them
correctly.

Throughout, the spatial coordinate is `z` (through-thickness): `z = 0` is the
bottom/substrate, `z = L` is the top/exposed surface (see Technical Reference §1).
The primary species state variable is the **mass concentration** `ξ_j` in kg·m⁻³
(not mass fraction). All inputs are SI: densities in kg·m⁻³, heat capacities in
J·kg⁻¹·K⁻¹, conductivities in W·m⁻¹·K⁻¹, activation energies in J·mol⁻¹.

---

## 3.1 Overview: the component → material hierarchy

The build order is always the same:

1. **Components** (`SolidComponent`, `LiquidComponent`, `GaseousComponent`) — each
   carries its phase as a *type*, so the solver dispatches on phase at compile
   time. Solids and liquids form the matrix; gases flow through it.
2. **Property functions** — every per-component thermophysical property (`ρ`, `c`,
   `k`, `λ`, `κ`, intrinsic density) is a callable of temperature. You may pass a
   bare number, a `(a, b)` tuple, a function, or an explicit property-function
   object; the constructors convert for you.
3. **Reactions** (`Reaction`) — single-reactant Arrhenius decomposition steps with
   product yields, activation energy, and heat of reaction. Components may be
   referenced by **symbol** (resolved automatically) or by integer index.
4. **Mixing rules** (`MixingSpec`, `MaterialMixing`) — pick how component thermal
   conductivity `k`, gas-transfer coefficient `λ`, and permeability `κ` are
   homogenized.
5. **Material** (`Material`) — the assembled object, optionally checked with the
   internal helper `Pyrolysis.Materials.validate_material` (or simply validated by
   `solve` when you build the problem; see §3.7).

A minimal end-to-end example (one solid, one gas, one reaction) appears in the
package's top-level docstring and in Chapter 2; this chapter expands every option.

> **A note on what is exported.** `using Pyrolysis` brings the *constructors* into
> scope: `SolidComponent`, `LiquidComponent`, `GaseousComponent`, `Material`,
> `Reaction`, `ReactionSet`, all property-function types, the `MixingRule` enum
> values, `MixingSpec`, and `MaterialMixing`. The various **helper and inspection
> functions** used in this chapter are *not* part of the curated top-level surface,
> so you reach them by qualifying through the submodule that defines them. This works
> the same way whether or not the function is exported from that submodule — Julia's
> `Module.name` syntax resolves any name defined in the module:
>
> - From the `Materials` submodule: `Pyrolysis.Materials.n_components`,
>   `Pyrolysis.Materials.n_reactions`, `Pyrolysis.Materials.dry_pyrolysis_gas_indices`,
>   `Pyrolysis.Materials.gas_component_indices`, and `Pyrolysis.Materials.has_darcy_flow`
>   are *exported* from `Materials`; the validation helpers
>   `Pyrolysis.Materials.validate_material`,
>   `Pyrolysis.Materials.verify_stoichiometry`, and
>   `Pyrolysis.Materials.check_material_mass_balance` are *internal* (not exported, not
>   part of the stable public API) but are still callable through the same qualified
>   path.
> - From the `Properties` submodule: the effective-property functions
>   (`Pyrolysis.Properties.mixture_density`, etc.).
>
> For convenience you may `using Pyrolysis.Internal`, which re-exports the union of
> every submodule's *exports* and so brings the exported helpers
> (`n_components`, `dry_pyrolysis_gas_indices`, `has_darcy_flow`, …) into scope
> unqualified — but **not** the internal validation helpers above, which still
> require the `Pyrolysis.Materials.` prefix. `Internal` is not the stable public
> surface; prefer qualified access in production code. If you want validation as part
> of the public API, build a `PyrolysisProblem` and let `solve` validate (it runs the
> same checks by default, `validate = true`); see §3.7.

---

## 3.2 Defining components

### 3.2.1 `SolidComponent`

A solid is the canonical matrix-forming component. The keyword constructor is:

```julia
SolidComponent(name::Symbol;
               ρ, c, k,
               ε = 0.9, α = 0.0, γ = 1.0,
               λ = 0.0, κ = 0.0,
               ρ_intrinsic = nothing, φ = nothing)
```

| Keyword | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `name` | — | Component identifier (positional, a `Symbol`) | — | required |
| `ρ` | `ρ_j` | **Bulk** (pure-phase) density | kg·m⁻³ | required |
| `c` | `c_{p,j}` | Specific heat capacity | J·kg⁻¹·K⁻¹ | required |
| `k` | `k_j` | Thermal conductivity | W·m⁻¹·K⁻¹ | required |
| `ε` | `ε_j` | Surface emissivity | – | `0.9` |
| `α` | `α_j` | Mass-basis absorption coefficient | m²·kg⁻¹ | `0.0` |
| `γ` | `γ_j` | Swelling factor (volume-change participation) | – | `1.0` |
| `λ` | `λ_j` | Gas-transfer (diffusion) coefficient contribution | m²·s⁻¹ | `0.0` |
| `κ` | `κ_j` | Darcy permeability | m² | `0.0` |
| `ρ_intrinsic` | `ρ_{i,j}` | Skeletal/intrinsic density (porosity only) | kg·m⁻³ | `= ρ` |
| `φ` | — | Pure-phase porosity (alternative to `ρ_intrinsic`) | – | `nothing` |

```julia
using Pyrolysis

# Virgin wood: temperature-rising heat capacity, opaque, tight matrix
virgin = SolidComponent(:virgin,
    ρ = 530.0,            # bulk density [kg/m³]
    c = (1500.0, 1.0),    # LinearProperty: cp(T) = 1500 + 1.0·T  [J/(kg·K)]
    k = 0.12,             # ConstantProperty 0.12 W/(m·K)
    ε = 0.88,             # emissivity
    λ = 2e-5,             # contributes to mixture gas-transfer coefficient
    κ = 1e-15)            # very low permeability (tight cellular structure)

# Char: lighter, more conductive, more emissive, absorbing, porous
char = SolidComponent(:char,
    ρ = 150.0, c = 1000.0, k = 0.10,
    ε = 0.95,
    α = 5.0,              # mass absorption coeff (used only by in-depth radiation)
    λ = 2e-5,
    κ = 1e-12)            # open, porous char conducts gas readily
```

Notes on the less-obvious fields:

- **`ε` and `α` are plain `Float64`s, not property functions.** Emissivity and the
  absorption coefficient do not carry temperature dependence at the component level.
- **`α` (absorption) is ignored unless an in-depth radiation model is active.** With
  `radiation_model = SURFACE_ABSORPTION` (or `NO_RADIATION`) the absorbed flux `α·Φ`
  is applied at the surface through the boundary energy balance and the component `α`
  values play no role — in
  fact `solve` will emit a warning if you set non-zero `α` together with
  `SURFACE_ABSORPTION`. Use `BEER_LAMBERT` (or P1, experimental) to make `α` matter.
  See §3.6 and Technical Reference §8.
- **`λ` on a solid** does not transport the solid (solids are stationary); it is the
  solid's *contribution* to the effective mixture gas-transfer coefficient through
  the mixing rule (see §3.5). ThermaKin-style setups typically use `λ ≈ 2e-5`.
- **`κ = 0.0` (default) means impermeable.** Setting any solid/liquid component's
  `κ > 0` switches the whole problem into Darcy–Fick (pressure-driven) gas transport
  (§3.6.2 and Technical Reference §9).
- **Swelling factor `γ`** is the fraction of the component's volume that counts as
  bulk (matrix) volume, with **one meaning everywhere**: the dilation rate θ, the
  mixture density, and the conductivity/emissivity volume weights all use the same
  `γ_j·ξ_j/ρ_j` weighting. `γ = 1.0` (default for solids) means the solid shrinks as
  it is consumed; `γ = 0.0` makes the component consistently *volumeless* — no volume
  change, no volume in the bulk density, and **no weight in conductivity/emissivity
  mixing** (use with care for a "rigid char"); intermediate values are for
  intumescent/partial-expansion cases. Values outside [0, 1] are rejected at
  construction. See Technical Reference §5.1–5.2 and §10.

### 3.2.2 `LiquidComponent`

Liquids are matrix-forming like solids (same field set) and are the standard way to
represent **bound/free moisture** in hygroscopic materials. The defaults differ from
solids: a liquid by default is strongly emissive (ε = 0.95, near water's ≈ 0.96),
non-swelling, and blocks gas flow. (For the default γ = 0 the emissivity carries no
weight in the mixture anyway — a γ = 0 liquid is volumeless in the emissivity mixing;
the 0.95 matters as soon as you give the liquid γ > 0.)

```julia
LiquidComponent(name::Symbol;
                ρ, c, k,
                λ = 2e-5, κ = 0.0,
                ε = 0.95, α = 0.0, γ = 0.0,
                ρ_intrinsic = nothing, φ = nothing)
```

| Keyword | Meaning | Units | Default |
|---|---|---|---|
| `ρ`, `c`, `k` | bulk density, heat capacity, conductivity | (as solid) | required |
| `λ` | gas-transfer contribution | m²·s⁻¹ | `2e-5` |
| `κ` | permeability | m² | `0.0` (liquid blocks pores) |
| `ε` | emissivity | – | `0.95` (water ≈ 0.96) |
| `α` | mass absorption | m²·kg⁻¹ | `0.0` |
| `γ` | swelling factor | – | `0.0` |
| `ρ_intrinsic` / `φ` | skeletal density / pure-phase porosity | kg·m⁻³ / – | `= ρ` / `nothing` |

```julia
# Liquid water (moisture) in wood
moisture = LiquidComponent(:moisture,
    ρ = 1000.0, c = 4180.0, k = 0.60)   # λ defaults to 2e-5, κ=0 (blocks flow)
```

The default `γ = 0.0` for liquids means moisture evaporation does **not** by itself
drive bulk shrinkage; set `γ = 1.0` if you want liquid loss to shrink the matrix.

### 3.2.3 `GaseousComponent`

Gases flow *through* the matrix; they are the products of decomposition reactions.
A gaseous component has **no permeability field** (`get_permeability` returns 0 for
gases) and carries a **molar mass** instead, which feeds the ideal-gas density and
the pore-pressure closure `P = Σ_g ξ_g R_g T/(M_g φ)`.

```julia
GaseousComponent(name::Symbol;
                 M, c, k, λ,
                 ε = 0.0, α = 0.0, γ = 0.0)
```

| Keyword | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `M` | `M_j` | Molar mass | kg·mol⁻¹ | required |
| `c` | `c_{p,j}` | Heat capacity | J·kg⁻¹·K⁻¹ | required |
| `k` | `k_j` | Thermal conductivity | W·m⁻¹·K⁻¹ | required |
| `λ` | `λ_j` | Gas-transfer (diffusion) coefficient | m²·s⁻¹ | required |
| `ε` | `ε_j` | Emissivity | – | `0.0` |
| `α` | `α_j` | Mass absorption | m²·kg⁻¹ | `0.0` |
| `γ` | `γ_j` | Swelling factor | – | `0.0` |

```julia
# Pyrolysis gas (modelled as a CO₂-like vapour) and water vapour
gas   = GaseousComponent(:gas, M = 0.044, c = 1800.0, k = 0.03, λ = 1e-4)
vapor = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025, λ = 2e-5)
```

**Density is automatic.** Unlike solids/liquids, you do *not* pass `ρ` to a gas. The
constructor installs the ideal-gas law

```
ρ_g(T) = P_ref · M / (R_g · T)
```

with `P_ref = 101325 Pa` and `R_g = 8.31446261815324 J·mol⁻¹·K⁻¹`. For example a
`M = 0.044` gas gives `ρ_g(300) ≈ 1.787 kg·m⁻³`. The ideal-gas default is correct for
almost all pyrolysis applications. See Technical Reference §4.3.

**Naming and water vapour.** Components named so that the name contains `H2O` or
`water` (case-insensitive, e.g. `:H2O`, `:water_vapor`) are treated as *moisture*
gases and are **excluded** from the dry-pyrolysis-gas set that drives lateral
shrinkage `χ` (Technical Reference §10.3). Name your dry-pyrolysis gas anything else
(`:gas`, `:MMA`, `:volatiles`). You can inspect the classification:

```julia
Pyrolysis.Materials.gas_component_indices(material)        # all gas indices
Pyrolysis.Materials.dry_pyrolysis_gas_indices(material)    # excludes :H2O / :water*
```

### 3.2.4 Bulk vs. intrinsic density and porosity

Two densities coexist, and getting them right matters for porous-media physics:

- **Bulk density `ρ_j`** is the pure-phase density used by *mixing rules*, *mixture
  density*, and *volume-change* kernels. This is what you pass as `ρ`.
- **Intrinsic (skeletal) density `ρ_{i,j}`** is the density of the cell-wall
  material itself, ignoring pores within that component. It is used **only** by the
  porosity formula

  ```
  φ = clamp(1 − Σ_{j∈solid,liquid} ξ_j / ρ_{i,j}, 0, 1)
  ```

  and, through `φ`, by the gas-pressure closure `P = N/φ` and the Carman–Kozeny
  permeability model (Technical Reference §5.3, §9.4).

If you never supply an intrinsic density, it **defaults to the bulk density**, which
makes the pure phase pore-free and gives `φ = 0` for a fully solid cell — the
historical, backward-compatible default. To model a genuinely porous solid you must
supply *either* `ρ_intrinsic` *or* `φ` (not both — passing both raises an
`ArgumentError`):

```julia
# Option A: give the skeletal density directly
porous_char = SolidComponent(:char, ρ = 150.0, c = 1000.0, k = 0.10,
                             ρ_intrinsic = 1450.0)   # cell-wall density of carbon

# Option B: give the pure-phase porosity (only valid with a constant bulk ρ)
porous = SolidComponent(:porous, ρ = 400.0, c = 1000.0, k = 0.10, φ = 0.5)
# → intrinsic_density = ρ / (1 − φ) = 400 / 0.5 = 800 kg/m³
```

The `φ` shortcut requires a **constant** bulk `ρ`; for a temperature-dependent bulk
density you must pass `ρ_intrinsic` explicitly. `φ` must satisfy `0 ≤ φ < 1`.

> **Practical implication.** If you intend to use Carman–Kozeny permeability or the
> gas-pressure formula and you leave intrinsic density at its default, porosity will
> be ~0 and those closures will behave degenerately. Always set `ρ_intrinsic`/`φ`
> on the solids when porosity matters.

### 3.2.5 Inspecting a component

Phase predicates and accessors dispatch statically on the component type:

```julia
is_solid(virgin)        # true
is_gas(gas)             # true
is_liquid(moisture)     # true
phase(virgin)           # SOLID  (a Phase enum value)
molar_mass(gas)         # 0.044  (0.0 for solids/liquids)

# Density / property accessors (qualify through Materials):
Pyrolysis.Materials.get_density(virgin, 300.0)            # bulk ρ at T
Pyrolysis.Materials.get_intrinsic_density(porous, 300.0)  # skeletal ρ
Pyrolysis.Materials.get_permeability(char, 500.0)         # κ at T (0 for gases)
```

---

## 3.3 Property functions

Every per-component thermophysical property is an `AbstractPropertyFunction`: a
callable `p(T)` (and, for some, `p(T, ξ)`). You rarely construct these directly —
the component constructors call `to_property` on whatever you pass — but knowing the
types lets you express any temperature dependence you need.

### 3.3.1 Shorthand conversions (`to_property`)

`to_property` (internal: `Pyrolysis.Materials.to_property`) maps convenient inputs to
property objects, and the component/reaction constructors apply it automatically to
`ρ`, `c`, `k`, `λ`, `κ`, `ρ_intrinsic`, and `h`:

| You pass | You get | Meaning |
|---|---|---|
| a `Real` number, e.g. `0.21` | `ConstantProperty(0.21)` | constant in `T` |
| a 2-tuple `(a, b)` | `LinearProperty(a, b)` | `a + b·T` |
| a 1-arg function `T -> …` | `CallableProperty` | arbitrary `f(T)` |
| an `AbstractPropertyFunction` | passes through unchanged | explicit object |

So `k = 0.21`, `k = (0.20, 1e-4)`, `k = T -> 0.20 + 1e-4*T`, and
`k = LinearProperty(0.20, 1e-4)` are four equivalent ways to specify a conductivity.

### 3.3.2 The property-function types

All are exported by `using Pyrolysis`.

**`ConstantProperty(value)`** — `p(T) = value`.

```julia
k = ConstantProperty(0.21)
```

**`LinearProperty(a, b)`** — `p(T) = a + b·T`. Equivalent to the `(a, b)` shorthand.

```julia
cp = LinearProperty(1420.0, 2.5)   # cp(T) = 1420 + 2.5·T
```

**`PolynomialProperty(coeffs::NTuple)`** — `p(T) = c₀ + c₁T + c₂T² + …`, evaluated
with `evalpoly` (numerically stable, Horner form). Pass the coefficients
low-order-first as a tuple:

```julia
cp = PolynomialProperty((1000.0, 2.0, 0.001))   # 1000 + 2T + 0.001T²
cp(300.0)   # = 1690.0
```

**`ArrheniusProperty(A, E[, T_ref])`** — `p(T) = A · exp(−E / (R_g·T))`, `E` in
J·mol⁻¹. `T_ref` defaults to 298.15 K. Mostly used for rate-constant-shaped
properties; reaction kinetics themselves are specified on the `Reaction` (§3.4), not
as an `ArrheniusProperty`.

```julia
p = ArrheniusProperty(8.5e12, 188e3)   # T_ref defaults to T_REF = 298.15
```

**`TableProperty(temperatures, values)`** — tabulated `(Tᵢ, vᵢ)` with linear
interpolation. Both arguments are tuples; **temperatures must be strictly
increasing** (the constructor errors otherwise). Queries outside the table
**extrapolate flat** (clamp to the nearest endpoint value). A binary search is used
for tables with more than 8 entries.

```julia
k = TableProperty((300.0, 400.0, 500.0), (0.20, 0.22, 0.25))
k(350.0)   # = 0.21  (linear interpolation between 0.20 and 0.22)
k(250.0)   # = 0.20  (flat extrapolation below the table)
```

**`CallableProperty(func)`** — wraps any `f(T) -> value`. Use for closed-form
dependences not covered by the other types:

```julia
k = CallableProperty(T -> 0.10 + 3.0e-4 * (T - 300.0))
```

**`StateDependentProperty(func)`** — a property that depends on both temperature and
the full cell concentration vector, `f(T, ξ) -> value`. This is the mechanism for
**composition-dependent heats of reaction**, most importantly the moisture heat of
sorption, where the sorption enthalpy rises as moisture content falls.

```julia
# Heat of sorption as a function of moisture content (%); ξ is the full
# mass-concentration tuple, so you index the moisture component yourself.
MOISTURE_IDX = 3
heat_sorption = StateDependentProperty((T, ξ) -> begin
    ξ_solid = ξ[1] + ξ[2]                              # dry-solid concentration
    MC = 100.0 * ξ[MOISTURE_IDX] / max(ξ_solid, 1.0)   # moisture content [%]
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3          # [J/kg], endothermic (h > 0)
end)
```

A `StateDependentProperty` **must** be called with state: `p(T, ξ)`. Calling it with
only `T` raises an error. When used as a reaction heat, the kinetics evaluation
passes `ξ` automatically. See §3.4.5 and Technical Reference §6.6.

### 3.3.3 Derivatives (for the implicit Jacobian)

You generally do not call these, but it is useful to know they exist and how they
behave. Temperature derivatives `dp/dT` are exact and closed-form for
`ConstantProperty` (0), `LinearProperty` (`b`), `PolynomialProperty` (unrolled), and
`ArrheniusProperty` (`p·E/(R_g T²)`); they are finite-difference (`δT = 0.1 K`) for
`TableProperty`, `CallableProperty`, and `StateDependentProperty`. Concentration
derivatives `∂p/∂ξ_j` are zero for temperature-only properties and adaptive
finite-difference (`δξ = max(|ξ|·1e-4, 1e-6)`) for `StateDependentProperty`. Because
table/callable/state-dependent properties use finite differences, prefer the
analytic types (constant/linear/polynomial/Arrhenius) where you can — they keep the
Jacobian exact and the solve crisp. See Technical Reference §4.5.

---

## 3.4 Defining reactions

Reactions convert one reactant component into one or two product components by
Arrhenius kinetics. They are usually attached to a material when you build it. The
full physics — smooth temperature gates, the depletion limiter, species sources — is
in **Chapter 4** and Technical Reference §6; here we cover the constructor surface
enough to build a material.

### 3.4.1 Constructor forms

There are four `Reaction` constructors. They differ only in how reactant/product are
named (integer index vs. symbol) and whether there are one or two products.

```julia
# (a) integer index, single product
Reaction(name, reactant::Int => product::Int;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0,
         validate_mass = true)

# (b) symbol, single product (resolved by Material)
Reaction(name, :reactant => :product;
         A, E, h, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0,
         validate_mass = true)

# (c) integer index, two products
Reaction(name, reactant::Int => (p1::Int, p2::Int);
         A, E, h, yields, n = 1.0, T_min = 0.0, T_max = Inf, gate_width = 1.0,
         validate_mass = true)

# (d) symbol, two products (resolved by Material)
Reaction(name, :reactant => (:p1, :p2);
         A, E, h, yields = (y1, y2), n = 1.0, T_min = 0.0, T_max = Inf,
         gate_width = 1.0, validate_mass = true)
```

| Argument | Symbol | Meaning | Units | Default |
|---|---|---|---|---|
| `name` | — | reaction identifier (`Symbol`, positional) | — | required |
| `A` | `A_i` | Arrhenius pre-exponential | s⁻¹ | required |
| `E` | `E_i` | activation energy | J·mol⁻¹ | required |
| `h` | `h` | heat of reaction (number or property fn) | J·kg⁻¹ | required |
| `yields` | `ν_{i,j}` | per-product mass yields (two-product forms) | – | required for (c)/(d) |
| `n` | `n_{i,j}` | reaction order in the reactant (must be ≥ 0) | – | `1.0` |
| `T_min` | `T_min,i` | lower temperature gate | K | `0.0` |
| `T_max` | `T_max,i` | upper temperature gate | K | `Inf` |
| `gate_width` | `w_gate` | tanh ramp width of both temperature gates (must be > 0 and finite) | K | `1.0` |
| `validate_mass` | — | enforce `Σ yields = 1` at construction | – | `true` |

Yields must each be ≥ 0 and the order must satisfy `n ≥ 0` — both are enforced
unconditionally at construction; `validate_mass = false` disables only the
`Σ yields = 1` sum check.

The temperature gates are smooth tanh ramps,
`0.5·(1 + tanh((T − T_min)/gate_width))` (mirrored for `T_max`), so the rate
and its Jacobian stay differentiable across the gate. The default 1 K width
suits kinetics that turn on over tens of kelvin; for gates that must act
nearly as hard switches (melting fronts, phase-change surrogates) use a
narrower width (e.g. `gate_width = 0.1`) and expect a stiffer Jacobian near
the gate. On a driving ΔT of only a few kelvin, a 1 K ramp visibly shifts
the effective switching temperature — narrow the gate accordingly.

### 3.4.2 Symbolic vs. index references

The **symbol** forms ((b), (d)) are the recommended style: you name components, and
the `Material` constructor resolves the names to indices using the `components`
tuple, so you never track integer positions by hand. Internally a symbolic reaction
is stored as a `PendingReaction` until `Material(...)` resolves it.

```julia
pyrolysis = Reaction(:pyrolysis, :virgin => (:char, :gas),
                     A = 1e10, E = 140e3, h = 500e3,   # h > 0: endothermic
                     yields = (0.25, 0.75),            # 25% char + 75% gas
                     T_min = 450.0)                    # gate below 450 K
```

The **index** forms ((a), (c)) reference components by their 1-based position in the
material's `components` tuple. They are handy in programmatically generated schemes:

```julia
rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3)
# component 1 → component 2, single product, default order/gates
```

### 3.4.3 Heat-of-reaction sign convention

Pyrolysis.jl uses the **storage convention `h > 0 = endothermic`** (cools the
material). The volumetric source is `Q_rxn = −Σ h_r r_r`, so an endothermic reaction
produces `Q_rxn < 0`. This is the *opposite* of the ThermaKin/Gpyro publication
convention (`h > 0 = exothermic`); when you transcribe literature parameters, flip
the sign. Examples:

```julia
endo = Reaction(:pyrolysis, :virgin => :gas, A = 1e10, E = 140e3, h =  500e3)  # cools
exo  = Reaction(:char_ox,   :char   => :ash, A = 1e6,  E = 120e3, h = -300e3)  # heats
```

See the overload note H1 in the nomenclature and Technical Reference §6.6.

### 3.4.4 Yields and mass balance

For single-product reactions the product yield is fixed at `1.0` (one unit of
reactant → one unit of product). For two-product reactions you supply
`yields = (y1, y2)`, and **mass conservation requires `y1 + y2 = 1.0`**. With the
default `validate_mass = true`, the constructor enforces this and errors on a
violation; you can bypass the check with `validate_mass = false` (not recommended).

```julia
# Charring: 1 kg virgin → 0.20 kg char + 0.80 kg gas
Reaction(:charring, :virgin => (:char, :gas),
         A = 1e15, E = 200e3, h = 500e3, yields = (0.20, 0.80))
```

Only **single-reactant** reactions are supported (the reactant is consumed with mass
stoichiometry 1.0). Multi-step schemes are modelled as *sequences* of single-reactant
reactions sharing intermediate components — see Chapter 4.

### 3.4.5 State-dependent heat (moisture sorption)

Pass a `StateDependentProperty` (§3.3.2) as `h` to make the heat depend on local
composition. This is the standard way to model the moisture heat of sorption:

```julia
heat_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[3] / max(ξ[1] + ξ[2], 1.0)
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3
end)

drying = Reaction(:drying, :moisture => :H2O,
                  A = 5.6e8, E = 43e3, h = heat_sorption)
```

### 3.4.6 Verifying reactions

Two internal helpers (not exported; qualify through `Materials`) check mass balance
after the fact:

```julia
Pyrolysis.Materials.verify_stoichiometry(rxn)                  # Bool, warns on failure
Pyrolysis.Materials.check_material_mass_balance(mat.reactions) # all reactions
```

Both compare the reactant mass (1.0) against the sum of product yields within a
tolerance of `1e-10` and emit a warning (returning `false`) on imbalance. They are
also invoked automatically by `validate_material` (§3.7) and by `solve`'s
problem-level validation, so for routine use you can rely on those rather than
calling these private helpers directly.

---

## 3.5 Mixing rules: from component to effective properties

The solver integrates **effective** (homogenized) bulk properties, not individual
component properties. Three properties are mixed per material — thermal conductivity
`k`, gas-transfer coefficient `λ`, and permeability `κ` — and **each may use a
different rule**. Heat capacity, mixture density, emissivity, and absorption are not
configurable: they always use fixed weighting rules (concentration- or
volume-weighted; see Technical Reference §5).

### 3.5.1 `MixingRule` and `MixingSpec`

The available rules are the `MixingRule` enum values (all exported):

| Rule | Formula | Valid channels |
|---|---|---|
| `PARALLEL` | `k = Σ_j v_j k_j` (volume-weighted, upper bound) | k, λ, κ |
| `SERIES` | `1/k = Σ_j v_j/k_j` (harmonic, lower bound) | k, λ, κ |
| `WEIGHTED` | `k = β·k_∥ + (1−β)·k_series` | k, λ, κ |
| `BRUGGEMAN` | EMT: `Σ_j v_j (k_j−k)/(k_j+2k) = 0` | **k only** |
| `CARMAN_KOZENY` | `κ = min(ℓ²φ³ / (180(1−φ)²), ℓ²)` (porosity-based) | **κ only** |

A `MixingSpec` bundles a rule with its parameters:

```julia
MixingSpec(rule::MixingRule, β::Real = 0.5; length_scale::Real = NaN)
```

- `β` is the blend weight, used **only** by `WEIGHTED`, and must lie in `[0, 1]`
  (default `0.5`). It is harmlessly carried for other rules.
- `length_scale` (`ℓ`, the characteristic particle/pore diameter in metres) is used
  **only** by `CARMAN_KOZENY` and is **required** for it — the constructor errors if
  you request `CARMAN_KOZENY` without a positive, finite `length_scale`.

```julia
MixingSpec(WEIGHTED, 0.6)                       # 60% parallel / 40% series for k
MixingSpec(PARALLEL)                            # β ignored
MixingSpec(CARMAN_KOZENY; length_scale = 5e-5)  # ℓ = 50 µm
MixingSpec(BRUGGEMAN)                            # EMT (k channel only)
```

### 3.5.2 `MaterialMixing`

`MaterialMixing` builds the `(k, λ, κ)` NamedTuple that `Material` expects. Each
channel accepts either a `MixingSpec` or a bare `MixingRule` (β defaults to 0.5 for
the bare form). The **defaults** are:

```julia
MaterialMixing(; k = MixingSpec(WEIGHTED, 0.5),   # conductivity
                 λ = MixingSpec(PARALLEL, 0.5),   # gas transfer
                 κ = MixingSpec(WEIGHTED, 0.5))   # permeability
```

```julia
mixing = MaterialMixing(
    k = MixingSpec(WEIGHTED, 0.5),                       # blended bounds
    λ = MixingSpec(PARALLEL),                            # parallel
    κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),  # porosity-based κ
)

# bare-rule shorthand (β = 0.5):
mixing2 = MaterialMixing(k = WEIGHTED, λ = PARALLEL, κ = SERIES)
```

If you omit `mixing` from `Material(...)`, the defaults above are used.

### 3.5.3 Choosing rules — guidance

- **Thermal conductivity `k`.** The default `WEIGHTED` with `β = 0.5` (midway between
  the parallel upper bound and series lower bound) is a sensible neutral choice. Use
  `BRUGGEMAN` for char-forming materials where the microstructure is random/
  self-similar — it captures percolation-like drops as a continuous phase breaks up
  (Technical Reference §5.5). `BRUGGEMAN` is solved by a warm-started Newton
  iteration and is valid **only for `k`**; requesting it for `λ` or `κ` raises an
  error at evaluation time.
- **Gas transfer `λ`.** `PARALLEL` (the default) is appropriate; EMT (`BRUGGEMAN`)
  is rejected for diffusivities because there is no physical EMT closure for them.
- **Permeability `κ`.** Two mutually exclusive approaches:
  1. **Component-based** (`PARALLEL`/`SERIES`/`WEIGHTED`): mix the per-component `κ`
     values you set on the solids/liquids.
  2. **Carman–Kozeny** (`CARMAN_KOZENY`): ignore component `κ` and compute
     permeability from **porosity** via `κ = min(ℓ²φ³/(180(1−φ)²), ℓ²)` — the
     packed-bed law capped at the free-flow scale `ℓ²` for very open structures
     (`φ ≳ 0.93`). This captures the steep rise in permeability as virgin → char
     raises porosity.

> **Important caveat for Carman–Kozeny.** Whether the solver runs in Darcy–Fick
> (pressure-driven) mode is decided by `has_darcy_flow(material)`, which checks
> whether any *component* carries `κ > 0` — it does **not** look at the κ mixing
> rule. So choosing `CARMAN_KOZENY` alone does **not** turn on Darcy transport: with
> all component `κ = 0`, `has_darcy_flow` returns `false` and gas transport stays
> Fickian. If you want pressure-driven flow *and* a Carman–Kozeny effective
> permeability, set a non-zero component `κ` (to trip `has_darcy_flow`) while using
> the `CARMAN_KOZENY` rule for the effective value, or use component-based κ mixing.
> Also remember that Carman–Kozeny needs a real porosity, so set `ρ_intrinsic`/`φ`
> on the solids (§3.2.4) — otherwise `φ ≈ 0` and `κ ≈ 0`.

### 3.5.4 Inspecting effective properties

The effective-property functions live in the `Properties` submodule (qualify them).
They take `(ξ, T, material)` where `ξ` is the mass-concentration tuple:

```julia
using Pyrolysis
P = Pyrolysis

ξ = (300.0, 50.0, 0.5)   # (virgin, char, gas) concentrations [kg/m³]
T = 600.0

P.Properties.mixture_density(ξ, T, material)         # ρ with swelling
P.Properties.effective_heat_capacity(ξ, T, material) # ρc_p (matrix only)
P.Properties.effective_conductivity(ξ, T, material)  # k_eff
P.Properties.effective_gas_transfer(ξ, T, material)  # λ_eff
P.Properties.effective_permeability(ξ, T, material)  # κ_eff
P.Properties.effective_emissivity(ξ, T, material)    # ε_eff
P.Properties.effective_absorption(ξ, material)       # α_eff (note: no T arg)
P.Properties.porosity(ξ, T, material)                # φ
```

Two physics notes worth remembering (both in Technical Reference §5):

- **Heat capacity excludes gas storage.** `effective_heat_capacity` sums only
  solid/liquid components: `ρc_p^eff = Σ_{j∈S,L} ξ_j c_{p,j}(T)`. Gas sensible heat
  is carried by the advective source term instead (the quasi-steady-gas
  approximation, ~0.2% of the matrix term at 1 atm). `total_heat_capacity` (which
  includes gas) exists for diagnostics only.
- **Conductivity and emissivity use the swelling-weighted (matrix) volume
  fractions** — γ = 0 components (canonical gases) carry no weight in these
  surface/matrix properties, matching the ThermaKin convention; a fully depleted
  cell falls back to the all-phase occupancy so k_eff degrades to the gas-mixture
  conductivity instead of 0. **Gas transfer uses the all-phase volume fractions**
  (gas included) — diffusion happens through the pore space the matrix basis
  excludes. Component-based permeability mixes on the matrix basis. Absorption uses
  concentration weighting; mixture density uses mass fractions with the swelling
  factor.

---

## 3.6 Assembling a `Material`

### 3.6.1 The constructor

```julia
Material(; name::Symbol,
           components::Tuple,
           reactions::Tuple = (),
           mixing = MaterialMixing(),
           lateral_shrinkage_law = nothing,
           depletion_limiter = DepletionLimiter())
```

| Keyword | Meaning | Default |
|---|---|---|
| `name` | material identifier (`Symbol`) | required |
| `components` | tuple of components | required |
| `reactions` | tuple of reactions (may be empty `()`) | `()` |
| `mixing` | `(k, λ, κ)` NamedTuple from `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `A(t)/A_0 = law(χ̄)` or `nothing` (identity) | `nothing` |
| `depletion_limiter` | `DepletionLimiter` applied to every reaction-rate evaluation (smooth `tanh(ξ/threshold)` roll-off near reactant depletion; §4.10) | `DepletionLimiter()` (`threshold = 1.0` kg·m⁻³, enabled) |

- **Component order defines the index space.** The first component is index 1, the
  second is 2, etc. The state vector's species blocks follow this order, and integer
  reaction references (§3.4.2) use it. Component names must be unique.
- **Symbolic reactions are resolved here.** Any `PendingReaction` built with symbol
  references is matched against the `components` tuple's name→index map. An unknown
  symbol raises a clear error naming the offending reaction and component.
- **`lateral_shrinkage_law`** is the variable-cross-section (intumescence/shrink)
  hook; it must be a `Function` of the column-average pyrolysis progress `χ̄`, or
  `nothing`. Most simulations leave it `nothing`. See Technical Reference §10.5 and
  Chapter 7.

### 3.6.2 Worked example: charring wood with moisture

A wood with virgin solid, char solid, liquid moisture, a dry pyrolysis gas, and water
vapour. The solids carry permeabilities, so the material will run in Darcy–Fick mode;
it uses `WEIGHTED` conductivity and `CARMAN_KOZENY` permeability, with porosity set
on the solids via `φ`.

```julia
using Pyrolysis

# --- Components -----------------------------------------------------------
virgin = SolidComponent(:virgin,
    ρ = 530.0, c = (1500.0, 1.0), k = 0.12,
    ε = 0.88, λ = 2e-5,
    κ = 1e-15,           # tight virgin matrix (trips Darcy mode)
    φ = 0.55)            # pure-phase porosity → intrinsic density set internally

char = SolidComponent(:char,
    ρ = 150.0, c = 1000.0, k = 0.10,
    ε = 0.95, α = 5.0, λ = 2e-5,
    κ = 1e-12,           # open, porous char
    φ = 0.85)

moisture = LiquidComponent(:moisture, ρ = 1000.0, c = 4180.0, k = 0.60)

gas   = GaseousComponent(:gas, M = 0.044, c = 1800.0, k = 0.03, λ = 1e-4)
vapor = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025, λ = 2e-5)

# --- Reactions ------------------------------------------------------------
pyrolysis = Reaction(:pyrolysis, :virgin => (:char, :gas),
                     A = 1e10, E = 140e3, h = 500e3,     # endothermic
                     yields = (0.25, 0.75), T_min = 450.0)

# Moisture evaporation with state-dependent heat of sorption
heat_sorption = StateDependentProperty((T, ξ) -> begin
    MC = 100.0 * ξ[3] / max(ξ[1] + ξ[2], 1.0)
    (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3
end)
drying = Reaction(:drying, :moisture => :H2O,
                  A = 5.6e8, E = 43e3, h = heat_sorption)

# --- Material -------------------------------------------------------------
wood = Material(
    name = :CharringWood,
    components = (virgin, char, moisture, gas, vapor),   # indices 1..5
    reactions  = (pyrolysis, drying),
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),
        λ = MixingSpec(PARALLEL),
        κ = MixingSpec(CARMAN_KOZENY; length_scale = 5e-5),
    ),
)

# --- Inspect / validate ---------------------------------------------------
# All reached via the Materials submodule (none are top-level Pyrolysis exports).
# n_components / n_reactions / has_darcy_flow / dry_pyrolysis_gas_indices are
# exported from Materials; validate_material is an internal helper — both kinds work
# through the qualified `Pyrolysis.Materials.` path.
Pyrolysis.Materials.n_components(wood)              # 5
Pyrolysis.Materials.n_reactions(wood)               # 2
Pyrolysis.Materials.has_darcy_flow(wood)            # true (virgin/char κ > 0)
Pyrolysis.Materials.dry_pyrolysis_gas_indices(wood) # (4,)  — :H2O (5) excluded
Pyrolysis.Materials.validate_material(wood)         # true (internal helper)
```

### 3.6.3 Worked example: a simple polymer (PMMA)

A non-charring polymer needs only a virgin solid and a single gas, one decomposition
reaction, and no permeability (so transport stays Fickian). Note the positive `h`:
PMMA depolymerization is endothermic (≈ 870 kJ·kg⁻¹ absorbed), so under the
`h > 0 = endothermic` convention it takes a positive value.

```julia
using Pyrolysis

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86),
        GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:depolymerization, :virgin => :MMA,
                 A = 8.5e12, E = 188e3, h = 870e3),   # h > 0 = endothermic
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)

Pyrolysis.Materials.validate_material(pmma)   # true (internal helper, qualified)
```

With `κ = 0` on every component, `has_darcy_flow(pmma)` is `false` and gas leaves by
Fickian diffusion only — the default transport mode. To make MMA absorb in depth you
would set `α` on the solid and choose `radiation_model = BEER_LAMBERT` at solve time.

---

## 3.7 Validating a material

`validate_material` is an **internal** helper (not exported; reach it as
`Pyrolysis.Materials.validate_material`). It is the one call to make if you want an
explicit material check before building a problem. It checks that:

- the material has at least one component,
- every reaction's reactant/product indices are in range `1:NC`, and
- every reaction is mass-balanced (via `verify_stoichiometry`).

```julia
Pyrolysis.Materials.validate_material(wood)   # returns true; warns on problems
```

Mass-balance failures are reported as **warnings** (the function still returns
`true`); index errors **throw**. For a standalone balance check use
`Pyrolysis.Materials.check_material_mass_balance(material.reactions)`. Note that the
two-product `Reaction` constructors already enforce `Σ yields = 1.0` at construction
time (unless `validate_mass = false`), so a material built normally will pass.

Because these helpers are private, the **public** way to validate is to build a
`PyrolysisProblem` and let `solve` validate by default (`validate = true` on
`solve`); the problem-level validation runs the same stoichiometry and index checks.
Calling the internal `validate_material` as you build the material is still useful for
the clearest, earliest error messages during development.

---

## 3.8 Quick reference

**Component defaults** (the values you get if you omit the keyword):

| Field | Solid | Liquid | Gas |
|---|---|---|---|
| `ε` (emissivity) | `0.9` | `0.95` | `0.0` |
| `α` (absorption) | `0.0` | `0.0` | `0.0` |
| `γ` (swelling) | `1.0` | `0.0` | `0.0` |
| `λ` (gas transfer) | `0.0` | `2e-5` | required |
| `κ` (permeability) | `0.0` | `0.0` | n/a (no field) |
| `ρ_intrinsic` | `= ρ` | `= ρ` | n/a |
| density | required `ρ` | required `ρ` | ideal gas from `M` |

**Reaction defaults:** `n = 1.0`, `T_min = 0.0`, `T_max = Inf`,
`validate_mass = true`; single-product yield fixed at `1.0`, two-product `yields`
must sum to `1.0`. Sign: `h > 0` endothermic.

**Mixing defaults:** `k = WEIGHTED (β=0.5)`, `λ = PARALLEL (β=0.5)`,
`κ = WEIGHTED (β=0.5)`. `BRUGGEMAN` is k-only; `CARMAN_KOZENY` is κ-only and requires
`length_scale`.

**Physical constants used by the constructors:** `R_g = 8.31446261815324
J·mol⁻¹·K⁻¹`, `P_ref = 101325 Pa`, `T_ref = 298.15 K`, `σ = 5.670374419e-8
W·m⁻²·K⁻⁴` (all exported as `R_GAS`, `P_REF`, `T_REF`, `STEFAN_BOLTZMANN`).

**Where to go next:** reactions in depth and TGA parameter recovery → Chapter 4;
boundary and initial conditions → Chapter 5; meshes → Chapter 6; building and
running the problem → Chapters 7–8. The physics behind everything here is in
Technical Reference §4 (materials & property functions), §5 (effective properties &
mixing), §6 (kinetics), and §10 (volume change & swelling).
