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
5. **Material** (`Material`) — the assembled object, validated together with the
   rest of the problem by `solve` (see §3.7).

A minimal end-to-end example (one solid, one gas, one reaction) appears in the
package's top-level docstring and in Chapter 2; this chapter expands every option.

> **A note on what is exported.** `using Pyrolysis` brings the *constructors* into
> scope: `SolidComponent`, `LiquidComponent`, `GaseousComponent`, `Material`,
> `Reaction`, `DepletionLimiter`, all property-function types, the `MixingRule` enum
> values, `MixingSpec`, and `MaterialMixing`. The various **helper and inspection
> functions** used in this chapter are *not* part of the curated top-level surface,
> so you reach them by qualifying through the submodule that defines them. This works
> the same way whether or not the function is exported from that submodule — Julia's
> `Module.name` syntax resolves any name defined in the module:
>
> - From the `Materials` submodule: `Pyrolysis.Materials.n_components`,
>   `Pyrolysis.Materials.n_reactions`, `Pyrolysis.Materials.dry_pyrolysis_gas_indices`,
>   are *exported* from `Materials`; the
>   stoichiometry helpers `Pyrolysis.Materials.verify_stoichiometry` and
>   `Pyrolysis.Materials.mass_balance_error` are *internal* (not exported, not part
>   of the stable public API) but are still callable through the same qualified
>   path.
> - From the `Properties` submodule: the effective-property functions
>   (`Pyrolysis.Properties.effective_conductivity`,
>   `mixture_emissivity`, `mixture_absorptivity`, `effective_permeability`,
>   `skeletal_volume_fraction`, …; §3.5.4).
> - From the `Physics` submodule: `reaction_rate` and `one_sided_factor`
>   (Chapter 4).
>
> For convenience you may `using Pyrolysis.Internal`, which re-exports the union of
> every submodule's *exports* and so brings the exported helpers
> (`n_components`, `dry_pyrolysis_gas_indices`, …) into scope
> unqualified — but **not** the internal stoichiometry helpers above, which still
> require the `Pyrolysis.Materials.` prefix. `Internal` is not the stable public
> surface; prefer qualified access in production code. Material validation happens
> at the problem level: build a `PyrolysisProblem` and let `solve` validate it (the
> default, `validate = true`); see §3.7.

---

## 3.2 Defining components

### 3.2.1 `SolidComponent`

```julia
virgin = SolidComponent(:virgin; ρ=530., ρ_intrinsic=1500.,
    c=(1500.,1.), k=.12, ε=.88, λ=2e-5, κ=1e-15)
char = SolidComponent(:char; ρ=150., ρ_intrinsic=1950.,
    c=1000., k=.10, ε=.95, α=5., λ=2e-5, κ=1e-12)
```

Bulk density `ρ`, heat capacity `c` and conductivity `k` are required.
Participating lumped matrices also require explicit superficial transfer `λ`.
Zero `λ` closes wall transfer independently of `κ`; zero permeability closes
Darcy flow. Neither coefficient is inferred from gas loading.

`ε` accepts a temperature law; `absorptivity` sets a constant incident-radiation
absorptivity, with `nothing` applying Kirchhoff's relation. `α` is mass-based
in-depth absorption and is used under `BEER_LAMBERT`.

`γ` sets participation in the volume-change law and property weights. A
zero-`γ` condensed phase still occupies skeletal volume and must not carry
nonzero `λ` or `κ`. Specify intrinsic density or pure-phase porosity to
describe that skeletal volume. Full signatures and units are in §16.2.1.

---

### 3.2.2 `LiquidComponent`

Liquids have the condensed component fields. Their defaults are `γ=0`,
`ε=.95`, `κ=0` and `λ=nothing`.

```julia
moisture = LiquidComponent(:moisture; ρ=1000., c=4180., k=.60)
melt = LiquidComponent(:melt; ρ=900., φ=.02, c=2000., k=.2,
                       γ=1., λ=1e-6, κ=0.)
```

Bound moisture with `γ=0` contributes no direct strain or matrix transport
weight, but it occupies skeletal volume and stores heat. A matrix-forming
melt with `γ>0` supplies explicit lumped transfer. The initialized sum of
all condensed skeletal volumes must leave resolved positive gas pore volume.

---

### 3.2.3 `GaseousComponent`

```julia
air = GaseousComponent(:air; M=0.02897, c=1006., k=0.026)
product = GaseousComponent(:product; M=0.044, c=1200., k=0.03)
```

`M` is a positive molar mass [kg/mol]. Heat capacity `c`, conductivity `k`,
emissivity `ε`, absorption `α` and swelling `γ` have their usual property-law
meaning. Gases carry no matrix `λ` or `κ`. `LumpedTransfer` uses coefficients
on the condensed matrix; `DustyGas` derives gas-specific Knudsen transport
from molar mass and pore geometry. Gas loading does not enter matrix transport
weights, regardless of `γ`. Bulk concentration is stored per total cell volume;
the central pore EOS divides by skeletal porosity to obtain intrinsic density.

---

### 3.2.4 Bulk and intrinsic density, pore volume

The condensed bulk density `ρ` defines volume-change and mixing data. The
intrinsic density `ρ_intrinsic` defines skeletal occupancy. For constant bulk
density, `φ` is a convenient alternative: `ρ_intrinsic=ρ/(1-φ)`.

```julia
matrix = SolidComponent(:matrix; ρ=500., ρ_intrinsic=1500.,
    c=1500., k=.2, λ=1e-5, κ=1e-14)
```

For the actual initialized composition, pore volume is
`φ = 1 - sum(ξ[j]/ρ_intrinsic[j](T) for condensed j)`. It must be resolved and
positive in every gas-bearing cell, even for sealed gas and zero permeability.
A missing intrinsic density defaults to bulk density; a partially filled cell
can therefore still be admissible. A name search cannot determine admissibility.
Zero `γ` excludes a component from matrix transport but not skeletal occupancy.
For dense melts, an explicit nominal pore volume requires joint porosity,
mesh and tolerance qualification (§7.5); the numerical divisor floor is not a
physical porosity choice.

---

### 3.2.5 Inspecting a component

Phase predicates and accessors dispatch statically on the component type:

```julia
is_solid(virgin)        # true
is_gas(gas)             # true
is_liquid(moisture)     # true
phase(virgin)           # SOLID  (a Phase enum value)
molar_mass(gas)         # 0.044  (0.0 for solids/liquids)

# Density / property accessors (qualify through Materials):
virgin.density(300.0)                                     # bulk ρ at T (property field)
Pyrolysis.Materials.get_intrinsic_density(porous, 300.0)  # skeletal ρ
Pyrolysis.Materials.get_permeability(char, 500.0)         # κ at T (0 for gases)
Pyrolysis.Materials.get_emissivity(char, 800.0)           # ε_j(T)
Pyrolysis.Materials.get_absorptivity(char, 800.0)         # a_j: absorptivity, else ε_j(T)
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

**`SutherlandProperty(μ_ref, T_ref, S)`** — Sutherland's law
`μ_ref (T/T_ref)^{3/2} (T_ref + S)/(T + S)`; `SutherlandProperty()` is air
(1.716e-5 Pa·s at 273.15 K, `S = 110.4` K, so `μ(300 K) ≈ 1.846e-5` Pa·s), the
default pore-gas viscosity of a `Material` (§3.6.1).

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

You never supply derivatives. The Jacobian is colored sparse forward-mode
automatic differentiation of the residual (Chapter 11), so every property law is
evaluated with dual-number arguments and its exact derivative comes out of the
same evaluation: the slope of the piecewise-linear interpolant for a
`TableProperty`, the derivative of your function for a `CallableProperty`, and
both `dp/dT` and `∂p/∂ξ_j` for a `StateDependentProperty`. A callable must
therefore be generic Julia code — do not restrict its arguments to `Float64`.
Every law evaluates on the property bracket
`[200, 3000]` K (its argument is clamped to it, so a nonphysical Newton trial
temperature cannot reach a fractional power), and the derivatives are zero
outside the bracket. The kernels hand the laws the cell's trial temperature, the
cell temperature smoothly clamped into the bracket with 1 K transitions, which is
the cell temperature itself on every physical state. See Technical Reference
§4.4–§4.5.

**Heat capacity and the sensible enthalpy.** The energy rows also use each
component's sensible enthalpy `h_j(T) = ∫_{T_ref}^{T} c_{p,j} dT′`
(`T_ref = 298.15` K): in the enthalpy carried by the gas flux and, with
`energy_state = :enthalpy`, in the energy state itself. Constant, linear,
polynomial and tabulated heat capacities are integrated in closed form; any other
law (a function, an Arrhenius or Sutherland law) is integrated once when the
problem is set up, on a 1 K grid over the bracket, and interpolated by cubic
Hermite polynomials, so its derivative is exact to that interpolant.

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
produces `Q_rxn < 0`. This is the *opposite* of the ThermaKin publication
convention (`h > 0 = exothermic`); when you transcribe ThermaKin parameters, flip
the sign and normalize the reaction extent (Technical Reference §6.9.1).
Gpyro and FDS use the convention of this package but still require basis and
unit checks. Examples:

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

The `Reaction` constructor already enforces mass balance (the product yields must
sum to the reactant mass 1.0) and throws on an imbalance unless you pass
`validate_mass = false`. To re-check a reaction after the fact, use the internal
helper (not exported; qualify through `Materials`):

```julia
Pyrolysis.Materials.verify_stoichiometry(rxn)   # Bool, warns on failure
```

It compares the reactant mass (1.0) against the sum of product yields within a
tolerance of `1e-10` and emits a warning (returning `false`) on imbalance.
`validate_problem(problem; strict = true)` runs it on every reaction of the
material (§3.7).

---

## 3.5 Mixing rules: from component to effective properties

The solver integrates **effective** (homogenized) bulk properties, not individual
component properties. Three properties are mixed per material — thermal conductivity
`k`, gas-transfer coefficient `λ`, and permeability `κ` — and **each may use a
different rule**. Heat capacity, emissivity, and absorption are not
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

For thermal conductivity, `v_j` is the normalized thermal volume weight.
For `λ` and component-based `κ`, replace it with normalized condensed matrix
occupancy `a_j ∝ γ_j ξ_j/ρ_bulk,j`; gases do not contribute. Mixing within
a cell is separate from matching total species flux through adjacent half cells
(Technical Reference §9).

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
  (Technical Reference §5.5). `BRUGGEMAN` is solved by a Newton
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

Carman–Kozeny directly supplies the effective permeability, including when
all component permeabilities are zero. Both gas models use this permeability.
The material still needs resolved physical porosity and a qualified matrix
endpoint or termination domain (§7.5).

### 3.5.4 Inspecting effective properties

The effective-property functions live in the `Properties` submodule (qualify them,
or `using Pyrolysis.Internal`). They take `(ξ, T, material)` where `ξ` is the
mass-concentration tuple; `material` below is any three-component material
(virgin, char, gas), such as the one of §3.6.2 without its moisture:

```julia
using Pyrolysis
P = Pyrolysis.Properties

ξ = (300.0, 50.0, 0.5)   # (virgin, char, gas) concentrations [kg/m³]
T = 600.0

sum(ξ)                                   # total bulk mass density
P.effective_conductivity(ξ, T, material)   # k_eff
P.effective_gas_transfer(ξ, T, material)   # λ_eff
φ = 1 - P.skeletal_volume_fraction(ξ, T, material)   # pore fraction (porosity)
P.effective_permeability(ξ, T, material, φ)          # κ_eff (CARMAN_KOZENY reads φ)
P.mixture_emissivity(ξ, T, material)       # ε_eff
P.mixture_absorptivity(ξ, T, material)     # a_eff (= ε_eff by Kirchhoff unless set)
P.effective_absorption(ξ, material)        # α_eff [1/m] (note: no T argument)

# Volumetric heat capacity ρc_p = Σ_j ξ_j c_p,j(T): every component, gas included
sum(ξ[j] * material.components[j].heat_capacity(T) for j in eachindex(ξ))
```

What the solver evaluates is the same functions, called by its state mapping at
each cell's trial temperature and at the positive part of each concentration
(a smooth one-sided guard that is exactly the concentration on physical states),
with smooth floors where a trial state could make them vanish: the conductivity
at 1e-6 W·m⁻¹·K⁻¹, the pore fraction at 1e-12, the gas viscosity at `MU_FLOOR`
(§3.6.1), and the cell's heat capacity `Σ_j m_j c_{p,j}(T)` at `RHO_CP_FLOOR`
(1e3 J·m⁻³·K⁻¹) times the cell's initial volume.

Two physics notes worth remembering (both in Technical Reference §5):

- **Heat capacity includes the pore gas.** The cell's heat capacity sums every
  component: `ρc_p^eff = Σ_j ξ_j c_{p,j}(T)`. The gas is at the matrix
  temperature and stores sensible heat like every other phase; the enthalpy its
  flux carries enters as the gas enthalpy flux `Σ_j J_j h_j` through each face
  (Technical Reference §7).
- **Conductivity and emissivity use the swelling-weighted (matrix) volume
  fractions** — γ = 0 components (canonical gases) carry no weight in these
  surface/matrix properties. Gas exclusion agrees with ThermaKin at zero gas
  swelling, but its liquid and mixture swelling rules differ (Technical Reference
  §5.14). As a cell's
  matrix volume fraction drops through `1e-3 … 1e-4` the weights blend smoothly
  into the all-phase occupancy, so an emptied cell degrades to the gas-mixture
  conductivity instead of 0. Transport uses a separate, strictly condensed
  γ-weighted bulk-volume basis. Both superficial transfer and component-based
  permeability use those matrix weights, with an explicit endpoint and domain
  when that matrix disappears. Skeletal pore occupancy includes every condensed
  phase regardless of γ. Absorption retains concentration weighting.


---

## 3.6 Assembling a `Material`

### 3.6.1 The constructor

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
| `name` | material identifier (`Symbol`) | required |
| `components` | tuple of components | required |
| `reactions` | tuple of reactions (may be empty `()`) | `()` |
| `mixing` | `(k, λ, κ)` NamedTuple from `MaterialMixing` | `MaterialMixing()` |
| `lateral_shrinkage_law` | `A(t)/A_0 = law(χ̄)` or `nothing` (identity) | `nothing` |
| `depletion_limiter` | `DepletionLimiter` applied to every reaction-rate evaluation (a C¹ rate factor that is zero at and below zero concentration and joins `tanh(ξ/threshold)` above it; §4.10) | `DepletionLimiter()` (`threshold = 1.0` kg·m⁻³, `join = 1.0`, `base_scale = 1e-8` kg·m⁻³, enabled) |
| `gas_transport` | declared pore-gas model and endpoint domain | `LumpedTransfer()` for gas-bearing materials |
| `gas_viscosity` | pore-gas dynamic viscosity μ(T) [Pa·s]; number, `(p0, p1)` tuple, `Function`, or property function | `SutherlandProperty()` (air) |

- **Component order defines the index space.** The first component is index 1, the
  second is 2, etc. The per-cell component masses of the state follow this order,
  and integer reaction references (§3.4.2) use it. Component names must be unique.
- **Symbolic reactions are resolved here.** Any `PendingReaction` built with symbol
  references is matched against the `components` tuple's name→index map. An unknown
  symbol raises a clear error naming the offending reaction and component.
- **`lateral_shrinkage_law`** is the variable-cross-section (intumescence/shrink)
  hook; it must be a `Function` of the column-average pyrolysis progress `χ̄`, or
  `nothing`. Most simulations leave it `nothing`. A material that carries a law
  must be solved on a moving mesh (`use_ale = true`) in planar coordinates; `setup`
  refuses it otherwise. See Technical Reference §10.5 and Chapters 6 and 7.
- **`gas_viscosity`** matters only in lumped gas transport mode, where it enters as the
  mobility `κ/μ` — so getting μ wrong by a factor is exactly the same as getting
  permeability wrong by the reciprocal factor. The default is air, which is fine
  for an N₂-swept pore space but overstates the viscosity of wood pyrolyzates by
  roughly 1.7× at 300 K falling to 1.3× at 1000 K. Supply the blend correlation
  for the gas your material actually produces:

  ```julia
  # Cantera mixture-averaged fit, ascending coefficients, valid 280–1200 K
  Material(name = :WhitePine, components = (...), reactions = (...),
           gas_viscosity = PolynomialProperty(mu_coeffs))
  ```

  It is **one curve per material**, so it averages over the pore gas shifting
  from background air to volatiles as pyrolysis proceeds; it cannot resolve that
  transition. The constructor rejects a curve that is non-positive or non-finite
  at any of 200, 300, 500, 800, 1200 and 1800 K, and the solver floors μ smoothly
  at `Pyrolysis.Properties.MU_FLOOR` (1e-12 Pa·s) as the backstop for trial
  states — both because a μ that reaches zero reverses the gas velocity rather
  than merely degrading it. Like every property law it is evaluated at the
  temperature clamped to `[200, 3000]` K. See Technical Reference §5.9.

### 3.6.2 Worked example: charring wood with moisture

A wood with virgin solid, char solid, liquid moisture, a dry pyrolysis gas, and water
vapour. The solids carry permeabilities, so the material will run in lumped gas transport mode;
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

gas   = GaseousComponent(:gas, M = 0.044, c = 1800.0, k = 0.03,)
vapor = GaseousComponent(:H2O, M = 0.018, c = 2000.0, k = 0.025,)

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

# --- Inspect ---------------------------------------------------------------
# All reached via the Materials submodule (none are top-level Pyrolysis exports).
Pyrolysis.Materials.n_components(wood)              # 5
Pyrolysis.Materials.n_reactions(wood)               # 2

Pyrolysis.Materials.dry_pyrolysis_gas_indices(wood) # (4,)  — :H2O (5) excluded
```

### 3.6.3 A polymer with specified pore volume

A gas-bearing polymer needs a physical pore volume even when the Darcy
permeability is zero. The PMMA cone examples declare nominal porosity and
lumped transfer separately from thermal and kinetic properties. Set
`PYROLYSIS_PMMA_POROSITY` to vary their pore-volume input. The transport
coefficient is superficial; changing porosity changes storage and the
pressure potential. Qualify this input with mesh and pressure-tolerance
refinement, rather than treating a trial divisor floor as a material property.

A positive permeability adds Darcy flow to either gas model. `κ=0` retains
a legal wall-only material; `λ=0, κ>0` is a legal Darcy-only lumped material.
A closed material can set both to zero and still has a pore EOS.

---

## 3.7 Validating a material

A material is checked in two places: its components and reactions validate
themselves when they are constructed, and the assembled material is validated as
part of the problem.

Each two-product `Reaction` constructor enforces `Σ yields = 1.0` at construction
time (unless `validate_mass = false`), and a symbolic component reference that
names no component of the material throws when the `Material` is built.

Property laws are checked when the **component** is built: every temperature-
dependent `ρ`, `c`, `k`, `λ`, `κ` and `ε` is sampled on the property bracket
`[200, 3000]` K at 1 K steps, and the constructor throws if a density, heat
capacity or conductivity is not positive there (a gas may have `c = 0`), if
`λ` or `κ` is negative, or if `ε` leaves `[0, 1]`. The message names the first
offending temperature — a linear fit `k = 0.3 − 4.2e-4·T` is rejected at
715 K, where it would make a cell a perfect insulator. Hold a fit constant
outside its measured range (`T -> law(clamp(T, T_lo, T_hi))`), or pass
`validate_properties = false` when the law comes from an optimizer's search
parameters and an unphysical value should reach the objective instead of
throwing.

Component and material construction validate transport inputs, including explicit
lumped matrix coefficients, named gas species and endpoint geometry. The
problem check validates mesh geometry, boundary tags and complete reservoir
states, initial-condition dimensions, the time span and initial step. Every
gas-bearing initial state must leave positive skeletal pore volume, including
sealed materials. Setup then checks the pore margin against the actual solver
tolerances, reachable matrix support and endpoint assumptions. Proposed and
post-event states undergo the same physical-domain checks. An explicitly zero
transport channel is valid; an omitted required coefficient is an error.

`validate_problem(problem; strict=true)` also rechecks reaction mass balance.
Errors are collected into one diagnostic when possible.

---

## 3.8 Quick reference

**Component defaults** (the values you get if you omit the keyword):

| Field | Solid | Liquid | Gas |
|---|---|---|---|
| `ε` (emissivity) | `0.9` | `0.95` | `0.0` |
| `absorptivity` (surface absorptivity) | `nothing` (Kirchhoff) | `nothing` (Kirchhoff) | n/a |
| `α` (absorption) | `0.0` | `0.0` | `0.0` |
| `γ` (swelling) | `1.0` | `0.0` | `0.0` |
| `λ` (superficial matrix transfer) | explicit under lumped transport | explicit for participating liquid | not a gas input |
| `κ` (permeability) | `0.0` | `0.0` | n/a (no field) |
| `ρ_intrinsic` | `= ρ` | `= ρ` | n/a |
| density | required `ρ` | required `ρ` | ideal gas from `M` |

**Reaction defaults:** `n = 1.0`, `T_min = 0.0`, `T_max = Inf`,
`gate_width = 1.0` K, `validate_mass = true`; single-product yield fixed at `1.0`, two-product `yields`
must sum to `1.0`. Sign: `h > 0` endothermic.

**Mixing defaults:** `k = WEIGHTED (β=0.5)`, `λ = PARALLEL (β=0.5)`,
`κ = WEIGHTED (β=0.5)`. `BRUGGEMAN` is k-only; `CARMAN_KOZENY` is κ-only and requires
`length_scale`.

**Physical constants used by the constructors:** `R_g = 8.31446261815324
J·mol⁻¹·K⁻¹`, `P_ref = 101325 Pa`, `T_ref = 298.15 K`, `σ = 5.670374419e-8
W·m⁻²·K⁻⁴`, exported by the `Materials` submodule as `R_GAS`, `P_REF`, `T_REF`,
`STEFAN_BOLTZMANN` (with the property bracket `T_PROPERTY_MIN = 200` K,
`T_PROPERTY_MAX = 3000` K): `Pyrolysis.Materials.R_GAS`, or unqualified after
`using Pyrolysis.Internal`.

**Limiter defaults:** `DepletionLimiter(threshold = 1.0, join = 1.0,
base_scale = 1e-8, enabled = true)`.

**Where to go next:** reactions in depth and TGA parameter recovery → Chapter 4;
boundary and initial conditions → Chapter 5; meshes → Chapter 6; building and
running the problem → Chapters 7–8. The physics behind everything here is in
Technical Reference §4 (materials & property functions), §5 (effective properties &
mixing), §6 (kinetics), and §10 (volume change & swelling).
