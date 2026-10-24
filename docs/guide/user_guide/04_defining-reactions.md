# 4. Defining Reactions

Reactions tell Pyrolysis.jl how the solid (and any liquid) phases of a material
decompose into char, intermediates, and gases as the temperature rises. Each
reaction is an Arrhenius rate law with a stoichiometry, a heat of reaction, and
an optional temperature window. This chapter shows how to build reactions with
the `Reaction` constructor, how to specify products and yields, how the heat-of-
reaction sign convention works, how to model multi-step schemes (the single-
reactant restriction), how the depletion limiter affects rates near zero
concentration, and how to validate your stoichiometry. The physics and
derivations behind every formula here are in Technical Reference §6 (Reaction
Kinetics and Stoichiometry).

In this chapter the index `i` denotes a reaction and `j` denotes a component, per
the kinetics convention in the nomenclature.

---

## 4.1 The reaction rate law

Every reaction in Pyrolysis.jl evaluates the same volumetric rate (kg·m⁻³·s⁻¹):

```
r_i = A_i · exp(−E_i / (R_g · T)) · Π_j b(ξ_j)^{n_{i,j}} · g(ξ_j/ξ_th) · gate(T)
```

where

- `A_i` is the Arrhenius pre-exponential factor (units depend on order; see §4.4),
- `E_i` is the activation energy (J·mol⁻¹),
- `R_g = 8.31446261815324` J·mol⁻¹·K⁻¹ is the universal gas constant,
- `T` is the cell's trial temperature (K): the cell temperature smoothly
  clamped into the property bracket `[200, 3000]` K with 1 K transitions, which
  is the cell temperature on every physical state,
- `ξ_j` is the mass concentration of reactant `j` (kg·m⁻³) — note
  **mass concentration**, not mole fraction or mass fraction,
- `n_{i,j}` is the reaction order with respect to reactant `j` (must be ≥ 0;
  negative orders are rejected at construction),
- `b(ξ_j)` is the finite power base of the depletion limiter (equal to `ξ_j`
  to within `10⁻⁸` kg·m⁻³) and `g(ξ_j/ξ_th)` its one-sided onset, exactly zero
  at and below zero concentration (§4.10), and
- `gate(T)` is a smooth temperature window built from `T_min` / `T_max` (§4.6).

The rate is on a **mass-concentration** basis, which is what makes the species
sources directly additive into the finite-volume residual. See Technical
Reference §6.1 for the derivation.

---

## 4.2 The `Reaction` constructor

`Reaction` is exported by `using Pyrolysis`. There are four call forms,
distinguished by whether you name components by **symbol** or by **integer index**,
and whether the reaction has **one product** or **two products**.

### Symbolic vs. numeric component references

You can refer to components either by the `Symbol` you gave them in their
constructor, or by their 1-based position in the material's `components` tuple.

```julia
using Pyrolysis

# Symbolic — components are referenced by name; resolved when the Material is built.
Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)

# Numeric — :virgin is component 1, :MMA is component 2 in the components tuple.
Reaction(:decomposition, 1 => 2; A = 8.5e12, E = 188e3, h = 870e3)
```

The symbolic form returns a lightweight `PendingReaction`; `Material(...)`
resolves the names to indices automatically using the `components` tuple, so you
never construct a `PendingReaction` yourself. If a reaction references a name that
is not in the components tuple, `Material(...)` raises a clear error. Prefer the
**symbolic** form: it is self-documenting and immune to component reordering.

### Single-product form

```
Reaction(name, reactant => product; A, E, h, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)
```

The product receives the full mass yield of 1.0 automatically; you do not pass
`yields`.

```julia
# virgin PMMA fully volatilizes to MMA gas, endothermic.
Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3, n = 1.0)
```

### Two-product form

```
Reaction(name, reactant => (product1, product2); A, E, h, yields, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)
```

Here `yields = (y1, y2)` gives the mass split between the two products and **must
sum to 1.0** (see §4.3).

```julia
# Charring: 1 kg virgin → 0.2 kg char + 0.8 kg gas, endothermic.
Reaction(:charring, :virgin => (:char, :gas); A = 1e15, E = 200e3, h = 500e3,
         yields = (0.2, 0.8))
```

### Keyword arguments and defaults

| Keyword | Meaning | Units | Default |
|---|---|---|---|
| `A` | Arrhenius pre-exponential factor (required) | s⁻¹ (first order) | — |
| `E` | Activation energy (required) | J·mol⁻¹ | — |
| `h` | Heat of reaction (required; number, `(a,b)` tuple, or property function) | J·kg⁻¹ | — |
| `yields` | Mass split for the two-product form (required there) | – | — |
| `n` | Reaction order w.r.t. the reactant | – | `1.0` |
| `T_min` | Lower temperature gate | K | `0.0` |
| `T_max` | Upper temperature gate | K | `Inf` |
| `gate_width` | Width of the tanh ramp of both gates (finite, > 0) | K | `1.0` |
| `validate_mass` | Enforce mass balance at construction | – | `true` |

`A`, `E`, and `h` are **required** keyword arguments — there are no defaults, and
omitting any of them is an error. `yields` is required only for the two-product
form. All arguments after the component pair are keyword arguments (note the `;`
in the signatures above), though the examples also work with `,` because Julia
accepts trailing keyword arguments either way.

---

## 4.3 Products, yields, and mass balance

Pyrolysis.jl uses a **single-reactant** convention (see §4.7). Per unit reaction
extent, the one reactant is consumed with mass stoichiometry `−1.0` (this is
implicit — you never type it), and the products are generated with the mass
stoichiometries you supply. Mass is conserved when the product yields sum to one:

```
Σ_products ν_{i,j} = 1.0
```

For the single-product form, the lone product automatically gets yield `1.0`, so
mass balance is trivially satisfied. For the two-product form, your `yields`
tuple must sum to `1.0`:

```julia
yields = (0.2, 0.8)   # 0.2 + 0.8 = 1.0  ✓  (mass conserved)
yields = (0.2, 0.7)   # 0.2 + 0.7 = 0.9  ✗  (validate_mass throws an error)
```

When `validate_mass = true` (the default), the constructor runs an internal
stoichiometry check and **throws an error at construction time** if the masses do
not balance. This catches typos before you ever solve. You can pass
`validate_mass = false` to skip the check (not recommended — it lets mass
conservation silently break).

The species source assembled per cell is the stoichiometric sum
`S_j = Σ_i ν_{i,j} r_i` (kg·m⁻³·s⁻¹), with the reactant carrying `−r_i` and each
product carrying `+yield · r_i`. See Technical Reference §6.5.

---

## 4.4 The Arrhenius parameters `A` and `E`

- `A` (pre-exponential factor) and `E` (activation energy, J·mol⁻¹) are stored as
  plain scalars (`Float64`). They are **temperature-independent**: a
  temperature-dependent `A(T)` is not supported directly (see §4.14).
- The units of `A` depend on the reaction order. For the standard first-order
  case (`n = 1.0`) the rate constant `A·exp(−E/R_gT)` has units s⁻¹ and the rate
  `r` comes out in kg·m⁻³·s⁻¹ because it is multiplied by a concentration. For a
  bimolecular convention `A` would be m³·kg⁻¹·s⁻¹. **The constructor does not
  check these units** — getting them right is your responsibility.

```julia
Reaction(:pyrolysis, :virgin => (:char, :gas);
         A = 1.0e13,    # s⁻¹
         E = 150e3,     # J/mol
         h = 300e3,     # J/kg
         yields = (0.20, 0.80))
```

§4.9 gives a practical recipe for recovering `A` and `E` from a TGA reference
temperature and peak rate, FDS-style.

---

## 4.5 Reaction order `n`

The keyword `n` sets the order with respect to the reactant concentration: the
rate carries the factor `ξ_reactant^n`. The default is first order (`n = 1.0`),
which is the usual choice for solid-state thermal decomposition. Non-integer
orders are allowed.

```julia
# Half-order decomposition
Reaction(:decomp, :virgin => :gas; A = 1e12, E = 180e3, h = 500e3, n = 0.5)
```

The rate evaluation has a fast path for orders 0 and 1 (avoiding a general
power), so first-order reactions are both the most common and the cheapest.

---

## 4.6 Temperature gates `T_min` and `T_max`

`T_min` and `T_max` define the temperature window in which the reaction is active.
They are **not hard cutoffs**: instead of step functions, Pyrolysis.jl multiplies
the rate by smooth `tanh` ramps,

```
r ← r · ½(1 + tanh((T − T_min)/w_gate)) · ½(1 + tanh((T_max − T)/w_gate))
```

with `w_gate = gate_width` (default 1 K), so the reaction turns on and off
smoothly. A narrower width approaches the hard step at the cost of a stiffer
Jacobian near the gate. The smoothness is deliberate — hard
steps create discontinuities in the Jacobian that break automatic differentiation
and Newton convergence (Technical Reference §6.2).

```julia
Reaction(:pyrolysis, :virgin => (:char, :gas);
         A = 1.0e13, E = 150e3, h = 300e3, yields = (0.20, 0.80),
         T_min = 400.0,    # ramps on near 400 K
         T_max = 1500.0)   # ramps off near 1500 K
```

Defaults are `T_min = 0.0` and `T_max = Inf`. The lower gate is applied only when
`T_min` is finite **and** strictly positive, and the upper gate only when `T_max`
is finite; with the defaults, neither gate fires and the bare Arrhenius rate is
used at all temperatures.

> **Caveat — gate leakage.** The ramp `tanh((T − T_min)/w_gate)` is sharp on the
> scale of the gate width but not instantaneous: the reaction is "weakly active"
> within roughly ±10 gate widths of each bound (±10 K at the default width). If you need a
> hard kinetic threshold, the gates are an approximation; for the AD- and
> Newton-stability reasons above, they are intentional. See Technical Reference
> §6.2.

---

## 4.7 The single-reactant restriction and multi-step schemes

The current implementation supports **single-reactant reactions only**
(`NRin == 1`). Each `Reaction` consumes exactly one component. Attempting to apply
a multi-reactant reaction is rejected by problem validation. This means you cannot
write a bimolecular reaction such as `char + O₂ → ...` as a single `Reaction` with
two reactants.

The way to model realistic decomposition is a **sequence of single-reactant
reactions**, each consuming one component and producing intermediates plus gases.
This is exactly how the shipped wood and pressure-treated-wood models work: a
chain `virgin → int1 → int2 → ... → char`, with a gas split released at every
step. Below is a three-step charring chain (symbolic names; numeric indices would
read identically):

```julia
using Pyrolysis

virgin = SolidComponent(:virgin, ρ = 530.0, c = 1500.0, k = 0.12)
int1   = SolidComponent(:int1,   ρ = 400.0, c = 1400.0, k = 0.11)
int2   = SolidComponent(:int2,   ρ = 250.0, c = 1300.0, k = 0.12)
char   = SolidComponent(:char,   ρ = 130.0, c = 1100.0, k = 0.10, ε = 0.95)
gas1   = GaseousComponent(:gas1, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
gas2   = GaseousComponent(:gas2, M = 0.040, c = 1500.0, k = 0.03, λ = 1e-5)
gas3   = GaseousComponent(:gas3, M = 0.050, c = 1500.0, k = 0.03, λ = 1e-5)

material = Material(
    name = :CharringWood,
    components = (virgin, int1, int2, char, gas1, gas2, gas3),
    reactions = (
        Reaction(:step1, :virgin => (:int1, :gas1);
                 A = 7.48e4,  E = 71.6e3,  h = -1.92e4, yields = (0.95, 0.05)),
        Reaction(:step2, :int1   => (:int2, :gas2);
                 A = 2.36e10, E = 140.9e3, h =  2.64e5, yields = (0.80, 0.20)),
        Reaction(:step3, :int2   => (:char, :gas3);
                 A = 6.40e0,  E = 56.0e3,  h = -3.3e2,  yields = (0.74, 0.26)),
    ),
)
```

Each step is mass-balanced on its own (every `yields` sums to 1.0), and the chain
as a whole conserves solid + gas mass.

> **Note.** O₂-dependent char oxidation (an exothermic surface reaction common in
> flaming-combustion models) cannot currently be expressed as a true bimolecular
> reaction. You can approximate it with a single-reactant exothermic char-→-ash
> step gated by temperature (§4.8), but the explicit oxygen dependence is a
> documented limitation (Technical Reference §6.10).

---

## 4.8 Heat of reaction `h` and its sign convention

The heat of reaction is supplied through the `h` keyword in J·kg⁻¹ (per kg of the
reactant). Pyrolysis.jl uses the **storage (enthalpy) sign convention**:

- **`h > 0` is endothermic** — the reaction absorbs heat and **cools** the
  material.
- **`h < 0` is exothermic** — the reaction releases heat and **heats** the
  material.

The volumetric heat source assembled into the energy equation is
`Q_rxn = −Σ_i h_i r_i` (W·m⁻³), so the minus sign converts the stored enthalpy
into a source: an endothermic reaction (`h > 0`) gives `Q_rxn < 0` (cooling), and
an exothermic reaction (`h < 0`) gives `Q_rxn > 0` (heating). See Technical
Reference §6.6.

> **Important — opposite to ThermaKin.** ThermaKin2Ds publishes heats of
> reaction with the *opposite* sign (there `h > 0` is exothermic; Gpyro and FDS
> use the sign of this package). When you port parameters from a ThermaKin input file,
> **flip the sign of every heat of reaction.** The shipped examples annotate this
> conversion explicitly; e.g. a ThermaKin "+2.78E6 J/kg" (exothermic there)
> becomes `h = -2.78e6` here (exothermic in this convention) — and a ThermaKin
> "−2.78E6 J/kg" (endothermic there) becomes `h = 2.78e6` here (endothermic in
> this convention). In both cases the rule is the same: flip the sign.

### Constant heat (the common case)

Pass a plain number; it is wrapped as a `ConstantProperty`:

```julia
# Endothermic pyrolysis (absorbs heat).
Reaction(:decomposition, :virgin => :gas; A = 1.5e14, E = 2.03e5, h = 820e3)

# Exothermic char-forming step (releases heat) — note the negative h.
Reaction(:exo_step, :int4 => (:char, :gas);
         A = 6.4e0, E = 56e3, h = -3.26e2, yields = (0.74, 0.26))
```

### Temperature-dependent heat

Pass a two-tuple `(a, b)` for an affine `h(T) = a + b·T` (a `LinearProperty`), or
any of the property-function types from Chapter 3:

```julia
# h(T) = 870e3 + 100*T  J/kg  (h > 0 = endothermic)
Reaction(:decomp, :virgin => :gas; A = 8.5e12, E = 188e3, h = (870e3, 100.0))
```

### State-dependent heat (moisture sorption)

For moisture evaporation, the heat of sorption depends on the current moisture
content, not only on temperature. Use a `StateDependentProperty`, whose closure
receives both the temperature `T` and the full concentration tuple `ξ`. The
kinetics evaluator automatically passes `ξ` when the heat is state-dependent.

```julia
const MOISTURE_IDX = 7        # index of moisture in the components tuple
const ρ_DRY_BASIS  = 530.0    # dry-wood density for the moisture-content ratio [kg/m³]

# HoS(MC) = 2450.4 + 875.6 * exp(-0.1122 * MC) [kJ/kg], MC = moisture content (% dry basis)
heat_of_sorption = StateDependentProperty(
    (T, ξ) -> begin
        MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS          # moisture content, % dry basis
        return (2450.420909987441 + 875.6441221735336 *
                exp(-0.1121684273316275 * MC)) * 1e3          # J/kg (endothermic, h > 0)
    end,
)

Reaction(:moisture_evaporation, :moisture => :water_vapor;
         A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0)
```

The closure indexes `ξ` by component position, so `MOISTURE_IDX` must match the
moisture component's slot in the `components` tuple. Evaporation is endothermic,
so the returned value is positive. See Technical Reference §6.6 and Chapter 3 for
the property-function types.

---

## 4.9 Recovering `A` and `E` from TGA data (FDS-style recipe)

A common practical task is fitting Arrhenius parameters to a single-step
thermogravimetric analysis (TGA) curve. The FDS approach uses two readouts from a
constant-heating-rate TGA run: the **peak-rate temperature** `T_p` (K, where the
mass-loss rate is maximal) and the **peak reaction rate** `r_p` (s⁻¹, the
maximum of the normalized mass-loss rate `−(1/m_0) dm/dt`; a curve reported per
kelvin is multiplied by the heating rate `β_h` in K·s⁻¹). For a first-order
reaction of a pure material the standard estimates are

```
E ≈ e · r_p · R_g · T_p² / β_h
A ≈ e · r_p · exp(E / (R_g · T_p))
```

where `e ≈ 2.71828` is Euler's number and `R_g = 8.314462618` J·mol⁻¹·K⁻¹
(Technical Reference §6.9.3 gives the general form with the initial mass
fraction). As a
runnable example, suppose the TGA peak is at `T_p = 640 K` with a heating rate of
10 K·min⁻¹ and a peak fractional mass-loss rate of 0.012 K⁻¹:

```julia
const R_g = 8.314462618    # J/(mol·K)

β_h = 10.0 / 60.0          # heating rate, K/s  (10 K/min)
T_p = 640.0                # peak-rate temperature, K
r_p = 0.012 * β_h          # peak reaction rate, s⁻¹  (0.012 K⁻¹ × β_h)

E = exp(1) * r_p * R_g * T_p^2 / β_h   # J/mol  (≈ 111.1 kJ/mol)
A = exp(1) * r_p * exp(E / (R_g * T_p)) # s⁻¹    (≈ 6.3e6 s⁻¹)

rxn = Reaction(:decomposition, :virgin => :gas; A = A, E = E, h = 500e3)
```

These are first-pass estimates (a first-order TGA integrated with these values
peaks at 640 K with a peak rate about 9 % above `r_p`, the error of the
approximation behind the formulas); for production fits, refine `A`, `E`, and `h`
against the measured curve using the sensitivity / inverse-analysis tooling (User
Guide Chapter 13). See Technical Reference §6.9 for the relationship to the FDS,
ThermaKin, and Gpyro kinetic formulations.

---

## 4.10 The depletion limiter

As a reactant concentration approaches zero, a bare Arrhenius rate can still try
to react at a finite speed, which makes implicit solvers stiff — and for a
zero-order fit it would keep consuming reactant straight through zero. An
implicit solver also evaluates trial states on which a concentration is
slightly negative. Pyrolysis.jl therefore multiplies the rate, per reactant, by
a one-sided factor and evaluates the order's power on a finite base:

```
r ∝ b(ξ)^n · g(ξ / threshold)

g(x) = 0                                   for x ≤ 0
g(x) = (s − 2y) t³ + (3y − s) t²           for 0 < x < join,  t = x/join,
                                           y = tanh(join), s = join·sech²(join)
g(x) = tanh(x)                             for x ≥ join

b(ξ) = base_scale · ln(1 + e^{ξ/base_scale})
```

The onset `g` has zero value and zero slope at `x = 0` and matches `tanh` and
its slope at `x = join`, so it is C¹ everywhere. The rate and every derivative
of it are exactly zero at and below zero concentration, for every order
`n ≥ 0`, including zero order: an absent or negative trial concentration is
never consumed. The base `b(ξ)` equals `ξ` to within `base_scale · ln 2` and is
finite for every real `ξ`, so the derivative of `ξⁿ` stays finite at
depletion for fractional orders. Far above the threshold the factor rounds to
one and the bare Arrhenius rate is recovered: the relative deficit is
`2 e^{−2ξ/threshold}`, below 4 % at `2·threshold` and below 5e-9 at
`10·threshold`, and the factor is exactly 1 from `20·threshold`.

The limiter is configured by the `DepletionLimiter` type (exported at the top
level, and by the `Physics` and `Materials` submodules):

```julia
DepletionLimiter(; threshold = 1.0, join = 1.0, base_scale = 1e-8, enabled = true)
```

| Field | Meaning | Units | Default |
|---|---|---|---|
| `threshold` | Concentration scale of the roll-off | kg·m⁻³ | `1.0` |
| `join` | Where the onset joins `tanh`, in units of `threshold` | – | `1.0` |
| `base_scale` | Scale of the finite power base | kg·m⁻³ | `1e-8` |
| `enabled` | Apply the onset factor | – | `true` |

`threshold`, `join` and `base_scale` must be positive (an `ArgumentError`
otherwise). With `enabled = false` the onset factor is one and only the finite
base remains: an order-zero reaction then runs at its full rate at any
concentration, and an absent reactant reacts at the rate of a `base_scale`
concentration; the integrator's admissibility check is then the only guard.

The limiter is a property of the **material**, not of individual reactions or
of `solve`: the default `DepletionLimiter()` is applied to every rate the
solver evaluates, and a custom one is passed with the `Material` constructor's
`depletion_limiter` keyword,

```julia
material = Material(name = …, components = …, reactions = …,
                    depletion_limiter = DepletionLimiter(threshold = 0.1))
```

The reaction kernel reads `material.depletion_limiter` for every rate
evaluation, and with an enabled limiter it skips a cell in which every reaction
lacks its reactant (that cell's rates and their derivatives are exactly zero).
**When to tune:** the default `threshold = 1.0` kg·m⁻³ is appropriate for most
condensed-phase materials, whose concentrations are hundreds to thousands of
kg·m⁻³. A lower threshold keeps the kinetics closer to the bare rate law near
depletion, which matters for genuinely dilute species (e.g. trace moisture at a
few kg·m⁻³); a higher one damps the stiffness of a vanishing reactant. See
Technical Reference §6.3.

You can evaluate a single rate with a custom or disabled limiter for testing,
and the onset itself with `one_sided_factor`:

```julia
using Pyrolysis
using Pyrolysis.Physics: reaction_rate, one_sided_factor   # exported by Physics only

rxn = Reaction(:decomp, 1 => 2; A = 1e10, E = 150e3, h = 5e5)
r      = reaction_rate(rxn, (0.5, 0.0), 700.0)   # the default limiter
r_bare = reaction_rate(rxn, (0.5, 0.0), 700.0; limiter = DepletionLimiter(enabled = false))
one_sided_factor(2.0, 1.0)                        # g(2) = tanh(2) ≈ 0.964
```

`reaction_rate(rxn, ξ, T; limiter)` takes a resolved (index-form) reaction, the
concentration `NTuple` of the cell and its temperature; without `limiter` it
uses `DepletionLimiter()`, not the material's.

---

## 4.11 Assembling reactions into a material

Reactions are passed to `Material(...)` as a tuple via the `reactions` keyword.
The components tuple comes first (so symbolic names can be resolved), and the
reaction tuple may mix heterogeneous reaction types (different orders, product
counts, heat-property types) without penalty:

```julia
using Pyrolysis

material = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86),
        SolidComponent(:char,   ρ = 5.0,  c = 1000, k = 0.10, ε = 0.95),
        GaseousComponent(:MMA,  M = 0.10012, c = 1500, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, :virgin => (:char, :MMA);
                 A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0,
                 yields = (0.002, 0.998)),   # 0.2% char, 99.8% gas
    ),
)
```

A material can also have **no reactions** — pass `reactions = ()` (the default) for
a pure heat-conduction problem. `Material(...)` resolves symbolic reactions (a
name that matches no component throws) into `Reaction`s, whose constructor
checks the mass balance; integer component indices are range-checked when
`solve` validates the problem (§3.7). Chapter 3 covers the rest of the
`Material` constructor (components, mixing rules, lateral-shrinkage law).

---

## 4.12 Validating reactions

One helper lets you check mass balance explicitly, beyond the automatic
check at construction. It is **not** part of the top-level `Pyrolysis`
exports — reach it through its submodule as
`Pyrolysis.Materials.verify_stoichiometry`, or bring it into scope with
`using Pyrolysis.Materials: verify_stoichiometry`.

**`verify_stoichiometry(rxn; tol=1e-10) -> Bool`** confirms that a single
reaction's product yields sum to 1.0 (within `tol`). It returns `true` when
balanced and issues a warning and returns `false` otherwise:

```julia
using Pyrolysis.Materials: verify_stoichiometry

rxn = Reaction(:charring, 1 => (2, 3);
               A = 1e15, E = 200e3, h = 500e3, yields = (0.25, 0.75))
verify_stoichiometry(rxn)   # true: 0.25 + 0.75 = 1.0

all(verify_stoichiometry, material.reactions)   # every reaction of a material
```

It takes a resolved `Reaction`: a reaction built with symbols is a
`PendingReaction` until `Material(...)` resolves it, so check the material's
`reactions` (or use the index form) rather than the symbolic object itself.

It is a diagnostic: mass conservation is enforced by the single-reactant
construction convention (§4.3), not by a projection during the solve, so the
check is for catching setup mistakes. See Technical Reference §6.8.

---

## 4.13 Worked example: a single-step endothermic decomposition

Putting it together — a one-reaction PMMA-like material, evaluated at a single
cell state to show the rate, heat source, and species sources:

```julia
using Pyrolysis
# `reaction_rate` lives in the Physics submodule and is exported from it, but is
# not re-exported at the top level — bring it into scope explicitly.
using Pyrolysis.Physics: reaction_rate
# `verify_stoichiometry` is not exported at all; reach it through the Materials
# submodule.
using Pyrolysis.Materials: verify_stoichiometry

# Three components: virgin solid, char residue, gas product.
virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
char   = SolidComponent(:char,   ρ = 5.0,  c = 1000, k = 0.10, ε = 0.95)
gas    = GaseousComponent(:MMA,  M = 0.10012, c = 1500, k = 0.02, λ = 1e-5)

rxn_spec = Reaction(:decomposition, :virgin => (:char, :MMA);
               A = 1.5e14,     # s⁻¹
               E = 2.03e5,     # J/mol
               h = 820e3,      # J/kg  (positive ⇒ endothermic ⇒ cools)
               n = 1.0,
               yields = (0.002, 0.998))

material = Material(name = :PMMA, components = (virgin, char, gas), reactions = (rxn_spec,))
rxn = material.reactions[1]                          # the resolved reaction

# Verify mass balance.
@assert verify_stoichiometry(rxn)                  # 0.002 + 0.998 = 1.0
@assert all(verify_stoichiometry, material.reactions)

# Evaluate the rate at one cell state (concentrations in kg/m³, T in K).
T_cell = 700.0
ξ_cell = (1000.0, 0.0, 0.0)                         # (virgin, char, MMA)
r = reaction_rate(rxn, ξ_cell, T_cell; limiter = material.depletion_limiter)  # kg·m⁻³·s⁻¹ (≈ 106.7)

Q_rxn = -rxn.heat(T_cell, ξ_cell) * r              # W·m⁻³  (= -h·r); endothermic ⇒ < 0
S = (-r, 0.002 * r, 0.998 * r)                     # species sources (virgin, char, gas)
@assert isapprox(S[2] + S[3], r; rtol = 1e-12)     # product mass = consumed mass
```

`r` is the volumetric reaction rate, `Q_rxn = −h·r` is the volumetric heat source
(negative here because the reaction is endothermic), and the species sources `S`
conserve mass (the products' generation balances the reactant's consumption).
The solver performs this computation in every cell with the material's
limiter at the cell's trial temperature and multiplies the sources by the cell
volume, since the component masses are the states; see Technical Reference §6.7.

---

## 4.14 Summary and limitations

- A `Reaction` is `A·exp(−E/R_gT)·ξ^n`, times the limiter's onset, with smooth
  `T_min`/`T_max` gates of width `gate_width` and a heat of reaction `h`.
- Reference components by **symbol** (recommended) or integer index; product
  `yields` must sum to **1.0**, enforced at construction unless
  `validate_mass = false`.
- **Sign convention: `h > 0` is endothermic, `h < 0` is exothermic** — opposite
  to ThermaKin; flip the sign when porting ThermaKin parameters.
- Only **single-reactant** reactions are supported; build multi-step schemes as
  chains of single-reactant reactions. Multi-reactant (e.g. O₂-dependent)
  kinetics, temperature-dependent `A`, and pressure-dependent rates are not
  supported (Technical Reference §6.10).
- The **depletion limiter** (default `threshold = 1.0` kg·m⁻³) is applied
  automatically: each reactant's rate factor is exactly zero at and below zero
  concentration and joins `tanh(ξ/threshold)` above it. It is rarely tuned;
  lower it only for genuinely dilute species and raise it only to relieve
  stiffness near full conversion.
- Re-check a reaction's mass balance with `verify_stoichiometry`.

For the mathematics behind the rate law, the smooth gates, the depletion limiter,
the stoichiometric sources, and the heat-source assembly, see Technical Reference
§6 (Reaction Kinetics and Stoichiometry). For defining the components a reaction
acts on, see Chapter 3 (Defining Materials, Components, and Properties).
