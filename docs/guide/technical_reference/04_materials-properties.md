# 4. Materials, Components, and Property Functions

This chapter specifies the *constitutive layer* of Pyrolysis.jl: how a pyrolyzing
substance is decomposed into discrete phase-typed **components**, how each
component carries its thermophysical **property functions**, how those functions
are evaluated and differentiated, and how components and reactions are assembled
into the static `Material` object that every physics kernel reads. The governing
PDE system (Chapter 3, *Governing Equations*) closes only once these laws are
supplied; the homogenization of per-component properties into mixture (effective)
properties is the subject of Chapter 5, *Effective Properties and Mixing Rules*,
to which this chapter forwards repeatedly. Throughout, the spatial coordinate is
$z$ (axial, through-thickness, $z=0$ substrate / bottom, $z=L$ exposed surface),
the primary species state variable is the **mass concentration** $\xi_j$ [kg·m⁻³],
and the heat-of-reaction sign convention is $h>0$ endothermic with the volumetric
source $Q_\mathrm{rxn} = -\sum_r h_r r_r$ (overload note H1 of Chapter 2). In the
kinetics passages of this chapter the index $i$ denotes a **reaction** and $j$ a
**component** (overload note I1); elsewhere $i$ is a cell index.

---

## 4.1 The phase model: solids, liquids, and gases

### 4.1.1 Why phase is a *type*, not a field

A pyrolyzing material is a heterogeneous mixture of condensed-phase species
(virgin polymer, char, moisture) and gaseous species (pyrolyzate, water vapor,
oxidizer). The three phases obey structurally different conservation laws
(Chapter 3): condensed species are **Lagrangian** — they do not transport
relative to the matrix, so the mass of a condensed component in a cell that
moves with the skeleton changes only through reaction — whereas gaseous species
flow through the porous matrix and therefore carry a flux divergence in addition
to a reaction source. The gas-vs-condensed dispatch appears in the species
balance, written in the frame of the skeleton ($\mathbf v_s$ its velocity, zero
on a fixed mesh; §3.2.3):

$$
\partial_t\xi_j + \nabla\!\cdot(\mathbf v_s\,\xi_j) =
\begin{cases}
-\,\nabla\!\cdot\mathbf{J}_j \;+\; \omega_j
   & j \text{ is gas},\\[1.2ex]
\omega_j
   & j \text{ solid or liquid,}
\end{cases}
\tag{4.1}
$$

with $\omega_j = \sum_r\nu_{r,j}r_r$ the net production by reaction and
$\mathbf J_j$ the gas flux relative to the skeleton.

Because this branch is taken for every face and every component, it must not
incur a runtime type check. Pyrolysis.jl encodes the phase **at the type
level**: there are three concrete subtypes of `AbstractComponent`, and the phase
predicate is a method that dispatches on the component's type. The predicates
are resolved once, when the problem is set up: the gas-transport kernel stores
the flags `map(is_gas, material.components)`, a tuple of `Bool`s, and forms face
fluxes only for the components it marks, so its inner loop reads a stored flag,
never a type. Where a predicate is called inside an unrolled loop over the
concrete component tuple (a `@generated` or recursive tuple function), the
compiler folds it to a constant. The runtime `Phase` enum (`SOLID=1`,
`LIQUID=2`, `GAS=3`) exists only for printing, logging, and comparison; it is
*not* consulted in hot paths.

```julia
@enum Phase begin SOLID = 1; LIQUID = 2; GAS = 3 end
```

> **Implementation note.** Phase types and predicates are defined in
> `SolidComponent`, `LiquidComponent`, `GaseousComponent`, and `is_solid` /
> `is_liquid` / `is_gas` / `phase` in `src/Materials/components.jl`. Eq. (4.1)
> is realized by `GasTransportKernel` (`src/Kernels/gas_transport.jl`), whose
> `gas` field holds the gas flags, and by the row assembly, which adds each
> component's face fluxes and reaction source to its mass row
> (`src/Assembly/residual.jl`). The matrix phases cross faces only with the
> skeleton, where an equidistribution correction moves the mesh relative to it
> (`RelativeFluxKernel`, §11.6).

### 4.1.2 Component storage

`SolidComponent` and `LiquidComponent` store bulk and intrinsic density,
heat capacity, conductivity, optional superficial matrix transfer,
permeability, emissivity, incident absorptivity, absorption and swelling.
Property laws have concrete type parameters. Participating lumped matrices
require explicit transfer; zero-`γ` condensed phases cannot have nonzero
transfer or permeability. Liquid defaults to zero swelling participation.

`GaseousComponent` stores molar mass, reference-density law, heat capacity,
conductivity, emissivity, absorption and swelling. It has no matrix transfer,
permeability or intrinsic-density field. Its local intrinsic pore density is
derived from the conserved inventory and pore volume, separately from its
reference-density law. Concrete property types preserve specialization through
the heterogeneous component tuple.

---

### 4.1.3 Compile-time phase accessors

All cross-phase property access goes through dispatched accessors that return
phase-appropriate defaults, so a kernel may ask any component for any quantity
without a runtime field-existence test:

| Accessor | Solid / Liquid | Gas |
|---|---|---|
| `is_solid(c)` | `true` / `false` | `false` |
| `is_liquid(c)` | `false` / `true` | `false` |
| `is_gas(c)` | `false` | `true` |
| `phase(c)` | `SOLID` / `LIQUID` | `GAS` |
| `molar_mass(c)` | `0.0` | `c.molar_mass` |
| `get_permeability(c, T)` | `c.permeability(T)` | `0.0` |
| `has_permeability(c)` | `true` | `false` |
| `get_intrinsic_density(c, T)` | `c.intrinsic_density(T)` | `0.0` |

Each is `@inline` and resolves at compile time on the component's concrete type; none branch
on a stored field. The gas's `get_permeability` and `get_intrinsic_density`
return exactly `0.0`, which is the algebraically correct value for the matrix
mixing rules (a gas contributes nothing to permeability or to the solid
skeleton).

> **Implementation note.** Accessors `get_intrinsic_density`, `get_permeability`,
> `has_permeability`, `get_emissivity`, `get_absorptivity`,
> `is_solid/is_liquid/is_gas`, `phase`, `molar_mass`, all in
> `src/Materials/components.jl`. Kernels read the remaining property fields
> (`density`, `heat_capacity`, …) directly.

### 4.1.4 Bulk versus intrinsic (skeletal) density

The single most important constitutive distinction at the component level is
between two densities (overload note A7 of Chapter 2):

- **Bulk (pure-phase) density** $\rho_j$ [kg·m⁻³] — the mass of pure phase $j$ per
  unit of the *bulk volume it would occupy in the mixture*. This is the density
  used to form thermal and matrix transport weights and by the volume-
  change kernels (Chapter 10). Field `density`.

- **Intrinsic / skeletal density** $\rho_{i,j}$ [kg·m⁻³] — the density of the
  cell-wall material itself, ignoring any pore volume internal to the component.
  This is used by **one formula only**, the porosity (Chapter 5, §5.3):

  $$
  \varphi = 1 - \sum_{j\in\{S,L\}} \frac{\xi_j}{\rho_{i,j}} .
  \tag{4.2}
  $$

  The formula carries no clamp: the state mapping floors it smoothly at
  $10^{-12}$ where the pressure divides by it, and the integrator's domain test
  reads it unguarded (§3.4.2, §4.2.3). Field `intrinsic_density`.

When a component is constructed without an explicit intrinsic density, the field
defaults to the bulk density, $\rho_{i,j}=\rho_j$: the pure phase is assumed to
have no internal pores. The resolution logic enforces that a
caller pass *either* `ρ_intrinsic` *or* a pure-phase porosity $\varphi$ (from which
$\rho_{i,j}=\rho_j/(1-\varphi)$ is computed, requiring a constant bulk $\rho_j$),
never both:

```
1. both ρ_intrinsic and φ given        → ArgumentError
2. ρ_intrinsic given                    → to_property(ρ_intrinsic)
3. φ given (requires ρ::Real, 0≤φ<1)    → ConstantProperty(ρ / (1 − φ))
4. neither                              → fall back to bulk ρ
```

> **Implementation note.** `_resolve_intrinsic_density`, `src/Materials/components.jl`.
> The `φ` path requires `raw_ρ isa Real` and `0 ≤ φ < 1`, throwing an
> `ArgumentError` otherwise. Gas components do not store an intrinsic density;
> `get_intrinsic_density(::GaseousComponent, ::Real) = 0.0`, which is consistent
> with Eq. (4.2) summing only over solid/liquid.

### 4.1.5 The swelling factor

The swelling factor $\gamma_j\in[0,1]$ controls whether changes in component $j$'s
concentration change the *bulk volume* of the cell. Physically it answers: when
this component is created or destroyed, does the material expand/shrink, or does
the mass simply enter/leave the existing pore space?

- $\gamma_j = 1$ (default for solids): full contribution. Depletion of the solid
  shrinks the bulk volume; this is the charring/recession behavior.
- $\gamma_j = 0$ (default for gases): no contribution. Gas fills existing pores
  and escapes without expanding the bulk.
- $\gamma_j \in (0,1)$: partial contribution, the intumescent case, where evolved
  gas creates new bulk volume.

$\gamma_j$ weights the volumetric strain rate
$\theta = \sum_j\gamma_j\rho_j^{-1}\mathrm{d}\xi_j/\mathrm{d}t$ (Chapter 10)
and the thermal matrix volume fractions. Transport weights use the same bulk
occupancy but exclude gases explicitly. Every condensed component contributes
to skeletal pore occupancy even when its $\gamma_j$ is zero (§9).

The supported range $0\le\gamma_j\le1$ is a constitutive restriction:
$\gamma_j$ represents the retained fraction of the declared bulk occupancy.
Normalized property weights still sum to one for any positive denominator;
normalization alone does not establish this range. Expansion can be represented
by a lower product bulk density: a complete conversion with unit swelling
factors gives $V_f/V_0=\sum_jy_j\rho_A/\rho_j$.

---

## 4.2 Component property fields and their roles

### 4.2.1 Inventory

Each property field of a component, its symbol, its SI units, and where it is
consumed downstream:

| Field | Symbol | Units | Phase | Primary consumer |
|---|---|---|---|---|
| `density` | $\rho_j$ | kg·m⁻³ | all | mixing rules, $\rho$ (§4.2.2), volume change (Ch. 10) |
| `intrinsic_density` | $\rho_{i,j}$ | kg·m⁻³ | S, L | porosity $\varphi$ (§4.2.3, Eq. 4.2) |
| `heat_capacity` | $c_{p,j}(T)$ | J·kg⁻¹·K⁻¹ | all | heat-capacity sum $C_i$ and sensible enthalpies $h_j(T)$ (§4.4.10): energy rows, gas enthalpy fluxes, deactivation transfer, energy ledger (Ch. 5, 7, 14, 16) |
| `thermal_conductivity` | $k_j(T)$ | W·m⁻¹·K⁻¹ | all | $k_\mathrm{eff}$ (Ch. 5), Fourier flux (Ch. 7) |
| `gas_transfer_coeff` | $\lambda_j(T)$ | m²·s⁻¹ | S, L | $\lambda_\mathrm{eff}$ (Ch. 5), gas diffusion (Ch. 9) |
| `permeability` | $\kappa_j(T)$ | m² | S, L | $\kappa_\mathrm{eff}$ (Ch. 5), Darcy velocity (Ch. 9) |
| `emissivity` | $\varepsilon_j(T)$ | – | all | $\varepsilon_\mathrm{eff}$ (Ch. 5), radiative BC (Ch. 12) |
| `surface_absorptivity` | $a_j$ (or `nothing` ⇒ $\varepsilon_j(T)$) | – | S, L | $a_\mathrm{eff}$ (§5.8.2), incident-flux BC (Ch. 12), Beer–Lambert incident intensity (Ch. 8) |
| `absorption_coeff` | $\alpha_j$ | m²·kg⁻¹ | all | $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ (Ch. 5, 8) |
| `swelling_factor` | $\gamma_j$ | – | all | $\rho$, strain rate $\theta$ (Ch. 10) |
| `molar_mass` | $M_j$ | kg·mol⁻¹ | G | $\rho_g$ (§4.3), $P$ (Ch. 9), Chapman–Enskog estimate (§5.10) |

Note that $\alpha_j$ is a **mass-basis** absorption coefficient [m²·kg⁻¹]; it is
multiplied by concentration to give the volumetric (extinction) absorption
coefficient $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ [m⁻¹] used by Beer–Lambert
radiation (Chapter 8), at the positive-part concentrations. Condensed
superficial transfer is explicit under the lumped model and uses condensed
matrix weights. Its value is constitutive data, independent of gas loading.

### 4.2.2 Bulk density, matrix occupancy and gas loading

The bulk mixture mass density is $\rho_{tot}=\sum_j\xi_j$. Matrix transport
uses $q_j=\gamma_j\xi_j/\rho_{b,j}$ only for condensed components and normalizes
by their sum. Gas loading cannot alter these weights, even if gas has nonzero
`γ`. Thermal mixing and the volume-change law retain their own stated bases.
Skeletal occupancy is a separate sum using intrinsic densities and no `γ`.

---

### 4.2.3 Pore volume

$\phi=1-\sum_{S\cup L}\xi_j/\rho_{i,j}(T)$ enters the pore EOS, intrinsic
densities, dusty-gas geometry and physical admissibility. For condensed
species, $\partial\phi/\partial\xi_j=-1/\rho_{i,j}$; for gas this derivative
is zero. Temperature dependence of skeletal density is differentiated too.
A condensed component with `γ=0` still occupies this skeletal volume.
Chapter 9 specifies the common trial extension and accepted-state margin.

---

## 4.3 Gas reference density and the local pore EOS

The gas component's reference density is
$\rho_{ref}(T)=P_{ref}M_j/(RT)$, with $P_{ref}=101325$ Pa. It supplies the
component's thermal-volume and swelling-property law. It is independent of
the local pressure. The actual pore density is $\rho_j^g=\xi_j/\phi$ and the
partial pressure is $p_j=\rho_j^gRT/M_j$ (Chapter 9).

Gas reference density and gas swelling never enter condensed matrix transport
weights. `GaseousComponent(name;M,c,k,ε=0,α=0,γ=0)` declares caloric and thermal
properties; matrix transfer belongs to the gas model and condensed components.

---

## 4.4 Property-function types

A *property function* is a callable object $p$ such that $p(T)$ (or, for
state-dependent properties, $p(T,\xi)$) returns the property value at the given
state. All eight concrete types subtype `AbstractPropertyFunction` and implement
both the `(p)(T)` and `(p)(T, ξ)` call signatures (the latter ignoring $\xi$ for
all but `StateDependentProperty`). This uniformity lets a single component field
hold any property law, and lets the kinetics evaluator pass $\xi$ to a heat-of-
reaction property without knowing whether it is composition-dependent.

**The property bracket.** Every temperature-dependent law evaluates on the
clamped argument $\tilde T = \operatorname{clamp}(T, T_{\min}, T_{\max})$ with
$[T_{\min}, T_{\max}] = [200, 3000]$ K (`property_temperature`;
`T_PROPERTY_MIN`, `T_PROPERTY_MAX` in `src/Materials/constants.jl`). Inside the
bracket the clamp is the identity with unit derivative, so nothing changes for
any physical state; outside it the law sees the nearest bracket end. The kernels
do not hand a law the raw cell temperature: they pass the **trial temperature**
$T_{\mathrm{trial}}$, the cell temperature smoothly clamped into the same bracket
with cubic joins 1 K wide (`trial_temperature`,
`src/Materials/guards.jl`), so the composite $p(T_{\mathrm{trial}}(T))$ is C¹ in the
state. A trial state with a nonphysical (even negative) temperature therefore
gives a finite residual with a finite derivative, not a `DomainError` out of a
fractional power or a fitted polynomial's far extrapolation, and the
integrator's domain test rejects the state (§15.5.4). The surface energy balance
evaluates its condition at the trial temperature of its Newton iterate and never
clamps the iterate itself (§12.2). No derivative is hand-coded: the colored
Jacobian evaluates each law on `ForwardDiff.Dual` arguments (§4.5), so the
derivatives quoted below are what automatic differentiation of the law returns —
and, through the clamp, zero outside the bracket.

The component constructors sample every law on the bracket at 1 K steps and
reject a component whose bulk or intrinsic density, heat capacity or thermal
conductivity is not positive there (a gas may have $c_p = 0$), whose gas
transfer coefficient or permeability is negative, or whose emissivity leaves
$[0, 1]$. The error names the first offending temperature: a linear fit
$k = 0.3 - 4.2\times10^{-4}\,T$ is rejected at 715 K because it turns a cell
into a perfect insulator once it reaches that temperature. Hold a fit constant
outside its measured range (`T -> law(clamp(T, T_lo, T_hi))`) or pass
`validate_properties = false` for laws built programmatically from search
parameters, where an unphysical value should reach the objective function
instead of throwing.

### 4.4.1 `ConstantProperty`

$$
p(T) = \texttt{value}, \qquad \frac{\mathrm{d}p}{\mathrm{d}T}=0.
$$

The workhorse for properties given a single number. It has a *fast path* in the
parameter-evaluation interface (§4.7) that bypasses the dynamic-dispatch hop.

### 4.4.2 `LinearProperty`

$$
p(T) = a + bT, \qquad \frac{\mathrm{d}p}{\mathrm{d}T}=b.
$$

Fields `a` (intercept at $T=0$) and `b` (slope). Constructed automatically from a
two-tuple via `to_property((a,b))`.

### 4.4.3 `PolynomialProperty`

$$
p(T) = \sum_{m=0}^{N-1} c_m\, T^{m} = c_0 + c_1 T + c_2 T^2 + \cdots,
$$

evaluated with `evalpoly` (Horner's scheme) for numerical stability and to avoid
recomputing powers. Forward-mode AD through Horner's scheme returns the exact
derivative

$$
\frac{\mathrm{d}p}{\mathrm{d}T} = \sum_{m=1}^{N-1} m\,c_m\,T^{m-1}
= c_1 + 2c_2 T + 3c_3 T^2 + \cdots .
\tag{4.6}
$$

### 4.4.4 `ArrheniusProperty`

$$
p(T) = A\,\exp\!\left(-\frac{E}{R_g\,T}\right),
\qquad
\frac{\mathrm{d}p}{\mathrm{d}T} = p(T)\,\frac{E}{R_g\,T^2}.
\tag{4.7}
$$

Fields `A` (pre-exponential), `E` (activation energy [J·mol⁻¹]), and `T_ref` (a
reference temperature, defaulting to $T_\mathrm{ref}=298.15$ K, carried for
normalization purposes but not used in the value itself). The law shares its
functional form with the reaction-rate constant (§4.5.3, Chapter 6).

### 4.4.5 `TableProperty`

Tabulated $(T_m, v_m)$ pairs with **linear interpolation** and **flat
extrapolation**. The constructor validates that the temperatures are strictly
increasing (`error` otherwise). For a query $T_q$:

$$
p(T_q) =
\begin{cases}
v_1 & T_q \le T_1,\\[1ex]
v_m + \dfrac{T_q - T_m}{T_{m+1}-T_m}\,(v_{m+1}-v_m)
   & T_m \le T_q < T_{m+1},\\[2ex]
v_N & T_q \ge T_N.
\end{cases}
\tag{4.8}
$$

The interval search is **linear for $N\le 8$** (fewer comparisons than a binary
search on small tables) and **binary for $N>8$** (logarithmic). Under AD the
derivative is the slope of the segment containing $T_q$ (the right-hand segment at
a knot) and zero in the flat extrapolation; it jumps at every knot.

### 4.4.6 `CallableProperty`

Wraps an arbitrary callable $f(T)$. Its value is $p(T)=f(\tilde T)$ on the
clamped argument of the property bracket; its derivative is whatever AD
propagates through $f$, so $f$ must accept `ForwardDiff.Dual` arguments (no
`Float64` annotations or conversions inside). This is the type produced by
`to_property(f)` for a function `f` and is how the ideal-gas density default
(Eq. 4.5) is stored.

### 4.4.7 `SutherlandProperty`

Sutherland's law for a gas transport property, carrying its own three constants:

$$
p(T) = \mu_{\text{ref}}
\left(\frac{T}{T_{\text{ref}}}\right)^{3/2}
\frac{T_{\text{ref}}+S}{T+S}.
$$

`SutherlandProperty()` is the canonical **air** triple
($\mu_{\text{ref}}=1.716\times10^{-5}$ Pa·s, $T_{\text{ref}}=273.15$ K,
$S=110.4$ K) and is the default `Material.gas_viscosity` (§5.9), so a material
that says nothing about its pore gas uses air. Its derivative, as AD returns it, is

$$
\frac{\mathrm{d}p}{\mathrm{d}T} = p(T)\left(\frac{3}{2T} - \frac{1}{T+S}\right),
$$

(the AD-cleanliness of the type is tested in
`test/materials/test_gas_viscosity_property.jl`).

Like every property law it evaluates on the clamped argument $\tilde T$ of the
property bracket. The positivity floor `MU_FLOOR` ($10^{-12}$ Pa·s) is applied
once, centrally, by the state mapping, which forms the Darcy mobility $\kappa/\mu$
with $\mu$ = `smooth_floor(material.gas_viscosity(T_trial), MU_FLOOR)` (§5.9), so
that every curve a material may carry is protected rather than just this one.

### 4.4.8 `StateDependentProperty`

A property that depends on **both** temperature and the full concentration tuple,
$p(T,\xi)=f(T,\xi)$. Calling it with $T$ alone raises an error (state is
mandatory). Its canonical use is a moisture-dependent heat of sorption: the heat
released/absorbed when water desorbs depends on the current moisture content, so

$$
h(T,\xi) = f\!\left(T,\; \mathrm{MC}(\xi)\right),
$$

where $\mathrm{MC}$ is, e.g., the moisture content reconstructed from the
concentration tuple. The kinetics evaluator always passes $\xi$ to the heat
property (§4.5, Chapter 6), so a `StateDependentProperty` heat is evaluated with
the live composition.

Its derivatives in $T$ and in each $\xi_j$ are propagated by AD through $f$
(§4.5), so $f$ must accept a `Dual` temperature and a tuple of `Dual`
concentrations.

### 4.4.9 `to_property` conversions

The convenience converter normalizes user inputs to property functions, so that
component constructors and reaction constructors accept either property objects
or bare values:

| Input | Result |
|---|---|
| an `AbstractPropertyFunction` | passes through unchanged |
| a `Real` number | `ConstantProperty(Float64(x))` |
| a `Tuple{Real,Real}` $(a,b)$ | `LinearProperty(a, b)` |
| a `Function` | `CallableProperty(f)` |

> **Implementation note.** All property types, their call methods,
> `linear_interpolate`, and `to_property` are in
> `src/Materials/property_functions.jl`.

### 4.4.10 Sensible enthalpy of a heat-capacity law

The energy rows, the enthalpy the face fluxes carry, the caloric relation of the
inventory form and the deactivation transfer all read the sensible enthalpy of
each component on one datum (Chapter 3, §3.3; Chapter 7),

$$
h_j(T) = \int_{T_\mathrm{ref}}^{T} c_{p,j}(\tilde T')\,\mathrm{d}T',
\qquad T_\mathrm{ref} = 298.15\ \mathrm{K},
\tag{4.9}
$$

the integral of the law as it evaluates, on the clamped argument, so that
$c_{p,j}$ is constant outside the property bracket and $h_j$ continues linearly
there. Each component's enthalpy law is built once per problem from its
heat-capacity law (`enthalpy_law`, collected in `CaloricLaws`):

- A `ConstantProperty`, `LinearProperty`, `PolynomialProperty` or
  `TableProperty` is integrated in closed form — the antiderivative of the
  clamped law, piecewise across the bracket ends and the table knots (the
  trapezoid rule is exact on each linear piece of a table) — so
  $\mathrm{d}h_j/\mathrm{d}T$ is the law itself to round-off.
- Any other law (`ArrheniusProperty`, `SutherlandProperty`, `CallableProperty`)
  is integrated once, when the problem is set up, onto the 1 K knots $T_k$ of the
  bracket by five-point Gauss–Legendre quadrature on each interval, and is
  evaluated between knots by the cubic Hermite interpolant of the integral $h_k$
  and the law $c_k$ at the knots (`TabulatedEnthalpy`):
  $$
  h(T) = (2s^3-3s^2+1)\,h_k + (s^3-2s^2+s)\,c_k + (-2s^3+3s^2)\,h_{k+1} + (s^3-s^2)\,c_{k+1},
  \qquad s = \frac{T - T_k}{1\ \mathrm{K}} \in [0,1),
  \tag{4.10}
  $$
  which is C¹ in $T$, its derivative given in closed form (`enthalpy_slope`). The
  interpolant reproduces exactly every heat-capacity law of polynomial degree two
  or less: wrapping the linear law $c_p = 1000 + T$ in a `CallableProperty`
  reproduces its closed-form enthalpy to within $3\times10^{-10}$ J·kg⁻¹, on
  enthalpies of order $10^6$ J·kg⁻¹.
  Both the interpolant and its slope are `NaN` at a `NaN` temperature, so a wild
  trial state gives a `NaN` residual, which the integrator rejects, rather than an
  error.
- A `StateDependentProperty` needs the composition and has no temperature-only
  form, so it cannot serve as a heat-capacity law.

> **Implementation note.** `sensible_enthalpy`, `enthalpy_slope`,
> `enthalpy_law`, `TabulatedEnthalpy`, `CaloricLaws` and `component_enthalpies`
> are in `src/States/caloric.jl`.

---

## 4.5 Property derivatives and their role in the Jacobian

Implicit time integration (Chapter 15) requires the Jacobian
$\partial\boldsymbol f/\partial\boldsymbol u$ of the residual. Because effective
properties depend on $T$ and $\xi$, and the residual depends on effective
properties, the per-component derivatives $\mathrm{d}p_j/\mathrm{d}T$ and
$\partial p_j/\partial \xi_k$ propagate through the chain rule into the Jacobian
(Chapter 5, §5.15 for the mixing-rule chain). The property types carry no
derivative methods, and no part of the Jacobian is hand-coded: the colored
Jacobian evaluates the whole residual — the state mapping, where the property
laws are evaluated, and every kernel — on `ForwardDiff.Dual` numbers, a few
columns at a time on a sparsity pattern built from the mesh connectivity
(Chapter 15). The property derivative is therefore the derivative of the law as
written, composed with the trial-temperature clamp and the positive part of the
masses through which the state mapping reaches it (§4.4). A dense forward-mode
Jacobian of the same residual is the reference the tests compare the colored one
against.

### 4.5.1 Temperature derivative $\mathrm{d}p/\mathrm{d}T$

What forward-mode AD returns for each type, inside the property bracket:

| Type | $\mathrm{d}p/\mathrm{d}T$ |
|---|---|
| `ConstantProperty` | $0$ |
| `LinearProperty` | $b$ |
| `PolynomialProperty` | $\sum_{m\ge 1} m c_m T^{m-1}$ (Eq. 4.6) |
| `ArrheniusProperty` | $p(T)\,E/(R_g T^2)$ (Eq. 4.7) |
| `TableProperty` | slope of the active segment; $0$ in the flat extrapolation |
| `CallableProperty` | $f'(T)$, as propagated through `f` |
| `SutherlandProperty` | $p(T)\,(3/(2T) - 1/(T+S))$ |
| `StateDependentProperty` | $\partial f/\partial T$ at fixed $\xi$ |

These values are exact derivatives of the evaluated law, with no truncation
error. A `TableProperty` has a discontinuous derivative at every knot, which the
integrator sees as a kink in the Jacobian; a smooth law (polynomial, Arrhenius,
Sutherland) avoids it.

### 4.5.2 Concentration derivative $\partial p/\partial \xi_j$

For every property type *except* `StateDependentProperty`, the two-argument call
`p(T, ξ)` ignores $\xi$, so $\partial p/\partial \xi_j = 0$ exactly. For a
`StateDependentProperty` the derivative is $\partial f/\partial \xi_j$, propagated by
AD through `f(T, ξ)` with the concentration tuple carrying `Dual` entries. The
kinetics evaluator passes the cell's concentrations to every heat-of-reaction
property, so a moisture-dependent heat of sorption contributes its composition
derivative to the rows of the cell's energy source.

> **Implementation note.** The call methods `(p)(T)` and `(p)(T, ξ)` are in
> `src/Materials/property_functions.jl`. The heats are evaluated with the rates by
> `evaluate_reactions!` (`src/Physics/kinetics.jl`), which `ReactionKernel`
> (`src/Kernels/reactions.jl`) calls per cell; the colored Jacobian
> (`src/ColoredJacobian/`) differentiates that kernel with the rest of the
> residual.

### 4.5.3 The reaction-rate connection

The Arrhenius derivative (Eq. 4.7) is identical in form to the
temperature-derivative of the reaction-rate constant
$k(T)=A\exp(-E/(R_g T))$, since $\mathrm{d}k/\mathrm{d}T = k\,E/(R_g T^2)$; the
Jacobian obtains it by forward-mode differentiation of the reaction kernel
(Chapter 15). The full reaction-rate derivative including the depletion limiter and
temperature gates is developed in Chapter 6.

---

## 4.6 Material assembly

### 4.6.1 The `Material` type

A `Material` is the static, fully type-specialized aggregate that every physics
kernel reads. It is **immutable** and **heterogeneously typed** so the compiler
can specialize every loop over its components and reactions:

```julia
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ,DL<:AbstractDepletionLimiter} <:
       AbstractMaterial
    name::Symbol
    components::C          # heterogeneous component tuple
    reactions::RT          # heterogeneous reaction tuple
    mixing::NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}
    lateral_shrinkage_law::LSL   # A(t)/A₀ = law(χ̄), or nothing
    depletion_limiter::DL  # one-sided rate roll-off near depletion (§6.3)
    gas_viscosity::Pμ      # pore-gas μ(T), default SutherlandProperty() (§5.9)
end
```

The type parameters carry compile-time information:

- `NC` — number of components, a static integer driving every `@generated`
  unroll (`Val(NC)`).
- `NR` — number of reactions, likewise static.
- `C`, `RT` — the *exact* concrete tuple types of the components and reactions.
  Because these are full tuple types (not `Vector{AbstractComponent}`), iterating
  over them is type-stable and allocation-free; the cost is that each distinct
  material *shape* is a distinct compiled type.
- `LSL` — the lateral-shrinkage-law type, a concrete `Function` subtype or
  `Nothing`. Parameterizing it here lets calls through
  `material.lateral_shrinkage_law` specialize statically rather than dispatch
  through a `Union{Function,Nothing}` (Chapter 10). The inner constructor
  enforces `LSL === Nothing || LSL <: Function`.
- `Pμ` — the type of the pore-gas viscosity property function, parameterized for
  the same reason so `material.gas_viscosity(T)` stays monomorphic. The inner
  constructor probes the curve at 200, 300, 500, 800, 1200 and 1800 K and rejects
  one that is not positive and finite at any of them (§5.9).
- `DL` — the concrete type of the depletion limiter, so that rate evaluation
  specializes on it.

The `mixing` field bundles three `MixingSpec`s — one each for thermal
conductivity ($k$), gas transfer ($\lambda$), and permeability ($\kappa$) — in a
fixed-order `NamedTuple` so that `material.mixing.k` compiles to a direct field
read. Each `MixingSpec` carries `(rule, β, ℓ)`: the mixing rule, the WEIGHTED
blend weight $\beta$ (default 0.5), and the Carman–Kozeny characteristic length
$\ell$ (NaN unless the rule is `CARMAN_KOZENY`). The default channel rules are

$$
k:\ \text{WEIGHTED},\ \beta=0.5;\qquad
\lambda:\ \text{PARALLEL};\qquad
\kappa:\ \text{WEIGHTED},\ \beta=0.5,
$$

built by `MaterialMixing()`. The `MixingSpec` constructor validates
$0\le\beta\le1$ and, for `CARMAN_KOZENY`, requires a finite positive
`length_scale`, failing fast otherwise. The semantics of each rule (PARALLEL,
SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY), the BRUGGEMAN rejection for $\lambda$
and $\kappa$, and the Newton solve for BRUGGEMAN, are all developed in Chapter 5.

The `depletion_limiter` field carries the material's own `DepletionLimiter`
(defaults `threshold = 1.0` kg·m⁻³, `join = 1.0`, `base_scale = 1e-8` kg·m⁻³,
`enabled = true`), set via the constructor's `depletion_limiter` keyword and read
by the reaction kernel for every reaction-rate evaluation. It is one-sided: each
reactant's rate factor is exactly zero at and below zero concentration and C¹
through it, and rises to $\tanh(\xi/\xi_\mathrm{th})$ above `join` times the
threshold (§6.3).

> **Implementation note.** `Material`, its inner/outer constructors,
> `_construct_material`, `MixingRule`, `MixingSpec`, `MaterialMixing`, in
> `src/Materials/materials.jl`. Counts via `n_components`, `n_reactions`.

### 4.6.2 `PendingReaction` and symbolic resolution

Reactions may be written against component **symbols** rather than integer
indices, which is far less error-prone for the user:

```julia
Reaction(:decomposition, :virgin => :MMA;          A=8.5e12, E=188e3, h=870e3)
Reaction(:charring,      :virgin => (:char, :gas); A=1e15,  E=2e5, h=5e5, yields=(0.2,0.8))
```

A symbolic constructor produces a `PendingReaction`, which stores the symbol
payloads plus the kinetic parameters. When `Material(...)` is called, it builds a
`name → index` map from the components tuple (erroring on a duplicate component
name) and resolves every `PendingReaction` into a numeric `Reaction`. The
resolution path is taken **only if** at least one reaction is pending — the common
all-numeric case skips the `Dict` allocation entirely:

```
needs_resolution = any(rxn isa PendingReaction for rxn in reactions)
needs_resolution || return reactions          # fast path, no Dict
```

Numeric reactions (built against integer indices directly) bypass resolution.

> **Implementation note.** `PendingReaction`, the symbolic `Reaction` constructors,
> `resolve_reaction`, and `_resolve_reactions` (called from the `Material`
> constructor) in `src/Materials/reactions.jl` and `src/Materials/materials.jl`.

### 4.6.3 The single-reactant mass-balance convention

The `Reaction` type is n-ary in both reactants and products in principle, but the
current model **supports only single-reactant reactions** ($N_\mathrm{in}=1$). The
mass-balance convention (developed fully in Chapter 6) is:

- The lone reactant is consumed with implicit mass stoichiometry $-1$ per unit
  reaction extent.
- Each product $j$ is produced with an explicit mass yield $\nu_j > 0$.
- Mass conservation requires $\sum_{\text{products}} \nu_j = 1$.

The `Reaction` inner constructor enforces this at build time — every
construction path, including the bare positional constructor, runs the
validators. Reaction orders must be $n \ge 0$ (`_validate_reaction_order`);
every product yield must be $\ge 0$ (`_validate_yield`, enforced
unconditionally — a negative yield consumes a species the rate law does not
see and drives its concentration negative without bound; note $(1.5, -0.5)$
sums to 1 and would pass a sum-only check); and for a single-reactant reaction
the product yields must sum to 1 (`_validate_stoichiometry`, relative tolerance
$\max(10^{-10}\cdot(1+\Sigma\nu), 10^{-14})$), raising an `error` on violation.
Passing `validate_mass=false` disables only the sum check (not recommended);
the non-negativity and order checks always run.

> **Implementation note.** `_validate_stoichiometry`, `verify_stoichiometry`
> (tol $10^{-10}$, returns `false` + `@warn` on failure), and `mass_balance_error`
> in `src/Materials/reactions.jl`.

### 4.6.4 Dry-pyrolysis-gas identification

For variable-cross-section (lateral shrinkage, Chapter 10) the cumulative gas
$G_i$ of the conversion-history pair, and with it the pyrolysis progress
$\chi_i = G_i/m_{\mathrm{dry},i}$, counts only **dry** pyrolysis gas — water vapor is
excluded because moisture evaporation should not drive radial shrinkage. The set
of dry-gas component indices is computed by a recursive tuple walk that includes a
component index iff it is a gas **and** its name is not a water-vapor name. The
water-vapor test is a substring/equality heuristic on the component name (it
contains `H2O`, contains `water` in any letter case, or is `h2o` in any letter
case):

```
is_gas(c) && !_is_water_vapor_name(c.name)  ⇒ include index
```

This name-based identification is a documented heuristic (open issue §4.9); a
gas whose name does not contain a recognized water token is treated as dry
pyrolyzate regardless of its actual chemistry.

> **Implementation note.** `dry_pyrolysis_gas_indices` →
> `_collect_dry_pyrolysis_gas_indices` and the heuristic `_is_water_vapor_name`,
> `src/Materials/materials.jl`. The phase test dispatches on the component types;
> the name test runs when the indices are computed, once, as the problem is set
> up: `LateralShrinkage` (`src/Closures/lateral.jl`) stores the result as the flag
> tuple `dry_gas`, and its history sources add the reaction mass sources of the
> flagged components into $G_i$.

### 4.6.5 Type-stable iteration over heterogeneous tuples

The defining performance pattern of the Materials subsystem is the use of
`@generated` functions and recursive tuple decomposition to **unroll** every
loop over the (heterogeneous) component and reaction tuples at compile time. A
naive `for c in material.components` over a heterogeneous tuple would box each
element and allocate; instead, e.g., the cell's heat-capacity sum
$C_i = \sum_j m_{j,i}\,c_{p,j}(T)$ is a recursion over the tuple of the cell's
masses and the tuple of the components' heat-capacity laws
(`src/States/caloric.jl`):

```julia
@inline _heat_capacity_sum(acc, m::Tuple, c::Tuple, T) = _heat_capacity_sum(
    iszero(first(m)) ? acc : acc + first(m) * first(c)(T),
    Base.tail(m),
    Base.tail(c),
    T,
)
@inline _heat_capacity_sum(acc, ::Tuple{}, ::Tuple{}, T) = acc
```

Each call peels the first element of both tuples, whose concrete type is known,
and recurses on the tails; the compiler inlines the recursion into $N_C$
statically dispatched terms with no runtime loop. A component whose mass is
exactly zero is skipped (its term would add an exact zero), and `iszero` of a
dual number reads its partials too, so a derivative is never dropped. The same
idiom — a recursion on tuples or a `@generated` unroll — appears in
`evaluate_reactions!`, the species sources, the partial-pressure and
skeletal-volume sums, the strain rate, and the mixing-rule sums (Chapter 5). The
net effect is that the state mapping and the kernels allocate nothing per
evaluation: all composition-dependent dispatch is resolved at compile time, all
intermediate tuples are statically sized, and there are no dictionary lookups in
hot paths. `test/jacobian/test_allocations.jl` asserts zero bytes for one
residual and one colored-Jacobian evaluation on every configuration it builds.

---

## 4.7 Time- and temperature-parameter evaluation

Several kernels accept "flexible" parameters that may be a constant, a function
of time, or a property function — for example a boundary heat flux $q(t)$, a
ramped convective coefficient, or a temperature-dependent property given as a
property object. Two thin, AD-preserving dispatch helpers normalize these.

**Time-dependent parameters** (`evaluate_param(param, t)`):

$$
\texttt{evaluate\_param}(p, t) =
\begin{cases}
\texttt{convert}(\mathrm{typeof}(t),\,p) & p::\text{Real},\\
p(t) & p::\text{Function}.
\end{cases}
$$

The `convert(typeof(t), p)` for a constant returns a value in the number type
of its argument, so the helper is type-stable for every argument type. In a
residual evaluation $t$ is a plain number: a dual-number time is not supported,
and every integrator policy takes the time derivative of the right-hand side by
finite differences (the docstring of `residual!`, Chapter 15). Where a
time-dependent parameter meets the state, the kernel promotes it explicitly
(the Beer–Lambert kernel adds `zero` of the state's element type to the incident
flux it evaluates), and parameter sensitivities are finite differences of whole
solves (Chapter 17).

**Temperature-dependent parameters** (`evaluate_param_T(param, T)`):

$$
\texttt{evaluate\_param\_T}(p, T) =
\begin{cases}
\texttt{convert}(\mathrm{typeof}(T),\,p) & p::\text{Real},\\
p(T) & p::\text{Function},\\
p(T) & p::\text{AbstractPropertyFunction}.
\end{cases}
$$

with the same type-preserving `convert` for constants. Here the conversion
matters: the temperature carries `ForwardDiff.Dual` numbers in every
colored-Jacobian pass — the surface temperature of the surface energy balance,
for example, at which a `RadiativeBC` evaluates an emissivity given as a law — so
a bare `Float64` constant is promoted to the dual type and the computation stays
in the AD number type end to end. A dedicated fast-path method short-circuits the
most common case:

```julia
@inline evaluate_param_T(param::ConstantProperty, T) = param.value
```

This bypasses the `param(T)` dynamic-dispatch hop and reads the value field
directly — measurable in tight loops where the same constant property is queried
per cell per step.

> **Implementation note.** `evaluate_param`, `evaluate_param_T`, and the
> `ConstantProperty` fast path in `src/Materials/property_functions.jl`.

---

## 4.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

### 4.8.1 Component / species representation

ThermaKin2Ds distinguishes solid, liquid and gas components. Its numerical
states are component masses and temperature (TK §3), as opposed to its
continuum concentration notation. Pyrolysis.jl likewise carries masses, with
temperature or enthalpy as its energy state.

Gpyro's $\bar\rho$ is the bulk density of the condensed mixture in vacuum,
excluding pore-gas mass (GP §§2.1, 4.1); it is not generally
$\sum_{j\in\mathrm{all}}\xi_j$. FDS similarly tracks condensed material
components separately from its CFD gas. Gpyro's reaction parameter $\chi_k$
and solid yield $\mathrm{SF}_k$ are not Pyrolysis.jl's component swelling
weights or lateral conversion variable (§10.11).

### 4.8.2 Heat capacity and the gas-storage convention

The common all-component storage coefficient of Pyrolysis.jl and TK Eq. 17 is

$$
\rho c_p^\mathrm{eff} = \sum_{j=1}^{N_C}\xi_j c_{p,j}(T).
\tag{4.11}
$$

Gpyro's default thermal-equilibrium exchange neglects gas storage and gas
conduction; `FULL_QSG` retains both (GP Eq. 65). FDS's standard solid model
has no pore-gas inventory. These are distinct energy closures (§7.8).

### 4.8.3 Reaction kinetics representation

The restrictions and heat-basis conversions are in §6.9. It is incorrect to
describe all reference heats as scalars: ThermaKin2Ds allows a temperature
function for reaction heat (TK Eq. 1 and §1.1), and FDS constructs
temperature-dependent reaction heats from its enthalpy accounting (FDS
mathematical model, Appendix P). Pyrolysis.jl accepts temperature-dependent
and state-dependent heat functions; matching their datum remains necessary.

### 4.8.4 Property functions and temperature dependence

TK Eq. 1 specifies $p_0+p_1T+p_2T^n$ for density, heat capacity, conductivity
and gas transfer, with constant component emissivity and absorption.
This is representable by a Pyrolysis.jl callable; it is not a documented
ThermaKin table input. GP §§2.1–2.5 specify power-law properties, optional
Gaussian melting heat capacity and a $T^3$ radiative conductivity term.
FDS supports temperature ramps for material properties. Each law can be
transcribed here, but interpolation, extrapolation, units and mixture weights
must be matched explicitly (§5.14).

---

## 4.9 Limitations and open issues

The following are documented limitations of the constitutive layer; several are
revisited where the relevant physics is developed.

1. **Single-reactant reactions only.** $N_\mathrm{in}=1$ is enforced. Multi-step
   schemes must be modeled as *sequences* of single-reactant reactions.
   Multi-reactant (e.g. genuine bimolecular char oxidation) would require
   extending the reactant tuple to 3-ary `(idx, order, mass_stoich)` and a
   revised mass-balance convention (Chapter 6).

2. **Intrinsic-density fallback may be unphysical.** When $\rho_{i,j}$ is not
   supplied it defaults to the bulk $\rho_j$, asserting the pure phase has no
   internal pores. For genuinely microporous solids this under-predicts porosity
   and over-predicts the matrix volume.

3. **Name-based dry-gas identification.** Water vapor is excluded from the
   pyrolysis-progress source by a substring match on the component name
   (`_is_water_vapor_name`). A water species named outside the recognized tokens
   would be (incorrectly) counted as dry pyrolyzate, and a non-water species named
   with a water token would be (incorrectly) excluded.

4. **Kinked and AD-dependent property laws.** A `TableProperty` is differentiated
   segment by segment, so its derivative jumps at every knot and the integrator
   sees a kink in the Jacobian there. `CallableProperty` and `StateDependentProperty`
   are differentiated by propagating `ForwardDiff.Dual` numbers through the user's
   function, which must therefore be generic in its argument types (§4.5).

5. **Heterogeneous-tuple compilation cost.** Preserving the exact component and
   reaction tuple types yields allocation-free, type-stable hot paths, but at the
   cost of recompiling the residual for each distinct material *shape*. Materials
   that differ only in numerical parameters share a compiled type; materials that
   differ in the *types* of their property functions do not.

6. **No temperature-dependent pre-exponential or pressure-dependent kinetics.**
   The Arrhenius property and `Reaction` carry constant $A$ and $E$ (no Kissinger/
   Ozawa $T$-dependent activation energy, no distributed-activation-energy
   tooling, no explicit pressure dependence). See Chapter 6.

7. **Heat-of-reaction sign mapping to reference codes.** Internally $h>0$ is
   endothermic ($Q_\mathrm{rxn}=-\sum_r h_r r_r$), whereas ThermaKin publishes
   $h>0$ as exothermic. Users porting parameters from ThermaKin must flip the
   sign; Gpyro and FDS use the same convention as Pyrolysis.jl, and their heats
   are entered as published (Chapter 2, overload note H1; Chapter 6, §6.9.2).
