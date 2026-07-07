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
relative to the matrix, so their concentration changes only through reaction —
whereas gaseous species are **Eulerian** — they flow through the porous matrix
and therefore carry a flux divergence in addition to a reaction source. The
gas-vs-condensed dispatch appears in the species residual:

$$
\frac{\partial \xi_j}{\partial t} =
\begin{cases}
-\,\nabla\!\cdot\mathbf{J}_j \;+\; \left.\dfrac{\partial \xi_j}{\partial t}\right|_\mathrm{rxn}
   & j \text{ is gas},\\[2.2ex]
\left.\dfrac{\partial \xi_j}{\partial t}\right|_\mathrm{rxn}
   & j \text{ solid or liquid.}
\end{cases}
\tag{4.1}
$$

Because this branch is taken inside the innermost residual loop, over every cell
and every component, it must not incur a runtime type check. Pyrolysis.jl
encodes the phase **at the type level**: there are three concrete subtypes of
`AbstractComponent`, and the phase predicate is a method that dispatches on the
singleton type, so the compiler eliminates the branch entirely at code
generation. The runtime `Phase` enum (`SOLID=1`, `LIQUID=2`, `GAS=3`) exists only
for printing, logging, and comparison; it is *not* consulted in hot paths.

```julia
@enum Phase begin SOLID = 1; LIQUID = 2; GAS = 3 end
```

> **Implementation note.** Phase types and predicates are defined in
> `SolidComponent`, `LiquidComponent`, `GaseousComponent`, and `is_solid` /
> `is_liquid` / `is_gas` / `phase` in `src/Materials/components.jl`. The
> compile-time-unrolled application of Eq. (4.1) is
> `apply_species_residual_unrolled!`, `src/Materials/materials.jl`, which emits
> one statically dispatched `is_gas(material.components[j])` branch per component
> — folded to a constant at codegen.

### 4.1.2 The three component types

**`SolidComponent`** is the matrix-forming condensed phase. It carries six
property-function fields and three scalar (`Float64`) fields:

```julia
struct SolidComponent{Pρ,Pc,Pk,Pλ,Pκ,Pρi} <: AbstractComponent
    name::Symbol
    density::Pρ                # bulk density ρ_j           [kg/m³]
    heat_capacity::Pc          # c_{p,j}(T)                 [J/(kg·K)]
    thermal_conductivity::Pk   # k_j(T)                     [W/(m·K)]
    gas_transfer_coeff::Pλ     # λ_j(T)                     [m²/s]
    permeability::Pκ           # κ_j(T)                     [m²]
    emissivity::Float64        # ε_j                        [-]
    absorption_coeff::Float64  # α_j (mass basis)           [m²/kg]
    swelling_factor::Float64   # γ_j                        [-]
    intrinsic_density::Pρi     # skeletal density ρ_{i,j}   [kg/m³]
end
```

The six property fields are each parameterized by their own type (`Pρ`, `Pc`, …),
so a component built with a constant conductivity and a polynomial heat capacity
is a *distinct concrete type* from one built with two constants. This preserves
type stability and allocation-free evaluation through the heterogeneous component
tuple. The three scalar fields (`ε_j`, `α_j`, `γ_j`) are stored as plain
`Float64` rather than property functions because they are not given temperature
dependence in the current model.

**`LiquidComponent`** has the *identical* field layout to `SolidComponent`
(bulk and intrinsic density, the four property functions, the three scalars). It
is the matrix-forming condensed phase used for moisture / liquid water in
hygroscopic materials. Its convenience constructor differs only in defaults
(see §4.1.5): $\lambda$ defaults to $2\times10^{-5}$ m²·s⁻¹, $\varepsilon$ to
0.95 (condensed liquids are strong emitters — water ≈ 0.96),
$\gamma$ to 0 (volumeless bound moisture, ThermaKin `L SWELLING: 0`), and
$\kappa$ to 0 (a liquid normally fills pore space and blocks gas flow).

**`GaseousComponent`** represents species that flow *through* the matrix:

```julia
struct GaseousComponent{Pρ,Pc,Pk,Pλ} <: AbstractComponent
    name::Symbol
    molar_mass::Float64        # M_j                        [kg/mol]
    density::Pρ                # ρ_g(T), default ideal gas  [kg/m³]
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    emissivity::Float64
    absorption_coeff::Float64
    swelling_factor::Float64
end
```

A gas has **no permeability field and no intrinsic density** — it is not part of
the matrix skeleton. It carries instead a `molar_mass` $M_j$, used by the ideal-
gas density law (§4.3), the ideal-gas pressure closure (Chapter 9), and the
Chapman–Enskog diffusivity estimate (§5.10).

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

Each is `@inline` and resolves at compile time on the singleton type; none branch
on a stored field. The gas's `get_permeability` and `get_intrinsic_density`
return exactly `0.0`, which is the algebraically correct value for the matrix
mixing rules (a gas contributes nothing to permeability or to the solid
skeleton).

> **Implementation note.** Accessors `get_density`, `get_intrinsic_density`,
> `get_heat_capacity`, `get_conductivity`, `get_gas_transfer`,
> `get_permeability`, `has_permeability`, `is_solid/is_liquid/is_gas`, `phase`,
> `molar_mass`, all in `src/Materials/components.jl`.

### 4.1.4 Bulk versus intrinsic (skeletal) density

The single most important constitutive distinction at the component level is
between two densities (overload note A7 of Chapter 2):

- **Bulk (pure-phase) density** $\rho_j$ [kg·m⁻³] — the mass of pure phase $j$ per
  unit of the *bulk volume it would occupy in the mixture*. This is the density
  used by every mixing rule, by the mixture-density formula, and by the volume-
  change kernels (Chapter 10). Field `density`.

- **Intrinsic / skeletal density** $\rho_{i,j}$ [kg·m⁻³] — the density of the
  cell-wall material itself, ignoring any pore volume internal to the component.
  This is used by **one formula only**, the porosity (Chapter 5, §5.3):

  $$
  \phi = \mathrm{clamp}\!\left(1 - \sum_{j\in\{S,L\}} \frac{\xi_j}{\rho_{i,j}},\; 0,\; 1\right).
  \tag{4.2}
  $$

  Field `intrinsic_density`.

When a component is constructed without an explicit intrinsic density, the field
defaults to the bulk density, $\rho_{i,j}=\rho_j$, which recovers the legacy
"pure phase has no internal pores" assumption — every example that predates the
intrinsic-density feature behaves unchanged. The resolution logic enforces that a
caller pass *either* `ρ_intrinsic` *or* a pure-phase porosity $\phi$ (from which
$\rho_{i,j}=\rho_j/(1-\phi)$ is computed, requiring a constant bulk $\rho_j$),
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

$\gamma_j$ enters two places: the mixture-density formula (§4.2.2, where a gas with
$\gamma_j=0$ contributes zero specific volume) and the volumetric strain rate
$\theta = \sum_j \gamma_j \rho_j^{-1}\,\mathrm{d}\xi_j/\mathrm{d}t$ that drives the
ALE mesh motion (Chapter 10). It also weights the **matrix volume fractions**
used for conductivity/emissivity (Chapter 5, §5.1): every component carries
weight $\gamma_j\,\xi_j/\rho_j$ there, so a $\gamma_j = 0$ component (the
canonical gas) carries no weight in the matrix radiative surface.

> **Why the cap $\gamma_j \le 1$ is a model decision, not a missing feature.**
> Construction rejects $\gamma_j > 1$ (`_validate_swelling_factor`) because
> the *same* $\gamma$ weight feeds three mixing rules with volume-fraction
> semantics — the strain kernel, the mixture density, and the $k/\varepsilon$
> matrix volume fractions — and a weight above 1 corrupts all three (volume
> fractions no longer sum consistently). Expansion beyond the reactant's
> volume is already expressible the physical way: give the product a **low
> bulk density** ($V_f/V_0 = \sum_j y_j\gamma_j\rho_A/\rho_j$ — the V-CPL-03
> P1 verification permutation demonstrates ×2 growth), with $\gamma \in
> (0,1]$ tuning how much of that specific volume the matrix keeps. A
> *decoupled* intumescent expansion factor (strain contribution independent
> of the mixing weight) would be a separate modeling design; no current use
> case demands it.

---

## 4.2 Component property fields and their roles

### 4.2.1 Inventory

Each property field of a component, its symbol, its SI units, and where it is
consumed downstream:

| Field | Symbol | Units | Phase | Primary consumer |
|---|---|---|---|---|
| `density` | $\rho_j$ | kg·m⁻³ | all | mixing rules, $\rho$ (§4.2.2), volume change (Ch. 10) |
| `intrinsic_density` | $\rho_{i,j}$ | kg·m⁻³ | S, L | porosity $\phi$ (§4.2.3, Eq. 4.2) |
| `heat_capacity` | $c_{p,j}(T)$ | J·kg⁻¹·K⁻¹ | all | $\rho c_p^\mathrm{eff}$, $\Delta h$, advection source (Ch. 5, 7) |
| `thermal_conductivity` | $k_j(T)$ | W·m⁻¹·K⁻¹ | all | $k_\mathrm{eff}$ (Ch. 5), Fourier flux (Ch. 7) |
| `gas_transfer_coeff` | $\lambda_j(T)$ | m²·s⁻¹ | all | $\lambda_\mathrm{eff}$ (Ch. 5), gas diffusion (Ch. 9) |
| `permeability` | $\kappa_j(T)$ | m² | S, L | $\kappa_\mathrm{eff}$ (Ch. 5), Darcy velocity (Ch. 9) |
| `emissivity` | $\varepsilon_j$ | – | all | $\varepsilon_\mathrm{eff}$ (Ch. 5), radiative BC (Ch. 12) |
| `absorption_coeff` | $\alpha_j$ | m²·kg⁻¹ | all | $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ (Ch. 5, 8) |
| `swelling_factor` | $\gamma_j$ | – | all | $\rho$, strain rate $\theta$ (Ch. 10) |
| `molar_mass` | $M_j$ | kg·mol⁻¹ | G | $\rho_g$ (§4.3), $P$ (Ch. 9), Chapman–Enskog estimate (§5.10) |

Note that $\alpha_j$ is a **mass-basis** absorption coefficient [m²·kg⁻¹]; it is
multiplied by concentration to give the volumetric (extinction) absorption
coefficient $\alpha_\mathrm{eff}=\sum_j\xi_j\alpha_j$ [m⁻¹] used in Beer–Lambert
and P1 radiation (Chapter 8). Even a solid carries a nonzero $\lambda_j$ default
in practice (e.g. $2\times10^{-5}$ m²·s⁻¹): although solids do not themselves
diffuse, their $\lambda_j$ contributes to the *mixture* gas-transfer coefficient
through the volume-weighted mixing rule, matching the ThermaKin convention.

### 4.2.2 Mixture density with per-component swelling

The mixture density combines per-component bulk densities through specific
volumes (volume per unit mass), every component weighted by its own swelling
factor (the universal γ semantics shared with θ and the mixing-rule volume
fractions, §5.1.3):

$$
\frac{1}{\rho}
= \sum_{j} \gamma_j\,\frac{Y_j}{\rho_j},
\qquad Y_j = \frac{\xi_j}{\rho_\mathrm{tot}},\quad
\rho_\mathrm{tot}=\sum_k \xi_k.
\tag{4.3}
$$

Here $Y_j$ is the mass fraction. A non-swelling gas ($\gamma_j=0$) contributes
**zero** specific volume — its mass is "free" in the pore space and does not
change the bulk volume — whereas a swelling gas ($\gamma_j=1$) contributes its
full specific volume like a solid. This is the discrete realization of the
swelling semantics of §4.1.5.

> **Important — correction to a common gloss.** The swelling weight in Eq. (4.3)
> is applied **per component**, with each component's *own* `swelling_factor` —
> not an averaged $\gamma_\mathrm{avg}$ multiplying the gas sum, and not gases
> only. The source evaluates `γ_j · Y_j/ρ_j` term by term for every component
> (a γ = 0 term contributes nothing).

> **Implementation note.** `mixture_density` and the inline, allocation-free
> kernel `_specific_vol_with_swelling_direct`,
> `src/Properties/mixture_rules.jl`. The
> inline form computes $Y_j=\xi_j\,/\rho_\mathrm{tot}$ on the fly to avoid
> allocating an intermediate mass-fraction tuple. Degeneracy guards: if
> $\rho_\mathrm{tot}<10^{-20}$ or the accumulated specific volume
> $<10^{-20}$, the function returns `0.0`.

### 4.2.3 Porosity and its single use

Porosity $\phi$ (Eq. 4.2) is the void fraction of the bulk *not* occupied by the
solid/liquid skeletal frame. It is computed from **intrinsic** densities and is
the canonical definition $\phi = \mathrm{clamp}(1-\sum_{j\in\{S,L\}}\xi_j/\rho_{i,j},0,1)$.
Its derivatives are exactly (in the unclamped interior; at an active clamp the
derivative is 0, and the AD path used in production returns exactly that)

$$
\frac{\partial \phi}{\partial \xi_j} =
\begin{cases}
-\,1/\rho_{i,j} & j\in\{S,L\},\\
0 & j\in G.
\end{cases}
\tag{4.4}
$$

A gas does not enter $\phi$ because it occupies the void that the solid frame
leaves; adding gas mass does not displace skeletal material. Porosity feeds (i)
the ideal-gas pressure closure $P=N(\xi,T)/\phi$ (Chapter 9) and (ii) the
Carman–Kozeny permeability closure (Chapter 5, §5.7). A *second*, volume-fraction
porosity $\phi=\sum_{j\in G} v_j$ also exists in the codebase
(`_porosity_sum`/`_porosity_sum_gen`) but the **intrinsic-density form Eq. (4.2)
is canonical** and is what `porosity` and `update_properties!` actually use
(overload note A8 of Chapter 2).

> **Implementation note.** Canonical porosity: `porosity` →
> `_porosity_void_frac`, in `src/Properties/mixture_rules.jl` (derivatives are
> taken by forward-mode AD through `porosity` itself; §5.3.3). In
> `_porosity_void_frac` a component
> contributes $\xi_j/\rho_{i,j}$ only when `is_solid(c) || is_liquid(c)` *and*
> $\rho_{i,j}>0$. Full development is in Chapter 5.

---

## 4.3 Ideal-gas density and the gaseous-component default

By default a `GaseousComponent` constructed via the convenience constructor is
given a density obeying the ideal-gas law at the **reference pressure**:

$$
\rho_g(T) = \frac{P_\mathrm{ref}\,M_j}{R_g\,T},
\tag{4.5}
$$

with $P_\mathrm{ref}=101325$ Pa, $R_g=8.31446261815324$ J·mol⁻¹·K⁻¹, and $M_j$
the gas molar mass. The constructor wraps this in a `CallableProperty`:

```julia
ρ_func = CallableProperty(T -> P_REF * M_float / (R_GAS * T))
```

Two modeling points deserve emphasis:

1. **The density used for mixing/volume is at $P_\mathrm{ref}$, not the local
   pressure $P$.** The component's `density(T)` is consulted for volume fractions,
   specific volumes, and the mixture density. The *local* pressure-dependent gas
   density (which matters for the pore-pressure closure $P=N/\phi$) is handled
   separately by the ideal-gas formula in the Darcy machinery (Chapter 9), not by
   the component's `density` field. The component-level density at $P_\mathrm{ref}$
   is what enters the **specific volume** in mixing rules. In the
   swelling-weighted bases (mixture density, conductivity/emissivity weights)
   the gas's contribution is usually gated to zero by $\gamma_j=0$ (Eq. 4.3);
   in the **all-phase** basis used for λ mixing (§5.6) the gas volume is
   present regardless of γ and is over-weighted by ~$P/P_\mathrm{ref}$ in a
   pressurized cell — acceptable at the ≤ atm-scale overpressures of the
   validated cases. The pressure closure $P=N/\phi$ and the porosity are
   unaffected (they never use this density).

2. **The user may override the gas density** by passing any property function
   (constant, polynomial, table, callable) — Eq. (4.5) is merely the default
   produced when one constructs a gas with `GaseousComponent(:name; M=…, c=…, k=…,
   λ=…)`.

> **Implementation note.** `GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0,
> γ=0.0)` constructor, `src/Materials/components.jl`. Defaults for a gas are
> $\varepsilon=0$, $\alpha=0$, $\gamma=0$ (a gas escapes without swelling and is
> not part of the radiating surface).

---

## 4.4 Property-function types

A *property function* is a callable object $p$ such that $p(T)$ (or, for
state-dependent properties, $p(T,\xi)$) returns the property value at the given
state. All seven concrete types subtype `AbstractPropertyFunction` and implement
both the `(p)(T)` and `(p)(T, ξ)` call signatures (the latter ignoring $\xi$ for
all but `StateDependentProperty`). This uniformity lets a single component field
hold any property law, and lets the kinetics evaluator pass $\xi$ to a heat-of-
reaction property without knowing whether it is composition-dependent.

### 4.4.1 `ConstantProperty`

$$
p(T) = \texttt{value}, \qquad \frac{\mathrm{d}p}{\mathrm{d}T}=0.
$$

The workhorse for properties given a single number. It has a *fast path* in the
parameter-evaluation interface (§4.8) that bypasses the dynamic-dispatch hop.

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
recomputing powers. The analytical derivative is computed by an unrolled loop
(no allocation, no `evalpoly` of a derived coefficient tuple):

$$
\frac{\mathrm{d}p}{\mathrm{d}T} = \sum_{m=1}^{N-1} m\,c_m\,T^{m-1}
= c_1 + 2c_2 T + 3c_3 T^2 + \cdots .
\tag{4.6}
$$

For $N=1$ (a degenerate constant polynomial) the derivative returns `0.0`.

### 4.4.4 `ArrheniusProperty`

$$
p(T) = A\,\exp\!\left(-\frac{E}{R_g\,T}\right),
\qquad
\frac{\mathrm{d}p}{\mathrm{d}T} = p(T)\,\frac{E}{R_g\,T^2}.
\tag{4.7}
$$

Fields `A` (pre-exponential), `E` (activation energy [J·mol⁻¹]), and `T_ref` (a
reference temperature, defaulting to $T_\mathrm{ref}=298.15$ K, carried for
normalization purposes but not used in the value itself). The closed-form
derivative in Eq. (4.7) is obtained by differentiating
$A\exp(-E/(R_gT))$: $\frac{\mathrm{d}}{\mathrm{d}T}[-E/(R_gT)] = E/(R_gT^2)$, hence
$p'(T)=p(T)\cdot E/(R_gT^2)$. This is the *only* property type with a nonzero
closed-form $T$-derivative beyond the polynomial/linear forms, and it shares its
functional form with the reaction-rate constant (§4.5, Chapter 6).

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
search on small tables) and **binary for $N>8$** (logarithmic). The derivative is
a *central finite difference* with a fixed step $\delta T = 0.1$ K,
$\mathrm{d}p/\mathrm{d}T \approx [p(T_q+\delta T) - p(T_q-\delta T)]/(2\delta T)$.
On a smooth segment this recovers the segment slope exactly; near a knot it blends
the two adjacent slopes.

### 4.4.6 `CallableProperty`

Wraps an arbitrary callable $f(T)$. Its value is $p(T)=f(T)$; its derivative is the
same central finite difference with $\delta T = 0.1$ K. This is the type produced
by `to_property(f)` for a function `f` and is how the ideal-gas density default
(Eq. 4.5) is stored.

### 4.4.7 `StateDependentProperty`

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

Its $T$-derivative is a central finite difference at fixed $\xi$,
$\partial p/\partial T \approx [p(T+\delta T,\xi)-p(T-\delta T,\xi)]/(2\delta T)$
with $\delta T = 0.1$ K. Its concentration derivative (§4.5.2) uses an *adaptive*
step.

### 4.4.8 `to_property` conversions

The convenience converter normalizes user inputs to property functions, so that
component constructors and reaction constructors accept either property objects
or bare values:

| Input | Result |
|---|---|
| an `AbstractPropertyFunction` | passes through unchanged |
| a `Real` number | `ConstantProperty(Float64(x))` |
| a `Tuple{Real,Real}` $(a,b)$ | `LinearProperty(a, b)` |
| a `Function` | `CallableProperty(f)` |

> **Implementation note.** All property types, their value/derivative methods,
> `linear_interpolate`, and `to_property` are in
> `src/Materials/property_functions.jl`.

---

## 4.5 Property derivatives and their role in the Jacobian

Implicit time integration (Chapter 15) requires the Jacobian
$\partial(\text{RHS})/\partial u$. Because effective properties depend on $T$ and
$\xi$, and the RHS depends on effective properties, the per-component derivatives
$\mathrm{d}p_j/\mathrm{d}T$ and $\partial p_j/\partial \xi_k$ propagate through the
chain rule into the structured Jacobian (Chapter 5, §5.15 for the mixing-rule
chain; Chapter 15 for assembly). Two derivative families are provided.

### 4.5.1 Temperature derivative $\mathrm{d}p/\mathrm{d}T$

The `derivative(p, T)` method dispatches on the property type:

| Type | $\mathrm{d}p/\mathrm{d}T$ | Method |
|---|---|---|
| `ConstantProperty` | $0$ | analytic |
| `LinearProperty` | $b$ | analytic |
| `PolynomialProperty` | $\sum_{m\ge 1} m c_m T^{m-1}$ | analytic, unrolled (Eq. 4.6) |
| `ArrheniusProperty` | $p(T)\,E/(R_g T^2)$ | analytic (Eq. 4.7) |
| `TableProperty` | central FD, $\delta T=0.1$ K | numerical |
| `CallableProperty` | central FD, $\delta T=0.1$ K | numerical |
| `StateDependentProperty` | central FD in $T$, $\delta T=0.1$ K | numerical |

The analytic derivatives are exact to machine precision; the finite-difference
derivatives carry an $O(\delta T^2)$ truncation error from the central scheme,
which at $\delta T=0.1$ K is negligible for the smooth thermophysical laws this
machinery serves. (The Jacobian itself may also be obtained by forward-mode
automatic differentiation, in which case these hand-coded derivatives are not the
source of truth; see §4.6 and Chapter 15. The hand-coded forms exist for the
*structured* Jacobian backend, which assembles per-operator contributions
analytically for speed.)

### 4.5.2 Concentration derivative $\partial p/\partial \xi_j$

For every property type *except* `StateDependentProperty`, the value does not
depend on $\xi$, so $\partial p/\partial \xi_j = 0$ exactly. For a
`StateDependentProperty` the derivative is an adaptive central finite difference,

$$
\frac{\partial p}{\partial \xi_j} \approx
\frac{p(T,\xi^{+}) - p(T,\xi^{-})}{2\,\delta\xi_j},
\qquad
\delta\xi_j = \max\!\bigl(|\xi_j|\cdot10^{-4},\;10^{-6}\bigr),
\tag{4.9}
$$

where $\xi^{\pm}$ are the concentration tuple with the $j$-th entry perturbed by
$\pm\delta\xi_j$ (using `Base.setindex` on the immutable tuple to avoid heap
allocation). The adaptive step is *relative* for non-small concentrations (the
$10^{-4}$ factor gives good FD accuracy) and switches to an *absolute floor*
$10^{-6}$ for vanishing $\xi_j$ to avoid catastrophic cancellation. At a boundary
where $\xi_j - \delta\xi_j < 0$ (which would push the concentration negative,
potentially out of the physical domain of $f$) the scheme falls back to a
**forward difference** $[p(T,\xi^{+})-p(T,\xi)]/\delta\xi_j$ to keep the
evaluation in the non-negative region; in the central branch the negative
perturbation is non-negative by construction.

> **Implementation note.** `derivative(p, T)` and
> `derivative_ξ(p, T, ξ, j)` in `src/Materials/property_functions.jl`. The default
> `derivative_ξ(::AbstractPropertyFunction, …) = 0.0` covers all $T$-only
> property types; only the `StateDependentProperty` method computes Eq. (4.9).

### 4.5.3 The reaction-rate connection

The Arrhenius derivative (Eq. 4.7) is identical in form to the
temperature-derivative of the reaction-rate constant
$k(T)=A\exp(-E/(R_g T))$, since $\mathrm{d}k/\mathrm{d}T = k\,E/(R_g T^2)$. This is
why the reaction-rate Jacobian energy block (Chapter 15) reuses the same
expression. The full reaction-rate derivative including the depletion limiter and
temperature gates is developed in Chapter 6.

---

## 4.6 The sensible-enthalpy difference (midpoint rule)

Several energy terms — the surface transpiration enthalpy (Chapter 12), the
gas-generation enthalpy sink $S_\mathrm{gen}$ (Chapter 7), and the energy ledgers
(Chapter 16) — require the **sensible-enthalpy difference** of a single component
between two temperatures, $\Delta h = \int_{T_\mathrm{lo}}^{T_\mathrm{hi}}
c_p(\tau)\,\mathrm{d}\tau$. Pyrolysis.jl evaluates this with the **midpoint rule**:

$$
\Delta h(T_\mathrm{hi}, T_\mathrm{lo})
= c_p\!\left(\frac{T_\mathrm{hi}+T_\mathrm{lo}}{2}\right)\,
  \bigl(T_\mathrm{hi}-T_\mathrm{lo}\bigr).
\tag{4.10}
$$

### 4.6.1 Why the midpoint rule

Three reasons justify Eq. (4.10) over alternatives such as trapezoidal,
endpoint, or analytic integration of a polynomial $c_p$:

1. **Second-order accuracy with a single evaluation.** The midpoint rule has
   truncation error $\tfrac{1}{24}c_p''(T_m)\,(\Delta T)^3$, i.e. it is exact for
   $c_p$ linear in $T$ and $O(\Delta T^2)$ for general $c_p$ — the same order as
   the trapezoidal rule but with **one** $c_p$ call instead of two. In a residual
   evaluated tens of thousands of times this halves the cost.

2. **Consistency with the cell-centered $c_p$.** The volumetric gas-advection
   source uses a *single* cell-centered evaluation $c_{p,g}(T_i)$ (Chapter 7,
   $S_\mathrm{conv} = -\sum_g c_{p,g}(T_i)\bar J_g (\partial T/\partial z)$).
   Using the midpoint rule for $\Delta h$ keeps the surface-transpiration BC at
   the *same discretization order* as the bulk gas-energy transport, so the two
   do not introduce inconsistent enthalpy accounting at the surface.

3. **AD-cleanliness.** Eq. (4.10) is a smooth composition of $c_p(\cdot)$ and
   multiplication, so both $T_\mathrm{hi}$ and $T_\mathrm{lo}$ may be
   `ForwardDiff.Dual`. ForwardDiff returns the derivative of *the midpoint-rule
   expression actually evaluated*, namely

   $$
   \frac{\partial(\Delta h)}{\partial T_\mathrm{hi}}
   = c_p(T_m) + \tfrac{1}{2}(T_\mathrm{hi}-T_\mathrm{lo})\,c_p'(T_m),
   \quad T_m=\tfrac{T_\mathrm{hi}+T_\mathrm{lo}}{2},
   $$

   which differs from the *true* integral derivative $c_p(T_\mathrm{hi})$ at
   $O(\Delta T)$ but is **exactly consistent** with the value the function
   returns. This consistency (value and derivative come from the same quadrature)
   is precisely what makes the Newton surface-energy solve and the sensitivity
   machinery well-behaved; an inconsistent value/derivative pair would corrupt
   the Newton convergence.

> **Implementation note.** `delta_enthalpy(comp, T_hi, T_lo)`,
> `src/Materials/components.jl`, which computes
> `comp.heat_capacity((T_hi+T_lo)*0.5) * (T_hi - T_lo)`. It dispatches on any
> `AbstractComponent` (uses only the `heat_capacity` field).

---

## 4.7 Material assembly

### 4.7.1 The `Material` type

A `Material` is the static, fully type-specialized aggregate that every physics
kernel reads. It is **immutable** and **heterogeneously typed** so the compiler
can specialize every loop over its components and reactions:

```julia
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL} <: AbstractMaterial
    name::Symbol
    components::C          # heterogeneous component tuple
    reactions::RT          # heterogeneous reaction tuple
    mixing::NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}
    lateral_shrinkage_law::LSL   # A(t)/A₀ = law(χ̄), or nothing
    depletion_limiter::DepletionLimiter   # tanh(ξ/threshold) rate roll-off (§6.3)
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
(default `threshold = 1.0` kg·m⁻³, `enabled = true`), set via the constructor's
`depletion_limiter` keyword and read by `evaluate_reactions!` for every
reaction-rate evaluation (§6.3).

> **Implementation note.** `Material`, its inner/outer constructors,
> `_construct_material`, `MixingRule`, `MixingSpec`, `MaterialMixing`, in
> `src/Materials/materials.jl`. Counts via `n_components`, `n_reactions`,
> `n_solids`, `n_liquids`, `n_gases`.

### 4.7.2 `PendingReaction` and symbolic resolution

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

### 4.7.3 The single-reactant mass-balance convention

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
> (tol $10^{-10}$, returns `false` + `@warn` on failure), `mass_balance_error`,
> and `check_material_mass_balance` in `src/Materials/reactions.jl`.

### 4.7.4 `validate_material`

After assembly, `validate_material(m)` performs whole-material consistency checks:

1. at least one component exists ($N_C\ge 1$, else `error`);
2. every reactant and product index in every reaction lies in $1\!:\!N_C$ (else
   `error`);
3. every reaction is mass-balanced (else `@warn` via `verify_stoichiometry`).

It returns `true` on success. This is a *build-time* invariant: by the time the
residual runs, indices are guaranteed valid and stoichiometry is (warned-)
balanced.

> **Implementation note.** `validate_material(m::Material)`,
> `src/Materials/materials.jl`.

### 4.7.5 Dry-pyrolysis-gas identification

For variable-cross-section (lateral shrinkage, Chapter 10) the pyrolysis-progress
state $\chi$ counts only **dry** pyrolysis gas — water vapor is excluded because
moisture evaporation should not drive radial shrinkage. The set of dry-gas
component indices is computed at compile time by a recursive tuple walk that
includes a component index iff it is a gas **and** its name is not a water-vapor
name. The water-vapor test is a substring/equality heuristic on the component
name (`"H2O"`, `"water"`/`"Water"`, or exactly `"h2o"` case-insensitively):

```
is_gas(c) && !_is_water_vapor_name(c.name)  ⇒ include index
```

This name-based identification is a documented heuristic (open issue §4.9); a
gas whose name does not contain a recognized water token is treated as dry
pyrolyzate regardless of its actual chemistry.

> **Implementation note.** `dry_pyrolysis_gas_indices` →
> `_collect_dry_pyrolysis_gas_indices` and the heuristic `_is_water_vapor_name`,
> `src/Materials/materials.jl`. The result is a statically known tuple of indices,
> constant-folded by the compiler.

### 4.7.6 Type-stable iteration over heterogeneous tuples

The defining performance pattern of the Materials subsystem is the use of
`@generated` functions and recursive tuple decomposition to **unroll** every
loop over the (heterogeneous) component and reaction tuples at compile time. A
naive `for c in material.components` over a heterogeneous tuple would box each
element and allocate; instead, e.g.,

```julia
@generated function foreach_component(f::F, material::Material{NC}) where {F,NC}
    return Expr(:block, [:(f(material, Val{$j}())) for j in 1:NC]...)
end
```

emits $N_C$ statically dispatched calls with no runtime loop. The same idiom
appears in `evaluate_reactions!`, the species-source and gas-flux kernels, the
mixing-rule sums (Chapter 5), and the porosity sums. The net effect is that the
entire property-evaluation loop in `update_properties!` is **allocation-free**
once the cache is sized (Chapter 5, §5.13): all composition-dependent dispatch is
generated, all intermediate tuples are statically sized, and there are no
dictionary lookups in hot paths.

---

## 4.8 Time- and temperature-parameter evaluation

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

The `convert(typeof(t), p)` for a constant is the crucial detail: when $t$ is a
`ForwardDiff.Dual` (during a sensitivity solve, Chapter 17), a bare `Float64`
constant is promoted to the Dual type so the computation stays in the AD number
type end to end. Returning the raw `Float64` would break the Dual propagation.

**Temperature-dependent parameters** (`evaluate_param_T(param, T)`):

$$
\texttt{evaluate\_param\_T}(p, T) =
\begin{cases}
\texttt{convert}(\mathrm{typeof}(T),\,p) & p::\text{Real},\\
p(T) & p::\text{Function},\\
p(T) & p::\text{AbstractPropertyFunction}.
\end{cases}
$$

with the same Dual-preserving `convert` for constants. A dedicated fast-path
method short-circuits the most common case:

```julia
@inline evaluate_param_T(param::ConstantProperty, T) = param.value
```

This bypasses the `param(T)` dynamic-dispatch hop and reads the value field
directly — measurable in tight loops where the same constant property is queried
per cell per step.

> **Implementation note.** `evaluate_param`, `evaluate_param_T`, and the
> `ConstantProperty` fast path in `src/Materials/property_functions.jl`.

---

## 4.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

Pyrolysis.jl's material model is a deliberate synthesis of the three reference
condensed-phase pyrolysis codes. We map their notation to ours on first use.

### 4.9.1 Component / species representation

**ThermaKin2Ds** represents a material as a set of *components* with mass
concentrations as the primary state — the convention Pyrolysis.jl adopts wholesale
($\xi_j$ in kg·m⁻³ is the species state, not mass fraction). ThermaKin distinguishes
condensed from gaseous components and applies the *solid-volume-fraction*
convention to conductivity and emissivity — exactly Pyrolysis.jl's
`solid_volume_fractions` (Chapter 5, §5.1), which excludes non-swelling gases
(ThermaKin: "when $\gamma=0$, gaseous components do not contribute to the thermal
conductivity"). ThermaKin's stoichiometric coefficient $\theta_i^{\,j}$ maps to our
$\nu_{i,j}$ (overload note A2).

**Gpyro** uses a control-volume formulation with $\bar\rho\,\Delta z$ (mass per
unit area) and mass fractions $Y_i$ of *condensed-phase species*; its gas species
are tracked separately. Gpyro's volume-fraction property averaging $X_\alpha$
maps to our $v_j$, and its power-law temperature dependence of properties is
representable here as a `PolynomialProperty` or `CallableProperty`. Gpyro's
solid-fraction $\mathrm{SF}_k$ / conversion $\chi_k$ mechanism for selecting
non-charring / charring / intumescent / shrinkage behavior corresponds to our
swelling-factor $\gamma_j$ plus the lateral-shrinkage law (Chapter 10).

**FDS** (solid phase) uses a layered material model with per-material density,
$c_p$, $k$, and emissivity, plus Arrhenius reaction parameters; its
pre-exponential $Z$ maps to our $A_i$ (overload note A3) and its solid absorption
coefficient $\kappa_s$ maps to our volumetric $\alpha$ (overload note A5). FDS's
material density-ratio swell/shrink rule corresponds to our $\gamma_j$-weighted
volume change (Chapter 10).

### 4.9.2 Heat capacity and the gas-storage convention

The most consequential modeling choice — and a genuine departure from one of the
three — concerns gas sensible storage. ThermaKin's energy equation (its Eq. 17)
sums the volumetric heat capacity over **all** components, gas included:
$\rho c_p^{\,\text{ThermaKin}} = \sum_{j\in\text{all}} \xi_j c_{p,j}$. Pyrolysis.jl
instead uses the **matrix-only** effective heat capacity in the temperature-
storage term (Formulation A, Chapter 3):

$$
\rho c_p^\mathrm{eff} = \sum_{j\in\{S,L\}} \xi_j\, c_{p,j}(T),
\tag{4.11}
$$

dropping the in-pore gas storage $\sum_{j\in G}(\rho c)_j\,\partial T/\partial t$
and carrying the gas sensible energy instead through the advective source
$S_\mathrm{conv}$ (the quasi-steady-gas approximation; Chapter 5, §5.4 and
Chapter 7). This is the **FDS/Gpyro-style** convention, *not* ThermaKin's. The
dropped term is $\sim0.2\%$ of the matrix term at 1 atm and scales with pressure;
the diagnostics tracker (Chapter 16) quantifies it as the *physical-form
residual* via the total heat capacity
$\rho c_p^\mathrm{total} = \sum_{j\in\text{all}} \xi_j c_{p,j}(T)$, which is
computed for diagnostics only and never enters the RHS.

> **Implementation note.** `effective_heat_capacity` →`_sum_heat_capacity`
> (gas-excluded, via the compile-time `is_gas` branch) and `total_heat_capacity`
> →`_sum_heat_capacity_total` (all components), `src/Properties/effective_properties.jl`.

### 4.9.3 Reaction kinetics representation

ThermaKin supports bimolecular rates $\propto \xi_{\text{COMP1}}\xi_{\text{COMP2}}$;
Gpyro uses an nth-order-in-conversion rate with optional O₂ sensitivity
$g(Y_{O_2})$ and a critical-temperature Heaviside; FDS uses
$\rho^n\cdot Y_{O_2}^{n}\cdot$ a temperature factor with kinetic-parameter recovery
from a TGA reference temperature and peak rate. Pyrolysis.jl restricts to
**single-reactant** mass-action kinetics in mass-concentration basis with smooth
tanh temperature gates and a smooth depletion limiter (Chapter 6). The
representation of heats of reaction as property functions — including the
state-dependent moisture-sorption heat via `StateDependentProperty` — generalizes
all three codes' scalar heats. The full kinetic comparison is in Chapter 6.

### 4.9.4 Property functions and temperature dependence

ThermaKin and Gpyro both permit temperature-dependent properties (tabulated or
power-law). Pyrolysis.jl's property-function abstraction subsumes both: a Gpyro
power law $p(T)=p_0(T/T_0)^m$ is a `CallableProperty`; a ThermaKin table is a
`TableProperty`; FDS ramp-style piecewise data is a `TableProperty` with linear
interpolation and flat extrapolation. The novel addition is the uniform
$(T,\xi)$ call signature, which allows composition-dependent properties (heat of
sorption) without a special-case code path.

---

## 4.10 Limitations and open issues

The following are documented limitations of the constitutive layer; several are
revisited where the relevant physics is developed.

1. **Single-reactant reactions only.** $N_\mathrm{in}=1$ is enforced. Multi-step
   schemes must be modeled as *sequences* of single-reactant reactions.
   Multi-reactant (e.g. genuine bimolecular char oxidation) would require
   extending the reactant tuple to 3-ary `(idx, order, mass_stoich)` and a
   revised mass-balance convention (Chapter 6).

2. **Two coexisting porosity definitions.** The intrinsic-density form (Eq. 4.2)
   is canonical and is what the pressure and permeability closures use, but a
   volume-fraction porosity also exists in the code. They differ for materials
   with internal pore structure; the choice is not user-selectable.

3. **Intrinsic-density fallback may be unphysical.** When $\rho_{i,j}$ is not
   supplied it defaults to the bulk $\rho_j$, asserting the pure phase has no
   internal pores. For genuinely microporous solids this under-predicts porosity
   and over-predicts the matrix volume.

4. **Name-based dry-gas identification.** Water vapor is excluded from the
   pyrolysis-progress source by a substring match on the component name
   (`_is_water_vapor_name`). A water species named outside the recognized tokens
   would be (incorrectly) counted as dry pyrolyzate, and a non-water species named
   with a water token would be (incorrectly) excluded.

5. **Finite-difference derivatives for table/callable/state-dependent properties.**
   `TableProperty`, `CallableProperty`, and the $T$- and $\xi$-derivatives of
   `StateDependentProperty` use fixed-step ($\delta T=0.1$ K) or adaptive-step
   (Eq. 4.9) finite differences. These carry truncation/round-off error and can
   degrade the structured-Jacobian accuracy relative to the analytic property
   types. (The AD-based Jacobian backend avoids this by differentiating the
   actual evaluation.)

6. **Heterogeneous-tuple compilation cost.** Preserving the exact component and
   reaction tuple types yields allocation-free, type-stable hot paths, but at the
   cost of recompiling the residual for each distinct material *shape*. Materials
   that differ only in numerical parameters share a compiled type; materials that
   differ in the *types* of their property functions do not.

7. **No temperature-dependent pre-exponential or pressure-dependent kinetics.**
   The Arrhenius property and `Reaction` carry constant $A$ and $E$ (no Kissinger/
   Ozawa $T$-dependent activation energy, no distributed-activation-energy
   tooling, no explicit pressure dependence). See Chapter 6.

8. **Heat-of-reaction sign mapping to reference codes.** Internally $h>0$ is
   endothermic ($Q_\mathrm{rxn}=-\sum_r h_r r_r$), whereas ThermaKin/Gpyro publish
   $h>0$ as exothermic. Users porting parameters from those codes must flip the
   sign (Chapter 2, overload note H1; Chapter 6).
