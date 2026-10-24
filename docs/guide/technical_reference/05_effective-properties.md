# 5. Effective (Mixture) Properties and Mixing Rules

This chapter derives, in full, how Pyrolysis.jl collapses the per-component
thermophysical and transport properties of a multi-phase pyrolyzing solid into
the *single* effective property each governing equation needs. It covers the
averaging bases (mass fractions and volume fractions, the two flavours of the
latter, and the blend between them in an emptying cell), the swelling-aware
mixture density, the intrinsic-density pore fraction, the mixture heat capacity
and the sensible enthalpies the energy row reads, the four conductivity mixing
rules (parallel / series / weighted / Bruggeman) together with the Bruggeman
Newton solve, the effective gas-transfer and permeability closures (component
mixing versus Carman–Kozeny), the radiative properties (emissivity, surface
absorptivity and absorption), the gas viscosity (Sutherland) and Chapman–Enskog
diffusivity, the absence of per-species transport scaling, the guarded
distance-weighted harmonic face mean, the cell properties the state mapping
evaluates, and how the colored Jacobian obtains their derivatives. Every
equation, constant, guard, threshold, and sign below follows
`src/Properties/mixture_rules.jl`, `src/Properties/effective_properties.jl`,
`src/States/derived.jl`, `src/States/caloric.jl`, `src/States/guards.jl`,
`src/Materials/components.jl`, and `src/Materials/materials.jl`; where the
implementation departs from the textbook form (most importantly which basis each
property uses, the harmonic-mean face average, and which mixing rules are
rejected for which channel), **the code is authoritative** and the deviation is
called out explicitly.

Throughout this chapter `j` is the **component (species) index** running over the
`N_C` components of a `Material`, and `i` is the **cell index** running over the
`n` finite volumes (overload resolution I1, §2). The primary species variable is
the **mass concentration** `ξ_j` (kg·m⁻³), not the mass fraction (§3, Governing
Equations); the integrated state is the cell mass `m_{j,i} = ξ_{j,i} V_i`. Two
density symbols are used with care: the **bulk** (pure-phase) density `ρ_j`
drives all mixing rules and the mixture density, while the **intrinsic**
(skeletal, cell-wall) density `ρ_{i,j}` appears *only* in the pore-fraction
formula (overload resolution A7, §2). The three binding sign conventions of §2
are assumed without restatement.

The property functions below are written in `(ξ, T)`. The residual evaluates
every one of them once per cell, in the state mapping, at the positive-part
concentrations `ξ⁺_{j,i}` and the trial temperature `T_trial` (§5.13); both
equal `ξ_{j,i}` and `T_i` on every physical state the tolerances can
distinguish.

---

## 5.1 Averaging bases: mass fractions and volume fractions

A homogenized property must be a weighted combination of its component values,
and the *correct weight* depends on the physics. Pyrolysis.jl uses three distinct
weighting bases, each computed from the cell's mass-concentration tuple
`ξ = (ξ_1,…,ξ_{N_C})` and (for the volume bases) the temperature `T`.

### 5.1.1 Total density and mass fractions

The **total (bulk) density** of a cell is simply the sum of mass concentrations,

$$
\rho_{\text{tot}}(\boldsymbol\xi) \;=\; \sum_{j=1}^{N_C} \xi_j ,
\tag{5.1}
$$

This is the density in *bulk-density space* — the mass per unit bulk volume — before any
swelling correction (§5.2). It is the natural denominator for **per-unit-mass**
properties, because such a property averaged over the mixture must weight each
component by its share of the mass.

The **mass fraction** of component `j` is

$$
Y_j \;=\; \frac{\xi_j}{\rho_{\text{tot}}}, \qquad \sum_j Y_j = 1 ,
\tag{5.2}
$$

Mass fractions are the conceptual weights for the swelling mixture density
(§5.2), where `Y_j` is folded in inline (with an empty-cell guard
`ρ_tot < 1e-20`) so that no mass-fraction tuple is built; the heat-capacity
sum of the energy row weights by the cell masses `m_{j,i}` directly (§5.4).
They are *not* the right weights for transport properties, which act over
volumes.

### 5.1.2 All-phase volume fractions

A transport property such as thermal conductivity is governed by *how much
volume* each phase occupies along the heat-flow path, not by how much mass it
carries. The **bulk volume fraction** of component `j` is

$$
v_j \;=\; \frac{\xi_j/\rho_j(T)}{\displaystyle\sum_{k=1}^{N_C}\xi_k/\rho_k(T)},
\qquad \sum_j v_j = 1 ,
\tag{5.3}
$$

where `ρ_j(T)` is the component's **bulk** density evaluated through its property
function (constant, polynomial, table, or callable — §4). The quantity
`ξ_j/ρ_j` is the volume occupied by component `j` per unit bulk volume, so (5.3)
normalizes those volumes to a partition of unity. This is `volume_fractions(ξ, T,
material)` (`mixture_rules.jl`).

The function builds the denominator `v_total = Σ_k ξ_k/ρ_k` via
`_sum_specific_vol` and, if `v_total < 1e-20` (a vanished cell), returns the zero
tuple. Each per-component term is itself guarded `ξ_k > 0 && ρ_k > 0`, so a
component with zero concentration or a (non-physical) zero density contributes
nothing. The all-phase basis is used for **gas transfer** (§5.6), where every
phase — including the gas that occupies pore space — participates in the
volumetric averaging; component-based permeability (§5.7.1) uses the
swelling-weighted basis (5.4).

### 5.1.3 Swelling-weighted (matrix) volume fractions and why surface properties use them

Thermal conductivity and surface emissivity are determined by the *condensed
matrix*, not by gas transiently occupying the pores. Pyrolysis.jl therefore
defines a **swelling-weighted volume fraction** in which every component
contributes the γ-fraction of its intrinsic volume:

$$
v_j^{\text{solid}} \;=\;
\dfrac{\gamma_j\,\xi_j/\rho_j(T)}{\displaystyle\sum_{k}\gamma_k\,\xi_k/\rho_k(T)},
\tag{5.4}
$$

so the fractions sum to one over the contributing (γ > 0) components. This is
`solid_volume_fractions(ξ, T, material)` (`mixture_rules.jl`). The swelling
factor here has **one universal meaning** — the fraction of a component's
intrinsic volume contributing to the bulk volume — shared with the dilation
rate $\theta = \sum_j\gamma_j(\mathrm{d}\xi_j/\mathrm{d}t)/\rho_j$ (§3.6,
§10) and the swelling mixture density (5.6).

The physical justification is the **ThermaKin convention**: "when γ = 0,
gaseous components do not contribute to the thermal conductivity"
(ThermaKin2Ds Reference, Eqs. 8–10), reproduced verbatim in the source
docstring. A non-swelling gas (`γ = 0`, the default for `GaseousComponent`,
§4) fills *existing* pore space and escapes; it neither forms a continuous
solid conduction path nor presents a radiating surface, so it carries zero
weight. A gas given `γ > 0` (an intumescent foaming agent that *creates* new
volume) is retained with weight γ.

> **Subtlety (γ = 0 condensed phases).** Because the weighting is universal,
> a *condensed* component with `γ = 0` (the `LiquidComponent` default, and
> the "rigid char" pattern) is consistently **volumeless**: zero weight here,
> zero strain in θ, and zero volume in the mixture density. Such a component
> contributes its mass and heat capacity but not its k or ε.

> **Void-cell blend.** Conductivity and emissivity consumers go through
> `mixing_volume_fractions`, which returns the swelling-weighted matrix
> fractions while the matrix volume fraction $v_m = \sum_j \gamma_j\xi_j/\rho_j$
> (`matrix_volume_fraction`) is at or above $10^{-3}$ (`MIXING_BLEND_HI`), the
> all-phase occupancy (5.3) at or below $10^{-4}$ (`MIXING_BLEND_LO`), and a
> cubic smoothstep blend $w\,v_j^{\text{solid}} + (1-w)\,v_j$,
> $w = s^2(3-2s)$, $s = (v_m - 10^{-4})/(10^{-3}-10^{-4})$, in between. A fully
> converted cell thus degrades to the gas-mixture conductivity instead of 0
> (which would zero both adjacent harmonic face means and, via Kirchhoff, its
> absorption), and because the matrix fractions are normalised by $v_m$ their
> derivatives grow as $1/v_m^2$: the smooth blend keeps
> $\partial k_{\rm eff}/\partial\xi$ bounded while a cell empties. An emptying
> cell stays in the system until deactivation removes it (§14), and with
> deactivation off it stays to the end of the solve, so the bound matters on
> every mesh.

The three bases and their consumers are summarized below; this mapping is the
single most error-prone part of effective-property bookkeeping, so it is stated
once authoritatively from the code. The authority is the state mapping
`derive!` (`src/States/derived.jl`), which evaluates every cell property once
per residual (§5.13):

| Basis | Symbol | Built by | Consumed by |
|---|---|---|---|
| Mass fraction | `Y_j` | inline in `_specific_vol_with_swelling_direct` | `mixture_density` (§5.2); the heat-capacity sum weights by the cell masses (§5.4) |
| All-phase volume fraction | `v_j` | `volume_fractions` | `effective_gas_transfer` |
| Swelling-weighted volume fraction | `v_j^solid` | `solid_volume_fractions` (through `mixing_volume_fractions`, which adds the void blend) | `effective_conductivity`, `mixture_emissivity`, `mixture_absorptivity`; `effective_permeability` reads `solid_volume_fractions` directly |

Note in particular that **gas transfer uses the all-phase basis `v_j`** while
**conductivity and emissivity use the solid basis `v_j^solid`**: in `derive!`,
`effective_conductivity`, `mixture_emissivity` and `mixture_absorptivity` each
build `mixing_volume_fractions`, while `effective_gas_transfer` builds
`volume_fractions`.

---

## 5.2 Mixture density with swelling

The **mixture (bulk) density** `ρ` answers a different question than
`ρ_tot`: it is the mass per unit *bulk volume* once each phase is allowed to
contribute, or not contribute, to that bulk volume according to its swelling
factor `γ_j`. The constitutive relation is a swelling-weighted harmonic
(specific-volume) average with the **universal** γ weighting (every
component, not gases only — the same weighting as θ and (5.4)),

$$
\frac{1}{\rho}
\;=\; \sum_{j} \gamma_j\,\frac{Y_j}{\rho_j(T)},
\tag{5.6}
$$

with `Y_j = ξ_j/ρ_tot` the mass fraction. This is `mixture_density(ξ, T,
material)` (`mixture_rules.jl`).

The semantics of `γ_j` follow directly from (5.6):

- **`γ_j = 1` (default for solids, `SolidComponent`).** The
  component contributes its full specific volume. When such a component depletes,
  its specific volume disappears and `1/ρ` decreases, i.e. the bulk *shrinks*.
  This is the ordinary charring/recession behaviour.
- **`γ_j = 0` (default for gases, `GaseousComponent`, **and for liquids**,
  `LiquidComponent`).** The term drops out of (5.6) entirely. Gas generated by
  a reaction fills pre-existing pore space and escapes through the matrix
  without inflating the bulk; a γ = 0 liquid is moisture bound in cell walls
  (ThermaKin `L SWELLING: 0`) — consistently volumeless across (5.4), (5.6),
  and θ. Set `γ = 1` for a liquid that physically occupies bulk volume (a
  melt layer).
- **`γ_j ∈ (0,1)`.** Partial contribution — the model for an intumescent
  blowing agent that *does* create new bulk volume as it is produced. Setting a
  gas's `γ_j` near 1 makes the generated gas expand the solid.

The volumetric strain rate `θ` that converts (5.6) into actual cell-volume change
(and hence into mesh motion under ALE) is derived in §10 (Volume Change and
Lateral Shrinkage) and §11 (The ALE Framework); here we only fix the
constitutive density. In the residual the mixture density has one consumer: the
film flux of a mass boundary condition, which multiplies the boundary cell's
gas-transfer coefficient by `ρ` (the product `λρ`), both at the cell's
positive-part concentrations and trial temperature (`GasTransportKernel`,
`src/Kernels/gas_transport.jl`; §9, §12).

> **Implementation note.** `mixture_density` computes the mass fractions
> *inline* via `_specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps, Val(NC))`
> (`mixture_rules.jl`), folding `Y_j = ξ_j · (1/ρ_tot)` into the
> specific-volume accumulation so that no intermediate `Y` tuple is allocated.
> Two empty-cell guards bracket the computation: `ρ_tot < 1e-20 ⇒ 0`, and
> `v_specific < 1e-20 ⇒ 0` (`mixture_rules.jl`). The swelling factor is
> applied per component (one universal γ semantics, §5.1.3), so a mixture may
> contain swelling and non-swelling components of any phase simultaneously.

---

## 5.3 Porosity

### 5.3.1 The canonical intrinsic-density form

**Porosity** `φ` is the bulk-volume fraction *not* occupied by the solid/liquid
skeleton. It is the only place the **intrinsic (skeletal) density** `ρ_{i,j}` is
used. The canonical definition is

$$
\varphi \;=\; 1 - \sum_{j\in\{\text{solid,liquid}\}} \frac{\xi_j}{\rho_{i,j}(T)},
\tag{5.7}
$$

with **no clamp**. The sum is `skeletal_volume_fraction(ξ, T, material)`
(`mixture_rules.jl`), the bulk-volume fraction the skeleton occupies; the pore
fraction is one minus it. The sum runs over **solid and liquid** components
only — the generated code branches on `is_solid(c) || is_liquid(c)`, and each
term is taken only if `ρ_intr > 0`. Gases are *excluded from the sum* because
they occupy the pore space being measured, not the skeleton.

The distinction between bulk and intrinsic density is essential here. The
*intrinsic* density is the density of the cell-wall material itself, ignoring any
micro-pores inside the pure phase; the *bulk* density is the apparent density of
that phase as it sits in the composite. For a porous solid these differ, and
using the intrinsic density in (5.7) correctly recovers the void fraction. A
component constructor takes either `ρ_intrinsic` or the pure phase's own
porosity `φ` (which, with a constant bulk density, sets
`ρ_intrinsic = ρ/(1 − φ)`). When a component is constructed with neither,
`get_intrinsic_density` falls back to the bulk density (`components.jl`), which
makes the summand `ξ_j/ρ_j`: the pure phase is taken to have no internal pores.
For a gas, `get_intrinsic_density(::GaseousComponent, ::Real) = 0.0`, consistent
with its exclusion from the sum.

Because (5.7) carries no clamp, a trial state can give `φ < 0` (a skeleton
larger than the cell) or `φ > 1` (negative condensed masses). Each consumer
guards the pore fraction as its own use requires:

1. **The state mapping**, under Darcy–Fick transport only, evaluates
   $$
   \varphi_i = \operatorname{smooth\_floor}\bigl(1 - s(\boldsymbol\xi^+_i, T_{\text{trial},i}),\ 10^{-12}\bigr),
   $$
   with `s` the skeletal volume fraction and `POROSITY_FLOOR = 1e-12`
   (`src/States/derived.jl`; the smooth floor is defined in §5.13). This `φ_i`
   feeds the **ideal-gas pressure closure** `P = N(ξ,T)/φ`, with
   `N = Σ_{g} ξ_g R_g T / M_g` the partial-pressure numerator
   (`ideal_gas_pressure`, which additionally takes `max(φ, 1e-12)`; §9), the
   **Carman–Kozeny permeability** closure (§5.7.2), and the face porosity of the
   Darcy advection (the arithmetic mean of the two cells' floored values, §9).
   Under Fickian transport no porosity is evaluated.
2. **The domain test** (`DomainCheck`, §15.5.4), under Darcy–Fick transport,
   reads the *unguarded* pore fraction of every live cell from its masses
   `m_{j,i}/V_i` at the trial temperature and rejects a proposed step when
   `φ < −tol` or `φ > 1 + tol`, `tol` being the skeletal volume fraction of the
   cell's mass tolerances (`src/Integration/checks.jl`).
3. **`background_gas_concentration`** (the pore gas at a given pressure, used to
   initialize a Darcy–Fick problem; §9) uses `max(1 − s, 0)`.
4. **`validate_darcy_initial_porosity`**, which `solve` runs through
   `validate_problem`, refuses a Darcy–Fick problem in which a cell with a
   positive effective permeability starts with an unclamped pore fraction below
   `1e-6` — the usual cause being a solid constructed without `ρ_intrinsic`
   (`src/Problem/problem_validation.jl`).
5. **Deactivation on a moving mesh** adds the pore volume (cell volume less the
   skeletal volume of its masses) of the cells it removes to the receiving
   cell, so pore space is conserved through the transfer (§14.3.3).

### 5.3.2 Porosity derivatives

Differentiating (5.7) term by term gives

$$
\frac{\partial \varphi}{\partial \xi_j}
=
\begin{cases}
-\dfrac{1}{\rho_{i,j}(T)}, & j\in\{\text{solid, liquid}\}\ \text{and}\ \rho_{i,j}>0,\\[8pt]
0, & j\in\text{gas (or }\rho_{i,j}\le 0).
\end{cases}
\tag{5.8}
$$

The sign is negative for the skeleton: adding solid mass reduces the void
fraction. Gases contribute zero — they fill pore space already accounted for by
the solid complement, so changing a gas concentration does not change `φ`.
These derivatives are not hand-coded: the colored Jacobian (§15) differentiates
the whole residual, the state mapping included, with forward-mode dual numbers,
so (5.8) is carried through the smooth floor exactly as the residual evaluates
it — the floor's derivative is zero at and below `10⁻¹²`, one at and above
`2·10⁻¹²`, and the cubic join's slope in between.

---

## 5.4 Effective volumetric heat capacity

### 5.4.1 Mixture heat capacity

The temperature-storage term of the energy equation (§3, §7) uses the
**mixture** volumetric heat capacity,

$$
(\rho c_p)_{\text{eff}}
\;=\; \sum_{j=1}^{N_C} \xi_j\, c_{p,j}(T) ,
\tag{5.9}
$$

summed over every component with no phase branch, so **gas sensible storage is
part of the left-hand side of the energy equation**. The units check:
`[kg·m⁻³] · [J·kg⁻¹·K⁻¹] = [J·m⁻³·K⁻¹]`.

Because the cell masses are the states, the energy row does not divide by a
volume (§7.3) and reads the extensive form of (5.9), the **heat-capacity sum**
of the cell,

$$
C_i \;=\; \operatorname{smooth\_floor}\Bigl(\sum_{j=1}^{N_C} m_{j,i}\,c_{p,j}(T_{\text{trial},i}),\ \ \rho c_{p,\text{floor}}\,V_i(0)\Bigr)
\qquad [\text{J·K}^{-1}],
\tag{5.10}
$$

evaluated in the state mapping (`cp_sum`, `src/States/derived.jl`) with the
masses as stored and the trial temperature (§5.13), so that `C_i/V_i` is the
discrete `(ρc_p)_eff`. The sum (`heat_capacity_sum`, `src/States/caloric.jl`)
does not evaluate the heat capacity of a component whose mass is exactly zero:
its term would add an exact zero for every finite heat capacity, and the test
`iszero` of a dual number reads its partials too, so no derivative is lost.

The floor is `RHO_CP_FLOOR = 10³` J·m⁻³·K⁻¹ (`src/Materials/constants.jl`) times
the cell's initial volume (`CellReference.capacity_floor`,
`src/States/tolerances.jl`): the thermal mass of pore gas at ambient
(ρ_g c_{p,g} ≈ 1.2 × 1000; Gpyro floors the same quantity). Because every
component is summed, the sum reaches the floor only in a cell emptied of all
mass — a burned-out cell whose pore gas has also left; without it the `1/C_i`
factor in `dT_i/dt` would blow up there. Virgin solids carry ρc_p ~ 10⁵–10⁶
J·m⁻³·K⁻¹, so the floor never binds in a cell that holds material. The floor is
the C¹ smooth floor of §5.13, applied in the state mapping, so the residual and
the colored Jacobian, which differentiates the state mapping, see the same
value and its consistent derivative.

### 5.4.2 Why the gas is stored

The one-temperature model puts the pore gas at the matrix temperature (§3.1.2),
so the gas stores sensible heat `ξ_g c_{p,g} ∂T/∂t` exactly like every other
component. Its mass is an integrated state, so the term costs nothing; dropping
it — Gpyro's default *quasi-steady-gas approximation*, which its `FULL_QSG`
option reverses, and the implicit convention of FDS, which has no in-depth gas —
would be an approximation with no benefit. The in-pore gas heat capacity is
roughly `φ ρ_g c_{p,g}` against the matrix `(1−φ) ρ_s c_{p,s}`, with `ρ_s` the
intrinsic density of the matrix. At 1 atm and 700 K, with comparable specific
heats and `ρ_s = 1500` kg·m⁻³, the ratio is about 3·10⁻⁴ for air in a matrix of
porosity 0.5 (ρ_g ≈ 0.50 kg·m⁻³) and about 5·10⁻³ for volatiles of molar mass
0.1 kg·mol⁻¹ in a char of porosity 0.8 (ρ_g ≈ 1.74 kg·m⁻³): from a few hundredths
to about half a percent of the storage, scaling with pressure. Small, but kept
exact because it is free.

The gas *flux* relative to the skeleton enters the energy equation only through
the sensible enthalpy it carries across the faces (§7.4). This is the ThermaKin
Eq. 17 convention term for term (§3.3.2).

No mass-averaged specific heat of the mixture enters the residual; the
enthalpy-form tolerance of a cell uses the mass-averaged heat capacity of its
initial contents as a scale (§15.5.6).

### 5.4.3 Sensible enthalpies and the extended caloric relation

The energy row reads, besides `C_i`, the **sensible enthalpy** of every
component on the common datum `T_ref = 298.15` K,

$$
h_j(T) \;=\; \int_{T_{\text{ref}}}^{T} c_{p,j}(T')\,\mathrm dT'
\qquad [\text{J·kg}^{-1}].
\tag{5.11}
$$

A property law evaluates at its argument clamped to the property bracket
`[200, 3000]` K (§4.4), so `c_{p,j}` is constant outside the bracket and `h_j`
continues linearly there. Constant, linear, polynomial and tabulated
heat-capacity laws are integrated in closed form, piecewise across the bracket
ends and the table knots, so `dh_j/dT` is the law itself to round-off. Any other
law is integrated once, at setup, into a `TabulatedEnthalpy`: the integral at
1 K knots over the bracket, accumulated by five-point Gauss–Legendre quadrature
on each interval, and evaluated between the knots by the cubic Hermite
interpolant of the integral and the law, whose exact derivative is the slope
the residual uses (`enthalpy_slope`). The state mapping evaluates `h_j(T_i)` at
the cell temperature itself, for every component (`src/States/caloric.jl`).

`h_j` enters the energy row in three places (§7.3, §7.4): the temperature form
subtracts `Σ_j h_j(T_i) Ṁ_{j,i}`, each face carries `Σ_j F_{j,f} h_j(T_f)`, and
the inventory form adds the formation term `Σ_j ω_{j,i} V_i h_j(T_i)`.

In the inventory form (`energy_state = :enthalpy`) the energy state of cell `i`
is its enthalpy inventory, and the state mapping recovers `T_i` from it by
inverting the **extended caloric relation**

$$
H_i \;=\; \sum_{j} m^{+}_{j,i}\,h_j(T_i) \;+\; C_{\mathrm{reg},i}(\boldsymbol m_i)\,\bigl(T_i - T_{\text{ref}}\bigr),
\tag{5.12}
$$

with `m⁺_{j,i}` the one-sided positive part of the mass at the cell's scale
`δ_i` (§5.13) and `C_reg,i` a C¹ capacity switch: with
`C_g = Σ_j m⁺_{j,i} c_{j,min}` the guaranteed heat capacity (`c_{j,min}` half the
smallest value of the component's law sampled at 1 K over the bracket) and
`C_floor,i = ρc_{p,floor} V_i(0)`, it is zero while `C_g ≥ 2 C_floor,i` and
`C_floor,i (1 − C_g/(2C_floor,i))²` below. The relation is strictly increasing
in `T` with slope at least `C_floor,i`, so it has exactly one root for every
real inventory, and it equals the physical `Σ_j m_{j,i} h_j(T_i)` on every state
in which each component is absent or holds at least `δ_i` and the guaranteed
capacity is at least `2 C_floor,i` (`CaloricExtension`, `src/States/caloric.jl`).
The inversion is a bracketed Newton iteration started at the cell's previous
temperature (`invert_caloric`); a dual-number inventory is differentiated
through the defining relation, by one Newton correction of the primal root per
nesting level, which carries the implicit-function derivatives exactly. A
negative heat capacity is refused at setup; a gas with zero heat capacity is
admitted and adds nothing to `C_g`.

---

## 5.5 Thermal conductivity mixing

The effective thermal conductivity is computed from the swelling-weighted
volume fractions `v_j^solid` (5.4) and the per-component conductivities
`k_j(T)` under one of four mixing rules selected by `material.mixing.k.rule`:
`effective_conductivity(ξ, T, material; skip_absent)` builds
`v = mixing_volume_fractions(…)` (the (5.4) weights, with the void blend of
§5.1.3) then dispatches to `_conductivity_from_v` (`effective_properties.jl`).
The default rule is `WEIGHTED` with `β = 0.5` (`materials.jl`).

The state mapping evaluates `k_i = effective_conductivity(ξ⁺_i, T_trial,i,
material; skip_absent = Val(true))` and floors it smoothly at
`CONDUCTIVITY_FLOOR = 1e-6` W·m⁻¹·K⁻¹ (`src/States/derived.jl`), far below any
material's conductivity: a trial state with no mass left has none, and the
surface-temperature solve of a boundary face divides by the cell conductivity
(§12). With `skip_absent = Val(true)` the parallel, series and weighted rules
do not evaluate the conductivity law of a component whose volume fraction is
exactly zero (its partials included, under automatic differentiation); the
result is the same for every finite conductivity.

### 5.5.1 Parallel and series bounds

The two classical bounds correspond to the geometric extremes of phase
arrangement relative to the heat-flow direction. With

$$
k_\parallel = \sum_j v_j^{\text{solid}}\, k_j(T),
\qquad
\frac{1}{k_{\text{series}}} = \sum_j \frac{v_j^{\text{solid}}}{k_j(T)},
\tag{5.13}
$$

the **PARALLEL** rule returns `k = k_∥` and the **SERIES** rule returns
`k = k_series`. Physically:

- **Parallel (upper bound, Voigt/Wiener):** all phases span the slab in the
  heat-flow direction and conduct independently in parallel; thermal
  conductances add. Appropriate when every phase forms a continuous path along
  the gradient.
- **Series (lower bound, Reuss/Wiener):** phases are stacked perpendicular to
  the gradient so heat passes through each in turn; thermal *resistances* add,
  giving a harmonic mean. Appropriate for stratified layering.

The true effective conductivity of any two-phase mixture lies between these
Hashin–Shtrikman-bracketed bounds.

> **Implementation note.** Both sums are produced in one pass by the generated
> `_conductivity_sums(v, T, comps, Val(NC), Val(skip_absent))`
> (`effective_properties.jl`). The series accumulation is guarded:
> `inv_k_series += v[i]/k` is taken only when `v[i] > 1e-20 && k > 1e-20`, so a
> phase with zero volume or (pathologically) zero conductivity does not blow up
> the harmonic sum. After the loop
> `k_series = inv_k_series > 1e-20 ? 1/inv_k_series : 0`
> (`effective_properties.jl`).

### 5.5.2 Weighted blend

The **WEIGHTED** rule linearly interpolates between the bounds with a single
user-supplied parameter `β ∈ [0,1]`:

$$
k_{\text{eff}} \;=\; \beta\, k_\parallel + (1-\beta)\, k_{\text{series}} .
\tag{5.14}
$$

`β = 1` recovers the parallel upper bound, `β = 0` the series lower bound, and the
default `β = 0.5` (`MixingSpec` default, `materials.jl`) is the arithmetic mean
of the two bounds. The blend is the package's default conductivity model
(`materials.jl`) because, lacking microstructural information, the mean of the
rigorous bounds is the least-biased single estimate. The `β` value is validated to
lie in `[0,1]` at `MixingSpec` construction (`materials.jl`).

### 5.5.3 Bruggeman effective-medium theory

For random heterogeneous microstructures — and especially for the self-similar
pore structure of forming char — the **BRUGGEMAN** symmetric effective-medium
closure is more faithful than a fixed blend. It treats each phase as a spherical
inclusion embedded in the *effective* medium and demands self-consistency: the
volume-averaged polarization vanishes,

$$
\sum_j v_j^{\text{solid}}\,
\frac{k_j - k_{\text{eff}}}{k_j + 2 k_{\text{eff}}} \;=\; 0 .
\tag{5.15}
$$

Unlike the parallel/series/weighted forms, (5.15) is *implicit* in `k_eff` and
can capture percolation-like behaviour (the effective conductivity dropping
sharply when a continuous phase fractures) and genuinely non-linear composition
dependence. The factor `+2 k_eff` is the spherical-inclusion depolarization
factor for an isotropic 3-D medium; the symmetric form (no distinguished host
phase) reflects an isotropic, statistically self-similar microstructure.

#### Newton–Raphson solve

Equation (5.15) is solved by Newton's method in
`bruggeman_effective_property(v, k; tol = 1e-10, maxiter = 100)`
(`effective_properties.jl`). Writing the residual
`f(k_{\text{eff}}) = \sum_j v_j (k_j - k_{\text{eff}})/(k_j + 2 k_{\text{eff}})`,
the exact derivative is

$$
\frac{df}{dk_{\text{eff}}}
= \sum_j v_j\,\frac{-3\,k_j}{(k_j + 2 k_{\text{eff}})^2} ,
\tag{5.16}
$$

derived in the source comment line-by-line (`effective_properties.jl`):
the numerator of `d/dk_eff [v(k−k_eff)/(k+2k_eff)]` collapses
`−(k+2k_eff) − 2(k−k_eff)` to `−3k`. The Newton update is
`k_eff ← k_eff − f/df`, with the following exact guards and conventions taken from
the implementation:

- **Degenerate volume.** If `Σ v_j < 1e-20`, return 0
  (`effective_properties.jl`).
- **Single active phase.** Active phases are those with `v_i > 1e-20`; if exactly
  one is active, return its (volume-weighted) parallel mean `k_parallel`
  immediately (`effective_properties.jl`). This is exact, since a
  single-phase Bruggeman medium is that phase.
- **Seed.** Every call starts the Newton iteration from the parallel mean
  `k_parallel`, the volume-weighted arithmetic mean of the active phases
  (`effective_properties.jl`); no value is carried from one evaluation to the
  next, so the result is a function of `(v, k)` alone.
- **Per-term contributions.** A phase enters the residual when `v_i > 1e-20`.
  A conducting phase (`k_i > 1e-20`) contributes the standard term; an
  **insulating** phase (`k_i ≤ 1e-20`) contributes its exact k → 0 limit
  `−v_i/2` to `f` (zero `df` contribution). The insulating term is what
  produces the percolation threshold: a conducting phase dispersed in an
  insulator yields `k_eff = 0` until its volume fraction exceeds 1/3
  (e.g. `v = (0.5, 0.5)`, `k = (1, 0)` → `k_eff = 0.25`).
- **Below-threshold collapse.** Below the percolation threshold the root is
  `k_eff = 0` and `f < 0` for every `k_eff > 0`, so the step-capped Newton
  decays `k_eff` geometrically; once `k_eff < 1e-30` the solve returns 0
  immediately instead of burning `maxiter`.
- **Singular-derivative bail-out.** If `|df| < 1e-20` the iteration breaks
  (it cannot take a Newton step; the problem is locally ill-posed)
  (`effective_properties.jl`).
- **Step limiter (positivity).** The raw step `dk = −f/df` is clipped so the
  update never reduces `k_eff` by more than 90 %: if
  `k_eff + dk < 0.1·k_eff`, set `dk = −0.9·k_eff`
  (`effective_properties.jl`). This prevents Newton from overshooting
  into a negative conductivity.
- **Convergence and return.** Iterate up to `maxiter = 100`, stopping when
  `|f| < tol = 1e-10`; the return is `max(k_eff, 0.0)`
  (`effective_properties.jl`).

> **Implementation note.** The component conductivities are gathered into a tuple
> by the generated `_conductivity_tuple(T, comps, Val(NC))`
> (`effective_properties.jl`), so the Bruggeman path is allocation-free.
> The dispatch in `_conductivity_from_v` checks `mixing == BRUGGEMAN` first and
> only otherwise computes the parallel/series sums
> (`effective_properties.jl`); the Bruggeman path evaluates every component's
> conductivity, so `skip_absent` does not apply to it.

---

## 5.6 Effective gas-transfer coefficient

The diffusive gas flux (§9) uses an effective Fickian transfer coefficient
`λ_eff` [m²·s⁻¹], homogenized from the per-component `λ_j(T)`. Crucially — and
unlike conductivity — gas transfer is averaged over the **all-phase** volume
fractions `v_j` (5.3), *not* the solid basis:

$$
\lambda_{\text{eff}}(T,\boldsymbol\xi)
\;=\; \text{MixingRule}\bigl(v_j,\ \lambda_j(T)\bigr),
\tag{5.17}
$$

`effective_gas_transfer(ξ, T, material)` builds `v = volume_fractions(…)` and
dispatches to `_gas_transfer_from_v` (`effective_properties.jl`). The state
mapping evaluates `λ_i = effective_gas_transfer(ξ⁺_i, T_trial,i, material)`
only when gas can move (§5.13). Which transport a problem has is decided once,
at setup (`_gas_transport`, `src/Integration/setup.jl`): none without a gas
component; Darcy–Fick when some condensed component carries a permeability law
that is not identically zero (`has_darcy_flow`, §5.7); Fickian when some
component carries a gas-transfer coefficient that is not identically zero or a
boundary condition passes mass; none otherwise. The supported rules mirror
conductivity:

$$
\lambda_\parallel = \sum_j v_j\,\lambda_j,\qquad
\frac{1}{\lambda_{\text{series}}}=\sum_j\frac{v_j}{\lambda_j},\qquad
\lambda_{\text{eff}}=\beta\lambda_\parallel+(1-\beta)\lambda_{\text{series}}.
\tag{5.18}
$$

The default rule for `λ` is **PARALLEL** with `β = 0.5` (`materials.jl`).

**Degenerate transport.** `SolidComponent` defaults to $\lambda_j = 0$. A
condensed component with no gas-transfer coefficient contributes nothing to
(5.18), so a cell of that component holding no gas has $\lambda_{\text{eff}}
\approx v_{\rm gas}\lambda_{\rm gas} \to 0$ under PARALLEL (and exactly zero
under SERIES), and the guarded face mean (5.29) goes to zero continuously as
either side does (§5.12). Gas transport then stops where the gas fraction stops:
the gas field of a gas-producing material ends in a sharp front, a trial state
that drives the front cell's gas concentration to zero closes both of its faces
and flattens its residual row, and the integrator rejects step after step on
that front. A small coefficient on the condensed components keeps the transport
non-degenerate; the Landau ablation case V-CPL-02 gives its solid
$\lambda = 10^{-6}$ m²/s for this reason (`test/verification/cases/V-CPL-02.jl`).
`validate_problem`, which `solve` calls, therefore warns when a material with a
gas-producing reaction and Fickian transport (no permeable component) has a
condensed component with $\lambda_j = 0$
(`warn_degenerate_gas_transport`, `src/Problem/problem_validation.jl`): "the
integrator rejects step after step on it. Give the condensed components a small
λ (for example `λ = 1e-6` m²/s)."

### 5.6.1 Bruggeman is rejected for gas transfer

Applying the effective-medium closure (5.15) to a diffusivity is **rejected**:
`_gas_transfer_from_v` throws an `ArgumentError` if `mixing == BRUGGEMAN`
(`effective_properties.jl`), with the message that the EMT closure "has
no physical justification for diffusivities." The reasoning is that a Fickian
transfer coefficient is already an average over a heterogeneous pore network
(molecular and Knudsen regimes), so wrapping it in a *second* effective-medium
self-consistency relation derived for a thermodynamic-like coefficient is not
physically supported. (The `effective_gas_transfer` docstring
(`effective_properties.jl`) explicitly states that `BRUGGEMAN` is **not**
supported for gas transfer and lists only `PARALLEL`/`SERIES`/`WEIGHTED`,
matching the running code, which is authoritative.)

> **Implementation note.** The series guard for `λ` is more elaborate than for
> `k`: in `_gas_transfer_sums` (`effective_properties.jl`), once any
> active phase has `λ ≤ 1e-20` the running `inv_λ_series` is set to `Inf` and the
> rest of the sum is short-circuited, which drives `λ_series → 0`. This is the
> correct series limit: a single impermeable-to-diffusion layer blocks the
> series path entirely.

---

## 5.7 Effective permeability

Permeability `κ` [m²] is the resistance-to-flow property in Darcy's law
`v_g = −(κ/μ)∇P` (§9). Only **solid and liquid** components carry a permeability;
gases flow *through* the matrix and `get_permeability(::GaseousComponent, T) =
0.0` (`components.jl`). Pyrolysis.jl offers two **mutually exclusive**
closures, selected by `material.mixing.κ.rule`
(`effective_permeability(ξ, T, material, φ)`, `effective_properties.jl`; the
porosity argument `φ` is read by the Carman–Kozeny rule only). The default rule
for `κ` is **WEIGHTED** with `β = 0.5` (`materials.jl`).

Darcy–Fick transport, and with it the permeability, is active when some
condensed component carries a permeability law that is not identically zero
(`has_darcy_flow`, `src/Materials/materials.jl`), whatever the mixing rule: a
material whose components all have `κ = 0` is Fickian even under
`CARMAN_KOZENY`. The state mapping then evaluates the cell mobility from the
permeability at the positive-part concentrations, the trial temperature and the
smoothly floored pore fraction of §5.3.1 (§5.7.4, §5.13).

### 5.7.1 Component mixing

When the rule is `PARALLEL`, `SERIES`, or `WEIGHTED`, the effective permeability
is a volume-fraction mixture of the component permeabilities, exactly as for
conductivity: on the **renormalized swelling-weighted basis** `v_j`
(`solid_volume_fractions` (5.4) — γ = 0 components carry no weight and the
surviving weights sum to 1):

$$
\kappa_\parallel = \sum_j v_j\,\kappa_j,\qquad
\frac{1}{\kappa_{\text{series}}}=\sum_j\frac{v_j}{\kappa_j},\qquad
\kappa_{\text{eff}}=\beta\kappa_\parallel+(1-\beta)\kappa_{\text{series}},
\tag{5.19}
$$

`_permeability_from_v` (`effective_properties.jl`), with the parallel/series
sums from `_permeability_sums`. An **impermeable phase** (`κ ≈ 0`, the
`SolidComponent` default) with `v_j > 0` makes the series (Wiener lower) bound
**exactly zero** — gas cannot cross a series-stacked impermeable layer; the
accumulator signals this with `1/κ_series = ∞`, the same convention as the
gas-transfer channel. The series guard floor is tighter than for
conductivity — `κ > 1e-30` rather than `1e-20` — because permeabilities
themselves are tiny (virgin wood ~10⁻¹⁵ m², char ~10⁻¹² m²). Representative
component values from the source docstring: virgin wood `~10⁻¹⁵ m²` (tight
cellular structure), char `~10⁻¹² m²` (open channels), impermeable `0.0 m²`.

### 5.7.2 Carman–Kozeny porosity closure

When the rule is `CARMAN_KOZENY`, permeability is *not* a per-component mix but a
porosity-based closure:

$$
\kappa \;=\; \min\!\left(\frac{\ell^2\,\varphi^3}{180\,(1-\varphi)^2},\;\ell^2\right),
\tag{5.20}
$$

`_carman_kozeny_permeability(φ, ℓ)` (`effective_properties.jl`), where
`ℓ` is the characteristic particle/pore diameter [m] supplied through the
`MixingSpec` `length_scale` and `φ` is the porosity the caller passes: in the
residual, the pore fraction (5.7) smoothly floored at `10⁻¹²` (§5.3.1), so the
permeability is smooth in the state. The constant **180** is the Kozeny
constant for a packed bed of spheres. The
`φ³/(1−φ)²` dependence captures the steep growth of permeability as porosity
rises during virgin → char conversion: as `φ → 0` (fully solid), `κ → 0`
(impermeable). As `φ → 1` the packed-bed closure loses validity and diverges,
so `κ` is capped at the physical free-flow scale `ℓ²` (an open channel of width
`ℓ` has `κ ~ ℓ²/32`; `ℓ²` is a generous upper bound). The cap binds for
`φ > 0.9328`, where `180(1−φ)² < φ³`; a tiny `1e-30` denominator guard remains
solely to avoid `Inf`/`NaN` at exactly `φ = 1` before the `min` resolves
(`effective_properties.jl`).

Because Carman–Kozeny needs the porosity, it is evaluated through
`effective_permeability(ξ, T, material, φ)` — calling the volume-fraction
variant `_permeability_from_v` with this rule throws an `ArgumentError`
explaining that concentrations are required (`effective_properties.jl`).
`MixingSpec(CARMAN_KOZENY; length_scale=ℓ)` fails fast at construction if `ℓ`
is not a positive, finite number (`materials.jl`).

The closure derivative needed for the Jacobian (§5.15) is not hand-coded: the
colored Jacobian (§15) differentiates the state mapping, `effective_permeability`
and the floored pore fraction included, with dual numbers, so the Carman–Kozeny
chain through the porosity's composition dependence is carried as evaluated. On
the capped branch the `min` gives the value-consistent zero `dκ/dφ` — value and
derivative saturate at the same point (a kink of the derivative at the cap,
§5.15).

### 5.7.3 Bruggeman is rejected for permeability

As with gas transfer, `BRUGGEMAN` is rejected for permeability:
`_permeability_from_v` throws an `ArgumentError` if `mixing == BRUGGEMAN`
("EMT is only defined for scalar transport coefficients … use
PARALLEL/SERIES/WEIGHTED, or CARMAN_KOZENY", `effective_properties.jl`).
Permeability is a flow-resistance of the matrix geometry, not a polarizable
scalar property to which the symmetric EMT applies. Carman–Kozeny is the
physically motivated porosity-based alternative for the percolation-like
virgin → char transition.

### 5.7.4 Face-averaged permeability

For the Darcy velocity at a face, the resistance-in-series combination is
applied to the **mobility** κ/μ rather than to κ alone. The state mapping forms
each cell's mobility

$$
\text{mob}_i = \frac{\kappa_{\text{eff}}(\boldsymbol\xi^+_i, T_{\text{trial},i}, \varphi_i)}
{\operatorname{smooth\_floor}\bigl(\mu(T_{\text{trial},i}),\ \mu_{\text{floor}}\bigr)},
\qquad \mu_{\text{floor}} = 10^{-12}\ \text{Pa·s}
$$

(`mob`, `src/States/derived.jl`; `μ` is the material's `gas_viscosity`, §5.9,
and `MU_FLOOR` the floor), and the gas-transport kernel combines the two sides
of an interior face with the guarded distance-weighted harmonic mean (§5.12)
with the mobility guard `MOBILITY_GUARD = 1e-30` m²·Pa⁻¹·s⁻¹ in place of the
default `10⁻²⁰` (`src/Kernels/gas_transport.jl`):

$$
\text{mob}_f = \frac{(d_o + d_n)\,\text{mob}_o\,\text{mob}_n}
{d_o\,\text{mob}_n + d_n\,\text{mob}_o + 10^{-30}(d_o + d_n)} .
$$

The guard lies far below the mobility of any permeable solid (κ = 10⁻¹⁸ m² and
μ = 5·10⁻⁵ Pa·s give 2·10⁻¹⁴), and the mean goes to zero continuously as either
side's mobility vanishes: an impermeable cell is an infinite series resistance,
so gas cannot cross into it, and the face's derivative stays finite.

---

## 5.8 Radiative properties: emissivity and absorption

### 5.8.1 Effective emissivity

The surface emissivity is a **solid**-volume-fraction average of the component
emissivities, each evaluated at the cell temperature:

$$
\varepsilon_{\text{eff}}(T,\xi) \;=\; \sum_j v_j^{\text{solid}}\,\varepsilon_j(T),
\tag{5.22}
$$

`mixture_emissivity(ξ, T, material)` builds the mixing weights
(`mixing_volume_fractions`, the solid basis with the void blend of §5.1.3) and
sums `v_j · get_emissivity(comp, T)` over the components
(`effective_properties.jl`). **No clamp** is applied, and none is needed: the
weights are a convex combination (non-negative, summing to one, in the blend as
on either side of it), and every component's emissivity law is validated to lie
in `[0, 1]` over the property bracket when the component is built
(`src/Materials/components.jl`), so the sum stays in `[0, 1]` to round-off and
its derivative has no kink. The solid basis is used for the same physical reason
as conductivity (§5.1.3): the radiating surface is the solid/liquid material,
not trapped non-swelling gas. Component emissivities are property functions
(`ε` accepts the same forms as `k`: number, linear tuple, callable, property
function) with defaults of `0.9` for `SolidComponent`, `0.95` for
`LiquidComponent` (condensed liquids are strong emitters — water ≈ 0.96), and
`0.0` for `GaseousComponent` (`src/Materials/components.jl`); a cell holding only
gases of the default emissivity therefore has `ε_eff = a_eff = 0`. Because the
sum is evaluated generically over the scalar type, a dual-number temperature
carries $\partial\varepsilon_{\text{eff}}/\partial T$ into the colored Jacobian
(§15). The effective emissivity is the surface property of the re-radiation
boundary condition (§12): the surface balance of a boundary face reads the
emissivity of its cell.

### 5.8.2 Effective surface absorptivity

The absorptivity for *incident* (heater-spectrum) radiation is mixed with the
same weights,

$$
a_{\text{eff}}(T,\xi) \;=\; \sum_j v_j^{\text{solid}}\,a_j,\qquad
a_j = \begin{cases}
\texttt{surface\_absorptivity}_j & \text{if set}\\
\varepsilon_j(T) & \text{otherwise (Kirchhoff)}
\end{cases}
\tag{5.22b}
$$

via `mixture_absorptivity(ξ, T, material)`, which sums
`v_j · get_absorptivity(comp, T)`; like (5.22) it is unclamped, a component's
explicit absorptivity being validated in `[0, 1]` at construction. The
Kirchhoff branch resolves at compile time (`Nothing` vs `Float64` is a type
parameter of the component), and for a material in which no component sets a
surface absorptivity the two generated sums perform the same floating-point
operations in the same order, so
$a_{\text{eff}} \equiv \varepsilon_{\text{eff}}$ bit-for-bit — a material that
sets no absorptivity absorbs exactly as Kirchhoff's law predicts
(`test/materials/test_emissivity_absorptivity.jl`, "a_eff: Kirchhoff invariant
and per-component override"). A component with an explicit $a_j$ decouples the
absorbed heater flux from the temperature-dependent emission: an emissivity that
falls with temperature (integrating-sphere data for virgin wood, for instance)
does not drag the absorbed cone flux down with it. $a_{\text{eff}}$ is consumed
by a `RadiativeFluxBC` with `absorptivity = nothing`: in the surface balance of
its face under `NO_RADIATION` and `SURFACE_ABSORPTION` (§12.2.6), and under
`BEER_LAMBERT` as the fraction of the incident flux that enters the first live
cell, the absorptivity of that cell (`incident_intensity`,
`src/Kernels/radiation.jl`; §8.3). Surface emission always uses
$\varepsilon_{\text{eff}}$.

### 5.8.3 Effective absorption coefficient

The volumetric (extinction) absorption coefficient of the in-depth radiation
model (`BEER_LAMBERT`, §8) is a **concentration**-weighted sum of the component
*mass* absorption coefficients `α_j` [m²·kg⁻¹]:

$$
\alpha_{\text{eff}} \;=\; \sum_j \xi_j\, \alpha_j \quad [\text{m}^{-1}],
\tag{5.23}
$$

`effective_absorption(ξ, material)` → `_sum_absorption`
(`effective_properties.jl`), evaluated by the state mapping at the positive-part
concentrations, and only under `BEER_LAMBERT`. The dimensional bookkeeping is
the key point: a *mass* absorption coefficient `α_j` [m²·kg⁻¹] times a mass
concentration `ξ_j` [kg·m⁻³] yields a *volumetric* absorption coefficient
[m⁻¹], which is what the optical thickness of a cell, `τ_i = α_i Δz_i`,
requires. Because the weight is the concentration `ξ_j` (not a normalized
fraction), `α_eff` correctly scales down as material is consumed — a thinning,
depleting layer becomes optically thinner. Component absorption coefficients
default to `0.0` (no absorption) and carry no temperature dependence.

---

## 5.9 Gas viscosity (Sutherland's law)

The dynamic viscosity of the escaping gas, needed for Darcy's law (§9), is the
material's `gas_viscosity` property. Its default is **Sutherland's law** for air,

$$
\mu(T) \;=\; \mu_{\text{ref}}
\left(\frac{T}{T_{\text{ref}}}\right)^{3/2}
\frac{T_{\text{ref}}+S}{T+S},
\tag{5.24}
$$

carried by `SutherlandProperty()` (§4.4.7) with the **canonical Sutherland air
constants** (e.g. White, *Viscous Fluid Flow*):

| Parameter | Default | Source |
|---|---|---|
| `μ_ref` | `1.716 × 10⁻⁵` Pa·s | `SutherlandProperty()`, `property_functions.jl` |
| `T_ref` | `273.15` K | `SutherlandProperty()`, `property_functions.jl` |
| `S` (Sutherland constant) | `110.4` K | `SutherlandProperty()`, `property_functions.jl` |

These give `μ(300 K) = 1.846 × 10⁻⁵` Pa·s and `μ(600 K) = 3.016 × 10⁻⁵` Pa·s.

> **Caution (reference temperature).** The viscosity reference temperature here
> is **`T_ref = 273.15 K`**, the temperature at which `μ_ref` is quoted for air —
> this is a *local default of the viscosity correlation* and is distinct from the
> global enthalpy datum `T_REF = 298.15 K` (§2). Do not conflate the two;
> `SutherlandProperty()` carries `T_ref = 273.15` in its own field.

Sutherland's law is the Chapman–Enskog rigid-sphere result for a single gas and
captures the *increase* of gas viscosity with temperature (`μ ∝ T^{1/2}` in the
high-`T` limit), the opposite of liquid behaviour. It is valid roughly over
`200–2000 K` for air-like gases.

### Material-supplied viscosity

The air constants are only the **default**. Each `Material` carries its own
pore-gas viscosity curve in the `gas_viscosity` field, and under Darcy–Fick
transport the state mapping evaluates it at the cell's trial temperature and
floors the result smoothly at `MU_FLOOR`,
`μ_i = smooth_floor(material.gas_viscosity(T_trial,i), 10⁻¹² Pa·s)`, the
denominator of the cell mobility (§5.7.4; `src/States/derived.jl`). A material
that says nothing carries `SutherlandProperty()` — the air triple above — so
this reduces to (5.24).

The field accepts anything `to_property` understands (a number, a `(p0, p1)`
tuple, a `Function`, or a property function), so another Sutherland triple or a
fitted correlation — a `PolynomialProperty` of ascending coefficients, for
instance a mixture-averaged fit for the material's pyrolyzates — drops in
directly:

```julia
using Pyrolysis

virgin = SolidComponent(:virgin, ρ = 500.0, c = 1500.0, k = 0.15, κ = 1e-15,
                        ρ_intrinsic = 1500.0)
char = SolidComponent(:char, ρ = 150.0, c = 1100.0, k = 0.08, κ = 1e-12,
                      ρ_intrinsic = 1500.0)
gas = GaseousComponent(:gas, M = 0.03, c = 1500.0, k = 0.03, λ = 1e-5)
wood = Material(
    name = :wood,
    components = (virgin, char, gas),
    reactions = (Reaction(:pyrolysis, :virgin => (:char, :gas); A = 1e8, E = 125e3,
                          h = 300e3, yields = (0.3, 0.7)),),
    gas_viscosity = SutherlandProperty(2.0e-5, 350.0, 120.0),  # μ_ref, T_ref, S
)
wood.gas_viscosity(600.0)   # 2.930e-5 Pa·s
```

**Why the guards matter.** μ enters the model *only* as the Darcy mobility
κ/μ (§9), never on its own. An error in μ is therefore exactly the reciprocal
error in permeability — and a μ that reaches zero or changes sign does not
degrade the answer, it reverses or blows up the gas velocity. Hence:

- the curve is a property law, so it evaluates on the property bracket
  $[200, 3000]$ K (§4.4), and the state mapping evaluates it at the trial
  temperature, itself inside the bracket: a wild trial state cannot evaluate a
  fitted polynomial absurdly far outside its range (the Chapman–Enskog
  diffusivity formula, which is not a property law, keeps its own 100 K floor
  `T_GAS_PROPERTY_FLOOR`);
- the result is floored smoothly at `MU_FLOOR` (1e-12 Pa·s, seven orders below
  any real gas) so an extrapolated fit cannot go non-positive;
- `Material`'s constructor probes the curve at 200, 300, 500, 800, 1200 and
  1800 K and refuses one that is non-positive or non-finite at any of them, in
  the spirit of `MixingSpec`'s `CARMAN_KOZENY`/`length_scale` check
  (`src/Materials/materials.jl`).

Because the bracket belongs to every property law and the floor is applied in
the state mapping, every curve a material can carry is protected, not just the
built-in Sutherland form.

**Accuracy.** Air is within ~2 % for an N₂-swept pore space, but overstates the
viscosity of pyrolysis volatiles by roughly 1.3–1.7× over 300–1000 K (the
`Material` docstring) — a bias in one direction everywhere, growing as the pore
space fills with volatiles.

**Composition dependence.** The model carries **one curve per material**, a
function of the temperature alone: an average over the composition that material
is expected to produce. It cannot resolve the pore gas shifting from background
air to volatiles as pyrolysis proceeds. A composition-dependent rule (Wilke)
would be evaluated in the same place, the state mapping — see §5.16 item 1.

---

## 5.10 Chapman–Enskog binary diffusivity

For materials where an explicit per-component `λ_j` is not supplied, the package
can compute a first-principles binary diffusion coefficient from
**Chapman–Enskog kinetic theory** (Bird, Stewart & Lightfoot). For a species `A`
diffusing through a background gas `B` (taken as air),

$$
D_{AB} \;=\; 5.95\times10^{-4}\,
\frac{T^{3/2}\sqrt{\,1/M_A + 1/M_B\,}}{P\,\sigma_{AB}^2\,\Omega_D(T^*)}
\quad [\text{m}^2\text{s}^{-1}],
\tag{5.25}
$$

with `T` in K, `M_A, M_B` in **kg·mol⁻¹**, `P` in **Pa**, `σ_{AB}` in **Å**, and
`D_{AB}` in **m²·s⁻¹**.

`chapman_enskog_diffusivity(T, M; P = 101325, σ = nothing, ε_k = nothing)`
(`effective_properties.jl`). It is not a property law, so it carries its own
temperature floor: it evaluates at `max(T, 100 K)` (`T_GAS_PROPERTY_FLOOR`),
which keeps `T^{3/2}` real at any trial temperature. Besides serving as the
first-principles estimate of a component's `λ` (§5.10.3), it supplies the gas
diffusivity of the blowing models of a `ConvectiveBC` when a species carries
neither an explicit diffusivity nor Fuller diffusion volumes
(`src/BoundaryConditions/blowing.jl`; §12.3.3). The prefactor and unit handling deserve
care, because the *textbook* Chapman–Enskog constant is `1.8583 × 10⁻³` in
**cgs-mixed** units (`T` in K, `M` in g·mol⁻¹, `P` in atm, `σ` in Å, `D` in
cm²·s⁻¹), stored as `CHAPMAN_ENSKOG_PREFACTOR = 1.8583e-3`
(`effective_properties.jl`). The implementation converts inputs and output to
SI explicitly (`effective_properties.jl`), and **every** conversion
contributes to the lumped prefactor — not just the `D` output conversion:

- `M`: kg·mol⁻¹ → g·mol⁻¹ (`×1000`). This appears *inside* `√(1/M)`, so it
  introduces a factor `1/√1000` in the diffusivity.
- `P`: Pa → atm (`/101325`). Because `P` is in the denominator, using SI `P` in
  Pa introduces a factor `×101325`.
- `D`: cm²·s⁻¹ → m²·s⁻¹ (`×1e-4`).

The net SI prefactor is therefore

$$
1.8583\times10^{-3}\;\times\;10^{-4}\;\times\;\frac{101325}{\sqrt{1000}}
\;\approx\; 5.95\times10^{-4},
$$

which is the constant quoted in the SI form of (5.25) (`M` in kg·mol⁻¹, `P` in
Pa). Folding *only* the `×1e-4` output conversion into the cgs prefactor would
give `1.8583 × 10⁻⁷`, but that value is valid **only** when the inputs are kept
in their cgs-mixed units (`M` in g·mol⁻¹, `P` in atm); pairing `1.8583 × 10⁻⁷`
with the SI symbols defined for (5.25) is dimensionally inconsistent and
under-predicts `D` by a factor `101325/√1000 ≈ 3.2 × 10³`. The physical scalings
are `D ∝ T^{3/2}` (through the explicit `T^{1.5}` and weak `Ω_D` dependence),
`D ∝ 1/P`, and `D ∝ √(1/M_A + 1/M_B)` (lighter species diffuse faster).

### 5.10.1 Combining rules and the collision integral

The binary interaction parameters use the standard Lennard-Jones combining rules
(`effective_properties.jl`):

$$
\sigma_{AB} = \tfrac12(\sigma_A + \sigma_B), \qquad
\varepsilon_{AB} = \sqrt{\varepsilon_A\,\varepsilon_B},
$$

and the **reduced temperature** is `T* = T/(ε_AB/k_B)`
(`effective_properties.jl`). The **diffusion collision integral** `Ω_D(T*)`
is evaluated with the Neufeld et al. (1972) correlation,

$$
\Omega_D(T^*) =
\frac{A}{(T^*)^B} + \frac{C}{e^{D T^*}} + \frac{E}{e^{F T^*}} + \frac{G}{e^{H T^*}},
\tag{5.26}
$$

`collision_integral_diffusion(T*)` (`effective_properties.jl`), with the
coefficients verified against the source (`effective_properties.jl`):

| | | | |
|---|---|---|---|
| `A = 1.06036` | `B = 0.15610` | `C = 0.19300` | `D = 0.47635` |
| `E = 1.03587` | `F = 1.52996` | `G = 1.76474` | `H = 3.89411` |

The fit is valid for `0.3 ≤ T* ≤ 400` with <0.1 % error.

> **Note (Neufeld G).** The coefficient `G` is **`1.76474`**, as in the code
> and the original Neufeld, Janzen & Aziz (1972) correlation.

### 5.10.2 Heuristic Lennard-Jones defaults and their unvalidated status

When the user does not supply `σ` or `ε/k_B` for a species, the code estimates
them from the molar mass via rough correlations anchored to small organics
(`effective_properties.jl`):

$$
\sigma_{\text{species}} \approx 3.5 + 0.5\log_{10}\!\Bigl(\max(M/0.016,\,1)\Bigr)\ [\text{Å}],
\qquad
\frac{\varepsilon}{k_B} \approx 150 + 100\,\frac{M}{0.044}\ [\text{K}],
\tag{5.27}
$$

with the air background fixed at `σ_air = 3.617 Å`, `ε/k_B|_air = 97 K`,
`M_air = 0.029 kg·mol⁻¹` (`effective_properties.jl`). These are **rough
estimates** suitable for organic vapours; they have **not** been validated
against tabulated LJ parameters for the specific multi-component gas mixtures that
arise in pyrolysis, and precise values should be supplied when available (§5.16).

### 5.10.3 Conversion to a gas-transfer coefficient

The model's diffusive flux is written in terms of a transfer coefficient `λ`
[m²·s⁻¹] with the convention `J = −λ ∂ξ/∂x` (the gas-density factor is folded
into `λ`), so for gas-in-gas diffusion `λ ≈ D` is the consistent
identification: assign the `chapman_enskog_diffusivity` value (optionally
porosity/tortuosity-corrected, and optionally wrapped in a temperature-dependent
property function) as the gaseous component's `λ`. There is no separate
conversion routine — the identification is one-to-one (§9 derives the precise
volume-fraction-gradient flux that `λ` multiplies).

---

## 5.11 No per-species transport scaling

The effective gas-transfer coefficient `λ_eff` of §5.6 is a property of the
**mixture**, and the gas flux kernels apply the same face value to every gas
species (§9.2). There is no kinetic-theory per-species scaling of the diffusive
flux: a light gas and a heavy gas in the same simulation diffuse identically.
Species-faithful magnitudes enter only through the values the user assigns to
the component `λ_j(T)` properties — for which the Chapman–Enskog estimate of
§5.10 is the first-principles starting point. (The nomenclature symbol `M_ref`
appears only inside the Chapman–Enskog combining rules of §5.10.)

---

## 5.12 Face averaging: distance-weighted harmonic mean

The cell-centered transport properties (`k`, `λ`, and the Darcy mobility `κ/μ`)
must be transferred to **faces** to form the finite-volume fluxes (§13). For a
conserved flux crossing a material interface, the correct face value is the one
that preserves flux continuity, which makes the property behave as a
**resistance in series**. With distances `d_L, d_R` from the two cell centers to
the face, the **distance-weighted harmonic mean** is

$$
\varphi_{\text{face}}
\;=\; \frac{d_L + d_R}{\,d_L/\varphi_L + d_R/\varphi_R\,} .
\tag{5.28}
$$

The derivation: the conductance of each half-cell along the flux path is
`φ/d`, the two resistances `d_L/φ_L` and `d_R/φ_R` add, and the effective face
property over the total span `d_L + d_R` is (5.28). On a uniform mesh
`d_L = d_R` and (5.28) reduces to the standard harmonic mean
`2 φ_L φ_R/(φ_L + φ_R)`. On a stretched mesh (e.g. boundary refinement) the
distance weighting keeps the steady flux through a piecewise-constant property
exact whatever the two spacings, where the unweighted mean would not.

The kernels evaluate (5.28) in a guarded product form,
`guarded_harmonic_mean(a, b, d_a, d_b, η)` (`effective_properties.jl`),

$$
\varphi_{\text{face}}
\;=\; \frac{(d_L + d_R)\,\varphi_L\,\varphi_R}
{d_L\,\varphi_R + d_R\,\varphi_L + \eta\,(d_L + d_R)} ,
\tag{5.29}
$$

which is smooth in every argument for non-negative values and positive
distances. Writing (5.28) as `H`, (5.29) is `H/(1 + η/H)`: it equals the
resistance mean to a relative `η/H`, and it goes to zero continuously, with a
finite slope, as either value vanishes, instead of being cut off below a
threshold (for `d_R = 3 d_L` the slope with respect to `φ_L` at `φ_L = 0` is
`(d_L + d_R)/d_L = 4`). The guard `η` carries the property's units and lies far
below its physical values: `FACE_MEAN_GUARD = 1e-20` for the conductivity and
the gas-transfer coefficient, `MOBILITY_GUARD = 1e-30` m²·Pa⁻¹·s⁻¹ for the Darcy
mobility (§5.7.4).

> **Implementation note.** The face means are formed inside the kernels, on the
> interior faces between two live cells only (`ConductionKernel`,
> `src/Kernels/conduction.jl`, for `k`; `GasTransportKernel`,
> `src/Kernels/gas_transport.jl`, for `λ` and the mobility); nothing is cached
> between evaluations. The distances are the center-to-face distances of the
> current geometry along the face normal (so the averaging follows a moving
> mesh), each smoothly floored at `10⁻⁸` of its reference value
> (`DISTANCE_FLOOR_FRACTION`), so an inverted trial mesh keeps finite fluxes. A
> boundary face takes no mean: its flux reads the boundary cell's own property
> (§7.1.3, §9, §12). On a mesh of more than one dimension the conduction face
> tensor is the face mean `k_f` times the distance-weighted mean of the two
> cells' dimensionless conductivity tensors (§7.1.3). The harmonic mean is
> applied to `λ` itself, **not** to the product `λρ` that a
> mass-fraction-gradient discretization `J = −λρ ∂Y/∂x` would average: the
> volume-fraction-gradient flux used by the gas transport (§9) carries the
> density dependence in the `ξT` driving potential and averages `λ` alone.

The same resistance-in-series averaging is used for the mobility at faces
(§5.7.4).

---

## 5.13 Cell properties in the state mapping

The effective properties of this chapter are evaluated **once per cell per
residual evaluation**, in the state mapping `derive!` (`src/States/derived.jl`),
the only reader of the state vector, and stored in the derived state
(`DerivedState`) that every kernel reads. Nothing is carried from one evaluation
to the next. The full mapping — masses, volumes, temperatures, the moving-mesh
and radiation unknowns — is the subject of §13; this section states where and
how the cell properties are evaluated.

### 5.13.1 The evaluation point

Every property law is evaluated at the cell's **positive-part concentrations**
and **trial temperature**,

$$
\xi^+_{j,i} = \frac{m^+_{j,i}}{V_{\text{eval},i}},\qquad
m^+ = \begin{cases} 0, & m\le 0,\\ m\,(3t^2 - 2t^3),\ t = m/\delta_i, & 0 < m < \delta_i,\\ m, & m \ge \delta_i,\end{cases}
\qquad
T_{\text{trial},i} = \operatorname{smooth\_clamp}(T_i,\ 200,\ 3000,\ 1\ \text{K}),
$$

with

- `δ_i = reltol · ε_m · m^s_i(0)` the cell's positive-part scale, equal to its
  mass absolute tolerance (`reltol` and `threshold` of the solve, `m^s_i(0)` the
  cell's initial condensed mass; `src/States/tolerances.jl`);
- `Ṽ_i = smooth_floor(V_i, 10⁻⁸ V_i(0))` the smoothly floored cell volume (`ds.V_eval`);
- `smooth_clamp(x, lo, hi, w)` constant outside `[lo, hi]`, the identity on
  `[lo + w, hi − w]`, and the cubic join `lo + w t²(2 − t)`, `t = (x − lo)/w`
  (mirrored at `hi`) between; `T_trial` is therefore `T_i` on every physical
  state (`trial_temperature`, `src/States/guards.jl`);
- `smooth_floor(x, a)` equal to `a` at and below `a`, to `x` at and above `2a`,
  and to the cubic `a(1 + 2t² − t³)`, `t = x/a − 1`, between.

Each guard is C¹ and is exactly the identity on the physical states the
tolerances can distinguish, so a trial state outside the physical domain still
has finite properties with finite derivatives, while the domain test reads the
unguarded quantities (§15.5.4). The kernels read the signed concentrations
`ξ_{j,i} = m_{j,i}/Ṽ_i` where a sign matters (the reaction rates, §6; the
gas fluxes, §9); the property laws see only `ξ⁺`.

### 5.13.2 The cell properties

| Field | Symbol | Units | Evaluated as | Filled |
|---|---|---|---|---|
| `cp_sum` | `C_i` | J·K⁻¹ | (5.10): `Σ_j m_{j,i} c_{p,j}(T_trial)`, smoothly floored at `ρc_{p,floor} V_i(0)` | always |
| `h` | `h_j(T_i)` | J·kg⁻¹ | (5.11), at the cell temperature itself | always |
| `k` | `k_eff` | W·m⁻¹·K⁻¹ | `effective_conductivity(ξ⁺, T_trial; skip_absent = Val(true))` (§5.5), smoothly floored at `10⁻⁶` | always |
| `ε` | `ε_eff` | – | `mixture_emissivity(ξ⁺, T_trial)` (5.22) | always |
| `a` | `a_eff` | – | `mixture_absorptivity(ξ⁺, T_trial)` (5.22b) | always |
| `α` | `α_eff` | m⁻¹ | `effective_absorption(ξ⁺)` (5.23) | `BEER_LAMBERT` |
| `λ` | `λ_eff` | m²·s⁻¹ | `effective_gas_transfer(ξ⁺, T_trial)` (5.17) | Fickian or Darcy–Fick transport |
| `φ` | `φ` | – | `smooth_floor(1 − skeletal_volume_fraction(ξ⁺, T_trial), 10⁻¹²)` (§5.3) | Darcy–Fick transport |
| `P` | `P` | Pa | `ideal_gas_pressure(ξ⁺, T_trial, φ)` (§9) | Darcy–Fick transport |
| `mob` | `κ/μ` | m²·Pa⁻¹·s⁻¹ | `effective_permeability(ξ⁺, T_trial, φ) / smooth_floor(μ(T_trial), 10⁻¹²)` (§5.7.4) | Darcy–Fick transport |

The gas transport of a problem (`NoGasTransport`, `FickTransport` or
`DarcyFickTransport`, chosen at setup, §5.6) and the radiation model are type
parameters of the setup, so the selection costs nothing at run time; the
buffers of the conditional fields exist for every problem and are filled only
when their consumers exist. The mixture density is not stored: its one
consumer, the film flux of a mass boundary condition, evaluates it at the
boundary cell (§5.2).

### 5.13.3 Buffers and element types

The derived state is allocated once per element type: with `Float64` elements
for the residual, and as a twin for every dual-number type the colored Jacobian
evaluates the residual at, built on first use (`buffers`,
`src/Assembly/definition.jl`; §7.7, §15.2). Its sizes are fixed for the whole
solve: a deactivated cell keeps its slot and is marked dead in `active`, so no
buffer is ever resized (§14.3.1). After warm-up a residual evaluation allocates
nothing (`test/assembly/test_residual_jacobian.jl` asserts
`@allocated residual!(du, u0, def, 0.0, ws) == 0`): all composition-dependent
dispatch is generated and unrolled at compile time over the component tuple,
every intermediate tuple is statically sized, and the Bruggeman Newton
iteration works on scalars.

---

## 5.14 Comparison to ThermaKin2Ds, Gpyro, and FDS

The mixing-rule choices map cleanly onto the established condensed-phase pyrolysis
codes. Mapping their notation onto ours on first use: ThermaKin's stoichiometric
`θ_i^j` is our `ν_{i,j}`, Gpyro's volume fraction `X_α` is our `v_j`, Gpyro's
mixture average `\barρ` corresponds to our swelling mixture density `ρ`, and the
FDS/Gpyro pre-exponential `Z` is our `A_i` (§2).

**Effective heat capacity / gas storage.** ThermaKin (Eq. 17) keeps the gas
sensible-storage term — it sums `Σ_j ξ_j c_j` over *all* components, gas included
— and Pyrolysis.jl does the same in (5.9). Gpyro's default thermal-equilibrium
form drops gas storage from `ρc_p` (its `FULL_QSG` option restores it) and FDS
carries no in-depth gas; both are approximations of (5.9) that Pyrolysis.jl
does not make (§5.4.2).

**Conductivity / emissivity basis.** ThermaKin determines surface properties from
the *solid* matrix: "when γ = 0, gaseous components do not contribute to the
thermal conductivity" (Eqs. 8–10). Pyrolysis.jl reproduces this with the
**swelling-weighted volume-fraction basis** (5.4) for both `k_eff` and `ε_eff`;
the `γ_j·ξ_j/ρ_j` weighting reduces exactly to ThermaKin's `γ = 0` rule for
canonical gases (zero weight) and condensed phases (full weight).

**Conductivity blend.** ThermaKin blends the series and parallel limits with a
weighting parameter — our `WEIGHTED` rule (5.14) with `β` is the direct analogue
(`k = β k_∥ + (1−β) k_series`), default `β = 0.5`. Gpyro homogenizes by
volume-fraction averaging (its `X_α`-weighted properties) with power-law
temperature dependence built into the component properties; our `PARALLEL` rule
(5.13, the volume-weighted arithmetic mean) is the volume-fraction average, and
arbitrary `k_j(T)` (including power laws) is supported through the property
functions (§4). The **BRUGGEMAN** effective-medium option (5.15) is an addition
beyond the volume-average/blend menu of ThermaKin and Gpyro, motivated by the
self-similar char microstructure.

**Permeability.** Gpyro uses a porosity-dependent permeability for its Darcian
gas transport; our **Carman–Kozeny** closure (5.20) is the canonical
porosity-based form (Kozeny constant 180), and is offered alongside component
mixing as a mutually exclusive alternative.

**Radiation properties.** FDS works with a solid absorption coefficient `κ_s`
(our `α`, §2) and a surface emissivity/absorptivity tied by Kirchhoff's law; our
`α_eff = Σ ξ_j α_j` (5.23) and `ε_eff` (5.22) provide exactly these. The surface
absorptivity `a_eff` (5.22b) equals `ε_eff` bit-for-bit unless a component sets
its own `surface_absorptivity`, which decouples the absorbed heater flux from the
emission (§5.8.2).

---

## 5.15 Mixing derivatives for the Jacobian

The colored Jacobian (§15) differentiates the whole residual — the state
mapping, every kernel and the row assembly — with forward-mode dual numbers on
the sparsity pattern of the mesh connectivity. No mixing rule carries a
hand-coded derivative, and no face flux or transmission chain is composed by
hand: every derivative follows the closure as the residual evaluates it, at the
same evaluation point `(ξ⁺, T_trial)` and through the same guards (§5.13).

### 5.15.1 Permeability and porosity derivatives

The pore fraction and the cell mobility are evaluated on dual numbers inside the
state mapping, so their derivatives track every closure: the
renormalized-solid-basis PARALLEL/SERIES/WEIGHTED mixing, including the
impermeable-phase series hard zero (where the derivative is exactly 0 on the
zero branch); the CARMAN_KOZENY chain through the smoothly floored pore fraction
(5.7); any user-supplied `κ_j(T)` (so `∂κ/∂T` is exact for a
temperature-dependent κ); and the viscosity curve with its smooth floor.
BRUGGEMAN is rejected for permeability (§5.7.3).

### 5.15.2 Conductivity, gas-transfer, and other derivatives

For `k_eff`, `λ_eff`, `ε_eff`, `a_eff`, and `α_eff`, the temperature dependence
enters through the component laws `k_j(T)`, `λ_j(T)`, `ε_j(T)` (§4.5) and the
composition dependence through the volume fractions (5.3)/(5.4) and the void
blend (for `α_eff`, through `ξ` directly). Because the properties are evaluated
at `ξ⁺ = m⁺/Ṽ`, the derivatives with respect to the masses carry the slope
of the positive part, and on a moving mesh the derivatives with respect to the
node positions arrive through the cell volume. The Bruggeman `k_eff` is
differentiated through its Newton iteration. The iteration starts from
`k_parallel`, carrying the partials of the parallel mean, and its convergence
test reads only the value; but the test follows the update, so at least one
dual-propagating Newton update runs before it is passed, and by the
implicit-function property each update brings the partials to those of the
implicit closure (5.15) to `O(tol)`. The one exception is the `|df| < 1e-20`
break before the first update, which requires degenerate inputs (all contrasts
≈ 0) that do not arise from valid materials.

The remaining non-smooth points of the mixing rules are thresholds and
branches, at which the derivative is that of the branch the value takes:

- the Carman–Kozeny cap `min(·, ℓ²)` (a kink of `dκ/dφ` at `φ = 0.9328`);
- the `10⁻²⁰` thresholds of the conductivity and gas-transfer series sums and
  the `10⁻³⁰` threshold of the permeability series sum, with the switch to an
  infinite series resistance once an active phase falls below its threshold;
- the empty-cell guards of the volume fractions and of the mixture density
  (`< 10⁻²⁰`);
- the Bruggeman branch selections (active and insulating phases, the step
  limiter, the collapse below `10⁻³⁰`), locally constant at a converged
  interior solution.

Together with the smooth `tanh` temperature gates and the C¹ depletion limiter
(§6), the C¹ guards of the state mapping (§5.13) and the property bracket
(§4.4), this is what keeps the colored Jacobian the exact derivative of the
residual — it is checked against a dense forward-mode Jacobian (§15.13) — and
the finite-difference sensitivities of §17 well-defined.

---

## 5.16 Limitations and open issues

The following modelling limitations of the effective-property subsystem are
recorded from the source; several are flagged here for the
verification and validation discussion of §18.

1. **Single-curve gas viscosity.** The pore-gas viscosity is one curve per
   material, a function of the temperature alone (§5.9). The *material*
   dependence is resolved — a material supplies the blend correlation for the
   gas it produces — but the *composition* dependence within a material is not:
   one curve averages over the pore gas shifting from background air to
   volatiles as pyrolysis proceeds. A Wilke mixing rule over per-component
   μⱼ(T) would resolve that transition. It would be evaluated where the
   viscosity is evaluated, in the pore-gas properties of the state mapping
   (`_gas_properties!`, `src/States/derived.jl`), from `ξ⁺` and `T_trial`, and
   the colored Jacobian would carry its composition derivatives with no further
   change. A Darcy–Fick test with more than one gas component would be needed to
   exercise it, since a one-species mixing rule is the identity.

2. **Unvalidated Chapman–Enskog LJ defaults.** The heuristic `σ(M)` and
   `ε/k_B(M)` correlations (5.27) are rough estimates anchored to small organics
   and have not been validated against literature LJ parameters for specific
   pyrolysis gas species. Supply explicit `σ`, `ε_k` where known.

3. **Permeability microstructure evolution.** Real materials may exhibit
   permeability changes with thermal stress and cracking beyond what the
   `κ_j(T)` property functions and porosity-based closures can represent;
   the colored Jacobian differentiates whatever closure is supplied
   (including `κ_j(T)` and the Carman–Kozeny porosity chain), but the closure
   itself remains the modelling limitation.

4. **Intrinsic-density temperature dependence.** Intrinsic (skeletal) density
   shares the same functional form as bulk density. In reality the cell-wall
   material may have a different thermal-expansion behaviour, which would shift
   the porosity (5.7) and hence the pressure and Carman–Kozeny permeability.

5. **Bruggeman near percolation.** The Newton iteration (5.15)–(5.16) is
   safeguarded by the 90 % step limiter, the insulating-phase exact limit, and
   the below-threshold collapse return (§5.5.3). If a continuous phase
   approaches zero volume while its conductivity diverges relative to the
   others, convergence can still stall; the step limiter is heuristic, not a
   percolation-aware bound.

6. **Carman–Kozeny microstructure assumption.** The `φ³/(1−φ)²` closure with
   Kozeny constant 180 is empirical for a packed bed of spheres. Real porous char
   may have a fractal or anisotropic pore structure for which a different closure
   (or a measured `κ(φ)` curve) would be more accurate.

The matching subsystem discussions are: the energy assembly that consumes the
heat-capacity sum, the sensible enthalpies and the face conductivities — §7
(Heat Conduction, Gas Convection, and Energy Assembly); the gas fluxes that
consume `λ_eff`, `λ_j`, `κ_eff`, `μ`, and `φ` — §9 (Gas-Phase Mass Transport and
Pressure); the radiation models that consume `α_eff`, `ε_eff` and `a_eff` — §8
(Thermal Radiation); the swelling/volume-change use of `ρ`, `γ_j` — §10 (Volume
Change, Swelling, Intumescence, and Lateral Shrinkage); and the differentiation
of every mixing rule — §15.
