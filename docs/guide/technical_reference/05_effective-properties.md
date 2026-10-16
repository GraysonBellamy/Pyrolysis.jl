# 5. Effective (Mixture) Properties and Mixing Rules

This chapter derives, in full, how Pyrolysis.jl collapses the per-component
thermophysical and transport properties of a multi-phase pyrolyzing solid into
the *single* effective property each governing equation needs. It covers the two
averaging bases (mass fractions and volume fractions, and the two flavours of the
latter), the swelling-aware mixture density, the intrinsic-density porosity, the
mixture effective heat capacity, the four
conductivity mixing rules (parallel / series / weighted / Bruggeman) together
with the Bruggeman Newton solve, the effective gas-transfer and permeability
closures (component mixing versus Carman–Kozeny), the radiative properties
(emissivity and absorption), the gas viscosity (Sutherland) and Chapman–Enskog
diffusivity, the absence of per-species transport scaling, the distance-weighted
harmonic face averaging, the property cache, and the mixing derivatives the
Jacobian consumes. Every equation, constant, guard, threshold, and sign below has
been confirmed line-by-line against `src/Properties/mixture_rules.jl`,
`src/Properties/effective_properties.jl`, `src/Properties/property_cache.jl`,
`src/Materials/components.jl`, and `src/Materials/materials.jl`; where the
implementation departs from the textbook form (most importantly which basis each
property uses, the harmonic-mean face average, and which mixing rules are
rejected for which channel), **the code is authoritative** and the deviation is
called out explicitly.

Throughout this chapter `j` is the **component (species) index** running over the
`N_C` components of a `Material`, and `i` is the **cell index** running over the
`n` finite volumes (overload resolution I1, §2). The primary species state is the
**mass concentration** `ξ_j` (kg·m⁻³), not the mass fraction (§3, Governing
Equations). Two density symbols are used with care: the **bulk** (pure-phase)
density `ρ_j` drives all mixing rules and the mixture density, while the
**intrinsic** (skeletal, cell-wall) density `ρ_{i,j}` appears *only* in the
porosity formula (overload resolution A7, §2). The three binding sign conventions
of §2 are assumed without restatement.

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
storage term weights by the raw concentrations `ξ_j` directly (§5.4). They are
*not* the right weights for transport properties, which act over volumes.

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
> fractions while the matrix volume fraction $\varphi_m = \sum_j \gamma_j\xi_j/\rho_j$
> is above $10^{-3}$, the all-phase occupancy (5.3) below $10^{-4}$, and a cubic
> smoothstep blend of the two in between. A fully converted cell thus degrades
> to the gas-mixture conductivity instead of 0 (which would zero both adjacent
> harmonic face means and, via Kirchhoff, its absorption), and because the
> matrix fractions are normalised by $\varphi_m$ their derivatives grow as
> $1/\varphi_m^2$ — the smooth blend keeps $\partial k_{\rm eff}/\partial\xi$
> bounded while a cell empties on a fixed mesh, which parks only its top live
> cell (§14.3.2).

The three bases and their consumers are summarized below; this mapping is the
single most error-prone part of effective-property bookkeeping, so it is stated
once authoritatively from the code in `update_properties!`
(`property_cache.jl`):

| Basis | Symbol | Built by | Consumed by |
|---|---|---|---|
| Mass fraction | `Y_j` | inline in `_specific_vol_with_swelling_direct` | swelling mixture density (§5.2; the heat-capacity storage term weights by raw `ξ_j`, §5.4) |
| All-phase volume fraction | `v_j` | `volume_fractions` | `effective_gas_transfer` |
| Swelling-weighted volume fraction | `v_j^solid` | `solid_volume_fractions` (via `mixing_volume_fractions`, which adds the depleted-cell all-phase fallback) | `effective_conductivity`, `effective_emissivity`, `effective_absorptivity`; `effective_permeability` reads `solid_volume_fractions` directly |

Note in particular that **gas transfer uses the all-phase basis `v_j`** while
**conductivity and emissivity use the solid basis `v_j^solid`** — this is what
`update_properties!` does (`cache.λ[i] = _gas_transfer_from_v(v_full, …)` versus
`cache.k[i] = _conductivity_from_v(v_solid, …)`).

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
constitutive density.

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
\phi \;=\; \operatorname{clamp}\!\left(
1 - \sum_{j\in\{\text{solid,liquid}\}} \frac{\xi_j}{\rho_{i,j}(T)},
\;0,\;1\right),
\tag{5.7}
$$

implemented as `porosity(ξ, T, material)` →
`_porosity_void_frac` (`mixture_rules.jl`). The sum runs over **solid and
liquid** components only — the generated code branches on
`is_solid(c) || is_liquid(c)`, and each term is taken only if
`ρ_intr > 0`. Gases are *excluded from the sum* because they occupy the pore
space being measured, not the skeleton. The result is clamped to `[0,1]`; a
clamp to 0 occurs only in the (numerical-error) case where the computed solid
volume exceeds the bulk volume.

The distinction between bulk and intrinsic density is essential here. The
*intrinsic* density is the density of the cell-wall material itself, ignoring any
micro-pores inside the pure phase; the *bulk* density is the apparent density of
that phase as it sits in the composite. For a porous solid these differ, and
using the intrinsic density in (5.7) correctly recovers the void fraction. When a
component is constructed **without** an explicit intrinsic density,
`get_intrinsic_density` falls back to the bulk density (`components.jl`),
which makes the summand `ξ_j/ρ_j`: the pure phase is taken to have no internal
pores. For a gas, `get_intrinsic_density(::GaseousComponent,
::Real) = 0.0`, consistent with its exclusion from the sum.

Porosity feeds two downstream models:

1. the **ideal-gas pressure closure** `P = N(ξ,T)/φ`, with
   `N = Σ_{g} ξ_g R_g T / M_g` the partial-pressure numerator (§9, Gas-Phase
   Mass Transport and Pressure); and
2. the **Carman–Kozeny permeability** closure (§5.7).

### 5.3.2 Porosity derivatives

Differentiating (5.7) term by term, and recognizing that the clamp is inactive in
the physical interior `0 < φ < 1` (where derivatives are needed), gives

$$
\frac{\partial \phi}{\partial \xi_j}
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
These derivatives are not hand-coded: the Jacobian path that needs them
(the `DarcyTransportOperator`'s Carman–Kozeny permeability and pressure
chains, §5.15, §15.8.3) differentiates
`porosity` itself with forward-mode `Dual`s, which is clamp-consistent by
construction — at the `clamp(·, 0, 1)` bounds the propagated derivative is 0,
exactly matching the residual's clamped value.

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

implemented as `effective_heat_capacity(ξ, T, material)` →
`_sum_heat_capacity` (`effective_properties.jl`). The generated helper adds
`ξ[j] * comps[j].heat_capacity(T)` for every component with no phase branch,
so **gas sensible storage is part of the left-hand side of the energy
equation**. The units check: `[kg·m⁻³] · [J·kg⁻¹·K⁻¹] = [J·m⁻³·K⁻¹]`.

The result is floored at `RHO_CP_FLOOR = 10³` J·m⁻³·K⁻¹
(`src/Materials/constants.jl`), the thermal mass of pore gas at ambient
(ρ_g c_{p,g} ≈ 1.2 × 1000; Gpyro floors the same quantity). Because every
component is summed, ρc_p reaches the floor only in a cell emptied of all
mass — a burned-out cell whose pore gas has also left; without it the 1/ρc_p
factor in dT/dt would blow up there. Virgin solids carry ρc_p ~ 10⁵–10⁶
J·m⁻³·K⁻¹, so the floor never binds in a cell that holds material. The floor
is applied *inside* `effective_heat_capacity`, so every consumer — the
residual, the property cache and the
structured-Jacobian operators that extract ∂ρc_p/∂(T, ξ) through Duals —
sees the identical clamped value and the consistent (zero) derivative at the
floor.

### 5.4.2 Why the gas is stored

The one-temperature model puts the pore gas at the matrix temperature (§3.1.2),
so the gas stores sensible heat `ξ_g c_{p,g} ∂T/∂t` exactly like every other
component. Its concentration `ξ_g` is an integrated state, so the term costs
nothing; dropping it — Gpyro's default *quasi-steady-gas approximation*, which
its `FULL_QSG` option reverses, and the implicit convention of FDS, which has
no in-depth gas — would be an approximation with no benefit. The in-pore gas
heat capacity is roughly `φ ρ_g c_{p,g}` against the matrix `(1−φ) ρ_s c_{p,s}`;
at 1 atm and pyrolysis temperatures the ratio is a few tenths of a percent and
it scales with pressure. On white-pine cone decks the term is at most about
half a percent of ρc_p and shifts the response by a uniform time lag of about
a tenth of a percent — small, but a thousand times above the solver's noise
floor and therefore worth keeping exact.

The gas *flux* relative to the skeleton enters the energy equation only
through the advective source `S_conv = −Σ_g c_{p,g}(T) J̄_g (∂T/∂z)` (§7). This is
the ThermaKin Eq. 17 convention term for term (§3.3.2).

No separate mass-averaged specific heat is computed or cached; consumers that
need `c_p` per unit mass divide `(ρc_p)_eff` by `ρ_tot`.

---

## 5.5 Thermal conductivity mixing

The effective thermal conductivity is computed from the swelling-weighted
volume fractions `v_j^solid` (5.4) and the per-component conductivities
`k_j(T)` under one of four mixing rules selected by `material.mixing.k.rule`:
`effective_conductivity(ξ, T, material)` builds
`v = mixing_volume_fractions(…)` (the (5.4) weights, with the depleted-cell
all-phase fallback of §5.1.3) then dispatches to `_conductivity_from_v`
(`effective_properties.jl`). The default rule is `WEIGHTED` with
`β = 0.5` (`materials.jl`).

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
> `_conductivity_sums(v, T, comps, Val(NC))`
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
`bruggeman_effective_property(v, k; tol, maxiter, k_init)`
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
- **Warm start.** The Newton seed is `k_init` when a positive `k_init` is supplied
  (`k_init > 1e-20`), otherwise the parallel mean `k_parallel`
  (`effective_properties.jl`). In production `k_init` is the *previous time
  step's converged* `k_eff` for the cell (passed by `update_properties!` as
  `cache.k[i]`, `property_cache.jl`), so on smooth time-stepping the solve
  needs only 1–3 iterations.
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
> (`effective_properties.jl`).

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
dispatches to `_gas_transfer_from_v` (`effective_properties.jl`). This is
confirmed in `update_properties!` where the cached gas-transfer property is
`cache.λ[i] = _gas_transfer_from_v(v_full, T_i, material)`
(`property_cache.jl`). The supported rules mirror conductivity:

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
under SERIES), and the face harmonic mean (5.29) vanishes when either side does.
Gas transport then stops where the gas fraction stops: the gas field of a
gas-producing material ends in a sharp front, and a Newton iterate that drives
the front cell's gas concentration below zero closes both of its faces and
flattens its residual row. The implicit solve fails repeatedly on that front.
On the Landau ablation case V-CPL-02 a solid coefficient of
$10^{-6}$ m²/s halves the step count and removes more than half of the
Newton failures, while the measured regression rate moves by $10^{-5}$
relative.
`validate_problem`, which `solve` calls, therefore warns when a material with a
gas-producing reaction and Fickian transport (no permeable component) has a
condensed component with $\lambda_j = 0$
(`warn_degenerate_gas_transport`, `src/Problem/problem_validation.jl`).

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
(`effective_permeability(ξ, T, material)`, `effective_properties.jl`).
The default rule for `κ` is **WEIGHTED** with `β = 0.5` (`materials.jl`).

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
\kappa \;=\; \min\!\left(\frac{\ell^2\,\phi^3}{180\,(1-\phi)^2},\;\ell^2\right),
\tag{5.20}
$$

`_carman_kozeny_permeability(φ, ℓ)` (`effective_properties.jl`), where
`ℓ` is the characteristic particle/pore diameter [m] supplied through the
`MixingSpec` `length_scale` and `φ` is the *canonical* porosity (5.7). The
constant **180** is the Kozeny constant for a packed bed of spheres. The
`φ³/(1−φ)²` dependence captures the steep growth of permeability as porosity
rises during virgin → char conversion: as `φ → 0` (fully solid), `κ → 0`
(impermeable). As `φ → 1` the packed-bed closure loses validity and diverges,
so `κ` is capped at the physical free-flow scale `ℓ²` (an open channel of width
`ℓ` has `κ ~ ℓ²/32`; `ℓ²` is a generous upper bound). The cap binds for
`φ ≳ 0.93`, where `180(1−φ)² < φ³`; a tiny `1e-30` denominator guard remains
solely to avoid `Inf`/`NaN` at exactly `φ = 1` before the `min` resolves
(`effective_properties.jl`).

Because Carman–Kozeny needs the porosity (which depends on `ξ` through the
intrinsic densities), it must be called through `effective_permeability(ξ, T,
material)` — calling the cache-aware `_permeability_from_v` with this rule throws
an `ArgumentError` explaining that concentrations are required
(`effective_properties.jl`). `MixingSpec(CARMAN_KOZENY; length_scale=ℓ)`
fails fast at construction if `ℓ` is not a positive, finite number
(`materials.jl`).

The closure derivative needed for the Jacobian (§5.15) is not hand-coded: the
`DarcyTransportOperator` (§15) differentiates `effective_permeability` itself by
forward-mode AD, so the Carman–Kozeny chain
(through the porosity's composition dependence) is picked up automatically. On
the capped branch the `min` gives the value-consistent zero `dκ/dφ` — value and
derivative saturate at the same point.

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
applied to the **mobility** κ/μ rather than to κ alone: the production path
(`update_darcy_properties!`, src/Discretization/finite_volume.jl) forms each
cell's mobility mobᵢ = κᵢ/μ(Tᵢ) and combines the two sides with the
distance-weighted harmonic mean (§5.12),

    mob_face = (d_L + d_R) / (d_L/mob_L + d_R/mob_R),

with mob_face = 0 exactly when either side's mobility is below `1e-30` — an
impermeable cell is an infinite series resistance, so gas cannot cross into it.
The Jacobian operator's `_dt_mobility_face` applies the identical rule.

---

## 5.8 Radiative properties: emissivity and absorption

### 5.8.1 Effective emissivity

The surface emissivity is a **solid**-volume-fraction average of the component
emissivities, each evaluated at the cell temperature, clamped to the physical
range:

$$
\varepsilon_{\text{eff}}(T,\xi) \;=\; \operatorname{clamp}\!\left(
\sum_j v_j^{\text{solid}}\,\varepsilon_j(T),\ 0,\ 1\right),
\tag{5.22}
$$

`effective_emissivity(ξ, T, material)` builds the solid basis and calls
`_emissivity_from_v(v, T, material)`, which evaluates `get_emissivity(comp, T)`
per component and applies `clamp(ε_eff, 0, 1)` (`effective_properties.jl`). The
solid basis is used for the same physical reason as conductivity (§5.1.3): the
radiating surface is the solid/liquid material, not trapped non-swelling gas.
Component emissivities are property functions (`ε` accepts the same forms as
`k`: number, linear tuple, callable, property function) with defaults of `0.9`
for `SolidComponent`, `0.95` for `LiquidComponent` (condensed liquids are strong
emitters — water ≈ 0.96), and `0.0` for `GaseousComponent`
(`src/Materials/components.jl`). Because the sum is evaluated generically over
the scalar type, a `ForwardDiff.Dual` temperature returns
$\partial\varepsilon_{\text{eff}}/\partial T$, which the surface-emission
(§12) and Beer–Lambert (§8) Jacobians pick up without operator-level
bookkeeping. The effective emissivity is the surface property used by the
re-radiation boundary condition (§12) and as the gray-body factor in the
volumetric P1 emission source of the experimental P1 model (§8), which `solve`
does not accept.

### 5.8.2 Effective surface absorptivity

The absorptivity for *incident* (heater-spectrum) radiation is mixed with the
same weights,

$$
a_{\text{eff}}(T,\xi) \;=\; \operatorname{clamp}\!\left(
\sum_j v_j^{\text{solid}}\,a_j,\ 0,\ 1\right),\qquad
a_j = \begin{cases}
\texttt{surface\_absorptivity}_j & \text{if set}\\
\varepsilon_j(T) & \text{otherwise (Kirchhoff)}
\end{cases}
\tag{5.22b}
$$

via `effective_absorptivity(ξ, T, material)` → `_absorptivity_from_v` →
`get_absorptivity(comp, T)`. The Kirchhoff branch resolves at compile time
(`Nothing` vs `Float64` is a type parameter of the component), and for a
material in which no component sets a surface absorptivity the two generated
sums perform the same floating-point operations in the same order, so
$a_{\text{eff}} \equiv \varepsilon_{\text{eff}}$ bit-for-bit — a material that sets no absorptivity
absorbs exactly as Kirchhoff's law predicts. A component with an explicit
$a_j$ decouples the absorbed heater flux from the temperature-dependent
emission: an emissivity that falls with temperature (integrating-sphere data for
virgin wood, for instance) does not drag the absorbed cone flux down with it.
$a_{\text{eff}}$ is consumed by a `RadiativeFluxBC` with
`absorptivity = nothing` (§12.2.6) and as the Beer–Lambert surface
transmissivity (§8.3); surface emission always uses
$\varepsilon_{\text{eff}}$.

### 5.8.3 Effective absorption coefficient

The volumetric (extinction) absorption coefficient used by the in-depth radiation
models (Beer–Lambert, and the experimental P1 model, §8) is a **concentration**-weighted sum of the
component *mass* absorption coefficients `α_j` [m²·kg⁻¹]:

$$
\alpha_{\text{eff}} \;=\; \sum_j \xi_j\, \alpha_j \quad [\text{m}^{-1}],
\tag{5.23}
$$

`effective_absorption(ξ, material)` → `_sum_absorption`
(`effective_properties.jl`). The dimensional bookkeeping is the key
point: a *mass* absorption coefficient `α_j` [m²·kg⁻¹] times a mass concentration
`ξ_j` [kg·m⁻³] yields a *volumetric* absorption coefficient [m⁻¹], which is what
the Beer–Lambert optical thickness `τ = α Δz` and the P1 diffusion coefficient
`D_rad = 1/(3α)` require. Because the weight is the raw concentration `ξ_j` (not a
normalized fraction), `α_eff` correctly scales down as material is consumed — a
thinning, depleting layer becomes optically thinner. Component absorption
coefficients default to `0.0` (no absorption) and carry no temperature
dependence.

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
pore-gas viscosity curve in the `gas_viscosity` field, and the Darcy path reads
it through `effective_gas_viscosity(T, material)`, which evaluates
`material.gas_viscosity(T)` on the property bracket and floors the result at
`MU_FLOOR`. A material that says nothing
carries `SutherlandProperty()` — the air triple above — so this reduces to
(5.24).

The field accepts anything `to_property` understands, so a fitted correlation
drops in directly:

```julia
# Cantera mixture-averaged fit for wood pyrolyzates, ascending coefficients,
# valid 280–1200 K.
Material(name = :WhitePine, components = (...), reactions = (...),
         gas_viscosity = PolynomialProperty(mu_coeffs))
```

**Why the guards matter.** μ enters the model *only* as the Darcy mobility
κ/μ (§9), never on its own. An error in μ is therefore exactly the reciprocal
error in permeability — and a μ that reaches zero or changes sign does not
degrade the answer, it reverses or blows up the gas velocity. Hence:

- the curve is a property law, so it evaluates on the property bracket
  $[200, 3000]$ K (§4.4): diverging Newton trial states cannot kill the solve
  on `T^{3/2}` or evaluate a fitted polynomial absurdly far outside its range
  (the Chapman–Enskog diffusivity formula, which is not a property law, keeps
  its own 100 K floor `T_GAS_PROPERTY_FLOOR`);
- the result is floored at `MU_FLOOR` (1e-12 Pa·s, seven orders below any real
  gas) so an extrapolated fit cannot go non-positive;
- `Material`'s constructor probes the curve over 200–1800 K and refuses one that
  is already non-positive or non-finite there, in the spirit of `MixingSpec`'s
  `CARMAN_KOZENY`/`length_scale` check.

Because the bracket belongs to every property law and the floor is applied in
`effective_gas_viscosity(T, material)`, every curve a material can carry is
protected, not just the built-in Sutherland form.

**Accuracy.** Air is within ~2 % for an N₂-swept pore space, but overstates the
viscosity of wood pyrolyzates by roughly 1.7× at 300 K falling to 1.3× at
1000 K — a bias in one direction everywhere, growing as the pore space fills
with volatiles.

**Composition dependence.** `mixture_gas_viscosity(ξ, T, material)` takes `ξ`
and every Darcy call site passes it, but the function returns the
material's own curve regardless: the model carries **one curve per material**, an average over the
composition that material is expected to produce. It cannot resolve the pore gas
shifting from background air to volatiles as pyrolysis proceeds. A
composition-dependent rule (Wilke) would be introduced by changing that one
function — but see §5.16 item 1 for why that is not a local change.

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

`chapman_enskog_diffusivity(T, M; P, σ, ε_k)`
(`effective_properties.jl`). The prefactor and unit handling deserve
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

All cell-centered transport properties (`k`, `λ`, `κ`) must be transferred to
**faces** to form the finite-volume fluxes (§13). For a conserved flux crossing a
material interface, the correct face value is the one that preserves flux
continuity, which makes the property behave as a **resistance in series**. With
distances `d_L, d_R` from the two cell centers to the face, the
**distance-weighted harmonic mean** is

$$
\varphi_{\text{face}}
\;=\; \frac{d_L + d_R}{\,d_L/\varphi_L + d_R/\varphi_R\,},
\tag{5.29}
$$

`harmonic_mean_weighted(a, b, d_a, d_b)`
(`effective_properties.jl`). The derivation: the conductance of each
half-cell along the flux path is `φ/d`, the two resistances `d_L/φ_L` and
`d_R/φ_R` add, and the effective face property over the total span `d_L + d_R` is
(5.29). On a uniform mesh `d_L = d_R` and (5.29) reduces to the standard harmonic
mean `2 φ_L φ_R/(φ_L + φ_R)`. On a stretched mesh (e.g. boundary refinement) the
distance weighting corrects the averaging and improves boundary fluxes by ~1–5 %.

> **Implementation note.** `update_properties!` applies (5.29) to fill `k_face`
> and `λ_face` for every interior face; at a boundary face (where one
> `cell_id == 0`) it simply copies the adjacent interior cell's value
> (`property_cache.jl`). The distances are computed from the live state
> positions `d_L = |z_face − z^c_L|`, `d_R = |z^c_R − z_face|`, so the averaging
> is correct under ALE mesh motion. A small-value guard returns 0 if either input
> is below `1e-20` (`effective_properties.jl`), preventing `Inf`/`NaN`
> from a vanished cell. The harmonic mean is applied to `λ` itself, **not** to
> the product `λρ` that a mass-fraction-gradient discretization
> `J = −λρ ∂Y/∂x` would average: the volume-fraction-gradient flux used by the
> gas transport (§9) carries the density dependence in the `ξT` driving
> potential and averages `λ` alone.

The same resistance-in-series averaging is used for permeability at faces
(§5.7.4) and for the P1 radiation diffusion coefficient `D_rad` at faces (§8).

---

## 5.13 The property cache

Effective properties are not recomputed inside every inner flux loop; they are
evaluated **once per cell per residual** and stored in a `PropertyCache`, then
read by the conduction, transport, radiation, and boundary kernels.

### 5.13.1 Layout

`PropertyCache{Tu, DarcyT}` (`property_cache.jl`) is parameterized on the
live scalar type `Tu` so the entire cache holds `ForwardDiff.Dual` values
unchanged during automatic differentiation (§15), and on `DarcyT` so the
optional Darcy sub-cache is either present (`DarcyCache{Tu}`) or `Nothing`. The
per-cell arrays (length `n`) are

| Field | Symbol | Units | Built from |
|---|---|---|---|
| `ρ` | mixture density | kg·m⁻³ | `mixture_density` (5.6) |
| `ρcp` | `(ρc_p)_eff` | J·m⁻³·K⁻¹ | `effective_heat_capacity` (5.9) |
| `k` | `k_eff` | W·m⁻¹·K⁻¹ | `_conductivity_from_v` (§5.5) |
| `λ` | `λ_eff` | m²·s⁻¹ | `_gas_transfer_from_v` (§5.6) |
| `ε` | `ε_eff` | – | `_emissivity_from_v` (5.22) |
| `a` | `a_eff` | – | `_absorptivity_from_v` (5.22b) |
| `α` | `α_eff` | m⁻¹ | `effective_absorption` (5.23) |

and the per-face arrays (length `n+1`) are `k_face` and `λ_face` (5.29). The
**Darcy sub-cache** `DarcyCache{Tu}` (`property_cache.jl`) is allocated only
for pressure-driven modes (`with_darcy = true`) and adds per-cell `P, φ, κ, μ` and
per-face `v_gas, ∇P` — keeping the four extra per-cell and two per-face arrays out
of memory for the common Fickian-only configuration. Note `ρcp` caches the
mixture `(ρc_p)_eff` of (5.9) — the storage term the energy
RHS divides by (§7).

### 5.13.2 The update kernel

`update_properties!(cache, state, mesh, material)`
(`property_cache.jl`) is the single entry point. It (i) resizes the cache
if its length differs from the state's; (ii) loops over the live cells
`1:n_active_cells` (parked ghost slots, §14.3, are left untouched), building `v_full =
volume_fractions(ξ_i, T_i, …)` and `v_solid = mixing_volume_fractions(ξ_i, T_i,
…)` (the swelling-weighted basis with the depleted-cell fallback, §5.1.3)
once each and threading them into the effective-property helpers; and (iii)
loops over faces applying the distance-weighted harmonic mean. The cell loop is
exactly:

```julia
cache.ρ[i]   = mixture_density(ξ_i, T_i, material)
cache.ρcp[i] = effective_heat_capacity(ξ_i, T_i, material)
cache.k[i]   = _conductivity_from_v(v_solid, T_i, material; k_init = k_prev)
cache.λ[i]   = _gas_transfer_from_v(v_full, T_i, material)
cache.ε[i]   = _emissivity_from_v(v_solid, T_i, material)
cache.a[i]   = _absorptivity_from_v(v_solid, T_i, material)
cache.α[i]   = effective_absorption(ξ_i, material)
```

This is the authoritative statement of which basis each property uses:
**conductivity, emissivity, and surface absorptivity from `v_solid`, gas
transfer from `v_full`**
(§5.1.3). The previous step's conductivity `k_prev = cache.k[i]` is passed as the
Bruggeman warm-start `k_init` (§5.5.3). After sizing, the loop is **allocation
free**: all composition-dependent dispatch is generated (compile-time unrolled),
all intermediate tuples are statically sized, and the Bruggeman Newton iteration
reuses scalars.

`resize_cache!` grows all arrays (and the Darcy sub-cache, if
present) in place after a mesh change.

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
`α_eff = Σ ξ_j α_j` (5.23) and `ε_eff` (5.22) provide exactly these, with the
`α = ε_eff` Kirchhoff identity used at the boundary (§12).

---

## 5.15 Mixing derivatives for the Jacobian

The structured Jacobian (§15) needs the sensitivity of each effective property to
the state `(T, ξ)`. No mixing rule carries a hand-coded derivative: the Jacobian
operators evaluate the mixing functions on `ForwardDiff.Dual` arguments, so every
derivative follows the closure as implemented.

### 5.15.1 Permeability and porosity derivatives

The `DarcyTransportOperator` (§15) differentiates `effective_permeability` and
`porosity` per cell with a single stack-allocated `Dual` pass over `(T, ξ)`
(`NC+1` partials) and feeds the results into its analytical Darcy-flux Jacobian.
The mixing rule therefore has a single source of truth: the derivatives track
every closure — the renormalized-solid-basis PARALLEL/SERIES/WEIGHTED mixing
including the impermeable-phase series hard-zero (where the derivative is
exactly 0 on the zero branch), the CARMAN_KOZENY porosity chain, and any
user-supplied `κ_j(T)` (so `∂κ/∂T` is exact for temperature-dependent κ).
BRUGGEMAN is rejected for permeability (§5.7.3).

### 5.15.2 Conductivity, gas-transfer, and other derivatives

For `k_eff`, `λ_eff`, `ε_eff`, `a_eff`, and `α_eff`, the temperature dependence
enters through the component laws `k_j(T)`, `λ_j(T)`, `ε_j(T)` (§4.5) and the
composition dependence through the volume fractions (5.3)/(5.4) (for `α_eff`,
through ξ directly). The conduction, gas-transport and Beer–Lambert operators
evaluate these mixing functions once per cell on `Dual` arguments over `(T, ξ)`
and compose the face fluxes or the transmission chain by hand in their
analytical Jacobians; the reaction-source and boundary operators run ForwardDiff
over stencil-local kernels (§15). The Bruggeman `k_eff` is differentiated by AD
through its Newton solve (started from `k_parallel` in the Jacobian, warm-started
from the cached `k` in `update_properties!`),
which is differentiable because every operation in the iteration is smooth (the
step limiter is inactive and the per-term branch selection is locally constant
at a converged interior solution).
This AD-cleanliness — together with the smooth `tanh` kinetics gates (§6) and the
clamped property laws (§4.4) — is what keeps the structured Jacobian (§15) exact
and the finite-difference sensitivities of §17 well-defined.

**Warm-start caveat.** The Bruggeman iteration's convergence test is on the
*value* only, and in an AD pass the Dual-typed cache means the warm start
`k_init` carries the *previous* step's partials. This is benign in the normal
path: at least one Dual-propagating Newton update always executes before the
tolerance break, and by the implicit-function property a single update
annihilates stale partials to O(tol) — the returned derivatives are those of
the implicit closure, not of the stale seed. The one failure mode is the
pre-update `|df| < 1e-20` degenerate-input break, which can return before any
update has run and thus pass the stale partials through; it requires
degenerate inputs (all contrasts ≈ 0) that do not arise from valid materials.

---

## 5.16 Limitations and open issues

The following modelling limitations of the effective-property subsystem are
recorded from the source; several are flagged here for the
verification and validation discussion of §18.

1. **Single-value mixture gas viscosity.** `mixture_gas_viscosity` ignores
   composition and returns the material's single viscosity curve. The
   *material* dependence is resolved — a material supplies the blend correlation
   for the gas it produces — but the *composition* dependence within a material
   is not: one curve averages over the pore gas shifting from background air to
   volatiles as pyrolysis proceeds. A Wilke mixing rule over per-component
   μⱼ(T) would resolve that transition.

   Two things make this more than a one-function change, and both are worth
   knowing before attempting it:

   - The production Darcy Jacobian is the **analytical override** in
     `Jacobian/operators/darcy_transport.jl`, not the AD path (which is retained
     only as the rtol = 1e-12 equivalence reference). Its mobility chain rule
     assumes `dμ/dξ = 0` and takes the ξ-path through κ alone; every
     `dvg_dξ_*_j` expression would need a `− mob·dμ_dξ[j]` term.
   - That equivalence gate **cannot catch the omission** while every Darcy test
     fixture carries a single gas component, because a one-species mixture rule
     is the identity. A multi-gas Darcy fixture has to land first.

   `test/discretization/test_darcy.jl` asserts ξ-invariance of
   `mixture_gas_viscosity` as a tripwire, so the day this changes the failure is
   loud and points here.

2. **Unvalidated Chapman–Enskog LJ defaults.** The heuristic `σ(M)` and
   `ε/k_B(M)` correlations (5.27) are rough estimates anchored to small organics
   and have not been validated against literature LJ parameters for specific
   pyrolysis gas species. Supply explicit `σ`, `ε_k` where known.

3. **Permeability microstructure evolution.** Real materials may exhibit
   permeability changes with thermal stress and cracking beyond what the
   `κ_j(T)` property functions and porosity-based closures can represent;
   the Darcy Jacobian differentiates whatever closure is supplied
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

The matching subsystem discussions are: the energy assembly that consumes
`(ρc_p)_eff` and the face conductivities — §7 (Heat Conduction, Convection &
Energy Assembly); the gas fluxes that consume `λ_eff`, `λ_j`, `κ_eff`, `μ`, and
`φ` — §9 (Gas-Phase Mass Transport & Pressure); the radiation models that consume
`α_eff` and `ε_eff` — §8 (Thermal Radiation); the swelling/volume-change use of
`ρ`, `γ_j` — §10 (Volume Change and Lateral Shrinkage); and the structured
assembly of every mixing derivative — §15 (Temporal Integration and Structured
Jacobian).
