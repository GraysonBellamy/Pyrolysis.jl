# 6. Reaction Kinetics and Stoichiometry

This chapter develops the complete condensed-phase reaction model of Pyrolysis.jl:
the Arrhenius rate law on a mass-concentration basis, the smooth temperature gates
and depletion limiter that replace the hard step functions of classical pyrolysis
codes, the stoichiometric assembly of species sources, and the heat-of-reaction
energy source with its (deliberately storage-oriented) sign convention. Every
governing equation, default value, and sign is derived from first principles and
then matched line-for-line to the implementation in `src/Physics/kinetics.jl` and
`src/Materials/reactions.jl`. We close with a careful contrast against ThermaKin2Ds,
Gpyro, and FDS, and a statement of the model's present limitations.

Throughout this chapter the index `i` denotes the **reaction index** and `j` the
**component (species) index** — the kinetics convention of overload note I1 in
§2, Nomenclature. The cell index, where it appears in the integration discussion,
is written explicitly. Mass concentration `ξ_j` (kg·m⁻³) is the primary species
state variable; stoichiometric coefficients are `ν_{i,j}` (overload note A2); the
Arrhenius pre-exponential is `A_i` (overload note A3); and the heat of reaction
`h` uses the storage sign convention `h > 0` endothermic (overload note H1).

---

## 6.1 The Arrhenius Reaction Rate Law

### 6.1.1 Continuous form

A single condensed-phase decomposition reaction $i$ proceeds at a volumetric rate
$r_i$ (kg·m⁻³·s⁻¹) given by an Arrhenius rate constant multiplied by a power-law
dependence on the mass concentrations of its reactants and gated to a temperature
window:

$$
r_i(T,\boldsymbol{\xi}) \;=\; A_i \,
\exp\!\left(-\frac{E_i}{R_g\,T}\right)
\;\prod_{j\,\in\,\mathrm{reactants}(i)} \xi_j^{\,n_{i,j}}
\,\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr)
\;\cdot\; H(T - T_{\min,i})\,H(T_{\max,i} - T).
\tag{6.1}
$$

Here

- $A_i$ is the pre-exponential factor (s⁻¹ for a first-order reaction; more
  generally with units that make $r_i$ come out in kg·m⁻³·s⁻¹ — see §6.4),
- $E_i$ is the activation energy (J·mol⁻¹),
- $R_g = 8.31446261815324$ J·mol⁻¹·K⁻¹ is the universal gas constant
  (`R_GAS`, src/Materials/constants.jl:8),
- $n_{i,j}$ is the reaction order of reaction $i$ with respect to reactant $j$
  ($n_{i,j} \ge 0$, enforced at construction),
- $\tanh(\xi_j/\xi_{\mathrm{th}})$ is the smooth **depletion rate factor**
  (§6.3), and
- $H(\cdot)$ is the Heaviside step, implemented as a **smooth tanh ramp** (§6.2).

The factorisation of $r_i$ into a temperature-only Arrhenius constant
$k_i(T) \equiv A_i \exp(-E_i/R_g T)$ and a composition-only product
$\prod_j \xi_j^{n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ is exact and is exploited
directly by the code:

```julia
k = rxn.A * exp(-rxn.E / (R_GAS * T))     # kinetics.jl:60
```

The rate constant carries the classical Arrhenius physical interpretation: $A_i$
is the collision/attempt frequency (modulated by an entropic steric factor), and
$\exp(-E_i/R_g T)$ is the Boltzmann fraction of molecular configurations that
surmount the activation barrier $E_i$. Both $A_i$ and $E_i$ are stored as plain
`Float64` scalars on the `Reaction` struct (src/Materials/reactions.jl:61–62), so
they are temperature-independent constants — see §6.10 for the consequences.

**Mass-concentration (not molar) basis.** The composition dependence in (6.1) is
written in terms of the **mass concentration** $\xi_j$ (kg of component $j$ per m³
of mixture), not molar concentration. This is the decisive modelling choice that
distinguishes a pyrolysis solver from a gas-phase combustion code. In a
decomposing solid the molar mass of the "reactant" (e.g. virgin wood, a polymer)
is ill-defined — the condensed phase is a continuum of macromolecules, not a
well-mixed ideal gas of identifiable molecules. Working in $\xi_j$ sidesteps the
need for a molar mass of the solid, makes $r_i$ directly the mass conversion rate,
and lets the stoichiometric coefficients $\nu_{i,j}$ be pure **mass** yields that
sum to unity (§6.5). The price is that $A_i$ inherits whatever units make the
product dimensionally consistent (§6.4); this is the universal convention of the
condensed-phase pyrolysis literature (ThermaKin, Gpyro, FDS all do the same — see
§6.9).

### 6.1.2 Discrete (per-cell) form

The finite-volume discretisation (§13, Finite-Volume Discretization) collocates
all state at cell centres. Within `compute_rhs!` the reaction model is evaluated
once per cell from the cell-centred state $(T, \xi_1,\dots,\xi_{N_C})$
(src/Discretization/finite_volume.jl:279–296):

```julia
T = state.T[i]
ξ = ntuple(j -> state.ξ[j, i], Val(NC))
evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)
```

so the continuous law (6.1) becomes simply $r_i^{\,\mathrm{cell}} =
r_i(T_{\mathrm{cell}}, \boldsymbol{\xi}_{\mathrm{cell}})$ — there is **no spatial
operator inside the reaction term**. Reactions are local (point-wise in space);
all spatial coupling enters through transport (conduction §7, gas diffusion/Darcy
§9) and through the ALE mesh motion (§11). This is the cell-as-batch-reactor view:
each control volume behaves, kinetically, like an isolated TGA sample, identical
to the Gpyro and FDS interpretation (§6.9).

The single-reaction kernel is `reaction_rate` (src/Physics/kinetics.jl:75–111). It
is fully scalar-type-generic (`<:Real` arguments), so `ForwardDiff.Dual` numbers
propagate cleanly through it for Jacobian and sensitivity work (§15, §17).

> **Implementation note.** Single-reaction rate: `reaction_rate`,
> src/Physics/kinetics.jl:75. All-reaction evaluation: `evaluate_reactions!`
> dispatching to the `@generated` kernel `_evaluate_reactions_gen!`,
> src/Physics/kinetics.jl:109–142.

---

## 6.2 Smooth Temperature Gates

### 6.2.1 Why hard steps are unacceptable

The temperature window $[T_{\min,i}, T_{\max,i}]$ serves two purposes (following
ThermaKin's original motivation): it spares the solver from evaluating a reaction
where it is kinetically negligible, and it lets a reaction be used to represent a
thermodynamically-controlled phase transition that is sharply switched on at a
threshold. The textbook expression of this is the product of two Heaviside steps,
$H(T - T_{\min,i})\,H(T_{\max,i} - T)$, equal to 1 inside the window and 0 outside.

A hard step is, however, **non-differentiable** at the threshold and has a zero
derivative everywhere else. Both properties are poison for the implicit machinery
of Pyrolysis.jl:

1. The stiff time integrators solve nonlinear systems with Newton's method, whose
   convergence depends on a meaningful Jacobian $\partial r_i/\partial T$. At a
   hard step the Jacobian is a Dirac delta (numerically: undefined), and away from
   it the step contributes nothing — Newton sees a flat rate that suddenly jumps.
2. Automatic differentiation (`ForwardDiff`) of a discontinuous function returns a
   derivative of zero almost everywhere and is simply wrong at the breakpoint,
   corrupting both the structured Jacobian (§15) and forward-mode sensitivities
   (§17).

### 6.2.2 The tanh ramp

Pyrolysis.jl replaces each Heaviside by a hyperbolic-tangent ramp. The lower gate
is applied only when $T_{\min,i}$ is finite and strictly positive; the upper gate
only when $T_{\max,i}$ is finite:

$$
g_{\min,i}(T) =
\begin{cases}
1, & T_{\min,i}\ \text{not finite or}\ \le 0,\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh(T - T_{\min,i})\bigr), & \text{otherwise},
\end{cases}
\qquad
g_{\max,i}(T) =
\begin{cases}
1, & T_{\max,i}\ \text{not finite},\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh(T_{\max,i} - T)\bigr), & \text{otherwise},
\end{cases}
\tag{6.2}
$$

and the gated rate is $r_i \leftarrow r_i\, g_{\min,i}(T)\, g_{\max,i}(T)$. The
implementation is verbatim:

```julia
if isfinite(rxn.T_min) && rxn.T_min > 0.0
    r *= 0.5 * (1 + tanh(T - rxn.T_min))
end
if isfinite(rxn.T_max)
    r *= 0.5 * (1 + tanh(rxn.T_max - T))
end
```

(src/Physics/kinetics.jl:102–108). The default gates are $T_{\min}=0$ and
$T_{\max}=\infty$ (src/Materials/reactions.jl:128–129, 144–145), so an
unconstrained reaction incurs no gating cost at all — the `T_min > 0.0` guard
short-circuits the lower ramp and `isfinite(Inf) == false` short-circuits the
upper ramp.

The ramp $\tfrac12(1 + \tanh(\Delta T))$ rises from $\approx 0$ to $\approx 1$
over a width set by the argument scale of `tanh`. Because the argument is
$(T - T_{\min})$ **in Kelvin with unit coefficient**, the transition has a
characteristic half-width of $\mathcal O(1\text{–}3)$ K and reaches
$\tfrac12(1+\tanh 5)\approx 0.99995$ only about 5 K inside the window. Equivalently,
the gate is $\approx 4.5\times10^{-5}$ about 5 K **outside** the window.

### 6.2.3 The ~10 K leakage caveat

The smoothing is not free: the reaction is "weakly active" within roughly 10 K of
each boundary, where a hard step would have it identically off. At the lower edge
$\tfrac12(1+\tanh(T - T_{\min}))$ takes the value $0.5$ exactly at $T = T_{\min}$,
$\approx 0.018$ at $T_{\min} - 2$ K, and $\approx 4.5\times10^{-5}$ at
$T_{\min}-5$ K. For decomposition kinetics, where $T_{\min}$ is a numerical
convenience rather than a sharp physical threshold, this leakage is immaterial
because the Arrhenius factor is itself vanishingly small that far below the
reaction's effective onset. For a reaction used as a sharp phase-transition
switch the user should be aware that the transition is smeared over a few Kelvin;
if a crisper gate is needed, $T_{\min}$ should be offset to recentre the ramp.
There is currently no way to tune the ramp width (it is hard-coded to a unit
Kelvin scale); see §6.10.

> **Implementation note.** Temperature gates: `reaction_rate`,
> src/Physics/kinetics.jl:102–108. Default `T_min = 0.0`, `T_max = Inf` set in the
> convenience constructors, src/Materials/reactions.jl.

---

## 6.3 The Depletion Limiter

### 6.3.1 The stiffness problem near zero concentration

As a reactant is consumed, $\xi_j \to 0$. For a first-order reaction ($n_{i,j}=1$)
the bare rate $r_i = k_i(T)\,\xi_j$ approaches zero linearly, which is benign. But
the **Jacobian entry** $\partial r_i/\partial \xi_j = k_i(T)$ does not vanish: at
high temperature $k_i(T)$ can be enormous ($A_i$ ranges over $10^8$–$10^{20}$ s⁻¹,
src/Materials/constants.jl:28). A control volume holding a numerically tiny,
possibly slightly-negative-from-roundoff residue of reactant therefore presents
the implicit solver with a stiff, fast-relaxing mode whose eigenvalue is
$-k_i(T)$ even though essentially no mass remains to react. This manifests as
time-step throttling, Newton non-convergence, and — if $\xi_j$ overshoots below
zero and $n_{i,j}$ is non-integer — complex or `NaN` rates from $\xi_j^{n_{i,j}}$.

Sub-first-order fits make this strictly worse: for $0 < n_{i,j} < 1$ the bare
$\partial r_i/\partial\xi_j \sim \xi_j^{\,n-1}$ **diverges** at depletion, and for
a zero-order fit ($n_{i,j}=0$, common in TGA analysis) the bare rate does not
vanish at all — left unguarded it keeps consuming reactant at $k_i(T)$ straight
through $\xi_j = 0$ and integrates the concentration negative without bound.

### 6.3.2 The tanh rate factor

Pyrolysis.jl rolls the rate off smoothly as a reactant depletes. Each reactant
contributes a multiplicative **depletion factor to the rate** (not to the
power-law base):

$$
r_i \;=\; k_i(T)\prod_{j\,\in\,\mathrm{reactants}(i)}
\xi_j^{\,n_{i,j}}\;\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr),
\tag{6.3}
$$

with a hard floor at zero applied first: $\xi_j \leftarrow \max(\xi_j, 0)$. The
implementation:

```julia
ξ_raw = max(ξ[idx], zero(eltype(ξ)))
conc_factor *= _fast_pow(ξ_raw, order)
if limiter.enabled
    conc_factor *= _depletion_factor(ξ_raw, limiter.threshold)
end
```

(src/Physics/kinetics.jl:89–95). The `max(ξ, 0)` clip guarantees a non-negative
base for the power even if the integrator transiently produces $\xi_j < 0$, and
it pins the depletion factor at $\tanh 0 = 0$ there — an overshooting cell
consumes nothing. `_depletion_factor` short-circuits to `one(ξ)` above
$20\,\xi_{\mathrm{th}}$, which is **bit-identical** to evaluating the `tanh`
(Float64 `tanh(x)` rounds to exactly `1.0` for $x \gtrsim 19.1$), so the common
path — bulk concentrations orders of magnitude above the threshold — pays no
transcendental (src/Physics/kinetics.jl:68–73).

Applying the factor to the rate rather than substituting an effective
concentration into the power-law base matters for low orders: a base-applied
limiter leaves a zero-order reactant unguarded ($\xi_{\mathrm{eff}}^{\,0}\equiv 1$)
and, for $0 < n < 0.5$, behaves as $(\xi^2/\xi_{\mathrm{th}})^n$ near zero with a
still-singular $\partial r/\partial\xi \sim \xi^{\,2n-1}$. The multiplicative
form removes the singular eigenvalue for **every** order $n \ge 0$.

### 6.3.3 Behaviour and default values

The limiter is configured by the `DepletionLimiter` struct
(defined in src/Materials/reactions.jl; re-exported from `Pyrolysis.Physics`),
whose defaults are

| field | default | meaning |
|---|---|---|
| `threshold` ($\xi_{\mathrm{th}}$) | $1.0$ kg·m⁻³ | rate roll-off scale |
| `enabled` | `true` | set `false` to recover the bare Arrhenius rate |

Each `Material` carries its own limiter in `material.depletion_limiter`, set
with the `Material` constructor's `depletion_limiter` keyword (default
`DepletionLimiter()`); the solver's kinetics path `evaluate_reactions!` reads it
for every rate evaluation. `DEFAULT_DEPLETION_LIMITER = DepletionLimiter()`
(src/Physics/kinetics.jl) backs only the standalone `reaction_rate` entry point
when no `limiter` keyword is passed.

The factor $\tanh(\xi/\xi_{\mathrm{th}})$ has the key properties:

- **Smooth everywhere.** There is no threshold branch: the rate is $C^\infty$ in
  $\xi > 0$ (the only non-smooth point left is the `max(ξ, 0)` clip exactly at
  zero). Neither the rate nor $\partial r/\partial\xi$ jumps anywhere.
- $\xi \gg \xi_{\mathrm{th}}$: the bare rate is recovered double-exponentially —
  the relative deficit is $2e^{-2\xi/\xi_{\mathrm{th}}}$, i.e. $\sim\!3.7\%$ at
  $2\,\xi_{\mathrm{th}}$, $4\times10^{-9}$ at $10\,\xi_{\mathrm{th}}$, and exactly
  zero in Float64 beyond $\approx 19.1\,\xi_{\mathrm{th}}$.
- $\xi \to 0$: the rate vanishes as $k\,\xi^{\,n+1}/\xi_{\mathrm{th}}$ with
  $\partial r/\partial\xi \to 0$ for $n > 0$ and $\to k/\xi_{\mathrm{th}}$
  (finite) for $n = 0$ — the stiff eigenvalue is removed for every order
  $n \ge 0$. For $n = 1$ the limited rate is exactly
  $k\,\xi\tanh(\xi/\xi_{\mathrm{th}})$.

The robustness is bought with a rate **bias near depletion**: the suppression is
$\sim\!24\%$ at $\xi = \xi_{\mathrm{th}}$ ($\tanh 1 \approx 0.762$) and grows as
the reactant runs out — $90\%$ at $0.1\,\xi_{\mathrm{th}}$. For first-order
kinetics the limited rate behaves as $k\,\xi^2/\xi_{\mathrm{th}}$ near zero,
which converts the exponential burnout tail into an algebraic
$\xi \approx \xi_{\mathrm{th}}/(kt)$ tail. Choose $\xi_{\mathrm{th}}$ well below
the concentrations that matter for the quantities of interest.

Reaction orders are validated at construction: $n < 0$ is rejected
(`_validate_reaction_order`, src/Materials/reactions.jl), since a negative-order
rate is singular at depletion and no limiter can repair it.

For validation studies that need the unbiased Arrhenius rate, the limiter is
disabled per call:

```julia
r = reaction_rate(rxn, ξ, T; limiter = DepletionLimiter(enabled = false))
```

or for an entire solve, per material:

```julia
Material(name = …, components = …, reactions = …,
         depletion_limiter = DepletionLimiter(enabled = false))
```

Disabling removes the only depletion guard: zero-order reactions then consume
reactant through $\xi = 0$ without bound, and orders $0 < n < 1$ have a singular
$\partial r/\partial\xi$ at depletion.

> **Implementation note.** Limiter struct: `DepletionLimiter`,
> src/Materials/reactions.jl (carried per material in
> `material.depletion_limiter`; re-exported from Pyrolysis.Physics).
> Per-reactant application inside the concentration product loop:
> src/Physics/kinetics.jl.

---

## 6.4 Reaction Order and the Concentration Product

The composition factor is the product over all reactants of the clipped
concentration raised to its per-reactant order, times the per-reactant
depletion factor of §6.3:

$$
\Phi_i(\boldsymbol{\xi}) = \prod_{j\,\in\,\mathrm{reactants}(i)}
\xi_j^{\,n_{i,j}}\,\tanh\!\bigl(\xi_j/\xi_{\mathrm{th}}\bigr).
\tag{6.4}
$$

Each reactant is stored as an `(component_index, reaction_order)` pair in the
`reactants::NTuple{NRin, Tuple{Int,Float64}}` field
(src/Materials/reactions.jl:59), so the order $n_{i,j}$ is a **`Float64`** and may
be non-integer (e.g. $n=0.5$, $n=1.5$, common in TGA fits and in nucleation/
diffusion reaction models). Orders must satisfy $n_{i,j} \ge 0$ — negative
orders are rejected at construction (§6.3.3). The tuple is iterated with the
loop fully unrolled by the compiler because `NRin` is a static type parameter.

The exponentiation routes through a fast-path helper that special-cases the two
most common orders:

```julia
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)
```

(src/Physics/kinetics.jl:66). For $n=1$ it returns $x$ directly (no `pow` call);
for $n=0$ it returns `one(x)` — a zero-order reactant contributes no
concentration dependence to the power law, but its depletion factor
$\tanh(\xi/\xi_{\mathrm{th}})$ still gates the rate to zero as it runs out
(§6.3.2); otherwise it falls back to the general
`x^n`. Because the `max(ξ,0)` clip (§6.3) guarantees a non-negative base,
$x^n$ is always real-valued even for fractional $n$. The accumulator is
initialised at the promoted scalar type so Duals carry through:

```julia
conc_factor = one(promote_type(eltype(ξ), typeof(T)))
```

(src/Physics/kinetics.jl:89).

**Units of the pre-exponential.** Dimensional consistency of (6.1) forces

$$
[A_i] = \text{kg}\,\text{m}^{-3}\,\text{s}^{-1}
\Big/ \Big(\text{kg}\,\text{m}^{-3}\Big)^{\sum_j n_{i,j}}
= \bigl(\text{kg}\,\text{m}^{-3}\bigr)^{1-\sum_j n_{i,j}}\,\text{s}^{-1}.
$$

For a first-order single reactant ($\sum_j n_{i,j}=1$) this reduces to s⁻¹; for a
second-order/bimolecular case it becomes m³·kg⁻¹·s⁻¹. The struct docstring records
this convention ("[s⁻¹] or [m³/(kg·s)]", src/Materials/reactions.jl:37) but the
code does **not** check it — the user is responsible for supplying $A_i$ in units
consistent with the chosen orders.

---

## 6.5 Stoichiometry and Species Source Terms

### 6.5.1 The species source

Summing the per-reaction mass conversion over all reactions gives the volumetric
chemical source for component $j$ (kg·m⁻³·s⁻¹):

$$
S_j(T,\boldsymbol{\xi}) \;=\;
\left.\frac{\partial \xi_j}{\partial t}\right|_{\mathrm{rxn}}
\;=\; \sum_{i=1}^{N_R} \nu_{i,j}\, r_i(T,\boldsymbol{\xi}),
\tag{6.5}
$$

where $\nu_{i,j}$ is the stoichiometric (mass) coefficient: **negative** when $j$
is a reactant of reaction $i$, **positive** when $j$ is a product. This $S_j$ is
the chemical contribution to the species conservation equation of §3, Governing
Equations; for gases it adds to the transport divergence, for solids/liquids it is
the entire right-hand side (§6.5.4).

### 6.5.2 The single-reactant mass convention

The current model supports exactly one reactant per reaction. Per unit reaction
extent, that reactant is consumed with **unit mass stoichiometry** (an implicit
$\nu = -1$), and the products are produced with their explicit per-product mass
yields. The kernel that applies one reaction to the source vector is:

```julia
@inline function _apply_reaction_source!(src, rxn::Reaction{NRin}, r) where {NRin}
    NRin == 1 || error("Multi-reactant reaction $(rxn.name) (NRin=$NRin) ...")
    @inbounds src[rxn.reactants[1][1]] -= r          # reactant: ν = -1
    @inbounds for (idx, ν) in rxn.products
        src[idx] += ν * r                            # products: ν = yield ≥ 0
    end
    return nothing
end
```

(src/Physics/kinetics.jl:182–191). The `NRin == 1` guard **errors at runtime** if
a multi-reactant reaction is ever applied — multi-reactant support would require
extending the reactant tuple to a 3-ary `(idx, order, mass_stoich)` form so that a
mass stoichiometry (not just a rate order) could be attached to each reactant
(src/Materials/reactions.jl:43–49). Note the asymmetry baked into the model: a
reactant carries a rate **order** $n_{i,j}$ (how its concentration scales the
rate) but its **mass** stoichiometry is fixed at $1$; a product carries a mass
**yield** but no rate order (products never appear in $\Phi_i$).

For a single-reactant reaction, then, the coefficients are

$$
\nu_{i,\text{reactant}} = -1, \qquad
\nu_{i,\text{product }p} = y_{i,p} \ge 0,
$$

and (6.5) becomes, for reaction $i$ acting on its reactant $a$ and products $p$:

$$
S_a \mathrel{-}= r_i, \qquad S_p \mathrel{+}= y_{i,p}\, r_i.
$$

### 6.5.3 Mass conservation as a build-time invariant

Total mass is conserved by reaction $i$ if and only if the consumed reactant mass
equals the produced product mass, i.e.

$$
\sum_{p\,\in\,\mathrm{products}(i)} y_{i,p} = 1.
\tag{6.6}
$$

This is enforced **at construction time**, not at run time — in the `Reaction`
inner constructor itself, so no construction path bypasses it.
`_validate_stoichiometry` (src/Materials/reactions.jl) sums the product yields
and requires $\sum_p y_{i,p} = 1$ against the implicit unit reactant mass of
the single-reactant convention. Each yield is additionally required to be
non-negative (`_validate_yield`, enforced even with `validate_mass = false` —
a pair like $(1.5, -0.5)$ sums to 1 but would silently drive the
negative-yield species' concentration negative without bound), and each
reaction order must satisfy $n \ge 0$ (`_validate_reaction_order`). The
tolerance is relative to the total coefficient mass with an absolute floor:

$$
\text{tol} = \max\!\bigl(10^{-10}\,(m_{\mathrm{react}} + m_{\mathrm{prod}}),\ 10^{-14}\bigr),
$$

(src/Materials/reactions.jl). A violation raises an error at construction unless
the user passes `validate_mass = false`, which disables only the sum check.

Worked example (the canonical charring reaction):

```julia
Reaction(:charring, 1 => (2, 3); A=1e15, E=200e3, h=500e3, yields=(0.2, 0.8))
```

assigns $\nu_{i,1} = -1$ (virgin consumed), $\nu_{i,2} = 0.2$ (char produced),
$\nu_{i,3} = 0.8$ (gas produced); the check $0.2 + 0.8 = 1.0$ passes. The species
sources at rate $r_i$ are then $S_1 = -r_i$, $S_2 = +0.2\,r_i$, $S_3 = +0.8\,r_i$,
whose sum is exactly $0$ — mass is conserved instantaneously.

Because conservation is guaranteed by the build-time check and the implicit
$\nu=-1$ convention, **no run-time projection or Lagrange multiplier is used**; the
species sum $\sum_j S_j$ is identically zero up to floating-point roundoff. The
validation utilities of §6.8 exist purely as diagnostics.

### 6.5.4 Gas-vs-solid dispatch

The species source $S_j$ is the same expression for every phase; the phase
distinction lives in the *conservation equation* into which $S_j$ is inserted, not
in the kinetics. At final RHS assembly (src/Discretization/finite_volume.jl:403)
the code calls `apply_species_residual_unrolled!`, which assembles

$$
\frac{d\xi_j}{dt} =
\begin{cases}
-(\nabla\!\cdot\! \mathbf J_j)_{\mathrm{cell}} + S_j, & j\ \text{gaseous},\\[4pt]
S_j, & j\ \text{solid or liquid},
\end{cases}
\tag{6.7}
$$

so solids/liquids have **no transport term** — they are Lagrangian, fixed to the
material point (in the ALE frame, §11) — while gases carry the divergence of their
mass flux $\mathbf J_j$ (§9). The dispatch is resolved at compile time from the
component's static phase type (`is_gas`, src/Materials/components.jl), with no
run-time branching. The full energy/species RHS assembly sequence is documented in
§7, Heat Conduction, Convection & Energy Assembly.

> **Implementation note.** Species source: `species_source_terms!` and the
> `@generated` kernel `_species_source_terms_gen!`, src/Physics/kinetics.jl:159–242;
> per-reaction application `_apply_reaction_source!`, src/Physics/kinetics.jl:182.
> Build-time balance: `_validate_stoichiometry`, src/Materials/reactions.jl:90.

---

## 6.6 Heat of Reaction

### 6.6.1 The reaction heat source and its sign convention

The volumetric heat liberated or absorbed by chemistry is

$$
Q_{\mathrm{rxn}}(T,\boldsymbol{\xi}) \;=\; -\sum_{i=1}^{N_R} h_i(T,\boldsymbol{\xi})\, r_i(T,\boldsymbol{\xi}),
\tag{6.8}
$$

(W·m⁻³), implemented as

```julia
@inline function reaction_heat_source(rates, heats)
    Tu = promote_type(eltype(rates), eltype(heats))
    Q = zero(Tu)
    @inbounds for i in eachindex(rates, heats)
        Q -= heats[i] * rates[i]          # Note the negative sign!
    end
    return Q
end
```

(src/Physics/kinetics.jl:257–264). The single minus sign is the entire content of
the sign convention, and it is worth stating precisely because it is the most
error-prone point in the whole model.

The field `rxn.heat` stores $h_i$ as an **enthalpy-of-reaction**: the heat
*absorbed* per kilogram of (first) reactant consumed. The storage convention is
(overload note H1):

$$
\boxed{\;h_i > 0 \Rightarrow \text{endothermic (absorbs heat, cools)};\quad
h_i < 0 \Rightarrow \text{exothermic (releases heat, heats)}.\;}
$$

The minus sign in (6.8) converts this enthalpy convention into the **heat-source**
convention required by the energy equation, where a positive source heats the
material:

- **Endothermic** ($h_i > 0$, $r_i > 0$): $Q_{\mathrm{rxn}} = -h_i r_i < 0$ — the
  cell is cooled, as a decomposition that consumes heat should.
- **Exothermic** ($h_i < 0$, $r_i > 0$): $Q_{\mathrm{rxn}} = -h_i r_i > 0$ — the
  cell is heated, as char oxidation should.

The source comment states this explicitly (src/Physics/kinetics.jl:251–255). In
the energy equation $Q_{\mathrm{rxn}}$ enters additively on the right-hand side
together with the conduction divergence, gas-advection source, and radiation
source, and the whole is divided by the matrix-only volumetric heat capacity
$\rho c_p^{\mathrm{eff}}$ (gas storage excluded — Formulation A; see §3, §5, §7):

$$
\frac{dT_{\mathrm{cell}}}{dt} =
\frac{-\,(\nabla\!\cdot\! q_{\mathrm{cond}}) + S_{\mathrm{conv}} + S_{\mathrm{gen}}
+ Q_{\mathrm{rad}} + Q_{\mathrm{rxn}}}{\rho c_p^{\mathrm{eff}}}
\tag{6.9}
$$

(src/Discretization/finite_volume.jl:393–400).

**Units check.** $[h_i] = $ J·kg⁻¹ (per kg of first reactant) and $[r_i] = $
kg·m⁻³·s⁻¹, so $[h_i r_i] = $ J·m⁻³·s⁻¹ $=$ W·m⁻³, the volumetric power density
demanded by (6.9). The "per kg of first reactant" basis is consistent with the
single-reactant mass convention of §6.5: the reactant is consumed at mass rate
$r_i$, and $h_i$ is the enthalpy per unit of that mass.

### 6.6.2 Temperature- and state-dependent heat

$h_i$ is not restricted to a constant. The `heat::Ph` field is typed by any
`AbstractPropertyFunction` (src/Materials/reactions.jl:63), and in
`evaluate_reactions!` the heat is evaluated with **both** $T$ and the full
concentration tuple $\boldsymbol{\xi}$:

```julia
heats[$i] = rxn.heat(T, ξ)   # Pass ξ for state-dependent heats
```

(src/Physics/kinetics.jl:136). The two-argument call dispatches as follows
(src/Materials/property_functions.jl:188–197):

- `ConstantProperty(h)` $\to h$ (the `(T,ξ)` method ignores both arguments);
- `LinearProperty(a,b)` $\to a + bT$ (affine in temperature, $\xi$ ignored);
- `PolynomialProperty` $\to \mathrm{evalpoly}(T,\text{coeffs})$;
- `StateDependentProperty(f)` $\to f(T, \boldsymbol{\xi})$ — the only type that
  actually consumes $\boldsymbol{\xi}$.

A constant heat is created by passing a number, an affine heat by passing a
`(a,b)` tuple, through `to_property` (src/Materials/property_functions.jl:387–390):

```julia
Reaction(:decomp, 1 => 2; A=8.5e12, E=188e3, h=870e3)          # ConstantProperty
Reaction(:decomp, 1 => 2; A=8.5e12, E=188e3, h=(870e3, 100.0)) # LinearProperty
```

(Positive `h` = endothermic, the storage convention of §6.6.1; PMMA-style
depolymerization absorbs ≈870 kJ·kg⁻¹.)

The principal use of `StateDependentProperty` is the **heat of moisture sorption**
in hygroscopic materials, where the differential heat depends on the current
moisture content $MC$. A representative closure (the form used in the moisture
examples) is

$$
h_{\mathrm{sorp}}(T,\boldsymbol{\xi}) =
\Bigl(2275 + 860\,e^{-0.120\,MC}\Bigr)\times 10^{3}\ \text{J·kg}^{-1},
\qquad MC = 100\,\frac{\xi_{\mathrm{moist}}}{\sum_{k}\xi_k},
$$

encoded as a closure over the moisture index (src/Materials/property_functions.jl:
158–162). Because the closure is passed $\boldsymbol{\xi}$ at evaluation, the
sorption heat self-consistently rises as the wood dries — a behaviour impossible
to represent with a temperature-only property. Calling a `StateDependentProperty`
with only $T$ raises an error (src/Materials/property_functions.jl:199–203),
guaranteeing the state argument is never silently dropped.

The midpoint-rule sensible-enthalpy helper $\Delta h(T_{\mathrm{hi}},T_{\mathrm{lo}})
= c_p(\tfrac{T_{\mathrm{hi}}+T_{\mathrm{lo}}}{2})(T_{\mathrm{hi}}-T_{\mathrm{lo}})$
(src/Materials/components.jl:515–518) is used for transpiration/gas-generation
enthalpy bookkeeping (§7, §12), **not** for the reaction heat itself. Whether
the reaction heat is "separate from" the sensible-enthalpy transport of the
released gases is precisely the **datum question** of §6.6.3 — the partition
between "chemical" and "sensible" enthalpy depends on the temperature at which
$h_i$ is referenced.

> **Implementation note.** Heat source: `reaction_heat_source`,
> src/Physics/kinetics.jl:257. Heat evaluation with state:
> `_evaluate_reactions_gen!`, src/Physics/kinetics.jl:136. Property dispatch:
> src/Materials/property_functions.jl:188–197.

### 6.6.3 The enthalpy datum: heats at reaction temperature vs the $T_{\text{ref}}$ datum

A heat of reaction is only defined relative to a **datum**: $h_i$ may be quoted

- **at the local reaction temperature** — $h_i(T) = \sum_j \nu_{ji}\,h_j(T)$,
  the heat a DSC measures at the decomposition temperature, and the usual
  meaning of inverse-calibrated (ThermaKin/Gpyro-style) effective heats; or
- **at the fixed ambient datum** $T_{\text{ref}} = 298.15$ K —
  $h_i = \sum_j \nu_{ji}\,h_j(T_{\text{ref}})$, the formation-enthalpy style.

The distinction is immaterial for the default `:standard` energy form: with
heats at the local reaction temperature, subtracting the species equations
from the total-enthalpy balance gives the **exact** temperature equation
$\rho c_p\,\partial_t T = -\nabla\cdot q + S_{\text{conv}} + Q_{\text{rad}} +
Q_{\text{rxn}}$ — no generation-sink term exists, because the gas products'
sensible enthalpy $\Delta h_g(T, T_{\text{ref}})$ is already *inside*
$h_i(T)$. The default form is therefore the consistent one for DSC-style
heats (its only approximation is the dropped gas heat storage, ≈0.2%).

The distinction is **load-bearing** for `energy_form = :with_generation_sink`
(§7.5): the sink $S_{\text{gen}} = -\sum_g \Delta h_g(T,T_{\text{ref}})\,
\dot m_g'''$ is exactly the conversion term between the two datums, so it is
consistent **only with $T_{\text{ref}}$-datum heats**. Enabling it with
DSC-style heats subtracts the gas sensible enthalpy a *second* time —
$\Delta h_g(900\,\text{K}, 298.15\,\text{K}) = 0.66$–$1.20$ MJ·kg⁻¹ for
$c_{p,g} = 1.1$–$2.0$ kJ·kg⁻¹·K⁻¹, comparable to or larger than typical
wood/polymer pyrolysis endotherms (0.1–0.9 MJ·kg⁻¹). `to_problem_def` emits a
one-time warning when the sink is selected (src/Problem/adapter.jl).

---

## 6.7 Allocation-Free Evaluation

The reaction model sits in the innermost loop of every residual and Jacobian
evaluation, so it is engineered to allocate nothing and to specialise completely
on the (compile-time-known) reaction and component counts $N_R$ and $N_C$. Three
`@generated` kernels achieve this.

**Rate and heat evaluation.** `evaluate_reactions!` forwards to
`_evaluate_reactions_gen!` (src/Physics/kinetics.jl:122–142), which the compiler
expands into a flat block — one explicit `rates[i] = reaction_rate(...)` and
`heats[i] = rxn.heat(T, ξ)` statement per reaction — with no run-time loop over
the heterogeneous `reactions` tuple. Heterogeneity (each `Reaction` may have a
different `NRin`, `NRout`, and heat property type) is therefore handled by static
dispatch with no `Union` boxing and no allocation.

**Species sources.** `species_source_terms!` zeroes the output and calls
`_species_source_terms_gen!` (src/Physics/kinetics.jl:201–242), again unrolled to
one `_apply_reaction_source!` per reaction. A second method accepts a plain `Int`
for `NC` and routes to the same `@generated` kernel
(src/Physics/kinetics.jl:213–222), which is the form called from `compute_rhs!`.
There is also a tuple-returning variant (src/Physics/kinetics.jl:159–175) that
builds the source in a stack-allocated `MVector` for local use.

**Heat sum.** `reaction_heat_source` is a single `@inline` reduction over the
pre-filled `rates`/`heats` vectors (src/Physics/kinetics.jl:257–264), promoting the
accumulator type so Duals carry through.

All three use the promoted scalar type so that, in a Dual-typed (AD) pass, the
caches hold Duals and derivatives propagate end-to-end. The per-cell scratch
(`rates`, `heats`, `dξ_rxn`, `Q_rxn`) lives in the RHS cache and is sized by
$N_R$, $N_C$, and the cell count (src/Physics/source_terms.jl:59–95;
src/Discretization/finite_volume.jl:123–154), so the hot loop reuses storage
across cells and time steps.

> **Implementation note.** Generated kernels: `_evaluate_reactions_gen!`
> (src/Physics/kinetics.jl:146), `_species_source_terms_gen!`
> (src/Physics/kinetics.jl:249). Compile-time constants $N_C, N_R$ are the
> `Material{NC,NR}` type parameters (src/Materials/materials.jl).

---

## 6.8 Validation Utilities

The kinetics module ships several diagnostics. They are not in the production hot
path; they exist to catch stoichiometry errors during model development.

**`verify_stoichiometry(rxn; tol=1e-10)`** (src/Materials/reactions.jl:518–525)
computes the mass-balance error of a single reaction via `mass_balance_error`
(src/Materials/reactions.jl:503–511):

$$
\varepsilon_{\mathrm{mb}}(i) = \Bigl|\,1 - \textstyle\sum_p y_{i,p}\,\Bigr|.
$$

For a single-reactant reaction it compares the implicit unit reactant mass against
the product-yield sum; for $N_{\mathrm{in}} > 1$ it returns `NaN` (unsupported).
The function warns and returns `false` if $\varepsilon_{\mathrm{mb}} > $ `tol`.

**`check_material_mass_balance(reactions; tol=1e-10)`**
(src/Materials/reactions.jl:547–557) iterates `verify_stoichiometry` over a
material's reaction tuple, returning `true` only if all balance.

**`validate_material(m)`** (src/Materials/materials.jl:437–457) is called as part
of material assembly: it checks $N_C \ge 1$, that every reactant/product index lies
in $1..N_C$, and that each reaction passes `verify_stoichiometry` (warning on
failure). Symbolic reactions are resolved to integer indices first by
`resolve_reaction` against the component name→index map
(src/Materials/reactions.jl:351–395), so a reaction that names a non-existent
component errors with a clear message before any solve.

**`validate_species_source_conservation(sources; rtol=1e-12)`**
(src/Physics/kinetics.jl:296–304) audits a computed source vector at run time:

$$
\Bigl|\textstyle\sum_j S_j\Bigr| \le \mathrm{rtol}\cdot \max_j |S_j|,
\qquad \max_j|S_j| \ \text{floored at}\ 1,
$$

with default $\mathrm{rtol} = 10^{-12}$ — a very tight relative tolerance
appropriate because the sum is mathematically zero by construction (§6.5.3) and
should fail only on a genuine bug. The checked wrapper
`species_source_terms_checked!` (src/Physics/kinetics.jl:337–360) computes the
sources, runs the audit, and either warns or (with `throw_on_error=true`) errors.
For production runs the unchecked `species_source_terms!` is used directly; the
audit's cost (one sum, one max) is negligible but pointless once a material is
known-good.

> **Implementation note.** `verify_stoichiometry`/`mass_balance_error`/
> `check_material_mass_balance`: src/Materials/reactions.jl:503–557.
> `validate_species_source_conservation`/`species_source_terms_checked!`:
> src/Physics/kinetics.jl:296–360. Material-level validation:
> `validate_material`, src/Materials/materials.jl:437.

---

## 6.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

All three reference codes solve essentially the same condensed-phase chemistry —
mass-concentration Arrhenius kinetics with mass-yield stoichiometry — and
Pyrolysis.jl deliberately reproduces their physics. The differences are in
generality, in the gating mechanism, and most consequentially in the **sign and
basis of the heat of reaction**. We map each code's notation to ours on first use.

### 6.9.1 ThermaKin2Ds

ThermaKin (Stoliarov & Lyon, *ThermaKin2Ds Reference v2*, Eq. 5) writes a reaction
as a generalised two-reactant, two-product step

$$
\theta_1\,\text{COMP1} + \theta_2\,\text{COMP2} \to
\theta_3\,\text{COMP3} + \theta_4\,\text{COMP4} + h,
$$

with stoichiometric coefficients $\theta_i^{\,j}$ (their notation) which map to our
$\nu_{i,j}$ (overload note A2 — $\theta$ is reserved for strain rate in our
nomenclature). The rate is **explicitly bimolecular** (ThermaKin Eq. 6):

$$
r = A\,\exp\!\Bigl(-\frac{E}{RT}\Bigr)\,\xi_{\mathrm{COMP1}}\,\xi_{\mathrm{COMP2}},
\qquad \xi = mf\cdot\rho,
$$

with the convention that "in the absence of a second reactant, $\xi_{\mathrm{COMP2}}$
is set to 1." Pyrolysis.jl differs in two ways. First, our rate is **n-ary in the
order but single-reactant in the mass**: we allow an arbitrary order $n_{i,j}$ on
the single supported reactant (a strict generalisation of ThermaKin's
unit-exponent factor for COMP1), but we do not yet support a genuine second
mass-bearing reactant (§6.10). ThermaKin's COMP2 mechanism is precisely what our
multi-reactant extension would need. Second, ThermaKin's temperature limit is a
**hard cutoff** ("if temperature decreases below the lower limit or increases above
the upper limit, the rate of reaction is set to 0"), whereas we use a smooth tanh
ramp (§6.2) for AD/Newton compatibility. ThermaKin's species/heat assembly
($\partial\xi_j/\partial t = \sum_i \theta_i^{\,j} r_i + \text{transport}$, their
Eqs. 15–16, and heat $= h\,r$) is identical in form to our (6.5) and (6.8).

The pivotal difference is the **heat-of-reaction sign**. ThermaKin states
explicitly: "Following a convention frequently adopted in pyrolysis studies, a
positive $h$ is used to represent an **exothermic** reaction." This is the
**opposite** of our storage convention ($h>0$ endothermic). The two are reconciled
by the minus sign in our $Q_{\mathrm{rxn}} = -\sum_i h_i r_i$: a ThermaKin input
of $h_{\mathrm{TK}} = +870$ kJ·kg⁻¹ (their exothermic) corresponds to our
$h = -870$ kJ·kg⁻¹, and both yield a positive (heating) source. A user porting a
ThermaKin material **must flip the sign of every heat of reaction.** (ThermaKin
also notes its reaction parameters are not constrained by mass or energy
conservation — the user owns that; Pyrolysis.jl instead enforces mass balance at
construction, §6.5.3.)

### 6.9.2 Gpyro

Gpyro (Lautenberger, *Gpyro Technical Reference*) tracks a **conversion** variable
$\alpha_{A_k}$ for condensed-phase reactant $A_k$ and writes the $n^{\text{th}}$
-order heterogeneous rate (their Eq. 39, with $Z_k \to A_i$, $E_k\to E_i$,
$n_k\to n_{i,j}$, and the Gpyro volumetric rate $\dot\omega''' \to r$):

$$
r_k = f(\alpha_{A_k})\,
\frac{(\bar\rho\,Y_{A_k}\Delta z)_\Sigma}{\Delta z}\,
Z_k\,e^{-E_k/RT}\,g(Y_{\mathrm{O_2}})\,H(T - T_{\mathrm{crit},k}),
\qquad f(\alpha_{A_k}) = (1-\alpha_{A_k})^{n_k}.
$$

Two Gpyro features have no analogue in Pyrolysis.jl: the oxygen-sensitivity factor
$g(Y_{\mathrm{O_2}})$ (their Eq. 40 — our kinetics have no $Y_{\mathrm{O_2}}$
dependence, §6.10), and the rich library of specialised reaction models
$f(\alpha)$ (nucleation, diffusion, Avrami, etc., their Table 1 — we offer only the
$n^{\text{th}}$-order power law of (6.4)). Gpyro's critical temperature
$H(T-T_{\mathrm{crit},k})$ is again a **hard** Heaviside, where we use a tanh ramp.

The conversion formalism is worth dwelling on because it reveals when the two codes
agree. Gpyro normalises by the cumulative species mass
$(\bar\rho Y_{A_k}\Delta z)_\Sigma$ (the initial plus all subsequently formed mass
of $A_k$), so that each cell "behaves as if it was a TGA sample." When $n_k=1$ and
$n_{\mathrm{O_2},k}=0$, Gpyro's rate collapses (their Eq. 41) to

$$
r_k = \bar\rho\,Y_{A_k}\,Z_k\,e^{-E_k/RT} = \xi_{A_k}\,A_i\,e^{-E_i/RT},
$$

since $\bar\rho Y_{A_k} = \xi_{A_k}$ is exactly our mass concentration. **For
first-order, oxygen-insensitive kinetics with no specialised model, our (6.1) and
Gpyro's rate are identical** — the conversion bookkeeping reduces to plain mass
concentration. Gpyro's stoichiometry (solid fraction $\mathrm{SF}_k$ and the
gas-yield split, their Eqs. 23–24) is a parameterisation of the same mass-balanced
yields we store directly as $\nu_{i,j}$; Gpyro's $\chi_k$ swelling parameter
belongs to the volume-change model (§10), not the kinetics.

On the **heat of reaction**, Gpyro (Eq. 43) splits the enthalpy into a solid part
$\Delta H_{\mathrm{sol},k}$ and a volatilisation part $\Delta H_{\mathrm{vol},k}$
and forms $\dot Q_{s,k}''' = -\dot\omega_{fB,k}'''\Delta H_{\mathrm{sol},k} -
\dot\omega_{fg,k}'''\Delta H_{\mathrm{vol},k}$. Setting
$\Delta H_{\mathrm{sol},k}=\Delta H_{\mathrm{vol},k}=\Delta H_k$ recovers (their
Eq. 45) the per-reactant form $\dot Q_{s,k}''' = -\dot\omega_{dA_k}'''\Delta H_k$,
which is **exactly our** (6.8) with $\Delta H_k \equiv h_i$ and
$\dot\omega_{dA_k}''' \equiv r_i$. Gpyro's sign convention agrees with ours:
"a positive value of $\Delta H_{\mathrm{vol},k}$ denotes an **endothermic**
reaction, and a negative value … an exothermic reaction." Our single heat $h_i$ is
thus Gpyro's per-reactant-consumed heat (the $\Delta H_{\mathrm{sol}} =
\Delta H_{\mathrm{vol}}$ special case), which Gpyro itself documents as the way to
specify a heat "per kg of reactant consumed" — precisely our per-first-reactant
basis (§6.6.1). Gpyro additionally separates sensible-enthalpy transport into the
$\sum(\dot\omega_f - \dot\omega_d)h_i$ term of its gas energy equation; in
Pyrolysis.jl that bookkeeping is the gas-advection / generation source machinery of
§7, kept strictly separate from the chemical $h_i$.

### 6.9.3 FDS

The FDS solid-phase pyrolysis model (FDS *Technical Reference*, Eqs. 7.34–7.35)
maps most closely term-for-term. With $\alpha,\beta$ the FDS material/reaction
indices and $\rho_{s,\alpha}\to\xi_j$, $A_{\alpha\beta}\to A_i$,
$E_{\alpha\beta}\to E_i$, FDS writes

$$
r_{\alpha\beta} = \rho_{s,\alpha}^{\,n_{s,\alpha\beta}}\,
A_{\alpha\beta}\,e^{-E_{\alpha\beta}/RT_s}\,
[X_{\mathrm{O_2}}(x)]^{\,n_{\mathrm{O_2},\alpha\beta}}\,
T_s^{\,n_{t,\alpha\beta}},
$$

a product of a reactant-concentration power, the Arrhenius factor, an oxygen
power, **and a temperature power** $T_s^{n_t}$. Pyrolysis.jl reproduces the first
two factors (with the same mass-concentration basis — FDS notes $A_{\alpha\beta}$
has units 1/s for $n_s=1$, exactly our convention) and omits the oxygen and
explicit-temperature powers (§6.10). FDS's species production
$S_\alpha = \sum v_{\alpha,\alpha'\beta} r_{\alpha'\beta}$ (Eq. 7.37) and gas
production $\dot m_\gamma''' = \sum v_{\gamma,\alpha\beta} r_{\alpha\beta}$
(Eq. 7.38) with yields $v$ summing to 1 are our (6.5)–(6.6). FDS's chemical heat
source $\dot q_{s,c}''' = -\sum r_{\alpha\beta} H_{r,\alpha\beta}$ (Eq. 7.42) is
our (6.8) with $H_{r,\alpha\beta}\to h_i$, and FDS uses the **same** sign as Gpyro
and Pyrolysis.jl ($H_r > 0$ endothermic). FDS interprets each in-depth cell as a
virtual TGA sample, matching our cell-as-batch-reactor view.

As Gpyro notes, FDS Version 5's kinetics are the special case of Gpyro's Eq. 42
with $n_{\mathrm{O_2}}=0$, conventional reaction order, and no volume change — and
therefore also the special case of our (6.1) for first/$n$-th-order
oxygen-insensitive kinetics.

**Kinetic-parameter recovery from TGA.** A practical FDS capability worth
reproducing in user workflows is the closed-form estimate of $A_i, E_i$ from a TGA
peak. Given a reference (peak) temperature $T_p$, the normalised peak mass-loss
rate $r_p/Y_s(0)$ (s⁻¹), a constant heating rate $\dot T$ (K·s⁻¹), and residue
yield $v_s$, FDS uses (FDS *User Guide*, Eqs. 9.4–9.5; the Lyon–Stoliarov
relations)

$$
E_i = \frac{e\,r_p}{Y_s(0)}\,\frac{R\,T_p^2}{\dot T},
\qquad
A_i = \frac{e\,r_p}{Y_s(0)}\,e^{E_i/(R\,T_p)},
\qquad
\frac{r_p}{Y_s(0)} = \frac{2\dot T}{\Delta T}(1 - v_s),
$$

with $e$ Euler's number and $\Delta T$ the (triangular) pyrolysis-range width
(default 80 K), $\dot T$ default 5 K·min⁻¹. These produce a first-order
$(A_i, E_i)$ pair that, fed into our `Reaction(...)` constructor, reproduces the
TGA peak — a recommended starting point for cone/TGA calibration (see §18,
Verification & Validation, and the User Guide).

### 6.9.4 Summary of correspondences

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| Rate basis | mass conc. $\xi_j$ | mass conc. $\xi$ | $\bar\rho Y_i$ (= $\xi$) | $\rho_{s,\alpha}$ (= $\xi$) |
| Reactant factor | $\xi_j^{n_{i,j}}\tanh(\xi_j/\xi_{\mathrm{th}})$ (1 reactant) | $\xi_{\mathrm{C1}}\xi_{\mathrm{C2}}$ (bimolec.) | $(1-\alpha)^{n_k}$ + models | $\rho_s^{n_s}$ |
| Extra factors | — | — | $g(Y_{\mathrm{O_2}})$ | $X_{\mathrm{O_2}}^{n_{\mathrm{O_2}}} T_s^{n_t}$ |
| Temp. gate | smooth tanh | hard cutoff | hard $H(T-T_{\mathrm{crit}})$ | (none / power) |
| Stoich. yields | $\nu_{i,j}$, $\sum y=1$ | $\theta_i^j$ | $\mathrm{SF}_k$ / yields | $v$, $\sum v=1$ |
| Heat sign ($h>0$) | **endothermic** | **exothermic** | endothermic | endothermic |
| Heat source | $-\sum h_i r_i$ | $+h\,r$ | $-\sum \Delta H_k r_k$ | $-\sum H_r r$ |

The lone sign reversal (ThermaKin's $h>0$ exothermic) is the single most important
porting hazard; Gpyro and FDS agree with us.

---

## 6.10 Limitations and Open Issues

The kinetics model is deliberately the simplest one that covers the target
applications (cone calorimetry, TGA/DSC/MCC idealisation, charring/non-charring/
intumescent solids, moisture). Its boundaries, with the source-of-truth guards:

1. **Single reactant only.** `_apply_reaction_source!` errors for $N_{\mathrm{in}}
   > 1$ (src/Physics/kinetics.jl:183–185), and `mass_balance_error` returns `NaN`
   for multi-reactant reactions (src/Materials/reactions.jl:504). True bimolecular
   chemistry (ThermaKin's COMP1·COMP2, char oxidation consuming a gaseous O₂
   reactant) is not representable. Extending it requires a 3-ary reactant tuple
   `(idx, order, mass_stoich)` so a mass stoichiometry can be attached to each
   reactant; the outer struct shape already anticipates this
   (src/Materials/reactions.jl:43–49). Multi-step schemes are instead built as
   **sequences** of single-reactant reactions.

2. **No oxygen or pressure dependence.** The rate (6.1) depends only on $T$ and
   condensed-phase mass concentrations. There is no $g(Y_{\mathrm{O_2}})$ (Gpyro)
   or $X_{\mathrm{O_2}}^{n_{\mathrm{O_2}}}$ (FDS) factor and no partial-pressure
   term, so oxidative pyrolysis and pressure-enhanced decomposition cannot be
   modelled directly. Pressure effects enter only through the separate gas
   transport/permeability block (§9), never the rate constant.

3. **Constant pre-exponential and activation energy.** $A_i$ and $E_i$ are scalar
   `Float64` fields (src/Materials/reactions.jl:61–62). There is no
   temperature-dependent $A(T)$, no explicit $T^{n_t}$ power (FDS's fourth factor),
   no kinetic compensation rule $A(E)=A_0 e^{\beta E}$, and no isokinetic shift. A
   temperature-dependent rate constant would require wrapping the entire constant
   in a `CallableProperty` and providing a custom evaluation path — not part of the
   public API.

4. **No distributed-activation-energy (DAEM) tooling.** A DAEM model can be
   approximated by stacking several parallel `Reaction` instances with different
   $(A_i, E_i)$, but there is no built-in quadrature over an activation-energy
   distribution or fitting helper for it.

5. **Tanh gate width is fixed.** The temperature ramp (6.2) has a hard-coded unit-
   Kelvin argument scale, giving the ~10 K leakage of §6.2.3. It is not currently
   configurable; a sharper or softer gate would need a width parameter on the
   reaction.

6. **Depletion-limiter bias.** The limiter (6.3) under-predicts the rate near
   and below $\xi_{\mathrm{th}} = 1$ kg·m⁻³ (~24% at the threshold, growing as
   the reactant runs out) and converts exponential burnout tails into algebraic
   ones (§6.3.3). For materials whose physically meaningful residue
   concentrations are comparable to the default threshold, the threshold should
   be lowered (at some cost in stiffness) or the limiter disabled for
   validation.

These limitations are consistent across the kinetics and materials knowledge
bases; where a richer model is required (bimolecular oxidation, oxygen-sensitive
char burnout, DAEM), the present architecture's chief virtue is that the
allocation-free, AD-clean, compile-time-specialised kernels of §6.7 provide a
clean foundation onto which the additional factors can be grafted without
disturbing the surrounding finite-volume and Jacobian machinery (§13, §15).
