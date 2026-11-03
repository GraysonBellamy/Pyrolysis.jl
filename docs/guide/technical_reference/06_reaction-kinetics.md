# 6. Reaction Kinetics and Stoichiometry

This chapter develops the complete condensed-phase reaction model of Pyrolysis.jl:
the Arrhenius rate law on a mass-concentration basis, the smooth temperature gates
and the one-sided depletion limiter that stand in for the hard step functions of
classical pyrolysis codes, the stoichiometric assembly of species sources, and the
heat-of-reaction energy source with its (deliberately storage-oriented) sign
convention. Every governing equation, default value, and sign is derived from
first principles and then matched line-for-line to the implementation in
`src/Physics/kinetics.jl`, `src/Materials/reactions.jl` and the reaction kernel
`src/Kernels/reactions.jl`. We close with a careful contrast against ThermaKin2Ds,
Gpyro V2, and FDS, and a statement of the model's limitations.

Throughout this chapter the index `i` denotes the **reaction index** and `j` the
**component (species) index** — the kinetics convention of overload note I1 in
§2, Nomenclature. The cell index, where it appears in the discretization
discussion, is written explicitly. Mass concentration $\xi_j$ (kg·m⁻³) is the
constitutive species variable (the integrated state of a cell is its mass
$\xi_j V$ of each component, §6.5.4); stoichiometric coefficients are $\nu_{i,j}$ (overload note
A2); the Arrhenius pre-exponential is `A_i` (overload note A3); and the heat of
reaction `h` uses the storage sign convention `h > 0` endothermic (overload note
H1).

---

## 6.1 The Arrhenius Reaction Rate Law

### 6.1.1 Continuous form

A single condensed-phase decomposition reaction $i$ of reactant $a$ proceeds at a
volumetric rate $r_i$ (kg·m⁻³·s⁻¹) given by an Arrhenius rate constant
multiplied by a dependence on how much of the reactant is left and gated to a
temperature window:

$$
r_i(T,\boldsymbol{\xi}) \;=\; A_i \,
\exp\!\left(-\frac{E_i}{R_g\,T}\right)
\Gamma_i(\xi_a)
\;\cdot\; H(T - T_{\min,i})\,H(T_{\max,i} - T),
\qquad \xi_a > 0,
\tag{6.1}
$$

$$
\Gamma_i(\xi_a) =
\begin{cases}
b(\xi_a)\,g\!\bigl(\xi_a/\xi_{\mathrm{th}}\bigr), & n_i = 1,\\[2pt]
\xi_{\mathrm{ref},i}\;b_x(x_i)^{\,n_i}\,g\!\bigl(x_i/x_{\mathrm{th}}\bigr),
\quad x_i = \xi_a/\xi_{\mathrm{ref},i}, & n_i \ne 1,
\end{cases}
$$

and $r_i = 0$ when the reactant has $\xi_a \le 0$. Here

- $A_i$ is the pre-exponential factor, in s⁻¹ for every order (§6.4),
- $E_i$ is the activation energy (J·mol⁻¹),
- $R_g = 8.31446261815324$ J·mol⁻¹·K⁻¹ is the universal gas constant
  (`R_GAS`, src/Materials/constants.jl),
- $n_i$ is the reaction order ($n_i \ge 0$, enforced at construction),
- $\xi_{\mathrm{ref},i}$ is the **reference concentration** of the reactant, the
  condensed mass its lineage started with at the material point over the cell
  volume (§6.4, §6.5.5), and $x_i$ its **conversion**,
- $b(\xi_a)$ and $b_x(x_i)$ are the **finite bases**, equal to $\xi_a$ to within
  $10^{-8}$ kg·m⁻³ and to $x_i$ to within $10^{-8}$ (§6.3.2),
- $g$ is the one-sided **depletion onset** (§6.3), and
- $H(\cdot)$ is the Heaviside step, implemented as a **smooth tanh ramp** (§6.2).

At first order the reference cancels: $\xi_{\mathrm{ref}}\,x = \xi_a$, and that
case is evaluated on the concentration directly.

The factorisation of $r_i$ into a temperature-only Arrhenius constant
$k_i(T) \equiv A_i \exp(-E_i/R_g T)$ and a composition factor $\Gamma_i$ is exact
and is exploited directly by the code:

```julia
k = rxn.A * exp(-rxn.E / (R_GAS * T))
```

The rate constant carries the classical Arrhenius physical interpretation: $A_i$
is the collision/attempt frequency (modulated by an entropic steric factor), and
$\exp(-E_i/R_g T)$ is the Boltzmann fraction of molecular configurations that
surmount the activation barrier $E_i$. Both $A_i$ and $E_i$ are stored as plain
`Float64` scalars on the `Reaction` struct (src/Materials/reactions.jl), so
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
sum to unity (§6.5). ThermaKin and FDS also use mass-based
rates, but their factors and stoichiometric bases differ; Gpyro V2's conversion
rate includes a cumulative-mass normalization (§6.9).

### 6.1.2 Discrete (per-cell) form

The finite-volume discretisation (§13, Finite-Volume Discretization) collocates
all state at cell centres. The reaction kernel (`ReactionKernel`,
src/Kernels/reactions.jl) evaluates the reaction model once per cell, from the
cell's concentrations and trial temperature in the derived state (§5.13):

```julia
ξ = ntuple(j -> ds.ξ[j, i], Val(NC))
kernel.limiter.enabled && _all_quenched(ξ, material.reactions) && continue
T = ds.T_trial[i]
V = ds.V_eval[i]
ξ_ref = ntuple(r -> kernel.reference[r, i] / V, Val(NR))
evaluate_reactions!(c.rates, c.heats, ξ, T, material, kernel.limiter, ξ_ref)
species_source_terms!(c.species, c.rates, material)
for j = 1:NC
    c.S_mass[j, i] += c.species[j] * V
end
c.S_energy[i] += reaction_heat_source(c.rates, c.heats) * V
```

`ξ_ref` holds the reference concentrations. A first-order reaction does not read
its entry, and a reaction of another order with a zero reference (§6.5.5) writes
zero rate and zero heat without reaching `reaction_rate`, whatever the limiter or
the trial state.

The rates read the signed concentrations $\xi_{j} = m_{j}/V_{\text{eval}}$ of the
cell (not the positive parts the property laws see: the limiter itself makes a
rate vanish at and below zero, §6.3) and the trial temperature
$T_{\text{trial}}$, the cell temperature smoothly clamped into the property
bracket $[200, 3000]$ K, which equals the cell temperature on every physical
state and keeps $\exp(-E/RT)$ finite, with a finite derivative, at any trial
temperature. Because the cell masses are the states, the kernel multiplies the
volumetric sources by the cell volume: the mass sources `S_mass` are in kg·s⁻¹,
the energy source `S_energy` in W. A cell in which every reaction has a reactant
at or below zero contributes exact zeros and, with the limiter enabled, is not
evaluated at all (`_all_quenched`).

The continuous law (6.1) thus becomes simply $r_i^{\,\mathrm{cell}} =
r_i(T_{\mathrm{trial,cell}}, \boldsymbol{\xi}_{\mathrm{cell}})$ with the cell's own
reference — there is **no spatial operator inside the reaction term**. Reactions are local (point-wise in
space); all spatial coupling enters through transport (conduction §7, gas
diffusion/Darcy §9) and through the mesh motion (§11). This is the
cell-as-batch-reactor view: each control volume behaves, kinetically, like an
isolated TGA sample. This local evaluation is also used in Gpyro V2 and FDS,
although their oxygen factors and rate normalization differ (§6.9).

The single-reaction function is `reaction_rate` (src/Physics/kinetics.jl). It
is fully scalar-type-generic (`<:Real` arguments), so the dual numbers of the
colored Jacobian (§15) propagate cleanly through it.

> **Implementation note.** Single-reaction rate: `reaction_rate`,
> src/Physics/kinetics.jl. All-reaction evaluation: `evaluate_reactions!`
> dispatching to the `@generated` kernel `_evaluate_reactions_gen!`,
> src/Physics/kinetics.jl. Per-cell sources: `ReactionKernel`,
> src/Kernels/reactions.jl.

---

## 6.2 Smooth Temperature Gates

### 6.2.1 Why hard steps are unacceptable

The temperature window $[T_{\min,i}, T_{\max,i}]$ lets a reaction represent a
thermodynamically controlled phase transition that is switched on sharply at a
threshold, as the Stefan verification case does (§6.2.3). It is not needed to spare
the solver work: outside the window of a decomposition step the Arrhenius factor is
already negligible, the default gates are off at no cost, and the example decks run
ungated. The textbook expression of a window is the product of two Heaviside steps,
$H(T - T_{\min,i})\,H(T_{\max,i} - T)$, equal to 1 inside the window and 0 outside.

A hard step is, however, **non-differentiable** at the threshold and has a zero
derivative everywhere else. Both properties are poison for the implicit machinery
of Pyrolysis.jl:

1. The stiff time integrators linearize the residual with its Jacobian: the
   default `Rodas5P`, a Rosenbrock method, solves one linear system with it per
   stage, and the `FBDF` and `KenCarp4` policies iterate Newton's method with it
   (§15.5). Their accuracy and their step acceptance depend on a meaningful
   $\partial r_i/\partial T$. At a hard step the Jacobian is a Dirac delta
   (numerically: undefined), and away from it the step contributes nothing — the
   linearization sees a flat rate that suddenly jumps.
2. Automatic differentiation of a discontinuous function returns a derivative of
   zero almost everywhere and is simply wrong at the breakpoint, corrupting the
   colored Jacobian (§15); a finite-difference sensitivity (§17) across the
   breakpoint measures the jump instead of a derivative.

### 6.2.2 The tanh ramp

Pyrolysis.jl realizes each Heaviside as a hyperbolic-tangent ramp. The lower gate
is applied only when $T_{\min,i}$ is finite and strictly positive; the upper gate
only when $T_{\max,i}$ is finite:

$$
H_{\min,i}(T) =
\begin{cases}
1, & T_{\min,i}\ \text{not finite or}\ \le 0,\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh\!\big((T - T_{\min,i})/w_{\text{gate},i}\big)\bigr), & \text{otherwise},
\end{cases}
\qquad
H_{\max,i}(T) =
\begin{cases}
1, & T_{\max,i}\ \text{not finite},\\[4pt]
\tfrac{1}{2}\bigl(1 + \tanh\!\big((T_{\max,i} - T)/w_{\text{gate},i}\big)\bigr), & \text{otherwise},
\end{cases}
\tag{6.2}
$$

and the gated rate is $r_i \leftarrow r_i\, H_{\min,i}(T)\, H_{\max,i}(T)$, with
$w_{\text{gate},i}$ the per-reaction **gate width** (`Reaction(...;
gate_width)`, default $1.0$ K, validated finite and $>0$). The implementation
is verbatim:

```julia
if isfinite(rxn.T_min) && rxn.T_min > 0.0
    r *= 0.5 * (1 + tanh((T - rxn.T_min) / rxn.gate_width))
end
if isfinite(rxn.T_max)
    r *= 0.5 * (1 + tanh((rxn.T_max - T) / rxn.gate_width))
end
```

(src/Physics/kinetics.jl). The default gates are $T_{\min}=0$ and
$T_{\max}=\infty$ (src/Materials/reactions.jl), so an
unconstrained reaction incurs no gating cost at all — the `T_min > 0.0` guard
short-circuits the lower ramp and `isfinite(Inf) == false` short-circuits the
upper ramp.

The ramp $\tfrac12(1 + \tanh(\Delta T/w_{\text{gate}}))$ rises from $\approx 0$
to $\approx 1$ over a width set by $w_{\text{gate}}$. At the default
$w_{\text{gate}} = 1$ K the transition has a characteristic half-width of
$\mathcal O(1\text{–}3)$ K and reaches $\tfrac12(1+\tanh 5)\approx 0.99995$
only about 5 K inside the window; equivalently, the gate is
$\approx 4.5\times10^{-5}$ about 5 K **outside** the window. All widths scale
linearly with $w_{\text{gate}}$.

### 6.2.3 The ~$10\,w_{\text{gate}}$ leakage caveat

The smoothing is not free: the reaction is "weakly active" within roughly
$10\,w_{\text{gate}}$ of each boundary, where a hard step would have it
identically off. At the lower edge
$\tfrac12(1+\tanh((T - T_{\min})/w_{\text{gate}}))$ takes the value $0.5$
exactly at $T = T_{\min}$, $\approx 0.018$ at $T_{\min} - 2\,w_{\text{gate}}$,
and $\approx 4.5\times10^{-5}$ at $T_{\min}-5\,w_{\text{gate}}$. For
decomposition kinetics, where $T_{\min}$ is a numerical convenience rather
than a sharp physical threshold, this leakage is immaterial because the
Arrhenius factor is itself vanishingly small that far below the reaction's
effective onset. For a reaction used as a sharp phase-transition switch
(melting fronts, T-gated latent heat as in V-CPL-01b), narrow the gate
(`gate_width = 0.1` K makes the smear a tenth of the default) — the colored
Jacobian differentiates `reaction_rate` with dual numbers, so the width threads
through it automatically, at the cost of stiffer gate-adjacent Jacobian
entries. On a driving ΔT of only a few kelvin the default 1 K ramp visibly
shifts the effective switching point; either narrow the gate or offset
$T_{\min}$ to recentre the ramp. The Stefan verification case V-CPL-01b, which
represents solidification as a $T_{\max}$-gated reaction, measures the effect:
at the default 1 K width a 5 K driving temperature difference shifts the
effective freezing point by about 20 % of the driving force (the front runs
19 % slow), so the case drives with ±50 K, where the gate is a 1–2 % effect
(`test/verification/cases/V-CPL-01.jl`).

> **Implementation note.** Temperature gates: `reaction_rate`,
> src/Physics/kinetics.jl. `gate_width` field, validation, and the default
> `T_min = 0.0`, `T_max = Inf` in the convenience constructors,
> src/Materials/reactions.jl.

---

## 6.3 The Depletion Limiter

### 6.3.1 The stiffness problem near zero concentration

As a reactant is consumed, $\xi_j \to 0$. For a first-order reaction ($n_{i,j}=1$)
the bare rate $r_i = k_i(T)\,\xi_j$ approaches zero linearly, which is benign. But
the **Jacobian entry** $\partial r_i/\partial \xi_j = k_i(T)$ does not vanish: at
high temperature $k_i(T)$ can be enormous ($A_i$ ranges over $10^8$–$10^{20}$ s⁻¹
in typical fits). A control volume holding a numerically tiny,
possibly slightly negative residue of reactant therefore presents
the implicit solver with a stiff, fast-relaxing mode whose eigenvalue is
$-k_i(T)$ even though essentially no mass remains to react. This manifests as
time-step throttling and rejected steps, and — if $\xi_j$ overshoots below zero
and $n_{i,j}$ is non-integer — complex or `NaN` rates from $\xi_j^{n_{i,j}}$.

Sub-first-order fits make this strictly worse: for $0 < n_{i,j} < 1$ the bare
$\partial r_i/\partial\xi_j \sim \xi_j^{\,n-1}$ **diverges** at depletion, and for
a zero-order fit ($n_{i,j}=0$, common in TGA analysis) the bare rate does not
vanish at all — left unguarded it keeps consuming reactant at $k_i(T)$ straight
through $\xi_j = 0$ and integrates the concentration negative without bound.

For an order other than one the same singularity sits at $x \to 0$: on a moving
mesh the concentration stays at its bulk density while the mass is consumed, so
a scale on $\xi$ never engages, and the derivative of $x^{n}$ at burnout is what
has to be removed. The limiter of such a reaction therefore acts on the
conversion.

### 6.3.2 The one-sided rate factor

Pyrolysis.jl rolls the rate off smoothly as a reactant depletes and stops it
exactly at zero. The reactant contributes a multiplicative factor to the rate,
the power of a finite base times a one-sided **onset** applied to the *rate*,
not to the power-law base:

$$
r_i \;=\; k_i(T)\,f(\xi_a;\,n_i),
\qquad
f(\xi;\,n) =
\begin{cases}
0, & \xi \le 0,\\[2pt]
b(\xi)^{\,n}\; g\!\bigl(\xi/\xi_{\mathrm{th}}\bigr), & \xi > 0 .
\end{cases}
\tag{6.3}
$$

written here for the first-order concentration form; a reaction of another
order carries the same factor on its conversion, (6.3b) below. The **finite
base** is a softplus at the scale $\delta_b$ (`base_scale`, default
$10^{-8}$ kg·m⁻³),

$$
b(\xi) \;=\; \xi + \delta_b \ln\!\bigl(1 + e^{-\xi/\delta_b}\bigr),
$$

which exceeds $\xi$ by at most $\delta_b\ln 2 \approx 6.9\times10^{-9}$ kg·m⁻³
(at $\xi \to 0^+$) and equals $\xi$ exactly once $e^{-\xi/\delta_b}$ rounds
away; being positive, it keeps fractional powers real. The **onset** is

$$
g(x) =
\begin{cases}
0, & x \le 0,\\[2pt]
(s - 2y)\,t^3 + (3y - s)\,t^2,\quad t = x/x_J, & 0 < x < x_J,\\[2pt]
\tanh x, & x_J \le x < 20,\\[2pt]
1, & x \ge 20,
\end{cases}
\qquad y = \tanh x_J,\quad s = x_J\,\mathrm{sech}^2 x_J ,
\tag{6.3a}
$$

with the join $x_J$ (`join`, default 1, in units of $\xi_{\mathrm{th}}$). The
cubic has zero value and zero slope at $x = 0$ and matches $\tanh$ and its slope
at $x = x_J$, so $g$ is C¹ everywhere (`one_sided_factor`, src/Physics/kinetics.jl);
at the default join its coefficients are $3y - s = 1.8648$ and
$s - 2y = -1.1032$. The cut to exactly one at $x = 20$ is bit-identical to
evaluating the `tanh`, which rounds to `1.0` in Float64 from $x \approx 19.06$,
so the common path — bulk concentrations orders of magnitude above the
threshold — pays no transcendental. The implementation:

```julia
if limiter.enabled
    ξ <= 0 && return zero(conc) * zero(ξ)
    conc *= _finite_base(ξ, limiter.base_scale)
    conc *= one_sided_factor(ξ / limiter.threshold, limiter.join)
    return conc
end
return conc * _finite_power(ξ, limiter.base_scale, 1.0)
```

(`_reactant_factor`, src/Physics/kinetics.jl). At and below zero the factor and
every derivative of it are exactly zero, so an absent or negative trial
concentration is never consumed, whatever the reaction order: an overshooting
cell consumes nothing, and the reaction kernel skips a cell in which every
reaction lacks a reactant (§6.1.2).

Applying the onset to the rate rather than substituting an effective
concentration into the power-law base matters for low orders: a base-applied
limiter leaves a zero-order reactant unguarded ($\xi_{\mathrm{eff}}^{\,0}\equiv 1$)
and, for $0 < n < 0.5$, behaves as $(\xi^2/\xi_{\mathrm{th}})^n$ near zero with a
still-singular $\partial r/\partial\xi \sim \xi^{\,2n-1}$. The multiplicative
form removes the singular eigenvalue for **every** order $n \ge 0$: near
depletion (with $\delta_b \ll \xi \ll \xi_{\mathrm{th}}$ and the default join)

$$
r_i \;\approx\; 1.865\,k_i(T)\,\frac{\xi_j^{\,n+2}}{\xi_{\mathrm{th}}^2},
$$

so $\partial r_i/\partial\xi_j \to 0$ for every $n \ge 0$.

**On the conversion.** A reaction of order $n \ne 1$ carries the factor on
$x = \xi/\xi_{\mathrm{ref}}$, with its own dimensionless scales and the same
join and the same $g$:

$$
r_i \;=\; k_i(T)\,\xi_{\mathrm{ref}}\;b_x(x)^{\,n}\;g\!\bigl(x/x_{\mathrm{th}}\bigr),
\qquad
b_x(x) = x + \delta_x\ln\!\bigl(1+e^{-x/\delta_x}\bigr),
\tag{6.3b}
$$

exactly zero at and below zero mass under an enabled limiter
(`_conversion_factor`, src/Physics/kinetics.jl):

```julia
if limiter.enabled
    ξ <= 0 && return zero(conc) * zero(ξ)
    x = ξ / ξ_ref
    conc *= ξ_ref * _fast_pow(_finite_base(x, limiter.conversion_base), order)
    conc *= one_sided_factor(x / limiter.conversion_threshold, limiter.join)
    return conc
end
return conc * ξ_ref * _finite_power(ξ / ξ_ref, limiter.conversion_base, order)
```

With $C = (3\tanh x_J - x_J\,\mathrm{sech}^2 x_J)/x_J^2$, 1.8648 at the default
join, the rate vanishes in two ranges,

$$
r_i \approx C\,k_i\,\xi_{\mathrm{ref}}\,\frac{x^{\,n+2}}{x_{\mathrm{th}}^2}
\quad (\delta_x \ll x \ll x_{\mathrm{th}}),
\qquad
r_i \approx C\,k_i\,\xi_{\mathrm{ref}}\,(\delta_x\ln 2)^{n}\,\frac{x^{2}}{x_{\mathrm{th}}^2}
\quad (x \ll \delta_x),
$$

and $\partial r_i/\partial x \to 0$ in both. The factor is exactly one for
$x \ge 20\,x_{\mathrm{th}}$ when $x_J \le 20$; $g(10) = 1 - 4.1\times10^{-9}$;
the base correction rounds away above $x \approx 3.6\times10^{-7}$.

### 6.3.3 Behaviour and default values

The limiter is configured by the `DepletionLimiter` struct
(defined in src/Materials/reactions.jl; re-exported from `Pyrolysis.Physics`),
whose defaults are

| field | default | meaning |
|---|---|---|
| `threshold` ($\xi_{\mathrm{th}}$) | $1.0$ kg·m⁻³ | concentration scale of the roll-off |
| `join` ($x_J$) | $1.0$ | where the cubic onset joins `tanh`, in units of `threshold` |
| `base_scale` ($\delta_b$) | $10^{-8}$ kg·m⁻³ | scale of the finite power base |
| `enabled` | `true` | apply the onset; `false` keeps only the finite base |
| `conversion_threshold` ($x_{\mathrm{th}}$) | $10^{-3}$ | conversion scale of the roll-off, orders other than one |
| `conversion_base` ($\delta_x$) | $10^{-8}$ | conversion scale of the finite power base, orders other than one |

`threshold`, `join`, `base_scale`, `conversion_threshold` and `conversion_base`
must be positive; the constructor refuses anything else with an `ArgumentError`.
The join is shared, and the statement that the conversion onset is exactly one
above $20\,x_{\mathrm{th}}$ holds at $x_J \le 20$. Each `Material` carries its own limiter
in `material.depletion_limiter`, set with the `Material` constructor's
`depletion_limiter` keyword (default `DepletionLimiter()`); the reaction kernel
is built from it and passes it to every rate evaluation. `DEFAULT_DEPLETION_LIMITER = DepletionLimiter()`
(src/Physics/kinetics.jl) backs only the standalone `reaction_rate` entry point
when no `limiter` keyword is passed.

The factor has the key properties:

- **One-sided and C¹.** It is exactly zero at and below zero concentration, C¹
  through zero and through the join, and smooth everywhere else. Neither the
  rate nor $\partial r/\partial\xi$ jumps anywhere.
- $\xi \gg \xi_{\mathrm{th}}$: the bare rate is recovered double-exponentially —
  the relative deficit $1 - \tanh x \approx 2e^{-2\xi/\xi_{\mathrm{th}}}$ is
  $3.6\%$ at $2\,\xi_{\mathrm{th}}$, $4.1\times10^{-9}$ at $10\,\xi_{\mathrm{th}}$,
  and exactly zero in Float64 beyond $\approx 19.06\,\xi_{\mathrm{th}}$; the base
  differs from $\xi$ by less than $\delta_b$.
- $\xi \to 0$: the rate vanishes as $1.865\,k\,\xi^{\,n+2}/\xi_{\mathrm{th}}^2$
  with $\partial r/\partial\xi \to 0$ for every $n \ge 0$ — the stiff eigenvalue
  is removed for every order.

The robustness is bought with a rate **bias near depletion**: the suppression is
$24\%$ at $\xi = \xi_{\mathrm{th}}$ ($g(1) = \tanh 1 \approx 0.762$), $67\%$ at
$0.5\,\xi_{\mathrm{th}}$ ($g = 0.328$) and $98\%$ at $0.1\,\xi_{\mathrm{th}}$
($g = 0.0175$). For first-order kinetics the limited rate behaves as
$1.865\,k\,\xi^3/\xi_{\mathrm{th}}^2$ near zero, which converts the exponential
burnout tail into an algebraic one: integrating
$\mathrm d\xi/\mathrm dt = -1.865\,k\,\xi^3/\xi_{\mathrm{th}}^2$ gives
$\xi \approx \xi_{\mathrm{th}}/\sqrt{3.73\,k\,t}$ at late times. Choose
$\xi_{\mathrm{th}}$ well below the concentrations that matter for the quantities
of interest.

The conversion limiter's bias has the same shape in $x$: choose
$x_{\mathrm{th}}$ well below the conversions that matter. A reactant whose fit
holds down to the end of its conversion needs $x_{\mathrm{th}}$ below the residue
of interest; the default leaves the rate within $4.1\times10^{-9}$ of the bare law
down to $x = 10^{-2}$.

Reaction orders are validated at construction: $n < 0$ is rejected
(`_validate_reaction_order`, src/Materials/reactions.jl), since a negative-order
rate is singular at depletion and no limiter can repair it.

For validation studies that need the bare Arrhenius rate, the onset is
disabled per call:

```julia
using Pyrolysis
using Pyrolysis.Physics: reaction_rate

rxn = Reaction(:decomp, 1 => 2; A = 8.5e12, E = 188e3, h = 870e3)
r = reaction_rate(rxn, (2.0, 0.0), 650.0; limiter = DepletionLimiter(enabled = false))
```

or for an entire solve, per material:

```julia
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190.0, φ=.01, c = 1420.0, k = 0.21, ε = 0.86, λ=1e-5),
        GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02,),
    ),
    reactions = (Reaction(:decomposition, :virgin => :MMA; A = 8.5e12, E = 188e3,
                          h = 870e3),),
    depletion_limiter = DepletionLimiter(enabled = false),
)
```

Disabling keeps only the finite base: for $\xi > 0$ the factor is $b(\xi)^n$, and
for $\xi \le 0$ the power is formed in log space, where the value and its slope
underflow together (`_finite_power`, src/Physics/kinetics.jl). An absent reactant
then reacts at the rate of a $\delta_b$ concentration, a zero-order reaction
runs at its full rate at any concentration, and orders $0 < n < 1$ have a steep
$\partial r/\partial\xi$ near depletion. The guard is then the solver's domain
test, which rejects a proposed step that takes a cell's mass more negative than
its absolute tolerance (§15.5.4).

The approach to depletion is tested for orders 0, 0.5 and 1: a single cell
decays monotonically without its concentration going below `-1e-9`, the two
components summing to their initial total throughout
(`test/kernels/test_kinetics.jl`, "approach to depletion"; the limiter's C¹ join
and its finite derivatives at absent reactants are tested in the same file and
in `test/physics/test_kinetics.jl`).

> **Implementation note.** Limiter struct: `DepletionLimiter`,
> src/Materials/reactions.jl (carried per material in
> `material.depletion_limiter`; re-exported from Pyrolysis.Physics).
> Per-reactant application inside the rate
> (`_reactant_factor`, `_finite_base`, `_finite_power`, `one_sided_factor`):
> src/Physics/kinetics.jl.

---

## 6.4 Reaction Order and the Conversion Basis

**The order acts on the conversion.** A first-order reaction is evaluated on its
reactant's concentration. Every other order acts on the reactant's conversion

$$
x_i = \frac{m_{a}}{M_{\mathrm{ref},i}} = \frac{\xi_a}{\xi_{\mathrm{ref},i}},
\qquad
r_i V = k_i(T)\,M_{\mathrm{ref},i}\;x_i^{\,n_i}\,g\!\bigl(x_i/x_{\mathrm{th}}\bigr),
\tag{6.4}
$$

with $M_{\mathrm{ref},i}$ the reactant's **lineage inventory**, the condensed mass
its lineage started with at the material point (§6.5.5). The reason is the mesh.
A fit to thermogravimetric data defines the order on $m/m_0$. On a mesh that
moves with the solid the cell volume follows the condensed volume (§10), so a
pure reactant's concentration stays at its bulk density while its mass is
consumed: an order on the concentration would be a constant factor, and every
order would behave as first order. On the conversion the law is the same on any
mesh,

$$
\frac{m}{m_0} = \frac{L}{L_0} = \bigl[1-(1-n)\,k\,t\bigr]^{1/(1-n)}
\quad\text{(uniform slab, one reaction, no residue)},
$$

which case V-RX-06 asserts on fixed and moving meshes at $n = 0$, $0.5$ and $2$;
at $n = 0$ a burning slab thins linearly, $L = L_0(1-kt)$.

**Units of the pre-exponential.** $A_i$ is in s⁻¹ for every order: a conversion is
dimensionless and the rate carries one reference concentration, so
$k_i\,\xi_{\mathrm{ref}}$ is in kg·m⁻³·s⁻¹ at any $n$. A rate constant fitted to
$\mathrm dx/\mathrm dt = -k\,x^n$ is the deck's $A$ as it stands.

**Evaluation.** The order is a `Float64` and may be non-integer. The exponentiation
routes through a fast-path helper that special-cases the two most common orders:

```julia
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)
```

(src/Physics/kinetics.jl). For $n=0$ it returns `one(x)`: a zero-order reaction
runs at $k_i\,\xi_{\mathrm{ref}}$, and its onset still gates the rate to zero as
the reactant runs out (§6.3.2). Because the finite base is positive, $x^n$ is
always real-valued even for fractional $n$. `reaction_rate(rxn, ξ, T, ξ_ref =
NaN; limiter)` branches on the order before it reads the reference: the first-order
path never forms a conversion, and any other order requires a finite positive
`ξ_ref` and is refused without one.

**The reference.** $M_{\mathrm{ref},r,i}$ is computed at `setup` for every reaction of
order other than one (§6.5.5) and stays a constant of the cell: at such an order the
cells keep their material, since `setup` refuses any correction (§6.10). The reference concentration is $\xi_{\mathrm{ref},r,i} = M_{\mathrm{ref},r,i}/V_i$.
The reference adds no row to the state and no entry to the Jacobian's pattern.

> **Implementation note.** Detection: `needs_lineage`, src/Materials/materials.jl.
> The table of references: `ReactionKernel.reference`, src/Kernels/reactions.jl, from
> `lineage_references`, src/Integration/setup.jl.

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
Equations; in the discrete mass rows it enters times the cell volume, beside the
face fluxes (§6.5.4).

### 6.5.2 The single-reactant mass convention

A reaction has exactly one reactant. Per unit reaction
extent, that reactant is consumed with **unit mass stoichiometry** (an implicit
$\nu = -1$), and the products are produced with their explicit per-product mass
yields. The function that applies one reaction to the source vector is:

```julia
@inline function _apply_reaction_source!(src, rxn::Reaction, r)
    @inbounds src[rxn.reactant] -= r
    @inbounds for (idx, ν) in rxn.products
        src[idx] += ν * r
    end
    return nothing
end
```

(src/Physics/kinetics.jl). Note the asymmetry baked into the model: the
reactant carries a rate **order** $n_i$ (how its conversion scales the rate) but
its **mass** stoichiometry is fixed at $1$; a product carries a mass **yield** but
no rate order (products never appear in the rate's composition factor).

The coefficients are, then,

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
constructor. `_validate_stoichiometry` (src/Materials/reactions.jl) sums the
product yields and requires $\sum_p y_{i,p} = 1$ against the implicit unit reactant
mass. Each yield is additionally required to be non-negative (`_validate_yield` —
a pair like $(1.5, -0.5)$ sums to 1 but would silently drive the
negative-yield species' concentration negative without bound), and each
reaction order must satisfy $n \ge 0$ (`_validate_reaction_order`). The
tolerance is relative to the total coefficient mass with an absolute floor:

$$
\text{tol} = \max\!\bigl(10^{-10}\,(1 + \textstyle\sum_p y_{i,p}),\ 10^{-14}\bigr),
$$

(src/Materials/reactions.jl). A violation throws an `ArgumentError`.

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
species sum $\sum_j S_j$ is identically zero up to floating-point roundoff.

### 6.5.4 Gas-vs-solid dispatch

The species source $S_j$ is the same expression for every phase; the phase
distinction lives in the *conservation equation* into which $S_j$ is inserted, not
in the kinetics. The integrated states are the cell masses, and the row
assembly (`assemble_rows!`, src/Assembly/residual.jl) writes, for every
component of every cell,

$$
\frac{dm_{j,\mathrm{cell}}}{dt} = S_j\,V + \dot M_j ,
\qquad
\dot M_j = -\sum_f s_f\,A_f\,F_{j,f},
\tag{6.7}
$$

with $S_j V$ the reaction kernel's mass source [kg·s⁻¹] (§6.1.2) and
$\dot M_j$ the face-flux part of the mass rate, scattered from the face fluxes
$F_{j,f}$ the kernels wrote. The gas-transport kernel writes face fluxes for the
gas components only (its `gas` flags, `map(is_gas, material.components)`), so on
a mesh that rides the skeleton the solids and liquids have **no transport
term** — they move with the material point (§11) — while the gases carry their
diffusive and Darcy fluxes (§9). Where a correction moves the mesh relative to
the solid, the relative-flux kernel carries *every* component across the faces
at the relative velocity (§11.6). Because the masses are the states and every
face flux is added to one cell and subtracted from its neighbour, mass is
conserved exactly; the energy row assembled beside (6.7) is documented in §7,
Heat Conduction, Gas Convection, and Energy Assembly.

> **Implementation note.** Species source: `species_source_terms!` and the
> `@generated` kernel `_species_source_terms_gen!`, src/Physics/kinetics.jl;
> per-reaction application `_apply_reaction_source!` (same file); the cell
> sources `S_mass = S_j V`: `ReactionKernel`, src/Kernels/reactions.jl.
> Build-time balance: `_validate_stoichiometry`, src/Materials/reactions.jl.

### 6.5.5 Lineages: the reference of a reactant

The reference of a reaction of order other than one is its reactant's **lineage
inventory**, the condensed mass that can become the reactant at that material
point. A reactant present at the start is its own lineage; an intermediate's
inventory is what its parents can yield. Water is a lineage like any other.

**The graph.** The nodes are the condensed components. A reaction with condensed
reactant $a$, products $p$ of yield $y_p > 0$ and self-return $y_{aa}$ consumes
nothing when $y_{aa} \ge 1$ and contributes no edge; otherwise it is a **channel**
of $a$ with edge weights $w_{a\to p} = y_p/(1-y_{aa})$ to each condensed $p \ne a$.
A product listed with yield zero is not an edge.

**The rules**, run only on a material with an order other than one, in this order:

1. *External feed, overridden or not.* A condensed product of a reaction whose
   reactant is a gas is externally fed, and so is every node reachable from it;
   a solve whose sources add mass makes every condensed node externally fed. A
   reactant of order other than one that is externally fed is refused, and
   `basis` is no remedy: its lineage has no initial inventory. A gas reactant of
   order other than one is refused for the same reason.
2. *Ancestors.* Without `basis`, the reactant's ancestors (the nodes from which it
   is reachable along positive-yield edges) must hold no cycle; otherwise the
   reaction is refused, naming the cycle, with `basis` as the remedy.
3. *Reach.* The mass of $a$ that one unit of $v$'s inventory can become, in reverse
   topological order: $\Phi_a(a) = 1$, and for $v$, over the channels of $v$ that
   reach $a$, each channel's value is $\sum_p w_{v\to p}\,\Phi_a(p)$. Two reaching
   channels that differ by more than $10^{-12}$ relative are refused, since the
   branching would depend on the rates, with `basis` as the remedy. A channel that
   produces nothing toward $a$ is ignored: a competitive branch is normalised to
   the whole of its origin's inventory, so $A\to B$ beside $A\to C$ gives
   $\Phi_B(A) = y_B$. This is a derivation from the balance convention of the
   Kinetics Neo documentation, in which each initial reactant starts at formal
   concentration one and may react in several ways and an intermediate's balance
   sums its producing and consuming steps; it is not a quoted rule for arbitrary
   graphs.
4. *Per cell.* With $Y_{v,i} = m_{v,i}(0)/M_{0,i}$ and $M_{0,i}$ the cell's initial
   condensed mass (exactly zero in a cell that starts empty), the coefficient is
   $c_{a,i} = \sum_v Y_{v,i}\,\Phi_a(v)$ and the reference $c_{a,i}M_{0,i}$. A cell
   that holds the reactant and also an ancestor it can be formed from is refused
   unless `basis` is given. The maximum the reactant could ever reach there, the mass
   present plus what the ancestors can yield, is well defined, but no fitting software
   is known to normalise a node both present and formed that way, so the reference is
   left to `basis` rather than assumed.
5. *`basis`.* A reaction with `basis` $= b$ has $c_{r,i} = b$ on the cells holding an
   ancestor of its reactant and zero elsewhere, with no cycle or channel check.
   `basis` is a fraction of the material point's initial condensed mass.
6. *Representability.* A positive coefficient below $10^{-12}$ and a `basis` below
   $10^{-6}$ are refused, so every positive reference and every scale derived from
   it is a normal float, and exact zero means structural absence, never underflow.
7. *Inactive reactions.* A reaction whose reference is zero in every cell has no
   lineage: its reactant cannot acquire mass, it gets no row, and the kernel masks
   it (§6.1.2).

The graph rules (1 to 3) run in `Material(...)` (`lineage_graph`,
src/Materials/materials.jl). The rest need the initial fields or the solve's
sources and run in `setup` (`lineage_references`, src/Integration/setup.jl). A
material with every order equal to one is never validated, so its cyclic and
competitive networks are unaffected.

**A worked lineage.** The Douglas fir chain with 94.3 mg of wood and 5.7 mg of
water in a 100 mg cell, yields 0.95, 0.80, 0.37, 0.79 and 0.74, has the references
94.3, 89.585, 71.668, 26.517, 20.949 and 15.502 mg along the chain, and 5.7 mg for
the water. A step of order 2 holding 35 mg of the second intermediate runs at
$k\,(35)^2/71.668 = 17.09\,k$ mg/s, against $35\,k$ mg/s at first order.

**Importing fitted parameters.** A Kinetics Neo fit is on these lineage
normalisations: $A$ carries over unchanged. A fit normalised to the total initial
mass (the Aria/Sierra members of a parallel or series scheme) takes `basis = 1.0`,
which reproduces the bare and the limited law; the alternative
$A \leftarrow A\,c_a^{\,n-1}$ reproduces the bare law only. A fit on wet mass of a
dry-lineage reactant takes `basis = 1.0`. A Gpyro single step carries $A$ over; a
Gpyro series scheme with $n \ne 1$ on an intermediate differs by convention (§6.9.2).
An order other than one on an intermediate is a fitting convention in every one of
these codes, not a mechanism: the fraction reacted of a material that is generated and
consumed at the same time has no single value, which is why the ICTAC multi-step
recommendations reserve the pseudo-order form for initial reactants. Such parameters
are valid only with the normalisation of the software that produced them.

> **Implementation note.** Graph rules: `lineage_graph`, `LineageGraph`,
> src/Materials/materials.jl. Per-cell rules: `lineage_references`,
> src/Integration/setup.jl.

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

(src/Physics/kinetics.jl). The single minus sign is the entire content of
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

The source comment states this explicitly (src/Physics/kinetics.jl). In
the energy row of a cell $Q_{\mathrm{rxn}}$ enters, times the cell volume, as part
of the cell's energy source together with the absorbed radiation, beside the
conduction and the sensible enthalpy the mass fluxes carry across its faces;
the temperature form divides by the heat-capacity sum $C$ of every component
(§5.4, §7.3). For a cell of volume $V$ whose faces $f$ carry the conductive flux
$q_f$ and the mass fluxes $F_{j,f}$ at the face temperature $T_f$ ($s_f = \pm1$
the orientation of face $f$ relative to the cell, §2):

$$
C\,\frac{dT_{\mathrm{cell}}}{dt} =
\bigl(Q_{\mathrm{rxn}} + Q_{\mathrm{rad}}\bigr)V
- \sum_{f} s_f A_f\Bigl[q_f + \sum_j F_{j,f}\bigl(h_j(T_f) - h_j(T_{\mathrm{cell}})\bigr)\Bigr]
\tag{6.9}
$$

(`assemble_rows!`, src/Assembly/residual.jl; the energy-state alternative, the
enthalpy inventory, adds the formation term $\sum_j S_j V\,h_j(T_{\mathrm{cell}})$
instead, §7.3.3).

**Units check.** $[h_i] = $ J·kg⁻¹ (per kg of reactant) and $[r_i] = $
kg·m⁻³·s⁻¹, so $[h_i r_i] = $ J·m⁻³·s⁻¹ $=$ W·m⁻³, the volumetric power density
demanded by (6.9). The "per kg of reactant" basis is consistent with the
single-reactant mass convention of §6.5: the reactant is consumed at mass rate
$r_i$, and $h_i$ is the enthalpy per unit of that mass.

### 6.6.2 Temperature- and state-dependent heat

$h_i$ is not restricted to a constant. The `heat::Ph` field is typed by any
`AbstractPropertyFunction` (src/Materials/reactions.jl), and in
`evaluate_reactions!` the heat is evaluated with **both** $T$ and the full
concentration tuple $\boldsymbol{\xi}$:

```julia
heats[$i] = rxn.heat(T, ξ)
```

(src/Physics/kinetics.jl). The two-argument call dispatches as follows
(src/Materials/property_functions.jl):

- `ConstantProperty(h)` $\to h$ (the $(T,\xi)$ method ignores both arguments);
- `LinearProperty(a,b)` $\to a + bT$ (affine in temperature, $\xi$ ignored);
- `PolynomialProperty` $\to \mathrm{evalpoly}(T,\text{coeffs})$;
- `StateDependentProperty(f)` $\to f(T, \boldsymbol{\xi})$ — the only type that
  actually consumes $\boldsymbol{\xi}$.

A constant heat is created by passing a number, an affine heat by passing a
`(a,b)` tuple, through `to_property` (src/Materials/property_functions.jl):

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
\Bigl(2275 + 860\,e^{-0.120\,MC}\Bigr)\times 10^{3}\ \text{J}\,\text{kg}^{-1},
\qquad MC = 100\,\frac{\xi_{\mathrm{moist}}}{\sum_{k}\xi_k},
$$

encoded as a closure over the moisture index (src/Materials/property_functions.jl).
Because the closure is passed $\boldsymbol{\xi}$ at evaluation, the
sorption heat self-consistently rises as the wood dries — a behaviour impossible
to represent with a temperature-only property. A heat that depends on the composition
makes the enthalpy of the liquid depend on the drying path; the accounting stays exact
because the reaction runs one way, so the path is monotone. Under the local-temperature
convention of §6.6.3 the base of an evaporation heat is the latent heat at the evaporation
temperature, about 2.26 MJ·kg⁻¹ at 373 K, with the differential heat of sorption added to it. Calling a `StateDependentProperty`
with only $T$ raises an error (src/Materials/property_functions.jl),
guaranteeing the state argument is never silently dropped.

Whether the reaction heat is "separate from" the sensible-enthalpy transport of
the released gases is precisely the **datum question** of §6.6.3 — the partition
between "chemical" and "sensible" enthalpy depends on the temperature at which
$h_i$ is referenced.

> **Implementation note.** Heat source: `reaction_heat_source`,
> src/Physics/kinetics.jl. Heat evaluation with state:
> `_evaluate_reactions_gen!` (same file). Property dispatch:
> src/Materials/property_functions.jl.

### 6.6.3 The enthalpy datum: heats at reaction temperature vs the $T_{\text{ref}}$ datum

A heat of reaction is only defined relative to a **datum**: $h_i$ may be quoted

- **at the local reaction temperature** — $h_i(T) = \sum_j \nu_{i,j}\,h_j(T)$,
  the reaction-temperature convention used here; an imported DSC or fitted
  heat needs its units, reaction basis and datum checked explicitly (§6.9); or
- **at the fixed ambient datum** $T_{\text{ref}} = 298.15$ K —
  $h_i = \sum_j \nu_{i,j}\,h_j(T_{\text{ref}})$, the formation-enthalpy style.

The energy equation (§3.3) uses the first, as do ThermaKin (its temperature-form
balance), Gpyro (a heat of volatilization "at whatever temperature the reaction occurs")
and FDS (a heat of reaction at temperature): heats fitted in those codes, and heats read
from DSC peak areas, are supplied unchanged. A constant heat then assumes that the
reactant's and the products' heat capacities match,
$h_i(T) = h_i(T_{\text{peak}}) + \int_{T_{\text{peak}}}^{T}\Delta c_p\,\mathrm{d}T'$; the difference is of
order 50 kJ·kg⁻¹ for $\Delta c_p$ = 0.5 kJ·kg⁻¹·K⁻¹ over 100 K, and a temperature-dependent
heat (§6.6.2) carries it where it matters. With heats at the local reaction
temperature, subtracting the species equations from the total-enthalpy balance
gives the **exact** temperature equation
$\rho c_p\,\partial_t T = -\nabla\cdot q - \sum_g c_{p,g}\,\mathbf J_g\cdot\nabla T
+ Q_{\text{rad}} + Q_{\text{rxn}}$, whose second term is the temperature-form
reading of the gas enthalpy flux at the faces (§7.4) — no further term exists,
because the gas products' sensible enthalpy
$\Delta h_g(T, T_{\text{ref}}) = \int_{T_{\text{ref}}}^{T} c_{p,g}\,\mathrm{d}T'$ is
already *inside* $h_i(T)$. The reaction kernel evaluates the heats at the local
temperature (§6.1.2), consistent with this datum.

A heat quoted at the $T_{\text{ref}}$ datum is converted by Kirchhoff's law,
$h_i(T) = h_i(T_{\text{ref}}) + \sum_j \nu_{i,j}\,\Delta h_j(T, T_{\text{ref}})$,
with the sum running over the reactant ($\nu = -1$) and the products; the
converted heat can be supplied as a temperature-dependent heat property (§6.6.2).
The conversion is the integral of the heat-capacity change across the reaction,
$\int_{T_{\text{ref}}}^{T}\bigl(\sum_p y_{i,p}c_{p,p} - c_{p,\text{reactant}}\bigr)\,\mathrm{d}T'$.
The gas products' share alone, $\Delta h_g(900\,\text{K}, 298.15\,\text{K}) =
0.66$–$1.20$ MJ·kg⁻¹ for $c_{p,g} = 1.1$–$2.0$ kJ·kg⁻¹·K⁻¹, is largely offset by
the reactant's, but a heat-capacity change of 0.5 kJ·kg⁻¹·K⁻¹ still shifts the heat
by 0.3 MJ·kg⁻¹ at 900 K, comparable to typical wood/polymer pyrolysis endotherms
(0.1–0.9 MJ·kg⁻¹). A $T_{\text{ref}}$-datum heat used unconverted mis-states the
reaction endotherm by that integral.

---

## 6.7 Allocation-Free Evaluation

The reaction model sits in the innermost loop of every residual and Jacobian
evaluation, so it is engineered to allocate nothing and to specialise completely
on the (compile-time-known) reaction and component counts $N_R$ and $N_C$. Two
`@generated` kernels and one inline reduction achieve this.

**Rate and heat evaluation.** `evaluate_reactions!(rates, heats, ξ, T, material,
limiter = material.depletion_limiter, ξ_ref = NaNs)` forwards to
`_evaluate_reactions_gen!` (src/Physics/kinetics.jl), which the compiler expands
into a flat block — one explicit `rates[i] = reaction_rate(rxn, ξ, T, ξ_ref[i];
limiter)` and `heats[i] = rxn.heat(T, ξ)` statement per reaction — with no run-time
loop over the heterogeneous `reactions` tuple. A reaction of order other than one
with a zero reference writes zero rate and heat, and with no reference given it is
refused. Heterogeneity (each `Reaction` may have a different number of products and
heat property type) is therefore handled by static dispatch with no `Union` boxing
and no allocation.

**Species sources.** `species_source_terms!(sources, rates, material)` zeroes the
output and calls `_species_source_terms_gen!`, again unrolled to one
`_apply_reaction_source!` per reaction. It is the single method, called by the
reaction kernel (§6.1.2).

**Heat sum.** `reaction_heat_source` is a single `@inline` reduction over the
pre-filled `rates`/`heats` vectors (src/Physics/kinetics.jl), promoting the
accumulator type so dual numbers carry through.

All three use the promoted scalar type, so the colored Jacobian (§15) runs the
same code on dual numbers and derivatives propagate end-to-end. The per-cell
scratch (`rates` and `heats`, sized by $N_R$; `species`, sized by $N_C$) lives in
the kernel contributions (`Contributions.rates`, `.heats`, `.species`,
src/Kernels/contributions.jl) and is reused across cells and evaluations; the
cell sources go to `S_mass` and `S_energy` of the same buffer. A dual-number
evaluation uses the dual twin of the buffer (§7.7).

> **Implementation note.** Generated kernels: `_evaluate_reactions_gen!` and
> `_species_source_terms_gen!` (src/Physics/kinetics.jl). Compile-time constants
> $N_C, N_R$ are the `Material{NC,NR}` type parameters
> (src/Materials/materials.jl).

---

## 6.8 Validation

Stoichiometry errors are caught when a reaction is built (§6.5.3).
`validate_problem`, which `solve` always runs first, checks among other
things (the mesh, the boundary conditions, the time span, the initial
conditions, the initial porosity of a gas problem) that the material has
at least one component and that every reactant and product index lies in
$1..N_C$ (numeric reactions are not index-checked when the `Material` is built),
and checks initial physical pore volume and gas-boundary data (§9).

Reactions given by name are resolved to integer indices by `resolve_reaction`
against the names of the components when the `Material` is built, so a reaction
that names a non-existent component errors with a clear message before any solve.

> **Implementation note.** `_validate_stoichiometry`, `resolve_reaction`:
> src/Materials/reactions.jl. `validate_problem`:
> src/Problem/problem_validation.jl.

---

## 6.9 Comparison with reference solvers

Reaction parameters can be transferred only after the reaction rate and heat
have been put on the same mass basis. The versions in §1.6.5 differ in
reactant factors, conversion history, temperature gates and enthalpy
conventions. The table gives the main porting checks; the equations below
state the restricted cases in which a mapping exists.

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro V2 | FDS 6.11.0 solid wall |
|:---|:---|:---|:---|:---|
| Bare reactant factor | $\xi$ at first order; $\xi_{\mathrm{ref}}x^n$ on the lineage inventory otherwise | one or two concentrations | conversion model times cumulative available mass/volume | $\rho_{s,\alpha}^{n_s}=\xi_\alpha^{n_s}$ |
| Extra rate factors | smooth gates and depletion limiter | hard temperature limits | oxygen dependence and critical temperature | oxygen and temperature powers, rate limits |
| Positive heat | endothermic | exothermic | endothermic | endothermic |
| Match before comparing | unit reactant consumption | normalized stoichiometry and flipped heat sign | conversion normalization and heat basis | units, rate options and enthalpy datum |

### 6.9.1 ThermaKin2Ds

TK Eqs. 5–6 define $\theta_1 C_1+\theta_2 C_2\to\theta_3 C_3+
\theta_4 C_4+h$ and

$$r_{\mathrm{TK}}=A_{\mathrm{TK}}e^{-E/(RT)}\xi_{C_1}\xi_{C_2},$$

where the second concentration factor is set to one for a single reactant, so
the single-reactant TK rate is first order on the concentration and maps onto
ours directly. A TK second-order step in one species,
$A_{\mathrm{TK}}\xi^2$, is the conversion form with $A = A_{\mathrm{TK}}\xi_{\mathrm{ref}}$
at a fixed reference; the two laws then agree on a fixed mesh and differ on a
moving one, where our reference concentration follows the mesh volume and the
TK law has none. Pyrolysis.jl does not support the general bimolecular source. A sequence of single-reactant steps is not in general equivalent to
a bimolecular rate law.

TK applies hard lower/upper temperature limits. Pyrolysis.jl uses smooth
tanh gates and a concentration depletion limiter. Even when the bare rates
match, these modifiers can produce different near-onset and near-depletion
behavior. Smoothness is our implementation choice, not a prerequisite for
all implicit solvers or all uses of automatic differentiation.

TK's heat source is $+h_{\mathrm{TK}}r_{\mathrm{TK}}$, positive for an
exothermic reaction; ours is $-h r$. For a single-reactant TK step with
$\theta_1=1$, the conversion is $h=-h_{\mathrm{TK}}$. If
$\theta_1\ne1$, normalize the rate and yields too:
$r=\theta_1r_{\mathrm{TK}}$, $A=\theta_1A_{\mathrm{TK}}$,
$\nu_j=\theta_j/\theta_1$, and $h=-h_{\mathrm{TK}}/\theta_1$.
This is an algebraic mapping of TK Eqs. 5–6 to our unit-consumption source,
not an implementation of two-reactant chemistry. Pyrolysis.jl validates the
product-yield sum when a reaction is built.

### 6.9.2 Gpyro V2

For the default reaction-order model with cumulative reactant accounting,
V2 uses the mass $m_\Sigma$ initially present plus that subsequently formed
in the cell, and $1-\alpha_A=m_A/m_\Sigma$. In cell-volume notation,

$$
r_k=(1-\alpha_A)^{n_k}\frac{m_\Sigma}{V}
Z_k e^{-E_k/(RT)}g(Y_{O_2})\,s_k(T).
$$

This is the rate in GP Eqs. 35–40 and V2 user guide §11.4, with the
implemented temperature gate
$s_k(T)=\min[\max((T-T_{\mathrm{crit},k})/(1\,\mathrm{K}),0),1]$.
The manual states a cutoff below `TCRIT`; the
[`REACTION_RATE_T` implementation](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_mass.f90#L68-L125)
also ramps the rate over the next kelvin. For $n_k=1$, oxygen-independent
kinetics and a fully open gate, the bare rate reduces to
$\xi_A Z_k e^{-E_k/(RT)}$. It matches our unregularized first-order law,
subject to our own limiter and gate settings. For $n_k\ne1$ on a reactant
present from the start, $m_\Sigma$ is our reference and $1-\alpha_A$ our
conversion, so $Z_k$ is our $A$. On an intermediate the two differ by
convention: $m_\Sigma$ accumulates what has been formed so far, our reference is
what the lineage can yield, known from the start (§6.5.5). That cumulative form
is not implemented, and a Gpyro series scheme with $n_k\ne1$ on an intermediate
does not transfer. Other conversion models and the
`CONVENTIONAL_RXN_ORDER` option change this accounting
([`REACTION_RATE_Y`](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_mass.f90#L129-L259)).
Local oxygen dependence and gas-consuming yields cannot generally be
represented by our single-reactant source.

GP Eqs. 43–45 and V2 user guide §11.5 split reaction heat into
condensed-product and volatilized-mass parts. With `CONSTANT_DHVOL=.TRUE.`,
the equivalent per-reactant heat is
$h_k=\mathrm{SF}_k\Delta H_{\mathrm{sol},k}+
(1-\mathrm{SF}_k)\Delta H_{\mathrm{vol},k}$, after matching the sensible
enthalpy convention. Setting the two heats equal is one special case.
Positive heat is endothermic in both codes.

With `CONSTANT_DHVOL=.FALSE.`, the
[`CHEMICAL_SOURCE_TERM` implementation](https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/src/common/gpyro_thermal.f90#L717-L814)
also subtracts $\dot\omega_{fg,k}[h_g(T)-h_A(T)]$, using
$h_g(T)=c_{p,g}T$. V2 user guide Eq. 44 prints the enthalpy difference
in the opposite order under the same minus sign. Use the implemented
correction and its enthalpy datum when converting this option; the fixed
heat formula alone is insufficient.

### 6.9.3 FDS

The FDS 6.11.0 solid-phase chapter and `PYROLYSIS` routine use

$$
r_{\alpha\beta}=\rho_{s,\alpha}^{n_s}A_{\alpha\beta}
e^{-E_{\alpha\beta}/(RT_s)}X_{O_2}(x)^{n_{O_2}}T_s^{n_t}.
$$

Here $\rho_{s,\alpha}=\xi_\alpha=m_\alpha/V$ is the reactant's mass
concentration (§4.8.3). At $n_s=1$ this is our bare single-reactant law. At
$n_s\ne1$ the FDS order acts on the concentration and ours on the conversion,
so on a fixed mesh $A = A_{\alpha\beta}\,\xi_{\mathrm{ref}}^{\,n_s-1}$, with $A$ in
s⁻¹. Pyrolysis.jl omits
the explicit oxygen and temperature powers. FDS's oxygen profile is a
prescribed exponential penetration model
with surface oxygen determined by a mass-transfer balance, not a pore-gas
species-transport solution. FDS's heat source has our endothermic-positive
sign, but its temperature-dependent reaction enthalpy accounting must also
be matched (FDS Appendix P). Older FDS 5 rate normalization must not be
substituted for this release's formula.

**Kinetic-parameter recovery from TGA.** The FDS User's Guide's
`REFERENCE_TEMPERATURE` / `REFERENCE_RATE` prescription is a first-order
estimate. With reference rate $r_p/Y_s(0)$ [s⁻¹], temperature $T_p$ [K]
and heating rate $\dot T$ [K/s], its formulas are

$$
E=\frac{e r_p}{Y_s(0)}\frac{RT_p^2}{\dot T},\qquad
A=\frac{e r_p}{Y_s(0)}e^{E/(RT_p)}.
$$

In the 6.11.0 `READ_MATL` implementation, an omitted `REFERENCE_RATE`
is replaced by $2\dot T/\Delta T$ using `PYROLYSIS_RANGE`; this branch
does not multiply by a residue-yield factor. Defaults are 80 K and 5 K/min
([source](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90#L7335-L7351)). FDS input temperatures are in °C; these equations use K.
These are starting estimates, not an exact fit to a measured TGA peak.
User Guide §4.9 shows their use to construct a Pyrolysis.jl `Reaction`;
there is no automatic FDS input importer.

### 6.9.4 Sources and numerical evidence

TK Eqs. 1, 5–6 describe the documented rate and heat basis. GP §§3.1–3.5
and V2 user guide §11 are checked against the linked V2 rate and heat-source
routines. The FDS solid-phase reaction equations and Appendix P are checked
against [`PYROLYSIS`](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90#L3184-L3580) and the parameter
recovery path linked above. Our rate and source assembly are in
`src/Physics/kinetics.jl`.

V-RX-01 and V-RX-02 test specified reaction networks; the retained external
overlays have the version scope in §18.1.3. They do not establish agreement
for oxygen-dependent kinetics or arbitrary reaction-heat conventions.
User Guide §§4.4–4.9 covers the practical parameter and heat conversions.

---

## 6.10 Limitations and Open Issues

The kinetics model is deliberately the simplest one that covers the target
applications (cone calorimetry, TGA/DSC/MCC idealisation, charring/non-charring/
intumescent solids, moisture). Its boundaries, with the source-of-truth guards:

1. **Single reactant only.** A `Reaction` holds one reactant index and one order
   (src/Materials/reactions.jl). True bimolecular
   chemistry (ThermaKin's COMP1·COMP2, char oxidation consuming a gaseous O₂
   reactant) is not representable. Extending it requires each reactant to carry
   its own order and mass stoichiometry. Multi-step schemes are instead built as
   **sequences** of single-reactant reactions.

2. **No oxygen or pressure dependence.** The rate (6.1) depends only on $T$ and
   condensed-phase mass concentrations. There is no $g(Y_{\mathrm{O_2}})$ (Gpyro V2)
   or $X_{\mathrm{O_2}}^{n_{\mathrm{O_2}}}$ (FDS) factor and no partial-pressure
   term, so oxidative pyrolysis and pressure-enhanced decomposition cannot be
   modelled directly. Pressure effects enter only through the separate gas
   transport/permeability block (§9), never the rate constant.

3. **Constant pre-exponential and activation energy.** $A_i$ and $E_i$ are scalar
   `Float64` fields (src/Materials/reactions.jl). There is no
   temperature-dependent $A(T)$, no explicit $T^{n_t}$ power (FDS's fourth factor),
   no kinetic compensation rule $A(E)=A_0 e^{\beta E}$, and no isokinetic shift. A
   temperature-dependent rate constant would require wrapping the entire constant
   in a `CallableProperty` and providing a custom evaluation path — not part of the
   public API.

4. **No distributed-activation-energy (DAEM) tooling.** A DAEM model can be
   approximated by stacking several parallel `Reaction` instances with different
   $(A_i, E_i)$, but there is no built-in quadrature over an activation-energy
   distribution or fitting helper for it.

5. **Depletion-limiter bias.** At first order the limiter (6.3) under-predicts the rate near
   and below $\xi_{\mathrm{th}} = 1$ kg·m⁻³ (24 % at the threshold, 98 % at
   $0.1\,\xi_{\mathrm{th}}$) and converts exponential burnout tails into
   algebraic ones, $\xi \approx \xi_{\mathrm{th}}/\sqrt{3.73\,k\,t}$ for
   first-order kinetics (§6.3.3). For materials whose physically meaningful residue
   concentrations are comparable to the default threshold, the threshold should
   be lowered (at some cost in stiffness) or the limiter disabled for
   validation. The conversion limiter (6.3b) carries the same shape of bias below
   $x_{\mathrm{th}}$.

6. **One conversion per cell.** The rate at an order other than one reads one
   conversion per cell, so a cell must hold one material: one conversion cannot
   reproduce the sum of the rates of material at different conversions, and a heated
   column holds a range of them. At such an order `setup` therefore refuses any
   correction, which moves material between cells. Measured on `examples/06` at 72
   cells against the same deck without one, the correction alone moves the mass-loss
   rate by 2.2 % of the peak at $n = 0.5$, 8.8 % at $n = 0.25$ and 32 % at $n = 0$.
   No test at a merge can bound what follows a merge of live material: below order one
   a relative difference $d$ grows as
   $\mathrm{d}\ln(|d|/x)/\mathrm{d}t = k(1-n)x^{n-1}$. Deactivation is allowed, on a
   line, for a cell that vanishes: the threshold reads the reactive inventory
   (§14.3.2), the dead cell's reference leaves with it, and the tail burns at the
   receiver's conversion. A cell that would keep 1 % or more of its volume stays. The
   threshold alone bounds neither what the transfer does to the rate nor what it does
   to the column: at a receiver near the limiter's onset the net change of the rate is
   several times the dead cell's own (3.2 × at $n=0$, 7.4 × at $n=0.5$, at conversions
   $9.9\times10^{-5}$ and $2\times10^{-4}$), so accuracy here is a measurement on
   each deck. On the exact solution of V-BO-09 the rate
   error falls as $\varepsilon^n$ and as $1/N$, to $4.4\times10^{-4}$ of the peak at
   $n = 0.5$, 40 cells and $\varepsilon_m = 10^{-4}$. On `examples/06` at the default
   threshold, against the same deck without deactivation, $n = 0.5$ lands within 0.06 %
   of the peak at 72 cells and 0.03 % at 144, 0.01 % integrated at both. At $n = 0$ the
   mass-loss rate is a staircase of single-cell burnouts: the deck's steps lead its
   reference by up to 2 s, 15 % of the peak pointwise, 1.0 % integrated, with the mass
   curve $3\times10^{-4}$ of the initial mass away. That is inside what the mesh costs:
   4.3 % integrated and $6.9\times10^{-4}$ between 72 and 144 cells without
   deactivation (V-BO-12a, `burnout_order_zero.jl`).

7. **Networks with no fixed lineage.** A lineage fed through the gas phase, a
   solve with external mass sources, gas reactants at orders other than one and
   multi-dimensional meshes are refused for any order other than one
   (§6.5.5, §15.1).

8. **Prefer first-order steps.** A multi-step scheme of first-order reactions needs no
   reference, no extra state, and no refusal, and transfers between codes without a
   normalisation convention. Where a fit allows it, more first-order steps are the
   better choice over a smaller number of steps of other orders.

9. **Yields are those of the fit.** The mass yields are constants (§6.5.3), so a char
   yield measured at thermogravimetric heating rates is applied at the surface of a
   cone sample, where the material heats a decade faster. Wood char yields fall with
   heating rate and final temperature, 20–35 % at TGA rates against 10–15 % in fast
   pyrolysis, while secondary charring of volatiles retained in a thick sample raises
   them; the net uncertainty on a TGA-calibrated residue at cone scale is of the order
   of ±5 % of initial mass, above the repeatability of residue measurements. Only a
   competitive scheme can make the yield depend on conditions, at the cost of two more
   kinetic pairs per component; the yields should be read as valid at the conditions
   of the fit.

10. **The triplet is extrapolated.** Kinetic parameters fitted at TGA heating rates,
   where the steps peak at 600–750 K, are applied at cone surface temperatures up to
   900 K and at heating rates of 1–10 K s⁻¹. The rate's sensitivity to temperature is
   $\partial \ln r/\partial T = E/R_gT^2 = 0.03$–$0.05$ K⁻¹ for $E$ = 150–200 kJ mol⁻¹
   at 650–750 K, so a 10 K error in the temperature a cell is fed changes its rate by
   30–50 %: within the decomposition window the kinetics are only as good as the
   temperature.

11. **Evaporation has no equilibrium limit.** A liquid modelled as a reactant
   evaporates one way at $k(T)\,\xi_w$ whatever the vapour partial pressure (§1.8.1
   item 4). The step is admissible only while $k(T_{\rm ambient})\,t_{\rm test} \ll 1$ and
   the pore pressure stays near the one at which its onset was placed ($T_{\rm sat}$ is
   384 K at 1.5 bar and 393 K at 2 bar). The wood examples do not meet the first
   condition: their parameters give $k(300\ {\rm K})\,t = 0.35$ (`examples/07`) and
   $0.25$ (`examples/08`) over a 600 s test, so the cold interior dries before it
   reaches 330 K. Measured on `examples/07` against the same deck with the water step
   gated at 373 K, the back-face temperature differs by up to 20 K during the heating
   phase, while the mass loss rate and the surface temperature move by less than 1 %
   and 4 K. A parameter set with a negligible ambient rate and an onset near
   $T_{\rm sat}$ removes two thirds of the difference; the factor
   $\max(0, 1 - p_v/p_{\rm sat}(T))$ on the water step would remove the rest and hold at
   any pore pressure without a new material parameter.

These limitations are consistent across the kinetics and materials chapters; where a richer model is required (bimolecular oxidation, oxygen-sensitive
char burnout, DAEM), the present architecture's chief virtue is that the
allocation-free, AD-clean, compile-time-specialised kernels of §6.7 provide a
clean foundation onto which the additional factors can be grafted without
disturbing the surrounding finite-volume and Jacobian machinery (§13, §15).
