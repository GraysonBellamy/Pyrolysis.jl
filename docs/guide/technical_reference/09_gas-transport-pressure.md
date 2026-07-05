# 9. Gas-Phase Mass Transport and Pressure

This chapter derives the constitutive and discrete laws by which gaseous
decomposition products move through the condensed-phase matrix in Pyrolysis.jl.
Two transport mechanisms are treated: ordinary (Fickian) diffusion driven by a
volume-fraction / partial-pressure gradient, and Darcy advection driven by the
gas-phase pressure field in a porous medium; together they form the *Darcy–Fick*
flux. We establish the ideal-gas pressure closure in pore space, the Darcy
velocity and its analytic pressure derivatives, the surface mass-flux boundary
treatments, the gas-only finite-volume divergence that closes the species
equations, and the analytic Jacobian contributions of every flux term. Every
equation in this chapter has been checked against the implementing source and,
where the textbook form and the code differ, the chapter documents what the code
actually does.

Throughout, the spatial coordinate is the through-thickness axis $z$, with
$z=0$ the substrate (bottom, boundary id 2) and $z=L$ the exposed surface (top,
boundary id 1). The face-flux sign convention is **positive = transport in the
$+z$ direction** (bottom $\to$ top), and the discrete divergence
$(F_R-F_L)A/V$ is positive when the quantity flows *out* of a cell. In this
chapter the index $i$ denotes a **cell**, $j$ a **component/species**, and $f$ a
**face**; $L,R$ denote the left/right cells bounding a face. The primary species
state is the mass concentration $\xi_j$ [kg·m⁻³], not the mass fraction
$Y_j=\xi_j/\rho$.

---

## 9.1 The volume-fraction-gradient diffusive flux

### 9.1.1 Continuous form and physical origin

Gaseous decomposition products in a hot porous solid diffuse down their
concentration gradient. The natural driving potential for an ideal gas is not the
mass concentration $\xi_j$ but the *mass fraction* $Y_j=\xi_j/\rho_g$, because
Fick's first law for a dilute species in a carrier of density $\rho_g$ is

$$
J_j = -\rho_g\,\lambda\,\frac{\partial Y_j}{\partial z}
    = -\rho_g\,\lambda\,\frac{\partial}{\partial z}\!\left(\frac{\xi_j}{\rho_g}\right),
$$

where $\lambda$ [m²·s⁻¹] is the effective gas-transfer (diffusion) coefficient of
the porous medium (a tortuosity- and porosity-corrected molecular diffusivity;
see §5 *Effective Properties and Mixing Rules*, and §9.2 below). For an ideal gas
in the pore space, the gas density at fixed pressure scales as

$$
\rho_g \propto \frac{P}{T},
$$

so along an isobar $\rho_g \propto 1/T$ and the volume fraction (equivalently, the
mass-fraction-like ratio) of species $j$ obeys $\xi_j/\rho_g \propto \xi_j T$. Up
to the pressure scale, then,

$$
J_j = -\rho_g\,\lambda\,\frac{\partial}{\partial z}\!\left(\frac{\xi_j}{\rho_g}\right)
    = -\frac{\lambda}{T}\,\frac{\partial(\xi_j T)}{\partial z}.
$$

This is the **volume-fraction-gradient formulation** documented at the head of
`src/Physics/gas_transport.jl`. Expanding the product derivative,

$$
\boxed{\;
J_j = -\lambda\left(\frac{\partial \xi_j}{\partial z}
       + \frac{\xi_j}{T}\,\frac{\partial T}{\partial z}\right)
\;}
$$

exposes the two physical contributions:

- **Ordinary diffusion** $-\lambda\,\partial\xi_j/\partial z$: transport down the
  composition gradient.
- **Thermal coupling** $-\lambda(\xi_j/T)\,\partial T/\partial z$: a Soret-like
  term arising *not* from genuine thermodiffusion but from the temperature
  dependence of the ideal-gas density. A hotter cell holds the same mass
  concentration in a more dilute gas (lower $\rho_g$), so equal $\xi_j$ across a
  temperature gradient still corresponds to a volume-fraction gradient and hence a
  flux.

Casting the driving potential in terms of $\xi_j T$ rather than $Y_j$ avoids ever
computing the gas density $\rho_g$ at the face, which would require evaluating $P$
and $T$ at an intermediate point. The formulation is well behaved in the
low-density (high-temperature, low-$P$) regime where $\rho_g\to0$ would otherwise
appear in a denominator.

### 9.1.2 Discrete face flux

Discretizing $-\frac{\lambda}{T}\,\partial(\xi T)/\partial z$ at an interior face
$f$ between cells $L$ and $R$, with the temperature in the prefactor taken at the
arithmetic-mean face value $T_{\text{face}}=\tfrac12(T_L+T_R)$, gives the
implemented kernel:

$$
\boxed{\;
J_j = -\,\lambda_{\text{face}}\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\;d_{LR}}
\;}
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

where $d_{LR}=|z_R^c - z_L^c|$ is the distance between the two cell centers and
$\lambda_{\text{face}}$ is the distance-weighted harmonic mean of the cell
transfer coefficients (§9.1.4). This is exactly

```julia
@inline function gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)
    T_face = (T_L + T_R) / 2
    return -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
end
```

*Implementation note.* Implemented in `gas_flux`, `src/Physics/gas_transport.jl`
(lines 36–47), and assembled face-by-face by `compute_gas_fluxes!` (with the
compile-time-unrolled inner helper `compute_gas_fluxes_unrolled!`).

### 9.1.3 The exact product-rule identity

The composition / thermal split of §9.1.1 is not an approximation of the discrete
kernel — it is algebraically identical to it when arithmetic-mean face values are
used. The identity

$$
\xi_R T_R - \xi_L T_L
   = T_{\text{face}}\,(\xi_R - \xi_L) + \xi_{\text{face}}\,(T_R - T_L),
\qquad \xi_{\text{face}}=\tfrac12(\xi_L+\xi_R),
$$

is proved by expanding both terms on the right:

$$
\begin{aligned}
T_{\text{face}}(\xi_R-\xi_L) + \xi_{\text{face}}(T_R-T_L)
&= \tfrac12(T_L+T_R)(\xi_R-\xi_L) + \tfrac12(\xi_L+\xi_R)(T_R-T_L) \\
&= \tfrac12\big(2\,\xi_R T_R - 2\,\xi_L T_L\big)
 = \xi_R T_R - \xi_L T_L .
\end{aligned}
$$

Dividing by $T_{\text{face}}\,d_{LR}$ and multiplying by $-\lambda_{\text{face}}$
recovers the discrete composition + thermal split

$$
J_j = -\lambda_{\text{face}}\!\left(\frac{\xi_R-\xi_L}{d_{LR}}
       + \frac{\xi_{\text{face}}}{T_{\text{face}}}\frac{T_R-T_L}{d_{LR}}\right).
$$

This equivalence is verified to machine precision (`atol = 1e-14`) in
`test/discretization/test_gas_flux.jl` (lines 99–113).

### 9.1.4 Why the harmonic mean is applied to $\lambda$, not $\lambda\rho$

At a material interface, *flux continuity* — not property continuity — is the
physical constraint. Treating the two half-cells as transport resistances in
series, the face transfer coefficient is the distance-weighted harmonic mean

$$
\boxed{\;
\lambda_{\text{face}} = \frac{d_L + d_R}{\dfrac{d_L}{\lambda_L} + \dfrac{d_R}{\lambda_R}}
\;}
$$

where $d_L,d_R$ are the distances from the two cell centers to the face. On a
uniform mesh ($d_L=d_R$) this reduces to the standard harmonic mean
$2\lambda_L\lambda_R/(\lambda_L+\lambda_R)$. The distance weighting makes the
scheme correct on stretched meshes (e.g. surface-refined grids); see §13.6 for
the general resistance-in-series argument and §7.1 for the identical treatment of
thermal conductivity. A small-value guard returns $0$ if either $\lambda_L$ or
$\lambda_R$ is below $10^{-20}$.

A subtle but important modeling choice: the code takes the harmonic mean of
$\lambda$ **alone**, not of the product $\lambda\rho_g$ that would appear in a
true mass-fraction-gradient discretization $J=-\lambda\rho_g\,\partial Y/\partial z$.
This is consistent with — indeed *required by* — the volume-fraction formulation,
in which $\rho_g$ never appears explicitly: the density dependence is folded into
the $\xi T$ driving potential, which uses *cell-center* concentrations $\xi_L,\xi_R$
and the arithmetic-mean face temperature, with no face density. No
$\lambda\rho$-averaging helper exists in the package; the distinction is
exercised in `test/discretization/test_gas_flux.jl`.

*Implementation note.* `λ_face` is precomputed once per step in the
`PropertyCache` via `harmonic_mean_weighted`,
`src/Properties/effective_properties.jl` (lines 1154–1159); per-cell $\lambda_i$
comes from `effective_gas_transfer` (§9.2 and §5).

### 9.1.5 Partial-pressure equivalence

Because the formulation is anchored to the ideal-gas relation, the same flux can
be written as a partial-pressure gradient. With $P_j = \xi_j R_g T / M_j$ the
partial pressure of species $j$,

$$
J_j = -\frac{\lambda_{\text{face}}\,M_j}{R_g\,T_{\text{face}}}\,\frac{\partial P_j}{\partial z}.
$$

To see the equivalence, note $\xi_j T = P_j M_j / R_g$, so
$\xi_R T_R - \xi_L T_L = (M_j/R_g)(P_{j,R}-P_{j,L})$, and substituting into the
discrete kernel gives the boxed pressure form directly. This is verified to
`atol = 1e-12` in `test/discretization/test_gas_flux.jl` (lines 119–142). The
partial-pressure reading reinforces the physical picture: gas diffuses down its
own partial-pressure gradient, with mobility $\lambda M_j/(R_g T)$.

---

## 9.2 One transfer coefficient for all species

Every gas component diffuses with the **same** effective transfer coefficient
$\lambda_{\text{eff}}$: the face coefficient $\lambda_{\text{face}}$ entering
the kernel of §9.1.2 is a single mixture property (§5.6), shared by all species
at that face. The model therefore does not discriminate between light and heavy
gases within one simulation — a light product (H₂O) and a heavy one (MMA)
diffuse at the same rate through the same matrix. The per-component `λ_j(T)`
property controls how each *condensed/gas phase* contributes to the mixture
coefficient as composition evolves, not how individual gas species move
relative to one another.

Where a species-faithful magnitude is wanted, the Chapman–Enskog helper
`chapman_enskog_diffusivity(T, M; P, σ, ε_k)` (§5.10) estimates the binary
diffusivity of a given vapor in air; the user can fold that estimate (with a
porosity/tortuosity correction) into the component `λ` values. Genuine
multicomponent (Stefan–Maxwell) diffusion, or even a per-species kinetic-theory
scaling of the diffusive flux, is not implemented (§9.12).

---

## 9.3 The Darcy–Fick combined flux

When the material has nonzero permeability, pressure-driven (Darcy) advection
augments diffusion. The total face flux is the sum of an advective and a diffusive
part:

$$
\boxed{\;
J_j^{\text{DF}} = \underbrace{\frac{\xi_j^{\text{upwind}}}{\varphi_{\text{face}}}\,v_g}_{J_j^{\text{adv}}}
   \;+\; \underbrace{\Big(\!-\lambda_{\text{face}}\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d_{LR}}\Big)}_{J_j^{\text{diff}}}
\;}
$$

where $v_g$ is the Darcy velocity at the face (§9.5, positive in $+z$). The
$1/\varphi_{\text{face}}$ is forced by the bookkeeping: $\xi_j$ is mass per
**bulk** volume (the pressure closure $P=\sum\xi_jR_gT/(M_j\varphi)$ divides by
$\varphi$ to recover the intrinsic pore concentration) while the Darcy velocity
is **superficial** (volumetric flux per total cross-section), so the mass flux
per bulk area is intrinsic density × superficial velocity,
$J^{\text{adv}}=(\xi/\varphi)\,v_g$ — equivalently $\xi\,v_i$ with the
interstitial velocity $v_i=v_g/\varphi$. The face porosity is the arithmetic
mean $\varphi_{\text{face}}=\tfrac12(\varphi_L+\varphi_R)$, floored at
$10^{-12}$ consistently with the pressure closure. The discrete kernel is

```julia
@inline function gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)
    T_face = (T_L + T_R) / 2
    J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
    φ_face = max((φ_L + φ_R) / 2, _PHI_FLOOR)
    ξ_upwind = v_g >= 0 ? ξ_L : ξ_R
    J_adv = v_g * ξ_upwind / φ_face
    return J_adv + J_diff
end
```

### 9.3.1 First-order upwind selection

The advected concentration is taken from the **upwind** cell — the cell from which
the flow originates:

$$
\xi_j^{\text{upwind}} = \begin{cases}
\xi_{j,L}, & v_g \ge 0 \quad (\text{flow toward }+z), \\
\xi_{j,R}, & v_g < 0 \quad (\text{flow toward }-z).
\end{cases}
$$

This first-order upwind rule is unconditionally stable and conservative (the same
face value enters the two cells sharing the face with opposite divergence signs).
It is total-variation-diminishing but only first-order accurate; no TVD/MUSCL
slope limiting is applied to the Eulerian Darcy–Fick advective term. (Higher-order
limited advection *is* available in the ALE relative-velocity flux; see §11.6–11.7,
*The ALE Framework*, and §13.4–13.5 for the limiter library.) A practical
consequence: when upwinding selects a cell with zero concentration adjacent to a
cell with positive concentration, the advective contribution can momentarily
produce a "wrong-side" zero flux; this is conservative but, over many steps,
relies on reaction and heat coupling to keep states physical.

### 9.3.2 Advection carries all species equally

The Darcy–Fick flux uses the bulk velocity $v_g$ in the advective term for every
component, and the diffusive term likewise uses one shared $\lambda_{\text{face}}$
for all species (§9.2). Physically, bulk pore flow is a mixture motion; the model
applies no molecular-weight discrimination in either term (§9.2, §9.12).

In the pure-Fickian limit $v_g\to0$, the Darcy–Fick flux reduces exactly to the
diffusion-only flux of §9.1 — the two assembly paths agree, which is checked in
the gas-flux test suite.

*Implementation note.* `gas_flux_darcy_fick`,
`src/Physics/gas_transport.jl` (lines 74–92); assembled by
`compute_gas_fluxes_darcy_fick!`.

---

## 9.4 Ideal-gas pressure in a porous medium

The pressure that drives Darcy flow is the gas-phase pressure inside the pore
space. Summing partial pressures of the gaseous components and dividing by the
porosity (the fraction of bulk volume actually occupied by gas) gives the
porous-medium ideal-gas law:

$$
\boxed{\;
P = \frac{1}{\varphi}\sum_{j:\,\text{is\_gas}(j)} \frac{\xi_j R_g T}{M_j}
  \equiv \frac{N(\xi,T)}{\varphi}
\;},
\qquad N \equiv \sum_{j:\,\text{gas}} \frac{\xi_j R_g T}{M_j}.
$$

The numerator $N$ is the sum of partial pressures *referred to the bulk volume*
(mass per bulk volume times $R_g T/M_j$). Dividing by the porosity $\varphi$
converts to the pressure *in the gas phase*, which occupies only the fraction
$\varphi$ of the bulk: the same gas mass confined to a smaller void volume exerts
a proportionally higher pressure. Without the $1/\varphi$ divisor the pressure
would be underestimated in nearly-solid cells.

Two implementation details are load-bearing:

- **Only true gases contribute to $N$.** A component contributes iff
  `is_gas(j)` is true *and* `molar_mass(j) > 0`; solids and liquids are skipped.
  This is a compile-time-eliminated branch in the generated sum
  `_gas_partial_pressure_sum`.
- **Porosity floor.** The divisor is clamped from below at
  $\varphi_{\min}=10^{-12}$ via `max(φ, _PHI_FLOOR)`, guaranteeing a finite
  pressure even if a cell becomes fully solid ($\varphi\to0$). This matches the
  clamping in the porosity routine (porosity is itself clamped to $[0,1]$).

The porosity $\varphi$ is computed from *intrinsic* (skeletal) component densities
$\rho_{i,j}$, not bulk densities:
$\varphi = \mathrm{clamp}\!\big(1-\sum_{j\in\text{S/L}}\xi_j/\rho_{i,j},\,0,\,1\big)$
(see §5.3). It is the same porosity used by the Carman–Kozeny permeability closure
(§9.5.2), so the pressure and the permeability share a single consistent void
fraction.

*Implementation note.* `ideal_gas_pressure`,
`src/Physics/pressure.jl` (lines 60–68); field assembly via `compute_pressure!`
(lines 77–89). Constant $R_g = 8.314462618$ J·mol⁻¹·K⁻¹.

---

## 9.5 Pressure gradient and Darcy velocity

### 9.5.1 Pressure gradient

At an interior face the pressure gradient is a plain central difference of
cell-centered pressures:

$$
(\nabla P)_f = \frac{P_R - P_L}{d_{LR}}.
$$

Boundary faces (one neighbor id $=0$) are assigned a zero gradient here; the
physical pressure boundary conditions are applied separately to the Darcy velocity
(§9.5.3, §12.5). The sign convention is the usual one: $(\nabla P)_f>0$ when
pressure increases toward $+z$.

*Implementation note.* `compute_pressure_gradient!`,
`src/Physics/pressure.jl` (lines 112–125).

### 9.5.2 Darcy velocity

Darcy's law gives the superficial (bulk) gas velocity at a face from the pressure
gradient, the permeability, and the gas viscosity:

$$
\boxed{\;
v_g = -\frac{\kappa_{\text{face}}}{\mu_{\text{face}}}\,(\nabla P)_f
    = -\frac{\kappa_{\text{face}}}{\mu_{\text{face}}}\,\frac{P_R - P_L}{d_{LR}}
\;}
$$

with $v_g>0$ for flow toward $+z$ (the surface). The face properties combine the
two cells consistently with their physical roles:

- **Mobility — distance-weighted harmonic mean.** Permeability is a flow
  *resistance* in series across the interface. The production path
  (`update_darcy_properties!`) forms the per-cell mobility
  $\text{mob}=\kappa/\mu$ and combines the two sides with the
  distance-weighted harmonic mean

  $$
  \text{mob}_{\text{face}} = \frac{d_L+d_R}{d_L/\text{mob}_L + d_R/\text{mob}_R},
  $$

  which reduces to the unweighted harmonic mean $2ab/(a+b)$ on uniform meshes.
  The harmonic mean is conservative: it tends toward the *smaller* of the two
  mobilities, so a single tight (low-$\kappa$) cell correctly throttles the
  interface flow. If **either** side's mobility is below $10^{-30}$ the face
  mobility is **exactly 0** — an impermeable cell is an infinite series
  resistance, and gas cannot cross into it. The Jacobian operator's
  `_dt_mobility_face` mirrors this rule.

- **Viscosity — per-cell temperature.** Viscosity is an intensive,
  smoothly varying gas property, so it is evaluated per cell at the cell
  temperature, $\mu_i=\mu_g(T_i)$, and enters the face through the per-cell
  mobilities $\kappa_i/\mu_i$ that the distance-weighted harmonic mean
  combines, using Sutherland's law (§5.9):

  $$
  \mu_g(T) = \mu_{\text{ref}}\left(\frac{T}{T_{\text{ref}}^{\mu}}\right)^{3/2}
             \frac{T_{\text{ref}}^{\mu}+S}{T+S},
  $$

  with the canonical air constants $\mu_{\text{ref}}=1.716\times10^{-5}$ Pa·s,
  $S=110.4$ K. **Caveat on the reference temperature:** the Sutherland reference
  in the code is $T_{\text{ref}}^{\mu}=273.15$ K (the temperature at which
  $\mu_{\text{ref}}$ is quoted for air), which is *distinct* from the enthalpy
  datum $T_{\text{ref}}=298.15$ K used elsewhere in the model. The two should not
  be conflated.

- **Guards.** A cell with $\kappa_i\le0$ or $\mu_i\le10^{-20}$ Pa·s has zero
  mobility, which zeroes the mobility of both adjacent faces (a fully solid or
  numerically degenerate cell carries no Darcy flow).

*Implementation note.* The production face mobility and velocity are
assembled in `update_darcy_properties!`
(`src/Discretization/finite_volume.jl`); `effective_gas_viscosity`,
`src/Properties/effective_properties.jl`. The `darcy_velocity_face` /
`harmonic_mean_permeability` family in `darcy_velocity.jl` is a test-only
alternative path — it never runs in production.

The effective permeability $\kappa$ entering the cell-center field is either a
component mixing rule (PARALLEL/SERIES/WEIGHTED on the renormalized solid
basis, with the series bound exactly 0 when any condensed phase is
impermeable; §5.7.1) or the porosity-based **Carman–Kozeny** closure

$$
\kappa = \min\!\left(\frac{\ell^2\,\varphi^3}{180\,(1-\varphi)^2},\;\ell^2\right),
$$

with $\ell$ the characteristic pore/particle diameter and $180$ the Kozeny
constant for a packed bed of spheres; as $\varphi\to1$ the packed-bed form
diverges and $\kappa$ is capped at the free-flow scale $\ell^2$ (binding for
$\varphi\gtrsim0.93$), with a $10^{-30}$ denominator guard solely against
$\mathrm{Inf}$ at exactly $\varphi=1$. Carman–Kozeny is
evaluated on the *canonical* porosity $\varphi(\xi,T)$ — the same porosity the
pressure formula uses — not on bulk-density volume fractions (see §5.7). The two
closures (component mixing vs. Carman–Kozeny) are mutually exclusive, and
BRUGGEMAN EMT is rejected for permeability (no physical basis for a transport
resistance; an `ArgumentError` is thrown). As virgin material converts to porous
char, $\varphi$ rises and $\kappa$ increases, easing gas escape.

### 9.5.3 Pressure boundary conditions and the Darcy venting flux

Pressure BCs do not enter the diffusive surface flux; they set the boundary-face
Darcy velocity, which drives the **advective venting flux into the species
equations**. The chain is: `boundary_darcy_velocity` (one shared scalar kernel
per BC type) → `apply_pressure_bc_to_darcy_velocity!` writes $v_b$ at the
boundary faces → `compute_gas_fluxes_darcy_fick!` writes the boundary face flux

$$
J_j = v_b\,\frac{\xi_j^{\text{int}}}{\varphi_{\text{int}}}
\quad\text{(outflow only)},
$$

upwinding against a zero-concentration ambient ghost: only outflow advects the
interior gas, inflow advects nothing (the ambient carries none of the tracked
species). The diffusive film-law part of the boundary flux is then **added** by
`apply_mass_bc_flux!`. Venting therefore depends on permeability and the
pressure gradient — the quantity Darcy–Fick mode exists to model.

With $d$ the cell-center-to-face distance,
$\kappa_{\text{int}},\mu_{\text{int}}$ the adjacent interior cell's
properties, and the tag-aware sign (kernel methods in
`boundary_darcy_velocity`):

- **Dirichlet** ($P=P_{\text{bc}}$):
  $\;v_b = -\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,(\nabla P)$, with
  $(\nabla P)=(P_{\text{bc}}-P_{\text{int}})/d$ at the top and
  $(P_{\text{int}}-P_{\text{bc}})/d$ at the bottom (default
  $P_{\text{bc}}=101325$ Pa).
- **Neumann** ($\partial P/\partial n$ given along the **outward normal**):
  $\;v_b = -\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\partial P/\partial n$
  at the top, and $+\dfrac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\partial
  P/\partial n$ at the bottom ($\hat n=-\hat z$ there, so
  $\partial P/\partial z=-\partial P/\partial n$);
  zero gradient $\Rightarrow$ impermeable.
- **Convective (Robin)**: $\;v_b = \pm\,h_P\,(P_{\text{int}}-P_{\text{ext}})$,
  with the $+$ sign at the top boundary and $-$ at the bottom (so positive
  outflow is consistent with the $+z$ convention), and $h_P$ [m·Pa⁻¹·s⁻¹] the
  pressure-transfer coefficient — gated on $\kappa_{\text{int}}>0$ (an
  impermeable surface cell cannot vent regardless of $h_P$). As
  $h_P\to\infty$ this approaches Dirichlet; as $h_P\to0$ it approaches a
  zero-gradient (impermeable) Neumann condition.

These are detailed in §12.5; the shared kernel is `boundary_darcy_velocity`
(`src/Discretization/darcy_velocity.jl`), used identically by the residual, the
`BCMassOperator` Jacobian, and the saved-output MLR diagnostics. Each pressure
BC also re-applies the $\kappa>0,\ \mu>10^{-20}$ guards, and Dirichlet
additionally zeroes the velocity if $d<10^{-15}$ m.

### 9.5.4 Initializing the pore gas: background air

Darcy–Fick problems must **not** start with an evacuated pore space. With all
gas $\xi_j = 0$ the closure gives $P \approx 0$ inside the solid against a
101325 Pa boundary value, so early product gas is spuriously advected *inward*
down the absolute-pressure slope instead of being resisted by ~1 atm of pore
air (peer codes initialize $P = P_\infty$ with a background species or solve
gauge pressure). Two ingredients make the initial state physical:

1. **A real pore space.** The condensed components must declare an intrinsic
   (skeletal) density via `ρ_intrinsic` (or the `φ` porosity kwarg) on
   `SolidComponent`/`LiquidComponent`. Without it, `get_intrinsic_density`
   falls back to the bulk density — the legacy "no-pores" model — and a cell
   at full bulk density has $\varphi = 0$: there is no pore volume to hold
   air, and any generated gas hits the $\varphi$ floor ($10^{-12}$) in
   $P = N/\varphi$, producing absurd pressures.
2. **An inert background air component** initialized at ambient pressure:

   ```julia
   air  = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
   ξ_air = background_gas_concentration(material, :AIR, ξ_condensed₀, T₀)
   # → use ξ_air as the :AIR entry of the initial composition
   ```

   `background_gas_concentration` (src/Physics/pressure.jl) is the
   single-species inverse of the pressure closure,
   $\xi_j = P\,M_j\,\varphi(\xi,T)/(R\,T)$, so the initial cell pressure is
   exactly `P` (101325 Pa by default).

The shipped DarcyFick examples (`07_douglas_fir_cone.jl`,
`08_darcy_fick_demo.jl`, `08_pressure_treated_wood_cone.jl`) follow this
pattern. Caveat: the diffusive film BC drives every species toward a
zero-concentration ambient, so the background air slowly leaks out over an
internal-diffusion timescale (~minutes for cone-scale slabs); the
initialization matters most early in heating, before significant gas
generation. Advective air ingress from the ambient is not
modeled (boundary Darcy inflow upwinds against a zero-concentration ghost,
§9.5.3).

---

## 9.6 Pressure derivatives in the Jacobian

The structured-operator Jacobian (§15) requires the sensitivity of
$P=N(\xi,T)/\varphi(\xi,T)$ to the state. In the clamp-inactive interior
($0<\varphi<1$, $\varphi$ above its $10^{-12}$ floor) the quotient rule gives

$$
\frac{\partial P}{\partial T} = \frac{1}{\varphi}\,\frac{\partial N}{\partial T},\qquad
\frac{\partial P}{\partial \xi_k} = \frac{1}{\varphi}\frac{\partial N}{\partial \xi_k}
- \frac{N}{\varphi^{2}}\frac{\partial\varphi}{\partial\xi_k},
$$

with $\partial N/\partial\xi_k = R_gT/M_k$ for a gas and $0$ otherwise, and
$\partial\varphi/\partial\xi_k = -1/\rho_{i,k}$ for a condensed phase and $0$ for
a gas. Two physical readings follow: a **gas** component raises $P$ directly
through the partial-pressure sum; a **solid/liquid** component raises $P$
indirectly — adding condensed mass shrinks the void fraction, confining the same
gas to a smaller volume.

These derivatives are not hand-coded. The `DarcyTransportOperator`
(src/Jacobian/operators/darcy_transport.jl) seeds the stencil state as
`ForwardDiff.Dual`s and evaluates the *same* functions the residual calls —
`porosity` → `ideal_gas_pressure` → `effective_permeability` →
`effective_gas_viscosity` — extracting
$(P,\partial P/\partial T,\partial P/\partial\xi)$, $(\kappa,\dots)$,
$(\varphi,\dots)$ and $(\mu,d\mu/dT)$ from the Dual partials. Differentiating
through the executed code makes the Jacobian consistent with every guard by
construction: at the porosity floor, AD through `max(φ, 1e-12)` returns a zero
$\varphi$-partial, exactly matching the residual's clamped value (a hand-written
$N/\varphi^2$ would blow up as $10^{24}N$ there). The chain also picks up
temperature dependence the interior formulas above elide — intrinsic densities
$\rho_{i,j}(T)$ inside $\varphi$, and user-supplied $\kappa_j(T)$ — with no
separate code path to keep in sync.

---

## 9.7 The Darcy property cache and the velocity Jacobian

Darcy quantities are recomputed **once per residual evaluation** and stored in
the conditional `DarcyCache` (per-cell $P,\varphi,\kappa,\mu$; per-face
$v_g,\nabla P$), which is allocated only when Darcy flow is enabled (§5.13).
`update_darcy_properties!` (src/Discretization/finite_volume.jl) fills the
per-cell fields from the state ($\varphi$ → $P$ → $\kappa$ → $\mu$), forms each
interior face's mobility as the distance-weighted harmonic mean of the per-cell
$\kappa_i/\mu_i$ (exactly 0 if either side is $\le10^{-30}$; §9.5.2), writes
$v_g[f]=-\text{mob}_f\,(\nabla P)_f$, and finishes by applying the pressure BCs
to the boundary faces (§9.5.3).

Holding the mobility fixed, the velocity's pressure sensitivity is exact and
constant in $P$:

$$
\frac{\partial v_g[f]}{\partial P_L} = +\frac{\text{mob}_f}{d_{LR}},\qquad
\frac{\partial v_g[f]}{\partial P_R} = -\frac{\text{mob}_f}{d_{LR}}.
$$

The structured Jacobian does not consume a cached coefficient table for this:
the `DarcyTransportOperator` (§15.8.3) reconstructs $v_g$ from the stencil's
cell-local $(P,\kappa,\mu)$ inside its own closure (`_dt_mobility_face`,
src/Jacobian/operators/darcy_transport.jl — the same face rule as the residual)
and obtains both the direct $\partial v_g/\partial P$ and the full chain to
$(T,\xi)$ — through $P$, $\varphi$, $\kappa(T,\xi)$, and $\mu(T)$ — from the
forward-mode pass of §9.6, with an analytical override gated against the AD path
at `rtol = 1e-12`.

(`compute_darcy_mobility!`, `velocity_from_mobility!`,
`DarcyVelocityJacobianCoeffs`, and `compute_darcy_velocity_jacobian!` in
src/Discretization/darcy_velocity.jl are exported scalar utilities for
standalone analysis; none of them runs in the production residual or Jacobian.)

---

## 9.8 Surface diffusive mass flux

At the exposed surface gas escapes to the environment. The surface diffusive flux
is obtained from the internal-diffusion-only form $J=-\lambda\,\partial\xi/\partial z$
with the environment assumed to hold *zero* concentration of decomposition
products, so the boundary gradient is approximated as $(\xi-0)/d=\xi/d$ with $d$
the distance from the adjacent cell center to the boundary face.

### 9.8.1 Infinite vs. finite external resistance

- **Pure internal diffusion** (no external film resistance, $h_m\to\infty$):

  $$
  J_j^{\text{surface}} = \frac{\lambda\,\xi_j}{d}.
  $$

  The implementation carries the surface coefficient as $\lambda\rho$ and writes
  $J=\lambda\rho\,(\xi/\rho)/d=\lambda\xi/d$; the two forms are identical because
  $\rho$ cancels (consistent with the volume-fraction convention of §9.1.4). The
  effective $\lambda$ is recovered by exact division $\lambda=\lambda\rho/\rho$
  (the kernel returns zero flux for $\rho \le 0$; no floor is applied — the
  callers pass $\lambda\rho = \lambda\cdot\rho$, so the division is exact for any
  positive density).

- **Combined internal + external resistance** (finite film coefficient
  $h_m$ [m·s⁻¹]):

  $$
  \boxed{\;
  J_j^{\text{surface}} = \frac{\lambda\,h_m\,\xi_j}{\lambda + h_m\,d}
  \;}
  $$

  This is two resistances in series. Writing the flux as
  $J=\xi_j/(R_{\text{int}}+R_{\text{ext}})$ with internal resistance
  $R_{\text{int}}=d/\lambda$ and external resistance $R_{\text{ext}}=1/h_m$,

  $$
  J_j = \frac{\xi_j}{\dfrac{d}{\lambda}+\dfrac{1}{h_m}}
      = \frac{\lambda h_m \xi_j}{h_m d + \lambda},
  $$

  recovering the boxed form. As $h_m\to\infty$, $R_{\text{ext}}\to0$ and
  $J\to\lambda\xi_j/d$ (pure internal limit); as $h_m\to0$, $J\to h_m\xi_j$
  (external-controlled).

A purely **convective** mass BC (external transfer only) gives
$J_j^{\text{surface}}=h_m(\xi_j-\xi_{j,\infty})$, with the ambient concentration
$\xi_{j,\infty}$ typically zero for decomposition products venting to atmosphere.
Pressure and impermeable BCs return zero surface *diffusive* flux (a pressure BC
acts through the Darcy velocity, §9.5.3).

### 9.8.2 Two sign conventions, reconciled

The surface-flux helpers `_compute_surface_flux` (in `Physics/gas_transport.jl`)
return **positive flux for gas escaping** — the natural sense for ALE surface
regression (positive efflux $\Rightarrow$ the surface recedes). This is the
*opposite* of the general BC convention used by `mass_bc_flux`
(`Discretization/boundary_fluxes.jl`), which returns **negative for gas escaping**
because its convention is "positive = into the domain."

The two are reconciled when the flux is written into the face array. Because the
finite-volume divergence is $(J_R-J_L)A/V$ and the top face is the *right* face of
the surface cell, the BC flux is negated at the top boundary and applied directly
at the bottom:

$$
J[f_{\text{top}}] = -J_{\text{bc}}, \qquad J[f_{\text{bottom}}] = +J_{\text{bc}}.
$$

The sign flip at the top is **not** a bug; it is the consistent realization of the
signed divergence, so that gas escaping at $z=L$ produces a positive flux in the
$+z$ direction and a positive (outflow) contribution to the surface cell's species
balance. The full sign-convention summary is:

| Quantity | Convention |
|---|---|
| Interior face flux $J_j[f]$ | $+$ = transport in $+z$ |
| `_compute_surface_flux` (ALE) | $+$ = gas escaping |
| `mass_bc_flux` (BC) | $-$ = gas escaping |
| Stored top-face flux | $J[f]=-J_{\text{bc}}$ |
| Stored bottom-face flux | $J[f]=+J_{\text{bc}}$ |

*Implementation note.* `_compute_surface_flux` dispatch,
`src/Physics/gas_transport.jl` (lines 369–401, with the sign-convention block at
lines 350–363); `mass_bc_flux(::DiffusiveMassBC, …)`,
`src/Discretization/boundary_fluxes.jl` (lines 248–264); face application in
`apply_mass_bc_flux!`, `src/Discretization/finite_volume.jl`. See §12.4 for the
full mass-BC catalogue (Impermeable, MassFlux, Convective, Diffusive, Combined).

---

## 9.9 Species conservation and the gas-only divergence

The species mass-conservation PDE for a gaseous component $j$ is

$$
\frac{\partial \xi_j}{\partial t} = -\nabla\!\cdot J_j + S_{j}^{\text{rxn}},
$$

while a solid or liquid component has no transport flux and evolves by reaction
(and, in ALE, mesh dilation) only:

$$
\frac{\partial \xi_j}{\partial t} = S_{j}^{\text{rxn}}
\qquad (j \text{ solid/liquid}),
$$

(see §3.2 for the gas-vs-solid dispatch and §11.2 for the ALE corrections). The
reaction source $S_j^{\text{rxn}}=\sum_r \nu_{j,r}\,r_r$ is treated in §6.

### 9.9.1 Discrete finite-volume divergence

Integrating the flux divergence over cell $i$ and applying the divergence theorem
in 1D yields the face-difference operator

$$
\boxed{\;
\left.\frac{d\xi_j}{dt}\right|_{\text{transport},\,i}
   = -\frac{A}{V_i}\Big(J_j^{R} - J_j^{L}\Big)
\;}
$$

where $A$ is the cross-sectional area and $V_i$ the cell volume. The discrete
operator uses **one** area for both faces (`divergence_species_1d!` reads
`state.A_face[f_left]`): all faces share the global scalar $A_{\text{section}}$
(also under `WithChi`, where $A(t)$ is column-global), so this is exactly
conservative — the operator does not support per-face areas. The
multi-component form is assembled by `divergence_species_1d!` over all components,
but only gas components feed the species residual — the dispatch
`apply_species_residual_unrolled!` writes

$$
\frac{d\xi_{j,i}}{dt} = -\,\text{div}_J[j,i] + S_{j,i}^{\text{rxn}}
   \quad(j\text{ gas}),
\qquad
\frac{d\xi_{j,i}}{dt} = S_{j,i}^{\text{rxn}}
   \quad(j\text{ solid/liquid}).
$$

This gas-vs-solid branch is a compile-time-eliminated `is_gas` test (zero runtime
cost). Solid components have $J_j\equiv0$ at all faces, so even if their divergence
were formed it would vanish; the explicit dispatch makes the intent — and the
mass-thickness locking of solids to the mesh (§11) — unambiguous.

The gas flux that feeds this divergence is the diffusion-only flux (Eulerian
Fickian mode) or the Darcy–Fick flux (when permeability is present), selected by
the simulation-mode trait (§15.1). In ALE mode the same lab-frame flux feeds the
divergence; the mesh-motion frame change is the pointwise $+w\,\partial\xi/\partial z$
correction of §11.6 (the relative-velocity flux of §11.7 is a test-only
formulation, not the production path).

### 9.9.2 Coupling to the energy equation

Gas transport feeds back on temperature through two volumetric sources in the
Formulation-A energy equation (derived in detail in §7.4–7.5):

- **Always-present advective source** (the advection piece of the gas-enthalpy
  divergence):

  $$
  S_{\text{conv},i} = -\sum_{j:\,\text{gas}} c_{p,j}(T_i)\,\bar{J}_{j,i}\,
      \left(\frac{\partial T}{\partial z}\right)_i,
  \qquad \bar{J}_{j,i}=\tfrac12\big(J_j^{L}+J_j^{R}\big),
  $$

  with $(\partial T/\partial z)_i$ the opposite-distance-weighted gradient in the
  interior; boundary cells blend the interior one-sided slope with the
  boundary-face slope at the solved $T_s$ (§7.4.2).

- **Optional generation-enthalpy sink** (the divergence-completion piece, added
  only when `energy_form = :with_generation_sink`):

  $$
  S_{\text{gen},i} = -\sum_{j:\,\text{gas}} \Delta h_j(T_i,T_{\text{ref}})\,
      (\nabla\!\cdot J_j)_i,
  \qquad
  \Delta h_j(T_i,T_{\text{ref}}) = \int_{T_{\text{ref}}}^{T_i} c_{p,j}(\tau)\,d\tau,
  $$

  using the *true* sensible-enthalpy integral (midpoint quadrature via
  `delta_enthalpy`, §4.6) rather than a single-point product, so that the discrete
  energy ledger telescopes exactly with the surface transpiration BC, which uses
  the identical quadrature and $T_{\text{ref}}$.

These appear in the assembled temperature RHS
$\rho c_p^{\text{eff}}\,dT_i/dt = -(\nabla\!\cdot q)_i + S_{\text{conv},i}
+ S_{\text{gen},i} + Q_{\text{rad},i} + Q_{\text{rxn},i}$, in which the gas
sensible-heat *storage* term has been deliberately removed from
$\rho c_p^{\text{eff}}$ (the quasi-steady-gas approximation; §5.4, §7.3) and is
carried instead by $S_{\text{conv}}$. Crucially, $S_{\text{conv}}$ and the surface
transpiration BC are **mutually exclusive** — enabling both double-counts the
advected enthalpy by exactly $\int S_{\text{conv}}\,dz$, so the transpiration BC is
off by default and rejected at problem construction if combined (see §7.6, §12.3).

*Implementation note.* `divergence_species_1d!`,
`src/Discretization/finite_volume.jl` (lines 81–97);
`apply_species_residual_unrolled!`, `src/Materials/materials.jl` (lines 890–909);
`compute_gas_convective_source!` and `compute_gas_generation_sink!`,
`src/Physics/convection.jl` (lines 67–139, 178–220). The 8-step `compute_rhs!`
ordering (interior conductive fluxes $\to$ reactions $\to$ gas transport $\to$
thermal/mass BCs $\to$ radiation $\to$ $S_{\text{conv}}$ $\to$ divergences $\to$
$S_{\text{gen}}$ $\to$ final RHS) is described in §7.8.

---

## 9.10 Direct gas-flux Jacobian derivatives

The `GasTransportOperator` contributes the analytic Jacobian of the gas-transport
terms over a 3-cell tridiagonal stencil ($i-1,i,i+1$). The core building block is
the set of direct partial derivatives of the single-face diffusive flux

$$
J = -\lambda\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d},
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

with respect to each scalar argument. Writing $N_{\text{num}}=\xi_R T_R-\xi_L T_L$
and $\text{denom}=T_{\text{face}}\,d$:

$$
\frac{\partial J}{\partial \lambda} = -\frac{N_{\text{num}}}{\text{denom}},
\qquad
\frac{\partial J}{\partial \xi_L} = \frac{\lambda T_L}{\text{denom}},
\qquad
\frac{\partial J}{\partial \xi_R} = -\frac{\lambda T_R}{\text{denom}},
$$

$$
\frac{\partial J}{\partial T_L} = \frac{\lambda \xi_L}{\text{denom}} - \frac{J}{2\,T_{\text{face}}},
\qquad
\frac{\partial J}{\partial T_R} = -\frac{\lambda \xi_R}{\text{denom}} - \frac{J}{2\,T_{\text{face}}}.
$$

The temperature derivatives have two pieces: a direct effect through the numerator
($\pm\lambda\xi/\text{denom}$) and an indirect effect through the
$T_{\text{face}}$ prefactor in the denominator
($-J/(2T_{\text{face}})$, since $\partial T_{\text{face}}/\partial T_{L,R}=\tfrac12$
and $\partial(1/T_{\text{face}})/\partial T_{\text{face}}=-1/T_{\text{face}}^2$).
A guard sets $\partial J/\partial\lambda=0$ when $N_{\text{num}}=0$ exactly (a flat
$\xi T$ field has no diffusive flux and no $\lambda$-sensitivity).

Beyond these direct partials, the operator composes two additional chains:

- **Mixing-rule chain (via $\lambda_{\text{face}}$).** Since
  $\lambda_{\text{face}}$ is a harmonic mean of the cell coefficients, and each
  $\lambda_{L/R}$ is itself a composition-dependent mixing rule, the flux
  sensitivity to a concentration also includes

  $$
  \frac{\partial J}{\partial \xi_{L/R,k}}
  \;\supset\;
  \frac{\partial J}{\partial \lambda}\cdot
  \frac{\partial \lambda_{\text{face}}}{\partial \lambda_{L/R}}\cdot
  \frac{\partial \lambda_{L/R}}{\partial \xi_{L/R,k}},
  $$

  with the last two factors supplied by the harmonic-mean derivative helper and
  the effective-gas-transfer mixing derivatives (§5.6, §5.15).

- **Convective-source chain (via $S_{\text{conv}}$).** The energy-equation
  coupling contributes terms such as

  $$
  \frac{\partial S_{\text{conv}}}{\partial T_C} \;\supset\;
  \left(\frac{\partial c_p}{\partial T}\Big|_{T_C}\right)\bar{J}\,
  \left(\frac{\partial T}{\partial z}\right),
  $$

  i.e. the explicit $c_p$-chain through the temperature dependence of the gas heat
  capacity, plus the dependence of $\bar J$ and of the gradient on the stencil
  temperatures. The species rows likewise receive $-\partial(\text{div}_J)/\partial(\cdot)$
  via the divergence $\text{div}_J=(J_R A_R-J_L A_L)/V_C$, written into the species
  block as $\text{out}[1+j]=-\text{div}_J$ for gas components and $0$ for
  solids/liquids.

The analytic path (`_contribute_to_jacobian_analytical!`) is the default; an
AD-default path (`_contribute_to_jacobian_ad!`, ForwardDiff over the full
`apply_cell!` closure) is preserved for validation and is gated against the
analytic version at `rtol = 1e-12` in
`test/jacobian/test_operator_gas_transport_analytical.jl`. Boundary-face
contributions are zeroed in the operator — the mass BCs and pressure BCs are
handled by `BCMassOperator`, and the boundary cell's $S_{\rm conv}$ energy row
by `BCConvectionOperator` (§15.8.8) — keeping
the stencil pure and avoiding double-counting.

*Implementation note.* `_gt_J_direct_partials`,
`src/Jacobian/operators/gas_transport.jl` (lines 532–551);
`_gt_species_div_unrolled!` (lines 299–313) and the analytic operator
`_contribute_to_jacobian_analytical!` (from line 565). The Darcy-velocity pressure
derivatives that close the advective term's Jacobian are in §9.7. See §15.7–15.9
for how these operator contributions assemble into the structured Jacobian.

---

## 9.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The Pyrolysis.jl gas-transport model is a faithful descendant of the ThermaKin2Ds
and Gpyro condensed-phase formulations, with the FDS solid-phase model as the
minimal-transport reference. The notation map below is used on first contact;
internally all chapters use the unified symbols of §2.

| External | Pyrolysis.jl | Meaning |
|---|---|---|
| ThermaKin $\theta_i^{\,j}$ | $\nu_{i,j}$ | stoichiometric yield |
| Gpyro $K$ | $\kappa$ | permeability |
| Gpyro $\bar\rho$, $X_\alpha$ | $\rho$, $v_j$ | mixture density, volume fraction |
| Gpyro $\nu$ (Darcy) | $\mu/\rho_g$ | kinematic viscosity |
| FDS $Z$ | $A_i$ | Arrhenius pre-exponential |

### ThermaKin2Ds

ThermaKin writes the gas mass flux as a volume-fraction-gradient diffusion law,
$J = -\rho_g\,\lambda\,\partial(\xi_g/\rho_g)/\partial x$, which — for an ideal gas
with $\rho_g\propto P/T$ — is identical to the Pyrolysis.jl form
$J=-(\lambda/T)\,\partial(\xi T)/\partial x$ (§9.1; the header of
`src/Physics/gas_transport.jl` cites ThermaKin §1.2 explicitly). The two codes
therefore share the same diffusive driving potential and the same
composition + thermal split. Pyrolysis.jl additionally offers the Darcy advective
augmentation (§9.3) for permeable/pressure-driven problems; in the rigid,
non-swelling limit ($\gamma=0$ for gases) ThermaKin's Boyle's-law Darcy treatment
maps onto the porous-medium ideal-gas pressure (§9.4) with $\kappa$ controlling the
flow resistance. A genuine difference is the energy-equation convention:
Pyrolysis.jl (like FDS/Gpyro) *excludes* gas sensible-heat storage from
$\rho c_p^{\text{eff}}$ and carries it advectively in $S_{\text{conv}}$ (§9.9.2,
§5.4), whereas ThermaKin Eq. 17 sums $\sum_j \xi_j c_j$ over *all* components,
including gas. The cost of the exclusion is ~0.2% of the matrix term at 1 atm and
is exposed by the physical-form energy residual in the diagnostics (§16.2, §16.4).

### Gpyro

Gpyro solves pressure-driven gas transport with a Darcian mass flux
$\dot m'' = -(K/\nu)(\partial P/\partial z - g\rho_g)$ together with a
pressure-evolution equation derived from gas continuity and the ideal-gas law.
Pyrolysis.jl uses the same Darcy constitutive law,
$v_g=-(\kappa/\mu)\nabla P$ (§9.5.2), with the distance-weighted harmonic mean
of the per-cell mobilities $\kappa/\mu$ (Sutherland viscosity at the cell
temperature); the gravity/buoyancy term $g\rho_g$ is omitted (negligible
for thin condensed-phase samples in the through-thickness direction). The pressure
is closed algebraically by the porous-medium ideal-gas law $P=N/\varphi$ (§9.4)
rather than time-integrated as a separate state. This is **not** an additional
approximation relative to Gpyro's pressure ODE: gas storage $\partial\xi_g/\partial t$
is fully retained in the species ODEs, and $P=N/\varphi$ is an exact EOS
evaluation of that evolving state — summing the species equations reproduces gas
continuity with no spurious source. (The model's quasi-steady approximation lives
in the *energy* treatment — neglected gas heat storage, §1.4.1.) Gpyro's
control-volume conservation of $\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, and
gas mass corresponds to the cell-integrated species divergence of §9.9.1 on the
moving (ALE) mesh of §11; the harmonic interface interpolation of Gpyro's Patankar
scheme is the same resistance-in-series face averaging used here (§9.1.4, §13.11).

### FDS

The FDS solid-phase model treats the condensed phase with 1D Fourier conduction
and Arrhenius kinetics but, in its standard form, does **not** resolve in-depth
gaseous mass transport — volatiles are assumed to escape instantaneously to the
gas phase. The Pyrolysis.jl `NO_RADIATION`/Fickian configuration with negligible
permeability reduces toward this minimal-transport limit; the explicit
diffusive/Darcy machinery of this chapter is the generalization FDS lacks. The
energy-accounting convention (gas storage excluded, advective enthalpy carried
separately) is, however, aligned with FDS practice.

---

## 9.12 Limitations and open issues

- **First-order upwind advection only.** The Eulerian Darcy–Fick advective term
  uses unlimited first-order upwinding (§9.3.1). It is stable and conservative but
  diffusive; sharp composition fronts are smeared. Higher-order TVD limiting is
  implemented only for the ALE relative-velocity advection (§11.6–11.7) and is not
  exposed for the Eulerian gas flux.

- **No species discrimination in diffusion.** All gas species share one
  effective transfer coefficient $\lambda_{\text{eff}}$ (§9.2): light and heavy
  vapors diffuse identically within a simulation. This is a simplified
  surrogate for a full multicomponent (Stefan–Maxwell) diffusion model; the
  Chapman–Enskog Lennard-Jones parameter estimates available for sizing
  $\lambda$ (§5.10) are heuristic and unvalidated for pyrolysis vapors.

- **No Wilke mixing for viscosity.** The Darcy viscosity uses a single Sutherland
  value for an air-like gas (§9.5.2); a Wilke mixing rule would be required to
  capture mixtures of very different molecular weights (e.g. H₂O + CO₂ at high
  temperature; §5 open questions). The Sutherland reference temperature (273.15 K)
  differs from the enthalpy datum (298.15 K) and must not be conflated.

- **Zero-exterior-concentration surface assumption.** Surface diffusive fluxes
  assume the environment holds no decomposition products ($\xi_\infty=0$ implicit;
  §9.8.1). Re-entry or accumulation of products in the boundary layer is not
  modeled.

- **Quasi-steady gas energy.** Gas sensible-heat storage is excluded from the
  effective heat capacity and carried advectively by $S_{\text{conv}}$ — the
  quasi-steady-gas assumption (gas residence $\ll$ thermal diffusion time),
  whose validity degrades at high pressure or under extremely fast heating. The
  pressure itself is *not* approximated: $P=N/\varphi$ is an exact
  equation-of-state evaluation with gas mass storage fully retained in the
  species ODEs (§3.4.2). The diagnostics expose the energy cost but do not flag
  when it becomes significant (§16.9).

- **Negative-advective-flux robustness.** When upwinding selects a zero-concentration
  donor cell, the advective term can transiently produce a "wrong-side" flux
  (§9.3.1); physicality over long runs relies on reaction/heat coupling rather
  than an explicit positivity-preserving limiter.
