# 9. Gas-Phase Mass Transport and Pressure

This chapter derives the constitutive and discrete laws by which gaseous
decomposition products move through the condensed-phase matrix in Pyrolysis.jl.
Two transport mechanisms are treated: ordinary (Fickian) diffusion driven by a
volume-fraction / partial-pressure gradient, and Darcy advection driven by the
gas-phase pressure field in a porous medium; together they form the *Darcy–Fick*
flux. We establish the ideal-gas pressure closure in pore space, the Darcy
velocity, the boundary venting and film fluxes, the flux form of the species
rows, the sensible enthalpy the gas fluxes carry across faces, and the structure
of their Jacobian. Every equation in this chapter has been checked against the
implementing source and, where the textbook form and the code differ, the
chapter documents what the code actually does.

Throughout, the spatial coordinate is the through-thickness axis $z$, with
$z=0$ the back face (tag `:bottom`) and $z=L$ the exposed surface (tag `:top`).
A face flux is stored **per unit face area along the face normal**, which points
from the face's owner cell to its neighbour and out of the domain on a boundary
face: on a line every interior face and the `:top` face have their normal along
$+z$, and the `:bottom` face along $-z$, so a boundary-face flux is positive when
gas leaves through it, at either end (§2). The assembler scatters
$-F_fA_f$ into the owner's rows and $+F_fA_f$ into the neighbour's, so a positive
flux removes the quantity from the owner (§7.2, §13.3). In this chapter the index
$i$ denotes a **cell**, $j$ a **component**, and $f$ a **face**; $L,R$ denote the
two cells bounding an interior face, $L$ its owner and $R$ its neighbour (below
and above it on a line), and $d_L,d_R$ their center-to-face distances. The
integrated species state is the cell mass $m_{j,i}$ [kg]; every kernel reads the
mass concentration $\xi_{j,i}=m_{j,i}/V_i$ [kg·m⁻³], not the mass fraction
$Y_j=\xi_j/\rho$.

`setup` chooses the transport from the material and the mass conditions
(§9.7): none without a gas component; **Darcy–Fick** when some condensed
component is permeable; **Fickian** when gas can move at all (a component carries
a non-zero gas-transfer coefficient, or a mass condition is not `ImpermeableBC`);
none otherwise. One `GasTransportKernel` (`src/Kernels/gas_transport.jl`) then
writes every gas flux of the residual.

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
appear in a denominator. Its isobaric assumption is discussed in §9.11.

### 9.1.2 Discrete face flux

Discretizing $-\frac{\lambda}{T}\,\partial(\xi T)/\partial z$ at an interior face
$f$ between cells $L$ and $R$, with the temperature in the prefactor taken at the
arithmetic-mean face value $T_{\text{face}}=\tfrac12(T_L+T_R)$, gives the
implemented flux:

$$
\boxed{\;
J_j = -\,\lambda_{\text{face}}\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\;d_{LR}}
\;}
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

where $d_{LR}=d_L+d_R$ is the distance between the two cell centers along the
face normal and $\lambda_{\text{face}}$ is the distance-weighted harmonic mean of
the cell transfer coefficients (§9.1.4). This is exactly

```julia
@inline function gas_flux(λ_face::Real, ξ_L::Real, ξ_R::Real, T_L::Real, T_R::Real,
                          d_LR::Real)
    T_face = (T_L + T_R) / 2
    return -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)
end
```

The kernel evaluates it at every interior face whose two cells are live, for every
gas component, with three choices that are load-bearing:

- **Trial temperatures.** $T_L$, $T_R$ are the cells' trial temperatures
  $T_{\text{trial}}$ (the temperature smoothly clamped into the property bracket
  $[200, 3000]$ K, §5.13), so a trial state far outside the bracket still gives a
  finite, smooth flux.
- **Raw concentrations.** $\xi_L$, $\xi_R$ are the stored masses over the volume,
  $m/V$ (with the volume smoothly floored), not their positive parts: the flux is
  linear in the masses, and the transported quantity is exactly what the mass rows
  hold.
- **Floored distances.** $d_L$ and $d_R$ are each smoothly floored at $10^{-8}$ of
  their reference value (`DISTANCE_FLOOR_FRACTION`), which keeps a collapsing cell
  on a moving mesh from dividing by zero.

*Implementation note.* `gas_flux`, `src/Physics/gas_transport.jl`, called face by
face by `contribute!` of the `GasTransportKernel`, `src/Kernels/gas_transport.jl`.

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

This equivalence is verified to `atol = 1e-14` in
`test/physics/test_gas_closures.jl`.

### 9.1.4 Why the harmonic mean is applied to $\lambda$, not $\lambda\rho$

At a material interface, *flux continuity* — not property continuity — is the
physical constraint. Treating the two half-cells as transport resistances in
series, the face transfer coefficient is the distance-weighted harmonic mean

$$
\lambda_{\text{face}} = \frac{d_L + d_R}{\dfrac{d_L}{\lambda_L} + \dfrac{d_R}{\lambda_R}} ,
$$

which on a uniform mesh ($d_L=d_R$) reduces to the standard harmonic mean
$2\lambda_L\lambda_R/(\lambda_L+\lambda_R)$. The distance weighting makes the
scheme correct on stretched meshes (e.g. surface-refined grids); §7.1 gives the
identical treatment of thermal conductivity. The kernel forms it in the guarded
product form (`guarded_harmonic_mean`, `src/Properties/effective_properties.jl`;
§5.12)

$$
\boxed{\;
\lambda_{\text{face}} = \frac{(d_L + d_R)\,\lambda_L\lambda_R}{d_L\lambda_R + d_R\lambda_L + \eta\,(d_L+d_R)}
\;},\qquad \eta = 10^{-20}\ \mathrm{m^2\,s^{-1}}\ (\texttt{FACE\_MEAN\_GUARD}),
$$

which equals the resistance-in-series mean to a relative $\eta/\lambda_{\text{face}}$,
is smooth in every argument, and goes to zero continuously, with a finite
derivative, as either cell coefficient vanishes: a cell with $\lambda=0$ is an
infinite resistance and passes nothing, without a threshold that would cut the
flux off discontinuously.

A subtle but important modeling choice: the code takes the harmonic mean of
$\lambda$ **alone**, not of the product $\lambda\rho_g$ that would appear in a
true mass-fraction-gradient discretization $J=-\lambda\rho_g\,\partial Y/\partial z$.
This is consistent with — indeed *required by* — the volume-fraction formulation,
in which $\rho_g$ never appears explicitly: the density dependence is folded into
the $\xi T$ driving potential, which uses *cell-center* concentrations $\xi_L,\xi_R$
and the arithmetic-mean face temperature, with no face density. No
$\lambda\rho$-averaging helper exists in the package.

*Implementation note.* The per-cell $\lambda_i$ is evaluated once per cell and
residual in the state mapping, `ds.λ[i] = effective_gas_transfer(ξ⁺, T_trial, material)`
(`src/States/derived.jl`; §5.6, §5.13), at the positive-part concentrations and
the trial temperature. The face mean is formed inside the kernel, face by face, in
every evaluation; nothing is cached between evaluations.

### 9.1.5 Partial-pressure equivalence

Because the formulation is anchored to the ideal-gas relation, the same flux can
be written as a partial-pressure gradient. With $P_j = \xi_j R_g T / M_j$ the
partial pressure of species $j$,

$$
J_j = -\frac{\lambda_{\text{face}}\,M_j}{R_g\,T_{\text{face}}}\,\frac{\partial P_j}{\partial z}.
$$

To see the equivalence, note $\xi_j T = P_j M_j / R_g$, so
$\xi_R T_R - \xi_L T_L = (M_j/R_g)(P_{j,R}-P_{j,L})$, and substituting into the
discrete kernel gives the pressure form directly. This is verified to
`atol = 1e-12` in `test/physics/test_gas_closures.jl`. The partial-pressure
reading reinforces the physical picture: gas diffuses down its own
partial-pressure gradient, with mobility $\lambda M_j/(R_g T)$.

---

## 9.2 One transfer coefficient for all species

Every gas component diffuses with the **same** effective transfer coefficient
$\lambda_{\text{eff}}$: the face coefficient $\lambda_{\text{face}}$ entering
the kernel of §9.1.2 is a single mixture property (§5.6), shared by all species
at that face. The model therefore does not discriminate between light and heavy
gases within one simulation — a light product (H₂O) and a heavy one (MMA)
diffuse at the same rate through the same matrix. The per-component `λ_j(T)`
property controls how each *condensed or gas phase* contributes to the mixture
coefficient as composition evolves, not how individual gas species move
relative to one another.

Where a species-faithful magnitude is wanted, the Chapman–Enskog helper
`chapman_enskog_diffusivity(T, M; P, σ, ε_k)` (§5.10) estimates the binary
diffusivity of a given vapor in air; the user can fold that estimate (with a
porosity/tortuosity correction) into the component `λ` values. Genuine
multicomponent (Stefan–Maxwell) diffusion, or even a per-species kinetic-theory
scaling of the diffusive flux, is not implemented (§9.11).

---

## 9.3 The Darcy–Fick combined flux

When the material has a permeable condensed component, pressure-driven (Darcy)
advection augments diffusion. The total face flux is the sum of an advective and
a diffusive part:

$$
\boxed{\;
J_j^{\text{DF}} = \underbrace{\frac{\xi_j^{\text{upwind}}}{\varphi_{\text{face}}}\,v_g}_{J_j^{\text{adv}}}
   \;+\; \underbrace{\Big(\!-\lambda_{\text{face}}\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d_{LR}}\Big)}_{J_j^{\text{diff}}}
\;}
$$

where $v_g$ is the superficial Darcy velocity at the face along its normal
(§9.5, positive toward $R$, i.e. $+z$ on a line). The $1/\varphi_{\text{face}}$ is
forced by the bookkeeping: $\xi_j$ is mass per **bulk** volume (the pressure
closure $P=\sum\xi_jR_gT/(M_j\varphi)$ divides by $\varphi$ to recover the
intrinsic pore concentration) while the Darcy velocity is **superficial**
(volumetric flux per total cross-section), so the mass flux per bulk area is
intrinsic density × superficial velocity, $J^{\text{adv}}=(\xi/\varphi)\,v_g$ —
equivalently $\xi\,v_i$ with the interstitial velocity $v_i=v_g/\varphi$. The face
porosity is the arithmetic mean $\varphi_{\text{face}}=\tfrac12(\varphi_L+\varphi_R)$
of the two cell porosities, each already smoothly floored at $10^{-12}$ in the
state mapping (§9.4), so the divisor is positive. The kernel's two helpers:

```julia
@inline function _darcy_face(::DarcyFickTransport, ds, o, nb, d_o, d_n)
    m_f = guarded_harmonic_mean(ds.mob[o], ds.mob[nb], d_o, d_n, MOBILITY_GUARD)
    v = -m_f * (ds.P[nb] - ds.P[o]) / (d_o + d_n)
    return v, (ds.φ[o] + ds.φ[nb]) / 2
end

@inline _advective_flux(v, φ_f, ξ_o, ξ_n) = v * (value(v) >= 0 ? ξ_o : ξ_n) / φ_f
```

with `o` and `nb` the face's owner and neighbour ($L$ and $R$). Under Fickian
transport `_darcy_face` returns nothing and the advective part is zero.

### 9.3.1 First-order upwind selection

The advected concentration is taken from the **upwind** cell — the cell from which
the flow originates:

$$
\xi_j^{\text{upwind}} = \begin{cases}
\xi_{j,L}, & v_g \ge 0 \quad (\text{flow from } L \text{ to } R), \\
\xi_{j,R}, & v_g < 0 \quad (\text{flow from } R \text{ to } L).
\end{cases}
$$

This first-order upwind rule is stable and conservative (the same face flux
enters the two cells sharing the face with opposite signs). It is monotone: the
advective outflow of a cell is proportional to its own concentration and vanishes
with it, so advection alone cannot drive a gas inventory negative in the
semi-discrete equations (a time step that overshoots is rejected by the domain
test, §15.5.4). It is only first-order accurate; no slope reconstruction is applied
to the Darcy–Fick advective term. (Second-order reconstruction with a smooth van
Albada slope exists only for the fluxes relative to a moving mesh under a
correction, §11.6.2.)

The side is read from the **value** of $v_g$ alone (`value(v) >= 0`). Under the
colored Jacobian $v_g$ is a dual number, and a dual of zero value compares by its
partials, which differ between seeds; reading the value makes every derivative
pass take the derivative of the same side, the $v_g\ge0$ side at $v_g=0$. The flux
is continuous across $v_g=0$, with a kink there. The kernel tests check that the
side does not depend on the seed (`test/kernels/test_gas_transport.jl`).

### 9.3.2 Advection carries all species equally

The Darcy–Fick flux uses the bulk velocity $v_g$ in the advective term for every
component, and the diffusive term likewise uses one shared $\lambda_{\text{face}}$
for all species (§9.2). Physically, bulk pore flow is a mixture motion; the model
applies no molecular-weight discrimination in either term (§9.2, §9.11).

In the pure-Fickian limit $v_g\to0$, the Darcy–Fick flux reduces exactly to the
diffusion-only flux of §9.1; the kernel tests check an interior Darcy–Fick face
against the Fick flux plus the upwind advection of $\xi/\varphi$, and that no
advection crosses into an impermeable cell.

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

The state mapping evaluates the pressure once per cell and residual, under
Darcy–Fick transport only (`_gas_properties!`, `src/States/derived.jl`; §5.13):

$$
\varphi_i = \mathrm{smooth\_floor}\!\big(1-s(\xi^+_i,T_{\text{trial},i}),\,10^{-12}\big),
\qquad
P_i = \frac{N(\xi^+_i,T_{\text{trial},i})}{\max(\varphi_i,10^{-12})},
$$

with $s=\sum_{j\in\text{S/L}}\xi_j/\rho_{i,j}$ the skeletal volume fraction
(`skeletal_volume_fraction`). Three implementation details are load-bearing:

- **Only true gases contribute to $N$.** A component contributes iff
  `is_gas(j)` is true *and* `molar_mass(j) > 0`; solids and liquids are skipped.
  This is a compile-time-eliminated branch in the generated sum
  `_gas_partial_pressure_sum` (`ideal_gas_pressure`, `src/Physics/pressure.jl`).
- **The porosity floor is smooth.** `smooth_floor(x, a)` is $a$ for $x\le a$, $x$
  for $x\ge 2a$ and a cubic join between (§5.3.2), so $\varphi_i\ge10^{-12}$ with a
  continuous derivative, and the pressure stays finite even if a cell becomes
  fully solid; the `max` inside `ideal_gas_pressure` is then inactive. No upper
  clamp is needed: the skeletal fraction of non-negative concentrations is
  non-negative, so $\varphi_i\le1$.
- **Positive parts and trial temperatures.** The closure reads the positive-part
  concentrations $\xi^+$ and the trial temperature, so a trial state with a
  slightly negative mass or an out-of-bracket temperature still gives a finite,
  physical pressure.

The porosity is computed from *intrinsic* (skeletal) component densities
$\rho_{i,j}$, not bulk densities (§5.3). It is the same porosity the Carman–Kozeny
permeability closure (§9.5.2) and the face porosity of the advective flux (§9.3)
read, so the pressure, the permeability and the advection share one void fraction.
The domain test reads the *unguarded* pore fraction $1-s$ from the raw masses and
rejects a proposed step in which it leaves $[0,1]$ by more than the skeletal volume
of the mass tolerances (§15.5.4). `validate_problem` refuses a Darcy–Fick problem
whose initial pore fraction falls below $10^{-6}$ at a cell centre with a positive
effective permeability (`validate_darcy_initial_porosity`; §9.5.4).

The constant is $R_g = 8.31446261815324$ J·mol⁻¹·K⁻¹ (`R_GAS`).

---

## 9.5 Pressure gradient and Darcy velocity

### 9.5.1 Pressure gradient

At an interior face between two live cells the pressure gradient along the face
normal is a plain central difference of the cell pressures:

$$
(\nabla P)_f = \frac{P_R - P_L}{d_{LR}},\qquad d_{LR}=d_L+d_R ,
$$

with the smoothly floored distances of §9.1.2, formed inside the kernel at each
face (`_darcy_face`). $(\nabla P)_f>0$ when pressure increases toward $R$ ($+z$ on
a line). Boundary faces take their velocity from the pressure condition instead
(§9.5.3, §12.5).

### 9.5.2 Darcy velocity

Darcy's law gives the superficial (bulk) gas velocity at a face from the pressure
gradient, the permeability, and the gas viscosity:

$$
\boxed{\;
v_g = -\,\text{mob}_{\text{face}}\,(\nabla P)_f
    = -\,\text{mob}_{\text{face}}\,\frac{P_R - P_L}{d_{LR}},
\qquad \text{mob}=\frac{\kappa}{\mu}
\;}
$$

with $v_g>0$ for flow toward $R$ (toward the surface on a line). The face
properties combine the two cells consistently with their physical roles:

- **Mobility — guarded distance-weighted harmonic mean.** Permeability is a flow
  *resistance* in series across the interface. The state mapping forms the
  per-cell mobility $\text{mob}_i=\kappa_i/\mu_i$, and the kernel combines the two
  sides with the guarded harmonic mean of §9.1.4,

  $$
  \text{mob}_{\text{face}} = \frac{(d_L+d_R)\,\text{mob}_L\,\text{mob}_R}{d_L\,\text{mob}_R + d_R\,\text{mob}_L + \eta_m(d_L+d_R)},
  \qquad \eta_m = 10^{-30}\ \mathrm{m^2\,Pa^{-1}\,s^{-1}}\ (\texttt{MOBILITY\_GUARD}),
  $$

  far below the mobility of any permeable solid ($\kappa=10^{-18}$ m² and
  $\mu = 5\times10^{-5}$ Pa·s give $2\times10^{-14}$). It tends toward the
  *smaller* of the two mobilities, so a single tight (low-$\kappa$) cell throttles
  the interface flow, and it is **exactly 0** when either side's mobility is zero —
  an impermeable cell is an infinite series resistance, and gas cannot cross into
  it — while staying smooth, with a finite slope, as a mobility goes to zero.

- **Viscosity — per-cell trial temperature.** Viscosity is an intensive,
  smoothly varying gas property, evaluated per cell at the trial temperature and
  smoothly floored, $\mu_i=\mathrm{smooth\_floor}\big(\mu_g(T_{\text{trial},i}),\,10^{-12}\ \mathrm{Pa\,s}\big)$
  (`MU_FLOOR`); it enters the face through the per-cell mobilities. $\mu_g$ is the
  material's `gas_viscosity` (§5.9), one curve per material, independent of the
  pore-gas composition; the default `SutherlandProperty()` is Sutherland's law for
  air,

  $$
  \mu_g(T) = \mu_{\text{ref}}\left(\frac{T}{T_{\text{ref}}^{\mu}}\right)^{3/2}
             \frac{T_{\text{ref}}^{\mu}+S}{T+S},
  $$

  with the canonical air constants $\mu_{\text{ref}}=1.716\times10^{-5}$ Pa·s,
  $S=110.4$ K. **Caveat on the reference temperature:** the Sutherland reference
  is $T_{\text{ref}}^{\mu}=273.15$ K (the temperature at which $\mu_{\text{ref}}$ is
  quoted for air), which is *distinct* from the enthalpy datum
  $T_{\text{ref}}=298.15$ K used elsewhere in the model. The two should not be
  conflated.

- **Guards.** A cell's mobility is zero exactly when its permeability is zero; the
  viscosity floor keeps the division finite without a threshold.

*Implementation note.* `ds.mob[i] = effective_permeability(ξ⁺, T_trial, material, φ) / μ`
in `_gas_properties!` (`src/States/derived.jl`); the face mobility and velocity
in `_darcy_face` (`src/Kernels/gas_transport.jl`).

The effective permeability $\kappa$ of a cell is either a component mixing rule
(PARALLEL/SERIES/WEIGHTED on the renormalized solid basis, with the series bound
exactly 0 when any condensed phase present is impermeable; §5.7.1) or the
porosity-based **Carman–Kozeny** closure

$$
\kappa = \min\!\left(\frac{\ell^2\,\varphi^3}{180\,(1-\varphi)^2},\;\ell^2\right),
$$

with $\ell$ the characteristic pore/particle diameter and $180$ the Kozeny
constant for a packed bed of spheres; as $\varphi\to1$ the packed-bed form
diverges and $\kappa$ is capped at the free-flow scale $\ell^2$ (binding for
$\varphi\gtrsim0.93$), with a $10^{-30}$ denominator guard solely against
$\mathrm{Inf}$ at exactly $\varphi=1$. Carman–Kozeny is evaluated on the
smoothly floored porosity $\varphi_i$ the state mapping passes — the same porosity
the pressure formula uses — not on bulk-density volume fractions (§5.7). The two
closures (component mixing vs. Carman–Kozeny) are mutually exclusive, and
BRUGGEMAN is rejected for permeability (no physical basis for a transport
resistance; an `ArgumentError` is thrown). As virgin material converts to porous
char, $\varphi$ rises and $\kappa$ increases, easing gas escape.

### 9.5.3 Pressure boundary conditions and the Darcy venting flux

Pressure conditions do not enter the diffusive film flux; they set the outward
Darcy velocity of the boundary face, which drives the **advective venting flux
into the species rows**. Each tag's mass condition applies on the face that bounds
the live cells on its side (`Meshes.live_face`: the boundary face itself, or after
cells beyond it are deactivated the face of the outermost live cell; §14). With $i$
the boundary cell, $d$ its smoothly floored center-to-face distance (half the cell
on a line), and the face normal pointing **out of the domain at both ends**, the
stored face flux of gas component $j$ is

$$
F_{j,f} \;=\; \underbrace{\frac{v_{\text{out}}^{+}}{\varphi_i}\,\xi_{j,i}}_{\text{venting}}
\;-\;\underbrace{\texttt{mass\_bc\_flux}(bc,\,\xi_{j,i},\,\lambda_i\rho_i,\,\rho_i,\,d,\,t,\,j)}_{\text{film law, positive into the domain}},
\qquad
v_{\text{out}}^{+} = \begin{cases} v_{\text{out}}, & v_{\text{out}}\ge 0,\\ 0, & v_{\text{out}}<0,\end{cases}
$$

with $\rho_i$ the mixture density at $\xi^+_i$ and $T_{\text{trial},i}$
(`_boundary_mass_flux!`, `src/Kernels/gas_transport.jl`). Only outflow advects the
interior gas; inflow advects nothing, the ambient carrying none of the tracked
gases. The outflow test reads the value of $v_{\text{out}}$, as the interior
upwind side does (§9.3.1). Venting therefore depends on permeability and the
pressure difference — the quantity Darcy–Fick transport exists to model. Under
Fickian transport no velocity is formed and only the film law acts.

`outward_darcy_velocity(bc, P_i, mob_i, d)` gives the superficial velocity along
the outward normal, by condition:

- **Dirichlet** ($P=P_{\text{bc}}$, default 101325 Pa):
  $\;v_{\text{out}} = \text{mob}_i\,(P_i-P_{\text{bc}})/d$, the one-sided gradient:
  an overpressure vents.
- **Neumann** ($\partial P/\partial n$ given along the **outward normal**):
  $\;v_{\text{out}} = -\,\text{mob}_i\,\partial P/\partial n$: a pressure falling
  outward drives outflow; zero gradient $\Rightarrow$ impermeable.
- **Convective (Robin)**:
  $\;v_{\text{out}} = h_P\,(P_i-P_{\text{ext}})\,\dfrac{\text{mob}_i}{\text{mob}_i+\eta_m}$,
  with $h_P$ [m·Pa⁻¹·s⁻¹] the pressure-transfer coefficient and the gate one on a
  permeable cell (exactly, in floating point, for any physical mobility), zero on
  an impermeable one and smooth between: an impermeable surface cell cannot vent
  regardless of $h_P$. As $h_P\to\infty$ this approaches Dirichlet; as $h_P\to0$ it
  approaches a zero-gradient (impermeable) Neumann condition.
- Every other mass condition drives none; a `CombinedMassBC` sums its members.

Because both boundary normals point outward, the same formula holds at `:top` and
`:bottom` with no tag-dependent sign; the kernel tests check the closed forms and
that a pressure condition on `:bottom` vents the bottom cell along its outward
normal (`test/kernels/test_gas_transport.jl`). These conditions are detailed in
§12.5. The residual, its Jacobian (the derivative of the same residual) and the
saved-output diagnostics (the solution extraction re-evaluates the residual at
every saved state) all read these fluxes from the one kernel.

### 9.5.4 Initializing the pore gas: background air

Darcy–Fick problems must **not** start with an evacuated pore space. With all
gas $\xi_j = 0$ the closure gives $P \approx 0$ inside the solid against a
101325 Pa boundary value, so early product gas is spuriously advected *inward*
down the absolute-pressure slope instead of being resisted by ~1 atm of pore
air (peer codes initialize $P = P_\infty$ with a background species or solve
gauge pressure). Two ingredients make the initial state physical:

1. **A real pore space.** The condensed components must declare an intrinsic
   (skeletal) density via `ρ_intrinsic` (or the `φ` porosity keyword) on
   `SolidComponent`/`LiquidComponent`. Without it, `get_intrinsic_density`
   returns the bulk density — a phase with no pores of its own — and a cell
   at full bulk density has $\varphi = 0$: there is no pore volume to hold
   air, and any generated gas meets the $\varphi$ floor ($10^{-12}$) in
   $P = N/\varphi$, producing absurd pressures. `validate_problem` refuses such a
   problem when the cell is permeable (§9.4).
2. **An inert background air component** initialized at ambient pressure.
   `background_gas_concentration(material, species, ξ, T; P = 101325.0)`
   (`src/Physics/pressure.jl`) is the single-species inverse of the pressure
   closure, $\xi_j = P\,M_j\,\varphi(\xi,T)/(R_g\,T)$ with
   $\varphi=\max(1-s,0)$ evaluated from the condensed entries of $\xi$, so the
   initial cell pressure is exactly `P`:

   ```julia
   using Pyrolysis

   wood = SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.12, κ = 1e-14,
                         ρ_intrinsic = 1500.0)
   volatiles = GaseousComponent(:volatiles, M = 0.05, c = 1500.0, k = 0.03, λ = 2e-5)
   air = GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026, λ = 2e-5)
   material = Material(name = :wood, components = (wood, volatiles, air), reactions = ())

   ξ_air = background_gas_concentration(material, :AIR, (500.0, 0.0, 0.0), 300.0)
   ξ_initial = [z -> 500.0, z -> 0.0, z -> ξ_air]   # ξ_air = 0.78455 kg/m³
   ```

   Here the skeleton fills a third of the bulk ($500/1500$), so $\varphi = 2/3$,
   and the pressure of this composition at 300 K is 101325 Pa (checked with
   `ideal_gas_pressure` in the documentation scratch environment).

The shipped Darcy–Fick examples (`07_douglas_fir_cone.jl`,
`08_darcy_fick_demo.jl`, `08_pressure_treated_wood_cone.jl`) follow this
pattern. Caveat: the diffusive film condition (`DiffusiveMassBC`) drives every
gas toward a zero-concentration ambient, so the background air slowly leaks out
over an internal-diffusion time scale; the initialization matters most early in
heating, before significant gas generation. A `ConvectiveMassBC` with a
per-species ambient (`ξ_inf` an `NTuple`) can hold the air by film exchange
(§12.4.3). Advective air ingress from the ambient is not modeled (boundary Darcy
inflow carries nothing, §9.5.3).

---

## 9.6 The gas-transport Jacobian

The Jacobian of the gas-transport terms is part of the colored forward-mode
derivative of the whole residual (`ColoredJacobian`, §15.6): the state is seeded
with dual numbers, one partial per color, and the duals propagate through the
state mapping (the concentrations, their positive parts, the trial temperatures,
$\lambda_i$, $\varphi_i$, $P_i$, $\text{mob}_i$) and through the kernel (the face
means, the fluxes, the enthalpy they carry). No derivative of the gas transport is
hand-coded. Differentiating the executed code makes the Jacobian consistent with
every guard by construction: the smooth floors of the porosity, the viscosity and
the distances, the guarded face means, the positive-part concentrations and the
trial clamp all carry continuous derivatives, and the derivative is that of the
value the residual uses. This section states what the colored derivative
contains.

**Pressure.** In the interior of the guards ($\varphi$ above twice its floor,
masses above their positive-part scale, temperatures inside the bracket) the
quotient rule gives

$$
\frac{\partial P}{\partial T} = \frac{1}{\varphi}\,\frac{\partial N}{\partial T}
- \frac{N}{\varphi^{2}}\frac{\partial\varphi}{\partial T},\qquad
\frac{\partial P}{\partial \xi_k} = \frac{1}{\varphi}\frac{\partial N}{\partial \xi_k}
- \frac{N}{\varphi^{2}}\frac{\partial\varphi}{\partial\xi_k},
$$

with $\partial N/\partial T = N/T$, $\partial N/\partial\xi_k = R_gT/M_k$ for a gas
and $0$ otherwise, $\partial\varphi/\partial\xi_k = -1/\rho_{i,k}$ for a condensed
phase and $0$ for a gas, and $\partial\varphi/\partial T$ non-zero only for
temperature-dependent intrinsic densities $\rho_{i,j}(T)$. Two physical readings
follow: a **gas** component raises $P$ directly through the partial-pressure sum;
a **solid/liquid** component raises $P$ indirectly — adding condensed mass shrinks
the void fraction, confining the same gas to a smaller volume. At the porosity
floor the smooth floor's slope goes continuously to zero, so the
$\varphi$-partial vanishes there exactly as the value saturates (a hand-written
$N/\varphi^2$ would grow as $10^{24}N$). The state enters through
$\partial\xi^+_k/\partial m_{k,i}$ ($1/V_i$ above the positive-part scale) and
$\mathrm{d}T_{\text{trial}}/\mathrm{d}T$ (1 inside the bracket).

**Darcy velocity.** At fixed mobility the velocity's pressure sensitivity is exact
and constant in $P$,

$$
\frac{\partial v_g}{\partial P_L} = +\frac{\text{mob}_{\text{face}}}{d_{LR}},\qquad
\frac{\partial v_g}{\partial P_R} = -\frac{\text{mob}_{\text{face}}}{d_{LR}},
$$

and the mobility adds its own chain to $(T,\xi)$ of both cells, through
$\kappa(\xi^+,T,\varphi)$, $\mu(T)$ and the derivative of the guarded mean, which
stays finite as a mobility goes to zero.

**Diffusive flux.** The direct partial derivatives of the single-face flux

$$
J = -\lambda\,\frac{\xi_R T_R - \xi_L T_L}{T_{\text{face}}\,d},
\qquad T_{\text{face}}=\tfrac12(T_L+T_R),
$$

with respect to each scalar argument are, writing $N_{\text{num}}=\xi_R T_R-\xi_L T_L$
and $\text{denom}=T_{\text{face}}\,d$,

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
A flat $\xi T$ field has no diffusive flux and no $\lambda$-sensitivity. The
temperatures here are the trial temperatures, so each temperature partial is
multiplied by $\mathrm{d}T_{\text{trial}}/\mathrm{d}T$.

**Chains.** Beyond these direct partials, the colored derivative composes:

- **The mixing-rule chain (via $\lambda_{\text{face}}$).** Since
  $\lambda_{\text{face}}$ is a guarded mean of the cell coefficients, and each
  $\lambda_{L/R}$ is itself a composition-dependent mixing rule, the flux
  sensitivity to a concentration also includes

  $$
  \frac{\partial J}{\partial \xi_{L/R,k}}
  \;\supset\;
  \frac{\partial J}{\partial \lambda}\cdot
  \frac{\partial \lambda_{\text{face}}}{\partial \lambda_{L/R}}\cdot
  \frac{\partial \lambda_{L/R}}{\partial \xi_{L/R,k}}
  $$

  (§5.6, §5.15).
- **The advective chain.** $\partial J^{\text{adv}}/\partial\xi_{\text{up}} = v_g/\varphi_{\text{face}}$
  for the upwind cell, plus $\xi_{\text{up}}/\varphi_{\text{face}}$ times the
  velocity's chain and $-J^{\text{adv}}/\varphi_{\text{face}}$ times the face
  porosity's; at $v_g=0$ the derivative of the $v_g\ge0$ side (§9.3.1).
- **The enthalpy chain.** The face energy flux $\sum_j F_j\,h_j(T_f)$ (§9.9.2)
  differentiates to $\sum_j[\partial F_j\,h_j(T_f) + F_j\,c_{p,j}(T_f)\,\partial T_f]$,
  with $\partial T_f/\partial T_L = d_R/d_{LR}$ and $\partial T_f/\partial T_R = d_L/d_{LR}$,
  and the cell term $-\sum_j h_j(T_i)\dot M_{j,i}$ of the temperature row adds
  $-\sum_j c_{p,j}(T_i)\dot M_{j,i}$ to its own temperature column.
- **The boundary faces.** The film law contributes, for a `DiffusiveMassBC`,
  $\partial F_{j,f}/\partial\xi_{j,i} = \lambda_i h_m/(\lambda_i + h_m d)$ plus the
  chain through $\lambda_i$; the venting flux contributes $v_{\text{out}}/\varphi_i$
  and the chain of $v_{\text{out}}$ through $P_i$, $\text{mob}_i$ and $\varphi_i$
  while gas flows out, and nothing while $v_{\text{out}}<0$ (a kink at zero, read
  from the value). The enthalpy the boundary fluxes carry at $T_s$ adds the chain
  through the surface temperature of the thermal boundary kernel (§12.2.3).

The Jacobian's pattern is built from the mesh connectivity, not from the kernels
(§15.11): every cell row already reaches every unknown of its neighbours, so gas
transport adds no entries, and the coloring needs $3(N_C+1)$ colors on a fixed
mesh. The kernel tests (`test/kernels/test_gas_transport.jl`) compare the colored
Jacobian with the dense forward-mode reference (§15.13) entry by entry to a
relative $10^{-12}$ on five states (the initial state, a hot state, a perturbed
state, the zero state and a negative state), in both energy forms, on 6 and 8
cells, for Fickian and Darcy–Fick transport with film, convective, prescribed,
pressure and combined mass conditions and the Couette, property-ratio and
pooling film blowing models; they assert the
color count, that the pattern contains every non-zero of the reference, and that
the residual allocates nothing.

---

## 9.7 Transport modes and the cell fields

`setup` chooses one gas transport for the problem (`_gas_transport`,
`src/Integration/setup.jl`) and builds the kernel and the state mapping for it:

| Mode | Selected when | Cell fields of the state mapping | Face quantities in the kernel |
|---|---|---|---|
| none (`NoGasTransport`) | no gas component; or no component carries a non-zero $\lambda$ and every mass condition is `ImpermeableBC` | — | — (no kernel) |
| Fickian (`FickTransport`) | gas can move and no condensed component is permeable | $\lambda_i$ | $\lambda_{\text{face}}$, $J^{\text{diff}}$, the film fluxes, $T_f$ |
| Darcy–Fick (`DarcyFickTransport`) | some solid or liquid component has a permeability law not identically zero (`has_darcy_flow`) | $\lambda_i$, $\varphi_i$, $P_i$, $\text{mob}_i$ | the above, and $\text{mob}_{\text{face}}$, $v_g$, $\varphi_{\text{face}}$, $J^{\text{adv}}$, $v_{\text{out}}$ |

The cell fields are evaluated **once per cell and residual**, in the state mapping
(`derive!`, `src/States/derived.jl`; §5.13), before any kernel runs, at the
positive-part concentrations and the trial temperature:

$$
\lambda_i = \lambda_{\text{eff}}(\xi^+_i,T_{\text{trial},i}),\quad
\varphi_i = \mathrm{smooth\_floor}(1-s_i,10^{-12}),\quad
P_i = \frac{N_i}{\varphi_i},\quad
\text{mob}_i = \frac{\kappa(\xi^+_i,T_{\text{trial},i},\varphi_i)}{\mathrm{smooth\_floor}(\mu_g(T_{\text{trial},i}),10^{-12})}.
$$

The derived state holds the vectors for every problem and fills each only under
the transport that reads it: a Fickian problem evaluates no pressure, porosity or
mobility, and a problem without gas transport not even $\lambda_i$. The face
quantities are formed inside the kernel, face by face, in every evaluation;
nothing is cached between evaluations, and the dual twin of every buffer serves
the colored Jacobian (§15.2). The kernels read only the derived state (§15.3), so
the residual, its Jacobian and the saved-output diagnostics see one set of fields.

Under Fickian transport the cell coefficient is the volume-fraction mix of the
component values (§5.6), so a condensed component with $\lambda = 0$ stops all
transport in a cell that holds no gas. The gas field then ends in a sharp front;
`validate_problem` warns about the configuration (`warn_degenerate_gas_transport`:
"… the integrator rejects step after step on it. Give the condensed components a
small λ"), and a small $\lambda$ on the condensed components (for example
$10^{-6}$ m²·s⁻¹) keeps the transport non-degenerate.

---

## 9.8 Surface diffusive mass flux

At a boundary face gas escapes to the environment. The diffusive part of the
boundary flux is obtained from the internal-diffusion-only form
$J=-\lambda\,\partial\xi/\partial z$ with the environment assumed to hold *zero*
concentration of decomposition products, so the boundary gradient is approximated
as $(\xi-0)/d=\xi/d$, with $d$ the smoothly floored distance from the boundary
cell's center to the face (half the cell on a line) and $\xi$ the cell's raw
concentration $m/V$.

### 9.8.1 Infinite vs. finite external resistance

- **Pure internal diffusion** (no external film resistance, $h_m\to\infty$,
  `DiffusiveMassBC(Inf)`):

  $$
  J_j^{\text{surface}} = \frac{\lambda\,\xi_j}{d}.
  $$

  The kernel passes the surface coefficient as $\lambda\rho$, with $\lambda=\lambda_i$
  and $\rho$ the mixture density of the boundary cell at $\xi^+$ and $T_{\text{trial}}$,
  and the condition writes $J=\lambda\rho\,(\xi/\rho)/d=\lambda\xi/d$; the two forms
  are identical because $\rho$ cancels (consistent with the volume-fraction
  convention of §9.1.4). The condition recovers $\lambda$ by exact division
  $\lambda=\lambda\rho/\rho$ and returns zero flux for $\rho \le 0$; no floor is
  applied, so the division is exact for any positive density.

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

A purely **convective** mass condition (`ConvectiveMassBC`, external transfer only)
gives $J_j^{\text{surface}}=h_m(\xi_j-\xi_{j,\infty})$, with the ambient
concentration $\xi_{j,\infty}$ a scalar for every gas, a function of $t$, or a
per-species tuple (§12.4.3); it is typically zero for decomposition products
venting to atmosphere. A `MassFluxBC` prescribes the flux into the domain. Pressure
and impermeable conditions return zero *film* flux (a pressure condition acts
through the Darcy velocity, §9.5.3). The catalogue is §12.4.

### 9.8.2 Two sign conventions, reconciled

The condition functions share the convention of the boundary-condition layer,
**positive into the domain**: `mass_bc_flux(bc, ξ_surface, λρ, ρ, d, t, j)`
(`src/BoundaryConditions/boundary_fluxes.jl`) returns a **negative** flux for gas
escaping. The face fluxes of the residual follow the face-normal convention, and
every boundary normal points out of the domain, so the kernel stores the negated
film flux, `F[j, f] -= mass_bc_flux(...)`, at **both** boundaries, plus the
venting flux of §9.5.3. A boundary-face flux is therefore positive when gas
escapes, at `:top` and at `:bottom` alike, and it removes mass from the boundary
cell through the same scatter as any other face. The diagnostic
`compute_surface_mass_flux` returns $-$`mass_bc_flux`, positive for gas escaping.

| Quantity | Convention |
|---|---|
| Interior face flux $F_{j,f}$ | $+$ = transport along the face normal, owner → neighbour ($+z$ on a line) |
| Boundary face flux $F_{j,f}$ (`:top` and `:bottom`) | $+$ = gas leaving the domain (outward normal) |
| `mass_bc_flux` (condition) | $+$ = into the domain; $-$ = gas escaping |
| `compute_surface_mass_flux` | $+$ = gas escaping |
| `sol.boundary_mass_flux[tag]`, `sol.extras.surface_mass_flux` | $+$ = outward, summed over the components |

*Implementation note.* `mass_bc_flux` and `compute_surface_mass_flux`,
`src/BoundaryConditions/boundary_fluxes.jl`; the face application in
`_boundary_mass_flux!`, `src/Kernels/gas_transport.jl`. See §12.4 for the full
mass-condition catalogue (Impermeable, MassFlux, Convective, Diffusive, Combined)
and §16.6 for how the ledgers read the boundary fluxes.

---

## 9.9 Species conservation in flux form

The species mass-conservation PDE for a gaseous component $j$ is

$$
\frac{\partial \xi_j}{\partial t} = -\nabla\!\cdot J_j + S_j ,
$$

while a solid or liquid component has no transport flux and, on a fixed mesh,
evolves by reaction only:

$$
\frac{\partial \xi_j}{\partial t} = S_j
\qquad (j \text{ solid/liquid}),
$$

(see §3.2 for the gas-vs-solid distinction and §11.2.3 for the balances on a
moving cell). The reaction source $S_j=\sum_r \nu_{r,j}\,r_r$ is treated in §6.

### 9.9.1 The mass rows

The state holds the cell masses $m_{j,i}$, and the mass row of cell $i$ is the
cell-integrated balance (`assemble_rows!`, `src/Assembly/residual.jl`):

$$
\boxed{\;
\frac{\mathrm{d}m_{j,i}}{\mathrm{d}t} = \omega_{j,i}\,V_i + \dot M_{j,i},
\qquad
\dot M_{j,i} = -\sum_{f\in\partial i} s_{i,f}\,A_f\,F_{j,f}
\;}
$$

with $s_{i,f}=+1$ when cell $i$ owns face $f$ and $-1$ when it is the neighbour.
On a line, with $R$ the face above cell $i$ and $L$ the face below,
$\dot M_{j,i} = -(A_RF_{j,R} - A_LF_{j,L})$ on an interior cell, and a boundary
face, owned by its cell, removes $A_fF_{j,f}$ from it. Each face flux is
multiplied by its own face area, so the rows are exact for areas that differ
between faces (cylindrical and spherical coordinates); on a planar line every face
has the column's area $A$ (the variable area $A(t)$ under a lateral-shrinkage law,
§10). Every interior face flux leaves one cell and enters its neighbour exactly,
so the face fluxes telescope: the total mass of a component changes only by the
boundary fluxes and by reaction, and mass is conserved exactly because the masses
are the states (§16.3).

The gas-transport kernel writes face fluxes for the gas components only (it skips
every component that is not `is_gas`); the kernel tests check that the condensed
components' rows receive no gas-transport flux. The flux that feeds the rows is
the diffusion-only flux under Fickian transport or the Darcy–Fick flux under
Darcy–Fick transport (§9.7).

On a moving mesh (`use_ale = true`) the gas is mesh-bound: the mesh moves with the
solid, $w = v_s$, unless a correction is installed, and the gas flux of this
chapter is the flux relative to the skeleton, so it feeds the rows unchanged. The
volume change of a moving cell needs no term of its own: the rows integrate the
masses, and the concentrations $\xi=m/V$ follow the geometric volume, so the
dilation is implicit (§11.3.4). Under a correction ($w = v_s + c$), the
`RelativeFluxKernel` adds the flux $u\,\xi_{j,f}$, $u = -c\cdot n$, of **every**
component across the moving faces (§11.6).

### 9.9.2 Coupling to the energy equation: the gas enthalpy flux at faces

Gas transport feeds back on temperature through the sensible enthalpy its fluxes
carry across faces (§7.4). At every interior face the kernel adds to the face's
energy flux the enthalpy of the same mass fluxes the mass rows carry,

$$
F^{E}_f \mathrel{+}= \sum_{j:\,\text{gas}} F_{j,f}\,h_j(T_f),
\qquad T_f = \frac{d_R\,T_L + d_L\,T_R}{d_L + d_R},
$$

with $T_f$ the distance-weighted face temperature of the cell temperatures (not
the trial temperatures). At a boundary face the thermal boundary kernel adds
$\sum_j F_{j,f}\,h_j(T_s)$ to the outward energy flux once the surface balance has
fixed $T_s$ (`_carried_enthalpy`, `src/Kernels/thermal_boundary.jl`): gas leaves at
the surface temperature. The temperature row (§7.3) subtracts the enthalpy of the
cell's net mass change at its own temperature,

$$
C_i\,\frac{\mathrm{d}T_i}{\mathrm{d}t} = E_i - \sum_j h_j(T_i)\,\dot M_{j,i},
$$

so each face contributes, net,

$$
-\,A_f\,s_{i,f}\sum_{j:\,\text{gas}} F_{j,f}\,\big(h_j(T_f)-h_j(T_i)\big)
$$

to cell $i$: the flux form of $-\sum_j \mathbf J_j\cdot\nabla h_j$. For a cell with
equal fluxes $F$ through both faces this is
$-A F\,[h(T_{R,f})-h(T_{L,f})] \approx -V\,F\,c_p\,\partial T/\partial z$, the
advective term of the energy equation. The form is conservative (each face carries
one enthalpy flux, which leaves one cell and enters the other), exact for
temperature-dependent $c_{p,j}$, and keeps a uniform temperature exactly at rest:
with $T_f = T_i$ everywhere every face term vanishes, whatever the mass fluxes. The
kernel tests check this to $10^{-12}$ in both energy forms, with Darcy flow
through a non-uniform gas (`test/kernels/test_gas_transport.jl`). In the enthalpy
form the same face fluxes enter the inventory row, which needs no cell term: the
inventory $H_i=\sum_j m_{j,i}h_j(T_i)$ changes with the masses by itself (§7.3).

The storage coefficient $C_i=\sum_j m_{j,i}c_{p,j}$ sums every component, the gas
included (§5.4); the gas *flux* enters only through the enthalpy it carries across
faces. The surface energy balance therefore carries no gas-enthalpy term: one on
top of the boundary-face enthalpy flux would count the escaping enthalpy a second
time (§7.4.3, §12.2). The blowing models of a convective condition read the
outward gas-transport fluxes of the face (§12.3.3).

*Implementation note.* The mass and enthalpy face fluxes in `contribute!` of the
`GasTransportKernel`; the scatter and the rows in `assemble_rows!`,
`src/Assembly/residual.jl`. The order of the kernels (reactions, gas transport,
in-depth radiation, conduction, thermal boundary, relative fluxes) is described in
§7.6; gas transport runs before the thermal boundary kernel so that the surface
balance and the carried enthalpy see its boundary fluxes.

---

## 9.10 Comparison to ThermaKin2Ds, Gpyro, and FDS

The Pyrolysis.jl gas-transport model follows the ThermaKin2Ds and Gpyro
condensed-phase formulations, with the FDS solid-phase model as the
minimal-transport reference. The notation map below is used on first contact;
internally all chapters use the unified symbols of §2.

| External | Pyrolysis.jl | Meaning |
|---|---|---|
| ThermaKin $\theta_i^{\,j}$ | $\nu_{i,j}$ | stoichiometric yield |
| Gpyro $K$ | $\kappa$ | permeability |
| Gpyro $\bar\rho$, $X_\alpha$ | $\rho$, $v_j$ | mixture density, volume fraction |
| Gpyro $\nu$ (Darcy) | $\mu/\rho_g$ | kinematic viscosity |
| Gpyro, FDS $Z$ | $A_i$ | Arrhenius pre-exponential |

### ThermaKin2Ds

ThermaKin writes the gas mass flux as a volume-fraction-gradient diffusion law,
$J = -\rho_g\,\lambda\,\partial(\xi_g/\rho_g)/\partial x$, which — for an ideal gas
with $\rho_g\propto P/T$ — is identical to the Pyrolysis.jl form
$J=-(\lambda/T)\,\partial(\xi T)/\partial x$ (§9.1; the header of
`src/Physics/gas_transport.jl` cites ThermaKin §1.2). The two codes therefore
share the same diffusive driving potential and the same composition + thermal
split. Pyrolysis.jl additionally offers the Darcy advective augmentation (§9.3) for
permeable, pressure-driven problems; in the rigid, non-swelling limit ($\gamma=0$
for gases) ThermaKin's Boyle's-law Darcy treatment maps onto the porous-medium
ideal-gas pressure (§9.4) with $\kappa$ controlling the flow resistance. The
energy-equation convention is the same as ThermaKin's: its Eq. 17 sums
$\sum_j \xi_j c_j$ over *all* components including gas and carries the gas flux
through $c_g(\mathbf J_g\cdot\nabla T)$, which is $\rho c_p^{\text{eff}}$ (§5.4)
and, discretely, the enthalpy the gas fluxes carry across faces (§9.9.2), a
conservative form of the same term. ThermaKin's Eq. 15 applies its deformation term
to every component, gas included — the skeleton-bound frame of §11.1.2.

### Gpyro

Gpyro solves pressure-driven gas transport with a Darcian mass flux
$\dot m'' = -(K/\nu)(\partial P/\partial z - g\rho_g)$ together with a
pressure-evolution equation derived from gas continuity and the ideal-gas law.
Pyrolysis.jl uses the same Darcy constitutive law, $v_g=-(\kappa/\mu)\nabla P$
(§9.5.2), with the guarded distance-weighted harmonic mean of the per-cell
mobilities $\kappa/\mu$ (the material's viscosity curve at the cell temperature);
the gravity/buoyancy term $g\rho_g$ is omitted (negligible for thin
condensed-phase samples in the through-thickness direction). The pressure is
closed algebraically by the porous-medium ideal-gas law $P=N/\varphi$ (§9.4)
rather than time-integrated as a separate state. This is **not** an additional
approximation relative to Gpyro's pressure equation: gas storage is fully retained
in the mass rows of the gas components, and $P=N/\varphi$ is an exact evaluation of
the equation of state at that evolving state — summing the species rows reproduces
gas continuity with no spurious source. Gpyro's control-volume conservation of
$\bar\rho\,\Delta z$, $\bar\rho Y_i\,\Delta z$, and gas mass — the gas balance
written on control volumes that deform with the solid, with only $\dot m''$
crossing their faces — corresponds to the cell mass rows of §9.9.1 on the moving
mesh of §11, where every component rides the mesh; the harmonic interface
interpolation of Gpyro's Patankar scheme is the same resistance-in-series face
averaging used here (§9.1.4, §13.11).

### FDS

The FDS solid-phase model treats the condensed phase with 1D Fourier conduction
and Arrhenius kinetics but, in its standard form, does **not** resolve in-depth
gaseous mass transport — volatiles are assumed to escape instantaneously to the
gas phase. A Pyrolysis.jl problem without gas transport (no gas component, or no
transfer coefficient and impermeable mass conditions) is this minimal-transport
limit; the diffusive and Darcy machinery of this chapter is the generalization
FDS lacks, and with no in-depth gas FDS has neither gas storage nor the enthalpy
the gas carries through the solid.

---

## 9.11 Limitations and open issues

- **First-order upwind advection only.** The Darcy–Fick advective term uses
  unlimited first-order upwinding (§9.3.1). It is stable, conservative and
  monotone but diffusive; sharp composition fronts are smeared. Second-order
  reconstruction with a smooth van Albada slope is implemented only for the
  fluxes relative to a moving mesh (§11.6.2) and is not used for the Darcy flux.
  The upwind face state and the outflow-only vent switch on the value of the
  velocity, so the flux has a kink at $v_g = 0$, where the Jacobian takes the
  $v_g\ge0$ side.

- **No species discrimination in diffusion.** All gas species share one
  effective transfer coefficient $\lambda_{\text{eff}}$ (§9.2): light and heavy
  vapors diffuse identically within a simulation. This is a simplified
  surrogate for a full multicomponent (Stefan–Maxwell) diffusion model; the
  Chapman–Enskog Lennard-Jones parameter estimates available for sizing
  $\lambda$ (§5.10) are heuristic and unvalidated for pyrolysis vapors.

- **Isobaric diffusive kernel.** $-(\lambda/T)\,\partial(\xi T)/\partial z$ equals
  Fick's law $-\rho_g D\,\nabla Y$ only when $\nabla(\varphi P)=0$; the dropped term
  is proportional to $-Y\,\nabla(\varphi P)/(\varphi P)$. Under Fickian transport,
  which resolves no pressure field, the convention is consistent. Under Darcy–Fick
  transport the code resolves $\nabla P$ (and sharp char fronts carry
  $\nabla\varphi$), so the diffusive term is then an approximation: the
  pressure-driven part of the transport is carried by the Darcy term while the
  diffusive term keeps the isobaric form (`gas_flux` docstring).

- **No Wilke mixing for viscosity.** The Darcy viscosity is one curve per material
  (Sutherland's law for air by default; §5.9, §9.5.2), independent of the pore-gas
  composition; a Wilke mixing rule would be required to capture mixtures of very
  different molecular weights (e.g. H₂O + CO₂ at high temperature; §5.16). The
  Sutherland reference temperature (273.15 K) differs from the enthalpy datum
  (298.15 K) and must not be conflated.

- **Zero-exterior-concentration film law.** The diffusive film condition assumes
  the environment holds no decomposition products ($\xi_\infty=0$ implicit;
  §9.8.1); only `ConvectiveMassBC` carries an ambient concentration. Re-entry or
  accumulation of products in the boundary layer is not modeled, and boundary
  Darcy inflow carries no ambient gas (§9.5.3).

- **No gravity in Darcy's law.** The buoyancy term $g\rho_g$ is omitted (§9.10).

- **One-temperature gas energy.** The pore gas is stored at the matrix
  temperature and its fluxes carry enthalpy at the face temperature (§9.9.2);
  there is no two-temperature gas energy equation, so regimes where the gas does
  not equilibrate with the pore walls within a cell are outside the model. The
  pressure itself is *not* approximated: $P=N/\varphi$ is an exact
  equation-of-state evaluation with gas mass storage fully retained in the
  mass rows (§3.4.2).
