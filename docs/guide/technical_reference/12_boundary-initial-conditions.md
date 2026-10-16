# 12. Boundary and Initial Conditions

This chapter specifies how Pyrolysis.jl closes the coupled energy, species, and
gas-pressure equations at the two physical surfaces of the 1D column ($z=0$
substrate / `:bottom`, $z=L$ exposed surface / `:top`) and how it seeds the
initial state. It derives the surface energy balance and the Newton iteration
that solves for the true surface temperature $T_s$ in a cell-centered
finite-volume scheme; documents every concrete thermal, mass, pressure, and
mesh-motion boundary condition with its exact flux law, temperature derivative,
and sign convention as implemented; explains the type-stable
`BoundaryConditionSet` container and the tag$\to$id bridge; and treats the
initial-condition closures and their auto-defaults. Throughout, the equations are
checked term-by-term against the source so that the math and the code agree
exactly — including the boundary-face sign flips, the Kirchhoff emissivity
coupling, the neglected $d\varepsilon/dT$ term, and the absence of a surface
gas-enthalpy term alongside the volumetric advection source $S_{\text{conv}}$.

This chapter conforms to the unified nomenclature of §2, Nomenclature, and uses
the three binding sign conventions: positive flux is transport in the $+z$
direction (bottom$\to$top); the discrete divergence $(F_R-F_L)A/V$ is positive
when the quantity leaves a cell; and for heats of reaction $h>0$ is endothermic
with $Q_{\text{rxn}}=-\sum_r h_r r_r$. In this chapter the index $i$ denotes the
cell index and $j$ the component index, consistent with the discretization
chapters.

---

## 12.1 Overview and design philosophy

### 12.1.1 What a boundary condition must supply

The continuum model of §3, Governing Equations, is a system of conservation laws
for the condensed-phase temperature $T(z,t)$ and the per-component mass
concentrations $\xi_j(z,t)$ (kg·m⁻³), with the gas-phase pressure $P$ as an
algebraic closure in Darcy mode. In flux (divergence) form the spatial operators
are

$$
\rho c_p^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\mathbf{q}_{\text{cond}} + S_{\text{conv}}
   + Q_{\text{rad}} + Q_{\text{rxn}},
\qquad
\frac{\partial \xi_j}{\partial t}
= -\nabla\!\cdot\mathbf{J}_j + S_{j,\text{rxn}} .
$$

A cell-centered finite-volume discretization (see §13, Finite-Volume
Discretization) integrates each equation over a cell $i$ and replaces the
divergence by a difference of face fluxes,

$$
(\nabla\!\cdot \mathbf{F})_i = \frac{(F_{f_R}-F_{f_L})\,A}{V_i}.
$$

For an interior face the flux $F_f$ is a constitutive law evaluated from the two
adjacent cells (Fourier, Fick, Darcy–Fick). For a **boundary face** there is no
neighboring cell, so the flux must instead be supplied by a boundary condition.
The boundary-condition subsystem therefore has exactly one job: for each boundary
face, given the interior state and the time $t$, produce the heat flux
$q_{\text{BC}}$, the species fluxes $J_j$, and (in Darcy mode) the gas velocity
$v_g$ that enter the divergence operator with the correct sign.

There is one subtlety unique to the **thermal** boundary: because temperature is
stored at cell centers and radiative and convective fluxes depend nonlinearly on
the *surface* temperature $T_s$ (which is not a state variable), the code does not
simply evaluate $q_{\text{BC}}$ at the adjacent cell-center temperature. Doing so
produces an $\sim$20% energy imbalance for $T^4$ radiation. Instead it solves a
small surface energy balance for the true $T_s$ (§12.2). The species and pressure
boundaries have no such nonlinearity and are evaluated directly.

### 12.1.2 Type hierarchy and dispatch

Boundary conditions are organized by physical kind, each rooted in an abstract
supertype declared in the `Materials` module:

- `AbstractThermalBC` — energy/heat boundary conditions;
- `AbstractMassBC` — species and pressure-driven (Darcy) boundary conditions;
- `AbstractRadiationBC` — P1 radiative-transfer boundary conditions (in the
  `Experimental` module; see §8, Thermal Radiation).

Concrete types live in the `BoundaryConditions` module
(`src/BoundaryConditions/boundary_conditions.jl`); the thermal flux laws with their
$T_s$ derivatives (`evaluate_bc_flux_and_derivative`) and the surface-temperature
Newton solve live in `src/Discretization/boundary_temperature.jl`, and the mass-BC
flux laws (`mass_bc_flux`) and the diagnostic heat-flux decomposition
(`compute_boundary_heat_fluxes`) in `src/Discretization/boundary_fluxes.jl`. The whole set is stored in
a type-stable, tag-keyed `BoundaryConditionSet` (§12.7) so that the residual
specializes statically with zero dispatch overhead.

> **Implementation note.** Concrete BC types and the `BoundaryConditionSet`
> container are defined in `src/BoundaryConditions/boundary_conditions.jl`; flux
> laws in `src/Discretization/boundary_fluxes.jl` (mass BCs) and, for thermal BCs
> with their $T_s$ derivatives, `src/Discretization/boundary_temperature.jl`; BC application into the flux
> arrays in `src/Discretization/finite_volume.jl` (`apply_thermal_bc_flux!`,
> `apply_mass_bc_flux!`); pressure BCs in `src/Discretization/darcy_velocity.jl`
> (`apply_pressure_bc_to_darcy_velocity!`).

### 12.1.3 The boundary-face sign convention

Two boundaries exist in 1D, with the **binding** id/tag/position mapping
(`src/Geometry/topology.jl`, `src/Geometry/mesh_generation.jl`):

| Tag | `boundary_id` | Position | Face id |
|---|---|---|---|
| `:top` | 1 | $z=L$, exposed surface | last node ($n+1$) |
| `:bottom` | 2 | $z=0$, substrate | first node ($1$) |

`tag_to_boundary_id(:top)=1`, `tag_to_boundary_id(:bottom)=2`,
`boundary_id_to_tag(1)=:top`.

Because the discrete divergence is $(F_{f_R}-F_{f_L})A/V_i$ with the convention
that a positive face flux points in $+z$, the BC application must orient the
computed *into-the-domain* flux correctly for each boundary. A flux
$q_{\text{cond}}$ that carries heat **into** the slab is in the $+z$ direction at
the bottom face (left face of cell 1) but in the $-z$ direction at the top face
(right face of cell $n$). Hence the implementation stores

$$
q[f_{\text{bottom}}] = +q_{\text{conducted}}, \qquad
q[f_{\text{top}}] = -q_{\text{conducted}},
$$

and identically for mass fluxes (`apply_thermal_bc_flux!`, `apply_mass_bc_flux!`,
`finite_volume.jl`). This is not a bug; it is the consistent realization of the signed
divergence. We return to it in §12.2.5 and §12.4.4.

---

## 12.2 The surface energy balance and the Newton solve

The surface temperature $T_s$ is the keystone of the thermal boundary. It is the
temperature of the exposed face itself, distinct from the adjacent cell-center
temperature $T_{\text{cell}}$, and it is determined by an instantaneous balance
between the externally applied flux and the heat conducted into the slab.

### 12.2.1 Continuous statement

Consider the exposed face. An energy balance on an infinitesimally thin control
volume at the surface (zero thermal mass) equates the net flux applied by the
environment to the flux conducted into the material:

$$
q_{\text{BC}}(T_s) \;=\; q_{\text{cond}}(T_s, T_{\text{cell}}),
$$

with, for the most general thermal BC,

$$
q_{\text{BC}}(T_s) = q_{\text{ext}}\;-\;\varepsilon_{\text{eff}}\sigma T_s^4
\;-\;h_{\text{conv}}\,(T_s-T_\infty)
$$

(absorbed external flux, re-radiation loss, convective loss; all signed positive
into the domain) and the conducted flux written as a one-sided Fourier difference
from the surface to the first cell center,

$$
q_{\text{cond}}(T_s,T_{\text{cell}})
= k_{\text{eff}}\,\frac{T_s - T_{\text{cell}}}{\Delta z_{1/2}},
$$

where $\Delta z_{1/2}$ is the distance from the boundary face to the first cell
center and $k_{\text{eff}}$ is the effective conductivity of the boundary cell
(§5, Effective Properties & Mixing Rules). Both sides are positive when heat
flows into the material. The gas enthalpy carried by the flux relative to the
skeleton is transported by the volumetric source $S_{\text{conv}}$ (§7), not by
a surface term; under local thermal equilibrium the gas crosses the surface at
$T_s$ on both sides, so the balance has no gas-enthalpy term (§1.4.3, §7.4.3). The
outward gas flux enters the balance only through the blowing model of a
`ConvectiveBC`, which scales the convective coefficient (§12.3.3).

The conductive distance is exactly half the boundary cell's thickness, computed
from state-dependent geometry (AD- and ALE-clean) as

$$
\Delta z_{1/2}
= \frac{1}{2}\,\frac{V_{\text{cell}}}{A_{\text{face}}}
\qquad\text{(`boundary\_temperature.jl`)}.
$$

Using the half-width (rather than the full cell width) places the surface at the
physical face and the datum $T_{\text{cell}}$ at the cell center, which is the
correct one-sided gradient for a cell-centered scheme.

### 12.2.2 Residual and Newton iteration

Define the surface residual

$$
f(T_s) \;=\; q_{\text{BC}}(T_s)\;-\;q_{\text{cond}}(T_s).
$$

The root $f(T_s)=0$ is found by Newton iteration

$$
T_s^{(k+1)} = T_s^{(k)} - \frac{f(T_s^{(k)})}{f'(T_s^{(k)})},
\qquad
f'(T_s) = \frac{dq_{\text{BC}}}{dT_s} - \frac{k_{\text{eff}}}{\Delta z_{1/2}}
$$

(`solve_boundary_temperature!`, `boundary_temperature.jl`).
The initial guess is the adjacent cell temperature clamped to the property
bracket, $T_s^{(0)} = \min\bigl(\max(T_{\text{cell}}, 200\ \text{K}),
3000\ \text{K}\bigr)$, so the solve is a pure function of its arguments. The
per-boundary `BoundaryState` receives the converged $T_s$ and fluxes for
diagnostics and is never read back; a non-finite final iterate is not stored.
A residual evaluation therefore never depends on the ones before it: a trial
state of a diverging outer Newton, with a cell temperature far outside the
bracket, fails that trial and leaves every later evaluation untouched. The
conduction derivative
$dq_{\text{cond}}/dT_s = k_{\text{eff}}/\Delta z_{1/2}$ is exact and always
present; the BC derivatives are given per-type in §12.2.4.

**Convergence.** The residual is compared against a relative tolerance scaled by
the dominant physical flux,

$$
|f(T_s)| < \max\!\Bigl(\texttt{tol}\cdot
\max\!\bigl(|q_{\text{BC}}|,\,|q_{\text{cond}}|,\,1\bigr),\ f_{\text{floor}}\Bigr),
\qquad \texttt{tol}=10^{-10},
$$

with a maximum of 50 iterations. The floor
$f_{\text{floor}} = 32\,\epsilon_{\text{mach}}\max(|T_s|,1)\,\bigl(k_{\text{eff}}/\Delta z_{1/2}+|dq_{\text{BC}}/dT_s|\bigr)$
is the flux change from a few ulps of $T_s$: the conducted flux is a difference of
two temperatures, so a face with near-zero net flux cannot resolve its residual
below it. Scaling by
the dominant flux makes the criterion meaningful whether the surface flux is
$\sim 10^4$ W·m⁻² (cone) or near zero (near-adiabatic). **Nothing aborts on
non-convergence**: on max-iterations or a vanishing Jacobian
($|f'|<10^{-30}$ — e.g. an adiabatic boundary with no temperature dependence
anywhere) the loop emits a throttled `@warn` and falls through to the final
flux evaluation at the unconverged iterate; the ODE integrator's step-error
control is the backstop that rejects a step built on a bad surface state.

**Clamping.** Each Newton iterate is clamped to the physical window
$[200,\,3000]$ K to prevent unphysical
excursions during the transient. The clamp constructs the clamped value with
**zero ForwardDiff partials** (ForwardDiff `convert` semantics) — derivatives
do *not* ride along through a clamped iterate. This self-heals in the normal
path because a subsequent un-clamped Newton update restores the implicit-
function derivatives; only when the *final* iterate is pinned at 200/3000 K
(root outside the window) are the returned $T_s$/$q$ sensitivities zero —
which is then the mathematically correct clamp derivative, and the solver
warns (throttled) with the AD consequence spelled out.

For the physical BCs, convergence is fast and robust: convection contributes a
constant negative slope $-h_{\text{conv}}$ (one Newton step is exact), and
radiation contributes a monotone cubic slope $-4F\varepsilon_{\text{eff}}\sigma
T_s^3$ that Newton handles in 2–5 iterations. Combined BCs sum monotone
contributions, so the residual is monotone decreasing in $T_s$ and the root is
unique.

### 12.2.3 Automatic differentiation through the solve

The solver is written to carry whatever scalar type its inputs carry. It promotes
to the widest type among $(T_{\text{cell}}, k_{\text{eff}}, \Delta z_{1/2},
t, \varepsilon_{\text{eff}}, a_{\text{eff}}, \dot m_g)$; under sensitivity analysis (§17, Sensitivity
Analysis) this is a `ForwardDiff.Dual`. The Newton *convergence test* uses
`ForwardDiff.value(f)` (the primal residual), so the iteration count is identical
with or without AD, while the returned tuple $(T_s,\,q_{\text{applied}},\,
q_{\text{conducted}})$ is fully Dual-typed and carries sensitivities into the
Jacobian. One guard is needed: on the very first iteration under AD the
perturbed residual may pass the primal tolerance while the dual part has not yet
propagated, so convergence is suppressed for `dual_mode && i==1`. The diagnostic
`BoundaryState` stores plain `Float64` (it is bookkeeping, never differentiated).

### 12.2.4 Per-type flux and derivative

The Newton kernel calls `evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff,
radiation_model, face)` (the `BlowingFace` at the iterate, read only by convective
BCs, §12.3.3), which returns $(q_{\text{BC}}, dq_{\text{BC}}/dT_s)$ for each
concrete thermal BC. The exact pairs (verified against
`boundary_temperature.jl`) are collected below; the
$\varepsilon_{\text{eff}}$ argument is the composition-weighted effective surface
emissivity from the property cache and the keyword `a_eff` (default
$\varepsilon_{\text{eff}}$) the effective surface absorptivity for incident
radiation (§12.2.6).

| BC | $q_{\text{BC}}(T_s)$ | $dq_{\text{BC}}/dT_s$ |
|---|---|---|
| `AdiabaticBC` | $0$ | $0$ |
| `HeatFluxBC` | $q(t)$ (const or $q(t)$) | $0$ |
| `ConvectiveBC` | $h_{\text{conv}}(t)\,(T_\infty(t)-T_s)$ | $-h_{\text{conv}}(t)$ |
| `RadiativeBC` | $F\,\varepsilon\,\sigma\,(T_\infty^4-T_s^4)$ | $-4F\,\varepsilon\,\sigma\,T_s^3$ |
| `RadiativeBC` (under `P1_QUASI_STEADY`) | $0$ | $0$ |
| `RadiativeFluxBC` (surface modes) | $\alpha\,q_{\text{ext}}(t)$ | $0$ |
| `RadiativeFluxBC` (volumetric modes) | $0$ | $0$ |
| `DirichletThermalBC` (penalty) | $h_{\text{large}}\,(T_{\text{bc}}-T_s)$ | $-h_{\text{large}}$ |
| `CombinedThermalBC` | $\sum_m q_{\text{BC},m}$ | $\sum_m dq_{\text{BC},m}/dT_s$ |

Each BC is derived and motivated in §12.3.

**Neglected emissivity slope.** For a temperature-dependent emissivity
$\varepsilon=\varepsilon(T_s)$ the exact radiative derivative would be

$$
\frac{dq_{\text{BC}}}{dT_s}
= -4F\varepsilon(T_s)\sigma T_s^3
  \;+\; F\,\frac{d\varepsilon}{dT_s}\,\sigma\,(T_\infty^4 - T_s^4),
$$

but the code keeps only the first term. The
omitted term is part of the Newton *Jacobian* of the iteration, not the residual
itself — the residual still uses the true $\varepsilon(T_s)$ — so the converged
$T_s$ is unbiased; only the convergence rate degrades from quadratic to
super-linear when $d\varepsilon/dT$ is large. This is documented as a known
limitation (§12.12).

### 12.2.5 Application into the flux array and the sign flip

After the Newton solve, `_apply_thermal_bc_static!` writes the **conducted** flux
$q_{\text{conducted}}$ (not the applied flux) into the energy face-flux array,
with the boundary-orientation sign

$$
q[f] =
\begin{cases}
+q_{\text{conducted}}, & \text{bottom (}z=0\text{)},\\[2pt]
-q_{\text{conducted}}, & \text{top (}z=L\text{)},
\end{cases}
$$

(`finite_volume.jl`). Writing the conducted flux (rather than the applied
flux) is what makes the scheme energy-conserving: the divergence operator then
sees the flux that actually enters the first interior cell, and the surface energy
balance $q_{\text{BC}}=q_{\text{cond}}$ guarantees that the energy delivered at
the surface equals the energy conducted inward. The diagnostic check
`check_boundary_energy_balance` returns the relative imbalance
$|q_{\text{applied}}-q_{\text{conducted}}|/|q_{\text{applied}}|$ (the absolute
difference when $|q_{\text{applied}}|\le10^{-10}$)
per boundary, which is $\approx$ machine precision when Newton has converged
(`boundary_temperature.jl`).

Recall (§7, §13) that the interior conductive fluxes are explicitly **zeroed at
boundary faces** in `compute_conductive_fluxes!` before BC application, so the
overwrite here is the sole source of the boundary face flux; there is no
double-counting.

### 12.2.6 Effective surface emissivity and absorptivity (composition-weighted)

When a `RadiativeBC` is constructed with $\varepsilon=\texttt{nothing}$, the
re-radiation uses the boundary cell's effective emissivity
$\varepsilon_{\text{eff}}$, computed in the property cache as the volume-weighted
mean over the solid/liquid components with each component emissivity evaluated
at the cell temperature,

$$
\varepsilon_{\text{eff}} = \sum_j v_j\,\varepsilon_j(T_{\text{cell}}),
$$

clamped to $[0,1]$ (see §5.8.1). When a `RadiativeFluxBC` is constructed with
$\alpha=\texttt{nothing}$, the absorbed incident flux uses the boundary cell's
effective surface absorptivity

$$
a_{\text{eff}} = \sum_j v_j\,a_j,\qquad
a_j = \texttt{surface\_absorptivity}_j \text{ if set, else } \varepsilon_j(T_{\text{cell}}),
$$

(§5.8.2). Two physical principles are encoded here. First, **Kirchhoff's law**
$\alpha=\varepsilon$ at thermal equilibrium is the default: with no component
absorptivity set, $a_{\text{eff}}\equiv\varepsilon_{\text{eff}}$ and the same
quantity serves as both emissivity (re-radiation) and absorptivity (absorbed
incident flux). A component-level `absorptivity` relaxes this where the heater
spectrum differs from the surface's own emission spectrum, keeping the two
mixtures on the same volume weights. Second, surface emissivity is a
*geometric* property of the exposed mixture, so a volume-weighted (not
mass-weighted) average is appropriate and it **automatically tracks composition
changes** during pyrolysis (e.g. polymer $\to$ char as $v_j$ shifts). The
boundary solver reads these as `cache.ε[cell_id]` and `cache.a[cell_id]` and
threads them into `solve_boundary_temperature!` as the positional
$\varepsilon_{\text{eff}}$ and the keyword `a_eff` (`boundary_temperature.jl`);
the Jacobian operators recompute both on the Dual-seeded $(T,\xi)$, so
$\partial\varepsilon_{\text{eff}}/\partial T_{\text{cell}}$ and
$\partial a_{\text{eff}}/\partial T_{\text{cell}}$ from a temperature-dependent
component emissivity enter the surface energy row exactly.

### 12.2.7 `BoundaryState` storage

The per-boundary diagnostic record (`src/Geometry/mesh.jl`) holds the
solved surface temperature, the applied flux, the conducted flux, the inputs
$(T_{\text{cell}}, k_{\text{eff}}, \Delta z_{1/2})$, convergence metadata
(`converged`, `iterations`), and the blowing diagnostics at the converged $T_s$
(`blowing_mcp`, `blowing_ratio`, and the pooling-model fields `pooling_factor`,
`face_fraction`, `face_buoyancy`; §12.3.3). The fields are plain `Float64`
(`converged` a `Bool`, `iterations` an `Int`);
they are diagnostics and are never differentiated. The live (possibly Dual-typed)
quantities are returned from the solve and flow into the residual.

---

## 12.3 Thermal boundary conditions

This section derives each concrete thermal BC, its physical meaning, and how it
plugs into §12.2. All flux laws are signed **positive into the domain**.

### 12.3.1 `AdiabaticBC`

Perfect insulation: $q_{\text{BC}}=0$, $dq_{\text{BC}}/dT_s=0$. The surface energy
balance reduces to $q_{\text{cond}}=0$, hence $T_s=T_{\text{cell}}$ — the surface
floats to the adjacent cell temperature with zero conduction across the face.
Typical for the substrate face in a cone-calorimeter idealization with a
well-insulated back. (`boundary_fluxes.jl`,
`boundary_temperature.jl`.)

### 12.3.2 `HeatFluxBC`

A prescribed flux $q_{\text{BC}}=q(t)$, constant (`HeatFluxBC{<:Real}`) or a
function of time (`HeatFluxBC{<:Function}`), independent of $T_s$ so
$dq_{\text{BC}}/dT_s=0$. The Newton solve still runs to find the $T_s$ that
balances this against conduction. Sign: positive $q$ heats the material.
(`boundary_fluxes.jl`, `boundary_temperature.jl`.)

> **$T_s$ is a half-cell extrapolation.** For a pure prescribed-flux boundary
> the surface Newton reduces to the linear extrapolation
> $T_s = T_{\text{cell}} + q''\,\Delta z_{1/2}/k_{\text{eff}}$, which assumes
> the *entire* prescribed flux crosses the boundary half-cell by conduction.
> When a volumetric consumer sits inside that half-cell (a reaction zone under
> MW-scale flux, in-depth absorption), the reported $T_s$ is biased by up to
> tens of K at coarse resolution — a real but $n$-shrinking effect. It enters
> the Landau traveling-wave identity of V-CPL-02 through
> $v_{ss}(\text{measured}\ T_s)$: with the reaction zone under-resolved
> ($\Delta z \approx 2\delta_r$) the identity misses by ≈1.3 %, collapsing to
> ≈5×10⁻⁴ once $\Delta z \lesssim \delta_r$. (A *mesh-independent* residual
> of up to about 1 % has a different source, the wave-development transient of
> a uniform-$T_0$ start; the case measures after at least six development time
> constants to exclude it.) Interpret surface-temperature outputs under `HeatFluxBC`
> accordingly, and refine the near-surface mesh when $T_s$ itself is a target
> quantity.

### 12.3.3 `ConvectiveBC` (Newton's law of cooling)

$$
q_{\text{BC}} = h_{\text{conv}}(t)\,\bigl(T_\infty(t) - T_s\bigr),
\qquad
\frac{dq_{\text{BC}}}{dT_s} = -h_{\text{conv}}(t).
$$

Both $h_{\text{conv}}$ and $T_\infty$ may be constants or functions of $t$
(evaluated through `evaluate_param`). The flux is positive (heating) when
$T_\infty>T_s$ and negative (cooling) otherwise, so for a hot surface in a cool
ambient this is a convective loss. Because the derivative is the constant
$-h_{\text{conv}}$, the contribution to the Newton residual is linear and exactly
resolved in a single step. (`boundary_temperature.jl`.)

**Blowing.** Every `ConvectiveBC` carries a blowing model (`bc.blowing`, an
`AbstractBlowingModel`, default `NoBlowing()`), and the surface balance applies
it at every Newton iterate to the face state: the **signed outward** gas mass
flux of every component, read from the same `J_gas[:, f]` the mass-flux pass
populates (film-law diffusion plus, in Darcy–Fick mode, the pressure-BC venting
flux), the surface-temperature iterate and the material (`BlowingFace`,
`Discretization/blowing.jl`). The convective flux is $h\,\chi\,(T_\infty - T_s)$
with $\chi = h_{\text{eff}}/h$ from the model (`blowing_state`). It modifies
only the convective *coefficient*, not the gas enthalpy accounting, so it
composes exactly with the $S_{\text{conv}}$-only convention (§7.4.3) — no
double-count is possible. Zero outflow gives $\chi = 1$ exactly for every
shipped model, so the back face is unaffected automatically.

*Film law with a fitted constant* (`CouetteFilm(λ)`):

$$
\chi \;=\; \frac{\phi}{e^{\phi}-1},
\qquad
\phi \;=\; \frac{2\lambda\,\dot m''\,c_{p,g}}{h},
\qquad
\dot m''\,c_{p,g} \;=\; \sum_{g}\dot m_g\,c_{p,g}(T_s).
$$

$\lambda = 1/2$ is the classical Couette-film (Spalding) factor, $\phi = B =
\dot m'' c_{p,g}/h$; a smaller $\lambda$ fitted to a resolved flow weakens the
reduction at a given outflow and holds for that gas only. Limits in $B$:
$B = 0$ recovers the bare $h$ exactly; $B = 4$ (peak-PMMA outgassing against
canonical cone $h$) gives $h_{\text{eff}}/h = 0.075$; $B < 0$ (net inflow) is
the smooth suction enhancement.

*Property-ratio film for any light gas* (`PropertyRatioFilm(free_stream; n_Γ)`),
anchored on the Couette film of the free stream:

$$
\phi \;=\; \frac{\dot m''\,c_{p,\infty}(T_s)}{h}\left(\frac{\Gamma_\infty}{\Gamma_g}\right)^{n_\Gamma},
\qquad \Gamma = \frac{k}{c_p},
$$

with $\Gamma_\infty$ of the free stream and $\Gamma_g$ of the mixture leaving the
face (heat capacity weighted by mass fraction, conductivity by Wilke's rule
over the components with outward flux, each species' viscosity from Eucken's
relation when the component carries none). $n_\Gamma = 0.6$ is Kays's turbulent
exponent on the heat-capacity ratio, reproduced by the chamber CFD; with
$n_\Gamma = 0$ this is the Couette film in the free-stream heat capacity. Under
net inflow the arriving gas is the free stream and the ratio is 1.

*Pooling* (`PoolingFilm(free_stream; n_Γ, n_D, ψ0, Λc, Λw)`): the same
$\phi$ with a pooling factor $\psi$, the share of the no-outflow exchange the
face keeps,

$$
\chi \;=\; \frac{\phi}{e^{\phi/\psi}-1} \;=\; \psi\,f(\phi/\psi),
\qquad
\psi \;=\; \psi_0 + \frac{1-\psi_0}{1+\exp[-(\Lambda-\Lambda_c)/\Lambda_w]},
\qquad
\Lambda \;=\; \frac{1-(M_{\text{face}}/T_s)(T_\infty/M_\infty)}{1-T_\infty/T_s},
$$

where $f$ is the film factor and $\Lambda$ the face buoyancy relative to the
hot free stream, from the face molar mass $M_{\text{face}}$ of the species
mass fractions

$$
Y_i \;=\; \omega_i\left[1-\exp\!\left(-\frac{\dot m''}{\psi\,g\,s_i^{\,n_D}}\right)\right],
\qquad g = \frac{h}{\bar c_{p,\infty}},
\qquad s_i = \frac{D_{i,\infty}}{\alpha_\infty},
$$

with $\omega_i$ the species' share of the outflow, $g$ the gas-exchange
conductance without outflow (the film analogy: $h$ over the free-stream heat
capacity averaged between $T_\infty$ and $T_s$, a 5-point Gauss–Legendre mean),
$D_{i,\infty}$ the binary diffusivity into the free stream (Fuller, from the
component's diffusion volume on the model, else the Chapman–Enskog estimate
from its molar mass; a lumped `volatiles` gas states its own) and $n_D = 1/2$
the penetration-theory exponent. Components whose name matches the free
stream count in $\dot m''$ but not in the face composition. $\psi$ and the
$Y_i$ are solved together: $G(\psi) = \Psi(\Lambda(Y(\psi))) - \psi$ has
$G(\psi_0) > 0$ and $G(1) < 0$, so a 16-interval scan of $[\psi_0, 1]$ on the
primal values brackets every root, each refined by safeguarded Newton; where
two states are self-consistent the mean of the lowest and highest is taken
(the resolved flow alternates between them). One Newton step in the working
scalar type at the primal root then gives the implicit-function derivative,
so `ForwardDiff` partials ride through the solve. A face no hotter than the
free stream keeps $\psi = 1$. The three switch constants belong to the
apparatus flow, as $h$ does.

*User rule* (`FunctionBlowing(f)`): $\chi = f(\text{face}, h, T_\infty)$.

The factor is evaluated by the numerically stable `blowing_factor` (Taylor
branch at $|\phi|<10^{-5}$, overflow-safe $\phi\,e^{-\phi}/(1-e^{-\phi})$ form for
large $\phi$). The Newton derivative uses the frozen-$\chi$ form $-h\,\chi$ —
first-order consistent, so convergence is marginally sub-quadratic and the
converged root is exact. Each convective member of a `CombinedThermalBC` applies its own model
against its own $h$ and $T_\infty$; all other BC types are untouched. At the
converged $T_s$ the `BoundaryState` records the blowing conductance
`blowing_mcp`, the applied ratio `blowing_ratio` and, for the pooling model,
`pooling_factor`, `face_fraction` and `face_buoyancy` (the first convective
member's), and the solution extras carry the ratio and the pooling fields per
output time (`surface_blowing_ratio`, `surface_pooling_factor`,
`surface_face_fraction`, `surface_face_buoyancy`). The structured
Jacobian's `BCThermalOperator` and `BCConvectionOperator` rebuild $\dot m''$
from the stencil state via the shared boundary-face kernel, so the
$\partial(h\chi)/\partial(T_b,\xi_b)$ chain is exact (support-complete strict
verification passes with the film and pooling models in Fickian, Darcy–Fick
and ALE modes). Giving a convective BC a model on a case calibrated with the
bare coefficient will raise $T_s$ and MLR; introduce it together with flame
heat-flux feedback before any recalibration pass (the two partially cancel;
§1.8.1).
(`Materials/blowing_models.jl` types and gas data,
`Discretization/blowing.jl` evaluation.)

### 12.3.4 `RadiativeBC` (gray surface exchange with surroundings)

$$
q_{\text{BC}} = F\,\varepsilon\,\sigma\,\bigl(T_\infty^4 - T_s^4\bigr),
\qquad
\frac{dq_{\text{BC}}}{dT_s} = -4\,F\,\varepsilon\,\sigma\,T_s^3,
$$

with the Stefan–Boltzmann constant $\sigma=5.670374419\times10^{-8}$ W·m⁻²·K⁻⁴
(`Materials/constants.jl`, exported as `STEFAN_BOLTZMANN`). The fields are: the
**view factor** $F\in[0,1]$ (fraction of the surface's hemisphere that exchanges
with surroundings at $T_\infty$; default $1$), the **emissivity** $\varepsilon$
(either a supplied constant/function or, when `nothing`, the composition-weighted
$\varepsilon_{\text{eff}}$ of §12.2.6), and the surroundings temperature
$T_\infty$. This term is the surface's own re-radiation to and absorption from a
gray environment; it is the dominant cooling mechanism for a hot pyrolyzing
surface. The cubic temperature derivative is the source of the $\sim$20% error
that would result from naively evaluating radiation at the cell-center
temperature — and the reason the Newton solve exists. (`boundary_temperature.jl`.)

The neglected $d\varepsilon/dT$ term (§12.2.4) applies here when $\varepsilon$ is
temperature-dependent.

### 12.3.5 `RadiativeFluxBC` (external incident radiation; model-dependent)

This BC represents incident radiation from an *external* source (cone heater,
furnace, radiant panel) — not the surface's own emission. The absorbed fraction is
$\alpha\,q_{\text{ext}}$, where $q_{\text{ext}}=$ `flux` (constant or $q(t)$) and
the absorptivity $\alpha$ is either supplied or, when `nothing`, defaults to
the mixture surface absorptivity $a_{\text{eff}}$, which equals $\varepsilon_{\text{eff}}$
by Kirchhoff's law unless a component sets its own `absorptivity` (§12.2.6).

Its behavior is **dispatched on the global `radiation_model`** passed to
`solve()` (§8, Thermal Radiation). Note the defaults diverge: `solve()`
defaults to `NO_RADIATION`, while several internal helpers
(`apply_thermal_bc_flux!`, `update_boundary_state!`, `solve_boundary_temperature!`
and the solution extractors) default to
`SURFACE_ABSORPTION` — a trap for direct callers of the internals; always
pass `radiation_model` explicitly when bypassing `solve()`. The crucial
design point is that incident radiation must be deposited exactly once:

- **`NO_RADIATION` and `SURFACE_ABSORPTION`**: the radiation is
  treated as a *surface* heat flux. The Newton solve receives
  $$q_{\text{BC}} = \alpha\,q_{\text{ext}}(t), \qquad dq_{\text{BC}}/dT_s = 0,$$
  so the absorbed flux enters the surface energy balance and correctly raises
  $T_s$ above $T_{\text{cell}}$ when there is net heating
  (`boundary_temperature.jl`).
- **`BEER_LAMBERT` and `P1_QUASI_STEADY`**: the radiation is deposited
  *volumetrically* by `compute_radiation_source!` (Beer–Lambert) or by the P1
  residual path. To avoid double-counting, the boundary flux is forced to zero,
  $$q_{\text{BC}} = 0, \qquad dq_{\text{BC}}/dT_s = 0$$
  (`boundary_temperature.jl`, `boundary_fluxes.jl`).

A note on consistency. The *Newton-solve* evaluator
(`evaluate_bc_flux_and_derivative`) applies $\alpha q_{\text{ext}}$ as a boundary
flux for `SURFACE_ABSORPTION`, while `compute_radiation_source!` deposits *zero*
volumetric radiation for `SURFACE_ABSORPTION`. These are consistent:
`SURFACE_ABSORPTION` absorbs at the boundary face only, so the energy enters via
the boundary flux and the volumetric source is zero. Only `BEER_LAMBERT` (and the
experimental P1 path) routes energy through the volume, and only then is the
boundary flux suppressed. The diagnostic decomposition
`compute_boundary_heat_fluxes` follows the *same* dispatch: its `net` is the
flux the Newton balance actually applies — $\alpha q_{\text{ext}}$ is inside
`net` under `NO_RADIATION`/`SURFACE_ABSORPTION` and absent under
`BEER_LAMBERT`/P1, where the volumetric deposit is reported separately as
`extras.volumetric_absorption` (§12.9.1).

> **Implementation note.** The incident intensity is extracted at compile time
> from the (possibly combined) BC by `extract_radiative_flux_from_bcs(bc_set,
> Val(tag), t)`, which returns `(I, has_explicit_absorptivity)` — the flux with
> an explicit absorptivity already folded in, or the bare incident flux whose
> $a_{\text{eff}}$ scaling happens at injection
> (`boundary_conditions.jl`, §8.3.4). The intensity is passed as `I_radiation`
> to `compute_rhs!` together with the `radiation_kirchhoff` flag.

### 12.3.6 `DirichletThermalBC` (prescribed surface temperature)

A prescribed surface temperature $T=T_{\text{bc}}(t)$ is enforced by one of two
mechanisms, selected by how the BC appears, rather than a hard algebraic
constraint:

1. **In the Newton solve** (`boundary_temperature.jl`): a penalty "flux"
   with a very large conductance,
   $$q_{\text{BC}} = h_{\text{large}}\,(T_{\text{bc}} - T_s),\qquad
   dq_{\text{BC}}/dT_s = -h_{\text{large}}, \qquad h_{\text{large}}=10^{8}
   \;\text{W·m}^{-2}\text{·K}^{-1}.$$
   The large slope drives $T_s\to T_{\text{bc}}$ to high accuracy in one step. The
   value $10^8$ is large enough to pin $T_s$ to $T_{\text{bc}}$ within Newton
   tolerance while keeping the iteration well conditioned (a still-larger penalty
   would degrade conditioning with no accuracy benefit).
2. **Standalone, in the RHS assembly** (`finite_volume.jl`): a bare `DirichletThermalBC`
   is special-cased and bypasses the Newton path entirely. The surface temperature is the
   prescribed value, $T_s = T_{\text{bc}}(t)$, and the conducted flux is written directly,
   $$q_{\text{bc}} = k_{\text{cell}}\,\frac{T_{\text{bc}} - T_{\text{cell}}}{\Delta z_{1/2}},
   \qquad
   q[f]= \begin{cases}+q_{\text{bc}},&\text{bottom},\\ -q_{\text{bc}},&\text{top},\end{cases}$$
   with $\Delta z_{1/2}=\tfrac12\,V_{\text{cell}}/A_{\text{face}}$ — the boundary face's own
   area, identical to the Newton path's distance, so the two mechanisms agree under a
   variable cross-section. The bypass also updates the diagnostic `BoundaryState`
   ($T_{\text{surface}}=T_{\text{bc}}$, $q_{\text{applied}}=q_{\text{conducted}}=q_{\text{bc}}$,
   `converged = true`), so surface-temperature and flux diagnostics report the prescribed
   state. The penalty mechanism of item 1 is used only when the Dirichlet BC is *nested
   inside* a `CombinedThermalBC` (the composite goes through the Newton path; the two
   mechanisms agree to $\sim10^{-4}$ K).

The penalty approach is acknowledged as less clean than a direct constraint (a
genuine algebraic constraint would require a DAE index reduction); it is robust in
practice because the prescribed temperature pins the boundary tightly.

### 12.3.7 `CombinedThermalBC` (additive superposition)

Real surfaces experience several mechanisms at once. The `+` operator builds an
additive composite, e.g. the canonical cone-calorimeter top BC

```julia
bc = RadiativeFluxBC(flux = 50e3)            # absorbed incident cone flux
   + RadiativeBC(ε = nothing, T_inf = 300)   # surface re-radiation loss
   + ConvectiveBC(h = 10.0, T_inf = 300)     # convective loss
```

Fluxes and derivatives sum termwise,

$$
q_{\text{BC}} = \sum_m q_{\text{BC},m}(T_s,t),\qquad
\frac{dq_{\text{BC}}}{dT_s} = \sum_m \frac{dq_{\text{BC},m}}{dT_s}.
$$

The sum is built by a `@generated` function that unrolls the tuple at compile time
for type stability (`boundary_temperature.jl`). Because the radiation model can affect a
`RadiativeFluxBC` member, the combined evaluator threads `radiation_model` through
to every child. The `+` operator overloads handle every nesting case (BC+BC,
combined+BC, BC+combined, combined+combined) so associativity holds
(`boundary_conditions.jl`).

---

## 12.4 Mass (species) boundary conditions

The species boundaries supply the gas mass flux $J_j$ at each boundary face.
Unlike the thermal boundary, there is no surface-temperature nonlinearity, so the
flux is evaluated directly from the adjacent cell state. The internal diffusive
flux model is the volume-fraction (partial-pressure) gradient form of §9,
Gas-Phase Mass Transport & Pressure; here it is closed at the surface by an
external resistance.

### 12.4.1 `ImpermeableBC` (default)

$J_j = 0$ for all components — no gas crosses the boundary. This is the default
for both faces (`boundary_fluxes.jl`) and is the physically correct choice for
the substrate (sealed back face) and any non-venting surface. It is also the
implicit closure for solid/liquid components everywhere, since only gas components
are transported.

### 12.4.2 `MassFluxBC`

A prescribed species flux $J_j = J(t)$ (constant or function), positive into the
domain (`boundary_fluxes.jl`). Used to impose a known gas influx/efflux,
typically zero or negative (gas escaping).

### 12.4.3 `ConvectiveMassBC` (external film transfer)

$$
J_j = h_m(t)\,\bigl(\xi_{j,\infty} - \xi_{j,s}\bigr),
$$

the mass-transfer analog of Newton's cooling law, with the film coefficient $h_m$
(m·s⁻¹, velocity-based) and the ambient concentration $\xi_{j,\infty}$ (typically
$0$ for pyrolysis gas escaping to atmosphere) (`boundary_fluxes.jl`). The
coefficient may be estimated by the heat–mass analogy,

$$
h_m \approx \frac{h_{\text{conv}}}{\rho_{\text{air}}\,c_{p,\text{air}}},
$$

e.g. $h_{\text{conv}}\approx10$ W·m⁻²·K⁻¹ gives $h_m\approx0.008$ m·s⁻¹. With
$\xi_{j,\infty}=0$ and $\xi_{j,s}>0$, $J_j<0$ (escaping), consistent with the
flux-into-domain sign.

The ambient $\xi_{j,\infty}$ accepts three forms, resolved per gas species $j$ by
`ambient_concentration(bc, t, j)`: a `Real` (the same ambient for every gas species —
typically 0, the atmosphere carrying no pyrolyzate), a `Function` of $t$ (time-varying
scalar, broadcast to every gas species), or an `NTuple{NC}` of per-component ambient
concentrations indexed by the material's component order (condensed-phase entries are
ignored; entries are validated non-negative at construction and the tuple length against
$N_C$ at problem validation). The per-species form is what represents ambient air ingress
without injecting pyrolyzate: a scalar $\xi_\infty>0$ would feed that concentration of
*every* gas species into the domain.

### 12.4.4 `DiffusiveMassBC` (combined internal/external resistance)

This BC couples the **internal** Fickian diffusion resistance from the surface
cell to the face with an **external** film resistance, in series. Working in the
volume-fraction formulation $J=-\lambda\,\partial\xi/\partial x$ with the surface
gradient approximated as $(\xi_s-0)/d$ (zero exterior concentration), the two
resistances are

$$
R_{\text{internal}} = \frac{d}{\lambda},
\qquad
R_{\text{external}} = \frac{1}{h_m},
$$

where $d=\Delta z_{1/2}$ is the cell-center-to-face distance and $\lambda$ is the
effective gas-transfer coefficient. The combined flux is

$$
J_j = -\frac{\xi_{j,s}}{R_{\text{internal}}+R_{\text{external}}}
    = -\frac{\lambda\,h_m\,\xi_{j,s}}{\lambda + h_m\,d}.
$$

The minus sign places escaping gas at $J_j<0$ (flux-into-domain convention). Two
limits follow:

- **No external resistance** ($h_m\to\infty$, the default `DiffusiveMassBC()` with
  `h_m = Inf`): pure internal diffusion,
  $$J_j = -\frac{\lambda\,\xi_{j,s}}{d};$$
- **Impermeable** ($h_m\to0$): $J_j\to0$.

> **Implementation note.** The flux law is in `mass_bc_flux(bc::DiffusiveMassBC,
> ξ_surface, λρ, ρ, d, t, j)` (`boundary_fluxes.jl`) — every mass-BC flux method carries
> the trailing species index `j` so per-species ambients resolve correctly. The
> discretization passes the *combined transport coefficient* $\lambda\rho$ (consistent
> with the rigorous gas-transport closure) and the density $\rho$; for $\rho \le 0$ the
> flux is zero, and otherwise the function recovers $\lambda = \lambda\rho/\rho$ by exact
> division (no floor — the callers construct $\lambda\rho$ as the product, so division is
> exact for any positive density), because $h_m$ is defined per concentration, not per
> mass fraction. The surface-cell quantities $(\xi_{j,s},\,\lambda,\,\rho,\,d)$ are read
> in `_apply_mass_bc_static!` (`finite_volume.jl`).

**Sign reconciliation.** Three sign conventions coexist and must be reconciled
(`boundary_fluxes.jl`):

| Quantity | Convention |
|---|---|
| `mass_bc_flux(...)` | **negative** = gas escaping (flux-into-domain) |
| `compute_surface_mass_flux(...)` (diagnostic) | **positive** = gas escaping |
| `J[face]` stored in the flux array | **positive** = flux in $+z$ |

The discretization bridges these by orienting per boundary
(`_apply_mass_bc_static!`, `finite_volume.jl`): at the bottom
$J[j,f]\mathrel{+}=J_{\text{bc}}$, at the top $J[j,f]\mathrel{-}=J_{\text{bc}}$.
With $J_{\text{bc}}<0$ for escaping gas, the top-face stored flux is positive —
gas leaving in $+z$ at the exposed surface — exactly what the divergence
operator needs. The application **accumulates** rather than overwrites: in
Darcy–Fick mode the boundary face already carries the pressure-BC advective
venting flux (§9.5.3), and the film-law flux adds to it; in Fickian mode the
boundary faces arrive zero so the accumulation is equivalent to assignment.

### 12.4.5 `CombinedMassBC`

Like the thermal composite, $J_j=\sum_m J_{j,m}$, built by a `@generated` unrolled
sum (`boundary_fluxes.jl`), e.g. `DiffusiveMassBC(Inf) +
ConvectiveMassBC(h_m=0.01, ξ_inf=0.0)`. Pressure BCs (§12.5) participate in the
composite by contributing zero *diffusive* flux (they act on the Darcy velocity,
not the diffusive flux), so they can be combined freely with `DiffusiveMassBC`
(`boundary_fluxes.jl`).

---

## 12.5 Pressure boundary conditions (Darcy flow)

In Darcy–Fick transport mode (§9), gas advection is driven by the Darcy velocity
$v_g = -(\kappa/\mu)\nabla P$, and the pressure boundaries set $v_g$ at the
boundary faces. These BCs are subtypes of `AbstractMassBC` so they live in the
mass slot of the `BoundaryConditionSet`; they do **not** contribute a diffusive
flux (returning $0$ from `mass_bc_flux`, §12.4.5) — they act through
`apply_pressure_bc_to_darcy_velocity!` (`darcy_velocity.jl`), which runs *after*
the interior Darcy velocities are computed and sets the boundary face velocity
$v_b$. That velocity then enters the species equations as the boundary advective
venting flux $J_j=v_b\,\xi_j^{\text{int}}/\varphi_{\text{int}}$ (outflow
only; §9.5.3), written by `compute_gas_fluxes_darcy_fick!` before the diffusive
mass BC is added.

### 12.5.1 `DirichletPressureBC` (prescribed pressure)

Prescribes the boundary pressure $P_{\text{bc}}$ (default $101325$ Pa,
atmospheric). The boundary face velocity follows from Darcy's law with a one-sided
pressure gradient between the boundary and the interior cell:

$$
v_g = -\frac{\kappa_{\text{int}}}{\mu_{\text{int}}}\,\nabla P,
\qquad
\nabla P =
\begin{cases}
(P_{\text{bc}} - P_{\text{int}})/d, & \text{top},\\[2pt]
(P_{\text{int}} - P_{\text{bc}})/d, & \text{bottom},
\end{cases}
$$

with $d$ the cell-center-to-face distance, $\mu_{\text{int}}$ the material's
pore-gas viscosity at the cell temperature (§5.9), and the velocity clamped to zero if
$\kappa_{\text{int}}\le0$ or $\mu_{\text{int}}<10^{-20}$ Pa·s
(`darcy_velocity.jl`). The orientation of the gradient ensures the correct
outflow direction at each face: at the top a higher interior pressure
($P_{\text{int}}>P_{\text{bc}}$) yields $v_g>0$ (outflow in $+z$); at the bottom
the same drives $v_g<0$ (outflow in $-z$).

### 12.5.2 `NeumannPressureBC` (prescribed gradient)

Prescribes the pressure gradient $\partial P/\partial n$ along the **outward
normal** (default $0$ = impermeable). For zero gradient the boundary velocity is
set to zero; for a nonzero prescribed gradient,
$v_g=-(\kappa_{\text{int}}/\mu_{\text{int}})\,\partial P/\partial n$ at the
top and $v_g=+(\kappa_{\text{int}}/\mu_{\text{int}})\,\partial P/\partial n$
at the bottom (where $\hat n=-\hat z$, so $\partial P/\partial z=-\partial
P/\partial n$), with the same $\kappa,\mu$ guards (`boundary_darcy_velocity`,
`darcy_velocity.jl`). The zero-gradient default makes a `NeumannPressureBC()`
equivalent to an impermeable wall for Darcy flow.

### 12.5.3 `ConvectivePressureBC` (Robin / surface resistance)

A Robin-type law for a surface with finite resistance to gas escape (char layer,
protective coating, blowing resistance):

$$
v_g\cdot\hat n = h_P\,\bigl(P_{\text{int}} - P_{\text{ext}}\bigr),
$$

with $h_P$ the pressure-transfer coefficient (m·Pa⁻¹·s⁻¹) and $P_{\text{ext}}$ the
exterior pressure, gated on $\kappa_{\text{int}}>0$ — an impermeable surface
cell cannot vent regardless of $h_P$. The implementation orients the outward
normal per boundary (`boundary_darcy_velocity`, `darcy_velocity.jl`):

$$
v_g[f] =
\begin{cases}
+h_P\,(P_{\text{int}}-P_{\text{ext}}), & \text{top},\\[2pt]
-h_P\,(P_{\text{int}}-P_{\text{ext}}), & \text{bottom}.
\end{cases}
$$

The limits bracket the other two pressure BCs: $h_P\to\infty$ pins the boundary
pressure (Dirichlet) and $h_P\to0$ seals it (zero-gradient Neumann / impermeable).

> **Implementation note.** All three pressure BCs are scalar methods of the
> shared kernel `boundary_darcy_velocity(bc, P_int, κ_int, μ_int, d, Val(tag))`
> (`darcy_velocity.jl`); a `@generated` overload sums the contributions of a
> `CombinedMassBC`, and the `AbstractMassBC` fallback contributes $0$ (e.g. a
> `DiffusiveMassBC` in the mass slot does not set a Darcy velocity). The
> face loop `_apply_pressure_bc_faces!` writes the kernel value per boundary
> face; the orchestration entry `apply_pressure_bc_to_darcy_velocity!` applies
> top then bottom. The same kernel is evaluated (with `Dual` arguments) inside
> the `BCMassOperator` Jacobian and (with the saved state) in the MLR
> diagnostics, so the three paths cannot drift apart.

---

## 12.6 Mesh-motion boundary conditions (ALE)

In ALE mode (§11, The ALE Framework) the node positions $z_i$ are state variables
and the mesh-velocity boundary treatment is **fixed by the formulation, not
user-configurable**: the bottom node is pinned ($w_1 = 0$, hard-coded in the
residual's mesh-velocity integration) and the velocity at every other node
follows cumulatively from the species dilation $\theta_i$ and the column-global
radial strain $\theta_A$, $w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i$
(§11.3, §15.3). The surface node is therefore automatically material-following
(Lagrangian) — it recedes or expands with the integrated dilation.

There are no user-facing mesh-motion BC types: `BoundaryConditionSet` carries
only thermal and mass slots, and the ALE bottom anchor ($w_1 = 0$) is
hard-coded in the residual's mesh-velocity assembly.
See §11 for the full ALE mesh-velocity derivation and the geometric
conservation law.

---

## 12.7 The `BoundaryConditionSet` container

### 12.7.1 Structure and tag keying

All boundaries are bundled into a single type-stable container
(`boundary_conditions.jl`):

```julia
struct BoundaryConditionSet{TT<:NamedTuple, MT<:NamedTuple}
    thermal::TT   # NamedTuple{tags, Tuple{<:AbstractThermalBC,...}}
    mass::MT      # NamedTuple{tags, Tuple{<:AbstractMassBC,...}}
end
```

Each boundary tag maps to exactly one thermal BC and one mass BC. The 1D
convention uses `:top` ($z=L$) and `:bottom` ($z=0$); the field shape is designed
so that a future 2D/3D constructor can add tags without changing the runtime
contract. Storing the BCs in a typed `NamedTuple` (rather than, say, a `Dict`)
means a symbol lookup compiles to a direct field read and the residual specializes
statically.

### 12.7.2 Constructor, defaults, and validation

The keyword constructor (`boundary_conditions.jl`) merges the
user-supplied partial NamedTuples with the 1D defaults

$$
\texttt{thermal} = (\,\texttt{top}=\texttt{AdiabaticBC()},\ \texttt{bottom}=\texttt{AdiabaticBC()}\,),
\qquad
\texttt{mass} = (\,\texttt{top}=\texttt{ImpermeableBC()},\ \texttt{bottom}=\texttt{ImpermeableBC()}\,),
$$

so a caller may specify only the slots it cares about. After merging it
**validates** that every thermal payload is an `AbstractThermalBC` and every mass
payload is an `AbstractMassBC`, throwing an `ArgumentError` if a BC was placed in
the wrong slot. The default — adiabatic and
impermeable on both faces — is the inert closure: a column with no defined BCs
neither gains energy nor loses mass.

### 12.7.3 Compile-time accessors

Access is through `Val{tag}` dispatch:

```julia
@inline get_thermal_bc(bc_set, ::Val{tag}) where {tag} = bc_set.thermal[tag]
@inline get_mass_bc(bc_set,    ::Val{tag}) where {tag} = bc_set.mass[tag]
```

At a call site such as `get_thermal_bc(bc_set, Val(:top))`, the `Val(:top)`
resolves at compile time to a direct NamedTuple field read with zero overhead
(`boundary_conditions.jl`). The RHS assembly uses these accessors so that
each boundary's BC type is known statically and the per-type flux law is inlined.
The tag$\to$id bridge `tag_to_boundary_id` (§12.1.3) connects the user-facing
symbol API to the integer-keyed `Mesh1D.boundary_faces` / `boundary_states`
dictionaries.

### 12.7.4 Constructor parameter validation

Individual BC constructors optionally validate their physical parameters at build
time (`validate_bc_params`, `boundary_conditions.jl`): $h\ge0$, $T>0$ K,
$\varepsilon\in[0,1]$, $F\in[0,1]$, $\alpha\in[0,1]$, $h_m\ge0$, $P>0$.
Function-valued parameters cannot be checked at construction and are skipped
(`_get_scalar_value` returns `nothing` for functions and property functions).

---

## 12.8 Initial conditions

### 12.8.1 The state seed

The initial state is specified by a temperature closure $T_{\text{initial}}(z)$
and per-component concentration closures $\xi_{j,\text{initial}}(z)$, evaluated at
each cell center $z_i^c$ (`initial_conditions.jl`):

$$
T_i(0) = T_{\text{initial}}(z_i^c), \qquad
\xi_{j,i}(0) = \xi_{j,\text{initial}}(z_i^c), \quad j=1,\dots,N_C .
$$

The closures are evaluated and stored in the mutable mesh's cell states, then
packed into the flat state vector (§15, Temporal Integration & Structured
Jacobian). Cell volumes are not altered by initialization.

The single problem-level entry point is `apply_initial_conditions!(problem)`,
which performs the cell-state write above **and** populates the dry-solid
denominator of §12.8.3 (the populator reads the freshly written cell states,
so the two steps are ordered). `solve()` and the exported composition helper
`initialize_state(problem; …) -> (def, ws, u0)` both call it before packing
`u0`; the mesh-level utility `apply_initial_conditions!(mesh, material; T, ξ)`
writes only uniform cell states and deliberately does **not** touch
`m_dry_solid_initial`. Composing the pure-residual layer by hand without the
problem-level IC step packs the mesh constructor's placeholder state instead
of the problem's initial conditions (§15.2) — a resolution-independent error
with a perfectly plausible-looking solution.

### 12.8.2 Auto-defaults

If the user supplies no initial conditions, the `PyrolysisProblem` constructor
synthesizes them (`pyrolysis_problem.jl`):

- $T_{\text{initial}}(z)=300$ K uniformly;
- $\xi_{1,\text{initial}}(z)=\rho_1(300\,\text{K})$ if the first component is a
  solid (the virgin solid at its reference density), and $\xi_{j,\text{initial}}=0$
  for all other components.

This "first solid at reference density, everything else zero" default models a
fresh, uniform virgin slab at room temperature — the standard cone/TGA starting
point. The uniform-IC variant `apply_initial_conditions!(mesh, material; T, ξ)`
follows the same rule when $\xi=\texttt{nothing}$, with the first solid's density
evaluated at the given $T$ rather than at 300 K (`initial_conditions.jl`). Profile helpers (`step_profile`,
`linear_profile`, `gaussian_profile`, `tanh_profile`) are provided for non-uniform
seeds (e.g. a pre-heated layer).

### 12.8.3 The dry-solid mass denominator $m_{\text{dry},i}$

A second, BC-adjacent quantity is fixed at initialization: the per-cell initial
dry-solid mass, which is the constant denominator of the pyrolysis-progress state
$\chi_i = m_{\text{released,dry},i}/m_{\text{dry},i}$ used by the variable
cross-section law (§10, Volume Change & Lateral Shrinkage). It is populated by
`populate_m_dry_solid_initial!` (`initial_conditions.jl`), wired into
`apply_initial_conditions!(problem)` after the cell-state write:

$$
m_{\text{dry},i} = V_i(0)\,\sum_{j\in\text{solid}} \xi_{j,i}(0),
$$

summing only the solid components (the `is_solid` predicate — not liquids, not
gases) at the initial state. This is computed once and stored on the mutable mesh
as `Float64`; under AD the $\chi$ numerator is Dual-typed but this denominator
stays `Float64` (see §15, §17).

---

## 12.9 Boundary-flux diagnostics

### 12.9.1 Flux decomposition

For auditing and post-processing, `compute_boundary_heat_fluxes(bc, T_s, t,
ε_eff, radiation_model[, face])` decomposes the boundary flux into named
components (`boundary_fluxes.jl`). The optional trailing `face` is the
`BlowingFace` at the converged $T_s$ (outward species fluxes, surface
temperature, material); when supplied, the convective component is reported
with the same ratio $\chi$ its blowing model gave the surface Newton
(§12.3.3), keeping the decomposition consistent with `net`. Components:

- `absorbed` — incident radiation absorbed *at the boundary face* [W·m⁻²];
- `reradiation` — re-radiation loss (positive = leaving);
- `convection` — convective loss (positive = leaving);
- `prescribed` — prescribed flux;
- `net` — the flux the surface energy balance actually applies at the face
  (positive into the domain).

For `ConvectiveBC` and `RadiativeBC`, the loss components are reported with the
*opposite* sign to the into-domain convention (so a hot surface re-radiating shows
a positive `reradiation`). Per component, `net` matches
`evaluate_bc_flux_and_derivative` — the Newton-applied flux — for every
radiation model, `DirichletThermalBC` excepted (its decomposition is all zero; for a
bare Dirichlet BC the saved extras report the solved $q_{\text{conducted}}$ as the
net flux), so the ledger identity
$\texttt{net} = q_{\text{cond}}$ holds at convergence
(§12.9.2, §16.6). For combined BCs the components **and** `net` are summed over
children; `net` is never recomputed from the components. Consequences of the
model dispatch:

- `NO_RADIATION` / `SURFACE_ABSORPTION`: a `RadiativeFluxBC` contributes
  `absorbed = net = α·q_ext` (α the explicit BC absorptivity or a_eff) — the
  flux is applied at the face.
- `BEER_LAMBERT` / P1: a `RadiativeFluxBC` contributes zero `absorbed` and zero
  `net` — the radiation enters volumetrically (§8.3) and is reported separately
  as `extras.volumetric_absorption` $= \alpha_{\text{surf}}\,q_{\text{ext}}\,
  (1-e^{-\tau_{\text{tot}}})$, the exact telescoped sum of the Beer–Lambert
  deposit on the saved state.

### 12.9.2 Energy-balance check

`check_boundary_energy_balance(mesh)` returns, per boundary, the relative
imbalance $|q_{\text{applied}}-q_{\text{conducted}}|/|q_{\text{applied}}|$ (the
absolute difference when $|q_{\text{applied}}|\le10^{-10}$), which should be
$\approx$ machine precision ($\sim10^{-15}$) when the
Newton solve has converged (`boundary_temperature.jl`). This is a direct
continuous self-check that the surface energy balance $q_{\text{BC}}=q_{\text{cond}}$
holds. The accessors `get_surface_temperature(mesh)` and
`get_back_face_temperature(mesh)` return the solved $T_s$ at the top (id 1) and
bottom (id 2) boundaries (`boundary_temperature.jl`).

---

## 12.10 Worked example: the cone-calorimeter top boundary

Assembling the canonical top boundary makes the pieces concrete. A 50 kW·m⁻² cone
with surface re-radiation and natural convection to a 300 K ambient, on a
sealed-but-venting slab, is

```julia
bc_set = BoundaryConditionSet(
    thermal = (
        top    = RadiativeFluxBC(flux = 50e3)            # absorbed incident, α = ε_eff
               + RadiativeBC(ε = nothing, T_inf = 300)   # re-radiation loss
               + ConvectiveBC(h = 10.0, T_inf = 300),    # convective loss
        bottom = AdiabaticBC(),                          # insulated back
    ),
    mass = (
        top    = DiffusiveMassBC(),                      # gas vents, internal diffusion limit
        bottom = ImpermeableBC(),                        # sealed back
    ),
)
```

With `radiation_model = SURFACE_ABSORPTION`, the surface energy
balance solved by Newton at the top is

$$
\underbrace{\alpha\,(50{,}000)}_{\text{absorbed}}
\;-\;\underbrace{\varepsilon_{\text{eff}}\,\sigma\,(T_s^4 - 300^4)}_{\text{re-radiation}}
\;-\;\underbrace{10\,(T_s - 300)}_{\text{convection}}
\;=\;
k_{\text{eff}}\,\frac{T_s - T_{\text{cell}}}{\Delta z_{1/2}},
$$

with $\alpha=\varepsilon_{\text{eff}}$ (Kirchhoff, both `nothing`). As the surface
chars, $\varepsilon_{\text{eff}}$ tracks the composition automatically, raising
both the absorbed and re-radiated terms. The Newton derivative is
$f'(T_s) = -4\varepsilon_{\text{eff}}\sigma T_s^3 - 10 - k_{\text{eff}}/\Delta
z_{1/2}$ (the absorbed term and convection-ambient are $T_s$-independent), which is
strictly negative — guaranteeing a unique root reached in a few iterations. The
conducted flux $q_{\text{conducted}}$ is written to the top face as
$q[f_{\text{top}}]=-q_{\text{conducted}}$, and the vented gas leaves the top with
$\dot m_g$ accumulated by `DiffusiveMassBC` and oriented to $+z$ via the top-face
sign flip. If `radiation_model = BEER_LAMBERT` is selected instead, the
`RadiativeFluxBC` contributes zero to this balance and the 50 kW·m⁻² enters the
energy equation as a volumetric Beer–Lambert source $q'''=\alpha I$ (§8).

---

## 12.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The three reference solvers express the same surface physics in their own
notation; the subsections below map their symbols to ours on first use and
contrast the modeling choices. (ThermaKin/Gpyro publish $h>0=$ exothermic; our
internal convention is $h>0=$ endothermic with $Q_{\text{rxn}}=-\sum h_r r_r$, but
this sign issue is confined to the volumetric reaction heat, not the surface BCs.)

### 12.11.1 ThermaKin2Ds

ThermaKin2Ds applies piecewise external/flame boundary conditions and supports an
external radiant source with surface re-radiation. The exposed-face energy
balance, in ThermaKin's form, is

$$
-k\frac{\partial T}{\partial x}\Big|_{\text{surf}}
= \alpha\,\dot q''_{\text{ex}}
  - \varepsilon\sigma\bigl(T_s^4 - T_\infty^4\bigr)
  - h_{\text{conv}}\bigl(T_s - T_\infty\bigr),
$$

with $\dot q''_{\text{ex}}$ the incident external flux, $\alpha$ the absorptivity
and $\varepsilon$ the emissivity (Kirchhoff $\alpha=\varepsilon$). Pyrolysis.jl's
default cone BC `RadiativeFluxBC + RadiativeBC + ConvectiveBC` is exactly this
balance, with our $q_{\text{ext}}=\dot q''_{\text{ex}}$ and our
$\varepsilon_{\text{eff}}$ playing the role of ThermaKin's surface $\varepsilon$.
ThermaKin also tracks a reflected component $I_{\text{ex}}^0$; our
composition-weighted $\varepsilon_{\text{eff}}$ supplies
$\alpha=\varepsilon_{\text{eff}}$ and the reflected fraction $(1-\alpha)$ is
implicitly discarded (opaque-surface assumption), consistent with the
SURFACE_ABSORPTION model. The shared idea — that both the absorbed external flux
and the surface re-radiation depend on the *true* surface temperature — is
realized by ThermaKin through its Crank–Nicolson surface node and by Pyrolysis.jl
through the explicit Newton surface-balance solve of §12.2. ThermaKin carries gas
enthalpy via a convective term $c_g(J_g\cdot\partial T/\partial x)$ (its Eq. 17),
which is precisely our $S_{\text{conv}}$ and the reason the surface balance carries
no gas-enthalpy term (§7.4.3).

### 12.11.2 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) writes the exposed-surface balance as a
radiative–convective balance with an impermeable back face by default. In Gpyro
notation,

$$
-\bar k\frac{\partial T}{\partial z}\Big|_{\text{surf}}
= \dot q''_e - \varepsilon\sigma\bigl(T_s^4 - T_\infty^4\bigr)
  - h_c\bigl(T_s - T_\infty\bigr),
$$

and gas escape at the surface uses a heat–mass transfer analogy for the film
coefficient, $h_m \sim h_c/(\rho c_p)$ — the same estimate we recommend for
`ConvectiveMassBC` (§12.4.3). Gpyro's control-volume formulation places the
surface temperature at the face via its harmonic interface interpolation; our
Newton solve achieves the equivalent by solving the one-sided Fourier balance for
$T_s$. Gpyro permits a non-zero back-face pressure/mass boundary in its
pressure-evolution model; our `DirichletPressureBC` / `NeumannPressureBC` /
`ConvectivePressureBC` (§12.5) are the analogous Darcy boundary closures, with
Gpyro's permeability $K$ mapping to our $\kappa$. Gpyro's default impermeable back
face corresponds to our default `ImpermeableBC` + `NeumannPressureBC()` (zero
gradient).

### 12.11.3 FDS

The FDS solid-phase model imposes the 1D surface boundary condition

$$
-k\frac{\partial T_s}{\partial x}\Big|_{x=0}
= \dot q''_c + \dot q''_r,
$$

i.e. the conducted flux equals the sum of the convective and net radiative surface
fluxes, with $\dot q''_r$ including absorbed incident and emitted re-radiation.
This is the same content as our $q_{\text{cond}}=q_{\text{BC}}$ with
$q_{\text{BC}}=q_{\text{ext}}-\varepsilon\sigma T_s^4-h_{\text{conv}}(T_s-
T_\infty)$. FDS additionally performs a *surface oxygen* Newton iteration for
oxidative pyrolysis (coupling surface $O_2$ to the gas phase); Pyrolysis.jl is a
condensed-phase-only solver with no gas-phase $O_2$ feedback, so our surface
Newton iterates on $T_s$ alone. FDS recovers Arrhenius parameters from TGA
reference temperature/rate (a kinetics, not a BC, concern; see §6, §18). The FDS
in-depth radiation (a two-flux Schuster–Schwarzschild model) corresponds to our
BEER_LAMBERT / P1 volumetric handling, under which the `RadiativeFluxBC` boundary
flux is suppressed to avoid double-counting (§12.3.5), just as FDS deposits the
in-depth absorbed fraction volumetrically rather than at the face.

---

## 12.12 Limitations and open issues

1. **Temperature-dependent emissivity derivative.** For a `RadiativeBC` whose
   own `ε` is a function of $T_s$, the radiative Newton derivative omits the
   $F\,(d\varepsilon/dT)\,\sigma\,(T_\infty^4-T_s^4)$ term (§12.2.4). The
   converged $T_s$ is unbiased, but convergence can become sub-quadratic for
   strongly temperature-sensitive $\varepsilon(T_s)$. This does not affect
   component-level $\varepsilon_j(T)$: the mixture $\varepsilon_{\text{eff}}$
   is evaluated at the *cell* temperature and is a constant of the surface
   Newton iteration, so that derivative is exact.

2. **Dirichlet thermal BC via penalty.** `DirichletThermalBC` uses a large
   conductance $h_{\text{large}}=10^8$ in the Newton solve and a special-cased
   direct flux in the RHS (§12.3.6). A genuine algebraic constraint (DAE) would be
   cleaner and would remove the residual conditioning sensitivity of the penalty.

3. **View-factor meaning in 1D.** The `RadiativeBC` view factor $F$ is a scalar
   surrogate for a 3D geometric quantity; in a strictly 1D model it is best read as
   a tuning factor for partial enclosure rather than a rigorously computed exchange
   factor.

4. **Mass-transfer coefficient estimation.** The film coefficient $h_m$ for
   `ConvectiveMassBC` / `DiffusiveMassBC` is supplied by the user, typically via
   the heat–mass analogy $h_m\approx h_{\text{conv}}/(\rho_{\text{air}}
   c_{p,\text{air}})$; there is no built-in correlation, so its accuracy depends on
   the user's estimate of $h_{\text{conv}}$ and the flow regime.

5. **Pressure-BC scope.** The three pressure BCs act only under Darcy–Fick
   transport, i.e. when a condensed component is permeable (`has_darcy_flow`); on a
   Fickian material they are accepted and contribute nothing. Their venting is
   outflow-only: the ambient ghost carries none of the tracked gas species, so an
   exterior pressure above the interior pressure draws no gas into the slab
   (§9.5.3). The residual, the Jacobian and the saved-output diagnostics evaluate
   the one kernel `boundary_darcy_velocity` (§12.5).
