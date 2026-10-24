# 12. Boundary and Initial Conditions

This chapter specifies how Pyrolysis.jl closes the coupled energy, species, and
gas-pressure equations at the boundary faces of the 1D column ($z=0$
substrate / `:bottom`, $z=L$ exposed surface / `:top`) and how it seeds the
initial state. It derives the surface energy balance and the Newton iteration
that solves for the surface temperature $T_s$ in a cell-centered
finite-volume scheme, and the closed forms of the adiabatic, prescribed-flux and
prescribed-temperature faces; documents every concrete thermal, mass, pressure, and
mesh-motion boundary condition with its exact flux law, temperature derivative,
and sign convention as implemented; explains the tag-keyed
`BoundaryConditionSet` and the per-tag boundary patches of the kernels; and treats the
initial fields at the cell centres and their defaults. Throughout, the equations are
checked term-by-term against the source so that the math and the code agree
exactly — including the outward boundary normals, the surface balance evaluated at
the trial clamp of the surface temperature, the Kirchhoff emissivity coupling, the
neglected $d\varepsilon/dT$ term, and why the surface balance carries no
gas-enthalpy term (the gas leaving through the face takes its enthalpy with the
face's energy flux).

This chapter conforms to the unified nomenclature of §2, Nomenclature, and uses
the binding sign conventions: an interior face flux is positive in the $+z$
direction (bottom$\to$top), and a boundary face flux is positive along the
face's **outward** normal, at both ends of the column; each face flux leaves its
owner cell and enters its neighbour, so a positive boundary flux removes the
quantity from the domain; and for heats of reaction $h>0$ is endothermic
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
= -\nabla\!\cdot\mathbf{q}_{\text{cond}} - \sum_j \mathbf J_j\cdot\nabla h_j
   + Q_{\text{rad}} + Q_{\text{rxn}},
\qquad
\frac{\partial \xi_j}{\partial t}
= -\nabla\!\cdot\mathbf{J}_j + S_{j},
$$

the second term of the energy equation being the sensible enthalpy the gas
carries relative to the skeleton (§7.4). The cell-centered finite-volume
discretization (Chapter 13) integrates each equation over a cell $i$: the rows are
extensive (the cell masses, and the cell temperature times the cell's heat
capacity), and every face flux $F_f$, stored per unit area along the face normal,
is scattered as $-F_fA_f$ into the face's owner and $+F_fA_f$ into its neighbour.
For an interior face the flux is a constitutive law evaluated from the two
adjacent cells (Fourier, Fick, Darcy–Fick). A **boundary face** has no neighbour,
so its flux must be supplied by a boundary condition. The boundary-condition
subsystem therefore has exactly one job: for each live boundary face, given the
state of its cell and the time $t$, produce the outward energy flux, the outward
gas mass fluxes, and (in Darcy mode) the outward venting the pressure condition
drives.

There is one subtlety unique to the **thermal** boundary: because temperature is
stored at cell centers and radiative and convective fluxes depend nonlinearly on
the *surface* temperature $T_s$ (which is not a state variable), the code does not
simply evaluate the applied flux at the adjacent cell-center temperature. Doing so
produces an $\sim$20 % energy imbalance for $T^4$ radiation. Instead it solves a
small surface energy balance for $T_s$ (§12.2). Three conditions need no iteration: an adiabatic face, a prescribed
flux and a prescribed temperature fix the face flux and $T_s$ in closed form
(§12.3.1, §12.3.2, §12.3.6). The species and pressure boundaries have no such
nonlinearity and are evaluated directly.

### 12.1.2 Type hierarchy and dispatch

Boundary conditions are organized by physical kind, each rooted in an abstract
supertype declared in the `Materials` module:

- `AbstractThermalBC` — energy/heat boundary conditions;
- `AbstractMassBC` — species and pressure-driven (Darcy) boundary conditions.

Concrete types live in the `BoundaryConditions` module
(`src/BoundaryConditions/boundary_conditions.jl`). The thermal flux laws with their
$T_s$ derivatives (`evaluate_bc_flux_and_derivative`) and the surface-temperature
Newton solve (`solve_boundary_temperature!`) live in
`src/BoundaryConditions/surface_balance.jl`, the mass-condition flux laws
(`mass_bc_flux`) and the diagnostic heat-flux decomposition
(`compute_boundary_heat_fluxes`) in `src/BoundaryConditions/boundary_fluxes.jl`, the
diagnostic record `BoundaryState` in `src/BoundaryConditions/boundary_state.jl`, and
the blowing models in `src/BoundaryConditions/blowing.jl`. The whole set is stored
in a type-stable, tag-keyed `BoundaryConditionSet` (§12.7).

Two kernels apply the conditions in the residual. `ThermalBoundaryKernel`
(`src/Kernels/thermal_boundary.jl`) applies the thermal conditions;
`GasTransportKernel` (`src/Kernels/gas_transport.jl`) applies the mass and pressure
conditions as the boundary part of the gas transport. Each holds one
`BoundaryPatch` per boundary tag of the mesh — the tag, the faces it names, and the
condition the set assigns to it — in a tuple that the kernel peels recursively, so
every condition is a concrete type and the residual specializes statically on it.

> **Implementation note.** `setup` builds both kernels from the problem's
> `bc_set`. A tag in the set that the mesh does not have is refused with an
> `ArgumentError` ("the boundary conditions name the tag :…, which the mesh does
> not have (mesh tags: …)"), and `validate_problem` (run by `solve` unless
> `validate = false`) requires the tags `:top` and `:bottom` among the mesh's and a
> thermal condition for every tag of the mesh. The mass conditions are read only
> when the solve transports gas: `setup` runs no gas transport for a material
> without a gas component, or when nothing can move the gas (no component carries a
> gas-transfer coefficient, no condensed component is permeable and every mass
> condition is `ImpermeableBC`). The kernels run in the order reactions, gas
> transport, in-depth radiation, conduction, thermal boundary, relative fluxes
> (§7.6), so the surface balance sees the gas-transport fluxes of its face.

### 12.1.3 The boundary-face sign convention

Every face of the mesh has an owner and a neighbour, and its normal points from
the owner to the neighbour; a boundary face has no neighbour, so its normal points
out of the domain (`src/Meshes/topology.jl`). On a line, cell $i$ lies between
nodes $i$ and $i+1$, face $f$ sits at node $f$, and the generator tags the two end
faces (`line_mesh`, `src/Meshes/generators.jl`):

| Tag | Position | Face | Its cell | Outward normal |
|---|---|---|---|---|
| `:top` | $z=L$, exposed surface | last node ($n+1$) | $n$ | $+\hat z$ |
| `:bottom` | $z=0$, substrate | first node ($1$) | $1$ | $-\hat z$ |

Every boundary flux — energy and mass — is stored per unit area along the outward
normal, so it is positive when the quantity leaves the domain, at both ends; the
row assembly removes $F_fA_f$ from the face's cell like any other owner. A condition
whose law gives the flux *into* the domain, $q_{\text{in}}$, therefore stores
$-q_{\text{in}}$, at the top and at the bottom alike; there is no per-boundary sign
flip and no integer boundary identifier. Fluxes in this chapter are written as the
condition laws define them (positive into the domain for the thermal laws and for
`mass_bc_flux`) and converted to the outward flux where they are stored.

A patch applies on the **live face** of its tag (`Meshes.live_face`): the boundary
face itself while its cell is live, and after cells beyond it are deactivated
(Chapter 14), the face of the outermost live cell — an interior face whose normal
points out of the live cells, so the same outward convention holds there.

---

## 12.2 The surface energy balance and the Newton solve

The surface temperature $T_s$ is the keystone of the thermal boundary. It is the
temperature of the exposed face itself, distinct from the adjacent cell-center
temperature $T_{\text{cell}}$, and it is determined by an instantaneous balance
between the externally applied flux and the heat conducted into the slab. Every
thermal condition goes through this balance except a bare `AdiabaticBC`,
`HeatFluxBC` or `DirichletThermalBC`, whose flux and face temperature follow in
closed form (§12.3.1, §12.3.2, §12.3.6); any combination of conditions, including a
combination of those three, goes through it.

### 12.2.1 Continuous statement

Consider the exposed face. An energy balance on an infinitesimally thin control
volume at the surface (zero thermal mass) equates the net flux applied by the
environment to the flux conducted into the material:

$$
q_{\text{BC}}(T_s) \;=\; q_{\text{cond}}(T_s, T_{\text{cell}}),
$$

with, for the most general thermal BC,

$$
q_{\text{BC}}(T_s) = a_{\text{eff}}\,q_{\text{inc}}\;-\;\varepsilon_{\text{eff}}\sigma \bigl(T_s^4-T_\infty^4\bigr)
\;-\;h_{\text{conv}}\,(T_s-T_\infty)
$$

(absorbed incident flux, re-radiation exchange, convective exchange; all signed
positive into the domain) and the conducted flux written as a one-sided Fourier
difference from the surface to the first cell center,

$$
q_{\text{cond}}(T_s,T_{\text{cell}})
= k\,\frac{T_s - T_{\text{cell}}}{d},
$$

where $d$ is the distance from the boundary face to the center of its cell and $k$
is the effective conductivity of the boundary cell (§5, Effective Properties &
Mixing Rules). Both sides are positive when heat flows into the material.

The balance carries no gas-enthalpy term. The gas leaving through the face takes
its sensible enthalpy at the face temperature with the face's energy flux,
$\sum_j J_j\,h_j(T_s)$ along the outward normal, which the thermal boundary kernel
adds once the balance has fixed $T_s$ (§7.4.2): under local thermal equilibrium the
gas crosses the surface at $T_s$, so a surface term
$\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$ on top of it would count that enthalpy
twice (§1.4.3, §7.4.3). The outward gas flux enters the balance only through the
blowing model of a `ConvectiveBC`, which scales the convective coefficient
(§12.3.3).

The conductive distance is the distance from the cell center to the face along the
face normal, computed from the state-dependent geometry (so automatic
differentiation and a moving mesh both reach it),

$$
d = \texttt{center\_to\_face\_distance}(\text{face}, \text{cell}),
$$

which on a line is exactly half the boundary cell's thickness (the cell center is
the midpoint of its two nodes in every coordinate system). It is smoothly floored at
$10^{-8}$ of its reference value, so a collapsing trial cell keeps it positive, and
$k$ is the cell conductivity smoothly floored at $10^{-6}$ W·m⁻¹·K⁻¹. Using the
half-width places the surface at the physical face and the datum $T_{\text{cell}}$
at the cell center, which is the one-sided gradient of a cell-centered scheme.

### 12.2.2 Residual and Newton iteration at the trial clamp

The condition is evaluated at the **trial clamp** of the surface temperature,
$\tilde T_s=\mathrm{trial}(T_s)$: $T_s$ smoothly clamped into the property bracket
$[200,\,3000]$ K with 1 K joins, the identity inside (`trial_temperature`,
`src/States/guards.jl`). Define the surface residual

$$
f(T_s) \;=\; q_{\text{BC}}\big(\tilde T_s\big)\;-\;q_{\text{cond}}(T_s).
$$

The root $f(T_s)=0$ is found by Newton iteration

$$
T_s^{(k+1)} = T_s^{(k)} - \frac{f(T_s^{(k)})}{f'(T_s^{(k)})},
\qquad
f'(T_s) = \frac{dq_{\text{BC}}}{dT}\Big|_{\tilde T_s}\,\mathrm{trial}'(T_s) - \frac{k}{d}
$$

(`solve_boundary_temperature!`, `surface_balance.jl`). The iterate itself is
never clamped. Inside the bracket the trial clamp is the identity and the balance is
the one of §12.2.1; outside it the condition is constant, the balance reduces to
conduction against a constant applied flux, and $f$ has the root
$T_s = T_{\text{cell}} + q_{\text{BC}}(T_{\text{bound}})\,d/k$. So the balance has a
root for any cell temperature it receives, with finite derivatives; when the root
lies beyond the bracket the face follows the cell one for one
($\partial T_s/\partial T_{\text{cell}}=1$). A test drives the surface with an
incident flux of 10 MW·m⁻² that the hottest admissible surface cannot re-radiate and
checks exactly this root and derivative. In the residual the cell temperature passed
to the balance is itself the cell's trial temperature, so a trial state of the
integrator with a cell temperature far outside the bracket (a second test uses
$\pm10^5$ K) gives a finite face flux and neither clamps nor warns.

The initial guess is the adjacent cell temperature clamped to the property
bracket, $T_s^{(0)} = \min\bigl(\max(T_{\text{cell}}, 200\ \text{K}),
3000\ \text{K}\bigr)$, taken on the plain value, so the solve is a pure function of
its arguments: it never starts from a stored state, and a residual evaluation never
depends on the ones before it. The cell temperature it receives is the cell's trial
temperature, the conductivity and emissivity are the cell's (§12.2.6), and the
outward gas fluxes are the face's gas-transport fluxes (§12.3.3).

**Convergence.** The residual is compared against a relative tolerance scaled by
the dominant physical flux,

$$
|f(T_s)| < \max\!\Bigl(\texttt{tol}\cdot
\max\!\bigl(|q_{\text{BC}}|,\,|q_{\text{cond}}|,\,1\bigr),\ f_{\text{floor}}\Bigr),
\qquad \texttt{tol}=10^{-10},
$$

with a maximum of 50 iterations. The floor
$f_{\text{floor}} = 32\,\epsilon_{\text{mach}}\max(|T_s|,1)\,\bigl(k/d+|dq_{\text{BC}}/dT|\,\mathrm{trial}'(T_s)\bigr)$
is the flux change from a few ulps of $T_s$: the conducted flux is a difference of
two temperatures, so a face with near-zero net flux cannot resolve its residual
below it. Scaling by the dominant flux makes the criterion meaningful whether the
surface flux is $\sim 10^4$ W·m⁻² (cone) or near zero (near-adiabatic).
**Nothing aborts and nothing warns on non-convergence**: on max-iterations or a
vanishing derivative ($|f'|<10^{-30}$) the flux is evaluated at the last iterate,
and the face's diagnostic record notes the outcome (`converged`, `iterations`;
§12.2.7, asserted by a test that limits the budget to one iteration). The
integrator's error control and domain test (§15.5.4) are the backstop that rejects
a step built on a bad surface state. A non-finite cell temperature gives a
non-finite face temperature, which is returned and not stored.

For the physical BCs, convergence is fast and robust: convection contributes a
constant negative slope $-h_{\text{conv}}$ (one Newton step is exact), and
radiation contributes a monotone cubic slope $-4F\varepsilon_{\text{eff}}\sigma
T_s^3$ that Newton resolves in a few iterations. Combined BCs sum monotone
contributions, and the trial clamp is monotone, so without blowing the residual is
strictly decreasing in $T_s$ and the root is unique.

### 12.2.3 Automatic differentiation through the solve

The solver is written to carry whatever scalar type its inputs carry. It promotes
to the widest type among $(T_{\text{cell}}, k, d, t, \varepsilon_{\text{eff}},
a_{\text{eff}})$ and the outward gas fluxes; when the colored Jacobian evaluates the
residual on dual numbers (Chapter 15) this is a `ForwardDiff.Dual`. The Newton
*convergence test* uses the plain value of the residual, so the iteration count is
identical with or without dual numbers, while the returned tuple
$(T_s,\,q_{\text{applied}},\,q_{\text{conducted}})$ is fully dual-typed and carries
the implicit-function derivatives into the Jacobian. One guard is needed: on the
very first iteration with dual numbers the residual may pass the plain tolerance
while the dual part has not yet propagated, so convergence is not accepted at the
first iteration in that case. The diagnostic `BoundaryState` stores plain `Float64`
(it is bookkeeping, never differentiated).

### 12.2.4 Per-type flux and derivative

The Newton kernel calls `evaluate_bc_flux_and_derivative(bc, T, t, ε_eff,
radiation_model, face; a_eff)` at $T=\tilde T_s$ (the `BlowingFace` at the iterate,
read only by convective BCs, §12.3.3), which returns $(q_{\text{BC}},
dq_{\text{BC}}/dT)$ for each concrete thermal BC. The pairs (from
`surface_balance.jl`) are collected below; $\varepsilon_{\text{eff}}$ is the
composition-weighted surface emissivity of the boundary cell and the keyword `a_eff`
(default $\varepsilon_{\text{eff}}$) its surface absorptivity for incident radiation
(§12.2.6).

| BC | $q_{\text{BC}}(T)$ | $dq_{\text{BC}}/dT$ |
|---|---|---|
| `AdiabaticBC` (inside a combination) | $0$ | $0$ |
| `HeatFluxBC` (inside a combination) | $q(t)$ (const or $q(t)$) | $0$ |
| `ConvectiveBC` | $h_{\text{conv}}(t)\,R_b\,(T_\infty(t)-T)$ | $-h_{\text{conv}}(t)\,R_b$ |
| `RadiativeBC` | $F\,\varepsilon\,\sigma\,(T_\infty^4-T^4)$ | $-4F\,\varepsilon\,\sigma\,T^3$ |
| `RadiativeFluxBC` (`NO_RADIATION`, `SURFACE_ABSORPTION`) | $a_s\,q_{\text{inc}}(t)$ | $0$ |
| `RadiativeFluxBC` (`BEER_LAMBERT`) | $0$ | $0$ |
| `DirichletThermalBC` (penalty, inside a combination) | $h_{\text{large}}\,(T_{\text{bc}}-T)$ | $-h_{\text{large}}$ |
| `CombinedThermalBC` | $\sum_m q_{\text{BC},m}$ | $\sum_m dq_{\text{BC},m}/dT$ |

Here $R_b = h_{\text{eff}}/h$ is the blowing ratio of the convective condition's
model ($R_b=1$ without outflow or with `NoBlowing()`), and $a_s$ is the
condition's own absorptivity or $a_{\text{eff}}$. Each BC is derived and motivated
in §12.3.

**Neglected emissivity slope.** For a `RadiativeBC` whose own emissivity is a
function of the surface temperature, $\varepsilon=\varepsilon(T_s)$, the exact
radiative derivative would be

$$
\frac{dq_{\text{BC}}}{dT_s}
= -4F\varepsilon(T_s)\sigma T_s^3
  \;+\; F\,\frac{d\varepsilon}{dT_s}\,\sigma\,(T_\infty^4 - T_s^4),
$$

but the code keeps only the first term. The
omitted term is part of the Newton *derivative* of the iteration, not the residual
itself — the residual still uses the true $\varepsilon(T_s)$ — so the converged
$T_s$ is unbiased; only the convergence rate degrades from quadratic to
super-linear when $d\varepsilon/dT$ is large. This is documented as a known
limitation (§12.12). The Jacobian of the system is not affected: it differentiates
the converged solve on dual numbers.

### 12.2.5 Application into the flux array

After the solve, the thermal boundary kernel writes the **conducted** flux
$q_{\text{conducted}}$ (not the applied flux) into the energy face flux, along the
outward normal, together with the enthalpy of the outward gas fluxes at the face
temperature, and records the face temperature:

$$
F^{E}_f \mathrel{+}= -\,q_{\text{conducted}} + \sum_j J_{j,f}\,h_j(T_s),
\qquad
T_f = T_s
$$

(`ThermalBoundaryKernel`, `src/Kernels/thermal_boundary.jl`; the enthalpy term only
when gas transport runs). Writing the conducted flux (rather than the applied flux)
is what makes the scheme energy-conserving: the row assembly then removes from the
boundary cell the flux that actually enters it, and the surface energy balance
$q_{\text{BC}}=q_{\text{cond}}$ guarantees that the energy delivered at the surface
equals the energy conducted inward. The saved solution reports both, and their
difference is at round-off when the balance has converged (7×10⁻¹² W·m⁻² in the
worked example of §12.10).

No other kernel writes a boundary flux: the conduction kernel visits only the
interior faces between two live cells (§7.1.3), and the boundary faces — and, after
deactivation, the live face of each tag — belong to the thermal boundary kernel, so
there is no double-counting.

### 12.2.6 Effective surface emissivity and absorptivity (composition-weighted)

When a `RadiativeBC` is constructed with $\varepsilon=\texttt{nothing}$, the
re-radiation uses the boundary cell's effective emissivity
$\varepsilon_{\text{eff}}$, the mean over the components weighted by the mixing
(matrix) volume fractions, each component emissivity evaluated at the cell's trial
temperature,

$$
\varepsilon_{\text{eff}} = \sum_j v_j\,\varepsilon_j(T_{\text{cell}}),
$$

a convex combination of values in $[0,1]$, so it stays in $[0,1]$ without a clamp
(§5.8.1). When a `RadiativeFluxBC` is constructed with
`absorptivity = nothing`, the absorbed incident flux uses the boundary cell's
effective surface absorptivity

$$
a_{\text{eff}} = \sum_j v_j\,a_j,\qquad
a_j = \texttt{surface\_absorptivity}_j \text{ if set, else } \varepsilon_j(T_{\text{cell}}),
$$

(§5.8.2). Two physical principles are encoded here. First, **Kirchhoff's law**
$a=\varepsilon$ at thermal equilibrium is the default: with no component
absorptivity set, $a_{\text{eff}}\equiv\varepsilon_{\text{eff}}$ and the same
quantity serves as both emissivity (re-radiation) and absorptivity (absorbed
incident flux). A component-level `absorptivity` relaxes this where the heater
spectrum differs from the surface's own emission spectrum, keeping the two
mixtures on the same volume weights. Second, surface emissivity is a
*geometric* property of the exposed mixture, so a volume-weighted (not
mass-weighted) average is appropriate and it **automatically tracks composition
changes** during pyrolysis (e.g. polymer $\to$ char as $v_j$ shifts). The state
mapping evaluates both per cell at the positive-part concentrations and the trial
temperature (`ds.ε[i]`, `ds.a[i]`, `src/States/derived.jl`), and the thermal boundary
kernel passes them to `solve_boundary_temperature!` as the positional
$\varepsilon_{\text{eff}}$ and the keyword `a_eff`. They are evaluated at the *cell*
temperature, so they are constants of the surface Newton iteration; the colored
Jacobian evaluates them on dual numbers, so
$\partial\varepsilon_{\text{eff}}/\partial T_{\text{cell}}$ and
$\partial a_{\text{eff}}/\partial T_{\text{cell}}$ from a temperature-dependent
component emissivity enter the energy rows exactly.

### 12.2.7 `BoundaryState` storage

The thermal boundary kernel holds one diagnostic record per face
(`kernel.states[f]`, `BoundaryState`, `src/BoundaryConditions/boundary_state.jl`):
the solved surface temperature, the applied flux, the conducted flux, the inputs
$(T_{\text{cell}}, k, d)$, convergence metadata (`converged`, `iterations`), and the
blowing diagnostics at the converged $T_s$ (`blowing_mcp`, `blowing_ratio`, and the
pooling-model fields `pooling_factor`, `face_fraction`, `face_buoyancy`; §12.3.3).
The fields are plain `Float64` (`converged` a `Bool`, `iterations` an `Int`); the
surface balance writes them as output only and never reads them back, so the
residual carries no state between evaluations (a test evaluates a wild trial state —
a cell temperature of $-9.1\times10^{40}$ K, or `NaN` — and then the original state,
and recovers the residual of a fresh setup to $10^{-12}$, on a fixed and a moving
mesh in both energy forms; another starts the solve from a record poisoned with
`NaN` and gets the fresh result). A face whose condition is a
closed form does not run the balance and does not write its record. The live
(possibly dual-typed) quantities are returned from the solve and flow into the
residual; the solution's per-tag blowing ratio reads the records after evaluating the
residual at each saved state.

---

## 12.3 Thermal boundary conditions

This section derives each concrete thermal BC, its physical meaning, and how it
plugs into §12.2. All flux laws are signed **positive into the domain**; the kernel
stores their negative along the outward normal (§12.1.3).

### 12.3.1 `AdiabaticBC`

Perfect insulation: $q_{\text{BC}}=0$, $dq_{\text{BC}}/dT_s=0$. The surface energy
balance reduces to $q_{\text{cond}}=0$, hence $T_s=T_{\text{cell}}$ — the surface
floats to the adjacent cell temperature with zero conduction across the face.
A bare `AdiabaticBC` is applied in closed form, without the Newton solve: the
outward energy flux is zero and the face temperature is the cell temperature.
Typical for the substrate face in a cone-calorimeter idealization with a
well-insulated back. (`boundary_energy_flux`, `src/Kernels/thermal_boundary.jl`.)

### 12.3.2 `HeatFluxBC`

A prescribed flux $q_{\text{BC}}=q(t)$, constant (`HeatFluxBC{<:Real}`) or a
function of time (`HeatFluxBC{<:Function}`), independent of $T_s$ so
$dq_{\text{BC}}/dT_s=0$. Sign: positive $q$ heats the material. A bare `HeatFluxBC`
is applied in closed form: the balance $q(t) = k\,(T_s-T_{\text{cell}})/d$ is linear,
so the outward energy flux is $-q(t)$ and the face temperature is

$$
T_s = T_{\text{cell}} + \frac{q(t)\,d}{k}.
$$

(On a 1 cm slab of 10 cells with $k=0.2$ W·m⁻¹·K⁻¹ and the surface cell at 490 K,
`HeatFluxBC(3e4)` gives $T_s = 565$ K $= 490 + 3\times10^4\cdot 5\times10^{-4}/0.2$
and an outward flux of $-3\times10^4$ W·m⁻².) Inside a `CombinedThermalBC` it enters
the Newton solve with the pair $(q(t), 0)$.

> **$T_s$ is a half-cell extrapolation.** For a pure prescribed-flux boundary
> the face temperature is the linear extrapolation
> $T_s = T_{\text{cell}} + q''\,d/k$, which assumes
> the *entire* prescribed flux crosses the boundary half-cell by conduction.
> When a volumetric consumer sits inside that half-cell (a reaction zone under
> MW-scale flux, in-depth absorption), the reported $T_s$ is biased by up to
> tens of K at coarse resolution — a real but $n$-shrinking effect. It enters
> the Landau traveling-wave identity of V-CPL-02 through
> $v_{ss}(\text{measured}\ T_s)$: with the reaction zone under-resolved
> ($\Delta z \approx 2\delta_r$) the identity misses by ≈1.3 %, falling to
> ≈7×10⁻⁴ at $\Delta z \approx \delta_r$. (A *mesh-independent* residual
> of about 0.9 % has a different source, the wave-development transient of
> a uniform-$T_0$ start; the case measures after at least six development time
> constants to exclude it.) Interpret surface-temperature outputs under `HeatFluxBC`
> accordingly, and refine the near-surface mesh when $T_s$ itself is a target
> quantity.

---

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
resolved in a single step. (`surface_balance.jl`.)

**Blowing.** Every `ConvectiveBC` carries a blowing model (`bc.blowing`, an
`AbstractBlowingModel`, default `NoBlowing()`), and the surface balance applies
it at every Newton iterate to the face state: the **signed outward** gas mass
flux of every component, read from the face's gas-transport fluxes (`c.F_mass[:, f]`
oriented out of the boundary cell: the film law of the mass condition plus, in
Darcy–Fick mode, the pressure-condition venting; §12.4, §12.5), the trial clamp of
the surface-temperature iterate and the material (`BlowingFace`,
`src/BoundaryConditions/blowing.jl`). The convective flux is
$h\,R_b\,(T_\infty - T_s)$ with $R_b = h_{\text{eff}}/h$ the blowing ratio of the
model (`blowing_state`; the code calls it `χ`, a letter the nomenclature reserves for
pyrolysis progress). It modifies only the convective *coefficient*, not the gas
enthalpy accounting, so it composes exactly with the gas enthalpy flux of the face
(§7.4.3) — no double-count is possible. Zero outflow gives $R_b = 1$ exactly for
every shipped model, so a face without outflow, such as a sealed back face, is
unaffected automatically. The gas-transport kernel runs before the thermal
boundary kernel and the relative-flux kernel after it, so the blowing models see the
gas-transport fluxes alone.

*Film law with a fitted constant* (`CouetteFilm(λ)`):

$$
R_b \;=\; \frac{B_\lambda}{e^{B_\lambda}-1},
\qquad
B_\lambda \;=\; \frac{2\lambda\,\dot m''\,c_{p,g}}{h},
\qquad
\dot m''\,c_{p,g} \;=\; \sum_{g}\dot m_g\,c_{p,g}(T_s).
$$

$\lambda = 1/2$ is the classical Couette-film (Spalding) factor, $B_\lambda = B =
\dot m'' c_{p,g}/h$; a smaller $\lambda$ fitted to a resolved flow weakens the
reduction at a given outflow and holds for that gas only. Limits in $B$:
$B = 0$ recovers the bare $h$ exactly; $B = 4$ (peak-PMMA outgassing against
canonical cone $h$) gives $h_{\text{eff}}/h = 0.075$; $B < 0$ (net inflow) is
the smooth suction enhancement.

*Property-ratio film for any light gas* (`PropertyRatioFilm(free_stream; n_Γ)`),
anchored on the Couette film of the free stream:

$$
B_\Gamma \;=\; \frac{\dot m''\,c_{p,\infty}(T_s)}{h}\left(\frac{\Gamma_\infty}{\Gamma_g}\right)^{n_\Gamma},
\qquad \Gamma = \frac{k}{c_p},
\qquad R_b = \frac{B_\Gamma}{e^{B_\Gamma}-1},
$$

with $\Gamma_\infty$ of the free stream and $\Gamma_g$ of the mixture leaving the
face (heat capacity weighted by mass fraction, conductivity by Wilke's rule
over the components with outward flux, each species' viscosity from Eucken's
relation when the component carries none). $n_\Gamma = 0.6$ is Kays's turbulent
exponent on the heat-capacity ratio, reproduced by the chamber CFD; with
$n_\Gamma = 0$ this is the Couette film in the free-stream heat capacity. Under
net inflow the arriving gas is the free stream and the ratio is 1.

*Pooling* (`PoolingFilm(free_stream; n_Γ, n_D, ψ0, Λc, Λw)`): the same
$B_\Gamma$ with a pooling factor $\psi$, the share of the no-outflow exchange the
face keeps,

$$
R_b \;=\; \frac{B_\Gamma}{e^{B_\Gamma/\psi}-1} \;=\; \psi\,F_b(B_\Gamma/\psi),
\qquad
\psi \;=\; \psi_0 + \frac{1-\psi_0}{1+\exp[-(\Lambda-\Lambda_c)/\Lambda_w]},
\qquad
\Lambda \;=\; \frac{1-(M_{\text{face}}/T_s)(T_\infty/M_\infty)}{1-T_\infty/T_s},
$$

where $F_b(x)=x/(e^x-1)$ is the film factor and $\Lambda$ the face buoyancy relative
to the hot free stream, from the face molar mass $M_{\text{face}}$ of the species
mass fractions

$$
Y_k \;=\; \frac{\dot m''_k}{\dot m''}\left[1-\exp\!\left(-\frac{\dot m''}{\psi\,\mathcal G\,s_k^{\,n_D}}\right)\right],
\qquad \mathcal G = \frac{h}{\bar c_{p,\infty}},
\qquad s_k = \frac{D_{k,\infty}}{\alpha_{\text{th},\infty}},
$$

with $\dot m''_k/\dot m''$ the species' share of the outflow, $\mathcal G$ the
gas-exchange conductance without outflow (the film analogy: $h$ over the free-stream
heat capacity averaged between $T_\infty$ and $T_s$, a 5-point Gauss–Legendre mean),
$D_{k,\infty}$ the binary diffusivity into the free stream (Fuller, from the
component's diffusion volume on the model, else the Chapman–Enskog estimate
from its molar mass; a lumped `volatiles` gas states its own),
$\alpha_{\text{th},\infty}$ the thermal diffusivity of the free stream, and $n_D = 1/2$
the penetration-theory exponent. Components whose name matches the free
stream count in $\dot m''$ but not in the face composition. $\psi$ and the
$Y_k$ are solved together: $G(\psi) = \Psi(\Lambda(Y(\psi))) - \psi$ has
$G(\psi_0) > 0$ and $G(1) < 0$, so a 16-interval scan of $[\psi_0, 1]$ on the
plain values brackets every root, each refined by safeguarded Newton; where
two states are self-consistent the mean of the lowest and highest is taken
(the resolved flow alternates between them). One Newton step in the working
scalar type at the plain root then gives the implicit-function derivative,
so `ForwardDiff` partials ride through the solve. A face no hotter than the
free stream keeps $\psi = 1$. The three switch constants belong to the
apparatus flow, as $h$ does.

*User rule* (`FunctionBlowing(f)`): $R_b = f(\text{face}, h, T_\infty)$.

The film factor is evaluated by the numerically stable `blowing_factor` (Taylor
branch at $|x|<10^{-5}$, overflow-safe $x\,e^{-x}/(1-e^{-x})$ form for
large $x$). The Newton derivative uses the frozen-ratio form $-h\,R_b$ —
first-order consistent, so convergence is marginally sub-quadratic and the
converged root is exact; the Jacobian of the system, which differentiates the
converged solve on dual numbers, carries the full derivative. Each convective member
of a `CombinedThermalBC` applies its own model against its own $h$ and $T_\infty$;
all other BC types are untouched. At the converged $T_s$ the `BoundaryState`
records the blowing conductance `blowing_mcp`, the applied ratio `blowing_ratio`
and, for the pooling model, `pooling_factor`, `face_fraction` and `face_buoyancy`
(the first convective member's), and the solution carries the ratio and the pooling
fields per output time (`sol.extras.surface_blowing_ratio`, `surface_pooling_factor`,
`surface_face_fraction`, `surface_face_buoyancy`; per tag,
`sol.boundary_blowing_ratio`). The colored Jacobian evaluates the face's outward
mass fluxes and the blowing model on dual numbers, so the
$\partial(hR_b)/\partial(T,\xi)$ chain through the boundary cell's state is exact.
Giving a convective BC a model on a case calibrated with the bare coefficient will
raise $T_s$ and MLR; introduce it together with a flame heat-flux feedback, with
which it partially cancels (§1.8.1).
(`src/Materials/blowing_models.jl` types and gas data,
`src/BoundaryConditions/blowing.jl` evaluation.)

### 12.3.4 `RadiativeBC` (gray surface exchange with surroundings)

$$
q_{\text{BC}} = F\,\varepsilon\,\sigma\,\bigl(T_\infty^4 - T_s^4\bigr),
\qquad
\frac{dq_{\text{BC}}}{dT_s} = -4\,F\,\varepsilon\,\sigma\,T_s^3,
$$

with the Stefan–Boltzmann constant $\sigma=5.670374419\times10^{-8}$ W·m⁻²·K⁻⁴
(`src/Materials/constants.jl`, exported as `STEFAN_BOLTZMANN`). The fields are: the
**view factor** $F\in[0,1]$ (fraction of the surface's hemisphere that exchanges
with surroundings at $T_\infty$; default $1$), the **emissivity** $\varepsilon$
(either a supplied constant or function of the surface temperature, evaluated at the
trial clamp of the iterate, or, when `nothing`, the composition-weighted
$\varepsilon_{\text{eff}}$ of §12.2.6), and the surroundings temperature
$T_\infty$. This term is the surface's own re-radiation to and absorption from a
gray environment; it is the dominant cooling mechanism for a hot pyrolyzing
surface. Its quartic temperature dependence is why the face temperature, not the
cell-center temperature, must enter it (§12.1.1). The condition stays at the
surface under every radiation model. (`surface_balance.jl`.)

The neglected $d\varepsilon/dT$ term (§12.2.4) applies here when $\varepsilon$ is
a function of the surface temperature.

### 12.3.5 `RadiativeFluxBC` (external incident radiation; model-dependent)

This BC represents incident radiation from an *external* source (cone heater,
furnace, radiant panel) — not the surface's own emission. The absorbed fraction is
$a_s\,q_{\text{inc}}$, where $q_{\text{inc}}=$ `flux` (constant or $q(t)$) and
the absorptivity $a_s$ (keyword `absorptivity`) is either supplied (a constant or a function of $t$) or,
when `nothing`, the mixture surface absorptivity $a_{\text{eff}}$, which equals
$\varepsilon_{\text{eff}}$ by Kirchhoff's law unless a component sets its own
`absorptivity` (§12.2.6).

Its behavior is **dispatched on the `radiation_model`** passed to `solve` or
`setup` (default `NO_RADIATION`; Chapter 8). The crucial design point is that
incident radiation must be deposited exactly once:

- **`NO_RADIATION` and `SURFACE_ABSORPTION`** treat the radiation identically, as a
  *surface* heat flux (a test asserts that the two give identical residuals). The
  Newton solve receives
  $$q_{\text{BC}} = a_s\,q_{\text{inc}}(t), \qquad dq_{\text{BC}}/dT_s = 0,$$
  so the absorbed flux enters the surface energy balance and raises $T_s$ above
  $T_{\text{cell}}$ when there is net heating. A direct call of
  `solve_boundary_temperature!` defaults its `radiation_model` keyword to
  `SURFACE_ABSORPTION`, which is the same treatment.
- **`BEER_LAMBERT`**: the radiation is absorbed *in depth* by `BeerLambertKernel`,
  through face intensities that are algebraic unknowns of the state, the exposed
  face's intensity being $a_s\,q_{\text{inc}}(t)$ with $a_s$ the condition's own
  absorptivity or that of the outermost live cell (§8.3, §8.3.4). To avoid
  double-counting, the condition contributes nothing to the surface balance,
  $$q_{\text{BC}} = 0, \qquad dq_{\text{BC}}/dT_s = 0 .$$

Under `BEER_LAMBERT` the beam is fed by the **first** `RadiativeFluxBC` member of
the `:top` condition (`get_radiative_flux_bc`, which returns the first
`RadiativeFluxBC` it finds in a `CombinedThermalBC`), while every member contributes
zero to the surface balance: a second radiative-flux member is dropped silently. On
a slab of emissivity 0.9, `RadiativeFluxBC(flux = 30e3) + RadiativeFluxBC(flux =
20e3)` gives an exposed-face intensity of 27 000 W·m⁻² ($0.9\times30$ kW·m⁻²) under
`BEER_LAMBERT` and a conducted surface flux of 45 000 W·m⁻² ($0.9\times50$ kW·m⁻²)
under `SURFACE_ABSORPTION`. Combine the incident fluxes into one condition when the
model is `BEER_LAMBERT` (§12.12). `setup` refuses a `RadiativeFluxBC` on any other tag
under `BEER_LAMBERT` (in-depth radiation enters through the `:top` faces only), and
refuses `BEER_LAMBERT` on cylindrical or spherical coordinates.

A note on consistency. The surface balance applies $a_s q_{\text{inc}}$ as a
boundary flux under `SURFACE_ABSORPTION`, and no volumetric source exists under that
model, so the energy enters once, through the face; only `BEER_LAMBERT` routes it
through the volume, and only then is the boundary flux suppressed. The diagnostic
decomposition `compute_boundary_heat_fluxes` follows the *same* dispatch: its `net`
is the flux the balance actually applies — $a_s q_{\text{inc}}$ is inside `net`
under `NO_RADIATION`/`SURFACE_ABSORPTION` and absent under `BEER_LAMBERT`, where the
in-depth deposit is reported separately as `extras.volumetric_absorption`
(§12.9.1).

> **Implementation note.** `BeerLambertKernel(topology, geometry, bc_set)`
> (`src/Kernels/radiation.jl`) takes its source from
> `get_radiative_flux_bc(bc_set.thermal[:top])` (a zero flux when the condition has
> no radiative-flux member) and evaluates the incident intensity at the time of each
> residual evaluation (`incident_intensity`, §8.3.4). The surface balance receives
> the radiation model from the thermal boundary kernel, which `setup` builds with
> the solve's `radiation_model`.

### 12.3.6 `DirichletThermalBC` (prescribed surface temperature)

A prescribed surface temperature $T=T_{\text{bc}}(t)$ is enforced by one of two
mechanisms, selected by how the BC appears, rather than by an algebraic
constraint:

1. **Bare, in closed form** (`boundary_energy_flux`,
   `src/Kernels/thermal_boundary.jl`): a bare `DirichletThermalBC` bypasses the
   Newton solve. The surface temperature is the prescribed value,
   $T_s = T_{\text{bc}}(t)$, and the conducted flux follows directly,
   $$q_{\text{cond}} = k\,\frac{T_{\text{bc}} - T_{\text{cell}}}{d},
   \qquad F^E_f = -q_{\text{cond}} = k\,\frac{T_{\text{cell}} - T_{\text{bc}}}{d}$$
   along the outward normal, with the same distance $d$ as the Newton path. (A face
   held at 600 K over a 490 K cell, $k = 0.2$ W·m⁻¹·K⁻¹, $d = 0.5$ mm, conducts
   $44\,000$ W·m⁻² inward: outward flux $-44\,000$ W·m⁻².) The saved solution reports
   $T_{\text{bc}}$ as the surface temperature and the conducted flux as the net
   surface flux.
2. **Nested inside a `CombinedThermalBC`**, in the Newton solve
   (`surface_balance.jl`): a penalty "flux" with a very large conductance,
   $$q_{\text{BC}} = h_{\text{large}}\,(T_{\text{bc}} - T_s),\qquad
   dq_{\text{BC}}/dT_s = -h_{\text{large}}, \qquad h_{\text{large}}=10^{8}
   \;\text{W·m}^{-2}\text{·K}^{-1}.$$
   The large slope drives $T_s\to T_{\text{bc}}$ in one step, with the face missing
   the prescribed value by $q_{\text{cond}}/h_{\text{large}}$: in the example above,
   with a convective member of zero coefficient added, $T_s = 599.99956$ K, a miss of
   $4.4\times10^{-4}$ K at 44 kW·m⁻². The value $10^8$ pins $T_s$ to $T_{\text{bc}}$
   within that margin while keeping the iteration well conditioned (a still-larger
   penalty would degrade conditioning with no accuracy benefit).

The penalty is less clean than a direct constraint; it is used only where the
prescribed temperature is combined with other conditions, and it is robust in
practice because the prescribed temperature pins the boundary tightly.

### 12.3.7 `CombinedThermalBC` (additive superposition)

Real surfaces experience several mechanisms at once. The `+` operator builds an
additive composite, e.g. the canonical cone-calorimeter top BC

```julia
bc = RadiativeFluxBC(flux = 50e3) +              # absorbed incident cone flux
     RadiativeBC(ε = nothing, T_inf = 300.0) +   # surface re-radiation loss
     ConvectiveBC(h = 10.0, T_inf = 300.0)       # convective loss
```

Fluxes and derivatives sum termwise,

$$
q_{\text{BC}} = \sum_m q_{\text{BC},m}(T_s,t),\qquad
\frac{dq_{\text{BC}}}{dT_s} = \sum_m \frac{dq_{\text{BC},m}}{dT_s}.
$$

The sum is built by a `@generated` function that unrolls the tuple at compile time
for type stability (`surface_balance.jl`). Because the radiation model can affect a
`RadiativeFluxBC` member, the combined evaluator threads `radiation_model` through
to every child, and each convective child applies its own blowing model. The `+`
operator overloads handle every nesting case (BC+BC, combined+BC, BC+combined,
combined+combined) so associativity holds (`boundary_conditions.jl`). A combination
always goes through the Newton solve, even when every member is one of the
closed-form conditions of §12.3.1, §12.3.2 and §12.3.6.

---

## 12.4 Mass (species) boundary conditions

The species boundaries supply the gas mass flux $J_j$ at each live boundary face.
Unlike the thermal boundary, there is no surface-temperature nonlinearity, so the
flux is evaluated directly from the state of the boundary cell. The internal
diffusive flux model is the volume-fraction (partial-pressure) gradient form of §9,
Gas-Phase Mass Transport & Pressure; here it is closed at the surface by an
external resistance. The conditions act only when the solve transports gas
(§12.1.2), and only on the gas components; solids and liquids never cross a boundary
face (on a moving mesh the correction vanishes at the end nodes, §11.7.5).

### 12.4.1 `ImpermeableBC` (default)

$J_j = 0$ for all components — no gas crosses the boundary. This is the default
for both faces (`BoundaryConditionSet`) and is the physically correct choice for
the substrate (sealed back face) and any non-venting surface. It is also the
implicit closure for solid/liquid components everywhere, since only gas components
are transported.

### 12.4.2 `MassFluxBC`

A prescribed species flux $J_j = J(t)$ (constant or function), positive into the
domain (`boundary_fluxes.jl`), applied to every gas component. Used to impose a
known gas influx/efflux, typically zero or negative (gas escaping).

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
*every* gas species into the domain. A test lets an ambient gas
(`ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.5))`) enter a gas-free column: the
outward flux of the top face is $-0.025$ kg·m⁻²·s⁻¹, the surface cell gains exactly
$A\cdot0.025$ kg/s, and the saved boundary mass flux reports the inflow on a fixed and
on a moving mesh.

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

where $d$ is the cell-center-to-face distance and $\lambda$ is the
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

> **Implementation note.** The flux law is `mass_bc_flux(bc::DiffusiveMassBC,
> ξ_surface, λρ, ρ, d, t, j)` (`boundary_fluxes.jl`) — every mass-BC flux method carries
> the trailing species index `j` so per-species ambients resolve correctly. The
> gas-transport kernel (`_boundary_mass_flux!`, `src/Kernels/gas_transport.jl`) passes
> the boundary cell's concentration $\xi_{j,s}=m_{j,s}/V_s$, the mixture density
> $\rho$ at its positive-part concentrations and trial temperature, the *combined
> transport coefficient* $\lambda\rho$ (the cell's gas-transfer coefficient times
> $\rho$), and the center-to-face distance smoothly floored at $10^{-8}$ of its
> reference value. For $\rho \le 0$ the flux is zero, and otherwise the function
> recovers $\lambda = \lambda\rho/\rho$ by exact division (no floor — the caller
> constructs $\lambda\rho$ as the product, so division is exact for any positive
> density), because $h_m$ is defined per concentration, not per mass fraction.

**Sign reconciliation.** Three sign conventions coexist and must be reconciled
(`boundary_fluxes.jl`, `gas_transport.jl`):

| Quantity | Convention |
|---|---|
| `mass_bc_flux(...)` | **negative** = gas escaping (flux-into-domain) |
| `compute_surface_mass_flux(...)` (diagnostic) | **positive** = gas escaping |
| boundary face flux `F[j, f]`, and `sol.boundary_mass_flux[tag]` | **positive** = along the outward normal (gas escaping), at both ends |

The kernel bridges these by negating the condition's flux into the face's outward
flux, `F[j, f] += v_out ξ_j/φ − mass_bc_flux(...)`, at the top and the bottom alike:
with $J_{\text{bc}}<0$ for escaping gas, the stored flux is positive — gas leaving
along the outward normal — which the row assembly removes from the boundary cell.
The application **accumulates**: in Darcy–Fick mode the face also carries the
pressure-condition venting flux $v_{\text{out}}\xi_j/\varphi$ (§12.5), and the film
law adds to it; in Fickian mode the venting term is zero. A test vents a symmetric
column through `DiffusiveMassBC(0.05)` at both faces: the bottom flux is outward and
equal to the top one, the integral of both closes against the lost inventory to
$10^{-3}$, and the mass-loss rate counts both faces.

### 12.4.5 `CombinedMassBC`

Like the thermal composite, $J_j=\sum_m J_{j,m}$, built by a `@generated` unrolled
sum (`boundary_fluxes.jl`), e.g. `DiffusiveMassBC(Inf) +
ConvectiveMassBC(h_m = 0.01, ξ_inf = 0.0)`. Pressure BCs (§12.5) participate in the
composite by contributing zero *film* flux (they act through the outward Darcy
velocity, not the film law), so they can be combined freely with `DiffusiveMassBC`
(`boundary_fluxes.jl`); the outward Darcy velocity of a composite is the sum of its
members' (§12.5).

---

## 12.5 Pressure boundary conditions (Darcy flow)

In Darcy–Fick transport mode (§9; selected when a condensed component is
permeable), gas advection is driven by the Darcy velocity
$v_g = -(\kappa/\mu)\nabla P$, and the pressure boundaries set the velocity at the
boundary faces. These BCs are subtypes of `AbstractMassBC` so they live in the
mass slot of the `BoundaryConditionSet`; they contribute **no** film flux
(returning $0$ from `mass_bc_flux`, §12.4.5) — they act through
`outward_darcy_velocity(bc, P, mob, d)` (`src/Kernels/gas_transport.jl`), the
superficial velocity along the face's **outward** normal from the boundary cell's
pressure $P$ and mobility $\mathrm{mob}=\kappa/\mu$ (the viscosity smoothly floored at
$10^{-12}$ Pa·s, §5.9), its center lying $d$ from the face. That velocity enters the
species rows as the boundary venting flux

$$
F_{j,f} \mathrel{+}= \frac{v_{\text{out}}}{\varphi_s}\,\xi_{j,s}
\quad (v_{\text{out}} \ge 0),
$$

outflow only: the ambient beyond the face carries none of the tracked gases, so an
inward velocity draws no gas in (§9.5.3). $\varphi_s$ is the boundary cell's smoothly
floored porosity, and the outflow test reads the plain value of $v_{\text{out}}$, as
the interior upwind side does. The orientation is the same at both ends, so one
formula serves the top and the bottom.

### 12.5.1 `DirichletPressureBC` (prescribed pressure)

Prescribes the boundary pressure $P_{\text{bc}}$ (default $101325$ Pa,
atmospheric). The outward velocity follows from Darcy's law with a one-sided
pressure gradient between the interior cell and the face:

$$
v_{\text{out}} = \mathrm{mob}\,\frac{P - P_{\text{bc}}}{d},
$$

with $d$ the cell-center-to-face distance and $\mathrm{mob}$ the cell's mobility. A higher
interior pressure ($P>P_{\text{bc}}$) yields $v_{\text{out}}>0$, outflow, at either
face; an impermeable cell ($\kappa = 0$) has zero mobility and vents nothing.

### 12.5.2 `NeumannPressureBC` (prescribed gradient)

Prescribes the pressure gradient $\partial P/\partial n$ along the **outward
normal** (default $0$ = impermeable):

$$
v_{\text{out}} = -\mathrm{mob}\,\frac{\partial P}{\partial n}.
$$

The zero-gradient default makes a `NeumannPressureBC()` equivalent to an impermeable
wall for Darcy flow.

### 12.5.3 `ConvectivePressureBC` (Robin / surface resistance)

A Robin-type law for a surface with finite resistance to gas escape (char layer,
protective coating, blowing resistance):

$$
v_{\text{out}} = h_P\,\bigl(P - P_{\text{ext}}\bigr)\,\frac{\mathrm{mob}}{\mathrm{mob} + \eta_m},
$$

with $h_P$ the pressure-transfer coefficient (m·Pa⁻¹·s⁻¹), $P_{\text{ext}}$ the
exterior pressure (`ConvectivePressureBC(h_P, P_ext)`), and the gate $\mathrm{mob}/(\mathrm{mob}+\eta_m)$, $\eta_m=10^{-30}$ m²·Pa⁻¹·s⁻¹,
equal to one on any permeable cell (exactly, in floating point, for any physical
mobility) and zero on an impermeable one, smooth between — an impermeable surface cell
cannot vent regardless of $h_P$. The limits bracket the other two pressure BCs:
$h_P\to\infty$ pins the boundary pressure (Dirichlet) and $h_P\to0$ seals it
(zero-gradient Neumann / impermeable).

> **Implementation note.** All three pressure BCs are methods of
> `outward_darcy_velocity(bc, P, mobility, d)` (`src/Kernels/gas_transport.jl`); a
> `CombinedMassBC` sums the contributions of its members, and every other mass
> condition contributes $0$ (e.g. a `DiffusiveMassBC` in the mass slot does not set
> a Darcy velocity). The gas-transport kernel evaluates it per live boundary face;
> the colored Jacobian evaluates the same code on dual numbers and the solution
> assembly on the saved states, so the three cannot drift apart.

---

## 12.6 Mesh-motion boundary conditions (ALE)

On a moving mesh (`use_ale = true`; Chapter 11) the node positions $z_n$ are
differential unknowns and the node velocities algebraic ones, and the mesh-velocity
boundary treatment is **fixed by the formulation, not user-configurable**: the
mesh's pinned node — on a line the back node, $z=0$ — is at rest, $0=w_1$, and
every other node velocity follows from the balance rows of the cells,
$0 = A_{i+1}(w-c)_{i+1} - A_i(w-c)_i - (\theta_i-\theta_A)V_i$ (`BalanceRowClosure`,
§11.3.2), which on a planar line without a correction is the recursion
$w_{i+1} = w_i + (\theta_i - \theta_A)\,\Delta z_i$ (Chapter 10). The exposed surface
node is an end node, where the correction vanishes (§11.7.5), so it is always
material-following (Lagrangian) — it recedes or expands with the integrated strain.

There are no user-facing mesh-motion BC types: `BoundaryConditionSet` carries only
thermal and mass slots. What a user chooses is `use_ale` (without it every node is
pinned and the geometry never changes), the mesh itself (`line_mesh` pins its first
node; `fixed = true` pins every node), and the `lagrangian` nodes that always move
with the solid (layer interfaces; §11.7.5); `w_override` prescribes every node's
velocity for verification (§11.7.8). The balance rows need exactly one pinned node,
at an end of the line. See Chapter 11 for the full mesh-velocity derivation and the
geometric conservation law.

---

## 12.7 The `BoundaryConditionSet` container

### 12.7.1 Structure and tag keying

All boundaries are bundled into a single type-stable container
(`boundary_conditions.jl`):

```julia
struct BoundaryConditionSet{TT<:NamedTuple,MT<:NamedTuple}
    thermal::TT   # NamedTuple{tags, Tuple{<:AbstractThermalBC,...}}
    mass::MT      # NamedTuple{tags, Tuple{<:AbstractMassBC,...}}
end
```

Each boundary tag maps to exactly one thermal BC and one mass BC. A line uses
`:top` ($z=L$) and `:bottom` ($z=0$); meshes of more dimensions, at the residual
level, carry their own tags, and the same container serves them. Storing the BCs in a
typed `NamedTuple` (rather than, say, a `Dict`) means a symbol lookup compiles to a
direct field read and every condition's type is known statically.

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
neither gains energy nor loses mass. The set is checked against the mesh when the
problem is set up (§12.1.2): a tag the mesh does not have is refused, and
`validate_problem` requires `:top` and `:bottom` and a thermal condition for every
tag of the mesh, and checks the length of a per-species ambient tuple against $N_C$.

### 12.7.3 Accessors and patches

Access by tag is through `Val{tag}` dispatch:

```julia
@inline get_thermal_bc(bc_set, ::Val{tag}) where {tag} = bc_set.thermal[tag]
@inline get_mass_bc(bc_set,    ::Val{tag}) where {tag} = bc_set.mass[tag]
```

At a call site such as `get_thermal_bc(bc_set, Val(:top))`, the `Val(:top)`
resolves at compile time to a direct NamedTuple field read with zero overhead
(`boundary_conditions.jl`). The kernels hold the conditions in their patches: one
`BoundaryPatch(tag, faces, bc)` per boundary tag of the mesh, with `faces =
boundary_faces(topology, tag)`, in a tuple peeled recursively so each patch's
condition is a concrete type. `boundary_condition(kernel, tag)` returns the thermal
condition the kernel applies on a tag. The faces are those of the mesh's topology;
there is no integer boundary identifier anywhere.

### 12.7.4 Constructor parameter validation

Individual BC constructors optionally validate their physical parameters at build
time (`validate_bc_params`, `boundary_conditions.jl`): $h\ge0$, $T>0$ K,
$\varepsilon\in[0,1]$, $F\in[0,1]$, $a_s\in[0,1]$, $h_m\ge0$, $P>0$.
Function-valued parameters cannot be checked at construction and are skipped
(`_get_scalar_value` returns `nothing` for functions and property functions).

---

## 12.8 Initial conditions

### 12.8.1 The state seed

The initial state is specified by a temperature closure $T_{\text{initial}}(z)$
and per-component concentration closures $\xi_{j,\text{initial}}(z)$ on the problem,
evaluated at each cell center $z_i^c$ of the mesh's reference positions by
`initial_fields(problem, mesh)` (`src/Problem/initial_conditions.jl`):

$$
T_i(0) = T_{\text{initial}}(z_i^c), \qquad
\xi_{j,i}(0) = \xi_{j,\text{initial}}(z_i^c), \quad j=1,\dots,N_C .
$$

On a line $z_i^c$ is the scalar midpoint of the cell (on meshes of more dimensions
the centroid, an `SVector`). The mesh holds no state: `setup` (which `solve` calls)
evaluates the fields of the problem every time, forms the initial inventories
$m_{j,i}(0) = \xi_{j,i}(0)\,V_i(0)$ from the reference cell volumes and the energy
state from $T_i(0)$ (the temperature, or the inventory $H_i=\sum_j m_{j,i}h_j(T_i)$
under `energy_state = :enthalpy`), writes the history pair under a lateral law
(§12.8.3) and the node positions on a moving mesh, and then sets the algebraic
unknowns consistent with them (`consistent_state!`: the face intensities, the area
and lateral rate, the node velocities; §15.4). For example:

```julia
using Pyrolysis
using Pyrolysis.Problem: initial_fields, linear_profile

slab = SolidComponent(:slab, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, 4),
    material = Material(name = :inert, components = (slab,), reactions = ()),
    bc_set = BoundaryConditionSet(thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC())),
    T_initial = z -> linear_profile(z; z_min = 0.0, z_max = 0.01,
                                    value_min = 300.0, value_max = 400.0),
    tspan = (0.0, 1.0),
)
T0, ξ0 = initial_fields(problem, problem.mesh)
```

gives $T_0 = (312.5,\ 337.5,\ 362.5,\ 387.5)$ K at the four cell centres, and the
default concentration $\xi_0 = 1000$ kg·m⁻³ of §12.8.2 in every cell.

### 12.8.2 Auto-defaults

If the user supplies no initial conditions, the `PyrolysisProblem` constructor
synthesizes them (`pyrolysis_problem.jl`):

- $T_{\text{initial}}(z)=300$ K uniformly;
- $\xi_{1,\text{initial}}(z)=\rho_1(300\,\text{K})$ if the first component is a
  solid (the virgin solid at its reference density), and $\xi_{j,\text{initial}}=0$
  for all other components.

This "first solid at reference density, everything else zero" default models a
fresh, uniform virgin slab at room temperature — the standard cone/TGA starting
point. A single function given as `ξ_initial` is used for every component. When no
`output_times` are given, the constructor saves at
$\min(101,\max(11,\mathrm{round}((t_{\text{end}}-t_{\text{start}})/(100\,\Delta t_{\text{initial}}))))$
equally spaced times from $t_{\text{start}}$ to $t_{\text{end}}$. Profile helpers
(`step_profile`, `linear_profile`, `gaussian_profile`, `tanh_profile`, in
`Pyrolysis.Problem`, not exported at the top level) are provided for non-uniform
seeds (e.g. a pre-heated layer).

### 12.8.3 The dry-solid mass denominator $m_{\text{dry},i}$

A second quantity is fixed at initialization under a lateral shrinkage law: the
per-cell initial dry-solid mass, the second member of the conversion-history pair
and the denominator of the pyrolysis progress $\chi_i = G_i/m_{\text{dry},i}$ used by
the variable cross-section law (§10.3). `setup` writes it into the cell's history row
of the initial state, with $G_i(0)=0$:

$$
m_{\text{dry},i} = V_i(0)\,\sum_{j\in\text{solid}} \xi_{j,i}(0),
$$

summing only the solid components (the `is_solid` predicate — not liquids, not
gases) at the initial state. It is a differential unknown with no source: it changes
only when the mesh moves relative to the solid (the relative fluxes carry it with the
solid) and at deactivation events (it goes to the receiver), and its column total is
conserved (§10.3.2, §10.3.4). Like every unknown it is differentiated by the colored
Jacobian.

---

## 12.9 Boundary-flux diagnostics

### 12.9.1 Flux decomposition

For auditing and post-processing, `compute_boundary_heat_fluxes(bc, T_s, t,
ε_eff, radiation_model[, face]; a_eff)` decomposes the boundary flux into named
components (`boundary_fluxes.jl`). The optional trailing `face` is the
`BlowingFace` at the converged $T_s$ (outward species fluxes, surface
temperature, material); when supplied, the convective component is reported
with the same ratio $R_b$ its blowing model gave the surface Newton
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
`evaluate_bc_flux_and_derivative` — the applied flux — for every
radiation model, `DirichletThermalBC` excepted (its decomposition is all zero; the
saved extras report the conducted flux as its net flux), so the ledger identity
$\texttt{net} = q_{\text{cond}}$ holds at convergence (§12.9.2, §16.6). For combined
BCs the components **and** `net` are summed over children; `net` is never recomputed
from the components. Consequences of the model dispatch:

- `NO_RADIATION` / `SURFACE_ABSORPTION`: a `RadiativeFluxBC` contributes
  `absorbed = net` $= a_s q_{\text{inc}}$ ($a_s$ the explicit BC absorptivity or $a_{\text{eff}}$) — the
  flux is applied at the face.
- `BEER_LAMBERT`: a `RadiativeFluxBC` contributes zero `absorbed` and zero `net` —
  the radiation enters in depth (§8.3) and is reported separately: per unit exposed
  area, `extras.volumetric_absorption` $= I_e - I_bA_b/A_e$, the intensity entering
  the exposed face less what leaves through the last face of the beam, and
  `extras.radiation_transmitted` $= I_bA_b/A_e$ (§8.3.5 works the numbers for a
  slab).

The solution evaluates the decomposition at every saved time
(`src/Outputs/solution.jl`): at the face temperature of each live `:top` face, with
the boundary cell's emissivity and absorptivity and the face's outward mass fluxes,
area-weighted over the faces, into `sol.extras.surface_heat_flux_net`,
`surface_heat_flux_absorbed`, `surface_heat_flux_reradiation` and
`surface_heat_flux_convection`, and the conducted flux $k(T_s - T_{\text{cell}})/d$
(`conducted_flux`) into `surface_heat_flux_conduction`. The face temperatures are
`sol.surface_T` (the `:top` faces) and `sol.boundary_T[tag]` per tag (after
deactivation, of the live face), with `sol.extras.back_temperature` and
`sol.extras.back_heat_flux` for the `:bottom` condition.

### 12.9.2 Energy-balance check

The saved surface fluxes are a direct self-check of the surface energy balance
$q_{\text{BC}}=q_{\text{cond}}$: at every saved time
`surface_heat_flux_net − surface_heat_flux_conduction` is at round-off for a condition
that goes through the balance and zero by construction for the closed forms (in the
worked example of §12.10, $7\times10^{-12}$ W·m⁻² against fluxes of
$4.5\times10^4$ W·m⁻²), and
`surface_heat_flux_absorbed − surface_heat_flux_reradiation −
surface_heat_flux_convection = surface_heat_flux_net` for a condition built from
those members. The per-face `BoundaryState` records whether each Newton solve
converged (§12.2.7). The balance of the whole column over a solve — the boundary
fluxes integrated over each accepted step against the inventories and the events —
is the conservation audit (`ConservationRecorder`, `conservation_audit`; §16.8).

---

## 12.10 Worked example: the cone-calorimeter top boundary

Assembling the canonical top boundary makes the pieces concrete. A 50 kW·m⁻² cone
with surface re-radiation and natural convection to a 300 K ambient, on a
sealed-but-venting slab, is

```julia
using Pyrolysis

bc_set = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 50e3) +              # absorbed incident flux, α = a_eff
              RadiativeBC(ε = nothing, T_inf = 300.0) +   # re-radiation loss
              ConvectiveBC(h = 10.0, T_inf = 300.0),      # convective loss
        bottom = AdiabaticBC(),                           # insulated back
    ),
    mass = (
        top = DiffusiveMassBC(),                          # gas vents, internal diffusion limit
        bottom = ImpermeableBC(),                         # sealed back
    ),
)
```

With `radiation_model = SURFACE_ABSORPTION`, the surface energy
balance solved by Newton at the top is

$$
\underbrace{a_{\text{eff}}\,(50{,}000)}_{\text{absorbed}}
\;-\;\underbrace{\varepsilon_{\text{eff}}\,\sigma\,(T_s^4 - 300^4)}_{\text{re-radiation}}
\;-\;\underbrace{10\,(T_s - 300)}_{\text{convection}}
\;=\;
k\,\frac{T_s - T_{\text{cell}}}{d},
$$

with $a_{\text{eff}}=\varepsilon_{\text{eff}}$ (Kirchhoff, both `nothing`). As the
surface chars, $\varepsilon_{\text{eff}}$ tracks the composition automatically,
changing both the absorbed and re-radiated terms. Inside the property bracket the
Newton derivative is $f'(T_s) = -4\varepsilon_{\text{eff}}\sigma T_s^3 - 10 - k/d$
(the absorbed term and the ambient parts are $T_s$-independent), which is strictly
negative — guaranteeing a unique root reached in a few iterations. The conducted flux
$q_{\text{conducted}}$ is written to the top face as the outward flux
$-q_{\text{conducted}}$, and, for a gas-producing material, the vented gas leaves the
top along its outward normal with the flux of `DiffusiveMassBC`, carrying its
enthalpy at $T_s$. If `radiation_model = BEER_LAMBERT` is selected instead, the
`RadiativeFluxBC` contributes zero to this balance and the 50 kW·m⁻² enters the
energy rows as in-depth absorption through the face intensities (§8.3).

On an inert slab ($\rho = 1000$ kg·m⁻³, $c = 1500$ J·kg⁻¹·K⁻¹,
$k = 0.2$ W·m⁻¹·K⁻¹, $\varepsilon = 0.9$; 1 cm in 20 cells, starting at 300 K; the
slab has no gas, so the mass conditions are not read):

```julia
slab = SolidComponent(:slab, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, 20),
    material = Material(name = :inert, components = (slab,), reactions = ()),
    bc_set = bc_set,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1000.0],
    tspan = (0.0, 60.0),
    output_times = [0.0, 30.0, 60.0],
)
sol = solve(problem; radiation_model = SURFACE_ABSORPTION)
x = sol.extras
```

the solve returns `:Success` with

| $t$ [s] | `surface_T` [K] | absorbed [W·m⁻²] | re-radiation [W·m⁻²] | convection [W·m⁻²] |
|---|---|---|---|---|
| 0 | 355.06 | 45 000 | 397.75 | 550.65 |
| 30 | 684.77 | 45 000 | 10 807.43 | 3 847.66 |
| 60 | 758.41 | 45 000 | 16 470.06 | 4 584.06 |

and `x.surface_heat_flux_net .- x.surface_heat_flux_conduction` of
$\pm7.3\times10^{-12}$ W·m⁻² at every saved time. At $t=0$ the cell is at 300 K, yet
the face is already at 355 K: the surface balance places the absorbed flux at the
face, not at the cell center.

---

## 12.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The three reference solvers express the same surface physics in their own
notation; the subsections below map their symbols to ours on first use and
contrast the modeling choices. (Pyrolysis.jl stores the heat of reaction with
$h>0$ endothermic, $Q_{\text{rxn}}=-\sum h_r r_r$, as Gpyro and FDS do; ThermaKin
publishes $h>0$ exothermic (§6.9.2). The sign concerns the volumetric reaction heat,
not the surface BCs.)

### 12.11.1 ThermaKin2Ds

ThermaKin2Ds applies piecewise external/flame boundary conditions and supports an
external radiant source with surface re-radiation. The exposed-face energy
balance, in ThermaKin's form, is

$$
-k\frac{\partial T}{\partial x}\Big|_{\text{surf}}
= a\,\dot q''_{\text{ex}}
  - \varepsilon\sigma\bigl(T_s^4 - T_\infty^4\bigr)
  - h_{\text{conv}}\bigl(T_s - T_\infty\bigr),
$$

with $\dot q''_{\text{ex}}$ the incident external flux, $a$ the absorptivity
(ThermaKin's $\alpha$) and $\varepsilon$ the emissivity (Kirchhoff $a=\varepsilon$). Pyrolysis.jl's
cone BC `RadiativeFluxBC + RadiativeBC + ConvectiveBC` is exactly this
balance, with our $q_{\text{inc}}=\dot q''_{\text{ex}}$ and our
$\varepsilon_{\text{eff}}$ playing the role of ThermaKin's surface $\varepsilon$.
ThermaKin also tracks a reflected component $I_{\text{ex}}^0$; our
composition-weighted $a_{\text{eff}}$ supplies the absorptivity and the reflected
fraction $(1-a_{\text{eff}})$ is implicitly discarded (opaque-surface assumption),
consistent with the `SURFACE_ABSORPTION` model. The shared idea — that both the
absorbed external flux and the surface re-radiation depend on the *true* surface
temperature — is realized by ThermaKin through its Crank–Nicolson surface node and by
Pyrolysis.jl through the surface-balance Newton solve of §12.2. ThermaKin carries gas
enthalpy via a convective term $c_g(J_g\cdot\partial T/\partial x)$ (its Eq. 17),
which corresponds to our gas enthalpy flux at the faces (§7.4) and is the reason the
surface balance carries no gas-enthalpy term (§7.4.3).

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
$q_{\text{BC}}=a_{\text{eff}}q_{\text{inc}}-\varepsilon\sigma (T_s^4-T_\infty^4)-h_{\text{conv}}(T_s-
T_\infty)$. FDS additionally performs a *surface oxygen* Newton iteration for
oxidative pyrolysis (coupling surface $O_2$ to the gas phase); Pyrolysis.jl is a
condensed-phase-only solver with no gas-phase $O_2$ feedback, so our surface
Newton iterates on $T_s$ alone. FDS recovers Arrhenius parameters from TGA
reference temperature/rate (a kinetics, not a BC, concern; see §6, §18). The FDS
in-depth radiation (a two-flux Schuster–Schwarzschild model) corresponds to our
`BEER_LAMBERT` in-depth absorption, under which the `RadiativeFluxBC` boundary
flux is suppressed to avoid double-counting (§12.3.5), just as FDS deposits the
in-depth absorbed fraction in the volume rather than at the face.

---

## 12.12 Limitations and open issues

1. **Temperature-dependent emissivity derivative.** For a `RadiativeBC` whose
   own `ε` is a function of $T_s$, the radiative Newton derivative omits the
   $F\,(d\varepsilon/dT)\,\sigma\,(T_\infty^4-T_s^4)$ term (§12.2.4). The
   converged $T_s$ is unbiased, but convergence can become sub-quadratic for
   strongly temperature-sensitive $\varepsilon(T_s)$. This does not affect
   component-level $\varepsilon_j(T)$: the mixture $\varepsilon_{\text{eff}}$
   is evaluated at the *cell* temperature and is a constant of the surface
   Newton iteration, and the Jacobian differentiates the converged solve exactly.

2. **Dirichlet thermal BC via penalty inside a combination.** A bare
   `DirichletThermalBC` is exact (closed form, §12.3.6). Nested in a
   `CombinedThermalBC` it uses a large conductance $h_{\text{large}}=10^8$ in the
   Newton solve, and the face misses the prescribed temperature by
   $q_{\text{cond}}/h_{\text{large}}$ ($4.4\times10^{-4}$ K at 44 kW·m⁻²). An
   algebraic surface-temperature row would remove the miss and the conditioning
   sensitivity of the penalty.

3. **One radiative-flux member under in-depth radiation.** Under `BEER_LAMBERT`
   only the first `RadiativeFluxBC` member of the `:top` condition feeds the beam,
   and every member contributes zero to the surface balance, so a second
   radiative-flux member is dropped silently (§12.3.5). Give the incident flux as a
   single condition (a function of time can combine several sources).

4. **Half-cell surface temperature.** The face temperature is the one-sided
   extrapolation over the boundary half-cell; a volumetric consumer inside that
   half-cell (a reaction zone, in-depth absorption) biases it at coarse resolution
   (§12.3.2). Refine the near-surface mesh when $T_s$ is a target quantity.

5. **View-factor meaning in 1D.** The `RadiativeBC` view factor $F$ is a scalar
   surrogate for a 3D geometric quantity; in a strictly 1D model it is best read as
   a tuning factor for partial enclosure rather than a rigorously computed exchange
   factor.

6. **Mass-transfer coefficient estimation.** The film coefficient $h_m$ for
   `ConvectiveMassBC` / `DiffusiveMassBC` is supplied by the user, typically via
   the heat–mass analogy $h_m\approx h_{\text{conv}}/(\rho_{\text{air}}
   c_{p,\text{air}})$; there is no built-in correlation, so its accuracy depends on
   the user's estimate of $h_{\text{conv}}$ and the flow regime.

7. **Pressure-BC scope.** The three pressure BCs act only under Darcy–Fick
   transport, i.e. when a condensed component is permeable (`has_darcy_flow`); on a
   Fickian material they are accepted and contribute nothing. Their venting is
   outflow-only: the ambient beyond the face carries none of the tracked gas species,
   so an exterior pressure above the interior pressure draws no gas into the slab
   (§9.5.3). The residual, the Jacobian and the saved-output diagnostics evaluate
   the one function `outward_darcy_velocity` (§12.5).
