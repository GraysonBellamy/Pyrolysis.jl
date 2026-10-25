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
material emissivity evaluated from the boundary cell, and the common energy
control surface for incoming/outgoing gas and thermal-film heat (§9.9).

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

### 12.1.1 Boundary control volumes

The finite-volume rows share one outward mass vector and one outward energy
flux at every live boundary. Thermal flux laws use the inward-positive
convention, and assembly converts them to outward energy. A boundary may
couple pore transport, external mass transfer, hydraulic resistance and a
thermal film. Their interface state is solved jointly (§12.2); adiabatic and
prescribed-heat conditions use the direct donor-enthalpy balance.

Surface temperature is an interface variable distinct from the adjacent
cell temperature. Conduction through the current half-cell distance matches
the supplied heat law. This captures nonlinear radiation and temperature-
dependent blowing without evaluating the heat law at the cell center.

---

### 12.1.2 Type hierarchy and dispatch

Thermal conditions subtype `AbstractThermalBC` and may be combined additively.
Mass conditions subtype `AbstractMassBC`: `ImpermeableBC`, `MassFluxBC`,
`AmbientGasBC`, `NeumannPressureBC` and `ConvectivePressureBC`. Each mass
condition supplies one complete closure and is not additive. `GasReservoir`
and `MassTransferFilm` supply its thermodynamic state and optional film data.
Typed boundary bindings retain gas indices, molar masses and species names.

---

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
condition laws define them: inward for heat and outward for mass. Assembly
converts heat and the live-face orientation to the stored outward convention.

A patch applies on the **live face** of its tag (`Meshes.live_face`): the boundary
face itself while its cell is live, and after cells beyond it are deactivated
(Chapter 14), the face of the outermost live cell — an interior face whose stored orientation may differ from the outward
normal of the remaining live cell. Assembly applies that orientation explicitly.

---

## 12.2 Coupled surface balance

### 12.2.1 Signed energy balance

Heat $q_s$ is positive into the material and species fluxes are positive
outward. Conduction uses $k(T_s-T_{cell})/d$ with current face distance.
For a wall thermal film, its heat correlation equals this conductive flux.
Boundary assembly converts to the direct-reservoir energy control surface:
$q_r=q_s+H_r-H_s$ and $E_{out}=H_r-q_r$. Inflow carries the specified incoming
enthalpy even under adiabatic or prescribed thermal conditions (§9.9).

![Radiation, convection and conduction at the surface, with donor carried enthalpy and the conversion from thermal-film heat to the reservoir energy control surface.](../../assets/figures/06_surface_energy.svg)

### 12.2.2 Local residuals and qualification

The unknowns are $T_s$, $N_g-1$ independent surface fractions when a finite
mass film is present, and $P_s$ for a hydraulic resistance. Coefficients,
Darcy donor density and blowing corrections are reevaluated together. The
species residual matches pore and film flux; the hydraulic residual matches
pressure drop and the selected flow; the thermal residual matches conduction
and applied wall heat, or the prescribed surface temperature.

Residuals and unknowns are scaled. Newton trials preserve nonnegative pressure
and the composition simplex, with line search and a relative trust bound.
Temperature-dependent properties use the finite trial-temperature extension;
physical accepted states are checked separately. A converged residual, resolved
correction and nonsingular local Jacobian are required. Failure is explicit;
no exhausted Newton iterate supplies a flux. The initial guess is independent
of diagnostic `BoundaryState` contents.

### 12.2.3 Derivatives and output

First derivatives are obtained from the implicit function theorem at the
qualified primal root. They include species/temperature coupling and blowing
properties, with deterministic primal branch choices. Degenerate channels and
vacuum branches have their stated explicit continuations (§9.5–9.6).
`get_boundary_exchange` reports the assembled fluxes, energy and local residual
and iteration counts; `BoundaryState` records thermal diagnostics.

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

**Blowing.** Each `ConvectiveBC` carries a model, default `NoBlowing()`.
The coupled boundary solve supplies signed outward species fluxes, surface
temperature and, when available, the solved surface composition. The model
returns $R_b=h_{eff}/h$, giving wall-film heat
$q_s=hR_b(T_\infty-T_s)$. Reservoir-control-surface heat then follows the
enthalpy conversion in §9.9. Relative mesh transport is handled separately.

For `CouetteFilm(λ)`,

$$
C=\sum_gJ_gc_{p,g}(T_s),\qquad
B=2\lambda C/h,\qquad R_b=\beta(B),\qquad
\beta(B)=\frac{B}{e^B-1},\quad\beta(0)=1.
$$

The classical factor is $\lambda=1/2$. The signed sum remains well defined
for counterflow, and suction enhances heat transfer. Zero net mass flow
need not give $C=0$ when species heat capacities differ. Surface heat
capacities give a frozen-property approximation for temperature-dependent
properties. The stable Taylor/exponential branches of `blowing_factor`
avoid cancellation at zero and overflow under strong blowing.

`PropertyRatioFilm(free_stream; n_Γ)` modifies the free-stream Couette
parameter using $(\Gamma_\infty/\Gamma_g)^{n_\Gamma}$, where $\Gamma=k/c_p$.
The outgoing mixture uses mass-weighted heat capacity and Wilke conductivity.
`PoolingFilm(free_stream; n_Γ,n_D,ψ0,Λc,Λw)` additionally applies a retention
factor $\psi$, returning $\psi\beta(B_\Gamma/\psi)$. It computes buoyancy from
the solved surface composition when one is supplied. Without that state it
uses its empirical outflow/exchange composition closure and requires zero
free-stream species flux. Both empirical models support outward products or
pure free-stream suction; they diagnose unsupported mixed counterflow.
Their fitted parameters describe the external apparatus, not pore transport.
The signed Couette model supplies the counterflow closure without those
empirical composition assumptions.

`FunctionBlowing(f)` evaluates `f(face,h,T_inf)` with signed species fluxes,
surface temperature, material properties and the available surface composition.
Each convective member of a combined condition applies its own model. The
coupled root's Jacobian and implicit derivatives include the full dependence
of mass flux and film response on the solved state (§12.2).

`BoundaryState` records `blowing_mcp`, `blowing_ratio`, `pooling_factor`,
`face_fraction` and `face_buoyancy`. Saved solution fields include
`surface_blowing_ratio`, `surface_pooling_factor`, `surface_face_fraction`,
`surface_face_buoyancy` and per-tag `boundary_blowing_ratio`. Inspect
`get_boundary_exchange` for the actual transmitted species and energy fluxes.
The types and gas data are in `src/Materials/blowing_models.jl`; evaluation
is in `src/BoundaryConditions/blowing.jl`.

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
$\varepsilon_{\text{eff}}$ of §5.8.2), and the surroundings temperature
$T_\infty$. This term is the surface's own re-radiation to and absorption from a
gray environment; it is the dominant cooling mechanism for a hot pyrolyzing
surface. Its quartic temperature dependence is why the face temperature, not the
cell-center temperature, must enter it (§12.1.1). The condition stays at the
surface under every radiation model. (`surface_balance.jl`.)

The displayed slope holds emissivity fixed. If emissivity depends on surface
temperature, the full derivative additionally contains
$F\sigma\varepsilon'(T_s)(T_\infty^4-T_s^4)$. The surface-root Jacobian and
its implicit derivative include that term by differentiating the flux residual.
The low-level flux helper also returns a frozen-coefficient slope; this slope
is not used to differentiate the qualified surface root.

### 12.3.5 `RadiativeFluxBC` (external incident radiation; model-dependent)

This BC represents incident radiation from an *external* source (cone heater,
furnace, radiant panel) — not the surface's own emission. The absorbed fraction is
$a_s\,q_{\text{inc}}$, where $q_{\text{inc}}=$ `flux` (constant or $q(t)$) and
the absorptivity $a_s$ (keyword `absorptivity`) is either supplied (a constant or a function of $t$) or,
when `nothing`, the mixture surface absorptivity $a_{\text{eff}}$, which equals
$\varepsilon_{\text{eff}}$ by Kirchhoff's law unless a component sets its own
`absorptivity` (§5.8.2).

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

## 12.4 Mass and reservoir boundary conditions

`ImpermeableBC()` gives zero species flux. `MassFluxBC(fluxes=(gas=...,...),
inflow_temperature=(gas=...,...))` prescribes named outward mass fluxes in
kg/(m² s); each negative flux needs its incoming temperature.

`GasReservoir(P=...,T=...,mole_fractions=(gas=...,...))` defines pressure [Pa],
incoming temperature [K] and normalized named gas mole fractions. Time-dependent
values are permitted and checked when evaluated. Only an explicit vacuum
`P=0` may omit composition. `AmbientGasBC(reservoir;film=nothing)` connects
the complete state to the pore path and permits both ingress and escape.

`MassTransferFilm(g)` adds the isobaric common-diffusivity film of §9.8, with
conductance in kg/(m² s), a number or temperature law. Zero conductance permits
only donor bulk flow; infinite conductance fixes composition. Pressure, species
and thermal equations are solved consistently. These conditions do not combine
with `+`; one closure is assigned to each tag.

---

## 12.5 Pressure-gradient and hydraulic-resistance boundaries

`NeumannPressureBC(gradient=...,T=...,mole_fractions=...)` sets the outward
pressure gradient [Pa/m] and provides the species composition and incoming
energy state. The surface pressure is the cell pressure plus gradient times
the current half-cell distance, and physical accepted pressures must be
nonnegative.

`ConvectivePressureBC(reservoir;resistance=...,flow=:darcy_velocity,film=nothing)`
solves $P_s-P_\infty=R_h f$. The flow is Darcy velocity [m/s], total mass flux
[kg/(m² s)] for `:mass`, or total molar flux [mol/(m² s)] for `:molar`.
Resistance is finite and nonnegative, with units pressure divided by the chosen
flow. A zero resistance fixes surface pressure. The complete reservoir and
optional film provide species and energy closures in either direction.
These pressure conditions are valid with either transport model.

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
tag of the mesh, and binds reservoir species names to the material gas components.

### 12.7.3 Typed bindings and patches

`bc_set.thermal[tag]` and `bc_set.mass[tag]` provide the assigned conditions.
Setup binds their gas species names and molar masses once. Boundary patches
follow the live face of each tag after deactivation; flux orientation is
converted to the stored mesh-face orientation before the row scatter.

---

### 12.7.4 Constructor parameter validation

Individual BC constructors optionally validate their physical parameters at build
time (`validate_bc_params`, `boundary_conditions.jl`): $h\ge0$, $T>0$ K,
$\varepsilon\in[0,1]$, $F\in[0,1]$, $a_s\in[0,1]$. Gas boundary constructors
validate their own nonnegative pressure and conductance, positive temperature,
and normalized composition; positive-pressure reservoirs require composition.
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
converged (§12.2.3). The balance of the whole column over a solve — the boundary
fluxes integrated over each accepted step against the inventories and the events —
is the conservation audit (`ConservationRecorder`, `conservation_audit`; §16.8).

---

## 12.10 Example: an atmospheric heat and mass film

A material with an `:air` gas can exchange with this complete reservoir:

```julia
reservoir = GasReservoir(P=101325., T=300., mole_fractions=(air=1.,))
bcs = BoundaryConditionSet(
    thermal=(top=HeatFluxBC(30e3) + ConvectiveBC(h=10., T_inf=300., blowing=CouetteFilm()),
             bottom=AdiabaticBC()),
    mass=(top=AmbientGasBC(reservoir; film=MassTransferFilm(.06)),
          bottom=ImpermeableBC()))
```

The mass conductance is in kg/(m² s). The common incoming stream temperature
is 300 K in both films. Other gases in the material have zero reservoir mole
fraction and can escape; air can enter. Initialize the pore gas inventories
from the declared initial pressure, temperature and physical porosity.
See the executable cone-calorimeter examples for complete material data.

---

## 12.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

### 12.11.1 ThermaKin2Ds

TK §2 defines time- and position-dependent external heating and empirical
flame feedback. Eq. 19 permits a linear gas-transfer boundary law and an
exponential surface-loss law. These are broader than the boundary models
used in our usual cone example. TK's $\alpha_i$ is a mass-specific absorption
coefficient, not dimensionless surface absorptivity. $I_{ex}^0$ denotes the
incident flux **after reflection**, not the reflected component.

TK Eqs. 12 and 23 and its MAX/RAND algorithm use absorption-weighted or
selected-element radiation properties. Only an appropriate surface-absorption
limit resembles our `RadiativeFluxBC + RadiativeBC + ConvectiveBC` balance;
the general TK model is not exactly that surface equation (§8.7.1).
The manual does not establish that TK uses our separate surface-temperature
Newton solve or an equivalent surface-node construction.

### 12.11.2 Gpyro

GP §5 describes convective/radiative thermal boundaries and alternative gas
pressure, flux and species conditions. Its 1D discretization puts temperature
nodes at the exposed and back boundaries in half cells (GP §4.1.1); a surface
temperature is not obtained by harmonic interpolation from a cell center.
Our boundary temperature is recovered by a one-sided surface balance.

The entering radiant flux includes absorptivity, and surface versus in-depth
deposition must be selected consistently. An impermeable back is one
configuration, not a general characterization of every Gpyro deck.
Pyrolysis.jl likewise provides separate species and pressure boundary
conditions; their matching requires checking the actual deck.

### 12.11.3 FDS

The FDS solid-phase chapter imposes conducted flux equal to convective plus
net radiative surface flux for opaque material. It sets the surface radiation
source to zero when in-depth radiation carries that exchange. Our
Beer–Lambert mode also avoids depositing the incident beam twice, but its
one-way absorption model differs from FDS's emitting two-flux model (§8.7.3).
FDS's Newton iteration for surface oxygen in oxidative reactions is separate
from our Newton solve for surface temperature.

For FDS 6.11.0, `READ_SURF` rejects `NET_HEAT_FLUX` or
`CONVECTIVE_HEAT_FLUX` combined with `MATL_ID` (error 306). This explains
the retained V-HC-03 overlay choice; it is a version-specific input
restriction, not a claim that all possible FDS configurations or later
versions cannot impose a net flux.

Sources: TK §§2–3; GP §§4.1.1, 5; [FDS boundary input validation](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90#L8290-L8300)
and solid-phase chapter (§1.6.5).

---

## 12.12 Boundary-model domain and numerical qualification

Mass films assume isobaric transport with one common diffusivity. Their
conductance is explicit data, supplied as a constant or temperature law.
A hydraulic resistance specifies which of Darcy velocity, mass flux or molar
flux it relates to pressure drop. A single surface temperature represents
LTE at the wall; incoming enthalpy uses the complete reservoir state.

The local nonlinear solve requires a physical, nonsingular root. Its
residual, correction and implicit derivative are qualified numerically;
there is no global uniqueness claim for arbitrary constitutive laws.
Failures report `TransportFailure`. Vacuum and closed channels use explicit
continuations, and empirical blowing laws reject mixtures outside their
stated composition domain. Refine the boundary half-cell when a reaction
zone or in-depth absorption makes its temperature variation significant.
