# 5. Effective (Mixture) Properties and Mixing Rules

This chapter defines the property bases used in the conservation equations.
Thermal conductivity, emissivity and absorptivity use thermal volume weights;
gas transfer and permeability use condensed bulk occupancy. Skeletal density
sets pore volume independently of both. Heat capacity and sensible enthalpy
are assembled from component masses and caloric laws.

The transport closures and their interface matching are specified here and
in Chapter 9. Property evaluation is implemented in `src/Properties`, with
trial-state extensions and caloric maps in `src/States` and `src/Materials`.

Throughout this chapter `j` is the **component (species) index** running over the
`N_C` components of a `Material`, and `i` is the **cell index** running over the
`n` finite volumes (overload resolution I1, §2). The primary species variable is
the **mass concentration** `ξ_j` (kg·m⁻³), not the mass fraction (§3, Governing
Equations); the integrated state is the cell mass `m_{j,i} = ξ_{j,i} V_i`. Two
density symbols are used with care: the **bulk** (pure-phase) density `ρ_j`
sets declared bulk occupancy in mixing and strain, while the **intrinsic**
(skeletal, cell-wall) density `ρ_{i,j}` appears *only* in the pore-fraction
formula (overload resolution A7, §2). The three binding sign conventions of §2
are assumed without restatement.

The property functions below are written in `(ξ, T)`. The residual evaluates
every one of them once per cell, in the state mapping, at the positive-part
concentrations `ξ⁺_{j,i}` and the trial temperature `T_trial` (§5.13); both
equal `ξ_{j,i}` and `T_i` on every physical state the tolerances can
distinguish.

---

## 5.1 Inventory and property-weighting bases

### 5.1.1 Mass fractions

$Y_j=\xi_j/\sum_k\xi_k$ weights properties expressed per unit total mixture
mass. Heat capacity is assembled directly as $\sum_jm_jc_{p,j}$; the
conservative energy row does not need a normalized mixture heat capacity.

### 5.1.2 Thermal volume fractions

`volume_fractions` normalizes all-phase $\xi_j/\rho_{b,j}$.
`solid_volume_fractions` normalizes all-phase
$\gamma_j\xi_j/\rho_{b,j}$. Thermal conductivity, emissivity and absorptivity
use `mixing_volume_fractions`: the swelling-weighted basis above matrix
occupancy $10^{-3}$, the all-phase basis below $10^{-4}$, and a cubic
smoothstep blend between them. Gas with nonzero `γ` participates in this
thermal basis. The trial extension at a vanishing total occupancy is stated
in the property implementation.

### 5.1.3 Gas matrix weights and skeletal occupancy

Gas transport has a distinct basis: only condensed components enter
$q_j=\gamma_j\xi_j/\rho_{b,j}$, and their normalized weights determine
$\Lambda$ and component-mixed $\kappa$. Gas loading never changes these
weights. At absent matrix support, the specified transport endpoint and
its skeletal-occupancy domain apply (§9.2).

A zero-`γ` condensed component has no matrix weight or direct contribution
to the strain law. It still has mass, heat capacity and skeletal volume
$\xi_j/\rho_{i,j}$. Zero weighted matrix occupancy is therefore insufficient
to identify an unobstructed void.

| Quantity | Basis | Consumer |
|---|---|---|
| Bulk density | $\sum_j\xi_j$ | inventory |
| Thermal fractions | all-phase bulk volume, with thermal blend | conductivity and surface optical properties |
| Gas matrix weights | condensed $\gamma_j\xi_j/\rho_{b,j}$ | $\Lambda$, component-mixed $\kappa$ |
| Skeletal occupancy | condensed $\xi_j/\rho_{i,j}$ | physical porosity and endpoint validity |
| Gas fractions | intrinsic gas inventories | friction, film and reservoir exchange |

---

## 5.2 Volume-change coefficients

The strain law uses $\gamma_j/\rho_{b,j}(T)$ for each component and its
reaction source (§3.6, Chapter 10). A nonzero gas swelling factor can change
bulk geometry without entering gas matrix weights. The pore EOS follows
the actual cell volume minus condensed skeletal volume; it does not replace
that geometric volume by a mixture-specific-volume average.

---

## 5.3 Porosity

Porosity is the complement of condensed skeletal occupancy,
$\phi=1-\sum_{S\cup L}\xi_j/\rho_{i,j}(T)$, independent of `γ` and gas loading.
It is evaluated for every gas-bearing material. The pressure, intrinsic
species density and gas fractions use this same pore state (§9.1).

### 5.3.1 Evaluation and admissibility

Trial concentrations are nonnegative extensions of raw inventories and the
porosity divisor has a C¹ floor at $10^{-12}$. A physical state must lie above
twice that floor plus its mass/temperature/geometry uncertainty margin.
Initial validation, step acceptance and post-deactivation checks enforce this
condition; smoothing does not legalize a filled pore volume.

### 5.3.2 Porosity derivatives

The derivative with respect to condensed concentration is $-1/\rho_{i,j}$,
and the temperature derivative is
$\sum_{S\cup L}\xi_j\rho'_{i,j}/\rho_{i,j}^2$. Inventory and geometry
derivatives follow $\xi_j=m_j/V$. These derivatives enter both the colored
Jacobian and proposed-state pressure-error control.

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
(ρ_g c_{p,g} ≈ 1.2 × 1000). Because every
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
component. Its mass is an integrated state, and no separate gas-temperature
unknown is needed. Gpyro's default thermal-equilibrium exchange neglects both
gas storage and gas conduction, restored by `FULL_QSG`; the standard FDS wall
model instead has no pore-gas inventory. The in-pore gas heat capacity is
roughly `φ ρ_g c_{p,g}` against the matrix `(1−φ) ρ_s c_{p,s}`, with `ρ_s` the
intrinsic density of the matrix. At 1 atm and 700 K, with comparable specific
heats and `ρ_s = 1500` kg·m⁻³, the ratio is about 3·10⁻⁴ for air in a matrix of
porosity 0.5 (ρ_g ≈ 0.50 kg·m⁻³) and about 5·10⁻³ for volatiles of molar mass
0.1 kg·mol⁻¹ in a char of porosity 0.8 (ρ_g ≈ 1.74 kg·m⁻³): from a few hundredths
to about half a percent of the storage, scaling with pressure. These estimates are small,
but our model retains this storage term.

The gas *flux* relative to the skeleton enters the energy equation only through
the sensible enthalpy it carries across the faces (§7.4). This corresponds to the gas temperature-advection term of ThermaKin
Eq. 17 under the matching assumptions in §7.8.

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

## 5.6 Effective superficial gas transfer

`effective_gas_transfer(ξ,T,material)` returns the lumped superficial
coefficient $\Lambda$ [m²/s], including a declared endpoint blend. Only
condensed $\gamma$-weighted bulk occupancies enter its normalized matrix
weights. Participating condensed components require explicit `λ`; gases
have no transfer-coefficient input. The coefficient multiplies the intrinsic
partial-pressure gradient, $J_j^D=-\Lambda M_j\nabla p_j/(RT)$.

`PARALLEL`, `SERIES` and `WEIGHTED` use the exact matrix rules of §9.2.
Positive-weight zero coefficients close the series channel. Support changes
at disappearing constituents require qualification; arbitrary positive floors
are not used. `BRUGGEMAN` is a thermal closure and is rejected for gas transfer.
A matrix-free limit requires an explicit endpoint and occupancy domain.

---

## 5.7 Effective permeability

`effective_permeability(ξ,T,material,φ)` returns $\kappa$ [m²]. Condensed
component mixing uses the same matrix basis as $\Lambda$. The default
`WEIGHTED` rule combines parallel and exact series values. Gases and
zero-`γ` constituents carry no matrix permeability weight.

### 5.7.1 Component mixing

For normalized condensed weights $w_j$, $\kappa_\parallel=\sum_jw_j\kappa_j$,
$\kappa_{series}=(\sum_jw_j/\kappa_j)^{-1}$ and
$\kappa_{weighted}=\beta\kappa_\parallel+(1-\beta)\kappa_{series}$.
A positive-weight zero closes the series value exactly. `κ=0` remains a legal
mechanism switch independently of wall transfer.

### 5.7.2 Carman–Kozeny closure

`MixingSpec(CARMAN_KOZENY;length_scale=ℓ)` instead gives
$\kappa=\min(\ell^2\phi^3/[180(1-\phi)^2],\ell^2)$, with a finite divisor
extension at unit porosity. The length scale is positive physical data.
The cap and its branch derivative are part of this constitutive closure.
It is evaluated even when component permeabilities are zero.

### 5.7.3 Constitutive and interface derivatives

`BRUGGEMAN` is rejected for permeability. State differentiation includes
matrix weights, property laws, pore geometry and endpoint blending. Interior
faces match total flux from the two half paths; independently averaging
permeability and transfer would not preserve mixed-mechanism paths (§9.6).

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
its face under `NO_RADIATION` and `SURFACE_ABSORPTION` (§5.8.2), and under
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
pore-gas viscosity curve in the `gas_viscosity` field, and under gas
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

virgin = SolidComponent(:virgin, ρ = 500.0, c = 1500.0, k = 0.15, λ=2e-5, κ = 1e-15,
                        ρ_intrinsic = 1500.0)
char = SolidComponent(:char, ρ = 150.0, c = 1100.0, k = 0.08, λ=2e-5, κ = 1e-12,
                      ρ_intrinsic = 1500.0)
gas = GaseousComponent(:gas, M = 0.03, c = 1500.0, k = 0.03,)
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

The package provides a binary molecular diffusion estimate from
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
estimate for explicit molecular reference data (§5.10.3), it supplies the gas
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

### 5.10.3 Transport-model inputs

A molecular diffusivity estimate can inform the explicitly supplied `D_ref`
of `DustyGas`, with its `T_ref` and `P_ref`. That model assumes one common
molecular diffusivity and derives species Knudsen coefficients from molar
mass and pore diameter. A fitted lumped coefficient belongs to the condensed
matrix and is superficial; its numerical value must respect that definition.

---

## 5.11 Species dependence

`LumpedTransfer` uses one superficial matrix coefficient for all partial
pressure potentials. `DustyGas` has $K_j\propto M_j^{-1/2}$ and couples all
gas fluxes through a common molecular-friction law. Species mass flux is
$M_j$ times molar flux. Neither model derives matrix weights from gas loading.

---

## 5.12 Conductive averaging and gas interface continuity

The conductive face coefficient uses the distance-weighted harmonic mean,
$k_f=(d_L+d_R)/(d_L/k_L+d_R/k_R)$, with its finite trial extension. This
combines two conductive resistances in series. Pore gas faces instead match
the total species flux through two half-cell paths (§9.6). For an isolated
linear wall channel this reduces to the series resistance
$R(d_LT_L/K_{j,L}+d_RT_R/K_{j,R})$. The coupled wall/molecular/Darcy law
requires the interface pressure state and cannot generally be represented
by separate harmonic means of its mechanisms. All distances use current
geometry, including ALE motion.

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

- `δ_i = reltol · ε_m · m^s_i(0)` for condensed components, equal to their
  mass absolute tolerance (`reltol` and `threshold` of the solve, `m^s_i(0)` the
  cell's initial condensed mass; `src/States/tolerances.jl`). For gas properties
  and transport, use the current pressure-budget scale `δ_{j,i}` of §9.10;
  the separate caloric extension keeps its reference-mass scales;
- `Ṽ_i = smooth_floor(V_i, 10⁻⁸ V_i(0))` the smoothly floored cell volume (`ds.V_eval`);
- `smooth_clamp(x, lo, hi, w)` constant outside `[lo, hi]`, the identity on
  `[lo + w, hi − w]`, and the cubic join `lo + w t²(2 − t)`, `t = (x − lo)/w`
  (mirrored at `hi`) between; `T_trial` is therefore `T_i` on every physical
  state (`trial_temperature`, `src/Materials/guards.jl`);
- `smooth_floor(x, a)` equal to `a` at and below `a`, to `x` at and above `2a`,
  and to the cubic `a(1 + 2t² − t³)`, `t = x/a − 1`, between.

Each guard is C¹ and is exactly the identity on the physical states the
tolerances can distinguish, so a trial state outside the physical domain still
has finite properties with finite derivatives, while the domain test reads the
unguarded quantities (§15.5.4). The kernels read the signed concentrations
`ξ_{j,i} = m_{j,i}/Ṽ_i` in the reaction rates (§6) and relative solid-motion
fluxes (§11); pore-gas transport and the property laws read `ξ⁺`.

### 5.13.2 Cell properties

`derive!` evaluates conductivity, emissivity, absorptivity, caloric properties
and, when present, radiation absorption. Every gas-bearing material also
gets the common physical/trial pore state: porosity, intrinsic gas densities,
partial and total pressures, and fractions. The declared transport model
supplies gas-only wall coefficients, molecular resistance and Darcy mobility.
A sealed gas material retains its EOS. Gas-free materials use the internal
`NoGasTransport` specialization. Properties and derivatives use the same
state and constitutive branches.

---

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

**Heat capacity.** TK Eqs. 2 and 17 include gas sensible storage.
GP Eqs. 11–12 mass-weight the condensed heat capacity and enthalpy; gas
storage is separate and enters thermal equilibrium through `FULL_QSG`
(Eq. 65). FDS's solid heat capacity sums condensed component contributions.
Pyrolysis.jl sums all component masses times their heat capacities.

**Conductivity and emissivity.** TK Eqs. 8–10 blend volume-weighted parallel
and series conductivity limits, corresponding to our `WEIGHTED` rule when
the volume fractions match. TK's mixture gas swelling factor (Eqs. 3–4)
is not generally our independent $\gamma_j$. In particular, `L SWELLING: 0`
does not remove liquid volume from TK Eq. 3. GP Eqs. 7–8 volume-average
condensed conductivities, including an optional pore-radiation term.
Our `PARALLEL` rule matches that averaging only on the same volume basis.
Additional mixing choices do not establish a more accurate physical model.

**Permeability.** GP Eq. 15 volume-averages prescribed component
permeabilities. It does not specify our Carman–Kozeny closure. In Pyrolysis.jl
that closure is an alternative to component mixing (§5.7).

**Absorption units.** TK Eq. 11 and Pyrolysis.jl use mass-specific component
coefficients $\alpha_j$ [m²/kg], with
$\alpha_{\mathrm{eff}}=\sum_j\xi_j\alpha_j$ [m⁻¹]. Gpyro's $\kappa_i$
(GP Eqs. 13–14) and FDS's $\kappa_{s,i}$ are already per length and are
volume-averaged. Match the **effective optical depth**. For a single pure
component at a specified density, $\alpha_j=\kappa_i/\rho_i$; this does not
establish a universal constant conversion for a changing mixture.
The radiation models themselves still differ (§8.7).

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

For `k_eff`, `ε_eff`, `a_eff`, and `α_eff`, temperature dependence enters
through the component laws (§4.5), and composition dependence through the
thermal volume fractions and their void blend (for `α_eff`, through `ξ`
directly). Gas transfer and permeability use condensed bulk-occupancy weights
and the explicitly configured transport endpoint window. Their derivatives
follow those closures, including temperature-dependent pore geometry. Because the properties are evaluated
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
- the numerical thresholds of thermal conductivity mixing;
- exact zero coefficients in a participating transport series channel, whose
  reachable support changes are validated at setup (§9);
- the empty-cell guards of thermal volume fractions (`< 10⁻²⁰`);
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
   change. A gas test with more than one gas component would be needed to
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
