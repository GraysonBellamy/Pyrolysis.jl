# 7. Heat Conduction, Gas Convection, and Energy Assembly

This chapter derives the condensed-phase energy equation as Pyrolysis.jl solves it: the
non-conservative temperature equation, in which the volumetric thermal storage sums the
heat capacity of every component, pore gas included, while the sensible heat carried by the
gas flux relative to the skeleton is transported by an explicit volumetric advection source.
We develop each term from its continuous form to its finite-volume discretization — Fourier
conduction with distance-weighted harmonic-mean face conductivity, the discrete divergence
operator, the gas convective source $S_{\text{conv}}$, the reaction-heat source
$Q_{\text{rxn}}$, and the radiation source $Q_{\text{rad}}$ — and we trace the exact
eight-step assembly sequence in `compute_rhs!` that produces $dT_i/dt$ and
$d\xi_{j,i}/dt$. Particular care is taken with the three sign conventions of §2
(Nomenclature), with the precise set of terms the code includes and excludes, and with why
the surface energy balance carries no gas-enthalpy term alongside $S_{\text{conv}}$. Where
useful, the Pyrolysis.jl forms are mapped to and
contrasted with ThermaKin2Ds, Gpyro, and FDS.

Index conventions for this chapter follow §2 overload note I1: **$i$ is the cell index**
($1\le i\le n$) and **$r$ is the reaction index**; the component/species index is $j$
($1\le j\le N_C$). The spatial coordinate is $z$ with $z=0$ the bottom/substrate (boundary
id 2) and $z=L$ the top/exposed surface (boundary id 1); heat enters from the top. Positive
flux is transport in the $+z$ direction. The symbol $x$ appears in a few places (e.g.
$\partial T/\partial x$) only because the implementation uses it as a synonym for the
through-thickness coordinate $z$; the two are identical in this 1D model.

> **Frame of this chapter.** The Eulerian energy balance developed here has
> no solid-convection term — it assumes the condensed phase carrying
> $\rho c_p^{\text{eff}}$ is **stationary in the lab frame**. On a fixed
> (Eulerian) mesh that is exact. Under ALE the matrix *moves* (it subsides
> with the reaction-driven mesh), and the correct time derivative of the
> matrix-carried temperature is the *mesh-frame* one, with no additional
> advection term — see §11.2.2 for the relative-velocity argument.

---

## 7.1 Fourier Conduction and the Face Heat Flux

### 7.1.1 Continuous law

Conductive transport of sensible heat in the condensed phase obeys Fourier's law,

$$
\mathbf{q}_{\text{cond}} = -k(T,\xi)\,\nabla T ,
\qquad
q_{\text{cond}}(z) = -k(T,\xi)\,\frac{\partial T}{\partial z}
\quad [\text{W}\,\text{m}^{-2}],
$$

with the minus sign expressing the second law: heat flows down the temperature gradient, in
the direction of decreasing $T$. The effective conductivity $k=k_{\text{eff}}(T,\xi)$ is the
composition- and temperature-dependent mixture conductivity built from the component
conductivities $k_j(T)$ by the chosen mixing rule (PARALLEL / SERIES / WEIGHTED / BRUGGEMAN)
on the **solid** volume-fraction basis; its construction is the subject of §5 (Effective
Properties & Mixing Rules), and only its *face value* enters here.

### 7.1.2 Discrete face flux: distance-weighted harmonic mean

In the cell-centered finite-volume layout (§13, Finite-Volume Discretization) the unknown
$T$ is stored at cell centers $z_i^c$, and fluxes are required at the faces that separate
adjacent cells. For an interior face $f$ separating a left cell $L$ and a right cell $R$,
with cell-center temperatures $T_L,T_R$ and center-to-face distances $d_L,d_R$, the
conductive flux is

$$
\boxed{\;
q_f \;=\; -\,k_f\,\frac{T_R-T_L}{d_L+d_R}
\;}
\tag{7.1}
$$

where the **face conductivity** $k_f$ is the **distance-weighted harmonic mean** of the two
cell-center conductivities:

$$
\boxed{\;
k_f \;=\; \frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
\;}
\tag{7.2}
$$

The denominator $d_L+d_R$ is exactly the center-to-center distance $d_{LR}=|z_R^c-z_L^c|$, so
on a mesh where each cell center sits midway between its faces, $d_L+d_R=d_{LR}$ and (7.1)
reduces to $q_f=-k_f(T_R-T_L)/d_{LR}$.

**Derivation of the harmonic mean from flux continuity.** Treat the half-cell on either side
of the face as a thermal resistance. With cross-sectional area $A$ the resistance of the left
half is $R_L=d_L/(k_L A)$ and of the right half $R_R=d_R/(k_R A)$. At steady state the same
flux $q_f$ passes through both half-cells in series, so the temperature drop from $T_L$ to
$T_R$ across the combined resistance is

$$
T_L-T_R \;=\; q_f\,A\,(R_L+R_R)
\;=\; q_f\left(\frac{d_L}{k_L}+\frac{d_R}{k_R}\right).
$$

Solving for $q_f$ and writing it as an effective conductivity $k_f$ acting over the total
distance $d_L+d_R$,

$$
q_f \;=\; -\,\frac{T_R-T_L}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}}
\;\equiv\; -\,k_f\,\frac{T_R-T_L}{d_L+d_R},
\qquad
k_f=\frac{d_L+d_R}{\dfrac{d_L}{k_L}+\dfrac{d_R}{k_R}},
$$

which is (7.1)–(7.2). The harmonic (resistance-in-series) mean is the *only* face averaging
that makes the discrete flux exact for a piecewise-constant conductivity field (a sharp
material interface at the face): an arithmetic mean would over-predict the flux through the
more-resistive layer. This is the same series-resistance averaging applied to the
gas-transfer coefficient $\lambda$ and the permeability $\kappa$ (§9, Gas-Phase Mass
Transport).

**Uniform-mesh limit.** When the cell centers are equidistant from the face, $d_L=d_R=d$,
(7.2) collapses to the standard (unweighted) harmonic mean,

$$
k_f \;\xrightarrow{\,d_L=d_R\,}\; \frac{2d}{\dfrac{d}{k_L}+\dfrac{d}{k_R}}
\;=\; \frac{2k_Lk_R}{k_L+k_R}.
$$

On **stretched** meshes — for example, boundary-refined grids used to resolve a steep
near-surface thermal gradient — $d_L\ne d_R$, and the distance weighting corrects the
averaging. Empirically the distance-weighted form recovers $1$–$5\%$ accuracy in the
near-boundary fluxes relative to the unweighted harmonic mean on such meshes, at no extra
cost.

**Small-value guard.** The implementation of (7.2) returns $0$ if either $k_L<10^{-20}$ or
$k_R<10^{-20}\;\text{W}\,\text{m}^{-1}\text{K}^{-1}$, preventing a division blow-up when a
cell has been driven to a vanishing conductivity (e.g. a near-fully-depleted cell). This is a
numerical safety net; physically realistic conductivities are many orders above the
threshold.

> **Implementation note.** The single-face flux is `conductive_flux(k_face, T_L, T_R, d_LR)`
> in `src/Physics/heat_conduction.jl`; it consumes the face conductivity $k_f$ that
> `update_properties!` pre-computes with `harmonic_mean_weighted`
> (`src/Properties/effective_properties.jl`), where the `<1e-20` guard lives.

### 7.1.3 Interior flux pass and the boundary-face zeroing convention

All interior face fluxes are computed in one pass before any boundary condition is applied,
by `compute_conductive_fluxes!` (`src/Physics/heat_conduction.jl`). For each face $f$ the
routine reads the left/right cell ids from the mesh connectivity; for an **interior** face it
reads $T_L,T_R$ from the state, the center-to-center distance $d_{LR}=|z_R^c-z_L^c|$ from the
(possibly ALE-moving) cell-center array, and the cached face conductivity $k_f$, then writes
$q_f$ via (7.1). For a **boundary** face — one whose left or right cell id is $0$ — it writes

$$
q_f \;=\; 0 .
\tag{7.3}
$$

The zeroing in (7.3) is deliberate, not a skip: the boundary face flux is **overwritten**
later (Step 4 of §7.6) by the surface energy-balance solver, which computes the true
conductive flux into the boundary cell from the solved surface temperature $T_s$ (see §12,
Boundary & Initial Conditions). Zeroing first guarantees that if the BC step did not touch a
face, that face contributes nothing rather than a stale value, and it keeps the interior
conduction operator and the boundary operator cleanly separated for Jacobian assembly (the
`ConductionOperator` fills only interior tridiagonal entries; the `BCThermalOperator` fills
the boundary self-coupling — see §15, Temporal Integration & Structured Jacobian).

The face conductivities $k_f$ themselves are pre-computed once per RHS evaluation in
`update_properties!` and stored in the `PropertyCache` field `k_face`, so the flux pass does
no mixing-rule work. The property cache is parameterized on the live scalar type $T_u$
(`Float64` in normal operation, `ForwardDiff.Dual` during automatic differentiation), so the
entire conduction computation carries forward-mode sensitivities without re-allocation.

---

## 7.2 The Discrete Divergence Operator

The conductive flux enters the energy balance through its divergence. Integrating
$\nabla\!\cdot\mathbf{q}_{\text{cond}}$ over cell $i$ and applying the divergence theorem
gives the net outflow through the two bounding faces; dividing by the cell volume $V_i$ yields
the cell-averaged divergence

$$
\boxed{\;
(\nabla\!\cdot\mathbf{q})_i \;=\; \frac{1}{V_i}\bigl(q_{f_R}-q_{f_L}\bigr)\,A
\;}
\tag{7.4}
$$

where $q_{f_L},q_{f_R}$ are the fluxes at the left and right faces of cell $i$, $A$ is the
face area (read from the state as `A_face[f_left]`, which may vary in time in ALE mode), and
$V_i$ is the cell volume.

**Sign convention (outflow-positive).** With the global rule "positive flux = transport in
$+z$," the right face flux $q_{f_R}$ is an *outflow* (it carries heat out of cell $i$ in the
$+z$ direction) and the left face flux $q_{f_L}$ is an *inflow*. Hence the $(q_{f_R}-q_{f_L})$
ordering makes $(\nabla\!\cdot\mathbf{q})_i>0$ when there is a *net loss* of heat from the
cell. Because the energy equation carries $-(\nabla\!\cdot\mathbf{q})_i$ (see §7.3), a positive
divergence (net outflow) correctly *reduces* $dT_i/dt$: a cell losing more heat than it gains
cools.

The same operator applies component-wise to the gas species fluxes,

$$
(\nabla\!\cdot\mathbf{J}_j)_i \;=\; \frac{1}{V_i}\bigl(J_{j,f_R}-J_{j,f_L}\bigr)\,A ,
\qquad j=1,\dots,N_C ,
\tag{7.5}
$$

with units $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ when $J$ is in $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$ — i.e.
$(\nabla\!\cdot\mathbf{J}_j)_i$ is the volumetric net mass-efflux rate of species $j$.

> **Implementation note.** `divergence_1d!` (scalar heat flux) and `divergence_species_1d!`
> (matrix species flux, `[N_C × n_cells]`) are in `src/Discretization/finite_volume.jl`. Both
> use the *current* state volumes and face areas, so in ALE mode they automatically pick up
> the geometric volume change; in Eulerian mode the volume is constant. The species variant
> uses the identical formula component-wise inside an inner `@simd` loop.

A constant flux field gives zero divergence everywhere, and a telescoping sum
$\sum_i (\nabla\!\cdot\mathbf{F})_i V_i = F_{\text{boundary}}^{R}-F_{\text{boundary}}^{L}$
recovers the global balance — the discrete conservation property exploited by the diagnostics
ledgers (§16, Conservation Diagnostics).

---

## 7.3 Effective Heat Capacity and the Energy Equation

### 7.3.1 The non-conservative temperature equation

Pyrolysis.jl solves the energy balance in **non-conservative (temperature-form)** shape:

$$
\boxed{\;
\rho c_{p}^{\text{eff}}(T,\xi)\,\frac{\partial T}{\partial t}
\;=\;
-\,\nabla\!\cdot\mathbf{q}_{\text{cond}}
\;+\; S_{\text{conv}}
\;+\; Q_{\text{rad}}
\;+\; Q_{\text{rxn}}
\;}
\tag{7.6}
$$

with, term by term,

| term | expression | role | § |
|---|---|---|---|
| $\mathbf{q}_{\text{cond}}$ | $-k\,\nabla T$ | Fourier conduction | 7.1 |
| $S_{\text{conv}}$ | $-\sum_g c_{p,g}(T)\,\mathbf{J}_g\!\cdot\!\nabla T$ | gas sensible-enthalpy advection | 7.4 |
| $Q_{\text{rad}}$ | model-dependent | volumetric radiation absorption | 7.5, §8 |
| $Q_{\text{rxn}}$ | $-\sum_r h_r r_r$ | reaction heat | 7.5 |

The subscript $g$ runs over the **gas** components only (the predicate `is_gas`). The storage
coefficient $\rho c_p^{\text{eff}}$ is the **mixture** volumetric heat capacity, defined
next.

### 7.3.2 Mixture heat capacity: every component stores sensible heat

The volumetric thermal storage sums every component, gas included:

$$
\boxed{\;
\rho c_p^{\text{eff}}(T,\xi)
\;=\!\!\sum_{j=1}^{N_C}\!\! \xi_j\,c_{p,j}(T)
\;}
\qquad [\text{J}\,\text{m}^{-3}\text{K}^{-1}]
\tag{7.7}
$$

The pore gas is at the matrix temperature (one-temperature model, §3.1.2) and its
concentration $\xi_g$ is an integrated state, so its storage term costs nothing and is
kept. Dropping it — the *quasi-steady-gas approximation* of Gpyro's default
thermal-equilibrium form (restored there by `FULL_QSG`) and the implicit convention of
FDS, which has no in-depth gas — would be an approximation with no benefit. At
$1\,\text{atm}$ the in-pore gas heat capacity is a few tenths of a percent of the matrix
heat capacity (the gas-to-solid density ratio is $\mathcal{O}(10^{-3})$ and the specific
heats are comparable), scaling with pressure; on white-pine cone decks it is at most
about half a percent and shifts the response by a uniform time lag of about a tenth
of a percent.

The gas enters the energy equation in two places and two places only: its storage in
(7.7) and its *flux relative to the skeleton* through the advective source
$S_{\text{conv}}$ (§7.4). The derivation of §3.3.2 shows these two terms, plus the
reaction heat with heats quoted at the reaction temperature, to be the complete
equation — ThermaKin Eq. 17 term for term. In ALE the whole $\rho c_p^{\text{eff}}$
rides the mesh with the skeleton, so there is no mesh-motion term in the energy
equation (§11.2.2).

> **Implementation note.** $\rho c_p^{\text{eff}}$ is `effective_heat_capacity` in
> `src/Properties/effective_properties.jl`, a `@generated` unrolled sum over the
> component tuple with no phase branch. The per-cell value used by the assembly is
> cached as `prop_cache.ρcp[i]`; the same function feeds the structured-Jacobian
> operators (through Dual evaluation) and the energy ledger; the park-and-mask
> temperature merge (`merge_temperatures_conservatively`, `src/Adaptivity/cell_merge.jl`)
> sums the same every-component heat capacities directly.

### 7.3.3 The discrete temperature update

Substituting the discrete divergence (7.4) into (7.6) and dividing by the cached mixture heat
capacity gives the cell temperature time derivative used by the integrator:

$$
\boxed{\;
\frac{dT_i}{dt}
\;=\;
\frac{1}{\rho c_{p,i}^{\text{eff}}}
\Bigl(
-\,(\nabla\!\cdot\mathbf{q})_i
\;+\; S_{\text{conv},i}
\;+\; Q_{\text{rad},i}
\;+\; Q_{\text{rxn},i}
\Bigr)
\;}
\tag{7.9}
$$

Equivalently, expanding the divergence,

$$
\frac{dT_i}{dt}
=
\frac{1}{V_i\,\rho c_{p,i}^{\text{eff}}}
\Bigl[
-\bigl(q_{f_R}-q_{f_L}\bigr)A
+ V_i\bigl(S_{\text{conv},i}+Q_{\text{rad},i}+Q_{\text{rxn},i}\bigr)
\Bigr].
$$

The three volumetric source terms are *cell-centered* (units $\text{W}\,\text{m}^{-3}$); only
the conduction term arrives as a flux divergence. The exact code is

```julia
ρcp = prop_cache.ρcp[i]
dT[i] =
    (
        -rhs_cache.div_q[i] +
        rhs_cache.S_conv[i] +
        rhs_cache.Q_rad[i] +
        rhs_cache.Q_rxn[i]
    ) / ρcp
```

so the storage coefficient enters at exactly one place — the divisor `ρcp`, the
mixture value (7.7).

### 7.3.4 Species equations (companion to the energy equation)

The same RHS pass assembles the species derivatives. With the gas-vs-solid dispatch (transport
for gases only — see §3, Governing Equations, and §9),

$$
\frac{d\xi_{j,i}}{dt}
=
\begin{cases}
-\,(\nabla\!\cdot\mathbf{J}_j)_i + S_{j,\text{rxn},i}, & \text{gas } j,\\[4pt]
S_{j,\text{rxn},i}, & \text{solid/liquid } j,
\end{cases}
\tag{7.10}
$$

where $S_{j,\text{rxn},i}=\sum_r \nu_{j,r}\,r_r$ is the stoichiometric reaction source (§6,
Reaction Kinetics). Solids and liquids have no transport flux ($\mathbf{J}_j\equiv 0$); they
move only through reactions (and, in ALE, with the mesh — §10, §11). This dispatch is
performed by the `@generated` `apply_species_residual_unrolled!` (`src/Materials/materials.jl`),
which emits `dξ[j,i] = -div_J[j,i] + dξ_rxn[j,i]` for gas components and
`dξ[j,i] = dξ_rxn[j,i]` otherwise, with the gas/solid branch resolved at compile time.

---

## 7.4 The Gas Convective (Advection) Source $S_{\text{conv}}$

### 7.4.1 Origin: the advection piece of the gas-enthalpy divergence

In a fully conservative energy equation, the sensible enthalpy carried by the gas appears as
the divergence of an enthalpy flux,

$$
\nabla\!\cdot\!\Bigl(\sum_g \Delta h_g(T,T_{\text{ref}})\,\mathbf{J}_g\Bigr),
\qquad
\Delta h_g(T,T_{\text{ref}}) = \int_{T_{\text{ref}}}^{T} c_{p,g}(\tau)\,d\tau .
$$

Expanding this divergence with the product rule (exact even for temperature-dependent
$c_{p,g}$, since $\nabla\Delta h_g = c_{p,g}(T)\,\nabla T$),

$$
\nabla\!\cdot\!\Bigl(\sum_g \Delta h_g(T,T_{\text{ref}})\,\mathbf{J}_g\Bigr)
=
\underbrace{\sum_g c_{p,g}\,\mathbf{J}_g\!\cdot\!\nabla T}_{\text{advection}}
\;+\;
\underbrace{\sum_g \Delta h_g(T,T_{\text{ref}})\,\nabla\!\cdot\!\mathbf{J}_g}_{\text{generation/divergence}} .
\tag{7.11}
$$

The energy equation keeps the **advection** piece as a volumetric source on the right-hand side of
the temperature equation, with a leading minus sign (because the conservative term sits on the
storage/transport side; moving it to the source side flips the sign):

$$
\boxed{\;
S_{\text{conv}} \;=\; -\sum_g c_{p,g}(T)\,\mathbf{J}_g\!\cdot\!\nabla T
\;}
\qquad [\text{W}\,\text{m}^{-3}]
\tag{7.12}
$$

The second (generation) piece of (7.11) needs no term of its own: with heats of reaction
quoted at the reaction temperature it is already part of $Q_{\text{rxn}}$ (§3.3.2, §6.6.3).
Physically, $S_{\text{conv}}$ is the rate at which flowing gas deposits or removes
sensible heat as it crosses a temperature gradient: gas moving from hot toward cold heats the
cell it passes through, and vice versa. It carries the enthalpy of the gas *flux* relative to
the skeleton; the gas *storage* sits in $\rho c_p^{\text{eff}}$ (§7.3.2).

### 7.4.2 Discretization

The cell-centered discrete form is

$$
\boxed{\;
S_{\text{conv},i}
\;=\;
-\sum_{g:\,\text{is\_gas}} c_{p,g}(T_i)\,\overline{J}_{g,i}\,
\left(\frac{\partial T}{\partial z}\right)_i
\;}
\tag{7.13}
$$

with two cell-centered ingredients:

**(a) Average gas flux $\overline{J}_{g,i}$.** The flux that "passes through" cell $i$ is the
arithmetic mean of the fluxes at its two bounding faces,

$$
\overline{J}_{g,i} \;=\; \tfrac12\bigl(J_{g,f_L}+J_{g,f_R}\bigr),
\tag{7.14}
$$

evaluated component-wise from the gas-flux matrix populated in Step 3 of the assembly.

**(b) Cell-centered temperature gradient $(\partial T/\partial z)_i$.** For an interior cell
the two one-sided differences are combined with **opposite-distance weights**
(`opposite_distance_gradient`, src/Geometry/geometric_quantities.jl),

$$
\left(\frac{\partial T}{\partial z}\right)_i
=
\frac{d_R\,\dfrac{T_i-T_L}{d_L} + d_L\,\dfrac{T_R-T_i}{d_R}}{d_L+d_R},
\tag{7.15}
$$

where $T_L,T_R$ are the neighbor cell-center temperatures and $d_L,d_R$ the respective
center-to-center distances. Each one-sided slope is a second-order estimate of $T'$ at its
segment midpoint ($z_i - d_L/2$, resp. $z_i + d_R/2$); weighting by the *opposite* distance
interpolates them linearly back to $z_i$, cancelling the leading $T''$ truncation term, so
(7.15) is second-order on arbitrary nonuniform meshes.

A **boundary cell** applies the same interpolation with the *surface* as the missing
neighbor: the solved surface temperature $T_s$ (Step 4 of the assembly) lives on the
boundary face at distance $\Delta z_{1/2} = V_i/(2A_f)$ from the center, giving the
face-side slope $s_{\rm face} = \pm(T_s - T_i)/\Delta z_{1/2}$ alongside the interior
one-sided slope $s_{\rm int}$, and

$$
\left(\frac{\partial T}{\partial z}\right)_i
= \frac{d_{\rm nb}\,s_{\rm face} + \Delta z_{1/2}\,s_{\rm int}}{d_{\rm nb}+\Delta z_{1/2}}
\tag{7.15b}
$$

— second order up to the boundary (exact on quadratics), and sensitive to the
surface-side gradient, where the outflowing gas is hottest, that a plain interior
one-sided difference misses. The implementation recovers the face slope AD-cleanly as
$s_{\rm face} = -q_{\rm total}[f]/k_i$, which is *algebraically* identical because the
thermal-BC pass writes the boundary flux as exactly
$q_{\rm total}[f] = \mp k_i (T_s - T_i)/\Delta z_{1/2}$ at the solved $T_s$ (the $k_i$
cancels). This is why $S_{\rm conv}$ is computed *after* the thermal BC pass. A fully
isolated cell (single-cell mesh, which does not occur in a well-formed problem) gives
$S_{\text{conv},i}=0$.

On a uniform mesh ($d_L=d_R=\Delta z$), (7.15) reduces to the standard central difference
$(T_R-T_L)/(2\Delta z)$. The unweighted mean of one-sided slopes and the wide-stencil form
$(T_R-T_L)/(z_R^c-z_L^c)$ both carry a first-order $\mathcal{O}(|d_R-d_L|)\,T''$ error on
stretched meshes and are not used.

> **Implementation note.** `compute_gas_convective_source!` in `src/Physics/convection.jl`
> builds $\overline{J}_{g,i}$ and the gradient per cell; the inner sum over gas components is
> the `@generated` `_gas_convective_source_unrolled`, which emits, for each gas $j$,
> `S -= cp_j * J_avg * dTdx` with `cp_j = comp.heat_capacity(T_i)` and `J_avg` from (7.14).
> Solid/liquid components are skipped by the compile-time `is_gas` branch. The accumulator is
> initialized to `zero(promote_type(typeof(T_i), typeof(dTdx), eltype(J_gas)))` so it is
> type-stable under AD.

### 7.4.3 The surface carries no gas-enthalpy term

$S_{\text{conv}}$ is always active and is the sole carrier of gas sensible-enthalpy
transport. It is also the reason the surface boundary condition is a *plain*
convective + radiative balance, $q_{\text{BC}}(T_s) = q_{\text{cond}}(T_s)$, with no
gas-enthalpy term: the gas energy leaving the domain has already been accounted for
volumetrically as it advected up to the surface.

Integrating (7.11) over the column shows the bookkeeping. With heats of reaction at the
reaction temperature, the generation piece of (7.11)
is inside $Q_{\text{rxn}}$, and the remaining advection piece, $\int_0^L S_{\text{conv}}\,dz$,
is the heating of the gas from its formation temperature to $T_s$ along its path. Under local
thermal equilibrium the gas crosses the surface at $T_s$ on both sides, so its enthalpy flux
is continuous there and the surface balance is conducted flux equals applied flux (§1.4.3).
A surface term $\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$ added to that balance would apply
the outflow a second time, on top of $\int S_{\text{conv}}\,dz$ — a spurious sink of order
$5$–$10\,\text{kW}\,\text{m}^{-2}$ that under-predicts the surface temperature and the
mass-loss rate. The convention is ThermaKin Eq. 17; Gpyro's thermal-equilibrium form uses the
same term. Parameter sets calibrated against a solver that applies a surface gas-enthalpy term
*in addition to* a volumetric advection source predict lower $T_s$ and MLR there than here;
re-validate imported calibrations (§18).

**The mass-flux/convection coupling is the blowing model of each `ConvectiveBC`.** The
*momentum-side* effect of the gas efflux — boundary-layer thickening that suppresses the
convective exchange coefficient, and the pooling of a dense gas on the face — is independent
of the enthalpy bookkeeping above and is the `blowing` model a `ConvectiveBC` carries
(§12.3.3): the coefficient becomes $h\,\chi$ with $\chi$ from the model and the outward gas
flux at the face. Because it modifies only the coefficient, it composes exactly with the
$S_{\text{conv}}$-only convention; no double-count arises. This is the only use of the gas
flux in the surface balance.

> **Implementation note.** The surface Newton solve is in
> `src/Discretization/boundary_temperature.jl`; it takes the outward gas flux at the face as
> its only gas-flux argument and passes it to the blowing models.

---

## 7.5 Reaction-Heat and Radiation Sources

### 7.5.1 Reaction heat $Q_{\text{rxn}}$ and its sign

The volumetric reaction-heat source is

$$
\boxed{\;
Q_{\text{rxn},i} \;=\; -\sum_r h_r\,r_r
\;}
\qquad [\text{W}\,\text{m}^{-3}]
\tag{7.16}
$$

where $r_r$ is the rate of reaction $r$ (units $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$; §6)
and $h_r$ is its heat of reaction per unit mass of the first reactant (units
$\text{J}\,\text{kg}^{-1}$). The **internal storage sign convention** (§2, overload note H1)
is

- $h_r > 0$ **endothermic** (the reaction *absorbs* heat). Then $Q_{\text{rxn}}<0$: the cell
  **cools**.
- $h_r < 0$ **exothermic** (the reaction *releases* heat). Then $Q_{\text{rxn}}>0$: the cell
  **heats**.

The leading minus sign in (7.16) converts the enthalpy-of-reaction convention into the
heat-source convention. This sign is the most error-prone in the whole assembly, so it is
worth stating precisely against the implementation: the accumulator computes
`Q -= heats[i] * rates[i]` — i.e. $Q=-\sum_r h_r r_r$ — with rates $r_r\ge 0$ for forward
reactions, so an exothermic reaction ($h_r<0$) produces $Q>0$.

> **Caution on external conventions.** ThermaKin and Gpyro publish $h>0=$ *exothermic*, the
> opposite of the internal storage convention here. The mapping is absorbed by the minus sign
> in (7.16): with their published $h$ one would write $Q_{\text{rxn}}=+\sum h r$.
> Pyrolysis.jl stores $h$ with $h>0=$ endothermic and forms $Q_{\text{rxn}}=-\sum h r$, so a
> user porting a ThermaKin heat of reaction must flip its sign (see §6).

> **Implementation note.** `reaction_heat_source(rates, heats)` in `src/Physics/kinetics.jl`
> performs the negated sum. The per-cell reaction rates and heats are evaluated by
> `evaluate_reactions!` (§6) at the cell's local $(T_i,\xi_i)$, and the stoichiometric species
> sources by `species_source_terms!`; `compute_rhs!` evaluates them inline, cell by cell.

### 7.5.2 Radiation source $Q_{\text{rad}}$ (dispatch only)

The volumetric radiation absorption source is computed by `compute_radiation_source!` and
dispatched on the active radiation model:

- **`NO_RADIATION`**: $Q_{\text{rad},i}=0$.
- **`SURFACE_ABSORPTION`**: $Q_{\text{rad},i}=0$ everywhere — the absorbed flux
  $\alpha\,q_{\mathrm{ext}}$ is applied through the surface energy balance
  (boundary Newton, §12.3.5), not as an $I/\Delta z$ volumetric deposit (opaque
  surface, better physics: the absorbed power raises $T_s$ directly).
- **`BEER_LAMBERT`**: in-depth absorption, $q'''(z)=\alpha(z)I(z)$ with
  $I(z)=I_{\text{surface}}\exp(-\!\int\!\alpha\,dz')$, assembled cell-exact via a backward
  sweep from the surface.
- **`P1_QUASI_STEADY`**: $Q_{\text{rad},i}=0$ here; P1 lives in `Pyrolysis.Experimental` with its own residual wrapper, and `solve` rejects it (§8.12).

The full derivation of each model — Beer–Lambert cell-averaged form, P1 diffusion, Marshak
boundary conditions, the energy-equation coupling, and the regime diagnostics — is the subject
of **§8 (Thermal Radiation)**. For the purposes of the energy assembly here, $Q_{\text{rad},i}$
is simply a cell-centered source term ($\text{W}\,\text{m}^{-3}$) added in (7.9). Note that the
boundary-incident radiative flux can also be applied through the *surface energy balance* (the
`RadiativeFluxBC` / boundary Newton solve, §12); which path a given model takes is part of §8's
coupling discussion.

---

## 7.6 The `compute_rhs!` Eight-Step Assembly Sequence

The complete energy + species RHS is assembled in `compute_rhs!`
(`src/Discretization/finite_volume.jl`) in a **fixed eight-step order**. The ordering is not
arbitrary: several steps consume results produced by earlier steps, and getting the order
wrong would break the surface energy balance, the boundary-cell $S_{\text{conv}}$, or the heat
divergence. One function serves
both transport laws: the law follows the property cache's Darcy sub-cache (`nothing` ⇔
Fickian), which selects the Darcy property refresh in Step 1 and the gas-flux form in Step 3
(§9). The sequence is:

**Step 1 — Transport properties and interior conductive fluxes.** In Darcy–Fick mode
`update_darcy_properties!` refreshes the cell pressure, porosity, permeability, and viscosity
and the face Darcy velocity (§9); in Fickian mode this is a no-op. `compute_conductive_fluxes!`
then fills $q_f$ at all interior faces via (7.1) and zeroes boundary faces (7.3). The result is
copied into the working flux vector `q_total` (`copyto!`), which will be modified at the
boundaries in Step 4.

**Step 2 — Reaction rates, heat source, and species sources.** For each cell,
`evaluate_reactions!` gives the rates $r_r$ and heats $h_r$; `reaction_heat_source` forms
$Q_{\text{rxn},i}$ (7.16); `species_source_terms!` forms $S_{j,\text{rxn},i}$, which is stored
in `dξ_rxn`. Step 8 reads both; in the ALE residual the strain rate and the mesh velocity
are built from `dξ_rxn` after `compute_rhs!` returns (§10.1.2, §11.3).

**Step 3 — Gas transport fluxes (interior + boundary).** `compute_gas_fluxes!` (Fickian
volume-fraction-gradient form) or `compute_gas_fluxes_darcy_fick!` (Darcy–Fick, with the
pressure-BC venting flux on the boundary faces) fills the gas-flux matrix `J_gas` (§9), then
`apply_mass_bc_flux!` adds the mass boundary conditions at the top and bottom faces. *This
must complete before the thermal-BC solve* (Step 4): the blowing model of each
`ConvectiveBC` reads the outward gas flux at its face from `J_gas` (§7.4.3, §12.3.3).

**Step 4 — Thermal boundary conditions.** `apply_thermal_bc_flux!` solves, at each
non-Dirichlet boundary, the surface energy balance for the surface temperature $T_s$ by Newton
iteration, and **overwrites** the boundary entries of `q_total` with the resulting conductive
flux $k_{\text{eff}}(T_s-T_{\text{cell}})/\Delta z_{\text{half}}$ (with the boundary sign
convention: $+$ at the bottom, $-$ at the top, so the divergence sees the correct
inflow/outflow). Dirichlet BCs set $T_s$ directly and evaluate the flux without iteration. The
full surface-balance derivation — the Newton residual $f(T_s)=q_{\text{BC}}(T_s)-
k_{\text{eff}}(T_s-T_{\text{cell}})/\Delta z_{\text{half}}$, its derivative, the
$[200,3000]\,\text{K}$ clamp, and the convergence criterion — is in §12.

**Step 5 — Radiation source.** `compute_radiation_source!` fills $Q_{\text{rad}}$ per the
active model (§7.5.2, §8).

**Step 6 — Gas convective source.** `compute_gas_convective_source!` fills $S_{\text{conv}}$
(7.13) from `J_gas`, the cell-centered temperature gradients, and — at the two boundary
cells — the boundary entries of `q_total` written in Step 4, from which the boundary-face
slope at the solved $T_s$ is recovered as $-q_{\rm total}[f]/k_i$ (7.15b). (It uses the
full face-averaged gas fluxes, including the boundary-applied values, now in `J_gas`.)

**Step 7 — Divergences.** `divergence_1d!` forms $(\nabla\!\cdot\mathbf{q})_i$ from
`q_total` (now including the boundary fluxes from Step 4), and `divergence_species_1d!` forms
$(\nabla\!\cdot\mathbf{J}_j)_i$.

**Step 8 — Final RHS assembly.** For each cell, the temperature derivative (7.9) is formed by
combining $-(\nabla\!\cdot\mathbf{q})_i$, $S_{\text{conv},i}$,
$Q_{\text{rad},i}$, $Q_{\text{rxn},i}$ and dividing by `prop_cache.ρcp[i]`; the species
derivative (7.10) is formed by `apply_species_residual_unrolled!`.

### 7.6.1 The ordering dependency graph

The dependencies that *fix* the order are:

1. **Gas fluxes (Step 3) → thermal BC solve (Step 4):** the blowing models of the
   convective BCs need the outward gas flux at the boundary.
2. **Thermal BC solve (Step 4) → gas convective source (Step 6) and heat divergence
   (Step 7):** both read the boundary fluxes that Step 4 wrote into `q_total` — the
   divergence directly, and $S_{\text{conv}}$ through the boundary cells' face slope
   $-q_{\rm total}[f]/k_i$ at the solved $T_s$ (7.15b).
3. **Reaction sources (Step 2) → everything downstream:** $Q_{\text{rxn}}$ and the species
   reaction sources are needed in Step 8; in the *ALE* residual (§11, §15) the reaction species
   sources also feed the strain rate $\theta_i$ and hence the cumulative mesh velocity $w$ — so
   reaction rates must be available *before* the mesh velocity is computed.

The full-residual entry point `residual!`, its two trait-dispatched paths (Eulerian and ALE),
and the ALE-specific extensions to this sequence (strain rate, mesh velocity, dilation
correction, $\chi$ dynamics) are described in §15
(Temporal Integration & Structured Jacobian) and §11 (ALE Framework). This chapter's
`compute_rhs!` is the Eulerian core that those paths build on.

---

## 7.7 Caches and Type Parameterization for AD

### 7.7.1 `RHSCache{Tu}`

All RHS intermediates are stored in a pre-allocated `RHSCache{Tu}`
(`src/Discretization/finite_volume.jl`), eliminating per-call allocations (verified by the
zero-allocation RHS test). Its fields:

| field | shape | meaning | units |
|---|---|---|---|
| `q_cond` | $n_{\text{faces}}$ | conductive face fluxes | $\text{W}\,\text{m}^{-2}$ |
| `q_total` | $n_{\text{faces}}$ | total energy face fluxes ($q_{\text{cond}}$ + BC) | $\text{W}\,\text{m}^{-2}$ |
| `S_conv` | $n_{\text{cells}}$ | gas convective source | $\text{W}\,\text{m}^{-3}$ |
| `J_gas` | $N_C\times n_{\text{faces}}$ | gas species face fluxes | $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$ |
| `Q_rxn` | $n_{\text{cells}}$ | reaction heat source | $\text{W}\,\text{m}^{-3}$ |
| `Q_rad` | $n_{\text{cells}}$ | radiation source | $\text{W}\,\text{m}^{-3}$ |
| `dξ_rxn` | $N_C\times n_{\text{cells}}$ | species reaction sources | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |
| `rates`, `heats` | $N_R$ | per-cell reaction workspace | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$, $\text{J}\,\text{kg}^{-1}$ |
| `species_src` | $N_C$ | per-cell species-source workspace | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |
| `div_q` | $n_{\text{cells}}$ | $\nabla\!\cdot\mathbf{q}$ | $\text{W}\,\text{m}^{-3}$ |
| `div_J` | $N_C\times n_{\text{cells}}$ | $\nabla\!\cdot\mathbf{J}_j$ | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |

The cache keeps a fixed size for the whole solve: depleted-cell parking (§14.3) leaves
ghost slots in place instead of changing the cell/face counts, and `compute_rhs!` loops
over the live cells `1:n_active_cells` only. It is constructed sized to a mesh and material
via `RHSCache(mesh, material)`, defaulting to `Tu = Float64`.

### 7.7.2 The scalar type parameter `Tu`

`RHSCache{Tu}`, `PropertyCache{Tu,DarcyT}`, and the typed state cache are all parameterized on
a single live scalar type `Tu`. In ordinary integration `Tu = Float64`. During forward-mode
automatic differentiation — used by the structured-Jacobian operators' stencil-local
ForwardDiff kernels and by the dense reference Jacobian of the verifier —
`Tu = ForwardDiff.Dual`. Because every flux, source, property, and geometric quantity is held
at type `Tu`, the entire RHS assembly propagates dual sensitivities without a single
re-allocation or type conversion. The compile-time component unrolling (`@generated` functions
for the heat-capacity sums, $S_{\text{conv}}$, and the species residual) is
what makes this type stability achievable over the heterogeneous component tuple: the gas/solid
branches are resolved at code generation, so there is no runtime tuple indexing (which would
allocate) and the compiler specializes each arithmetic line to the concrete `Tu`.

> **Implementation note.** The production assembly uses `RHSCache` exclusively. The
> `PropertyCache` (§5) supplies `ρ`, `ρcp`,
> `k`, `λ`, `ε`, `a`, `α` per cell and `k_face`, `λ_face` per face (plus an optional Darcy sub-cache
> of `P`, `φ`, `κ`, `μ` per cell and `v_gas`, `∇P` per face); it must be refreshed by `update_properties!` before `compute_rhs!` is
> called.

---

## 7.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to the unified nomenclature on first use and contrast the
energy assembly. (The reference codes are condensed-phase pyrolysis / fire solvers:
ThermaKin2Ds by Stoliarov et al.; Gpyro by Lautenberger & Fernandez-Pello; the FDS solid-phase
model by McGrattan et al.)

### 7.8.1 The energy-equation skeleton

All four codes solve a 1D (through-thickness) Fourier-conduction energy balance with
volumetric reaction heat and a gas-enthalpy coupling. In our symbols the common skeleton is

$$
(\rho c_p)\,\frac{\partial T}{\partial t}
=
\frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right)
+ Q_{\text{rxn}} + Q_{\text{rad}} + (\text{gas enthalpy term}),
$$

and the codes differ mainly in (i) whether gas heat capacity is in the storage term, (ii) how
the gas enthalpy coupling is written, and (iii) the radiation closure.

### 7.8.2 ThermaKin2Ds

ThermaKin writes the through-thickness energy equation with the all-component storage
$\sum_j \xi_j c_{p,j}\,\partial T/\partial t$, the gas-advection term
$c_g\,(\mathbf{J}_g\!\cdot\!\partial T/\partial x)$ (its Eq. 17) — *exactly* the
$S_{\text{conv}}$ of (7.12) — and a deformation term that is the skeleton-following part
of the time derivative, which the ALE mesh frame absorbs (§11.2.2). Pyrolysis.jl's energy
equation is therefore ThermaKin Eq. 17 term for term. ThermaKin maps to our symbols as:
ThermaKin's stoichiometric $\theta_i^{\,j}\to \nu_{j,r}$ (the strain rate $\theta$ is
reserved here for dilation; §2 note A2); ThermaKin's bimolecular concentration product
$\xi_{\text{COMP1}}\xi_{\text{COMP2}}$ has no analog because Pyrolysis.jl enforces a single
reactant (§6); ThermaKin's gas flux $J=-\rho_g\lambda\,\partial(\xi_g/\rho_g)/\partial x$ is
the volume-fraction-gradient form implemented here (§9).

### 7.8.3 Gpyro

Gpyro's energy equation (Lautenberger) uses a control-volume conservative formulation on
$\bar\rho\,\Delta z$, with the bulk (averaged) density $\bar\rho$ (our $\rho$),
volume-fraction ($X_\alpha\to v_j$) property averaging, and power-law temperature-dependent
properties. Its thermal-equilibrium form carries the gas through the solid-to-gas heat
transfer $\dot Q_{s-g}''' \approx \dot m'' c_{pg}\,\partial T/\partial z$ (its Eq. 65) — our
$S_{\text{conv}}$ with the Darcy mass flux — and drops the gas storage
$\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default; its `FULL_QSG` option restores it,
which makes the equation (7.6) with (7.7). The mapping: Gpyro's pre-exponential
$Z\to A_r$; its conversion variable $\alpha_{\text{conv}}$ has no single analog
(Pyrolysis.jl tracks mass concentrations $\xi_j$ directly, with $\chi$ as a separate
progress variable for cross-section, §10); its permeability $K\to\kappa$; its mass flux
$\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ is the Darcy term of §9 (Pyrolysis.jl
omits the gravity/buoyancy term $g\rho_g$ as negligible in the thin condensed phase).
Numerically Gpyro uses Patankar control-volume discretization with Peclet-based scheme
selection, harmonic interface interpolation, and a TDMA solve; Pyrolysis.jl uses the same
harmonic interface interpolation (7.2) but a method-of-lines + stiff-ODE integrator with a
structured Jacobian (§15) rather than a hand-rolled implicit TDMA sweep.

### 7.8.4 FDS

The FDS solid-phase model solves a 1D Fourier conduction equation with the surface balance
$-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$, Arrhenius reaction kinetics
($\rho^n O_2^n T^n$ rate form), and no in-depth gas: volatiles leave instantaneously, so
there is no gas storage and no in-depth gas-enthalpy coupling; the gas enthalpy is closed
at the surface. FDS maps as: pre-exponential $Z\to A_r$; solid absorption coefficient
$\kappa_s\to\alpha$; the surface energy balance with $\dot q_c''$ (convective) and
$\dot q_r''$ (radiative) is the boundary Newton solve of §12. The radiation closure differs
(FDS uses a two-flux Schuster–Schwarzschild model in-depth; Pyrolysis.jl offers
SURFACE_ABSORPTION / BEER_LAMBERT / P1, §8). Pyrolysis.jl, following ThermaKin, resolves
the gas in depth and distributes its enthalpy transport volumetrically as $S_{\text{conv}}$
rather than as a surface outflow term, and the surface balance carries no gas-enthalpy term
on top of it (§7.4.3).

### 7.8.5 Summary table

| feature | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| storage $\rho c_p$ | all components (7.7) | all components (incl. gas) | matrix by default; all with `FULL_QSG` | matrix (no in-depth gas) |
| gas enthalpy coupling | $S_{\text{conv}}$ volumetric | $S_{\text{conv}}$ volumetric (Eq. 17) | $\dot m'' c_{pg}\,\partial T/\partial z$ (thermal equilibrium) | surface outflow |
| face conductivity | distance-weighted harmonic mean | harmonic | harmonic (Patankar) | — |
| conduction discretization | cell-centered FV, MOL + stiff ODE | Crank–Nicolson (x) | Patankar CV + TDMA | implicit FD |
| radiation | none / surface / Beer–Lambert / P1 | Beer–Lambert + re-radiation | in-depth exp absorption | two-flux S–S |
| heat-of-reaction sign | $h>0$ endothermic, $Q_{\text{rxn}}=-\sum hr$ | $h>0$ exothermic (flip) | $h>0$ exothermic (flip) | per-reaction |

---

## 7.9 Limitations and Open Issues

The energy assembly inherits a small number of known limitations, drawn from the source notes
and discussed where relevant above:

1. **One temperature.** The pore gas is stored at the matrix temperature (§7.3.2); there is
   no two-temperature gas energy equation, so regimes where the gas does not equilibrate
   with the pore walls within a cell are outside the model (§3.9).

2. **One global radiation model.** $Q_{\text{rad}}$ is dispatched on a single model for the
   whole domain; per-cell or temperature-dependent model switching is not supported (§8).

3. **Discrete energy-balance closure for composition-varying mixtures.** The energy ledger
   (§16) audits bookkeeping closure to integrator tolerance but does not carry the
   composition-change enthalpy $\sum_j h_j(T)\,\partial_t\xi_j$ as a conserved quantity
   (§3.3.2) or the ALE volume-change work, so the *physical* balance (as opposed to the
   discrete-bookkeeping balance) carries those as part of the quantified approximation cost
   rather than as an exactly-closed term.

These do not affect the correctness of the discrete update (7.9)–(7.10); they bound the
*modeling* fidelity of the energy equation and are exposed, not hidden, through the
conservation diagnostics of §16.
