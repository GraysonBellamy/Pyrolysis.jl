# 7. Heat Conduction, Gas Convection, and Energy Assembly

This chapter derives the condensed-phase energy balance as Pyrolysis.jl assembles
it: a cell balance of sensible enthalpy in flux form, advanced in one of two energy
states (the cell temperature by default, the cell's enthalpy inventory as the
alternative), in which the storage sums the heat capacity of every component, pore
gas included, and the sensible enthalpy carried by the gas flux relative to the
skeleton crosses the cell faces together with that flux. We develop each term from
its continuous form to its finite-volume discretization — Fourier conduction with
the guarded distance-weighted harmonic-mean face conductivity, the face scatter that
forms the discrete divergence, the gas enthalpy flux and its face temperatures, the
reaction-heat source $Q_{\text{rxn}}$, and the radiation source $Q_{\text{rad}}$ —
and we trace the residual sequence that produces the energy and mass rows. Particular
care is taken with the sign conventions of §2 (Nomenclature), with the precise set of
terms the code includes and excludes, and with why the surface energy balance carries
no gas-enthalpy term. Where useful, the Pyrolysis.jl forms are mapped to and
contrasted with ThermaKin2Ds, Gpyro, and FDS.

Index conventions for this chapter follow §2 overload note I1: **$i$ is the cell index**
($1\le i\le n$) and **$r$ is the reaction index**; the component/species index is $j$
($1\le j\le N_C$). The spatial coordinate is $z$ with $z=0$ the bottom/substrate
(boundary tag `:bottom`) and $z=L$ the top/exposed surface (tag `:top`); heat enters
from the top. Every face flux $F_f$ is a transport per unit face area **along the face
normal**, which points from the face's owner cell to its neighbour and out of the
domain on a boundary face: on a line the interior faces and the `:top` face have their
normal along $+z$, the `:bottom` face along $-z$. The orientation of face $f$ relative
to cell $i$ is $s_{if} = +1$ when the normal points out of $i$ (cell $i$ owns the face)
and $-1$ otherwise.

> **Frame of this chapter.** The balance developed here is written on cells that ride
> the skeleton — a fixed mesh, or a moving mesh without a correction. The condensed
> phase carrying $\rho c_p^{\text{eff}}$ is then at rest relative to the cell, and the
> energy row has no solid-convection term: the mesh-frame rate of the matrix-carried
> temperature is its material rate (§11.2.2). Where a correction moves the mesh
> relative to the solid, every face additionally carries the inventories at the
> relative velocity and the sensible enthalpy they carry (§11.6); that flux enters the
> same energy row as one more face energy flux (§7.3.3, §7.4.2).

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
Properties & Mixing Rules). The state mapping evaluates it once per cell, at the cell's
positive-part concentrations and trial temperature, and floors it smoothly at
$10^{-6}$ W·m⁻¹·K⁻¹ (§5.5, §5.13); only its *face value* enters here.

### 7.1.2 Discrete face flux: distance-weighted harmonic mean

In the cell-centered finite-volume layout (§13, Finite-Volume Discretization) the
temperature $T$ is stored at cell centers $z_i^c$, and fluxes are required at the faces
that separate adjacent cells. For an interior face $f$ with owner cell $o$ and neighbour
cell $n$ (on a line, the cell below the face and the cell above it), cell temperatures
$T_o, T_n$ and center-to-face distances $d_o, d_n$ along the normal, the conductive flux
is

$$
\boxed{\;
q_f \;=\; -\,k_f\,\frac{T_n-T_o}{d_o+d_n}
\;}
\tag{7.1}
$$

where the **face conductivity** $k_f$ is the **distance-weighted harmonic mean** of the two
cell conductivities, in the guarded form of §5.12:

$$
\boxed{\;
k_f \;=\; \frac{(d_o+d_n)\,k_o\,k_n}{d_o\,k_n + d_n\,k_o + \eta\,(d_o+d_n)},
\qquad \eta = 10^{-20}\ \text{W·m}^{-1}\text{K}^{-1}
\;}
\tag{7.2}
$$

Without the guard, (7.2) is the resistance-in-series mean
$(d_o+d_n)/(d_o/k_o + d_n/k_n)$, to which it is equal to a relative $\eta/k_f$. The
denominator $d_o+d_n$ of (7.1) is the center-to-center distance along the normal; on a
line it is $|z_n^c - z_o^c|$.

**Derivation of the harmonic mean from flux continuity.** Treat the half-cell on either side
of the face as a thermal resistance. With face area $A_f$ the resistance of the owner's
half is $R_o=d_o/(k_o A_f)$ and of the neighbour's half $R_n=d_n/(k_n A_f)$. At steady state
the same flux $q_f$ passes through both half-cells in series, so the temperature drop from
$T_o$ to $T_n$ across the combined resistance is

$$
T_o-T_n \;=\; q_f\,A_f\,(R_o+R_n)
\;=\; q_f\left(\frac{d_o}{k_o}+\frac{d_n}{k_n}\right).
$$

Solving for $q_f$ and writing it as an effective conductivity $k_f$ acting over the total
distance $d_o+d_n$,

$$
q_f \;=\; -\,\frac{T_n-T_o}{\dfrac{d_o}{k_o}+\dfrac{d_n}{k_n}}
\;\equiv\; -\,k_f\,\frac{T_n-T_o}{d_o+d_n},
\qquad
k_f=\frac{d_o+d_n}{\dfrac{d_o}{k_o}+\dfrac{d_n}{k_n}},
$$

which is (7.1)–(7.2) without the guard. The harmonic (resistance-in-series) mean is the
*only* face averaging that makes the discrete flux exact for a piecewise-constant
conductivity field (a sharp material interface at the face): an arithmetic mean would
over-predict the flux through the more-resistive layer. The same series-resistance
averaging is applied to the gas-transfer coefficient $\lambda$ and, with a smaller guard,
to the Darcy mobility $\kappa/\mu$ (§5.12, §9).

**Uniform-mesh limit.** When the cell centers are equidistant from the face, $d_o=d_n=d$,
(7.2) collapses to the standard (unweighted) harmonic mean,

$$
k_f \;\xrightarrow{\,d_o=d_n\,}\; \frac{2k_ok_n}{k_o+k_n}
$$

(up to the guard). On **stretched** meshes — for example, boundary-refined grids used to
resolve a steep near-surface thermal gradient — $d_o\ne d_n$, and the distance weighting
keeps the steady flux through a sharp interface at the face exact, at no extra cost.

**Vanishing conductivity.** The guard makes (7.2) go to zero continuously, with a finite
derivative, as either cell's conductivity vanishes, instead of cutting the face off below a
threshold; the cell conductivity itself is floored at $10^{-6}$ W·m⁻¹·K⁻¹ in the state
mapping (§5.13), so the guard only matters on trial states.

The kernel evaluates (7.1) at the cell temperatures $T_i$ themselves, not at their trial
clamp, and with the center-to-face distances of the current geometry, each smoothly floored
at $10^{-8}$ of its reference value so that an inverted trial mesh keeps a finite flux.

> **Implementation note.** The interior face fluxes are added by `ConductionKernel`
> (`src/Kernels/conduction.jl`) to the face energy fluxes of the contributions (§7.7),
> with the face conductivity from `guarded_harmonic_mean`
> (`src/Properties/effective_properties.jl`) and the distance floors
> `DISTANCE_FLOOR_FRACTION = 1e-8`.

### 7.1.3 Interior faces and boundary faces

The conduction kernel visits **only the interior faces between two live cells**. A face is
skipped when it is a boundary face (no neighbour) or when either of its cells has been
deactivated (§14). The heat flux through a boundary face is written by a different kernel,
the thermal boundary kernel (`ThermalBoundaryKernel`, `src/Kernels/thermal_boundary.jl`):
for each boundary tag it applies the tag's condition on the face that bounds the live cells
on that side (`live_face`; the boundary face itself while no cell beyond it has been
deactivated, otherwise the face of the outermost live cell), and it writes that face's
energy flux and its face temperature from the surface energy balance at the solved surface
temperature $T_s$ (§7.4.3, §12). Every kernel *adds* to the shared face fluxes; no kernel
overwrites another's, and a face beyond the live cells carries nothing. The interior
conduction and the boundary conditions are thus separated by face, and the colored Jacobian
(§15) differentiates both with the rest of the residual.

On a mesh of more than one dimension (residual level only, not a supported path for full
solves; §13) the conductivity of a cell is a tensor $k_i\mathbf A_i$ with a dimensionless
symmetric positive-definite $\mathbf A_i$, and the face tensor is
$\mathbf K_f = k_f\,(d_n\mathbf A_o + d_o\mathbf A_n)/(d_o + d_n)$. The flux along the normal
$\hat{\mathbf n}$ splits into a two-point normal part and a cross-diffusion part,

$$
q_f = -(\hat{\mathbf n}\cdot\mathbf K_f\hat{\mathbf n})\,\frac{T_n - T_o}{d_o + d_n}
- \mathbf t\cdot\overline{\nabla T}_f ,
\qquad
\mathbf t = \mathbf K_f\hat{\mathbf n} - (\hat{\mathbf n}\cdot\mathbf K_f\hat{\mathbf n})\,\frac{\mathbf d}{d_o + d_n},
$$

with $\mathbf d$ the vector between the cell centers and $\overline{\nabla T}_f$ the
distance-weighted mean of the two cells' least-squares gradients over their live face
neighbours, so a linear temperature field's flux is exact wherever the neighbours span every
direction; where they do not, a Tikhonov term ($10^{-8}$ of the trace of the normal matrix)
makes the unresolved component of the gradient zero. On a line $\mathbf d$ is along the
normal and the cross part vanishes identically. Element death on such a mesh exposes the
faces between live and dead cells, which take the condition of the exposed tag (§14.5).

---

## 7.2 The Discrete Divergence Operator

The face fluxes enter the balances through their divergence. Integrating
$\nabla\!\cdot\mathbf{F}$ over cell $i$ and applying the divergence theorem gives the net
outflow through the bounding faces,

$$
\boxed{\;
\int_{\Omega_i}\nabla\!\cdot\mathbf F\,dV \;=\; \sum_{f\in\partial i} s_{if}\,A_f\,F_f
\;}
\tag{7.3}
$$

with $A_f$ the area of face $f$ in the current geometry ($A$ on a planar line; $2\pi rA$ and
$Ar^2$ on cylindrical and spherical meshes, §13) and $(1/V_i)\sum_f s_{if}A_fF_f$ the
cell-averaged divergence. The rows are **extensive**: because the integrated states are the
cell masses and the temperature row divides by the extensive heat-capacity sum $C_i$ (§7.3.3),
the assembly forms the face sums (7.3) themselves and never divides by a volume.

**The scatter.** The row assembly (`assemble_rows!`, `src/Assembly/residual.jl`) loops over the
faces once: each face adds $-F_f A_f$ to its owner's row and $+F_f A_f$ to its neighbour's (a
boundary face has no neighbour). The energy row of cell $i$ thus receives
$-\sum_f s_{if}A_f F^E_f$ from its face energy fluxes (§7.3.3), and the mass rows receive the
**face-flux part of the mass rate**

$$
\dot M_{j,i} \;=\; -\sum_{f\in\partial i} s_{if}\,A_f\,F_{j,f}
\qquad [\text{kg}\,\text{s}^{-1}],
\qquad j=1,\dots,N_C ,
\tag{7.4}
$$

which the temperature row also reads (§7.3.3).

**Sign convention (outflow-positive).** With the normal pointing out of the owner, a positive
$F_f$ carries the quantity out of the owner and into the neighbour. On a line, cell $i$ owns
its top face and is the neighbour of its bottom face, so (7.3) reads
$A\,(F_{\text{top}} - F_{\text{bot}})$: positive when there is a *net loss* from the cell.
Because the rows carry the face sums with a minus sign, a net outflow correctly *reduces*
$dT_i/dt$ and $dm_{j,i}/dt$: a cell losing more heat than it gains cools.

**Conservation.** Every interior face flux is added to one cell and subtracted from the other,
so the sum of the rows over the live cells telescopes to the boundary faces alone: a constant
flux field gives zero net rate in every cell, and
$\sum_i \int_{\Omega_i}\nabla\!\cdot\mathbf F\,dV = \sum_{f\in\text{boundary}}A_fF_f$
recovers the global balance exactly. This is the discrete conservation property the ledgers
of §16 rest on (`test/assembly/test_residual_jacobian.jl`, "row assembly conserves mass and
scatters fluxes antisymmetrically"). The face areas and the cell geometry are those of the
current state, so on a moving mesh the scatter picks up the geometric change of every cell
(§11.5).

---

## 7.3 Effective Heat Capacity and the Energy Equation

### 7.3.1 The energy balance and its temperature form

The quantity Pyrolysis.jl conserves is the **sensible enthalpy** of the material,
$\sum_j \xi_j h_j(T)$, with $h_j(T) = \int_{T_{\text{ref}}}^{T} c_{p,j}\,dT'$ on the common
datum $T_{\text{ref}} = 298.15$ K (§5.4.3). On cells that ride the skeleton its balance is
(§1.4.1, §3.3)

$$
\boxed{\;
\frac{\partial}{\partial t}\Big(\sum_j \xi_j h_j\Big)
+ \nabla\!\cdot\Big(\mathbf q_{\text{cond}} + \sum_g h_g\,\mathbf J_g\Big)
\;=\; Q_{\text{rad}} + Q_{\text{rxn}} + \sum_j S_j\,h_j
\;}
\tag{7.5}
$$

with $\mathbf J_g$ the mass flux of gas component $g$ relative to the skeleton,
$Q_{\text{rxn}} = -\sum_r h_r r_r$ the reaction heat with the heats quoted at the local
temperature (§6.6), and $S_j = \sum_r \nu_{r,j} r_r$ the reaction source of component $j$;
the last term converts the heat of reaction to the sensible-enthalpy datum (§6.6.3).
Subtracting $h_j$ times each species balance
$\partial_t\xi_j + \nabla\!\cdot\mathbf J_j = S_j$ ($\mathbf J_j \equiv 0$ for condensed
components) and using $\nabla h_g = c_{p,g}\nabla T$ gives the **temperature form**

$$
\boxed{\;
\rho c_{p}^{\text{eff}}(T,\xi)\,\frac{\partial T}{\partial t}
\;=\;
-\,\nabla\!\cdot\mathbf{q}_{\text{cond}}
\;-\; \sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T
\;+\; Q_{\text{rad}}
\;+\; Q_{\text{rxn}}
\;}
\tag{7.6}
$$

with, term by term,

| term | expression | role | § |
|---|---|---|---|
| $\mathbf{q}_{\text{cond}}$ | $-k\,\nabla T$ | Fourier conduction | 7.1 |
| $-\sum_g c_{p,g}\,\mathbf J_g\!\cdot\!\nabla T$ | temperature-form reading of $-\nabla\!\cdot(\sum_g h_g\mathbf J_g)$ | gas sensible-enthalpy transport | 7.4 |
| $Q_{\text{rad}}$ | model-dependent | volumetric radiation absorption | 7.5, §8 |
| $Q_{\text{rxn}}$ | $-\sum_r h_r r_r$ | reaction heat | 7.5 |

The subscript $g$ runs over the **gas** components only. The storage coefficient
$\rho c_p^{\text{eff}}$ is the **mixture** volumetric heat capacity, defined next. The code
does not assemble (7.6) pointwise: it assembles the cell integral of (7.5) in flux form, and
its temperature row is the cell integral of (7.6) obtained from it exactly (§7.3.3).

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

and its discrete, extensive form is the cell's **heat-capacity sum**
$C_i = \sum_j m_{j,i}\,c_{p,j}(T_{\text{trial},i})$, smoothly floored at
$\rho c_{p,\text{floor}}V_i(0)$ with $\rho c_{p,\text{floor}} = 10^3$ J·m⁻³·K⁻¹ (§5.4.1,
equation (5.10)), so that a cell emptied of all mass keeps a finite rate.

The pore gas is at the matrix temperature (one-temperature model, §3.1.2) and its mass is an
integrated state, so its storage term costs nothing and is kept. Dropping it — the
*quasi-steady-gas approximation* of Gpyro's default thermal-equilibrium form (restored there
by `FULL_QSG`) and the implicit convention of FDS, which has no in-depth gas — would be an
approximation with no benefit. At $1\,\text{atm}$ the in-pore gas heat capacity is a small
fraction of the matrix heat capacity, from a few hundredths of a percent to about half a
percent (the estimate is in §5.4.2), scaling with pressure.

The gas enters the energy equation in two places and two places only: its storage in (7.7)
and the sensible enthalpy its *flux relative to the skeleton* carries across the faces
(§7.4). The derivation of §3.3.2 shows these two terms, plus the reaction heat with heats
quoted at the reaction temperature, to be the complete equation — ThermaKin Eq. 17 term for
term. On a moving mesh the whole $\rho c_p^{\text{eff}}$ rides the mesh with the skeleton, so
there is no mesh-motion term in the energy equation (§11.2.2, §11.2.4).

> **Implementation note.** $C_i$ is the field `cp_sum` of the derived state, evaluated by the
> state mapping (`src/States/derived.jl`) as the unrolled sum `heat_capacity_sum` over the
> component tuple with no phase branch (`src/States/caloric.jl`). The same masses and
> sensible enthalpies define the energy inventory of the ledger (§16.4), and a deactivation
> event merges cells at their summed sensible enthalpy (§14.3.3).

### 7.3.3 The discrete energy row

Integrating (7.5) over cell $i$ with the divergence theorem (7.3) gives the cell's
**energy right-hand side**

$$
E_i \;=\; S_{E,i} \;-\; \sum_{f\in\partial i} s_{if}\,A_f\,F^{E}_{f},
\qquad
F^{E}_{f} \;=\; q_f + \sum_j F_{j,f}\,h_j(T_f)
\quad [\text{W}\,\text{m}^{-2}],
\tag{7.8}
$$

with $S_{E,i}$ [W] the cell's energy source — the reaction heat $Q_{\text{rxn},i}V_i$ and the
absorbed radiation (§7.5) — and $F^E_f$ the face's **energy flux**: conduction plus the
sensible enthalpy its mass fluxes $F_{j,f}$ carry at the face temperature $T_f$ (§7.4.2). The
energy row is then written in one of two energy states (`energy_state`):

$$
\textbf{temperature form:}\quad
\frac{dT_i}{dt} \;=\; \frac{E_i - \sum_j h_j(T_i)\,\dot M_{j,i}}{C_i},
\tag{7.9}
$$

$$
\textbf{inventory form:}\quad
\frac{dH_i}{dt} \;=\; E_i + \sum_j \omega_{j,i}\,V_i\,h_j(T_i),
\tag{7.10}
$$

with $\dot M_{j,i}$ the face-flux part of the mass rate (7.4), $\omega_{j,i}V_i$ the reaction
mass source [kg·s⁻¹] and $H_i$ the cell's enthalpy inventory (the extended caloric relation
(5.12), equal to $\sum_j m_{j,i}h_j(T_i)$ away from empty cells). Substituting (7.4) and (7.8)
into (7.9),

$$
C_i\,\frac{dT_i}{dt}
\;=\; S_{E,i}
\;-\; \sum_{f\in\partial i} s_{if}\,A_f\Big[q_f + \sum_j F_{j,f}\big(h_j(T_f) - h_j(T_i)\big)\Big],
\tag{7.11}
$$

the cell integral of (7.6): the conductive divergence, the gas enthalpy transport as the face
enthalpy relative to the cell's own (§7.4.2), and the sources. The two rows are one balance:
with the mass rows $dm_{j,i}/dt = \omega_{j,i}V_i + \dot M_{j,i}$ (§7.3.4) the product rule
gives $dH_i/dt = C_i\,dT_i/dt + \sum_j h_j(T_i)\,dm_{j,i}/dt$, which turns (7.9) into (7.10).
In the inventory form the inventory less its integrated rate is a linear invariant of the
semi-discrete system, conserved to round-off by every Runge–Kutta step; in the temperature
form energy is conserved to the integrator's order (§15.3.2, §16.5.2). The code is exactly
(`_energy_row`, `src/Assembly/residual.jl`):

```julia
@inline function _energy_row(::TemperatureState, E, S_T, c::Contributions{NC}, ds, i) where {NC}
    for j = 1:NC
        E -= ds.h[j, i] * c.M_flux[j, i]
    end
    return E / ds.cp_sum[i] + S_T
end

@inline function _energy_row(::EnthalpyState, E, S_T, c::Contributions{NC}, ds, i) where {NC}
    for j = 1:NC
        E += c.S_mass[j, i] * ds.h[j, i]
    end
    return E + ds.cp_sum[i] * S_T
end
```

with `E` the cell's energy source plus its scattered face energy fluxes, `c.M_flux` the
face-flux part of the mass rates, `c.S_mass` the reaction mass sources and `ds.h` the
sensible enthalpies at the cell temperature. `S_T` is an external temperature-rate source
[K·s⁻¹] (`VolumetricSources`, for verification by manufactured solutions; zero otherwise),
added as $S_T$ to the temperature row and as $C_iS_T$ to the inventory row. The storage
coefficient enters at exactly one place — the divisor $C_i$ of the temperature form, the
mixture value (7.7).

### 7.3.4 Species equations (companion to the energy equation)

The same row assembly writes the mass rows, the cell masses being the unknowns:

$$
\frac{dm_{j,i}}{dt}
\;=\; \omega_{j,i}\,V_i \;+\; \dot M_{j,i},
\qquad
\omega_{j,i} = \sum_r \nu_{r,j}\,r_{r,i},
\tag{7.12}
$$

with the reaction mass source written by the reaction kernel (§6.1.2) and the face fluxes
scattered by (7.4) (and, for verification, an external mass source). The gas-transport
kernel writes face fluxes for the gas components only (§9), so on a mesh that rides the
skeleton the solids and liquids move only through reactions; where a correction moves the mesh
relative to the solid, the relative-flux kernel carries every component across the faces
(§11.6). Because the masses are the states and each face flux enters two cells with opposite
signs, mass is conserved exactly.

---

## 7.4 The Gas Enthalpy Flux

### 7.4.1 Origin: the gas-enthalpy divergence

In the conservative balance (7.5) the sensible enthalpy carried by the gas appears as the
divergence of an enthalpy flux,

$$
\nabla\!\cdot\!\Bigl(\sum_g h_g(T)\,\mathbf{J}_g\Bigr),
\qquad
h_g(T) = \int_{T_{\text{ref}}}^{T} c_{p,g}(\tau)\,d\tau .
$$

Expanding this divergence with the product rule (exact even for temperature-dependent
$c_{p,g}$, since $\nabla h_g = c_{p,g}(T)\,\nabla T$),

$$
\nabla\!\cdot\!\Bigl(\sum_g h_g\,\mathbf{J}_g\Bigr)
=
\underbrace{\sum_g c_{p,g}\,\mathbf{J}_g\!\cdot\!\nabla T}_{\text{advection}}
\;+\;
\underbrace{\sum_g h_g\,\nabla\!\cdot\!\mathbf{J}_g}_{\text{generation/divergence}} .
\tag{7.13}
$$

The temperature form (7.6) keeps the **advection** piece, with a leading minus sign; the
generation piece cancels against $h_g$ times the gas species balance, with the reaction
heat quoted at the reaction temperature (§3.3.2, §6.6.3). Physically, the advection piece is
the rate at which flowing gas deposits or removes sensible heat as it crosses a temperature
gradient: gas moving from hot toward cold heats the cell it passes through, and vice versa. It
carries the enthalpy of the gas *flux* relative to the skeleton; the gas *storage* sits in
$\rho c_p^{\text{eff}}$ (§7.3.2).

Pyrolysis.jl does not evaluate the advection piece as a volumetric source. It assembles the
conservative flux itself: every face carries $\sum_j F_{j,f}\,h_j(T_f)$, the sensible enthalpy
of the **same** mass fluxes $F_{j,f}$ that the mass rows carry, and the temperature row (7.9)
subtracts the enthalpy the cell's net mass inflow brings at the cell temperature,
$\sum_j h_j(T_i)\dot M_{j,i}$. Together they are the face sum of (7.11), whose bracket
$h_j(T_f) - h_j(T_i)$ is the discrete advection (§7.4.2).

### 7.4.2 Discretization: the face temperatures

The face energy flux (7.8) gains $\sum_j F_{j,f}\,h_j(T_f)$ from every kernel that writes a
mass flux, each with the face temperature of its own flux:

- **Interior face, gas transport.** The gas-transport kernel (`GasTransportKernel`,
  `src/Kernels/gas_transport.jl`) evaluates the sensible enthalpies of the gas components at
  the **distance-weighted face temperature**

  $$
  T_f \;=\; \frac{d_n\,T_o + d_o\,T_n}{d_o + d_n},
  \tag{7.14}
  $$

  the linear interpolation of the two cell temperatures to the face, at the cell
  temperatures themselves. The mass flux it multiplies is computed separately, at the trial
  temperatures, with the arithmetic mean $\bar T = (T_o + T_n)/2$ in its driving potential
  (§9).
- **Boundary face.** The thermal boundary kernel adds $\sum_j F_{j,f}\,h_j(T_s)$ with the
  outward mass fluxes the gas-transport kernel wrote on that face (the film flux of the mass
  condition and, under Darcy–Fick transport, the vent the pressure condition drives) at the
  surface temperature $T_s$ of the face's energy balance, once that balance has fixed it
  (`_carried_enthalpy`, `src/Kernels/thermal_boundary.jl`; §7.4.3). Gas leaves at the surface
  temperature.
- **Relative-flux face** (a moving mesh under a correction). The relative-flux kernel
  (`RelativeFluxKernel`, `src/Kernels/relative_flux.jl`) evaluates the enthalpies at the face
  temperature reconstructed from the same upwind stencil as the face's concentrations (the
  upwind cell's value plus its smooth van Albada slope times its distance to the face,
  §11.6.2), so an isothermal composition profile advected across the mesh stays isothermal.

The construction has three properties:

1. **Conservative.** Each face's enthalpy flux leaves one cell and enters the next, so the
   column's sensible enthalpy changes only through the boundary faces and the sources.
2. **Exact for temperature-dependent $c_p$.** The enthalpies are the exact integrals of the
   heat-capacity laws (§5.4.3); no heat capacity is evaluated at an intermediate temperature.
3. **A uniform temperature stays at rest.** On an isothermal field $T_f = T_i$ on every face,
   so the bracket $h_j(T_f) - h_j(T_i)$ of (7.11) vanishes identically, whatever the mass
   fluxes. This is tested on a stretched mesh with Darcy–Fick transport, a heat capacity
   linear in $T$ and non-zero gas fluxes: the temperature rates vanish to $10^{-12}$ relative
   in the temperature form, and the inventory rates equal $\sum_j h_j\,dm_j/dt$ to the same
   accuracy in the inventory form (`test/kernels/test_gas_transport.jl`, "isothermal transport
   leaves the temperature at rest").

**Consistency.** On a line, let cell $i$ of thickness $\Delta z$ carry a uniform gas flux $F$
along $+z$ through its bottom face $L$ (where it is the neighbour, $s = -1$) and its top face
$R$ (where it is the owner, $s = +1$). With $T_{R,f} - T_i \approx \tfrac12\Delta z\,\partial_z T$
and $T_{L,f} - T_i \approx -\tfrac12\Delta z\,\partial_z T$, the gas part of the face sum of
(7.11) is

$$
-A\,F\bigl[h(T_{R,f}) - h(T_i)\bigr] + A\,F\bigl[h(T_{L,f}) - h(T_i)\bigr]
\;\approx\; -A\,\Delta z\,F\,c_p\,\frac{\partial T}{\partial z}
\;=\; -V_i\,c_p\,F\,\frac{\partial T}{\partial z},
$$

the cell integral of $-c_{p,g}\mathbf J_g\cdot\nabla T$ in (7.6). On a uniform mesh the
interpolation (7.14) is the central one; on a stretched mesh its distance weighting keeps it
exact for a linear temperature field.

### 7.4.3 The surface carries no gas-enthalpy term

The surface boundary condition is a *plain* convective + radiative balance,
$q_{\text{BC}}(T_s) = k\,(T_s - T_{\text{cell}})/d$, with no gas-enthalpy term
(`solve_boundary_temperature!`, `src/BoundaryConditions/surface_balance.jl`; §12.2). The gas
energy leaving the domain is carried by the boundary face's energy flux instead: its outward
mass flux carries its sensible enthalpy at $T_s$ (§7.4.2), and the heating of the gas on its
way to the surface was carried by the interior face fluxes as it crossed the cells.

Summing (7.11) over the column shows the bookkeeping. With heats of reaction at the reaction
temperature, the gas leaves the reaction at the local temperature with its sensible enthalpy
inside $h_r(T)$; the interior face terms $h(T_f) - h(T_i)$ heat it from its formation
temperature to the temperature of the cells it crosses; and the boundary face carries it out
at $T_s$. Under local thermal equilibrium the gas crosses the surface at $T_s$ on both sides,
so its enthalpy flux is continuous there and the surface balance is conducted flux equals
applied flux (§1.4.3). A surface term $\sum_g \dot m_g\,\Delta h_g(T_s,T_\infty)$ added to
that balance would apply the outflow a second time — a spurious sink of order
$\dot m''\,c_{p,g}\,(T_s - T_\infty)$, about $7.5\,\text{kW}\,\text{m}^{-2}$ for
$10\,\text{g}\,\text{m}^{-2}\text{s}^{-1}$ of gas with $c_{p,g} = 1.5$ kJ·kg⁻¹·K⁻¹ heated
by 500 K — that under-predicts the surface temperature and the mass-loss rate. The convention
is ThermaKin Eq. 17; Gpyro's thermal-equilibrium form uses the same term. Parameter sets
calibrated against a solver that applies a surface gas-enthalpy term *in addition to* the
in-depth gas enthalpy transport predict lower $T_s$ and MLR there than here; re-validate
imported calibrations (§18).

**The mass-flux/convection coupling is the blowing model of each `ConvectiveBC`.** The
*momentum-side* effect of the gas efflux — boundary-layer thickening that suppresses the
convective exchange coefficient, and the pooling of a dense gas on the face — is independent
of the enthalpy bookkeeping above and is the `blowing` model a `ConvectiveBC` carries
(§12.3.3): the coefficient becomes $h\,R_b$ with the blowing ratio $R_b$ from the model and the outward gas
flux at the face. Because it modifies only the coefficient, it composes exactly with the
face enthalpy flux; no double-count arises. The surface balance reads the gas-transport
fluxes of the face alone: the relative-flux kernel, which runs after the thermal boundary
kernel (§7.6), adds its fluxes only afterwards. This is the only use of the gas flux in the
surface balance.

> **Implementation note.** The thermal boundary kernel calls the surface Newton solve
> `solve_boundary_temperature!` (`src/BoundaryConditions/surface_balance.jl`) with the outward
> mass fluxes of the face (`J_out`), which feed only the blowing models, and the emissivity
> and absorptivity of the boundary cell; it then adds the conducted flux and the carried
> enthalpy at the solved $T_s$ to the face's energy flux (`src/Kernels/thermal_boundary.jl`).

---

## 7.5 Reaction-Heat and Radiation Sources

### 7.5.1 Reaction heat $Q_{\text{rxn}}$ and its sign

The volumetric reaction-heat source is

$$
\boxed{\;
Q_{\text{rxn},i} \;=\; -\sum_r h_r\,r_r
\;}
\qquad [\text{W}\,\text{m}^{-3}]
\tag{7.15}
$$

where $r_r$ is the rate of reaction $r$ (units $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$; §6)
and $h_r$ is its heat of reaction per unit mass of the first reactant (units
$\text{J}\,\text{kg}^{-1}$), quoted at the local temperature. The **internal storage sign
convention** (§2, overload note H1) is

- $h_r > 0$ **endothermic** (the reaction *absorbs* heat). Then $Q_{\text{rxn}}<0$: the cell
  **cools**.
- $h_r < 0$ **exothermic** (the reaction *releases* heat). Then $Q_{\text{rxn}}>0$: the cell
  **heats**.

The leading minus sign in (7.15) converts the enthalpy-of-reaction convention into the
heat-source convention. This sign is the most error-prone in the whole assembly, so it is
worth stating precisely against the implementation: the accumulator computes
`Q -= heats[i] * rates[i]` — i.e. $Q=-\sum_r h_r r_r$ — with rates $r_r\ge 0$ for forward
reactions, so an exothermic reaction ($h_r<0$) produces $Q>0$.

> **Caution on external conventions.** ThermaKin publishes $h>0=$ *exothermic*, the opposite
> of the internal storage convention here; Gpyro and FDS use the same sign as Pyrolysis.jl
> ($h>0$ endothermic, §6.9). With ThermaKin's published $h$ one would write
> $Q_{\text{rxn}}=+\sum h r$; Pyrolysis.jl stores $h$ with $h>0=$ endothermic and forms
> $Q_{\text{rxn}}=-\sum h r$, so a user porting a ThermaKin heat of reaction must flip its
> sign (§6.9.1).

> **Implementation note.** `reaction_heat_source(rates, heats)` in `src/Physics/kinetics.jl`
> performs the negated sum. The reaction kernel (`ReactionKernel`, `src/Kernels/reactions.jl`)
> evaluates the rates and heats of every cell with `evaluate_reactions!` (§6) at the cell's
> signed concentrations and trial temperature, with the material's depletion limiter, and adds
> $Q_{\text{rxn},i}V_i$ to the cell's energy source and $\omega_{j,i}V_i$ to its mass sources
> (§6.1.2).

### 7.5.2 Radiation source $Q_{\text{rad}}$ (dispatch only)

The absorbed radiation depends on the radiation model of the solve (`radiation_model`):

- **`NO_RADIATION`** and **`SURFACE_ABSORPTION`**: no volumetric source. The absorbed flux
  $a\,q_{\text{inc}}$ of a `RadiativeFluxBC` is applied through the surface energy balance of
  its face (§12.3.5), not as a volumetric deposit (opaque surface: the absorbed power raises
  $T_s$ directly). The two models give identical residuals.
- **`BEER_LAMBERT`**: in-depth absorption, $q'''(z)=\alpha(z)I(z)$ with
  $I(z)=I_{\text{inc}}\exp(-\!\int\!\alpha\,dz')$. The face intensities are algebraic unknowns
  of the state, with the rows $0 = I_{\text{down}} - I_{\text{up}}\,e^{-\alpha_i\Delta_i}$ for
  every live cell along the beam and $0 = I_e - a\,q_{\text{inc}}$ at the face the beam enters
  (the upstream face of the first live cell). The beam kernel (`BeerLambertKernel`,
  `src/Kernels/radiation.jl`) adds to each live cell's energy source
  $I_{\text{up}}A_{\text{up}} - I_{\text{down}}A_{\text{down}}$ [W], the exact cell average of
  the attenuated beam, and the radiative member of the boundary condition contributes zero to
  the surface balance.

The full derivation of each model — the Beer–Lambert cell-averaged form, the intensity rows,
the energy-equation coupling and the regime diagnostics — is the subject of **§8 (Thermal
Radiation)**. For the purposes of the energy assembly here, the absorbed radiation is part of
the cell's energy source $S_{E,i}$ in (7.8).

---

## 7.6 The Residual Sequence

The complete right-hand side is assembled by `residual!(du, u, def, t, ws)`
(`src/Assembly/residual.jl`) in a **fixed order**. The ordering is not arbitrary: several
steps consume results produced by earlier ones, and getting the order wrong would break the
surface energy balance, the enthalpy the gas carries out, or the mesh closure. One function
serves every configuration: the physics options are type parameters of the definition (the
kernel tuple, the gas transport, the closure), so the steps that do not apply are absent
rather than skipped (§15.2). The sequence is:

1. **The state mapping** (`derive!`, §5.13, §13): the only reader of the state vector. It
   produces the masses, volumes, concentrations and positive parts, the temperatures (the
   stored temperature, or the inversion of the enthalpy inventory) and their trial clamp, the
   heat-capacity sum and the sensible enthalpies, the cell properties, and the pore-gas
   properties the gas transport needs; on a moving mesh, the geometry from the node positions.
2. **The node corrections** (`correction!`): zero without a correction, the node
   velocities themselves under a prescribed mesh velocity, the solid being at rest
   (§11.7.8), and the relaxed
   equidistribution of each layer under the equidistribution correction (§11.7).
3. **Reset** (`reset!`): every face flux and cell source of the contributions is zeroed.
4. **The kernels**, in assembler order (`_kernels`, `src/Integration/setup.jl`), each adding
   its face fluxes and cell sources:
   1. the **reactions** (if the material has any): the mass sources $\omega_{j,i}V_i$ and the
      reaction heat $Q_{\text{rxn},i}V_i$ (§7.5.1, §6);
   2. the **gas transport** (if gas can move): the interior gas mass fluxes and the enthalpy
      they carry at the face temperature (7.14), and the outward mass fluxes of the boundary
      faces (§9, §7.4.2);
   3. the **in-depth radiation** (`BEER_LAMBERT`): the intensity residuals and the absorbed
      power (§7.5.2);
   4. **conduction**: the interior conductive fluxes (7.1) (and, in more dimensions, the
      cells' least-squares gradients, §7.1.3);
   5. the **thermal boundary**: the surface energy balance of every live boundary face, its
      conducted flux and the enthalpy its outward mass fluxes carry at $T_s$ (§7.4.3, §12);
   6. under a correction, the **relative fluxes**: every component's mass flux at the
      relative velocity, the enthalpy it carries, and the conversion-history pair (§11.6).
5. **The mesh closure** (`closure_rows!`, on a moving mesh): the cell strain rates
   $\theta_i$ from the reaction mass sources (`volume_change_rate` at the trial temperature),
   the history sources under a variable cross section, and the residuals of the velocity rows
   (§10, §11.3).
6. **The external sources** (`add_sources!`): a manufactured mass source, for verification.
7. **The row assembly** (`assemble_rows!`): the scatter (7.4) of every face flux, the mass
   rows (7.12), the energy row (7.9) or (7.10), the intensity, motion, history and
   cross-section rows, and last the dead rows of deactivated cells, which become identity rows
   (zero derivative on a differential row, the unknown itself on an algebraic one, §14.3.1).
   It is the only writer of `du`.

### 7.6.1 The ordering dependency graph

The dependencies that *fix* the order are:

1. **State mapping → everything.** Every kernel reads only the derived state (and the
   contributions of the kernels before it), never the state vector, so the mapping runs first
   and once.
2. **Corrections → kernels.** The relative-flux kernel reads the node corrections.
3. **Gas transport → thermal boundary.** The surface balance of a boundary face reads the
   face's outward mass fluxes, which feed the blowing models of its convective conditions, and
   the face's energy flux carries the sensible enthalpy of those fluxes at the solved $T_s$
   (§7.4.3).
4. **Thermal boundary → relative fluxes.** The blowing models must see the gas-transport
   fluxes of the face alone, not the relative fluxes of the mesh motion (§7.4.3).
5. **Conduction → thermal boundary** (more dimensions). The boundary kernel reads the
   least-squares cell gradients the conduction kernel computes, for the cross part of the
   boundary flux (§7.1.3).
6. **Reactions → closure.** The strain rates are the reactions' alone, read from the mass
   sources the reaction kernel wrote; a manufactured mass source, added after the closure,
   does not move the mesh.
7. **Assembly last.** It is the only writer of the derivative, so every kernel only adds to
   the contributions.

The residual as a function of the whole state — the algebraic rows, the mass matrix and the
integration — is described in §15 (Temporal Integration, the DAE, and the Colored Jacobian),
and the mesh closure and corrections in §11 (The ALE Framework).

---

## 7.7 Buffers and Element Types

### 7.7.1 The contributions

The kernels write into a pre-allocated buffer of contributions (`Contributions{NC,T}`,
`src/Kernels/contributions.jl`), which the row assembly reads. Its fields:

| field | shape | meaning | units |
|---|---|---|---|
| `F_mass` | $N_C\times n_{\text{faces}}$ | mass flux of every component per unit face area along the face normal | $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$ |
| `F_energy` | $n_{\text{faces}}$ | energy flux per unit face area along the normal: conduction and the enthalpy carried by the mass fluxes | $\text{W}\,\text{m}^{-2}$ |
| `S_mass` | $N_C\times n_{\text{cells}}$ | reaction mass source of every component, $\omega_{j,i}V_i$ | $\text{kg}\,\text{s}^{-1}$ |
| `S_energy` | $n_{\text{cells}}$ | cell energy source: reaction heat and absorbed radiation, without the formation term | W |
| `θ` | $n_{\text{cells}}$ | cell strain rate (mesh closure) | $\text{s}^{-1}$ |
| `M_flux` | $N_C\times n_{\text{cells}}$ | face-flux part of every mass rate, $\dot M_{j,i}$ (7.4), filled by the assembler | $\text{kg}\,\text{s}^{-1}$ |
| `T_face` | $n_{\text{faces}}$ | face temperatures; on a boundary face the surface temperature of its energy balance | K |
| `r_intensity` | $n_{\text{faces}}$ | residuals of the face-intensity rows (`BEER_LAMBERT`) | $\text{W}\,\text{m}^{-2}$ |
| `r_velocity` | $D\cdot n_{\text{nodes}}$ | residuals of the node-velocity rows (moving mesh) | those of the rows (§11.3) |
| `r_multiplier` | $n_{\text{cells}}$ | residuals of the closure multipliers (potential closure, more dimensions) | those of the rows (§11.3.5) |
| `F_carried`, `S_carried` | $2\times n_{\text{faces}}$, $2\times n_{\text{cells}}$ | flux and source of the conversion-history pair (variable cross section) | $\text{kg}\,\text{m}^{-2}\text{s}^{-1}$, $\text{kg}\,\text{s}^{-1}$ |
| `gradient` | $D\times n_{\text{cells}}$ | least-squares cell temperature gradients (more dimensions) | $\text{K}\,\text{m}^{-1}$ |
| `rates`, `heats`, `species` | $N_R$, $N_R$, $N_C$ | per-cell scratch of the reaction kernel | $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$, $\text{J}\,\text{kg}^{-1}$, $\text{kg}\,\text{m}^{-3}\text{s}^{-1}$ |

The cell properties the kernels read are the fields of the derived state (§5.13.2).

### 7.7.2 The workspace and the element type

The derived state and the contributions are held in the workspace (`CoreWorkspace`,
`src/Assembly/definition.jl`) with `Float64` elements, together with a twin of both for every
dual-number element type at which automatic differentiation evaluates the residual, built on
first use and selected by the element type of the state (`buffers`). In ordinary
integration the residual runs in the `Float64` buffers; the colored Jacobian (§15) evaluates
it on dual numbers in their twins, and the dense forward-mode Jacobian that serves as the test
reference does the same. Because every flux, source, property and geometric quantity of an
evaluation is held at the state's element type, the entire assembly propagates dual
sensitivities without a single conversion.

The buffers keep a fixed size for the whole solve: a deactivated cell keeps its slot and is
masked (its rows become identity rows, and the kernels skip every face beyond the live cells),
so no buffer is resized (§14.3.1). After warm-up a residual evaluation allocates nothing
(`test/assembly/test_residual_jacobian.jl` and `test/kernels/test_gas_transport.jl` assert
`@allocated residual!(du, u0, def, 0.0, ws) == 0`). The compile-time unrolling over the
component tuple — the `@generated` mixing sums and reaction evaluation, the recursive
heat-capacity sum — is what makes this type stability achievable over a heterogeneous
component tuple: the phase branches are resolved at code generation, so there is no runtime
tuple indexing (which would allocate) and the compiler specializes each arithmetic line to the
concrete element type.

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
temperature-form reading $-\sum_g c_{p,g}\mathbf J_g\cdot\nabla T$ of the face enthalpy flux
of §7.4 — and a deformation term that is the skeleton-following part of the time derivative,
which the moving mesh absorbs (§11.2.2). Pyrolysis.jl's energy equation is therefore ThermaKin
Eq. 17 term for term. ThermaKin maps to our symbols as: ThermaKin's stoichiometric
$\theta_i^{\,j}\to \nu_{r,j}$ (the strain rate $\theta$ is reserved here for dilation; §2 note
A2); ThermaKin's bimolecular concentration product $\xi_{\text{COMP1}}\xi_{\text{COMP2}}$ has
no analog because Pyrolysis.jl enforces a single reactant (§6); ThermaKin's gas flux
$J=-\rho_g\lambda\,\partial(\xi_g/\rho_g)/\partial x$ is the volume-fraction-gradient form
implemented here (§9).

### 7.8.3 Gpyro

Gpyro's energy equation (Lautenberger) uses a control-volume conservative formulation on
$\bar\rho\,\Delta z$, with the bulk (averaged) density $\bar\rho$ (our $\rho$),
volume-fraction ($X_\alpha\to v_j$) property averaging, and power-law temperature-dependent
properties. Its thermal-equilibrium form carries the gas through the solid-to-gas heat
transfer $\dot Q_{s-g}''' \approx \dot m'' c_{pg}\,\partial T/\partial z$ (its Eq. 65) — the
temperature-form reading of our face enthalpy flux with the Darcy mass flux — and drops the gas
storage $\rho_g\bar\psi c_{pg}\,\partial T/\partial t$ by default; its `FULL_QSG` option
restores it, which makes the equation (7.6) with (7.7). The mapping: Gpyro's pre-exponential
$Z\to A_r$; its conversion variable $\alpha_{\text{conv}}$ has no single analog (Pyrolysis.jl
tracks the component masses directly, with $\chi$ as a separate progress variable for the
cross section, §10); its permeability $K\to\kappa$; its mass flux
$\dot m''=-(K/\nu)(\partial P/\partial z - g\rho_g)$ is the Darcy term of §9 (Pyrolysis.jl
omits the gravity/buoyancy term $g\rho_g$ as negligible in the thin condensed phase).
Numerically Gpyro uses Patankar control-volume discretization with Peclet-based scheme
selection, harmonic interface interpolation, and a TDMA solve; Pyrolysis.jl uses the same
harmonic interface interpolation (7.2) but a method-of-lines discretization advanced by a
stiff integrator with the colored Jacobian of the whole residual (§15), rather than a
hand-rolled implicit TDMA sweep.

### 7.8.4 FDS

The FDS solid-phase model solves a 1D Fourier conduction equation with the surface balance
$-k\,\partial T/\partial x = \dot q_c'' + \dot q_r''$, Arrhenius reaction kinetics
($\rho^n O_2^n T^n$ rate form), and no in-depth gas: volatiles leave instantaneously, so
there is no gas storage and no in-depth gas-enthalpy coupling; the gas enthalpy is closed
at the surface. FDS maps as: pre-exponential $Z\to A_r$; solid absorption coefficient
$\kappa_s\to\alpha$; the surface energy balance with $\dot q_c''$ (convective) and
$\dot q_r''$ (radiative) is the boundary Newton solve of §12. The radiation closure differs
(FDS uses a two-flux Schuster–Schwarzschild model in depth; Pyrolysis.jl offers
`NO_RADIATION`, `SURFACE_ABSORPTION` and `BEER_LAMBERT`, §8). Pyrolysis.jl, following
ThermaKin, resolves the gas in depth and carries its enthalpy across every face with its mass
flux rather than as a surface outflow term, and the surface balance carries no gas-enthalpy
term on top of it (§7.4.3).

### 7.8.5 Summary table

| feature | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS |
|---|---|---|---|---|
| storage $\rho c_p$ | all components (7.7) | all components (incl. gas) | matrix by default; all with `FULL_QSG` | matrix (no in-depth gas) |
| gas enthalpy coupling | $\sum_j F_{j,f}h_j(T_f)$ at every face (temperature form $-\sum_g c_{p,g}\mathbf J_g\cdot\nabla T$) | $c_g\,\mathbf J_g\cdot\partial T/\partial x$ volumetric (Eq. 17) | $\dot m'' c_{pg}\,\partial T/\partial z$ (thermal equilibrium) | surface outflow |
| face conductivity | distance-weighted harmonic mean | harmonic | harmonic (Patankar) | — |
| conduction discretization | cell-centered FV, method of lines + stiff integrator | Crank–Nicolson (x) | Patankar CV + TDMA | implicit FD |
| radiation | none / surface / Beer–Lambert | Beer–Lambert + re-radiation | in-depth exp absorption | two-flux S–S |
| heat-of-reaction sign | $h>0$ endothermic, $Q_{\text{rxn}}=-\sum hr$ | $h>0$ exothermic (flip) | $h>0$ endothermic | $h>0$ endothermic |

---

## 7.9 Limitations and Open Issues

The energy assembly inherits a small number of known limitations, drawn from the source notes
and discussed where relevant above:

1. **One temperature.** The pore gas is stored at the matrix temperature (§7.3.2); there is
   no two-temperature gas energy equation, so regimes where the gas does not equilibrate
   with the pore walls within a cell are outside the model (§3.9).

2. **One global radiation model.** The radiation model is chosen once for the whole domain;
   per-cell or temperature-dependent model switching is not supported (§8).

3. **Sensible-enthalpy datum, no pressure work.** The conserved inventory is the sensible
   enthalpy $\sum_j m_{j,i}h_j(T_i)$ on the $T_{\text{ref}}$ datum, with the reaction heat
   entering as a source at the local temperature (§7.3.1, §16.4); no formation enthalpies are
   carried, so a heat of reaction quoted at the $T_{\text{ref}}$ datum must be converted before
   use (§6.6.3). No pressure–volume work enters the energy row: the pore-gas pressure drives
   the Darcy flux (§9) and does no work on the matrix in the model.

These do not affect the correctness of the discrete update (7.9)–(7.12): the energy ledger of
§16 closes the discrete balance of either energy form, whatever the composition — to
round-off in the inventory form and at the integrator's order in the temperature form. They
bound the *modeling* fidelity of the energy equation and are exposed, not hidden, through the
conservation ledgers and the audit of §16.
