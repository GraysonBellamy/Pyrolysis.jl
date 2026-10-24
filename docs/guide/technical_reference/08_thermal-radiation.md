# 8. Thermal Radiation

This chapter derives and documents the thermal-radiation models implemented in
Pyrolysis.jl and the way each couples into the energy equation. The solver offers
three mutually exclusive treatments of the incident radiation — `NO_RADIATION`,
`SURFACE_ABSORPTION` and `BEER_LAMBERT` — each with its own state, Jacobian
structure and physical regime of validity. We give the continuous radiative
transfer equation (RTE) and its reduction to the Beer–Lambert law, the exact
cell-averaged discretization, the face intensities that carry the beam as
algebraic unknowns of the state, and the banded Jacobian those unknowns give.
Every governing form, discrete update, constant, sign, factor and guard below
has been checked against the implementation; where a docstring or textbook form
diverges from the code, **the code is authoritative** and the divergence is
called out explicitly.

Throughout, the spatial coordinate is $z$ (axial, through-thickness), with $z=0$
the back face (tag `:bottom`) and $z=L$ the exposed surface (tag `:top`). Heat
enters from the top. A face flux is stored per unit face area along the face
normal, which points from the face's owner cell to its neighbour and out of the
domain on a boundary face: on a line the interior faces and the `:top` face have
their normal along $+z$, the `:bottom` face along $-z$ (§2). The incident beam
travels in $-z$; its intensity $I_f$ at face $f$ is stored as a non-negative
magnitude, and the radiative source of a cell is **positive = heating**. The cell
index is $i$ ($1\le i\le n$), the component index $j$ ($1\le j\le N_C$). The
Stefan–Boltzmann constant is $\sigma = 5.670374419\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}$
(`STEFAN_BOLTZMANN`, `src/Materials/constants.jl`).

---

## 8.1 The Radiation Model Family

### 8.1.1 The three models

The radiation treatment is selected at `solve` (or `setup`) time by one
enumeration value passed as `radiation_model`; the default is `NO_RADIATION`. The
enum is defined in `src/Materials/radiation_types.jl`:

```julia
@enum RadiationModel begin
    NO_RADIATION = 0
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
end
```

The models are **mutually exclusive** — exactly one is active for the whole domain
in a given run — and they differ in where the incident external radiation of a
`RadiativeFluxBC` is deposited:

1. **`NO_RADIATION`** and **`SURFACE_ABSORPTION`.** The incident radiation is
   absorbed at the exposed face. It is **not** a volumetric source: the absorbed
   flux $a\,q_{\text{inc}}(t)$ enters the energy balance as a boundary heat flux
   through the surface energy balance and its Newton solve for $T_s$ (§8.2, §12.2).
   Component absorption coefficients $\alpha_j$ are not read. The two names select
   **the same residual**: the `RadiativeFluxBC` dispatch treats them identically
   (`evaluate_bc_flux_and_derivative`, `src/BoundaryConditions/surface_balance.jl`),
   and the kernel tests assert identical residuals for the two
   (`test/kernels/test_boundary_conditions.jl`). `NO_RADIATION` names the intent
   of a run in which radiation does not matter (no `RadiativeFluxBC`, or a
   transparent sample, $\tau \ll 0.1$); `SURFACE_ABSORPTION` names the intent of
   an opaque solid whose first cell is optically thick, $\tau_n\gg 1$. Surface
   re-radiation (`RadiativeBC`) acts at the face under every model.

2. **`BEER_LAMBERT`.** In-depth exponential absorption following the Beer–Lambert
   law, $I(z)=I_e\exp\!\big(-\int_z^L\alpha_{\text{eff}}\,\mathrm{d}z'\big)$, with
   volumetric source $q'''=\alpha_{\text{eff}}I$. Radiation penetrates over the
   absorption length $1/\alpha_{\text{eff}}$. The intensity of every face is an
   **algebraic unknown** of the state (§8.3.3); each of its rows reads one cell and
   that cell's two faces, so the Jacobian stays banded (§8.6). Appropriate for
   semi-transparent polymers and foams ($0.1\lesssim\tau\lesssim 10$).

There is no spherical-harmonics (diffusion) or two-flux model of the radiation
in depth; §8.8 states what one would take.

### 8.1.2 What the model changes

The model-dependent decisions are made where they arise, all of them by `setup`
(`src/Integration/setup.jl`):

- **The state.** Under `BEER_LAMBERT` the layout carries the intensity of every
  face as an algebraic unknown (`StateLayout(...; intensities = true)`): on a line
  of $n$ cells, $n+1$ rows, one in the node block of each face's node
  (`intensity_index(layout, f)`, §15.1). The other two models add nothing.
- **The kernels.** `BEER_LAMBERT` adds a `BeerLambertKernel`
  (`src/Kernels/radiation.jl`) after the gas-transport kernel and before
  conduction; the other two add none. The thermal boundary kernel reads the model
  for the `RadiativeFluxBC` dispatch (§8.4).
- **The tolerances and the domain test.** The intensity rows take the absolute
  tolerance $\mathrm{reltol}\cdot S_I$, with $S_I$ the largest incident intensity
  over the time span (§8.3.3); the domain test checks their residuals (§15.5.4).
- **The refusals.** `setup` refuses `BEER_LAMBERT` with an `ArgumentError` off
  planar coordinates, on a mesh of more than one dimension, when the `:top` tag
  has more than one face, and when a `RadiativeFluxBC` sits on any tag other than
  `:top` (§8.3.4).

No check compares the material with the model: `SURFACE_ABSORPTION` with non-zero
$\alpha_j$ (not read) and `BEER_LAMBERT` with every $\alpha_j = 0$ (nothing
absorbed in depth, everything transmitted, §8.3.3) both run without a warning.

| Model | Extra state | Kernel | Jacobian pattern |
|---|---|---|---|
| `NO_RADIATION` | — | — | banded (cell blocks and their neighbours) |
| `SURFACE_ABSORPTION` | — | — | the same |
| `BEER_LAMBERT` | $n+1$ face intensities (algebraic) | `BeerLambertKernel` | banded; each intensity row reads three intensities and the two cells beside its face (§8.6) |

---

## 8.2 Surface Absorption

### 8.2.1 Surface absorption is a boundary flux, not a volumetric source

Under `NO_RADIATION` and `SURFACE_ABSORPTION` no kernel deposits radiation in
the cells. The incident radiation is folded into the surface energy balance of
the face carrying the `RadiativeFluxBC` (`evaluate_bc_flux_and_derivative`,
`src/BoundaryConditions/surface_balance.jl`):

$$
q_{\text{BC}}(T_s) = a\,q_{\text{inc}}(t),\qquad
\frac{\mathrm{d}q_{\text{BC}}}{\mathrm{d}T_s}=0,
\qquad\text{for }\texttt{NO\_RADIATION}\text{ and }\texttt{SURFACE\_ABSORPTION},
$$

where $q_{\text{inc}}(t)$ [W·m⁻²] is the condition's incident flux and $a$ the
surface absorptivity. With no explicit `absorptivity` on the condition, $a$ is the
effective surface absorptivity of the boundary cell (§5.8.2, §12.2.6),

$$
a_{\text{eff}}=\sum_j v_j\,a_j ,
$$

with $v_j$ the mixing volume fractions and $a_j$ the component's
`absorptivity` when set and otherwise its emissivity $\varepsilon_j(T)$; with no
component absorptivity set this is Kirchhoff's law, $a_{\text{eff}}=\varepsilon_{\text{eff}}$.
The state mapping evaluates it once per cell and residual,
`ds.a[i] = mixture_absorptivity(ξ⁺, T_trial, material)` (`src/States/derived.jl`),
at the positive-part concentrations and the cell's trial temperature. No clamp is
applied: the weights are convex and every component value lies in $[0,1]$
(validated on the property bracket when the component is built), so the sum stays
in $[0,1]$. An explicit `absorptivity` on the condition (a constant or a function
of $t$) overrides $a_{\text{eff}}$. The absorbed incident flux does not depend on
$T_s$. Surface re-radiation is supplied separately by a `RadiativeBC` term,
$q=F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$, whose $T_s$-dependence
enters the Newton derivative (§12.3.4).

The textbook alternative — depositing the absorbed flux as a volumetric source
$a\,q_{\text{inc}}/\Delta z_n$ in the top cell — is deliberately **not**
implemented. The boundary-flux route lets the surface Newton solve resolve the
(possibly large) difference between $T_s$ and the cell temperature $T_n$
consistently with convection and re-radiation; for an opaque material with a
well-resolved first cell the two constructions are numerically close anyway.

### 8.2.2 Assumptions and validity

Surface absorption presumes the material is optically thick over the first cell,
$\tau_n=\alpha_{\text{eff}}\Delta z_n\gg 1$, so that essentially no radiation
penetrates beyond cell $n$. It adds no unknowns and no Jacobian entries beyond the
boundary cell's own, and it is the recommended choice for charring solids and fast
cone-calorimeter analysis.

---

## 8.3 The Beer–Lambert In-Depth Model

### 8.3.1 Continuous form

Beer–Lambert follows from the 1D radiative transfer equation for a collimated beam
with absorption coefficient $\alpha_{\text{eff}}(z)$ and negligible scattering and
emission along the *incident* beam. With the beam entering at the exposed surface
$z=L$ and travelling into the depth, the depth $s=L-z$ is the path length, and

$$
\frac{\mathrm{d}I}{\mathrm{d}s} = -\alpha_{\text{eff}}\, I
\qquad\Longleftrightarrow\qquad
\frac{\mathrm{d}I}{\mathrm{d}z} = +\alpha_{\text{eff}}\, I .
$$

Integration from the surface gives exponential attenuation,

$$
I(z) = I_e\,\exp\!\left(-\int_{z}^{L}\alpha_{\text{eff}}(z')\,\mathrm{d}z'\right),
$$

with $I_e$ the intensity that enters through the exposed face (§8.3.2). The local
energy deposited per unit volume is the absorbed power density,

$$
q'''(z) = \alpha_{\text{eff}}(z)\,I(z)\qquad[\mathrm{W\,m^{-3}}].
$$

The effective absorption coefficient is the composition-weighted volumetric
extinction,

$$
\alpha_{\text{eff}} = \sum_j \xi^+_j\,\alpha_j,
$$

where $\xi^+_j$ [kg·m⁻³] is the positive-part mass concentration and $\alpha_j$
[m²·kg⁻¹] the component's constant mass absorption coefficient (keyword `α`;
`effective_absorption`, `src/Properties/effective_properties.jl`). The product
$\xi_j\alpha_j$ has units m⁻¹, so $\alpha_{\text{eff}}$ is a *volumetric*
coefficient. The state mapping evaluates it once per cell and residual, and only
under `BEER_LAMBERT` (`ds.α[i]`, `src/States/derived.jl`; §5.8, §5.13). The
optical thickness of the whole layer is $\tau = \int_0^L \alpha_{\text{eff}}\,\mathrm{d}z$.

### 8.3.2 The exact cell-averaged discretization

A naive evaluation $q'''_i \approx \alpha_i I(z_i^c)$ at the cell center is only
first-order accurate and can over- or under-deposit energy depending on $\tau_i$.
Pyrolysis.jl instead uses the **exact cell-integrated** source: integrating
$q'''=\alpha I$ over a cell of constant $\alpha_i$ and entrance intensity
$I_{\text{up}}$ (the intensity at the cell's face nearer the exposed surface),

$$
\int_{\text{cell }i} \alpha_i\,I\,\mathrm{d}s
= \int_0^{\Delta z_i}\!\alpha_i\, I_{\text{up}}\, e^{-\alpha_i s}\,\mathrm{d}s
= I_{\text{up}}\big(1 - e^{-\alpha_i \Delta z_i}\big),
$$

so the *cell-averaged* volumetric source is the absorbed power divided by the cell
thickness:

$$
\boxed{\;q^{\text{rad}}_i = I_{\text{up}}\,\frac{1 - e^{-\tau_i}}{\Delta z_i},
\qquad \tau_i \equiv \alpha_i\,\Delta z_i\;}
$$

This is exact for piecewise-constant $\alpha$ and conserves energy by
construction: the power absorbed in cell $i$ per unit area is
$q^{\text{rad}}_i\,\Delta z_i = I_{\text{up}}(1-e^{-\tau_i}) = I_{\text{up}}-I_{\text{down}}$,
exactly the drop in beam intensity across the cell, with
$I_{\text{down}} = I_{\text{up}}e^{-\tau_i}$ the intensity at its far face. The
implementation carries exactly these face intensities (§8.3.3).

The intensity that crosses the exposed face is $I_e = a\,q_{\text{inc}}(t)$, with
$a$ the surface absorptivity: on the default **Kirchhoff path** (the
`RadiativeFluxBC` carries no explicit `absorptivity`) it is the effective surface
absorptivity $a_{\text{eff}}$ of the first live cell from the exposed face
(`ds.a`, §8.2.1, §5.8.2) — identical to its $\varepsilon_{\text{eff}}$ unless a
component sets its own `absorptivity`, so by default absorption and surface
emission use the same coefficient (Kirchhoff's law, detailed balance) — and the
remaining $(1-a)\,q_{\text{inc}}$ is reflected and never enters the medium. With
an **explicit absorptivity** on the condition, $a(t)$ is that value at the time of
the evaluation. The entrance intensity of cell $i$ is then the cumulative
attenuation from the surface down to its upstream face,

$$
I_{\text{up}}(i) = a\,q_{\text{inc}}\prod_{k>i} e^{-\tau_k}
= a\,q_{\text{inc}}\,\exp\!\left(-\sum_{k>i}\alpha_k\Delta z_k\right),
$$

the product running over the live cells between cell $i$ and the exposed face.
The transmitted (un-absorbed) intensity at the back face is
$I_b = a\,q_{\text{inc}}\,e^{-\tau}$; it leaves through the back face (§8.3.3).

### 8.3.3 The face intensities as algebraic rows

Eliminating the recurrence $I_{\text{down}} = I_{\text{up}}e^{-\tau_i}$ would make
the source of every cell depend on the state of every cell above it. The
`BeerLambertKernel` (`src/Kernels/radiation.jl`) keeps the recurrence local
instead: the intensity $I_f$ of every face is an algebraic unknown, and each face
carries one row. With the cells ordered from the exposed face into the depth
(`order`), and `up[i]`, `down[i]` the faces of cell $i$ nearer and farther from
the exposed face,

$$
\begin{aligned}
0 &= I_e - a\,q_{\text{inc}}(t) &&\text{(the exposed face $e$)},\\
0 &= I_{\text{down}(i)} - I_{\text{up}(i)}\,\exp(-\alpha_i\Delta_i) &&\text{(every other face, the far face of cell $i$)},
\end{aligned}
$$

and cell $i$ absorbs the power

$$
S^{\text{rad}}_i = I_{\text{up}(i)}\,A_{\text{up}(i)} - I_{\text{down}(i)}\,A_{\text{down}(i)}\qquad[\mathrm{W}],
$$

added to its energy source (`S_energy`). At intensities that satisfy the rows,
$S^{\text{rad}}_i/V_i$ is exactly the cell average $q^{\text{rad}}_i$ of §8.3.2. The
kernel, verbatim:

```julia
function contribute!(c::Contributions, kernel::BeerLambertKernel, ds::DerivedState,
                     topo::MeshTopology{1}, material, t)
    g = ds.geometry
    s = _first_live(kernel, ds)
    s == 0 && return c
    e = kernel.up[kernel.order[s]]
    c.r_intensity[e] = ds.I[e] - incident_intensity(kernel, ds, t)
    @inbounds for k = s:length(kernel.order)
        i = kernel.order[k]
        up, down = kernel.up[i], kernel.down[i]
        I_up, I_down = ds.I[up], ds.I[down]
        c.r_intensity[down] = I_down - I_up * exp(-_optical_thickness(kernel, ds, topo, i))
        c.S_energy[i] += I_up * g.face_area[up] - I_down * g.face_area[down]
    end
    return c
end
```

The assembler writes `r_intensity[f]` into the intensity row of face $f$
(`assemble_rows!`, `src/Assembly/residual.jl`). Several details are load-bearing:

- **Cell thickness is geometric.** $\Delta_i = d_{\text{up}} + d_{\text{down}}$, the
  sum of the cell center's distances to its two faces in the current geometry
  (the node positions on a moving mesh), each smoothly floored at $10^{-8}$ of its
  reference value (`DISTANCE_FLOOR_FRACTION`, `_optical_thickness`). On a line
  $\Delta_i = \Delta z_i$, so the attenuation follows the cell as the mesh moves
  and shrinks (§10, §11). Under a variable cross section every face of a planar
  column has the same area $A(t)$, so $A_{\text{up}} = A_{\text{down}}$.
- **The beam's energy telescopes for any intensities.** Summed over the live
  cells, $\sum_i S^{\text{rad}}_i = I_e A_e - I_b A_b$ identically, whether or not
  the intensities satisfy their rows: the power deposited in the cells is the
  power that entered less the power that left, also on a trial state off the
  algebraic manifold.
- **Transparent cells.** A cell with $\alpha_i = 0$ has $e^{-\tau_i}=1$: its row
  passes the intensity unchanged and it absorbs exactly nothing. No threshold
  guard is needed.
- **Zero flux.** A zero incident flux, or no `RadiativeFluxBC` on `:top` (the
  kernel then uses `RadiativeFluxBC(0.0)`), gives zero intensities through the
  rows; nothing is short-circuited.
- **Transmitted flux.** Whatever intensity reaches the last face, $I_b$, leaves
  through the back face: no row deposits it, and the back-face condition
  describes only the surface it models. That is the right answer for a slab whose
  back face is open and a modelling error for one that rests on an opaque holder.
  The solution reports it per saved time (`sol.extras.radiation_transmitted`,
  §8.5); no warning is issued.
- **Planar coordinates only.** The rows attenuate intensities while the cell
  sources read face areas; on a cylinder or a sphere the two areas of a cell
  differ, so a transparent cell would absorb $I\,(A_{\text{up}}-A_{\text{down}})$.
  `setup` refuses `BEER_LAMBERT` off planar coordinates
  (`"in-depth radiation needs planar coordinates"`).
- **Deactivation.** After cells are deactivated from the exposed side (§14), the
  beam enters the first live cell of `order` through its upstream face
  (`exposed_face`), the exposed-face row moves there, and $a$ is that cell's
  absorptivity. The intensity rows of the faces upstream of it are masked with
  the dead cells (§14.3.1): they become identity rows, and their intensities are
  zero.
- **A consistent start.** `setup` writes intensities that satisfy every row into
  the initial state (`consistent_state!`, which calls `sweep_intensities!`: the
  incident intensity at the exposed face, attenuated cell by cell, zero upstream
  of the first live cell), and the deactivation event does the same after every
  transfer (§15.4). The integrator runs no initialization of its own.
- **Tolerances and the domain test.** Each intensity row takes the absolute
  tolerance $\mathrm{reltol}\cdot S_I$, where $S_I$ is the largest
  $|I_{\text{inc}}(t)|=|a\,q_{\text{inc}}(t)|$ over 101 equally spaced times of the
  time span, at least 1 W·m⁻² (`_intensity_scale`, `src/Integration/setup.jl`;
  `abstol_vector`, `src/States/tolerances.jl`). The domain test rejects a proposed
  step whose live intensity rows leave a residual above
  $\mathrm{rtol}\cdot S_I + \mathrm{abstol}$ (§15.5.4).

The source is **positive = heating**. The absorbed power over the live cells is
$\sum_i S^{\text{rad}}_i = A\,a\,q_{\text{inc}}\,(1-e^{-\tau}) \le A\,a\,q_{\text{inc}}$,
with equality only if the layer is fully opaque ($e^{-\tau}\to 0$) — consistent
with `SURFACE_ABSORPTION`, which applies $a\,q_{\text{inc}}$ at the face with the
same $a$ (the absorptivity of the same cell), so switching models does not change
the absorbed power in the opaque limit.

### 8.3.4 Where the incident intensity comes from

The incident flux fed to Beer–Lambert is the flux of the `RadiativeFluxBC` on the
`:top` tag: the condition itself, or the **first** `RadiativeFluxBC` member of a
`CombinedThermalBC` (`get_radiative_flux_bc`,
`src/BoundaryConditions/boundary_conditions.jl`). `incident_intensity(kernel, ds, t)`
returns

$$
I_{\text{inc}}(t) = a\,q_{\text{inc}}(t),\qquad
a = \begin{cases}
\text{the condition's absorptivity at } t, & \text{when it has one},\\
a_{\text{eff}} \text{ of the first live cell}, & \text{otherwise (Kirchhoff path)}.
\end{cases}
$$

The choice between the two paths is a property of the condition and is made by
the same code in every evaluation; there is no separate flag to keep in step
between the residual and its Jacobian, which differentiates the same residual
(§8.6).

Beer–Lambert models **top-incident** radiation only. Under this model a
`RadiativeFluxBC` contributes zero to the surface balance of *every* face (§8.4),
so a `RadiativeFluxBC` on any other tag could not be honored: `setup` refuses it
with an `ArgumentError` (`"in-depth radiation enters through the :top faces
only; the condition on :bottom carries an incident radiative flux"`) rather than
silently dropping the back-side heating. Use `SURFACE_ABSORPTION` (face deposit
at either boundary) or a `HeatFluxBC` for back-side radiant heating. The `:top`
tag must have exactly one face, and a mesh of more than one dimension is refused
(`"in-depth radiation is a line's"`, §8.8).

A second `RadiativeFluxBC` member of the `:top` condition is **dropped
silently**: only the first feeds the beam, and every member contributes zero to
the surface balance. With the slab of §8.3.5 and a top condition
`RadiativeFluxBC(flux = 50e3) + RadiativeFluxBC(flux = 20e3) + RadiativeBC(T_inf = 300.0)`,
the incident intensity is $0.9\times 50 = 45\ \mathrm{kW\,m^{-2}}$ under
`BEER_LAMBERT`, while `SURFACE_ABSORPTION` absorbs $0.9\times 70 = 63\ \mathrm{kW\,m^{-2}}$
at the face (run in the documentation scratch environment). Combine the sources
into one `RadiativeFluxBC` (a function of $t$ can sum them).

The incident flux is evaluated inside **every** residual evaluation, at that
evaluation's own time $t$, through `incident_intensity`. The colored Jacobian
differentiates the same residual at its own time, so it sees the flux at the time
it is evaluated. There is no constant-flux detection and no snapshot: a ramped or
modulated cone flux can never be frozen at its $t=0$ value. The kernel tests run
a sinusoidally modulated flux and check the exposed-face intensity against
$a\,q_{\text{inc}}(t)$ at the residual's time (`test/kernels/test_radiation.jl`).

### 8.3.5 Worked numbers

For a 10 mm semi-transparent slab with $\alpha_{\text{eff}}=100\ \mathrm{m^{-1}}$
($\rho = 1000\ \mathrm{kg\,m^{-3}}$, $\alpha_j = 0.1\ \mathrm{m^2\,kg^{-1}}$),
uniform $\Delta z=1\ \mathrm{mm}$, $q_{\text{inc}}=50\ \mathrm{kW\,m^{-2}}$ and
$\varepsilon_{\text{eff}}=0.9$ (Kirchhoff path), the intensity entering the
medium is $I_e = 0.9\times50 = 45\ \mathrm{kW\,m^{-2}}$. With $\tau_i=0.1$ the top
cell ($i=n=10$) absorbs $q_{10}=45\,000\,(1-e^{-0.1})/0.001 \approx 4.28\ \mathrm{MW\,m^{-3}}$;
the intensity at its far face is $45\,000\,e^{-0.1}=40.7\ \mathrm{kW\,m^{-2}}$,
giving $q_9\approx 3.88\ \mathrm{MW\,m^{-3}}$, then $q_8\approx 3.51\ \mathrm{MW\,m^{-3}}$,
and so on. The total optical thickness is $\tau=\alpha L=1$, so
$e^{-1}\approx 0.37$ of the entering intensity, $16.55\ \mathrm{kW\,m^{-2}}$, is
transmitted through the slab — Beer–Lambert is the right tool in this regime,
whereas surface absorption would deposit everything at the surface.

The run below reproduces these numbers:

```julia
using Pyrolysis

slab = SolidComponent(:slab, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.1)
material = Material(name = :slab, components = (slab,), reactions = ())
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, 10),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (top = RadiativeFluxBC(flux = 50e3), bottom = AdiabaticBC()),
    ),
    ξ_initial = [z -> 1000.0],
    tspan = (0.0, 1.0),
    output_times = [0.0, 1.0],
)
sol = solve(problem; radiation_model = BEER_LAMBERT)
sol.extras.volumetric_absorption        # 28445.4 W/m² absorbed in depth
sol.extras.radiation_transmitted        # 16554.6 W/m² through the back face
sol.extras.surface_heat_flux_absorbed   # 0.0: nothing is absorbed at the face
```

The initial face intensities of this problem, from the exposed face down, are
45000.0, 40717.7, 36842.9, … , 16554.6 W·m⁻², and the cell sources of the top three
cells 4.282, 3.875 and 3.506 MW·m⁻³. With `radiation_model = SURFACE_ABSORPTION` the
same problem absorbs 45000 W·m⁻² at the face, all of it conducted into the slab, and
nothing in depth.

---

## 8.4 Energy-Equation Coupling for Each Model

The condensed-phase energy equation (§3.3, §7.3) carries a volumetric radiation
source $Q_{\text{rad}}$ on its right-hand side:

$$
\rho c_p^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\!\mathbf{q}_{\text{cond}} - \sum_j \mathbf{J}_j\!\cdot\!\nabla h_j
+ Q_{\text{rad}} + Q_{\text{rxn}} .
$$

Discretely, every source and face flux of cell $i$ is collected in its energy
right-hand side $E_i$ [W], and the energy row is (`assemble_rows!`,
`src/Assembly/residual.jl`; §7.3)

$$
C_i\,\frac{\mathrm{d}T_i}{\mathrm{d}t} = E_i - \sum_j h_j(T_i)\,\dot M_{j,i}
\quad\text{(temperature form)},\qquad
\frac{\mathrm{d}H_i}{\mathrm{d}t} = E_i + \sum_j \omega_{j,i}V_i\,h_j(T_i)
\quad\text{(enthalpy form)},
$$

with $C_i$ the heat-capacity sum of the cell (§5.4) and $\dot M_{j,i}$ the
face-flux part of its mass rate. The three models populate $Q_{\text{rad}}$ — or
bypass it for the surface route — as follows.

- **`NO_RADIATION` and `SURFACE_ABSORPTION`.** $Q_{\text{rad}}=0$: no kernel adds a
  radiative source. A `RadiativeFluxBC` enters the surface energy balance of its
  face as $q_{\text{BC}}=a\,q_{\text{inc}}$ (§8.2); the balance solves for $T_s$,
  and the flux conducted from the face into the boundary cell, $k(T_s-T_n)/d$,
  enters that cell's $E_n$ through the boundary face (§12.2). A `RadiativeBC`
  acts in the same balance.

- **`BEER_LAMBERT`.** Cell $i$'s energy source gains the absorbed power
  $S^{\text{rad}}_i = I_{\text{up}}A_{\text{up}} - I_{\text{down}}A_{\text{down}}$
  [W] of §8.3.3, the same term in both energy forms; in the temperature form it is
  divided by $C_i$ with the rest of the row, so
  $(\mathrm{d}T_i/\mathrm{d}t)_{\text{rad}} = S^{\text{rad}}_i/C_i =
  q^{\text{rad}}_i/\rho c_{p,i}^{\text{eff}}$ at consistent intensities. The
  `RadiativeFluxBC` contributes **zero** to the surface balance under this model
  (`evaluate_bc_flux_and_derivative` returns $(0, 0)$), so the incident energy is
  not double-counted: it is already deposited in the cells. A `RadiativeBC` stays
  active at the face: Beer–Lambert deposits only the *incident* beam, and the
  emission of the material remains a surface process.

**The radiation-model dispatch of `RadiativeFluxBC` is the linchpin.** It keeps
the incident energy single-counted across the models: deposited at the face for
`NO_RADIATION` and `SURFACE_ABSORPTION`, in the cells (and therefore zero at the
face) for `BEER_LAMBERT`.

| Condition | `NO_RADIATION`, `SURFACE_ABSORPTION` | `BEER_LAMBERT` |
|---|---|---|
| `RadiativeFluxBC` at `:top` | $a\,q_{\text{inc}}$ in the surface balance | zero in the surface balance; the first member feeds the beam |
| `RadiativeFluxBC` at another tag | $a\,q_{\text{inc}}$ in that face's balance | refused by `setup` |
| `RadiativeBC` | $F\varepsilon\sigma(T_\infty^4-T_s^4)$ in the surface balance | the same |
| `ConvectiveBC`, `HeatFluxBC`, … | unchanged | unchanged |

The same dispatch threads through `CombinedThermalBC`, so a cone-calorimeter
condition such as
`RadiativeFluxBC(flux = 50e3) + RadiativeBC(T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)`
behaves correctly under every model: under `BEER_LAMBERT` its radiative-flux
member contributes exactly zero at the face and feeds the beam, while the
re-radiative and convective members act at the face. The diagnostics follow the
same split: `sol.extras.surface_heat_flux_absorbed` is the flux absorbed at the
face (zero under `BEER_LAMBERT`), and `sol.extras.volumetric_absorption` the flux
absorbed in depth (zero under the other two models).

---

## 8.5 Regime Diagnostics

The validity bands of the models are stated in terms of the optical thickness of
the layer, $\tau=\alpha_{\text{eff}}L=\int_0^L\alpha_{\text{eff}}\,\mathrm{d}z$
(`RadiationModel` docstring, `src/Materials/radiation_types.jl`):

| $\tau=\alpha_{\text{eff}} L$ | regime | recommendation |
|---|---|---|
| $\tau>10$ | optically thick | surface absorption is adequate |
| $0.1<\tau<10$ | semi-transparent | `BEER_LAMBERT` |
| $\tau<0.1$ | transparent | `NO_RADIATION` |

The per-cell optical thickness $\tau_i=\alpha_i\Delta z_i$ is a resolution
measure, not a validity criterion: the exact cell average of §8.3.2 is correct at
any $\tau_i$, and refining the mesh shrinks every $\tau_i$ without changing the
physics. What it decides is whether `SURFACE_ABSORPTION` is close to the in-depth
answer, which needs $\tau_n\gg1$ in the exposed cell (§8.2.2).

No function computes $\tau$ for a material; it follows from the component values,
$\tau = L\sum_j\xi_j\alpha_j$ for a uniform layer. A run diagnoses itself through
two series of `sol.extras`, both per unit exposed area at every saved time:
`volumetric_absorption`, $I_e - I_bA_b/A_e$, the intensity the cells absorbed, and
`radiation_transmitted`, $I_bA_b/A_e$, the intensity leaving through the back face
(`src/Outputs/solution.jl`). Their sum is the entering intensity $I_e$ by
construction, and by the telescoping of §8.3.3 the first is exactly the power the
cells absorbed per unit exposed area. Their ratio gives the optical thickness the
beam actually crossed,

$$
\tau = \ln\frac{I_e}{I_b} = \ln\frac{\text{volumetric\_absorption}+\text{radiation\_transmitted}}{\text{radiation\_transmitted}},
$$

which for the slab of §8.3.5 is $\ln(45\,000/16\,554.6)=1.000$. A transmitted
fraction near one signals a layer that absorbs little in depth (every $\alpha_j=0$
transmits the whole absorbed share of the incident flux); a fraction that stays
near zero with most of the absorption in the exposed cell signals that surface
absorption would do.

---

## 8.6 The Beer–Lambert Jacobian

### 8.6.1 Banded, not triangular

Eliminated, the attenuation recurrence makes the source of cell $i$ depend on the
absorption coefficient of every cell between it and the exposed face, through
$I_{\text{up}}(i)= a\,q_{\text{inc}}\exp(-\sum_{k>i}\alpha_k\Delta z_k)$: a dense
upper-triangular block in the Jacobian. With the face intensities as unknowns,
every row is local (§8.3.3), and the pattern built from the mesh connectivity
(`jacobian_pattern`, `src/ColoredJacobian/pattern.jl`; §15.11) is banded:

- the intensity row of face $f$ reaches the intensities of nodes $f-1$, $f$ and
  $f+1$, the masses and energy of the cells on either side of the face (the
  attenuation of the cell the beam leaves through it, and for the exposed face
  the absorptivity of the exposed cell), on a moving mesh the positions of nodes
  $f-1\dots f+1$, and under a variable cross section the area;
- every cell row reads the intensities of its two faces.

Every row holds its diagonal, and the coloring stays bounded independently of the
cell count: the kernel tests assert at most $3(N_C+1)+2$ colors on a fixed mesh
(`test/kernels/test_radiation.jl`). The Jacobian is the colored forward-mode
derivative of the residual itself (§15.6); no derivative of the radiation terms is
hand-coded.

### 8.6.2 The entries

For a live cell $i$ with upstream face $u$ and downstream face $d$, the colored
derivative contains

$$
\frac{\partial r_d}{\partial I_d}=1,\qquad
\frac{\partial r_d}{\partial I_u}=-e^{-\tau_i},\qquad
\frac{\partial r_d}{\partial x_i}=I_u\,e^{-\tau_i}\,\Delta_i\,\frac{\partial\alpha_i}{\partial x_i},
\qquad
\frac{\partial S^{\text{rad}}_i}{\partial I_u}=A_u,\qquad
\frac{\partial S^{\text{rad}}_i}{\partial I_d}=-A_d,
$$

with $x_i$ a mass of cell $i$: $\partial\alpha_i/\partial m_{j,i} = \alpha_j\,
(\mathrm{d}m^+_{j,i}/\mathrm{d}m_{j,i})/V_i$, which is $\alpha_j/V_i$ above the
cell's positive-part scale. The mass absorption coefficients are constants, so
$\partial\alpha_i/\partial T_i = 0$. On a moving mesh the same pass gives the
position columns, through $\Delta_i$, through the volume inside $\xi^+$, and
through the face areas. The exposed-face row has $\partial r_e/\partial I_e = 1$
and, on the Kirchhoff path,

$$
\frac{\partial r_e}{\partial x_s} = -\,q_{\text{inc}}\,\frac{\partial a_{\text{eff}}}{\partial x_s},
$$

with $x_s$ the masses and the energy state of the exposed cell: $a_{\text{eff}}$
depends on the composition through the mixing volume fractions, and on the
temperature through $\varepsilon_j(T)$ (a component without its own absorptivity)
and through temperature-dependent densities inside the volume fractions. With an
explicit absorptivity these entries vanish. The tests check exactly this pattern
of zeros: under Kirchhoff's law the exposed cell's mass column reaches the
exposed-face row and under an explicit absorptivity it does not; its temperature
column reaches that row only when the cell's absorptivity depends on temperature
and the condition sets none of its own.

### 8.6.3 The eliminated sensitivities

Eliminating the intensity block (its Schur complement) recovers the
upper-triangular sensitivities of the eliminated model, so the two descriptions
are the same derivative. With
$q^{\text{rad}}_i=I_{\text{up}}(i)\,(1-e^{-\alpha_i\Delta z_i})/\Delta z_i$ and a
state variable $x_k$ of cell $k$ (a temperature or a mass),

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_k} =
\begin{cases}
\big(I_{\text{up}}(i) - q^{\text{rad}}_i\,\Delta z_i\big)\,\dfrac{\partial \alpha_i}{\partial x_i}, & k=i,\\[2.0ex]
-\,q^{\text{rad}}_i\,\Delta z_k\,\dfrac{\partial \alpha_k}{\partial x_k}, & k \text{ between } i \text{ and the exposed face},\\[1.0ex]
0, & k \text{ deeper than } i .
\end{cases}
$$

The **diagonal** term has a clean identity. From
$q^{\text{rad}}_i=I_{\text{up}}(1-e^{-\alpha_i\Delta z_i})/\Delta z_i$,

$$
\frac{\partial q^{\text{rad}}_i}{\partial\alpha_i}
= I_{\text{up}}\,e^{-\alpha_i\Delta z_i}
= I_{\text{up}} - q^{\text{rad}}_i\,\Delta z_i ,
$$

so the diagonal coefficient is exactly
$\big(I_{\text{up}}-q^{\text{rad}}_i\Delta z_i\big)\partial\alpha_i/\partial x_i$.
(The form $I_{\text{up}}(1-(1-e^{-\tau})/\tau)$ is sometimes quoted; it is *not*
equal to this and is incorrect.) The **off-diagonal** term comes from
$\partial I_{\text{up}}(i)/\partial\alpha_k = -I_{\text{up}}(i)\,\Delta z_k$ for
every cell $k$ the beam crosses before cell $i$.

**The Kirchhoff surface-absorptivity chain.** On the Kirchhoff path every
$q^{\text{rad}}_i$ is proportional to $a_{\text{eff}}$ of the exposed cell $s$, so
each eliminated source also depends on that cell's state,

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_s} \mathrel{+}=
\frac{q^{\text{rad}}_i}{a_{\text{eff}}}\,
\frac{\partial a_{\text{eff}}}{\partial x_s}\qquad(i\neq s).
$$

In the row form this whole chain is the single entry $\partial r_e/\partial x_s$
of §8.6.2, carried down the column by the attenuation rows. On a moving mesh the
eliminated sources also depend on every node position above them, through the
$\Delta z_k$ of the transmission; in the row form these are the position entries
of each row's own cell.

### 8.6.4 Verification

The kernel tests (`test/kernels/test_radiation.jl`) compare the colored Jacobian
with the dense forward-mode reference (§15.13) entry by entry to a relative
$10^{-12}$, on the initial state, a hot state, a perturbed state, the zero state
and a negative state, in both energy forms, for a Kirchhoff and an explicit
absorptivity, two absorber laws (a component absorptivity, and an absorptivity
that depends on temperature), with and without gas transport, and on a fixed and a
moving mesh; the pattern contains every non-zero of the dense reference. The
exposed cell's temperature and mass columns are compared with central differences
of the residual, an oracle independent of automatic differentiation. A further
test checks the telescoping and the exact cell averages at consistent
intensities: each cell absorbs $I_{\text{up}}(1-e^{-\tau_i})$ and the cells
together absorb the entering intensity less the transmitted one, to a relative
$10^{-12}$. The residual allocates nothing.

---

## 8.7 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to ours on first use. All three reference
solvers are condensed-phase pyrolysis models; their radiation treatments span the
same idioms Pyrolysis.jl supports, plus a two-flux variant in FDS.

### 8.7.1 ThermaKin2Ds

ThermaKin (Stoliarov & Lyon) uses **in-depth exponential (Beer–Lambert)
absorption** of the external flux together with explicit surface **re-radiation**.
In ThermaKin notation the external radiant flux is $\dot q''_{ex}$ (our
$q_{\text{inc}}$), the absorption coefficient is per-component and forms a
volumetric coefficient (our $\alpha_{\text{eff}}=\sum_j\xi_j\alpha_j$), and the
in-depth source is $\dot q''_{ex}\,\alpha\,e^{-\int\alpha\,\mathrm{d}x}$ — identical
in form to our $q'''=\alpha I$ with $I=I_e e^{-\int\alpha\,\mathrm{d}z}$ (our
`BEER_LAMBERT`). ThermaKin additionally models re-radiation $\varepsilon\sigma T_s^4$
at the surface (often with reflection of a fraction of $\dot q''_{ex}$), which
Pyrolysis.jl supplies via a top `RadiativeBC` term
$F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$ in the surface balance and the
reflected share $(1-a)\,q_{\text{inc}}$ (§8.3.2, §12). ThermaKin also offers a
radiative-conductivity augmentation $k\propto T^3$ to mimic radiative diffusion in
char; Pyrolysis.jl has no separate radiative-diffusion model, and the same
augmentation is a component conductivity law with a $T^3$ term (a function or
`PolynomialProperty` for `k`, §5.5). ThermaKin publishes heats of reaction with
the opposite sign convention ($h>0$ exothermic); this does not touch the radiation
terms (§6.9.2).

### 8.7.2 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) deposits the external radiant flux
**in-depth** with an exponential absorption profile and treats the surface energy
balance with combined radiative–convective exchange. In Gpyro notation the
pre-exponential of any property is $Z$ (our $A_i$), the volume fraction is
$X_\alpha$ (our $v_j$), and the in-depth absorbed flux follows
$\dot q''_e\,\kappa\,e^{-\kappa x}$ — again our `BEER_LAMBERT` form, with Gpyro's
$\kappa$ (the absorption coefficient, distinct from the permeability $K$ in §9)
mapping to our volumetric $\alpha_{\text{eff}}$. Gpyro's effective radiative
conductivity in optically thick char (a temperature-cubed conductivity
contribution) plays the role a $T^3$ term in our $k_j(T)$ would play; Pyrolysis.jl
carries it, if at all, in the mixture conductivity of §5 and §7. Gpyro's surface
emissivity enters the boundary balance just as our
$\varepsilon_{\text{eff}}=\sum_j v_j\varepsilon_j$ does.

### 8.7.3 FDS

The FDS solid-phase model uses a **two-flux (Schuster–Schwarzschild)** in-depth
radiation model for the condensed phase: forward and backward hemispheric
intensities $I^+,I^-$ are propagated with absorption coefficient $\kappa_s$ (our
$\alpha$), and the radiative source in the 1D Fourier conduction equation is the
divergence of the net radiative flux $I^+-I^-$. Pyrolysis.jl's `BEER_LAMBERT`
corresponds to the single-flux limit: only the inward-going beam, with no internal
re-emission. FDS applies surface absorptivity and emissivity through Kirchhoff's
law at the gas–solid interface — the default ($a=\varepsilon_{\text{eff}}$)
Pyrolysis.jl uses for `RadiativeFluxBC`. FDS's gray-band absorption maps to our
gray (wavelength-independent) $\alpha_{\text{eff}}$. The principal difference is
the explicit two-stream bookkeeping: FDS resolves $I^\pm$ with emission in depth,
whereas Pyrolysis.jl resolves the inward beam only (§8.8).

### 8.7.4 Summary table

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS (solid) |
|---|---|---|---|---|
| In-depth absorption | Beer–Lambert $q'''=\alpha I$, face intensities as unknowns | Beer–Lambert | exp. in-depth | two-flux $I^\pm$ |
| Radiative diffusion in depth | none (a $T^3$ term in $k_j(T)$ can carry it) | radiative $k\propto T^3$ | radiative $k$ | two-flux divergence |
| Surface absorptivity | $a_{\text{eff}}$ (Kirchhoff $\varepsilon_{\text{eff}}$ unless a component sets $a_j$) | $\varepsilon$ + reflection | $\varepsilon$ | Kirchhoff |
| Spectral model | gray | gray | gray | gray-band |
| Surface re-radiation | `RadiativeBC` $F\varepsilon\sigma(T_\infty^4-T_s^4)$ | $\varepsilon\sigma T_s^4$ | radiative–convective | $\varepsilon\sigma T_s^4$ |

---

## 8.8 Limitations and Open Issues

The radiation subsystem makes several deliberate simplifications; the following
are the open issues most likely to matter in practice.

1. **No scattering.** There is no scattering coefficient. Strongly scattering
   chars or foams would need $\sigma_s>0$ and, for anisotropy, a higher-order
   closure.

2. **Gray (wavelength-independent) properties.** $\alpha_j$, $\varepsilon_j$ and
   $a_j$ are gray; real polymers and chars are spectrally selective. A band
   extension would need per-band absorption coefficients and one beam per band.

3. **No emission in depth.** `BEER_LAMBERT` carries the incident beam only: the
   medium's own emission, and the radiative transport it drives inside a hot,
   optically thick char, are not modeled; emission is the surface process of
   `RadiativeBC`. Neither the first-order spherical-harmonics approximation (an
   incident-radiation field $G$) nor a two-flux pair $I^\pm$ is implemented.
   Either would be one more block of algebraic unknowns of the state with
   tridiagonal rows (a quasi-steady diffusion equation for $G$, or two opposed
   attenuation-and-emission recurrences), banded like the intensity rows of
   §8.3.3 and differentiated by the same colored automatic differentiation.

4. **One beam, through `:top`, on a planar line.** Only the first
   `RadiativeFluxBC` member of the `:top` condition feeds the beam; a second member
   is dropped silently (§8.3.4). The beam enters through the single `:top` face;
   `BEER_LAMBERT` is refused off planar coordinates and on meshes of more than one
   dimension.

5. **The transmitted beam leaves.** The intensity reaching the back face is lost
   through it (§8.3.3). For a sample on an opaque holder that is a modelling
   error; `sol.extras.radiation_transmitted` measures its size, and no warning is
   issued. Nor is any warning issued for a model that does not match the material
   (§8.1.2).

6. **The surface Newton derivative of a condition-supplied $\varepsilon(T)$.** A
   `RadiativeBC` with its own emissivity law evaluates $\varepsilon(T_s)$ at the
   surface temperature, and the derivative of the surface-balance Newton
   iteration drops $\mathrm{d}\varepsilon/\mathrm{d}T_s$ (§12.2.4). This affects
   only the convergence of the iteration, not its root. The default
   $\varepsilon_{\text{eff}}$ and $a_{\text{eff}}$ are those of the boundary cell
   at its composition and trial temperature, not at $T_s$; they are constants
   within the iteration.

7. **`SURFACE_ABSORPTION` assumes an opaque first cell** (§8.2). Absorbed incident
   radiation enters through the *boundary balance*; the model is accurate only
   when the first cell is optically thick ($\tau_n\gg1$). If the first cell is
   coarse or only marginally opaque, refine the surface cell or switch to
   `BEER_LAMBERT`.
