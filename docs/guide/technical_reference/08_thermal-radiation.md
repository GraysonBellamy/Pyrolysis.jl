# 8. Thermal Radiation

This chapter derives and documents every thermal-radiation model implemented in
Pyrolysis.jl and the way each couples into the energy equation. The solver offers
four mutually exclusive volumetric-radiation treatments — `NO_RADIATION`,
`SURFACE_ABSORPTION`, `BEER_LAMBERT`, and a quasi-steady P1 (first-order
spherical-harmonics) approximation — each with a distinct cost, Jacobian
structure, and physical regime of validity. We give the continuous radiative
transfer equation (RTE), its reduction to the Beer–Lambert and P1 forms, the exact
cell-averaged discretizations, the Marshak/Reflective/Dirichlet boundary closures,
the allocation-free Newton/Thomas solver for the P1 field, and the non-local
(upper-triangular) Jacobian that Beer–Lambert induces. Every governing form,
discrete update, constant, sign, factor, guard, and limiter below has been verified
against the implementation; where a docstring, textbook form, or knowledge-base
note diverges from the code, **the code is authoritative** and the divergence is
called out explicitly.

Throughout, the spatial coordinate is $z$ (axial, through-thickness), with $z=0$
the bottom/substrate (boundary id 2, `:bottom`) and $z=L$ the exposed/top surface
(boundary id 1, `:top`). Heat enters from the top. The flux sign convention is
**positive flux = transport in the $+z$ direction**; the discrete divergence
$(F_R-F_L)A/V$ is positive when the quantity flows *out* of a cell. The cell index
is $i$ ($1\le i\le n$), the component index is $j$ ($1\le j\le N_C$). The
Stefan–Boltzmann constant is $\sigma = 5.670374419\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}$
(`STEFAN_BOLTZMANN`, src/Materials/constants.jl:9).

---

## 8.1 The Radiation Model Family and its Traits

### 8.1.1 The four models

Volumetric radiation in Pyrolysis.jl is selected at `solve` time by a single
enumeration value passed as `radiation_model`. The enum is defined in
src/Materials/radiation_types.jl:82–87:

```julia
@enum RadiationModel begin
    NO_RADIATION = 0
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
    P1_QUASI_STEADY = 3
end
```

The models are **mutually exclusive** — exactly one is active for the whole domain
in a given run — and they differ in where the incident external radiation is
deposited and how the material's own re-emission is treated:

1. **`NO_RADIATION`.** No volumetric radiation source ($Q_{\text{rad}}=0$).
   External radiative *boundary* conditions (`RadiativeBC` surface re-radiation,
   `RadiativeFluxBC` incident flux) still act at the surface through the boundary
   energy balance (§8.8, and §12, Boundary & Initial Conditions). Use for
   optically transparent media ($\tau \ll 0.1$) or when in-depth absorption is
   negligible.

2. **`SURFACE_ABSORPTION`.** The incident external radiation is absorbed
   entirely at the exposed surface. The absorbed flux is **not** deposited as a
   volumetric source; it enters the energy balance as a *boundary heat flux*
   through the surface-temperature Newton solve (§8.2, §8.8). Component absorption
   coefficients $\alpha_j$ are ignored. Cost $O(1)$; tridiagonal Jacobian. This is
   the recommended choice for opaque solids with first-cell optical thickness
   $\tau\gg 1$.

3. **`BEER_LAMBERT`.** In-depth exponential absorption following the Beer–Lambert
   law, $I(z)=I_{\text{surface}}\exp\!\big(-\int_L^z\alpha_{\text{eff}}\,\mathrm{d}z'\big)$,
   with volumetric source $q'''=\alpha_{\text{eff}}I$. Radiation penetrates over the
   absorption length $1/\alpha_{\text{eff}}$. Cost $O(n)$; the intensity at cell $i$
   depends on all cells above it, producing a **non-local, upper-triangular**
   Jacobian coupling (§8.10). Appropriate for semi-transparent polymers and foams
   ($0.1\lesssim\tau\lesssim 10$).

4. **`P1_QUASI_STEADY`.** The P1 diffusion approximation, solving an auxiliary
   incident-radiation field $G$ from a quasi-steady diffusion equation, including
   *volumetric re-emission* $4\sigma\alpha T^4$. Cost $\approx O(10n)$ (a small
   Newton loop per step); the $G$ equation is tridiagonal. Recommended when
   internal re-radiation matters ($T\gtrsim 800$ K, optically thick char). P1 lives
   in the `Experimental` module and runs through a *separate residual* that
   re-solves the quasi-steady $G$ field into its radiation cache at every
   evaluation — $G$ is not part of the ODE state (§8.12).

### 8.1.2 Model traits

Three trait queries (src/Materials/radiation_types.jl:109–149) let the assembly and
solver configuration branch at compile time without inspecting the model in hot
loops:

| Trait | True for | Meaning |
|---|---|---|
| `requires_absorption_coefficients(model)` | `BEER_LAMBERT`, `P1_QUASI_STEADY` | the model needs the per-component mass absorption $\alpha_j$ to form $\alpha_{\text{eff}}$ |
| `has_nonlocal_jacobian_coupling(model)` | `BEER_LAMBERT` only | the Jacobian is no longer tridiagonal (upper-triangular fill) |
| `requires_radiation_cache(model)` | `P1_QUASI_STEADY` only | a `RadiationCache` for $G$, $D_{\mathrm{rad}}$, and Newton work arrays is allocated |

The performance/structure summary (radiation_types.jl:74–80):

| Model | Loop cost | Jacobian pattern | Extra memory |
|---|---|---|---|
| `NO_RADIATION` | $O(1)$ | tridiagonal | baseline |
| `SURFACE_ABSORPTION` | $O(1)$ | tridiagonal | baseline |
| `BEER_LAMBERT` | $O(n)$ | upper-triangular | $+O(n^2)$ (dense fill) |
| `P1_QUASI_STEADY` | $O(10n)$ | tridiagonal | $+O(n)$ |

> **Implementation note.** The dispatch `compute_radiation_source!`
> (`q_rad, state, mesh, cache, I_surface, model, apply_kirchhoff_absorptivity`,
> src/Physics/radiation.jl:118–141)
> *zeros* `q_rad` for `NO_RADIATION`, `P1_QUASI_STEADY`, **and**
> `SURFACE_ABSORPTION`, and only calls `compute_beer_lambert_source!` for
> `BEER_LAMBERT`. P1 is added on a separate residual path (§8.12). The boundary
> energy balance accounts for the surface-absorbed radiation in `NO_RADIATION` and
> `SURFACE_ABSORPTION` modes.

---

## 8.2 Surface Absorption

### 8.2.1 Surface absorption is a boundary flux, not a volumetric source

In the assembled RHS, the `SURFACE_ABSORPTION` model contributes *zero*
volumetric source — `compute_radiation_source!` zeroes `q_rad`
(src/Physics/radiation.jl) — and the incident radiation is folded into the
surface energy balance through the `RadiativeFluxBC` dispatch
(boundary_temperature.jl):

$$
q_{\text{BC}}(T_s) = \Phi_{\text{incident}}\,\alpha
\quad\text{for }\texttt{NO\_RADIATION}\text{ or }\texttt{SURFACE\_ABSORPTION},
$$

where $\Phi_{\text{incident}}=I_{\text{ext}}(t)$ is the incident intensity
[W·m⁻²] and $\alpha$ is the surface absorptivity. By Kirchhoff's law at radiative
equilibrium the default absorptivity equals the effective surface emissivity,
$\alpha=\varepsilon_{\text{eff}}=\sum_j v_j^{\text{solid}}\varepsilon_j$ (clamped
to $[0,1]$); an explicit `absorptivity` overrides this. The temperature derivative
of this contribution is $\mathrm{d}q_{\text{BC}}/\mathrm{d}T_s=0$ (absorbed
incident flux is independent of $T_s$). Surface re-radiation loss is supplied
separately by a `RadiativeBC` term, $q=F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$,
whose $T_s$-dependence does enter the Newton derivative (§12).

The textbook alternative — depositing the absorbed flux as a volumetric source
$I/\Delta z_N$ in the top cell — is deliberately **not** implemented. The
boundary-flux route is the better construction: it lets the surface-temperature
Newton solve resolve the (possibly large) jump between $T_s$ and the cell-center
temperature $T_N$ consistently with convection and re-radiation; for an opaque
material with a well-resolved first cell the two constructions would be
numerically close anyway.

### 8.2.2 Assumptions and validity

Surface absorption presumes the material is optically thick over the first cell,
$\tau_N=\alpha_{\text{eff}}\Delta z_N\gg 1$, so that essentially no radiation
penetrates beyond cell $N$. It is $O(1)$, adds no Jacobian non-locality (the
tridiagonal conduction stencil is preserved), and is the recommended choice for
charring solids and fast cone-calorimeter analysis.

---

## 8.3 The Beer–Lambert In-Depth Model

### 8.3.1 Continuous form

Beer–Lambert follows from the 1D radiative transfer equation for a collimated beam
along $z$ with absorption coefficient $\alpha_{\text{eff}}(z)$ and negligible
scattering/emission of the *incident* beam:

$$
\frac{\partial I}{\partial z} = -\alpha_{\text{eff}}\, I .
$$

With the beam entering at the top surface $z=L$ and travelling into depth (toward
$z=0$), integration gives exponential attenuation:

$$
I(z) = I_{\text{surface}}\,\exp\!\left(-\int_{L}^{z}\alpha_{\text{eff}}(z')\,\mathrm{d}z'\right).
$$

The local energy deposited per unit volume is the absorbed power density,

$$
q'''(z) = \alpha_{\text{eff}}(z)\,I(z)\qquad[\mathrm{W\,m^{-3}}].
$$

The effective absorption coefficient is the composition-weighted volumetric
extinction,

$$
\alpha_{\text{eff}} = \sum_j \xi_j\,\alpha_j,
$$

where $\xi_j$ [kg·m⁻³] is the mass concentration and $\alpha_j$ [m²·kg⁻¹] is the
per-component mass absorption coefficient (`effective_absorption`,
src/Properties/effective_properties.jl:661–672). The product $\xi_j\alpha_j$ has
units m⁻¹, so $\alpha_{\text{eff}}$ is a *volumetric* coefficient; it is stored in
the property cache as `cache.α[i]` (see §5, Effective Properties & Mixing Rules).
The optical thickness of the whole layer is
$\tau = \int_0^L \alpha_{\text{eff}}\,\mathrm{d}z$.

### 8.3.2 The exact cell-averaged discretization

A naive evaluation $q'''_i \approx \alpha_i I(z_i^c)$ at the cell center is only
first-order accurate and can over- or under-deposit energy depending on
$\tau_i$. Pyrolysis.jl instead uses the **exact cell-integrated** source: integrate
$q'''=\alpha I$ over a cell of constant $\alpha_i$ and entrance intensity
$I_{\text{in}}$,

$$
\int_{\text{cell }i} \alpha_i\,I\,\mathrm{d}z
= \int_0^{\Delta z_i}\!\alpha_i\, I_{\text{in}}\, e^{-\alpha_i s}\,\mathrm{d}s
= I_{\text{in}}\big(1 - e^{-\alpha_i \Delta z_i}\big),
$$

so the *cell-averaged* volumetric source is the absorbed power divided by the cell
thickness:

$$
\boxed{\;q^{\text{rad}}_i = I_{\text{in}}\,\frac{1 - e^{-\tau_i}}{\Delta z_i},
\qquad \tau_i \equiv \alpha_i\,\Delta z_i\;}
$$

This is exact for piecewise-constant $\alpha$ and conserves energy by
construction: the power absorbed in cell $i$ is $q^{\text{rad}}_i\,\Delta z_i =
I_{\text{in}}(1-e^{-\tau_i})$, exactly the drop in beam power across the cell.

The flux that actually crosses the surface is $\alpha_{\text{surf}}\,
I_{\text{surface}}$, where $\alpha_{\text{surf}}$ is the surface
transmissivity: on the default **Kirchhoff path** (the `RadiativeFluxBC`
carries no explicit `absorptivity`) it is the *live* surface cell's effective
emissivity $\varepsilon_{\text{eff}}(T_N,\xi_{\cdot,N})$ — absorption and
surface emission then use the same coefficient (Kirchhoff's law
$\alpha=\varepsilon$, detailed balance), and the remaining
$(1-\varepsilon_{\text{eff}})I_{\text{surface}}$ is reflected and never
enters the medium. With an **explicit BC absorptivity**, $\alpha(t)$ is folded
into the extracted intensity upstream and the kernel does not rescale
(`apply_kirchhoff_absorptivity = false`). The entrance intensity is then the
cumulative attenuation from the surface down to the top face of cell $i$:

$$
I_{\text{in}}(i) = \alpha_{\text{surf}}\,I_{\text{surface}}\prod_{j=i+1}^{N} e^{-\tau_j}
= \alpha_{\text{surf}}\,I_{\text{surface}}\,\exp\!\left(-\sum_{j>i}\alpha_j\Delta z_j\right).
$$

The transmitted (un-absorbed) fraction leaving the bottom is
$\alpha_{\text{surf}}\,I_{\text{surface}}\,e^{-\tau}$; the diagnostic
`transmitted_radiation` (radiation.jl) reports the attenuation factor for the
intensity it is given.

### 8.3.3 Backward sweep and guards

Because $I_{\text{in}}(i)$ depends only on cells *above* $i$ (larger index, nearer
the surface), the implementation sweeps from the top cell $N$ down to cell 1,
carrying a running intensity `I_current`. Forward iteration would require
not-yet-computed information; the backward sweep is the natural causal order.

```julia
I_current = oneunit(eltype(q_rad)) * I_surface
if apply_kirchhoff_absorptivity
    I_current *= cache.ε[n_cells_val]   # Kirchhoff: α_surf = ε_eff (live surface cell)
end
A = state.A_section[]
@inbounds for i = n_cells_val:-1:1
    α  = cache.α[i]
    Δz = state.volume[i] / A
    if α > 1e-20
        τ = α * Δz
        absorption_factor = 1 - exp(-τ)
        q_rad[i] = I_current * absorption_factor / Δz
        I_current *= exp(-τ)
    else
        q_rad[i] = zero(eltype(q_rad))
    end
end
```
(`compute_beer_lambert_source!`, src/Physics/radiation.jl.)

Three details are load-bearing:

- **Cell thickness is geometric.** $\Delta z_i = V_i / A$ from `state.volume[i]`
  and the live cross-section `state.A_section[]`, so the source remains exact under
  ALE mesh motion and variable cross-section (the lateral-shrinkage area cancels
  between $V_i$ and $A$; see §10, Volume Change & Lateral Shrinkage, and §11, The
  ALE Framework).
- **Transparency guard.** Cells with $\alpha\le 10^{-20}$ get $q^{\text{rad}}_i=0$
  and do not attenuate the beam, preventing `exp(0)`-degenerate work and spurious
  absorption in radiation-free regions.
- **Zero-flux short-circuit.** If $I_{\text{surface}}\le 0$, the whole array is
  zeroed and the routine returns (radiation.jl:60–63).

The source is **positive = heating**: it is added to the energy RHS as
$+q^{\text{rad}}_i/\rho c_p^{\text{eff}}$ (§8.8). Energy conservation holds in the
strong sense $\sum_i q^{\text{rad}}_i\,V_i = \alpha_{\text{surf}}\,
I_{\text{surface}}\,(1-e^{-\tau_{\text{tot}}})\,A \le \alpha_{\text{surf}}\,
I_{\text{surface}}\,A$, with equality only if the layer is fully opaque
($e^{-\tau}\to 0$) — consistent with `SURFACE_ABSORPTION`, which applies
$\varepsilon_{\text{eff}}\,I_{\text{surface}}$ at the face (switching models
therefore does not change the absorbed power in the opaque limit).

### 8.3.4 Where $I_{\text{surface}}$ comes from

The surface intensity fed to Beer–Lambert is the incident flux carried by the
top-boundary `RadiativeFluxBC`. `extract_radiative_flux_from_bcs`
(src/BoundaryConditions/boundary_conditions.jl) returns
`(I, has_explicit_absorptivity)`: with an explicit BC `absorptivity` it returns
$I_{\text{ext}}\alpha(t)$ and `true`; otherwise the bare incident value and
`false` — the live surface-cell $\varepsilon_{\text{eff}}$ is then applied at
injection inside `compute_beer_lambert_source!` (the Kirchhoff path). `solve()`
derives the static `radiation_kirchhoff` flag from the BC once and threads it
through `to_problem_def` to both the residual and the structured Jacobian, so
the two conventions can never mix.

Beer–Lambert models **top-incident** radiation only: the sweep starts at the
exposed surface, and the per-BC face flux of a `RadiativeFluxBC` is zeroed at
*every* boundary under this model. A `RadiativeFluxBC` on the **bottom**
boundary therefore cannot be honored, and `solve()` rejects the combination
with an `ArgumentError` rather than silently dropping the back-side heating
(src/Solver/solve.jl); use `SURFACE_ABSORPTION` (face deposit at either
boundary) or a `HeatFluxBC` for back-side radiant heating.

The solver feeds this flux to the ODE by re-sampling it at **every** RHS
evaluation — `ws.I_radiation = radiation_flux(t)` — whenever the active model
consumes it (`SURFACE_ABSORPTION` or `BEER_LAMBERT`; src/Solver/solve.jl).
There is no constant-flux detection: one user-function call per RHS is
negligible against the residual cost, and unconditional re-sampling guarantees
that a ramped or modulated cone flux can never be frozen at its $t=0$ value.

> **Implementation note.** Snapshotting a time-varying incident flux at $t=0$
> would silently freeze a cone-heater ramp — a silent error. The unconditional
> per-RHS re-sampling forecloses this failure mode by construction; a regression
> with $q_0(1+\sin 2\pi t)$ guards it (§18).

### 8.3.5 Worked numbers

For a 10 mm semi-transparent slab with $\alpha_{\text{eff}}=100\ \mathrm{m^{-1}}$,
uniform $\Delta z=1\ \mathrm{mm}$, $I_{\text{surface}}=50\ \mathrm{kW\,m^{-2}}$,
$\varepsilon_{\text{eff}}=0.9$ (Kirchhoff path): the flux entering the medium
is $0.9\times50 = 45\ \mathrm{kW\,m^{-2}}$. With $\tau_i=0.1$ the top cell
absorbs $q_1=45\,000\,(1-e^{-0.1})/0.001 \approx 4.28\ \mathrm{MW\,m^{-3}}$;
the beam entering the second cell is $45\,000\,e^{-0.1}=40.7\ \mathrm{kW\,m^{-2}}$,
giving $q_2\approx 3.88\ \mathrm{MW\,m^{-3}}$, and so on. The total optical
thickness is $\tau=\alpha L=1$, so $e^{-\tau}\approx 0.37$ of the entering flux
is transmitted through the slab — Beer–Lambert is the right tool in this regime,
whereas surface absorption would (incorrectly) deposit everything at the surface.

---

## 8.4 The P1 (First-Order Spherical-Harmonics) Approximation

### 8.4.1 From the RTE to a diffusion equation for $G$

For an absorbing–emitting (and, in principle, isotropically scattering) gray
medium, expanding the directional intensity in spherical harmonics and truncating
at first order (P1) yields a diffusion equation for the **incident radiation**
$G=\int_{4\pi} I\,\mathrm{d}\Omega$ [W·m⁻²]:

$$
\boxed{\;\nabla\!\cdot\!\big(D_{\mathrm{rad}}\,\nabla G\big) = \alpha\,G - 4\sigma\alpha T^4\;}
$$

with the P1 radiation diffusion coefficient

$$
D_{\mathrm{rad}} = \frac{1}{3(\alpha+\sigma_s)}\;\xrightarrow{\sigma_s\to0}\;\frac{1}{3\alpha}
\qquad[\mathrm{m}],
$$

where $\alpha$ is the volumetric absorption coefficient, $\sigma_s$ is the
scattering coefficient (taken $\approx 0$ for pyrolysis), and $\sigma$ is the
Stefan–Boltzmann constant. The left side is radiative *diffusion* (transport by
repeated absorption/re-emission); the right side balances removal from the field
by absorption ($\alpha G$) against blackbody emission ($4\sigma\alpha T^4$).

> **Implementation note.** `compute_radiation_diffusion_coefficient`
> (src/Experimental/radiation_p1.jl) returns
> $1/\big(3\max(\alpha+\sigma_s,\,10^{-10})\big)$ with $\sigma_s=0$ by default —
> the transparent limit is clamped to the large finite value
> $1/(3\cdot 10^{-10})$ rather than zeroed, because $D=0$ would delete the $G$
> row of the quasi-steady Newton system wherever $\alpha=0$ (singular
> tridiagonal matrix). This is the **single** canonical guard: the one
> `update_radiation_properties!` method routes every $D_{\mathrm{rad}}$ through
> it (see open issue 8 in §8.13).

### 8.4.2 The net volumetric source and its sign

Integrating the RTE over solid angle gives the radiative-flux divergence
$\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}} = \alpha(4\sigma T^4 - G)$. The energy
equation carries $-\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}}$, so the **net
volumetric radiation source** that heats the matrix is

$$
\boxed{\;S_{\mathrm{rad}} = \alpha\,(G - 4\sigma T^4)\;}
$$

with the convention

- $S_{\mathrm{rad}}>0$: net absorption ($G>4\sigma T^4$) — the material **heats**;
- $S_{\mathrm{rad}}<0$: net emission ($T^4$ exceeds $G$) — the material **cools**;
- $S_{\mathrm{rad}}=0$: local radiative equilibrium.

This is what the active P1 residual uses: `compute_radiation_flux_divergence!`
(src/Experimental/radiation_fluxes.jl) computes the **algebraic** source
$S_{\mathrm{rad}}[i] = \alpha_i(G_i - 4\sigma T_i^4)$ cell-by-cell and stores it in
`rad_cache.S_rad`, and the residual adds $\mathrm{d}T_i/\mathrm{d}t \mathrel{+}=
S_{\mathrm{rad}}[i]/\rho c_p^{\text{eff}}_i$ (`compute_p1_residual!`,
RadiationP1.jl). This algebraic identity is equivalent to $-\nabla\!\cdot\!
\mathbf{q}_{\mathrm{rad}}$ *only when the discrete $G$ field exactly satisfies the
P1 equation*; the Newton solve below drives it there.

> **Two divergence forms — only one is on the active path.** An alternative
> `apply_radiation_divergence!` (radiation_fluxes.jl) instead forms the
> *flux* divergence $S=(q_{\mathrm{rad}}[R]-q_{\mathrm{rad}}[L])A/V$ from the
> diffusive face fluxes $q_{\mathrm{rad}}=-D_{\mathrm{rad}}\nabla G$ and adds
> $\mathrm{d}T_i/\mathrm{d}t \mathrel{+}= -S/\rho c_p$. The minus sign there is
> correct because $\nabla\!\cdot\!(D\nabla G) = -\nabla\!\cdot\!\mathbf{q}_{\mathrm{rad}}$.
> The production P1 residual (`compute_p1_residual!`, RadiationP1.jl) calls
> the *algebraic* `compute_radiation_flux_divergence!`, not this
> flux-divergence variant; the two agree at the converged $G$.

### 8.4.3 The quasi-steady assumption

Radiation equilibrates orders of magnitude faster than conduction or reaction
($\tau_{\mathrm{rad}}\sim 10^{-12}$ s vs. $\tau_{\mathrm{cond}}\sim 10^{2}$ s), so
the field is treated as instantaneous:

$$
\frac{\partial G}{\partial t}\approx 0 .
$$

Operationally, $G$ is **not** a state variable at all: it lives in the
`RadiationCache` and is re-solved as a nonlinear algebraic system by the Newton
iteration of §8.7 at *every* residual evaluation (§8.12). The temperature
non-linearity $T^4$ makes the $G$ system nonlinear even though it is linear in
$G$ alone at fixed $T$.

---

## 8.5 P1 Volumetric Source Terms

The two pieces of the right-hand side of the P1 equation are computed by
`compute_radiation_source_terms!` (src/Experimental/radiation_p1.jl:98–130):

**Emission** [W·m⁻³]:

$$
Q_{\mathrm{em}}[i] = 4\sigma\,\alpha_i\,T_i^4 .
$$

**Absorption** [W·m⁻³]:

$$
Q_{\mathrm{abs}}[i] = \alpha_i\,G_i .
$$

The net source used in the energy equation is $S_{\mathrm{rad}} = Q_{\mathrm{abs}} -
Q_{\mathrm{em}} = \alpha(G - 4\sigma T^4)$ (§8.4.2), consistent with
`net_radiation_source` (src/Properties/radiation_cache.jl:255–257).

> **No emissivity factor in the volumetric emission — the code is authoritative.**
> Several references (and even the docstrings of `compute_radiation_source_terms!`
> and the `RadiationCache.Q_emission` field) write emission as $4\sigma
> \varepsilon\alpha T^4$. The executed code (radiation_p1.jl:122) uses **blackbody**
> emission $4\sigma\alpha T^4$ with **no** $\varepsilon$ factor. The physical
> justification is that in the P1 approximation the interior is a gray body in local
> thermodynamic equilibrium ($\varepsilon=1$ locally): a control volume that absorbs
> $\alpha G$ must, at equilibrium $G=4\sigma T^4$, emit exactly $4\sigma\alpha T^4$,
> giving $S_{\mathrm{rad}}=0$. Surface emissivity $\varepsilon<1$ enters **only** at
> the boundaries through the Marshak condition (§8.6), where it represents partial
> reflection of the diffuse field, not reduced volumetric emission. Writing
> $\varepsilon$ in the interior emission would break the equilibrium identity.

The $\alpha$ used here is read from `rad_cache.α` when populated (and positive),
otherwise from `prop_cache.α` (radiation_p1.jl:112–114); both hold the same
volumetric $\alpha_{\text{eff}}=\sum_j\xi_j\alpha_j$. Caching $\alpha$ on the
radiation cache avoids recomputing it in the source, diffusion-coefficient, and
Jacobian kernels.

---

## 8.6 P1 Boundary Conditions

P1 boundary conditions act on the radiation field $G$, not on $T$, and are stored
in a type-stable `RadiationBCSet{B1,B2}` with `top` and `bottom` slots
(src/Experimental/radiation_boundary.jl:259–277), accessed by
`get_radiation_bc(set, Val(1))` (top) / `Val(2)` (bottom). The **default** set is
a Marshak condition at the top and a reflective condition at the bottom
(`default_radiation_bcs`, radiation_boundary.jl:232–245).

### 8.6.1 Marshak condition (exposed surface)

The Marshak (mixed) condition closes the P1 half-range fluxes at a gray,
diffuse surface of emissivity $\varepsilon$. With $q_{\text{in}} = G/4 + q/2$
and $q_{\text{out}} = G/4 - q/2$ (q = net flux along the inward normal), the
inward half-range flux just inside the surface equals the transmitted external
radiosity plus the reflected outgoing half-range flux,

$$
\frac{G}{4} + \frac{q}{2}
= \varepsilon\big(I_{\text{ext}} + \sigma T_{\text{amb}}^4\big)
+ (1-\varepsilon)\Big(\frac{G}{4} - \frac{q}{2}\Big).
$$

Solving for the net flux into the domain gives the relation implemented by
`radiation_bc_flux` (src/Experimental/radiation_p1.jl):

$$
\boxed{\;q_{\text{BC}} = c\,\big(G_{\text{ext}} + 4 I_{\text{ext}} -
G_{\text{surf}}\big),\qquad c = \frac{\varepsilon}{2(2-\varepsilon)}\;}
$$

with

- $I_{\text{ext}}=I_{\text{external}}(t)$ — external incident radiation [W·m⁻²]
  (cone heater, radiant panel, flame). It carries the half-range weight
  $4c = 2\varepsilon/(2-\varepsilon)$ — **the same** weight the ambient term
  carries ($c\,G_{\text{ext}} = c\cdot 4\sigma T_{\text{amb}}^4 \equiv
  4c\cdot\sigma T_{\text{amb}}^4$): 2.0 at $\varepsilon=1$, 1.64 at
  $\varepsilon=0.9$. A bare weight of 1 on $I_{\text{ext}}$ under-injects the
  external flux by that factor (for a cold optically thick absorber at
  $\varepsilon = 0.9$ the self-consistent absorbed fraction is
  $4c/(1+\sqrt{3}c) \approx 0.96$ with the correct weight versus
  $1/(1+\sqrt{3}c) \approx 0.59$ with weight 1);
- $G_{\text{ext}} = 4\sigma T_{\text{amb}}^4$ — equilibrium incident field of the
  surroundings at the ambient temperature $T_{\text{amb}}(t)$;
- $G_{\text{surf}}$ — the P1 incident field at the surface;
- $\varepsilon=\varepsilon(T_s)$ — surface emissivity (constant, $\varepsilon(T)$,
  or BC-supplied).

```julia
c     = _marshak_exchange_coefficient(bc, T_surf)        # ε/(2(2−ε))
flux  = c * (_marshak_driving_field(bc, t) - G_surf)     # driving = 4σT_amb⁴ + 4I_ext
```

This formulation satisfies the required limits: at equilibrium
($G=4\sigma T^4$, $T_{\text{amb}}=T_{\text{surf}}$, $I_{\text{ext}}=0$) the flux
vanishes; the external and ambient sources enter with consistent half-range
weights; and a fraction $(1-\varepsilon)$ of the interior radiation is reflected
back into the domain.

**Half-cell-corrected application.** The solver applies the relation at the
boundary-cell *center*, eliminating the face value $G_{\text{face}}$ through
the discrete diffusive flux $q = D_{\mathrm{rad}}(G_{\text{face}} - G_N)/\delta$
($\delta$ = center-to-face distance). The result is the same relation with the
series-resistance coefficient

$$
q_{\text{BC}} = c_{\text{eff}}\,\big(G_{\text{ext}} + 4I_{\text{ext}} - G_N\big),
\qquad c_{\text{eff}} = \frac{1}{1/c + \delta/D_{\mathrm{rad}}},
$$

(`_marshak_effective_coefficient`), which removes the $O(\delta)$ error of a
cell-centered evaluation. The flux application (`_apply_marshak_bc!`) and the
quasi-steady Newton Jacobian use the identical $c_{\text{eff}}$.
$\varepsilon(T)$ and $\alpha$ are evaluated at the boundary-cell temperature —
the surface-Newton $T_s$ is not available inside the radiation solve.

> **The $\varepsilon/(2(2-\varepsilon))$ coefficient is the correct Marshak form**
> (Modest, *Radiative Heat Transfer*, 3rd ed., Eq. 16.48), **not** $\tfrac{1}{2}\varepsilon$.
> A common simplification, $0.5\varepsilon$, silently *underestimates* radiative
> exchange for $\varepsilon<1$ (the implemented coefficient is
> `_marshak_exchange_coefficient`, radiation_p1.jl). For $\varepsilon=1$ both
> reduce to $\tfrac12$; the discrepancy grows as $\varepsilon$ falls (e.g. for
> $\varepsilon=0.5$ the correct coefficient is $0.5/(2\cdot1.5)=1/6$ vs. $0.25$).

The flux is applied at the boundary face with the sign-orientation convention of
the divergence operator (radiation_p1.jl:372–400): at a **right** boundary
(`right_id==0`, the top) the stored face flux is $-q_{\text{BC}}$; at a **left**
boundary (the bottom) it is $+q_{\text{BC}}$. Because the `RadiationCache.q_rad`
array is `Float64`, the Marshak flux is stored via `ForwardDiff.value` to strip any
dual part during the Float64-only P1 solve (radiation_p1.jl:368–369, 391–392).

The Marshak BC is constructed as `MarshakRadiationBC(; I_external, T_ambient=300.0,
emissivity=0.9)` (radiation_boundary.jl:79–116) with validation that
$I_{\text{external}}\ge 0$, $T_{\text{ambient}}>0$, and $\varepsilon\in[0,1]$. Any of
the three fields may be a constant or a callable of $t$ (or $T$ for emissivity).

### 8.6.2 Reflective condition (insulated / symmetry surface)

`ReflectiveRadiationBC` (radiation_boundary.jl:153) imposes a zero normal gradient,

$$
\frac{\partial G}{\partial n} = 0,
$$

realized by setting the boundary face flux to zero (`_apply_reflective_bc!`,
radiation_p1.jl:402–408). It is the natural choice for an insulated back face, a
symmetry plane, or a reflective substrate. It is the default bottom BC.

### 8.6.3 Dirichlet condition (prescribed $G$)

`DirichletRadiationBC(; G)` (radiation_boundary.jl:201–206) prescribes the
boundary intensity directly. The current implementation is a placeholder: the
boundary face flux is set to zero (`_apply_dirichlet_bc!`, radiation_p1.jl:410–417;
`radiation_bc_flux` returns `zero(t)`), so $G$ is effectively pinned only through
the interior gradient. It is used mainly for validation against analytic solutions
and is rarely used in practice; Marshak is preferred because it couples the
external source to surface emission physically.

---

## 8.7 P1 Numerics: Diffusion Faces, Fluxes, and the Newton/Thomas Solve

### 8.7.1 Face diffusion coefficient (harmonic mean)

Interior face values of $D_{\mathrm{rad}}$ use a **distance-weighted harmonic
mean** (`update_radiation_properties!`, radiation_p1.jl:178–202):

$$
D_{\mathrm{rad,face}} = \frac{d_L + d_R}{\dfrac{d_L}{D_{\mathrm{rad},L}} + \dfrac{d_R}{D_{\mathrm{rad},R}}},
$$

where $d_L,d_R$ are the distances from the left/right cell centers to the face. The
harmonic mean enforces flux continuity across the interface (resistance in series),
exactly as for thermal conductivity (§7, Heat Conduction) and the transport
coefficients (§5). At a boundary face the value of the adjacent interior cell is
used.

### 8.7.2 Interior diffusive flux

The radiative flux at an interior face is the discrete diffusive flux
(radiation_p1.jl:245–257):

$$
q_{\mathrm{rad}}[f] = -D_{\mathrm{rad,face}}\,\frac{G_R - G_L}{\Delta x},
\qquad \Delta x = |z_R^c - z_L^c|,
$$

with the standard sign convention (positive flux toward $+z$). Boundary faces are
overwritten by the BC handlers of §8.6 in `_apply_all_radiation_bcs!`
(radiation_p1.jl:265–296).

### 8.7.3 The quasi-steady residual

`solve_quasi_steady_radiation!` (src/Experimental/radiation_fluxes.jl:208–414)
solves the discrete P1 equation for $G$ by Newton–Raphson. The per-cell residual
combines the algebraic source terms with the discrete flux divergence
(radiation_fluxes.jl:266–281):

$$
R_i = Q_{\mathrm{abs}}[i] - Q_{\mathrm{em}}[i] + S_{\mathrm{rad,div}}[i],
\qquad
S_{\mathrm{rad,div}}[i] = \big(q_{\mathrm{rad}}[R] - q_{\mathrm{rad}}[L]\big)\frac{A}{V_i}.
$$

This is the discrete form of $R = \alpha G - 4\sigma\alpha T^4 -
\nabla\!\cdot\!(D_{\mathrm{rad}}\nabla G)$: since $q_{\mathrm{rad}}=-D_{\mathrm{rad}}
\nabla G$, the flux divergence $S_{\mathrm{rad,div}} = -\nabla\!\cdot\!(D_{\mathrm{rad}}
\nabla G)$, and $R = (Q_{\mathrm{abs}}-Q_{\mathrm{em}}) + S_{\mathrm{rad,div}}$.

### 8.7.4 The tridiagonal Jacobian

Differentiating $R_i$ with respect to the unknowns $G$ gives a tridiagonal Jacobian
(radiation_fluxes.jl:290–378). With
$\mathrm{coeff}_L = \tfrac{A}{V_i}\tfrac{D_{\mathrm{rad},L}}{\Delta x_L}$ and
$\mathrm{coeff}_R = \tfrac{A}{V_i}\tfrac{D_{\mathrm{rad},R}}{\Delta x_R}$ for
interior faces,

$$
J_{ii} = \alpha_i + \mathrm{coeff}_L + \mathrm{coeff}_R,\qquad
J_{i,i-1} = -\mathrm{coeff}_L,\qquad
J_{i,i+1} = -\mathrm{coeff}_R .
$$

At a Marshak boundary the boundary flux contributes $\partial q_{\text{BC}}/\partial
G_N$ to the diagonal. The applied flux is the half-cell-corrected form
$q_{\text{BC}} = c_{\text{eff}}(G_{\text{ext}} + 4I_{\text{ext}} - G_N)$
(§8.6.1), so the assembled Jacobian uses the **same** effective coefficient as
the residual,

$$
c_{\text{eff}} = \frac{1}{1/c_M + \delta/D_{\mathrm{rad}}},\qquad
c_M(T)=\frac{\varepsilon(T)}{2(2-\varepsilon(T))},
$$

with $\delta$ the boundary-cell center-to-face distance, evaluated from the
boundary condition's emissivity law at the boundary-cell temperature
(`_marshak_effective_coefficient`, radiation_p1.jl):

$$
J_{NN}\mathrel{+}= \frac{A}{V_N}\,c_{\text{eff}}\quad\text{(top Marshak)},
\qquad
J_{11}\mathrel{+}= \frac{A}{V_1}\,c_{\text{eff}}\quad\text{(bottom Marshak)}.
$$

This keeps the Newton matrix exactly consistent with the Marshak residual for
non-blackbody boundaries. Temperature derivatives of emissivity are not included
because this Jacobian is with respect to the radiation unknowns $G$, while the
thermal state is held fixed during the quasi-steady P1 solve.

### 8.7.5 Allocation-free Thomas (TDMA) solve

The tridiagonal system $J\,\Delta G = -R$ is solved by the Thomas algorithm
(`thomas_solve!`, radiation_fluxes.jl:48–74) using pre-allocated diagonals
(`newton_dl`, `newton_d`, `newton_du`), a residual buffer `newton_R`, an update
buffer `newton_ΔG`, and a work array `newton_d_work` that preserves the original
diagonal across iterations (all on the `RadiationCache`). The solve is therefore
allocation-free and $O(n)$ per Newton step. Thomas is stable without pivoting for
the diagonally dominant P1 system.

### 8.7.6 The Newton loop, scaling, and clipping

Per iteration (radiation_fluxes.jl:248–406):

1. recompute the radiation fluxes from the current $G$ (`compute_radiation_fluxes!`);
2. recompute $Q_{\mathrm{em}}$, $Q_{\mathrm{abs}}$ (`compute_radiation_source_terms!`);
3. assemble $R$;
4. test convergence on the residual (below);
5. assemble the tridiagonal Jacobian;
6. solve $J\,\Delta G=-R$ via Thomas;
7. update $G_i \mathrel{+}= \text{damping}\cdot\Delta G_i$ and clip $G_i \leftarrow \max(G_i,0)$;
8. test convergence on the update.

**Defaults (the code is authoritative):** `tol = 1e-6`, `maxiter = 20`,
`damping = 1.0`, `initialize = true` (radiation_fluxes.jl:208–218). (A docstring
elsewhere says `tol=1e-10`; the executed signature uses `1e-6`.)

**Relative convergence.** $G$ is large ($10^4$–$10^6\ \mathrm{W\,m^{-2}}$), so an
absolute tolerance would be unphysical. The residual test scales by
$G_{\text{scale}}=\max(\|G\|_\infty, 1)$,

$$
\|R\|_\infty < \mathrm{tol}\cdot G_{\text{scale}}\ \Rightarrow\ \text{converged},
$$

and a second test on the relative update,
$\|\Delta G\|_\infty / \max(\|G\|_\infty,1) < \mathrm{tol}$, catches a
machine-zero plateau when $G$ is small. A warning is emitted only if the final
residual exceeds $100\,\mathrm{tol}\,G_{\text{scale}}$ (significant
non-convergence), avoiding noise on benign near-convergence.

**Non-negativity.** $G\ge0$ is enforced by `G[i] = max(G[i], 0.0)` after each
update, since the incident radiation cannot be negative.

### 8.7.7 Initialization

When $G$ is all-zero on entry (cold start) the field is initialized from local
thermodynamic equilibrium (`initialize_g_from_temperature!`,
radiation_fluxes.jl:435–448):

$$
G_i = 4\sigma T_i^4 ,
$$

which makes $S_{\mathrm{rad}}=\alpha(G-4\sigma T^4)$ vanish initially — a good first
guess that avoids cold-start oscillation. (This same blackbody equilibrium is
exposed by `initialize_radiation_intensity!`,
src/Properties/radiation_cache.jl:191–224.)

---

## 8.8 Energy-Equation Coupling for Each Model

The condensed-phase energy equation (Formulation A, §3 Governing Equations and §7)
carries a volumetric radiation source $Q_{\text{rad}}$ on its right-hand side:

$$
\rho c_p^{\text{eff}}\,\frac{\partial T}{\partial t}
= -\nabla\!\cdot\!\mathbf{q}_{\text{cond}} + S_{\text{conv}} + S_{\text{gen}}
+ Q_{\text{rad}} + Q_{\text{rxn}} .
$$

The four models populate $Q_{\text{rad}}$ — or bypass it for the surface route — as
follows.

- **`NO_RADIATION`.** $Q_{\text{rad}}=0$. Surface re-radiation/incident handled only
  by the boundary balance (a `RadiativeBC`/`RadiativeFluxBC`).

- **`SURFACE_ABSORPTION`.** $Q_{\text{rad}}=0$ in the dispatched RHS
  (`compute_radiation_source!` zeroes it). The absorbed incident flux enters through
  the surface-temperature Newton solve as the boundary flux $q_{\text{BC}}=I_{\text{ext}}
  \alpha$ (boundary_temperature.jl:414–423); see §8.2 and §12.

- **`BEER_LAMBERT`.** $Q_{\text{rad},i}=q^{\text{rad}}_i$ from §8.3, added to the
  RHS in step 5 of `compute_rhs!` (src/Discretization/finite_volume.jl:339–348) and
  divided by $\rho c_p^{\text{eff}}$ in the final assembly (finite_volume.jl:384–404):

  $$
  \Big(\frac{\mathrm{d}T_i}{\mathrm{d}t}\Big)_{\text{rad}}
  = \frac{q^{\text{rad}}_i}{\rho c_p^{\text{eff}}_i}.
  $$

  The `RadiativeFluxBC` boundary flux is set to **zero** in Beer–Lambert mode
  (boundary_temperature.jl:426–428) so that the incident energy is not double-counted
  (it is already deposited volumetrically). The structured-operator wrapper
  `RadiationOperator` (src/Jacobian/operators/radiation.jl) routes the same physics
  through `compute_beer_lambert_source!`, treating `dT` as the `q_rad` buffer before
  the $\rho c_p$ divide.

- **`P1_QUASI_STEADY`.** Handled on the separate P1 residual
  (`compute_p1_residual!`, src/Experimental/RadiationP1.jl): the unified
  `residual!` runs first with `radiation_model = P1_QUASI_STEADY` from
  `def.options` (the `P1Residual` constructor enforces this), which zeroes the
  core volumetric source, the boundary `RadiativeFluxBC` **and** the boundary
  `RadiativeBC`; then the quasi-steady $G$ is solved and
  $S_{\mathrm{rad},i}=\alpha_i(G_i-4\sigma T_i^4)$ is added as
  $\mathrm{d}T_i/\mathrm{d}t\mathrel{+}=S_{\mathrm{rad},i}/\rho c_p^{\text{eff}}_i$.
  All radiation energy is owned by P1: the external flux enters through the
  Marshak BC (§8.6.1), emission and ambient exchange through the volumetric
  $4\sigma\alpha T^4$ term plus the Marshak exchange.

**The radiation-model BC dispatch is the linchpin.** The model-dependent
behavior of `RadiativeFluxBC` (boundary_temperature.jl, boundary_fluxes.jl) is
what keeps the incident energy single-counted across all four models:
surface-deposited for `NO_RADIATION`/`SURFACE_ABSORPTION`, volumetric (and
hence boundary-zero) for `BEER_LAMBERT`/`P1_QUASI_STEADY`. `RadiativeBC` is
additionally zeroed under `P1_QUASI_STEADY` only — under P1 the medium emits
volumetrically and exchanges with the ambient through the Marshak BC, so a
surface re-radiation law on top would double-count both channels; under
`BEER_LAMBERT` (which deposits only the *incident* flux volumetrically) the
surface emission stays at the boundary and `RadiativeBC` remains active. The
same dispatch threads through `CombinedThermalBC`, so a cone-calorimeter BC
such as
`RadiativeFluxBC(50e3) + RadiativeBC(ε=0.9,T_inf=300) + ConvectiveBC(h=10,T_inf=300)`
behaves correctly under any radiation model — under P1 the two radiative
members contribute exactly zero and only the convective member survives at
the boundary.

---

## 8.9 Regime Diagnostics

Three diagnostics quantify whether a chosen model is appropriate.

**Optical thickness.** `optical_thickness(mesh, prop_cache)`
(radiation_p1.jl) returns $\tau=\sum_i\alpha_i\Delta z_i=\int_0^L\alpha\,
\mathrm{d}z$, and `optical_thickness_cell` gives the per-cell
$\tau_i=\alpha_i\Delta z_i$ (a resolution diagnostic, not a validity
criterion). The standard validity bands (radiation_types.jl:53–57) are:

| $\tau=\alpha L$ | regime | recommendation |
|---|---|---|
| $\tau>10$ | optically thick | P1 excellent; surface absorption acceptable |
| $1<\tau<10$ | intermediate | P1 good; Beer–Lambert recommended |
| $0.1<\tau<1$ | optically thin | Beer–Lambert recommended |
| $\tau<0.1$ | transparent | `NO_RADIATION` |

**P1 validity.** `check_p1_validity(mesh, prop_cache)` (radiation_p1.jl)
classifies the regime by the **domain** optical thickness and warns when P1 is
marginal (15–30 % error for $0.1<\tau<1$) or invalid ($\tau<0.1$, $>30$ % error —
consider discrete ordinates). Validity is a property of the domain $\tau$
only: a per-cell criterion ($\tau_i \ge 1$ in every cell) would be
mesh-dependent — refining the mesh shrinks every $\tau_i$ with no change in the
physics — so no per-cell validity test is applied.

**Radiation energy balance.** `radiation_energy_balance`
(radiation_p1.jl) integrates $\int Q_{\mathrm{em}}\,\mathrm{d}V$ and
$\int Q_{\mathrm{abs}}\,\mathrm{d}V$, reports the net source and the top/bottom
boundary fluxes, and forms the integral residual of the quasi-steady G
equation, $(q_{\text{top}}-q_{\text{bottom}}) + (\text{net source})$, which
vanishes at a converged quasi-steady state (the boundary-flux difference and
the net source *cancel* — they do not equal each other: `q_rad` stores the
physical $+z$ flux at every face, and $\int(Q_{\mathrm{abs}}-Q_{\mathrm{em}})
\,\mathrm{d}V + \oint \mathbf{q}\cdot\hat{\mathbf{n}}\,\mathrm{d}A = 0$). A
speed-of-light radiation energy density $E_{\text{rad}}=G/c$
($c=2.998\times10^8\ \mathrm{m\,s^{-1}}$) is also available as a diagnostic
(radiation_cache.jl:241–244).

---

## 8.10 The Beer–Lambert Non-Local Jacobian

### 8.10.1 Why it is upper-triangular

Beer–Lambert makes $q^{\text{rad}}_i$ depend on the absorption coefficient of every
cell **above** it, through the cumulative attenuation $I_{\text{in}}(i)=
I_{\text{surface}}\exp(-\sum_{j>i}\alpha_j\Delta z_j)$. Differentiating with respect
to a state variable $x_j$ (either $T_j$ or $\xi_{k,j}$) gives a structurally
**upper-triangular** coupling (src/Jacobian/operators/radiation.jl:22–39):

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_j} =
\begin{cases}
\big(I_{\text{in}}(i) - q^{\text{rad}}_i\,\Delta z_i\big)\,\dfrac{\partial \alpha_i}{\partial x_i}, & j=i,\\[2.0ex]
-\,q^{\text{rad}}_i\,\Delta z_j\,\dfrac{\partial \alpha_j}{\partial x_j}, & j>i,\\[1.0ex]
0, & j<i.
\end{cases}
$$

The **diagonal** ($j=i$) term has a clean identity. From
$q^{\text{rad}}_i=I_{\text{in}}(1-e^{-\alpha_i\Delta z_i})/\Delta z_i$,

$$
\frac{\partial q^{\text{rad}}_i}{\partial\alpha_i}
= I_{\text{in}}\,e^{-\alpha_i\Delta z_i}
= I_{\text{in}} - q^{\text{rad}}_i\,\Delta z_i ,
$$

so the diagonal coefficient is exactly $\big(I_{\text{in}}-q^{\text{rad}}_i\Delta z_i\big)\partial\alpha_i/\partial x_i$.
(An alternative form $I_{\text{in}}(1-(1-e^{-\tau})/\tau)$ appears in some notes; it
is *not* equal to this and is incorrect — the code form above is the verified
derivative.)

The **off-diagonal** ($j>i$) term comes from
$\partial I_{\text{in}}(i)/\partial\alpha_j = -I_{\text{in}}(i)\,\Delta z_j$ for
$j>i$, giving $\partial q^{\text{rad}}_i/\partial x_j = -q^{\text{rad}}_i\,\Delta z_j\,
\partial\alpha_j/\partial x_j$.

**The Kirchhoff surface-transmissivity chain.** On the default Kirchhoff path
(§8.3.2) the injected intensity is
$\varepsilon_{\text{eff}}(T_n,\xi_{\cdot,n})\,I_{\text{surface}}$, so every
$q^{\text{rad}}_i \propto \varepsilon_{\text{eff}}$ and each energy row gains a
dependence on the **surface cell's** state:

$$
\frac{\partial q^{\text{rad}}_i}{\partial x_n} \mathrel{+}=
q^{\text{rad}}_i\,\frac{1}{\varepsilon_{\text{eff}}}
\frac{\partial\varepsilon_{\text{eff}}}{\partial x_n}\qquad(i<n).
$$

In the structured factorization this folds into the same suffix coupling via
$v_n \mathrel{-}= (\partial\varepsilon_{\text{eff}}/\partial x_n)/\varepsilon_{\text{eff}}$;
the surface cell's own diagonal closure recomputes
$\varepsilon_{\text{eff}}(T,\xi)$ live so its
$\partial\varepsilon_{\text{eff}}/\partial(T_n,\xi_{\cdot,n})$ enters exactly;
and the sparsity pattern gains $T$-row($i$) × $T$-column($n$) entries for $i<n$
(the $\xi_{\cdot,n}$ columns already exist in the suffix block;
`_emit_radiation_pattern!`, src/Jacobian/sparsity.jl). With an explicit BC
absorptivity (`radiation_kirchhoff = false`) the scaling is folded into
$I_{\text{surface}}$ upstream and none of these terms appear. The chain is
verified entry-by-entry against the dense AD reference (§18.3).

### 8.10.2 Structured factorization (SuffixSumCoupling)

The $j>i$ block factors as $\mathrm{Diag}(-q^{\text{rad}})\cdot U\cdot
\mathrm{diag}(v)$ where $v_j=\Delta z_j\,\partial\alpha_j/\partial x_j -
\delta_{j,n}(\partial\varepsilon_{\text{eff}}/\partial x_n)/\varepsilon_{\text{eff}}$
(the last term on the Kirchhoff path only) and $U$ is
the strict upper-triangular all-ones matrix — i.e. a **suffix scan**
$\sum_{k>i}v_k$. This is the `SuffixSumCoupling(v)` primitive
(`build_suffix_sum_coupling_alpha`, radiation.jl:256–283; for the energy row the
per-row scale is $-q^{\text{rad}}_i/\rho c_p^{\text{eff}}_i$, embedded as
`EmbeddedScaledSuffixSumCoupling`, radiation.jl:416–470). It applies in $O(n)$
without ever materializing the dense $O(n^2)$ fill. The diagonal term — including
the $\rho c_p$ chain rule $-q^{\text{rad}}/(\rho c_p)^2\,\partial\rho c_p/\partial x_i$
— is captured by a cell-local ForwardDiff over the closure
$(T,\xi)\mapsto I_{\text{in}}(1-e^{-\alpha\Delta z})/(\Delta z\,\rho c_p^{\text{eff}})$
with $I_{\text{in}}$ held fixed (radiation.jl:376–407). For the standard linear
mixing $\alpha=\sum_j\xi_j\alpha_j$ the derivatives reduce to $\partial\alpha/\partial
T=0$ and $\partial\alpha/\partial\xi_j=\alpha_j$, so the suffix coupling is cheap;
the AD path keeps it generic for future $T$-dependent absorption laws. (See §15,
Temporal Integration & the Structured Jacobian, for how `SuffixSumCoupling` is
applied and solved.)

> **Solver impact.** `has_nonlocal_jacobian_coupling(BEER_LAMBERT) == true`, so the
> sparsity pattern is upper-triangular rather than tridiagonal. The other three
> models keep the tridiagonal pattern. Under the `:ale_reduced` assembly profile the
> embedded suffix couplings are dropped (radiation.jl:410–416), trading a small
> Jacobian inexactness for speed.

---

## 8.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to ours on first use. All three reference
solvers are condensed-phase pyrolysis models; their radiation treatments span the
same idioms Pyrolysis.jl supports, plus a two-flux variant in FDS.

### 8.11.1 ThermaKin2Ds

ThermaKin (Stoliarov & Lyon) uses **in-depth exponential (Beer–Lambert)
absorption** of the external flux together with explicit surface **re-radiation**.
In ThermaKin notation the external radiant flux is $\dot q''_{ex}$ (our
$I_{\text{ext}}$), the absorption coefficient is per-component and forms a
volumetric coefficient (our $\alpha_{\text{eff}}=\sum_j\xi_j\alpha_j$), and the
in-depth source is $\dot q''_{ex}\,\alpha\,e^{-\int\alpha\,\mathrm{d}x}$ — identical
in form to our $q'''=\alpha I$ with $I=I_{\text{surface}}e^{-\int\alpha\,\mathrm{d}z}$
(our `BEER_LAMBERT`). ThermaKin additionally models re-radiation $I_{rr}=
\varepsilon\sigma T_s^4$ at the surface (often with reflection of a fraction of
$\dot q''_{ex}$), which Pyrolysis.jl supplies via a top `RadiativeBC` re-radiation
term $F\varepsilon_{\text{eff}}\sigma(T_\infty^4-T_s^4)$ in the boundary balance
(§12). ThermaKin also offers a radiative-conductivity augmentation $k\propto T^3$
to mimic radiative diffusion in char — the physical content of our P1 diffusion
term $\nabla\!\cdot\!(D_{\mathrm{rad}}\nabla G)$, which Pyrolysis.jl models
explicitly through `P1_QUASI_STEADY` rather than by inflating $k$. ThermaKin
publishes heats of reaction with the opposite sign convention ($h>0$ exothermic);
this does not touch the radiation terms (overload note H1, §2, Nomenclature).

### 8.11.2 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) deposits the external radiant flux
**in-depth** with an exponential absorption profile and treats the surface energy
balance with combined radiative–convective exchange. In Gpyro notation the
pre-exponential of any property is $Z$ (our $A_i$), the volume fraction is $X_\alpha$
(our $v_j$), and the in-depth absorbed flux follows $\dot q''_e\,\kappa\,e^{-\kappa
x}$ — again our `BEER_LAMBERT` form with Gpyro's $\kappa$ (the absorption
coefficient, distinct from the permeability $K$ in §9) mapping to our volumetric
$\alpha_{\text{eff}}$. Gpyro's effective radiative conductivity in optically thick
char (a temperature-cubed conductivity contribution) plays the same role as our P1
diffusion term; Pyrolysis.jl keeps the two mechanisms separate (a genuine $k(T,\xi)$
mixture conductivity in §5/§7 and an explicit $G$ field in P1). Gpyro's surface
emissivity enters the boundary balance just as our $\varepsilon_{\text{eff}}=
\sum_j v_j^{\text{solid}}\varepsilon_j$ does.

### 8.11.3 FDS

The FDS solid-phase model uses a **two-flux (Schuster–Schwarzschild)** in-depth
radiation model for the condensed phase: forward and backward hemispheric
intensities $I^+,I^-$ are propagated with absorption coefficient $\kappa_s$ (our
$\alpha$), and the radiative source in the 1D Fourier conduction equation is the
divergence of the net radiative flux $I^+-I^-$. Pyrolysis.jl's `BEER_LAMBERT`
corresponds to the single-flux limit (only the inward-going beam, no internal
re-emission), while `P1_QUASI_STEADY` is the diffusion-limit analog of the two-flux
model (P1 and two-flux/Schuster–Schwarzschild are the two classical first-order
closures of the RTE; both reduce to a single second-order ODE for a combined field).
FDS applies surface absorptivity/emissivity through Kirchhoff's law at the gas–solid
interface — the same default ($\alpha=\varepsilon_{\text{eff}}$) Pyrolysis.jl uses
for `RadiativeFluxBC`. FDS's gray-band absorption maps to our gray (wavelength-
independent) $\alpha_{\text{eff}}$. The principal difference is dimensionality and
the explicit two-stream bookkeeping: FDS resolves $I^\pm$ separately, whereas
Pyrolysis.jl resolves either the inward beam ($I$, Beer–Lambert) or the isotropic
incident field ($G$, P1).

### 8.11.4 Summary table

| Aspect | Pyrolysis.jl | ThermaKin2Ds | Gpyro | FDS (solid) |
|---|---|---|---|---|
| In-depth absorption | Beer–Lambert $q'''=\alpha I$ | Beer–Lambert | exp. in-depth | two-flux $I^\pm$ |
| Re-radiation diffusion | explicit P1 ($G$) | radiative $k\propto T^3$ | radiative $k$ | two-flux divergence |
| Surface absorptivity | $\varepsilon_{\text{eff}}$ (Kirchhoff) | $\varepsilon$ + reflection | $\varepsilon$ | Kirchhoff |
| Spectral model | gray | gray | gray | gray-band |
| Surface re-radiation | `RadiativeBC` $\varepsilon\sigma(T_\infty^4-T_s^4)$ | $\varepsilon\sigma T_s^4$ | radiative–convective | $\varepsilon\sigma T_s^4$ |

---

## 8.12 Experimental Status of P1: the Separate Residual

P1 is **Experimental**: it is research-grade and out of scope for the unified
trait-dispatched residual that drives `SURFACE_ABSORPTION`/`BEER_LAMBERT` runs
(`RadiationP1.jl` header). Three structural consequences follow.

**No dedicated state layout.** $G$ is an algebraic (quasi-steady) field, not a
state variable: it lives in the `RadiationCache` and is re-solved from the
current $(T,\xi,z)$ state at **every** residual evaluation. The ODE state
vector is the ordinary unified `StateLayout` (T, ξ, z when ALE, χ when
WithChi); build `def`, `ws`, and `u0` with
`initialize_state(problem; radiation_model = P1_QUASI_STEADY)` (which applies
the problem's initial conditions before packing — a bare
`create_state_vector` call packs the mesh as-is, §15.2) and read
the field with `p.rad_cache.G` or `get_radiation_intensity(p::P1Residual, i)`
— always the field of the most recent residual evaluation, never a frozen
copy.

**A thin dedicated residual.** `compute_p1_residual!`
(src/Experimental/RadiationP1.jl) delegates to the unified `residual!` —
conduction, reactions, transport, ALE, χ all come from the production code
path — and then (1) mirrors `u` into the legacy `Mesh1D` for the
Float64-only P1 kernels, (2) updates $\alpha$/$D_{\mathrm{rad}}$ and solves
the quasi-steady $G$ system, and (3) adds
$S_{\mathrm{rad}}/\rho c_p = \alpha(G-4\sigma T^4)/\rho c_p$ to the energy
rows. The `P1Residual` constructor **requires**
`def.options.radiation_model == P1_QUASI_STEADY` (build the def with
`to_problem_def(problem; radiation_model = P1_QUASI_STEADY)`): the unified
residual then runs with zero core volumetric radiation source and zero
boundary `RadiativeFluxBC`/`RadiativeBC` flux, so the incident energy and the
surface emission are counted exactly once — through the Marshak BC and the
volumetric $4\sigma\alpha T^4$ emission (§8.8).

**Float64-only.** The `RadiationCache` arrays are `Float64`, and the P1 Newton
solve runs in Float64; AD `Dual` parts of the time variable are stripped via
`ForwardDiff.value` when storing boundary fluxes (radiation_p1.jl), and
state-vector `Dual`s are not supported. Consequently P1 sensitivities through
the $G$ field are not propagated by forward-mode AD the way the core models'
Jacobians are; P1 stays here until the trait system gains a first-class
`P1Radiation` mode (§15, §17 Sensitivity Analysis).

---

## 8.13 Limitations and Open Issues

The radiation subsystem makes several deliberate simplifications; the following are
the open issues most likely to matter in practice.

1. **Scattering is not exercised.** $\sigma_s$ is plumbed through
   $D_{\mathrm{rad}}=1/(3(\alpha+\sigma_s))$ but is taken $=0$ everywhere for
   pyrolysis. Strongly scattering chars or foams would need $\sigma_s>0$ and, for
   anisotropy, a higher-order closure (P3 or discrete ordinates).

2. **Gray (wavelength-independent) properties.** $\alpha_j$ and $\varepsilon_j$ are
   gray; real polymers and chars are spectrally selective. A gray-to-band extension
   would require per-band absorption coefficients and a banded $G$ field.

3. **P1 is quasi-steady only.** $\partial G/\partial t=0$ is built in: $G$ is
   an algebraic field in the `RadiationCache`, re-solved at every residual
   evaluation, and is not part of the ODE state at all. A fully transient P1
   would integrate $G$ as a genuine state variable with its own time scale —
   rarely needed given $\tau_{\mathrm{rad}}\ll\tau_{\mathrm{cond}}$, but
   unavailable today.

4. **The surface energy-balance Newton Jacobian neglects $\mathrm{d}\varepsilon/\mathrm{d}T$**
   (§8.8, §12.2). The Marshak coefficient itself is consistent — the half-cell
   effective coefficient $c_{\text{eff}} = 1/(1/c_M + \delta/D_{\mathrm{rad}})$,
   $c_M = \varepsilon/(2(2-\varepsilon))$, appears identically in the residual
   and in the Jacobian (§8.7.4). The remaining approximation is that,
   where a state-dependent emissivity is used, the temperature dependence of
   $\varepsilon$ is dropped from the surface-temperature derivative. For the
   near-constant emissivities used in practice this is negligible; a strongly
   $T$-dependent $\varepsilon$ can slow Newton convergence.

5. **No interior emissivity factor by design** (§8.5). This is physically correct
   for a gray body in LTE, but it means surface emissivity influences the interior
   field *only* through the Marshak boundary; users expecting a volumetric
   $\varepsilon\sigma T^4$ emission will not find one.

6. **P1 is AD-opaque (Float64-only).** Sensitivities through the $G$ field are not
   propagated by the forward-mode AD path used elsewhere; P1 remains Experimental and
   outside the unified trait dispatch (§8.12, §15, §17).

7. **`SURFACE_ABSORPTION` assumes an opaque first cell** (§8.2). Absorbed
   incident radiation enters through the *boundary balance*; the model is
   accurate only when the first cell is optically thick ($\tau_N\gg1$). If the
   first cell is coarse or only marginally opaque, refine the surface cell or
   switch to `BEER_LAMBERT`.

8. **Transparent-limit clamp on $D_{\mathrm{rad}}$.** Below
   $10^{-10}\ \mathrm{m^{-1}}$ extinction,
   `compute_radiation_diffusion_coefficient` — the single canonical guard, used
   by the one `update_radiation_properties!` method — clamps the diffusion
   coefficient to the large finite value $1/(3\cdot 10^{-10})$ rather than
   letting it diverge (or zeroing it, which would make the quasi-steady Newton
   system singular wherever $\alpha = 0$). A transparent cell thus behaves as a
   near-perfect conductor of $G$. P1 is not recommended for $\tau<0.1$
   regardless.
