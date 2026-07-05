# 2. Nomenclature

This is the **single, authoritative symbol table** for the Pyrolysis.jl Technical
Reference Guide and User Guide. Every writer **MUST** conform to these symbols,
meanings, and SI units verbatim. Where a symbol is used differently in different
source subsystems, the **conflict has been resolved here** and the resolution is
binding. Conflicts and overloads are called out explicitly in the
"Overloaded / reserved symbols" section at the end.

Conventions:
- The spatial coordinate is **z** (axial, through-thickness). `z = 0` is the
  **bottom / substrate** (boundary id 2, tag `:bottom`); `z = L` is the **top /
  exposed surface** (boundary id 1, tag `:top`). Heat enters from the top; the
  material shrinks downward. The transverse/lateral plane carries the
  cross-sectional area `A`.
- Sign convention for fluxes: **positive flux = transport in the +z direction**
  (bottom → top). Divergence `(F_R − F_L)·A/V` is positive when the quantity
  flows *out* of a cell.
- Heat-of-reaction sign convention: **h > 0 endothermic (cools)**,
  **h < 0 exothermic (heats)**; the volumetric heat source is `Q_rxn = −Σ h_r r_r`,
  so endothermic reactions give `Q_rxn < 0`. (This is the storage convention used
  internally; note ThermaKin/Gpyro publish `h > 0 = exothermic` — see overload
  note H1.)
- "Bulk" density `ρ_j` is the pure-phase density used by mixing rules and mixture
  density; "intrinsic/skeletal" density `ρ_{i,j}` (cell-wall density) is used
  ONLY in the porosity formula.
- Mass concentration `ξ_j` (kg/m³) is the **primary species state variable**, not
  mass fraction. `ξ_j = Y_j ρ`.

---

## 1. Latin letters

| Symbol | Meaning | SI units |
|---|---|---|
| `A` | Cross-sectional area of the 1D column (may vary in time via lateral-shrinkage law); also face area in the divergence operator | m² |
| `A_0` | Initial / calibration-reference cross-sectional area (fixed at setup) | m² |
| `A_{αβ}`, `A_i` | Arrhenius pre-exponential factor for reaction *i* (or *αβ*) | s⁻¹ (1st-order); m³·kg⁻¹·s⁻¹ (bimolecular) |
| `B` | Spalding mass-transfer number (liquid evaporation, comparison only) | – |
| `c_p`, `c_{p,j}` | Specific heat capacity (of component *j*) | J·kg⁻¹·K⁻¹ |
| `c` | Speed of light (radiation energy density diagnostic only) | m·s⁻¹ |
| `D_{AB}`, `D_j` | Binary (Chapman–Enskog) diffusion coefficient | m²·s⁻¹ |
| `D_rad` | P1 radiation diffusion coefficient, `1/(3(α+σ_s))` | m |
| `d`, `d_L`, `d_R`, `d_{LR}` | Distance from cell center to face / between adjacent cell centers | m |
| `E`, `E_{αβ}`, `E_i` | Activation energy for reaction *i* | J·mol⁻¹ |
| `E_matrix` | Matrix-only sensible energy of the domain (gas excluded) | J |
| `E_total` | Total sensible energy of the domain (matrix + gas) | J |
| `e` | Euler's number / base of natural log (literal) | – |
| `F` | Radiative view factor | – |
| `F_L`, `F_R` | Generic face flux (left/right) in the divergence operator | (flux units)·m⁻² |
| `f` | Generic field, or Newton-residual function `f(T_s)` | varies |
| `G` | Incident radiation intensity (P1 primary variable) | W·m⁻² |
| `G_ext` | Ambient radiation field at boundary, `4σT_∞⁴` | W·m⁻² |
| `g` | Gravitational acceleration (Darcy buoyancy, comparison only) | m·s⁻² |
| `H(·)` | Heaviside step (implemented as smooth tanh ramp) | – |
| `h` (kinetics) | Heat of reaction, per kg of first reactant | J·kg⁻¹ |
| `h_conv` | Convective heat-transfer coefficient | W·m⁻²·K⁻¹ |
| `h_m` | External mass-transfer (film) coefficient | m·s⁻¹ |
| `h_P` | Pressure-transfer coefficient (convective pressure BC) | m·Pa⁻¹·s⁻¹ |
| `Δh`, `Δh_g` | Sensible enthalpy difference (midpoint-rule integral) | J·kg⁻¹ |
| `I(z)`, `I_in` | Radiation intensity along path / at cell entrance (Beer–Lambert) | W·m⁻² |
| `I_surface`, `I_ext` | Incident radiation intensity at domain surface / external source | W·m⁻² |
| `J_j` | Mass flux of (gas) component *j* at a face | kg·m⁻²·s⁻¹ |
| `J_j^diff`, `J_j^adv` | Diffusive / advective part of the gas flux | kg·m⁻²·s⁻¹ |
| `K`, `κ` | Permeability — **use `κ`**; `K` reserved for Gpyro-comparison text only | m² |
| `k`, `k_j` | Thermal conductivity (of component *j*) | W·m⁻¹·K⁻¹ |
| `k_eff`, `k_f` | Effective (mixture) / face-averaged thermal conductivity | W·m⁻¹·K⁻¹ |
| `k_B` | Boltzmann constant (reduced temperature `T* = k_B T/ε`) | J·K⁻¹ |
| `L` | Domain thickness (total material length); also characteristic length in convection correlations | m |
| `M`, `M_j` | Molar mass (of component *j*) | kg·mol⁻¹ |
| `M_ref` | Reference molar mass in the Chapman–Enskog combining rules (= 0.029, air; §5.10) | kg·mol⁻¹ |
| `M_dry` | Total initial dry-solid mass in the column | kg |
| `m_j` | Total mass of component *j* in domain | kg |
| `m_{dry,i}` | Initial dry-solid mass in cell *i* (χ denominator) | kg |
| `ṁ_g` | Outward gas mass flux at surface (transpiration) | kg·m⁻²·s⁻¹ |
| `N` | Partial-pressure sum (numerator of ideal-gas pressure) | Pa |
| `N_C` (`NC`) | Number of chemical components (type parameter) | – |
| `N_R` (`NR`) | Number of reactions (type parameter) | – |
| `n` | Number of finite-volume cells | – |
| `n_nodes`, `n_faces` | Number of mesh nodes / faces (= n+1 in 1D) | – |
| `n_{i,j}`, `n_s` | Reaction order of reaction *i* w.r.t. reactant *j* | – |
| `Nu` | Nusselt number (convection BC) | – |
| `P` | Gas pressure in the porous medium | Pa |
| `P_ref` | Reference pressure (= 101325) | Pa |
| `Q_rxn` | Volumetric reaction heat source, `−Σ h_r r_r` | W·m⁻³ |
| `Q_rad` | Volumetric radiation absorption source | W·m⁻³ |
| `Q_em`, `Q_abs` | P1 volumetric emission `4σαT⁴` / absorption `αG` | W·m⁻³ |
| `q`, `q_cond` | (Conductive) heat flux | W·m⁻² |
| `q_BC` | Boundary-condition heat flux (positive into domain) | W·m⁻² |
| `q'''`, `q_rad'''` | Volumetric radiation source (Beer–Lambert) | W·m⁻³ |
| `R_g` (`R`, `R_gas`) | Universal gas constant (= 8.314462618) | J·mol⁻¹·K⁻¹ |
| `r`, `r_i` | Reaction rate (per unit volume) of reaction *i* | kg·m⁻³·s⁻¹ |
| `r` (limiters) | Slope ratio argument of a flux limiter `ψ(r)` | – |
| `S` | Sutherland constant (= 110.4 K, air) in the viscosity law | K |
| `S_j`, `S_{j,i}` | Species source term for component *j* (in cell *i*) | kg·m⁻³·s⁻¹ |
| `S_conv` | Volumetric gas-advection (convection) energy source | W·m⁻³ |
| `S_gen` | Volumetric gas-generation enthalpy sink (optional) | W·m⁻³ |
| `S_rad` | Net volumetric P1 radiation source, `α(G − 4σT⁴)` | W·m⁻³ |
| `S_dry,i` | Dry-pyrolysis-gas source term in cell *i* | kg·m⁻³·s⁻¹ |
| `T` | Absolute temperature | K |
| `T_s`, `T_w` | Surface (boundary face) temperature | K |
| `T_cell` | Temperature at the cell center adjacent to a boundary | K |
| `T_∞` | Ambient / surroundings temperature | K |
| `T_ref` | Reference temperature for enthalpy datum (= 298.15) | K |
| `T_face`, `T_f` | Arithmetic-mean temperature at a face | K |
| `T_min,i`, `T_max,i` | Lower / upper temperature gate for reaction *i* | K |
| `T*` | Reduced temperature `k_B T / ε_AB` (collision integral) | – |
| `t` | Time | s |
| `u` | Flat ODE/DAE state vector | mixed |
| `du` | Time derivative of state vector | mixed |
| `V`, `V_i` | Cell volume (of cell *i*) | m³ |
| `v_g` | Bulk (Darcy) gas velocity at a face | m·s⁻¹ |
| `v_j` | Volume fraction of component *j* | – |
| `v_j^solid` | Solid volume fraction (excludes non-swelling gases) | – |
| `W` | Implicit-RK "W" matrix, `I/(γΔt) − J` | – |
| `w`, `w_i` | ALE mesh node velocity (at node *i*) | m·s⁻¹ |
| `w_cell` | Cell-centered ALE mesh velocity | m·s⁻¹ |
| `X_α` | Volume fraction (FDS/Gpyro comparison notation; internally `v_j`) | – |
| `x`, `y` | Transverse/lateral coordinates (multi-D comparison only) | m |
| `Y_j` | Mass fraction of component *j*, `ξ_j/ρ` | – |
| `Z` | Pre-exponential factor (Gpyro/FDS comparison notation; internally `A`) | varies |
| `z` | Axial / through-thickness coordinate (primary spatial variable) | m |
| `z_i` | Position of mesh node *i* (ALE state) | m |
| `z_i^c` | Cell-center position of cell *i* | m |
| `Δz`, `Δz_i` | Cell thickness (of cell *i*) | m |

---

## 2. Greek letters

| Symbol | Meaning | SI units |
|---|---|---|
| `α` (`α_eff`) | Volumetric absorption (extinction) coefficient | m⁻¹ |
| `α_j` | Mass-basis absorption coefficient of component *j* | m²·kg⁻¹ |
| `α_conv` | Conversion / extent variable (Gpyro kinetics comparison) | – |
| `β` | Weighting parameter of the WEIGHTED mixing rule, `k = βk_∥ + (1−β)k_series` | – |
| `γ` (`γ_j`) | Swelling factor of component *j* (1 = solids, 0 = non-swelling gas) | – |
| `γ_RK` | Implicit-RK stage coefficient (in `W = I/(γΔt) − J`) | – |
| `Δ` | Forward-difference / change operator (prefix) | varies |
| `δ_{ij}` | Kronecker delta | – |
| `δT`, `δξ_j` | Finite-difference perturbation step | K; kg·m⁻³ |
| `ε`, `ε_j` | Emissivity (of component *j*) | – |
| `ε_eff` | Effective (composition-weighted) emissivity | – |
| `ε_AB` | Lennard-Jones well depth (collision integral) | J |
| `ζ` | Dummy integration variable along the optical path | m |
| `η_i` | AMR error-indicator value for cell *i* | varies |
| `θ`, `θ_i` | Volumetric strain rate (dilation) of cell *i*, `(1/V) dV/dt` | s⁻¹ |
| `θ_A` | Column-global **radial** volumetric strain rate, `(1/A) dA/dt` | s⁻¹ |
| `θ_L` | **Axial** volumetric strain rate (`θ − θ_A`) | s⁻¹ |
| `θ_i^j` | Stoichiometric coefficient (ThermaKin comparison notation; internally `ν_{i,j}`) | – |
| `κ` (`κ_j`, `κ_face`) | Permeability (Darcy) of component *j* / at a face | m² |
| `κ_s` | Solid absorption coefficient (FDS comparison notation; internally `α`) | m⁻¹ |
| `λ` (`λ_j`, `λ_face`, `λ_eff`) | Gas transfer (diffusion) coefficient | m²·s⁻¹ |
| `ℓ` | Characteristic particle/pore diameter (Carman–Kozeny) | m |
| `μ`, `μ_face` | Gas dynamic (Sutherland) viscosity | Pa·s |
| `μ_ref` | Reference viscosity (= 1.716×10⁻⁵ at 273.15 K, air) | Pa·s |
| `ν` | Kinematic viscosity (Gpyro Darcy comparison) | m²·s⁻¹ |
| `ν_{i,j}` | Stoichiometric coefficient (yield) of component *j* in reaction *i* (negative reactant, positive product) | – |
| `ξ_j`, `ξ_{j,i}` | Mass concentration of component *j* (in cell *i*) — **primary species state** | kg·m⁻³ |
| `tanh(ξ/ξ_th)` | Depletion rate factor (per reactant) | – |
| `ρ` | Mixture / bulk density (`Σ ξ_j` for total; mixture-density formula for effective) | kg·m⁻³ |
| `ρ_j` | Bulk (pure-phase) density of component *j* | kg·m⁻³ |
| `ρ_{i,j}` | Intrinsic / skeletal density of component *j* (porosity formula only) | kg·m⁻³ |
| `ρ_g` | Gas-phase density (ideal-gas law) | kg·m⁻³ |
| `ρc_p^eff` | Effective volumetric heat capacity (matrix only) | J·m⁻³·K⁻¹ |
| `ρc_p^total` | Total volumetric heat capacity (matrix + gas) | J·m⁻³·K⁻¹ |
| `σ` | Stefan–Boltzmann constant (= 5.670374419×10⁻⁸) | W·m⁻²·K⁻⁴ |
| `σ_s` | Scattering coefficient (≈ 0 for pyrolysis) | m⁻¹ |
| `σ_AB` | Lennard-Jones collision diameter | Å |
| `σ_mon` | Gaussian width of the interface monitor (r-refinement) | m |
| `τ`, `τ_i` | Optical thickness (total / cell-local), `α Δz` | – |
| `τ_ILU` | ILU drop tolerance | – |
| `φ` | Porosity (gas void fraction) | – |
| `χ`, `χ_i` | **Pyrolysis progress** in cell *i* (cumulative dry-gas mass released per unit initial dry-solid mass) | – |
| `χ̄` | Mass-weighted column-average pyrolysis progress | – |
| `ψ(r)` | Flux-limiter function | – |
| `Ω_D` | Diffusion collision integral (Neufeld) | – |
| `ω` (`ω(z)`) | Monitor function (r-refinement equidistribution) | varies |
| `ω̇'''` | Volumetric reaction rate (Gpyro comparison notation; internally `r`) | kg·m⁻³·s⁻¹ |

---

## 3. Subscripts

| Subscript | Meaning |
|---|---|
| `j` | Component / species index (1…N_C) |
| `i` | Cell index (1…n); also reaction index in kinetics context (see overload I1) |
| `r`, `αβ` | Reaction index |
| `g` | Gas phase / gaseous component |
| `s` | Solid / condensed phase |
| `ℓ`, `liq` | Liquid phase |
| `L`, `R` | Left / right of a face |
| `f` | Face-centered quantity |
| `face` | Face-averaged property |
| `eff` | Effective (mixture) property |
| `cell` | Cell-centered quantity |
| `s` (boundary) | Surface (boundary face) — e.g. `T_s`, `ξ_{j,s}` |
| `∞` | Ambient / surroundings |
| `ref` | Reference value |
| `in` | Cell-entrance (radiation) |
| `ext` | External source |
| `min`, `max` | Lower / upper bound |
| `dry` | Dry-pyrolysis-gas (excludes water vapor) |
| `init`, `0` | Initial value |
| `top`, `bot` | Top (z=L) / bottom (z=0) boundary |
| `rxn` | Reaction-source contribution |
| `cond`, `conv`, `rad`, `gen` | Conduction / convection / radiation / gas-generation contribution |
| `∥`, `series` | Parallel / series mixing limit |
| `matrix`, `total` | Matrix-only / matrix+gas (energy, heat capacity) |

## 4. Superscripts and accents

| Mark | Meaning |
|---|---|
| `°`, `^o`, `_init` | Value at previous time step / initial value |
| `+`, `−` (radiation) | Forward / backward hemispheric intensity (two-flux) |
| `solid` | Restricted to solid/liquid (non-swelling) components |
| `diff`, `adv` | Diffusive / advective flux component |
| `^c` | Cell-center (e.g. `z_i^c`) |
| `‾` (overbar) | Column / spatial average (e.g. `χ̄`, `J̄_g`) or Gpyro-style mixture average (`ρ̄`) |
| `·` (overdot) | Time rate / flux (e.g. `ṁ`, `ω̇`) |
| `'`, `''`, `'''` | Per-length / per-area / per-volume (FDS/Gpyro comparison notation) |

## 5. Operators

| Operator | Meaning |
|---|---|
| `∂/∂t` | Partial time derivative (Eulerian, fixed point) |
| `∂/∂t\|_χ` | ALE material-point (moving-mesh) time derivative |
| `d/dt` | Total time derivative |
| `∇`, `∇·` | Gradient / divergence (1D: `∂/∂z`, `(F_R−F_L)A/V`) |
| `∇²` | Laplacian / second derivative (curvature indicator) |
| `Σ`, `∫` | Sum / integral |
| `clamp(x,a,b)` | Clamp `x` to `[a,b]` |
| `tanh` | Hyperbolic tangent (smooth gates and depletion limiter) |
| `‖·‖_∞` | Infinity norm |
| `diag(v)` | Diagonal matrix from vector `v` |
| `J`, `A`, `C_j` | Jacobian / its sparse block / structured coupling *j* |

---

## 6. Physical constants (fixed values)

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | `R_g` | 8.314462618 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | `σ` | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference pressure | `P_ref` | 101325 | Pa |
| Reference temperature | `T_ref` | 298.15 | K |
| Reference molar mass (air) | `M_ref` | 0.029 | kg·mol⁻¹ |
| Reference viscosity (air, 273.15 K) | `μ_ref` | 1.716×10⁻⁵ | Pa·s |
| Sutherland constant (air) | `S` | 110.4 | K |
| Kozeny constant (spheres) | — | 180 | – |

---

## 7. Overloaded / reserved symbols — resolutions (BINDING)

These symbols are overloaded somewhere in the source subsystems. The resolution
below is authoritative; writers must use the disambiguating subscript/context.

- **χ (chi) — A1 [pyrolysis progress vs. cross-section ratio].** The
  `volume_change`, `problem_residual`, and `geometry` subsystems define `χ` as
  the **pyrolysis progress** (dimensionless, cumulative dry-gas mass released per
  unit initial dry-solid mass). The `jacobian` KB once glosses it as
  "cross-section area ratio." **Resolution:** `χ` is *always* pyrolysis progress;
  the cross-section ratio is `A/A_0 = law(χ̄)` and is never written as `χ`.

- **θ (theta) — A2 [strain rate vs. stoichiometry].** Internally `θ` is the
  **volumetric strain rate** (s⁻¹), with `θ_A` radial and `θ_L` axial. ThermaKin
  uses `θ_i^j` for stoichiometric coefficients. **Resolution:** use `ν_{i,j}` for
  stoichiometry everywhere; reserve `θ` for strain rate. `θ_i^j` appears only in
  the ThermaKin-comparison passage and is immediately mapped to `ν_{i,j}`.

- **A — A3 [area vs. pre-exponential factor].** `A` is the **cross-sectional
  area**; the Arrhenius pre-exponential is `A_i` / `A_{αβ}` (always subscripted by
  reaction). In Jacobian/linear-algebra context `A` also denotes the sparse
  Jacobian block. **Resolution:** bare `A` = area in physics chapters; `A_i` =
  pre-exponential (kinetics); `A` (linear algebra) only inside the Jacobian
  chapter where the meaning is unambiguous and stated. The FDS/Gpyro
  pre-exponential `Z` is mapped to `A_i` on first use.

- **κ vs. K — A4 [permeability].** Use **`κ`** for permeability throughout.
  `K` (Gpyro permeability) appears only in comparison text and is mapped to `κ`.

- **α — A5 [absorption vs. conversion vs. heat-transfer].** `α` (or `α_eff`) =
  **volumetric absorption coefficient** (m⁻¹); `α_j` = mass-basis absorption
  coefficient (m²·kg⁻¹). Conversion (Gpyro) is `α_conv`; FDS `κ_s` maps to `α`.
  The mass-transfer "absorptivity" of a surface is also `α` but always paired
  with a BC context and Kirchhoff `α = ε_eff`.

- **λ — A6 [gas transfer vs. anything else].** `λ` is reserved for the **gas
  transfer (diffusion) coefficient** (m²·s⁻¹). Do not use `λ` for thermal
  conductivity (that is `k`) or eigenvalues.

- **ρ subscripts — A7 [bulk vs. intrinsic].** `ρ_j` = bulk density (mixing,
  mixture density); `ρ_{i,j}` = intrinsic/skeletal density (porosity only). Never
  write a bare `ρ_i` for skeletal density without the explicit "intrinsic"
  qualifier; prefer `ρ_{i,j}`.

- **h — H1 [heat of reaction; sign].** `h` (kinetics) = heat of reaction per kg
  of first reactant, with **h > 0 endothermic** internally; ThermaKin/Gpyro
  publish `h > 0 = exothermic`. The conversion to the volumetric source carries
  the minus sign: `Q_rxn = −Σ h_r r_r`. The convective coefficient is `h_conv`,
  the mass coefficient `h_m`, the pressure coefficient `h_P`, all subscripted.

- **i — I1 [cell index vs. reaction index].** `i` is the **cell index** in
  discretization/geometry/ALE; in the kinetics chapter `i` is the **reaction
  index** (with `j` the component index). Each chapter states its convention in
  its first paragraph. Where both appear together, use `i` = cell, `r` =
  reaction.

- **φ — A8 [porosity only].** `φ` is **porosity**; do not reuse for view factor
  (that is `F`) or any potential. Two historical porosity definitions exist
  (volume-fraction `Σ v_j^gas` vs. intrinsic-density `1 − Σ ξ_j/ρ_{i,j}`); the
  intrinsic-density form is canonical and is what `φ` denotes.

- **S — S1 [source vs. Sutherland constant].** `S` with a subscript
  (`S_conv`, `S_gen`, `S_rad`, `S_j`) is a source term; bare `S` is the
  Sutherland constant in the viscosity law only.
