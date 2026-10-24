# 2. Nomenclature

This is the **single, authoritative symbol table** for the Pyrolysis.jl Technical
Reference Guide and User Guide. Every writer **MUST** conform to these symbols,
meanings, and SI units verbatim. Where a symbol is used differently in different
source subsystems, the **conflict has been resolved here** and the resolution is
binding. Conflicts and overloads are called out explicitly in the
"Overloaded / reserved symbols" section at the end.

Conventions:
- The spatial coordinate is **z** (axial, through-thickness). `z = 0` is the
  **back / bottom / substrate** (boundary tag `:bottom`, face 1, at the pinned
  node 1); `z = L` is the **exposed surface / top** (boundary tag `:top`, face
  `n + 1`). Heat enters from the top; the material shrinks downward. The
  transverse/lateral plane carries the cross-sectional area `A`. Boundary
  conditions are keyed by tag; there are no numeric boundary ids.
- Sign convention for fluxes: a face flux `F_f` is the transport **per unit face
  area along the face normal**, and the normal points from the face's owner cell
  to its neighbour, out of the domain on a boundary face. On a line every
  interior face and the `:top` face therefore have their normal along **+z**
  (positive flux = transport bottom → top), and the `:bottom` face along **−z**
  (positive flux = outflow through the back). The discrete divergence
  `(1/V_i) Σ_f s_if A_f F_f` is positive when the quantity flows *out* of cell
  `i`.
- Heat-of-reaction sign convention: **h > 0 endothermic (cools)**,
  **h < 0 exothermic (heats)**; the volumetric heat source is `Q_rxn = −Σ h_r r_r`,
  so endothermic reactions give `Q_rxn < 0`. (This is the storage convention used
  internally; note ThermaKin publishes `h > 0 = exothermic`; Gpyro and FDS use this convention — see overload
  note H1.)
- "Bulk" density `ρ_j` is the pure-phase density used by mixing rules and mixture
  density; "intrinsic/skeletal" density `ρ_{i,j}` (cell-wall density) is used
  ONLY in the porosity (skeletal volume fraction) formula.
- Mass concentration `ξ_j` (kg/m³) is the **constitutive species variable**, not
  mass fraction: `ξ_j = Y_j ρ`. The integrated discrete state is the cell
  inventory `m_{j,i} = ξ_{j,i} V_i` (kg); every kernel reads
  `ξ_{j,i} = m_{j,i}/Ṽ_i`, with `Ṽ_i` the smoothly floored cell volume.
- The chapters write their mathematics in LaTeX; this table gives each symbol in
  Unicode. Two pairs of glyphs are distinct symbols and must not be interchanged:
  porosity is **`φ`** (U+03C6, LaTeX `\varphi`) and the potential-closure
  multiplier is **`ϕ`** (U+03D5, LaTeX `\phi`), overload note A8; the
  emissivities and the deactivation fractions are **`ε`** (LaTeX `\varepsilon`)
  and the numerical slope and regularization scales are **`ϵ`** (LaTeX
  `\epsilon`), overload note E1.
- Calligraphic letters denote sets (`𝒟`, `𝓕_live`), the correction's monitor
  `𝓜`, the pooling conductance `𝒢`, the extended caloric relation `𝓔_i`, and
  the smooth guards `𝒫`, `𝒮`.

---

## 1. Latin letters

| Symbol | Meaning | SI units |
|---|---|---|
| `A` | Cross-sectional area of the 1D column (the mesh's cross-section measure; an algebraic unknown under a lateral-shrinkage law) | m² |
| `A_0` | Initial / calibration-reference cross-sectional area (fixed at setup) | m² |
| `A_f`, `A_n` | Area of face `f` (`A` in planar coordinates; `2π r A`, `A r²` on cylindrical and spherical meshes). On a line face `n` sits at node `n`, so the balance rows' `A_i`, `A_{i+1}` are the areas of the two faces of cell *i* (overload A3) | m² |
| `A_{αβ}`, `A_i` | Arrhenius pre-exponential factor for reaction *i* (or *αβ*) | s⁻¹ (1st-order); m³·kg⁻¹·s⁻¹ (bimolecular) |
| `a_j` | Surface absorptivity of component *j* for incident radiation (`surface_absorptivity`; `ε_j(T)` when unset, Kirchhoff) | – |
| `a_eff` | Effective (composition-weighted) surface absorptivity, `Σ_j v_j a_j` | – |
| `a`, `a_s` | Absorptivity applied to the incident flux of a `RadiativeFluxBC`: the condition's own `absorptivity` when it has one, otherwise `a_eff` of the outermost live cell (§8.3.2, §12.3.5) | – |
| `a_i`, `a_n` | Activation switch of cell *i* / node *n* of the equidistribution correction (§11; overload A9) | – |
| `B` | Spalding blowing number `ṁ''c_{p,g}/h` of the convective blowing models (§12.3.3) | – |
| `B_λ`, `B_Γ` | Film parameters of the blowing models: `B_λ = 2λ ṁ'' c_{p,g}/h` (`CouetteFilm`, `λ` its film constant, overload A6); `B_Γ = (ṁ'' c_{p,∞}/h)(Γ_∞/Γ_g)^{n_Γ}` (`PropertyRatioFilm`, `PoolingFilm`) (§12.3.3; overload B1) | – |
| `C_i` | Heat-capacity sum of cell *i*, `Σ_j m_{j,i} c_{p,j}(T̃_i)` at the trial temperature, smoothly floored at `C_{floor,i}` | J·K⁻¹ |
| `C_{floor,i}` | Capacity floor of cell *i*, `ρc_{p,floor} V_i(0)` (§5.4.3, §13.10.4) | J·K⁻¹ |
| `C_{reg,i}` | Regularizing capacity of the extended caloric relation: zero while the cell's guaranteed heat capacity is at least `2C_{floor,i}`, rising to `C_{floor,i}` as that capacity vanishes (§5.4.3, §13.10.4, §16.4) | J·K⁻¹ |
| `c_p`, `c_{p,j}` | Specific heat capacity (of component *j*) | J·kg⁻¹·K⁻¹ |
| `c_{p,j,min}` | Half the sampled minimum of the heat-capacity law of component *j*, the guaranteed capacity of the caloric extension (§13.10.4) | J·kg⁻¹·K⁻¹ |
| `c_{p,ref}` | Mass-weighted heat capacity of a cell at its initial temperature, the scale of its enthalpy-row tolerance (§13.10.5, §15.5.6) | J·kg⁻¹·K⁻¹ |
| `c`, `c_n` | Node correction of a moving mesh, `c = w − v_s`: the mesh velocity relative to the solid's (overload C1) | m·s⁻¹ |
| `D_{AB}`, `D_j`, `D_{k,∞}` | Binary diffusion coefficient (Chapman–Enskog, §5.10); `D_{k,∞}` of species *k* into the free stream (pooling film, §12.3.3) | m²·s⁻¹ |
| `𝒟` | Set of dead cells that merge into the receiver *k* in one deactivation event, outermost first (§14.3.3, §16.7) | – |
| `d`, `d_L`, `d_R`, `d_o`, `d_n`, `d_{LR}` | Distance from cell center to face (left/right, owner/neighbour) / between adjacent cell centers (overload D1) | m |
| `d_m`, `d_H`, `d_m^rel`, `d_H^rel` | Mass and energy defects of the conservation audit, and each over the largest term of its balance (§16.3.3, §16.5.4; overload D1) | kg; J; – |
| `𝐝_k` | Vector from a cell's centre to the centre of its neighbour *k* in the least-squares gradient (more dimensions, §13.6.4; overload D1) | m |
| `E`, `E_{αβ}`, `E_i` | Activation energy for reaction *i* | J·mol⁻¹ |
| `E_i` (assembly) | Energy right-hand side of cell *i*, `S_{E,i} − Σ_f s_if A_f F^E_f`: its sources plus its face energy fluxes (§1.4.4, §7.3.3, §13.3.4) | W |
| `e` | Euler's number / base of natural log (literal) | – |
| `𝓔_i(T; m)` | Extended caloric relation of cell *i*: its sensible-enthalpy inventory as a function of the temperature at the masses `m`, `Σ_j (m_j)⁺ h_j(T) + C_{reg,i}(m)(T − T_ref)` (§13.10.4, §14.3.3) | J |
| `F` | Radiative view factor (overload F1) | – |
| `F_b(x)` | Film factor `x/(e^x − 1)` of the blowing models (§12.3.3; overload F1) | – |
| `F_i` | Deformation gradient of cell *i* (meshes of more dimensions, residual level; §15.3.2; overload F1) | – |
| `F_f`, `F_{j,f}` | Face flux per unit area along the face normal (of component *j*: the mass flux the kernels write) | (flux units)·m⁻² |
| `F^E_f`, `F^E_{rel,f}` | Face energy flux per unit area along the normal, `q_f + Σ_j F_{j,f} h_j(T_f)`: conduction plus the sensible enthalpy the mass fluxes carry; `F^E_{rel,f}` its relative-flux part (§1.4.4, §7.3.3, §11.6.4, §16.5.1) | W·m⁻² |
| `𝓕_live` | Set of live boundary faces: for each boundary face of the mesh, the face that bounds the live cells on its side (§16.3.2) | – |
| `f` | Generic field, or Newton-residual function `f(T_s)` (`Φ` where the face index `f` is in the same formula) | varies |
| `G_i` | Cumulative dry pyrolysis gas produced by cell *i*: first member of the conversion-history pair (variable cross section) | kg |
| `𝒢` | Gas-exchange conductance of the pooling film without outflow, `h/c̄_{p,∞}` (§12.3.3) | kg·m⁻²·s⁻¹ |
| `g` | Gravitational acceleration (Darcy buoyancy, comparison only) | m·s⁻² |
| `g(x)` | One-sided onset of the depletion limiter, `x = ξ/ξ_th` (§6.3) | – |
| `H(·)` | Heaviside step (implemented as smooth tanh ramp; overload H2) | – |
| `H_{min,i}(T)`, `H_{max,i}(T)` | Smooth lower / upper temperature gates of reaction *i*, tanh ramps at `T_min,i`, `T_max,i` (§6.2.2; overload H2) | – |
| `H_i` | Sensible-enthalpy inventory of cell *i*, `Σ_j m_{j,i} h_j(T_i)`; the energy state of the enthalpy form (overload H2) | J |
| `H_live` | Sensible-enthalpy inventory of the live cells, `Σ_{i live} H_i` (§16.4) | J |
| `Ḣ_in` | Enthalpy rate of the live region, positive entering: `−F^E_f A_f` over the live boundary faces plus the live cells' energy sources and formation terms (§16.5.1) | W |
| `H_m` | Summed sensible enthalpy a deactivation receiver takes, less what the event loses (§16.7) | J |
| `H^lost_e`, `H^reg_e` | Sensible enthalpy deactivation event *e* removes from the live region without passing it to the receiver (vented gas, the last cell's inventory) / caloric regularization the merged cells held that the receiver does not (§16.7) | J |
| `h` (kinetics), `h_r` | Heat of reaction, per kg of first reactant | J·kg⁻¹ |
| `h_j(T)` | Sensible enthalpy of component *j* on the `T_ref` datum, `∫_{T_ref}^{T} c_{p,j} dT` (overload H1) | J·kg⁻¹ |
| `h_conv` | Convective heat-transfer coefficient (bare `h` inside §12.3.3); `h_eff = R_b h_conv` the coefficient under blowing | W·m⁻²·K⁻¹ |
| `h_m` | External mass-transfer (film) coefficient | m·s⁻¹ |
| `h_P` | Pressure-transfer coefficient (convective pressure BC) | m·Pa⁻¹·s⁻¹ |
| `Δh`, `Δh_g` | Sensible enthalpy difference `∫ c_p dT` between two temperatures, e.g. `Δh_g(T_s, T_∞)`; `Δh_c,eff` the effective heat of combustion (User Guide 2) | J·kg⁻¹ |
| `I(z)` | Radiation intensity along the beam (Beer–Lambert) | W·m⁻² |
| `I_f`, `I_up`, `I_down` | Intensity at face *f* (an algebraic unknown under Beer–Lambert) / at a cell's upstream and downstream face | W·m⁻² |
| `I_inc`, `I_e` | Intensity entering the medium through the exposed face *e*, `a q_inc` (§8.3.2, §8.3.4) | W·m⁻² |
| `I_b` | Intensity reaching the last face of the beam, transmitted out through the back (`radiation_transmitted`, §8.3.2) | W·m⁻² |
| `J` | Jacobian `∂f/∂u` of the residual (linear-algebra context only) | mixed |
| `J_j`, `J_{j,f}` | Mass flux of (gas) component *j* relative to the solid (at face *f*: the gas-transport flux, §9, §13.7) | kg·m⁻²·s⁻¹ |
| `J_j^diff`, `J_j^adv` | Diffusive / advective part of the gas flux | kg·m⁻²·s⁻¹ |
| `J_{bc,j}` | Flux of component *j* into the domain that a mass condition prescribes at a boundary face (§12.4, §13.7.5) | kg·m⁻²·s⁻¹ |
| `K`, `κ` | Permeability — **use `κ`**; `K` reserved for Gpyro-comparison text only | m² |
| `k`, `k_j` | Thermal conductivity (of component *j*) | W·m⁻¹·K⁻¹ |
| `k_eff`, `k_f` | Effective (mixture) / face-averaged thermal conductivity | W·m⁻¹·K⁻¹ |
| `k_B` | Boltzmann constant (reduced temperature `T* = k_B T/ε`) | J·K⁻¹ |
| `L` | Domain thickness (total material length); also characteristic length in convection correlations | m |
| `L_layer`, `L̇_layer` | Thickness of a layer of the equidistribution correction, from its inner bounding node to its outermost live node, and its rate of change (§11.7.6) | m; m·s⁻¹ |
| `M`, `M_j` | Molar mass (of component *j*) (overload M1) | kg·mol⁻¹ |
| `M_face`, `M_∞` | Molar mass of the gas at the face (from the pooling film's mass fractions) / of the free stream (§12.3.3) | kg·mol⁻¹ |
| `M` (linear algebra) | The integrator's diagonal mass matrix: one on the differential rows, zero on the algebraic rows (`𝖬` in §3 and §15) | – |
| `M_ref` | Reference molar mass in the Chapman–Enskog combining rules (= 0.029, air; §5.10) | kg·mol⁻¹ |
| `M_dry` | Total initial dry-solid mass in the column, `Σ_i m_{dry,i}` | kg |
| `𝓜`, `𝓜̃_i` | Monitor of the equidistribution correction, `𝓜 = α_M + β_M ρ_c/ρ_ref` (`ρ_c` a cell's condensed mass per volume, `ρ_ref` the largest initial one), and its value in cell *i* over its layer's initial mean (§11.7.4; overload M1) | – |
| `m_{j,i}` | Mass of component *j* in cell *i* — **the integrated state** | kg |
| `m_j` | Total mass of component *j* in the domain, `Σ_i m_{j,i}` | kg |
| `m_live` | Mass inventory of the live cells, `Σ_{i live} Σ_j m_{j,i}` (§16.3.1) | kg |
| `m^lost_{j,e}`, `m^lost_e` | Mass of component *j* (summed over the components) that deactivation event *e* removes from the live region without passing it to the receiver: the vented gas on a fixed mesh, the last cell's inventory (§16.3.3, §16.7); `m^vent_{j,d}` the gas mass vented from dead cell *d* (§14.3.3) | kg |
| `m^s_i(0)` | Reference condensed mass of cell *i*: its initial solid and liquid mass (the deactivation reference and the scale of its mass tolerances); `m^s_i` its current condensed mass | kg |
| `m_{dry,i}` | Initial dry-solid mass of cell *i*: second member of the conversion-history pair (code `M`), the χ denominator | kg |
| `m̂_{j,k}`, `V̂_k` | Mass of component *j* in, and volume of, the receiver *k* after a deactivation merge (§14.3.3) | kg; m³ |
| `ṁ_g`, `ṁ''`, `ṁ''_k` | Outward gas mass flux at the surface (`ṁ''` in the blowing models, `ṁ''_k` of species *k*; §12.3.3) | kg·m⁻²·s⁻¹ |
| `ṁ_out` | Outflow rate of the live region, positive leaving, `Σ_{f∈𝓕_live} Σ_j F_{j,f} A_f` (§16.3.2) | kg·s⁻¹ |
| `Ṁ_{j,i}` | Face-flux part of the mass rate of component *j* in cell *i*, `−Σ_f s_if A_f F_{j,f}` | kg·s⁻¹ |
| `mob_i`, `mob_f` | Gas mobility `κ/μ` of cell *i* (permeability over the floored viscosity) / its guarded harmonic face mean (§5.7.4, §9.5.2; overload M1) | m²·Pa⁻¹·s⁻¹ |
| `N` | Partial-pressure sum (numerator of ideal-gas pressure; overload N1) | Pa |
| `N_C` (`NC`) | Number of chemical components (type parameter) | – |
| `N_R` (`NR`) | Number of reactions (type parameter) | – |
| `N_state` | Length of the state vector (§3.1.2, §15.1.1) | – |
| `N_node`, `N_cell`, `N_blk` | Rows of a node block / of a cell block of the interleaved state on a line, and their sum, the stride (§15.1.1) | – |
| `𝖭_i` | Normal matrix of the least-squares gradient of cell *i*, `Σ_k w_k 𝐝_k 𝐝_kᵀ` (§13.6.4) | – |
| `n` | Number of finite-volume cells | – |
| `n_c` | Number of colors of the colored Jacobian (§15.7) | – |
| `n_layer` | Number of cells across a layer of the equidistribution correction; its gain is `n_layer²/τ_eq` (§11.7.2; overload N1) | – |
| `n_nodes`, `n_faces` | Number of mesh nodes / faces (= n+1 in 1D) | – |
| `n_{i,j}`, `n_s` | Reaction order of reaction *i* w.r.t. reactant *j* | – |
| `n_Γ`, `n_D` | Exponents of the property-ratio film (`n_Γ` on `Γ_∞/Γ_g`, default 0.6) and of the pooling film's face fractions (`n_D`, default 1/2) (§12.3.3) | – |
| `P` | Gas pressure in the porous medium | Pa |
| `P_ref` | Reference pressure (= 101325) | Pa |
| `Q_rxn` | Volumetric reaction heat source, `−Σ h_r r_r` | W·m⁻³ |
| `Q_rad` | Volumetric radiation absorption source | W·m⁻³ |
| `q`, `q_cond`, `q_f` | (Conductive) heat flux; at face *f*, along its normal | W·m⁻² |
| `q_BC` | Boundary-condition heat flux (positive into domain) | W·m⁻² |
| `q_inc` | Incident radiative flux of a `RadiativeFluxBC` | W·m⁻² |
| `q'''`, `q_rad'''`, `q^rad_i` | Volumetric radiation source (Beer–Lambert; `q^rad_i` its average over cell *i*, §8.3.2) | W·m⁻³ |
| `R_b` | Blowing ratio `h_eff/h` of a convective condition's blowing model, 1 at zero outflow (§12.3.3; the code's `χ`, overload B1) | – |
| `R_g` (`R`, `R_gas`) | Universal gas constant (= 8.31446261815324) | J·mol⁻¹·K⁻¹ |
| `r`, `r_i` | Reaction rate (per unit volume) of reaction *i* | kg·m⁻³·s⁻¹ |
| `r_hi`, `r_lo` | Activation bounds of the equidistribution correction, fractions of a cell's initial thickness (defaults 0.3, 0.1) | – |
| `S` | Sutherland constant (= 110.4 K, air) in the viscosity law | K |
| `S_j`, `S_{j,i}` | Species source term for component *j* (in cell *i*), `Σ_r ν_{r,j} r_r` (= `ω_{j,i}`); the rows carry the cell's reaction source `ω_{j,i} Ṽ_i` [kg·s⁻¹] | kg·m⁻³·s⁻¹ |
| `S_dry,i` | Dry-pyrolysis-gas source term in cell *i* (source of `G_i`) | kg·m⁻³·s⁻¹ |
| `S_{E,i}` | Energy source of cell *i*: its reaction heat `Q_{rxn,i} Ṽ_i` plus the in-depth radiation it absorbs (§7.3.3, §13.3.4, §16.5.1) | W |
| `S_T`, `S_{T,i}` | External temperature-rate source (manufactured solutions): `S_T` on the temperature row, `C_i S_T` on the enthalpy row (§7.3.3, §13.5.3) | K·s⁻¹ |
| `S^rad_i` | In-depth radiation absorbed by cell *i*, `I_up A_up − I_down A_down` (§8.3.3) | W |
| `s_if` | Orientation of face *f* relative to cell *i* (written `s_{if}` or `s_{i,f}`; `s_f` in §6, where *i* is the reaction index and the cell carries no index): +1 when its normal points out of *i* (*i* its owner), −1 otherwise | – |
| `s_k` | Diffusivity ratio `D_{k,∞}/α_{th,∞}` of species *k* in the pooling film (§12.3.3) | – |
| `T` | Absolute temperature | K |
| `T_s` | Surface (boundary face) temperature | K |
| `T_cell` | Temperature at the cell center adjacent to a boundary | K |
| `T_∞` | Ambient / surroundings temperature | K |
| `T_ref` | Reference temperature for enthalpy datum (= 298.15) | K |
| `T_f` | Face temperature of the enthalpy the gas fluxes carry: the distance-weighted interpolation `(d_n T_o + d_o T_n)/(d_o + d_n)` at an interior face, the surface temperature on a boundary face (§13.6.2) | K |
| `T̄_f`, `T_face` | Mean `(T̃_o + T̃_n)/2` of the two cells' trial temperatures inside the Fickian flux (written `T_face = (T_L + T_R)/2` in §3.4 and §9.1; §13.6.2, §13.7.2) | K |
| `T̃_f` | Face temperature of the relative fluxes: the upwind van Albada reconstruction on the stencil of the mass flux it carries (§11.6.4, §13.6.2) | K |
| `T_trial`, `T̃_i` | Trial temperature: `T` smoothly clamped into the property bracket `[200, 3000]` K with 1 K transitions | K |
| `T_m` | Merged temperature of a deactivation receiver: the root of the extended caloric relation at `H_m` and the merged masses (§14.3.3, §16.7) | K |
| `T_min,i`, `T_max,i` | Lower / upper temperature gate for reaction *i* | K |
| `T*` | Reduced temperature `k_B T / ε_AB` (collision integral) | – |
| `t` | Time | s |
| `u` | Flat DAE state vector | mixed |
| `du` | Time derivative of state vector | mixed |
| `V`, `V_i`, `V_i(0)` | Cell volume (of cell *i*; initial) | m³ |
| `Ṽ_i` | Smoothly floored cell volume, `smooth_floor(V_i, 10⁻⁸ V_i(0))`, the volume the concentrations and the reaction sources use (code `V_eval`; §5.13.1, §13.10.2, §13.10.3) | m³ |
| `v_g` | Bulk (Darcy, superficial) gas velocity at a face | m·s⁻¹ |
| `v_s` | Solid (skeleton) velocity; at the nodes `v_s = w − c` | m·s⁻¹ |
| `v_j` | Volume fraction of component *j* | – |
| `v_j^solid` | Solid volume fraction (excludes non-swelling gases) | – |
| `v_m` | Matrix volume fraction, `Σ_j γ_j ξ_j/ρ_j` (§5.1.3) | – |
| `v_out` | Outward Darcy velocity a pressure condition drives at a boundary face (§12.5, §13.7.5) | m·s⁻¹ |
| `v_{rel,f}` | Velocity of the solid relative to face *f* along its normal, `(v_s − w)·n̂_f = −c_f n̂_f` (§11.6.1) | m·s⁻¹ |
| `W` | The integrator's iteration matrix `J − M/(γ_RK Δt)` | – |
| `w`, `w_n` | Mesh node velocity (at node *n*), `w = v_s + c`; an algebraic unknown on a moving mesh (overload W1) | m·s⁻¹ |
| `w_k` | Weight `\|𝐝_k\|⁻²` of neighbour *k* in the least-squares gradient (§13.6.4); a quadrature weight in §16.5.3 (overload W1) | m⁻²; – |
| `X_α` | Volume fraction (FDS/Gpyro comparison notation; internally `v_j`) | – |
| `x`, `y` | Transverse/lateral coordinates (multi-dimensional meshes and comparison text) | m |
| `x_J` | Join of the depletion limiter's onset, in units of `ξ_th` (`join`, default 1; §6.3.2) | – |
| `Y_j` | Mass fraction of component *j*, `ξ_j/ρ` (`Y_k` of species *k* at the face in the pooling film, §12.3.3) | – |
| `Z` | Pre-exponential factor (Gpyro/FDS comparison notation; internally `A`) | varies |
| `z` | Axial / through-thickness coordinate (primary spatial variable) | m |
| `z_n` | Position of mesh node *n* (a differential unknown on a moving mesh) | m |
| `z_i^c` | Cell-center position of cell *i* | m |
| `Δz`, `Δz_i`, `Δz_i(0)` | Cell thickness (of cell *i*; initial) | m |

---

## 2. Greek letters

| Symbol | Meaning | SI units |
|---|---|---|
| `α` (`α_eff`) | Volumetric absorption (extinction) coefficient | m⁻¹ |
| `α_j` | Mass-basis absorption coefficient of component *j* | m²·kg⁻¹ |
| `α_conv` | Conversion / extent variable (Gpyro kinetics comparison) | – |
| `α_{th,∞}` | Thermal diffusivity of the free stream (pooling film, §12.3.3; overload A5) | m²·s⁻¹ |
| `α_M`, `β_M` | Monitor coefficients of the equidistribution correction, `monitor = (α_M, β_M)` (defaults 1, 0) | – |
| `β` | Weighting parameter of the WEIGHTED mixing rule, `k = βk_∥ + (1−β)k_series` | – |
| `γ` (`γ_j`) | Swelling factor of component *j* (1 = solids, 0 = non-swelling gas) | – |
| `γ_RK` | Stage coefficient of the integrator in `W = J − M/(γ_RK Δt)` | – |
| `Γ` | Thermal property ratio `k/c_p` of a gas (property-ratio and pooling films, §12.3.3; overload G1) | kg·m⁻¹·s⁻¹ |
| `Γ_i` | Volume rate the node corrections sweep through the faces of cell *i*, `A_{i+1} c_{i+1} − A_i c_i` on a line (§11.5.2; overload G1) | m³·s⁻¹ |
| `Δ` | Forward-difference / change operator (prefix) | varies |
| `Δ_i` | Thickness of cell *i* along the Beer–Lambert beam, the sum of its two floored centre-to-face distances (`Δz_i` away from the floors; §8.3.3, §13.5.2) | m |
| `δ_{ij}` | Kronecker delta | – |
| `δ_i` | Positive-part mass scale of cell *i*, `reltol · ε_m · m^s_i(0)` (its mass tolerance) | kg |
| `δ_b` | Scale of the softplus finite base of the depletion limiter (`base_scale`, default 10⁻⁸; §6.3.2) | kg·m⁻³ |
| `δ_mov`, `δ_BL`, `δ_χ` | Configuration flags of the state layout: one on a moving mesh / under Beer–Lambert radiation / under a lateral-shrinkage law, zero otherwise (§3.1.2, §15.1.1) | – |
| `ε`, `ε_j` | Emissivity (of component *j*) | – |
| `ε_eff` | Effective (composition-weighted) emissivity | – |
| `ε_m` | Deactivation threshold: fraction of a cell's reference condensed mass below which it deactivates (`threshold`, default 1e-2; overload E1) | – |
| `ε_band` | Near-threshold band of deactivation, `ε_m ≤ ε_band < 1` (`band`) | – |
| `ε_AB` | Lennard-Jones well depth (collision integral) | J |
| `ϵ_ξ`, `ϵ_T`, `ϵ_G` | Slope scales of the van Albada reconstruction for the concentrations, the temperature and the history values: the square of a tolerance over the upwind cell's reference thickness (§11.6.2; overload E1) | (field units)²·m⁻² |
| `ϵ_LS` | Relative regularization of the least-squares normal matrix, 10⁻⁸ (`GRADIENT_REGULARIZATION`; §13.6.4; overload E1) | – |
| `η` | Guard of the distance-weighted harmonic face means, 10⁻²⁰ in the property's units (`FACE_MEAN_GUARD`; §1.4.4, §5.12, §7.1.2, §9.1.4); in §17 the relative finite-difference step (overload E2) | property units |
| `η_m` | Guard of the face mobility mean, 10⁻³⁰ (`MOBILITY_GUARD`; §9.5.2; overload E2) | m²·Pa⁻¹·s⁻¹ |
| `θ`, `θ_i` | Volumetric strain rate (dilation) of cell *i*, `(1/V) dV/dt` (bold `θ` in §17 is the parameter vector; overload A2) | s⁻¹ |
| `θ_A` | Column-global **radial** volumetric strain rate, `(1/A) dA/dt` (an algebraic unknown under a lateral-shrinkage law) | s⁻¹ |
| `θ_L` | **Axial** volumetric strain rate (`θ − θ_A`) | s⁻¹ |
| `θ_i^j` | Stoichiometric coefficient (ThermaKin comparison notation; internally `ν_{i,j}`) | – |
| `κ` (`κ_j`, `κ_eff`) | Permeability (Darcy) of component *j* / of the cell's mixture | m² |
| `κ_s` | Solid absorption coefficient (FDS comparison notation; internally `α`) | m⁻¹ |
| `Λ` | Face buoyancy of the pooling film relative to the hot free stream, `[1 − (M_face/T_s)(T_∞/M_∞)]/(1 − T_∞/T_s)` (§12.3.3; overload B1) | – |
| `Λ_c`, `Λ_w` | Centre and width of the pooling factor's logistic switch in `Λ` (`PoolingFilm` keywords `Λc`, `Λw`, defaults 0.35, 0.10; §12.3.3) | – |
| `λ` (`λ_j`, `λ_face`, `λ_eff`) | Gas transfer (diffusion) coefficient | m²·s⁻¹ |
| `λ_τ` | Factor of the layer equidistribution time, `τ_eq = λ_τ L_layer/\|L̇_layer\|` (`ActivationTime`, default 1.8; overload A6) | – |
| `ℓ` | Characteristic particle/pore diameter (Carman–Kozeny); in §17 the scalar loss | m |
| `μ`, `μ_i` | Gas dynamic (Sutherland) viscosity (of cell *i*) | Pa·s |
| `μ_ref` | Reference viscosity (= 1.716×10⁻⁵ at 273.15 K, air) | Pa·s |
| `ν` | Kinematic viscosity (Gpyro Darcy comparison) | m²·s⁻¹ |
| `ν_{i,j}` | Stoichiometric coefficient (yield) of component *j* in reaction *i* (negative reactant, positive product) | – |
| `ξ_j`, `ξ_{j,i}` | Mass concentration of component *j* (in cell *i*), `m_{j,i}/Ṽ_i` — **constitutive species variable** | kg·m⁻³ |
| `ξ⁺_{j,i}` | Positive-part concentration `(m_{j,i})⁺/Ṽ_i` at which the property laws are evaluated | kg·m⁻³ |
| `ξ_th` | Concentration scale of the depletion limiter (`DepletionLimiter` `threshold`, default 1 kg·m⁻³) | kg·m⁻³ |
| `ξ_{tol,i}` | Mass tolerance of cell *i* over its initial volume, the positivity scale of the relative-flux face states (§11.6.2, §11.6.3) | kg·m⁻³ |
| `ρ` | Mixture / bulk density (`Σ ξ_j` for total; mixture-density formula for effective) | kg·m⁻³ |
| `ρ_j` | Bulk (pure-phase) density of component *j* | kg·m⁻³ |
| `ρ_{i,j}` | Intrinsic / skeletal density of component *j* (skeletal volume fraction only) | kg·m⁻³ |
| `ρ_g` | Gas-phase density (ideal-gas law) | kg·m⁻³ |
| `ρc_p^eff` | Effective (mixture) volumetric heat capacity, every component, `Σ_j ξ_j c_{p,j}` (discretely `C_i/V_i`) | J·m⁻³·K⁻¹ |
| `σ` | Stefan–Boltzmann constant (= 5.670374419×10⁻⁸) | W·m⁻²·K⁻⁴ |
| `σ_s` | Scattering coefficient (not modeled) | m⁻¹ |
| `σ_AB` | Lennard-Jones collision diameter | Å |
| `τ_eq` | Equidistribution time of a layer in the correction (keyword `τ`; overload T1) | s |
| `τ_i` | Optical thickness of cell *i*, `α_i Δz_i` (always subscripted or called optical; overload T1) | – |
| `φ` | Porosity (gas void fraction), `1 − Σ_{S,L} ξ_j/ρ_{i,j}` (U+03C6, LaTeX `\varphi`; overload A8) | – |
| `ϕ_i` | Multiplier of cell *i* in the potential closure (multi-dimensional meshes only; U+03D5, LaTeX `\phi`; overload A8) | m²·s⁻¹ |
| `χ`, `χ_i` | **Pyrolysis progress** in cell *i*, `G_i/m_{dry,i}` (cumulative dry-gas mass released per unit initial dry-solid mass) | – |
| `χ̄` | Column conversion, `Σ_i G_i / Σ_i m_{dry,i}` (the `m_dry`-weighted mean of `χ_i`) | – |
| `χ̄_c` | Column conversion smoothly clamped into `[−10⁻³, 1 + 10⁻³]` (`LATERAL_CLAMP_WIDTH`), the argument of the lateral-shrinkage law (§10.4) | – |
| `ψ`, `ψ_0` | Pooling factor of the pooling film, the share of the no-outflow exchange the face keeps, and its floor (`PoolingFilm` keyword `ψ0`, default 0.59; §12.3.3; overload B1) | – |
| `ψ_VA(a, b; ϵ)` | Smooth van Albada slope of the relative-flux reconstruction (§11.6.2) | (field units)·m⁻¹ |
| `Φ` | Generic cell field in the two-point gradient, where the face index `f` shares the formula (§13.6.1); `Φ_i(ξ)` the reactant concentration product of reaction *i* (§6.4); in §17 the terminal part of a loss | varies |
| `Ω(t)`, `Ω_i(t)` | Live region of the column (§16.2); region of cell *i* (§7.2, §11.2.1) | – |
| `Ω_D` | Diffusion collision integral (Neufeld) | – |
| `ω_{j,i}` | Net volumetric production rate of component *j* by reaction in cell *i*, `Σ_r ν_{r,j} r_r` (= `S_{j,i}`; the code's `ω`); the cell's reaction source is `ω_{j,i} Ṽ_i` | kg·m⁻³·s⁻¹ |
| `ω̇'''` | Volumetric reaction rate (Gpyro comparison notation; internally `r`) | kg·m⁻³·s⁻¹ |

---

## 3. Subscripts

| Subscript | Meaning |
|---|---|
| `j` | Component / species index (1…N_C) |
| `i` | Cell index (1…n); also reaction index in kinetics context (see overload I1) |
| `n` | Node index (1…n+1) in node quantities (`z_n`, `w_n`, `c_n`); bare `n` is the cell count |
| `r`, `αβ` | Reaction index; in §15, where stated, a row of the state vector (`f_r`, `s_r`, `J_{rc}`) |
| `k` | Layer or cell index inside a sweep (closure and correction context); the receiver of a merge (§14, §16.7); a neighbour in the least-squares gradient (§13.6.4); a species of the pooling film (§12.3.3); a quadrature node (§16.5.3) |
| `q` | Parameter index (§17); a node of a face in the multi-dimensional sweep vectors `s_{q,f}` (§11.6.5) |
| `e` | Exposed face of the Beer–Lambert beam (`I_e`, `A_e`, §8); a deactivation event in the ledger (`m^lost_e`, §16) |
| `b` | Last face of the beam, at the back (`I_b`, `A_b`, §8) |
| `g` | Gas phase / gaseous component |
| `s` | Solid / condensed phase |
| `ℓ`, `liq` | Liquid phase |
| `L`, `R` | Left / right of a face |
| `o`, `n` (face context) | Owner / neighbour cell of a face |
| `f` | Face-centered quantity |
| `face` | Face-averaged property |
| `eff` | Effective (mixture) property |
| `cell` | Cell-centered quantity |
| `s` (boundary) | Surface (boundary face) — e.g. `T_s`, `ξ_{j,s}` |
| `up`, `down` | Upstream / downstream face of a cell along the Beer–Lambert beam |
| `inc` | Incident (radiation entering the exposed face) |
| `∞` | Ambient / surroundings; the free stream of a blowing model |
| `ref` | Reference value |
| `in` | Cell-entrance (radiation); entering the live region (`Ḣ_in`) |
| `out` | Leaving (`ṁ_out`, `v_out`) |
| `live` | Restricted to the live cells or the live boundary faces (§16) |
| `layer` | A layer of the equidistribution correction (`n_layer`, `L_layer`) |
| `rel` | Relative to the moving mesh (`v_{rel,f}`, `F^E_{rel,f}`) |
| `tol` | Tolerance scale (`ξ_{tol,i}`) |
| `floor`, `reg` | Floor / regularizing part of the caloric extension (`C_{floor,i}`, `C_{reg,i}`) |
| `ext` | External source |
| `min`, `max` | Lower / upper bound |
| `dry` | Dry-pyrolysis-gas (excludes water vapor) |
| `init`, `0` | Initial value |
| `top`, `bot` | Top (z=L) / bottom (z=0) boundary |
| `rxn` | Reaction-source contribution |
| `cond`, `conv`, `rad` | Conduction / convection / radiation contribution |
| `∥`, `series` | Parallel / series mixing limit |
| `total` | All-component total (e.g. `ρ_total = Σ_j ξ_j`) |

## 4. Superscripts and accents

| Mark | Meaning |
|---|---|
| `°`, `^o`, `_init` | Value at previous time step / initial value |
| `(0)` | Value at the initial (reference) state, e.g. `V_i(0)`, `Δz_i(0)`, `m^s_i(0)` |
| `+`, `−` (radiation) | Forward / backward hemispheric intensity (two-flux, FDS comparison) |
| `⁺` | One-sided positive part, e.g. `ξ⁺`, `(m)⁺` (zero at and below zero, the identity above the scale `δ_i`); never a state after an event, which carries a hat |
| `^` (hat) | The receiver's state after a deactivation merge (`m̂_{j,k}`, `V̂_k`, §14.3.3); `n̂_f` the unit normal of face *f* |
| `solid` | Restricted to solid/liquid (non-swelling) components |
| `diff`, `adv` | Diffusive / advective flux component |
| `E` | Energy (the face energy flux `F^E_f`) |
| `rad` | In-depth radiation (`q^rad_i`, `S^rad_i`) |
| `lost`, `reg` | What a deactivation event removes from the live region / the caloric regularization it writes to its record (`m^lost_e`, `H^lost_e`, `H^reg_e`) |
| `rel` (on a defect) | Relative defect (`d_m^rel`, `d_H^rel`) |
| `^c` | Cell-center (e.g. `z_i^c`) |
| `‾` (overbar) | Column / spatial average (e.g. `χ̄`, `J̄_g`), a two-point face mean (`T̄_f`), or Gpyro-style mixture average (`ρ̄`); `c̄_{p,∞}` the free-stream heat capacity averaged between `T_∞` and `T_s` |
| `~` (tilde) | Normalized, clamped or guarded trial-state quantity (e.g. `𝓜̃`, `T̃_i` the trial temperature, `Ṽ_i` the floored volume); on a face temperature, `T̃_f`, the relative-flux reconstruction |
| `·` (overdot) | Time rate / flux (e.g. `ṁ`, `Ṁ`, `ω̇`) |
| `'`, `''`, `'''` | Per-length / per-area / per-volume (blowing-model fluxes `ṁ''`; FDS/Gpyro comparison notation) |

## 5. Operators

| Operator | Meaning |
|---|---|
| `∂/∂t` | Partial time derivative (Eulerian, fixed point) |
| `D_s/Dt` | Solid-following time derivative, `∂/∂t + v_s·∇` (the mesh-frame rate where the mesh moves with the solid) |
| `d/dt` | Total time derivative; of a cell inventory, following the cell |
| `∇`, `∇·` | Gradient / divergence (1D: `∂/∂z`, `(1/V_i) Σ_f s_if A_f F_f`) |
| `∇²` | Laplacian / second derivative |
| `Σ`, `∫` | Sum / integral |
| `clamp(x,a,b)` | Clamp `x` to `[a,b]` |
| `smooth_clamp`, `smooth_floor`, `(·)⁺`, `𝒫(x, δ)` | The C¹ guards of trial states (§13.10.3): clamp with cubic joins, floor, one-sided positive part (written `𝒫(x, δ)` with its scale `δ` in §11.6.3) |
| `𝒮(x)` | Smoothstep `x²(3 − 2x)` clamped to `[0, 1]` (the correction's activation, §11.7.3) |
| `tanh` | Hyperbolic tangent (smooth temperature gates; the depletion onset above its join) |
| `‖·‖_2`, `‖·‖_W` | Euclidean norm / norm weighted by the matrix `W` of a least-squares loss (§17, §18) |
| `diag(v)` | Diagonal matrix from vector `v` |
| `J`, `M`, `W` | Jacobian / mass matrix / iteration matrix (§15) |

---

## 6. Physical constants (fixed values)

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | `R_g` | 8.31446261815324 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | `σ` | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference pressure | `P_ref` | 101325 | Pa |
| Reference temperature | `T_ref` | 298.15 | K |
| Reference molar mass (air) | `M_ref` | 0.029 | kg·mol⁻¹ |
| Reference viscosity (air, 273.15 K) | `μ_ref` | 1.716×10⁻⁵ | Pa·s |
| Sutherland constant (air) | `S` | 110.4 | K |
| Kozeny constant (spheres) | — | 180 | – |

### 6.1 Numerical constants and defaults

| Quantity | Symbol | Value | Units | Source |
|---|---|---|---|---|
| Property bracket | `[T_min, T_max]` | [200, 3000] | K | `T_PROPERTY_MIN`, `T_PROPERTY_MAX` |
| Width of the trial clamp's joins | — | 1 | K | `TRIAL_TEMPERATURE_WIDTH` |
| Floor of the volumetric heat capacity | `ρc_{p,floor}` | 1.0×10³ | J·m⁻³·K⁻¹ | `RHO_CP_FLOOR` |
| Porosity floor (smooth) | — | 10⁻¹² | – | `POROSITY_FLOOR` |
| Viscosity floor (smooth) | — | 10⁻¹² | Pa·s | `MU_FLOOR` |
| Guard of the harmonic face means | `η` | 10⁻²⁰ | property units | `FACE_MEAN_GUARD` |
| Guard of the face mobility mean | `η_m` | 10⁻³⁰ | m²·Pa⁻¹·s⁻¹ | `MOBILITY_GUARD` |
| Regularization of the least-squares gradient | `ϵ_LS` | 10⁻⁸ | – | `GRADIENT_REGULARIZATION` |
| Half-width of the clamp band of `χ̄_c` | — | 10⁻³ | – | `LATERAL_CLAMP_WIDTH` |
| Deactivation threshold | `ε_m` | 10⁻² | – | `solve(...; threshold)` |
| Near-threshold band under a moving correction | `ε_band` | 0.1 | – | `CORRECTION_BAND` |
| Activation bounds | `(r_hi, r_lo)` | (0.3, 0.1) | – | `EquidistributionCorrection` |
| Layer-time factor | `λ_τ` | 1.8 | – | `ACTIVATION_TIME_FACTOR` |
| Fallback equidistribution time | `τ_eq` | 720 | s | `EQUIDISTRIBUTION_TIME` |
| Depletion-limiter scale | `ξ_th` | 1.0 | kg·m⁻³ | `DepletionLimiter()` |
| Join of the depletion onset | `x_J` | 1.0 | – (units of `ξ_th`) | `DepletionLimiter()` |
| Scale of the finite power base | `δ_b` | 10⁻⁸ | kg·m⁻³ | `DepletionLimiter()` |

---

## 7. Overloaded / reserved symbols — resolutions (BINDING)

These symbols are overloaded somewhere in the source subsystems. The resolution
below is authoritative; writers must use the disambiguating subscript/context.

- **χ (chi) — A1 [pyrolysis progress vs. cross-section ratio].** The variable
  cross section defines `χ_i = G_i/m_{dry,i}` as the **pyrolysis progress**
  (dimensionless, cumulative dry-gas mass released per unit initial dry-solid
  mass) and `χ̄ = Σ G_i/Σ m_{dry,i}` as the column conversion; it is easily
  misread as the cross-section area ratio. **Resolution:** `χ` is *always*
  pyrolysis progress; the cross-section ratio is `A/A_0 = law(χ̄)` and is never
  written as `χ`. No time-derivative operator carries a `χ` subscript. The
  blowing ratio, which the code calls `χ`, is written `R_b` (B1).

- **θ (theta) — A2 [strain rate vs. stoichiometry vs. parameters].** Internally
  `θ` is the **volumetric strain rate** (s⁻¹), with `θ_A` radial and `θ_L`
  axial. ThermaKin uses `θ_i^j` for stoichiometric coefficients.
  **Resolution:** use `ν_{i,j}` for stoichiometry everywhere; reserve `θ` for
  strain rate. `θ_i^j` appears only in the ThermaKin-comparison passage and is
  immediately mapped to `ν_{i,j}`. The one exception is §17 (sensitivity), where
  the bold `θ ∈ ℝ^p` is the vector of parameters and `θ_q` its *q*-th entry;
  the strain rate does not appear in that chapter.

- **A — A3 [area vs. pre-exponential factor].** `A` is the **cross-sectional
  area** and `A_f` the area of face `f`; the Arrhenius pre-exponential is `A_i` /
  `A_{αβ}` (always subscripted by reaction). **Resolution:** bare `A` = area in
  physics chapters; `A_f` = face area; `A_i` = pre-exponential (kinetics). The
  FDS/Gpyro pre-exponential `Z` is mapped to `A_i` on first use. In the
  mesh-motion rows on a line (§10, §11, §15) a face is named by its node, so
  `A_i` and `A_{i+1}` there are the areas of the faces at nodes `i` and `i + 1`,
  the two faces of cell `i`; a reaction subscript never appears in those rows.

- **κ vs. K — A4 [permeability].** Use **`κ`** for permeability throughout.
  `K` (Gpyro permeability) appears only in comparison text and is mapped to `κ`.

- **α — A5 [absorption vs. conversion vs. absorptivity].** `α` (or `α_eff`) =
  **volumetric absorption coefficient** (m⁻¹); `α_j` = mass-basis absorption
  coefficient (m²·kg⁻¹). Conversion (Gpyro) is `α_conv`; FDS `κ_s` maps to `α`.
  The absorptivity of a surface for incident radiation is **`a`** (`a_j`,
  `a_eff`, `a_s`), never `α`; it equals `ε_eff` under Kirchhoff's law unless a
  component sets its own `surface_absorptivity`. The correction's monitor
  coefficients are `α_M`, `β_M`, always subscripted, and the free stream's
  thermal diffusivity in the pooling film is `α_{th,∞}`.

- **λ — A6 [gas transfer vs. anything else].** `λ` is reserved for the **gas
  transfer (diffusion) coefficient** (m²·s⁻¹). Do not use `λ` for thermal
  conductivity (that is `k`) or eigenvalues. Two exceptions carry a code keyword
  `λ`: the factor of `ActivationTime` (default 1.8), written **`λ_τ`** wherever a
  formula needs it (prose may call it "the layer-time factor λ = 1.8" when
  `ActivationTime` is named in the same sentence); and the film constant of
  `CouetteFilm` (default 1/2), which appears only inside `B_λ = 2λ ṁ'' c_{p,g}/h`
  in §12.3.3 and the User Guide's blowing passage, always with `CouetteFilm`
  named. The layer-time factor is never written `Λ`; `Λ` is the pooling film's
  face buoyancy (B1).

- **ρ subscripts — A7 [bulk vs. intrinsic].** `ρ_j` = bulk density (mixing,
  mixture density); `ρ_{i,j}` = intrinsic/skeletal density (skeletal volume
  fraction only). Never write a bare `ρ_i` for skeletal density without the
  explicit "intrinsic" qualifier; prefer `ρ_{i,j}`.

- **h — H1 [heat of reaction; sign; sensible enthalpy].** `h` (kinetics) and
  `h_r` = heat of reaction per kg of first reactant, with **h > 0 endothermic**
  internally; ThermaKin publishes `h > 0 = exothermic` (Gpyro and FDS do not). The conversion to the
  volumetric source carries the minus sign: `Q_rxn = −Σ h_r r_r`. With a
  *component* subscript, `h_j(T)` is the component's **sensible enthalpy**
  `∫_{T_ref}^{T} c_{p,j} dT`; the two never share a subscript (`r` reaction, `j`
  component). The sensible enthalpy on the datum is always `h_j(T)`, never
  `Δh_j(T)`; `Δh` is reserved for a difference between two stated temperatures
  (`Δh_g(T_s, T_∞)`) and for the effective heat of combustion `Δh_c,eff`. The
  convective coefficient is `h_conv`, the mass coefficient `h_m`, the pressure coefficient
  `h_P`, all subscripted (inside §12.3.3 the convective coefficient of the
  blowing formulas is the bare `h`, and `h_eff = R_b h`). A cell thickness is
  `Δz_i`, never `h_i`, even where the code names it `h`. In §17 `h_q` is the
  finite-difference step of parameter `q`.

- **H — H2 [Heaviside vs. inventory].** `H(·)` with an argument is the
  Heaviside step (realized as a tanh ramp), and `H_{min,i}(T)`, `H_{max,i}(T)`
  are the smooth temperature gates of reaction `i` (§6.2.2); `H_i` with a cell
  subscript is the sensible-enthalpy inventory of cell `i`, and `H_live`, `H_m`,
  `H^lost_e`, `H^reg_e` are the ledger's and the deactivation event's
  sensible enthalpies (§16).

- **i — I1 [cell index vs. reaction index].** `i` is the **cell index** in
  discretization/geometry/mesh-motion chapters; in the kinetics chapter `i` is
  the **reaction index** (with `j` the component index). Each chapter states its
  convention in its first paragraph. Where both appear together, use `i` = cell,
  `r` = reaction.

- **φ — A8 [porosity only].** `φ` (U+03C6, LaTeX `\varphi`) is **porosity**; do
  not reuse it for view factor (that is `F`), a potential or a generic field
  (that is `f`, or `q` in a reconstruction). `φ` denotes the intrinsic-density
  form `1 − Σ_{S,L} ξ_j/ρ_{i,j}`, one minus the skeletal volume fraction; it is
  evaluated unclamped where admissibility is tested and floored smoothly at
  10⁻¹² where the pressure closure divides by it. The all-phase gas volume
  fraction `Σ v_j^gas` is a different quantity and is not written `φ`. The
  potential closure's per-cell multiplier (code `φ`) is written `ϕ_i` (U+03D5,
  LaTeX `\phi`); since the two glyphs are easily confused, `ϕ` appears only with
  its cell subscript and only in the multi-dimensional closure. The blowing
  models' film parameter, which the code calls `φ`, is written `B_λ` or `B_Γ`
  (B1).

- **S — S1 [source vs. Sutherland constant].** `S` with a subscript
  (`S_j`, `S_dry`, `S_{E,i}`, `S_T`, `S^rad_i`) is a source term; bare `S` is
  the Sutherland constant in the viscosity law only. The lowercase `s_if` is the
  face orientation sign; `s_r` (row subscript, §15.5.4) is the equation scale of
  a row and `s_k` (§12.3.3) the pooling film's diffusivity ratio. In §17,
  `S_{kq} = ∂y_k/∂θ_q` is the sensitivity matrix and `S_q^1`, `S_q^T` the
  first-order and total Sobol indices; no source term appears there.

- **ε — E1 [emissivity vs. thresholds vs. numerical scales].** `ε`, `ε_j`,
  `ε_eff` are emissivities; `ε_m` (deactivation threshold) and `ε_band` (band)
  are fractions of a cell's reference condensed mass; `ε_AB` is the
  Lennard-Jones well depth. The subscript always disambiguates. The numerical
  scales of the van Albada slope (`ϵ_ξ`, `ϵ_T`, `ϵ_G`, and the `ϵ` argument of
  `ψ_VA`) and of the least-squares regularization (`ϵ_LS`) use the lunate `ϵ`
  (LaTeX `\epsilon`); the physical quantities and the deactivation fractions use
  `ε` (LaTeX `\varepsilon`).

- **η — E2 [guards vs. finite-difference step].** `η` is the guard added to the
  denominator of a guarded harmonic face mean, in the units of the averaged
  property (10⁻²⁰); the mobility mean's guard is `η_m` (10⁻³⁰ m²·Pa⁻¹·s⁻¹). In
  §17 `η` is the relative finite-difference step and `η_abs` its absolute
  floor; no face mean appears there.

- **τ — T1 [optical thickness vs. equidistribution time].** `τ_i` with a cell
  subscript, or "optical thickness" in words, is `α_i Δz_i`; `τ_eq` in the
  mesh-correction and burnout chapters (§11, §14) is a layer's equidistribution
  time [s]. A time variable of integration is written `t`, never `τ`.

- **M — M1 [molar mass vs. mass matrix vs. history vs. monitor].** `M`, `M_j`
  (and `M_face`, `M_∞`) is molar mass; `M` in the linear-algebra context of §15
  is the integrator's mass matrix (sans-serif `𝖬` in §3 and §15). The second
  member of the conversion-history pair, which the code calls `M`, is written
  `m_{dry,i}`, its rate `ṁ_{dry,i}`; the correction's monitor is the
  calligraphic `𝓜` (`𝓜̃_i` normalized). The gas mobility `κ/μ` is written
  `mob_i` (`mob_f` at a face), never `m`, which is a mass.

- **c — C1 [heat capacity vs. correction].** The specific heat capacity is
  always `c_p` (`c_{p,j}`); bare `c` (or `c_n` at node `n`) is the node correction
  `w − v_s`.

- **a — A9 [absorptivity vs. activation].** `a_j`, `a_eff`, `a_s` (component,
  effective or surface subscript) are surface absorptivities; `a_i`, `a_n` (cell
  or node subscript, correction context only) are the correction's activation
  switches. A row's absolute tolerance is written `atol_r` (§15.5.4), never
  `a_r`.

- **N — N1 [partial-pressure sum vs. counts].** `N` in the pressure closure is
  the partial-pressure sum; every count carries a subscript: `N_C`, `N_R`,
  `N_state`, `N_node`, `N_cell`, `N_blk`. The number of cells across a layer of
  the equidistribution correction is `n_layer`, in the gain `n_layer²/τ_eq`.

- **Blowing — B1 [blowing ratio, film parameters, pooling].** The code's
  blowing ratio `χ` is written **`R_b`** (`h_eff/h`) and its film parameter `φ`
  is written **`B_λ`** (`CouetteFilm`) or **`B_Γ`** (property-ratio and pooling
  films), because `χ` is the pyrolysis progress (A1) and `φ` the porosity (A8);
  `F_b(x) = x/(e^x − 1)` is the film factor, `ψ` the pooling factor (`ψ_0` its
  floor) and `Λ` the face buoyancy (`Λ_c`, `Λ_w` the switch constants). Code
  keyword names (`λ`, `n_Γ`, `n_D`, `ψ0`, `Λc`, `Λw`) are kept as they are in
  code wherever a keyword is meant.

- **d — D1 [distances vs. defects vs. least-squares vectors].** `d` with a face
  or cell subscript (`d_L`, `d_o`, `d_{U,f}`, `d_{LR}`) is a distance;
  `d_m`, `d_H` (and their relative forms `d^rel`) are the conservation audit's
  mass and energy defects (§16); the bold `𝐝_k` is the centre-to-centre vector
  of the least-squares gradient (§13.6.4).

- **F — F1 [view factor vs. fluxes vs. film factor vs. deformation].** Bare `F`
  is the radiative view factor; `F_f`, `F_{j,f}`, `F^E_f` (face subscript) are
  face fluxes; `F_b(x)` is the blowing models' film factor; `F_i` (cell
  subscript, meshes of more dimensions) is the deformation gradient.

- **Γ — G1 [thermal property ratio vs. sweep].** `Γ = k/c_p` (with `∞` or `g`,
  free stream or leaving gas) belongs to the blowing models of §12.3.3; `Γ_i`
  with a cell subscript is the volume rate the node corrections sweep (§11.5.2).

- **w — W1 [node velocity vs. weights].** `w`, `w_n` is the mesh node velocity;
  `w_k` (neighbour or quadrature-node subscript) is a weight of the
  least-squares gradient (§13.6.4) or of a quadrature rule (§16.5.3).
