# 2. Nomenclature

This is the **single, authoritative symbol table** for the Pyrolysis.jl Technical
Reference Guide and User Guide. Every writer **MUST** conform to these symbols,
meanings, and SI units verbatim. Where a symbol is used differently in different
source subsystems, the **conflict has been resolved here** and the resolution is
binding. Conflicts and overloads are called out explicitly in the
"Overloaded / reserved symbols" section at the end.

Conventions:
- The spatial coordinate is **z** (axial, through-thickness). $z = 0$ is the
  **back / bottom / substrate** (boundary tag `:bottom`, face 1, at the pinned
  node 1); $z = L$ is the **exposed surface / top** (boundary tag `:top`, face
  $n + 1$). Heat enters from the top; the material shrinks downward. The
  transverse/lateral plane carries the cross-sectional area $A$. Boundary
  conditions are keyed by tag; there are no numeric boundary ids.
- Sign convention for fluxes: a face flux $F_f$ is the transport **per unit face
  area along the face normal**, and the normal points from the face's owner cell
  to its neighbour, out of the domain on a boundary face. On a line every
  interior face and the `:top` face therefore have their normal along **+z**
  (positive flux = transport bottom → top), and the `:bottom` face along **−z**
  (positive flux = outflow through the back). The discrete divergence
  $(1/V_i) \sum_f s_{if} A_f F_f$ is positive when the quantity flows *out* of cell
  $i$.
- Heat-of-reaction sign convention: **h > 0 endothermic (cools)**,
  **h < 0 exothermic (heats)**; the volumetric heat source is $Q_{\text{rxn}} = -\sum h_r r_r$,
  so endothermic reactions give $Q_{\text{rxn}} < 0$. (This is the storage convention used
  internally; note ThermaKin publishes $h > 0 = \text{exothermic}$; Gpyro V2 and FDS use this convention — see overload
  note H1.)
- "Bulk" density $\rho_j$ is the pure-phase density used by mixing rules and mixture
  density; "intrinsic/skeletal" density $\rho_{i,j}$ (cell-wall density) is used
  ONLY in the porosity (skeletal volume fraction) formula.
- Mass concentration $\xi_j$ (kg/m³) is the **constitutive species variable**, not
  mass fraction: $\xi_j = Y_j \rho$. The integrated discrete state is the cell
  inventory $m_{j,i} = \xi_{j,i} V_i$ (kg); every kernel reads
  $\xi_{j,i} = m_{j,i}/\tilde V_i$, with $\tilde V_i$ the smoothly floored cell volume.
- **Solver notation.** ThermaKin2Ds's mass concentration and FDS's component
  density describe the same quantity as $\xi_j$: component mass per whole
  current cell volume. For a condensed component, $\rho_{s,j} = \xi_j = m_j/V$.
  FDS's total condensed density is $\rho_s = \sum_{S,L} \xi_j$. These inventory
  densities are distinct from the material density properties $\rho_j$ and
  $\rho_{i,j}$ (§§3.8, 4.8).
- The chapters write their mathematics in LaTeX; this table gives each symbol in
  Unicode. Two pairs of glyphs are distinct symbols and must not be interchanged:
  porosity is **$\varphi$** (U+03C6, LaTeX `\varphi`) and the potential-closure
  multiplier is **$\phi$** (U+03D5, LaTeX `\phi`), overload note A8; the
  emissivities and the deactivation fractions are **$\varepsilon$** (LaTeX `\varepsilon`)
  and the numerical slope and regularization scales are **$\epsilon$** (LaTeX
  `\epsilon`), overload note E1.
- Calligraphic letters denote sets ($\mathcal{D}$, $\mathcal{F}_{\text{live}}$), the correction's monitor
  $\mathcal{M}$, the pooling conductance $\mathcal{G}$, the extended caloric relation $\mathcal{E}_i$, and
  the smooth guards $\mathcal{P}$, $\mathcal{S}$.

---

## 1. Latin letters

| Symbol | Meaning | SI units |
|---|---|---|
| $A$ | Cross-sectional area of the 1D column (the mesh's cross-section measure; an algebraic unknown under a lateral-shrinkage law) | m² |
| $A_0$ | Initial / calibration-reference cross-sectional area (fixed at setup) | m² |
| $A_f$, $A_n$ | Area of face $f$ ($A$ in planar coordinates; $2\pi r A$, $A r^{2}$ on cylindrical and spherical meshes). On a line face $n$ sits at node $n$, so the balance rows' $A_i$, $A_{i+1}$ are the areas of the two faces of cell *i* (overload A3) | m² |
| $A_{\alpha\beta}$, $A_i$ | Arrhenius pre-exponential factor for reaction *i* (or *αβ*) | s⁻¹ |
| $a_j$ | Surface absorptivity of component *j* for incident radiation (`surface_absorptivity`; $\varepsilon_j(T)$ when unset, Kirchhoff) | – |
| $a_{\text{eff}}$ | Effective (composition-weighted) surface absorptivity, $\sum_j v_j a_j$ | – |
| $a$, $a_s$ | Absorptivity applied to the incident flux of a `RadiativeFluxBC`: the condition's own `absorptivity` when it has one, otherwise $a_{\text{eff}}$ of the outermost live cell (§8.3.2, §12.3.5) | – |
| $a_i$, $a_n$ | Activation switch of cell *i* / node *n* of the equidistribution correction (§11; overload A9) | – |
| $B$ | Spalding blowing number $\dot m''c_{p,g}/h$ of the convective blowing models (§12.3.3) | – |
| $B_{\lambda}$, $B_{\Gamma}$ | Film parameters of the blowing models: $B_{\lambda} = 2\lambda \sum_j N_j c_{p,j}(T_s)/h$ (`CouetteFilm`, $\lambda$ its film constant, overload A6); $B_{\Gamma} = (\dot m'' c_{p,\infty}/h)(\Gamma_{\infty}/\Gamma_g)^{n_{\Gamma}}$ (`PropertyRatioFilm`, `PoolingFilm`) (§12.3.3; overload B1) | – |
| $C_i$ | Heat-capacity sum of cell *i*, $\sum_j m_{j,i} c_{p,j}(\tilde T_i)$ at the trial temperature, smoothly floored at $C_{\text{floor},i}$ | J·K⁻¹ |
| $C_{\text{floor},i}$ | Capacity floor of cell *i*, $\rho c_{p,\text{floor}} V_i(0)$ (§5.4.3, §13.10.4) | J·K⁻¹ |
| $C_{\text{reg},i}$ | Regularizing capacity of the extended caloric relation: zero while the cell's guaranteed heat capacity is at least $2C_{\text{floor},i}$, rising to $C_{\text{floor},i}$ as that capacity vanishes (§5.4.3, §13.10.4, §16.4) | J·K⁻¹ |
| $c_{r,i}$ | Lineage coefficient of reaction *r* in cell *i*: its reference as a fraction of the cell's initial condensed mass (§6.5.5) | – |
| $c_p$, $c_{p,j}$ | Specific heat capacity (of component *j*) | J·kg⁻¹·K⁻¹ |
| $c_{p,j,\text{min}}$ | Half the sampled minimum of the heat-capacity law of component *j*, the guaranteed capacity of the caloric extension (§13.10.4) | J·kg⁻¹·K⁻¹ |
| $c_{p,\text{ref}}$ | Mass-weighted heat capacity of a cell at its initial temperature, the scale of its enthalpy-row tolerance (§13.10.5, §15.5.6) | J·kg⁻¹·K⁻¹ |
| $c$, $c_n$ | Node correction of a moving mesh, $c = w - v_s$: the mesh velocity relative to the solid's (overload C1) | m·s⁻¹ |
| $D_{AB}$, $D_j$, $D_{k,\infty}$ | Binary diffusion coefficient (Chapman–Enskog, §5.10); $D_{k,\infty}$ of species *k* into the free stream (pooling film, §12.3.3) | m²·s⁻¹ |
| $D_f$, $a_o$, $a_n$ | Conductivity an interior face adds to its gas enthalpy flux, $\max(0, a_o - k_f, a_n - k_f)$; the heat-capacity flux leaving the owner or the neighbour through the face times its center-to-face distance (§7.4.2) | W·m⁻¹·K⁻¹ |
| $\mathcal{D}$ | Set of dead cells that merge into the receiver *k* in one deactivation event, outermost first (§14.3.3, §16.7) | – |
| $d$, $d_L$, $d_R$, $d_o$, $d_n$, $d_{LR}$ | Distance from cell center to face (left/right, owner/neighbour) / between adjacent cell centers (overload D1) | m |
| $d_m$, $d_H$, $d_m^{\text{rel}}$, $d_H^{\text{rel}}$ | Mass and energy defects of the conservation audit, and each over the largest term of its balance (§16.3.3, §16.5.4; overload D1) | kg; J; – |
| $\mathbf{d}_k$ | Transmission image of neighbour *k* in a cell's least-squares gradient: the vector between the cell centres, moved across a material interface for the jump of the normal slope (more dimensions, §13.6.4; overload D1) | m |
| $E$, $E_{\alpha\beta}$, $E_i$ | Activation energy for reaction *i* | J·mol⁻¹ |
| $E_i$ (assembly) | Energy right-hand side of cell *i*, $S_{E,i} - \sum_f s_{if} A_f F^E_f$: its sources plus its face energy fluxes (§1.4.4, §7.3.3, §13.3.4) | W |
| $e$ | Euler's number / base of natural log (literal) | – |
| $\mathcal{E}_i(T; m)$ | Extended caloric relation of cell *i*: its sensible-enthalpy inventory as a function of the temperature at the masses $m$, $\sum_j (m_j)^{+} h_j(T) + C_{\text{reg},i}(m)(T - T_{\text{ref}})$ (§13.10.4, §14.3.3) | J |
| $F$ | Radiative view factor (overload F1) | – |
| $F_b(x)$ | Film factor $x/(e^x - 1)$ of the blowing models (§12.3.3; overload F1) | – |
| $F_i$ | Deformation gradient of cell *i* (meshes of more dimensions, residual level; §15.3.2; overload F1) | – |
| $F_f$, $F_{j,f}$ | Face flux per unit area along the face normal (of component *j*: the mass flux the kernels write) | (flux units)·m⁻² |
| $F^E_f$, $F^E_{\text{rel},f}$ | Face energy flux per unit area along the normal, $q_f + \sum_j F_{j,f} h_j(T_f)$: conduction (at an interior face with the gas transport's added conduction $D_f (T_o - T_n)/(d_o + d_n)$) plus the sensible enthalpy the mass fluxes carry; $F^E_{\text{rel},f}$ its relative-flux part (§1.4.4, §7.3.3, §11.6.4, §16.5.1) | W·m⁻² |
| $\mathcal{F}_{\text{live}}$ | Set of live boundary faces: for each boundary face of the mesh, the face that bounds the live cells on its side (§16.3.2) | – |
| $f$ | Generic field, or Newton-residual function $f(T_s)$ ($\Phi$ where the face index $f$ is in the same formula) | varies |
| $G_i$ | Cumulative dry pyrolysis gas produced by cell *i*: first member of the conversion-history pair (variable cross section) | kg |
| $\mathcal{G}$ | Gas-exchange conductance of the pooling film without outflow, $h/\bar c_{p,\infty}$ (§12.3.3) | kg·m⁻²·s⁻¹ |
| $g$ | Gravitational acceleration (Darcy buoyancy, comparison only) | m·s⁻² |
| $g(x)$ | One-sided onset of the depletion limiter, $x = \xi/\xi_{\text{th}}$ (§6.3) | – |
| $H(\cdot)$ | Heaviside step (implemented as smooth tanh ramp; overload H2) | – |
| $H_{\text{min},i}(T)$, $H_{\text{max},i}(T)$ | Smooth lower / upper temperature gates of reaction *i*, tanh ramps at $T_{\text{min},i}$, $T_{\text{max},i}$ (§6.2.2; overload H2) | – |
| $H_i$ | Sensible-enthalpy inventory of cell *i*, $\sum_j m_{j,i} h_j(T_i)$; the energy state of the enthalpy form (overload H2) | J |
| $H_{\text{live}}$ | Sensible-enthalpy inventory of the live cells, $\sum_{i\ \text{live}} H_i$ (§16.4) | J |
| $\dot H_{\text{in}}$ | Enthalpy rate of the live region, positive entering: $-F^E_f A_f$ over the live boundary faces plus the live cells' energy sources and formation terms (§16.5.1) | W |
| $H_m$ | Summed sensible enthalpy a deactivation receiver takes, less what the event loses (§16.7) | J |
| $H^{\text{lost}}_e$, $H^{\text{reg}}_e$ | Sensible enthalpy deactivation event *e* removes from the live region without passing it to the receiver (vented gas, the last cell's inventory) / caloric regularization the merged cells held that the receiver does not (§16.7) | J |
| $h$ (kinetics), $h_r$ | Heat of reaction, per kg of first reactant | J·kg⁻¹ |
| $h_j(T)$ | Sensible enthalpy of component *j* on the $T_{\text{ref}}$ datum, $\int_{T_{\text{ref}}}^T c_{p,j} dT$ (overload H1) | J·kg⁻¹ |
| $h_{\text{conv}}$ | Convective heat-transfer coefficient (bare $h$ inside §12.3.3); $h_{\text{eff}} = R_b h_{\text{conv}}$ the coefficient under blowing | W·m⁻²·K⁻¹ |
| $h_m$ | External mass-transfer (film) coefficient | m·s⁻¹ |
| $h_P$ | Pressure-transfer coefficient (convective pressure BC) | m·Pa⁻¹·s⁻¹ |
| $\Delta h$, $\Delta h_g$ | Sensible enthalpy difference $\int c_p dT$ between two temperatures, e.g. $\Delta h_g(T_s, T_{\infty})$; $\Delta h_{c,\text{eff}}$ the effective heat of combustion (User Guide 2) | J·kg⁻¹ |
| $I(z)$ | Radiation intensity along the beam (Beer–Lambert) | W·m⁻² |
| $I_f$, $I_{\text{up}}$, $I_{\text{down}}$ | Intensity at face *f* (an algebraic unknown under Beer–Lambert) / at a cell's upstream and downstream face | W·m⁻² |
| $I_{\text{inc}}$, $I_e$ | Intensity entering the medium through the exposed face *e*, $a q_{\text{inc}}$ (§8.3.2, §8.3.4) | W·m⁻² |
| $I_b$ | Intensity reaching the last face of the beam, transmitted out through the back (`radiation_transmitted`, §8.3.2) | W·m⁻² |
| $J$ | Jacobian $\partial f/\partial u$ of the residual (linear-algebra context only) | mixed |
| $J_j$, $J_{j,f}$ | Mass flux of (gas) component *j* relative to the solid (at face *f*: the gas-transport flux, §9, §13.7) | kg·m⁻²·s⁻¹ |
| $J_j^{\text{diff}}$, $J_j^{\text{adv}}$ | Diffusive / advective part of the gas flux | kg·m⁻²·s⁻¹ |
| $J_{\text{bc},j}$ | Outward flux of component *j* prescribed at a boundary face (§12.4, §13.7.5) | kg·m⁻²·s⁻¹ |
| $\kappa$ | Effective condensed-matrix permeability | m² |
| $K_j$ | Species wall diffusivity in the gas-friction equations | m²·s⁻¹ |
| $k$, $k_j$ | Thermal conductivity (of component *j*) | W·m⁻¹·K⁻¹ |
| $k_{\text{eff}}$, $k_f$ | Effective (mixture) / face-averaged thermal conductivity | W·m⁻¹·K⁻¹ |
| $k_B$ | Boltzmann constant (reduced temperature $T^* = k_B T/\varepsilon$) | J·K⁻¹ |
| $L$ | Domain thickness (total material length); also characteristic length in convection correlations | m |
| $L_{\text{layer}}$, $\dot L_{\text{layer}}$ | Thickness of a layer of the equidistribution correction, from its inner bounding node to its outermost live node, and its rate of change (§11.7.6) | m; m·s⁻¹ |
| $M$, $M_j$ | Molar mass (of component *j*) (overload M1) | kg·mol⁻¹ |
| $M_{\text{face}}$, $M_{\infty}$ | Molar mass of the gas at the face (from the pooling film's mass fractions) / of the free stream (§12.3.3) | kg·mol⁻¹ |
| $M$ (linear algebra) | The integrator's diagonal mass matrix: one on the differential rows, zero on the algebraic rows ($\mathsf{M}$ in §3 and §15) | – |
| $M_{\text{ref}}$ | Reference molar mass in the Chapman–Enskog combining rules (= 0.029, air; §5.10) | kg·mol⁻¹ |
| $M_{\text{dry}}$ | Total initial dry-solid mass in the column, $\sum_i m_{\text{dry},i}$ | kg |
| $\mathcal{M}$, $\tilde{\mathcal{M}}_i$ | Monitor of the equidistribution correction, $\mathcal{M} = \alpha_M + \beta_M \rho_c/\rho_{\text{ref}}$ ($\rho_c$ a cell's condensed mass per volume, $\rho_{\text{ref}}$ the largest initial one), and its value in cell *i* over its layer's initial mean (§11.7.4; overload M1) | – |
| $m_{j,i}$ | Mass of component *j* in cell *i* — **the integrated state** | kg |
| $m_j$ | Total mass of component *j* in the domain, $\sum_i m_{j,i}$ | kg |
| $m_{\text{live}}$ | Mass inventory of the live cells, $\sum_{i\ \text{live}} \sum_j m_{j,i}$ (§16.3.1) | kg |
| $m^{\text{lost}}_{j,e}$, $m^{\text{lost}}_e$ | Mass of component *j* (summed over the components) that deactivation event *e* removes from the live region without passing it to the receiver: the vented gas on a fixed mesh, the last cell's inventory (§16.3.3, §16.7); $m^{\text{vent}}_{j,d}$ the gas mass vented from dead cell *d* (§14.3.3) | kg |
| $m^s_i(0)$ | Reference condensed mass of cell *i*: its initial solid and liquid mass (the deactivation reference and the scale of its mass tolerances); $m^s_i$ its current condensed mass | kg |
| $M_{\text{ref},r}$ | Reference of reaction *r*: the lineage inventory of its reactant at the material point (§6.4, §6.5.5) | kg |
| $m_{\text{dry},i}$ | Initial dry-solid mass of cell *i*: second member of the conversion-history pair (code `M`), the χ denominator | kg |
| $\hat m_{j,k}$, $\hat V_k$ | Mass of component *j* in, and volume of, the receiver *k* after a deactivation merge (§14.3.3) | kg; m³ |
| $\dot m_g$, $\dot m''$, $\dot m''_k$ | Outward gas mass flux at the surface ($\dot m''$ in the blowing models, $\dot m''_k$ of species *k*; §12.3.3) | kg·m⁻²·s⁻¹ |
| $\dot m_{\text{out}}$ | Outflow rate of the live region, positive leaving, $\sum_{f\in\mathcal{F}_{\text{live}}} \sum_j F_{j,f} A_f$ (§16.3.2) | kg·s⁻¹ |
| $\dot M_{j,i}$ | Face-flux part of the mass rate of component *j* in cell *i*, $-\sum_f s_{if} A_f F_{j,f}$ | kg·s⁻¹ |
| $\text{mob}_i$, $\text{mob}_f$ | Gas mobility $\kappa/\mu$ of cell *i*; interface fluxes use each half-cell mobility (§5.7.4, §9.6; overload M1) | m²·Pa⁻¹·s⁻¹ |
| $N_j$, $n_j$ | Species mass flux and molar flux relative to the skeleton | kg·m⁻²·s⁻¹; mol·m⁻²·s⁻¹ |
| $N_C$ (`NC`) | Number of chemical components (type parameter) | – |
| $N_R$ (`NR`) | Number of reactions (type parameter) | – |
| $N_{\text{state}}$ | Length of the state vector (§3.1.2, §15.1.1) | – |
| $N_{\text{node}}$, $N_{\text{cell}}$, $N_{\text{blk}}$ | Rows of a node block / of a cell block of the interleaved state on a line, and their sum, the stride (§15.1.1) | – |
| $\mathsf{N}_i$ | Normal matrix of the least-squares gradient of cell *i*, $\sum_k w_k \mathbf{d}_k \mathbf{d}_k^{T}$ (§13.6.4) | – |
| $n$ | Number of finite-volume cells | – |
| $n_c$ | Number of colors of the colored Jacobian (§15.7) | – |
| $n_{\text{layer}}$ | Number of cells across a layer of the equidistribution correction; its gain is $n_{\text{layer}}^{2}/\tau_{\text{eq}}$ (§11.7.2; overload N1) | – |
| $n_{\text{nodes}}$, $n_{\text{faces}}$ | Number of mesh nodes / faces (= n+1 in 1D) | – |
| $n_{i,j}$, $n_s$ | Reaction order of reaction *i* w.r.t. reactant *j*, an order on the conversion for $n \ne 1$ | – |
| $n_{\Gamma}$, $n_D$ | Exponents of the property-ratio film ($n_{\Gamma}$ on $\Gamma_{\infty}/\Gamma_g$, default 0.6) and of the pooling film's face fractions ($n_D$, default 1/2) (§12.3.3) | – |
| $P$ | Gas pressure in the porous medium | Pa |
| $P_{\text{ref}}$ | Reference pressure (= 101325) | Pa |
| $Q_{\text{rxn}}$ | Volumetric reaction heat source, $-\sum h_r r_r$ | W·m⁻³ |
| $Q_{\text{rad}}$ | Volumetric radiation absorption source | W·m⁻³ |
| $q$, $q_{\text{cond}}$, $q_f$ | (Conductive) heat flux; at face *f*, along its normal | W·m⁻² |
| $q_{\text{BC}}$ | Boundary-condition heat flux (positive into domain) | W·m⁻² |
| $q_{\text{inc}}$ | Incident radiative flux of a `RadiativeFluxBC` | W·m⁻² |
| $q'''$, $q_{\text{rad}}'''$, $q^{\text{rad}}_i$ | Volumetric radiation source (Beer–Lambert; $q^{\text{rad}}_i$ its average over cell *i*, §8.3.2) | W·m⁻³ |
| $R_b$ | Blowing ratio $h_{\text{eff}}/h$ of a convective condition's blowing model, 1 at zero species fluxes (§12.3.3; the code's $\chi$, overload B1) | – |
| $R_g$ ($R$, $R_{\text{gas}}$) | Universal gas constant (= 8.31446261815324) | J·mol⁻¹·K⁻¹ |
| $R_{s,i}$ | Reference field *s* of cell *i*, a differential unknown where some order differs from one (§6.4) | kg |
| $r$, $r_i$ | Reaction rate (per unit volume) of reaction *i* | kg·m⁻³·s⁻¹ |
| $r_{\text{hi}}$, $r_{\text{lo}}$ | Activation bounds of the equidistribution correction, fractions of a cell's initial thickness (defaults 0.3, 0.1) | – |
| $S$ | Sutherland constant (= 110.4 K, air) in the viscosity law | K |
| $S_j$, $S_{j,i}$ | Species source term for component *j* (in cell *i*), $\sum_r \nu_{r,j} r_r$ (= $\omega_{j,i}$); the rows carry the cell's reaction source $\omega_{j,i} \tilde V_i$ [kg·s⁻¹] | kg·m⁻³·s⁻¹ |
| $S_{\text{dry},i}$ | Dry-pyrolysis-gas source term in cell *i* (source of $G_i$) | kg·m⁻³·s⁻¹ |
| $S_{E,i}$ | Energy source of cell *i*: its reaction heat $Q_{\text{rxn},i} \tilde V_i$ plus the in-depth radiation it absorbs (§7.3.3, §13.3.4, §16.5.1) | W |
| $S_T$, $S_{T,i}$ | External temperature-rate source (manufactured solutions): $S_T$ on the temperature row, $C_i S_T$ on the enthalpy row (§7.3.3, §13.5.3) | K·s⁻¹ |
| $S^{\text{rad}}_i$ | In-depth radiation absorbed by cell *i*, $I_{\text{up}} A_{\text{up}} - I_{\text{down}} A_{\text{down}}$ (§8.3.3) | W |
| $s_{if}$ | Orientation of face *f* relative to cell *i* (written $s_{if}$ or $s_{i,f}$; $s_f$ in §6, where *i* is the reaction index and the cell carries no index): +1 when its normal points out of *i* (*i* its owner), −1 otherwise | – |
| $s_k$ | Diffusivity ratio $D_{k,\infty}/\alpha_{\text{th},\infty}$ of species *k* in the pooling film (§12.3.3) | – |
| $T$ | Absolute temperature | K |
| $T_s$ | Surface (boundary face) temperature | K |
| $T_{\text{cell}}$ | Temperature at the cell center adjacent to a boundary | K |
| $T_{\infty}$ | Ambient / surroundings temperature | K |
| $T_{\text{ref}}$ | Reference temperature for enthalpy datum (= 298.15) | K |
| $T_f$ | Face temperature of the enthalpy the gas fluxes carry: the distance-weighted interpolation $(d_n T_o + d_o T_n)/(d_o + d_n)$ at an interior face, the surface temperature for boundary outflow and specified incoming temperature for inflow (§13.6.2) | K |
| $\bar T_f$, $T_{\text{face}}$ | Face temperature; gas transport uses the distance-weighted interface state $T^*=(d_RT_L+d_LT_R)/(d_L+d_R)$, reducing to the arithmetic mean for equal half-cell distances (§9.6) | K |
| $\tilde T_f$ | Face temperature of the relative fluxes: the upwind van Albada reconstruction on the stencil of the mass flux it carries (§11.6.4, §13.6.2) | K |
| $T_{\text{trial}}$, $\tilde T_i$ | Trial temperature: $T$ smoothly clamped into the property bracket $[200, 3000]$ K with 1 K transitions | K |
| $T_m$ | Merged temperature of a deactivation receiver: the root of the extended caloric relation at $H_m$ and the merged masses (§14.3.3, §16.7) | K |
| $T_{\text{min},i}$, $T_{\text{max},i}$ | Lower / upper temperature gate for reaction *i* | K |
| $T^*$ | Reduced temperature $k_B T / \varepsilon_{AB}$ (collision integral) | – |
| $t$ | Time | s |
| $u$ | Flat DAE state vector | mixed |
| $du$ | Time derivative of state vector | mixed |
| $V$, $V_i$, $V_i(0)$ | Cell volume (of cell *i*; initial) | m³ |
| $\tilde V_i$ | Smoothly floored cell volume, $\operatorname{smooth\_floor}(V_i, 10^{-8} V_i(0))$, the volume the concentrations and the reaction sources use (code `V_eval`; §5.13.1, §13.10.2, §13.10.3) | m³ |
| $v_g$ | Bulk (Darcy, superficial) gas velocity at a face | m·s⁻¹ |
| $v_s$ | Solid (skeleton) velocity; at the nodes $v_s = w - c$ | m·s⁻¹ |
| $v_j$ | Volume fraction of component *j* | – |
| $v_j^{\text{solid}}$ | Condensed volume fraction, over the condensed (solid and liquid) components | – |
| $v_m$ | Matrix volume fraction, $\sum_j \xi_j/\rho_j$ over the condensed components (§5.1.3) | – |
| $v_{\text{out}}$ | Outward Darcy velocity a pressure condition drives at a boundary face (§12.5, §13.7.5) | m·s⁻¹ |
| $v_{\text{rel},f}$ | Velocity of the solid relative to face *f* along its normal, $(v_s - w)\cdot\hat n_f = -c_f \hat n_f$ (§11.6.1) | m·s⁻¹ |
| $W$ | The integrator's iteration matrix $J - M/(\gamma_{\text{RK}} \Delta t)$ | – |
| $w$, $w_n$ | Mesh node velocity (at node *n*), $w = v_s + c$; an algebraic unknown on a moving mesh (overload W1) | m·s⁻¹ |
| $w_k$ | Weight $\lvert\mathbf{d}_k\rvert^{-2}$ of neighbour *k* in the least-squares gradient (§13.6.4); a quadrature weight in §16.5.3 (overload W1) | m⁻²; – |
| $X_{\alpha}$ | Volume fraction (FDS/Gpyro V2 comparison notation; internally $v_j$) | – |
| $x$, $y$ | Transverse/lateral coordinates (multi-dimensional meshes and comparison text) | m |
| $x_a$ | Conversion of reactant $a$ against its reference, $\xi_a/\xi_{\text{ref}}$ (§6.4) | – |
| $x_J$ | Join of the depletion limiter's onset, in units of its threshold (`join`, default 1; §6.3.2) | – |
| $x_{\text{th}}$, $\delta_x$ | Conversion scales of the depletion limiter's onset and finite base (`conversion_threshold`, `conversion_base`; §6.3.2) | – |
| $Y_j$ | Mass fraction of component *j*, $\xi_j/\rho$ ($Y_k$ of species *k* at the face in the pooling film, §12.3.3) | – |
| $Z$ | Pre-exponential factor (Gpyro V2 notation; FDS uses $A_{\alpha\beta}$; internally $A$) | varies |
| $z$ | Axial / through-thickness coordinate (primary spatial variable) | m |
| $z_n$ | Position of mesh node *n* (a differential unknown on a moving mesh) | m |
| $z_i^c$ | Cell-center position of cell *i* | m |
| $\Delta z$, $\Delta z_i$, $\Delta z_i(0)$ | Cell thickness (of cell *i*; initial) | m |

---

## 2. Greek letters

| Symbol | Meaning | SI units |
|---|---|---|
| $\alpha$ ($\alpha_{\text{eff}}$) | Volumetric absorption (extinction) coefficient | m⁻¹ |
| $\alpha_j$ | Mass-basis absorption coefficient of component *j* | m²·kg⁻¹ |
| $\alpha_{\text{conv}}$ | Conversion / extent variable (Gpyro V2 kinetics comparison) | – |
| $\alpha_{\text{th},\infty}$ | Thermal diffusivity of the free stream (pooling film, §12.3.3; overload A5) | m²·s⁻¹ |
| $\alpha_M$, $\beta_M$ | Monitor coefficients of the equidistribution correction, `monitor = (α_M, β_M)` (defaults 1, 0) | – |
| $\beta$ | Weighting parameter of the WEIGHTED mixing rule, $k = \beta k_{\parallel} + (1-\beta)k_{\text{series}}$ | – |
| $\gamma_{\text{RK}}$ | Stage coefficient of the integrator in $W = J - M/(\gamma_{\text{RK}} \Delta t)$ | – |
| $\Gamma$ | Thermal property ratio $k/c_p$ of a gas (property-ratio and pooling films, §12.3.3; overload G1) | kg·m⁻¹·s⁻¹ |
| $\Gamma_i$ | Volume rate the node corrections sweep through the faces of cell *i*, $A_{i+1} c_{i+1} - A_i c_i$ on a line (§11.5.2; overload G1) | m³·s⁻¹ |
| $\Delta$ | Forward-difference / change operator (prefix) | varies |
| $\Delta_i$ | Thickness of cell *i* along the Beer–Lambert beam, the sum of its two floored centre-to-face distances ($\Delta z_i$ away from the floors; §8.3.3, §13.5.2) | m |
| $\delta_{ij}$ | Kronecker delta | – |
| $\delta_i$ | Positive-part mass scale of cell *i*, $\text{reltol} \cdot \varepsilon_m \cdot m^s_i(0)$ (its mass tolerance) | kg |
| $\delta_b$ | Scale of the softplus finite base of the depletion limiter (`base_scale`, default 10⁻⁸; §6.3.2) | kg·m⁻³ |
| $\delta_{\text{mov}}$, $\delta_{\text{BL}}$, $\delta_{\chi}$ | Configuration flags of the state layout: one on a moving mesh / under Beer–Lambert radiation / under a lateral-shrinkage law, zero otherwise (§3.1.2, §15.1.1) | – |
| $\varepsilon$, $\varepsilon_j$ | Emissivity (of component *j*) | – |
| $\varepsilon_{\text{eff}}$ | Effective (composition-weighted) emissivity | – |
| $\varepsilon_m$ | Deactivation threshold: fraction of a cell's reference below which it deactivates (`threshold`, default 1e-2, at an order other than one a tenth of the limiter's `conversion_threshold`; overload E1) | – |
| $\varepsilon_{\text{band}}$ | Near-threshold band of deactivation, $\varepsilon_m \le \varepsilon_{\text{band}} < 1$ (`band`) | – |
| $\varepsilon_{AB}$ | Lennard-Jones well depth (collision integral) | J |
| $\epsilon_{\xi}$, $\epsilon_T$, $\epsilon_G$ | Slope scales of the van Albada reconstruction for the concentrations, the temperature and the history values: the square of a tolerance over the upwind cell's reference thickness (§11.6.2; overload E1) | (field units)²·m⁻² |
| $\epsilon_{\text{LS}}$ | Relative regularization of the least-squares normal matrix, 10⁻⁸ (`GRADIENT_REGULARIZATION`), switched on smoothly by the fit's conditioning (§13.6.4; overload E1) | – |
| $\eta$ | Guard of the distance-weighted harmonic face means, 10⁻²⁰ in the property's units (`FACE_MEAN_GUARD`; §1.4.4, §5.12, §7.1.2); in §17 the relative finite-difference step (overload E2) | property units |
| $\theta$, $\theta_i$ | Volumetric strain rate (dilation) of cell *i*, $(1/V) dV/dt$ (bold $\theta$ in §17 is the parameter vector; overload A2) | s⁻¹ |
| $\theta_A$ | Column-global **radial** volumetric strain rate, $(1/A) dA/dt$ (an algebraic unknown under a lateral-shrinkage law) | s⁻¹ |
| $\theta_L$ | **Axial** volumetric strain rate ($\theta - \theta_A$) | s⁻¹ |
| $\theta_i^j$ | Stoichiometric coefficient (ThermaKin comparison notation; internally $\nu_{i,j}$) | – |
| $\kappa$ ($\kappa_j$, $\kappa_{\text{eff}}$) | Permeability (Darcy) of component *j* / of the cell's mixture | m² |
| $\kappa_s$ | Solid absorption coefficient (FDS; compare effective $\alpha_{\text{eff}}$, not mass-specific $\alpha_j$) | m⁻¹ |
| $\Lambda$ | Face buoyancy of the pooling film relative to the hot free stream, $[1 - (M_{\text{face}}/T_s)(T_{\infty}/M_{\infty})]/(1 - T_{\infty}/T_s)$ (§12.3.3; overload B1) | – |
| $\Lambda_c$, $\Lambda_w$ | Centre and width of the pooling factor's logistic switch in $\Lambda$ (`PoolingFilm` keywords `Λc`, `Λw`, defaults 0.35, 0.10; §12.3.3) | – |
| $\lambda$ ($\lambda_j$, $\lambda_{\text{face}}$, $\lambda_{\text{eff}}$) | Gas transfer (diffusion) coefficient | m²·s⁻¹ |
| $\lambda_{\tau}$ | Factor of the layer equidistribution time, $\tau_{\text{eq}} = \lambda_{\tau} L_{\text{layer}}/\vert \dot L_{\text{layer}}\vert $ (`ActivationTime`, default 1.8; overload A6) | – |
| $\ell$ | Characteristic particle/pore diameter (Carman–Kozeny); in §17 the scalar loss | m |
| $\mu$, $\mu_i$ | Gas dynamic (Sutherland) viscosity (of cell *i*) | Pa·s |
| $\mu_{\text{ref}}$ | Reference viscosity (= 1.716×10⁻⁵ at 273.15 K, air) | Pa·s |
| $\nu$ | Kinematic viscosity (Gpyro V2 Darcy comparison) | m²·s⁻¹ |
| $\nu_{i,j}$ | Stoichiometric coefficient (yield) of component *j* in reaction *i* (negative reactant, positive product) | – |
| $\xi_j$, $\xi_{j,i}$ | Mass concentration of component *j* (in cell *i*), $m_{j,i}/\tilde V_i$ — **constitutive species variable** | kg·m⁻³ |
| $\xi^{+}_{j,i}$ | Positive-part concentration $(m_{j,i})^{+}/\tilde V_i$ at which the property laws are evaluated | kg·m⁻³ |
| $\xi_{\text{ref},i}$ | Reference concentration of reaction *i*, its reference over the cell volume (§6.1, §6.4) | kg·m⁻³ |
| $\xi_{\text{th}}$ | Concentration scale of the depletion limiter (`DepletionLimiter` `threshold`, default 1 kg·m⁻³) | kg·m⁻³ |
| $\xi_{\text{tol},i}$ | Mass tolerance of cell *i* over its initial volume, the positivity scale of the relative-flux face states (§11.6.2, §11.6.3) | kg·m⁻³ |
| $\rho$ | Bulk density ($\sum \xi_j$ for the total inventory density) | kg·m⁻³ |
| $\rho_j$ | Bulk (pure-phase) density of component *j* | kg·m⁻³ |
| $\rho_{i,j}$ | Intrinsic / skeletal density of component *j* (skeletal volume fraction only) | kg·m⁻³ |
| $\rho_g$ | Gas-phase density (ideal-gas law) | kg·m⁻³ |
| $\rho_{s,j}$ ($\rho_{s,\alpha}$) | FDS component density: the condensed-component mass concentration $\xi_j = m_j/V$ | kg·m⁻³ |
| $\rho_s$ (FDS comparisons) | Total condensed-material density, $\sum_{S,L} \xi_j$ | kg·m⁻³ |
| $\rho c_p^{\text{eff}}$ | Effective (mixture) volumetric heat capacity, every component, $\sum_j \xi_j c_{p,j}$ (discretely $C_i/V_i$) | J·m⁻³·K⁻¹ |
| $\sigma$ | Stefan–Boltzmann constant (= 5.670374419×10⁻⁸) | W·m⁻²·K⁻⁴ |
| $\sigma_s$ | Scattering coefficient (not modeled) | m⁻¹ |
| $\sigma_{AB}$ | Lennard-Jones collision diameter | Å |
| $\tau_{\text{eq}}$ | Equidistribution time of a layer in the correction (keyword `τ`; overload T1) | s |
| $\tau_i$ | Optical thickness of cell *i*, $\alpha_i \Delta z_i$ (always subscripted or called optical; overload T1) | – |
| $\varphi$ | Porosity (gas void fraction), $1 - \sum_{S,L} \xi_j/\rho_{i,j}$ (U+03C6, LaTeX `\varphi`; overload A8) | – |
| $\phi_i$ | Multiplier of cell *i* in the potential closure (multi-dimensional meshes only; U+03D5, LaTeX `\phi`; overload A8) | m²·s⁻¹ |
| $\chi$, $\chi_i$ | **Pyrolysis progress** in cell *i*, $G_i/m_{\text{dry},i}$ (cumulative dry-gas mass released per unit initial dry-solid mass) | – |
| $\bar{\chi}$ | Column conversion, $\sum_i G_i / \sum_i m_{\text{dry},i}$ (the $m_{\text{dry}}$-weighted mean of $\chi_i$) | – |
| $\bar{\chi}_c$ | Column conversion smoothly clamped into $[-10^{-3}, 1 + 10^{-3}]$ (`LATERAL_CLAMP_WIDTH`), the argument of the lateral-shrinkage law (§10.4) | – |
| $\psi$, $\psi_0$ | Pooling factor of the pooling film, the share of the no-outflow exchange the face keeps, and its floor (`PoolingFilm` keyword `ψ0`, default 0.59; §12.3.3; overload B1) | – |
| $\psi_{\text{VA}}(a, b; \epsilon)$ | Smooth van Albada slope of the relative-flux reconstruction (§11.6.2) | (field units)·m⁻¹ |
| $\Phi$ | Generic cell field in the two-point gradient, where the face index $f$ shares the formula (§13.6.1); $\Phi_a(v)$ the reach of reactant *a* from node *v* of a reaction network, the mass of *a* one unit of *v*'s inventory can become (§6.5.5); in §17 the terminal part of a loss | varies |
| $\Omega(t)$, $\Omega_i(t)$ | Live region of the column (§16.2); region of cell *i* (§7.2, §11.2.1) | – |
| $\Omega_D$ | Diffusion collision integral (Neufeld) | – |
| $\omega_{j,i}$ | Net volumetric production rate of component *j* by reaction in cell *i*, $\sum_r \nu_{r,j} r_r$ (= $S_{j,i}$; the code's $\omega$); the cell's reaction source is $\omega_{j,i} \tilde V_i$ | kg·m⁻³·s⁻¹ |
| $\dot{\omega}'''$ | Volumetric reaction rate (Gpyro V2 comparison notation; internally $r$) | kg·m⁻³·s⁻¹ |

---

## 3. Subscripts

| Subscript | Meaning |
|---|---|
| $j$ | Component / species index (1…$N_C$) |
| $i$ | Cell index (1…n); also reaction index in kinetics context (see overload I1) |
| $n$ | Node index (1…n+1) in node quantities ($z_n$, $w_n$, $c_n$); bare $n$ is the cell count |
| $r$, $\alpha\beta$ | Reaction index; in §15, where stated, a row of the state vector ($f_r$, $s_r$, $J_{rc}$) |
| $k$ | Layer or cell index inside a sweep (closure and correction context); the receiver of a merge (§14, §16.7); a neighbour in the least-squares gradient (§13.6.4); a species of the pooling film (§12.3.3); a quadrature node (§16.5.3) |
| $q$ | Parameter index (§17); a node of a face in the multi-dimensional sweep vectors $s_{q,f}$ (§11.6.5) |
| $e$ | Exposed face of the Beer–Lambert beam ($I_e$, $A_e$, §8); a deactivation event in the ledger ($m^{\text{lost}}_e$, §16) |
| $b$ | Last face of the beam, at the back ($I_b$, $A_b$, §8) |
| $g$ | Gas phase / gaseous component |
| $s$ | Solid / condensed phase |
| $\ell$, $\text{liq}$ | Liquid phase |
| $L$, $R$ | Left / right of a face |
| $o$, $n$ (face context) | Owner / neighbour cell of a face |
| $f$ | Face-centered quantity |
| $\text{face}$ | Face-averaged property |
| $\text{eff}$ | Effective (mixture) property |
| $\text{cell}$ | Cell-centered quantity |
| $s$ (boundary) | Surface (boundary face) — e.g. $T_s$, $\xi_{j,s}$ |
| $\text{up}$, $\text{down}$ | Upstream / downstream face of a cell along the Beer–Lambert beam |
| $\text{inc}$ | Incident (radiation entering the exposed face) |
| $\infty$ | Ambient / surroundings; the free stream of a blowing model |
| $\text{ref}$ | Reference value |
| $\text{in}$ | Cell-entrance (radiation); entering the live region ($\dot H_{\text{in}}$) |
| $\text{out}$ | Leaving ($\dot m_{\text{out}}$, $v_{\text{out}}$) |
| $\text{live}$ | Restricted to the live cells or the live boundary faces (§16) |
| $\text{layer}$ | A layer of the equidistribution correction ($n_{\text{layer}}$, $L_{\text{layer}}$) |
| $\text{rel}$ | Relative to the moving mesh ($v_{\text{rel},f}$, $F^E_{\text{rel},f}$) |
| $\text{tol}$ | Tolerance scale ($\xi_{\text{tol},i}$) |
| $\text{floor}$, $\text{reg}$ | Floor / regularizing part of the caloric extension ($C_{\text{floor},i}$, $C_{\text{reg},i}$) |
| $\text{ext}$ | External source |
| $\text{min}$, $\text{max}$ | Lower / upper bound |
| $\text{dry}$ | Dry pyrolysis gas (the gases the solids give off) |
| $\text{init}$, $0$ | Initial value |
| $\text{top}$, $\text{bot}$ | Top (z=L) / bottom (z=0) boundary |
| $\text{rxn}$ | Reaction-source contribution |
| $\text{cond}$, $\text{conv}$, $\text{rad}$ | Conduction / convection / radiation contribution |
| $\parallel$, $\text{series}$ | Parallel / series mixing limit |
| $\text{total}$ | All-component total (e.g. $\rho_{\text{total}} = \sum_j \xi_j$) |

## 4. Superscripts and accents

| Mark | Meaning |
|---|---|
| $x^{\circ}$, $x^{o}$, $x_{\text{init}}$ | Value at previous time step / initial value |
| $x(0)$ | Value at the initial (reference) state, e.g. $V_i(0)$, $\Delta z_i(0)$, $m^s_i(0)$ |
| $x^{+}$, $x^{-}$ (radiation) | Forward / backward hemispheric intensity (two-flux, FDS comparison) |
| $x^{+}$ | One-sided positive part, e.g. $\xi^{+}$, $(m)^{+}$ (zero at and below zero, the identity above the scale $\delta_i$); never a state after an event, which carries a hat |
| $\hat x$ (hat) | The receiver's state after a deactivation merge ($\hat m_{j,k}$, $\hat V_k$, §14.3.3); $\hat n_f$ the unit normal of face *f* |
| $x^{\text{solid}}$ | Restricted to the condensed (solid and liquid) components |
| $x^{\text{diff}}$, $x^{\text{adv}}$ | Diffusive / advective flux component |
| $x^{E}$ | Energy (the face energy flux $F^E_f$) |
| $x^{\text{rad}}$ | In-depth radiation ($q^{\text{rad}}_i$, $S^{\text{rad}}_i$) |
| $x^{\text{lost}}$, $x^{\text{reg}}$ | What a deactivation event removes from the live region / the caloric regularization it writes to its record ($m^{\text{lost}}_e$, $H^{\text{lost}}_e$, $H^{\text{reg}}_e$) |
| $x^{\text{rel}}$ (on a defect) | Relative defect ($d_m^{\text{rel}}$, $d_H^{\text{rel}}$) |
| $x^{c}$ | Cell-center (e.g. $z_i^c$) |
| $\bar x$ (overbar) | Column / spatial average (e.g. $\bar{\chi}$, $\bar J_g$), a two-point face mean ($\bar T_f$), or Gpyro V2-style mixture average ($\bar{\rho}$); $\bar c_{p,\infty}$ the free-stream heat capacity averaged between $T_{\infty}$ and $T_s$ |
| $\tilde x$ (tilde) | Normalized, clamped or guarded trial-state quantity (e.g. $\tilde{\mathcal{M}}$, $\tilde T_i$ the trial temperature, $\tilde V_i$ the floored volume); on a face temperature, $\tilde T_f$, the relative-flux reconstruction |
| $\dot x$ (overdot) | Time rate / flux (e.g. $\dot m$, $\dot M$, $\dot{\omega}$) |
| $x'$, $x''$, $x'''$ | Per-length / per-area / per-volume (blowing-model fluxes $\dot m''$; FDS/Gpyro V2 comparison notation) |

## 5. Operators

| Operator | Meaning |
|---|---|
| $\partial/\partial t$ | Partial time derivative (Eulerian, fixed point) |
| $D_s/Dt$ | Solid-following time derivative, $\partial/\partial t + v_s\cdot\nabla$ (the mesh-frame rate where the mesh moves with the solid) |
| $d/dt$ | Total time derivative; of a cell inventory, following the cell |
| $\nabla$, $\nabla\cdot$ | Gradient / divergence (1D: $\partial/\partial z$, $(1/V_i) \sum_f s_{if} A_f F_f$) |
| $\nabla^{2}$ | Laplacian / second derivative |
| $\sum$, $\int$ | Sum / integral |
| $\operatorname{clamp}(x,a,b)$ | Clamp $x$ to $[a,b]$ |
| `smooth_clamp`, `smooth_floor`, $(\cdot)^{+}$, $\mathcal{P}(x, \delta)$ | The C¹ guards of trial states (§13.10.3): clamp with cubic joins, floor, one-sided positive part (written $\mathcal{P}(x, \delta)$ with its scale $\delta$ in §11.6.3) |
| $\mathcal{S}(x)$ | Smoothstep $x^{2}(3 - 2x)$ clamped to $[0, 1]$ (the correction's activation, §11.7.3) |
| $\tanh$ | Hyperbolic tangent (smooth temperature gates; the depletion onset above its join) |
| $\Vert \cdot\Vert _2$, $\Vert \cdot\Vert _W$ | Euclidean norm / norm weighted by the matrix $W$ of a least-squares loss (§17, §18) |
| $\operatorname{diag}(v)$ | Diagonal matrix from vector $v$ |
| $J$, $M$, $W$ | Jacobian / mass matrix / iteration matrix (§15) |

---

## 6. Physical constants (fixed values)

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Universal gas constant | $R_g$ | 8.31446261815324 | J·mol⁻¹·K⁻¹ |
| Stefan–Boltzmann constant | $\sigma$ | 5.670374419×10⁻⁸ | W·m⁻²·K⁻⁴ |
| Reference pressure | $P_{\text{ref}}$ | 101325 | Pa |
| Reference temperature | $T_{\text{ref}}$ | 298.15 | K |
| Reference molar mass (air) | $M_{\text{ref}}$ | 0.029 | kg·mol⁻¹ |
| Reference viscosity (air, 273.15 K) | $\mu_{\text{ref}}$ | 1.716×10⁻⁵ | Pa·s |
| Sutherland constant (air) | $S$ | 110.4 | K |
| Kozeny constant (spheres) | — | 180 | – |

### 6.1 Numerical constants and defaults

| Quantity | Symbol | Value | Units | Source |
|---|---|---|---|---|
| Property bracket | $[T_{\text{min}}, T_{\text{max}}]$ | [200, 3000] | K | `T_PROPERTY_MIN`, `T_PROPERTY_MAX` |
| Width of the trial clamp's joins | — | 1 | K | `TRIAL_TEMPERATURE_WIDTH` |
| Floor of the volumetric heat capacity | $\rho c_{p,\text{floor}}$ | 1.0×10³ | J·m⁻³·K⁻¹ | `RHO_CP_FLOOR` |
| Porosity floor (smooth) | — | 10⁻¹² | – | `POROSITY_FLOOR` |
| Viscosity floor (smooth) | — | 10⁻¹² | Pa·s | `MU_FLOOR` |
| Guard of the harmonic face means | $\eta$ | 10⁻²⁰ | property units | `FACE_MEAN_GUARD` |
| Regularization of the least-squares gradient | $\epsilon_{\text{LS}}$ | 10⁻⁸ | – | `GRADIENT_REGULARIZATION` |
| Conditioning of the gradient fit at and below which the regularization is full | — | 10⁻⁶ | – | `GRADIENT_SPAN` |
| Conditioning of the gradient fit at and above which the regularization vanishes | — | 10⁻⁴ | – | `GRADIENT_SUPPORT` |
| Relative normal conductivity of an insulating neighbour in the gradient fit | — | 10⁻¹² | – | `INSULATING_NEIGHBOUR` |
| Half-width of the clamp band of $\bar{\chi}_c$ | — | 10⁻³ | – | `LATERAL_CLAMP_WIDTH` |
| Deactivation threshold | $\varepsilon_m$ | 10⁻² | – | `solve(...; threshold)` |
| Near-threshold band under a moving correction | $\varepsilon_{\text{band}}$ | 0.1 | – | `CORRECTION_BAND` |
| Activation bounds | $(r_{\text{hi}}, r_{\text{lo}})$ | (0.3, 0.1) | – | `EquidistributionCorrection` |
| Layer-time factor | $\lambda_{\tau}$ | 1.8 | – | `ACTIVATION_TIME_FACTOR` |
| Fallback equidistribution time | $\tau_{\text{eq}}$ | 720 | s | `EQUIDISTRIBUTION_TIME` |
| Depletion-limiter scale | $\xi_{\text{th}}$ | 1.0 | kg·m⁻³ | `DepletionLimiter()` |
| Join of the depletion onset | $x_J$ | 1.0 | – (units of $\xi_{\text{th}}$) | `DepletionLimiter()` |
| Scale of the finite power base | $\delta_b$ | 10⁻⁸ | kg·m⁻³ | `DepletionLimiter()` |

---

## 7. Overloaded / reserved symbols — resolutions (BINDING)

These symbols are overloaded somewhere in the source subsystems. The resolution
below is authoritative; writers must use the disambiguating subscript/context.

- **χ (chi) — A1 [pyrolysis progress vs. cross-section ratio].** The variable
  cross section defines $\chi_i = G_i/m_{\text{dry},i}$ as the **pyrolysis progress**
  (dimensionless, cumulative dry-gas mass released per unit initial dry-solid
  mass) and $\bar{\chi} = \sum G_i/\sum m_{\text{dry},i}$ as the column conversion; it is easily
  misread as the cross-section area ratio. **Resolution:** $\chi$ is *always*
  pyrolysis progress; the cross-section ratio is $A/A_0 = \text{law}(\bar\chi)$ and is never
  written as $\chi$. No time-derivative operator carries a $\chi$ subscript. The
  blowing ratio, which the code calls `χ`, is written $R_b$ (B1).

- **θ (theta) — A2 [strain rate vs. stoichiometry vs. parameters].** Internally
  $\theta$ is the **volumetric strain rate** (s⁻¹), with $\theta_A$ radial and $\theta_L$
  axial. ThermaKin uses $\theta_i^j$ for stoichiometric coefficients.
  **Resolution:** use $\nu_{i,j}$ for stoichiometry everywhere; reserve $\theta$ for
  strain rate. $\theta_i^j$ appears only in the ThermaKin-comparison passage and is
  immediately mapped to $\nu_{i,j}$. The one exception is §17 (sensitivity), where
  the bold $\theta \in \mathbb{R}^p$ is the vector of parameters and $\theta_q$ its *q*-th entry;
  the strain rate does not appear in that chapter.

- **A — A3 [area vs. pre-exponential factor].** $A$ is the **cross-sectional
  area** and $A_f$ the area of face $f$; the Arrhenius pre-exponential is $A_i$ /
  $A_{\alpha\beta}$ (always subscripted by reaction). **Resolution:** bare $A$ = area in
  physics chapters; $A_f$ = face area; $A_i$ = pre-exponential (kinetics). The
  Gpyro V2 pre-exponential $Z$ and FDS $A_{\alpha\beta}$ are compared with $A_i$; copying values requires matching the rate basis (§6.9). In the
  mesh-motion rows on a line (§10, §11, §15) a face is named by its node, so
  $A_i$ and $A_{i+1}$ there are the areas of the faces at nodes $i$ and $i + 1$,
  the two faces of cell $i$; a reaction subscript never appears in those rows.

- **κ vs. K — A4 [permeability].** Use **$\kappa$** for permeability throughout.
  $K$ (Gpyro V2 permeability) appears only in comparison text and is mapped to $\kappa$.

- **α — A5 [absorption vs. conversion vs. absorptivity].** $\alpha$ (or $\alpha_{\text{eff}}$) =
  **volumetric absorption coefficient** (m⁻¹); $\alpha_j$ = mass-basis absorption
  coefficient (m²·kg⁻¹). Conversion (Gpyro V2) is $\alpha_{\text{conv}}$; FDS $\kappa_s$ has the units of $\alpha_{\text{eff}}$, not $\alpha_j$; matching attenuation also requires the angular closure (§8.7.3).
  The absorptivity of a surface for incident radiation is **$a$** ($a_j$,
  $a_{\text{eff}}$, $a_s$), never $\alpha$; it equals $\varepsilon_{\text{eff}}$ under Kirchhoff's law unless a
  component sets its own `surface_absorptivity`. The correction's monitor
  coefficients are $\alpha_M$, $\beta_M$, always subscripted, and the free stream's
  thermal diffusivity in the pooling film is $\alpha_{\text{th},\infty}$.

- **λ — A6 [gas transfer vs. anything else].** $\lambda$ is reserved for the **gas
  transfer (diffusion) coefficient** (m²·s⁻¹). Do not use $\lambda$ for thermal
  conductivity (that is $k$) or eigenvalues. Two exceptions carry a code keyword
  $\lambda$: the factor of `ActivationTime` (default 1.8), written **$\lambda_{\tau}$** wherever a
  formula needs it (prose may call it "the layer-time factor λ = 1.8" when
  `ActivationTime` is named in the same sentence); and the film constant of
  `CouetteFilm` (default 1/2), which appears only inside $B_{\lambda} = 2\lambda \dot m'' c_{p,g}/h$
  in §12.3.3 and the User Guide's blowing passage, always with `CouetteFilm`
  named. The layer-time factor is never written $\Lambda$; $\Lambda$ is the pooling film's
  face buoyancy (B1).

- **ρ subscripts — A7 [bulk vs. intrinsic].** $\rho_j$ = bulk density (mixing,
  mixture density); $\rho_{i,j}$ = intrinsic/skeletal density (skeletal volume
  fraction only). Never write a bare $\rho_i$ for skeletal density without the
  explicit "intrinsic" qualifier; prefer $\rho_{i,j}$.

- **h — H1 [heat of reaction; sign; sensible enthalpy].** $h$ (kinetics) and
  $h_r$ = heat of reaction per kg of first reactant, with **h > 0 endothermic**
  internally; ThermaKin publishes $h > 0 = \text{exothermic}$ (Gpyro V2 and FDS do not). The conversion to the
  volumetric source carries the minus sign: $Q_{\text{rxn}} = -\sum h_r r_r$. With a
  *component* subscript, $h_j(T)$ is the component's **sensible enthalpy**
  $\int_{T_{\text{ref}}}^T c_{p,j} dT$; the two never share a subscript ($r$ reaction, $j$
  component). The sensible enthalpy on the datum is always $h_j(T)$, never
  $\Delta h_j(T)$; $\Delta h$ is reserved for a difference between two stated temperatures
  ($\Delta h_g(T_s, T_{\infty})$) and for the effective heat of combustion $\Delta h_{c,\text{eff}}$. The
  convective coefficient is $h_{\text{conv}}$, the mass coefficient $h_m$, the pressure coefficient
  $h_P$, all subscripted (inside §12.3.3 the convective coefficient of the
  blowing formulas is the bare $h$, and $h_{\text{eff}} = R_b h$). A cell thickness is
  $\Delta z_i$, never $h_i$, even where the code names it $h$. In §17 $h_q$ is the
  finite-difference step of parameter $q$.

- **H — H2 [Heaviside vs. inventory].** $H(\cdot)$ with an argument is the
  Heaviside step (realized as a tanh ramp), and $H_{\text{min},i}(T)$, $H_{\text{max},i}(T)$
  are the smooth temperature gates of reaction $i$ (§6.2.2); $H_i$ with a cell
  subscript is the sensible-enthalpy inventory of cell $i$, and $H_{\text{live}}$, $H_m$,
  $H^{\text{lost}}_e$, $H^{\text{reg}}_e$ are the ledger's and the deactivation event's
  sensible enthalpies (§16).

- **i — I1 [cell index vs. reaction index].** $i$ is the **cell index** in
  discretization/geometry/mesh-motion chapters; in the kinetics chapter $i$ is
  the **reaction index** (with $j$ the component index). Each chapter states its
  convention in its first paragraph. Where both appear together, use $i$ = cell,
  $r$ = reaction.

- **φ — A8 [porosity only].** $\varphi$ (U+03C6, LaTeX `\varphi`) is **porosity**; do
  not reuse it for view factor (that is $F$), a potential or a generic field
  (that is $f$, or $q$ in a reconstruction). $\varphi$ denotes the intrinsic-density
  form $1 - \sum_{S,L} \xi_j/\rho_{i,j}$, one minus the skeletal volume fraction; it is
  evaluated unclamped where admissibility is tested and floored smoothly at
  10⁻¹² where the pressure closure divides by it. The all-phase gas volume
  fraction $\sum v_j^{\text{gas}}$ is a different quantity and is not written $\varphi$. The
  potential closure's per-cell multiplier (code `φ`) is written $\phi_i$ (U+03D5,
  LaTeX `\phi`); since the two glyphs are easily confused, $\phi$ appears only with
  its cell subscript and only in the multi-dimensional closure. The blowing
  models' film parameter, which the code calls `φ`, is written $B_{\lambda}$ or $B_{\Gamma}$
  (B1).

- **S — S1 [source vs. Sutherland constant].** $S$ with a subscript
  ($S_j$, $S_{\text{dry}}$, $S_{E,i}$, $S_T$, $S^{\text{rad}}_i$) is a source term; bare $S$ is
  the Sutherland constant in the viscosity law only. The lowercase $s_{if}$ is the
  face orientation sign; $s_r$ (row subscript, §15.5.4) is the equation scale of
  a row and $s_k$ (§12.3.3) the pooling film's diffusivity ratio. In §17,
  $S_{kq} = \partial y_k/\partial\theta_q$ is the sensitivity matrix and $S_q^1$, $S_q^T$ the
  first-order and total Sobol indices; no source term appears there.

- **ε — E1 [emissivity vs. thresholds vs. numerical scales].** $\varepsilon$, $\varepsilon_j$,
  $\varepsilon_{\text{eff}}$ are emissivities; $\varepsilon_m$ (deactivation threshold) and $\varepsilon_{\text{band}}$ (band)
  are fractions of a cell's reference condensed mass; $\varepsilon_{AB}$ is the
  Lennard-Jones well depth. The subscript always disambiguates. The numerical
  scales of the van Albada slope ($\epsilon_{\xi}$, $\epsilon_T$, $\epsilon_G$, and the $\epsilon$ argument of
  $\psi_{\text{VA}}$) and of the least-squares regularization ($\epsilon_{\text{LS}}$) use the lunate $\epsilon$
  (LaTeX `\epsilon`); the physical quantities and the deactivation fractions use
  $\varepsilon$ (LaTeX `\varepsilon`).

- **η — E2 [guards vs. finite-difference step].** $\eta$ is the guard added to the
  denominator of a guarded harmonic face mean, in the units of the averaged
  property (10⁻²⁰), for thermal conduction. In
  §17 $\eta$ is the relative finite-difference step and $\eta_{\text{abs}}$ its absolute
  floor; no face mean appears there.

- **τ — T1 [optical thickness vs. equidistribution time].** $\tau_i$ with a cell
  subscript, or "optical thickness" in words, is $\alpha_i \Delta z_i$; $\tau_{\text{eq}}$ in the
  mesh-correction and burnout chapters (§11, §14) is a layer's equidistribution
  time [s]. A time variable of integration is written $t$, never $\tau$.

- **M — M1 [molar mass vs. mass matrix vs. history vs. monitor].** $M$, $M_j$
  (and $M_{\text{face}}$, $M_{\infty}$) is molar mass; $M$ in the linear-algebra context of §15
  is the integrator's mass matrix (sans-serif $\mathsf{M}$ in §3 and §15). The second
  member of the conversion-history pair, which the code calls `M`, is written
  $m_{\text{dry},i}$, its rate $\dot m_{\text{dry},i}$; the correction's monitor is the
  calligraphic $\mathcal{M}$ ($\tilde{\mathcal{M}}_i$ normalized). The gas mobility $\kappa/\mu$ is written
  $\text{mob}_i$ ($\text{mob}_f$ at a face), never $m$, which is a mass.

- **c — C1 [heat capacity vs. correction].** The specific heat capacity is
  always $c_p$ ($c_{p,j}$); bare $c$ (or $c_n$ at node $n$) is the node correction
  $w - v_s$.

- **a — A9 [absorptivity vs. activation].** $a_j$, $a_{\text{eff}}$, $a_s$ (component,
  effective or surface subscript) are surface absorptivities; $a_i$, $a_n$ (cell
  or node subscript, correction context only) are the correction's activation
  switches. A row's absolute tolerance is written $\text{atol}_r$ (§15.5.4), never
  $a_r$.

- **N — N1 [species fluxes vs. counts].** $N_j$ is species mass flux and
  $n_j$ is species molar flux. Total pressure is $P = \sum_j p_j$.
  Every count carries a subscript: $N_C$, $N_R$,
  $N_{\text{state}}$, $N_{\text{node}}$, $N_{\text{cell}}$, $N_{\text{blk}}$. The number of cells across a layer of
  the equidistribution correction is $n_{\text{layer}}$, in the gain $n_{\text{layer}}^{2}/\tau_{\text{eq}}$.

- **Blowing — B1 [blowing ratio, film parameters, pooling].** The code's
  blowing ratio $\chi$ is written **$R_b$** ($h_{\text{eff}}/h$) and its film parameter $\varphi$
  is written **$B_{\lambda}$** (`CouetteFilm`) or **$B_{\Gamma}$** (property-ratio and pooling
  films), because $\chi$ is the pyrolysis progress (A1) and $\varphi$ the porosity (A8);
  $F_b(x) = x/(e^x - 1)$ is the film factor, $\psi$ the pooling factor ($\psi_0$ its
  floor) and $\Lambda$ the face buoyancy ($\Lambda_c$, $\Lambda_w$ the switch constants). Code
  keyword names (`λ`, `n_Γ`, `n_D`, `ψ0`, `Λc`, `Λw`) are kept as they are in
  code wherever a keyword is meant.

- **d — D1 [distances vs. defects vs. least-squares vectors].** $d$ with a face
  or cell subscript ($d_L$, $d_o$, $d_{U,f}$, $d_{LR}$) is a distance;
  $d_m$, $d_H$ (and their relative forms $d^{\text{rel}}$) are the conservation audit's
  mass and energy defects (§16); the bold $\mathbf{d}_k$ is the centre-to-centre vector
  of the least-squares gradient (§13.6.4).

- **F — F1 [view factor vs. fluxes vs. film factor vs. deformation].** Bare $F$
  is the radiative view factor; $F_f$, $F_{j,f}$, $F^E_f$ (face subscript) are
  face fluxes; $F_b(x)$ is the blowing models' film factor; $F_i$ (cell
  subscript, meshes of more dimensions) is the deformation gradient.

- **Γ — G1 [thermal property ratio vs. sweep].** $\Gamma = k/c_p$ (with $\infty$ or $g$,
  free stream or leaving gas) belongs to the blowing models of §12.3.3; $\Gamma_i$
  with a cell subscript is the volume rate the node corrections sweep (§11.5.2).

- **w — W1 [node velocity vs. weights].** $w$, $w_n$ is the mesh node velocity;
  $w_k$ (neighbour or quadrature-node subscript) is a weight of the
  least-squares gradient (§13.6.4) or of a quadrature rule (§16.5.3).
