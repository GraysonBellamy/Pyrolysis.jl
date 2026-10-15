# 13. Finite-Volume Spatial Discretization

This chapter derives the cell-centered finite-volume (FVM) spatial discretization that converts the continuous coupled PDE system of §3 (Governing Equations) into the spatially discrete semi-discrete system of ordinary differential / differential-algebraic equations integrated in §15 (Temporal Integration & Structured Jacobian). It documents the 1D mesh layout (nodes, cells, faces), the geometric quantities that enter every flux, the divergence operator and its sign conventions, the conductive / gas-diffusive / Darcy–Fick face-value formulas, the TVD flux limiters and gradient-reconstruction machinery, the moving-mesh (ALE) advection discretization, mesh generation (uniform and geometric-stretch), and the block-major state-vector packing with allocation-free geometry-from-state. Every equation, sign, factor, guard, and limiter below has been confirmed against the implementation; where the code deviates from the textbook form (or from the design-document prose), the code is reported as ground truth and the deviation is called out.

The unifying design principle is that **all transport is expressed as fluxes evaluated at cell faces**, and the time derivative of every conserved or quasi-conserved cell quantity is the discrete divergence of those face fluxes plus volumetric sources. This guarantees discrete conservation up to the boundary-flux treatment (see §16, Conservation Diagnostics) and produces a sparse, banded coupling structure that the structured Jacobian of §15 exploits.

Throughout this chapter the index `i` is the **cell index** (overload resolution I1: `i` = cell in discretization; the reaction index is `r`), `j` is the **component index** (`1…N_C`), and `f` is the **face index**. The spatial coordinate is `z` with `z = 0` at the bottom/substrate (boundary id 2, tag `:bottom`) and `z = L` at the top/exposed surface (boundary id 1, tag `:top`); heat enters from the top and the material recedes downward. The three sign conventions of §2 (Nomenclature) are in force: positive flux is transport in `+z`; the divergence `(F_R − F_L)A/V` is positive for net outflow; and the heat-of-reaction convention `h > 0` endothermic, `Q_rxn = −Σ_r h_r r_r`.

---

## 13.1 Cell-Centered Finite-Volume Layout in 1D

### 13.1.1 Topology: nodes, cells, faces

Pyrolysis.jl stores its primary unknowns — temperature `T` and the mass concentrations `ξ_j` — at **cell centers**, not at faces or nodes. This is the classic cell-centered (collocated) finite-volume arrangement. The 1D mesh consists of three topological entity types whose counts are rigidly linked:

$$
n_{\text{nodes}} = n_{\text{faces}} = n + 1, \qquad n \equiv n_{\text{cells}}.
$$

- A **node** (`Node1D`) is a mesh vertex at position `z_i`. Node 1 sits at `z = 0` (bottom, `boundary_id = 2`), node `n+1` sits at `z = L` (top, `boundary_id = 1`); the rest are interior (`boundary_id = 0`). In ALE mode the node positions become part of the ODE state (§13.10).
- A **cell** (`Cell1D`) is the finite volume. Cell `i` spans the two adjacent nodes `(i, i+1)` and is bounded by the two faces `(i, i+1)`. Cells carry the unknowns through their `CellState{NC}` (`T`, the `NTuple{NC}` of `ξ_j`, and a diagnostic `volume`). Cells are immutable in connectivity; refinement/coarsening rebuilds the mesh (§14, Adaptive Mesh Refinement).
- A **face** (`Face1D`) is, in 1D, a single point co-located with a node. Face `f` carries its `cell_ids = (left, right)` with `0` marking a domain boundary: the bottom face (id 1, at `z = 0`) has `cell_ids = (0, 1)`; the top face (id `n+1`, at `z = L`) has `cell_ids = (n, 0)`; interior face `i` has `cell_ids = (i−1, i)`.

"Left" always denotes the lower-`z` (toward substrate) side of a face and "Right" the higher-`z` (toward surface) side. Cell indices increase from bottom to top, so node `i` is the left endpoint of cell `i` and the right endpoint of cell `i−1`; face `i` co-locates with node `i`.

> **Implementation note.** The data structures are defined in `Node1D` (src/Geometry/nodes.jl), `Cell1D`/`CellState` (src/Geometry/cells.jl), `Face1D` (src/Geometry/faces.jl), and the container `Mesh1D{NC}` (src/Geometry/mesh.jl). The immutable connectivity view `MeshTopology` and the mutable per-step `MeshGeometry` are kept separate (src/Geometry/topology.jl, src/Geometry/mesh_geometry.jl) so that the same connectivity supports many geometry snapshots — this is what makes pure (mesh-non-mutating) residual and Jacobian evaluation possible.

### 13.1.2 Boundary identifier convention

The internal representation uses integer boundary ids; the user-facing API uses symbols. The bridge is

$$
\texttt{tag\_to\_boundary\_id}: \;\; :\texttt{top} \mapsto 1, \quad :\texttt{bottom} \mapsto 2, \quad \text{otherwise} \mapsto 0,
$$

with the inverse `boundary_id_to_tag` (src/Geometry/topology.jl). The boundary-face map `mesh.boundary_faces` stores `1 => [n+1]` (top) and `2 => [1]` (bottom). This convention propagates into the boundary-temperature solver, the mass/pressure/mesh boundary conditions of §12, and the diagnostics of §16.

### 13.1.3 What a "cell-centered" scheme means for accuracy

Because unknowns live at cell centers but fluxes are needed at faces, every flux requires an *interpolation/reconstruction* from cell values to faces. For diffusion this is a face-property average plus a two-point gradient (§13.6, §13.7); for advection it is an upwind or limited reconstruction (§13.4, §13.5). On a uniform mesh the central two-point gradient at a face is second-order; the cell-centered gradient (§13.6) is second-order in the interior and degrades to first-order one-sided at boundaries. The harmonic-mean face conductivity (§13.7) preserves flux continuity exactly across material interfaces, which is the physically correct behavior for a layered/heterogeneous medium and the main reason the cell-centered scheme is preferred over a vertex-centered one here.

---

## 13.2 Geometric Quantities

All discrete operators are built from a small set of geometric quantities derived from node positions and the (possibly time-varying) cross-sectional area `A(t)`.

### 13.2.1 Cell volume and cell center

For cell `i` spanning nodes `{z_i, z_{i+1}}`:

$$
V_i = |z_{i+1} - z_i|\, A, \qquad z_i^c = \tfrac{1}{2}\,(z_i + z_{i+1}).
$$

The absolute value guards against transiently inverted nodes during ALE motion; the cross-sectional area `A` is the **current** column area (it equals `A_0` unless a lateral-shrinkage law is active — see §10, Volume Change). All faces in 1D carry the same area `A_f = A`.

> **Implementation note.** `update_geometry!` (src/Geometry/geometric_quantities.jl) writes `cell_centers[i]`, `cell_volumes[i] = abs(z_right − z_left) * A`, and `face_areas[f] = A`. The mesh-resident path is used between integrator steps; the state-vector path `compute_geometry_from_state!` is used *during* residual/Jacobian assembly (§13.10).

### 13.2.2 Interior face distance

Diffusive fluxes use the **center-to-center distance** across a face, not the cell width. For an interior face `f` with neighbors `L` and `R`:

$$
d_{LR} = |z_R^c - z_L^c|.
$$

This is the denominator of every two-point gradient in §13.7 and is the correct choice for second-order accuracy on non-uniform meshes (using cell width instead would introduce a zeroth-order multiplicative bias).

### 13.2.3 Ghost-cell boundary distance

At a boundary face the "other side" cell is absent. The code supplies a reflected ghost cell whose center sits at the mirror position `z^c_{\text{ghost}} = 2 z_f − z^c_{\text{adj}}`, so the effective boundary-face distance is

$$
d_f^{\text{bnd}} = 2\,|z_f - z^c_{\text{adj}}|.
$$

The factor of 2 encodes a symmetric (zero-gradient / reflective) extrapolation across the domain edge, turning a Dirichlet value at the face into an equivalent Neumann stencil. This is implemented in `face_distance` and `compute_face_distances_from_state!` (src/Geometry/geometric_quantities.jl), with `left_id == 0` for the bottom face and `right_id == 0` for the top.

### 13.2.4 The ALE advection boundary distance

The ALE advection operator (§13.8) needs an *upwind* distance, computed by `_upwind_gradient_distance` (src/Discretization/ale_fluxes.jl). For an interior upwind neighbor it returns the center-to-center distance, consistent with §13.2.2. For a boundary (`neighbor_id == 0`) it returns

$$
d^{\text{bnd}}_{\text{upwind}} = \frac{V_i}{A} = \Delta z_i,
$$

i.e. the **full cell width**. This is precisely `face_distance` (§13.2.3) evaluated at a boundary: the center-to-mirror-ghost spacing, twice the center-to-face distance — the same footing as the interior center-to-center distance. At a boundary the value is in any case inert: the zero-gradient ghost makes the first-order gradient identically zero in that direction, so no residual or Jacobian entry depends on it (§11.6.3).

> **Implementation note.** The ALE state-vector geometry (`GeometryCache`, src/Geometry/geometric_quantities.jl) caches `volumes`, `centers`, `face_positions`, and `face_distances`; the boundary entries of `face_distances` use the same factor-2 ghost-cell form of §13.2.3.

---

## 13.3 The Divergence Operator

### 13.3.1 Discrete divergence of a scalar flux field

The control-volume integral of a divergence reduces, via the divergence theorem, to a sum of outward face fluxes. In 1D with a left face (outward normal `−ẑ`) and a right face (outward normal `+ẑ`):

$$
(\nabla\!\cdot\!\mathbf F)_i \;=\; \frac{1}{V_i}\oint_{\partial V_i}\mathbf F\cdot \hat{\mathbf n}\,\mathrm dS
\;\xrightarrow{\text{1D}}\;
\boxed{\;\text{div}[i] = \frac{(F_{f_R} - F_{f_L})\,A}{V_i}\;}
$$

where `F_{f_L}`, `F_{f_R}` are the face fluxes at the left and right faces of cell `i`, `A` is the (common) face area, and `V_i` is the cell volume. Because positive flux is `+z` transport, the right-face flux is outflow `(+)` and the left-face flux is inflow `(−)`; a positive divergence therefore means **net efflux from the cell**, which in the energy equation enters with a leading minus sign so that net efflux cools the cell.

> **Implementation note.** `divergence_1d!` (src/Discretization/finite_volume.jl) reads `V_i = state.volume[i]` and `A = state.A_face[f_left]` and computes `(fluxes[f_right] − fluxes[f_left]) * A / V_i`. Both `V_i` and `A` carry the live scalar type, so in ALE/variable-cross-section mode the divergence is differentiable through the geometry (essential for the Jacobian and for the sensitivities of §17).

### 13.3.2 Multi-component (species) divergence

For the gas mass fluxes, the same operator is applied component-wise to the `NC × n_faces` flux matrix `J`, producing the `NC × n_cells` divergence matrix:

$$
\text{div}_j[i] = \frac{(J_{j,f_R} - J_{j,f_L})\,A}{V_i}, \qquad j = 1,\dots,N_C.
$$

> **Implementation note.** `divergence_species_1d!` (src/Discretization/finite_volume.jl) loops over cells and components with the identical formula. Solid/liquid components carry zero face flux (they are mesh-bound, §13.7.4), so their divergence is identically zero and the species residual reduces to reactions only (§13.3.4).

### 13.3.3 Discrete conservation property

Summing the scalar divergence over all cells telescopes:

$$
\sum_{i=1}^{n} V_i\,(\nabla\!\cdot\!\mathbf F)_i = A\sum_{i=1}^{n}(F_{f_R(i)} - F_{f_L(i)}) = A\,(F_{\text{top face}} - F_{\text{bottom face}}),
$$

so the integral of the divergence equals the net boundary flux exactly — the discrete analog of Gauss's theorem. This identity is the foundation of the conservation diagnostics in §16 and is exercised by the unit tests (constant flux → zero divergence; linear flux → constant divergence; integral closure to boundary fluxes; test/discretization/test_discretization.jl, test/discretization/test_conservation.jl).

### 13.3.4 Sign assembly in the final residual

The cell-wise temperature derivative assembled at the end of `compute_rhs!` is

$$
\frac{dT_i}{dt} = \frac{1}{\rho c_{p,\text{eff},i}}\Bigl[\,-(\nabla\!\cdot\!\mathbf q)_i + S_{\text{conv},i} + S_{\text{gen},i} + Q_{\text{rad},i} + Q_{\text{rxn},i}\,\Bigr],
$$

and the species derivative (per cell, per component) is

$$
\frac{d\xi_{j,i}}{dt} =
\begin{cases}
-(\nabla\!\cdot\!\mathbf J_j)_i + S_{j,\text{rxn},i}, & \text{component } j \text{ is gas},\\[4pt]
S_{j,\text{rxn},i}, & \text{component } j \text{ is solid/liquid}.
\end{cases}
$$

The leading `1/ρc_{p,\text{eff}}` uses the mixture volumetric heat capacity (every component; see §7, Heat Conduction & Energy Assembly, and §13.7.5). The gas/solid dispatch is exactly the `is_gas`-branch in `apply_species_residual_unrolled!` (src/Materials/materials.jl): gas → `-div_J + dξ_rxn`, solid/liquid → `dξ_rxn`. The full eight-step assembly order is given in §7.8; this chapter is concerned only with the spatial operators that produce `div_q`, `div_J`, and the face fluxes feeding them.

---

## 13.4 Flux Limiters (TVD Schemes)

Advective transport with a centered reconstruction is second-order but oscillatory near steep fronts; a pure upwind reconstruction is monotone but only first-order and excessively diffusive. Total-Variation-Diminishing (TVD) flux limiters blend the two: they recover second order in smooth regions and degrade gracefully to first order near extrema, provably not increasing the total variation `Σ|φ_{i+1}−φ_i|`. Pyrolysis.jl uses TVD limiting in the ALE advection (§13.8) and in second-order conservative remapping (§11, the ALE Framework).

### 13.4.1 The limiter enumeration

A small `@enum Limiter` selects the scheme (src/Discretization/flux_limiters.jl):

| Enum | `ψ(r)` | Character |
|---|---|---|
| `MINMOD` | `max(0, min(1, r))` | Most diffusive; strictly TVD |
| `VANLEER` | `(r + \|r\|)/(1 + \|r\|)` | Balanced; **default in ALE** |
| `SUPERBEE` | `max(0, min(2r, 1), min(r, 2))` | Least diffusive; compressive |
| `NONE` | `0` | First-order upwind (no limiting) |
| `CENTRAL` | `1` | Second-order central (no limiting) |

The dispatcher `apply_limiter(r, limiter)` returns the corresponding `ψ(r)`, defaulting to `minmod` for any unrecognized value.

### 13.4.2 Limiter functions and the TVD region

The slope ratio `r` measures the smoothness of the field across the upwind direction (defined in §13.5). Each limiter is a function `ψ(r)` of that ratio:

$$
\psi_{\text{minmod}}(r) = \max(0,\min(1,r)), \qquad
\psi_{\text{vanLeer}}(r) = \frac{r + |r|}{1 + |r|}, \qquad
\psi_{\text{superbee}}(r) = \max\bigl(0,\min(2r,1),\min(r,2)\bigr).
$$

All three satisfy `ψ(r) = 0` for `r ≤ 0` (extremum → first-order upwind), `ψ(1) = 1` (smooth → second order), and lie within Sweby's TVD region `0 ≤ ψ(r) ≤ min(2r, 2)`. Minmod hugs the lower edge of the region (maximal dissipation, maximal robustness); superbee hugs the upper edge (minimal dissipation, mild steepening/compression that can clip smooth extrema); van Leer is the smooth symmetric interior choice — it is the analytic harmonic mean of `r` and `1`, is differentiable for `r ≠ 0` (which matters for the AD-based Jacobian and sensitivities of §15/§17), and is the reason it is the package default.

> **Implementation note.** `minmod`, `vanleer`, `superbee`, `no_limiter`, `central_limiter`, and the `apply_limiter` dispatcher are all `@inline` scalar functions (src/Discretization/flux_limiters.jl). The hard `max`/`min` in minmod/superbee make those two non-differentiable at the kink points; van Leer is preferred when smooth derivatives are required end-to-end.

---

## 13.5 Limited Reconstruction: What Is Live

The TVD machinery in `src/Discretization/flux_limiters.jl` consists of the
limiter functions $\psi(r)$ (minmod, van Leer, superbee; §13.8.2) and the
`Limiter` enum / `apply_limiter` dispatcher. Their consumer is the mesh-relative
temperature advection's cell-centered limited gradient `_limited_gradient`
(§11.6.2, §13.8.2), which reconstructs $\partial f/\partial z$ donor-side with the
distance-correct blend $\beta = d_0/(d_0+d_1)$ and runs only under the
`w_override` verification hook.

There is no limited *face-value* reconstruction in production: the Eulerian
gas-advection face value is first-order upwind on the Darcy velocity
(`gas_flux_darcy_fick`, §9; a TVD upgrade is a documented accuracy-debt item),
and conduction/diffusion use the exact two-point face gradients of §13.7. The
exported `face_value_*` helpers of §13.6.3 (including the TVD-limited
`face_value_limited`, a fixed $\tfrac12$ blend valid on uniform meshes) are
uncalled building blocks.

---

## 13.6 Gradient Reconstruction

### 13.6.1 Cell-centered gradient on non-uniform meshes

Several volumetric source terms (notably the gas-advection source `S_conv` of §7.4) and the error/refinement indicators (§14) require a cell-centered gradient. The interior formula combines the two one-sided center-to-center slopes with **opposite-distance weights** (`opposite_distance_gradient`, src/Geometry/geometric_quantities.jl):

$$
\nabla\varphi_i = \frac{d_R\,s_L + d_L\,s_R}{d_L+d_R},
\qquad
s_L = \frac{\varphi_i - \varphi_{i-1}}{d_L},\quad
s_R = \frac{\varphi_{i+1} - \varphi_i}{d_R},
$$

with $d_L = z_i^c - z_{i-1}^c$ and $d_R = z_{i+1}^c - z_i^c$. Each one-sided slope estimates $\varphi'$ at its segment midpoint; the opposite-distance weighting interpolates them back to $z_i^c$, cancelling the leading $\varphi''$ truncation term — second order on arbitrary nonuniform meshes. On a uniform mesh this reduces to the familiar central difference $(\varphi_{i+1}-\varphi_{i-1})/(2\Delta z)$. The wide-stencil form $(\varphi_{i+1}-\varphi_{i-1})/(z_{i+1}^c-z_{i-1}^c)$ — which weights each slope by its *own* distance — carries a first-order $\tfrac12(d_R-d_L)\varphi''$ error on stretched meshes and is not used.

### 13.6.2 One-sided boundary gradients

At domain edges, one neighbor is missing, so the gradient degrades to a one-sided difference:

$$
\nabla\varphi_i =
\begin{cases}
\dfrac{\varphi_{i+1} - \varphi_i}{z_{i+1}^c - z_i^c}, & \text{bottom cell (left\_id = 0)},\\[10pt]
\dfrac{\varphi_i - \varphi_{i-1}}{z_i^c - z_{i-1}^c}, & \text{top cell (right\_id = 0)},\\[10pt]
0, & \text{isolated cell (both ids = 0; degenerate)}.
\end{cases}
$$

> **Implementation note.** `compute_gradients!` and the temperature-specialized `compute_temperature_gradients!` (src/Discretization/gradient_reconstruction.jl) implement these branches exactly, reading `state.cell_centers` so that the gradients carry AD sensitivities through the geometry. The interior `opposite_distance_gradient` reappears in the `S_conv` temperature-gradient computation (src/Physics/convection.jl) — which however upgrades the boundary cells to the surface-slope blend of §7.4.2 (the thermal BC supplies a solved $T_s$ there; a generic field has no boundary value, so the one-sided fallback above remains correct for `compute_gradients!`) — and in the AMR gradient indicators (src/Adaptivity/error_indicators.jl, r_refinement.jl).

### 13.6.3 Face-value interpolation schemes

For completeness the discretization module exposes a family of face-value building blocks (src/Discretization/gradient_reconstruction.jl); none are called by the production fluxes:

- **Central:** `face_value_central(φ_L, φ_R) = ½(φ_L + φ_R)` — second order, unlimited.
- **First-order upwind:** `face_value_upwind(φ_L, φ_R, v) = (v ≥ 0 ? φ_L : φ_R)`.
- **Distance-weighted:** `face_value_weighted(φ_L, φ_R, w_L) = w_L φ_L + (1−w_L) φ_R`, with `w_L = d_{f,R}/d_{LR}`.
- **Second-order via gradient:** `face_value_second_order` returns `φ_L + ∇φ_L d_{Lf}` for `v ≥ 0` and `φ_R − ∇φ_R d_{Rf}` for `v < 0` (MUSCL-style reconstruction to the face from the upwind cell center).
- **TVD-limited:** `face_value_limited(φ_LL, φ_L, φ_R, limiter) = φ_L + ½ψ(r)(φ_R − φ_L)` with `r = (φ_L − φ_LL)/(φ_R − φ_L)` (`compute_slope_ratio`) — a uniform-mesh ½ blend; exported but not called by any production flux.
- **Face gradient:** `face_gradient(φ_L, φ_R, d_{LR}) = (φ_R − φ_L)/d_{LR}`, with boundary faces zeroed in `compute_face_gradients!`.

These are building blocks; the conduction and gas-diffusion fluxes (§13.7) do not call them directly but use the equivalent closed-form face expressions.

---

## 13.7 Conductive, Gas-Diffusive, and Darcy–Fick Face Values

This section states the discrete face fluxes as they appear in the FVM. The physical derivations are in §7 (conduction, convection) and §9 (gas transport); here the focus is on the *discrete* form, the face averaging, and the signs.

### 13.7.1 Conductive heat flux

Fourier's law at a face, discretized with a distance-weighted harmonic-mean conductivity:

$$
q_f = -k_f\,\frac{T_R - T_L}{d_L + d_R}, \qquad
k_f = \frac{d_L + d_R}{\dfrac{d_L}{k_L} + \dfrac{d_R}{k_R}}.
$$

The harmonic mean is the **resistance-in-series** average: across an interface the thermal resistances `d_L/k_L` and `d_R/k_R` add, and continuity of heat flux then fixes `k_f` as above. On a uniform mesh `d_L = d_R` and this collapses to the standard harmonic mean `2k_Lk_R/(k_L+k_R)`; on stretched meshes the distance weighting improves near-boundary flux accuracy by roughly 1–5%. The minus sign is Fourier's law: `T_R > T_L` gives `q_f < 0` (heat flows toward `−z`, down the gradient).

> **Implementation note.** `conductive_flux(k_L, k_R, T_L, T_R, d_L, d_R)` computes `k_face = harmonic_mean_weighted(...)` then `-k_face*(T_R−T_L)/(d_L+d_R)`; the precomputed-face variant is `conductive_flux(k_face, T_L, T_R, d_LR) = -k_face*(T_R−T_L)/d_LR` (src/Physics/heat_conduction.jl). The assembly `compute_conductive_fluxes!` loops all faces using the precomputed `cache.k_face[f]` and `d_LR = |z_R^c − z_L^c|`, and **explicitly zeros boundary faces** (`left_id == 0 || right_id == 0 → q_cond[f] = 0`). The zero is then overwritten by the surface-energy-balance Newton solve (§12); zeroing first prevents double-counting (the interior operator must not see a boundary flux that the BC solver will also supply).

### 13.7.2 Gas-diffusive (volume-fraction-gradient) flux

The gas diffusive flux is written in **volume-fraction (equivalently partial-pressure) gradient** form. Starting from `J_g = −ρ_g λ ∂(ξ_g/ρ_g)/∂x` with the ideal-gas relation `ρ_g ∝ P/T`, this collapses to `J = −(λ/T)·∂(ξT)/∂x`, whose discrete face form is

$$
\boxed{\;J_f = -\lambda_f\,\frac{\xi_R T_R - \xi_L T_L}{T_f\,d_{LR}}\;}, \qquad T_f = \tfrac12(T_L + T_R).
$$

The product `ξT` is proportional to partial pressure (`P_j = ξ_j R_g T / M_j`), so this is diffusion down the partial-pressure gradient. The algebraic identity

$$
\xi_R T_R - \xi_L T_L = T_f\,(\xi_R - \xi_L) + \xi_f\,(T_R - T_L)
$$

shows the flux contains both an ordinary-diffusion term (`∝ ∂ξ/∂x`) and a thermal-diffusion / Soret-like term (`∝ ξ/T·∂T/∂x`); the single combined expression captures both without ever forming `ρ_g` at the face.

Crucially, the **harmonic mean is taken of `λ` alone, not of `λρ`**. The density dependence is carried implicitly by the mass concentrations `ξ` in the numerator, so averaging `λ` (a transport resistance) is the correct interface treatment:

$$
\lambda_f = \frac{d_L + d_R}{\dfrac{d_L}{\lambda_L} + \dfrac{d_R}{\lambda_R}}.
$$

> **Implementation note.** `gas_flux(λ_face, ξ_L, ξ_R, T_L, T_R, d_LR)` (src/Physics/gas_transport.jl) returns exactly `-λ_face*(ξ_R*T_R − ξ_L*T_L)/(T_face*d_LR)` with `T_face = (T_L+T_R)/2`. The face coefficient `λ_face` is the distance-weighted harmonic mean from `harmonic_mean_weighted` (src/Properties/effective_properties.jl). All gas species share the same `λ_face`; there is no per-species scaling of the diffusive term (§9.2).

### 13.7.3 Darcy–Fick combined advection–diffusion flux

When pressure-driven Darcy flow is active, the face flux adds a first-order-upwinded advective term to the diffusive term:

$$
J_f^{\text{DF}} = \underbrace{\frac{\xi^{\text{upwind}}}{\varphi_f}\,v_g}_{\text{advection}} \;+\; \underbrace{\Bigl(-\lambda_f\,\frac{\xi_R T_R - \xi_L T_L}{T_f\,d_{LR}}\Bigr)}_{\text{diffusion}},
\qquad
\xi^{\text{upwind}} =
\begin{cases}
\xi_L, & v_g \ge 0,\\
\xi_R, & v_g < 0,
\end{cases}
$$

with $\varphi_f=\max(\tfrac12(\varphi_L+\varphi_R),10^{-12})$ the face porosity: $\xi$ is per bulk volume and $v_g$ is superficial, so the advected scalar is the intrinsic density $\xi/\varphi$ (§9.3). The upwinding chooses the cell *from which* the bulk gas flows, which is the monotone first-order choice. The bulk Darcy velocity at the face is

$$
v_g = -\text{mob}_f\,\frac{P_R - P_L}{d_{LR}},\qquad
\text{mob}_f = \frac{d_L+d_R}{d_L\mu_L/\kappa_L + d_R\mu_R/\kappa_R},
$$

the **distance-weighted harmonic mean of the cell mobilities** $\kappa/\mu$, with `mob_f = 0` exactly when either side's mobility is below $10^{-30}$ (an impermeable side is an infinite series resistance) and `μ > 10^{-20}` guards per cell. The ideal-gas pressure closure `P = (Σ_g ξ_j R_g T/M_j)/φ` (porosity floor `10^{-12}`) and the pressure gradient are derived in §9; here they enter only as the velocity driving the advective face flux. All gas species share the same `v_g` (advection carries every species equally) and the same `λ_face` (§9.2).

> **Implementation note.** `gas_flux_darcy_fick(v_g, ξ_L, ξ_R, φ_L, φ_R, λ_face, T_L, T_R, d_LR)` (src/Physics/gas_transport.jl) computes `J_diff` exactly as above, sets `ξ_upwind = v_g ≥ 0 ? ξ_L : ξ_R`, and returns `J_adv + J_diff` with `J_adv = v_g·ξ_upwind/φ_face`. The face velocity and mobility come from `update_darcy_properties!` (src/Discretization/finite_volume.jl); the Darcy RHS path is `compute_rhs_darcy_fick!`, which differs from the Fickian path in computing pressure/velocity first, calling the Darcy–Fick flux instead of the diffusion-only flux, and writing the boundary-face advective venting flux from the pressure BCs (§9.5.3).

### 13.7.4 Solid/liquid components carry no transport flux

Solid and liquid components are mesh-bound (Lagrangian): they are not transported by diffusion or Darcy flow. Their face fluxes are identically zero, so `divergence_species_1d!` produces zero for them and the species residual reduces to reactions only (§13.3.4). This is what locks solid mass-thickness `ρ_s·Δz` to the mesh and makes the ALE motion (§13.8) the mechanism by which solids are advected. The gas fluxes of §13.7.2–§13.7.3 are defined relative to that skeleton, so under ALE the gas rides the mesh too and the same divergence is the whole of its transport (§11.1.2).

### 13.7.5 The energy-storage convention (Formulation A)

The conductive divergence feeds the energy equation with the mixture heat capacity `ρc_{p,\text{eff}} = Σ_j ξ_j c_{p,j}(T)` in the denominator — every component, gas included, stores sensible heat. The enthalpy carried by the gas flux relative to the skeleton is the volumetric advection source `S_conv` (always on in `:standard` energy form) and, optionally, `S_gen` (in `:with_generation_sink`). This is the ThermaKin Eq. 17 convention. The discretization layer is agnostic to this choice — it simply assembles `div_q`, `S_conv`, `S_gen`, `Q_rad`, `Q_rxn` and divides by `prop_cache.ρcp[i]`. See §7.3–§7.6 for the full derivation, and §13.3.4 for the assembled formula.

---

## 13.8 ALE Flux Discretization (Recap)

This section summarizes how the moving mesh modifies the spatial operators; the full ALE theory (mesh-velocity construction, GCL, remapping) is in §11, The ALE Framework.

### 13.8.1 The ALE advection correction

The arbitrary Lagrangian–Eulerian time-derivative identity for a moving sample point is

$$
\frac{\partial f}{\partial t}\bigg|_{\chi} = \frac{\partial f}{\partial t}\bigg|_{x} + w\,\frac{\partial f}{\partial z},
$$

so the moving-mesh residual equals the fixed-mesh (Eulerian) residual plus a mesh-motion advection term `w·∇f`. In the discretization this correction is computed cell-by-cell and **added** to the assembled Eulerian RHS:

$$
\frac{dT_i}{dt}\bigg|_{\text{ALE}} = \frac{dT_i}{dt}\bigg|_{\text{Euler}} + w_{\text{cell},i}\,\Bigl(\frac{\partial T}{\partial z}\Bigr)_i,
\qquad
w_{\text{cell},i} = \tfrac12\,(w_{n_L} + w_{n_R}).
$$

The cell velocity is the average of its two node velocities; node velocities are built by bottom-up integration of the per-cell dilation (`w_1 = 0` pinned at the substrate; §11.3).

### 13.8.2 Upwind, TVD-limited gradient for the advection term

Because the correction `w·∇φ` is *added* to the Eulerian RHS, the mesh-frame equation is `∂φ/∂τ + c·∂φ/∂z = RHS` with apparent transport velocity `c = −w` (the lab-frame field drifts at `−w` relative to the mesh). The advected gradient `∂φ/∂z` is therefore reconstructed **donor-side with respect to `c = −w`** — the right neighbor for `w > 0`, the left for `w < 0` — with an optional TVD limiter (default `VANLEER`). The first-order donor-side gradient, always written in the `+z` sense, is

$$
\text{grad}_1 =
\begin{cases}
(T_{\text{donor}} - T_i)/d, & w_{\text{cell}} > 0 \;(c<0,\ \text{donor on the right}),\\[4pt]
(T_i - T_{\text{donor}})/d, & w_{\text{cell}} < 0 \;(c>0,\ \text{donor on the left}),
\end{cases}
$$

and the second-order TVD correction, using the next-donor-side gradient `grad_0`, the **near/far** ratio `r = grad_1/grad_0` (see §11.6.2 for why the inverse ratio is unbounded), and the distance-correct blend `β = d₀/(d₀+d₁)` (½ on uniform meshes; §11.6.2), is

$$
\text{grad}_{\text{TVD}} = \text{grad}_1 + \beta\,\psi(r)\,(\text{grad}_1 - \text{grad}_0).
$$

The correction is dropped (first order only) when there is no donor-donor cell, when `limiter == NONE`, or when `|grad_0| < 10^{-20}` (flat far slope — ratio undefined). Boundary neighbors use zero-gradient (ghost = cell value) extrapolation; in the standard shrinking case (`w < 0`) the surface cell's donor is the interior neighbor, so the receding surface receives a real gradient rather than the boundary ghost's zero. A negligible mesh velocity `|w_cell| < 10^{-20}` short-circuits the advection to zero.

> **Implementation note.** `compute_ale_advection_term!` (energy) and `_limited_gradient` (src/Discretization/ale_fluxes.jl) implement these formulas; the assignment is `advection[i] = w_cell * dT_dz` and this is added to `dT` in `_apply_ale_T_advection!` (src/Residual/residual.jl) — **only under `Workspace.w_override`**, with `w = w_rel = w_applied − w_physical` (§11.2.2): production ALE has no energy advection because the matrix carrying the thermal mass moves with the mesh.

### 13.8.3 Species: dilation only

Every component rides the mesh (§11.1.2), so the species equations carry no mesh-motion advection at all. The Eulerian species RHS (`-div_J + dξ_rxn` for gases, `dξ_rxn` for solids/liquids) is used unchanged, and the residual adds the dilation `−ξ_{j,i}·θ_i` to every component (§11.2.3, §11.3.4): per cell, `d(ξ_j V)/dt = S_j V − [J_j A]_faces`, so each component's mass changes only by reaction and by its skeleton-relative flux through the faces, and the column total changes only through the boundary faces. In a sealed collapsing column the pore gas compresses cell by cell (V-CPL-03 P4 pins exact conservation); at an open face the gas leaves only through the mass BC's diffusive or Darcy flux.

> **Implementation note.** `_apply_ale_dilation_correction!` (src/Residual/residual.jl): `dξ[j, i] -= state.ξ[j, i] * θ[i]` for every component, with `θ` the cell's reaction-driven strain rate. Its Jacobian is the cell-local block of §15.8.4.

---

## 13.9 Mesh Generation

### 13.9.1 Uniform mesh

The simplest mesh divides `[0, L]` into `n` equal cells:

$$
\Delta z = \frac{L}{n}, \qquad z_i = (i-1)\,\Delta z \;\; (i = 1,\dots,n+1), \qquad V_i = \Delta z\,A.
$$

Node 1 (`z = 0`) is the bottom boundary (`boundary_id = 2`); node `n+1` (`z = L`) is the top (`boundary_id = 1`). Faces, cells, the active-cell mask (all `true`), the per-cell dry-solid-mass denominator `m_dry_solid_initial` (zeroed, populated later at problem setup), and the boundary-face map (`1 => [n+1]`, `2 => [1]`) are all initialized, and `update_geometry!` computes the cached centers/volumes/areas.

> **Implementation note.** `create_uniform_mesh(L, n_cells, Val(NC); T_initial, ξ_initial, cross_section_area)` (src/Geometry/mesh_generation.jl).

### 13.9.2 Geometric-stretch mesh

To resolve the steep thermal and reaction fronts near the heated surface without paying for fine resolution everywhere, a geometric stretch concentrates cells near one boundary. With a stretch ratio `s` (`Δz_{i+1}/Δz_i = s`), the cell widths form a geometric series whose sum is `L`:

$$
L = \sum_{i=1}^{n}\Delta z_i = \Delta z_0\,\frac{s^{\,n} - 1}{s - 1}
\;\;\Longrightarrow\;\;
\Delta z_0 = L\,\frac{s - 1}{s^{\,n} - 1}, \qquad \Delta z_i = \Delta z_0\,s^{\,i-1}.
$$

Node positions are accumulated `z_{i} = z_{i-1} + Δz_{i-1}` and then **renormalized** by `L/z_{n+1}` so that `z_{n+1} = L` exactly (removing floating-point drift in the partial sums). For `|s − 1| < 10^{-10}` the code falls back to the uniform `range(0, L, n+1)`.

The `stretch_from` keyword names the boundary that gets the **finest** cells. The series is always built fine-end-first from `z = 0` (smallest cell `Δz_0` at the bottom, cells growing toward `z = L`); with the default `stretch_from = :top` the positions are then reflected so the finest cells sit at the exposed surface `z = L`:

$$
\tilde z_i = L - z_{n+2-i}.
$$

`stretch_from = :bottom` keeps the unreflected series (finest cells at the substrate). The constructor requires `stretch ≥ 1` and `stretch_from ∈ (:top, :bottom)`; orientation is selected only by `stretch_from`, never by a sub-unity ratio.

> **Implementation note.** `create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from=:top, ...)` (src/Geometry/mesh_generation.jl). The geometric series uses `Δz_i = Δz_0 * stretch^(i-2)` in the node-accumulation loop (so the first interior gap is `Δz_0 * s^0 = Δz_0`), followed by the `positions .*= L/positions[end]` renormalization, and the `:top` reflection `positions = L .- reverse(positions)`. For convergence studies on stretched meshes see §18 (Verification & Validation), e.g. `convergence_douglas_fir`.

---

## 13.10 State-Vector Packing and Geometry-from-State

### 13.10.1 Block-major layout `[T | ξ | z | χ]`

The flat ODE/DAE state vector `u` is packed **block-major**: all temperatures first, then each component's concentration block, then (in ALE) the node positions, then (in WithChi) the pyrolysis-progress block. With `n = n_cells`:

$$
\underbrace{T_1,\dots,T_n}_{\text{T block: }1..n}\;\bigl|\;
\underbrace{\xi_{1,1},\dots,\xi_{1,n}}_{\xi_1}\,,\dots,\,\underbrace{\xi_{N_C,1},\dots,\xi_{N_C,n}}_{\xi_{N_C}}\;\bigl|\;
\underbrace{z_1,\dots,z_{n+1}}_{\text{z block (ALE)}}\;\bigl|\;
\underbrace{\chi_1,\dots,\chi_n}_{\text{$\chi$ block (WithChi)}}.
$$

The exact ranges (1-based) are:

$$
T:\;1..n, \qquad
\xi_j:\;n + (j-1)n + (1..n), \qquad
z:\;n(1+N_C) + (1..n_{\text{nodes}}), \qquad
\chi:\;\text{trailing }n.
$$

Block presence is determined **strictly** by the `SimulationMode` trait: `T` and `ξ` always present; `z` iff `mesh_mode(mode) isa ALE`; `χ` iff `chi_mode(mode) isa WithChi`. The total length is `n(1+N_C)` plus `n_nodes` (ALE) plus `n` (WithChi).

> **Implementation note.** `StateLayout{NC,Mode}` (src/Problem/state_layout.jl) is the single source of truth for every pack/unpack/view. Accessors: `temperature_range`, `concentration_range(l, j)`, `z_range`, `chi_range`, `z_offset`, `chi_offset`, `state_length`, plus the compile-time `block_view(u, l, Val{:T}|Val{:ξ}|Val{:z}|Val{:χ})`. `pack_state!` writes `mesh.cell_states` and (ALE) node positions into `u`; the χ block is initialized to zero because χ is reconstructed each residual rather than integrated as raw state (§10.3).

### 13.10.2 Geometry-from-state during assembly

In ALE mode the geometry is a function of the `z` block, so it must be recomputed from `u` *without mutating the mesh* (the residual and Jacobian must be pure functions of `u`). The geometry kernel reads node positions directly from `u` at the `z_offset`:

$$
V_i = |u_{z_{\text{off}}+i+1} - u_{z_{\text{off}}+i}|\,A, \quad
z_i^c = \tfrac12\,(u_{z_{\text{off}}+i} + u_{z_{\text{off}}+i+1}), \quad
z_f = u_{z_{\text{off}}+f},
$$

with face distances assembled afterward (interior: center-to-center; boundary: factor-2 ghost form of §13.2.3). Because `A` and the `u`-derived positions carry the live scalar type, the entire geometry — and every flux and divergence built on it — is differentiable end-to-end.

> **Implementation note.** `compute_geometry_from_state!`, `compute_face_distances_from_state!`, and the unified `update_geometry_cache!` populate a preallocated `GeometryCache` (src/Geometry/geometric_quantities.jl), called once per RHS evaluation in ALE mode. Between *accepted* integrator steps, `update_mesh_from_state!` (src/Problem/state_layout.jl) mirrors the state back into `mesh.nodes[i].z` for AMR, depletion checks, and post-processing. The `StateCache` (src/Problem/state_cache.jl) holds `T`, `ξ`, `z`, `χ`, the geometry quantities `volume`/`A_face`/`A_section`, and `cell_centers`, all parameterized on the scalar type `Tu`.

### 13.10.3 RHS cache and resizing

The face fluxes, sources, and divergences are stored in a preallocated `RHSCache{Tu}` (src/Discretization/finite_volume.jl) carrying `q_cond`, `q_total`, `S_conv`, `S_gen`, `J_gas` (`NC × n_faces`), `Q_rxn`, `Q_rad`, `dξ_rxn` (`NC × n_cells`), the reaction workspaces `rates`/`heats`/`species_src`, and the divergences `div_q`/`div_J`. The type parameter `Tu` is `Float64` for ordinary evaluation and `ForwardDiff.Dual` for AD-based Jacobian/sensitivity passes, so no reallocation occurs when switching modes. After an AMR resize the cache is resized in lockstep via `resize_rhs_cache!`, and the geometry cache via `resize_geometry_cache!` (§14).

---

## 13.11 Comparison to Gpyro and ThermaKin2Ds Discretizations

We map the reference codes' notation onto the unified nomenclature of §2 on first use, then contrast the discretization choices.

### 13.11.1 Gpyro (Lautenberger & Fernandez-Pello)

Gpyro is a 1D (optionally pseudo-2D) generalized condensed-phase pyrolysis solver using a **Patankar-style finite-volume** discretization (control volumes about cell-centered nodes). Mapping its symbols to ours: Gpyro's pre-exponential `Z → A_i`, permeability `K → κ`, kinematic viscosity `ν` (in `ṁ'' = −(K/ν)(∂P/∂z − gρ_g)`), conversion `α_conv`, and the overbarred mixture density `ρ̄`.

- **Conservation form vs. non-conservative.** Gpyro writes integral control-volume balances of `ρ̄Δz`, `ρ̄Y_iΔz`, gas mass, and energy on volumes that deform with the solid. Pyrolysis.jl uses the same cell-centered FVM divergence operator (§13.3) on the ALE mesh, with every component riding the mesh, but solves the **non-conservative (Formulation A)** temperature equation, storing every component in `ρc_p` and carrying the gas-flux enthalpy via `S_conv` (§13.7.5). Both telescope to the same boundary-flux balance (§13.3.3).
- **Interface interpolation.** Patankar's hallmark is the **harmonic interface conductivity** for diffusion coefficients, justified by flux continuity across a control-volume face — identical in spirit to the harmonic-mean `k_f`, `λ_f`, `κ_f` of §13.7.1–§13.7.3. Pyrolysis.jl additionally distance-weights the harmonic mean for non-uniform meshes.
- **Convection scheme selection.** Patankar/Gpyro selects the face scheme by the cell **Péclet number** (central below ≈2, upwind/power-law above), a built-in switching. Pyrolysis.jl does not auto-switch on Péclet: diffusion is always the volume-fraction-gradient flux and gas advection is first-order upwind (Darcy–Fick) or TVD-limited (ALE). Choosing resolution to keep the cell Péclet number moderate is left to the user (see §18.7).
- **Linear solver.** Gpyro is fully implicit (backward Euler) with Patankar linearization solved by **TDMA (Thomas algorithm)** for the tridiagonal system. Pyrolysis.jl delegates time integration to stiff OrdinaryDiffEq integrators with a structured/sparse Jacobian (§15); a TDMA appears internally only in the P1 radiation solve (§8.7).
- **Pressure/Darcy.** Gpyro's Darcian gas flux `ṁ'' = −(K/ν)(∂P/∂z − gρ_g)` includes buoyancy; Pyrolysis.jl's `v_g = −(κ/μ)∇P` (§13.7.3) omits the gravity term (1D through-thickness, gravity negligible) and uses dynamic viscosity `μ` rather than kinematic `ν`.

### 13.11.2 ThermaKin2Ds (Stoliarov et al.)

ThermaKin is a 1D (with 2D extensions) finite-difference pyrolysis solver. Mapping its symbols: stoichiometric coefficient `θ_i^j → ν_{i,j}` (overload A2; ThermaKin's `θ` is stoichiometry, ours is strain rate), and its published heat-of-reaction sign is `h > 0 = exothermic` (overload H1; internally we use `h > 0` endothermic with `Q_rxn = −Σh_r r_r`).

- **Gas transport.** Pyrolysis.jl's volume-fraction-gradient flux `J = −(λ/T)∂(ξT)/∂x` (§13.7.2) is taken directly from ThermaKin's Section 1.2 gas-transport formulation `J = −ρ_g λ ∂(ξ_g/ρ_g)/∂x`. For a rigid non-swelling matrix (`γ = 0`) ThermaKin recovers a Boyle's-law Darcy form; Pyrolysis.jl provides the explicit Darcy–Fick path (§13.7.3) for that regime.
- **Convective energy term.** ThermaKin's Eq. 17 carries gas sensible enthalpy as `c_g (J_g·∂T/∂x)`, which is exactly the `S_conv` source of §13.7.5 (the default S_conv-only convention), and sums every component in the storage term — both as here.
- **Time discretization.** ThermaKin uses Crank–Nicolson in the primary (through-thickness, `x`) direction and explicit updates in the transverse `(y, z)` directions; the spatial operators are second-order finite differences on a structured grid. Pyrolysis.jl is finite-**volume** (flux-based, conservative-divergence) rather than finite-difference, and offloads temporal integration to adaptive stiff solvers (§15).
- **Moving mesh.** ThermaKin performs dynamic mesh adjustment around a single reference point as the material deforms; Pyrolysis.jl recasts the same deformation through the ALE framework (§11) with cumulative mesh velocity built from per-cell dilation, and (optionally) a separate lateral-shrinkage `χ` channel (§10).

The shared lineage (harmonic interface averaging, S_conv gas-enthalpy convention, volume-fraction-gradient gas flux) is intentional — Pyrolysis.jl is designed to reproduce ThermaKin/Gpyro physics while adding type-stable, allocation-free, AD-compatible assembly and a structured Jacobian.

---

## 13.12 Limitations and Open Issues

The following discretization-specific issues are open (drawn from the source annotations and the verification against the implementation):

1. **ALE advection limiter is not exposed through `solve`.** `VANLEER` is the default limiter in `compute_ale_advection_term!` (the mesh-relative temperature advection, override-only). The limiter is selectable only at the workspace level (`build_workspace(...; ale_limiter = MINMOD)`); there is no `solve`-level keyword to switch it (minmod/superbee/none), which would aid verification studies.

2. **First-order Darcy–Fick advection only.** The pressure-driven advective face flux uses first-order upwinding (§13.7.3); no TVD/MUSCL option exists for Darcy advection (the TVD machinery is wired only into the ALE advection and the conservative remapping). Sharp pressure-driven composition fronts are therefore smeared at first order unless the mesh is refined.

3. **`T_ref` for the volumetric generation sink under asymmetric ambients.** When `energy_form = :with_generation_sink` and the two boundaries have different ambient temperatures with a permeable back face, the single-scalar `T_ref` (resolved from one boundary) breaks exact energy closure at the opposite boundary by `Σ_g ṁ_g c_{p,g}(T_{∞,\text{top}} − T_{∞,\text{bottom}})` (§7.5). This is harmless for the typical impermeable back face but would require a per-cell datum to fully resolve.

4. **Multi-dimensional generalization.** All geometry, the divergence operator, the GCL, and the face structures are specialized to 1D (one face per node, scalar face areas, all faces sharing `A`). The `MeshTopology{D,...}` and `MeshGeometry{D}` types are parameterized by dimension `D` but only `D = 1` is implemented; 2D/3D would require vector node positions, per-face normals and variable areas, and a full divergence-theorem GCL.

For the temporal integration and Jacobian that consume the operators of this chapter, see §15 (Temporal Integration & Structured Jacobian); for adaptive refinement that resizes the mesh and these caches, see §14 (Adaptive Mesh Refinement); for the conservation audits that exploit the telescoping divergence of §13.3.3, see §16 (Conservation Diagnostics).
