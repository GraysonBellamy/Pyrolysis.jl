# 15. Temporal Integration, DAE/Residual Formulation, and the Structured Jacobian

This chapter describes how the spatially-discretized governing equations of §13
(Finite-Volume Spatial Discretization) are assembled into a single flat
ODE/DAE right-hand side, integrated in time by stiff implicit Runge–Kutta
methods, and linearized through a structured-operator Jacobian. We derive the
block-major state-vector layout and the four dispatched residual paths; the
unified `residual!` entry point and the immutable-definition / mutable-workspace
split that makes it pure and thread-safe; the full ALE residual assembly
sequence; the link to OrdinaryDiffEq's implicit-RK $W$-matrix; the
three-layer structured-Jacobian architecture (operators → couplings →
linear-solve strategies); the exact per-operator Jacobian formulas with their
sign conventions; the O(N) prefix/suffix/rank-1 couplings that represent the
non-local terms; the linear-solve strategies and iterative solvers; the
two-pattern sparsity construction with greedy coloring; and the verification
machinery that anchors every backend to a dense AD reference. Every discrete
rule is grounded in the implementing function so the mathematics is traceable to
the code; where the code makes an approximation, a guard, or a sign choice, we
state it and trust the code over any idealized form.

Conventions used in this chapter (binding, from §2, Nomenclature):

- The spatial coordinate is $z$ (axial, through-thickness). **$z=0$ is the
  bottom/substrate** (boundary id 2, tag `:bottom`); **$z=L$ is the top/exposed
  surface** (boundary id 1, tag `:top`). Cell indices increase from bottom to
  top.
- **Flux sign convention:** positive flux is transport in the $+z$ direction.
- **Divergence convention:** the discrete divergence $(F_R-F_L)A/V$ is positive
  when the quantity flows *out* of the cell; the energy RHS carries
  $-(\nabla\cdot\mathbf{q})$.
- **Heat-of-reaction sign convention:** $h>0$ endothermic (cools), $h<0$
  exothermic (heats); the volumetric source is $Q_{\rm rxn}=-\sum_r h_r r_r$
  (overload resolution H1 of §2).
- In this chapter $i$ is the **cell index** ($1,\dots,n$), $f$ a face index, $j$
  a component index, $r$ a reaction index, $k$ either a node index or a generic
  loop index where stated. In the linear-algebra sections, $\mathsf{A}$ denotes
  the sparse Jacobian block (not the cross-sectional area $A$); the disambiguation
  follows overload resolution A3 of §2 and is stated where both appear.

The state vector is denoted $\boldsymbol{u}$; its time derivative is
$\boldsymbol{du}=d\boldsymbol{u}/dt$. The residual function (in the
method-of-lines sense, an *explicit* RHS, not an implicit DAE residual
$F(\boldsymbol{u},\boldsymbol{du},t)=0$) is `residual!`, which writes
$\boldsymbol{du}=\mathrm{RHS}(\boldsymbol{u},t)$.

---

## 15.1 State layout and the simulation-mode trait

### 15.1.1 Block-major packing

All time-dependent unknowns are packed into a single flat vector
$\boldsymbol{u}\in\mathbb{R}^{N_{\rm state}}$ in **block-major** order, with the
block set determined by the active `SimulationMode`. The descriptor is
`StateLayout`:

$$
\boldsymbol{u} =
\big[\,\underbrace{T_1,\dots,T_n}_{T\ \text{block}}\ \big|\
\underbrace{\xi_{1,1},\dots,\xi_{1,n}}_{\xi_1\ \text{block}},\dots,
\underbrace{\xi_{N_C,1},\dots,\xi_{N_C,n}}_{\xi_{N_C}\ \text{block}}\ \big|\
\underbrace{z_1,\dots,z_{n_{\rm nodes}}}_{z\ \text{block (ALE)}}\ \big|\
\underbrace{\chi_1,\dots,\chi_n}_{\chi\ \text{block (WithChi)}}\,\big].
$$

The explicit index map (used verbatim by every pack/unpack routine and by the
Jacobian's `FastSparsityMap`) is

$$
\begin{aligned}
\text{T index}        &: \quad i, &&i=1,\dots,n,\\[2pt]
\text{$\xi_{j,i}$ index} &: \quad n + (j-1)\,n + i, &&j=1,\dots,N_C,\ i=1,\dots,n,\\[2pt]
\text{$z_k$ index}     &: \quad n(1+N_C) + k, &&k=1,\dots,n_{\rm nodes}\ \text{(ALE only)},\\[2pt]
\text{$\chi_i$ index}  &: \quad n(1+N_C) + \delta_{\rm ALE}\,n_{\rm nodes} + i, &&i=1,\dots,n\ \text{(WithChi only)},
\end{aligned}
$$

and the total length is

$$
N_{\rm state} = n\,(1+N_C) + \delta_{\rm ALE}\,n_{\rm nodes} + \delta_{\rm WithChi}\,n,
$$

with $\delta_X=1$ if the mode includes block $X$ and $0$ otherwise. In 1D
$n_{\rm nodes}=n+1$. The $\chi$ block is **appended after** the $z$ block, not
interleaved; its offset arithmetic is therefore dense and is the single source
of subtle indexing bugs the code guards against by always going through the
range accessors `temperature_range`, `concentration_range(l,j)`, `z_range(l)`,
and `chi_range(l)`.

> **Implementation note.** The layout, accessors, and the
> `pack_state!`/`update_mesh_from_state!`/`create_state_vector` operations live
> in `StateLayout`, src/Problem/state_layout.jl. The same packing is mirrored in
> the Jacobian sparsity module's index macros `_temp_idx`, `_xi_idx`, `_z_idx`,
> `_chi_idx`, src/Jacobian/sparsity.jl, and in solution extraction,
> src/Solver/solve.jl. The block order $[T\mid\xi\mid z\mid\chi]$ is binding: any
> code that reads or writes raw state must respect it.

The block-major (rather than cell-major / interleaved) choice is deliberate. It
makes each physical field a contiguous slice, which (i) lets `block_view` expose
a strided sub-vector for diagnostics and callbacks without
copying, (ii) keeps the conduction/transport stencils as small offsets within a
block, and (iii) makes the ALE geometry coupling a clean rectangular sub-block
(the $z$-columns), which the two-pattern sparsity construction of §15.11
exploits.

The experimental P1 radiation path (quarantined to
`Pyrolysis.Experimental.RadiationP1`) uses this *same* layout: the
incident-radiation field $G$ is an algebraic (quasi-steady) field stored in
the `RadiationCache` and re-solved at every residual evaluation, never a
state block (see §8.12). The layout above is therefore the canonical one for
every mode, production and experimental alike.

### 15.1.2 The bundled `SimulationMode` trait

Rather than carry four orthogonal trait axes and suffer combinatorial method
explosion, Pyrolysis.jl bundles all four into a single type parameter:

$$
\texttt{SimulationMode}\{\,\underbrace{\text{Mesh}}_{\rm Eulerian/ALE},\ \underbrace{\text{Rad}}_{\rm NoRadiation/SurfaceAbsorption/BeerLambert},\ \underbrace{\text{Trans}}_{\rm Fickian/DarcyFick},\ \underbrace{\text{Chi}}_{\rm NoChi/WithChi}\,\}.
$$

The four axes are:

| Axis | Type | Values | Effect |
|---|---|---|---|
| Mesh motion | `MeshMode` | `Eulerian`, `ALE` | $z$ block present iff ALE; geometry constant iff Eulerian |
| Radiation transport | `RadiationMode` | `NoRadiation`, `SurfaceAbsorption`, `BeerLambert` | volumetric source $Q_{\rm rad}$ and Jacobian non-locality |
| Gas transport | `TransportMode` | `Fickian`, `DarcyFick` | diffusion-only vs. pressure-driven Darcy + diffusion |
| Cross-section | `ChiMode` | `NoChi`, `WithChi` | $\chi$ block present iff `WithChi`; $A(t)=A_0\,\mathrm{law}(\bar\chi)$ |

Accessors extract the type parameters at compile time:

```julia
@inline mesh_mode(::SimulationMode{M,R,T,C})      where {M,R,T,C} = M()
@inline radiation_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = R()
@inline transport_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = T()
@inline chi_mode(::SimulationMode{M,R,T,C})       where {M,R,T,C} = C()
```

Because the mode is a type parameter, `mesh_mode(def.mode) isa ALE` and friends
resolve to compile-time constants, and the residual's dispatch chain
specializes fully — no runtime branching on the mode. The four *combinations
that actually matter for the residual entry* are the cross product of mesh mode
(Eulerian/ALE) with transport mode (Fickian/DarcyFick); radiation and $\chi$
are further dispatched inside those paths. This yields the **four dispatched
residual paths** referenced in the scope:

1. Eulerian × Fickian, 2. Eulerian × DarcyFick, 3. ALE × Fickian,
4. ALE × DarcyFick — with `NoChi`/`WithChi` and the radiation model
specializing each.

> **Implementation note.** `SimulationMode`, the `MeshMode`/`RadiationMode`/
> `TransportMode`/`ChiMode` abstract hierarchies, and the accessors are in
> src/Problem/traits.jl. The mode is *derived*, not user-set directly:
> `to_problem_def(problem; use_ale, radiation_model, energy_form, …)`
> (src/Problem/adapter.jl) assembles it from the `use_ale` flag, the radiation
> model, and the material's Darcy/$\chi$ traits (`has_darcy_flow`, the presence
> of a `lateral_shrinkage_law`). The `RadiationModel` enum value
> `P1_QUASI_STEADY` is rejected at solve setup (it lives in `Experimental`), so
> the production `RadiationMode` axis spans only the first three values.

---

## 15.2 The problem definition / workspace separation

The single most important structural decision in the temporal layer is that the
residual is a **pure function** of state and time, parameterized by an
*immutable* problem definition and a *mutable* scratch workspace:

```julia
residual!(du, u, def::ProblemDef, t, ws::Workspace)  ->  nothing
```

- `def::ProblemDef` is **immutable** and shareable across threads / sweep
  workers. It strips the mutable mesh down to its connectivity (`MeshTopology`)
  and retains only the *unchanging* specification: the material, BC set, initial
  conditions, time span, the `SimulationMode`, and a NamedTuple of runtime
  options (`radiation_model`, `energy_form`, `use_transpiration_bc`, …).
- `ws::Workspace` is **mutable** and owned by exactly one worker. It holds the
  typed state mirror (`StateCache`), property cache (`PropertyCache`), RHS work
  arrays (`RHSCache`), and (in ALE mode) the ALE work arrays. It also caches the
  dry-gas component indices and the boundary face/cell IDs.
- `t` is a plain `Real`; the kernel accepts whatever scalar type $\boldsymbol{u}$
  and $\boldsymbol{du}$ carry — `Float64` for ordinary evaluation, or
  `ForwardDiff.Dual` for AD Jacobian assembly and forward sensitivity (§17).

This split is what makes parameter sweeps (§ User Guide UG-12) and pluggable
AD Jacobians possible without re-allocating problem metadata: every worker shares
one `def` and owns one `ws`; the AD Jacobian builds a *parallel* `Dual`-typed
workspace (`build_ad_workspace`) that shares the mesh, layout, and boundary
indices but carries independent `Dual`-typed caches.

A defining invariant is that **the mutable `Mesh1D` is read-only during residual
and Jacobian assembly** — connectivity only. All state and all state-derived
geometry (cell centers, volumes, face areas, cross-section) live in the typed
`StateCache`, which is populated once at residual entry by `populate!` and
carries `eltype(u)` end-to-end. The mutable mesh is re-synchronized from the
accepted integrator state *only between steps*, by the `StepUpdateCallback`
(§15.5). This guarantees that the many residual evaluations inside one implicit
step (Newton iterations, finite-difference Jacobian coloring) all see a
consistent geometry, and that AD sensitivities propagate cleanly through
geometry because the geometry is a function of `cache.z` and `cache.A_section`,
both `Dual`-typed in the AD pass.

> **Implementation note.** `ProblemDef` and the runtime-option NamedTuple:
> src/Problem/problem_def.jl. `Workspace`, `build_workspace`, and
> `build_ad_workspace`: src/Problem/workspace.jl. `StateCache`, `populate!`,
> `update_state_geometry!`, `update_cross_section_area!`:
> src/Problem/state_cache.jl. The user-facing `PyrolysisProblem` and the
> `to_problem_def` adapter: src/Problem/pyrolysis_problem.jl and
> src/Problem/adapter.jl.

### 15.2.1 The `energy_form` / transpiration guard

The `to_problem_def` adapter enforces the Formulation-A energy convention of §7
at construction time. Two options interact:

- `energy_form = :standard` (default): the gas-advection source
  $S_{\rm conv}=-\sum_g c_{p,g}\,\bar{J}_{g}\,\partial T/\partial z$ is the
  *only* gas-enthalpy carrier; the gas-generation enthalpy sink $S_{\rm gen}$ is
  zeroed. This matches FDS/Gpyro/ThermaKin.
- `energy_form = :with_generation_sink`: adds the cell source
  $S_{\rm gen,i}=-\sum_g \Delta h_g(T_i,T_{\rm ref})(\nabla\cdot\mathbf{J}_g)_i$,
  a *partial* move toward conservative form (full conservation would also
  restore gas storage to $\rho c_p^{\rm eff}$).

The surface transpiration enthalpy term $\sum_g \dot{m}_g\,\Delta h_g(T_s,T_\infty)$
overlaps $S_{\rm conv}$ exactly by $\int S_{\rm conv}$, so combining both
double-counts the gas enthalpy and produces a spurious $\sim5\text{–}10\,\mathrm{kW/m^2}$
surface sink. The adapter therefore **rejects** `use_transpiration_bc=true` at
construction (and `solve` rejects it in `_validate_solve_kwargs`). The supported
production path is $S_{\rm conv}$-only; see §7.6 and §16.4 for the energy-ledger
consequences.

---

## 15.3 The unified residual and its dispatch

The top-level entry is a one-line dispatcher on the mesh mode:

```julia
@inline function residual!(du, u, def::ProblemDef, t, ws::Workspace)
    _residual!(du, u, def, t, ws, mesh_mode(def.mode))
    return nothing
end
```

### 15.3.1 Eulerian path

For a fixed mesh the geometry is constant (mirrored at workspace construction),
so the residual is a pure function of $(\boldsymbol{u},t)$ in four steps:

1. `populate!(ws.state_cache, u, layout)` — unpack $T,\xi$ from $\boldsymbol{u}$
   into the typed cache (single source of truth).
2. `Properties.update_properties!` — evaluate all state-dependent properties
   ($k_{\rm eff}$, $\rho c_p^{\rm eff}$, $\lambda$, $\varepsilon$, $\alpha$, and
   face-averaged $k_{\rm face}$, $\lambda_{\rm face}$; plus the Darcy sub-cache
   $P,\varphi,\kappa,\mu$ if transport is `DarcyFick`).
3. `_compute_physics_rhs!` — dispatch on `transport_mode(def.mode)` to either
   `Discretization.compute_rhs!` (Fickian) or
   `Discretization.compute_rhs_darcy_fick!` (DarcyFick). Both populate
   `ws.dT_work[i]`, `ws.dξ_work[j,i]`, and `ws.rhs_cache.dξ_rxn[j,i]`.
4. Pack $\boldsymbol{du}$: $du[i]=dT_i/dt$, $du[n+(j-1)n+i]=d\xi_{j,i}/dt$.

The discrete updates assembled in step 3 are exactly those derived in §7 and
§13 (repeated here for completeness):

$$
\frac{dT_i}{dt} = \frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i + S_{{\rm conv},i} + S_{{\rm gen},i} + Q_{{\rm rad},i} + Q_{{\rm rxn},i}\Big],
$$

$$
\frac{d\xi_{j,i}}{dt} = -(\nabla\cdot\mathbf{J}_j)_i + S_{j,{\rm rxn},i}\quad(\text{gas}),\qquad
\frac{d\xi_{j,i}}{dt} = S_{j,{\rm rxn},i}\quad(\text{solid/liquid}),
$$

with the face fluxes, harmonic-mean conductivity, smooth-gated Arrhenius rates,
and the radiation model all as specified in their respective chapters. The
reaction heat source carries the storage sign convention,
$Q_{{\rm rxn},i}=-\sum_r h_r\,r_{r}(T_i,\xi_i)$ (endothermic $h>0$ gives
$Q_{\rm rxn}<0$); the conduction divergence enters with a leading minus so that
outflow cools the cell.

> **Implementation note.** Eulerian path: `_residual!(…, ::Eulerian)`,
> src/Residual/residual.jl. The physics RHS dispatch and the
> 8-step `compute_rhs!`/`compute_rhs_darcy_fick!` assembly are in
> src/Discretization/finite_volume.jl (see §7.8 for the strict step ordering).

### 15.3.2 Order dependencies inside the physics RHS

The 8-step `compute_rhs!` sequence (conductive fluxes → reactions → gas
transport → thermal/mass BCs → radiation → $S_{\rm conv}$ → divergences /
$S_{\rm gen}$ → final RHS) is *ordered* for correctness, not convenience:

- **Reactions before gas BCs.** The boundary-temperature Newton solve (step 4)
  needs the surface gas mass flux $\dot{m}_g$ (for transpiration when enabled),
  which requires the interior gas fluxes (step 3), which in turn need
  `dξ_rxn` populated (step 2).
- **Gas fluxes before thermal BCs.** Same chain: mass BCs (step 3) must finish
  before thermal BCs (step 4) so the surface enthalpy term has its $\dot{m}_g$.
- **Reaction rates before mesh velocity (ALE).** The ALE mesh velocity is built
  from the per-cell strain rate $\theta_i$, which is a function of the species
  reaction sources `dξ_rxn`; hence the physics RHS (which fills `dξ_rxn`) runs
  before the mesh-velocity integration.

---

## 15.4 The full ALE residual assembly

When `mesh_mode(def.mode) isa ALE`, mesh motion is part of the coupled system,
and the residual additionally evolves the node positions $z_k$ and (under
`WithChi`) the pyrolysis-progress $\chi_i$. The assembly is a twelve-step
sequence; we present it as the canonical pipeline and give the governing
equation for each step. Cross-references: the ALE kinematics and GCL are in §11;
the volume-change / $\chi$ model is in §10.

1. **Populate state cache.** Unpack $T,\xi,z,\chi$ from $\boldsymbol{u}$ into the
   typed cache.

2. **Update $\chi$-driven cross-section** (`update_cross_section_area!`,
   dispatched on `chi_mode`). Under `WithChi`, compute the mass-weighted
   column-average progress and set the column area:
   $$
   \bar\chi = \frac{\sum_i \chi_i\,m_{{\rm dry},i}}{\sum_i m_{{\rm dry},i}},\qquad
   A = A_0\,\mathrm{law}(\bar\chi),
   $$
   storing $A$ in `cache.A_section[]` carrying `eltype(u)` (AD-safe). The
   numerator accumulates `Dual` types; the denominator $M_{\rm dry}$ is a
   `Float64` constant stored on the mesh. Under `NoChi` this is a no-op
   returning $A_0$.

3. **Update state-derived geometry** (`update_state_geometry!(…, ALE())`).
   Recompute face positions $z_{f}=z_{{\rm node}(f)}$, face areas $A_f=A$, cell
   centers $z_i^c=\tfrac12(z_L+z_R)$, and cell volumes $V_i=|z_R-z_L|\,A$, all
   from `cache.z` and `cache.A_section`, preserving the scalar type.

4. **Update properties.** As in the Eulerian path, but now on the moving
   geometry.

5. **Physics RHS.** Same `compute_rhs!`/`compute_rhs_darcy_fick!` as Eulerian;
   importantly it also fills `ws.rhs_cache.dξ_rxn`, which the ALE steps reuse.

6. **Column-global radial strain and cumulative mesh velocity.** First compute
   $\theta_A$ (the column-global radial channel; $0$ under `NoChi`), then
   integrate the axial strain bottom-up to obtain the node velocities. The
   per-cell volumetric strain rate is the volume-change law of §10,
   $$
   \theta_i = \frac{1}{V_i}\frac{dV_i}{dt} = \sum_{j}\gamma_j\,\frac{1}{\rho_j}\,\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn},
   $$
   (only $\gamma_j>0$ components contribute; gases with $\gamma=0$ escape without
   changing matrix volume). The axial part removes the radial channel,
   $\theta_{L,i}=\theta_i-\theta_A$, and the cumulative integration with pinned
   bottom node $w_1=0$ gives
   $$
   w_{k+1} = w_k + (\theta_i-\theta_A)\,\Delta z_i,\qquad \Delta z_i = \frac{V_i}{A},
   $$
   i.e. $w_k=\sum_{i<k}\theta_{L,i}\,\Delta z_i$ (an exclusive cumulative sum;
   this is exactly the structure represented by the `PrefixSumCoupling` in the
   Jacobian, §15.9). The radial channel under `WithChi` is, by the chain rule,
   $$
   \theta_A = \frac{1}{A}\frac{dA}{dt} = \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\frac{d\bar\chi}{dt},\qquad
   \frac{d\bar\chi}{dt} = \frac{1}{M_{\rm dry}}\sum_i S_{{\rm dry},i}\,V_i,
   $$
   with $\mathrm{law}'(\bar\chi)$ obtained by `ForwardDiff.derivative(law, χ̄)`
   and $S_{{\rm dry},i}=\sum_{j\in\text{dry-gas}}(d\xi_{j,i}/dt)|_{\rm rxn}$. The
   guard `law(χ̄) > ε` prevents the radial channel from diverging as the
   cross-section collapses.

7. **Interpolate node velocities to faces** (`_interpolate_node_to_face!`): in
   1D the face sits at a node, and the cell-centered velocity used by advection
   is $w_{{\rm cell},i}=\tfrac12(w_L+w_R)$.

8. **ALE advection** (`compute_ale_advection_term!`,
   `compute_ale_species_advection!`). The advection correction is driven by the
   mesh velocity $w$ alone — the gas Darcy flux is already in the lab-frame RHS
   from step 4, and the $+w\,\partial_z$ correction converts it to the moving
   frame (see §11.6, §13.8).

9. **Apply ALE corrections** to $dT$ and $d\xi$ (`_apply_ale_advection!`). The
   ALE chain-rule identity (§11.2) adds $w\cdot\nabla f$ for every field whose
   lab-frame rate the Eulerian RHS already computes — $T$, and gas $\xi_g$
   (whose full lab-frame transport flux entered in step 4, so the same
   $+w\cdot\nabla$ frame conversion applies); solids are mesh-bound and get no
   advection (their transport *is* the mesh motion, handled by dilation in
   step 10):
   $$
   \frac{dT_i}{dt}\bigg|_{\rm ALE} \mathrel{+}= w_{{\rm cell},i}\Big(\frac{\partial T}{\partial z}\Big)_i,\qquad
   \frac{d\xi_{g,i}}{dt}\bigg|_{\rm ALE} \mathrel{+}= w_{{\rm cell},i}\Big(\frac{\partial \xi_g}{\partial z}\Big)_i,
   $$
   with the gradients computed by a TVD-limited 3-point donor-side stencil (Van
   Leer default), upwinded on the sign of the apparent transport velocity
   $c=-w$ (§11.6.1, §13.8).

10. **Dilation correction** for mesh-bound (non-gas) species
    (`_apply_ale_dilation_correction!`):
    $$
    \frac{d\xi_{j,i}}{dt}\bigg|_{\rm dilation} \mathrel{-}= \xi_{j,i}\,\theta_i\quad(j\ \text{non-gas}),
    $$
    which makes the solids Lagrangian (their *mass thickness* is locked to the
    mesh; concentration rises as the cell contracts). Gases are excluded because
    they flow through the pores and are not mesh-bound.

11. **Pack the $T,\xi,z$ blocks.** $du[i]=dT_i/dt$ (with the ALE correction),
    $du[n+(j-1)n+i]=d\xi_{j,i}/dt$ (with advection + dilation), and
    $du[n(1+N_C)+k]=w_k$ — i.e. $dz_k/dt=w_k$ (Euler-forward node motion
    realized by the integrator).

12. **Pack the $\chi$ block** (`_pack_chi_derivative_dispatched!`, `WithChi`
    only). The progress variable obeys
    $$
    \frac{d\chi_i}{dt} = \frac{S_{{\rm dry},i}\,V_i}{m_{{\rm dry},i}},\qquad
    S_{{\rm dry},i}=\sum_{j\in\text{dry-gas}}\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn},
    $$
    with a `max(m_init, eps())` guard against division by zero. Note $\chi$ is
    *reconstructed* from `dξ_rxn` each residual rather than carried as an
    independent integral — it is still a genuine ODE state (the integrator
    advances it), but its derivative is recomputed from first principles every
    evaluation, keeping it exactly in sync with cumulative dry-gas production.
    Under identity-law materials the $\chi$ derivatives are zeroed.

The complete ALE discrete updates therefore read

$$
\frac{dz_k}{dt}=w_k,
$$
$$
\frac{dT_i}{dt}=\frac{-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+S_{{\rm gen},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}}{\rho c_p^{\rm eff}{}_{,i}} + w_{{\rm cell},i}\Big(\frac{\partial T}{\partial z}\Big)_i,
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i} + \underbrace{w_{{\rm cell},i}\Big(\tfrac{\partial\xi_j}{\partial z}\Big)_i}_{\text{advection (gas only)}} - \underbrace{\xi_{j,i}\,\theta_i}_{\text{dilation (non-gas)}}.
$$

> **Implementation note.** ALE path: `_residual!(…, ::ALE)`,
> `_compute_ale_mesh_velocity!`, `_compute_theta_A` / `_compute_theta_A_dispatched`,
> `_apply_ale_advection!`, `_apply_ale_dilation_correction!`, and
> `_pack_chi_derivative!` are all in src/Residual/residual.jl. The strain-rate
> law `volume_change_rate` is in src/Physics/volume_change.jl. The ALE flux
> kernels are in src/Discretization/ale_fluxes.jl.

### 15.4.1 Conservation-tracker hook and energy forms

The two energy forms (`:standard` vs `:with_generation_sink`) enter at step 5
inside `compute_rhs!`/`compute_rhs_darcy_fick!` exactly as in the Eulerian path:
`:with_generation_sink` resolves a single volumetric reference temperature
$T_{\rm ref}$ from the boundary ambient and adds $S_{\rm gen}$ using the
midpoint-rule sensible enthalpy $\Delta h_g$, while `:standard` zeroes
$S_{\rm gen}$. The diagnostics layer (§16) attaches a per-substep accumulator
that audits the matrix and total energy ledgers; this is a *callback*, not part
of the residual, so it does not perturb the integration. The conservation
audit is what quantifies the quasi-steady-gas cost (the gas-storage gap, $<1\%$
at 1 atm).

---

## 15.5 Time integration

### 15.5.1 The method-of-lines ODE and the integrator

The packed residual turns the PDE system into the autonomous-in-form ODE
$\boldsymbol{du}/dt=\mathrm{RHS}(\boldsymbol{u},t)$, which `solve` hands to
OrdinaryDiffEq.jl. Pyrolysis kinetics are *stiff* (Arrhenius source terms span
many decades over a narrow temperature window; conduction adds the usual
parabolic stiffness $\sim\Delta t/\Delta z^2$), so the default integrator is an
L-stable singly-diagonally-implicit Runge–Kutta (SDIRK) method:

$$
\texttt{KenCarp4}\big(\text{autodiff}=\text{AutoFiniteDiff()},\ \text{linsolve}=\text{SparspakFactorization()}\big).
$$

KenCarp4 is a 4th-order ESDIRK pair with an embedded error estimator for
adaptive step-size control; its singly-diagonal structure means every implicit
stage shares the same Newton iteration matrix, which we exploit below. The
default linear solver is `SparspakFactorization` (reported $\approx28\%$ faster
than KLU on the block-tridiagonal 1D pattern); KLU and the ILU-Krylov solvers of
§15.12 are drop-in alternatives.

> **Implementation note.** Integrator resolution: `_resolve_integrator`,
> src/Solver/solve.jl. A user may pass a constructor (`integrator=Rodas5`) plus
> a `linear_solver=`, or a fully-built algorithm. The RHS closure re-samples
> the incident radiation at every evaluation (`ws.I_radiation = radiation_flux(t)`)
> whenever the active model consumes it (`SURFACE_ABSORPTION` / `BEER_LAMBERT`)
> — constant and time-varying fluxes alike (§8.3.4). A one-time sample at
> $t=0$ additionally seeds the structured-Jacobian cache's
> `radiation_intensity`.

### 15.5.2 The implicit-RK $W$ matrix and warm start

Each implicit SDIRK stage solves a nonlinear algebraic system for the stage
increment by Newton's method. The Newton iteration matrix — OrdinaryDiffEq's
"$W$" matrix — is

$$
\boxed{\,W = \frac{1}{\gamma_{\rm RK}\,\Delta t}\,I - J\,,}\qquad J = \frac{\partial\,\mathrm{RHS}}{\partial\boldsymbol{u}},
$$

where $\gamma_{\rm RK}$ is the (common, by the singly-diagonal property) stage
coefficient and $J$ is the residual Jacobian. (Equivalently $W$ may be scaled by
$\gamma_{\rm RK}\Delta t$ to $W=I-\gamma_{\rm RK}\Delta t\,J$; the scaling is an
internal detail of the integrator. The factor that matters for our purposes is
that $W$ is a *diagonal shift of $-J$*, so the diagonal of $J$ adds to the
$1/(\gamma_{\rm RK}\Delta t)$ scale, which is what makes the structured Jacobian's
sparse block $\mathsf{A}$ diagonally dominant and the Richardson iteration of
§15.10 converge in a few steps.) Each Newton step solves $W\,\delta\boldsymbol{u}=-r$
for the stage residual $r$; the linear solve is delegated to the chosen
`LinearSolveStrategy` (DirectSolve / StructuredSolve / ApproxSolve /
MatrixFreeSolve, §15.10) or, when no pluggable Jacobian is supplied, to the
integrator's colored finite-difference Jacobian (`AutoFiniteDiff`) with the
configured `linsolve`.

The Jacobian $J$ is supplied to the integrator through an `ODEFunction` carrying
both `jac = jac_fn!` and a `jac_prototype` sparsity pattern (the union of the
physics and geometry patterns of §15.11):

```julia
jac_fn! = (J, u, p, t) -> compute_jacobian!(J, jac_spec, u, def, t, ws)
ode_fn  = ODEFunction(rhs!; jac = jac_fn!, jac_prototype = Float64.(S_phys .| S_geom))
```

**Warm start.** A prior solution's final state can seed a new solve's initial
condition (`warm_start=prev_solution`), provided the component count, mesh size,
ALE mode, and $\chi$ mode all agree (enforced implicitly by a state-vector
length check). The integrator's `tspan[1]` still defines the start time, so
warm-starting treats the prior final state as fresh initial data — useful for
continuation studies and for sensitivity sweeps, with the caveat (§17.7) that
warm-start artifacts can perturb finite-difference sensitivities.

> **Implementation note.** ODE construction, Jacobian wiring, and warm-start
> application: src/Solver/solve.jl (`solve`, `_apply_warm_start!`). The
> per-step linear-solve protocol is OrdinaryDiffEq's; Pyrolysis supplies only
> $J$ and the strategy.

### 15.5.3 Output-time handling and callbacks

Snapshots are saved at `saveat=problem.output_times` (11–101 evenly spaced
samples by default). Callbacks are assembled in a fixed order and combined into
a `CallbackSet`:

| Order | Callback | Type | Role |
|---|---|---|---|
| 1 | `StepUpdateCallback` | discrete (every accepted step) | mirror integrator $\boldsymbol{u}\to$ mutable mesh; sync geometry & cross-section |
| 2 | `ThicknessTerminationCallback` | continuous (ALE) | terminate when $z_{\rm top}-z_{\rm bot}<$ `min_thickness` |
| 3 | `DepletionCallback` | discrete (ALE, opt-in) | merge cells below their per-cell thickness threshold; resize state |
| 4 | `ProgressCallback` | periodic | progress bar |
| 5 | `DiagnosticsCallback` | discrete (every `diagnostics_stride`) | snapshot energy/mass ledgers |
| 6 | `VerifyJacobianCallback` | discrete (at sample times) | compare active Jacobian vs. dense AD reference |
| 7 | user callback | — | custom |

The **`StepUpdateCallback`** is the linchpin of the read-only-mesh design: after
each accepted step it sets `mesh.cross_section_area = cross_section_area_from_state(…)`
and calls `update_mesh_from_state!(mesh, u, layout)`, so all between-step
consumers (AMR, depletion detection, diagnostics, the surface BC solve in
post-processing) see a mesh consistent with the integrator state.

The **`DepletionCallback`** (§14.3) detects cells whose current thickness has
fallen below `depletion_threshold × Δz_init,i` (default $5\%$ of each cell's
own initial thickness), merges them into a
neighbor with the condensed-enthalpy-conserving temperature merge of §14.3.2

$$
\sum_{j\notin{\rm gas}} m_j\!\int_{T_{\rm ref}}^{T_{\rm merged}}\! c_{p,j}\,dT'
= H_a(T_a) + H_b(T_b)
\ \text{(bracketed Newton; gas heat capacity excluded)},\qquad
\xi_{{\rm merged},j} = \frac{V_a\,\xi_{a,j}+V_b\,\xi_{b,j}}{V_a+V_b},
$$

then resizes the workspace caches and the integrator state vector
(`resize!(integrator, new_len)`; `u_modified!(integrator, true)`). Because the
state vector changes length mid-solve, depletion is **incompatible with a
pluggable `JacobianSpec`** (whose caches are sized at setup) — when
`handle_depletion=true`, `solve` sets `jac_spec=nothing` and falls back to the
integrator's colored finite-difference Jacobian (§15.13).

> **Implementation note.** Callbacks: src/Solver/callbacks.jl
> (`create_step_update_callback`, `create_thickness_termination_callback`,
> `create_depletion_callback`, `create_merged_state_vector`) and src/Solver/solve.jl
> (assembly, diagnostics, verify-Jacobian). The depletion↔pluggable-Jacobian
> guard is in `_validate_solve_kwargs`, src/Solver/solve.jl.

---

## 15.6 Jacobian backends

Pyrolysis.jl exposes exactly two production backend families, unified behind a
single user-facing handle `JacobianSpec`:

1. **`DenseAD_ForwardDiff`** — runs `ForwardDiff` over the *entire* residual
   vector and returns a dense Jacobian. It is correct by construction (it
   differentiates the same `residual!` the integrator uses) but $O(N^2)$ in
   storage and slow; it is the **verification anchor**, not a production
   backend.

2. **`Structured{S}`** — the production backend, parameterized by a
   `LinearSolveStrategy` $S$ (DirectSolve / StructuredSolve / ApproxSolve /
   MatrixFreeSolve). It assembles the Jacobian *operator by operator*, exploiting
   the block-tridiagonal physics structure and representing the non-local terms
   as O(N) couplings rather than dense matrix entries.

The auto-default when the user passes no `jacobian` kwarg is
`Structured(strategy = DirectSolve())` — exact assembly plus a sparse direct
solve. This is the production baseline; large ALE problems can opt into
`ApproxSolve` for speed (§15.10).

`JacobianSpec` bundles the backend with everything the integrator needs:

```julia
struct JacobianSpec{B<:JacobianBackend, S, C}
    backend::B
    sparsity_phys::S      # SparseMatrixCSC{Bool,Int} or nothing
    sparsity_geom::S      # geometry (z-column) pattern; nothing for Eulerian
    cache::C              # backend-private workspace (not shareable)
    colorvec_phys::Vector{Int}   # column coloring for FD/sparse-AD
    colorvec_geom::Vector{Int}
end
```

A single entry point `compute_jacobian!(J, spec, u, def, t, ws)` dispatches to
`_compute_jacobian!(J, spec.backend, …)`.

> **Implementation note.** Backend contract and `JacobianSpec`:
> src/Jacobian/interface.jl. Dense reference: src/Jacobian/dense_ad_reference.jl.
> Structured backend assembly: src/Jacobian/structured_backend.jl. The
> `is_approximate(backend)` trait (default `false`; `true` for
> `Structured{ApproxSolve}`) selects strict vs. masked verification (§15.13).

---

## 15.7 The structured-operator architecture

The structured Jacobian is organized into three orthogonal layers.

### 15.7.1 Layer 1 — Operators

Each physics block is an `Operator` subtype: a lightweight stateless tag that
owns *both* its residual contribution (`apply!`) and its Jacobian contribution
(`contribute_to_jacobian!`, optionally `contribute_to_structured!`,
`contribute_to_implicit!`). The default `contribute_to_jacobian!` differentiates
the operator's own `apply!` with `ForwardDiff` over its declared stencil, so an
operator is *correct by construction*; analytical overrides are performance
optimizations, not correctness requirements. Each operator also declares:

- `reads(op)` — derived fields it consumes (e.g. `:w_cell`, `:theta_A`),
- `writes(op)` — derived fields it produces,
- `stencil(op)` — its cell/face neighborhood (e.g. `(cells = -1:1)` for
  conduction).

The fixed assembly walk is

```
1. ReactionSourceOperator
2. ConductionOperator
3. GasTransportOperator | DarcyTransportOperator
4. RadiationOperator
5. BCThermalOperator
6. BCMassOperator
7. BCConvectionOperator                    (boundary-cell S_conv energy row)
8. (ALE)      ALEMeshVelocityOperator      → contribute_to_structured!
9. (ALE)      ALEAdvectionOperator
10.(WithChi)  ChiOperator
11.(ALE & !ale_reduced) GeometryOperator   (full-residual AD over z-columns)
12.(WithChi)  ChiOperator                  → contribute_to_structured! (rank-1)
```

This order guarantees that providers run before consumers: properties are
updated first; the ALE mesh velocity is computed before advection reads it; and
the $\chi$ derivatives are available to the column-global rank-1 coupling.

The dependency order is made explicit by the **`ResidualContext`**, which
topologically sorts operators by their `reads`/`writes` declarations and holds a
per-step `derived_fields` map. A provider pushes a *bare* `StructuredCoupling`
primitive into `derived_fields[:w_cell]`; a consumer reads it and `compose_into!`
wraps it at its own row/column offsets. Operators thus only know their *local*
math; the framework handles the global embedding.

> **Implementation note.** Operator framework: src/Jacobian/operators/Operator.jl;
> dependency graph: src/Jacobian/residual_context.jl; per-operator
> implementations under src/Jacobian/operators/. Per-operator `ForwardDiff` chunk
> sizes: reaction `1+N_C`; conduction $3(1+N_C)$; gas transport $3(1+N_C)$ (or 5
> with the limiter window); BC thermal $1+N_C$; BC mass $2+N_C$ (the energy row
> reads the interior neighbor's $T$); ALE advection $5(1+N_C)$;
> $\chi$ $1+N_C$.

### 15.7.2 Layer 2 — `StructuredJacobian`

The assembled Jacobian is *not* a monolithic sparse matrix. It is

```julia
struct StructuredJacobian{Mode}
    A::SparseMatrixCSC{Float64,Int}        # local block: block-tridiag physics + sparse augments
    structured::Vector{StructuredCoupling} # typed O(N)-apply non-local couplings
    implicit::Vector{ImplicitTerm}         # IFT-resolved couplings (rarely used)
    smap::FastSparsityMap                  # O(1) sparse-write helper
    n::Int                                 # full state dimension
end
```

The split is the architectural heart of the design: the sparse block
$\mathsf{A}$ carries everything *local* (block-tridiagonal mass / conduction /
reaction / transport plus the cell-local self-terms), and the `structured`
vector carries everything *non-local* (the cumulative ALE mesh velocity, the
Beer–Lambert optical-depth suffix, the column-global $\theta_A$ and $\chi$
back-reactions) as operators with O(N) apply methods. A solve never has to
materialize the dense non-local blocks unless the chosen strategy is
`DirectSolve`.

`FastSparsityMap` is an O(1) sparse-write helper: a dense lookup table
`indices[row + (col-1)*n_rows] = nz_idx` (zero meaning "not in pattern"). It
trades quadratic memory ($\sim2.5\,\mathrm{MB}$ at 560 DOFs, $\sim90\,\mathrm{MB}$ at
3360 DOFs) for O(1) scattered writes via `fast_sparse_write!` /
`fast_sparse_add!` / `fast_sparse_sub!`.

### 15.7.3 Layer 3 — `LinearSolveStrategy`

A strategy is the rule for inverting $\mathsf{A}+\sum_j C_j$. It is orthogonal to
the assembling backend; at solve time the integrator calls
`ldiv!(out, strategy, jac, b)`. The four strategies are detailed in §15.10.

---

## 15.8 Per-operator Jacobian formulas

We give the exact discrete derivative each operator contributes, with the sign
conventions verified against the code. Throughout, the residual writes
$dT_i\mathrel{+}=(\text{source})/\rho c_p^{\rm eff}$ and the conduction
divergence enters with a leading minus, so signs are tracked *per operator* and
accumulated additively.

### 15.8.1 Reaction source

The cell-local contribution is

$$
\frac{dT_i}{dt}\bigg|_{\rm rxn}=\frac{Q_{\rm rxn}(T_i,\xi_i)}{\rho c_p^{\rm eff}(T_i,\xi_i)},\qquad
\frac{d\xi_{j,i}}{dt}\bigg|_{\rm rxn}=S_{j,{\rm rxn}}(T_i,\xi_i),
$$
$$
Q_{\rm rxn}=-\sum_r h_r\,r_r,\qquad S_{j,{\rm rxn}}=\sum_r \nu_{r,j}\,r_r.
$$

The Jacobian is purely cell-local — the $(1+N_C)\times(1+N_C)$ block coupling
$(T_i,\xi_{1,i},\dots,\xi_{N_C,i})$ to itself — obtained by `ForwardDiff` over
`apply_cell!` (which evaluates the reactions, the heat source *including the
$1/\rho c_p$ division*, and the species sources in place). Note the energy row
differentiates the *quotient* $Q_{\rm rxn}/\rho c_p^{\rm eff}$, so it captures
both the explicit $T,\xi$ dependence of the rates *and* the $\xi$-dependence of
the heat capacity. The smooth tanh temperature gates and the multiplicative
$\tanh(\xi/\text{threshold})$ depletion rate factor (§6) make these
derivatives well-defined and free of the spurious spikes a hard Heaviside or
$\min$ would produce — this is the reason the kinetics are AD-clean.

> **Implementation note.** `ReactionSourceOperator`, `apply_cell!`,
> src/Jacobian/operators/reaction_source.jl. Chunk size $1+N_C$.

### 15.8.2 Conduction (tridiagonal divergence)

The energy row depends on the three cells $\{i-1,i,i+1\}$ through the divergence
of the harmonic-mean face flux. With the distance-weighted harmonic-mean
conductivity (§7.1)

$$
q_f=-k_f\frac{T_R-T_L}{d_L+d_R},\qquad k_f=\frac{d_L+d_R}{d_L/k_L+d_R/k_R},
$$

the operator computes (verified in src/Jacobian/operators/conduction.jl)

$$
q_{\rm right}=-k_{f,R}\frac{T_R-T_C}{d_{CR}},\quad
q_{\rm left}=-k_{f,L}\frac{T_C-T_L}{d_{LC}},
$$
$$
\frac{dT_i}{dt}\bigg|_{\rm cond}=-\frac{(q_{\rm right}A_{\rm right}-q_{\rm left}A_{\rm left})/V_C}{\rho c_p^{\rm eff}(T_C,\xi_C)}.
$$

The leading minus realizes "outflow on the right cools the cell," matching the
discrete divergence $\mathrm{div}[i]=(q_R-q_L)A/V_i$ and the residual's
$-\mathrm{div}[i]/\rho c_p$. In ALE with $\chi$-shrinkage, $A_{\rm left}$ and
$A_{\rm right}$ can differ slightly (the face areas carry the geometry), so the
code keeps them distinct rather than factoring out a single $A$. The Jacobian is
the tridiagonal block from `ForwardDiff` over the 3-cell stencil with chunk size
$3(1+N_C)$: each of the three cells contributes $(T,\xi_1,\dots,\xi_{N_C})$
inputs because $k_f$ and $\rho c_p^{\rm eff}$ both depend on composition.

**Boundary-cell self-fill.** At a boundary cell (cell 1 has no left neighbor,
cell $n$ no right), the missing neighbor's $(T,\xi)$ in the AD stencil is filled
with the cell's *own* value, and the missing face flux is zero (the BC operator
owns the boundary flux). AD over a zero face produces structural-zero columns,
which the sparse scatter respects. The genuine boundary flux contribution is the
`BCThermalOperator`'s responsibility, computed through the surface-temperature
Newton solve (§12.2).

> **Implementation note.** `ConductionOperator`, `apply_cell!`, `_pack_stencil!`,
> src/Jacobian/operators/conduction.jl.

### 15.8.3 Gas transport / Darcy

The Fickian gas-transport operator differentiates the volume-fraction-gradient
flux (§9.1) over a $\le\pm1$-cell stencil (chunk $3(1+N_C)$), giving the direct
derivatives $\partial J/\partial\lambda_f$, $\partial J/\partial\xi_L$,
$\partial J/\partial\xi_R$ and their $T$-coupling through $T_{\rm face}$. The
Darcy operator (`DarcyTransportOperator`) additionally carries the
pressure-driven advection $(\xi_{\rm upwind}/\varphi_{\rm face})\,v_g$: it
reconstructs the face velocity inside its own stencil closure from the
cell-local state — $P_i$ from the ideal-gas closure, $\kappa_i$ from the mixing
rule, $\mu_i$ from Sutherland, the face mobility by the residual's
distance-weighted harmonic rule (`_dt_mobility_face`, zero when either side
$\le10^{-30}$) — so the chain $\partial J/\partial\xi$ and
$\partial J/\partial T$ through $P$, $\varphi$, $\kappa$, $\mu$ is
differentiated through the same formulas the residual executes (including the
$\xi/\varphi_{\rm face}^2$ porosity chain, zero when the $\varphi$ floor binds).
First-order upwinding keeps the advective stencil local; the diffusive part uses
harmonic-mean $\lambda_f$ as in the residual. Boundary faces are zeroed here —
the boundary venting flux is `BCMassOperator`'s row, and the boundary cell's
whole $S_{\rm conv}$ energy row (whose gradient runs through the surface-Newton
$T_s$) is `BCConvectionOperator`'s (§15.8.8); the transport operators keep the
boundary cell's ξ-rows and $S_{\rm gen}$ interior-face slice but contribute no
boundary-cell $S_{\rm conv}$. The hand-coded
analytical override is the default and is gated against the local-AD fallback at
`rtol = 10^{-12}`.

> **Implementation note.** `GasTransportOperator`,
> src/Jacobian/operators/gas_transport.jl; `DarcyTransportOperator`,
> src/Jacobian/operators/darcy_transport.jl (`_dt_props` Dual pass,
> `_dt_mobility_face`, `_dt_J_face_partials`, analytical override
> `_contribute_to_jacobian_analytical!`).

### 15.8.4 ALE mesh velocity and the solid dilation correction

The mesh velocity operator does **not** AD a single residual; it builds the
analytical chain-rule pieces and emits a structured coupling. Two derivative
sweeps populate workspace buffers:

- `compute_strain_rate_derivatives!` — cell-local `ForwardDiff` over
  $\theta(T,\xi)=\big(\sum_j\gamma_j S_j(T,\xi)\big)/\rho(T)$, producing
  $d\theta/dT_i$ and $d\theta/d\xi_{j,i}$.
- `compute_theta_A_derivatives!` (WithChi) — analytical $d\theta_A/d\chi_j$ and
  $dA/d\chi_j$.

The cumulative velocity $w_k=\sum_{i<k}\theta_{L,i}\Delta z_i$ has the
*exclusive prefix-sum* structure: $\partial w_k/\partial x_i=\Delta z_i\,
\partial\theta_i/\partial x_i$ for $i<k$ (and $0$ for $i\ge k$). The per-cell
scaling $v[i]=\Delta z_i\,\partial\theta_i/\partial x_i$ is captured as a bare
`PrefixSumCoupling{Float64}` per source channel ($T$ and each $\xi_j$) and
consumed twice:

- pushed into `derived_fields[:w_cell]` for the advection consumer, which
  wraps it with its local row scaling (§15.9);
- embedded **at the z-row offsets** as `EmbeddedNodePrefixSumCoupling` blocks —
  the z-row equations $dz_k/dt = w_k$ are themselves consumers of the
  cumulative coupling, with $\partial(dz_k/dt)/\partial x_j =
  \mathbb{1}[j\le k-1]\,v[j]$. The geometry AD pass (§15.11.4) covers z
  *columns* only and plays no part in these rows; the z-row × $(T,\xi)$ block
  is owned entirely by these embedded couplings.

Under `WithChi`, the column-global $\theta_A$ enters every
$w_k = \sum_{i<k}(\theta_i-\theta_A)\Delta z_i$, so each source channel
$y\in\{\chi,\,T,\,\xi_1..\xi_{N_C}\}$ contributes a genuinely rank-1 block
$\partial w_k/\partial y_j = -\text{cumZ}_{\rm node}[k]\cdot
\partial\theta_A/\partial y_j$ with

$$
\frac{\partial\theta_A}{\partial\chi_j} = G'\,\frac{m_{{\rm init},j}}{M}\,
\frac{d\bar\chi}{dt} + G\,\frac{\partial A/\partial\chi_j}{M}\,
\frac{\sum_i S_{{\rm dry},i}V_i}{A},\qquad
\frac{\partial\theta_A}{\partial T_j} = \frac{G\,V_j}{M}\,
\frac{\partial S_{{\rm dry},j}}{\partial T_j},\qquad
\frac{\partial\theta_A}{\partial \xi_{k,j}} = \frac{G\,V_j}{M}\,
\frac{\partial S_{{\rm dry},j}}{\partial \xi_{k,j}},
$$

where $G=\mathrm{law}'(\bar\chi)/\mathrm{law}(\bar\chi)$. The bare rank-1s
(shared left factor $u[k]=-\text{cumZ}_{\rm node}[k]$) are pushed to
`derived_fields[:w_node_chi]` ($\chi$ channel) and
`derived_fields[:w_node_theta_A]` ($T,\xi$ channels) for the advection
consumer, and embedded directly at the z-rows as `EmbeddedRank1Coupling`
blocks.

The **solid dilation correction** $d\xi_{j,i}\mathrel{-}=\xi_{j,i}\theta_i$ adds
the cell-local Jacobian (non-gas $j$ only)

$$
\frac{\partial}{\partial T_i}:\ -\xi_{j,i}\frac{\partial\theta_i}{\partial T_i},\qquad
\frac{\partial}{\partial\xi_{k,i}}:\ \begin{cases}
-\theta_i-\xi_{j,i}\dfrac{\partial\theta_i}{\partial\xi_{j,i}} & k=j,\\[6pt]
-\xi_{j,i}\dfrac{\partial\theta_i}{\partial\xi_{k,i}} & k\ne j,
\end{cases}
$$

reusing the cached $d\theta/dT$, $d\theta/d\xi$ from the strain-rate sweep
(`_apply_ale_dilation_jacobian!`).

> **Implementation note.** `ALEMeshVelocityOperator`,
> `compute_strain_rate_derivatives!`, `compute_theta_A_derivatives!`,
> src/Jacobian/operators/ale_mesh_velocity.jl; dilation Jacobian in
> src/Jacobian/structured_backend.jl. A documented sharp edge: the
> $d\theta/dT$, $d\theta/d\xi$ buffers persist in the workspace and are reused by
> both the prefix-sum builder and the dilation Jacobian; assembly is atomic per
> Newton step, so the reuse is safe in practice.

### 15.8.5 ALE advection (TVD)

The advection energy row depends on the 5-cell TVD window
$[i-2,i+2]$ through the limited gradient and on *all upstream cells* through the
mesh velocity. The direct part is the chain rule
$\partial dT_i/\partial T_j$ for $j$ in the limiter window; the indirect part is
$\partial(\,(\partial T/\partial z)_i\,w_{{\rm cell},i})/\partial x_k$ for all
$k\le i$, realized by composing the provider's prefix-sum coupling with the local
row scaling $\text{row\_scale}[i]=(\partial T/\partial z)_i$ (and
$(\partial\xi_j/\partial z)_i$ for gas species). Chunk size $5(1+N_C)$.

Under the `:ale_reduced` profile (used by `ApproxSolve`), the advection operator
**skips its `compose_into!` calls entirely** — it emits no structured couplings.
This is the defining feature of the reduced profile: it is not merely a sparsity
choice, it changes which operators contribute to `jac.structured`.

> **Implementation note.** `ALEAdvectionOperator`,
> src/Jacobian/operators/ale_advection.jl.

### 15.8.6 Variable cross-section ($\chi$)

Under `WithChi` (and only for non-identity laws; identity-law materials have
$d\chi/dt\equiv0$), the $\chi$ row is

$$
\frac{d\chi_i}{dt}=\frac{S_{{\rm dry},i}(T_i,\xi_i)\,V_i}{m_{{\rm init},i}}.
$$

Its Jacobian has three parts:

- **Cell-local $(T,\xi)$**: `ForwardDiff` over $S_{\rm dry}\,V/m_{\rm init}$ gives
  $\partial(d\chi_i/dt)/\partial T_i$ and $\partial(d\chi_i/dt)/\partial\xi_{j,i}$.
- **$z$-columns (ALE)**: because $V_i=A|z_{i+1}-z_i|$ is linear in the two
  bounding node positions, the analytical entries are
  $\partial(d\chi_i/dt)/\partial z_i=-A\,S_{\rm dry}/m_{\rm init}$ and
  $\partial(d\chi_i/dt)/\partial z_{i+1}=+A\,S_{\rm dry}/m_{\rm init}$.
- **$\chi$-columns (rank-1)**: the area $A=A_0\,\mathrm{law}(\bar\chi)$ depends on
  *every* $\chi_j$ through $\bar\chi$, giving a dense column coupling factored as
  an `EmbeddedRank1Coupling` with
  $u[i]=(S_{{\rm dry},i}/m_{{\rm init},i})\,\text{cumZ}[i]$ and $v[j]=dA/d\chi_j$.

> **Implementation note.** `ChiOperator`, `apply_cell!`,
> src/Jacobian/operators/chi.jl. Chunk size $1+N_C$. The $\chi$ block is appended
> after $(T,\xi,z)$; the offset is $\chi_{\rm off}=n(1+N_C)+\delta_{\rm ALE}n_{\rm nodes}$.

### 15.8.7 Beer–Lambert non-local radiation

In `BeerLambert` mode the absorbed flux in cell $i$ depends on the optical depth
of *all cells above it*,
$q_{\rm rad}[i]=I_{\rm in}\,(1-e^{-\tau_i})/\Delta z_i$ with
$\tau_i=\alpha_i\Delta z_i$ and $I_{\rm in}$ the cumulative-attenuated intensity
from the top. Differentiating (verified in src/Jacobian/operators/radiation.jl):

$$
\frac{\partial q_{\rm rad}[i]}{\partial x_j}=
\begin{cases}
I_{\rm in}\,e^{-\tau_i}\,\dfrac{\partial\alpha_i}{\partial x_j}=\big(I_{\rm in}-q_{\rm rad}[i]\,\Delta z_i\big)\dfrac{\partial\alpha_i}{\partial x_j} & j=i,\\[10pt]
-\,q_{\rm rad}[i]\,\Delta z_j\,\dfrac{\partial\alpha_j}{\partial x_j} & j>i,\\[6pt]
0 & j<i.
\end{cases}
$$

The energy row $T[i]$ therefore couples to the $\xi$-columns of every cell
$j\ge i$ — an **upper-triangular** (strictly, a *suffix*) pattern. The diagonal
$j=i$ term is added analytically by the operator; the strict-suffix $j>i$
couplings are captured as an `EmbeddedScaledSuffixSumCoupling`. The split
($j=i$ analytical, $j>i$ structured) avoids double-counting the diagonal.

On the Kirchhoff default path (`radiation_kirchhoff = true`, §8.3.2) the
injected intensity is $\varepsilon_{\rm eff}(T_n,\xi_{\cdot,n})\,I_{\rm surface}$,
so every energy row $i<n$ additionally couples to the **surface cell**: the
suffix vector carries
$v_n \mathrel{-}= (\partial\varepsilon_{\rm eff}/\partial x_n)/\varepsilon_{\rm eff}$,
the surface cell's diagonal closure recomputes $\varepsilon_{\rm eff}(T,\xi)$
live, and the sparsity pattern gains $T$-row($i$) × $T$-column($n$) slots for
$i<n$ (`_emit_radiation_pattern!`). With an explicit BC absorptivity the factor
is folded into $I_{\rm surface}$ upstream and these terms vanish.

One deliberate approximation remains: under ALE the operator fills **T and ξ
columns only** — the dependence of the deposit on the live geometry
($\Delta z_i = V_i/A$ from the z block) is omitted from the Jacobian (the
residual is exact). This degrades only Newton convergence of the implicit
integrator, never the converged solution (src/Jacobian/operators/radiation.jl).

> **Implementation note.** `RadiationOperator`,
> src/Jacobian/operators/radiation.jl. Active only in `BeerLambert` mode;
> `SurfaceAbsorption`/`NoRadiation` leave a tridiagonal pattern.

### 15.8.8 Boundary conditions

`BCThermalOperator` overwrites the boundary-cell conductive flux with the
BC-specified surface flux (convective + radiative + fixed-$T$), obtained from the
surface-temperature Newton solve $f(T_s)=q_{\rm BC}(T_s)-k_{\rm eff}(T_s-T_{\rm cell})/\Delta z_{1/2}=0$
(§12.2). Its Jacobian is taken by AD through that Newton solve, producing the
boundary-cell row's dependence on $(T_{\rm cell},\xi_{\rm cell})$. The outward
mass flux fed to the Newton solve is **hard-zeroed**, mirroring the residual:
transpiration is rejected at construction (`to_problem_def` throws on
`use_transpiration_bc = true`) and the residual zeroes $\dot m_g$ before its
surface solve, so the operator must differentiate the same balance. Chunk size
$1+N_C$.

`BCMassOperator` owns the species rows of the boundary-face gas flux $F_g$
(film law plus, in Darcy mode, the pressure-BC advective venting),
$d\xi_{g,\rm bdy} \mathrel{+}= F_g A_f/V$, and the energy row's boundary-face
slice of the (opt-in) generation sink,

$$
\Delta S_{\rm gen} = +\sum_g \Delta h_g(T_{\rm bdy},T_{\rm ref})\,F_g\,A_f/V,
$$

because the residual's $\nabla\!\cdot\!J_{\rm gas}$ is formed *after*
`apply_mass_bc_flux!` fills the boundary faces, while the transport operators
zero those faces. $S_{\rm gen}$ is linear in the face fluxes, so this
boundary-face slice composes exactly with the transport operators'
interior-face slice. Everything is cell-local on the boundary cell; chunk
size $1+N_C$. The shared per-component flux kernel is
`_bc_boundary_face_fluxes` (film law + venting, positive into the domain).

`BCConvectionOperator` is the **single owner of the boundary cell's
$S_{\rm conv}$ energy row**. At a boundary cell the residual's gradient
blends the interior one-sided slope with the boundary-face slope at the
solved surface temperature (§7.4.2), so the row is the product

$$
S_{\rm conv,b} = -\sum_g c_{p,g}(T_b)\,
\underbrace{\tfrac12\bigl(J_{g,f_{\rm bc}} + J_{g,f_{\rm int}}\bigr)}_{\bar J_g}
\,\underbrace{\frac{d_{\rm nb}\,s_{\rm face}(T_s) + \Delta z_{1/2}\,s_{\rm int}}
{d_{\rm nb}+\Delta z_{1/2}}}_{(\partial T/\partial z)_b},
$$

whose factors are owned by three different operators — the BC face flux
(mass-BC kernel), the interior face flux (transport kernel), and the $T_s$
implicit chain (thermal-BC Newton). A split assembly would need
$\bar J \times \partial(\partial T/\partial z)/\partial x$ cross-products
spanning owners, so the whole row lives here instead: the operator re-runs
the Dual-aware surface Newton (outward mass flux hard-zeroed, mirroring the
residual; $s_{\rm face} = \pm(T_s-T_b)/\Delta z_{1/2}$ is algebraically
identical to the residual's $-q_{\rm total}[f]/k_b$ because the applied flux
is exactly $\mp k_b(T_s-T_b)/\Delta z_{1/2}$), evaluates both face-flux
kernels mode-appropriately (Fickian `gas_flux` or Darcy-Fick
`gas_flux_darcy_fick` for the interior face; `_bc_boundary_face_fluxes` for
the BC face), and ADs the whole product in one piece. The transport
operators zero their boundary-cell $S_{\rm conv}$ (keeping the boundary
cell's ξ-rows and the linear $S_{\rm gen}$ interior slice). Stencil: the
boundary cell and its interior neighbor's $(T,\xi)$ — inside the physics
tridiagonal; chunk size $2(1+N_C)$
(src/Jacobian/operators/bc_convection.jl).

---

## 15.9 Structured couplings

The non-local terms are represented by typed coupling objects, each implementing
`apply!(y,c,x)` (the action $y\mathrel{+}=Cx$), `apply_inverse!(y,c,x)`,
`densify!(M,c)` (for `DirectSolve`), and `size(c)`. All applies are O(N).

### 15.9.1 `PrefixSumCoupling` (exclusive prefix scan)

Represents $C=L\,\mathrm{diag}(v)$ with $L$ lower-triangular all-ones, so

$$
(C\,x)[k]=\sum_{i\le k-1}v[i]\,x[i]\qquad(\text{exclusive prefix sum}).
$$

The apply is a single forward pass accumulating $v[k]x[k]$ *after* writing
$y[k]$:

```julia
acc = 0
for k = 1:N
    y[k] += acc
    acc  += v[k]*x[k]
end
```

This is the cumulative ALE mesh velocity: $w_k$ depends on all *lower* cells'
strain rates. The matrix $L\,\mathrm{diag}(v)$ is **singular** — its first row is
zero and its last column is zero. The inverse (used by `StructuredSolve` when
composing with $\mathsf{A}$) is a first-difference recursion
$y[k]=(x[k+1]-x[k])/v[k]$ for $k<N$, with $y[N]=0$ by convention (the
right-nullspace direction is pinned by the rest of $\mathsf{A}$). It errors if
any $v[k]=0$.

### 15.9.2 `SuffixSumCoupling` (backward scan)

The mirror image, $C=U\,\mathrm{diag}(v)$ with $U$ strict upper-triangular:

$$
(C\,x)[k]=\sum_{i\ge k+1}v[i]\,x[i]\qquad(\text{exclusive suffix sum}),
$$

applied by a single backward pass. This is the Beer–Lambert non-local
absorption: cell $i$'s absorbed flux depends on all cells *above* it ($j>i$),
i.e. a strict suffix. The matrix is singular (last row zero, first column zero);
the inverse is the forward first-difference $y[k]=(x[k-1]-x[k])/v[k]$ for
$k\ge2$, $y[1]=0$.

### 15.9.3 `Rank1Coupling`

Represents $C=u\,v^{\mathsf T}$:

$$
(C\,x)[i]=u[i]\sum_j v[j]\,x[j],
$$

an inner product followed by a scaled scatter — O(N) apply. This is the
column-global $\theta_A$ / $\chi$ back-reaction in `WithChi`: the cross-section
$A=A_0\,\mathrm{law}(\bar\chi)$ depends on *all* $\chi_j$, so the affected rows
get a dense column dependence captured as a single outer product. Its inverse is
not formed elementwise; in a direct solve it is densified, and in
`StructuredSolve` it is handled via the Richardson loop (a Sherman–Morrison
formula could in principle be used, but the production path applies it through the
preconditioned Richardson iteration rather than an explicit rank-1 update).

### 15.9.4 Embedded variants

The bare couplings carry only their local math. The *embedded* variants wrap them
in a global embedding — row/column offsets plus per-cell scalings — so they can
be scattered into the full state vector:

- **`EmbeddedScaledPrefixSumCoupling`** — wraps a `PrefixSumCoupling` with a
  cell-centering average and a per-row scale, capturing
  $M[\text{row}_i,\text{col}_j]=\text{row\_scale}[i]\cdot\tfrac12(\mathbb{1}[j\le n_L-1]+\mathbb{1}[j\le n_R-1])\,v[j]$.
  Its apply computes
  $\text{out}[i]\mathrel{+}=\text{row\_scale}[i]\cdot\tfrac12(S[n_L]+S[n_R])$ with
  $S[k]=\sum_{j<k}v[j]x[j]$ — the cell-centered mesh velocity entering advection.
- **`EmbeddedScaledSuffixSumCoupling`** — the Beer–Lambert suffix at the right
  row/column offsets.
- **`EmbeddedRank1Coupling`** — the column-global $\chi$ back-reaction at its
  offsets; its `apply_inverse!` is intentionally undefined (handled by the
  solver).

The `compose_into!` machinery (called by consumers) is what lifts a provider's
bare coupling into the appropriate embedded form at the consumer's offsets.

> **Implementation note.** All coupling types, their `apply!`/`apply_inverse!`/
> `densify!`, and the `compose_into!` wrappers are in
> src/Jacobian/structured_couplings.jl. The rank-deficiency of the prefix/suffix
> couplings is documented and *load-bearing*: `apply_inverse!` sets the null
> direction to zero, and the surrounding $\mathsf{A}$ pins it — a subtlety easy to
> miss when reasoning about exactness.

---

## 15.10 Linear-solve strategies

A strategy inverts $\mathsf{A}+\sum_j C_j$ for the Newton step. All four operate
on the same assembled `StructuredJacobian`.

### 15.10.1 `DirectSolve` (anchor)

Materialize $M=\mathsf{A}+\sum_j\mathrm{densify}(C_j)$ and factor it directly
(KLU/UMFPACK):

```julia
M = Matrix(jac.A)
for c in jac.structured; densify!(M, c); end
out .= M \ b
```

Exact, robust, and the verification anchor; the default backend strategy. It
pays $O(N^2)$ to densify the non-local couplings, which is why the other
strategies exist.

### 15.10.2 `StructuredSolve` (preconditioned Richardson)

Solve $(\mathsf{A}+\sum_j C_j)\,y=b$ by Richardson iteration preconditioned with
$\mathsf{A}^{-1}$, applying each coupling through its fast `apply!`:

$$
F=\mathrm{lu}(\mathsf{A}),\quad y_0=F\backslash b,\qquad
\begin{aligned}
r_k &= b - \mathsf{A}y_k - \textstyle\sum_j C_j y_k,\\
\text{stop if } \|r_k\|&\le \text{tol}\,(\|b\|+1),\\
y_{k+1} &= y_k + F\backslash r_k.
\end{aligned}
$$

The default tolerance is $10^{-12}$ and the cap is 50 iterations; if it does not
converge it *throws* rather than returning a stale iterate. Convergence is fast —
**2–5 iterations to machine precision** — precisely because the implicit-RK $W$
shift makes $\mathsf{A}$ strongly diagonally dominant: the cumulative-coupling
entries are small relative to the $1/(\gamma_{\rm RK}\Delta t)$ + conduction +
reaction diagonal. Cost is $O(n_{\rm iter}\cdot n_{\rm structured})$ apply passes
plus the reused $\mathsf{A}$ factorization. When `jac.structured` is empty it
degenerates to the single direct $\mathsf{A}$-solve.

### 15.10.3 `ApproxSolve` (Schur-style, drop couplings)

Solve $\mathsf{A}^{-1}b$ only, dropping every structured coupling:

```julia
out .= jac.A \ b
```

Newton then sees a Jacobian that omits the dense ALE $C$ block — it takes a few
extra outer iterations, but each linear solve is far cheaper. The KB benchmark
reports roughly a $7.8\times$ Jacobian-side speedup on large ALE problems
(treat this as an order-of-magnitude figure, not a guaranteed constant; it
depends on $N$ and the coupling magnitudes). `ApproxSolve` uses the
`:ale_reduced` sparsity profile and is the **only** strategy carrying the
`is_approximate` trait, which switches the verifier to *masked* comparison
(in-support only) since the dropped couplings are *expected* omissions.

### 15.10.4 `MatrixFreeSolve` (GMRES over the operator)

Solve with GMRES over the `StructuredJacobian`'s action
$x\mapsto \mathsf{A}x+\sum_j C_j x$, never densifying the couplings. Experimental;
useful for very large $N$ where even forming $\mathsf{A}+C$ is expensive.
Production callers should prefer `DirectSolve` unless a problem-specific benchmark
shows the matrix-free path winning.

| Strategy | Inverts | Exact? | Sparsity profile | Verification |
|---|---|---|---|---|
| `DirectSolve` | $\mathsf{A}+\sum C_j$ densified | yes | `:full` | strict |
| `StructuredSolve` | $\mathsf{A}+\sum C_j$ via Richardson | yes (to tol) | `:full` | strict |
| `ApproxSolve` | $\mathsf{A}$ only | no | `:ale_reduced` | masked |
| `MatrixFreeSolve` | $\mathsf{A}+\sum C_j$ via GMRES | yes (to tol) | `:full` | strict |

> **Implementation note.** Strategy types: src/Jacobian/strategies.jl;
> `sparsity_profile(::Strategy)` maps `DirectSolve`/`StructuredSolve`/
> `MatrixFreeSolve` → `:full`, `ApproxSolve` → `:ale_reduced`. The four `ldiv!`
> methods live in src/Jacobian/structured.jl.

---

## 15.11 Sparsity construction and coloring

### 15.11.1 Two patterns: physics and geometry

The ALE geometry coupling motivates a **two-pattern** model:

- **$S_{\rm phys}$** — every coupling with $T$/$\xi$/$\chi$ columns:
  block-tridiagonal physics ($T\!\leftrightarrow\!T$, $T\!\leftrightarrow\!\xi$,
  $\xi\!\leftrightarrow\!\xi$ at each cell and its neighbors), the Beer–Lambert
  suffix (row $i$ depends on $\xi$-columns $j\ge i$, plus — on the Kirchhoff
  default path — the surface temperature column $T_{\rm col}(n)$ for every row
  $i<n$, the $\varepsilon_{\rm eff}(T_n,\xi_{\cdot,n})$ transmissivity chain of
  §15.8.7), the $\chi$-row block
  (cell-local plus the dense $\chi$-column), and the ALE advection limiter +
  cumulative widening.
- **$S_{\rm geom}$** — every coupling with $z$ columns (ALE only; empty for
  Eulerian): physics$\to z$ (cells read nodes; reach $\le i+3$ without $\chi$,
  full $1{:}n_{\rm nodes}$ with `WithChi` via the $\theta_A$ path), $z\to$physics
  (node $k$ depends on cells $i<k$, lower-triangular), and $z\to z$
  (lower-triangular including the diagonal, the cumulative nature of mesh
  velocity).

Five emitter functions compose these patterns, invoked in order by
`build_sparsity(layout, mode; profile)`:

1. `_emit_physics_pattern!` — block-tridiagonal;
2. `_emit_geometry_pattern!` — $z$-columns (ALE);
3. `_emit_chi_pattern!` — $\chi$-row and dense $\chi$-column (WithChi);
4. `_emit_radiation_pattern!` — Beer–Lambert suffix (BeerLambert);
5. `_emit_ale_advection_pattern!` — ALE limiter + cumulative (ALE).

The integrator's `jac_prototype` is the union $S_{\rm phys}\cup S_{\rm geom}$.
Patterns are memoized by the fingerprint
$(\mathrm{typeof}(\text{mode}),N_C,n_{\rm cells},n_{\rm nodes},\text{profile})$ in
`_SPARSITY_CACHE`.

### 15.11.2 Profiles

- `:full` (default) — all entries assembled; structured couplings emitted; the
  geometry AD pass runs.
- `:ale_reduced` (`ApproxSolve` only) — drops every structured coupling (the
  cumulative ALE widening, the z-row prefix-sum blocks, the $\theta_A$ rank-1
  channels, the Beer–Lambert suffix, and the $\chi$-column rank-1) and the
  geometry block. The masked verifier expects exactly this set of omissions.

### 15.11.3 Greedy coloring

For finite-difference / sparse-AD Jacobian construction (the integrator's
`AutoFiniteDiff` path and any colored AD), `matrix_colors` runs a greedy graph
coloring over a CSR representation of the pattern: columns sharing a color have
no common nonzero rows and can therefore be perturbed *simultaneously*, so the
number of RHS evaluations needed equals the number of colors rather than $N$. The
colorings are stored on the `JacobianSpec` (`colorvec_phys`, `colorvec_geom`).

### 15.11.4 The geometry operator (z-column AD)

Because $z$ couples through many indirect paths (mesh velocity, dilation,
$\chi$-dependent area), the $z$-columns are filled by a *separate* full-residual
`ForwardDiff` pass with respect to the $z$ subvector (`GeometryOperator`,
`GeometryADCache`), run *after* the operator walk and skipped under
`:ale_reduced`. This decouples the $z$-column computation from the per-operator
framework, at the cost that its entries must be consistent with what the
operators scattered into the $T/\xi/\chi$ rows independently — consistency that
verification (§15.13) checks.

The geometry pass covers **columns only**: it says nothing about the z-*rows*'
dependence on $(T,\xi,\chi)$. Those entries are structured couplings emitted by
`ALEMeshVelocityOperator` (§15.8.4) — the `EmbeddedNodePrefixSumCoupling`
prefix-sum blocks and, under `WithChi`, the $\theta_A$ rank-1 channels.

> **Implementation note.** `build_sparsity`, the five emitters, `matrix_colors`:
> src/Jacobian/sparsity.jl. `GeometryOperator`: src/Jacobian/operators/geometry.jl.
> `FastSparsityMap`: src/Jacobian/fast_sparsity_map.jl.

---

## 15.12 ILU-preconditioned Krylov solvers and the MKL ipiv caveat

For very large ALE problems where the direct factorization of $W$ dominates, two
iterative `LinearSolve.jl` factorizations are provided. They operate on the
*concrete* sparse $W=I/(\gamma_{\rm RK}\Delta t)-J$ matrix the integrator forms,
not on the `StructuredJacobian` (these are integrator-level linear solvers,
orthogonal to the §15.10 strategies):

- **`ILUGMRESFactorization(τ=0.1, gmres_rtol=1e-8, restart=30)`** — build an
  ILU($\tau$) preconditioner $P=\mathrm{ilu}(W,\tau)$ and solve $Wx=b$ by
  preconditioned GMRES. If ILU fails on a near-singular transient, it falls back
  to diagonal scaling $P=\mathrm{Diagonal}(1/\max(|{\rm diag}(W)|,10^{-12}))$.
- **`ILUBiCGSTABFactorization(τ=0.1, rtol=1e-8)`** — same preconditioner, BiCGSTAB
  iteration (fixed memory, no restart).

The drop tolerance default $\tau=0.1$ trades fill for preconditioner quality;
inexact solves can increase Newton rejections on Arrhenius transients, so the
iterative path is recommended only when the direct factorization is the
bottleneck.

**The MKL ipiv resize caveat.** When `handle_depletion=true`, the state vector is
resized mid-solve via `resize!(integrator, …)`. `LinearSolve.jl`'s default solver
can switch to `MKLLUFactorization` for large systems ($N\gtrsim500$), whose pivot
array `ipiv` is not always resized correctly, producing
`DimensionMismatch: ipiv has length X, but needs to be Y`. Julia's native
`LUFactorization` handles `resize!` correctly, and it is what
the depletion path uses: with `handle_depletion=true` there is no sparse
`jac_prototype` (the colored-FD path), so `_resolve_integrator` defaults the
linear solver to the dense, resize-safe `LUFactorization()` and rebuilds an
explicitly requested sparse factorization (Sparspak/KLU) into the user's
algorithm as `LUFactorization`, with a warning. The MKL pitfall therefore
survives only when a fully constructed integrator is passed with its
`linsolve` left unset — LinearSolve's auto-selection can then pick the MKL
LU — and `solve` warns when it detects that combination on a large depleting
state.

> **Implementation note.** `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`,
> and the Krylov integration: src/Solver/ilu_gmres.jl. The MKL ipiv warning:
> src/Solver/solve.jl (the `MKL_RESIZE_THRESHOLD = 500` block).

---

## 15.13 Verification

Every backend is anchored to `DenseAD_ForwardDiff` — `ForwardDiff` over the whole
residual, which is correct by definition because it differentiates the same
`residual!` the integrator integrates. `verify_jacobian_at(spec, def, ws, u, t;
reference=DenseAD_ForwardDiff(), rtol, atol)` compares the active backend's
Jacobian against this reference at a state point and returns
`(passed, max_abs_dev, max_rel_dev, worst_index, …, J_active, J_ref)`.

Two comparison modes, selected by the `is_approximate` trait:

- **Strict / support-complete** (exact backends, `is_approximate==false`):
  compares the *full* Jacobian under a **per-entry** tolerance with a
  row-magnitude floor,

  $$
  |\Delta J_{ij}| \le \text{atol} + \text{rtol}\cdot
  \max\big(|J_{\rm ref}[i,j]|,\ \max_j |J_{\rm ref}[i,:]|\big).
  $$

  Pyrolysis Jacobians span ~10 orders of magnitude across rows (Arrhenius
  ξ-rows ~$10^6$ vs ALE z-rows ~$10^{-2}$), so a global max-$|J_{\rm ref}|$
  normalization would let a wholly missing small-magnitude block hide behind
  the matrix maximum, with pass/fail depending on problem size rather than
  correctness. Under the per-entry rule, a missing entry contributes to its
  own row's reference scale and cannot hide, while honest roundoff at true
  zeros of a stiff row stays within $\text{rtol}\cdot\text{row max}$. Any missing
  entry the reference has nonzero (beyond its per-entry tolerance) is a
  failure.
- **Masked** (`is_approximate==true`, i.e. `Structured{ApproxSolve}`): compares
  only inside the declared support, under the coarse global-scale rule
  ($\max|\Delta J| \le \text{atol}$ or
  $\max|\Delta J|/\max|J_{\rm ref}| \le \text{rtol}$); entries the strategy
  *intentionally* omits (the dropped ALE/radiation/$\chi$ couplings of
  `:ale_reduced`) are informational, not failures. The in-support values of an
  approximate backend legitimately deviate by the dropped coupling magnitude,
  so masked mode is a sanity check, not an exactness gate.

Verification can run as a setup-time check or as a `VerifyJacobianCallback`
firing at chosen `sample_times` during a live solve. The callback *warns* on
deviation rather than aborting — a single transient outlier should not kill an
otherwise-good production run; the user reads `max_abs_dev`/`max_rel_dev` to
judge severity.

A reported incompatibility worth restating: **surface depletion and a pluggable
`JacobianSpec` cannot coexist** (the cache is sized at setup and cannot survive a
mid-solve `resize!`); with `handle_depletion=true` the spec is forced to
`nothing` and the integrator uses its colored finite-difference Jacobian, which
*is* verified implicitly by the integrator's own error control but is not
checked against the dense reference.

> **Implementation note.** `verify_jacobian_at`, strict/masked compare:
> src/Jacobian/verify.jl; `is_approximate`: src/Jacobian/interface.jl;
> the verify callback: src/Solver/solve.jl.

---

## 15.14 Comparison to ThermaKin2Ds, Gpyro, and FDS

We map each reference code's notation to ours on first use.

### 15.14.1 Gpyro

Gpyro (Lautenberger & Fernandez-Pello) solves the condensed phase by a
**fully-implicit backward-Euler** time discretization with **Patankar-style
linearization** of the source terms and a **tridiagonal (TDMA / Thomas)** linear
solve, sweeping the coupled energy/species/gas/pressure equations to convergence
each step. Mapping: Gpyro's volume fractions $X_\alpha$ are our $v_j$; its
pre-exponential $Z$ is our $A_i$; its permeability $K$ is our $\kappa$; its
volumetric reaction rate $\dot\omega'''$ is our $r$; its mixture average
$\bar\rho$ corresponds to our $\rho=\sum_j\xi_j$.

Contrasts with Pyrolysis.jl:

- **Time scheme.** Gpyro is first-order backward-Euler with an outer
  source-linearization sweep; Pyrolysis.jl uses an adaptive higher-order SDIRK
  (KenCarp4) with a true Newton solve per stage. Both are L-stable and implicit;
  Pyrolysis.jl gains temporal order and step-size adaptivity, at the cost of a
  more elaborate $W=I/(\gamma_{\rm RK}\Delta t)-J$ inversion than Gpyro's single
  TDMA pass.
- **Jacobian.** Gpyro never forms a global Jacobian; its Patankar linearization
  produces a tridiagonal coefficient matrix per equation, solved by TDMA. The
  coupling between equations is handled by the outer iteration. Pyrolysis.jl
  forms the *full coupled* Jacobian as a structured operator (block-tridiagonal
  $\mathsf{A}$ plus O(N) couplings), which captures cross-equation coupling
  directly inside the Newton step. The block-tridiagonal core of our
  $\mathsf{A}$ is the multi-field analog of Gpyro's per-equation tridiagonal
  matrix; our prefix/suffix couplings have no Gpyro analog because Gpyro's
  fixed-grid control volumes do not carry a cumulative-mesh-velocity term.
- **Mesh.** Gpyro uses fixed (or prescribed variable-$\Delta z$) control volumes
  with a solid-fraction / $\chi_k$ mechanism for swell/shrink; Pyrolysis.jl moves
  the mesh by ALE and represents the cumulative volume-change kinematics as a
  prefix-sum Jacobian coupling (§15.9).

### 15.14.2 ThermaKin2Ds

ThermaKin2Ds (Stoliarov et al.) solves the 1D (and 2D) condensed phase with a
**Crank–Nicolson** scheme in the primary (through-thickness, here $x$/our $z$)
direction and explicit treatment in the transverse directions, on a **dynamically
adjusted mesh**. Mapping: ThermaKin's stoichiometric coefficients $\theta_i^j$
are our $\nu_{i,j}$ (we reserve $\theta$ for strain rate, overload resolution
A2); its component concentrations are our $\xi_j$.

Contrasts:

- **Time scheme.** Crank–Nicolson is second-order and A-stable (not L-stable),
  trapezoidal in time; KenCarp4 is fourth-order, L-stable, with embedded error
  control. L-stability matters for the stiff Arrhenius transients at ignition,
  where Crank–Nicolson can ring; the L-stable SDIRK damps the fast modes
  cleanly. ThermaKin's fixed-$\theta$ trapezoidal rule contrasts with our
  adaptive SDIRK $\gamma_{\rm RK}$.
- **Energy convention.** ThermaKin (Eq. 17) carries gas enthalpy by the
  convective term $c_g(\mathbf{J}_g\cdot\partial T/\partial x)$ — exactly our
  $S_{\rm conv}$ — with gas storage excluded from the storage term, which is the
  Formulation-A convention Pyrolysis.jl adopts as canonical. The sign convention
  for heat of reaction differs in *publication* ($h>0$ exothermic in
  ThermaKin/Gpyro vs. $h>0$ endothermic internally here), but the volumetric
  source $Q_{\rm rxn}=-\sum h_r r_r$ is consistent once the convention is mapped
  (overload H1).
- **Mesh.** ThermaKin's dynamic mesh adjustment is the closest analog to our ALE
  + remapping; both track a receding surface. Pyrolysis.jl's distinction is the
  *Jacobian-aware* treatment: the moving-mesh kinematics enter the Newton matrix
  through the geometry-AD $z$-columns and the prefix-sum couplings, so the
  implicit step "knows" the mesh is moving.

### 15.14.3 FDS

The FDS (McDermott et al.) solid-phase pyrolysis model integrates a 1D Fourier
conduction equation with Arrhenius source terms; the surface boundary couples to
the gas-phase CFD. FDS uses a node-based 1D solid grid and an implicit/explicit
split appropriate to its multiphysics coupling. Pyrolysis.jl is a *standalone*
condensed-phase solver: it does not solve a gas-phase flow field, but it provides
a far richer in-solid model (multi-component effective properties, four radiation
models, Darcy gas transport, ALE with variable cross-section). FDS shares the
gas-storage-excluded energy convention and the Arrhenius kinetic form; its
two-flux (Schuster–Schwarzschild) radiation maps onto our Beer–Lambert / P1
choices (§8). On the *temporal/Jacobian* axis specifically, FDS does not expose a
structured-operator Jacobian of the kind described here — the relevant comparison
is methodological (implicit conduction + explicit/implicit coupling) rather than
a direct mapping of Jacobian internals.

---

## 15.15 Worked summary of the discrete update

Collecting the pieces, the production residual writes, for each cell $i$ and
component $j$:

$$
\frac{dT_i}{dt}=\frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+S_{{\rm gen},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}\Big] + \underbrace{w_{{\rm cell},i}(\partial T/\partial z)_i}_{\text{ALE only}},
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i}+\underbrace{w_{{\rm cell},i}(\partial\xi_j/\partial z)_i}_{\text{ALE, gas}}-\underbrace{\xi_{j,i}\theta_i}_{\text{ALE, non-gas}},
$$
$$
\frac{dz_k}{dt}=w_k=\sum_{i<k}(\theta_i-\theta_A)\Delta z_i\quad(\text{ALE}),\qquad
\frac{d\chi_i}{dt}=\frac{S_{{\rm dry},i}V_i}{m_{{\rm dry},i}}\quad(\text{WithChi}).
$$

The Jacobian of this RHS, $J=\partial\,\mathrm{RHS}/\partial\boldsymbol{u}$, enters
the implicit-RK $W=I/(\gamma_{\rm RK}\Delta t)-J$ and is assembled as the
structured operator $J=\mathsf{A}+\sum_j C_j$, where $\mathsf{A}$ is the
block-tridiagonal local block (reaction, conduction, transport, BC, dilation,
$\chi$-local, $z$-columns) and the $C_j$ are the prefix-sum (mesh velocity),
suffix-sum (Beer–Lambert), and rank-1 ($\theta_A$/$\chi$) couplings. The linear
solve $W\delta\boldsymbol{u}=-r$ is performed by the chosen strategy, and the
whole construction is verified against the dense AD reference.

---

## 15.16 Limitations and open issues

- **Geometry operator under `:ale_reduced`.** The $z$-columns are dropped
  entirely under `ApproxSolve`. Whether *any* residual $z$-coupling survives the
  profile (i.e. whether the profile is guaranteed to leave exactly zero $z$-rows)
  is asserted by the masked verifier rather than proven; the design intent is
  zero, and verification is the safeguard.
- **Sherman–Morrison for `Rank1Coupling`.** `StructuredSolve` handles the rank-1
  $\theta_A$/$\chi$ update through the preconditioned Richardson loop rather than
  an explicit Sherman–Morrison correction. An explicit rank-1 inverse could in
  principle accelerate convergence when the rank-1 magnitude is large, but the
  current path relies on $\mathsf{A}$'s diagonal dominance to converge.
- **Chunk-size tuning.** Per-operator `ForwardDiff` chunk sizes are hardcoded
  (e.g. $3(1+N_C)$ for conduction). They are the informationally-necessary
  stencil widths, not adaptively tuned for cache behavior; no auto-tuning is
  present.
- **`ImplicitTerm` path.** The IFT-resolved coupling type exists but the default
  boundary-Jacobian path is AD-through-Newton; explicit `ImplicitTerm` overrides
  are not exercised in the production walk.
- **WithChi × Darcy.** `WithChi` is tested mainly with Fickian transport; the
  interaction of the $\chi$ rank-1 coupling with the full Darcy operator set is
  not exhaustively exercised.
- **Depletion ↔ pluggable Jacobian / sensitivity.** Mid-solve `resize!` is
  incompatible with the sized `JacobianSpec` cache (forcing the colored-FD
  fallback) and, transitively, with adjoint sensitivity through depletion. Making
  the merge callback resize-aware for the structured backend is an open item.
- **ILU drop tolerance.** The default $\tau=0.1$ for the ILU-Krylov solvers has
  not been swept for optimality across pyrolysis regimes; inexact solves can
  raise Newton rejections on Arrhenius transients.
- **Warm-start artifacts.** Reusing a prior final state as fresh initial data can
  perturb finite-difference sensitivities and continuation studies; its effect is
  unquantified (§17.7).
- **Diagnostics cadence.** Energy/mass ledger snapshots are stride-based (every
  `diagnostics_stride` steps), not interpolated to output times, so the audited
  ledger is sampled rather than continuous (§16.8).
