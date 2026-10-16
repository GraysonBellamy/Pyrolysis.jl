# 15. Temporal Integration, DAE/Residual Formulation, and the Structured Jacobian

This chapter describes how the spatially-discretized governing equations of §13
(Finite-Volume Spatial Discretization) are assembled into a single flat
ODE/DAE right-hand side, integrated in time by stiff implicit Runge–Kutta
methods, and linearized through a structured-operator Jacobian. We derive the
block-major state-vector layout and the two dispatched residual paths; the
unified `residual!` entry point and the immutable-definition / mutable-workspace
split that makes it pure and thread-safe; the full ALE residual assembly
sequence; the link to OrdinaryDiffEq's implicit-RK $W$-matrix; the
three-layer structured-Jacobian architecture (operators → couplings →
linear-solve strategies); the exact per-operator Jacobian formulas with their
sign conventions; the O(N) prefix/suffix/rank-1 couplings that represent the
non-local terms; the linear-solve strategies and iterative solvers; the
two-pattern sparsity construction; and the verification
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
> the Jacobian sparsity module's index helpers `_xi_idx` and `_z_idx`,
> src/Jacobian/sparsity.jl, and in solution extraction,
> src/Solver/solve.jl. The block order $[T\mid\xi\mid z\mid\chi]$ is binding: any
> code that reads or writes raw state must respect it.
>
> **IC contract.** `create_state_vector` packs the mesh **as-is** — it does not
> consult the problem's initial-condition closures (which live on the
> `PyrolysisProblem`, not the mesh) and does not populate
> `mesh.m_dry_solid_initial`. The canonical setup path is the exported
> `initialize_state(problem; …) -> (def, ws, u0)`
> (src/Problem/adapter.jl), which mirrors `solve()` exactly:
> `apply_initial_conditions!(problem)` → `to_problem_def` → `build_workspace` →
> `create_state_vector`. Hand-rolled compositions that skip the IC step
> integrate from the mesh constructor's placeholder state — a
> resolution-independent error that still returns `retcode == :Success` —
> and, under `WithChi`, leave the χ denominator at zero (§12.8.3).

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
specializes fully — no runtime branching on the mode. The residual entry
dispatches on the mesh mode only, giving **two residual paths**, Eulerian and
ALE. Inside each, the transport law follows the property cache (a Darcy sub-cache
exists only in `DarcyFick` mode, and `compute_rhs!` selects the Darcy–Fick flux
from it), the $\chi$ block is included iff `WithChi`, and the radiation model
specializes the source.

> **Implementation note.** `SimulationMode`, the `MeshMode`/`RadiationMode`/
> `TransportMode`/`ChiMode` abstract hierarchies, and the accessors are in
> src/Problem/traits.jl. The mode is *derived*, not user-set directly:
> `to_problem_def(problem; use_ale, radiation_model, radiation_intensity,
> radiation_kirchhoff)` (src/Problem/adapter.jl) assembles it from the `use_ale`
> flag, the radiation
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
  workers. It carries the mesh connectivity (`MeshTopology`, a placeholder the
  1D residual does not read) and the *unchanging* specification: the material,
  the BC set, the `SimulationMode`, and a NamedTuple of runtime options
  (`radiation_model`, `radiation_intensity`, `radiation_kirchhoff`).
- `ws::Workspace` is **mutable** and owned by exactly one worker. It holds the
  typed state mirror (`StateCache`), property cache (`PropertyCache`), RHS work
  arrays (`RHSCache`), and (in ALE mode) the ALE work arrays. It also caches the
  dry-gas component indices and the boundary face IDs.
- `t` is a plain `Real`; the kernel accepts whatever scalar type $\boldsymbol{u}$
  and $\boldsymbol{du}$ carry — `Float64` for ordinary evaluation, or
  `ForwardDiff.Dual` for AD Jacobian assembly.

This split is what makes parameter sweeps (§ User Guide UG-12) and the
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
step (Newton iterations, Jacobian evaluations) all see a
consistent geometry, and that AD sensitivities propagate cleanly through
geometry because the geometry is a function of `cache.z` and `cache.A_section`,
both `Dual`-typed in the AD pass.

> **Implementation note.** `ProblemDef` and the runtime-option NamedTuple:
> src/Problem/problem_def.jl. `Workspace`: src/Problem/workspace.jl;
> `build_workspace` and `build_ad_workspace`: src/Problem/adapter.jl. `StateCache`, `populate!`,
> `update_state_geometry!`, `update_cross_section_area!`:
> src/Problem/state_cache.jl. The user-facing `PyrolysisProblem` and the
> `to_problem_def` adapter: src/Problem/pyrolysis_problem.jl and
> src/Problem/adapter.jl. The canonical `(def, ws, u0)` front door is
> `initialize_state(problem; …)` (same file) — it applies the problem's
> initial conditions before packing and is the path `solve()` itself takes;
> `to_problem_def`/`build_workspace`/`create_state_vector` are public as
> the advanced layer (AD workspaces, mid-simulation repacking, warm starts)
> with the IC step on the caller (§15.1.1).

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
   face-averaged $k_{\rm face}$, $\lambda_{\rm face}$).
3. `_compute_physics_rhs!` — calls `Discretization.compute_rhs!`, which picks
   the transport law from the property cache's Darcy sub-cache (`nothing` ⇔
   Fickian), refreshes that sub-cache ($P,\varphi,\kappa,\mu$ and the Darcy
   velocity, `update_darcy_properties!`) in its first step under `DarcyFick`, and
   populates `ws.dT_work[i]`, `ws.dξ_work[j,i]`, and
   `ws.rhs_cache.dξ_rxn[j,i]`.
4. Pack $\boldsymbol{du}$: $du[i]=dT_i/dt$, $du[n+(j-1)n+i]=d\xi_{j,i}/dt$.

The discrete updates assembled in step 3 are exactly those derived in §7 and
§13 (repeated here for completeness):

$$
\frac{dT_i}{dt} = \frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i + S_{{\rm conv},i} + Q_{{\rm rad},i} + Q_{{\rm rxn},i}\Big],
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
> src/Residual/residual.jl. The 8-step `compute_rhs!` assembly is in
> src/Discretization/finite_volume.jl (see §7.6 for the strict step ordering).

### 15.3.2 Order dependencies inside the physics RHS

The 8-step `compute_rhs!` sequence (transport properties and conductive fluxes →
reactions → gas transport → thermal BCs → radiation → $S_{\rm conv}$ →
divergences → final RHS) is *ordered* for correctness, not convenience:

- **Gas fluxes before thermal BCs.** The mass BCs (step 3) must finish before
  the boundary-temperature Newton solve (step 4), because the blowing model of
  each convective BC reads the outward gas mass flux $\dot{m}_g$ at the face.
- **Thermal BCs before $S_{\rm conv}$.** The boundary cells' temperature
  gradient in $S_{\rm conv}$ (step 6) blends in the boundary-face slope at the
  solved $T_s$, read from the flux that step 4 wrote.
- **Reaction rates before mesh velocity (ALE).** The ALE mesh velocity is built
  from the per-cell strain rate $\theta_i$, which is a function of the species
  reaction sources `dξ_rxn`; hence the physics RHS (which fills `dξ_rxn`) runs
  before the mesh-velocity integration.

---

## 15.4 The full ALE residual assembly

When `mesh_mode(def.mode) isa ALE`, mesh motion is part of the coupled system,
and the residual additionally evolves the node positions $z_k$ and (under
`WithChi`) the pyrolysis-progress $\chi_i$. The assembly is the following
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

5. **Physics RHS.** Same `compute_rhs!` as Eulerian; importantly it also fills
   `ws.rhs_cache.dξ_rxn`, which the ALE steps reuse.

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
   this is exactly the structure of the prefix-sum couplings in the Jacobian,
   §15.9). The radial channel under `WithChi` is, by the chain rule,
   $$
   \theta_A = \frac{1}{A}\frac{dA}{dt} = \frac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\frac{d\bar\chi}{dt},\qquad
   \frac{d\bar\chi}{dt} = \frac{1}{M_{\rm dry}}\sum_i S_{{\rm dry},i}\,V_i,
   $$
   with $\mathrm{law}'(\bar\chi)$ obtained by `ForwardDiff.derivative(law, χ̄)`
   and $S_{{\rm dry},i}=\sum_{j\in\text{dry-gas}}(d\xi_{j,i}/dt)|_{\rm rxn}$. The
   guard `law(χ̄) > ε` prevents the radial channel from diverging as the
   cross-section collapses.

7. **Mesh-relative temperature advection** (`compute_ale_advection_term!`,
   override-only). Every component rides the mesh, so the thermal mass
   $\rho c_p^{\rm eff}$ moves with it and $T$ is advected only by the mesh
   velocity relative to the material, $w_{\rm rel} = w_{\rm applied} −
   w_{\rm physical}$ (§11.2.2): in production $w_{\rm rel} \equiv 0$ and nothing
   fires; under `Workspace.w_override` the TVD-limited kernel adds
   $w_{\rm rel}\,\partial T/\partial z$ to $dT$.

8. **Dilation correction** for every component
   (`_apply_ale_dilation_correction!`):
   $$
   \frac{d\xi_{j,i}}{dt}\bigg|_{\rm dilation} \mathrel{-}= \xi_{j,i}\,\theta_i\quad(\text{every } j),
   $$
   which keeps each component's cell mass on its reaction-and-flux budget as
   the cell volume changes: $d(\xi_jV)/dt = S_jV − [J_jA]_{\rm faces}$. The gas
   is included — it is carried with the pores it occupies, and its flux
   $J_g$ is relative to the skeleton (§11.2.3).

9. **Pack the $T,\xi,z$ blocks.** $du[i]=dT_i/dt$ (with the ALE correction),
    $du[n+(j-1)n+i]=d\xi_{j,i}/dt$ (with dilation), and
    $du[n(1+N_C)+k]=w_k$ — i.e. $dz_k/dt=w_k$ (Euler-forward node motion
    realized by the integrator).

10. **Pack the $\chi$ block** (`_pack_chi_derivative_dispatched!`, `WithChi`
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
\frac{dT_i}{dt}=\frac{-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}}{\rho c_p^{\rm eff}{}_{,i}}
\;\Bigl[\,+ w_{{\rm rel},i}\Big(\frac{\partial T}{\partial z}\Big)_i \text{ under \texttt{w\_override} only; } w_{\rm rel}\equiv 0 \text{ in production}\Bigr],
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i}
- \underbrace{\xi_{j,i}\,\theta_i}_{\text{dilation (every component)}}.
$$

> **Implementation note.** ALE path: `_residual!(…, ::ALE)`,
> `_compute_ale_mesh_velocity!`, `_compute_theta_A` / `_compute_theta_A_dispatched`,
> `_apply_ale_T_advection!`, `_apply_ale_dilation_correction!`, and
> `_pack_chi_derivative!` are all in src/Residual/residual.jl. The strain-rate
> law `volume_change_rate` is in src/Physics/volume_change.jl. The ALE flux
> kernels are in src/Discretization/ale_fluxes.jl.

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
\texttt{KenCarp4}\big(\text{autodiff}=\text{AutoFiniteDiff()}\big).
$$

KenCarp4 is a 4th-order ESDIRK pair with an embedded error estimator for
adaptive step-size control; its singly-diagonal structure means every implicit
stage shares the same Newton iteration matrix, which we exploit below.

**Default linear solver.** `_resolve_integrator` picks the sparse direct
factorization from the Jacobian's pattern: `UMFPACKFactorization` when the
Jacobian materializes the ALE geometry block (`DirectSolve` under ALE,
`materializes_geometry_block`, §15.11.2), and `KLUFactorization` for every
other sparse Jacobian. The rule depends on the mode and the strategy, never on
the mesh size. It rests on end-to-end KenCarp4 runs of 100-cell PMMA and
charring-wood decks (Eulerian and ALE, surface and Beer–Lambert absorption,
with and without a lateral shrinkage law) and of V-BO-01:

- On the local patterns (every `ApproxSolve` Jacobian, Eulerian `DirectSolve`
  with surface absorption) KLU is 2.6–4.4× faster than Sparspak and 1.1–1.4×
  faster than UMFPACK, allocates 8–83 MB per solve against 19–534 MB for
  Sparspak, and takes the same steps to within 1%.
- On the exact ALE Jacobian UMFPACK is 1.5–1.9× faster than Sparspak and
  2.6–2.7× faster than KLU, with identical step and Newton counts.
- On the Eulerian Beer–Lambert pattern of `DirectSolve` KLU is 1.7× faster
  than Sparspak. UMFPACK's threshold pivoting (relative pivot tolerance 0.1)
  produced a factor with growth near $10^{50}$ on a $W$ of condition number
  about 30; the resulting Newton steps failed, and the run took 50 times the
  steps. That pattern therefore uses KLU.
- All three factorizations are backward stable on these systems (normwise
  backward error below $3\times10^{-15}$). Sparspak does not pivot: on the
  Darcy decks its componentwise backward error reaches 1 in individual rows,
  where KLU's and UMFPACK's stay below $10^{-13}$. On the exact ALE Jacobian the
  relative residuals $\lVert Wx-b\rVert/\lVert b\rVert$ reach
  $6\times10^{-11}$ for UMFPACK, $2\times10^{-11}$ for KLU and
  $3\times10^{-9}$ for Sparspak, without effect on the Newton iteration.

Sparspak, KLU, UMFPACK and the ILU-Krylov solvers of §15.12 are interchangeable
through `linear_solver`.

> **Implementation note.** Integrator resolution: `_resolve_integrator`,
> src/Solver/solve.jl. A user may pass a constructor (`integrator=Rodas5P`) plus
> a `linear_solver=`, or a fully-built algorithm. The RHS closure re-samples
> the incident radiation at every evaluation (`ws.I_radiation = radiation_flux(t)`)
> whenever the active model consumes it (`SURFACE_ABSORPTION` / `BEER_LAMBERT`)
> — constant and time-varying fluxes alike (§8.3.4). The Jacobian closure
> refreshes `ws.I_radiation` at its own $t$ the same way, so a time-varying flux
> enters the Beer–Lambert Jacobian at the time the Jacobian is evaluated.

**Problem-independent integrator types.** OrdinaryDiffEq compiles its stepper,
nonlinear solver and factorization for the types of the functions it is handed.
The RHS and Jacobian closures capture the `ProblemDef` and workspace, whose
types carry the material, boundary-condition and mode types. `solve` therefore
hands them to the integrator behind fixed-signature `FunctionWrapper`s:
`ErasedRHS` for `(du, u, p, t)` with `Vector{Float64}` state and `Float64`
time, and `ErasedJacobian{M}` for the Jacobian matrix type `M`. The callbacks
do the same (§15.5.3). The integrator's type then depends only on the solver
configuration: the mesh mode, depletion handling, the installed optional
callbacks, the integrator and its factorization. A problem with new material or
boundary-condition types compiles only the Pyrolysis kernels (initialization,
residual, Jacobian, solution extraction). The precompile workloads
(`src/precompile.jl`) cache the integrator for the default configurations. An
algorithm that calls the RHS with other argument types, such as dual numbers
for a ForwardDiff time derivative, reaches the closure by dynamic dispatch. A
call of the Jacobian with a different matrix type is a `MethodError`, because
a converted copy would silently drop the entries.

> **Implementation note.** src/Solver/integrator_functions.jl (`ErasedRHS`,
> `ErasedJacobian`); the `ODEFunction` is built in `solve`, src/Solver/solve.jl.
> `solve` also converts the tolerances, step bounds and `saveat` to `Float64`,
> because the solve keywords enter the integrator's type.
> test/solver/test_integrator_type.jl checks that two decks differing in every
> problem type share one integrator type.

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
$1/(\gamma_{\rm RK}\Delta t)$ scale.) Each Newton step solves
$W\,\delta\boldsymbol{u}=-r$ for the stage residual $r$ with the sparse
factorization of §15.5.1; the strategy (`DirectSolve` / `ApproxSolve`, §15.10)
decides which entries of $J$ are assembled.

The Jacobian $J$ is supplied to the integrator through an `ODEFunction` carrying
both `jac = jac!` and a `jac_prototype` sparsity pattern (the union of the
physics and geometry patterns of §15.11). `solve` builds the spec for its own
`(def, ws)`:

```julia
jac_spec = JacobianSpec(backend, def, ws)
jac! = (J, u, p, t) -> compute_jacobian!(J, jac_spec, u, def, t, ws)
jac_prototype = Float64.(jac_spec.sparsity_phys .| jac_spec.sparsity_geom)
ode_fn = ODEFunction{true,FullSpecialize}(ErasedRHS(rhs!);
    jac = ErasedJacobian{typeof(jac_prototype)}(jac!),
    jac_prototype = jac_prototype)
```

(with the incident flux refreshed at the Jacobian's own $t$ inside `jac!`,
§15.5.1).

**Warm start.** A prior solution's final state can seed a new solve's initial
condition (`warm_start=prev_solution`), provided the component count, mesh size,
ALE mode, and $\chi$ mode all agree (enforced implicitly by a state-vector
length check) and the prior solve parked no cells. The integrator's `tspan[1]` still defines the start time, so
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
| 2 | thickness termination | continuous (ALE) | terminate when the live thickness $z_{n_{\rm live}+1}-z_1<$ `min_thickness` |
| 3 | `DepletionCallback` | discrete (opt-in) | park cells below their depletion threshold; the state keeps its size |
| 4 | Newton-failure cap | discrete (every accepted step; Newton integrators) | after a step whose Newton needed retries, hold the proposed step at or below 1.25× that step for 20 steps |
| 5 | progress bar | discrete (every accepted step) | progress bar |
| 6 | `DiagnosticsCallback` | discrete (every `diagnostics_stride`) | log a sensible-energy snapshot row (ledger fields NaN; §16.8.3) |
| 7 | Jacobian verification | discrete (at sample times) | compare active Jacobian vs. dense AD reference |
| 8 | user callback | — | custom |

Every callback of rows 1–7 keeps its problem-typed work behind a fixed-signature
`FunctionWrapper` on the state vector and the time: the mesh sync, the
live-thickness evaluation, the depletion check and park, the diagnostics row and
the Jacobian comparison. Only the actions on the integrator itself stay in the
callback: `terminate!`, `set_proposed_dt!`, `derivative_discontinuity!` and the
reads of its step and Newton-failure counters. The `CallbackSet`'s type
therefore depends on which callbacks are installed, not on the problem (§15.5.1).

The **`StepUpdateCallback`** is the linchpin of the read-only-mesh design: after
each accepted step it sets `mesh.cross_section_area = cross_section_area_from_state(…)`
and calls `update_mesh_from_state!(mesh, u, layout)`, so all between-step
consumers (depletion detection, diagnostics, the surface BC solve in
post-processing) see a mesh consistent with the integrator state.

The **`DepletionCallback`** (§14.3) detects cells whose current thickness has
fallen below `depletion_threshold × Δz_init,i` (default $5\%$ of each cell's
own initial thickness; on a fixed mesh, the top live cell whose condensed
concentration has fallen below the same fraction of its initial value), merges
it into a neighbour with the sensible-enthalpy-conserving temperature merge of
§14.3.3,

$$
\sum_{j} m_j\!\int_{T_{\rm ref}}^{T_{\rm merged}}\! c_{p,j}\,dT'
= H_a(T_a) + H_b(T_b)
\ \text{(bracketed Newton; every component, gas included)},\qquad
\xi_{{\rm merged},j} = \frac{V_a\,\xi_{a,j}+V_b\,\xi_{b,j}}{V_{\rm keep}},
$$

and parks the freed slot above the live surface
(`derivative_discontinuity!(integrator, true)`). The state keeps its size, so
the Jacobian's sparsity pattern and caches stay valid and `saveat`
interpolation is unaffected (§14.3.5–§14.3.6).

The **Newton-failure cap** (`create_newton_failure_cap_callback`) is a step cap with a memory of Newton failures. After
every accepted step it reads the integrator's count of failed nonlinear solves.
When the step just accepted needed retries, the proposed next step is held at or
below $g\,\Delta t_{\rm acc}$ with $g = 1.25$ (`FAILURE_MEMORY_GROWTH`) for the
next $M = 20$ accepted steps (`FAILURE_MEMORY_STEPS`); a further failure
restarts the memory from the new step. The error controller alone regrows the
step to its accuracy limit after a failure, and on a strongly nonlinear front
that limit lies beyond the range in which the stage Newton converges, so the
integrator settles into a fail–halve–regrow rhythm. On the Landau ablation case
V-CPL-02 the Newton fails near $\Delta t = 5$ ms and converges at 2.5 ms while
the controller keeps proposing 5 ms: the cap takes the case from 1087 steps and
892 Newton failures to 1032 steps and 78 failures. A run without Newton
failures never triggers the cap and keeps its trajectory exactly. It is
installed for every integrator that solves its stages by Newton (an
algorithm with an `nlsolve` field).

> **Implementation note.** Callbacks: src/Solver/callbacks.jl
> (`create_step_update_callback`, `create_thickness_termination_callback`,
> `create_depletion_callback`, `park!`, `create_newton_failure_cap_callback`)
> and src/Solver/solve.jl (assembly, diagnostics, verify-Jacobian). The
> cross-keyword checks are in `_validate_solve_kwargs`, src/Solver/solve.jl.

---

## 15.6 Jacobian backends

The production backend is **`Structured{S}`**, parameterized by a
`LinearSolveStrategy` $S$ (`DirectSolve` or `ApproxSolve`). It assembles the
Jacobian *operator by operator*, exploiting the block-tridiagonal physics
structure and representing the non-local terms as O(N) couplings until they are
written into the sparse matrix. A second backend, **`DenseAD_ForwardDiff`**,
runs `ForwardDiff` over the *entire* residual vector and returns a dense
Jacobian. It is correct by construction (it differentiates the same `residual!`
the integrator uses) but $O(N^2)$ in storage and slow; it is the
**verification reference** (§15.13), not a solve backend.

`solve` takes the bare backend: `jacobian = Structured(strategy =
DirectSolve())` or `Structured(strategy = ApproxSolve())`. With `jacobian =
nothing` (the default) it picks `DirectSolve` on a fixed mesh — exact assembly
plus a sparse direct solve — and `ApproxSolve` under ALE, which drops the
structured couplings (§15.10.2). `solve` builds the spec for its own mesh copy,
mode, and options; the `jacobian` keyword is typed `Union{Structured,Nothing}`,
so passing a prebuilt spec or `DenseAD_ForwardDiff()` raises a `TypeError`.

`JacobianSpec` bundles the backend with everything the integrator needs:

```julia
struct JacobianSpec{B<:JacobianBackend, S, C}
    backend::B
    sparsity_phys::S      # SparseMatrixCSC{Bool,Int}, or nothing (dense reference)
    sparsity_geom::S      # geometry (z-column) pattern; nothing for Eulerian
    cache::C              # backend-private workspace (not shareable)
end
```

It is built by `JacobianSpec(backend, def, ws)` — the constructor `solve` calls,
and the one the verifier and other advanced callers use with a `(def, ws)` pair
from `initialize_state`. A single entry point
`compute_jacobian!(J, spec, u, def, t, ws)` dispatches to
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

Each physics block is an `Operator` subtype: a lightweight stateless tag with two
methods, its residual contribution `apply!` (a thin wrapper over the residual
kernel it differentiates) and its Jacobian contribution
`contribute_to_jacobian!(jac, op, u, t, ws, args...)`, where `ws` is the
operator's workspace from the backend cache and the remaining arguments are
operator-specific. There is no default fallback: each operator implements its
own contribution, as `ForwardDiff` over a stencil-local kernel (the
`GeometryOperator` differentiates the full residual with respect to the $z$
block), or as an analytical override (conduction, gas transport, Darcy transport)
tested against a ForwardDiff reference of the same kernel at `rtol = 1e-12`. Local entries are added into the
sparse block `jac.A`; non-local blocks are pushed onto `jac.structured` as
`StructuredCoupling`s. The ALE mesh-velocity operator emits its non-local
blocks through a separate `contribute_to_structured!`; the $\chi$ operator pushes
its rank-1 block from its own `contribute_to_jacobian!`.

The backend walks the operators in a fixed order:

```
1. (ALE)      ALEMeshVelocityOperator      apply! + strain-rate (and θ_A) derivatives,
                                           then contribute_to_structured! (z-row prefix sums, θ_A rank-1)
2. ReactionSourceOperator
3. ConductionOperator
4. GasTransportOperator | DarcyTransportOperator
5. RadiationOperator
6. BCThermalOperator
7. BCMassOperator
8. BCConvectionOperator                    (boundary-cell S_conv energy row)
9. (ALE)      dilation rows                (cell-local, from the cached ∂θ)
10.(ALE & :full) GeometryOperator          (full-residual AD over z-columns)
11.(WithChi)  ChiOperator                  (cell-local rows; χ-column rank-1)
```

The order is fixed so that the data each step reads exists: the state cache,
the ALE geometry, and the properties are refreshed first (mirroring the
residual), the mesh-velocity workspace and its strain-rate derivatives are filled
before the dilation rows read them, and the $\chi$ operator reads the
cross-section derivatives the mesh-velocity workspace holds. The radiation,
mesh-velocity and $\chi$ operators build their embedded couplings at the global
row/column offsets of the blocks they fill; `compose_into!` embeds the prefix-sum
factor vectors at the node rows.

> **Implementation note.** Operator interface: src/Jacobian/operators/Operator.jl;
> per-operator implementations under src/Jacobian/operators/; the walk is
> `_compute_jacobian!` in src/Jacobian/structured_backend.jl. Per-operator
> `ForwardDiff` chunk sizes: reaction `1+N_C`; conduction $3(1+N_C)$; gas
> transport $3(1+N_C)$; BC thermal $1+N_C$; BC mass $1+N_C$; BC convection
> $2(1+N_C)$; $\chi$ $1+N_C$. The ALE dilation derivatives come from the
> mesh-velocity operator's strain-rate sweep.

### 15.7.2 Layer 2 — `StructuredJacobian`

The assembled Jacobian is *not* a monolithic sparse matrix. It is

```julia
struct StructuredJacobian
    A::SparseMatrixCSC{Float64,Int}        # local block: block-tridiag physics + sparse augments
    structured::Vector{StructuredCoupling} # non-local couplings in factored form
    smap::FastSparsityMap                  # O(1) sparse-write helper
    n::Int                                 # full state dimension
end
```

The split is the architectural heart of the design: the sparse block
$\mathsf{A}$ carries everything *local* (block-tridiagonal mass / conduction /
reaction / transport plus the cell-local self-terms), and the `structured`
vector carries everything *non-local* (the cumulative ALE mesh velocity, the
Beer–Lambert optical-depth suffix, the column-global $\theta_A$ and $\chi$
back-reactions) in factored form. The assembly writes them into the sparse
Jacobian with `densify!` only under the exact `:full` profile.

`FastSparsityMap` is an O(1) sparse-write helper: a dense lookup table
`indices[row + (col-1)*n_rows] = nz_idx` (zero meaning "not in pattern"). It
trades quadratic memory ($\sim2.5\,\mathrm{MB}$ at 560 DOFs, $\sim90\,\mathrm{MB}$ at
3360 DOFs) for O(1) scattered additions via `fast_sparse_add!`.

### 15.7.3 Layer 3 — `LinearSolveStrategy`

A strategy selects which couplings the assembled Jacobian carries; it is
orthogonal to the assembling backend. `DirectSolve` writes every structured
coupling into the sparse matrix (the exact Jacobian), and `ApproxSolve` drops
them (the reduced pattern). The integrator then factors the resulting sparse $W$
with KLU or UMFPACK (§15.5.1). The two strategies are detailed in §15.10.

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
the tridiagonal block, assembled by a hand-coded analytical override
(`_contribute_to_jacobian_analytical!`) that is tested against a `ForwardDiff`
reference over the 3-cell stencil (`_contribute_to_jacobian_ad!`, chunk size
$3(1+N_C)$) at `rtol = 1e-12`; each of the three cells contributes
$(T,\xi_1,\dots,\xi_{N_C})$ inputs because $k_f$ and $\rho c_p^{\rm eff}$ both
depend on composition.

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
rule, $\mu_i$ from the material's viscosity curve, the face mobility by the residual's
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
boundary cell's ξ-rows but contribute no boundary-cell $S_{\rm conv}$. Each
operator's `contribute_to_jacobian!` is a hand-coded analytical override, gated
against a ForwardDiff reference of the same kernel (`_contribute_to_jacobian_ad!`)
at `rtol = 1e-12`.

> **Implementation note.** `GasTransportOperator`,
> src/Jacobian/operators/gas_transport.jl; `DarcyTransportOperator`,
> src/Jacobian/operators/darcy_transport.jl (`_dt_props` Dual pass,
> `_dt_mobility_face`, `_dt_J_face_partials`, analytical override
> `_contribute_to_jacobian_analytical!`).

### 15.8.4 ALE mesh velocity and the dilation correction

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
scaling $v[i]=\Delta z_i\,\partial\theta_i/\partial x_i$ is held in a
`PrefixSumCoupling{Float64}` (a plain carrier of the per-cell factors) per source
channel ($T$ and each $\xi_j$), which `compose_into!` embeds at the z-row offsets
as `EmbeddedNodePrefixSumCoupling` blocks. The z-row equations $dz_k/dt = w_k$
are themselves consumers of the cumulative coupling, with $\partial(dz_k/dt)/\partial x_j =
\mathbb{1}[j\le k-1]\,v[j]$. The geometry AD pass (§15.11.3) covers z
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

where $G=\mathrm{law}'(\bar\chi)/\mathrm{law}(\bar\chi)$. The rank-1s (shared
left factor $u[k]=-\text{cumZ}_{\rm node}[k]$) are embedded directly at the
z-rows as `EmbeddedRank1Coupling` blocks.

The **dilation correction** $d\xi_{j,i}\mathrel{-}=\xi_{j,i}\theta_i$ adds the
cell-local Jacobian (every component $j$)

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

### 15.8.5 No ALE advection operator

Every component rides the mesh and the thermal mass rides with it (§11.1.2,
§11.2.2), so neither the ξ rows nor the T row carry a mesh-motion coupling
through $w$ in production: the only mesh-motion term in the species rows is the
cell-local dilation of §15.8.4, and the mesh-relative temperature kernel runs
only under the `w_override` verification hook, which the structured Jacobian
does not support (override runs use finite-difference Jacobians). The cumulative
mesh velocity therefore couples only the z rows, through the prefix-sum blocks
the mesh-velocity operator embeds itself. The block-tridiagonal physics
pattern is identical under the `:full` and `:ale_reduced` profiles; the profiles
differ in the geometry block and the structured couplings (§15.11.2).

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
  bounding node positions, the entries are
  $\partial(d\chi_i/dt)/\partial z_i=-A\,S_{\rm dry}/m_{\rm init}$ and
  $\partial(d\chi_i/dt)/\partial z_{i+1}=+A\,S_{\rm dry}/m_{\rm init}$. The
  `ChiOperator` does not write them: the `GeometryOperator`'s AD pass over the
  $z$ subvector fills them under `:full` (§15.11.3), and `:ale_reduced` drops them.
- **$\chi$-columns (rank-1)**: the area $A=A_0\,\mathrm{law}(\bar\chi)$ depends on
  *every* $\chi_j$ through $\bar\chi$, giving a dense column coupling factored as
  an `EmbeddedRank1Coupling` with
  $u[i]=S_{{\rm dry},i}\,\Delta z_i/m_{{\rm init},i}$ ($\Delta z_i=V_i/A$) and
  $v[j]=dA/d\chi_j$, pushed only under `:full`.

> **Implementation note.** `ChiOperator`, `apply!`, `_chi_residual_kernel!`,
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
$j=i$ term is added by the operator through a cell-local `ForwardDiff` pass (with
$I_{\rm in}$ held fixed); the strict-suffix $j>i$
couplings are captured as an `EmbeddedScaledSuffixSumCoupling`. The split
($j=i$ cell-local, $j>i$ structured) avoids double-counting the diagonal.

On the Kirchhoff default path (`radiation_kirchhoff = true`, §8.3.2) the
injected intensity is $a_{\rm eff}(T_n,\xi_{\cdot,n})\,I_{\rm surface}$ (the
effective surface absorptivity of §5.8.2, equal to $\varepsilon_{\rm eff}$ unless a
component sets its own `absorptivity`), so every energy row $i<n$ additionally
couples to the **surface cell**: the suffix vector carries
$v_n \mathrel{-}= (\partial a_{\rm eff}/\partial x_n)/a_{\rm eff}$,
the surface cell's diagonal closure recomputes $a_{\rm eff}(T,\xi)$
live, and the sparsity pattern gains $T$-row($i$) × $T$-column($n$) slots for
$i<n$ (`_emit_radiation_pattern!`). With an explicit BC absorptivity the factor
is folded into $I_{\rm surface}$ upstream and these terms vanish.

Under ALE the operator fills **T and ξ columns only**. The dependence of the
deposit on the live geometry ($\Delta z_i = V_i/A$ from the z block) lands in the
z columns, which the `GeometryOperator`'s AD pass over the full residual writes
under the exact `:full` profile (§15.11.3); the `:ale_reduced` profile drops
every z column (src/Jacobian/operators/radiation.jl).

> **Implementation note.** `RadiationOperator`,
> src/Jacobian/operators/radiation.jl. Active only in `BeerLambert` mode;
> `SurfaceAbsorption`/`NoRadiation` leave a tridiagonal pattern.

### 15.8.8 Boundary conditions

`BCThermalOperator` overwrites the boundary-cell conductive flux with the
BC-specified surface flux (convective + radiative + fixed-$T$), obtained from the
surface-temperature Newton solve $f(T_s)=q_{\rm BC}(T_s)-k_{\rm eff}(T_s-T_{\rm cell})/\Delta z_{1/2}=0$
(§12.2). Its Jacobian is taken by AD through that Newton solve, producing the
boundary-cell row's dependence on $(T_{\rm cell},\xi_{\rm cell})$. The outward
gas mass flux fed to the Newton solve is rebuilt from the stencil state by the
shared boundary-face kernel `_bc_boundary_face_fluxes` (film law plus, in Darcy
mode, the venting flux) — the same flux the residual passes to the blowing models
— so the $(T_b,\xi_b)$ derivatives of a blowing model's coefficient $h\chi$ ($\chi$
the blowing factor of §12.3.3, not the pyrolysis progress) are exact
(§12.3.3). Chunk size $1+N_C$.

`BCMassOperator` owns the species rows of the boundary-face gas flux $F_g$
(film law plus, in Darcy mode, the pressure-BC advective venting),
$d\xi_{g,\rm bdy} \mathrel{+}= F_g A_f/V$, because the residual's
$\nabla\!\cdot\!J_{\rm gas}$ is formed *after* `apply_mass_bc_flux!` fills the
boundary faces, while the transport operators zero those faces. Everything is
cell-local on the boundary cell; its block is $N_C\times(1+N_C)$ and the chunk
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
the Dual-aware surface Newton (fed the outward flux from the shared
boundary-face kernel, so each convective member's blowing model is differentiated
exactly; $s_{\rm face} = \pm(T_s-T_b)/\Delta z_{1/2}$ is algebraically
identical to the residual's $-q_{\rm total}[f]/k_b$ because the applied flux
is exactly $\mp k_b(T_s-T_b)/\Delta z_{1/2}$), evaluates both face-flux
kernels mode-appropriately (Fickian `gas_flux` or Darcy-Fick
`gas_flux_darcy_fick` for the interior face; `_bc_boundary_face_fluxes` for
the BC face), and ADs the whole product in one piece. The transport
operators zero their boundary-cell $S_{\rm conv}$ (keeping the boundary
cell's ξ-rows). Stencil: the
boundary cell and its interior neighbor's $(T,\xi)$ — inside the physics
tridiagonal; chunk size $2(1+N_C)$
(src/Jacobian/operators/bc_convection.jl).

---

## 15.9 Structured couplings

The non-local terms are represented by typed coupling objects, subtypes of
`StructuredCoupling`, each carrying O(N) data — a per-cell factor vector and its
global row/column offsets — and implementing `densify!(M, c)`, the addition
$M \mathrel{+}= C$ that the `:full` profile uses to write the block into the
sparse Jacobian. All three shapes are *embedded*: they know where they sit in
the full state vector.

### 15.9.1 `EmbeddedNodePrefixSumCoupling` (exclusive prefix sum at the node rows)

The cumulative ALE mesh velocity: node $k$'s velocity depends on the strain rates
of all *lower* cells,

$$
M[\text{row}_k,\text{col}_j] \mathrel{+}= \mathbb{1}[j\le k-1]\,v[j],
\qquad v[j]=\Delta z_j\,\frac{\partial\theta_j}{\partial x_j},
$$

i.e. $C=L\,\mathrm{diag}(v)$ with $L$ strictly lower-triangular all-ones
(exclusive prefix sum). Node 1, the pinned bottom node, gets no entries. The
per-cell factors arrive in a `PrefixSumCoupling` — a plain carrier of the vector
$v$, not itself a `StructuredCoupling` — which `compose_into!` embeds at the node
rows.

### 15.9.2 `EmbeddedScaledSuffixSumCoupling` (strict suffix sum)

The mirror image, with a per-row scale:

$$
M[\text{row}_i,\text{col}_j] \mathrel{+}= \text{row\_scale}[i]\;\mathbb{1}[j\ge i+1]\;v[j],
$$

i.e. $C=\mathrm{diag}(\text{row\_scale})\,U\,\mathrm{diag}(v)$ with $U$ strictly
upper-triangular all-ones. This is the Beer–Lambert non-local absorption: cell
$i$'s absorbed flux depends on all cells *above* it ($j>i$), with
$\text{row\_scale}[i]=-q_{\rm rad}[i]/\rho c_{p,i}^{\rm eff}$ for the energy row
(§15.8.7). The diagonal $j=i$ is not part of the coupling; the radiation operator
adds it into $\mathsf{A}$ through a cell-local `ForwardDiff` pass.

### 15.9.3 `EmbeddedRank1Coupling`

A rank-1 block $C=u\,v^{\mathsf T}$ at given offsets:

$$
M[\text{row}_i,\text{col}_j] \mathrel{+}= u[i]\,v[j].
$$

This is the column-global $\theta_A$ / $\chi$ back-reaction in `WithChi`: the
cross-section $A=A_0\,\mathrm{law}(\bar\chi)$ depends on *all* $\chi_j$, so the
affected rows get a dense column dependence captured as a single outer product
(the $\theta_A$ channels at the z rows, §15.8.4, and the $\chi$-column block of
the $\chi$ rows, §15.8.6).

> **Implementation note.** `PrefixSumCoupling`, the three embedded coupling types,
> their `densify!`, and `compose_into!` are in src/Jacobian/structured_couplings.jl.
> The couplings alias their factor vectors (no copies) and are rebuilt at every
> Jacobian evaluation; `zero!` empties `jac.structured` at the start of each
> assembly.

---

## 15.10 Linear-solve strategies

A strategy decides which couplings the assembled Jacobian carries. Both operate on
the same assembled `StructuredJacobian`; the integrator then factors the sparse
$W$ with the factorization of §15.5.1.

### 15.10.1 `DirectSolve` (anchor)

Write every structured coupling into the sparse matrix:
$J=\mathsf{A}+\sum_j\mathrm{densify}(C_j)$ under the `:full` sparsity profile.
The linear solve then sees the exact Jacobian. It is the verification anchor and
the default on a fixed mesh. Under ALE the prefix-sum and geometry blocks couple
every node to every cell below it, which makes the factorized $W$ effectively
dense — the reason for the other strategy.

### 15.10.2 `ApproxSolve` (Schur-style, drop couplings)

Keep the local block only, $J=\mathsf{A}$, under the `:ale_reduced` sparsity
profile: every structured coupling is omitted at its producer and the geometry
block is not assembled. Newton then sees a Jacobian that omits the dense ALE
block, while each linear solve stays local and cheap. `ApproxSolve` is the
**only** strategy carrying the `is_approximate` trait, which switches the
verifier to *masked* comparison (in-support only) since the dropped couplings are
*expected* omissions.

It is the default under ALE. With $W = I-\gamma\Delta t\,J$ the modified
Newton iteration that factors the reduced matrix has the error operator
$E = I - W_{\rm red}^{-1}W_{\rm full}$, whose spectral radius is its
asymptotic contraction per iteration. It was estimated at every accepted step
of KenCarp4 runs, at $\gamma$ times the larger of the step taken and the step
proposed, by power iteration in the integrator's weighted norm and checked
against dense eigenvalues, on plain-ALE decks without parking at 100 cells:
PMMA at 25, 50 and 100 kW/m² (solid $\lambda = 10^{-6}$, and 0 at 50 kW/m²)
run to 5 % of its thickness; a charring wood (Darcy, Beer–Lambert) whose char
shrinks the cells to 0.83 of their volume under a lateral shrinkage law; an
intumescent coating that swells 3.6× in volume under a widening lateral law;
the Landau ablation deck of V-CPL-02; and 80-s cone decks of PMMA and charring
wood. It never exceeded 0.19, and it peaked before any step collapse, so the
iteration contracts by at least a factor of five per iteration. $E$ is
strongly non-normal: because every $z$ column is dropped, its one-iteration
weighted norm reaches $10^2$ on the cone decks and far more once cells
collapse, with no visible cost in steps or Newton failures.

The two strategies follow the same trajectories:

- The decks that run to their end take the same steps to within two (54/54,
  979/979, 82/83, 80/80, 337/336, 286/288 reduced/exact) and agree to
  $3\times10^{-4}$ K in surface temperature and $5\times10^{-6}$ of the peak
  in gas generation. PMMA at 25 kW/m² reaches 5 % of its thickness at the same
  time, 1617.5 s, in 21,036 and 22,607 steps.
- Without parking, surface cells burn down to $10^{-10}$ to $10^{-13}$ m. On
  the faster burns the step then collapses, and both strategies stall alike:
  at equal step counts they reach the same simulated time. PMMA at 50 kW/m²
  with solid $\lambda = 0$ ends in `DtLessThanMin` at 225 s (reduced) and
  232.5 s (exact), the ablation deck at 0.25 s with both, and PMMA at 50 and
  100 kW/m² had not reached 5 % of its thickness after an hour of either. In
  that regime the surface mass flux changes sign between outputs with either
  Jacobian, while the gas generation and the thickness stay smooth and agree
  to $10^{-7}$. Parking (§14.3) is what carries a run through burnout; the
  exact Jacobian does not.

Per unit of work the reduced strategy with KLU is 8–17× faster than
`DirectSolve` with UMFPACK on these decks (median CPU time over five
interleaved rounds: 0.057–0.121 of it where the decks run to their end, and
0.107–0.120 for the stalling decks over an equal budget of 4000 integrator
iterations), and allocates 20–8000 times less.

| Strategy | Jacobian assembled | Exact? | Sparsity profile | Verification |
|---|---|---|---|---|
| `DirectSolve` | $\mathsf{A}+\sum C_j$ densified | yes | `:full` | strict |
| `ApproxSolve` | $\mathsf{A}$ only | no | `:ale_reduced` | masked |

> **Implementation note.** Strategy types: src/Jacobian/strategies.jl;
> `sparsity_profile(::Strategy)` maps `DirectSolve` → `:full` and `ApproxSolve` →
> `:ale_reduced`. The emit step `_emit_jacobian!` (src/Jacobian/structured_backend.jl)
> copies $\mathsf{A}$ into the integrator's matrix and densifies the couplings only
> under `:full`.

---

## 15.11 Sparsity construction

### 15.11.1 Two patterns: physics and geometry

The ALE geometry coupling motivates a **two-pattern** model:

- **$S_{\rm phys}$** — every coupling with $T$/$\xi$/$\chi$ columns:
  block-tridiagonal physics ($T\!\leftrightarrow\!T$, $T\!\leftrightarrow\!\xi$,
  $\xi\!\leftrightarrow\!\xi$ at each cell and its neighbors), the Beer–Lambert
  suffix (row $i$ depends on $\xi$-columns $j\ge i$, plus — on the Kirchhoff
  default path — the surface temperature column $T_{\rm col}(n)$ for every row
  $i<n$, the $a_{\rm eff}(T_n,\xi_{\cdot,n})$ transmissivity chain of
  §15.8.7), and the $\chi$-row block
  (cell-local plus the dense $\chi$-column). The ALE dilation is cell-local and
  lies inside the block-tridiagonal physics.
- **$S_{\rm geom}$** — every coupling with $z$ columns (ALE only; empty for
  Eulerian): physics$\to z$ (cells read nodes; reach $\le i+3$ without $\chi$,
  full $1{:}n_{\rm nodes}$ with `WithChi` via the $\theta_A$ path, and full for
  the energy rows under Beer–Lambert, whose optical depth
  $\tau_i=\sum_{k>i}\alpha_k\Delta z_k$ depends on every node above the cell),
  $z\to$physics (node $k$ depends on cells $i<k$, lower-triangular), and $z\to z$
  (lower-triangular including the diagonal, the cumulative nature of mesh
  velocity).

Four emitter functions compose these patterns. `build_sparsity(layout, mode;
profile)` builds $S_{\rm phys}$ from the first three, in this order, and
$S_{\rm geom}$ from the fourth:

1. `_emit_physics_pattern!` — block-tridiagonal;
2. `_emit_radiation_pattern!` — Beer–Lambert suffix (BeerLambert);
3. `_emit_chi_pattern!` — $\chi$-row and dense $\chi$-column (WithChi);
4. `_emit_radiation_pattern!` — Beer–Lambert suffix (BeerLambert).

The integrator's `jac_prototype` is the union $S_{\rm phys}\cup S_{\rm geom}$.
Every row of the union stores its diagonal in every profile, so the shifted
matrix $W$ of §15.5.2 has exactly the pattern of $J$; where no operator writes
the diagonal (the node rows and the $\chi$ rows under `:ale_reduced`) it holds
a structural zero. Patterns are memoized by the fingerprint
$(\mathrm{typeof}(\text{mode}),N_C,n_{\rm cells},n_{\rm nodes},\text{profile})$ in
`_SPARSITY_CACHE`.

### 15.11.2 Profiles

- `:full` (default) — all entries assembled; structured couplings emitted; the
  geometry AD pass runs.
- `:ale_reduced` (`ApproxSolve` only) — drops every structured coupling (the
  z-row prefix-sum blocks, the $\theta_A$ rank-1
  channels, the Beer–Lambert suffix, and the $\chi$-column rank-1) and the
  geometry block, of which only the node rows' structural diagonal remains.
  The masked verifier expects exactly this set of omissions, and counts the
  stored but never-written diagonals of the node and $\chi$ rows
  (`structural_zero_diagonal`) among them.

Only `:full` couples cells an arbitrary distance apart, through the couplings
it materializes: the ALE geometry block and prefix sums, the Beer–Lambert
suffix, and the $\chi$-column rank-1. Each adds $O(N^2)$ entries; every other
pattern couples a cell only to its immediate neighbours.
`materializes_geometry_block(mode, profile)` (`:full` under ALE) selects the
default factorization (§15.5.1).

### 15.11.3 The geometry operator (z-column AD)

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

> **Implementation note.** `build_sparsity` and the four emitters:
> src/Jacobian/sparsity.jl. `GeometryOperator`: src/Jacobian/operators/geometry.jl.
> `FastSparsityMap`: src/Jacobian/fast_sparsity_map.jl.

---

## 15.12 ILU-preconditioned Krylov solvers

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

**Linear-solver family.** Every `solve` attaches a sparse `jac_prototype`, which
pairs with a sparse direct default chosen from its pattern (`KLUFactorization` or
`UMFPACKFactorization`, §15.5.1). An explicit `linear_solver` — one of the
iterative solvers above, `SparspakFactorization()`, or any other LinearSolve
backend, such as `MKLPardisoFactorize()` (with `using LinearSolve` and the Pardiso
/ MKL packages loaded) — always wins.

> **Implementation note.** `ILUGMRESFactorization`, `ILUBiCGSTABFactorization`,
> and the Krylov integration: src/Solver/ilu_gmres.jl.

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

Verification can run as a setup-time check or as a callback (`verify_jacobian = (rtol=, atol=, sample_times=)`)
firing at chosen `sample_times` during a live solve. The callback *warns* on
deviation rather than aborting — a single transient outlier should not kill an
otherwise-good production run; the user reads `max_abs_dev`/`max_rel_dev` to
judge severity.

Depletion handling (§14.3) keeps the state size and the Jacobian's pattern and
caches, so `verify_jacobian` runs through park events; the parked states of every
mode are part of the verifier's own matrix
(`test/jacobian/test_park_jacobian.jl`).

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
  $S_{\rm conv}$ — with every component in the storage term, which is the
  energy-equation convention Pyrolysis.jl adopts. The sign convention
  for heat of reaction differs in *publication* ($h>0$ exothermic in
  ThermaKin/Gpyro vs. $h>0$ endothermic internally here), but the volumetric
  source $Q_{\rm rxn}=-\sum h_r r_r$ is consistent once the convention is mapped
  (overload H1).
- **Mesh.** ThermaKin's dynamic mesh adjustment is the closest analog to our ALE
  with depletion parking; both track a receding surface. Pyrolysis.jl's distinction is the
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
models, Darcy gas transport, ALE with variable cross-section). FDS carries no
in-depth gas, so its energy equation is the condensed-phase limit of ours; it
shares the Arrhenius kinetic form, and its
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
\frac{dT_i}{dt}=\frac{1}{\rho c_p^{\rm eff}{}_{,i}}\Big[-(\nabla\cdot\mathbf{q})_i+S_{{\rm conv},i}+Q_{{\rm rad},i}+Q_{{\rm rxn},i}\Big],
$$
$$
\frac{d\xi_{j,i}}{dt}=-(\nabla\cdot\mathbf{J}_j)_i+S_{j,{\rm rxn},i}-\underbrace{\xi_{j,i}\theta_i}_{\text{ALE, every component}},
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
suffix-sum (Beer–Lambert), and rank-1 ($\theta_A$/$\chi$) couplings, written into
the sparse matrix under `DirectSolve` and dropped under `ApproxSolve`. The linear
solve $W\delta\boldsymbol{u}=-r$ is a sparse factorization of that matrix, and the
whole construction is verified against the dense AD reference.

---

## 15.16 Limitations and open issues

- **Geometry operator under `:ale_reduced`.** The $z$-columns are dropped
  entirely under `ApproxSolve`. Whether *any* residual $z$-coupling survives the
  profile (i.e. whether the profile is guaranteed to leave exactly zero $z$-rows)
  is asserted by the masked verifier rather than proven; the design intent is
  zero, and verification is the safeguard.
- **Chunk-size tuning.** Per-operator `ForwardDiff` chunk sizes are hardcoded
  (e.g. $3(1+N_C)$ for conduction). They are the informationally-necessary
  stencil widths, not adaptively tuned for cache behavior; no auto-tuning is
  present.
- **WithChi × Darcy.** `WithChi` is tested mainly with Fickian transport; the
  interaction of the $\chi$ rank-1 coupling with the full Darcy operator set is
  not exhaustively exercised.
- **Sensitivities across park events.** A park is a discrete rewrite of
  the state (§14.3.4); finite-difference sensitivities are valid only between
  events.
- **ILU drop tolerance.** The default $\tau=0.1$ for the ILU-Krylov solvers has
  not been swept for optimality across pyrolysis regimes; inexact solves can
  raise Newton rejections on Arrhenius transients.
- **Warm-start artifacts.** Reusing a prior final state as fresh initial data can
  perturb finite-difference sensitivities and continuation studies; its effect is
  unquantified (§17.7).
- **Diagnostics cadence.** The diagnostics snapshots are stride-based (every
  `diagnostics_stride` steps), not interpolated to output times, and record the
  sensible energy only; the ledger fields are NaN (§16.8.3).
